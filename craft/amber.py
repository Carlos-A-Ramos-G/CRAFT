"""
craft.amber
Run the post-Gaussian AMBER parameterization pipeline:

  espgen      – extract ESP grid from Gaussian HF log
  resp        – fit RESP charges
  antechamber – assign AMBER atom types and write .ac file
  prepgen     – build residue topology (.prepin)
  parmchk2    – generate missing force-field parameters (.frcmod)
               twice: once for GAFF, once for ff14SB

All tools must be available in $PATH (i.e. AMBER or AmberTools installed).
"""

import os
import re
import shutil
import subprocess
from pathlib import Path


def _run(cmd, cwd):
    """Run a shell command, raise on failure, print output."""
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout.rstrip())
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {result.returncode}):\n"
            f"  {' '.join(cmd)}\n"
            f"{result.stderr}"
        )


def remap_ac_atom_names(ac_path, capped_pdb_path):
    """
    Rename atoms in an antechamber .ac file to match a capped PDB.

    Builds a positional mapping (ac_name[i] → pdb_name[i]) from the order
    of ATOM records and applies it to all ATOM and BOND lines.  Overwrites
    ac_path in place.

    Parameters
    ----------
    ac_path         : str | Path — antechamber .ac file to update
    capped_pdb_path : str | Path — capped PDB whose atom names are the target
    """
    ac_path         = Path(ac_path)
    capped_pdb_path = Path(capped_pdb_path)

    ac_names = [
        ln.split()[2]
        for ln in ac_path.read_text().splitlines()
        if ln.startswith('ATOM')
    ]
    pdb_names = [
        ln.split()[2]
        for ln in capped_pdb_path.read_text().splitlines()
        if ln.startswith('ATOM')
    ]

    if len(ac_names) != len(pdb_names):
        raise ValueError(
            f"Atom count mismatch: {ac_path.name} has {len(ac_names)} atoms, "
            f"{capped_pdb_path.name} has {len(pdb_names)} atoms."
        )

    name_map = dict(zip(ac_names, pdb_names))

    out_lines = []
    for line in ac_path.read_text().splitlines(keepends=True):
        nl = '\n' if line.endswith('\n') else ''
        s  = line.rstrip('\n')

        if s.startswith('ATOM'):
            # Groups: (ATOM + serial) | name | (resname + rest)
            m = re.match(r'^(ATOM\s+\d+)\s+(\S+)\s+(\S+\s+.*)$', s)
            if m:
                new_name = name_map.get(m.group(2), m.group(2))
                # Left-justify name in 4-char field; 2 spaces before name
                s = f"{m.group(1)}  {new_name:<4}{m.group(3)}"

        elif s.startswith('BOND'):
            # Groups: (BOND + indices) | name1 | (spaces) | name2 | (rest)
            m = re.match(r'^(BOND\s+\d+\s+\d+\s+\d+\s+\d+\s+)(\S+)(\s+)(\S+)(.*)$', s)
            if m:
                new1 = name_map.get(m.group(2), m.group(2))
                new2 = name_map.get(m.group(4), m.group(4))
                s = m.group(1) + new1 + m.group(3) + new2 + m.group(5)

        out_lines.append(s + nl)

    ac_path.write_text(''.join(out_lines))
    print(f"  Remapped {len(name_map)} atom names in {ac_path.name}")
    return name_map


def run_amber_pipeline(hf_log, resname, charge, mc_file,
                       workdir='.', atom_type='amber', ff14sb=True,
                       capped_pdb=None):
    """
    Run the full post-Gaussian parameterization pipeline.

    Parameters
    ----------
    hf_log     : str | Path — Gaussian HF/ESP log (MEO_hf.log)
    resname    : str        — three-letter residue code (e.g. 'MEO')
    charge     : int        — net molecular charge of the capped model
    mc_file    : str | Path — prepgen main-chain file (.mc)
    workdir    : str | Path — directory where all files are written (default: '.')
    atom_type  : str        — antechamber -at flag: 'amber', 'gaff', or 'gaff2'
    ff14sb     : bool       — also run parmchk2 with ff14SB parameters
    capped_pdb : str | Path | None — capped PDB used to rename atoms in the .ac
                 file after antechamber; defaults to {resname}_capped.pdb in
                 workdir; pass None to skip renaming
    """
    hf_log  = str(Path(hf_log).resolve())
    mc_file = str(Path(mc_file).resolve())
    wd      = str(Path(workdir).resolve())

    ac_file      = f"{resname}.ac"
    prepin_file  = f"{resname}.prepin"
    gaff_frcmod  = f"{resname}_gaff.frcmod"
    ff14sb_frcmod= f"{resname}_ff14SB.frcmod"

    amberhome = os.environ.get('AMBERHOME', '')
    parm10    = str(Path(amberhome) / 'dat/leap/parm/parm10.dat') if amberhome else 'parm10.dat'

    print("\n-- espgen ------------------------------------------------------------")
    _run(['espgen', '-i', hf_log, '-o', 'esp.dat'], cwd=wd)

    print("\n-- resp --------------------------------------------------------------")
    _run([
        'resp', '-O',
        '-i', 'resp.in',
        '-o', 'resp.out',
        '-p', 'resp.pch',
        '-t', 'resp.chg',
        '-q', 'resp.qin',
        '-e', 'esp.dat',
    ], cwd=wd)

    print("\n-- antechamber -------------------------------------------------------")
    _run([
        'antechamber',
        '-fi', 'gout',
        '-i',  hf_log,
        '-bk', resname,
        '-fo', 'ac',
        '-o',  ac_file,
        '-c',  'rc',
        '-cf', 'resp.chg',
        '-at', atom_type,
        '-nc', str(charge),
    ], cwd=wd)

    print("\n-- remap atom names -------------------------------------------------")
    _capped = Path(capped_pdb) if capped_pdb else Path(wd) / f"{resname}_capped.pdb"
    if _capped.exists():
        remap_ac_atom_names(Path(wd) / ac_file, _capped)
    else:
        print(f"  {_capped.name} not found — atom names in {ac_file} not remapped.")

    # prepgen has a short fixed-length path buffer (~256 chars); long absolute
    # paths are silently truncated. Copy the .mc file into workdir and pass
    # only the basename to avoid this.
    mc_local = Path(mc_file).name
    mc_dst = Path(wd) / mc_local
    if mc_dst.resolve() != Path(mc_file).resolve():
        shutil.copy(mc_file, mc_dst)

    print("\n-- prepgen -----------------------------------------------------------")
    _run([
        'prepgen',
        '-i',  ac_file,
        '-o',  prepin_file,
        '-m',  mc_local,
        '-rn', resname,
    ], cwd=wd)

    print("\n-- parmchk2 (GAFF) ---------------------------------------------------")
    _run([
        'parmchk2',
        '-i', ac_file,
        '-f', 'ac',
        '-o', gaff_frcmod,
    ], cwd=wd)

    if ff14sb:
        print("\n-- parmchk2 (ff14SB) -------------------------------------------------")
        _run([
            'parmchk2',
            '-i', ac_file,
            '-f', 'ac',
            '-o', ff14sb_frcmod,
            '-a', 'Y',
            '-p', parm10,
        ], cwd=wd)

    output_files = [ac_file, prepin_file, gaff_frcmod]
    if ff14sb:
        output_files.append(ff14sb_frcmod)

    print(f"\nDone. Output files in {wd}:")
    for f in output_files:
        p = Path(wd) / f
        status = "✓" if p.exists() else "✗ MISSING"
        print(f"  {status}  {f}")
