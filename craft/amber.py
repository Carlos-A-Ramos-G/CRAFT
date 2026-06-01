"""
craft.amber
Run the post-Gaussian AMBER parameterization pipeline:

  espgen      - extract ESP grid from Gaussian HF log
  resp        - fit RESP charges
  antechamber - assign AMBER atom types and write .ac file
  prepgen     - build residue topology (.prepin)
  parmchk2    - generate missing force-field parameters (.frcmod)
               twice: once for GAFF, once for the target protein FF

All tools must be available in $PATH (i.e. AMBER or AmberTools installed).

Output file names carry a position prefix for terminal variants:
  middle : {resname}.ac, {resname}.prepin, {resname}_gaff.frcmod, ...
  cterm  : C{resname}.ac, C{resname}.prepin, C{resname}_gaff.frcmod, ...
  nterm  : N{resname}.ac, N{resname}.prepin, N{resname}_gaff.frcmod, ...
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

    Builds a positional mapping (ac_name[i] -> pdb_name[i]) from the order
    of ATOM records and applies it to all ATOM and BOND lines.  Overwrites
    ac_path in place.

    Parameters
    ----------
    ac_path         : str | Path -- antechamber .ac file to update
    capped_pdb_path : str | Path -- capped PDB whose atom names are the target
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
            m = re.match(r'^(ATOM\s+\d+)\s+(\S+)\s+(\S+\s+.*)$', s)
            if m:
                new_name = name_map.get(m.group(2), m.group(2))
                s = f"{m.group(1)}  {new_name:<4}{m.group(3)}"

        elif s.startswith('BOND'):
            m = re.match(r'^(BOND\s+\d+\s+\d+\s+\d+\s+\d+\s+)(\S+)(\s+)(\S+)(.*)$', s)
            if m:
                new1 = name_map.get(m.group(2), m.group(2))
                new2 = name_map.get(m.group(4), m.group(4))
                s = m.group(1) + new1 + m.group(3) + new2 + m.group(5)

        out_lines.append(s + nl)

    ac_path.write_text(''.join(out_lines))
    print(f"  Remapped {len(name_map)} atom names in {ac_path.name}")
    return name_map


PARM_FILES = {
    'ff14SB': 'parm10.dat',
    'ff19SB': 'parm19.dat',
}


def run_amber_pipeline(hf_log, resname, charge, mc_file,
                       workdir='.', atom_type='amber', forcefield='ff14SB',
                       capped_pdb=None, position='middle'):
    """
    Run the full post-Gaussian parameterization pipeline.

    Parameters
    ----------
    hf_log      : str | Path -- Gaussian HF/ESP log (MEO_hf.log)
    resname     : str        -- three-letter residue code (e.g. 'MEO')
    charge      : int        -- net molecular charge of the capped model
    mc_file     : str | Path -- prepgen main-chain file (.mc)
    workdir     : str | Path -- directory where all files are written
    atom_type   : str        -- antechamber -at flag: 'amber', 'gaff', or 'gaff2'
    forcefield  : str | None -- protein FF for the second parmchk2 run:
                                'ff14SB' (parm10.dat), 'ff19SB' (parm19.dat),
                                or None to skip
    capped_pdb  : str | Path | None -- capped PDB for atom name remapping;
                  defaults to {base}_capped.pdb in workdir; pass None to skip
    position    : str        -- 'middle', 'cterm', or 'nterm'; determines the
                                output file name prefix (C/N/none)
    """
    if position not in ('middle', 'cterm', 'nterm'):
        raise ValueError(
            f"position must be 'middle', 'cterm', or 'nterm'; got {position!r}")

    prefix = {'middle': '', 'cterm': 'C', 'nterm': 'N'}[position]
    base   = f"{prefix}{resname}"

    hf_log  = str(Path(hf_log).resolve())
    mc_file = str(Path(mc_file).resolve())
    wd      = str(Path(workdir).resolve())

    ac_file     = f"{base}.ac"
    prepin_file = f"{base}.prepin"
    gaff_frcmod = f"{base}_gaff.frcmod"
    ff_frcmod   = f"{base}_{forcefield}.frcmod" if forcefield else None

    amberhome = os.environ.get('AMBERHOME', '')
    if forcefield:
        parm_name = PARM_FILES.get(forcefield)
        if parm_name is None:
            raise ValueError(
                f"Unknown forcefield {forcefield!r}. "
                f"Supported: {list(PARM_FILES)}"
            )
        parm_file = str(Path(amberhome) / 'dat/leap/parm' / parm_name) if amberhome else parm_name

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
    if capped_pdb:
        _capped = Path(capped_pdb)
    else:
        _capped = Path(wd) / f"{base}_capped.pdb"
        if not _capped.exists():
            _capped = Path(wd) / f"{resname}_capped.pdb"

    if _capped.exists():
        remap_ac_atom_names(Path(wd) / ac_file, _capped)
    else:
        print(f"  {_capped.name} not found -- atom names in {ac_file} not remapped.")

    # prepgen has a short fixed-length path buffer (~256 chars); long absolute
    # paths are silently truncated. Copy the .mc file into workdir and pass
    # only the basename to avoid this.
    mc_local = Path(mc_file).name
    mc_dst   = Path(wd) / mc_local
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

    if forcefield:
        print(f"\n-- parmchk2 ({forcefield}) {'-' * (51 - len(forcefield))}")
        _run([
            'parmchk2',
            '-i', ac_file,
            '-f', 'ac',
            '-o', ff_frcmod,
            '-a', 'Y',
            '-p', parm_file,
        ], cwd=wd)

    output_files = [ac_file, prepin_file, gaff_frcmod]
    if forcefield:
        output_files.append(ff_frcmod)

    print(f"\nDone. Output files in {wd}:")
    for f in output_files:
        p = Path(wd) / f
        status = "ok" if p.exists() else "MISSING"
        print(f"  [{status}]  {f}")
