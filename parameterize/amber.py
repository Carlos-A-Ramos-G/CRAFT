"""
parameterize.amber
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


def run_amber_pipeline(hf_log, resname, charge, mc_file,
                       workdir='.', atom_type='amber', ff14sb=True):
    """
    Run the full post-Gaussian parameterization pipeline.

    Parameters
    ----------
    hf_log    : str | Path — Gaussian HF/ESP log (MEO_hf.log)
    resname   : str        — three-letter residue code (e.g. 'MEO')
    charge    : int        — net molecular charge of the capped model
    mc_file   : str | Path — prepgen main-chain file (.mc)
    workdir   : str | Path — directory where all files are written (default: '.')
    atom_type : str        — antechamber -at flag: 'amber', 'gaff', or 'gaff2'
    ff14sb    : bool       — also run parmchk2 with ff14SB parameters
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

    print("\n── espgen ────────────────────────────────────────────────────────────")
    _run(['espgen', '-i', hf_log, '-o', 'esp.dat'], cwd=wd)

    print("\n── resp ──────────────────────────────────────────────────────────────")
    _run([
        'resp', '-O',
        '-i', 'resp.in',
        '-o', 'resp.out',
        '-p', 'resp.pch',
        '-t', 'resp.chg',
        '-q', 'resp.qin',
        '-e', 'esp.dat',
    ], cwd=wd)

    print("\n── antechamber ───────────────────────────────────────────────────────")
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

    # prepgen has a short fixed-length path buffer (~256 chars); long absolute
    # paths are silently truncated. Copy the .mc file into workdir and pass
    # only the basename to avoid this.
    mc_local = Path(mc_file).name
    mc_dst = Path(wd) / mc_local
    if mc_dst.resolve() != Path(mc_file).resolve():
        shutil.copy(mc_file, mc_dst)

    print("\n── prepgen ───────────────────────────────────────────────────────────")
    _run([
        'prepgen',
        '-i',  ac_file,
        '-o',  prepin_file,
        '-m',  mc_local,
        '-rn', resname,
    ], cwd=wd)

    print("\n── parmchk2 (GAFF) ───────────────────────────────────────────────────")
    _run([
        'parmchk2',
        '-i', ac_file,
        '-f', 'ac',
        '-o', gaff_frcmod,
    ], cwd=wd)

    if ff14sb:
        print("\n── parmchk2 (ff14SB) ─────────────────────────────────────────────────")
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
