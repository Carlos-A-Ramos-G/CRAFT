"""
craft.cli
Console-script entry points registered in pyproject.toml.

Each function here is also called by its matching top-level convenience
script (run.py, make_hf_input.py, ...) so both invocation styles work:

    python run.py                 # clone-and-run (no install needed)
    craft-run                     # after pip install .
"""

import sys


# ---------------------------------------------------------------------------
# craft-run   (Phase 1)
# ---------------------------------------------------------------------------

def run():
    """Cap termini and generate all pre-Gaussian inputs."""
    import argparse
    import yaml
    from pathlib import Path
    from craft import cap, get_resname, write_com, write_resp_in, write_resp_qin, write_mc
    from craft.gaussian import NPROC_DEFAULT, MEM_DEFAULT, ROUTE_DEFAULT

    parser = argparse.ArgumentParser(
        description="Cap termini and generate all pre-Gaussian inputs (Phase 1)",
    )
    parser.add_argument('--config', default='config.yaml',
                        help='Config file (default: config.yaml)')
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    res       = cfg['residue']
    input_pdb = res['input_pdb']
    charge    = res.get('charge', 0)
    mult      = res.get('multiplicity', 1)
    position  = res.get('position', 'middle')
    resname   = get_resname(input_pdb)

    if position not in ('middle', 'cterm', 'nterm'):
        sys.exit(f"Error: residue.position must be 'middle', 'cterm', or 'nterm'; "
                 f"got {position!r}")

    prefix = {'middle': '', 'cterm': 'C', 'nterm': 'N'}[position]
    base   = f"{prefix}{resname}"

    workdir = Path(resname) / position
    workdir.mkdir(parents=True, exist_ok=True)

    cap_cfg    = cfg.get('cap') or {}
    capped_pdb = cap_cfg.get('output_pdb') or str(workdir / f"{base}_capped.pdb")

    print("=" * 60)
    print(f"Step 1 -- Cap termini  [{position}]")
    print("=" * 60)
    cap(input_pdb, capped_pdb, position=position)

    gauss_cfg = cfg.get('gaussian_opt') or cfg.get('gaussian') or {}
    opt_com   = gauss_cfg.get('output_com') or str(workdir / f"{base}_opt.com")

    print()
    print("=" * 60)
    print("Step 2 -- Geometry-optimisation Gaussian input")
    print("=" * 60)
    write_com(
        capped_pdb, opt_com,
        charge=charge, mult=mult,
        nproc=gauss_cfg.get('nproc', NPROC_DEFAULT),
        mem=gauss_cfg.get('mem', MEM_DEFAULT),
        route=gauss_cfg.get('route', ROUTE_DEFAULT),
    )

    print()
    print("=" * 60)
    print("Steps 3-5 -- RESP and prepgen inputs")
    print("=" * 60)
    write_resp_in( capped_pdb, charge, resname, str(workdir / 'resp.in'),
                   position=position)
    write_resp_qin(capped_pdb, str(workdir / 'resp.qin'),
                   position=position)
    write_mc(      capped_pdb, charge, str(workdir / f"{base}.mc"),
                   position=position)

    print()
    print(f"All Phase 1 outputs written to: {workdir}/")
    print()
    print("Next steps:")
    print(f"  1. Submit {workdir}/{base}_opt.com to HPC")
    print(f"  2. Copy {base}_opt.log into {workdir}/, then:")
    print(f"       craft-hf-input {workdir}/{base}_opt.log")
    print(f"  3. Submit {workdir}/{base}_hf.com to HPC")
    print(f"  4. Copy {base}_hf.log into {workdir}/, then:")
    print(f"       craft-amber {workdir}/{base}_hf.log")
    print(f"\n  Residue '{resname}'  position '{position}'  read from {input_pdb}")


# ---------------------------------------------------------------------------
# craft-hf-input   (Phase 2b)
# ---------------------------------------------------------------------------

def hf_input():
    """Extract optimised geometry and write HF/6-31G(d) single-point input."""
    import argparse
    import yaml
    from pathlib import Path
    from craft import parse_opt_log, write_hf_com, get_resname
    from craft.gaussian import NPROC_DEFAULT, MEM_DEFAULT, HF_ROUTE_DEFAULT

    parser = argparse.ArgumentParser(
        description="Gaussian opt log -> HF/6-31G(d) single-point .com",
    )
    parser.add_argument('log',            help='Gaussian optimisation log (e.g. KME3_opt.log)')
    parser.add_argument('com', nargs='?', help='Output .com  (default: <base>_hf.com)')
    parser.add_argument('--config',       default='config.yaml',
                        help='Config file (default: config.yaml)')
    parser.add_argument('-c', '--charge', type=int, default=None)
    parser.add_argument('-m', '--mult',   type=int, default=None)
    parser.add_argument('-n', '--nproc',  type=int, default=None)
    parser.add_argument('--mem',          default=None)
    args = parser.parse_args()

    try:
        cfg = yaml.safe_load(open(args.config))
    except FileNotFoundError:
        cfg = {}
    res_cfg = cfg.get('residue', {})
    g_cfg   = cfg.get('gaussian_hf', {}) or cfg.get('gaussian', {}) or {}

    charge   = args.charge if args.charge is not None else res_cfg.get('charge', 0)
    mult     = args.mult   if args.mult   is not None else res_cfg.get('multiplicity', 1)
    nproc    = args.nproc  if args.nproc  is not None else g_cfg.get('nproc', NPROC_DEFAULT)
    mem      = args.mem    if args.mem    is not None else g_cfg.get('mem', MEM_DEFAULT)
    position = res_cfg.get('position', 'middle')

    log_path = Path(args.log)
    workdir  = log_path.parent

    if res_cfg.get('input_pdb'):
        resname = get_resname(res_cfg['input_pdb'])
    else:
        stem = log_path.stem
        for tag in ('_opt', '_hf'):
            stem = stem.replace(tag, '')
        if position == 'nterm' and stem.startswith('N'):
            stem = stem[1:]
        elif position == 'cterm' and stem.startswith('C'):
            stem = stem[1:]
        resname = stem

    prefix   = {'middle': '', 'cterm': 'C', 'nterm': 'N'}[position]
    base     = f"{prefix}{resname}"
    com_path = args.com or str(workdir / f"{base}_hf.com")

    print(f"Parsing optimised geometry from {args.log} ...")
    atoms_xyz = parse_opt_log(args.log)
    print(f"  Found {len(atoms_xyz)} atoms in final Standard orientation block")

    print(f"Writing HF/6-31G(d) single-point input ...")
    write_hf_com(atoms_xyz, com_path, base,
                 charge=charge, mult=mult, nproc=nproc, mem=mem,
                 route=g_cfg.get('route') or HF_ROUTE_DEFAULT)


# ---------------------------------------------------------------------------
# craft-amber   (Phase 3)
# ---------------------------------------------------------------------------

def amber():
    """Run espgen -> resp -> antechamber -> prepgen -> parmchk2."""
    import argparse
    import yaml
    from pathlib import Path
    from craft import run_amber_pipeline, get_resname

    parser = argparse.ArgumentParser(
        description="Run AMBER parameterization pipeline from Gaussian HF log",
    )
    parser.add_argument('log',             help='Gaussian HF/ESP log')
    parser.add_argument('--config',        default='config.yaml',
                        help='Config file (default: config.yaml)')
    parser.add_argument('--resname',       default=None)
    parser.add_argument('-c', '--charge',  type=int, default=None)
    parser.add_argument('--workdir',       default='.')
    args = parser.parse_args()

    try:
        cfg = yaml.safe_load(open(args.config))
    except FileNotFoundError:
        cfg = {}
    res_cfg = cfg.get('residue', {})
    amb_cfg = cfg.get('amber',   {}) or {}
    cap_cfg = cfg.get('cap',     {}) or {}

    charge      = args.charge if args.charge is not None else res_cfg.get('charge', 0)
    position    = res_cfg.get('position', 'middle')
    forcefield  = amb_cfg.get('forcefield', 'ff14SB')

    workdir = (args.workdir if args.workdir != '.'
               else (amb_cfg.get('workdir') or str(Path(args.log).resolve().parent)))

    if args.resname:
        resname = args.resname
    elif res_cfg.get('input_pdb'):
        resname = get_resname(res_cfg['input_pdb'])
    else:
        stem = Path(args.log).stem
        for tag in ('_hf', '_opt'):
            stem = stem.replace(tag, '')
        if position == 'nterm' and stem.startswith('N'):
            stem = stem[1:]
        elif position == 'cterm' and stem.startswith('C'):
            stem = stem[1:]
        resname = stem

    prefix  = {'middle': '', 'cterm': 'C', 'nterm': 'N'}[position]
    base    = f"{prefix}{resname}"
    mc_file = Path(workdir) / f"{base}.mc"
    if not mc_file.exists():
        sys.exit(f"Error: {mc_file} not found -- run 'craft-run' first.")

    capped_pdb = cap_cfg.get('output_pdb') or None

    print(f"Residue  : {resname}  (charge {charge:+d})")
    print(f"Position : {position}")
    print(f"Log      : {args.log}")
    print(f"MC file  : {mc_file}")
    print()

    run_amber_pipeline(
        hf_log     = args.log,
        resname    = resname,
        charge     = charge,
        mc_file    = str(mc_file),
        workdir    = workdir,
        atom_type   = amb_cfg.get('atom_type', 'amber'),
        forcefield  = forcefield,
        capped_pdb = capped_pdb,
        position   = position,
    )


# ---------------------------------------------------------------------------
# craft-slurm   (SLURM script generator)
# ---------------------------------------------------------------------------

def slurm():
    """Generate a SLURM batch script for the full pipeline."""
    import argparse
    import yaml
    from pathlib import Path
    from craft import write_slurm, get_resname

    parser = argparse.ArgumentParser(
        description="Generate a SLURM batch script for the full pipeline",
    )
    parser.add_argument('--config', default='config.yaml',
                        help='Config file (default: config.yaml)')
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    cfg = yaml.safe_load(open(config_path))

    proj_root = Path.cwd().resolve()
    res_cfg   = cfg.get('residue', {})
    input_pdb = res_cfg.get('input_pdb', '')
    position  = res_cfg.get('position', 'middle')
    resname   = get_resname(input_pdb) if input_pdb else 'residue'
    prefix    = {'middle': '', 'cterm': 'C', 'nterm': 'N'}[position]
    base      = f"{prefix}{resname}"

    workdir = proj_root / resname / position
    workdir.mkdir(parents=True, exist_ok=True)

    output = str(workdir / f"{base}_craft.sh")

    write_slurm(cfg, output, proj_root, workdir, position, config_path=config_path)


# ---------------------------------------------------------------------------
# craft-check   (environment checker)
# ---------------------------------------------------------------------------

def check():
    """Check that all required tools and packages are available."""
    from craft.check import check_env
    check_env()
