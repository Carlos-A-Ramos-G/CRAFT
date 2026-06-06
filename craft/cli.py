"""
craft.cli
Console-script entry points registered in pyproject.toml.

Each command detects from the config file whether to run the single-residue
or two-residue reaction workflow:
  - Single residue : config has a `residue` key
  - Reaction       : config has `residue1`, `residue2`, and `reaction` keys

Both invocation styles work after pip install .:
    craft-run --config config.yaml
    craft-run --config config_react.yaml
"""

import sys


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _load_config(path):
    """Load and return a YAML config, or exit with a clear message if not found."""
    import yaml
    from pathlib import Path
    p = Path(path)
    if not p.exists():
        sys.exit(
            f"Error: config file '{path}' not found.\n"
            f"\n"
            f"  All craft-* commands must be run from the directory that contains\n"
            f"  your config.yaml, or you must supply an explicit path:\n"
            f"\n"
            f"    craft-run --config /path/to/config.yaml\n"
            f"    craft-hf-input <opt.log> --config /path/to/config.yaml\n"
            f"    craft-amber   <hf.log>  --config /path/to/config.yaml\n"
            f"    craft-slurm            --config /path/to/config.yaml\n"
            f"\n"
            f"  Copy config.yaml (single residue) or config_react.yaml (two bonded\n"
            f"  residues) from the CRAFT repository into your working directory and\n"
            f"  fill in the values for your system."
        )
    with open(p) as f:
        return yaml.safe_load(f)


def _resnames_from_pdb(pdb_path):
    """Return non-cap residue names from pdb_path in order of first appearance."""
    from craft.cap import parse_pdb
    atoms = parse_pdb(str(pdb_path))
    seen = []
    for a in atoms:
        rn = a['resName']
        if rn not in ('ACE', 'NME') and rn not in seen:
            seen.append(rn)
    return seen


def _resolve_react_resnames(cfg, workdir=None):
    """
    Return (resname1, resname2) from a reaction config.

    Resolution order:
      1. reaction.combined_pdb  (takes priority; auto-detect residue names from file)
      2. residue1.input_pdb / residue2.input_pdb  (both must be set)
      3. glob *_react_combined.pdb in workdir  (phase-3 fallback; requires workdir)

    A warning is printed when both combined_pdb and input_pdb fields are set,
    since they are mutually exclusive — combined_pdb wins.
    """
    from craft import get_resname
    from pathlib import Path

    res1_cfg = cfg['residue1']
    res2_cfg = cfg['residue2']
    pdb1 = res1_cfg.get('input_pdb')
    pdb2 = res2_cfg.get('input_pdb')
    combined_pdb = cfg.get('reaction', {}).get('combined_pdb')

    if combined_pdb:
        if pdb1 or pdb2:
            print("  Warning: reaction.combined_pdb is set — "
                  "residue1.input_pdb and residue2.input_pdb will be ignored.")
        rns = _resnames_from_pdb(combined_pdb)
        if len(rns) != 2:
            sys.exit(
                f"Error: expected exactly 2 residue names in {combined_pdb}, "
                f"found {len(rns)}: {rns}. "
                f"Set residue1.input_pdb and residue2.input_pdb instead."
            )
        return tuple(rns)

    if pdb1 and pdb2:
        return get_resname(pdb1), get_resname(pdb2)

    if workdir is not None:
        candidates = list(Path(workdir).glob('*_react_combined.pdb'))
        if len(candidates) == 1:
            rns = _resnames_from_pdb(candidates[0])
            if len(rns) != 2:
                sys.exit(
                    f"Error: expected 2 residue names in {candidates[0]}, "
                    f"found {len(rns)}: {rns}."
                )
            return tuple(rns)

    sys.exit(
        "Error: residue1.input_pdb and residue2.input_pdb are required "
        "when reaction.combined_pdb is not set."
    )


# ---------------------------------------------------------------------------
# craft-run   (Phase 1)
# ---------------------------------------------------------------------------

def _run_react(cfg, args):
    """Phase 1 body for a two-residue reaction."""
    import json
    from pathlib import Path
    from craft import cap
    from craft.gaussian import NPROC_DEFAULT, MEM_DEFAULT
    from craft.react import (assemble_react_pdb, _split_combined_pdb,
                              write_react_com, write_react_resp_in,
                              write_react_resp_qin, write_react_mc)

    res1_cfg = cfg['residue1']
    res2_cfg = cfg['residue2']
    rxn_cfg  = cfg['reaction']
    g_cfg    = cfg.get('gaussian_opt') or cfg.get('gaussian') or {}

    pdb1      = res1_cfg.get('input_pdb')
    pdb2      = res2_cfg.get('input_pdb')
    charge1   = res1_cfg.get('charge', 0)
    charge2   = res2_cfg.get('charge', 0)
    position1 = res1_cfg.get('position', 'middle')
    position2 = res2_cfg.get('position', 'middle')
    combined_input = rxn_cfg.get('combined_pdb')

    if position1 != 'middle' or position2 != 'middle':
        sys.exit(
            f"Error: the reaction workflow currently only supports position='middle' "
            f"for both residues (cterm/nterm will be added in a future release). "
            f"Got residue1={position1!r}, residue2={position2!r}."
        )

    resname1, resname2 = _resolve_react_resnames(cfg)

    atom1        = rxn_cfg['atom1']
    atom2        = rxn_cfg['atom2']
    bond_length  = rxn_cfg.get('bond_length', 1.5)
    total_charge = charge1 + charge2

    base    = f"{resname1}_{resname2}_react"
    workdir = Path(f"{resname1}_{resname2}")
    sub1    = workdir / resname1
    sub2    = workdir / resname2
    workdir.mkdir(parents=True, exist_ok=True)
    sub1.mkdir(exist_ok=True)
    sub2.mkdir(exist_ok=True)

    capped1         = str(sub1 / f"{resname1}_capped.pdb")
    capped2         = str(sub2 / f"{resname2}_capped.pdb")
    combined_pdb    = str(workdir / f"{base}_combined.pdb")
    com_path        = str(workdir / f"{base}_opt.com")
    mc1_path        = str(sub1 / f"{resname1}.mc")
    mc2_path        = str(sub2 / f"{resname2}.mc")
    resp_in_path    = str(workdir / 'resp.in')
    resp_qin_path   = str(workdir / 'resp.qin')
    rename_map_path = str(workdir / 'rename_map.json')

    if combined_input:
        print("=" * 60)
        print("Steps 1+2 -- Cap residues from pre-assembled structure")
        print("=" * 60)
        split1 = str(sub1 / f"{resname1}_split.pdb")
        split2 = str(sub2 / f"{resname2}_split.pdb")
        _split_combined_pdb(combined_input, resname1, resname2, split1, split2)
        cap(split1, capped1, position=position1)
        cap(split2, capped2, position=position2)

        print()
        print("=" * 60)
        print("Step 3 -- Assemble with caps, preserving pre-assembled geometry")
        print("=" * 60)
        _, _, rename_map = assemble_react_pdb(
            capped1, capped2, atom1, atom2,
            skip_reposition=True, output=combined_pdb)
    else:
        print("=" * 60)
        print(f"Step 1 -- Cap {resname1}  [{position1}]")
        print("=" * 60)
        cap(pdb1, capped1, position=position1)

        print()
        print("=" * 60)
        print(f"Step 2 -- Cap {resname2}  [{position2}]")
        print("=" * 60)
        cap(pdb2, capped2, position=position2)

        print()
        print("=" * 60)
        print(f"Step 3 -- Assemble combined model  ({atom1}—{atom2}, {bond_length} Å)")
        print("=" * 60)
        _, _, rename_map = assemble_react_pdb(
            capped1, capped2, atom1, atom2, bond_length, combined_pdb)

    with open(rename_map_path, 'w') as f:
        json.dump(rename_map, f, indent=2)
    if rename_map:
        print(f"  rename_map.json : {len(rename_map)} atom(s) renamed in {resname2}")

    print()
    print("=" * 60)
    print("Step 4 -- Frozen-backbone Gaussian opt input")
    print("=" * 60)
    write_react_com(
        combined_pdb, com_path,
        charge=total_charge, mult=1,
        nproc=g_cfg.get('nproc', NPROC_DEFAULT),
        mem=g_cfg.get('mem', MEM_DEFAULT),
        route=g_cfg.get('route', "#P b3lyp/6-31g* opt(modredundant)"),
        freeze_backbone=g_cfg.get('freeze_backbone', True),
    )

    print()
    print("=" * 60)
    print("Step 5 -- MC and RESP input files")
    print("=" * 60)
    write_react_mc(combined_pdb, 2, charge1, mc1_path, position1)
    write_react_mc(combined_pdb, 5, charge2, mc2_path, position2)
    write_react_resp_in(combined_pdb, total_charge, resname1, resname2,
                        resp_in_path, position1, position2)
    write_react_resp_qin(combined_pdb, resp_qin_path, position1, position2)

    print()
    print(f"All Phase 1 outputs written to: {workdir}/")
    print(f"\nNext steps:")
    print(f"  1. Submit {workdir}/{base}_opt.com to HPC")
    print(f"  2. Copy {base}_opt.log into {workdir}/, then:")
    print(f"       craft-hf-input {workdir}/{base}_opt.log "
          f"--charge {total_charge} --config {args.config}")
    print(f"  3. Submit {workdir}/{base}_hf.com to HPC")
    print(f"  4. Copy {base}_hf.log into {workdir}/, then:")
    print(f"       craft-amber {workdir}/{base}_hf.log --config {args.config}")


def run():
    """Cap residue(s) and generate all pre-Gaussian inputs (Phase 1)."""
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description="Cap residue(s) and generate all pre-Gaussian inputs (Phase 1). "
                    "Detects single-residue or reaction mode from the config file.",
    )
    parser.add_argument('--config', default='config.yaml',
                        help='Config file (default: config.yaml)')
    args = parser.parse_args()

    cfg = _load_config(args.config)

    if 'residue1' in cfg:
        _run_react(cfg, args)
        return

    # -- single-residue path ---------------------------------------------------
    from craft import cap, get_resname, write_com, write_resp_in, write_resp_qin, write_mc
    from craft.gaussian import NPROC_DEFAULT, MEM_DEFAULT, ROUTE_DEFAULT

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
# craft-hf-input   (Phase 2b) — shared, unchanged
# ---------------------------------------------------------------------------

def hf_input():
    """Extract optimised geometry and write HF/6-31G(d) single-point input."""
    import argparse
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

    cfg = _load_config(args.config)
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

def _amber_react(cfg, args):
    """Phase 3 body for a two-residue reaction."""
    import json
    from pathlib import Path
    from craft.react import run_react_amber_pipeline

    res1_cfg = cfg['residue1']
    res2_cfg = cfg['residue2']
    amb_cfg  = cfg.get('amber', {}) or {}

    charge1      = res1_cfg.get('charge', 0)
    charge2      = res2_cfg.get('charge', 0)
    total_charge = args.charge if args.charge is not None else charge1 + charge2

    workdir = (args.workdir if args.workdir != '.'
               else str(Path(args.log).resolve().parent))

    resname1, resname2 = _resolve_react_resnames(cfg, workdir=workdir)

    base         = f"{resname1}_{resname2}_react"
    mc1_file     = str(Path(workdir) / resname1 / f"{resname1}.mc")
    mc2_file     = str(Path(workdir) / resname2 / f"{resname2}.mc")
    combined_pdb = str(Path(workdir) / f"{base}_combined.pdb")

    rename_map      = {}
    rename_map_path = Path(workdir) / 'rename_map.json'
    if rename_map_path.exists():
        with open(rename_map_path) as f:
            rename_map = json.load(f)
    else:
        print("  Warning: rename_map.json not found -- residue2 atom names "
              "may not be restored in the .prepin. Run craft-run first.")

    print(f"Residues     : {resname1} (charge {charge1:+d})"
          f" + {resname2} (charge {charge2:+d})")
    print(f"Total charge : {total_charge:+d}")
    print(f"Log          : {args.log}")
    if rename_map:
        print(f"Rename map   : {len(rename_map)} atom(s) will be restored "
              f"in {resname2}.prepin")
    print()

    run_react_amber_pipeline(
        hf_log       = args.log,
        resname1     = resname1,
        resname2     = resname2,
        total_charge = total_charge,
        mc_file1     = mc1_file,
        mc_file2     = mc2_file,
        workdir      = workdir,
        atom_type    = amb_cfg.get('atom_type', 'amber'),
        forcefield   = amb_cfg.get('forcefield', 'ff14SB'),
        combined_pdb = combined_pdb,
        rename_map   = rename_map,
    )


def amber():
    """Run espgen -> resp -> antechamber -> prepgen -> parmchk2 (Phase 3)."""
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description="Run AMBER parameterization pipeline from Gaussian HF log. "
                    "Detects single-residue or reaction mode from the config file.",
    )
    parser.add_argument('log',            help='Gaussian HF/ESP log')
    parser.add_argument('--config',       default='config.yaml',
                        help='Config file (default: config.yaml)')
    parser.add_argument('--resname',      default=None,
                        help='Override residue name (single-residue mode only)')
    parser.add_argument('-c', '--charge', type=int, default=None)
    parser.add_argument('--workdir',      default='.')
    args = parser.parse_args()

    cfg = _load_config(args.config)

    if 'residue1' in cfg:
        _amber_react(cfg, args)
        return

    # -- single-residue path ---------------------------------------------------
    from craft import run_amber_pipeline, get_resname

    res_cfg = cfg.get('residue', {})
    amb_cfg = cfg.get('amber',   {}) or {}
    cap_cfg = cfg.get('cap',     {}) or {}

    charge     = args.charge if args.charge is not None else res_cfg.get('charge', 0)
    position   = res_cfg.get('position', 'middle')
    forcefield = amb_cfg.get('forcefield', 'ff14SB')

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
        atom_type  = amb_cfg.get('atom_type', 'amber'),
        forcefield = forcefield,
        capped_pdb = capped_pdb,
        position   = position,
    )


# ---------------------------------------------------------------------------
# craft-slurm   (SLURM script generator)
# ---------------------------------------------------------------------------

def _slurm_react(cfg, config_path):
    """SLURM script generator body for a two-residue reaction."""
    from pathlib import Path
    from craft.slurm import write_react_slurm

    res1_cfg = cfg['residue1']
    res2_cfg = cfg['residue2']

    charge1      = res1_cfg.get('charge', 0)
    charge2      = res2_cfg.get('charge', 0)
    total_charge = charge1 + charge2

    resname1, resname2 = _resolve_react_resnames(cfg)

    proj_root = Path.cwd().resolve()
    workdir   = proj_root / f"{resname1}_{resname2}"
    workdir.mkdir(parents=True, exist_ok=True)

    base   = f"{resname1}_{resname2}_react"
    output = str(workdir / f"{base}_craft.sh")

    write_react_slurm(cfg, output, proj_root, workdir,
                      resname1, resname2, total_charge,
                      config_path=config_path)


def slurm():
    """Generate a SLURM batch script for the full pipeline."""
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description="Generate a SLURM batch script for the full pipeline. "
                    "Detects single-residue or reaction mode from the config file.",
    )
    parser.add_argument('--config', default='config.yaml',
                        help='Config file (default: config.yaml)')
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    cfg = _load_config(args.config)

    if 'residue1' in cfg:
        _slurm_react(cfg, config_path)
        return

    # -- single-residue path ---------------------------------------------------
    from craft import write_slurm, get_resname

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
