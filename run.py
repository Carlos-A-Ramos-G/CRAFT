#!/usr/bin/env python3
"""
run.py  --  Phase 1: pre-Gaussian preparation
Reads config.yaml and produces everything needed before submitting to HPC.
All outputs are written to <resname>/<position>/ (created automatically).

  Step 1 -- Cap the residue PDB according to its terminal position
  Step 2 -- Write the geometry-optimisation Gaussian input (B3LYP/6-31G*)
  Step 3 -- Write resp.in  (RESP charge-fitting control file)
  Step 4 -- Write resp.qin (initial RESP charges from AMBER ff14SB)
  Step 5 -- Write <base>.mc (prepgen main-chain definition)

Terminal positions (residue.position in config.yaml):
  middle : ACE---residue---NME  ->  <resname>/middle/
  cterm  : ACE---residue        ->  <resname>/cterm/
  nterm  : residue---NME        ->  <resname>/nterm/

Usage:
    python run.py                     # uses config.yaml in current directory
    python run.py path/to/config.yaml
"""

import sys
import yaml
from pathlib import Path
from craft import cap, get_resname, write_com, write_resp_in, write_resp_qin, write_mc
from craft.gaussian import NPROC_DEFAULT, MEM_DEFAULT, ROUTE_DEFAULT, HF_ROUTE_DEFAULT


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else 'config.yaml'
    cfg      = load_config(cfg_path)

    res       = cfg['residue']
    input_pdb = res['input_pdb']
    charge    = res.get('charge', 0)
    mult      = res.get('multiplicity', 1)
    position  = res.get('position', 'middle')
    resname   = get_resname(input_pdb)

    if position not in ('middle', 'cterm', 'nterm'):
        sys.exit(f"Error: residue.position must be 'middle', 'cterm', or 'nterm'; "
                 f"got {position!r}")

    suffix = '' if position == 'middle' else f'_{position}'
    base   = f"{resname}{suffix}"

    # -- Create output directory -----------------------------------------------
    workdir = Path(resname) / position
    workdir.mkdir(parents=True, exist_ok=True)

    # -- Step 1: cap -----------------------------------------------------------
    cap_cfg    = cfg.get('cap') or {}
    capped_pdb = cap_cfg.get('output_pdb') or str(workdir / f"{base}_capped.pdb")

    print("=" * 60)
    print(f"Step 1 -- Cap termini  [{position}]")
    print("=" * 60)
    cap(input_pdb, capped_pdb, position=position)

    # -- Step 2: geometry-optimisation Gaussian input --------------------------
    gauss_cfg = cfg.get('gaussian_opt') or cfg.get('gaussian') or {}
    opt_com   = gauss_cfg.get('output_com') or str(workdir / f"{base}_opt.com")

    print()
    print("=" * 60)
    print("Step 2 -- Geometry-optimisation Gaussian input")
    print("=" * 60)
    write_com(
        capped_pdb, opt_com,
        charge=charge,
        mult=mult,
        nproc=gauss_cfg.get('nproc', NPROC_DEFAULT),
        mem=gauss_cfg.get('mem',   MEM_DEFAULT),
        route=gauss_cfg.get('route', ROUTE_DEFAULT),
    )

    # -- Steps 3-5: RESP and prepgen inputs -----------------------------------
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
    print(f"       python make_hf_input.py {workdir}/{base}_opt.log")
    print(f"  3. Submit {workdir}/{base}_hf.com to HPC")
    print(f"  4. Copy {base}_hf.log into {workdir}/, then:")
    print(f"       python amber_pipeline.py {workdir}/{base}_hf.log")
    print(f"\n  Residue '{resname}'  position '{position}'  read from {input_pdb}")


if __name__ == '__main__':
    main()
