#!/usr/bin/env python3
"""
run.py  —  Phase 1: pre-Gaussian preparation
Reads config.yaml and produces everything needed before submitting to HPC:
  Step 1 — Cap the residue PDB with ACE / NME termini
  Step 2 — Write the geometry-optimisation Gaussian input (B3LYP/6-31G*)
  Step 3 — Write resp.in  (RESP charge-fitting control file)
  Step 4 — Write resp.qin (initial RESP charges from AMBER ff14SB)
  Step 5 — Write <resname>.mc (prepgen main-chain definition)

After HPC jobs finish, run:
  python make_hf_input.py   <MEO_opt.log>   → MEO_hf.com  (submit to HPC)
  python amber_pipeline.py  <MEO_hf.log>    → .ac / .prepin / .frcmod files

Usage:
    python run.py                     # uses config.yaml in current directory
    python run.py path/to/config.yaml
"""

import sys
import yaml
from pathlib import Path
from parameterize import cap, get_resname, write_com, write_resp_in, write_resp_qin, write_mc
from parameterize.gaussian import NPROC_DEFAULT, MEM_DEFAULT, ROUTE_DEFAULT, HF_ROUTE_DEFAULT


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
    base      = Path(input_pdb).stem          # used for output filenames
    resname   = get_resname(input_pdb)        # residue name from PDB records

    # ── Step 1: cap ──────────────────────────────────────────────────────────
    cap_cfg    = cfg.get('cap') or {}
    capped_pdb = cap_cfg.get('output_pdb') or f"{resname}_capped.pdb"

    print("=" * 60)
    print("Step 1 — Cap termini")
    print("=" * 60)
    cap(input_pdb, capped_pdb)

    # ── Step 2: geometry-optimisation Gaussian input ──────────────────────────
    gauss_cfg = cfg.get('gaussian_opt') or cfg.get('gaussian') or {}
    opt_com   = gauss_cfg.get('output_com') or f"{resname}_opt.com"

    print()
    print("=" * 60)
    print("Step 2 — Geometry-optimisation Gaussian input")
    print("=" * 60)
    write_com(
        capped_pdb, opt_com,
        charge=charge,
        mult=mult,
        nproc=gauss_cfg.get('nproc', NPROC_DEFAULT),
        mem=gauss_cfg.get('mem',   MEM_DEFAULT),
        route=gauss_cfg.get('route', ROUTE_DEFAULT),
    )

    # ── Steps 3-5: RESP and prepgen inputs ───────────────────────────────────
    print()
    print("=" * 60)
    print("Steps 3-5 — RESP and prepgen inputs")
    print("=" * 60)
    write_resp_in( capped_pdb, charge, resname, 'resp.in')
    write_resp_qin(capped_pdb, 'resp.qin')
    write_mc(      capped_pdb, charge, f"{resname}.mc")

    print()
    print("Next steps:")
    print(f"  1. Submit {opt_com} to HPC (sbatch or equivalent)")
    print(f"  2. When {resname}_opt.log arrives:")
    print(f"       python make_hf_input.py {resname}_opt.log")
    print(f"  3. Submit {resname}_hf.com to HPC")
    print(f"  4. When {resname}_hf.log arrives:")
    print(f"       python amber_pipeline.py {resname}_hf.log")
    print(f"\n  Residue name '{resname}' read from {input_pdb}")


if __name__ == '__main__':
    main()
