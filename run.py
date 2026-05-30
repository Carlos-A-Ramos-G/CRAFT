#!/usr/bin/env python3
"""
run.py
Run the full parameterization pipeline from config.yaml:
  Step 1 — Cap the residue PDB with ACE / NME termini
  Step 2 — Write the Gaussian input file for geometry optimisation + RESP

Usage:
    python run.py                    # uses config.yaml in current directory
    python run.py path/to/config.yaml
"""

import sys
import yaml
from pathlib import Path
from parameterize import cap, write_gjf
from parameterize.gaussian import NPROC_DEFAULT, MEM_DEFAULT, ROUTE_DEFAULT


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else 'config.yaml'
    cfg = load_config(cfg_path)

    res        = cfg['residue']
    input_pdb  = res['input_pdb']
    charge     = res.get('charge', 0)
    mult       = res.get('multiplicity', 1)

    # ── Step 1: cap ──────────────────────────────────────────────────────────
    cap_cfg    = cfg.get('cap') or {}
    capped_pdb = cap_cfg.get('output_pdb') or (Path(input_pdb).stem + '_capped.pdb')

    print("=" * 60)
    print("Step 1 — Cap termini")
    print("=" * 60)
    cap(input_pdb, capped_pdb)

    # ── Step 2: Gaussian input ───────────────────────────────────────────────
    gauss_cfg = cfg.get('gaussian') or {}
    base      = Path(capped_pdb).stem.replace('_capped', '')
    gjf       = gauss_cfg.get('output_gjf') or f"{base}_hf.gjf"

    print()
    print("=" * 60)
    print("Step 2 — Gaussian input")
    print("=" * 60)
    write_gjf(
        capped_pdb, gjf,
        charge=charge,
        mult=mult,
        nproc=gauss_cfg.get('nproc', NPROC_DEFAULT),
        mem=gauss_cfg.get('mem',   MEM_DEFAULT),
        route=gauss_cfg.get('route', ROUTE_DEFAULT),
    )


if __name__ == '__main__':
    main()
