#!/usr/bin/env python3
"""
make_hf_input.py  —  Phase 2 intermediate step
Extract the final optimised geometry from a Gaussian opt log and write
the HF/6-31G(d) single-point .com for ESP/RESP charge fitting.

Usage:
    python make_hf_input.py MEO_opt.log
    python make_hf_input.py MEO_opt.log MEO_hf.com   # explicit output name
    python make_hf_input.py MEO_opt.log -c 1          # for cationic residue
"""

import sys
import argparse
import yaml
from pathlib import Path
from craft import parse_opt_log, write_hf_com, get_resname
from craft.gaussian import NPROC_DEFAULT, MEM_DEFAULT, HF_ROUTE_DEFAULT


def load_config(path='config.yaml'):
    try:
        with open(path) as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        return {}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Gaussian opt log → HF/6-31G(d) single-point .com"
    )
    parser.add_argument('log',             help='Gaussian optimisation log (e.g. MEO_opt.log)')
    parser.add_argument('com', nargs='?',  help='Output .com  (default: <base>_hf.com)')
    parser.add_argument('-c', '--charge',  type=int, default=None,
                        help='Net molecular charge (default: from config.yaml or 0)')
    parser.add_argument('-m', '--mult',    type=int, default=None,
                        help='Spin multiplicity  (default: from config.yaml or 1)')
    parser.add_argument('-n', '--nproc',   type=int, default=None)
    parser.add_argument('--mem',           default=None)
    args = parser.parse_args()

    # Fall back to config.yaml for any unspecified values
    cfg     = load_config()
    res_cfg = cfg.get('residue', {})
    g_cfg   = cfg.get('gaussian_hf', {}) or cfg.get('gaussian', {}) or {}

    charge   = args.charge if args.charge is not None else res_cfg.get('charge', 0)
    mult     = args.mult   if args.mult   is not None else res_cfg.get('multiplicity', 1)
    nproc    = args.nproc  if args.nproc  is not None else g_cfg.get('nproc', NPROC_DEFAULT)
    mem      = args.mem    if args.mem    is not None else g_cfg.get('mem', MEM_DEFAULT)
    position = res_cfg.get('position', 'middle')

    log_path = Path(args.log)
    workdir  = log_path.parent   # write output alongside the log

    if res_cfg.get('input_pdb'):
        resname = get_resname(res_cfg['input_pdb'])
    else:
        resname = log_path.stem.replace('_opt', '').replace('_hf', '')
        for pos in ('_cterm', '_nterm'):
            resname = resname.replace(pos, '')

    suffix   = '' if position == 'middle' else f'_{position}'
    base     = f"{resname}{suffix}"
    com_path = args.com or str(workdir / f"{base}_hf.com")

    print(f"Parsing optimised geometry from {args.log} ...")
    atoms_xyz = parse_opt_log(args.log)
    print(f"  Found {len(atoms_xyz)} atoms in final Standard orientation block")

    print(f"Writing HF/6-31G(d) single-point input ...")
    write_hf_com(atoms_xyz, com_path, base,
                 charge=charge, mult=mult, nproc=nproc, mem=mem)
