#!/usr/bin/env python3
"""
make_gaussian_input.py  —  standalone entry point (delegates to parameterize.gaussian)
Generates the B3LYP geometry-optimisation .com from a capped PDB.

Usage:
    python make_gaussian_input.py MEO_capped.pdb
    python make_gaussian_input.py MEO_capped.pdb -c 1 -m 1
    python make_gaussian_input.py MEO_capped.pdb MEO_opt.com
"""

import sys
import argparse
from pathlib import Path
from parameterize.gaussian import write_com, NPROC_DEFAULT, MEM_DEFAULT, ROUTE_DEFAULT

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert a capped PDB to a Gaussian optimisation input (.com)"
    )
    parser.add_argument('pdb',            help='Capped PDB file')
    parser.add_argument('com', nargs='?', help='Output .com  (default: <base>_opt.com)')
    parser.add_argument('-c', '--charge', type=int, default=0,
                        help='Net molecular charge (default: 0)')
    parser.add_argument('-m', '--mult',   type=int, default=1,
                        help='Spin multiplicity  (default: 1)')
    parser.add_argument('-n', '--nproc',  type=int, default=NPROC_DEFAULT)
    parser.add_argument('--mem',          default=MEM_DEFAULT)
    args = parser.parse_args()

    base = Path(args.pdb).stem.replace('_capped', '')
    com  = args.com or f"{base}_opt.com"

    write_com(args.pdb, com, args.charge, args.mult,
              nproc=args.nproc, mem=args.mem, route=ROUTE_DEFAULT)
