#!/usr/bin/env python3
"""
cap_termini.py  —  standalone entry point (delegates to craft.cap)

Usage:
    python cap_termini.py MEO.pdb
    python cap_termini.py MEO.pdb MEO_capped.pdb
"""

import sys
from pathlib import Path
from craft import cap

if __name__ == '__main__':
    inp = sys.argv[1] if len(sys.argv) > 1 else 'KME3.pdb'
    out = sys.argv[2] if len(sys.argv) > 2 else Path(inp).stem + '_capped.pdb'
    cap(inp, out)
