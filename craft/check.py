"""
craft.check
Environment checker -- verifies that all required external tools and
Python packages are present before submitting a pipeline run.
"""

import os
import shutil
import sys
from pathlib import Path


def check_env():
    """
    Check external tools and Python packages required by CRAFT.
    Exits with code 1 if any required component is missing.
    """
    all_ok = True

    # -- Python packages -------------------------------------------------------
    print("Python packages")
    print("-" * 40)
    for pkg, import_name in [('numpy', 'numpy'), ('pyyaml', 'yaml')]:
        try:
            __import__(import_name)
            print(f"  [ok]      {pkg}")
        except ImportError:
            print(f"  [MISSING] {pkg}  ->  pip install {pkg}")
            all_ok = False

    try:
        import rdkit
        print(f"  [ok]      rdkit  (full symmetry-based RESP equivalence)")
    except ImportError:
        print(f"  [--]      rdkit  not installed (optional; geometry-only"
              " equivalence will be used)")

    # -- Gaussian --------------------------------------------------------------
    print()
    print("Gaussian")
    print("-" * 40)
    found_gaussian = False
    for exe in ('g16', 'g09'):
        path = shutil.which(exe)
        if path:
            print(f"  [ok]      {exe} -> {path}")
            found_gaussian = True
            break
    if not found_gaussian:
        print("  [MISSING] g16 or g09 -- not found in $PATH")
        all_ok = False

    # -- AmberTools ------------------------------------------------------------
    print()
    print("AmberTools")
    print("-" * 40)
    for tool in ('espgen', 'resp', 'antechamber', 'prepgen', 'parmchk2'):
        path = shutil.which(tool)
        if path:
            print(f"  [ok]      {tool} -> {path}")
        else:
            print(f"  [MISSING] {tool} -- not found in $PATH")
            all_ok = False

    amberhome = os.environ.get('AMBERHOME', '')
    print()
    if amberhome:
        print(f"  $AMBERHOME = {amberhome}")
        parm10 = Path(amberhome) / 'dat/leap/parm/parm10.dat'
        if parm10.exists():
            print(f"  [ok]      parm10.dat found")
        else:
            print(f"  [warn]    parm10.dat not found at expected path:")
            print(f"            {parm10}")
            print(f"            ff14SB frcmod generation (parmchk2 -a Y) may fail")
    else:
        print("  $AMBERHOME is not set -- ff14SB frcmod generation will be skipped")

    # -- Summary ---------------------------------------------------------------
    print()
    print("=" * 40)
    if all_ok:
        print("All required tools found. CRAFT is ready to run.")
    else:
        print("One or more required tools are missing. See above.")
        sys.exit(1)
