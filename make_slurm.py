"""
Generate a single SLURM batch script that runs the entire parameterization
pipeline end-to-end.

Usage:
    python make_slurm.py              # reads config.yaml, writes <base>_param.sh
    python make_slurm.py my.yaml      # alternate config file
    python make_slurm.py my.yaml out.sh  # explicit output path

Submit the generated script with:
    sbatch <base>_param.sh
"""

import sys
import yaml
from pathlib import Path
from parameterize import write_slurm, get_resname


def main():
    config_path = sys.argv[1] if len(sys.argv) > 1 else 'config.yaml'
    cfg = yaml.safe_load(Path(config_path).read_text())

    input_pdb = cfg.get('residue', {}).get('input_pdb', '')
    resname   = get_resname(input_pdb) if input_pdb else 'residue'
    output    = sys.argv[2] if len(sys.argv) > 2 else f"{resname}_param.sh"

    write_slurm(cfg, output)


if __name__ == '__main__':
    main()
