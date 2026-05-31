"""
Generate a single SLURM batch script that runs the entire parameterization
pipeline end-to-end. The script is written into <resname>/<position>/
alongside all other pipeline outputs.

Usage:
    python make_slurm.py              # reads config.yaml
    python make_slurm.py my.yaml      # alternate config file
    python make_slurm.py my.yaml out.sh  # explicit output path

After generation, submit with:
    cd <resname>/<position>
    sbatch <base>_craft.sh
"""

import sys
import yaml
from pathlib import Path
from craft import write_slurm, get_resname


def main():
    config_path = sys.argv[1] if len(sys.argv) > 1 else 'config.yaml'
    cfg = yaml.safe_load(Path(config_path).read_text())

    proj_root = Path.cwd().resolve()
    res_cfg   = cfg.get('residue', {})
    input_pdb = res_cfg.get('input_pdb', '')
    position  = res_cfg.get('position', 'middle')
    resname   = get_resname(input_pdb) if input_pdb else 'residue'
    suffix    = '' if position == 'middle' else f'_{position}'
    base      = f"{resname}{suffix}"

    workdir = proj_root / resname / position
    workdir.mkdir(parents=True, exist_ok=True)

    output = sys.argv[2] if len(sys.argv) > 2 else str(workdir / f"{base}_craft.sh")

    write_slurm(cfg, output, proj_root, workdir, position)


if __name__ == '__main__':
    main()
