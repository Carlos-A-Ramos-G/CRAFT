"""
craft.slurm
Generate a single SLURM batch script that runs the entire parameterization
pipeline end-to-end without manual intervention between phases:

  Phase 1  – python run.py
  Phase 2a – g16 < <base>_opt.com > <base>_opt.log
  Phase 2b – python make_hf_input.py <base>_opt.log
  Phase 2c – g16 < <base>_hf.com  > <base>_hf.log
  Phase 3  – python amber_pipeline.py <base>_hf.log
"""

from pathlib import Path
from .cap import get_resname


_TEMPLATE = """\
#!/bin/bash
{sbatch_directives}
# ── Environment ───────────────────────────────────────────────────────────────
{module_lines}
export GAUSS_SCRDIR=$(pwd)
{conda_block}
# ── Phase 1: cap termini, generate Gaussian inputs and RESP files ─────────────
echo "[$(date '+%H:%M:%S')] Phase 1 — run.py"
python run.py
[ $? -ne 0 ] && echo "ERROR in Phase 1" && exit 1

# ── Phase 2a: geometry optimisation ──────────────────────────────────────────
echo "[$(date '+%H:%M:%S')] Phase 2a — geometry optimisation ({opt_com})"
g16 < {opt_com} > {opt_log}
[ $? -ne 0 ] && echo "ERROR in Phase 2a (Gaussian opt)" && exit 1

# ── Phase 2b: build HF/ESP input from optimised geometry ─────────────────────
echo "[$(date '+%H:%M:%S')] Phase 2b — make_hf_input.py"
python make_hf_input.py {opt_log}
[ $? -ne 0 ] && echo "ERROR in Phase 2b (make_hf_input)" && exit 1

# ── Phase 2c: HF/6-31G(d) single-point for ESP/RESP ─────────────────────────
echo "[$(date '+%H:%M:%S')] Phase 2c — HF single-point ({hf_com})"
g16 < {hf_com} > {hf_log}
[ $? -ne 0 ] && echo "ERROR in Phase 2c (Gaussian HF)" && exit 1

# ── Phase 3: AMBER parameterization ──────────────────────────────────────────
echo "[$(date '+%H:%M:%S')] Phase 3 — amber_pipeline.py"
python amber_pipeline.py {hf_log}
[ $? -ne 0 ] && echo "ERROR in Phase 3 (AMBER pipeline)" && exit 1

echo "[$(date '+%H:%M:%S')] Done."
"""


def write_slurm(cfg, output):
    """
    Generate the SLURM batch script from config dict.

    Parameters
    ----------
    cfg    : dict — full parsed config.yaml
    output : str | Path — path of the generated script (e.g. 'KME3_craft.sh')
    """
    res_cfg  = cfg.get('residue', {})
    g_opt    = cfg.get('gaussian_opt', {}) or {}
    g_hf     = cfg.get('gaussian_hf',  {}) or {}
    sl       = cfg.get('slurm', {}) or {}

    input_pdb = res_cfg.get('input_pdb', '')
    resname   = get_resname(input_pdb) if input_pdb else Path(input_pdb or 'residue.pdb').stem
    job_name  = sl.get('job_name') or f"{resname}_craft"

    opt_com  = g_opt.get('output_com') or f"{resname}_opt.com"
    opt_log  = Path(opt_com).stem + '.log'
    hf_com   = g_hf.get('output_com')  or f"{resname}_hf.com"
    hf_log   = Path(hf_com).stem  + '.log'

    # ── #SBATCH directives ────────────────────────────────────────────────────
    # Always emit job-name first (may be derived, not literally in config).
    # Then iterate every key in the slurm section in config order, skipping
    # non-directive keys and blank values — so users can freely add/remove
    # fields (mem, time, nodes, account, …) without touching this code.
    _SKIP = {'job_name', 'modules', 'conda_env'}

    directives = [f'#SBATCH --job-name={job_name}']
    for key, val in sl.items():
        if key in _SKIP or val is None or val == '':
            continue
        directives.append(f'#SBATCH --{key.replace("_", "-")}={val}')

    # ── module load lines ─────────────────────────────────────────────────────
    modules = sl.get('modules') or []
    module_lines = (f"module load {' '.join(modules)}"
                    if modules else '# (no modules configured)')

    # ── conda activation ──────────────────────────────────────────────────────
    conda_env = sl.get('conda_env', '').strip()
    if conda_env:
        conda_block = (
            f'source "$(conda info --base)/etc/profile.d/conda.sh"\n'
            f'conda activate {conda_env}\n'
        )
    else:
        conda_block = ''

    script = _TEMPLATE.format(
        sbatch_directives='\n'.join(directives),
        module_lines=module_lines,
        conda_block=conda_block,
        opt_com=opt_com,
        opt_log=opt_log,
        hf_com=hf_com,
        hf_log=hf_log,
    )

    Path(output).write_text(script)
    print(f"SLURM script : {output}")
    print(f"  Residue    : {resname}")
    for d in directives:
        print(f"  {d.lstrip('#')}")
    if modules:
        print(f"  module load : {', '.join(modules)}")
    if conda_env:
        print(f"  conda env   : {conda_env}")
    print(f"  Workflow   : {opt_com} → {opt_log} → {hf_com} → {hf_log}")
    print(f"\nSubmit with:  sbatch {output}")
