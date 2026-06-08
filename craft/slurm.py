"""
craft.slurm
Generate a single SLURM batch script that runs the entire parameterization
pipeline end-to-end without manual intervention between phases:

  Phase 1  – craft-run
  Phase 2a – g16 < <base>_opt.com > <base>_opt.log
  Phase 2b – craft-hf-input <base>_opt.log
  Phase 2c – g16 < <base>_hf.com  > <base>_hf.log
  Phase 3  – craft-amber <base>_hf.log
"""

from pathlib import Path
from .cap import get_resname


_TEMPLATE = """\
#!/bin/bash
{sbatch_directives}
# -- Paths (baked in at generation time) --------------------------------------
PROJ_ROOT="{proj_root}"
CRAFT_WD="{workdir}"
CRAFT_CONFIG="{config_path}"

# -- Environment ---------------------------------------------------------------
{module_lines}
export GAUSS_SCRDIR="$CRAFT_WD"
{conda_block}
cd "$PROJ_ROOT"

# -- Phase 1: cap termini, generate Gaussian inputs and RESP files -------------
echo "[$(date '+%H:%M:%S')] Phase 1 -- craft-run"
craft-run --config "$CRAFT_CONFIG"
[ $? -ne 0 ] && echo "ERROR in Phase 1" && exit 1

# -- Phase 2a: geometry optimisation ------------------------------------------
echo "[$(date '+%H:%M:%S')] Phase 2a -- geometry optimisation ({opt_com})"
g16 < "$CRAFT_WD/{opt_com}" > "$CRAFT_WD/{opt_log}"
[ $? -ne 0 ] && echo "ERROR in Phase 2a (Gaussian opt)" && exit 1

# -- Phase 2b: build HF/ESP input from optimised geometry ---------------------
echo "[$(date '+%H:%M:%S')] Phase 2b -- craft-hf-input"
craft-hf-input "$CRAFT_WD/{opt_log}" --config "$CRAFT_CONFIG"
[ $? -ne 0 ] && echo "ERROR in Phase 2b (craft-hf-input)" && exit 1

# -- Phase 2c: HF/6-31G(d) single-point for ESP/RESP -------------------------
echo "[$(date '+%H:%M:%S')] Phase 2c -- HF single-point ({hf_com})"
g16 < "$CRAFT_WD/{hf_com}" > "$CRAFT_WD/{hf_log}"
[ $? -ne 0 ] && echo "ERROR in Phase 2c (Gaussian HF)" && exit 1

# -- Phase 3: AMBER parameterization ------------------------------------------
echo "[$(date '+%H:%M:%S')] Phase 3 -- craft-amber"
craft-amber "$CRAFT_WD/{hf_log}" --config "$CRAFT_CONFIG"
[ $? -ne 0 ] && echo "ERROR in Phase 3 (craft-amber)" && exit 1

echo "[$(date '+%H:%M:%S')] Done. Output files in $CRAFT_WD"
"""


def write_slurm(cfg, output, proj_root, workdir, position='middle', config_path=None):
    """
    Generate the SLURM batch script from config dict.

    Parameters
    ----------
    cfg         : dict       -- full parsed config.yaml
    output      : str | Path -- path of the generated script
    proj_root   : str | Path -- absolute path to the working directory
    workdir     : str | Path -- absolute path to the variant output directory
                                (e.g. <resname>/<position>/)
    position    : str        -- 'middle', 'cterm', or 'nterm'
    config_path : str | Path -- absolute path to the config file; baked into
                                the script so craft-* commands find it regardless
                                of the file name (default: proj_root/config.yaml)
    """
    res_cfg  = cfg.get('residue', {})
    g_opt    = cfg.get('gaussian_opt', {}) or {}
    g_hf     = cfg.get('gaussian_hf',  {}) or {}
    sl       = cfg.get('slurm', {}) or {}

    input_pdb = res_cfg.get('input_pdb', '')
    resname   = get_resname(input_pdb) if input_pdb else 'residue'
    prefix    = {'middle': '', 'cterm': 'C', 'nterm': 'N'}[position]
    base      = f"{prefix}{resname}"
    job_name  = sl.get('job_name') or f"{base}_craft"

    opt_com  = Path(g_opt.get('output_com') or f"{base}_opt.com").name
    hf_com   = Path(g_hf.get('output_com')  or f"{base}_hf.com").name
    opt_log  = Path(opt_com).stem + '.log'
    hf_log   = Path(hf_com).stem  + '.log'

    # -- #SBATCH directives ----------------------------------------------------
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

    # -- module load lines -----------------------------------------------------
    modules = sl.get('modules') or []
    module_lines = (f"module load {' '.join(modules)}"
                    if modules else '# (no modules configured)')

    # -- conda activation ------------------------------------------------------
    conda_env = (sl.get('conda_env') or '').strip()
    if conda_env:
        conda_block = (
            f'source "$(conda info --base)/etc/profile.d/conda.sh"\n'
            f'conda activate {conda_env}\n'
        )
    else:
        conda_block = ''

    if config_path is None:
        config_path = Path(proj_root) / 'config.yaml'

    script = _TEMPLATE.format(
        sbatch_directives='\n'.join(directives),
        module_lines=module_lines,
        conda_block=conda_block,
        proj_root=str(proj_root),
        workdir=str(workdir),
        config_path=str(config_path),
        opt_com=opt_com,
        opt_log=opt_log,
        hf_com=hf_com,
        hf_log=hf_log,
    )

    output = Path(output)
    output.write_text(script)
    print(f"SLURM script : {output}")
    print(f"  Residue    : {resname}")
    for d in directives:
        print(f"  {d.lstrip('#')}")
    if modules:
        print(f"  module load : {', '.join(modules)}")
    if conda_env:
        print(f"  conda env   : {conda_env}")
    print(f"  Workflow   : {opt_com} → {opt_log} → {hf_com} → {hf_log}")
    print(f"\nTo submit:")
    print(f"  cd {workdir}")
    print(f"  sbatch {output.name}")


# ---------------------------------------------------------------------------
# Bond workflow SLURM template
# ---------------------------------------------------------------------------

_BOND_TEMPLATE = """\
#!/bin/bash
{sbatch_directives}
# -- Paths (baked in at generation time) --------------------------------------
PROJ_ROOT="{proj_root}"
CRAFT_WD="{workdir}"
CRAFT_CONFIG="{config_path}"

# -- Environment ---------------------------------------------------------------
{module_lines}
export GAUSS_SCRDIR="$CRAFT_WD"
{conda_block}
cd "$PROJ_ROOT"

# -- Phase 1: prepare combined model and generate Gaussian/RESP inputs ---------
echo "[$(date '+%H:%M:%S')] Phase 1 -- craft-run"
craft-run --config "$CRAFT_CONFIG"
[ $? -ne 0 ] && echo "ERROR in Phase 1" && exit 1

# -- Phase 2a: frozen-backbone geometry optimisation --------------------------
echo "[$(date '+%H:%M:%S')] Phase 2a -- geometry optimisation ({opt_com})"
g16 < "$CRAFT_WD/{opt_com}" > "$CRAFT_WD/{opt_log}"
[ $? -ne 0 ] && echo "ERROR in Phase 2a (Gaussian opt)" && exit 1

# -- Phase 2b: build HF/ESP input from optimised geometry ---------------------
echo "[$(date '+%H:%M:%S')] Phase 2b -- craft-hf-input"
craft-hf-input "$CRAFT_WD/{opt_log}" --charge {total_charge} --config "$CRAFT_CONFIG"
[ $? -ne 0 ] && echo "ERROR in Phase 2b (craft-hf-input)" && exit 1

# -- Phase 2c: HF/6-31G(d) single-point for ESP/RESP -------------------------
echo "[$(date '+%H:%M:%S')] Phase 2c -- HF single-point ({hf_com})"
g16 < "$CRAFT_WD/{hf_com}" > "$CRAFT_WD/{hf_log}"
[ $? -ne 0 ] && echo "ERROR in Phase 2c (Gaussian HF)" && exit 1

# -- Phase 3: AMBER parameterization ------------------------------------------
echo "[$(date '+%H:%M:%S')] Phase 3 -- craft-amber"
craft-amber "$CRAFT_WD/{hf_log}" --config "$CRAFT_CONFIG"
[ $? -ne 0 ] && echo "ERROR in Phase 3 (craft-amber)" && exit 1

echo "[$(date '+%H:%M:%S')] Done. Output files in $CRAFT_WD"
"""


def write_bond_slurm(cfg, output, proj_root, workdir,
                     resname1, resname2, total_charge, config_path=None):
    """
    Generate the SLURM batch script for the two-residue covalent bond workflow.

    Parameters
    ----------
    cfg          : dict       -- full parsed config yaml
    output       : str | Path -- path of the generated script
    proj_root    : str | Path -- absolute path to the project root
    workdir      : str | Path -- absolute path to <resname1>_<resname2>/
    resname1     : str        -- first residue name
    resname2     : str        -- second residue name
    total_charge : int        -- sum of residue1.charge + residue2.charge
    config_path  : str | Path -- absolute config path baked into the script
    """
    sl      = cfg.get('slurm', {}) or {}
    g_opt   = cfg.get('gaussian_opt', {}) or {}
    g_hf    = cfg.get('gaussian_hf',  {}) or {}

    base     = f"{resname1}_{resname2}_bond"
    job_name = sl.get('job_name') or f"{base}_craft"

    opt_com = Path(g_opt.get('output_com') or f"{base}_opt.com").name
    hf_com  = Path(g_hf.get('output_com')  or f"{base}_hf.com").name
    opt_log = Path(opt_com).stem + '.log'
    hf_log  = Path(hf_com).stem  + '.log'

    _SKIP = {'job_name', 'modules', 'conda_env'}
    directives = [f'#SBATCH --job-name={job_name}']
    for key, val in sl.items():
        if key in _SKIP or val is None or val == '':
            continue
        directives.append(f'#SBATCH --{key.replace("_", "-")}={val}')

    modules      = sl.get('modules') or []
    module_lines = (f"module load {' '.join(modules)}"
                    if modules else '# (no modules configured)')

    conda_env = (sl.get('conda_env') or '').strip()
    if conda_env:
        conda_block = (
            f'source "$(conda info --base)/etc/profile.d/conda.sh"\n'
            f'conda activate {conda_env}\n'
        )
    else:
        conda_block = ''

    if config_path is None:
        config_path = Path(proj_root) / 'config.yaml'

    script = _BOND_TEMPLATE.format(
        sbatch_directives='\n'.join(directives),
        module_lines=module_lines,
        conda_block=conda_block,
        proj_root=str(proj_root),
        workdir=str(workdir),
        config_path=str(config_path),
        opt_com=opt_com,
        opt_log=opt_log,
        hf_com=hf_com,
        hf_log=hf_log,
        total_charge=total_charge,
    )

    output = Path(output)
    output.write_text(script)
    print(f"SLURM script : {output}")
    print(f"  Residues   : {resname1} + {resname2}  (charge {total_charge:+d})")
    for d in directives:
        print(f"  {d.lstrip('#')}")
    if modules:
        print(f"  module load : {', '.join(modules)}")
    if conda_env:
        print(f"  conda env   : {conda_env}")
    print(f"  Workflow   : {opt_com} → {opt_log} → {hf_com} → {hf_log}")
    print(f"\nTo submit:")
    print(f"  cd {workdir}")
    print(f"  sbatch {output.name}")
