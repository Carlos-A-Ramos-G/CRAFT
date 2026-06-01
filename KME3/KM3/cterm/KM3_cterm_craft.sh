#!/bin/bash
#SBATCH --job-name=KM3_cterm_craft
#SBATCH --output=param_%j.out
#SBATCH --error=param_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
# -- Paths (baked in at generation time) --------------------------------------
PROJ_ROOT="/path/to/project"
CRAFT_WD="/path/to/project/KM3/cterm"
CRAFT_CONFIG="/path/to/project/config_cterm.yaml"

# -- Environment ---------------------------------------------------------------
module load apps/gaussian/g16 apps/amber/24
export GAUSS_SCRDIR="$CRAFT_WD"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate YOUR_ENV

cd "$PROJ_ROOT"

# -- Phase 1: cap termini, generate Gaussian inputs and RESP files -------------
echo "[$(date '+%H:%M:%S')] Phase 1 -- craft-run"
craft-run "$CRAFT_CONFIG"
[ $? -ne 0 ] && echo "ERROR in Phase 1" && exit 1

# -- Phase 2a: geometry optimisation ------------------------------------------
echo "[$(date '+%H:%M:%S')] Phase 2a -- geometry optimisation (KM3_cterm_opt.com)"
g16 < "$CRAFT_WD/KM3_cterm_opt.com" > "$CRAFT_WD/KM3_cterm_opt.log"
[ $? -ne 0 ] && echo "ERROR in Phase 2a (Gaussian opt)" && exit 1

# -- Phase 2b: build HF/ESP input from optimised geometry ---------------------
echo "[$(date '+%H:%M:%S')] Phase 2b -- craft-hf-input"
craft-hf-input "$CRAFT_WD/KM3_cterm_opt.log" --config "$CRAFT_CONFIG"
[ $? -ne 0 ] && echo "ERROR in Phase 2b (craft-hf-input)" && exit 1

# -- Phase 2c: HF/6-31G(d) single-point for ESP/RESP -------------------------
echo "[$(date '+%H:%M:%S')] Phase 2c -- HF single-point (KM3_cterm_hf.com)"
g16 < "$CRAFT_WD/KM3_cterm_hf.com" > "$CRAFT_WD/KM3_cterm_hf.log"
[ $? -ne 0 ] && echo "ERROR in Phase 2c (Gaussian HF)" && exit 1

# -- Phase 3: AMBER parameterization ------------------------------------------
echo "[$(date '+%H:%M:%S')] Phase 3 -- craft-amber"
craft-amber "$CRAFT_WD/KM3_cterm_hf.log" --config "$CRAFT_CONFIG"
[ $? -ne 0 ] && echo "ERROR in Phase 3 (craft-amber)" && exit 1

echo "[$(date '+%H:%M:%S')] Done. Output files in $CRAFT_WD"
