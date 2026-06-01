# CRAFT

**Custom Residue AMBER Forcefield Toolkit**

CRAFT automates the parameterization of non-standard amino acid residues for AMBER molecular dynamics simulations, following the standard RESP charge-fitting protocol (Bayly et al., *J. Phys. Chem.* 97, 10269, 1993).

---

## Overview

Given a single-residue PDB file, CRAFT produces the `.prepin` topology and `.frcmod` force-field parameter files needed to simulate the residue in AMBER, using AMBER's standard ff14SB backbone charges and RESP-fitted sidechain charges.

Three terminal positions are supported:

| `position` | Capping | Use case |
|---|---|---|
| `middle` | ACE---residue---NME | Interior position in a peptide chain |
| `cterm`  | ACE---residue       | C-terminal residue |
| `nterm`  |       residue---NME | N-terminal residue |

Each variant is parameterized independently and written to its own subdirectory (`<resname>/<position>/`).

```
raw PDB
  |
  v  Phase 1 (local)
cap termini (ACE, NME, or both depending on position)
generate Gaussian opt input (.com)
generate RESP files (resp.in, resp.qin)
generate prepgen main-chain file (.mc)
  |
  v  Phase 2a (HPC)
Gaussian B3LYP/6-31G* geometry optimisation
  |
  v  Phase 2b (local)
extract optimised geometry -> HF/6-31G(d) input (.com)
  |
  v  Phase 2c (HPC)
Gaussian HF/6-31G(d) single-point (ESP)
  |
  v  Phase 3 (local / HPC)
espgen -> resp -> antechamber -> prepgen -> parmchk2
  |
  v
<resname>/<position>/<base>.prepin        -- residue topology
<resname>/<position>/<base>_gaff.frcmod  -- GAFF missing parameters
<resname>/<position>/<base>_ff14SB.frcmod -- ff14SB missing parameters
```

Where `<base>` is `<resname>` for `middle`, `<resname>_cterm` for `cterm`, and `<resname>_nterm` for `nterm`.

---

## Installation

**Python dependencies**

```bash
pip install -r requirements.txt
```

| Package | Required | Notes |
|---------|----------|-------|
| `numpy` | yes | |
| `pyyaml` | yes | |
| `rdkit` | no | Full symmetry-based RESP equivalence detection; falls back to geometry-only without it |

**Install CRAFT as a package**

```bash
pip install .
```

This registers the `craft-*` commands in the active conda environment. Run from any directory once installed.

**Verify the environment**

```bash
craft-check
```

This checks that all required tools (Gaussian, AmberTools) and Python packages are available before you submit a job.

**External tools** (must be in `$PATH`)

| Tool | Source |
|------|--------|
| `g16` (or `g09`) | Gaussian 16 |
| `espgen`, `resp`, `antechamber`, `prepgen`, `parmchk2` | AmberTools / AMBER |

---

## Quick start

**1. Prepare your input**

Place a single-residue PDB in your working directory and create a `config.yaml` (copy from the reference template):

```yaml
residue:
  input_pdb: KME3/KME3.pdb
  charge: +1
  multiplicity: 1
  position: middle        # middle | cterm | nterm
```

The PDB can be a free amino acid (zwitterionic), a residue cut from a PDB chain (bare N, single amide H, no OXT), or anything in between. CRAFT inspects the termini and handles all cases. For `cterm` and `nterm` variants, the free terminus must be correctly protonated in the input PDB.

**2. Run the pipeline**

*If CRAFT is installed as a package (`pip install .`):*

```bash
# Phase 1 -- local; creates <resname>/<position>/
craft-run

# Phase 2a -- submit geometry optimisation to HPC
#   Input:  <resname>/<position>/<base>_opt.com
#   Output: <resname>/<position>/<base>_opt.log  (copy back from HPC)

# Phase 2b -- local, after opt log arrives
craft-hf-input <resname>/<position>/<base>_opt.log

# Phase 2c -- submit HF/ESP single-point to HPC
#   Input:  <resname>/<position>/<base>_hf.com
#   Output: <resname>/<position>/<base>_hf.log  (copy back from HPC)

# Phase 3 -- local, after HF log arrives
craft-amber <resname>/<position>/<base>_hf.log
```

*Automated (single SLURM job):*

```bash
craft-slurm                              # generates <resname>/<position>/<base>_craft.sh
cd <resname>/<position>
sbatch <base>_craft.sh
```

*If running directly from source (no install):*

```bash
python run.py
python make_hf_input.py <resname>/<position>/<base>_opt.log
python amber_pipeline.py <resname>/<position>/<base>_hf.log
# or
python make_slurm.py
cd <resname>/<position> && sbatch <base>_craft.sh
```

The SLURM script has absolute paths baked in at generation time, so it runs correctly regardless of where SLURM sets the working directory.

---

## Commands

| Command | Phase | Description |
|---------|-------|-------------|
| `craft-check` | — | Verify all required tools and packages are available |
| `craft-run` | 1 | Cap termini, generate all pre-Gaussian inputs |
| `craft-hf-input` | 2b | Extract optimised geometry, write HF `.com` |
| `craft-amber` | 3 | Run espgen → resp → antechamber → prepgen → parmchk2 |
| `craft-slurm` | — | Generate a single SLURM script for the full pipeline |

Equivalent convenience scripts (`run.py`, `make_hf_input.py`, `amber_pipeline.py`, `make_slurm.py`) are also provided for users who prefer to run directly from source without installing.

---

## Configuration reference (`config.yaml`)

```yaml
residue:
  input_pdb: KME3/KME3.pdb   # raw single-residue PDB
  charge: +1                  # net molecular charge
  multiplicity: 1             # spin multiplicity (1 = singlet)
  position: middle            # middle | cterm | nterm

cap:
  output_pdb:                 # leave blank -> <base>_capped.pdb

gaussian_opt:
  nproc: 16
  mem: 512MB
  route: "#P b3lyp/6-31g* opt"
  output_com:                 # leave blank -> <base>_opt.com

gaussian_hf:
  nproc: 16
  mem: 512MB
  route: "#p hf/6-31g(d) SCF=Tight Pop=MK IOp(6/33=2)"
  output_com:                 # leave blank -> <base>_hf.com

amber:
  atom_type: amber            # amber | gaff | gaff2
  ff14sb_frcmod: true         # generate ff14SB frcmod (requires $AMBERHOME)
  workdir:                    # leave blank -> derived from HF log path

slurm:
  job_name:                   # leave blank -> <base>_craft
  output: param_%j.out
  error:  param_%j.err
  ntasks: 1
  cpus_per_task: 16
  account: YOUR_ACCOUNT
  partition: YOUR_PARTITION
  modules:
    - apps/gaussian/g16
    - apps/amber/24
  conda_env: your_env         # conda env where CRAFT is installed (pip install .)
                              # and that has numpy/pyyaml; leave blank if not needed
```

---

## Package structure

```
craft/
  __init__.py    -- public API
  cap.py         -- ACE/NME capping, PDB I/O, geometry utilities
  gaussian.py    -- Gaussian .com writers, opt log parser
  resp.py        -- resp.in / resp.qin generation, RESP equivalence detection
  mc.py          -- prepgen main-chain (.mc) file writer
  amber.py       -- antechamber, prepgen, parmchk2 runner; atom name remapping
  slurm.py       -- SLURM batch script generator
  cli.py         -- craft-* command entry points
  check.py       -- environment checker (craft-check)
```

---

## RESP charge protocol

Fixed and free atoms depend on the terminal position:

| Position | Fixed (ff14SB values) | Free (RESP-fitted) |
|---|---|---|
| `middle` | ACE + backbone N, H, C, O + NME | sidechain |
| `cterm` | ACE + backbone N, H | sidechain + C-terminal C, O, OXT |
| `nterm` | NME + backbone C, O | sidechain + N-terminal N, H atoms |

Symmetry-equivalent sidechain atoms (e.g. the three NZ-methyl carbons of trimethyllysine) are constrained equal using RDKit canonical Morgan ranking. Without RDKit, only H atoms bonded to the same heavy atom are constrained (a warning is printed).

---

## Reference

Bayly, C. I.; Cieplak, P.; Cornell, W. D.; Kollman, P. A. *A well-behaved electrostatic potential based method using charge restraints for deriving atomic charges: the RESP model.* J. Phys. Chem. **1993**, 97, 10269–10280.
