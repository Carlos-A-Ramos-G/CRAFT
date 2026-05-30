# CRAFT

**Custom Residue AMBER Forcefield Toolkit**

CRAFT automates the parameterization of non-standard amino acid residues for AMBER molecular dynamics simulations, following the standard RESP charge-fitting protocol (Bayly et al., *J. Phys. Chem.* 97, 10269, 1993).

---

## Overview

Given a single-residue PDB file, CRAFT produces the `.prepin` topology and `.frcmod` force-field parameter files needed to simulate the residue in AMBER, using AMBER's standard ff14SB backbone charges and RESP-fitted sidechain charges.

```
raw PDB
  │
  ▼  Phase 1 (local)
cap termini (ACE/NME)
generate Gaussian opt input (.com)
generate RESP files (resp.in, resp.qin)
generate prepgen main-chain file (.mc)
  │
  ▼  Phase 2a (HPC)
Gaussian B3LYP/6-31G* geometry optimisation
  │
  ▼  Phase 2b (local)
extract optimised geometry → HF/6-31G(d) input (.com)
  │
  ▼  Phase 2c (HPC)
Gaussian HF/6-31G(d) single-point (ESP)
  │
  ▼  Phase 3 (local / HPC)
espgen → resp → antechamber → prepgen → parmchk2
  │
  ▼
<resname>.prepin   — residue topology
<resname>_gaff.frcmod   — GAFF missing parameters
<resname>_ff14SB.frcmod — ff14SB missing parameters
```

---

## Requirements

**Python packages**

```
numpy
pyyaml
rdkit        # optional — enables full symmetry-based RESP equivalence detection
```

Install with:

```bash
pip install -r requirements.txt
```

Without RDKit, CRAFT falls back to geometry-only equivalence detection (H atoms bonded to the same heavy atom are constrained equal, but equivalent heavy atoms such as symmetric methyl groups are not). A warning is printed at runtime.

**External tools** (must be available in `$PATH`)

| Tool | Source |
|------|--------|
| `g16` (or `g09`) | Gaussian 16 |
| `espgen`, `resp`, `antechamber`, `prepgen`, `parmchk2` | AmberTools / AMBER |

---

## Quick start

**1. Prepare your input**

Place a single-residue PDB in your working directory and fill in `config.yaml`:

```yaml
residue:
  input_pdb: KME3/KME3.pdb
  charge: +1
  multiplicity: 1
```

The PDB can be a free amino acid (zwitterionic), a residue cut from a PDB chain (bare N, single amide H, no OXT), or anything in between. CRAFT inspects the termini and handles all cases.

**2. Run the pipeline**

Phase 1 creates a `<resname>/` subdirectory and writes all outputs there, keeping the project root clean.

*Manual (step-by-step):*

```bash
# Phase 1 — local; creates <resname>/
python run.py

# Phase 2a — submit geometry optimisation to HPC
#   Input:  <resname>/<resname>_opt.com
#   Output: <resname>/<resname>_opt.log  (copy back from HPC)

# Phase 2b — local, after opt log arrives
python make_hf_input.py <resname>/<resname>_opt.log

# Phase 2c — submit HF/ESP single-point to HPC
#   Input:  <resname>/<resname>_hf.com
#   Output: <resname>/<resname>_hf.log  (copy back from HPC)

# Phase 3 — local, after HF log arrives
python amber_pipeline.py <resname>/<resname>_hf.log
```

*Automated (single SLURM job):*

```bash
python make_slurm.py     # generates <resname>/<resname>_craft.sh
cd <resname>
sbatch <resname>_craft.sh
```

The SLURM script has absolute paths baked in at generation time, so it can be submitted from any directory and runs correctly regardless of where SLURM sets the working directory.

---

## Scripts

| Script | Phase | Description |
|--------|-------|-------------|
| `run.py` | 1 | Cap termini, generate all pre-Gaussian inputs |
| `make_hf_input.py` | 2b | Extract optimised geometry, write HF `.com` |
| `amber_pipeline.py` | 3 | Run espgen → resp → antechamber → prepgen → parmchk2 |
| `cap_termini.py` | — | Standalone ACE/NME capping utility |
| `make_gaussian_input.py` | — | Standalone Gaussian opt `.com` writer |
| `make_slurm.py` | — | Generate a single SLURM script for the full pipeline |

---

## Configuration reference (`config.yaml`)

```yaml
residue:
  input_pdb: KME3/KME3.pdb   # raw single-residue PDB
  charge: +1                  # net molecular charge
  multiplicity: 1             # spin multiplicity (1 = singlet)

cap:
  output_pdb:                 # leave blank → <stem>_capped.pdb

gaussian_opt:
  nproc: 16
  mem: 512MB
  route: "#P b3lyp/6-31g* opt"
  output_com:                 # leave blank → <base>_opt.com

gaussian_hf:
  nproc: 16
  mem: 512MB
  route: "#p hf/6-31g(d) SCF=Tight Pop=MK IOp(6/33=2)"
  output_com:                 # leave blank → <base>_hf.com

amber:
  atom_type: amber            # amber | gaff | gaff2
  ff14sb_frcmod: true         # generate ff14SB frcmod (requires $AMBERHOME)
  workdir:                    # leave blank → <resname>/ (derived from HF log path)

slurm:
  job_name:                   # leave blank → <base>_craft
  output: param_%j.out
  error:  param_%j.err
  ntasks: 1
  cpus_per_task: 16
  account: YOUR_ACCOUNT
  partition: YOUR_PARTITION
  modules:
    - apps/gaussian/g16
    - apps/amber/24
  conda_env: your_env         # conda env with numpy/pyyaml; leave blank if not needed
```

---

## Package structure

```
craft/
  __init__.py    — public API
  cap.py         — ACE/NME capping, PDB I/O, geometry utilities
  gaussian.py    — Gaussian .com writers, opt log parser
  resp.py        — resp.in / resp.qin generation, RESP equivalence detection
  mc.py          — prepgen main-chain (.mc) file writer
  amber.py       — antechamber, prepgen, parmchk2 runner; atom name remapping
  slurm.py       — SLURM batch script generator
```

---

## RESP charge protocol

- **ACE and NME cap atoms** — fixed to AMBER ff14SB charges.
- **Backbone N, H (amide), C, O** — fixed to ff14SB values.
- **Sidechain atoms** — free to be fit by RESP; symmetry-equivalent atoms (e.g. the three NZ-methyl carbons of trimethyllysine) are constrained equal using RDKit canonical Morgan ranking.

---

## Roadmap

CRAFT is an early iteration, N-terminal and C-terminal forms of non-standard residues will be added in the future.

---

## Reference

Bayly, C. I.; Cieplak, P.; Cornell, W. D.; Kollman, P. A. *A well-behaved electrostatic potential based method using charge restraints for deriving atomic charges: the RESP model.* J. Phys. Chem. **1993**, 97, 10269–10280.
