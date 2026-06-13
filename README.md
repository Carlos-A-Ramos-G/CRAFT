# CRAFT

**Custom Residue AMBER Forcefield Toolkit**

CRAFT automates the parameterization of non-standard amino acid residues for AMBER molecular dynamics simulations, following the standard RESP charge-fitting protocol.

---

## Overview

Given a single-residue PDB file, CRAFT produces the `.prepin` topology and `.frcmod` force-field parameter files needed to simulate the residue in AMBER, using AMBER's standard backbone charges and RESP-fitted sidechain charges. Both ff14SB and ff19SB are supported.

CRAFT also handles **two-residue covalent parameterization** (disulfide bonds, isopeptide bonds, NOS bonds, acylation of Ser/Cys, and similar cross-links). RESP charges are fitted jointly on the bonded system so mutual polarization across the new bond is captured, and `antechamber` sees the full bonded geometry to assign correct atom types at the bond interface. Separate `.prepin` files are produced for each residue; the cross-link bond is declared in tleap.

If you managed to parameterize any residue following my <a href="https://carlosramosg.com/amber-custom-residue-parameterization" target="_blank">previous tutorial</a>, CRAFT will make the process feel like a walk in the park. I tried to automate and simplify everything as much as possible.


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
espgen -> resp -> antechamber -> resolve DU atom types -> prepgen -> parmchk2
post-process frcmod: remove ATTN placeholders; warn on penalty score > 100
  |
  v
<resname>/<position>/<base>.prepin              -- residue topology
<resname>/<position>/<base>_gaff.frcmod        -- GAFF missing parameters
<resname>/<position>/<base>_{forcefield}.frcmod -- FF-specific missing parameters
```

Where `<base>` is `<resname>` for `middle`, `C<resname>` for `cterm`, and `N<resname>` for `nterm`.

As you will see below, the entire protocol can be run in an end-to-end fashion on a cluster machine if Gaussian is available.

---

## Installation (Recommended approach)

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

This will register the `craft-*` commands in the active conda environment. Run from any directory once installed.

```bash
pip install .
```


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

All commands must be run from the **project root** — the directory that contains `config.yaml`. CRAFT creates all outputs under `<resname>/<position>/` relative to that directory. Running from a subdirectory will cause `craft-run` and `craft-slurm` to fail (config not found), and will cause `craft-hf-input` and `craft-amber` to write outputs to the wrong location.

Expected layout after Phase 1:

```
project_root/          ← run all commands from here
├── config.yaml
├── KME3/
│   └── KME3.pdb       ← input_pdb
└── KME3/middle/       ← created by craft-run
    ├── KME3_capped.pdb
    ├── KME3_opt.com
    ├── resp.in
    ├── resp.qin
    └── KME3.mc
```

*If CRAFT is installed as a package (`pip install .`):*

```bash
# Phase 1 -- local; creates <resname>/<position>/
craft-run --config config.yaml

# Phase 2a -- submit geometry optimisation to HPC
#   Input:  <resname>/<position>/<base>_opt.com
#   Output: <resname>/<position>/<base>_opt.log  (copy back from HPC)

# Phase 2b -- local, after opt log arrives
craft-hf-input <resname>/<position>/<base>_opt.log --config config.yaml

# Phase 2c -- submit HF/ESP single-point to HPC
#   Input:  <resname>/<position>/<base>_hf.com
#   Output: <resname>/<position>/<base>_hf.log  (copy back from HPC)

# Phase 3 -- local, after HF log arrives
craft-amber <resname>/<position>/<base>_hf.log --config config.yaml
```

*Automated (single SLURM job):*

```bash
craft-slurm --config config.yaml        # generates <resname>/<position>/<base>_craft.sh
cd <resname>/<position>
sbatch <base>_craft.sh
```

*If running directly from source (no install):*

The root-level scripts (`run.py`, `make_hf_input.py`, `amber_pipeline.py`, `make_slurm.py`) are thin wrappers that delegate directly to `craft.cli.*`. They accept the same flags as the `craft-*` commands and produce identical output. The `craft` package must still be importable — run from the repository root or add it to `PYTHONPATH`.

```bash
python run.py --config config.yaml
python make_hf_input.py <resname>/<position>/<base>_opt.log --config config.yaml
python amber_pipeline.py <resname>/<position>/<base>_hf.log --config config.yaml
# or
python make_slurm.py --config config.yaml
cd <resname>/<position> && sbatch <base>_craft.sh
```

The SLURM script has absolute paths baked in at generation time, so it runs correctly regardless of where SLURM sets the working directory.

---

## Two-residue covalent parameterization

`craft-run` and `craft-amber` handle pairs of residues that form a covalent bond through their side chains — disulfide bonds, isopeptide bonds, NOS bonds (Lys–Cys), acylation of Ser or Cys, and similar cross-links — when the config file contains `residue1`, `residue2`, and `bond` keys. Both residues are parameterized jointly: RESP charges are fitted on the assembled bonded system so that mutual polarization across the new bond is captured, and `antechamber` sees the full bonded geometry to assign correct atom types at the reactive interface. Separate `.prepin` files are produced for each residue; the cross-link bond itself is declared in tleap.
When `bond.combined_pdb` is set, `atom1`, `atom2`, and `bond_length` are ignored. `input_pdb` fields become optional (used only for residue name resolution if present).

For the full annotated config template see `config_bond.yaml` in the repository root.

### Folder structure

All outputs are written under `<resname1>/<resname2>/`. Per-residue files are placed in dedicated subdirectories; shared QM and AMBER files stay at the pair level.

```
<resname1>/<resname2>/
  <r1>_<r2>_bond_combined.pdb     ← assembled model compound (4 cap groups)
  <r1>_<r2>_bond_opt.com          ← frozen-backbone Gaussian opt input
  resp.in, resp.qin, rename_map.json
  <r1>_<r2>_bond.ac
  <r1>_<r2>_bond_gaff.frcmod      ← shared for both residues
  <r1>_<r2>_bond_ff14SB.frcmod    ← shared for both residues
  <resname1>/
    <resname1>_capped.pdb
    <resname1>.mc
    <resname1>.prepin
  <resname2>/
    <resname2>_capped.pdb
    <resname2>.mc
    <resname2>.prepin
```

### Workflow

```
[Standard path]                        [Pre-assembled path]
residue1.pdb  residue2.pdb             combined.pdb (bonded, already capped)
      |               |                        |
   cap()           cap()              rename conflicting atoms (block 2)
      |               |               renumber resSeq to canonical 1–6
      +-- assemble --+                         |
   angle + torsion geometry                    |
   (RDKit UFF or numpy torsion scan)           |
              |                               |
              +------------- merge, unique atom names -----------+
                                              |
                    combined.pdb  (2 to 4 ACE/NME cap groups, unique names)
                                              |
                     write frozen-backbone opt .com
                       (caps + backbone N/CA/C/O fixed; side chains free)
                     write resp.in / resp.qin / two .mc files
                                              |
                              Phase 2a (HPC)
                     Gaussian B3LYP/6-31G* frozen-backbone optimisation
                                              |
                              Phase 2b (local)
                     craft-hf-input <opt.log> --charge <total> --config config.yaml
                                              |
                              Phase 2c (HPC)
                     Gaussian HF/6-31G(d) single-point (ESP)
                                              |
                              Phase 3 (local)
                     espgen → resp → antechamber (combined) → resolve DU atom types → prepgen × 2 → parmchk2
                     post-process frcmod: remove ATTN placeholders; warn on penalty score > 100
```

### Commands

```bash
# Phase 1 -- local
craft-run --config config.yaml

# Phase 2b -- after opt log arrives from HPC
craft-hf-input <r1>_<r2>/<r1>_<r2>_bond_opt.log --charge <total> --config config.yaml

# Phase 3 -- after HF log arrives from HPC
craft-amber <r1>_<r2>/<r1>_<r2>_bond_hf.log --config config.yaml
```

### tleap

Each residue's `.prepin` describes it in isolation. The cross-link bond is declared in tleap following the same structure used while describing disulfide bonds:

```
loadAmberPrep <resname1>.prepin
loadAmberPrep <resname2>.prepin
loadAmberParams <r1>_<r2>_bond_ff14SB.frcmod
mol = loadPdb system.pdb
bond mol.X.<atom1> mol.Y.<atom2>
saveAmberParm mol mol.prmtop mol.inpcrd
```

### Pre-assembled geometry

If you already have a PDB with both side chains bonded **and already capped with ACE/NME** — from your own QM workflow or molecular modelling — set `bond.combined_pdb` in the config. The PDB must contain the full capped model: ACE and/or NME groups on each residue according to the requested `position` values (e.g. `middle`+`middle` requires 6 residue groups total). CRAFT renames any conflicting atoms in the second residue block, renumbers resSeq to its canonical scheme, and proceeds directly to writing the frozen-backbone opt input. The CRAFT assembly step is skipped; the Gaussian geometry optimisation is not. Any resSeq numbering is accepted.

```yaml
bond:
  atom1: SG
  atom2: CG
  combined_pdb: my_bonded_structure.pdb
```

The residue names in the combined PDB must match those inferred from `residue1.input_pdb` and `residue2.input_pdb`.

---

## Commands

| Command | Phase | Description |
|---------|-------|-------------|
| `craft-check` | — | Verify all required tools and packages are available |
| `craft-run` | 1 | Cap termini, generate all pre-Gaussian inputs |
| `craft-hf-input` | 2b | Extract optimised geometry, write HF `.com` |
| `craft-amber` | 3 | Run espgen → resp → antechamber → resolve DU types → prepgen → parmchk2; clean frcmod |
| `craft-slurm` | — | Generate a single SLURM script for the full pipeline |

All commands read from `config.yaml` by default; pass `--config <path>` for a different file. Bond mode is detected automatically when the config contains `residue1`/`residue2`/`bond` keys. CLI flags override config values — run any command with `--help` for the full list.

If invoking from outside the project root, pass an absolute config path:

```bash
craft-hf-input KME3/middle/KME3_opt.log --config /abs/path/to/config.yaml
```

Equivalent convenience scripts (`run.py`, `make_hf_input.py`, `amber_pipeline.py`, `make_slurm.py`) are thin wrappers that accept the same flags as the corresponding `craft-*` commands. The `craft` package must be importable (run from the repository root, or set `PYTHONPATH`).

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
  forcefield: ff14SB          # ff14SB | ff19SB
  workdir:                    # leave blank -> derived from HF log path
  atom_type_overrides:        # optional: manually assign AMBER/GAFF types for atoms
    # ATOMNAME: TYPE          #   that antechamber could not resolve (DU type)
    #   e.g.  SG: S

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

For the bond config template see `config_bond.yaml` in the repository root — it is fully annotated.

---

## Package structure

```
craft/
  __init__.py    -- public API
  cap.py         -- ACE/NME capping, PDB I/O, geometry utilities
  gaussian.py    -- Gaussian .com writers, opt log parser
  resp.py        -- resp.in / resp.qin generation, RESP equivalence detection
  mc.py          -- prepgen main-chain (.mc) file writer
  amber.py       -- antechamber, prepgen, parmchk2 runner; atom name remapping; frcmod post-processing
  bond.py        -- two-residue covalent parameterization; angle-correct assembly geometry
  slurm.py       -- SLURM batch script generator
  cli.py         -- craft-* command entry points
  check.py       -- environment checker (craft-check)
```

---

## Examples

Ready-to-use input PDBs and config files are included at the repository root:

```
single_AA/             -- single-residue example (trimethyllysine KM3, all three positions)
two_bonded_AA/         -- two-residue covalent example (CYA–ASA thioether bond)
```

Each directory contains the PDB(s) and `config.yaml` needed to run the full pipeline from `craft-run` through `craft-amber`.

---

## RESP charge protocol

Fixed and free atoms depend on the terminal position:

| Position | Fixed (ff14SB/ff19SB values) | Free (RESP-fitted) |
|---|---|---|
| `middle` | ACE + backbone N, H, CA, HA, C, O + NME | sidechain |
| `cterm` | ACE + backbone N, H, CA, HA | sidechain + C-terminal C, O, OXT |
| `nterm` | NME + backbone CA, HA, C, O | sidechain + N-terminal N, H atoms |

ff14SB and ff19SB share identical backbone charges (N−0.4157, H+0.2719, CA+0.0337, HA+0.0823, C+0.5973, O−0.5679), so the same fixed constraints apply to both. The HA charge is adjusted by −0.0015 to make the six backbone atoms sum to exactly zero, preventing charge transfer artifacts at QM/MM boundaries. For glycine the HA charge is split equally between HA2 and HA3; for α-substituted residues with no alpha-H it is folded onto CA. The `forcefield` setting only affects which parm file is passed to `parmchk2` for the FF-specific `.frcmod`.

Symmetry-equivalent sidechain atoms (e.g. the three NZ-methyl carbons of trimethyllysine) are constrained equal using RDKit canonical Morgan ranking. Without RDKit, only H atoms bonded to the same heavy atom are constrained (a warning is printed).

---

## Reference

Bayly, C. I.; Cieplak, P.; Cornell, W. D.; Kollman, P. A. *A well-behaved electrostatic potential based method using charge restraints for deriving atomic charges: the RESP model.* J. Phys. Chem. **1993**, 97, 10269–10280.
