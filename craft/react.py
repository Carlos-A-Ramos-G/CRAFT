"""
craft.react
Two-residue side-chain reaction parameterization.

Workflow
--------
Phase 1  – craft-run --config config.yaml
             cap each residue independently (ACE/NME on backbone)
             assemble combined model compound with unique atom names
             write frozen-backbone Gaussian opt input
             write resp.in / resp.qin / two .mc files
Phase 2a – g16 < <r1>_<r2>_react_opt.com > <r1>_<r2>_react_opt.log
Phase 2b – craft-hf-input <r1>_<r2>_react_opt.log --charge <total> --config config.yaml
Phase 2c – g16 < <r1>_<r2>_react_hf.com  > <r1>_<r2>_react_hf.log
Phase 3  – craft-amber <r1>_<r2>_react_hf.log --config config.yaml

ResSeq conventions in the combined PDB
---------------------------------------
  1 = ACE(res1)   2 = res1   3 = NME(res1)
  4 = ACE(res2)   5 = res2   6 = NME(res2)

Atom naming
-----------
All atom names in the combined PDB are globally unique. Residue2 atoms that
conflict with residue1 names get a numeric suffix (N→N2, CA→CA2, etc.).
The rename_map {renamed→original} is saved to rename_map.json so that
craft-amber can restore original names in residue2's .prepin after prepgen.

antechamber is run on the combined HF log (not per-residue) so that the reactive
atom at the bond interface receives the correct AMBER atom type for its bonded
state. prepgen is then run twice with per-residue .mc OMIT lists.
"""

import re
import numpy as np
from pathlib import Path

from .cap import parse_pdb, pdb_line, _elem
from .gaussian import NPROC_DEFAULT, MEM_DEFAULT, _write_gjf
from .resp import (
    _find_equiv, ATOMIC_NUMBERS,
    ACE_CHARGES, NME_CHARGES,
    BACKBONE_N_CHARGES, BACKBONE_CA_CHARGE,
    BACKBONE_HA_TOTAL, BACKBONE_C_CHARGES,
)


_REACT_ROUTE_DEFAULT = "#P b3lyp/6-31g* opt(modredundant)"
_CAP_RESSEQS         = {1, 3, 4, 6}
_RES_RESSEQS         = {2, 5}
_BACKBONE_NAMES      = {'N', 'CA', 'C', 'O'}


def _pos(a):
    return np.array([a['x'], a['y'], a['z']])


def _base_name(name):
    """Strip trailing digits: 'CA2'→'CA', 'N2'→'N', 'CG1'→'CG'.
    Safe for backbone identification because no standard sidechain atom
    (ND, CG, OE, …) degrades to a backbone name (N, CA, C, O)."""
    return name.rstrip('0123456789')


# -- PDB assembly --------------------------------------------------------------

def _make_unique_names(atoms1, atoms2):
    """
    Return (atoms1_copy, atoms2_renamed, rename_map).
    atoms1 keeps its names; conflicting names in atoms2 get a numeric suffix.
    rename_map : {new_name → original_name} for every atom renamed in atoms2.
    """
    used = {a['name'] for a in atoms1}
    renamed    = []
    rename_map = {}
    for a in atoms2:
        name = a['name']
        if name in used:
            k = 2
            while f"{name}{k}" in used:
                k += 1
            name = f"{name}{k}"
            rename_map[name] = a['name']
        used.add(name)
        renamed.append({**a, 'name': name})
    return list(atoms1), renamed, rename_map


def _split_combined_pdb(combined_pdb, resname1, resname2, out1, out2):
    """
    Split a user-supplied combined PDB into two single-residue PDB files.

    Atoms are identified by resName.  Both output files are written with
    resSeq=1 so that cap() treats them as fresh single-residue inputs.
    Raises ValueError if either resName is not found in the combined PDB.
    """
    atoms = parse_pdb(combined_pdb)
    for resname, out in [(resname1, out1), (resname2, out2)]:
        res_atoms = [a for a in atoms if a['resName'] == resname]
        if not res_atoms:
            raise ValueError(
                f"No atoms with resName {resname!r} found in {combined_pdb}. "
                f"Residue names in file: "
                f"{sorted({a['resName'] for a in atoms})}")
        lines = []
        for i, a in enumerate(res_atoms, start=1):
            lines.append(pdb_line(i, a['name'], a['resName'], 1,
                                  (a['x'], a['y'], a['z']),
                                  a.get('occupancy', 1.0), a.get('tempFactor', 0.0)))
        lines.append(f'TER   {len(res_atoms) + 1:5d}\n')
        lines.append('END\n')
        Path(out).write_text(''.join(lines))


def assemble_react_pdb(capped_pdb1, capped_pdb2, atom1, atom2,
                       bond_length=1.5, output=None, skip_reposition=False):
    """
    Combine two capped single-residue PDBs into one model compound.

    By default (skip_reposition=False) residue2 is translated so that atom2
    sits at bond_length from atom1, along the extension of the CA1→atom1
    vector.  When skip_reposition=True the coordinates of both PDBs are used
    as-is; atom1, atom2, and bond_length are ignored.  Use skip_reposition
    when the user has supplied a pre-assembled combined structure and the bond
    geometry is already correct — only capping and atom-name uniqueness are
    needed.

    Parameters
    ----------
    capped_pdb1    : str | Path -- capped PDB for residue1 (output of cap())
    capped_pdb2    : str | Path -- capped PDB for residue2 (output of cap())
    atom1          : str        -- reactive atom name in residue1 (ignored when
                                   skip_reposition=True)
    atom2          : str        -- reactive atom name in residue2 (ignored when
                                   skip_reposition=True)
    bond_length    : float      -- initial bond length in Å; ignored when
                                   skip_reposition=True
    output         : str | Path | None -- path for combined PDB; None = no file
    skip_reposition: bool       -- if True, keep coordinates from both input
                                   PDBs unchanged (pre-assembled geometry path)

    Returns
    -------
    (atoms1_out, atoms2_out, rename_map)
    rename_map : {renamed_name → original_name} for residue2 atoms that were
                 renamed to ensure global uniqueness.
    """
    atoms1 = parse_pdb(capped_pdb1)
    atoms2 = parse_pdb(capped_pdb2)

    if skip_reposition:
        atoms2_shifted = list(atoms2)
    else:
        ca1   = next(a for a in atoms1 if a['resSeq'] == 2 and a['name'] == 'CA')
        bond1 = next(a for a in atoms1 if a['resSeq'] == 2 and a['name'] == atom1)
        bond2 = next(a for a in atoms2 if a['resSeq'] == 2 and a['name'] == atom2)

        direction = _pos(bond1) - _pos(ca1)
        direction /= np.linalg.norm(direction)
        delta = (_pos(bond1) + bond_length * direction) - _pos(bond2)

        atoms2_shifted = [{**a, 'x': a['x'] + delta[0],
                                'y': a['y'] + delta[1],
                                'z': a['z'] + delta[2]} for a in atoms2]

    atoms1_out, atoms2_out, rename_map = _make_unique_names(atoms1, atoms2_shifted)

    remap = {1: 4, 2: 5, 3: 6}
    for a in atoms2_out:
        a['resSeq'] = remap[a['resSeq']]

    if output is not None:
        ser   = 1
        lines = []
        for a in atoms1_out + atoms2_out:
            lines.append(pdb_line(ser, a['name'], a['resName'], a['resSeq'],
                                  (a['x'], a['y'], a['z']),
                                  a.get('occupancy', 1.0), a.get('tempFactor', 0.0)))
            ser += 1
        lines.append(f'TER   {ser:5d}\n')
        lines.append('END\n')
        Path(output).write_text(''.join(lines))
        print(f"Combined   : {output}  ({ser - 1} atoms, {len(rename_map)} renamed)")

    return atoms1_out, atoms2_out, rename_map


# -- Gaussian frozen-backbone opt input ----------------------------------------

def _frozen_indices(atoms):
    """
    1-based indices of atoms frozen during backbone-frozen optimisation.
    Frozen: all cap atoms (resSeq 1,3,4,6) and backbone heavy atoms of both
    residues (resSeq 2,5): N/CA/C/O identified via _base_name so renamed atoms
    (N2, CA2 …) are correctly frozen.
    """
    frozen = []
    for i, a in enumerate(atoms, start=1):
        if a['resSeq'] in _CAP_RESSEQS:
            frozen.append(i)
        elif a['resSeq'] in _RES_RESSEQS and _base_name(a['name']) in _BACKBONE_NAMES:
            frozen.append(i)
    return frozen


def write_react_com(combined_pdb, com_path, charge, mult,
                    nproc=NPROC_DEFAULT, mem=MEM_DEFAULT,
                    route=_REACT_ROUTE_DEFAULT):
    """
    Write a frozen-backbone geometry-optimisation .com for the combined system.
    Cap and backbone atoms are frozen via ModRedundant; only side chains move.
    """
    atoms = parse_pdb(combined_pdb)
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records in {combined_pdb}")

    if 'modredundant' not in route.lower():
        route = route.replace(' opt', ' opt(modredundant)', 1)

    base   = Path(combined_pdb).stem
    header = [f"%nprocshared={nproc}", f"%mem={mem}", route]
    atoms_xyz = [(_elem(a['name']), a['x'], a['y'], a['z']) for a in atoms]
    _write_gjf(com_path, header,
               f"{base}  reaction complex backbone-frozen optimisation",
               charge, mult, atoms_xyz)

    frozen = _frozen_indices(atoms)
    with open(com_path, 'a') as f:
        f.write('\n'.join(f"X {i} F" for i in frozen) + '\n\n')

    print(f"Input  : {combined_pdb}  ({len(atoms)} atoms, {len(frozen)} frozen)")
    print(f"Output : {com_path}")
    print(f"  Charge / mult : {charge} / {mult}  |  nproc / mem : {nproc} / {mem}")


# -- RESP classification -------------------------------------------------------

def _classify_react(atoms, position1='middle', position2='middle'):
    """
    Classify every atom in the combined two-residue system.

    resSeq 1,4 → 'ace'   resSeq 3,6 → 'nme'
    resSeq 2   → backbone/sidechain per position1
    resSeq 5   → backbone/sidechain per position2

    _base_name is used for all backbone name comparisons so that renamed atoms
    (N2, CA2, C2, O2 …) are correctly classified as backbone rather than sidechain.
    """
    def _classify_res(res_indices, position):
        amide_h = res_indices[1]
        ca = next((atoms[i] for i in res_indices
                   if _base_name(atoms[i]['name']) == 'CA'), None)
        alpha_h = set()
        if ca is not None:
            ca_pos = _pos(ca)
            for i in res_indices:
                if (_elem(atoms[i]['name']) == 'H'
                        and np.linalg.norm(_pos(atoms[i]) - ca_pos) < 1.15):
                    alpha_h.add(i)
        grp = {}
        for i in res_indices:
            bn = _base_name(atoms[i]['name'])
            if position == 'middle':
                if   bn == 'N':      grp[i] = 'bb_N'
                elif i == amide_h:   grp[i] = 'bb_H'
                elif bn == 'CA':     grp[i] = 'bb_CA'
                elif i in alpha_h:   grp[i] = 'bb_HA'
                elif bn == 'C':      grp[i] = 'bb_C'
                elif bn == 'O':      grp[i] = 'bb_O'
                else:                grp[i] = 'sidechain'
            elif position == 'cterm':
                if   bn == 'N':      grp[i] = 'bb_N'
                elif i == amide_h:   grp[i] = 'bb_H'
                elif bn == 'CA':     grp[i] = 'bb_CA'
                elif i in alpha_h:   grp[i] = 'bb_HA'
                else:                grp[i] = 'sidechain'
            else:  # nterm
                if   bn == 'CA':     grp[i] = 'bb_CA'
                elif i in alpha_h:   grp[i] = 'bb_HA'
                elif bn == 'C':      grp[i] = 'bb_C'
                elif bn == 'O':      grp[i] = 'bb_O'
                else:                grp[i] = 'sidechain'
        return grp

    res1_grp = _classify_res(
        [i for i, a in enumerate(atoms) if a['resSeq'] == 2], position1)
    res2_grp = _classify_res(
        [i for i, a in enumerate(atoms) if a['resSeq'] == 5], position2)

    groups = []
    for i, a in enumerate(atoms):
        if   a['resSeq'] in (1, 4): groups.append('ace')
        elif a['resSeq'] in (3, 6): groups.append('nme')
        elif a['resSeq'] == 2:      groups.append(res1_grp[i])
        else:                       groups.append(res2_grp[i])
    return groups


# -- RESP input files ----------------------------------------------------------

def write_react_resp_in(combined_pdb, charge, resname1, resname2, output,
                         position1='middle', position2='middle'):
    """Write RESP control file (resp.in) for the combined two-residue system."""
    atoms  = parse_pdb(combined_pdb)
    groups = _classify_react(atoms, position1, position2)
    equiv  = _find_equiv(atoms, groups)

    natoms = len(atoms)
    n_sc   = sum(1 for g in groups if g == 'sidechain')

    lines = [
        "capped-resp run #1",
        " &cntrl",
        " nmol=1,",
        " ihfree=1,",
        " qwt=0.0005,",
        " iqopt=2,",
        " /",
        "    1.00000",
        f"{resname1}-{resname2}",
        f"{charge:5d}{natoms:5d}",
    ]
    for g1, (a, grp) in enumerate(zip(atoms, groups), start=1):
        anum       = ATOMIC_NUMBERS[_elem(a['name'])]
        constraint = equiv.get(g1, 0) if grp == 'sidechain' else -1
        lines.append(f"{anum:5d}{constraint:5d}")
    lines += ["", ""]

    Path(output).write_text('\n'.join(lines) + '\n')
    print(f"  resp.in  : {output}  ({natoms} atoms, {n_sc} sidechain free)")


def write_react_resp_qin(combined_pdb, output, position1='middle', position2='middle'):
    """Write initial charges file (resp.qin) for the combined two-residue system."""
    atoms  = parse_pdb(combined_pdb)
    groups = _classify_react(atoms, position1, position2)

    def _ha_charges(resseq):
        n  = sum(1 for i, g in enumerate(groups)
                 if g == 'bb_HA' and atoms[i]['resSeq'] == resseq)
        ca = BACKBONE_CA_CHARGE + (BACKBONE_HA_TOTAL if n == 0 else 0.0)
        ha = BACKBONE_HA_TOTAL / n if n > 0 else 0.0
        return ca, ha

    ca1, ha1 = _ha_charges(2)
    ca2, ha2 = _ha_charges(5)

    ace_ctr = {1: 0, 4: 0}
    nme_ctr = {3: 0, 6: 0}

    charges = []
    for a, grp in zip(atoms, groups):
        rs = a['resSeq']
        if   grp == 'ace':   charges.append(ACE_CHARGES[ace_ctr[rs]]); ace_ctr[rs] += 1
        elif grp == 'nme':   charges.append(NME_CHARGES[nme_ctr[rs]]); nme_ctr[rs] += 1
        elif grp == 'bb_N':  charges.append(BACKBONE_N_CHARGES[0])
        elif grp == 'bb_H':  charges.append(BACKBONE_N_CHARGES[1])
        elif grp == 'bb_CA': charges.append(ca1 if rs == 2 else ca2)
        elif grp == 'bb_HA': charges.append(ha1 if rs == 2 else ha2)
        elif grp == 'bb_C':  charges.append(BACKBONE_C_CHARGES[0])
        elif grp == 'bb_O':  charges.append(BACKBONE_C_CHARGES[1])
        else:                charges.append(0.0)

    n_sc  = sum(1 for g in groups if g == 'sidechain')
    lines = [''.join(f"{q:10.6f}" for q in charges[i:i + 8])
             for i in range(0, len(charges), 8)]

    Path(output).write_text('\n'.join(lines) + '\n')
    print(f"  resp.qin : {output}  ({len(charges)} charges, {n_sc} sidechain free)")


# -- MC files ------------------------------------------------------------------

def write_react_mc(combined_pdb, target_resseq, charge, output, position='middle'):
    """
    Write prepgen .mc file for one residue in the combined system.

    target_resseq : 2 for residue1, 5 for residue2.

    HEAD/TAIL/MAIN_CHAIN use the actual (possibly renamed) atom names so they
    match the combined .ac file produced by antechamber.  _base_name is used to
    locate CA and C unambiguously even after renaming (CA2, C2).

    OMIT includes all cap atoms (resSeq 1,3,4,6) and all atoms of the other
    residue so prepgen extracts only the target residue from the combined .ac.

    After prepgen, call _restore_prepin_names() to convert renamed atoms back
    to their original names in residue2's .prepin.
    """
    atoms    = parse_pdb(combined_pdb)
    res_idx  = [i for i, a in enumerate(atoms) if a['resSeq'] == target_resseq]
    omit_idx = [i for i, a in enumerate(atoms) if a['resSeq'] != target_resseq]

    head_i = res_idx[0]
    ca_i   = next(i for i in res_idx if _base_name(atoms[i]['name']) == 'CA')
    tail_i = next(i for i in res_idx if _base_name(atoms[i]['name']) == 'C')

    lines = []
    if position in ('middle', 'cterm'):
        lines.append(f"HEAD_NAME {atoms[head_i]['name']}")
    if position in ('middle', 'nterm'):
        lines.append(f"TAIL_NAME {atoms[tail_i]['name']}")
    lines.append(f"MAIN_CHAIN {atoms[ca_i]['name']}")
    for i in omit_idx:
        lines.append(f"OMIT_NAME {atoms[i]['name']}")
    if position in ('middle', 'cterm'):
        lines.append("PRE_HEAD_TYPE C")
    if position in ('middle', 'nterm'):
        lines.append("POST_TAIL_TYPE N")
    lines.append(f"CHARGE {float(charge):.1f}")

    Path(output).write_text('\n'.join(lines) + '\n')
    print(f"  {Path(output).name:<20s}: HEAD={atoms[head_i]['name']}  "
          f"CA={atoms[ca_i]['name']}  TAIL={atoms[tail_i]['name']}  "
          f"omit {len(omit_idx)} atoms")


# -- prepin post-processing ----------------------------------------------------

def _restore_prepin_names(prepin_path, rename_map):
    """
    Reverse atom renaming in a prepgen .prepin file.

    rename_map : {renamed_name → original_name}
    Applied longest-key-first so 'CA2' is replaced before any shorter key that
    might be a substring. Word boundaries prevent 'N2' from matching inside
    longer tokens like 'ND2'.
    """
    content = Path(prepin_path).read_text()
    for new, orig in sorted(rename_map.items(), key=lambda kv: -len(kv[0])):
        content = re.sub(r'\b' + re.escape(new) + r'\b', orig, content)
    Path(prepin_path).write_text(content)
    print(f"  Restored {len(rename_map)} atom name(s) in {Path(prepin_path).name}")


# -- AMBER pipeline ------------------------------------------------------------

def run_react_amber_pipeline(hf_log, resname1, resname2, total_charge,
                              mc_file1, mc_file2, workdir='.',
                              atom_type='amber', forcefield='ff14SB',
                              combined_pdb=None, rename_map=None):
    """
    Post-Gaussian AMBER pipeline for the combined two-residue system.

    espgen → resp → antechamber (combined .ac) → prepgen × 2 → parmchk2

    antechamber runs on the combined HF log so that the reactive atom at the
    bond interface is typed correctly for its bonded state.  prepgen is run
    twice using per-residue .mc files with full OMIT lists.  rename_map is
    applied to residue2's .prepin to restore original atom names after prepgen.

    parmchk2 runs once on the combined .ac; the resulting frcmod covers missing
    parameters for both residues and is loaded once in tleap.
    """
    import shutil, os
    from .amber import _run, remap_ac_atom_names, PARM_FILES

    hf_log   = str(Path(hf_log).resolve())
    mc_file1 = str(Path(mc_file1).resolve())
    mc_file2 = str(Path(mc_file2).resolve())
    wd       = str(Path(workdir).resolve())

    base        = f"{resname1}_{resname2}_react"
    ac_file     = f"{base}.ac"
    gaff_frcmod = f"{base}_gaff.frcmod"
    ff_frcmod   = f"{base}_{forcefield}.frcmod" if forcefield else None

    if forcefield:
        parm_name = PARM_FILES.get(forcefield)
        if parm_name is None:
            raise ValueError(
                f"Unknown forcefield {forcefield!r}. Supported: {list(PARM_FILES)}")
        amberhome = os.environ.get('AMBERHOME', '')
        parm_file = (str(Path(amberhome) / 'dat/leap/parm' / parm_name)
                     if amberhome else parm_name)

    print("\n-- espgen ------------------------------------------------------------")
    _run(['espgen', '-i', hf_log, '-o', 'esp.dat'], cwd=wd)

    print("\n-- resp --------------------------------------------------------------")
    _run([
        'resp', '-O',
        '-i', 'resp.in',  '-o', 'resp.out',
        '-p', 'resp.pch', '-t', 'resp.chg',
        '-q', 'resp.qin', '-e', 'esp.dat',
    ], cwd=wd)

    print("\n-- antechamber -------------------------------------------------------")
    _run([
        'antechamber',
        '-fi', 'gout', '-i', hf_log,
        '-bk', resname1,
        '-fo', 'ac',   '-o', ac_file,
        '-c',  'rc',   '-cf', 'resp.chg',
        '-at', atom_type,
        '-nc', str(total_charge),
    ], cwd=wd)

    print("\n-- remap atom names --------------------------------------------------")
    if combined_pdb and Path(combined_pdb).exists():
        remap_ac_atom_names(Path(wd) / ac_file, combined_pdb)
    else:
        print("  combined PDB not found -- atom names in .ac not remapped.")

    for resname, mc_file in [(resname1, mc_file1), (resname2, mc_file2)]:
        sub_dir  = Path(wd) / resname
        sub_dir.mkdir(exist_ok=True)
        mc_local = Path(mc_file).name
        mc_dst   = sub_dir / mc_local
        if mc_dst.resolve() != Path(mc_file).resolve():
            shutil.copy(mc_file, mc_dst)
        prepin_rel = f"{resname}/{resname}.prepin"
        mc_rel     = f"{resname}/{mc_local}"
        print(f"\n-- prepgen ({resname}) {'-' * (51 - len(resname))}")
        _run(['prepgen', '-i', ac_file, '-o', prepin_rel,
              '-m', mc_rel, '-rn', resname], cwd=wd)

    if rename_map:
        print("\n-- restore atom names in residue2 prepin -------------------------")
        _restore_prepin_names(Path(wd) / resname2 / f"{resname2}.prepin", rename_map)

    print("\n-- parmchk2 (GAFF) ---------------------------------------------------")
    _run(['parmchk2', '-i', ac_file, '-f', 'ac', '-o', gaff_frcmod], cwd=wd)

    if forcefield:
        print(f"\n-- parmchk2 ({forcefield}) {'-' * (51 - len(forcefield))}")
        _run(['parmchk2', '-i', ac_file, '-f', 'ac', '-o', ff_frcmod,
              '-a', 'Y', '-p', parm_file], cwd=wd)

    output_files = [ac_file,
                    f"{resname1}/{resname1}.prepin",
                    f"{resname2}/{resname2}.prepin",
                    gaff_frcmod]
    if forcefield:
        output_files.append(ff_frcmod)
    print(f"\nDone. Output files in {wd}:")
    for fname in output_files:
        status = 'ok' if (Path(wd) / fname).exists() else 'MISSING'
        print(f"  [{status}]  {fname}")
