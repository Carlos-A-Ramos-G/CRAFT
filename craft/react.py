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

# Ideal bond angles (degrees) at the reactive atom, keyed by element symbol.
# Used when assembling two residues from individual PDBs.
_IDEAL_ANGLE = {'C': 109.5, 'N': 109.5, 'O': 109.5, 'S': 103.0, 'P': 109.5}

# Van der Waals radii (Å) used for inter-residue clash scoring.
_VDW_RADII = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'H': 1.20}


def _pos(a):
    return np.array([a['x'], a['y'], a['z']], dtype=float)


def _base_name(name):
    """Strip trailing digits: 'CA2'→'CA', 'N2'→'N', 'CG1'→'CG'.
    Safe for backbone identification because no standard sidechain atom
    (ND, CG, OE, …) degrades to a backbone name (N, CA, C, O)."""
    return name.rstrip('0123456789')


def _vdw_r(atom_name):
    return _VDW_RADII.get(_elem(atom_name), 1.70)


def _rodrigues(v, axis, theta):
    """Rotate vector v around unit-vector axis by angle theta (radians)."""
    c, s = np.cos(theta), np.sin(theta)
    return v * c + np.cross(axis, v) * s + axis * np.dot(axis, v) * (1.0 - c)


def _sidechain_parent(atoms, reactive_name, resseq, max_bond=2.1):
    """Return the nearest heavy non-backbone atom in resseq bonded to reactive_name."""
    target = next(a for a in atoms if a['resSeq'] == resseq and a['name'] == reactive_name)
    tp = _pos(target)
    exclude = {'N', 'CA', 'C', 'O', reactive_name}
    best_d, best_a = np.inf, None
    for a in atoms:
        if a['resSeq'] != resseq or a['name'] in exclude or _elem(a['name']) == 'H':
            continue
        d = np.linalg.norm(_pos(a) - tp)
        if d < max_bond and d < best_d:
            best_d, best_a = d, a
    return best_a


def _angle_placement(atoms1, atoms2, a1_name, a2_name, bl):
    """
    Translate atoms2 rigidly so that atom a2_name sits at distance bl from
    a1_name and the parent1–atom1–atom2 angle equals the ideal bond angle for
    atom1's element.  A perpendicular direction is chosen using the CA1 atom
    to give a consistent (though arbitrary) initial torsion.

    Returns a new atoms2 list with updated coordinates.
    """
    bond1  = next(a for a in atoms1 if a['resSeq'] == 2 and a['name'] == a1_name)
    bond2  = next(a for a in atoms2 if a['resSeq'] == 2 and a['name'] == a2_name)
    ca1    = next(a for a in atoms1 if a['resSeq'] == 2 and a['name'] == 'CA')
    parent = _sidechain_parent(atoms1, a1_name, 2) or ca1

    b1 = _pos(bond1)
    u  = _pos(parent) - b1
    u /= np.linalg.norm(u)                     # parent → atom1 unit vector

    # perpendicular direction anchored to CA1 (Gram-Schmidt)
    ref  = _pos(ca1) - b1
    ref -= np.dot(ref, u) * u
    if np.linalg.norm(ref) < 1e-6:             # CA1 is collinear: use global z
        ref = np.array([0., 0., 1.]) - np.dot([0., 0., 1.], u) * u
    perp = ref / np.linalg.norm(ref)

    theta  = np.radians(_IDEAL_ANGLE.get(_elem(a1_name), 109.5))
    d_new  = np.cos(theta) * u + np.sin(theta) * perp  # atom1 → atom2 direction
    target = b1 + bl * d_new
    delta  = target - _pos(bond2)

    return [{**a, 'x': a['x'] + delta[0],
                  'y': a['y'] + delta[1],
                  'z': a['z'] + delta[2]} for a in atoms2]


def _torsion_scan(atoms1, atoms2_placed, a1_name, a2_name, n_steps=36):
    """
    Rotate atoms2_placed around the atom1→atom2 bond axis to find the
    orientation that minimises inter-residue Van der Waals clashes.

    Clash score: sum of cubic overlaps (r1 + r2 - distance)³ for all
    heavy-atom pairs between the two residues.

    Returns the best atoms2 list (lowest clash score).
    """
    b1   = _pos(next(a for a in atoms1        if a['resSeq'] == 2 and a['name'] == a1_name))
    b2   = _pos(next(a for a in atoms2_placed if a['resSeq'] == 2 and a['name'] == a2_name))
    axis = (b2 - b1) / np.linalg.norm(b2 - b1)

    # residue1 heavy atoms (vectorised for speed)
    r1_pos = np.array([[a['x'], a['y'], a['z']]
                       for a in atoms1 if _elem(a['name']) != 'H'])
    r1_rad = np.array([_vdw_r(a['name'])
                       for a in atoms1 if _elem(a['name']) != 'H'])

    best_score, best_atoms = np.inf, atoms2_placed

    for step in range(n_steps):
        theta   = 2.0 * np.pi * step / n_steps
        rotated = []
        for a in atoms2_placed:
            v   = np.array([a['x'], a['y'], a['z']]) - b1
            v_r = _rodrigues(v, axis, theta) + b1
            rotated.append({**a, 'x': v_r[0], 'y': v_r[1], 'z': v_r[2]})

        score = 0.0
        for a in rotated:
            if _elem(a['name']) == 'H':
                continue
            p2  = np.array([a['x'], a['y'], a['z']])
            ov  = (r1_rad + _vdw_r(a['name'])) - np.linalg.norm(r1_pos - p2, axis=1)
            score += float(np.sum(ov[ov > 0] ** 3))

        if score < best_score:
            best_score, best_atoms = score, rotated

    return best_atoms


def _geometry_torsion(atoms1, atoms2, a1_name, a2_name, bl):
    """Angle-corrected placement + torsion-scan clash minimisation (numpy only)."""
    placed = _angle_placement(atoms1, atoms2, a1_name, a2_name, bl)
    return _torsion_scan(atoms1, placed, a1_name, a2_name)


def _geometry_rdkit(capped_pdb1, capped_pdb2, atoms1, atoms2, a1_name, a2_name, bl):
    """
    RDKit UFF geometry: angle-corrected initial placement, then UFF minimisation
    with all cap atoms and backbone heavy atoms of both residues frozen.
    Only the sidechains near the new bond are free to relax.

    Raises ImportError if RDKit is not installed.
    Raises RuntimeError if UFF setup fails (falls back to torsion scan in caller).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    atoms2_init = _angle_placement(atoms1, atoms2, a1_name, a2_name, bl)

    mol1 = Chem.MolFromPDBFile(str(capped_pdb1), removeHs=False, sanitize=False)
    mol2 = Chem.MolFromPDBFile(str(capped_pdb2), removeHs=False, sanitize=False)
    if mol1 is None or mol2 is None:
        raise RuntimeError("RDKit could not parse the capped PDB files")

    # update mol2 conformer with initial-placement coordinates (match by atom name)
    init_xyz = {a['name']: (a['x'], a['y'], a['z']) for a in atoms2_init}
    conf2 = mol2.GetConformer()
    for i in range(mol2.GetNumAtoms()):
        info = mol2.GetAtomWithIdx(i).GetMonomerInfo()
        if info:
            xyz = init_xyz.get(info.GetName().strip())
            if xyz:
                conf2.SetAtomPosition(i, xyz)

    n1      = mol1.GetNumAtoms()
    combined = Chem.RWMol(Chem.CombineMols(mol1, mol2))

    def _find_by_name(name, start, stop):
        for i in range(start, stop):
            info = combined.GetAtomWithIdx(i).GetMonomerInfo()
            if info and info.GetName().strip() == name:
                return i
        return None

    idx1 = _find_by_name(a1_name, 0, n1)
    idx2 = _find_by_name(a2_name, n1, combined.GetNumAtoms())
    if idx1 is None or idx2 is None:
        raise RuntimeError(f"Reactive atoms {a1_name!r}/{a2_name!r} not found in RDKit mol")

    combined.AddBond(idx1, idx2, Chem.BondType.SINGLE)

    # partial sanitization — don't fail on unusual residue atom types
    Chem.SanitizeMol(combined,
        Chem.SanitizeFlags.SANITIZE_FINDRADICALS  |
        Chem.SanitizeFlags.SANITIZE_SETAROMATICITY |
        Chem.SanitizeFlags.SANITIZE_SETCONJUGATION |
        Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION |
        Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
        catchErrors=True)

    mol_f = combined.GetMol()
    ff = AllChem.UFFGetMoleculeForceField(mol_f)
    if ff is None:
        raise RuntimeError("UFF setup failed — atom types not recognized")

    # freeze caps (ACE, NME) and backbone heavy atoms (N, CA, C, O) of both residues
    _CAPS = {'ACE', 'NME'}
    _BB   = {'N', 'CA', 'C', 'O'}
    for i in range(mol_f.GetNumAtoms()):
        info = mol_f.GetAtomWithIdx(i).GetMonomerInfo()
        if info and (info.GetResidueName().strip() in _CAPS or
                     info.GetName().strip() in _BB):
            ff.AddFixedPoint(i)

    ff.Minimize(maxIts=500)

    # extract updated residue2 coordinates (match mol2 atom-name → rdkit index)
    conf = mol_f.GetConformer()
    rdkit_idx = {}
    for i in range(mol2.GetNumAtoms()):
        info = mol2.GetAtomWithIdx(i).GetMonomerInfo()
        if info:
            rdkit_idx[info.GetName().strip()] = n1 + i

    result = list(atoms2)
    for j, a in enumerate(atoms2):
        ri = rdkit_idx.get(a['name'])
        if ri is not None:
            p = conf.GetAtomPosition(ri)
            result[j] = {**a, 'x': p.x, 'y': p.y, 'z': p.z}
    return result


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
                       bond_length=None, output=None, skip_reposition=False):
    """
    Combine two capped single-residue PDBs into one model compound.

    By default (skip_reposition=False) residue2 is repositioned so that atom2
    bonds to atom1 at the correct distance AND bond angle.  If RDKit is
    available a UFF minimisation (backbone frozen) is run to relax the
    sidechain geometry; otherwise a numpy torsion scan finds the orientation
    with the lowest inter-residue VdW clash.  When skip_reposition=True the
    coordinates of both PDBs are used as-is (pre-assembled geometry path) and
    atom1, atom2, and bond_length are ignored.

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

    geom_method = None
    if skip_reposition:
        atoms2_shifted = list(atoms2)
    else:
        bl = bond_length if bond_length is not None else 1.5
        try:
            atoms2_shifted = _geometry_rdkit(capped_pdb1, capped_pdb2,
                                             atoms1, atoms2, atom1, atom2, bl)
            geom_method = 'RDKit UFF'
        except ImportError:
            atoms2_shifted = _geometry_torsion(atoms1, atoms2, atom1, atom2, bl)
            geom_method = 'torsion scan'
        except RuntimeError as exc:
            print(f"  Warning: RDKit geometry failed ({exc}) — using torsion scan")
            atoms2_shifted = _geometry_torsion(atoms1, atoms2, atom1, atom2, bl)
            geom_method = 'torsion scan'

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
        geom_note = f'  geometry: {geom_method}' if geom_method else ''
        print(f"Combined   : {output}  ({ser - 1} atoms, {len(rename_map)} renamed){geom_note}")

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
    from .amber import _run, remap_ac_atom_names, PARM_FILES, _postprocess_frcmod

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
    _postprocess_frcmod(Path(wd) / gaff_frcmod)

    if forcefield:
        print(f"\n-- parmchk2 ({forcefield}) {'-' * (51 - len(forcefield))}")
        _run(['parmchk2', '-i', ac_file, '-f', 'ac', '-o', ff_frcmod,
              '-a', 'Y', '-p', parm_file], cwd=wd)
        _postprocess_frcmod(Path(wd) / ff_frcmod)

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
