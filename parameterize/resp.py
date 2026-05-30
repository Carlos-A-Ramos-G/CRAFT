"""
parameterize.resp
Generate RESP charge-fitting input files (resp.in and resp.qin).

Constraints follow the standard AMBER RESP protocol:
  - ACE and NME cap atoms are FIXED to AMBER ff14SB charges.
  - Backbone N, H(amide), C, O of the residue are FIXED to ff14SB values.
  - Sidechain atoms are FREE to be fit; symmetry-equivalent atoms are
    constrained equal using RDKit canonical ranking (falls back to geometry-only
    H-on-same-heavy-atom detection if RDKit is unavailable).

Backbone atoms are identified by exact name ('N', 'C', 'O') and by position
(the amide H is always the second residue atom in cap.py output). This is
robust regardless of atom ordering within the residue.

Reference: Bayly et al., J. Phys. Chem. 97, 10269 (1993).
"""

import numpy as np
from pathlib import Path
from .cap import parse_pdb, _elem


# ── AMBER ff14SB charges for fixed atoms ──────────────────────────────────────
# Order matches the capped-PDB atom ordering produced by cap.py:
#   ACE: methyl-C  methyl-H1  methyl-H2  methyl-H3  carbonyl-C  carbonyl-O
#   Backbone N-terminus side: N  H(amide)
#   Backbone C-terminus side: C  O
#   NME: amide-N  amide-H  methyl-C  methyl-H1  methyl-H2  methyl-H3
# (Actual atom names in the PDB are unique per-molecule; charges map by position.)

ACE_CHARGES        = [-0.3662,  0.1123,  0.1123,  0.1123,  0.5972, -0.5679]
BACKBONE_N_CHARGES = [-0.4157,  0.2719]   # backbone N, amide H
BACKBONE_C_CHARGES = [ 0.5972, -0.5679]   # backbone C (carbonyl), O
NME_CHARGES        = [-0.4157,  0.2719, -0.1490,  0.0976,  0.0976,  0.0976]

ATOMIC_NUMBERS = {
    'H': 1, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'P': 15, 'S': 16, 'Cl': 17, 'Br': 35,
}


# ── Atom classification ───────────────────────────────────────────────────────

def _classify(atoms):
    """
    Classify every atom as one of:
      'ace' | 'bb_N' | 'bb_H' | 'sidechain' | 'bb_C' | 'bb_O' | 'nme'

    Backbone N  → exact name 'N' in residue.
    Backbone H  → second residue atom (cap.py always writes N then H).
    Backbone C  → exact name 'C' in residue (carbonyl C, not CA/CB/etc.).
    Backbone O  → exact name 'O' in residue (carbonyl O, not OG/OD1/etc.).
    Everything else in the residue is sidechain.
    """
    res_indices = [i for i, a in enumerate(atoms) if a['resSeq'] == 2]
    amide_h_global = res_indices[1]   # 0-based index in full atom list

    groups = []
    for i, a in enumerate(atoms):
        if a['resSeq'] == 1:
            groups.append('ace')
        elif a['resSeq'] == 3:
            groups.append('nme')
        else:
            if a['name'] == 'N':
                groups.append('bb_N')
            elif i == amide_h_global:
                groups.append('bb_H')
            elif a['name'] == 'C':
                groups.append('bb_C')
            elif a['name'] == 'O':
                groups.append('bb_O')
            else:
                groups.append('sidechain')

    return groups


# ── Equivalence detection ─────────────────────────────────────────────────────

def _pos(a):
    return np.array([a['x'], a['y'], a['z']])


def _find_equiv_rdkit(atoms, sc_indices):
    """
    Symmetry-equivalent sidechain atoms via RDKit canonical ranking.

    Builds the full molecular graph from covalent-radii distances (all single
    bonds), runs the Morgan algorithm with breakTies=False, and groups sidechain
    atoms that share the same rank.  This correctly identifies equivalent heavy
    atoms (e.g. the three NZ-methyl carbons in trimethyllysine) AND propagates
    the equivalence to their hydrogens.

    Returns {global_1based: ref_global_1based}.
    """
    from rdkit import Chem

    mol = Chem.RWMol()
    for a in atoms:
        mol.AddAtom(Chem.Atom(_elem(a['name'])))

    conf = Chem.Conformer(len(atoms))
    for i, a in enumerate(atoms):
        p = _pos(a)
        conf.SetAtomPosition(i, (float(p[0]), float(p[1]), float(p[2])))
    mol.AddConformer(conf, assignId=True)

    pt = Chem.GetPeriodicTable()
    n = mol.GetNumAtoms()
    for i in range(n):
        ri = pt.GetRcovalent(_elem(atoms[i]['name']))
        for j in range(i + 1, n):
            rj = pt.GetRcovalent(_elem(atoms[j]['name']))
            if np.linalg.norm(_pos(atoms[i]) - _pos(atoms[j])) < 1.3 * (ri + rj):
                mol.AddBond(i, j, Chem.BondType.SINGLE)

    # Skip valence check so unusual charge states (e.g. N+ without formal charge)
    # do not abort sanitization before hybridization / ring detection runs.
    Chem.SanitizeMol(
        mol,
        Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
    )

    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))

    rank_groups = {}
    for i in sc_indices:
        rank_groups.setdefault(ranks[i], []).append(i)

    equiv = {}
    for idxs in rank_groups.values():
        if len(idxs) < 2:
            continue
        idxs = sorted(idxs)
        ref = idxs[0] + 1          # 1-based global index
        for later in idxs[1:]:
            equiv[later + 1] = ref
    return equiv


def _find_equiv_geom(atoms, sc_indices):
    """
    Geometry-only fallback: constrain H atoms bonded to the same heavy atom.
    Returns {global_1based: ref_global_1based}.
    """
    sc_atoms = [atoms[i] for i in sc_indices]

    h_to_heavy = {}
    for li, a in enumerate(sc_atoms):
        if _elem(a['name']) != 'H':
            continue
        for lj, b in enumerate(sc_atoms):
            if li == lj or _elem(b['name']) not in ('C', 'N', 'O', 'S', 'P'):
                continue
            if np.linalg.norm(_pos(a) - _pos(b)) < 1.35:
                h_to_heavy[li] = lj
                break

    heavy_to_hs = {}
    for h_idx, heavy_idx in h_to_heavy.items():
        heavy_to_hs.setdefault(heavy_idx, []).append(h_idx)

    equiv = {}
    for _, hs in heavy_to_hs.items():
        hs.sort()
        ref_global = sc_indices[hs[0]] + 1
        for later_local in hs[1:]:
            equiv[sc_indices[later_local] + 1] = ref_global
    return equiv


def _find_equiv(atoms, groups):
    """
    Dispatcher: try RDKit-based symmetry detection, fall back to geometry-only.
    Returns {global_1based: ref_global_1based} for all constrained sidechain atoms.
    """
    sc_indices = [i for i, g in enumerate(groups) if g == 'sidechain']
    try:
        return _find_equiv_rdkit(atoms, sc_indices)
    except ImportError:
        print("  Warning: RDKit not found — falling back to geometry-only equivalence "
              "detection (H atoms on the same heavy atom only). Equivalent heavy atoms "
              "such as symmetric methyl groups will not be constrained. "
              "Install RDKit for full symmetry detection.")
        return _find_equiv_geom(atoms, sc_indices)
    except Exception as e:
        print(f"  Warning: RDKit equivalence detection failed ({e}) — "
              "falling back to geometry-only detection.")
        return _find_equiv_geom(atoms, sc_indices)


# ── resp.in ───────────────────────────────────────────────────────────────────

def write_resp_in(capped_pdb, charge, resname, output):
    """
    Write the RESP control file (resp.in).

    Parameters
    ----------
    capped_pdb : str | Path  — capped PDB produced by cap.py
    charge     : int         — net molecular charge of the capped model
    resname    : str         — residue name used as title in the file
    output     : str | Path  — output path (e.g. 'resp.in')
    """
    atoms  = parse_pdb(capped_pdb)
    groups = _classify(atoms)

    sc_equiv = _find_equiv(atoms, groups)

    # Build per-atom constraint table
    table = []
    for global_1, (a, grp) in enumerate(zip(atoms, groups), start=1):
        anum = ATOMIC_NUMBERS[_elem(a['name'])]
        if grp == 'sidechain':
            constraint = sc_equiv.get(global_1, 0)
        else:
            constraint = -1
        table.append((anum, constraint))

    natoms = len(table)
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
        resname,
        f"{charge:5d}{natoms:5d}",
    ]
    for anum, constraint in table:
        lines.append(f"{anum:5d}{constraint:5d}")
    lines += ["", ""]   # two trailing blank lines required by resp

    Path(output).write_text('\n'.join(lines) + '\n')
    print(f"  resp.in  : {output}  ({natoms} atoms, {n_sc} sidechain)")


# ── resp.qin ──────────────────────────────────────────────────────────────────

def write_resp_qin(capped_pdb, output):
    """
    Write the initial charges file (resp.qin).

    Fixed atoms receive their AMBER ff14SB values; free sidechain atoms get 0.0.
    """
    atoms  = parse_pdb(capped_pdb)
    groups = _classify(atoms)

    ace_q = list(ACE_CHARGES)
    nme_q = list(NME_CHARGES)
    ace_i = nme_i = 0

    charges = []
    for grp in groups:
        if   grp == 'ace':       charges.append(ace_q[ace_i]); ace_i += 1
        elif grp == 'nme':       charges.append(nme_q[nme_i]); nme_i += 1
        elif grp == 'bb_N':      charges.append(BACKBONE_N_CHARGES[0])
        elif grp == 'bb_H':      charges.append(BACKBONE_N_CHARGES[1])
        elif grp == 'bb_C':      charges.append(BACKBONE_C_CHARGES[0])
        elif grp == 'bb_O':      charges.append(BACKBONE_C_CHARGES[1])
        else:                    charges.append(0.0)   # sidechain

    n_sc = sum(1 for g in groups if g == 'sidechain')

    lines = []
    for i in range(0, len(charges), 8):
        lines.append(''.join(f"{q:10.6f}" for q in charges[i:i + 8]))

    Path(output).write_text('\n'.join(lines) + '\n')
    print(f"  resp.qin : {output}  ({len(charges)} charges, {n_sc} sidechain free)")
