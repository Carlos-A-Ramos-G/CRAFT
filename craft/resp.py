"""
craft.resp
Generate RESP charge-fitting input files (resp.in and resp.qin).

Constraints follow the standard AMBER RESP protocol adapted for each
terminal position:

  middle : fix ACE + backbone N/H/CA/HA/C/O + NME; sidechain atoms free
  cterm  : fix ACE + backbone N/H/CA/HA; C-terminal C/O/OXT and sidechain free
  nterm  : fix NME + backbone CA/HA/C/O; N-terminal N/H atoms and sidechain free

CA and HA are fixed such that the six backbone atoms (N, H, CA, HA, C, O)
sum to exactly zero net charge, which prevents charge transfer artifacts at
QM/MM boundaries.  Alpha-H atoms are detected geometrically (distance to CA
< 1.15 Å) so the correct number is found automatically for glycine (two) or
alpha-substituted residues with no alpha-H.

Backbone atoms are identified by exact name ('N', 'H', 'CA', 'C', 'O') for
the fixed end, and by residue sequence number for cap atoms.

Reference: Bayly et al., J. Phys. Chem. 97, 10269 (1993).
"""

import numpy as np
from pathlib import Path
from .cap import parse_pdb, _elem


# -- AMBER ff14SB charges for fixed atoms --------------------------------------
# Order matches the capped-PDB atom ordering produced by cap.py:
#   ACE: methyl-C  methyl-H1  methyl-H2  methyl-H3  carbonyl-C  carbonyl-O
#   Backbone N-terminus side: N  H(amide)
#   Backbone C-terminus side: C  O
#   NME: amide-N  amide-H  methyl-C  methyl-H1  methyl-H2  methyl-H3
# (Actual atom names in the PDB are unique per-molecule; charges map by position.)

ACE_CHARGES        = [-0.3662,  0.1123,  0.1123,  0.1123,  0.5972, -0.5679]
BACKBONE_N_CHARGES = [-0.4157,  0.2719]   # backbone N, amide H
BACKBONE_CA_CHARGE =  0.0337             # alpha carbon (ff14SB)
BACKBONE_HA_TOTAL  =  0.0808             # total alpha-H charge, split equally
                                         # among all H bonded to CA; 0.0015
                                         # below ff14SB (0.0823) so that
                                         # N+H+CA+HA+C+O = 0.0000 exactly
BACKBONE_C_CHARGES = [ 0.5972, -0.5679]  # backbone C (carbonyl), O
NME_CHARGES        = [-0.4157,  0.2719, -0.1490,  0.0976,  0.0976,  0.0976]

ATOMIC_NUMBERS = {
    'H': 1, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'P': 15, 'S': 16, 'Cl': 17, 'Br': 35,
}


# -- Atom classification -------------------------------------------------------

def _classify(atoms, position='middle'):
    """
    Classify every atom as one of:
      'ace' | 'bb_N' | 'bb_H' | 'bb_CA' | 'bb_HA' |
      'sidechain' | 'bb_C' | 'bb_O' | 'nme'

    CA and alpha-H atoms are always fixed regardless of position so that the
    six backbone atoms (N, H, CA, HA, C, O) sum to exactly zero net charge.
    Alpha-H atoms are detected geometrically (H in resSeq=2 within 1.15 Å of
    CA), which handles glycine (two alpha-H) and alpha-substituted residues
    with no alpha-H without special cases.

    Which other backbone atoms are fixed depends on position:
      middle : N, H(amide), CA, HA, C, O all fixed
      cterm  : N, H(amide), CA, HA fixed; C-terminus (C, O, OXT) free
      nterm  : CA, HA, C, O fixed; N-terminus (N and its H atoms) free

    The amide H is identified as the second residue atom in the output PDB
    (cap.py always writes N then the amide H for middle/cterm modes).
    """
    res_indices    = [i for i, a in enumerate(atoms) if a['resSeq'] == 2]
    amide_h_global = res_indices[1]   # only meaningful for middle and cterm

    # Detect alpha-H: H in resSeq=2 within 1.15 Å of CA
    ca_atom = next((a for a in atoms if a['resSeq'] == 2 and a['name'] == 'CA'), None)
    alpha_h_indices = set()
    if ca_atom is not None:
        ca_pos = _pos(ca_atom)
        for i, a in enumerate(atoms):
            if a['resSeq'] == 2 and _elem(a['name']) == 'H':
                if np.linalg.norm(_pos(a) - ca_pos) < 1.15:
                    alpha_h_indices.add(i)

    groups = []
    for i, a in enumerate(atoms):
        if a['resSeq'] == 1:
            groups.append('ace')
        elif a['resSeq'] == 3:
            groups.append('nme')
        else:
            # residue atom -- classification depends on position
            if position == 'middle':
                if a['name'] == 'N':
                    groups.append('bb_N')
                elif i == amide_h_global:
                    groups.append('bb_H')
                elif a['name'] == 'CA':
                    groups.append('bb_CA')
                elif i in alpha_h_indices:
                    groups.append('bb_HA')
                elif a['name'] == 'C':
                    groups.append('bb_C')
                elif a['name'] == 'O':
                    groups.append('bb_O')
                else:
                    groups.append('sidechain')
            elif position == 'cterm':
                # ACE-capped N-terminus is fixed; C-terminus is free
                if a['name'] == 'N':
                    groups.append('bb_N')
                elif i == amide_h_global:
                    groups.append('bb_H')
                elif a['name'] == 'CA':
                    groups.append('bb_CA')
                elif i in alpha_h_indices:
                    groups.append('bb_HA')
                else:
                    groups.append('sidechain')
            else:  # nterm
                # NME-capped C-terminus is fixed; N-terminus is free
                if a['name'] == 'CA':
                    groups.append('bb_CA')
                elif i in alpha_h_indices:
                    groups.append('bb_HA')
                elif a['name'] == 'C':
                    groups.append('bb_C')
                elif a['name'] == 'O':
                    groups.append('bb_O')
                else:
                    groups.append('sidechain')

    return groups


# -- Equivalence detection -----------------------------------------------------

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
        print("  Warning: RDKit not found -- falling back to geometry-only equivalence "
              "detection (H atoms on the same heavy atom only). Equivalent heavy atoms "
              "such as symmetric methyl groups will not be constrained. "
              "Install RDKit for full symmetry detection.")
        return _find_equiv_geom(atoms, sc_indices)
    except Exception as e:
        print(f"  Warning: RDKit equivalence detection failed ({e}) -- "
              "falling back to geometry-only detection.")
        return _find_equiv_geom(atoms, sc_indices)


# -- resp.in -------------------------------------------------------------------

def write_resp_in(capped_pdb, charge, resname, output, position='middle'):
    """
    Write the RESP control file (resp.in).

    Parameters
    ----------
    capped_pdb : str | Path  -- capped PDB produced by cap.py
    charge     : int         -- net molecular charge of the capped model
    resname    : str         -- residue name used as title in the file
    output     : str | Path  -- output path (e.g. 'resp.in')
    position   : str         -- 'middle', 'cterm', or 'nterm'
    """
    atoms  = parse_pdb(capped_pdb)
    groups = _classify(atoms, position)

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
    print(f"  resp.in  : {output}  ({natoms} atoms, {n_sc} sidechain free)")


# -- resp.qin ------------------------------------------------------------------

def write_resp_qin(capped_pdb, output, position='middle'):
    """
    Write the initial charges file (resp.qin).

    Fixed atoms receive their AMBER ff14SB values; free atoms get 0.0.
    """
    atoms  = parse_pdb(capped_pdb)
    groups = _classify(atoms, position)

    # Per-HA charge: split BACKBONE_HA_TOTAL equally; if no alpha-H exists
    # (e.g. alpha-substituted residues), fold the full amount onto CA.
    n_bb_ha  = groups.count('bb_HA')
    ca_charge = BACKBONE_CA_CHARGE + (BACKBONE_HA_TOTAL if n_bb_ha == 0 else 0.0)
    ha_charge = BACKBONE_HA_TOTAL / n_bb_ha if n_bb_ha > 0 else 0.0

    ace_q = list(ACE_CHARGES)
    nme_q = list(NME_CHARGES)
    ace_i = nme_i = 0

    charges = []
    for grp in groups:
        if   grp == 'ace':   charges.append(ace_q[ace_i]); ace_i += 1
        elif grp == 'nme':   charges.append(nme_q[nme_i]); nme_i += 1
        elif grp == 'bb_N':  charges.append(BACKBONE_N_CHARGES[0])
        elif grp == 'bb_H':  charges.append(BACKBONE_N_CHARGES[1])
        elif grp == 'bb_CA': charges.append(ca_charge)
        elif grp == 'bb_HA': charges.append(ha_charge)
        elif grp == 'bb_C':  charges.append(BACKBONE_C_CHARGES[0])
        elif grp == 'bb_O':  charges.append(BACKBONE_C_CHARGES[1])
        else:                charges.append(0.0)   # sidechain / free terminal

    n_sc = sum(1 for g in groups if g == 'sidechain')

    lines = []
    for i in range(0, len(charges), 8):
        lines.append(''.join(f"{q:10.6f}" for q in charges[i:i + 8]))

    Path(output).write_text('\n'.join(lines) + '\n')
    print(f"  resp.qin : {output}  ({len(charges)} charges, {n_sc} sidechain free)")
