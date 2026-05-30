"""
parameterize.resp
Generate RESP charge-fitting input files (resp.in and resp.qin).

Constraints follow the standard AMBER RESP protocol:
  - ACE and NME cap atoms are FIXED to AMBER ff14SB charges.
  - Backbone N, H(amide), C, O of the residue are FIXED to ff14SB values.
  - Sidechain atoms are FREE to be fit; equivalent H atoms (same bonded heavy
    atom) are constrained equal to each other.

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
#   ACE: CH3(C)  H1  H2  H3  C  O
#   Backbone N-terminus side: N  H(amide)
#   Backbone C-terminus side: C  O
#   NME: N  H  C  H1  H2  H3

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


def _find_equiv(sidechain_atoms):
    """
    Return {local_idx: ref_local_idx} for H atoms in the sidechain that are
    equivalent to (bonded to the same heavy atom as) an earlier H.
    Only constrained atoms appear as keys; the reference H does not.
    """
    h_to_heavy = {}
    for i, a in enumerate(sidechain_atoms):
        if _elem(a['name']) != 'H':
            continue
        for j, b in enumerate(sidechain_atoms):
            if i == j:
                continue
            if _elem(b['name']) in ('C', 'N', 'O', 'S', 'P'):
                if np.linalg.norm(_pos(a) - _pos(b)) < 1.35:
                    h_to_heavy[i] = j
                    break

    heavy_to_hs = {}
    for h_idx, heavy_idx in h_to_heavy.items():
        heavy_to_hs.setdefault(heavy_idx, []).append(h_idx)

    equiv = {}
    for _, hs in heavy_to_hs.items():
        hs.sort()
        for later_h in hs[1:]:
            equiv[later_h] = hs[0]
    return equiv


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

    # Collect sidechain atoms (with their 1-based global indices) for equiv detection
    sc_pairs = [(i + 1, atoms[i])
                for i, g in enumerate(groups) if g == 'sidechain']
    sc_global  = [idx for idx, _ in sc_pairs]
    sc_dicts   = [a   for _, a   in sc_pairs]

    equiv_local = _find_equiv(sc_dicts)
    # Convert local (sidechain-list) indices to global 1-based indices
    sc_equiv = {sc_global[loc]: sc_global[ref]
                for loc, ref in equiv_local.items()}

    # Build per-atom constraint table
    table = []
    for global_1, (a, grp) in enumerate(zip(atoms, groups), start=1):
        anum = ATOMIC_NUMBERS[_elem(a['name'])]
        if grp == 'sidechain':
            constraint = sc_equiv.get(global_1, 0)
        else:
            constraint = -1
        table.append((anum, constraint))

    natoms  = len(table)
    n_sc    = sum(1 for g in groups if g == 'sidechain')

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
