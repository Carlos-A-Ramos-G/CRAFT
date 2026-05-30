"""
parameterize.cap
Add ACE (N-terminal) and NME (C-terminal) capping groups to a single-residue PDB.

Works regardless of how the residue was obtained:
  - Free amino acid in zwitterionic form (NH3+ / COO-)
  - Cut from a PDB chain (bare N, single amide H, no OXT)
  - Anything in between (NH2, partial H, etc.)

Output order: ACE (resSeq 1) → residue (resSeq 2) → NME (resSeq 3)
Atom naming matches AMBER ff14SB conventions.
"""

import numpy as np
from pathlib import Path


# ── PDB I/O ───────────────────────────────────────────────────────────────────

def parse_pdb(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            try:
                occ  = float(line[54:60])
                bfac = float(line[60:66])
            except (ValueError, IndexError):
                occ, bfac = 1.0, 0.0
            atoms.append({
                'name':    line[12:16].strip(),
                'resName': line[17:20].strip(),
                'chainID': line[21],
                'resSeq':  int(line[22:26]),
                'x': float(line[30:38]),
                'y': float(line[38:46]),
                'z': float(line[46:54]),
                'occupancy':  occ,
                'tempFactor': bfac,
            })
    return atoms


def _fmt_name(name):
    """PDB columns 13-16: 4-char atom name field."""
    return name if len(name) == 4 else f" {name:<3s}"


def _elem(name):
    """Derive element symbol from atom name."""
    return name.lstrip('0123456789')[0].upper()


def pdb_line(serial, name, res, seq, xyz, occ=1.0, bfac=0.0):
    x, y, z = xyz
    return (f"ATOM  {serial:5d} {_fmt_name(name)} {res:<3s}  "
            f"{seq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{occ:6.2f}{bfac:6.2f}          {_elem(name):>2s}\n")


# ── Geometry ──────────────────────────────────────────────────────────────────

def nerf(A, B, C, length, angle_deg, dihedral_deg):
    """
    Place D bonded to C (|CD|=length), angle B-C-D=angle_deg,
    dihedral A-B-C-D=dihedral_deg.

    Sign convention: uses n = cross(A-B, C-B), the negative of the standard
    n1 = cross(B-A, C-B). Input dihedral_deg = standard dihedral + 180°.
    """
    a = np.radians(angle_deg)
    d = np.radians(dihedral_deg)

    bc = C - B
    bc_hat = bc / np.linalg.norm(bc)

    n = np.cross(A - B, bc)
    if np.linalg.norm(n) < 1e-10:
        perp = np.array([1., 0., 0.]) if abs(bc_hat[0]) < 0.9 else np.array([0., 1., 0.])
        n = np.cross(bc_hat, perp)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc_hat)

    return C + length * (
        -np.cos(a) * bc_hat
        + np.sin(a) * np.cos(d) * m
        + np.sin(a) * np.sin(d) * n
    )


def sp2_third(center, p1, p2, bond_length):
    """Third sp2 substituent: opposite to the resultant of the other two bonds."""
    u1 = p1 - center;  u1 /= np.linalg.norm(u1)
    u2 = p2 - center;  u2 /= np.linalg.norm(u2)
    v = -(u1 + u2)
    return center + bond_length * v / np.linalg.norm(v)


def methyl_H(C_pos, bonded_to, bond_length=1.090):
    """Three tetrahedral H positions on a –CH3 group."""
    axis = bonded_to - C_pos
    axis /= np.linalg.norm(axis)
    perp = np.array([1., 0., 0.]) if abs(axis[0]) < 0.9 else np.array([0., 1., 0.])
    perp -= np.dot(perp, axis) * axis
    perp /= np.linalg.norm(perp)
    perp2 = np.cross(axis, perp)
    cos_a = np.cos(np.radians(109.471))
    sin_a = np.sin(np.radians(109.471))
    return [
        C_pos + bond_length * (cos_a * axis
                               + sin_a * (np.cos(np.radians(t)) * perp
                                          + np.sin(np.radians(t)) * perp2))
        for t in (0, 120, 240)
    ]


# ── Terminus inspection ───────────────────────────────────────────────────────

CTERM_ONAMES = {'O1', 'OXT', 'OT1', 'OT2', 'O2'}

def inspect_termini(atoms, N_pos, C_pos):
    """
    Examine the actual state of the N- and C-termini.

    Returns
    -------
    h_on_N   : list of atom indices — H atoms bonded to backbone N
    oxt_idx  : list of atom indices — extra carboxylate O on C
    n_label  : human-readable description of the N-terminus state
    c_label  : human-readable description of the C-terminus state
    """
    h_on_N  = []
    oxt_idx = []

    for i, a in enumerate(atoms):
        p = np.array([a['x'], a['y'], a['z']])
        if a['name'].startswith('H') and np.linalg.norm(p - N_pos) < 1.3:
            h_on_N.append(i)
        if a['name'] in CTERM_ONAMES and np.linalg.norm(p - C_pos) < 1.6:
            oxt_idx.append(i)

    n = len(h_on_N)
    n_label = {
        0: 'bare N — no H attached (cut from PDB, heavy atoms only)',
        1: 'N–H — single amide H (cut from PDB with H, or peptide fragment)',
        2: 'NH2 — neutral free N-terminus',
        3: 'NH3+ — zwitterionic free amino acid',
    }.get(n, f'N with {n} H atoms (unusual)')

    if oxt_idx:
        oxt_name = atoms[oxt_idx[0]]['name']
        c_label = f'COO⁻ — has {oxt_name} ({len(oxt_idx)} extra O), free amino acid C-terminus'
    else:
        c_label = 'C=O only — cut from PDB or already a single carbonyl'

    return h_on_N, oxt_idx, n_label, c_label


# ── Main logic ────────────────────────────────────────────────────────────────

def cap(input_pdb, output_pdb):
    atoms = parse_pdb(input_pdb)
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {input_pdb}")

    pos = {}
    for a in atoms:
        pos.setdefault(a['name'], np.array([a['x'], a['y'], a['z']]))

    missing = [n for n in ('N', 'CA', 'C', 'O') if n not in pos]
    if missing:
        raise ValueError(
            f"Missing backbone atom(s): {missing}\n"
            "The input PDB must contain N, CA, C, and O backbone atoms."
        )

    N, CA, C, O = pos['N'], pos['CA'], pos['C'], pos['O']
    resName = atoms[0]['resName']

    h_on_N, oxt_idx, n_label, c_label = inspect_termini(atoms, N, C)

    # Single amide H (PDB-cut): keep it in place and use it to anchor ACE.
    # For 0 H (bare N) or 2+ H (free amino acid), remove all and recompute.
    if len(h_on_N) == 1:
        h_remove = []
        existing_H_pos = np.array([atoms[h_on_N[0]]['x'],
                                   atoms[h_on_N[0]]['y'],
                                   atoms[h_on_N[0]]['z']])
        inject_amide_H = False
    else:
        h_remove = h_on_N
        existing_H_pos = None
        inject_amide_H = True

    remove_idx = set(h_remove) | set(oxt_idx)

    print(f"Input : {input_pdb}  ({len(atoms)} atoms, residue {resName})")
    print(f"  N-terminus : {n_label}")
    print(f"  C-terminus : {c_label}")
    h_action = ("keep existing amide-H, use it to anchor ACE"
                if existing_H_pos is not None
                else f"remove {len(h_remove)} H(s) on N, add amide-H")
    print(f"  Action     : {h_action}"
          + (f", remove {atoms[oxt_idx[0]]['name']}" if oxt_idx else "")
          + ", add ACE / add NME")

    # ── Place ACE ─────────────────────────────────────────────────────────────
    if existing_H_pos is not None:
        # Derive ACE_C as the third sp2 bond from N, opposite N→H and N→CA.
        ACE_C = sp2_third(N, existing_H_pos, CA, 1.335)
    else:
        # nerf convention: input dihedral = standard dihedral + 180°
        # want standard C–CA–N–ACE_C = 180° (trans ω)  →  pass 0°
        ACE_C = nerf(C, CA, N, 1.335, 121.0, 0.0)
    # want standard CA–N–ACE_C–O  = 0°  (cis, trans amide)  →  pass 180°
    ACE_O   = nerf(CA, N, ACE_C,  1.229, 120.5, 180.0)
    ACE_CH3 = sp2_third(ACE_C, ACE_O, N, 1.522)
    ACE_Hs  = methyl_H(ACE_CH3, ACE_C)

    # ── Backbone amide H ─────────────────────────────────────────────────────
    if inject_amide_H:
        # want standard CA–ACE_C–N–H = 180°  →  pass 0°
        N_H = nerf(CA, ACE_C, N, 1.010, 118.0, 0.0)
    else:
        N_H = None

    # ── Place NME ─────────────────────────────────────────────────────────────
    # C is sp2: NME_N is the third substituent opposite the CA–C–O bisector.
    # Using sp2_third (not NERF) ensures NME_N is placed from the *actual*
    # O position, not a backbone dihedral that may not match where O sits.
    NME_N   = sp2_third(C, CA, O, 1.335)
    # want standard CA–C–NME_N–H  = 0°  (cis to CA)  →  pass 180°
    NME_H   = nerf(CA, C, NME_N,   1.010, 119.0, 180.0)
    NME_CH3 = sp2_third(NME_N, C, NME_H, 1.449)
    NME_Hs  = methyl_H(NME_CH3, NME_N)

    # ── Write output PDB ──────────────────────────────────────────────────────
    ser = 1
    out = []

    # ACE (resSeq 1): CH3 H1 H2 H3 C O
    for name, xyz in [('CH3', ACE_CH3), ('H1', ACE_Hs[0]),
                      ('H2', ACE_Hs[1]), ('H3',  ACE_Hs[2]),
                      ('C',  ACE_C),     ('O',   ACE_O)]:
        out.append(pdb_line(ser, name, 'ACE', 1, xyz)); ser += 1

    # Residue (resSeq 2): write N, then inject new amide H if needed, then rest
    n_written = 0
    for i, a in enumerate(atoms):
        if i in remove_idx:
            continue
        out.append(pdb_line(ser, a['name'], resName, 2,
                            (a['x'], a['y'], a['z']),
                            a['occupancy'], a['tempFactor']))
        ser += 1; n_written += 1
        if a['name'] == 'N' and inject_amide_H:   # inject amide H right after N
            out.append(pdb_line(ser, 'H', resName, 2, N_H))
            ser += 1; n_written += 1

    # NME (resSeq 3): N H C H1 H2 H3
    for name, xyz in [('N', NME_N), ('H',  NME_H),  ('C',  NME_CH3),
                      ('H1', NME_Hs[0]), ('H2', NME_Hs[1]), ('H3', NME_Hs[2])]:
        out.append(pdb_line(ser, name, 'NME', 3, xyz)); ser += 1

    out.append(f'TER   {ser:5d}      NME     3\n')
    out.append('END\n')

    Path(output_pdb).write_text(''.join(out))

    print(f"Output: {output_pdb}")
    print(f"  ACE  resSeq=1 :  6 atoms")
    print(f"  {resName:<4s} resSeq=2 : {n_written} atoms")
    print(f"  NME  resSeq=3 :  6 atoms")
    print(f"  Total         : {ser - 1} atoms")
