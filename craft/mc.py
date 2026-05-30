"""
craft.mc
Generate the prepgen main-chain (.mc) file for a capped residue PDB.

The .mc file tells prepgen:
  - Which atom is the N-terminal connection point (HEAD_NAME)
  - Which atom is the C-terminal connection point (TAIL_NAME)
  - Which atom is the alpha-carbon (MAIN_CHAIN)
  - Which atoms to omit from the residue definition (the ACE / NME caps)
  - What AMBER atom type connects on each side (PRE_HEAD_TYPE / POST_TAIL_TYPE)
  - The net charge of the residue in the peptide

Atom names in the .mc file are taken directly from the capped PDB, which
matches the names written into the .ac file by remap_ac_atom_names().
"""

from pathlib import Path
from .cap import parse_pdb


def write_mc(capped_pdb, charge, output):
    """
    Write the prepgen main-chain file (.mc).

    Parameters
    ----------
    capped_pdb : str | Path — capped PDB produced by cap.py
    charge     : int        — net charge of the residue (= total capped charge,
                              since ACE and NME are neutral)
    output     : str | Path — output path (e.g. 'MEO.mc')
    """
    atoms    = parse_pdb(capped_pdb)
    name_map = {i: atoms[i]['name'] for i in range(len(atoms))}  # index → pdb name

    ace_idx  = [i for i, a in enumerate(atoms) if a['resSeq'] == 1]
    res_idx  = [i for i, a in enumerate(atoms) if a['resSeq'] == 2]
    nme_idx  = [i for i, a in enumerate(atoms) if a['resSeq'] == 3]

    # Backbone atoms in residue: N (first), CA (named 'CA'), C (named 'C')
    head_i    = res_idx[0]                                      # N
    ca_i      = next(i for i in res_idx if atoms[i]['name'] == 'CA')
    tail_i    = next(i for i in res_idx if atoms[i]['name'] == 'C')

    omit_idx  = ace_idx + nme_idx

    lines = [
        f"HEAD_NAME {name_map[head_i]}",
        f"TAIL_NAME {name_map[tail_i]}",
        f"MAIN_CHAIN {name_map[ca_i]}",
    ]
    for i in omit_idx:
        lines.append(f"OMIT_NAME {name_map[i]}")
    lines += [
        "PRE_HEAD_TYPE C",
        "POST_TAIL_TYPE N",
        f"CHARGE {float(charge):.1f}",
    ]

    Path(output).write_text('\n'.join(lines) + '\n')
    print(f"  {Path(output).name:<12s}: HEAD={name_map[head_i]}  "
          f"CA={name_map[ca_i]}  TAIL={name_map[tail_i]}  "
          f"omit {len(omit_idx)} cap atoms")
