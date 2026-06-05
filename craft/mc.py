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

For terminal variants:
  middle : HEAD_NAME N, TAIL_NAME C, omit ACE + NME
  cterm  : HEAD_NAME N only (no TAIL), omit ACE only
  nterm  : TAIL_NAME C only (no HEAD), omit NME only

Atom names in the .mc file are taken directly from the capped PDB, which
matches the names written into the .ac file by remap_ac_atom_names().
"""

from pathlib import Path
from .cap import parse_pdb


def write_mc(capped_pdb, charge, output, position='middle'):
    """
    Write the prepgen main-chain file (.mc).

    Parameters
    ----------
    capped_pdb : str | Path -- capped PDB produced by cap.py
    charge     : int        -- net charge of the residue (= total capped charge,
                               since ACE and NME are neutral)
    output     : str | Path -- output path (e.g. 'MEO.mc')
    position   : str        -- 'middle', 'cterm', or 'nterm'
    """
    atoms    = parse_pdb(capped_pdb)

    ace_idx  = [i for i, a in enumerate(atoms) if a['resSeq'] == 1]
    res_idx  = [i for i, a in enumerate(atoms) if a['resSeq'] == 2]
    nme_idx  = [i for i, a in enumerate(atoms) if a['resSeq'] == 3]

    head_i   = res_idx[0]
    ca_i     = next(i for i in res_idx if atoms[i]['name'] == 'CA')
    tail_i   = next(i for i in res_idx if atoms[i]['name'] == 'C')

    if position == 'cterm':
        omit_idx = ace_idx
    elif position == 'nterm':
        omit_idx = nme_idx
    else:
        omit_idx = ace_idx + nme_idx

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

    info_parts = []
    if position in ('middle', 'cterm'):
        info_parts.append(f"HEAD={atoms[head_i]['name']}")
    info_parts.append(f"CA={atoms[ca_i]['name']}")
    if position in ('middle', 'nterm'):
        info_parts.append(f"TAIL={atoms[tail_i]['name']}")
    info_parts.append(f"omit {len(omit_idx)} cap atoms")
    print(f"  {Path(output).name:<16s}: {' '.join(info_parts)}")
