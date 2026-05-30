"""
parameterize.gaussian
Generate a Gaussian input (.com) from a capped PDB file.

The checkpoint name is derived from the PDB stem by stripping a trailing
'_capped' suffix, e.g. MEO_capped.pdb → MEO_opt.chk / MEO_opt.com.
"""

from pathlib import Path


NPROC_DEFAULT = 16
MEM_DEFAULT   = "512MB"
ROUTE_DEFAULT = "#P b3lyp/6-31g* opt"


def _elem(name):
    """Element symbol from a PDB atom name."""
    return name.lstrip('0123456789')[0].upper()


def parse_pdb(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            atoms.append({
                'name': line[12:16].strip(),
                'x':    float(line[30:38]),
                'y':    float(line[38:46]),
                'z':    float(line[46:54]),
            })
    return atoms


def write_com(pdb_path, com_path, charge, mult,
              nproc=NPROC_DEFAULT, mem=MEM_DEFAULT, route=ROUTE_DEFAULT):
    atoms = parse_pdb(pdb_path)
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records in {pdb_path}")

    stem     = Path(pdb_path).stem
    base     = stem.replace('_capped', '')
    chk_name = f"{base}_opt.chk"
    title    = f"{base}  b3lyp/6-31g* geometry optimisation"

    blocks = [
        f"%nprocshared={nproc}",
        f"%mem={mem}",
        f"%chk={chk_name}",
        route,
        "",
        title,
        "",
        f"{charge} {mult}",
    ]

    for a in atoms:
        elem = _elem(a['name'])
        blocks.append(f" {elem:<2s}  {a['x']:12.6f}  {a['y']:12.6f}  {a['z']:12.6f}")

    blocks.append("")   # Gaussian requires a trailing blank line

    Path(com_path).write_text('\n'.join(blocks) + '\n')

    print(f"Input  : {pdb_path}  ({len(atoms)} atoms)")
    print(f"Output : {com_path}")
    print(f"  Charge / mult : {charge} / {mult}")
    print(f"  Checkpoint    : {chk_name}")
    print(f"  nproc / mem   : {nproc} / {mem}")
