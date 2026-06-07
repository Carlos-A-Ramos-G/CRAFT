"""
craft.amber
Run the post-Gaussian AMBER parameterization pipeline:

  espgen      - extract ESP grid from Gaussian HF log
  resp        - fit RESP charges
  antechamber - assign AMBER atom types and write .ac file
  prepgen     - build residue topology (.prepin)
  parmchk2    - generate missing force-field parameters (.frcmod)
               twice: once for GAFF, once for the target protein FF

All tools must be available in $PATH (i.e. AMBER or AmberTools installed).

Output file names carry a position prefix for terminal variants:
  middle : {resname}.ac, {resname}.prepin, {resname}_gaff.frcmod, ...
  cterm  : C{resname}.ac, C{resname}.prepin, C{resname}_gaff.frcmod, ...
  nterm  : N{resname}.ac, N{resname}.prepin, N{resname}_gaff.frcmod, ...
"""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


_GAFF2_TO_AMBER = {
    'c':  'C',   'c1': 'CX', 'c2': 'CA', 'c3': 'CT', 'ca': 'CA',
    'cb': 'CB',  'cc': 'CA', 'cd': 'CA', 'ce': 'CA', 'cf': 'CA',
    'cp': 'CA',  'cq': 'CA',
    'n':  'N',   'n1': 'N1', 'n2': 'N2', 'n3': 'N3', 'n4': 'N3',
    'na': 'NA',  'nb': 'NA', 'nh': 'N2',
    'o':  'O',   'oh': 'OH', 'os': 'OS', 'op': 'OS', 'oq': 'OS',
    'sh': 'SH',  'ss': 'S',  's2': 'S',  's4': 'S',  's6': 'S',
    'sx': 'S',   'sy': 'S',
    'p2': 'P',   'p3': 'P',  'p4': 'P',  'p5': 'P',
    'hc': 'HC',  'ha': 'HA', 'hn': 'H',  'ho': 'HO', 'hs': 'HS',
    'hp': 'HP',  'h1': 'H1', 'h2': 'H2', 'h3': 'H3', 'h4': 'H4', 'h5': 'H5',
    'f':  'F',   'cl': 'Cl', 'br': 'Br', 'i':  'I',
}


def _run(cmd, cwd):
    """Run a shell command, raise on failure, print output."""
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout.rstrip())
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {result.returncode}):\n"
            f"  {' '.join(cmd)}\n"
            f"{result.stderr}"
        )


def _postprocess_frcmod(frcmod_path):
    """
    Post-process a parmchk2 frcmod in place:
      - Remove lines containing 'ATTN, need revision' (zero placeholder parameters
        that parmchk2 writes when it cannot locate any suitable parameter).
      - Warn and list any parameters whose penalty score exceeds 100.
    """
    path  = Path(frcmod_path)
    lines = path.read_text().splitlines(keepends=True)

    clean        = []
    high_penalty = []

    for line in lines:
        if 'ATTN' in line:
            continue
        m = re.search(r'penalty score=\s*([\d.]+)', line)
        if m and float(m.group(1)) > 100:
            high_penalty.append(line.rstrip())
        clean.append(line)

    path.write_text(''.join(clean))

    if high_penalty:
        print(f"\n  Warning: it's worth noticing that there are parameters with a "
              f"penalty score above 100 in {path.name}:")
        for ln in high_penalty:
            print(f"    {ln.strip()}")


def remap_ac_atom_names(ac_path, capped_pdb_path):
    """
    Rename atoms in an antechamber .ac file to match a capped PDB.

    Builds a positional mapping (ac_name[i] -> pdb_name[i]) from the order
    of ATOM records and applies it to all ATOM and BOND lines.  Overwrites
    ac_path in place.

    Parameters
    ----------
    ac_path         : str | Path -- antechamber .ac file to update
    capped_pdb_path : str | Path -- capped PDB whose atom names are the target
    """
    ac_path         = Path(ac_path)
    capped_pdb_path = Path(capped_pdb_path)

    ac_names = [
        ln.split()[2]
        for ln in ac_path.read_text().splitlines()
        if ln.startswith('ATOM')
    ]
    pdb_names = [
        ln.split()[2]
        for ln in capped_pdb_path.read_text().splitlines()
        if ln.startswith('ATOM')
    ]

    if len(ac_names) != len(pdb_names):
        raise ValueError(
            f"Atom count mismatch: {ac_path.name} has {len(ac_names)} atoms, "
            f"{capped_pdb_path.name} has {len(pdb_names)} atoms."
        )

    name_map = dict(zip(ac_names, pdb_names))

    out_lines = []
    for line in ac_path.read_text().splitlines(keepends=True):
        nl = '\n' if line.endswith('\n') else ''
        s  = line.rstrip('\n')

        if s.startswith('ATOM'):
            m = re.match(r'^(ATOM\s+\d+)\s+(\S+)\s+(\S+\s+.*)$', s)
            if m:
                new_name = name_map.get(m.group(2), m.group(2))
                s = f"{m.group(1)}  {new_name:<4}{m.group(3)}"

        elif s.startswith('BOND'):
            m = re.match(r'^(BOND\s+\d+\s+\d+\s+\d+\s+\d+\s+)(\S+)(\s+)(\S+)(.*)$', s)
            if m:
                new1 = name_map.get(m.group(2), m.group(2))
                new2 = name_map.get(m.group(4), m.group(4))
                s = m.group(1) + new1 + m.group(3) + new2 + m.group(5)

        out_lines.append(s + nl)

    ac_path.write_text(''.join(out_lines))
    print(f"  Remapped {len(name_map)} atom names in {ac_path.name}")
    return name_map


def _read_ac_types(ac_path):
    """Return {atom_name: atom_type} for every ATOM record in an .ac file."""
    types = {}
    for line in Path(ac_path).read_text().splitlines():
        if line.startswith('ATOM'):
            parts = line.split()
            if len(parts) >= 10:
                types[parts[2]] = parts[-1]
    return types


def _patch_ac_types(ac_path, corrections):
    """Replace atom types in-place for the specified atom names."""
    ac_path = Path(ac_path)
    out = []
    for line in ac_path.read_text().splitlines(keepends=True):
        stripped = line.rstrip()
        if stripped.startswith('ATOM'):
            parts = stripped.split()
            if len(parts) >= 10 and parts[2] in corrections:
                trailing = line[len(stripped):]
                stripped  = re.sub(r'\S+$', corrections[parts[2]], stripped)
                line      = stripped + trailing
        out.append(line)
    ac_path.write_text(''.join(out))


def resolve_du_atom_types(ac_path, hf_log, wd, bk_resname, total_charge,
                          atom_type, user_overrides=None):
    """
    Scan ac_path for DU atom types, suggest corrections via a gaff2 probe run,
    merge with user_overrides (overrides take priority), and patch the .ac file
    in place.  Hard-stops if any DU atom cannot be resolved.
    """
    user_overrides = user_overrides or {}
    ac_path        = Path(ac_path)

    ac_types = _read_ac_types(ac_path)
    du_names = [name for name, t in ac_types.items() if t == 'DU']
    if not du_names:
        return

    print("\n-- resolving DU atom types -------------------------------------------")
    print(f"  antechamber assigned DU to: {', '.join(du_names)}")

    gaff2_types = {}
    suggestions  = {}

    if atom_type != 'gaff2':
        print("  re-running with -at gaff2 to generate type suggestions ...")
        probe_name = '_du_probe.ac'
        probe_path = Path(wd) / probe_name
        try:
            _run([
                'antechamber',
                '-fi', 'gout', '-i', hf_log,
                '-bk', bk_resname,
                '-fo', 'ac',   '-o', probe_name,
                '-c',  'rc',   '-cf', 'resp.chg',
                '-at', 'gaff2',
                '-nc', str(total_charge),
            ], cwd=wd)
            gaff2_types = _read_ac_types(probe_path)
        except RuntimeError:
            print("  Warning: gaff2 probe failed — auto-suggestions unavailable.")
        finally:
            probe_path.unlink(missing_ok=True)

        for name in du_names:
            g2 = gaff2_types.get(name, 'DU')
            if g2 == 'DU':
                continue
            if atom_type == 'amber':
                amber = _GAFF2_TO_AMBER.get(g2.lower())
                if amber:
                    suggestions[name] = amber
            else:
                suggestions[name] = g2
    else:
        print("  atom_type is gaff2 — probe would be identical; skipping auto-suggest.")

    # Resolution table
    w = max((len(n) for n in du_names), default=4) + 2
    print()
    if gaff2_types:
        print(f"  {'Atom':<{w}} {'GAFF2':<8} {'Applied ({})'.format(atom_type):<22} Source")
        print(f"  {'-'*w} {'-'*8} {'-'*22} {'-'*28}")
        for name in du_names:
            g2 = gaff2_types.get(name, '--')
            if name in user_overrides:
                auto = suggestions.get(name, '--')
                apl  = user_overrides[name]
                src  = f"config override  (auto: {auto})"
            elif name in suggestions:
                apl  = suggestions[name]
                src  = "auto-suggestion"
            else:
                apl  = "--"
                src  = "UNRESOLVED"
            print(f"  {name:<{w}} {g2:<8} {apl:<22} {src}")
    else:
        print(f"  {'Atom':<{w}} {'Applied ({})'.format(atom_type):<22} Source")
        print(f"  {'-'*w} {'-'*22} {'-'*28}")
        for name in du_names:
            if name in user_overrides:
                print(f"  {name:<{w}} {user_overrides[name]:<22} config override")
            else:
                print(f"  {name:<{w}} {'--':<22} UNRESOLVED")

    # Determine corrections
    corrections = {}
    unresolved  = []
    for name in du_names:
        if name in user_overrides:
            corrections[name] = user_overrides[name]
        elif name in suggestions:
            corrections[name] = suggestions[name]
        else:
            unresolved.append(name)

    if unresolved:
        print()
        print(f"  Error: unresolved DU type(s): {', '.join(unresolved)}")
        print(f"  Add atom_type_overrides to your config:")
        print(f"    amber:")
        print(f"      atom_type_overrides:")
        for name in unresolved:
            g2_hint = gaff2_types.get(name, '?')
            print(f"        {name}: <type>  # gaff2 assigned: {g2_hint}")
        sys.exit(1)

    auto_applied = {n: t for n, t in corrections.items() if n not in user_overrides}
    if auto_applied:
        print()
        print("  Auto-suggestions applied. To override, add to config:")
        print("    amber:")
        print("      atom_type_overrides:")
        for name, atype in auto_applied.items():
            print(f"        {name}: {atype}")

    _patch_ac_types(ac_path, corrections)
    applied = ',  '.join(f"{n}: DU → {t}" for n, t in corrections.items())
    print(f"\n  Patched {ac_path.name}: {applied}")


PARM_FILES = {
    'ff14SB': 'parm10.dat',
    'ff19SB': 'parm19.dat',
}


def run_amber_pipeline(hf_log, resname, charge, mc_file,
                       workdir='.', atom_type='amber', forcefield='ff14SB',
                       capped_pdb=None, position='middle',
                       atom_type_overrides=None):
    """
    Run the full post-Gaussian parameterization pipeline.

    Parameters
    ----------
    hf_log      : str | Path -- Gaussian HF/ESP log (MEO_hf.log)
    resname     : str        -- three-letter residue code (e.g. 'MEO')
    charge      : int        -- net molecular charge of the capped model
    mc_file     : str | Path -- prepgen main-chain file (.mc)
    workdir     : str | Path -- directory where all files are written
    atom_type   : str        -- antechamber -at flag: 'amber', 'gaff', or 'gaff2'
    forcefield  : str | None -- protein FF for the second parmchk2 run:
                                'ff14SB' (parm10.dat), 'ff19SB' (parm19.dat),
                                or None to skip
    capped_pdb  : str | Path | None -- capped PDB for atom name remapping;
                  defaults to {base}_capped.pdb in workdir; pass None to skip
    position    : str        -- 'middle', 'cterm', or 'nterm'; determines the
                                output file name prefix (C/N/none)
    """
    if position not in ('middle', 'cterm', 'nterm'):
        raise ValueError(
            f"position must be 'middle', 'cterm', or 'nterm'; got {position!r}")

    prefix = {'middle': '', 'cterm': 'C', 'nterm': 'N'}[position]
    base   = f"{prefix}{resname}"

    hf_log  = str(Path(hf_log).resolve())
    mc_file = str(Path(mc_file).resolve())
    wd      = str(Path(workdir).resolve())

    ac_file     = f"{base}.ac"
    prepin_file = f"{base}.prepin"
    gaff_frcmod = f"{base}_gaff.frcmod"
    ff_frcmod   = f"{base}_{forcefield}.frcmod" if forcefield else None

    amberhome = os.environ.get('AMBERHOME', '')
    if forcefield:
        parm_name = PARM_FILES.get(forcefield)
        if parm_name is None:
            raise ValueError(
                f"Unknown forcefield {forcefield!r}. "
                f"Supported: {list(PARM_FILES)}"
            )
        parm_file = str(Path(amberhome) / 'dat/leap/parm' / parm_name) if amberhome else parm_name

    print("\n-- espgen ------------------------------------------------------------")
    _run(['espgen', '-i', hf_log, '-o', 'esp.dat'], cwd=wd)

    print("\n-- resp --------------------------------------------------------------")
    _run([
        'resp', '-O',
        '-i', 'resp.in',
        '-o', 'resp.out',
        '-p', 'resp.pch',
        '-t', 'resp.chg',
        '-q', 'resp.qin',
        '-e', 'esp.dat',
    ], cwd=wd)

    print("\n-- antechamber -------------------------------------------------------")
    _run([
        'antechamber',
        '-fi', 'gout',
        '-i',  hf_log,
        '-bk', resname,
        '-fo', 'ac',
        '-o',  ac_file,
        '-c',  'rc',
        '-cf', 'resp.chg',
        '-at', atom_type,
        '-nc', str(charge),
    ], cwd=wd)

    print("\n-- remap atom names -------------------------------------------------")
    if capped_pdb:
        _capped = Path(capped_pdb)
    else:
        _capped = Path(wd) / f"{base}_capped.pdb"
        if not _capped.exists():
            _capped = Path(wd) / f"{resname}_capped.pdb"

    if _capped.exists():
        remap_ac_atom_names(Path(wd) / ac_file, _capped)
    else:
        print(f"  {_capped.name} not found -- atom names in {ac_file} not remapped.")

    resolve_du_atom_types(
        Path(wd) / ac_file, hf_log, wd,
        resname, charge, atom_type, atom_type_overrides,
    )

    # prepgen has a short fixed-length path buffer (~256 chars); long absolute
    # paths are silently truncated. Copy the .mc file into workdir and pass
    # only the basename to avoid this.
    mc_local = Path(mc_file).name
    mc_dst   = Path(wd) / mc_local
    if mc_dst.resolve() != Path(mc_file).resolve():
        shutil.copy(mc_file, mc_dst)

    print("\n-- prepgen -----------------------------------------------------------")
    _run([
        'prepgen',
        '-i',  ac_file,
        '-o',  prepin_file,
        '-m',  mc_local,
        '-rn', resname,
    ], cwd=wd)

    print("\n-- parmchk2 (GAFF) ---------------------------------------------------")
    _run([
        'parmchk2',
        '-i', ac_file,
        '-f', 'ac',
        '-o', gaff_frcmod,
    ], cwd=wd)
    _postprocess_frcmod(Path(wd) / gaff_frcmod)

    if forcefield:
        print(f"\n-- parmchk2 ({forcefield}) {'-' * (51 - len(forcefield))}")
        _run([
            'parmchk2',
            '-i', ac_file,
            '-f', 'ac',
            '-o', ff_frcmod,
            '-a', 'Y',
            '-p', parm_file,
        ], cwd=wd)
        _postprocess_frcmod(Path(wd) / ff_frcmod)

    output_files = [ac_file, prepin_file, gaff_frcmod]
    if forcefield:
        output_files.append(ff_frcmod)

    print(f"\nDone. Output files in {wd}:")
    for f in output_files:
        p = Path(wd) / f
        status = "ok" if p.exists() else "MISSING"
        print(f"  [{status}]  {f}")
