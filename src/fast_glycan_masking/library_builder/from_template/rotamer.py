#!/usr/bin/env python3
"""
rotamer_library_generator.py

Build a rotamer library from a single template molecule (e.g. the R1 / MTSSL
spin label) by systematically scanning a set of torsion angles, discarding any
conformer that self-clashes.

This is a standalone Python re-implementation of the idea behind the MATLAB
``RX_labeling.m`` rotamer-building step: instead of driving PyMOL, everything is
done with numpy so it runs headless.

Inputs
------
1. A PDB file holding ONE molecule/residue (the label template). Atom names must
   be unique because the torsions are addressed by atom name. Hydrogens are
   stripped by default (``--keep-hydrogens`` to keep them).
2. A YAML config with two sections:

     torsions:                       # part 1 -- what to sample
       - name: chi1
         atoms: [N, CA, CB, SG]      # the 4 atom names that define the dihedral
         range: [0, 360]             # absolute dihedral value, degrees
         step: 30                    # spacing, degrees
       - ...
     clash:                          # part 2 -- self-clash rule
       cutoff: 2.5                   # angstrom; a pair closer than this = clash
       exclude_bonds: 3              # skip pairs separated by <= this many bonds

Semantics
---------
* For each torsion the sampled values are the ABSOLUTE dihedral targets in
  ``range`` spaced by ``step`` (0-360 space). For a full turn (span == 360) the
  duplicate 360==0 endpoint is dropped automatically.
* Every combination of torsion values (Cartesian product) is generated. For each
  combination the template is taken, each dihedral is rotated to its target (the
  downstream atoms move; the rest stay put), then the conformer is clash-checked.
* Self-clash: any pair of atoms whose bond-graph separation is greater than
  ``exclude_bonds`` and whose distance is below ``cutoff`` marks the rotamer as
  clashed, and it is discarded.
* Output frame == input frame. Provide the template already positioned in the
  reference frame you want the library in (e.g. CA at the origin).

Outputs (both written by default)
---------------------------------
* ``<prefix>.pdb``  : multi-model PDB, one MODEL per surviving rotamer.
* ``<prefix>.npz``  : coords (R, n_atoms, 3), atom_names, elements, resname,
                      dihedrals (R, n_torsion) and torsion_names.

Usage
-----
    python3 rotamer_library_generator.py --pdb R1_template.pdb \\
            --config r1_config.yaml --out-prefix R1_library

    python3 rotamer_library_generator.py --self-test     # verify the engine
"""

import argparse
import itertools
import sys

import numpy as np

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


# --------------------------------------------------------------------------- #
# Chemistry helpers
# --------------------------------------------------------------------------- #

# Cordero et al. (2008) covalent radii, angstrom. Enough for organic labels.
COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07,
    "F": 0.57, "CL": 1.02, "BR": 1.20, "I": 1.39, "SE": 1.20,
}


def infer_element(atom_name):
    """Infer the element from a PDB atom name.

    Macromolecular/label atom names encode the element as the leading letter
    (CA/CB/CE/C3 -> C, SG/SD -> S, N1 -> N, O1 -> O, HB2/H91 -> H). We therefore
    take the first alphabetic character rather than trusting columns 77-78,
    which are garbled in the MMM "pseudo-PDB" files.
    """
    s = atom_name.strip().lstrip("0123456789")
    return s[0].upper() if s else "C"


def covalent_radius(element):
    return COVALENT_RADII.get(element.upper(), 0.77)


# --------------------------------------------------------------------------- #
# PDB parsing / writing
# --------------------------------------------------------------------------- #

class Atom:
    __slots__ = ("serial", "name", "resname", "chain", "resseq", "element", "xyz")

    def __init__(self, serial, name, resname, chain, resseq, element, xyz):
        self.serial = serial
        self.name = name
        self.resname = resname
        self.chain = chain
        self.resseq = resseq
        self.element = element
        self.xyz = np.asarray(xyz, dtype=float)


def read_pdb(path, keep_hydrogens=False):
    """Return (atoms, conect) where atoms is a list[Atom] and conect maps a
    serial to a list of bonded serials (empty if no CONECT records)."""
    atoms = []
    conect = {}
    with open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec in ("ATOM", "HETATM"):
                name = line[12:16].strip()
                element = line[76:78].strip()
                if not element or not element.isalpha():
                    element = infer_element(name)
                element = element.upper()
                if element == "H" and not keep_hydrogens:
                    continue
                try:
                    serial = int(line[6:11])
                except ValueError:
                    serial = len(atoms) + 1
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    resseq = 1
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                atoms.append(Atom(serial, name, line[17:20].strip(),
                                  line[21:22].strip() or "A", resseq, element, xyz))
            elif rec == "CONECT":
                nums = [line[i:i + 5] for i in range(6, len(line.rstrip()), 5)]
                nums = [int(x) for x in nums if x.strip().lstrip("-").isdigit()]
                if nums:
                    conect.setdefault(nums[0], []).extend(nums[1:])
    return atoms, conect


def _format_atom_name(name):
    """PDB atom-name field (cols 13-16). Names <=3 chars get a leading space."""
    if len(name) >= 4:
        return name[:4]
    return " " + name.ljust(3)


def write_multimodel_pdb(path, atoms, coords_all, dihedrals=None, torsion_names=None):
    """coords_all: (R, n_atoms, 3)."""
    with open(path, "w") as fh:
        fh.write("REMARK  Rotamer library generated by rotamer_library_generator.py\n")
        fh.write("REMARK  %d rotamers, %d atoms each\n"
                 % (coords_all.shape[0], len(atoms)))
        if torsion_names:
            fh.write("REMARK  torsions: %s\n" % ", ".join(torsion_names))
        for r in range(coords_all.shape[0]):
            fh.write("MODEL %8d\n" % (r + 1))
            if dihedrals is not None and torsion_names:
                pairs = ", ".join("%s=%.1f" % (n, v)
                                  for n, v in zip(torsion_names, dihedrals[r]))
                fh.write("REMARK  %s\n" % pairs)
            for i, at in enumerate(atoms):
                x, y, z = coords_all[r, i]
                fh.write("HETATM%5d %s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
                         % (i + 1, _format_atom_name(at.name), "",
                            at.resname[:3], at.chain[:1], at.resseq,
                            x, y, z, 1.0, 0.0, at.element.rjust(2)))
            fh.write("ENDMDL\n")
        fh.write("END\n")


# --------------------------------------------------------------------------- #
# Bond graph
# --------------------------------------------------------------------------- #

def build_bonds(atoms, conect=None, tolerance=1.3):
    """Return adjacency as list[set[int]] over atom indices.

    Uses CONECT records if provided, otherwise perceives bonds by distance:
    a pair is bonded when d < tolerance * (r_i + r_j) covalent radii.
    """
    n = len(atoms)
    adj = [set() for _ in range(n)]
    if conect:
        serial_to_idx = {at.serial: i for i, at in enumerate(atoms)}
        used = False
        for a_serial, partners in conect.items():
            if a_serial not in serial_to_idx:
                continue
            i = serial_to_idx[a_serial]
            for b_serial in partners:
                j = serial_to_idx.get(b_serial)
                if j is not None and j != i:
                    adj[i].add(j)
                    adj[j].add(i)
                    used = True
        if used:
            return adj  # trust explicit connectivity

    coords = np.array([at.xyz for at in atoms])
    radii = np.array([covalent_radius(at.element) for at in atoms])
    for i in range(n):
        d = np.linalg.norm(coords - coords[i], axis=1)
        cutoff = tolerance * (radii + radii[i])
        for j in np.where((d > 0.4) & (d < cutoff))[0]:
            if j != i:
                adj[i].add(int(j))
                adj[int(j)].add(i)
    return adj


def bond_graph_distances(adj):
    """All-pairs shortest path length (#bonds) via BFS. Unreachable -> n+1."""
    n = len(adj)
    INF = n + 1
    dist = np.full((n, n), INF, dtype=int)
    for src in range(n):
        dist[src, src] = 0
        frontier = [src]
        d = 0
        seen = {src}
        while frontier:
            d += 1
            nxt = []
            for x in frontier:
                for nb in adj[x]:
                    if nb not in seen:
                        seen.add(nb)
                        dist[src, nb] = d
                        nxt.append(nb)
            frontier = nxt
    return dist


def moving_atoms(adj, b, c):
    """Indices that must rotate to change a dihedral about bond b-c: the
    connected component containing c once the b-c bond is removed (c stays on
    the axis, so including it is harmless). Raises if b-c lies in a ring."""
    seen = {c}
    stack = [c]
    while stack:
        x = stack.pop()
        for nb in adj[x]:
            if {x, nb} == {b, c}:      # never traverse the rotation bond
                continue
            if nb not in seen:
                seen.add(nb)
                stack.append(nb)
    if b in seen:
        raise ValueError(
            "bond %d-%d is inside a ring; it is not a free torsion" % (b, c))
    return seen


# --------------------------------------------------------------------------- #
# Geometry
# --------------------------------------------------------------------------- #

def dihedral_deg(p0, p1, p2, p3):
    """Signed dihedral p0-p1-p2-p3 in degrees, IUPAC convention, [-180,180)."""
    b1 = p1 - p0
    b2 = p2 - p1
    b3 = p3 - p2
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    b2u = b2 / np.linalg.norm(b2)
    m1 = np.cross(n1, b2u)
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return np.degrees(np.arctan2(y, x))


def rotation_matrix(axis_unit, theta):
    """Rodrigues rotation matrix for angle theta (rad) about a unit axis."""
    u = axis_unit
    c, s = np.cos(theta), np.sin(theta)
    ux, uy, uz = u
    return np.array([
        [c + ux * ux * (1 - c),      ux * uy * (1 - c) - uz * s, ux * uz * (1 - c) + uy * s],
        [uy * ux * (1 - c) + uz * s, c + uy * uy * (1 - c),      uy * uz * (1 - c) - ux * s],
        [uz * ux * (1 - c) - uy * s, uz * uy * (1 - c) + ux * s, c + uz * uz * (1 - c)],
    ])


def set_dihedral(coords, a, b, c, d, target_deg, movers):
    """Rotate `movers` about the b-c axis so that dihedral a-b-c-d == target."""
    current = dihedral_deg(coords[a], coords[b], coords[c], coords[d])
    # A right-hand rotation by +alpha about the B->C axis changes the dihedral
    # by -alpha, so to reach `target` we rotate by (current - target).
    delta = np.radians(current - target_deg)
    axis = coords[c] - coords[b]
    axis /= np.linalg.norm(axis)
    R = rotation_matrix(axis, delta)
    pivot = coords[c]
    idx = list(movers)
    coords[idx] = (coords[idx] - pivot) @ R.T + pivot
    return coords


# --------------------------------------------------------------------------- #
# Core generation
# --------------------------------------------------------------------------- #

def sample_values(lo, hi, step):
    """Absolute dihedral targets from lo to hi spaced by step. A full 360 span
    drops the duplicate endpoint (360 == 0)."""
    vals = []
    v = lo
    while v <= hi + 1e-9:
        vals.append(round(v, 6))
        v += step
    if abs((hi - lo) - 360.0) < 1e-9 and len(vals) > 1:
        vals.pop()
    return vals


class Torsion:
    def __init__(self, name, atom_idx, quad_names, values, movers):
        self.name = name
        self.atom_idx = atom_idx      # (a, b, c, d) indices
        self.quad_names = quad_names
        self.values = values
        self.movers = movers


def prepare_torsions(config_torsions, name_to_idx, adj):
    torsions = []
    for t in config_torsions:
        quad = t["atoms"]
        if len(quad) != 4:
            raise ValueError("torsion %r needs exactly 4 atom names, got %r"
                             % (t.get("name"), quad))
        try:
            a, b, c, d = (name_to_idx[nm] for nm in quad)
        except KeyError as e:
            raise ValueError("torsion %r: atom %s not found in PDB"
                             % (t.get("name"), e))
        if c not in adj[b]:
            raise ValueError("torsion %r: %s-%s are not bonded (cannot rotate "
                             "about a non-bond)" % (t.get("name"), quad[1], quad[2]))
        movers = moving_atoms(adj, b, c)
        lo, hi = t["range"]
        vals = sample_values(float(lo), float(hi), float(t["step"]))
        torsions.append(Torsion(t.get("name", "%s-%s-%s-%s" % tuple(quad)),
                                (a, b, c, d), quad, vals, movers))
    return torsions


def generate_rotamers(atoms, adj, torsions, cutoff, exclude_bonds, verbose=True):
    """Return (coords_all (R,n,3), dihedrals (R,T)). Discards self-clashers."""
    template = np.array([at.xyz for at in atoms])
    n = len(atoms)

    gd = bond_graph_distances(adj)
    iu, ju = np.triu_indices(n, k=1)
    check = gd[iu, ju] > exclude_bonds
    Iarr, Jarr = iu[check], ju[check]
    cutoff2 = cutoff * cutoff

    total = 1
    for t in torsions:
        total *= len(t.values)
    if verbose:
        print("[info] atoms=%d  torsions=%d  grid=%d combinations  "
              "checked-pairs=%d" % (n, len(torsions), total, len(Iarr)),
              file=sys.stderr)

    kept_coords, kept_dih = [], []
    for k, combo in enumerate(itertools.product(*[t.values for t in torsions])):
        coords = template.copy()
        for t, target in zip(torsions, combo):
            a, b, c, d = t.atom_idx
            set_dihedral(coords, a, b, c, d, target, t.movers)
        if len(Iarr):
            diff = coords[Iarr] - coords[Jarr]
            if np.einsum("ij,ij->i", diff, diff).min() < cutoff2:
                continue  # self-clash -> drop
        kept_coords.append(coords)
        kept_dih.append(combo)
        if verbose and total >= 20000 and (k + 1) % 20000 == 0:
            print("[info]   scanned %d/%d, kept %d"
                  % (k + 1, total, len(kept_coords)), file=sys.stderr)

    if verbose:
        print("[info] kept %d / %d rotamers after clash filtering"
              % (len(kept_coords), total), file=sys.stderr)
    if not kept_coords:
        return np.empty((0, n, 3)), np.empty((0, len(torsions)))
    return np.array(kept_coords), np.array(kept_dih, dtype=float)


# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #

def load_config(path):
    if yaml is None:
        raise RuntimeError("PyYAML is required to read the config file")
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    if "torsions" not in cfg or "clash" not in cfg:
        raise ValueError("config must have 'torsions' and 'clash' sections")
    clash = cfg["clash"]
    return cfg["torsions"], float(clash.get("cutoff", 2.5)), \
        int(clash.get("exclude_bonds", 3))


# --------------------------------------------------------------------------- #
# Self-test (engine verification, no external files)
# --------------------------------------------------------------------------- #

def _build_test_molecule():
    """A-B-C-D-E chain so that dihedral A-B-C-D is well defined and atom E
    (4 bonds from A) can be swung toward or away from A to create a clash."""
    specs = [
        ("A", "C", (-0.5, 1.40, 0.0)),
        ("B", "C", (0.0, 0.0, 0.0)),
        ("C", "C", (1.5, 0.0, 0.0)),
        ("D", "C", (2.0, 1.40, 0.0)),
        ("E", "C", (1.55, 2.75, 0.0)),
    ]
    atoms = [Atom(i + 1, nm, "LIG", "A", 1, el, xyz)
             for i, (nm, el, xyz) in enumerate(specs)]
    adj = [set() for _ in atoms]
    for i, j in [(0, 1), (1, 2), (2, 3), (3, 4)]:   # A-B-C-D-E
        adj[i].add(j)
        adj[j].add(i)
    return atoms, adj


def self_test():
    atoms, adj = _build_test_molecule()
    name_to_idx = {at.name: i for i, at in enumerate(atoms)}
    ok = True

    # --- bond-graph distances / exclusion ---
    gd = bond_graph_distances(adj)
    assert gd[0, 1] == 1 and gd[0, 2] == 2 and gd[0, 3] == 3 and gd[0, 4] == 4
    print("[test] bond-graph distances A..E = "
          "%d %d %d %d  -> with exclude_bonds=3 only A-E is checked  OK"
          % (gd[0, 1], gd[0, 2], gd[0, 3], gd[0, 4]))

    # --- dihedral is set to the requested absolute target ---
    tcfg = [{"name": "t1", "atoms": ["A", "B", "C", "D"],
             "range": [0, 360], "step": 30}]
    tors = prepare_torsions(tcfg, name_to_idx, adj)
    assert moving_atoms(adj, 1, 2) == {2, 3, 4}, "moving set must be {C,D,E}"
    coords, dih = generate_rotamers(atoms, adj, tors, cutoff=0.0,
                                    exclude_bonds=3, verbose=False)
    max_err = 0.0
    for r in range(coords.shape[0]):
        measured = dihedral_deg(coords[r, 0], coords[r, 1], coords[r, 2], coords[r, 3])
        want = dih[r, 0]
        err = abs(((measured - want + 180) % 360) - 180)
        max_err = max(max_err, err)
    assert max_err < 1e-6, "dihedral target not met, max err=%g" % max_err
    print("[test] set_dihedral reproduces absolute target (max err %.2e deg)  OK"
          % max_err)

    # --- A, B, C never move; D and E do ---
    assert np.allclose(coords[:, 0], coords[0, 0]) \
        and np.allclose(coords[:, 1], coords[0, 1]) \
        and np.allclose(coords[:, 2], coords[0, 2]), "upstream atoms moved!"
    assert not np.allclose(coords[:, 4], coords[0, 4]), "atom E should move"
    print("[test] upstream atoms (A,B,C) fixed, downstream (D,E) move  OK")

    # --- clash filtering matches an independent recomputation ---
    cutoff = 2.5
    coords_c, dih_c = generate_rotamers(atoms, adj, tors, cutoff=cutoff,
                                        exclude_bonds=3, verbose=False)
    # independent brute-force count of clash-free grid points
    grid = sample_values(0, 360, 30)
    indep_ok = 0
    for target in grid:
        cc = np.array([at.xyz for at in atoms])
        set_dihedral(cc, 0, 1, 2, 3, target, moving_atoms(adj, 1, 2))
        dAE = np.linalg.norm(cc[0] - cc[4])   # only A-E is checkable
        if dAE >= cutoff:
            indep_ok += 1
    assert coords_c.shape[0] == indep_ok, \
        "clash count mismatch: tool kept %d, brute force %d" % (coords_c.shape[0], indep_ok)
    # every kept rotamer must really be clash-free
    for r in range(coords_c.shape[0]):
        assert np.linalg.norm(coords_c[r, 0] - coords_c[r, 4]) >= cutoff
    assert 0 < coords_c.shape[0] < len(grid), \
        "test geometry should reject some but not all rotamers"
    print("[test] clash filter: kept %d/%d, matches brute force, and every kept "
          "rotamer has A-E >= %.1f A  OK" % (coords_c.shape[0], len(grid), cutoff))

    # --- ring bond is rejected ---
    ring_adj = [set(s) for s in adj]
    ring_adj[0].add(4)   # close A-E -> now B-C sits in a ring A-B-C-D-E-A
    ring_adj[4].add(0)
    try:
        moving_atoms(ring_adj, 1, 2)
        ok = False
        print("[test] ring detection FAILED (should have raised)")
    except ValueError:
        print("[test] ring bond correctly rejected as non-rotatable  OK")

    print("\nSELF-TEST %s" % ("PASSED" if ok else "FAILED"))
    return ok


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", help="template PDB (single molecule/residue)")
    ap.add_argument("--config", help="YAML config (torsions + clash)")
    ap.add_argument("--out-prefix", default="rotamer_library",
                    help="output prefix (writes <prefix>.pdb and <prefix>.npz)")
    ap.add_argument("--keep-hydrogens", action="store_true",
                    help="keep hydrogens (default: strip to heavy atoms)")
    ap.add_argument("--bond-tolerance", type=float, default=1.3,
                    help="distance bond-perception tolerance factor (default 1.3)")
    ap.add_argument("--self-test", action="store_true",
                    help="run the built-in engine verification and exit")
    args = ap.parse_args(argv)

    if args.self_test:
        return 0 if self_test() else 1

    if not args.pdb or not args.config:
        ap.error("--pdb and --config are required (or use --self-test)")

    atoms, conect = read_pdb(args.pdb, keep_hydrogens=args.keep_hydrogens)
    if not atoms:
        ap.error("no atoms read from %s" % args.pdb)
    names = [at.name for at in atoms]
    if len(set(names)) != len(names):
        dup = sorted({n for n in names if names.count(n) > 1})
        ap.error("atom names must be unique; duplicates: %s. Provide a single "
                 "residue/molecule as the template." % ", ".join(dup))
    name_to_idx = {n: i for i, n in enumerate(names)}

    adj = build_bonds(atoms, conect, tolerance=args.bond_tolerance)
    cfg_torsions, cutoff, exclude_bonds = load_config(args.config)
    torsions = prepare_torsions(cfg_torsions, name_to_idx, adj)

    coords_all, dihedrals = generate_rotamers(
        atoms, adj, torsions, cutoff=cutoff, exclude_bonds=exclude_bonds)

    if coords_all.shape[0] == 0:
        print("[warn] no rotamers survived clash filtering; nothing written",
              file=sys.stderr)
        return 2

    tnames = [t.name for t in torsions]
    pdb_out = args.out_prefix + ".pdb"
    npz_out = args.out_prefix + ".npz"
    write_multimodel_pdb(pdb_out, atoms, coords_all, dihedrals, tnames)
    np.savez(npz_out,
             coords=coords_all,
             atom_names=np.array(names),
             elements=np.array([at.element for at in atoms]),
             resname=atoms[0].resname,
             dihedrals=dihedrals,
             torsion_names=np.array(tnames))
    print("[done] wrote %s (%d models) and %s"
          % (pdb_out, coords_all.shape[0], npz_out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
