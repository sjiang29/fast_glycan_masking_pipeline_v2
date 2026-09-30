#!/usr/bin/env python3
"""Generate a reusable glycan conformer library from a Rosetta glycoprotein PDB.

N-linked Man5 conformer-library generator:
  * follows LINK records from a selected ASN to recover the glycan tree;
  * samples the ASN->first-GlcNAc attachment torsion CG-ND2-C1-O5;
  * identifies phi/psi and omega (for 1->6) glycosidic torsions;
  * samples torsions from linkage-specific distributions in YAML;
  * rotates the complete downstream glycan subtree;
  * keeps attachment bond lengths/angles fixed while sampling only the torsion;
  * rejects glycan self-clashes;
  * writes a glycan-only multi-model PDB and NPZ library.

The source ASN CB/CG/ND2 frame is retained so every sampled conformer can be
rigidly transferred to a target ASN by placement/placer.py.
Angle distributions are intentionally data/config driven because published angle signs and
ranges depend on atom ordering and convention.
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple

import numpy as np
import yaml

from .rotamer import (
    Atom,
    bond_graph_distances,
    build_bonds,
    dihedral_deg,
    infer_element,
    rotation_matrix,
    write_multimodel_pdb,
)

AtomKey = Tuple[str, int, str]  # chain, residue number, atom name
ResidueKey = Tuple[str, int]
SUGAR_NAMES = {"GLC", "MAN", "BMA", "NAG", "NDG", "GLCNAC"}


@dataclass(frozen=True)
class Link:
    atom1: AtomKey
    resname1: str
    atom2: AtomKey
    resname2: str


@dataclass
class GlycoTorsion:
    name: str
    kind: str
    linkage: str
    atoms: Tuple[int, int, int, int]
    movers: np.ndarray
    template_angle: float


def parse_atom_line(line: str) -> Atom:
    name = line[12:16].strip()
    element = line[76:78].strip().upper()
    if not element or not element.isalpha():
        element = infer_element(name)
    return Atom(
        int(line[6:11]), name, line[17:20].strip(),
        line[21:22].strip() or "A", int(line[22:26]), element,
        (float(line[30:38]), float(line[38:46]), float(line[46:54])),
    )


def parse_link_line(line: str) -> Link:
    a1 = line[12:16].strip(); r1 = line[17:20].strip()
    c1 = line[21:22].strip() or "A"; n1 = int(line[22:26])
    a2 = line[42:46].strip(); r2 = line[47:50].strip()
    c2 = line[51:52].strip() or "A"; n2 = int(line[52:56])
    return Link((c1, n1, a1), r1, (c2, n2, a2), r2)


def read_rosetta_pdb(path: str, keep_hydrogens: bool = False):
    atoms: List[Atom] = []
    links: List[Link] = []
    hetnam: Dict[ResidueKey, str] = {}
    for line in open(path):
        if line.startswith(("ATOM  ", "HETATM")):
            at = parse_atom_line(line)
            if not keep_hydrogens and at.element == "H":
                continue
            atoms.append(at)
        elif line.startswith("LINK"):
            links.append(parse_link_line(line))
        elif line.startswith("HETNAM"):
            # Rosetta example: HETNAM     Man A 273  ->4)-alpha-D-Manp
            m = re.match(r"HETNAM\s+(\S+)\s+(\S)\s*(\d+)\s+(.*)", line)
            if m:
                hetnam[(m.group(2), int(m.group(3)))] = m.group(4).strip()
    if not atoms:
        raise ValueError(f"No atoms read from {path}")
    return atoms, links, hetnam


def is_sugar(resname: str) -> bool:
    return resname.upper() in SUGAR_NAMES or resname.upper().startswith(("GLC", "MAN"))


def find_attached_glycan(links: Sequence[Link], chain: str, resid: int):
    """Return root, glycan residues, oriented sugar links, and ASN-root LINK."""
    root: ResidueKey | None = None
    attachment_link: Link | None = None
    sugar_edges: Dict[ResidueKey, List[Tuple[ResidueKey, Link]]] = {}
    for lk in links:
        r1 = lk.atom1[:2]; r2 = lk.atom2[:2]
        if r1 == (chain, resid) and lk.atom1[2] == "ND2" and is_sugar(lk.resname2):
            root = r2
            # Orient attachment as ASN ND2 -> root-sugar C1.
            attachment_link = lk
        elif r2 == (chain, resid) and lk.atom2[2] == "ND2" and is_sugar(lk.resname1):
            root = r1
            attachment_link = Link(lk.atom2, lk.resname2, lk.atom1, lk.resname1)
        if is_sugar(lk.resname1) and is_sugar(lk.resname2):
            # Rosetta LINK records in this file list parent O first, child C1 second.
            if lk.atom1[2].startswith("O") and lk.atom2[2] == "C1":
                sugar_edges.setdefault(r1, []).append((r2, lk))
            elif lk.atom2[2].startswith("O") and lk.atom1[2] == "C1":
                rev = Link(lk.atom2, lk.resname2, lk.atom1, lk.resname1)
                sugar_edges.setdefault(r2, []).append((r1, rev))
    if root is None or attachment_link is None:
        raise ValueError(f"No N-linked glycan found at {chain}{resid}")
    if attachment_link.atom1[2] != "ND2" or attachment_link.atom2[2] != "C1":
        raise ValueError(
            f"Expected ASN ND2 -> root C1 attachment at {chain}{resid}; "
            f"found {attachment_link.atom1[2]} -> {attachment_link.atom2[2]}"
        )
    residues: Set[ResidueKey] = set()
    oriented: List[Tuple[ResidueKey, ResidueKey, Link]] = []
    stack = [root]
    while stack:
        parent = stack.pop()
        if parent in residues:
            continue
        residues.add(parent)
        for child, lk in sugar_edges.get(parent, []):
            oriented.append((parent, child, lk))
            stack.append(child)
    return root, residues, oriented, attachment_link


def component_after_cut(adj: Sequence[Set[int]], start: int, cut_a: int, cut_b: int) -> Set[int]:
    seen = {start}; stack = [start]
    while stack:
        x = stack.pop()
        for nb in adj[x]:
            if {x, nb} == {cut_a, cut_b}:
                continue
            if nb not in seen:
                seen.add(nb); stack.append(nb)
    return seen


def normalize_angle(x: float) -> float:
    return ((x + 180.0) % 360.0) - 180.0


def set_dihedral(coords: np.ndarray, quad: Tuple[int, int, int, int], target: float,
                 movers: np.ndarray) -> None:
    a, b, c, d = quad
    current = dihedral_deg(coords[a], coords[b], coords[c], coords[d])
    delta = math.radians(current - target)
    axis = coords[c] - coords[b]
    norm = np.linalg.norm(axis)
    if norm < 1e-8:
        raise ValueError("Degenerate torsion axis")
    R = rotation_matrix(axis / norm, delta)
    pivot = coords[c]
    coords[movers] = (coords[movers] - pivot) @ R.T + pivot


def set_attachment_dihedral(coords: np.ndarray, asn_cg: np.ndarray, asn_nd2: np.ndarray,
                            root_c1_idx: int, root_o5_idx: int, target: float) -> None:
    """Set CG(ASN)-ND2(ASN)-C1(root)-O5(root) by rotating the whole glycan.

    Only the torsion changes.  The ASN coordinates and the ND2-C1 bond geometry
    are not rebuilt or stretched; C1 lies on the rotation axis and remains fixed.
    """
    current = dihedral_deg(asn_cg, asn_nd2, coords[root_c1_idx], coords[root_o5_idx])
    delta = math.radians(current - target)
    axis = coords[root_c1_idx] - asn_nd2
    norm = np.linalg.norm(axis)
    if norm < 1e-8:
        raise ValueError("Degenerate ASN-GlcNAc attachment axis")
    R = rotation_matrix(axis / norm, delta)
    pivot = coords[root_c1_idx]
    coords[:] = (coords - pivot) @ R.T + pivot


def infer_anomer(child: ResidueKey, hetnam: Mapping[ResidueKey, str]) -> str:
    text = hetnam.get(child, "").lower()
    if "alpha" in text:
        return "alpha"
    if "beta" in text:
        return "beta"
    return "unknown"


def build_torsions(atoms: Sequence[Atom], adj: Sequence[Set[int]], oriented_links,
                   hetnam: Mapping[ResidueKey, str]) -> List[GlycoTorsion]:
    key_to_idx: Dict[AtomKey, int] = {(a.chain, a.resseq, a.name): i for i, a in enumerate(atoms)}
    coords = np.array([a.xyz for a in atoms])
    torsions: List[GlycoTorsion] = []

    def idx(res: ResidueKey, atom: str) -> int:
        try:
            return key_to_idx[(res[0], res[1], atom)]
        except KeyError as exc:
            raise ValueError(f"Missing atom {res[0]}:{res[1]}:{atom}") from exc

    for parent, child, lk in oriented_links:
        acceptor_o = lk.atom1[2]
        m = re.fullmatch(r"O(\d+)", acceptor_o)
        if not m:
            raise ValueError(f"Cannot infer linkage position from {acceptor_o}")
        pos = int(m.group(1))
        anomer = infer_anomer(child, hetnam)
        linkage = f"{anomer}1-{pos}"

        # Reversed atom order preserves the same dihedral value but places the
        # downstream child side on atom c, matching set_dihedral's rotation sign.
        phi_quad = (idx(parent, f"C{pos}"), idx(parent, acceptor_o),
                    idx(child, "C1"), idx(child, "O5"))
        phi_movers = np.array(sorted(component_after_cut(adj, phi_quad[2], phi_quad[1], phi_quad[2])), int)
        torsions.append(GlycoTorsion(
            f"{parent[1]}_{child[1]}_phi", "phi", linkage, phi_quad,
            phi_movers, dihedral_deg(*(coords[i] for i in phi_quad))))

        prev = "C5" if pos == 6 else f"C{pos-1}"
        psi_quad = (idx(parent, prev), idx(parent, f"C{pos}"),
                    idx(parent, acceptor_o), idx(child, "C1"))
        psi_movers = np.array(sorted(component_after_cut(adj, psi_quad[2], psi_quad[1], psi_quad[2])), int)
        torsions.append(GlycoTorsion(
            f"{parent[1]}_{child[1]}_psi", "psi", linkage, psi_quad,
            psi_movers, dihedral_deg(*(coords[i] for i in psi_quad))))

        if pos == 6:
            omega_quad = (idx(parent, "C4"), idx(parent, "C5"),
                          idx(parent, "C6"), idx(parent, "O6"))
            omega_movers = np.array(sorted(component_after_cut(adj, omega_quad[2], omega_quad[1], omega_quad[2])), int)
            torsions.append(GlycoTorsion(
                f"{parent[1]}_{child[1]}_omega", "omega", linkage, omega_quad,
                omega_movers, dihedral_deg(*(coords[i] for i in omega_quad))))
    return torsions


def circular_normal(rng: np.random.Generator, mean: float, sd: float) -> float:
    return normalize_angle(float(rng.normal(mean, sd)))


def sample_target(rng: np.random.Generator, torsion: GlycoTorsion, cfg: Mapping) -> float:
    link_cfg = cfg.get("linkages", {}).get(torsion.linkage, {})
    spec = link_cfg.get(torsion.kind, cfg.get("default", {}).get(torsion.kind, {}))
    mode = spec.get("mode", "template_jitter")
    if mode == "normal":
        return circular_normal(rng, float(spec["mean"]), float(spec["sd"]))
    if mode == "template_jitter":
        return circular_normal(rng, torsion.template_angle + float(spec.get("offset", 0.0)),
                               float(spec.get("sd", 15.0)))
    if mode == "discrete":
        vals = np.asarray(spec["values"], float)
        weights = np.asarray(spec.get("weights", np.ones(len(vals))), float)
        weights /= weights.sum()
        return float(rng.choice(vals, p=weights))
    raise ValueError(f"Unknown sampling mode {mode!r} for {torsion.linkage}/{torsion.kind}")


def sample_attachment_target(rng: np.random.Generator, template_angle: float,
                             cfg: Mapping) -> float:
    """Sample CG-ND2-C1-O5 from the YAML ``attachment`` section."""
    spec = cfg.get("attachment", {})
    mode = spec.get("mode", "template_jitter")
    if mode == "normal":
        return circular_normal(rng, float(spec["mean"]), float(spec["sd"]))
    if mode == "template_jitter":
        return circular_normal(
            rng, template_angle + float(spec.get("offset", 0.0)),
            float(spec.get("sd", 20.0))
        )
    if mode == "discrete":
        vals = np.asarray(spec["values"], float)
        weights = np.asarray(spec.get("weights", np.ones(len(vals))), float)
        weights /= weights.sum()
        return float(rng.choice(vals, p=weights))
    if mode == "fixed":
        return normalize_angle(float(spec.get("value", template_angle)))
    raise ValueError(f"Unknown sampling mode {mode!r} for attachment torsion")


def self_clash_pairs(atoms: Sequence[Atom], adj: Sequence[Set[int]], exclude_bonds: int):
    gd = bond_graph_distances(adj)
    i, j = np.triu_indices(len(atoms), 1)
    mask = gd[i, j] > exclude_bonds
    return i[mask], j[mask]


def generate(atoms: Sequence[Atom], adj: Sequence[Set[int]], torsions: Sequence[GlycoTorsion],
             cfg: Mapping, n_conformers: int, max_attempts: int, seed: int,
             clash_cutoff: float, exclude_bonds: int, asn_cg: np.ndarray,
             asn_nd2: np.ndarray, root_c1_idx: int, root_o5_idx: int,
             attachment_template_angle: float):
    rng = np.random.default_rng(seed)
    template = np.array([a.xyz for a in atoms])
    ii, jj = self_clash_pairs(atoms, adj, exclude_bonds)
    cutoff2 = clash_cutoff ** 2
    kept, angles = [], []
    attempts = 0
    while len(kept) < n_conformers and attempts < max_attempts:
        attempts += 1
        xyz = template.copy(); sampled = []

        # Level 1: sample ASN -> first GlcNAc orientation around ND2-C1.
        attachment_target = sample_attachment_target(rng, attachment_template_angle, cfg)
        set_attachment_dihedral(
            xyz, asn_cg, asn_nd2, root_c1_idx, root_o5_idx, attachment_target
        )
        sampled.append(attachment_target)

        # Level 2: sample downstream glycan-glycan phi/psi/omega torsions.
        for tor in torsions:
            target = sample_target(rng, tor, cfg)
            set_dihedral(xyz, tor.atoms, target, tor.movers)
            sampled.append(target)
        if len(ii):
            d = xyz[ii] - xyz[jj]
            if np.any(np.einsum("ij,ij->i", d, d) < cutoff2):
                continue
        kept.append(xyz); angles.append(sampled)
    n_sampled_torsions = 1 + len(torsions)
    if not kept:
        return (np.empty((0, len(atoms), 3)),
                np.empty((0, n_sampled_torsions)), attempts)
    return np.asarray(kept), np.asarray(angles), attempts


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--resid", required=True, type=int, help="ASN residue carrying the reference glycan")
    ap.add_argument("--config", required=True)
    ap.add_argument(
        "--n-conformers",
        type=int,
        default=10000,
        help="Number of clash-free glycan conformers to generate (default: 10000)",
    )
    # Backward-compatible alias used by earlier versions.
    ap.add_argument(
        "--n-models",
        dest="n_conformers",
        type=int,
        help=argparse.SUPPRESS,
    )
    ap.add_argument("--max-attempts", type=int, default=100000)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--clash-cutoff", type=float, default=2.0)
    ap.add_argument("--exclude-bonds", type=int, default=3)
    ap.add_argument("--out-prefix", default="Man5_library")
    ap.add_argument("--keep-hydrogens", action="store_true")
    ap.add_argument("--report-only", action="store_true")
    args = ap.parse_args(argv)

    all_atoms, links, hetnam = read_rosetta_pdb(args.pdb, args.keep_hydrogens)
    root, residues, oriented, attachment_link = find_attached_glycan(
        links, args.chain, args.resid
    )
    atoms = [a for a in all_atoms if (a.chain, a.resseq) in residues]
    # Add explicit LINK connectivity to distance-perceived intra-residue bonds.
    adj = build_bonds(atoms, conect=None, tolerance=1.25)
    key_to_idx = {(a.chain, a.resseq, a.name): i for i, a in enumerate(atoms)}
    for _, _, lk in oriented:
        i = key_to_idx[lk.atom1]; j = key_to_idx[lk.atom2]
        adj[i].add(j); adj[j].add(i)
    torsions = build_torsions(atoms, adj, oriented, hetnam)

    # Source ASN coordinates are external to the glycan-only atom array.
    all_key_to_xyz = {(a.chain, a.resseq, a.name): np.asarray(a.xyz, float) for a in all_atoms}
    try:
        asn_cg = all_key_to_xyz[(args.chain, args.resid, "CG")]
        asn_nd2 = all_key_to_xyz[(args.chain, args.resid, "ND2")]
        root_c1_idx = key_to_idx[(root[0], root[1], "C1")]
        root_o5_idx = key_to_idx[(root[0], root[1], "O5")]
    except KeyError as exc:
        raise ValueError(
            f"Cannot build ASN-GlcNAc attachment torsion at {args.chain}{args.resid}; "
            f"missing atom {exc.args[0]}"
        ) from exc
    template_coords = np.array([a.xyz for a in atoms])
    attachment_template_angle = dihedral_deg(
        asn_cg, asn_nd2, template_coords[root_c1_idx], template_coords[root_o5_idx]
    )

    print(f"[info] root sugar: {root[0]}{root[1]}; residues={len(residues)} atoms={len(atoms)}")
    print(
        f"[attachment] {args.chain}{args.resid} CG-ND2-{root[0]}{root[1]} C1-O5 "
        f"template={attachment_template_angle:8.2f} deg"
    )
    for t in torsions:
        print(f"[torsion] {t.name:18s} {t.linkage:10s} {t.kind:5s} template={t.template_angle:8.2f}")
    if args.report_only:
        return 0

    cfg = yaml.safe_load(open(args.config)) or {}
    coords, angle_values, attempts = generate(
        atoms, adj, torsions, cfg, args.n_conformers, args.max_attempts, args.seed,
        args.clash_cutoff, args.exclude_bonds, asn_cg, asn_nd2,
        root_c1_idx, root_o5_idx, attachment_template_angle)
    if len(coords) < args.n_conformers:
        print(
            f"[warn] generated {len(coords)}/{args.n_conformers} conformers "
            f"after {attempts} attempts",
            file=sys.stderr,
        )
    if len(coords) == 0:
        return 2
    out = Path(args.out_prefix)
    all_torsion_names = ["asn_glcnac_attachment"] + [t.name for t in torsions]
    write_multimodel_pdb(str(out) + ".pdb", atoms, coords, angle_values,
                         all_torsion_names)
    # Store the source ASN side-chain frame so the library can be rigidly
    # transferred to any other ASN without rebuilding the attachment geometry.
    asn_atoms = {(a.chain, a.resseq, a.name): a.xyz for a in all_atoms}
    try:
        attachment_frame = np.array([
            asn_atoms[(args.chain, args.resid, "CB")],
            asn_atoms[(args.chain, args.resid, "CG")],
            asn_atoms[(args.chain, args.resid, "ND2")],
        ], dtype=float)
    except KeyError as exc:
        raise ValueError(f"Source ASN {args.chain}{args.resid} lacks CB/CG/ND2") from exc

    np.savez(str(out) + ".npz", coords=coords,
             atom_names=np.array([a.name for a in atoms]),
             residue_numbers=np.array([a.resseq for a in atoms]),
             residue_names=np.array([a.resname for a in atoms]),
             chains=np.array([a.chain for a in atoms]),
             elements=np.array([a.element for a in atoms]),
             attachment_frame=attachment_frame,
             root_residue_number=int(root[1]),
             attachment_parent_atom=np.array(attachment_link.atom1[2]),
             attachment_child_atom=np.array(attachment_link.atom2[2]),
             attachment_torsion_name=np.array("asn_glcnac_attachment"),
             attachment_template_angle=np.array(attachment_template_angle),
             sampled_attachment_angles=angle_values[:, 0],
             linkage_parent_resid=np.array([p[1] for p, c, lk in oriented], dtype=int),
             linkage_child_resid=np.array([c[1] for p, c, lk in oriented], dtype=int),
             linkage_parent_atom=np.array([lk.atom1[2] for p, c, lk in oriented]),
             linkage_child_atom=np.array([lk.atom2[2] for p, c, lk in oriented]),
             torsion_names=np.array(all_torsion_names),
             linkage_types=np.array(["ASN-N-linked"] + [t.linkage for t in torsions]),
             torsion_kinds=np.array(["attachment"] + [t.kind for t in torsions]),
             template_angles=np.array(
                 [attachment_template_angle] + [t.template_angle for t in torsions]
             ),
             sampled_angles=angle_values,
             source_pdb=args.pdb, source_chain=args.chain, source_resid=args.resid,
             seed=args.seed)
    print(f"[done] wrote {out}.pdb and {out}.npz; {len(coords)} conformers from {attempts} attempts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
