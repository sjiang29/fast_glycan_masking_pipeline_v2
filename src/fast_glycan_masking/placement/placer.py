#!/usr/bin/env python3
"""Place a reusable Man5 conformer library at one or more N-glycosylation sites.

The input protein must already contain valid N-X-S/T sequons at every requested
site. Use fast_glycan_masking.pipeline when new sites need to be created by Rosetta
FastDesign first.

For each output model this script independently chooses one library conformer per
site, aligns it by the ASN CB/CG/ND2 frame, checks protein-glycan clashes and
between-glycan clashes, and writes one complete multi-glycosylated PDB.

By default all carbohydrate residues already present in the input PDB are removed
before placement. This mirrors the old Rosetta workflow with strip_existing="1":
native sites and newly engineered sites are all rebuilt from the reusable library.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np

from ..library_builder.from_template.rotamer import Atom, infer_element

WATER_NAMES = {"HOH", "WAT", "H2O", "DOD", "SOL", "TIP", "TIP3", "TIP4", "TIP5", "SPC"}
SUGAR_NAMES = {
    "GLC", "MAN", "BMA", "NAG", "NDG", "GLCNAC", "FUC", "FUL", "GAL", "SIA",
    "NEU", "NAN", "A2G", "BGC", "G6D", "XYS", "XYP", "ARA", "RIB"
}
AA3_TO_1 = {
    "ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C", "GLN":"Q",
    "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I", "LEU":"L", "LYS":"K",
    "MET":"M", "PHE":"F", "PRO":"P", "SER":"S", "THR":"T", "TRP":"W",
    "TYR":"Y", "VAL":"V"
}
VDW = {"H":1.20, "C":1.70, "N":1.55, "O":1.52, "S":1.80, "P":1.80,
       "F":1.47, "CL":1.75, "BR":1.85, "I":1.98, "SE":1.90}


@dataclass(frozen=True)
class Site:
    chain: str
    resid: int

    @property
    def label(self) -> str:
        return f"{self.chain}{self.resid}"


@dataclass
class PDBRecord:
    atom: Atom
    record: str
    altloc: str = ""
    occupancy: float = 1.0
    bfactor: float = 0.0


def parse_site(text: str) -> Site:
    m = re.fullmatch(r"\s*([^:\s]+)\s*:\s*(-?\d+)\s*", text)
    if not m:
        raise argparse.ArgumentTypeError(f"Invalid site {text!r}; use CHAIN:RESID, e.g. A:200")
    return Site(m.group(1), int(m.group(2)))


def parse_atom_line(line: str) -> PDBRecord:
    name = line[12:16].strip()
    element = line[76:78].strip().upper()
    if not element or not element.isalpha():
        element = infer_element(name)
    chain = line[21:22].strip() or "A"
    try: serial = int(line[6:11])
    except ValueError: serial = 0
    try: resid = int(line[22:26])
    except ValueError: resid = 0
    try: occ = float(line[54:60])
    except ValueError: occ = 1.0
    try: b = float(line[60:66])
    except ValueError: b = 0.0
    at = Atom(serial, name, line[17:20].strip(), chain, resid, element,
              (float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return PDBRecord(at, line[:6].strip() or "ATOM", line[16:17].strip(), occ, b)


def read_structure(path: str, keep_hydrogens: bool = False, keep_water: bool = False,
                   keep_existing_glycans: bool = False) -> List[PDBRecord]:
    records: List[PDBRecord] = []
    for line in open(path):
        rec = line[:6].strip()
        if rec == "ENDMDL":
            break
        if rec not in {"ATOM", "HETATM"}:
            continue
        r = parse_atom_line(line)
        rn = r.atom.resname.upper()
        if not keep_hydrogens and r.atom.element == "H":
            continue
        if not keep_water and rn in WATER_NAMES:
            continue
        if not keep_existing_glycans and rn in SUGAR_NAMES:
            continue
        if r.altloc not in {"", "A"}:
            continue
        records.append(r)
    if not records:
        raise ValueError(f"No atoms read from {path}")
    return records


def residue_name(records: Sequence[PDBRecord], chain: str, resid: int) -> str | None:
    names = {r.atom.resname.upper() for r in records
             if r.atom.chain == chain and r.atom.resseq == resid}
    if not names:
        return None
    if len(names) > 1:
        raise ValueError(f"Multiple residue names at {chain}{resid}: {sorted(names)}")
    return next(iter(names))


def validate_sequon(records: Sequence[PDBRecord], site: Site) -> Tuple[bool, str]:
    r0 = residue_name(records, site.chain, site.resid)
    r1 = residue_name(records, site.chain, site.resid + 1)
    r2 = residue_name(records, site.chain, site.resid + 2)
    if r0 is None or r1 is None or r2 is None:
        return False, f"{site.label}: residues i, i+1, and i+2 must all exist in the PDB"
    aa0, aa1, aa2 = (AA3_TO_1.get(x, "?") for x in (r0, r1, r2))
    ok = aa0 == "N" and aa1 != "P" and aa2 in {"S", "T"}
    return ok, f"{site.label}: {aa0}-{aa1}-{aa2} ({r0}-{r1}-{r2})"


def kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = Qc - R @ Pc
    return R, t


def target_asn_frame(records: Sequence[PDBRecord], site: Site) -> np.ndarray:
    found = {r.atom.name: r.atom.xyz for r in records
             if r.atom.chain == site.chain and r.atom.resseq == site.resid}
    for name in ("CB", "CG", "ND2"):
        if name not in found:
            raise ValueError(f"Target {site.label} lacks ASN attachment atom {name}")
    return np.array([found["CB"], found["CG"], found["ND2"]], float)


def load_library(path: str) -> Dict[str, np.ndarray]:
    """Load and validate the NPZ fields required for glycan placement."""
    with np.load(path, allow_pickle=False) as z:
        required = {"coords", "atom_names", "residue_numbers", "residue_names",
                    "chains", "elements", "attachment_frame"}
        missing = required.difference(z.files)
        if missing:
            raise ValueError(f"Library {path} lacks fields {sorted(missing)}")
        return {k: z[k].copy() for k in z.files}


def vdw_radius(element: str) -> float:
    return VDW.get(str(element).upper(), 1.70)


def clash_between(a_xyz: np.ndarray, a_elem: Sequence[str], b_xyz: np.ndarray,
                  b_elem: Sequence[str], mode: str, cutoff: float, scale: float) -> bool:
    d2 = np.sum((a_xyz[:, None, :] - b_xyz[None, :, :]) ** 2, axis=-1)
    if mode == "distance":
        return bool(np.any(d2 < cutoff * cutoff))
    ar = np.array([vdw_radius(e) for e in a_elem])
    br = np.array([vdw_radius(e) for e in b_elem])
    threshold = scale * (ar[:, None] + br[None, :])
    return bool(np.any(d2 < threshold ** 2))


def protein_clash_mask(placed: np.ndarray, gly_elements: Sequence[str],
                       protein: Sequence[PDBRecord], all_sites: Sequence[Site],
                       mode: str, cutoff: float, scale: float,
                       exclude_target_residues: bool = True,
                       chunk_size: int = 32) -> np.ndarray:
    site_keys = {(s.chain, s.resid) for s in all_sites}
    prot = [r for r in protein if not (exclude_target_residues and
            (r.atom.chain, r.atom.resseq) in site_keys)]
    if not prot:
        return np.ones(len(placed), dtype=bool)
    pxyz = np.array([r.atom.xyz for r in prot])
    pelem = np.array([r.atom.element for r in prot])
    pr = np.array([vdw_radius(e) for e in pelem])
    gr = np.array([vdw_radius(e) for e in gly_elements])
    keep = np.ones(len(placed), dtype=bool)
    for start in range(0, len(placed), chunk_size):
        stop = min(start + chunk_size, len(placed))
        xyz = placed[start:stop]
        d2 = np.sum((xyz[:, :, None, :] - pxyz[None, None, :, :]) ** 2, axis=-1)
        if mode == "distance":
            clash = np.any(d2 < cutoff * cutoff, axis=(1, 2))
        else:
            threshold = scale * (gr[:, None] + pr[None, :])
            clash = np.any(d2 < threshold[None, :, :] ** 2, axis=(1, 2))
        keep[start:stop] = ~clash
    return keep


def _format_atom_name(name: str) -> str:
    return name[:4] if len(name) >= 4 else " " + name.ljust(3)


def format_atom(record: str, serial: int, name: str, resname: str, chain: str,
                resid: int, xyz: Sequence[float], element: str,
                occupancy: float = 1.0, bfactor: float = 0.0) -> str:
    x, y, z = xyz
    return (f"{record:<6}{serial:5d} {_format_atom_name(name)} {resname[:3]:>3s} "
            f"{chain[:1]}{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{occupancy:6.2f}{bfactor:6.2f}          {element:>2s}\n")


def unique_residue_order(residue_numbers: Sequence[int]) -> List[int]:
    seen: List[int] = []
    for x in residue_numbers:
        x = int(x)
        if x not in seen:
            seen.append(x)
    return seen


def build_numbering(protein: Sequence[PDBRecord], sites: Sequence[Site],
                    lib_resids: np.ndarray, glycan_chain: str | None,
                    residue_start: int | None) -> Dict[Site, Tuple[str, Dict[int, int]]]:
    order = unique_residue_order(lib_resids)
    used_max: Dict[str, int] = {}
    for r in protein:
        used_max[r.atom.chain] = max(used_max.get(r.atom.chain, 0), r.atom.resseq)
    result = {}
    next_override = residue_start
    for site in sites:
        chain = glycan_chain or site.chain
        start = next_override if next_override is not None else used_max.get(chain, 0) + 1
        mapping = {old: start + i for i, old in enumerate(order)}
        used_max[chain] = max(mapping.values())
        if next_override is not None:
            next_override = max(mapping.values()) + 1
        result[site] = (chain, mapping)
    return result


def write_model(path: Path, protein: Sequence[PDBRecord], placements: Mapping[Site, np.ndarray],
                lib: Mapping[str, np.ndarray], numbering: Mapping[Site, Tuple[str, Dict[int,int]]]) -> None:
    atom_names = lib["atom_names"]; resids = lib["residue_numbers"]
    resnames = lib["residue_names"]; elements = lib["elements"]
    root_old = int(resids[0])
    with open(path, "w") as fh:
        fh.write("REMARK  Fast multi-site Man5 conformer placement\n")
        fh.write("REMARK  Sites: " + ", ".join(s.label for s in placements) + "\n")
        serial = 1
        for r in protein:
            a = r.atom
            fh.write(format_atom(r.record, serial, a.name, a.resname, a.chain, a.resseq,
                                 a.xyz, a.element, r.occupancy, r.bfactor))
            serial += 1
        for site, xyz_all in placements.items():
            out_chain, rmap = numbering[site]
            for i, xyz in enumerate(xyz_all):
                fh.write(format_atom("HETATM", serial, str(atom_names[i]), str(resnames[i]),
                                     out_chain, rmap[int(resids[i])], xyz, str(elements[i])))
                serial += 1
        # LINK records
        for site in placements:
            out_chain, rmap = numbering[site]
            root_idx = int(np.where(resids == root_old)[0][0])
            root_name = str(resnames[root_idx])[:3]
            fh.write(f"LINK         ND2 ASN {site.chain}{site.resid:4d}                 C1 {root_name:>3s} {out_chain}{rmap[root_old]:4d}\n")
            if "linkage_parent_resid" in lib:
                for p, c, pa, ca in zip(lib["linkage_parent_resid"], lib["linkage_child_resid"],
                                        lib["linkage_parent_atom"], lib["linkage_child_atom"]):
                    p, c = int(p), int(c)
                    pi = int(np.where(resids == p)[0][0]); ci = int(np.where(resids == c)[0][0])
                    fh.write(f"LINK        {str(pa):>4s} {str(resnames[pi])[:3]:>3s} {out_chain}{rmap[p]:4d}                "
                             f"{str(ca):>4s} {str(resnames[ci])[:3]:>3s} {out_chain}{rmap[c]:4d}\n")
        fh.write("END\n")



def glycan_direction_vectors(placed: np.ndarray, target_frame: np.ndarray) -> np.ndarray:
    """Represent each placed glycan by the unit vector ND2 -> glycan geometric center."""
    """Unit vectors from target ASN ND2 to each glycan geometric center."""
    nd2 = np.asarray(target_frame, float)[2]
    centers = np.asarray(placed, float).mean(axis=1)
    vec = centers - nd2[None, :]
    norms = np.linalg.norm(vec, axis=1)
    good = norms > 1e-12
    out = np.zeros_like(vec)
    out[good] = vec[good] / norms[good, None]
    return out

def farthest_point_order(
    indices: np.ndarray,
    directions: np.ndarray,
    rng: np.random.Generator,
    n_select: int | None = None,
) -> np.ndarray:

    indices = np.asarray(indices, dtype=int)

    if n_select is None:
        n_select = len(indices)

    n_select = min(n_select, len(indices))

    if n_select <= 0:
        return np.array([], dtype=int)

    if len(indices) <= 1:
        return indices.copy()

    cand_dirs = directions[indices]

    # Random starting conformer
    first_pos = int(rng.integers(len(indices)))

    selected_pos = [first_pos]

    selected_mask = np.zeros(len(indices), dtype=bool)
    selected_mask[first_pos] = True

    # Similarity of every candidate to its most similar
    # already-selected conformer.
    max_dot = cand_dirs @ cand_dirs[first_pos]

    while len(selected_pos) < n_select:

        scores = max_dot.copy()
        scores[selected_mask] = np.inf

        # Lowest maximum similarity =
        # farthest from the current selected set.
        next_pos = int(np.argmin(scores))

        selected_pos.append(next_pos)
        selected_mask[next_pos] = True

        # Only one new vector comparison is required.
        new_dot = cand_dirs @ cand_dirs[next_pos]
        max_dot = np.maximum(max_dot, new_dot)

    return indices[np.asarray(selected_pos, dtype=int)]


def run(args: argparse.Namespace) -> int:
    sites: List[Site] = list(dict.fromkeys(args.site))
    if not sites:
        raise ValueError("At least one --site CHAIN:RESID is required")
    lib = load_library(args.library)
    protein = read_structure(args.protein, args.keep_hydrogens, args.keep_water,
                             args.keep_existing_glycans)

    if not args.skip_sequon_validation:
        errors = []
        for site in sites:
            ok, msg = validate_sequon(protein, site)
            print(f"[sequon] {msg}")
            if not ok:
                errors.append(msg)
        if errors:
            raise ValueError("Invalid N-glycosylation sequon(s). Run FastDesign first:\n  " + "\n  ".join(errors))

    coords = np.asarray(lib["coords"], float)
    source_frame = np.asarray(lib["attachment_frame"], float)
    site_data = {}
    for site in sites:
        target = target_asn_frame(protein, site)
        R, t = kabsch(source_frame, target)
        rmsd = float(np.sqrt(np.mean(np.sum((source_frame @ R.T + t - target) ** 2, axis=1))))
        placed = coords @ R.T + t
        keep = protein_clash_mask(placed, lib["elements"], protein, sites,
                                  args.clash_mode, args.clash_cutoff, args.vdw_scale,
                                  exclude_target_residues=not args.include_target_residues_in_clash)
        candidates = np.flatnonzero(keep)
        if len(candidates) == 0:
            raise RuntimeError(f"No library conformers pass protein clash filtering at {site.label}")
        directions = glycan_direction_vectors(placed, target)
        site_data[site] = {"placed": placed, "keep": keep, "candidates": candidates,
                           "rmsd": rmsd, "target_frame": target, "directions": directions}
        print(f"[info] {site.label}: attachment RMSD={rmsd:.4f} A; {len(candidates)}/{len(coords)} pass protein clash filter")

    rng = np.random.default_rng(args.seed)
    accepted: List[Dict[Site, int]] = []
    seen = set()
    attempts = 0

    if args.selection == "diverse" and len(sites) == 1:
        site = sites[0]

        ordered = farthest_point_order(site_data[site]["candidates"],
                                    site_data[site]["directions"],
                                    rng,
                                    n_select=args.n_models,)
        for idx in ordered:
            accepted.append({site: int(idx)})
        
        attempts = len(accepted)
        print(f"[info] {site.label}: selected {len(accepted)} conformers by max-min angular diversity")
    else:
        # For multi-site placement, build clash-free joint combinations first.  In
        # diverse mode, each site's candidate stream follows a farthest-point order;
        # random mode preserves the original behavior.
        if args.selection == "diverse":
            orders = {s: farthest_point_order(site_data[s]["candidates"],
                                              site_data[s]["directions"], rng)
                      for s in sites}
            cursors = {s: 0 for s in sites}
        while len(accepted) < args.n_models and attempts < args.max_attempts:
            attempts += 1
            if args.selection == "random":
                choice = {s: int(rng.choice(site_data[s]["candidates"])) for s in sites}
            else:
                choice = {}
                for offset, s in enumerate(sites):
                    order = orders[s]
                    pos = (cursors[s] + attempts - 1 + offset) % len(order)
                    choice[s] = int(order[pos])
                # advance one site each attempt to avoid repeatedly pairing the same ranks
                cursors[sites[(attempts - 1) % len(sites)]] += 1
            key = tuple(choice[s] for s in sites)
            if key in seen:
                continue
            seen.add(key)
            ok = True
            for i, s1 in enumerate(sites):
                xyz1 = site_data[s1]["placed"][choice[s1]]
                for s2 in sites[i+1:]:
                    xyz2 = site_data[s2]["placed"][choice[s2]]
                    if clash_between(xyz1, lib["elements"], xyz2, lib["elements"],
                                     args.clash_mode, args.clash_cutoff, args.vdw_scale):
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                accepted.append(choice)

    if not accepted:
        raise RuntimeError("No joint multi-glycan combinations passed clash filtering")
    if len(accepted) < args.n_models:
        print(f"[warn] requested {args.n_models}, obtained {len(accepted)} after {attempts} attempts", file=sys.stderr)

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    numbering = build_numbering(protein, sites, lib["residue_numbers"],
                                args.glycan_chain, args.residue_start)
    for rank, choice in enumerate(accepted, 1):
        placements = {s: site_data[s]["placed"][choice[s]] for s in sites}
        write_model(out_dir / f"{args.prefix}_{rank:04d}.pdb", protein, placements, lib, numbering)

    selected_matrix = np.array([[choice[s] for s in sites] for choice in accepted], dtype=int)
    site_labels = np.array([s.label for s in sites])
    placed_payload = dict(
        selected_library_indices=selected_matrix,
        site_labels=site_labels,
        attachment_rmsd=np.array([site_data[s]["rmsd"] for s in sites]),
        coords=np.stack([[site_data[s]["placed"][choice[s]] for s in sites]
                         for choice in accepted]),
        atom_names=lib["atom_names"], residue_numbers=lib["residue_numbers"],
        residue_names=lib["residue_names"], elements=lib["elements"],
        selected_direction_vectors=np.stack([[site_data[s]["directions"][choice[s]] for s in sites]
                                             for choice in accepted]),
    )
    # New libraries generated by glycan_rotamer_generator_updated.py carry one
    # sampled ASN->GlcNAc attachment angle per library conformer.  Preserve the
    # angles selected at each site for troubleshooting and analysis.
    selected_attachment_angles = None
    if "sampled_attachment_angles" in lib:
        source_attachment_angles = np.asarray(lib["sampled_attachment_angles"], float)
        selected_attachment_angles = source_attachment_angles[selected_matrix]
        placed_payload["selected_attachment_angles"] = selected_attachment_angles
        if "attachment_template_angle" in lib:
            placed_payload["library_attachment_template_angle"] = np.asarray(
                lib["attachment_template_angle"]
            )

    # Multi-reference libraries can retain provenance for every conformer.
    # Preserve the source model/path selected at each target site so that a
    # successful placement can be traced back to the Rosetta reference that
    # contributed the glycan geometry.
    selected_source_model_ids = None
    selected_source_pdbs = None
    if "source_model_ids" in lib:
        source_model_ids = np.asarray(lib["source_model_ids"])
        selected_source_model_ids = source_model_ids[selected_matrix]
        placed_payload["selected_source_model_ids"] = selected_source_model_ids
    if "source_pdbs" in lib:
        source_pdbs = np.asarray(lib["source_pdbs"])
        selected_source_pdbs = source_pdbs[selected_matrix]
        placed_payload["selected_source_pdbs"] = selected_source_pdbs
    if "source_reference_index" in lib:
        source_reference_index = np.asarray(lib["source_reference_index"], int)
        placed_payload["selected_source_reference_index"] = source_reference_index[selected_matrix]

    np.savez(out_dir / f"{args.prefix}_placed.npz", **placed_payload)

    manifest = {
        "library": str(Path(args.library).resolve()),
        "protein": str(Path(args.protein).resolve()),
        "sites": [s.label for s in sites],
        "sequons": {s.label: validate_sequon(protein, s)[1] for s in sites},
        "strip_existing_glycans": not args.keep_existing_glycans,
        "library_models": int(len(coords)),
        "protein_clash_free_by_site": {s.label: int(len(site_data[s]["candidates"])) for s in sites},
        "written_models": len(accepted),
        "joint_sampling_attempts": attempts,
        "selected_library_indices": selected_matrix.tolist(),
        "selection": args.selection,
        "orientation_definition": "unit vector from target ASN ND2 to glycan geometric center",
        "selected_direction_vectors": placed_payload["selected_direction_vectors"].tolist(),
        "clash_mode": args.clash_mode,
        "clash_cutoff_A": args.clash_cutoff,
        "vdw_scale": args.vdw_scale,
    }
    if selected_attachment_angles is not None:
        manifest["selected_attachment_angles_deg"] = selected_attachment_angles.tolist()
        if "attachment_template_angle" in lib:
            manifest["library_attachment_template_angle_deg"] = float(
                np.asarray(lib["attachment_template_angle"]).reshape(())
            )
    if selected_source_model_ids is not None:
        manifest["selected_source_model_ids"] = selected_source_model_ids.tolist()
    if selected_source_pdbs is not None:
        manifest["selected_source_pdbs"] = selected_source_pdbs.tolist()
    if "library_kind" in lib:
        manifest["library_kind"] = str(np.asarray(lib["library_kind"]).reshape(()))

    (out_dir / f"{args.prefix}_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"[done] wrote {len(accepted)} multi-glycosylated models to {out_dir}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--library", required=True)
    ap.add_argument("--protein", required=True)
    ap.add_argument("--site", action="append", type=parse_site, required=True,
                    help="site as CHAIN:RESID; repeat for multiple glycans")
    ap.add_argument("--n-models", type=int, default=100)
    ap.add_argument("--max-attempts", type=int, default=100000,
                    help="maximum random joint combinations to test")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--selection", choices=["random", "diverse"], default="random",
                    help="random preserves legacy selection; diverse spreads selected glycans over available COM directions")
    ap.add_argument("--out-dir", default="placed_glycans")
    ap.add_argument("--prefix", default="Glyc_fast")
    ap.add_argument("--clash-mode", choices=["distance", "vdw"], default="vdw")
    ap.add_argument("--clash-cutoff", type=float, default=2.0)
    ap.add_argument("--vdw-scale", type=float, default=0.70)
    ap.add_argument("--include-target-residues-in-clash", action="store_true")
    ap.add_argument("--skip-sequon-validation", action="store_true",
                    help="not recommended; bypass N-X-S/T validation")
    ap.add_argument("--keep-existing-glycans", action="store_true",
                    help="retain carbohydrate atoms in input; default strips all and rebuilds requested sites")
    ap.add_argument("--keep-hydrogens", action="store_true")
    ap.add_argument("--keep-water", action="store_true")
    ap.add_argument("--glycan-chain", default=None)
    ap.add_argument("--residue-start", type=int, default=None)
    return ap


def main(argv=None) -> int:
    try:
        return run(build_parser().parse_args(argv))
    except (ValueError, RuntimeError, FileNotFoundError) as exc:
        print(f"[error] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
