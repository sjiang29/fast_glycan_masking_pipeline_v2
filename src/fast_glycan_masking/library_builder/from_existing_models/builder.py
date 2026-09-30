"""Construct an empirical N-glycan conformer library from Rosetta models.

The builder treats ``pos_0`` as the wild-type score baseline, selects ``pos_*``
folders using the mean Rosetta ``total_score`` criterion, extracts every matched
ASN-linked glycan, maps it into a common ASN attachment frame, and writes a
reusable NPZ library. Optional local torsion sampling can expand the empirical
conformer pool without changing covalent bond lengths or bond angles.

This module intentionally preserves provenance for every conformer so later
analysis can trace a library member back to its source position and PDB model.
"""
from __future__ import annotations
import csv, json
from pathlib import Path
from typing import List
import numpy as np
from .scoring import analyze_positions, discover_glyc_pdbs, parse_rosetta_score_file, write_position_report
from .glycan import find_n_glycans, canonicalize_tree, Atom
from .geometry import transform
from .sampling import expand_conformer

def write_multimodel_pdb(path, atoms, coords):
    with Path(path).open("w") as fh:
        for m,xyz in enumerate(coords,1):
            fh.write(f"MODEL     {m:4d}\n")
            for i,(a,p) in enumerate(zip(atoms,xyz),1):
                fh.write(
                    f"HETATM{i:5d} {a.name:>4s} {a.resname:>3s} G{a.resid:4d}    "
                    f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00          {a.element:>2s}\n"
                )
            fh.write("ENDMDL\n")

def build_library(input_dir, out_prefix, reu_delta=0.0, include_wt_glycans=False,
                  native_positions=None, samples_per_glycan=0, seed=2026,
                  attachment_sigma_deg=20.0, glycosidic_sigma_deg=20.0,
                  max_library_size=None, write_pdb=False, strict=False):
    """Build an empirical glycan library from a ``pos_*`` Rosetta dataset.

    Each selected source model may contribute one or more ASN-linked glycans.
    Conformers are transformed into a shared attachment frame before storage.
    The returned manifest summarizes selection, extraction, failures, and the
    number of empirical/expanded conformers written to ``<out_prefix>.npz``.
    """
    input_dir=Path(input_dir); out=Path(out_prefix)
    out.parent.mkdir(parents=True,exist_ok=True)
    native_positions=set(native_positions or [])

    pos_rows=analyze_positions(input_dir,reu_delta)
    write_position_report(pos_rows,Path(str(out)+"_positions.csv"))
    wt_mean=next(x.wt_mean_total_score for x in pos_rows if x.position==0)

    selected=[x for x in pos_rows if x.selected or (include_wt_glycans and x.position==0)]
    if not selected:
        raise RuntimeError("No position folders passed the REU filter")

    canonical_frame=None; canonical_atoms=None; expected_sig=None; expected_edges=None
    coords=[]; provenance=[]; failures=[]
    rng=np.random.default_rng(seed)

    for prow in selected:
        folder=Path(prow.folder)
        pdbs=discover_glyc_pdbs(folder)
        scores=parse_rosetta_score_file(folder/"Glyc_score.sc")
        matched=sorted(set(pdbs)&set(scores))
        for stem in matched:
            pdb=pdbs[stem]
            try:
                trees=find_n_glycans(pdb)
                if not trees:
                    raise ValueError("no ASN-linked glycans found")
                for tree_index,tree in enumerate(trees):
                    source_atoms,catoms,sig,edges=canonicalize_tree(tree)
                    raw=np.array([a.xyz for a in source_atoms],float)
                    if canonical_frame is None:
                        canonical_frame=tree.asn_frame.copy()
                        canonical_atoms=catoms
                        expected_sig=sig; expected_edges=edges
                    if sig != expected_sig or edges != expected_edges:
                        raise ValueError(
                            f"glycan topology/atom signature differs from canonical library"
                        )
                    norm=transform(raw,tree.asn_frame,canonical_frame)
                    role="engineered" if (prow.position!=0 and tree.attachment_resid==prow.position) else "native"
                    if native_positions and role=="native" and tree.attachment_resid not in native_positions:
                        role="other_n_linked"

                    base_meta=dict(
                        source_folder=folder.name,
                        source_position=prow.position,
                        source_pdb=pdb.name,
                        source_model=stem,
                        source_total_score=scores[stem],
                        position_mean_total_score=prow.mean_total_score,
                        wt_mean_total_score=wt_mean,
                        position_delta_reu=prow.delta_reu,
                        attachment_chain=tree.attachment_chain,
                        attachment_resid=tree.attachment_resid,
                        glycan_role=role,
                        source_tree_index=tree_index,
                    )
                    coords.append(norm)
                    provenance.append({**base_meta,"conformer_origin":"rosetta_empirical",
                                       "parent_empirical_index":len(coords)-1})

                    if samples_per_glycan>0:
                        atom_names=np.array([a.name for a in canonical_atoms])
                        resnums=np.array([a.resid for a in canonical_atoms],int)
                        variants=expand_conformer(
                            norm,atom_names,resnums,expected_edges,canonical_frame,
                            samples_per_glycan,rng,attachment_sigma_deg,glycosidic_sigma_deg
                        )
                        parent_idx=len(coords)-1
                        for v in variants:
                            coords.append(v)
                            provenance.append({**base_meta,"conformer_origin":"torsion_expanded",
                                               "parent_empirical_index":parent_idx})
            except Exception as exc:
                failures.append((str(pdb),str(exc)))
                if strict: raise

    if not coords:
        raise RuntimeError("No glycans were successfully extracted")

    coords=np.asarray(coords,float)
    if max_library_size is not None and len(coords)>max_library_size:
        keep=np.sort(rng.choice(len(coords),size=max_library_size,replace=False))
        coords=coords[keep]
        provenance=[provenance[i] for i in keep]

    atom_names=np.array([a.name for a in canonical_atoms])
    residue_numbers=np.array([a.resid for a in canonical_atoms],int)
    residue_names=np.array([a.resname for a in canonical_atoms])
    elements=np.array([a.element for a in canonical_atoms])
    # Atom-level chain identifiers.
    # The downstream glycan_library_placer expects this field.
    # "G" denotes the standalone glycan chain in the conformer library.
    chains = np.full(len(atom_names), "G", dtype="<U1")
    # String provenance arrays are stored explicitly for easy downstream use.
    np.savez(
        str(out)+".npz",
        coords=coords,
        atom_names=atom_names,
        residue_numbers=residue_numbers,
        residue_names=residue_names,
        chains=chains,
        elements=elements,
        attachment_frame=canonical_frame,
        linkage_edges=np.array(expected_edges,dtype=str),
        source_folder=np.array([x["source_folder"] for x in provenance]),
        source_position=np.array([x["source_position"] for x in provenance],int),
        source_pdb=np.array([x["source_pdb"] for x in provenance]),
        source_model=np.array([x["source_model"] for x in provenance]),
        source_total_score=np.array([x["source_total_score"] for x in provenance],float),
        position_mean_total_score=np.array([x["position_mean_total_score"] for x in provenance],float),
        wt_mean_total_score=np.array([x["wt_mean_total_score"] for x in provenance],float),
        position_delta_reu=np.array([x["position_delta_reu"] for x in provenance],float),
        attachment_chain=np.array([x["attachment_chain"] for x in provenance]),
        attachment_resid=np.array([x["attachment_resid"] for x in provenance],int),
        glycan_role=np.array([x["glycan_role"] for x in provenance]),
        conformer_origin=np.array([x["conformer_origin"] for x in provenance]),
        parent_empirical_index=np.array([x["parent_empirical_index"] for x in provenance],int),
        library_kind=np.array("existing_rosetta_models"),
    )

    prov_path=Path(str(out)+"_provenance.csv")
    with prov_path.open("w",newline="") as fh:
        fields=["library_index"]+list(provenance[0].keys())
        w=csv.DictWriter(fh,fieldnames=fields); w.writeheader()
        for i,row in enumerate(provenance):
            w.writerow({"library_index":i,**row})

    if write_pdb:
        write_multimodel_pdb(str(out)+".pdb",canonical_atoms,coords)

    manifest=dict(
        input_dir=str(input_dir.resolve()), reu_delta=reu_delta,
        include_wt_glycans=include_wt_glycans,
        native_positions=sorted(native_positions),
        selected_positions=[x.position for x in selected],
        empirical_and_expanded_conformers=len(coords),
        samples_per_glycan=samples_per_glycan,
        attachment_sigma_deg=attachment_sigma_deg,
        glycosidic_sigma_deg=glycosidic_sigma_deg,
        failed_pdbs=len(failures), seed=seed,
    )
    Path(str(out)+"_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    if failures:
        with Path(str(out)+"_skipped.tsv").open("w") as fh:
            fh.write("pdb\treason\n")
            for p,e in failures: fh.write(f"{p}\t{e}\n")
    return manifest
