#!/usr/bin/env python3
"""Create new N-X-S/T sites with Rosetta FastDesign, then place Man5 at all sites.

This wrapper keeps the scientifically important part of the old pipeline:
FastDesign is still responsible for introducing a valid sequon and building a
reasonable ASN side chain. The slow GlycanTreeModeler ensemble step is replaced
by the reusable-library placement module.

New masking sites are supplied with --new-site. Native sites that should also be
resampled are supplied with --native-site. All sites are placed simultaneously.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List

from .placement.placer import (Site, parse_site, read_structure, residue_name,
                               validate_sequon)

AA3_TO_1 = {
    "ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C", "GLN":"Q",
    "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I", "LEU":"L", "LYS":"K",
    "MET":"M", "PHE":"F", "PRO":"P", "SER":"S", "THR":"T", "TRP":"W",
    "TYR":"Y", "VAL":"V"
}


def required_mutations(records, sites: List[Site]):
    mutations = {}
    for site in sites:
        r0 = residue_name(records, site.chain, site.resid)
        r1 = residue_name(records, site.chain, site.resid + 1)
        r2 = residue_name(records, site.chain, site.resid + 2)
        if None in (r0, r1, r2):
            raise ValueError(f"Cannot create sequon at {site.label}: i/i+1/i+2 not all present")
        if r0 != "ASN":
            mutations[(site.chain, site.resid)] = "N"
        if r1 == "PRO":
            mutations[(site.chain, site.resid + 1)] = "A"
        if r2 not in {"SER", "THR"}:
            mutations[(site.chain, site.resid + 2)] = "T"
    return mutations


def write_resfile(path: Path, mutations) -> None:
    lines = ["NATRO", "EX 1 EX 2", "start"]
    for (chain, resid), aa in sorted(mutations.items(), key=lambda x: (x[0][0], x[0][1])):
        lines.append(f"{resid} {chain} PIKAA {aa}")
    path.write_text("\n".join(lines) + "\n")


def make_fastdesign_xml(template: Path, output: Path, resfile: Path) -> None:
    tree = ET.parse(template)
    root = tree.getroot()
    movers = list(root.iter("ReadResfile"))
    if not movers:
        raise ValueError(f"No <ReadResfile> element found in {template}")
    for mover in movers:
        mover.set("filename", str(resfile))
    tree.write(output)


def infer_designed_pdb(input_pdb: Path) -> Path:
    return input_pdb.with_name(f"Des_{input_pdb.stem}_0001.pdb")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--protein", required=True)
    ap.add_argument("--library", required=True)
    ap.add_argument("--new-site", action="append", type=parse_site, default=[],
                    help="engineered site CHAIN:RESID; repeat as needed")
    ap.add_argument("--native-site", action="append", type=parse_site, default=[],
                    help="native site CHAIN:RESID to resample; repeat as needed")
    ap.add_argument("--n-models", type=int, default=100)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--max-attempts", type=int, default=100000)
    ap.add_argument("--out-dir", default="fast_glycan_models")
    ap.add_argument("--prefix", default="Glyc_Des_fast")
    ap.add_argument("--clash-mode", choices=["vdw", "distance"], default="vdw")
    ap.add_argument("--vdw-scale", type=float, default=0.70)
    ap.add_argument("--clash-cutoff", type=float, default=2.0)
    ap.add_argument("--fastdesign-template", default="template_FastDesign.xml")
    ap.add_argument("--fastdesign-script", default="fast_design.sh")
    ap.add_argument("--work-dir", default="fastdesign_work")
    ap.add_argument("--skip-fastdesign", action="store_true",
                    help="require all new sites to already be valid sequons")
    args = ap.parse_args(argv)

    sites = list(dict.fromkeys(args.new_site + args.native_site))
    if not sites:
        ap.error("Supply at least one --new-site or --native-site")

    protein_path = Path(args.protein).resolve()
    records = read_structure(str(protein_path), keep_existing_glycans=False)
    mutations = required_mutations(records, args.new_site)
    designed_pdb = protein_path

    if mutations:
        if args.skip_fastdesign:
            msgs = [validate_sequon(records, s)[1] for s in args.new_site]
            raise SystemExit("New sites are not valid sequons and --skip-fastdesign was used:\n" + "\n".join(msgs))
        work = Path(args.work_dir).resolve(); work.mkdir(parents=True, exist_ok=True)
        local_pdb = work / protein_path.name
        shutil.copy2(protein_path, local_pdb)
        resfile = work / "fast_man5_sites.resfile"
        xml = work / "FastDesign.xml"
        write_resfile(resfile, mutations)
        make_fastdesign_xml(Path(args.fastdesign_template).resolve(), xml, resfile)
        script = Path(args.fastdesign_script).resolve()
        if not script.exists():
            raise SystemExit(f"FastDesign script not found: {script}")
        print("[info] FastDesign mutations:")
        for (chain, resid), aa in sorted(mutations.items()):
            print(f"  {chain}{resid} -> {aa}")
        subprocess.run(["bash", str(script), local_pdb.name], cwd=work, check=True)
        designed_pdb = work / f"Des_{local_pdb.stem}_0001.pdb"
        if not designed_pdb.exists():
            raise SystemExit(f"FastDesign completed but expected output was not found: {designed_pdb}")
    else:
        print("[info] All requested new sites already satisfy N-X-S/T; FastDesign not needed")

    # Validate all sites in the structure that will be placed.
    designed_records = read_structure(str(designed_pdb), keep_existing_glycans=False)
    invalid = []
    for site in sites:
        ok, msg = validate_sequon(designed_records, site)
        print(f"[sequon] {msg}")
        if not ok:
            invalid.append(msg)
    if invalid:
        raise SystemExit("Sequon validation failed after FastDesign:\n" + "\n".join(invalid))

    # Invoke the placer as a package module so this wrapper remains valid after
    # source-tree reorganization and when installed with pip.
    cmd = [sys.executable, "-m", "fast_glycan_masking.placement.placer",
           "--library", str(Path(args.library).resolve()),
           "--protein", str(designed_pdb), "--n-models", str(args.n_models),
           "--max-attempts", str(args.max_attempts), "--seed", str(args.seed),
           "--out-dir", args.out_dir, "--prefix", args.prefix,
           "--clash-mode", args.clash_mode, "--vdw-scale", str(args.vdw_scale),
           "--clash-cutoff", str(args.clash_cutoff)]
    for site in sites:
        cmd += ["--site", f"{site.chain}:{site.resid}"]
    print("[info] Running multi-site fast Man5 placement")
    subprocess.run(cmd, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
