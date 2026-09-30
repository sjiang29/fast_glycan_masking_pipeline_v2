"""Command-line interface for empirical library analysis and construction."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
from .scoring import analyze_positions, write_position_report
from .builder import build_library

def parse_positions(text):
    if not text: return []
    return [int(x.strip()) for x in text.split(",") if x.strip()]

def parser():
    ap=argparse.ArgumentParser(
        prog="fast_glycan_masking.library_builder.from_existing_models",
        description="Build N-glycan libraries from pre-generated Rosetta glycoprotein models."
    )
    sub=ap.add_subparsers(dest="cmd",required=True)

    a=sub.add_parser("analyze",help="Compare pos_* mean total_score against pos_0")
    a.add_argument("--input-dir",required=True)
    a.add_argument("--reu-delta",type=float,default=0.0,
                   help="select when mean(pos)-mean(pos_0) < this value; default 0")
    a.add_argument("--out",default="position_scores.csv")

    b=sub.add_parser("build",help="Analyze positions, extract glycans, and build library")
    b.add_argument("--input-dir",required=True)
    b.add_argument("--out-prefix",required=True)
    b.add_argument("--reu-delta",type=float,default=0.0)
    b.add_argument("--include-wt-glycans",action="store_true")
    b.add_argument("--native-positions",default="",
                   help="optional comma-separated known native ASN positions, e.g. 141,200")
    b.add_argument("--samples-per-glycan",type=int,default=0,
                   help="additional torsion-expanded conformers per empirical glycan")
    b.add_argument("--attachment-sigma-deg",type=float,default=20.0)
    b.add_argument("--glycosidic-sigma-deg",type=float,default=20.0)
    b.add_argument("--max-library-size",type=int,default=None)
    b.add_argument("--seed",type=int,default=2026)
    b.add_argument("--write-pdb",action="store_true",
                   help="also write a potentially very large multi-model PDB")
    b.add_argument("--strict",action="store_true")
    return ap

def main(argv=None):
    args=parser().parse_args(argv)
    try:
        if args.cmd=="analyze":
            rows=analyze_positions(Path(args.input_dir),args.reu_delta)
            write_position_report(rows,Path(args.out))
            print(f"[done] wrote {args.out}")
            for r in rows:
                print(f"pos_{r.position}: n={r.n_matched} mean={r.mean_total_score} "
                      f"delta={r.delta_reu} status={r.status}")
        else:
            m=build_library(
                args.input_dir,args.out_prefix,args.reu_delta,args.include_wt_glycans,
                parse_positions(args.native_positions),args.samples_per_glycan,args.seed,
                args.attachment_sigma_deg,args.glycosidic_sigma_deg,
                args.max_library_size,args.write_pdb,args.strict
            )
            print(f"[done] library: {args.out_prefix}.npz")
            print(f"[done] conformers: {m['empirical_and_expanded_conformers']}")
            print(f"[done] selected positions: {m['selected_positions']}")
        return 0
    except Exception as exc:
        print(f"[error] {exc}",file=sys.stderr)
        return 2

if __name__=="__main__":
    raise SystemExit(main())
