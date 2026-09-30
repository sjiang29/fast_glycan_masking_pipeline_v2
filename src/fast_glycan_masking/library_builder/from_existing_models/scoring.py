"""Analyze Rosetta position folders before empirical library construction.

Expected input folders are named ``pos_0``, ``pos_1``, ... . ``pos_0`` provides
the wild-type mean ``total_score`` baseline. Candidate positions are evaluated
from the PDB models that have matching rows in ``Glyc_score.sc``.
"""
from __future__ import annotations
import csv, re
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional

POS_RE = re.compile(r"^pos_(\d+)$")
PDB_RE = re.compile(r"^Glyc.*\d+\.pdb$", re.IGNORECASE)

@dataclass
class PositionScore:
    position: int
    folder: str
    n_pdb: int
    n_score_rows: int
    n_matched: int
    mean_total_score: Optional[float]
    wt_mean_total_score: Optional[float] = None
    delta_reu: Optional[float] = None
    selected: bool = False
    status: str = ""

def discover_positions(input_dir: Path):
    """Return numerically sorted ``(position, folder)`` pairs for ``pos_*`` directories."""
    found = []
    for p in Path(input_dir).iterdir():
        if p.is_dir():
            m = POS_RE.match(p.name)
            if m:
                found.append((int(m.group(1)), p))
    return sorted(found)

def discover_glyc_pdbs(folder: Path) -> Dict[str, Path]:
    """Map Rosetta glycosylated-model stems to PDB paths in one position folder."""
    out = {}
    for p in folder.iterdir():
        if p.is_file() and PDB_RE.match(p.name):
            out[p.stem] = p
    return out

def parse_rosetta_score_file(path: Path):
    """Return description -> total_score from Rosetta SCORE table."""
    if not path.exists():
        return {}
    header = None
    rows = {}
    with path.open() as fh:
        for line in fh:
            if not line.startswith("SCORE:"):
                continue
            parts = line.split()
            vals = parts[1:]
            if "total_score" in vals and "description" in vals:
                header = vals
                continue
            if header is None or len(vals) < len(header):
                continue
            rec = dict(zip(header, vals))
            try:
                rows[Path(rec["description"]).stem] = float(rec["total_score"])
            except (KeyError, ValueError):
                continue
    return rows

def analyze_positions(input_dir: Path, reu_delta: float = 0.0) -> List[PositionScore]:
    """Compare each position mean score with the ``pos_0`` wild-type baseline."""
    input_dir = Path(input_dir)
    positions = discover_positions(input_dir)
    if not any(i == 0 for i,_ in positions):
        raise ValueError(f"{input_dir} has no pos_0 WT folder")

    raw = []
    for pos, folder in positions:
        pdbs = discover_glyc_pdbs(folder)
        scores = parse_rosetta_score_file(folder / "Glyc_score.sc")
        matched = sorted(set(pdbs) & set(scores))
        mean = (sum(scores[k] for k in matched) / len(matched)) if matched else None
        if not pdbs:
            status = "NO_GLYCAN_MODELS"
        elif not scores:
            status = "NO_SCORE_FILE_OR_ROWS"
        elif not matched:
            status = "NO_MATCHED_SCORES"
        else:
            status = "OK"
        raw.append(PositionScore(
            pos, str(folder), len(pdbs), len(scores), len(matched), mean,
            status=status
        ))

    wt = next(x for x in raw if x.position == 0)
    if wt.mean_total_score is None:
        raise ValueError("pos_0 has no Glyc*.pdb models matched to Glyc_score.sc")
    wt_mean = wt.mean_total_score

    for x in raw:
        x.wt_mean_total_score = wt_mean
        if x.mean_total_score is not None:
            x.delta_reu = x.mean_total_score - wt_mean
        if x.position == 0:
            x.selected = False
            x.status = "WT"
        elif x.mean_total_score is not None and x.delta_reu < reu_delta:
            x.selected = True
            x.status = "SELECTED"
        elif x.mean_total_score is not None:
            x.status = "REJECTED_REU"
    return raw

def write_position_report(rows: List[PositionScore], path: Path):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(asdict(rows[0]).keys()) if rows else []
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(asdict(row))
