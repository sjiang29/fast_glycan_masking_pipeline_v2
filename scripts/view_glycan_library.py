#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Export glycan conformers from an NPZ library to PDB files."
    )

    parser.add_argument(
        "--library",
        type=Path,
        required=True,
        help="Path to glycan conformer library (.npz)"
    )

    parser.add_argument(
        "--n-models",
        type=int,
        default=100,
        help="Number of conformers to export (default: 100)"
    )

    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: <library_name>_visualization)"
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible sampling (default: 42)"
    )

    return parser.parse_args()


def export_pdb(
    pdb_path,
    coordinates,
    atom_names,
    residue_names,
    residue_numbers,
    elements,
):
    with open(pdb_path, "w") as f:
        for j, xyz in enumerate(coordinates):
            x, y, z = xyz

            atom_name = str(atom_names[j])
            residue_name = str(residue_names[j])
            residue_number = int(residue_numbers[j])
            element = str(elements[j])

            f.write(
                f"HETATM{j + 1:5d} "
                f"{atom_name:4s} "
                f"{residue_name:3s} "
                f"A{residue_number:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{1.00:6.2f}{0.00:6.2f}"
                f"          {element:>2s}\n"
            )

        f.write("END\n")


def main():

    """
    python view_glycan_library.py \
    --library library/Man5_library.npz \
    --n-models 200 \
    --out-dir outputs/man5_visualization
    """
    args = parse_args()
    if not args.library.is_file():
        raise FileNotFoundError(
            f"Library not found: {args.library}"
        )

    if args.n_models <= 0:
        raise ValueError("--n-models must be positive")

    out_dir = (
        args.out_dir
        if args.out_dir is not None
        else Path(f"{args.library.stem}_visualization")
    )

    out_dir.mkdir(parents=True, exist_ok=True)

    required = [
        "coords",
        "atom_names",
        "residue_names",
        "residue_numbers",
        "elements",
    ]

    with np.load(args.library, allow_pickle=False) as library:
        missing = [key for key in required if key not in library.files]

        if missing:
            raise ValueError(
                f"Missing required library fields: {missing}"
            )

        coords = library["coords"]
        atom_names = library["atom_names"]
        residue_names = library["residue_names"]
        residue_numbers = library["residue_numbers"]
        elements = library["elements"]

    n_total = len(coords)
    n_export = min(args.n_models, n_total)

    if n_total == 0:
        raise ValueError("Library contains no conformers")

    if coords.ndim != 3 or coords.shape[2] != 3:
        raise ValueError(
            f"Expected coords shape (N, atoms, 3), got {coords.shape}"
        )

    n_atoms = coords.shape[1]

    for name, values in [
        ("atom_names", atom_names),
        ("residue_names", residue_names),
        ("residue_numbers", residue_numbers),
        ("elements", elements),
    ]:
        if len(values) != n_atoms:
            raise ValueError(
                f"{name} has {len(values)} entries; expected {n_atoms}"
            )

    print(f"Library: {args.library}")
    print(f"Total conformers: {n_total}")
    print(f"Atoms per conformer: {n_atoms}")
    print(f"Exporting: {n_export}")
    print(f"Output directory: {out_dir}")

    rng = np.random.default_rng(args.seed)

    indices = rng.choice(
        n_total,
        size=n_export,
        replace=False
    )

    for i, idx in enumerate(indices, start=1):
        pdb_path = out_dir / f"{args.library.stem}_{i:04d}.pdb"

        export_pdb(
            pdb_path,
            coords[idx],
            atom_names,
            residue_names,
            residue_numbers,
            elements,
        )

    # Save source indices for reproducibility
    np.savetxt(
        out_dir / "selected_indices.txt",
        indices,
        fmt="%d",
    )

    print(f"\nDone! Exported {n_export} PDB files.")


if __name__ == "__main__":
    main()
