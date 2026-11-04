# Glycan Library Merger

## Overview

`merge_library.py` merges two glycan conformer libraries stored in NumPy `.npz` format into a single library.

The script is designed for the **Fast Glycan Masking Pipeline**, particularly for combining Man5 conformers generated from different sources, such as Rosetta empirical models and torsion-sampled conformers.

The merger supports:

- Heavy-atom normalization
- Glycan topology validation
- Automatic atom reordering
- ASN attachment-frame alignment using Kabsch superposition
- RMSD-based conformer deduplication
- Preservation of per-conformer provenance
- JSON summary generation

The resulting library can be used for downstream glycan placement and orientation analysis, provided those tools support heavy-atom-only conformers.

## Requirements

- Python 3.9+
- NumPy
- SciPy

Install dependencies:

```bash
pip install numpy scipy
```

## Usage

### Basic merge

```bash
python merge_library.py \
    --library-a Man5_flu.npz \
    --library-b Man5_library_new.npz \
    --out Man5_merged.npz
```

By default, conformers are deduplicated using a **0.5 Å RMSD threshold**.

### Specify deduplication threshold

```bash
python merge_library.py \
    --library-a Man5_flu.npz \
    --library-b Man5_library_new.npz \
    --out Man5_merged.npz \
    --dedup-rmsd 0.5
```

### Merge without deduplication

```bash
python merge_library.py \
    --library-a Man5_flu.npz \
    --library-b Man5_library_new.npz \
    --out Man5_merged.npz \
    --no-dedup
```

## Command-Line Arguments

| Argument | Description | Default |
|---|---|---|
| `--library-a` | First input NPZ library | Required |
| `--library-b` | Second input NPZ library | Required |
| `--out` | Output merged NPZ file | Required |
| `--dedup-rmsd` | RMSD threshold in Å for duplicate removal | `0.5` |
| `--no-dedup` | Disable deduplication | `False` |

`--dedup-rmsd` and `--no-dedup` are mutually exclusive.

## Input Library Requirements

Both input libraries must contain:

- `coords`: conformer coordinates, shape `(N, A, 3)`
- `atom_names`: atom names
- `residue_numbers`: glycan residue indices
- `residue_names`: residue names
- `elements`: chemical elements
- `attachment_frame`: ASN attachment coordinates, shape `(3, 3)`

The attachment frame must contain the coordinates of ASN atoms in the following order:

```text
CB
CG
ND2
```

Each library must also provide glycan linkage information through either `linkage_edges` or the supported parent/child linkage arrays.

### Different library formats

The script supports libraries with different atom counts due to the presence or absence of hydrogens.

For example:

| Property | Man5_flu.npz | Man5_library_new.npz |
|---|---|---|
| Conformers | 64,125 | 10,000 |
| Total atoms | 160 | 83 |
| Heavy atoms | 83 | 83 |
| Hydrogens | 77 | 0 |

The merger removes hydrogens from the coordinate representation, normalizes residue numbering, and matches corresponding heavy atoms.

**Important:** The output library contains heavy atoms only. Hydrogen coordinates are not retained.

## Merge Workflow

### 1. Validate input libraries

The script checks that:

- Required NPZ fields exist
- Coordinate and attachment-frame dimensions are valid
- Coordinates are finite
- Heavy-atom identities match
- Glycan linkage topology matches
- ASN attachment frames are geometrically compatible

If the heavy-atom sets or linkage topologies differ, the merge is rejected.

### 2. Normalize atom ordering

Heavy atoms are identified using:

```text
(residue index, residue name, atom name, element)
```

Residue indices are normalized according to their order of first appearance.

Library B is reordered to match Library A's normalized heavy-atom ordering.

This allows libraries with different residue numbering conventions to be combined.

### 3. Align ASN attachment frames

Each library may have been generated in a different global coordinate frame.

The script uses the Kabsch algorithm to calculate a rigid-body transformation between the ASN attachment frames.

Library B's coordinates are transformed into Library A's frame:

```text
Library B coordinates
        |
        v
Kabsch alignment using ASN CB, CG, ND2
        |
        v
Coordinates in Library A's reference frame
```

The alignment preserves the glycan's orientation relative to its ASN attachment site.

The merge is rejected if the attachment-frame RMSD exceeds **0.3 Å**.

### 4. Combine conformers

After normalization and alignment, the two coordinate arrays are concatenated.

For the example libraries:

```text
Man5_flu.npz          64,125 conformers
Man5_library_new.npz  10,000 conformers
                     -----------------
Total                 74,125 conformers
```

This is the total before deduplication.

### 5. Remove duplicate conformers

By default, the script uses a 0.5 Å heavy-atom RMSD cutoff.

For two conformers with N corresponding heavy atoms:

\[
\mathrm{RMSD} =
\sqrt{\frac{1}{N}\sum_{i=1}^{N}
\|\mathbf{x}_i-\mathbf{y}_i\|^2}
\]

**No additional glycan superposition is performed during deduplication.**

This is important because the goal is to preserve conformational and orientation diversity relative to the ASN attachment site.

Two conformers with similar internal glycan structures but different attachment orientations are therefore not automatically considered duplicates.

Deduplication uses a greedy, first-kept strategy with a SciPy KD-tree to accelerate comparisons against previously retained conformers.

Library A is processed before Library B, so Library A has priority when duplicates occur.

The method is order-dependent and is not a global clustering algorithm.

## Provenance Preservation

The script preserves available per-conformer metadata from both input libraries.

Two additional fields identify the origin of every retained conformer:

| Field | Description |
|---|---|
| `source_library` | Absolute path to the original NPZ library |
| `original_library_index` | Zero-based conformer index in the original library |

Other available provenance fields, such as `source_pdb`, `source_position`, `source_total_score`, `conformer_origin`, and `parent_empirical_index`, are also retained.

If a metadata field exists in only one library, missing values are filled according to data type:

- Floating-point values: `NaN`
- Integer values: `-1`
- Strings: empty string
- Boolean values: `False`

If `parent_empirical_index` exists, the script attempts to remap parent indices to the merged library. Parents removed by deduplication are represented by `-1`.

The original parent indices are also stored in `original_parent_empirical_index`.

Static metadata that cannot be merged as per-conformer arrays is recorded in the JSON summary.

## Output Files

For:

```bash
--out Man5_merged.npz
```

the script generates:

```text
Man5_merged.npz
Man5_merged_merge_summary.json
```

### Merged NPZ

The output contains:

- Merged heavy-atom coordinates
- Normalized atom and residue metadata
- ASN attachment frame
- Glycan linkage information
- Per-conformer provenance
- Source library identifiers
- Original conformer indices

The `library_kind` field is set to:

```text
merged_man5_heavy_atom
```

### JSON Summary

The summary includes:

| Field | Description |
|---|---|
| `input_a` | Number of conformers in Library A |
| `input_b` | Number of conformers in Library B |
| `total` | Combined conformers before deduplication |
| `retained` | Conformers retained after deduplication |
| `duplicates_removed` | Number of removed conformers |
| `retained_a` | Retained conformers from Library A |
| `retained_b` | Retained conformers from Library B |
| `heavy_atoms` | Number of heavy atoms per conformer |
| `dedup_rmsd_angstrom` | Applied RMSD cutoff |
| `asn_frame_alignment_rmsd_angstrom` | ASN-frame alignment RMSD |
| `original_static_metadata` | Non-conformer metadata from both inputs |

The exact number of retained conformers depends on the similarity between the input libraries.

## Downstream Orientation Analysis

The merged library can be analyzed using the Fast Glycan Masking Pipeline's orientation analysis tool.

For example:

```bash
python -m fast_glycan_masking.analysis.orientation \
    --library Man5_merged.npz \
    --coordinate-frame asn-local \
    --out-prefix merged_orientation
```

The ASN-local coordinate system allows glycan orientations to be compared independently of the original global coordinate frames.

Useful comparisons include:

- Root C1 azimuth distributions
- Root C1 polar-angle distributions
- Whole-glycan center orientations
- Orientation-space occupancy
- Angular diversity before and after merging

These analyses can help determine whether combining the libraries improves orientation diversity.

## Limitations and Considerations

1. **Matching topology is required.** The script does not convert between different glycan types or linkage topologies.
2. **Heavy-atom-only output.** Hydrogen coordinates are discarded.
3. **Frame consistency is required.** Both attachment frames must represent ASN CB, CG, and ND2 in that order.
4. **Order-dependent deduplication.** Library A conformers are preferred over duplicate conformers from Library B.
5. **No additional RMSD alignment.** Deduplication measures differences in the common ASN frame.
6. **Large libraries may require substantial memory and computation.** The KD-tree is rebuilt periodically as retained conformers accumulate.
7. **Downstream compatibility should be checked.** In particular, confirm that glycan placement tools accept the normalized residue numbering and heavy-atom-only representation.

## Recommended Workflow

```text
Rosetta empirical library ──────┐
                                │
                                v
                         Validate topology
                                ^
                                │
Torsion-sampled library ────────┘
                                |
                                v
                      Normalize heavy atoms
                                |
                                v
                       Align ASN frames
                                |
                                v
                       Merge conformers
                                |
                                v
                     Deduplicate (0.5 Å)
                                |
                                v
                    Preserve provenance
                                |
                                v
                      Merged Man5 library
                                |
                   +------------+------------+
                   |                         |
                   v                         v
          Orientation analysis       Glycan placement
```
