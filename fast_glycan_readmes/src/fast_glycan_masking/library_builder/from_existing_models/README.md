# Existing-Model Library Builder

Builds an empirical glycan conformer library from previously generated glycosylated protein models.

This workflow is intended to reuse expensive modeling calculations rather than regenerating glycan conformations from scratch.

## Expected Input

A source dataset is organized into position directories:

```text
input/
├── pos_0/
├── pos_1/
├── pos_2/
└── ...
```

`pos_0` is the reference/wild-type calculation. A `pos_i` directory represents models generated for candidate position `i`.

The builder uses `Glyc_score.sc` and the `total_score` column for score-based position filtering. Glycosylated structures are identified from `Glyc*.pdb` files. Positions with no matching glycosylated models are skipped.

## Workflow

```text
read WT scores from pos_0
        ↓
calculate mean total_score
        ↓
compare each candidate position with WT
        ↓
retain accepted positions
        ↓
extract glycans from accepted models
        ↓
align glycans into a common attachment frame
        ↓
optionally expand conformers
        ↓
write NPZ library + manifest
```

## Example

```bash
python -m fast_glycan_masking.library_builder.from_existing_models build \
    --input-dir work/regression_test \
    --out-prefix work/regression_test/output/Man5_test \
    --samples-per-glycan 0 \
    --seed 42
```

Use `--samples-per-glycan 0` for a purely empirical library without additional torsional expansion.

For all options:

```bash
python -m fast_glycan_masking.library_builder.from_existing_models build --help
```

## Outputs

The primary outputs are:

```text
<out-prefix>.npz
<out-prefix>_manifest.json
```

The NPZ stores glycan coordinates, atom metadata, attachment-frame information, and source provenance used by downstream placement and analysis.

For detailed NPZ fields and attachment-frame conventions, see `docs/LIBRARY_FORMAT.md` and `docs/ATTACHMENT_GEOMETRY.md`.
