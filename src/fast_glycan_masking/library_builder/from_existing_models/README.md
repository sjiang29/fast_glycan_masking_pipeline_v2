# Existing-Model Library Builder

Builds an empirical glycan conformer library from previously generated glycosylated protein models.

This workflow reuses expensive modeling calculations rather than regenerating glycan conformations from scratch. It supports optional torsional expansion and RMSD-based conformer deduplication.

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
align glycans into a common ASN attachment frame
        ↓
optionally expand conformers
        ↓
optionally remove duplicate conformers (RMSD)
        ↓
write NPZ library + manifest
```

## Example

Build a purely empirical library without torsional expansion or deduplication:

```bash
python -m fast_glycan_masking.library_builder.from_existing_models build \
    --input-dir work/regression_test \
    --out-prefix work/regression_test/output/Man5_test \
    --samples-per-glycan 0 \
    --seed 42
```

Build a library with RMSD-based deduplication:

```bash
python -m fast_glycan_masking.library_builder.from_existing_models build \
    --input-dir work/regression_test \
    --out-prefix work/regression_test/output/Man5_dedup \
    --samples-per-glycan 0 \
    --dedup-rmsd 0.5 \
    --seed 42
```

Build a library with torsional expansion and deduplication:

```bash
python -m fast_glycan_masking.library_builder.from_existing_models build \
    --input-dir work/regression_test \
    --out-prefix work/regression_test/output/Man5_expanded \
    --samples-per-glycan 10 \
    --dedup-rmsd 0.5 \
    --seed 42
```

Use `--samples-per-glycan 0` for a purely empirical library without additional torsional expansion.

## Conformer Deduplication

The optional `--dedup-rmsd` argument removes identical or near-identical glycan conformers based on their heavy-atom RMSD.

All extracted glycans are first aligned into a common attachment frame using the ASN side-chain atoms **CB, CG, and ND2**. RMSD is then calculated directly between corresponding glycan heavy atoms without additional structural superposition.

For conformers A and B:

\[
\mathrm{RMSD}(A,B)=\sqrt{\frac{1}{N}\sum_{i=1}^{N}\|\mathbf r_i^A-\mathbf r_i^B\|^2}
\]

where \(N\) is the number of corresponding glycan heavy atoms.

The deduplication procedure:

1. Processes conformers in their existing order.
2. Retains the first conformer.
3. Discards subsequent conformers whose RMSD is below the specified threshold relative to an already retained conformer.
4. Preserves the coordinates and provenance of retained conformers.

Because no additional glycan alignment is performed, conformers with different orientations relative to ASN are treated as structurally distinct.

Suggested thresholds for testing:

| Threshold | Interpretation |
|---|---|
| `0.1` Å | Removes nearly identical conformers |
| `0.25` Å | Conservative near-duplicate removal |
| `0.5` Å | Moderate deduplication |
| `1.0` Å | More aggressive filtering |

Deduplication is disabled by default. The `--dedup-rmsd` option can be used with or without torsional expansion.

The procedure is greedy and order-dependent. The retained conformer set may depend on the input order.

## Command-Line Options

| Argument | Description | Default |
|---|---|---|
| `--samples-per-glycan` | Number of additional torsion-expanded conformers per empirical glycan | `0` |
| `--attachment-sigma-deg` | Standard deviation for attachment torsion sampling | `20.0` |
| `--glycosidic-sigma-deg` | Standard deviation for glycosidic torsion sampling | `20.0` |
| `--dedup-rmsd` | Heavy-atom RMSD threshold for deduplication (Å) | Disabled |
| `--max-library-size` | Maximum number of conformers retained | No limit |
| `--seed` | Random seed for reproducibility | `2026` |

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

When deduplication is enabled, the NPZ contains only retained conformers. The manifest records the deduplication threshold and conformer counts before and after filtering.

For detailed NPZ fields and attachment-frame conventions, see `docs/LIBRARY_FORMAT.md` and `docs/ATTACHMENT_GEOMETRY.md`.
