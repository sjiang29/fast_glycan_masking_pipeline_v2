# Fast Glycan Masking Pipeline

A Python package for building reusable N-linked glycan conformer libraries, placing glycans onto target ASN sites, and analyzing glycan orientation diversity.

The pipeline is designed to reuse previously generated glycan conformations so that glycan masking can be evaluated without repeating expensive glycan modeling for every target site.

## Main Workflows

### 1. Build a glycan library

Two approaches are supported:

- **From a template** — generate conformers by sampling a glycosylated template structure.
- **From existing models** — extract glycans from previously generated glycosylated structures, filter source positions using model scores, and align the glycans into a common attachment frame.

### 2. Place glycans

Library conformers are aligned to a target ASN attachment site, filtered for protein clashes, and selected for output. Both normal and orientation-diverse selection are supported.

### 3. Analyze orientation diversity

Orientation analysis can be applied to a complete conformer library or to conformers produced after placement.

## Project Structure

```text
src/fast_glycan_masking/
├── library_builder/
│   ├── from_template/
│   └── from_existing_models/
├── placement/
├── analysis/
└── pipeline.py

tests/       automated tests and small test fixtures
examples/    usage examples
scripts/     convenience scripts
configs/     reusable configuration files
docs/        detailed scientific and technical documentation
work/        local calculations and large regression data (not tracked by Git)
```

## Installation

From the repository root:

```bash
python -m pip install -e .
```

For development:

```bash
python -m pip install -e ".[dev]"
```

## Quick Examples

Build an empirical library from existing models:

```bash
python -m fast_glycan_masking.library_builder.from_existing_models build \
    --input-dir /path/to/models \
    --out-prefix output/Man5_library \
    --samples-per-glycan 0 \
    --seed 42
```

Place library conformers on a target site:

```bash
python -m fast_glycan_masking.placement.placer \
    --library output/Man5_library.npz \
    --protein target.pdb \
    --site A:581 \
    --n-models 10 \
    --out-dir outputs/A581
```

Analyze library orientation diversity:

```bash
python -m fast_glycan_masking.analysis.orientation \
    --library output/Man5_library.npz \
    --out-prefix results/Man5_library
```

## Testing

Run the automated test suite with:

```bash
pytest -v
```

See the module documentation and `docs/` directory for implementation details, library representation, attachment geometry, orientation analysis, and validation.
