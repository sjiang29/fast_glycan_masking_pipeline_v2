# Glycan Placement

The placement module transfers conformers from a glycan library to one or more target ASN sites.

## Workflow

```text
glycan library + target protein
        ↓
construct target ASN attachment frame
        ↓
align library conformers to target site
        ↓
protein clash filtering
        ↓
select output conformers
        ↓
write glycosylated models
```

## Example

```bash
python -m fast_glycan_masking.placement.placer \
    --library output/Man5_library.npz \
    --protein target.pdb \
    --site A:581 \
    --n-models 10 \
    --out-dir outputs/A581 \
    --prefix Glyc_A581
```

## Selection

The placer supports normal selection and orientation-diverse selection.

Diverse selection uses glycan orientation information to choose conformers that are separated in orientation space. It cannot recover orientations that were already removed by protein clash filtering.

The placer also writes a `*_placed.npz` file that can be used by the orientation-analysis module.

For details on attachment alignment and orientation-diverse selection, see `docs/ATTACHMENT_GEOMETRY.md` and `docs/ORIENTATION_DIVERSITY.md`.
