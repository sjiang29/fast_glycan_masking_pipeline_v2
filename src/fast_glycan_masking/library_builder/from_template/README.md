# Template-Based Library Builder

Builds a glycan conformer library starting from a glycosylated template structure.

## Workflow

```text
glycosylated template
        ↓
extract glycan
        ↓
sample attachment/internal glycan degrees of freedom
        ↓
store conformers in a reusable library
```

This approach is useful when a large collection of pre-generated glycosylated models is not available.

Because all conformers originate from one template, the resulting library may retain orientation or conformational bias from that starting structure. This limitation motivated the complementary `from_existing_models` workflow.

## Usage

Check the current command-line interface with:

```bash
python -m fast_glycan_masking.library_builder.from_template.builder --help
```

For detailed attachment-frame and torsion concepts, see the documents under `docs/`.
