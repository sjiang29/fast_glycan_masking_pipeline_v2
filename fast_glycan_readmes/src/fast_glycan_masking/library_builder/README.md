# Library Builders

This package contains the two supported approaches for constructing reusable glycan conformer libraries.

## Builders

- `from_template/` — generates conformers from a glycosylated template by sampling glycan degrees of freedom.
- `from_existing_models/` — extracts conformers from previously generated glycosylated models and aligns them into a common attachment frame.

Both builders produce libraries that can be consumed by the placement module.

For detailed scientific descriptions of attachment geometry and the library representation, see `docs/ATTACHMENT_GEOMETRY.md` and `docs/LIBRARY_FORMAT.md`.
