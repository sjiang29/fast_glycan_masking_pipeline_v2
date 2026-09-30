# Documentation

This directory contains the detailed scientific and technical documentation for the Fast Glycan Masking Pipeline.

The documentation is separated from the root README and component READMEs so that each level has a clear purpose:

- **Root `README.md`** — project overview, installation, and quick-start commands.
- **Component `README.md` files** — how to use a specific part of the package.
- **`docs/`** — deeper explanations of algorithms, coordinate conventions, design decisions, validation, and scientific assumptions.

## Planned Documentation

- `ARCHITECTURE.md` — overall data flow and relationships between package components.
- `ATTACHMENT_GEOMETRY.md` — ASN attachment frame, rigid-body alignment, Kabsch transformation, and attachment RMSD.
- `LIBRARY_FORMAT.md` — NPZ schema, structural fields, provenance fields, and compatibility requirements.
- `ORIENTATION_DIVERSITY.md` — glycan direction vectors, angular distance, orientation bins, entropy, and farthest-point sampling.
- `MULTI_SOURCE_LIBRARIES.md` — combining libraries generated from multiple source proteins while preserving provenance and controlling source bias.
- `TESTING.md` — unit, integration, and regression testing strategy.

Detailed documentation should describe *why* an algorithm or representation is used. Usage instructions that are specific to one component should remain in that component's README.
