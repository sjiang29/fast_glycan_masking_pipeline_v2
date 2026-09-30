"""Glycan conformer-library construction workflows.

Two builders are provided:

* :mod:`from_template` generates conformers by torsion sampling from one
  glycosylated reference structure.
* :mod:`from_existing_models` extracts empirical conformers from collections
  of pre-generated Rosetta glycoprotein models.

Both workflows write reusable NPZ libraries for the placement module.
"""
