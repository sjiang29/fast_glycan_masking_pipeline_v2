"""Smoke tests for the package structure.

These tests catch broken module paths and relative imports after code
reorganization. They intentionally do not test scientific behavior.
"""


def test_package_imports():
    import fast_glycan_masking
    import fast_glycan_masking.library_builder
    import fast_glycan_masking.placement
    import fast_glycan_masking.analysis


def test_template_builder_imports():
    import fast_glycan_masking.library_builder.from_template.builder
    import fast_glycan_masking.library_builder.from_template.rotamer


def test_existing_models_builder_imports():
    import fast_glycan_masking.library_builder.from_existing_models.builder
    import fast_glycan_masking.library_builder.from_existing_models.scoring
    import fast_glycan_masking.library_builder.from_existing_models.glycan
    import fast_glycan_masking.library_builder.from_existing_models.sampling
    import fast_glycan_masking.library_builder.from_existing_models.geometry


def test_placement_imports():
    import fast_glycan_masking.placement.placer


def test_analysis_imports():
    import fast_glycan_masking.analysis.orientation


def test_pipeline_imports():
    import fast_glycan_masking.pipeline
