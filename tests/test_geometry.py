"""Unit tests for geometry utilities used by empirical library construction.

The geometry module aligns coordinates between two reference frames using
Kabsch rigid-body alignment. These tests verify that the transformation
preserves coordinates correctly for identity, translation, and rotation.
"""

import numpy as np

from fast_glycan_masking.library_builder.from_existing_models.geometry import (
    kabsch,
    transform,
)


def rmsd(a, b):
    """Calculate RMSD between corresponding Cartesian coordinates."""
    return np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1)))


def test_transform_identity():
    """Transforming between identical frames should leave coordinates unchanged."""

    frame = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )

    coords = np.array(
        [
            [0.2, 0.3, 0.0],
            [0.7, 0.1, 0.5],
            [1.2, 0.4, -0.3],
        ]
    )

    aligned = transform(coords, frame, frame)

    assert np.allclose(aligned, coords, atol=1e-8)


def test_transform_translation():
    """Transform should correctly map coordinates between translated frames."""

    source_frame = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )

    translation = np.array([5.0, -3.0, 2.0])

    target_frame = source_frame + translation

    coords = np.array(
        [
            [0.2, 0.3, 0.0],
            [0.7, 0.1, 0.5],
            [1.2, 0.4, -0.3],
        ]
    )

    expected = coords + translation

    aligned = transform(coords, source_frame, target_frame)

    assert np.allclose(aligned, expected, atol=1e-8)
    assert rmsd(aligned, expected) < 1e-8


def test_transform_rotation_and_translation():
    """Transform should correctly map coordinates through rotation and translation."""

    source_frame = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )

    # 90-degree rotation around z.
    rotation = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0,  0.0, 0.0],
            [0.0,  0.0, 1.0],
        ]
    )

    translation = np.array([5.0, -3.0, 2.0])

    target_frame = source_frame @ rotation.T + translation

    coords = np.array(
        [
            [0.2, 0.3, 0.0],
            [0.7, 0.1, 0.5],
            [1.2, 0.4, -0.3],
        ]
    )

    expected = coords @ rotation.T + translation

    aligned = transform(coords, source_frame, target_frame)

    assert np.allclose(aligned, expected, atol=1e-8)
    assert rmsd(aligned, expected) < 1e-8


def test_kabsch_returns_rotation_and_translation():
    """Kabsch should return a 3x3 rotation matrix and 3-vector translation."""

    source = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )

    target = source + np.array([5.0, -3.0, 2.0])

    R, t = kabsch(source, target)

    assert R.shape == (3, 3)
    assert t.shape == (3,)
