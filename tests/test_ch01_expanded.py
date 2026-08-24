"""Numerical checks for the expanded image-basics chapter."""

from contextlib import contextmanager
from pathlib import Path
import os
import runpy

import matplotlib
import numpy as np
import pytest


matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPOSITORY = Path(__file__).resolve().parents[1]
BOOK = REPOSITORY / "book"


@contextmanager
def working_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


@pytest.fixture(scope="module")
def chapter():
    with working_directory(BOOK):
        namespace = runpy.run_path(str(BOOK / "01_image_basics.py"))
    yield namespace
    plt.close("all")


def test_roi_mask_and_color_round_trip(chapter):
    assert chapter["roi"].shape == (100, 100)
    assert chapter["circular_mask"].dtype == bool
    assert np.all(chapter["masked_coins"][~chapter["circular_mask"]] == 0)
    assert np.allclose(chapter["rgb_from_hsv"], chapter["rgb"], atol=1e-12)
    assert np.array_equal(chapter["rgb_restored"], chapter["rgb_uint8"])


def test_lossless_and_lossy_codec_behavior(chapter):
    round_trip = chapter["codec_round_trip"]
    image = np.arange(64 * 64, dtype=np.uint16).reshape(64, 64)

    png, _ = round_trip(image, ".png")
    tiff, _ = round_trip(image, ".tiff")
    assert png.dtype == np.uint16 and np.array_equal(png, image)
    assert tiff.dtype == np.uint16 and np.array_equal(tiff, image)

    color = chapter["sample_bgr"]
    jpeg, _ = round_trip(
        color, ".jpg", [chapter["cv2"].IMWRITE_JPEG_QUALITY, 80]
    )
    assert jpeg.shape == color.shape
    assert not np.array_equal(jpeg, color)


def test_histogram_cdf_is_normalized_and_monotone(chapter):
    histogram, _, cdf = chapter["histogram_and_cdf"](
        np.array([[0.0, 0.25], [0.5, 1.0]])
    )
    assert histogram.sum() == pytest.approx(1.0)
    assert np.all(np.diff(cdf) >= 0)
    assert cdf[-1] == pytest.approx(1.0)


def test_euclidean_transform_preserves_distance(chapter):
    transform = chapter["forward"]
    points = np.array([[4.0, 7.0], [19.0, -3.0]])
    transformed = transform(points)
    assert np.linalg.norm(transformed[1] - transformed[0]) == pytest.approx(
        np.linalg.norm(points[1] - points[0])
    )
    assert chapter["wide_canvas"].shape == (320, 320)
    assert not np.allclose(chapter["reflected_boundary"], chapter["bilinear"])


def test_antialiasing_suppresses_unrepresentable_checkerboard(chapter):
    without_antialiasing = chapter["downsampled_no_aa"]
    with_antialiasing = chapter["downsampled_aa"]
    assert with_antialiasing.std() < without_antialiasing.std() / 5
