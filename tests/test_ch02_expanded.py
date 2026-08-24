"""Numerical checks for the expanded filtering and segmentation chapter."""

from pathlib import Path
import runpy

import matplotlib
import numpy as np
import pytest


matplotlib.use("Agg")
import matplotlib.pyplot as plt


CHAPTER = Path(__file__).resolve().parents[1] / "book" / "02_filter_segment.py"


@pytest.fixture(scope="module")
def chapter():
    namespace = runpy.run_path(str(CHAPTER))
    plt.close("all")
    return namespace


def test_noise_model_changes_filter_ranking(chapter):
    gaussian_errors = chapter["gaussian_errors"]
    impulse_errors = chapter["impulse_errors"]

    assert gaussian_errors["Gaussian"] < gaussian_errors["median"]
    assert impulse_errors["median"] < impulse_errors["Gaussian"]
    assert gaussian_errors["Gaussian"] < gaussian_errors["noisy"]
    assert impulse_errors["median"] < impulse_errors["noisy"]


def test_local_threshold_handles_the_example_illumination_ramp(chapter):
    scores = chapter["threshold_scores"]

    assert scores["local"] > scores["fixed"]
    assert scores["local"] > scores["Otsu"]
    assert chapter["binary_iou"](
        np.zeros((3, 3), dtype=bool), np.zeros((3, 3), dtype=bool)
    ) == 1.0


def test_morphology_restores_three_known_components(chapter):
    labels = chapter["component_labels"]
    properties = chapter["component_properties"]

    assert labels.max() == 3
    assert len(properties) == 3
    assert chapter["clean_mask"][62, 55]


def test_gradient_pyramid_and_hog_geometry(chapter):
    assert chapter["gradient_magnitude"].shape == chapter["edge_input"].shape
    assert chapter["canny_smoothed"].sum() < chapter["canny_without_smoothing"].sum()
    assert len(chapter["gaussian_pyramid"]) == 4
    assert chapter["features"].size == chapter["calculated_hog_length"]
    assert chapter["features"].size == 14_076
