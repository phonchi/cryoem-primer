"""Numerical checks for the expanded Fourier and multiscale lessons."""

from __future__ import annotations

import os
import runpy
from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cryoem-tests")

BOOK = Path(__file__).resolve().parents[1] / "book"


@pytest.fixture(scope="module")
def fourier():
    namespace = runpy.run_path(str(BOOK / "03_fourier.py"))
    namespace["plt"].close("all")
    return namespace


@pytest.fixture(scope="module")
def wavelet():
    namespace = runpy.run_path(str(BOOK / "04_wavelet.py"))
    namespace["plt"].close("all")
    return namespace


def test_translation_changes_phase_but_not_magnitude(fourier):
    original_fft = fourier["camera_fft"]
    shifted_fft = np.fft.fft2(fourier["shifted_camera"])
    assert np.allclose(np.abs(shifted_fft), np.abs(original_fft))
    assert np.allclose(shifted_fft, fourier["predicted_shift_fft"], atol=1e-9)


def test_spatial_and_fft_linear_convolution_agree(fourier):
    assert np.allclose(fourier["spatial_full"], fourier["fft_full"], atol=1e-12)


def test_frequency_filters_and_regularization_are_numerically_meaningful(fourier):
    assert fourier["ideal_impulse"].min() < -1e-4
    assert fourier["gaussian_impulse"].min() > -1e-10
    assert fourier["mse_values"].min() < fourier["observed_mse"]
    assert np.allclose(
        fourier["psnr_values"],
        10 * np.log10(1.0 / fourier["mse_values"]),
    )


def test_pyramid_and_multilevel_dwt_reconstruct(wavelet):
    assert np.allclose(
        wavelet["pyramid_reconstruction"], wavelet["pyramid_image"], atol=1e-12
    )
    assert np.allclose(
        wavelet["reconstruction_2d"], wavelet["image"], atol=1e-12
    )


def test_wavelet_denoising_improves_psnr_and_ssim(wavelet):
    quality = wavelet["quality"]
    assert quality["BayesShrink"][0] > quality["noisy"][0]
    assert quality["BayesShrink"][1] > quality["noisy"][1]
    assert abs(wavelet["sigma_mad"] - wavelet["sigma_true"]) < 0.03


def test_swt_is_shift_equivariant_for_periodized_example(wavelet):
    for (approximation, detail), (shifted_approximation, shifted_detail) in zip(
        wavelet["swt_original"], wavelet["swt_shifted"]
    ):
        assert np.allclose(shifted_approximation, np.roll(approximation, 1), atol=1e-10)
        assert np.allclose(shifted_detail, np.roll(detail, 1), atol=1e-10)

