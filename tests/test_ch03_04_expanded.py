"""Numerical invariants of the restored Fourier and wavelet examples."""

from __future__ import annotations

import numpy as np
import pywt
from scipy import ndimage, signal


def test_original_phase_exchange_preserves_donor_magnitude_and_phase(restored_chapter):
    chapter = restored_chapter("03_fourier.py")
    lucario = np.fft.fft2(chapter["lucario"])
    pikachu = np.fft.fft2(chapter["pikachu"])
    exchanges = [
        ("lucario_magnitude_pikachu_phase", lucario, pikachu),
        ("pikachu_magnitude_lucario_phase", pikachu, lucario),
    ]
    for name, magnitude_donor, phase_donor in exchanges:
        observed = np.fft.fft2(chapter[name])
        np.testing.assert_allclose(np.abs(observed), np.abs(magnitude_donor), atol=1e-9)
        # Phase is informative only where both spectra have nonzero amplitude.
        valid = (np.abs(magnitude_donor) > 1e-7) & (np.abs(phase_donor) > 1e-7)
        assert np.count_nonzero(valid) > 0
        relative_phase = observed[valid] * phase_donor[valid].conj()
        np.testing.assert_allclose(
            relative_phase / np.abs(relative_phase), 1, atol=1e-7
        )


def test_lucario_direct_and_fft_convolutions_use_same_normalized_kernel(restored_chapter):
    chapter = restored_chapter("03_fourier.py")
    # Independent spatial reference for the original 3 x 3, sigma=3 timing example.
    axis = np.arange(-1, 2, dtype=float)
    weights = np.exp(-(axis**2) / (2 * 3**2))
    kernel = np.outer(weights, weights)
    kernel /= kernel.sum()
    expected = signal.convolve2d(
        chapter["lucario"], kernel, mode="same", boundary="fill", fillvalue=0
    )
    np.testing.assert_allclose(chapter["im_blurred1"], expected, atol=1e-12)
    np.testing.assert_allclose(chapter["im_blurred2"], expected, atol=1e-12)


def test_moonlanding_peak_suppression_preserves_protected_frequencies(restored_chapter):
    chapter = restored_chapter("03_fourier.py")
    original_fft = chapter["F"]
    filtered_fft = chapter["F_dim"]
    removed = chapter["remove_peaks"]
    protected = chapter["protected_center"]
    assert removed.any() and (~removed).any()
    assert not removed[protected].any()
    np.testing.assert_array_equal(filtered_fft[~removed], original_fft[~removed])
    np.testing.assert_array_equal(filtered_fft[removed], 0)
    assert filtered_fft[0, 0] == original_fft[0, 0]

    rows = (-np.arange(removed.shape[0])) % removed.shape[0]
    cols = (-np.arange(removed.shape[1])) % removed.shape[1]
    np.testing.assert_array_equal(removed, removed[np.ix_(rows, cols)])
    reconstructed = np.fft.ifft2(filtered_fft)
    assert np.max(np.abs(reconstructed.imag)) < 1e-9
    np.testing.assert_allclose(chapter["image_filtered"], reconstructed.real, atol=1e-10)


def test_pikachu_laplacian_api_returns_image_minus_same_scale_smoothing(restored_chapter):
    chapter = restored_chapter("03_fourier.py")
    image = chapter["pikachu_rgb"]
    # pyramid_laplacian defaults to sigma = 2 * downscale / 6, downscale=2.
    # This is not the residual against an expanded image from the next level.
    smoothed = ndimage.gaussian_filter(image, sigma=(2 / 3, 2 / 3, 0), mode="reflect")
    first_detail = chapter["laplacian_pyramid"][0]
    np.testing.assert_allclose(first_detail, image - smoothed, atol=1e-12)
    assert first_detail.min() < 0 < first_detail.max()
    gaussian_levels = chapter["pyramid"]
    np.testing.assert_array_equal(gaussian_levels[0], image)
    assert all(
        coarse.shape[:2] == tuple((size + 1) // 2 for size in fine.shape[:2])
        for fine, coarse in zip(gaussian_levels, gaussian_levels[1:])
    )


def test_ecg_soft_threshold_changes_details_but_preserves_coarse_approximation(restored_chapter):
    chapter = restored_chapter("04_wavelet.py")
    noisy = chapter["x_noisy"]
    before = pywt.wavedec(noisy, "db4", mode="periodization")
    after = pywt.wavedec(chapter["rec"], "db4", mode="periodization")
    np.testing.assert_allclose(after[0], before[0], atol=1e-10)
    threshold = 0.1 * np.max(noisy)
    assert threshold > 0
    for original, shrunk in zip(before[1:], after[1:]):
        assert np.all(np.abs(shrunk) <= np.abs(original) + 1e-10)
        np.testing.assert_allclose(shrunk[np.abs(original) <= threshold], 0, atol=1e-10)
        large = np.abs(original) > threshold
        np.testing.assert_allclose(
            np.abs(original[large]) - np.abs(shrunk[large]), threshold, atol=1e-10
        )
    np.testing.assert_allclose(chapter["lowpassfilter"](noisy, 0), noisy, atol=1e-10)


def test_rgb_crop_psnr_uses_the_same_unit_data_range_for_all_methods(restored_chapter):
    chapter = restored_chapter("04_wavelet.py")
    reference = chapter["original"]
    variants = [
        ("noisy", "psnr_noisy"),
        ("im_bayes", "psnr_bayes"),
        ("im_visushrink", "psnr_visushrink"),
        ("im_visushrink2", "psnr_visushrink2"),
        ("im_visushrink4", "psnr_visushrink4"),
    ]
    for image_name, score_name in variants:
        image = chapter[image_name]
        assert image.shape == reference.shape
        mse = np.mean((image - reference) ** 2)
        assert mse > 0
        np.testing.assert_allclose(chapter[score_name], 10 * np.log10(1 / mse), atol=1e-10)
