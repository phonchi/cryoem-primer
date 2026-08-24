"""Executable checks for the four Jupytext foundation chapters.

The chapter sources contain their own numerical assertions.  Running them here
keeps prose examples and tested computations on the same code path.
"""

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


def run_chapter(filename):
    with working_directory(BOOK):
        namespace = runpy.run_path(str(BOOK / filename))
    plt.close("all")
    return namespace


@pytest.mark.parametrize(
    "filename",
    [
        "01_image_basics.py",
        "02_filter_segment.py",
        "03_fourier.py",
        "04_wavelet.py",
    ],
)
def test_foundation_chapter_executes(filename):
    run_chapter(filename)


def test_rectangular_dft_and_parseval():
    namespace = run_chapter("03_fourier.py")
    image = np.arange(24, dtype=float).reshape(4, 6)
    direct = namespace["direct_dft2"](image)
    transformed = np.fft.fft2(image)

    assert np.allclose(direct, transformed)
    assert np.allclose(np.fft.ifft2(transformed).real, image)
    assert transformed[0, 0] == pytest.approx(image.sum())
    assert np.sum(image**2) == pytest.approx(
        np.sum(np.abs(transformed) ** 2) / image.size
    )


def test_watershed_returns_individual_objects_and_zero_background():
    namespace = run_chapter("02_filter_segment.py")
    yy, xx = np.mgrid[:100, :140]
    mask = ((yy - 50) ** 2 + (xx - 45) ** 2 < 28**2) | (
        (yy - 50) ** 2 + (xx - 88) ** 2 < 28**2
    )
    _, _, labels, marker_count = namespace["split_touching_objects"](
        mask, min_distance=24
    )

    assert marker_count == 2
    assert labels.max() == 2
    assert np.all(labels[~mask] == 0)


def test_dwt_round_trip_for_odd_length_signal():
    namespace = run_chapter("04_wavelet.py")
    signal = np.random.default_rng(10).normal(size=259)

    for mode in ("zero", "symmetric", "periodization"):
        reconstructed, _ = namespace["dwt_round_trip"](signal, mode=mode)
        assert np.allclose(reconstructed, signal, atol=1e-10)


def test_phase_flip_preserves_observed_magnitude_away_from_zeros():
    namespace = run_chapter("03_fourier.py")
    observed = np.arange(16).reshape(4, 4) + 1j
    ctf = np.linspace(-1, 1, 16).reshape(4, 4)
    corrected = namespace["phase_flip_spectrum"](observed, ctf)
    nonzero = ctf != 0

    assert np.allclose(np.abs(corrected[nonzero]), np.abs(observed[nonzero]))
