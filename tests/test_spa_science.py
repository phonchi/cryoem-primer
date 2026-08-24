"""Focused scientific-contract tests for the SPA teaching chapters."""

from __future__ import annotations

import ast
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SYNTHETIC_SOURCE = ROOT / "book" / "07_synthetic_data.py"


def _load_functions(*names):
    tree = ast.parse(SYNTHETIC_SOURCE.read_text(encoding="utf-8"))
    selected = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in names
    ]
    namespace = {"np": np}
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(SYNTHETIC_SOURCE), "exec"), namespace)
    return [namespace[name] for name in names]


def test_isotropic_sampler_is_reproducible_and_uniform_in_cos_tilt():
    sample_isotropic_view_angles, = _load_functions("sample_isotropic_view_angles")
    rot_a, tilt_a = sample_isotropic_view_angles(np.random.default_rng(7), 100_000)
    rot_b, tilt_b = sample_isotropic_view_angles(np.random.default_rng(7), 100_000)

    np.testing.assert_array_equal(rot_a, rot_b)
    np.testing.assert_array_equal(tilt_a, tilt_b)
    assert np.all((0 <= rot_a) & (rot_a < 2 * np.pi))
    assert abs(np.mean(np.cos(tilt_a))) < 0.01
    assert abs(np.mean(np.cos(tilt_a) ** 2) - 1 / 3) < 0.01


def test_view_expansion_preserves_image_count_and_balances_repeats():
    expand_view_angles, = _load_functions("expand_view_angles")
    rot = np.linspace(0, 1, 7)
    tilt = np.linspace(0.2, 2.8, 7)
    angles, view_ids = expand_view_angles(rot, tilt, 100)

    assert angles.shape == (100, 3)
    assert view_ids.shape == (100,)
    counts = np.bincount(view_ids)
    assert counts.max() - counts.min() <= 1
    assert np.all((0 <= angles[:, 2]) & (angles[:, 2] < 2 * np.pi))


def test_ctf_dc_and_200_kv_wavelength_match_chapter_convention():
    electron_wavelength, ctf_1d = _load_functions("electron_wavelength", "ctf_1d")

    assert abs(electron_wavelength(200) - 0.02508) < 2e-4
    assert ctf_1d(np.array([0.0]), 15_000, w=0.1)[0] == -0.1


def test_global_awgn_realizes_requested_variance_snr():
    add_global_awgn, = _load_functions("add_global_awgn")
    clean = np.random.default_rng(11).normal(size=(512, 64)).astype(np.float32)
    noisy, sigma, standardized_noise = add_global_awgn(
        clean, 0.1, np.random.default_rng(12)
    )
    noise = noisy - clean
    realized = np.var(clean, dtype=np.float64) / np.var(noise, dtype=np.float64)

    assert sigma > 0
    assert standardized_noise.shape == clean.shape
    assert abs(realized - 0.1) < 0.003


def test_synthetic_chapter_declares_reproducibility_contract():
    source = SYNTHETIC_SOURCE.read_text(encoding="utf-8")
    required_terms = {
        "CRYOEM_NUM_IMAGES",
        "CRYOEM_OUTPUT_DIR",
        "ground_truth.jsonl",
        "run_manifest.json",
        "schema_version",
        "euler_convention",
        "rotation_matrix",
        "rotation_semantics",
        "volume_sha256",
        "software_versions",
        "get_metadata(as_dict=True)",
    }
    assert not {term for term in required_terms if term not in source}
    assert "relion_source._metadata" not in source
    assert "預設 ±L/16" not in source


def test_new_spa_pages_state_key_claim_boundaries():
    formation = (ROOT / "book" / "05_image_formation.md").read_text(encoding="utf-8")
    validation = (ROOT / "book" / "06_reconstruction_validation.md").read_text(encoding="utf-8")

    for term in ("弱相位物體", "reference wave", "DQE", "Dose fractionation", "SSNR"):
        assert term in formation
    assert r"\mathrm{SNR}_{\mathrm{out}}^2" not in formation
    assert r"\frac{\mathrm{SNR}_{\mathrm{out}}(s)}{\mathrm{SNR}_{\mathrm{in}}(s)}" in formation
    for term in ("Fourier slice theorem", "soft assignment", "half-set", "FSC", "Mask"):
        assert term in validation
    assert "不直接衡量真實性" in validation
