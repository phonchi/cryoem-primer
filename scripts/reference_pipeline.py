#!/usr/bin/env python3
"""Build and validate the public reference metadata for the primer.

The full extracted text remains in ``References/document_cache``.  This module
only publishes metadata, routing decisions, claim provenance, and an asset
inventory.  Paths are resolved relative to the repository unless overridden;
there are deliberately no machine-specific absolute defaults.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = 1
PIPELINE_VERSION = "1.1"
ROUTING_STATUSES = {"full_digest", "route_note", "out_of_scope"}


def _chapter_key(year: str, number: int) -> str:
    return f"mie{year}-ch{number:02d}"


# Project-specific routing is intentionally kept here, separate from the
# generic document-folder-ingest skill.
CHAPTERS: dict[tuple[str, str], dict[str, str]] = {}


def _chapter(
    year: str,
    number: int,
    filename_token: str,
    title: str,
    status: str,
    digest: str,
    cite_key: str = "",
) -> None:
    CHAPTERS[(year, filename_token)] = {
        "source_id": _chapter_key(year, number),
        "title": title,
        "coverage_status": status,
        "digest_path": digest,
        "citation_key": cite_key,
    }


_chapter("2010", 1, "Chapter-One", "Fundamentals of Three-Dimensional Reconstruction from Projections", "full_digest", "notes/ref_digest/mie2010-ch1-3d-reconstruction.md", "penczek2010")
_chapter("2010", 2, "Chapter-Two", "Image Restoration in Cryo-Electron Microscopy", "full_digest", "notes/ref_digest/mie2010-ch2-image-restoration.md", "penczek2010restoration")
_chapter("2010", 3, "Chapter-Three", "Resolution Measures in Molecular Electron Microscopy", "full_digest", "notes/ref_digest/mie2010-core-addenda.md#resolution-measures", "penczek2010resolution")
_chapter("2010", 4, "Chapter-Four", "3D Reconstruction from 2D Crystal Image and Diffraction Data", "out_of_scope", "notes/ref_digest/scope-boundaries.md#electron-crystallography")
_chapter("2010", 5, "Chapter-Five", "Fourier-Bessel Reconstruction of Helical Assemblies", "out_of_scope", "notes/ref_digest/scope-boundaries.md#helical-reconstruction")
_chapter("2010", 6, "Chapter-Six", "Reconstruction of Helical Filaments and Tubes", "out_of_scope", "notes/ref_digest/scope-boundaries.md#helical-reconstruction")
_chapter("2010", 7, "Chapter-Seven", "Three-Dimensional Asymmetric Reconstruction of Tailed Bacteriophage", "route_note", "notes/ref_digest/route-notes.md#asymmetric-virus-reconstruction")
_chapter("2010", 8, "Chapter-Eight", "Single Particle Analysis at High Resolution", "full_digest", "notes/ref_digest/mie2010-core-addenda.md#high-resolution-spa", "cong2010")
_chapter("2010", 9, "Chapter-Nine", "The Orthogonal Tilt Reconstruction Method", "route_note", "notes/ref_digest/route-notes.md#tilt-based-initial-models")
_chapter("2010", 10, "Chapter-Ten", "An Introduction to Maximum-Likelihood Methods in Cryo-EM", "full_digest", "notes/ref_digest/mie2010-ch10-maximum-likelihood.md", "scheres2010")
_chapter("2010", 11, "Chapter-Eleven", "Classification of Structural Heterogeneity by Maximum-Likelihood Methods", "full_digest", "notes/ref_digest/mie2010-core-addenda.md#maximum-likelihood-classification", "scheres2010heterogeneity")
_chapter("2010", 12, "Chapter-Twelve", "Methods for Three-Dimensional Reconstruction of Heterogeneous Assemblies", "full_digest", "notes/ref_digest/mie2010-core-addenda.md#heterogeneous-reconstruction", "leschziner2010heterogeneity")
_chapter("2010", 13, "Chapter-Thirteen", "Alignment of Cryo-Electron Tomograms without Fiducial Markers", "out_of_scope", "notes/ref_digest/scope-boundaries.md#cryo-electron-tomography")
_chapter("2010", 14, "Chapter-Fourteen", "Correcting for the Ewald Sphere in High-Resolution Single-Particle Reconstructions", "full_digest", "notes/ref_digest/mie2010-core-addenda.md#ewald-sphere", "leong2010")
_chapter("2010", 15, "Chapter-Fifteen", "Software Tools for Molecular Microscopy", "route_note", "notes/ref_digest/route-notes.md#software-landscape")

_chapter("2016", 1, "Chapter-One", "Direct Electron Detectors", "full_digest", "notes/ref_digest/mie2016-core-addenda.md#direct-electron-detectors", "mcmullan2016")
_chapter("2016", 2, "Chapter-Two", "Specimen Behavior in the Electron Beam", "full_digest", "notes/ref_digest/mie2016-core-addenda.md#specimen-behavior", "glaeser2016")
_chapter("2016", 3, "Chapter-Three", "Specimen Preparation for High-Resolution Cryo-EM", "full_digest", "notes/ref_digest/mie2016-core-addenda.md#specimen-preparation", "passmore2016")
_chapter("2016", 4, "Chapter-Four", "Strategies for Automated Cryo-EM Data Collection", "full_digest", "notes/ref_digest/mie2016-core-addenda.md#automated-data-collection", "cheng2016")
_chapter("2016", 5, "Chapter-Five", "Processing of Cryo-EM Movie Data", "full_digest", "notes/ref_digest/mie2016-ch5-movie-processing.md", "rubinstein2016")
_chapter("2016", 6, "Chapter-Six", "Processing of Structurally Heterogeneous Cryo-EM Data in RELION", "full_digest", "notes/ref_digest/mie2016-core-addenda.md#relion-and-heterogeneity", "scheres2016")
_chapter("2016", 7, "Chapter-Seven", "Single-Particle Refinement and Variability Analysis in EMAN2.1", "full_digest", "notes/ref_digest/mie2016-core-addenda.md#refinement-and-variability", "ludtke2016")
_chapter("2016", 8, "Chapter-Eight", "FREALIGN: An Exploratory Tool for Single-Particle Cryo-EM", "full_digest", "notes/ref_digest/mie2016-core-addenda.md#frealign", "grigorieff2016")
_chapter("2016", 9, "Chapter-Nine", "Testing the Validity of Single-Particle Maps at Low and High Resolution", "full_digest", "notes/ref_digest/mie2016-core-addenda.md#map-validation", "rosenthal2016")
_chapter("2016", 10, "Chapter-Ten", "Tools for Model Building and Optimization into Electron Cryo-Microscopy Density Maps", "out_of_scope", "notes/ref_digest/scope-boundaries.md#atomic-model-building")
_chapter("2016", 11, "Chapter-Eleven", "Refinement of Atomic Structures Against Cryo-EM Maps", "out_of_scope", "notes/ref_digest/scope-boundaries.md#atomic-model-building")
_chapter("2016", 12, "Chapter-Twelve", "Cryo-EM Structure Determination Using Segger", "out_of_scope", "notes/ref_digest/scope-boundaries.md#helical-reconstruction")
_chapter("2016", 13, "Chapter-Thirteen", "Cryo-Electron Tomography and Subtomogram Averaging", "out_of_scope", "notes/ref_digest/scope-boundaries.md#cryo-electron-tomography")
_chapter("2016", 14, "Chapter-Fourteen", "High-Resolution Macromolecular Structure Determination by MicroED", "out_of_scope", "notes/ref_digest/scope-boundaries.md#electron-crystallography")
_chapter("2016", 15, "Chapter-Fifteen", "Databases and Archiving for Cryo-EM", "route_note", "notes/ref_digest/route-notes.md#databases-and-archiving")


STANDALONE: dict[str, dict[str, str]] = {
    "Computer Vision - A Modern Approach.pdf": {
        "source_id": "forsyth2012", "title": "Computer Vision: A Modern Approach", "coverage_status": "route_note", "digest_path": "notes/ref_digest/cv-modern-approach-index.md", "citation_key": "forsyth2012"
    },
    "Computer Vision algorithms and applications.pdf": {
        "source_id": "szeliski2022", "title": "Computer Vision: Algorithms and Applications", "coverage_status": "route_note", "digest_path": "notes/ref_digest/szeliski-cv-index.md", "citation_key": "szeliski2022"
    },
    "Learning OpenCV 4 Computer Vision with Python 3.pdf": {
        "source_id": "opencv2020", "title": "Learning OpenCV 4 Computer Vision with Python 3", "coverage_status": "route_note", "digest_path": "notes/ref_digest/learning-opencv4-index.md", "citation_key": "howse2020"
    },
    "[2016] Principles of cryo-EM single-particle image.pdf": {
        "source_id": "sigworth2016", "title": "Principles of Cryo-EM Single-Particle Image Processing", "coverage_status": "full_digest", "digest_path": "notes/ref_digest/2015-principles-spa.md", "citation_key": "sigworth2016"
    },
    "[2020] IEEE_Review.pdf": {
        "source_id": "bendory2020", "title": "Single-Particle Cryo-EM: Mathematical Theory, Computational Challenges, and Opportunities", "coverage_status": "full_digest", "digest_path": "notes/ref_digest/2020-ieee-review.md", "citation_key": "bendory2020"
    },
    "[2020] Math_review.pdf": {
        "source_id": "singer2020", "title": "Computational Methods for Single-Particle Electron Cryomicroscopy", "coverage_status": "full_digest", "digest_path": "notes/ref_digest/2020-math-review.md", "citation_key": "singer2020"
    },
    "[2020] Math_review_sup.pdf": {
        "source_id": "singer2020-supp", "title": "Supplementary Material: Computational Methods for Single-Particle Electron Cryomicroscopy", "coverage_status": "full_digest", "digest_path": "notes/ref_digest/2020-math-review.md", "citation_key": "singer2020"
    },
}

# Explicitly reviewed additions; alternate PDFs retain distinct hashes.
STANDALONE['[1-1] The EM algorithm.pdf'] = {'source_id': 'ma2019em', 'title': 'The EM algorithm', 'coverage_status': 'full_digest', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#ma2019em', 'citation_key': 'ma2019em'}
STANDALONE['[1998] ML_RFA.pdf'] = {'source_id': 'sigworth1998', 'title': 'A Maximum-Likelihood Approach to Single-Particle Image Refinement', 'coverage_status': 'full_digest', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#sigworth1998', 'citation_key': 'sigworth1998'}
STANDALONE['[2005] ML2D.pdf'] = {'source_id': 'scheres2005', 'title': 'Maximum-likelihood Multi-reference Refinement for Electron Microscopy Images', 'coverage_status': 'full_digest', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#scheres2005', 'citation_key': 'scheres2005'}
STANDALONE['[2010] CL2D.pdf'] = {'source_id': 'sorzano2010', 'title': 'A clustering approach to multireference alignment of single-particle projections in electron microscopy', 'coverage_status': 'full_digest', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#sorzano2010', 'citation_key': 'sorzano2010'}
STANDALONE['[2012] ISAC.pdf'] = {'source_id': 'yang2012', 'title': 'Iterative Stable Alignment and Clustering of 2D Transmission Electron Microscope Images', 'coverage_status': 'full_digest', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#yang2012', 'citation_key': 'yang2012'}
STANDALONE['[2012] ISAC-sup.pdf'] = {'source_id': 'yang2012-supp', 'title': 'Supplemental Experimental Procedures: ISAC', 'coverage_status': 'full_digest', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#yang2012-supp', 'citation_key': 'yang2012'}
STANDALONE['[2012] Relion.pdf'] = {'source_id': 'scheres2012bayes', 'title': 'A Bayesian View on Cryo-EM Structure Determination', 'coverage_status': 'full_digest', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#scheres2012bayes', 'citation_key': 'scheres2012bayes'}
STANDALONE['[2014] gSUP.pdf'] = {'source_id': 'chen2014', 'title': 'gamma-SUP: A Clustering Algorithm for Cryo-Electron Microscopy Images of Asymmetric Particles', 'coverage_status': 'full_digest', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#chen2014', 'citation_key': 'chen2014'}
STANDALONE['[2015] A Primer to cryoEM.pdf'] = {'source_id': 'cheng2015primer', 'title': 'A Primer to Single-Particle Cryo-Electron Microscopy', 'coverage_status': 'full_digest', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#cheng2015primer', 'citation_key': 'cheng2015primer'}
STANDALONE['[2017] CryoSparc.pdf'] = {'source_id': 'punjani2017', 'title': 'cryoSPARC: algorithms for rapid unsupervised cryo-EM structure determination', 'coverage_status': 'full_digest', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#punjani2017', 'citation_key': 'punjani2017'}
STANDALONE['[2017] CryoSparc_Sup.pdf'] = {'source_id': 'punjani2017-supp', 'title': 'Supplementary Information: cryoSPARC', 'coverage_status': 'full_digest', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#punjani2017-supp', 'citation_key': 'punjani2017'}
STANDALONE['[2021] 3DVA.pdf'] = {'source_id': 'punjani2021', 'title': '3D variability analysis: Resolving continuous flexibility and discrete heterogeneity from single particle cryo-EM', 'coverage_status': 'full_digest', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#punjani2021', 'citation_key': 'punjani2021'}
STANDALONE['[2022] On bias_var_overfitting_spa.pdf'] = {'source_id': 'sorzano2022', 'title': 'On bias, variance, overfitting, gold standard and consensus in single-particle analysis by cryo-electron microscopy', 'coverage_status': 'full_digest', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#sorzano2022', 'citation_key': 'sorzano2022'}
STANDALONE['[5] 2SDR.pdf'] = {'source_id': 'chung2020', 'title': 'Two-Stage Dimension Reduction for Noisy High-Dimensional Images and Application to Cryogenic Electron Microscopy', 'coverage_status': 'full_digest', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#chung2020', 'citation_key': 'chung2020'}
STANDALONE['[2] ML-EM in cryo-EM.pdf'] = {'source_id': 'mie2010-ch10-alt', 'title': 'An Introduction to Maximum-Likelihood Methods in Cryo-EM (alternate PDF)', 'coverage_status': 'route_note', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#mie2010-ch10-alt', 'citation_key': 'scheres2010'}
STANDALONE['[2] ML-EM review (Section4).pdf'] = {'source_id': 'singer2020-alt', 'title': 'Computational Methods for Single-Particle Electron Cryomicroscopy (alternate PDF)', 'coverage_status': 'route_note', 'digest_path': 'notes/ref_digest/spa-core-20260911.md#singer2020-alt', 'citation_key': 'singer2020'}


CLAIMS: list[dict[str, Any]] = [{'claim_id': 'image-forward-model',
  'chapter': '07_synthetic_data',
  'anchor': 'forward-model',
  'claim': 'A particle image is modeled as a rotated 3D density projected to 2D, shifted, filtered by the '
           'point-spread function or CTF, and corrupted by noise.',
  'source_id': 'bendory2020',
  'locator': 'pp. 62-63, Eqs. (1)-(2)',
  'scope': 'Single-particle cryo-EM under the projection approximation.'},
 {'claim_id': 'fourier-slice-theorem',
  'chapter': '06_reconstruction_validation',
  'anchor': 'fourier-coverage',
  'claim': 'The 2D Fourier transform of a projection is a central plane through the 3D Fourier transform of the '
           'object.',
  'source_id': 'sigworth2016',
  'locator': 'p. 59',
  'scope': 'Projection approximation; Ewald-sphere curvature is neglected.'},
 {'claim_id': 'ctf-single-image-zero',
  'chapter': '05_image_formation',
  'anchor': 'ctf-correction-example',
  'claim': 'Information at a CTF zero cannot be recovered from that image alone; multiple defocus groups provide '
           'complementary frequency coverage.',
  'source_id': 'sigworth2016',
  'locator': 'p. 60',
  'scope': 'Idealized single-image transfer; recovery relies on other observations and regularization.'},
 {'claim_id': 'ctf-phase-flipping',
  'chapter': '05_image_formation',
  'anchor': 'ctf-correction-example',
  'claim': 'Phase flipping multiplies Fourier coefficients by sign(CTF), correcting sign inversions but not '
           'restoring attenuated amplitudes.',
  'source_id': 'bendory2020',
  'locator': 'p. 69, Eqs. (16)-(17)',
  'scope': 'A basic CTF correction, not Wiener amplitude correction.'},
 {'claim_id': 'ctf-model',
  'chapter': '07_synthetic_data',
  'anchor': 'ctf',
  'claim': 'The CTF phase depends on wavelength, defocus, spherical aberration, spatial frequency, and amplitude '
           'or phase contrast terms.',
  'source_id': 'singer2020',
  'locator': 'pp. 172-174',
  'scope': 'Conventional weak-phase, thin-specimen CTF model.'},
 {'claim_id': 'dft-dc',
  'chapter': '03_fourier',
  'anchor': 'dft',
  'claim': 'For the unnormalized forward DFT, the zero-frequency coefficient equals the sum of image samples.',
  'source_id': 'szeliski2022',
  'locator': 'Ch. 3.4, pp. 143-146',
  'scope': 'Normalization conventions change the scale, not the zero-frequency interpretation.'},
 {'claim_id': 'convolution-theorem',
  'chapter': '03_fourier',
  'anchor': 'convolution',
  'claim': 'Spatial convolution corresponds to multiplication in the Fourier domain.',
  'source_id': 'szeliski2022',
  'locator': 'Ch. 3.2, p. 122',
  'scope': 'Boundary assumptions distinguish circular from linear discrete convolution.'},
 {'claim_id': 'nyquist-limit',
  'chapter': '01_image_basics',
  'anchor': 'sampling',
  'claim': 'The Nyquist frequency is one half-cycle per sample, so a pixel size a implies an ideal sampling limit '
           'of resolution 2a.',
  'source_id': 'szeliski2022',
  'locator': 'Ch. 2.3.1, p. 84',
  'scope': 'A sampling bound, not a guarantee of achieved experimental resolution.'},
 {'claim_id': 'padding-boundaries',
  'chapter': '02_filter_segment',
  'anchor': 'filtering',
  'claim': 'Padding choice changes filter behavior at image boundaries and can create edge artifacts.',
  'source_id': 'szeliski2022',
  'locator': 'Ch. 3.2, pp. 123-124',
  'scope': 'Finite digital images.'},
 {'claim_id': 'dog-picking-bias',
  'chapter': '06_workflow',
  'anchor': 'workflow',
  'claim': 'Difference-of-Gaussians can provide a low-bias generic particle-picking template, whereas detailed '
           'references can introduce model bias.',
  'source_id': 'sigworth2016',
  'locator': 'p. 62',
  'scope': 'A teaching baseline, not a claim that DoG is universally optimal.'},
 {'claim_id': 'wavelet-perfect-reconstruction',
  'chapter': '04_wavelet',
  'anchor': 'dwt-perfect-reconstruction',
  'claim': 'A matched analysis and synthesis discrete wavelet filter bank can reconstruct its input exactly up to '
           'numerical precision.',
  'source_id': 'szeliski2022',
  'locator': 'Ch. 3.5, pp. 159-163',
  'scope': 'Compatible boundary mode and complete coefficient set are required.'},
 {'claim_id': 'spa-low-snr',
  'chapter': '05_background',
  'anchor': 'low-dose-snr',
  'claim': 'Low dose and weak contrast make single-particle cryo-EM images extremely noisy, so information is '
           'recovered statistically across many particles.',
  'source_id': 'bendory2020',
  'locator': 'pp. 60-63',
  'scope': 'Single-particle cryo-EM.'},
 {'claim_id': 'ded-movies',
  'chapter': '05_image_formation',
  'anchor': 'detectors',
  'claim': 'Direct electron detectors enabled dose-fractionated movies and computational correction of '
           'beam-induced motion.',
  'source_id': 'mie2016-ch01',
  'locator': 'pp. 1-13',
  'scope': 'Historical and practical detector context as of 2016.'},
 {'claim_id': 'motion-translation-limit',
  'chapter': '06_workflow',
  'anchor': 'motion-correction',
  'claim': 'The compared movie-alignment methods model in-plane translations and do not fully model rotation or '
           'out-of-plane motion.',
  'source_id': 'mie2016-ch05',
  'locator': 'pp. 103-122',
  'scope': 'Methods reviewed in the 2016 chapter.'},
 {'claim_id': 'workflow-stages',
  'chapter': '06_workflow',
  'anchor': 'workflow',
  'claim': 'The SPA workflow separates micrograph-level motion, CTF, and picking steps from particle-level '
           'classification and reconstruction.',
  'source_id': 'bendory2020',
  'locator': 'pp. 60-61, Fig. 4',
  'scope': 'Conceptual workflow; software packages may order or combine stages differently.'},
 {'claim_id': 'hard-soft-assignment',
  'chapter': '06_reconstruction_validation',
  'anchor': 'soft-assignment',
  'claim': 'Refinement may use hard best-match assignments or soft likelihood-based weights over orientations and '
           'shifts.',
  'source_id': 'bendory2020',
  'locator': 'pp. 64-66',
  'scope': 'Conceptual distinction; implementations add priors and regularization.'},
 {'claim_id': 'ml-latent-marginalization',
  'chapter': '06_reconstruction_validation',
  'anchor': 'soft-assignment',
  'claim': 'Maximum-likelihood refinement marginalizes nuisance variables such as orientation and translation '
           'rather than estimating a single fixed value first.',
  'source_id': 'mie2010-ch10',
  'locator': 'pp. 263-294',
  'scope': 'Likelihood-based cryo-EM image processing.'},
 {'claim_id': 'gold-standard-halves',
  'chapter': '06_reconstruction_validation',
  'anchor': 'validation',
  'claim': 'Gold-standard refinement keeps two independently refined half-sets to reduce overfitting and assesses '
           'agreement with FSC.',
  'source_id': 'mie2016-ch06',
  'locator': 'pp. 129-131',
  'scope': 'Independence must be preserved through refinement.'},
 {'claim_id': 'fsc-mask-caveat',
  'chapter': '06_reconstruction_validation',
  'anchor': 'validation',
  'claim': 'Masks can inflate FSC, so map validation must test masking and high-resolution noise behavior.',
  'source_id': 'mie2016-ch09',
  'locator': 'pp. 239-245',
  'scope': 'Half-map resolution validation.'},
 {'claim_id': 'orientation-bias',
  'chapter': '06_reconstruction_validation',
  'anchor': 'fourier-coverage',
  'claim': 'Preferred orientations leave anisotropic Fourier coverage and therefore anisotropic resolution.',
  'source_id': 'sigworth2016',
  'locator': 'pp. 59-60',
  'scope': 'Projection data with nonuniform viewing directions.'},
 {'claim_id': 'heterogeneity-identifiability',
  'chapter': '06_heterogeneity',
  'anchor': 'discrete-heterogeneity',
  'claim': 'Low SNR can make separation of distinct conformational populations unreliable with finite data.',
  'source_id': 'sigworth2016',
  'locator': 'pp. 63-64',
  'scope': 'Finite-sample discrimination; this does not assert that every Gaussian mixture becomes theoretically '
           'unidentifiable at low SNR.'},
 {'claim_id': 'simulation-noise-boundary',
  'chapter': '07_synthetic_data',
  'anchor': 'noise',
  'claim': 'Cryo-EM image noise is often modeled as additive Gaussian after frame integration, but its empirical '
           'spectrum need not be white.',
  'source_id': 'bendory2020',
  'locator': 'pp. 62-63',
  'scope': "The chapter's global iid AWGN is an explicit teaching simplification, not a realistic detector-noise "
           'claim.'},
 {'claim_id': 'spa-sample-optimization',
  'chapter': '05_background',
  'anchor': 'sample-preparation',
  'claim': 'Specimen quality and processing inform each other; ice and orientation distribution affect '
           'recoverable information.',
  'source_id': 'cheng2015primer',
  'locator': 'pp.439–440, Fig.1 and Specimen Preparation',
  'scope': 'Teaching overview grounded in the 2015 primer; not a current software capability claim.',
  'verification_status': 'verified',
  'verification_basis': 'Manual reading of the cited local PDF passages, 2026-09-11'},
 {'claim_id': 'spa-colored-noise',
  'chapter': '05_image_formation',
  'anchor': 'noise-model',
  'claim': 'Frequency-dependent Gaussian variances permit colored noise while maintaining conditional '
           'independence assumptions in Fourier coordinates.',
  'source_id': 'scheres2012bayes',
  'locator': 'p.410, Eq.(6)',
  'scope': 'Statistical noise approximation, not proof of empirical independence.',
  'verification_status': 'verified',
  'verification_basis': 'Manual reading of the cited local PDF passages, 2026-09-11'},
 {'claim_id': 'spa-em-lower-bound',
  'chapter': '05_statistical_inference',
  'anchor': 'spa-em',
  'claim': 'The posterior over latent variables makes the Jensen lower bound tight at the current parameters.',
  'source_id': 'ma2019em',
  'locator': 'PDF pp.4–5, Eqs.(6–9)',
  'scope': 'Exact EM increases likelihood under a valid M-step; no global optimum guarantee.',
  'verification_status': 'verified',
  'verification_basis': 'Manual reading of the cited local PDF passages, 2026-09-11'},
 {'claim_id': 'spa-motion-frequency-evidence',
  'chapter': '06_workflow',
  'anchor': 'motion-correction',
  'claim': 'Motion correction may improve Thon rings even when raw and corrected real-space image sums look '
           'similar.',
  'source_id': 'cheng2015primer',
  'locator': 'p.440, Fig.2',
  'scope': 'Evidence from the example in the primer, not a guarantee for every movie.',
  'verification_status': 'verified',
  'verification_basis': 'Manual reading of the cited local PDF passages, 2026-09-11'},
 {'claim_id': 'spa-soft-2d-alignment',
  'chapter': '06_alignment_classification',
  'anchor': 'soft-alignment',
  'claim': 'The 2D ML update averages transformed images using normalized weights over possible in-plane poses.',
  'source_id': 'sigworth1998',
  'locator': 'pp.331–332, Eqs.(12),(16),(20)',
  'scope': 'Single underlying 2D view with white Gaussian noise and pose prior.',
  'verification_status': 'verified',
  'verification_basis': 'Manual reading of the cited local PDF passages, 2026-09-11'},
 {'claim_id': 'spa-robust-cl2d',
  'chapter': '06_alignment_classification',
  'anchor': 'cl2d',
  'claim': 'CL2D combines correntropy with divisive clustering and a class-relative robust assignment criterion.',
  'source_id': 'sorzano2010',
  'locator': 'pp.198–201, Sections 2.1–2.4',
  'scope': 'Correntropy is a similarity measure, not a posterior probability.',
  'verification_status': 'verified',
  'verification_basis': 'Manual reading of the cited local PDF passages, 2026-09-11'},
 {'claim_id': 'spa-isac-stability',
  'chapter': '06_alignment_classification',
  'anchor': 'isac',
  'claim': 'Similar repeated class averages and FRC do not by themselves establish within-class homogeneity.',
  'source_id': 'yang2012-supp',
  'locator': 'PDF p.2, Section S2, Figs.S8–S10',
  'scope': 'The supplemental experiment supports using per-particle alignment stability; not a universal purity '
           'guarantee.',
  'verification_status': 'verified',
  'verification_basis': 'Manual reading of the cited local PDF passages, 2026-09-11'},
 {'claim_id': 'spa-gsup-updates',
  'chapter': '06_alignment_classification',
  'anchor': 'gamma-sup',
  'claim': 'Gamma-SUP updates both centers and their representative data through a blurring self-updating process '
           'with local influence weights.',
  'source_id': 'chen2014',
  'locator': 'Sections 2.3,3.2, Eqs.(10–19)',
  'scope': 'The compact-support regime assumes q<1 and gamma>1-q; it is not ordinary Gaussian-mixture EM.',
  'verification_status': 'verified',
  'verification_basis': 'Manual reading of the cited local PDF passages, 2026-09-11'},
 {'claim_id': 'spa-2sdr-ranks',
  'chapter': '06_dimension_reduction',
  'anchor': 'spa-2sdr',
  'claim': '2SDR combines MPCA rank reduction with PCA of the compressed coefficients, using SURE then GIC for '
           'rank selection.',
  'source_id': 'chung2020',
  'locator': 'pp.293–297, Section 2, Eqs.(2.1–2.5), Fig.1',
  'scope': 'Rank consistency depends on model and theorem assumptions; hand-chosen ranks are a teaching '
           'simplification.',
  'verification_status': 'verified',
  'verification_basis': 'Manual reading of the cited local PDF passages, 2026-09-11'},
 {'claim_id': 'spa-map-regularization',
  'chapter': '06_reconstruction_validation',
  'anchor': 'reconstruction',
  'claim': 'The MAP reconstruction update adds prior precision to a CTF/noise-weighted reconstruction '
           'denominator.',
  'source_id': 'scheres2012bayes',
  'locator': 'p.410, Eqs.(8–12)',
  'scope': 'Gaussian prior regularizes uncertain frequencies; no recovery guarantee at unobserved frequencies.',
  'verification_status': 'verified',
  'verification_basis': 'Manual reading of the cited local PDF passages, 2026-09-11'},
 {'claim_id': 'spa-sgd-marginalization',
  'chapter': '06_reconstruction_validation',
  'anchor': 'soft-assignment',
  'claim': 'The cryoSPARC ab initio objective marginalizes class and pose before stochastic gradient updates.',
  'source_id': 'punjani2017-supp',
  'locator': 'PDF pp.1–4, Supplementary Note 1, Eqs.(1–9)',
  'scope': 'Uniform class probabilities and pose prior are settings of this paper, not universal assumptions.',
  'verification_status': 'verified',
  'verification_basis': 'Manual reading of the cited local PDF passages, 2026-09-11'},
 {'claim_id': 'spa-3dva-fixed-poses',
  'chapter': '06_heterogeneity',
  'anchor': 'linear-heterogeneity',
  'claim': '3DVA fits linear variability with a supplied consensus map, poses and CTFs; the implementation uses '
           'the PPCA limiting least-squares updates.',
  'source_id': 'punjani2021',
  'locator': 'PDF p.12, Section 5.2, Eqs.(4–8)',
  'scope': 'Latent coordinates are not time; pose errors may affect variability estimates.',
  'verification_status': 'verified',
  'verification_basis': 'Manual reading of the cited local PDF passages, 2026-09-11'},
 {'claim_id': 'spa-shared-bias-fsc',
  'chapter': '06_resolution_validation',
  'anchor': 'fsc',
  'claim': 'Shared systematic errors can inflate agreement between half maps, so FSC does not by itself certify '
           'structural correctness.',
  'source_id': 'sorzano2022',
  'locator': 'pp.419–420, Section 5.1, Eqs.(3–4)',
  'scope': 'This motivates complementary bias checks, not abandoning independent half-set validation.',
  'verification_status': 'verified',
  'verification_basis': 'Manual reading of the cited local PDF passages, 2026-09-11'}]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def write_jsonl(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n" for row in rows)
    path.write_text(text, encoding="utf-8")


def live_pdf_state(source_root: Path) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for path in sorted(source_root.rglob("*.pdf")):
        if "document_cache" in path.parts:
            continue
        rel = path.relative_to(source_root).as_posix()
        result[rel] = {"sha256": sha256_file(path), "size_bytes": path.stat().st_size}
    return result


def cache_state(index_path: Path) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    rows = read_jsonl(index_path)
    return rows, {row["relative_path"]: row for row in rows}


def _chapter_metadata(relative_path: str) -> dict[str, str] | None:
    parts = Path(relative_path).parts
    if len(parts) != 2 or parts[0] not in {"2010", "2016"}:
        return None
    for (year, token), metadata in CHAPTERS.items():
        if parts[0] == year and f"{token}---" in parts[1]:
            return metadata
    return None


def classify(relative_path: str, sha: str, canonical_for_sha: dict[str, str]) -> dict[str, Any]:
    if canonical_for_sha[sha] != relative_path:
        return {
            "source_id": "duplicate-" + sha[:12],
            "title": Path(relative_path).stem,
            "content_role": "duplicate_alias",
            "coverage_status": None,
            "digest_path": None,
            "citation_key": None,
        }
    metadata = _chapter_metadata(relative_path) or STANDALONE.get(relative_path)
    if metadata:
        return {**metadata, "content_role": "substantive"}
    # Only known volume paratext is front matter. Unknown PDFs must be routed.
    name = Path(relative_path).name
    is_front = len(Path(relative_path).parts) == 2 and Path(relative_path).parts[0] in {"2010", "2016"} and any(
        name.startswith(token + "_") for token in
        ("Author-Index", "Contributors", "Copyright", "Preface", "Subject-Index", "Title-Page", "Volume-in-Series", "Series-Page")
    )
    return {
        "source_id": ("frontmatter-" if is_front else "unrouted-") + sha[:12],
        "title": Path(relative_path).stem.replace("_", " "),
        "content_role": "front_matter" if is_front else "unrouted",
        "coverage_status": None,
        "digest_path": None,
        "citation_key": None,
    }


def public_registry(index_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_sha: dict[str, list[str]] = defaultdict(list)
    for row in index_rows:
        by_sha[row["sha256"]].append(row["relative_path"])
    # Prefer the volume copy over top-level duplicates; otherwise lexical order.
    canonical = {sha: sorted(paths, key=lambda p: ("/" not in p, p))[0] for sha, paths in by_sha.items()}
    records = []
    for row in sorted(index_rows, key=lambda r: r["relative_path"]):
        routing = classify(row["relative_path"], row["sha256"], canonical)
        aliases = sorted(p for p in by_sha[row["sha256"]] if p != row["relative_path"])
        records.append({
            "schema_version": SCHEMA_VERSION,
            "source_id": routing["source_id"],
            "title": routing["title"],
            "relative_path": row["relative_path"],
            "sha256": row["sha256"],
            "pages": row.get("pages"),
            "size_bytes": row.get("size_bytes"),
            "extraction_status": row.get("status"),
            "duplicate_aliases": aliases,
            "content_role": routing["content_role"],
            "coverage_status": routing["coverage_status"],
            "digest_path": routing["digest_path"],
            "citation_key": routing["citation_key"],
            "topic_hints": sorted(row.get("keyword_hits", [])),
        })
    return records


def seed_claims(registry: list[dict[str, Any]], path: Path, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        return
    hashes = {r["source_id"]: r["sha256"] for r in registry if r["content_role"] == "substantive"}
    rows = []
    for claim in CLAIMS:
        row = dict(claim)
        row.update({"schema_version": SCHEMA_VERSION, "source_sha256": hashes[claim["source_id"]], "verification_status": claim.get("verification_status", "needs_review")})
        rows.append(row)
    write_jsonl(path, rows)


def write_snapshot(path: Path, registry: list[dict[str, Any]], extraction_action: str) -> None:
    unique = {r["sha256"] for r in registry}
    roles = Counter(r["content_role"] for r in registry)
    coverage = Counter(r["coverage_status"] for r in registry if r["content_role"] == "substantive")
    digest = hashlib.sha256("\n".join(sorted(unique)).encode()).hexdigest()
    text = f"""# Reference ingestion snapshot

This is a public metadata snapshot. Extracted full text remains in the local
`References/document_cache/` and is intentionally not published.

- Schema version: {SCHEMA_VERSION}
- Pipeline version: {PIPELINE_VERSION}
- Input PDF records: {len(registry)}
- Unique SHA-256 values: {len(unique)}
- Exact duplicate groups: {len({r['sha256'] for r in registry if r['duplicate_aliases']})}
- Successful extractions: {sum(r['extraction_status'] == 'ok' for r in registry)}
- Substantive unique sources: {roles['substantive']}
- Front matter records: {roles['front_matter']}
- Duplicate aliases: {roles['duplicate_alias']}
- Coverage: {coverage['full_digest']} `full_digest`, {coverage['route_note']} `route_note`, {coverage['out_of_scope']} `out_of_scope`
- Extraction action: `{extraction_action}`
- Snapshot fingerprint: `{digest}`

The registry contains titles, relative source names, hashes, extraction state,
and routing decisions only. It excludes extracted prose and machine-specific
absolute paths.
"""
    path.write_text(text, encoding="utf-8")


def write_asset_manifest(path: Path, repo_root: Path) -> None:
    image_root = repo_root / "book" / "images"
    content_files = list((repo_root / "book").glob("*.md")) + list((repo_root / "book").glob("*.py"))
    content = "\n".join(p.read_text(encoding="utf-8", errors="replace") for p in content_files)
    lines = ["schema_version: 1", "policy: public_metadata_only", "assets:"]
    for image_path in sorted(p for p in image_root.rglob("*") if p.is_file()):
        rel = image_path.relative_to(repo_root).as_posix()
        book_rel = image_path.relative_to(repo_root / "book").as_posix()
        used = book_rel in content or image_path.name in content
        lines.extend([
            f"  - path: {json.dumps(rel)}",
            f"    sha256: {sha256_file(image_path)}",
            f"    size_bytes: {image_path.stat().st_size}",
            f"    referenced: {'true' if used else 'false'}",
            "    provenance: project_provided",
            "    license_status: needs_review",
            "    alt_text_status: needs_review",
        ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def resolve_defaults(args: argparse.Namespace) -> None:
    repo_root = Path(args.repo_root).resolve() if args.repo_root else Path(__file__).resolve().parents[1]
    source_root = Path(args.source_root).resolve() if args.source_root else repo_root.parent / "References"
    args.repo_root = repo_root
    args.source_root = source_root
    args.cache_root = Path(args.cache_root).resolve() if args.cache_root else source_root / "document_cache"
    args.catalog_dir = Path(args.catalog_dir).resolve() if args.catalog_dir else repo_root / "notes" / "reference_catalog"
    args.digest_dir = Path(args.digest_dir).resolve() if args.digest_dir else repo_root / "notes" / "ref_digest"
    args.bib = Path(args.bib).resolve() if args.bib else repo_root / "book" / "references.bib"


def upstream_ingest(args: argparse.Namespace) -> None:
    script = args.ingest_script or os.environ.get("DOCUMENT_INGEST_SCRIPT") or shutil.which("ingest_documents.py")
    if not script:
        raise RuntimeError("Source/cache SHA mismatch. Pass --ingest-script or set DOCUMENT_INGEST_SCRIPT to refresh extracted text.")
    # The generic CLI scans all extensions, including an existing in-tree cache.
    # Reuse its extraction API with this project's explicit PDF-only source list.
    import importlib.util
    spec = importlib.util.spec_from_file_location("primer_document_ingest", script)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load document ingestion script: {script}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    text_dir = args.cache_root / "text"
    text_dir.mkdir(parents=True, exist_ok=True)
    records = [module.ingest_one(args.source_root / relative, args.source_root, text_dir,
               ["cryo-em", "likelihood", "classification", "heterogeneity"])
               for relative in sorted(live_pdf_state(args.source_root))]
    module.write_outputs(records, args.source_root, args.cache_root)


def ingest(args: argparse.Namespace) -> int:
    resolve_defaults(args)
    index_path = args.cache_root / "references" / "document_index.jsonl"
    if not index_path.exists():
        upstream_ingest(args)
    live = live_pdf_state(args.source_root)
    index_rows, cached = cache_state(index_path)
    mismatch = set(live) != set(cached) or any(live[p]["sha256"] != cached[p]["sha256"] for p in set(live) & set(cached))
    action = "skipped_up_to_date"
    if mismatch or args.force:
        upstream_ingest(args)
        index_rows, cached = cache_state(index_path)
        action = "refreshed"
    args.catalog_dir.mkdir(parents=True, exist_ok=True)
    registry = public_registry(index_rows)
    write_jsonl(args.catalog_dir / "source_registry.jsonl", registry)
    seed_claims(registry, args.catalog_dir / "claim_map.jsonl", overwrite=args.overwrite_claims)
    write_snapshot(args.catalog_dir / "ingestion_snapshot.md", registry, action)
    write_asset_manifest(args.catalog_dir / "asset_manifest.yml", args.repo_root)
    report = validate_paths(args, emit=False)
    print(json.dumps({"extraction_action": action, **report}, ensure_ascii=False, sort_keys=True))
    return 0 if not report["errors"] else 1


def _digest_file(repo_root: Path, digest_path: str) -> Path:
    return repo_root / digest_path.split("#", 1)[0]


def validate_paths(args: argparse.Namespace, emit: bool = True) -> dict[str, Any]:
    index_path = args.cache_root / "references" / "document_index.jsonl"
    registry_path = args.catalog_dir / "source_registry.jsonl"
    claims_path = args.catalog_dir / "claim_map.jsonl"
    errors: list[str] = []
    warnings: list[str] = []
    needs_reverify: list[str] = []
    if not index_path.exists() or not registry_path.exists() or not claims_path.exists():
        missing = [str(p) for p in (index_path, registry_path, claims_path) if not p.exists()]
        result = {"errors": ["missing required files: " + ", ".join(missing)], "warnings": [], "needs_reverify": [], "pdf_records": 0, "unique_sha256": 0, "substantive_sources": 0, "claims": 0, "assets": 0}
        if emit:
            print(json.dumps(result, ensure_ascii=False, indent=2))
        return result
    registry = read_jsonl(registry_path)
    claims = read_jsonl(claims_path)
    live = live_pdf_state(args.source_root)
    registered = {r["relative_path"]: r for r in registry}
    if set(live) != set(registered):
        errors.append("source_registry paths differ from live PDF paths")
    for relative_path in sorted(set(live) & set(registered)):
        if live[relative_path]["sha256"] != registered[relative_path]["sha256"]:
            errors.append(f"source SHA changed: {relative_path}")
    for row in registry:
        if row["content_role"] == "unrouted":
            errors.append(f"unrouted PDF requires explicit review: {row['relative_path']}")
    substantive = [r for r in registry if r["content_role"] == "substantive"]
    ids = [r["source_id"] for r in substantive]
    if len(ids) != len(set(ids)):
        errors.append("substantive source_id values are not unique")
    for row in substantive:
        if row["coverage_status"] not in ROUTING_STATUSES:
            errors.append(f"invalid coverage status: {row['source_id']}")
        if not row["digest_path"] or not _digest_file(args.repo_root, row["digest_path"]).exists():
            errors.append(f"missing digest: {row['source_id']} -> {row['digest_path']}")
    hashes = {r["source_id"]: r["sha256"] for r in substantive}
    claim_ids: set[str] = set()
    for claim in claims:
        cid = claim.get("claim_id", "<missing>")
        if cid in claim_ids:
            errors.append(f"duplicate claim_id: {cid}")
        claim_ids.add(cid)
        sid = claim.get("source_id")
        if sid not in hashes:
            errors.append(f"claim {cid} references unknown source_id {sid}")
        elif claim.get("source_sha256") != hashes[sid]:
            needs_reverify.append(cid)
        for field in ("chapter", "anchor", "claim", "locator", "scope"):
            if not claim.get(field):
                errors.append(f"claim {cid} missing {field}")
    if needs_reverify:
        warnings.append(f"{len(needs_reverify)} claims need re-verification because a source SHA changed")
    if args.bib.exists():
        bibkeys = set(re.findall(r"@[A-Za-z]+\{([^,]+),", args.bib.read_text(encoding="utf-8")))
        used = {r["citation_key"] for r in substantive if r.get("citation_key")}
        missing_keys = sorted(used - bibkeys)
        if missing_keys:
            errors.append("missing BibTeX keys: " + ", ".join(missing_keys))
    else:
        errors.append(f"missing bibliography: {args.bib}")
    asset_path = args.catalog_dir / "asset_manifest.yml"
    assets = 0
    if asset_path.exists():
        assets = sum(1 for line in asset_path.read_text(encoding="utf-8").splitlines() if line.startswith("  - path:"))
    else:
        errors.append(f"missing asset manifest: {asset_path}")
    result = {
        "errors": errors,
        "warnings": warnings,
        "needs_reverify": needs_reverify,
        "pdf_records": len(registry),
        "unique_sha256": len({r["sha256"] for r in registry}),
        "substantive_sources": len(substantive),
        "coverage": dict(sorted(Counter(r["coverage_status"] for r in substantive).items())),
        "claims": len(claims),
        "assets": assets,
    }
    if emit:
        print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    return result


def validate(args: argparse.Namespace) -> int:
    resolve_defaults(args)
    result = validate_paths(args)
    return 0 if not result["errors"] and not result["needs_reverify"] else 1


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=("ingest", "validate"))
    p.add_argument("--repo-root")
    p.add_argument("--source-root")
    p.add_argument("--cache-root")
    p.add_argument("--catalog-dir")
    p.add_argument("--digest-dir")
    p.add_argument("--bib")
    p.add_argument("--ingest-script", help="Path to the generic document-folder-ingest script")
    p.add_argument("--force", action="store_true", help="Refresh the generic extraction cache even when SHA state matches")
    p.add_argument("--overwrite-claims", action="store_true", help="Re-seed claim hashes after claims have been manually verified")
    return p


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        return ingest(args) if args.command == "ingest" else validate(args)
    except (OSError, RuntimeError, subprocess.CalledProcessError, ValueError) as exc:
        print(f"reference_pipeline: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
