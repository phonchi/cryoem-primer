# Original-notebook restoration, 2026-09-13

Approved: restore chapters 01–04 to the four parent-directory original notebooks, their original example images, parameter values and teaching sequence. Correct verified errors and APIs, annotate changes. Never overwrite originals. Never change 07_synthetic_data.py/ipynb or existing image/data assets. No SPA content or SPA cross-links in general-image chapters; root moves required material into SPA. Preserve current page URLs. Complete prose in Traditional Chinese, use speak-human-tw; core explanations/body, longer derivations/necessary extra material dropdowns. No comprehension-check sections. Link-only original topics become readable concise teaching sections, not new giant labs. Restore all original diagrams and original links with contextual readings.

Website interaction is client-side JavaScript/Canvas/SVG, no Python server and no paid service. Existing book typography/light-dark theme, blue inputs / orange active kernel / green output. Original illustrations and data preserved. Python is numerical source of truth. A source-cell-to-section/example coverage ledger is mandatory for each chapter.

Shared code contract owned by root:
```
from pathlib import Path
import sys
BOOK = Path('book') if Path('book').is_dir() else Path('.')
sys.path.insert(0, str(BOOK.resolve()))
from _support import image_path, show_images, lab
```
`image_path(name)` returns local original asset path (filename names unchanged, parent download URLs via requests).
`show_images(*images, titles=None, cmap=None, figsize=None, shared=False)` makes/shows a matplotlib figure, per-image normalization by default, preserve RGB colors. Returns figure.
`lab(kind, **data)` displays a self-contained iframe widget, accepting NumPy arrays and scalar/list data. root implements engine; all widgets have unique IDs, play/pause/reset and mobile layout. Static Python figures remain present and provide fallback.
Kinds:
- geometry: image=venusaurf (RGBA uint8), rgb=venusaurf2, text_image=ski.data.text(); UI original presets and matrix/center/angle/translation/scale/shear/projective controls. order=5 comparisons stay exact Python images in chapter.
- convolution1d: signal=original noisy_signal, kernels={"3-point":list(...),"11-point":list(...)}. Steps valid convolution, padding comparison separate Python original.
- convolution2d: image=original 7x7 array, kernel=3x3 mean array. exact scan.
- fourier_basis: n=128, kx=5, ky=2, amplitude=1; bounds [-32,32], sine/cosine; original FFT coefficient normalization, DC handled.
- filter_sweep: original=normalized grayscale Gengar, cutoffs=list(range(1,25)), lowpass=list of real float result arrays, highpass=list of signed result arrays; source retains exact per-cutoff original/corrected mask computation; iframe shows cutoff mask and spectrum if supplied `spectrum=log1p(abs(fftshift(fft2(original))))`.
Root adds comparison/diagram viewers where useful. Do not invent other lab kinds without message.

Agents own individual chapter .py and its ledger only; root owns shared support, assets, notebook syncing, tests, SPA migration, diagrams and frontend. Do not sync notebooks/build globally from agent. No broad search/code-review gates; independent science and source checks then two reader reviews, browser QA, commit/push/deployment verification authorized. Tests must be adapted to original examples, not preserve substitute camera experiments.
