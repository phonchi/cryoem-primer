# ASPIRE API errata: `Simulation_data.ipynb` vs aspire==0.14.3

Verified against the installed `aspire 0.14.3` in the `cryoem-book` conda env,
by reading the installed package source directly (not from memory), and
cross-checked against `/home/phonchi/2sdr_research/synthetic_validation/generate_synthetic.py`,
which pins `ASPIRE_VERSION = "0.14.3"` and is a known-working reference for
this exact pipeline (project, volume, CTF, simulate, save).

A full 7-stage smoke test with these fixes applied ran clean end-to-end
(`notes/smoke_test.py`, n=100, img_size=130). **Result: PASSED.**

Each entry: original notebook code → what breaks under 0.14.3 → the fix used
in `smoke_test.py` → why.

---

## 1. Dead import: `uniform_random_angles` no longer exists

```python
# Notebook (cell 1), breaks immediately:
from aspire.utils.coor_trans import uniform_random_angles
```

`ImportError: cannot import name 'uniform_random_angles' from 'aspire.utils.coor_trans'`.

This function was never actually called anywhere in the notebook (the
notebook builds its Euler angles by hand in cell 9) — it was just dead
scaffolding for a comment. In 0.14.3, random-orientation generation lives on
`aspire.utils.Rotation.generate_random_rotations(...)`, used internally by
`Simulation._init_angles()` when `angles=None`.

**Fix:** delete the import line entirely. No replacement needed since the
notebook never calls it.

---

## 2. `RadialCTFFilter(pixel_size, voltage, defocus=d, Cs=2.0, alpha=0.1)` — wrong signature, raises TypeError

```python
# Notebook (cell 10):
filters = [RadialCTFFilter(pixel_size, voltage, defocus=d, Cs=2.0, alpha=0.1)
           for d in np.linspace(defocus_min, defocus_max, defocus_ct)]
```

0.14.3's signature is:

```python
RadialCTFFilter(voltage=200, defocus=15000, Cs=2.26, alpha=0.07, B=0)
```

There is **no `pixel_size` parameter at all** anymore. The notebook's two
positional args (`pixel_size`, `voltage`) land on `voltage` and `defocus`
respectively; then the explicit `defocus=d` keyword collides with the
already-filled `defocus` positional slot →
`TypeError: RadialCTFFilter.__init__() got multiple values for argument 'defocus'`.

**Fix:**

```python
filters = [RadialCTFFilter(voltage=voltage, defocus=d, Cs=Cs, alpha=alpha)
           for d in np.linspace(defocus_min, defocus_max, defocus_ct)]
```

`pixel_size` moved out of the per-filter constructor and onto the
`Simulation`/`ImageSource` itself — see item 3.

---

## 3. `Simulation(...)` without `pixel_size=` does NOT crash — it silently uses the wrong pixel size (1.0 A instead of 2.82 A)

**Correction (verified empirically, see below): this is not a crash bug.**
An earlier pass of this errata claimed omitting `pixel_size` raises a
`RuntimeError` from `CTFFilter.evaluate()`. That's wrong — checked directly
against the installed package and it is silent, which is actually worse for
a teaching notebook.

The notebook's `Simulation(...)` call (cell 11) never sets `pixel_size`.
Reading `aspire/source/simulation.py` (lines 114-123):

```python
# Infer pixel_size from volume. Otherwise default to 1.0 angstrom.
if pixel_size is None:
    if self.vols.pixel_size is None:
        pixel_size = 1.0
    else:
        pixel_size = self.vols.pixel_size
```

The notebook builds its `Volume` as `Volume(vols)` (cell 7) from a bare
`mrcfile` array — no `pixel_size` is ever attached, even though the source
`.mrc` header may carry voxel spacing. So `vols.pixel_size` is `None`, and
`Simulation` silently falls back to **`pixel_size=1.0`** angstrom instead of
the intended `2.82`. No exception anywhere — `sim.pixel_size` is simply
`1.0`, and every downstream CTF evaluation
(`aspire/operators/filters.py`, `RadialCTFFilter._evaluate`, which computes
bandwidth as `BW = 1/(pixel_size/10)`) uses that wrong value.

**Verified impact:** built two otherwise-identical `Simulation` objects
(same seed, same 5 `RadialCTFFilter`s, `n=20`), one with no `pixel_size`
(→ 1.0 default) and one with `pixel_size=2.82`, and diffed
`sim.clean_images[:20]`. Result: **relative L2 difference of ~152%,
max abs pixel difference 0.22** — the CTF envelope is completely different
between the two, not a subtle rounding difference. The notebook as literally written never reaches this bug — it crashes first
at cell 10 (item 2, `RadialCTFFilter` signature). But with items 1, 2, 4, 5,
6 fixed and this one left alone, the pipeline runs clean end-to-end and
silently generates CTF-corrupted images that do not correspond to the
stated `pixel_size = 2.82` experimental parameter, with no error or warning
to flag it. That's the dangerous case: anyone who patches the crashes
without knowing about this one ships wrong figures.

**Fix:** add `pixel_size=pixel_size` to the `Simulation(...)` call, matching
`generate_synthetic.py`'s `simulation = api["Simulation"](..., pixel_size=config.pixel_size)`.

```python
sim = Simulation(
    L=img_size, n=num_imgs, vols=v, C=num_maps,
    unique_filters=filters, angles=my_angles,
    offsets=0.0, amplitudes=1.0,
    pixel_size=pixel_size,   # silently defaults to 1.0 A if omitted (see above)
)
```

(Passing `pixel_size=pixel_size` to `Volume(...)` at construction, cell 7,
would also fix it via inheritance — either the Volume or the Simulation
needs it, but the notebook currently sets neither.)

---

## 4. `sim.projections(start=0, num=num_imgs)` — no longer callable

```python
# Notebook (cell 13):
imgs_clean = sim.projections(start=0, num=num_imgs)
```

`projections` is now a `@property` returning a subscriptable `_ImageAccessor`,
not a bound method. Calling it raises:

```
TypeError: '_ImageAccessor' object is not callable
```

**Fix:** slice instead of calling:

```python
imgs_clean = sim.projections[:num_imgs]
```

(`sim.projections[:]` also works and defaults to all `n` images.) The same
`_ImageAccessor` pattern applies to `.images` and `.clean_images` (see items
5 and 6) — this is a systematic API shift, not a one-off.

---

## 5. `sim.eval_filters(imgs_clean)` — method removed entirely

```python
# Notebook (cell 15):
imgs_ctf_clean = sim.eval_filters(imgs_clean)
```

`Simulation.eval_filters` does not exist anywhere in the installed 0.14.3
package (confirmed by grepping the whole package tree — zero hits). This
raises `AttributeError: 'Simulation' object has no attribute 'eval_filters'`.

**Fix:** use the `clean_images` accessor, which returns projections with
filters/shifts/amplitudes applied but no noise — exactly what
`eval_filters(projections())` used to produce:

```python
imgs_ctf_clean = sim.clean_images[:num_imgs]
```

No separate call to `sim.projections` is even needed if only the CTF'd
images are wanted; the notebook computes `imgs_clean` (stage 6a) purely for
the mid-pipeline visualization in cell 14, so both accessors are kept in the
smoke test for parity but `clean_images` alone would suffice for the final
output.

---

## 6. `relion_src._images()` — private method now requires `indices`

```python
# Notebook (cell 24):
imgs = relion_src._images()
```

`_images` is the internal method backing the public `.images` accessor and
now requires an explicit 1-D index array:

```
TypeError: _images() missing 1 required positional argument: 'indices'
```

**Fix:** use the public, subscriptable accessor instead:

```python
imgs = relion_src.images[:]
```

---

## 7. `sim._metadata` / `relion_src._metadata` — now a `dict`, not a pandas DataFrame

```python
# Notebook cells 12 and 22 (exploratory, not required for the pipeline):
sim._metadata
relion_src._metadata
```

These don't crash, but in 0.14.3 `_metadata` is a plain `dict` of NumPy
arrays keyed by RELION column names (e.g. `_rlnAngleRot`, `_rlnDefocusU`),
not the pandas `DataFrame` the notebook likely rendered as a table under the
older ASPIRE version. If the book reproduces these cells for illustration,
expect a dict repr rather than a DataFrame table — cosmetic only, no fix
required unless the chapter wants a DataFrame view (`pd.DataFrame(sim._metadata)`
works fine as a drop-in for display purposes).

---

## 8. Things that did *not* break (confirmed working as-is)

- `from aspire.source.simulation import Simulation` — deep submodule import
  still resolves in 0.14.3 (though `generate_synthetic.py` prefers the
  shorter `from aspire.source import Simulation`; either works).
- `Volume(vols)` — still accepts a raw 3-D `(L, L, L)` NumPy array directly
  (auto-expands to a 1-volume stack); no need to route through `Volume.load()`.
- `from aspire.utils.matrix import anorm`, `aspire.utils.types.utest_tolerance`,
  `aspire.utils.random.Random`, `IdentityFilter` — all import cleanly.
- `sim.save(starfile_filepath, batch_size=num_imgs, overwrite=True)` —
  signature unchanged; output `.mrcs` filename convention unchanged
  (`{stem}_{start}_{end-1}.mrcs`, e.g. `simulate_0_4999.mrcs` for
  `num_imgs=5000`), so the notebook's hardcoded filename in cell 23 is still
  correct **only because it matches `num_imgs=5000` exactly** — if the
  chapter parameterizes `num_imgs`, that filename must be built dynamically
  (`f"simulate_0_{num_imgs-1}.mrcs"`) rather than hardcoded.
- `Simulation(..., angles=my_angles)` — still accepts an `(n, 3)` NumPy array
  of ZYZ Euler angles directly; no need to wrap in a `Rotation` object.

---

## Smoke test: timing (n=100, img_size=130) and extrapolation to n=5000

Measured on the `cryoem-book` env (`notes/smoke_test.py`, single run, no
warmup, wall-clock via `time.perf_counter()`):

| stage | n=100 time (s) | linear extrapolation to n=5000 (s) |
|---|---:|---:|
| load_volume (fixed cost, not per-image) | 0.010 | ~0.01 (unchanged) |
| build_angles | 0.0001 | ~0.005 |
| build_filters (50 filters, fixed) | 0.0001 | ~0.0001 (unchanged) |
| build_simulation | 0.002 | ~0.002 (unchanged) |
| **projections** (`sim.projections[:n]`) | 0.367 | **~18.3** |
| **eval_filters replacement** (`sim.clean_images[:n]`) | 0.459 | **~22.9** |
| add_noise (NumPy, in-memory) | 0.022 | ~1.1 |
| **save** (`sim.save(...)`, STAR + `.mrcs` write) | 0.463 | **~23.1** |
| overwrite_mrcs (re-write with noisy stack) | 0.005 | ~0.24 |
| reload (`RelionSource` + `.images[:]`) | 0.008 | ~0.4 |
| **total** | **~1.34** | **~66 s (about 1.1 minutes)** |

Extrapolation method: the three dominant stages (`projections`,
`clean_images`, `save`) scale essentially linearly with `num_imgs` at fixed
`img_size` (per-image projection/CTF work, per-image STAR rows + per-image
`.mrcs` frames), so each is scaled by 5000/100 = 50×. Fixed-cost stages
(volume load, filter construction, Simulation object build) are left
unscaled. This is a rough estimate, not a benchmark guarantee — actual
`save` time for the full run may scale slightly worse than linear once the
~338 MB output `.mrcs` file (5000 × 130 × 130 × 4 bytes) no longer fits
comfortably in page cache, and disk I/O contention on the day of class could
add variance. Budget **roughly 1–2 minutes** of wall time for the full
`num_imgs=5000` pipeline on comparable hardware, dominated by
`projections` + `clean_images` (CTF evaluation) + `save` (disk I/O) in
roughly equal thirds.

Smoke test outputs (not committed to the book) were written to
`/tmp/claude-1000/-home-phonchi-cryo-em-general/e51677a5-94a4-486b-89b8-bf3b7f14845c/scratchpad/smoke/`
(`smoke.star`, `smoke_0_99.mrcs`, `timings.json`).
