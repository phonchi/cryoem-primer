"""
Smoke test: reproduce the 7-stage ASPIRE synthetic-data pipeline from
Simulation_data.ipynb under aspire==0.14.3, at reduced scale (n=100).

Purpose: verify the pipeline runs end-to-end on the installed ASPIRE
version and measure per-stage wall time so it can be extrapolated to the
full n=5000 run used in the notebook. See notes/api_notes.md for the
full list of API differences this test uncovered.

Run with:
    conda run -n cryoem-book python notes/smoke_test.py
"""

import json
import os
import time

import mrcfile
import numpy as np

# --- Fixed API imports for aspire 0.14.3 -----------------------------------
# NOTE: `from aspire.utils.coor_trans import uniform_random_angles` (used by
# the notebook) no longer exists in 0.14.3 and is dropped entirely -- it was
# unused dead code in the notebook anyway (angles are built by hand).
from aspire.operators import RadialCTFFilter
from aspire.source import Simulation           # was aspire.source.simulation
from aspire.source.relion import RelionSource
from aspire.volume import Volume
from aspire.utils.matrix import anorm  # noqa: F401  (kept for parity; unused below)

OUT_DIR = "/tmp/claude-1000/-home-phonchi-cryo-em-general/e51677a5-94a4-486b-89b8-bf3b7f14845c/scratchpad/smoke"
VOL_PATH = "/home/phonchi/cryo-em_general/dataset/70S_Conform1.mrc"
os.makedirs(OUT_DIR, exist_ok=True)

timings = {}


def tic():
    return time.perf_counter()


def toc(t0):
    return time.perf_counter() - t0


# --- Stage 1: parameters (n=100 for the smoke test; everything else as in
# the notebook) -------------------------------------------------------------
img_size = 130
num_imgs = 100          # notebook uses 5000; scaled down for the smoke test
duplicate = 100
num_maps = 1
sn_ratio = 0.1

pixel_size = 2.82
voltage = 200
defocus_min = 1.5e4
defocus_max = 2.0e4
defocus_ct = 50
Cs = 2.0
alpha = 0.1

np.random.seed(0)

# --- Stage 2: load volume ---------------------------------------------------
t0 = tic()
with mrcfile.open(VOL_PATH) as infile:
    vols = infile.data
v = Volume(vols)
timings["load_volume"] = toc(t0)
print(f"[stage2] Volume loaded: shape={v.shape}, dtype={v.dtype}  ({timings['load_volume']:.3f}s)")

# --- Stage 3: Euler angles ---------------------------------------------------
t0 = tic()
n_dirs = max(1, int(num_imgs / duplicate))
rot = np.repeat(np.random.random(n_dirs) * 2 * np.pi, duplicate)[:num_imgs]
tilt = np.repeat(np.arccos(2 * np.random.random(n_dirs) - 1), duplicate)[:num_imgs]
psi = np.tile(np.linspace(0, 2 * np.pi, duplicate, endpoint=False), n_dirs)[:num_imgs]
my_angles = np.column_stack((rot, tilt, psi))
timings["build_angles"] = toc(t0)
print(f"[stage3] angles shape={my_angles.shape}  ({timings['build_angles']:.4f}s)")

# --- Stage 4: CTF filters -----------------------------------------------------
# FIX: RadialCTFFilter no longer takes pixel_size as its first positional
# argument (0.14.3 signature is (voltage=200, defocus=15000, Cs=2.26,
# alpha=0.07, B=0) -- no pixel_size param at all). The notebook's
# `RadialCTFFilter(pixel_size, voltage, defocus=d, Cs=2.0, alpha=0.1)` would
# raise "got multiple values for argument 'defocus'" (pixel_size lands in
# `voltage`, voltage lands in `defocus` positionally, and then `defocus=d`
# collides). pixel_size now belongs on the Simulation/ImageSource instead
# (see stage 5).
t0 = tic()
filters = [
    RadialCTFFilter(voltage=voltage, defocus=d, Cs=Cs, alpha=alpha)
    for d in np.linspace(defocus_min, defocus_max, defocus_ct)
]
timings["build_filters"] = toc(t0)
print(f"[stage4] built {len(filters)} RadialCTFFilter objects  ({timings['build_filters']:.4f}s)")

# --- Stage 5: Simulation ------------------------------------------------------
# FIX: must pass pixel_size= explicitly. This does NOT raise an error if
# omitted -- aspire/source/simulation.py silently falls back to
# `pixel_size = self.vols.pixel_size or 1.0`. Since `v = Volume(vols)` above
# never set a pixel_size either, omitting it here means every CTFFilter is
# silently evaluated at pixel_size=1.0 A instead of 2.82 A -- verified this
# changes sim.clean_images by ~150% relative L2 norm (not a rounding error).
# See notes/api_notes.md item 3 for the measurement.
t0 = tic()
sim = Simulation(
    L=img_size,
    n=num_imgs,
    vols=v,
    C=num_maps,
    unique_filters=filters,
    angles=my_angles,
    offsets=0.0,
    amplitudes=1.0,
    pixel_size=pixel_size,   # NEW: required in 0.14.3, absent from notebook
)
timings["build_simulation"] = toc(t0)
print(f"[stage5] Simulation built, n={sim.n}  ({timings['build_simulation']:.4f}s)")

# --- Stage 6: projections -> CTF -> noise ------------------------------------
# FIX: `sim.projections` is now a property returning a subscriptable
# `_ImageAccessor`, not a callable. The notebook's
# `sim.projections(start=0, num=num_imgs)` raises
# "'_ImageAccessor' object is not callable". Use `sim.projections[:num_imgs]`.
t0 = tic()
imgs_clean = sim.projections[:num_imgs]
timings["projections"] = toc(t0)
print(f"[stage6a] projections: {type(imgs_clean).__name__}, shape={imgs_clean.shape}  ({timings['projections']:.3f}s)")

# FIX: `Simulation.eval_filters()` no longer exists anywhere in aspire
# 0.14.3 (removed; grep of the installed package returns zero hits). The
# equivalent replacement is the `sim.clean_images` accessor, which returns
# projections with filters/shifts/amplitudes applied but no noise -- exactly
# what eval_filters(projections) used to produce.
t0 = tic()
imgs_ctf_clean = sim.clean_images[:num_imgs]
timings["eval_filters"] = toc(t0)
print(f"[stage6b] clean_images (CTF-applied): shape={imgs_ctf_clean.shape}  ({timings['eval_filters']:.3f}s)")

imgs_ctf_clean = imgs_ctf_clean.asnumpy()
print(f"[stage6b] clean var = {imgs_ctf_clean.var():.6g}")

t0 = tic()
noise_var = imgs_ctf_clean.var() / sn_ratio
imgs_noise = imgs_ctf_clean + np.sqrt(noise_var) * np.random.randn(num_imgs, img_size, img_size)
timings["add_noise"] = toc(t0)
print(f"[stage6c] noise added, noise_var={noise_var:.6g}  ({timings['add_noise']:.4f}s)")

# --- Stage 7: save + reload ---------------------------------------------------
star_path = os.path.join(OUT_DIR, "smoke.star")
t0 = tic()
sim.save(star_path, batch_size=num_imgs, overwrite=True)
timings["save"] = toc(t0)
print(f"[stage7a] sim.save() done  ({timings['save']:.3f}s)")

mrcs_path = os.path.join(OUT_DIR, f"smoke_0_{num_imgs - 1}.mrcs")
assert os.path.exists(mrcs_path), f"expected mrcs at {mrcs_path}, contents: {os.listdir(OUT_DIR)}"

t0 = tic()
with mrcfile.open(mrcs_path, mode="r+") as mrc:
    mrc.set_data(np.array(imgs_noise, dtype=np.float32))
timings["overwrite_mrcs"] = toc(t0)
print(f"[stage7b] overwrote mrcs with noisy stack  ({timings['overwrite_mrcs']:.3f}s)")

# FIX: RelionSource._images() is an internal method requiring an explicit
# `indices` array argument; the notebook's `relion_src._images()` (no args)
# raises "missing 1 required positional argument: 'indices'". Use the public
# `.images` accessor instead: `relion_src.images[:]`.
t0 = tic()
relion_src = RelionSource(star_path)
imgs_reloaded = relion_src.images[:]
timings["reload"] = toc(t0)
print(f"[stage7c] reloaded via RelionSource: shape={imgs_reloaded.shape}  ({timings['reload']:.3f}s)")

print("\n=== timings (seconds) ===")
for k, val in timings.items():
    print(f"  {k:20s} {val:.4f}")

with open(os.path.join(OUT_DIR, "timings.json"), "w") as f:
    json.dump({"num_imgs": num_imgs, "img_size": img_size, "timings": timings}, f, indent=2)

print("\nSMOKE TEST PASSED")
