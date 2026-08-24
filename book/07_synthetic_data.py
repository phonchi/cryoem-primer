# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: cryoem-book
# ---

# %% [markdown]
# # 合成資料：把 SPA 前向模型寫成程式
#
# ```{admonition} 讀完本章，你應該能
# :class: important
# - 在實空間與 Fourier 空間正確寫出 projection、CTF、平移與雜訊。
# - 說明有限個隨機方向只是球面均勻分布的樣本，不是「完美均勻覆蓋」。
# - 說清楚本章 global AWGN 的 SNR 定義與限制。
# - 產生影像時，同步保存 `ground_truth.jsonl` 與 `run_manifest.json`。
# ```
#
# ```{figure} images/pptx/s22_1.png
# :width: 40%
# :name: fig-70s-map
#
# 本章使用 EMD-1056 的 70S ribosome map，產生帶 CTF 與雜訊的模擬粒子影像。
# ```

# %% [markdown]
# (forward-model)=
# ## 前向模型有兩種等價寫法
#
# 第 $i$ 張影像在實空間寫成
#
# $$
# X_i(\mathbf x)=h_i * [P_{R_i}V](\mathbf x-\mathbf t_i)+N_i(\mathbf x),
# $$
#
# 在 Fourier 空間則是
#
# $$
# \widehat X_i(\mathbf k)=H_i(\mathbf k)\widehat{P_{R_i}V}(\mathbf k)
# e^{-2\pi\mathrm{i}\mathbf k\cdot\mathbf t_i}+\widehat N_i(\mathbf k).
# $$
#
# $V$ 是 3D density，$P_{R_i}$ 是取向 $R_i$ 下的投影，$h_i$ 是 point-spread function，
# $H_i=\mathcal F\{h_i\}$ 才是 CTF。也就是說：**實空間與 PSF 卷積，Fourier 空間與 CTF 相乘**。
# 這個線性模型是薄樣品 SPA 的教學近似；它沒有模擬多重散射、空間非平穩背景、逐 frame 輻射損傷，
# 也沒有模擬 beam-induced motion（{cite}`singer2020`, Eq. 10；詳見 {doc}`05_image_formation`）。
#
# 合成資料的用途不是做出「像真的」圖片，而是保存已知答案，讓方向、CTF、平移、雜訊與構形的誤差可以分開量化。

# %% [markdown]
# ## 參數與輸出位置
#
# 預設產生 5,000 張。網站快速建置可設 `CRYOEM_NUM_IMAGES=100`；大型 MRCS 與 provenance 檔
# 寫到 `CRYOEM_OUTPUT_DIR`，未設定時放在 Jupyter Book 的 `_build/synthetic_data`，不寫進原始資料目錄。

# %%
import hashlib
import importlib.metadata
import json
import logging
import os
import platform
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import mrcfile
import numpy as np

from aspire.operators import RadialCTFFilter
from aspire.source import Simulation
from aspire.source.relion import RelionSource
from aspire.volume import Volume

logging.getLogger("aspire").setLevel(logging.WARNING)

# %%
schema_version = "1.0.0"
seed = 0
rng = np.random.default_rng(seed)

img_size = 130
num_imgs = int(os.environ.get("CRYOEM_NUM_IMAGES", "5000"))
if num_imgs < 1:
    raise ValueError("CRYOEM_NUM_IMAGES 必須是正整數")

num_views = min(50, num_imgs)
num_maps = 1
sn_ratio = 0.1

pixel_size = 2.82       # Å / pixel
voltage = 200           # kV
defocus_min = 1.5e4     # Å = 1.5 µm
defocus_max = 2.0e4     # Å = 2.0 µm
defocus_ct = 50
Cs = 2.0                # mm
alpha = 0.1             # amplitude contrast ratio
envelope_B = 0.0        # 本章不額外模擬 envelope decay

volume_path = Path("data/70S_Conform1.mrc")
out_dir = Path(os.environ.get("CRYOEM_OUTPUT_DIR", "_build/synthetic_data"))
out_dir.mkdir(parents=True, exist_ok=True)

print(f"n={num_imgs}, views={num_views}, output={out_dir.resolve()}")

# %% [markdown]
# ## 載入 3D density 並記錄來源雜湊
#
# MRC header 能記錄 voxel size，但把 NumPy array 傳給 `Volume` 時，ASPIRE 不會自動沿用該 header。
# 因此稍後仍要把 `pixel_size` 明確傳給 `Simulation`。SHA-256 用來確認本次 ground truth 對應哪一個輸入檔。

# %%
def sha256_file(path, chunk_size=1024 * 1024):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


volume_sha256 = sha256_file(volume_path)
with mrcfile.open(volume_path) as infile:
    volume_array = np.asarray(infile.data).copy()

volume = Volume(volume_array)
print(f"volume={volume_array.shape}, dtype={volume_array.dtype}, sha256={volume_sha256[:12]}…")

# %%
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for axis in range(3):
    axes[axis].contourf(np.sum(volume_array, axis=axis))
    axes[axis].set_title(f"sum along axis {axis}")
    axes[axis].set_aspect("equal")
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 取向：從等向分布抽出有限樣本
#
# ASPIRE 接受 ZYZ Euler angles `(rot, tilt, psi)`，單位是 radians。`rot` 與 `tilt` 決定 viewing direction，
# `psi` 是平面內旋轉。本章讓 $\alpha\sim U(0,2\pi)$、$\cos\beta\sim U(-1,1)$，所以每個 viewing direction
# 都來自球面上的等向分布。不過 50 個方向只是有限的蒙地卡羅樣本；固定 seed 會得到同一批方向，
# 不能宣稱每次都不同，也不能說它們形成完全均勻的球格。

# %%
def sample_isotropic_view_angles(random_generator, n_views):
    """Sample viewing directions with alpha uniform and cos(beta) uniform."""
    rot_view = random_generator.uniform(0.0, 2 * np.pi, size=n_views)
    cos_tilt = random_generator.uniform(-1.0, 1.0, size=n_views)
    tilt_view = np.arccos(cos_tilt)
    return rot_view, tilt_view


def expand_view_angles(rot_view, tilt_view, n_images):
    """Repeat finite views and spread each view's in-plane angles over [0, 2pi)."""
    n_views = len(rot_view)
    counts = np.full(n_views, n_images // n_views, dtype=int)
    counts[: n_images % n_views] += 1
    view_ids = np.repeat(np.arange(n_views), counts)
    rot = np.repeat(rot_view, counts)
    tilt = np.repeat(tilt_view, counts)
    psi = np.concatenate(
        [np.linspace(0.0, 2 * np.pi, count, endpoint=False) for count in counts]
    )
    return np.column_stack((rot, tilt, psi)), view_ids


rot_view, tilt_view = sample_isotropic_view_angles(rng, num_views)
angles, view_ids = expand_view_angles(rot_view, tilt_view, num_imgs)

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].scatter(np.degrees(rot_view), np.cos(tilt_view), s=25)
axes[0].set_xlabel("rot (degree)")
axes[0].set_ylabel("cos(tilt)")
axes[0].set_title(f"{num_views} isotropic random draws (finite sample)")
axes[1].hist(np.degrees(angles[:, 2]), bins=min(50, num_imgs))
axes[1].set_xlabel("psi (degree)")
axes[1].set_title("In-plane rotations")
plt.tight_layout()
plt.show()

# %% [markdown]
# 對照錯誤取法：若讓 $\beta\sim U(0,\pi)$，球面面積元素中的 $\sin\beta$ 沒有被補償，樣本會偏向兩極。

# %%
n_demo = 1500
alpha_demo = rng.uniform(0.0, 2 * np.pi, size=n_demo)
u_demo = rng.uniform(0.0, 1.0, size=n_demo)
beta_naive = np.pi * u_demo
beta_isotropic = np.arccos(2 * u_demo - 1)

fig = plt.figure(figsize=(11, 5))
for panel, (beta_demo, title) in enumerate(
    [(beta_naive, "tilt uniform: pole-biased"),
     (beta_isotropic, "cos(tilt) uniform: isotropic draws")]
):
    x = np.sin(beta_demo) * np.cos(alpha_demo)
    y = np.sin(beta_demo) * np.sin(alpha_demo)
    z = np.cos(beta_demo)
    ax = fig.add_subplot(1, 2, panel + 1, projection="3d")
    ax.scatter(x, y, z, s=3, alpha=0.5)
    ax.set_title(title)
    ax.set_box_aspect((1, 1, 1))
plt.tight_layout()
plt.show()

# %% [markdown]
# (ctf)=
# ## CTF：Fourier 空間的轉移函數
#
# 本章採徑向對稱、無散光的教學模型：
#
# $$
# H(s)=\sqrt{1-w^2}\sin\chi(s)-w\cos\chi(s),
# \qquad
# \chi(s)=2\pi\left(-\frac12\Delta f\lambda s^2+\frac14 C_s\lambda^3s^4\right).
# $$
#
# 不同離焦讓 CTF 零點錯開，但「50 組離焦」不保證任何有限頻帶都沒有共同弱點；是否互補要看實際參數、
# envelope 與取樣。本章也沒有模擬散光或額外 envelope，不能拿來驗證依賴這些效應的方法。

# %%
def electron_wavelength(voltage_kv):
    """Relativistic electron wavelength in angstrom; voltage_kv is in kV."""
    voltage_v = voltage_kv * 1e3
    return 12.2639 / np.sqrt(voltage_v + 0.97845e-6 * voltage_v**2)


def ctf_1d(spatial_frequency, defocus_A, voltage_kv=200, cs_mm=2.0, w=0.1):
    """Radial CTF using the sign convention used in this chapter."""
    wavelength = electron_wavelength(voltage_kv)
    cs_A = cs_mm * 1e7
    phase = 2 * np.pi * (
        -0.5 * defocus_A * wavelength * spatial_frequency**2
        + 0.25 * cs_A * wavelength**3 * spatial_frequency**4
    )
    return np.sqrt(1 - w**2) * np.sin(phase) - w * np.cos(phase)


nyquist = 1 / (2 * pixel_size)
frequency = np.linspace(0, nyquist, 800)
plt.figure(figsize=(9, 4))
for defocus, style in [(defocus_min, "-"), (defocus_max, "--")]:
    plt.plot(
        frequency,
        ctf_1d(frequency, defocus, voltage, Cs, alpha),
        style,
        label=f"defocus={defocus / 1e4:.1f} µm",
    )
plt.axhline(0, color="gray", linewidth=0.8)
plt.xlabel("spatial frequency (1/Å)")
plt.ylabel("CTF")
plt.title(f"Radial CTF at {voltage} kV; Nyquist={nyquist:.3f} 1/Å")
plt.legend()
plt.show()

# %%
defocus_values = np.linspace(defocus_min, defocus_max, defocus_ct)
ctf_filters = [
    RadialCTFFilter(
        voltage=voltage,
        defocus=defocus,
        Cs=Cs,
        alpha=alpha,
        B=envelope_B,
    )
    for defocus in defocus_values
]
filter_indices = rng.integers(0, defocus_ct, size=num_imgs)

# %% [markdown]
# ## 組裝 `Simulation`：明確關掉平移與亮度變化
#
# ASPIRE 0.14.3 若省略 `offsets`，會從**無界的常態分布**抽樣，每一軸的標準差是 $L/16$；
# $L/16$ 不是硬性範圍。這裡用 `offsets=0.0` 關閉平移，並用 `amplitudes=1.0` 關閉亮度縮放。
# 若省略 `pixel_size`，此處的 `Volume` 沒有 pixel-size metadata，ASPIRE 會靜默採用 1.0 Å；
# 因此必須明確傳入 2.82 Å。

# %%
state_indices_aspire = np.ones(num_imgs, dtype=int)  # ASPIRE 的 volume state 是 1-based
sim = Simulation(
    L=img_size,
    n=num_imgs,
    vols=volume,
    C=num_maps,
    states=state_indices_aspire,
    unique_filters=ctf_filters,
    filter_indices=filter_indices,
    angles=angles,
    offsets=0.0,
    amplitudes=1.0,
    pixel_size=pixel_size,
    seed=seed,
)
print(f"Simulation ready: n={sim.n}, L={sim.L}, pixel_size={sim.pixel_size} Å")

assigned_defocus = defocus_values[sim.filter_indices]
plt.figure(figsize=(8, 3.5))
plt.hist(assigned_defocus / 1e4, bins=defocus_ct)
plt.xlabel("defocus (µm)")
plt.ylabel("count")
plt.title("Random CTF-group assignment")
plt.show()

# %% [markdown]
# (noise)=
# ## 三層影像與 global AWGN
#
# - `sim.projections[...]`：純投影 $P_RV$。
# - `sim.clean_images[...]`：套用 CTF、平移與 amplitude 後，尚未加噪。
# - `noisy_images`：本章另加的白高斯雜訊。
#
# 本章把整疊 CTF-filtered images 視為一個母體，先扣除全域平均，以
#
# $$
# \sigma_N^2=\operatorname{Var}(X_{\mathrm{CTF}})/\mathrm{SNR}
# $$
#
# 設定單一雜訊變異數。這是 **global, mean-centered, variance SNR**。它不保證每張影像、每個 mask 或每個
# frequency shell 都是 SNR 0.1；真實 micrograph 的背景也不一定是獨立、同分布、白高斯雜訊。

# %%
def add_global_awgn(clean_images, snr, random_generator):
    """Add zero-mean global AWGN using Var(clean)/Var(noise)=snr."""
    if snr <= 0:
        raise ValueError("snr must be positive")
    signal_variance = float(np.var(clean_images, dtype=np.float64))
    noise_sigma = np.sqrt(signal_variance / snr)
    noise = random_generator.standard_normal(clean_images.shape, dtype=np.float32)
    noisy = np.asarray(clean_images, dtype=np.float32) + noise_sigma * noise
    return noisy.astype(np.float32), noise_sigma, noise


clean_projections = sim.projections[:num_imgs].asnumpy().astype(np.float32)
ctf_clean_images = sim.clean_images[:num_imgs].asnumpy().astype(np.float32)
noisy_images, noise_sigma, standardized_noise = add_global_awgn(
    ctf_clean_images, sn_ratio, rng
)
realized_noise = noise_sigma * standardized_noise
achieved_snr = float(
    np.var(ctf_clean_images, dtype=np.float64)
    / np.var(realized_noise, dtype=np.float64)
)
print(f"target global SNR={sn_ratio}; realized global SNR={achieved_snr:.4f}")

# %%
n_show = min(8, num_views)
show_indices = np.array([np.flatnonzero(view_ids == view_id)[0] for view_id in range(n_show)])
layers = [
    (clean_projections, "projection"),
    (ctf_clean_images, "+ CTF"),
    (noisy_images, f"+ global AWGN (target SNR={sn_ratio})"),
]
fig, axes = plt.subplots(3, n_show, figsize=(1.8 * n_show, 5.6), squeeze=False)
for row, (stack, label) in enumerate(layers):
    for column, image_index in enumerate(show_indices):
        axes[row, column].imshow(stack[image_index], cmap="gray")
        axes[row, column].set_xticks([])
        axes[row, column].set_yticks([])
    axes[row, 0].set_ylabel(label)
plt.suptitle("Forward model: projection → CTF → noise")
plt.tight_layout()
plt.show()

# %%
same_view_indices = np.flatnonzero(view_ids == 0)[:8]
fig, axes = plt.subplots(1, len(same_view_indices), figsize=(1.8 * len(same_view_indices), 2), squeeze=False)
for column, image_index in enumerate(same_view_indices):
    axes[0, column].imshow(clean_projections[image_index], cmap="gray")
    axes[0, column].set_title(f"ψ={np.degrees(angles[image_index, 2]):.1f}°", fontsize=8)
    axes[0, column].set_xticks([])
    axes[0, column].set_yticks([])
plt.suptitle("One viewing direction, different in-plane rotations")
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 輸出影像、ground truth 與 run manifest
#
# STAR 保存儀器與處理軟體可使用的 metadata；`ground_truth.jsonl` 另存模擬時才知道的答案，
# 避免評估程式無意間讀到 label。`run_manifest.json` 記錄全域參數、版本、SNR 定義與輸出路徑。

# %%
star_path = out_dir / "simulate.star"
sim.save(str(star_path), batch_size=num_imgs, overwrite=True)
noisy_mrcs_path = out_dir / f"simulate_0_{num_imgs - 1}.mrcs"
with mrcfile.open(noisy_mrcs_path, mode="r+") as mrc:
    mrc.set_data(noisy_images)
    mrc.voxel_size = pixel_size

projection_path = out_dir / "clean_projection.mrcs"
with mrcfile.new(projection_path, overwrite=True) as mrc:
    mrc.set_data(clean_projections)
    mrc.voxel_size = pixel_size

ctf_clean_path = out_dir / "clean_ctf.mrcs"
with mrcfile.new(ctf_clean_path, overwrite=True) as mrc:
    mrc.set_data(ctf_clean_images)
    mrc.voxel_size = pixel_size


def package_version(distribution_name):
    try:
        return importlib.metadata.version(distribution_name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


software_versions = {
    "python": platform.python_version(),
    "aspire": package_version("aspire"),
    "numpy": np.__version__,
    "mrcfile": package_version("mrcfile"),
    "matplotlib": package_version("matplotlib"),
    "platform": platform.platform(),
}
euler_convention = (
    "ASPIRE ZYZ Euler angles (rot, tilt, psi), radians; "
    "passed directly to Simulation(angles=...)"
)

ground_truth_path = out_dir / "ground_truth.jsonl"
with ground_truth_path.open("w", encoding="utf-8") as stream:
    for image_index in range(num_imgs):
        filter_index = int(sim.filter_indices[image_index])
        record = {
            "schema_version": schema_version,
            "seed": seed,
            "image_index": image_index,
            "view_id": int(view_ids[image_index]),
            "euler_convention": euler_convention,
            "euler_angles_rad": {
                "rot": float(angles[image_index, 0]),
                "tilt": float(angles[image_index, 1]),
                "psi": float(angles[image_index, 2]),
            },
            "rotation_matrix": np.asarray(sim.rotations[image_index]).tolist(),
            "rotation_semantics": (
                "ASPIRE image-plane-to-volume Fourier coordinate mapping: "
                "q_volume = R @ [k_x, k_y, 0]; passed to Volume.project"
            ),
            "filter_index": filter_index,
            "ctf": {
                "model": "ASPIRE RadialCTFFilter; no astigmatism",
                "voltage_kv": voltage,
                "defocus_u_A": float(defocus_values[filter_index]),
                "defocus_v_A": float(defocus_values[filter_index]),
                "defocus_angle_rad": 0.0,
                "Cs_mm": Cs,
                "amplitude_contrast": alpha,
                "B": envelope_B,
            },
            "offset_px": [
                float(sim.offsets[image_index, 0]),
                float(sim.offsets[image_index, 1]),
            ],
            "state": {"aspire_index": 1, "label": "70S_Conform1"},
            "pixel_size_A": pixel_size,
            "volume_sha256": volume_sha256,
        }
        stream.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")

manifest = {
    "schema_version": schema_version,
    "created_by": "book/07_synthetic_data.py",
    "command_environment": {
        "CRYOEM_NUM_IMAGES": os.environ.get("CRYOEM_NUM_IMAGES"),
        "CRYOEM_OUTPUT_DIR": os.environ.get("CRYOEM_OUTPUT_DIR"),
    },
    "seed": seed,
    "num_images": num_imgs,
    "image_size_px": img_size,
    "num_views": num_views,
    "euler_convention": euler_convention,
    "volume": {
        "path": str(volume_path),
        "sha256": volume_sha256,
        "state_label": "70S_Conform1",
    },
    "pixel_size_A": pixel_size,
    "ctf": {
        "filter_count": defocus_ct,
        "defocus_min_A": defocus_min,
        "defocus_max_A": defocus_max,
        "voltage_kv": voltage,
        "Cs_mm": Cs,
        "amplitude_contrast": alpha,
        "B": envelope_B,
        "astigmatism": False,
    },
    "offset_model": {"enabled": False, "offset_px": [0.0, 0.0]},
    "amplitude_model": {"enabled": False, "amplitude": 1.0},
    "noise_model": {
        "family": "global additive white Gaussian noise",
        "definition": "Var(mean-centered CTF-clean stack) / Var(noise)",
        "target_snr": sn_ratio,
        "realized_snr": achieved_snr,
        "noise_sigma": float(noise_sigma),
        "limitations": [
            "SNR is global, not per-image or per-shell",
            "noise is iid and spatially stationary",
            "clean target is CTF-filtered projection",
        ],
    },
    "outputs": {
        "star": str(star_path.resolve()),
        "noisy_mrcs": str(noisy_mrcs_path.resolve()),
        "clean_projection_mrcs": str(projection_path.resolve()),
        "clean_ctf_mrcs": str(ctf_clean_path.resolve()),
        "ground_truth_jsonl": str(ground_truth_path.resolve()),
    },
    "software_versions": software_versions,
}
manifest_path = out_dir / "run_manifest.json"
with manifest_path.open("w", encoding="utf-8") as stream:
    json.dump(manifest, stream, ensure_ascii=False, indent=2, sort_keys=True)
    stream.write("\n")

print(f"saved noisy images: {noisy_mrcs_path}")
print(f"saved ground truth: {ground_truth_path}")
print(f"saved run manifest: {manifest_path}")

# %% [markdown]
# 用公開的 `get_metadata()` 讀取 STAR metadata；不要依賴 `_metadata` 這類私有屬性。

# %%
relion_source = RelionSource(str(star_path))
reloaded_images = relion_source.images[: min(8, num_imgs)].asnumpy()
public_metadata = relion_source.get_metadata(as_dict=True)
print("metadata columns:")
print(sorted(public_metadata))
print(f"round-trip images: {reloaded_images.shape}")

# %% [markdown]
# ## 文獻如何約束模擬協定
#
# | 來源 | 可支持的模擬設計 | 不能從該來源延伸出的主張 |
# |---|---|---|
# | {cite}`sigworth2016`, p. 58, Fig. 1 | projection → CTF → noise 的教學骨架 | 不指定本章的 50 個方向、50 組離焦或 SNR 0.1 |
# | {cite}`singer2020`, Eq. 10 | pose、投影、PSF／CTF 與加成雜訊的線性模型 | 不保證真實雜訊是 global iid AWGN |
# | {cite}`penczek2010`, pp. 5–8 | ZYZ Euler 慣例、方向覆蓋與重建幾何 | 不表示有限隨機方向能達到均勻球格的確定性覆蓋 |
# | {cite}`scheres2010`, pp. 273–286 | 以高斯雜訊建立 ML 模型及其 hidden variables | 作者也明確提醒加成、獨立、白高斯等假設可能不符合真實資料 |
#
# 這張表只連結「來源實際支持的 claim」。本章數值是可重現的教學設定，不冒充任何一篇論文的原始 protocol。

# %% [markdown]
# ## 這份資料能驗證什麼
#
# 本章只有一個 70S 構形，而且平移、亮度變化、散光與 colored noise 都被關閉。因此它適合測試 I/O、
# CTF-aware forward model、已知 view labels 的分群或去雜訊基線；它不能單獨支持方法已能處理真實資料、
# 異質性、preferred orientation 或 alignment error 的主張。
#
# 教材目前沒有附 `70S_Conform2.mrc`，也沒有預先定義的分群輸出，所以不會假裝提供可直接執行的 purity 練習。
# 若要加入第二構形，必須先保存兩個 volume hashes 與每張影像的 state label，再明確定義如何把 cluster labels
# 對應到 state labels；purity 只衡量分群與已知狀態的一致性，不衡量姿態或重建是否正確。

# %% [markdown]
# ## 理解檢查
#
# 1. 為什麼 `h * image` 與 `H * FFT(image)` 不能在沒有 Fourier transform 的情況下視為同一個運算？
# 2. 固定 seed 後，50 個方向的「隨機」與「可重現」如何同時成立？為什麼仍不能稱為完美均勻球格？
# 3. ASPIRE 預設 offset 的 $L/16$ 是標準差還是範圍？無界常態分布對模擬可能造成什麼邊界情況？
# 4. `realized_snr` 接近 0.1 能證明哪些事？為什麼不能推出每張影像或每個 frequency shell 都是 0.1？
# 5. 若要設計兩構形 purity 實驗，除了第二個 volume，還必須新增哪些 ground-truth 欄位與評估約定？
# 6. 比較 phase flipping 與 Wiener-style correction 時，兩個輸出的目標與所需先驗資訊有何不同？
#
# 更多 SPA 成像與驗證限制見 {doc}`05_image_formation` 與 {doc}`06_reconstruction_validation`。
