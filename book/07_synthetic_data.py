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
# - 說明有限個隨機方向如何近似球面均勻分布，並辨認有限取樣造成的疏密起伏。
# - 說清楚本章全域加性白高斯雜訊的 SNR 定義與限制。
# - 產生影像時，同步輸出 `ground_truth.jsonl`，記錄每張影像的已知真值。
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
# 合成資料保留已知真值，讓方向、CTF、平移、雜訊與構形的誤差可以分開量化。

# %% [markdown]
# ## 設定影像數量與輸出資料夾
#
# 第一次執行時，可把 `num_imgs_default` 改成 100，先確認完整流程與輸出格式；正式產生資料時再使用 5,000。
# `output_dir_default` 指定輸出資料夾，STAR、MRCS 與 `ground_truth.jsonl` 都會寫到該處。

# %%
import json
import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import mrcfile
import numpy as np

from aspire.operators import RadialCTFFilter
from aspire.image import Image
from aspire.source import Simulation
from aspire.source.relion import RelionSource
from aspire.volume import Volume

logging.getLogger("aspire").setLevel(logging.WARNING)

# %%
seed = 0
rng = np.random.default_rng(seed)

img_size = 130
num_imgs_default = 5000
output_dir_default = Path("synthetic_data_output")
num_imgs = int(os.environ.get("CRYOEM_NUM_IMAGES", str(num_imgs_default)))
if num_imgs < 1:
    raise ValueError("CRYOEM_NUM_IMAGES 必須是正整數")

num_views = min(50, num_imgs)
num_maps = 1
sn_ratio = 0.1
max_shift_px = 4.0
shift_seed = seed + 1

pixel_size = 2.82       # Å / pixel
voltage = 200           # kV
defocus_min = 1.5e4     # Å = 1.5 µm
defocus_max = 2.0e4     # Å = 2.0 µm
defocus_ct = 50
Cs = 2.0                # mm
amplitude_contrast = 0.15  # RELION classification example: rlnAmplitudeContrast=0.15
envelope_B = 0.0        # 本章不額外模擬 envelope decay

volume_path = Path("data/70S_Conform1.mrc")
out_dir = Path(os.environ.get("CRYOEM_OUTPUT_DIR", str(output_dir_default)))
out_dir.mkdir(parents=True, exist_ok=True)

print(f"n={num_imgs}, views={num_views}")

# %% [markdown]
# ## 載入 3D density map
#
# MRC header 能記錄 voxel size。把 NumPy array 傳給 `Volume` 後，這項資訊不會跟著陣列進入 ASPIRE，
# 因此稍後仍要把 `pixel_size` 明確傳給 `Simulation`。

# %%
with mrcfile.open(volume_path) as infile:
    volume_array = np.asarray(infile.data).copy()

volume = Volume(volume_array)
print(f"volume={volume_array.shape}, dtype={volume_array.dtype}")

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
# `psi` 是平面內旋轉。本章讓 $\mathrm{rot}\sim U(0,2\pi)$、$\cos(\mathrm{tilt})\sim U(-1,1)$，所以每個 viewing direction
# 都來自球面上的等向分布。50 個方向是有限的蒙地卡羅樣本，仍可看到隨機的疏密起伏；固定 seed
# 會重現同一批方向。規則球格則由設計好的節點與間距構成，性質不同。
#
# 不同軟體的 Euler angle 順序與 active／passive rotation 定義可能不同。
# [3DEM conventions](https://github.com/azazellochg/3DEM-conventions) 整理了常見軟體的慣例；跨軟體交換角度時，
# 應同時確認角度順序、單位與 rotation matrix 的作用方向。
#
# | 本章欄位 | 定義 |
# |---|---|
# | `rot, tilt, psi` | ASPIRE 使用的 ZYZ Euler angles |
# | 單位 | radians |
# | `rotation_matrix` | $\mathbf q_{volume}=R[k_x,k_y,0]^{\mathsf T}$ 的 image-plane-to-volume mapping |

# %%
def sample_isotropic_view_angles(random_generator, n_views):
    """Sample viewing directions with rot uniform and cos(tilt) uniform."""
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


def sample_uniform_offsets(random_generator, n_images, max_abs_shift_px):
    """Sample independent continuous x/y translations in pixel units."""
    if max_abs_shift_px < 0:
        raise ValueError("max_abs_shift_px must be non-negative")
    return random_generator.uniform(
        -max_abs_shift_px, max_abs_shift_px, size=(n_images, 2)
    )


rot_view, tilt_view = sample_isotropic_view_angles(rng, num_views)
angles, view_ids = expand_view_angles(rot_view, tilt_view, num_imgs)
offsets = sample_uniform_offsets(
    np.random.default_rng(shift_seed), num_imgs, max_shift_px
)

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
# 對照錯誤取法：若讓 $\mathrm{tilt}\sim U(0,\pi)$，球面面積元素中的 $\sin(\mathrm{tilt})$ 沒有被補償，樣本會偏向兩極。

# %%
n_demo = 1500
rot_demo = rng.uniform(0.0, 2 * np.pi, size=n_demo)
u_demo = rng.uniform(0.0, 1.0, size=n_demo)
tilt_naive = np.pi * u_demo
tilt_isotropic = np.arccos(2 * u_demo - 1)

fig = plt.figure(figsize=(11, 5))
for panel, (tilt_demo, title) in enumerate(
    [(tilt_naive, "tilt uniform: pole-biased"),
     (tilt_isotropic, "cos(tilt) uniform: isotropic draws")]
):
    x = np.sin(tilt_demo) * np.cos(rot_demo)
    y = np.sin(tilt_demo) * np.sin(rot_demo)
    z = np.cos(tilt_demo)
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
# 不同離焦會讓 CTF 零點錯開。各組 CTF 在某段頻率是否互補，仍取決於實際參數、envelope 與取樣。
# 本章省略散光與額外 envelope；相關方法仍需使用更完整的合成資料或真實資料評估。
# [RELION classification example](https://www3.mrc-lmb.cam.ac.uk/relion/index.php?title=Classification_example)
# 將這組 70S benchmark 的 `rlnAmplitudeContrast` 設為 0.15；本章採用相同的 amplitude contrast ratio。

# %%
def electron_wavelength(voltage_kv):
    """Relativistic electron wavelength in angstrom; voltage_kv is in kV."""
    voltage_v = voltage_kv * 1e3
    return 12.2639 / np.sqrt(voltage_v + 0.97845e-6 * voltage_v**2)


def ctf_1d(spatial_frequency, defocus_A, voltage_kv=200, cs_mm=2.0, w=0.15):
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
        ctf_1d(frequency, defocus, voltage, Cs, amplitude_contrast),
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
        alpha=amplitude_contrast,
        B=envelope_B,
    )
    for defocus in defocus_values
]
filter_indices = rng.integers(0, defocus_ct, size=num_imgs)

# %% [markdown]
# ## 組裝 `Simulation`：加入平移並固定亮度
#
# 原始 2SDR 合成實驗交代了影像大小、方向、CTF 與雜訊，沒有列出 2D 平移分布
# （{cite}`chung2020`）。本章另為每張影像抽取連續平移：
#
# $$
# t_x,t_y\overset{\mathrm{iid}}{\sim}U(-4,4)\ \text{pixels}.
# $$
#
# 對 130 像素寬、2.82 Å／pixel 的影像而言，單軸最大位移約佔寬度 3.1%，相當於 11.28 Å。
# 平移使用獨立的隨機種子，因此改動平移設定時，取向仍保持原來的抽樣結果。`amplitudes=1.0` 讓所有影像維持相同亮度。
# `pixel_size=2.82` 也明確傳入 ASPIRE，避免套用與密度圖不符的預設值。

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
    offsets=offsets,
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
# ## 三層影像與全域加性白高斯雜訊
#
# - `sim.projections[...]`：純投影 $P_RV$。
# - `sim.clean_images[...]`：套用 CTF、平移與 amplitude 後，尚未加入雜訊。
# - `noisy_images`：本章另加的白高斯雜訊。
#
# 本章把整疊經 CTF 處理的影像視為一個母體，先扣除全域平均，以
#
# $$
# \sigma_N^2=\operatorname{Var}(X_{\mathrm{CTF}})/\mathrm{SNR}
# $$
#
# 設定單一雜訊變異數。這裡的 SNR 是用整批影像、扣除平均後的變異數定義。`realized_snr` 接近 0.1，
# 只檢查整疊影像的變異數比是否符合設定；單張影像、遮罩內區域與各頻率殼層會有不同的訊雜比。
# 真實 micrograph 的背景往往還帶有空間相關與非定常成分。

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

# %% [markdown]
# ## 已知平移的置中對照
#
# 平移會讓粒子的共同特徵在平均影像中變寬。以下取前 500 張影像，另外產生一份平移為零的置中對照，
# 再用已知的 $(t_x,t_y)$ 將含平移影像移回中心。這項比較只使用無雜訊的 CTF 影像，以便單獨觀察
# 平移造成的差異。校正前後都與同一份置中對照比較；MSE 越小、相關係數越接近 1，表示置中越準確。

# %%
diagnostic_count = min(500, num_imgs)
centered_sim = Simulation(
    L=img_size,
    n=diagnostic_count,
    vols=volume,
    C=num_maps,
    states=state_indices_aspire[:diagnostic_count],
    unique_filters=ctf_filters,
    filter_indices=filter_indices[:diagnostic_count],
    angles=angles[:diagnostic_count],
    offsets=0.0,
    amplitudes=1.0,
    pixel_size=pixel_size,
    seed=seed,
)
centered_control = centered_sim.clean_images[:diagnostic_count].asnumpy()
shifted_subset = ctf_clean_images[:diagnostic_count]
recentered_subset = Image(shifted_subset, pixel_size=pixel_size).shift(
    -offsets[:diagnostic_count]
).asnumpy()


def normalized_correlation(first, second):
    """Return Pearson correlation after subtracting each array's mean."""
    first_centered = np.asarray(first, dtype=np.float64) - np.mean(first)
    second_centered = np.asarray(second, dtype=np.float64) - np.mean(second)
    denominator = np.linalg.norm(first_centered) * np.linalg.norm(second_centered)
    return float(np.sum(first_centered * second_centered) / denominator)


centered_average = np.mean(centered_control, axis=0)
shifted_average = np.mean(shifted_subset, axis=0)
recentered_average = np.mean(recentered_subset, axis=0)

for label, average in [
    ("校正前", shifted_average),
    ("用已知平移校正後", recentered_average),
]:
    mse = float(np.mean((average - centered_average) ** 2))
    correlation = normalized_correlation(average, centered_average)
    print(f"{label}: MSE={mse:.6g}, correlation={correlation:.6f}")

fig, axes = plt.subplots(1, 3, figsize=(10, 3.2))
for axis, (image, title) in zip(
    axes,
    [
        (shifted_average, "校正前平均"),
        (recentered_average, "已知平移校正後"),
        (centered_average, "零平移對照"),
    ],
):
    axis.imshow(image, cmap="gray")
    axis.set_title(title)
    axis.axis("off")
plt.tight_layout()
plt.show()

# %%
n_show = min(8, num_views)
show_indices = np.array([np.flatnonzero(view_ids == view_id)[0] for view_id in range(n_show)])
layers = [
    (clean_projections, "projection"),
    (ctf_clean_images, "+ CTF"),
    (noisy_images, f"+ 全域 AWGN（目標 SNR={sn_ratio}）"),
]
fig, axes = plt.subplots(3, n_show, figsize=(1.8 * n_show, 5.6), squeeze=False)
for row, (stack, label) in enumerate(layers):
    for column, image_index in enumerate(show_indices):
        axes[row, column].imshow(stack[image_index], cmap="gray")
        axes[row, column].set_xticks([])
        axes[row, column].set_yticks([])
    axes[row, 0].set_ylabel(label)
plt.suptitle("前向模型：投影 → CTF → 雜訊")
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
# ## 輸出影像與已知真值
#
# STAR 儲存分析軟體可讀取的中繼資料；`ground_truth.jsonl` 另存模擬時才知道的方向、CTF、位移與狀態標籤。
# 將已知真值和觀測資料分開，可避免後續分析程式無意間讀到答案。

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

euler_convention = (
    "ASPIRE ZYZ Euler angles (rot, tilt, psi), radians; "
    "passed directly to Simulation(angles=...)"
)

ground_truth_path = out_dir / "ground_truth.jsonl"
with ground_truth_path.open("w", encoding="utf-8") as stream:
    for image_index in range(num_imgs):
        filter_index = int(sim.filter_indices[image_index])
        record = {
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
                "amplitude_contrast": amplitude_contrast,
                "B": envelope_B,
            },
            "offset_px": [
                float(sim.offsets[image_index, 0]),
                float(sim.offsets[image_index, 1]),
            ],
            "state": {"aspire_index": 1, "label": "70S_Conform1"},
            "pixel_size_A": pixel_size,
        }
        stream.write(json.dumps(record, ensure_ascii=False) + "\n")

print(f"saved noisy images: {noisy_mrcs_path.name}")
print(f"saved ground truth: {ground_truth_path.name}")

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
# ## 模型與參數的參考來源
#
# | 來源 | 提供的概念 | 本章的使用方式 |
# |---|---|---|
# | {cite}`sigworth2016`, p. 58, Fig. 1 | projection → CTF → noise | 建立合成影像的基本順序 |
# | {cite}`singer2020`, Eq. 10 | 姿態、投影、PSF／CTF 與加性雜訊 | 寫成實空間與 Fourier 空間的前向模型 |
# | {cite}`penczek2010`, pp. 5–8；[3DEM conventions](https://github.com/azazellochg/3DEM-conventions) | ZYZ Euler angles 與重建幾何 | 定義 `rot, tilt, psi` 與 rotation matrix |
# | {cite}`scheres2010`, pp. 273–286 | 高斯雜訊模型 | 建立可控制訊雜比的基準資料 |
# | {cite}`chung2020` | 130 × 130 pixels、50 個方向、5,000 張影像與 50 組 CTF 的 2SDR 合成實驗 | 延續影像數量與有限方向的設定；原文未指定 2D 平移分布 |
# | [RELION classification example](https://www3.mrc-lmb.cam.ac.uk/relion/index.php?title=Classification_example) | 70S benchmark 的 `rlnAmplitudeContrast=0.15` | 將 amplitude contrast ratio 設為 0.15 |
#
# 每軸 $U(-4,4)$ pixels 的平移與 SNR 0.1 是本章的教學設定。改變方向數、離焦組數、平移範圍或
# 訊雜比，都會直接改變資料難度。

# %% [markdown]
# ## 這組合成資料的適用範圍
#
# 本章固定使用一個 70S 核糖體密度圖，為每張影像加入不超過 4 pixels 的 x、y 平移，亮度維持為一，
# CTF 採無散光的徑向模型，背景則加入全域加性白高斯雜訊。學生可從已知真值逐項核對檔案讀寫、
# 投影方向、CTF、平移校正與去雜訊結果。
#
# 構形異質性、偏好取向、姿態估計誤差、散光與空間相關雜訊均未納入本例。研究這些效應時，應另外設計
# 含有對應變因的合成資料，並以真實 micrograph 檢查方法在實際背景與成像誤差下的表現。

# %% [markdown]
# ## 理解檢查
#
# 1. 請從實空間模型推導 Fourier 空間模型，並指出平移 $\mathbf t_i$ 對頻譜造成的變化。
#
#    ```{dropdown} 參考答案
#    對
#    $X_i(\mathbf x)=h_i*[P_{R_i}V](\mathbf x-\mathbf t_i)+N_i(\mathbf x)$
#    取 Fourier transform。卷積定理把 $h_i*$ 變成 $H_i(\mathbf k)$ 的乘法；平移定理把
#    $f(\mathbf x-\mathbf t_i)$ 變成
#    $\widehat f(\mathbf k)e^{-2\pi\mathrm{i}\mathbf k\cdot\mathbf t_i}$。因此
#
#    $$
#    \widehat X_i(\mathbf k)=H_i(\mathbf k)\widehat{P_{R_i}V}(\mathbf k)
#    e^{-2\pi\mathrm{i}\mathbf k\cdot\mathbf t_i}+\widehat N_i(\mathbf k).
#    $$
#
#    平移只改變相位，頻譜振幅保持不變。常見誤解是把 $h_i$ 和 $H_i$ 當成同一個陣列直接相乘；前者在
#    實空間與影像卷積，後者才在頻域與頻譜相乘。離散影像使用的頻率座標也必須和 pixel size、FFT
#    排列及平移單位一致。
#    ```
#
# 2. 本章如何抽取方向與平移？固定隨機種子後，為何仍可稱為隨機抽樣？
#
#    ```{dropdown} 參考答案
#    方向以 $\mathrm{rot}\sim U(0,2\pi)$、$\cos(\mathrm{tilt})\sim U(-1,1)$ 抽取，才能讓球面上的每一小塊
#    面積具有相同機率。有限的 50 個方向仍會有疏密起伏；規則球格則由特別設計的節點與間距構成。
#    平面內角 `psi` 在每個方向的重複影像之間鋪在 $[0,2\pi)$，而 x、y 平移分別由
#    $U(-4,4)$ pixels 抽取。
#
#    偽隨機產生器從固定種子出發會重現同一串數值；「隨機」描述抽樣模型，「可重現」描述程式在相同
#    初始狀態下的結果，兩者可以同時成立。本章讓平移使用 `seed + 1`，使調整平移程式時不會改掉既有方向。
#    校正結果可和零平移對照比較；MSE 降低且相關係數升高，表示影像更接近正確中心。
#    ```
#
# 3. `realized_snr` 接近 0.1 能檢查什麼？它對單張影像與各頻率殼層有何限制？
#
#    ```{dropdown} 參考答案
#    程式計算
#    $\operatorname{Var}(X_{\mathrm{CTF}})/\operatorname{Var}(N)$；結果接近 0.1 表示整疊影像的全域變異數比
#    符合設定。訊號能量會隨方向、離焦與影像內容改變，同一個雜訊變異數因此會產生不同的單張 SNR。
#    CTF 又會依頻率增強、削弱或反轉訊號，各頻率殼層的 SNR 也會不同。
#
#    全域 SNR 無法代替單張、遮罩區域或頻率殼層的 SNR。比較方法時必須先確認 SNR 的估計範圍、
#    是否扣除平均，以及採用功率比或振幅比。
#    ```
#
# 4. `rlnAmplitudeContrast=0.15` 如何影響本章 CTF 在零頻率附近的值？
#
#    ```{dropdown} 參考答案
#    本章使用
#    $H(s)=\sqrt{1-w^2}\sin\chi(s)-w\cos\chi(s)$。零頻率時 $\chi(0)=0$，所以
#    $H(0)=-w=-0.15$。若 $w=0$，同一慣例下的 $H(0)$ 會等於零；加入 amplitude contrast 後，
#    CTF 在零頻率附近多了 cosine 分量，零點位置與低頻對比也跟著改變。
#
#    `0.15` 表示振幅對比比例 $w$；若把整條 CTF 再乘以 0.15，會得到另一個錯誤模型。CTF 正負號還會受
#    defocus、Fourier transform 與影像對比慣例影響；跨軟體比較時要連同完整公式核對。本章的數值取自
#    RELION 70S classification example。
#    ```
#
# 5. `rot, tilt, psi` 與 `rotation_matrix` 各記錄什麼？跨軟體交換時需要核對哪些慣例？
#
#    ```{dropdown} 參考答案
#    本章使用 ASPIRE 的 ZYZ Euler angles：`rot` 與 `tilt` 決定 viewing direction，`psi` 表示平面內旋轉，
#    三者都以 radians 儲存。`rotation_matrix` 則直接記錄同一旋轉的矩陣形式；本章的語意是
#    $\mathbf q_{\mathrm{volume}}=R[k_x,k_y,0]^{\mathsf T}$，把影像平面的 Fourier 座標映到體積座標。
#
#    跨軟體時須一起核對 Euler 順序、角度單位、active／passive 定義、矩陣左右乘、座標軸方向及影像到
#    體積或體積到影像的映射。三個數字相同，只能說欄位值相同；上述慣例一致後，才表示相同的物理取向。
#    可用單位向量或已知非對稱體積做一次投影測試，直接檢查轉換結果。
#    ```
#
# 更多 SPA 成像與驗證限制見 {doc}`05_image_formation` 與 {doc}`06_reconstruction_validation`。
