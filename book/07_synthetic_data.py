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
# - 用已知平移把粒子移回中心，並比較校正前後的平均影像。
# - 從加入雜訊前後的變異數算出 SNR，解釋為何每張影像的清晰程度仍會不同。
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
# 從 3D 密度圖裁出一張乾淨投影，只完成了前向模型的第一步。粒子的取向與位置會改變，顯微鏡再以
# CTF 調變不同空間頻率，最後才疊上雜訊。本章保留每一層影像，方便逐步比較這些因素造成的變化。
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
# 先從薄樣品的線性模型開始，投影、CTF、平移與雜訊的效果比較容易分開觀察。每張模擬影像另有
# 已知的取向與平移，可直接拿來檢查估計誤差 {cite}`singer2020`。

# %% [markdown]
# ## 設定影像數量與輸出資料夾
#
# 第一次執行時，先把 `num_imgs_default` 改成 100，跑完所有圖形並看看輸出檔案。熟悉流程後再使用
# 預設的 5,000 張。`output_dir_default` 指定輸出資料夾，STAR、MRCS 與 `ground_truth.jsonl` 都會寫到該處。

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
# ## 讀取 3D density map 與像素大小
#
# CTF 的空間頻率以 Å⁻¹ 計算，因此像素大小必須正確。MRC 讀成 NumPy array 後，ASPIRE 不會從陣列取得
# 原本的 voxel size；建立 `Simulation` 時要明確傳入 `pixel_size=2.82`。

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
# 不同軟體可能採用不同的 Euler angle 順序與 active／passive rotation 定義。跨軟體交換角度時，
# 要同時確認角度順序、單位與 rotation matrix 的作用方向。
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
# 試跑 100 張時，50 個 viewing directions 各有兩張影像，所以每個方向只分到 $\psi=0^\circ$ 與
# $180^\circ$。使用預設 5,000 張時，每個方向有 100 個平面內旋轉，直方圖會變得更密。
#
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
# ## CTF：先看離焦如何移動零點
#
# 先比較 1.5 與 2.0 µm 離焦的徑向 CTF：
#
# $$
# H(s)=\sqrt{1-w^2}\sin\chi(s)-w\cos\chi(s),
# \qquad
# \chi(s)=2\pi\left(-\frac12\Delta f\lambda s^2+\frac14 C_s\lambda^3s^4\right).
# $$
#
# 不同離焦會把 CTF 零點移到不同頻率，因此多組離焦能互相補足部分缺口。這裡取振幅對比比例
# $w=0.15$；把 $s=0$ 代入公式可得 $H(0)=-0.15$。`envelope_B=0` 讓我們先專心觀察 CTF 振盪；加入
# envelope 後，高頻振幅還會逐漸衰減。徑向模型畫出的 Thon rings 是圓形，散光則會讓零點隨方向改變。

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
# ## 讓粒子離開影像中心
#
# 從 micrograph 裁切粒子時，挑選座標通常不會剛好落在粒子中心。這裡讓每張影像在 $x$、$y$ 方向
# 各自平移一個連續的隨機距離：
#
# $$
# t_x,t_y\overset{\mathrm{iid}}{\sim}U(-4,4)\ \text{pixels}.
# $$
#
# 對 130 pixels 寬、2.82 Å／pixel 的影像而言，單軸最大位移約佔寬度 3.1%，相當於 11.28 Å。
# 這個範圍足以讓未置中的平均影像看出模糊，又不至於把粒子大幅移出裁切框。`amplitudes=1.0`
# 固定所有影像的亮度，讓接下來的比較只多出平移、CTF 與雜訊。

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
# ## 把雜訊調到 SNR = 0.1
#
# - `sim.projections[...]`：純投影 $P_RV$。
# - `sim.clean_images[...]`：套用 CTF、平移與 amplitude 後，尚未加入雜訊。
# - `noisy_images`：在上一層加入白高斯雜訊。
#
# 把整疊經 CTF 處理的影像先扣除全域平均，再用
#
# $$
# \sigma_N^2=\operatorname{Var}(X_{\mathrm{CTF}})/\mathrm{SNR}
# $$
#
# 算出要加入的雜訊變異數。這裡的 SNR 是整疊影像、扣除平均後的功率比。所有影像加入相同的
# $\sigma_N$ 後，訊號較強的影像仍會比較清楚；CTF 對各頻率的傳遞也不同，所以單張影像與各頻率殼層
# 不會同時等於 SNR 0.1。

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
# ## 把粒子移回中心
#
# 先平均含平移的影像，再用已知的 $(t_x,t_y)$ 把最多取前 500 張影像移回中心後重做平均。第三張圖使用同一批
# 取向與 CTF，但一開始就把平移設為零。三張圖都不加雜訊，平移造成的模糊會更容易看見。校正後若
# MSE 降低、相關係數升高，數值與影像外觀便朝同一方向改變。

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
    ("before recentering", shifted_average),
    ("recentered with known shifts", recentered_average),
]:
    mse = float(np.mean((average - centered_average) ** 2))
    correlation = normalized_correlation(average, centered_average)
    print(f"{label}: MSE={mse:.6g}, correlation={correlation:.6f}")

fig, axes = plt.subplots(1, 3, figsize=(10, 3.2))
for axis, (image, title) in zip(
    axes,
    [
        (shifted_average, "Uncorrected average"),
        (recentered_average, "Recentered with known shifts"),
        (centered_average, "Zero-shift reference"),
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
# ## 輸出影像與已知真值
#
# MRCS 存放影像，STAR 存放分析軟體會用到的中繼資料，`ground_truth.jsonl` 則記錄模擬時使用的方向、
# CTF、位移與狀態標籤。練習估計參數時先只讀影像，完成後再打開真值計算誤差。

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
# 接著用 `get_metadata()` 讀回 STAR 檔，確認影像索引、離焦與像素大小都有成功儲存。

# %%
relion_source = RelionSource(str(star_path))
reloaded_images = relion_source.images[: min(8, num_imgs)].asnumpy()
public_metadata = relion_source.get_metadata(as_dict=True)
preview_columns = (
    "_rlnImageName",
    "_rlnDefocusU",
    "_rlnDefocusV",
    "_rlnImagePixelSize",
    "_rlnAngleRot",
    "_rlnAngleTilt",
    "_rlnAnglePsi",
    "_rlnOriginXAngst",
    "_rlnOriginYAngst",
)
print("first STAR record:")
for column in preview_columns:
    if column in public_metadata:
        value = np.asarray(public_metadata[column]).reshape(-1)[0]
        print(f"  {column}: {value}")
print(f"round-trip images: {reloaded_images.shape}")

# %% [markdown]
# ## 這次練習使用的設定
#
# | 參數 | 數值 | 改動後會先看到什麼 |
# |---|---:|---|
# | 影像大小 | $130\times130$ pixels | 方框太小時，平移後的粒子容易碰到邊界 |
# | 像素大小 | 2.82 Å／pixel | Nyquist frequency 與 CTF 的實體頻率會改變 |
# | 影像數 | 預設 5,000；試跑可用 100 | 影像少時，方向與離焦分布的隨機起伏較明顯 |
# | viewing directions | 50 個等向隨機方向 | 方向少時，球面覆蓋會變得稀疏 |
# | 平移 | 每軸 $U(-4,4)$ pixels | 範圍愈大，未校正平均影像愈模糊 |
# | 離焦 | 1.5–2.0 µm，共 50 組 | CTF 零點的位置與互補程度會改變 |
# | 振幅對比比例 | 0.15 | CTF 的低頻值與零點會改變 |
# | 全域功率 SNR | 0.1 | 數值降低時，粒子在雜訊中更難辨認 |

# %% [markdown]
# ## 把題目再變難
#
# 目前先固定一種 70S 構形與等向方向，使用徑向 CTF 和白高斯雜訊。確認平移校正、CTF 與 SNR 的圖形
# 都符合預期後，可以一次改一個條件：
#
# - 讓 $x$、$y$ 方向有不同離焦，觀察 Thon rings 如何由圓形變成橢圓。
# - 把白雜訊換成含低頻背景的 colored noise，比較 whitening 前後的功率頻譜。
# - 讓 viewing directions 集中在球面的一部分，觀察偏好取向如何改變方向覆蓋。
# - 準備兩個不同的 3D density maps 並儲存 state label，再測試分類能否把兩種構形分開。

# %% [markdown]
# ## 延伸閱讀
#
# - SPA 前向模型與計算流程可接著讀 {cite}`sigworth2016,singer2020`。
# - 2SDR 的合成實驗提供另一個 130 × 130 pixels、5,000 張粒子影像的去雜訊例子 {cite}`chung2020`。
# - [RELION classification example](https://www3.mrc-lmb.cam.ac.uk/relion/index.php?title=Classification_example)
#   示範相同 70S benchmark 與 `rlnAmplitudeContrast=0.15` 的設定。
# - 跨軟體交換角度時，可查 [3DEM conventions](https://github.com/azazellochg/3DEM-conventions) 的 Euler angle 對照。

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
#    平移只改變相位，頻譜振幅保持不變。$h_i$ 在實空間與影像卷積，$H_i$ 才在頻域與頻譜相乘。
#    離散影像使用的頻率座標還要和 pixel size、FFT 排列及平移單位一致。
#    ```
#
# 2. 為什麼等向抽取 viewing direction 時要讓 $\cos(\mathrm{tilt})$ 均勻分布？
#
#    ```{dropdown} 參考答案
#    方向以 $\mathrm{rot}\sim U(0,2\pi)$、$\cos(\mathrm{tilt})\sim U(-1,1)$ 抽取，才能讓球面上的每一小塊
#    面積具有相同機率。有限的 50 個方向仍會有疏密起伏；規則球格則由特別設計的節點與間距構成。
#    若直接讓 $\mathrm{tilt}\sim U(0,\pi)$，球面面積元素中的 $\sin(\mathrm{tilt})$ 沒有被補償，樣本會
#    過度集中在兩極。固定 seed 只負責重現同一批隨機樣本，不會改變所抽取的分布。
#    ```
#
# 3. 用已知平移把粒子移回中心後，平均影像、MSE 與相關係數應如何改變？
#
#    ```{dropdown} 參考答案
#    未置中的粒子在不同方向偏離中心，平均後的共同特徵會被攤寬。逐張套用反向平移後，平均影像應靠近
#    零平移參考：輪廓變清楚、MSE 降低，與參考影像的相關係數則升高。
#
#    MSE 會受整體強度尺度影響，相關係數則先扣除平均，較著重圖形的相對變化。兩個數值和影像外觀一起看，
#    可以避免只憑一張較銳利的圖下判斷。這裡使用已知平移；真實資料還要先從含雜訊影像估計位移。
#    ```
#
# 4. 為什麼整疊影像的 `realized_snr` 接近 0.1，單張粒子的清晰程度仍會不同？
#
#    ```{dropdown} 參考答案
#    程式計算
#    $\operatorname{Var}(X_{\mathrm{CTF}})/\operatorname{Var}(N)$；結果接近 0.1 表示整疊影像的全域變異數比
#    符合設定。訊號能量會隨方向、離焦與影像內容改變，同一個雜訊變異數因此會產生不同的單張 SNR。
#    CTF 又會依頻率增強、削弱或反轉訊號，各頻率殼層的 SNR 也會不同。
#
#    因此 0.1 描述整疊影像的整體難度。若要比較另一個 SNR 數值，還要先看它使用單張或整疊影像、
#    是否扣除平均，以及採用功率比或振幅比。
#    ```
#
# 5. `rlnAmplitudeContrast=0.15` 如何影響本章 CTF 在零頻率附近的值？
#
#    ```{dropdown} 參考答案
#    本章使用
#    $H(s)=\sqrt{1-w^2}\sin\chi(s)-w\cos\chi(s)$。零頻率時 $\chi(0)=0$，所以
#    $H(0)=-w=-0.15$。若 $w=0$，同一慣例下的 $H(0)$ 會等於零；加入 amplitude contrast 後，
#    CTF 在零頻率附近多了 cosine 分量，零點位置與低頻對比也跟著改變。
#
#    `0.15` 表示振幅對比比例 $w$；把整條 CTF 再乘以 0.15 會得到不同且錯誤的模型。CTF 正負號還會受
#    defocus、Fourier transform 與影像對比慣例影響，因此跨軟體比較時要連同完整公式核對。
#    ```
#
# 想把這個模擬接到完整 SPA 流程，可接著閱讀 {doc}`05_image_formation` 與 {doc}`06_reconstruction_validation`。
