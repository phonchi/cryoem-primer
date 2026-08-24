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
# 合成資料保留已知答案，讓方向、CTF、平移、雜訊與構形的誤差可以分開量化。

# %% [markdown]
# ## 設定影像數量與輸出資料夾
#
# 第一次執行時，可把 `num_imgs_default` 改成 100，先確認完整流程與輸出格式；正式生成時再使用 5,000。
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

print(f"n={num_imgs}, views={num_views}, output={out_dir.resolve()}")

# %% [markdown]
# ## 載入 3D density map
#
# MRC header 能記錄 voxel size，但把 NumPy array 傳給 `Volume` 時，ASPIRE 不會自動沿用該 header。
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
# 都來自球面上的等向分布。不過 50 個方向只是有限的蒙地卡羅樣本；固定 seed 會得到同一批方向，
# 不能宣稱每次都不同，也不能說它們形成完全均勻的球格。
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
# 不同離焦讓 CTF 零點錯開，但「50 組離焦」不保證任何有限頻帶都沒有共同弱點；是否互補要看實際參數、
# envelope 與取樣。本章也沒有模擬散光或額外 envelope，不能拿來驗證依賴這些效應的方法。
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
# - `sim.clean_images[...]`：套用 CTF、平移與 amplitude 後，尚未加入雜訊。
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

print(f"saved noisy images: {noisy_mrcs_path}")
print(f"saved ground truth: {ground_truth_path}")

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
# | {cite}`singer2020`, Eq. 10 | pose、投影、PSF／CTF 與加成雜訊 | 寫成實空間與 Fourier 空間的前向模型 |
# | {cite}`penczek2010`, pp. 5–8；[3DEM conventions](https://github.com/azazellochg/3DEM-conventions) | ZYZ Euler angles 與重建幾何 | 定義 `rot, tilt, psi` 與 rotation matrix |
# | {cite}`scheres2010`, pp. 273–286 | 高斯雜訊模型 | 建立可控制訊雜比的基準資料 |
# | [RELION classification example](https://www3.mrc-lmb.cam.ac.uk/relion/index.php?title=Classification_example) | 70S benchmark 的 `rlnAmplitudeContrast=0.15` | 將 amplitude contrast ratio 設為 0.15 |
#
# 50 個方向、50 組離焦與 SNR 0.1 是本章的教學設定。改變這些值會直接改變方向覆蓋、CTF 互補程度與資料難度。

# %% [markdown]
# ## 這組合成資料的適用範圍
#
# 本章使用單一 70S 核糖體狀態，將平移固定為零、亮度縮放固定為一，並關閉散光；雜訊採全域白高斯模型。
# 這組資料可用來檢查檔案讀寫、已知取向下的投影、CTF 模擬，以及在上述理想條件下比較去雜訊結果。
#
# 這組資料不含構形異質性、偏好取向、姿態估計誤差、散光或空間相關雜訊，因此不能用來評估這些問題，
# 也不能取代真實 micrograph 的驗證。

# %% [markdown]
# ## 理解檢查
#
# 1. 為什麼 `h * image` 與 `H * FFT(image)` 不能在沒有 Fourier transform 的情況下視為同一個運算？
# 2. 固定 seed 後，50 個方向的「隨機」與「可重現」如何同時成立？為什麼仍不能稱為完美均勻球格？
# 3. ASPIRE 預設 offset 的 $L/16$ 是標準差還是範圍？無界常態分布對模擬可能造成什麼邊界情況？
# 4. `realized_snr` 接近 0.1 能證明哪些事？為什麼不能推出每張影像或每個 frequency shell 都是 0.1？
# 5. `rlnAmplitudeContrast=0.15` 如何改變 CTF 在零頻率附近的值？
# 6. 跨軟體交換 ZYZ Euler angles 時，為什麼不能只複製三個角度數字？
# 7. 比較 phase flipping 與 Wiener-style correction 時，兩個輸出的目標與所需先驗資訊有何不同？
#
# 更多 SPA 成像與驗證限制見 {doc}`05_image_formation` 與 {doc}`06_reconstruction_validation`。
