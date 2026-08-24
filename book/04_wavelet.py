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
# # 小波與多尺度分析
#
# 全域 Fourier spectrum 告訴我們整段訊號的頻率組成，但不保留事件發生的位置。小波把分析函數限制在局部，再以不同尺度掃描訊號。本章從一般影像的 Gaussian／Laplacian pyramids 出發，分清 short-time Fourier transform（STFT）、continuous wavelet transform（CWT）與 decimated discrete wavelet transform（DWT）；它們都能描述局部尺度，輸出與用途卻不同。多尺度表示的影像處理背景可參考 {cite}`szeliski2022,forsyth2012`。
#
# ```{admonition} 學習目標
# :class: important
#
# - 說明有限觀測下 Fourier spectrum 的解析度與 spectral leakage。
# - 比較 STFT 固定視窗與 CWT 可變尺度的解析度。
# - 寫出 CWT 與 dyadic wavelet family，說明位移網格如何隨尺度改變。
# - 將 DWT 解釋成 analysis filter bank 加 downsampling。
# - 驗證 pyramid 與 DWT perfect reconstruction，並比較 boundary modes。
# - 比較 hard、soft、universal 與 BayesShrink 去雜訊，說明它們依賴的假設。
# ```

# %%
import matplotlib.pyplot as plt
import numpy as np
import pywt
import skimage as ski
from scipy import ndimage as ndi
from scipy import signal
from skimage.metrics import peak_signal_noise_ratio, structural_similarity
from skimage.restoration import denoise_wavelet, estimate_sigma

plt.rcParams["image.cmap"] = "gray"
plt.rcParams["figure.figsize"] = (6, 5)


def imshow_all(*images, titles=None, size=4, **kwargs):
    if titles is None:
        titles = [""] * len(images)
    fig, axes = plt.subplots(1, len(images), figsize=(size * len(images), size))
    axes = np.atleast_1d(axes)
    for ax, image, title in zip(axes, images, titles):
        ax.imshow(image, **kwargs)
        ax.set_title(title)
        ax.axis("off")
    fig.tight_layout()
    return fig, axes


# %% [markdown]
# ## 全域頻譜看得到頻率，看不到位置
#
# 有限長度 DFT 的 frequency-bin spacing 由觀測時間決定；非整數週期還會產生 spectral leakage。因此 Fourier spectrum 描述整段訊號在各個頻率 bins 的能量，解析精度仍受觀測長度與 window 影響。
#
# 下例中訊號 A 讓四個頻率全程同時存在，訊號 B 則讓它們依序出現。兩者都在相同頻率附近有能量，但 B 的每個成分只持續四分之一秒，峰值較低，分段邊界也造成明顯 leakage；兩張 magnitude spectra 的形狀與幅度因此有清楚差異。

# %%
sample_rate = 200
t = np.arange(sample_rate) / sample_rate
component_frequencies = [4, 30, 60, 90]

signal_constant = sum(np.sin(2 * np.pi * frequency * t)
                      for frequency in component_frequencies)
signal_sequential = np.zeros_like(t)
segment_length = len(t) // len(component_frequencies)
for index, frequency in enumerate(component_frequencies):
    segment = slice(index * segment_length, (index + 1) * segment_length)
    signal_sequential[segment] = np.sin(2 * np.pi * frequency * t[segment])

magnitude_constant = np.abs(np.fft.rfft(signal_constant))
magnitude_sequential = np.abs(np.fft.rfft(signal_sequential))
frequency_axis = np.fft.rfftfreq(t.size, d=1 / sample_rate)

assert magnitude_constant.max() > 3 * magnitude_sequential.max()

fig, axes = plt.subplots(2, 2, figsize=(11, 6))
axes[0, 0].plot(t, signal_constant)
axes[0, 0].set_title("A: four frequencies throughout")
axes[1, 0].plot(t, signal_sequential)
axes[1, 0].set_title("B: one frequency per interval")
axes[0, 1].plot(frequency_axis, magnitude_constant)
axes[0, 1].set_title("global spectrum of A")
axes[1, 1].plot(frequency_axis, magnitude_sequential)
axes[1, 1].set_title("global spectrum of B: lower peaks + leakage")
for ax in axes[:, 0]:
    ax.set_xlabel("time [s]")
for ax in axes[:, 1]:
    ax.set_xlabel("frequency [Hz]")
fig.tight_layout()


# %% [markdown]
# ## STFT 與 CWT 採用不同的解析度策略
#
# Short-time Fourier transform（STFT）用固定長度 window 切出局部頻譜，因此時間／頻率解析度由同一個 window 決定。短 window 能較精準定位變化時間，但 frequency bins 較粗；長 window 能分開較接近的頻率，卻會把短暫事件拉寬。下例使用 SciPy `ShortTimeFFT`；window、overlap、padding 與 scaling 都是轉換定義的一部分，解讀時必須一起報告。完整參數見 [SciPy `ShortTimeFFT` 文件](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.ShortTimeFFT.html)。

# %%
stft_short = signal.ShortTimeFFT.from_window(
    "hann", fs=sample_rate, nperseg=24, noverlap=18, scale_to="magnitude"
)
stft_long = signal.ShortTimeFFT.from_window(
    "hann", fs=sample_rate, nperseg=80, noverlap=60, scale_to="magnitude"
)
short_coefficients = stft_short.stft(signal_sequential)
long_coefficients = stft_long.stft(signal_sequential)
short_times = stft_short.t(signal_sequential.size)
long_times = stft_long.t(signal_sequential.size)

assert short_coefficients.shape == (stft_short.f.size, short_times.size)
assert long_coefficients.shape == (stft_long.f.size, long_times.size)
assert stft_long.delta_f < stft_short.delta_f

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for ax, transform, coefficients, times, title in (
    (axes[0], stft_short, short_coefficients, short_times, "short window"),
    (axes[1], stft_long, long_coefficients, long_times, "long window"),
):
    ax.pcolormesh(times, transform.f, np.abs(coefficients), shading="auto")
    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(0, sample_rate / 2)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("frequency [Hz]")
    ax.set_title(title)
fig.tight_layout()

# %% [markdown]
# CWT 改用可縮放的母小波 $\psi$：
#
# $$W(a,b)=\frac{1}{\sqrt{|a|}}\int h(t)\psi^*\!\left(\frac{t-b}{a}\right)dt.$$
#
# 大尺度通常對應較低 pseudo-frequency 與較寬的時間支撐；小尺度對應較高 pseudo-frequency 與較窄的時間支撐。CWT 對每個尺度保留密集的位移取樣，結果可畫成尺度圖（scalogram）。母小波與取樣週期共同決定尺度和 Hz 的換算，尺度的倒數只提供粗略直覺。PyWavelets 的 `cwt` 會依 `sampling_period` 和母小波的 central frequency 回傳 pseudo-frequency，參數意義可查閱 [CWT 官方文件](https://pywavelets.readthedocs.io/en/stable/ref/cwt.html)。

# %%
scales = np.arange(1, 65)
cwt_coefficients, cwt_frequencies = pywt.cwt(
    signal_sequential, scales, "morl", sampling_period=1 / sample_rate
)
assert cwt_coefficients.shape == (len(scales), len(signal_sequential))
assert np.all(np.diff(cwt_frequencies) < 0)

fig, ax = plt.subplots(figsize=(10, 4))
extent = [t[0], t[-1], cwt_frequencies[-1], cwt_frequencies[0]]
ax.imshow(np.abs(cwt_coefficients), aspect="auto", extent=extent, origin="upper")
ax.set_xlabel("time [s]")
ax.set_ylabel("pseudo-frequency [Hz]")
ax.set_title("CWT scalogram")
fig.tight_layout()

# %% [markdown]
# 訊號邊界附近放不下完整小波，程式會藉由 padding 或 extension 處理。尺度越大，受邊界影響的時間範圍越寬；解讀 cone-of-influence 區域時要比中央係數保守。若事件只出現在訊號起點，宜取得更長的觀測，並比較不同 extension 設定，再判斷它是否為穩定的物理訊號。


# %% [markdown]
# ## Dyadic family 與 DWT
#
# 一個常見的 dyadic wavelet family 可寫成
#
# $$\psi_{j,k}(t)=2^{-j/2}\psi(2^{-j}t-k).$$
#
# 以 CWT 參數理解，就是 $a=2^j$、$b=k2^j$：尺度變大時，位移網格也跟著變粗。正交／雙正交 DWT 再利用特定 scaling function 與 wavelet，形成可逆、非冗餘或低冗餘的離散表示。
#
# 實作上，一階 DWT 將訊號分別通過 analysis low-pass 與 high-pass filters，再 downsample by 2，得到 approximation coefficients $cA_1$ 與 detail coefficients $cD_1$。下一階只繼續分解 $cA_1$。Inverse DWT 則 upsample、通過 synthesis filters 並相加。

# %%
wavelet = pywt.Wavelet("db4")
print("analysis low-pass:", np.round(wavelet.dec_lo, 4))
print("analysis high-pass:", np.round(wavelet.dec_hi, 4))
print("synthesis low-pass:", np.round(wavelet.rec_lo, 4))
print("synthesis high-pass:", np.round(wavelet.rec_hi, 4))
assert len(wavelet.dec_lo) == len(wavelet.dec_hi)


# %% [markdown]
# ## `wavefun(level=...)` 與 decomposition level
#
# `Wavelet.wavefun(level=r)` 的 `level` 是 cascade／refinement 的取樣精細度：增加它會用更多點近似同一個 scaling function 與 mother wavelet。`wavedec(..., level=L)` 的 `L` 才是資料的 decomposition level；兩個參數同名，控制的運算不同。

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
for level in (4, 7):
    _, psi, x_axis = wavelet.wavefun(level=level)
    axes[0].plot(x_axis, psi, label=f"refinement={level}")
axes[0].set_title("same db4 mother wavelet")
axes[0].legend()

ecg = pywt.data.ecg().astype(float)
for level in (1, 3, 5):
    coeffs = pywt.wavedec(ecg, wavelet, level=level, mode="symmetric")
    axes[1].plot([len(coefficient) for coefficient in coeffs], "o-", label=f"decomposition={level}")
axes[1].set_title("coefficient lengths by decomposition")
axes[1].legend()
fig.tight_layout()


# %% [markdown]
# ## 小波家族：support、symmetry 與 vanishing moments 的取捨
#
# 離散小波的 order、vanishing moments、symmetry、filter length 與支撐區間（support）會一起影響係數的局部性與邊緣表現。Haar（`haar`，也是 `db1`）最局部且完全對稱，卻只有一個 vanishing moment。高 order Daubechies（`db`）可消去較高階多項式，但 filter 變長且不對稱。Symlets（`sym`）追求近似對稱；Coiflets（`coif`）同時要求 scaling function 與 wavelet 的 moments；biorthogonal（`bior`）使用不同 analysis／synthesis bases 換取線性相位與對稱性。
#
# CWT 另使用連續小波：Mexican hat（`mexh`）是 real-valued second derivative of Gaussian；Morlet（`morl`）用正弦調變 Gaussian，適合呈現局部振盪。選擇小波家族時，除了係數圖，也要檢查 noise model、boundary effect，以及後續分析結果。

# %%
discrete_families = ["haar", "db4", "sym4", "coif2", "bior2.2"]
continuous_families = ["mexh", "morl"]
fig, axes = plt.subplots(2, 4, figsize=(13, 6))
for ax, name in zip(axes.ravel(), discrete_families + continuous_families):
    if name in continuous_families:
        psi, x_family = pywt.ContinuousWavelet(name).wavefun(level=8)
    else:
        functions = pywt.Wavelet(name).wavefun(level=7)
        psi, x_family = functions[1], functions[-1]
    ax.plot(x_family, psi)
    ax.set_title(name)
    ax.axhline(0, color="black", linewidth=0.5)
axes.ravel()[-1].axis("off")
fig.tight_layout()

assert pywt.Wavelet("db4").vanishing_moments_psi == 4
assert pywt.Wavelet("bior2.2").biorthogonal


# %% [markdown]
# ## 一維 multilevel DWT：每層對應不同尺度
#
# `wavedec(x, level=L)` 回傳 `[cA_L, cD_L, ..., cD_1]`。$cD_1$ 來自第一次 high-pass，對應最細尺度；$cD_L$ 是對 low-pass branch 重複分解後的較粗細節。這是 dyadic filter bank 的尺度分工，不表示每層都是固定頻寬的「理想頻帶」；真實頻率響應由 wavelet filters 決定。PyWavelets 的 [multilevel DWT 文件](https://pywavelets.readthedocs.io/en/stable/ref/dwt-discrete-wavelet-transform.html) 列出係數順序與長度規則。

# %%
ecg_segment = ecg[:512]
multilevel_coefficients = pywt.wavedec(ecg_segment, "db4", level=4, mode="symmetric")
multilevel_reconstruction = pywt.waverec(
    multilevel_coefficients, "db4", mode="symmetric"
)[:ecg_segment.size]
assert np.allclose(multilevel_reconstruction, ecg_segment, atol=1e-10)

coefficient_names = ["cA4", "cD4", "cD3", "cD2", "cD1"]
fig, axes = plt.subplots(len(multilevel_coefficients), 1, figsize=(10, 8))
for ax, coefficient, name in zip(axes, multilevel_coefficients, coefficient_names):
    ax.plot(coefficient)
    ax.set_ylabel(name, rotation=0, labelpad=18)
    ax.set_xlim(0, len(coefficient) - 1)
axes[-1].set_xlabel("coefficient index")
fig.tight_layout()


# %% [markdown]
# (dwt-perfect-reconstruction)=
# ## Perfect reconstruction 與 boundary modes
#
# 在 analysis／synthesis filters 與 extension mode 相容的情況下，不修改係數就能在浮點誤差內重建原訊號。有限訊號仍需延拓邊界；`zero`、`symmetric`、`periodization` 等 mode 會改變係數長度與邊緣係數。分解與重建必須使用相容設定。

# %%
def dwt_round_trip(signal_in, wavelet_name="db4", level=4, mode="symmetric"):
    """分解後不改係數，重建並裁回原長度。"""
    signal_in = np.asarray(signal_in, dtype=float)
    coefficients = pywt.wavedec(signal_in, wavelet_name, level=level, mode=mode)
    reconstructed = pywt.waverec(coefficients, wavelet_name, mode=mode)
    return reconstructed[: signal_in.size], coefficients


rng = np.random.default_rng(3)
test_signal = rng.normal(size=257)
for boundary_mode in ("zero", "symmetric", "periodization"):
    reconstructed, coefficients = dwt_round_trip(test_signal, mode=boundary_mode)
    error = np.max(np.abs(reconstructed - test_signal))
    print(boundary_mode, "coefficient lengths:", [len(c) for c in coefficients], "error:", error)
    assert error < 1e-10

# %% [markdown]
# ```{admonition} Perfect reconstruction 與去雜訊是兩項不同檢查
# :class: note
#
# Perfect reconstruction 只驗證 transform pair 沒有在「不改係數」時丟資料。只要 threshold、truncate 或量化係數，輸出便帶有偏差；是否更接近未知真訊號需要另外評估。
# ```


# %% [markdown]
# ## Gaussian 與 Laplacian pyramids
#
# Gaussian pyramid 在每次 downsampling 前先 low-pass，避免高頻 alias 到低頻。它適合快速瀏覽多種尺度，但單獨保留縮小影像無法完整重建原圖。Laplacian pyramid 改存每層 Gaussian image 與下一層 upsampled image 之差，再加上最粗層；在使用同一個 resize operator 的情況下可完整重建。

# %%
def gaussian_laplacian_pyramid(image_in, levels=4):
    """Build matching Gaussian and Laplacian pyramids."""
    gaussian_levels = [np.asarray(image_in, dtype=float)]
    for _ in range(levels):
        blurred = ndi.gaussian_filter(gaussian_levels[-1], sigma=1.0, mode="reflect")
        gaussian_levels.append(blurred[::2, ::2])
    laplacian_levels = []
    for fine, coarse in zip(gaussian_levels[:-1], gaussian_levels[1:]):
        expanded = ski.transform.resize(
            coarse, fine.shape, order=1, mode="reflect", anti_aliasing=False,
            preserve_range=True,
        )
        laplacian_levels.append(fine - expanded)
    laplacian_levels.append(gaussian_levels[-1])
    return gaussian_levels, laplacian_levels


def reconstruct_laplacian_pyramid(laplacian_levels):
    """Reconstruct using the same interpolation used to build the pyramid."""
    reconstructed = laplacian_levels[-1]
    for detail in reversed(laplacian_levels[:-1]):
        reconstructed = ski.transform.resize(
            reconstructed, detail.shape, order=1, mode="reflect",
            anti_aliasing=False, preserve_range=True,
        ) + detail
    return reconstructed


pyramid_image = ski.transform.resize(
    ski.util.img_as_float(ski.data.camera()), (256, 256), anti_aliasing=True
)
gaussian_levels, laplacian_levels = gaussian_laplacian_pyramid(pyramid_image, levels=3)
pyramid_reconstruction = reconstruct_laplacian_pyramid(laplacian_levels)
assert np.allclose(pyramid_reconstruction, pyramid_image, atol=1e-12)

fig, axes = plt.subplots(2, 4, figsize=(12, 6))
for column, level_image in enumerate(gaussian_levels):
    axes[0, column].imshow(level_image)
    axes[0, column].set_title(f"Gaussian {column}")
for column, level_image in enumerate(laplacian_levels):
    limit = np.max(np.abs(level_image)) if column < 3 else None
    axes[1, column].imshow(level_image, vmin=-limit if limit else None,
                           vmax=limit if limit else None)
    axes[1, column].set_title(f"Laplacian {column}")
for ax in axes.ravel():
    ax.axis("off")
fig.tight_layout()


# %% [markdown]
# ### Multiband blending：低頻平滑過渡，高頻保留局部邊緣
#
# 直接用 binary mask 拼接兩張影像會產生一條高頻接縫。Multiband blending 在每個 Laplacian level 使用對應尺度的 Gaussian mask：粗尺度採寬廣過渡，細尺度則保留局部細節。這項經典應用把「尺度」落實為不同空間頻帶的分工，內容比單純縮放顯示圖更完整。

# %%
left_image = pyramid_image
right_image = ski.transform.resize(
    ski.util.img_as_float(ski.data.coins()), pyramid_image.shape, anti_aliasing=True
)
mask = np.zeros_like(pyramid_image)
mask[:, :mask.shape[1] // 2] = 1.0
_, left_laplacian = gaussian_laplacian_pyramid(left_image, levels=3)
_, right_laplacian = gaussian_laplacian_pyramid(right_image, levels=3)
mask_gaussian, _ = gaussian_laplacian_pyramid(mask, levels=3)
blended_levels = [
    weight * left_level + (1 - weight) * right_level
    for weight, left_level, right_level
    in zip(mask_gaussian, left_laplacian, right_laplacian)
]
multiband_blend = reconstruct_laplacian_pyramid(blended_levels)
direct_blend = mask * left_image + (1 - mask) * right_image
assert multiband_blend.shape == pyramid_image.shape
assert np.isfinite(multiband_blend).all()
imshow_all(direct_blend, multiband_blend,
           titles=["hard seam", "multiband blend"])


# %% [markdown]
# ## 二維 DWT：四個子帶
#
# Separable 2D DWT 先沿一軸、再沿另一軸套 filter bank，一階得到 approximation `cA` 與 horizontal、vertical、diagonal details。不同套件對 `LH/HL` 命名方向可能不同；使用 PyWavelets 時應依 `(cH, cV, cD)` API 語意與測試影像確認方向。

# %%
image = ski.util.img_as_float(ski.data.camera())
cA, (cH, cV, cD) = pywt.dwt2(image, "db4", mode="symmetric")
reconstructed_image = pywt.idwt2((cA, (cH, cV, cD)), "db4", mode="symmetric")
reconstructed_image = reconstructed_image[: image.shape[0], : image.shape[1]]
assert np.allclose(reconstructed_image, image, atol=1e-12)

imshow_all(cA, cH, cV, cD,
           titles=["approximation", "horizontal detail", "vertical detail", "diagonal detail"])


# %% [markdown]
# ### 二維 multilevel DWT
#
# `wavedec2` 回傳 `[cA_L, (cH_L,cV_L,cD_L), ..., (cH_1,cV_1,cD_1)]`。最粗的 approximation 和三個 detail subbands 都受所選 filters、downsampling 與 extension mode 影響，無法只當成原圖縮圖與三張一般「邊緣圖」。未修改係數時，`waverec2` 應在數值誤差內還原影像。

# %%
coefficients_2d = pywt.wavedec2(image, "db4", level=3, mode="symmetric")
reconstruction_2d = pywt.waverec2(coefficients_2d, "db4", mode="symmetric")
reconstruction_2d = reconstruction_2d[:image.shape[0], :image.shape[1]]
assert np.allclose(reconstruction_2d, image, atol=1e-12)
print("2D coefficient shapes:", [
    coefficient.shape if isinstance(coefficient, np.ndarray)
    else tuple(detail.shape for detail in coefficient)
    for coefficient in coefficients_2d
])


# %% [markdown]
# ## Thresholding 需要 sparse-signal 與 noise assumptions
#
# Wavelet shrinkage 假設目標在選定 wavelet domain 中相對 sparse，且許多小 detail coefficients 主要來自可建模雜訊。Hard threshold 將小係數設為零、保留大係數，在 threshold 處不連續；soft threshold 再將保留的大係數往零收縮，較平滑但會產生 bias。
#
# 對長度 $N$ 的白 Gaussian noise，universal threshold 常寫成 $\tau=\hat\sigma\sqrt{2\log N}$；$\hat\sigma$ 可用最細層 diagonal coefficients 的 median absolute deviation（MAD）估計：$\hat\sigma=\operatorname{median}(|cD_1|)/0.6745$。BayesShrink 改為每個子帶估計 threshold。這些公式仍需配合資料調整；遇到空間相關雜訊、CTF-shaped spectrum 或非平穩變異時，簡單白雜訊假設可能失準，弱小但真實的高解析訊號也可能被刪掉。

# %%
original = ski.util.img_as_float(ski.data.camera())
rng = np.random.default_rng(8)
sigma_true = 0.12
noisy = np.clip(original + rng.normal(scale=sigma_true, size=original.shape), 0, 1)
sigma_estimate = float(estimate_sigma(noisy, channel_axis=None))
noisy_coefficients = pywt.wavedec2(noisy, "db4", level=3, mode="symmetric")
finest_diagonal = noisy_coefficients[-1][2]
sigma_mad = np.median(np.abs(finest_diagonal)) / 0.6745
universal_threshold = sigma_mad * np.sqrt(2 * np.log(noisy.size))


def threshold_detail_coefficients(coefficients, threshold, mode):
    """Threshold every detail subband, leaving the approximation unchanged."""
    return [coefficients[0]] + [
        tuple(pywt.threshold(detail, threshold, mode=mode) for detail in level)
        for level in coefficients[1:]
    ]


denoised_hard = pywt.waverec2(
    threshold_detail_coefficients(noisy_coefficients, universal_threshold, "hard"),
    "db4", mode="symmetric",
)[:original.shape[0], :original.shape[1]]
denoised_soft = pywt.waverec2(
    threshold_detail_coefficients(noisy_coefficients, universal_threshold, "soft"),
    "db4", mode="symmetric",
)[:original.shape[0], :original.shape[1]]
denoised_bayes = denoise_wavelet(
    noisy,
    method="BayesShrink",
    mode="soft",
    wavelet="db4",
    rescale_sigma=True,
    channel_axis=None,
)

denoising_results = {
    "noisy": noisy,
    "universal hard": np.clip(denoised_hard, 0, 1),
    "universal soft": np.clip(denoised_soft, 0, 1),
    "BayesShrink": np.clip(denoised_bayes, 0, 1),
}
quality = {
    name: (
        peak_signal_noise_ratio(original, result, data_range=1),
        structural_similarity(original, result, data_range=1),
    )
    for name, result in denoising_results.items()
}

print(f"sigma: estimate_sigma={sigma_estimate:.3f}, MAD={sigma_mad:.3f}")
for name, (psnr_value, ssim_value) in quality.items():
    print(f"{name:>14s}: PSNR={psnr_value:.2f} dB, SSIM={ssim_value:.3f}")
assert abs(sigma_mad - sigma_true) < 0.03
assert quality["BayesShrink"][0] > quality["noisy"][0]
assert quality["BayesShrink"][1] > quality["noisy"][1]

fig, axes = plt.subplots(2, 3, figsize=(12, 8))
display_names = ["noisy", "universal hard", "universal soft", "BayesShrink"]
axes[0, 0].imshow(original)
axes[0, 0].set_title("reference")
for ax, name in zip(axes.ravel()[1:5], display_names):
    psnr_value, ssim_value = quality[name]
    ax.imshow(denoising_results[name])
    ax.set_title(f"{name}\n{psnr_value:.1f} dB, SSIM {ssim_value:.3f}")
residual = denoising_results["BayesShrink"] - original
residual_limit = np.max(np.abs(residual))
axes[1, 2].imshow(residual, cmap="coolwarm", vmin=-residual_limit, vmax=residual_limit)
axes[1, 2].set_title("BayesShrink residual")
for ax in axes.ravel():
    ax.axis("off")
fig.tight_layout()


# %% [markdown]
# 峰值訊雜比（PSNR）是 MSE 的對數量尺，SSIM 比較局部亮度、對比與結構；殘差圖（residual）則顯示去雜訊方法從參考影像刪掉或加入了哪些內容。三者適合用在有已知真值的教學評估，實際分析仍要加上任務相關的檢查。本例使用固定隨機數產生器，所以不同方法看到同一份雜訊。`estimate_sigma` 與 `denoise_wavelet` 的參數定義見 [scikit-image restoration 文件](https://scikit-image.org/docs/stable/api/skimage.restoration.html)。


# %% [markdown]
# ## Decimated DWT 的 shift variance 與 SWT
#
# DWT 每層都 downsample，因此輸入只平移一點，也可能大幅改變哪些係數被保留。Stationary wavelet transform（SWT）不做 downsampling，改為逐層在 filters 中插零；它保留每層與輸入同長的係數，在 periodized 情況下對整數平移具 equivariance，代價是冗餘儲存與計算量。

# %%
shift_signal = signal_sequential[:128]
shifted_signal = np.roll(shift_signal, 1)
dwt_original = pywt.wavedec(shift_signal, "db2", level=3, mode="periodization")
dwt_shifted = pywt.wavedec(shifted_signal, "db2", level=3, mode="periodization")
swt_original = pywt.swt(shift_signal, "db2", level=3)
swt_shifted = pywt.swt(shifted_signal, "db2", level=3)

assert not np.allclose(dwt_shifted[-1], np.roll(dwt_original[-1], 1))
for (approximation, detail), (shifted_approximation, shifted_detail) in zip(
    swt_original, swt_shifted
):
    assert np.allclose(shifted_approximation, np.roll(approximation, 1), atol=1e-10)
    assert np.allclose(shifted_detail, np.roll(detail, 1), atol=1e-10)

fig, axes = plt.subplots(2, 2, figsize=(11, 5))
axes[0, 0].plot(dwt_original[-1], label="original")
axes[0, 0].plot(dwt_shifted[-1], "--", label="shifted input")
axes[0, 0].set_title("DWT finest detail")
axes[0, 0].legend()
axes[0, 1].plot(np.roll(dwt_original[-1], 1) - dwt_shifted[-1])
axes[0, 1].set_title("DWT: not a one-index roll")
axes[1, 0].plot(swt_original[-1][1], label="original")
axes[1, 0].plot(swt_shifted[-1][1], "--", label="shifted input")
axes[1, 0].set_title("SWT finest detail")
axes[1, 0].legend()
axes[1, 1].plot(np.roll(swt_original[-1][1], 1) - swt_shifted[-1][1])
axes[1, 1].set_title("SWT: aligned difference")
fig.tight_layout()


# %% [markdown]
# ```{admonition} 與 cryo-EM 的連結
# :class: caution
#
# Gaussian／Laplacian pyramids 可幫助分開粒子的整體輪廓與較局部細節；小波可作多尺度視覺化、特定稀疏先驗下的正則化，或經驗證的去雜訊元件。實際資料中的結構未必集中在少數大係數，雜訊係數也可能彼此相關。CTF、運動、冰層與偵測器響應都會造成相關且頻率相依的訊號／雜訊。
#
# 對顯示用影像做增強可以幫助人眼檢查，但未經獨立 half-set 或已知真值驗證，不要把增強後的粒子影像直接當成定量重建輸入，也不要把「看起來更銳利」解讀成解析度提升。
# ```
#
# ## 理解檢查
#
# 1. STFT 與 CWT 如何在時間定位和頻率解析度之間取捨？
#
#    ```{dropdown} 參考答案
#    STFT 在每個時間位置使用同一長度的 window。若 window 含 64 個 samples，時間定位約受這 64 點寬度限制；改成 256 點後，frequency-bin spacing 約縮小為原來的四分之一，但短暫事件會在較寬的時間範圍內出現。所有頻率都共用這組取捨。
#
#    CWT 會縮放 mother wavelet。小尺度支撐範圍窄，適合定位快速變化；大尺度支撐範圍寬，可分辨較慢的振盪。Pseudo-frequency 還要由母小波的 central frequency 與 `sampling_period` 換算，$1/a$ 只提供尺度 $a$ 與頻率的反比趨勢。常見誤解是把 scalogram 上每個亮點都視為精確頻率；window／wavelet、取樣率與邊界效應都會影響亮點的寬度和位置。
#    ```
#
# 2. Gaussian／Laplacian pyramid、STFT、CWT 與 DWT 各自如何描述位置與尺度？
#
#    ```{dropdown} 參考答案
#    Gaussian pyramid 反覆做 low-pass filtering 與 downsampling：
#
#    $$
#    G_{j+1}=\downarrow_2\,(g*G_j).
#    $$
#
#    每層保留影像平面上的位置，但解析度與取樣密度逐層降低。Laplacian pyramid 用相鄰 Gaussian levels 的差異保留頻帶細節，配合最粗層可重建原圖。以一張同時含小亮點與寬廣亮斑的影像為例，小亮點主要出現在細層 Laplacian bands，寬廣亮斑會延續到較粗層。
#
#    STFT 用固定 window 取得局部頻譜，所有頻率共用相同的位置／頻率解析度。CWT 縮放 mother wavelet，對每個尺度密集取樣位置，因此高頻事件定位較細、低頻事件的頻率分辨較好。DWT 以 dyadic filters 與 downsampling 產生緊密的 multiresolution coefficients，儲存量較省，但 decimation 會帶來 shift sensitivity。這四種方法的 level 或 scale 都受 filter、取樣率與邊界設定影響，不能直接當成物理解析度；要依目標是影像縮放、局部頻率分析、事件偵測或可逆壓縮來選擇。
#    ```
#
# 3. DWT 可以 perfect reconstruction，為何 thresholding 仍可能刪掉訊號？輸入平移一個 sample 又會發生什麼事？
#
#    ```{dropdown} 參考答案
#    Perfect reconstruction 指 analysis 與 synthesis filters 滿足重建條件。係數未修改時，`waverec(wavedec(x))` 可在浮點誤差內還原 $x$。Thresholding 主動把部分 detail coefficients 設成零或縮小，重建輸入已經改變；若弱邊緣與雜訊都落在小係數區，兩者會一起被刪除。重建誤差此時來自係數修改，與 filter bank 是否可逆是兩回事。
#
#    Decimated DWT 每層 filtering 後取每隔一點的樣本。輸入平移一個 sample 會改變 downsampling 選到的奇偶位置，許多係數都可能改變，結果未必是原係數平移一格。SWT 省略 downsampling，shift sensitivity 通常較低，代價是係數量增加。常見誤解是把 perfect reconstruction 當成去雜訊品質保證；它只檢查未修改係數的可逆性。
#    ```
#
# 4. 如何判斷 cryo-EM 去雜訊方法是否真的改善分析結果？
#
#    ```{dropdown} 參考答案
#    有已知真值的模擬資料可先量 MSE、PSNR、SSIM，並查看殘差圖是否含有粒子結構。假設 clean image 的資料範圍是 $[0,1]$，MSE 從 0.010 降到 0.0025，則
#
#    $$
#    \Delta\mathrm{PSNR}=10\log_{10}(0.010/0.0025)\approx6.02\ \mathrm{dB}.
#    $$
#
#    這個數字說明像素誤差下降，還要確認解析度相關資訊有沒有被平滑掉。真實資料沒有逐張 clean target，兩個 half-sets 應分開去雜訊與估計參數，並檢查訓練、調參或共享模型是否造成跨 half-set 資訊洩漏。接著可比較對位穩定度、2D class averages、重建 FSC，以及同一生物結構特徵能否在兩個 half-maps 重現。只看影像變銳利或 PSNR 上升，無法回答姿態估計與重建是否改善；評估項目要配合方法預計改善的分析步驟。
#    ```
