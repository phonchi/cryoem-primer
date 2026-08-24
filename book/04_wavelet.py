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
# 全域 Fourier spectrum 告訴我們整段訊號的頻率組成，但不保留事件發生的位置。小波把分析函數限制在局部，再以不同尺度掃描訊號。本章分清 continuous wavelet transform（CWT）與 decimated discrete wavelet transform（DWT）；兩者都叫小波轉換，輸出與用途卻不同。
#
# ```{admonition} 學習目標
# :class: important
#
# - 說明有限觀測下 Fourier spectrum 的解析度與 spectral leakage。
# - 寫出 CWT 與 dyadic wavelet family，說明位移網格如何隨尺度改變。
# - 將 DWT 解釋成 analysis filter bank 加 downsampling。
# - 驗證 DWT perfect reconstruction，並比較 boundary modes。
# - 說明小波閾值去雜訊依賴哪些假設，以及何時不該用於定量重建。
# ```

# %%
import matplotlib.pyplot as plt
import numpy as np
import pywt
import skimage as ski
from skimage.metrics import peak_signal_noise_ratio
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
# 有限長度 DFT 的 frequency-bin spacing 由觀測時間決定；非整數週期還會產生 spectral leakage。因此 Fourier spectrum 提供的是整段訊號的頻率組成，不是無條件「精確列出所有頻率」。
#
# 下例中訊號 A 讓四個頻率全程同時存在，訊號 B 則讓它們依序出現。兩者都在相同頻率附近有能量，但 B 的每個成分只持續四分之一秒，峰值較低，分段邊界也造成明顯 leakage；不能說兩張 magnitude spectra 幾乎相同。

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
# Short-time Fourier transform（STFT）用固定長度 window 切出局部頻譜，因此時間／頻率解析度由同一個 window 決定。CWT 改用可縮放的母小波 $\psi$：
#
# $$W(a,b)=\frac{1}{\sqrt{|a|}}\int h(t)\psi^*\!\left(\frac{t-b}{a}\right)dt.$$
#
# 大尺度通常對應較低 pseudo-frequency 與較寬的時間支撐；小尺度對應較高 pseudo-frequency 與較窄的時間支撐。CWT 對每個尺度保留密集的位移取樣，結果可畫成尺度圖（scalogram）。尺度和 Hz 的換算依母小波與取樣週期而定，不能只寫成「尺度就是頻率倒數」。

# %%
scales = np.arange(1, 65)
cwt_coefficients, cwt_frequencies = pywt.cwt(
    signal_sequential, scales, "morl", sampling_period=1 / sample_rate
)
assert cwt_coefficients.shape == (len(scales), len(signal_sequential))

fig, ax = plt.subplots(figsize=(10, 4))
extent = [t[0], t[-1], cwt_frequencies[-1], cwt_frequencies[0]]
ax.imshow(np.abs(cwt_coefficients), aspect="auto", extent=extent, origin="upper")
ax.set_xlabel("time [s]")
ax.set_ylabel("pseudo-frequency [Hz]")
ax.set_title("CWT scalogram")
fig.tight_layout()


# %% [markdown]
# ## Dyadic family 與 DWT
#
# 一個常見的 dyadic wavelet family 可寫成
#
# $$\psi_{j,k}(t)=2^{-j/2}\psi(2^{-j}t-k).$$
#
# 以 CWT 參數理解，就是 $a=2^j$、$b=k2^j$：尺度變大時，位移網格也變粗，不是每個尺度都任取同一組整數位移。正交／雙正交 DWT 再利用特定 scaling function 與 wavelet，形成可逆、非冗餘或低冗餘的離散表示。
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
# ## `wavefun(level=...)` 不是 decomposition level
#
# `Wavelet.wavefun(level=r)` 的 `level` 是 cascade／refinement 的取樣精細度：增加它會用更多點近似同一個 scaling function 與 mother wavelet，並不表示對某筆資料做了 $r$ 階 DWT，也不是把母小波「拉伸」。真正的 decomposition level 是 `wavedec(..., level=L)` 的 `L`。

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
# 小波家族的 order、vanishing moments、symmetry 與支撐區間（support）會一起影響係數的局部性與邊緣表現。沒有一個家族對所有訊號都最好；選擇時要連同 noise model 與下游任務驗證。


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
# ```{admonition} Perfect reconstruction 不等於去雜訊成功
# :class: note
#
# Perfect reconstruction 只驗證 transform pair 沒有在「不改係數」時丟資料。只要 threshold、truncate 或量化係數，輸出便帶有偏差；是否更接近未知真訊號需要另外評估。
# ```


# %% [markdown]
# ## 二維 DWT：四個子帶
#
# Separable 2D DWT 先沿一軸、再沿另一軸套 filter bank，一階得到 approximation `cA` 與 horizontal、vertical、diagonal details。不同套件對 `LH/HL` 命名方向可能不同；使用 PyWavelets 時應依 `(cH, cV, cD)` API 語意，不只看字母猜方向。

# %%
image = ski.util.img_as_float(ski.data.camera())
cA, (cH, cV, cD) = pywt.dwt2(image, "db4", mode="symmetric")
reconstructed_image = pywt.idwt2((cA, (cH, cV, cD)), "db4", mode="symmetric")
reconstructed_image = reconstructed_image[: image.shape[0], : image.shape[1]]
assert np.allclose(reconstructed_image, image, atol=1e-12)

imshow_all(cA, cH, cV, cD,
           titles=["approximation", "horizontal detail", "vertical detail", "diagonal detail"])


# %% [markdown]
# ## Thresholding 需要 sparse-signal 與 noise assumptions
#
# Wavelet shrinkage 假設目標在選定 wavelet domain 中相對 sparse，且許多小 detail coefficients 主要來自可建模雜訊。Soft threshold 會把大係數也往零收縮，因此會產生 bias。若雜訊有空間相關性、CTF-shaped spectrum 或非平穩變異，簡單的白 Gaussian noise threshold 未必合適；弱小但真實的高解析訊號也可能被刪掉。

# %%
original = ski.util.img_as_float(ski.data.camera())
rng = np.random.default_rng(8)
sigma_true = 0.12
noisy = np.clip(original + rng.normal(scale=sigma_true, size=original.shape), 0, 1)
sigma_estimate = float(estimate_sigma(noisy, channel_axis=None))
denoised = denoise_wavelet(
    noisy,
    method="BayesShrink",
    mode="soft",
    wavelet="db4",
    rescale_sigma=True,
    channel_axis=None,
)

psnr_noisy = peak_signal_noise_ratio(original, noisy, data_range=1)
psnr_denoised = peak_signal_noise_ratio(original, denoised, data_range=1)
print(f"sigma estimate={sigma_estimate:.3f}")
print(f"PSNR noisy={psnr_noisy:.2f} dB, denoised={psnr_denoised:.2f} dB")
assert psnr_denoised > psnr_noisy

imshow_all(original, noisy, denoised,
           titles=["reference", f"noisy ({psnr_noisy:.1f} dB)", f"denoised ({psnr_denoised:.1f} dB)"])


# %% [markdown]
# ```{admonition} 與 cryo-EM 的連結
# :class: caution
#
# 小波可作多尺度視覺化、特定稀疏先驗下的 regularization，或經驗證的去雜訊元件；不能預設真實 cryo-EM 結構必然集中在少數大係數，也不能預設雜訊必然是許多獨立小係數。CTF、motion、ice 與 detector response 都會造成相關且頻率相依的訊號／雜訊。
#
# 對顯示用影像做增強可以幫助人眼檢查，但未經 half-set 或 ground-truth 驗證，不要把增強後 particles 直接當成定量重建輸入，也不要把「看起來更銳利」解讀成解析度提升。
# ```
#
# ## 理解檢查
#
# 1. 為什麼訊號 A 與 B 在相同頻率附近有能量，卻不能說 magnitude spectra 幾乎相同？
# 2. 在 dyadic family 中，尺度變成兩倍時，允許的位移網格如何改變？
# 3. `wavefun(level=7)` 與 `wavedec(level=7)` 的 `level` 各代表什麼？
# 4. DWT 能 perfect reconstruction，為什麼 thresholding 後仍可能刪掉真實訊號？
# 5. 在 cryo-EM 中，什麼證據才能支持「去雜訊改善了分析」而不只是改善視覺外觀？
