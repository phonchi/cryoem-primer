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
# # 頻域影像處理：從 DFT 到 CTF correction
#
# 傅立葉轉換把「影像裡哪裡亮」改寫成「各空間頻率有多少振幅與相位」。Cryo-EM 的 projection、CTF、平移與 reconstruction 都在這套語言中變得比較直接；前提是索引、正規化與 boundary convention 沒有混用。
#
# ```{admonition} 學習目標
# :class: important
#
# - 寫出矩形影像的 2D DFT／inverse DFT，核對 DC、Parseval 與 Hermitian symmetry。
# - 解釋有限視窗造成的 spectral leakage，以及 padding 為何影響 convolution。
# - 分清 inverse filtering、regularized inversion、phase flipping 與 Wiener-style CTF correction。
# - 用 Fourier slice theorem 連結 2D projections 與 3D reconstruction。
# - 說出簡化 CTF 省略了哪些物理因素。
# ```

# %%
import matplotlib.pyplot as plt
import numpy as np
import skimage as ski
from scipy import ndimage as ndi
from scipy import signal

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
# ## 一維 DFT：複數係數保留振幅與相位
#
# 對長度 $N$ 的訊號 $h[n]$，NumPy 採用未正規化 forward DFT：
#
# $$H[k]=\sum_{n=0}^{N-1}h[n]e^{-i2\pi kn/N},\qquad
# h[n]=\frac1N\sum_{k=0}^{N-1}H[k]e^{i2\pi kn/N}.$$
#
# 實數訊號的頻譜具有 Hermitian symmetry：$H[-k]=H[k]^*$。因此 `rfft` 只需儲存非負頻率半邊。

# %%
sample_rate = 100.0
t = np.arange(200) / sample_rate
x = np.sin(2 * np.pi * 10 * t + 0.3)
X = np.fft.fft(x)
frequencies = np.fft.fftfreq(x.size, d=1 / sample_rate)

assert np.allclose(np.fft.ifft(X).real, x)
assert np.allclose(X[1:], np.conj(X[:0:-1]))

fig, axes = plt.subplots(1, 2, figsize=(11, 3.5))
axes[0].plot(t[:60], x[:60])
axes[0].set_xlabel("time [s]")
axes[0].set_title("10 Hz signal")
axes[1].stem(frequencies, np.abs(X))
axes[1].set_xlim(-50, 50)
axes[1].set_xlabel("frequency [Hz]")
axes[1].set_title("magnitude spectrum")
fig.tight_layout()


# %% [markdown]
# ## 矩形影像的二維 DFT
#
# 令影像 $h[y,x]$ 有 $N$ 列、$M$ 欄。兩個軸必須各用自己的長度：
#
# $$H[k_y,k_x]=\sum_{y=0}^{N-1}\sum_{x=0}^{M-1}
# h[y,x]e^{-i2\pi(k_yy/N+k_xx/M)},$$
#
# $$h[y,x]=\frac1{NM}\sum_{k_y=0}^{N-1}\sum_{k_x=0}^{M-1}
# H[k_y,k_x]e^{i2\pi(k_yy/N+k_xx/M)}.$$
#
# 這和陣列 `[row, column]` 對應：`k_y` 走列方向，`k_x` 走欄方向。

# %%
def direct_dft2(image):
    """教學用矩陣式 2D DFT；只適合小陣列。"""
    image = np.asarray(image, dtype=complex)
    n_rows, n_cols = image.shape
    y = np.arange(n_rows)
    x = np.arange(n_cols)
    basis_y = np.exp(-2j * np.pi * np.outer(y, y) / n_rows)
    basis_x = np.exp(-2j * np.pi * np.outer(x, x) / n_cols)
    return basis_y @ image @ basis_x


rectangular = np.arange(15, dtype=float).reshape(3, 5)
H_direct = direct_dft2(rectangular)
H_fft = np.fft.fft2(rectangular)
round_trip = np.fft.ifft2(H_fft).real

assert np.allclose(H_direct, H_fft)
assert np.allclose(round_trip, rectangular)
assert np.isclose(H_fft[0, 0], rectangular.sum())
assert np.isclose(H_fft[0, 0] / rectangular.size, rectangular.mean())


# %% [markdown]
# ### DC、Parseval 與 Hermitian symmetry
#
# 在這個慣例下，`H[0,0]` 是所有像素的**總和**；除以 $NM$ 才是平均值。Parseval identity 則把兩域的能量連起來：
#
# $$\sum_{y,x}|h[y,x]|^2=\frac1{NM}\sum_{k_y,k_x}|H[k_y,k_x]|^2.$$
#
# 對實數 2D 影像，頻率座標的對稱關係是 $H[-k_y,-k_x]=H[k_y,k_x]^*$，索引要以各軸長度取模。

# %%
energy_space = np.sum(np.abs(rectangular) ** 2)
energy_frequency = np.sum(np.abs(H_fft) ** 2) / rectangular.size
assert np.isclose(energy_space, energy_frequency)

n_rows, n_cols = rectangular.shape
for ky in range(n_rows):
    for kx in range(n_cols):
        assert np.allclose(H_fft[-ky % n_rows, -kx % n_cols], np.conj(H_fft[ky, kx]))
print("rectangular DFT checks passed")


# %% [markdown]
# ## `fftshift` 只改顯示順序
#
# FFT 輸出的 DC 位於 `[0,0]`。`fftshift()` 把它移到畫面中央，方便看低頻到高頻的空間排列；它不會新增資訊，也不會改變頻率係數。

# %%
camera = ski.util.img_as_float(ski.data.camera())
camera_fft = np.fft.fft2(camera)
log_magnitude = np.log1p(np.abs(np.fft.fftshift(camera_fft)))
phase = np.angle(np.fft.fftshift(camera_fft))
imshow_all(camera, log_magnitude, phase,
           titles=["image", "log magnitude", "phase"])


# %% [markdown]
# ## 有限觀測與 spectral leakage
#
# DFT 把有限長度訊號視為週期重複。若觀測窗內不是整數週期，首尾接合會產生不連續，能量便散到鄰近 frequency bins，稱為 spectral leakage。加窗能降低遠端 leakage，但會拓寬主峰；zero-padding 只把頻譜取樣畫得更密，不會提高由觀測長度決定的真實解析力。

# %%
n = 128
sample_index = np.arange(n)
non_bin_tone = np.sin(2 * np.pi * 10.35 * sample_index / n)
rect_spectrum = np.abs(np.fft.rfft(non_bin_tone))
hann_spectrum = np.abs(np.fft.rfft(non_bin_tone * np.hanning(n)))

fig, ax = plt.subplots(figsize=(8, 3.5))
ax.semilogy(rect_spectrum / rect_spectrum.max(), label="rectangular window")
ax.semilogy(hann_spectrum / hann_spectrum.max(), label="Hann window")
ax.set_xlabel("frequency bin")
ax.set_ylabel("normalized magnitude")
ax.legend()
fig.tight_layout()


# %% [markdown]
# (convolution)=
# ## Linear convolution 與 circular convolution
#
# 直接把同尺寸 FFT 相乘再 inverse FFT，得到的是 circular convolution：超出右端的訊號會繞回左端。要取得長度 $N+K-1$ 的 linear convolution，兩個輸入都要 zero-pad 到至少這個長度。

# %%
def fft_linear_convolve_1d(x, h):
    """以足夠 zero-padding 計算 full linear convolution。"""
    output_length = len(x) + len(h) - 1
    return np.fft.ifft(np.fft.fft(x, output_length) * np.fft.fft(h, output_length)).real


x_small = np.array([1.0, 0.0, 0.0, 2.0])
h_small = np.array([1.0, 1.0, 1.0])
linear = fft_linear_convolve_1d(x_small, h_small)
circular = np.fft.ifft(np.fft.fft(x_small) * np.fft.fft(h_small, len(x_small))).real
assert np.allclose(linear, np.convolve(x_small, h_small, mode="full"))
assert not np.allclose(circular, linear[: len(circular)])
print("linear:", linear, "circular:", circular)


# %% [markdown]
# ## Inverse filtering 是病態問題
#
# 若 $G=HF+N$，沒有雜訊時可以形式上寫 $F=G/H$。但 $H$ 接近零時，任何微小的 $N$ 都會被放大。程式裡直接寫 `G / (H + eps)` 已經不是 exact inverse；`eps` 改變了 transfer function，且可能連相位也一起偏移。
#
# 較透明的示範是 truncated inverse：只在 $|H|\ge\tau$ 的頻率做除法，其餘設為零。$\tau$ 是明確的 regularization parameter；結果不應稱為「完美復原」。

# %%
def psf_to_otf(psf, shape):
    """將以中心表示的 PSF padding 並移到 FFT origin。"""
    padded = np.zeros(shape, dtype=float)
    slices = tuple(slice(0, size) for size in psf.shape)
    padded[slices] = psf
    for axis, size in enumerate(psf.shape):
        padded = np.roll(padded, -(size // 2), axis=axis)
    return np.fft.fft2(padded)


def truncated_inverse(observed, transfer, threshold):
    """只反轉 transfer magnitude 不低於 threshold 的頻率。"""
    observed_fft = np.fft.fft2(observed)
    inverse = np.zeros_like(transfer, dtype=complex)
    stable = np.abs(transfer) >= threshold
    inverse[stable] = 1.0 / transfer[stable]
    return np.fft.ifft2(observed_fft * inverse).real, stable


original = ski.util.img_as_float(ski.data.camera())[::2, ::2]
g1 = signal.windows.gaussian(11, std=2)
psf = np.outer(g1, g1)
psf /= psf.sum()
transfer = psf_to_otf(psf, original.shape)
blurred = np.fft.ifft2(np.fft.fft2(original) * transfer).real
rng = np.random.default_rng(0)
noisy_blurred = blurred + rng.normal(scale=0.01 * blurred.std(), size=blurred.shape)
restored_truncated, retained = truncated_inverse(noisy_blurred, transfer, threshold=0.08)

assert 0 < retained.mean() < 1
print("retained frequency fraction:", retained.mean())
imshow_all(original, noisy_blurred, restored_truncated,
           titles=["original", "blurred + noise", "truncated inverse"])


# %% [markdown]
# ## Wiener regularization
#
# 在簡化的白雜訊／平穩訊號模型下，常見的 Wiener-style estimator 可寫成
#
# $$\widehat F=\frac{H^*}{|H|^2+K}G.$$
#
# $K$ 代表 noise-to-signal power 的近似；$K$ 越大，越不願意在 transfer 弱的頻率放大觀測值。這是 bias–variance trade-off，不會創造零點處沒有被量到的資訊。

# %%
def wiener_frequency(observed, transfer, regularization):
    if regularization <= 0:
        raise ValueError("regularization must be positive")
    observed_fft = np.fft.fft2(observed)
    estimator = np.conj(transfer) / (np.abs(transfer) ** 2 + regularization)
    return np.fft.ifft2(estimator * observed_fft).real


restored_wiener = wiener_frequency(noisy_blurred, transfer, regularization=2e-3)
assert np.isfinite(restored_wiener).all()
imshow_all(noisy_blurred, restored_truncated, restored_wiener,
           titles=["observed", "truncated inverse", "Wiener-style"])


# %% [markdown]
# ## Fourier slice theorem
#
# 對 2D 物件沿 $y$ 積分得到 1D projection $p[x]=\sum_y f[y,x]$，$p$ 的 1D Fourier transform 等於 $f$ 的 2D Fourier transform 中 $k_y=0$ 的 central line。3D 情況同理：每張 2D projection 的 Fourier transform 對應 3D Fourier volume 中一個通過原點的平面；粒子取向決定平面的方向。

# %%
phantom = ski.data.shepp_logan_phantom()
phantom = ski.transform.resize(phantom, (96, 128), anti_aliasing=True)
projection = phantom.sum(axis=0)
projection_fft = np.fft.fft(projection)
central_slice = np.fft.fft2(phantom)[0, :]
assert np.allclose(projection_fft, central_slice, atol=1e-10)

fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
axes[0].imshow(phantom)
axes[0].set_title("2D object")
axes[0].axis("off")
axes[1].plot(np.abs(np.fft.fftshift(projection_fft)), label="FFT of projection")
axes[1].plot(np.abs(np.fft.fftshift(central_slice)), "--", label="central slice")
axes[1].legend()
fig.tight_layout()

# %% [markdown]
# ```{dropdown} 從 theorem 到 SPA reconstruction 還缺什麼？
#
# 真實 projection 有 CTF、平移、雜訊、有限 box、未知取向與可能的構型差異。Fourier slice theorem 說明資料幾何，沒有自動解決 pose estimation、interpolation、regularization 或 heterogeneity。
# ```


# %% [markdown]
# ## 從簡化 CTF 到較完整的參數
#
# 教學上常先用無散光、無振幅對比、無 envelope 的純相位模型
#
# $$\mathrm{CTF}(s)=-\sin\chi(s),\qquad
# \chi(s)=\pi\lambda\Delta f\,s^2-\frac{\pi}{2}C_s\lambda^3s^4+\phi.$$
#
# 符號會隨 defocus 與 Fourier convention 改寫；真正重要的是全章採同一慣例。較完整的模型會混合 amplitude contrast，並以 envelope 描述 temporal／spatial coherence、motion 與偵測器造成的高頻衰減。Astigmatism 讓 defocus 隨方位角改變；phase plate 會加入額外 phase shift。

# %%
def electron_wavelength_A(voltage_kv):
    """Relativistic electron wavelength in Å for accelerating voltage in kV。"""
    if voltage_kv <= 0:
        raise ValueError("voltage_kv must be positive")
    voltage = voltage_kv * 1_000.0
    return 12.2639 / np.sqrt(voltage * (1.0 + 0.97845e-6 * voltage))


def ctf_1d(spatial_frequency, defocus_A, voltage_kv=300.0, cs_mm=2.7,
           amplitude_contrast=0.0, phase_shift_rad=0.0, b_factor_A2=0.0):
    """Isotropic CTF；不含 astigmatism，envelope 以可選 B factor 近似。"""
    s = np.asarray(spatial_frequency, dtype=float)
    if not 0 <= amplitude_contrast < 1:
        raise ValueError("amplitude_contrast must be in [0, 1)")
    wavelength = electron_wavelength_A(voltage_kv)
    cs_A = cs_mm * 1e7
    chi = (
        np.pi * wavelength * defocus_A * s**2
        - 0.5 * np.pi * cs_A * wavelength**3 * s**4
        + phase_shift_rad
    )
    phase_weight = np.sqrt(1.0 - amplitude_contrast**2)
    envelope = np.exp(-0.25 * b_factor_A2 * s**2)
    return -envelope * (phase_weight * np.sin(chi) + amplitude_contrast * np.cos(chi))


assert np.isclose(electron_wavelength_A(300), 0.01969, rtol=2e-3)
s = np.linspace(0, 0.5, 4000)
ctf_simple = ctf_1d(s, defocus_A=15_000, voltage_kv=300, cs_mm=2.7)
ctf_with_terms = ctf_1d(
    s, defocus_A=15_000, voltage_kv=300, cs_mm=2.7,
    amplitude_contrast=0.1, phase_shift_rad=0.15, b_factor_A2=40,
)

fig, ax = plt.subplots(figsize=(9, 3.5))
ax.plot(s, ctf_simple, label="pure phase, no envelope")
ax.plot(s, ctf_with_terms, label="amplitude contrast + phase shift + envelope")
ax.axhline(0, color="black", linewidth=0.6)
ax.set_xlabel("spatial frequency [1/Å]")
ax.set_ylabel("CTF")
ax.legend()
fig.tight_layout()


# %% [markdown]
# 若固定電壓、$C_s$、振幅對比、phase shift 與散光設定，改變離焦（defocus）會移動 CTF 零點。但零點不由離焦單獨決定；電子波長、球差、振幅對比與額外 phase 都可能改變它。散光則使零點依方位角而異。


# %% [markdown]
# ## Phase flipping 與 Wiener-style CTF correction 是兩件事
#
# 若觀測頻譜 $G=C F+N$：
#
# - phase flipping 使用 $\operatorname{sign}(C)G$，只修正 sign／phase reversal，不補償 $|C|$ 的衰減。
# - Wiener-style correction 使用 $C^*/(|C|^2+K)G$，同時做帶正則化的振幅校正；$K$ 決定零點附近抑制程度。
#
# 以下以 2D isotropic CTF 示範。同一張影像的 CTF 零點附近無法可靠復原；多個離焦組可以提供互補觀測，但合併仍需正確的 pose、noise weighting 與 validation。

# %%
def radial_frequency_grid(shape, pixel_size_A):
    fy = np.fft.fftfreq(shape[0], d=pixel_size_A)
    fx = np.fft.fftfreq(shape[1], d=pixel_size_A)
    yy, xx = np.meshgrid(fy, fx, indexing="ij")
    return np.sqrt(xx**2 + yy**2)


def phase_flip_spectrum(observed_fft, ctf_values):
    signs = np.sign(ctf_values)
    return signs * observed_fft


def wiener_ctf_spectrum(observed_fft, ctf_values, regularization):
    if regularization <= 0:
        raise ValueError("regularization must be positive")
    return np.conj(ctf_values) * observed_fft / (np.abs(ctf_values) ** 2 + regularization)


true_projection = ski.transform.resize(ski.data.camera(), (192, 192), anti_aliasing=True)
true_projection = ski.util.img_as_float(true_projection)
frequency_radius = radial_frequency_grid(true_projection.shape, pixel_size_A=1.5)
ctf_grid = ctf_1d(
    frequency_radius, defocus_A=15_000, voltage_kv=300,
    cs_mm=2.7, amplitude_contrast=0.1, b_factor_A2=30,
)
true_fft = np.fft.fft2(true_projection)
rng = np.random.default_rng(7)
noise = rng.normal(scale=0.03 * true_projection.std(), size=true_projection.shape)
observed = np.fft.ifft2(ctf_grid * true_fft).real + noise
observed_fft = np.fft.fft2(observed)

phase_flipped = np.fft.ifft2(phase_flip_spectrum(observed_fft, ctf_grid)).real
wiener_corrected = np.fft.ifft2(
    wiener_ctf_spectrum(observed_fft, ctf_grid, regularization=0.03)
).real

nonzero = np.abs(ctf_grid) > 1e-8
assert np.allclose(
    np.abs(phase_flip_spectrum(observed_fft, ctf_grid)[nonzero]),
    np.abs(observed_fft[nonzero]),
)
assert np.isfinite(wiener_corrected).all()

imshow_all(true_projection, observed, phase_flipped, wiener_corrected,
           titles=["true projection", "CTF + noise", "phase flipping", "Wiener-style"])


# %% [markdown]
# ```{admonition} 主張範圍
# :class: caution
#
# 這個 CTF demo 假設 isotropic defocus，並用單一 B factor 代表 envelope；沒有完整建模 astigmatism、anisotropic magnification、beam tilt、higher-order aberrations、DQE 或 colored noise。它用來區分 correction operators，不是 production CTF estimation pipeline。
# ```
#
# ## 理解檢查
#
# 1. 對 shape `(N, M)` 的影像，為什麼 `k_y` 必須除以 `N`、`k_x` 必須除以 `M`？
# 2. `H[0,0]` 與影像平均值差一個什麼因子？
# 3. 加 Hann window 與 zero-padding 分別改變頻譜的哪個部分？
# 4. 為什麼 `G/(H+eps)` 不能稱為 exact inverse filter？
# 5. Phase flipping 是否會復原 CTF 壓低的振幅？Wiener regularization 又付出什麼代價？
