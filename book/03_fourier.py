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
# 傅立葉轉換把「影像裡哪裡亮」改寫成「各空間頻率有多少振幅與相位」。這是一般影像濾波、復原與取樣的共通語言，也讓 cryo-EM 的投影、對比傳遞函數（contrast transfer function, CTF）、平移與重建變得比較直接。本章先建立通用頻域工具，再把它們用於 cryo-EM。索引、正規化與邊界慣例必須全章一致。頻域濾波與影像復原的更完整討論可參考 {cite}`szeliski2022,forsyth2012,penczek2010restoration`。
#
# ```{admonition} 學習目標
# :class: important
#
# - 寫出矩形影像的 2D DFT／inverse DFT，核對 DC、Parseval 與 Hermitian symmetry。
# - 解釋頻率軸、振幅、相位、有限視窗與 zero-padding。
# - 比較空間域與頻域卷積，並設計低通、高通、帶通與 notch filters。
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
# FFT 輸出的 DC 位於 `[0,0]`。`fftshift()` 只重新排列頻率係數，把 DC 移到畫面中央，方便觀察低頻到高頻的空間排列。

# %%
camera = ski.util.img_as_float(ski.data.camera())
camera_fft = np.fft.fft2(camera)
log_magnitude = np.log1p(np.abs(np.fft.fftshift(camera_fft)))
phase = np.angle(np.fft.fftshift(camera_fft))
imshow_all(camera, log_magnitude, phase,
           titles=["image", "log magnitude", "phase"])


# %% [markdown]
# ## 二維 Fourier basis、頻率軸與方向
#
# 二維 DFT 的每個 basis image 是
#
# $$b_{k_y,k_x}[y,x]=e^{i2\pi(k_y y/N+k_x x/M)}.$$
#
# 實數的 cosine 需要一對共軛頻率來表示。`fftfreq` 會依取樣間隔產生正確頻率軸；對 pixel size 為 $p$ 的影像，單位為 cycles per physical length，Nyquist frequency 為 $1/(2p)$。條紋的空間方向與頻率向量 $(f_y,f_x)$ 垂直；頻譜亮點指出頻率向量的方向，與條紋延伸方向相差 $90^\circ$。

# %%
n_rows, n_cols = 96, 144
pixel_size = 1.5
fy0, fx0 = 8 / (n_rows * pixel_size), 15 / (n_cols * pixel_size)
y_A = np.arange(n_rows)[:, None] * pixel_size
x_A = np.arange(n_cols)[None, :] * pixel_size
oriented_grating = np.cos(2 * np.pi * (fy0 * y_A + fx0 * x_A) + 0.4)
grating_fft = np.fft.fftshift(np.fft.fft2(oriented_grating))
fy_axis = np.fft.fftshift(np.fft.fftfreq(n_rows, d=pixel_size))
fx_axis = np.fft.fftshift(np.fft.fftfreq(n_cols, d=pixel_size))

peak_indices = np.argpartition(np.abs(grating_fft).ravel(), -2)[-2:]
peak_coordinates = np.array(np.unravel_index(peak_indices, grating_fft.shape)).T
measured_peaks = {(round(fy_axis[i], 8), round(fx_axis[j], 8))
                  for i, j in peak_coordinates}
expected_peaks = {(round(fy0, 8), round(fx0, 8)),
                  (round(-fy0, 8), round(-fx0, 8))}
assert measured_peaks == expected_peaks

fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
axes[0].imshow(oriented_grating)
axes[0].set_title("oriented cosine basis")
axes[0].axis("off")
axes[1].imshow(np.log1p(np.abs(grating_fft)), extent=[fx_axis[0], fx_axis[-1],
                                                     fy_axis[-1], fy_axis[0]])
axes[1].set_xlabel("$f_x$ [cycles/Å]")
axes[1].set_ylabel("$f_y$ [cycles/Å]")
axes[1].set_title("conjugate frequency pair")
fig.tight_layout()


# %% [markdown]
# ## Magnitude、phase 與 translation theorem
#
# 頻譜的 magnitude $|H|$ 表示每個 basis 的強度，phase $\arg H$ 決定這些 basis 在空間中如何對齊。將影像平移 $(\Delta y,\Delta x)$ 時，magnitude 不變，但頻譜會乘上線性 phase ramp：
#
# $$h[y-\Delta y,x-\Delta x]\ \longleftrightarrow\
# H[k_y,k_x]e^{-i2\pi(k_y\Delta y/N+k_x\Delta x/M)}.$$
#
# 所以「magnitude 比 phase 重要」或「phase 比 magnitude 重要」都過度簡化：它們承載不同資訊，完整影像需要兩者。在 cryo-EM 中，particle translation 正是以這個 phase ramp 來處理。

# %%
translation = (13, -9)
shifted_camera = np.roll(camera, translation, axis=(0, 1))
fy = np.fft.fftfreq(camera.shape[0])[:, None]
fx = np.fft.fftfreq(camera.shape[1])[None, :]
phase_ramp = np.exp(-2j * np.pi * (fy * translation[0] + fx * translation[1]))
predicted_shift_fft = camera_fft * phase_ramp

assert np.allclose(np.abs(np.fft.fft2(shifted_camera)), np.abs(camera_fft))
assert np.allclose(np.fft.fft2(shifted_camera), predicted_shift_fft, atol=1e-9)
imshow_all(camera, shifted_camera,
           titles=["original", f"circular shift {translation}"])


# %% [markdown]
# ## 有限觀測與 spectral leakage
#
# DFT 把有限長度訊號視為週期重複。若觀測窗內含有非整數個週期，首尾接合會產生不連續，能量便散到鄰近 frequency bins，稱為 spectral leakage。加窗能降低遠端 leakage，但會拓寬主峰；zero-padding 會增加頻譜的取樣點，真正的頻率解析力仍由觀測長度決定。

# %%
n = 128
sample_index = np.arange(n)
non_bin_tone = np.sin(2 * np.pi * 10.35 * sample_index / n)
rect_spectrum = np.abs(np.fft.rfft(non_bin_tone))
hann_spectrum = np.abs(np.fft.rfft(non_bin_tone * np.hanning(n)))
zero_padded_spectrum = np.abs(np.fft.rfft(non_bin_tone, n=8 * n))

fig, ax = plt.subplots(figsize=(8, 3.5))
ax.semilogy(rect_spectrum / rect_spectrum.max(), label="rectangular window")
ax.semilogy(hann_spectrum / hann_spectrum.max(), label="Hann window")
ax.set_xlabel("frequency bin")
ax.set_ylabel("normalized magnitude")
ax.legend()
fig.tight_layout()

assert zero_padded_spectrum.size == 4 * n + 1
assert np.isclose(np.fft.rfftfreq(n)[1], 1 / n)
assert np.isclose(np.fft.rfftfreq(8 * n)[1], 1 / (8 * n))

# %% [markdown]
# Hann window 以較寬的 main lobe 換取較低的 side lobes，因此適合減少遠離主峰的 leakage；頻率解析度則受變寬的主瓣限制。將 128 點補零成 1024 點，只是在同一個有限視窗的 discrete-time Fourier transform 上取得更密的樣本；真正分開相近頻率的能力仍由原始觀測長度與 window 決定 {cite}`szeliski2022`。


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
# ### 二維空間卷積與 FFT 卷積
#
# 卷積定理在二維一樣成立：$g=h*f$ 對應 $G=HF$。空間域直接卷積和 FFT 卷積可以得到相同的 linear-convolution 結果，但必須具有相同的 padding、裁切與 boundary 定義。小 kernel 的直接法常常較快；大 kernel 或重複使用同一 transfer function 時，FFT 法才逐漸有利。切換點受影像大小、kernel 大小、硬體與實作影響，不存在「FFT 永遠較快」的通則。

# %%
small_image = ski.util.img_as_float(ski.data.camera())[::8, ::8]
kernel_2d = np.outer(signal.windows.gaussian(9, 1.5),
                     signal.windows.gaussian(9, 1.5))
kernel_2d /= kernel_2d.sum()
spatial_full = signal.convolve2d(small_image, kernel_2d, mode="full")
fft_full = signal.fftconvolve(small_image, kernel_2d, mode="full")
assert np.allclose(spatial_full, fft_full, atol=1e-12)


# %% [markdown]
# ## 頻域濾波：從 transfer function 讀出空間效應
#
# 頻域濾波器直接乘在影像頻譜上。Low-pass 保留緩慢變化，high-pass 保留快速變化，band-pass 只保留一段尺度，notch filter 則壓低特定方向與頻率。這些名稱描述 transfer function，可用不同方法實作。例如 ideal low-pass 的截止邊界不連續，其 impulse response 具有長距離振盪，因而在強邊緣附近產生 ringing；Gaussian 轉換平滑，空間域的 impulse response 保持正值。Butterworth 的 order 則提供兩者之間的過渡 {cite}`szeliski2022,forsyth2012`。

# %%
def frequency_radius(shape):
    """Return an unshifted radial-frequency grid in cycles/pixel."""
    fy = np.fft.fftfreq(shape[0])[:, None]
    fx = np.fft.fftfreq(shape[1])[None, :]
    return np.sqrt(fy**2 + fx**2)


def lowpass_transfer(shape, cutoff, kind="gaussian", order=3):
    """Construct ideal, Gaussian, or Butterworth radial low-pass transfer."""
    if not 0 < cutoff <= 0.5:
        raise ValueError("cutoff must be in (0, 0.5]")
    radius = frequency_radius(shape)
    if kind == "ideal":
        return (radius <= cutoff).astype(float)
    if kind == "gaussian":
        return np.exp(-0.5 * (radius / cutoff) ** 2)
    if kind == "butterworth":
        return 1.0 / (1.0 + (radius / cutoff) ** (2 * order))
    raise ValueError("kind must be ideal, gaussian, or butterworth")


def apply_transfer(image, transfer_function):
    """Apply a same-shape transfer function under circular boundary conditions."""
    if image.shape != transfer_function.shape:
        raise ValueError("image and transfer_function must have the same shape")
    return np.fft.ifft2(np.fft.fft2(image) * transfer_function).real


filter_image = ski.transform.resize(camera, (192, 192), anti_aliasing=True)
lowpasses = {kind: lowpass_transfer(filter_image.shape, 0.08, kind, order=4)
             for kind in ("ideal", "gaussian", "butterworth")}
filtered = {kind: apply_transfer(filter_image, transfer_function)
            for kind, transfer_function in lowpasses.items()}

fig, axes = plt.subplots(3, 3, figsize=(10, 9))
for row, kind in enumerate(("ideal", "gaussian", "butterworth")):
    transfer_function = lowpasses[kind]
    impulse_response = np.fft.fftshift(np.fft.ifft2(transfer_function).real)
    axes[row, 0].imshow(np.fft.fftshift(transfer_function), vmin=0, vmax=1)
    axes[row, 0].set_title(f"{kind}: transfer")
    center = tuple(size // 2 for size in impulse_response.shape)
    crop = impulse_response[center[0]-16:center[0]+17,
                            center[1]-16:center[1]+17]
    axes[row, 1].imshow(crop)
    axes[row, 1].set_title("impulse response")
    axes[row, 2].imshow(filtered[kind])
    axes[row, 2].set_title("filtered image")
for ax in axes.ravel():
    ax.axis("off")
fig.tight_layout()

ideal_impulse = np.fft.fftshift(np.fft.ifft2(lowpasses["ideal"]).real)
gaussian_impulse = np.fft.fftshift(np.fft.ifft2(lowpasses["gaussian"]).real)
assert ideal_impulse.min() < 0
assert gaussian_impulse.min() > -1e-10


# %% [markdown]
# ### Low-pass、high-pass、band-pass 與 notch 是可組合的工具
#
# High-pass 可以用 $1-L$ 從 low-pass $L$ 建立；band-pass 可由兩個不同截止頻率的 low-pass 相減。Notch filter 必須成對放在 $(f_y,f_x)$ 與 $(-f_y,-f_x)$，才會對應實數輸出。頻率濾波會一起壓低同頻率的物件細節與干擾；notch 移除的資訊也會留在結果中。

# %%
radius = frequency_radius(filter_image.shape)
butter_low = lowpass_transfer(filter_image.shape, 0.08, "butterworth", order=4)
butter_high = 1.0 - butter_low
butter_band = (
    lowpass_transfer(filter_image.shape, 0.16, "butterworth", order=4)
    - lowpass_transfer(filter_image.shape, 0.04, "butterworth", order=4)
)

fy_grid = np.fft.fftfreq(filter_image.shape[0])[:, None]
fx_grid = np.fft.fftfreq(filter_image.shape[1])[None, :]
notch_center = (0.12, 0.18)
notch_width = 0.015
notch = np.ones(filter_image.shape)
for sign_value in (-1, 1):
    distance_squared = ((fy_grid - sign_value * notch_center[0]) ** 2
                        + (fx_grid - sign_value * notch_center[1]) ** 2)
    notch *= 1.0 - np.exp(-distance_squared / (2 * notch_width**2))

transfer_examples = [butter_low, butter_high, butter_band, notch]
transfer_titles = ["low-pass", "high-pass", "band-pass", "paired notch"]
fig, axes = plt.subplots(2, 4, figsize=(13, 6))
for column, (transfer_function, title) in enumerate(zip(transfer_examples, transfer_titles)):
    axes[0, column].imshow(np.fft.fftshift(transfer_function), vmin=0, vmax=1)
    axes[0, column].set_title(title)
    axes[1, column].imshow(apply_transfer(filter_image, transfer_function))
    axes[1, column].set_title("output")
for ax in axes.ravel():
    ax.axis("off")
fig.tight_layout()

notch_complex_output = np.fft.ifft2(np.fft.fft2(filter_image) * notch)
notch_output = notch_complex_output.real
assert np.max(np.abs(notch_complex_output.imag)) < 1e-12


# %% [markdown]
# ## Inverse filtering 是病態問題
#
# 若 $G=HF+N$，在 $N=0$ 且 $H\ne0$ 時可寫成 $F=G/H$。當 $H$ 接近零，微小的 $N$ 也會被大幅放大。程式裡的 `G / (H + eps)` 屬於帶偏差的近似；`eps` 改變了 transfer function，也可能造成相位偏移。
#
# 較透明的示範是 truncated inverse：只在 $|H|\ge\tau$ 的頻率做除法，其餘設為零。$\tau$ 是明確的 regularization parameter，輸出屬於帶正則化的近似復原。

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
# Point-spread function（PSF）描述系統對單一點的空間域響應；它的 Fourier transform 是 optical transfer function（OTF）。程式中 PSF 的中心通常畫在陣列中央，但 FFT 的 origin 在 `[0,0]`，因此 `psf_to_otf` 必須先 padding，再把 PSF 中心移到 origin。少了這步會多出一個 phase ramp，造成復原影像平移。影像復原的前向模型、邊界與正則化必須配套解讀 {cite}`penczek2010restoration`。


# %% [markdown]
# ## Wiener regularization
#
# 在簡化的白雜訊／平穩訊號模型下，常見的 Wiener-style estimator 可寫成
#
# $$\widehat F=\frac{H^*}{|H|^2+K}G.$$
#
# $K$ 代表 noise-to-signal power 的近似；$K$ 越大，在 transfer 弱的頻率給予的增益越小。這是 bias–variance trade-off，零點處未被量到的資訊仍然缺失。

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
# ### Regularization 強度決定 bias–variance trade-off
#
# $K$ 太小時，弱 transfer 處的雜訊被大幅放大；$K$ 太大時，復原過度平滑。下例因為有已知參考影像，可以計算均方誤差（mean squared error, MSE）與峰值訊雜比（peak signal-to-noise ratio, PSNR）。真實資料缺少未模糊的已知真值，選擇 $K$ 時要配合雜訊模型與獨立驗證。PSNR 衡量像素誤差，科學結構是否正確還需要其他證據。

# %%
regularization_values = np.array([1e-6, 1e-4, 2e-3, 5e-2])
restoration_candidates = [
    wiener_frequency(noisy_blurred, transfer, value)
    for value in regularization_values
]
mse_values = np.array([
    np.mean((candidate - original) ** 2)
    for candidate in restoration_candidates
])
psnr_values = np.array([
    ski.metrics.peak_signal_noise_ratio(original, candidate, data_range=1.0)
    for candidate in restoration_candidates
])
observed_mse = np.mean((noisy_blurred - original) ** 2)

assert np.isfinite(mse_values).all()
assert np.allclose(psnr_values, 10 * np.log10(1.0 / mse_values))
assert mse_values.min() < observed_mse

fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
axes[0].semilogx(regularization_values, mse_values, "o-")
axes[0].set_xlabel("regularization $K$")
axes[0].set_ylabel("MSE")
axes[0].set_title("lower is better")
axes[1].semilogx(regularization_values, psnr_values, "o-")
axes[1].set_xlabel("regularization $K$")
axes[1].set_ylabel("PSNR [dB]")
axes[1].set_title("higher is better")
fig.tight_layout()


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
    """Isotropic CTF；envelope 以可選 B factor 近似，astigmatism 留待二維模型。"""
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
# ```{admonition} 這段示範包含的條件
# :class: caution
#
# 這個 CTF 示範採用 isotropic defocus，並用單一 B factor 近似 envelope，適合比較 phase flipping 與 Wiener-style correction 的運算差異。完整的 CTF estimation 還要處理 astigmatism、anisotropic magnification、beam tilt、higher-order aberrations、DQE 與 colored noise。
# ```
#
# ## 理解檢查
#
# 1. 對 shape `(N, M)` 的影像，二維 DFT 的正規化因子如何決定？`H[0,0]` 又和影像平均值差多少？
#
#    ```{dropdown} 參考答案
#    令影像為 $h[n,m]$，其中 $0\le n<N$、$0\le m<M$。正向 DFT 的指數項是
#
#    $$
#    e^{-i2\pi(kn/N+\ell m/M)},
#    $$
#
#    因此垂直頻率索引 $k$ 要除以列數 $N$，水平頻率索引 $\ell$ 要除以欄數 $M$。若交換兩個分母，矩形影像的基底頻率便會算錯。採用 NumPy 預設慣例時，正向轉換不乘正規化常數，反向轉換乘 $1/(NM)$。零頻係數為
#
#    $$
#    H[0,0]=\sum_{n=0}^{N-1}\sum_{m=0}^{M-1}h[n,m]=NM\,\bar h.
#    $$
#
#    例如 $2\times3$ 影像的平均值是 4，則 `H[0,0]` 是 $6\times4=24$。常見誤解是把 DC coefficient 直接當成平均值；只有使用 `norm="ortho"` 或自行正規化時，比例才會不同。回答這類問題時要先寫明採用的 DFT 慣例。
#    ```
#
# 2. 影像平移後，Fourier magnitude 與 phase 會如何改變？
#
#    ```{dropdown} 參考答案
#    若 $g(y,x)=f(y-\Delta y,x-\Delta x)$，translation theorem 給出
#
#    $$
#    G(f_y,f_x)=F(f_y,f_x)e^{-i2\pi(f_y\Delta y+f_x\Delta x)}.
#    $$
#
#    右側新增的因子絕對值為 1，所以 $|G|=|F|$；phase 則多出一個隨頻率線性變化的斜坡。以向右平移 3 pixels 為例，水平頻率 $f_x$ 的相位改變 $-2\pi(3f_x)$。這項性質說明 magnitude spectrum 無法單獨決定影像位置。離散影像若採 circular shift，公式可直接套用；實際裁切、補零或邊界截斷會改變影像內容，此時 magnitude 也可能跟著改變。
#    ```
#
# 3. Hann window 與 zero-padding 各自改變了什麼？為什麼 ideal low-pass 容易產生 ringing？
#
#    ```{dropdown} 參考答案
#    Hann window 先讓訊號兩端平順降到接近零，可減輕週期接合處的不連續，因此遠離主峰的 spectral leakage 會下降；代價是主瓣變寬，鄰近頻率更難分開。Zero-padding 只在既有 DFT 取樣點之間加入更密的頻率取樣，畫出的頻譜較平滑，觀測時間與真正的頻率解析力都沒有增加。
#
#    Ideal low-pass 在截止頻率處突然從 1 跳到 0。其 inverse Fourier transform 含有延伸很遠的振盪旁瓣，強邊緣經過卷積後便出現正負交替的波紋。Gaussian low-pass 的 transfer function 平滑，impulse response 也保持 Gaussian 形狀，因此同類 ringing 較弱。常見誤解是把 zero-padding 後較密的頻譜當成新增解析度，或把 ringing 當成影像原有的細節。
#    ```
#
# 4. 已知 $G=HF+N$ 時，為何 `G/(H+eps)` 仍會產生偏差？正則化如何在雜訊與細節之間取捨？
#
#    ```{dropdown} 參考答案
#    理想反濾波 $F=G/H$ 需要 $N=0$ 且每個頻率的 $H\ne0$。加入 `eps` 後，對無雜訊資料得到
#
#    $$
#    \widehat F=\frac{HF}{H+\varepsilon}
#    =F\frac{H}{H+\varepsilon},
#    $$
#
#    乘上的比例通常偏離 1，所以估計帶有偏差；若 $H$ 為複數，直接加實數 `eps` 還可能改變相位。Wiener-style filter 常寫成
#
#    $$
#    \widehat F=\frac{H^*}{|H|^2+K}G.
#    $$
#
#    在 $|H|$ 很小的頻率，分母中的 $K$ 可抑制雜訊放大；$K$ 太大會抹去可恢復的細節，太小則接近不穩定的 inverse filter。例如 $H=0.01$、$K=0.01$ 時，Wiener gain 約為 $0.01/(0.0001+0.01)\approx0.99$，遠低於 inverse gain 100。`eps` 或 $K$ 都需要配合雜訊水準與任務評估，沒有一個值能通用於所有影像。
#    ```
#
# 5. Phase flipping 與 Wiener-style CTF correction 分別改變哪些 Fourier 資訊？
#
#    ```{dropdown} 參考答案
#    Phase flipping 把觀測頻譜乘上 $\operatorname{sign}(\mathrm{CTF})$。在 $\mathrm{CTF}\ne0$ 的頻率，乘數的絕對值是 1，所以觀測振幅保持不變；CTF 為負時，相位會翻回來。依本章實作，零點的乘數設為 0。CTF 已壓低的振幅與零點處遺失的資訊仍維持原狀。
#
#    Wiener-style correction 使用
#
#    $$
#    \frac{\mathrm{CTF}^*}{|\mathrm{CTF}|^2+K},
#    $$
#
#    同時調整相位與振幅。$K$ 防止 CTF 接近零時增益暴增，也會讓估計產生平滑偏差。舉例來說，CTF 為 $-0.5$、$K=0.05$ 時，gain 為 $-0.5/(0.25+0.05)\approx-1.67$；phase flipping 的 gain 則是 $-1$。兩者都無法憑單張影像補回 CTF 零點，實務上會利用不同離焦影像互補頻率資訊，並配合更完整的 CTF 與雜訊模型。
#    ```
