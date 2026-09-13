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
# # 傅立葉轉換、影像金字塔與特徵偵測
#
# 同一張影像可以用像素的位置與亮度描述，也可以用不同頻率的波描述。本章先從 10 Hz 的正弦波認識 DFT，再用 Gengar、Lucario、Pikachu 和 Dragonite 觀察影像頻譜，操作濾波、復原與去除週期干擾。最後回到影像空間，建立金字塔，找出角點與斑點。
#
# ## 從一個 10 Hz 正弦波開始


# %% [markdown]
# 先執行下列設定，再依章節順序操作。互動圖可直接在網頁調整，Python 圖則列出完整的比較結果。

# %% tags=["hide-input"]
from pathlib import Path
import sys
import math
import timeit

import matplotlib.pyplot as plt
import numpy as np
import skimage as ski
from scipy import fft as fp, signal
from scipy.signal import convolve2d as conv2
from IPython.display import Image, display

BOOK = Path('book') if Path('book').is_dir() else Path('.')
sys.path.insert(0, str(BOOK.resolve()))
from _support import image_path, show_images, lab

plt.rcParams['image.cmap'] = 'gray'


def read_rgb(name):
    image = ski.io.imread(image_path(name))
    if image.ndim == 2:
        return ski.color.gray2rgb(ski.util.img_as_float(image))
    if image.shape[-1] == 4:
        return ski.color.rgba2rgb(image)
    return ski.util.img_as_float(image)


def read_gray(name):
    return ski.color.rgb2gray(read_rgb(name))


def real_ifft2(coefficients):
    """確認共軛對稱與虛部殘差，再取實部。"""
    rows = (-np.arange(coefficients.shape[0])) % coefficients.shape[0]
    cols = (-np.arange(coefficients.shape[1])) % coefficients.shape[1]
    partner = np.conj(coefficients[np.ix_(rows, cols)])
    scale = max(1.0, float(np.max(np.abs(coefficients))))
    assert np.allclose(coefficients, partner, atol=1e-10 * scale, rtol=1e-10)
    result = fp.ifft2(coefficients)
    assert np.max(np.abs(result.imag)) < 1e-9 * max(1.0, np.max(np.abs(result.real)))
    return result.real


def frequency_indices(shape):
    """FFT儲存順序下各軸的有號整數頻率索引。"""
    return (np.rint(fp.fftfreq(shape[0]) * shape[0]).astype(int)[:, None],
            np.rint(fp.fftfreq(shape[1]) * shape[1]).astype(int)[None, :])


def square_lowpass(shape, cutoff):
    ky, kx = frequency_indices(shape)
    return (np.abs(ky) <= cutoff) & (np.abs(kx) <= cutoff)


def gaussian_psf(shape, std):
    psf = np.outer(signal.windows.gaussian(shape[0], std),
                   signal.windows.gaussian(shape[1], std))
    return psf / psf.sum()


def coefficient_basis(n, kx, ky, amplitude=1.0, wave_type='cosine'):
    """amplitude是單一FFT係數大小，非空間波振幅。"""
    spectrum = np.zeros((n, n), dtype=complex)
    pos, neg = (ky % n, kx % n), ((-ky) % n, (-kx) % n)
    if pos == neg:
        spectrum[pos] = amplitude if wave_type == 'cosine' else 0.0
    else:
        value = amplitude if wave_type == 'cosine' else -1j * amplitude
        spectrum[pos], spectrum[neg] = value, np.conj(value)
    return spectrum, real_ifft2(spectrum)


def draw_spectrum_panels(image, title):
    coefficients = np.fft.fft2(image)
    reconstructed = np.fft.ifft2(coefficients).real
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    panels = [image, reconstructed,
              np.log1p(np.abs(np.fft.fftshift(coefficients))),
              np.angle(np.fft.fftshift(coefficients))]
    for ax, shown, name in zip(axes.flat, panels,
                              ['Original', 'IFFT reconstruction', 'Log magnitude', 'Phase']):
        ax.imshow(shown, cmap='twilight' if name == 'Phase' else 'gray')
        ax.set_title(f'{title}: {name}')
        ax.axis('off')
    fig.tight_layout()
    plt.show()
    assert np.allclose(reconstructed, image, atol=1e-12)
    return coefficients


def montage(images, titles, signed=False, ncols=4, figsize=(12, 18)):
    nrows = math.ceil(len(images) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    for ax, image, title in zip(axes.flat, images, titles):
        options = {}
        if signed:
            limit = max(1e-12, float(np.max(np.abs(image))))
            options = dict(cmap='RdBu_r', vmin=-limit, vmax=limit)
        ax.imshow(image, **options)
        ax.set_title(title)
        ax.axis('off')
    for ax in list(axes.flat)[len(images):]:
        ax.axis('off')
    fig.tight_layout()
    plt.show()


# %% [markdown]
# 以每秒 100 次的頻率取樣兩秒，得到 200 個取樣點。訊號 $h(t)=\sin(10\times2\pi t)$ 每秒重複 10 次；先看它在時間軸上的樣子。

# %%
f = 10  # Frequency, in cycles per second, or Hertz
f_s = 100  # Sampling rate, or number of measurements per second

t = np.linspace(0, 2, 2*f_s, endpoint=False)
x = np.sin(f*2*np.pi*t)

fig, ax = plt.subplots()
ax.plot(t, x)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Signal amplitude');

# %% [markdown]
# 這是一個頻率為 10 Hz 的正弦波，其在時間域（time domain）中呈現出規律的波動。

# %% [markdown]
# DFT 將這 200 個數改寫成不同頻率的複數係數。某個係數的幅度大，表示對應的波在訊號中占有較大分量。對實數正弦波，正、負頻率需要成對出現。

# %% [markdown]
# ### 用波表示訊號
#
# 傅立葉級數以離散頻率描述週期函數；傅立葉轉換則以連續頻率描述適當的非週期訊號。對電腦裡長度為 $N$ 的陣列，DFT 是一個可逆的線性座標轉換：使用 $N$ 個複指數基底，就能表示這 $N$ 個取樣值。這不表示原本的連續訊號已被完整量到，取樣與觀測時間仍會限制可恢復的資訊。

# %% [markdown]
# ```{dropdown} 從連續積分到 DFT 的取樣格點
# 以 $e^{i\omega t}$ 為基底，連續傅立葉轉換採負指數：
#
# $$H(\omega)=\int_{-\infty}^{\infty}h(t)e^{-i\omega t}\,dt.$$
#
# 有限時間 $T=N\Delta t$ 的觀測，可在取樣夠密且適合數值積分時近似為
#
# $$H_T(\omega)\approx\Delta t\sum_{n=0}^{N-1}h[n]e^{-i\omega n\Delta t}.$$
#
# DFT 選用 $\omega_k=2\pi k/(N\Delta t)$ 這組正交格點，並省去共同的 $\Delta t$ 縮放，得到
#
# $$H[k]=\sum_{n=0}^{N-1}h[n]e^{-i2\pi kn/N},\qquad
# h[n]=\frac1N\sum_{k=0}^{N-1}H[k]e^{i2\pi kn/N}.$$
#
# 有限序列也可在其他頻率計算其 DTFT；DFT 的 $N$ 個格點已足以表示這個陣列。格點間距為 $1/T$；能否分辨鄰近真實頻率，還受觀測長度、窗函數與雜訊影響。
#
# 對實數 $h[n]$，實部是與 cosine 的內積，虛部是與 sine 內積的**負值**。係數幅度與相位分別是 $|H[k]|$、$\arg H[k]$；正負頻率滿足 $H[-k]=H[k]^*$，索引按 $N$ 取模。
# ```

# %%
X = fp.fft(x)
freqs = fp.fftfreq(len(x), d=1/f_s) # d: Sample spacing (inverse of the sampling rate)
print(freqs, len(x), len(X))
# The first component is np.mean(x) * N
fig, ax = plt.subplots()

ax.stem(freqs, np.abs(X))
ax.set_xlabel('Frequency in Hertz [Hz]')
ax.set_ylabel('Frequency Domain (Spectrum) Magnitude')
ax.set_xlim(-f_s / 2, f_s / 2)
ax.set_ylim(-5, 110);

# %% [markdown]
# 頻譜在 +10 Hz 與 −10 Hz 各有一個峰。理想的數學結果只有這兩個非零係數；數值 FFT 還會留下浮點誤差尺度的小值。因為 forward DFT 沒有除以 $N$，振幅為 1 的 sine 在兩個峰上的係數幅度各是 $N/2=100$。

# %% [markdown]
# > 有關一維討論，請參考以下資源：
#
# > [https://github.com/elegant-scipy/notebooks/blob/master/notebooks/ch4.ipynb](https://github.com/elegant-scipy/notebooks/blob/master/notebooks/ch4.ipynb)
#
# > [https://scipy-lectures.org/intro/scipy.html?highlight=fft#fast-fourier-transforms-scipy-fftpack](https://scipy-lectures.org/intro/scipy.html?highlight=fft#fast-fourier-transforms-scipy-fftpack)
#
# > [https://dsp.stackexchange.com/questions/23994/meaning-of-real-and-imaginary-part-of-fourier-transform-of-a-signal](https://dsp.stackexchange.com/questions/23994/meaning-of-real-and-imaginary-part-of-fourier-transform-of-a-signal)

# %% [markdown]
# 本章使用 [`scipy.fft`](https://docs.scipy.org/doc/scipy/reference/fft.html)：`fft`／`ifft` 處理一維，`fft2`／`ifft2` 處理二維，`fftn`／`ifftn` 處理多維。`fftfreq` 建立頻率軸，`fftshift` 與 `ifftshift` 在顯示順序和 FFT 儲存順序之間轉換。NumPy 的同名 FFT 函式採用相同的預設正規化。

# %% [markdown]
# ```{dropdown} 補充：實數 FFT、餘弦轉換與窗函數
# 實數訊號的頻譜具有共軛對稱性，`rfft` 只儲存非負頻率部分；反轉換使用 `irfft`，並指定原訊號長度。DCT／DST 分別以餘弦／正弦基底表示資料，反轉換是 `idct`／`idst`。JPEG 的轉換編碼使用 DCT，可把平滑區塊的主要變化集中在少數係數中。
#
# 有限長度的截取可能讓訊號兩端接合時產生跳躍，使能量散到其他頻率，稱為頻譜洩漏（spectral leakage）。`np.hanning`、`np.hamming`、`np.bartlett`、`np.blackman`、`np.kaiser` 等窗函數讓邊界逐漸衰減；代價是主瓣變寬，鄰近頻率較難分開。
#
# 原稿提到的 [FFTPACK](https://www.netlib.org/fftpack/) 是歷史實作背景；[FFTW](https://www.fftw.org/) 是另一套 FFT 函式庫。使用本章的 `scipy.fft` 時，請以其文件為準，勿把舊 `fftpack.rfft` 的儲存格式直接套用到新 API。
# ```

# %% [markdown]
# #### 頻率與其排序方式（Frequencies and their ordering）

# %% [markdown]
# FFT 的輸出先放零頻率，再放正頻率與負頻率。例如十個 1 組成的陣列沒有起伏，只有零頻率係數不為零。

# %% [markdown]
# (dft)=
#
# **DC 係數**是零頻率係數。依本章的正規化，
#
# $$H[0]=\sum_{n=0}^{N-1}h[n]=N\bar h.$$
#
# 它除以 $N$ 才是訊號平均值。

# %%
N = 10
print(fp.fft(np.ones(N)))

# %% [markdown]
# 接著用九個數觀察共軛對稱。對實數序列，正負頻率的實部相同、虛部符號相反；這裡的「對稱」指頻率 $k$ 與 $-k$，不是直接把未排序陣列左右對折。

# %%
x = np.array([1, 5, 12, 7, 3, 0, 4, 3, 2])
X = fp.fft(x)

with np.printoptions(precision=2):
    print("Real part:     ", X.real)
    print("Imaginary part:", X.imag)

# %% [markdown]
# `fftfreq` 函數可以告訴我們所對應的頻率為何，也就是每個頻譜成分實際代表的頻率:

# %%
print(fp.fftfreq(len(x), d=1))

# %%
# Sort the Fourier coefficients by frequency
freqs = fp.fftfreq(len(x), d=1)
sorted_indices = np.argsort(freqs)
sorted_X = X[sorted_indices]

# Display the sorted real and imaginary parts
with np.printoptions(precision=2):
    print("Real part sorted:     ", sorted_X.real)
    print("Imaginary part sorted:", sorted_X.imag)

# %%
# Use fftshift to Sort the Fourier coefficients by frequency
sorted_X = fp.fftshift(X)

# Display the sorted real and imaginary parts
with np.printoptions(precision=2):
    print("Real part sorted:     ", sorted_X.real)
    print("Imaginary part sorted:", sorted_X.imag)

# %% [markdown]
# 請參考這個連結了解更多細節：
# [https://numpy.org/doc/stable/reference/generated/numpy.fft.fftfreq.html](https://numpy.org/doc/stable/reference/generated/numpy.fft.fftfreq.html)

# %% [markdown]
# ### 為什麼我們需要離散傅立葉轉換（DFT）？

# %% [markdown]
# 頻域表示有兩個直接用途。它能顯示影像有哪些尺度與方向的變化，也能把卷積改寫成頻譜的逐元素相乘。許多自然影像的低頻能量較強，但重要邊緣、紋理或小物件也可能出現在高頻；係數較小的成分仍可能對辨認特徵有用。

# %% [markdown]
# ## 二維 DFT：影像也是波的組合
#
# 影像的每一個像素值，可由不同方向、不同頻率的二維波加總而成。二維係數仍是複數，保留了幅度與相位。以下先操作單一基底，再看真實圖像的頻譜。

# %% [markdown]
# 對灰階影像使用 `fft2()` 與 `ifft2()`。它們等價於沿每一軸各做一次一維轉換；兩個軸的長度不必相同。

# %% [markdown]
# 二維的連續與離散轉換可對照下列三張原始示意圖。對 $M$ 列、$N$ 欄的影像，各軸長度分別寫成
#
# $$H[k_y,k_x]=\sum_{y=0}^{M-1}\sum_{x=0}^{N-1}
# h[y,x]e^{-2\pi i(k_y y/M+k_x x/N)}.$$
#
# 反轉換使用正指數，並除以 $MN$。影像的 DC 係數為像素總和。
#
# ```{dropdown} 將各個二維基底加回像素值
# 逆轉換為
#
# $$h[y,x]=\frac1{MN}\sum_{k_y=0}^{M-1}\sum_{k_x=0}^{N-1}
# H[k_y,k_x]e^{2\pi i(k_y y/M+k_x x/N)}.$$
#
# 每個 $H[k_y,k_x]$ 指定一個基底的係數。改變 $k_x,k_y$ 會改變波的頻率向量，條紋延伸方向與這個向量垂直。真實影像的頻譜滿足 $H[-k_y,-k_x]=H[k_y,k_x]^*$，所以一對共軛基底加總後的虛部抵消。
# ```
#
# 下面操作 128×128 的基底。控制項的 amplitude 表示**頻譜係數大小**：非自共軛頻率放入一對大小為 $a$ 的係數時，空間波的振幅是 $2a/128^2$。DC 只有一個獨立係數，cosine 產生常數 $a/128^2$，sine 在 DC 為零。

# %%
# 圖片來源： https://drive.google.com/uc?id=1ArBNhQHt4OHQX6fGRsQP80J9I8e8GpyI
display(Image(filename=str(image_path('ch03-cell28-1.png'))))

# %%
# 圖片來源： https://drive.google.com/uc?id=1AyaTNoPLbMAmpS3eruaZiA9JsiuQtQYk
display(Image(filename=str(image_path('ch03-cell28-2.png'))))

# %%
# 圖片來源： https://upload.wikimedia.org/wikipedia/commons/f/fa/2D_Fourier_Transform_and_Base_Images.png
display(Image(filename=str(image_path('ch03-cell28-3.png'))))

# %%
basis_spectrum, basis_image = coefficient_basis(128, 5, 2, 1.0, 'cosine')
show_images(fp.fftshift(np.abs(basis_spectrum)), basis_image,
            titles=['Two conjugate FFT coefficients', '128 x 128 cosine basis'])
for kind in ('cosine', 'sine'):
    _, dc_image = coefficient_basis(128, 0, 0, 1.0, kind)
    assert np.allclose(dc_image, 1 / 128**2 if kind == 'cosine' else 0)
lab('fourier_basis', n=128, kx=5, ky=2, amplitude=1.0)

# %% [markdown]
# ### Gengar：轉換後還能回到原圖嗎？
#
# 將 RGBA 圖片疊在白底上，轉成灰階，再計算 FFT 與 IFFT。比較重建誤差，確認改變表示方式並沒有刪掉資訊。

# %%
gengar = read_gray('gengar.png')
freq = fp.fft2(gengar)
gengar_reconstructed = real_ifft2(freq)
print('maximum reconstruction error:', np.max(np.abs(gengar_reconstructed - gengar)))
assert np.allclose(gengar, gengar_reconstructed, atol=1e-12)
show_images(gengar, gengar_reconstructed, titles=['Gengar', 'FFT then IFFT'])

# %% [markdown]
# #### 繪製頻率頻譜（Plotting the frequency spectrum）

# %% [markdown]
# 幅度圖顯示每個頻率係數的大小。DC 經常比其他係數大很多，直接用線性色階顯示時，較小的係數會看不清楚。

# %% [markdown]
# 使用 `log1p(abs(H))` 壓縮顯示範圍，再用 `fftshift` 把 DC 移到中心。取對數與位移只用於顯示；反轉換仍使用原本的複數係數。

# %%
show_images(np.log1p(np.abs(fp.fftshift(freq))), titles=['Gengar: centered log spectrum'])

# %% [markdown]
# ### Lucario 與 Pikachu：幅度、相位、重建
#
# 對另外兩張圖做相同操作，並排觀察原圖、重建、對數幅度與相位。這一段改用 `numpy.fft`，也可與 SciPy 的結果互相比對。

# %%
lucario = read_gray('lucario.png')
freq1 = draw_spectrum_panels(lucario, 'Lucario')

# %% [markdown]
# Lucario 的頻譜中，低頻通常較強；不同方向的紋理也會留下對應的頻率分布。相位圖雖然看起來雜亂，重建時仍需要它決定各頻率如何疊合。

# %% [markdown]
# 用同樣的色彩轉換與 FFT 流程處理 Pikachu，準備比較兩張圖的頻譜。

# %%
pikachu = read_gray('pikachu.png')
freq2 = draw_spectrum_panels(pikachu, 'Pikachu')

# %% [markdown]
# ### 交換幅度與相位
#
# 保留 Lucario 的幅度，換入 Pikachu 的相位，再反向交換一次。操作的是 $|H|$ 與 $\arg H$，不是實部與虛部。兩張圖必須有相同尺寸，才能逐係數組合。

# %%
assert lucario.shape == pikachu.shape, 'Phase exchange requires matching original image dimensions.'
combined = np.abs(freq1) * np.exp(1j * np.angle(freq2))
lucario_magnitude_pikachu_phase = real_ifft2(combined)
show_images(lucario_magnitude_pikachu_phase,
            titles=['Lucario magnitude + Pikachu phase'])

# %%
combined = np.abs(freq2) * np.exp(1j * np.angle(freq1))
pikachu_magnitude_lucario_phase = real_ifft2(combined)
show_images(pikachu_magnitude_lucario_phase,
            titles=['Pikachu magnitude + Lucario phase'])

# %% [markdown]
# 在這兩張圖上，交換後的輪廓明顯受到相位來源影響。相位控制各個波如何在空間中對齊；幅度仍保留尺度、方向與能量分布，也會改變重建結果。完整重建需要兩者。
#
# ```{dropdown} 為什麼整張圖平移後，幅度不變？
# 若將影像平移，頻譜只多出隨頻率變化的相位因子，其絕對值為 1。因此幅度對整體平移不敏感。這不表示幅度沒有空間資訊：例如水平與垂直條紋的幅度峰會落在不同方向。
# ```

# %% [markdown]
# (convolution)=
#
# ## 卷積定理：在頻域做 Gaussian 平滑

# %% [markdown]
# 讓我們從 **卷積定理（convolution theorem）** 開始，看看在頻率域中，卷積運算是如何變得更簡單的。

# %%
# 圖片來源： https://drive.google.com/uc?id=1B5QZeuQ3NiujL7GZjN69DMWR7NoBwuXq
display(Image(filename=str(image_path('ch03-cell46-1.png'))))

# %% [markdown]
# 下圖說明了在頻率域中進行濾波的基本步驟：

# %%
# 圖片來源： https://drive.google.com/uc?id=1BAhUknFnpeQeRdqO8LaviZksTv7l_nQ-
display(Image(filename=str(image_path('ch03-cell47-1.png'))))

# %% [markdown]
# ### Dragonite：直接相乘兩個頻譜
#
# 建立與影像同尺寸、標準差為 1 的 Gaussian 核，先讓權重總和為 1，再把核中心移到 FFT 原點。影像頻譜與核頻譜逐元素相乘，IFFT 後得到平滑影像。這裡沒有擴大 FFT 陣列，對應週期邊界的 circular convolution。

# %%
dragonite = read_gray('dragonite.png')
gauss_kernel = gaussian_psf(dragonite.shape, std=1)
freq = fp.fft2(dragonite)
freq_kernel = fp.fft2(fp.ifftshift(gauss_kernel))
convolved = freq * freq_kernel
im1 = real_ifft2(convolved)
assert np.isclose(im1.mean(), dragonite.mean())
fig, axes = plt.subplots(2, 3, figsize=(13, 8))
for ax, image, title in zip(axes.flat,
        [dragonite, gauss_kernel, im1,
         np.log1p(np.abs(fp.fftshift(freq))),
         np.log1p(np.abs(fp.fftshift(freq_kernel))),
         np.log1p(np.abs(fp.fftshift(convolved)))],
        ['Dragonite', 'Normalized Gaussian: std=1', 'Circular convolution',
         'Image spectrum', 'Kernel spectrum', 'Output spectrum']):
    ax.imshow(image)
    ax.set_title(title)
    ax.axis('off')
fig.tight_layout()
plt.show()

# %% [markdown]
# ### Lucario：用 `fftconvolve` 計算線性卷積
#
# 改用 11×11、標準差為 1 的 Gaussian 核。`fftconvolve` 會處理足夠的 padding，計算線性卷積，再依 `mode="same"` 裁成輸入尺寸。它包含 FFT、頻譜逐元素相乘與 IFFT，並非只做一次乘法。

# %%
im = lucario
gauss_kernel = gaussian_psf((11, 11), std=1)
im_blurred = signal.fftconvolve(im, gauss_kernel, mode='same')
show_images(im, gauss_kernel, im_blurred,
            titles=['Lucario', '11 x 11 Gaussian: std=1', 'Linear FFT convolution'])

# %% [markdown]
# 下列程式碼區塊示範了如何在卷積後繪製原始影像與模糊影像的頻譜（spectrum）：

# %%
show_images(np.log1p(np.abs(fp.fftshift(fp.fft2(im)))),
            np.log1p(np.abs(fp.fftshift(fp.fft2(im_blurred)))),
            titles=['Original spectrum', 'Blurred spectrum'])

# %% [markdown]
# ### 空間卷積與 FFT 卷積的時間比較
#
# 使用 Python 的 [`timeit`](https://docs.python.org/zh-tw/3/library/timeit.html) 量測執行時間。使用 3×3、標準差為 3 的核，每種方法量測 100 次。空間法明確指定 `method="direct"`，避免 `signal.convolve` 自動改選 FFT。先確認兩者輸出一致，再比較時間分布。小核的直接法未必比 FFT 慢；可在完成本例後自行改變核大小。

# %%
im = lucario
gauss_kernel = gaussian_psf((3, 3), std=3)
def direct_convolution():
    return signal.convolve(im, gauss_kernel, mode='same', method='direct')
def fft_convolution():
    return signal.fftconvolve(im, gauss_kernel, mode='same')
im_blurred1 = direct_convolution()
im_blurred2 = fft_convolution()
assert np.allclose(im_blurred1, im_blurred2, atol=1e-12)
times1 = timeit.repeat(direct_convolution, number=1, repeat=100)
times2 = timeit.repeat(fft_convolution, number=1, repeat=100)
show_images(im, im_blurred1, im_blurred2,
            titles=['Lucario', 'Direct convolution', 'FFT convolution'])

# %%
fig, ax = plt.subplots(figsize=(8, 5))
boxes = ax.boxplot([times1, times2], patch_artist=True,
                  tick_labels=['direct', 'FFT'])
for patch, color in zip(boxes['boxes'], ['#5aa9df', '#f2a65a']):
    patch.set_facecolor(color)
ax.set_ylabel('Time per call [s]')
ax.set_title('Same image, same normalized 3 x 3 kernel')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 頻率域濾波：高通與低通

# %% [markdown]
# ### 高通濾波器（High-Pass Filter, HPF）

# %% [markdown]
# 高通濾波器移除低頻，保留較快速的空間變化，因此常顯示邊緣與細節，也會保留同頻帶的雜訊。以下使用 Gengar，將中心的一塊方形低頻區設為零。

# %% [markdown]
# 計算 `fft2`，以 `fftshift` 找到中央低頻區，套上遮罩，再用 `ifftshift` 和 `ifft2` 回到影像。遮罩必須讓正負共軛頻率成對處理，才能得到實數輸出。

# %%
show_images(gengar, titles=['Gengar'])

# %%
gengar_fft = fp.fft2(gengar)
show_images(np.log1p(np.abs(fp.fftshift(gengar_fft))), titles=['Original spectrum'])

# %%
low_mask20 = square_lowpass(gengar.shape, 20)
high_mask20 = ~low_mask20
high_fft20 = gengar_fft * high_mask20
assert low_mask20.sum() == 41**2
show_images(fp.fftshift(high_mask20), np.log1p(np.abs(fp.fftshift(high_fft20))),
            titles=['High-pass mask: 41 x 41 center removed', 'Filtered spectrum'])

# %% [markdown]
# 高通結果含有正、負響應。使用以零為中心的對称色階，就能同時看見邊緣兩側；把負值截成零會改變結果。

# %%
high20 = real_ifft2(high_fft20)
limit = np.max(np.abs(high20))
fig, ax = plt.subplots(figsize=(5, 5))
ax.imshow(high20, cmap='RdBu_r', vmin=-limit, vmax=limit)
ax.set_title('Signed high-pass output: F=20')
ax.axis('off')
plt.show()

# %% [markdown]
# ### 從截止索引 1 看到 24
#
# 用同一張圖逐步加大被移除的中央方形。$F$ 是各軸保留／排除的最大頻率索引，方形邊長為 $2F+1$，不是一個徑向截止頻率。下方完整列出 24 個高通結果。

# %%
cutoffs = list(range(1, 25))
lowpass_results = [real_ifft2(gengar_fft * square_lowpass(gengar.shape, f)) for f in cutoffs]
highpass_results = [real_ifft2(gengar_fft * ~square_lowpass(gengar.shape, f)) for f in cutoffs]
for low, high in zip(lowpass_results, highpass_results):
    assert np.allclose(low + high, gengar, atol=1e-12)
montage(highpass_results, [f'High-pass F={f}' for f in cutoffs], signed=True)

# %% [markdown]
# ### 低通濾波器（Low-Pass Filter, LPF）

# %% [markdown]
# 低通保留中央的低頻區，移除外圍高頻。它能保留較大尺度的輪廓，但細節與銳利邊緣會一起減少。

# %% [markdown]
# 對同一組 FFT 係數套上互補遮罩：高通排除的中央方形，就是低通保留的區域。

# %% [markdown]
# 先看 $F=20$ 的結果。中央區塊包含每一軸索引 −20 到 20，因此是 **41×41**。

# %%
low_fft20 = gengar_fft * low_mask20
low20 = real_ifft2(low_fft20)
show_images(low20, titles=['Low-pass F=20'])

# %% [markdown]
# 查看套過低通遮罩的頻譜，確認中心保留、外圍變成零。

# %%
show_images(np.log1p(np.abs(fp.fftshift(low_fft20))), titles=['Low-pass spectrum: F=20'])

# %% [markdown]
# 用相同的 24 個截止索引計算低通結果；與前面的高通形成一組互補比較。

# %%
montage(lowpass_results, [f'Low-pass F={f}' for f in cutoffs])

# %% [markdown]
# $F$ 愈大，低通保留的頻帶愈寬，圖像細節逐漸增加。互動圖可以逐一切換截止索引，同時比較遮罩、低通和保留正負值的高通結果。

# %%
lab('filter_sweep', original=gengar, cutoffs=cutoffs,
    lowpass=lowpass_results, highpass=highpass_results,
    spectrum=np.log1p(np.abs(fp.fftshift(gengar_fft))))

# %% [markdown]
# ## 去卷積：已知模糊核時如何復原影像？

# %% [markdown]
# 先用標準差為 3 的 Gaussian 核模糊 Gengar，再使用同一個核的 Fourier transform 嘗試反演。這一次仍採同尺寸 FFT，模糊與復原使用一致的週期邊界。

# %%
im = gengar
gauss_kernel = gaussian_psf(im.shape, std=3)
freq_kernel = fp.fft2(fp.ifftshift(gauss_kernel))
im_blur = real_ifft2(fp.fft2(im) * freq_kernel)
show_images(im, gauss_kernel, im_blur,
            titles=['Gengar', 'Normalized Gaussian: std=3', 'Blurred image'])

# %% [markdown]
# ### 反濾波與差值圖
# 理想的無雜訊模型為 $G=HF$，在 $H\ne0$ 時可以除以 $H$。實際數值中，接近零的 $H$ 會放大捨入誤差。設定 $\epsilon=10^{-12}$，計算 $G/(H+\epsilon)$，並把重建與原圖的差值畫出來。

# %%
inverse_epsilon = 1e-12
inverse_transfer = 1 / (freq_kernel + inverse_epsilon)
im_restored = real_ifft2(fp.fft2(im_blur) * inverse_transfer)
inverse_difference = im_restored - im
print('inverse epsilon:', inverse_epsilon)
print('maximum absolute restoration error:', np.max(np.abs(inverse_difference)))
print('restoration MSE:', np.mean(inverse_difference**2))
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
for ax, image, title in zip(axes.flat,
        [im, im_blur, im_restored, inverse_difference],
        ['Gengar', 'Blurred', 'Approximate inverse', 'Restored minus original']):
    if title == 'Restored minus original':
        limit = max(1e-12, np.max(np.abs(image)))
        ax.imshow(image, cmap='RdBu_r', vmin=-limit, vmax=limit)
    else:
        ax.imshow(image)
    ax.set_title(title)
    ax.axis('off')
fig.tight_layout()
plt.show()

# %%
print('difference array shape:', inverse_difference.shape)
print(inverse_difference)

# %% [markdown]
# 差值可能來自浮點誤差，以及分母加入 $\epsilon$ 後的偏差。當某些頻率的核響應很弱，近似反演就不會完全還原原係數；不能從這張差值圖推論每次無雜訊去卷積都必然丟失資訊。

# %% [markdown]
# ```{dropdown} 為什麼含雜訊時除以 $H$ 特別危險？
# 若 $G=HF+N$，直接反演得到 $G/H=F+N/H$。當 $|H|$ 很小，雜訊會被大幅放大。加入 $\epsilon$ 後，訊號項也改成 $H/(H+\epsilon)F$，所以這是帶偏差的近似；它不是適用於所有模糊核的最佳正則化。一般反濾波器也不能直接等同於高通濾波器。
# ```

# %% [markdown]
# ### 帶雜訊的影像：`unsupervised_wiener`

# %% [markdown]
# 改用 2×2 box PSF，以線性卷積產生模糊圖，再加入標準差為模糊圖標準差一半的 Gaussian 雜訊。`unsupervised_wiener` 會從模型與資料估計相關的精度參數。固定亂數種子方便比較，並觀察它恢復的輪廓與殘留誤差。
#
# [`scikit-image restoration 文件`](https://scikit-image.org/docs/stable/api/skimage.restoration.html#skimage.restoration.unsupervised_wiener)說明此函式的模型、回傳值與參數。產生影像的 `same` 線性卷積與復原器的頻域模型在邊界處可能不同，判讀邊界時要把這個差別算進去。

# %%
im = gengar
n = 2
psf = np.ones((n, n)) / n**2
im1 = conv2(im, psf, 'same')
rng = np.random.default_rng(42)
im1 = im1 + 0.5 * im1.std() * rng.standard_normal(im1.shape)
im2, chains = ski.restoration.unsupervised_wiener(im1, psf, rng=42)
show_images(im, im1, im2,
            titles=['Gengar', 'Noisy blurred image', 'Self-tuned restoration'])

# %% [markdown]
# ## Moon landing：去除週期性干擾

# %% [markdown]
# 登月影像有規律條紋。先看原圖與頻譜，再比較「直接縮窄頻帶」及「針對強峰抑制」兩種做法。這張 `moonlanding.png` 本身已有干擾，直接分析它的頻譜。

# %%
image = ski.io.imread(image_path('moonlanding.png'))
assert image.ndim == 2
M, N = image.shape
show_images(image, titles=['Moon landing: original image'])
print(image.shape, image.dtype)

# %% [markdown]
# 二維傅立葉轉換（2-D FFT）等價於先對影像的每一列進行一維傅立葉轉換，再對每一欄進行一維傅立葉轉換（或反過來亦可）。

# %%
F = fp.fft2(image)
F_magnitude = fp.fftshift(np.abs(F))

# %% [markdown]
# 同樣地，在顯示之前，我們先對頻譜取對數，以壓縮數值範圍：

# %%
show_images(np.log1p(F_magnitude), titles=['Moon landing: log spectrum'], cmap='viridis')

# %%
keep_fraction = 0.1
ky_moon, kx_moon = frequency_indices(image.shape)
# Retain signed frequencies within 10% of each axis length; paired at both ends.
moon_low_mask = ((np.abs(ky_moon) < M * keep_fraction)
                 & (np.abs(kx_moon) < N * keep_fraction))
im_fft2 = F * moon_low_mask
show_images(fp.fftshift(moon_low_mask),
            np.log1p(np.abs(fp.fftshift(im_fft2))),
            titles=['Symmetric low-pass mask', 'Filtered spectrum'])

# %%
im_new = real_ifft2(im_fft2)
show_images(im_new, titles=['Moon landing: low-pass reconstruction'])

# %% [markdown]
# 低通會一起移除干擾與高頻細節。接著改用頻譜峰值找候選干擾：先保護中心低頻區，再將其餘幅度超過第 98 百分位的係數成對壓掉。強峰也可能來自真實紋理，因此要比較輸出影像，不能把所有高頻峰都當成雜訊。

# %%
F = fp.fft2(image)
K = 40
# A symmetric 81 x 81 protected center replaces the asymmetric [-40, 39] slice.
protected_center = square_lowpass(image.shape, K)
peak_scores = np.abs(F).copy()
peak_scores[protected_center] = 0
threshold98 = np.percentile(peak_scores, 98)
remove_peaks = (peak_scores >= threshold98) & ~protected_center
neg_rows = (-np.arange(M)) % M
neg_cols = (-np.arange(N)) % N
remove_peaks |= remove_peaks[np.ix_(neg_rows, neg_cols)]
F_dim = F * ~remove_peaks
image_filtered = real_ifft2(F_dim)
print('98th-percentile amplitude:', threshold98)
print('removed Fourier coefficients:', np.count_nonzero(remove_peaks))
show_images(np.log1p(np.abs(fp.fftshift(F_dim))),
            titles=['Spectrum after paired peak suppression'], cmap='viridis')
show_images(image, im_new, image_filtered,
            titles=['Original moon landing', 'Low-pass', 'Peak suppression'])

# %% [markdown]
# ## 影像金字塔：在多個尺度看同一張圖

# %% [markdown]
# 影像中的物件可能有不同大小。尋找人臉時，可以先建立一組由大到小的影像，在各尺度中搜尋；編輯影像時，也可以分開處理大尺度明暗與局部細節。這就是影像金字塔的用途。
#
# Gaussian pyramid 從原圖開始，每次先平滑，再取較少的像素，逐層得到較小的影像。下方兩張原始示意圖展示金字塔的層級與尺度關係。
#
# 可重建的 residual Laplacian pyramid 儲存「某一層 Gaussian image 減去下一個粗層放大後的影像」，再保留最粗層；重建時按同一個放大操作逐層加回差值。稍後呼叫的 scikit-image `pyramid_laplacian` 採用另一種每層高通定義，兩者的輸出不能直接混用。

# %%
# 圖片來源： https://drive.google.com/uc?id=1BLZ-Chch9k3r9EHNsnkbi75maq_a7n-y
display(Image(filename=str(image_path('ch03-cell101-1.png'))))

# %%
# 圖片來源： https://drive.google.com/uc?id=1BOq0GkhwEvPCRG6YGpZCClzvcVryP1RQ
display(Image(filename=str(image_path('ch03-cell101-2.png'))))

# %% [markdown]
# ### Pikachu 的 Gaussian pyramid
#
# `pyramid_gaussian(image, downscale=2, channel_axis=-1)` 每次先平滑，再縮小。保留各層尺寸與完整 montage，觀察細節如何隨尺度消失。

# %%
pikachu_rgb = read_rgb('pikachu.png')
image = pikachu_rgb
pyramid = tuple(ski.transform.pyramid_gaussian(image, downscale=2, channel_axis=-1))
show_images(*pyramid, titles=[f'{p.shape[0]} x {p.shape[1]}' for p in pyramid],
            figsize=(20, 4))
nrows, ncols = image.shape[:2]
right_height = sum(p.shape[0] for p in pyramid[1:])
right_width = max((p.shape[1] for p in pyramid[1:]), default=0)
composite_gaussian = np.zeros((max(nrows, right_height), ncols + right_width, 3))
composite_gaussian[:nrows, :ncols] = pyramid[0]
i_row = 0
for p in pyramid[1:]:
    nr, nc = p.shape[:2]
    composite_gaussian[i_row:i_row + nr, ncols:ncols + nc] = p
    i_row += nr
show_images(composite_gaussian, titles=['Pikachu: complete Gaussian pyramid'], figsize=(9, 8))

# %% [markdown]
# ### `pyramid_laplacian` 的各尺度細節
#
# 這個 API 的每層是「當層影像減去其平滑版本」，接著從平滑、縮小後的影像繼續計算。因此它會凸顯各尺度的細節，但不是上一段提到的「相鄰 Gaussian 層－expand 粗層」residual 定義。本例顯示各層及 montage，不宣稱用這些回傳值即可按另一套公式精確重建。

# %%
laplacian_pyramid = tuple(ski.transform.pyramid_laplacian(
    pikachu_rgb, downscale=2, channel_axis=-1))
laplacian_gray = [ski.color.rgb2gray(p) for p in laplacian_pyramid]
montage(laplacian_gray, [f'{p.shape[0]} x {p.shape[1]}' for p in laplacian_gray],
        signed=True, figsize=(12, 8))
nrows, ncols = laplacian_gray[0].shape
right_height = sum(p.shape[0] for p in laplacian_gray[1:])
right_width = max((p.shape[1] for p in laplacian_gray[1:]), default=0)
composite_laplacian = np.zeros((max(nrows, right_height), ncols + right_width))
composite_laplacian[:nrows, :ncols] = laplacian_gray[0]
i_row = 0
for p in laplacian_gray[1:]:
    nr, nc = p.shape
    composite_laplacian[i_row:i_row + nr, ncols:ncols + nc] = p
    i_row += nr
limit = max(1e-12, np.max(np.abs(composite_laplacian)))
fig, ax = plt.subplots(figsize=(9, 8))
ax.imshow(composite_laplacian, cmap='RdBu_r', vmin=-limit, vmax=limit)
ax.set_title('Pikachu: complete pyramid_laplacian output')
ax.axis('off')
plt.show()

# %% [markdown]
# Laplacian 層含正、負值，以零為中心顯示才能看清楚局部亮暗差異；Gaussian 層則保留平滑後的影像。兩者都保留空間位置，只是呈現的尺度與內容不同。

# %% [markdown]
# ### 金字塔也能用來融合影像
#
# 若直接把兩張圖各切一半接起來，接縫通常很明顯。多尺度融合讓低頻的過渡較寬，高頻的細節較局部；可沿著下列範例查看如何建立兩張圖與遮罩的金字塔。
#
# - [Image blending using Laplacian pyramids](https://becominghuman.ai/image-blending-using-laplacian-pyramids-2f8e9982077f)：查看不同層的融合與合成順序。
# - [OpenCV Image Pyramids](https://docs.opencv.org/4.x/dc/dff/tutorial_py_pyramids.html)：由 Gaussian pyramid 讀到 Laplacian pyramid 的融合範例。

# %% [markdown]
# ## 特徵偵測：角點與斑點

# %% [markdown]
# 前面用頻率與尺度描述整張影像；特徵偵測則選出值得注意的位置，例如角點或特定大小的斑點。這些位置可以供後續比對、辨識或量測使用。

# %% [markdown]
# ### Harris 角點：棋盤與足球

# %% [markdown]
# 邊緣主要沿一個方向有強烈變化；在角點附近，視窗往不同方向移動都會遇到明顯的強度變化。Harris response 就用這個差別找候選角點。設定 `k=0.001`，將大於最大response 1%的區域塗紅。這是response閾值圖，尚未使用局部極大值抑制把每一區縮成一個點。
#
# [CS131範例notebook](https://github.com/mikucy/CS131/blob/master/hw3_release/hw3.ipynb)提供角點與特徵比對的延伸練習。

# %%
image = read_rgb('chess_football.png')
image_gray = ski.color.rgb2gray(image)
harris_response = ski.feature.corner_harris(image_gray, k=0.001)
harris_regions = harris_response > 0.01 * harris_response.max()
harris_overlay = image.copy()
harris_overlay[harris_regions] = [1.0, 0.0, 0.0]
show_images(image, harris_response, harris_overlay,
            titles=['Chess and football', 'Harris response', 'Response above 1% of maximum'])

# %% [markdown]
# ### Hubble 影像：LoG、DoG 與 DoH 斑點偵測

# %% [markdown]
# 高斯模糊會將影像「平滑化」。在對比度變化非常小的區域，即使一張影像的模糊程度比另一張強，這些區域看起來仍然相似。當某區域變化很小時，兩張影像相減的結果接近零（黑色）。而在高對比度的區域（例如邊緣或斑點），模糊強度的影響較大。

# %%
# 圖片來源： https://drive.google.com/uc?id=1C-bXNQH9ad1HNGMHvqh4XUJjHAoZ8FcZ
display(Image(filename=str(image_path('ch03-cell114-1.png'))))

# %%
# 圖片來源： https://drive.google.com/uc?id=1BPel1EjTj2gqHeHNlT2uz80CaXhd2vlr
display(Image(filename=str(image_path('ch03-cell114-2.png'))))

# %% [markdown]
# https://medium.com/@vad710/cv-for-busy-devs-improving-features-df20c3aa5887

# %%
image = ski.data.hubble_deep_field()[0:500, 0:500]
image_gray = ski.color.rgb2gray(image)

# Note there is a parameter sigma_ratio
blobs_log = ski.feature.blob_log(image_gray, max_sigma=30, num_sigma=10, threshold=.1)
# Compute radii in the 3rd column.
blobs_log[:, 2] = blobs_log[:, 2] * math.sqrt(2)

blobs_dog = ski.feature.blob_dog(image_gray, max_sigma=30, threshold=.1)
blobs_dog[:, 2] = blobs_dog[:, 2] * math.sqrt(2)

blobs_doh = ski.feature.blob_doh(image_gray, max_sigma=30, threshold=.01)

blobs_list = [blobs_log, blobs_dog, blobs_doh]
colors = ['yellow', 'lime', 'red']
titles = ['Laplacian of Gaussian', 'Difference of Gaussian',
          'Determinant of Hessian']
sequence = zip(blobs_list, colors, titles)

fig, axes = plt.subplots(1, 3, figsize=(9, 3), sharex=True, sharey=True)
ax = axes.ravel()

for idx, (blobs, color, title) in enumerate(sequence):
    ax[idx].set_title(title)
    ax[idx].imshow(image)
    for blob in blobs:
        y, x, r = blob
        c = plt.Circle((x, y), r, color=color, linewidth=2, fill=False)
        ax[idx].add_patch(c)
    ax[idx].set_axis_off()

plt.tight_layout()
plt.show();

# %% [markdown]
# 比較三種方法畫出的圓：中心與尺度是否相近？某些模糊小點是否只被一種方法找出？LoG／DoG 的二維半徑以 $\sqrt2\sigma$ 近似，DoH 的輸出尺度則直接當作半徑。三者閾值的意義不同，要將中心位置、尺度與實際物件一起比較，才能判斷結果。
#
# - [scikit-image feature API](https://scikit-image.org/docs/stable/api/skimage.feature.html)：查看三種blob函式的尺度與閾值定義。
# - [CS131特徵實作](https://github.com/mikucy/CS131/blob/master/hw3_release/hw3.ipynb)：延伸到特徵描述與比對。

# %% [markdown]
# ## 延伸閱讀與參考資料

# %% [markdown]
# - [OpenCV Image Processing](https://docs.opencv.org/4.x/d2/d96/tutorial_py_table_of_contents_imgproc.html)：從濾波、頻域操作到影像金字塔，依本章單元對照實作。
# - [OpenCV Feature Detection and Description](https://docs.opencv.org/4.x/db/d27/tutorial_py_table_of_contents_feature2d.html)：從角點偵測接到描述子與影像匹配。
# - 課程指定教科書第 3.4、3.5 與 7.1 節：依課堂使用的版本查閱，對照濾波、金字塔與特徵偵測。
# - [Hands-On Image Processing with Python 原始碼](https://github.com/PacktPublishing/Hands-On-Image-Processing-with-Python)：本章moonlanding、chess_football影像與頻域實驗的來源之一。
# - [SciPy FFT文件](https://docs.scipy.org/doc/scipy/reference/fft.html)：核對頻率排序、正規化與多維轉換。
# - [scikit-image金字塔API](https://scikit-image.org/docs/stable/api/skimage.transform.html#skimage.transform.pyramid_laplacian)：特別比較本章兩種Laplacian定義。
