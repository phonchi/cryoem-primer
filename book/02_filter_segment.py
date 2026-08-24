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
# # 濾波、匹配與分割
#
# 本章把三類常被混在一起的工作拆開：濾波改變像素，匹配產生候選位置，分割則替像素或物件編號。它們都可能讓圖「比較好看」，但只有明確定義輸入、邊界與評估方式，結果才可重現。
#
# ```{admonition} 學習目標
# :class: important
#
# - 比較 convolution 與 correlation，並說明 boundary mode 的影響。
# - 驗證 2D Gaussian filter 可以拆成兩次 1D filtering。
# - 使用 local normalization、normalized cross-correlation（NCC）與 Difference of Gaussians（DoG）找候選物件。
# - 用 distance transform、local maxima 與 masked watershed 分開相鄰物件。
# - 說明 template bias 與非線性前處理對 cryo-EM 分析的限制。
# ```

# %%
import matplotlib.pyplot as plt
import numpy as np
import skimage as ski
from scipy import ndimage as ndi

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
# (filtering)=
# ## Convolution、correlation 與邊界
#
# 對離散訊號 $x$ 與 kernel $h$，convolution 會翻轉 kernel；correlation 不翻轉。kernel 若對稱，兩者相同。影像是有限陣列，kernel 到邊緣時還需要指定陣列外的值：補零、重複邊緣、鏡射或週期延拓會得到不同答案。

# %%
step = np.zeros(20)
step[7:15] = 1
asymmetric_kernel = np.array([1.0, 0.0, -1.0])

conv = ndi.convolve1d(step, asymmetric_kernel, mode="reflect")
corr = ndi.correlate1d(step, asymmetric_kernel, mode="reflect")
assert np.allclose(conv, -corr)

edge_signal = np.array([1.0, 0.0, 0.0, 0.0])
mean_kernel = np.full(3, 1 / 3)
boundary_results = {
    mode: ndi.convolve1d(edge_signal, mean_kernel, mode=mode)
    for mode in ("constant", "nearest", "reflect", "wrap")
}
for mode, values in boundary_results.items():
    print(f"{mode:>8}:", np.round(values, 3))

# %% [markdown]
# ```{admonition} 邊界不是小細節
# :class: caution
#
# FFT multiplication naturally corresponds to circular convolution；一般空間域的 `same` convolution 則常假設陣列外為零。若沒有先 padding 與 crop，兩者不能直接比較。粒子盒邊緣的訊號尤其容易因 boundary mode 產生假影。
# ```


# %% [markdown]
# ## Gaussian filter 為什麼可以拆成兩次 1D filtering？
#
# Isotropic 2D Gaussian kernel 是兩個 1D Gaussian 的 outer product，因此是 separable filter。先沿列方向，再沿欄方向濾波，和直接套 2D kernel 應在浮點誤差內一致；計算量則從每像素約 $K^2$ 次乘加降到 $2K$。

# %%
radius = 4
sigma = 1.5
axis = np.arange(-radius, radius + 1)
gaussian_1d = np.exp(-(axis**2) / (2 * sigma**2))
gaussian_1d /= gaussian_1d.sum()
gaussian_2d = np.outer(gaussian_1d, gaussian_1d)

rng = np.random.default_rng(42)
image = rng.normal(size=(48, 64))
direct_2d = ndi.convolve(image, gaussian_2d, mode="reflect")
separable = ndi.convolve1d(image, gaussian_1d, axis=0, mode="reflect")
separable = ndi.convolve1d(separable, gaussian_1d, axis=1, mode="reflect")
assert np.allclose(direct_2d, separable, atol=1e-12)
print("max |2D - separable| =", np.max(np.abs(direct_2d - separable)))


# %% [markdown]
# ## 線性與非線性平滑
#
# Gaussian filter 是線性低通；median filter 是非線性排序操作，對 salt-and-pepper outliers 很有效，但會改變細小結構與像素統計。它適合當作教學比較或處理已確認的脈衝型污染，不是所有 cryo-EM micrograph 的預設前處理。

# %%
coins = ski.util.img_as_float(ski.data.coins())
coins_impulse = ski.util.random_noise(coins, mode="s&p", amount=0.03, rng=0)
gaussian = ski.filters.gaussian(coins_impulse, sigma=1.2)
median = ndi.median_filter(coins_impulse, size=3, mode="reflect")
imshow_all(coins_impulse, gaussian, median,
           titles=["impulse noise", "Gaussian", "median"])


# %% [markdown]
# ## Local normalization 與 whitening 解決的問題不同
#
# Local normalization 以局部平均與標準差校正緩慢變化的背景：
#
# $$z(\mathbf{x})=\frac{I(\mathbf{x})-\mu_{\mathrm{local}}(\mathbf{x})}
# {\sqrt{\sigma^2_{\mathrm{local}}(\mathbf{x})+\epsilon}}.$$
#
# Whitening 則在頻域依背景 power spectrum 重新加權頻率。兩者都常出現在 particle picking 的前處理，但不能互相冒充；以下只示範 local normalization。

# %%
def local_normalize(image, sigma=12.0, eps=1e-6):
    """以 Gaussian-weighted local moments 做正規化。"""
    image = np.asarray(image, dtype=float)
    mean = ndi.gaussian_filter(image, sigma=sigma, mode="reflect")
    mean_square = ndi.gaussian_filter(image**2, sigma=sigma, mode="reflect")
    variance = np.maximum(mean_square - mean**2, 0.0)
    return (image - mean) / np.sqrt(variance + eps)


rows, cols = coins.shape
background_ramp = np.linspace(0, 0.7, cols)[None, :]
uneven = coins + background_ramp
normalized = local_normalize(uneven)
assert np.isfinite(normalized).all()
imshow_all(uneven, normalized, titles=["uneven background", "local normalization"])


# %% [markdown]
# ## Normalized cross-correlation（NCC）
#
# Raw correlation 會偏好整體較亮或能量較大的區域。NCC 先扣除局部與模板平均，再除以兩者能量；理想情況下分數落在 $[-1,1]$，較能比較背景亮度不同的位置。`match_template()` 實作的就是這類正規化匹配。

# %%
search_image = ski.util.img_as_float(ski.data.coins())
template = search_image[170:215, 75:120]
ncc = ski.feature.match_template(search_image, template, pad_input=True)
peak_row, peak_col = np.unravel_index(np.argmax(ncc), ncc.shape)
assert ncc.max() <= 1 + 1e-12

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].imshow(search_image)
axes[0].plot(peak_col, peak_row, "r+")
axes[0].set_title("best NCC location")
axes[1].imshow(ncc, vmin=-1, vmax=1)
axes[1].set_title("NCC score")
for ax in axes:
    ax.axis("off")
fig.tight_layout()

# %% [markdown]
# ```{admonition} Template bias
# :class: caution
#
# NCC 只會找「像模板」的區域。模板來源、低通截止、取向覆蓋與 score threshold 都會改變 picks，還可能偏向預先期待的結構。實務上要用獨立資料、negative controls、2D classification 與多樣性檢查評估，而不能把高 NCC 分數直接當成粒子真值。
# ```


# %% [markdown]
# (particle-picking)=
# ## Difference of Gaussians（DoG）blob detection
#
# 兩個尺度的 Gaussian-smoothed images 相減，形成 band-pass-like response。`blob_dog()` 在位置與尺度上找 extrema，適合產生近圓形物件的候選點；它不理解分子種類，也不會自動排除冰污染。

# %%
blobs = ski.feature.blob_dog(search_image, min_sigma=5, max_sigma=18,
                             sigma_ratio=1.4, threshold=0.08, overlap=0.5)
assert blobs.ndim == 2 and blobs.shape[1] == 3

fig, ax = plt.subplots(figsize=(7, 5))
ax.imshow(search_image)
for row, col, sigma_found in blobs[:80]:
    circle = plt.Circle((col, row), np.sqrt(2) * sigma_found,
                        fill=False, color="tab:red", linewidth=0.8)
    ax.add_patch(circle)
ax.set_title(f"DoG candidates (showing {min(len(blobs), 80)} of {len(blobs)})")
ax.axis("off")
fig.tight_layout()


# %% [markdown]
# ## 用 watershed 分開相鄰物件
#
# 二值 mask 只回答前景／背景。要讓每個相鄰物件得到獨立 label，可以對前景做 Euclidean distance transform，以局部 maxima 當作物件 seeds，再在 `-distance` 上執行 masked watershed。背景維持 label 0，`regionprops()` 就不會把背景當成一個物件。

# %%
def split_touching_objects(mask, min_distance=12):
    """回傳 distance map、seed markers 與 instance labels。"""
    mask = np.asarray(mask, dtype=bool)
    distance = ndi.distance_transform_edt(mask)
    coordinates = ski.feature.peak_local_max(
        distance, labels=mask, min_distance=min_distance, exclude_border=False
    )
    seeds = np.zeros(mask.shape, dtype=bool)
    seeds[tuple(coordinates.T)] = True
    markers, marker_count = ndi.label(seeds)
    labels = ski.segmentation.watershed(-distance, markers, mask=mask)
    return distance, markers, labels, marker_count


yy, xx = np.mgrid[:160, :220]
touching_mask = (
    (yy - 80) ** 2 + (xx - 75) ** 2 < 42**2
) | (
    (yy - 80) ** 2 + (xx - 137) ** 2 < 42**2
)
distance, markers, instance_labels, marker_count = split_touching_objects(
    touching_mask, min_distance=35
)
props = ski.measure.regionprops(instance_labels)

assert marker_count == 2
assert len(props) == 2
assert instance_labels[~touching_mask].max(initial=0) == 0
print("instance count:", len(props), "areas:", [p.area for p in props])

imshow_all(touching_mask, distance, markers, instance_labels,
           titles=["binary mask", "distance", "seeds", "watershed labels"])


# %% [markdown]
# ## HOG：方向分布的描述子
#
# Histogram of Oriented Gradients（HOG）把每個 cell 的梯度方向累積成直方圖，再用相鄰 cells 組成 block 做正規化。這裡使用 `(2,2)` cells per block，確實包含跨 cell 的 block normalization；HOG 是通用電腦視覺示範，不是現代 cryo-EM picking 的完整模型。

# %%
features, hog_image = ski.feature.hog(
    search_image,
    orientations=9,
    pixels_per_cell=(16, 16),
    cells_per_block=(2, 2),
    visualize=True,
)
assert features.ndim == 1 and features.size > 0
imshow_all(search_image, hog_image, titles=["input", "HOG visualization"])


# %% [markdown]
# ## 與 cryo-EM 的連結
#
# Particle picking 常會使用 band-pass、whitening 或 local normalization 改善候選偵測的條件；Gaussian 與 median filter 在本章主要是教學比較，median filter 只適合有明確脈衝型污染的特定情境。任何非線性濾波都可能改變微弱高解析訊號，不應在未驗證下套到 reconstruction-bound particles。
#
# 模板匹配與 blob detection 可以產生候選座標，但候選不等於真實粒子。模板偏差、preferred orientation、污染與 ice thickness 都會影響被挑出的分布。形態學與 watershed 適合清理污染遮罩或示範 instance segmentation；它們不是把重疊投影「物理解混」的工具。
#
# ## 理解檢查
#
# 1. 為什麼同一個 mean kernel 在 `constant` 與 `reflect` 邊界會得到不同結果？
# 2. Local normalization 與 whitening 分別校正哪一種不均勻？
# 3. NCC 分數高，為什麼仍不能斷言該位置是真實粒子？
# 4. Watershed 若只有「背景 marker 1」與「前景 marker 2」，為什麼不能分開兩個相鄰物件？
