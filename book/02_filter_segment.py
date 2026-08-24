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
# # 濾波、邊緣、分割與匹配
#
# 本章從局部濾波開始，依序處理邊緣、閾值化、形態學、物件分割與特徵匹配。這幾類操作常被串在同一條流程裡，任務卻不相同：濾波改變像素，邊緣偵測估計強度變化，分割替像素或物件編號，匹配則產生候選位置。讀完後，應能根據雜訊種類與分析目標選擇方法，而不只是挑一張看起來最平滑的結果。基本觀念可參考 {cite}`szeliski2022,forsyth2012`，實作細節則以 [SciPy `ndimage`](https://docs.scipy.org/doc/scipy/reference/ndimage.html) 與 [scikit-image API](https://scikit-image.org/docs/stable/api/api.html) 為準。
#
# ```{admonition} 學習目標
# :class: important
#
# - 比較 convolution 與 correlation，並說明 boundary mode 的影響。
# - 依雜訊種類比較 box、Gaussian、median 與 bilateral filtering。
# - 從 $G_x$、$G_y$ 建立梯度大小與方向，解釋 Sobel、Gaussian derivative、Laplacian、LoG、DoG、unsharp masking 與 Canny。
# - 比較固定、Otsu 與局部閾值，並用 IoU 量化分割結果。
# - 使用侵蝕、膨脹、開運算、閉運算、connected components 與 `regionprops()` 整理二值遮罩。
# - 使用 local normalization、normalized cross-correlation（NCC）與 Difference of Gaussians（DoG）找候選物件。
# - 用 distance transform、local maxima 與 masked watershed 分開相鄰物件。
# - 推導 HOG 特徵向量的維度，並說明 template bias 與非線性前處理對 cryo-EM 分析的限制。
# ```

# %%
import matplotlib.pyplot as plt
import numpy as np
import skimage as ski
from scipy import ndimage as ndi

plt.rcParams["image.cmap"] = "gray"
plt.rcParams["figure.figsize"] = (6, 5)


def imshow_all(*images, titles=None, size=4, **kwargs):
    """並排顯示影像，供教學比較使用。"""
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
# ## 常用的平滑濾波器
#
# Box filter 對視窗內像素給相同權重；Gaussian filter 讓鄰近中心的像素占較大權重。兩者都是線性低通，會一併削弱雜訊與細節。Median filter 以鄰域中位數取代中心值，是非線性排序操作；bilateral filter 則同時依空間距離與像素差異加權，較不容易跨越強邊緣平均。這些差異來自運算定義，不能簡化成某一種濾波器「最好」。相關推導見 {cite}`szeliski2022,forsyth2012`；函式參數可查 [scikit-image restoration 文件](https://scikit-image.org/docs/stable/api/skimage.restoration.html)。

# %%
coins = ski.util.img_as_float(ski.data.coins())
coins_impulse = ski.util.random_noise(coins, mode="s&p", amount=0.03, rng=0)
box_filtered = ndi.uniform_filter(coins_impulse, size=5, mode="reflect")
gaussian_filtered = ndi.gaussian_filter(coins_impulse, sigma=1.2, mode="reflect")
median_filtered = ndi.median_filter(coins_impulse, size=3, mode="reflect")
bilateral_filtered = ski.restoration.denoise_bilateral(
    coins_impulse, sigma_color=0.08, sigma_spatial=3, channel_axis=None
)
imshow_all(
    coins_impulse,
    box_filtered,
    gaussian_filtered,
    median_filtered,
    bilateral_filtered,
    titles=["impulse noise", "box", "Gaussian", "median", "bilateral"],
    size=3,
)


# %% [markdown]
# ### 雜訊模型會改變比較結果
#
# Additive Gaussian noise 會輕微擾動多數像素；salt-and-pepper noise 則把少數像素推到動態範圍兩端。下面使用同一張乾淨影像與固定亂數種子，再以 mean squared error（MSE）比較結果。這個例子裡，Gaussian smoothing 對 Gaussian noise 的 MSE 較低，而 median filtering 對脈衝雜訊較有利。結論只適用於這組雜訊強度與參數；改變細節尺度、kernel 大小或評估指標，排序也可能改變。

# %%
clean = ski.util.img_as_float(ski.data.camera())[80:336, 80:336]
rng = np.random.default_rng(2026)
gaussian_noise = np.clip(clean + rng.normal(0, 0.08, clean.shape), 0, 1)
impulse_noise = ski.util.random_noise(
    clean, mode="s&p", amount=0.04, rng=2026
)


def compare_gaussian_median(clean_image, noisy_image):
    """回傳 Gaussian 與 median filtering 的結果和 MSE。"""
    smooth_gaussian = ndi.gaussian_filter(noisy_image, sigma=1, mode="reflect")
    smooth_median = ndi.median_filter(noisy_image, size=3, mode="reflect")
    errors = {
        "noisy": ski.metrics.mean_squared_error(clean_image, noisy_image),
        "Gaussian": ski.metrics.mean_squared_error(clean_image, smooth_gaussian),
        "median": ski.metrics.mean_squared_error(clean_image, smooth_median),
    }
    return smooth_gaussian, smooth_median, errors


gaussian_smooth, gaussian_median, gaussian_errors = compare_gaussian_median(
    clean, gaussian_noise
)
impulse_smooth, impulse_median, impulse_errors = compare_gaussian_median(
    clean, impulse_noise
)

assert gaussian_errors["Gaussian"] < gaussian_errors["median"]
assert impulse_errors["median"] < impulse_errors["Gaussian"]
print("Gaussian noise MSE:", {k: round(v, 5) for k, v in gaussian_errors.items()})
print("Impulse noise MSE:", {k: round(v, 5) for k, v in impulse_errors.items()})

fig, axes = plt.subplots(2, 4, figsize=(13, 7))
rows_to_show = [
    (gaussian_noise, gaussian_smooth, gaussian_median, "Gaussian noise"),
    (impulse_noise, impulse_smooth, impulse_median, "impulse noise"),
]
for row_axes, (noisy, smooth_g, smooth_m, noise_name) in zip(axes, rows_to_show):
    for ax, shown, title in zip(
        row_axes,
        (clean, noisy, smooth_g, smooth_m),
        ("clean", noise_name, "Gaussian", "median"),
    ):
        ax.imshow(shown)
        ax.set_title(title)
        ax.axis("off")
fig.tight_layout()

# %% [markdown]
# ```{admonition} 與 cryo-EM 的連結
# :class: note
#
# 真實 micrograph 的背景通常不是獨立同分布的白色 Gaussian noise，還會受偵測器、冰層、beam-induced motion 與 CTF 影響。上面的比較用來理解濾波器，不足以替 cryo-EM 前處理選定唯一方法。Median filter 只有在已確認脈衝型污染時才有明確對應；任何非線性濾波都可能改變微弱的高解析訊號。
# ```


# %% [markdown]
# ## 從一階導數找邊緣
#
# 對影像 $I(x,y)$，水平與垂直偏導數可寫成 $G_x=\partial I/\partial x$、$G_y=\partial I/\partial y$。梯度大小與方向分別為
#
# $$
# \lVert\nabla I\rVert=\sqrt{G_x^2+G_y^2},\qquad
# \theta=\operatorname{atan2}(G_y,G_x).
# $$
#
# Sobel kernel 同時做局部差分與少量平滑。另一種做法是以 Gaussian derivative 直接估計平滑後影像的導數；`sigma` 越大，對雜訊越不敏感，但定位也會變得較粗。邊緣偵測與尺度選擇的關係可參考 {cite}`szeliski2022,forsyth2012`。

# %%
edge_input = clean
gx = ndi.sobel(edge_input, axis=1, mode="reflect") / 8
gy = ndi.sobel(edge_input, axis=0, mode="reflect") / 8
gradient_magnitude = np.hypot(gx, gy)
gradient_orientation = np.arctan2(gy, gx)

gaussian_gx = ndi.gaussian_filter(
    edge_input, sigma=1.2, order=(0, 1), mode="reflect"
)
gaussian_gy = ndi.gaussian_filter(
    edge_input, sigma=1.2, order=(1, 0), mode="reflect"
)
gaussian_gradient = np.hypot(gaussian_gx, gaussian_gy)

assert gradient_magnitude.shape == edge_input.shape
assert np.isfinite(gradient_orientation).all()
assert np.all(gradient_magnitude + 1e-15 >= np.abs(gx))

fig, axes = plt.subplots(2, 3, figsize=(11, 7))
edge_panels = [
    (edge_input, "input", {}),
    (gx, "$G_x$ (Sobel)", {"cmap": "coolwarm"}),
    (gy, "$G_y$ (Sobel)", {"cmap": "coolwarm"}),
    (gradient_magnitude, "gradient magnitude", {}),
    (gradient_orientation, "orientation", {"cmap": "twilight", "vmin": -np.pi, "vmax": np.pi}),
    (gaussian_gradient, "Gaussian derivative", {}),
]
for ax, (shown, title, options) in zip(axes.flat, edge_panels):
    ax.imshow(shown, **options)
    ax.set_title(title)
    ax.axis("off")
fig.tight_layout()


# %% [markdown]
# ## 二階導數、DoG 與銳化
#
# Laplacian $\nabla^2 I$ 是二階導數，會在快速變化處產生正負響應，也會強烈放大高頻雜訊。Laplacian of Gaussian（LoG）先用 Gaussian 控制尺度，再計算 Laplacian。Difference of Gaussians（DoG）以兩個相近尺度的 Gaussian 結果相減，可近似尺度正規化後的 LoG；它既可用來找 blob，也可看成帶通響應。Unsharp masking 則把原圖與低通影像之差加回原圖：
#
# $$I_{\mathrm{sharp}}=I+\alpha\left(I-G_\sigma*I\right).$$
#
# 銳化增加局部對比，也會放大雜訊與 ringing，不能視為新增解析資訊。

# %%
smoothed_for_second_order = ndi.gaussian_filter(edge_input, sigma=1.2)
laplacian_response = ndi.laplace(smoothed_for_second_order, mode="reflect")
log_response = ndi.gaussian_laplace(edge_input, sigma=1.2, mode="reflect")
dog_response = (
    ndi.gaussian_filter(edge_input, sigma=1.0, mode="reflect")
    - ndi.gaussian_filter(edge_input, sigma=1.6, mode="reflect")
)
unsharp = ski.filters.unsharp_mask(edge_input, radius=1.5, amount=1.2)

assert np.isfinite(log_response).all()
assert abs(dog_response.mean()) < 1e-3

imshow_all(
    laplacian_response,
    log_response,
    dog_response,
    unsharp,
    titles=["Laplacian after smoothing", "LoG", "DoG", "unsharp masking"],
    size=3.2,
)


# %% [markdown]
# ## Canny：從梯度到單像素邊緣
#
# Canny detector 串起四個步驟：Gaussian smoothing、梯度估計、non-maximum suppression，以及由高低兩個閾值控制的 hysteresis。高閾值先保留可靠邊緣；低閾值只接受與可靠邊緣相連的弱響應。這比直接對梯度大小切一個閾值多了邊緣細化與連通性判斷。完整參數定義見 [scikit-image Canny 文件](https://scikit-image.org/docs/stable/api/skimage.feature.html#skimage.feature.canny)。

# %%
canny_without_smoothing = ski.feature.canny(
    gaussian_noise, sigma=0, low_threshold=0.05, high_threshold=0.15
)
canny_smoothed = ski.feature.canny(
    gaussian_noise, sigma=2, low_threshold=0.05, high_threshold=0.15
)
assert canny_smoothed.sum() < canny_without_smoothing.sum()
imshow_all(
    gaussian_noise,
    canny_without_smoothing,
    canny_smoothed,
    titles=["noisy input", "Canny, sigma=0", "Canny, sigma=2"],
)


# %% [markdown]
# ## Local normalization 與 whitening 解決的問題不同
#
# Local normalization 以局部平均與標準差校正緩慢變化的背景：
#
# $$z(\mathbf{x})=\frac{I(\mathbf{x})-\mu_{\mathrm{local}}(\mathbf{x})}
# {\sqrt{\sigma^2_{\mathrm{local}}(\mathbf{x})+\epsilon}}.$$
#
# Whitening 則在頻域依背景 power spectrum 重新加權頻率。兩者都可能出現在粒子挑選（particle picking）的前處理，但不能互相代替；以下只示範 local normalization。

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
# NCC 只會找「像模板」的區域。模板來源、低通截止、取向覆蓋與分數閾值都會改變候選座標，還可能偏向預先期待的結構。實務上要用獨立資料、negative controls、2D classification 與多樣性檢查評估，不能把高 NCC 分數直接當成粒子真值。
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
# ## 固定、Otsu 與局部閾值
#
# 閾值化只根據強度把像素分成前景與背景。固定閾值容易解釋，卻依賴影像尺度；Otsu 方法從全域直方圖選出能最大化類別間變異的閾值；局部閾值則讓閾值隨位置改變，較能處理緩慢變化的背景。它們都不會自動理解「物件」是什麼。方法背景見 {cite}`szeliski2022,forsyth2012,howse2020`，API 行為見 [scikit-image `threshold_otsu` 文件](https://scikit-image.org/docs/stable/api/skimage.filters.html#skimage.filters.threshold_otsu)。
#
# 已知真值時，可用 intersection over union（IoU）量化二值分割：
#
# $$\operatorname{IoU}(A,B)=\frac{|A\cap B|}{|A\cup B|}.$$
#
# 下例刻意加入由左至右變亮的背景。局部方法在這組固定參數下表現較好；這不代表它在所有影像上都會勝出。

# %%
def binary_iou(prediction, truth):
    """計算兩個二值遮罩的 intersection over union。"""
    prediction = np.asarray(prediction, dtype=bool)
    truth = np.asarray(truth, dtype=bool)
    union = np.logical_or(prediction, truth).sum()
    if union == 0:
        return 1.0
    return np.logical_and(prediction, truth).sum() / union


height, width = 192, 256
yy_threshold, xx_threshold = np.mgrid[:height, :width]
threshold_truth = (
    (yy_threshold - 62) ** 2 + (xx_threshold - 55) ** 2 < 25**2
) | (
    (yy_threshold - 130) ** 2 + (xx_threshold - 140) ** 2 < 30**2
) | (
    (np.abs(yy_threshold - 62) < 18) & (np.abs(xx_threshold - 205) < 22)
)

rng_threshold = np.random.default_rng(2027)
illumination = 0.12 + 0.48 * xx_threshold / (width - 1)
uneven_objects = np.clip(
    illumination
    + 0.25 * threshold_truth
    + rng_threshold.normal(0, 0.025, threshold_truth.shape),
    0,
    1,
)

fixed_mask = uneven_objects > 0.52
otsu_value = ski.filters.threshold_otsu(uneven_objects)
otsu_mask = uneven_objects > otsu_value
local_surface = ski.filters.threshold_local(
    uneven_objects, block_size=91, method="gaussian", offset=-0.03
)
local_mask = uneven_objects > local_surface

threshold_scores = {
    "fixed": binary_iou(fixed_mask, threshold_truth),
    "Otsu": binary_iou(otsu_mask, threshold_truth),
    "local": binary_iou(local_mask, threshold_truth),
}
assert threshold_scores["local"] > threshold_scores["fixed"]
assert threshold_scores["local"] > threshold_scores["Otsu"]
print("segmentation IoU:", {k: round(v, 3) for k, v in threshold_scores.items()})

imshow_all(
    uneven_objects,
    threshold_truth,
    fixed_mask,
    otsu_mask,
    local_mask,
    titles=["uneven input", "known truth", "fixed", "Otsu", "local"],
    size=3,
)


# %% [markdown]
# ## 形態學：用結構元素描述鄰域
#
# 形態學運算以 structuring element（scikit-image 稱為 footprint）定義鄰域形狀。侵蝕要求 footprint 覆蓋處都屬於前景，因此會縮小前景；膨脹只要鄰域碰到前景就擴張。開運算是先侵蝕再膨脹，常用來移除比 footprint 小的亮點；閉運算是先膨脹再侵蝕，可接合窄縫或填平小缺口。圓盤、菱形與矩形對方向的偏好不同，footprint 必須配合要保留的結構。詳細定義見 [scikit-image morphology 文件](https://scikit-image.org/docs/stable/api/skimage.morphology.html)。

# %%
disk_footprint = ski.morphology.disk(4)
diamond_footprint = ski.morphology.diamond(4)
rectangle_footprint = ski.morphology.footprint_rectangle((7, 11))
imshow_all(
    disk_footprint,
    diamond_footprint,
    rectangle_footprint,
    titles=["disk", "diamond", "rectangle"],
    size=2.7,
)

morphology_demo = np.zeros((96, 128), dtype=bool)
morphology_demo[24:72, 35:93] = True
morphology_demo[43:53, 58:68] = False
morphology_demo[12:15, 14:17] = True
morphology_demo[79:82, 110:113] = True
small_disk = ski.morphology.disk(3)
eroded = ski.morphology.erosion(morphology_demo, small_disk)
dilated = ski.morphology.dilation(morphology_demo, small_disk)
opened = ski.morphology.opening(morphology_demo, small_disk)
closed = ski.morphology.closing(morphology_demo, small_disk)
assert eroded.sum() < morphology_demo.sum() < dilated.sum()

imshow_all(
    morphology_demo,
    eroded,
    dilated,
    opened,
    closed,
    titles=["input", "erosion", "dilation", "opening", "closing"],
    size=2.8,
)


# %% [markdown]
# ## Connected components、孔洞清理與 `regionprops`
#
# 二值遮罩仍只記錄前景與背景。Connected-component labeling 會依指定連通性替每個物件編號；`regionprops()` 再從 label image 計算面積、質心、bounding box、eccentricity 等幾何量。`remove_small_objects()` 與 `remove_small_holes()` 的面積參數都以像素計，換成另一個 pixel size 時，代表的物理尺度也會改變。

# %%
mask_with_defects = threshold_truth.copy()
rng_morphology = np.random.default_rng(2028)
noise_rows = rng_morphology.integers(0, height, size=180)
noise_cols = rng_morphology.integers(0, width, size=180)
mask_with_defects[noise_rows, noise_cols] = True
mask_with_defects[(yy_threshold - 62) ** 2 + (xx_threshold - 55) ** 2 < 5**2] = False

without_small_objects = ski.morphology.remove_small_objects(
    mask_with_defects, max_size=250
)
clean_mask = ski.morphology.remove_small_holes(
    without_small_objects, max_size=150
)
component_labels = ski.measure.label(clean_mask, connectivity=2)
component_properties = ski.measure.regionprops(component_labels)

assert len(component_properties) == 3
assert clean_mask[62, 55]
for prop in component_properties:
    row, col = prop.centroid
    print(
        f"label {prop.label}: area={prop.area:.0f} px, "
        f"centroid=({row:.1f}, {col:.1f}), eccentricity={prop.eccentricity:.3f}"
    )

imshow_all(
    mask_with_defects,
    without_small_objects,
    clean_mask,
    component_labels,
    titles=["defects", "remove small objects", "fill small holes", "labels"],
    size=3.2,
)


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
# ## Gaussian pyramid：同一影像的多個尺度
#
# Gaussian pyramid 在每次縮小前先做低通，降低 decimation 造成的 aliasing。較粗層級保留大尺度輪廓，細節與像素數則逐層減少。這是 scale-space 與多尺度偵測的入門；Laplacian pyramid 與小波會在第 4 章展開。若物件在 micrograph 中的直徑範圍很廣，多尺度搜尋可以減少單一模板尺寸造成的遺漏，但每一層的座標都要正確換回原影像。

# %%
gaussian_pyramid = tuple(
    ski.transform.pyramid_gaussian(
        search_image, downscale=2, max_layer=3, preserve_range=True
    )
)
assert len(gaussian_pyramid) == 4
assert all(
    coarse.shape[0] < fine.shape[0] and coarse.shape[1] < fine.shape[1]
    for fine, coarse in zip(gaussian_pyramid, gaussian_pyramid[1:])
)

fig, axes = plt.subplots(1, len(gaussian_pyramid), figsize=(13, 3.5))
for level, (ax, pyramid_image) in enumerate(zip(axes, gaussian_pyramid)):
    ax.imshow(pyramid_image)
    ax.set_title(f"level {level}\n{pyramid_image.shape}")
    ax.axis("off")
fig.tight_layout()


# %% [markdown]
# ## HOG：從梯度到固定長度的特徵向量
#
# Histogram of Oriented Gradients（HOG）先計算梯度，再於每個 cell 累積方向直方圖，最後把相鄰 cells 組成 block 做正規化。對 $H\times W$ 影像、`pixels_per_cell=(p_r,p_c)`、`cells_per_block=(b_r,b_c)` 與 $K$ 個方向 bins，完整 cell 的數量是
#
# $$n_r=\left\lfloor H/p_r\right\rfloor,\qquad
# n_c=\left\lfloor W/p_c\right\rfloor,$$
#
# 因而 feature length 為
#
# $$(n_r-b_r+1)(n_c-b_c+1)b_rb_cK.$$
#
# 下面使用 `(2,2)` cells per block，因此正規化確實跨越相鄰 cells。實作選項見 [scikit-image HOG 文件](https://scikit-image.org/docs/stable/api/skimage.feature.html#skimage.feature.hog)；HOG 是通用電腦視覺描述子，不是現代 cryo-EM picking 的完整模型。

# %%
def expected_hog_length(
    image_shape, orientations, pixels_per_cell, cells_per_block
):
    """依 cell 與 block 幾何計算灰階 HOG 特徵長度。"""
    cell_rows = image_shape[0] // pixels_per_cell[0]
    cell_cols = image_shape[1] // pixels_per_cell[1]
    block_rows = cell_rows - cells_per_block[0] + 1
    block_cols = cell_cols - cells_per_block[1] + 1
    return (
        block_rows
        * block_cols
        * cells_per_block[0]
        * cells_per_block[1]
        * orientations
    )


hog_orientations = 9
hog_pixels_per_cell = (16, 16)
hog_cells_per_block = (2, 2)
features, hog_image = ski.feature.hog(
    search_image,
    orientations=hog_orientations,
    pixels_per_cell=hog_pixels_per_cell,
    cells_per_block=hog_cells_per_block,
    visualize=True,
)
calculated_hog_length = expected_hog_length(
    search_image.shape,
    hog_orientations,
    hog_pixels_per_cell,
    hog_cells_per_block,
)
assert features.ndim == 1
assert features.size == calculated_hog_length
print("HOG feature length:", features.size)
imshow_all(search_image, hog_image, titles=["input", "HOG visualization"])


# %% [markdown]
# ## 與 cryo-EM 的連結
#
# 單粒子 cryo-EM 會把本章的工具放進不同階段，而不是把它們當成同一種「去雜訊」。粒子挑選可能使用 band-pass、whitening 或 local normalization 改善候選偵測條件；Gaussian 與 median filter 在本章主要用來比較濾波定義，median filter 只適合已確認的脈衝型污染。任何非線性濾波都可能改變微弱的高解析訊號；若處理結果要進入三維重建，必須先以獨立資料與頻域指標確認沒有引入偏差。
#
# Edge detector 與 thresholding 可用來建立碳膜邊緣、厚冰或污染區域的初始遮罩，再以 opening、closing 與 connected components 整理碎片。這類遮罩的 footprint 與面積門檻都對應物理尺度；micrograph 經過 binning 後，參數也要跟著換算。
#
# 模板匹配與 blob detection 可以產生候選座標，但候選不等於真實粒子。模板偏差、偏好取向、污染與冰層厚度都會影響挑選結果的分布。Gaussian pyramid 能支援多尺度搜尋，卻不能消除模板覆蓋不足。形態學與 watershed 適合清理污染遮罩或示範 instance segmentation，無法把物理上重疊的粒子投影分解成各自的訊號。
#
# ## 理解檢查
#
# 1. 為什麼同一個 mean kernel 在 `constant` 與 `reflect` 邊界會得到不同結果？
# 2. Gaussian noise 與 impulse noise 的比較中，為什麼不能只看哪張結果最平滑？
# 3. Canny 比直接對 gradient magnitude 切閾值多了哪些步驟？
# 4. 全域 Otsu 在不均勻背景下為何可能失效？局部閾值又新增了哪些參數？
# 5. Opening 與 `remove_small_objects()` 都能移除小區域，兩者的判斷依據有何不同？
# 6. Local normalization 與 whitening 分別校正哪一種不均勻？
# 7. NCC 分數高，為什麼仍不能斷言該位置是真實粒子？
# 8. Watershed 若只有「背景 marker 1」與「前景 marker 2」，為什麼不能分開兩個相鄰物件？
# 9. HOG 改用 `(3,3)` cells per block 時，特徵維度與邊界可用 block 數會如何變化？
