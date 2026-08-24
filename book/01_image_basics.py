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
# # 影像處理基礎：從像素到幾何變換
#
# 數位影像是陣列，但光知道陣列裡的數字還不夠。檔案格式、資料型別、色彩模型、取樣間距與座標變換，都會改變我們對這些數字的解讀。本章從 NumPy 陣列出發，一路處理影像讀寫、色彩與對比，最後進入取樣和幾何變換。這些是一般影像處理的基礎，也是正確讀取 cryo-EM 資料的前提。{cite}`szeliski2022,forsyth2012`
#
# ```{admonition} 學習目標
# :class: important
#
# - 使用 `(row, column)` 索引影像，並正確轉成 `(x, y)` 座標。
# - 以 ROI 與布林遮罩選取像素，並辨識通道維度。
# - 分清 dtype 的可表示範圍與資料實際採用的數值慣例。
# - 比較 PNG、TIFF 與 JPEG 的重要差異，並處理 RGB、RGBA、HSV、灰階與 BGR。
# - 使用 gamma、百分位對比伸展、直方圖等化與 CLAHE，同時說明它們不會產生新資訊。
# - 從 pixel size、取樣率與 Nyquist frequency 判斷 aliasing。
# - 讀取 MRC 的陣列與 header，不把檔案格式誤當成物理單位。
# - 比較 Euclidean、Similarity、Affine 與 Projective 變換，並說明插值、輸出畫布、邊界與反向座標映射。
# ```

# %%
from pathlib import Path

import matplotlib.pyplot as plt
import cv2
import mrcfile
import numpy as np
import skimage as ski

plt.rcParams["image.cmap"] = "gray"
plt.rcParams["figure.figsize"] = (6, 5)


def imshow_all(*images, titles=None, size=4, **kwargs):
    """並排顯示影像；不在此函式偷偷改變數值範圍。"""
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
# ## 影像是帶有慣例的陣列
#
# 灰階影像通常是 `(row, column)` 的二維陣列；彩色影像再加一個通道（channel）維度。`shape` 只告訴我們有多少格，沒有告訴我們一格是幾 Å，也沒有說數值是電子計數、正規化強度或其他量。

# %%
coins = ski.data.coins()
print("shape:", coins.shape, "dtype:", coins.dtype)
print("min/max/mean:", coins.min(), coins.max(), coins.mean())
assert coins.ndim == 2


# %% [markdown]
# ## 陣列索引、繪圖座標與物理座標
#
# `image[row, column]` 的第一個索引沿畫面向下，第二個索引沿畫面向右。若把欄當成 $x$、列當成 $y$，陣列索引 `[r, c]` 對應座標 `(x=c, y=r)`。
#
# `matplotlib.pyplot.imshow()` 預設 `origin="upper"`，所以 `(0, 0)` 和陣列一樣畫在左上角。只有設定 `origin="lower"` 時，顯示座標的垂直方向才會反轉。一般笛卡兒座標常把原點畫在左下；那是另一套慣例。

# %%
demo = np.zeros((5, 8), dtype=float)
demo[1, 6] = 1.0

fig, axes = plt.subplots(1, 2, figsize=(9, 3))
axes[0].imshow(demo, origin="upper")
axes[0].set_title("imshow default: origin='upper'")
axes[1].imshow(demo, origin="lower")
axes[1].set_title("origin='lower'")
fig.tight_layout()

row, col = np.unravel_index(np.argmax(demo), demo.shape)
assert (row, col) == (1, 6)
print(f"array index [row, col] = [{row}, {col}]; plot coordinate (x, y) = ({col}, {row})")


# %% [markdown]
# 若像素大小是 $p$ Å/pixel，忽略額外的原點偏移時，像素中心可寫成 $(x,y)=(cp,rp)$ Å。真實 MRC／STAR 處理流程可能另有 origin、binning 與 crop offset；換座標系時必須一起追蹤。

# %%
pixel_size_A = 1.2
x_A, y_A = col * pixel_size_A, row * pixel_size_A
print(f"physical coordinate: ({x_A:.1f}, {y_A:.1f}) Å")


# %% [markdown]
# ## ROI、切片與布林遮罩
#
# 矩形感興趣區域（region of interest, ROI）可以直接用切片取出。不規則區域則適合布林遮罩：遮罩和影像的前兩個維度相同，`True` 表示選中的像素。對灰階影像使用布林索引會得到一維像素集合；要保留原本形狀，可用 `np.where()` 將遮罩外的值改掉。{cite}`szeliski2022`

# %%
roi = coins[35:135, 90:190]
rows, cols = np.ogrid[: coins.shape[0], : coins.shape[1]]
center = (105, 155)
radius = 48
circular_mask = (rows - center[0]) ** 2 + (cols - center[1]) ** 2 <= radius**2
masked_coins = np.where(circular_mask, coins, 0)

assert roi.shape == (100, 100)
assert circular_mask.dtype == bool
assert np.count_nonzero(masked_coins[~circular_mask]) == 0
imshow_all(coins, roi, masked_coins, titles=["original", "rectangular ROI", "boolean mask"])

# %% [markdown]
# ```{admonition} ROI 的座標不會自動跟著走
# :class: caution
#
# `coins[35:135, 90:190]` 取出後，新陣列的左上角索引是 `(0, 0)`，但它在原影像中的位置是 `(35, 90)`。若 ROI 要與其他座標對應，必須另外記錄這個 offset。
# ```


# %% [markdown]
# ## dtype 與數值範圍是兩件事
#
# 整數 dtype 有固定的可表示範圍，例如 `uint8` 是 0–255。浮點 dtype 只規定精度與可表示範圍，**不保證**影像落在 `[0, 1]` 或 `[-1, 1]`。許多 `scikit-image` 函式以這些區間作為 float 影像慣例，但資料本身仍可能是電子計數、z-score 或任意實數。
#
# `img_as_float()` 會把整數依 dtype 範圍縮放；輸入若已是 float，通常保留原值。因此 `.astype(float)` 與 `img_as_float()` 對整數的結果不同，對 float 則都不會自動替你正規化。

# %%
integer_image = np.array([0, 64, 255], dtype=np.uint8)
float_outside_unit = np.array([-2.0, 0.5, 3.0], dtype=np.float32)

print("astype(float):", integer_image.astype(float))
print("img_as_float(uint8):", ski.util.img_as_float(integer_image))
print("img_as_float(float):", ski.util.img_as_float(float_outside_unit))

assert np.allclose(ski.util.img_as_float(integer_image), integer_image / 255)
assert np.array_equal(ski.util.img_as_float(float_outside_unit), float_outside_unit)

# %% [markdown]
# ```{admonition} 轉型前先問
# :class: caution
#
# 這個函式期待哪個數值區間？負值有沒有意義？輸出要拿來顯示，還是進入定量重建？若答案不清楚，不要先 clip 到 `[0,1]`；clip 會永久丟掉資訊。
# ```


# %% [markdown]
# ## PNG、TIFF 與 JPEG 讀寫
#
# `ski.io.imread()` 會把影像讀成 NumPy 陣列，`ski.io.imsave()` 則將陣列寫回檔案。副檔名會決定編碼器、是否採用有損壓縮，以及可支援的 bit depth 與中繼資料（metadata）。{cite}`howse2020`
#
# | 格式 | 壓縮與常見用途 | 定量資料的注意事項 |
# |---|---|---|
# | PNG | 無損；網頁圖、遮罩、標註 | 常見 8/16-bit 整數；不適合儲存任意浮點實數 |
# | TIFF | 可無損也可壓縮；顯微鏡與出版 | 支援多頁、多種 bit depth 與較豐富的中繼資料，但要確認讀寫程式支援的子類型 |
# | JPEG | 有損；一般照片與網路預覽 | 重新讀回後像素通常不會逐點相同，不可當作定量分析的中間格式 |
#
# 下面讓同一張影像分別通過 PNG、TIFF 與 JPEG 的編碼器，再解碼回 NumPy 陣列。比較讀回的像素與檔案大小，就能直接看出無損與有損壓縮的差別。

# %%
def codec_round_trip(image, extension, params=None):
    """用 OpenCV codec 編碼後再解碼；輸入與輸出都是 NumPy 陣列。"""
    params = [] if params is None else params
    ok, encoded = cv2.imencode(extension, image, params)
    if not ok:
        raise ValueError(f"cannot encode {extension}")
    decoded = cv2.imdecode(encoded, cv2.IMREAD_UNCHANGED)
    if decoded is None:
        raise ValueError(f"cannot decode {extension}")
    return decoded, encoded.nbytes


codec_sample = ski.data.astronaut()[80:208, 80:208]
sample_bgr = cv2.cvtColor(codec_sample, cv2.COLOR_RGB2BGR)
png_bgr, png_bytes = codec_round_trip(sample_bgr, ".png")
tiff_bgr, tiff_bytes = codec_round_trip(sample_bgr, ".tiff")
jpeg_bgr, jpeg_bytes = codec_round_trip(
    sample_bgr, ".jpg", [cv2.IMWRITE_JPEG_QUALITY, 85]
)

assert np.array_equal(png_bgr, sample_bgr)
assert np.array_equal(tiff_bgr, sample_bgr)
assert not np.array_equal(jpeg_bgr, sample_bgr)
jpeg_mae = np.mean(np.abs(jpeg_bgr.astype(float) - sample_bgr.astype(float)))
print(f"PNG/TIFF exact round trip; JPEG mean absolute error = {jpeg_mae:.2f}")
print(f"encoded bytes: PNG={png_bytes}, TIFF={tiff_bytes}, JPEG={jpeg_bytes}")

# %%
gradient16 = np.linspace(0, 65535, 128, dtype=np.uint16)[None, :].repeat(32, axis=0)
gradient16_back, _ = codec_round_trip(gradient16, ".png")
assert gradient16_back.dtype == np.uint16
assert np.array_equal(gradient16_back, gradient16)
print("16-bit PNG dtype:", gradient16_back.dtype, "range:", gradient16_back.min(), gradient16_back.max())

# %% [markdown]
# JPEG 的檔案往往較小，代價是區塊狀壓縮痕跡與高頻細節的改變。兩張影像視覺上可能很接近，像素值卻已經改變；上面的逐像素相等檢查與 mean absolute error 分別檢查這兩個層次。


# %% [markdown]
# (sampling)=
# ## 取樣、Nyquist frequency 與 aliasing
#
# 若 pixel size 是 $p$ Å/pixel，取樣率為 $1/p$ cycles/Å，而 Nyquist frequency 是
#
# $$f_N=\frac{1}{2p}\quad\text{cycles/Å}.$$
#
# 相對應的 Nyquist resolution 是 $2p$ Å。這是取樣網格能表示的上限，不保證資料真的含有該解析度的可靠訊號。高於 Nyquist 的頻率會折回較低頻率，形成 aliasing；降採樣前先低通濾波，就是為了移除會折回的頻率。{cite}`szeliski2022`

# %%
def nyquist_frequency(pixel_size):
    """由正的 pixel size 計算 Nyquist frequency。"""
    if pixel_size <= 0:
        raise ValueError("pixel_size must be positive")
    return 1.0 / (2.0 * pixel_size)


assert np.isclose(nyquist_frequency(1.2), 1 / 2.4)
print("pixel size 1.2 Å/pixel -> Nyquist resolution 2.4 Å")

# %%
n = 64
index = np.arange(n)
# 0.60 cycles/pixel 超過 0.5 Nyquist，離散樣本與 0.40 cycles/pixel 無法區分（相位符號除外）。
above_nyquist = np.cos(2 * np.pi * 0.60 * index)
aliased = np.cos(2 * np.pi * 0.40 * index)
assert np.allclose(above_nyquist, aliased)

fig, ax = plt.subplots(figsize=(8, 3))
ax.plot(index[:20], above_nyquist[:20], "o-", label="sampled 0.60 cycles/pixel")
ax.plot(index[:20], aliased[:20], "x--", label="0.40 cycles/pixel alias")
ax.set_xlabel("sample index")
ax.legend()
fig.tight_layout()


# %% [markdown]
# ## MRC：陣列之外還有 header
#
# MRC 常用來儲存 micrograph、粒子堆疊與 3D 密度圖。Header 會記錄陣列尺寸、儲存 mode、cell dimensions 等欄位；常用套件則將 voxel size 提供為較容易讀取的屬性。Header 的值可能缺漏或被舊程式寫錯，分析前仍要和擷取與處理的中繼資料交叉核對。

# %%
mrc_path = Path("data/70S_Conform1.mrc")
with mrcfile.mmap(mrc_path, permissive=True, mode="r") as volume_mrc:
    volume_shape = volume_mrc.data.shape
    volume_dtype = volume_mrc.data.dtype
    voxel_size = tuple(float(v) for v in volume_mrc.voxel_size.tolist())
    mrc_mode = int(volume_mrc.header.mode)

print("MRC shape:", volume_shape)
print("MRC dtype/mode:", volume_dtype, mrc_mode)
print("voxel size [Å]:", voxel_size)
assert len(volume_shape) == 3
assert all(length > 0 for length in volume_shape)

# %% [markdown]
# ```{dropdown} MRC 的軸順序為什麼特別容易出錯？
#
# `mrcfile` 回傳的 3D NumPy 陣列通常以 `(section, row, column)` 存取，也就是常寫的 `(z, y, x)`；MRC header 另有 `mapc/mapr/maps` 描述檔案欄、列、切片對應哪個空間軸。只看陣列 shape 不足以判定生物結構的方向。
# ```


# %% [markdown]
# ## RGB、RGBA、HSV、灰階與 BGR
#
# RGB 用紅、綠、藍三個通道表示色彩；RGBA 再加上 alpha 通道表示不透明度。這裡先把 alpha 與白色背景合成，再轉成灰階。`rgb2gray()` 依亮度加權三個通道；目前 `scikit-image` 使用的係數是 0.2125、0.7154、0.0721。{cite}`szeliski2022,howse2020`

# %%
rgba = ski.io.imread("images/blastoise.png")
rgb = ski.color.rgba2rgb(rgba)
gray = ski.color.rgb2gray(rgb)
gray_manual = rgb @ np.array([0.2125, 0.7154, 0.0721])
assert np.allclose(gray, gray_manual)
assert rgba.shape[-1] == 4 and rgb.shape[-1] == 3

fig, axes = plt.subplots(1, 5, figsize=(15, 3))
channel_names = ["red", "green", "blue", "alpha"]
for index_channel, (axis, name) in enumerate(zip(axes[:4], channel_names)):
    axis.imshow(rgba[..., index_channel], vmin=0, vmax=255)
    axis.set_title(name)
    axis.axis("off")
axes[4].imshow(rgb)
axes[4].set_title("RGBA on white")
axes[4].axis("off")
fig.tight_layout()

# %% [markdown]
# HSV 把色相（hue）、飽和度（saturation）與明度（value）分開，方便依顏色而非單純依亮度選取區域。但當飽和度接近零時，色相對微小雜訊很敏感，不適合單獨拿來分類。

# %%
hsv = ski.color.rgb2hsv(rgb)
rgb_from_hsv = ski.color.hsv2rgb(hsv)
assert np.allclose(rgb_from_hsv, rgb, atol=1e-12)
imshow_all(
    rgb,
    gray,
    hsv[..., 0],
    hsv[..., 1],
    hsv[..., 2],
    titles=["RGB", "grayscale", "hue", "saturation", "value"],
)

# %% [markdown]
# OpenCV 的 `imread()` 預設回傳 BGR，而 `matplotlib` 與 `scikit-image` 使用 RGB。這個差異不會改變 row 與 column，但會交換紅、藍通道。使用 `cv2.cvtColor()` 明確轉換，比在程式中預設通道順序安全。{cite}`howse2020`

# %%
rgb_uint8 = ski.util.img_as_ubyte(rgb)
bgr = cv2.cvtColor(rgb_uint8, cv2.COLOR_RGB2BGR)
rgb_restored = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
assert np.array_equal(rgb_restored, rgb_uint8)
imshow_all(rgb_uint8, bgr, rgb_restored, titles=["RGB", "BGR shown as RGB", "converted back"])


# %% [markdown]
# ## 強度變換、直方圖與局部對比
#
# 直方圖計數各強度區間出現的次數，累積分布函數（cumulative distribution function, CDF）則回答「不大於某強度的像素佔多少」。兩者都不保留像素的空間位置，因此直方圖相同的兩張影像，結構可能完全不同。{cite}`szeliski2022,forsyth2012`
#
# 對比伸展將指定強度區間線性映射到輸出範圍；百分位數可避免少數離群值主導整個映射。Gamma 變換 $s=r^\gamma$ 是非線性的：$\gamma<1$ 抬高暗部，$\gamma>1$ 壓低暗部。全域直方圖等化改寫整張影像的 CDF；CLAHE 在局部區塊處理，並以 clip limit 限制雜訊被放大的程度。

# %%
def histogram_and_cdf(image, bins=256):
    """回傳正規化直方圖、bin centers 與 CDF。"""
    histogram, centers = ski.exposure.histogram(image, nbins=bins, normalize=True)
    cdf = np.cumsum(histogram)
    cdf /= cdf[-1]
    return histogram, centers, cdf


moon = ski.util.img_as_float(ski.data.moon())
p2, p98 = np.percentile(moon, (2, 98))
moon_stretched = ski.exposure.rescale_intensity(moon, in_range=(p2, p98))
moon_gamma = ski.exposure.adjust_gamma(moon, gamma=0.6)
moon_equalized = ski.exposure.equalize_hist(moon)
moon_clahe = ski.exposure.equalize_adapthist(moon, clip_limit=0.02)

histogram, bin_centers, cdf = histogram_and_cdf(moon)
assert np.isclose(histogram.sum(), 1.0)
assert np.all(np.diff(cdf) >= 0) and np.isclose(cdf[-1], 1.0)
assert moon_stretched.min() == 0 and moon_stretched.max() == 1

# %%
contrast_images = [moon, moon_stretched, moon_gamma, moon_equalized, moon_clahe]
contrast_titles = ["original", "2–98% stretch", "gamma = 0.6", "global equalization", "CLAHE"]
fig, axes = plt.subplots(2, len(contrast_images), figsize=(16, 6))
for column, (contrast_image, title) in enumerate(zip(contrast_images, contrast_titles)):
    axes[0, column].imshow(contrast_image, vmin=0, vmax=1)
    axes[0, column].set_title(title)
    axes[0, column].axis("off")
    hist, centers, local_cdf = histogram_and_cdf(contrast_image)
    axes[1, column].plot(centers, hist, color="black", label="histogram")
    cdf_axis = axes[1, column].twinx()
    cdf_axis.plot(centers, local_cdf, color="tab:red", label="CDF")
    axes[1, column].set_xlim(0, 1)
    axes[1, column].set_yticks([])
    cdf_axis.set_yticks([])
fig.tight_layout()

# %% [markdown]
# ```{admonition} 顯示清晰度與資訊量
# :class: caution
#
# Gamma、對比伸展和等化可讓特徵更容易看見，但不會恢復原本沒有取樣到的頻率，也不會提高實驗本身的訊雜比。在 cryo-EM 中，顯示用的強度映射應和後續定量處理分開；否則很容易把「更好看」誤解為「解析度更高」。
# ```


# %% [markdown]
# ## 幾何變換的階層
#
# 平面上的點可以用齊次座標 $(x,y,1)^\mathsf{T}$ 表示，讓平移也能寫進 $3\times3$ 矩陣。從 Euclidean 到 Projective，模型可表示的形變愈多，自由度也愈高。可以擬合更複雜的形變，不表示就應該無條件使用它：控制點不足或座標有雜訊時，額外自由度反而容易擬合誤差。{cite}`szeliski2022,forsyth2012`
#
# | 變換 | 2D 自由度 | 保留的幾何性質 | 最少點對（一般位置） |
# |---|---:|---|---:|
# | Euclidean（剛體） | 3 | 距離、角度、平行性 | 2 |
# | Similarity（相似） | 4 | 角度、形狀、平行性；距離可統一縮放 | 2 |
# | Affine（仿射） | 6 | 直線、平行性、同一線上的長度比 | 3 |
# | Projective（投影） | 8 | 直線與交比；平行線可交於消失點 | 4 |
#
# 自由度不包含齊次矩陣整體乘上一個非零常數的重複表示。例如 Projective 矩陣雖有 9 個數，整體尺度不改變映射，所以只有 8 個自由度。
#
# ### `warp()` 問的是反向映射
#
# 幾何變換常以 forward map 表示「輸入點會移到哪裡」。但要產生每一個輸出像素，`warp()` 必須回到輸入影像查值，因此參數 `inverse_map` 的方向是 **output → input**。這裡的 inverse 指反向座標映射；把矩陣乘以 $-1$ 會得到錯誤的變換。

# %%
image = ski.util.img_as_float(ski.data.camera())[::2, ::2]
forward = ski.transform.EuclideanTransform(rotation=np.deg2rad(8), translation=(15, -8))
warped = ski.transform.warp(image, inverse_map=forward.inverse)

similarity = ski.transform.SimilarityTransform(
    scale=0.82, rotation=np.deg2rad(-6), translation=(18, 14)
)
affine = ski.transform.AffineTransform(scale=(0.95, 0.8), shear=np.deg2rad(12))
warped_similarity = ski.transform.warp(image, inverse_map=similarity.inverse)
warped_affine = ski.transform.warp(image, inverse_map=affine.inverse)

# 若已知輸出矩形座標 `output_xy` 對應輸入影像中的 `input_xy`，直接估計的就是 output → input。
output_xy = np.array([[0, 0], [0, 80], [160, 80], [160, 0]])
input_xy = np.array([[35, 10], [15, 100], [175, 90], [155, 5]])
output_to_input = ski.transform.ProjectiveTransform.from_estimate(output_xy, input_xy)
assert output_to_input is not None
rectified = ski.transform.warp(image, inverse_map=output_to_input, output_shape=(80, 160))

imshow_all(
    image,
    warped,
    warped_similarity,
    warped_affine,
    rectified,
    titles=["input", "Euclidean", "Similarity", "Affine", "Projective rectification"],
)

# %% [markdown]
# ### 插值、輸出畫布與邊界
#
# 反向映射通常不會剛好落在整數像素上，因此需要插值。Nearest-neighbor 只取最近樣本，適合類別標籤；bilinear 使用 $2\times2$ 鄰域，速度與平滑度折衷良好；bicubic 使用 $4\times4$ 鄰域，通常更平滑，卻可能在銳利邊緣附近產生 overshoot。對 segmentation label 使用高階插值會創造本來不存在的類別值，這時應用 nearest-neighbor。{cite}`szeliski2022`
#
# `output_shape` 決定輸出畫布。畫布太小會裁掉變換後的內容，畫布變大也不會自動調整座標原點。反向映射超出輸入範圍時，`mode` 決定要使用常數、邊界值、反射或環繞；這項選擇會直接改變邊緣附近的數值。

# %%
small_patch = image[70:125, 85:140]
subpixel_shift = ski.transform.EuclideanTransform(translation=(0.45, 0.45))
nearest = ski.transform.warp(
    small_patch, subpixel_shift.inverse, order=0, mode="constant", cval=0
)
bilinear = ski.transform.warp(
    small_patch, subpixel_shift.inverse, order=1, mode="constant", cval=0
)
bicubic = ski.transform.warp(
    small_patch, subpixel_shift.inverse, order=3, mode="constant", cval=0
)
wide_canvas = ski.transform.warp(
    image, forward.inverse, output_shape=(320, 320), mode="constant", cval=0
)
reflected_boundary = ski.transform.warp(
    small_patch, subpixel_shift.inverse, order=1, mode="reflect"
)
assert not np.allclose(nearest, bilinear)
assert not np.allclose(bilinear, bicubic)
assert wide_canvas.shape == (320, 320)
assert not np.allclose(reflected_boundary, bilinear)
imshow_all(small_patch, nearest, bilinear, bicubic, titles=["input", "nearest", "bilinear", "bicubic"])

fig, axes = plt.subplots(1, 3, figsize=(10, 3))
axes[0].imshow(wide_canvas)
axes[0].set_title("larger output canvas")
axes[1].imshow(bilinear)
axes[1].set_title("constant boundary")
axes[2].imshow(reflected_boundary)
axes[2].set_title("reflected boundary")
for axis in axes:
    axis.axis("off")
fig.tight_layout()

# %% [markdown]
# ### 縮小前的 anti-aliasing
#
# 影像縮小也是重新取樣。新網格的 Nyquist frequency 較低，所以縮小前應先低通濾波。只改輸出 shape 而沒有 anti-aliasing，高頻棋盤格會折成低頻假紋理。

# %%
grid_row, grid_col = np.indices((192, 192))
fine_checkerboard = ((grid_row // 3 + grid_col // 3) % 2).astype(float)
downsampled_no_aa = ski.transform.resize(
    fine_checkerboard,
    (48, 48),
    order=1,
    anti_aliasing=False,
    preserve_range=True,
)
downsampled_aa = ski.transform.resize(
    fine_checkerboard,
    (48, 48),
    order=1,
    anti_aliasing=True,
    preserve_range=True,
)
assert downsampled_aa.std() < downsampled_no_aa.std() / 5
imshow_all(
    fine_checkerboard,
    downsampled_no_aa,
    downsampled_aa,
    titles=["fine grid", "resize without anti-aliasing", "resize with anti-aliasing"],
)


# %% [markdown]
# ## 與 cryo-EM 的連結
#
# 一套單粒子分析資料會經過幾個層次：偵測器先記錄電影式影格（movie frames）；完成位移校正後，將影格加權或未加權加總，才得到 micrograph；接著從 micrograph 擷取粒子影像，後續對位、分類與重建再使用這些影像。名稱不同，陣列的意義也不同。
#
# Pixel size、crop origin、binning 與座標慣例必須隨資料一起記錄。例如 2× binning 後，陣列每個方向的像素數約減半，pixel size 則變成原來的兩倍；若只改 shape 卻忘了更新 pixel size，所有以 Å 表示的距離與頻率都會跟著錯。
#
# MRC 是一種檔案容器。同樣的格式可以儲存 2D micrograph、粒子堆疊或 3D 密度圖，資料的含義要連同 shape、header 與外部中繼資料解讀。與一般影像相同，幾何重新取樣會改變頻率內容；對粒子影像做旋轉或位移時，插值方法、邊界模式與輸出畫布仍然會影響後續對位。
#
# ```{admonition} 換成低劑量影像時
# :class: note
#
# 本章先用灰階照片練習陣列操作。讀取低劑量 cryo-EM 影像時，還要把 shot noise、gain correction、motion、ice 與 DQE 放進解讀；直接電子偵測器的響應也無法用一般相機常見的 `[0,255]` 範圍完整描述。遇到負值或大於 1 的浮點數時，應先查清楚處理步驟與單位，再決定如何顯示或轉型。
# ```
#
# ## 延伸閱讀
#
# - 影像取樣、幾何變換與色彩表示可接著讀 {cite}`szeliski2022,forsyth2012`。
# - OpenCV 的色彩順序與檔案讀寫範例可參考 {cite}`howse2020`；函式參數則以目前安裝版本的官方文件為準。
#
# ## 理解檢查
#
# 1. 一個 `float32` 陣列的最大值是 12，為什麼還不能判定它是否經過正規化？
#
#    ```{dropdown} 參考答案
#    資料不足。`float32` 只規定數值的儲存方式與可表示範圍，正規化則是資料的數值慣例。最大值 12 可能來自原始量測值，也可能是某種自訂尺度的正規化結果。還需要查讀檔案說明、數值範圍與處理程式，才能判定 12 的意義。
#
#    例如資料若以電子計數表示，12 可以是正常量測值；若定義為零平均、單位變異數，12 又可能是極端值。
#    判讀時要一起查看單位、預期範圍、是否扣除平均，以及前一步做過哪些縮放或轉型。
#    ```
#
# 2. 灰階區域的 saturation 很低時，用 hue 做分類會遇到什麼問題？
#
#    ```{dropdown} 參考答案
#    Hue 描述色相，saturation 則描述顏色偏離灰階的程度。當 $R\approx G\approx B$ 時，通道最大值與最小值的差（chroma）接近零，色相的計算對很小的像素誤差都會很敏感。理想的完全灰階像素其 hue 甚至沒有明確的物理意義。
#
#    分類時可先設 saturation 下限。只有高於下限的像素才使用 hue，其餘像素可依 brightness、RGB 距離或其他特徵處理。問題來自色相在低飽和度時的不穩定性，與 HSV 轉換程式是否正確無關。
#    ```
#
# 3. 百分位對比伸展讓粒子顯示得更清楚時，能從圖上確定什麼？訊雜比與解析度又該如何判斷？
#
#    ```{dropdown} 參考答案
#    只能得到「在這組顯示參數下，人眼較容易看出強度差」。若將下百分位 $p_l$ 與上百分位 $p_h$ 間的像素線性映射至 $[0,1]$，中間區間可寫成
#
#    $$I'(x)=\frac{I(x)-p_l}{p_h-p_l}.$$
#
#    這個映射同時放大訊號與雜訊的局部差異，也會將區間外的數值截斷。它沒有自動增加樣品傳遞到偵測器的資訊。
#
#    訊雜比需要明確的 signal/noise 定義與估計法；解析度則要根據獨立資料的頻域指標，例如 half-map FSC。可以把原圖與顯示後的圖並排，再分別計算訊雜比與 FSC，避免用「看起來更清楚」代替這兩項數值。
#    ```
#
# 4. pixel size 為 1.5 Å/pixel 時，Nyquist frequency 與 Nyquist resolution 各是多少？
#
#    ```{dropdown} 參考答案
#    取樣間距為 $\Delta x=1.5$ Å/pixel，Nyquist frequency 是取樣頻率的一半：
#
#    $$f_{\mathrm{Nyq}}=\frac{1}{2\Delta x}=\frac{1}{3.0}\approx0.333\ \text{cycles/Å}.$$
#
#    對應的空間週期，也就是 Nyquist resolution，為
#
#    $$d_{\mathrm{Nyq}}=\frac{1}{f_{\mathrm{Nyq}}}=2\Delta x=3.0\ \text{Å}.$$
#
#    頻率的單位是 cycles/Å，解析度則以 Å 表示。數值較大的頻率對應較小的結構尺度。Nyquist resolution 是離散取樣允許的理論上限，實際可達解析度還會受訊雜比、CTF、motion 與重建流程影響。
#    ```
#
# 5. `warp(image, forward.inverse)` 為什麼要傳入 `inverse`？
#
#    ```{dropdown} 參考答案
#    `forward` 描述輸入座標如何移到輸出畫布。建立輸出影像時，程式逐一走訪輸出像素 $\mathbf{x}_{out}$，再透過
#
#    $$\mathbf{x}_{in}=T^{-1}\mathbf{x}_{out}$$
#
#    找到輸入影像的取樣位置。這種 output-to-input 的 pull sampling 可避免 forward mapping 在輸出畫布留下空洞。這裡的 `inverse` 只指座標映射方向，與像素值的反向變化或矩陣乘上 $-1$ 無關。Euclidean、Similarity、Affine 與 Projective 的保留性質可回到本節表格比較，不會改變 `warp()` 所需的映射方向。
#    ```
