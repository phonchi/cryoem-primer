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
# # 影像處理基礎：像素、取樣與 MRC
#
# 一張數位影像不只是「一堆數字」。要正確解讀它，至少要同時知道陣列形狀、資料型別、像素值的意義，以及每個像素對應的實際長度。本章先把這四件事釐清。
#
# ```{admonition} 學習目標
# :class: important
#
# - 使用 `(row, column)` 索引影像，並正確轉成 `(x, y)` 座標。
# - 分清 dtype 的可表示範圍與資料實際採用的數值慣例。
# - 從 pixel size、取樣率與 Nyquist frequency 判斷 aliasing。
# - 讀取 MRC 的陣列與 header，不把檔案格式誤當成物理單位。
# - 說明 `warp()` 為什麼需要 output → input 的反向座標映射。
# ```

# %%
from pathlib import Path

import matplotlib.pyplot as plt
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
# 灰階影像通常是 `(row, column)` 的二維陣列；彩色影像再加一個 channel 維度。`shape` 只告訴我們有多少格，沒有告訴我們一格是幾 Å，也沒有說數值是電子計數、正規化強度或其他量。

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
# 若像素大小是 $p$ Å/pixel，忽略額外的原點偏移時，像素中心可寫成 $(x,y)=(cp,rp)$ Å。真實 MRC／STAR workflow 可能另有 origin、binning 與 crop offset；換座標系時必須一起追蹤。

# %%
pixel_size_A = 1.2
x_A, y_A = col * pixel_size_A, row * pixel_size_A
print(f"physical coordinate: ({x_A:.1f}, {y_A:.1f}) Å")


# %% [markdown]
# ## dtype 不等於數值範圍慣例
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
# (sampling)=
# ## 取樣、Nyquist frequency 與 aliasing
#
# 若 pixel size 是 $p$ Å/pixel，取樣率為 $1/p$ cycles/Å，而 Nyquist frequency 是
#
# $$f_N=\frac{1}{2p}\quad\text{cycles/Å}.$$
#
# 相對應的 Nyquist resolution 是 $2p$ Å。這是取樣網格能表示的上限，不保證資料真的含有該解析度的可靠訊號。高於 Nyquist 的頻率會折回較低頻率，形成 aliasing；降採樣前先低通濾波，就是為了移除會折回的頻率。

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
# MRC 常用來存 micrograph、particle stack 與 3D density map。header 會記錄陣列尺寸、儲存 mode、cell dimensions 等欄位；常用套件把 voxel size 暴露成較易讀的屬性。header 的值可能缺漏或被舊程式寫錯，分析前仍要和 acquisition／processing metadata 交叉核對。

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
# ## RGB、灰階與顯示
#
# 一般照片適合用來練習 channel，但 cryo-EM 強度影像通常是灰階。`rgb2gray()` 採用版本所定義的亮度權重；目前 `scikit-image` 使用的係數是 0.2125、0.7154、0.0721，和常見 Rec. 709 四捨五入係數略有差異。

# %%
rgba = ski.io.imread("images/blastoise.png")
rgb = ski.color.rgba2rgb(rgba)
gray = ski.color.rgb2gray(rgb)
gray_manual = rgb @ np.array([0.2125, 0.7154, 0.0721])
assert np.allclose(gray, gray_manual)
imshow_all(rgb, gray, titles=["RGB", "grayscale"])


# %% [markdown]
# ## 幾何變換：`warp()` 問的是反向映射
#
# 幾何變換常以 forward map 表示「輸入點會移到哪裡」。但要產生每一個輸出像素，`warp()` 必須回到輸入影像查值，因此參數 `inverse_map` 的方向是 **output → input**。這裡的 inverse 是反向座標映射，不是把矩陣乘以 $-1$。

# %%
image = ski.util.img_as_float(ski.data.camera())[::2, ::2]
forward = ski.transform.EuclideanTransform(rotation=np.deg2rad(8), translation=(15, -8))
warped = ski.transform.warp(image, inverse_map=forward.inverse)

# 若已知輸出矩形座標 `output_xy` 對應輸入影像中的 `input_xy`，直接估計的就是 output → input。
output_xy = np.array([[0, 0], [0, 80], [160, 80], [160, 0]])
input_xy = np.array([[35, 10], [15, 100], [175, 90], [155, 5]])
output_to_input = ski.transform.ProjectiveTransform.from_estimate(output_xy, input_xy)
rectified = ski.transform.warp(image, inverse_map=output_to_input, output_shape=(80, 160))

imshow_all(image, warped, rectified, titles=["input", "Euclidean warp", "projective rectification"])


# %% [markdown]
# ## 與 cryo-EM 的連結
#
# 偵測器先輸出 movie frames；經 motion correction 的加權或未加權 sum 才通常稱為 micrograph。從 micrograph 擷取 particle images 後，pixel size、crop origin、binning 與座標慣例都要跟著資料走。MRC 只是一種容器：同樣的 MRC 格式可以裝 2D micrograph、particle stack 或 3D density map，真正含義來自 shape、header 與外部 metadata 的組合。
#
# ```{admonition} 主張範圍
# :class: caution
#
# 本章的灰階照片用來教陣列操作，不代表直接電子偵測器的完整 response。低劑量 cryo-EM 的雜訊還受 shot noise、gain correction、motion、ice 與 DQE 等因素影響，不能只用一般相機的 `[0,255]` 直覺解讀。
# ```
#
# ## 理解檢查
#
# 1. 一個 `float32` 陣列的最大值是 12，能否僅憑 dtype 判定它「沒有正規化」？
# 2. pixel size 為 1.5 Å/pixel 時，Nyquist frequency 與 Nyquist resolution 各是多少？
# 3. 為什麼 `warp(image, forward.inverse)` 中的 `inverse` 不是「把影像做反變換」的口語說法？
# 4. 只看到一個 MRC 的 shape 是 `(256, 256, 256)`，還缺哪些資訊才能解讀三個軸？
