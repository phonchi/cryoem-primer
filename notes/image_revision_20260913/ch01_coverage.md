# 第 1 章原始教材覆蓋紀錄

來源：`../01_image_processing_1_zh.ipynb`，171 cells。依原順序逐一恢復；原notebook未修改。`source-cells` 使用零起算。所有原程式圖片/幾何參數保留，舊安裝/plugin操作改為現有讀寫介面，兩處先使用後定義的變數合併到賦值後。

| 原cell | 類型 | 新source行 | 教學去向／理由 |
|---:|---|---:|---|
| 0 | code | 18 | 環境與helper；原單元：# @title 🛠️ Setup %pip install imread -qq import os import reques；保留helper功能並改用shared support；移除Colab安裝與網路下載 |
| 1 | markdown | 73 | # 影像處理基礎：從像素、色彩到幾何變換；原單元：# 📘 什麼是電腦視覺？ ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 2 | markdown | 83 | ## 什麼是電腦視覺？；原單元：<u>電腦視覺</u>是一個處理電腦如何從數位影像或影片中獲得高階理解的領域。從工程角度來看，*它尋求理解並自動化人類視覺系統能夠；原圖／參數／操作保留；必要語句修正，順序不變 |
| 3 | markdown | 89 | ## 影像表示法；原單元：# 📘 影像表示法；原圖／參數／操作保留；必要語句修正，順序不變 |
| 4 | markdown | 95 | ## 影像表示法；原單元：<u>點陣圖</u>在我們需要儲存類比影像時使用。類比影像會被**取樣與量化**。資料的強度（顏色）在空間中變化，透過取樣來獲取我；原圖／參數／操作保留；必要語句修正，順序不變 |
| 5 | markdown | 106 | ## 影像表示法；原單元：<center><img src="https://drive.google.com/uc?id=170NOZD-HtG4jt-P；原圖／參數／操作保留；必要語句修正，順序不變 |
| 6 | markdown | 111 | ## 影像表示法；原單元：<center><img src="https://drive.google.com/uc?id=170igILkg6b76ofZ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 7 | markdown | 122 | ## 影像是 NumPy 陣列；原單元：# 📘 影像是 `NumPy` 陣列；原圖／參數／操作保留；必要語句修正，順序不變 |
| 8 | markdown | 128 | ## 影像是 NumPy 陣列；原單元：`scikit-image` 是一個與 `NumPy` 陣列一起使用的影像處理 `Python` 套件。其他可用於影像處理並與 `；原圖／參數／操作保留；必要語句修正，順序不變 |
| 9 | markdown | 134 | ## 影像是 NumPy 陣列；原單元：此外，強大的影像處理函式庫也有 `Python` 綁定： * [OpenCV](https://docs.opencv.org/4；原圖／參數／操作保留；必要語句修正，順序不變 |
| 10 | markdown | 140 | ## 影像是 NumPy 陣列；原單元：在 ``scikit-image`` 中，影像使用標準的 `NumPy` 陣列來表示。這使得與科學 `Python` 生態系統中的；原圖／參數／操作保留；必要語句修正，順序不變 |
| 11 | markdown | 146 | ## 影像是 NumPy 陣列；原單元：／ 影像: ／ np.ndarray ／ ／-----------------／-------------------------；原圖／參數／操作保留；必要語句修正，順序不變 |
| 12 | markdown | 160 | ## 影像是 NumPy 陣列；原單元：讓我們看看如何建立一個灰階影像作為二維陣列：；原圖／參數／操作保留；必要語句修正，順序不變 |
| 13 | code | 166 | ## 影像是 NumPy 陣列；原單元：random_image.shape；shape查詢合併到cell14建立陣列後 |
| 14 | code | 170 | ## 影像是 NumPy 陣列；原單元：random_image = np.random.random([512, 512]) plt.imshow(random_ima；原圖／參數／操作保留；必要語句修正，順序不變 |
| 15 | markdown | 177 | ## 影像是 NumPy 陣列；原單元：大多數的 `skimage` 函式都集中在其子模組中。`skimage.data` 子模組提供一組可返回範例影像的函式，可協助使用；原圖／參數／操作保留；必要語句修正，順序不變 |
| 16 | code | 183 | ## 影像是 NumPy 陣列；原單元：coins = ski.data.coins() print('Type:', type(coins)) print('dtype；原圖／參數／操作保留；必要語句修正，順序不變 |
| 17 | markdown | 193 | ## 影像是 NumPy 陣列；原單元：當然，也可以使用 `skimage.io.imread()` 從影像檔案載入自己的圖片，並以 `NumPy` 陣列的形式讀取進來。；原圖／參數／操作保留；必要語句修正，順序不變 |
| 18 | code | 199 | ## 影像是 NumPy 陣列；原單元：download_from_pokemondb('https://img.pokemondb.net/sprites/red-bl；原圖／參數／操作保留；必要語句修正，順序不變 |
| 19 | markdown | 204 | ## 影像是 NumPy 陣列；原單元：取得影像的幾何資訊與像素數量：；原圖／參數／操作保留；必要語句修正，順序不變 |
| 20 | code | 210 | ## 影像是 NumPy 陣列；原單元：venusaur.shape, venusaur.size, venusaur.dtype；原圖／參數／操作保留；必要語句修正，順序不變 |
| 21 | markdown | 214 | ## 影像是 NumPy 陣列；原單元：取得影像強度值的統計資訊：；原圖／參數／操作保留；必要語句修正，順序不變 |
| 22 | code | 220 | ## 影像是 NumPy 陣列；原單元：venusaur.min(), venusaur.max(), venusaur.mean()；原圖／參數／操作保留；必要語句修正，順序不變 |
| 23 | markdown | 224 | ### 索引、座標與局部修改；原單元：## `NumPy` 索引（座標）；原圖／參數／操作保留；必要語句修正，順序不變 |
| 24 | markdown | 230 | ### 索引、座標與局部修改；原單元：可以使用 `NumPy` 索引來查看像素值，也可以用來修改像素值： ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 25 | code | 236 | ### 索引、座標與局部修改；原單元：# Get the value of the pixel at the 10th row and 20th column venu；原圖／參數／操作保留；必要語句修正，順序不變 |
| 26 | code | 241 | ### 索引、座標與局部修改；原單元：# Set to black the pixel at the 3rd row and 10th column venusaur[；原圖／參數／操作保留；必要語句修正，順序不變 |
| 27 | markdown | 246 | ### 索引、座標與局部修改；原單元：📌 請注意！在 `NumPy` 索引中，第一個維度（`venusaur.shape[0]`）對應到列，第二個維度（`venusau；原圖／參數／操作保留；必要語句修正，順序不變 |
| 28 | markdown | 252 | ### 索引、座標與局部修改；原單元：二維（2D）灰階影像（例如前面的 `venusaur`）是以列與欄來索引的，簡寫為 `(row, col)` 或 `(r, c)`；原圖／參數／操作保留；必要語句修正，順序不變 |
| 29 | markdown | 258 | ### 索引、座標與局部修改；原單元：在文獻中，可以看到不同的影像數值表示方式： ``` 0 - 255 其中 0 表示黑色，255 表示白色 0 - 1 其中 0 表；原圖／參數／操作保留；必要語句修正，順序不變 |
| 30 | code | 264 | ### 索引、座標與局部修改；原單元：## Slicing: Set the first ten lines to "black" (0) venusaur[:10] ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 31 | code | 269 | ### 索引、座標與局部修改；原單元：plt.imshow(venusaur, cmap='gray');；原圖／參數／操作保留；必要語句修正，順序不變 |
| 32 | markdown | 273 | ### 索引、座標與局部修改；原單元：遮罩（使用 `boolean` 遮罩進行索引）在需要選取一組像素以進行操作時非常有用。遮罩可以是任何與影像形狀相同（或可廣播為影像；原圖／參數／操作保留；必要語句修正，順序不變 |
| 33 | code | 279 | ### 索引、座標與局部修改；原單元：np.mgrid[0:nrows, 0:ncols]；mgrid示意合併到cell34建立nrows/ncols後 |
| 34 | code | 283 | ### 索引、座標與局部修改；原單元：nrows, ncols = venusaur.shape row, col = np.mgrid[0:nrows, 0:ncol；原圖／參數／操作保留；必要語句修正，順序不變 |
| 35 | markdown | 292 | ## 影像資料型別與數值範圍；原單元：## [影像資料型別及其意義](https://scikit-image.org/docs/stable/user_guide/d；原圖／參數／操作保留；必要語句修正，順序不變 |
| 36 | markdown | 300 | ## 影像資料型別與數值範圍；原單元：在 `skimage` 中，影像就是一般的 `NumPy` 陣列，而 `NumPy` 支援多種資料型別（即 "dtypes"）。為；原圖／參數／操作保留；必要語句修正，順序不變 |
| 37 | markdown | 316 | ## 影像資料型別與數值範圍；原單元：**請注意，雖然浮點數影像的資料型別本身可以超出範圍，但其像素值應限制在 0 到 1 之間；而所有整數型別的像素強度則可以涵蓋整個；原圖／參數／操作保留；必要語句修正，順序不變 |
| 38 | code | 324 | ## 影像資料型別與數值範圍；原單元：image = np.arange(0, 50, 10, dtype=np.uint8) print(image.astype(f；原圖／參數／操作保留；必要語句修正，順序不變 |
| 39 | markdown | 330 | ### 輸入型別與轉換工具；原單元：### 輸入型別 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 40 | markdown | 336 | ### 輸入型別與轉換工具；原單元：雖然我們的目標是保留輸入影像的數值範圍與資料型別，但某些函式僅支援特定的資料型別。在這種情況下，輸入影像將會（如果可行的話）自動轉；原圖／參數／操作保留；必要語句修正，順序不變 |
| 41 | code | 349 | ### 輸入型別與轉換工具；原單元：# These functions convert images to the desired dtype and properl；原圖／參數／操作保留；必要語句修正，順序不變 |
| 42 | markdown | 356 | ### 輸入型別與轉換工具；原單元：📌 此外，有些函式提供 `preserve_range` 參數，用於在範圍轉換是方便但非必要的情況下使用。然而，在某些情況下，**；原圖／參數／操作保留；必要語句修正，順序不變 |
| 43 | code | 364 | ### 輸入型別與轉換工具；原單元：coin = ski.data.coins() coin.dtype, coin.min(), coin.max(), coin.；原圖／參數／操作保留；必要語句修正，順序不變 |
| 44 | code | 369 | ### 輸入型別與轉換工具；原單元：rescaled = ski.transform.rescale(coin, 0.5) rescaled.dtype, np.ro；原圖／參數／操作保留；必要語句修正，順序不變 |
| 45 | code | 374 | ### 輸入型別與轉換工具；原單元：rescaled = ski.transform.rescale(coin, 0.5, preserve_range=True) ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 46 | markdown | 379 | ### 輸入型別與轉換工具；原單元：我們建議使用浮點數格式，因為 `scikit-image` 在內部主要也是採用此格式。 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 47 | markdown | 385 | ### 輸出型別；原單元：### 輸出型別 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 48 | markdown | 391 | ### 輸出型別；原單元：函式的輸出型別是由函式的作者所決定，並會在文件中說明以方便使用者查閱。若使用者需要特定的輸出型別（例如用於顯示用途），可以這樣撰寫；原圖／參數／操作保留；必要語句修正，順序不變 |
| 49 | code | 397 | ### 輸出型別；原單元：out = ski.util.img_as_uint(ski.filters.sobel(coin)) plt.imshow(ou；原圖／參數／操作保留；必要語句修正，順序不變 |
| 50 | markdown | 403 | ## 影像輸入與輸出；原單元：## 🔍 [影像輸入／輸出（Image I/O）](https://scikit-image.org/docs/stable/ap；原圖／參數／操作保留；必要語句修正，順序不變 |
| 51 | markdown | 411 | ## 影像輸入與輸出；原單元：通常，我們不會使用 scikit-image 範例資料集中的輸入影像。這些影像通常是以 JPEG 或 PNG 格式儲存。由於 sc；原圖／參數／操作保留；必要語句修正，順序不變 |
| 52 | code | 417 | ## 影像輸入與輸出；原單元：ski.io.find_available_plugins()；舊plugin清單併入cell57收合說明 |
| 53 | code | 421 | ## 影像輸入與輸出；原單元：ski.io.plugin_info('tifffile')；舊plugin_info併入cell57收合說明 |
| 54 | code | 425 | ## 影像輸入與輸出；原單元：ski.io.find_available_plugins(loaded=True)；舊loaded plugin清單併入cell57收合說明 |
| 55 | code | 429 | ## 影像輸入與輸出；原單元：ski.io.use_plugin('imageio', 'imread')；淘汰plugin設定移除，使用預設讀圖 |
| 56 | code | 433 | ## 影像輸入與輸出；原單元：coins = ski.data.coins() coins；原圖／參數／操作保留；必要語句修正，順序不變 |
| 57 | markdown | 438 | ## 影像輸入與輸出；原單元：```python io.use_plugin('pil', 'imread') # Use only the imread ca；原圖／參數／操作保留；必要語句修正，順序不變 |
| 58 | markdown | 446 | ## 色彩、通道與透明度；原單元：## [色彩（Colors）](https://scikit-image.org/docs/stable/api/skimage.；原圖／參數／操作保留；必要語句修正，順序不變 |
| 59 | markdown | 454 | ## 色彩、通道與透明度；原單元：<center><img src="https://drive.google.com/uc?id=171RCIjgsQYhQTRv；原圖／參數／操作保留；必要語句修正，順序不變 |
| 60 | code | 465 | ## 色彩、通道與透明度；原單元：download_from_pokemondb('https://img.pokemondb.net/sprites/diamon；原圖／參數／操作保留；必要語句修正，順序不變 |
| 61 | code | 470 | ## 色彩、通道與透明度；原單元：ski.io.use_plugin('imread', 'imread')；淘汰imread plugin移除，cell60/62保留RGBA讀取 |
| 62 | code | 474 | ## 色彩、通道與透明度；原單元：venusaurf = ski.io.imread("venusaur-f.png") print("Shape:", venus；原圖／參數／操作保留；必要語句修正，順序不變 |
| 63 | code | 480 | ## 色彩、通道與透明度；原單元：ski.io.imsave('venusaur-f.jpg', venusaurf[:,:,:3])；原RGB JPEG輸出改記憶體roundtrip，不覆寫資產；新增PNG對照 |
| 64 | markdown | 493 | ## 色彩、通道與透明度；原單元：與先前一樣，我們可以取得與設定像素值： ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 65 | code | 499 | ## 色彩、通道與透明度；原單元：venusaurf[10, 20], venusaurf[40, 75]；原圖／參數／操作保留；必要語句修正，順序不變 |
| 66 | markdown | 503 | ## 色彩、通道與透明度；原單元：在 alpha 通道中，0 表示完全透明，255 表示完全不透明。 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 67 | code | 509 | ## 色彩、通道與透明度；原單元：# Set the pixel at (40th row, 75th column) to "black" venusaurf[4；原圖／參數／操作保留；必要語句修正，順序不變 |
| 68 | code | 514 | ## 色彩、通道與透明度；原單元：# set the pixel at (40th row, 76st column) to "green" venusaurf[4；原圖／參數／操作保留；必要語句修正，順序不變 |
| 69 | code | 519 | ## 色彩、通道與透明度；原單元：plt.imshow(venusaurf);；原圖／參數／操作保留；必要語句修正，順序不變 |
| 70 | code | 523 | ## 色彩、通道與透明度；原單元：# https://stackoverflow.com/questions/52632718/display-image-only；原圖／參數／操作保留；必要語句修正，順序不變 |
| 71 | code | 548 | ## 色彩、通道與透明度；原單元：without_red = rgb_decomposition(venusaurf2, 'R') without_green = ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 72 | code | 554 | ## 色彩、通道與透明度；原單元：#https://stackoverflow.com/questions/49643907/clipping-input-data；三個去通道陣列皆顯示，保留原without_red/255.0 |
| 73 | markdown | 560 | ## 色彩、通道與透明度；原單元：`skimage.data` 將通道資訊儲存在最後一個軸上。`scikit-image` 函式庫中支援彩色影像的函式都有一個 `c；原圖／參數／操作保留；必要語句修正，順序不變 |
| 74 | markdown | 566 | ## 色彩、通道與透明度；原單元：在多通道影像的情況下，任意一個維度（陣列軸）都可以用來表示色彩通道，並以 `channel` 或 `ch` 表示。需要多通道資料的；原圖／參數／操作保留；必要語句修正，順序不變 |
| 75 | markdown | 572 | ## 色彩、通道與透明度；原單元：最後，對於體積式（3D）影像，例如影片、磁振造影（MRI）掃描、共軛焦顯微影像等，我們將最前面的維度稱為平面（plane），簡寫為；原圖／參數／操作保留；必要語句修正，順序不變 |
| 76 | markdown | 585 | # im3d 已是 (plane, row, column) 陣列；原單元：在 `scikit-image` 中，許多函式可以直接處理三維影像： ```python rng = np.random.defa；原三維操作保留收合code；watershed模組更正，不配置巨量隨機資料 |
| 77 | markdown | 605 | ### RGBA 轉成 RGB：和背景混合；原單元：### 從 RGBA 轉換為 RGB —— 透過 alpha 混合移除 alpha 通道 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 78 | markdown | 611 | ### RGBA 轉成 RGB：和背景混合；原單元：透過將 RGBA 影像與背景進行 alpha 混合，可使用 `rgba2rgb()` 將其轉換為 RGB 影像。(預設backgr；原圖／參數／操作保留；必要語句修正，順序不變 |
| 79 | code | 621 | ### RGBA 轉成 RGB：和背景混合；原單元：img_rgb = ski.color.rgba2rgb(venusaurf) plt.imshow(img_rgb);；原圖／參數／操作保留；必要語句修正，順序不變 |
| 80 | code | 626 | ### RGBA 轉成 RGB：和背景混合；原單元：img_rgb.shape；原圖／參數／操作保留；必要語句修正，順序不變 |
| 81 | markdown | 630 | ### 彩色與灰階之間的轉換；原單元：### 彩色與灰階值之間的轉換 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 82 | markdown | 636 | ### 彩色與灰階之間的轉換；原單元：在 [skimage.color](https://scikit-image.org/docs/stable/api/skimag；原圖／參數／操作保留；必要語句修正，順序不變 |
| 83 | markdown | 642 | ### 彩色與灰階之間的轉換；原單元：影像的*相對亮度*（relative luminance）指的是每個點所發出的光強度。由於人眼對不同顏色的敏感程度不同，因此各色彩；保留原亮度公式，明示線性RGB前提與兩套權重/量尺差異 |
| 84 | code | 652 | ### 彩色與灰階之間的轉換；原單元：# https://en.wikipedia.org/wiki/Grayscale#Converting_color_to_gra；原圖／參數／操作保留；必要語句修正，順序不變 |
| 85 | markdown | 664 | ### 彩色與灰階之間的轉換；原單元：使用 `gray2rgb()` 將灰階影像轉換為 RGB 時，會簡單地將灰階值複製到三個色彩通道上。 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 86 | markdown | 670 | ### HSV 與 CIELAB；原單元：### 色彩模型之間的轉換 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 87 | markdown | 676 | ### HSV 與 CIELAB；原單元：彩色影像可以使用不同的色彩空間來表示。其中最常見的是 RGB 色彩空間，影像包含紅色、綠色與藍色通道。然而，也有其他廣泛使用的色彩；原圖／參數／操作保留；必要語句修正，順序不變 |
| 88 | code | 682 | ### HSV 與 CIELAB；原單元：# bright saturated red red_pixel_rgb = np.array([[[255, 0, 0]]], ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 89 | markdown | 688 | ### HSV 與 CIELAB；原單元：現在我們將使用<u>CIELAB 色彩空間</u>來分離 chelsea 圖像中的眼睛。我們首先繪製 `L`、`a`、`b` 三個；原圖／參數／操作保留；必要語句修正，順序不變 |
| 90 | code | 694 | ### HSV 與 CIELAB；原單元：color_image = ski.util.img_as_float(ski.data.chelsea()) lab_image；原圖／參數／操作保留；必要語句修正，順序不變 |
| 91 | markdown | 702 | ### HSV 與 CIELAB；原單元：> 來描述人眼可見的所有顏色的最完備的色彩模型。它是為這個特殊目的而由國際照明委員會（Commission Internation；原圖／參數／操作保留；必要語句修正，順序不變 |
| 92 | code | 708 | ### HSV 與 CIELAB；原單元：plt.imshow(color_image);；原圖／參數／操作保留；必要語句修正，順序不變 |
| 93 | markdown | 712 | ### HSV 與 CIELAB；原單元：這張圖像顯示錯誤，因為它不是 RGB 影像，但直方圖本身仍然是有效的。請注意，我們是以 RGB 的順序繪製直方圖，因此這裡的 "A；原圖／參數／操作保留；必要語句修正，順序不變 |
| 94 | code | 718 | ### HSV 與 CIELAB；原單元：imshow_with_histogram(lab_image);；保留故意LAB錯色顯示與未截斷直方圖 |
| 95 | markdown | 722 | ### HSV 與 CIELAB；原單元：為了取得眼睛和鼻子的區域，可以注意到它在 "A" 分量中相當明亮，因此我們可以選取綠色直方圖右側的尾端（> 30）。 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 96 | code | 728 | ### HSV 與 CIELAB；原單元：eyes = a < 0 nose = a > 30 mask = eyes／nose mask = np.expand_dims；原圖／參數／操作保留；必要語句修正，順序不變 |
| 97 | code | 736 | ### HSV 與 CIELAB；原單元：plt.imshow(color_image*mask);；原圖／參數／操作保留；必要語句修正，順序不變 |
| 98 | markdown | 740 | ### 影像反相；原單元：### 影像反相（Image Inversion） ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 99 | markdown | 746 | ### 影像反相；原單元：反相影像也稱為互補影像（complementary image）。對於二值影像，`True` 值會變成 `False`，反之亦然；；原圖／參數／操作保留；必要語句修正，順序不變 |
| 100 | code | 752 | ### 影像反相；原單元：inverted_img = ski.util.invert(venusaur) plt.imshow(inverted_img,；原圖／參數／操作保留；必要語句修正，順序不變 |
| 101 | markdown | 757 | ## OpenCV 的 BGR 與 scikit-image 的 RGB；原單元：# 📘 使用 OpenCV 進行影像處理 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 102 | markdown | 763 | ## OpenCV 的 BGR 與 scikit-image 的 RGB；原單元：有時你可能需要在 `skimage` 與 OpenCV 之間互相使用影像資料。OpenCV 的影像資料可以直接（不需複製）透過 `；原圖／參數／操作保留；必要語句修正，順序不變 |
| 103 | code | 769 | ## OpenCV 的 BGR 與 scikit-image 的 RGB；原單元：ski.io.use_plugin('matplotlib', 'imread')；移除全域plugin切換 |
| 104 | code | 773 | ## OpenCV 的 BGR 與 scikit-image 的 RGB；原單元：from google.colab.patches import cv2_imshow download_from_pokemon；原Blastoise用OpenCV讀取；取代僅Colab可用的cv2_imshow |
| 105 | code | 780 | ## OpenCV 的 BGR 與 scikit-image 的 RGB；原單元：blastoise.shape, blastoise.dtype；原圖／參數／操作保留；必要語句修正，順序不變 |
| 106 | markdown | 784 | ### BGR 與 RGB 的轉換；原單元：## BGR 與 RGB 之間的轉換 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 107 | markdown | 790 | ### BGR 與 RGB 的轉換；原單元：在 skimage 和 OpenCV 中，彩色影像具有三個維度：寬度、高度與色彩。不過 RGB 和 BGR 使用的是相同的色彩空間；原圖／參數／操作保留；必要語句修正，順序不變 |
| 108 | code | 796 | ### BGR 與 RGB 的轉換；原單元：plt.imshow(blastoise);；原圖／參數／操作保留；必要語句修正，順序不變 |
| 109 | code | 800 | ### BGR 與 RGB 的轉換；原單元：blastoise = blastoise[:, :, ::-1] plt.imshow(blastoise);；原圖／參數／操作保留；必要語句修正，順序不變 |
| 110 | markdown | 805 | ### 兩套函式庫之間的型別轉換；原單元：## 在 `skimage` 中使用來自 OpenCV 的影像 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 111 | markdown | 811 | ### 兩套函式庫之間的型別轉換；原單元：如果 `cv_image` 是一個無號位元組（unsigned byte）陣列，`skimage` 預設就能正確處理它。若你偏好使；原圖／參數／操作保留；必要語句修正，順序不變 |
| 112 | markdown | 817 | ## 用 Matplotlib 並排顯示影像；原單元：# 📘 [使用 matplotlib 顯示影像](https://scikit-image.org/docs/stable/use；原圖／參數／操作保留；必要語句修正，順序不變 |
| 113 | code | 825 | ## 用 Matplotlib 並排顯示影像；原單元：download_from_pokemondb('https://img.pokemondb.net/sprites/go/nor；原圖／參數／操作保留；必要語句修正，順序不變 |
| 114 | code | 830 | ## 用 Matplotlib 並排顯示影像；原單元：f, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 5)) ax0.imshow(ch；原圖／參數／操作保留；必要語句修正，順序不變 |
| 115 | markdown | 839 | ## 用 Matplotlib 並排顯示影像；原單元：[cv2.imshow](https://stackoverflow.com/questions/69140016/graysca；原圖／參數／操作保留；必要語句修正，順序不變 |
| 116 | markdown | 845 | ## 對比與曝光調整；原單元：# 📘 [對比與曝光調整（Contrast and Exposure）](https://scikit-image.org/doc；原圖／參數／操作保留；必要語句修正，順序不變 |
| 117 | markdown | 853 | ## 對比與曝光調整；原單元：影像像素的數值範圍取決於影像的 `dtype`，例如 `uint8` 影像的範圍是 `[0, 255]`，而浮點數影像的範圍則是 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 118 | markdown | 859 | ## 對比與曝光調整；原單元：第一類方法是針對影像強度套用非線性函數，且這些函數與特定影像的像素值無關。這類方法常用於校正感測器或人眼等接收器的已知非線性反應。；原圖／參數／操作保留；必要語句修正，順序不變 |
| 119 | markdown | 865 | ## 對比與曝光調整；原單元：其他方法則會根據影像的直方圖重新分佈像素值。像素值的直方圖可以使用 `skimage.exposure.histogram()` ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 120 | code | 871 | ## 對比與曝光調整；原單元：image = np.array([[1, 3], [1, 1]]) ski.exposure.histogram(image)；原圖／參數／操作保留；必要語句修正，順序不變 |
| 121 | markdown | 876 | ## 對比與曝光調整；原單元：`histogram()` 會回傳每個數值區間（bin）內的像素數量，以及該區間的**中心值**。因此，`histogram()`；原圖／參數／操作保留；必要語句修正，順序不變 |
| 122 | markdown | 884 | ### 重設強度範圍；原單元：## 重設強度值範圍（Rescaling Intensity Values） ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 123 | markdown | 890 | ### 重設強度範圍；原單元：在可行的情況下，函式應避免盲目地拉伸影像強度值（例如將浮點數影像的 `min` 和 `max` 強度強制調整為 0 和 1），因為；原圖／參數／操作保留；必要語句修正，順序不變 |
| 124 | markdown | 896 | ### 重設強度範圍；原單元：不過有時候，你**已確定的影像理應**涵蓋整個強度範圍，卻實際上沒有。例如，有些相機會以每個像素 10、12 或 14 位元的深度；原圖／參數／操作保留；必要語句修正，順序不變 |
| 125 | markdown | 906 | ### 重設強度範圍；原單元：> 此處的 `in_range` 參數設定為 10 位元影像的最大範圍。`rescale_intensity` 預設會將 `in_；原圖／參數／操作保留；必要語句修正，順序不變 |
| 126 | markdown | 912 | ### 重設強度範圍；原單元：另外一個使用情境是，許多人經常使用帶符號的資料型別（signed dtypes）來表示影像，即使實際上只操作影像中的正值（例如在 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 127 | markdown | 918 | # img_int32 是需要保留正負資訊的陣列；原單元：換句話說，當從帶符號的資料型別轉換為無符號型別時，負值會被截斷為 0。（若在帶符號型別之間轉換，則會保留負值。）為了避免這種截斷行；修正signed→unsigned浮點範圍，再img_as_ubyte |
| 128 | markdown | 932 | # img_int32 是需要保留正負資訊的陣列；原單元：即使影像已使用整個數值範圍，有時其實際像素值在該範圍的兩端所佔比例非常小。在這種情況下，透過使用影像的百分位數來截斷像素值，可以提；原圖／參數／操作保留；必要語句修正，順序不變 |
| 129 | code | 938 | # img_int32 是需要保留正負資訊的陣列；原單元：# Load an example image img = ski.data.moon() # Contrast stretchi；原圖／參數／操作保留；必要語句修正，順序不變 |
| 130 | code | 953 | # img_int32 是需要保留正負資訊的陣列；原單元：# Display results fig = plt.figure(figsize=(16, 10)) axes = np.ze；原moon圖/直方圖/CDF保留，CDF標籤修正 |
| 131 | markdown | 981 | # img_int32 是需要保留正負資訊的陣列；原單元：函式 `equalize_hist()` 會將像素值的累積分布函數（cdf）**映射到線性的 cdf** 上（[直方圖均衡化](h；原圖／參數／操作保留；必要語句修正，順序不變 |
| 132 | markdown | 989 | ## 影像的幾何變換；原單元：# 📘 [影像的幾何轉換（Geometrical Transformations of Images）](https://scik；原圖／參數／操作保留；必要語句修正，順序不變 |
| 133 | markdown | 997 | ### 裁切、重縮放與調整尺寸；原單元：## 裁切、調整大小與重縮影像（Cropping, Resizing, and Rescaling Images） ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 134 | markdown | 1003 | ### 裁切、重縮放與調整尺寸；原單元：由於影像是 `NumPy` 陣列，因此可以透過簡單的切片操作來裁切影像。以下範例中，我們裁切出一個 `100x100` 的正方形區；原圖／參數／操作保留；必要語句修正，順序不變 |
| 135 | code | 1009 | ### 裁切、重縮放與調整尺寸；原單元：img = charizardy top_left = img[:100, :100] plt.imshow(top_left);；原圖／參數／操作保留；必要語句修正，順序不變 |
| 136 | markdown | 1015 | ### 裁切、重縮放與調整尺寸；原單元：為了改變影像的形狀，`skimage.transform` 提供了多個函式來進行重縮放（rescale）、調整大小（resize）；原圖／參數／操作保留；必要語句修正，順序不變 |
| 137 | markdown | 1021 | ### 裁切、重縮放與調整尺寸；原單元：* **Rescale**（重縮放）：根據給定的縮放因子調整影像大小。縮放因子可以是單一浮點數（套用到所有軸），也可以是對每個軸分；原圖／參數／操作保留；必要語句修正，順序不變 |
| 138 | code | 1034 | ### 裁切、重縮放與調整尺寸；原單元：image = ski.color.rgb2gray(ski.color.rgba2rgb(charizardx)) image_；原圖／參數／操作保留；必要語句修正，順序不變 |
| 139 | code | 1043 | ### 裁切、重縮放與調整尺寸；原單元：fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 8)) ax = a；保留三種縮小對照，no aliasing改為具體處理名 |
| 140 | markdown | 1063 | ### 射影變換與齊次座標；原單元：## 投影轉換（同調變換，Projective Transforms / Homographies） ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 141 | markdown | 1069 | ### 射影變換與齊次座標；原單元：同調變換（homographies）是對歐幾里得空間的轉換，會保留點的**對齊關係**。某些特定類型的同調變換會額外保留更多幾何性；原圖／參數／操作保留；必要語句修正，順序不變 |
| 142 | markdown | 1077 | ### 射影變換與齊次座標；原單元：在二維歐幾里得空間中的同調變換（例如用於 2D 灰階或多通道影像），可以透過一個 `3x3 矩陣` 定義。所有類型的同調變換都可以；原圖／參數／操作保留；必要語句修正，順序不變 |
| 143 | markdown | 1083 | ### 射影變換與齊次座標；原單元：<center><img src="https://drive.google.com/uc?id=173_GZa3-rLcI6-8；原圖／參數／操作保留；必要語句修正，順序不變 |
| 144 | markdown | 1093 | ### 歐幾里得（剛體）變換；原單元：## 歐幾里得（剛性）變換（Euclidean / Rigid Transformation） ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 145 | markdown | 1099 | ### 歐幾里得（剛體）變換；原單元：歐幾里得變換，又稱為剛性變換（rigid transformation），**會保留點與點之間的歐幾里得距離**。它可以表示為**；原圖／參數／操作保留；必要語句修正，順序不變 |
| 146 | code | 1105 | ### 歐幾里得（剛體）變換；原單元：tform = ski.transform.EuclideanTransform( rotation = np.pi / 12.,；原圖／參數／操作保留；必要語句修正，順序不變 |
| 147 | code | 1113 | ### 歐幾里得（剛體）變換；原單元：matrix = np.array([[np.cos(np.pi/12), -np.sin(np.pi/12), 30], [np；原圖／參數／操作保留；必要語句修正，順序不變 |
| 148 | markdown | 1122 | ### 歐幾里得（剛體）變換；原單元：轉換物件的轉換矩陣可透過其 `tform.params` 屬性取得。若要組合多個變換，可以使用 `@` 運算子（矩陣乘法）來相乘。；說清warp反向取樣，只在tform是forward時傳inverse |
| 149 | code | 1130 | ### 歐幾里得（剛體）變換；原單元：img = ski.util.img_as_float(venusaurf) tf_img = ski.transform.war；原圖／參數／操作保留；必要語句修正，順序不變 |
| 150 | code | 1136 | ### 歐幾里得（剛體）變換；原單元：tform = ski.transform.EuclideanTransform( rotation = np.pi / 12.,；原圖／參數／操作保留；必要語句修正，順序不變 |
| 151 | code | 1147 | ### 內插階數與旋轉中心；原單元：tform = ski.transform.EuclideanTransform( rotation = np.pi / 12.,；保留原order5；補原點順時針/中心逆時針差異 |
| 152 | code | 1165 | ### 內插階數與旋轉中心；原單元：venusaurf2.dtype；原圖／參數／操作保留；必要語句修正，順序不變 |
| 153 | code | 1169 | ### 內插階數與旋轉中心；原單元：img = ski.transform.rotate(venusaurf2, 15, mode='constant', cval=；原圖／參數／操作保留；必要語句修正，順序不變 |
| 154 | code | 1174 | ### 內插階數與旋轉中心；原單元：img = ski.transform.rotate(venusaurf2, -50, mode='constant', cval；原圖／參數／操作保留；必要語句修正，順序不變 |
| 155 | markdown | 1179 | ### 相似變換；原單元：## 🔍 [相似變換（Similarity Transformation）](https://en.wikipedia.org/w；原圖／參數／操作保留；必要語句修正，順序不變 |
| 156 | markdown | 1187 | ### 相似變換；原單元：相似變換會保留物體的形狀，並結合了**縮放、平移與旋轉**等幾何操作。 ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 157 | code | 1193 | ### 相似變換；原單元：tform = ski.transform.SimilarityTransform( scale=0.5, rotation=np；原圖／參數／操作保留；必要語句修正，順序不變 |
| 158 | markdown | 1203 | ### 仿射變換；原單元：## 🔍 仿射變換（Affine Transform） ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 159 | markdown | 1209 | ### 仿射變換；原單元：仿射變換（Affine Transformation）會**保留直線**，並保留直線之間的**平行關係**。它可以分解為一個**相；仿射明示非等向縮放，修正過窄的分解說法 |
| 160 | code | 1217 | ### 仿射變換；原單元：tform = ski.transform.AffineTransform( shear=np.pi/6, ) print(tfo；原圖／參數／操作保留；必要語句修正，順序不變 |
| 161 | markdown | 1226 | ### 一般射影變換；原單元：## 🔍 射影變換（Projective Transformation / Homographies） ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 162 | markdown | 1232 | ### 一般射影變換；原單元：同調變換（homography），亦稱為射影變換（projective transformation），**會保留直線，但不一定保；原圖／參數／操作保留；必要語句修正，順序不變 |
| 163 | code | 1238 | ### 一般射影變換；原單元：matrix = np.array([[1, -0.5, 40], [0.1, 0.9, 20], [0.0015, 0.0015；原圖／參數／操作保留；必要語句修正，順序不變 |
| 164 | markdown | 1247 | ### 用對應點估計變換；原單元：## 🔍 參數估計（Parameter Estimation） ；原圖／參數／操作保留；必要語句修正，順序不變 |
| 165 | markdown | 1253 | ### 用對應點估計變換；原單元：除了上述的基本功能外，你還可以使用 **最小平方法（least-squares method**來估計幾何變換的參數。這在影像配準；原圖／參數／操作保留；必要語句修正，順序不變 |
| 166 | code | 1261 | ### 直接調整幾何參數；原單元：text = ski.data.text() src = np.array([[0, 0], [0, 50], [300, 50]；原控制點與warp方向保留；estimate更換目前from_estimate；新增geometry互動 |
| 167 | markdown | 1293 | ### 直接調整幾何參數；原單元：> 上述參數估計依賴於**點位置的精確資訊以及對應點的正確選取**。若點位置存在不確定性，可以提供權重，使得最終的變換優先符合權重；原圖／參數／操作保留；必要語句修正，順序不變 |
| 168 | markdown | 1299 | ### 直接調整幾何參數；原單元：變換估計(Transformation Estimation) 的應用範例如下： https://github.com/deepa；原圖／參數／操作保留；必要語句修正，順序不變 |
| 169 | markdown | 1307 | ## 延伸閱讀；原單元：# 📘 延伸閱讀與參考；原圖／參數／操作保留；必要語句修正，順序不變 |
| 170 | markdown | 1313 | ## 延伸閱讀；原單元：- OpenCV: 請參考 https://docs.opencv.org/4.x/d2/d96/tutorial_py_tabl；原圖／參數／操作保留；必要語句修正，順序不變 |

## 數學與行為修正

- 索引標籤改為零起算對應的實際列欄；`random_image.shape`、`mgrid` 移到相關變數建立後。
- dtype儲存範圍和影像慣例分開；astype不一概禁用。signed→unsigned先映射成0–1再轉型。
- RGBA是三維陣列四通道；alpha丟棄和背景混合分開。HSV座標不同不意味統計獨立。
- 原灰階權重與套件權重、0–255/0–1量尺差异明示，兩個原計算保留。LAB直接imshow的錯色示範仍保留並註明，Chelsea閾值 a<0/a>30保留。
- moon百分位2/98與CLAHE .03保留，CDF縱軸修成累積像素比例，不保證離散直方圖完全均勻。
- 縮小對照保留Charizard X與三種原方法，標題移除保證no aliasing。
- Venusaur幾何維持cell67/68修改後RGBA；Euclidean15°/(30,-20)、手寫同矩陣、origin純15°/order5、rotate中心15°/白255與-50°/白1、similarity.5/15°/(30,-20)、shearπ/6、原projective矩陣與text四點50×300均保留。
- warp需要output→input，前三類示範傳inverse，text原src→dst已output→input故直接傳tform；原點正角順時針和rotate中心正角逆時針分開解釋。
- 三維watershed舊模組更正；大型示意保留收合讀法，不配置100×1000×1000陣列。
- 一般仿射可含非等向縮放，不宣稱只需相似變換加單一剪切。

## 驗證限制

已做Python AST語法檢查；完整執行等root建立shared support及原圖後統一進行。未同步ipynb、未build。原圖由root提供image_path，包括所有Google Drive示意圖。
## 原圖來源

- cell 18: `venusaur.png` ← https://img.pokemondb.net/sprites/red-blue/normal/venusaur.png
- cell 60: `venusaur-f.png` ← https://img.pokemondb.net/sprites/diamond-pearl/normal/venusaur-f.png
- cell 104: `blastoise.png` ← https://img.pokemondb.net/sprites/bank/normal/blastoise-mega.png
- cell 113: `charizard-mega-x.png` ← https://img.pokemondb.net/sprites/go/normal/charizard-mega-x.png
- cell 113: `charizard-mega-y.png` ← https://img.pokemondb.net/sprites/go/normal/charizard-mega-y.png

原cell標記放在`remove-cell` code註解，避免將工程定位渲染成正文。JPEG現用imageio直接編碼，避開skimage舊plugin參數警告；PNG無損round-trip與JPEG有損對照均使用原Venusaur。

原始notebook SHA-256：`1e6a3d87ee5b5912c8a2d4eb05dd68876ae567dd595218343d57fa01cb3c3310`。本次只讀原始檔。
