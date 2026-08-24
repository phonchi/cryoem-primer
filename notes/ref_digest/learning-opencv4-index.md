# [2020] Learning OpenCV 4 Computer Vision with Python 3 (3rd Edition)

作者：Joseph Howse、Joe Minichino｜Packt，2020 年 2 月｜正文 p.7–342，Index 起於 p.345

印刷頁 ↔ PDF 頁：**PDF 頁 = 印刷頁 + 12**（已驗證：PDF 14 = p.2，PDF 20 = p.8）

## 一句話定位

以 OpenCV 4 + Python 3 為主軸的實作型入門書，把「影像即 NumPy 陣列」的觀念一路帶到濾波、邊緣、輪廓、分割與機器學習；本站 01/02 章可直接對照，03 章只有概念性鋪陳，04 章幾乎無對應。

## 關鍵主張與內容

- 影像在 Python 端就是 `numpy.array`，逐像素存取與切片操作是所有後續處理的基礎（p.30–34）。
- 讀寫檔案、影像 ↔ raw bytes、視訊與攝影機串流構成完整 I/O 骨架（p.26–40）。
- **Fourier transform 在本書是「動機」而非「工具」**：p.55 只提到 NumPy 的 `numpy.fft` 有 `fft2` 可算 DFT、p.56 說明 magnitude spectrum，之後 HPF/LPF 全部改用**空間域 kernel 卷積**實作。**本書沒有 `cv2.dft` 或頻域遮罩濾波的完整範例**——03 章可引它建立「濾波 = 頻率選擇」的直覺，不要引為頻域實作來源。
- 高通濾波以「原圖減去模糊圖」（`img - cv2.GaussianBlur(...)`）與 3×3 / 5×5 差分 kernel 呈現；低通即 Gaussian blur（p.56–59）。
- 邊緣 → 輪廓 → 幾何形狀是一條完整的實作鏈：Canny（p.66）→ threshold + findContours（p.68）→ 外接框/最小面積矩形/最小外接圓（p.69–72）→ 凸包與 Douglas-Peucker（p.73）→ Hough 直線與圓（p.75–78）。
- 分割獨立成 Ch4：GrabCut 前景抽取（p.98）與 Watershed 分水嶺（p.102–105），morphology 前處理夾在 Watershed 段落內。
- **本書完全沒有 wavelet / multiresolution 內容**；唯一近親是 image pyramid（見下）。

## 一頁式章節索引

- **Ch 1 Setting Up OpenCV（p.7）** — 各平台安裝 OpenCV 4 / Python 3 / NumPy / SciPy → 不適用（環境說明可略提）。
- **Ch 2 Handling Files, Cameras, and GUIs（p.25）** — 影像與視訊 I/O、把影像當陣列操作、視窗顯示 → **對應本站 01**。
  - ★ Reading/writing an image file（p.26–29）：`cv2.imread` / `cv2.imwrite` 與各種 `IMREAD_*` 旗標
  - ★ Converting between an image and raw bytes（p.29–32）：bytearray ↔ `numpy.array` reshape，理解 dtype 與 shape 的關鍵段
  - ★ Accessing image data with numpy.array（p.32–34）：切片、ROI、通道存取；**01 章「影像即陣列」的最佳引用點**
  - Reading/writing a video file（p.34）、Capturing camera frames（p.36）、Displaying an image in a window（p.38）
- **Ch 3 Processing Images with OpenCV（p.53）** — 色彩模型、濾波、邊緣、輪廓、形狀偵測 → **對應本站 02，部分 03**。
  - ★ Converting images between different color models（p.54–55）：BGR/HSV/灰階與 `cv2.cvtColor`
  - ★ Exploring the Fourier transform（p.55）+ HPFs and LPFs（p.56–60）：**03 章的直接對照**，但僅概念 + 空間域實作（見上方警語）
  - ★ Edge detection（p.60）：`cv2.Laplacian`、`cv2.medianBlur`、`split`/`merge`（p.61）
  - ★ Custom kernels – getting convoluted（p.62–64）：`cv2.filter2D` 與自訂 kernel，**02 章卷積小節的核心引用**
  - ★ Edge detection with Canny（p.66）：`cv2.Canny` 完整最小範例
  - ★ Contour detection（p.68）：`cv2.threshold` + `cv2.findContours` + `cv2.drawContours`
  - ★ Bounding box / min-area rect / min enclosing circle（p.69–72）
  - ★ Convex contours and Douglas-Peucker（p.73–75）：`cv2.convexHull`、`cv2.approxPolyDP`
  - ★ Detecting lines, circles, and other shapes（p.75–78）：`cv2.HoughLinesP`、`cv2.HoughCircles`
- **Ch 4 Depth Estimation and Segmentation（p.80）** — 深度相機、視差圖、前景偵測與分割 → **對應本站 02（分割段）**。
  - ★ Foreground detection with GrabCut（p.98–102）：`cv2.grabCut` 與遮罩後處理
  - ★ Image segmentation with Watershed（p.102–105）：`cv2.threshold`(OTSU) → `cv2.morphologyEx`/`cv2.dilate`（p.103）→ `cv2.distanceTransform` → `cv2.connectedComponents` → `cv2.watershed`（p.104）。**注意：morphology 範例在此，不在 Ch3**
  - Depth estimation with a normal camera（p.91–98，StereoSGBM）— 不適用
- **Ch 5 Detecting and Recognizing Faces（p.106）** — Haar cascade 人臉偵測與 Eigen/Fisher/LBPH 辨識 → 不適用。
- **Ch 6 Retrieving Images and Searching Using Image Descriptors（p.133）** — Harris/SIFT/SURF/ORB 特徵與 FLANN 比對 → 大致不適用；唯 **p.141「Anatomy of a keypoint」內含 image pyramid 概念說明**，可作為本站 04 的 multiresolution 側面引用。
- **Ch 7 Building Custom Object Detectors（p.168）** — HOG、SVM、BoW、NMS → 不適用；p.181–190 有 image pyramid + sliding window 的實作，04 章若談尺度空間可備用。
- **Ch 8 Tracking Objects（p.201）** — 背景相減、MeanShift/CamShift、Kalman filter → 不適用。
- **Ch 9 Camera Models and Augmented Reality（p.244）** — 相機與鏡頭參數、solvePnPRansac、3D 追蹤 → 不適用。
- **Ch 10 Introduction to Neural Networks with OpenCV（p.288）** — ANN、MNIST、第三方 DNN → 不適用。
- **Appendix A Bending Color Space with the Curves Filter（p.333）** — SciPy 插值做色彩曲線 → 不適用。

## 可引用的程式範例 / API（此書偏實作，列 API 與頁碼）

| 主題 | API / 範例 | 印刷頁 |
| --- | --- | --- |
| 建立空白影像 | `numpy.zeros`、`numpy.uint8`、`cv2.cvtColor(COLOR_GRAY2BGR)` | p.26–27 |
| 讀寫影像 | `cv2.imread` / `cv2.imwrite`、`IMREAD_GRAYSCALE`、`IMREAD_UNCHANGED`、`IMREAD_REDUCED_*` | p.28–29 |
| 影像 ↔ bytes | `bytearray` ↔ `numpy.array`、`numpy.random.randint` 造隨機影像 | p.30–31 |
| 陣列存取 | `numpy.array` 索引、ROI 切片、`item`/`itemset` | p.32–34 |
| 視訊 I/O | `cv2.VideoCapture`、`cv2.VideoWriter`、`VideoWriter_fourcc`、`CAP_PROP_*` | p.35–37 |
| 顯示與事件 | `cv2.imshow`、`cv2.waitKey`、`cv2.namedWindow`、`cv2.setMouseCallback` | p.38–40 |
| **DFT 概念** | 文字提及 `numpy.fft` 的 `fft2`；magnitude spectrum 說明 | **p.55–56** |
| **高/低通濾波** | `cv2.GaussianBlur` + `img - blurred` 得 g_hpf；`scipy.ndimage.convolve` 搭 3×3 / 5×5 差分 kernel | **p.57** |
| 自訂卷積 | `cv2.filter2D`、`cv2.sepFilter2D`（文中提及）、`numpy.array` kernel | p.57–58、p.62–63 |
| 邊緣 | `cv2.Laplacian`、`cv2.medianBlur`、`cv2.split` / `cv2.merge` | p.61 |
| Canny | `cv2.Canny(img, 200, 300)` | p.66 |
| 輪廓 | `cv2.threshold`、`cv2.findContours`（`RETR_TREE` / `RETR_EXTERNAL`、`CHAIN_APPROX_SIMPLE`）、`cv2.drawContours` | p.68–70 |
| 輪廓幾何 | `cv2.boundingRect`、`cv2.minAreaRect` + `cv2.boxPoints`、`cv2.minEnclosingCircle` | p.70–72 |
| 輪廓近似 | `cv2.arcLength`、`cv2.approxPolyDP`、`cv2.convexHull` | p.73–75 |
| Hough | `cv2.HoughLinesP`、`cv2.HoughCircles(HOUGH_GRADIENT)` | p.76–77 |
| 降採樣 | `cv2.pyrDown` | p.70、p.74 |
| GrabCut | `cv2.grabCut`、`GC_INIT_WITH_RECT`、`np.where` 造遮罩 | p.100–101 |
| Morphology | `cv2.morphologyEx(MORPH_OPEN)`、`cv2.dilate`、`cv2.threshold(THRESH_OTSU)` | p.103 |
| Watershed | `cv2.distanceTransform(DIST_L2)`、`cv2.connectedComponents`、`cv2.watershed` | p.104 |

圖表方面：p.56 有 magnitude spectrum 的說明性圖示，p.59 有 Gaussian blur 前後對照，p.67、p.72、p.74 有輪廓／外接框／凸包的結果圖，p.102、p.105 有 GrabCut 與 Watershed 分割結果圖——皆為結果截圖，非可重繪的示意圖。

## 適用章節

- **01 陣列基礎** — 主力來源。Ch2 p.26–34（imread/imwrite、bytes ↔ array、numpy.array 存取）。
- **02 濾波與分割** — 主力來源。Ch3 p.57–78（卷積、Laplacian、Canny、threshold、contour、Hough）＋ Ch4 p.98–105（GrabCut、morphology、Watershed）。
- **03 傅立葉** — 僅可引 Ch3 p.55–60 建立動機與 HPF/LPF 直覺；**頻域實作需另尋來源**。
- **04 小波** — 幾乎無對應。僅 image pyramid：p.141（概念）、p.181–190（sliding window 實作）、`cv2.pyrDown` p.70/74。
- **05 / 06 / 07 cryo-EM 相關** — 不適用；本書無電子顯微鏡、無合成資料生成內容。
