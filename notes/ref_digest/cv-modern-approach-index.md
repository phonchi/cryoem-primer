# [2012] Computer Vision: A Modern Approach（第 2 版，Forsyth & Ponce，Pearson/Prentice Hall）

印刷頁 ↔ PDF 頁：PDF 頁 = 印刷頁 + 32
（已於印刷 p.101 / p.107 / p.253 / p.661 四點驗證，全書 793 PDF 頁）

## 一句話定位

以幾何與統計為骨幹的電腦視覺經典教科書；對本站最有用的是 Part II「Early Vision: Just One Image」的線性濾波、影像梯度與紋理三章，把 convolution、Fourier、sampling、edge、texture filter bank 講得完整且有數學細節。

## 關鍵主張與內容

- 線性濾波的核心是 shift-invariant linear system，其輸出即 convolution；離散與連續 convolution、邊界效應在 4.1–4.2（p.107–118）。
- 頻域是理解濾波的自然語言：convolution theorem（訊號域摺積 = 頻域相乘）明確出現在 p.120；Fourier transform 一節在 4.3（p.118–121）。
- Sampling 與 aliasing 必須一起談：sampling theorem / Nyquist 條件在 p.126，並說明降採樣前務必先平滑（4.4.3 Smoothing and Resampling, p.126）。
- Gaussian 是預設的低通核，但作者指出它並非理想低通濾波器；smoothing 範例與 σ 的意義在 p.109，2D Gaussian kernel 圖在同節。
- 影像梯度用 derivative of Gaussian filter 一次完成平滑與微分（5.1.1, p.142）；edge 由 gradient magnitude 加 nonmaximum suppression 與 hysteresis 門檻取得（5.2.1, p.145）。
- Texture 用 filter bank（spots and bars）的輸出向量描述局部結構（6.1, p.166），再以 vector quantization / k-means 聚成 textons（6.2, p.171）。
- 去雜訊自成一節：non-local means、BM3D、learned sparse coding（6.4, p.182–186），對 cryo-EM 低訊噪比討論很有參考價值。
- Segmentation 被統一視為 clustering 問題：pixel clustering、graph-based（normalized cut）、watershed 等（Ch 9, p.255–287）。
- 模型擬合與穩健估計：Hough transform（p.290）、line/curve fitting（p.293–297）、robustness 與 RANSAC（p.299）、EM 機率模型（p.306）。
- 本書**沒有 wavelet 章節**；multiresolution 只以 image pyramid / Gaussian pyramid 形式出現（4.7, p.134–136）。

### 一頁式章節索引

**Part I — Image Formation（p.1）**

- Ch 1 Geometric Camera Models（p.3）— pinhole/weak perspective、intrinsic & extrinsic、projection matrix、camera calibration → 07（投影幾何背景）
- Ch 2 Light and Shading（p.32）— 反射模型、photometric stereo、shape from shading → —
- Ch 3 Color（p.68）— 色彩感知、linear/non-linear color space、color constancy → —

**Part II — Early Vision: Just One Image（p.105）**

- ★ Ch 4 Linear Filters（p.107）— 全書對本站最關鍵一章：convolution（4.1, p.107）、shift invariant systems（4.2, p.112）、Fourier transform（4.3, p.118）、sampling & aliasing（4.4, p.121）、filters as templates 與 normalized correlation（4.5–4.6, p.131）、image pyramids（4.7, p.134）→ 02 + 03 + 04
- ★ Ch 5 Local Image Features（p.141）— image gradient（5.1, p.141）、gradient-based edge detector（5.2, p.144）、corner detection（5.3, p.148）、SIFT/HOG（5.4, p.155）→ 02
- ★ Ch 6 Texture（p.164）— filter bank 紋理表示（6.1, p.166）、textons 與 k-means（6.2, p.171）、texture synthesis / inpainting（6.3, p.176）、image denoising（6.4, p.182）、shape from texture（6.5, p.187）→ 02（去雜訊段落亦可支援 07）

**Part III — Early Vision: Multiple Images（p.195）**

- Ch 7 Stereopsis（p.197）— epipolar geometry、binocular fusion → —
- Ch 8 Structure from Motion（p.221）— 由多視角回推 3D 結構與相機運動 → 06（多視角重建的類比，僅概念層）

**Part IV — Mid-Level Vision（p.253）**

- ★ Ch 9 Segmentation by Clustering（p.255）— Gestalt 分群原則（9.1, p.256）、pixel clustering / k-means / watershed（9.3, p.268）、graph-based segmentation（9.4, p.277）、實務建議（9.5, p.285）→ 02
- ★ Ch 10 Grouping and Model Fitting（p.290）— Hough transform（10.1, p.290）、line & curve fitting（10.2–10.3, p.293）、robustness/RANSAC（10.4, p.299）、EM 機率擬合（10.5, p.306）、model selection（10.7, p.319）→ 02（形狀擬合與穩健估計）
- Ch 11 Tracking（p.326）— Kalman filter（11.3, p.339）、particle filter（11.5, p.350）→ —

**Part V — High-Level Vision（p.365）**

- Ch 12 Registration（p.367）— ICP（12.1.1, p.368）、model-based registration（12.2, p.375）、deformable registration 與醫學影像應用（12.3.3, p.383）→ 06（對位/alignment 概念背景）
- Ch 13 Smooth Surfaces and Their Outlines（p.391）— 微分幾何與輪廓 → —
- Ch 14 Range Data（p.422）— 深度感測、range image segmentation 與 registration、Kinect → —
- Ch 15 Learning to Classify（p.457）— loss、overfitting、cross-validation、SVM/kernel/boosting → 06（particle picking 的分類器背景）
- Ch 16 Classifying Images（p.482）— visual words、spatial pyramid kernel、PCA 降維（16.1.5, p.493）→ 06（PCA 段落可對照 2D classification）
- Ch 17 Detecting Objects in Images（p.519）— sliding window（17.1, p.519）、deformable part model（17.2, p.530）→ 06（particle picking 的偵測框架）
- Ch 18 Topics in Object Recognition（p.540）— 辨識任務的開放問題討論 → —

**Part VI — Applications and Topics（p.557）**

- Ch 19 Image-Based Modeling and Rendering（p.559）— visual hull、PMVS、light field → —
- Ch 20 Looking at People（p.590）— HMM、人體姿態與動作辨識 → —
- Ch 21 Image Search and Retrieval（p.627）— 檢索與標註預測 → —

**Part VII — Background Material（p.661）**

- Ch 22 Optimization Techniques（p.663）— linear least squares（22.1, p.663）、nonlinear least squares（22.2, p.669）、sparse coding & dictionary learning（22.3, p.672）、min-cut/max-flow（22.4, p.675）→ 背景數學參考

## 可引用的公式/圖表

- Convolution 定義與離散/連續形式：4.1–4.2（p.107–118）。
- 2D symmetric Gaussian kernel 公式與示意圖（Figure 4.2）、σ 的效果：Example: Smoothing with a Gaussian（p.109）。
- Convolution theorem（空間域摺積 = 頻域相乘）：p.120。
- Fourier transform 定義與影像頻譜：4.3.1（p.119）。
- Sampling theorem / Nyquist 條件與縮圖前先平滑：4.4（p.121–126，Nyquist 敘述在 p.126）。
- Gaussian pyramid 建構與應用：4.7.1–4.7.2（p.135–136）。
- Derivative of Gaussian filter：5.1.1（p.142）。
- Edge detection 流程（gradient magnitude → nonmaximum suppression → hysteresis 雙門檻），含 Figure 5.5：5.2.1（p.145）。
- Corner detection（Harris 型）：5.3.1（p.149）。
- Texture filter bank：spots and bars（6.1.1, p.167）；textons 與 k-means（6.2.1–6.2.2, p.172）。
- Denoising：non-local means（p.183）、BM3D（p.183）、learned sparse coding（p.184）、結果比較（p.186）。
- Graph-based segmentation / normalized cut：9.4（p.277）。
- Hough transform 累積器（10.1, p.290）；RANSAC（10.4, p.299）。
- PCA 降維：16.1.5（p.493）。

## 適用章節

- **02 濾波與分割**（最主要）：Ch 4（convolution、filter、smoothing）、Ch 5（gradient、edge）、Ch 6（texture filter bank、denoising）、Ch 9（segmentation）、Ch 10（Hough、RANSAC）。
- **03 傅立葉**：Ch 4.3–4.4（Fourier transform p.118、convolution theorem p.120、sampling/Nyquist p.126）。
- **04 小波**：本書無 wavelet 內容，只能引 Ch 4.7 image pyramid / Gaussian pyramid（p.134–136）作為 multiresolution 的入口；wavelet 本身需另尋來源（如 Szeliski）。
- **01 陣列基礎**：可借 Ch 4 開頭「影像即像素權重和」的敘述（p.107），但本書不談 NumPy。
- **06 cryo-EM workflow**：Ch 12 registration（p.367）、Ch 15 分類（p.457）、Ch 17 sliding-window 偵測（p.519）可作 alignment 與 particle picking 的方法論類比。
- **07 合成資料生成**：Ch 1 projection 幾何（p.3）、Ch 6.4 雜訊模型與去雜訊（p.182）。
