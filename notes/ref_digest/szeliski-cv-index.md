# [2022] Computer Vision: Algorithms and Applications（第 2 版，Szeliski；final draft 2021-09，Springer）

`印刷頁 ↔ PDF 頁：PDF 頁 = 印刷頁 + 26`（適用於正文阿拉伯頁碼；前言羅馬頁碼不適用。全書 1232 PDF 頁）

## 一句話定位

本站 01–04 章（陣列、濾波、傅立葉、小波）的**主要教科書級參考**：把卷積、邊界處理、非線性濾波、DFT/FFT、金字塔與 wavelet 放在同一套符號下講清楚，且每個算子都附實作考量與圖例，可直接對應 NumPy/SciPy 練習。

## 關鍵主張與內容

### 章節索引（一行一章）

- Ch 1 Introduction（p.1–32）— 電腦視覺概觀與課程建議；本站不常用
- Ch 2 Image formation（p.33–106）— 幾何/光度成像、數位相機、取樣與 aliasing（2.3.1, p.84）→ 站 05/07
- Ch 3 Image processing（p.107–190）— **本書對本站最核心的一章**：point operator、線性/非線性濾波、Fourier、金字塔與小波、幾何變換 → 站 01/02/03/04
- Ch 4 Model fitting and optimization（p.191–234）— 內插、變分法與正則化、MRF → 站 02（進階分割）
- Ch 5 Deep Learning（p.235–342）— 監督/非監督學習、DNN、CNN；含 PCA（5.2.3, p.262）→ 選用
- Ch 6 Recognition（p.343–416）— 分類、偵測、語意分割；6.4.1 醫學影像分割（p.390）→ 站 02 延伸
- Ch 7 Feature detection and matching（p.417–500）— 關鍵點、邊緣與輪廓、線與消失點、分割 → 站 02
- Ch 8 Image alignment and stitching（p.501–554）— 配準、RANSAC、拼接、混色 → 選用
- Ch 9 Motion estimation（p.555–606）— 平移對齊、Fourier-based alignment（9.1.2, p.563）、光流 → 站 06（movie/motion correction 類比）
- Ch 10 Computational photography（p.607–680）— HDR、去噪、超解析、去模糊 → 站 07 雜訊模擬可參照
- Ch 11 SfM and SLAM（p.681–748）／Ch 12 Depth estimation（p.749–804）／Ch 14 Image-based rendering（p.861–918）— 本站不用
- Ch 13 3D reconstruction（p.805–860）— 表面/體積表示法 → 站 06 重建概念類比
- 附錄 A 線性代數與數值方法（p.919）／B 貝氏建模與推論（p.939）／C 資料集與軟體（p.953）；References p.973；Index p.1179

### 3.2–3.5 重點（本 PDF 實際標題：3.2 Linear filtering／3.3 More neighborhood operators／3.4 Fourier transforms／3.5 Pyramids and wavelets）

- 鄰域算子定義與線性濾波開場（p.119）；correlation `g=f⊗h` 與 convolution `g=f*h` 只差偏移正負號，convolution 的 h 稱 impulse response（p.120）。
- 卷積具交換律、結合律，且**兩影像卷積的 Fourier transform 等於各自 transform 的乘積**（convolution theorem，p.122）；LSI 算子滿足疊加與位移不變（p.122）；卷積可寫成稀疏矩陣乘法 `g=Hf`（p.122）。
- 邊界 padding 六種模式：zero / constant / clamp / wrap / mirror / extend，並示範各自模糊後的邊緣假影（p.123–124）。
- Separable filter：K×K 卷積由 K² 降到 2K 次運算；可用 SVD 判斷可分離性，只有 σ₀≠0 即可分離（p.124–125）。
- 常用核比較：box / bilinear（tent）/「Gaussian」binomial / Sobel / corner（p.125–127）；unsharp masking 用模糊差值銳化（p.126）。
- Band-pass 與 steerable filter：Gaussian 一階/二階導數、Laplacian、LoG；二階方向導數只需 3 個基底函數即可任意轉向（p.127–129）。
- Summed area table（integral image）：兩次遞迴掃描建表，任意矩形和只需 4 次取樣（p.129–130）；IIR/遞迴濾波與 FIR 的區別（p.130）。
- 非線性濾波：median filter 對 shot noise 遠優於 Gaussian；α-trimmed mean、weighted median 為折衷（p.132–133）。
- Bilateral filter = domain kernel（空間高斯）× range kernel（灰階相似度），能保邊去噪；迭代版等價於 anisotropic diffusion（p.133–136）；guided image filter 用引導影像的局部仿射模型（p.136–137）。
- 二值形態學：以 structuring element 卷積後閾值化，統一定義 dilation / erosion / majority / opening / closing（p.138–139）。
- Distance transform（city-block 兩次 raster 掃描；Euclidean 較複雜）、signed distance、skeleton/MAT（p.139–141）；connected components 與區域統計（面積、周長、質心、二階矩）（p.141）。
- Fourier：以正弦輸入量測 gain 與 phase 定義 transform；連續式與 DFT（p.143–144）；FFT 把 O(N²) 降為 O(N log₂N）（p.144）。
- **注意**：第 2 版**刪去**了 Fourier 性質表與常用 transform pair 表（叫讀者查第 1 版 Szeliski 2010, Tables 3.1–3.3），也刪去 Wiener filter 推導（p.144, p.147）——本站要引用這兩項時不要引第 2 版。
- 小核的頻率響應（Table 3.1, p.145）：box 對高頻壓抑不均勻會 ringing；binomial 分離高低頻尚可但仍留高頻導致 downsample 後 aliasing；Sobel 先線性放大再衰減。
- 2D DFT 公式與性質全部由 1D 直接推廣（p.146）；DCT 為分塊壓縮用變體，近似 KL/PCA 最佳解耦（p.147–148）。
- 去噪/銳化實務：今日多用非線性濾波或 DNN；評估用 PSNR，但更建議 SSIM 或感知相似度（p.148–149）。
- Interpolation：linear/tent 會有摺痕、cubic B-spline 偏軟、bicubic a=−1 最接近 sinc 但有 ringing、a=−0.5 可精確重現線性與二次函數；windowed sinc 品質最高（p.150–153）。
- Decimation：先低通再取樣；列出 linear/binomial/cubic/windowed sinc/QMF-9/JPEG2000 9/7 六組 2× 降採樣濾波器與頻率響應（p.153–154）。
- Gaussian pyramid：Burt–Adelson 五抽頭 `[1 4 6 4 1]/16` binomial 核，每層邊長減半（octave pyramid）（p.154–155）。
- **Laplacian pyramid**：把上採樣重建的低通版從原圖減去，得 band-pass 層；具 perfect reconstruction（p.156–157）。名稱其實是 DoG；LoG 才是真正的二階導數，兩者形狀相近（p.158）。
- Half-octave / quincunx pyramid，以及特徵偵測常用的 quarter-octave（p.159）。
- Wavelet 與 pyramid 之別：傳統金字塔 overcomplete，wavelet 是 tight frame（係數數量等同像素數）；更關鍵的差別是 wavelet **方向選擇性**更好（p.159–160）。2D wavelet 以可分離方式做，產生 HH/HL/LH 三個 detail band，HL/LH 強調水平與垂直邊緣（p.161）。
- Lifting scheme（second-generation wavelet）：先拆奇偶序列，再做可就地、可逆的 predict + update，等效低通為 `[−1/8, 1/4, 3/4, 1/4, −1/8]`（p.162–163）；steerable pyramid 用徑向弧分割頻域，overcomplete 但自逆且具方向選擇性（p.164–165）。
- Laplacian pyramid blending（蘋果/橘子經典圖）：兩張圖各建 Laplacian pyramid，遮罩建 Gaussian pyramid，逐層加權相加（p.165–167）。

### 7.1–7.2 重點（本 PDF 實際標題：7.1 Points and patches／7.2 Edges and contours）

- 關鍵點流程四階段：detection → description → matching →（或）tracking（p.421）。
- Aperture problem：無紋理區無法定位，單一方向邊緣只能沿法線對齊，需**兩個以上方向的梯度**才好定位（p.422–423）。
- Auto-correlation surface 泰勒展開得 `E≈Δuᵀ A Δu`，A = 梯度外積再與權重核卷積（p.424）；A 的特徵值分析給出定位不確定度橢圓（p.425）。
- 關鍵點度量：Shi–Tomasi 最小特徵值、Harris `det(A)−α·trace(A)²`（α=0.06）、Triggs `λ₀−αλ₁`、Brown 等人的調和平均 `det A / tr A`（p.425–426）。
- ANMS（adaptive non-maximal suppression）讓特徵點分布更均勻（p.426–427）；repeatability 為評估指標（p.428）。
- 尺度不變：Lindeberg 用 LoG 極值選尺度；Lowe 的 SIFT 用 sub-octave DoG 金字塔找 3D（空間+尺度）極值，再以二次擬合求次像素位置，並用 Hessian 比值 `Tr²/Det > 10` 剔除邊緣響應（p.429–430）。
- 方向估計：36-bin 梯度方向直方圖取主峰（p.431）；仿射不變用橢圓擬合或 MSER（p.431–432）。
- 描述子：MOPS（8×8 去偏去增益）、**SIFT（16×16 窗、4×4×8 = 128 維、正規化後截斷 0.2 再正規化）**、PCA-SIFT、RootSIFT、GLOH（log-polar bin）、BRIEF/ORB/BRISK/FREAK 位元串（p.435–439）。
- 配對策略與錯誤率：TP/FP/FN/TN 混淆矩陣、TPR/FPR/PPV/ACC、ROC 與 AUC（p.441–444）；**NNDR（最近鄰／次近鄰距離比）**比固定閾值穩健（p.444）。
- 高效配對：LSH、k-d tree 與 best-bin-first、FLANN/Faiss（p.445–447）；再用 RANSAC 幾何驗證剔除 outlier（p.447）。
- 7.2 邊緣：梯度 `J=∇I` 的模長是邊緣強度、方向垂直於等值線；微分會放大高頻雜訊，故先用 Gaussian 平滑——Gaussian 是**唯一可分離且圓對稱**的濾波器，所以幾乎所有邊緣偵測都用它（p.456–457）。
- `∇(G*I) = (∇G)*I`；再取一次散度得 LoG，並可拆成兩組可分離濾波（p.457）；實務常以 DoG 取代 LoG，若已算過 Laplacian pyramid 更省事（p.458）。
- 邊緣細化 = 沿梯度方向找強度極大值 = 找 LoG/DoG 的 zero crossing；次像素位置用線性內插求 x-intercept（p.458）。
- 尺度選擇：Elder & Zucker 依雜訊水準逐像素算出「最小可靠尺度」，再用二階導數零交越挑邊（p.459–460）。
- 彩色邊緣：直接加總各通道梯度會互相抵消（純紅到綠邊），較好的做法是各通道算 oriented energy 再合併，或比較局部色彩統計（p.460）。
- 輪廓連接：edgel 串成 chain，套用 **hysteresis 雙閾值**（Canny）去弱段（p.463）；chain code 與 arc-length 參數化 `x(s)`；正規化後取 Fourier transform 即得可比對的輪廓描述子（p.463–464）；直接平滑會使曲線收縮，需加補償項（p.464–465）。
- Martin/Arbeláez 等人以機器學習合併亮度+色彩+紋理三種 cue，勝過單純梯度法（p.461）。

## 可引用的公式/圖表

- 式 (3.12)–(3.15) correlation vs convolution，Fig 3.10 卷積示意（p.120）
- Fig 3.13 六種 padding 及其模糊結果（p.123）
- 式 (3.20)–(3.21) 可分離性與 SVD 判準（p.124–125）
- **Fig 3.14 五種可分離核（box / bilinear / Gaussian / Sobel / corner）的 2D 核、1D 核與濾波結果**（p.125）— 站 02 最佳教學圖
- 式 (3.24) Gaussian、(3.25) Laplacian、(3.26) LoG（p.127）
- Fig 3.17 + 式 (3.30)–(3.32) summed area table（p.129–130）
- Fig 3.18 median / bilateral 對 Gaussian noise 與 shot noise 的比較（p.131）
- 式 (3.34)–(3.37) bilateral filter 完整權重（p.134）
- Fig 3.22 形態學六圖（原圖/膨脹/侵蝕/majority/opening/closing）（p.138）
- 式 (3.44)–(3.48) 閾值、形態學、距離變換定義（p.138–140）
- 式 (3.56)–(3.57) 連續/離散 Fourier transform；(3.59)–(3.60) 2D 版（p.143–146）
- **Table 3.1 六個小核的解析頻率響應與圖示**（p.145）— 站 03 講「濾波器的頻域長相」直接用
- 式 (3.62)–(3.63) 1D/2D DCT，Fig 3.25 DCT 基底（p.147–148）
- Fig 3.27 四種內插結果（bilinear / bicubic a=−1 / a=−0.5 / windowed sinc）+ Fig 3.28 其頻譜（p.151–152）
- Table 3.2 + Fig 3.30 六組 2× 降採樣濾波器係數與頻率響應（p.154–155）
- **Fig 3.31 影像金字塔、Fig 3.32 Gaussian pyramid 訊號流、Fig 3.33 Laplacian pyramid 完整分析/合成圖**（p.156–157）— 站 04 主圖
- 式 (3.70)–(3.72) DoG 與 LoG 定義，Fig 3.34 兩者在空間與頻域的比較（p.158）
- Fig 3.36 影像 wavelet 分解（PyWavelet 產生，含 HH/HL/LH）+ Fig 3.37 2D 可分離 wavelet 流程圖（p.160–161）
- Fig 3.38/3.39 lifting scheme 分析與合成訊號流（p.162–163）
- **Fig 3.41 蘋果+橘子 Laplacian pyramid blending 與 Fig 3.42 逐層細節**（p.166–167）
- 式 (7.1)–(7.9) WSSD、auto-correlation 矩陣 A、Harris 響應（p.423–425）；Fig 7.5 三種 auto-correlation surface（p.423）
- Fig 7.11 SIFT 的 sub-octave DoG 尺度空間與 3D 極值偵測（p.430）
- Fig 7.16 SIFT 描述子構造示意（p.436）
- 式 (7.14)–(7.18) TPR/FPR/PPV/ACC 與 NNDR；Fig 7.22 ROC 曲線（p.442–444）
- 式 (7.19)–(7.25) 梯度、∇Gσ、LoG 及其可分離拆解、零交越次像素定位（p.456–458）
- Fig 7.34 BG / CG / TG 與合併邊界偵測器輸出對照（p.462）
- Fig 7.35 chain code 與 arc-length 參數化（p.463）

## 適用章節

- **01 陣列基礎** — 弱關聯：p.120 的 `g=Hf` 稀疏矩陣觀點、p.122 raster-order 向量化，可當「影像即陣列」的理論註腳。
- **02 濾波與分割** — **核心**：3.2（線性濾波、可分離性、Sobel/LoG）、3.3（median/bilateral、形態學、距離變換、connected components）、7.2（邊緣偵測、Canny hysteresis、輪廓連接）。7.5 分割（p.483–490）可補充。
- **03 傅立葉** — **核心**：3.4 全節（DFT/FFT/2D DFT/DCT）+ p.122 convolution theorem + Table 3.1（p.145）。**但 Fourier 性質表與 Wiener filter 已從第 2 版刪除，需改引第 1 版。**
- **04 小波** — **核心**：3.5 全節（interpolation/decimation/Gaussian & Laplacian pyramid/wavelet/lifting/steerable pyramid/blending）。
- **05 cryo-EM 背景** — 僅 2.3.1 取樣與 aliasing（p.84）可作為「解析度上限」的類比。
- **06 cryo-EM workflow** — 弱關聯：9.1.1 階層式運動估計、9.1.2 Fourier-based alignment（p.562–563）可與 motion correction 對照。
- **07 合成資料生成** — 弱關聯：p.148 PSNR/SSIM 評估、p.122 的 spatially varying kernel 模糊模型（可類比 CTF 空間變化）；p.146 頻譜隨機噪聲場生成法（先在頻域填高斯、再逆 FFT）可直接用於雜訊模擬。
