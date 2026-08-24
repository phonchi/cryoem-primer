# Cryo-EM_Overview.pptx 抽取對照表（manifest）

- 來源：`/home/phonchi/cryo-em_general/Cryo-EM_Overview.pptx`（26 張投影片）
- 圖片輸出：`/home/phonchi/cryoem-primer/book/images/pptx/`，命名 `s{NN}_{k}.png`
- 統計：**共 88 張圖成功輸出**（72 張直接抽出、16 張 EMF 經 LibreOffice 轉 PNG，全部成功、無 SKIPPED）；另有 4 張過小裝飾圖略過、7 個 rId 在 pptx 內部即是壞連結（Target="NULL"，原檔就沒有圖，非抽取失敗）
- 輸出總量約 15 MB
- 各投影片整頁 render（佈局參考用）：`/home/phonchi/cryoem-primer/notes/pptx_slides/slide-01.jpg` … `slide-26.jpg`

## 轉檔與重複圖備註

- **EMF 轉檔**：PIL 無法解 EMF/WMF，改用 `soffice --convert-to pdf` + `pdftoppm -r 200` 再自動裁白邊，皆成功。受影響檔案：s04_1、s04_2、s04_3、s05_1、s06_1、s06_3、s10_4、s10_5、s12_7、s12_8、s13_3、s15_2、s15_3、s15_7、s15_8、s22_12。
- **同一張圖在多張投影片重複出現**（blob 內容相同，網站引用時擇一即可）：
  - 紅色 70S ribosome 3D map：`s04_3` = `s10_4` = `s12_7` = `s13_3`（另一視角 `s10_5` = `s12_8`）
  - CTF 同心圓環（2D CTF pattern）：`s04_2` = `s22_12`
  - 粒子堆疊（particle stack 示意）：`s11_3` = `s13_1`
  - 2D class average 直欄（3 張一欄）：`s11_5` = `s13_4`、`s11_6` = `s13_5`
- **黑色矩形圖**：`s06_4`、`s06_5`、`s06_6`、`s11_1`、`s11_2`、`s11_4`、`s13_2`、`s15_5` 是內嵌影片的黑色 poster frame 或黑色佔位框，**網站不需使用**。
- **深色簡報截圖**：`s07_1`、`s10_1`、`s15_1`、`s20_1` 是從另一份深色簡報剪來的文字截圖（bullet list），內容已在下方逐字列出，網站可直接重打文字、不必用圖。
- 壞連結（原檔即遺失）：slide 4 的 rId6/rId10/rId14、slide 13 的 rId4/rId5/rId6/rId8。從 render 看不影響版面主圖。

## 章節對應總表

| 投影片 | 網站章節 |
|---|---|
| 1–4 | 05_background |
| 5–13 | 06_workflow（主流程） |
| 14–15、25 | 06_workflow 課後補（heterogeneity） |
| 16–21、23–24、26 | 06_workflow 課後補（工具連結／CTF／picking／旋轉表示／motion correction 補充） |
| 22 | 07_synthetic_data 引言圖（同時標 06 課後補） |

---

## Slide 1 — Cryo-EM Image Processing Overview

**文字**：標題頁。「Cryo-EM Image Processing Overview / Szu-Chi Chung」。
**Notes**：無。
**圖片**：無。
**建議章節**：無（標題頁）。

## Slide 2 — The 3D structure of the Proteins

**文字**：Protein performs a vast array of functions within organisms。連結：https://deepmind.com/research/case-studies/alphafold
**Speaker notes**（重點）：蛋白質在生物體內負責催化代謝反應、DNA 複製、對刺激反應、提供結構、運輸分子。數十年來科學家嘗試僅從胺基酸序列可靠地決定蛋白質結構，即「protein folding problem」。AlphaFold 2 的成功可能提供洞見；但模型以 PDB 訓練、偏向易結晶的蛋白，且以 MSA 為輸入，在 MSA 淺或不具資訊性時（蛋白設計、突變序列、抗體）仍待驗證。AlphaFold 具內部信心指標可標示預測可靠區段。（Notes 亦含 CASP14 相關評論連結。）

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s02_1.png | 837×513 | DeepMind 風格資訊圖：1 個胺基酸 → 20 種胺基酸 → 蛋白質由數百個組成 → 人體 2 萬種 → 地球已知 2 億種（綠底插畫） |
| s02_2.png | 161×189 | 胺基酸鏈（珠鏈狀）示意圖，標註 "Amino acids" |
| s02_3.png | 241×243 | 摺疊後的蛋白質 ribbon 圖（藍綠色 α-helix 結構） |

**建議章節**：05_background（蛋白質結構背景；三張可組成「序列→鏈→摺疊結構」敘事）。

## Slide 3 — Cryo-electron microscopy (Cryo-EM) in structural biology

**文字**：Insight for biological function [1] and the strength of cryo-EM compared with x-ray crystallography：Efficiently determine the 3D structure；Recently, computational methods combined with cryo-EM have also been studied in various works [3]；A general method that can analyze conformations；Provide information for developing vaccines/drug。連結：https://cryosparc.com/blog/2019-nCoV
**Speaker notes**（重點）：結構生物學的基本信念——能夠足夠細緻地直接觀察巨分子，就能理解 3D 結構如何決定生物功能。NMR 難解大型蛋白。目前已知約 2 億種蛋白、每年新增 3 千萬，但已知確切結構者只占極小部分；預測未知蛋白結構有助於對抗疾病、開發新藥。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s03_1.png | 600×450 | 紫色 cryo-EM 密度圖（SARS-CoV-2 spike 蛋白，出自 cryoSPARC 2019-nCoV blog） |

**建議章節**：05_background（cryo-EM 動機、與 X-ray 的比較；spike 圖適合當疫苗/藥物開發應用例）。

## Slide 4 — Cryo-electron microscopy (Cryo-EM) in structural biology [4]

**文字**：連結 https://cryoem101.org/chapter-5/ 、http://biomachina.org/courses/structures/091.pdf ；標註「Contrast transfer function (CTF)」。
**Speaker notes**（重點）：單顆粒子影像的高解析內容無法直接量測；單粒子資料品質只能靠影像集合的統計（2D class averages、3D 重建、power spectrum、variance/covariance）衡量，比晶體學的繞射點間接得多。CTF 效應在多張影像資訊合併時可以處理。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s04_1.png | 1351×740 | 成像原理示意（EMF 轉檔）：電子束（上方箭頭）穿過冰層中隨機取向的藍色粒子，下方形成 2D 投影影像 |
| s04_2.png | 166×155 | CTF 2D 同心圓環 pattern（EMF 轉檔；與 s22_12 同圖） |
| s04_3.png | 342×297 | 紅色 3D 密度圖（70S ribosome；EMF 轉檔；與 s10_4 等同圖） |
| s04_4.png | 1079×543 | 左右兩欄 micrograph 對照：Raw image vs Denoised image（4096×4096 標註） |
| s04_5.png | 453×358 | ribosome 原子模型 ribbon 圖（多色 RNA/蛋白鏈） |
| （略過） | 55×61 | 小裝飾圖（<10KB 且 <100px） |

**建議章節**：05_background（成像模型與資料特性；s04_1 是「電子束穿過冰中粒子→投影」的關鍵教學圖）。

## Slide 5 — Stage-wise computational workflow [1]

**文字**：Data source：Movies → Micrograph → Particle stacks → 3D Volumes。
**Notes**：無。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s05_1.png | 561×603 | **主流程流程圖**（EMF 轉檔）：Movie Alignment → CTF Estimation → Particle Picking → 2D Classification → Initial Model → 3D Classification → 3D Refinement，含中間產物標籤（micrographs、coordinates、particles、2D averages、initial 3D volume…） |
| s05_2.png | 226×703 | Scipion（ASCEP_SPA）protocol tree 軟體截圖：Import Data / Movies / Micrographs / Particles / 2D / 3D / Utilities |

**建議章節**：06_workflow（總覽節的主圖，s05_1 為全章骨架圖）。

## Slide 6 — (1,2) Preprocessing - Motion correction [1]

**文字**：Motion Correction：The invention of DDD；The movie mode allows Motion estimation to be performed；工具：Unblur (MPI)、MotionCor2 (GPU)、Xmipp (MPI)。右欄流程標籤：Motion Correction → CTF Estimation → Particle Picking。
**Speaker notes**：The CTF (Contrast Transfer Function) distorts the image. Estimate the parameters for each micrograph.

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s06_1.png | 1399×493 | movie frames 對齊示意（EMF 轉檔）：多張歪斜堆疊的 frame → 箭頭 → 對齊後堆疊 |
| s06_2.png | 255×255 | 單張原始 frame（極高雜訊、幾乎純噪點） |
| s06_3.png | 329×336 | motion correction 後平均的 micrograph（EMF 轉檔，可見粒子） |
| s06_4.png / s06_5.png / s06_6.png | ~133×60 | 黑色矩形（影片 poster/佔位框）——**不用** |

**建議章節**：06_workflow（motion correction 節；s06_2 vs s06_3 是「校正前後」好對比）。

## Slide 7 — Cryo-EM Image Processing (2/2)

**文字**：（render 顯示為 "Features of Cryo-EM data"）High dimensional data (p)；Low Signal to Noise Ratio (SNR) (≤ 0.1)；Lots of outlier and artifact；Particles are densely connected。標註 4096（micrograph 尺寸）。
**Speaker notes**：The acquired image is further modulated by the CTF function. The data is mixture of conformation. Heterogeneity is another problem in Cryo-EM 3D determination.

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s07_1.png | 925×819 | 深色簡報截圖（bullet 文字，同上）——可重打文字，不必用圖 |
| s07_2.png | 988×973 | 完整 4096×4096 micrograph（含座標軸的顯示截圖） |
| s07_3.png | 988×967 | micrograph 放大區域，紅圈標出粒子位置 |

**建議章節**：06_workflow（資料特性/挑戰節；s07_2+s07_3 呈現「整張圖 vs 放大找粒子」）。

## Slide 8 — (1,2) Preprocessing – CTF Estimation [1]

**文字**：工具：CTFFind4 (Intel MKL)、Gctf (GPU)、ASPIRE-CTF。連結：https://github.com/wjiang/ctfsimulation
**Speaker notes**（重點）：CTF 對 micrograph 2D power spectrum 的擬合：左半是實驗 power spectrum、右上是圓周平均、右下是擬合的 CTF²(|k|)；擬合出 defocus 2.66 μm、astigmatism 49 nm。1D PSD（黑線）、理論 CTF²（藍線）與逐週期相關係數（紅線）顯示資訊存在於 |k| < 0.3 Å⁻¹。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s08_1.png | 971×581 | 教科書圖 a–f：3D 結構→無 CTF 投影→不同 defocus（1/2/4 μm）之 2D 投影 + 底部 Information vs Spatial frequency 的 CTF 曲線圖（不同 defocus 三色曲線） |
| s08_2.png | 613×541 | 實驗 power spectrum 與擬合 CTF 的四象限對照圖（Thon rings），座標軸為 Spatial frequency (Å⁻¹) |
| s08_3.png | 1179×471 | CTF 擬合診斷曲線圖：CTF fit / Quality of fit / Power spectrum 三色曲線 vs spatial frequency |

**建議章節**：06_workflow（CTF estimation 節主圖；s08_1 解釋 CTF 概念、s08_2 解釋 Thon ring 擬合）。

## Slide 9 — (3) Particle picking

**文字**：Particle of interests / Contaminates。方法分類：Reference-based（Template matching using cross correlation；Disk, Gaussian function）；Unsupervised（Distinguish particle from background；Particle will have higher variance；KLT Picker）；Discriminative（Train a classifier，近年多為神經網路；Topaz (Pytorch)、CrYolo (Tensorflow)）。
**Speaker notes**（重點）：計算旋轉後 reference 與 micrograph 各位置的 cross-correlation，之後用較快的 peak-detection（需調 threshold 與 minimum distance 參數），程式會輸出 FOM maps。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s09_1.png | 987×1028 | micrograph 上以綠色小標記標出 picking 結果的軟體截圖 |
| s09_2.png | 989×1023 | micrograph 上藍圈（template matching 候選）+ 綠色區塊（選定區）的 picking 對照截圖 |

**建議章節**：06_workflow（particle picking 節）。

## Slide 10 — Cryo-EM Image Processing（Features of Cryo-EM data）

**文字**：（render）Features of Cryo-EM data：High dimensional data (p)；A large number of samples (n) (unlabeled)；Unknown orientation and conformation；Low SNR (≤ 0.1)。標註：100×100、…、>100,000。
**Speaker notes**：同 slide 7（CTF 調變、conformation 混合、heterogeneity 問題）。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s10_1.png | 926×819 | 深色簡報截圖（bullet 文字）——可重打文字 |
| s10_2.png | 1300×812 | 3×5 粒子小圖網格（100×100 的 noisy particle images） |
| s10_3.png | 448×91 | 一列 5 張粒子小圖（particle stack 示意） |
| s10_4.png | 342×297 | 紅色 ribosome 3D map（=s04_3） |
| s10_5.png | 328×285 | 紅色 ribosome 3D map 另一視角（=s12_8） |

**建議章節**：06_workflow（資料特性節；s10_2 適合展示「十萬張未標註低 SNR 粒子」）。

## Slide 11 — (4) 2D classification

**文字**：2D MultiReference Alignment：ISAC (MPI, GPU)；Maximum likelihood methods (EM-ML)：Relion (GPU)；Hybrid：ROME (Intel MKL)。流程標籤：Particle images → Orientation parameters（50~500）→ PCA+Kmeans → 2D models。
**Speaker notes**（重點）：演算法不保證收斂到全域最佳；outlier 影響嚴重；可能發生 group collapsing；若 K 猜錯且群不可分（雜訊大時必然），結果強烈依賴初始化。K 個 template 的初始化方式：隨機挑 K 個物件當 template，或隨機均分 K 群取組內平均。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s11_1.png / s11_2.png / s11_4.png | 黑色 | 黑色矩形（佔位/影片框）——**不用** |
| s11_3.png | 155×110 | 粒子影像堆疊（cascade 疊圖示意；=s13_1） |
| s11_5.png | 80×232 | 一欄 3 張 2D class averages（白底黑框；=s13_4） |
| s11_6.png | 75×232 | 一欄 3 張 2D class averages（黑底；=s13_5） |

**建議章節**：06_workflow（2D classification 節；流程示意建議直接用 slide render 或重繪，因原圖多為小碎圖）。

## Slide 12 — Cryo-EM Image Processing（Alignment and Clustering）

**文字**：Alignment and Clustering (2D Classification)；Noise reduction for data cleaning。（render 中有三列「粒子strip → class average → ✓/✗ → 3D map」的資料清理示意，✗ 列被剔除。）
**Notes**：無。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s12_1.png | 508×99 | 一列 5 張 noisy 粒子（含冰污染/邊緣 artifact 的壞粒子列） |
| s12_2.png | 199×199 | 對應的 2D class average（模糊、帶黑色條紋 = 壞 class） |
| s12_3.png | 659×130 | 一列 5 張 noisy 粒子 |
| s12_4.png | 259×261 | 對應的 2D class average（清楚的粒子輪廓 = 好 class） |
| s12_5.png | 659×131 | 一列 5 張 noisy 粒子 |
| s12_6.png | 258×259 | 對應的 2D class average（清楚 = 好 class） |
| s12_7.png | 342×297 | 紅色 ribosome 3D map（=s04_3） |
| s12_8.png | 328×285 | 紅色 ribosome 3D map 另一視角（=s10_5） |

**建議章節**：06_workflow（2D classification/資料清理節；strip+average 配對是「用 2D class 篩粒子」的好素材）。

## Slide 13 — (5) 3D Refinement/Classification

**文字**：Multiple initial 3D models will be generated through perturbation of the consensus 3D models。流程標籤：Initial model → Small perturbation → 3D models → 2D projections → Particle images（Calculate CC in different (φ,x,y)）→ 3D maps → Filtering（迭代迴圈）。
**Speaker notes**（重點）：找最可能產生觀測資料的參數集合；離散指派改為對所有可能指派的機率加權積分；每張影像以權重被指派到所有 orientation 與 class。角度範圍 0≤φ≤360、0≤θ≤180、0≤ψ≤360。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s13_1.png | 155×110 | 粒子堆疊 cascade（=s11_3） |
| s13_2.png | 211×106 | 黑色矩形——**不用** |
| s13_3.png | 342×297 | 紅色 ribosome 3D map（=s04_3） |
| s13_4.png / s13_5.png | 80×232 / 75×232 | 2D class average 直欄（=s11_5 / s11_6） |

**建議章節**：06_workflow（3D refinement 節；迭代流程圖建議用 slide render 或重繪）。

## Slide 14 —（無標題：Heterogeneity 方法總覽）

**文字**：僅連結 https://www.sciencedirect.com/science/article/pii/S1047847722000909
**Notes**：無。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s14_1.png | 1279×868 | Heterogeneity 重建方法總覽圖（出自上述 JSB 論文）：上半是 Encoder→Hidden Variables→Rendering Pipeline 的生成式框架示意（volume 參數化：VG/NN/GM）；下半表格比較 CryoGAN、RELION、MultiCryoGAN、E2GMM、CryoPoseNet、CryoSPARC、FSTDiff、CryoFold、CryoAI、CryoDRGN2、AtomVAE、3DFlex 等方法的變數組合 |

**建議章節**：06_workflow 課後補（heterogeneity 總覽圖，也可當進階閱讀入口）。

## Slide 15 —（無標題：Analysis of Z_i / 3DVA）

**文字**：Simple mode (Consensus structure plus a fraction of eigenvolume)；Intermediate mode (Separate weighted (kernel) reconstruction along a reaction coordinate)；Generate energy landscape and find Minimum Energy Path (MEP)。（render 另含自由能公式 π(y₁,y₂)=e^(−ΔG/kBT)。）
**Speaker notes**：連結 cryoSPARC 3D Variability Analysis (3DVA) 筆記（notion）。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s15_1.png | 1825×194 | 深色截圖標題列「Analysis of Z_i」——**不用** |
| s15_2.png | 452×420 | 3DVA latent 散佈圖：Variability Component 1 vs 2（EMF 轉檔） |
| s15_3.png | 1139×407 | 3DVA simple mode 結果：三個 ribosome 密度圖沿 variability component 漸變（紅→藍，標註 60S/40S、40S head missing；EMF 轉檔） |
| s15_4.png | 758×739 | 能量地形圖（heatmap）與白色虛線 Minimum Energy Path |
| s15_5.png | 627×101 | 黑色矩形——**不用** |
| s15_6.png | 433×303 | latent 空間分群散佈圖（5 群、5 色，K-means/UMAP 風格） |
| s15_7.png | 295×247 | 灰色 eigenvolume/密度圖（EMF 轉檔） |
| s15_8.png | 289×171 | 灰色細長密度圖（另一 eigenvolume；EMF 轉檔） |

**建議章節**：06_workflow 課後補（heterogeneity/3DVA 進階）。

## Slide 16 — References

**文字**：
[1] Singer & Sigworth (2020). Computational Methods for Single-Particle Electron Cryomicroscopy. Annual Review of Biomedical Data Science 3.
[2] Demo data - https://drive.google.com/drive/folders/12LuCFg3SdRQ9RII844wiAyG-kqp4cjJZ?usp=sharing
[3] https://github.com/geoffwoollard/learn_cryoem_math
[4] https://github.com/phonchi/Computational-CryoEM
**圖片**：無。
**建議章節**：06_workflow 參考文獻／全站 references。

## Slide 17 — Thanks for your listening!

**圖片**：s17_1.png（1013×686）— 裝飾用蛋白質密度＋原子模型混合渲染圖（藍/橘漸層，左密度右原子）。可留作裝飾素材，非必要。
**建議章節**：無（結尾頁）。

## Slide 18 — Appendix

**文字**：分節頁「Appendix」。無圖。
**建議章節**：無（分節頁）。

## Slide 19 — Notes（工具連結）

**文字**：
- Metadata and data manipulation：https://github.com/asarnow/pyem 、https://github.com/sami-chaaban/starparser
- 3D data manipulation：https://github.com/zhonge/cryodrgn 、Chimera、https://github.com/MaximilianBeckers/SPOC
- Helical：https://github.com/wjiang/HI3D
- MTF, DQE：https://github.com/M4I-nanoscopy/mtf-nps-dqe
**Speaker notes**：https://github.com/Guillawme/cryoEM-scripts
**圖片**：無。
**建議章節**：06_workflow 課後補（工具連結清單）。

## Slide 20 — (1,2) Preprocessing – CTF Estimation [1]（補充）

**文字**：（render）Find the three CTF parameters：DefocusU (Δf₁)、DefocusV (Δf₂)、DefocusAngle (α_ast)；Frequency-dependent phase shift is a function of electron wavelength λ、object defocus、spherical aberration C_s and additional phase shift ΔΦ。
**Notes**：無。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s20_1.png | 1809×819 | 深色簡報截圖（同上文字）——可重打文字 |
| s20_2.png | 363×82 | 公式：CTF = −sin(χ_φ(g)) |
| s20_3.png | 592×62 | 公式：χ_φ(g) = πλ|g|²Δf − (π/2)λ³|g|⁴C_s + Δφ |
| s20_4.png | 797×69 | 公式：Δf = ½(Δf₁+Δf₂+(Δf₁−Δf₂)cos(2(α−α_ast))) |
| s20_5.png | 665×540 | astigmatism 橢圓示意圖：x-y 平面上的橢圓、向量 g、角度 α_g/α_ast、Δf₁/Δf₂ 方向 |

**建議章節**：06_workflow 課後補（CTF 數學細節；公式建議在網站以 LaTeX 重排，s20_5 可直接用）。

## Slide 21 — (3) Particle picking（補充：denoising）

**文字**：Noise2Noise：Topaz、Restore。連結：https://github.com/eugenepalovcak/restore （render 上方另有 frame 平均分半示意：奇偶 frame 各自平均產生兩張獨立雜訊影像）。
**Notes**：無。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s21_1.png | 1847×897 | 左右對照：raw micrograph vs 還原（denoised）後的 proteasome micrograph（"Restored proteasomes"） |
| s21_2.png | 173×319 | 直向 frame 堆疊示意（多張 frame） |
| s21_3.png | 175×172 | 較小的 frame 堆疊示意 |

**建議章節**：06_workflow 課後補（picking/denoising 補充；s21_1 是 Noise2Noise 效果好素材）。

## Slide 22 — Synthetic Data

**文字**：步驟：1. Uniform project 3D map to 50 orientations to obtain 50 projected images；2. Duplicate each projected images 100 times（Convoluted each images with CTF functions）；3. Add i.i.d Gaussian Noise to the 5000 images。標籤：3D structure、Projection/Slice Operator、2D particle image、Independent, Gaussian distributed noise。
**Speaker notes**：defocus 1.5–2 μm、voltage 300ev、astigmatism angle 0.2–1.4 radian、Cs=2。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s22_1.png | 400×400 | 藍色 3D 密度圖（EMD-1056，70S ribosome，合成資料的 ground-truth map） |
| s22_2.png | 122×121 | 乾淨 2D 投影（無雜訊 projection） |
| s22_3.png | 119×119 | 加噪後的粒子影像（幾乎純噪） |
| s22_4.png | 122×120 | 乾淨 2D 投影（另一 orientation） |
| s22_5.png | 120×120 | 加噪後的粒子影像 |
| s22_6.png | 118×121 | CTF 卷積後的投影（灰階、對比反轉感） |
| s22_7.png | 119×120 | 加噪後的粒子影像 |
| s22_8.png | 124×120 | 乾淨投影（黑底白粒子） |
| s22_9.png | 121×119 | 乾淨投影（黑底白粒子，另一視角） |
| s22_10.png | 119×121 | CTF 卷積後投影 |
| s22_11.png | 377×117 | 生成模型公式：X_ij = CTF_ij Σ_{l=1}^{L} P^φ_{jl} V_l + N_ij |
| s22_12.png | 166×155 | CTF 2D 同心圓環（=s04_2） |
| （略過×3） | ≤76px | 小箭頭/裝飾（<10KB 且 <100px） |

註：此頁在原始 pptx 中同一張圖被重複貼放多次（表現「複製 100 份」的視覺效果），抽取時已依 rId 去重，故檔案數（12+3 略過）少於投影片上的圖形物件數（約 31 個）。

**建議章節**：**07_synthetic_data 引言圖**（整張 slide 的 pipeline：map→投影→CTF→加噪，s22_11 公式即前向模型；亦標 06 課後補）。合成流程可參考整頁 render `notes/pptx_slides/slide-22.jpg` 的佈局重繪，或依「s22_1 → s22_2/4/8/9 → s22_6/10 → s22_3/5/7」順序排版。

## Slide 23 — (4) 2D Processing

**文字**：Refer to https://github.com/geoffwoollard/learn_cryoem_math/blob/master/nb/2d_class_fourier.ipynb
**圖片**：無。
**建議章節**：06_workflow 課後補（2D 分類數學 notebook 連結）。

## Slide 24 — Description of spatial rotation

**文字**：To describe spatial rotations, three approaches widely used：Rotation matrix；Euler angles；Unit quaternion (Metric learning)。連結：https://github.com/geoffwoollard/learn_cryoem_math/blob/master/nb/fsc.ipynb 、https://github.com/azazellochg/3DEM-conventions
**Speaker notes**（重點）：Rotation matrix——最基本、可直接旋轉 3D 向量/物體、難以估計與最佳化。Euler angles——用三個角度表示繞三個基本軸的連續旋轉、旋轉描述不均勻、定義與軸慣例繁多。Unit quaternion——可執行基本旋轉、可在旋轉空間產生均勻/準均勻取樣、較難理解。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s24_1.png | 695×695 | 球面上的 Euler 角示意圖（z-y-z：θ、φ、ψ 與向量 r，黑底球體） |
| s24_2.png | 817×597 | FSC 曲線圖（Fourier Shell Correlation - unmasked，0.143 門檻，標註 3.74 Å） |

**建議章節**：06_workflow 課後補（旋轉表示法；s24_2 亦可用於解析度/FSC 說明）。

## Slide 25 — Motion-based methods

**文字**：Conformation variability tends to conserve mass and preserve local geometry。兩路線：(1) Represent 3D density map as 3D Gaussian mixture model (five degrees of freedom)：3D coordinates, amplitude and width；Motion is represented by the location change of 3D coordinate；**e2gmm**。(2) Represent 3D flow field as the vertices of tetrahedral mesh；The deformation field is parameterized by a 3D flow vector at each vertex of the tetrahedral mesh；**3DFlex**。
**Speaker notes**（重點）：低解析度下 Gaussian 函數是建模電子密度的自然方式；若推展到原子解析度可能需納入 atomic form factors。3DFlex 的 deformation field 以線性 FEM shape function 在每個 mesh element 內插值。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s25_1.png | 446×532 | 四面體 mesh（tetrahedral mesh）3D 線框示意（3DFlex 用） |
| s25_2.png | 219×214 | ribosome 原子模型細線 wireframe 圖 |
| s25_3.png | 229×209 | 灰色表面密度圖（GMM 表示的 map） |
| s25_4.png | 241×209 | 灰色表面密度圖（另一構型，與 s25_3 成對表示構型變化） |
| s25_5.png | 250×118 | GMM 公式：Σ_{i=1}^{N} ω_i K_σ(‖r − r_i‖) |

**建議章節**：06_workflow 課後補（heterogeneity 的 motion-based 方法：e2gmm、3DFlex）。

## Slide 26 — (1,2) Preprocessing - Motion correction [1]（補充）

**文字**：連結 https://nramm.nysbc.org/wp-content/seminars/2017/slides/nramm_2017_Rubinstein.pdf
**Speaker notes**（重點）：迭代式對齊——每次把一張 raw frame 對到「其餘所有 frame 的目前最佳總和」（leave-one-out 避免 frame「找到自己」）；每輪後對 X/Y shifts 擬合 spline 以降低雜訊敏感度；frame 總和可再套 exposure filter。

| 檔案 | 尺寸 | 內容描述 |
|---|---|---|
| s26_1.png | 446×446 | beam-induced motion 漂移軌跡圖（X-Y 平面上的藍點紅線軌跡，單位像素） |
| s26_2.png | 1352×706 | leave-one-out 對齊流程圖：Take movie → Remove a frame → Align removed frame to sum → Fit trajectories to a (smooth) spline → Calculate sum → repeat until convergence |

**建議章節**：06_workflow 課後補（motion correction 演算法細節；s26_1 也可提前用在 06 主節展示漂移現象）。
