# 第 4 章原始教材覆蓋紀錄

來源：`../04_images_processing_4_zh.ipynb`，全部 38 cells 已閱讀並恢復。原檔未修改；SHA-256：`c289b403b2c2779172f8157ec9283d436a5947ca28712443c31619e76f8f8f90`。

| 原cell | 新source行 | 教學去向與修正 |
|---:|---:|---|
| 0 | 18 | 套件由既有環境供應；移除舊imread安裝 |
| 1 | 22 | 原setup改shared support與必要imports；原下載但未用之Gengar/Pikachu/Dragonite/moonlanding/chess不增新實驗 |
| 2 | 41 | # 小波分析：從時間與尺度到影像去雜訊 上一章用 Fourier 表示頻率，並用金字塔分開影像尺度。本章把分析函數放到訊號的不同位置，再調整它的尺度，觀察局部變化。先比較時間與頻率的取捨，再 |
| 3 | 49 | 原四頻率圖保留；修正Fourier不可區分/零時間資訊與第二區間索引 |
| 4 | 66 | 原時間頻率圖保留；明示STFT/CWT取捨 |
| 5 | 81 | ## 小波如何分析訊號？ |
| 6 | 87 | 原正弦與小波圖保留；scalogram/尺度頻率/卷積共軛定義修正 |
| 7 | 100 | 原CWT公式图+可讀公式；dyadic位移隨尺度改變 |
| 8 | 119 | ## 小波家族 |
| 9 | 125 | 小波的支撐長度、平滑性、對稱性與消失動差會影響它對局部形狀的反應。[PyWavelets 文件](https://pywavelets.readthedocs.io/en/latest/)整 |
| 10 | 131 | 原families與全wavelist列出 |
| 11 | 146 | 有限能量、零平均與admissibility分開；後者收合 |
| 12 | 160 | 下面第一列是 `db5`、`sym5`、`coif5`、`bior2.4`，第二列是 `mexh`、`morl`、`cgau5`、`gaus5`。離散小波的 `wavefun()` 同時回傳 |
| 13 | 166 | 原8家族全部保留；psi取值更正，bior分析合成與cgau複數分開 |
| 14 | 198 | ### 同一家族的階數與取樣精細度 濾波器係數數目、消失動差數目與分解層數是不同量。以 Daubechies 的 `dbN` 為例，N 表示消失動差數，濾波器長度為 2N。下面用 `db1` |
| 15 | 206 | db1..5×wavefun level1..5全25圖；phi/psi更正 |
| 16 | 226 | wavefun取樣level與資料DWT level分開；係數數目/消失動差不同 |
| 17 | 234 | ## CWT：閱讀時間與尺度 |
| 18 | 240 | 原Kaggle閱讀保留，增加CWT尺度圖讀法；不加新實驗 |
| 19 | 248 | ## DWT：用濾波器分解訊號 |
| 20 | 254 | DWT 以分析濾波器組把訊號分成近似與細節，並搭配降採樣。保留全部係數時可以重建；若量化、捨棄或縮小部分係數，就能用於壓縮或去雜訊。[這則 DWT 說明](https://dsp.stack |
| 21 | 260 | 原1000Hz分頻例保留但標理想化；filterbank alias cancellation與邊界上限修正 |
| 22 | 270 | 兩張原濾波樹圖保留；dwt/wavedec係數順序 |
| 23 | 288 | 原ECG db1 smooth單階分解並驗證還原 |
| 24 | 299 | 原ECG整段重建對照 |
| 25 | 306 | 原db1 level8多階與前1000點圖，驗證還原 |
| 26 | 318 | 原閾值圖保留；soft縮小大係數而非刪大係數 |
| 27 | 339 | 原振幅比例閾值保留，periodization等價mode，保留輸入長度 |
| 28 | 352 | 原ecg/256與sigma.05，註明SD並固定seed42 |
| 29 | 359 | 原啟發式.1與BayesShrink db4 level4 |
| 30 | 367 | 原三條曲線加圖例/軸單位 |
| 31 | 378 | ## 二維 DWT：影像的四個子帶 |
| 32 | 384 | 原2D圖與Medium閱讀保留；子帶軸命名說明 |
| 33 | 397 | 原Lucario灰階+bior1.3四圖，增加idwt2還原確認 |
| 34 | 426 | 原LucarioRGB crop40:150/80:160，SD.15；五種结果+原圖2×3面板；sigma_est,/2,/4完整；seed42+PSNRrange1 |
| 35 | 497 | 原文章影片保留，補閾值選擇與PSNR判讀 |
| 36 | 509 | ## 延伸閱讀 |
| 37 | 515 | 原課本節號與Packt repo保留，補完整書名與annotated reading |

## 限制

未同步ipynb、未build、未執行全章；root負責統一執行與原圖顯示驗收。原CWT時間頻率示意全部保留，解釋以正確的解析度取捨為準。移除現有非原章主線的pyramid/blending與SPA/理解檢查。