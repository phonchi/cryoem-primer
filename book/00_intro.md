# Cryo-EM 單顆粒分析入門

這本教材寫給第一次接觸冷凍電子顯微鏡（cryo-electron microscopy, cryo-EM）單顆粒分析（single-particle analysis, SPA）的專題生。主線不是背軟體按鈕，而是先弄清楚一張影像代表什麼，再把每個計算步驟接回成像模型。

```{admonition} 讀完後，你應該能做什麼？
:class: important

- 分清楚 movie frame、motion-corrected micrograph、particle image 與 3D density map。
- 用陣列、濾波、傅立葉轉換與多尺度分析解釋常見的 SPA 操作。
- 說明投影、CTF 與雜訊如何構成簡化的 forward model。
- 生成一組保留 ground truth 的合成粒子影像，並知道它不能取代真實資料驗證。
```

## 一筆 SPA 資料怎麼走到 3D 密度圖？

直接電子偵測器先記錄一段 **dose-fractionated movie**。每個 frame 的電子劑量低，樣品也會在照射期間移動；將 frames 做 motion correction 並加總後，才得到後續處理常用的 **micrograph**。一張 micrograph 內含許多低訊噪比、經 CTF 調制的二維粒子投影，連同冰層、碳膜、污染物與偵測器雜訊。particle picking 找出座標並切出 particle images；對齊、分類與 3D refinement 再從大量不同取向的粒子影像估計 3D 密度圖。

```text
Movie frames → motion correction → micrograph
                                      ├─ CTF estimation ─┐
                                      └─ particle picking ┴→ particle images
                                              → 2D classification
                                              → 3D reconstruction／refinement
                                              → half-map validation
```

這條資料鏈很重要：movie／micrograph 是實驗觀測，particle image 是從 micrograph 擷取的分析單位，projection 則是 forward model 裡的理想化物理量。三者不能互換。

## 閱讀順序

| 章節 | 問題 | 讀完後的最低標準 |
|---|---|---|
| 1　影像基礎 | 像素、取樣與檔案到底存了什麼？ | 能解釋座標、dtype、Nyquist 與 MRC header |
| 2　濾波與分割 | 怎麼找局部特徵，又不把前處理當成真訊號？ | 能比較濾波邊界、NCC、DoG 與 watershed |
| 3　傅立葉轉換 | 為什麼 CTF 與卷積要在頻域思考？ | 能核對 DFT 慣例並區分 phase flipping 與 Wiener correction |
| 4　小波 | 固定解析度的頻譜不夠時怎麼辦？ | 能分清 CWT、DWT、filter bank 與閾值假設 |
| 5–6　SPA 背景與流程 | 各計算步驟如何接回物理與統計模型？ | 能畫出完整流程並指出主要失敗模式 |
| 7　合成資料 | 已知取向、平移、CTF 與雜訊時能驗證什麼？ | 能重跑生成流程並保存 ground truth |
| 8　延伸資源 | 下一步該查教科書、綜述還是軟體文件？ | 能依問題選對來源 |

建議依序閱讀第 1–4 章；已熟悉數位影像處理的讀者，可以先做各章末的理解檢查，再跳到第 5 章。第 7 章會用到前面所有觀念，不適合作為第一章。

## 先備知識與執行環境

讀者應熟悉 Python 基本語法、`NumPy` 陣列與高中程度的三角函數。複數、線性代數與機率若不熟，可以先讀主線；較深的推導放在可展開的進階區塊。

程式以 `cryoem-book` conda 環境執行，主要套件包括 `NumPy`、`SciPy`、`scikit-image`、`mrcfile`、`PyWavelets` 與 ASPIRE-Python。每章的 `.py` 是 Jupytext single source of truth，同步後可下載 `.ipynb` 執行。

```{admonition} 本書的主張範圍
:class: caution

教材聚焦 SPA 的計算觀念。程式中的小影像、簡化 CTF 與白雜訊是教學模型，不代表完整的顯微鏡、樣品或偵測器。視覺上「變清楚」也不等於解析度提高；任何會改變像素的增強法，都要用獨立資料與合適的驗證指標確認。
```

## 理解檢查

1. 為什麼不能把單張 movie frame、motion-corrected micrograph 與 particle image 都叫作「原始投影」？
2. particle picking 改變的是觀測資料，還是從觀測資料取出的分析單位？
3. 合成資料有 ground truth，為什麼仍不能取代真實資料的驗證？
