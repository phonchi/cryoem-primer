# 從影像處理到 Cryo-EM 單粒子分析

一張影像能告訴我們什麼，取決於它如何形成、哪些變化來自訊號，以及我們如何處理誤差。一般影像處理教我們表示、轉換與分析像素；Cryo-EM 單粒子分析（single-particle analysis, SPA）則進一步問：如何從大量含雜訊、方向未知的二維影像，推估共同的三維結構？

本書分成三個部分。**一般影像處理**四章建立通用工具；**SPA** 篇沿著實驗、成像、統計推論與結構驗證展開；最後的**合成資料**章提供可執行的例子，讓取向、CTF、平移與雜訊的作用可以逐一觀察。下面依閱讀順序列出各章要回答的問題。

## 一般影像處理：四章完整基礎

| 章節 | 主要問題 |
|---|---|
| {doc}`01_image_basics` | 像素、色彩與強度如何表示，旋轉、縮放、剪切與射影校正如何操作？ |
| {doc}`02_filter_segment` | 如何濾波、找邊緣、偵測與分割物件？ |
| {doc}`03_fourier` | 如何用原圖觀察幅度／相位、頻域濾波、金字塔與特徵偵測？ |
| {doc}`04_wavelet` | 如何同時描述位置與尺度，並用多尺度係數處理雜訊？ |

這四章可以獨立閱讀。程式以小型例子呈現方法，適合逐格執行並觀察參數效果。已熟悉這些概念的讀者，也可以直接進入 SPA 篇，在需要時再回查相關工具。

## SPA：從實驗問題走到統計解釋

SPA 篇面向修過基礎統計、微積分與線性代數的大學生。機率密度、平均數與變異數是起點；概似、後驗、EM、降維與模型驗證會隨著問題逐步介紹。讀者可先掌握正文的動機、重要公式與完整例子，再展開推導查看數學細節。

| 閱讀順序 | 核心問題 |
|---|---|
| {doc}`05_background` | 為何要用許多低劑量影像研究結構，製樣如何影響可用資料？ |
| {doc}`05_image_formation` | 三維結構如何經投影、CTF 與雜訊成為粒子影像？ |
| {doc}`05_statistical_inference` | 未知角度與位移如何用機率描述，EM 與 MAP 在估計什麼？ |
| {doc}`06_workflow` | Movie、micrograph 與粒子影像如何逐步形成分析資料？ |
| {doc}`06_alignment_classification` | 如何同時處理對齊、分類與平均？ |
| {doc}`06_dimension_reduction` | 高維影像如何用共用的低維表示描述？ |
| {doc}`06_reconstruction_validation` | 不同方向的二維資訊如何合成三維結構？ |
| {doc}`06_heterogeneity` | 多種構形如何與角度、位移及 CTF 造成的差異區分？ |
| {doc}`06_resolution_validation` | 重建中的哪些細節受到資料支持，又有哪些可能來自共同偏差？ |

閱讀每一章時，都可以問：觀測到的是什麼？哪些量還不知道？估計結果有哪些證據支持？例如，影像匹配使用的雜訊權重，也會用在三維重建；分析粒子資料時考慮的獨立性，也會影響 half-map 驗證。

## 先分清楚資料的層級

電子偵測器把一次曝光分成多個 movie frames。估計影格間的移動並加權合成後，得到 micrograph；找到粒子中心並裁切後，得到 particle images。這些影像仍帶有 CTF 調變與雜訊，需要藉由對齊、分類及重建來估計三維結構。

```text
movie frames → 位移校正與劑量加權 → micrograph
                                  ↓ CTF 估計與粒子挑選
                              particle images
                                  ↓ 對齊、分類與重建
                              3D density map
                                  ↓ 獨立資料與模型檢查
                              可解讀的結構資訊
```

前向模型中的「投影」是理想化物理量；movie frame 記錄一小段曝光時間內的觀測，micrograph 是處理後的整張視野，particle image 則是裁切後的分析單位。它們保留的資訊與處理歷程不同，解讀結果時要知道自己正在比較哪一層資料。

## 最後用合成資料串起來

{doc}`07_synthetic_data`從一張三維密度圖產生不同取向的投影，再加入 CTF、平移與雜訊。這裡有已知真值，可以比較處理前後的影像及誤差。合成模型讓因素容易分開觀察，真實資料則帶來更複雜的背景、樣品移動與構形變化；比較兩者，可以看出簡化模型解釋了哪些現象，還漏掉哪些因素。

SPA 篇的敘述、圖解與公式可以直接閱讀。要執行程式時，再設定 Python 環境。環境與 notebook 的執行方式見網站 repository 的 README。閱讀時遇到單位、旋轉或檔案術語，可查{doc}`appendix_conventions`；各章的延伸閱讀與{doc}`08_resources`提供下一步的來源。

## 延伸閱讀

- **Cheng, Y., Grigorieff, N., Penczek, P. A. 與 Walz, T.（2015）．*A Primer to Single-Particle Cryo-Electron Microscopy*. Cell, 161(3), 438–449。** [DOI](https://doi.org/10.1016/j.cell.2015.03.050)。先看 Figure 1 的整體流程，再讀製樣與影像處理段落，建立全書的實驗背景 {cite}`cheng2015primer`。
- **Singer, A. 與 Sigworth, F. J.（2020）．*Computational Methods for Single-Particle Electron Cryomicroscopy*. Annual Review of Biomedical Data Science, 3, 163–190。** [開放全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8412055/)。先讀第 1–2 節掌握資料與成像，再隨本站章節閱讀 maximum likelihood 與異質性內容 {cite}`singer2020`。
