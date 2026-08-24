# Cryo-EM 單粒子分析：從 movie 到 3D 密度圖

```{admonition} 讀完本章，你應該能
:class: important
- 分清楚 movie、micrograph、particle stack 與 3D volume 四個資料層級。
- 說明 motion correction、CTF estimation、particle picking、classification 與 refinement 各自在估計什麼。
- 分辨簡化的 hard projection matching 與 maximum-likelihood／Bayesian soft assignment。
- 知道清楚的 class average 或漂亮的 3D 圖不是單獨成立的驗證證據。
```

(workflow)=
## SPA 主線與資料層級

單粒子分析（single-particle analysis, SPA）的資料通常沿著這條路徑前進：

**movie frames → motion-corrected micrograph → particle stack → 2D／3D estimates → 3D density map**。

Movie 是一次曝光拆成的多個時間 frame；micrograph 是 frame 對齊並加權合成後的整張視野；particle stack 則是從 micrograph 依座標裁出的粒子小圖。三者不能混稱為「原始影像」。

```{figure} images/pptx/s05_1.png
:width: 55%
:name: fig-workflow-overview
SPA 主流程：Movie Alignment → CTF Estimation → Particle Picking → 2D Classification → Initial Model → 3D Classification → 3D Refinement。
```

[Computational-CryoEM](https://github.com/phonchi/Computational-CryoEM) 整理了各處理階段的論文、軟體與教學資源。本章只走 SPA；tomographic reconstruction 屬於 cryo-electron tomography 的分支，不是一般 SPA 的必經階段。

一個資料集可含數百到數千部 movie，以及 $10^5$ 到 $10^6$ 張粒子影像。除了檔案量大，每張粒子的三維取向、平面內旋轉、位移、CTF 與構形也可能未知，且單張影像的訊雜比很低 {cite}`singer2020`（pp. 171–180）。

```{figure} images/pptx/s10_2.png
:width: 70%
:name: fig-particle-grid
裁切後的粒子影像。肉眼看不清單張高頻細節，不代表整批資料沒有可估計的共同訊號。
```

## Motion correction：先處理曝光期間的移動

直接電子偵測器（direct electron detector, DED）可以把一次曝光記成 movie。電子束會引起整體漂移與局部、非均勻的 beam-induced motion；若直接加總 frames，高頻資訊會因位移而模糊。

常見流程先以互相關估計全域或 patch 位移，再對空間與時間上的軌跡做平滑，最後對齊加總。高解析訊號通常在曝光早期受損較少，但最早的 frames 又可能有較劇烈的運動，因此實務上還會依空間頻率做 dose weighting。不同軟體的運動模型與正則化不同，motion-corrected micrograph 不是唯一、無誤差的答案 {cite}`rubinstein2016,singer2020`（Ripstein & Rubinstein, pp. 103–124；Singer & Sigworth, pp. 170–172）。

```{figure} images/pptx/s06_2.png
:width: 30%
:name: fig-frame-raw
單一 movie frame 的電子計數有限，訊雜比通常很低。
```

```{figure} images/pptx/s06_3.png
:width: 40%
:name: fig-frame-aligned
frames 經對齊與加權後形成 micrograph。可見度改善不等於所有局部運動都已完全校正。
```

**前置概念**：{doc}`02_filter_segment` 的互相關與 {doc}`03_fourier` 的 Fourier shift theorem。

## CTF estimation：估計每張影像如何轉移對比

離焦（defocus）讓弱相位物體產生可見對比，也使 CTF 隨空間頻率振盪、反相並出現零點。CTF 估計通常從 micrograph 的 2D periodogram／power spectrum 擬合 defocus、散光與其他參數。無散光近似或初始估計可看徑向平均；要估散光大小與方向，必須保留二維的非圓對稱資訊 {cite}`singer2020`（p. 172）。

三種處理角色要分開：

- **Phase flipping** 乘上 $\operatorname{sign}(H)$，只修正對比正負號，不恢復振幅。
- **Wiener-style correction** 以 SSNR 正則化振幅校正，目標與假設都不同。
- **多離焦合併** 利用不同影像的 CTF 零點錯開，補上單張影像在零點沒有量到的頻率。

沒有任何單張影像校正能恢復 CTF 零點處已遺失的資訊。較大的離焦通常提高低頻對比，但也讓高頻 CTF 振盪更快、零點更密；高頻包絡衰減還受到 temporal／spatial coherence、運動與輻射損傷等因素影響，不能全歸因於離焦。

```{figure} images/pptx/s08_1.png
:width: 75%
:name: fig-ctf-theory
不同離焦下的對比轉移示意。曲線的振盪與 envelope attenuation 是不同機制。
```

```{figure} images/pptx/s08_2.png
:width: 55%
:name: fig-ctf-fit
實驗 2D power spectrum 與 CTF 擬合的四象限對照；非圓對稱的 Thon rings 可攜帶散光資訊。
```

**前置概念**：{doc}`03_fourier` 的卷積定理、phase flipping 與 Wiener filter；成像模型見 {doc}`05_image_formation`。

## Particle picking：找候選座標，不是宣告真實粒子

Particle picking 需要在 micrograph 找出候選中心並避開污染。常見方法包括簡單 blob／DoG template、以資料統計量驅動的無監督方法，以及 Topaz、crYOLO 等監督式模型。詳細參考模板若來自同一批低訊雜比資料，可能把 template bias 帶回結果；候選數量也會受閾值、最小距離、冰厚與粒子大小影響 {cite}`sigworth2016`（p. 62）。

```{figure} images/pptx/s09_1.png
:width: 60%
:name: fig-picking
軟體標出的候選粒子位置仍須經裁切、分類與資料品質檢查。
```

## 2D classification：同時估計對位與類別

2D classification 會把影像在平面內旋轉、平移並分群，再以加權平均估計 class averages。平均能提升共同訊號的相對強度，但前提是分到同一類的粒子確實可比較。

模糊或條紋狀 class average 可能來自污染與誤挑，也可能是錯誤對位、低 occupancy、preferred orientation 或構形異質性。輪廓清楚的 class 也不能單憑外觀保證正確；挑選或剔除應搭配粒子數、方向覆蓋、重複分析與下游驗證。

```{figure} images/pptx/s12_2.png
:width: 22%
:name: fig-class-bad
品質可疑的 2D class average。外觀指出需要檢查，不直接指定唯一原因。
```

```{figure} images/pptx/s12_4.png
:width: 22%
:name: fig-class-good
輪廓較清楚的 2D class average。是否納入 3D 分析仍需看粒子數與其他診斷。
```

{doc}`04_wavelet` 的多尺度觀察可幫助理解不同頻帶的結構，但 class averaging 的 SNR 提升來自對位後的統計平均，不是小波閾值去雜訊。

(refinement)=
(heterogeneity)=
## Initial model、3D classification 與 refinement

Fourier slice theorem 說明了已知取向的 2D 投影如何約束 3D Fourier volume。SPA 的取向未知，因此這是一個非線性反問題。最簡化的 projection matching 會從目前模型產生許多參考投影，為每張影像挑出相關最高的方向與位移，再重建新模型。這個 hard assignment 適合建立直覺，卻不是現代 maximum-likelihood／Bayesian refinement 的完整描述。

以 RELION 類方法為例，演算法會對候選方向、位移與類別計算後驗權重，更新時用機率加權，不必把每張低訊雜比影像過早鎖定在單一姿態。先驗與正則化會抑制資料不足頻帶中的不穩定解；它們也會影響結果，因此不能把收斂等同於真實 {cite}`scheres2010,singer2020`（Sigworth et al., pp. 273–277；Singer & Sigworth, pp. 175–180）。

3D classification 用離散類別描述構形混合，但類別數要事先選擇，小族群可能因訊號不足而漏掉。連續異質性方法屬於進階延伸；本書不展開 cryo-ET、helical reconstruction 或 atomic fitting。

```{figure} images/pptx/s13_3.png
:width: 35%
:name: fig-3d-map
70S ribosome 的 3D 密度圖示例。地圖外觀、解析度估計與生物詮釋是三個不同層次的問題。
```

完整的 Fourier slice、重建、gold-standard half-set 與 FSC 限制見 {doc}`06_reconstruction_validation`；可控的前向模型見 {doc}`07_synthetic_data`。

## 理解檢查

1. 為什麼徑向平均後的 power spectrum 不足以估計散光方向？
2. Phase flipping、Wiener-style correction 與多離焦合併各自處理什麼問題？
3. Hard projection matching 在低訊雜比資料中會忽略哪一類不確定性？
4. 一個 class average 看起來模糊時，至少列出三種可能原因，以及可用來區分它們的證據。
