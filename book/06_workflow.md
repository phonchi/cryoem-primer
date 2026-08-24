# Cryo-EM 單粒子分析：從 movie 到 3D 密度圖

```{admonition} 讀完本章，你應該能
:class: important
- 分清楚 movie、micrograph、particle stack 與 3D volume 四個資料層級。
- 說明 motion correction、CTF estimation、particle picking、classification 與 refinement 各自在估計什麼。
- 分辨簡化的 hard projection matching 與 maximum-likelihood／Bayesian soft assignment。
- 說明為何還要配合粒子數、方向分布與獨立驗證，才能解讀清楚的類別平均影像或 3D 密度圖。
```

(workflow)=
## SPA 主線與資料層級

單粒子分析（single-particle analysis, SPA）的資料通常依序經過以下階段：

**movie frames → motion-corrected micrograph → particle stack → 2D／3D estimates → 3D density map**。

Movie 是一次曝光拆成的多個時間 frame；micrograph 是 frame 對齊並加權合成後的整張視野；particle stack 則是從 micrograph 依座標裁出的粒子小圖。請保留這三個名稱，才能看出資料位於哪個處理階段。

```{figure} images/pptx/s05_1.png
:width: 55%
:name: fig-workflow-overview
SPA 主流程：Movie Alignment → CTF Estimation → Particle Picking → 2D Classification → Initial Model → 3D Classification → 3D Refinement。
```

[Computational-CryoEM](https://github.com/phonchi/Computational-CryoEM) 整理了各處理階段的論文、軟體與教學資源。以下沿著 SPA 的資料流往下看；冷凍電子斷層掃描使用另一條斷層重建（tomographic reconstruction）路徑。

一個資料集可含數百到數千部 movie，以及 $10^5$ 到 $10^6$ 張粒子影像。除了檔案量大，每張粒子的三維取向、平面內旋轉、位移、CTF 與構形也可能未知，且單張影像的訊雜比很低 {cite}`singer2020`。

```{figure} images/pptx/s10_2.png
:width: 70%
:name: fig-particle-grid
裁切後的粒子影像。單張影像的高頻細節很難用肉眼辨認；大量粒子經正確對位後，仍可估計共同訊號。
```

## Motion correction：先處理曝光期間的移動

直接電子偵測器（direct electron detector, DED）可以把一次曝光記成 movie。電子束會引起整體漂移與局部、非均勻的 beam-induced motion；若直接加總 frames，高頻資訊會因位移而模糊。

常見流程先以互相關估計全域或區塊位移，再對空間與時間上的軌跡做平滑，最後對齊加總。高解析訊號通常在曝光早期受損較少，但最早的 frames 又可能有較劇烈的運動，因此實務上還會依空間頻率做 dose weighting。若校正後的 micrograph 仍有方向一致的拖影，或不同區域的細節清晰度差很多，可以回頭檢查局部運動模型與軌跡是否合適 {cite}`rubinstein2016,singer2020`。

```{figure} images/pptx/s06_2.png
:width: 30%
:name: fig-frame-raw
單一 movie frame 的電子計數有限，訊雜比通常很低。
```

```{figure} images/pptx/s06_3.png
:width: 40%
:name: fig-frame-aligned
frames 經對齊與加權後形成 micrograph。可見度改善後，仍要檢查局部運動是否留下殘差。
```

**前置概念**：{doc}`02_filter_segment` 的互相關與 {doc}`03_fourier` 的 Fourier shift theorem。

## CTF estimation：估計每張影像如何轉移對比

離焦（defocus）讓弱相位物體產生可見對比，也使 CTF 隨空間頻率振盪、反相並出現零點。CTF 估計通常從 micrograph 的 2D periodogram／power spectrum 擬合 defocus、散光與其他參數。無散光近似或初始估計可看徑向平均；要估散光大小與方向，必須保留二維的非圓對稱資訊 {cite}`singer2020`。

以下三種處理的目的不同：

- **Phase flipping** 乘上 $\operatorname{sign}(H)$，作用限於修正對比正負號；CTF 壓低的振幅仍維持原狀。
- **Wiener-style correction** 以 SSNR 正則化振幅校正，目標與假設都不同。
- **多離焦合併** 利用不同影像的 CTF 零點錯開，補上單張影像在零點沒有量到的頻率。

單張影像在 CTF 零點處沒有可用資訊，校正方法無從恢復。較大的離焦通常提高低頻對比，同時也讓高頻 CTF 振盪更快、零點更密。高頻包絡衰減還受到 temporal／spatial coherence、運動與輻射損傷等因素影響，離焦只是其中一項。

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

## Particle picking：找出候選座標

Particle picking 需要在 micrograph 找出候選中心並避開污染。常見方法包括 blob／DoG 偵測、模板匹配、依資料統計特徵建立的無監督方法，以及 Topaz、crYOLO 等監督式模型。詳細參考模板若來自同一批低訊雜比資料，可能把 template bias 帶回結果；候選數量也會受閾值、最小距離、冰厚與粒子大小影響 {cite}`sigworth2016`。

把挑選座標疊回 micrograph，是最快的初步檢查。若座標大量落在碳膜、冰晶或污染上，應先調整閾值或訓練資料；若同一顆粒子附近出現多個座標，則要檢查最小距離。接著再看裁切框，確認粒子大致位於中央，而且 box size 足以容納完整投影。

```{figure} images/pptx/s09_1.png
:width: 60%
:name: fig-picking
軟體標出的候選粒子位置仍須經裁切、分類與資料品質檢查。
```

## 2D classification：同時估計對位與類別

2D classification 會把影像在平面內旋轉、平移並分群，再以加權平均估計 class averages。平均能提升共同訊號的相對強度，但前提是分到同一類的粒子確實可比較。

模糊或帶有條紋的類別平均影像，可能來自污染、誤挑、錯誤對位、該類粒子數太少、偏好取向或構形異質性。輪廓清楚時，仍要搭配粒子數、方向分布、重複分析與後續驗證，才能決定是否保留該類粒子。

```{figure} images/pptx/s12_2.png
:width: 22%
:name: fig-class-bad
品質可疑的 2D 類別平均影像。這種外觀提醒我們檢查資料，成因仍要靠其他資訊判斷。
```

```{figure} images/pptx/s12_4.png
:width: 22%
:name: fig-class-good
輪廓較清楚的 2D 類別平均影像。是否納入 3D 分析，還要查看粒子數與其他檢查結果。
```

{doc}`04_wavelet` 的多尺度觀察可幫助理解不同頻帶的結構。類別平均影像的訊雜比提升來自對位後的統計平均；小波閾值去雜訊採用另一套假設與運算。

(refinement)=
(heterogeneity)=
## Initial model、3D classification 與 refinement

Fourier slice theorem 說明了已知取向的 2D 投影如何約束 3D Fourier volume。SPA 的取向未知，因此這是一個非線性反問題。最簡化的 projection matching 會從目前模型產生許多參考投影，為每張影像挑出相關最高的方向與位移，再重建新模型。這個 hard assignment 適合建立直覺；現代 maximum-likelihood／Bayesian refinement 會保留多個候選姿態的不確定性。

以 RELION 這類方法為例，演算法會對候選方向、位移與類別計算後驗權重，更新時用機率加權，讓每張低訊雜比影像保留多個可能姿態。先驗與正則化會抑制資料不足頻帶中的不穩定解，也會影響最後結果。目標函數趨於穩定時，演算法便會停止更新。這時應繼續查看方向分布、half-map 與 FSC；它們回答的問題和收斂條件不同 {cite}`scheres2010,singer2020`。

3D classification 用離散類別描述構形混合，但類別數要事先選擇，小族群可能因訊號不足而漏掉。連續異質性、cryo-ET、helical reconstruction 與 atomic fitting 可作為後續主題。

```{figure} images/pptx/s13_3.png
:width: 35%
:name: fig-3d-map
70S ribosome 的 3D 密度圖示例。地圖外觀、解析度估計與生物詮釋是三個不同層次的問題。
```

下一章會從 Fourier slice 接到重建、gold-standard half-set 與 FSC：{doc}`06_reconstruction_validation`。若想逐項改變取向、平移、CTF 與雜訊，可接著操作 {doc}`07_synthetic_data`。

## 延伸閱讀

- 想了解整條 SPA 分析流程，可讀 Singer 與 Sigworth 的綜述 {cite}`singer2020`。
- 想進一步了解 maximum-likelihood refinement，可讀 Scheres 等人的方法說明 {cite}`scheres2010`。
- 實際操作可參考 [RELION](https://relion.readthedocs.io/en/release-5.0/) 或 [cryoSPARC](https://guide.cryosparc.com/) 的官方教學。

## 理解檢查

1. 為什麼徑向平均後的 power spectrum 不足以估計散光方向？

```{dropdown} 參考答案
散光會讓離焦量隨 Fourier 平面的方位角改變，因此 Thon rings 會呈現橢圓或其他非圓對稱形狀。徑向平均把同一半徑上的方位角資訊合併，只能保留接近平均離焦的徑向變化，散光方向也隨之消失。

因此，清楚的徑向峰谷只能用於無散光近似或初始估計；估計散光大小與方向仍需使用二維 power spectrum。
```

2. Phase flipping、Wiener-style correction 與多離焦合併各自處理什麼問題？

```{dropdown} 參考答案
Phase flipping 乘上 $\operatorname{sign}(H)$，修正 CTF 造成的對比反相。Wiener-style correction 常寫成

$$
\widehat F_{\mathrm W}=\frac{H^*}{|H|^2+K}\widehat I,
$$

它同時校正對比方向與振幅，並用 $K$ 抑制 CTF 接近零時的雜訊放大；$K$ 的選擇需要訊號與雜訊的先驗資訊。多離焦合併則利用不同影像的 CTF 零點互相錯開，讓某張影像缺少的頻率由其他影像提供。

判讀時先問要修正的是正負號、振幅，還是零點處缺少的頻率。單張影像在 CTF 零點遺失的資訊，只能由其他離焦影像補充。
```

3. Hard projection matching 在低訊雜比資料中會忽略哪一類不確定性？

```{dropdown} 參考答案
低訊雜比下，相關最高與次高的方向、位移或類別可能得到很接近的分數。Hard projection matching 只保留最高分候選，因此丟掉「這張影像其實也可能來自其他姿態」的不確定性。

機率加權方法會保留多個候選，例如

$$
p(z\mid I,V)\propto p(I\mid z,V)p(z),
$$

再以後驗機率更新模型。最高相關分數只標出目前候選中的最佳者；Soft assignment 的權重也會隨初始模型、候選網格、雜訊模型與先驗設定而改變。
```

4. 一個 class average 看起來模糊時，至少列出三種可能原因，以及可用來區分它們的證據。

```{dropdown} 參考答案
可能原因包括污染或誤挑、平面內對位錯誤、該類粒子數太少、同一類混入不同構形，以及偏好取向造成的資料分布不均。可依序查看原始粒子框與挑選座標、旋轉與位移分布、每類粒子數、不同初始值下的分類穩定性，以及後續 3D 的方向分布。

一張模糊的類別平均影像通常無法單獨指出成因。輪廓清楚時，表示該類具有一致的二維訊號；是否把這些粒子帶入 3D 分析，則要合併上述檢查結果判斷。
```
