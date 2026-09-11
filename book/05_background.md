# Cryo-EM 背景：為什麼要拍這麼多張低訊雜比影像

蛋白質的三維結構能幫助我們解釋分子如何辨識配體、催化反應或改變構形。冷凍電鏡單粒子分析（single-particle analysis, SPA）透過大量低訊雜比影像研究分子結構：先拍下數十萬張難以直接辨認細節的影像，再從中估計共同的三維結構。本章從樣品製備與電子劑量說明，為什麼需要這麼多影像，以及量測有哪些限制。

這些量測限制也決定了統計分析要處理的問題。製樣決定了哪些粒子進入資料、它們的取向如何分布；低劑量決定了每張影像能帶多少訊號。後面幾章的概似函數、對位演算法與解析度評估，處理的都是這裡產生的限制。知道限制從哪裡來，才能理解各個處理步驟為何需要這些假設。

## 三維結構圖背後的量測

結構預測可以提出可檢驗的模型；實際樣品還可能含有配體、修飾、多個組成或多種構形。模型信心描述預測本身的可信程度，實驗量測則用來判斷樣品在特定條件下呈現的狀態。兩者回答的問題不同，因此常常互相補充。

```{figure} images/pptx/s02_1.png
:width: 65%
:name: fig-protein-scale
從胺基酸序列到蛋白質三維結構的尺度示意；各序列是否已有實驗結構，仍須查閱結構資料庫。
```

```{figure} images/pptx/s02_3.png
:width: 35%
:name: fig-protein-fold
蛋白質 ribbon 表示法。模型由結構資料建立；顯微鏡量到的密度仍需經過重建與模型詮釋。
```

## 三種方法回答的問題不完全相同

| 方法 | 主要量測 | 常見優勢 | 需要留意的限制 |
|---|---|---|---|
| X 光晶體學（X-ray crystallography） | 晶體的繞射強度 | 成熟且可達原子尺度 | 樣品必須形成合適晶體；量到的強度沒有直接給出 Fourier 相位 |
| 核磁共振（NMR） | 原子核的頻譜與核間約束 | 可在溶液中研究結構與動態 | 傳統溶液 NMR 通常最適合較小蛋白質；尺寸增加時，譜線重疊與弛豫會使分析迅速變難，實際可行範圍也受標定與實驗設計影響 |
| 冷凍電鏡 SPA | 玻璃態冰中大量單粒子的低劑量投影 | 不需結晶；可用分類研究混合樣品中的離散狀態 | 單張粒子影像的訊雜比很低，取向、位移、CTF 與構形也未知，必須用統計方法共同估計 |

三種方法量到的物理量不同，也各自需要不同的樣品條件 {cite}`singer2020`。晶體學量的是整個晶格的繞射強度，因此每個量測都已經是巨大數量分子的平均；SPA 量的是單一分子的投影，平均這件事被搬到了電腦裡。這讓 SPA 能在平均前處理不同取向與構形，也增加了對位的難度：晶格原本提供的規則排列，在 SPA 中必須改由含雜訊的影像估計。

選擇方法時，解析度只是其中一項條件，還要考慮樣品大小、穩定性、可取得的量、是否能結晶，以及研究問題是否關心溶液動態或構形混合。例如，同樣想觀察蛋白質的構形變化，NMR 可追蹤溶液中的動態，SPA 則可從大量粒子中分出較穩定的離散狀態。

```{figure} images/pptx/s03_1.png
:width: 45%
:name: fig-spike
SARS-CoV-2 spike 的冷凍電鏡密度圖示例。從密度圖到生物機制或藥物設計，仍需要樣品、生化與功能證據一起判讀。
```

(sample-preparation)=
## 從溶液到玻璃態冰

製樣時，先把水溶液中的複合體放到載網上，減薄成一層薄液膜，再以足夠快的速度冷凍，使水來不及結晶。薄層厚度通常落在 100 到 800 Å 之間，視分子大小而定 {cite}`passmore2016`。水分子來不及排成晶格，形成玻璃態冰（vitreous ice），因此可避免結晶冰造成的強烈繞射與樣品擾動，也能在短時間尺度上固定分子構形。

冰的厚度同時牽動兩件事。太厚會加入額外的散射背景，壓低對比；太薄則可能把粒子排擠出視野，在冰洞中央留下沒有粒子的圓形區域。實務上希望冰只比粒子直徑略厚一些 {cite}`passmore2016`。

### 濃度算出來的粒子密度

冰層的厚度與溶液濃度一起決定了每平方微米能看到幾顆粒子。由這個估計，可以把樣品條件換算成後面需要收集的資料量。取 800 Å 厚的冰、2 mg/mL 的蛋白質溶液，分子量 1 MDa。每 1 μm² 的視野對應的體積是

$$
1\ \mathrm{\mu m^2}\times 0.08\ \mathrm{\mu m}=0.08\ \mathrm{\mu m^3}
=8\times10^{-17}\ \mathrm{L},
$$

其中的蛋白質質量為 $2\ \mathrm{g/L}\times 8\times10^{-17}\ \mathrm{L}=1.6\times10^{-16}\ \mathrm g$。除以莫耳質量 $10^{6}\ \mathrm{g/mol}$ 再乘上 Avogadro 常數，得到

$$
\frac{1.6\times10^{-16}}{10^{6}}\times 6.022\times10^{23}\approx 96
$$

顆粒子。也就是大約 100 particles/μm²；同樣濃度下換成 250 kDa 的分子，粒子數變成四倍，約 400 particles/μm² {cite}`passmore2016`。若載網上的洞直徑是 1.2 μm，單一洞的面積為 $\pi(0.6)^2\approx1.13\ \mathrm{\mu m^2}$，因此一個洞裡約有 110 顆 1 MDa 的粒子。想累積 $2\times10^5$ 張粒子影像，就需要接近兩千個洞的資料；這正是自動化資料收集成為標準配備的原因 {cite}`cheng2016`。

### 界面會挑選粒子

實際觀察到的粒子密度常和上面的計算差距很大，而差距本身就是診斷資訊。玻璃化過程中，分子會接觸許多表面：疏水的氣–水界面、非晶碳支撐膜、金屬載網，甚至吸走多餘液體的濾紙。薄冰層中的粒子在冷凍前很容易接觸氣–水界面，因此需要特別考慮它的影響 {cite}`passmore2016`。

表面吸附可能把洞裡的粒子濃縮得比溶液還密，也可能把粒子整批黏走而使洞裡幾乎沒有東西。吸附在界面上的粒子常有較低的對比、看起來比預期大，而且可能已經部分或完全變性 {cite}`passmore2016`。這些粒子仍然會被挑選程式選中，仍然會進入後面的分類與重建。

吸附也會改變取向分布。分子以特定的面貼附界面時，那個取向在資料中的比例會遠高於其他方向，形成偏好取向（preferred orientation）。{doc}`06_reconstruction_validation`會說明，取向分布直接決定三維 Fourier 空間哪些方向被取樣到；某些方向資料稀少時，重建圖在該方向會較模糊或被拉長。改變緩衝液、加入界面活性劑、調整電漿處理條件，或改用石墨烯等支撐膜，都是為了改變分子與表面的交互作用，進而改變取向分布 {cite}`passmore2016`。

因此，「接近原生狀態」應理解為樣品在這套製備條件下保留下來的狀態。製樣不只影響影像好不好看，它決定了我們手上這批粒子是母體的哪一個子集合；統計推論的所有結論，都是針對這個子集合而言。

```{figure} images/pptx/s04_1.png
:width: 55%
:name: fig-imaging-principle
SPA 的幾何近似：不同取向的同類粒子各自形成一張 2D 投影。真實影像還會受到 CTF、位移、偵測器與雜訊影響。
```

穿透式電子顯微鏡（transmission electron microscope, TEM）讓電子束穿過冰層。對薄的生物樣品，可把樣品近似為弱相位物體：樣品主要改變電子波的相位。未散射波是參考波，會與彈性散射波經物鏡後干涉，形成相位對比；有限電子計數則帶來 shot noise。更完整的推導見 {doc}`05_image_formation`。

(low-dose-snr)=
## 低劑量保護樣品，也限制單張影像

電子是游離輻射。入射電子與樣品發生非彈性散射時會把能量留在樣品裡，破壞化學鍵、產生自由基，分子結構隨曝光累積而改變 {cite}`glaeser2016`。這個效應強到什麼程度，可以直接換算出來。

對玻璃態冰與 300 keV 的電子，單位質量的線性能量轉移約為 $2.4\ \mathrm{MeV\,cm^2/g}$ {cite}`glaeser2016`。$10\ \mathrm{e^-/\mathring A^2}$ 的曝光相當於 $10^{17}$ 個電子每平方公分，若先忽略二次電子把能量帶出樣品，可估計能量轉移的量級：

$$
2.4\ \mathrm{\frac{MeV\,cm^2}{g}}\times10^{17}\ \mathrm{cm^{-2}}
=2.4\times10^{17}\ \mathrm{MeV/g}
\approx3.8\times10^{7}\ \mathrm{J/kg}.
$$

也就是約 $3.8\times10^{7}\ \mathrm{Gy}$，或 $3.8\times10^{9}\ \mathrm{rad}$。這是忽略二次電子逸出的量級估計；薄樣品有部分能量被帶走，Glaeser 指出實際沉積劑量可依厚度下降一半或更多。即使考慮這個修正，生物分子的高解析度資訊仍會隨曝光損失，玻璃態冰本身也會受到游離輻射破壞 {cite}`glaeser2016`。

於是電子劑量變成一種預算。總 fluence 決定整次曝光收到多少電子；提高它可以改善每張影像的統計品質，也會加深輻射損傷。早期的低劑量成像例子約使用 $20\ \mathrm{e^- / \mathring{A}^2}$ {cite}`sigworth2016`；現代 SPA 常把數十 $\mathrm{e^- / \mathring{A}^2}$ 的總 fluence 分散到 movie frames。把同樣的總劑量拆成多張 frames，不會增加電子數，但增加了時間軸，讓我們能估計曝光期間的移動，也能依累積劑量調整各頻率的權重。每個 frame 的劑量描述時間分配，dose weighting 則決定合併影格時各頻率占多少權重；細節見 {doc}`05_image_formation` 與 {doc}`06_workflow`。

以下一律用功率比定義純量訊雜比（signal-to-noise ratio, SNR）：

$$
\mathrm{SNR}=\frac{P_{\mathrm{signal}}}{P_{\mathrm{noise}}},
\qquad
\mathrm{SNR}_{\mathrm{dB}}=10\log_{10}(\mathrm{SNR}).
$$

例如 $-20\ \mathrm{dB}$ 對應 $\mathrm{SNR}=0.01$，也就是雜訊功率為訊號功率的一百倍。估計區域、是否扣除平均與所選頻帶都會改變數值，計算時要採用相同設定。另外要區分兩種常見寫法：若某份資料說「SNR = 0.1」指的是振幅標準差比，換算成功率比是 $0.1^2=0.01$，相差一個數量級。比較數字以前，先確認公式、估計區域與 clean target 的定義。

頻率相依的 spectral SNR（SSNR）則寫成 $\mathrm{SSNR}(s)$。一張影像可能保有清楚的低頻輪廓，高頻訊號卻已接近雜訊水準；單一 SNR 數字只概括指定範圍內的功率比，把這個頻率差異壓成一個數。

```{figure} images/pptx/s04_4.png
:width: 80%
:name: fig-raw-denoise
低訊雜比 micrograph 與處理後示例。處理可改善觀看效果，影像中未被量到的頻率資訊仍無法由此找回。
```

## 平均為什麼有效，以及什麼時候效果打折

SPA 利用平均累積微弱訊號。平均能改善多少訊雜比，取決於影像是否對齊、訊號是否相同，以及雜訊是否獨立；以下用統計模型說明。

設每張影像有 $d$ 個像素，第 $i$ 張已正確對位的影像為 $y_i=s+\varepsilon_i$，其中 $s$ 是所有粒子共享的訊號，$\varepsilon_i$ 的各像素是零平均、變異數 $\sigma^2$ 的雜訊。以下變異數式逐像素成立；訊號的平均像素功率定義為 $P_s=\|s\|^2/d$。平均之後

$$
\bar y=\frac1N\sum_{i=1}^{N}y_i=s+\frac1N\sum_{i=1}^{N}\varepsilon_i .
$$

訊號原封不動；雜訊項的期望為零。若各張影像的雜訊互相獨立，

$$
\operatorname{Var}\!\left(\frac1N\sum_i\varepsilon_i\right)=\frac{\sigma^2}{N}.
$$

雜訊的標準差因此降為 $\sigma/\sqrt N$。訊號振幅與雜訊標準差的比值改善 $\sqrt N$ 倍，而功率比（本書採用的 SNR 定義）改善 $N$ 倍：

$$
\mathrm{SNR}(\bar y)=\frac{P_s}{\sigma^2/N}=N\cdot\mathrm{SNR}(y_i).
$$

兩個倍率描述同一件事，只是一個用振幅、一個用功率。看到「平均可讓訊雜比改善 $\sqrt N$」或「改善 $N$ 倍」時，先確認對方用哪一種比值。

### 雜訊帶有共同成分時

上面的變異數推導採用了影像間雜訊獨立的假設。實際資料裡，同一張 micrograph 的粒子共享冰層背景、同一部 movie 的粒子共享殘餘運動，偵測器與內插也會在鄰近像素之間造成相關。把這件事寫進模型最簡單的方式，是讓每個雜訊含有一個共同成分：

$$
\varepsilon_i=\sqrt{\rho}\,u+\sqrt{1-\rho}\,v_i,
$$

其中 $0\leq\rho\leq1$，$u,v_1,\ldots,v_N$ 彼此獨立，都是零平均、變異數 $\sigma^2$ 的隨機量。$u$ 是所有影像共用的成分。任兩張影像的雜訊相關係數為 $\rho$，而平均的變異數變成

$$
\operatorname{Var}(\bar\varepsilon)
=\rho\sigma^2+\frac{(1-\rho)\sigma^2}{N}
=\frac{1+(N-1)\rho}{N}\,\sigma^2 .
$$

取 $N=10^4$、$\rho=0.01$。獨立假設會預測 $\operatorname{Var}=10^{-4}\sigma^2$；實際值是

$$
\frac{1+9999\times0.01}{10^4}\sigma^2
=0.010099\,\sigma^2,
$$

大了一百倍。換句話說，這一萬張影像的效果只相當於

$$
N_{\mathrm{eff}}=\frac{N}{1+(N-1)\rho}=\frac{10^4}{100.99}\approx 99
$$

張獨立影像。$N\to\infty$ 時 $N_{\mathrm{eff}}\to1/\rho=100$，再加影像也無法突破。即使相關係數只有 1%，增加資料量的效益也會受到明確限制。這解釋了為什麼獨立性在 SPA 裡是反覆出現的主題：從 half-set 的切分方式，到粒子是否來自同一張 micrograph，都需要考慮如何讓 $\rho$ 盡量小。

### 其他讓增益打折的因素

即使雜訊真的獨立，平均的增益還需要「所有影像共享同一個 $s$」。粒子若混有不同構形，平均得到的是混合物；取向或位移估計錯誤，等於在平均前先把影像各自推開一點，高頻細節會互相抵消；CTF 沒有處理好，不同影像在同一頻率的正負號可能相反，加起來反而變小。{doc}`06_alignment_classification` 會把對位誤差寫成一個隨頻率衰減的因子，{doc}`06_heterogeneity` 則處理構形混合。

```{figure} images/pptx/s04_5.png
:width: 40%
:name: fig-atomic-model
大量低訊雜比投影經過對位與重建後，可用來建立 3D 密度圖和原子模型。密度圖的解析度描述可分辨的空間尺度，原子模型是否放對則要另外檢查。
```

## 接下來

從這些例子可看出，樣品製備決定了資料的母體與取向分布，電子劑量決定了單張影像的統計品質，大量平均則在獨立性與同質性成立時把訊號累積起來。下一步要把「一張粒子影像長什麼樣」寫成精確的數學式，包括電子波如何產生對比、CTF 如何調變各個空間頻率，以及偵測器保留了多少資訊。請接著閱讀 {doc}`05_image_formation`。

## 延伸閱讀

- **Cheng, Y., Grigorieff, N., Penczek, P. A. 與 Walz, T.（2015）．*A Primer to Single-Particle Cryo-Electron Microscopy*. Cell, 161(3), 438–449。** [DOI](https://doi.org/10.1016/j.cell.2015.03.050)。先看 Figure 1 的整體流程圖與製樣段落，再讀影像處理概述；適合在本章之後補上生物學動機與完整名詞對照 {cite}`cheng2015primer`。
- **Passmore, L. A. 與 Russo, C. J.（2016）．*Specimen Preparation for High-Resolution Cryo-EM*. Methods in Enzymology, 579, 51–86。** [DOI](https://doi.org/10.1016/bs.mie.2016.04.011)。重點讀第 2.3 節 Diagnostic Cryo-EM（頁 56–59），其中有本章使用的粒子密度計算、氣–水界面吸附與變性的討論；第 7 節 Vitrification（頁 76–80）說明玻璃化步驟；第 2.4 節 Initial Cryo-EM Data Collection（頁 59–60）談如何用初步資料檢查取向分布是否足以支撐各向同性的三維重建，以及被大量捨棄的粒子多半在製樣階段就已受損 {cite}`passmore2016`。
- **Glaeser, R. M.（2016）．*Specimen Behavior in the Electron Beam*. Methods in Enzymology, 579, 19–50。** [DOI](https://doi.org/10.1016/bs.mie.2016.04.010)。第 2.3 節（頁 24–25）給出本章使用的 $2.4\ \mathrm{MeV\,cm^2/g}$ 與劑量換算，第 3 節（頁 25–30）說明結構損傷的證據與敏感殘基，第 6 節（頁 37–40）解釋 beam-induced motion 的來源 {cite}`glaeser2016`。
- **Sigworth, F. J.（2016）．*Principles of Cryo-EM Single-Particle Image Processing*. Microscopy, 65(1), 57–67。** [Oxford Academic 免費全文](https://doi.org/10.1093/jmicro/dfv370)。開頭兩頁（頁 57–58）用 TRPV1 的資料說明單張粒子影像的 spectral SNR 為何在 16 Å 以下就低於 1，可與本章的平均論證並讀；這是一篇篇幅短、適合第一次接觸的公開文獻 {cite}`sigworth2016`。
- **Cheng, A., Tan, Y. Z., Dandey, V. P., Potter, C. S. 與 Carragher, B.（2016）．*Strategies for Automated CryoEM Data Collection Using Direct Detectors*. Methods in Enzymology, 579, 87–102。** [DOI](https://doi.org/10.1016/bs.mie.2016.04.008)。說明要收集本章估算出的資料量時，自動化如何安排載網位置、對焦與曝光；讀第 1 與第 2 節即可建立概念 {cite}`cheng2016`。
- **Singer, A. 與 Sigworth, F. J.（2020）．*Computational Methods for Single-Particle Electron Cryomicroscopy*. Annual Review of Biomedical Data Science, 3, 163–190。** [DOI](https://doi.org/10.1146/annurev-biodatasci-021020-093826)。第 1.2 節比較 X 光晶體學、NMR 與冷凍電鏡的量測對象，第 7.1 節討論製樣與蛋白質變性如何限制目前的方法 {cite}`singer2020`。
