# Cryo-EM 單粒子分析：從 movie 到 3D 密度圖

{doc}`05_statistical_inference` 把 SPA 寫成一個統計問題：觀測影像 $y_i$、未知角度、位移與可能類別 $z_i$、共享結構 $V$，以及一組雜訊假設。但那組概似函數並非直接作用在顯微鏡輸出的原始檔案上。在統計推論開始之前，資料要先經過幾道處理，把一部 movie 變成一疊可以互相比較的粒子影像。

本章依資料的處理順序，說明各步驟估計的量、輸入與輸出的資料層級，以及誤差或偏差如何傳到後續分析。前處理決定了統計模型實際分析的資料；{doc}`06_resolution_validation` 提到的許多偏差，源頭都在這裡。

(workflow)=
## SPA 主線與資料層級

單粒子分析（single-particle analysis, SPA）的資料通常依序經過以下階段：

**movie frames → motion-corrected micrograph → particle stack → 2D／3D estimates → 3D density map**。

Movie 是一次曝光拆成的多個時間 frame；micrograph 是 frame 對齊並加權合成後的整張視野；particle stack 則是從 micrograph 依座標裁出的粒子小圖。區分這三個名稱，就能看出資料位於哪個處理階段，並判斷應在哪一層檢查問題。

```{figure} images/pptx/s05_1.png
:width: 55%
:name: fig-workflow-overview
SPA 主流程：Movie Alignment → CTF Estimation → Particle Picking → 2D Classification → Initial Model → 3D Classification → 3D Refinement。
```

一個資料集可含數百到數千部 movie，以及 $10^5$ 到 $10^6$ 張粒子影像。除了檔案量大，每張粒子的三維取向、平面內旋轉、位移、CTF 與構形也可能未知，且單張影像的訊雜比很低 {cite}`singer2020`。資料量大、未知量多且訊雜比低，使得所有步驟都必須能在大量資料上自動執行，而且每一步都要假設自己的輸入含有大量雜訊與一部分錯誤。

```{figure} images/pptx/s10_2.png
:width: 70%
:name: fig-particle-grid
裁切後的粒子影像。單張影像的高頻細節很難用肉眼辨認；大量粒子經正確對位後，仍可估計共同訊號。
```

(motion-correction)=
## 運動估計：三個層級的同一個問題

直接電子偵測器把一次曝光記成 movie，因為電子束會引起樣品的整體漂移與局部、非均勻的移動。若直接加總 frames，高頻資訊會因位移而模糊。運動校正的工作是估計每個 frame（或每個區塊、每顆粒子）相對於參考位置的平移，再把它補回去。

這三個層級都可以利用影像間的一致性估計位移；光流或直接最佳化合併目標也是可用的途徑，以下先用互相關說明。兩張影像若只差一個平移 $\mathbf t$，它們的互相關函數會在 $\mathbf t$ 處出現峰值。實作上通常先把一張影像的頻譜乘上另一張的共軛，再做逆 Fourier 變換，一次得到所有候選平移的分數。{doc}`02_filter_segment` 與 {doc}`03_fourier` 已建立這兩項工具。

### Frame 層級

最早的做法把整張 frame 當成剛體。以 $T$ 個 frames 為例，可以計算所有不重複配對的互相關函數，共 $T(T-1)/2$ 個，每個給出一個量測到的相對位移 $\mathbf m_{mn}$。若真正的逐 frame 位移為 $\mathbf s_{m,m+1}$，則

$$
\mathbf m_{mn}=\sum_{l=m}^{n-1}\mathbf s_{l,l+1}
$$

在無雜訊時成立。這是一組線性方程式：$T=30$ 時有 435 條方程式、29 個未知量，因此系統高度超定，可用最小平方求解 {cite}`rubinstein2016`。單一配對的互相關峰值可能不穩定；利用幾百條彼此一致的約束求解，可降低個別量測誤差的影響。

計算時還要處理偵測器圖樣與峰值定位的問題。計算互相關以前通常先套一個抑制高頻的濾波器，避免偵測器的固定圖樣（fixed-pattern noise）主導峰值；相鄰或接近相鄰的 frame 配對，因位移太小而容易被原點附近的殘餘峰值淹沒，可以先排除；峰值位置可用 Fourier 補零內插到次像素精度 {cite}`rubinstein2016`。

### Patch 層級

整張 frame 的剛體位移只是近似。冰層在電子束下會像鼓面一樣彎曲，視野不同位置的移動方向與大小可以明顯不同 {cite}`glaeser2016`。把 frame 切成區塊各自對齊，可以描述這種空間變化；但區塊愈小，訊號愈弱，估計愈不穩定。

常見做法是加入平滑假設：軌跡在時間上與空間上都應該連續。實作上把量測到的位移擬合成時間的樣條或多項式，同時要求鄰近區塊的軌跡相似，再把平滑後的軌跡套回去，反覆迭代直到位移更新小於某個門檻（例如 0.5 pixel）{cite}`rubinstein2016`。平滑參數因此是一個取捨：太強會壓掉真實的局部運動，太弱則把雜訊當成運動。

### Particle 層級

這是取得初步三維重建後，回到原始 frames 改善粒子軌跡的精修步驟。第一次處理資料可先做 frame／patch 校正，再經挑選、分類與重建取得參考模型；流程因此會回頭迭代。下面的 FSC 與 B factor 用來比較各影格保留的資訊，其定義與 Guinier plot 的讀法見{doc}`06_resolution_validation`。

最細的層級是逐顆粒子估計軌跡。這時可以利用三維模型的投影當作參考，計算粒子影像與投影的互相關來得到每個 frame 的位移。單一 frame 的訊雜比太低，所以通常改用數個 frames 的滑動平均，並假設同一顆粒子的軌跡接近線性、鄰近粒子的軌跡相似 {cite}`rubinstein2016`。

這個層級也需要比較各 frame 的資訊量。先用某一 frame 的粒子資料分成兩半，各自重建並計算 half-map FSC；再對全部 frames 合成的粒子資料做同樣的獨立半組分析。比較這兩組訊號估計，才能從 Guinier plot 取得各 frame 的相對 B factor。這裡並非直接把單一 frame map 與包含它的 all-frame map 做 FSC。相對 B factor 為負時降低高頻權重，為正時相對提高；這套符號與一般 $\exp(-Bs^2/4)$ 的衰減 B 因子相反 {cite}`rubinstein2016`。

### Dose weighting 與 motion correction 的關係

頻率權重可以從資料估計，也可以根據累積劑量的經驗曲線設定。前者可同時反映運動與輻射損傷，後者主要描述曝光造成的衰減；兩者都與{doc}`05_image_formation`的 dose weighting 有關，採用的估計模型則有所不同。

時間軸的兩端因此有不同的問題。曝光最初 1–2 $\mathrm{e^-/\mathring A^2}$ 期間的 beam-induced motion 通常最劇烈，使最早的幾個 frames 難以使用；高解析度資訊主要來自其後的早期 frames，而低解析度資訊（正是挑粒子與定向所倚重的）可以從整部 movie 累積，因此常使用遠高於 $20\ \mathrm{e^-/\mathring A^2}$ 的總劑量 {cite}`sigworth2016`。因此設定各 frame 的權重時，需要同時考慮殘餘運動與輻射損傷。

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

若校正後的 micrograph 仍有方向一致的拖影，或不同區域的細節清晰度差很多，可以回頭檢查局部運動模型與平滑參數是否合適。

## CTF estimation：從 power spectrum 認出振盪

{doc}`05_image_formation` 已給出 CTF 的形式與參數。這一步要從資料把 $\Delta f_1$、$\Delta f_2$ 與散光方位角估出來。

估計的依據是 micrograph 的二維 power spectrum。把它寫成三項的組合：

$$
\mathrm{PSD}_I(\mathbf k)=H^2(\mathbf k)\,\mathrm{PSD}_f(\mathbf k)+\mathrm{PSD}_b(\mathbf k),
$$

其中 $\mathrm{PSD}_f$ 是粒子、碳膜與污染物的平均功率頻譜，$\mathrm{PSD}_b$ 是背景雜訊的功率頻譜 {cite}`singer2020`。

這個式子能否用來估計 CTF，取決於各項隨頻率變化的方式。等號右邊有兩個完全未知的函數與一組 CTF 參數，單看方程式本身是不定的。讓問題變得可解的，是三者在頻率上的行為差異：$\mathrm{PSD}_f$ 與 $\mathrm{PSD}_b$ 都隨頻率緩慢且大致單調地變化，而 $H^2$ 快速振盪 {cite}`singer2020`。因此可以先用局部平均或低次多項式估出緩慢變化的背景，扣掉之後剩下的快速振盪就幾乎只受 CTF 參數控制。這是一種靠尺度分離達成的可識別性（identifiability）：兩組未知量若在同一個頻率尺度上變化，就沒有辦法只從一條曲線分開。

扣掉背景後，可以用局部極小值的位置直接讀出零點，也可以用理論 $H^2$ 與扣背景、適當正規化後的 power spectrum 做比對，取最佳的參數組合。後者在雜訊較大時較穩定 {cite}`singer2020`。

```{figure} images/pptx/s08_3.png
:width: 80%
:name: fig-ctf-radial
徑向平均的 power spectrum（綠）、CTF 功率模型的擬合曲線（橘）與逐頻率的擬合品質（藍）。快速振盪的成分帶有 CTF 資訊；擬合品質在高頻下降，表示該頻段已缺乏可辨認的 Thon rings。
```

散光讓離焦隨方位角改變，Thon rings 於是變成橢圓。徑向平均把同一半徑上的方位角資訊合併，只保留接近平均離焦的徑向變化，方向資訊隨之消失。因此散光的大小與方向必須從二維頻譜擬合，做法與一維相同，只是把互相關推廣到二維 {cite}`singer2020`。

```{figure} images/pptx/s08_2.png
:width: 55%
:name: fig-ctf-fit
實驗 2D power spectrum 與 CTF 擬合的四象限對照；非圓對稱的 Thon rings 攜帶散光資訊。
```

估好 CTF 之後，有三種用途不同的處理常被混在一起：

**Phase flipping** 乘上 $\operatorname{sign}(H)$，只修正對比的正負號，被 CTF 壓低的振幅維持原狀。它的優點是保持雜訊的統計性質 {cite}`singer2020`。**Wiener 型校正**同時處理正負號與振幅，並用訊雜比資訊抑制 $H$ 接近零時的放大；它需要 SSNR 的估計。**多離焦合併**則是在 2D 類別平均或 3D 重建時，讓不同影像的零點互相錯開，補上單張影像缺少的頻率。前兩者作用在單張影像上，無法救回零點處遺失的資訊；只有第三種能真正填補頻率缺口。

較大的離焦能加強部分低頻對比，也使高頻的 CTF 振盪更密集，增加頻率取樣與尺度校準的要求，並加強 delocalization。離焦值本身需要多準，則由目標頻率與電子波長決定：$|\delta\chi|=\pi\lambda s^2|\delta f|$ 在固定 $s,\lambda$ 下不依賴原先的離焦量。提高目標解析度時，必須降低離焦估計誤差，才能控制相位誤差。

**前置概念**：{doc}`03_fourier` 的卷積定理、phase flipping 與 Wiener filter；成像模型見 {doc}`05_image_formation`。

## Particle picking 與選樣偏差

Particle picking 要在 micrograph 找出候選中心並避開污染。常見方法包括 blob／difference-of-Gaussians 偵測、以參考影像做模板匹配、依資料統計特徵建立的無監督方法，以及以標註資料訓練的監督式模型 {cite}`sigworth2016,singer2020`。

從統計角度看，這一步做的是取樣：它從 micrograph 的所有位置中選出一個子集合，而後面所有推論都只針對這個子集合。因此挑選規則裡的任何偏好，都會變成資料本身的性質。

以模板匹配為例，若用一張細節豐富的參考影像當模板，被選中的就是那些和模板相關最高的位置；即使 micrograph 裡完全沒有粒子，被選中的純雜訊小圖平均起來也會浮現模板的形狀。這個現象在文獻中被稱為 “Einstein from noise”：用愛因斯坦的照片當模板去純雜訊影像中挑選，平均後會得到近似愛因斯坦的圖案 {cite}`singer2020`。以資料訓練的深度模型同樣可能帶有這類偏差。

**閾值決定了兩種錯誤的比例**。閾值放寬會納入更多真粒子，也納入更多冰晶、碳膜邊緣、聚集物與純雜訊；收緊則相反，但可能系統性地漏掉對比較低的取向或較小的粒子。由於不同取向的粒子投影對比本來就不同，閾值會直接改變取向分布。

誤挑的粒子會繼續影響後續估計。它們會參與 2D 與 3D 的平均，稀釋共同訊號，也可能被分類演算法聚成看似合理的類別。{doc}`06_resolution_validation` 會說明，這類系統性的錯誤成分屬於偏差而非隨機誤差，增加資料量無法消除它。

實務上最快的初步檢查是把挑選座標疊回 micrograph 看一遍。座標大量落在碳膜或污染上，表示閾值或訓練資料需要調整；同一顆粒子附近出現多個座標，表示最小距離設定太小。接著檢查裁切框：粒子應大致位於中央，而且 box size 要足以容納 CTF 造成的 delocalization（{doc}`05_image_formation` 給出的估計）。

```{figure} images/pptx/s07_3.png
:width: 60%
:name: fig-picking-overlay
把候選座標疊回 micrograph。圈選的位置是否對應真實粒子、是否漏掉某些取向，需要逐張檢視才看得出來。
```

```{figure} images/pptx/s09_1.png
:width: 60%
:name: fig-picking
軟體標出的候選粒子位置仍須經裁切、分類與資料品質檢查。
```

## 交給統計模型之前

走到這裡，資料流已經把 movie 變成一疊粒子影像，並附上每張影像的 CTF 參數與來源 micrograph。剩下的未知量是每顆粒子的取向、平移與可能的構形類別，也就是 {doc}`05_statistical_inference` 裡的 $z_i$。

後面的驗證還需要保留前處理的來源資訊：每顆粒子屬於哪張 micrograph（影響 half-set 切分的獨立性）、屬於哪個 CTF group、來自哪一批載網。把這些中繼資料和影像一起保存，是 STAR 等格式存在的原因。

為了能追蹤各步驟與重現結果，資料與參數需要一起保存。

處理鏈的每一步都應該只做一件事，並把輸出寫成下一步可讀的中繼資料，這樣才能在不重跑整條流程的情況下更換其中一步。每次執行使用的參數、隨機種子與軟體版本要跟著結果一起記錄，否則重現性檢查會失去依據。互相獨立的檢查應該分開保存，例如取向分布、每類粒子數與局部解析度分別回答不同的問題。最後，任何把資料分成兩半的規則必須在高解析度資訊產生之前就決定；{doc}`06_resolution_validation` 會說明為什麼。

[Computational-CryoEM](https://github.com/phonchi/Computational-CryoEM) 依處理階段整理了論文、軟體與教學資源，可作為查找各步驟現行方法的入口。實際操作時，依正在使用的版本查閱 [RELION](https://relion.readthedocs.io/en/release-5.0/)、[cryoSPARC](https://guide.cryosparc.com/) 或 [Scipion](https://scipion-em.github.io/docs/) 的官方文件。

(refinement)=
(heterogeneity)=
## 從流程接到估計問題

處理鏈後半段的 2D 分類、初始模型、3D 分類與 refinement，直接估計角度、位移與結構。後續各章會分別說明這些方法，先看它們如何接在一起。

**2D 對位與分類**同時估計每張影像的平面內旋轉、平移與所屬類別，再用加權平均得到 class averages。對位與分類互相依賴：影像要先對齊才好分群，但要分群正確才知道該和誰對齊。處理這個耦合的方式，以及硬指派與機率加權的差別，見 {doc}`06_alignment_classification`。

**Refinement** 把同一套想法推到三維。演算法從目前模型產生候選投影，為每張影像評估各組候選角度與位移的合理程度，再用這些權重重建新模型，反覆迭代。最簡化的版本只保留相關最高的候選（hard projection matching），現代方法則對角度、位移與類別做邊際化，保留多個候選的不確定性。重建本身如何把加權後的投影合成三維密度，以及初始模型偏差、對稱性與手性的問題，見 {doc}`06_reconstruction_validation`。

**異質性**指的是同一批粒子並非全部來自同一個結構。離散的 3D classification 用有限個類別描述它，本書介紹的連續方法 3DVA 則用少數幾個線性變化模式描述。兩種描述各有適用範圍，也各有可能把取向差異或成像條件誤讀成構形差異的風險；完整討論見 {doc}`06_heterogeneity`。

```{figure} images/pptx/s13_3.png
:width: 35%
:name: fig-3d-map
70S ribosome 的 3D 密度圖示例。密度圖外觀、解析度估計與生物詮釋是三個不同層次的問題。
```

## 接下來

下一章從二維開始：如何把一疊未對齊、未分類的粒子影像整理成可以互相比較的群組，以及為什麼這件事和統計推論密不可分。請接著閱讀 {doc}`06_alignment_classification`。若想動手逐項改變取向、平移、CTF 與雜訊，可操作 {doc}`07_synthetic_data`。

## 延伸閱讀

- **Ripstein, Z. A. 與 Rubinstein, J. L.（2016）．*Processing of Cryo-EM Movie Data*. Methods in Enzymology, 579, 103–124。** [DOI](https://doi.org/10.1016/bs.mie.2016.04.009)。本章運動估計三個層級的主要來源：第 2 節（頁 107–109）是全域 frame 對齊與最小平方解，第 3 節（頁 109–112）與第 5 節（頁 114–116）是區塊與光流做法，第 4 節（頁 112–114）說明樣條平滑與曝光加權，第 6 節（頁 116–121）是逐粒子軌跡與相對 B factor 權重。全章使用一致的符號，適合逐節比較不同方法的假設 {cite}`rubinstein2016`。
- **Singer, A. 與 Sigworth, F. J.（2020）．*Computational Methods for Single-Particle Electron Cryomicroscopy*. Annual Review of Biomedical Data Science, 3, 163–190。** [DOI](https://doi.org/10.1146/annurev-biodatasci-021020-093826)。第 3 節（頁 171–174）依序處理 motion correction、CTF estimation 與 particle picking；本章使用的 $\mathrm{PSD}_I=H^2\mathrm{PSD}_f+\mathrm{PSD}_b$ 與可識別性論證出自第 3.2 節，“Einstein from noise” 的討論出自第 3.3 節 {cite}`singer2020`。
- **Sigworth, F. J.（2016）．*Principles of Cryo-EM Single-Particle Image Processing*. Microscopy, 65(1), 57–67。** [Oxford Academic 免費全文](https://doi.org/10.1093/jmicro/dfv370)。頁 62–63 的 Particle picking 與 CTF determination 兩節簡短說明 difference-of-Gaussians 模板、模板偏差與離焦精度的要求；頁 65 給出本章引用的傳統低劑量設定（約 $20\ \mathrm{e^-/\mathring A^2}$）、曝光最初 1–2 $\mathrm{e^-/\mathring A^2}$ 的劇烈運動，以及高低頻資訊各自來自曝光哪一段的說明。公開取用，適合作為入門對照 {cite}`sigworth2016`。
- **Glaeser, R. M.（2016）．*Specimen Behavior in the Electron Beam*. Methods in Enzymology, 579, 19–50。** [DOI](https://doi.org/10.1016/bs.mie.2016.04.010)。第 6 節（頁 37–40）說明冷凍樣品為何會出現鼓面般的彎曲與難以預測的集體位移，這是 patch 層級運動模型存在的物理理由；第 7 節（頁 40–42）討論可能的機制 {cite}`glaeser2016`。
- **Cheng, A., Tan, Y. Z., Dandey, V. P., Potter, C. S. 與 Carragher, B.（2016）．*Strategies for Automated CryoEM Data Collection Using Direct Detectors*. Methods in Enzymology, 579, 87–102。** [DOI](https://doi.org/10.1016/bs.mie.2016.04.008)。說明資料如何被自動收集、每個位置的曝光與對焦如何安排；理解 movie 與 micrograph 的來源後，前處理的許多設計選擇會變得比較清楚 {cite}`cheng2016`。
- **[Computational-CryoEM](https://github.com/phonchi/Computational-CryoEM) 資源彙整。** 開放維護的清單，依處理階段列出論文、開放原始碼工具與教學。想知道某一步目前有哪些做法時，從對應章節的條目往外找，比直接搜尋關鍵字有效率。
