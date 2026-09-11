# 二維對位與分類：把可比較的影像放在一起

{doc}`06_workflow` 交給我們一疊粒子影像。每張影像都含有同一種分子的某個投影，但它們的平面內旋轉、平移各不相同，投影方向也不同，其中還混著誤挑的雜物與可能的不同構形。在做任何平均之前，必須先決定兩件事：每張影像該怎麼轉正，以及哪些影像可以放在一起平均。

這兩個問題互相依賴。要把影像轉正，得先知道它們該和誰對齊；要判斷哪些影像相似，又得先把它們轉到同一個座標系。本章說明這個耦合如何處理，如何量化對位誤差造成的高頻訊號衰減，以及幾種實際使用的做法各自做了什麼假設。

(rigid-alignment)=
## 剛體對位：把變換寫成一個群

一張二維影像的平面內對位，可以用兩個參數描述：旋轉角 $\psi$ 與平移 $\mathbf t$。把這對參數寫成 $g=(\psi,\mathbf t)$，它對影像的作用是

$$
(g\cdot a)(\mathbf x)=a\!\left(R_\psi^{-1}(\mathbf x-\mathbf t)\right),
$$

其中 $R_\psi$ 是二維旋轉矩陣。這些變換構成一個群：兩個變換接連作用仍是同一形式的變換，每個變換都有反變換，恆等變換對應 $(0,\mathbf 0)$。在冷凍電鏡裡通常還要加上鏡射，因為粒子可能以兩個面之一朝向偵測器；加入鏡射後群的結構也隨之改變。

把變換視為群，有一個直接的後果：**對位的解只在一個整體變換之下確定**。若把所有影像的 $g_i$ 都右乘同一個 $g_0$，同時把參考影像左乘 $g_0^{-1}$，所有殘差完全不變。因此「粒子在類別平均中的絕對方位」沒有意義，需要另外指定一個約定（例如把平均影像的質心放在中央）。同樣地，比較兩次獨立分析的對位參數時，必須先把兩者對到同一個約定上。

```{figure} images/pptx/s24_1.png
:width: 38%
:name: fig-angles
三個角度描述一顆粒子：$\phi$ 與 $\theta$ 決定投影方向在球面上的位置，$\psi$ 是影像平面內的旋轉。二維對位處理的是 $\psi$ 與平移，三維 refinement 才同時估計 $\phi$、$\theta$。
```

## 最小平方為什麼變成相關

對位的目標函數通常寫成平方誤差。給定參考影像 $a$，找出使殘差最小的變換：

$$
\hat g_i=\arg\min_{g}\left\|y_i-g\cdot a\right\|^2 .
$$

把平方展開，

$$
\left\|y_i-g\cdot a\right\|^2
=\|y_i\|^2-2\left\langle y_i,\;g\cdot a\right\rangle+\|g\cdot a\|^2 .
$$

第一項與 $g$ 無關。若變換也保持參考影像的能量，使第三項與 $g$ 無關，那麼最小化平方誤差就等價於**最大化內積**，也就是互相關。此處的等價關係以變換保持參考影像能量為前提 {cite}`sorzano2010`。

這個條件什麼時候會失效？旋轉在連續、圓對稱支撐的情形下保持能量，但離散影像的內插會改變它；平移在週期延拓下保持能量，實際裁切時卻會把內容移出方框。套用 mask 之後，被 mask 保留的能量也隨 $g$ 改變。因此實作上常改用正規化互相關，或先把影像做圓形 mask 再比較。

這個目標函數也採用了特定的雜訊假設。 {doc}`05_statistical_inference` 說明，未加權的平方誤差對應「白雜訊、各像素變異數相同」的高斯概似。冷凍電鏡的雜訊在低頻遠強於高頻，所以直接用未加權平方距離比較，等於讓低頻主導對位分數。實務上會先做 whitening 或帶通濾波，把各頻率的貢獻調整到與其訊雜比相稱。

## 對位誤差如何變成高頻的衰減

對位通常會留下誤差。指定誤差模型後，就能計算它對平均訊號的影響。設所有影像的訊號相同，但第 $i$ 張影像被套上一個估計錯誤造成的殘餘平移 $\boldsymbol\delta_i$。在 Fourier 空間，平移只改變相位：

$$
\widehat{y_i}(\mathbf k)=\widehat s(\mathbf k)\,e^{-2\pi\mathrm i\,\mathbf k\cdot\boldsymbol\delta_i}+\widehat\varepsilon_i(\mathbf k).
$$

平均 $N$ 張影像時，訊號項變成 $\widehat s(\mathbf k)\cdot\frac1N\sum_i e^{-2\pi\mathrm i\mathbf k\cdot\boldsymbol\delta_i}$。$N$ 夠大時，括號內收斂到期望值。若 $\boldsymbol\delta_i$ 的兩個分量獨立、各為平均零、標準差 $\delta_0$ 的常態分布，這個期望值是常態分布的特徵函數：

$$
\mathbb E\left[e^{-2\pi\mathrm i\,\mathbf k\cdot\boldsymbol\delta}\right]
=\exp\!\left(-2\pi^2\delta_0^2\,\lvert\mathbf k\rvert^2\right).
$$

也就是說，**隨機的對位誤差在平均影像上留下一個隨頻率遞減的包絡**，形式和 CTF 的 envelope 或晶體學的溫度因子完全相同。把它寫成 $\exp(-Bs^2/4)$ 的形式，可得 $B=8\pi^2\delta_0^2$；$\delta_0=1\ \mathrm{\mathring A}$ 對應 $B\approx79\ \mathrm{\mathring A^2}$。

例如，取 $\delta_0=1$ pixel：

在 $\lvert\mathbf k\rvert=1/6$ pixel⁻¹（對應 6 pixels 的尺度），保留比例為 $\exp(-2\pi^2/36)\approx0.58$；在 $\lvert\mathbf k\rvert=1/3$ pixel⁻¹（3 pixels），只剩 $\exp(-2\pi^2/9)\approx0.11$。

一個像素的對位誤差幾乎不影響粗輪廓，卻讓三個像素尺度的細節只剩下約一成。低頻看起來清楚的類別平均影像，完全可能在高頻已被對位誤差抹平。

小的平面旋轉誤差會在離中心較遠的地方造成較大位移。半徑 $r$ 處，切向位移約為 $r\,\Delta\psi$。若角度誤差為零平均 Gaussian，標準差為 $1.7^\circ$，在 $r=65$ pixels 的外緣，切向位移標準差約為 $65\times0.0297\approx1.9$ pixels。這時沿外緣切向、6 pixels 尺度的細節，在小角度局部近似下保留約 0.13 的振幅；徑向細節的影響則不同。固定的角度偏轉只會旋轉影像，隨機角度誤差的平均才造成模糊。

```{dropdown} 小角度旋轉為何造成方向性的模糊？
令 $J=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$，位置 $\mathbf x$ 的一階位移為 $\Delta\psi J\mathbf x$。若 $\Delta\psi\sim N(0,\sigma_\psi^2)$，局部 Fourier 模式的平均相位因子為

$$
\exp\![-2\pi^2\sigma_\psi^2(\mathbf k\cdot J\mathbf x)^2].
$$

它依賴位置與頻率方向。外緣切向頻率受影響最大；當 $\mathbf k$ 沿徑向，內積為零，一階近似沒有此衰減。它不能當成整張影像的各向同性平移包絡。
```

角度精度本身又由訊雜比決定。以模擬資料測試 projection matching 時，當粒子影像的訊雜比與某組實驗資料相當，角度誤差的標準差約為 $1.7^\circ$；訊雜比減半時升到 $3^\circ$，再減半則常出現 $10^\circ$ 以上的誤差，低解析度重建都變得困難 {cite}`sigworth2016`。誤差對訊雜比的這種急遽反應，是 SPA 需要盡量壓低每一環節額外雜訊的理由。

提高角度精度也會增加計算成本。窮舉式的 projection matching 要為每張影像測試所有候選方向、平面內旋轉與平移。把影像邊長記為 $n$、平移不確定度記為 $t$ pixels、影像數記為 $m$、對稱階數記為 $n_{\mathrm{sym}}$，運算量大致是

$$
n_c\approx \frac{t^2\pi^2n^5m}{n_{\mathrm{sym}}},
$$

若粒子直徑為 $D$、目標解析度為 $d$，在相應取樣下約有 $n\simeq2D/d$。這個估算描述逐一比較候選的做法，顯示成本如何隨粒子大小與目標解析度上升 {cite}`sigworth2016`。這解釋了為什麼實作上會用多解析度搜尋、局部搜尋與各種加速結構，而這些近似本身又會改變估計的性質。

(alignment-classification-coupling)=
## 對位與分類為何需要一起估計

前面的推導假設參考影像已知，但實際分析往往還需要從資料建立參考。資料又含有多個投影方向，把方向差異很大的影像強行對齊平均，會使結構細節模糊。

因此對位與分類必須一起做。這個耦合的困難在文獻中被明確指出：同一群內的影像必須夠相似才能正確對齊，但要用分群技術找出相似的影像，又需要影像已經對齊 {cite}`yang2012`。

實作上的共同架構是迭代：先指定一組初始參考（或直接用隨機挑出的影像），把所有影像對到最相似的參考並分群，用各群的平均更新參考，重複。這個架構稱為 multireference alignment。它有效，但也帶來幾個可預期的問題。

初始參考會影響結果。 低訊雜比時，互相關會出現許多接近的局部極大；初始參考中的偶然特徵會在迭代中被放大，最後的分群強烈依賴初始種子 {cite}`scheres2005`。

類別大小也會影響指派。 比較影像與各類平均時，粒子多的類別其平均影像雜訊較低，容易在分數上佔優。可以把這件事算出來：設影像 $y=s+\varepsilon$，類別平均 $a_k=s_k+\eta_k$，其中 $\eta_k$ 的每像素變異數為 $\sigma^2/N_k$。此例假設參考由獨立影像建立，或比較時採 leave-one-out，使 $\varepsilon$ 與 $\eta_k$ 不相關。則

$$
\mathbb E\left\|y-a_k\right\|^2
=\|s-s_k\|^2+d\sigma^2+\frac{d\sigma^2}{N_k}.
$$

取 $d=16{,}900$ 個像素、$\sigma^2=1$。由 1,000 張影像組成的大類別 A，最後一項為 16.9；由 50 張影像組成的小類別 B，這一項是 338。

若測試影像真正屬於小類別 B，當 $\|s_B-s_A\|^2<338-16.9=321.1$ 時，大類別 A 的期望平方距離反而較低。這顯示參考雜訊水準可能影響指派；單次含雜訊的比較仍會波動，不能由期望直接斷定每張影像的歸屬。若小類別因此流失成員，平均影像的雜訊會進一步增加 {cite}`sorzano2010,yang2012`。

參考與待對位影像共用雜訊時，還會出現人為的一致性。 把影像對齊到由同一批影像算出的平均，等於讓每張影像的雜訊參與決定自己的對位參數。這會讓對齊後的資料看起來比實際更一致，文獻中稱為 noise alignment {cite}`penczek2010resolution`。{doc}`06_resolution_validation` 會說明這件事如何直接影響解析度估計。

(soft-alignment)=
## 硬指派與機率加權

低訊雜比下，多個候選可能得到相近分數。對位方法需要決定，要只保留最佳候選，還是同時考慮其他可能。

**硬指派**為每張影像選出分數最高的變換與類別，其餘候選丟棄。它容易實作、容易理解，但在最佳與次佳分數幾乎相同時，這個選擇幾乎是隨機的，而下一輪迭代會把這個隨機選擇當成事實。

**機率加權**改為對所有候選計算後驗權重。以 $z$ 表示類別、平面內旋轉與平移的組合，$\Theta$ 表示類別平均、雜訊等模型參數，

$$
r_i(z)=\frac{p(y_i\mid z,\Theta)\,p(z)}
{\sum_{z'}p(y_i\mid z',\Theta)\,p(z')},
\qquad \sum_z r_i(z)=1,
$$

再用這些權重更新類別平均。連續的旋轉與平移需要積分；以離散網格近似時，加權和應包含網格的體積元素，等權相加只適用於相應的離散先驗或等測度網格。完整推導、EM 的單調性，以及權重如何隨雜訊變異數變化，見 {doc}`05_statistical_inference`。

機率加權保留了其他角度、位移與類別的可能性。在模擬資料上，這個做法得到的模型偏差比傳統互相關對位小，在低訊雜比時差異更明顯 {cite}`scheres2005`。

## 四種實際做法的比較

下面比較聯合對齊／分類、穩定性檢查，以及對齊後的穩健分群，說明各方法處理流程中的哪一個問題。

**ML2D**（maximum-likelihood multireference refinement）把類別標籤 $k$、平面內旋轉與平移一起當成隱藏變數，對它們做邊際化，並把雜訊標準差與平移的先驗寬度也列為待估參數，用 EM 更新 {cite}`scheres2005`。它將雜訊大小列入待估參數，讓權重隨資料更新；代價是每次迭代都要對整個角度與位移網格求和，計算量大。

(cl2d)=
**CL2D** 改動的是相似度與分群策略。它用 correntropy 取代相關：

$$
V_\sigma(X,Y)=\mathbb E\left[\kappa_\sigma(X-Y)\right],
\qquad
\kappa_\sigma(x)=\exp\!\left(-\frac{x^2}{2\sigma^2}\right).
$$

這是一個相似度：兩張影像完全相同時 $V_\sigma$ 取到上界 1，差距愈大值愈小。把核函數的指數展開可以看出，$V_\sigma$ 包含 $X-Y$ 的所有偶次動差，而相關只用到二次動差 {cite}`sorzano2010`。額外的高次項讓最佳點附近的曲面更尖銳，因此對小的錯位與小的類別差異更敏感。核寬度 $\sigma$ 由背景區域估出的雜訊變異數決定。分群方面，CL2D 採用分裂式的向量量化：先把最大的類別依相似度高低對半切開，再重新分群，反覆把類別數加倍，使每個類別平均的影像數維持在較小的規模 {cite}`sorzano2010`。它另外設計了一個穩健的指派準則，處理上一節算過的「大類別佔優」問題。

(isac)=
**ISAC**（iterative stable alignment and clustering）把重點放在驗證。它用等大小分群（EQK-means）強制所有類別成員數相同，直接避免小類別塌陷 {cite}`yang2012`。它還使用兩層檢查：**對位穩定性**以隨機初始化重跑 $L$ 次無參考對位（原文取 $L=5$），比較同一張影像在各次得到的旋轉、平移與鏡射，用像素誤差判斷它是否穩定；**可重現性**則把整個分群流程在近乎獨立的條件下重跑多次，用多分割比對找出成員一致的類別，隨迭代把要求從兩次一致逐步提高到四次一致。通過檢查的影像與類別被取出，剩下的影像再跑一次，如此進行約十輪；原文的示例中，第一輪約取出一半影像，後續每輪取出的比例下降到一成以下 {cite}`yang2012`。ISAC 交出的因此是一個經過驗證的子集合，而非全部影像的分群結果。

(gamma-sup)=
**γ-SUP** 對已對齊的影像做穩健分群，處理前面對齊步驟留下的錯位離群值；它本身不估計旋轉與平移來修正影像。它用 q-Gaussian 混合描述資料，以最小 γ-divergence 估計，再用 self-updating 的迭代求數值解 {cite}`chen2014`。在原文的 $q<1$、$\gamma>1-q$ 條件下，這個組合有兩個性質：最小 γ-divergence 會對偏離的點自動降權（soft rejection），原文選用 $q<1$ 的 q-Gaussian 具有有限支撐，因此偏離超過一定範圍的點會被完全排除。由於每個資料點一開始都是自己的群代表，演算法既不需要事先指定群數，也不需要隨機初始中心。作者以刻意錯位的模擬影像測試，γ-SUP 能把大部分錯位影像標為離群值 {cite}`chen2014`。

這四種方法的著重點不同：ML2D 改的是**統計模型**，CL2D 改的是**相似度與分群結構**，ISAC 改的是**驗證程序**，γ-SUP 改的是**估計的穩健性**。前三者處理對齊與分類流程中的不同問題，γ-SUP 則在已有對齊結果上處理離群值與穩健分群，在不同資料上的表現取決於哪個弱點是主要瓶頸。

(class-stability)=
## 類別穩定嗎？類別是什麼？

清楚的類別平均顯示，分到同一類的影像在二維上互相一致。要判斷這種一致性是否穩定、是否來自相同構形，還需要其他證據。

**穩定性檢查**是最直接的。換一組初始種子、換隨機數種子、把資料重新切分後重跑，比較成員是否大致相同。ISAC 把這件事做成流程的一部分，但任何方法都可以用同樣的方式自行檢查 {cite}`yang2012`。若某個類別在重跑後散掉，它多半反映的是某次迭代的偶然結果。

**離群值與污染**會以兩種方式出現。少量嚴重偏離的影像（冰晶、聚集物、誤挑的雜訊）可能自成一類，也可能散進各類別拉偏平均。穩健方法把它們標為離群值，一般方法則需要人工檢視類別平均並剔除 {cite}`chen2014`。

**剔除本身是一種取樣。** 每次丟掉「看起來不好」的類別，都在改變剩下粒子的母體。若被丟掉的類別系統性地對應某些取向或某些構形，取向分布與構形比例都會跟著改變。檢查 2D 類別時把每類的粒子數一併記錄，後續才能追蹤資料量與分布怎麼變化。

影像相似可能有多種原因。 同一個類別裡的影像之所以相似，可能因為投影方向接近，也可能因為離焦相近、冰厚相近、位於同一批載網，或只是含有相似的對位誤差。把一組類別直接讀成一組構形，需要先排除這些成像因素。{doc}`06_heterogeneity` 會用三維的觀測模型處理這個問題，{doc}`06_dimension_reduction` 則說明如何把主要變化方向和已知的成像參數對照。

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

```{figure} images/pptx/s12_6.png
:width: 28%
:name: fig-class-masked
套上圓形 mask 的類別平均。Mask 讓比較集中在粒子所在區域，同時也改變了「能量隨變換保持不變」這個假設。
```

模糊的類別平均通常無法單獨指出成因。可能是污染或誤挑、平面內對位錯誤、該類粒子數太少、同一類混入不同構形，或偏好取向造成的分布不均。可依序查看原始粒子框與挑選座標、旋轉與位移的分布、每類粒子數、不同初始值下的分類穩定性，以及後續 3D 的方向分布，分別判斷資料、對位與分類可能出了什麼問題。

{doc}`04_wavelet` 的多尺度觀察可幫助理解不同頻帶的結構。類別平均影像的訊雜比提升來自對位後的統計平均；小波閾值去雜訊採用另一套假設與運算。

## 接下來

本章把影像整理成可比較的群組，但每張影像仍然是一個上萬維的向量，而我們實際關心的變化可能只有少數幾個方向。下一章討論如何從整批資料學出共用的低維表示，以及這個表示保留了什麼、遺失了什麼：{doc}`06_dimension_reduction`。

## 延伸閱讀

- **Scheres, S. H. W., Valle, M., Nuñez, R., Sorzano, C. O. S., Marabini, R., Herman, G. T. 與 Carazo, J.-M.（2005）．*Maximum-likelihood Multi-reference Refinement for Electron Microscopy Images*. Journal of Molecular Biology, 348, 139–149。** [DOI](https://doi.org/10.1016/j.jmb.2005.02.031)。Introduction（頁 140）列出互相關式多參考對位的三個弱點：局部極大、對初始種子的依賴，以及初始偏好被放大；Mathematical Background（頁 140–142）給出把類別、旋轉與平移一起邊際化的模型與 EM 更新式，並說明 PSF 造成的像素相關如何違反獨立雜訊假設 {cite}`scheres2005`。
- **Sorzano, C. O. S., Bilbao-Castro, J. R., Shkolnisky, Y. 等（2010）．*A Clustering Approach to Multireference Alignment of Single-Particle Projections in Electron Microscopy*. Journal of Structural Biology, 171, 197–206。** [DOI](https://doi.org/10.1016/j.jsb.2010.03.011)。第 2.1 節（頁 198–199）定義 correntropy 並說明它為何比相關對小錯位敏感；第 2.2 與 2.3 節（頁 199–200）介紹分裂式向量量化與處理不同雜訊水準類別的穩健準則。閱讀時把它的相似度與本章的最小平方展開並排比較 {cite}`sorzano2010`。
- **Yang, Z., Fang, J., Chittuluru, J., Asturias, F. J. 與 Penczek, P. A.（2012）．*Iterative Stable Alignment and Clustering of 2D Transmission Electron Microscope Images*. Structure, 20(2), 237–247。** [DOI](https://doi.org/10.1016/j.str.2011.12.007)。Introduction 與 Results 前段（頁 237–240）說明對位與分群互相依賴的困境、K-means 的類別塌陷，以及 EQK-means 的等大小約束；接著的段落定義像素誤差的穩定性測試（$L=5$ 次無參考對位）與跨次執行的可重現性測試，並描述以「pass」逐步取出穩定子集合的流程 {cite}`yang2012`。
- **Chen, T.-L., Hsieh, D.-N., Hung, H., Tu, I-P., Wu, P.-S., Wu, Y.-M., Chang, W.-H. 與 Huang, S.-Y.（2014）．*γ-SUP: A Clustering Algorithm for Cryo-Electron Microscopy Images of Asymmetric Particles*. The Annals of Applied Statistics, 8(1), 259–285。** [開放預印本 arXiv:1205.2034](https://arxiv.org/abs/1205.2034)。第 2 節回顧 γ-divergence 與 q-Gaussian 分布，第 3 節建立 γ-SUP 並說明 soft rejection 與拒絕區如何共同帶來穩健性，第 4.2 節示範以刻意錯位的模擬影像檢驗離群值偵測。這篇原始方法論文可公開取得全文，適合細讀估計式的推導 {cite}`chen2014`。
- **Sigworth, F. J.（2016）．*Principles of Cryo-EM Single-Particle Image Processing*. Microscopy, 65(1), 57–67。** [Oxford Academic 免費全文](https://doi.org/10.1093/jmicro/dfv370)。頁 59–60 推導 projection matching 的運算量公式，頁 60–61 的 Figure 2 給出本章引用的角度誤差對訊雜比的關係（$1.7^\circ$、$3^\circ$、$10^\circ$ 三個水準），並以模擬影像顯示相近投影之間的差異有多細微 {cite}`sigworth2016`。
- **Sigworth, F. J., Doerschuk, P. C., Carazo, J.-M. 與 Scheres, S. H. W.（2010）．*An Introduction to Maximum-Likelihood Methods in Cryo-EM*. Methods in Enzymology, 482, 263–294。** [DOI](https://doi.org/10.1016/S0076-6879(10)82011-7)。從二維對位的概似開始，逐步推到分類與三維情形，可作為本章機率加權一節的完整推導 {cite}`scheres2010`。
