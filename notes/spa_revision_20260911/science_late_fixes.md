# 三章科學修正紀錄

## 1. `book/06_reconstruction_validation.md` 第 5 行
原句：

第一層是幾何：即使取向完全已知，從有限個投影恢復三維物體仍是一個病態的反問題，需要說清楚什麼被資料決定、什麼由正則化決定。

原因：有限投影或長方矩陣不必病態；改用識別與穩定性條件。

改後：

第一層是幾何：即使角度完全已知，投影方向、CTF 與離散取樣仍可能使部分密度無法識別，或讓估計對雜訊很敏感。要分清楚哪些部分能由資料決定，哪些部分需要正則化。

## 2. `book/06_reconstruction_validation.md` 第 57 行
原句：

## 離散化之後，反演本身就是病態的

上面的論述用的是連續函數。實際資料是有限張、取樣過的投影，而待求的是有限個 voxel。把三維密度排成向量 $V$、所有投影排成向量 $y$、投影運算排成矩陣 $A$，問題變成：已知 $y$ 與 $A$，求 $\widehat V$ 使 $A\widehat V\approx y$。

$A$ 幾乎總是長方形的，因此反演必然是病態的：依投影數目而定，系統可能超定或欠定 {cite}`penczek2010`。更嚴重的是，離散化帶來一個無法靠增加影像解決的問題——對任何一組有限的投影方向，都存在非零的物體 $\widehat V_0\neq0$，使得它在所有量測方向上的投影恰為零，也就是 $A\widehat V_0=0$ {cite}`penczek2010`。這類物體稱為 ghost。它們的存在表示重建結果可以和原物體差很多，卻仍與資料完全相符；離散的投影變換即使沒有雜訊也不可逆。

實務上 ghost 多半是很高頻的物體，因此適當的正則化可以壓住它們的影響 {cite}`penczek2010`。但這句話要正確讀：正則化把答案從「無窮多個同樣符合資料的解」縮到一個，靠的是資料以外的假設。{doc}`05_statistical_inference` 的 MAP 觀點把這個假設寫成先驗；哪一部分來自資料、哪一部分來自先驗，是後面驗證章節反覆要問的問題。



原因：改正長方矩陣必病態、有限voxel必ghost及增加影像無用的錯誤推論。

改後：

## 離散反問題：什麼時候能解，什麼時候穩定？

實際資料是有限張、取樣過的投影，待求的是有限個 voxel。把三維密度排成向量 $V$、所有投影排成向量 $y$、成像運算排成矩陣 $A$，問題變成：已知 $y$ 與 $A$，求 $\widehat V$ 使 $A\widehat V\approx y$。

矩陣是否長方形，不能決定反問題是否病態。例如

$$
A=\begin{pmatrix}1&0\\0&1\\1&1\end{pmatrix}
$$

雖然是長方矩陣，卻有兩個非零奇異值 $\sqrt3$ 與 $1$，最小平方解唯一，而且對小擾動相當穩定。真正要看的是矩陣的秩與奇異值：

- 若 $A$ 欠秩，就存在非零 $V_0$ 滿足 $AV_0=0$。因此 $V$ 與 $V+V_0$ 產生完全相同的量測，資料無法區分它們。
- 若 $A$ 滿秩，但某些奇異值很小，對應的結構變化只會在影像裡留下微弱痕跡。反演這些方向會放大雜訊，估計因而不穩定。

投影文獻中的 **ghost** 指在所選量測方向下投影為零的物體，也就是觀測算子的 null space 成分 {cite}`penczek2010`。有限方向的連續物體問題可以存在這類歧義；限定到特定的有限 voxel 模型後，是否仍有非零 null space，必須看實際矩陣，不能只由「投影有限」斷言一定有 ghost。

正則化會抑制不穩定方向，並在不可識別的部分選擇符合額外假設的解。增加互補方向或不同 CTF 的影像，有機會增加秩或改善奇異值；只重複相同量測，則不會補足 null space。{doc}`05_statistical_inference` 的 MAP 觀點把額外假設寫成先驗，讓我們能追問每個結構細節是由資料還是先驗決定。



## 3. `book/06_reconstruction_validation.md` 第 88 行
原句：

其中 $W$ 是量測權重。右邊正是（加權的）反投影；左邊的 $A^{\mathsf T}WA$ 描述每個 Fourier 位置被取樣了多少次。反投影只完成了分子；要得到估計值，還要除以取樣密度。

原因：一般normal operator有非對角插值耦合，不能當成逐voxel取樣次數。

改後：

其中 $W$ 是量測權重。右邊是加權反投影；左邊的 $A^{\mathsf T}WA$ 則描述成像與反投影如何共同作用在密度上。一般的內插會耦合不同 voxel，因此這是一個需要求解的矩陣方程。只有採用對角或 binning 近似時，才可將更新理解為逐個 Fourier 位置除以累積權重。

## 4. `book/06_reconstruction_validation.md` 第 97 行
原句：

考慮最簡單的情形：某個 Fourier 位置有三個觀測

原因：明定實數例子，使無共軛CTF公式的適用範圍清楚。

改後：

先用實數係數示範：某個 Fourier 位置有三個觀測

## 5. `book/06_reconstruction_validation.md` 第 105 行
原句：

估計是無偏的：$\widehat V=1$。

原因：無偏是重複取樣期望，不是每次估計恰等於真值。

改後：

估計是無偏的：$\mathbb E[\widehat V]=1$，但每次含雜訊的估計值仍會波動。

## 6. `book/06_reconstruction_validation.md` 第 107 行
原句：

在無雜訊的資料下 $\widehat V=1.21/2.21\approx0.548$，被朝零收縮了。

原因：將mean、variance及MSE全部放在同一重複取樣模型。

改後：

重複收集含雜訊的資料時，估計的期望為 $\mathbb E[\widehat V]=1.21/2.21\approx0.548$，比真值更接近零。

## 7. `book/06_reconstruction_validation.md` 第 129 行
原句：

把這個權重和上一節的 CTF 加權合併接起來，就得到 RELION 所採用的 MAP 更新的結構。對每個 Fourier 位置 $\mathbf k$：

$$
\widehat V^{(n+1)}(\mathbf k)=
\frac{\displaystyle\sum_i\sum_z r_i^{(n)}(z)\,
\frac{H_i(\mathbf k)\,\widehat y_i(\mathbf k)}{\sigma_i^2(\mathbf k)}}
{\displaystyle\sum_i\sum_z r_i^{(n)}(z)\,
\frac{H_i^2(\mathbf k)}{\sigma_i^2(\mathbf k)}
+\frac{1}{\tau^2(\mathbf k)}}
$$

（式中的 $\widehat y_i$ 與 $H_i$ 依姿態 $z$ 對應到該 voxel 的切片位置）。分子是姿態權重與 CTF 加權後的反投影，分母是累積的資料精度加上先驗精度 {cite}`scheres2012bayes`。

這個式子的三個特點值得記住。第一，它同時完成 CTF 校正、取樣密度補償與低通濾波，這些在早期做法裡是分開的、各帶參數的步驟。第二，$\sigma_i^2(\mathbf k)$ 與 $\tau^2(\mathbf k)$ 都隨頻率變化，而且由 EM 從資料本身更新，所以「該相信資料到多高的頻率」是估計出來的，而非人為設定 {cite}`scheres2012bayes`。第三，這個特性也帶來風險：若某一輪不小心把雜訊當成訊號，先驗會在下一輪允許更多雜訊。{doc}`06_resolution_validation` 會說明為什麼這使得獨立的資料分割成為必要。



原因：以精確operator式取代缺少z-dependent投影/shift的分式；補足超參數與2012經驗設定界線。

改後：

把候選權重加入成像模型，可以寫出精確的 MAP 正規方程式。令 $A_i(z)$ 包含候選角度的投影、CTF，以及位移對應的 Fourier 相位；$\Sigma_i$ 是雜訊共變異數，$\Lambda$ 是零平均 Gaussian 結構先驗的精度矩陣。固定這一輪的權重及精度後，更新滿足

$$
\left[\sum_{i,z}r_i^{(n)}(z)A_i(z)^*\Sigma_i^{-1}A_i(z)+\Lambda\right]V^{(n+1)}
=\sum_{i,z}r_i^{(n)}(z)A_i(z)^*\Sigma_i^{-1}y_i .
$$

星號表示共軛轉置；在實數座標下就是轉置。右邊把每張影像按候選角度、位移、CTF 與雜訊精度反投影；左邊同時納入成像算子的耦合與先驗。若模型含多個構形，$V$ 可堆疊各構形的密度，$A_i(z)$ 依候選類別選取對應部分。

上一節的單係數分式，是這個方程在各係數可分開估計時的簡化。一般的 gridding 與內插會耦合鄰近係數；使用逐 voxel 的「資料分子除以資料精度加先驗精度」時，必須說明採用的對角或 binning 近似，也不能漏掉候選角度的切片位置與位移相位 {cite}`scheres2012bayes`。

Scheres 2012 的 Bayesian 實作把 CTF 校正、取樣權重與先驗收縮放進同一估計過程，並從資料逐輪更新頻率相依的雜訊與訊號功率。不過，模型仍需先驗形式、角度取樣與參數設定；原文也使用經驗倍率 $T=4$ 調整先驗功率，不能把它理解為完全不需人工設定的方法 {cite}`scheres2012bayes`。

從重建本身更新先驗也有風險：若某輪把雜訊當成訊號，下一輪就可能保留更多雜訊。{doc}`06_resolution_validation` 會說明獨立 half-set 如何協助檢查這個問題。



## 8. `book/06_reconstruction_validation.md` 第 146 行
原句：

上述迭代是局部最佳化，它會收斂到起點附近的解。因此初始模型會影響結果——這件事有一個專門的名字：model bias，指最終結構受到初始模型的影響

原因：局部最佳化不保證參數解位於起點附近。

改後：

上述迭代屬於局部最佳化；不同初始模型可能導向不同的局部解，不能保證找到全域最佳解，也不保證解在起點附近。初始參考中的特徵可能因此保留在最終結構裡，形成 model bias

## 9. `book/06_reconstruction_validation.md` 第 154 行
原句：

每張投影同時約束 $|G|$ 個等價的切片方向，等效的資料量因此放大 $|G|$ 倍，Fourier 空間的覆蓋也更完整。

原因：對稱複製不產生獨立資料。

改後：

每張投影可約束 $|G|$ 個對稱等價的切片方向，並減少待估的自由度。這會改善對稱模型下的覆蓋與估計精度，但不等於增加 $|G|$ 倍獨立觀測：由同一張影像複製出的對稱資訊共享雜訊。

## 10. `book/06_reconstruction_validation.md` 第 181 行
原句：

RELION 與 cryoSPARC 是目前最常被比較的兩套系統。把它們的差別放在正確的層次上，比記住誰比較快有用。

原因：限定歷史實作，避免將2012/2017當現行能力。

改後：

以下比較 Scheres 2012 與 Punjani 等人 2017 兩篇論文描述的實作，分別看統計模型與最佳化策略。這些歷史版本不能代表兩套軟體目前的全部功能。

## 11. `book/06_reconstruction_validation.md` 第 148 行
原句：

cryoSPARC 採用後者：

原因：限定ab initio實作版本。

改後：

Punjani 等人 2017 的 cryoSPARC 方法採用後者：

## 12. `book/06_reconstruction_validation.md` 第 148 行
原句：

用這些帶雜訊的步伐廣泛探索結構空間，因此可以從隨機初始化開始

原因：SGD隨機梯度不保證廣泛探索或跳脫局部解。

改後：

用這些近似梯度更新結構，演算法可從隨機初始化開始

## 13. `book/06_reconstruction_validation.md` 第 185 行
原句：

這些帶雜訊的步伐讓演算法能從隨機初始化廣泛探索，因而支援 ab initio。

原因：不將SGD噪聲直接推為廣泛搜尋保證。

改後：

該文以隨機初始化配合這些更新做 ab initio 重建，但仍需檢查初始化敏感性。

## 14. `book/06_reconstruction_validation.md` 第 213 行
原句：

Mask 可以排除大面積溶劑、提高數值穩定性，也可能因為兩張 half-map 套上同一個帶有結構細節的 mask 而人為提高相關性。

原因：區分固定mask作用與共同訊號頻率洩漏/data-dependent bias。

改後：

Mask 可以排除大面積溶劑、提高數值穩定性，但也會混合頻率；共同的低頻訊號可能混入高頻，而由完整資料細節製作的 mask 還可能引入共同偏差。事先固定的同一個 mask 並不會自行讓原本獨立的兩組雜訊產生交叉相關。

## 15. `book/06_reconstruction_validation.md` 第 225 行
原句：

branch-and-bound 對位的細節在 Supplementary Note 1。

原因：已讀Note1為SGD，不能誤標branch-and-bound定位。

改後：

Supplementary Note 1 詳述 SGD，branch-and-bound 對位見主文相應方法說明。

## 16. `book/06_heterogeneity.md` 第 11 行

原句：

把所有粒子放在一起重建，得到的是兩者的加權平均。

原因：β平均需要已知角度及一致覆蓋權重，非任意joint reconstruction。

改後：

先假設角度已知，兩個狀態的方向覆蓋及重建權重相同，並使用線性重建。在這個簡化下，把所有粒子放在一起會得到兩個密度的加權平均；若角度與分類也要聯合估計，結果不一定等於這個平均。

## 17. `book/06_heterogeneity.md` 第 21 行

原句：

也就是說，一個只有 5 Å 的構形差異，就足以讓活動區域在 10 Å 以上的細節從平均圖中消失，而分子的其餘部分依然清晰。

原因：兩狀態cos絕對值振盪，不是單調高頻抹除。

改後：

但在 **5 Å（$k=0.2$）因子又回到 $|\cos\pi|=1$**。這兩個固定狀態造成的是隨頻率振盪的干涉：某些頻率抵消，其他頻率恢復，不能說所有高於 10 Å 所對應頻率的細節都消失。分子的固定部分則不受這個平移混合影響。

## 18. `book/06_heterogeneity.md` 第 25 行

原句：

若位移近似服從標準差 $\sigma_d$ 的常態分布

原因：向量characteristic function需要isotropic covariance假設。

改後：

若三維位移服從零平均、各方向等變異的常態分布 $\mathbf d\sim N(0,\sigma_d^2I)$

## 19. `book/06_heterogeneity.md` 第 32 行

原句：

和 {doc}`06_alignment_classification` 的對位誤差包絡形式相同。$\sigma_d=2\ \mathrm{\mathring A}$ 時，10 Å 尺度保留 0.45，5 Å 尺度只剩 0.04。

原因：用既有新圖明確區分離散干涉與連續衰減。

改後：

和 {doc}`06_alignment_classification` 的對位誤差包絡形式相同。$\sigma_d=2\ \mathrm{\mathring A}$ 時，10 Å 尺度保留 0.45，5 Å 尺度只剩 0.04。

```{figure} images/spa/heterogeneity_attenuation.png
:name: fig-heterogeneity-attenuation
:alt: 兩個固定平移狀態的振幅因子隨頻率振盪，Gaussian 連續位移的振幅包絡則隨頻率下降。
:width: 90%

兩種混合模型造成不同的頻率效應。兩個固定狀態可在某些頻率完全抵消後又恢復；Gaussian 位移分布則形成隨頻率下降的包絡。判讀局部模糊時，要先分清楚使用的是哪一種模型。
```

## 20. `book/06_heterogeneity.md` 第 77 行

原句：

**實務上通常從每類的粒子數反推。** 在 RELION 的使用經驗裡

原因：標明2016經驗值不是當前硬門檻。

改後：

**每類粒子數可協助選擇起始設定。** Scheres 2016 整理當時 RELION 的使用經驗

## 21. `book/06_heterogeneity.md` 第 77 行

原句：

這些數字是經驗法則，用意是讓每個類別的訊號足以支撐有意義的密度，而非某種理論上的最佳值。

原因：避免歷史類別數及粒子數標準化。

改後：

這些數字是 2016 的經驗法則，用意是避免每類訊號太少；實際所需粒子數仍隨分子大小、訊雜比、方向覆蓋、目標解析度與方法而變，不能當作目前軟體的硬門檻。

## 22. `book/06_heterogeneity.md` 第 79 行

原句：

一般 2D classification 使用 $T=2$，3D classification 使用約 $T=4$

原因：歷史T值不能無條件變成一般預設。

改後：

該 2016 章節的經驗設定是 2D classification 使用 $T=2$，3D classification 使用約 $T=4$

## 23. `book/06_heterogeneity.md` 第 81 行

原句：

**可識別性本身就受訊雜比限制。**

原因：Gaussian mixtures在低SNR仍可能理論可識別；有限資料分辨困難不等於不識別。

改後：

**有限資料下，能否分辨狀態也受訊雜比限制。**

## 24. `book/06_heterogeneity.md` 第 81 行

原句：

分類能否成功，因此先受限於資料本身的訊雜比，而非分類方法的巧妙程度。

原因：區分有限样本統計分辨與母體識別。

改後：

在有限樣本下，較低訊雜比會讓分類更難，也可能需要更多粒子；不能僅因直方圖看不出雙峰，就判定兩個族群在數學上永遠不可識別。

## 25. `book/06_heterogeneity.md` 第 83 行

原句：

那就是 4,000 顆——依上面的經驗法則，勉強夠重建一個低解析度的密度。

原因：不能由4000個粒子推定低解析度必足夠。

改後：

那就是 4,000 顆。僅憑這個數量，還無法判定能否得到有用的低解析度密度；也要看這些粒子的訊雜比與方向覆蓋。

## 26. `book/06_heterogeneity.md` 第 94 行

原句：

3D variability analysis（3DVA）把上面的線性模型直接放進成像模型。在 Fourier 空間，第 $i$ 張影像的預測是

$$
\widehat y_i=\alpha_i\,H_i\,P(\phi_i)
\left[V_0+\sum_{l=1}^{L}z_{il}V_l\right]+\eta_i ,
$$

其中 $P(\phi_i)$ 是姿態 $\phi_i$ 的投影算子、$H_i$ 是該影像的 CTF、$\alpha_i$ 是逐粒子的強度比例 {cite}`punjani2021`。這個寫法的關鍵是：**線性變化發生在三維，比較發生在二維**。三維的變化模式先經過投影與 CTF，才和觀測影像比較。因此 CTF 的振盪與各粒子不同的視角都被正確計入，而非事後修正。

演算法把潛在座標 $\mathbf z_i$ 視為零平均、單位共變異數的高斯隱藏變數，成為機率主成分分析（probabilistic PCA）的一個變體，再用 EM 求解：E-step 更新每張影像潛在座標的後驗平均，M-step 更新比例參數 $\alpha_i$ 與變化成分 $V_{1:L}$ {cite}`punjani2021`。雜訊被建模成在 Fourier 空間各環帶內變異數固定的 colored noise，並在計算前把影像與 CTF 一起做 whitening，把問題化簡成等向雜訊的最小平方 {cite}`punjani2021`——這與 {doc}`05_image_formation` 的 whitening 是同一個操作。

這個方法的能力來自幾個明確的假設，理解它們比記住演算法名稱更重要。

**它假設 $V_0$、每顆粒子的姿態 $\phi_i$ 與 CTF 都已知。** 實際流程是先用全部粒子做一次 consensus refinement，再把得到的姿態拿來做 3DVA {cite}`punjani2021`。於是 consensus refinement 的所有偏差——包括上一節提到的、活動區域在平均圖中被抹平——都會傳遞進來。若某個構形的粒子在 consensus 階段就被指派了錯誤的取向，它在潛在空間的位置也會跟著錯。



原因：修正Fourier投影算子定义、PPCA零噪聲計算極限、固定參數、全部粒子貢獻與矩陣可逆條件。

改後：

3D variability analysis（3DVA）把線性變化放進成像模型。沿用全書的 $A_i$ 記號，在 Fourier 空間寫成

$$
\widehat y_i=\alpha_i A_i
\left[V_0+\sum_{l=1}^{L}z_{il}V_l\right]+\eta_i .
$$

$A_i$ 依已估計的角度從三維 Fourier 密度擷取中心切片，再乘上 CTF 與位移相位；$\alpha_i$ 是逐粒子的強度比例 {cite}`punjani2021`。線性變化發生在三維，預測與觀測的比較則發生在二維。這樣才能同時考慮 CTF 的振盪及各粒子的不同視角。

原文先以 probabilistic PCA 描述零平均、單位共變異數的 Gaussian 潛在座標，並以 Fourier 環帶內變異數固定的 colored noise 描述雜訊。計算時先將觀測影像與成像算子一起 whitening；接著採用 PPCA 的零雜訊**計算極限**，得到交替最小平方更新。這個極限是簡化求解方式，不是在宣稱實驗影像沒有雜訊，也不是每轮保留完整潛在後驗共變異數的一般 PPCA EM {cite}`punjani2021`。

令 $\widehat y_i'$ 與 $A_i'$ 表示 whitening 後的觀測及算子，目標為

$$
\min_{V_{1:L},\{\alpha_i,\mathbf z_i\}}
\frac12\sum_i\left\|\widehat y_i'-\alpha_i A_i'
\left(V_0+\sum_{l=1}^{L}z_{il}V_l\right)\right\|^2 .
$$

解每張影像的 $\mathbf z_i$ 時，先固定變化成分、比例與成像算子；接著固定座標，更新比例及變化成分。因此所有粒子都參與各成分的估計，只是各自的係數不同。

```{dropdown} 固定變化成分後，如何求一張影像的座標？
令 $W_i=[A_i'V_1,\ldots,A_i'V_L]$，並記 $d_i=\widehat y_i'-\alpha_i A_i'V_0$。固定非零 $\alpha_i$ 時，要解的是

$$
\min_{\mathbf z_i}\|d_i-\alpha_iW_i\mathbf z_i\|^2.
$$

若 $W_i^*W_i$ 可逆，解為

$$
\mathbf z_i=\frac1{\alpha_i}(W_i^*W_i)^{-1}W_i^*d_i.
$$

若某個視角讓不同三維模式投影成相同或接近相同的二維變化，這個矩陣可能欠秩或病態，需使用廣義逆或適當正則化。單張影像因此不一定能可靠地區分所有模式。
```

**$V_0$、角度、位移與 CTF 都先固定。** 實際流程先做 consensus refinement，再把估得的平均密度與粒子參數交給 3DVA {cite}`punjani2021`。若某些粒子在 consensus 階段被指派錯誤角度，這些誤差可能進入潛在座標與變化模式。把 consensus map 當作已知，也表示這一階段沒有重新估計它的不確定性。



## 27. `book/06_heterogeneity.md` 第 135 行

原句：

各成分彼此正交，每個解釋一種不同的變化模式；某個成分的潛在座標在整批影像上的變異數愈大，表示它解釋的變異愈多

原因：mode尺度任意；需正交單位長度與主方向旋轉才能比较score variance。

改後：

這個比較需要先把成分正交化並正規化成單位長度，再依估得的座標共變異數旋轉到主方向。完成這些步驟後，某個座標在整批影像上的經驗變異數愈大，才表示相應成分解釋的變異愈多

## 28. `book/06_heterogeneity.md` 第 135 行

原句：

這是一個關於資料中觀測變異的陳述，把它讀成物理上的重要性需要額外的論證。

原因：區分empirical mode-score variance與N0I prior。

改後：

這裡比較的是估計座標的經驗變異數，不是前面單位共變異數先驗的固定變異數；也不能直接把它當成生物功能上的重要性。

## 29. `book/06_heterogeneity.md` 第 154 行

原句：

那麼「該構形」與「該視角」在資料中就是混淆的，任何只看影像的方法都無法把兩者分開。這時需要改變製樣條件並看結論是否穩定。

原因：視角與構形相關不等於永久不可識別。

改後：

兩者在資料中就會相關。相關不代表一定無法區分；若各狀態仍有足夠且互補的角度覆蓋，成像模型可能提供區分線索。若缺乏跨狀態、跨角度的覆蓋，資料才可能不足以分開這兩種解釋。改變製樣條件或補充角度，有助於檢查結論。

## 30. `book/06_heterogeneity.md` 第 168 行

原句：

**Masked refinement** 在每一輪對參考套上三維 mask，讓對位只依據 mask 內的質量進行，於是 mask 外的變化不影響取向估計

原因：參考mask不能刪除觀測中的重疊外部投影。

改後：

**Masked refinement** 在每一輪對參考套上三維 mask，讓匹配更著重 mask 內的質量。不過，觀測影像仍含有外部結構的投影，外部訊號可能與局部訊號重疊，所以只遮參考並不能完全消除外部變化對角度估計的影響

## 31. `book/06_heterogeneity.md` 第 170 行

原句：

這個做法有一個規模下限。要能對某個區域獨立估計取向，該區域本身必須含有足夠的質量——經驗上大約需要和一個可以獨立解結構的複合體相當，也就是至少 150–200 kDa

原因：150–200kDa是2016經驗，不是現行物理門檻。

改後：

要對某個區域獨立估計角度，該區域本身必須含有足夠的訊號。Scheres 2016 的章節提出約 150–200 kDa 的經驗尺度，反映當時資料與方法的條件，並非現在通用的質量下限

## 32. `book/06_heterogeneity.md` 第 172 行

原句：

Mask 一律需要柔邊，否則會在 Fourier 空間造成振鈴；

原因：新增局部訊號扣除及model/angles條件，避免mask萬能敘述。

改後：

若已有外部區域的可靠模型 $V_{\mathrm{out}}$，可進一步做 partial signal subtraction：

$$
y_i^{\mathrm{sub}}=y_i-A_i\widehat V_{\mathrm{out}}.
$$

用同一組角度、位移、CTF 與強度慣例預測外部投影，再從影像扣除，剩餘資料才更接近局部區域的觀測。外部模型或角度估計不準時，殘留與過度扣除也會形成誤差，因此要檢查扣除後的影像及其重建是否合理。

Mask 通常使用柔邊，以降低 Fourier 空間的振鈴；

## 33. `book/06_heterogeneity.md` 第 184 行

原句：

異質性分析讓這個問題更尖銳，因為每個類別或每個變化成分都只用到一部分資料，而類別的劃分本身就是從資料估出來的。

原因：連續成分不是每個只用部分影像。

改後：

異質性分析增加了需要估計的量：離散類別分攤粒子，各類的有效資料量減少；連續成分則可由全部影像共同估計，但必須同時估計模式與各粒子的座標。分類和座標本身的不確定性也要納入判讀。

## 34. `book/06_heterogeneity.md` 第 189 行

原句：

及 150–200 kDa 的規模下限

原因：延伸閱讀同步歷史範圍。

改後：

及當時 150–200 kDa 的經驗尺度

## 35. `book/06_heterogeneity.md` 第 193 行

原句：

同一組雙族群資料在雜訊較大時無法與單族群分布區分

原因：圖示分辨困難不等於理論不可識別。

改後：

有限樣本的雙族群資料在雜訊較大時難以與單族群分布區分

## 36. `book/06_resolution_validation.md` 第 12 行

原句：

\mathbb E\left[\left(\widehat V-V\right)^2\right]
=\underbrace{\left(\mathbb E[\widehat V]-V\right)^2}_{\text{偏差平方}}
+\underbrace{\operatorname{Var}\left(\widehat V\right)}_{\text{變異}} .

原因：複數Fourier係數須使用模平方。

改後：

\mathbb E\left[\left|\widehat V-V\right|^2\right]
=\underbrace{\left|\mathbb E[\widehat V]-V\right|^2}_{\text{偏差平方}}
+\underbrace{\operatorname{Var}\left(\widehat V\right)}_{\text{變異}} .

## 37. `book/06_resolution_validation.md` 第 17 行

原句：

第一項描述「如果重複收集很多次資料，平均的估計偏離真值多少」；第二項描述不同批資料之間的波動

原因：明定complex variance，保留實數例子。

改後：

此處 $\operatorname{Var}(\widehat V)=\mathbb E|\widehat V-\mathbb E\widehat V|^2$，適用於實數或複數係數。第一項描述重複收集資料後，平均估計偏離真值的程度；第二項描述不同批資料之間的波動

## 38. `book/06_resolution_validation.md` 第 19 行

原句：

這個分解之所以關鍵，是因為兩項對資料量的反應完全不同。變異隨影像數增加而下降；偏差不會。

原因：偏差可以隨N改變，不能泛稱永不下降。

改後：

增加資料常能降低變異，但偏差如何改變，要看估計方法與誤差來源；有些偏差會隨資料增加而縮小，有些則會持續存在。

## 39. `book/06_resolution_validation.md` 第 21 行

原句：

收集更多資料把變異壓到接近零，偏差原封不動地留在那裡。

原因：限定固定b示例。

改後：

在這個固定 $b$ 的例子裡，增加資料會把變異壓低，卻不會消除預設的偏差。

## 40. `book/06_resolution_validation.md` 第 23 行

原句：

Sorzano 等人主張，冷凍電鏡裡被稱為 overfitting 的現象，主要屬於偏差而非變異：以 100,000 顆 $200\times200$ 的粒子重建一個 $200^3$ 的體，量測數遠多於未知數，體本身沒有多少自由度去擬合雜訊。真正造成錯誤結構的，是模型與資料的不匹配、姿態與類別參數的錯誤估計，以及目標函數或演算法本身帶有的偏差 {cite}`sorzano2022`。這個判斷改變了驗證的重點：要檢查的對象是參數估計是否正確，而非最後那張圖是否「太銳利」。



原因：保留作者bias觀點但不接受量測數多必無過擬合。

改後：

Sorzano 等人強調，SPA 中一些被稱為 overfitting 的現象來自系統性偏差：例如以 100,000 顆 $200\times200$ 粒子重建一個 $200^3$ 的體，表面量測數雖遠多於 voxel 數，角度與類別估計錯誤、成像模型不符，以及演算法偏差仍可能產生錯誤結構 {cite}`sorzano2022`。這是作者用來提醒讀者檢查參數與模型的觀點。

量測數多，並不能單獨保證可識別或排除過擬合：投影可能冗餘、方向覆蓋可能不足，角度與類別也需從含雜訊資料估計。驗證因此要同時檢查參數、資料獨立性與重建結果，不能只憑圖像銳利或量測數較多就判定可靠。



## 41. `book/06_resolution_validation.md` 第 38 行

原句：

同一個物理粒子的重複觀測、同一部 movie 或同一張 micrograph 的粒子，彼此的雜訊相關。

原因：同一micrograph不必然所有noise相關，條件化說明。

改後：

同一個物理粒子的重複觀測會共享資訊；同一部 movie 或同一張 micrograph 中的粒子，也可能因共同背景或校正誤差而相關。

## 42. `book/06_resolution_validation.md` 第 54 行

原句：

FSC 之所以能換算成訊雜比，靠的是一個簡單的雜訊模型。

原因：正文呈現核心換算，區分母體rho與樣本FSC。

改後：

要把 FSC 換算成訊雜比，需假設兩張 half-map 共享訊號，而各自的雜訊零平均、彼此獨立，且具有相同功率。令 $P_S$ 為訊號功率、$P_N$ 為單張 half-map 的雜訊功率，母體相關為

$$
\rho(s)=\frac{P_S}{P_S+P_N}.
$$

FSC 是用有限 shell 係數估計 $\rho$ 的樣本量，並非每次都等於這個功率比。有效樣本足夠時，才可用 FSC 近似 $\rho$，得到

$$
\mathrm{SSNR}_{\mathrm{half}}(s)=\frac{\rho}{1-\rho},
\qquad
\mathrm{SSNR}_{\mathrm{full}}(s)=\frac{2\rho}{1-\rho}.
$$

第二式對應兩張 half-map 的平均：獨立且等功率的雜訊經平均後，功率減半。有限樣本的 FSC 可能為負，這不代表真實 SSNR 為負，而是該頻帶的相關估計受抽樣波動影響，不宜直接套入訊雜比換算 {cite}`penczek2010resolution`。

## 43. `book/06_resolution_validation.md` 第 86 行

原句：

其中 $n_s$ 是 shell 內的係數個數；所有交叉項因獨立與零平均而消失。代入定義得

原因：不能將分子分母期望代入而宣稱隨機比值恆等。

改後：

其中 $n_s$ 是 shell 內的係數個數；交叉項的期望因獨立與零平均而消失。用這些母體功率定義的相关為

## 44. `book/06_resolution_validation.md` 第 89 行

原句：

\rho\equiv\mathrm{FSC}(s)=\frac{P_S}{P_S+P_N}.

原因：移除sample FSC與population rho恆等。

改後：

\rho(s)=\frac{P_S}{P_S+P_N}.

## 45. `book/06_resolution_validation.md` 第 92 行

原句：

定義單張 half-map 的頻譜訊雜比

原因：補足ratio-estimator推論界線。

改後：

FSC 的分子與分母在有效樣本增加時可接近各自期望，因而近似上式；這不是有限樣本下「比值的期望等於期望的比值」。定義單張 half-map 的頻譜訊雜比

## 46. `book/06_resolution_validation.md` 第 111 行

原句：

定義完整 map 與無誤差的參考結構之間的相關

原因：Cref同樣需population模型。

改後：

在同一訊號加雜訊模型下，完整 map 與無誤差參考結構之間的母體相關為

## 47. `book/06_resolution_validation.md` 第 114 行

原句：

C_{\mathrm{ref}}(s)=\sqrt{\frac{2\,\mathrm{FSC}(s)}{1+\mathrm{FSC}(s)}},

原因：Cref以rho寫精確式，FSC僅plug-in估計。

改後：

C_{\mathrm{ref}}(s)=\sqrt{\frac{2\rho(s)}{1+\rho(s)}},

## 48. `book/06_resolution_validation.md` 第 117 行

原句：

這個量相當於晶體學的 figure of merit。代入 $\mathrm{FSC}=0.143$：$2\times0.143=0.286$，除以 $1.143$ 得 $0.2502$，開根號得 $0.5002$。**half-map FSC 取 0.143，正好對應完整 map 與真實結構的相關為 0.5**，也對應該 shell 的平均相位誤差 $60^\circ$（因為 $\cos60^\circ=0.5$）{cite}`rosenthal2016`。

原因：1/7非精確.143；相關.5不能推出平均相位60度。

改後：

以樣本 $\mathrm{FSC}=0.143$ 近似 $\rho$：$2\times0.143=0.286$，除以 $1.143$ 得 $0.2502$，開根號得 $0.5002$。精確令 $C_{\mathrm{ref}}=0.5$ 時，對應的是 $\rho=1/7\approx0.143$。因此慣用的 0.143 門檻，約對應完整 map 與真實訊號的相關為 0.5 {cite}`rosenthal2016`。這是相關量的換算，不能由它推出平均角度或平均相位誤差。

## 49. `book/06_resolution_validation.md` 第 142 行

原句：

**共同偏差在 FSC 裡的表現和真實訊號完全一樣**：它讓兩張圖更一致，因此讓 FSC 上升、讓報告的解析度變好。

原因：共bias可建設或抵消，不能說必增FSC。

改後：

FSC 無法把共同成分中的真實訊號與偏差分開。因此即使偏差很大，FSC 仍可能很高。不過，偏差不必然提高 FSC：若 $b=-S_{\mathrm{true}}$，共同訊號反而被抵消。

## 50. `book/06_resolution_validation.md` 第 146 行

原句：

這個結論的用途是把 half-set 放回正確的位置：它有效地防止兩邊互相學習對方的雜訊，對兩邊共有的系統性錯誤則無能為力。

原因：保留gold standard用途，作者批判不寫成一般替代共識。

改後：

上述是作者對偏差的論點，不代表應放棄 half-set。保持兩組獨立仍能降低共享高頻雜訊造成的樂觀估計；共同偏差則需要改變初始模型、方法或外部證據來檢查。

## 51. `book/06_resolution_validation.md` 第 156 行

原句：

更危險的情形是 mask 本身帶有結構資訊。

原因：區分same fixed mask與data-dependent mask。

改後：

事先固定的同一個 mask，並不會使原本獨立的兩組零平均雜訊自動產生交叉相關；前面的高頻上升主要來自共同訊號的頻率混合。另一個問題是 mask 的形狀由資料本身決定。

## 52. `book/06_resolution_validation.md` 第 158 行

原句：

校正的方式是把 mask 造成的相關估出來再扣掉。標準流程是：選一個明顯低於 unmasked 解析度的截止頻率，把兩張 half-map 在該頻率以上的相位隨機化；對這兩張相位隨機化的圖套上同一個 mask，重新計算 FSC。若沒有卷積效應，這條曲線在截止頻率以上應為零，因此任何非零的相關都來自 mask；用它去校正未做相位隨機化的 masked FSC，得到的 corrected 曲線就只反映移除溶劑雜訊帶來的改善 {cite}`scheres2016`。

這個檢查有一個細節值得知道：原始構想是在**實驗影像**層級做相位隨機化，而常見的實作是在**體**的層級做，後者測到的主要是所施加操作本身的邏輯結果，驗證力較弱 {cite}`sorzano2022`。判讀時把 masked、unmasked 與 corrected 三條曲線並排看，比只讀一個數字有用得多。



原因：加入正確mask normalized correction、Hermitian/independent phases、頻段及用途區分；不用未讀Chen2013全文作已驗來源。

改後：

可用高頻相位隨機化估計 mask 造成的頻率混合。選定截止頻率後，對兩張 half-map 在該頻率以上的相位**分別獨立隨機化**，並保留 Hermitian 對稱，使反轉換仍得到實數密度。再套用同一個 mask，得到 $\mathrm{FSC}_{\mathrm{randomized}}$。在截止頻率以上，尚未隨機化的共同低頻訊號仍可經 mask 混入高頻，因此這條曲線可估計該效應；有限樣本也會使曲線波動。

校正式不是兩條 FSC 直接相減，而是

$$
\mathrm{FSC}_{\mathrm{corrected}}
=\frac{\mathrm{FSC}_{\mathrm{masked}}-\mathrm{FSC}_{\mathrm{randomized}}}
{1-\mathrm{FSC}_{\mathrm{randomized}}}.
$$

例如 masked FSC 為 0.5、randomized FSC 為 0.2，校正後是 $(0.5-0.2)/(1-0.2)=0.375$。[Rich Hite 的公開講義第 11 頁](https://semc.nysbc.org/wp-content/uploads/2018/08/20190311_SEMCcourse_SPAIII-RH.pdf)列出這個公式。它應用在相位已隨機化且避開截止頻率過渡影響的頻帶；較低頻仍使用未隨機化的結果。若 randomized FSC 接近 1，分母很小，校正會不穩定，應重新檢查 mask 與截止頻率。

還要分清楚兩種檢查的目的：在**粒子影像**層級做高頻 noise substitution 並重跑 refinement，可以測試演算法是否把無訊號的頻率學成共同結構；在 **half-map** 層級做相位隨機化，主要用來估計後處理 mask 的頻率混合。後者不能取代前者的流程驗證，但對 mask 校正仍有明確用途 {cite}`sorzano2022,scheres2016`。

將 masked、unmasked、randomized 與 corrected 曲線一起看，才能區分去除溶劑雜訊後的改善與 mask 引入的效應。校正並不保證消除所有共同偏差。



## 53. `book/06_resolution_validation.md` 第 196 行

原句：

因此銳化要和頻率相依的加權一起使用。最小均方誤差意義下的最佳線性濾波增益，可以直接用 half-map FSC 表示：

原因：Wiener目標為已衰減map訊號，不能混同envelope前真值。

改後：

因此銳化通常配合頻率相依的加權，並限制在有可靠訊號的範圍。在 $F=S+N$ 模型下，若目標是估計 map 中的訊號 $S$，最小均方誤差的線性增益為：

## 54. `book/06_resolution_validation.md` 第 203 行

原句：

{cite}`penczek2010resolution`。這個增益在 $\rho\to0$ 的頻帶趨近於零，正好抵消銳化在無訊號區域的放大。順帶一提，它恰好等於上面 $C_{\mathrm{ref}}$ 的平方——最佳增益與「與真值的相關」之間有這層直接關係。

原因：Wiener不正好抵消任意sharpening；補envelope-aware增益。

改後：

{cite}`penczek2010resolution`。實作時可用 FSC 近似 $\rho$；弱訊號頻帶的權重較小，但它不會正好抵消任意銳化因子。這個增益也等於同一模型下的 $C_{\mathrm{ref}}^2$。

如果目標改成包絡衰減前的真實係數 $T$，模型就應寫成 $F=E T+N$。已知包絡 $E$，且 $T,N$ 不相關時，估計 $T$ 的線性增益為

$$
\frac{E^*P_T}{|E|^2P_T+P_N}.
$$

它還需要包絡與衰減前的訊號功率，不能只由 $2\rho/(1+\rho)$ 決定。這也說明為何銳化、弱頻帶收縮及可靠頻率上限要一起判讀。

## 55. `book/06_resolution_validation.md` 第 217 行

原句：

若大多數估計彼此接近，較可能處於隨機誤差的情形，此時把這些估計平均起來可以降低變異；若它們分散，則指出某處出了問題

原因：多次一致不能推定較可能無bias；平均降變異需相關性條件。

改後：

若估計彼此接近，表示結果在這些執行條件下穩定，但也可能共享偏差；若估計分散，則需要檢查初始化敏感性、資料不足或模型不符。將相容的估計平均可能降低變異，效果還取決於估計間的相關性

## 56. `book/06_resolution_validation.md` 第 231 行

原句：

並指出 half-map FSC 的 0.143 與 $C_{\mathrm{ref}}$ 的 0.5 及 $60^\circ$ 相位誤差對應同一個解析度——本章 0.143 的換算即取自此處

原因：延伸閱讀刪除不成立的60度推論。

改後：

並將 half-map FSC 的約 0.143 對照到 $C_{\mathrm{ref}}$ 的約 0.5；閱讀時要分清相關係數與平均相位誤差，不能直接把相關值取反餘弦當成平均誤差

## 57. `book/06_resolution_validation.md` 第 234 行

原句：

每筆密度圖都附有 FSC 曲線、使用的 mask、報告的解析度與處理流程描述。

原因：EMDB非每筆必附完整mask/FSC。

改後：

先查看條目實際提供哪些 FSC 曲線、mask、解析度與處理流程中繼資料；不同條目的資料完整度可能不同。

## 58. `06_reconstruction_validation.md` 第 69 行（回讀追加）

原句：若 $A$ 欠秩，就存在非零 $V_0$

原因：欠定系統可滿行秩仍有nullspace，必須比較列向量數即未知係數數。

改後：若 $A$ 的秩小於未知係數數量，就存在非零 $V_0$

## 59. `06_reconstruction_validation.md` 第 70 行（回讀追加）

原句：若 $A$ 滿秩，但某些奇異值很小

原因：唯一性需要full column rank。

改後：若 $A$ 的秩等於未知係數數量，但某些奇異值很小

## 60. `06_reconstruction_validation.md` 第 97 行（回讀追加）

原句：觀測模型是 $y_j=H_jV+\varepsilon_j$。

原因：重算無偏variance需明定noise零均與獨立。

改後：觀測模型是 $y_j=H_jV+\varepsilon_j$，各 $\varepsilon_j$ 零平均且互相獨立。

## 61. `06_reconstruction_validation.md` 第 148 行（回讀追加）

原句：兩張投影的中心切片必定在三維 Fourier 空間相交於一條直線

原因：同方向切片可能重合，不只共線。

改後：兩張不同投影方向的中心切片，在三維 Fourier 空間相交於一條直線

## 62. `06_reconstruction_validation.md` 第 225 行（回讀追加）

原句：說明隨機梯度下降為何能從任意初始化出發

原因：不能宣稱任何初始化都成功。

改後：說明該演算法如何從隨機初始化出發

## 63. `06_reconstruction_validation.md` 第 202 行（回讀追加）

原句：對半徑為 $s$ 的 Fourier shell，Fourier shell correlation 定義為

原因：保留完整shell公式並明示實數密度共軛對稱使交叉和為實數。

改後：取包含正負共軛頻率的完整對稱 shell，半徑記為 $s$。Fourier shell correlation 定義為

## 64. `06_heterogeneity.md` 第 11 行（回讀追加）

原句：並使用線性重建。在這個簡化下

原因：精確β平均还需要忽略重建誤差，線性不等於精確反演。

改後：並使用線性重建，暫不考慮雜訊與重建誤差。在這個簡化下

## 65. `06_heterogeneity.md` 第 56 行（回讀追加）

原句：類別標籤 $c_i$ 是另一個隱藏變數

原因：避免z原章包含類別而此式重複c，並區分後續連續boldz。

改後：這一式的 $z_i$ 只表示角度與位移；類別另外用 $c_i$ 表示。類別標籤 $c_i$ 是另一個隱藏變數

## 66. `06_heterogeneity.md` 第 75 行（回讀追加）

原句：**$K$ 增大時，估計本身也會變差。**

原因：增加正確類別可改善估計，不應泛稱K增必變差。

改後：**$K$ 增大時，每類可用的資料可能變少。**

## 67. `06_heterogeneity.md` 第 168 行（回讀追加）

原句：讓匹配更著重 mask 內的質量

原因：匹配影像訊號是density，後文物理質量門檻保留。

改後：讓匹配更著重 mask 內的密度

## 68. `06_heterogeneity.md` 第 178 行（回讀追加）

原句：用同一組角度、位移、CTF 與強度慣例預測外部投影

原因：與前文顯式alpha模型一致，subtraction算子需包含scale。

改後：這裡的 $A_i$ 也納入已估計的強度比例，並使用與影像相同的角度、位移和 CTF 預測外部投影

## 69. `06_resolution_validation.md` 第 43 行（回讀追加）

原句：對半徑為 $s$ 的 Fourier shell，Fourier shell correlation 定義為

原因：完整對稱shell保證此未取實部的公式為實數。

改後：取包含正負共軛頻率的完整對稱 shell，半徑記為 $s$。Fourier shell correlation 定義為

## 70. `06_resolution_validation.md` 第 60 行（回讀追加）

原句：FSC 是用有限 shell 係數估計 $\rho$ 的樣本量

原因：sample statistic不是sample size。

改後：FSC 是用有限 shell 係數估計 $\rho$ 的統計量

## 71. `06_resolution_validation.md` 第 21 行（回讀追加）

原句：兩項相當

原因：0.04是0.01四倍，說相當不準確。

改後：偏差已大於變異

## 72. `06_resolution_validation.md` 第 229 行（回讀追加）

原句：給出 $\mathrm{SSNR}=\mathrm{FSC}/(1-\mathrm{FSC})$ 與 half-set 版本的 $2\mathrm{FSC}/(1-\mathrm{FSC})$

原因：延伸閱讀half/full語意清楚，不逆轉原式係數。

改後：給出單張 half-map 的 $\mathrm{SSNR}=\mathrm{FSC}/(1-\mathrm{FSC})$ 與兩半平均之 full-map 的 $2\mathrm{FSC}/(1-\mathrm{FSC})$

## 73. `06_resolution_validation.md` 第 229 行（回讀追加）

原句：以 Fisher $z$ 變換建立顯著性檢定的必要性

原因：Fisher transform是方法之一，非所有檢定都必須。

改後：以 Fisher $z$ 變換近似相關係數抽樣分布的做法

## 74. `06_resolution_validation.md` 第 119 行（回讀追加）

原句：通常要先做 Fisher 的 $z$ 變換

原因：避免Fisher普遍必要。

改後：可用 Fisher 的 $z$ 變換作近似

## 75. `06_resolution_validation.md` 第 230 行（回讀追加）

原句：本章關於共同偏差的討論全部出自此文，且可公開取得全文

原因：後續加入獨立反例，不能再稱全部出自此文。

改後：此文是本章共同偏差討論的主要參考，可公開取得全文

## 最終算式與檔案核驗

以下為固定參數的直接計算，無隨機模擬。MAP 以兩個固定 2×2 複數算子驗证非對角耦合及共軛轉置；LS使用固定3×2矩陣與scale=1.7。

```json
{
  "rectangular_singular_values": [
    1.7320508075688772,
    1.0
  ],
  "ctf_map": {
    "mean": 0.5475113122171945,
    "variance": 0.24774267521140025,
    "bias_squared": 0.20474601257140523,
    "mse": 0.4524886877828055,
    "unregularized_mse": 0.8264462809917356
  },
  "two_state_amplitude": {
    "30": 0.8660254037844387,
    "20": 0.7071067811865476,
    "10": 6.123233995736766e-17,
    "5": 1.0
  },
  "gaussian_amplitude_sigma2": {
    "10": 0.4540407387272451,
    "5": 0.04249905628536256
  },
  "mask_correction": 0.37499999999999994,
  "correlation_exact_one_seventh": 0.5,
  "correlation_rounded_fsc": 0.5002186748409814,
  "fixed_bias_mse": {
    "100": 1.04,
    "10000": 0.05000000000000001,
    "1000000": 0.04010000000000001
  },
  "map_normal_gradient_norm": 1.249000902703301e-16,
  "map_normal_offdiagonal_abs": 0.3490902462114918,
  "latent_ls_error": 0.0,
  "source_structure": "3 chapters: original anchors, figures, math/fence pairing, no pose term or comprehension blocks pass"
}
```

科學修正75處、純語句潤稿49處。已逐章回讀；未建置網站或執行notebook，交由主代理整合驗收。Hite講義公式沿用主代理已核實的第11頁；未聲稱讀取Chen2013全文。
