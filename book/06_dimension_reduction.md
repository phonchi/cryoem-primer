# SPA 高維影像與降維：PCA、MPCA 與 2SDR

上一章把粒子影像對齊、分群，讓共同訊號能透過平均顯現。接著會遇到另一個問題：一張 $130\times130$ 的小影像，仍有 16,900 個像素座標。比較成千上萬張低訊雜比影像時，大量雜訊也一起進入距離與共變異數。能否用少量係數描述真正反覆出現的變化？

降維的想法是從整批資料學出共用的表示，再用每張影像的少量係數描述差異。本章先從 PCA 的幾何意義出發，再保留影像的矩陣結構，引入 MPCA 與兩階段降維（2SDR）。這些方法可用來表示影像及去雜訊；用在 SPA 時，還要分辨取向、CTF、位移與構形所造成的變化。

(spa-high-dimension)=
## 像素很多，距離中也會累積很多雜訊

將影像向量寫成 $y_i=s_i+\varepsilon_i\in\mathbb R^d$。若不同影像的雜訊互相獨立，且 $\varepsilon_i\sim N(0,\sigma^2 I_d)$，則

$$
\mathbb E\|y_i-y_j\|^2
=\|s_i-s_j\|^2+2d\sigma^2.
$$

這裡把兩張乾淨訊號 $s_i,s_j$ 固定，只對雜訊取期望。當 $d=16{,}900$、$\sigma^2=1$，即使兩張乾淨影像完全相同，觀測平方距離的期望仍為 33,800。訊號造成的距離差，要在這個背景上辨認。

若雜訊具有共變異數 $\Sigma$，上述雜訊項改為 $2\operatorname{tr}(\Sigma)$；其中的方向性還會影響鄰近關係。因此，相似度的選擇、whitening、對齊與降維，都會改變分群時影像之間的距離與鄰近關係。這些步驟的順序與估計方式，需要配合實際模型。

## PCA 找的是整批影像共同變動的方向

令 $\bar y=n^{-1}\sum_i y_i$，以中心化影像建立樣本共變異數

$$
S=\frac1n\sum_{i=1}^{n}(y_i-\bar y)(y_i-\bar y)^{\mathsf T}.
$$

此處採 $1/n$，方便和最小平方目標對照；改用 $1/(n-1)$ 會等比例縮放特徵值，特徵向量保持相同。對任何長度為 1 的向量 $u$，投影係數 $u^{\mathsf T}(y_i-\bar y)$ 的樣本變異數為 $u^{\mathsf T}Su$。使它最大的方向就是最大特徵值對應的特徵向量。

保留前 $r$ 個正交方向，排成 $U_r\in\mathbb R^{d\times r}$，每張影像的係數與重建為

$$
a_i=U_r^{\mathsf T}(y_i-\bar y),\qquad
\widehat y_i=\bar y+U_r a_i.
$$

$U_r$ 的每一欄都可以排回一張影像，稱為 eigenimage。它描述一種帶正負號的變化模式；$a_i$ 指出第 $i$ 張影像含有多少這種變化。基底的正負號可以整欄反轉，同時反轉係數而不改變重建。

考慮只有兩個像素的例子：

$$
S=\begin{pmatrix}5&4\\4&5\end{pmatrix}.
$$

特徵值為 9 與 1，對應方向是 $(1,1)/\sqrt2$ 與 $(1,-1)/\sqrt2$。兩個像素主要一起變亮或變暗；相反方向的變化較小。只留第一個方向，可保留這批資料 90% 的總變異。這個比例描述的是**觀測變異**；要把剩下的 10% 解讀成雜訊，還需要雜訊模型。

```{figure} images/spa/pca_geometry.png
:name: fig-spa-pca
:alt: 二像素資料的主要變動沿對角線，PCA 將資料投影至第一主成分；旁圖比較各秩的重建偏差與保留雜訊。
:width: 92%

左圖以兩個像素表示一張影像，橘線為第一主成分。右圖使用已知基底與訊號能量的簡化例子，顯示增加秩時，遺失訊號減少，保留的雜訊增加。
```

(pca-reconstruction)=
## 最大變異與最小重建誤差如何連在一起？

PCA 也可以由重建問題定義：

$$
\min_{U^{\mathsf T}U=I_r,\{a_i\}}
\sum_i\|(y_i-\bar y)-Ua_i\|^2.
$$

固定 $U$ 後，最佳係數是正交投影 $a_i=U^{\mathsf T}(y_i-\bar y)$。每張影像的能量分成保留部分與垂直於子空間的殘差，因此最小化殘差，相當於最大化投影後的總變異。這讓「主要變化方向」與「最佳低秩平方重建」成為同一個問題。

把中心化影像排成 $Y\in\mathbb R^{n\times d}$，使用 thin SVD：

$$
Y=L D U^{\mathsf T},\qquad S=U(D^2/n)U^{\mathsf T}.
$$

可以從資料矩陣直接求奇異向量，不必先建立 $d\times d$ 共變異數矩陣。當 $n<d$，中心化資料的秩至多為 $n-1$，這也說明可估計方向數量受到樣本限制。比較 PCA 與矩陣方法的計算量時，要說清楚採用完整共變異數、thin SVD 或近似分解；不同實作的成本差很多。

```{dropdown} PCA 最小重建誤差的推導

令 $P=UU^{\mathsf T}$，其中 $U^{\mathsf T}U=I_r$。因 $P=P^{\mathsf T}=P^2$，

$$
\frac1n\sum_i\|(I-P)(y_i-\bar y)\|^2
=\operatorname{tr}(S)-\operatorname{tr}(U^{\mathsf T}SU).
$$

對稱半正定矩陣 $S$ 的前 $r$ 個特徵向量最大化右側第二項。若特徵值由大到小為 $\lambda_1,\ldots,\lambda_d$，最小的每張影像平均殘差平方為 $\sum_{j>r}\lambda_j$。這個量針對訓練觀測影像；相對於未知乾淨影像的誤差，則要另外評估。
```

## 低秩為何可能去雜訊？又會遺失什麼？

先假設乾淨影像位於一個已知的 $r$ 維子空間，平均數也已知。將白雜訊影像投影到該子空間，留下的雜訊能量為

$$
\mathbb E\|U_rU_r^{\mathsf T}\varepsilon_i\|^2=r\sigma^2,
$$

原本則為 $d\sigma^2$。共用的訊號方向得到保留，垂直於子空間的雜訊被移除。若真正訊號有一部分在子空間外，投影也會一併刪掉它。

對固定子空間及固定乾淨訊號 $s_i$，兩者可寫在同一式中：

$$
\mathbb E\|P_r(s_i+\varepsilon_i)-s_i\|^2
=\|(I-P_r)s_i\|^2+r\sigma^2.
$$

例如沿四個已知正交方向的訊號能量為 $(9,4,0.2,0)$，$\sigma^2=1$。保留 $r=1,2,3,4$ 時，期望誤差分別是 $5.2,2.2,3,4$。第三個方向只有 0.2 的訊號能量，加入它卻多保留 1 的雜訊能量，因此這個例子的最佳秩為 2。

實際 PCA 的基底由同一批含雜訊資料估計，會有額外的子空間估計誤差；上式是固定基底下的教學分解。高維且樣本不足時，最大樣本特徵值也可能受雜訊推高。增加秩本來就會降低訓練殘差，因此秩選擇還需要評估保留雜訊的代價。

(spa-mpca)=
## MPCA 保留影像的列與行結構

PCA 將 $p\times q$ 影像排列成 $pq$ 維向量。MPCA 則用兩組共用基底表示整批矩陣：

$$
X_i=M+A Z_i B^{\mathsf T}+E_i,
$$

其中 $A\in\mathbb R^{p\times r_p}$、$B\in\mathbb R^{q\times r_q}$，且 $A^{\mathsf T}A=I$、$B^{\mathsf T}B=I$；$M$ 是平均影像，$Z_i\in\mathbb R^{r_p\times r_q}$ 是第 $i$ 張影像的壓縮係數。固定基底後，

$$
Z_i=A^{\mathsf T}(X_i-M)B,\qquad
\widehat X_i=M+A Z_iB^{\mathsf T}.
$$

$A$ 與 $B$ 分別沿影像的兩個維度壓縮。它們由整批資料共同估計，讓所有影像的係數具有相同意義；對每張影像各自做一次 SVD，得到的個別基底通常不同。

估計 $A,B$ 時，可最小化所有影像的總重建誤差。固定 $B$ 後，以

$$
\sum_i (X_i-M)BB^{\mathsf T}(X_i-M)^{\mathsf T}
$$

的前 $r_p$ 個特徵向量更新 $A$；固定 $A$ 後，交換矩陣方向更新 $B$，如此交替。這是耦合的子空間問題，單次分別計算列、行共變異數通常尚未完成這個最佳化 {cite}`chung2020`。

矩陣表示把基底限制在可分離的列、行子空間，降低表示成本，同時加入結構假設。影像若在各方向有複雜旋轉，可能需要較大的 $(r_p,r_q)$ 才能保留變化；比較方法時，要同時看參數數量與保留下來的訊號。

```{dropdown} 用 Kronecker product 看 MPCA 與 PCA 的關係

採用依欄堆疊的 $\operatorname{vec}$ 定義，有

$$
\operatorname{vec}(A Z_i B^{\mathsf T})=(B\otimes A)\operatorname{vec}(Z_i).
$$

因此 MPCA 使用的向量子空間由 $B\otimes A$ 張成，維度為 $r_pr_q$。一般 PCA 的同維基底可以在更大的候選集合中搜尋；MPCA 以可分離結構限制基底。$Z_i$ 的不同元素仍可相關，這正是第二階段還有進一步壓縮空間的原因。換用依列堆疊時，Kronecker 因子的順序也要跟著調整。
```

(spa-2sdr)=
## 2SDR：先縮小矩陣，再學係數之間的共同變化

兩階段降維先由 MPCA 取得 $Z_i$，再對 $u_i=\operatorname{vec}(Z_i)$ 做 PCA。令 $G_r\in\mathbb R^{r_pr_q\times r}$ 為第二階段的正交基底，則

$$
a_i=G_r^{\mathsf T}u_i,\qquad
\widehat X_i=M+A\operatorname{unvec}(G_r a_i)B^{\mathsf T}.
$$

第一階段利用影像的列、行結構；第二階段再移除壓縮係數中冗餘的方向。譬如 $130\times130$ 影像先壓成 $20\times20$ 矩陣，再用 40 個主成分表示，資料維度依序為 $16{,}900\to400\to40$。這組數字說明每一步的維度如何改變；實際使用的秩，應由資料與方法準則決定。

Chung 等人的 2SDR 原始方法包含完整的秩選擇程序：先估計雜訊、以 SURE 選第一階段的矩陣秩，再用 GIC 選第二階段 PCA 秩。只指定兩個手選秩並串接 MPCA、PCA，是兩階段表示的示範；要重現論文方法，還需實作其雜訊估計與兩種秩準則 {cite}`chung2020`。

## 如何選秩：從訓練誤差走向估計風險

增加係數數量幾乎總能改善訓練重建，因此秩選擇需要衡量額外自由度的代價。對 $y=x+\varepsilon$、$\varepsilon\sim N(0,\sigma^2I_d)$ 及適當可微的估計 $f(y)$，Stein's unbiased risk estimate（SURE）使用

$$
\operatorname{SURE}(y)
=\|f(y)-y\|^2+2\sigma^2\operatorname{div}f(y)-d\sigma^2.
$$

其中散度（divergence）定義為 $\operatorname{div}f(y)=\sum_{j=1}^{d}\partial f_j(y)/\partial y_j$，累加輸入各座標的小變化對相應估計值的影響。若 $f(y)=Py$ 是固定的正交投影，便有 $\operatorname{div}f=\operatorname{tr}(P)=r$，接回前面保留 $r$ 個方向的例子。基底由資料估計時，$P$ 也隨輸入改變，微分需要把這部分算進去。

SURE 的期望等於相對於乾淨訊號的平方風險。第一項是觀測殘差，第二項補償模型跟隨資料的靈活程度，最後一項扣除觀測本身的雜訊能量。當基底由資料估計，微分必須包含基底對資料的依賴，不能直接把自由度替換成「最後留下幾個係數」。2SDR 的第一階段準則正是依 MPCA 模型處理這個問題。

GIC 則比較模型的擬合程度與估計複雜度。原文的第二階段以壓縮係數的共變異數特徵值建立準則，並討論模型下的秩選擇性質。這些結果依賴論文的雜訊、子空間與特徵值條件；將其套到未校正 CTF 或有色雜訊資料時，需要重新檢視雜訊與子空間條件，才能判斷準則是否適用。

也可以使用{doc}`07_synthetic_data`的合成粒子，利用已知真值分開量測保留的訊號與雜訊。真實資料則可用獨立資料檢查表示穩定性，並把預處理和基底估計納入資料分割的範圍。

## 在 SPA 中，主成分究竟代表什麼？

假設同一結構的兩組影像只差一點平移。對小位移 $\delta$，一維截面的近似為

$$
f(x-\delta)\approx f(x)-\delta f'(x).
$$

即使所有粒子的構形完全相同，位移變動仍會產生近似導數形狀的主成分。旋轉、離焦、亮度、冰厚與背景也會產生大幅變異。主成分能指出資料中主要的變化方向；解釋為構形前，要把它與這些成像因素對照。

實際判讀可把主成分係數對上已知的離焦、取向、位移與拍攝批次，並檢查沿主成分改變時的影像。如果係數主要跟離焦一起變化，這個方向首先反映的是成像條件。將已對齊且條件相近的粒子分組，有助於縮小干擾，但對齊誤差仍會留在資料中。

同樣地，對 raw particle images 做 PCA，得到的是二維影像變異；{doc}`06_heterogeneity`中的 3DVA 則依已估計的角度投影三維線性變化，再套用位移與 CTF，直接和粒子影像比較。兩者的觀測模型不同，座標的物理意義也需要各自建立 {cite}`punjani2021`。

## 如何判讀二維嵌入圖？

將 PCA 或 2SDR 係數交給 t-SNE 等方法，可以在平面上查看鄰近關係。嵌入圖的群集外觀還受到距離定義、鄰域參數與最佳化影響；圖上的群間距離、面積及空白區，需要回到原表示空間檢查。

若圖上出現一個小群，可先看該群原始粒子、類別平均與離焦分布，再比較不同參數下是否保留相同成員。清楚的分群提示可能存在結構差異；要確認這些差異，還需檢查資料模型，並做獨立重現與三維分析。這也把降維接回上一章的分類穩定性問題。

接下來的{doc}`06_reconstruction_validation`把目標從二維影像表示推進到三維結構，會使用成像算子來說明資料如何約束未知密度。

## 延伸閱讀

- **Chung, S.-C., Wang, S.-H., Niu, P.-Y., Huang, S.-Y., Chang, W.-H. 與 Tu, I.-P.（2020）．*Two-Stage Dimension Reduction for Noisy High-Dimensional Images and Application to Cryogenic Electron Microscopy*. Annals of Mathematical Sciences and Applications, 5(2), 283–316。** [DOI](https://doi.org/10.4310/AMSA.2020.v5.n2.a4)。依序讀第 1.3 節的 PCA／MPCA、第 2 節與第 2.1 節的四步程序及 Figure 1（頁 288–297），再讀第 3 節實驗。SURE 與 GIC 的角色要和兩階段表示一起閱讀 {cite}`chung2020`。
- **Chen, T.-L. 等（2014）．*γ-SUP: A Clustering Algorithm for Cryo-Electron Microscopy Images of Asymmetric Particles*. The Annals of Applied Statistics, 8(1), 259–285。** [開放預印本](https://arxiv.org/abs/1205.2034)；讀第 2.3 節與第 3 節的估計及 self-updating 機制，對照降維之後如何處理高雜訊與離群影像。γ-SUP 的更新有自己的穩健估計目標 {cite}`chen2014`。
- **Singer, A. 與 Sigworth, F. J.（2020）．*Computational Methods for Single-Particle Electron Cryomicroscopy*. Annual Review of Biomedical Data Science, 3, 163–190。** [開放全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8412055/)；閱讀異質性相關段落，辨認觀測影像的共變異數如何受到不同投影方向影響 {cite}`singer2020`。
- **Punjani, A. 與 Fleet, D. J.（2021）．*3D Variability Analysis: Resolving Continuous Flexibility and Discrete Heterogeneity from Single Particle Cryo-EM*. Journal of Structural Biology, 213, 107702。** [DOI](https://doi.org/10.1016/j.jsb.2021.107702)。第 5.2 節的資料模型與最小平方更新，把共用線性子空間放進三維成像問題；可先看模型式，再接續下一篇異質性教學 {cite}`punjani2021`。
