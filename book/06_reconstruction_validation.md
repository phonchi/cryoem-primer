# 3D 重建：從 Fourier 切片到加權最小平方

{doc}`06_dimension_reduction` 把粒子影像壓成少量係數，處理的仍是二維表示。本章跨到三維：一組不同方向的二維投影，如何共同決定一個三維密度？

這個問題有兩個層次。第一層是幾何：即使角度完全已知，投影方向、CTF 與離散取樣仍可能使部分密度無法識別，或讓估計對雜訊很敏感。要分清楚哪些部分能由資料決定，哪些部分需要正則化。第二層是 SPA 特有的：取向本身也是未知量，必須和結構一起估計。先分清這兩層問題，再來討論初始模型偏差、對稱性與解析度。

(fourier-coverage)=
## Fourier slice theorem 是重建的幾何核心

沿用 {doc}`05_image_formation` 的定義，取向 $R$ 下的投影是

$$
P_RV(\mathbf u)=\int V\!\left(R\,[u_x,u_y,z]^{\mathsf T}\right)\mathrm dz .
$$

Fourier slice theorem 指出，這張投影的二維 Fourier transform 等於三維 Fourier transform 在一個通過原點的平面上的取值：

$$
\mathcal F_2\{P_RV\}(\mathbf k)
=\widehat V\!\left(R\,[k_x,k_y,0]^{\mathsf T}\right).
$$

一張投影因此提供 $\widehat V$ 的一個中心切片，切片的法向量就是投影方向。收集方向分散的切片，可以擴大三維 Fourier 空間的覆蓋；缺少切片的方向，則缺乏資料約束。

```{dropdown} Fourier slice theorem 的推導
直接計算二維 Fourier transform：

$$
\mathcal F_2\{P_RV\}(\mathbf k)
=\int\!\!\int P_RV(\mathbf u)\,e^{-2\pi\mathrm i(k_xu_x+k_yu_y)}\,\mathrm du
=\int\!\!\int\!\!\int V\!\left(R[u_x,u_y,z]^{\mathsf T}\right)
e^{-2\pi\mathrm i(k_xu_x+k_yu_y)}\,\mathrm du\,\mathrm dz .
$$

換變數 $\mathbf r=R[u_x,u_y,z]^{\mathsf T}$。因為 $R$ 是旋轉矩陣，$R^{-1}=R^{\mathsf T}$ 且 $|\det R|=1$，所以 $[u_x,u_y,z]^{\mathsf T}=R^{\mathsf T}\mathbf r$，體積元素不變。指數中的相位可以寫成

$$
k_xu_x+k_yu_y
=[k_x,k_y,0]\cdot\left(R^{\mathsf T}\mathbf r\right)
=\left(R\,[k_x,k_y,0]^{\mathsf T}\right)\cdot\mathbf r ,
$$

最後一步用到 $[k_x,k_y,0]R^{\mathsf T}=\left(R[k_x,k_y,0]^{\mathsf T}\right)^{\mathsf T}$。代回去得到

$$
\mathcal F_2\{P_RV\}(\mathbf k)
=\int V(\mathbf r)\,
e^{-2\pi\mathrm i\left(R[k_x,k_y,0]^{\mathsf T}\right)\cdot\mathbf r}\,\mathrm d\mathbf r
=\widehat V\!\left(R[k_x,k_y,0]^{\mathsf T}\right).
$$

積分沿 $z$ 的方向做，對應的頻率分量就是 $k_z=0$，這正是「投影對應中心切片」的來源。注意這裡的 $R$ 與投影定義中的 $R$ 是同一個矩陣，沒有轉置；若把投影定義改成 $V(R^{\mathsf T}[\cdot])$，切片公式中的矩陣也要跟著改。跨軟體比較時，必須把兩式一起核對。
```

切片的分布直接決定重建的品質。取向若集中在少數區域（偏好取向），對應的切片也集中在少數方向，其他方向留下稀疏取樣的空洞，重建圖在該方向會較模糊或被拉長，形成各向異性解析度。增加更多相同取向的粒子，只會提高已取樣方向的訊雜比；缺少的方向仍然空白。正則化能讓數值解穩定，卻無法補出未被量到的資訊。

````{dropdown} 二維物件的中央切線：可直接重算的例子
沿 $y$ 加總二維物件 $f[y,x]$，得到 $p[x]=\sum_y f[y,x]$。重排求和可得

$$\operatorname{FFT}_1(p)[k_x]=\operatorname{FFT}_2(f)[0,k_x].$$

下圖使用 96×128 的 Shepp–Logan 物件，同時計算投影的 FFT 與二維 FFT 的中央切線。兩條曲線重合；三維結構的二維投影則對應一個中央平面。

```{figure} images/spa/slice_example_migrated.png
:alt: 二維 Shepp–Logan 物件，以及投影 FFT 與中央切線相互重合的頻譜曲線。
兩種計算採用相同 FFT 正規化，差異只剩浮點誤差。這裡檢查投影與切線的關係，沒有估計未知角度。
```
````

## 離散反問題：什麼時候能解，什麼時候穩定？

實際資料是有限張、取樣過的投影，待求的是有限個 voxel。把三維密度排成向量 $V$、所有投影排成向量 $y$、成像運算排成矩陣 $A$，問題變成：已知 $y$ 與 $A$，求 $\widehat V$ 使 $A\widehat V\approx y$。

矩陣是否長方形，不能決定反問題是否病態。例如

$$
A=\begin{pmatrix}1&0\\0&1\\1&1\end{pmatrix}
$$

雖然是長方矩陣，卻有兩個非零奇異值 $\sqrt3$ 與 $1$，最小平方解唯一，而且對小擾動相當穩定。真正要看的是矩陣的秩與奇異值：

- 若 $A$ 的秩小於未知係數數量，就存在非零 $V_0$ 滿足 $AV_0=0$。因此 $V$ 與 $V+V_0$ 產生完全相同的量測，資料無法區分它們。
- 若 $A$ 的秩等於未知係數數量，但某些奇異值很小，對應的結構變化只會在影像裡留下微弱痕跡。反演這些方向會放大雜訊，估計因而不穩定。

投影文獻中的 **ghost** 指在所選量測方向下投影為零的物體，也就是觀測算子的 null space 成分 {cite}`penczek2010`。有限方向的連續物體問題可以存在這類歧義；限定到特定的有限 voxel 模型後，是否仍有非零 null space，需要檢查實際矩陣的秩；投影數有限這項條件，本身不足以確定有非零 ghost。

正則化會抑制不穩定方向，並在不可識別的部分選擇符合額外假設的解。增加互補方向或不同 CTF 的影像，有機會增加秩或改善奇異值；只重複相同量測，則不會補足 null space。{doc}`05_statistical_inference` 的 MAP 觀點把額外假設寫成先驗，讓我們能追問每個結構細節是由資料還是先驗決定。

## Backprojection 如何把投影放回三維空間

最直觀的重建方式是把每張投影沿它的方向「塗回」三維空間再相加。這個運算就是 $A$ 的伴隨（adjoint）$A^{\mathsf T}$。它容易計算，但它不是 $A$ 的反矩陣。

差別可以從 Fourier 空間看清楚。每張投影貢獻一個通過原點的切片；把所有切片直接相加時，接近原點的區域被許多切片重複覆蓋，遠離原點的區域則覆蓋稀疏。結果是低頻被過度加權，重建圖看起來像被一個模糊核卷積過。傳統的 filtered backprojection 因此在反投影前先對每張投影套一個隨頻率上升的濾波器，補償這個不均勻的取樣密度。

把同一件事寫成統計問題會更清楚。加權最小平方的解滿足正規方程式

$$
\left(A^{\mathsf T}WA\right)\widehat V=A^{\mathsf T}Wy,
$$

其中 $W$ 是量測權重。右邊是加權反投影；左邊的 $A^{\mathsf T}WA$ 則描述成像與反投影如何共同作用在密度上。一般的內插會耦合不同 voxel，因此這是一個需要求解的矩陣方程。只有採用對角或 binning 近似時，才可將更新理解為逐個 Fourier 位置除以累積權重。

在 Fourier 空間直接做這件事的方法稱為 direct Fourier inversion。它的困難在於切片上的取樣點落在不規則的位置，必須重新內插到規則格點，同時正確地計入不均勻的取樣密度。理論上恰當的權重，是讓加權和近似對取樣分布的積分，也就是各取樣點所佔的 Voronoi 胞體積；實際的 gridding 方法用各種近似來避免直接計算三維 Voronoi 圖的成本 {cite}`penczek2010`。至於代數迭代法（ART、SIRT 一類），它們直接求解線性系統，取樣密度不會顯式出現在公式裡。即使公式裡沒有單獨的取樣密度項，方向覆蓋的缺口仍會影響解 {cite}`penczek2010`。

(reconstruction)=
## CTF 加權合併：一個可以手算的例子

真實重建還要處理 CTF。同一個 Fourier voxel 會被許多張影像的切片經過，而每張影像在該頻率的 $H_i$ 大小與正負號都不同。正確的合併方式是把每個觀測依它的可信度加權。

先用實數係數示範：某個 Fourier 位置有三個觀測，它們的 CTF 值分別是 $H_1=0.9$、$H_2=0.2$、$H_3=-0.6$，雜訊變異數都是 $\sigma^2=1$，而真值 $V=1$。觀測模型是 $y_j=H_jV+\varepsilon_j$，各 $\varepsilon_j$ 零平均且互相獨立。加入高斯先驗（訊號功率 $\tau^2$）後的 MAP 估計為

$$
\widehat V=\frac{\sum_j H_j\,y_j}{\sum_j H_j^2+\sigma^2/\tau^2}.
$$

先算分母的資料項：$\sum_j H_j^2=0.81+0.04+0.36=1.21$。

**不加先驗時**（$\sigma^2/\tau^2=0$），估計是無偏的：$\mathbb E[\widehat V]=1$，但每次含雜訊的估計值仍會波動。它的變異數為 $\sigma^2/\sum H_j^2=1/1.21\approx0.826$，這就是均方誤差。

**取 $\tau^2=1$ 時**，分母變成 $1.21+1=2.21$。重複收集含雜訊的資料時，估計的期望為 $\mathbb E[\widehat V]=1.21/2.21\approx0.548$，比真值更接近零。它的變異數是 $\sigma^2\sum H_j^2/(\sum H_j^2+\sigma^2/\tau^2)^2=1.21/4.8841\approx0.248$，偏差平方為 $(0.548-1)^2\approx0.205$，兩者相加得到均方誤差 $\approx0.452$。

收縮讓估計變得有偏，均方誤差卻從 0.826 降到 0.452。這是{doc}`05_statistical_inference`的偏差－變異取捨在三維重建裡的具體樣子。

先驗的價值在資料稀少時更明顯。若某個頻率只被一張接近 CTF 零點的影像涵蓋，$H=0.02$，則無先驗估計的變異數是 $1/0.0004=2500$，估計容易被雜訊主導。加入同樣的先驗後，分母成為 $0.0004+1\approx1$，估計值被壓到接近零。這個頻率的資料精度遠低於先驗精度，因此估計更接近先驗的零平均。

(soft-assignment)=
## 角度與位移未知時：邊際化的 refinement

以上都假設取向已知。醫學斷層掃描確實如此，投影角度由儀器控制；SPA 沒有這個條件。

教學上可把 refinement 寫成三步：由目前模型產生候選投影；為每張粒子影像找最相似的方向與位移；以估計的角度與位移重建新模型，再重複。這是 hard projection matching。它建立直覺，但在單張影像訊雜比很低時，最佳與次佳角度、位移組合的分數可能幾乎相同，只保留最高分等於把一個近乎隨機的選擇當成事實，並在下一輪迭代中放大它。

機率加權的做法則對候選組合做邊際化。以 $z$ 表示角度、位移與可能的類別，每張影像對每個候選的後驗權重是

$$
r_i(z)=\frac{p(y_i\mid z,V)\,p(z)}{\sum_{z'}p(y_i\mid z',V)\,p(z')},
\qquad \sum_z r_i(z)=1 .
$$

角度與位移是連續變數時，需要積分；以離散網格近似時，加權和要包含網格的體積元素。更新模型時使用 $r_i(z)$ 加權所有候選，讓每張影像保留多組可能的角度與位移 {cite}`scheres2010,sigworth1998`。

把候選權重加入成像模型，可以寫出精確的 MAP 正規方程式。令 $A_i(z)$ 包含候選角度的投影、CTF，以及位移對應的 Fourier 相位；$\Sigma_i$ 是雜訊共變異數，$\Lambda$ 是零平均 Gaussian 結構先驗的精度矩陣。固定這一輪的權重及精度後，更新滿足

$$
\left[\sum_{i,z}r_i^{(n)}(z)A_i(z)^*\Sigma_i^{-1}A_i(z)+\Lambda\right]V^{(n+1)}
=\sum_{i,z}r_i^{(n)}(z)A_i(z)^*\Sigma_i^{-1}y_i .
$$

星號表示共軛轉置；在實數座標下就是轉置。右邊把每張影像按候選角度、位移、CTF 與雜訊精度反投影；左邊同時納入成像算子的耦合與先驗。若模型含多個構形，$V$ 可堆疊各構形的密度，$A_i(z)$ 依候選類別選取對應部分。

上一節的單係數分式，是這個方程在各係數可分開估計時的簡化。一般的 gridding 與內插會耦合鄰近係數；使用逐 voxel 的「資料分子除以資料精度加先驗精度」時，必須說明採用的對角或 binning 近似，也不能漏掉候選角度的切片位置與位移相位 {cite}`scheres2012bayes`。

Scheres 2012 的 Bayesian 實作把 CTF 校正、取樣權重與先驗收縮放進同一估計過程，並從資料逐輪更新頻率相依的雜訊與訊號功率。不過，模型仍需先驗形式、角度取樣與參數設定；原文也使用經驗倍率 $T=4$ 調整先驗功率，不能把它理解為完全不需人工設定的方法 {cite}`scheres2012bayes`。

從重建本身更新先驗也有風險：若某輪把雜訊當成訊號，下一輪就可能保留更多雜訊。{doc}`06_resolution_validation` 會說明獨立 half-set 如何協助檢查這個問題。

## 初始模型、ab initio 與參考偏差

上述迭代屬於局部最佳化；不同初始模型可能導向不同的局部解，不能保證找到全域最佳解，也不保證解在起點附近。初始參考中的特徵可能因此保留在最終結構裡，形成 model bias {cite}`sigworth2016`。資料品質高時，即使從平滑的橢球或隨機密度出發也常能收斂到正確結構；訊雜比低或偏好取向嚴重時，偏差可能很嚴重 {cite}`sigworth2016`。

從資料本身建立初始模型的做法稱為 ab initio。常見路線包括利用 common lines——兩張不同投影方向的中心切片，在三維 Fourier 空間相交於一條直線，這條共線關係可用來同時決定一組取向 {cite}`penczek2010`；以及隨機初始化加上隨機性優化。Punjani 等人 2017 的 cryoSPARC 方法採用後者：把整批資料的邊際概似目標寫下來，每次迭代只用一個隨機抽取的影像子集合估計梯度，用這些近似梯度更新結構，演算法可從隨機初始化開始，而且能在不假設各狀態彼此相似的情況下做 ab initio 的三維分類 {cite}`punjani2017`。

無論用哪條路線，參考偏差的檢查方式是一致的：換初始模型、換隨機種子重跑，看主要密度特徵是否重現。這和 {doc}`06_workflow` 提到的模板挑選偏差是同一類問題——參考裡的特徵會透過「挑出最像參考的東西」這個動作回到結果裡。

## 對稱假設與手性歧義

**施加對稱**是常見的加速手段。若結構具有點群對稱 $G$，每張投影可約束 $|G|$ 個對稱等價的切片方向，並減少待估的自由度。這會改善對稱模型下的覆蓋與估計精度；不過，由同一張影像複製出的對稱資訊共享雜訊，不能當成新增的 $|G|$ 倍獨立觀測。代價是：若真實結構並不具有該對稱，強制對稱等於把不對稱的特徵平均掉。這種錯誤在兩個 half-set 上會以完全相同的方式發生，因此半數資料之間的一致性檢查看不到它。

**手性**是更根本的歧義。把三維密度鏡射，同時對所有取向做對應的變換，得到的投影集合與原來完全相同。

```{dropdown} 為什麼投影資料本身無法決定手性
設 $M$ 是一個鏡射矩陣（$\det M=-1$，且 $M^{-1}=M$），定義鏡射後的體 $V'(\mathbf r)=V(M\mathbf r)$。再取 $S=\mathrm{diag}(1,1,-1)$，並令

$$
R'=M\,R\,S .
$$

因為 $\det R'=(-1)(1)(-1)=1$，$R'$ 仍是合法的旋轉。代入投影定義：

$$
P_{R'}V'(\mathbf u)
=\int V\!\left(M R'[u_x,u_y,z]^{\mathsf T}\right)\mathrm dz
=\int V\!\left(R\,S\,[u_x,u_y,z]^{\mathsf T}\right)\mathrm dz
=\int V\!\left(R\,[u_x,u_y,-z]^{\mathsf T}\right)\mathrm dz .
$$

把積分變數換成 $z'=-z$，右邊就是 $P_RV(\mathbf u)$。也就是說，對每一個取向 $R$，鏡射後的體都能在取向 $R'$ 下產生**完全相同**的投影。兩個互為鏡像的結構因此對投影資料給出相同的概似，無論收集多少影像都無法分辨。

要定出手性，需要資料以外的資訊：已知的傾斜角度關係、可辨認的右手螺旋二級結構，或與已知結構的比對。實務上也要留意，軟體慣例（Euler 角順序、矩陣作用方向、影像的軸向）不一致時，會讓輸出的密度整體鏡射。
```

## 兩套實作各自最佳化什麼

以下比較 Scheres 2012 與 Punjani 等人 2017 兩篇論文描述的實作，分別看統計模型與最佳化策略。這些歷史版本不能代表兩套軟體目前的全部功能。

**統計模型層**：兩者寫下的目標函數形式非常接近。RELION 的出發點是 MAP 估計——在 Fourier 空間對訊號加上零平均高斯先驗，讓平滑性成為明確的統計假設，再用 EM 同時更新結構、雜訊功率與先驗功率 {cite}`scheres2012bayes`。cryoSPARC 把結構決定明白寫成一個 Bayesian 概似的最佳化問題，對每張影像的角度與位移積分、對多個結構求和，並加上模型的先驗項 {cite}`punjani2017`。兩者都以含隱藏變數的機率模型估計結構。

**最佳化策略**：RELION 的 EM 每一輪使用全部影像計算期望步驟，更新穩定但每輪昂貴，而且作為局部最佳化，它對初始模型敏感 {cite}`scheres2012bayes`。cryoSPARC 用隨機梯度下降，每次迭代只取一個隨機影像子集合來近似目標函數；該文以隨機初始化配合這些更新做 ab initio 重建，但仍需檢查初始化敏感性。高解析度階段則改用 branch-and-bound 的對位搜尋，藉由排除不可能的候選來減少重複計算 {cite}`punjani2017`。

**該怎麼讀論文裡的數字**：兩篇論文報告的速度與解析度，是在各自當時的實作、資料集與硬體條件下取得的。把它們讀成「某某方法一律較好」會超出證據範圍；要比較，必須在同一批資料、同一組前處理與同一套評估標準下進行，並把隨機種子與參數設定一併記錄。

```{figure} images/pptx/s12_8.png
:width: 35%
:name: fig-reconstructed-map
由大量投影重建出的三維密度圖。密度的外觀、解析度估計與生物詮釋是三個不同層次的問題，需要分開評估。
```

概似值趨於穩定、相鄰兩輪的 map 很相似，或角度與位移不再大幅改變時，演算法通常會停止迭代。這些現象表示更新幅度已經很小，也可能發生在局部最佳解附近。接下來要換初始模型或合理的參數重跑，再比較方向分布、密度特徵與獨立重建是否一致。

(validation)=
## 驗證的入口：從一致性到解析度

重建結束後最常被問的問題是「這張圖的解析度是多少」。SPA 回答這個問題的方式，是檢查兩份獨立處理的資料能否得到一致的結果。

標準做法是在 refinement 之前把粒子隨機分成兩個 half-set，兩組從此獨立估計角度與位移並各自重建，得到 half-map A 與 half-map B。取包含正負共軛頻率的完整對稱 shell，半徑記為 $s$。Fourier shell correlation 定義為

$$
\mathrm{FSC}(s)=
\frac{\sum_{\mathbf k\in s}F_A(\mathbf k)F_B(\mathbf k)^*}
{\sqrt{\sum_{\mathbf k\in s}|F_A(\mathbf k)|^2
\sum_{\mathbf k\in s}|F_B(\mathbf k)|^2}} .
$$

FSC 衡量兩張 half-map 的一致性：數值為 1 表示兩組 Fourier 係數完全線性一致；若其中一組只是另一組乘上共同的正比例常數，FSC 仍為 1，所以它不要求振幅數值完全相同。接近 0 則表示該頻帶缺少線性一致性。

Mask 可以排除大面積溶劑、提高數值穩定性，但也會混合頻率；共同的低頻訊號可能混入高頻，而由完整資料細節製作的 mask 還可能引入共同偏差。事先固定的同一個 mask 並不會自行讓原本獨立的兩組雜訊產生交叉相關。同樣地，若兩邊共用了利用高解析度資料估得的角度與位移資訊或來自同一部 movie 的重複觀測，half-set 的獨立性就被破壞。

這些議題——FSC 與頻譜訊雜比的關係、0.143 這個慣用門檻從哪裡來、mask 造成的頻率混合如何校正、局部與方向解析度、銳化與 B factor，以及兩份資料共有的偏差為什麼檢查不出來——構成完整的一章：{doc}`06_resolution_validation`。

## 接下來

本章假設所有粒子來自同一個結構。實際樣品常含有多種構形或組成狀態，而把它們混在一起重建，會讓變動的部分在平均中模糊掉。下一章處理這件事：{doc}`06_heterogeneity`。

## 延伸閱讀

- **Penczek, P. A.（2010）．*Fundamentals of Three-Dimensional Reconstruction from Projections*. Methods in Enzymology, 482, 1–33。** [DOI](https://doi.org/10.1016/S0076-6879(10)82001-4)。第 1 節（頁 2–5）以矩陣形式說明反演的病態性與 ghost 的存在，是本章離散反問題一節的來源；第 2 節（頁 5–8）給出 central section theorem；第 5 節（頁 13–18）是代數迭代法，第 6 節（頁 18–22）是 filtered backprojection，第 7 節（頁 22–24）說明 gridding 與以 Voronoi 胞體積作為取樣權重的理由 {cite}`penczek2010`。
- **Scheres, S. H. W.（2012）．*A Bayesian View on Cryo-EM Structure Determination*. Journal of Molecular Biology, 415, 406–418。** [DOI](https://doi.org/10.1016/j.jmb.2011.11.010)。頁 408–411 是核心：先用 Wiener filter 的常見寫法說明 CTF 加權合併與正則化常數的角色，再導出以 Gaussian prior 為基礎的 MAP 更新式，並說明雜訊功率與訊號功率如何由 EM 逐輪從資料估出。本章的更新式結構即取自此處 {cite}`scheres2012bayes`。
- **Punjani, A., Rubinstein, J. L., Fleet, D. J. 與 Brubaker, M. A.（2017）．*cryoSPARC: Algorithms for Rapid Unsupervised Cryo-EM Structure Determination*. Nature Methods, 14(3), 290–296。** [DOI](https://doi.org/10.1038/nmeth.4169)。頁 290–291 寫下與上一篇同一族的 Bayesian 目標函數，並以 Figure 1 說明該演算法如何從隨機初始化出發；Supplementary Note 1 詳述 SGD，branch-and-bound 對位見主文相應方法說明。閱讀時把「統計模型」與「最佳化策略」分開記 {cite}`punjani2017`。
- **Sigworth, F. J., Doerschuk, P. C., Carazo, J.-M. 與 Scheres, S. H. W.（2010）．*An Introduction to Maximum-Likelihood Methods in Cryo-EM*. Methods in Enzymology, 482, 263–294。** [DOI](https://doi.org/10.1016/S0076-6879(10)82011-7)。把角度與位移的邊際化及 EM 更新從二維推到三維，補上本章 $r_i(z)$ 一節省略的推導步驟 {cite}`scheres2010`。
- **Sigworth, F. J.（2016）．*Principles of Cryo-EM Single-Particle Image Processing*. Microscopy, 65(1), 57–67。** [Oxford Academic 免費全文](https://doi.org/10.1093/jmicro/dfv370)。頁 61–62 集中討論 model bias：它如何隨訊雜比與偏好取向惡化，以及 ab initio 模型、傾斜對分析等檢查方式。公開取用，適合作為本章參考偏差一節的對照 {cite}`sigworth2016`。
- **Singer, A. 與 Sigworth, F. J.（2020）．*Computational Methods for Single-Particle Electron Cryomicroscopy*. Annual Review of Biomedical Data Science, 3, 163–190。** [DOI](https://doi.org/10.1146/annurev-biodatasci-021020-093826)。第 4 節（頁 174–182）依序處理最大概似、EM、Bayesian 推論與隨機梯度下降，可與 {doc}`05_statistical_inference` 的概似及 EM 公式對照 {cite}`singer2020`。
- **[ASPIRE-Python](https://github.com/ComputationalCryoEM/ASPIRE-Python) 的 `aspire.reconstruction` 模組。** 開放原始碼實作了本章描述的加權最小平方重建與正則化項，可直接對照公式與程式中的權重、內插與正規化步驟；{doc}`07_synthetic_data` 產生的模擬資料可直接餵給它。
