# SPA 的統計推論：從含雜訊影像到 EM

上一章的成像模型告訴我們：給定結構、取向、位移與 CTF，可以預測一張粒子影像。分析實驗資料時，問題的方向倒過來了。我們看到的是影像，要由它推估結構；偏偏取向與位移也未知。一張模糊影像常同時支持好幾種解釋，統計推論的工作就是比較這些解釋，並把不確定性帶進結構估計。

本章從基礎統計的常態分布與加權平均出發，先建立概似，再用一個可以手算的混合模型學習 EM，最後把同一套思路接回 SPA。後面的對齊、分類與重建會反覆使用這些概念。

(spa-likelihood)=
## 觀測、未知量與假設要先分清楚

令 $y_i\in\mathbb R^d$ 表示第 $i$ 張粒子影像排列成的向量，$d$ 是像素數。用 $z_i$ 表示這張影像未知的角度、位移與可能的類別，用 $\theta$ 表示整批資料共享的參數，例如三維結構、類別比例與雜訊變異數。模型寫成

$$
y_i=m_i(z_i;\theta)+\varepsilon_i.
$$

$m_i$ 是把候選結構經過投影、CTF 與平移後得到的預測影像。CTF 已由前處理估計時，可先把每張影像的 CTF 當作已知參數；此時模型沒有另外描述 CTF 的估計誤差。把哪些量固定、哪些量一起估計，會影響最後的不確定性。

先採用 $\varepsilon_i\sim N(0,\sigma^2I_d)$，並假設不同粒子的雜訊獨立。給定一個候選 $z$，觀測影像的機率密度為

$$
p(y_i\mid z,\theta)
=(2\pi\sigma^2)^{-d/2}
\exp\!\left[-\frac{\|y_i-m_i(z;\theta)\|^2}{2\sigma^2}\right].
$$

固定已觀測的 $y_i$，把右側視為候選參數的函數，就是**概似（likelihood）**。它回答「若使用這組參數，觀測到這張影像有多合理」。連續資料的單點機率為零，因此這裡比較的是密度與密度比。

在相同 $\sigma^2$ 下，較小的平方殘差有較高概似。這說明最小平方法為何會自然出現在影像匹配；若不同位置的雜訊精度不同，則需要改變每個殘差的權重 {cite}`scheres2010,sigworth1998`。

## 雜訊有相關性時，距離也要改變

插值、濾波與共同背景可能讓鄰近像素的誤差相關。若 $\varepsilon_i\sim N(0,\Sigma_i)$，負對數概似包含

$$
\frac12\bigl[y_i-m_i(z;\theta)\bigr]^{\mathsf T}
\Sigma_i^{-1}\bigl[y_i-m_i(z;\theta)\bigr]
+\frac12\log\det\Sigma_i.
$$

第一項是以雜訊共變異數調整的距離。誤差本來就大的方向，獲得的權重較低；高度相關的像素，也不宜全都當成新的獨立證據。當 $\Sigma_i$ 已知且固定時，第二項不影響同一張影像在不同角度與位移下的比較；估計雜訊參數時，它仍是概似的一部分。

例如兩個像素的雜訊變異數分別為 1 與 9，且互相獨立。候選 A 的殘差為 $(2,0)$，候選 B 為 $(0,3)$。普通平方距離為 4 與 9，偏向 A；依雜訊變異數加權後，距離成為 4 與 1，轉而偏向 B。第二個像素偏離 3，只是一個標準差；第一個像素偏離 2，已經是兩個標準差。

**Whitening** 可用 $W_i\Sigma_iW_i^{\mathsf T}=I$ 的轉換，把這個加權距離改寫為 $\|W_i(y_i-m_i)\|^2$。觀測與預測需要一起轉換。SPA 常在 Fourier 空間用雜訊功率頻譜近似各頻率的變異數；這種近似依賴平穩性等假設；背景不平穩或視窗有限時，近似可能不夠準確。

(spa-posterior)=
## 從兩組候選角度與位移看 Bayes 更新

概似描述影像和候選的符合程度；先驗 $p(z\mid\theta)$ 描述看到這張影像以前對候選的權重。Bayes 公式把兩者相乘再正規化：

$$
r_i(z)=p(z\mid y_i,\theta)
=\frac{p(y_i\mid z,\theta)p(z\mid\theta)}
{\sum_{z'}p(y_i\mid z',\theta)p(z'\mid\theta)}.
$$

$r_i(z)$ 是後驗權重，也稱 responsibility。分母收集所有候選的證據，因此每張影像的權重總和為 1。角度與位移是連續變數時，需要積分；數值網格的加權和應包含網格體積，等權候選只適用於相應的離散先驗或等測度網格。

考慮兩組先驗相等的候選角度與位移 A、B。若它們的白雜訊平方殘差分別為 $E_A$、$E_B$，則

$$
\frac{r_i(A)}{r_i(B)}
=\exp\!\left[\frac{E_B-E_A}{2\sigma^2}\right].
$$

取 $E_A=10$、$E_B=12$。當 $\sigma^2=1$，A 的權重約為 0.731；當 $\sigma^2=4$，A 的權重約為 0.562。同樣的殘差差距，在較大的雜訊下只提供微弱的區別。只選分數最高的角度與位移組合時，兩種情況都會留下 A；保留機率權重，則能讓重建更新反映兩者不同的把握程度。

```{figure} images/spa/em_uncertainty.png
:name: fig-spa-posterior
:alt: 相同的候選殘差差距，在較大雜訊變異數下產生較接近二分之一的後驗權重。
:width: 88%

左圖比較兩個雜訊水準下，候選角度與位移的權重；右圖畫出殘差差距與後驗機率的關係。曲線越平緩，越需要較大的分數差距才能明確區分候選。
```

這些權重是**給定目前結構與假設後**的條件機率。如果候選結構錯誤、雜訊變異數低估，或真實角度與位移根本沒有放進候選集合，分布仍可能很集中。判讀後驗時，也要檢查產生它的模型。

(spa-marginalization)=
## 未知角度與位移如何進入整批資料的概似？

同一張影像可能對應多組角度與位移。因此觀測機率要把這些互斥的可能性加起來：

$$
p(y_i\mid\theta)=\sum_z p(y_i,z\mid\theta),\qquad
\ell(\theta)=\sum_{i=1}^{n}\log\!\left[\sum_z p(y_i,z\mid\theta)\right].
$$

這個動作稱為**邊際化（marginalization）**。它保留多個候選對同一影像的解釋能力。取最大值只採用最佳候選；求和則把其他合理候選的證據也納入。當後驗非常集中，兩種處理可能接近；低訊雜比下，差異往往更明顯。

困難在於對數裡面還有求和，候選角度與位移的權重又隨結構改變。若角度與位移已知，重建常能整理成加權最小平方法；這些量未知時，就要和結構一起估計。EM 就利用「如果隱藏資料已知，問題比較容易」這個特性，反覆更新權重與模型。

(spa-em)=
## 先用兩群一維資料做完一次 EM

暫時把影像縮成一個數 $x_i$，以兩個常態成分描述資料：

$$
p(x_i)=\pi_1\phi(x_i;\mu_1,\sigma^2)+\pi_2\phi(x_i;\mu_2,\sigma^2),
\qquad \pi_1+\pi_2=1.
$$

$\phi$ 是常態密度，隱藏類別為 $c_i\in\{1,2\}$。為了看清一輪更新，固定 $\sigma^2=1$，使用 $x=(0,1,3,4)$，初始平均數 $(\mu_1,\mu_2)=(0,4)$，初始比例各為一半。

**E-step：先計算每筆資料的類別權重。**

$$
r_{ik}^{(t)}=
\frac{\pi_k^{(t)}\phi(x_i;\mu_k^{(t)},1)}
{\sum_{h=1}^{2}\pi_h^{(t)}\phi(x_i;\mu_h^{(t)},1)}.
$$

| $x_i$ | 第 1 群權重 | 第 2 群權重 |
|---:|---:|---:|
| 0 | 0.999665 | 0.000335 |
| 1 | 0.982014 | 0.017986 |
| 3 | 0.017986 | 0.982014 |
| 4 | 0.000335 | 0.999665 |

靠近 0 的資料主要支持第一群，靠近 4 的資料主要支持第二群。每筆資料都保留一點分到另一群的可能性。

**M-step：固定剛才的權重，更新模型。**

$$
N_k=\sum_i r_{ik}^{(t)},\qquad
\pi_k^{(t+1)}=N_k/n,\qquad
\mu_k^{(t+1)}=\frac{\sum_i r_{ik}^{(t)}x_i}{N_k}.
$$

此例 $N_1=N_2=2$，新的平均數約為 $(0.519, 3.481)$。第一群的中心向 0 與 1 的中間移動，第二群則向 3 與 4 的中間移動。接著用新中心重算權重，再更新平均數，直到更新變小。這裡固定變異數；若把變異數也列為未知量，M-step 還要更新相應的加權殘差平方平均。

## 把這兩步換回粒子影像

對一般參數 $\theta$，EM 的 M-step 最大化

$$
Q(\theta\mid\theta^{(t)})=
\sum_i\sum_z r_i^{(t)}(z)\log p(y_i,z\mid\theta).
$$

先用舊參數得到 $r_i^{(t)}$，再把它們當作固定權重更新 $\theta$。E-step 的 E 是對隱藏變數取條件期望；M-step 的 M 是最大化這個期望後的完整資料對數概似。

在二維分類中，$z$ 可以同時包括類別、平面旋轉與平移，$\theta$ 包括各類的平均影像。在三維重建中，$z$ 包括三維取向與平移，$m_i(z;\theta)$ 是候選結構經 CTF 調變的投影。每張粒子的多組可能角度與位移，都依權重參與更新 {cite}`sigworth1998,scheres2005,scheres2010`。

白高斯雜訊且變異數固定時，模型更新的主要部分可寫成

$$
\min_V\sum_{i,z}\frac{r_i^{(t)}(z)}{2\sigma_i^2}
\|y_i-A_i(z)V\|^2.
$$

$A_i(z)$ 包含投影、CTF 與位移；在本式中權重與操作算子都已由 E-step 固定。後續的{doc}`06_reconstruction_validation`會說明如何利用這個加權反問題更新結構。

```{dropdown} EM 為何能使概似不下降？

為每張影像選一個機率分布 $q_i(z)$，定義

$$
\mathcal L(q,\theta)=\sum_{i,z}q_i(z)
\log\frac{p(y_i,z\mid\theta)}{q_i(z)}.
$$

由 Bayes 公式可以整理出

$$
\ell(\theta)=\mathcal L(q,\theta)+
\sum_i D_{\mathrm{KL}}\!\left(q_i\,\|\,p(z\mid y_i,\theta)\right).
$$

KL divergence 非負，所以 $\mathcal L$ 是對數概似的下界。E-step 選 $q_i=p(z\mid y_i,\theta^{(t)})$，使下界在舊參數處等於概似。M-step 固定 $q$ 並提高下界，於是

$$
\ell(\theta^{(t+1)})\geq\mathcal L(q,\theta^{(t+1)})
\geq\mathcal L(q,\theta^{(t)})=\ell(\theta^{(t)}).
$$

因此，只要 E-step 精確，M-step 又提高同一下界，概似就不會下降。若要再推出概似值收斂，需要上界；若要討論參數收斂或全域最佳解，還需要額外條件。自由變異數的 Gaussian mixture 甚至可能讓某個成分縮到單一資料點，造成概似無界。網格截斷、近似搜尋與不同更新規則，也要另外檢查其目標函數 {cite}`ma2019em`。
```

(spa-map)=
## 資料不足時，先驗如何穩定估計？

CTF 零點或少見取向使某些結構方向幾乎沒有資料支持。最大概似估計只按資料符合程度選模型；MAP 估計再加入結構先驗：

$$
\widehat V_{\mathrm{MAP}}
=\arg\max_V\{\log p(y\mid V)+\log p(V)\}.
$$

先看一個實數係數的例子。令 $y=h x+\varepsilon$、$\varepsilon\sim N(0,\sigma^2)$，並給 $x\sim N(0,\tau^2)$ 的先驗。MAP 目標等價於最小化

$$
\frac{(y-hx)^2}{2\sigma^2}+\frac{x^2}{2\tau^2},
$$

其解為

$$
\widehat x=\frac{h}{h^2+\sigma^2/\tau^2}\,y.
$$

這正是 Wiener 型收縮。當 $h$ 很小，直接除以 $h$ 會放大雜訊；先驗讓估計朝零收縮，反映這個係數缺乏充分證據。取 $h=0.1$、$\sigma^2/\tau^2=0.04$，逆濾波增益為 10，MAP 增益為 2。訊號與雜訊在輸出中都受到這個增益影響。

在 Gaussian 線性例子中，後驗平均與後驗眾數相同，因此 MAP 也給出平方損失下的 Bayes 估計；更一般的分布需要分開看待兩者。RELION 的 Bayesian 觀點把這種資料與先驗的平衡帶入頻率相依的三維結構估計 {cite}`scheres2012bayes`。

```{dropdown} 從求導看收縮解與後驗變異數

對 $x$ 微分並令結果為零：

$$
-\frac{h(y-hx)}{\sigma^2}+\frac{x}{\tau^2}=0,
\qquad
\left(\frac{h^2}{\sigma^2}+\frac1{\tau^2}\right)x=\frac{hy}{\sigma^2}.
$$

因此後驗變異數為 $(h^2/\sigma^2+1/\tau^2)^{-1}$。資料精度與先驗精度相加；當 $h=0$，觀測沒有提供這個係數的資訊，後驗回到先驗。可見估計值即使穩定，也可能主要由先驗決定。
```

## 平滑帶來的取捨，如何用偏差與變異描述？

固定真實訊號 $x$，想像重新收集很多批資料，每批都得到估計 $\widehat x$。平方誤差的期望為

$$
\mathbb E[(\widehat x-x)^2]
=\bigl(\mathbb E[\widehat x]-x\bigr)^2+
\operatorname{Var}(\widehat x).
$$

第一項描述平均估計偏離真值的程度；第二項描述不同資料批次之間的波動。收縮通常能降低變異，但也可能把真實的小特徵一起壓低。對 $y=x+\varepsilon$ 使用 $\widehat x=a y$ 時，MSE 為 $(a-1)^2x^2+a^2\sigma^2$，直接看得出兩項如何隨 $a$ 變化。

若 $x=1$、$\sigma^2=4$，不收縮的 $a=1$ 給出 MSE 4；取 $a=0.2$，MSE 為 $0.64+0.16=0.8$。這個改善取決於指定的訊號與雜訊條件，不能由訓練影像看起來更平滑就直接推得。真實資料沒有完整真值，後面的{doc}`06_resolution_validation`會用獨立 half-set、模型假設與共同偏差進一步討論。

有了概似、權重與先驗，下一章可以沿著{doc}`06_workflow`追蹤每一步究竟在估計什麼，以及前處理如何改變後續統計分析的輸入。

## 延伸閱讀

- **Ma, T. 與 Ng, A.（2019）．*The EM algorithm*, CS229 lecture notes, Part IX。** [Stanford 開放講義](https://cs229.stanford.edu/notes2020spring/cs229-notes8.pdf)。本地獨立講義第 2 節、頁 2–8 從邊際概似講到下界與單調性，適合在做完本章一輪數值更新後閱讀 {cite}`ma2019em`。
- **Sigworth, F. J.（1998）．*A maximum-likelihood approach to single-particle image refinement*. Journal of Structural Biology, 122, 328–339。** [DOI](https://doi.org/10.1006/jsbi.1998.4014)。閱讀其影像模型與 maximum-likelihood 推導，看二維旋轉、平移如何成為隱藏變數 {cite}`sigworth1998`。
- **Sigworth, F. J., Doerschuk, P. C., Carazo, J.-M. 與 Scheres, S. H. W.（2010）．*An Introduction to Maximum-Likelihood Methods in Cryo-EM*. Methods in Enzymology, 482, 263–294。** 從影像概似讀到 EM 更新，再看分類與重建的應用；完整書目見引用 {cite}`scheres2010`，可接著使用 [Singer 與 Sigworth 的開放綜述](https://pmc.ncbi.nlm.nih.gov/articles/PMC8412055/)對照三維情境。
- **Scheres, S. H. W.（2012）．*A Bayesian View on Cryo-EM Structure Determination*. Journal of Molecular Biology, 415, 406–418。** [DOI](https://doi.org/10.1016/j.jmb.2011.11.010)。重點讀 Gaussian prior 與 MAP 更新的部分，對照本章單一係數的收縮解 {cite}`scheres2012bayes`。
- **Singer, A. 與 Sigworth, F. J.（2020）．*Computational Methods for Single-Particle Electron Cryomicroscopy*. Annual Review of Biomedical Data Science, 3, 163–190。** [開放全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8412055/)；第 4 節 Maximum Likelihood 一併說明角度與位移的邊際化及計算代價，適合銜接完整 SPA 模型 {cite}`singer2020`。
