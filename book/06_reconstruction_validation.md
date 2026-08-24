# 3D 重建與驗證：從 Fourier 切片到 half-map FSC

```{admonition} 讀完本章，你應該能
:class: important
- 用 Fourier slice theorem 解釋 2D 投影如何約束 3D 結構。
- 分辨 hard projection matching 與機率加權的 pose estimation。
- 說明 gold-standard half-set、FSC、mask 與 overfitting 之間的關係。
- 分清「演算法收斂」「FSC 通過閾值」與「結構正確」各自能說明什麼。
```

(fourier-coverage)=
## Fourier slice theorem 是重建的幾何核心

令 $P_RV$ 為三維結構 $V$ 在取向 $R$ 下的投影。Fourier slice theorem 表示

$$
\mathcal F_2\{P_RV\}(\mathbf k)
=\widehat V\!\left(R^{\mathsf T}[k_x,k_y,0]^{\mathsf T}\right).
$$

一張投影的 2D Fourier transform，對應到 3D Fourier volume 中通過原點的一個平面。若取向已知且覆蓋足夠，可把許多切片插值到 3D 格點，再做 inverse FFT；也可在實空間用 weighted backprojection，或以 ART／SIRT 等迭代法解離散線性系統 {cite}`penczek2010`。

有限方向、非均勻取樣、CTF 零點與插值都讓反演不適定。偏好取向（preferred orientation）會使 Fourier 空間某些方向取樣較弱，造成各向異性解析度。方向分布圖若集中在少數區域，重建圖上常可看到某些方向較模糊或被拉長。正則化能穩定數值解，但無法補出未取樣方向的資訊。

(soft-assignment)=
## SPA 多了一個未知量：姿態

醫學 CT 的投影角度由儀器控制；SPA 的粒子取向與平面位移通常未知。教學上可把 refinement 寫成三步：

1. 由目前 3D 模型產生候選投影。
2. 為每張粒子影像找最相似的方向與位移。
3. 以估計姿態重建新模型，再重複。

這是 hard projection matching 的近似。當單張影像的訊雜比很低，最佳與次佳姿態可能幾乎同樣合理；只保留最高相關者會丟掉不確定性，也可能放大初始模型偏差。

Maximum-likelihood 方法改為對姿態、位移與類別做邊際化。若 $z$ 表示這些 latent variables，一張影像的 responsibility 可概念化為

$$
r_i(z)\propto p(I_i\mid z,V)\,p(z),
\qquad \sum_z r_i(z)=1.
$$

更新模型時使用 $r_i(z)$ 加權所有候選，讓影像保留多個可能姿態。Bayesian／MAP 方法再加入結構先驗或頻率相依正則化。後驗權重若集中在單一候選附近，表示目前模型下的姿態較明確；權重分散時，影像仍支持多個方向或位移。權重形狀也會隨概似函數、雜訊模型、先驗分布、離散網格與初始模型而改變 {cite}`scheres2010,singer2020`。

(reconstruction)=
## 從加權切片到 3D volume

實際重建會同時處理 CTF、取樣密度與正則化。對同一 Fourier voxel，來自不同粒子的觀測可帶有不同 CTF 正負號與振幅；合併時需保留這些權重，並在所有影像都缺乏資訊的頻率避免不穩定除法。不同軟體可能採用 gridding、backprojection、iterative least squares 或數種方法的組合，實作內容遠超過單一步驟的反投影。

概似值趨於穩定、相鄰兩輪 map 很相似，或方向不再大幅改變時，演算法通常會停止迭代。這些現象表示更新幅度已經很小，也可能發生在局部最佳解附近。接下來可更換初始模型或合理參數重跑，再比較方向分布、half-map 與密度特徵是否一致。

(validation)=
## Gold-standard half-set：把獨立性放進流程

常見做法是在 refinement 前把粒子隨機分成兩半，兩組從此獨立 refinement，得到 half-map A 與 half-map B。兩邊可共享事先確定的處理規則；高頻 map、姿態，以及由完整資料高頻細節製作的 mask 必須維持分離，否則 half-map correlation 會被人為提高。

實際資料流可寫成：**particle stack → 分成兩組 particle IDs → 各自估計姿態並重建 → half-map A／B → FSC**。若兩張 half-map 的低頻輪廓相近，但其中一張在局部顯得較模糊，可先比較兩組的粒子數、方向分布與局部解析度。

「Gold-standard」指的是用資料分割降低高頻過擬合的原則。若同一個物理粒子、同一段 movie 的重複觀測或高度相關資料被分到兩邊，兩組資料的獨立性就會降低。方向分布、局部解析度與生化證據則從其他角度檢查重建結果。

## FSC 衡量兩張 half-map 的一致性

對半徑為 $s$ 的 Fourier shell，half-map FSC 常寫成

$$
\mathrm{FSC}(s)=
\frac{\sum_{\mathbf k\in s}F_A(\mathbf k)F_B(\mathbf k)^*}
{\sqrt{\sum_{\mathbf k\in s}|F_A(\mathbf k)|^2
\sum_{\mathbf k\in s}|F_B(\mathbf k)|^2}}.
$$

FSC = 1 表示兩組 Fourier 係數完全線性一致；若其中一組只是另一組乘上共同的正比例常數，FSC 仍為 1，
所以它不要求每個係數的振幅數值完全相同。接近 0 則表示兩張圖在該頻帶缺少線性一致性。Half-map FSC
常用 0.143 交點估計整體解析度；這個閾值建立在獨立 half-map 與特定訊號—雜訊模型的前提上。局部解析度圖
與方向性 FSC 可進一步顯示整體曲線隱藏的空間和方向差異 {cite}`singer2020`。

```{figure} images/pptx/s24_2.png
:width: 72%
:name: fig-fsc-unmasked

這條 unmasked FSC 曲線在約 $0.267\ \mathrm{\mathring A^{-1}}$ 穿過 0.143，對應解析度
$1/0.267\approx3.74\ \mathrm{\mathring A}$。先找交點，再檢查高頻是否異常上升；加入 mask 後，還要把
masked 與 unmasked 曲線並排比較。
```

讀 FSC 曲線時，先確認 half-set 如何建立，再查看 map 是否套用 mask、mask 如何產生，以及曲線是否校正 mask effect。若曲線在高頻突然上升或 masked 與 unmasked 結果差異很大，應回頭檢查 mask 與 half-set 的獨立性。

## Mask 能提高穩定性，也能製造相關

Mask 排除大面積溶劑，可降低雜訊並改善數值穩定性；邊界太硬會在 Fourier 空間產生 ringing。更危險的是用 full map 的高頻細節製作緊貼結構的 mask，再把同一個 mask 套回兩張 half-maps：mask 本身攜帶的共同形狀可能提高 FSC。實務上需要柔邊、獨立或低解析來源的 mask，並比較 masked、unmasked 或經校正的曲線。

過度擬合也可能來自把一半資料估出的高頻姿態資訊帶到另一半。High-resolution noise substitution、phase randomization 或獨立重跑可協助檢查；每種方法都有自己的適用條件。

## 從幾張圖一起讀重建結果

一條 FSC 曲線無法呈現所有細節。把下列資訊放在一起看，較容易找出結果中的弱點：

- **重現性**：更換 half-set 的分法、隨機種子或合理處理設定後，主要密度特徵是否仍然出現？
- **方向覆蓋**：是否有 preferred orientation 與各向異性解析度？
- **局部品質**：只看整體 FSC 時，是否掩蓋了局部解析度與構形彈性的差異？
- **資料—模型一致性**：投影與反投影的殘差若帶有系統性結構，通常表示模型仍漏掉影像中的規律。
- **處理偏差**：粒子篩選、初始模型、對稱性與 mask 是否把預期答案帶進結果？
- **生物背景**：密度中的組成與構形可與生化或功能實驗相互比較。

原子模型擬合屬於後續主題。進行擬合後，還要查看模型幾何、map-to-model agreement 與交叉驗證；這些量測針對原子模型，與 half-map FSC 的問題不同。

## 延伸閱讀

- Fourier slice、姿態估計與重建的完整推導可參考 Singer 與 Sigworth 的綜述 {cite}`singer2020`。
- Maximum-likelihood 與 Bayesian refinement 的方法背景可參考 Scheres 等人的工作 {cite}`scheres2010`。
- 實際查看 half-map、mask 與 FSC 時，可搭配 [RELION post-processing 文件](https://relion.readthedocs.io/en/release-5.0/)閱讀。

## 理解檢查

1. 偏好取向會在 3D Fourier 空間留下什麼取樣問題？

```{dropdown} 參考答案
每張 2D 投影會在 3D Fourier 空間提供一個通過原點的平面。粒子集中在少數取向時，這些平面也會集中在少數方向，其他方向便留下稀疏取樣區域，造成各向異性解析度與方向性偽影。

增加更多相同取向的粒子，可以提高已取樣方向的訊雜比。缺少的方向仍然空白；實際缺口形狀取決於取向分布、對稱性與 CTF，應搭配方向分布圖與局部或方向性解析度評估。
```

2. 為什麼機率加權在低訊雜比資料中比單一最佳方向更能表達不確定性？

```{dropdown} 參考答案
低訊雜比影像常有多個方向與位移得到相近的概似值。機率加權會計算

$$
r_i(z)\propto p(I_i\mid z,V)p(z),\qquad \sum_z r_i(z)=1,
$$

並讓所有候選依權重參與模型更新。權重的分散程度直接反映姿態的不確定性；單一最佳方向只留下最大值的位置。

因此，讀後驗權重時也要記下候選網格、先驗、雜訊模型與初始結構；這些設定都會改變權重分布。
```

3. 哪些資訊若跨 half-set 傳遞，會讓 FSC 過度樂觀？

```{dropdown} 參考答案
高頻 map、由其中一半估出的高頻姿態，以及從完整資料高頻細節製作的緊密 mask，都可能把共同結構帶到另一半。重複觀測同一個物理粒子或同一段 movie，若被分到兩邊，也會造成資料洩漏。

Half-set 獨立的重點是讓兩邊的高頻雜訊無法互相學習。事先固定的低解析初始模型與共同處理規則通常可以共享；使用兩邊高頻資料重新調整規則，則會破壞獨立性。
```

4. FSC 的 0.143 交點能說明什麼？Masked FSC 升高時，還要檢查哪些資訊？

```{dropdown} 參考答案
在兩張 half-map 維持獨立、FSC 定義與 mask 校正都清楚的前提下，0.143 交點可作為整體解析度的慣用估計。FSC 衡量兩張 half-map 在各 Fourier shell 的一致性；局部解析度圖與方向性 FSC 會顯示整體數值中看不到的空間與方向差異。

Masked FSC 升高可能來自合理排除溶劑雜訊，也可能來自 mask 讓兩張 half-map 帶有相同邊界。判讀時要查看 mask 的來源、柔邊寬度、masked 與 unmasked 曲線、phase randomization 或 noise substitution 檢查，以及局部解析度和方向分布。

因此，0.143 交點提供的是整體一致性尺度，需和 mask 資訊、局部解析度與方向分布一起閱讀。
```
