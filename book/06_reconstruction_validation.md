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

一張投影的 2D Fourier transform，對應到 3D Fourier volume 中通過原點的一個平面。若取向已知且覆蓋足夠，可把許多切片插值到 3D 格點，再做 inverse FFT；也可在實空間用 weighted backprojection，或以 ART／SIRT 等迭代法解離散線性系統 {cite}`penczek2010`（pp. 5–24）。

有限方向、非均勻取樣、CTF 零點與插值都讓反演不適定。偏好取向（preferred orientation）會使 Fourier 空間某些方向取樣較弱，造成各向異性解析度。正則化可以穩定數值解；缺少觀測的方向仍然缺乏資訊。

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

更新模型時使用 $r_i(z)$ 加權所有候選，讓影像保留多個可能姿態。Bayesian／MAP 方法再加入結構先驗或頻率相依正則化。這些方法比較能表達不確定性；結果仍取決於概似函數、雜訊模型、先驗分布、離散網格與初始模型 {cite}`scheres2010,singer2020`（Sigworth et al., pp. 273–277；Singer & Sigworth, pp. 175–180）。

(reconstruction)=
## 從加權切片到 3D volume

實際重建會同時處理 CTF、取樣密度與正則化。對同一 Fourier voxel，來自不同粒子的觀測可帶有不同 CTF 正負號與振幅；合併時需保留這些權重，並在所有影像都缺乏資訊的頻率避免不穩定除法。不同軟體可能採用 gridding、backprojection、iterative least squares 或數種方法的組合，實作內容遠超過單一步驟的反投影。

演算法停止迭代，只能說數值已經收斂。概似值趨於穩定、相鄰兩輪 map 很相似，或方向不再大幅改變，都可能發生在局部最佳解附近；初始模型偏差與錯誤的雜訊模型也可能留在結果中。

(validation)=
## Gold-standard half-set：把獨立性放進流程

常見做法是在 refinement 前把粒子隨機分成兩半，兩組從此獨立 refinement，得到 half-map A 與 half-map B。兩邊可共享事先確定的處理規則；高頻 map、姿態，以及由完整資料高頻細節製作的 mask 必須維持分離，否則 half-map correlation 會被人為提高。

「Gold-standard」指的是用資料分割降低高頻過擬合的原則。結構正確性還需要方向分布、局部解析度與生化證據等檢查。若同一個物理粒子、同一段 movie 的重複觀測或高度相關資料被分到兩邊，兩組資料的獨立性也會被削弱。

## FSC 衡量兩張 half-map 的一致性

對半徑為 $s$ 的 Fourier shell，half-map FSC 常寫成

$$
\mathrm{FSC}(s)=
\frac{\sum_{\mathbf k\in s}F_A(\mathbf k)F_B(\mathbf k)^*}
{\sqrt{\sum_{\mathbf k\in s}|F_A(\mathbf k)|^2
\sum_{\mathbf k\in s}|F_B(\mathbf k)|^2}}.
$$

FSC = 1 表示該 shell 的 Fourier 係數完全一致，接近 0 表示缺少線性一致性。FSC 不直接衡量真實性。Half-map FSC 常用 0.143 交點報告整體解析度；這個閾值建立在獨立 half-map 與特定訊號—雜訊模型的前提上。局部解析度、方向各向異性與模型正確性需要其他檢查 {cite}`singer2020`（supplement pp. 5–6）。

報告 FSC 時至少要交代 half-set 如何建立、map 是否套用 mask、mask 如何產生、曲線是否校正 mask effect，以及交點前後是否出現不符合物理預期的變化。

## Mask 能提高穩定性，也能製造相關

Mask 排除大面積溶劑，可降低雜訊並改善數值穩定性；邊界太硬會在 Fourier 空間產生 ringing。更危險的是用 full map 的高頻細節製作緊貼結構的 mask，再把同一個 mask 套回兩張 half-maps：mask 本身攜帶的共同形狀可能提高 FSC。實務上需要柔邊、獨立或低解析來源的 mask，並比較 masked、unmasked 或經校正的曲線。

過度擬合也可能來自把一半資料估出的高頻姿態資訊帶到另一半。High-resolution noise substitution、phase randomization 或獨立重跑可協助檢查；每種方法都有自己的適用條件。

## 驗證要回答多個問題

判讀 SPA 結果時，至少要檢查以下幾點：

- **重現性**：更換 half-set 的分法、隨機種子或合理處理設定後，是否得到相容結果？
- **方向覆蓋**：是否有 preferred orientation 與各向異性解析度？
- **局部品質**：只看整體 FSC 時，是否掩蓋了局部解析度與構形彈性的差異？
- **資料—模型一致性**：投影與反投影的殘差是否帶有系統性結構？
- **處理偏差**：粒子篩選、初始模型、對稱性與 mask 是否把預期答案帶進結果？
- **生物合理性**：密度中的組成與構形是否有獨立的生化或功能證據可相互印證？

原子模型擬合屬於後續主題。整體 half-map FSC 很高時，仍要另外檢查模型幾何、map-to-model agreement 與交叉驗證。

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

常見誤解是把 soft assignment 當成自動排除 model bias 的保證。候選網格、先驗、雜訊模型與初始結構仍會改變後驗權重，結果仍需獨立驗證。
```

3. 哪些資訊若跨 half-set 傳遞，會讓 FSC 過度樂觀？

```{dropdown} 參考答案
高頻 map、由其中一半估出的高頻姿態，以及從完整資料高頻細節製作的緊密 mask，都可能把共同結構帶到另一半。重複觀測同一個物理粒子或同一段 movie，若被分到兩邊，也會造成資料洩漏。

Half-set 獨立的重點是讓兩邊的高頻雜訊無法互相學習。事先固定的低解析初始模型與共同處理規則通常可以共享；使用兩邊高頻資料重新調整規則，則會破壞獨立性。
```

4. FSC 的 0.143 交點能說明什麼？Masked FSC 升高時，還要檢查哪些資訊？

```{dropdown} 參考答案
在兩張 half-map 維持獨立、FSC 定義與 mask 校正都清楚的前提下，0.143 交點可作為整體解析度的慣用估計。FSC 衡量的是兩張 half-map 在各 Fourier shell 的一致性；結構正確性、局部解析度與方向各向異性需要另外評估。

Masked FSC 升高可能來自合理排除溶劑雜訊，也可能來自 mask 讓兩張 half-map 帶有相同邊界。判讀時要查看 mask 的來源、柔邊寬度、masked 與 unmasked 曲線、phase randomization 或 noise substitution 檢查，以及局部解析度和方向分布。

常見誤解是把單一 0.143 數值當成整個結構的品質證明。這個交點只在上述前提下提供整體的一致性尺度。
```
