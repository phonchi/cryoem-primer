# 3D 重建與驗證：從 Fourier 切片到 half-map FSC

```{admonition} 讀完本章，你應該能
:class: important
- 用 Fourier slice theorem 解釋 2D 投影如何約束 3D 結構。
- 分辨 hard projection matching 與機率加權的 pose estimation。
- 說明 gold-standard half-set、FSC、mask 與 overfitting 之間的關係。
- 把「演算法收斂」「FSC 通過閾值」與「結構正確」分成不同主張。
```

(fourier-coverage)=
## Fourier slice theorem 是重建的幾何核心

令 $P_RV$ 為三維結構 $V$ 在取向 $R$ 下的投影。Fourier slice theorem 表示

$$
\mathcal F_2\{P_RV\}(\mathbf k)
=\widehat V\!\left(R^{\mathsf T}[k_x,k_y,0]^{\mathsf T}\right).
$$

一張投影的 2D Fourier transform，對應到 3D Fourier volume 中通過原點的一個平面。若取向已知且覆蓋足夠，可把許多切片插值到 3D 格點，再做 inverse FFT；也可在實空間用 weighted backprojection，或以 ART／SIRT 等迭代法解離散線性系統 {cite}`penczek2010`（pp. 5–24）。

有限方向、非均勻取樣、CTF 零點與插值都讓反演不適定。Preferred orientation 會使 Fourier 空間某些方向取樣較弱，造成各向異性解析度；正則化能穩定解，不能創造沒有量到的資訊。

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

更新模型時使用 $r_i(z)$ 加權所有候選，不必先把影像硬指派到唯一姿態。Bayesian／MAP 方法再加入結構先驗或頻率相依正則化。這些方法比較能表達不確定性，但結果仍取決於 likelihood、noise model、prior、離散網格與初始模型 {cite}`scheres2010,singer2020`（Sigworth et al., pp. 273–277；Singer & Sigworth, pp. 175–180）。

(reconstruction)=
## 從加權切片到 3D volume

實際重建會同時處理 CTF、取樣密度與正則化。對同一 Fourier voxel，來自不同粒子的觀測可帶有不同 CTF 正負號與振幅；合理的合併會保留這些權重，並避免在共同資訊不足處做不穩定除法。不同軟體採用 gridding、backprojection、iterative least squares 或其組合，不能只用一個「反投影」詞概括所有實作。

迭代停止也不是證明。Likelihood 不再明顯改善、兩輪 map 很相似，或方向不再大幅改變，都只代表某個目標函數附近的數值收斂；局部最佳解、model bias 與錯誤 noise model 仍可能存在。

(validation)=
## Gold-standard half-set：把獨立性放進流程

常見做法是在 refinement 前把粒子隨機分成兩半，兩組從此獨立 refinement，得到 half-map A 與 half-map B。兩邊可共享事先確定的處理規則，但不能在高頻反覆交換 map、姿態或 mask 所帶出的資料特定資訊，否則 half-map correlation 會被人為提高。

「gold-standard」描述的是防止高頻過擬合的資料分割原則，不表示整個結構已經通過所有驗證。若同一個物理粒子、同一段 movie 的重複觀測或高度相關資料被分到兩邊，獨立性也會被削弱。

## FSC 衡量一致性，不直接衡量真實性

對半徑為 $s$ 的 Fourier shell，half-map FSC 常寫成

$$
\mathrm{FSC}(s)=
\frac{\sum_{\mathbf k\in s}F_A(\mathbf k)F_B(\mathbf k)^*}
{\sqrt{\sum_{\mathbf k\in s}|F_A(\mathbf k)|^2
\sum_{\mathbf k\in s}|F_B(\mathbf k)|^2}}.
$$

FSC = 1 表示該 shell 的 Fourier 係數完全一致，接近 0 表示缺少線性一致性。Half-map FSC 常用 0.143 crossing 報告 nominal resolution；這是建立在特定獨立 half-map 與訊號—雜訊模型下的慣例，不是所有資料都適用的自然常數，也不等於局部解析度或模型正確性 {cite}`singer2020`（supplement pp. 5–6）。

報告 FSC 時至少要交代 half-set 如何建立、map 是否 masked、mask 如何產生、曲線是否校正 mask effect，以及 crossing 前後是否有非物理解讀。

## Mask 能提高穩定性，也能製造相關

Mask 排除大面積溶劑，可降低雜訊並改善數值穩定性；邊界太硬會在 Fourier 空間產生 ringing。更危險的是用 full map 的高頻細節製作緊貼結構的 mask，再把同一個 mask 套回兩張 half-maps：mask 本身攜帶的共同形狀可能提高 FSC。實務上需要柔邊、獨立或低解析來源的 mask，並比較 masked、unmasked 或經校正的曲線。

Overfitting 也可能來自把一半資料估出的高頻姿態資訊帶到另一半。High-resolution noise substitution、phase randomization 或獨立重跑可協助檢查，但每種診斷都有適用條件。

## 驗證要回答多個問題

一份可辯護的 SPA 結果通常要分開檢查：

- **重現性**：不同 half-sets、seeds 或合理處理設定是否得到相容結果？
- **方向覆蓋**：是否有 preferred orientation 與各向異性解析度？
- **局部品質**：局部解析度與構形彈性是否讓 global FSC 過度簡化？
- **資料—模型一致性**：projection／back-projection residual 是否有系統性結構？
- **處理偏差**：particle selection、initial model、symmetry 與 mask 是否把預期答案帶進結果？
- **生物合理性**：密度中的組成與構形是否得到獨立生化或功能證據支持？

本章不展開 atomic model fitting；即使 global half-map FSC 很好，模型幾何、map-to-model agreement 與交叉驗證仍是另一組問題。

## 理解檢查

1. Preferred orientation 會在 3D Fourier 空間留下什麼取樣問題？
2. 為什麼 soft assignment 在低訊雜比資料中比單一最佳方向更能表達不確定性？
3. 哪些資訊若跨 half-set 傳遞，會讓 FSC 過度樂觀？
4. Masked FSC 比 unmasked FSC 高時，還需要哪些資訊才能判斷這是合理去除溶劑，還是 mask-induced correlation？
5. FSC 0.143 crossing 能支持什麼主張？不能單獨支持什麼主張？
