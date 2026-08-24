# 影像形成：電子波如何變成像素值

```{admonition} 讀完本章，你應該能
:class: important
- 從弱相位物體近似解釋相位對比的來源。
- 在實空間卷積與 Fourier 空間相乘兩種寫法間正確切換。
- 說明 DED、DQE、dose fractionation 與 dose weighting 的角色。
- 分辨純量 SNR 與頻率相依的 SSNR。
```

## 從投影開始，但不能停在投影

令 $V(\mathbf r)$ 表示三維散射位能的教學近似，$P_RV$ 表示將旋轉 $R$ 後的體沿電子束方向積分。理想投影只是前向模型的第一步。第 $i$ 張粒子影像可寫成

$$
I_i(\mathbf x)
=h_i * \left[P_{R_i}V\right](\mathbf x-\mathbf t_i)+\varepsilon_i(\mathbf x),
$$

或在 Fourier 空間寫成

$$
\widehat I_i(\mathbf k)
=H_i(\mathbf k)\,\widehat{P_{R_i}V}(\mathbf k)
e^{-2\pi\mathrm{i}\mathbf k\cdot\mathbf t_i}
+\widehat\varepsilon_i(\mathbf k).
$$

$h_i$ 是 point-spread function，$H_i=\mathcal F\{h_i\}$ 是 CTF，$\mathbf t_i$ 是平面位移。第一式是實空間卷積，第二式才是逐頻率相乘；兩者不能混在同一個等號裡。這個線性模型忽略多重散射、厚樣品與部分儀器非理想性，適合薄的 SPA 樣品與本書的合成資料 {cite}`singer2020`（Eqs. 3–5, 10）。

## 弱相位物體為什麼需要離焦

對薄、弱散射的生物樣品，exit wave 可近似為

$$
\psi_{\mathrm{exit}}(\mathbf x)\approx 1+\mathrm{i}\,\sigma V_p(\mathbf x),
$$

其中 $V_p$ 是沿束流方向積分的位能，$\sigma$ 是交互作用常數。常數項可視為未散射的 reference wave，第二項是樣品造成相位改變的散射波。如果物鏡沒有把兩者的相位差轉成強度差，感光元件量到的 $|\psi|^2$ 幾乎沒有線性於 $V_p$ 的對比。

離焦與球面像差改變散射波在各空間頻率的相位，讓它與未散射波干涉後形成可見的 phase contrast。這也是為什麼「未散射電子只產生背景」不夠準確：它同時參與對比形成，而有限計數仍會造成 shot noise。

## CTF、振幅對比與 envelope

忽略散光時，一個常見的徑向 CTF 模型是

$$
H(s)=\left[\sqrt{1-w^2}\sin\chi(s)-w\cos\chi(s)\right]E(s),
$$

$$
\chi(s)=2\pi\left(-\frac{1}{2}\Delta f\lambda s^2
+\frac{1}{4}C_s\lambda^3s^4\right).
$$

$s$ 是空間頻率、$\Delta f$ 是離焦、$\lambda$ 是電子波長、$C_s$ 是球面像差、$w$ 是振幅對比比例，$E(s)$ 是包絡。不同軟體對離焦正負號、Fourier 指數與 phase shift 的慣例可能不同；比較參數時要連同 convention 一起看。

CTF 的正負號決定對比方向，零點代表該張影像在該頻率沒有線性傳遞訊號。Envelope 則描述 temporal／spatial coherence、樣品移動、輻射損傷與其他因素造成的高頻衰減。零點變密與 envelope 下降是兩種效應。

散光時，離焦會隨 Fourier 平面的方位角改變，Thon rings 因而不再是同心圓。徑向平均適合無散光近似；完整估計需要二維擬合。

(detectors)=
## DED 與 DQE：偵測器也會改變可用訊號

直接電子偵測器（DED）直接記錄電子事件，省去傳統 scintillator 與光學耦合造成的一部分模糊，並能快速讀出 movie frames。偵測量子效率（detective quantum efficiency, DQE）是頻率相依的量，常寫成

$$
\mathrm{DQE}(s)=\frac{\mathrm{SNR}_{\mathrm{out}}(s)}{\mathrm{SNR}_{\mathrm{in}}(s)}.
$$

這裡的 SNR 沿用本書的功率比定義；若文獻把 S/N 定義成振幅比，分子與分母會各自平方，兩種寫法等價。DQE 越接近 1，代表輸入電子所帶的 SNR 保留得越多；它不是「偵測到多少電子」的單一百分比 {cite}`mcmullan2016`（pp. 2–3）。Counting mode 可降低 readout noise，但電子事件重疊時會有 coincidence loss，因此曝光率仍需控制。DED 促成了 movie-based motion correction，但不會消除輻射損傷或 CTF 零點。

(dose-weighting)=
## Dose fractionation 與 dose weighting

Dose fractionation 把一次總曝光拆成多個 frames。若讀出雜訊夠低，拆分本身不會憑空增加總電子數；它提供時間解析度，讓我們估計 beam-induced motion 與不同曝光階段的資訊衰減。

曝光早期通常保留較多高頻資訊，累積劑量增加後高頻先受輻射損傷；低頻訊號可從較多 frames 受益。因此 dose weighting 會依 frame 與空間頻率調整權重，而不是只把後半段 frames 全部丟掉。最早 frames 又可能有較大的運動，實際權重需由資料與方法共同決定 {cite}`sigworth2016,rubinstein2016`（Sigworth, pp. 65–66；Ripstein & Rubinstein, pp. 103–124）。

```{admonition} 三個容易混淆的量
:class: warning
- **總 fluence**：整部 movie 累積的 $\mathrm{e^- / \mathring{A}^2}$。
- **每 frame fluence**：總 fluence 在時間上的分配。
- **exposure rate**：單位時間或 detector area 的電子事件率，會影響 counting coincidence。
```

## SNR 與 SSNR 要帶著定義使用

純量 SNR 把某個區域或頻帶內的訊號與雜訊功率壓成一個數：

$$
\mathrm{SNR}=\frac{P_{\mathrm{signal}}}{P_{\mathrm{noise}}}.
$$

若以影像變異數估計，還要說明是否先扣平均、mask 取在哪裡，以及 clean target 是純投影還是 CTF-filtered projection。Spectral SNR 則保留頻率軸：

$$
\mathrm{SSNR}(s)=\frac{\mathrm{PSD}_{\mathrm{signal}}(s)}
{\mathrm{PSD}_{\mathrm{noise}}(s)}.
$$

Wiener filter 與解析度評估需要的是頻率相依資訊。把不同定義、不同頻帶或不同 clean target 的「SNR 0.1」直接並列，通常沒有可比性。

## 模型邊界

本章模型適合建立 SPA 前向模型的第一層直覺。以下效應沒有被完整納入：厚樣品的多重散射、非彈性散射背景、anisotropic magnification、beam tilt、Ewald sphere curvature、粒子逐 frame 的非剛體運動，以及空間非平穩雜訊。合成資料若沒有模擬這些效應，只能驗證方法在較理想條件下的行為。

## 理解檢查

1. 為什麼正焦的弱相位物體可能幾乎沒有可見對比？
2. 在前向模型中，平移在實空間與 Fourier 空間分別如何出現？
3. 為什麼 DQE 必須寫成頻率的函數？
4. Dose fractionation 提供了什麼可估計資訊？它沒有增加什麼？
5. 比較兩個模擬資料集的 SNR 前，至少要核對哪些定義？
