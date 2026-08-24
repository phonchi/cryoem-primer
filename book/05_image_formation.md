# 影像形成：電子波如何變成像素值

```{admonition} 讀完本章，你應該能
:class: important
- 從弱相位物體近似解釋相位對比的來源。
- 在實空間卷積與 Fourier 空間相乘兩種寫法間正確切換。
- 說明 DED、DQE、dose fractionation 與 dose weighting 的角色。
- 分辨純量 SNR 與頻率相依的 SSNR。
```

## 從投影擴充成成像模型

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

$h_i$ 是 point-spread function，$H_i=\mathcal F\{h_i\}$ 是 CTF，$\mathbf t_i$ 是平面位移。第一式用實空間卷積描述，第二式經 Fourier transform 後成為逐頻率相乘。這個線性模型適合薄的 SPA 樣品與本書的合成資料；多重散射、厚樣品與部分儀器非理想性需要更完整的模型 {cite}`singer2020`（Eqs. 3–5, 10）。

## 弱相位物體為什麼需要離焦

對薄、弱散射的生物樣品，exit wave 可近似為

$$
\psi_{\mathrm{exit}}(\mathbf x)\approx 1+\mathrm{i}\,\sigma V_p(\mathbf x),
$$

其中 $V_p$ 是沿束流方向積分的位能，$\sigma$ 是交互作用常數。常數項可視為未散射的 reference wave，第二項是樣品造成相位改變的散射波。如果物鏡沒有把兩者的相位差轉成強度差，感光元件量到的 $|\psi|^2$ 幾乎沒有線性於 $V_p$ 的對比。

離焦與球面像差改變散射波在各空間頻率的相位，讓它與未散射波干涉後形成可見的 phase contrast。未散射電子參與對比形成，也因有限計數而帶來 shot noise。

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

散光時，離焦會隨 Fourier 平面的方位角改變，Thon rings 會呈現橢圓或更複雜的方向差異。徑向平均適合無散光近似；完整估計需要二維擬合。

(detectors)=
## DED 與 DQE：偵測器也會改變可用訊號

直接電子偵測器（DED）直接記錄電子事件，省去傳統 scintillator 與光學耦合造成的一部分模糊，並能快速讀出 movie frames。偵測量子效率（detective quantum efficiency, DQE）是頻率相依的量，常寫成

$$
\mathrm{DQE}(s)=\frac{\mathrm{SNR}_{\mathrm{out}}(s)}{\mathrm{SNR}_{\mathrm{in}}(s)}.
$$

這裡的 SNR 沿用本書的功率比定義；若文獻把 S/N 定義成振幅比，分子與分母會各自平方，兩種寫法等價。DQE 越接近 1，代表輸入電子所帶的 SNR 保留得越多；它是輸出與輸入 SNR 的比值，涵蓋偵測器的模糊與雜訊特性 {cite}`mcmullan2016`（pp. 2–3）。Counting mode 可降低 readout noise，但電子事件重疊時會有 coincidence loss，因此曝光率仍需控制。DED 促成了 movie-based motion correction；輻射損傷與 CTF 零點仍需分別處理。

(dose-weighting)=
## Dose fractionation 與 dose weighting

Dose fractionation 把一次總曝光拆成多個 frames。總 fluence 固定時，拆分前後的總電子數相同；新增的時間解析度可用來估計 beam-induced motion 與不同曝光階段的資訊衰減。

曝光早期通常保留較多高頻資訊，累積劑量增加後高頻先受輻射損傷；低頻訊號可從較多 frames 受益。因此 dose weighting 會依 frame 與空間頻率調整權重。最早 frames 又可能有較大的運動，實際權重需由資料與方法共同決定 {cite}`sigworth2016,rubinstein2016`（Sigworth, pp. 65–66；Ripstein & Rubinstein, pp. 103–124）。

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

Wiener filter 與解析度評估需要頻率相依資訊。比較「SNR 0.1」時，兩份資料須採用相同定義、頻帶與 clean target，數值才有共同基準。

## 模型邊界

本章模型適合建立 SPA 前向模型的第一層直覺，採用薄樣品、線性成像與簡化雜訊等條件。厚樣品的多重散射、非彈性散射背景、anisotropic magnification、beam tilt、Ewald sphere curvature、粒子逐 frame 的非剛體運動，以及空間非平穩雜訊，均留待更完整的模型處理。用本章設定產生的合成資料，可觀察較理想成像條件下的行為。

## 理解檢查

1. 為什麼正焦的弱相位物體可能幾乎沒有可見對比？

   ```{dropdown} 參考答案
   對實值 projected potential $V_p$，弱相位近似為

   $$
   \psi_{\mathrm{exit}}\approx1+\mathrm{i}\sigma V_p.
   $$

   正焦且沒有額外 phase shift 時，一階展開只允許我們判斷

   $$
   |\psi_{\mathrm{exit}}|^2
   =1+O\!\left((\sigma V_p)^2\right).
   $$

   線性於 $V_p$ 的交叉項相消，所以一階模型預測影像接近均勻強度。二階項的係數要從更高階 exit-wave 模型推導；若採精確的純相位波 $e^{\mathrm{i}\sigma V_p}$，其強度恆為 1。離焦與球面像差讓不同空間頻率累積相位差，未散射 reference wave 與散射波干涉後便產生一階 phase contrast。增加離焦量可提高部分低頻對比，也會讓 CTF 零點變密並加重高頻衰減，實際設定需在兩者之間取捨。
   ```

2. 平移在實空間與 Fourier 空間分別如何出現？

   ```{dropdown} 參考答案
   若未平移投影為 $p_i(\mathbf x)=P_{R_i}V(\mathbf x)$，平移 $\mathbf t_i$ 後寫成 $p_i(\mathbf x-\mathbf t_i)$；再與 point-spread function 卷積得到

   $$
   I_i(\mathbf x)=h_i*p_i(\mathbf x-\mathbf t_i)+\varepsilon_i(\mathbf x).
   $$

   Fourier translation theorem 把同一個平移轉成 phase ramp：

   $$
   \widehat I_i(\mathbf k)=H_i(\mathbf k)\widehat p_i(\mathbf k)
   e^{-2\pi\mathrm{i}\mathbf k\cdot\mathbf t_i}+\widehat\varepsilon_i(\mathbf k).
   $$

   例如沿 $x$ 方向平移 2 pixels，若 $k_x$ 的單位是 cycles/pixel，相位會多出 $-2\pi(2k_x)$；若 $k_x$ 是寬度 $M$ 的 DFT 整數索引，則寫成 $-2\pi(2k_x/M)$。兩種寫法的 magnitude 都保持不變。這個結論假設採用相同 Fourier convention 且平移沒有因裁切丟失內容；zero padding、邊界與 pixel／Å 單位都要寫清楚。
   ```

3. DQE 為什麼依空間頻率而變？Dose fractionation 又多提供了哪些資訊？

   ```{dropdown} 參考答案
   偵測器的 point-spread function 會讓高頻訊號比低頻訊號更容易衰減，電子散射、readout noise 與 counting coincidence 對不同頻率的影響也不一樣。因此

   $$
   \mathrm{DQE}(s)=\frac{\mathrm{SNR}_{\mathrm{out}}(s)}
   {\mathrm{SNR}_{\mathrm{in}}(s)}
   $$

   必須保留頻率 $s$；單一百分比會藏起低頻與高頻表現的差異。舉例來說，某偵測器可在低頻有 DQE 0.9，在 Nyquist 附近降到 0.3，表示兩個頻帶保留的輸入 SNR 比例不同。

   Dose fractionation 在總 fluence 固定時，把電子事件分配到多個 frames，總電子數維持不變。新增的時間軸可估計 beam-induced motion，並觀察高頻資訊隨累積劑量衰減；dose weighting 才據此調整各 frame、各頻率的權重。若每個 frame 的 readout noise 很高或事件率造成 coincidence loss，切得更細仍可能付出代價。
   ```

4. 比較兩個模擬資料集的 SNR 前，至少要核對哪些定義？

   ```{dropdown} 參考答案
   應核對 SNR 是功率比或振幅比、是否換算成 dB，以及 signal 與 noise 各指什麼。Clean target 可能是純投影，也可能是 CTF-filtered projection；兩者的訊號功率不同。還要列出是否扣平均、使用哪個 mask、資料正規化方式、估計頻帶，以及雜訊是白雜訊或 colored noise。

   若功率 SNR 為 0.1，dB 值是 $10\log_{10}(0.1)=-10\ \mathrm{dB}$，對應的 RMS 振幅比為 $\sqrt{0.1}\approx0.316$。若另一篇文章把振幅比 0.1 直接稱為 SNR，它的功率比其實是 0.01，也就是 $-20\ \mathrm{dB}$。單一 SNR 也會壓掉頻率差異，因此 CTF 或濾波工作最好再比較 $\mathrm{SSNR}(s)$。
   ```
