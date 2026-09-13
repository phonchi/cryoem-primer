# 影像形成：電子波如何變成像素值

{doc}`05_background` 說明了為什麼 SPA 需要大量低劑量影像。要從這些影像估計結構，先要建立一個模型：給定三維結構、粒子的取向與位置，以及顯微鏡的設定，就能預測感光元件會記錄到什麼。這稱為前向模型（forward model）。

前向模型提供了比較候選結構的依據。有了它，「從影像推測結構」就變成一個明確的數學問題：找出一組結構與參數，使預測影像和實際影像盡量相符。本章依序說明投影、平移、對比轉移、偵測與雜訊，以及各步驟的參數與假設。

## 前向模型有哪些量

令 $V(\mathbf r)$ 表示三維散射位能的教學近似。取向 $R$ 下的投影定義為沿電子束方向的積分：

$$
P_RV(\mathbf u)=\int V\!\left(R\,[u_x,u_y,z]^{\mathsf T}\right)\mathrm dz,
\qquad \mathbf u=(u_x,u_y).
$$

這裡的 $R$ 把影像平面的座標送進體的座標系，因此稱為 image-plane-to-volume 旋轉。{doc}`07_synthetic_data` 的模擬程式與 {doc}`06_reconstruction_validation` 的 Fourier slice theorem 都採用同一個約定；跨軟體比較角度時，必須連同矩陣的作用方向一起核對。

第 $i$ 張粒子影像在實空間寫成

$$
I_i(\mathbf x)
=h_i * \left[P_{R_i}V\right](\mathbf x-\mathbf t_i)+\varepsilon_i(\mathbf x),
$$

在 Fourier 空間則是

$$
\widehat I_i(\mathbf k)
=H_i(\mathbf k)\,\widehat{P_{R_i}V}(\mathbf k)
e^{-2\pi\mathrm{i}\mathbf k\cdot\mathbf t_i}
+\widehat\varepsilon_i(\mathbf k).
$$

$h_i$ 是 point-spread function（PSF），$H_i=\mathcal F\{h_i\}$ 是對比轉移函數（contrast transfer function, CTF），$\mathbf t_i$ 是平面位移。第一式用實空間卷積描述，第二式經 Fourier transform 後成為逐頻率相乘。

後面幾章會把整條運算鏈縮寫成一個線性算子。把這張影像的所有未知量收進 $z_i=(R_i,\mathbf t_i,\dots)$，並把影像排成向量，前向模型就是

$$
y_i=A_i(z_i)\,V+\varepsilon_i,
$$

其中 $A_i(z_i)$ 依序執行投影、平移與 CTF 調變。這個寫法強調兩件事：對固定的 $z_i$，影像對結構 $V$ 是線性的；但 $z_i$ 未知時，整個問題對未知量而言是非線性的。{doc}`05_statistical_inference` 的概似與 EM、{doc}`06_reconstruction_validation` 的加權最小平方，都建立在這個區分上 {cite}`singer2020`。

下表整理各層參數的物理來源與估計方式。

| 層 | 符號 | 由什麼決定 | 通常如何取得 |
|---|---|---|---|
| 投影方向 | $R_i$ | 粒子在冰中的取向 | 由資料估計（未知） |
| 平面位移 | $\mathbf t_i$ | 裁切框與粒子中心的偏差 | 由資料估計（未知） |
| 對比轉移 | $H_i$ | 離焦、散光、$C_s$、電壓、振幅對比 | 由 power spectrum 擬合 |
| 雜訊 | $\varepsilon_i$ | 電子計數、冰層、偵測器、前處理 | 由背景區域或頻譜估計 |

## 弱相位物體為什麼需要離焦

對薄、弱散射的生物樣品，離開樣品的電子波（exit wave）可近似為

$$
\psi_{\mathrm{exit}}(\mathbf x)\approx 1+\mathrm{i}\,\sigma V_p(\mathbf x),
$$

其中 $V_p$ 是沿束流方向積分的位能，$\sigma$ 是交互作用常數。這就是弱相位物體近似。常數項 1 可視為未散射的 reference wave，第二項是樣品造成相位改變的散射波。

感光元件只能量到強度。在理想正焦且忽略像差時，直接計算強度，一階項會消失：

$$
|\psi_{\mathrm{exit}}|^2
=\left|1+\mathrm{i}\sigma V_p\right|^2
=1+\sigma^2V_p^2
=1+O\!\left((\sigma V_p)^2\right).
$$

線性於 $V_p$ 的交叉項來自 $2\,\mathrm{Re}\{1\cdot\overline{\mathrm i\sigma V_p}\}$，而 $\mathrm i\sigma V_p$ 對實數 $V_p$ 是純虛數，實部為零。因此正焦的弱相位物體幾乎沒有可見對比。若採精確的純相位波 $e^{\mathrm i\sigma V_p}$，其強度恆為 1，結論相同。

離焦與球面像差的作用，是在散射波上加入一個隨空間頻率變化的額外相位。這個相位讓散射波的一部分轉到與 reference wave 同相的方向，兩者干涉後就在強度上留下線性於 $V_p$ 的成分。相位對比的來源因此是干涉，未散射電子在其中扮演參考波的角色，同時也因有限計數帶來 shot noise。

```{dropdown} 從 exit wave 推到線性相位對比
本節的離焦與球差相位 $\chi(\mathbf k)$ 為實偶函數。物鏡的作用可寫成在 Fourier 空間乘上一個相位因子 $e^{-\mathrm i\chi(\mathbf k)}$。對

$$
\widehat\psi_{\mathrm{exit}}(\mathbf k)=\delta(\mathbf k)+\mathrm i\sigma\widehat V_p(\mathbf k),
$$

像平面的波為

$$
\widehat\psi_{\mathrm{img}}(\mathbf k)=\delta(\mathbf k)+\mathrm i\sigma\widehat V_p(\mathbf k)e^{-\mathrm i\chi(\mathbf k)} .
$$

回到實空間並保留一階項，

$$
|\psi_{\mathrm{img}}|^2\approx 1+2\,\mathrm{Re}
\left\{\overline{\mathrm i\sigma \left[V_p * \mathcal F^{-1}\{e^{-\mathrm i\chi}\}\right]}\right\}.
$$

把 $e^{-\mathrm i\chi}=\cos\chi-\mathrm i\sin\chi$ 代入，只有 $\sin\chi$ 的部分會與常數參考波產生實值的一階項，於是強度的一階成分在 Fourier 空間正比於 $+2\sigma\widehat V_p(\mathbf k)\sin\chi(\mathbf k)$。$\chi\equiv0$ 時這一項為零，這正是正焦沒有對比的原因。加入振幅對比比例 $w$ 後，$\cos\chi$ 的部分也會貢獻，得到本章使用的完整 $H(s)$。不同教材對 $\chi$ 的正負號與 Fourier 指數慣例不同，實際比較時要連同約定一起看。
```

## CTF：振盪、零點與包絡

忽略散光時，一個常見的徑向 CTF 模型是

$$
H(s)=\left[\sqrt{1-w^2}\sin\chi(s)-w\cos\chi(s)\right]E(s),
$$

$$
\chi(s)=2\pi\left(-\frac{1}{2}\Delta f\lambda s^2
+\frac{1}{4}C_s\lambda^3s^4\right).
$$

$s$ 是空間頻率、$\Delta f$ 是離焦、$\lambda$ 是電子波長、$C_s$ 是球面像差、$w$ 是振幅對比比例，$E(s)$ 是包絡。採用 $E(0)=1$ 的包絡正規化。零頻率處 $\chi(0)=0$，所以 $H(0)=-w$；$w=0.15$ 時 $H(0)=-0.15$。{doc}`07_synthetic_data` 使用同一組公式與慣例。不同軟體對離焦正負號、Fourier 指數與額外 phase shift 的慣例可能不同；比較參數時要連同 convention 一起看。

電子波長由加速電壓決定。以相對論修正的近似式計算，200 kV 對應 $\lambda\approx0.0251\ \mathrm{\mathring A}$，300 kV 對應 $\lambda\approx0.0197\ \mathrm{\mathring A}$。這兩個數字會反覆出現在下面的估算中。

### 零點的位置可以直接算

先忽略 $C_s$ 與 $w$，$H$ 的第一個非零頻率的零點出現在 $|\chi|=\pi$，也就是 $\Delta f\,\lambda\,s^2=1$。解出對應的實空間尺度：

$$
d_1=\frac1{s_1}=\sqrt{\Delta f\,\lambda}.
$$

取 300 kV、$\Delta f=1\ \mathrm{\mu m}=10^4\ \mathrm{\mathring A}$，得到 $d_1=\sqrt{10^4\times0.0197}=\sqrt{197}\approx14\ \mathrm{\mathring A}$。更高頻的零點依 $|\chi|/\pi=\Delta f\lambda s^2$ 計數：在 7 Å 處這個量是 4.0，在 3.5 Å 處是 16.1。也就是說，同一張 1 μm 離焦的影像，訊號在到達 3.5 Å 以前已經反轉了十六次 {cite}`sigworth2016`。離焦加倍，反轉次數也加倍。

在 3.5 Å 處，$C_s$ 項的貢獻為 $\tfrac12 C_s\lambda^3s^4$。取 $C_s=2.7\ \mathrm{mm}=2.7\times10^7\ \mathrm{\mathring A}$，這個量約為 $0.69$（以 $\pi$ 為單位），和離焦造成的 16.1 相比不大，但足以在高頻把零點位置移開，因此高解析度分析必須把 $C_s$ 寫進模型。

### 離焦要估得多準

既然 CTF 的正負號決定訊號方向，離焦估計錯誤會讓不同影像在同一頻率互相抵消。可以量化這個要求：離焦誤差 $\delta$ 造成的相位改變為 $\Delta\chi=-\pi\delta\lambda s^2$。當 $|\Delta\chi|=\pi$ 時，該頻率的對比完全反相，足以造成這種完全反相的離焦誤差尺度為

$$
\delta=\frac{1}{\lambda s^2}=\frac{d^2}{\lambda}.
$$

在 300 kV、$d=3.5\ \mathrm{\mathring A}$ 時，$\delta=12.25/0.0197\approx620\ \mathrm{\mathring A}$，也就是約 62 nm。這時對比已完全反相；實際估計希望相位誤差遠小於 $\pi$。目標解析度愈高，造成同樣相位誤差的離焦誤差尺度依 $d^2$ 縮小。離焦與散光需要由每張 micrograph 的二維頻譜估計，並依目標頻率判斷其精度是否足夠。

### 散光、包絡與 PSF 的延伸

散光時，離焦會隨 Fourier 平面的方位角改變，通常以 $\Delta f_1$、$\Delta f_2$ 與散光方位角三個參數描述。Thon rings 於是從圓形變成橢圓。徑向平均適合無散光近似或初始估計；要估出散光大小與方向，必須保留二維的非圓對稱資訊。

包絡 $E(s)$ 描述部分相干性、殘餘移動與輻射損傷等因素造成的衰減。CTF 的零點可藉由不同離焦提供互補資料；更換離焦卻不會消除共同的包絡衰減。對非零的 $E(s)$ 做補償仍可估計訊號，但也會放大雜訊；訊號已弱到無法可靠辨認時，補償不會增加原先量測到的資訊。

CTF 在實空間的對應是 PSF。振盪的 $H$ 會把一個點的訊號攤到周圍，稱為 delocalization，其半徑約為

$$
r\approx\frac{\lambda\,\Delta f}{d}.
$$

取 $\lambda=0.0197\ \mathrm{\mathring A}$、$\Delta f=2.2\ \mathrm{\mu m}=2.2\times10^4\ \mathrm{\mathring A}$、$d=3.5\ \mathrm{\mathring A}$，得到約 124 Å {cite}`sigworth2016`。裁切粒子時，方框在粒子外各側需預留約這個距離，總寬約比粒子直徑多 $2r$，否則被攤出去的高解析度訊號會被切掉。

```{figure} images/pptx/s08_1.png
:width: 75%
:name: fig-ctf-theory
不同離焦下的對比轉移示意。曲線的振盪由 $\chi(s)$ 決定，整體高頻衰減則來自包絡；兩者是不同機制。
```

(ctf-correction-example)=
````{dropdown} 數值例：CTF 參數、phase flipping 與 Wiener 校正
第一張圖比較純相位 CTF 與加入振幅對比、phase shift、包絡後的曲線，參數為 300 kV、離焦 1.5 µm、球差 2.7 mm。此例的程式採

$$
\gamma(s)=\pi\lambda\Delta f s^2-\tfrac\pi2 C_s\lambda^3s^4+\phi,
\qquad H=-E\,[\sqrt{1-w^2}\sin\gamma+w\cos\gamma].
$$

它與本章採用的 $\chi=-\gamma$ 寫法相容。比較不同公式時，需同時核對相位定義與 CTF 前面的符號。

```{figure} images/spa/ctf_terms_migrated.png
:alt: 純相位 CTF 與加入振幅對比、相位位移及包絡後的曲線。
曲線比較使用振幅對比 0.1、額外相位 0.15 rad 與包絡 B=40 Å²；每個參數改變的效應要分開觀察。
```

再取 70S 三維密度沿一個固定方向的投影，以 2.82 Å/pixel 取樣，加入 CTF 與白高斯雜訊。此圖使用振幅對比 0.1、B=30 Å²、雜訊標準差為乾淨投影標準差的 0.03 倍；Wiener 項取 $K=0.03$。

```{figure} images/spa/ctf_correction_migrated.png
:alt: 70S 理想投影、經 CTF 加雜訊的觀測、phase flipping 及 Wiener 校正的四圖比較。
Phase flipping 乘上 $\operatorname{sign}(H)$，保留各非零頻率的振幅；Wiener 校正乘上 $H^*/(|H|^2+K)$，同時調整振幅與抑制不穩定頻帶。各圖採各自顯示範圍，對比強弱要連同公式解讀。
```

這些運算使用同一張觀測影像，零點附近的訊號仍受到資料限制。完整成像與已知真值的產生方式見{doc}`07_synthetic_data`。
````

(detectors)=
## DED 與 DQE：偵測器保留了多少資訊

直接電子偵測器（direct electron detector, DED）直接記錄電子事件，省去傳統 scintillator 與光學耦合造成的一部分模糊，並能快速讀出 movie frames。描述偵測器表現的量是偵測量子效率（detective quantum efficiency, DQE）：

$$
\mathrm{DQE}(s)=\frac{\mathrm{SNR}_{\mathrm{out}}(s)}{\mathrm{SNR}_{\mathrm{in}}(s)}.
$$

這裡的 SNR 沿用本書的功率比定義。若資料把 S/N 定義成振幅比，分子與分母會各自平方，兩種寫法等價 {cite}`mcmullan2016`。DQE 越接近 1，代表輸入電子所帶的資訊保留得越多；完美偵測器的 DQE 恆為 1。

DQE 必須畫成頻率的函數，因為偵測器的模糊（由 modulation transfer function 描述）與各種雜訊來源對不同頻率的影響不同。McMullan 等人（2016）比較的 300 keV 直接偵測器在 integrating mode 下，半 Nyquist 頻率的 DQE 約在 40%–60%，到 Nyquist 通常降到 25% 左右；counting mode 把每個電子事件換成一個理想計數，可讓零頻率的 DQE 提高到 80% 上下。作為對照，早期的 phosphor／fiber-optic／CCD 相機在 300 keV、半 Nyquist 處的 DQE 只有 7%–10%，底片約為 30%–35% {cite}`mcmullan2016`。

### 從電子計數算一次 SNR

先從沒有偵測器損失的理想計數，算出輸入 SNR，再看 DQE 如何影響輸出。設某個像素平均收到 $n$ 個電子。電子到達是 Poisson 過程，所以計數的變異數也是 $n$。若粒子區域與背景的相對對比為 $c$，訊號的大小是 $cn$ 個電子，於是功率訊雜比為

$$
\mathrm{SNR}_{\mathrm{in}}=\frac{(cn)^2}{n}=c^2n .
$$

取 $40\ \mathrm{e^-/\mathring A^2}$ 的總 fluence 與 $1\ \mathrm{\mathring A^2}$ 的像素，$n=40$；低頻對比取 $c=5\%$，得到 $\mathrm{SNR}_{\mathrm{in}}=0.0025\times40=0.1$。這個 $-10\ \mathrm{dB}$ 是指定像素、已知背景與弱對比近似下的功率比；{doc}`07_synthetic_data`另以整疊影像的變異數比定義 SNR，兩者的估計範圍不同。

再用另一個例子看頻率相依的 DQE：假設兩台偵測器在同一指定頻率的輸入功率 SSNR 都為 0.1，DQE 分別為 0.5 與 0.08，輸出 SSNR 就是 0.05 與 0.008。這個輸入頻帶功率是另外給定的條件，無法直接從前面的單像素對比推得。在其餘成像條件相同、影像可以正確對齊且雜訊獨立時，要達到相同的輸出 SSNR，後者需要 $0.5/0.08=6.25$ 倍粒子數。

若同一部偵測器在半 Nyquist 的 DQE 是 0.5、Nyquist 是 0.25，它在後者保留的輸入 SSNR 比例較低。判斷高頻有多少可用資訊時，還要把各頻率原本的訊號與雜訊功率納入，才能由 DQE 算出實際輸出 SSNR。

Counting mode 可大幅降低 readout noise，但電子事件在時間或空間上重疊時會被算成一個，造成 coincidence loss，因此曝光率仍須控制。DED 提供的快速讀出也使 movie-based motion correction 成為可能；輻射損傷與 CTF 零點則要在後續步驟處理。

(dose-weighting)=
## Dose fractionation 與 dose weighting

Dose fractionation 把一次總曝光拆成多個 frames。總 fluence 固定時，拆分前後的總電子數相同；新增的是時間解析度。時間資訊可用來估計 beam-induced motion，也可用來比較各曝光階段保留的高頻訊號。

曝光早期通常保留較多高頻資訊，累積劑量增加後高頻先受輻射損傷；低頻訊號則可從較多 frames 受益。因此 dose weighting 會依 frame 與空間頻率調整權重，讓每個頻率主要採用它還可靠的那一段曝光。最早的 frames 又可能帶有較大的運動，所以 motion correction 與 dose weighting 需要搭配使用 {cite}`sigworth2016,rubinstein2016`。實作細節與三個層級的運動估計見 {doc}`06_workflow`。

```{admonition} 三個容易混淆的量
:class: warning
- **總 fluence**：整部 movie 累積的 $\mathrm{e^- / \mathring{A}^2}$。
- **每 frame fluence**：總 fluence 在時間上的分配。
- **exposure rate**：單位時間或 detector area 的電子事件率，會影響 counting coincidence。
```

(noise-model)=
## 雜訊：從 Poisson 計數到影像的高斯近似

$\varepsilon_i$ 的分布與相關性決定了概似函數的形式。要選擇合適的雜訊模型，需要分清電子計數、影像前處理與像素間相關性。

先看電子計數。到達某個像素的電子數服從 Poisson 分布，其變異數等於平均值。這個分布依賴平均計數，因此嚴格說來雜訊的大小隨影像內容改變。

接著看處理過的影像。一部 movie 經過對齊、加權、相加、內插與正規化之後，每個像素值是許多電子事件的線性組合。計數夠多時，中央極限定理讓結果接近高斯分布，於是 SPA 的模型通常寫成

$$
\varepsilon_i\sim N(0,\Sigma_i).
$$

這是一個工作假設：它把 Poisson 的內容相依變異數換成一個固定的共變異數，換取可解析處理的概似函數。

最後還要考慮像素間的相關結構。$\Sigma_i$ 是否為對角矩陣，取決於處理過程。內插、後處理濾波與偵測器事件擴散都可能讓鄰近像素的誤差相關；冰層厚度變化與非彈性散射背景則帶來大尺度的低頻起伏 {cite}`scheres2005`。$\Sigma_i=\sigma^2 I$ 的白雜訊模型因此只是近似。前面成像 CTF 的 PSF 改變的是預測的平均訊號；理想 Poisson 計數在不重疊像素間仍可獨立。

### 共變異數、功率頻譜與 whitening

把雜訊視為平穩隨機場時，它的統計可用自共變異數

$$
C_\varepsilon(\Delta\mathbf x)=\mathbb E\left[\varepsilon(\mathbf x)\,\varepsilon(\mathbf x+\Delta\mathbf x)\right]
$$

描述。自共變異數的 Fourier transform 就是功率頻譜密度（power spectral density, PSD）：

$$
\mathrm{PSD}_{\mathrm{noise}}(\mathbf k)=\mathcal F\{C_\varepsilon\}(\mathbf k).
$$

白雜訊的自共變異數集中在原點，PSD 是常數；colored noise 的 PSD 隨頻率變化，冷凍電鏡資料通常在低頻遠高於高頻。平穩性讓 Fourier 基底近似對角化 $\Sigma_i$，於是各頻率可以分開處理。

雜訊的頻率分布也會影響影像相似度的計算。以未加權的平方距離比較影像時，低頻的大雜訊會主導結果。**Whitening** 是先套一個轉換 $W_i$ 使 $W_i\Sigma_iW_i^{\mathsf T}=I$；在平穩近似下，這相當於把每個 Fourier 係數除以 $\sqrt{\mathrm{PSD}_{\mathrm{noise}}}$。觀測與預測必須一起轉換。

例如，設某頻帶 A（低頻）的雜訊功率為 100、訊號功率為 10；頻帶 B（高頻）的雜訊功率為 1、訊號功率為 0.5（單位任意）。未加權比較時，A 的雜訊功率是 B 的一百倍，殘差幾乎完全由 A 決定。Whitening 之後，兩個頻帶的貢獻依各自的訊雜比排序：A 為 $10/100=0.1$，B 為 $0.5/1=0.5$。振幅較小的 B，其每係數 SSNR 反而是 A 的五倍。是否更能區分候選結構，還要看候選在這個頻帶的預測差異。{doc}`05_statistical_inference` 會說明這個加權如何從高斯概似自然導出。

### SNR 與 SSNR 要帶著定義使用

純量 SNR 把某個區域或頻帶內的訊號與雜訊功率壓成一個數：

$$
\mathrm{SNR}=\frac{P_{\mathrm{signal}}}{P_{\mathrm{noise}}} .
$$

Spectral SNR（SSNR）則保留頻率軸：

$$
\mathrm{SSNR}(s)=\frac{\mathrm{PSD}_{\mathrm{signal}}(s)}
{\mathrm{PSD}_{\mathrm{noise}}(s)}.
$$

例如功率 SNR 為 0.1 時，$10\log_{10}(0.1)=-10\ \mathrm{dB}$，RMS 振幅比為 $\sqrt{0.1}\approx0.316$；若 0.1 指的是振幅比，對應的功率比則是 0.01，也就是 $-20\ \mathrm{dB}$。比較兩份資料以前，先確認 signal 指的是純投影或 CTF-filtered projection，再確認是否扣平均、mask 位置、正規化方式、估計頻帶，以及雜訊屬於白雜訊或 colored noise。

單一 SNR 數字會把頻率差異壓掉，因此 Wiener filter、CTF 校正與解析度評估都需要 $\mathrm{SSNR}(s)$ 這種頻率相依的量。{doc}`06_resolution_validation` 會說明如何從兩張 half-map 的相關性估計它。

## 影像出現其他效應時

薄樣品與線性成像近似讓我們能分開描述投影、CTF 與雜訊。樣品變厚時，多重散射與非彈性散射背景會更明顯；Thon rings 的非圓對稱形狀可提示散光或非等向放大。Beam tilt 造成的主要彗形像差屬於奇像差，會改變相位，通常無法只靠單張 power spectrum 的環紋形狀辨認，需要額外資訊或 refinement 檢查。高解析度資料還可能需要考慮 Ewald sphere curvature。加入這些因素時，要說清楚它改變的是訊號、相位、座標還是雜訊，才能決定新增哪些待估參數。

## 接下來

本章把一張粒子影像寫成 $y_i=A_i(z_i)V+\varepsilon_i$，並說明了 $A_i$ 的每一層與 $\varepsilon_i$ 的統計結構。真實資料裡 $V$ 與 $z_i$ 都是未知的，而且 $z_i$ 逐張不同。下一章用概似、後驗權重與 EM 比較同一張模糊影像的多種可能解釋，進而估計未知結構與參數：{doc}`05_statistical_inference`。

## 延伸閱讀

- **Sigworth, F. J.（2016）．*Principles of Cryo-EM Single-Particle Image Processing*. Microscopy, 65(1), 57–67。** [Oxford Academic 免費全文](https://doi.org/10.1093/jmicro/dfv370)。頁 57–58 的 Figure 1 逐層展示投影、CTF 與雜訊，並標出本章使用的 delocalization 估計；頁 62 的 CTF determination 一節給出 1 μm 離焦在 14 Å、7 Å、3.5 Å 的零點與反轉次數。篇幅短且公開取用，適合作為本章的第一份對照閱讀 {cite}`sigworth2016`。
- **McMullan, G., Faruqi, A. R. 與 Henderson, R.（2016）．*Direct Electron Detectors*. Methods in Enzymology, 579, 1–17。** [DOI](https://doi.org/10.1016/bs.mie.2016.05.056)。頁 2–3 給出 DQE 與 MTF 的定義及 Nyquist 慣例，頁 8–9 列出本章引用的 integrating 與 counting mode DQE 數值以及底片、CCD 的對照；閱讀時注意該文以振幅比的平方定義 DQE {cite}`mcmullan2016`。
- **Singer, A. 與 Sigworth, F. J.（2020）．*Computational Methods for Single-Particle Electron Cryomicroscopy*. Annual Review of Biomedical Data Science, 3, 163–190。** [DOI](https://doi.org/10.1146/annurev-biodatasci-021020-093826)。第 2.1 與 2.2 節（頁 166–170）完整推導電子顯微鏡的成像模型與相位對比，第 2.3 節（頁 170–171）說明直接偵測器與 beam-induced motion；可作為本章前向模型的正式版本 {cite}`singer2020`。
- **Ripstein, Z. A. 與 Rubinstein, J. L.（2016）．*Processing of Cryo-EM Movie Data*. Methods in Enzymology, 579, 103–124。** [DOI](https://doi.org/10.1016/bs.mie.2016.04.009)。第 1 節（頁 103–107）說明 DDD 如何讓一次曝光成為 movie，以及 dose fractionation 帶來的時間資訊；第 6 節（頁 116–121）說明依頻率調整 frame 權重的做法 {cite}`rubinstein2016`。
- **Penczek, P. A.（2010）．*Image Restoration in Cryo-Electron Microscopy*. Methods in Enzymology, 482, 35–72。** [DOI](https://doi.org/10.1016/S0076-6879(10)82002-6)。從 CTF 模型讀到 Wiener filter 與 SSNR 的估計，補上本章 whitening 與 PSD 段落的數學細節 {cite}`penczek2010restoration`。
- **[ASPIRE-Python](https://github.com/ComputationalCryoEM/ASPIRE-Python) 的 `aspire.operators` 原始碼。** 開放原始碼的 `RadialCTFFilter` 與 `CTFFilter` 直接實作本章的 $H(s)$ 與散光版本，可對照參數名稱、單位與正負號慣例；{doc}`07_synthetic_data` 使用同一套介面。
