# [2020] Computational Methods for Single-Particle Electron Cryomicroscopy

- 作者：Amit Singer（Princeton）、Fred J. Sigworth（Yale）
- 出處：*Annual Review of Biomedical Data Science* 2020, 3: 163–190，doi:10.1146/annurev-biodatasci-021020-093826
- 頁碼對照：正文印刷頁 163–190 ＝ PDF 頁 1–28（`印刷頁 = PDF 頁 + 162`）；正文止於 p.185，p.186–190 為參考文獻。
- 補充材料 `[2020] Math_review_sup.pdf` 自成一套頁碼（p.1–11），本檔以 **sup p.N** 標記。
- 原檔：`References/[2020] Math_review.pdf`、`References/[2020] Math_review_sup.pdf`

## 一句話定位
全套參考文獻中**數學最完整、公式編號最清楚**的一篇：CTF 從繞射物理一路推導到工程模型、ML/MAP-EM 全式、FSC 與 B-factor 的定量關係，是 06/07 章公式的主要出處，也是 05 章成像原理的權威依據。

## 關鍵主張與內容
- 電子波長 λ = h/(m_e v)：100–300 keV 對應約 0.037–0.020 Å；理論上可達次埃解析度，實務因物鏡像差限制在約 1 Å（p.166）。
- 玻璃態冰的三個作用：避免結晶冰破壞分子、避免結晶冰的高對比淹沒微弱訊號、快速冷卻使分子來不及變形（p.166–167）。
- 成像物理：多數電子直接穿透只貢獻 shot noise；非彈性散射造成輻射損傷、不貢獻對比（可用 energy filter 濾除，否則變成 shot noise），且隨冰厚增加；**彈性散射才產生 tomographic projection**，主要成像模式是 phase contrast。內電位：水/冰約 5 V、脂雙層約 6.5 V、蛋白質約 7.5 V——**我們成像的其實是蛋白質與冰的對比**（p.167–168）。
- 解析度的意義：1.5 Å 才能直接定出個別原子；3–4 Å 已足以配出有意義的原子模型（靠化學約束補足）。典型胺基酸殘基約 120 Da，1 MDa 複合體尺度約 200 Å（p.165）。
- cryo-EM vs X-ray vs NMR（05 章可直接用）：X-ray 的難點在純化與結晶，繞射強度只給 Fourier 係數的模平方、相位須另求；NMR 靠核偶極耦合估原子對距離，限於約 200 個殘基以下的小蛋白；**電子每次有用散射事件造成的損傷比 X-ray 小數個數量級**，故理論上只需約 10⁴ 顆粒，相較 X-ray 需要約 10⁸ 個分子（p.165–166）。
- DED 與 resolution revolution：CCD 時代 SPR 幾乎都差於 6 Å（blob-ology）；2012 起 DED 直接偵測電子、DQE 更高且可錄 movie（例如 2 s 錄 40 frames）；實驗發現 **beam-induced motion 可達 20 Å**，且最大位移集中在前幾個 frame；motion correction 後 B-factor 變小、Thon rings 延伸到更高頻（p.170–171）。
- pipeline 分兩段：micrograph 層級（motion correction、CTF 估計、particle picking，屬經典訊號/影像處理）與 particle 層級（2D 分類、3D 重建、異質性分析，屬統計估計與最佳化）（p.171）。
- Motion correction：早期做全域 frame 對位的最小平方；因位移場在視野中不均勻，改以小 patch 互相關估局部位移，再用低次多項式對位置與時間內插擬合，最後（可加 dose weighting）對齊平均（p.171–172）。
- CTF 估計：由式 9 `PSD_I = H²·PSD_f + PSD_b` 出發。此式表面上 ill-posed（PSD_f、PSD_b 都未知），但可解的關鍵是 **PSD_f 與 PSD_b 隨徑向頻率緩變，而 H² 快速震盪**。做法：徑向平均 → 扣背景（局部平均／低次多項式／線性規劃）→ 由局部極小值定 zero crossing 得 defocus；雜訊大時改用背景扣除後 PSD 與理論 CTF 的互相關；散光參數以 2D 互相關求得（p.172）。
- CTF 校正：**phase flipping**（套用 sign(H) 濾波）只修相位不修振幅，但保留雜訊統計、簡單；**Wiener filter** 用 SSNR = PSD_f/PSD_b 兼修振幅，是所有線性濾波中 MSE 最小者，但兩者都救不回 zero crossing 處已完全消失的資訊——真正的解法是合併不同 defocus 的影像（2D 類平均與 3D 重建時自然發生）（p.172–173）。
- Particle picking：偵測可靠度取決於顆粒大小，目前實務下限約 100 kDa，defocus contrast 對約 40 kDa 以下無能為力。**大 defocus 提升低頻對比、利於 picking，但壓抑高頻**，故高解析仍需低 defocus 影像。Volta phase plate 可提升對比但有 phase drift 爭議。模板可用 Gaussian blob 或其差（DoG）、負染投影、或手挑 1,000–10,000 顆算出 10–100 個 2D 類平均；風險是 **"Einstein from noise"** 的模板偏差（p.173）。
- SPR vs 醫學 CT 的兩個關鍵差異：CT 的 pose 已知（線性反問題、有成熟解法），SPR 的 pose 未知（非線性反問題）；SPR 的雜訊高出許多，需大量影像平均（p.174）。
- MLE 必須對 nuisance parameter（旋轉、平移、每張 defocus、雜訊功率譜、normalization）做**邊際化**，否則參數數隨影像數無限增長，估計不一致（Neyman–Scott paradox），表現為過擬合（p.175）。
- ML-EM 的計算複雜度（06 章可引用的硬數字）：SO(3) 離散成 S = O(L³) 個旋轉（O(L²) 個視角 × O(L) 個平面內旋轉）；E 步與 M 步各為 **O(nL⁴ log L)**。典型 n = 10⁵、L = 200 → nL⁴ = 1.6×10¹⁴，再乘上 10×10 平移搜尋（×100）與 NUFFT 常數（×10），**輕易超過 10¹⁷ flops**；瓶頸來自 L⁴ 隨解析度的成長（p.176–177）。
- 加速手段：responsibility 稀疏化（極端情況即 hard assignment）、以 gridding 取代 NUFFT、branch-and-bound 限縮搜尋空間、adaptive EM（粗網格保留累積機率 99.9% 的點再細化，收斂附近可快一個數量級）、Subspace-EM（用 PCA 降維並預算相關）、GPU 實作（p.177–178）。
- Bayesian / MAP-EM：無正則化的 ML-EM 會在高頻長出假訊號（過擬合）。RELION 式先驗假設每個 Fourier 體素為零均值複高斯、變異數 τ² 於迭代中更新（式 23、24），τ 取為 1D 徑向函數並乘上經驗常數 T（建議 2–4，因先驗忽略了 Fourier 係數間的相關性）。此高斯先驗等價於在 Fourier 域做 **Tikhonov 正則化（強度 ∝ 1/τ²）**，也可視為 Wiener 濾波中的 SNR；初始為低解析時 τ 只在低頻大，高頻被強力壓制，因而自然產生 **frequency marching**（p.179–180）。
- Masking 是另一種先驗：可壓縮參數空間、推高解析度，但先跑一輪再用結果定 mask 屬於 **inverse crime**，會給出過度樂觀的結果（p.180）。
- ML-EM 2D 分類的兩個病症：K 有限與視角連續分布不符；以及 **rich-get-richer**——高 SNR 的類別吸走大部分影像，結果只產出少數低解析類別（p.180）。
- SGD：目標函數可拆成每張影像一項，每輪只用 mini-batch 近似梯度；隨機性有助跳出局部極小，但迭代永不收斂到全域最佳，故**只適合低解析 ab initio（約 10 Å），不適合高解析重建**（p.181–182）。
- 異質性是目前公認的最大計算挑戰。離散 3D classification 在 K 小時有效，K 變大時每類影像太少而崩壞、複雜度隨 K 線性成長，實務上大量人工介入（先小 K、丟壞類、再階層分類＋masking）使可重現性存疑（p.182）。
- 連續異質性的四條路線：共變異數估計（用 covariance 版的 Fourier slice theorem——Fourier 空間中任兩個非共線頻率加原點唯一決定一個中心切片，但**需要視角全覆蓋**，比重建的要求嚴格得多；矩陣有 O(L⁶) 項故只能低解析估計）、normal mode analysis、manifold learning／diffusion maps、multibody refinement 與深度學習（p.182–183）。
- 未來挑戰：air–water interface 造成的變性與聚集（PaaZ 以 graphene oxide 支撐膜解決）、preferred orientation（可用傾斜收資料但加劇 beam-induced motion）；**B-factor 是高解析的終極限制**，目前好資料約 B ≈ 100 Å²，代表 5 Å 處訊號已衰到 1/e、2.5 Å 處衰到 1/e²；最大貢獻來自未修正的 beam-induced motion（尤其前 1–2 e/Å²）與輻射損傷（p.184–185）。

## 可引用的公式/圖表
- **式 1–2（p.168）**：弱相位近似下 `|ψ|² ≈ |ψ_u|²[1 + 2Φ(k)cos φ]`，φ = π/2 + πzλ|k|²，故 **`H(k) = cos φ = sin(−πzλ|k|²)`**。注意 z = 0（正焦）時對比消失——這是「為什麼一定要 defocus」的第一原理解釋。
- **式 3–5（p.168–169）**：`I = h * P f`、`P f = ∫ f dz`、`F I = H · F P f`。PSF h 即 CTF 的逆傅立葉轉換。
- **式 6（p.169）**：`|F I|² = H² · |F P f|²` → power spectrum 上的 **Thon rings** 由此而來。
- **式 8（p.170）——07 章 CTF 實作的主公式**：`H(k) = sin(−πzΔλ|k|² + (π/2)C_s λ³|k|⁴ − w) · E(|k|)`，其中 w 為 amplitude contrast 比例、C_s 為球差、E 為包絡；`E(|k|) = exp(−B|k|²/4)`，好的資料 B ≈ 60–100 Å²。散光模型 `Δz = ½(z₁+z₂) + ½(z₁−z₂)cos(2(α_k − α_z))`（p.169）。此模型在約 2.5 Å 以內有效。
- **式 9（p.172）**：`PSD_I = H²·PSD_f + PSD_b`，CTF 估計的出發點。
- **式 10（p.174）**：影像形成模型 `I_i = h_i * P R_i∘φ + ε_i`。
- **式 11–16（p.175–176）**：似然函數、對 SO(3) 的邊際化、邊際對數似然與 MLE 定義。
- **式 17–19（p.176–177）**：E 步的 responsibility γ、M 步的加權最小平方、含 back-projection 算子 Pᵀ 的顯式解。
- 式 20–22（p.179）Bayes/MAP；式 23–24（p.179–180）高斯先驗與 τ² 更新；式 25–27（p.181–182）SGD。
- **Figure 1（p.167）**：PaaZ micrograph（434×434 nm、1.06 Å/px、約 220 顆）、模板比對挑出的顆粒（含一個誤挑的冰球）、2D power spectrum 與擬合 CTF²（defocus 2.66 µm、astigmatism 49 nm）、以及 1D PSD 與理論 H² 疊圖。**06 章講 CTF 估計、07 章驗證模擬結果都用得上。**
- **Figure 2（p.175）——07 章的核心對照圖**：(a) PaaZ 密度圖的投影 → (b) 最佳匹配取向的投影 → (c) 經 CTF 濾波後（defocus 2.6 µm 造成明顯資訊離域與環狀條紋）→ (d) 真實實驗顆粒影像。四格正好是「投影 → CTF → 雜訊」的教學序列。
- **Figure 3（p.181）**：118,000 顆粒的原始影像 gallery 與 RELION 2.0 算出的完整 2D 類平均集（附各類佔比 %）。
- **Figure 4（p.184）**：把原子模型配進 2.9 Å 密度圖（α-helix pitch 5.4 Å 可作尺標）。
- **補充材料的公式**：
  - `sup p.2` 式 1：common lines 的同步化最小平方 `min Σ‖R_i c_ij − R_j c_ji‖²`（SDP 鬆弛）；同頁有 handedness 的共軛關係 `R̃ = J R J⁻¹`, J = diag(1,1,−1)。
  - `sup p.3`：Kam 的 autocorrelation analysis——均勻視角下二階矩等價於 SO(3) 上的自相關，展開係數只定到一組正交矩陣，需三階矩補足。
  - `sup p.4`：steerable PCA 的區塊對角結構（最大區塊 O(L×L)，原共變異數為 L²×L²）、bispectrum 旋轉不變表示、vector diffusion maps。
  - **`sup p.5` 式 2：FSC 定義** `FSC(k) = Σ U_s V_s* / sqrt(Σ|U_s|² · Σ|V_s|²)`；隨機分半重建、**0.143 判準**（與 X-ray 判準對齊）；`sup p.6` 式 3–4 給出 FSC 與 SSNR 的關係，map-to-model FSC 用 0.5 判準。
  - **`sup p.7` 式 5–6：B-factor 與顆粒數的定量關係** `SSNR(k) = N c e^{−Bk²/2}`，取 SSNR(k₀) = 1/3 得 **`k₀² = 2 ln(3Nc)/B`**——即 Wilson plot（k₀² 對 ln N 為直線，PaaZ 資料 B ≈ 100 Å²）。**這是「要多少顆粒才能到某解析度」的可計算公式，06 章講解析度時極有價值。**

## 適用章節
- **07（主要）**：式 8 的 CTF（含 C_s、amplitude contrast w、B-factor 包絡）、式 1–6 的推導、式 10 的成像模型、Figure 2 的四格序列、雜訊來源（shot noise / 非彈性散射）。
- **06（主要）**：兩段式 pipeline、motion correction、CTF 估計與校正（phase flipping vs Wiener）、particle picking、ML/MAP-EM 全式與 O(nL⁴ log L) 複雜度、frequency marching、2D/3D 分類的病症、FSC 0.143 與 Wilson plot。
- **05（主要）**：電子波長、玻璃態冰、phase contrast 與內電位數值、解析度的生物意義、cryo-EM vs X-ray vs NMR、DED 與 resolution revolution、樣品製備（air–water interface、preferred orientation）。
- **03**：Fourier slice theorem 及其 covariance 版本、convolution theorem（式 3→5）、NUFFT/gridding、Wiener filter 與 Tikhonov 正則化。
- 02（次要）：template matching / DoG picking（p.173）、PCA 去雜訊與 steerable PCA（sup p.4）。
