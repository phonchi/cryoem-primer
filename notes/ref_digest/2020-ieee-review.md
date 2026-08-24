# [2020] Single-Particle Cryo-EM: Mathematical Theory, Computational Challenges, and Opportunities

- 作者：Tamir Bendory, Alberto Bartesaghi, Amit Singer
- 出處：*IEEE Signal Processing Magazine* 37(2): 58–76, March 2020，doi:10.1109/MSP.2019.2957822
- 頁碼對照：印刷頁 58–76 ＝ PDF 頁 1–19（`印刷頁 = PDF 頁 + 57`）；本檔頁碼一律用印刷頁。
- 原檔：`References/[2020] IEEE_Review.pdf`

## 一句話定位
從訊號處理視角寫的 cryo-EM 全景 review：有完整的 12 階段 pipeline 流程圖、明確的影像形成數學模型與 CTF 公式，是 06 章 workflow 與 07 章模擬的主要引用來源。

## 關鍵主張與內容
- 背景：2015 Nature Methods「Method of the Year」、2017 諾貝爾化學獎（Dubochet / Frank / Henderson）；2013 前解析度優於 7 Å 者寥寥，被戲稱 blob-ology，DED 出現後才有 resolution revolution（p.58–60）。
- 解析度的生物意義：9 Å 可分辨 α-helix、4.8 Å 可分辨 β-strand、3.5 Å 可分辨多數側鏈（p.59）。
- **cryo-EM vs X-ray vs NMR**（05 章直接可用）：X-ray 三大弱點——難結晶（尤其膜蛋白）、crystal contact 扭曲構形、繞射訊號是所有分子的集體平均因而看不到結構變異；cryo-EM 樣品製備簡單、玻璃態冰保留近生理狀態、每顆粒獨立成像故可解多個功能態；NMR 可在生理條件下量測但僅限約 50 kDa 以下（p.60–61）。
- SNR 極低的兩個原因：無重金屬染色故對比低、必須低電子劑量以免輻射損傷。SNR 常低於 0 dB，可低到 −20 dB（雜訊功率為訊號的 100 倍）（p.60, 62）。
- **12 階段 pipeline**（Fig.4）：raw movie → motion correction／movie frame alignment → CTF estimation → particle picking（＋pruning）→ 2D classification → ab initio model building → 3D classification（選用）→ 3D refinement → per-particle refinement（選用）→ postprocessing → final structure（p.60–61）。
- 雜訊模型：frame 層級主要是 Poisson shot noise；平均成 micrograph 後慣例假設 Gaussian。雜訊頻譜非白（非彈性散射、冰厚變化、PSF 所致），但假設為 1D 徑向函數，參數由重建過程或無訊號區域估計（p.63）。**這段是 07 章雜訊模擬的直接依據。**
- 重建有三個本質不確定性：整體 3D 旋轉、分子中心位置、手性（handedness / chirality）。反射造成的旋轉共軛關係 R̃ = JRJ⁻¹，J = diag(1,1,−1)（p.62）。
- 低 SNR 下取向估計必然失敗：以簡單的 1D 旋轉估計為例，σ→∞ 時相關峰位置在圓周上均勻分布；Cramér–Rao bound 正比於 σ² 且與樣本數 N 無關（p.63）。因此應考慮繞過旋轉估計的統計方法（ML、method of moments）。
- 重建方法：FBP 不適合 cryo-EM（視角覆蓋不均、Fourier 取樣非均勻）；ART/Kaczmarz 因無法用 FFT 加速而少用；現代做法用 gridding 與 NUFFT（p.63–64）。
- Refinement 分兩類：**hard assignment**（template matching，每張影像指派單一視角，如 EMAN2、cisTEM）與 **soft assignment**（給每個影像—模板對一個 responsibility 權重，即 ML-EM／Bayesian，如 RELION、cryoSPARC）（p.64–65）。
- ML-EM 要點：對 nuisance variable（旋轉、平移）做邊際化，否則參數數隨樣本數增長而不一致（Neyman–Scott paradox）；E-step 算 responsibility 權重、M-step 更新結構（p.65）。RELION 假設 SO(3) 均勻先驗、每個 Fourier 係數獨立高斯，變異數逐輪由同心頻率殼平均更新，故低頻先收斂、高頻逐步放行，即 **frequency marching**；cryoSPARC 用每體素 Poisson 先驗（p.65–66）。主要瓶頸是 E-step 的 5D 搜尋，可用多尺度或 branch-and-bound 加速（p.66）。
- Model bias 與 **"Einstein from noise"**：純高斯雜訊影像與參考影像對齊後平均，會收斂到參考影像本身而非零。教訓：不審慎設計演算法，重建出的會是研究者的預設而非資料（p.66–67）。
- Ab initio 三條路：SGD 於負對數後驗（cryoSPARC，可達約 10 Å）（p.67）；**common lines / angular reconstitution**——由 Fourier slice theorem，任兩張投影的 FT 必共享一條中心線，可定出相對旋轉的兩個尤拉角，第三張影像定出第三個角，穩健版本為 synchronization（譜方法或 SDP）（p.67）；**method of moments**（Kam 1980）——均勻旋轉分布時三階矩唯一決定結構，非均勻分布時二階矩即足夠（p.67–68）。
- Motion correction：MotionCor2 把 micrograph 切 patch、各自以互相關估計位移，再擬合成隨時間變化的 2D 多項式（空間二次、時間三次），最後（可加 dose weighting）加總；更新的做法是 per-particle motion correction（p.69）。
- Particle picking：template matching（RELION 用手挑數百顆 → 2D 分類 → 類平均當模板）、edge detection、機器學習；APPLE picker 用局部均值/變異數挑參考視窗再訓練 SVM 以避免模板偏差；常見 pruning 用多個 picker 取共識或高斯混合模型（p.70）。
- 2D classification：類平均用途是建 ab initio 模型、當 picking 模板、快速檢視資料、剔除壞顆粒、偵測對稱性。RELION 的 ML-EM 版有三個弱點：計算量大、假設只有 K 個類別（實際取向連續分布於 SO(3)）、winner-takes-all 導致只產出少數低解析類別。替代方案：以 bispectrum（旋轉不變特徵）做最近鄰平均，並以 steerable PCA 壓縮（p.70）。
- Steerable PCA：把影像展開在 steerable basis（Fourier–Bessel、2D prolate spheroidal、Zernike）上，共變異數矩陣呈區塊對角，複雜度由 O(NL⁴+L⁶) 降到 O(NL³+L⁴)；高維下用 spiked covariance model 的 eigenvalue shrinkage 去雜訊（p.71）。
- 兩個抽象數學模型：**MRA**（multireference alignment，未知群作用下的訊號回復）與 **MTD**（multitarget detection，直接從 micrograph 重建、不做 particle picking）。低 SNR 下 sample complexity 由「第一個能區分訊號的矩」決定：旋轉均勻時 ∝ 1/SNR³、非均勻時 ∝ 1/SNR²（p.71–72）。
- 異質性：離散模型（K 個 volume，易併入 ML-EM，但 K 大時不 scale）vs 連續模型（構形落在低維子空間或 manifold，用 PCA、diffusion maps）（p.72–73）。
- 驗證與解析度：tilt-pair validation；主流是把資料隨機分兩半各自重建，取兩者仍一致的最高頻率（gold-standard FSC），但兩半用同一流程可能產生相關的系統性錯誤（p.73）。
- 深度學習目前多用於 particle picking、pruning、驗證與去雜訊；在極低 SNR 與監督式 model bias 兩點上仍有疑慮，對 3D 重建與異質性的影響有限（p.74）。

## 可引用的公式/圖表
- **式 (1)(2)（p.62）**：影像形成模型 `Iᵢ = hᵢ * T_{tᵢ} P R_{ωᵢ} z + noise`——旋轉 → 沿 z 軸積分（X-ray transform）→ 平移 → 與 PSF 卷積 → 加雜訊。**07 章合成資料的權威式子。**
- **式 (4)（p.63）**：Fourier slice theorem 的算子形式 `S R F₃ z = F₂ P R z`。
- **式 (5)–(10)（p.65）**：ML-EM 的一般模型、後驗、E-step 權重 w、M-step。
- **式 (15)（p.69）**：CTF 模型 `CTF(k) = −sin(π λ|k|²Δf + (π/2) λ³|k|⁴C_s + α) · E(|k|)`，含 defocus Δf、球差 C_s、相位偏移 α、B-factor 包絡 E。**07 章 CTF 實作直接對照。**
- **式 (16)(17)（p.69）**：phase flipping，`I_corrected = sign(CTF)·I = |CTF|·I`。
- 式 (18)（p.70）2D 分類的生成模型；式 (19)（p.71）MRA；式 (20)–(22)（p.72）MTD 與 micrograph 生成模型；式 (23)(24)（p.73）異質性模型。
- **Fig.1（p.59）**：8 個代表性結構的解析度畫廊（4.5 Å GPCR → 1.6 Å apoferritin）——05 章開場圖。
- **Fig.2（p.60）**：β-galactosidase 的 micrograph、power spectrum 上的 Thon rings 與擬合 CTF、1.9 Å 重建圖與側鏈特寫。**06/07 章解釋 CTF 估計的最佳圖。**
- **Fig.3（p.60）**：EMDB 累積結構數與平均解析度的年度曲線，2013 年出現拐點。
- **Fig.4（p.60–61）**：完整 pipeline 流程圖。**06 章 12 階段架構直接照這張。**
- **Fig.5（p.63）**：模擬 3D 結構與十餘張無雜訊投影——07 章「投影」步驟的示意。
- **Fig.6（p.69）**：global / semi-local / per-particle 三種 motion correction 的軌跡對照。
- **Fig.7（p.71）、Fig.8（p.72）**：MRA 與 MTD 在 σ = 0 / 0.2 / 1.2 下的對照——**極佳的「雜訊如何摧毀對齊與偵測」教學圖，可對應 07 章加雜訊實驗。**

## 適用章節
- **06（主要）**：Fig.4 的 pipeline、各階段（motion correction、CTF、picking、2D 分類、ab initio、refinement、3D 分類、驗證）的原理與代表軟體。
- **07（主要）**：式 (1)(2) 的影像形成模型、式 (15) 的 CTF、雜訊模型段落（p.63）、Fig.5/Fig.7/Fig.8。
- **05（主要）**：解析度的生物意義（p.59）、cryo-EM vs X-ray vs NMR 的三點比較（p.60–61）、DED 與 resolution revolution。
- **03**：式 (4) Fourier slice theorem、FBP/gridding/NUFFT、common lines（p.63, 67）。
- 02（次要）：template matching、edge detection 式的 particle picking（p.70）；PCA 去雜訊（p.71）。
- 01（次要）：資料規模與維度（micrograph 數十至 100 megapixel、資料集達數 TB、200³ volume ≈ 8×10⁶ 體素）（p.64）。
