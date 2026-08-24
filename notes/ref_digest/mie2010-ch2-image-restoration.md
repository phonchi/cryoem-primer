# [2010] Image Restoration in Cryo-Electron Microscopy（Methods in Enzymology Vol. 482, Ch.2）

- 作者：Pawel A. Penczek（University of Texas, Houston）
- 出處：*Methods in Enzymology*, 2010, Vol. 482, pp. 35–72，doi:10.1016/S0076-6879(10)82002-6
- 頁碼對照：印刷頁 35–72 ＝ PDF 頁 1–38（`PDF 頁 = 印刷頁 − 34`）
- 原檔：`References/2010/Chapter-Two---Image-Restoration-in-Cryo-Electron-Mic_2010_Methods-in-Enzymol.pdf`

## 一句話定位
CTF 與影像復原的權威教科書式導論：從線性成像模型、CTF 相位函數推導，到 pseudoinverse／chi-square／Wiener filter 三大復原準則的完整數學脈絡，以及各主流軟體（SPIDER/IMAGIC/EMAN/SPARX）CTF 校正策略的比較，是 05/06 章講「為什麼要做 CTF correction」的核心引用來源。

## 關鍵主張與內容
- **影像復原 vs 影像增強**：復原是有明確數學目標（在雜訊限制下盡量還原原始物體）、線性、盡量不放大雜訊的過程；增強是為了視覺辨識而做的啟發式非線性處理，兩者目標與方法完全不同——SPR 中應優先考慮復原（p.37–38）。
- **弱相位物件近似下的成像模型**（式 2.1，p.38）：2D 投影是 3D 庫倫電位沿光軸積分後，與顯微鏡 psf 卷積，並加上兩種噪聲——溶劑/support film 造成的背景噪聲 mᴮ 與（如使用）碳膜造成的 mˢ（p.38–39）。
- **CTF 與相位擾動函數**（式 2.9–2.10，p.41）：`CTF(s;Δz) = √(1−A²)·sin(γ) − A·cos(γ)`，`γ(s;Δz) = 2π(−½Δzλs² + ¼Csλ³s⁴)`；A 為振幅對比比例（cryo 資料常設 6–12%），Cs 為球差常數，λ 由加速電壓決定。忽略像散時 CTF 為 1D（僅依賴 |s|）函數（p.41）。
- **B-factor 包絡函數**（式 2.6–2.7，p.40）：`E(s) = exp(−Bs²/4)`，`σ_B = √(2/B)` 為振幅衰減至 ~60% 的空間頻率；3σ_B 處振幅已衰減至 ~1%。B-factor 單位 Å²，是**面積**單位這件事常被誤解。
- **SSNR 的三種定義**（p.42–43）：理想 SSNR（式 2.15，僅信號本身）、SSNR^Data（式 2.16，把 CTF 與包絡函數計入有效訊號，較常用）、SSNR^Ali（式 2.17，額外考慮 alignment 誤差造成的等效包絡）——三者關係不直觀，需視資料處理協定而定。
- **CTF 參數估計流程**（p.44–46）：加速電壓與 Cs 視為已知精確值，振幅對比常設為固定值；待估的是 defocus 與像散。做法有二：對整張顯微照片估計 2D power spectrum（優點：統計量大、可用於改善粒子挑選）或對已挑選粒子窗口估計（可額外用於估計包絡函數/SSNR）。背景噪聲 PW 需先從資料 PW 中減去才能擬合純 CTF 曲線（式 2.18，Fig.2.2）。
- **蛋白質理論 PW 模型**（3.1 節，p.46–47）：低頻部分近似高斯（標準差關聯蛋白有效直徑），高頻部分趨於常數；SPARX 採用的經驗式（式 2.19）僅有效依賴兩參數（t₂, t₃），對多數大分子複合體適用良好（Fig.2.3 展示 20S proteasome、GroEL、CaMKII、70S ribosome 四例）。
- **Pseudoinverse 復原**（4.1 節，p.47–53）：忽略雜訊時單張影像的復原即簡單除以 CTF·E（式 2.21–2.22），但 CTF 過零點處資訊不可逆失去，且除法會放大雜訊造成嚴重實空間假影——**單張影像的 CTF 校正本質上有資訊缺口，必須靠收集多種 defocus 的資料互補**。多影像的 chi-square 準則（MCE，式 2.23）給出加權合併公式（式 2.29–2.31），此即許多套件實作「多 defocus 合併」的理論基礎。
- **CMCE（constrained pseudoinverse）**（4.2 節，p.53–54）：以高通濾波器 U(s) 對解施加平滑性約束，避免除以極小值的問題（式 2.32–2.34）；γ=0 時退化為 MCE，γ→∞ 時退化為簡單加權求和。
- **迭代解法**（4.3 節，p.54–56）：式 2.36 為 Richardson 型梯度下降迭代，優點是不需顯式除以 CTF（可自然處理過零點）、可加入正則化（式 2.37），且是**唯一在 3D 重建中同時做 CTF 校正在數學上站得住腳的方法**（direct Fourier inversion 與 filtered backprojection 都沒有一致的 3D CTF 校正版本）。
- **POCS 框架**（4.4 節，p.56–58）：以凸集交集表達先驗知識約束（非負性、緊支撐、預定 Fourier 振幅），迭代投影序列收斂到滿足所有約束的解；在 Random Conical Tilt 的 missing cone 問題上曾嘗試應用但效果有限。
- **Wiener filter（MMSE 準則）**（第 5 節，p.58–62）：最小化 `E{|F−F̂|²}`（式 2.39），解出的權重（式 2.45）與 CTF/包絡/SSNR 直接相關；與 MCE/CMCE 數學形式幾乎一致，**但分母多了常數 1**，這一項使 Wiener filter 自動在 CTF 或 SSNR 趨近零處把振幅抑制到零，避免病態放大——這是 Wiener filter 優於 pseudoinverse 的關鍵性質。高 SSNR 極限下退化為 MCE 解；低 SSNR 極限下退化为对影像做低通濾波（式 2.47），此性質是 MMSE 獨有、MSE/MCE 框架沒有對應物。
- **Sharpening（振幅校正）**（5.2 節，p.62–64）：對已完成的 2D 平均/3D 重建做事後濾波，需獨立估計 P_f（目標結構的理論 PW，可來自 X-ray/SAXS 或經驗模型）與 P_m（背景噪聲 PW），並利用 FSC 得到 SSNR_v，才能做出視覺上可與 X-ray 密度圖比擬的「sharpened」地圖（式 2.53）。
- **異質樣本（heterogeneous data）的問題**（第 6 節，p.64–67）：標準 Wiener filter 假設所有影像訊號相同；若用於異質樣本，低頻會被過度放大（式 2.54，典型放大約 1/A ≈ 10 倍），導致模糊的平均影像。兩種啟發式解法：
  - **Phase flipping**（式 2.55，p.65–66）：僅校正 CTF 正負號、不做振幅校正，優點是可用不變量做高效對齊、避免低頻過度放大，缺點是忽略振幅加權導致平均後 SNR 次優。
  - **Reciprocal space adaptive Wiener filter**（式 2.56–2.60，p.66–67）：引入使用者指定的「同質子集數」參數 N_g，在純同質（N_g=1，退化為標準 Wiener/MCE）與純異質（N_g=N，退化為簡單 CTF 加權平均，式 2.58）之間連續過渡，兼顧振幅加權與異質性穩健度。
- **各軟體 CTF 校正策略比較**（第 7 節，p.67–69）：SPIDER 假設 SSNR 為常數、按 defocus 分組各自重建後用 Wiener filter 合併；IMAGIC 先用 phase-flipped 資料做 2D 對齊分類，只在 3D 重建步驟才校正振幅；EMAN1 用擬合的解析函數；EMAN2 改為直接由資料估計影像形成特徵（更穩健、更廣泛適用）；SPARX 的 2D/3D 都用 reciprocal space adaptive Wiener filter 變體。**目前沒有任何套件被證明 SSNR/包絡函數估計的精細度能顯著改善最終解析度**，方法選擇仍高度啟發式。

## 可引用的公式/圖表
- **式 (2.1)（p.38）**：線性成像模型 `g_n(x) = psf_n * e_n * [∫d(T_n r)dz + m^S_n] + m^B_n`。
- **式 (2.8)–(2.10)（p.41）——CTF 完整定義**：`G_n(s) = CTF(s;Δz_n)E_n(s)[(D(Ts))_{s_z=0} + M^S_n(s)] + M^B_n(s)`；`CTF(s;Δz)=√(1−A²)sin(γ)−Acos(γ)`；`γ(s;Δz)=2π(−½Δzλs²+¼Csλ³s⁴)`。**05 章解釋 CTF 過零點/相位反轉的必用式。**
- **式 (2.6)–(2.7)（p.40）**：B-factor 包絡 `E(s)=exp(−Bs²/4)`，`σ_B=√(2/B)`。
- 式 (2.15)–(2.17)（p.42–43）：SSNR、SSNR^Data、SSNR^Ali 三種定義。
- **式 (2.19)（p.47）**：SPARX 經驗 1D PW 模型，兩參數（t₂, t₃）決定形狀；Fig.2.3（p.48）四種蛋白複合體的擬合實例。
- **Fig.2.2（p.46）**：PW 估計流程圖（micrograph PW / 背景 PW / 背景相減後擬合 CTF+包絡），e2ctf.py 實際輸出範例。**05 章講 CTF 估計流程的最佳配圖。**
- 式 (2.21)–(2.22)（p.49）：單張影像 pseudoinverse 復原。
- 式 (2.29)/(2.31)（p.51–52）：多 defocus MCE 合併公式。
- 式 (2.34)（p.54）：CMCE 解。
- **式 (2.45)（p.59）——Wiener filter 完整式**：`F̂(s) = Σ CTF_n E_n SSNR_n G_n / (Σ CTF_n² E_n² SSNR_n + 1)`。**06 章「多 defocus 合併」與 05 章「CTF correction 方法」的核心式，與 mie2010-ch1 式(1.27) 對照引用效果最佳。**
- **Fig.2.4（p.50）**：pseudoinverse filter（A 圖，可見無窮多尖峰）vs Wiener filter（C 圖，SSNR 不同時的平滑抑制曲線）直接對比。**02/05 章講「除以 CTF 為何危險」的最佳圖。**
- 式 (2.53)（p.64）：sharpening 公式。
- 式 (2.55)（p.65）：phase flipping，`G̃_n(s)=sign(CTF_n(s))G_n(s)`。
- 式 (2.59)–(2.60)（p.67）：reciprocal space adaptive Wiener filter 與插值參數 α(s)。

## 適用章節
- **05（主要）**：CTF 定義與相位函數推導（式 2.8–2.10）、B-factor 包絡函數、SSNR 三種定義、CTF 估計流程（Fig.2.2）、多 defocus 資料收集的理論理由。
- **06（主要）**：Wiener filter（式 2.45）、pseudoinverse/CMCE 三種復原準則比較、phase flipping vs adaptive Wiener filter 的異質樣本處理、各主流軟體 CTF 校正策略對照表（第 7 節）——與 mie2010-ch1 的式 (1.27) 直接呼應，適合並列引用。
- **03（次要）**：CTF 校正在 Fourier 空間的除法/濾波操作、convolution theorem 的實際應用案例（式 2.1 的成像模型即為卷積模型）。
- **02（次要）**：Fig.2.4 的「除以接近零的濾波器會放大雜訊」是講 deconvolution/inverse filtering 危險性的絕佳教材，可與影像復原/去模糊章節呼應。
- **07（次要）**：模擬合成資料時加入 CTF、包絡函數與可控 SSNR 的理論依據（式 2.1、式 2.6）；B-factor 與雜訊模型可直接用於合成資料的雜訊注入設計。
