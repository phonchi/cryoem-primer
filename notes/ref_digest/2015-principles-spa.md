# [2016] Principles of cryo-EM single-particle image processing

- 作者：Fred J. Sigworth（Yale）
- 出處：*Microscopy* 2016, 65(1): 57–67，doi:10.1093/jmicro/dfv370（線上發表 2015-12-24）
- 頁碼對照：印刷頁 57–67 ＝ PDF 頁 1–11（`PDF 頁 = 印刷頁 − 56`）；本檔頁碼一律用印刷頁。
- 原檔：`References/[2016] Principles of cryo-EM single-particle image.pdf`

## 一句話定位
最精簡好讀的 single-particle reconstruction (SPR) 原理導論，適合當 06 章 workflow 的骨架，也提供「為什麼要模擬雜訊與 CTF」的量化理由（SNR 決定一切）給 07 章。

## 關鍵主張與內容
- SPR 能成立靠兩個現象：(1) 急速冷凍的分子構形一致到數 Å 內，資訊可跨顆粒疊加；(2) 極低訊噪的影像本身仍含足以判定取向的資訊（p.57）。
- **單張顆粒影像的高解析內容無法被量測**（原文黑體結論）。TRPV1 資料集中，個別顆粒影像在細於約 16 Å 處 spectral SNR < 1；資料品質只能靠 2D class average、3D 重建、power spectrum、變異數等統計量間接評估（p.58）。
- 影像形成模型：3D 密度沿視線投影 → 乘上 defocus 造成的 CTF → 加上低劑量電子計數雜訊（p.58, Fig.1b–d）。這正是 07 章合成資料的三步驟。
- **Fourier slice theorem**：投影影像的 2D FT 等於 3D 密度 FT 中通過原點、垂直於投影方向的一個切片。由此得直接 Fourier 反演的五步驟流程（取投影 → FFT → 定切片平面 → 填值 → 逆 FFT）（p.59）。
- 直接 Fourier 反演的四個實務難題：插值、切片交會處數值不一致（需加權合併）、Fourier 空間覆蓋不足（missing wedge / missing cone / preferential orientation → 各向異性解析度）、CTF 校正（p.59–60）。
- **CTF 不可能從單張影像「移除」**（零點處資訊全失），但可在 3D 重建時以帶正負號的權重合併多張不同 defocus 的影像來處理（p.60）。
- 取向決定用 projection matching：由初始 3D model 產生大量參考投影，逐張比對取最佳者，再重建、迭代數十次（p.60）。
- 計算複雜度估計 `n_c ≈ t²π²n⁵m/s`（n = 2D/d 為影像邊長、t 為平移不確定度、m 顆粒數、s 對稱階數）；**複雜度隨顆粒尺寸與解析度極陡上升**。D=120 Å、d=3 Å、m=30000、s=4、t²=25 → n_c ≈ 6×10¹⁵ ≈ 1600 CPU-hours（p.60–61）。可用 polar FT、spherical harmonics、steerable basis 把 πn³ 降到 n²log n（p.61）。
- 取向精度需求：要保住 d=3 Å、D=120 Å 的資訊，角度誤差須小於 sin⁻¹(d/D) ≈ 1.5°。模擬顯示 SNR=1 時誤差 1.7°；SNR 減半誤差變 3°（解析度上限約 6 Å）；再減半則誤差 >10°，連低解析重建都做不出來（p.61, Fig.2）。
- 統計加權（maximum-likelihood，如 RELION）不做硬性取向指派，改用取向機率做「模糊」加權，因此對取向誤差較耐受（p.61）。
- **Model bias：最終結構會受初始模型影響**。資料好時從橢球或亂數起始也能收斂；SNR 低或取向偏好嚴重時 bias 明顯。對策：common lines / 隨機聚類做 ab initio 初始模型、random conical tilt、subtomogram averaging，以及事後的 tilt-pair analysis 驗證（p.61–62）。
- Particle picking：difference of Gaussians (DoG) 是出奇有效的通用 2D 模型；用細節豐富的參考會造成 2D model bias（純雜訊也能平均出參考的樣子，即 "Einstein from noise"）。結論金句：**看不到顆粒，多半就是沒有顆粒**（p.62）。
- CTF 判定：300 keV、defocus 1 µm 時第一個零點在 (14 Å)⁻¹；到 7 Å 已翻轉 4 次、3.5 Å 已翻轉 16 次；defocus 4 µm 則翻轉密度再乘四。高 defocus 有利於 picking 與定向（保留 200–20 Å 訊號）但要求更精確的 CTF 模型；**124 nm 的 defocus 誤差就會讓 3.5 Å 處的 CTF 極性完全反轉**；astigmatism 必須一併建模（p.62）。
- 2D classification（RELION 的 ML 版）與 SPR 幾乎同構，差別只在 SPR 的參考來自同一個 3D volume，因而多了自洽性約束；2D 分類無此約束，故不會收斂到唯一結果，主要用途是看資料內容與剔除壞類（p.62–63）。
- 3D classification 與異質性：可做「in silico 純化」，但需事先給定類別數、少數構形常漏掉（p.63）。分辨極限是統計上的 **identifiability problem**——雜訊夠大時兩個母體的分布與單一母體無法區分，SNR 加倍才分得開（p.63–64, Fig.3）。
- 異質性分析工具鏈：3D 變異數 → 共變異數（10⁶ voxel 的共變異數矩陣有 10¹² 項，不可行）→ 取共變異數的主成分（少數自由度即可描述 hinge 等運動）→ 非線性 manifold embedding（p.64）。
- SNR 是 SPR 的根本限制。單顆粒影像 SNR 大致正比於分子量；早年高解析僅限 >1 MDa，電子計數相機出現後 TRPV1（400 kDa）成功，當時最小近原子解析結構為 gamma-secretase（120 kDa）（p.65）。
- 提升 SNR 的途徑：計數相機 DQE 由 0.7/0.4 往 1 提升、薄冰、zero-loss energy filter（非彈性散射率約為彈性的三倍，既降對比又添雜訊）、劑量加權（低劑量約 20 e/Å²；前 1–2 e/Å² 因 beam-induced motion 常須丟棄，高解析資訊集中在早期 frame、低解析資訊可累積整段 movie）、Volta phase plate（把平均平方振幅 ½ 的振盪 CTF 換成接近 1 的常數傳遞，理論上 SNR ×2）。總結預測 SNR 可再提升 3–4 倍，使 <50 kDa 蛋白質成為可能（p.65–66）。

## 可引用的公式/圖表
- **Fig.1（p.58）**：TRPV1 micrograph、boxed particle（256 px、1.22 Å/px）、對應的 3D map 投影、加了 CTF 的無雜訊模擬影像、以及平均 spectral SNR 曲線。**這張圖就是 07 章「投影 → CTF → 雜訊」的最佳教學對照圖**。
- 訊號離域公式：λδ/d（p.58 圖說）。300 keV 時 λ = 2 pm，δ = 2.2 µm，d = 3.5 Å → 離域約 120 Å，這是 box size 必須放大的原因。
- **Fig.2（p.61）**：SNR = 1 / 0.5 / 0.25 三個雜訊水準下的模擬顆粒影像與角度誤差等高線（1.7° / 3° / >10°）。直接支撐「SNR 決定可達解析度」的教學論點。
- **Fig.3（p.64）**：identifiability problem 的 1D/2D 示意（雙母體在高雜訊下與單母體無法區分；SNR 加倍後 PC1 顯現）。
- 複雜度公式 `n_c = t²π²n⁵m/s`（p.60），以及 n = 2D/d、角度精度 sin⁻¹(d/D)（p.61）。
- CTF 零點與翻轉次數的數值例（p.62），可直接拿來設計 07 章的 CTF 參數練習。

## 適用章節
- **06（主要）**：整篇即 SPR pipeline 的原理骨架——picking、CTF 估計、2D 分類、3D 重建、3D 分類、異質性。
- **07（主要）**：Fig.1 的影像形成模型、λδ/d 離域、CTF 零點數值、Fig.2 的 SNR 對照，全部可直接對應 ASPIRE 投影＋CTF＋雜訊的模擬設計。
- **03**：Fourier slice theorem 的最簡陳述與五步驟重建流程（p.59），可當 03 章傅立葉的 cryo-EM 應用收尾。
- **05**：p.65 關於分子量、劑量、輻射損傷、energy filter 的段落，可支援成像原理與「為何 cryo-EM 難」的敘述。
- 02（次要）：DoG particle picking 可當濾波/樣板比對的實例（p.62）。
