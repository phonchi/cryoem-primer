> 這是 Claude 寫作工作結束前的初稿來源紀錄，保留作協作證據。主代理其後已修正科學、術語與銜接問題，最終狀態見 `validation.md`、`science_late_fixes.md` 與 `reader_review.md`。其中 PMCID 衝突等提醒來自較早快照：Singer 2020 正確為 PMC8412055，Bendory 2020 為 PMC7213211，網站連結已核實。

# Claude 撰寫紀錄：SPA 七章來源與核算（2026-09-11）

本檔記錄我實際讀過的來源、每個數值範例的獨立核算，以及尚未解決的問題。學生端頁面不含任何佔位符或待辦標記；所有未確認事項都留在這裡。

## 已寫入的檔案（僅這七個 + 本檔）

| 檔案 | 狀態 |
|---|---|
| `book/05_background.md` | 重寫 |
| `book/05_image_formation.md` | 重寫 |
| `book/06_workflow.md` | 重寫 |
| `book/06_alignment_classification.md` | 新增 |
| `book/06_reconstruction_validation.md` | 重寫 |
| `book/06_heterogeneity.md` | 新增 |
| `book/06_resolution_validation.md` | 新增 |

未觸碰：01–04、07、`_toc.yml`、`_config.yml`、`references.bib`、`00_intro.md`、`08_resources.md`、`appendix_conventions.md`、任何 `.py`／`.ipynb`／`tests/`／`images/`（僅引用既有圖檔，未新增或修改資產）。

## 實際讀過的來源

全部取自本機 `References/`，多數透過既有的 `References/document_cache/text/*.txt` 抽取文字閱讀；未下載任何外部檔案，未使用付費服務。

| 書目鍵 | 快取檔（`References/document_cache/text/`） | 讀到的段落 |
|---|---|---|
| `passmore2016` | `2016-chapter-three-specimen-preparation-...-ed16b046.txt` | 目錄、第 1 節、第 2 節、2.3（頁 56–59）、2.4（頁 59–60） |
| `glaeser2016` | `2016-chapter-two-specimen-behavior-...-5b96ee6e.txt` | 目錄、第 1 節、2.2–2.3（頁 22–25）、Figure 2 說明、第 6 節標題 |
| `mcmullan2016` | `2016-chapter-one-direct-electron-detectors-...-6ba89b4a.txt` | DQE／MTF 定義（頁 2–3）、偵測器效能比較（頁 8–9） |
| `rubinstein2016` | `2016-chapter-five-processing-of-cryo-em-movie-data-...-a351adea.txt` | 目錄、第 1 節、第 2–6 節的演算法描述與加權式 |
| `cheng2016` | `2016-chapter-four-strategies-...-845e6f5a.txt` | 首頁 DOI 與章節定位 |
| `scheres2016` | `2016-chapter-six-processing-of-structurally-heterogeneous-...-0d35407a.txt` | 2.2、2.3（頁 129–133）、頁 139、頁 142–143、4.4–4.5（頁 144–146） |
| `rosenthal2016` | `2016-chapter-nine-testing-the-validity-...-dfd4fbce.txt` | 3.2–3.3、第 4 節（頁 244–246） |
| `penczek2010` | `2010-chapter-one-fundamentals-...-fe0a4ebc.txt` | 目錄、第 1 節（頁 2–5）、第 2 節 central section theorem、第 3 節、第 5 節、第 7 節 gridding／Voronoi |
| `penczek2010resolution` | `2010-chapter-three-resolution-measures-...-9dffd482.txt` | 第 3 節（FRC／DPR／SSNR 定義、兩種 half-set 建法）、第 4 節（式 3.14／3.15、門檻、Wiener 增益式 3.20–3.22） |
| `scheres2005` | `2005-ml2d-pdf-4245e350.txt` | Abstract、Introduction（頁 140）、Mathematical Background（頁 140–142） |
| `sorzano2010` | `2010-cl2d-pdf-496ae98b.txt` | 第 2 節、2.1 correntropy、2.2 divisive VQ、2.3 穩健準則（頁 198–200） |
| `yang2012` | `2012-isac-pdf-aad8de60.txt` | SUMMARY、Introduction、EQK-means／stability／reproducibility／pass 流程 |
| `chen2014` | `2014-gsup-pdf-6114cc27.txt` | Abstract、第 1 節、節次結構（2、2.1–2.3、3、3.1–3.3、4.1–4.2、5、6） |
| `scheres2012bayes` | `2012-relion-pdf-a3756c86.txt` | 頁 407–411：Wiener filter 寫法（式 3、4）、Bayesian 觀點、MAP／EM 更新（式 9–12） |
| `punjani2017` | `2017-cryosparc-pdf-597d0f73.txt` | 頁 290–291：式 (1) 的 Bayesian 目標、SGD 說明、Figure 1 圖說、branch-and-bound 段落 |
| `punjani2021` | `2021-3dva-pdf-3154b03b.txt` | 第 1 節、第 2 節（式 1、2 與 PCA 的困難）、第 3 節前段、5.2 節（式 4、6–10） |
| `sorzano2022` | `2022-on-bias-var-overfitting-spa-pdf-52ec4cbf.txt` | 第 2 節（MSE 分解、bias–variance）、第 4 節片段、5.2、5.3、第 6 節 |
| `sigworth2016` | `2016-principles-of-cryo-em-single-particle-image-pdf-38aaccdc.txt` | 頁 57–58、59–60（複雜度公式）、60–61（Figure 2 角度誤差）、61–62（model bias）、62–63（picking、CTF）、64–65（Figure 3 可識別性、變異／共變異數、SNR 限制與劑量） |
| `singer2020` | `2020-math-review-pdf-53188edd.txt` | 節次結構全覽、3.2（頁 172）、3.3（頁 173）、第 5 節（頁 182–183） |

僅由既有書目資料引用、未新讀全文者：`scheres2010`、`scheres2010heterogeneity`、`leschziner2010heterogeneity`、`penczek2010restoration`、`cheng2015primer`、`sigworth1998`。這些條目的延伸閱讀說明只寫到我能從書目與章節標題確認的層級（章名、頁碼區間、主題），沒有宣稱具體段落內容。

DOI 全部取自各快取檔首頁的 DOI 行，未憑記憶填寫：
`bs.mie.2016.04.008 / .009 / .010 / .011 / .012`、`bs.mie.2016.05.056`、`bs.mie.2016.06.004`、
`S0076-6879(10)82001-4 / 82002-6 / 82003-8 / 82011-7 / 82012-9 / 82013-0`。
期刊論文的 DOI 沿用 `references.bib` 既有欄位。

## 數值與公式的獨立核算

每一項都自己重算過，下面列出算式與結果。

1. **粒子密度**（`05_background`）：$2\ \mathrm{g/L}\times 8\times10^{-17}\ \mathrm L=1.6\times10^{-16}\ \mathrm g$；$\div10^{6}\ \mathrm{g/mol}\times6.022\times10^{23}=96.4$ → 約 100 particles/μm²；250 kDa 為四倍 → 386 ≈ 400。與來源給的「約 100」「約 400」一致。1.2 μm 洞的面積 $\pi(0.6)^2=1.131\ \mathrm{\mu m^2}$ → 約 113 顆（文中寫「約 110」）。
2. **輻射劑量**：$2.4\ \mathrm{MeV\,cm^2/g}\times10^{17}\ \mathrm{cm^{-2}}=2.4\times10^{17}\ \mathrm{MeV/g}$；$\times1.602\times10^{-13}\ \mathrm{J/MeV}=3.845\times10^{4}\ \mathrm{J/g}=3.845\times10^{7}\ \mathrm{Gy}=3.845\times10^{9}\ \mathrm{rad}$。與 `glaeser2016` 的 $3.8\times10^9$ rad 相符。
3. **相關雜訊下的平均變異**：$\operatorname{Var}(\bar\varepsilon)=[1+(N-1)\rho]\sigma^2/N$。$N=10^4$、$\rho=0.01$ → $0.010099\sigma^2$；$N_{\mathrm{eff}}=10^4/100.99=99.02$；$N\to\infty$ 極限 $1/\rho=100$。
4. **電子波長**：以 $\lambda=12.2639/\sqrt{V+0.97845\times10^{-6}V^2}$（$V$ 為伏特）計算，200 kV → 0.025079 Å（與 `07_synthetic_data` 測試的 0.02508 一致），300 kV → 0.019687 Å。
5. **CTF 第一零點**：$d_1=\sqrt{\Delta f\lambda}=\sqrt{10^4\times0.0197}=14.04\ \mathrm{\mathring A}$，與 `sigworth2016` 頁 62 的「1 μm 離焦、300 keV，第一零點在 14 Å」相符。反轉次數 $\Delta f\lambda s^2$：7 Å → 4.02；3.5 Å → 16.08，與同頁的 4 與 16 相符。
6. **$C_s$ 的高頻貢獻**：$\tfrac12C_s\lambda^3s^4$，取 $C_s=2.7\times10^7\ \mathrm{\mathring A}$、$\lambda=0.0197$、$s=1/3.5$ → 0.688（以 $\pi$ 為單位）。
7. **離焦容許誤差**：$\Delta\chi=-\pi\delta\lambda s^2$，$|\Delta\chi|=\pi$ → $\delta=d^2/\lambda$；$d=3.5\ \mathrm{\mathring A}$、300 kV → 622 Å ≈ 62 nm。**注意**：`sigworth2016` 頁 62 寫「a defocus error of 124 nm yields a complete reversal in the polarity of the CTF at 3.5 Å」，恰為我算出值的兩倍。在本書採用的 $\chi=2\pi(-\tfrac12\Delta f\lambda s^2+\cdots)$ 慣例下，$\pi$ 相移對應 62 nm。因此章內只寫我自己推導的式子與 62 nm，並未引用 124 nm 這個數字。差異可能來自該文對「complete reversal」的定義或不同的 $\chi$ 慣例，尚未查證。
8. **Delocalization**：$\lambda\Delta f/d=0.0197\times2.2\times10^4/3.5=123.8\ \mathrm{\mathring A}$，與 `sigworth2016` 頁 57 圖說的「about 120 Å」相符。
9. **Poisson 計數的 SNR**：$\mathrm{SNR}=(cn)^2/n=c^2n$；$c=0.05$、$n=40$ → 0.1。DQE 換算：$0.5/0.08=6.25$ 倍粒子數。
10. **Whitening 的頻帶例子**：未加權時低頻雜訊功率 100 vs 高頻 1；whitening 後每係數證據為 $10/100=0.1$ 與 $0.5/1=0.5$，相差五倍。
11. **對位誤差包絡**：二維等向高斯位移的特徵函數 $\exp(-2\pi^2\delta_0^2|\mathbf k|^2)$。$\delta_0=1$ px：$|\mathbf k|=1/6$ → $\exp(-0.5483)=0.578$；$|\mathbf k|=1/3$ → $\exp(-2.193)=0.1115$。等效 $B=8\pi^2\delta_0^2$，$\delta_0=1\ \mathrm{\mathring A}$ → 78.96 Å²。
12. **角度誤差換算**：$1.7^\circ=0.02967$ rad，$\times65$ px $=1.93$ px；代回包絡在 $|\mathbf k|=1/6$ 得 $\exp(-19.739\times3.725/36)=\exp(-2.043)=0.130$。
13. **大類別吸引效應**：$\mathbb E\|y-a_k\|^2=\|s-s_k\|^2+d\sigma^2+d\sigma^2/N_k$。$d=16{,}900$、$\sigma^2=1$：$N=1000$ → 16.9；$N=50$ → 338；門檻差 321.1。
14. **CTF 加權合併**：$\sum H_j^2=0.81+0.04+0.36=1.21$。無先驗：$\widehat V=1$、$\operatorname{Var}=1/1.21=0.826$。$\tau^2=1$：分母 2.21，$\mathbb E[\widehat V]=0.5475$，$\operatorname{Var}=1.21/4.8841=0.2477$，$\mathrm{bias}^2=0.2047$，MSE $=0.4525$。CTF 零點附近 $H=0.02$：無先驗變異數 $1/0.0004=2500$。
15. **手性歧義**：取 $M$ 為鏡射（$\det M=-1$、$M^{-1}=M$）、$S=\mathrm{diag}(1,1,-1)$、$R'=MRS$，則 $\det R'=1$ 且 $P_{R'}V'(\mathbf u)=P_RV(\mathbf u)$（換 $z'=-z$）。推導自行完成，未引用來源。
16. **Fourier slice theorem**：在 $P_RV(\mathbf u)=\int V(R[u_x,u_y,z]^{\mathsf T})\mathrm dz$ 的定義下，切片為 $\widehat V(R[k_x,k_y,0]^{\mathsf T})$，**沒有轉置**。舊版 `06_reconstruction_validation.md` 寫成 $R^{\mathsf T}[k_x,k_y,0]^{\mathsf T}$，與凍結的 `07_synthetic_data` 慣例不一致，本次已改正並附推導。
17. **兩構形平均的抵消**：$\lvert\beta+(1-\beta)e^{-\mathrm i\varphi}\rvert$，$\beta=\tfrac12$ 時 $=\lvert\cos(\varphi/2)\rvert$。$d=5\ \mathrm{\mathring A}$：30 Å → 0.866；20 Å → 0.707；10 Å → 0。$\varphi=\pi$、$\beta=0.8$ → $\lvert2\beta-1\rvert=0.6$。連續位移：$\sigma_d=2\ \mathrm{\mathring A}$，10 Å → 0.454，5 Å → 0.0425。
18. **FSC 與 SSNR**：$\rho=P_S/(P_S+P_N)$ → $\mathrm{SSNR}_{\mathrm{half}}=\rho/(1-\rho)$；完整 map 雜訊功率減半 → $\mathrm{SSNR}_{\mathrm{full}}=2\rho/(1-\rho)$。與 `penczek2010resolution` 式 (3.14)／(3.15) 一致。門檻核對：$\rho=1/3\Rightarrow\mathrm{SSNR}_{\mathrm{full}}=1$；$\rho=0.5\Rightarrow2$（與該章列出的對照相同）。
19. **0.143 的來源**：$C_{\mathrm{ref}}=\sqrt{2\rho/(1+\rho)}$，$\rho=0.143$ → $0.286/1.143=0.25022$，開根號 $=0.5002$；$\cos60^\circ=0.5$。與 `rosenthal2016` 頁 245 的敘述一致。
20. **Wiener 增益**：$\mathrm{SSNR}_{\mathrm{full}}/(\mathrm{SSNR}_{\mathrm{full}}+1)=2\rho/(1+\rho)=C_{\mathrm{ref}}^2$。與 `penczek2010resolution` 式 (3.20)–(3.22) 一致（該處 OCR 的分母被截斷，已自行補回並驗證）。
21. **偏差－變異數字**：$b=0.2$、$\sigma^2=100$：$N=100$ → $0.04+1=1.04$；$N=10^4$ → $0.05$；$N=10^6$ → $0.0401$。
22. **銳化倍率**：$\exp(-Bs^2/4)$，$B=-100\ \mathrm{\mathring A^2}$：$s=0.1$ → $\exp(0.25)=1.284$；$s=1/3$ → $\exp(2.778)=16.09$。
23. **FSC 圖的讀數**：$1/0.267=3.745\ \mathrm{\mathring A}$（沿用原圖說的 3.74 Å）。
24. **Motioncorr 方程式數**：$T=30$ → $T(T-1)/2=435$ 條方程式、29 個未知逐 frame 位移。

## 尚未解決或需要他人確認的事項

1. **`tests/test_site_sources.py::test_understanding_checks_stay_focused_and_collapsible` 會失敗。** 該測試的 `expected` 字典仍要求 `05_background.md`、`05_image_formation.md`、`06_workflow.md`、`06_reconstruction_validation.md` 各含 3／4／4／4 題「理解檢查」。本次任務明確要求移除問答段落並把答案內容併入正文，因此這四個項目需要從字典中移除（保留 `00_intro.md`、01–04、`07_synthetic_data.py`、`appendix_conventions.md`）。測試檔不在我的可寫清單內，請由擁有 `tests/` 的一方調整。
2. **`_toc.yml` 尚未列入三個新章。** 目前 `parts` 只有 `05_background`、`05_image_formation`、`06_workflow`、`06_reconstruction_validation`。需要補上 `05_statistical_inference`、`06_alignment_classification`、`06_dimension_reduction`、`06_heterogeneity`、`06_resolution_validation`，順序依 `implementation.md`。`_toc.yml` 不在我的可寫清單內。
3. **`singer2020` 的開放全文連結在既有頁面中不一致**：`05_statistical_inference.md` 與 `06_dimension_reduction.md` 用 `PMC7213211`，`08_resources.md` 用 `PMC8412055`，而 `08_resources.md` 又把 `PMC7213211` 指給 `bendory2020`。兩者至多有一個正確。我在自己的七章裡一律只用 DOI，避免傳播未經確認的 PMCID。請文獻代理確認後統一。
4. **`sigworth2016` 的 124 nm 離焦容許值**（見上面第 7 點）與本書慣例下的推導相差兩倍，尚未釐清；章內未引用該數字。
5. **`chen2014` 的頁碼定位**：手上的快取是 IMS 電子重印本，頁碼自 1 起算（期刊頁碼為 259–285）。為避免換算錯誤，延伸閱讀只寫節次（第 2、3、4.2 節），未標頁碼。
6. **圖檔用途未經原作者確認**：我為新章引用了先前未使用的既有圖檔 `s24_1`（球面角度）、`s12_6`（圓形 mask 類別平均）、`s13_4`（數個 3D 類別）、`s15_8`（破碎拉長的密度）、`s07_3`（挑選座標疊圖）、`s08_3`（徑向 power spectrum 與 CTF 擬合），並把 `s08_1`（`fig-ctf-theory`）由 `06_workflow` 移到 `05_image_formation`、`s24_2`（`fig-fsc-unmasked`）由 `06_reconstruction_validation` 移到 `06_resolution_validation`。圖說都寫成不宣稱來源資料細節的敘述。若原始投影片對這些圖有更精確的說明，圖說可再收緊。未修改任何圖檔本身。
7. **未執行任何測試或建置**：本次工作環境沒有 shell 工具，因此 `pytest`、Jupyter 建置與 `pdftotext` 都沒有執行；所有來源均透過既有的文字快取閱讀，所有檢查以人工核算與字串比對完成。與測試相關的字串（`05_image_formation.md` 的 `\frac{\mathrm{SNR}_{\mathrm{out}}(s)}{\mathrm{SNR}_{\mathrm{in}}(s)}`、`06_reconstruction_validation.md` 的「FSC 衡量兩張 half-map 的一致性」與「共同的正比例常數，FSC 仍為 1」等）已逐一以搜尋確認存在，禁用詞彙與 `不是…而是` 一類句式也已全檔掃描。

## 保留的錨點

| 檔案 | 錨點 |
|---|---|
| `05_background.md` | `low-dose-snr` |
| `05_image_formation.md` | `detectors`、`dose-weighting` |
| `06_workflow.md` | `workflow`、`motion-correction`（新增）、`refinement`、`heterogeneity` |
| `06_reconstruction_validation.md` | `fourier-coverage`、`reconstruction`、`soft-assignment`、`validation` |
| `06_alignment_classification.md` | `rigid-alignment`、`alignment-classification-coupling`、`class-stability` |
| `06_heterogeneity.md` | `discrete-heterogeneity`、`linear-heterogeneity`、`confounding` |
| `06_resolution_validation.md` | `halfset`、`fsc`、`sharpening` |

`refinement`、`heterogeneity` 與 `validation` 三個舊錨點的位置都保留了具實質內容的橋接段落，並以 `{doc}` 連到對應的新章。
