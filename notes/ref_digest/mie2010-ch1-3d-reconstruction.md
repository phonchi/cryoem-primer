# [2010] Fundamentals of Three-Dimensional Reconstruction from Projections（Methods in Enzymology Vol. 482, Ch.1）

- 作者：Pawel A. Penczek（University of Texas, Houston）
- 出處：*Methods in Enzymology*, 2010, Vol. 482, pp. 1–33，doi:10.1016/S0076-6879(10)82001-4
- 頁碼對照：印刷頁 1–33 ＝ PDF 頁 1–33（`印刷頁 = PDF 頁`，位移為 0）
- 原檔：`References/2010/Chapter-One---Fundamentals-of-Three-Dimensional-Reconstr_2010_Methods-in-Enz.pdf`

## 一句話定位
3D 重建演算法的權威教科書式導論：從 ray transform、central section theorem 到 ART/SIRT、filtered backprojection、gridding 直接傅立葉反演一路推到底，是 06 章「3D 重建」階段與 03 章傅立葉應用的主要引用來源。

## 關鍵主張與內容
- 離散的投影反演本質上是 **ill-posed 反問題**：投影矩陣 P 為長方形，故必然過定或欠定；且解不一定連續（g 相近不保證重建結果相近）（p.3–4）。
- **Ghost 現象**：對任一組有限的投影方向，都存在非零物體 d̂₀ 使 P d̂₀ = 0——即**離散投影轉換即使無雜訊也不可逆**。實務上 ghost 屬極高頻成分，可用正則化壓下（p.4）。
- EM 重建與一般斷層掃描的五個差異（p.4–5）：(1) SPR 與 double-tilt ET 是罕見的「真 3D」問題，多數 CT 可拆成一系列 2D 切片重建；(2) SPR 的取樣幾何無法控制、投影分布隨機；(3) 分布極不均勻造成 SSNR 在物體內不均，且資料 SNR 極低；(4) 取向參數只是近似值，誤差兼有隨機與系統性，故重建演算法必須放在 refinement 流程中一起看；(5) 投影數遠多於影像邊長（嚴重過取樣），造成極大的計算負擔。
- **Ray transform vs Radon transform**：ray transform 沿直線積分，把 nD 函數變成 (n−1)D；Radon transform 是對 (n−1)D 超平面積分。2D 時兩者只差記號，3D 時不同——cryo-EM 處理的是 3D ray transform 的反演（p.5）。
- 尤拉角採 **ZYZ 慣例**（SPIDER、IMAGIC、MRC、FREALIGN、EMAN2、SPARX 都用）。角 ψ 為平面內旋轉，不影響 ray transform 的資訊量。φ、θ、ψ 加上平移 t_x, t_y 合稱 **projection orientation parameters**（p.5–6）。
- **Orlov's condition**：ray transform 可逆的最小投影分布，是單位向量沿任一條連接球面對徑點的連續曲線移動所掃出的集合。ET 的 single-axis tilt（全角度範圍）滿足此條件；SPR 因取向隨機不可控，Orlov 條件的實務意義有限。SPR 反而希望投影方向**盡量均勻**，以獲得均勻的 SSNR 分布（p.7–8）。
- 重建方法二分法（p.8–9）：
  - **代數法**（algebraic）：先離散化再解線性方程組。優點是可直接用線性代數工具、離散化/插值自然內建於 P、不需為特定投影分布設計加權；缺點是幾乎必為迭代法、計算量極大、迭代次數難以事先決定。
  - **轉換法**（transform）：立基於 central section theorem。又分兩支——(a) 在 Fourier 空間設計線性濾波器後在實空間做簡單 backprojection（即 FBP，好寫、夠快，但不均勻分布時濾波函數只是粗略近似）；(b) 完全在 Fourier 空間操作的 **direct Fourier inversion**（實作困難但效率極高，其中 gridding 演算法目前公認最準確）。
- 取樣與插值（p.10–12）：像素大小 a ≤ 1/(2 s_max)；實務經驗法則是**取像素大小約為目標解析度的 1/3**，以壓低插值的不良影響。Backprojection 可用 voxel-driven 或 ray-driven 實作，兩者插值假影不同，迭代演算法中必須保持兩步驟的矩陣互為轉置。兩者複雜度皆為 O(K³)。
- **插值造成的高頻損失（03 章可直接用的量化結論）**：nearest-neighbor 插值等價於與寬 1 pixel 的矩形窗卷積，頻域即乘上 sinc(s)；bilinear 插值等價於與三角函數卷積，頻域為 sinc²(s)。在 Nyquist 頻率處，NN 損失 36%、bilinear 損失 60%。常見的三倍過取樣可把 bilinear 的高頻壓抑降到 10%，但 2D 資料量與計算時間都會變成 9 倍；gridding 法則幾乎在全頻段完美且不需過取樣（p.11–13, Fig.1.3）。
- **SIRT**：以梯度下降最小化 |Pd − g|²，全部體素同時更新；在 EM 各種條件下表現優異，且投影分布有角度缺口時假影最不擾人。可用 steepest descent 或 conjugate gradient 加速（CG 約 10 次迭代收斂，但須加正則項 |Bd|²，B 為 Laplacian 或高階導數的離散近似）（p.14–15）。
- **ART（= Kaczmarz 法）**：每次只針對單一投影的單一像素做修正並立即更新結構，因此收斂遠快於 SIRT。若選取順序使相鄰兩像素的投影方向互相垂直可再加速，隨機順序效果也幾乎一樣好。Xmipp 的「ART with blobs」以 Kaiser–Bessel 球對稱基底函數表示結構，大幅減少所需迭代數（p.15–16）。
- 對代數迭代法的重要警語：**ART 與 SIRT 都不會用有意義的資訊填補 Fourier 空間的缺口**；更激進的正則化（如 Maximum Entropy）能減少假影但不會增加資訊。文中直指此領域曾有「誇大的宣稱與過度美化的結果」。根本原因之一是取向參數的誤差使方程組 (1.9) 本身不自洽（p.16–17）。
- 迭代法可用 a priori 知識做正則化（與量測資料相似、密度為正、有界支撐等），形式上以 **POCS**（projections onto convex sets）描述（p.17）。
- **FBP** 三步驟（2D 版）：對每個投影做 1D FFT → 乘上 ramp filter |R| → 逆 1D FFT → 實空間 backprojection。3D 的加權函數設計才是難題：Radermacher 的 general weighting function（GWF）、Harauz–van Heel 的 "exact filter"（以 Fourier 中心切片的幾何重疊倒數為權重）。三種加權函數都依賴物體半徑 r_max，可藉調整 r_max 針對特定取樣幾何最佳化。**加權函數要算兩兩投影間的距離，複雜度隨投影數平方成長，大資料集下效率不佳**（p.18–21）。
- **BFP（backprojection-filtering）不能取代 FBP**：在離散形式下，Jacobian 權重的套用順序不可反轉——「加權投影之和」不等於「投影之和再加權」。均勻分布時 BFP 尚可，所有有趣的 3D 非均勻情形都不行（p.19）。
- **GDFR（gridding-based direct Fourier reconstruction）四步驟**（p.22–24, Fig.1.4）：(1) reverse gridding 把 2D 投影重取樣到 2D 極座標；(2) gridding——以 2D 單位球面 Voronoi 圖的胞元面積作為 gridding weight（避開昂貴的 3D Voronoi），與卷積核 F[w] 做卷積累積到 Cartesian 格點；(3) 3D 逆 FFT；(4) 實空間除以 w 移除權重。建議窗函數為可分離的 Kaiser–Bessel 窗、Fourier 支撐設為 6 個 voxel。**GDFR 在無雜訊測試中全頻段幾乎完美，優於 SIRT 再優於 weighted backprojection；且比 weighted backprojection 快 6 倍、比 SIRT 快 17 倍。**
- 各軟體的 CTF 處理策略（p.24–27）：SPIDER 把資料依 defocus 分組、各組獨立重建，最後用 Wiener filter 合併（低顆粒產率的資料無法處理）；IMAGIC/EMAN 先在 2D 做含 CTF 校正的 Wiener filter 得類平均再重建（但類平均的 SSNR 不均勻並未被後續考慮）；FREALIGN 用改良的 trilinear 插值直接實作式 (1.27)；SPARX 用 NN direct inversion（先補零到 4 倍、2D FFT、以 NN 插值累積並套 Wiener filter，最後套上 Bracewell 式的「local density」加權）。**含 CTF 校正的迭代重建（式 1.14）是目前唯一理論上乾淨的做法，但只在 SPIDER 實作且未普及。**
- 歷史脈絡（可作 05/06 章的敘事素材）：1956 Bracewell 由地球多角度視圖重建太陽黑子並提出 FBP；1968 DeRosier & Klug 利用二十面體對稱重建病毒衣殼；1969 Hoppe 提出非週期樣品的 3D 高解析 EM；1972 Hounsfield & Cormack 的 CAT。演算法重心的演變：1980 年代 FBP（general / exact filter）→ 1990 年代初 cryo 製備帶來低對比資料，轉向 ART/SIRT → 1990 年代末資料量暴增（70S ribosome 從 303 張影像變成超過 500,000 張），轉向快速的 direct Fourier inversion；目前除 IMAGIC 外幾乎所有主流套件都實作了某種版本（p.2–3, 27–28）。

## 可引用的公式/圖表
- **式 (1.1)（p.5）**：ray transform `g(x_t) = ∫ d(r) dt`，t ⊥ x。
- **式 (1.2)（p.6）**：ZYZ 尤拉角的三個旋轉矩陣乘積（φ, θ, ψ）。
- **式 (1.3)（p.6）——central section theorem**：`F[g(x_t)] = D(s_t)`，s_t ⊥ t。**這是 03 章與 06 章最該引用的一式。**
- 式 (1.4)–(1.6)（p.11）：Shannon 取樣定理式的 Fourier 空間插值與 Crowther 的最小平方解。
- **式 (1.7)（p.12）**：`sinc(s) = sin(πs)/(πs)`——NN 插值的頻域效應；bilinear 對應 sinc²(s)。
- **式 (1.8)–(1.10)（p.13–14）**：離散化為 `g = P d`，最小平方解 `d̂ = (PᵀP)⁻¹Pᵀg`。
- **式 (1.11)–(1.12)（p.14）**：SIRT 目標函數與更新式 `d_{i+1} = d_i − λ_i(PᵀP d_i − Pᵀg)`。
- **式 (1.13)–(1.14)（p.15）**：含 CTF 與正則化的目標 `L(d) = (1−α)|SPd − g|² + α|Bd|²` 及其迭代式；S 為 PSF 的代數表示。
- **式 (1.15)–(1.16)（p.15–16）**：ART / Kaczmarz 更新式。
- **式 (1.17)（p.18）**：FBP 反演式，`d(x) = Backprojection[Filtration_{|R|}(g)]`（ramp filter |R|）。
- 式 (1.18)–(1.19)（p.18–19）：Riemann-sum 加權 `c(R_j, ψ_j) = R_j Δψ_j / 2π` 與取樣充分性條件 `Δψ_j ≤ 1/r_max`。
- 式 (1.20)（p.19）Bracewell 的 local-density 加權；式 (1.21)（p.19）BFP；式 (1.22)（p.20）Radermacher GWF；式 (1.23)–(1.26)（p.20–21）Harauz–van Heel 的 exact filter。
- **式 (1.27)（p.25）——Wiener filter 合併多 defocus 資料**：`D = Σ_n CTF_n·SSNR_n·G_n / (Σ_n CTF_n²·SSNR_n + 1)`。**06 章解釋「為何要收集不同 defocus」的關鍵式。**
- 式 (1.28)–(1.29)（p.26）SPARX 的 local-density 加權函數（參數 β = 0.2、n = 3）。
- **Fig.1.1（p.6）**：投影球面與 g(x_t) 的幾何示意（含尤拉角定義）。
- **Fig.1.2（p.7）**：真實資料（Thermus thermophilus ribosome + EF-Tu，6.4 Å，322,688 個角度）的投影方向半球分布圖，色階表示各方向的顆粒數。**06 章講 angular distribution / preferred orientation 的最佳實例圖。**
- **Fig.1.3（p.13）**：把 1024² 高斯白雜訊影像旋轉 45° 後，旋轉平均功率譜的比值曲線——bilinear（最差）、quadratic（中）、gridding（幾乎為 1）。**02/03 章講插值品質時直接可用的量化圖。**
- **Fig.1.4（p.23）**：GDFR 五格流程圖（2D FFT → reverse gridding 到極座標 → 球面 Voronoi 權重 → gridding + 3D 逆 FFT → 實空間除權重）。

## 適用章節
- **06（主要）**：3D 重建階段的完整方法學——ART/SIRT vs FBP vs direct Fourier inversion 的取捨、各主流軟體實作、CTF 校正策略、angular distribution 的影響（Fig.1.2）、式 (1.27) 的多 defocus 合併。
- **03（主要）**：central section theorem（式 1.3）、ramp filter 與 FBP 推導（式 1.17）、Fourier 空間插值、gridding/NUFFT、sinc 與 sinc² 的插值頻域效應。
- **07（次要）**：ray transform 的定義（式 1.1）與 ZYZ 尤拉角慣例（式 1.2）——ASPIRE 產生投影時的座標與角度約定；像素大小取解析度 1/3 的經驗法則。
- **02（次要）**：插值方法（NN / bilinear / quadratic / gridding）的品質比較與 Fig.1.3，可作為重取樣與旋轉操作的延伸閱讀。
- **05（次要）**：p.2–3 與 p.27–28 的歷史脈絡（Bracewell 1956 → DeRosier & Klug 1968 → Hoppe 1969 → CAT 1972；資料集規模從 303 張到 50 萬張的演變）。
