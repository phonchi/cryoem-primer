# [2016] Processing of Cryo-EM Movie Data（Methods in Enzymology Vol. 579, Ch.5）

- 作者：Z.A. Ripstein, J.L. Rubinstein（Hospital for Sick Children / University of Toronto）
- 出處：*Methods in Enzymology*, 2016, Vol. 579, pp. 103–124，doi:10.1016/bs.mie.2016.04.009
- 頁碼對照：印刷頁 103–124 ＝ PDF 頁 1–22（`PDF 頁 = 印刷頁 − 102`）
- 原檔：`References/2016/Chapter-Five---Processing-of-Cryo-EM-Movie-Data_2016_Methods-in-Enzymology.pdf`

## 一句話定位
以統一數學記號系統性比較五種 motion correction 演算法（Motioncorr、alignframes_lmbfgs/alignparts_lmbfgs、Unblur、Optical Flow、RELION particle polishing）的權威綜述，是 06 章「movie processing / beam-induced motion correction」段落的核心引用來源。

## 關鍵主張與內容
- **DDD 相機改善影像 SNR 的四種機制**（p.104）：(1) monolithic active pixel sensor 的 DQE 優於底片/CCD；(2) CMOS 高幀率可讓電子計數（electron counting），進一步提升 DQE；(3) 高幀率允許把單次曝光拆成多幀 movie，藉由對齊各幀來部分校正 beam-induced movement；(4) 可用低劑量記錄高解析特徵、高劑量記錄低解析特徵（後者對輻射損傷較耐受）。
- **三家商用 DDD 相機比較**（p.104–105）：Direct Electron LP 的 DE 系列（6.0–6.5 μm 像素，30–60 fps）；FEI Falcon（14 μm 像素，高 DQE 但僅 18 fps，尚不適合電子計數）；Gatan K2（5 μm 像素，內部 400 fps 電子計數、可記錄至 40 fps）。不同相機因像素尺寸與幀率不同，最佳使用策略也不同：積分型相機（DE/Falcon）宜用短且高幀率的 movie 減少載台漂移；K2 因電子計數需限制曝光速率避免 coincidence loss，故載台漂移在單一 movie 內可能較顯著。
- **大顆粒 vs 小顆粒的處理策略差異**（p.105）：大型病毒顆粒因訊號足夠，可對每個顆粒單獨的 movie frame 做剛體平移+旋轉校正；小顆粒訊號不足以單獨追蹤，需要利用整個 movie 的資訊（局部鄰近顆粒的相關性）才能估計軌跡。
- **本文所有方法的共同限制**（p.105）：目前的 whole-frame 與 individual-particle motion correction 方法都**忽略樣品的面外運動（out-of-plane motion）與旋轉**，將樣品/局部樣品在各幀之間的運動簡化為平移。
- **Motioncorr**（第 2 節，p.107–109）：對 movie 中所有幀兩兩計算 Fourier 空間互相關函數（T 幀共得 T(T−1)/2 個唯一幀對），每對互相關峰值給出位移估計向量 m̄_mn；用固定圖樣噪聲抑制濾波 `exp(−B/4d²)`（B 為溫度因子、d 為解析度）預處理影像；因互相關向量數遠多於待求的 T−1 個真實幀間位移 s̄_mn，構成 overdetermined 線性方程組（式 1），以最小平方法求解。可用峰值次像素定位（Fourier padding 插值）與剔除殘差過大的方程來改善精度；原始實作需 GPU，且**不適合對小於 2000×2000 像素的影像區域做對齊**（即無法直接用於個別顆粒層級的形變校正）。
- **alignframes_lmbfgs 與 alignparts_lmbfgs**（第 3 節，p.109–112）：把對齊問題表述為一個可解析求導的目標函數最佳化問題，用 Fourier 空間相位偏移（式 2：相位變化 `φ_jt = k_x(j)·x_t·2π/N + k_y(j)·y_t·2π/N`）取代影像域內插；目標函數（式 3）為所有已位移幀與位移幀總和的未正規化互相關之和的負值，可解析求梯度（式 4–5），以 L-BFGS 最佳化——這是程式名稱由來。額外加入**二階平滑懲罰項**（式 6–8，懲罰軌跡的二階差分，即防止突然轉向/加速），視為使用者可調參數 λ。alignparts_lmbfgs 額外對個別顆粒軌跡做**高斯加權局部平均**（式 9–10：權重 `w_mn = exp(−d²_mn/2σ²)`，d 為顆粒間距），並依 icosahedral 病毒顆粒研究得出的 critical exposure 曲線做輻射損傷加權。alignparts_lmbfgs 建議只用於**未對齊的原始 movie**（因二階平滑約束若疊加在已做過 whole-frame 對齊的資料上，可能扭曲真實顆粒軌跡）。
- **Unblur**（第 4 節，p.112–114）：同樣用互相關法，但不建立單一聯合目標函數，而是**反覆將每一幀對齊到「不含該幀」的其餘幀平均影像**（正規化互相關，式 11），迭代直到位移收斂（如低於 0.5 像素閾值）。以三次雲線（cubic spline）擬合位移軌跡以平滑掉雜訊造成的假轉向，等價於最小化含平滑懲罰項的泛函（式 12），且**以交叉驗證自動決定平滑權重 λ**（不需使用者手動設定）。同時對疊加平均做**曝光加權（exposure weighting）**以補償輻射損傷，權重來自直接測量的 critical exposure 曲線。
- **Optical Flow**（第 5 節，p.114–116）：Abrishami et al. 方法不用互相關，而用像素強度差異的一階泰勒展開估計位移（式 13–14b：`I(x+Δx,y+Δy,t+Δt) ≈ I(x,y,t) + (∂I/∂x)Δx + (∂I/∂y)Δy + (∂I/∂t)Δt`，假設總強度不變）；每個像素窗口都給出局部超定線性方程組（形式同 Motioncorr 的 A·s̄=B），求解得到整張影像的向量場 OF。因單步不夠穩健，採用**金字塔式疊代方案**（p.115–116）：先平均全部幀，再把 movie 拆成 2^k 組分別平均並與全域平均算 optical flow（式 15），逐層細分直到每組只剩一幀。
- **RELION Particle Polishing**（第 6 節，p.116–120）：與其他方法根本不同之處在於**於 3D map 精修「之後」執行**（其他方法都在初始 2D 影像分析階段之前），利用 3D 參考圖的 2D 投影而非原始影像資訊來對齊個別顆粒。流程：先用 whole-frame 對齊後的平均影像做 3D refinement 得到每顆粒的 Euler angle；固定角度後，用 map projection 與各幀（常為 running average）顆粒影像算互相關估計每幀的平移 −x̄_t；假設軌跡為**線性、等速**（式 16：對每顆粒 n 最小化 `O(Θ)_n = Σ_m w_m Σ_t (x̄_mt − (ᾱ_n + β̄_n t))²`，w_mn 同 alignparts_lmbfgs 的高斯鄰近權重），只估計 x、y 平移（不含旋轉）。另外對每幀做**資料驅動的頻率相關加權**：對每幀單獨算 3D map 並與整體平均 map 的 FSC 比較，得相對溫度因子 B_t 與截距 C_t（式 17），再算加權函數 `w_t(1/d) = exp(B_t/4d² + C_t) / Σ_t' w_t'(1/d)`（式 18，正規化使加權平均總功率不變、不做額外銳化）。**注意此處相對 B-factor 正負號與一般包絡函數相反**：正的相對 B-factor 代表該幀的高頻訊號優於平均、應被放大。實務上前幾幀常因初始劇烈運動被降權，之後曲線形狀與輻射損傷測量的降權曲線（Unblur/alignparts_lmbfgs 所用）相似（Fig.3B）。
- **各方法優劣總結**（第 7 節，p.121–122）：Motioncorr 的互相關高度超定但雜訊幀對雜訊幀比對雜訊大、計算量大（需 GPU）；alignframes_lmbfgs/Unblur/Optical Flow 都改成與（更高 SNR 的）平均影像比對，對低對比影像可能表現較好；Unblur 是三者中最簡單（重複對齊到平均值）；Optical Flow 用金字塔式平均與光流取代互相關；alignframes_lmbfgs 靠解析梯度的最佳化與軌跡平滑懲罰項快速求解。**不同演算法對同一 micrograph 給出的 whole-frame 位移不完全一致**，因為 whole-frame 軌跡本質上只是各局部區域運動的某種折衷，並不存在唯一正解（Fig.2B）——這正是需要 individual-particle motion correction 的理由。
- **alignparts_lmbfgs vs RELION polishing 對比**（p.122）：根本差異在於 alignparts_lmbfgs 只參照平均影像本身做對齊（影像分析流程「前段」），polishing 用（SNR 更高的）map projection 在角度已知後做對齊（流程「後段」）；polishing 的軌跡估計依賴精確的 Euler angle/CTF/顆粒均質性等假設；alignparts_lmbfgs 保留任意方向速度變化的軌跡，polishing 目前版本強制**線性等速**軌跡（不含旋轉）；polishing 的加權是資料驅動（不需輻射損傷實測值），但計算成本遠高於 alignparts_lmbfgs。
- **未來展望**（第 8 節，p.122）：目前所有方法只校正 xy 平面內平移，尚未處理樣品的**旋轉與 z 方向運動**（曝光下樣品的複雜彎曲變形）；此外可設計不同的資料收集方案以更好地處理曝光初期的劇烈運動。

## 可引用的公式/圖表
- **式 (1)（p.108）——Motioncorr 超定線性系統**：`A·s̄_mn = m̄_mn`（4 幀範例的矩陣形式），最小平方求 T−1 個真實幀間位移。
- **式 (2)（p.109）——Fourier 空間平移相位公式**：`φ_jt = k_x(j)·x_t·2π/N + k_y(j)·y_t·2π/N`，是所有基於 Fourier shift theorem 的 motion correction 演算法共同基礎。**03 章可直接引用作 shift theorem 的實例。**
- **式 (3)–(5)（p.109–110）——alignframes_lmbfgs 目標函數與解析梯度**。
- **式 (6)–(8)（p.110–111）——二階平滑懲罰項**（軌跡二階差分懲罰，`λ[(x_t−2x_{t−1}+x_{t−2})² + (γ_t−2γ_{t−1}+γ_{t−2})²]`）。
- **式 (9)–(10)（p.111–112）——高斯鄰近平滑**：`w_mn = exp(−d²_mn/2σ²)`，用於 alignparts_lmbfgs 與 RELION polishing 共用的局部軌跡平均。
- **式 (11)（p.112）——Unblur 正規化互相關**（單幀對其餘幀平均之互相關）。
- **式 (12)（p.114）——cubic spline 平滑泛函**：`Σ(x̄_t−f(t))² + λ∫f''(t)dt`。
- **式 (13)–(14b)（p.115）——Optical Flow 的一階泰勒展開光流方程**。
- **式 (16)（p.117）——RELION polishing 線性軌跡擬合目標函數**。
- **式 (17)–(18)（p.118–119）——RELION 資料驅動的幀加權**：相對 B-factor `B_t` 與 `C_t` 來自逐幀 FSC，加權 `w_t(1/d) = exp(B_t/4d²+C_t)/Σw_t'`。
- **Fig.1（p.104）**：whole-frame alignment 示意圖 + Thon rings 校正前後對比，展示高解析度資訊如何因對齊而恢復。**05/06 章講「為什麼需要 motion correction」的最佳配圖。**
- **Fig.2（p.113）**：alignparts_lmbfgs 個別顆粒軌跡圖、三種演算法對同一區域軌跡的比較、RELION polishing 的線性軌跡疊圖。**06 章講「不同演算法給出不同 whole-frame 解」的關鍵圖。**
- **Fig.3（p.120）**：RELION 逐幀 B_t/C_t 測量、最佳曝光加權曲線（不同解析度目標對應不同最佳總曝光量）、motion correction 前後的地圖品質改善實例（β-galactosidase、Complex-I、20S proteasome）。

## 適用章節
- **06（主要）**：整章都是 movie processing / motion correction 的方法學核心——Motioncorr 互相關法、alignframes/alignparts_lmbfgs 的最佳化框架、Unblur 的迭代對齊+CV 平滑、Optical Flow 的光流金字塔法、RELION particle polishing 的後驗軌跡估計，以及最終的方法比較表（第 7 節）。
- **05（次要）**：DDD 相機原理與 DQE 優勢（p.104）、輻射損傷與 critical exposure 曲線（Fig.3B）、beam-induced motion 作為 cryo-EM 特有的影像劣化來源。
- **03（次要）**：式 (2) 的 Fourier shift theorem 是所有互相關式對齊方法的數學基礎，可與傅立葉章節的 shift theorem 直接呼應；互相關求位移即為傅立葉域的相位相關（phase correlation）應用範例。
- **07（次要）**：合成資料若要模擬 movie（多幀累積、beam-induced motion、逐幀輻射損傷加權），本章式 (1)–(18) 的運動模型與加權方案可作為生成器設計依據。
