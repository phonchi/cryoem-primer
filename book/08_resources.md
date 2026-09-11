# 延伸閱讀與資料入口

各章的延伸閱讀依「先建立直覺，再讀核心方法」排列。本頁依 SPA 的問題整理原始論文，列出適合接續教材閱讀的段落。先選自己想理解的問題，再讀指定段落；讀完後回到章內公式，核對論文的觀測量、未知量與雜訊假設。

## 第一輪：讀懂整個問題

- **Cheng, Y., Grigorieff, N., Penczek, P. A. 與 Walz, T.（2015）．*A Primer to Single-Particle Cryo-Electron Microscopy*. Cell, 161(3), 438–449。** [DOI](https://doi.org/10.1016/j.cell.2015.03.050)。先看 Figure 1 與製樣段落，再讀影像處理流程，適合從生物問題進入計算方法 {cite}`cheng2015primer`。
- **Sigworth, F. J.（2016）．*Principles of Cryo-EM Single-Particle Image Processing*. Microscopy, 65(1), 57–67。** [DOI](https://doi.org/10.1093/jmicro/dfv370)。依投影、CTF、對齊、分類的順序閱讀，篇幅較短，可配合成像章中的圖形建立模型 {cite}`sigworth2016`。
- **Bendory, T., Bartesaghi, A. 與 Singer, A.（2020）．*Single-Particle Cryo-Electron Microscopy: Mathematical Theory, Computational Challenges, and Opportunities*. IEEE Signal Processing Magazine, 37(2), 58–76。** [開放全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7213211/)。從訊號處理角度讀前向模型、角度與位移估計，以及反問題；初讀可先略過較抽象的數學方法 {cite}`bendory2020`。
- **Singer, A. 與 Sigworth, F. J.（2020）．*Computational Methods for Single-Particle Electron Cryomicroscopy*. Annual Review of Biomedical Data Science, 3, 163–190。** [開放全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8412055/)。可作為整個 SPA 篇的共同參考，第 4 節尤其適合和 EM 章一起讀 {cite}`singer2020`。

## 第二輪：依統計問題深入

| 想理解的問題 | 本站章節 | 核心來源與讀法 |
|---|---|---|
| 為何低劑量需要統計平均？ | {doc}`05_background`、{doc}`05_image_formation` | Cheng primer；McMullan 等人的偵測器章；先比較訊號與雜訊如何隨劑量改變 |
| 為何對多組角度與位移加權？ | {doc}`05_statistical_inference` | Ma／Ng EM 講義；Sigworth 1998 的概似，以及對角度與位移的積分；先用一輪數值例子理解責任權重 |
| 對齊與分類為何互相影響？ | {doc}`06_alignment_classification` | ML2D、CL2D、ISAC、γ-SUP；分別比較目標函數、相似度、穩定性與離群值處理 |
| 低維表示保留了什麼？ | {doc}`06_dimension_reduction` | 2SDR 第 1.3 節及第 2 節；分清資料表示、秩選擇與最終重建任務 |
| CTF 很小時，為何需要收縮？ | {doc}`05_statistical_inference`、{doc}`06_reconstruction_validation` | Scheres 2012 Bayesian 論文；比較資料精度與先驗精度在分母的作用 |
| 快速重建如何減少計算？ | {doc}`06_reconstruction_validation` | cryoSPARC 2017 主文及 Supplementary Note 1；分開看統計模型與最佳化策略 |
| 二維影像如何描述三維變動？ | {doc}`06_heterogeneity` | 3DVA 第 5.2 節；先寫出模式依角度投影，再套用位移與 CTF 後的預測影像 |
| 兩張 half-map 一致能說明什麼？ | {doc}`06_resolution_validation` | Penczek 解析度章、Rosenthal 驗證章、Sorzano 2022；一起看獨立性、mask 與共同偏差 |

### EM、分類與穩定性

- **Ma, T. 與 Ng, A.（2019）．*The EM algorithm*, CS229 Part IX。** [Stanford 開放講義](https://cs229.stanford.edu/notes2020spring/cs229-notes8.pdf)。第 2 節由 Jensen 下界推到 E/M-step，適合在本站例子之後閱讀 {cite}`ma2019em`。
- **Sigworth, F. J.（1998）．*A Maximum-Likelihood Approach to Single-Particle Image Refinement*. Journal of Structural Biology, 122(3), 328–339。** [DOI](https://doi.org/10.1006/jsbi.1998.4014)。頁 329–332 的模型與更新式說明如何保留多組角度與位移，並依機率加權 {cite}`sigworth1998`。
- **Scheres, S. H. W. 等（2005）．*Maximum-likelihood Multi-reference Refinement for Electron Microscopy Images*. Journal of Molecular Biology, 348, 139–149。** [DOI](https://doi.org/10.1016/j.jmb.2005.02.031)。讀 Theory 與 Discussion，理解加入類別後如何共同估計參考影像、角度與位移 {cite}`scheres2005`。
- **Sorzano, C. O. S. 等（2010）．*A Clustering Approach to Multireference Alignment of Single-Particle Projections in Electron Microscopy*. Journal of Structural Biology, 171, 197–206。** [DOI](https://doi.org/10.1016/j.jsb.2010.03.011)。第 2 節比較 correntropy、分裂式分類及穩健指派，留意不同大小群集的平均影像雜訊 {cite}`sorzano2010`。
- **Yang, Z., Fang, J., Chittuluru, J., Asturias, F. J. 與 Penczek, P. A.（2012）．*Iterative Stable Alignment and Clustering of 2D Transmission Electron Microscope Images*. Structure, 20(2), 237–247。** [DOI](https://doi.org/10.1016/j.str.2011.12.007)。閱讀 Experimental Procedures 的 alignment stability 與 cluster reproducibility，再看補充材料的反例 {cite}`yang2012`。
- **Chen, T.-L. 等（2014）．*γ-SUP: A Clustering Algorithm for Cryo-Electron Microscopy Images of Asymmetric Particles*. The Annals of Applied Statistics, 8(1), 259–285。** [開放預印本](https://arxiv.org/abs/1205.2034)。第 2.3 節與第 3 節把穩健估計接到 self-updating；與 EM 比較時要核對各自的目標 {cite}`chen2014`。

### 降維、重建與異質性

- **Chung, S.-C. 等（2020）．*Two-Stage Dimension Reduction for Noisy High-Dimensional Images and Application to Cryogenic Electron Microscopy*. Annals of Mathematical Sciences and Applications, 5(2), 283–316。** [DOI](https://doi.org/10.4310/AMSA.2020.v5.n2.a4)。第 2 節及 Figure 1 串起 MPCA、SURE、PCA 與 GIC，之後再比較第 3 節實驗 {cite}`chung2020`。
- **Scheres, S. H. W.（2012）．*A Bayesian View on Cryo-EM Structure Determination*. Journal of Molecular Biology, 415, 406–418。** [DOI](https://doi.org/10.1016/j.jmb.2011.11.010)。Theory 的 Gaussian prior 與 MAP 更新解釋頻率相依收縮，Discussion 說明模型與估計的限制 {cite}`scheres2012bayes`。
- **Punjani, A., Rubinstein, J. L., Fleet, D. J. 與 Brubaker, M. A.（2017）．*cryoSPARC: Algorithms for Rapid Unsupervised Cryo-EM Structure Determination*. Nature Methods, 14(3), 290–296。** [DOI](https://doi.org/10.1038/nmeth.4169)。搭配 Supplementary Note 1 讀 marginal likelihood、SGD 與搜尋策略；文中的效能數值依當時資料與硬體解讀 {cite}`punjani2017`。
- **Punjani, A. 與 Fleet, D. J.（2021）．*3D Variability Analysis: Resolving Continuous Flexibility and Discrete Heterogeneity from Single Particle Cryo-EM*. Journal of Structural Biology, 213(2), 107702。** [DOI](https://doi.org/10.1016/j.jsb.2021.107702)。第 5.2 節是模型與更新的核心，第 4 節有助於理解線性表示的適用範圍 {cite}`punjani2021`。
- **Sorzano, C. O. S. 等（2022）．*On Bias, Variance, Overfitting, Gold Standard and Consensus in Single-Particle Analysis by Cryo-Electron Microscopy*. Acta Crystallographica D, 78, 410–423。** [DOI](https://doi.org/10.1107/S2059798322001978)。讀第 2 節誤差模型、第 5.1 節 FSC 與第 6 節討論；把作者對共同偏差的分析與 half-set 的用途一起理解 {cite}`sorzano2022`。

## 通用影像處理與軟體文件

一般影像處理四章的主要參考為 Szeliski 的 *Computer Vision: Algorithms and Applications*、Forsyth 與 Ponce 的 *Computer Vision: A Modern Approach*，以及 Howse 與 Minichino 的 *Learning OpenCV 4 Computer Vision with Python 3* {cite}`szeliski2022,forsyth2012,howse2020`。可先使用 [Szeliski 作者網站](https://szeliski.org/Book/)提供的書籍入口，再依章末指引閱讀卷積、頻域與多尺度章節。

需要操作軟體時，依正在使用的版本選官方教學：[RELION](https://relion.readthedocs.io/en/release-5.0/)、[cryoSPARC](https://guide.cryosparc.com/)、[Scipion](https://scipion-em.github.io/docs/)。本站 SPA 篇著重方法與統計解釋；實際操作步驟、輸入輸出與參數設定，請查軟體文件。

{doc}`07_synthetic_data`保留可執行的 ASPIRE 範例，產生帶有已知取向、平移、CTF 與雜訊的影像；其環境與參數以該章及 repository 的 README 為準。[ASPIRE-Python](https://github.com/ComputationalCryoEM/ASPIRE-Python)提供原始碼與文件入口。

## 公開資料與單位慣例

[EMPIAR](https://www.ebi.ac.uk/empiar/)提供 movie、micrograph 或粒子資料等實驗影像；[EMDB](https://www.ebi.ac.uk/emdb/)提供重建密度圖與相關資訊。挑選資料前先看資料層級、檔案大小、像素大小及原作者的處理說明，才能判斷它適合重跑哪一段流程。

跨軟體交換資料時，可將 [3DEM conventions](https://github.com/azazellochg/3DEM-conventions)與{doc}`appendix_conventions`並排閱讀。角度的軸順序、矩陣作用方向及 Å／pixel 的換算，會直接改變前向模型的預測。

## 本站引用文獻

以下集中列出教材引用的完整書目。各章末的閱讀說明提供較具體的閱讀順序與定位。

```{bibliography}
:style: plain
```
