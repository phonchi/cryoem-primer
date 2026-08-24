# 延伸資源與查閱路徑

先確定目前要解決的問題，再選對應的教科書、論文或軟體文件。

## 第一次讀這個領域

- Sigworth 的 *Principles of cryo-EM single-particle image processing* 適合在讀完成像模型後回頭看物理與統計假設 {cite}`sigworth2016`。
- Singer 與 Sigworth 的 *Computational Methods for Single-Particle Electron Cryomicroscopy* 適合補 Fourier、取向估計、最大概似與重建 {cite}`singer2020`。
- Bendory、Bartesaghi 與 Singer 的訊號處理綜述提供完整 SPA 問題地圖 {cite}`bendory2020`。
- [Computational-CryoEM](https://github.com/phonchi/Computational-CryoEM) 依處理階段整理方法、論文與軟體入口；方法細節可繼續閱讀其中連結的原始論文與官方文件。

## 依學習任務查找

| 目前要解決的問題 | 先讀本站 | 再查 |
|---|---|---|
| 影像陣列、取樣與 MRC | 第 1 章、慣例附錄 | MRC2014 規格、scikit-image 文件 |
| 濾波、模板匹配與 picking | 第 2 章 | Szeliski、APPLE／KLT picker 原始論文 |
| Fourier、CTF 與校正 | 第 3 章、成像模型 | Sigworth、CTFFIND 方法與官方文件 |
| SPA 流程 | 第 6 章 | RELION、cryoSPARC、Scipion 官方教學 |
| 重建與 FSC | 重建與驗證章 | gold-standard FSC 與 validation 原始文獻 |
| 合成資料與 ground truth | 第 7 章 | ASPIRE API、Simulation 實作與版本說明 |

## 軟體與官方文件

- [RELION 5.0 文件](https://relion.readthedocs.io/en/release-5.0/)：Bayesian refinement、classification 與 post-processing。
- [cryoSPARC](https://guide.cryosparc.com/)：互動式 SPA workflow 與 heterogeneous refinement。
- [Scipion 文件](https://scipion-em.github.io/docs/)：跨套件分析流程。
- [EMAN2](https://blake.bcm.edu/emanwiki/doku.php?id=eman2)：單粒子分析工具與教學。
- [ASPIRE-Python](https://github.com/ComputationalCryoEM/ASPIRE-Python)：第 7 章使用的數學與模擬程式庫。

版本會影響 API 與預設參數。重現本站結果時，以 `environment.yml` 鎖定的 ASPIRE 0.14.3 為準，不把最新版文件的行為直接套回舊版程式。

## 資料、格式與慣例

- [EMPIAR](https://www.ebi.ac.uk/empiar/) 典藏原始 movie、micrograph 與 particle data。
- [EMDB](https://www.ebi.ac.uk/emdb/) 典藏重建密度圖與相關中繼資料。
- [3DEM conventions](https://github.com/azazellochg/3DEM-conventions) 整理 Euler 角與跨軟體慣例；實作時仍要核對所用版本。
- 本站的 {doc}`appendix_conventions` 固定教材內使用的座標、頻率、術語與資料層級。

## 可以繼續做的觀察

這些活動不另附解答。每一項都先寫下預期，再改一個變因並記錄實際輸出。

1. 在第 7 章加入已知平移，比較未對齊平均與用已知 offset 校正後的平均。
2. 對同一批 CTF-modulated 影像分別做 phase flipping 與 Wiener-style correction，比較零點附近的行為。
3. 固定 clean projections，只改全域白高斯雜訊的訊雜比，量測實際訊雜比與類別平均影像的變化。
4. 選一個小型公開資料集，把其 movie、micrograph、particle stack 與 final map 對回第 6 章的資料層級。

## 本站引用文獻

```{bibliography}
:style: plain
```
