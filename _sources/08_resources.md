# 延伸資源與下一步

```{admonition} 本章重點
:class: important
- 這一章是「地圖」：告訴你每個主題可以去哪裡深入。
- 最完整的資源索引是老師維護的 [Computational-CryoEM](https://github.com/phonchi/Computational-CryoEM) 清單，本章按照本書章節對應整理其入口。
```

## 按主題的延伸資源

### 入門與綜述（對應第 5–6 章）

- **[Computational-CryoEM](https://github.com/phonchi/Computational-CryoEM)** — 依 workflow 12 階段整理的論文/軟體/教學索引，本書第 6 章的每個階段都能在這裡找到對應小節。
- Sigworth (2016), *Principles of cryo-EM single-particle image processing* — 單顆粒方法的物理與統計原理短篇綜述。
- Singer & Sigworth (2020), *Computational Methods for Single-Particle Electron Cryomicroscopy* — 數學視角的完整 review。
- Bendory, Bartesaghi & Singer (2020), IEEE Signal Processing Magazine review — 訊號處理視角。

### 主流軟體（對應第 6 章）

| 軟體 | 特色 |
|---|---|
| [RELION](https://relion.readthedocs.io/) | 貝葉斯最大後驗框架，學界標準 |
| [CryoSPARC](https://cryosparc.com/) | 商用、速度快、UI 友善 |
| [Scipion](http://scipion.i2pc.es/) | 整合多套軟體的 workflow 引擎 |
| [EMAN2](https://blake.bcm.edu/emanwiki/EMAN2) / [SPHIRE](http://sphire.mpg.de/) | 老牌全功能套件 |
| [ASPIRE-Python](https://github.com/ComputationalCryoEM/ASPIRE-Python) | 演算法研究導向，本書第 7 章使用 |

### 資料格式與慣例（對應第 7 章）

- MRC/MRCS 格式規格與 STAR 檔（RELION metadata）——第 7 章實作過讀寫。
- [3DEM Conventions](https://github.com/azazellochg/3DEM-conventions)：各軟體 Euler 角慣例對照——跨軟體轉檔時的頭號陷阱。
- 公開資料集：[EMPIAR](https://www.ebi.ac.uk/empiar/)（原始影像）、[EMDB](https://www.ebi.ac.uk/emdb/)（密度圖）。

### 影像處理教科書(對應第 1–4 章)

- Szeliski, *Computer Vision: Algorithms and Applications*（[免費電子版](https://szeliski.org/Book/)）— 第 3 章（濾波/傅立葉）、第 7 章（特徵）對應本書 1–4 章。
- [scikit-image 官方範例集](https://scikit-image.org/docs/stable/auto_examples/)。

## 建議的下一步練習

1. **異質性資料集**：`dataset` 內還有第二個構型 `70S_Conform2.mrc`。把兩個構型各生成一半影像混合，用第 7 章的 `purity_score` 檢驗 2D 分群能否分開兩種構型。
2. **打開被關掉的變因**：第 7 章刻意設 `offsets=0`（無平移）。打開隨機平移後，觀察對齊與分類的難度變化。
3. **CTF 修正實作**：結合第 3 章的 Wiener 去卷積與第 7 章的 CTF 參數，對模擬影像做 phase flipping，比較修正前後的 2D 平均品質。
4. **真實資料**：從 EMPIAR 下載一組小型資料集，用 RELION 或 CryoSPARC 跑完整 pipeline，對照第 6 章的每個階段。
