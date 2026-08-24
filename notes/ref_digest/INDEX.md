# 文獻 Digest 索引

本目錄彙整 cryo-EM 教學網站（`cryoem-primer`）用到的參考文獻精讀筆記。每份 digest 對應 `References/` 下一篇論文或教科書章節，格式統一為「一句話定位／關鍵主張／可引用公式圖表／適用章節」。網站章節代號：**01** 陣列基礎、**02** 濾波與分割、**03** 傅立葉、**04** 小波、**05** cryo-EM 背景、**06** cryo-EM workflow、**07** 合成資料生成、**08** 資源。

> 覆蓋狀態與 SHA 的機器可讀權威檔是 `notes/reference_catalog/source_registry.jsonl`。本索引只提供人類閱讀入口；完整抽取文字保留在本機 `References/document_cache/`，不隨網站發布。

## 覆蓋摘要

- `full_digest`：SPA 核心來源，包含既有獨立 digest，以及 [2010 核心補充](mie2010-core-addenda.md)與 [2016 核心補充](mie2016-core-addenda.md)。
- `route_note`：相鄰方法與會變動的工具入口，見 [路由筆記](route-notes.md)及三本 computer vision 索引。
- `out_of_scope`：cryo-ET、helical reconstruction、electron crystallography、atomic model building，見 [範圍界線](scope-boundaries.md)。

| 檔名 | 一句話定位 | 適用章節 |
|---|---|---|
| [2015-principles-spa.md](2015-principles-spa.md) | 最精簡好讀的 SPR（single-particle reconstruction）原理導論，適合當 06 章 workflow 骨架，也給 07 章「為何要模擬雜訊與 CTF」的量化理由。 | 06（主）、07（主）、03、05、02 |
| [2020-ieee-review.md](2020-ieee-review.md) | 訊號處理視角的 cryo-EM 全景 review，附完整 12 階段 pipeline 流程圖與明確的影像形成/CTF 數學模型。 | 06（主）、07（主）、05（主）、03、02、01 |
| [2020-math-review.md](2020-math-review.md) | 全套文獻中數學最完整、公式編號最清楚的一篇：CTF 從繞射物理推導到工程模型、ML/MAP-EM 全式、FSC 與 B-factor 定量關係。 | 07（主）、06（主）、05（主）、03、02 |
| [cv-modern-approach-index.md](cv-modern-approach-index.md) | 以幾何與統計為骨幹的 CV 經典教科書，Part II「Early Vision」把 convolution、Fourier、sampling、edge、texture filter bank 講得完整且有數學細節。 | 02（主）、03、04、01、06、07 |
| [learning-opencv4-index.md](learning-opencv4-index.md) | OpenCV 4 + Python 3 實作型入門書，把「影像即 NumPy 陣列」的觀念一路帶到濾波、邊緣、輪廓、分割；01/02 章主力來源，03 僅概念性鋪陳，04/05/06/07 幾乎無對應。 | 01（主）、02（主）、03、04 |
| [mie2010-ch1-3d-reconstruction.md](mie2010-ch1-3d-reconstruction.md) | 3D 重建演算法的權威教科書式導論：central section theorem、ART/SIRT、FBP、gridding direct Fourier inversion 一路推到底。 | 06（主）、03（主）、07、02、05 |
| [mie2010-ch2-image-restoration.md](mie2010-ch2-image-restoration.md) | CTF 與影像復原的權威導論：從線性成像模型、CTF 相位函數，到 pseudoinverse／chi-square／Wiener filter 三大復原準則，以及各軟體 CTF 校正策略比較。 | 05（主）、06（主）、03、02、07 |
| [mie2010-ch10-maximum-likelihood.md](mie2010-ch10-maximum-likelihood.md) | 把 cryo-EM 的對位、分類與重建重新表述成統計估計問題：用 maximum likelihood 取代 cross-correlation／least-squares，以 EM 迭代求解。 | 06、05、07 |
| [mie2016-ch5-movie-processing.md](mie2016-ch5-movie-processing.md) | 以統一數學記號系統性比較五種 DDD movie motion correction 演算法（Motioncorr、alignframes/alignparts_lmbfgs、Unblur、Optical Flow、RELION particle polishing）。 | 06（主）、05、03、07 |
| [szeliski-cv-index.md](szeliski-cv-index.md) | 本站 01–04 章的主要教科書級參考：把卷積、邊界處理、非線性濾波、DFT/FFT、金字塔與 wavelet 放在同一套符號下講清楚；另含 Ch7.1–7.2 關鍵點/邊緣偵測重點。 | 01–04（主）、05、06、07 |

## 依章節反查

- **01 陣列基礎**：learning-opencv4-index（主）、szeliski-cv-index、2020-ieee-review、cv-modern-approach-index
- **02 濾波與分割**：cv-modern-approach-index（主）、learning-opencv4-index（主）、szeliski-cv-index（主）、mie2010-ch2-image-restoration、mie2010-ch1-3d-reconstruction、2020-ieee-review、2020-math-review、2015-principles-spa
- **03 傅立葉**：mie2010-ch1-3d-reconstruction（主）、szeliski-cv-index（主）、2020-math-review、2020-ieee-review、mie2016-ch5-movie-processing、mie2010-ch2-image-restoration、cv-modern-approach-index、learning-opencv4-index、2015-principles-spa
- **04 小波**：szeliski-cv-index（主）、cv-modern-approach-index（次要，僅 pyramid）
- **05 cryo-EM 背景**：mie2010-ch2-image-restoration（主）、2020-math-review（主）、2020-ieee-review（主）、mie2010-ch1-3d-reconstruction、mie2010-ch10-maximum-likelihood、mie2016-ch5-movie-processing、2015-principles-spa、szeliski-cv-index
- **06 cryo-EM workflow**：2015-principles-spa（主）、2020-ieee-review（主）、2020-math-review（主）、mie2010-ch1-3d-reconstruction（主）、mie2010-ch2-image-restoration（主）、mie2016-ch5-movie-processing（主）、mie2010-ch10-maximum-likelihood、cv-modern-approach-index、szeliski-cv-index
- **07 合成資料生成**：2015-principles-spa（主）、2020-ieee-review（主）、2020-math-review（主）、mie2010-ch10-maximum-likelihood、mie2010-ch1-3d-reconstruction、mie2010-ch2-image-restoration、mie2016-ch5-movie-processing、cv-modern-approach-index、szeliski-cv-index
- **08 資源**：目前無專屬 digest（本目錄十份皆為 01–07 內容型文獻）
