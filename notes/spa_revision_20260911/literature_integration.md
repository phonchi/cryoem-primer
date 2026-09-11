# SPA 文獻整合交付（2026-09-11）

## 完成項目

- 本地 `References/document_cache` 更新為 68 PDF、68 成功擷取、66 SHA、2 組 binary duplicates；52→68 的新增及 Sigworth 檔名更正已納入。
- `source_registry.jsonl` 有 53 個 substantive PDF 記錄，含 2 份同作品 alternate PDF；它不是 53 篇獨立作品。13 份為卷冊前後資料、2 份為 binary aliases，沒有 unrouted 文件。
- 完整保留原 bibliography 位元內容作為前綴；只新增 11 keys。原 23 claims 不變，新增 13 manually verified claims，含 source hash／方法定位／適用界線。
- 深入筆記 `notes/ref_digest/spa-core-20260911.md` 包含主文與已提供的 ISAC、cryoSPARC 補充材料；各核心主張回溯至方法段，沒有把摘要擷取標成精讀。
- `scripts/reference_pipeline.py` 不再把未知 PDF 預設為 front matter；新增文件需要明確路由。手動未確認的 seed claim 預設 needs_review。
- 強制重建快取使用 generic skill 的 PDF extraction API，來源明確限定 `.pdf` 並排除 cache。這修正 generic CLI 會把既有 cache 當輸入的遞迴膨脹問題，沒有修改全域 skill。

## 新章來源對照

| SPA 章 | 優先來源 | 方法定位 |
|---|---|---|
| Background | cheng2015primer；既有 passmore2016 | Fig.1、pp.439–440 |
| Image formation | scheres2012bayes；既有 sigworth2016 | pp.408–410 Eqs.(1),(6) |
| Statistical inference | ma2019em、sigworth1998、scheres2005 | lecture pp.3–8；1998 pp.331–332；2005 Theory |
| Workflow | cheng2015primer；既有 rubinstein2016 | Fig.2、movie-processing 既有digest |
| Alignment/classification | sigworth1998、scheres2005、sorzano2010、yang2012、chen2014 | CL2D §2；ISAC Experimental Procedures/S2；gamma-SUP §§2.3,3.2 |
| Dimension reduction | chung2020 | §§1.3–1.5、2–2.1，pp.293–297 |
| Reconstruction | scheres2012bayes、punjani2017 | MAP Eqs.(8–12)；cryoSPARC Supplementary Note 1 |
| Heterogeneity | punjani2021；既有 scheres2010heterogeneity | 3DVA §5.2 Eqs.(4–8) |
| Resolution/validation | sorzano2022；既有 rosenthal2016 | §§4–5，FSC Eqs.(3–4) |

## 書目與來源界線

- `[2012] Relion.pdf` 是 JMB *A Bayesian View…*，key `scheres2012bayes`；不可當作同年的 JSB RELION implementation 論文。
- `[2] ML-EM in cryo-EM.pdf` 與 volume Ch.10 同作品但不同 SHA；`[2] ML-EM review (Section4).pdf` 是完整 Singer2020，並非 section-only PDF。對應關係另存 `alternate_editions.jsonl`。
- Ma & Ng 講義日期核對為 2019-05-13，官方公開頁：[Stanford 版本](https://cs229.stanford.edu/notes2020spring/cs229-notes8.pdf)。沒有捏造 DOI。
- 1998 DOI 由[出版社](https://www.sciencedirect.com/science/article/pii/S104784779894014X)核實。cryoSPARC 的卷頁由[Nature](https://www.nature.com/articles/nmeth.4169)核實；Cheng2015 卷頁由[Cell/Elsevier](https://www.sciencedirect.com/science/article/pii/S0092867415003700)核實。其他新增 DOI 與作者／刊名主要核對本地原論文首頁，不宣稱逐一完成外網 HTTP 檢查。
- γ-SUP 本地是 reprint，PDF頁不同於刊頁；請用 section／equation locator。3DVA 的 least-squares limit 不等於完整 PPCA posterior covariance 更新。
- Sorzano2022 supplementary experiments 不在資料夾中，本次沒有下載或驗證；只使用正文可支持的 bias/FSC 論述。不要把作者的 gold-standard 批判寫成取代 half-set 的普遍共識。
- 新 claim anchors 是概念名稱，主代理應在章節整合時對齊真實錨點；文獻管線只驗來源／hash／字段，不驗 HTML 錨點。

## 驗證

`literature_validation.md` 保留 force-ingest、validate、7項單元測試結果。確認原 bibliography prefix 和23項既有claim完整保留；沒有修改五個保護章、notebook或圖片。本次未建置book、未執行notebook、未發布。
