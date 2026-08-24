# 科學與教學複驗 ledger

來源：2026-08-24 Codex 複驗報告、本輪第 7 章審查與完成後的獨立科學複驗。`fixed` 均經公式／原始碼對照、針對性測試、notebook 執行或來源定位確認；沒有用「網站可 build」代替科學結案。

| ID | 嚴重度 | 問題 | 狀態 | 結案證據 |
|---|---|---|---|---|
| R00-01 | major | 原始資料誤寫成單張低 SNR 投影 | fixed | 00 明列 movie → micrograph → particle image → 3D map |
| R00-02 | major | 第 7 章不存在 | stale | 第 7 章與 TOC 已存在 |
| R00-03 | minor | 去雜訊與模板術語不一致 | fixed | glossary 與語言測試 |
| R01-01 | major | `imshow` 與笛卡兒座標原點混淆 | fixed | 01 座標段與 foundations test |
| R01-02 | major | float dtype 與數值範圍混淆 | fixed | float passthrough assertion |
| R01-03 | minor | Rec. 709 與 `rgb2gray` 係數未區分 | fixed | 版本實作係數 assertion |
| R01-04 | major | `warp` mapping 方向錯誤 | fixed | output-to-input mapping 示範 |
| R01-05 | major | micrograph 定義衝突 | fixed | 00、01、glossary 定義一致 |
| R02-01 | major | watershed 無法分開相鄰 instance | fixed | distance-transform markers；2-instance assertion |
| R02-02 | minor | HOG normalization 與參數不符 | fixed | 2×2 cells per block |
| R02-03 | major | 把 Gaussian／median 說成 picking 常規前處理 | fixed | 改為有條件教學工具並補 local normalization／whitening |
| R02-04 | minor | 中國用語 | fixed | 語言測試 |
| R03-01 | major | inverse DFT 索引與維度錯置 | fixed | 非方陣直接 DFT 對照 NumPy |
| R03-02 | major | DC 係數誤寫成平均 | fixed | sum 與 mean 各自 assertion |
| R03-03 | major | 正則化除法冒充精確 inverse | fixed | truncated inverse 與 Wiener 分開 |
| R03-04 | minor | 簡化 CTF 假設未列出 | fixed | CTF 主張範圍 admonition |
| R03-05 | major | CTF 零點誤稱由 defocus 唯一決定 | fixed | 公式明列 voltage、Cs、w、phase 等參數 |
| R03-06 | minor | 演算法、離焦、粒子術語混用 | fixed | glossary 與語言測試 |
| R04-01 | major | 兩訊號頻譜「幾乎相同」與數值不符 | fixed | 重設示範並量化 leakage |
| R04-02 | minor | 有限觀測頻率解析說得過度絕對 | fixed | window length／leakage 限制 |
| R04-03 | major | DWT 尺度與位移網格錯誤 | fixed | dyadic basis 公式 |
| R04-04 | major | refinement 與 decomposition level 混淆 | fixed | 兩者拆段示範 |
| R04-05 | major | 過度保證 cryo-EM 小波稀疏與去雜訊效果 | fixed | 降為有假設的視覺化／多尺度類比 |
| R04-06 | minor | scalogram、支撐區間與台灣用語 | fixed | 語言與拼字複核 |
| R05-01 | major | AlphaFold 訓練資料縮成晶體結構 | fixed | 移除不必要且過窄的訓練資料敘述 |
| R05-02 | minor | 「兩億種蛋白質」缺定義與日期 | fixed | 移除無法穩定支持的數字 |
| R05-03 | minor | 50 kDa 被寫成 NMR 硬上限 | fixed | 改寫為條件式方法比較 |
| R05-04 | major | 忽略 unscattered reference wave | fixed | 新成像章 weak-phase／interference 段 |
| R05-05 | minor | 20 e⁻/Å² 被寫成無條件典型劑量 | fixed | 標成較早期例子並區分現代總 fluence |
| R05-06 | major | SNR 未定義功率與估計範圍 | fixed | power-SNR、dB、SSNR 分開定義 |
| R05-07 | minor | 冷凍電顯、顆粒與翻譯腔 | fixed | glossary 與語言測試 |
| R06-01 | major | tomography 混入 SPA 主線 | fixed | workflow 明確排除 cryo-ET |
| R06-02 | major | DDD 中文名稱與縮寫不對應 | fixed | 統一 DED |
| R06-03 | major | phase flipping、Wiener 與多離焦混淆 | fixed | 03 與 06 各自分開說明 |
| R06-04 | major | 圓周平均後仍宣稱可估 astigmatism | fixed | 明列 2D fitting 與徑向近似界線 |
| R06-05 | major | CTF 振盪與 envelope attenuation 混淆 | fixed | 新成像章分開兩種效應 |
| R06-06 | minor | 2D class quality 被過度確定地解讀 | fixed | 列出多重原因與交叉診斷 |
| R06-07 | minor | class averaging 誤對應 wavelet thresholding | fixed | 移除錯誤直接對應 |
| R06-08 | major | refinement 被教成 hard assignment | fixed | 新重建章加入 soft likelihood weights |
| R06-09 | major | 第 7 章不存在 | stale | 第 7 章與 TOC 已存在 |
| R06-10 | minor | movie 量詞、離焦、閾值、Ab initio 術語 | fixed | glossary 與語言測試 |
| R08-01 | major | 練習依賴不存在的 Conform2／purity | fixed | 移除壞練習，改成自足觀察活動 |
| R08-02 | major | Wiener 與 phase flipping 寫成同一方法 | fixed | 改列兩種獨立觀察 |
| R08-03 | minor | 半形標點、review 與翻譯腔 | fixed | `speak-human-tw` 複核 |
| N07-01 | major | transfer function 與空間卷積核混淆 | fixed | 實／頻域 forward model 分式 |
| N07-02 | major | ASPIRE offset 分布說錯 | fixed | 明列預設為無界常態、σ=L/16；本章固定 0 |
| N07-03 | major | 缺檔與未定義 purity 練習 | fixed | 移除不可執行承諾 |
| N07-04 | major | 研究協定表沒有可追溯來源 | fixed | 改成可解析 citations 與 claim boundaries |
| N07-05 | major | 宣稱 ground truth 但未輸出 | fixed | E2E 產生 JSONL 與 run manifest |
| N07-06 | minor | 固定 seed 卻稱每次方向不同 | fixed | 改為可重現抽樣 |
| N07-07 | minor | 分布均勻與有限樣本均勻覆蓋混淆 | fixed | 文字與 sampler test 分開兩種主張 |
| N07-08 | major | global AWGN 邊界未說明 | fixed | manifest 記錄定義、realized SNR 與 limitations |
| N07-09 | major | 使用 ASPIRE 私有 `_metadata` | fixed | public `get_metadata(as_dict=True)` E2E 通過 |
| A01 | major | DQE 對 power-SNR 重複平方 | fixed | 改為 SNR_out／SNR_in，並說明 amplitude S/N 寫法 |
| A02 | major | ground truth 缺 rotation matrix 語意 | fixed | 每筆保存 matrix 與 image-plane-to-volume 方程 |
| A03 | major | claim map 指向 stale chapter／anchor | fixed | 刪除 stale claim，23 筆重綁顯式 MyST labels |
| A04 | minor | truncated inverse 的 `\tau` 字元損壞 | fixed | 公式原始碼與 build 複驗 |

## 獨立複驗摘要

- 44 個舊 finding：42 `fixed`、2 `stale`，0 open。
- ASPIRE 0.14.3 `RadialCTFFilter.evaluate` 對照本站 `ctf_1d`：500 個頻點 max abs diff `1.98e-14`，correlation `1.0`。
- 快速 E2E：STAR／MRCS shape、voxel size、metadata round-trip 一致；realized global SNR 約 `0.1`。
- Foundations 與 SPA 科學測試、全文 reference 驗證及零警告 Jupyter Book build 皆須在發布前再次執行。
