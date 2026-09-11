# 全文讀者審查與銜接修正

使用者要求最後由多位讀者全文閱讀，再修正並推送。兩位原生 agent 各自按 `_toc.yml` 完整讀完 17 個來源頁面，沒有只看目錄或抽樣，也沒有互相交換結論。Reader A 以首次學 SPA 的統計本科生視角閱讀；Reader B 以已有一般影像處理背景的本科生視角閱讀。

## 實際閱讀範圍

`00_intro.md`、`01_image_basics.py`、`02_filter_segment.py`、`03_fourier.py`、`04_wavelet.py`、`05_background.md`、`05_image_formation.md`、`05_statistical_inference.md`、`06_workflow.md`、`06_alignment_classification.md`、`06_dimension_reduction.md`、`06_reconstruction_validation.md`、`06_heterogeneity.md`、`06_resolution_validation.md`、`07_synthetic_data.py`、`08_resources.md`、`appendix_conventions.md`。

包含 `.py` 的全部教學段落、可見程式、圖說及收合答案，不重讀內容相同的 `.ipynb`。這是全文連貫性審查；程式執行、公式與圖片渲染另外驗證。

## 讀者 findings 與實際修正

| 問題 | 讀者 | 修正位置與內容 |
|---|---|---|
| Movie 初次處理卻突然需要三維參考，似乎循環依賴 | A、B | `06_workflow` Particle 層級明示：取得初步重建後回到 frames 精修；初次可採 frame／patch 校正，再進入分類與重建 |
| Guinier plot 尚未交代兩軸與斜率 | A | `06_resolution_validation` 說明橫軸 s²、縱軸 log 振幅、衰減斜率 −B_damp/4 與銳化 B=−B_damp，並區分對數功率的2倍斜率 |
| 同一θ在相鄰段落分別代表平面旋轉與Euler傾角 | B | `06_alignment_classification` 把平面旋轉及其誤差統一記為ψ，與圖中平面內旋轉ψ一致；模型參數Θ保留 |
| 一般Fourier章與SPA的CTF正負號似乎矛盾 | A、B | `05_image_formation` 明示無額外phase shift時兩章χ互為負值，sin奇／cos偶，使相同設定產生相同CTF |
| SURE式中的div尚未定義 | A、B | `06_dimension_reduction` 補散度偏導加總，並以固定投影div f=tr P=r接回前例；保留資料依存基底需另算導數的提醒 |
| 合成章多角度影像的置中平均容易被誤當class average | A、B | `06_resolution_validation` 的合成章導讀說明這是同一批角度／CTF下隔離位移效應的對照，尚未估旋轉、分類或重建 |

六項均在SPA可修改範圍內完成；前四章及合成章沒有修改。

## 確認無誤的主線

兩位讀者都認為無需重排章序。一般影像的像素／取樣→局部操作→Fourier→多尺度可連續閱讀；首頁將一般影像處理與SPA清楚分篇。SPA的平均→成像→Gaussian概似→Bayes／EM→MAP→對齊分類→降維→重建→異質性→FSC與共同偏差逐步增加問題與假設。重複公式是在新問題中應用，有教學用途；SNR定義與CTF-filtered target的差異已有交代。

Reader A 曾暫記02→03的whitening承諾，全文對照後確認03的頻譜濾波可合理承接，因此撤回，不以此修改保護章。
