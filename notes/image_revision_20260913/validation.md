# 原 notebook 還原驗證

- 兩位讀者完整對照四章與原稿 171／134／120／38 cells；修正結果見 readers.md。
- 四章 clean kernel 執行成功，code cells 178／56／52／48，零例外；後續僅補閱讀錨點、標題層級、用詞及收合設定，Jupytext 保留輸出。
- pytest 全部 42 項通過；非互動 matplotlib 的 92 則警告不影響數值檢查。
- Jupyter Book 全站 warnings-as-errors 建置成功；HTML 連結、錨點、圖片 alt 檢查通過。
- Chromium 四章各 1440／390 px：圖片載入完整、頁面無橫向溢出、7 個互動／圖解框正常、pageerror 為 0。高低通 1／12／24 截止操作已驗證。
- 原始幾何 8 組 preset 對 skimage：矩陣誤差小於 1e-9，顯示差小於 2/255；FFT 與卷積互動對 NumPy／SciPy 誤差小於 1e-12。詳見 labs_acceptance.json。
- 99 個受保護檔案與四份原 notebook SHA256 全部一致；合成資料章與既有圖／資料沒有改動。
- 36 個原素材皆已本地保存，來源及雜湊見 assets.json；SPA 搬遷例子的參數見 spa_migration.json。
- 瀏覽器發現 _static 中的 iframe CSS／JS 會被全站自動載入，已移到 _support/image-lab；重新驗證零錯誤。
- 完整執行與建置 log、本機截圖保留在 logs/image_revision_20260913/（不發布）；可分享的數值與驗收摘要存於本目錄。
- FFT API 補充核對 SciPy https://docs.scipy.org/doc/scipy/reference/fft.html 與 NumPy https://numpy.org/doc/stable/reference/routines.fft.html。
