# 原稿索引移出教學內容

移除四章共 486 個 source-cells 標記（包括 HTML 註解）。原稿對照改存 source_cell_map.json；舊 coverage 紀錄中的標記指向此表及 baseline commit，無須出現在教學來源。source-image-url 改成「圖片來源」，保留出處。

直接同步修改 py／ipynb 來源，逐格比對一致；所有執行輸出保持不變，包含圖片及互動。

驗證：9 項來源測試通過；全站 warnings-as-errors 建置與 HTML 檢查通過；Chromium 逐章確認 HTML 及可見正文皆無索引標記。99 個受保護檔案及四份原 notebook 雜湊不變。
