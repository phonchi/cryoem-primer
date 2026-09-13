# 第 1 章必要檢查

使用本地 `/home/phonchi/miniconda3/envs/cryoem-book/bin/python`。

- Python AST、所有 Jupytext code cell compile：通過。
- Jupytext 解析為 278 cells，原 source-cell ID 0–170 全覆蓋；markdown 中沒有 source-cells 工程標記洩漏。
- Euclidean(rotation=π/12, translation=(30,-20)) 與原手寫矩陣逐項相等：通過。
- 原 text 四點以 from_estimate 建立映射；最大點對殘差 1.7053e-13；直接 warp 輸出 (50,300)：通過。
- JPEG 直接 imageio.v3 記憶體編碼/讀回 (10,10,3)：通過；使用小陣列只測codec API，不作替代教材圖片。
- 一般教材無SPA/MRC/理解檢查；保留的kernel metadata cryoem-book 是環境名稱。
- 未執行全章，待root完成 `_support` 與原圖後統一執行；未同步ipynb，未建置。

## 整章執行發現並修正
原 venusaur-f.png 為80×80 palette PNG，透明索引0；現有skimage reader回傳RGB，無法執行原四通道像素操作。改以Pillow明確解碼RGBA，保留原圖透明資訊及原cell67–70操作，不修改原圖。
