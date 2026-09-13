# 從影像處理到 Cryo-EM 單粒子分析

這是一套給修過基礎統計、微積分與線性代數的大學生使用的繁體中文教材。一般影像處理四章、Cryo-EM 單粒子分析（SPA）篇與合成資料章分別編排。SPA 篇從成像與統計推論走到分類、降維、重建、異質性與驗證，正文以動機、圖解、核心公式及完整例子教學，詳細推導可展開閱讀。

- 公開網站：https://phonchi.github.io/cryoem-primer/
- 教材來源：`book/`
- 每個可執行章節皆提供 Jupyter notebook。

一般影像處理四章依原始 notebook 的教學順序、圖片與參數編排：像素與幾何 → 局部濾波與分割 → Fourier、金字塔與特徵 → 小波與去雜訊。原有數值例子和靜態圖保留；已確認的科學敘述及淘汰 API 另作修正。逐 cell 對照保存在 `notes/image_revision_20260913/ch01_coverage.md` 至 `ch04_coverage.md`。

幾何變換、一維／二維卷積、Fourier 基底與濾波掃描，可直接在瀏覽器操作。互動以本機 JavaScript／Canvas 執行，不需要 Python 伺服器、API key 或付費服務；Python 靜態圖提供相同實驗的數值對照。第 7 章仍是獨立的 SPA 合成資料練習。


## 學生：建立環境與執行教材

```bash
conda env create -f environment.yml
conda activate cryoem-book
python -m ipykernel install --user --name cryoem-book
jupyter notebook
```

圖片存於 `book/images/original_notebooks/`，共用顯示與互動程式位於 `book/_support/` 及 `book/_support/image-lab/`。執行下載的 notebook 時，請保留這些相對目錄。

第一次執行第 7 章時，可先把影像數設為 100；確認流程後再產生完整的 5,000 張影像。輸出資料夾內會包含 STAR、MRCS 與每張影像的已知真值。

## 引用與版權

教材使用自行撰寫的解釋與可執行範例，引用來源列在各章與書末。原始教科書、論文與資料集的權利仍屬原作者及發布單位。

## 維護者：測試與發布

```bash
make sync
make qa-fast
make build-full
make linkcheck
```

`.py` 是可執行章節的 Jupytext 來源，修改後需同步相應 `.ipynb`。只改特定章節時，使用 `jupytext --sync book/章節名稱.py`，避免重新同步未修改章節。

原始圖片與示意圖的來源、尺寸及雜湊記於 `notes/image_revision_20260913/assets.json`。`scripts/fetch_original_assets.py` 可依專案父目錄的四份原 notebook 重新擷取素材；原 notebook 是來源，不應被教材同步流程覆寫。引用路由由 `scripts/reference_pipeline.py validate` 核對。推送 `main` 後，GitHub Actions 會建立並部署網站。
