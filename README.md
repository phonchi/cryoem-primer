# Cryo-EM 單顆粒分析入門

這是一套給大三以上專題生使用的繁體中文教材。主線從數位影像與 Fourier 分析出發，接到 cryo-EM 的成像模型、單顆粒分析流程、3D 重建與驗證，最後用 ASPIRE 產生一套具有 ground truth 的合成資料。

- 公開網站：https://phonchi.github.io/cryoem-primer/
- 教材來源：`book/`
- 文獻摘要與 provenance：`notes/ref_digest/`、`notes/reference_catalog/`
- `.py` 是可執行章節的 single source of truth；`.ipynb` 由 Jupytext 同步產生。

## 建立環境

```bash
conda env create -f environment.yml
conda activate cryoem-book
python -m ipykernel install --user --name cryoem-book
```

## 常用命令

```bash
make sync          # 將 book/*.py 同步到 paired notebooks
make qa-fast       # 100 張模擬影像：測試、同步檢查與零警告 build
make build-full    # 5,000 張模擬影像：發布前完整 build
make linkcheck     # 外部連結檢查
make ingest        # 更新本機 References/document_cache（不會發布全文）
```

第 7 章的完整模擬會寫出約數百 MB 的 MRCS 與 provenance 檔。預設輸出位於 repo 外的暫存目錄；可用 `CRYOEM_OUTPUT_DIR` 指定位置。生成資料不會進入網站 artifact。

## 引用與版權

原始 PDF 與抽取全文只保留在本機 `References/document_cache/`。公開 repo 僅保存來源 registry、自行撰寫的摘要、claim map、短引文與 BibTeX。教材中的數值、公式、方法比較與限制性主張，均應能回溯到原始頁碼、式號或官方文件。

## 發布

推送 `main` 後，GitHub Actions 會建立完整網站並部署同一份已驗收的 artifact 至 GitHub Pages。不要手動執行 `ghp-import`，以免發布內容與 CI 驗證內容不同。
