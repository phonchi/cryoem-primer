# 從影像處理到 Cryo-EM 單粒子分析

這是一套給大三以上專題生使用的繁體中文教材。第 1–4 章完整介紹數位影像處理，第 5–7 章再把這些工具用於 Cryo-EM 單粒子分析與合成資料。

- 公開網站：https://phonchi.github.io/cryoem-primer/
- 教材來源：`book/`
- 每個可執行章節皆提供 Jupyter notebook。

## 學生：建立環境與執行教材

```bash
conda env create -f environment.yml
conda activate cryoem-book
python -m ipykernel install --user --name cryoem-book
jupyter notebook
```

第一次執行第 7 章時，可先把影像數設為 100；確認流程後再生成完整的 5,000 張影像。輸出資料夾內會包含 STAR、MRCS、已知真值與執行參數紀錄。

## 引用與版權

教材使用自行撰寫的解釋與可執行範例，引用來源列在各章與書末。原始教科書、論文與資料集的權利仍屬原作者及發布單位。

## 維護者：測試與發布

```bash
make sync
make qa-fast
make build-full
make linkcheck
```

`.py` 是可執行章節的 Jupytext 來源，修改後需同步 `.ipynb`。推送 `main` 後，GitHub Actions 會建立並部署網站。
