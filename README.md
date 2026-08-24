# Cryo-EM 單顆粒分析入門（教學網站）

專題生教材：從一般影像處理到 cryo-EM 單顆粒分析，重點是**合成資料的生成**（ASPIRE、70S 核糖體、5,000 張模擬粒子）。

- 網站：https://phonchi.github.io/cryoem-primer/
- 教材本體：`book/`（Jupyter Book，8 章；各章以 jupytext py:percent 格式撰寫，`.py` 是 single source of truth）
- 文獻摘要：`notes/ref_digest/`（供備課引用）

## 建置

```bash
conda create -n cryoem-book --clone aspire-2sdr   # 或依 environment.yml
conda run -n cryoem-book pip install jupyter-book jupytext ipykernel \
  opencv-python-headless PyWavelets ghp-import
conda run -n cryoem-book python -m ipykernel install --user --name cryoem-book

# .py 修改後同步回 .ipynb
conda run -n cryoem-book jupytext --sync book/*.py

conda run -n cryoem-book jupyter-book build book/
# 開啟 book/_build/html/index.html
```

首次建置會實際執行所有 notebook（第 7 章生成 5,000 張 130×130 模擬影像，需數分鐘；輸出寫到 `book/data/output/`，不進版控）。

## 部署

```bash
conda run -n cryoem-book ghp-import -n -p -f book/_build/html
```
