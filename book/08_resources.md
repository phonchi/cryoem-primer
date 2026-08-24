# 延伸閱讀與實作入口

遇到陌生術語時，不必從頭讀完一整本書。先確認問題落在分析流程、成像物理、程式實作或公開資料，再從下面的入口往下找。

## 先看懂整體流程

想先建立 SPA 的全貌，可從第 6 章的 movie、micrograph、particle stack 與 3D map 四個資料層級開始，再閱讀 Bendory、Bartesaghi 與 Singer 的訊號處理綜述 {cite}`bendory2020`。若想沿著每個處理步驟繼續找方法與軟體，[Computational-CryoEM](https://github.com/phonchi/Computational-CryoEM) 整理了 motion correction、CTF estimation、particle picking、classification 與 reconstruction 的入口。

第一次操作完整資料時，可選一套官方教學跟著做：

- [RELION 5.0 文件](https://relion.readthedocs.io/en/release-5.0/)：Bayesian refinement、classification 與 post-processing。
- [cryoSPARC guide](https://guide.cryosparc.com/)：互動式 SPA 分析流程與異質性 refinement。
- [Scipion 文件](https://scipion-em.github.io/docs/)：在同一個流程中串接多套軟體。
- [EMAN2](https://blake.bcm.edu/emanwiki/doku.php?id=eman2)：單粒子分析工具與教學。

## 補成像物理與數學

第 3 章介紹 Fourier transform，第 5 章把投影、平移、CTF、偵測器與雜訊連成前向模型。讀完這兩章後，可用下列資料補上推導與方法背景：

- Sigworth 的 *Principles of cryo-EM single-particle image processing* 深入說明成像物理與統計假設 {cite}`sigworth2016`。
- Singer 與 Sigworth 的 *Computational Methods for Single-Particle Electron Cryomicroscopy* 涵蓋 Fourier、取向估計、最大概似與重建 {cite}`singer2020`。
- 第 6 章的 half-map 與 FSC 可再對照 [RELION post-processing 文件](https://relion.readthedocs.io/en/release-5.0/)，實際查看 mask 對曲線的影響。

查 CTF 公式時，請連同離焦正負號、Fourier 指數、振幅對比與額外 phase shift 的定義一起核對。Euler 角跨軟體交換則可搭配 [3DEM conventions](https://github.com/azazellochg/3DEM-conventions) 與本站的 {doc}`appendix_conventions`。

## 動手實作

前四章使用 NumPy、SciPy、scikit-image 與 PyWavelets 練習一般影像處理；第 7 章使用 [ASPIRE-Python](https://github.com/ComputationalCryoEM/ASPIRE-Python) 產生已知取向、平移、CTF 與雜訊的粒子影像。

執行第 7 章前，先在 Python 檢查版本：

```python
import aspire

print(aspire.__version__)
```

本章的程式以 ASPIRE 0.14.3 測試。若顯示其他版本，請改查相同版本的 API 文件，並留意 `Simulation`、`RadialCTFFilter` 與旋轉慣例是否相同。

可依序做以下觀察，每次只改一個變因：

1. 把第 7 章的 `max_shift_px` 設成 0、2、4、8 pixels，比較未校正平均與已知平移校正後的平均。
2. 對同一批經 CTF 調變的影像分別做 phase flipping 與 Wiener-style correction，觀察零點附近的差異。
3. 固定無雜訊投影，只改全域白高斯雜訊的目標訊雜比，量測實際訊雜比與平均影像的變化。
4. 從同一張 micrograph 取一小塊背景與一組粒子框，比較兩者的 power spectrum，找出冰、碳膜或漂移留下的方向性特徵。

每次先寫下預期，再執行程式並比較輸出。若結果和預期不同，先檢查單位、軸順序與亂數種子，再回到公式找原因。

## 找公開資料

- [EMPIAR](https://www.ebi.ac.uk/empiar/) 典藏原始 movie、micrograph 與 particle data，適合練習從早期處理步驟開始的分析。
- [EMDB](https://www.ebi.ac.uk/emdb/) 典藏重建密度圖與相關中繼資料，適合練習讀取 MRC volume、像素大小與對稱性資訊。

下載資料前先查看檔案大小、像素大小、顯微鏡參數與作者提供的處理說明。開始分析後，可把每個檔案對回 movie、micrograph、particle stack 或 3D map，避免把不同層級的資料當成同一種影像。

## 依問題快速查找

| 想解決的問題 | 先讀本站 | 接著查 |
|---|---|---|
| 影像陣列、取樣與 MRC | 第 1 章、慣例附錄 | MRC2014 規格、scikit-image 文件 |
| 濾波、模板匹配與 picking | 第 2 章 | Szeliski、APPLE／KLT picker 原始論文 |
| Fourier、CTF 與校正 | 第 3 章、成像模型 | Sigworth、CTFFIND 方法與官方文件 |
| SPA 分析流程 | 第 6 章 | RELION、cryoSPARC、Scipion 官方教學 |
| 重建、half-map 與 FSC | 重建與驗證章 | Singer 與 Sigworth、RELION post-processing 文件 |
| 合成粒子影像 | 第 7 章 | ASPIRE API 與 `Simulation` 文件 |

## 本站引用文獻

```{bibliography}
:style: plain
```
