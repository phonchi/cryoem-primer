# 從影像處理到 Cryo-EM 單粒子分析

這本教材有兩條可以獨立閱讀、也能前後銜接的主線。第 1–4 章從數位影像的表示、濾波、分割、傅立葉轉換一路講到多尺度分析；第 5–7 章再把這些工具用於 Cryo-EM 單粒子分析（single-particle analysis, SPA）。

```{admonition} 讀完後，你應該能做什麼？
:class: important

- 正確讀寫、顯示、轉換與分析一般數位影像。
- 比較空間域濾波、分割、頻域處理與多尺度方法的假設和限制。
- 分清楚 movie frame、motion-corrected micrograph、particle image 與 3D density map。
- 說明投影、CTF 與雜訊如何構成簡化的前向模型。
- 產生一組取向與平移都已知的合成粒子影像，直接比較校正前後的圖形與誤差。
```

## 兩條閱讀路徑

| 路徑 | 適合誰 | 建議順序 |
|---|---|---|
| 一般影像處理 | 想建立完整數位影像處理基礎的讀者 | 01 → 02 → 03 → 04 |
| Cryo-EM 單粒子分析 | 已熟悉影像處理，想了解 SPA 的讀者 | 05 → 05A → 06 → 06A → 07 |

讀者可以只讀其中一條路徑，也可以從第 1 章一路讀到第 7 章。第 1–4 章本身是一套完整的影像處理教材，每章也會穿插 Cryo-EM 案例；第 5–7 章則會反覆連回前面的通用方法。

## 一套 SPA 資料怎麼走到 3D 密度圖？

直接電子偵測器先記錄一段 **dose-fractionated movie**。每個 frame 的電子劑量低，樣品也會在照射期間移動；將 frames 做 motion correction 並加總後，才得到後續處理常用的 **micrograph**。一張 micrograph 內含許多低訊雜比、經 CTF 調變的二維粒子投影，連同冰層、碳膜、污染物與偵測器雜訊。Particle picking 找出座標並切出 particle images；對齊、分類與 3D refinement 再從大量不同取向的粒子影像估計 3D 密度圖。

```text
Movie frames → motion correction → micrograph
                                      ├─ CTF estimation ─┐
                                      └─ particle picking ┴→ particle images
                                              → 2D classification
                                              → 3D reconstruction／refinement
                                              → half-map validation
```

Movie 與 micrograph 是不同層次的實驗觀測，particle image 是從 micrograph 擷取的分析單位，projection 則是前向模型裡的理想化物理量。這些名詞各自對應不同的資料來源與處理階段。

## 章節地圖

| 章節 | 核心問題 | 學完後可以做什麼 |
|---|---|---|
| 1　影像基礎 | 像素、色彩、曝光、取樣與幾何變換如何表示？ | 能正確處理 dtype、色彩、直方圖、幾何變換與常用檔案 |
| 2　濾波與分割 | 如何平滑、找邊緣、分割物件與描述局部特徵？ | 能依雜訊與任務選擇濾波、閾值、形態學與分割方法 |
| 3　傅立葉轉換 | 空間頻率如何連結濾波、卷積與影像復原？ | 能核對 DFT 慣例、設計頻域濾波器並理解正則化 |
| 4　多尺度分析 | 固定解析度不足時，如何同時觀察不同尺度？ | 能分清金字塔、STFT、CWT、DWT 與閾值假設 |
| 5–6　SPA 背景與流程 | 各計算步驟如何接回物理與統計模型？ | 能畫出完整流程並指出主要失敗模式 |
| 7　合成資料 | 粒子影像如何從 3D density 逐步產生？ | 能調整取向、平移、CTF 與雜訊，並比較每一層影像 |
| 8　延伸資源 | 下一步該查教科書、綜述還是軟體文件？ | 能依問題選對來源 |

## 先備知識與執行方式

讀者應熟悉 Python 基本語法、`NumPy` 陣列與高中程度的三角函數。複數、線性代數與機率若不熟，可以先讀主線；較深的推導放在可展開的進階區塊。

程式使用 `NumPy`、`SciPy`、`scikit-image`、OpenCV、`mrcfile`、`PyWavelets` 與 ASPIRE-Python。下載各章 notebook 後，可以逐格執行並改動參數；每次只改一個量，最容易看出它如何影響圖形與數值。

```{admonition} 讀圖時記得對照數值
:class: tip

教材先用小影像、簡化 CTF 與白雜訊，讓一項參數的作用能直接從圖上看出來。進入真實資料後，再把顯微鏡、樣品與偵測器特性逐一加回來。影像看起來更清楚時，先找出究竟是哪個像素值改變，再用訊雜比或頻域指標回答資訊是否增加；顯示效果和解析度回答的是不同問題。
```

## 理解檢查

1. 一般影像處理路徑中的哪些概念，也會出現在 Cryo-EM 流程？

   ```{dropdown} 參考答案
   取樣、座標與幾何變換會出現在 motion correction、粒子擷取與對齊；濾波與頻域分析用於估計 CTF、比較頻率內容與改善偵測條件；分割可用來建立污染或碳膜遮罩；多尺度表示則有助於比較不同大小的結構。

   這些連結可從資料流程看出來：movie frames 經過位移估計與重新取樣後得到 micrograph，micrograph 再經過座標偵測與裁切得到 particle images，後續的對齊與重建則大量使用傅立葉表示。

   但通用工具的名稱相同，並不保證參數也能照搬。例如濾波尺度要換算成實際長度，邊界處理也要配合粒子方框與後續分析。
   ```

2. 把單張 movie frame、motion-corrected micrograph 與 particle image 統稱為「原始投影」，會混淆哪些資料層次？

   ```{dropdown} 參考答案
   三者的產生方式和處理程度不同。Movie frame 是偵測器在一小段曝光時間內的記錄；motion-corrected micrograph 是將多張 frames 估計位移後對齊、加權並加總的結果；particle image 則是根據挑選座標，從 micrograph 裁下的局部方框。

   前向模型中的 projection 指三維密度在某一取向下的二維投影。實驗影像還包含 CTF、位移、雜訊、冰層背景與偵測器效應，可寫成簡化形式

   $$I_i = T_{\Delta_i}\,H_i P_{R_i}V + \varepsilon_i.$$

   把實驗記錄直接稱為原始投影，會漏掉它的處理歷程與成像效應。Particle image 已經歷位移校正、影格加總和裁切；偵測器直接輸出的是 movie frames。
   ```

3. 學習前向模型時，為什麼先使用已知真值的合成影像？

   ```{dropdown} 參考答案
   合成影像讓每個成像因素都能單獨調整，而且正確答案已知。例如，已知真實取向 $R_i$ 時，可比較估計取向 $\hat R_i$ 與 $R_i$；已知無雜訊影像 $x_i$ 時，可用 MSE

   $$\operatorname{MSE}=\frac{1}{N}\sum_{p=1}^{N}\bigl(\hat x_i[p]-x_i[p]\bigr)^2$$

   比較去雜訊結果。若改變一項已知參數後誤差突然升高，就能進一步追查取向、CTF、平移或雜訊是哪一環出了問題。

   解讀數值時，要先查看產生資料時放入了哪些條件。白高斯雜訊的練習可以回答這種雜訊下的誤差；真實 micrograph 還會出現複雜背景、樣品移動、污染、粒子異質性與成像模型誤差。接著換用真實資料，才能觀察方法遇到這些現象時的表現。
   ```
