# 慣例與術語速查

這一頁整理最容易讓程式「算得動，答案卻轉了九十度」的細節。跨軟體交換資料前，先確認這裡列出的慣例。陣列形狀相同時，軸向與座標語意仍可能不同。

## 本書固定用語

| 英文 | 本書用語 | 備註 |
|---|---|---|
| cryo-electron microscopy | 冷凍電子顯微鏡；簡稱冷凍電鏡 | 全書統一使用「冷凍電鏡」 |
| particle | 粒子 | 指裁切影像時仍用「粒子影像」 |
| defocus | 離焦 | 欠焦與過焦需另外說明正負號慣例 |
| denoising | 去雜訊 | 增強後影像進入定量重建前仍需驗證偏差 |
| visualization | 視覺化 | 影像增強與定量分析必須分開 |
| threshold | 閾值 | 程式參數可保留原英文名稱 |
| template matching | 模板匹配 | 必須交代正規化、搜尋角度與 template bias |
| review article | 綜述 | 第一次出現時可附英文 |

## 陣列與座標

- NumPy 影像以 `array[row, column]` 索引；`row` 對應垂直方向，`column` 對應水平方向。
- `imshow()` 預設把第 0 列畫在上方。笛卡兒座標的 \(y\) 軸方向與它不同，除非明確指定 `origin="lower"`。
- Fourier 陣列沿用 `[k_y, k_x]`；本書的公式也寫成 \(h[y,x]\)，讓陣列索引與平面座標維持一致。
- `fftshift()` 只重排頻率格點，Fourier transform 的正規化維持原值。

## 頻率、像素大小與解析度

- 像素大小以 Å/pixel 表示；空間頻率以 Å⁻¹ 表示。
- Nyquist 頻率是 \(1/(2\Delta x)\)。超過 Nyquist 的連續訊號會 alias 到可觀察頻帶內。
- 「解析度 \(d\)」常對應空間頻率 \(s=1/d\)，但 FSC 報告的解析度還依賴 half-set、mask 與閾值定義。
- CTF 的離焦、球差、振幅對比、額外 phase shift 與 Fourier 指數正負號，在不同軟體可能採不同慣例；比較參數時必須連同公式一起核對。

## Euler 角與旋轉

Euler 角有多套慣例。RELION、cryoSPARC、EMAN2 與程式庫可能使用不同的軸順序、主動／被動旋轉及座標手性。本書輸出已知真值時，同時記錄：

- 角度數值與單位；
- Euler convention 名稱；
- 對應的 \(3\times3\) rotation matrix；
- matrix 的作用方向與明確座標方程；
- image-plane-to-volume 或 volume-to-image 的方向。

可靠重現投影需要角度數值、單位、Euler 慣例、旋轉矩陣及其作用方向。

## Movie、micrograph 與 particle stack

1. **Movie frames** 是直接電子偵測器連續記錄的低劑量影格。
2. **Micrograph** 通常指經 gain correction、motion correction 與加權整合後的大幅影像。
3. **Particle stack** 是依 picking 座標從 micrograph 裁出的粒子影像集合。
4. **3D map／volume** 是由大量未知方向的粒子影像估計出的三維密度。

這四個資料層級具有不同的雜訊、中繼資料與校正方式；保留各自名稱可避免混淆。

## MRC、MRCS 與 STAR

- MRC header 的 `mode` 決定資料型別；`nx, ny, nz` 是檔案維度，生物學軸向仍需另查 header 與中繼資料。
- `.mrcs` 通常表示 2D 影像 stack；完整解讀仍需格式規範與中繼資料。
- STAR 儲存表格型中繼資料。影像位置常寫成 `index@stack.mrcs`，而 index 是否從 0 或 1 起算必須依格式規格確認。
- 讀回測試至少核對 shape、dtype、pixel size、影像順序、CTF 群組與抽樣像素值。

## 理解檢查

1. 一個 shape 為 `(128, 128)` 的陣列，還需要哪些資訊才能判斷第一軸是 \(x\) 或 \(y\)？

```{dropdown} 參考答案
Shape 記錄每個軸的長度；軸的物理意義要從索引慣例、header、中繼資料與讀取程式的文件判斷。第一軸可能是 row、$y$、section 或其他量。

常見誤解是把 `(128, 128)` 的第一個數字直接當成 $x$。NumPy 影像通常使用 `[row, column]`，也就是 `[y,x]`；顯示時還要確認 `origin` 設定。
```

2. 兩套軟體都輸出 `(rot, tilt, psi)`，交換前還要核對哪些慣例？

```{dropdown} 參考答案
欄位名稱相同，只能確認三個數字的標籤相同。可靠交換還要核對 Euler 軸順序、角度單位、主動或被動旋轉、座標手性，以及旋轉矩陣究竟把 image-plane 座標映到 volume，或採相反方向。

最安全的做法是用一個已知向量或已知投影，比較兩套軟體產生的 $3\times3$ rotation matrix 與實際作用結果。直接複製三個角度，常會得到旋轉方向、手性或平面內角度錯誤。
```

3. 同樣 shape 為 `(100, 128, 128)` 的 MRC 資料，為何可能代表 100 張粒子影像，也可能代表一個 3D 體積？讀取時應如何判斷？

```{dropdown} 參考答案
Shape 只記錄三個陣列軸的長度。對 particle stack 而言，第一軸通常列舉 100 張互相獨立的 2D 影像；對 3D map 而言，第一軸通常列舉同一體積的 100 個切片。兩者的陣列外形可以完全相同，物理意義卻不同。

判讀時要一起查看副檔名慣例、MRC header、voxel size、軸對應、中繼資料與產生檔案的軟體文件。`.mrcs` 常用於影像堆疊，但副檔名本身仍不足以確定內容。最可靠的檢查包括讀回數張索引影像、核對 STAR 的 `index@stack.mrcs` 對應，以及確認 3D 檢視器中的三個正交切面是否形成連續體積。
```
