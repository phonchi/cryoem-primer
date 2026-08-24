# 慣例與術語速查

這一頁處理最容易讓程式「算得動，答案卻轉了九十度」的細節。跨軟體交換資料前，先確認這裡列出的慣例；不要只看陣列形狀相同就假設語意相同。

## 本書固定用語

| 英文 | 本書用語 | 備註 |
|---|---|---|
| cryo-electron microscopy | 冷凍電子顯微鏡；簡稱冷凍電鏡 | 不使用「冷凍電顯」 |
| particle | 粒子 | 指裁切影像時仍用「粒子影像」 |
| defocus | 離焦 | 欠焦與過焦需另外說明正負號慣例 |
| denoising | 去雜訊 | 不把增強後影像直接當成無偏的重建輸入 |
| visualization | 視覺化 | 影像增強與定量分析必須分開 |
| threshold | 閾值 | 程式參數可保留原英文名稱 |
| template matching | 模板匹配 | 必須交代正規化、搜尋角度與 template bias |
| review article | 綜述 | 第一次出現時可附英文 |

## 陣列與座標

- NumPy 影像以 `array[row, column]` 索引；`row` 對應垂直方向，`column` 對應水平方向。
- `imshow()` 預設把第 0 列畫在上方。笛卡兒座標的 \(y\) 軸方向與它不同，除非明確指定 `origin="lower"`。
- Fourier 陣列沿用 `[k_y, k_x]`；本書的公式也寫成 \(h[y,x]\)，避免把陣列索引與平面座標混在一起。
- `fftshift()` 只重排頻率格點，不改變 Fourier transform 的正規化。

## 頻率、像素大小與解析度

- 像素大小以 Å/pixel 表示；空間頻率以 Å⁻¹ 表示。
- Nyquist 頻率是 \(1/(2\Delta x)\)。超過 Nyquist 的連續訊號會 alias 到可觀察頻帶內。
- 「解析度 \(d\)」常對應空間頻率 \(s=1/d\)，但 FSC 報告的解析度還依賴 half-set、mask 與閾值定義。
- CTF 的離焦、球差、振幅對比、額外 phase shift 與 Fourier 指數正負號，在不同軟體可能採不同慣例；只比較數值而不比較公式並不安全。

## Euler 角與旋轉

Euler 角不是單一標準。RELION、cryoSPARC、EMAN2 與程式庫可能使用不同的軸順序、主動／被動旋轉及座標手性。本書輸出已知真值時，同時記錄：

- 角度數值與單位；
- Euler convention 名稱；
- 對應的 \(3\times3\) rotation matrix；
- matrix 的作用方向與明確座標方程；
- image-plane-to-volume 或 volume-to-image 的方向。

只記錄三個角度、卻沒有記錄 convention，之後通常無法可靠重現投影。

## Movie、micrograph 與 particle stack

1. **Movie frames** 是直接電子偵測器連續記錄的低劑量影格。
2. **Micrograph** 通常指經 gain correction、motion correction 與加權整合後的大幅影像。
3. **Particle stack** 是依 picking 座標從 micrograph 裁出的粒子影像集合。
4. **3D map／volume** 是由大量未知方向的粒子影像估計出的三維密度。

這四個資料層級的雜訊、中繼資料與可做的校正不同，不應都簡稱為「原始影像」。

## MRC、MRCS 與 STAR

- MRC header 的 `mode` 決定資料型別；`nx, ny, nz` 是檔案維度，不自動保證生物學軸向。
- `.mrcs` 通常表示 2D 影像 stack，但副檔名本身不是完整 schema。
- STAR 儲存表格型中繼資料。影像位置常寫成 `index@stack.mrcs`，而 index 是否從 0 或 1 起算必須依格式規格確認。
- 讀回測試至少核對 shape、dtype、pixel size、影像順序、CTF 群組與抽樣像素值。

## 理解檢查

1. 一個 shape 為 `(128, 128)` 的陣列，能否單靠 shape 判斷第一軸是 \(x\) 還是 \(y\)？
2. 兩套軟體都輸出 `(rot, tilt, psi)`，為什麼仍不能直接交換？
3. 一張經過去雜訊的 micrograph 看起來更清楚，是否代表它一定適合送進 3D refinement？
