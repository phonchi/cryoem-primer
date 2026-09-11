# 慣例與術語速查

這一頁整理影像處理常用的座標、單位與檔案慣例。跨軟體交換資料前，先確認這裡列出的慣例。陣列形狀相同時，軸向與座標語意仍可能不同。

## 本站固定用語

| 英文 | 本站用語 | 備註 |
|---|---|---|
| cryo-electron microscopy | 冷凍電子顯微鏡；簡稱冷凍電鏡 | 全書統一使用「冷凍電鏡」 |
| particle | 粒子 | 指裁切影像時仍用「粒子影像」 |
| defocus | 離焦 | 欠焦與過焦需另外說明正負號慣例 |
| denoising | 去雜訊 | 定量重建前要比較去雜訊影像與原始資料 |
| visualization | 視覺化 | 觀看增強影像時，另保留未增強版本供定量分析 |
| threshold | 閾值 | 程式參數可保留原英文名稱 |
| template matching | 模板匹配 | 使用時記下正規化方式、搜尋角度與 template bias |
| review article | 綜述 | 第一次出現時可附英文 |

## 陣列與座標

- NumPy 影像以 `array[row, column]` 索引；`row` 對應垂直方向，`column` 對應水平方向。
- `imshow()` 預設把第 0 列畫在上方。笛卡兒座標的 \(y\) 軸方向與它不同，除非明確指定 `origin="lower"`。
- Fourier 陣列沿用 `[k_y, k_x]`；本站的公式也寫成 \(h[y,x]\)，讓陣列索引與平面座標維持一致。
- `fftshift()` 只重排頻率格點，Fourier transform 的正規化維持原值。

## 頻率、像素大小與解析度

- 像素大小以 Å/pixel 表示；空間頻率以 Å⁻¹ 表示。
- Nyquist 頻率是 \(1/(2\Delta x)\)。超過 Nyquist 的連續訊號會 alias 到可觀察頻帶內。
- 「解析度 \(d\)」常對應空間頻率 \(s=1/d\)，但 FSC 報告的解析度還依賴 half-set、mask 與閾值定義。
- CTF 的離焦、球差、振幅對比、額外 phase shift 與 Fourier 指數正負號，在不同軟體可能採不同慣例；比較參數時必須連同公式一起核對。

## Euler 角與旋轉

Euler 角有多套慣例。RELION、cryoSPARC、EMAN2 與程式庫可能使用不同的軸順序、主動／被動旋轉及座標手性。交換 Euler 角時，應一起記錄：

- 角度數值與單位；
- Euler convention 名稱；
- 對應的 \(3\times3\) rotation matrix；
- matrix 的作用方向與明確座標方程；
- image-plane-to-volume 或 volume-to-image 的方向。

合成資料章的 `ground_truth.jsonl` 提供一個具體例子：`rot, tilt, psi` 使用 ASPIRE 的 ZYZ Euler angles，單位是 radians；`rotation_matrix` 則記錄同一個旋轉，並採用

$$
\mathbf q_{\mathrm{volume}}=R[k_x,k_y,0]^{\mathsf T}
$$

的 image-plane-to-volume 方向。把這組資料交給另一套軟體前，可先用一個已知方向產生投影，確認觀察方向、平面內旋轉與手性都相同，再交換整批角度。

## Movie、micrograph 與 particle stack

1. **Movie frames** 是直接電子偵測器連續記錄的低劑量影格。
2. **Micrograph** 通常指經 gain correction、motion correction 與加權整合後的大幅影像。
3. **Particle stack** 是依 picking 座標從 micrograph 裁出的粒子影像集合。
4. **3D map／volume** 是由大量未知方向的粒子影像估計出的三維密度。

這四個資料層級的雜訊、中繼資料與校正方式不同，討論時要分清楚。

## MRC、MRCS 與 STAR

- MRC header 的 `mode` 決定資料型別；`nx, ny, nz` 是檔案維度，生物學軸向仍需另查 header 與中繼資料。
- `.mrcs` 通常表示 2D 影像 stack；完整解讀仍需格式規範與中繼資料。
- STAR 儲存表格型中繼資料。影像位置常寫成 `index@stack.mrcs`；index 從 0 或 1 起算，請依格式規格確認。
- 讀取後先核對資料用途，再查看 shape、dtype、pixel size、影像順序與 CTF 群組。顯示第一張、中間一張與最後一張影像，也能快速發現索引偏移或軸順序錯誤。

## SPA 公式中的統一記號

| 記號 | 本篇意義 |
|---|---|
| $V$、$V_k$ | 三維密度與第 $k$ 個構形的密度 |
| $y_i$、$X_i$ | 第 $i$ 張影像的向量或矩陣表示 |
| $R_i$ | image-plane-to-volume 的三維旋轉 |
| $\mathbf t_i$ | 平面內位移，須配合頻率的單位 |
| $P_{R_i}$ | 在 $R_i$ 方向的投影算子 |
| $H_i$ | Fourier 空間的 CTF；實空間對應 PSF |
| $z_i$ | 統計推論章中未知的角度、位移與可能的類別；異質性章會另外定義連續座標 |
| $r_i(z)$ | 給定影像與目前模型後的正規化候選權重 |
| $\sigma^2$、$\Sigma$ | 雜訊變異數、雜訊共變異數矩陣 |

在 SPA 篇的旋轉定義下，投影與 Fourier slice 一起寫成

$$
(P_RV)(x,y)=\int V(R[x,y,z]^{\mathsf T})\,dz,
\qquad
\widehat{P_RV}(k_x,k_y)=\widehat V(R[k_x,k_y,0]^{\mathsf T}).
$$

若某個來源把旋轉定義為 volume-to-image，該來源的矩陣 $R_{\mathrm{v2i}}$ 與本篇關係為 $R_{\mathrm{v2i}}=R^{\mathsf T}$。先比較座標方程，再比較矩陣數值，就能看出公式為何在不同位置使用轉置。

## 像素、頻率與統計量的換算例子

像素大小為 2.82 Å/pixel 時，Nyquist 頻率為 $1/(2\times2.82)\approx0.1773$ Å⁻¹，對應的理想取樣長度是 5.64 Å。實驗可達解析度還受劑量、CTF、取向與估計誤差影響。

以 Å⁻¹ 表示的頻率 $\mathbf k$，要搭配以 Å 表示的位移 $\mathbf t$，使平移因子 $e^{-2\pi\mathrm{i}\mathbf k\cdot\mathbf t}$ 的指數無單位。若程式使用 cycles/pixel，位移就應使用 pixels。

本篇 SNR 採訊號功率除以雜訊功率，SSNR 用各頻率的功率比；採振幅比的來源需要平方後才能比較。共變異數估計、遮罩、扣平均與頻帶範圍也要一併說明。全影像的單一 SNR 和局部頻帶 SNR 分別描述不同尺度的資訊。

## 檔案之間如何對照？

MRCS 存放影像堆疊，STAR 記錄影像索引、CTF 與其他參數。只有 STAR 而沒有對應影像，無法重建粒子像素；只有 MRCS 而缺少必要的 CTF 與像素資訊，則容易用錯成像模型。合成資料另提供已知真值，讀取時應逐筆以影像索引對應，避免將真值順序與重新排序的粒子堆疊混用。

## 延伸閱讀

- **Singer, A. 與 Sigworth, F. J.（2020）．*Computational Methods for Single-Particle Electron Cryomicroscopy*. Annual Review of Biomedical Data Science, 3, 163–190。** [開放全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8412055/)；第 2 節的成像模型提供投影、CTF 與座標的脈絡 {cite}`singer2020`。
- **3DEM conventions：** [旋轉慣例對照](https://github.com/azazellochg/3DEM-conventions)。交換 Euler angles 時對照軸順序、主動／被動旋轉與矩陣作用方向，並使用一個已知投影確認。
- **MRC2014：** [格式規格](https://www.ccpem.ac.uk/mrc_format/mrc2014.php)。讀 header 的軸向、取樣與 voxel size 欄位，搭配合成資料章的讀取例子。
