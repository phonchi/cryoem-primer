# [2010] An Introduction to Maximum-Likelihood Methods in Cryo-EM（Methods in Enzymology Vol. 482, Chapter Ten）
- 作者：Fred J. Sigworth, Peter C. Doerschuk, Jose-Maria Carazo, Sjors H. W. Scheres
- 出處：Methods in Enzymology, 2010, Vol. 482, pp. 263–294
- 頁碼對照：印刷頁 263–294 ＝ PDF 頁 1–32（`印刷頁 = PDF 頁 + 262`）
- 原檔：`References/2010/Chapter-Ten---An-Introduction-to-Maximum-Likelihood-Me_2010_Methods-in-Enzym.pdf`
- **選讀文獻：本 digest 為半頁摘要，非完整精讀。**

## 一句話定位
把 cryo-EM 的對位、分類與重建重新表述成統計估計問題的入門教材：用 maximum likelihood 取代傳統的 cross-correlation／least-squares，並以 expectation-maximization 迭代求解。

## 關鍵主張與內容
- 由 Bayes' rule 出發：後驗 ∝ likelihood × prior；若 prior 視為常數，maximum a posteriori 就退化成 maximum likelihood。likelihood 被看成模型 Θ 的函數，而非資料的函數（p.265）。
- 影像已對齊且雜訊為獨立同分布高斯時，MLE 就是單純的平均影像，等同 least-squares。作者強調這個特例常被誤用來否定 ML 的價值（p.272–273）。
- 未知的旋轉／平移屬於 hidden variable（nuisance parameter）。傳統做法把它們當成待估參數，ML 做法則是對它們積分（marginalization），只留下真正要估的 Θ（p.273–274）。
- prior 可以自然嵌入：例如把 particle picking 殘餘位移建模成高斯分布，等於對大位移做加權壓抑，比傳統「限定搜尋範圍」更有原則（p.274）。
- EM–ML 交替執行 E-step（算 hidden variable 的機率分布）與 M-step（更新模型），每次迭代保證 log-likelihood 不下降，收斂到局部極大（p.275–276）。
- 參考影像更新是「機率加權平均」而非只取單一最佳方位——這是 ML 與 cross-correlation 在實作上最直觀的差異（p.277）。
- 應用系譜：Sigworth (1998) 的 ML2D、kerdenSOM 分類、ML3D、含 scale/offset 的 MLn3D，以及 icosahedral virus、2D crystallography、tomography 的 subtomogram averaging（p.278–286，表 10.1 在 p.279）。
- 統計模型的假設（加成性、獨立、高斯白雜訊）是整套方法的成敗關鍵；真實 cryo-EM 影像因 CTF 與偵測器而不完全符合（p.286）。
- 主要瓶頸是計算量：single-particle 的積分橫跨 5–6 維空間。加速手段包括 domain reduction、grid interpolation、GPU、分塊式 partial E-step（p.289–290）。
- 展望明確指向 MAP 估計：加入 prior（例如懲罰粒子邊界外的溶劑密度，即 solvent flattening 的統計版本）可在資料量有限時改善模型品質（p.291）。

## 可引用的公式/圖表
- Eq. 10.1（p.265）Bayes' rule：`P(Θ|X) = P(X|Θ)P(Θ)/P(X)`——ML 與 MAP 的分野。
- Eq. 10.9（p.271）log-likelihood：`L(Θ) = Σ_i log P(X_i|Θ)`；Eq. 10.13（p.272）在對齊+高斯雜訊下 `A_MLE = (1/N)Σ_i X_i`。
- Eq. 10.14 / 10.16（p.273–274）marginal likelihood：`L(Θ) = Σ_i log ∫_φ P(X_i|φ,Θ)P(φ|Θ)dφ`——對 hidden variable φ（角度＋平移）積分。
- Eq. 10.19–10.21（p.276）EM 兩步：E-step 求 `Q(Θ,Θ⁽ⁿ⁾) = ∫ P(Y|X,Θ⁽ⁿ⁾) log P(X|Y,Θ)P(Y|Θ)dY`，M-step 取 `Θ⁽ⁿ⁺¹⁾ = argmax Q`；Eq. 10.24（p.277）機率加權重建 `A⁽ⁿ⁺¹⁾ = (1/N)Σ_i ∫ P(φ|X_i,Θ⁽ⁿ⁾) R_φ⁻¹ X_i dφ`。
- Fig. 10.4（p.275）機率隨相對方位的變化；Fig. 10.5（p.277）機率加權平均示意；Table 10.1（p.279）各 ML 方法的 data model、Θ 與 hidden variable 對照。

## 適用章節
- 06 cryo-EM workflow：2D 分類、3D refinement 階段的理論依據，可用 Table 10.1 說明各階段的模型參數與 hidden variable。
- 05 cryo-EM 背景：解釋為何低訊噪比下要用機率加權而非硬性指派方位。
- 07 合成資料生成：Eq. 10.10 / 10.15 的加成性高斯雜訊資料模型，正好對應投影＋雜訊的合成流程。
