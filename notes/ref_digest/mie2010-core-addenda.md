# Methods in Enzymology 2010：SPA 核心章節補充 digest

本檔補足尚未有獨立長篇 digest、但直接支撐本站 SPA 主線的章節。頁碼採各章印刷頁；它是閱讀與引用路由，不取代回到 PDF 核對原文、公式與圖說。

## Resolution measures

來源：Penczek，〈Resolution Measures in Molecular Electron Microscopy〉，pp. 73–100。

- 解析度是「兩個獨立重建在某個空間頻率仍有多少一致訊號」的操作性判準。FSC 必須連同門檻、遮罩與資料切分方式一起報告。
- 單一全域數字會掩蓋方向性與局部差異；粒子取向不均、柔性或局部佔有率不同時尤其明顯。
- 本站用途：`06_workflow` 的 half-set／FSC／masking 限制；不能把 FSC crossing 寫成顯微鏡的硬體解析度。

## High-resolution SPA

來源：Cong 與 Ludtke，〈Single Particle Analysis at High Resolution〉，pp. 211–235。

- 章節把 specimen、取向覆蓋、影像擷取、particle boxing、CTF／envelope correction、reference-free 2D analysis、initial model、refinement 與結果評估串成完整 SPA 流程。
- 高解析度需要足夠且均勻的取向、可靠的 CTF／包絡模型，以及不受壞粒子或錯誤對齊主導的重建；這些條件彼此耦合，不能只靠增加粒子數補救。
- 本站用途：`05_background`、`06_workflow`；軟體操作細節有年代性，應引用方法限制而非把舊版 EMAN 步驟當成現行標準。

## Maximum-likelihood classification

來源：Scheres，〈Classification of Structural Heterogeneity by Maximum-Likelihood Methods〉，pp. 295–320。

- 2D／3D classification 將 class、orientation 與 shift 視為未知變數，以各種可能狀態的機率權重取代一次性的硬指派。
- 正規化、粒子擷取尺度、初始模型與類別數都會改變分類結果；分類不是自動揭示「真實構型數」的黑盒子。
- 本站用途：`06_workflow` 的 hard／soft assignment 與 heterogeneity；不把 class average 的視覺差異直接解讀成生物構型。

## Heterogeneous reconstruction

來源：〈Methods for Three-Dimensional Reconstruction of Heterogeneous Assemblies〉，pp. 321–341。

- 結構異質性會破壞「所有粒子是同一 3D 結構之投影」的基本假設，因此需要先分群、局部分析或多模型重建。
- 可靠性取決於訊號量、族群比例、取向分布及起始參考；少數族群在低 SNR 下可能不可辨識。
- 本站用途：`06_workflow` 的異質性章節；僅建立問題與限制，不把任一演算法列為普遍最佳解。

## Ewald sphere

來源：Leong 等人，〈Correcting for the Ewald Sphere in High-Resolution Single-Particle Reconstructions〉，pp. 369–380。

- 平面 central slice 是投影近似；解析度提高、粒子變厚時，Ewald sphere 曲率使不同深度的相位差不可再忽略。
- 該章介紹 paraboloid／prec correction，顯示曲率修正是高解析度延伸，而非入門 forward model 每次都必須開啟的項目。
- 本站用途：作為 Fourier slice theorem 的明確適用範圍；正文維持平面投影模型，將曲率列為進階限制。
