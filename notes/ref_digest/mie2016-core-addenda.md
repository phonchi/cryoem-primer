# Methods in Enzymology 2016：SPA 核心章節補充 digest

本檔涵蓋偵測器、樣品、資料收集、refinement 與 validation。頁碼採各章印刷頁；涉及現行軟體選項時仍須另查官方文件。

## Direct electron detectors

來源：McMullan、Faruqi 與 Henderson，〈Direct Electron Detectors〉，pp. 1–17。

- 直接偵測器讓高能電子直接進入感測層，避開閃爍體與光學耦合造成的額外模糊；評估不能只看像素數，還要看 MTF、NPS 與 DQE 隨空間頻率的變化。
- Electron counting 可改善低劑量效率，但高事件密度會產生 coincidence loss，因此 dose rate 與 frame rate 是互相關聯的操作條件。
- 本站用途：`05_background` 的 DED／DQE 與 movie；歷史機型數字只作時代背景，不當成現行採購建議。

## Specimen behavior

來源：Glaeser，〈Specimen Behavior in the Electron Beam〉，pp. 19–50。

- 彈性散射提供結構對比；非彈性散射沉積能量並驅動輻射損傷。玻璃態冰只能減緩部分後續反應，不能消除損傷。
- 照射也會造成冰層變薄、起泡、集體位移與 doming；因此「低劑量」和 motion correction 是物理限制導出的設計，不只是軟體慣例。
- 本站用途：`05_background` 的 dose／damage／beam-induced motion，並提醒不同空間頻率的可用曝光量不同。

## Specimen preparation

來源：Passmore 與 Russo，〈Specimen Preparation for High-Resolution Cryo-EM〉，pp. 51–86。

- 高解析資料需要薄而均勻的玻璃態冰、完整且分散的粒子、足夠取向覆蓋，以及乾淨且機械穩定的 support。
- 樣品製備是逐步診斷流程：先確認生化品質與 negative stain，再測 vitrification、濃度、冰厚、表面吸附及資料收集條件。
- 本站用途：`05_background` 的樣品限制；不把「影像處理能修好」當成冰太厚、聚集或取向偏好的替代方案。

## Automated data collection

來源：Cheng 等人，〈Strategies for Automated CryoEM Data Collection Using Direct Detectors〉，pp. 87–102。

- 自動化資料收集必須平衡 throughput、漂移穩定、視野／孔洞選擇、dose fractionation 與資料品質監測。
- 更多 micrographs 不等於更多有效資訊；污染、厚冰、錯誤 defocus 或嚴重 motion 會使後段可用粒子率下降。
- 本站用途：`06_workflow` 的 acquisition 上游脈絡；具體控制軟體與硬體設定視為版本相關資訊。

## RELION and heterogeneity

來源：Scheres，〈Processing of Structurally Heterogeneous Cryo-EM Data in RELION〉，pp. 125–157。

- Empirical Bayesian regularization 以資料估計不同頻率的訊號強度，避免高頻雜訊在 refinement 中被無限制擬合。
- Gold-standard 做法把資料分成兩半並維持獨立 refinement；postprocessing 的 masking、FSC 與 sharpening 必須在此獨立性下解讀。
- 3D classification、focused mask 與 signal subtraction 可探查異質性，但結果依賴 mask、角度搜尋、粒子數及模型可辨識度。
- 本站用途：`06_workflow` 的 soft assignment、half-set、heterogeneity 與 postprocessing 限制。

## Refinement and variability

來源：Ludtke，〈Single-Particle Refinement and Variability Analysis in EMAN2.1〉，pp. 159–189。

- Refinement 是 projection matching 與 reconstruction 的迭代閉環；resolution、map accuracy、symmetry 與 variability 必須分開診斷。
- 對稱性可提高有效平均數，也可能抹去真正的不對稱訊號；局部柔性與成分變異會同時限制全域解析度。
- 本站用途：`06_workflow` 的 refinement 直覺與 variability；EMAN2.1 操作只作歷史方法案例。

## Frealign

來源：〈FREALIGN: An Exploratory Tool for Single-Particle Cryo-EM〉，pp. 191–226。

- FREALIGN 以帶 CTF 與雜訊權重的似然／相關目標精修取向、平移與其他參數，展示 forward model 如何進入實際 refinement。
- 參數搜尋、解析度範圍、mask 與 reference 都可能造成局部最優或 bias；逐輪提升解析度必須由獨立驗證約束。
- 本站用途：`06_workflow` 的 projection matching 實例，不把工具專屬參數提升成一般物理定律。

## Map validation

來源：Rosenthal，〈Testing the Validity of Single-Particle Maps at Low and High Resolution〉，pp. 227–252。

- 低解析度可用 tilt-pair 檢查取向與 handedness；高解析度需看獨立 half-map FSC、high-resolution noise substitution、mask 影響及 map-to-model agreement。
- 柔邊 mask 仍可能引入相關，FSC 門檻本身也不能證明局部特徵正確；應同時報局部解析度與可重現性檢查。
- 本站用途：`06_workflow` validation；任何「達到 X Å」的敘述都需附 FSC protocol 與限制。

