# 第03章檢查紀錄

此agent只負責03.py及兩份ledger，不sync notebook、不build、不下載素材、不執行全章。

已實作的檢查：DFT/IFFT誤差、共軛對稱與IFFT虛部、DC cosine/sine、phase交換原尺寸、Gaussian核歸一化、direct/FFT相等、41×41遮罩、各cutoff low+high重建。這些完整執行由root整合後驗證，尚不可標成通過。

必要改正：DFT負sine、magnitude/phase非real/imag、FFT卷積成本、direct明設、signed高通、inverse epsilon偏差、對稱notch、pyramid_laplacian非跨層residual。保留原重要圖像、操作、數字與threshold；增加固定seed42以便重現。

文字依speak-human-tw直接潤稿：移除籠統的「相位才有空間資訊」「FFT一定較快」等句；補明操作與假設；原完整來源links保留。原稿不修改。

## 已執行的最低必要檢查

- `ast.parse` 語法檢查通過。
- 僅抽出 `real_ifft2/frequency_indices/square_lowpass/gaussian_psf/coefficient_basis`，不載入資產、不呼叫plot、不執行全章。
- 128×128的cosine/sine，測(5,2)、DC、三種Nyquist自共軛點及(-32,32)：FFT往返、實數性通過。
- 256×256與251×320各測cutoff1..24：low/high共軛對稱、虛部容差、low+high回到原圖均通過。
- 256×256 Gaussian std3與epsilon1e-12的反演數值：虛部／共軛檢查通過。
- 未執行素材圖像、unsupervised_wiener、瀏覽器或完整notebook；交root整合驗證。

最終static檢查：120/120 source-cells均有映射；9個原示意圖URL以source-image-url保留；無camera／CTF／SPA／理解檢查／Fourier slice／Colab絕對路徑／執行時下載命令；AST通過。
