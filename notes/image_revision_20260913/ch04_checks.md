# 第 4 章必要檢查

以 `/home/phonchi/miniconda3/envs/cryoem-book/bin/python` 在記憶體進行。

- 原38-cell source ID全部存在；8張原示意圖均有image_path，唯一實驗圖片為原Lucario。
- Python AST與Jupytext解析通過；72個新cells的全部code cells可編譯。
- db5/sym5/coif5 wavefun回傳3項，psi為第2項；bior2.4回傳5項，使用psi_d與psi_r。cgau5為複數，因此分開顯示實虛部。
- db1..db5 × wavefun level1..5共25組，全部取得函數值。
- ECG單階db1 smooth還原最大誤差5.6843e-14；db1 level8還原最大誤差1.7053e-13。
- Soft threshold [-2,-.5,0,.5,2]，threshold1得到[-1,0,0,0,1]。
- ECG/256、Gaussian SD .05、rng42：啟發式閾值 .1*max(noisy)=0.1062745786；MSE noisy=0.002419444、manual db4=0.000924711、BayesShrink db4 level4=0.000642114。此為固定例子的檢查，不宣稱普遍排序。
- 沒有SPA/MRC/理解檢查及指定禁用句式。

未執行全章或Lucario完整figure pipeline；原圖實際顯示、shared support、notebook同步與建置交root統一驗證。原notebook未改。
