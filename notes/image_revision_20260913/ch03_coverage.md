# 第03章原始notebook對照
原稿：`../03_images_processing_3_zh.ipynb`，120 cells，索引0-based。原稿未改；原次序保留。`03_fourier.py`各cell有source-cells標記。0–1安裝／下載／ipywidgets改成本地資產、共用support及前端lab；2–119逐cell均有對應。
| 原cells | 主題／原圖 | 本站處置 |
|---|---|---|
| 2–23 | 10Hz／100Hz／2秒；DC、小9元素陣列、fftfreq排序 | 保留原數字、各操作；修正負sine內積，CFT推導收合 |
| 24–29 | 2D公式圖與128×128基底 | 保留所有3张原图；Python静态＋fourier_basis，係數大小與DC修正 |
| 30–40 | Gengar／Lucario／Pikachu FFT、幅度相位重建 | 原圖恢复，不使用camera |
| 41–44 | Lucario與Pikachu雙向幅度相位交換 | 恢復；不取abs扭曲空間signed值，修正文義 |
| 45–53 | Dragonite全尺寸Gaussian／Lucario11×11卷積 | 恢復原std1；核歸一化、circular與linear邊界區分；原示意圖 |
| 54–56 | Lucario 3×3 std3 timing及boxplot | 100次原比較；direct明確設定 |
| 57–77 | Gengar方形HPF／LPF，1..24完整montage | 48結果全部保留；signed HPF；41×41、共軛mask；filter_sweep |
| 78–85 | Gengar std3 Gaussian inverse與difference | epsilon1e-12；保持正負差值，說明偏差與浮點誤差 |
| 86–88 | Gengar 2×2box＋noise、unsupervised_wiener | 原noise比0.5；固定seed42；非手寫固定K替代 |
| 89–99 | moonlanding低通及98percentile峰抑制 | 兩種原操作；修共軛對稱及DC保護，K40中心改81×81 |
| 100–107 | Pikachu Gaussian／Laplacian各層與montage | 回到03原次序；明確API定義；blending只導讀原links |
| 108–112 | chess_football Harris | 原k=.001、1%threshold、紅色區域 |
| 113–117 | Hubble 500×500 LoG／DoG／DoH | 原3方法、threshold、sigma與圈色半徑 |
| 118–119 | 延伸閱讀 | 原OpenCV兩入口、書節、Packt保留，補API定位 |

## 每cell覆蓋
| cell | 類型 | 狀態 |
|---|---|---|
| 0–1 | setup | 本地資產與共用support取代安裝／下載；移除未使用helper |
| 2 | markdown | restored with correction |
| 3 | markdown | restored with correction |
| 4 | code | retained original operation/prose |
| 5 | markdown | retained original operation/prose |
| 6 | markdown | restored with correction |
| 7 | markdown | restored with correction |
| 8 | markdown | restored with correction |
| 9 | code | retained original operation/prose |
| 10 | markdown | restored with correction |
| 11 | markdown | retained original operation/prose |
| 12 | markdown | restored with correction |
| 13 | markdown | retained original operation/prose |
| 14 | markdown | restored with correction |
| 15 | markdown | restored with correction |
| 16 | code | restored with correction |
| 17 | markdown | restored with correction |
| 18 | code | retained original operation/prose |
| 19 | markdown | retained original operation/prose |
| 20 | code | restored with correction |
| 21 | code | retained original operation/prose |
| 22 | code | retained original operation/prose |
| 23 | markdown | retained original operation/prose |
| 24 | markdown | retained original operation/prose |
| 25 | markdown | restored with correction |
| 26 | markdown | restored with correction |
| 27 | markdown | restored with correction |
| 28 | markdown | restored with correction |
| 29 | code | restored with correction |
| 30 | markdown | restored with correction |
| 31 | code | restored with correction |
| 32 | markdown | retained original operation/prose |
| 33 | markdown | restored with correction |
| 34 | markdown | restored with correction |
| 35 | code | restored with correction |
| 36 | markdown | restored with correction |
| 37 | code | restored with correction |
| 38 | markdown | restored with correction |
| 39 | markdown | restored with correction |
| 40 | code | restored with correction |
| 41 | markdown | restored with correction |
| 42 | code | restored with correction |
| 43 | code | restored with correction |
| 44 | markdown | restored with correction |
| 45 | markdown | restored with correction |
| 46 | markdown | retained original operation/prose |
| 47 | markdown | retained original operation/prose |
| 48 | markdown | restored with correction |
| 49 | code | restored with correction |
| 50 | markdown | restored with correction |
| 51 | code | restored with correction |
| 52 | markdown | retained original operation/prose |
| 53 | code | restored with correction |
| 54 | markdown | restored with correction |
| 55 | code | restored with correction |
| 56 | code | restored with correction |
| 57 | markdown | restored with correction |
| 58 | markdown | retained original operation/prose |
| 59 | markdown | restored with correction |
| 60 | markdown | restored with correction |
| 61 | code | restored with correction |
| 62 | code | restored with correction |
| 63 | code | restored with correction |
| 64 | markdown | restored with correction |
| 65 | code | restored with correction |
| 66 | markdown | restored with correction |
| 67 | code | restored with correction |
| 68 | markdown | retained original operation/prose |
| 69 | markdown | restored with correction |
| 70 | markdown | restored with correction |
| 71 | markdown | restored with correction |
| 72 | code | restored with correction |
| 73 | markdown | restored with correction |
| 74 | code | restored with correction |
| 75 | markdown | restored with correction |
| 76 | code | restored with correction |
| 77 | markdown | restored with correction |
| 78 | markdown | restored with correction |
| 79 | markdown | restored with correction |
| 80 | code | restored with correction |
| 81 | markdown | restored with correction |
| 82 | code | restored with correction |
| 83 | code | restored with correction |
| 84 | markdown | restored with correction |
| 85 | markdown | restored with correction |
| 86 | markdown | restored with correction |
| 87 | markdown | restored with correction |
| 88 | code | restored with correction |
| 89 | markdown | restored with correction |
| 90 | markdown | restored with correction |
| 91 | code | restored with correction |
| 92 | markdown | retained original operation/prose |
| 93 | code | restored with correction |
| 94 | markdown | retained original operation/prose |
| 95 | code | restored with correction |
| 96 | code | restored with correction |
| 97 | code | restored with correction |
| 98 | markdown | restored with correction |
| 99 | code | restored with correction |
| 100 | markdown | restored with correction |
| 101 | markdown | restored with correction |
| 102 | markdown | restored with correction |
| 103 | code | restored with correction |
| 104 | markdown | restored with correction |
| 105 | code | restored with correction |
| 106 | markdown | restored with correction |
| 107 | markdown | restored with correction |
| 108 | markdown | restored with correction |
| 109 | markdown | restored with correction |
| 110 | markdown | restored with correction |
| 111 | markdown | restored with correction |
| 112 | code | restored with correction |
| 113 | markdown | restored with correction |
| 114 | markdown | retained original operation/prose |
| 115 | markdown | retained original operation/prose |
| 116 | code | retained original operation/prose |
| 117 | markdown | restored with correction |
| 118 | markdown | restored with correction |
| 119 | markdown | restored with correction |

## 既有網站額外內容
原例以外的camera實驗、CTF、Fourier slice、SPA連結與理解檢查全部移出此檔；SPA遷移由root負責。原本沒有的固定K掃描不取代unsupervised_wiener。示意圖以image_path載入ch03-cell28-1..3、46-1、47-1、101-1..2、114-1..2。
