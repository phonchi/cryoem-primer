# 第02章局部核驗

未執行整章、同步ipynb或建置網站；只計算原固定參數的獨立子實驗。原notebook沒有修改。

```json
{
  "valid_lengths": [
    98,
    90
  ],
  "mean_square_center_and_window22": [
    1.0000000000000002,
    0.4444444444444444
  ],
  "binary_pattern_identical": true,
  "area_lt21_semantics": true,
  "watershed": {
    "labels": [
      1,
      2
    ],
    "areas": [
      77442.0,
      38910.0
    ],
    "perimeters": [
      5082.625863363164,
      3577.8613613753473
    ]
  },
  "page": {
    "yen": 121,
    "sauvola_shape": [
      191,
      384
    ],
    "black_tophat_exact": true
  },
  "Snorlax": {
    "rgba_shape": [
      256,
      256,
      4
    ],
    "pixelated": [
      128,
      128
    ],
    "hog_length": 2048
  },
  "retained_original_links": 28,
  "diamond_radius3_diagonal": true,
  "syntax_scope": "pass",
  "source_cell_ids": 133,
  "full_chapter_executed": false,
  "original_sha256": "31ee0a4c0d15a158273638017e1de2bdc4325cef04ff3aeea1fec6c2c2f0b1a9"
}
```

最初URL比較regex誤含Markdown分隔符，以及SPA搜尋誤中spacing；修正檢查器後全部通過，未因此改動教材。source-cells涵蓋133個保留cells，cell0僅安裝magic由專案環境取代，ledger完整記錄134cells。原圖渲染與browser lab由root整合驗收。

Snorlax局部核驗還核對原notebook cell125的saved text輸出：((128,128),(2048,))，與復原輸入/特徵維度完全一致。最後補入原code-comment對應的十一點平均、valid/same正文說明，未改數值。
