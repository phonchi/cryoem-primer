# 教學關係圖檢查

## 範圍與設計依據

新增geometry-family.html與wavelet-filterbank.html，不重繪或取代原始圖。實際讀取diagram-design SKILL.md、style-guide.md與type-nested.md、type-tree.md、type-process.md；geometry選Nested，wavelet選Tree（合成為反向樹）。

專案已明確核准白／淺底、ink #192c43、blue #1769aa、accent #cc6517、green #237858及system-ui/Noto Sans TC。按這份專案規格生成，沒有修改全域plugin/style-guide，沒有載入Google Fonts或其他外部資產。

640 viewBox，SVG主字20px／標題24px；手機SVG最小顯示寬384px，必要時在figure內橫向捲動，使字不小於12px。Geometry四層4nodes0arrows；wavelet analysis9nodes8arrows，synthesis7nodes6arrows；只遞迴cA1，回傳次序[cA2,cD2,cD1]。全部靜態，不增加動畫runtime。

## 驗證完成

兩圖 self_check.py 均通過。讀者代理完成 672／390 px 瀏覽器檢查與四張截圖：外部請求 0、SVG 文字溢出 0，字級至少 20／12 px；瀏覽器程序 exit 0。
