# 對齊分類章的逐項科學修正

共 13 項具體修正。前面三章的科學修正、數學推導與來源核對見 validation.md。

## 1. 06_alignment_classification.md：修正成本式s是對稱階數

原文：目標解析度對應的角度間隔記為 $s$

修正：對稱階數記為 $n_{\mathrm{sym}}$

## 2. 06_alignment_classification.md：避免與空間頻率s混淆

原文：\frac{t^2\pi^2n^5m}{s}

修正：\frac{t^2\pi^2n^5m}{n_{\mathrm{sym}}}

## 3. 06_alignment_classification.md：補齊解析度與n的關係

原文：也就是隨粒子大小與目標解析度急速上升

修正：若粒子直徑為 $D$、目標解析度為 $d$，在相應取樣下約有 $n\simeq2D/d$。這個估算描述逐一比較候選的做法，顯示成本如何隨粒子大小與目標解析度上升

## 4. 06_alignment_classification.md：旋轉衰減是局部且具方向性

原文：角度誤差的效果可以換算成同一個尺度。半徑 $r$ 處，角度誤差 $\Delta\theta$ 造成的位移是 $r\,\Delta\theta$。對一顆半徑 65 pixels 的粒子，$1.7^\circ$ 的誤差在外緣相當於 $65\times0.0297\approx1.9$ pixels 的位移；把這個值代入上面的包絡，6 pixels 尺度的訊號只剩約 0.13。粒子愈大，同樣的角度誤差造成的邊緣位移愈大，對角度精度的要求也愈嚴。

修正：小的平面旋轉誤差會在離中心較遠的地方造成較大位移。半徑 $r$ 處，切向位移約為 $r\,\Delta\theta$。若角度誤差為零平均 Gaussian，標準差為 $1.7^\circ$，在 $r=65$ pixels 的外緣，切向位移標準差約為 $65\times0.0297\approx1.9$ pixels。這時沿外緣切向、6 pixels 尺度的細節，在小角度局部近似下保留約 0.13 的振幅；徑向細節的影響則不同。固定的角度偏轉只會旋轉影像，隨機角度誤差的平均才造成模糊。

## 5. 06_alignment_classification.md：補方向性近似的收合推導

原文：角度精度本身又由訊雜比決定。

修正：```{dropdown} 小角度旋轉為何造成方向性的模糊？
令 $J=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$，位置 $\mathbf x$ 的一階位移為 $\Delta\theta J\mathbf x$。若 $\Delta\theta\sim N(0,\sigma_\theta^2)$，局部 Fourier 模式的平均相位因子為

$$
\exp\![-2\pi^2\sigma_\theta^2(\mathbf k\cdot J\mathbf x)^2].
$$

它依賴位置與頻率方向。外緣切向頻率受影響最大；當 $\mathbf k$ 沿徑向，內積為零，一階近似沒有此衰減。它不能當成整張影像的各向同性平移包絡。
```

角度精度本身又由訊雜比決定。

## 6. 06_alignment_classification.md：補獨立參考假設

原文：其中 $\eta_k$ 的每像素變異數為 $\sigma^2/N_k$。

修正：其中 $\eta_k$ 的每像素變異數為 $\sigma^2/N_k$。此例假設參考由獨立影像建立，或比較時採 leave-one-out，使 $\varepsilon$ 與 $\eta_k$ 不相關。

## 7. 06_alignment_classification.md：期望距離不等於必然分派

原文：現在設某張影像真正屬於 B。它與 $s_B$ 的訊號差距為零，所以殘差期望是 $d\sigma^2+338$；它與 $s_A$ 的差距為 $\|s-s_A\|^2$，殘差期望是 $\|s-s_A\|^2+d\sigma^2+16.9$。只有在 $\|s-s_A\|^2>321$ 時，這張影像才會被正確指派到 B。換言之，兩個類別的真實訊號差距若不夠大，影像會系統性地流向粒子較多的那一類。於是小類別繼續失去成員、平均更吵、更難吸引成員，最終塌陷 {cite}`sorzano2010,yang2012`。

修正：若測試影像真正屬於小類別 B，當 $\|s_B-s_A\|^2<338-16.9=321.1$ 時，大類別 A 的期望平方距離反而較低。這顯示參考雜訊水準可能影響指派；單次含雜訊的比較仍會波動，不能由期望直接斷定每張影像的歸屬。若小類別因此流失成員，平均影像的雜訊會進一步增加 {cite}`sorzano2010,yang2012`。

## 8. 06_alignment_classification.md：區分模型參數與旋轉角

原文：p(y_i\mid z,\theta)

修正：p(y_i\mid z,\Theta)

## 9. 06_alignment_classification.md：明定模型參數記號

原文：以 $z$ 表示類別、平面內旋轉與平移的組合，

修正：以 $z$ 表示類別、平面內旋轉與平移的組合，$\Theta$ 表示類別平均、雜訊等模型參數，

## 10. 06_alignment_classification.md：區分分群與聯合對齊

原文：**γ-SUP** 從穩健統計切入。

修正：**γ-SUP** 對已對齊的影像做穩健分群，處理前面對齊步驟留下的錯位離群值；它本身不估計旋轉與平移來修正影像。

## 11. 06_alignment_classification.md：有限支撐的q條件

原文：q-Gaussian 的有限支撐

修正：原文選用 $q<1$ 的 q-Gaussian 所具有的有限支撐

## 12. 06_alignment_classification.md：不將gamma-SUP說成角度估計方法

原文：它們處理的是同一個耦合問題的不同弱點

修正：前三者處理對齊與分類流程中的不同問題，γ-SUP 則在已有對齊結果上處理離群值與穩健分群

## 13. 06_alignment_classification.md：移除不實的唯一公開全文說法

原文：這是本章唯一可公開取得全文的原始方法論文，

修正：這篇原始方法論文可公開取得全文，
