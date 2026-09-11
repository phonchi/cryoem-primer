# 解析度與驗證：哪些細節真的被資料支持

經過 {doc}`06_heterogeneity` 的分類與多輪重建後，我們還要判斷：密度圖的哪些細節有資料支持，哪些可能來自處理流程的偏差？

冷凍電鏡沒有辦法直接回答這個問題，因為真實結構未知。我們可以先檢查兩份獨立處理的資料能否得到一致的結果，並把一致性換算成一個頻率尺度。本章說明這個換算的推導、它成立的條件，以及哪些誤差無法靠一致性發現。

## 用偏差與變異分解估計誤差

從估計的角度看，重建結果 $\widehat V$ 與真實結構 $V$ 的差距可以分解成兩項。對每個 voxel（或每個 Fourier 係數），

$$
\mathbb E\left[\left|\widehat V-V\right|^2\right]
=\underbrace{\left|\mathbb E[\widehat V]-V\right|^2}_{\text{偏差平方}}
+\underbrace{\operatorname{Var}\left(\widehat V\right)}_{\text{變異}} .
$$

此處 $\operatorname{Var}(\widehat V)=\mathbb E|\widehat V-\mathbb E\widehat V|^2$，適用於實數或複數係數。第一項描述重複收集資料後，平均估計偏離真值的程度；第二項描述不同批資料之間的波動 {cite}`sorzano2022`。

增加資料常能降低變異，但偏差如何改變，要看估計方法與誤差來源；有些偏差會隨資料增加而縮小，有些則會持續存在。舉例來說，設每個 Fourier 係數的估計有固定偏差 $b=0.2$ 與變異 $\sigma^2/N$，其中 $\sigma^2=100$：

$N=100$ 時均方誤差為 $0.04+1=1.04$，幾乎全部來自變異；$N=10^4$ 時是 $0.04+0.01=0.05$，偏差已大於變異；$N=10^6$ 時是 $0.04+0.0001\approx0.04$，幾乎全部來自偏差。在這個固定 $b$ 的例子裡，增加資料會把變異壓低，卻不會消除預設的偏差。

Sorzano 等人強調，SPA 中一些被稱為 overfitting 的現象來自系統性偏差：例如以 100,000 顆 $200\times200$ 粒子重建一個 $200^3$ 的體，表面量測數雖遠多於 voxel 數，角度與類別估計錯誤、成像模型不符，以及演算法偏差仍可能產生錯誤結構 {cite}`sorzano2022`。這是作者用來提醒讀者檢查參數與模型的觀點。

量測數多，並不能單獨保證可識別或排除過擬合：投影可能冗餘、方向覆蓋可能不足，角度與類別也需從含雜訊資料估計。驗證因此要同時檢查參數、資料獨立性與重建結果；圖像銳利或量測數較多，都不足以判定結果可靠。

(halfset)=
## 獨立 half-set：把獨立性放進流程

標準做法是在 refinement 之前把粒子隨機分成兩個 half-set，兩組從此各自獨立估計角度與位移並重建，得到 half-map A 與 half-map B。資料流可以寫成：

**particle stack → 分成兩組 particle IDs → 各自估計角度與位移並重建 → half-map A／B → FSC**。

這個安排要處理的具體問題是 {doc}`06_alignment_classification` 提過的 noise alignment：把影像對齊到由同一批影像算出的參考時，每張影像的雜訊參與決定自己的對位參數，使資料看起來比實際更一致 {cite}`penczek2010resolution`。兩組資料若各自獨立地經歷這個過程，它們各自「學到」的是各自的雜訊，這些雜訊不會互相一致。

「Gold-standard」指的就是用資料分割降低高頻過擬合的這個原則。RELION 的實作在每輪迭代中維持兩個獨立的重建，並用兩者的 FSC 決定下一輪要使用到多高的頻率；早期把整批資料對同一個參考精修、事後才切半計算 FSC 的做法，會讓 FSC 曲線被人為抬高 {cite}`scheres2016`。

要維持獨立性，必須注意哪些東西可以共享、哪些不行。事先固定的低解析度初始模型與共同的處理規則通常可以共享；高解析度的 map、由其中一半的高解析度資料估出的角度與位移，以及從完整資料的高解析度細節製作的緊密 mask，都會把共同結構帶到另一半。資料層面也要小心：同一個物理粒子的重複觀測會共享資訊；同一部 movie 或同一張 micrograph 中的粒子，也可能因共同背景或校正誤差而相關。{doc}`05_background` 已經算過，雜訊之間只要有 1% 的相關，一萬張影像的效果就只相當於約一百張獨立影像。雜訊相關的程度，會直接影響平均能改善多少。

(fsc)=
## FSC 與頻譜訊雜比

取包含正負共軛頻率的完整對稱 shell，半徑記為 $s$。Fourier shell correlation 定義為

$$
\mathrm{FSC}(s)=
\frac{\sum_{\mathbf k\in s}F_A(\mathbf k)F_B(\mathbf k)^*}
{\sqrt{\sum_{\mathbf k\in s}|F_A(\mathbf k)|^2\;
\sum_{\mathbf k\in s}|F_B(\mathbf k)|^2}} .
$$

分子是兩張圖的交叉功率，分母把各自的功率正規化掉。因此 FSC 量的是線性一致性：若兩組係數只差一個共同的正比例常數，FSC 仍等於 1。

要把 FSC 換算成訊雜比，需假設兩張 half-map 共享訊號，而各自的雜訊零平均、彼此獨立，且具有相同功率。令 $P_S$ 為訊號功率、$P_N$ 為單張 half-map 的雜訊功率，母體相關為

$$
\rho(s)=\frac{P_S}{P_S+P_N}.
$$

FSC 是用有限 shell 係數估計 $\rho$ 的統計量，並非每次都等於這個功率比。有效樣本足夠時，才可用 FSC 近似 $\rho$，得到

$$
\mathrm{SSNR}_{\mathrm{half}}(s)=\frac{\rho}{1-\rho},
\qquad
\mathrm{SSNR}_{\mathrm{full}}(s)=\frac{2\rho}{1-\rho}.
$$

第二式對應兩張 half-map 的平均：獨立且等功率的雜訊經平均後，功率減半。真實 SSNR 由非負功率定義。有限樣本的 FSC 卻可能因抽樣波動而為負；這些頻帶不宜直接套入訊雜比換算 {cite}`penczek2010resolution`。

```{dropdown} 從「訊號加雜訊」推出 FSC 與 SSNR 的關係
把某個 shell 內的係數寫成訊號加雜訊：

$$
F_A=S+N_A,\qquad F_B=S+N_B .
$$

假設 $S$ 是兩張圖共有的真實訊號，$N_A$ 與 $N_B$ 平均為零、互相獨立、也與 $S$ 獨立，且在該 shell 內具有相同的功率 $P_N$。記訊號功率為 $P_S$。取期望值：

$$
\mathbb E\left[\sum F_AF_B^*\right]=n_s P_S,
\qquad
\mathbb E\left[\sum|F_A|^2\right]
=\mathbb E\left[\sum|F_B|^2\right]=n_s\left(P_S+P_N\right),
$$

其中 $n_s$ 是 shell 內的係數個數；交叉項的期望因獨立與零平均而消失。用這些母體功率定義的相關為

$$
\rho(s)=\frac{P_S}{P_S+P_N}.
$$

FSC 的分子與分母在有效樣本增加時可接近各自期望，因而近似上式；這不是有限樣本下「比值的期望等於期望的比值」。定義單張 half-map 的頻譜訊雜比 $\mathrm{SSNR}_{\mathrm{half}}=P_S/P_N$，把上式反解：

$$
\mathrm{SSNR}_{\mathrm{half}}(s)=\frac{\rho}{1-\rho}
$$

完整的 map 是兩張 half-map 的平均，$F=\tfrac12(F_A+F_B)=S+\tfrac12(N_A+N_B)$。雜訊功率變成 $\tfrac14(P_N+P_N)=P_N/2$，訊號功率不變，所以

$$
\mathrm{SSNR}_{\mathrm{full}}(s)=\frac{2\rho}{1-\rho}
$$

也就是完整 map 的頻譜訊雜比恰為單張 half-map 的兩倍——正是把兩組獨立、等功率雜訊平均的結果 {cite}`penczek2010resolution`。

這個推導用到四個假設：兩張圖共享同一個訊號、雜訊零平均、雜訊互相獨立、兩邊雜訊功率相同。任何一項不成立，上面的換算都要重新檢視。特別是「共享同一個訊號」這一條——下面會看到，它比表面上寬鬆得多。
```

可以用這個關係解讀常見門檻。$\rho=0.5$ 對應 $\mathrm{SSNR}_{\mathrm{full}}=2$；$\rho=1/3$ 對應 $\mathrm{SSNR}_{\mathrm{full}}=1$，也就是完整 map 在該 shell 的訊號功率恰等於雜訊功率 {cite}`penczek2010resolution`。

那 0.143 從哪裡來？它來自另一個量。在同一訊號加雜訊模型下，完整 map 與無誤差參考結構之間的母體相關為

$$
C_{\mathrm{ref}}(s)=\sqrt{\frac{2\rho(s)}{1+\rho(s)}},
$$

以樣本 $\mathrm{FSC}=0.143$ 近似 $\rho$：$2\times0.143=0.286$，除以 $1.143$ 得 $0.2502$，開根號得 $0.5002$。精確令 $C_{\mathrm{ref}}=0.5$ 時，對應的是 $\rho=1/7\approx0.143$。因此慣用的 0.143 門檻，約對應完整 map 與真實訊號的相關為 0.5 {cite}`rosenthal2016`。這是相關量的換算，不能由它推出平均角度或平均相位誤差。

這說明 0.143 是一個經過換算的慣例，它把「half-map 之間的一致性」對應到「與真值的一致性」這個比較容易解釋的尺度。它本身沒有構成統計檢定。若要問某個 FSC 值是否顯著，需要知道 FSC 的抽樣分布，而相關係數的分布並非常態，可用 Fisher 的 $z$ 變換作近似，而且 shell 內獨立係數的個數會因 mask、對稱與內插而遠少於 voxel 數 {cite}`penczek2010resolution`。同樣重要的是，門檻的適用性建立在前面那四個假設上；當 half-set 不完全獨立或 mask 引入共同結構時，不管挑哪個數字都會偏樂觀。

```{figure} images/pptx/s24_2.png
:width: 72%
:name: fig-fsc-unmasked

這條 unmasked FSC 曲線在約 $0.267\ \mathrm{\mathring A^{-1}}$ 穿過 0.143，對應解析度
$1/0.267\approx3.74\ \mathrm{\mathring A}$。讀曲線時先找交點，再檢查高頻是否異常上升；加入 mask 後，還要把
masked 與 unmasked 曲線並排比較。
```

## 共同偏差：高 FSC 仍可能有錯誤

回到上面推導中的第一個假設：兩張圖共享同一個訊號 $S$。推導只要求 $S$ 是兩邊共有的成分，並沒有要求它等於真實結構。

把共同偏差明確寫出來。設

$$
F_A=S_{\mathrm{true}}+b+N_A,
\qquad
F_B=S_{\mathrm{true}}+b+N_B,
$$

其中 $b$ 是兩邊共有的系統性誤差。整個推導完全不變，只是原來的 $P_S$ 現在變成 $S_{\mathrm{true}}+b$ 的功率。FSC 無法把共同成分中的真實訊號與偏差分開。因此即使偏差很大，FSC 仍可能很高。不過，偏差不必然提高 FSC：若 $b=-S_{\mathrm{true}}$，共同訊號反而被抵消。

什麼會造成共同偏差？兩個 half-set 若從同一個有偏的初始模型出發、被施加了錯誤的對稱、含有同一批誤挑的粒子、用同一個不當的 mask，或使用同一個帶有偏差的目標函數，兩邊就會得到相同或相似的偏差 {cite}`sorzano2022`。Sorzano 等人因此主張 gold standard 對於排除偏差**既不充分也不必要**：不充分，因為兩半很容易落入同一種偏差；不必要，因為在迭代中系統性交換兩半影像的流程也未必出現過擬合的徵兆 {cite}`sorzano2022`。

上述是作者對偏差的論點；half-set 在檢查雜訊過度擬合時仍有用途。保持兩組獨立仍能降低共享高頻雜訊造成的樂觀估計；共同偏差則需要改變初始模型、方法或外部證據來檢查。

## Mask 的兩面

Mask 排除大面積溶劑，可降低雜訊功率、改善數值穩定性，也讓 FSC 更能反映結構所在區域的品質。但它在 Fourier 空間的作用需要留意。

實空間的相乘對應 Fourier 空間的卷積：把 map 乘上 mask $m$，等於把 $\widehat V$ 與 $\widehat m$ 卷積。於是原本分屬不同 shell 的係數被混在一起，$\mathrm{FSC}(s)$ 不再只描述半徑 $s$ 的資訊。

混合的方向很不對稱，這正是問題所在。低頻的係數強度大、兩張 half-map 之間本來就高度相關；高頻的係數弱、相關性低。卷積把強而相關的低頻成分摻進弱而不相關的高頻 shell，於是 FSC 在高頻上升，而這個上升並不反映該頻帶真正的訊雜比 {cite}`scheres2016`。mask 的細節愈多，$\widehat m$ 的能量愈分散，這個效應愈強；柔邊（例如寬度數個 pixel 的餘弦邊界）把 $\widehat m$ 集中起來，可以壓低混合與振鈴 {cite}`scheres2016`。

事先固定的同一個 mask，並不會使原本獨立的兩組零平均雜訊自動產生交叉相關；前面的高頻上升主要來自共同訊號的頻率混合。另一個問題是 mask 的形狀由資料本身決定。用完整資料的高解析度細節做出緊貼分子的 mask，再把同一個 mask 套回兩張 half-map，等於給兩邊加上同一個帶有高頻細節的形狀——按照上一節的分析，這是一個標準的共同偏差。

可用高頻相位隨機化估計 mask 造成的頻率混合。選定截止頻率後，對兩張 half-map 在該頻率以上的相位**分別獨立隨機化**，並保留 Hermitian 對稱，使反轉換仍得到實數密度。再套用同一個 mask，得到 $\mathrm{FSC}_{\mathrm{randomized}}$。在截止頻率以上，尚未隨機化的共同低頻訊號仍可經 mask 混入高頻，因此這條曲線可估計該效應；有限樣本也會使曲線波動。

校正式還需要用分母正規化：

$$
\mathrm{FSC}_{\mathrm{corrected}}
=\frac{\mathrm{FSC}_{\mathrm{masked}}-\mathrm{FSC}_{\mathrm{randomized}}}
{1-\mathrm{FSC}_{\mathrm{randomized}}}.
$$

例如 masked FSC 為 0.5、randomized FSC 為 0.2，校正後是 $(0.5-0.2)/(1-0.2)=0.375$。[Rich Hite 的公開講義第 11 頁](https://semc.nysbc.org/wp-content/uploads/2018/08/20190311_SEMCcourse_SPAIII-RH.pdf)列出這個公式。它應用在相位已隨機化且避開截止頻率過渡影響的頻帶；較低頻仍使用未隨機化的結果。若 randomized FSC 接近 1，分母很小，校正會不穩定，應重新檢查 mask 與截止頻率。

還要分清楚兩種檢查的目的：在**粒子影像**層級做高頻 noise substitution 並重跑 refinement，可以測試演算法是否把無訊號的頻率學成共同結構；在 **half-map** 層級做相位隨機化，主要用來估計後處理 mask 的頻率混合。後者不能取代前者的流程驗證，但對 mask 校正仍有明確用途 {cite}`sorzano2022,scheres2016`。

將 masked、unmasked、randomized 與 corrected 曲線一起看，才能區分去除溶劑雜訊後的改善與 mask 引入的效應。校正並不保證消除所有共同偏差。

## 一個數字描述不了整張圖

整體 FSC 給出一個數字，但密度圖各處的品質通常差很多。

**局部解析度**在小的局部區域內重複做同樣的一致性估計，得到逐位置的解析度圖。分子核心通常比週邊柔性區域好，而 {doc}`06_heterogeneity` 說明的局部模糊會直接顯示在這張圖上。局部估計的代價是每個區域可用的 Fourier 係數變少，估計本身的不確定性上升。

**方向解析度**沿不同方向分別計算一致性。{doc}`06_reconstruction_validation` 說明過，偏好取向會讓三維 Fourier 空間某些方向取樣稀疏；球面平均的 FSC 把這個各向異性平均掉。頻譜訊雜比原則上可以逐 voxel 估計，據此描述各向異性 {cite}`penczek2010resolution`。方向分布圖、方向性 FSC 與局部解析度圖回答的是三個不同的問題，應該分別保存。

(sharpening)=
## 銳化：誰的振幅被放大了

重建圖的高頻振幅通常低於真實結構，原因包括輻射損傷、殘餘運動、對位誤差與包絡衰減。銳化把這一段補回來，方式是乘上一個隨頻率上升的因子：

$$
\widehat V_{\mathrm{sharp}}(s)=\exp\!\left(-\frac{B\,s^2}{4}\right)\widehat V(s),
\qquad B<0 .
$$

需要這個步驟，是因為高頻成分在成像、偵測與影像處理三個環節都被壓低了 {cite}`scheres2016`。Guinier plot 用平方空間頻率 $s^2$ 作橫軸、徑向平均振幅的自然對數 $\log a(s)$ 作縱軸。在選定的可靠頻帶內，若以 $a(s)\approx a_0\exp(-B_{\mathrm{damp}}s^2/4)$ 描述衰減，直線斜率就是 $-B_{\mathrm{damp}}/4$；補償時，銳化式採 $B=-B_{\mathrm{damp}}$。先檢查擬合頻帶、雜訊及結構本身的頻率變化，才能把斜率解讀為衰減。比較不同 frames 時也可畫對數相對振幅，這與{doc}`06_workflow`的相對 B factor 相接；若縱軸改用對數功率，斜率會多一個 2 的因子。取 $B=-100\ \mathrm{\mathring A^2}$，在 10 Å（$s=0.1$）處放大倍率是 $\exp(100\times0.01/4)=\exp(0.25)\approx1.28$；在 3 Å（$s\approx0.333$）處是 $\exp(100\times0.1111/4)=\exp(2.78)\approx16.1$。

**銳化因子同時乘在該 shell 的訊號與雜訊上**。3 Å 的訊號被放大 16 倍，同一個 shell 的雜訊也被放大 16 倍。銳化改變的是振幅的頻率分布，訊雜比逐頻率維持不變。它讓已經存在的高頻特徵變得容易看見，也同樣讓已經存在的高頻雜訊變得容易看見。

因此銳化通常配合頻率相依的加權，並限制在有可靠訊號的範圍。在 $F=S+N$ 模型下，若目標是估計 map 中的訊號 $S$，最小均方誤差的線性增益為：

$$
\frac{\mathrm{SSNR}_{\mathrm{full}}}{\mathrm{SSNR}_{\mathrm{full}}+1}
=\frac{2\rho}{1+\rho}
$$

{cite}`penczek2010resolution`。實作時可用 FSC 近似 $\rho$；弱訊號頻帶的權重較小，但它不會正好抵消任意銳化因子。這個增益也等於同一模型下的 $C_{\mathrm{ref}}^2$。

如果目標改成包絡衰減前的真實係數 $T$，模型就應寫成 $F=E T+N$。已知包絡 $E$，且 $T,N$ 不相關時，估計 $T$ 的線性增益為

$$
\frac{E^*P_T}{|E|^2P_T+P_N}.
$$

它還需要包絡與衰減前的訊號功率，單一 $2\rho/(1+\rho)$ 因子不足以決定這個增益。這也說明為何銳化、弱頻帶收縮及可靠頻率上限要一起判讀。

## 重現性與共識

最後回到一開始的問題：如何檢查那些 FSC 看不見的偏差？

單一次執行給出單一組參數估計，我們無法判斷它的誤差屬於「圍繞真值的隨機波動」還是「系統性地偏離真值」；即使拿到多組估計，也無法個別判斷哪一組屬於哪一種。但有多組估計時（至少兩組），可以做一件事：比較它們是否在參數空間的同一個區域取得共識。若估計彼此接近，表示結果在這些執行條件下穩定，但也可能共享偏差；若估計分散，則需要檢查初始化敏感性、資料不足或模型不符。將相容的估計平均可能降低變異，效果還取決於估計間的相關性 {cite}`sorzano2022`。

實際分析時可以這樣檢查。換一個初始模型、換隨機種子、換 half-set 的切分方式重跑，比較主要密度特徵是否重現。用不同的軟體或不同的演算法處理同一批資料，比較逐粒子的角度指派是否一致——在參數層級比較，比在最終 map 層級比較更靈敏，因為偏差在被投影進三維之前更容易被辨認 {cite}`sorzano2022`。此外，投影與反投影的殘差若帶有系統性結構，通常表示模型還漏掉影像中的某些規律。

作者也提醒，兩個一致的估計仍可能同時有偏 {cite}`sorzano2022`。共識提高信心，卻無法提供保證。最後一層檢查來自模型之外——密度中的組成與構形能否與生化、功能或既有結構資訊互相印證。

## 接下來

本篇的 SPA 內容到此完成：從樣品與劑量的限制、影像形成、統計推論，到對位、分類、降維、重建、異質性與驗證。可以用已知真值的資料逐項改變條件，觀察每個估計步驟如何反應。{doc}`07_synthetic_data` 提供可執行的模擬：它按照 {doc}`05_image_formation` 的前向模型產生帶有已知取向、平移、CTF 與雜訊的粒子影像，並保存全部真值供比對。

合成章的置中平均，是在同一批角度與 CTF 下，比較已知位移校正與零位移對照，藉此單獨觀察平移的效果。它尚未估計旋轉、完成分類或做三維重建；讀圖時要把這項診斷與前面完成角度對齊的 class average 分開看。

## 延伸閱讀

- **Penczek, P. A.（2010）．*Resolution Measures in Molecular Electron Microscopy*. Methods in Enzymology, 482, 73–100。** [DOI](https://doi.org/10.1016/S0076-6879(10)82003-8)。本章 FSC 與 SSNR 關係的主要來源：第 3 節（頁 78–83）比較 Q-factor、DPR、FRC 與 SSNR 的定義並說明兩種 half-set 建立方式（頁 81–82 討論 noise alignment）；第 4 節（頁 83–86）給出單張 half-map 的 $\mathrm{SSNR}=\mathrm{FSC}/(1-\mathrm{FSC})$ 與兩半平均之 full-map 的 $2\mathrm{FSC}/(1-\mathrm{FSC})$、各種門檻對應的 SSNR 值、以 Fisher $z$ 變換近似相關係數抽樣分布的做法，以及以 FSC 表示的最佳線性濾波增益 {cite}`penczek2010resolution`。
- **Sorzano, C. O. S., Jiménez-Moreno, A., Maluenda, D. 等（2022）．*On Bias, Variance, Overfitting, Gold Standard and Consensus in Single-Particle Analysis by Cryo-Electron Microscopy*. Acta Crystallographica Section D, 78, 410–423。** [開放全文（IUCr）](https://doi.org/10.1107/S2059798322001978)。第 2 節（頁 410–411）建立偏差－變異的框架並論證 EM 的 overfitting 主要是偏差；第 3 與第 4 節（頁 411–419）逐項列出實驗與演算法造成偏差的來源；第 5.1 節討論以 FSC 偵測 overfitting 的限制，第 5.2 節（頁 419–420）給出「gold standard 既不充分也不必要」的論證與共同偏差的清單，第 5.3 節指出相位隨機化在影像層級與體層級實作的差別；第 6 節（頁 420）說明多次估計取得共識的策略及其極限。此文是本章共同偏差討論的主要參考，可公開取得全文 {cite}`sorzano2022`。
- **Rosenthal, P. B.（2016）．*Testing the Validity of Single-Particle Maps at Low and High Resolution*. Methods in Enzymology, 579, 227–252。** [DOI](https://doi.org/10.1016/bs.mie.2016.06.004)。第 3.2 與 3.3 節（頁 244–245）說明相位隨機化與其他雜訊檢定的流程；第 4 節（頁 245–246）定義 $C_{\mathrm{ref}}=\sqrt{2\mathrm{FSC}/(1+\mathrm{FSC})}$，並將 half-map FSC 的約 0.143 對照到 $C_{\mathrm{ref}}$ 的約 0.5；閱讀時要分清相關係數與平均相位誤差，不能直接把相關值取反餘弦當成平均誤差 {cite}`rosenthal2016`。
- **Scheres, S. H. W.（2016）．*Processing of Structurally Heterogeneous Cryo-EM Data in RELION*. Methods in Enzymology, 579, 125–157。** [DOI](https://doi.org/10.1016/bs.mie.2016.04.012)。第 2.2 節（頁 129–131）說明 Bayesian 方法中先驗功率由重建自身估計為何可能累積雜訊、gold-standard refinement 如何處理它，以及 0.143 門檻在獨立重建下才適用的來由；第 2.3 節（頁 131–133）是本章 mask 一節的來源，說明卷積如何把強而相關的低頻摻進高頻、相位隨機化校正的完整步驟，以及銳化為何被放在同一個 postprocessing 階段；柔邊寬度 3–10 pixels 的實用範圍出現在第 4.4 節（頁 145）{cite}`scheres2016`。
- **Sigworth, F. J.（2016）．*Principles of Cryo-EM Single-Particle Image Processing*. Microscopy, 65(1), 57–67。** [Oxford Academic 免費全文](https://doi.org/10.1093/jmicro/dfv370)。頁 64–65 說明訊雜比為何是 SPA 的根本限制：它同時決定取向估計的精度與能否分辨混合樣品中的不同族群，而更大的資料量與統計方法只能部分補償；同一段也給出當時電子計數相機的 DQE 範圍（低頻約 0.7、在常用解析度極限附近約 0.4）。篇幅短、公開取用 {cite}`sigworth2016`。
- **[EMDB](https://www.ebi.ac.uk/emdb/) 的公開條目。** 先查看條目實際提供哪些 FSC 曲線、mask、解析度與處理流程中繼資料；不同條目的資料完整度可能不同。挑幾筆自己感興趣的結構，對照本章的檢查項目逐一閱讀它們的驗證資訊，可以練習判讀實際資料的驗證結果。

- **Hite, R.（2019-03-11）．*Single-particle analysis (Part IV): Cryo-EM map interpretation*，SEMC course。** [公開講義](https://semc.nysbc.org/wp-content/uploads/2018/08/20190311_SEMCcourse_SPAIII-RH.pdf)。第 11 頁列出以相位隨機化 FSC 校正 masked FSC 的公式；搭配本章的 0.5、0.2 數值例子閱讀。
