# FSJ 定義鏈查證（只讀；未改程式與論文）

證據來源縮寫：
- `circ.py` = server `~/miniconda3/envs/ciriquant/lib/python3.7/site-packages/CIRIquant/circ.py`（**實際被 import 的版本，version.py = 1.1.3**）
- `analysis.R` / `merge_counts.py` = server `~/circRNA_agent/scripts/…`（offset 區塊自 2026-05-26 起未變：`git log -S fsj_pseudo` → 20c206e；2026-10-02 的 2f27d7b 只改了 filterByExpr 的參數寫法）
- 重算腳本：`audit/f1.R f3.py f4.R f5.py f6.py f7.py f8.R`（輸出只寫 `~/audit_tmp/`）

> **先更正我自己先前的一個發現（C040 / C277 / C328）**：我之前說「CIRIquant 實際是 1.1.2、論文 §2.3 寫 1.1.3 是錯的」。這是錯的。
> server 的 conda env 同時裝了兩份：`CIRIquant-1.1.2-py3.7.egg`（Python 2 語法，如 `izip_longest`、`cPickle`，在 py3.7 下不能執行）與 conda 安裝的 `ciriquant-1.1.3.dist-info` + `site-packages/CIRIquant/`（`version.py` = '1.1.3'）。
> `CIRIquant --version` 輸出 `CIRIquant 1.1.3`；`import CIRIquant` 解析到 1.1.3 目錄；`pip show` 也是 1.1.3。`conda list` 顯示 `1.1.2=pypi_0` 是 egg 殘留造成的誤導（CLAUDE.md 的「更正」也是被它誤導）。
> 因此：**§2.3 的 1.1.3 正確；Table 1 與 Fig 6 圖說的 1.1.2 才是錯的。** 另外 1.1.3 目錄的 `direct_url.json` 指向 `file:///opt/conda/conda-bld/ciriquant_1749565571344/work`（本地 conda-build，時間戳換算約 2025-06-10），程式碼含 Python 3 移植修改與中文註解（`circ.py` 內「將 chrom_start 和 chrom_length 轉換為整數」）——它是**本地建置、含 patch 的 1.1.3**，不是原封不動的上游發行版（論文未交代，S3）。`claims.tsv` 已據此更正。

---

## 1. CIRIquant 端：fsj 到底是什麼

### 1.1 計的是哪些 reads（`circ.py:530-550`，`genome_worker`）
對每個 circRNA（起點 `parser.start`、終點 `parser.end`）分別掃兩個邊界：
- 起點：`sam.fetch(region=chrom:start-start)`（:532），保留 read 若：未 unmapped、非 supplementary（:533）、`mapping_quality > 10`（:535）、`read.get_overlap(start-1, start+THRESHOLD-1) >= THRESHOLD`（:537，即 read 必須覆蓋起點起算的 5 個鹼基；THRESHOLD = `-a/--anchor`，預設 5，`main.py:35-36`）、且 read 兩端 CIGAR 皆「mapped」（`is_mapped`，:539；`circ.py:569-586`：軟剪裁 ≤10 nt 才算 mapped）。
- 終點：同上，視窗 `get_overlap(end-THRESHOLD, end)`（:542-550）。
- 兩個邊界的檢查**沒有任何剪接條件、也沒有 strand 條件**（整段 :530-550 不引用 `read.is_reverse`、`parser.strand`）。所以它計的是「**在基因組比對 BAM 中，覆蓋 circRNA 邊界外顯子端點 5 nt 以上、且末端沒被大幅剪裁的線性 reads**」——可以是跨剪接的 reads，也可以是未剪接的（內含子保留/pre-mRNA）reads。名稱「forward-splice junction」是 CIRIquant 沿用的簡稱，程式並未要求 read 真的跨過剪接點。

### 1.2 兩個邊界怎麼合併：既不是平均、也不是總和/最大值——是「去重後的聯集」
- 兩個邊界的 reads 都 append 到同一個 list `fsj_reads`（:540、:550）。
- `proc_genome_bam`（`circ.py:444-445`）：`fsj_reads[pair_id][mate_id] = circ_id`——每個 (pair, mate) **只保留一個 circ_id（後寫覆蓋前寫）**。
- `circ.py:456-460`：`circ_fsj[circ_id].update({query_prefix(pair_id): 1})`——以「去掉 /1 /2 的 pair id」為鍵，所以同一個 fragment 的兩個 mate、或同時碰到起點與終點，只算 1。
- `circ.py:704`：`fsj = len(fsj_reads[circ_id])`。
→ **fsj = 在起點或終點邊界上有線性 reads 的相異 fragment 數（聯集，fragment 層級）**，且排除同時是 BSJ 候選 read 的 reads（:458-459）。
→ **重要副作用（本次新發現）**：同一條 read 若同時覆蓋多個 circRNA（共享邊界或邊界相近）的錨點，只會被記給其中一個 circ（:445 覆蓋）。共享起點/終點的同基因 circRNA 異構體，其 FSJ 會被任意地「分配」。見 §5。

### 1.3 在哪個階段算
- 在 `proc()`（`circ.py:628-671`）內，順序是：產生偽環狀參考（:650）→ 建 HISAT2 index（:653）→ 偽環狀參考重比對得 BSJ 候選（:657-661）→ **`proc_genome_bam(hisat_bam, …)`（:662）：用步驟前期 HISAT2 對基因組的比對 BAM** 掃 FSJ。
- 所以：**BSJ = 偽環狀參考重比對（再以基因組 BAM 複核 false-positive，:510-524）；FSJ = 基因組比對 BAM**。兩者來自不同的比對。

### 1.4 fsj = 0 的條件
- `expression_level`（`circ.py:697-709`）：`fsj = len(fsj_reads[circ_id]) if circ_id in fsj_reads else 0`；只有 `bsj==0 且 fsj==0` 才整列跳過（:705）；因此 GTF 內有「bsj>0、fsj=0」的列。
- fsj=0 發生於：(a) 邊界處沒有通過 MAPQ>10/錨長/剪裁條件的線性 read；(b) read 被 :445 的單值覆蓋「讓給」共享邊界的另一個 circ；(c) 使用 `--no-fsj`（`main.py:71`；管線**未**使用，`circrna.smk:50-59`）。
- 輸出到 GTF：`bsj {:.3f}; fsj {:.3f}; junc_ratio …`（`circ.py:786-791`）；`junc_ratio = 2·bsj/(2·bsj+fsj)`（:708）——管線不用這個比值。

## 2. 管線端

### 2.1 `merge_counts.py`（server `scripts/merge_counts.py`）
- `parse_gtf`（:17-43）：用正規表示式 `FSJ\s+([\d.]+)`（:35，不分大小寫）原樣讀 GTF 的 fsj 值；沒有正規化、取整、座標合併。
- `build_matrix`（:62-95）：以 `pivot_table(aggfunc="sum", fill_value=0)`（:87-90）轉成 circRNA × sample 矩陣。**該樣本的 GTF 沒有此 circ（CIRIquant 認定 bsj=0 且 fsj=0 而略過）→ BSJ 與 FSJ 都被補 0**（所以矩陣中 FSJ=0 的格子有兩種來源：真的 0，或根本沒被該樣本偵測到）。
- 矩陣僅保留 consensus 的 high-confidence BED 內的 circRNA（:79-84）。

### 2.2 `analysis.R`：FSJ → offset（逐步，`scripts/analysis.R`）
1. :192-198 讀 FSJ 矩陣，`round()` 取整，與 BSJ 矩陣取共同 circRNA。
2. :200-203 `filterByExpr(bsj, …)`——**只用 BSJ 決定保留哪些 circ**；FSJ 矩陣只是跟著取同一子集。
3. :207 `fsj_pseudo <- fsj + 0.5`（**fsj=0 的處理：加 0.5 prior count，不排除**）。
4. :208-209 `calcNormFactors(DGEList(counts = fsj), method="TMM")`——TMM 是對**未加偽計數的原始 fsj**、**僅在已通過 BSJ filterByExpr 的子集**上算。
5. :210 `fsj_scaled = fsj_pseudo × norm.factor`（每個樣本一個因子）。
6. :211 `fsj_cpm = fsj_scaled / colSums(fsj_scaled) × 1e6`——分母是**該子集內**的 FSJ 總和（不是全庫存）。
7. :212 `offset_mat <- log(fsj_cpm + 0.5)`——自然對數，**再加一次 0.5**。
8. :215-216 `dge$offset <- offset_mat`：**offset 是 per-locus × per-sample 的矩陣**（不是 per-sample 純量）。edgeR 的 `getOffset`（我在 server 上印出函式本體確認）在 `y$offset` 非空時直接回傳它，**不再加 log(library size)**——也就是 BSJ 的 library size 完全沒有進模型，模型比較的是「BSJ count ÷ 該 locus 該樣本的 FSJ-CPM」。
9. :217-219 `estimateDisp` → `glmQLFit(robust=TRUE)` → `glmQLFTest(coef = cond_coef)`。
10. **獨立 FSJ 檢定另外處理**（:223-228）：`fsj_safe <- fsj + 1L`（加整數 1，非 0.5）、獨立的 `DGEList`、`calcNormFactors` 與一般 library-size 正規化、`estimateDisp`、`glmQLFit`、`glmQLFTest`——**與 offset 用的不是同一份轉換**（不同偽計數、不含 CPM/offset）。其 `logFC_fsj`、`PValue_fsj` 供 Type I/II 判斷（`analysis.R:168-170`）。

重現驗證：我把 :192-212 與 :215-219 原樣在 R 中重算（GSE133998，唯讀、輸出只寫 `~/audit_tmp/`），650 個 circRNA 的 log2FC 與 p 值和 `de_results.tsv` 的最大差異分別是 4.9e-15 與 1.8e-15——上述鏈與實際產出檔完全一致。

## 3. 真實數值範例 → `audit/fsj_worked_example.tsv`
GSE133998（6T+6N，unpaired 設計）。TSV 含三個 Table 2 前段的 circRNA × 12 樣本：BSJ 原始計數、FSJ 原始計數、TMM 因子、偽計數後 FSJ、縮放後、FSJ-CPM、offset（ln）、檔案中的 log2FC/p、重擬合的 log2FC/p。
**兩個邊界各自的接點計數無法提供**：CIRIquant 的 GTF 只輸出合併後的 fsj；中間的 `circ/`、`align/`、`gene/` 目錄被 `circrna.smk:60` 在每個樣本跑完時刪除，HISAT2 BAM 也已清理，需重跑 CIRIquant（約 11–12 h/樣本）才能取得。

**範例 A：N4BP2L2，chr13:33091994|33101669（Table 2 第 1 名，log2FC +6.65，p=0.0080，Type I）**
TMM 因子（12 樣本，順序 329…340）：1.3211, 0.6810, 0.9534, 1.1871, 1.2932, 0.8862, 1.4020, 0.7976, 1.0301, 1.0201, 0.9221, 0.7909。

| 樣本 | 組別 | BSJ | FSJ | FSJ-CPM（含 0.5 偽計數 + TMM） | offset = ln(CPM+0.5) |
|---|---|---|---|---|---|
| SRR11600329 | normal | 100 | 235 | 2331.5 | 7.755 |
| SRR11600330 | tumor | 60 | 278 | 1315.6 | 7.183 |
| SRR11600331 | normal | 181 | 325 | 2170.4 | 7.683 |
| SRR11600332 | tumor | 58 | 221 | 1953.0 | 7.577 |
| SRR11600333 | normal | 73 | 171 | 1856.5 | 7.527 |
| SRR11600334 | tumor | 32 | 257 | 1532.0 | 7.335 |
| SRR11600335 | normal | 128 | 351 | 4021.3 | 8.300 |
| SRR11600336 | tumor | 77 | 333 | 2051.3 | 7.627 |
| SRR11600337 | normal | 70 | 275 | 2502.2 | 7.825 |
| **SRR11600338** | **tumor** | **126** | **0** | **4.0** | **1.505** |
| SRR11600339 | normal | 84 | 261 | 1823.9 | 7.509 |
| SRR11600340 | tumor | 76 | 519 | 2489.4 | 7.820 |

（CIRIquant 原始 GTF：`GSE133998_results/circRNA/SRR11600338/SRR11600338.gtf` 該列為 `bsj 126.000; fsj 0.000; junc_ratio 1.000`，所以 0 來自 CIRIquant，不是 merge_counts。）

**這個範例同時是一個嚴重問題**：除 SRR11600338 外，腫瘤的 BSJ（60, 58, 32, 77, 76）反而低於正常（100, 181, 73, 128, 70, 84），FSJ 相近；只有那一個 FSJ=0 的腫瘤樣本把 offset 壓到 1.5（其餘約 7.5，差 6 個 ln 單位 ≈ 8.7 個 log2 單位）、使其 BSJ/FSJ 比值極大，單一樣本把腫瘤組的 log2FC 拉到 +6.65。我另外重擬合（**這是我的重算，不是管線輸出**）：
| circRNA | 全部 12 樣本 | 排除 FSJ≈0 的樣本後 | 純 BSJ（CPM）粗略 log2 比 |
|---|---|---|---|
| N4BP2L2（rank 1） | log2FC +6.65，p=0.0080 | 排除 SRR11600338：−0.46，p=0.31 | +0.04 |
| PTK2 chr8:141856359\|141900868（rank 2） | −7.25，p=0.0038 | 排除 SRR11600331、SRR11600337（正常組 FSJ=0）：−0.28，p=0.66 | +0.01 |
| ATP5C1 chr10:7839010\|7844817（rank 3） | −5.17，p=0.0043 | 排除 SRR11600335（正常組 FSJ=2）：−1.85，p=0.025 | +0.17 |
→ Table 2 前三名中，前兩名的「顯著」完全由 1–2 個 FSJ≈0 的樣本造成，第三名縮小到約 1/3 幅度。**不建議用範例 A 當論文的示範範例**（它示範的是失效模式）；若要放範例，建議挑 BSJ、FSJ 皆無 0 的 locus（`fsj_worked_example.tsv` 目前無此類；可由我再挑）。

## 4. 邊界情況

### 4.1 intronic / intergenic / antisense 的 fsj
`circ.py:530-550` 的 FSJ 計算**完全是位置式**，不需要 host transcript（`circRNA_attr` 只用於標註 `circ_type`，:876-）。實測（所有樣本 GTF 列，`f7.py`）：
| 資料集 | circ_type | 列數 | FSJ=0 比例 | FSJ 中位數 |
|---|---|---|---|---|
| GSE133998 | exon | 114,971 | 3.6% | 101 |
| | intron | 18,258 | 7.3% | 10 |
| | intergenic | 2,427 | 24.2% | 3 |
| | antisense | 3,507 | 8.4% | 7 |
| GSE113230 | exon | 146,304 | 19.5% | 6 |
| | intron | 22,347 | 34.9% | 2 |
| | intergenic | 3,583 | 37.0% | 1 |
| GSE108735 | exon | 29,929 | 2.1% | 41 |
| | intron | 6,948 | 8.8% | 10 |
| | intergenic | 4,180 | 15.4% | 6 |
→ 這些 circ 的 FSJ 不是「host transcript 的線性丰度」，而是邊界處任何線性 reads（且不分股向）。論文把 Type II 解讀為 host-gene 轉錄變化，對 intronic/intergenic/antisense 的 circ 缺乏依據（GSE133998 的 intergenic 中位數 3）。

### 4.2 GSE108735 的 57/14,262：你的推測（FSJ 趨近於零）被資料否定
- `filterByExpr` 只看 BSJ（`analysis.R:200`）。若把同一函式套到 FSJ 矩陣，通過數是 **332（2.33%）**，不是更少；BSJ 通過 57（0.40%），兩者交集 52。
- FSJ 並不稀少：GSE108735 每樣本 FSJ 總和中位數 307,941，BSJ 僅 7,133（FSJ/BSJ = 43.2，在 19 個資料集中屬中間偏低：GSE58135 235、GSE148036 210、GSE221107 124 更高）；只有 FSJ 欄為 0 的格子 87.5%，與 BSJ 為 0 的格子 87.0% 幾乎相同——表示「0」主要是該樣本沒偵測到此 circ（merge_counts 補 0），而非 FSJ 耗竭。
- 通過率與 BSJ 的 library 大小強相關：19 個資料集，通過率 vs 每樣本 BSJ 總和中位數 Spearman ρ=0.86（p=2e-6）；vs FSJ 總和 ρ=0.73；vs FSJ/BSJ 比值 ρ=−0.63（FSJ 相對越多反而越低，與你的推測方向相反）。GSE108735 的 BSJ 總和（中位 7,133）使 `filterByExpr` 的 CPM 門檻（約 5 reads 對應 CPM）極高，7 vs 7 樣本需在 ≥7 個樣本達標，所以僅 57 個過關。
→ 結論：**否定**。§2.8/§3.4.2/§4.5/§4.7 把 57/14,262 歸因於「FSJ 訊號不足」不成立；成因是 BSJ 深度/樣本間 BSJ 總量（腫瘤 7.4k–15.4k vs 正常 1.9k–3.0k 的 3–5 倍差距）。

## 5. 本次查證新增的兩項發現（會影響論文）
1. **FSJ 零值與 offset 失真（影響範圍）**：19 個資料集的 4,357 個 DEC 中，3,383 個（77.6%）至少有一個樣本滿足「BSJ ≥ 5 且 FSJ ≤ 2」；有 427 個 DEC 的 log2FC 方向與純 BSJ（CPM）方向相反，其中 359 個同時有此類樣本（`f5.py`；GSE192410 2,346 個中 1,960 個、116 個方向相反）。這個比例遠高於可歸因於生物變異者。
2. **成因（程式）**：`circ.py:445` 的單值覆蓋：一條 read 只被記給一個 circ。在 GSE133998 的 (BSJ≥5, FSJ=0) 格子中，98.1% 在同一樣本 GTF 中有與它共享起點或終點的兄弟 circ（97.1% 的兄弟 FSJ>0），對照 FSJ>0 格子只有 61.9%；GSE113230 為 84.1%（FSJ>0 對照 65.7%）；GSE108735 為 41.5%（對照 26.9%）（`f6.py`）。PTK2 就是例子：SRR11600331 的 GTF 中共享終點 141900868 的 9 個 circ，fsj 為 130, 99, 0, 0, 100, **0（PTK2 本尊，bsj 31）**, 201, 0, 433。這與覆蓋機制一致，但屬關聯性證據；N4BP2L2 在 SRR11600338 沒有座標完全相同的兄弟 circ，其 fsj=0 的原因無法由現存檔案判定（BAM 已刪）。

**建議（供你決定，我沒有動任何東西）**：§2.8 應照實寫出 FSJ 的定義（邊界線性 reads 的 fragment 聯集、MAPQ>10、5 nt 錨、無 strand/剪接要求）、offset 為 `log(CPM(含 0.5 與 TMM) + 0.5)` 的 locus×sample 矩陣，以及 FSJ=0 的 locus 不排除；Type II 解讀與 Table 2/§3.2 的 N4BP2L2 要重新評估。
