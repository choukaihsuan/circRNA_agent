# 程式碼凍結紀錄（2026-10-02）

## `results-2026-07` tag — 產出全部19個已發表資料集結果的程式狀態

```
commit 49b122dfcbcb0e501faf033db3c014b9b2a0085d
tag    results-2026-07
```

**論文所有行號引用應以此tag為準**，不是後續commit。這是2026-10-02參數開放與計分修正工作開始前的最後一個commit。

## `2f27d7b` — 參數開放 + biomarker score 缺值修正（尚未用於任何已發表結果）

```
commit 2f27d7b739ce1ea4d05012c41ec09aa47bf26ba7
```

內容：
- `config.yaml` schema 擴充（`tools.ciriquant.*`、`tools.dcc.*`、`consensus.min_confidence`、`de.filter_expr_*`、`de.deseq2_size_factor_type`、`de.paired_design_override`、`de.biomarker.*`）
- 新增 `scripts/validate_config.py`
- `workflow/Snakefile`、`workflow/rules/circrna.smk`、`workflow/rules/de.smk`：接上新參數（全部預設值與修改前行為逐字元等價，見 commit message）
- `scripts/rank_biomarkers.py`：修正 interaction（miRNA/RBP）缺值處理——未被 `predict_interactions.py` 查詢過的候選，原本被當成0分並納入min-max正規化範圍，現在改為排除於平均之外（新增 `interaction_queried` 欄位）；新增六維分數權重設定（預設等權，與論文現況一致）；新增 `--fully-evaluated-only` 獨立輸出模式
- `scripts/generate_report.py`：Biomarker區塊新增interaction覆蓋率揭露（低於80%變色示警）

**⚠️ 影響幅度（已用GSE113230驗證，非估計）**：平均score變化0.097、83%候選排名變動超過5名、Top-10有4個換人（見本目錄`new_vs_old_top10_GSE113230.md`或對話紀錄）。**此修正尚未套用到任何一個已發表資料集的`biomarker_candidates.tsv`**；要不要、何時對19組重跑這一步，是獨立待決定事項。

## `e085bc0` — compute_cost.py 三個bug修正（**Figure 6(c) / 892 / 118.9 / 30.3% / 24.5% 這組數字的唯一可重現來源**）

```
commit e085bc0
```

論文引用的 Figure 6(c) compute cost 數字（Our pipeline Total_wall_min=891.8、CPU-h=118.9、相對 CirComPara2_4tools 節省 30.3%、相對 nf-core/circrna 節省 24.5%）**在此commit之前無法由repo重現**——`results-2026-07`（`49b122d`）當時的`compute_cost.py`含三個bug（CirComPara2_4tools的STAR mate1/mate2時間寫死、DE_min錯誤併入Total、nf-core/circrna漏算STAR時間），算出的是811.1/1280.8/1084.2這組舊數字。**此commit之後，`benchmark/Snakefile`+`benchmark/scripts/compute_cost.py`才是這組benchmark數值的權威來源**；完整修正前後逐欄對照見`audit/compute_cost_before_after.tsv`。

重跑指令與log/tsv mtime核對方式見CLAUDE.md「compute_cost.py 三處修正」段落。

## `5dfc55d` — `/login` CSRF豁免（獨立修正，與參數開放/compute_cost無關）

```
commit 5dfc55d
```

CLAUDE.md已知問題表記載的既有bug修正（CSRFProtect全域啟用導致登入前無法取得有效session、CSRF永遠驗證失敗的死鎖）。調查確認此檔案改動僅此一行，與本次凍結範圍內其他工作無關聯，獨立commit。

## `798f696` — cross_dataset.html circ_id顯示0-based化（獨立修正，純前端顯示層）

```
commit 798f696
```

調查確認此檔案改動僅新增`_dispCoord()`做顯示層轉換（1-based內部座標→0-based顯示，比照circBase慣例），不影響任何查詢/join用的circ_id實際值，與本次凍結範圍內其他工作無關聯，獨立commit。

## `a8dfd79` — Snakefile_multi WIP（ENA FTP路徑泛化，刻意不完成）

```
commit a8dfd79
```

`_ena_fastq_base()` helper，屬於另一個已部分提交（`49b122d`）的「3-dataset RNase R benchmark expansion」功能，尚未接上、未測試，commit訊息明確標註「不可執行」。刻意以WIP狀態提交而非繼續放在工作目錄——未提交檔案留在工作目錄正是這次compute_cost provenance問題的根源之一，而這個檔案是獨立的Snakefile、不被主流程引用，提交沒有風險。完成與提交時機是獨立待決定事項，排在論文定稿與其他修訂之後。

## `c7e0559` — CLAUDE.md資料集事實校正 + Table2/compute_cost稽核紀錄

```
commit c7e0559
```

單一commit，涵蓋兩類先前未提交的文件更新（commit訊息已逐項列出）：
1. **資料集事實校正**（依`audit/per_dataset_summary.tsv`、`audit/thesis_dataset_table_19ds.tsv`，非沿用文件舊敘述）：GSE229705讀長100bp→51bp；GSE58135補揭露poly(A)-selected建庫法；新增19組資料集讀長/library prep/strandedness線索/有效單工具比例/15組跨資料集分析納入與否的稽核校正表。
2. **Table 2可重現性 + compute_cost修正紀錄**（本次session較早的調查成果，先前因「橫跨多主題」被擱置未提交）：版本號更正（CIRIquant/R/MultiQC）、CirComPara2/nf-core未實際執行的揭露、基因體來源查證、資料庫access date、compute_cost三個bug的修正說明與前後數字。

## 本次凍結（2026-10-02～03）最終狀態：無遺留未提交檔案

原列「來源不明」的`scripts/web_ui.py`、`scripts/templates/cross_dataset.html`已查明並分別提交於`5dfc55d`／`798f696`；原「多主題擱置」的`CLAUDE.md`已提交於`c7e0559`；原「無關但未完成」的`benchmark/Snakefile_multi`已以WIP狀態提交於`a8dfd79`。本次對話開始時的6個未提交檔案，現皆已入庫。

**⚠️ 798f696（cross_dataset.html的circ_id顯示改0-based）與2f27d7b內混入的`_disp_coord()`（generate_report.py多處0-based顯示）仍待覆核**：調查確認論文Table 2/Figure 9不印`chr:start|end`格式座標（只用circBase ID或基因名），故不構成「論文vs網頁」的直接衝突；但這兩處改動本身互相之間、以及跟pipeline資料層（TSV、consensus_filter、rank_biomarkers，全部1-based）、跟18/19已發表report.html的主要顯示（表格/JSON，1-based）不一致。建議revert兩處的0-based轉換，回到顯示原始1-based circ_id，待使用者決定後執行。

## 歷史改寫紀錄（2026-10-03，推送前；本機 main 的 origin/main..HEAD 這段，共改寫兩次）

原因：舊 `5b38d5b` 誤把整個 `audit/` 帶進 commit，含論文稿與其抽取檔。已從歷史移出論文稿、其全文抽取、論文圖、逐條論文稽核文件及其產生工具，以及所有產生出來的圖（規則見 `.gitignore`；判準：拿掉論文後仍有意義者留，成品圖表只追蹤腳本）。未推送過。

**tag 與基準 commit 不受影響**：`results-2026-07`（annotated tag `88f54bd`）指向 `49b122d`，即 `origin/main`，在改寫範圍之外，未變動。`e085bc0`（compute_cost 修正）在範圍內但早於第一個被改動的 commit，hash 不變。

舊（改寫前原 hash）→ 新（目前）對照：

| 舊 hash | 新 hash | 說明 |
|---|---|---|
| 7154a69 | decfd7b | 跨資料集數字對照表與腳本（圖已移出追蹤） |
| b403018 | ac74e99 | organ/cancer-type 說明文字修正 |
| 5b38d5b | a438230 | audit 分析產物；訊息已改為實際內容 |
| 116fc0b | 28d8742 | 論文風格 Fig 7/8/9 腳本（成品圖不追蹤） |
| a8dfd79 | f71b361 | Snakefile_multi WIP（`_ena_fastq_base()` 未完成，不可執行） |
| c7e0559 | 41a8c59 | CLAUDE.md 資料集事實校正 |
| afa0e27 | d9852cf | 本檔先前版本 |
| 49b122d、2f27d7b、eda5424、e085bc0、5dfc55d、798f696 | 不變 | |

本檔以上文字中的 `5b38d5b`、`a8dfd79`、`c7e0559`、`116fc0b`、`afa0e27` 皆指改寫前的 hash，請依此表對照。
