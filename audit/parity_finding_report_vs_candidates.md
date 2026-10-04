# 稽核發現：已發表報告的生物標記分數 vs 候選檔（2026-10-04）

**性質**：記錄既成事實。兩條計分路徑在統一之前已經分歧。已發表的報告與論文數字**未回頭修改**；
報告改為讀取候選檔後，15 份報告會以候選檔數值重新產生（封存於 `published_reports_pre_v2/`）。

## 兩條計分路徑

1. **候選檔**（`scripts/rank_biomarkers.py`）：`biomarker_candidates.tsv` 與 v2 的
   `biomarker_candidates_v2_<method>.tsv`。
2. **報告內嵌**（`scripts/generate_report.py` 的 `_compute_bm_table_data` 與 `_bm_lookup`）：
   依各方法的 DE 表自行計分，並在報告端重算 `mirna_n`／`rbp_n`，再以 `interactions.json` 覆寫 `n_mirna`。

## 比對結果（15 組，同一 circRNA 的 biomarker_score）

| 方法 | 比對對象 | 結果 |
|---|---|---|
| edgeR | 舊 v1 候選檔 | 15 組全部一致（≤ 0.0005，即小數點第三位的捨入） |
| edgeR | 現行 v2 候選檔 | 9 組有差異，最大差值 0.30（GSE113230）；6 組一致 |
| DESeq2 | v2 候選檔 | 15 組全部有差異，最大差值約 0.23；miRNA 計數差值最高 2,968 |
| limma | v2 候選檔 | 15 組全部有差異，最大差值約 0.23；miRNA 計數差值最高 2,968 |

edgeR 9 組的差異來源是 2f27d7b 之前的零填補（未查詢 miRNA/RBP 的 circRNA 原本以 0 計入 6 維平均，
之後改為以其餘 4 維平均）。報告內嵌的 edgeR 分數與 v1 檔一致，即報告沿用了舊的零填補。

DESeq2 與 limma 的差異來源未在本次追溯到單一程式行；報告內嵌路徑對這兩個方法的 `n_mirna` 來源與候選檔不同，
這是兩條路徑分歧的直接證據。

## 處置

- 已發表數字不回改。封存既有 15 份報告後，以候選檔重新產生。
- 論文方法段需寫明：分數正規化母體為 edgeR 顯著集合；三方法一致與效應量穩定性於排序後套用。
- 此分歧列入 §4.10 的限制或稽核附錄，依作者決定。
