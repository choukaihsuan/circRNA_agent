# 第 8 項：已發表 v2 候選檔的終端重新產生（最終範圍與結果）

## 範圍（修訂）

原計畫是全量 Snakemake 重新產生。A 節的 dry-run 證明這不可行：17 組資料集的 BAM 與 trimmed FASTQ 已在磁碟清理時刪除，Snakemake 會規劃重建 STAR 與 CIRIquant，而 CIRIquant 約需 11–12 小時／樣本。B 節的預測顯示五項變更中只有 92a6713 tie-break 會改變輸出，且只改同分列的 rank。因此改為終端重新產生，不使用 Snakemake：

1. `rank_biomarkers_v2`（三個方法）— 直接呼叫 `scripts/rank_biomarkers.py`，參數與 `workflow/rules/de.smk` 相同，LOO 目錄、cross-method 與 stability 門檻取自各資料集 status.json（三個方法一致，已核對）。加上 `--no-run-params`，不產生 run_params.yaml（留待 B5）。
2. `generate_report` — 以 mock 物件呼叫 `scripts/generate_report.py`，輸入與參數同 `de.smk`。

工具：`audit/regen_terminal.py`、`audit/verify_regen_candidates.py`、事後比對用的 manifest（sha256，未 commit）。

## 前置修正

- **sig_cap / fc_cap / stability_max_candidates 的 v2 接線**（commit 4c5d41a）：rank_biomarkers_v2 原本不輸出這三個 flag，config 鍵對已發表的候選檔無效。現在偏離預設值才輸出（與 v1 相同寫法），預設命令不變。
- **dry-run 斷言測試** `audit/test_biomarker_wiring.sh`（同 commit）：7/7 通過。每個偏離預設值的 `de.biomarker` 鍵都出現在命令列中；控制組不含任何可選 flag。

## 執行前後的證據

- 執行前：15 組的完整檔案 sha256 manifest，與封存時（`~/archive_published_20261005/`）逐檔一致（15/15 SAME）。
- 執行：rank 45 次（15 組 × 3 方法）、generate_report 15 次，全部 exit 0。
- 執行後：完整 manifest 比對，變動檔案只有三類：
  - `de/biomarker_candidates_v2_{edgeR,DESeq2,limma}.tsv`（edgeR 7 個檔、DESeq2 10 個、limma 15 個，合計 32 個）
  - `report.html`（15 個）
  - 其他所有檔案逐位元不變，包括 `de_results.tsv`、`count_matrix.tsv`、`high_confidence.bed`、status.json、plots，且沒有產生 run_params.yaml。
- 欄位級驗證（`verify_regen_candidates.py`，45/45 檔案）：
  - 所有欄位（除 `rank`）逐值相同；`rank_default`、`pass_default`、列數、circ_id 集合全部相同。
  - `rank` 變動共 15,991 列，**全部**落在完全同分（biomarker_score 精確相等）的群組內；無法解釋的變動為 0。
  - 與 scratch 實驗及 B 節的檔案推算一致。

## Table 2 核對（GSE133998，edgeR，重新產生後）

| Gene | circ_id | score（檔案） | 論文 | rank_default |
|---|---|---|---|---|
| DENND4C | chr9:19286767\|19305525 | 0.5665 | 0.567 | 1 |
| SCAF11 | chr12:46319925\|46322642 | 0.5222 | 0.522 | 2 |
| FAT1 | chr4:187627717\|187630999 | 0.5192 | 0.519 | 3 |
| RNU7-61P | chr16:30495148\|30495584 | 0.5051 | 0.505 | 4 |
| CHD9 | chr16:53288350\|53308214 | 0.3976 | 0.398 | 5 |
| UBXN7 | chr3:196118684\|196129890 | 0.3272 | 0.327 | 6 |
| MAPK8 | chr10:49609655\|49618211 | 0.3210 | 0.321 | 7 |
| RP11-74K11.2 | chr12:94562929\|94580249 | 0.2062 | 0.206 | 8 |

- 候選總數 79、通過 8，與論文一致。八列的 rank（20, 29, 30, 34, 52, 66, 68, 75）與重新產生前相同。
- 變動的兩列（MAP3K4、AFF4，同分 0.385）都未通過，不在 Table 2 中。
- 「Rank of 79」以 `rank` 欄（79 個候選中的名次）解讀，這個解讀請你確認。

## 其他數值的影響

- 候選檔的列數、`pass_default`、`rank_default`、分數全部不變。
- 論文中其餘數值（Table 4 的 tested／significant、§3.3 的 benchmark、§4.2 權重敏感度、§4.10 的 nested 統計）來自未被重新產生的檔案或獨立目錄，預期零變動；本輪以 sha256 證明它們的輸入檔逐位元不變。

## 四組排除資料集（GSE58135、GSE192410、GSE192849、GSE229705）

不重新產生，停留在舊 commit，明確註明：

- GSE229705 的 dry-run 因 `IncompleteFilesException` 無法規劃，未執行 `--cleanup-metadata`。它不供應任何論文數字，直接排除。
- **GSE192410 與 GSE192849 是 `--dcc-count-threshold` 修正唯一會咬到的資料集**。這兩組的 9 個 DCC 樣本是 header-only（0 條資料列），fallback 從未產生任何 locus。因此這個修正在資料上從未生效，它們維持舊版是可接受的，但必須寫在 audit 中，不能默默留著。

## 其他查證結果

- **CIRIquant 的 provenance 觸發**（`audit/ciriquant_provenance_check.md`）：來源是 2f27d7b（2026-10-02）對 rule 的 params 變更（新增 anchor_arg、no_gene）。所有資料集的 ciriquant metadata 都只記錄 2 個 params，早於該 commit。不是無法解釋的 hash 變動；預設命令逐位元不變，不影響輸出。
- **GSE133998 的 CIRIquant GTF 寫入保護**：保留，不解除（CLAUDE.md 已知問題表的規程）。

## 未做的事（依指示）

- 不重跑 CIRIquant、STAR、DCC、consensus、merge_counts、de_analysis。
- 不解除 GSE133998 GTF 的寫入保護。
- 不對 GSE229705 執行 `--cleanup-metadata`。
- 不修改任何 `.docx`。
