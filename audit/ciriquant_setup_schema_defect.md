# 缺陷紀錄：generate_ciriquant_config() 產生 CIRIquant 無法讀取的 schema（代辦第 13 項，先記錄不修）

## 現象

`scripts/utils.py` 的 `generate_ciriquant_config()` 寫出的 YAML 以 `database:` 為頂層鍵，其中索引欄位是
`database.hisat2_genome`。CIRIquant 1.1.3 讀的是 `reference.hisat_index`（見
`CIRIquant/utils.py:108` 的 `if 'hisat_index' not in config['reference']`）。用這個函式產生的檔案，CIRIquant
會因為缺少 `reference` 區段而拒絕執行（CLAUDE.md 已知問題表有記：`KeyError: 'reference'`）。

## 觸發路徑

- `scripts/agent.py --setup-ciriquant`（`cmd_setup_ciriquant`，`agent.py:93`）呼叫此函式。
- `scripts/agent.py` 的用法說明列出該子命令（第 8 行）。README.md 中我沒有找到這個子命令，所以它不在使用者文件裡；
  只在程式的 `--help` 中可見。

## 影響

- 產生的檔案無法直接使用。
- 目前本研究的所有執行都使用手寫的 `config/ciriquant.yaml`（正確 schema，絕對路徑，指定 java 與 perl），
  所以已發表結果不受影響。
- 但「他人可重現」的主張依賴這個入口。任何人照程式的 `--help` 跑 `--setup-ciriquant` 都會得到壞檔案。

## 與 genome.hisat2_index 的關係

`genome.hisat2_index`（專案快照）只被此函式消費，寫入的正是上面那個錯誤的鍵。因此快照中的這個值從未決定
CIRIquant 使用的索引。實際的索引位置是共用檔 `config/ciriquant.yaml` 的 `reference.hisat_index`，
目前指向 `/home3/choukaihsuan/reference/hg19/hg19_hisat2_index`（8 個 `.ht2` 索引檔存在）。

## 修正方向（未執行）

- 讓 `generate_ciriquant_config()` 產生 `reference:` 區段（`fasta`、`gtf`、`bwa_index`、`hisat_index`），
  並加入 java 與 perl 等工具路徑；或改為不產生、只驗證手寫檔。
- 加入一個測試：產生的檔案必須能被 CIRIquant 的 `check_config` 通過。
- 在 README 補上這個子命令的正確使用方式，或移除它。

## 與 CIRIquant provenance 的同一缺口

`config/ciriquant.yaml` 是單一共用檔，手動維護，沒有逐次執行紀錄（GSE97239 在 2026-07-29 有修改過）。
這與 `audit/ciriquant_provenance_check.md` 記錄的 provenance 缺口是同一件事，不另外列為第二條。
