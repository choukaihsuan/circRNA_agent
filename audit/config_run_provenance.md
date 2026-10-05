# 執行所用 config 的來源查證（15 組論文資料集，唯讀）

問題：已發表的執行是否都帶了 `--configfile config/projects/{GSE}.yaml`？還是有的用全域 `config.yaml`？

## 證據與限制

1. **Snakemake 的 config 訊息**：`Config file config.yaml is extended by additional config specified via the
   command line` 只在有傳入 `--configfile` 時出現（對照測試：無參數 0 次，帶專案快照 1 次）。但這則訊息不會說明
   是哪一個檔案，所以它只能證明「有傳入額外 config」，不能證明是專案快照。
2. **CIRIquant metadata 的 `-t`**：`.snakemake/metadata` 的 shellcmd 記錄了產生 GTF 的那次 CIRIquant 的實際
   thread 數。全域值有兩個：server 上的副本是 18，已 commit 的版本是 8。若某組的快照值既不是 18 也不是 8，
   而 `-t` 等於它，就能確定當時讀了專案快照。
3. **bash 歷史**：只有 13 行 `--configfile config.yaml`，沒有記錄專案快照的執行，所以無法用它判定。

## 結果（15 組）

| 資料集 | 快照 threads | 已發表 GTF 的 CIRIquant `-t` | 判定 |
|---|---|---|---|
| GSE121842 | 24 | 24 | 專案快照（確定） |
| GSE130078 | 12 | 12 | 專案快照（確定） |
| GSE133998 | 12 | 12 | 專案快照（確定） |
| GSE143797 | 24 | 24 | 專案快照（確定） |
| GSE148036 | 24 | 24 | 專案快照（確定） |
| GSE221107 | 12 | 12 | 專案快照（確定） |
| GSE248612 | 12 | 12 | 專案快照（確定） |
| GSE77509 | 12 | 12 | 專案快照（確定） |
| PRJNA553289 | 12 | 12 | 專案快照（確定） |
| GSE108735 | 18 | 18 | 無法區分（與全域 18 相同） |
| GSE136569 | 18 | 18 | 無法區分（與全域 18 相同） |
| GSE171011 | 18 | 18 | 無法區分（與全域 18 相同） |
| GSE97239 | 18 | 18 | 無法區分（與全域 18 相同） |
| SRP156355 | 8 | 8 | 無法區分（與已 commit 全域 8 相同） |
| **GSE113230** | **12** | **8** | **不一致**：GTF 由 threads 8 產生，快照的 12 不是產生該 GTF 的值 |

### GSE113230 的不一致

- 快照寫 12，但 6 個樣本的 CIRIquant GTF 都是以 `-t 8` 產生。
- CLAUDE.md 記錄 GSE113230 的 GTF 沿用舊檔（`-Nr` 重跑時未重跑 CIRIquant）。舊 GTF 產生於快照設為 12 之前，
  所以快照的 threads 與論文使用的 GTF 不一致。
- 附錄 A2 的 threads 欄以 metadata 的實際 `-t` 為準（GSE113230 = 8），並在表內註明。

### 沒有 pipeline log 或 log 不含 configfile 訊息的組別

- GSE113230：`logs/` 中沒有名為 `pipeline_GSE113230*` 的 log（只有 `reNr2_GSE113230_*` 階段 log）。
- GSE143797：最近的 `pipeline_GSE143797-rerun.log` 不含 configfile 訊息（較早的 `pipeline_GSE143797.log` 有）。
  最近一次 rerun 很可能只重跑了下游，且 CIRIquant GTF 仍來自先前的執行。

## 結論

- 9 組：有直接證據，已發表的 CIRIquant GTF 是以專案快照的 threads 產生。
- 5 組：無法從 threads 區分，但效果上與全域值相同，沒有觀察到不一致。
- 1 組（GSE113230）：不一致，已如實記錄。
- 「15 組全部都帶了專案快照」這個說法 **無法被完整證明**，只能陳述上述範圍。

## 與其他紀錄的關係

- `config.yaml` 在 server 上未與 repo 同步（threads 18 對 8，且缺 `tools.dcc` 區段）。這是 repo 衛生問題，
  不影響已發表執行的 threads 來源（上表證據已排除全域值的影響，除了 5 組無法區分者）。
- 完整的 provenance 需要 B5 的 run_params.yaml；目前只能依上述間接證據。
