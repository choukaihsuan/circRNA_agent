# README 實測報告（clean 環境照 README 走一次）

日期：2026-10-09。環境：雲端 session（Linux、Python 3.13、無 conda／R／任何生資工具、有 docker CLI 但**無 daemon**、
無 Singularity、可連 PyPI／conda channel／Docker Hub／GitHub，**連不到 NCBI**，也連不到 `172.16.0.178`）。
原則：只做唯讀或在 scratchpad 內可丟棄的操作；沒有建立 conda 環境（專案規則：安裝／建環境前須先回報）。

## 結果總表

| # | README 步驟 | 結果 | 說明 |
|---|-------------|------|------|
| 1 | `git clone` | ✅ | |
| 2 | Requirements 表 | ⚠️ | 原文寫「Python ≥ 3.7」，但 `envs/circrna.yaml` 解出 **Python 3.9、R 4.3.1、edgeR 4.0.16、DESeq2 1.42**；已發表結果是在 server（Python 3.7、R 4.2.2、**edgeR 3.40**）產生。新環境的 DE 數字可能略有差異 → 已寫入 README「Known limitations」 |
| 3 | `mamba env create -f envs/circrna.yaml` | ✅（只做 solve） | 以 micromamba `create --dry-run` 驗證：431 個套件、下載 913 MB、solve 約 68 秒，CIRIquant 1.1.3、DCC 0.5.0、STAR 2.7.10b、HISAT2 2.2.1、BWA 0.7.17、samtools 1.17、snakemake 7.32.4 皆可解。**未實際安裝**。環境名稱確為 `ciriquant`（與 README 一致）|
| 4 | 編輯 `config.yaml`／`config/ciriquant.yaml` | ❌ 原本走不通 → ✅ 已補模板（2026-10-09）：新增 `config.example.yaml`、`config/ciriquant.example.yaml`（`/path/to/...` 佔位、有測試確保鍵與 live 檔一致且不含個人路徑），README 改為 `cp` 模板 | 兩檔是維護者的工作副本：`/mnt/d/ref/hg19/...`、`/home/choukaihsuan/miniforge3/envs/circrna/bin/...`、`download.sra_cache_dir: /home/choukaihsuan/sra_cache`。沒有 `config.example.yaml`／`ciriquant.example.yaml` 模板，審閱者必須自己逐欄改。原因（保留供查）|
| 5 | `python scripts/web_ui.py --host 0.0.0.0 --port 5000` | ✅ 可啟動 | 在乾淨 venv（flask、flask-wtf、pyyaml、pandas）下啟動成功 |
| 6 | 登入（README：「沒設 email 時 console 會印出 magic link」）| ❌ → ✅ 已修 | **原本走不通**：程式的 console fallback 只印收件者與主旨，**不印連結**，所以沒有郵件服務就無法登入 Web UI。已新增明確 opt-in 的 `PIPELINE_DEV_PRINT_LINK=1`（預設關閉；公開部署不可開）並更新 README。端對端實測：`POST /login` → console 取得連結 → `GET /auth/<token>` 302 → 首頁 200；未登入存取首頁 302 到登入頁 |
| 7 | Method 1：GEO 一鍵（`agent.py --gse`）| ❌ 雲端無法驗證 | 需要 `pysradb`（`envs/circrna.yaml` 的 pip 區段有；本 venv 沒有 → `ImportError: pysradb not installed`，訊息清楚）與 NCBI 連線（本環境連不到）|
| 8 | Method 2／3（手動 SRR、本地 FASTQ）| ✅ 輸入驗證與佇列行為 | 以 Flask test client 驗證（見 `tests/test_web_security.py`）：合法輸入入佇列、非法輸入 400。**未執行 pipeline** |
| 9 | 實際跑一次完整分析 | ❌ 雲端不可行 | 需要約 0.9 GB 的環境＋參考基因體＋索引＋數小時至數天運算＋數百 GB 磁碟 |
| 10 | CLI：`snakemake --snakefile workflow/Snakefile --configfile config.yaml ...` | ⚠️ 僅 dry-run | 在暫存 fixture 內 `snakemake -n` 成功建 DAG（`tests/test_config_wiring.py`）。`agent.py --gse ... --dry-run` 同樣需要 pysradb／NCBI |
| 11 | Docker | ⚠️ 部分 | Docker Hub 上 `choukaihsuan/circrna-pipeline` 的 **1.0.1 與 1.0.0 都存在**（各 9 層、約 1.7 GB 壓縮）。`docker pull/run` **未執行**（無 daemon）。原 README 只有 Singularity 指令、**沒有純 Docker 用法**；而 Dockerfile 的映像**不含本 repo 程式碼**（只裝工具鏈），直接 `docker run` 跑不起來 → 已補 README「Plain Docker」並說明需掛載 repo |
| 12 | Singularity `--use-singularity` | ✅ 不一致已修（Snakefile 改 `:1.0.1`）；實跑未測 | `workflow/Snakefile` 的 `singularity:` 寫 `:1.0.0`，但 README／Dockerfile／CLAUDE.md 都寫 `1.0.1`。照 README 用 `--use-singularity` 會拉到**舊版映像**。已依你的決定改為 `1.0.1`。實際 `--use-singularity` 執行仍未測（無 Singularity）|
| 13 | 最小可重現範例 | 原本 ❌ 沒有 → ✅ 已新增 | README 的 Quick Start 直接要求完整 server＋基因體。新增 `examples/minimal/`（僅需 Python 3、<1 秒）：雙工具 consensus 步驟＋預期輸出＋CI 測試 |
| 14 | 輸入與預期輸出 | 原本 ⚠️ 零散 → ✅ 已補 | 新增「Inputs and expected outputs」（metadata 格式、輸出檔樹）|
| 15 | 已知限制 | 原本 ⚠️ 散落各處 → ✅ 已彙整 | 新增「Known limitations」|
| 16 | Configuration Reference | ⚠️ 已修 | README 寫 `de_sig_by: pvalue`，`config.yaml`／rules 實際預設是 `auto`（Storey q<0.2 → 否則 nominal p<0.05），`rank_biomarkers.py` 還接受 `qvalue`。已改。**CLAUDE.md 也有同樣的過時敘述（範例 config 寫 `pvalue`）**，未動 |

## 仍未解決（需要你決定）

1. （已解決）config 模板與映像版本。維護者的 `config.yaml` 仍在追蹤中，是否改名／移出追蹤由你決定。
2. **版本漂移的誠實揭露**（#2）：README 現在寫明「新環境可能與已發表數字略有差異」。若論文 §2.14 宣稱可重現，建議審閱者拿到的是 pin 住 R／edgeR 版本的環境檔（目前 `bioconductor-edger>=3.40` 不 pin，解出 4.x）。
3. **無法在雲端驗證的項目**：實際 pipeline 執行、Docker／Singularity 實跑、GEO 一鍵啟動。這些需要在 server 或有 daemon 的機器上走一次。
