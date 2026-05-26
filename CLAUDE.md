# circRNA Analysis Pipeline — Project Context for Claude

## 專案概述

本專案是一個以 **Snakemake** 驅動的 circRNA（環狀 RNA）全流程分析管線，
從 GEO/SRA 原始數據下載，到差異表現分析（DE）與 HTML 報告輸出。

- **目標數據集**：GSE113230（肝癌 tumor vs. normal，6 個 sample）
- **主要工具**：CIRIquant（circRNA 偵測）+ DCC（輔助偵測，雙工具共識）
- **執行環境**：基因體中心 HPC server（`172.16.0.178`，CentOS 7，96 cores，377 GB RAM）
- **本機開發**：Windows 11 + WSL2（Ubuntu），程式碼在 `/mnt/c/Users/User/develop/circRNA_agent/`
- **Server 路徑**：`~/circRNA_agent/`（即 `/home3/choukaihsuan/circRNA_agent/`，`/home/choukaihsuan` 是 symlink）

---

## 目錄結構

```
circRNA_agent/
├── config.yaml                  # 主設定檔（路徑、參數、工具選擇）
├── config/
│   ├── ciriquant.yaml           # CIRIquant 工具路徑設定（server 版本另存在 server 上）
│   └── .ciriquant_ready         # touch 檔，驗證 ciriquant.yaml 存在後建立
├── metadata/
│   ├── library_info.csv         # SRR ID、配對資訊（srr_id 欄位為必要）
│   └── sample_groups.csv        # 樣本分組（srr_id, condition 欄位）tumor/normal
├── workflow/
│   ├── Snakefile                # 主 Snakefile，載入 rules，設定 target
│   └── rules/
│       ├── download.smk         # SRA 下載（prefetch + fasterq-dump）
│       ├── qc.smk               # FastQC + fastp + MultiQC
│       ├── circrna.smk          # circRNA 偵測主規則（見下方詳述）
│       └── de.smk               # 差異表現 + 圖表 + HTML 報告
├── scripts/
│   ├── agent.py                 # CLI 入口（--gse, --setup-ciriquant 等）
│   ├── prepare_metadata.py      # 從 GEO RunInfo 建立 library_info.csv
│   ├── download_geo.py          # SRA 下載輔助
│   ├── consensus_filter.py      # CIRIquant + DCC 共識過濾（輸出 BED + confidence score）
│   ├── merge_counts.py          # 從 CIRIquant GTF 建立 BSJ + FSJ count matrix
│   ├── annotate_circbase.py     # circBase hg19 座標比對注釋（自動下載或讀本地檔）
│   ├── rank_biomarkers.py       # Biomarker 候選排序（composite score）
│   ├── assign_isoforms.py       # circRNA BSJ 座標對應 host gene（依 GTF）
│   ├── isoform_switching.R      # 計算 IUI、Wilcoxon rank-sum 測試 isoform switching
│   ├── analysis.R               # DE 分析（edgeR_ciriquant / deseq2 / limma）
│   ├── generate_report.py       # 輸出 HTML 報告（互動式 Plotly + Type I/II + biomarker）
│   ├── utils.py                 # 共用工具函數
│   ├── web_ui.py                # Flask Web UI（GEO 一鍵啟動 + 進度視覺化）
│   └── templates/
│       ├── index.html           # 主設定頁面（GEO 入口 + Step 1-3）
│       └── status.html          # Pipeline 狀態頁（進度條 + rule 狀態格 + log）
├── envs/
│   └── circrna.yaml             # Conda 環境定義
└── logs/                        # Snakemake 各 rule 的 log 檔
```

---

## Pipeline 流程

```
SRA/GEO
  │
  ▼
[download] prefetch + fasterq-dump
  → raw_dir/{srr}_1.fastq.gz, {srr}_2.fastq.gz
  │
  ▼
[qc] FastQC（raw）+ fastp（trim）+ MultiQC
  → trimmed_dir/{srr}_1.fastq.gz, {srr}_2.fastq.gz
  → results/qc/fastp/{srr}.json, multiqc_report.html
  │
  ├──────────────────────────────────────────┐
  ▼                                          ▼
[ciriquant]                         [star_align]（paired-end）
HISAT2 + BWA-MEM                    STAR chimeric junctions
→ {srr}/{srr}.gtf                   → Chimeric.out.junction
  (BSJ + FSJ counts)                → Aligned.sortedByCoord.out.bam
                                           │
                                    ┌──────┴──────┐
                                    ▼             ▼
                             [star_align_mate1] [star_align_mate2]
                             STAR R1 only      STAR R2 only
                             → mate1/Chimeric  → mate2/Chimeric
                                    │
                                    ▼
                               [dcc]
                               DCC 0.5.0
                               -mt1 mate1/Chimeric
                               -mt2 mate2/Chimeric
                               → DCC/CircCoordinates
  │                                 │
  └──────────────────┬──────────────┘
                     ▼
             [consensus_filter]
             共識過濾（min_tools=2, slop=10, min_bsj=2）
             pseudo-circ QC（BSJ/FSJ > max_junction_ratio 過濾）
             confidence_score = Σ[log2(bsj+1)×(1−dist/slop)] / n_supporting
             → high_confidence.bed（7 欄，含 confidence_score）
             → consensus_summary.tsv
             （支援單工具模式：USE_CIRIQUANT / USE_DCC）
                     │
                     ▼
             [merge_counts]
             parse GTF → BSJ + FSJ count matrix
             → circRNA/count_matrix.tsv
             → circRNA/fsj_count_matrix.tsv
                     │
                     ├─────────────────────────┐
                     ▼                         ▼
             [de_analysis]              [annotate_circbase]
             edgeR_ciriquant /          比對 circBase hg19
             deseq2 / limma             → circRNA/circbase_annotated.tsv
             → de/de_results.tsv              │
             → plots/*.pdf                    │
                     │                        │
                     └──────────┬─────────────┘
                                ▼
                        [rank_biomarkers]
                        composite score：
                        significance + |log2FC|
                        + confidence + circBase bonus
                        → de/biomarker_candidates.tsv
                                │
                                ▼
                        [generate_report]
                        → report.html（含 Type I/II、biomarker 表）
```

---

## 關鍵技術細節

### circRNA 偵測工具選擇（config `consensus.tools`）

`Snakefile` 在啟動時讀取 `config.consensus.tools`（list），設定以下全域變數：

```python
TOOLS        = set(config["consensus"]["tools"])  # {"ciriquant", "dcc"} 或子集
USE_CIRIQUANT = "ciriquant" in TOOLS
USE_DCC       = "dcc"       in TOOLS
```

- **兩者並行**（預設）：`consensus_filter` 取交集，min_tools=2
- **CIRIquant only**：跳過 STAR 三次比對 + DCC；`consensus_filter` 單工具模式（min_tools=1）
- **DCC only**：跳過 CIRIquant；`consensus_filter` 單工具模式

`consensus_filter.py` 的 `--cirique` 和 `--dcc` 都是 **optional**，至少提供一個即可。

### DCC 的正確呼叫方式（DCC 0.5.0 重要）

DCC 0.5.0 必須同時提供三個 junction 檔案，缺少 `-mt2` 會出現：
`TypeError: object of type 'NoneType' has no len()`

```bash
DCC {paired_junction} \
    -mt1 {mate1_junction} \     # R1 單端 STAR 比對
    -mt2 {mate2_junction} \     # R2 單端 STAR 比對
    -D -an {gtf} \
    -Pi -F -M -Nr 5 1 \
    -G -B {bam} \
    -O {outdir} -T 4
```

- 命令是大寫 `DCC`（不是 `dcc`）
- BAM 參數是 `-B`（不是 `-A`）
- `-mt1`/`-mt2` 直接傳檔案路徑（不支援 `@filelist` 格式）
- **`-fg` 已移除**：觸發 CircSkip 解析，GTF exon_number 屬性尾部引號造成 `ValueError`

### Consensus 過濾與 Confidence Score

`consensus_filter.py` 的過濾流程（每個 sample 獨立執行）：

1. **min_bsj 過濾**：各工具輸出中，BSJ < min_bsj 的 circRNA 直接丟棄
2. **Pseudo-circ QC**（CIRIquant only）：解析 GTF 中的 FSJ，若 BSJ/FSJ > `max_junction_ratio`（預設 1.0）則丟棄。真實 circRNA 幾乎都有 BSJ < FSJ；比值 > 1 是誤比對或 repeat 假陽性的警訊。FSJ = 0 的 locus 不做此過濾（缺乏對應線性轉錄本資訊）
3. **座標共識投票**：各工具的座標在 slop 範圍內（預設 10 bp）視為一致
4. **Confidence score** = `Σ[log2(bsj+1) × (1−dist/slop)] / n_supporting_tools`
   - 分母是**實際支持該 circRNA 的工具數**（非所有工具數），代表「各支持工具的平均每工具信心」
   - 論文應標明此為 weighted scoring heuristic，**非機率值**，引用 CirComPara2 + Hansen (2018) 作為 consensus 正當性基礎

### circBase 注釋（`annotate_circbase.py`）

- 從 `http://www.circbase.org/download/hsa_hg19_circRNA.txt` 自動下載（或讀 `--circbase-file` 本地檔案，預設快取於 `/tmp/circbase_hg19.txt`）
- 以 max(|Δstart|, |Δend|) ≤ slop 判斷為已知 circRNA
- 輸出欄位：`circbase_id`（或 `"novel"`）、`circbase_gene`、`in_circbase`（0/1）

### Biomarker 排序（`rank_biomarkers.py`）

四維 composite score，每維度在 significant set 內 min-max 標準化後平均：

```
biomarker_score = (sig_norm + fc_norm + conf_norm + known_bonus) / 4
  sig_norm  = -log10(padj),  上限 10，標準化
  fc_norm   = |log2FC|,      上限 5，標準化
  conf_norm = confidence_score 標準化
  known_bonus = 1 若 in_circbase，否則 0（不做標準化）
```

輸入：`de_results.tsv` + `circbase_annotated.tsv`，只對 FDR 顯著的 circRNA 排序。

### CIRIquant GTF 輸出格式

```
chr1  CIRIquant  circRNA  100  200  .  +  .  gene_id "G1"; BSJ 5; FSJ 20; circ_id "chr1:100|200";
```

- `BSJ`：back-splice junction reads（circRNA 特異性）
- `FSJ`：forward-splice junction reads（linear RNA 同位點）
- 兩者都由 `merge_counts.py` 的 `parse_gtf()` 提取，輸出兩個獨立矩陣

### DE 分析方法（`config de.method`）

| 值 | 說明 |
|----|------|
| `edgeR_ciriquant`（預設）| 複製 `CIRI_DE_replicate`：edgeR GLM + FSJ offset（測 BSJ/FSJ 比值），輸出 Type I/II 分類 |
| `deseq2` | DESeq2 RLE normalization + Wald test |
| `limma` | limma-voom |

**`edgeR_ciriquant` 核心邏輯**：
FSJ counts → TMM normalization → `log(lib.size × norm.factors)` 作為 GLM offset
→ 效果等同於測試 BSJ/FSJ 比值是否在 tumor vs. normal 之間改變
→ 同時對 FSJ 跑獨立 QLFTest，若 FSJ 也顯著且方向相同 → Type II（基因層次），否則 → Type I（circRNA 專一性）

`analysis.R` 有 **backward-compatible fallback**：若 `snakemake@input[["fsj_matrix"]]` 不存在（舊 DAG），自動 fallback 到 deseq2。

---

## 設定檔說明（config.yaml）

```yaml
project_id: GSE113230
metadata:   metadata/library_info.csv
groups:     metadata/sample_groups.csv

# 路徑（server 版本使用 /home3/choukaihsuan/...）
raw_dir:     /mnt/d/circRNA_data/raw_fastq
trimmed_dir: /mnt/d/circRNA_data/trimmed
results_dir: /mnt/d/circRNA_results/GSE113230

genome:
  fasta:       /mnt/d/ref/hg19/hg19.fa
  gtf:         /mnt/d/ref/hg19/hg19.gtf
  bwa_index:   /mnt/d/ref/hg19/bwa_index/hg19
  hisat2_index: /mnt/d/ref/hg19/hisat2_index/hg19
  star_index:  /mnt/d/ref/hg19/star_index
  species:     hg19

ciriquant_config: config/ciriquant.yaml

consensus:
  tools:             [ciriquant, dcc]  # 工具選擇
  min_tools:         2                # 共識閾值
  slop:              10               # 座標容忍 bp
  min_bsj_reads:     2
  max_junction_ratio: 1.0             # pseudo-circ QC：BSJ/FSJ 上限（CIRIquant only）

circbase_file: ""   # 留空 = 自動下載；或填本地路徑

de:
  method:       edgeR_ciriquant    # edgeR_ciriquant / deseq2 / limma
  fdr_cutoff:   0.05
  log2fc_cutoff: 1.0
  tumor_label:  tumor
  normal_label: normal

threads: 8
```

**Server 版 config**（`~/circRNA_agent/config.yaml`）路徑前綴不同：
- `raw_dir: /home3/choukaihsuan/GSE113230/raw`
- `trimmed_dir: /home3/choukaihsuan/GSE113230/trimmed`
- `results_dir: /home3/choukaihsuan/GSE113230_results`
- `star_index: /home3/choukaihsuan/reference/hg19/star_index`

---

## HPC Server 環境

| 項目 | 值 |
|------|----|
| 主機 | `172.16.0.178`（也可用 `172.16.0.179`） |
| 使用者 | `choukaihsuan` |
| Home | `/home/choukaihsuan` → symlink 到 `/home3/choukaihsuan` |
| OS | CentOS 7 |
| Python | 3.7.12（conda env `ciriquant`） |
| CPU | 96 cores |
| RAM | 377 GB |
| 磁碟 | /home3：596 GB 可用 |
| Conda env | `ciriquant`（CIRIquant 1.1.3, DCC 0.5.0, STAR, HISAT2, BWA, samtools, snakemake） |
| Java | `/usr/bin/java`（不在 conda env 內，ciriquant.yaml 必須指定此路徑） |

**SSH 連線**：
```bash
ssh choukaihsuan@172.16.0.178
# 或用設定好的 alias
ssh genomics
```

**執行 pipeline**：
```bash
cd ~/circRNA_agent
conda activate ciriquant
nohup snakemake \
    --snakefile workflow/Snakefile \
    --configfile config.yaml \
    --cores 36 \
    --resources mem_gb=300 \
    --keep-going \
    --rerun-incomplete \
    > logs/pipeline_run.log 2>&1 &
echo "PID: $!"
```

**同步本機到 server**（config.yaml 不同步，路徑不同）：
```bash
rsync -avz /mnt/c/Users/User/develop/circRNA_agent/scripts/ choukaihsuan@172.16.0.178:~/circRNA_agent/scripts/
rsync -avz /mnt/c/Users/User/develop/circRNA_agent/workflow/ choukaihsuan@172.16.0.178:~/circRNA_agent/workflow/
```

---

## Web UI

Flask 應用，讓使用者在瀏覽器選擇分析方法：

```bash
cd ~/circRNA_agent
conda activate ciriquant
python scripts/web_ui.py --host 0.0.0.0 --port 5000
# 瀏覽器開 http://172.16.0.178:5000
```

**功能**：
- **GEO 一鍵啟動**（頂部卡片）：輸入 GSE ID + cores → POST `/run_gse` → 呼叫 `agent.py --gse {gse_id}` → 跳轉狀態頁
- Step 1：circRNA 工具選擇（CIRIquant / DCC / 兩者），自動顯示共識模式說明
- Step 2：DE 方法選擇（edgeR_ciriquant 推薦 / DESeq2 / limma-voom）
- Step 3：進階參數（min_bsj, slop, **max_junction_ratio**, FDR, log2FC, threads）
- 儲存設定 → 更新 `config.yaml`
- 儲存並執行 → 更新 config 後啟動 Snakemake subprocess
- 狀態頁（`/status`）：每 5 秒 auto-refresh log，顯示 pipeline 是否執行中

**Web UI routes**：
- `GET /` — 主設定頁
- `POST /update` — 儲存設定（+ 可選執行 Snakemake）
- `POST /run_gse` — GEO 一鍵啟動
- `GET /status` — 狀態頁（進度條 + 18 stage 格 + collapsible log）
- `GET /api/log` — log JSON（前端 polling 用）
- `GET /api/progress` — Snakemake log 解析 JSON（stages 陣列 + finished/total count + running bool）

---

## 已知問題與解決方法

| 問題 | 原因 | 解決方法 |
|------|------|----------|
| DCC `TypeError: NoneType has no len()` | DCC 0.5.0 不支援缺少 `-mt2` | 加入 `star_align_mate1` 和 `star_align_mate2` rule，分別比對 R1/R2 |
| DCC `command not found` | 指令為大寫 `DCC` | shell command 改為 `DCC` |
| DCC `-A` flag error | DCC 0.5.0 BAM 參數是 `-B` 不是 `-A` | 改為 `-B {bam}` |
| CIRIquant 輸出 `.bed` 而非 `.bsj` | CIRIquant 1.1.3 bioconda 版本差異 | rule output 只宣告 `.gtf`，忽略 `.bed` |
| multiqc ImportError TypedDict | multiqc 1.17 + markdown 3.6 不支援 Python 3.7 | `pip install markdown==3.3.7` |
| Snakemake LockException | 前次執行被 kill 留下 lock | `snakemake --unlock` 後重跑 |
| samtools sort "File exists" | 多個 CIRIquant 進程同時寫同一 sample 的暫存檔 | `pkill -f CIRIquant`，`rm -rf results/circRNA/`，重新跑 |
| `nohup: failed to run command 'snakemake'` | conda env 未啟動 | 先 `conda activate ciriquant` |
| wildcard ambiguity（star_align vs mate1/mate2） | `{srr}` wildcard 匹配到 `SRR7012366/mate1` | 在 `circrna.smk` 加 `wildcard_constraints: srr = r"[A-Z]+\d+"` |
| star_align temp dir 硬編路徑 | 原本 `/home/choukaihsuan/star_tmp/{srr}` 只適用本機 | 改為 `RESULTS_DIR + "/circRNA/{srr}/star_tmp"` |
| DCC `ValueError: invalid literal for int(): '4"'` | `-fg` 觸發 CircSkip 解析，GTF exon_number 屬性殘留尾部引號 | 移除 `-fg` flag，CircSkip 計數非必要 |
| DCC 六個樣本全部失敗（IndexError: list index out of range） | 多個 DCC 並行共用工作目錄下的 `_tmp_DCC/`，競爭條件 + 上次失敗殘留的 partial 資料 | 改為 `(cd {params.outdir} && DCC ...)` subshell，每個 sample 的 `_tmp_DCC/` 獨立在各自 outdir 內 |
| DCC log 路徑失敗（`logs/dcc/SRR.log: No such file or directory`） | `cd {outdir}` 後 `> {log}` 的相對路徑從 outdir 解析 | 同上，subshell 讓 `> {log}` 在 parent shell 的 CWD（`~/circRNA_agent/`）執行 |
| multiqc numpy 版本衝突（`numpy 1.16.4 < 1.17`） | conda env 安裝的 numpy 過舊，matplotlib 要求 ≥1.17 | `pip install 'numpy>=1.17'` 在 ciriquant env |
| `consensus_filter.py` `TypeError: 'type' object is not subscriptable` | Python 3.7：module-level 的 `CoordMap = dict[tuple[...], ...]` 賦值在執行時求值，不能用內建 `dict`/`tuple` 做 subscript；`from __future__ import annotations` 只延遲 annotation 求值，**不影響普通賦值** | 改為 `from typing import Dict, Tuple; CoordMap = Dict[Tuple[str,int,int], float]` |
| `parse_ciriquant` 回傳 0 筆 circRNA | CIRIquant 1.1.3 GTF 屬性欄用小寫 `bsj`/`fsj`，regex 搜尋大寫 `BSJ`/`FSJ` 全部未匹配 | 兩個 `re.search()` 加 `re.IGNORECASE` flag |
| `parse_dcc` 回傳 0 筆 circRNA | `CircCoordinates` 只有座標（8 欄，無 count）；舊程式讀 col 3 得到 Gene 欄（字串），`float()` 失敗後全部 skip | `parse_dcc()` 改為優先讀同目錄的 `CircRNACount`（col 3 = junction count）；不存在時 fallback count=5 |

---

## 重要的 Snakemake 行為

- **DAG 在啟動時編譯**：修改 `.smk` 規則不影響正在跑的進程
- **`script:` 腳本在執行時讀取**：修改 `analysis.R` 或 `.py` 腳本對正在跑的進程立即生效
- **`--rerun-incomplete`**：重跑上次被中斷的 job（不自動重跑已完成的）
- **`--keep-going`**：某個 rule 失敗時繼續跑其他獨立的 rule
- **資源控制**：`--cores 36 --resources mem_gb=300`（保留給其他 server 使用者）

---

## HTML 報告內容（`generate_report.py`）

報告為自包含 HTML（PDF 圖表以 base64 內嵌）：

| 區塊 | 說明 |
|------|------|
| Summary stat-boxes | 樣本數、total circRNAs、顯著數、Up/Down |
| **Type I/II 分類** | edgeR_ciriquant 模式才顯示；橫向進度條 + 各自數量；Type I = circRNA 專一性，Type II = 基因層次 |
| **Biomarker 候選表** | top 30，欄位：rank, circ_id, log2FC, padj, biomarker_score, in_circbase, circbase_id, circbase_gene, Type |
| Top DE table | FDR 顯著的 circRNA，含 Type 欄（若存在） |
| Volcano plot | **Plotly 互動式**（hover 顯示 circ_id / log2FC / padj / Type）；fallback to PDF embed if Plotly unavailable |
| PCA | **Plotly 互動式**（hover 顯示 SRR ID / condition，tumor/normal 顏色區分）；numpy SVD |
| Heatmap | **Plotly 互動式**（top 50 DE，hover 顯示 circRNA ID，RdBu_r colorscale，z-score 標準化）；fallback to PDF embed |

報告標頭顯示使用的 DE 方法（`method-tag` badge）。
Plotly 依賴：`plotly`、`numpy`；若兩者未安裝則自動 fallback 到靜態 PDF embed。

---

## 目前執行進度（2026-05-26）

| 步驟 | 狀態 |
|------|------|
| fastp QC/trim | ✅ 6/6 完成 |
| FastQC raw | ✅ 6/6 完成 |
| CIRIquant | ✅ 6/6 完成（GTF 在 `~/GSE113230_results/circRNA/SRRxxxxxxx/`） |
| STAR paired-end | ✅ 6/6 完成 |
| STAR mate1 | ✅ 6/6 完成 |
| STAR mate2 | ✅ 6/6 完成 |
| DCC | ✅ 5/6 完成；SRR7012368 重跑中（12:15 重啟） |
| consensus_filter（含 pseudo-circ QC） | ✅ 5/6 完成（1,594–3,157 circRNAs / sample）；等 SRR7012368 |
| merge_counts | ⏳ 等所有 consensus_filter 完成 |
| annotate_circbase | ⏳ |
| rank_biomarkers | ⏳ |
| DE analysis | ⏳ |
| isoform switching | ⏳ |
| report | ⏳ |

Pipeline 以 `nohup` 在背景執行，log 在 `~/circRNA_agent/logs/pipeline_run5.log`。
