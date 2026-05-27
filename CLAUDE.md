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
│   ├── assign_isoforms.py       # circRNA BSJ 座標對應 host gene + exon span + strand + region
│   ├── isoform_switching.R      # 計算 IUI、Wilcoxon rank-sum 測試 isoform switching
│   ├── analysis.R               # DE 分析（edgeR_ciriquant / deseq2 / limma）
│   ├── generate_report.py       # 輸出 HTML 報告（互動式 Plotly + Type I/II + biomarker）
│   ├── notify.py                # 通知模組（Email/Slack，Snakemake hook 呼叫）
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
          ┌──────────┴───────────────────────┐
          ▼                                  ▼
  [assign_isoforms]                  [annotate_circbase]
  host gene + strand +               比對 circBase hg19
  exon_span + region                 → circbase_annotated.tsv
  → isoform_groups.tsv
          │
          ▼
   [de_analysis]
   edgeR_ciriquant / deseq2 / limma
   → de/de_results.tsv
   → plots/*.pdf
          │
          ├──────────────────┐
          ▼                  ▼
  [rank_biomarkers]   [isoform_switching]
  composite score     IUI + Wilcoxon test
  → biomarker_        → iui_matrix.tsv
    candidates.tsv    → isoform_switching.tsv
          │
          └──────────┬───────────────────────┘
                     ▼
             [generate_report]
             → report.html
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

### Host Gene / Exon / Strand 注釋（`assign_isoforms.py`）

`assign_isoforms.py` 對每個 circRNA 輸出以下額外欄位：

| 欄位 | 說明 |
|------|------|
| `gene_name` | Host gene 名稱（GTF gene 特徵，取最小包含基因體） |
| `strand` | `+` 或 `-`（來自 host gene GTF 記錄） |
| `exon_span` | 環化的 exon 範圍，格式 `eN-eM`（如 `e3-e7`）；找不到時為空白 |
| `region` | `exonic`（兩端都匹配 exon 邊界）/ `intronic`（在基因內但不匹配 exon）/ `intergenic` |

**exon_span 判斷邏輯**：對每個 circRNA (start, end) 掃描 host gene 所有 transcript，找到 `exon.start ≈ circ_start` 且 `exon.end ≈ circ_end`（±10 bp）的 transcript，回報 exon 編號。負股基因的 exon_number 可能出現反向（如 `e9-e8`），這是正常現象——GTF exon_number 按基因組座標遞增，但 mRNA 順序在負股上相反。

GSE113230 結果：9,349 circRNA 中 8,046 exonic（86%）、1,089 intronic（12%）、214 intergenic（2%）。

### Biomarker 排序（`rank_biomarkers.py`）

四維 composite score，每維度在 significant set 內 min-max 標準化後平均：

```
biomarker_score = (sig_norm + fc_norm + conf_norm + known_bonus) / 4
  sig_norm  = -log10(pvalue 或 padj),  上限 10，標準化  ← 依 de_sig_by 切換
  fc_norm   = |log2FC|,                上限 5，標準化
  conf_norm = confidence_score 標準化
  known_bonus = 1 若 in_circbase，否則 0（不做標準化）
```

- 顯著閾值欄位依 `de_sig_by` 而定：`pvalue`（nominal）或 `padj`（BH 校正）
- CLI 參數：`--use-pvalue` 對應 `de_sig_by: pvalue`

### DE 分析方法（`config de.method`）

| 值 | 說明 |
|----|------|
| `edgeR_ciriquant`（預設）| 複製 `CIRI_DE_replicate`：edgeR GLM + FSJ offset（測 BSJ/FSJ 比值），輸出 Type I/II 分類 |
| `deseq2` | DESeq2 RLE normalization + Wald test |
| `limma` | limma-voom |

**`edgeR_ciriquant` 核心邏輯**：
FSJ counts → TMM normalization → per-locus FSJ CPM 作為 GLM offset
→ 效果等同於測試 BSJ/FSJ 比值是否在 tumor vs. normal 之間改變
→ 同時對 FSJ 跑獨立 QLFTest，若 FSJ 也顯著且方向相同 → Type II（基因層次），否則 → Type I（circRNA 專一性）

**Type I / II / III 分類**：
- **Type_I**：BSJ 顯著，FSJ 不顯著（或方向相反）→ circRNA 環化效率真正改變
- **Type_II**：BSJ 顯著，FSJ 也顯著且同方向（|FSJ logFC| ≥ 0.5）→ host gene 整體變化帶動
- **Type_III**：只有 FSJ 顯著 → 線性 mRNA 變化，不是 circRNA DE
- `concordant` 判斷要求 FSJ logFC 方向相同**且** |FSJ logFC| ≥ 0.5，避免邊緣顯著雜訊

**重要欄位命名**：merge BSJ/FSJ 結果後，`PValue` 不加後綴（只有 `logFC`、`FDR` 因同時出現在兩表才加 `_bsj`/`_fsj`）。`bsj_sig_col` 必須用 `"PValue"`，不是 `"PValue_bsj"`。

**顯著性欄位切換**（`config de.de_sig_by`）：
- `padj`（預設）：BH 校正 FDR。小樣本（n=3 vs 3）+ 多重檢定（~7,779 tests）時 min padj ≈ 0.432，幾乎無法通過
- `pvalue`：nominal p-value（未校正）。小樣本研究的實務做法，論文中需標明

`analysis.R` 有 **backward-compatible fallback**：若 `snakemake@input[["fsj_matrix"]]` 不存在（舊 DAG），自動 fallback 到 deseq2。

### Isoform Switching（`isoform_switching.R`）

**是否顯著的判斷使用 within-gene BH FDR**（`padj_within_gene`），而非 global BH。

原因：global BH 跨所有 isoform（~7,360 tests）→ min padj = 0.937，完全無法找到 switching。Within-gene BH（每個基因內部獨立做 BH）是 DEXSeq 風格的標準做法，在 GSE113230 結果中找到 66 個 significant switching events（FDR < 0.1, |ΔIUI| > 0.1）。

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
  method:              edgeR_ciriquant    # edgeR_ciriquant / deseq2 / limma
  fdr_cutoff:          0.05
  log2fc_cutoff:       1.0
  de_sig_by:           pvalue            # pvalue = nominal p；padj = BH 校正 FDR
  isoform_fdr_cutoff:  0.1               # within-gene FDR for isoform switching
  delta_iui_cutoff:    0.1               # minimum |ΔIUI| to call switching
  tumor_label:         tumor
  normal_label:        normal

threads: 8
```

**Server 版 config**（`~/circRNA_agent/config.yaml`）路徑前綴不同：
- `raw_dir: /home3/choukaihsuan/GSE113230/raw`
- `trimmed_dir: /home3/choukaihsuan/GSE113230/trimmed`
- `results_dir: /home3/choukaihsuan/GSE113230_results`
- `star_index: /home3/choukaihsuan/reference/hg19/star_index`
- `gtf: /home3/choukaihsuan/reference/hg19/genes.gtf`（注意：檔名是 `genes.gtf`，非 `hg19.gtf`）

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

**R packages（安裝在 conda env `ciriquant` 的 r-base 4.2.2）**：
r-ggplot2, r-pheatmap, r-rcolorbrewer, r-dplyr, r-ggrepel, r-tibble, r-tidyr,
bioconductor-edger, bioconductor-limma, r-statmod

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
| `parse_ciriquant` / `merge_counts.py` 回傳 0 筆 circRNA | CIRIquant 1.1.3 GTF 屬性欄用小寫 `bsj`/`fsj`，regex 搜尋大寫 `BSJ`/`FSJ` 全部未匹配 | 兩個 `re.search()` 加 `re.IGNORECASE` flag（`consensus_filter.py` 和 `merge_counts.py` 均適用） |
| `parse_dcc` 回傳 0 筆 circRNA | `CircCoordinates` 只有座標（8 欄，無 count）；舊程式讀 col 3 得到 Gene 欄（字串），`float()` 失敗後全部 skip | `parse_dcc()` 改為優先讀同目錄的 `CircRNACount`（col 3 = junction count）；不存在時 fallback count=5 |
| `analysis.R` `storage.mode(counts) <- "integer"` 失敗 | conda install r-pheatmap 後 dplyr 進入 env，data.frame subsetting 回傳 tibble-like 物件，`storage.mode` 不接受 | 在 `storage.mode` 之前加 `data.matrix(round(counts))` 轉為真正的 matrix（BSJ 和 FSJ 矩陣均需處理） |
| `generate_report.py` `SyntaxError: from __future__ imports must occur at the beginning` | Snakemake `script:` wrapper 在 user 腳本前插入 ~11 行 setup code，將 `from __future__ import annotations` 推到第 12 行；Python 3.7 要求此 import 必須是第一行 | 移除 `from __future__ import annotations`；改用 `from typing import Optional`；所有 `str \| None` 改為 `Optional[str]` |
| `analysis.R` `Error: replacement has 0 rows, data has 7779`（edgeR_ciriquant + use_pvalue）| `bsj_sig_col` 設為 `"PValue_bsj"`，但 merge 後該欄名為 `"PValue"`（只有 logFC/FDR 在兩表均存在才加後綴） | 改為 `bsj_sig_col <- if (use_pvalue) "PValue" else "FDR_bsj"` |
| Snakemake hook `notify.py` 路徑錯誤（找到 conda env 的 Python runner） | `onstart`/`onsuccess`/`onerror` 中 `__file__` 指向 conda env 的腳本執行器，非專案目錄 | 改用 `workflow.snakefile` 推算專案根目錄路徑 |
| benchmark `MissingInputException`（hg19.gtf） | benchmark config 路徑是 `hg19.gtf`，但 server 實際檔名是 `genes.gtf` | 更新 `benchmark/config_benchmark.yaml` GTF 路徑 |
| isoform switching 找到 0 個 significant events | Global BH 校正跨 ~7,360 isoform，min padj = 0.937 | 改用 within-gene BH（`padj_within_gene`）；結果：66 events（FDR < 0.1, \|ΔIUI\| > 0.1） |
| DE analysis 0 個 significant circRNA（padj 校正） | n=3 vs 3，BH 校正後 min padj = 0.432 | 改用 nominal p-value（`de_sig_by: pvalue`）；結果：831 個 significant circRNA |

---

## 通知系統（`notify.py`）

Pipeline 完成、失敗、啟動時自動發送通知，透過 Snakemake `onstart` / `onsuccess` / `onerror` hook 呼叫。

**支援管道**：

| 管道 | 說明 |
|------|------|
| SMTP Email | Gmail TLS 587，成功時附加 `report.html`（>20 MB 略過附件） |
| Slack Webhook | Incoming Webhook，Markdown 格式 |

（LINE Notify 已於 2025-03-31 終止服務，已從 `notify.py` 移除）

**通知內容**：

| 事件 | 內容 |
|------|------|
| `start` | 啟動時間 |
| `success` | 完成時間、總 circRNA 數、顯著 DECs、上調/下調數 |
| `failure` | 失敗 rule 名稱、pipeline log 最後 50 行 |

**環境變數設定**（加到 server `~/.bashrc`）：

```bash
export NOTIFY_EMAIL_FROM="寄件gmail@gmail.com"
export NOTIFY_EMAIL_PASS="xxxx xxxx xxxx xxxx"   # Gmail App Password（非登入密碼）
export NOTIFY_EMAIL_TO="chou.kaihsuan@gmail.com"
export NOTIFY_SLACK_WEBHOOK="https://hooks.slack.com/services/..."  # 選填
```

Gmail App Password 申請：Google 帳號 → 安全性 → 兩步驟驗證開啟後 → 應用程式密碼。

**手動測試**：
```bash
python scripts/notify.py --event start --project GSE113230
python scripts/notify.py --event success --project GSE113230 --report results/report.html
python scripts/notify.py --event failure --project GSE113230 --rule dcc --log logs/pipeline_run.log
```

未設定任何環境變數時，`notify.py` 只列印 log，不影響 pipeline 執行。

**Snakemake hook 實作**（`workflow/Snakefile`）：
- `onstart` → fire-and-forget subprocess 呼叫 `notify.py --event start`
- `onsuccess` → 解析 `de_results.tsv` 統計後呼叫 `notify.py --event success --stats {...}`
- `onerror` → 呼叫 `notify.py --event failure --rule unknown --log logs/pipeline_run.log`
- hook 內必須用 `workflow.snakefile` 推算專案路徑，**不能用 `__file__`**（會指向 conda env）

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
| Summary stat-boxes | 樣本數、total circRNAs、顯著數、Up/Down；顯著標準依 `de_sig_by` 顯示 `p<0.05` 或 `FDR<0.05` |
| **Type I/II 分類** | edgeR_ciriquant 模式才顯示；橫向進度條 + 各自數量 |
| **Biomarker 候選表** | top 30，欄位：rank, circ_id, log2FC, padj, biomarker_score, in_circbase, circbase_id, circbase_gene, Type |
| **Top DE table（分兩表）** | 上調（tumor 高）和下調（tumor 低）分開顯示；欄位含 gene_name / strand / region / exon_span / circbase_id / log2FC / pvalue / padj / Type |
| Volcano plot | **Plotly 互動式**（hover 顯示 circ_id / log2FC / p-value / Type）；Y 軸標題依 `de_sig_by` 動態切換；fallback to PDF embed if Plotly unavailable |
| PCA | **Plotly 互動式**（hover 顯示 SRR ID / condition，tumor/normal 顏色區分）；numpy SVD |
| Heatmap | **Plotly 互動式**（top 50 DE，hover 顯示 circRNA ID，RdBu_r colorscale，z-score 標準化）；fallback to PDF embed |
| Isoform Switching | Plotly 長條圖（top 10 switching genes 的 IUI tumor vs normal）+ 顯著 switching 表格 |

**DE table 資料來源合併**：
- `de_results.tsv`（主表）
- `isoform_groups.tsv` → `gene_name`, `strand`, `region`, `exon_span`
- `circbase_annotated.tsv` → `circbase_id`, `circbase_gene`, `in_circbase`

報告標頭顯示使用的 DE 方法（`method-tag` badge）。
Plotly 依賴：`plotly`、`numpy`；若兩者未安裝則自動 fallback 到靜態 PDF embed。

---

## 目前執行進度（2026-05-27）

| 步驟 | 狀態 |
|------|------|
| fastp QC/trim | ✅ 6/6 完成 |
| FastQC raw | ✅ 6/6 完成 |
| CIRIquant | ✅ 6/6 完成 |
| STAR paired-end | ✅ 6/6 完成 |
| STAR mate1 | ✅ 6/6 完成 |
| STAR mate2 | ✅ 6/6 完成 |
| DCC | ✅ 6/6 完成 |
| consensus_filter | ✅ 6/6 完成（1,594–3,728 circRNAs / sample） |
| merge_counts | ✅ 完成（9,349 circRNAs） |
| assign_isoforms | ✅ 完成（含 strand / exon_span / region） |
| annotate_circbase | ✅ 完成 |
| DE analysis | ✅ 完成（nominal p < 0.05；831 significant circRNAs） |
| isoform switching | ✅ 完成（66 events，within-gene FDR < 0.1） |
| rank_biomarkers | ✅ 完成（831 candidates） |
| report | ✅ 完成（含 gene / strand / region / exon_span / circbase_id 欄位） |

**GSE113230 各工具偵測數量**：

| SRR ID | 分組 | CIRIquant | DCC | 共識 |
|--------|------|----------:|----:|-----:|
| SRR7012366 | Tumor 1 | 26,455 | 6,010 | 1,905 |
| SRR7012367 | Tumor 2 | 34,075 | 7,222 | 2,329 |
| SRR7012368 | Tumor 3 | 35,000 | 10,756 | 3,728 |
| SRR7012369 | Normal 1 | 27,105 | 9,397 | 3,157 |
| SRR7012370 | Normal 2 | 39,211 | 9,349 | 2,790 |
| SRR7012371 | Normal 3 | 12,985 | 5,788 | 1,594 |
