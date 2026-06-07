# circRNA Analysis Pipeline — Project Context for Claude

## 專案概述

本專案是一個以 **Snakemake** 驅動的 circRNA（環狀 RNA）全流程分析管線，
從 GEO/SRA 原始數據下載，到差異表現分析（DE）與 HTML 報告輸出。

- **目標數據集**：GSE113230（三陰性乳癌 tumor vs. normal，6 個 sample）；GSE58135（乳癌，進行中）
- **主要工具**：CIRIquant（circRNA 偵測）+ DCC（輔助偵測，雙工具共識）
- **執行環境**：基因體中心 HPC server（`172.16.0.178`，CentOS 7，96 cores，377 GB RAM）
- **本機開發**：Windows 11 + WSL2（Ubuntu 26.04），程式碼在 `/mnt/c/Users/User/develop/circRNA_agent/`
- **Server 路徑**：`~/circRNA_agent/`（即 `/home3/choukaihsuan/circRNA_agent/`，`/home/choukaihsuan` 是 symlink）
- **Container**：Docker image `choukaihsuan/circrna-pipeline:1.0.0`；HPC 用 Singularity 拉取

---

## 目錄結構

```
circRNA_agent/
├── config.yaml                  # 主設定檔（路徑、參數、工具選擇）
├── Dockerfile                   # 容器化環境定義（mamba + conda env circrna.yaml）
├── .dockerignore                # Docker 建置排除清單
├── config/
│   ├── ciriquant.yaml           # CIRIquant 工具路徑設定（server 版本另存在 server 上）
│   ├── .ciriquant_ready         # touch 檔，驗證 ciriquant.yaml 存在後建立
│   └── projects/                # 各 GSE 專案獨立設定快照（由 web_ui.py 自動建立）
│       └── {GSE_ID}.yaml
├── containers/
│   └── build_and_deploy.sh      # Docker image 建置與推送腳本
├── metadata/
│   ├── library_info.csv         # 目前啟用專案的 SRR/配對資訊
│   ├── sample_groups.csv        # 目前啟用專案的分組
│   └── {GSE_ID}/                # 各 GSE 專屬 metadata 目錄（由 web_ui.py 自動建立）
│       ├── library_info.csv
│       └── sample_groups.csv
├── workflow/
│   ├── Snakefile                # 主 Snakefile，載入 rules，設定 target
│   └── rules/
│       ├── download.smk         # SRA 下載（aria2c S3 > ascp > prefetch + fasterq-dump）
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
│   ├── predict_interactions.py  # CircInteractome miRNA/RBP interaction 查詢
│   ├── notify.py                # 通知模組（Email/Slack，Snakemake hook 呼叫）
│   ├── utils.py                 # 共用工具函數
│   ├── web_ui.py                # Flask Web UI（GEO 一鍵啟動 + 進度視覺化）
│   └── templates/
│       ├── index.html           # 主設定頁面（GEO 入口 + Step 1-3）
│       └── status.html          # Pipeline 狀態頁（進度條 + rule 狀態格 + log）
├── envs/
│   └── circrna.yaml             # Conda 環境定義（mamba 建置；已移除 defaults channel）
└── logs/                        # Snakemake 各 rule 的 log 檔
```

---

## Pipeline 流程

```
SRA/GEO
  │
  ▼
[download] aria2c (S3, 16連線) > ascp (Aspera) > prefetch (HTTPS fallback)
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

**Adaptive fallback**（`--adaptive` flag）：
若兩工具偵測數量嚴重失衡（`min(counts) / max(counts) < --adaptive-ratio`，預設 0.1），
自動將 `min_tools` 從 2 降為 1，避免共識過濾後近乎零回收率。
失衡時在 stderr 印出警告。新 sample 類型的探索性分析時建議開啟。

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

六維 composite score，每維度在 significant set 內 min-max 標準化後平均：

```
biomarker_score = (sig_norm + fc_norm + conf_norm + known_bonus + mirna_norm + rbp_norm) / 6
  sig_norm   = −log10(pvalue 或 padj), 上限 10，標準化  ← 依 de_sig_by 切換
  fc_norm    = |log2FC|,               上限 5，標準化
  conf_norm  = confidence_score 標準化
  known_bonus = 1 若 in_circbase，否則 0（不標準化）
  mirna_norm = distinct miRNA binders 數，min-max 標準化（無 interaction data 者為 0）
  rbp_norm   = distinct RBP binders 數，min-max 標準化（無 interaction data 者為 0）
```

- 顯著閾值欄位依 `de_sig_by` 而定：`pvalue`（nominal）或 `padj`（BH 校正）
- CLI 參數：`--use-pvalue` 對應 `de_sig_by: pvalue`
- miRNA/RBP interaction 資料來自 `predict_interactions.py`（CircInteractome 查詢，top 50 circRNA）

### DE 分析方法（`config de.method`）

`analysis.R` **預設三種方法全跑**，各輸出獨立 TSV：
- `de_results_edgeR_ciriquant.tsv`（同時作為主 `de_results.tsv`）
- `de_results_deseq2.tsv`
- `de_results_limma.tsv`

| 值 | 說明 |
|----|------|
| `edgeR_ciriquant`（預設主方法）| 複製 `CIRI_DE_replicate`：edgeR GLM + FSJ offset（測 BSJ/FSJ 比值），輸出 Type I/II 分類 |
| `deseq2` | DESeq2 RLE normalization + Wald test |
| `limma` | limma-voom（TMM normalization）|

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
| Conda env | `ciriquant`（CIRIquant 1.1.3, DCC 0.5.0, STAR, HISAT2, BWA, samtools, snakemake, **aria2c 1.36.0**） |
| Java | `/usr/bin/java`（不在 conda env 內，ciriquant.yaml 必須指定此路徑） |

**SRA 下載優先順序**（`workflow/rules/download.smk`）：

| 優先 | 方法 | 速度 | 說明 |
|------|------|------|------|
| 1 | **aria2c + S3**（預設）| ~25 MB/s（16 連線）| `srapath --location s3` 取 S3 URL → aria2c 多連線下載 |
| 2 | ascp（Aspera）| ~50 MB/s | 需 Aspera key，目前 server 未安裝 |
| 3 | prefetch（HTTPS）| ~0.5 MB/s | NCBI 單連線，最慢，S3 失敗時的 fallback |

`_find_tool("aria2c")` 搜尋優先順序：`sra_env` → `circrna` → **`ciriquant`**（已加入）→ `which()`。
若 `srapath returned no S3 URL`（NCBI API 暫時失敗），重啟 pipeline 即可；S3 URL 通常幾分鐘後恢復。

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
- Step 2：**三種 DE 方法全部執行**（固定）；選擇「報告預設顯示方法（主方法）」；各方法附適用情境說明（edgeR=circRNA 特異性、DESeq2=樣本多/保守、limma=小樣本穩定）
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

## Container 部署（Docker + Singularity HPC）

### Docker image

```bash
# 本機建置並推送（WSL2）
cd /mnt/c/Users/User/develop/circRNA_agent
bash containers/build_and_deploy.sh   # 自動執行 docker build + docker push
```

image name：`choukaihsuan/circrna-pipeline:1.0.0`
`Dockerfile` 使用 `mamba env create --yes`（非互動式），比 conda 快 3–5×。

### `envs/circrna.yaml` 重要設定

- channels：`bioconda`, `conda-forge`（**已移除 `defaults`**，Anaconda 商業授權問題）
- `bioconductor-edger>=3.40`（不鎖定 patch 版本，避免 PackagesNotFoundError）
- `bioconductor-deseq2>=1.40`
- `bioconductor-qvalue`（Storey q-value，analysis.R 需要）
- `dcc=0.5.0` 在 conda 依賴中（bioconda），**不在 pip**（PyPI 的 DCC 是另一個無關套件，版本 0.7+）

### Singularity / Apptainer（HPC server）

Docker image 已上傳至 Docker Hub，任何有 Singularity 或 Apptainer 且開啟 user namespace 的 HPC server，執行以下一行即可取得完整環境：

```bash
singularity pull circrna-pipeline.sif docker://choukaihsuan/circrna-pipeline:1.0.0
```

或使用 Apptainer（語法相同）：

```bash
apptainer pull circrna-pipeline.sif docker://choukaihsuan/circrna-pipeline:1.0.0
```

拉取後用容器執行 pipeline：

```bash
snakemake \
    --snakefile workflow/Snakefile \
    --configfile config.yaml \
    --cores 36 \
    --use-singularity \
    --singularity-args "--bind /home3/choukaihsuan:/home3/choukaihsuan" \
    --keep-going --rerun-incomplete
```

`workflow/Snakefile` 頂部已設定：`singularity: "docker://choukaihsuan/circrna-pipeline:1.0.0"`

**注意**：目前使用的 server（172.16.0.178，CentOS 7）預設關閉 user namespace，conda 版 apptainer 無法執行。需請管理員執行：
```bash
echo 10000 > /proc/sys/user/max_user_namespaces
```
或安裝 setuid 版本的 Singularity。在此之前繼續使用 conda env `ciriquant`。

---

## Per-project Config 系統

每個 GSE 分析自動使用獨立設定，避免多個分析互蓋 `config.yaml` 和 `metadata/`。

### 機制

- **`save_project_snapshot(cfg)`**（`web_ui.py`）：每次儲存設定時，同時寫入：
  1. `config.yaml`（全域啟用，目前執行的專案）
  2. `config/projects/{GSE_ID}.yaml`（該專案的永久快照）
- **`_configfile_for(gse_id)`**：若 `config/projects/{GSE_ID}.yaml` 存在，Snakemake 的 `--configfile` 自動指向它；否則 fallback 到 `config.yaml`
- **`run_manual()`**（手動填入 SRR 啟動）：metadata 同時儲存到 `metadata/library_info.csv` 和 `metadata/{GSE_ID}/library_info.csv`；config 的 `metadata`/`groups` 欄位更新為專案路徑

### 切換專案

```bash
# Web UI 輸入新 GSE ID → run_gse 自動載入對應 config/projects/{GSE_ID}.yaml
# 若要手動切換：
cp config/projects/GSE113230.yaml config.yaml
```

### Server 端獨立重跑腳本（Mock snakemake 物件）

當 server 上 config.yaml 已切換到其他專案時，用以下 wrapper 直接執行，繞過 Snakemake DAG：

**`/tmp/run_generate_report.py`**（Python mock，`snakemake` = `types.SimpleNamespace`）：
```bash
# 必須用 conda run，否則 base python 缺少 plotly，會生成靜態 PDF 報告
conda run -n ciriquant python /tmp/run_generate_report.py
```

**`/tmp/run_de_analysis.R`**（R mock，S4 class with list slots）：
```bash
conda run -n ciriquant Rscript /tmp/run_de_analysis.R
```
R mock 格式：
```r
setClass('Snakemake', representation(input='list', output='list', params='list', log='list'))
snakemake <- new('Snakemake',
  input  = list(matrix="...", fsj_matrix="...", groups="...", circbase_annot="..."),
  output = list(de="...", de_edger="...", de_deseq="...", de_limma="...", volcano="...", ...),
  params = list(de_method="edgeR_ciriquant", fdr=0.05, lfc=1.0, ...),
  log    = list("/path/to/log")
)
source('/home/choukaihsuan/circRNA_agent/scripts/analysis.R')
```

兩個 wrapper 都 hardcode GSE113230 的 input/output/params 路徑。此模式可複用於任何需要獨立重跑 terminal rule 的情境。

---

## Benchmark 設計（`benchmark/`）

評估本 pipeline 對抗三個已發表方法的偵測準確率與 DE 品質，輸出自包含 HTML 報告。

### 執行方式

```bash
cd ~/circRNA_agent
conda activate ciriquant
snakemake \
    --snakefile benchmark/Snakefile \
    --configfile benchmark/config_benchmark.yaml \
    --cores 8 \
    --resources mem_gb=60 \
    --keep-going --rerun-incomplete
```

只跑準確率（跳過 DE quality，適合 GSE113230 尚未完成時）：
```bash
snakemake --snakefile benchmark/Snakefile \
    --configfile benchmark/config_benchmark.yaml \
    --cores 8 --until accuracy_benchmark compute_cost
```

### Ground Truth 建立（`rnaser_ground_truth.py`）

使用 **GSE55872**（Hs68 cell line, hg19）的 RNase R enrichment 實驗作為 ground truth：

| SRR ID | 角色 |
|--------|------|
| SRR444655 | Total RNA（偵測基準）|
| SRR444974 | RNase R replicate 1 |
| SRR445016 | RNase R replicate 2 |

RNase R 消化線性 RNA 並富集環狀 RNA，是目前 circRNA 偵測 benchmark 的標準方法。

**Enrichment Ratio（ER）判斷**：
- ER = RNase R BSJ counts / Total RNA BSJ counts
- ER > 1.5 → **True Positive**（真實 circRNA）
- ER < 0.5 → **True Negative**（假陽性候選）
- 中間值排除於評估之外

GSE55872 的 FASTQs 由 `bench_download` rule 從 **EBI FTP** 自動下載，無需手動準備。

### Task 1 – 偵測準確率比較（`accuracy_benchmark.py`）

對 SRR444655（total RNA）執行多工具共識策略，再對照 RNase R ground truth：

| 策略 | 工具組合 | slop | pseudo-circ QC | 對應論文 |
|------|----------|------|----------------|----------|
| **Our_adaptive** | CIRIquant + DCC | 10 bp | ✅ selective（BSJ<5）+ adaptive | — |
| **Our_no_QC**（消融）| CIRIquant + DCC | 10 bp | ❌ | — |
| **CirComPara2_4tools** | CIRIquant + DCC + CIRCexplorer2 + find_circ | 10 bp | ❌ | Gaffo et al. 2022 |
| **nfcore_3tools** | CIRIquant + CIRCexplorer2 + find_circ | 0（精確匹配）| ❌ | Digby-Bell et al. 2023 |
| CirComPara2_sim（ablation）| CIRIquant + DCC only | 10 bp | ❌ | — |

**重要設計決策**：CirComPara2_4tools 不含獨立 CIRI2，因為 CIRIquant 內部已呼叫 CIRI2 做 BSJ 偵測；同時納入兩者等於讓同一演算法投兩票，破壞 consensus 獨立性。4 個工具分別使用 HISAT2+BWA / STAR / STAR / Bowtie2，是真正獨立的偵測策略。

**find_circ 輸出過濾**：`parse_find_circ()` 必須只保留 `category` 欄含 `CIRCULAR` 的行，否則 LINEAR（~233K）和 AMBIGUOUS（~190K）junction 會進入 consensus，導致座標匹配極慢（測試發現 55 min 仍未完成）。真正 CIRCULAR 只有 1,895 個。

**評估指標**：Precision、Recall、F1、**Specificity**、**TN**、AUC-PR

**分層分析**：依 Total RNA 中的 BSJ count 分三層：
- Low：1–4 RPM
- Mid：5–19 RPM
- High：≥ 20 RPM

輸出：`benchmark/accuracy_summary.tsv`（含 TN、Specificity 欄）、`benchmark/stratified_f1.tsv`

### Task 2 – DE 分析品質比較（`de_quality_benchmark.py`）

使用 **GSE113230**（三陰性乳癌，6 samples）的 count matrix，比較三種 DE 方法：

| 方法 | DE 演算法 | 特色 |
|------|----------|------|
| **Our method** | edgeR_ciriquant（BSJ/FSJ ratio）| Type I/II 分類；FSJ offset |
| DESeq2 baseline | DESeq2 on BSJ counts | 無 FSJ offset |
| **limma-voom** | limma-voom（TMM + voom weights）| 無 FSJ offset |

**前置條件**：需先執行 `/tmp/run_de_analysis.R`（或 Snakemake de_analysis rule）產生三個 TSV；以及 `/tmp/run_de_baseline.R` 執行 benchmark baseline。

**評估指標**：
- 各方法顯著 DE circRNA 數量
- 兩兩 Jaccard similarity（Our vs DESeq2 / Our vs limma / DESeq2 vs limma）
- Type I circRNA 中，僅我們方法偵測到的比例
- Top 20 DE circRNA 中 circBase 已知的比例

輸出：`benchmark/de_quality_summary.tsv`、`benchmark/de_jaccard.tsv`

### Task 3 – 計算資源成本（`compute_cost.py`）

整合各步驟實測時間（`/usr/bin/time -v` log），輸出 `benchmark/compute_cost.tsv`。

**nf-core 實測 time logs**（benchmark Snakefile 現在記錄以下步驟）：
- `logs/bench/time_ciriquant_{srr}.log`：CIRIquant（與 Our pipeline 共用）
- `logs/bench/time_circexplorer2_{srr}.log`：CIRCexplorer2（約 6 秒）
- `logs/bench/time_find_circ_map_{srr}.log`：bowtie2 unmapped（在 HPC NFS 上約 3h 44min）
- `logs/bench/time_find_circ_{srr}.log`：find_circ detection

**注意**：find_circ_map wall time 受 HPC NFS 磁碟 I/O 速度影響顯著（寫入 28GB unmapped reads），比 AWS SSD 環境慢。報告 Source 欄標注「measured on HPC with NFS storage」。若 time logs 不存在，自動 fallback 使用 Digby-Bell 2023 文獻值。

**已知實測時間**（SRR444655，8 cores，HPC NFS）：
- CIRCexplorer2：0:05.99（6 秒，0.13 GB RAM）
- find_circ_map：3:44:29（3.54 GB RAM）
- find_circ detect：2:41:05（5.14 GB RAM）
- CIRIquant：**11:41:07（47.2 GB RAM）** ← NFS I/O 瓶頸（HISAT2 5h8m + BWA-MEM 3h5m + CIRI2 56m + de novo quant 1h22m）

**CIRIquant 步驟分解**（NFS 環境，SAM/BAM 寫入放大效應）：
- HISAT2 genome alignment：00:01 → 05:09（5h 8min；unmapped.sam 124 GB 寫入 NFS）
- Gene abundance：05:09 → 05:22（13 min）
- BWA-MEM：05:22 → 08:27（3h 5min）
- CIRI2.pl detection：08:27 → 09:23（56 min）
- Build circular index：09:23 → 09:31（8 min）
- De novo HISAT2 alignment：09:31 → 10:15（44 min）
- BSJ/FSJ detection & quantification：10:15 → 11:42（87 min）

### 輸出報告

`benchmark/comparison_report.html`：自包含 HTML，整合三個面向（準確率、DE 品質、資源成本）的比較表與圖表。

**獨立重跑 benchmark 腳本**（繞過 Snakemake，適合 config 指向其他專案時）：
```bash
RESULTS='/home3/choukaihsuan/GSE113230_results'
BENCH=$RESULTS/benchmark
SCRIPTS=~/circRNA_agent/benchmark/scripts
PY=/home/choukaihsuan/miniconda3/envs/ciriquant/bin/python
GS55=/home3/choukaihsuan/GSE55872_results/circRNA/SRR444655

# 0. nf-core 3-tool re-filter（需 find_circ CIRCULAR-only fix）
$PY ~/circRNA_agent/scripts/consensus_filter.py \
    --cirique $GS55/SRR444655.gtf --circexplorer2 $GS55/CIRCexplorer2/known_circ.txt \
    --find-circ $GS55/find_circ/splice_sites.bed \
    --output $BENCH/detection/nfcore_3tools.bed --summary $BENCH/detection/nfcore_3tools_summary.tsv \
    --min-tools 2 --slop 0 --min-bsj 2 --max-junction-ratio 999

# 1. CirComPara2 4-tool consensus（CIRIquant + DCC + CIRCexplorer2 + find_circ）
$PY ~/circRNA_agent/scripts/consensus_filter.py \
    --cirique $GS55/SRR444655.gtf --dcc $GS55/DCC/CircCoordinates \
    --circexplorer2 $GS55/CIRCexplorer2/known_circ.txt \
    --find-circ $GS55/find_circ/splice_sites.bed \
    --output $BENCH/detection/circompara2_4tools.bed --summary $BENCH/detection/circompara2_4tools_summary.tsv \
    --min-tools 2 --slop 10 --min-bsj 2 --max-junction-ratio 999

# 2. accuracy benchmark
$PY $SCRIPTS/accuracy_benchmark.py \
    --ground-truth $BENCH/rnaser/ground_truth.tsv \
    --our-bed $BENCH/detection/our_method.bed --our-summary $BENCH/detection/our_method_summary.tsv \
    --our-no-qc-bed $BENCH/detection/our_no_qc.bed --our-no-qc-summary $BENCH/detection/our_no_qc_summary.tsv \
    --circompara2-bed $BENCH/detection/circompara2_sim.bed --circompara2-summary $BENCH/detection/circompara2_sim_summary.tsv \
    --circompara2-4tools-bed $BENCH/detection/circompara2_4tools.bed --circompara2-4tools-summary $BENCH/detection/circompara2_4tools_summary.tsv \
    --nfcore-bed $BENCH/detection/nfcore_3tools.bed --nfcore-summary $BENCH/detection/nfcore_3tools_summary.tsv \
    --output-summary $BENCH/accuracy_summary.tsv --output-stratified $BENCH/stratified_f1.tsv \
    --output-fp-comparison $BENCH/fp_score_comparison.tsv --slop 10 --min-bsj 2

# 3. DE baseline（mock snakemake，見 /tmp/run_de_baseline.R）
conda run -n ciriquant Rscript /tmp/run_de_baseline.R

# 4. DE quality
$PY $SCRIPTS/de_quality_benchmark.py \
    --our-de $RESULTS/de/de_results.tsv \
    --nfcore-de $BENCH/de/nfcore_deseq2_results.tsv --limma-de $BENCH/de/nfcore_limma_results.tsv \
    --circbase-annot $RESULTS/circRNA/circbase_annotated.tsv \
    --fdr 0.05 --lfc 1.0 \
    --output-summary $BENCH/de_quality_summary.tsv --output-jaccard $BENCH/de_jaccard.tsv

# 5. comparison report
$PY $SCRIPTS/generate_comparison_report.py \
    --accuracy $BENCH/accuracy_summary.tsv --stratified $BENCH/stratified_f1.tsv \
    --compute $BENCH/compute_cost.tsv \
    --de-quality $BENCH/de_quality_summary.tsv --de-jaccard $BENCH/de_jaccard.tsv \
    --fp-comparison $BENCH/fp_score_comparison.tsv \
    --output $BENCH/comparison_report.html
```

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
| DE analysis 0 個 significant circRNA（padj 校正） | n=3 vs 3，BH 校正後 min padj = 0.432 | 改用 nominal p-value（`de_sig_by: pvalue`）；結果：482 個 significant circRNA |
| `bioconductor-edger=3.44.0` PackagesNotFoundError | conda bioconductor channel 無此 patch 版本 | 改為 `bioconductor-edger>=3.40`（不鎖 patch） |
| `DCC==0.5.0` pip install 失敗 | PyPI 上的 DCC 是不相干套件（版本 0.7+）；circRNA DCC 在 bioconda | 移至 conda 依賴 `dcc=0.5.0`，從 pip section 移除 |
| Docker build `mamba env create` 等待確認 | `mamba env create` 預設互動式 `[Y/n]` | 加 `--yes` flag：`mamba env create --yes` |
| Anaconda `repo.anaconda.com` 商業授權警告 | `defaults` channel 要求商業許可 | 移除 `defaults`，僅保留 `bioconda` + `conda-forge` |
| `generate_report` 重跑觸發上游 rules（download/align） | server config.yaml 指向其他專案，`--forcerun` 重建整個 DAG | 改用 `/tmp/run_generate_report.py` mock snakemake 物件直接執行 |
| SVG badge 距離弧段過遠（舊 angular de-overlap） | `deOverlap` 只往順時針推移，多個 badge 連鎖推移後遠離原弧段 | 改為 radial staggering：badge 固定在弧段正上方角度，僅在 angular gap < 0.14 rad 時改變半徑（`MI_R0=+7` vs `MI_R1=+19`）；虛線放在 `<g>` 內隨 badge 隱藏 |
| `analysis.R` heatmap 錯誤（`pheatmap` 少於 2 行） | DE circRNA 數量太少時 `pheatmap` 會報錯 | 加 `if (nrow(mat) >= 2)` 判斷，否則 `plot.new()` 顯示提示文字 |
| `analysis.R` `slice_min` 版本警告 | dplyr 新版棄用 `slice_min` 的部分用法 | 改為 `arrange(.data[[col]]) %>% head(n)` |
| `generate_report` Plotly 圖表沒有互動（靜態 PDF embed） | 用 base conda Python（3.13）執行，無 plotly | 改用 `conda run -n ciriquant python /tmp/run_generate_report.py` |
| modal 分頁無法點擊（整個 JS 失效） | `priTitle` 字串中的 Python `\n` 變成 JS 字串裡的實際換行 → SyntaxError | 改為 `\\n`，讓 JS 收到合法逸脫序列 |
| Isoform switching bar chart x-axis label 重疊 | 10 基因 × 2 條件 = 20 個 label 太密 | 改用 Plotly multi-level x-axis（`[[gene,gene],["Normal","Tumor"]]`） |
| `analysis.R` limma `object 'logFC' not found` | dplyr 1.1.x `rename(log2FC = logFC)` NSE 在 ciriquant env 下找不到欄位 `logFC`；DESeq2 和 edgeR 用 base R 賦值故無問題 | limma 區塊改用 base R：`names(res_df)[names(res_df) == "logFC"] <- "log2FC"` |
| benchmark `--forcerun comparison_report` 觸發 bowtie2 重跑 | `comparison_report` 依賴 `accuracy_summary` → `nfcore_sim.bed` → `bench_find_circ` → `bench_find_circ_map`；任何 `--forcerun` 包含 `comparison_report` 都會觸發整個 dependency chain | 改用直接腳本執行（見 Benchmark 獨立重跑腳本）；或只 `--forcerun de_deseq2_baseline de_quality_benchmark` |
| benchmark `de_deseq2_baseline` groups 路徑失效 | `benchmark/config_benchmark.yaml` 的 `groups` 指向全域 `metadata/sample_groups.csv`；config.yaml 切換到 GSE58135 後該檔案是 GSE58135 的 sample list | 更新 `benchmark/config_benchmark.yaml` 的 `groups` 改為 `metadata/GSE113230/sample_groups.csv` |
| Snakemake `--forcerun de_analysis` 觸發完整 DAG 重建 | `--forcerun` 標記 de_analysis 需重跑，但 Snakemake 評估 transitive inputs；若 `raw_dir` FASTQ 不存在則排入 download_fastq | 改用 `/tmp/run_de_analysis.R` mock snakemake 直接執行，或確保 `--configfile` 指向正確專案且所有上游 output 都存在 |
| `FULL_HEATMAP_DATA.conditions` 全為灰色（wrong SRR IDs）| Server 上 `config.yaml` 切換到 GSE58135，`metadata/sample_groups.csv` 被蓋成 GSE58135 SRR IDs；`/tmp/run_generate_report.py` 的 `groups` 參數指向此被蓋掉的檔案 | 改為指向 `metadata/GSE113230/sample_groups.csv`（固定路徑不受全域 config 影響）|
| 報告中 circRNA 以座標 ID 儲存，搜尋 circBase 名稱找不到 | 管線使用基因組座標 `chr:start\|end` 作為 circ_id，不是 circBase ID | 用 `circbase_annotated.tsv` 做 ID 對應；`annotate_circbase.py --slop 10` 容許 ±10 bp 誤差（1-based/0-based 轉換造成的 1 bp 差異在容許範圍內）|
| find_circ.py Python 3 三個 bug（`NameError: _np`、`mismatches` 縮排錯誤、`complement()` KeyError: `'\n'`）| Python 2→3 遷移不完整：`mismatches()` 內部的 `import numpy as _np` 縮排在 inner function，外層程式碼無法存取；`complement()` 用 `COMPLEMENT[x]` dict lookup 但基因組序列讀取包含 `\n`；`fromstring(a,dtype=byte)` Python 3 棄用 | 修 `/home3/choukaihsuan/tools/find_circ/find_circ.py`：(1) 模組層級加 `import numpy as _np`；(2) `mismatches()` 改用純 Python `sum(x!=y for x,y in zip(a,b))`；(3) `complement()` 加 `s = s.replace(chr(10),"").replace(chr(13),"")`  |
| `parse_find_circ()` 載入 235K 行造成 consensus_filter 卡住 55+ 分鐘 | find_circ 輸出含 LINEAR（233K）、AMBIGUOUS（190K）、CIRCULAR（僅 1,895）三類 junction，舊 parser 未過濾 category 欄 | `parse_find_circ()` 加判斷：`if len(parts)>=18 and "CIRCULAR" not in parts[17]: continue`；真正 circRNA 從 235K 降至 1,895 個，執行時間從 55min+ 降至秒級 |
| RBP 分頁 Binding Seq 欄大多顯示 "N/A" | ENCORI 的 `circ_pos` 是絕對座標格式（`chr6:148390208-148390208`），CircInteractome 是相對座標（`156–191`）；`_absPos()` 把 ENCORI 絕對座標再加 `chromStart`，產生超出染色體長度的位置 → UCSC API 回傳無效 → "N/A" | `_absPos()` 和 `_seq_logo` 加偵測：若 `circ_pos` 以 `chr` 開頭則視為絕對座標直接使用，否則才加 `chromStart`（`generate_report.py`） |
| benchmark `config_benchmark.yaml` BWA index 路徑錯誤 | `bwa_index: /home3/choukaihsuan/reference/hg19/bwa_index/hg19` 目錄不存在；實際 BWA index prefix 是 `hg19.fa` | 改為 `bwa_index: /home3/choukaihsuan/reference/hg19/hg19.fa` |
| find_circ_map timing bowtie2 卡住 13+ 小時 | timing script 將 28GB unmapped reads 寫入 `/tmp`（掛載在 `/` 分區，98% 使用率），磁碟 I/O 極慢 | 改寫到 `/home3/choukaihsuan/timing_tmp/`（348GB 可用）；`/` 分區的 timing_tmp 清理後釋放 11GB |

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
| **DE 方法切換器** | 頁面頂部三個按鈕（edgeR_ciriquant / DESeq2 / limma-voom）；切換時即時更新 stat-boxes、Volcano、Heatmap、**Top DE 表格**、**Biomarker 灰化**（`Plotly.react` 原地更新，無頁面 reload）|
| Summary stat-boxes | 樣本數、total circRNAs、顯著數、Up/Down；id="stat-n-sig/up/dn"，供方法切換器更新 |
| **Type I/II 分類** | edgeR_ciriquant 模式才顯示；橫向進度條 + 各自數量 |
| **3-method Venn diagram** | SVG 三圓 Venn；圓心 A=(170,128) B=(290,128) C=(230,210) r=90；交集數字加白色 halo（`paint-order="stroke"`）；高度 345px |
| **Biomarker 候選表** | top 30 表格；方法切換時非顯著行灰化（opacity 0.25）|
| **Biomarker Score 分布圖** | 兩圖並排：① Ranked scatter（x=rank, y=score，Top 30 紅點）② Histogram + Normal fit + Shapiro-Wilk 檢定；垂直線 μ / μ±σ / μ±2σ；說明文字在圖下方 |
| **Top DE table（分兩表）** | 方法切換時完全重繪（`_renderDETables()`）；欄位含 gene_name / strand / region / exon_span / circbase_id / log2FC / p-value / Type |
| Volcano plot | **Plotly 互動式**；方法切換時 Plotly.react 更新 |
| PCA | **Plotly 互動式**（tumor/normal 顏色區分）；numpy SVD |
| **Heatmap top-N 控制** | 預設 10+10，最高 50+50；tumor=紅、normal=**綠**（`#2CA02C`） |
| Heatmap | **Plotly 互動式**（top N DE，z-score 標準化）；方法切換時同步更新；Plotly 內建 title 移除，由 HTML h2 顯示 |
| Isoform Switching | Plotly 長條圖 + 顯著 switching 表格；**不隨 DE 方法切換更新**（IUI 計算固定）；section 下方有說明文字 |
| **SVG Circular Diagram** | 每個 circRNA 的環狀圖（exon 結構 + miRNA/RBP binding site 弧段 + 流水號 badge）|

**JS 全域狀態變數**：
- `const ALL_DE_METHODS`：三方法的完整 volcano + stats + heatmap + **de_table + sig_ids** 資料
- `const FULL_HEATMAP_DATA`：pool=50 up + 50 down，每 circRNA 含 `{z, pval, log2fc, label}`
- `let _HEATMAP_DATA_CACHE`：目前顯示方法的 heatmap data（方法切換時更新）；若 `conditions` 為空則 fallback 到 `FULL_HEATMAP_DATA.conditions`

**`switchDEMethod(method)`**：更新 stat-boxes → Plotly.react volcano → updateMainHeatmap → **`_renderDETables()`** → **`_updateBiomarkerHighlight()`**

**`_renderDETables(method, md)`**：從 `md.de_table` 重建 up/down HTML 表格，插入 `#de-tables-section`；清除靜態 `table, h3, .tbl-dl-bar`。

**`_updateBiomarkerHighlight(sigIds)`**：在 `#biomarker-section` 中，對每個 `<tr>` 讀取 `circ-link` onclick 的 circ_id，不在 sigIds 內者設 `opacity:0.25`。

**列印排版（`@media print`）**：
- `h2, h3 { break-after: avoid }` — 標題後不換頁
- `.plotly-graph-div { break-inside: avoid }` — 圖表不跨頁
- `#de-tables-section { break-before: page }` — Top DE 表格從新頁開始
- flex 區塊改 block（雙欄分布圖垂直排列）

**circRNA 詳細 modal（點擊任意 circ_id 開啟）**：
- **⬛ Circular Structure**：SVG 環狀圖；底部「⬇ SVG」下載按鈕
- **📺 miRNA Sponge**：互動表格（Priority 排序）；「⬇ CSV」下載；Binding Seq 欄自動從 UCSC hg19 REST API 獲取序列
- **🧬 RBP Binding**：同上
- **📈 Volcano / 🔥 Heatmap**：Plotly mini-chart；「⬇ PNG」下載
- **Priority Score**（miRNA）：seed type（8mer=+4…6mer=+1）+ CLIP>0（+3）+ ENCORI（+2）+ in_circ（+1）
- **Priority Score**（RBP）：bindingSites log₂×2 max3 + internal（+2）+ CLIP>0（+3）+ ENCORI（+2）+ in_circ（+1）
- 所有 header 可點擊排序（▲/▼）
- Chr Position 欄顯示 chr 絕對座標（circRNA start + circ_pos 1-based offset）

**Volcano toggle**：右上角「○ Heatmap circles: OFF / ON」按鈕，heatmap top circRNA 標記預設隱藏。

**全域下載功能**：
- 頁面頂部 sticky 列：「🖨 列印 / 存為 PDF」（`window.print()`）
- 各表格右上角「⬇ CSV」按鈕：DE 上調/下調、Biomarker、Isoform switching

**Plotly 安裝要求**：`pip install plotly` 在 `ciriquant` conda env（已安裝 5.18.0）。
必須用 `conda run -n ciriquant python` 執行 generate_report.py，否則 fallback 靜態 PDF。

**SVG Circular Diagram — Badge 放置邏輯（radial staggering）**：

- miRNA badge（外圈）：預設放在弧段外 `+7` px（`MI_R0`）；若與前一 badge 角度差 < 0.14 rad，改放 `+19` px（`MI_R1`）
- RBP badge（內圈）：同邏輯，往內 `−8` px（`RBP_R0`）或 `−20` px（`RBP_R1`）
- **不做 angular de-overlap**：badge 角度永遠等於弧段中心角，只在徑向上交錯
- 虛線 connector（`stroke-dasharray="2,2"`）放在 `<g id="_mib_N">` 內，toggle 時隨 badge 一起隱藏
- 虛線只在 badge 在外圈（`MI_R1`）時才畫，距離弧段緊貼時不畫

**DE table 資料來源合併**：
- `de_results.tsv`（主表）
- `isoform_groups.tsv` → `gene_name`, `strand`, `region`, `exon_span`
- `circbase_annotated.tsv` → `circbase_id`, `circbase_gene`, `in_circbase`

報告標頭顯示使用的 DE 方法（`method-tag` badge）。
Plotly 依賴：`plotly`、`numpy`；若兩者未安裝則自動 fallback 到靜態 PDF embed。

---

## 目前執行進度（2026-06-05 更新）

### GSE113230（三陰性乳癌）

**所有步驟已完成（含三方法 DE + 新版報告 + benchmark）。** 報告位置：`~/GSE113230_results/report.html`（server）

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
| merge_counts | ✅ 完成（9,349 circRNAs；filterByExpr 後 4,630） |
| assign_isoforms | ✅ 完成（含 strand / exon_span / region） |
| annotate_circbase | ✅ 完成 |
| DE analysis | ✅ 完成（三方法全跑：edgeR 482 / DESeq2 409 / limma 736 significant）|
| predict_interactions | ✅ 完成（top 50；CircInteractome；interactions.json） |
| isoform switching | ✅ 完成（66 events，within-gene FDR < 0.1） |
| rank_biomarkers | ✅ 完成（482 candidates；**6D score**：sig+FC+conf+circbase+miRNA+RBP） |
| report | ✅ 完成 v2（動態 DE 表格切換；Biomarker 分布圖 + 常態檢定；Venn diagram 修正；列印排版）|
| benchmark accuracy | ✅ 完成（4-method + Our_no_QC ablation；report.html 更新）|
| benchmark compute cost | ✅ 完成（CIRIquant 實測 11:41:07 on HPC NFS；compute_cost.tsv + comparison_report.html 已更新）|

**主要數值結果**：
- 偵測：9,349 consensus circRNAs → filterByExpr 後 4,630 tested
- DE（edgeR_ciriquant）：482 significant（nominal p < 0.05，|log2FC| > 1）；min Storey q = 0.384（underpowered）
- DE（DESeq2 baseline）：409 significant；DE（limma-voom）：736 significant
- Isoform switching：66 events（within-gene FDR < 0.1，|ΔIUI| > 0.1）
- Biomarker score：6D（sig, FC, confidence, circBase, #miRNA, #RBP），87 circRNAs 有 interaction data（interactions.json May 29 更新）
- Top 1 biomarker：chr10:5836848|5842668（hsa_circ_0002665，GDI2；score=0.8202，83 miRNA，118 RBP binders，log2FC=7.48，Type_I）；interactions.json 於 May 29 重跑後更新，118 RBP binders 為資料集最大值（rbp_n=1.0）
- Benchmark（含 CirComPara2_4tools + Our_no_QC 消融）：

| Method | Precision | Recall | F1 | Specificity | AUC-PR |
|--------|-----------|--------|----|-------------|--------|
| Our_adaptive | 0.877 | 0.171 | 0.286 | 0.959 | 0.946 |
| Our_no_QC | 0.879 | 0.173 | 0.290 | 0.959 | 0.946 |
| CirComPara2_sim（2-tool ablation）| 0.879 | 0.173 | 0.290 | 0.959 | 0.946 |
| **CirComPara2_4tools** | 0.852 | **0.235** | **0.368** | 0.930 | 0.921 |
| nfcore_3tools | 0.873 | 0.182 | 0.301 | 0.955 | 0.943 |

CirComPara2_4tools Recall 最高但 Specificity 最低；Our pipeline AUC-PR 最優（0.946）。

**Biomarker score 公式（6D）**：
```
score = (sig_norm + fc_norm + conf_norm + known_bonus + mirna_norm + rbp_norm) / 6
  sig_norm   = −log10(pvalue), 上限 10，min-max 標準化
  fc_norm    = |log2FC|, 上限 5，min-max 標準化
  conf_norm  = confidence_score, min-max 標準化
  known_bonus = 1 若 in_circbase，否則 0（不標準化）
  mirna_norm = distinct miRNA binders 數，min-max 標準化（無 interaction data 者為 0）
  rbp_norm   = distinct RBP binders 數，min-max 標準化（無 interaction data 者為 0）
```

**Cascade 顯著性邏輯（de_sig_by: auto）**：
1. filterByExpr(min.count=5) → 4,630 tests（原 9,349）
2. Storey q-value < 0.2（bioconductor-qvalue，已安裝）
3. min q = 0.384 → fallback to nominal p < 0.05
4. 482 significant circRNAs（論文 Methods 需說明樣本量限制）

**GSE113230 各工具偵測數量**：

| SRR ID | 分組 | CIRIquant | DCC | 共識 |
|--------|------|----------:|----:|-----:|
| SRR7012366 | Tumor 1 | 26,455 | 6,010 | 1,905 |
| SRR7012367 | Tumor 2 | 34,075 | 7,222 | 2,329 |
| SRR7012368 | Tumor 3 | 35,000 | 10,756 | 3,728 |
| SRR7012369 | Normal 1 | 27,105 | 9,397 | 3,157 |
| SRR7012370 | Normal 2 | 39,211 | 9,349 | 2,790 |
| SRR7012371 | Normal 3 | 12,985 | 5,788 | 1,594 |

---

### GSE58135（乳癌）

**進行中。** SRA 下載中（10 個 SRR）；fasterq-dump 完成後繼續跑 QC → circRNA 偵測。

| 步驟 | 狀態 |
|------|------|
| prefetch + fasterq-dump | 🔄 進行中（10/10 SRR 已 prefetch） |
| fastp QC/trim | ⏳ 待執行 |
| CIRIquant | ⏳ 待執行 |
| STAR / DCC | ⏳ 待執行 |
| consensus → DE → report | ⏳ 待執行 |

**Server config**（`config/projects/GSE58135.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE58135/raw`
- `results_dir: /home3/choukaihsuan/GSE58135_results`
