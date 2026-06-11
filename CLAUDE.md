# circRNA Analysis Pipeline — Project Context for Claude

## 專案概述

本專案是一個以 **Snakemake** 驅動的 circRNA（環狀 RNA）全流程分析管線，
從 GEO/SRA 原始數據下載，到差異表現分析（DE）與 HTML 報告輸出。

- **目標數據集**：GSE113230（三陰性乳癌 tumor vs. normal，6 samples，✅ 完成）；GSE58135（乳癌，10 samples，✅ 完成）；GSE323364（TNBC cell line EZH2 inhibitor，6 samples，✅ 完成）；GSE133998（乳癌 tumor vs. normal，12 samples，✅ 完成）
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
│   ├── prepare_metadata.py      # 從 GEO/SRA RunInfo 建立 library_info.csv；T/N/APN 尾碼自動偵測
│   ├── download_geo.py          # SRA metadata 抓取（GSE→pysradb；PRJNA/SRP→NCBI eUtils API）
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
→ 同時對 FSJ 跑獨立 QLFTest，若 FSJ 也顯著（不做方向一致性檢查）→ Type II（兩層調控同時發生），否則 → Type I（circRNA 專一性）

**Type I / II / III 分類**：
- **Type_I**：BSJ 顯著，FSJ 不顯著 → circRNA 環化效率真正改變（circRNA-specific regulation）
- **Type_II**：BSJ/FSJ ratio 顯著（BSJ offset test）且 FSJ 獨立測試也顯著 → 兩種調控同時發生
- **Type_III**：只有 FSJ 顯著 → 線性 mRNA 變化，不是 circRNA DE
- **不做方向一致性檢查**：offset approach 造成 log2FC（ratio 方向）與 logFC_fsj（FSJ 方向）天然反向相關（FSJ↑ → offset↑ → ratio log2FC 傾向負值），強制要求 sign 一致等同於把幾乎所有 Type II 歸入 Type I。正確判斷：`sig_bsj AND sig_fsj`（無方向條件）
- `fsj_concordance_lfc`（config `de.fsj_concordance_lfc`，預設 `0.0`）：可設定最小 |logFC_fsj| 作為 FSJ 顯著的過濾條件（0 = 僅看統計顯著性）

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

| 優先 | 方法 | 實測速度 | 說明 |
|------|------|----------|------|
| 1 | **aria2c + S3**（預設）| **~25–30 MB/s**（6–8 連線，NFS 環境）| `srapath --location s3` 取 S3 URL → aria2c 多連線下載；`-x 6 -s 6` 適合 NFS（過多連線反而限速）|
| 2 | ascp（Aspera）| ~50 MB/s | 需 Aspera key，目前 server 未安裝 |
| 3 | prefetch（HTTPS）| **~20–100 KB/s**（實測 22 KB/s）| NCBI 單連線，最慢；656MB 需 8+ 小時；S3 失敗時的 fallback |

**後處理速度（6GB SRA file，NFS）**：

| 步驟 | 工具 | 實測速度 | 說明 |
|------|------|----------|------|
| SRA → FASTQ | fasterq-dump `-e 6` | ~10–15 min/6GB | 多執行緒解壓縮；`--split-files` 輸出 PE |
| FASTQ 壓縮 | **pigz** `-p 3` | **~15 min/15GB** | 多核 gzip，較快 |
| FASTQ 壓縮 | gzip（單執行緒）| ~2h/15GB | pigz 未安裝時的 fallback |

**aria2c 並行下載建議**：
- **S3 per-IP 總連線數安全閾值**：≤ 30 連線（實測 32 連線可接受）；超過 60 會觸發 DNS SERVFAIL（AWS S3 throttle），持續約 2 小時
- **當前 download.smk 設定**：`-x 8 -s 8`；Snakemake `threads: 8`，max 4 parallel → 4 × 8 = 32 連線（安全範圍上緣）
- **勿手動腳本同時跑多個 SRR**：10 × 16 = 160 連線曾觸發 DNS 封鎖（GSE133998 教訓）
- S3 間歇性不可用（NCBI API）：`srapath --location s3` 回傳空字串，等待恢復後手動重啟
- 下載時若 kill 中斷：aria2c 會留下 `.sra.aria2` resume 檔 → 重啟可斷點續傳；但若多次中斷導致 `.aria2` 狀態損壞，需刪除 `.sra` 和 `.aria2` 重新下載，否則 fasterq-dump 報 `rcBlob,rcCorrupt`

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
- **手動 SRR 清單**：逐筆輸入 SRR ID + condition，或上傳 CSV（srr_id, condition）
- **本地 FASTQ**（新功能）：指定 server 上 FASTQ 目錄路徑 → 自動偵測配對檔（支援 `_1/_2`、`_R1/_R2`、`_R1_001/_R2_001` 等格式）→ 使用者在表格中指定各 sample 的 condition → 系統在 `raw_dir/` 建立 symlink（`{name}_1.fastq.gz` → 原始路徑），讓 Snakemake 自動跳過 download 步驟直接從 QC 開始
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
- `POST /run_manual` — 手動 SRR 清單或 CSV 上傳啟動
- `GET /api/scan_fastq?path=...` — 掃描 server 目錄，回傳偵測到的 FASTQ 配對 JSON
- `POST /run_local` — 本地 FASTQ 建立 symlink 後啟動 pipeline
- `GET /status` — 狀態頁（進度條 + 18 stage 格 + collapsible log）
- `GET /api/log` — log JSON（前端 polling 用）
- `GET /api/progress` — Snakemake log 解析 JSON（stages 陣列 + finished/total count + running bool）
- `GET /api/detect_labels?gse=...` — 自動偵測 case/control label

**GEO 資料集選擇指引**（`templates/index.html` 獨立折疊區塊）：
Web UI 主頁新增常駐卡片，使用者點標題行即可展開；內容包含：
- ✅/⚠️ 六項 checklist（Read length / Library type / replicates / case+control / 深度 / 樣本類型）
- 五行因素表（讀長 / RNA-Seq 方式 / 樣本類型 / 樣本數 / 設計）+ 顏色標記（綠/黃/紅）
- 已知資料集摘要（GSE113230/GSE58135/GSE323364/GSE133998 一覽）
- GEO 查詢提示（Library Strategy / avgLength 欄位位置）

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

**AUC-PR 計算方式（重要）**：
- 二元偵測（每個 circRNA 只有 detected/not detected，無連續分數）直接套用 `sklearn.average_precision_score` 會嚴重虛高 AUC-PR。
- **根本原因**：Our_adaptive 偵測到的 circRNA 全部 score=1，未偵測的全部 score=0；ground truth 中 87.8% 的 TP 在 score=0 的大池（即 FN）；樂觀排序（label=1 排在 label=0 前面）把 2,953 個 FN 全部排在 1,997 個 TN 前面，產生長段假精確度，AUC-PR = 0.946（虛高）。
- **正確方法**：門檻掃描（threshold sweep）。對 CIRI2 output 和 DCC `CircRNACount` 各取 `min_bsj` 門檻（1–50），重新建立共識 → 計算每個門檻下的 (Precision, Recall) → 繪製真實 PR 曲線 → trapezoid AUC。
- **實測結果**：Our_adaptive **誠實 AUC-PR = 0.155**（最大 Recall 17.3%，Precision 維持 90.5%）；其他方法的 AUC-PR 仍為二元偵測值（尚未做門檻掃描，待後續補充）。
- **新增 CLI 參數**：`accuracy_benchmark.py --ciri2-file <CIRI2 output> --dcc-count-file <CircRNACount> --output-pr-curve <pr_curve.tsv>`；`generate_comparison_report.py --pr-curve <pr_curve.tsv>`（在報告中插入 SVG PR 曲線圖）。
- **SRR444655 無 CIRIquant GTF**：benchmark total RNA 樣本（SRR444655）只有 RNase R replicates 才有 CIRIquant GTF；門檻掃描改用 `SRR444655.ciri2` + `DCC/CircRNACount`（CIRI2 輸出 col 4 = BSJ count，DCC CircRNACount col 3 = junction count）。

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
- CIRIquant：**11:50:25（49.1 GB RAM）** ← NFS I/O 瓶頸（HISAT2 5h8m + BWA-MEM 3h3m + CIRI2 + de novo quant）；2026-06-10 重跑確認值（disk-full 期間仍完成）

**CIRIquant 步驟分解**（NFS 環境，SAM/BAM 寫入放大效應）：
- HISAT2 genome alignment：00:01 → 05:09（5h 8min；unmapped.sam 124 GB 寫入 NFS）
- Gene abundance：05:09 → 05:22（13 min）
- BWA-MEM：05:22 → 08:25（3h 3min）
- CIRI2.pl detection：08:25 → 09:21（56 min）
- Build circular index：09:21 → 09:29（8 min）
- De novo HISAT2 alignment：09:29 → 10:13（44 min）
- BSJ/FSJ detection & quantification：10:13 → 11:50（97 min，含 disk-full 延遲）

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

# 2. accuracy benchmark（含門檻掃描 PR curve，需提供 CIRI2 output 和 DCC CircRNACount）
$PY $SCRIPTS/accuracy_benchmark.py \
    --ground-truth $BENCH/rnaser/ground_truth.tsv \
    --our-bed $BENCH/detection/our_method.bed --our-summary $BENCH/detection/our_method_summary.tsv \
    --our-no-qc-bed $BENCH/detection/our_no_qc.bed --our-no-qc-summary $BENCH/detection/our_no_qc_summary.tsv \
    --circompara2-bed $BENCH/detection/circompara2_sim.bed --circompara2-summary $BENCH/detection/circompara2_sim_summary.tsv \
    --circompara2-4tools-bed $BENCH/detection/circompara2_4tools.bed --circompara2-4tools-summary $BENCH/detection/circompara2_4tools_summary.tsv \
    --nfcore-bed $BENCH/detection/nfcore_3tools.bed --nfcore-summary $BENCH/detection/nfcore_3tools_summary.tsv \
    --output-summary $BENCH/accuracy_summary.tsv --output-stratified $BENCH/stratified_f1.tsv \
    --output-fp-comparison $BENCH/fp_score_comparison.tsv --slop 10 --min-bsj 2 \
    --ciri2-file $GS55/SRR444655.ciri2 \
    --dcc-count-file $GS55/DCC/CircRNACount \
    --output-pr-curve $BENCH/pr_curve.tsv

# 3. DE baseline（mock snakemake，見 /tmp/run_de_baseline.R）
conda run -n ciriquant Rscript /tmp/run_de_baseline.R

# 4. DE quality
$PY $SCRIPTS/de_quality_benchmark.py \
    --our-de $RESULTS/de/de_results.tsv \
    --nfcore-de $BENCH/de/nfcore_deseq2_results.tsv --limma-de $BENCH/de/nfcore_limma_results.tsv \
    --circbase-annot $RESULTS/circRNA/circbase_annotated.tsv \
    --fdr 0.05 --lfc 1.0 \
    --output-summary $BENCH/de_quality_summary.tsv --output-jaccard $BENCH/de_jaccard.tsv

# 5. comparison report（加 --pr-curve 可在報告中插入 SVG PR 曲線圖）
$PY $SCRIPTS/generate_comparison_report.py \
    --accuracy $BENCH/accuracy_summary.tsv --stratified $BENCH/stratified_f1.tsv \
    --compute $BENCH/compute_cost.tsv \
    --de-quality $BENCH/de_quality_summary.tsv --de-jaccard $BENCH/de_jaccard.tsv \
    --fp-comparison $BENCH/fp_score_comparison.tsv \
    --pr-curve $BENCH/pr_curve.tsv \
    --output $BENCH/comparison_report.html
```

---

## 已知問題與解決方法

| 問題 | 原因 | 解決方法 |
|------|------|----------|
| PRJNA/SRP metadata 抓取失敗（server 上）| Server 防火牆封鎖 NCBI eUtils HTTP；GSE 透過 pysradb 走不同 endpoint 可成功 | 見下方詳述 |
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
| `--adaptive` 沒有傳入 `consensus_filter.py` | `circrna.smk` 的 consensus_filter shell 指令從來沒帶 `--adaptive` flag；adaptive 只在腳本內部實作但未啟用 | `circrna.smk` 加入 `adaptive_flag = "--adaptive" if config["consensus"].get("adaptive", True) else ""`，預設開啟；同時加 `--adaptive-ratio` 參數 |
| `vst()` 在少於 1000 circRNA 時失敗（GSE58135）| DESeq2 的 `vst()` 要求輸入行數 ≥ nsub（預設 1000）；小資料集（如 28 circRNA）呼叫 `vst()` 報錯：`less than 'nsub' rows` | `analysis.R` 改為 `tryCatch(vst(dds, blind=FALSE), error=function(e) varianceStabilizingTransformation(dds, blind=FALSE))`，自動 fallback |
| GSE58135（50bp reads）consensus 只有 28 個 circRNA | DCC 在 50bp 讀長下每個 sample 只偵測到 3–17 個 circRNA（STAR chimeric 短讀限制），但 `--adaptive` 未傳入，min_tools=2 要求兩工具交集，結果全部 10 sample 只有 28 個 | 修復 `--adaptive` 傳入問題後，ratio ≈ 0.002 << 0.1，adaptive 自動降級 min_tools=1，結果 1,607 個 circRNA |
| `_biomarker_normality_plot` Shapiro-Wilk ValueError（n < 3）| GSE323364 只有 2 個 biomarker candidates（p < 0.05），`scipy.stats.shapiro` 要求 n ≥ 3，直接呼叫會 raise ValueError | `generate_report.py` 加 `if n < 3` guard：顯示「樣本數不足（n=N），無法進行常態檢定」而不執行 Shapiro-Wilk |
| Web UI 手動 SRR 表單 Project ID 預設填入當前專案 ID | `value="{{ config.get('project_id', 'CUSTOM') }}"` 讓使用者未改 ID 就送出，蓋掉現有專案的 metadata | 改為 `value="" required`，強制使用者填入正確 GSE ID |
| `config/projects/GSE133998.yaml` sra_cache_dir/tmp_dir 繼承舊專案路徑 | Web UI 從前一專案（GSE323364）快照複製 config，路徑殘留 `GSE323364/sra_cache`；Snakemake DAG 正常但中間檔案放錯位置 | 用 `sed -i` 替換 config 中的路徑（`GSE323364` → `GSE133998`）|
| gzip 壓縮 SRA 轉換 FASTQ 太慢（~2h/15GB）| fasterq-dump 後 gzip 單執行緒壓縮，NFS 寫入慢 | 在 ciriquant env 安裝 `pigz`（`conda install -y -c conda-forge pigz`）；`download.smk` 已使用 pigz（若 pigz 不存在，`find_tool("pigz")` fallback 到 gzip）|
| NCBI S3 URL 間歇性不可用（`srapath` 回傳空）| NCBI API 暫時失敗，`srapath --location s3` 回傳空字串；download.smk fallback 到 prefetch HTTPS（~22 KB/s）| 等待 S3 恢復後手動執行 aria2c；一般幾分鐘至數小時內恢復 |
| S3 aria2c DNS SERVFAIL（`DNS server returned general failure`）| 並行下載連線總數過多（例如 10 samples × 16 connections = 160）觸發 AWS S3 per-IP throttle，DNS 回傳 SERVFAIL；持續約 2 小時後自動恢復 | 限制總連線數 ≤ 30：`-x 8` × 4 parallel = 32（當前設定，可接受）；`-x 4` × 4 = 16（保守）；早上手動腳本 10 × 16 = 160 是觸發根因；恢復後重跑 pipeline 即可 |
| `kill <PID>`（SIGTERM）殺掉 Snakemake 時刪除 output 檔 | Snakemake 收到 SIGTERM 執行 graceful shutdown，自動刪除所有「正在執行 job 的 output 檔案」以防止不完整輸出；若早先下載完成的 FASTQ 恰好被列為「執行中 job 的 output」（因未設 protected），就會被刪除 | **一律用 `kill -9 <PID>`（SIGKILL）終止 Snakemake**，SIGKILL 無法被捕捉，不觸發 cleanup；永遠不用 bare `kill` 或 `pkill` 不加 `-9`；GSE133998 因此損失 9 個 FASTQ（~63 GB），需重新下載 |
| 修改 `.smk` 規則觸發所有 jobs 重跑（"Code has changed"）| Snakemake 預設追蹤 rule 的 code hash；修改任何 `.smk` 檔案後，Snakemake 在下次執行時標記該 rule 所有 jobs 為「Code has changed since last execution」→ 全部重跑，即使 output 已存在 | **只在 Snakemake 未執行時修改 `.smk`**；若需 mid-run 修改，加 `--rerun-triggers mtime` 旗標讓 Snakemake 只依時間戳判斷（忽略 code 變更）；GSE133998 因將 `-x 4` 改為 `-x 8` 而觸發所有 download_fastq 重跑 |
| `--rerun-incomplete` + 未 protected output + SIGTERM = output 消失 | `--rerun-incomplete` 把被 kill 前的「未完成 job」全部重排；若這些 job 的 output（早先其實已完成）未被 Snakemake 的 `protected()` 標記（原始 run 被 kill -9 前來不及 chmod），Snakemake 視為「正在執行」；再次 SIGTERM 時 cleanup 刪除這些 output | 避免在有執行中 job 的情況下用 `--rerun-incomplete` 重啟；若必須重啟，先手動 `chmod 444` 保護重要 output 再重跑，**並且**確認 Snakemake metadata（`.snakemake/metadata/`）也記錄該檔案為 protected（光 chmod 不足，需 Snakemake 自己完成的 job 才記錄）|
| Circular Structure 第一個分頁空白（無弧段）| interactions.json 的 `spliced_length` 欄位固定為 0；JS `if(totalLen>0)` guard 阻擋所有弧段繪製 | JS 改為：`const totalLen = parseInt(info.spliced_length) \|\| (exonBds.length>0 ? exonBds[exonBds.length-1].cum_end \|\| 0 : 0)`；從 `exon_boundaries[last].cum_end` 推算真實長度 |
| Venn detail 「方法」欄顯示三個方法（包含不顯著的）| `_venn_3_svg` 建立 `circ_info` 時對每個方法的**全部** DE 行都加入 `mlabel`，不論該 circRNA 在該方法中是否顯著 | 加入 `is_sig = cid in m_sig_set` 判斷；只有當 circRNA 在 `sig_sets[m]`（顯著集合）中才加入 `methods` 清單 |
| Heatmap tumor/normal 欄位混排 | 樣本欄依 SRR ID 字母順序排列，奇數偶數交錯（normal→tumor→normal…）| `_plotly_heatmap()` 和 `FULL_HEATMAP_DATA` 建立時，讀取 `sample_groups.csv` 按 condition 重排：tumor 欄在左、normal 欄在右 |
| `analysis.R` `coef=2` 在配對設計下指向 patient 係數而非 condition | `model.matrix(~patient+condition)` 有多個 patient dummy；`coef=2` 指向第一個 patient，非 condition | 改為 `cond_coef = ncol(design)`（condition 永遠在最後一列）；`glmQLFTest(fit, coef=cond_coef)`、`topTable(fit, coef=cond_coef)` |
| `download.smk` Python 3.7 不支援 `unlink(missing_ok=True)` | `missing_ok` 參數在 Python 3.8 才加入；CentOS 7 server 的 conda Python 3.7 執行時 TypeError | 改用 `try: lock.unlink() except OSError: pass` 相容寫法 |
| Type II DE circRNA 幾乎全為零 | 兩個根本原因：(1) `abs(logFC_fsj) >= 0.5` 閾值過嚴；(2) offset approach 造成 log2FC 與 logFC_fsj **反向相關**，`sign(log2FC)==sign(logFC_fsj)` 幾乎從不成立（73 sig_fsj 中只有 1 個同方向） | 移除方向一致性檢查，改為 `sig_bsj & sig_fsj`（BSJ ratio 顯著 AND FSJ 獨立顯著）；結果：GSE113230 從 0 → 73 Type_II（15.1%），GSE58135 → 2（13.3%），GSE133998 → 2（2.4%） |
| benchmark AUC-PR 0.946 虛高（二元偵測假象）| `sklearn.average_precision_score` 對 score ∈ {0,1} 做樂觀排序：ground truth 87.8% 的 TP 都在 score=0 大池（FN）；樂觀排序把 2,953 個 FN 全部排在 1,997 個 TN 前，產生假精確度長段，AUC-PR = 0.946 | 改用**門檻掃描**（threshold sweep）：對 CIRI2 + DCC CircRNACount 各取 min_bsj ∈ {1,2,3,...,50}，每個門檻計算真實 (Precision, Recall)，trapezoid AUC；誠實值 = **0.155**；新增 `--ciri2-file / --dcc-count-file / --output-pr-curve` CLI 參數至 `accuracy_benchmark.py`，`--pr-curve` 至 `generate_comparison_report.py` |
| CirComPara2 工具數誤標（5→4 工具）| `generate_comparison_report.py` 的說明文字列出 CIRI2 + CIRIquant + DCC + find_circ + CIRCexplorer2 = 5 工具；但 CIRIquant 內部已呼叫 CIRI2 做 BSJ 偵測，兩者是同一演算法 | 修正為 4 工具（CIRIquant + DCC + find_circ + CIRCexplorer2）；Consensus 門檻改為「≥2/4 tools」；報告說明文字同步更新 |

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
| **3-method Venn diagram** | SVG 三圓 Venn；圓心 A=(170,128) B=(290,128) C=(230,210) r=90；交集數字加白色 halo（`paint-order="stroke"`）；高度 345px；**7 個區域均可點擊**，點擊後在 `div#venn-detail` 顯示對應 circRNA 清單（含 gene_name / log2FC / p-value / circbase_id / **方法欄**）；「⬇ CSV」下載各區域清單 |
| **Biomarker 候選表** | top 30 表格；方法切換時非顯著行灰化（opacity 0.25）|
| **Biomarker Score 分布圖** | 兩圖並排：① Ranked scatter（x=rank, y=score，Top 30 紅點）② Histogram + Normal fit + Shapiro-Wilk 檢定；垂直線 μ / μ±σ / μ±2σ；說明文字在圖下方 |
| **Top DE table（分兩表）** | 方法切換時完全重繪（`_renderDETables()`）；欄位含 gene_name / strand / region / exon_span / circbase_id / log2FC / p-value / Type |
| Volcano plot | **Plotly 互動式**；方法切換時 Plotly.react 更新 |
| PCA | **Plotly 互動式**（tumor/normal 顏色區分）；numpy SVD |
| **Heatmap top-N 控制** | 預設 10+10，最高 50+50；tumor=紅、normal=**綠**（`#2CA02C`） |
| Heatmap | **Plotly 互動式**（top N DE，z-score 標準化）；方法切換時同步更新；Plotly 內建 title 移除，由 HTML h2 顯示；**欄位順序：tumor 在左、normal 在右**（按 `sample_groups.csv` condition 分類排列）|
| Isoform Switching | Plotly 長條圖 + 顯著 switching 表格；**不隨 DE 方法切換更新**（IUI 計算固定）；section 下方有說明文字 |
| **SVG Circular Diagram** | 每個 circRNA 的環狀圖（exon 結構 + miRNA/RBP binding site 弧段 + 流水號 badge）|

**JS 全域狀態變數**：
- `const ALL_DE_METHODS`：三方法的完整 volcano + stats + heatmap + **de_table + sig_ids** 資料
- `const FULL_HEATMAP_DATA`：pool=50 up + 50 down，每 circRNA 含 `{z, pval, log2fc, label}`；欄位順序 tumor 在左、normal 在右
- `let _HEATMAP_DATA_CACHE`：目前顯示方法的 heatmap data（方法切換時更新）；若 `conditions` 為空則 fallback 到 `FULL_HEATMAP_DATA.conditions`
- `const VENN_REGION_DATA`：7 個 Venn 區域的 circRNA 清單（`{label, circs:[{id, gene, lfc, pval, cb, m}]}`）

**`switchDEMethod(method)`**：更新 stat-boxes → Plotly.react volcano → updateMainHeatmap → **`_renderDETables()`** → **`_updateBiomarkerHighlight()`**

**`_renderDETables(method, md)`**：從 `md.de_table` 重建 up/down HTML 表格，插入 `#de-tables-section`；清除靜態 `table, h3, .tbl-dl-bar`。

**`_updateBiomarkerHighlight(sigIds)`**：在 `#biomarker-section` 中，對每個 `<tr>` 讀取 `circ-link` onclick 的 circ_id，不在 sigIds 內者設 `opacity:0.25`。

**列印排版（`@media print`）**：
- `h2, h3 { break-after: avoid }` — 標題後不換頁
- `.plotly-graph-div { break-inside: avoid }` — 圖表不跨頁
- `#de-tables-section { break-before: page }` — Top DE 表格從新頁開始
- flex 區塊改 block（雙欄分布圖垂直排列）

**circRNA 詳細 modal（點擊任意 circ_id 開啟）**：
- **⬛ Circular Structure**：SVG 環狀圖；底部「⬇ SVG」下載按鈕；`totalLen` 由 `exon_boundaries[last].cum_end` 推算（interactions.json 的 `spliced_length` 欄位固定為 0，不可用）
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

## 目前執行進度（2026-06-10 更新）

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
| report | ✅ 完成 v3（動態 DE 表格切換；Biomarker 分布圖 + 常態檢定；Venn diagram 修正；列印排版；**Venn 可點擊區域**；Heatmap tumor 在左；Circular Structure totalLen 修正）|
| benchmark accuracy | ✅ 完成（4-method + Our_no_QC ablation；report.html 更新）|
| benchmark compute cost | ✅ 完成（CIRIquant 實測 **11:50:25** on HPC NFS；2026-06-10 重跑確認；compute_cost.tsv + comparison_report.html 已更新）|

**主要數值結果**：
- 偵測：9,349 consensus circRNAs → filterByExpr 後 4,630 tested
- DE（edgeR_ciriquant）：482 significant（nominal p < 0.05，|log2FC| > 1）；**409 Type_I (84.9%) / 73 Type_II (15.1%)**；min Storey q = 0.384（underpowered）
- DE（DESeq2 baseline）：409 significant；DE（limma-voom）：736 significant
- Isoform switching：66 events（within-gene FDR < 0.1，|ΔIUI| > 0.1）
- Biomarker score：6D（sig, FC, confidence, circBase, #miRNA, #RBP），87 circRNAs 有 interaction data（interactions.json May 29 更新）
- Top 1 biomarker：chr10:5836848|5842668（hsa_circ_0002665，GDI2；score=0.8202，83 miRNA，118 RBP binders，log2FC=7.48，Type_I）；interactions.json 於 May 29 重跑後更新，118 RBP binders 為資料集最大值（rbp_n=1.0）
- Benchmark（含 CirComPara2_4tools + Our_no_QC 消融）：

| Method | Precision | Recall | F1 | Specificity | AUC-PR |
|--------|-----------|--------|----|-------------|--------|
| **Our_adaptive** | 0.877 | 0.171 | 0.286 | 0.959 | **0.155**（門檻掃描，誠實值）|
| Our_no_QC | 0.879 | 0.173 | 0.290 | 0.959 | 0.946†|
| CirComPara2_sim（2-tool ablation）| 0.879 | 0.173 | 0.290 | 0.959 | 0.946†|
| **CirComPara2_4tools** | 0.852 | **0.235** | **0.368** | 0.930 | 0.921†|
| nfcore_3tools | 0.873 | 0.182 | 0.301 | 0.955 | 0.943†|

†：二元偵測樂觀 AUC-PR（score ∈ {0,1}，虛高；Our_adaptive 誠實值 = 0.155，見 AUC-PR 計算說明）

CirComPara2_4tools Recall 最高但 Specificity 最低；Our pipeline Specificity 最優（0.959）。

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

**完成。** 報告位置：`~/GSE58135_results/report.html`（server）

**特殊情況：50bp reads + CIRIquant/DCC 嚴重失衡**
- Read length 50bp → STAR chimeric junction 偵測極差（DCC 僅 3–17 circRNA/sample）
- CIRIquant vs DCC 比例 ≈ 0.002，遠低於 adaptive_ratio=0.1 閾值
- `--adaptive` flag 先前未傳入 circrna.smk（bug 已修正）；修正後 adaptive fallback 觸發，min_tools 從 2 降為 1（CIRIquant-only 模式）
- 最終共識：1,607 circRNAs；filterByExpr 後測試數量較小，DESeq2 vst() 失敗已以 varianceStabilizingTransformation() fallback 修正

| 步驟 | 狀態 |
|------|------|
| prefetch + fasterq-dump | ✅ 10/10 完成 |
| fastp QC/trim | ✅ 10/10 完成 |
| CIRIquant | ✅ 10/10 完成 |
| STAR / DCC | ✅ 10/10 完成（DCC 偵測數極少，adaptive 模式僅用 CIRIquant） |
| consensus_filter（--adaptive）| ✅ 完成（1,607 circRNAs；CIRIquant-only 模式） |
| merge_counts / assign_isoforms | ✅ 完成 |
| DE analysis | ✅ 完成（edgeR 15 / DESeq2 122 / limma 508 significant）|
| report | ✅ 完成 |

**注意**：edgeR_ciriquant 顯著數極少（15）是因為 50bp reads → circRNA 偵測數量有限 + 樣本間差異較大；**13 Type_I (86.7%) / 2 Type_II (13.3%)**；limma-voom 在小樣本較穩定（508 significant）。

**Server config**（`config/projects/GSE58135.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE58135/raw`
- `results_dir: /home3/choukaihsuan/GSE58135_results`

**中間檔案已清理**（保留報告、count matrix、DE 結果）；raw FASTQ 已刪除以釋放磁碟空間。

---

### GSE323364（TNBC cell line，EZH2 inhibitor）

**完成。** 報告位置：`~/GSE323364_results/report.html`（server）

MDA-MB-436 TNBC cell line，EZH2 抑制劑 EPZ-6438 vs. DMSO，150bp PE，total RNA，各 3 replicates。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 6/6 完成 |
| fastp QC/trim | ✅ 6/6 完成 |
| CIRIquant | ✅ 6/6 完成 |
| STAR / DCC | ✅ 6/6 完成 |
| consensus_filter（--adaptive）| ✅ 完成 |
| merge_counts / assign_isoforms | ✅ 完成 |
| DE analysis | ✅ 完成（edgeR 15 / DESeq2 122 / limma 508 significant）|
| report | ✅ 完成 |

**設定**：
- case/control label：`EPZ6438` / `DMSO`
- SRR 清單：SRR37484804–SRR37484809（6 個 sample；3 EPZ6438 + 3 DMSO）
- genome：hg19（同 GSE113230）

**主要數值結果**：
- edgeR_ciriquant：2 significant circRNAs（nominal p < 0.05）；**2 Type_I (100%) / 0 Type_II**；EZH2 抑制劑對 circRNA 影響極有限（細胞株 + 藥物處理）
- Biomarker candidates：2 個（p < 0.05 篩選極嚴）→ `_biomarker_normality_plot` 需 n ≥ 3 的 Shapiro-Wilk 保護已加入

**中間檔案已清理**（raw FASTQ + trimmed + sra_cache + 中間 BAM 已刪除，釋放 77GB）；
報告 1.3MB，保留於 `~/GSE323364_results/report.html`。

**Server config**（`config/projects/GSE323364.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE323364/raw`
- `results_dir: /home3/choukaihsuan/GSE323364_results`

**Condition list CSV 格式**（6 行）：
```
srr_id,condition,gsm_id,description
SRR37484809,EPZ6438,GSM9564370,MDA-MB-436 EPZ-6438 rep1
SRR37484808,DMSO,GSM9564371,MDA-MB-436 DMSO rep1
SRR37484807,EPZ6438,GSM9564372,MDA-MB-436 EPZ-6438 rep2
SRR37484806,DMSO,GSM9564373,MDA-MB-436 DMSO rep2
SRR37484805,EPZ6438,GSM9564374,MDA-MB-436 EPZ-6438 rep3
SRR37484804,DMSO,GSM9564375,MDA-MB-436 DMSO rep3
```

---

### GSE133998（乳癌，paired tumor/normal）

**完成。** 報告位置：`~/GSE133998_results/report.html`（server）

乳癌手術切除組織，cancer tissue vs. adjacent normal，150bp PE，Illumina HiSeq X Ten，各 6 replicates（H36–H42，共 12 samples）。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 12/12 完成（aria2c S3 + prefetch fallback） |
| fastp QC/trim | ✅ 12/12 完成 |
| CIRIquant | ✅ 12/12 完成 |
| STAR paired-end | ✅ 12/12 完成 |
| STAR mate1 | ✅ 12/12 完成 |
| STAR mate2 | ✅ 12/12 完成 |
| DCC | ✅ 12/12 完成 |
| consensus_filter（--adaptive）| ✅ 完成（10,979 circRNAs 總計） |
| merge_counts / assign_isoforms | ✅ 完成 |
| DE analysis | ✅ 完成（edgeR 84 / DESeq2 194 / limma 674 significant）|
| report | ✅ 完成 v3（Venn 可點擊區域；Heatmap tumor 在左；Circular Structure 修正；**Venn 方法欄修正**）|

**主要數值結果**：
- 偵測：10,979 consensus circRNAs；filterByExpr 後 640 tested
- DE（edgeR_ciriquant）：84 significant（nominal p < 0.05，|log2FC| ≥ 1）；上調 34 / 下調 50；**82 Type_I (97.6%) / 2 Type_II (2.4%)**
- DE（DESeq2）：194 significant；DE（limma-voom）：674 significant
- Isoform switching：8,624 rows（within-gene BH FDR 分析）
- Biomarker candidates：84 個

**各 sample 共識 circRNA 數量**：

| SRR ID | 條件 | 共識 circRNA |
|--------|------|----------:|
| SRR11600329 | normal | 871 |
| SRR11600330 | tumor | 2,892 |
| SRR11600331 | normal | 1,765 |
| SRR11600332 | tumor | 597 |
| SRR11600333 | normal | 3,089 |
| SRR11600334 | tumor | 505 |
| SRR11600335 | normal | 4,013 |
| SRR11600336 | tumor | 808 |
| SRR11600337 | normal | 2,975 |
| SRR11600338 | tumor | 1,051 |
| SRR11600339 | normal | 3,307 |
| SRR11600340 | tumor | 2,848 |

**注意**：
- 配對設計（同患者 tumor+normal），嘗試加入 `design = ~patient + condition` 後：edgeR 變差（min padj=0.996，patient dummy 消耗 5 個 df，FSJ offset 已吸收個體差異）；limma 進步（51 padj<0.05）；**決定維持 unpaired 設計**，主要看 edgeR_ciriquant 結果。
- `analysis.R` 已加入向後相容的配對設計支援：若 `sample_groups.csv` 含 `patient_id` 欄則自動啟用 `~patient+condition`，否則維持 `~condition`；`cond_coef=ncol(design)` 確保 condition 係數位置正確。
- edgeR filterByExpr 後僅 640 circRNA（10,979 中）是因為 tumor/normal 樣本間表現量分佈差異較大。

**Server config**（`config/projects/GSE133998.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE133998/raw`
- `results_dir: /home3/choukaihsuan/GSE133998_results`

**設定**：
- case/control label：`tumor` / `normal`
- SRR 清單：SRR11600329–SRR11600340（12 samples；6 tumor + 6 normal）
- genome：hg19（同 GSE113230）
- Condition list CSV：`/mnt/c/Users/User/Desktop/GSE133998_condition_list.csv`

**SRR → 樣本對應**：

| SRR ID | 條件 | 樣本 |
|--------|------|------|
| SRR11600329 | normal | H36-adjacent normal |
| SRR11600330 | tumor | H36-cancer tissue |
| SRR11600331 | normal | H37-adjacent normal |
| SRR11600332 | tumor | H37-cancer tissue |
| SRR11600333 | normal | H38-adjacent normal |
| SRR11600334 | tumor | H38-cancer tissue |
| SRR11600335 | normal | H39-adjacent normal |
| SRR11600336 | tumor | H39-cancer tissue |
| SRR11600337 | normal | H40-adjacent normal |
| SRR11600338 | tumor | H40-cancer tissue |
| SRR11600339 | normal | H42-adjacent normal |
| SRR11600340 | tumor | H42-cancer tissue |

**Server config**（`config/projects/GSE133998.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE133998/raw`
- `results_dir: /home3/choukaihsuan/GSE133998_results`

**下載注意事項**：
- NCBI S3 間歇性不可用（`srapath --location s3` 回傳空），直接用 aria2c 手動觸發（S3 URL 恢復後）
- SRR11600334 初次下載失敗（HTTPS fallback ~22 KB/s），需重試
- `pigz` 已安裝於 ciriquant env（`conda install -y -c conda-forge pigz`），gzip 壓縮從 ~2h 降至 ~15min/15GB

---

## GEO / SRA BioProject 資料集選擇指引

### 支援的 Accession 格式

Pipeline 支援三種 accession 格式，輸入 `--gse` 或 Web UI 的「GEO 資料集」欄位：

| 格式 | 範例 | 來源 | 查詢方式 |
|------|------|------|---------|
| `GSE*` | `GSE113230` | NCBI GEO Series | pysradb（GSE → SRP → SRR） |
| `PRJNA*` | `PRJNA808398` | NCBI SRA BioProject | NCBI eUtils API |
| `SRP*` | `SRP156355` | NCBI SRA Study | NCBI eUtils API |

**HPC 網路限制**：Server（172.16.0.178）防火牆封鎖 NCBI eUtils HTTP；PRJNA/SRP 需在本機執行 `python scripts/download_geo.py --gse PRJNA808398` 取得 metadata CSV，再透過 Web UI 的「手動 SRR 清單 → 上傳 CSV」路徑匯入。GSE 透過 pysradb 可在 server 上直接執行。

**Condition 自動偵測**（`prepare_metadata.py`）：

| 樣本命名慣例 | 偵測邏輯 | 範例 |
|------------|---------|------|
| SRA T/N 尾碼 | `_detect_condition_by_suffix()`（優先）：名稱結尾為 `T` → tumor；`N`/`APN`/`CN` → normal | `M269T`→tumor, `87APN`→normal |
| GEO 關鍵字 | `_detect_condition()`（fallback）：掃描 disease_state/source_name/tissue 欄位 | `tumor tissue`→tumor |
| 無法偵測 | 留空（如 DCIS `D` 尾碼）| `1086D`→空，需手動指定 |

---

### 已分析/已調查資料集

| 編號 | 樣本類型 | 對比組 | Library | Read length | 樣本數 | 定序深度 | circRNA 偵測 | DE 顯著數（edgeR）| 評估 |
|------|----------|--------|---------|-------------|--------|----------|-------------|-------------------|------|
| **GSE113230** | 組織（TNBC tumor vs. normal）| tumor vs. normal | Total RNA（rRNA-depleted）| 150bp PE | 6（3T+3N）| ~100M reads/sample | 9,349（consensus）| **482** | ✅ 最佳示範 |
| **GSE58135** | 組織（乳癌 tumor vs. normal）| tumor vs. normal | Total RNA | **50bp PE** | 10（5T+5N）| ~50M reads/sample | 1,607（CIRIquant-only）| **15** | ⚠️ 50bp 讀長，DCC 失效→adaptive fallback |
| **GSE323364** | **細胞株**（MDA-MB-436 TNBC）| EPZ6438 vs. DMSO | Total RNA | 150bp PE | 6（3+3）| ~60M reads/sample | 中等 | **15** | ⚠️ 細胞株+藥物，DE 極少 |
| **GSE133998** | 組織（乳癌 tumor vs. normal）| tumor vs. normal | Total RNA | 150bp PE | 12（6T+6N）| ~80M reads/sample | 10,979 | **84** | ✅ 配對設計，樣本數最多 |
| **SRP156355** | 組織（早期乳癌 IDC/DCIS）| tumor vs. normal | rRNA-depleted（`other`）| **100bp PE** | 23（6T+6N+DCIS+APN）| avg 48M spots（~96M reads）| 待跑 | 待跑 | ✅ 待分析；6 對配對 T/N；LibrarySelection=other 適合 circRNA；Cancer Institute WIA 印度 |
| **PRJNA808398** | 組織（TNBC tumor vs. normal）| tumor vs. normal | **cDNA（poly-A）** | 150bp PE | 50（25T+25N）| avg 18M spots（低）| ❌ | ❌ | ❌ 不建議：poly-A selection，circRNA 接近零；A.C. CAMARGO 巴西 |

### 影響分析品質的關鍵因素

#### 1. Read Length（讀長）

| 讀長 | circRNA 偵測 | 說明 |
|------|-------------|------|
| **≥ 100bp PE**（建議）| ✅ 良好 | BSJ 需橫跨 back-splice junction，讀長越長越容易偵測；150bp 是目前標準 |
| 75bp PE | ⚠️ 尚可 | CIRIquant 可用，DCC 偵測能力下降 |
| **50bp PE** | ❌ 不佳 | STAR chimeric 幾乎無法偵測 BSJ；DCC 每 sample 只得 3–17 個 circRNA；需強制 CIRIquant-only 模式 |

**判斷方法**：GEO Series → Library Strategy 欄位；或下載 SRA run info 查 `avgLength` 欄。

#### 2. 樣本類型：組織 vs. 細胞株

| 類型 | circRNA 數量 | DE 結果 | 適合問題 |
|------|-------------|---------|---------|
| **組織**（tumor vs. normal）| 較多（整體轉錄組豐富）| 較多顯著 circRNA | 腫瘤生物標記、臨床相關性 |
| **細胞株**（drug treatment）| 較少（基因組背景單純）| 少量顯著（藥物效果有限）| 機制研究、pathway 分析 |
| **細胞株**（KO/OE）| 視目標基因而定 | 通常較集中 | 單一分子機制 |

**注意**：細胞株資料的 circRNA biomarker 臨床轉化性低，適合做機制驗證而非 biomarker 篩選。

#### 3. RNA-Seq 方式

| 方式 | circRNA 偵測 | 說明 |
|------|-------------|------|
| **Total RNA（rRNA-depleted）**（建議）| ✅ 最佳 | 保留 circRNA；rRNA-depleted 維持全轉錄組代表性 |
| **Total RNA（RNase R enriched）**| ✅✅ 最高靈敏度 | 消化線性 RNA，富集 circRNA；用於 benchmark ground truth（GSE55872）|
| poly-A selection | ❌ 不建議 | circRNA 無 poly-A tail → 幾乎偵測不到；若數據集用 poly-A，circRNA 數量會極少 |
| strand-specific | ✅ 加分 | 可推斷 circRNA strand；非必要 |

**判斷方法**：GEO → Library Selection 欄（`POLY_A`/`cDNA`/`other`）；Supplementary 找 library prep protocol。

#### 4. 定序深度（sequencing depth）

| 深度 | 說明 |
|------|------|
| **≥ 80M reads/sample**（建議）| circRNA 偵測靈敏度高；低豐度 circRNA 也能偵測 |
| 50–80M | 可用，主流研究水準 |
| < 30M | circRNA 偵測靈敏度明顯下降；少數 high-abundance circRNA 可用，但 DE 功率不足 |

**查詢方式**：SRA run info 的 `spots` 欄（paired-end: spots × 2 = total reads）。

#### 5. 樣本數與統計功率

| 組別樣本數 | DE 方法建議 | 預期顯著數 | 說明 |
|-----------|------------|-----------|------|
| n ≥ 5 vs. 5 | edgeR / DESeq2 | FDR 校正可用 | BH 校正有效；padj < 0.05 有意義 |
| n = 3 vs. 3（小樣本）| edgeR + nominal p | 需用 pvalue（未校正）| BH 校正後幾乎全不顯著（min padj ≈ 0.4）；論文需說明限制 |
| **配對設計**（tumor+normal same patient）| edgeR paired / limma blocking | 功率最高 | 消除個體差異；GSE133998 可加入配對項 |

**GSE133998 特別說明**：H36–H42 為同一患者的 cancer tissue 與 adjacent normal，是**嚴格配對設計**。已測試 `design = ~patient + condition`：limma 改善明顯（51 padj<0.05），但 edgeR 反而變差（FSJ offset 已吸收個體差異 + patient dummy 消耗 5 df）。最終決定維持 unpaired design，主要以 edgeR_ciriquant 為主方法。若要啟用配對設計，在 `metadata/GSE133998/sample_groups.csv` 加入 `patient_id` 欄即可（`analysis.R` 已支援自動偵測）。

### SRP156355 樣本清單（6 對 T/N，已生成 condition.csv）

| Patient | Normal SRR | Tumor SRR |
|---------|-----------|-----------|
| P138 | SRR7645073 (138N) | SRR7645074 (138T) |
| P148 | SRR7645077 (148N) | SRR7645078 (148T) |
| P204 | SRR7645071 (204N) | SRR7645072 (204T) |
| P272 | SRR7645075 (272N) | SRR7645076 (272T) |
| P690 | SRR7645087 (690N) | SRR7645088 (690T) |
| P1123 | SRR7645079 (1123N) | SRR7645080 (1123T) |

Condition CSV 位置：`/mnt/c/Users/User/Desktop/SRP156355_condition.csv`（含 `patient_id` 欄，可啟用配對設計）。
排除的樣本：DCIS（66D/299D/712D/803D/1102D/1151D）和 APN（719APN/768APN/87APN/91311APN/93277APN）留待未來三組比較分析。

---

### 選擇新資料集的 Checklist

在 GEO / SRA BioProject 找新 dataset 時確認以下條件（✅ = 必要，⚠️ = 建議）：

- ✅ Read length ≥ 100bp（PE）
- ✅ RNA-Seq library = Total RNA（rRNA-depleted or RNase R）
- ✅ 每組 ≥ 3 replicates（建議 ≥ 5）
- ✅ 有明確的 case vs. control（tumor/normal、treatment/vehicle）
- ⚠️ 定序深度 ≥ 50M reads/sample
- ⚠️ 組織樣本優先（而非細胞株），若要 biomarker 研究尤其重要
- ⚠️ hg19 或 hg38 人類基因組（本 pipeline 以 hg19 為主）
- ⚠️ 配對設計（paired tumor/normal）功率最高

**排除條件**：
- poly-A selection → circRNA 接近零
- single-end reads → BSJ 偵測困難
- < 30bp reads → 完全不可用
