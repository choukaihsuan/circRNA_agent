# CIRIquant `-l/--library-type` 的作用範圍（原始碼查證）

來源：server `ciriquant` conda env 的 `site-packages/CIRIquant/`（`version.py` = 1.1.3；19 組 182 個輸出 GTF 標頭皆為 `##version: 1.1.3`）。只讀查證，未執行 CIRIquant。

## 參數怎麼流動（`grep -n lib_type *.py` 的全部出現處）
| 位置 | 內容 |
|---|---|
| `main.py:37-39` | `-l/--library-type`，預設 0（0 unstranded，1 read1 同義股，2 read1 反義股） |
| `main.py:86-93` | 轉 int 並驗證 ∈ {0,1,2}；`141-146` 僅用於 log 顯示 library 名稱 |
| `main.py:197` | **唯一傳遞點**：`circ.proc(..., anchor, lib_type, args.no_fsj, args.bsj_read_file)` |
| `circ.py:628` | `proc(..., lib_type, ...)` |
| `circ.py:661` | `cand_bsj = proc_denovo_bam(denovo_bam, thread, circ_info, anchor, lib_type)` |
| `circ.py:296, 317` | `proc_denovo_bam` → 把 `lib_type` 交給 `denovo_worker` |
| `circ.py:343-390` | `denovo_worker`：**實際使用處**。`lib_type==0`（375）不做任何過濾；`==1`（377-383）與 `==2`（384-390）依 read1/read2 與 host circRNA 的 `parser.strand` 丟棄 strand 不符的 read |
| `circ.py:358` | 程式內註解：`#TODO: add support for stranded specific RNA_seq libraries` |

## 沒有接收 `lib_type` 的步驟
- **FSJ**：`circ.py:662` 呼叫 `proc_genome_bam(hisat_bam, thread, circ_info, cand_bsj, anchor, circ_dir, is_no_fsj)`（定義 `406`；`genome_worker` `484`）——參數列沒有 `lib_type`。`circ.py` 第 400–627 行（FSJ/genome-BAM 區段）內**沒有任何 `strand` 字樣**（awk 掃描為空）。
- **CIRI2 偵測、HISAT2/BWA 比對、基因丰度（StringTie）**：`pipeline.py`、`main.py` 其他段落沒有任何 `lib_type` 引用。

## 結論（限於程式碼可證的部分）
1. `-l 0`（我們所有資料集的設定）對 stranded library 的**唯一**影響，是 **BSJ 候選 read 的 strand 過濾被略過**（`circ.py:375`）：host-strand 不符的 BSJ read 不會被丟棄，可能計入 BSJ。
2. **FSJ 計數不受 `-l` 影響**（`proc_genome_bam` 不接收、`circ.py` FSJ 區段不含 strand 邏輯）。
3. 影響方向只能是 BSJ 計數偏高（多收了 strand 不符的 read），幅度大小**本查證無法給出**：沒有實測（重跑 CIRIquant 11–12 小時/樣本），也不引用外部文獻數字。
