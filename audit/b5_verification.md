# B5 verification record (2026-10-05)

## 1. Rerun set before and after the run_params change (dry-run, read-only)

Both trees use the same metadata. "before" is the server working tree (scripts equal to commit 348bd06);
"after" is a scratch copy with the four instrumented scripts and run_params.py. Each run uses
`snakemake -n -r --configfile config/projects/{GSE}.yaml`, with and without `--rerun-triggers mtime`.

| dataset | trigger | before rules | after rules | added | removed | count changes | status before / after |
|---|---|---|---|---|---|---|---|
| GSE113230 | default | 20 | 20 | - | - | - | ok / ok |
| GSE113230 | mtime | 20 | 20 | - | - | - | ok / ok |
| GSE97239 | default | 20 | 20 | - | - | - | ok / ok |
| GSE97239 | mtime | 0 | 0 | - | - | - | ok / ok |
| GSE121842 | default | 20 | 20 | - | - | - | ok / ok |
| GSE121842 | mtime | 20 | 20 | - | - | - | ok / ok |
| GSE133998 | default | 20 | 20 | - | - | - | ProtectedOutput / ProtectedOutput |
| GSE133998 | mtime | 20 | 20 | - | - | - | ProtectedOutput / ProtectedOutput |

Three datasets (GSE113230, GSE97239, GSE121842) plan completely, with identical rule sets before and after.
GSE133998 aborts planning at the write-protected CIRIquant GTF in both trees (pre-existing; kept as is), so
its comparison is partial. The earlier stop flag in the comparison script wrongly counted that pre-existing
error as a change; the corrected rule compares added rules, removed rules, count changes and status changes.

## 2. Small-subset predict_interactions test (scratch)

Input: the first five rows of the GSE121842 DE table (scratch copy). Output: `run_params/predict_interactions.yaml`
written by the step itself. Result: provenance `executed`; `external_queries` holds `queried_at`, `n_circrnas: 5`,
CircInteractome and ENCORI hit counts for miRNA and RBP, and both source names. The record also holds the
effective parameters (`top_n: 50`), four input fingerprints, and the code commit, read with `cwd=` because the server's git does not support `-C`. That commit is the
server's own history (see audit/server_worktree_provenance.md).

## 3. Full end-to-end run (GSE121842, scratch)

consensus_filter (one sample), rank_biomarkers_v2 (three methods in parallel), generate_report (mock runner):
one record file each. The three rank records have identical start times, confirming parallel writes did
not collide. Consolidation merged the per-step files into `run_params.yaml` (provenance `executed`).
predict_interactions was not part of the full run; its record is verified by section 2.

## 4. Known limits

- The server git history is 246 commits behind the local repository (see audit/server_worktree_provenance.md).
  Records written on the server therefore carry the server's commit SHA (2c4b31f) with `git_dirty: true`;
  they do not identify the local code. No backfill has been done.
- CIRIquant, STAR, DCC, fastp, merge_counts, de_analysis, isoform_switching and assign_isoforms have no
  records yet.
