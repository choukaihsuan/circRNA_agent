# Server working tree: provenance against local git history (2026-10-05, read-only)

Method: sha256 of every .py/.R/.smk/Snakefile/.html/.sh/.yaml under scripts/ and workflow/ on the server; for each path, every local commit that touched it was checked with `git show <sha>:<path>` and the blob hash compared. Nothing was pulled, reset, checked out or written on either side.

## Summary

- files examined: 43
- matched exactly one commit: 32; of these, equal to the local working tree now: 28; differs from local working tree (local has later edits not on the server): 4
- matched no commit at all: 11 (see below)

## Files that differ from the local working tree (server = an earlier commit)

| path | server matches commit | date |
|---|---|---|
| scripts/consensus_filter.py | 348bd06 | 2026-10-05 |
| scripts/generate_report.py | 348bd06 | 2026-10-05 |
| scripts/predict_interactions.py | 4644cf0 | 2026-07-07 |
| scripts/rank_biomarkers.py | 348bd06 | 2026-10-05 |

These are the four files changed by the run_params work (uncommitted on the server by design). The server holds the B2 state (commit 348bd06) for them.

## Files that match no commit

| path | note |
|---|---|
| scripts/download.smk | not a repository file; differs from workflow/rules/download.smk; origin unknown |
| scripts/index.html | not in the repository (templates/index.html is) |
| scripts/login.html | not in the repository |
| scripts/mnt/c/Users/User/develop/circRNA_agent/scripts/templates/queue.html | nested path scripts/mnt/c/Users/User/develop/circRNA_agent/... : a rsync that used an absolute path (directory created 2026-05-14) |
| scripts/mnt/c/Users/User/develop/circRNA_agent/scripts/templates/status.html | nested path scripts/mnt/c/Users/User/develop/circRNA_agent/... : a rsync that used an absolute path (directory created 2026-05-14) |
| scripts/mnt/c/Users/User/develop/circRNA_agent/scripts/web_ui.py | nested path scripts/mnt/c/Users/User/develop/circRNA_agent/... : a rsync that used an absolute path (directory created 2026-05-14) |
| scripts/qc.smk | identical to workflow/rules/qc.smk: a misplaced copy |
| scripts/queue.html | not in the repository |
| scripts/sensitivity_GSE58135.py | not in the repository |
| scripts/status.html | not in the repository |
| scripts/templates/web_ui.py | not in the repository |

## Files that match one commit and equal the local working tree (28)

These are current on both sides. Listed by commit:

- 1483a79 2026-06-20: scripts/assign_isoforms.py
- 20c206e 2026-05-26: scripts/annotate_circbase.py
- 2c4b31f 2026-05-14: scripts/utils.py
- 2d19cd7 2026-05-27: scripts/merge_counts.py
- 2fa4253 2026-07-02: scripts/gen_pipeline_diagram.py
- 32df542 2026-10-04: scripts/loo_stability.R
- 348bd06 2026-10-05: scripts/analysis.R, scripts/validate_config.py, workflow/Snakefile, workflow/rules/circrna.smk
- 42453d4 2026-08-02: scripts/templates/login.html
- 4644cf0 2026-07-07: scripts/notify.py
- 478f2f1 2026-10-04: scripts/templates/index.html, scripts/web_ui.py
- 49b122d 2026-10-01: scripts/run_permutation.py
- 4c5d41a 2026-10-05: workflow/rules/de.smk
- 4dfa28a 2026-07-15: scripts/cross_dataset_biomarkers.py, scripts/prepare_metadata.py, scripts/validate_condition_labels.py
- 798f696 2026-10-03: scripts/templates/cross_dataset.html
- 79926ee 2026-07-01: workflow/rules/qc.smk
- 93d8739 2026-06-10: scripts/agent.py, scripts/download_geo.py
- 9c493de 2026-07-21: scripts/gen_pipeline_diagram2.py, workflow/rules/download.smk
- ac8b060 2026-08-10: scripts/templates/status.html
- e003930 2026-06-07: scripts/isoform_switching.R
- f03ef98 2026-08-22: scripts/templates/queue.html

## What this means

- The server is not a single commit. Its tracked files correspond to different commits (2026-05-14 to 2026-10-05), consistent with rsync-based deployment. No single commit reproduces the server tree.
- For the published results, the relevant question is whether each step's code matches an identifiable commit. For workflow/ and most scripts it does; none of the 11 unmatched files is referenced anywhere in the repository (checked by grep on workflow/, scripts/ and config.yaml). The Snakefile includes only workflow/rules/*.smk.
- Nothing here has been synchronised. A pull or reset would overwrite the server's working files, which are the files that produced the published results; that is why this comparison was done first.
