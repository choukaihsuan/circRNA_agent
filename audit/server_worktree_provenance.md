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

## Unreferenced server files (kept in place; not moved or deleted)

Verified unreferenced: none of these paths or basenames is named in workflow/, scripts/ or config.yaml.
Recorded so their identity is known without moving them.

| server path | sha256 |
|---|---|
| scripts/download.smk | 30a401415e8bf456af426a77463402e6d724e9c7783fa684f5b2f6e3f367971c |
| scripts/index.html | cb5f91613b5801913321f67ff3c82d71f4ac80e433c9ae4735d3ef3354420fe7 |
| scripts/login.html | c3555770f8404e6da158837e3b3eeff80339e6689d77c52b8339308fd24534f4 |
| scripts/mnt/c/Users/User/develop/circRNA_agent/scripts/templates/queue.html | f10e77e779d348afb74548e3aea4b1cb68f844e23f52de02fbc1f6d3dde06468 |
| scripts/mnt/c/Users/User/develop/circRNA_agent/scripts/templates/status.html | c9328a49a1f1aefcf8de82d9ec30950d58a623fe742476ca1d07e2b78c717953 |
| scripts/mnt/c/Users/User/develop/circRNA_agent/scripts/web_ui.py | e512ba1a6c00fbe924009dd33322901055e325da228df23de1a7ff1513bd5b59 |
| scripts/qc.smk | ea3cd285b373cb13584a1cd1e4e1d46ac7f0bbe8d496b6b749df4e3f59dae25c |
| scripts/queue.html | 982b78924abebabb8b0ec39e32d14f172f1b08f9e85fcd64cf31c5c7886a683b |
| scripts/sensitivity_GSE58135.py | 32c046d07edf83183a016a8cfbf2dba4b51cce63b6ea0ae6856abd8c29387db7 |
| scripts/status.html | 375f9745b2b2a7a57ca214a24d901893c183b18a3a785e58032aca85441a7f0a |
| scripts/templates/web_ui.py | 86cd94569fbc8371dd180e0a223a5356a9ea1584d6b5768b8ded7931e53a92f6 |
