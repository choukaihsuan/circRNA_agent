# FSJ floor sweep on the RNase R datasets
- `fsj_floor_sweep_paired_<GSE>.tsv` (GSE97239, GSE136569): produced by `~/audit_tmp/g3p.R` = the existing `g3.R` with the design switched to `~patient+condition` when `patient_id` exists (the production design of these two datasets). Baseline rows reproduce the published n_sig (140, 199).
- A first run with the unpaired `g3.R` was discarded: its baseline gave 488 / 780 significant, which does not match the published 140 / 199 (those two datasets are paired in the pipeline). Its outputs are not in this directory; copies named `fsj_floor_sweep_<GSE>.tsv` for these two datasets still sit in the server's `~/audit_tmp/` and should not be used.
- `fsj_floor_sweep_GSE113230.tsv`, `fsj_floor_sweep_GSE133998.tsv`: comparison sweeps produced by `g3.R` in another session (unpaired = their production design; baselines reproduce 240, 79).
- `retention_summary.tsv` is produced by `summarise_retention.py`.
