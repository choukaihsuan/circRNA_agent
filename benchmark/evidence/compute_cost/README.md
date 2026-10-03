# compute_cost provenance (captured 2026-10-03)

- `time_logs/` — raw `/usr/bin/time -v` logs copied read-only from the HPC server `~/circRNA_agent/logs/bench/` (mtimes 2026-06-05 … 2026-10-01), plus `compute_cost.log` from the 2026-10-01 11:08 run.
- `tsv/compute_cost_after_fix.tsv` — byte copy of the server's `GSE113230_results/benchmark/compute_cost.tsv` (mtime 2026-10-01 11:08:12), produced by the fixed script.
- `tsv/compute_cost_before_fix.REPRODUCED.tsv` — NOT the original. The original pre-fix tsv was overwritten on the server and was never committed or backed up. This file is the output of `compute_cost.py` at git HEAD (2f27d7b) run with the pre-fix Snakefile arguments (no `--our-star-mate1/2-log`) on the logs above. Its values match the numbers recorded for the pre-fix tsv (e.g. circDEX Total 811.1, nf-core 1084.2).
- `compute_cost_before_after.tsv` — cell-level diff with reasons.

The fixed `compute_cost.py` / `benchmark/Snakefile` were uncommitted when this evidence was captured.
