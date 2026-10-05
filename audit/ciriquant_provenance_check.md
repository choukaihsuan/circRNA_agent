# CIRIquant provenance trigger: explained (read-only check, 2026-10-05)

Question: default-trigger dry-runs showed `ciriquant` jobs with "Code has changed; Params have changed"
for GSE97239, GSE192410 and GSE192849. Is this an unexplained hash change?

Answer: no. The trigger is the rule change in commit 2f27d7b (2026-10-02), which added the
`anchor_arg` and `no_gene` params to `rule ciriquant`. The recorded metadata predates that commit.

Evidence (`.snakemake/metadata`, read-only; the summary script was run from scratch and is not committed):

- 197 ciriquant metadata records in the repo tree. 196 record 2 params (`outdir`, `sample`), which is
  the schema before 2f27d7b. 1 record records 4 params; it is a rerun from 2026-10-03, after the commit.
- Every published dataset's ciriquant metadata shows 2 params. Their recorded start times fall in
  May–July 2026 (GSE97239 about 2026-07-29; GSE192410 and GSE192849 about 2026-07-29 to 2026-08-02).
- The current rule records 4 params. The shell text is unchanged for default values:
  `-a` is emitted only when anchor differs from 5, and `--no-gene` only when `no_gene` is true, so the
  default CIRIquant command is byte-identical to the one that produced the published GTFs.

Consequences:

1. Under default triggers, `ciriquant` is planned for every published dataset (not only the three named),
   because the params changed on 2026-10-02. This is a provenance difference, not a change in output.
2. The regeneration in item 8 (terminal, no DAG) does not run CIRIquant, so the provenance difference
   does not affect the outputs. It remains a reason never to run Snakemake on published directories
   without `--rerun-triggers mtime`.
3. A future decision could refresh the metadata without rerunning (for example by a dedicated record step).
   That is not done here; the metadata is left untouched.

Excluded from this count: 8 ciriquant records under the benchmark directory (no `_results` path match).
