# Platform input validation — design notes

## SRR/accession format checks live at the input layer, not in Snakemake

`scripts/web_ui.py`'s `/run_manual` route validates every `srr_id` (from both
the uploaded-CSV path and the manual-entry form fields) against `SRR_RE`
(`^(SRR|ERR|DRR)\d+$`) **before** any row is written to `metadata/library_info.csv`.
One invalid row aborts the entire submission (400) and reports which row —
nothing is written and nothing is silently dropped.

`/run_gse` separately validates the dataset accession against `ACCESSION_RE`
(`^(GSE|PRJNA|SRP|SRX|SRS|SRR|ERP|DRP)\d+$`, matching exactly what
`download_geo.py`'s `_SRA_PREFIXES` + `GSE` actually dispatch on) before it is
used as a project id or passed to `agent.py --gse`.

**This validation is deliberately NOT relied on being redundant with
`circrna.smk`'s `wildcard_constraints: srr = r"[A-Z]+\d+"`.** That constraint
is real and, as far as we understand Snakemake's documented behavior, applies
workflow-wide once `circrna.smk` is included into the top-level `Snakefile` —
but this has not been verified by actually running Snakemake (not even
`--dry-run`) against a deliberately malformed `srr` value, because the task
that added this input-layer validation was explicitly scoped to avoid running
Snakemake at all. **Do not treat the wildcard constraint as the thing standing
between a malformed id and a scheduled job.** If the input-layer checks in
`web_ui.py` are ever refactored away on the assumption that "Snakemake already
blocks this," that assumption has never actually been exercised end-to-end —
verify it first (e.g. a `--dry-run` with a crafted `srr` value), or keep both
layers.

## `source_type` and the accession/local split

Every project config now carries a `source_type` key, set to `"accession"` in
`/run_gse` and `"local"` in `/run_manual` and `/run_local`. The one place in
the app where a project's own id (as opposed to a per-sample SRR id) is passed
as an argument to something that triggers a network download —
`/run_gse`'s fallback to `agent.py --gse {gse_id}` for a fresh dataset —
asserts `cfg.get("source_type") == "accession"` immediately before building
that command. This doesn't change behavior for any currently-reachable code
path (that branch is only ever reached from `/run_gse`, which always sets
`source_type="accession"` a few lines earlier) — it exists so that if a future
change ever lets a `CUSTOM`/`LOCAL` (manual-upload) project id reach that
branch, the assertion fails loudly instead of quietly building a download
command from an id that was never validated as a real accession.

## `PROJECT_ID_RE` vs `ACCESSION_RE`/`SRR_RE`: why there are three regexes, not one

- `PROJECT_ID_RE` (`^[A-Za-z0-9][A-Za-z0-9_-]{0,39}$`) is the generic,
  pre-existing path-safety check (`_validate_project_id`). It has to admit
  `CUSTOM`, `LOCAL`, and any other id a user types for a manual/local project,
  so it only guarantees "safe as a path component" — not "looks like a real
  accession."
- `ACCESSION_RE` is strictly narrower and is only applied in `/run_gse`, the
  one route that treats its id as a dataset accession to fetch over the
  network.
- `SRR_RE` is narrower still (run-level accessions only, not study/project/
  experiment-level ones) and is only applied to individual `srr_id` values in
  `/run_manual`.

Applying `ACCESSION_RE` to `/run_manual`'s `project_id` field, or vice versa,
would be wrong: `/run_manual`'s `project_id` is a free-form project label
(`CUSTOM` by default), never a value passed to a download.
