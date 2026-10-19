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

## Incident: `load_project_config()` leaked one project's real metadata into another

**What happened.** Testing the validation above against a brand-new,
never-before-seen accession (`GSE999999`, chosen specifically because it
should not exist) did not fail fast as expected — it downloaded real public
sequencing data belonging to a *different*, real, already-published project
(`GSE97239`) into a directory labeled `GSE999999`. A second test (a made-up
`SRR12345678`, intended to be a syntactically-valid-but-nonexistent run
accession) turned out to coincidentally be a real SRA accession and also
began downloading and processing real data before being caught. Both runs
were killed (`kill -9`) within roughly a minute of starting; all resulting
files, config snapshots, metadata, queue rows, and logs were deleted
afterward. Confirmed against all 20 published datasets' own
`config/projects/*.yaml`: every one of them correctly points at its own
metadata, before and after this incident — **this was not a published-result
problem**, only a problem for ids typed into the web form during testing.

**Root cause** (`scripts/web_ui.py`, pre-fix `load_project_config()`): when a
project's own `config/projects/{id}.yaml` snapshot didn't exist yet, the
function fell through to `load_config()` — the single, mutable, shared
`config.yaml` that reflects whichever project was most recently submitted by
*anyone*, through *any* route. `/run_gse` then only overwrote a handful of
fields on that inherited config (`project_id`, `source_type`, `threads`, and
the raw/trimmed/results directory names via now-removed
`_update_paths_for_project()`'s string substitution) — critically, it never
reset `metadata` or `groups`. So a brand-new id inherited the *previous*
project's real `metadata` path, that path existed (because the previous
project was real and already downloaded), and `/run_gse`'s "does metadata
already exist → treat this as a re-run, skip `agent.py`'s download step and
go straight to Snakemake" check fired on someone else's data.

**The fix** (this change): `load_project_config()` is gone, replaced by two
choke points:

- `build_project_config(project_id, source_type, threads)` — used by all
  three submit routes (`/run_gse`, `/run_manual`, `/run_local`) instead of
  each hand-rolling its own load-then-patch sequence (that per-route
  duplication is exactly how `/run_gse` ended up missing the `metadata`/
  `groups` reset that the other two routes happened to include). When a
  project has no snapshot of its own yet, it starts from `config/default.yaml`
  — a clean, committed template with **no per-project identity fields at
  all** — never from the shared `config.yaml`. In both the new-project and
  existing-project cases, identity/path fields (`raw_dir`, `trimmed_dir`,
  `results_dir`, `metadata`, `groups`, `download.sra_cache_dir`/`tmp_dir`)
  are then *unconditionally* (re)derived from `project_id` and `data_root`
  via plain f-string construction — never string-substituted from whatever
  the loaded config happened to already contain, which is what let a stale
  path "look" correct whenever the old id happened to already be embedded in
  it.
- `_load_existing_project_config_or_none()` / `_load_existing_project_config_or_404()`
  — the read-side counterpart, used by `/api/progress`, `/report`,
  `/download`, `/qc`. These never fall through to any other config either;
  a missing snapshot means "nothing to show" (404, or a graceful no-op for
  `/api/progress`'s best-effort polling), not "show whatever's currently
  active globally." This closes a narrower, lower-probability sibling of the
  same bug: if a job's own snapshot were ever missing, these routes could
  previously have served a *different* project's `report.html`/QC report to
  someone holding a valid `job_id` for the missing-snapshot one.

Also fixed as a direct consequence: previously, *every* field a route didn't
explicitly reset (`consensus.tools`, `circbase_file`, `de.method`,
`de.fdr_cutoff`, `genome.*`, etc.) silently inherited the previous project's
value too, on all three submit routes — not just `metadata`/`groups`, and
not just `/run_gse`. `build_project_config()`'s clean-template-or-own-
snapshot rule closes this for every field, not only the two that caused the
visible incident.

**Deferred, not fixed, in this round**: replacing "does `metadata/{id}/
library_info.csv` exist" as the re-run-vs-fresh-download signal in `/run_gse`
with an explicit status field recorded by the project itself. With
`build_project_config()`'s fix, that path is now *deterministically* always
`metadata/{project_id}/library_info.csv` — it can no longer resolve to a
different, foreign project's real metadata the way it could before — so the
check is safe again, just not as architecturally clean as an explicit status
field would be. Implementing that properly requires `agent.py` (a separate
process from `web_ui.py`, invoked as a subprocess) to write a completion
marker back into the project's own config after `prepare_metadata.py`
succeeds, which is a larger change than this round's scope. Tracked as a
follow-up, not silently dropped.

**Not investigated further in this round, flagged for awareness**: whether a
published dataset's config could have been *transiently* mis-pointed at some
earlier point and then overwritten by a later, correct submission before
anyone looked — the current state and file mtimes are clean and consistent
with "correct from the start," but mtimes alone cannot fully rule out a
since-overwritten transient error. No evidence of this was found (no
CLAUDE.md record, commit message, or leftover metadata directory suggests it
happened), but absence of evidence isn't proof it didn't.
