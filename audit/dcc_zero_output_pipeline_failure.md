# DCC zero-output hard failure (robustness defect, not fixed)

Tracking: 代辦清單 第 13 項（輸入驗證／規模閘門）。Not fixed in this round by instruction.
Manuscript: candidate limitation for §4.10 (Discussion).

## Observed

- Test: SRR444655 (SRP011042), first 2,000,000 read pairs (fastq-dump -X), STAR paired/mate1/mate2 chimeric junctions, DCC with `-Nr 2 1`.
- DCC log (`DCC-*.log`): `No circRNA passed the expression threshold filtering.`
- DCC exited non-zero. Snakemake treated the `dcc` rule as failed and stopped the run (no `--keep-going` fallback for this path).
- Date: 2026-10-05 (B3 threads test, scratch directory only).

## Why this is a defect, not a data result

- Zero circRNA passing the count threshold is a legitimate outcome for a shallow sample. The correct behaviour is for DCC to report zero circRNA and let `consensus_filter` handle it.
- `consensus_filter.parse_dcc` already handles an empty or header-only `CircCoordinates` (0 loci), and `--adaptive` degrades the consensus to the remaining tool. So the downstream logic already supports the outcome.
- The failure is therefore the non-zero exit code of DCC itself, which Snakemake reads as a hard error.

## Same root cause seen in production

- GSE192410 (9 samples) and GSE192849 (6 samples) have header-only `DCC/CircCoordinates` files. These were produced by hand (a header-only file) to bypass DCC's O(n²) duplicate-marking stall on high-junction samples. Those workarounds exist because DCC's failure modes are not handled by the pipeline.
- Both cases are the same defect seen from two sides: a legitimate DCC result of "nothing to report" or "too slow" has no first-class path through the pipeline. The current path is a non-zero exit (here) or a manual header-only file (there).

## Required behaviour (for the future fix)

- DCC returning zero circRNA is not an error. The rule should succeed with an empty result and let `consensus_filter` apply its adaptive fallback.
- Input validation (代辦 第 13 項): check junction counts before running DCC and report a scale warning instead of letting DCC fail or stall.
- Any manual header-only workaround must be replaced by the same code path, so the outcome is recorded rather than hand-edited.

## Limitations to disclose

- The 2M-read subsample yields no DCC circRNA, so it is not a usable benchmark input for any DCC-dependent comparison.
- The GSE192410 and GSE192849 DCC files were hand-edited to header-only; this must appear in their Methods as a deviation.
