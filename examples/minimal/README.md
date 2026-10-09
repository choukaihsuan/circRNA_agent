# Minimal reproducible example

A 1-second, dependency-free check that the **dual-tool consensus** step behaves as documented.
The inputs are hand-made toy data (not real sequencing data); `input/SCENARIOS.tsv` says what each
CIRIquant row is meant to test.

```bash
bash examples/minimal/run.sh
```

| Scenario | CIRIquant | DCC | Result with `--min-tools 2 --slop 10 --min-bsj 2` |
|----------|-----------|-----|-----|
| A  chr1:1000-2000 | BSJ 12 | 10 | **kept** (exact match) |
| B  chr1:5000-6000 | BSJ 9  | 8 (start +4 bp) | **kept** (within slop) |
| C  chr2:3000-4000 | BSJ 7  | 6 (start +25 bp) | dropped (outside slop) |
| D  chr3:800-1500  | BSJ 20 | – | dropped (one tool only) |
| E  chr4:100-900   | BSJ 1  | – | dropped (`BSJ < min-bsj`) |
| F  chr5:2000-2600 | BSJ 3, FSJ 1 | – | removed by pseudo-circRNA QC (`BSJ<5` and `BSJ/FSJ>1.0`) |
| G  chr6:7000-8000 | BSJ 30, FSJ 5 | 25 | **kept** (high BSJ is exempt from the ratio QC) |
| H  chr7:4000-4800 | BSJ 6 | 5 | **kept** |
| –  chr9:100-500   | – | 11 | dropped (DCC only) |

Expected output is in `expected/` (`high_confidence.bed`, `consensus_summary.tsv`); the script diffs against it,
and `tests/test_minimal_example.py` does the same in CI.

This covers one step of the pipeline. A full run (download → QC → CIRIquant/STAR/DCC → DE → report) needs the
conda environment or container, a reference genome with indices, and several hours; see the top-level README.
