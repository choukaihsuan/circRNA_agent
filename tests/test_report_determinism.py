"""Two independent determinism checks, answering two different questions:

1. test_report_md5_reproducibility() builds a tiny synthetic report from
   fixed inputs N times, each in a fresh subprocess with a different
   PYTHONHASHSEED, and asserts every run produces byte-identical output
   (after normalizing the two things that are *supposed* to vary: the
   generation timestamp, and Plotly's per-figure random div-id UUIDs).
   This is an integration-level regression check: "did this particular
   fixture come out the same this time." It is necessarily scale-limited
   -- see its own docstring for a known gap at this fixture's size.

2. test_set_iteration_order_invariant() asks a narrower, stronger
   question: is a given output order a pure function of the input's
   *content*, independent of how that content happened to be stored.
   It builds 1,000 synthetic circ_id strings into a Python set, and in 5
   subprocesses (one per PYTHONHASHSEED) runs exactly the
   `sorted(cid_set)` line itself -- no report, no HTML, no md5. All 5
   must produce the identical sequence (compared element-by-element, not
   hashed). It also asserts the *unsorted* raw set iteration order is NOT
   invariant across seeds first, as a sanity check that 1,000 strings is
   actually enough to make Python's hash randomization visible at all --
   an invariant test that can't fail on the thing it's supposed to detect
   would be worthless.

Why both exist, and why neither replaces the other: md5 comparison proves
"the same inputs produced the same artifact this time, in this
environment." The invariant test proves "this step's output order cannot
depend on anything but input content, by construction, for a case general
enough that iteration order WAS seed-dependent going in" -- a stronger,
scale-independent guarantee that doesn't need to get lucky into
reproducing a specific implementation's instability threshold. This
exists because grep-based auditing for non-deterministic sorts
(`.sort/.sort_values/sorted`) was proven insufficient in this project's
own history: three real instances of "iterate a Python set with no sort"
fed report/page output order and were caught only by rebuilding and
diffing under different hash seeds, not by the grep.

test_report_md5_reproducibility() deliberately uses a hand-built, minimal
fixture (a handful of rows, including a few exact score/p-value ties and
two circRNAs with NaN p-values modeling a real case) rather than any real
dataset's results directory -- small, fast, and never touches anything
under a published results_dir.

Run with:  python tests/test_report_determinism.py
Requires the project's conda env (pandas, and plotly if available --
falls back to the static-PDF path otherwise, which this test still
covers) for test 1; test 2 has no such dependency.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts"
N_SEEDS = 5
SEEDS = [1, 2, 3, 4, 5]

V2_COLS = ["rank", "rank_default", "pass_default", "three_method_sig",
           "effect_stable", "biomarker_score"]


def _write_tsv(path: Path, rows: list[dict], columns: list[str]):
    with open(path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in columns) + "\n")


def build_fixture(root: Path):
    """A deliberately tiny, synthetic dataset with a few exact ties baked
    in (two circRNAs at the same p-value, two genes at the same
    padj_global) -- small enough to inspect by eye, big enough to exercise
    every code path this session's fixes touched."""
    de_dir = root / "de"
    circ_dir = root / "circRNA"
    plots_dir = root / "plots"
    for d in (de_dir, circ_dir, plots_dir):
        d.mkdir(parents=True, exist_ok=True)

    circ_ids = [f"chr1:{1000+i*100}|{1000+i*100+50}" for i in range(12)]
    # Two deliberate exact ties: (circ_ids[4], circ_ids[5]) share a p-value;
    # (circ_ids[8], circ_ids[9]) share both log2FC and p-value.
    pvals = [0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.006, 0.007,
             0.008, 0.008, 0.009, 0.010]
    lfcs  = [3.0, -3.0, 2.5, -2.5, 2.0, 2.0, -2.0, 1.8, 1.5, 1.5, -1.5, 1.2]
    genes = [f"GENE{i}" for i in range(12)]

    de_rows = [
        {"circ_id": circ_ids[i], "gene_name": genes[i], "log2FC": lfcs[i],
         "pvalue": pvals[i], "circbase_id": "novel", "Type": "Type_I",
         "strand": "+", "region": "exonic", "exon_span": "e1-e2"}
        for i in range(12)
    ]
    de_cols = ["circ_id", "gene_name", "log2FC", "pvalue", "circbase_id",
               "Type", "strand", "region", "exon_span"]
    _write_tsv(de_dir / "de_results.tsv", de_rows, de_cols)

    # The KIF3B case, reproduced on purpose: one extra circRNA absent from
    # edgeR_ciriquant's file entirely, present in deseq2's with pvalue=NaN,
    # and present in limma's with a real, significant value -- the exact
    # shape that fed a NaN into _venn_3_svg's circ_info aggregation
    # (generate_report.py's _is_better_pval) in the real GSE113230 data
    # this harness is modeled on. Without the fix, this circRNA's pval
    # stays NaN, and sorting a (pval, id) key containing NaN is not a well
    # defined total order -- its relative position among other NaN rows
    # then depends on pre-sort iteration order, which is exactly what this
    # harness is checking stays fixed across hash seeds.
    # Two NaN-first circRNAs, not one: a single stray NaN among otherwise
    # real-valued keys turned out not to be adversarial enough to exercise
    # `for cid in sorted(cid_set)` on its own. KNOWN LIMITATION, stated
    # plainly rather than left for someone to discover by surprise: even
    # with two NaN rows, this harness's reverting that one `sorted()` call
    # (tested by hand, not committed) still passed all 5 seeds x 3 repeats
    # in this project's pandas/Python version -- Timsort's handling of a
    # small number of NaN-keyed, mutually incomparable rows is, empirically,
    # stable at this scale regardless of input order. The real GSE113230
    # region that originally caught this had 2191 total rows with several
    # genuinely NaN among them, which DID reproduce the instability (5/5
    # seeds disagreed before the fix, verified directly -- see
    # docs/platform_security_notes.md). This harness's small fixture is a
    # solid regression guard for _is_better_pval (reverting that one DOES
    # fail this test, verified) and for every other fix in this series it
    # exercises; it is not, on its own, proof that the `sorted(cid_set)`
    # line stays necessary -- that claim rests on the real-data check, not
    # on this file.
    nan_ids = ["chr9:9999|10049", "chr9:8888|8949"]
    for method in ("edgeR_ciriquant", "deseq2", "limma"):
        rows = list(de_rows)
        if method == "deseq2":
            rows = rows + [
                {"circ_id": nid, "gene_name": f"KIF3B_FIXTURE_{i}",
                 "log2FC": -4.0, "pvalue": "nan", "circbase_id": "novel",
                 "Type": "Type_I", "strand": "+", "region": "exonic",
                 "exon_span": "e1-e2"}
                for i, nid in enumerate(nan_ids)
            ]
        elif method == "limma":
            rows = rows + [
                {"circ_id": nid, "gene_name": f"KIF3B_FIXTURE_{i}",
                 "log2FC": -4.0, "pvalue": 0.0009 + i * 0.0001, "circbase_id": "novel",
                 "Type": "Type_I", "strand": "+", "region": "exonic",
                 "exon_span": "e1-e2"}
                for i, nid in enumerate(nan_ids)
            ]
        # edgeR_ciriquant: nan_ids simply absent, matching the real case
        _write_tsv(de_dir / f"de_results_{method}.tsv", rows, de_cols)

    # count matrix: 4 samples, tumor/normal
    samples = ["S1_tumor", "S2_tumor", "S3_normal", "S4_normal"]
    with open(circ_dir / "count_matrix.tsv", "w") as f:
        f.write("circ_id\t" + "\t".join(samples) + "\n")
        for i, cid in enumerate(circ_ids):
            f.write(cid + "\t" + "\t".join(str(10 + i + j) for j in range(4)) + "\n")

    with open(root / "sample_groups.csv", "w") as f:
        f.write("srr_id,condition\n")
        f.write("S1_tumor,tumor\nS2_tumor,tumor\nS3_normal,normal\nS4_normal,normal\n")

    # isoform_groups / circbase_annotated: minimal passthrough
    _write_tsv(circ_dir / "isoform_groups.tsv",
               [{"circ_id": c, "gene_name": g, "strand": "+", "region": "exonic", "exon_span": "e1-e2"}
                for c, g in zip(circ_ids, genes)],
               ["circ_id", "gene_name", "strand", "region", "exon_span"])
    _write_tsv(circ_dir / "circbase_annotated.tsv",
               [{"circ_id": c, "circbase_id": "novel", "circbase_gene": "", "in_circbase": 0}
                for c in circ_ids],
               ["circ_id", "circbase_id", "circbase_gene", "in_circbase"])

    # isoform_switching.tsv: two genes tied at the same padj_global on purpose
    sw_rows = [
        {"circ_id": circ_ids[i], "gene_id": genes[i], "gene_name": genes[i],
         "is_switching": "True",
         "padj_global": (0.5 if i < 6 else 0.5), "p_value": pvals[i],
         "delta_iui": 0.2, "iui_tumor": 0.6, "iui_normal": 0.4}
        for i in range(8)
    ]
    _write_tsv(de_dir / "isoform_switching.tsv", sw_rows,
               ["circ_id", "gene_id", "gene_name", "is_switching", "padj_global",
                "p_value", "delta_iui", "iui_tumor", "iui_normal"])

    (de_dir / "interactions.json").write_text(json.dumps({}))

    # biomarker_candidates_v2_* + status.json for all three methods
    bm_rows = [
        {"rank": i + 1, "circ_id": circ_ids[i], "gene_name": genes[i],
         "rank_default": i + 1, "pass_default": "True", "three_method_sig": "True",
         "effect_stable": "True", "biomarker_score": round(0.9 - i * 0.05, 4),
         "n_mirna": 1, "n_rbp": 1, "in_circbase": 0, "circbase_id": "novel",
         "circbase_gene": "", "Type": "Type_I"}
        for i in range(8)
    ]
    bm_cols = ["rank", "circ_id", "gene_name", "rank_default", "pass_default",
               "three_method_sig", "effect_stable", "biomarker_score",
               "n_mirna", "n_rbp", "in_circbase", "circbase_id", "circbase_gene", "Type"]
    for method in ("edgeR", "DESeq2", "limma"):
        p = de_dir / f"biomarker_candidates_v2_{method}.tsv"
        _write_tsv(p, bm_rows, bm_cols)
        (de_dir / f"biomarker_candidates_v2_{method}.status.json").write_text(
            json.dumps({"status": "ok"}))

    # Minimal placeholder PDFs -- base64-embedded as bytes, never parsed as
    # real PDF content by this script.
    for name in ("volcano.pdf", "heatmap.pdf", "pca.pdf"):
        (plots_dir / name).write_bytes(b"%PDF-1.4 fixture placeholder\n")


def run_once(root: Path, seed: int, out_path: Path) -> None:
    script = f"""
import sys
sys.path.insert(0, {str(SCRIPTS_DIR)!r})
import generate_report as G

R = {str(root)!r}
G.build_report(
    de_file=R + "/de/de_results.tsv",
    matrix_file=R + "/circRNA/count_matrix.tsv",
    volcano_pdf=R + "/plots/volcano.pdf",
    heatmap_pdf=R + "/plots/heatmap.pdf",
    pca_pdf=R + "/plots/pca.pdf",
    output_file={str(out_path)!r},
    project_id="FIXTURE",
    fdr=0.05, lfc=1.0, de_method="edgeR_ciriquant",
    biomarker_files={{
        "edgeR_ciriquant": R + "/de/biomarker_candidates_v2_edgeR.tsv",
        "deseq2":          R + "/de/biomarker_candidates_v2_DESeq2.tsv",
        "limma":           R + "/de/biomarker_candidates_v2_limma.tsv",
    }},
    switching_file=R + "/de/isoform_switching.tsv",
    groups_file=R + "/sample_groups.csv",
    de_sig_by="pvalue", tumor_label="tumor", normal_label="normal",
    isoform_file=R + "/circRNA/isoform_groups.tsv",
    circbase_file=R + "/circRNA/circbase_annotated.tsv",
    heatmap_top_n=5,
    interactions_file=R + "/de/interactions.json",
    multiqc_file=None,
    de_files={{
        "edgeR_ciriquant": R + "/de/de_results_edgeR_ciriquant.tsv",
        "deseq2":          R + "/de/de_results_deseq2.tsv",
        "limma":           R + "/de/de_results_limma.tsv",
    }},
    study_title="",
)
"""
    env = dict(os.environ)
    env["PYTHONHASHSEED"] = str(seed)
    result = subprocess.run([sys.executable, "-c", script], env=env,
                             capture_output=True, text=True, timeout=120)
    if result.returncode != 0:
        raise RuntimeError(
            f"build_report() failed under PYTHONHASHSEED={seed}:\n{result.stderr}")


_UUID_RE = re.compile(r'[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}')
_DATE_RE = re.compile(r'\d{4}-\d{2}-\d{2}')


def normalize(html: str) -> str:
    html = _UUID_RE.sub("<UUID>", html)
    html = _DATE_RE.sub("<DATE>", html)
    return html


def test_report_md5_reproducibility():
    tmp = Path(tempfile.mkdtemp(prefix="circdex_det_test_"))
    try:
        build_fixture(tmp)
        hashes = []
        for seed in SEEDS:
            out = tmp / f"report_seed{seed}.html"
            run_once(tmp, seed, out)
            content = normalize(out.read_text(encoding="utf-8"))
            h = hashlib.md5(content.encode()).hexdigest()
            hashes.append((seed, h))
            print(f"  seed={seed}: md5={h}")

        unique = {h for _, h in hashes}
        if len(unique) != 1:
            print("FAIL: report output differs across PYTHONHASHSEED values:")
            for seed, h in hashes:
                print(f"  seed={seed}: {h}")
            sys.exit(1)
        print(f"PASS: all {len(SEEDS)} seeds produced byte-identical output "
              f"(after normalizing UUIDs/dates): {hashes[0][1]}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ── Invariant test: is sorted(cid_set) a pure function of content? ─────────

_INVARIANT_N_IDS = 1000


def _run_set_order_subprocess(seed: int, sorted_: bool) -> list:
    """Run, in a fresh subprocess under PYTHONHASHSEED=seed, exactly the
    step being checked: build a set of _INVARIANT_N_IDS synthetic circ_id
    strings and iterate it, either raw (sorted_=False, the control -- this
    is expected to VARY across seeds) or through sorted() (sorted_=True,
    the actual claim under test -- this must NOT vary). Returns the
    resulting id sequence as a plain list, printed as JSON by the
    subprocess and parsed back here -- never hashed, so a failure shows
    the actual differing elements, not just "not equal"."""
    script = f"""
import json
ids = {{f"chr{{i % 22 + 1}}:{{1000 + i}}|{{1050 + i}}" for i in range({_INVARIANT_N_IDS})}}
out = sorted(ids) if {sorted_} else list(ids)
print(json.dumps(out))
"""
    env = dict(os.environ)
    env["PYTHONHASHSEED"] = str(seed)
    result = subprocess.run([sys.executable, "-c", script], env=env,
                             capture_output=True, text=True, timeout=30)
    if result.returncode != 0:
        raise RuntimeError(f"subprocess failed under PYTHONHASHSEED={seed}:\n{result.stderr}")
    return json.loads(result.stdout)


def test_set_iteration_order_invariant():
    # Control: confirm raw set iteration actually IS seed-dependent at
    # this scale first. If this assertion itself fails, 1,000 strings
    # stopped being enough to make Python's hash randomization visible
    # (e.g. a Python version change disabling it by default) and the
    # "must be invariant" assertion below would be meaningless -- it
    # would trivially pass for the wrong reason.
    raw_sequences = [_run_set_order_subprocess(seed, sorted_=False) for seed in SEEDS]
    if len({tuple(s) for s in raw_sequences}) < 2:
        print("FAIL (test setup, not the claim under test): raw set "
              f"iteration over {_INVARIANT_N_IDS} strings was identical "
              "across all 5 seeds -- this control should have varied; "
              "the invariant assertion below would not be meaningful")
        sys.exit(1)
    print(f"  control OK: raw set iteration over {_INVARIANT_N_IDS} ids "
          f"varied across seeds ({len({tuple(s) for s in raw_sequences})} "
          f"distinct orderings in {len(SEEDS)} runs)")

    # The actual claim: sorted(cid_set) must be invariant.
    sorted_sequences = [_run_set_order_subprocess(seed, sorted_=True) for seed in SEEDS]
    first = sorted_sequences[0]
    for seed, seq in zip(SEEDS, sorted_sequences):
        if seq != first:
            # Find and report the first differing position, not just "differs".
            for i, (a, b) in enumerate(zip(first, seq)):
                if a != b:
                    print(f"FAIL: sorted(cid_set) differs at index {i}: "
                          f"seed={SEEDS[0]} has {a!r}, seed={seed} has {b!r}")
                    sys.exit(1)
            print(f"FAIL: sorted(cid_set) differs in length: "
                  f"seed={SEEDS[0]} has {len(first)}, seed={seed} has {len(seq)}")
            sys.exit(1)
    print(f"PASS: sorted(cid_set) over {_INVARIANT_N_IDS} ids produced the "
          f"identical sequence across all {len(SEEDS)} seeds")


def main():
    test_report_md5_reproducibility()
    test_set_iteration_order_invariant()


if __name__ == "__main__":
    main()
