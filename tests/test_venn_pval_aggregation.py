"""_is_better_pval() decides whether a newly-seen p-value should replace
the "best p-value across methods" tracked for one circRNA in
_venn_3_svg()'s circ_info aggregation (generate_report.py). NaN means
"this method didn't test/report this circRNA," not "its p-value equals
NaN" -- the bare `candidate < current` comparison this replaced gets that
wrong under IEEE 754 semantics: any comparison involving NaN is False, so
once the first method visited for a circ_id left a NaN there, no later
method's real value could ever win.

Named case: chr20:30904021|30904677 (KIF3B) in GSE113230 -- absent from
edgeR_ciriquant's de_results entirely, pval=NaN in deseq2's (visited
first, in de_files dict order), pval=0.0445603018127046 in limma's
(visited second). Before the fix, circ_info's stored pval for this
circRNA stayed NaN forever despite limma's real value existing; the Venn
detail list showed "NaN" for a circRNA that does have a real p-value.

Run with:  python tests/test_venn_pval_aggregation.py
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import generate_report as G  # noqa: E402

NAN = float("nan")
KIF3B_LIMMA_PVAL = 0.0445603018127046  # real value from GSE113230's de_results_limma.tsv


def test_nan_then_real_real_wins():
    # deseq2 visited first (NaN), limma visited second (real) -- the exact
    # KIF3B case.
    current = NAN
    assert G._is_better_pval(KIF3B_LIMMA_PVAL, current) is True, (
        "a real p-value must replace a NaN 'current' -- KIF3B's limma "
        "value should not be discarded because deseq2 reported NaN first")
    print("NaN-then-real (KIF3B case): real value wins -- OK")


def test_real_then_nan_real_stays():
    # Same circRNA, methods visited in the opposite order: real value seen
    # first, then a later method reports NaN. This direction already
    # worked before the fix (comparing a NaN candidate against anything
    # was already always False) -- locked in here so a future refactor
    # can't flip it.
    current = KIF3B_LIMMA_PVAL
    assert G._is_better_pval(NAN, current) is False, (
        "a NaN candidate must never overwrite an already-real current value")
    print("real-then-NaN: real value stays -- OK")


def test_all_nan_stays_nan():
    current = NAN
    result_is_better = G._is_better_pval(NAN, current)
    assert result_is_better is False, (
        "when every method reports NaN for a circRNA, the aggregate must "
        "remain NaN (there is no real value to prefer)")
    # Simulate the actual accumulation: current never changes because no
    # candidate was ever judged "better".
    assert current != current  # still NaN, i.e. math.isnan(current)
    assert math.isnan(current)
    print("all-NaN: result stays NaN -- OK")


def test_two_real_values_normal_comparison_unaffected():
    # Sanity check: the NaN-handling rewrite must not change ordinary
    # real-vs-real comparisons.
    assert G._is_better_pval(0.01, 0.05) is True
    assert G._is_better_pval(0.05, 0.01) is False
    assert G._is_better_pval(0.05, 0.05) is False  # not strictly better -> no replace
    print("real-vs-real comparisons unaffected -- OK")


def main():
    test_nan_then_real_real_wins()
    test_real_then_nan_real_stays()
    test_all_nan_stays_nan()
    test_two_real_values_normal_comparison_unaffected()
    print("_is_better_pval() regression tests passed")


if __name__ == "__main__":
    main()
