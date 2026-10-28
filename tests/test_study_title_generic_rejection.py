"""_reject_if_generic() must catch NCBI's own content-free placeholder
titles (e.g. "Homo sapiens Transcriptome or Gene expression", seen verbatim
from PRJNA553289's own BioProject record during this feature's development)
without being so broad that it also discards a real, short, title-shaped
string -- that second failure mode would be silent data loss: a genuine
title quietly replaced with "", indistinguishable from an actual NCBI
lookup failure. Both directions are tested here, on fixed strings, with no
network call.

Requires Flask + Flask-WTF importable (web_ui.py's own dependencies):
    conda run -n ciriquant python tests/test_study_title_generic_rejection.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

try:
    import web_ui as W  # noqa: E402
except ModuleNotFoundError as exc:
    print(f"SKIP (dependency not installed in this interpreter: {exc}); "
          f"rerun with: conda run -n ciriquant python {__file__}")
    sys.exit(0)


def test_known_generic_placeholder_is_rejected():
    # Verbatim string this feature's own development run got back from
    # db=bioproject's project_title for PRJNA553289 -- a real NCBI response,
    # not a synthetic example, which is exactly the kind of case the regex
    # must keep catching after any future edit.
    title, source = W._reject_if_generic(
        "Homo sapiens Transcriptome or Gene expression", "bioproject_desc")
    assert (title, source) == ("", "generic_rejected"), (
        f"expected the generic placeholder to be rejected, got {(title, source)!r}")
    print("known NCBI generic placeholder rejected: OK")


def test_real_short_title_is_not_rejected():
    # Verbatim value this feature's own /run_gse dry-run test fetched for
    # PRJNA808398 (db=bioproject's project_description) -- a real title,
    # just a short one, and the exact shape of string a too-broad regex
    # would start eating: "Capitalized-word lowercase-word ..." is the same
    # surface pattern as the placeholder above, minus the specific trailing
    # noun phrase the regex keys on.
    title, source = W._reject_if_generic(
        "RNA-Seq for Triple Negative Breast Cancer", "bioproject_desc")
    assert title == "RNA-Seq for Triple Negative Breast Cancer", (
        f"a real title was wrongly discarded as generic: got {(title, source)!r}")
    assert source == "bioproject_desc"
    print("real short title not mistaken for a generic placeholder: OK")


def main():
    test_known_generic_placeholder_is_rejected()
    test_real_short_title_is_not_rejected()
    print("study_title generic-rejection tests passed")


if __name__ == "__main__":
    main()
