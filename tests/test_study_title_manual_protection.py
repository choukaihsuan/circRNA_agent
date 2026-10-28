"""_maybe_autofetch_study_title() must never overwrite a human-typed
("manual") study_title with an automatic fetch. This is the highest-stakes
property in the whole study_title feature: the seven manually-written
titles in config/projects/*.yaml (GSE108735, GSE143797, GSE192410,
GSE192849, GSE97239, PRJNA553289, SRP156355) have no backup source -- if a
future refactor breaks this guard and a GEO/SRA auto-fetch silently
overwrites one of them, the original text is gone for good. A regression
here is irreversible in a way the other two study_title bugs (missing
/run_gse call, GEO UID padding) are not: both of those just produce an
empty or stale title, which a re-fetch can always repair after the fact.

This test never touches the network -- it monkeypatches urllib.request's
urlopen to raise on any call, so if the manual-protection guard is ever
removed or weakened, this test fails on "urlopen was called" rather than on
a wrong answer that happened to still look plausible.

Requires Flask + Flask-WTF importable (web_ui.py's own dependencies):
    conda run -n ciriquant python tests/test_study_title_manual_protection.py
"""
import sys
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

try:
    import web_ui as W  # noqa: E402
except ModuleNotFoundError as exc:
    print(f"SKIP (dependency not installed in this interpreter: {exc}); "
          f"rerun with: conda run -n ciriquant python {__file__}")
    sys.exit(0)


def _urlopen_must_not_be_called(*a, **k):
    raise AssertionError(
        "urlopen() was called -- the manual-protection guard in "
        "_maybe_autofetch_study_title() failed to short-circuit before any "
        "network fetch, which means it would have overwritten a "
        "human-written title")


def test_manual_title_is_never_overwritten():
    cfg = {
        "study_title": "RNA-seq profiling of small cell lung cancer tumor and adjacent normal lung tissues (Nanjing Chest Hospital)",
        "study_title_source": "manual",
    }
    original = dict(cfg)
    with mock.patch.object(W.urllib.request, "urlopen", side_effect=_urlopen_must_not_be_called):
        W._maybe_autofetch_study_title(cfg, "PRJNA553289")
    assert cfg == original, (
        f"manual study_title was modified: before={original!r}, after={cfg!r}")
    print("manual title protected from overwrite: OK")


def test_geo_sourced_title_is_still_refreshable():
    """The guard must only block source=='manual'. A previously auto-fetched
    title (source in geo/sra_study/bioproject_desc/generic_rejected/
    fetch_failed) is not a deliberate human choice and must remain
    refreshable on a later call -- otherwise a dataset whose first fetch
    failed, or whose GEO record later gets corrected, could never recover."""
    cfg = {"study_title": "stale cached title", "study_title_source": "geo"}

    def fake_urlopen(url, timeout=8):
        class _Resp:
            def __enter__(self):
                return self
            def __exit__(self, *exc):
                return False
            def read(self):
                import json
                uid = url.split("id=")[1].split("&")[0]
                return json.dumps({"result": {uid: {"title": "fresh refetched title"}}}).encode()
        return _Resp()

    with mock.patch.object(W.urllib.request, "urlopen", side_effect=fake_urlopen):
        W._maybe_autofetch_study_title(cfg, "GSE113230")

    assert cfg["study_title"] == "fresh refetched title", (
        f"a source=='geo' title should have been refreshed, got {cfg!r}")
    assert cfg["study_title_source"] == "geo"
    print("geo-sourced title remains refreshable (not mistaken for manual): OK")


def main():
    test_manual_title_is_never_overwritten()
    test_geo_sourced_title_is_still_refreshable()
    print("study_title manual-protection tests passed")


if __name__ == "__main__":
    main()
