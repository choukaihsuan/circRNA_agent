"""Every figure under audit/figs/ must have a committed script that produces it.

Background: a figure that went into the manuscript once had no discoverable source and it took a
full round to reconstruct it. Rule: no orphan images.
"""
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT_SUFFIXES = {".py", ".R", ".r", ".sh", ".ipynb"}


def find_orphans(audit_dir: Path):
    """PNGs under ``audit_dir/figs`` that no script under ``audit_dir`` (outside figs/) mentions."""
    figs = audit_dir / "figs"
    if not figs.is_dir():
        return []
    scripts = [p for p in audit_dir.rglob("*")
               if p.is_file() and p.suffix in SCRIPT_SUFFIXES and figs not in p.parents]
    corpus = {p: p.read_text(errors="replace") for p in scripts}
    orphans = []
    for png in sorted(figs.rglob("*.png")):
        if not any(png.name in t or png.stem in t for t in corpus.values()):
            orphans.append(png.relative_to(audit_dir.parent).as_posix())
    return orphans


def test_no_orphan_figures_in_repo():
    figs = ROOT / "audit" / "figs"
    if not figs.is_dir() or not list(figs.rglob("*.png")):
        pytest.skip("no audit/figs/*.png in this checkout")
    assert find_orphans(ROOT / "audit") == [], \
        "figures without a producing script under audit/ (commit the script)"


# The checker itself, tested against a synthetic tree (the real tree has no figures yet).
def test_checker_flags_orphan_and_accepts_matched(tmp_path):
    audit = tmp_path / "audit"
    (audit / "figs").mkdir(parents=True)
    (audit / "figs" / "fig_ok.png").write_bytes(b"\x89PNG")
    (audit / "figs" / "fig_orphan.png").write_bytes(b"\x89PNG")
    (audit / "make_ok.py").write_text("fig.savefig('audit/figs/fig_ok.png')\n")
    assert find_orphans(audit) == ["audit/figs/fig_orphan.png"]


def test_checker_ignores_scripts_inside_figs_dir(tmp_path):
    audit = tmp_path / "audit"
    (audit / "figs").mkdir(parents=True)
    (audit / "figs" / "a.png").write_bytes(b"\x89PNG")
    (audit / "figs" / "a.py").write_text("# mentions a.png but lives in figs/ - not a real producer\n")
    assert find_orphans(audit) == ["audit/figs/a.png"]
