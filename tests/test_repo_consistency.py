"""Cheap guards against documentation / config drift that has bitten this project."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _tags(text):
    return set(re.findall(r"circrna-pipeline[:_]([0-9]+\.[0-9]+\.[0-9]+)", text))


def test_container_image_tag_is_consistent():
    """Snakefile `singularity:` pinned :1.0.0 while README/Dockerfile said 1.0.1 (found 2026-10-09)."""
    snake = _tags((ROOT / "workflow" / "Snakefile").read_text())
    readme = _tags((ROOT / "README.md").read_text())
    build = set(re.findall(r'IMAGE_TAG="([0-9.]+)"', (ROOT / "containers" / "build_and_deploy.sh").read_text()))
    label = set(re.findall(r'version="([0-9.]+)"', (ROOT / "Dockerfile").read_text()))
    assert len(snake) == 1 and snake == readme == build == label, (snake, readme, build, label)


def test_claude_md_never_tells_anyone_to_use_bare_kill():
    """Project rule: always kill -9 (SIGTERM makes Snakemake delete outputs)."""
    text = (ROOT / "CLAUDE.md").read_text()
    bad = [l for l in text.splitlines() if re.search(r"\bpkill\s+(?!-9)\S", l)]
    assert not bad, bad


def test_claude_md_de_sig_by_default_matches_config():
    import yaml
    default = yaml.safe_load((ROOT / "config.yaml").read_text())["de"]["de_sig_by"]
    text = (ROOT / "CLAUDE.md").read_text()
    assert "  de_sig_by:           %s" % default in text, "CLAUDE.md sample config out of date"
    for rule in ("de.smk",):
        assert 'get("de_sig_by", "%s")' % default in (ROOT / "workflow" / "rules" / rule).read_text()


def test_no_paper_content_in_repo():
    """The manuscript lives elsewhere: no Word/TeX/BibTeX files and no discussion-draft section."""
    tracked = [p for p in ROOT.rglob("*") if p.is_file() and ".git" not in p.parts
               and "node_modules" not in p.parts and p.suffix.lower() in {".docx", ".doc", ".tex", ".bib"}]
    assert not tracked, tracked
    assert "## 論文 Discussion 素材" not in (ROOT / "CLAUDE.md").read_text()
