"""The README's minimal example must keep producing exactly the documented output."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EX = ROOT / "examples" / "minimal"


def test_minimal_example_matches_expected(tmp_path):
    bed, summ = tmp_path / "o.bed", tmp_path / "s.tsv"
    r = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "consensus_filter.py"),
         "--cirique", str(EX / "input" / "SAMPLE.gtf"),
         "--dcc", str(EX / "input" / "DCC" / "CircCoordinates"),
         "--output", str(bed), "--summary", str(summ),
         "--min-tools", "2", "--slop", "10", "--min-bsj", "2",
         "--max-junction-ratio", "1.0", "--qc-bsj-threshold", "5",
         "--adaptive", "--adaptive-ratio", "0.1"],
        capture_output=True, text=True, cwd=str(ROOT), timeout=60)
    assert r.returncode == 0, r.stderr
    assert bed.read_text() == (EX / "expected" / "high_confidence.bed").read_text()
    assert summ.read_text() == (EX / "expected" / "consensus_summary.tsv").read_text()


def test_documented_scenarios_hold(tmp_path):
    """Spot-check the semantics promised in examples/minimal/README.md."""
    lines = (EX / "expected" / "consensus_summary.tsv").read_text().splitlines()[1:]
    ids = {l.split("\t")[0] for l in lines}
    assert ids == {"chr1:1000|2000", "chr1:5000|6000", "chr6:7000|8000", "chr7:4000|4800"}
    assert not ({"chr2:3000|4000", "chr3:800|1500", "chr4:100|900", "chr5:2000|2600", "chr9:100|500"} & ids)
