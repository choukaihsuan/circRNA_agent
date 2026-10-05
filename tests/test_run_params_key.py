"""canonical_key() must follow the program files a step runs, not the repository version.

Run with:  python tests/test_run_params_key.py
"""
import os
import sys
import tempfile
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import run_params as rp  # noqa: E402


def record_for(res_dir, step_file, started="2026-10-05T10:00:00+08:00"):
    path = rp.write_record(
        res_dir, "demo_step", wildcard="w1", started_at=started,
        effective_params={"alpha": 0.05, "beta": 2},
        inputs=[],
        script=step_file, program_files=[step_file],
    )
    return yaml.safe_load(open(path, encoding="utf-8"))


def main():
    with tempfile.TemporaryDirectory() as tmp:
        step = os.path.join(tmp, "step.py")
        other = os.path.join(tmp, "unrelated.py")
        Path(step).write_text("print('step v1')\n")
        Path(other).write_text("print('other v1')\n")
        res = os.path.join(tmp, "results")

        base = record_for(res, step)
        key0 = rp.canonical_key(base)

        # 1. an unrelated file changes: the key must not change
        Path(other).write_text("print('other v2 — unrelated edit')\n")
        after_unrelated = record_for(res, step, started="2026-10-06T10:00:00+08:00")
        assert rp.canonical_key(after_unrelated) == key0, "unrelated edit changed the key"

        # 2. the repository HEAD (information only) changes: the key must not change
        moved = dict(after_unrelated)
        moved["code"] = dict(after_unrelated["code"], repo_head="0" * 40)
        assert rp.canonical_key(moved) == key0, "repo_head changed the key"

        # 3. a program file the step runs changes: the key must change
        Path(step).write_text("print('step v2 — real edit')\n")
        after_step = record_for(res, step, started="2026-10-07T10:00:00+08:00")
        assert rp.canonical_key(after_step) != key0, "edit to the step's own program did not change the key"

        # 4. the key ignores timestamps and environment
        ts = dict(after_step)
        ts["started_at"], ts["finished_at"] = "2030-01-01T00:00:00+00:00", "2030-01-01T00:00:01+00:00"
        ts["environment"] = {"hostname": "elsewhere"}
        assert rp.canonical_key(ts) == rp.canonical_key(after_step)

        # 5. a changed parameter changes the key
        p = dict(after_step)
        p["effective_params"] = dict(after_step["effective_params"], beta=3)
        assert rp.canonical_key(p) != rp.canonical_key(after_step)

        # 6. the program-file list is recorded with its hash
        assert after_step["program_files"][0]["sha256"] and after_step["program_files"][0]["path"].endswith("step.py")
    print("run_params canonical_key tests passed")


if __name__ == "__main__":
    main()
