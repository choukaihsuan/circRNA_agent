"""Shared fixtures. Tests run without real data: everything lives in tmp dirs."""
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

# Must be set before web_ui is imported (module-level config reads them).
os.environ.setdefault("PIPELINE_SECRET_KEY", "test-secret-key-not-for-production")
os.environ.setdefault("PIPELINE_COOKIE_SECURE", "0")
os.environ["PIPELINE_PUBLIC_URL"] = "https://circdex.example.org"
os.environ.pop("PIPELINE_ALLOWED_EMAILS", None)


@pytest.fixture()
def web(tmp_path, monkeypatch):
    """web_ui module re-pointed at a throwaway BASE_DIR."""
    import web_ui

    base = tmp_path / "repo"
    (base / "config" / "projects").mkdir(parents=True)
    (base / "metadata").mkdir()
    (base / "logs").mkdir()
    (base / "jobs").mkdir()
    raw = tmp_path / "data" / "GSE1" / "raw"
    raw.mkdir(parents=True)
    cfg = {
        "project_id": "GSE1",
        "raw_dir": str(raw),
        "trimmed_dir": str(tmp_path / "data" / "GSE1" / "trimmed"),
        "results_dir": str(tmp_path / "data" / "GSE1_results"),
        "threads": 8,
        "consensus": {"tools": ["ciriquant", "dcc"], "min_tools": 2},
        "de": {"method": "edgeR_ciriquant", "tumor_label": "tumor", "normal_label": "normal"},
        "download": {},
    }
    import yaml
    (base / "config.yaml").write_text(yaml.safe_dump(cfg))

    monkeypatch.setattr(web_ui, "BASE_DIR", base)
    monkeypatch.setattr(web_ui, "CONFIG_PATH", base / "config.yaml")
    monkeypatch.setattr(web_ui, "REGISTRY_PATH", base / "jobs" / "registry.json")
    monkeypatch.setattr(web_ui, "QUEUE_DB", base / "jobs" / "queue.db")
    monkeypatch.setattr(web_ui, "LOG_PATH", base / "logs" / "pipeline_run.log")
    monkeypatch.setattr(web_ui, "_snake_bin", lambda: "snakemake")
    monkeypatch.setattr(web_ui, "_fetch_geo_title", lambda gse: "")
    web_ui.init_queue_db()
    web_ui.app.config.update(TESTING=True, WTF_CSRF_ENABLED=False)
    web_ui.TEST_RAW = raw
    return web_ui


@pytest.fixture()
def client(web):
    c = web.app.test_client()
    with c.session_transaction() as s:
        s["email"] = "tester@example.org"
        s["session_expires"] = "2999-01-01 00:00:00"
    return c
