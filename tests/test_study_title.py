"""study_title data flow: GEO response -> config snapshot -> report HTML."""
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import geo_title  # noqa: E402


@pytest.mark.parametrize("acc,uid", [
    ("GSE113230", "200113230"),     # 6 digits: old formula happened to be right
    ("GSE58135",  "200058135"),     # 5 digits: old formula gave 20058135 (another record)
    ("GSE77509",  "200077509"),
    ("gse97239",  "200097239"),
    ("GSE1234567", "201234567"),    # 7 digits: old formula gave 2001234567
])
def test_geo_uid(acc, uid):
    assert geo_uid_ok(acc) == uid


def geo_uid_ok(acc):
    return geo_title.geo_uid(acc)


@pytest.mark.parametrize("bad", ["SRP156355", "PRJNA553289", "GSE", "GSE1;x", "", "CUSTOM"])
def test_geo_uid_rejects_non_gse(bad):
    assert geo_title.geo_uid(bad) is None


class _Resp:
    def __init__(self, payload):
        self._p = payload

    def read(self):
        return self._p

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def test_fetch_refuses_title_of_a_different_accession(monkeypatch):
    import json
    payload = json.dumps({"result": {"200058135": {"accession": "GSE99999", "title": "SOMEONE ELSE"}}}).encode()
    monkeypatch.setattr(geo_title.urllib.request, "urlopen", lambda *a, **k: _Resp(payload))
    assert geo_title.fetch_geo_title("GSE58135") == ""          # never attach another series' title


def test_fetch_returns_matching_title(monkeypatch):
    import json
    payload = json.dumps({"result": {"200058135": {"accession": "GSE58135", "title": "Right one"}}}).encode()
    seen = {}
    monkeypatch.setattr(geo_title.urllib.request, "urlopen",
                        lambda url, **k: seen.setdefault("url", url) and _Resp(payload))
    assert geo_title.fetch_geo_title("GSE58135") == "Right one"
    assert "id=200058135" in seen["url"]


def test_new_project_does_not_inherit_previous_projects_title(web):
    cfg = web.load_config()
    cfg.update(project_id="GSE111", study_title="Title of the PREVIOUS project", notify={"email_to": "a@b.org"})
    web.save_config(cfg)
    fresh = web.load_project_config("GSE222")                    # no snapshot -> falls back to global
    assert "study_title" not in fresh and "notify" not in fresh
    same = web.load_project_config("GSE111")                     # same project: global config is its own
    assert same["study_title"] == "Title of the PREVIOUS project"


def test_run_gse_fetches_title_for_the_new_dataset(web, client, monkeypatch):
    cfg = web.load_config()
    cfg.update(project_id="GSE111", study_title="Title of the PREVIOUS project")
    web.save_config(cfg)
    monkeypatch.setattr(web, "_fetch_geo_title", lambda g: "Title of " + g)
    r = client.post("/run_gse", data={"gse_id": "GSE222", "cores": "4"})
    assert r.status_code == 302
    snap = yaml.safe_load((web.BASE_DIR / "config" / "projects" / "GSE222.yaml").read_text())
    assert snap["study_title"] == "Title of GSE222"


def test_run_gse_existing_snapshot_keeps_its_title(web, client, monkeypatch):
    snap = web.load_config()
    snap.update(project_id="GSE222", study_title="Curated title")
    (web.BASE_DIR / "config" / "projects" / "GSE222.yaml").write_text(yaml.safe_dump(snap))
    monkeypatch.setattr(web, "_fetch_geo_title", lambda g: "NETWORK TITLE")
    client.post("/run_gse", data={"gse_id": "GSE222", "cores": "4"})
    out = yaml.safe_load((web.BASE_DIR / "config" / "projects" / "GSE222.yaml").read_text())
    assert out["study_title"] == "Curated title"


def test_update_endpoint_stores_user_title_capped(web, client):
    client.post("/update", data={"study_title": "  My title  " + "x" * 900})
    t = web.load_config()["study_title"]
    assert t.startswith("My title") and len(t) <= 500


def test_report_escapes_study_title_and_labels(tmp_path):
    pytest.importorskip("pandas")
    import generate_report as gr
    groups = tmp_path / "g.csv"
    groups.write_text('srr_id,condition,patient_id\nSRR1,tumor,"<img src=x onerror=alert(2)>"\nSRR2,normal,P2\n')
    html = gr._sample_overview_section(str(groups), "tumor", "normal", None,
                                       study_title='<script>alert(1)</script> & "quoted"')
    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;alert(1)&lt;/script&gt; &amp; &quot;quoted&quot;" in html
    assert "<img src=x" not in html and "&lt;img src=x" in html


def test_audit_script_flags_mismatch_and_not_applicable(tmp_path):
    import audit_study_titles as a
    (tmp_path / "GSE1.yaml").write_text(yaml.safe_dump({"study_title": "Right"}))
    (tmp_path / "GSE2.yaml").write_text(yaml.safe_dump({"study_title": "Wrong"}))
    (tmp_path / "GSE3.yaml").write_text(yaml.safe_dump({}))
    (tmp_path / "SRP9.yaml").write_text(yaml.safe_dump({"study_title": "x"}))
    geo = {"GSE1": "Right", "GSE2": "Real", "GSE3": "Real3"}
    rows = {r["project"]: r["status"] for r in a.audit(tmp_path, fetch=lambda p: geo.get(p, ""))}
    assert rows == {"GSE1": "OK", "GSE2": "MISMATCH", "GSE3": "EMPTY", "SRP9": "N-A"}
