"""Security regression tests for scripts/web_ui.py (OWASP ASVS L1 hardening).

Every test sends a *malicious* input and asserts it is rejected / neutralised. They run with
a Flask test client against a throwaway BASE_DIR: no server, no real data, no snakemake.
"""
import ast
import json
import re
from pathlib import Path

import pytest

import security as sec

ROOT = Path(__file__).resolve().parent.parent
EVIL = ["../../etc/passwd", "GSE1;rm -rf /", "--cores", "GSE1\x00", "A" * 500,
        "<script>alert(1)</script>", "GSE1 && id", "$(id)", "GSE12345/../../x"]


# ── pure validators ──────────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("ok", ["GSE113230", "srp156355", "PRJNA553289"])
def test_accession_whitelist_accepts(ok):
    assert sec.validate_accession(ok) == ok.upper()


@pytest.mark.parametrize("bad", EVIL + ["", "GSE", "GSEabc", "SRR1", "gse 1", "GSE1\nGSE2"])
def test_accession_whitelist_rejects(bad):
    with pytest.raises(sec.ValidationError):
        sec.validate_accession(bad)


@pytest.mark.parametrize("bad", EVIL + ["srr", "SRR", "SRR12a", "ERR123", "SRR1\nSRR2"])
def test_srr_rejects(bad):
    with pytest.raises(sec.ValidationError):
        sec.validate_srr(bad)


def test_project_id_cannot_look_like_option_or_path():
    for bad in ["-x", "--gse", "../x", "a/b", ".hidden", "1ABC", ""]:
        with pytest.raises(sec.ValidationError):
            sec.validate_project_id(bad)
    assert sec.validate_project_id("custom") == "CUSTOM"


def test_safe_join_blocks_escape(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    assert sec.safe_join(root, "a", "b.txt") == (root / "a" / "b.txt").resolve()
    for parts in [("..", "x"), ("a/../../x",), ("/etc/passwd",), ("x\x00",)]:
        with pytest.raises(sec.ValidationError):
            sec.safe_join(root, *parts)


def test_safe_join_blocks_symlink_escape(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    (tmp_path / "outside").mkdir()
    (root / "evil").symlink_to(tmp_path / "outside")
    with pytest.raises(sec.ValidationError):
        sec.safe_join(root, "evil", "file")


def test_safe_target_is_component_based_not_string_prefix(tmp_path):
    good = tmp_path / "home3" / "x"
    sibling = tmp_path / "home3" / "xevil"       # shares the string prefix "…/x"
    good.mkdir(parents=True)
    sibling.mkdir()
    (good / "a.fq.gz").write_text("")
    (sibling / "b.fq.gz").write_text("")
    assert sec.safe_target(str(good / "a.fq.gz"), [good])
    with pytest.raises(sec.ValidationError):      # old code used str.startswith and let this pass
        sec.safe_target(str(sibling / "b.fq.gz"), [good])


def test_safe_target_blocks_system_dirs_and_symlinks(tmp_path):
    root = tmp_path / "r"
    root.mkdir()
    link = root / "passwd"
    link.symlink_to("/etc/passwd")
    with pytest.raises(sec.ValidationError):
        sec.safe_target(str(link), [root, Path("/")])   # even if "/" were allowed, /etc is not
    with pytest.raises(sec.ValidationError):
        sec.safe_target("/etc/shadow", [Path("/")])


def test_job_id_strength_and_uniqueness():
    ids = {sec.new_job_id("GSE1") for _ in range(2000)}
    assert len(ids) == 2000
    suffix = next(iter(ids)).split("-", 1)[1]
    assert len(suffix) >= 12 and re.fullmatch(r"[A-Za-z0-9_-]+", suffix)
    sec.validate_job_id(next(iter(ids)))


def test_mask_link_hides_token():
    link = "https://h/auth/" + "A" * 20 + "WXYZ?lang=en"
    masked = sec.mask_link(link)
    assert "A" * 20 not in masked and masked.endswith("?lang=en")


def test_email_validation_blocks_header_injection():
    for bad in ["a@b.org\r\nBcc: x@y.org", "a@b.org,c@d.org", "a b@c.org", "x", ""]:
        with pytest.raises(sec.ValidationError):
            sec.validate_email(bad)
    assert sec.validate_email("Ab@Example.org") == "ab@example.org"


# ── HTTP-level: input validation ─────────────────────────────────────────────────────────────

def _queued(web):
    with web._qdb() as c:
        return c.execute("SELECT COUNT(*) FROM jobs").fetchone()[0]


@pytest.mark.parametrize("bad", EVIL)
def test_run_gse_rejects_bad_accession_without_echo(web, client, bad):
    r = client.post("/run_gse", data={"gse_id": bad, "cores": "4"})
    assert r.status_code == 400
    body = r.get_data(as_text=True)
    assert "alert(1)" not in body and bad.strip() not in body or bad.strip() == ""
    assert _queued(web) == 0


def test_run_gse_valid_still_works(web, client):
    r = client.post("/run_gse", data={"gse_id": "gse123456", "cores": "4"})
    assert r.status_code == 302 and r.headers["Location"].endswith("/queue")
    assert _queued(web) == 1
    with web._qdb() as c:
        row = c.execute("SELECT gse_id, command FROM jobs").fetchone()
    assert row["gse_id"] == "GSE123456"
    cmd = json.loads(row["command"])
    assert isinstance(cmd, list) and "GSE123456" in " ".join(cmd)


def test_run_gse_cores_garbage_does_not_500(web, client):
    r = client.post("/run_gse", data={"gse_id": "GSE1", "cores": "lots"})
    assert r.status_code == 302


@pytest.mark.parametrize("pid", ["../../tmp/pwn", "--cores", "a/b", "GSE1;id", "x" * 100])
def test_run_manual_rejects_bad_project_id(web, client, pid):
    r = client.post("/run_manual", data={"project_id": pid, "srr_id[]": ["SRR1"],
                                         "condition[]": ["tumor"]})
    assert r.status_code == 400 and _queued(web) == 0
    assert not (web.BASE_DIR.parent / "tmp").exists()
    assert not list(web.BASE_DIR.glob("metadata/*/library_info.csv"))


@pytest.mark.parametrize("srr", ["SRR1; id", "../SRR1", "-SRR1", "ERR1", "SRR1\n--x"])
def test_run_manual_rejects_bad_srr(web, client, srr):
    r = client.post("/run_manual", data={"project_id": "CUSTOM", "srr_id[]": [srr],
                                         "condition[]": ["tumor"]})
    assert r.status_code == 400 and _queued(web) == 0


def test_run_manual_rejects_bad_srr_in_csv(web, client):
    import io
    csv_bytes = b"srr_id,condition\nSRR1,tumor\n../../x,normal\n"
    r = client.post("/run_manual", data={"project_id": "CUSTOM",
                                         "csv_file": (io.BytesIO(csv_bytes), "a.csv")},
                    content_type="multipart/form-data")
    assert r.status_code == 400 and _queued(web) == 0


def test_run_manual_valid_still_works(web, client):
    r = client.post("/run_manual", data={"project_id": "custom", "srr_id[]": ["SRR100", "SRR101"],
                                         "condition[]": ["tumor", "normal"]})
    assert r.status_code == 302 and _queued(web) == 1
    assert (web.BASE_DIR / "metadata" / "CUSTOM" / "library_info.csv").exists()


def test_run_manual_bad_notify_email_header_injection(web, client):
    r = client.post("/run_manual", data={"project_id": "CUSTOM", "srr_id[]": ["SRR1"],
                                         "condition[]": ["tumor"],
                                         "notify_email": "a@b.org\r\nBcc: evil@x.org"})
    assert r.status_code == 400 and _queued(web) == 0


def _local(client, samples, **extra):
    data = {"project_id": "LOCAL", "samples_json": json.dumps(samples)}
    data.update(extra)
    return client.post("/run_local", data=data)


def test_run_local_blocks_name_traversal(web, client, tmp_path):
    src = web.TEST_RAW.parent / "r1.fastq.gz"
    src.write_text("")
    r = _local(client, [{"name": "../../escape", "r1": str(src), "r2": str(src),
                         "condition": "tumor"}])
    assert r.status_code == 400
    assert not list(tmp_path.rglob("escape*"))


def test_run_local_blocks_symlink_to_sensitive_and_sibling_prefix(web, client, tmp_path, monkeypatch):
    fq_root = tmp_path / "fq"
    sibling = tmp_path / "fqevil"                     # shares the *string* prefix of fq_root
    fq_root.mkdir()
    sibling.mkdir()
    (sibling / "s.fq.gz").write_text("")
    good = fq_root / "ok_1.fastq.gz"
    good.write_text("")
    monkeypatch.setenv("PIPELINE_FASTQ_ROOTS", str(fq_root))
    raw = web.TEST_RAW.parent.parent / "LOCAL" / "raw"      # run_local re-roots raw_dir under project id
    r = _local(client, [
        {"name": "etc", "r1": "/etc/passwd", "r2": "/etc/shadow", "condition": "tumor"},
        {"name": "sib", "r1": str(sibling / "s.fq.gz"), "r2": "", "condition": "tumor"},
        {"name": "good", "r1": str(good), "r2": str(good), "condition": "normal"},
    ])
    assert r.status_code in (302, 200)
    assert not (raw / "etc_1.fastq.gz").exists() and not (raw / "etc_2.fastq.gz").exists()
    assert not (raw / "sib_1.fastq.gz").exists()          # old startswith() check let this through
    assert (raw / "good_1.fastq.gz").is_symlink()          # legitimate use still works
    assert (raw / "good_1.fastq.gz").resolve() == good.resolve()


def test_run_local_rejects_non_list_json(web, client):
    r = client.post("/run_local", data={"project_id": "LOCAL", "samples_json": '{"a":1}'})
    assert r.status_code == 400
    r = client.post("/run_local", data={"project_id": "LOCAL", "samples_json": "{not json"})
    assert r.status_code == 400


def test_update_rejects_unknown_method_and_tools_and_survives_garbage(web, client):
    assert client.post("/update", data={"de_method": "x; rm -rf /"}).status_code == 400
    assert client.post("/update", data={"tools": ["ciriquant", "evil"]}).status_code == 400
    r = client.post("/update", data={"min_bsj": "NaN", "slop": "x", "fdr": "inf", "threads": "-5"})
    assert r.status_code == 302                      # old code: ValueError -> HTTP 500
    cfg = web.load_config()
    assert cfg["threads"] == 1 and cfg["de"]["fdr_cutoff"] == 0.05


def test_detect_labels_rejects_bad_accession(web, client):
    r = client.get("/api/detect_labels", query_string={"gse": "../../x"})
    assert r.status_code == 400


def test_scan_fastq_cannot_list_arbitrary_dirs(web, client):
    for path in ["/etc", "/", "/root", "../../etc", "/proc/self"]:
        r = client.get("/api/scan_fastq", query_string={"path": path})
        assert r.status_code == 403, path
        assert path not in r.get_data(as_text=True) or path == "/"


def test_scan_fastq_allows_configured_root(web, client):
    d = web.TEST_RAW.parent / "fq"
    d.mkdir()
    (d / "s_1.fastq.gz").write_text("")
    (d / "s_2.fastq.gz").write_text("")
    r = client.get("/api/scan_fastq", query_string={"path": str(d)})
    assert r.status_code == 200 and r.get_json()["count"] == 1


# ── reflected XSS / error messages ───────────────────────────────────────────────────────────

XSS_JOB_IDS = ['";alert(1);"', '";alert(1)//', "</script><script>alert(1)</script>", "x'-alert(1)-'", "%3Cimg%20src=x%3E"]


@pytest.mark.parametrize("jid", XSS_JOB_IDS)
def test_status_page_does_not_reflect_job_id(client, jid):
    r = client.get("/status/" + jid)
    body = r.get_data(as_text=True)
    assert "alert(1)" not in body and "<img src=x" not in body
    assert "const JOB_ID = null;" in body or r.status_code == 404


def test_status_page_known_job_id_is_json_encoded(web, client):
    web.queue_add("GSE1-abcDEF_123", "GSE1", ["x"], str(web.BASE_DIR / "logs" / "a.log"), 4, "")
    r = client.get("/status/GSE1-abcDEF_123")
    assert 'const JOB_ID = "GSE1-abcDEF_123";' in r.get_data(as_text=True)


def test_qc_pending_page_does_not_leak_server_path(web, client):
    web.save_registry({"GSE1-abc": {"gse_id": "GSE1", "log": "logs/a.log"}})
    r = client.get("/qc/GSE1-abc")
    assert r.status_code == 202
    assert str(web.BASE_DIR) not in r.get_data(as_text=True) and "/qc/" not in r.get_data(as_text=True).split("<h2>")[0]


def test_json_for_script_neutralises_script_breakout(web):
    out = web._json_for_script({"gene": "</script><script>alert(1)</script>"})
    assert "</script>" not in out and "<" not in out
    assert json.loads(out)["gene"] == "</script><script>alert(1)</script>"


# ── login: Host-header poisoning, allow-list, enumeration, token hygiene ─────────────────────

@pytest.fixture()
def mail(web, monkeypatch):
    sent = []
    monkeypatch.setattr(web, "send_magic_link", lambda to, link, lang="zh": sent.append((to, link)) or True)
    return sent


def test_magic_link_ignores_spoofed_host(web, client, mail):
    r = client.post("/login", data={"email": "victim@example.org"},
                    headers={"Host": "evil.example.net", "X-Forwarded-Host": "evil.example.net"})
    assert r.status_code == 200 and len(mail) == 1
    assert mail[0][1].startswith("https://circdex.example.org/auth/")
    assert "evil" not in mail[0][1]


def test_forwarded_host_not_trusted_even_behind_proxy(web):
    # ProxyFix (if enabled) is configured with x_proto/x_for only, never x_host.
    assert "x_host" not in open(ROOT / "scripts" / "web_ui.py").read().split("ProxyFix(")[1].split(")")[0]


def test_login_rejects_malformed_email_without_sending(web, client, mail):
    for bad in ["a@b.org\r\nBcc: x@y.org", "no-at-sign", "a@b", "<x>@y.org"]:
        r = client.post("/login", data={"email": bad})
        assert r.status_code == 200
    assert mail == []


def test_allowlist_blocks_unlisted_but_looks_identical(web, client, mail, monkeypatch):
    monkeypatch.setattr(web, "ALLOWED_EMAILS", {"ok@example.org"})
    r1 = client.post("/login", data={"email": "stranger@example.org"})
    r2 = client.post("/login", data={"email": "ok@example.org"})
    assert [m[0] for m in mail] == ["ok@example.org"]          # stranger got nothing
    strip = lambda t, e: t.replace(e, "")
    assert strip(r1.get_data(as_text=True), "stranger@example.org") == \
        strip(r2.get_data(as_text=True), "ok@example.org")     # no account enumeration


def test_magic_link_token_is_strong_and_not_logged(web, client, mail, capsys):
    client.post("/login", data={"email": "a@example.org"})
    token = mail[0][1].split("/auth/")[1].split("?")[0]
    assert len(token) >= 43                                      # token_urlsafe(32)
    out = capsys.readouterr()
    assert token not in out.out and token not in out.err


def test_auth_open_redirect_blocked(web, client, mail):
    client.post("/login", data={"email": "a@example.org"})
    token = mail[0][1].split("/auth/")[1]
    r = client.get("/auth/" + token, query_string={"next": "https://evil.example.net/x"})
    assert r.status_code == 302 and "evil" not in r.headers["Location"]


# ── headers, CSRF, authn ─────────────────────────────────────────────────────────────────────

def test_security_headers_present(client):
    r = client.get("/", base_url="https://circdex.example.org")
    h = r.headers
    assert h["X-Content-Type-Options"] == "nosniff"
    assert h["X-Frame-Options"] == "DENY"
    assert h["Referrer-Policy"] == "no-referrer"
    assert "max-age=" in h["Strict-Transport-Security"]
    csp = h["Content-Security-Policy"]
    assert "frame-ancestors 'none'" in csp and "object-src 'none'" in csp


def test_generated_report_is_sandboxed(web, client):
    cfgdir = Path(web.load_config()["results_dir"])
    cfgdir.mkdir(parents=True)
    (cfgdir / "report.html").write_text("<html><script>1</script></html>")
    web.save_registry({"GSE1-abc": {"gse_id": "GSE1", "log": "logs/a.log"}})
    r = client.get("/report/GSE1-abc")
    assert r.status_code == 200
    assert r.headers["Content-Security-Policy"].startswith("sandbox allow-scripts")


def test_csrf_enforced_on_state_changing_posts(web):
    web.app.config["WTF_CSRF_ENABLED"] = True
    c = web.app.test_client()
    with c.session_transaction() as s:
        s["email"] = "t@example.org"
        s["session_expires"] = "2999-01-01 00:00:00"
    for url, data in [("/run_gse", {"gse_id": "GSE1"}), ("/update", {}),
                      ("/run_manual", {"project_id": "X"}), ("/run_local", {})]:
        assert c.post(url, data=data).status_code == 400, url
    assert _queued(web) == 0


@pytest.mark.parametrize("url", ["/", "/queue", "/status/x", "/report/x", "/download/x",
                                 "/api/scan_fastq?path=/", "/api/progress", "/cross_dataset"])
def test_everything_requires_login(web, url):
    c = web.app.test_client()
    r = c.get(url)
    assert r.status_code == 302 and "/login" in r.headers["Location"]


def test_session_cookie_flags(web):
    assert web.app.config["SESSION_COOKIE_HTTPONLY"] is True
    assert web.app.config["SESSION_COOKIE_SAMESITE"] == "Lax"
    assert web.app.config["MAX_CONTENT_LENGTH"] <= 10 * 1024 * 1024


# ── static guard: no shell-string subprocess calls ───────────────────────────────────────────

def test_no_shell_string_subprocess_calls():
    """Every subprocess.* call must pass a list/variable argv and never shell=True."""
    problems = []
    for py in [ROOT / "scripts" / "web_ui.py", ROOT / "scripts" / "security.py"]:
        tree = ast.parse(py.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            f = node.func
            name = None
            if isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name):
                name = f.value.id + "." + f.attr
            if name in ("os.system", "os.popen"):
                problems.append("{}:{} {}".format(py.name, node.lineno, name))
            if name and name.startswith("subprocess."):
                for kw in node.keywords:
                    if kw.arg == "shell" and not (isinstance(kw.value, ast.Constant) and kw.value.value is False):
                        problems.append("{}:{} shell=...".format(py.name, node.lineno))
                if name != "subprocess.PIPE" and node.args:
                    a0 = node.args[0]
                    if isinstance(a0, (ast.Constant, ast.JoinedStr, ast.BinOp)):
                        problems.append("{}:{} string argv".format(py.name, node.lineno))
    assert not problems, problems


def test_secure_cookie_default_follows_public_url(monkeypatch):
    """Secure cookie is on for an https public origin, but not forced for plain-HTTP intranet use."""
    import importlib
    import web_ui
    for env, public, expect in [(None, "https://x.org", True), (None, "", False),
                                ("0", "https://x.org", False), ("1", "", True)]:
        if env is None:
            monkeypatch.delenv("PIPELINE_COOKIE_SECURE", raising=False)
        else:
            monkeypatch.setenv("PIPELINE_COOKIE_SECURE", env)
        monkeypatch.setenv("PIPELINE_PUBLIC_URL", public)
        mod = importlib.reload(web_ui)
        assert mod.app.config["SESSION_COOKIE_SECURE"] is expect, (env, public)
    monkeypatch.setenv("PIPELINE_PUBLIC_URL", "https://circdex.example.org")
    monkeypatch.setenv("PIPELINE_COOKIE_SECURE", "0")
    importlib.reload(web_ui)


def test_dev_print_link_is_opt_in(web, client, mail, capsys, monkeypatch):
    client.post("/login", data={"email": "a@example.org"})
    assert "[dev] magic link" not in capsys.readouterr().out          # default: token never printed
    monkeypatch.setenv("PIPELINE_DEV_PRINT_LINK", "1")
    client.post("/login", data={"email": "a@example.org"})
    assert "[dev] magic link for a@example.org: https://circdex.example.org/auth/" in capsys.readouterr().out
