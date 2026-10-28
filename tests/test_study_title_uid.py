"""fetch_study_title()'s GEO UID formula must zero-pad the GSE number to 6
digits before prepending "200" -- a bare "200" + digits concatenation
produces a UID that matches no GEO record for any accession under 6 digits
(GSE1 .. GSE999999), and esummary.fcgi then returns an empty result, which
this function maps to ("", "fetch_failed") with no indication the UID
itself -- not the network, not the accession -- was wrong. That is exactly
how this bug hid for GSE58135/GSE77509/GSE97239: it failed silently, with a
return value indistinguishable from a genuine network/lookup failure.

This test imports the real fetch_study_title() from scripts/web_ui.py (not
a reimplementation of the formula) and intercepts urllib.request.urlopen to
read the UID it actually put in the request URL, without making a real
request -- so a future refactor that drops the zfill is caught here even
though no network call happens.

Requires Flask + Flask-WTF importable (web_ui.py's own dependencies) -- run
with the project's conda env, same as web_ui.py itself:
    conda run -n ciriquant python tests/test_study_title_uid.py
A bare `python3 tests/test_study_title_uid.py` outside that env will fail
on the `import web_ui` with ModuleNotFoundError: flask, which is a missing
dependency, not a test failure -- this script exits with a distinct message
in that case rather than reporting a pass or a silent skip.
"""
import io
import json
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


class _FakeResponse:
    """Minimal stand-in for the object urllib.request.urlopen() returns,
    used as a context manager the same way fetch_study_title() uses it."""
    def __init__(self, payload: bytes):
        self._payload = payload

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read(self):
        return self._payload


def _geo_uid_used_in_request(accession: str) -> str:
    """Calls the real fetch_study_title(), capturing the UID it embeds in
    the esummary.fcgi URL, by returning a canned (empty) GEO response so the
    function completes without ever touching the network."""
    captured = {}

    def fake_urlopen(url, timeout=8):
        captured["url"] = url
        # An esummary response with no matching uid -- shaped like what
        # real NCBI servers return for a UID that doesn't exist, so the
        # function's own "if not title: return fetch_failed" path is what
        # ends this call, not an exception from a malformed fake.
        return _FakeResponse(json.dumps({"result": {}}).encode())

    with mock.patch.object(W.urllib.request, "urlopen", side_effect=fake_urlopen):
        W.fetch_study_title(accession)

    assert "url" in captured, (
        f"fetch_study_title({accession!r}) never called urlopen -- "
        f"the accession likely failed the GSE regex before reaching the "
        f"UID computation this test exists to check")
    # id={uid}& -- uid is everything between "id=" and the next "&"
    url = captured["url"]
    marker = "id="
    start = url.index(marker) + len(marker)
    end = url.index("&", start)
    return url[start:end]


def main():
    cases = {
        "GSE58135":  "200058135",   # 5-digit GSE number -- the bug's signature case
        "GSE97239":  "200097239",   # 5-digit
        "GSE121842": "200121842",   # 6-digit, zfill is a no-op here
        "GSE323364": "200323364",   # 6-digit, zfill is a no-op here
    }
    for accession, expected_uid in cases.items():
        got = _geo_uid_used_in_request(accession)
        assert got == expected_uid, (
            f"{accession}: fetch_study_title() put UID {got!r} in its "
            f"request, expected {expected_uid!r} -- a bare \"200\"+digits "
            f"concatenation (no zfill) would silently compute the wrong "
            f"UID for any GSE number under 6 digits")
    print("fetch_study_title() GEO UID zero-padding test passed")


if __name__ == "__main__":
    main()
