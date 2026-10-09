"""Config -> command-line wiring tests.

Guards against *silent wiring gaps*: a config option exists (and the UI lets you change it) but
never reaches the script that should use it. This has happened three times in this project
(``components_arg``, ``sig_cap``/``fc_cap``, ``tools.ciriquant.anchor``); each time it was only
found by comparing code line by line.

What is checked
  1. Dynamic: every scalar under ``consensus``, ``de`` and ``tools`` in ``config.yaml`` (plus every
     key the rules read with a ``.get(key, default)`` but that config.yaml does not declare) is set
     to a clearly deviating value; the real shell commands are re-rendered with a Snakemake
     dry-run; the new value must appear literally in a command (or change the DAG for the
     keys that select samples/tools).
  2. Static, for ``script:`` rules (analysis.R, isoform_switching.R, generate_report.py), whose
     params are not visible in a dry-run: every ``params`` entry must be consumed by the script and
     every ``snakemake@params[[...]]`` / ``snakemake.params.x`` the script reads must be defined.
  3. Config keys that reach neither (and are not listed in KNOWN_GAPS) fail.

No snakemake execution: ``-n`` (dry-run) only, inside a throw-away fixture, no real data.
"""
from __future__ import annotations

import ast
import copy
import os
import re
from pathlib import Path

import pytest

import wiring as w

ROOT = w.ROOT
RULES = sorted((ROOT / "workflow" / "rules").glob("*.smk"))
SECTIONS = ("consensus", "de", "tools")

pytestmark = pytest.mark.skipif(
    w.snakemake_bin() is None and not os.environ.get("REQUIRE_SNAKEMAKE"),
    reason="snakemake not installed (set REQUIRE_SNAKEMAKE=1 to make this a failure)")


# ── what we know about each key ──────────────────────────────────────────────────────────────

# Keys whose effect is *which jobs exist* rather than a literal flag value.
# value to set -> predicate over (baseline_cmds, new_cmds, baseline_rules, new_rules)
SELECTORS = {
    ("consensus", "tools"): (
        ["ciriquant"],
        lambda b, n: "--dcc" in " ".join(b["consensus_filter"]) and
        "--dcc" not in " ".join(n["consensus_filter"]) and not n.get("dcc")),
    ("consensus", "adaptive"): (
        False,
        lambda b, n: "--adaptive " in " ".join(b["consensus_filter"]) + " " and
        not re.search(r"--adaptive(\s|$)", " ".join(n["consensus_filter"]))),
    ("de", "tumor_label"): ("tumor_X", lambda b, n: set(n) != set(b) or n != b),
    ("de", "normal_label"): ("normal_X", lambda b, n: set(n) != set(b) or n != b),
}

# Config keys that are consumed through the Snakefile's module-level aliases.
ALIASES = {"DE_METHOD": ("de", "method")}

# Keys we KNOW are not wired. Each entry needs a reason. strict xfail => if somebody wires one,
# the test flips to XPASS(strict) and fails until the entry is removed (so this list cannot rot).
KNOWN_GAPS: dict = {}


def _blocks(text: str) -> dict:
    """{rule_name: rule_text} for one .smk file."""
    out, name, buf = {}, None, []
    for line in text.splitlines():
        m = re.match(r"^rule (\w+):", line)
        if m:
            if name:
                out[name] = "\n".join(buf)
            name, buf = m.group(1), [line]
        elif name:
            buf.append(line)
    if name:
        out[name] = "\n".join(buf)
    return out


def _all_rules() -> dict:
    d = {}
    for f in RULES:
        d.update(_blocks(f.read_text()))
    return d


def _params_block(rule_text: str) -> str:
    m = re.search(r"^    params:\n((?:        .*\n|\n)+)", rule_text + "\n", re.M)
    return m.group(1) if m else ""


def _script_of(rule_text: str):
    m = re.search(r'^    script:\s*\n\s+"([^"]+)"', rule_text, re.M)
    return (ROOT / "workflow" / "rules" / m.group(1)).resolve() if m else None


def _config_refs(text: str) -> set:
    """(section, key) pairs read from config in a piece of Snakemake/Python source."""
    refs = set()
    for m in re.finditer(r'config\["(\w+)"\]\["(\w+)"\]', text):
        refs.add((m.group(1), m.group(2)))
    for m in re.finditer(r'config\["(\w+)"\]\.get\(\s*"(\w+)"', text):
        refs.add((m.group(1), m.group(2)))
    for m in re.finditer(r"config\.get\(\s*\"(\w+)\"\s*,\s*\{\}\s*\)\.get\(\s*\"(\w+)\"", text):
        refs.add((m.group(1), m.group(2)))
    for alias, key in ALIASES.items():
        if re.search(r"\b%s\b" % alias, text):
            refs.add(key)
    return refs


def _declared_defaults() -> dict:
    """Keys the rules read with ``.get(key, default)`` -> literal default (undeclared in config.yaml)."""
    out = {}
    text = "\n".join(f.read_text() for f in RULES)
    for m in re.finditer(r'config\["(\w+)"\]\.get\(\s*"(\w+)"\s*,\s*([^\)\n]+?)\s*\)', text):
        try:
            out[(m.group(1), m.group(2))] = ast.literal_eval(m.group(3))
        except (ValueError, SyntaxError):
            pass
    return out


def _candidate_keys():
    base = w.load_base_config()
    keys = {}
    for sect in SECTIONS:
        if sect in base and isinstance(base[sect], dict):
            for p in w.leaf_paths({sect: base[sect]}):
                keys[p] = w.get_path(base, p)
    for (sect, key), default in _declared_defaults().items():
        if sect in SECTIONS and (sect, key) not in keys:
            keys[(sect, key)] = default
    # list-valued keys are selectors, handled explicitly
    keys.setdefault(("consensus", "tools"), base.get("consensus", {}).get("tools", ["ciriquant", "dcc"]))
    return sorted(keys.items())


CANDIDATES = _candidate_keys()


@pytest.fixture(scope="module")
def fx(tmp_path_factory):
    return tmp_path_factory.mktemp("wiring_fx")


@pytest.fixture(scope="module")
def baseline(fx):
    return w.render(fx, w.load_base_config())


def test_dry_run_renders_expected_rules(baseline):
    """Sanity: the harness sees the shell commands we rely on."""
    for rule in ("consensus_filter", "rank_biomarkers", "ciriquant", "predict_interactions"):
        assert baseline.get(rule), "no shell command rendered for %s" % rule
    assert any("CIRIquant" in c for c in baseline["ciriquant"])


# ── 1+3. dynamic: each key must reach a command ───────────────────────────────────────────────

SCRIPT_RULES = {n: t for n, t in _all_rules().items() if _script_of(t)}
SCRIPT_REFS = {}
for _n, _t in SCRIPT_RULES.items():
    for _ref in _config_refs(_params_block(_t)):
        SCRIPT_REFS.setdefault(_ref, []).append(_n)


def _id(p):
    return ".".join(p[0])


@pytest.mark.parametrize("path,default", CANDIDATES, ids=[_id(c) for c in CANDIDATES])
def test_config_key_reaches_a_consumer(path, default, fx, baseline, request):
    dotted = ".".join(path)
    if dotted in KNOWN_GAPS:
        request.applymarker(pytest.mark.xfail(reason=KNOWN_GAPS[dotted], strict=True))

    cfg = copy.deepcopy(w.load_base_config())
    if path in SELECTORS:
        new, predicate = SELECTORS[path]
        w.set_path(cfg, path, new)
        after = w.render(fx, cfg)
        assert predicate(baseline, after), (
            "%s=%r did not change the job graph/commands as expected" % (dotted, new))
        return

    new = w.deviate(default)
    w.set_path(cfg, path, new)
    after = w.render(fx, cfg)

    hit_rules = [r for r, cmds in after.items()
                 if any(re.search(r"(?<![\w.])%s(?![\w.])" % re.escape(str(new)), c) for c in cmds)
                 and after[r] != baseline.get(r)]
    if hit_rules:
        return                                              # value reached a real command line
    if path in SCRIPT_REFS:
        return                                              # reaches a script rule's params (see 2.)
    pytest.fail(
        "SILENT WIRING GAP: %s was set to %r but appears in no rendered shell command and in no "
        "script-rule params. Either wire it into the rule or document it in KNOWN_GAPS." %
        (dotted, new))


# ── 2. static: script rules' params <-> script reads ──────────────────────────────────────────

def _defined_param_names(rule_text: str) -> set:
    return set(re.findall(r"^        (\w+)\s*=", _params_block(rule_text), re.M))


def _script_param_reads(script: Path) -> set:
    s = script.read_text()
    names = set(re.findall(r'snakemake@params\[\[\s*"(\w+)"\s*\]\]', s))
    names |= set(re.findall(r"snakemake@params\$(\w+)", s))
    names |= set(re.findall(r"snakemake\.params\.(\w+)", s))
    names |= set(re.findall(r'snakemake\.params\[\s*["\'](\w+)["\']\s*\]', s))
    names |= set(re.findall(r'(?:sn|snakemake)\.params,\s*["\'](\w+)["\']', s))
    return names


@pytest.mark.parametrize("rule", sorted(SCRIPT_RULES))
def test_script_rule_params_consumed_and_defined(rule):
    text = SCRIPT_RULES[rule]
    script = _script_of(text)
    assert script.exists(), "%s: script %s missing" % (rule, script)
    defined = _defined_param_names(text)
    read = _script_param_reads(script)
    never_read = sorted(defined - read)
    never_defined = sorted(read - defined)
    msgs = []
    if never_read:
        msgs.append("rule %s passes params the script never reads (dead wiring): %s" % (rule, never_read))
    if never_defined:
        msgs.append("script %s reads params the rule never defines (NULL/AttributeError at run time): %s"
                    % (script.name, never_defined))
    assert not msgs, "\n".join(msgs)


# ── sections the project believes exist ───────────────────────────────────────────────────────

EXPECTED_SECTIONS = ["de.biomarker", "tools.ciriquant", "rank_biomarkers_v2"]


@pytest.mark.xfail(reason="repo HEAD has none of de.biomarker.*, tools.ciriquant.*, rank_biomarkers_v2 "
                          "(components_arg / sig_cap / fc_cap / anchor live in an unpushed version). "
                          "This test XPASSes once they are merged; the loops above then cover them "
                          "automatically.", strict=False)
def test_expected_sections_present():
    base = w.load_base_config()
    rules = " ".join(_all_rules())
    missing = []
    for name in EXPECTED_SECTIONS:
        if name == "rank_biomarkers_v2":
            ok = "rank_biomarkers_v2" in rules or (ROOT / "scripts" / "rank_biomarkers_v2.py").exists()
        else:
            sect, sub = name.split(".")
            ok = isinstance(base.get(sect, {}).get(sub), dict)
        if not ok:
            missing.append(name)
    assert not missing, "missing in this checkout: %s" % missing
