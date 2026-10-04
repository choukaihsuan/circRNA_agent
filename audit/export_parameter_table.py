"""Parameter table for the thesis appendix: code default vs value used in the published runs.

Columns: kind | name | code_default | published_value | status | help | defined_in

  code_default     argparse default (scripts/*.py), the R fallback (scripts/analysis.R),
                   the .get(key, default) fallback in workflow/rules/*.smk, or the value
                   in config.yaml. This is the effective value when a key is absent.
  published_value  effective value across the 20 published per-project configs
                   (config/projects/*.yaml). "same" if every project agrees, otherwise
                   "varies: v (n); ..." listing each distinct value and its project count.
  status           same / DIFF / n/a (not a config-driven value)

Also writes audit/parameter_discrepancies.tsv with every DIFF row.

usage:
  python3 audit/export_parameter_table.py <projects_dir>
where <projects_dir> holds the 20 config/projects/*.yaml files copied from the HPC server.
"""
import ast
import collections
import glob
import os
import re
import sys

import yaml

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
PROJ_DIR = sys.argv[1] if len(sys.argv) > 1 else ""
if not PROJ_DIR or not glob.glob(os.path.join(PROJ_DIR, "*.yaml")):
    raise SystemExit("usage: export_parameter_table.py <projects_dir with *.yaml>")
rows = []


def rel(p):
    return os.path.relpath(p, ROOT).replace(os.sep, "/")


def flat(d, pre=""):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(flat(v, f"{pre}.{k}" if pre else str(k)))
    else:
        out[pre] = d
    return out


projects = {os.path.basename(f): flat(yaml.safe_load(open(f, encoding="utf-8")) or {})
            for f in sorted(glob.glob(os.path.join(PROJ_DIR, "*.yaml")))}


def published(path, code_default):
    """Effective value across projects: the project's value if present, else the code default."""
    vals = collections.Counter()
    for pf in projects.values():
        v = pf[path] if path in pf else code_default
        vals[str(v)] += 1
    if len(vals) == 1:
        return next(iter(vals)), "same"
    return "varies: " + "; ".join(f"{v} ({n})" for v, n in vals.most_common()), "DIFF"


# 1. argparse defaults --------------------------------------------------------
argparse_rows = []
for path in sorted(glob.glob(os.path.join(ROOT, "scripts", "*.py"))):
    try:
        tree = ast.parse(open(path, encoding="utf-8").read())
    except SyntaxError:
        continue
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "add_argument"):
            continue
        names = [a.value for a in node.args if isinstance(a, ast.Constant) and isinstance(a.value, str)]
        if not names:
            continue
        opt = [n for n in names if n.startswith("--")] or names
        kw = {k.arg: k.value for k in node.keywords if k.arg}

        def _val(n):
            if n is None:
                return None
            try:
                return ast.literal_eval(n)
            except Exception:
                return "(not static)"
        d = _val(kw.get("default"))
        if kw.get("required") is not None and _val(kw["required"]):
            default = "(required)"
        else:
            default = "(not given)" if d is None else str(d)
        argparse_rows.append(dict(flag=opt[0], code_default=default,
                                  help=str(_val(kw.get("help")) or "").replace("\n", " "),
                                  source=rel(path), line=node.lineno))

# 2. Snakemake params blocks: name -> (config path, .get fallback, file) ------
smk_map = {}
smk_files = sorted(glob.glob(os.path.join(ROOT, "workflow", "rules", "*.smk")))
smk_files.append(os.path.join(ROOT, "workflow", "Snakefile"))
for smk in smk_files:
    if not os.path.exists(smk):
        continue
    txt = open(smk, encoding="utf-8").read()
    for m in re.finditer(r"^\s+(\w+)\s*=\s*(config\[.*?)(?:,\s*$|$)", txt, re.M):
        name, expr = m.group(1), m.group(2).strip()
        keys = re.findall(r'\[\s*"([^"]+)"\s*\]', expr)  # config["a"]["b"] chain
        fb = re.search(r'\.get\(\s*"([^"]+)"\s*,\s*([^)]+)\)', expr)
        if fb:  # config["de"].get("de_sig_by", default): the key is the .get argument
            keys = keys + [fb.group(1)]
        path = ".".join(keys) if keys else None
        fallback = fb.group(2).strip().strip('"') if fb else None
        smk_map.setdefault(name, (path, fallback, rel(smk)))

# 3. config.yaml (parsed as YAML so lists and nested keys are read correctly) ----
cfg_path = os.path.join(ROOT, "config.yaml")
_cfg_text = open(cfg_path, encoding="utf-8").read().splitlines()
cfg_rows = []
for dotted, val in flat(yaml.safe_load(open(cfg_path, encoding="utf-8")) or {}).items():
    line_no = next((i for i, l in enumerate(_cfg_text, 1)
                    if re.match(r"^\s*" + re.escape(dotted.split(".")[-1]) + r"\s*:", l)), 0)
    comment = ""
    m = re.search(r"#\s*(.*)$", _cfg_text[line_no - 1]) if line_no else None
    if m:
        comment = m.group(1).strip()
    cfg_rows.append((dotted, str(val) if not isinstance(val, list) else str(val).replace("'", ""), comment, line_no))
# Snakefile fallback applied when a key is absent from config.yaml
SNAKEFILE_FALLBACK = {"consensus.tools": ("['ciriquant', 'dcc']", "workflow/Snakefile:80")}
cfg_default = {d: v for d, v, _, _ in cfg_rows}
_proj_only = sorted({k for pf in projects.values() for k in pf} - set(cfg_default))

# 4. R fallbacks in analysis.R -------------------------------------------------
r_path = os.path.join(ROOT, "scripts", "analysis.R")
r_src = open(r_path, encoding="utf-8").read().splitlines()
r_rows = {}
for i, line in enumerate(r_src, 1):
    for m in re.finditer(r'snakemake@params\[\["([^"]+)"\]\]', line):
        name = m.group(1)
        if name in r_rows:
            continue
        window = "\n".join(r_src[max(0, i - 1): i + 6])
        d = re.search(r'<-\s*(as\.[a-z.]+\()?\s*("[^"]*"|[0-9.]+L?|TRUE|FALSE)', window)
        r_rows[name] = (d.group(2).strip('"') if d else None, i)

# 5. build rows ----------------------------------------------------------------
for k in _proj_only:
    code_default = SNAKEFILE_FALLBACK[k][0] if k in SNAKEFILE_FALLBACK else "(not in config.yaml)"
    pub, status = published(k, code_default)
    rows.append(dict(kind="config", name="config." + k, code_default=code_default, published=pub,
                     status=status, help="key present in some published projects only" +
                     (f"; code fallback at {SNAKEFILE_FALLBACK[k][1]}" if k in SNAKEFILE_FALLBACK else ""),
                     source="config.yaml", line=0))
for dotted, val, comment, line_no in cfg_rows:
    pub, status = published(dotted, val)
    rows.append(dict(kind="config", name="config." + dotted, code_default=val, published=pub,
                     status=status, help=comment, source="config.yaml", line=line_no))

for r in argparse_rows:
    passed_by = None
    for smk in smk_files:
        if not os.path.exists(smk):
            continue
        mm = re.search(re.escape(r["flag"]) + r"\s+\{params\.(\w+)\}", open(smk, encoding="utf-8").read())
        if mm:
            passed_by = mm.group(1)
            break
    if passed_by and smk_map.get(passed_by, (None,))[0]:
        path, fb, _ = smk_map[passed_by]
        pub, status = published(path, cfg_default.get(path, fb))
        rows.append(dict(kind="argparse", name=r["flag"], code_default=r["code_default"],
                         published=pub + f" [config.{path}]", status=status,
                         help=r["help"], source=r["source"], line=r["line"]))
    else:
        rows.append(dict(kind="argparse", name=r["flag"], code_default=r["code_default"],
                         published="not passed by pipeline (code default)", status="n/a",
                         help=r["help"], source=r["source"], line=r["line"]))

smk_map["de_method"] = ("de.method", None, "workflow/Snakefile (DE_METHOD)")
for name, (rv, rl) in sorted(r_rows.items()):
    path, fb, smk_file = smk_map.get(name, (None, None, None))
    if path:
        # the rule always passes the config value, so the repo config is the code default;
        # analysis.R's own fallback applies only if the parameter is missing.
        code = cfg_default.get(path, fb if fb is not None else rv)
        pub, status = published(path, code)
        fb_note = f"; analysis.R fallback {rv}" if rv is not None and str(rv) != str(code) else ""
        rows.append(dict(kind="R", name=name, code_default=str(code), published=pub + f" [config.{path}]",
                         status=status, help="analysis.R parameter; passed from " + str(smk_file) + fb_note,
                         source=rel(r_path), line=rl))
    else:
        code = rv if rv is not None else "(not static)"
        rows.append(dict(kind="R", name=name, code_default=str(code), published="not from config",
                         status="n/a", help="analysis.R parameter", source=rel(r_path), line=rl))

IDENT = re.compile(r"^config\.(project_id|raw_dir|trimmed_dir|results_dir|groups|metadata|log_dir|study_title|"
                   r"download\.sra_cache_dir|download\.tmp_dir|notify\..*)$")
REF = re.compile(r"^config\.(genome\..*|circbase_file)$")
COMP = re.compile(r"^config\.(threads|download\.max_concurrent)$")
for r in rows:
    n = r["name"]
    if IDENT.match(n):
        r["category"] = "identifier/path"
    elif REF.match(n):
        r["category"] = "reference"
    elif COMP.match(n):
        r["category"] = "computational"
    else:
        r["category"] = "analysis"
    if n.startswith("config.notify.") and "email" in n:
        r["published"] = "<redacted: personal address>" if r["status"] == "DIFF" else r["published"]
        r["code_default"] = "<redacted>" if r["code_default"] else r["code_default"]
for r in rows:
    if r["status"] == "DIFF" and r["category"] == "analysis" and r["name"] == "config.de.method":
        pass

# 6. write ----------------------------------------------------------------------
rows.sort(key=lambda r: (r["kind"], r["name"]))
with open(os.path.join(ROOT, "audit", "parameter_table.tsv"), "w", encoding="utf-8") as fh:
    fh.write("kind\tcategory\tname\tcode_default\tpublished_value\tstatus\thelp\tdefined_in\tline\n")
    for r in rows:
        fh.write(f'{r["kind"]}\t{r["category"]}\t{r["name"]}\t{r["code_default"]}\t{r["published"]}\t{r["status"]}\t'
                 f'{r["help"]}\t{r["source"]}\t{r["line"]}\n')
with open(os.path.join(ROOT, "audit", "parameter_table.md"), "w", encoding="utf-8") as fh:
    fh.write("| 參數名 | 程式碼預設值 | 已發表執行使用值 | 狀態 | 作用 | 定義位置 |\n|---|---|---|---|---|---|\n")
    for r in rows:
        fh.write(f'| `{r["name"]}` | {r["code_default"]} | {r["published"]} | {r["status"]} | '
                 f'{r["help"]} | {r["source"]}:{r["line"]} |\n')
def _dup_of_config(r):
    """argparse/R rows that just restate a config-driven value are listed once, via the config row."""
    return r["kind"] != "config" and r["published"].endswith("]") and "[config." in r["published"]


disc = [r for r in rows if r["status"] == "DIFF" and not _dup_of_config(r)]
order = {"analysis": 0, "reference": 1, "computational": 2, "identifier/path": 3}
disc.sort(key=lambda r: (order.get(r["category"], 9), r["name"]))
with open(os.path.join(ROOT, "audit", "parameter_discrepancies.tsv"), "w", encoding="utf-8") as fh:
    fh.write("category\tkind\tname\tcode_default\tpublished_value\thelp\tdefined_in\tline\n")
    for r in disc:
        fh.write(f'{r["category"]}\t{r["kind"]}\t{r["name"]}\t{r["code_default"]}\t{r["published"]}\t{r["help"]}\t'
                 f'{r["source"]}\t{r["line"]}\n')
print(f"{len(rows)} parameters; {len(disc)} DIFF rows -> audit/parameter_table.tsv|md, audit/parameter_discrepancies.tsv")
