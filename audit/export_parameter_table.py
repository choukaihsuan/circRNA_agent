"""Export every tunable parameter of the pipeline from the code (not from memory).

Sources:
  1. argparse options in scripts/*.py        (name, default, help, file:line)
  2. R defaults in scripts/analysis.R        (snakemake@params[[...]] with a fallback default)
  3. keys and values in config.yaml          (top-level and nested, with file:line)

Output: audit/parameter_table.tsv (one row per parameter and source) and the same rows as
markdown in audit/parameter_table.md. Rows whose default cannot be read statically are kept
with default = "(not static)" rather than dropped.

usage: python3 audit/export_parameter_table.py   (run from the repo root)
"""
import ast
import glob
import os
import re

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
rows = []


def rel(p):
    return os.path.relpath(p, ROOT).replace(os.sep, "/")


# 1. argparse ----------------------------------------------------------------
for path in sorted(glob.glob(os.path.join(ROOT, "scripts", "*.py"))):
    src = open(path, encoding="utf-8").read()
    try:
        tree = ast.parse(src)
    except SyntaxError as e:
        rows.append(dict(name="(parse error)", default="", help=str(e), source=rel(path), line=0,
                         kind="argparse"))
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
        default = _val(kw.get("default"))
        help_ = _val(kw.get("help")) or ""
        rows.append(dict(name=opt[0], default=("(required)" if "required" in kw and _val(kw["required"]) else
                                             ("(not given)" if default is None else str(default))),
                         help=str(help_).replace("\n", " "), source=rel(path), line=node.lineno,
                         kind="argparse"))

# 2. R defaults --------------------------------------------------------------
r_path = os.path.join(ROOT, "scripts", "analysis.R")
r_src = open(r_path, encoding="utf-8").read().splitlines()
seen = set()
for i, line in enumerate(r_src, 1):
    for m in re.finditer(r'snakemake@params\[\["([^"]+)"\]\]', line):
        name = m.group(1)
        if name in seen:
            continue
        seen.add(name)
        # default: a later assignment like `if (...) x <- 5L` or `<- "..."` near the read
        default = "(not static)"
        window = "\n".join(r_src[max(0, i - 1): i + 6])
        d = re.search(r'<-\s*(as\.[a-z.]+\()?\s*("[^"]*"|[0-9.]+L?|TRUE|FALSE)', window)
        if d:
            default = d.group(2)
        rows.append(dict(name=name, default=default, help="snakemake params read in analysis.R",
                         source=rel(r_path), line=i, kind="R"))
# bare R constants that define a default
for i, line in enumerate(r_src, 1):
    m = re.match(r'\s*(\w+)\s*<-\s*(\d+(?:\.\d+)?L?)\s*(#.*)?$', line)
    if m and m.group(1).endswith(("_cutoff", "_thr", "_min", "_count")):
        rows.append(dict(name=m.group(1), default=m.group(2), help=(m.group(3) or "").lstrip("# "),
                         source=rel(r_path), line=i, kind="R-constant"))

# 3. config.yaml -------------------------------------------------------------
cfg_path = os.path.join(ROOT, "config.yaml")
stack = []
for i, raw in enumerate(open(cfg_path, encoding="utf-8").read().splitlines(), 1):
    if not raw.strip() or raw.lstrip().startswith("#"):
        continue
    m = re.match(r"^(\s*)([A-Za-z_][\w]*)\s*:\s*(.*)$", raw)
    if not m:
        continue
    indent, key, val = len(m.group(1)), m.group(2), m.group(3)
    while stack and stack[-1][0] >= indent:
        stack.pop()
    path_ = ".".join([k for _, k in stack] + [key])
    comment = ""
    if "#" in val:
        val, comment = val.split("#", 1)
        comment = comment.strip()
    val = val.strip()
    if val == "":
        stack.append((indent, key))
        continue
    rows.append(dict(name="config." + path_, default=val.strip('"'), help=comment,
                     source=rel(cfg_path), line=i, kind="config"))

# 4. R parameters that the Snakemake rule passes from config (fill "(not static)") --------
R_FROM_CONFIG = {"fdr": "config.de.fdr_cutoff", "lfc": "config.de.log2fc_cutoff",
                 "de_sig_by": "config.de.de_sig_by", "tumor_label": "config.de.tumor_label",
                 "normal_label": "config.de.normal_label", "heatmap_top_n": "config.de.heatmap_top_n",
                 "fsj_concordance_lfc": "config.de.fsj_concordance_lfc"}
cfg_vals = {r["name"]: r["default"] for r in rows if r["kind"] == "config"}
for r in rows:
    if r["kind"] == "R" and r["default"] == "(not static)" and r["name"] in R_FROM_CONFIG:
        key = R_FROM_CONFIG[r["name"]]
        if key in cfg_vals:
            r["default"] = f"{cfg_vals[key]} (from {key}; code fallback in rule)"
            r["help"] = (r["help"] + "; passed by workflow/rules/de.smk").strip("; ")

# discrepancy: config.yaml default vs the value in the published per-project configs.
# config/projects/*.yaml exists on the HPC server only (not tracked here), so the observation
# is recorded with its date and source rather than re-read at export time.
rows.append(dict(name="config.de.de_sig_by (published projects)", default="pvalue (20 of 20 files)",
                 help="observed 2026-10-04 on server: grep de_sig_by ~/circRNA_agent/config/projects/*.yaml; "
                      "config.yaml default is auto -- the two disagree",
                 source="config/projects/ (server)", line=0, kind="audit-note"))

# write ----------------------------------------------------------------------
rows.sort(key=lambda r: (r["kind"], r["name"]))
tsv = os.path.join(ROOT, "audit", "parameter_table.tsv")
with open(tsv, "w", encoding="utf-8") as fh:
    fh.write("kind\tname\tdefault\thelp\tdefined_in\tline\n")
    for r in rows:
        fh.write(f'{r["kind"]}\t{r["name"]}\t{r["default"]}\t{r["help"]}\t{r["source"]}\t{r["line"]}\n')
md = os.path.join(ROOT, "audit", "parameter_table.md")
with open(md, "w", encoding="utf-8") as fh:
    fh.write("| 來源 | 參數名 | 預設值 | 作用 | 定義位置 |\n|---|---|---|---|---|\n")
    for r in rows:
        fh.write(f'| {r["kind"]} | `{r["name"]}` | {r["default"]} | {r["help"]} | {r["source"]}:{r["line"]} |\n')
print(f"{len(rows)} parameters -> {rel(tsv)}, {rel(md)}")
