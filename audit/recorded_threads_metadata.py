"""Read-only: the threads value recorded in each dataset's CIRIquant metadata (the job that produced the
published GTFs). Written to a JSON file for appendix_params.py."""
import base64
import glob
import json
import os
import re
import sys
from collections import Counter

REPO = os.path.expanduser("~/circRNA_agent")
OUT = sys.argv[1]
GSES = ["GSE108735", "GSE113230", "GSE121842", "GSE130078", "GSE133998", "GSE136569", "GSE143797",
        "GSE148036", "GSE171011", "GSE221107", "GSE248612", "GSE77509", "GSE97239", "PRJNA553289", "SRP156355",
        "GSE58135", "GSE192410", "GSE192849", "GSE229705"]
from collections import defaultdict
recorded = defaultdict(Counter)
for path in glob.glob(os.path.join(REPO, ".snakemake", "metadata", "*")):
    key = os.path.basename(path)
    try:
        out = base64.b64decode(key + "=" * (-len(key) % 4)).decode("utf-8", "replace")
    except Exception:
        continue
    m = re.search(r"/(GSE\d+|SRP\d+|PRJNA\d+)_results/circRNA/[^/]+/[^/]+\.gtf$", out)
    if not m:
        continue
    meta = json.load(open(path))
    if meta.get("rule") != "ciriquant":
        continue
    t = re.search(r"-t (\d+)", meta.get("shellcmd", ""))
    if t:
        recorded[m.group(1)][int(t.group(1))] += 1
result = {g: {str(k): v for k, v in recorded[g].items()} for g in GSES}
json.dump(result, open(OUT, "w"), indent=1)
print("wrote", OUT, {g: dict(c) for g, c in recorded.items() if c})
