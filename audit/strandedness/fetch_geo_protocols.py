"""Fetch GSM-level extract_protocol text for every sample actually analysed (19 datasets).
SRR -> GSM via ENA portal API (sample_alias / experiment_alias); GSM text via GEO accession page (SOFT, view=full).
Writes raw/<GSM>.soft (GEO text, kept for provenance) and srr_gsm.tsv.  Needs network; read-only.
Run: python3 fetch_geo_protocols.py
"""
import os, re, sys, time, urllib.request, urllib.parse, csv
HERE = os.path.dirname(os.path.abspath(__file__)); RAW = os.path.join(HERE, "raw"); os.makedirs(RAW, exist_ok=True)
SAMPLES = {}
for line in open(os.path.join(HERE, "sample_srr_list.txt")):
    a = line.split()
    if a: SAMPLES[a[0]] = a[1:]
def get(url, tries=4):
    for i in range(tries):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "circDEX-audit"}), timeout=60).read().decode("utf8", "replace")
        except Exception as e:
            err = e; time.sleep(3*(i+1))
    raise err
rows = []
for d, srrs in SAMPLES.items():
    q = " OR ".join("run_accession=%s" % s for s in srrs)
    url = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv&fields=run_accession,sample_alias,experiment_alias,sample_title&query=" + urllib.parse.quote(q)
    for r in csv.DictReader(get(url).splitlines(), delimiter="\t"):
        gsm = next((x for x in (r["sample_alias"], r["experiment_alias"], r["sample_title"]) if re.match(r"GSM\d+", x or "")), "")
        gsm = re.match(r"GSM\d+", gsm).group(0) if gsm else ""
        rows.append((d, r["run_accession"], gsm, r["sample_alias"], r["experiment_alias"]))
    print(d, len(srrs), "samples queried", file=sys.stderr)
with open(os.path.join(HERE, "srr_gsm.tsv"), "w") as f:
    f.write("dataset\tsrr\tgsm\tena_sample_alias\tena_experiment_alias\n")
    for r in rows: f.write("\t".join(r) + "\n")
for gsm in sorted({r[2] for r in rows if r[2]}):
    p = os.path.join(RAW, gsm + ".soft")
    if os.path.exists(p): continue
    t = get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=full" % gsm)
    open(p, "w").write(t); time.sleep(0.4)
print("done", len({r[2] for r in rows if r[2]}), "GSMs")
