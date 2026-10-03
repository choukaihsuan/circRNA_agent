"""GSE148036: ENA sample_alias is a sample name (AD001...), not a GSM. Map SRR -> SRX (ENA) -> GSM (GEO series SOFT 'relation = SRA: ...SRX').
Appends to srr_gsm.tsv in place and downloads the GSM pages. Run after fetch_geo_protocols.py."""
import csv, re, os, urllib.request, urllib.parse, time, gzip, io
HERE = os.path.dirname(os.path.abspath(__file__))
def get(u):
    return urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "circDEX-audit"}), timeout=90).read()
rows = list(csv.reader(open(os.path.join(HERE, "srr_gsm.tsv")), delimiter="\t")); hdr, body = rows[0], rows[1:]
srrs = [r[1] for r in body if r[0] == "GSE148036"]
q = " OR ".join("run_accession=" + s for s in srrs)
ena = get("https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv&fields=run_accession,experiment_accession&query=" + urllib.parse.quote(q)).decode()
srr2srx = {l.split("\t")[0]: l.split("\t")[1] for l in ena.splitlines()[1:]}
soft = gzip.decompress(get("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148036/soft/GSE148036_family.soft.gz")).decode("utf8", "replace")
srx2gsm = {}
for blk in soft.split("^SAMPLE = ")[1:]:
    gsm = blk.split("\n")[0].strip()
    m = re.search(r"relation = SRA: .*?term=(SRX\d+)", blk)
    if m: srx2gsm[m.group(1)] = gsm
for r in body:
    if r[0] == "GSE148036": r[2] = srx2gsm.get(srr2srx.get(r[1], ""), "")
with open(os.path.join(HERE, "srr_gsm.tsv"), "w") as f:
    for r in [hdr] + body: f.write("\t".join(r) + "\n")
for g in sorted({r[2] for r in body if r[0] == "GSE148036" and r[2]}):
    p = os.path.join(HERE, "raw", g + ".soft")
    if not os.path.exists(p):
        open(p, "w").write(get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=full" % g).decode("utf8", "replace")); time.sleep(0.4)
print({r[1]: r[2] for r in body if r[0] == "GSE148036"})
