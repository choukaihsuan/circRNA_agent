"""Shared loaders/metadata for thesis cross-dataset numbers and figures.
Inputs (all in ./data, copied read-only from the HPC server 2026-10-03):
  <GSE>.de.tsv.gz        circ_id,log2FC,pvalue,Type  (columns of <GSE>_results/de/de_results.tsv, edgeR_ciriquant)
  gene_map.tsv.gz        circ_id -> gene_name (isoform_groups.tsv, all 19)
  per_dataset_summary.tsv  (audit/ file: n_tumor/n_normal, single_tool_rows_pct, library_prep ...)
Significance (thesis Sec 2.7): pvalue < 0.05 and |log2FC| > 1.
"""
import csv, gzip, os, collections
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "data")
# paper order of Table/figure grouping: alphabetical by organ, then cancer type
META = {  # gse: (organ [paper naming], cancer_type label [paper abbreviation], lib_class)
 "GSE113230":  ("breast","TNBC","rRNA depletion"),
 "GSE58135":   ("breast","breast cancer, bulk","poly(A) selection"),
 "GSE133998":  ("breast","breast cancer, bulk","rRNA depletion"),
 "SRP156355":  ("breast","IDC","not stated"),
 "GSE192849":  ("breast","node-positive breast cancer","RNase R"),
 "GSE77509":   ("liver","HCC","rRNA depletion"),
 "GSE130078":  ("esophagus","ESCC","rRNA depletion"),
 "GSE248612":  ("stomach","gastric cancer","not stated"),
 "GSE221107":  ("prostate","prostate cancer","rRNA depletion"),
 "PRJNA553289":("lung","SCLC","not stated"),
 "GSE229705":  ("lung","LUAD","SPIA/Trio rRNA-transcript removal"),
 "GSE148036":  ("lung","LUAD","rRNA depletion"),
 "GSE121842":  ("colorectal","CRC","not stated"),
 "GSE136569":  ("pancreas","PDAC","rRNA depletion + RNase R"),
 "GSE143797":  ("nasopharynx","NPC","ncRNA-seq size fractionation"),
 "GSE108735":  ("kidney","RCC","ncRNA-seq size fractionation"),
 "GSE171011":  ("thyroid","PTC","not stated"),
 "GSE97239":   ("bladder","bladder cancer","rRNA depletion + RNase R"),
 "GSE192410":  ("ovary","ovarian cancer","not stated"),
}
DS19 = ["GSE113230","GSE58135","GSE133998","SRP156355","GSE77509","GSE130078","GSE248612","GSE221107","PRJNA553289","GSE229705","GSE148036","GSE121842","GSE136569","GSE143797","GSE108735","GSE171011","GSE97239","GSE192410","GSE192849"]
EXCL = ["GSE58135","GSE229705","GSE192410","GSE192849"]
DS15 = [d for d in DS19 if d not in EXCL]

def _open(p): return gzip.open(p,"rt")
def load():
    summ = {r["accession"]: r for r in csv.DictReader(open(os.path.join(DATA,"per_dataset_summary.tsv")),delimiter="\t")}
    sig, tested, genes = {}, {}, {}
    for d in DS19:
        S = {}; n = 0
        for r in csv.DictReader(_open(os.path.join(DATA,d+".de.tsv.gz")),delimiter="\t"):
            n += 1; p = float(r["pvalue"]); l = float(r["log2FC"])
            if p < 0.05 and abs(l) > 1: S[r["circ_id"]] = (l, r["Type"])
        sig[d] = S; tested[d] = n
    for line in _open(os.path.join(DATA,"gene_map.tsv.gz")):
        a = line.rstrip("\n").split("\t")
        if len(a) == 2 and a[1]: genes[a[0]] = a[1]
    return summ, sig, tested, genes

def tiers(sig, ds):
    c = collections.Counter()
    for d in ds:
        for k in sig[d]: c[k] += 1
    return c
