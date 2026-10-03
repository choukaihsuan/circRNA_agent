"""Build audit/library_strandedness.tsv and audit/library_strandedness_mismatch.tsv from the GEO GSM pages in raw/.
Judgement rules (applied to the GSM-level !Sample_extract_protocol_ch1 text only; nothing from the paper or CLAUDE.md):
  stranded    : text explicitly says strand-specific / stranded library|kit / dUTP second-strand marking, and not 'non-stranded'
  unstranded  : text explicitly says non-stranded/unstranded  -> tier 'explicit'
                OR cDNA made with random primers and no strand-marking step described -> tier 'text-implied' (precedent: GSE97239)
  未載明      : anything else (incl. 'standard Illumina protocols', 'none provided', SRA-only submissions with no GSM)
Kit names that merely *imply* a mode from product knowledge (e.g. NEBNext Ultra RNA, Trio/SPIA) are NOT counted; a note is recorded.
"""
import csv, re, os, collections
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "..")
sg = list(csv.DictReader(open(os.path.join(HERE, "srr_gsm.tsv")), delimiter="\t"))
order = []; byds = collections.defaultdict(list)
for r in sg:
    if r["dataset"] not in byds: order.append(r["dataset"])
    byds[r["dataset"]].append(r)
EXCL = {"GSE58135", "GSE229705", "GSE192410", "GSE192849"}
STRANDED = re.compile(r"strand-specific|stranded (?:RNA|library|mRNA|Total)|TruSeq Stranded|KC\S* Stranded|dUTP", re.I)
NONSTR = re.compile(r"non-?stranded|unstranded", re.I)
RANDOM = re.compile(r"random (?:primer|hexamer)", re.I)
def sentence(text, m):
    s = max(text.rfind(".", 0, m.start()), text.rfind("||", 0, m.start())); e = text.find(".", m.end())
    return text[s+1 if s >= 0 else 0: e+1 if e >= 0 else len(text)].strip(" |")[:420]
rows = []; mism = []
for d in order:
    rs = byds[d]; gsms = sorted({r["gsm"] for r in rs if r["gsm"]})
    if not gsms:
        rows.append([d, "(none: SRA-only submission, no GEO series/GSM)", "SRA DESIGN_DESCRIPTION only (no extract_protocol exists): " + {"SRP156355": "rRNA depletion", "PRJNA553289": "Transcriptome sequencing"}.get(d, ""), "未載明", "未載明", "no GSM exists; SRA experiment XML (ENA browser API) has no library-construction/strand text", len(rs), "yes" if d in EXCL else "no"]); continue
    texts = set(); rep = gsms[0]
    for g in gsms:
        t = open(os.path.join(HERE, "raw", g + ".soft"), encoding="utf8", errors="ignore").read()
        texts.add(" || ".join(re.findall(r"^!Sample_extract_protocol_ch1 = (.*)$", t, re.M)))
    assert len(texts) == 1, (d, len(texts))
    text = texts.pop()
    ms, mn, mr = STRANDED.search(text), NONSTR.search(text), RANDOM.search(text)
    note = ""
    if mn: verdict, tier, snip, basis = "unstranded", "explicit", sentence(text, mn), "text states non-stranded"
    elif ms: verdict, tier, snip, basis = "stranded", "explicit", sentence(text, ms), "text states '%s'" % ms.group(0)
    elif mr: verdict, tier, snip, basis = "unstranded", "text-implied", sentence(text, mr), "cDNA synthesised with random primer; no dUTP/strand-marking step described (precedent GSE97239) - implied by absence, not stated"
    else:
        lib = re.search(r"librar|Ribo|TruSeq|NEBNext|Trio", text, re.I)
        snip = text[max(0, lib.start()-40) if lib else 0:][:300]; verdict, tier, basis = "未載明", "未載明", "protocol text does not mention strandedness"
        if d == "GSE148036": note = "kit named 'NEBNext Ultra RNA Library Prep Kit' (non-directional product per NEB catalogue) - not stated in text, not counted"
        if d == "GSE229705": note = "Trio RNA-Seq (SPIA) - text silent on strand; kit documentation not used as evidence"
    ex_note = ""
    if re.search(r"RNase ?R", text, re.I): ex_note = "also: text mentions RNase R digestion"
    rows.append([d, "%s (+%d more, identical text)" % (rep, len(gsms)-1) if len(gsms) > 1 else rep, snip, verdict, tier, "; ".join(x for x in (basis, note, ex_note) if x), len(rs), "yes" if d in EXCL else "no"])
hdr = ["dataset", "gsm", "protocol_snippet_verbatim", "verdict", "evidence_tier", "basis_and_notes", "n_samples_analysed", "excluded_from_15ds"]
with open(os.path.join(OUT, "library_strandedness.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow(hdr); w.writerows(rows)
V = {r[0]: r for r in rows}
with open(os.path.join(OUT, "library_strandedness_mismatch.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["list", "dataset", "verdict_of_real_library", "evidence_tier", "tool_setting_used", "excluded_from_15ds"])
    for d in order:
        r = V[d]
        if r[3] == "unstranded": w.writerow(["A: DCC run stranded (no -N) but library is unstranded", d, r[3], r[4], "DCC default = stranded (fr-firststrand)", r[7]])
    for d in order:
        r = V[d]
        if r[3] == "stranded": w.writerow(["B: CIRIquant run unstranded (-l 0) but library is stranded", d, r[3], r[4], "CIRIquant default -l 0", r[7]])
    for d in order:
        r = V[d]
        if r[3] == "未載明": w.writerow(["C: strandedness 未載明 (cannot be checked)", d, r[3], r[4], "-", r[7]])
print(open(os.path.join(OUT, "library_strandedness_mismatch.tsv")).read())
