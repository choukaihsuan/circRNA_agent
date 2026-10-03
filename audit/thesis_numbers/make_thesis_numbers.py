"""Recompute every cross-dataset number used in the thesis for the 19-dataset and 15-dataset (dual-tool) sets.
Outputs: ../thesis_numbers_15ds.tsv, ../thesis_dataset_table_19ds.tsv
Run: python3 make_thesis_numbers.py   (pure stdlib; inputs in ./data)
"""
import csv, collections, os
from thesis_common import *
summ, sig, tested, genes = load()
OUT = os.path.join(HERE, "..")
rows = []
def add(item, v19, v15, src, loc=""): rows.append([item, v19, v15, src, loc])
SRC_DE = "data/<GSE>.de.tsv.gz (=<GSE>_results/de/de_results.tsv; sig = pvalue<0.05 & |log2FC|>1)"
def nsig(d): return len(sig[d])
def typ(d):
    c = collections.Counter(t for _, t in sig[d].values()); return c
def pct(a, b): return "%.1f" % (100*a/b) if b else "NA"
def rng(vals, f="%.1f"): return f % min(vals) + "–" + f % max(vals)

# ---- totals
t19 = sum(nsig(d) for d in DS19); t15 = sum(nsig(d) for d in DS15)
add("total significant calls (sum over datasets)", t19, t15, SRC_DE, "Abstract; §3.4 (opening paragraph)")
add("dataset with most calls", "GSE192410=%d (%s%% of total)" % (nsig("GSE192410"), pct(nsig("GSE192410"), t19)),
    "%s=%d (%s%% of total)" % (max(DS15, key=nsig), nsig(max(DS15, key=nsig)), pct(nsig(max(DS15, key=nsig)), t15)), SRC_DE, "§3.4 (opening paragraph): GSE192410 share of all calls")
add("range of calls per dataset (min–max)", "%d–%d" % (min(map(nsig, DS19)), max(map(nsig, DS19))), "%d–%d" % (min(map(nsig, DS15)), max(map(nsig, DS15))), SRC_DE, "§3.4 (opening paragraph): range of calls per dataset")
add("dataset with fewest calls", "GSE229705=1", "%s=%d" % (min(DS15, key=nsig), nsig(min(DS15, key=nsig))), SRC_DE, "§3.4 (opening paragraph)")
add("n datasets", len(DS19), len(DS15), "thesis_common.py DS19/DS15", "throughout")

# ---- Type I/II
def tstats(ds):
    p = []; T1 = T2 = 0; other = 0
    for d in ds:
        c = typ(d); T1 += c["Type_I"]; T2 += c["Type_II"]; other += sum(v for k, v in c.items() if k not in ("Type_I", "Type_II"))
        p.append(100*c["Type_I"]/nsig(d))
    return p, T1, T2, other
p19, a1, a2, o19 = tstats(DS19); p15, b1, b2, o15 = tstats(DS15)
add("Type I share per dataset, range (%)", rng(p19), rng(p15), SRC_DE, "Abstract; §3.4.1; §4.1, §4.5; §5")
lo19 = min(DS19, key=lambda d: 100*typ(d)["Type_I"]/nsig(d)); lo15 = min(DS15, key=lambda d: 100*typ(d)["Type_I"]/nsig(d))
add("dataset with lowest Type I share", "%s (%.1f%%)" % (lo19, min(p19)), "%s (%.1f%%)" % (lo15, min(p15)), SRC_DE, "")
h19 = [d for d in DS19 if typ(d)["Type_II"] == 0]; h15 = [d for d in DS15 if typ(d)["Type_II"] == 0]
add("datasets with 100% Type I", "%d (%s)" % (len(h19), ";".join(h19)), "%d (%s)" % (len(h15), ";".join(h15)), SRC_DE,
    "§3.4.1 refers to two datasets with 100% Type I; this becomes ONE dataset in the 15 set")
add("pooled Type I (n, %)", "%d/%d (%s%%)" % (a1, t19, pct(a1, t19)), "%d/%d (%s%%)" % (b1, t15, pct(b1, t15)), SRC_DE, "")
add("pooled Type II (n, %)", "%d (%s%%)" % (a2, pct(a2, t19)), "%d (%s%%)" % (b2, pct(b2, t15)), SRC_DE, "")
add("sig calls with Type other than I/II (sanity, should be 0)", o19, o15, SRC_DE, "")
add("datasets where Type I > Type II", "%d/19" % sum(typ(d)["Type_I"] > typ(d)["Type_II"] for d in DS19), "%d/15" % sum(typ(d)["Type_I"] > typ(d)["Type_II"] for d in DS15), SRC_DE, "§3.4.1")
for d in DS19:
    c = typ(d); n = nsig(d)
    v = "TypeI=%d (%s%%); TypeII=%d (%s%%); n=%d" % (c["Type_I"], pct(c["Type_I"], n), c["Type_II"], pct(c["Type_II"], n), n)
    add("Type I/II: " + d, v, v if d in DS15 else "EXCLUDED", SRC_DE, "Figure 7")

# ---- direction
def dirs(d):
    up = sum(1 for l, _ in sig[d].values() if l > 0); return up, nsig(d) - up
def dom(d):
    u, n = dirs(d); return "up" if u > n else ("down" if n > u else "tie")
for d in DS19:
    u, n = dirs(d); v = "up=%d; down=%d; dominant=%s (up share %s%%)" % (u, n, dom(d), pct(u, u+n))
    add("direction: " + d, v, v if d in DS15 else "EXCLUDED", SRC_DE, "Figure 8; §3.4.2")
cu19 = collections.Counter(dom(d) for d in DS19); cu15 = collections.Counter(dom(d) for d in DS15)
add("datasets up-dominant / down-dominant / tie", "%d / %d / %d" % (cu19["up"], cu19["down"], cu19["tie"]), "%d / %d / %d" % (cu15["up"], cu15["down"], cu15["tie"]), SRC_DE, "Supports direction-heterogeneity statements (§3.4.2; §5)")
add("up-dominant datasets", ";".join(d for d in DS19 if dom(d) == "up"), ";".join(d for d in DS15 if dom(d) == "up"), SRC_DE, "")
add("down-dominant datasets", ";".join(d for d in DS19 if dom(d) == "down"), ";".join(d for d in DS15 if dom(d) == "down"), SRC_DE, "")
def strong(ds, th):  # >=th share in one direction
    return sum(1 for d in ds if max(dirs(d))/nsig(d) >= th)
add("datasets with >=80% of calls in one direction", "%d/19" % strong(DS19, .8), "%d/15" % strong(DS15, .8), SRC_DE, "")
add("pooled up / down calls", "%d / %d" % (sum(dirs(d)[0] for d in DS19), sum(dirs(d)[1] for d in DS19)), "%d / %d" % (sum(dirs(d)[0] for d in DS15), sum(dirs(d)[1] for d in DS15)), SRC_DE, "pooled is dominated by GSE192410 (down) in the 19 set")
add("§3.4.2 named examples hold in 15 set?", "GSE192410 down-dominant; SRP156355 down; GSE113230 up; GSE108735 100% up",
    "GSE192410 EXCLUDED (sentence must change); SRP156355 %s; GSE113230 %s; GSE108735 up=%d down=%d" % (dom("SRP156355"), dom("GSE113230"), *dirs("GSE108735")), SRC_DE, "§3.4.2")

# ---- organ / cancer type
def org(ds): return sorted({META[d][0] for d in ds})
def ctype(ds): return sorted({META[d][1] for d in ds})
o19, o15 = org(DS19), org(DS15)
add("organ sites: n distinct (paper naming)", len(o19), len(o15), "thesis_common.META[...][0]", "Abstract; §3.4; Figure 8 caption; §4.1, §4.5; §5")
add("organ sites: list", ", ".join(o19), ", ".join(o15), "", "§3.4 (opening paragraph) lists 13")
add("organ sites removed by 15-set", "-", ", ".join(sorted(set(o19)-set(o15))), "", "ovary disappears (only GSE192410); breast keeps 3 sets, lung keeps 2")
cnt19 = collections.Counter(META[d][0] for d in DS19); cnt15 = collections.Counter(META[d][0] for d in DS15)
add("datasets per organ site", "; ".join("%s=%d" % kv for kv in sorted(cnt19.items())), "; ".join("%s=%d" % kv for kv in sorted(cnt15.items())), "", "")
add("organ sites with >=2 datasets (within-organ comparisons)", "; ".join(k for k, v in sorted(cnt19.items()) if v > 1), "; ".join(k for k, v in sorted(cnt15.items()) if v > 1), "", "web_ui 'Compare within' (breast, lung)")
c19, c15 = ctype(DS19), ctype(DS15)
add("cancer-type labels (histology/subtype as in paper abbreviations): n distinct", len(c19), len(c15), "thesis_common.META[...][1]", "NOT in v19 text as a number; likely source of the '15 種癌症／組織類型' wording")
add("cancer-type labels: list", "; ".join(c19), "; ".join(c15), "", "")
c15merged = sorted({("breast cancer" if META[d][0] == "breast" else META[d][1]) for d in DS15})
add("cancer-type count if all breast subtypes merged to 'breast cancer'", len({("breast cancer" if META[d][0]=="breast" else META[d][1]) for d in DS19}), len(c15merged), "", "alt. definition: gives 13 for the 15 set")
add("NOTE definition difference", "organ site = anatomical site (breast, lung…); cancer-type label = histology/subtype (TNBC, IDC, SCLC, LUAD…)",
    "15 set: 12 organ sites vs 15 cancer-type labels vs 15 datasets; the 15=15 equality is coincidental (breast has 3 labels TNBC/IDC/bulk, lung has 2 SCLC/LUAD, other 10 organs 1 each)", "", "choose ONE: '12 organ sites' for the 15 set; do not reuse '13'")

# ---- sample size
def npairs(d): return min(int(summ[d]["n_tumor"]), int(summ[d]["n_normal"]))
add("samples per arm (min–max pairs)", "%d–%d" % (min(map(npairs, DS19)), max(map(npairs, DS19))), "%d–%d" % (min(map(npairs, DS15)), max(map(npairs, DS15))), "data/per_dataset_summary.tsv n_tumor/n_normal", "§4.5; §4.10 (Limitations)")
sz19 = collections.Counter(npairs(d) for d in DS19); sz15 = collections.Counter(npairs(d) for d in DS15)
add("n per arm distribution", "; ".join("%d:%d" % kv for kv in sorted(sz19.items())), "; ".join("%d:%d" % kv for kv in sorted(sz15.items())), "", "")
add("designs that are NOT same-patient pairs (wording 'pairs' imprecise)", "GSE148036 (5T+5N unrelated), GSE143797 (4 NPC + 4 nasopharyngitis)", "same two (both in 15 set)", "CLAUDE.md dataset notes (unverified against GEO)", "v19 wording is 'pairs' (§4.5, §4.10)")

# ---- library prep
lib19 = collections.Counter(META[d][2] for d in DS19); lib15 = collections.Counter(META[d][2] for d in DS15)
doc19 = [k for k in lib19 if k != "not stated"]; doc15 = [k for k in lib15 if k != "not stated"]
add("library-prep: n distinct classes documented in GEO/SRA", len(doc19), len(doc15), "per_dataset_summary.tsv library_prep (GEO/SRA metadata) mapped by META[...][2]", "v19 gives no number for library-preparation protocols (§4.5, §4.6)")
add("library-prep: classes (datasets)", "; ".join("%s=%d" % kv for kv in sorted(lib19.items())), "; ".join("%s=%d" % kv for kv in sorted(lib15.items())), "", "")
add("library-prep: classes if 'not stated' counted as a class", len(lib19), len(lib15), "", "")
add("poly(A)-selected datasets included", "GSE58135", "none", "per_dataset_summary.tsv", "contradicts §2.13 exclusion in the 19 set")
add("RNase R in GSM extract_protocol text", "; ".join(d for d in DS19 if "RNase R" in META[d][2]), "; ".join(d for d in DS15 if "RNase R" in META[d][2]), "audit/library_strandedness.tsv basis column + GEO GSM extract_protocol (GSE136569 added 2026-10-03)", "§4.10 (currently says only GSE97239 among the 15)")
add("reads/sample range (M, both mates, fastp pre-filter)", rng([x for d in DS19 for x in map(float, summ[d]["reads_per_sample_M(both_mates,fastp_before_filtering)"].split("-"))], "%.0f"),
    rng([x for d in DS15 for x in map(float, summ[d]["reads_per_sample_M(both_mates,fastp_before_filtering)"].split("-"))], "%.0f"), "per_dataset_summary.tsv", "")
add("read length range (bp)", rng([x for d in DS19 for x in map(float, summ[d]["read_len_mean_range_bp"].split("-"))], "%.0f"),
    rng([x for d in DS15 for x in map(float, summ[d]["read_len_mean_range_bp"].split("-"))], "%.0f"), "per_dataset_summary.tsv", "50/51 bp sets are GSE58135/GSE229705 (both excluded)")

# ---- recurrence
c19c, c15c, c18c = tiers(sig, DS19), tiers(sig, DS15), tiers(sig, [d for d in DS19 if d != "GSE192410"])
def tier_row(c): return {k: sum(1 for v in c.values() if v >= k) for k in (1, 2, 3, 4, 5)}
r19, r15, r18 = tier_row(c19c), tier_row(c15c), tier_row(c18c)
for k in (1, 2, 3, 4, 5):
    add("recurrence: circRNAs significant in >=%d datasets" % k, r19[k], r15[k], SRC_DE, "Abstract; §3.4.3; §4.1, §4.6; Figure 9a/b" if k > 1 else "")
add("recurrence: exactly 5 datasets", sum(1 for v in c19c.values() if v == 5), sum(1 for v in c15c.values() if v == 5), SRC_DE, "§4.6")
add("recurrence: max n datasets", max(c19c.values()), max(c15c.values()), SRC_DE, "")
add("recurrence: unique circRNAs sig in >=1 dataset", len(c19c), len(c15c), SRC_DE, "")
add("sensitivity: excluding only GSE192410 (>=2/>=3/>=4)", "%d/%d/%d" % (r18[2], r18[3], r18[4]), "n/a", SRC_DE, "§4.10 (sensitivity analysis excluding GSE192410); reproduced here for 18 datasets")
add("recurrence: 15-set tiers as shares of 15-set unique", "-", "; ".join("%d:%s%%" % (k, pct(r15[k], r15[1])) for k in (2, 3, 4, 5)), SRC_DE, "")
# ---- seven candidates (section 4.6)
seven = [k for k, v in c19c.items() if v >= 5]
for k in sorted(seven, key=lambda x: (-c15c[x], x)):
    ds19 = [d for d in DS19 if k in sig[d]]; ds15 = [d for d in ds19 if d in DS15]
    add("§4.6 candidate %s (%s): n datasets" % (k, genes.get(k, "NA")), len(ds19), len(ds15), SRC_DE, "§4.6 seven candidates; 19-set members: " + ";".join(ds19) + (" | removed: " + ";".join(d for d in ds19 if d not in DS15) if len(ds19) != len(ds15) else ""))
add("§4.6 candidates still 5-of-5 in the 15 set", "7/7", "%d/7" % sum(1 for k in seven if c15c[k] >= 5), SRC_DE, "")
add("§4.6 candidates still >=4 in the 15 set", "7/7", "%d/7" % sum(1 for k in seven if c15c[k] >= 4), SRC_DE, "")
add("circRNAs at 5 datasets in the 15 set", "-", ";".join("%s(%s)" % (k, genes.get(k, "NA")) for k, v in c15c.items() if v >= 5), SRC_DE, "none are new vs the 19-set seven" if all(k in seven for k, v in c15c.items() if v >= 5) else "NEW 5-of-5 members present")
add("circRNAs at >=4 datasets (Fig 9c heatmap rows)", r19[4], r15[4], SRC_DE, "Figure 9c caption (circRNA count)")

with open(os.path.join(OUT, "thesis_numbers_15ds.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow(["item", "value_19_datasets", "value_15_datasets", "source", "location_in_v19"]); w.writerows(rows)

# ---- dataset table (all 19 kept)
with open(os.path.join(OUT, "thesis_dataset_table_19ds.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["accession", "organ_site", "cancer_type_label", "n_tumor", "n_normal", "n_tested", "n_sig", "n_up", "n_down", "TypeI", "TypeII", "TypeI_pct", "library_class", "single_tool_rows_pct", "in_cross_dataset_analysis", "exclusion_reason"])
    for d in DS19:
        c = typ(d); u, n = dirs(d); sp = float(summ[d]["single_tool_rows_pct"])
        inc = "yes" if d in DS15 else "no"
        why = "" if d in DS15 else "effective single-tool: %.1f%% of consensus rows supported by 1 tool (config tools=%s, min_tools=%s); exclusion threshold >=97.8%%" % (sp, summ[d]["config_tools"], summ[d]["config_min_tools"])
        w.writerow([d, META[d][0], META[d][1], summ[d]["n_tumor"], summ[d]["n_normal"], tested[d], nsig(d), u, n, c["Type_I"], c["Type_II"], pct(c["Type_I"], nsig(d)), META[d][2], summ[d]["single_tool_rows_pct"], inc, why])
print(len(rows), "rows")
