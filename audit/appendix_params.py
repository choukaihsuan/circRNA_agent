"""Build appendix tables A1 (parameter catalogue) and A2 (per-dataset effective values).

Source of every value, in order of precedence:
  1. "命令列"        the expanded shell command of the rule (snakemake -n -p, audit/dryrun_shell_dump.sh,
                     run with --cores so threads are not capped)
  2. "config 快照"   the key is present in config/projects/{GSE}.yaml
  3. "程式預設"      neither of the above; the code default the pipeline applies when the key is absent
  "未取得"           no source gives the effective value; never guessed

This is reconstructed from config, code defaults and the dry-run expansion. It is NOT produced from
the execution record: run_params.yaml (item B5) is not complete yet. Once B5 lands, these tables should
be regenerated from run_params.yaml, which is the reproducible source.

usage: appendix_params.py <dump_dir> <projects_dir> <committed_config_yaml> <out_dir> <recorded_threads_json> <ciriquant_yaml>
  committed_config_yaml: the repo config.yaml from git HEAD (the code defaults). The server copy
  is not synced and differs (threads 18, no tools section), so it must not be used.
"""
import collections
import glob
import os
import re
import sys

import yaml

import json

DUMP, PROJ, CFG_COMMITTED, OUT, RECORDED_JSON, CIRIQUANT_YAML = sys.argv[1:7]
COHORT15 = ["GSE108735", "GSE113230", "GSE121842", "GSE130078", "GSE133998", "GSE136569", "GSE143797",
            "GSE148036", "GSE171011", "GSE221107", "GSE248612", "GSE77509", "GSE97239", "PRJNA553289", "SRP156355"]
EXCLUDED4 = ["GSE58135", "GSE192410", "GSE192849", "GSE229705"]
# GSE323364 (cell line with drug treatment, not tumour vs normal) is outside the paper's 19 datasets.
GSES = sorted(set(COHORT15 + EXCLUDED4))
RECORDED = json.load(open(RECORDED_JSON))
CIRIQUANT_REF = yaml.safe_load(open(CIRIQUANT_YAML))["reference"]["hisat_index"]

PREFIX = ("input:", "output:", "log:", "jobid:", "reason:", "wildcards:", "resources:", "rule ", "#",
          "Config file", "Building DAG", "Job stats", "Reasons", "### ")


def blocks(path):
    """rule -> list of expanded command lines (header lines removed)."""
    out = collections.defaultdict(list)
    cur = None
    for raw in open(path, encoding="utf-8", errors="replace"):
        line = raw.rstrip("\n")
        if line.startswith("### RULE "):
            cur = line.split()[2]
            continue
        if line.startswith("### EXIT"):
            cur = None
            continue
        s = line.strip()
        # Snakemake's "Job stats" rows look like "dcc   1   1   1" (rule name and counts); not commands
        if re.fullmatch(r"[A-Za-z_0-9]+\s+\d+\s+\d+\s+\d+", s) or s == cur:
            continue
        if cur and s and not s.startswith(PREFIX):
            out[cur].append(s)
    return out


def flag(lines, tok, script=None):
    """Value after `tok` in the first matching command line; '(flag)' if no value; None if absent."""
    for s in lines:
        if script and f"scripts/{script}" not in s:
            continue
        parts = s.split()
        for i, t in enumerate(parts):
            if t == tok:
                nxt = parts[i + 1] if i + 1 < len(parts) else ""
                if nxt and not nxt.startswith(("-", "|", "&", "<", ">", "2>", "(", ")")) \
                        and not any(c in nxt for c in "<>|;"):
                    return nxt.strip("'\"")
                return "(flag)"
    return None


def has_flag(lines, tok, script=None):
    return flag(lines, tok, script) is not None


def snap_flat(d, pre=""):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(snap_flat(v, f"{pre}.{k}" if pre else str(k)))
    else:
        out[pre] = d
    return out


cfg_default = snap_flat(yaml.safe_load(open(CFG_COMMITTED)) or {})

# Defaults the pipeline applies when a key is absent and it is not in config.yaml (rule .get fallbacks).
RULE_DEFAULT = {
    "de.biomarker.cross_method_min": 3, "de.biomarker.stability_min_ratio": 0.5,
    "de.biomarker.stability_max_candidates": 500, "de.biomarker.sig_cap": 10, "de.biomarker.fc_cap": 5,
    "de.isoform_fdr_cutoff": 0.1, "de.delta_iui_cutoff": 0.1,
    "de.filter_expr_min_count": 5, "de.filter_expr_min_total_count": 15, "de.fsj_concordance_lfc": 0.0,
    "de.deseq2_size_factor_type": "poscounts", "de.paired_design_override": "auto",
    "de.heatmap_top_n": 10, "de.interaction_top_n": 50, "de.encori_clip_exp_num": 1,
    "de.encori_program_num": 2, "consensus.min_confidence": 0.0, "consensus.adaptive": True,
    "consensus.adaptive_ratio": 0.1, "consensus.qc_bsj_threshold": 5,
    "consensus.max_junction_ratio": 1.0, "tools.ciriquant.anchor": 5, "tools.ciriquant.no_gene": False,
    "tools.dcc.filter_recommended": True, "tools.dcc.filter_mito": True, "tools.dcc.end_tol": 5,
    "tools.dcc.max_length": 1000000, "tools.dcc.min_length": 30, "de.method": "edgeR_ciriquant",
    "study_title": "", "download.max_concurrent": 2, "download.retry": 3, "download.ascp_speed": "500m",
    "de.tumor_label": "tumor", "de.normal_label": "normal",
}
COMPONENTS = ["sig", "fc", "conf", "known", "mirna", "rbp"]


def dflt(name):
    v = RULE_DEFAULT.get(name, cfg_default.get(name))
    if name == "de.biomarker.components":
        return "all enabled, weight 1.0"
    return v if v is not None else "未取得"



def fmt(v):
    if isinstance(v, bool):
        return "True" if v else "False"
    if isinstance(v, float) and v.is_integer():
        return str(v)
    return str(v)


def facts_for(g):
    """Effective values for one dataset. Each entry: name -> (value, source)."""
    cfg = yaml.safe_load(open(os.path.join(PROJ, f"{g}.yaml"), encoding="utf-8")) or {}
    snap = snap_flat(cfg)
    b = blocks(os.path.join(DUMP, f"{g}.txt"))
    F = {}

    def from_snap(name, default_key=None):
        if name in snap:
            return snap[name], "config 快照"
        d = RULE_DEFAULT.get(name, cfg_default.get(name))
        if d is None:
            return "未取得", "未取得"
        return d, "程式預設"

    def cmd(name, lines, tok, default, script=None):
        v = flag(lines, tok, script) if lines else None
        if v is None:
            return default, "程式預設（命令列未輸出）"
        return v, "命令列"

    cons = b.get("consensus_filter", [])
    # a consensus command was expanded only if it carries the consensus flags; otherwise (for example
    # a MissingInput error block) nothing about it is known, and the values must not default
    cons_ok = any("scripts/consensus_filter.py" in x and "--min-tools" in x for x in cons)
    if not cons_ok:
        cons = []
    F["consensus.min_tools"] = cmd("", cons, "--min-tools", dflt("consensus.min_tools"), "consensus_filter.py")
    F["consensus.tools"] = (("ciriquant,dcc" if has_flag(cons, "--cirique", "consensus_filter.py")
                             and has_flag(cons, "--dcc", "consensus_filter.py") else
                             "ciriquant" if has_flag(cons, "--cirique", "consensus_filter.py") else
                             "dcc" if has_flag(cons, "--dcc", "consensus_filter.py") else "未取得"),
                            "命令列" if cons else "未取得")
    F["consensus.slop"] = cmd("", cons, "--slop", dflt("consensus.slop"), "consensus_filter.py")
    F["consensus.min_bsj_reads"] = cmd("", cons, "--min-bsj", dflt("consensus.min_bsj_reads"), "consensus_filter.py")
    F["consensus.max_junction_ratio"] = cmd("", cons, "--max-junction-ratio", dflt("consensus.max_junction_ratio"), "consensus_filter.py")
    F["consensus.qc_bsj_threshold"] = cmd("", cons, "--qc-bsj-threshold", dflt("consensus.qc_bsj_threshold"), "consensus_filter.py")
    if cons:
        F["consensus.adaptive"] = (True if has_flag(cons, "--adaptive", "consensus_filter.py") else False,
                                   "命令列")
    else:
        F["consensus.adaptive"] = ("未取得", "未取得")
    F["consensus.adaptive_ratio"] = cmd("", cons, "--adaptive-ratio", dflt("consensus.adaptive_ratio"), "consensus_filter.py")
    F["consensus.min_confidence"] = cmd("", cons, "--min-confidence", 0.0, "consensus_filter.py")
    if not cons_ok:
        F["consensus.tools"] = ("未取得", "未取得（dry-run 無法展開 consensus）")
    if not cons_ok:
        for k in ["consensus.min_tools", "consensus.slop", "consensus.min_bsj_reads",
                  "consensus.max_junction_ratio", "consensus.qc_bsj_threshold", "consensus.adaptive",
                  "consensus.adaptive_ratio", "consensus.min_confidence"]:
            F[k] = ("未取得", "未取得（dry-run 無法展開 consensus）")
    uses_dcc = "dcc" in str(F["consensus.tools"][0]).split(",")
    dcc = b.get("dcc", [])
    if not uses_dcc:
        dcc = []
    if dcc:
        nr = None
        for s in dcc:
            m = re.search(r"-Nr (\d+) (\d+)", s)
            if m:
                nr = (m.group(1), m.group(2))
                break
        F["tools.dcc.count_threshold"] = (nr[0] if nr else "未取得", "命令列" if nr else "未取得")
        F["tools.dcc.replicate_threshold"] = (nr[1] if nr else "未取得", "命令列" if nr else "未取得")
        F["tools.dcc.filter_recommended"] = (has_flag(dcc, "-F") and True or False, "命令列")
        F["tools.dcc.filter_mito"] = (has_flag(dcc, "-M") and True or False, "命令列")
        F["tools.dcc.end_tol"] = cmd("", dcc, "-E", 5)
        F["tools.dcc.max_length"] = cmd("", dcc, "-m", 1000000)
        F["tools.dcc.min_length"] = cmd("", dcc, "-n", 30)
        F["tools.dcc.threads"] = cmd("", dcc, "-T", dflt("tools.dcc.threads"))
    else:
        if F["consensus.tools"][0] == "未取得":
            na = ("未取得", "未取得（dry-run 無法展開 consensus）")
        else:
            na = ("不適用", "不適用（consensus.tools 不含 dcc）") if not uses_dcc else ("未取得", "未取得（dry-run 無法展開）")
        for k in ["tools.dcc.count_threshold", "tools.dcc.replicate_threshold", "tools.dcc.filter_recommended",
                  "tools.dcc.filter_mito", "tools.dcc.end_tol", "tools.dcc.max_length", "tools.dcc.min_length",
                  "tools.dcc.threads"]:
            F[k] = na
    cq = b.get("ciriquant", [])
    rec = RECORDED.get(g, {})
    if rec:
        # the thread count the published CIRIquant GTFs were produced with (Snakemake metadata)
        best = max(rec.items(), key=lambda kv: kv[1])[0]
        F["threads"] = (best if len(rec) == 1 else "混合：" + "/".join(sorted(rec)),
                        "Snakemake metadata（已發表 GTF 的實際 -t）")
    else:
        F["threads"] = ("未取得", "未取得（無 CIRIquant metadata）")
    F["tools.ciriquant.anchor"] = cmd("", cq, "-a", 5)
    F["tools.ciriquant.no_gene"] = ((True if has_flag(cq, "--no-gene") else False) if cq else "未取得",
                                    "命令列" if cq else "未取得")
    star = b.get("star_align", [])
    if star and not cq:
        F["threads"] = cmd("", star, "--runThreadN", dflt("threads"))
    rk = b.get("rank_biomarkers_v2", [])
    if rk:
        F["de.biomarker.cross_method_min"] = cmd("", rk, "--cross-method-min", dflt("de.biomarker.cross_method_min"), "rank_biomarkers.py")
        F["de.biomarker.stability_min_ratio"] = cmd("", rk, "--stability-min-ratio", dflt("de.biomarker.stability_min_ratio"), "rank_biomarkers.py")
        F["de.biomarker.stability_max_candidates"] = cmd("", rk, "--stability-max-candidates", 500,
                                                         "rank_biomarkers.py")
        F["de.biomarker.sig_cap"] = cmd("", rk, "--sig-cap", 10, "rank_biomarkers.py")
        F["de.biomarker.fc_cap"] = cmd("", rk, "--fc-cap", 5, "rank_biomarkers.py")
        comp = flag(rk, "--components-json", "rank_biomarkers.py")
        F["de.biomarker.components"] = ((comp if comp and comp != "(flag)" else
                                         "all enabled, weight 1.0"), "命令列" if comp else "程式預設（命令列未輸出）")
        F["de.de_sig_by"] = cmd("", rk, "--de-sig-by", dflt("de.de_sig_by"), "rank_biomarkers.py")
        F["de.fdr_cutoff"] = cmd("", rk, "--fdr", dflt("de.fdr_cutoff"), "rank_biomarkers.py")
        F["de.log2fc_cutoff"] = cmd("", rk, "--lfc", dflt("de.log2fc_cutoff"), "rank_biomarkers.py")
        F["de.biomarker.primary_method"] = ("edgeR / DESeq2 / limma (one candidate file each)",
                                            "rule wildcard (de.smk rank_biomarkers_v2)")
    pi = b.get("predict_interactions", [])
    F["de.interaction_top_n"] = cmd("", pi, "--top-n", 50, "predict_interactions.py")
    F["de.encori_clip_exp_num"] = cmd("", pi, "--clip-exp-num", 1, "predict_interactions.py")
    F["de.encori_program_num"] = cmd("", pi, "--program-num", 2, "predict_interactions.py")

    # Items the dry-run cannot show (R scripts and generate_report use script: directives).
    for name in ["de.method", "de.tumor_label", "de.normal_label", "de.isoform_fdr_cutoff",
                 "de.delta_iui_cutoff", "de.filter_expr_min_count", "de.filter_expr_min_total_count",
                 "de.fsj_concordance_lfc", "de.deseq2_size_factor_type", "de.paired_design_override",
                 "de.heatmap_top_n", "download.max_concurrent", "download.retry", "download.ascp_speed",
                 "study_title", "genome.hisat2_index"]:
        if name in F:
            continue
        if name == "genome.hisat2_index":
            F[name] = (snap.get("genome.hisat2_index", "未取得"), "config 快照（有效值未取得）")
            continue
        F[name] = from_snap(name)
    F.setdefault("de.biomarker.components", ("all enabled, weight 1.0", "程式預設（命令列未輸出）"))
    return F, snap


FACTS = {}
SNAPS = {}
for g in GSES:
    FACTS[g], SNAPS[g] = facts_for(g)

# A1 registry: name, tier, description
REG = [
    ("threads", 1, "CPU threads for CIRIquant and STAR (fastp uses the same value)"),
    ("tools.dcc.threads", 1, "DCC -T thread count"),
    ("download.max_concurrent", 1, "parallel SRA download jobs (operational only)"),
    ("download.retry", 1, "retry count for SRA download (operational only)"),
    ("download.ascp_speed", 1, "Aspera speed (not used on the current server)"),
    ("de.heatmap_top_n", 1, "heatmap top-N in the PDF output (display only)"),
    ("study_title", 1, "study title shown in the report (display only; text, varies by dataset)"),
    ("consensus.min_tools", 2, "minimum tools that must agree in consensus voting"),
    ("consensus.tools", 2, "tools used for consensus (ciriquant, dcc)"),
    ("consensus.slop", 2, "coordinate tolerance in bp when matching loci across tools"),
    ("consensus.min_bsj_reads", 2, "minimum BSJ reads per tool before voting"),
    ("consensus.max_junction_ratio", 2, "CIRIquant pseudo-circRNA QC: maximum BSJ/FSJ ratio"),
    ("consensus.qc_bsj_threshold", 2, "BSJ count below which the BSJ/FSJ ratio QC applies"),
    ("consensus.adaptive", 2, "fall back to one tool when detection counts are very unbalanced"),
    ("consensus.adaptive_ratio", 2, "min/max detection ratio that triggers the adaptive fallback"),
    ("consensus.min_confidence", 2, "minimum confidence score after voting"),
    ("tools.dcc.count_threshold", 2, "DCC -Nr first argument: minimum junction count (the published runs used 2; the consensus fallback count uses the same value)"),
    ("tools.dcc.replicate_threshold", 2, "DCC -Nr second argument: minimum replicate count"),
    ("tools.dcc.filter_recommended", 2, "DCC -F: recommended filtering"),
    ("tools.dcc.filter_mito", 2, "DCC -M: exclude mitochondrial circRNAs"),
    ("tools.dcc.end_tol", 2, "DCC -E: end tolerance (bp)"),
    ("tools.dcc.max_length", 2, "DCC -m: maximum circRNA length"),
    ("tools.dcc.min_length", 2, "DCC -n: minimum circRNA length"),
    ("tools.ciriquant.anchor", 2, "CIRIquant -a: anchor length for back-splice junction reads"),
    ("tools.ciriquant.no_gene", 2, "CIRIquant --no-gene: skip gene abundance"),
    ("de.method", 2, "primary DE method reported in the report"),
    ("de.de_sig_by", 2, "significance column: pvalue (nominal) or padj/qvalue/auto"),
    ("de.fdr_cutoff", 2, "significance threshold (nominal P in this study)"),
    ("de.log2fc_cutoff", 2, "minimum |log2FC| for a significant circRNA"),
    ("de.isoform_fdr_cutoff", 2, "within-gene FDR threshold for isoform switching"),
    ("de.delta_iui_cutoff", 2, "minimum |delta IUI| for isoform switching"),
    ("de.filter_expr_min_count", 2, "filterByExpr min.count"),
    ("de.filter_expr_min_total_count", 2, "filterByExpr min.total.count"),
    ("de.fsj_concordance_lfc", 2, "minimum |logFC_fsj| for Type II concordance"),
    ("de.deseq2_size_factor_type", 2, "DESeq2 estimateSizeFactors type"),
    ("de.paired_design_override", 2, "paired design: auto (detect patient_id), force_paired, force_unpaired"),
    ("de.tumor_label", 2, "label used for the tumour group"),
    ("de.normal_label", 2, "label used for the normal group (sign of log2FC)"),
    ("de.interaction_top_n", 2, "per method, the top-N significant circRNAs up-regulated and top-N down-regulated are queried for miRNA/RBP interactions"),
    ("de.encori_clip_exp_num", 2, "ENCORI minimum CLIP experiments"),
    ("de.encori_program_num", 2, "ENCORI minimum prediction programmes"),
    ("de.biomarker.components", 2, "six biomarker components: enabled flag and weight each"),
    ("de.biomarker.sig_cap", 2, "cap on -log10(p) in the significance component"),
    ("de.biomarker.fc_cap", 2, "cap on |log2FC| in the fold-change component"),
    ("de.biomarker.cross_method_min", 2, "minimum methods in which a candidate is significant"),
    ("de.biomarker.stability_min_ratio", 2, "LOO stability: minimum |log2FC| ratio to the full value"),
    ("de.biomarker.stability_max_candidates", 2, "maximum candidates tested in the LOO stability step"),
    ("de.biomarker.primary_method", 2, "method whose table is scored (one candidate file per method)"),
]

LIMIT = {"consensus.tools": "ciriquant,dcc", "de.biomarker.primary_method": "per method"}


def effective(name, g):
    v = FACTS[g].get(name)
    if v is None:
        return ("未取得", "未取得")
    return v


def canon(v):
    if isinstance(v, tuple):
        v = v[0]
    if isinstance(v, str) and v == "":
        return "未設定"
    return fmt(v) if not isinstance(v, str) else v


# A2 columns: varying parameters across the 20 projects, plus forced threads, plus hisat2 index.
ALL_NAMES = [n for n, _, _ in REG]
varying = []
for n in ALL_NAMES:
    vals = {canon(effective(n, g)) for g in GSES}
    if len(vals) > 1 or n == "threads":
        varying.append(n)


os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "appendix_A1_parameters.tsv"), "w", encoding="utf-8") as fh:
    fh.write("# 收錄範圍：Tier 1（資源與顯示參數）與 Tier 2（會改變分析結果的參數）。排除 Tier 3：路徑、專案身分、參考基因組與外部資源位置、個資。\n")
    fh.write("# 值的來源優先序：CIRIquant metadata（threads）→ dry-run 展開命令列 → 專案快照 → 程式預設。取不到者標「未取得」，不猜。\n")
    fh.write("# 本表由 config、程式預設與 dry-run 展開重建，不是由執行紀錄產生；B5（run_params.yaml）完成後應改由其產生。\n")
    fh.write("# CIRIquant 的索引與參數取自單一共用檔 config/ciriquant.yaml（手動維護，無逐次執行紀錄），各資料集執行當下的值無法逐一還原。\n")
    fh.write("參數名\t層級\t本研究使用值\t程式預設值\t本研究值＝預設？\t來源\t作用\n")
    fh.write(f"reference.hisat_index（ciriquant.yaml，共用）\t3（例外，依要求列出）\t{CIRIQUANT_REF}\t—\t—\tciriquant.yaml\tCIRIquant 實際使用的 HISAT2 索引（唯一位置）\n")
    for name, tier, desc in REG:
        vals = {canon(effective(name, g)) for g in COHORT15}
        srcs = sorted({effective(name, g)[1] for g in COHORT15})
        if len(vals) == 1:
            use = next(iter(vals))
        else:
            use = "見表 A2"
        default = dflt(name) if name in LIMIT or name == "de.biomarker.components" else dflt(name)
        if name == "de.biomarker.primary_method":
            default = "per method (wildcard)"
        same = "same" if len(vals) == 1 and use == fmt(default) else "DIFF"
        if use == "見表 A2":
            same = "見表 A2"
        if name == "de.biomarker.primary_method":
            same = "n/a（wildcard）"
        fh.write(f"{name}\t{tier}\t{use}\t{fmt(default)}\t{same}\t{'; '.join(srcs)}\t{desc}\n")

with open(os.path.join(OUT, "appendix_A2_per_dataset.tsv"), "w", encoding="utf-8") as fh:
    fh.write("# 每組資料集的有效參數（涵蓋論文 15 組與排除 4 組；只列實際有差異的欄位，threads 強制列出）。\n")
    fh.write("# threads 取自已發表 GTF 的 CIRIquant metadata。GSE113230 快照寫 12，但其 GTF 由 threads 8 產生（沿用舊檔）。\n")
    fh.write("dataset\t納入分析\t" + "\t".join(varying) + "\n")
    for g in GSES:
        inc = ("論文 15 組" if g in COHORT15 else "排除（不供應論文數字）" if g in EXCLUDED4
               else "未納入論文 19 組（待確認）")
        row = [g, inc] + [canon(effective(n, g)) for n in varying]
        fh.write("\t".join(row) + "\n")
print("A1 rows", len(REG), "A2 cols", len(varying), "datasets", len(GSES))
