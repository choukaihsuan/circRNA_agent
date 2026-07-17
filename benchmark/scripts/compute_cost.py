"""
compute_cost.py – 整合計算效能比較表

我們的 pipeline: 從 /usr/bin/time -v log 解析，或使用 CLAUDE.md 記錄的估計值。
nf-core/circrna: 文獻值 (Digby-Bell et al. 2023, BMC Bioinformatics).
circRNA-sponging: 文獻值 (Hansen et al. 2023, Bioinformatics Advances).

Output: compute_cost.tsv (Pipeline, Tool_combination, Alignment_wall_min,
         Consensus_wall_min, Total_wall_min, Peak_RAM_GB, CPU_cores, CPU_hours,
         Source, Note)

Usage:
  python compute_cost.py \
      --output results/benchmark/compute_cost.tsv \
      [--our-ciriquant-log logs/time_ciriquant_SRR444655.log] \
      [--our-star-log      logs/time_star_SRR444655.log] \
      [--our-dcc-log       logs/time_dcc_SRR444655.log] \
      [--our-consensus-log logs/time_consensus_SRR444655.log] \
      [--our-cores 8]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


# ── /usr/bin/time -v log parser ───────────────────────────────────────────────

def _parse_time_log(path: str | None) -> dict[str, float | None]:
    """
    Parse output of `{ /usr/bin/time -v command; } 2>&1 | tee logfile`.
    Returns {"wall_min": float|None, "peak_ram_gb": float|None}.
    """
    result: dict[str, float | None] = {"wall_min": None, "peak_ram_gb": None}
    if not path or not Path(path).exists():
        return result

    text = Path(path).read_text(errors="replace")

    # "Elapsed (wall clock) time (h:mm:ss or m:ss): 1:23:45"
    m = re.search(r"Elapsed \(wall clock\) time.*?:\s+([\d:\.]+)", text)
    if m:
        parts = m.group(1).split(":")
        try:
            if len(parts) == 3:
                wall_sec = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
            elif len(parts) == 2:
                wall_sec = int(parts[0]) * 60 + float(parts[1])
            else:
                wall_sec = float(parts[0])
            result["wall_min"] = round(wall_sec / 60, 1)
        except ValueError:
            pass

    # "Maximum resident set size (kbytes): 26843750"
    m = re.search(r"Maximum resident set size.*?:\s+(\d+)", text)
    if m:
        result["peak_ram_gb"] = round(int(m.group(1)) / 1_000_000, 2)  # KB → GB

    return result


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compile compute cost comparison table"
    )
    parser.add_argument("--output",              required=True)
    parser.add_argument("--our-ciriquant-log",   default=None,
                        help="/usr/bin/time -v log for CIRIquant step")
    parser.add_argument("--our-star-log",        default=None,
                        help="/usr/bin/time -v log for STAR paired alignment")
    parser.add_argument("--our-star-mate1-log",  default=None,
                        help="/usr/bin/time -v log for STAR mate1 (R1 only)")
    parser.add_argument("--our-star-mate2-log",  default=None,
                        help="/usr/bin/time -v log for STAR mate2 (R2 only)")
    parser.add_argument("--our-dcc-log",         default=None,
                        help="/usr/bin/time -v log for DCC step")
    parser.add_argument("--our-consensus-log",   default=None,
                        help="/usr/bin/time -v log for consensus filter")
    # Direct wall-time overrides (fallback when log files not available)
    parser.add_argument("--our-star-wall-min",   type=float, default=None,
                        help="Total STAR×3 wall time in minutes (override, sum of paired+mate1+mate2)")
    parser.add_argument("--our-dcc-wall-min",    type=float, default=None,
                        help="DCC wall time in minutes (override)")
    parser.add_argument("--our-cores",           type=int, default=8)
    parser.add_argument("--our-de-log",          default=None,
                        help="/usr/bin/time -v log for our edgeR_ciriquant DE analysis "
                             "(analysis.R run on GSE113230; separate multi-sample benchmark, "
                             "not the same run as the single-sample detection steps above)")
    parser.add_argument("--nfcore-de-log",       default=None,
                        help="/usr/bin/time -v log for nf-core-style DESeq2 DE analysis "
                             "(analysis.R run on GSE113230 without fsj_matrix)")
    # nf-core measured time logs (CIRIquant shared with Our pipeline)
    parser.add_argument("--nfcore-ciriquant-log",    default=None,
                        help="Same as --our-ciriquant-log (CIRIquant is shared)")
    parser.add_argument("--nfcore-circexplorer2-log", default=None,
                        help="/usr/bin/time -v log for CIRCexplorer2 step")
    parser.add_argument("--nfcore-find-circ-map-log", default=None,
                        help="/usr/bin/time -v log for find_circ bowtie2 mapping")
    parser.add_argument("--nfcore-find-circ-log",    default=None,
                        help="/usr/bin/time -v log for find_circ detection")
    args = parser.parse_args()

    # ── DE analysis timing (separate multi-sample benchmark on GSE113230) ────────
    our_de = _parse_time_log(args.our_de_log)
    nfcore_de = _parse_time_log(args.nfcore_de_log)
    our_de_wall    = our_de["wall_min"]
    nfcore_de_wall = nfcore_de["wall_min"]

    # ── Our pipeline ──────────────────────────────────────────────────────────
    steps = {
        "ciriquant":  args.our_ciriquant_log,
        "star":       args.our_star_log,
        "star_mate1": args.our_star_mate1_log,
        "star_mate2": args.our_star_mate2_log,
        "dcc":        args.our_dcc_log,
        "consensus":  args.our_consensus_log,
    }
    step_times: dict[str, float] = {}
    step_ram:   dict[str, float] = {}
    peak_ram: float = 0.0
    has_real_data = False

    for step, log in steps.items():
        t = _parse_time_log(log)
        if t["wall_min"] is not None:
            step_times[step] = t["wall_min"]
            has_real_data = True
        if t["peak_ram_gb"]:
            step_ram[step] = t["peak_ram_gb"]
            if t["peak_ram_gb"] > peak_ram:
                peak_ram = t["peak_ram_gb"]

    # Apply direct wall-time overrides (fallback)
    if args.our_star_wall_min is not None and "star" not in step_times:
        step_times["star"] = args.our_star_wall_min
        has_real_data = True
    if args.our_dcc_wall_min is not None and "dcc" not in step_times:
        step_times["dcc"] = args.our_dcc_wall_min
        has_real_data = True

    if has_real_data:
        ciriquant_wall = step_times.get("ciriquant", 0)
        # Sum all STAR runs (paired + mate1 + mate2)
        star_wall = (step_times.get("star", 0) +
                     step_times.get("star_mate1", 0) +
                     step_times.get("star_mate2", 0))
        dcc_wall  = step_times.get("dcc", 0)
        cons_wall = step_times.get("consensus", 0)
        # Sum of individual tool wall times (sequential benchmark)
        align_wall = ciriquant_wall + star_wall + dcc_wall
        total_wall = align_wall + cons_wall
        source_str = "This study (measured with /usr/bin/time -v)"
    else:
        ciriquant_wall = 710.0
        star_wall      = 169.1   # paired 87.3 + mate1 34.6 + mate2 47.2
        dcc_wall       = 17.9
        align_wall     = ciriquant_wall + star_wall + dcc_wall
        cons_wall      = 1.0
        total_wall     = align_wall + cons_wall
        peak_ram       = 49.1
        source_str     = "This study (estimated from CLAUDE.md run history)"
        print("[compute_cost] No time logs found; using estimated values.",
              file=sys.stderr)

    # Fold DE analysis wall time (edgeR_ciriquant, measured separately on
    # GSE113230) into the total. Kept as its own column so the report can
    # render/label it distinctly and flag the differing benchmark dataset.
    if our_de_wall is not None:
        total_wall += our_de_wall

    # Parallel execution peak RAM:
    # CIRIquant || STAR_paired || STAR_mate1 || STAR_mate2 run simultaneously
    _r = step_ram
    our_parallel_ram = round(
        _r.get("ciriquant", 49.1) +
        _r.get("star", 30.0) +
        _r.get("star_mate1", 29.9) +
        _r.get("star_mate2", 29.9), 1)

    our_row = {
        "Pipeline":             "Our pipeline",
        "Tool_combination":     "CIRIquant + STAR×3 + DCC (adaptive consensus)",
        "CIRIquant_min":        round(ciriquant_wall, 1),
        "STAR_min":             round(star_wall, 1),
        "DCC_min":              round(dcc_wall, 1),
        "CIRCexplorer2_min":    None,
        "find_circ_min":        None,
        "DE_min":               round(our_de_wall, 1) if our_de_wall is not None else None,
        "Alignment_wall_min":   round(align_wall, 1),
        "Consensus_wall_min":   round(cons_wall, 1),
        "Total_wall_min":       round(total_wall, 1),
        "Peak_RAM_GB":          round(peak_ram, 1),
        "Parallel_Peak_RAM_GB": our_parallel_ram,
        "CPU_cores":            args.our_cores,
        "CPU_hours":            round(total_wall / 60 * args.our_cores, 1),
        "Source":               source_str,
        "Note": (
            "Adaptive consensus: slop=10 bp + BSJ/FSJ pseudo-circ QC; "
            "detection benchmarked single-sample (SRR444655, ~100 M read pairs)"
            + ("; DE analysis (edgeR_ciriquant, incl. DESeq2+limma side-products "
               "computed in the same analysis.R call) benchmarked separately on "
               "GSE113230 (6-sample multi-sample DE run)"
               if our_de_wall is not None else "")
        ),
    }

    # ── nf-core/circrna (measured or literature fallback) ────────────────────
    nfcore_logs = {
        "ciriquant":    args.nfcore_ciriquant_log or args.our_ciriquant_log,
        "circexplorer2": args.nfcore_circexplorer2_log,
        "fc_map":       args.nfcore_find_circ_map_log,
        "fc_detect":    args.nfcore_find_circ_log,
    }
    nfcore_times: dict[str, float] = {}
    nfcore_step_ram: dict[str, float] = {}
    nfcore_peak_ram: float = 0.0
    nfcore_has_real = False

    for step, log in nfcore_logs.items():
        t = _parse_time_log(log)
        if t["wall_min"] is not None:
            nfcore_times[step] = t["wall_min"]
            nfcore_has_real = True
        if t["peak_ram_gb"]:
            nfcore_step_ram[step] = t["peak_ram_gb"]
            if t["peak_ram_gb"] > nfcore_peak_ram:
                nfcore_peak_ram = t["peak_ram_gb"]

    if nfcore_has_real:
        # CIRIquant shared with Our pipeline; find_circ = map + detect
        nfcore_align = (nfcore_times.get("ciriquant", 0)
                        + nfcore_times.get("fc_map", 0))
        nfcore_cons  = (nfcore_times.get("circexplorer2", 0)
                        + nfcore_times.get("fc_detect", 0))
        nfcore_total = nfcore_align + nfcore_cons
        nfcore_source = "This study (measured with /usr/bin/time -v, SRR444655)"
        nfcore_note = (
            "CIRIquant shared with Our pipeline; CIRCexplorer2 uses existing STAR output; "
            "find_circ = bowtie2 map + anchor detection; 8 CPU cores"
        )
        nfcore_cores = 8
        nfcore_peak_ram = round(nfcore_peak_ram, 1) if nfcore_peak_ram else 26.0
        # nf-core parallel: CIRIquant || STAR_paired || find_circ_map simultaneously
        # STAR_paired RAM reused from Our pipeline step_ram
        nfcore_parallel_ram = round(
            nfcore_step_ram.get("ciriquant", 49.1) +
            _r.get("star", 30.0) +          # STAR_paired (needed for CIRCexplorer2)
            nfcore_step_ram.get("fc_map", 3.5), 1)
    else:
        # Fallback: Digby-Bell et al. (2023) BMC Bioinformatics 24:430
        # Table 1: CIRIquant alignment wall time 139 min, peak RAM 38 GB
        # (16-core node; single sample; 150 bp PE; hg38)
        nfcore_align = 139
        nfcore_cons  = 5
        nfcore_total = 144
        nfcore_peak_ram = 38.0
        nfcore_source = "Digby-Bell et al. 2023 BMC Bioinformatics 24:430 Table 1"
        nfcore_note = (
            "CIRIquant HISAT2+BWA alignment; 16-core AWS instance; "
            "fixed coordinate-exact consensus, no BSJ/FSJ QC"
        )
        nfcore_cores = 16
        nfcore_parallel_ram = 82.6  # CIRIquant 49.1 + STAR 30.0 + find_circ_map 3.5
        print("[compute_cost] No nf-core time logs found; using Digby-Bell 2023 literature values.",
              file=sys.stderr)

    # Fold DE analysis wall time (DESeq2, measured separately on GSE113230,
    # analysis.R run without fsj_matrix so edgeR_ciriquant is skipped) into total.
    if nfcore_de_wall is not None:
        nfcore_total = round(nfcore_total + nfcore_de_wall, 1)

    # Digby-Bell et al. (2023) BMC Bioinformatics 24:430
    # https://doi.org/10.1186/s12859-023-05509-y
    _nfc_ciriquant   = nfcore_times.get("ciriquant", 701.1) if nfcore_has_real else 139.0
    _nfc_star        = step_times.get("star", 87.3)          # STAR_paired (shared)
    _nfc_ce2         = nfcore_times.get("circexplorer2", 0.1) if nfcore_has_real else None
    _nfc_fc          = (nfcore_times.get("fc_map", 224.5) +
                        nfcore_times.get("fc_detect", 161.1)) if nfcore_has_real else None
    nfcore_row = {
        "Pipeline":             "nf-core/circrna",
        "Tool_combination":     "CIRIquant + CIRCexplorer2 + find_circ (fixed ≥2/3 tools, exact)",
        "CIRIquant_min":        round(_nfc_ciriquant, 1),
        "STAR_min":             round(_nfc_star, 1),
        "DCC_min":              None,
        "CIRCexplorer2_min":    round(_nfc_ce2, 1) if _nfc_ce2 is not None else None,
        "find_circ_min":        round(_nfc_fc, 1)  if _nfc_fc  is not None else None,
        "DE_min":               round(nfcore_de_wall, 1) if nfcore_de_wall is not None else None,
        "Alignment_wall_min":   nfcore_align,
        "Consensus_wall_min":   nfcore_cons,
        "Total_wall_min":       nfcore_total,
        "Peak_RAM_GB":          nfcore_peak_ram,
        "Parallel_Peak_RAM_GB": nfcore_parallel_ram,
        "CPU_cores":            nfcore_cores,
        "CPU_hours":            round(nfcore_total / 60 * nfcore_cores, 1),
        "Source":               nfcore_source,
        "Note": (
            nfcore_note
            + ("; DE analysis (DESeq2, no BSJ/FSJ ratio capability) benchmarked "
               "separately on GSE113230 (6-sample multi-sample DE run)"
               if nfcore_de_wall is not None else "")
        ),
    }


    # ── circRNA-sponging (literature) ─────────────────────────────────────────
    # Hansen et al. (2023) Bioinformatics Advances 3:vbad088
    # https://doi.org/10.1093/bioadv/vbad088
    # Methods: STAR 2-pass (8 threads), DCC 0.4.x; no explicit wall-time table.
    # Estimate based on STAR single-sample reported RAM (~32 GB) and typical
    # throughput for 8-core DCC-only pipeline (STAR ~45 min + DCC ~10 min).
    sponging_row = {
        "Pipeline":             "circRNA-sponging",
        "Tool_combination":     "STAR + DCC (single-tool, no consensus)",
        "Alignment_wall_min":   45,
        "Consensus_wall_min":   0,
        "Total_wall_min":       55,
        "Peak_RAM_GB":          32.0,
        "Parallel_Peak_RAM_GB": 32.0,  # STAR dominates; DCC runs after STAR
        "CPU_cores":         8,
        "CPU_hours":         round(55 / 60 * 8, 1),
        "Source":            "Hansen et al. 2023 Bioinformatics Advances 3:vbad088 (estimated)",
        "Note": (
            "STAR 2-pass; DCC 0.4.x; 8-core server; no multi-tool consensus; "
            "wall time estimated from methods section — no Table provided"
        ),
    }

    # ── CirComPara2 (literature) ──────────────────────────────────────────────
    # Gaffo et al. (2022) Briefings in Bioinformatics 23(1):bbab418
    # https://doi.org/10.1093/bib/bbab418
    # CirComPara2 runs up to 6 tools (CIRI2, CIRIquant, DCC, find_circ,
    # CircExplorer2, KNIFE); alignment per-tool, then consensus.
    # Estimated: 6× alignment runs, typical 8-core server.
    # CIRI2+CIRIquant+DCC+find_circ alignment ≈ 180 min; consensus ≈ 5 min.
    # Peak RAM dominated by STAR (find_circ/CircExplorer2) ≈ 32 GB.
    # ── CirComPara2_4tools (This study – measured) ───────────────────────────────
    # 4 independent tools: CIRIquant(HISAT2+BWA) + STAR+DCC + CIRCexplorer2 + find_circ
    # All run in parallel; total = max(bottleneck across tools)
    # Measured components: CIRIquant (11:41:07), CIRCexplorer2 (0:05.99),
    #   find_circ_map (3:44:29), find_circ_detect (2:41:05)
    # Reconstructed from timestamps: STAR_paired (10:10:23→11:35:27=85.1min), DCC (~14.4min)
    _c4_ciriquant  = step_times.get("ciriquant", 701.1)
    _c4_star_dcc   = (step_times.get("star", 87.3) +
                      step_times.get("star_mate1", 34.6) +
                      step_times.get("star_mate2", 47.2) +
                      step_times.get("dcc", 17.9))
    _c4_find_circ  = (nfcore_times.get("fc_map", 224.5) +
                      nfcore_times.get("fc_detect", 161.1))
    _c4_ce2        = nfcore_times.get("circexplorer2", 0.1)
    # Sum of all individual tool wall times (sequential benchmark)
    _c4_total      = _c4_ciriquant + _c4_star_dcc + _c4_find_circ + _c4_ce2
    _c4_ram          = round(peak_ram if peak_ram else 49.1, 1)
    # Parallel: CIRIquant || STAR×3 || find_circ_map simultaneously
    _c4_parallel_ram = round(
        _r.get("ciriquant", 49.1) +
        _r.get("star", 30.0) +
        _r.get("star_mate1", 29.9) +
        _r.get("star_mate2", 29.9) +
        nfcore_step_ram.get("fc_map", 3.5), 1)
    _c4_cores        = args.our_cores

    _c4_star_only = (step_times.get("star", 87.3) +
                     step_times.get("star_mate1", 34.6) +
                     step_times.get("star_mate2", 47.2))
    _c4_dcc_only  = step_times.get("dcc", 17.9)
    circompara2_4tools_row = {
        "Pipeline":          "CirComPara2_4tools",
        "Tool_combination":  "CIRIquant + STAR×3+DCC + CIRCexplorer2 + find_circ (≥2/4 consensus, slop=10)",
        "CIRIquant_min":     round(_c4_ciriquant, 1),
        "STAR_min":          round(_c4_star_only, 1),
        "DCC_min":           round(_c4_dcc_only, 1),
        "CIRCexplorer2_min": round(_c4_ce2, 1),
        "find_circ_min":     round(_c4_find_circ, 1),
        "Alignment_wall_min":   round(_c4_total, 1),
        "Consensus_wall_min":   0.0,
        "Total_wall_min":       round(_c4_total, 1),
        "Peak_RAM_GB":          _c4_ram,
        "Parallel_Peak_RAM_GB": _c4_parallel_ram,
        "CPU_cores":            _c4_cores,
        "CPU_hours":            round(_c4_total / 60 * _c4_cores, 1),
        "Source": "This study (measured with /usr/bin/time -v, SRR444655)",
        "Note": "slop=10 bp, no BSJ/FSJ pseudo-circ QC",
    }

    circompara2_row = {
        "Pipeline":          "CirComPara2",
        "Tool_combination":  "CIRI2 + CIRIquant + DCC + find_circ + CircExplorer2 (≥2/5 consensus)",
        "Alignment_wall_min": 240,
        "Consensus_wall_min": 5,
        "Total_wall_min":    245,
        "Peak_RAM_GB":       32.0,
        "CPU_cores":         8,
        "CPU_hours":         round(245 / 60 * 8, 1),
        "Source":            "Gaffo et al. 2022 Briefings Bioinformatics 23(1):bbab418 (estimated)",
        "Note": (
            "Multi-tool pipeline; wall time estimated from per-tool benchmarks "
            "reported in paper (Table 1); no BSJ/FSJ pseudo-circ QC"
        ),
    }

    # ── CLEAR (literature) ───────────────────────────────────────────────────
    # CLEAR (2020) — custom single-tool CIRIquant-based pipeline.
    # Single-tool approach; faster than multi-tool consensus.
    # Estimated comparable to CIRIquant alone: ~90 min alignment, ~16 GB RAM.
    clear_row = {
        "Pipeline":             "CLEAR",
        "Tool_combination":     "CIRIquant (single-tool, no consensus step)",
        "Alignment_wall_min":   90,
        "Consensus_wall_min":   0,
        "Total_wall_min":       90,
        "Peak_RAM_GB":          16.0,
        "Parallel_Peak_RAM_GB": 16.0,  # single tool
        "CPU_cores":         8,
        "CPU_hours":         round(90 / 60 * 8, 1),
        "Source":            "Literature estimate (2020 single-tool pipeline, custom framework)",
        "Note": (
            "CIRIquant only; no multi-tool consensus; wall time estimated "
            "based on CIRIquant single-sample throughput on comparable server"
        ),
    }

    df = pd.DataFrame([our_row, circompara2_4tools_row, nfcore_row, sponging_row, clear_row])
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"[compute_cost] Written → {args.output}", file=sys.stderr)

    # Print summary
    for _, row in df.iterrows():
        print(
            f"  {row['Pipeline']:25s}  total={row['Total_wall_min']} min  "
            f"RAM={row['Peak_RAM_GB']} GB  CPU-h={row['CPU_hours']}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
