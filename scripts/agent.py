"""
agent.py – circRNA Agent CLI entry point.

Usage examples:
  python scripts/agent.py --gse GSE113230
  python scripts/agent.py --runinfo RunInfo.csv --interactive
  python scripts/agent.py --gse GSE113230 --cores 16 --dry-run
  python scripts/agent.py --setup-ciriquant
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Make scripts/ importable when called from project root
sys.path.insert(0, str(Path(__file__).parent))

from download_geo import prepare_metadata
from prepare_metadata import detect_group_labels, load_or_create
from utils import ensure_dirs, generate_ciriquant_config, load_config, run_snakemake, save_config


# ── Sub-commands ──────────────────────────────────────────────────────────────

def cmd_run(args: argparse.Namespace) -> None:
    """Full pipeline: metadata → groups → snakemake."""
    ensure_dirs()

    cfg = load_config("config.yaml")

    # 1. Resolve metadata
    print("\n── Step 1: Metadata ─────────────────────────────────", flush=True)
    metadata_df = prepare_metadata(
        gse_id      = args.gse,
        runinfo_csv = args.runinfo,
        output_file = cfg["metadata"],
    )
    print(f"   {len(metadata_df)} samples found", flush=True)

    # 2. Update project_id in config (works for GSE / PRJNA / SRP)
    if args.gse:
        cfg["project_id"] = args.gse
        cfg["results_dir"] = f"results/{args.gse}"
    cfg["threads"] = args.cores

    # 2b. Auto-detect case/control labels from GEO metadata
    print("\n── Step 1b: Detect group labels ─────────────────────", flush=True)
    case_label, control_label = detect_group_labels(metadata_df)
    existing_case    = cfg.get("de", {}).get("tumor_label", "tumor")
    existing_control = cfg.get("de", {}).get("normal_label", "normal")
    if case_label != existing_case or control_label != existing_control:
        cfg.setdefault("de", {})["tumor_label"]  = case_label
        cfg.setdefault("de", {})["normal_label"]  = control_label
        print(f"   [detect] case={case_label!r}  control={control_label!r}  → config.yaml updated", flush=True)
    else:
        print(f"   [detect] case={case_label!r}  control={control_label!r}  (unchanged)", flush=True)

    save_config(cfg)

    # 3. Assign case/control groups
    print("\n── Step 2: Sample grouping ──────────────────────────", flush=True)
    load_or_create(
        metadata_file  = cfg["metadata"],
        groups_file    = cfg["groups"],
        interactive    = args.interactive,
        case_label     = case_label,
        control_label  = control_label,
    )

    # 4. Generate CIRIquant config
    print("\n── Step 3: CIRIquant config ─────────────────────────", flush=True)
    generate_ciriquant_config(cfg, output=cfg["ciriquant_config"])

    # 5. Launch Snakemake
    print("\n── Step 4: Snakemake pipeline ───────────────────────", flush=True)
    extra = []
    if args.cluster:
        extra += ["--cluster", args.cluster, "--jobs", str(args.jobs)]
    if args.profile:
        extra += ["--profile", args.profile]

    rc = run_snakemake(
        snakefile  = "workflow/Snakefile",
        cores      = args.cores,
        extra_args = extra,
        dry_run    = args.dry_run,
    )
    sys.exit(rc)


def cmd_setup_ciriquant(_: argparse.Namespace) -> None:
    """Only regenerate the CIRIquant config YAML (useful after editing config.yaml)."""
    cfg = load_config("config.yaml")
    generate_ciriquant_config(cfg, output=cfg["ciriquant_config"])


def cmd_grouping(args: argparse.Namespace) -> None:
    """Re-run (or force-interactive) sample grouping without re-downloading."""
    cfg = load_config("config.yaml")
    # Remove existing groups file to force re-creation
    gf = Path(cfg["groups"])
    if gf.exists() and args.reset:
        gf.unlink()
    load_or_create(
        metadata_file = cfg["metadata"],
        groups_file   = cfg["groups"],
        interactive   = True,
    )


# ── Argument parser ───────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog        = "circRNA-agent",
        description = "Automated circRNA analysis pipeline (GEO / SRA BioProject → HTML report)",
        formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command")

    # ── run ──
    run_p = sub.add_parser("run", help="Run the full pipeline")
    src = run_p.add_mutually_exclusive_group(required=True)
    src.add_argument("--gse",     metavar="ACCESSION",
                     help="GEO or SRA accession: GSE113230 / PRJNA808398 / SRP156355")
    src.add_argument("--runinfo", metavar="CSV",    help="Path to NCBI SRA RunInfo.csv")
    run_p.add_argument("--cores",       type=int,  default=8,    help="CPU cores (default: 8)")
    run_p.add_argument("--interactive", action="store_true",     help="Manually assign tumor/normal labels")
    run_p.add_argument("--dry-run",     action="store_true",     help="Show Snakemake plan without executing")
    run_p.add_argument("--cluster",     metavar="CMD",           help="Snakemake --cluster submit command")
    run_p.add_argument("--jobs",        type=int,  default=50,   help="Max concurrent cluster jobs")
    run_p.add_argument("--profile",     metavar="DIR",           help="Snakemake profile directory")
    run_p.set_defaults(func=cmd_run)

    # ── setup-ciriquant ──
    sc_p = sub.add_parser("setup-ciriquant", help="Regenerate config/ciriquant.yaml")
    sc_p.set_defaults(func=cmd_setup_ciriquant)

    # ── grouping ──
    gr_p = sub.add_parser("grouping", help="(Re-)assign tumor/normal labels interactively")
    gr_p.add_argument("--reset", action="store_true", help="Remove existing groups file first")
    gr_p.set_defaults(func=cmd_grouping)

    return parser


def main() -> None:
    parser = build_parser()

    # Allow legacy positional-style: `python agent.py --gse GSE123`
    # by routing bare --gse/--runinfo flags to the 'run' subcommand
    argv = sys.argv[1:]
    if argv and argv[0] not in ("run", "setup-ciriquant", "grouping", "--help", "-h", "--version"):
        argv = ["run"] + argv

    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        sys.exit(1)
    args.func(args)


if __name__ == "__main__":
    main()
