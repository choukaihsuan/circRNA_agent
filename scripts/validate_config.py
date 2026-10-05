"""
validate_config.py -- Pre-flight validation and default-filling for config.yaml.

Called from workflow/Snakefile BEFORE any rule (including download) is
evaluated, so a bad value fails fast with a clear message instead of
crashing hours into a run.

Design contract (circDEX parameter-exposure project, 2026-10):
  1. Every default documented here is the value this pipeline has always
     used for its production runs (CIRIquant/DCC hardcoded values, or
     config.yaml values already shipped). Loading an old config.yaml /
     config/projects/*.yaml that is missing any of these keys MUST
     reproduce the exact same shell commands as before -- see
     tests/test_validate_config.py (added in Phase 5) for the
     byte-for-byte check.
  2. Strand-related parameters (DCC -N/-ss, CIRIquant -l/--library-type)
     are intentionally NOT exposed here yet -- see CLAUDE.md "strand
     assumption investigation" section. Do not add them without reading
     that section first.
  3. `tools.dcc.threads` default is 4 to match the value this pipeline's
     production DCC rule has always used (circrna.smk, previously
     hardcoded `threads: 4`, now config-driven but same default).
     benchmark/config_benchmark.yaml pins this to 4 explicitly and
     independently, because Figure 6(c)'s compute-cost numbers were
     measured at -T 4; do not let that file inherit this default blindly.

Usage (standalone CLI, same checks the Snakefile import uses):
  python scripts/validate_config.py config.yaml
  python scripts/validate_config.py config/projects/GSE113230.yaml
"""
from __future__ import annotations

import sys
from typing import Any, Dict, List, Tuple

import yaml


class ConfigError(Exception):
    """Raised for any invalid config value; message is meant to be read by a human."""


# ── Default values -- MUST match current production behaviour exactly ──────
# (tool, key) -> (default, type, validator)
# validator(value) -> None, raises ConfigError with a clear message on failure.

def _bool(v: Any) -> bool:
    if isinstance(v, bool):
        return v
    raise ConfigError(f"expected true/false, got {v!r}")


def _int_range(lo: int = None, hi: int = None):
    def f(v: Any) -> int:
        if isinstance(v, bool) or not isinstance(v, int):
            raise ConfigError(f"expected an integer, got {v!r}")
        if lo is not None and v < lo:
            raise ConfigError(f"must be >= {lo}, got {v}")
        if hi is not None and v > hi:
            raise ConfigError(f"must be <= {hi}, got {v}")
        return v
    return f


def _float_range(lo: float = None, hi: float = None, allow_eq_lo: bool = True):
    def f(v: Any) -> float:
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            raise ConfigError(f"expected a number, got {v!r}")
        v = float(v)
        if lo is not None:
            if allow_eq_lo and v < lo:
                raise ConfigError(f"must be >= {lo}, got {v}")
            if not allow_eq_lo and v <= lo:
                raise ConfigError(f"must be > {lo}, got {v}")
        if hi is not None and v > hi:
            raise ConfigError(f"must be <= {hi}, got {v}")
        return v
    return f


def _choice(*options: str):
    def f(v: Any) -> str:
        if v not in options:
            raise ConfigError(f"must be one of {options}, got {v!r}")
        return v
    return f


# Schema: section -> {key: (default, validator)}
# Only keys this project has chosen to expose. Any key not listed here is
# passed through untouched (e.g. genome.*, metadata, groups, project_id).
SCHEMA: Dict[str, Dict[str, Tuple[Any, Any]]] = {
    "tools.ciriquant": {
        "no_gene": (False, _bool),                 # --no-gene: skip stringtie gene abundance (speed only, no effect on circRNA calls)
        "anchor":  (5, _int_range(lo=1)),           # -a/--anchor: min anchor length for junction alignment
    },
    "tools.dcc": {
        "threads":             (4, _int_range(lo=1, hi=96)),   # -T; was hardcoded 4 in circrna.smk
        "count_threshold":     (2, _int_range(lo=1)),          # -Nr arg 1
        "replicate_threshold": (1, _int_range(lo=1)),          # -Nr arg 2
        "filter_mito":         (True, _bool),                  # -M
        "filter_recommended":  (True, _bool),                  # -F
        "end_tol":             (5, _int_range(lo=0, hi=9)),    # -E
        "max_length":          (1_000_000, _int_range(lo=1)), # -m
        "min_length":          (30, _int_range(lo=1)),         # -n
    },
    "consensus": {
        "min_tools":          (2, _int_range(lo=1, hi=2)),
        "slop":                (10, _int_range(lo=0)),
        "min_bsj_reads":       (2, _int_range(lo=1)),
        "max_junction_ratio":  (1.0, _float_range(lo=0.0, allow_eq_lo=False)),
        "qc_bsj_threshold":    (5, _int_range(lo=0)),
        # NOTE: default True to match production -- circrna.smk's previous
        # `.get("adaptive", True)` fallback is what every published dataset
        # actually ran with. consensus_filter.py's own argparse default of
        # False is NOT what was used in practice; do not "fix" this to False.
        "adaptive":            (True, _bool),
        "adaptive_ratio":      (0.1, _float_range(lo=0.0, hi=1.0)),
        "min_confidence":      (0.0, _float_range(lo=0.0)),
    },
    "de": {
        "method":                      ("edgeR_ciriquant", _choice("edgeR_ciriquant", "deseq2", "limma")),
        "fdr_cutoff":                  (0.05, _float_range(lo=0.0, hi=1.0)),
        "log2fc_cutoff":               (1.0, _float_range(lo=0.0)),
        "de_sig_by":                   ("pvalue", _choice("auto", "pvalue", "padj", "qvalue")),
        "heatmap_top_n":               (10, _int_range(lo=1)),
        "fsj_concordance_lfc":         (0.0, _float_range(lo=0.0)),
        "isoform_fdr_cutoff":          (0.1, _float_range(lo=0.0, hi=1.0)),
        "delta_iui_cutoff":            (0.1, _float_range(lo=0.0)),
        "interaction_top_n":           (50, _int_range(lo=1)),
        "encori_clip_exp_num":         (1, _int_range(lo=0)),
        "encori_program_num":          (2, _int_range(lo=0)),
        # New in this phase -- previously hardcoded in analysis.R
        "filter_expr_min_count":       (5, _int_range(lo=1)),
        "filter_expr_min_total_count": (15, _int_range(lo=1)),
        "deseq2_size_factor_type":     ("poscounts", _choice("poscounts", "ratio", "iterate")),
        "paired_design_override":      ("auto", _choice("auto", "force_paired", "force_unpaired")),
    },
    "de.biomarker": {
        "sig_cap": (10.0, _float_range(lo=0.0, allow_eq_lo=False)),
        "fc_cap":  (5.0, _float_range(lo=0.0, allow_eq_lo=False)),
    },
}

# Biomarker component weights live one level deeper (de.biomarker.components.<name>)
# and need per-component validation plus whole-set degenerate checks, so they
# are handled separately from SCHEMA rather than forced into its flat shape.
BIOMARKER_COMPONENTS = ("sig", "fc", "conf", "known", "mirna", "rbp")


def _get_section(config: Dict[str, Any], dotted: str) -> Dict[str, Any]:
    """Walk dotted path, creating empty dicts as needed (mutates config)."""
    node = config
    for part in dotted.split("."):
        if part not in node or node[part] is None:
            node[part] = {}
        if not isinstance(node[part], dict):
            raise ConfigError(f"config.{dotted} must be a mapping, got {type(node[part]).__name__}")
        node = node[part]
    return node


def validate_and_fill_defaults(config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate every key in SCHEMA that is present; fill in the documented
    default for every key that is absent (backward compatibility with old
    config.yaml / config/projects/*.yaml files that predate this schema).

    Mutates and returns `config`. Raises ConfigError on the first invalid
    value found (message identifies the exact dotted key).
    """
    errors: List[str] = []
    for section, keys in SCHEMA.items():
        node = _get_section(config, section)
        for key, (default, validator) in keys.items():
            full_key = f"{section}.{key}"
            if key not in node or node[key] is None:
                node[key] = default
                continue
            try:
                node[key] = validator(node[key])
            except ConfigError as e:
                errors.append(f"  config.{full_key} = {node[key]!r}: {e}")

    _validate_cross_field(config, errors)
    _validate_biomarker_components(config, errors)

    if errors:
        raise ConfigError(
            "Invalid config.yaml -- fix before running the pipeline:\n" + "\n".join(errors)
        )
    return config


def _validate_biomarker_components(config: Dict[str, Any], errors: List[str]) -> None:
    """
    Validate de.biomarker.components.<name> = {enabled, weight}.
    Fills in the equal-weight default (this paper's setting) for any missing
    component. Rejects degenerate configurations outright -- see CLAUDE.md
    "六維生物標記分數" section for the rationale behind each rejection.
    """
    comp = _get_section(config, "de.biomarker.components")
    for name in BIOMARKER_COMPONENTS:
        if name not in comp or comp[name] is None:
            comp[name] = {"enabled": True, "weight": 1.0}
            continue
        c = comp[name]
        if not isinstance(c, dict):
            errors.append(f"  config.de.biomarker.components.{name} must be a mapping, got {c!r}")
            continue
        c.setdefault("enabled", True)
        c.setdefault("weight", 1.0)
        try:
            c["enabled"] = _bool(c["enabled"])
        except ConfigError as e:
            errors.append(f"  config.de.biomarker.components.{name}.enabled = {c['enabled']!r}: {e}")
        try:
            # weight=0 is allowed (equivalent to disabling without unchecking
            # the box); negative weights are not.
            c["weight"] = _float_range(lo=0.0)(c["weight"])
        except ConfigError as e:
            errors.append(f"  config.de.biomarker.components.{name}.weight = {c['weight']!r}: {e}")

    if errors:
        return  # per-component errors already explain the problem; skip set-level checks on bad data

    enabled = {name: comp[name]["weight"] for name in BIOMARKER_COMPONENTS if comp[name]["enabled"]}
    if not enabled:
        errors.append(
            "  config.de.biomarker.components: all six components are disabled -- "
            "at least one must be enabled to compute a Biomarker Score"
        )
        return
    if sum(enabled.values()) <= 0:
        errors.append(
            "  config.de.biomarker.components: every enabled component has weight 0 "
            f"(enabled: {list(enabled)}) -- the score would be undefined (0/0)"
        )
        return
    if set(enabled) == {"known"}:
        errors.append(
            "  config.de.biomarker.components: only 'known' (in_circbase) is enabled. "
            "known_bonus is a binary 0/1 value -- ranking on it alone ties every "
            "circBase-listed candidate at the same score with no real ordering "
            "(the apparent order would just reflect row input order, not biology). "
            "Enable at least one more component, or use the circbase_id column "
            "directly instead of the Biomarker Score for this purpose."
        )


def _validate_cross_field(config: Dict[str, Any], errors: List[str]) -> None:
    """Checks that need more than one field at a time."""
    dcc = config.get("tools", {}).get("dcc", {})
    if dcc.get("min_length", 30) >= dcc.get("max_length", 1_000_000):
        errors.append(
            f"  tools.dcc.min_length ({dcc.get('min_length')}) must be < "
            f"tools.dcc.max_length ({dcc.get('max_length')})"
        )


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(f"Usage: python {sys.argv[0]} <config.yaml>")
    path = sys.argv[1]
    with open(path) as f:
        config = yaml.safe_load(f) or {}
    try:
        validate_and_fill_defaults(config)
    except ConfigError as e:
        sys.exit(f"[validate_config] FAILED: {e}")
    print(f"[validate_config] OK: {path}")


if __name__ == "__main__":
    main()
