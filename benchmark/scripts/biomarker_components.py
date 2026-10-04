"""Six biomarker components and the weighted score, reproduced from scripts/rank_biomarkers.py.

Components (same definitions and normalisation populations as _compute_weighted_score):
  S  = min(-log10 p, 10)            min-max over the significant set
  F  = min(|log2FC|, 5)             min-max over the significant set
  C  = max(confidence_score, 0)     min-max over the significant set
  K  = in_circbase (0/1)            not normalised
  M  = n_mirna                      min-max over interaction-queried rows only; NaN elsewhere
  R  = n_rbp                        min-max over interaction-queried rows only; NaN elsewhere

Missing values are NaN, not zero. The candidate file stores n_mirna = n_rbp = 0 for rows that were
never queried; the interaction_queried flag is what marks them missing. This module reads the flag.

min-max follows rank_biomarkers._minmax: a constant column maps to 0.5.

Weighted score for a weight scheme w (dict over S,F,C,K,M,R):
  score = sum_c(w_c * x_c) / sum_c(w_c)   summed over components c that are valid (non-NaN) for that row.
The denominator is the weight of the components that actually take part for that row, not a fixed
sum of six. A fixed denominator would make a missing M/R and a dropped component interact.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

COMPONENTS = ["S", "F", "C", "K", "M", "R"]
SIG_CAP = 10.0
FC_CAP = 5.0


def _minmax(s: pd.Series) -> pd.Series:
    valid = s.dropna()
    if valid.empty:
        return pd.Series(np.nan, index=s.index)
    lo, hi = valid.min(), valid.max()
    if hi == lo:
        out = pd.Series(np.nan, index=s.index)
        out.loc[valid.index] = 0.5
        return out
    return (s - lo) / (hi - lo)


def components(cand: pd.DataFrame) -> pd.DataFrame:
    """Normalised components for every row of a v2 candidate table (one primary method)."""
    p = cand["pvalue"].astype(float)
    s_raw = p.map(lambda x: min(-math.log10(max(x, 1e-10)), SIG_CAP))
    f_raw = cand["log2FC"].astype(float).abs().clip(upper=FC_CAP)
    c_raw = cand["confidence_score"].astype(float).clip(lower=0)
    queried = cand["interaction_queried"].astype(str).isin(["True", "true", "1"])
    m_raw = cand["n_mirna"].astype(float).where(queried)
    r_raw = cand["n_rbp"].astype(float).where(queried)
    out = pd.DataFrame(index=cand.index)
    out["S"] = _minmax(s_raw)
    out["F"] = _minmax(f_raw)
    out["C"] = _minmax(c_raw)
    out["K"] = cand["in_circbase"].astype(float)  # 0/1, never missing
    out["M"] = _minmax(m_raw)
    out["R"] = _minmax(r_raw)
    return out


def weighted_score(comp: pd.DataFrame, weights: dict) -> pd.Series:
    """Renormalised weighted mean over the components valid for each row (see module docstring)."""
    w = pd.Series({c: float(weights.get(c, 0.0)) for c in COMPONENTS})
    valid = comp[COMPONENTS].notna()
    num = (comp[COMPONENTS].fillna(0.0) * w).sum(axis=1)
    den = (valid * w).sum(axis=1)
    return (num / den.replace(0, np.nan))
