"""metaforge.pairwise — pairwise meta-analysis for all common effect measures.

Binary:     OR, RR, RD  (from 2x2 counts)
Continuous: MD, SMD     (from mean/sd/n per arm)
Generic:    any         (from a precomputed effect + 95% CI)

Ratio measures are pooled on the log scale and back-transformed.
"""
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Optional
from . import common

RATIO = {"OR", "RR", "HR"}
DIFF = {"MD", "SMD", "RD"}


@dataclass
class Effect:
    name: str
    yi: float
    vi: float
    n: Optional[int] = None
    raw: dict = None


def _cc(tE, tN, cE, cN):
    """Continuity correction: add 0.5 to all cells only if any cell is 0."""
    cells = [tE, tN - tE, cE, cN - cE]
    if min(cells) == 0:
        return tE + 0.5, tN + 1.0, cE + 0.5, cN + 1.0, True
    return tE, tN, cE, cN, False


def effect_OR(name, tE, tN, cE, cN) -> Optional[Effect]:
    if tN <= 0 or cN <= 0 or tE < 0 or cE < 0 or tE > tN or cE > cN:
        return None
    a, n1, c, n0, _ = _cc(tE, tN, cE, cN)
    b, d = n1 - a, n0 - c
    try:
        yi = math.log((a * d) / (b * c))
        vi = 1/a + 1/b + 1/c + 1/d
    except (ValueError, ZeroDivisionError):
        return None
    return Effect(name, yi, vi, n=tN + cN,
                  raw={"tE": tE, "tN": tN, "cE": cE, "cN": cN})


def effect_RR(name, tE, tN, cE, cN) -> Optional[Effect]:
    if tN <= 0 or cN <= 0 or tE < 0 or cE < 0 or tE > tN or cE > cN:
        return None
    a, n1, c, n0, _ = _cc(tE, tN, cE, cN)
    try:
        yi = math.log((a / n1) / (c / n0))
        vi = 1/a - 1/n1 + 1/c - 1/n0
    except (ValueError, ZeroDivisionError):
        return None
    return Effect(name, yi, vi, n=tN + cN,
                  raw={"tE": tE, "tN": tN, "cE": cE, "cN": cN})


def effect_RD(name, tE, tN, cE, cN) -> Optional[Effect]:
    if tN <= 0 or cN <= 0 or tE < 0 or cE < 0 or tE > tN or cE > cN:
        return None
    p1, p0 = tE / tN, cE / cN
    yi = p1 - p0
    vi = p1 * (1 - p1) / tN + p0 * (1 - p0) / cN
    if vi <= 0:
        vi = 1e-9
    return Effect(name, yi, vi, n=tN + cN,
                  raw={"tE": tE, "tN": tN, "cE": cE, "cN": cN})


def effect_MD(name, m1, sd1, n1, m0, sd0, n0) -> Optional[Effect]:
    if n1 < 2 or n0 < 2 or sd1 < 0 or sd0 < 0:
        return None
    yi = m1 - m0
    vi = sd1 ** 2 / n1 + sd0 ** 2 / n0
    if vi <= 0:
        return None
    return Effect(name, yi, vi, n=n1 + n0,
                  raw={"m1": m1, "sd1": sd1, "n1": n1, "m0": m0, "sd0": sd0, "n0": n0})


def effect_SMD(name, m1, sd1, n1, m0, sd0, n0) -> Optional[Effect]:
    """Hedges' g (small-sample corrected SMD)."""
    if n1 < 2 or n0 < 2:
        return None
    sp2 = ((n1 - 1) * sd1 ** 2 + (n0 - 1) * sd0 ** 2) / (n1 + n0 - 2)
    if sp2 <= 0:
        return None
    d = (m1 - m0) / math.sqrt(sp2)
    J = 1 - 3 / (4 * (n1 + n0 - 2) - 1)  # Hedges correction
    g = J * d
    vi = (n1 + n0) / (n1 * n0) + g ** 2 / (2 * (n1 + n0 - 2))
    return Effect(name, g, vi, n=n1 + n0,
                  raw={"m1": m1, "sd1": sd1, "n1": n1, "m0": m0, "sd0": sd0, "n0": n0})


def effect_generic(name, effect, ci_low, ci_high, measure) -> Optional[Effect]:
    if measure in RATIO:
        if effect <= 0 or ci_low <= 0 or ci_high <= 0:
            return None
        yi = math.log(effect)
        se = (math.log(ci_high) - math.log(ci_low)) / (2 * common.Z975)
    else:
        yi = effect
        se = (ci_high - ci_low) / (2 * common.Z975)
    if se <= 0 or not math.isfinite(se):
        return None
    return Effect(name, yi, se * se,
                  raw={"effect": effect, "ci_low": ci_low, "ci_high": ci_high})


def build_effects(trials, measure):
    """Build Effect objects from a list of trial dicts based on `measure`."""
    effects = []
    skipped = []
    for i, t in enumerate(trials):
        name = t.get("name") or t.get("nct") or f"Study {i+1}"
        e = None
        if measure == "OR" and all(k in t for k in ("tE", "tN", "cE", "cN")):
            e = effect_OR(name, t["tE"], t["tN"], t["cE"], t["cN"])
        elif measure == "RR" and all(k in t for k in ("tE", "tN", "cE", "cN")):
            e = effect_RR(name, t["tE"], t["tN"], t["cE"], t["cN"])
        elif measure == "RD" and all(k in t for k in ("tE", "tN", "cE", "cN")):
            e = effect_RD(name, t["tE"], t["tN"], t["cE"], t["cN"])
        elif measure == "MD" and all(k in t for k in ("m1", "sd1", "n1", "m0", "sd0", "n0")):
            e = effect_MD(name, t["m1"], t["sd1"], t["n1"], t["m0"], t["sd0"], t["n0"])
        elif measure == "SMD" and all(k in t for k in ("m1", "sd1", "n1", "m0", "sd0", "n0")):
            e = effect_SMD(name, t["m1"], t["sd1"], t["n1"], t["m0"], t["sd0"], t["n0"])
        elif all(k in t for k in ("effect", "ci_low", "ci_high")):
            e = effect_generic(name, t["effect"], t["ci_low"], t["ci_high"], measure)
        if e is None:
            skipped.append(name)
        else:
            e.raw = e.raw or {}
            e.raw["nct"] = t.get("nct")
            e.raw["year"] = t.get("year")
            effects.append(e)
    return effects, skipped


def analyze(trials, measure, tau2_method="PM"):
    """Run a pairwise meta-analysis. Returns (pool, effects, skipped) with
    display-scale fields attached to pool.extra."""
    effects, skipped = build_effects(trials, measure)
    if len(effects) < 1:
        return None, effects, skipped
    pool = common.pool_inverse_variance(
        [e.yi for e in effects], [e.vi for e in effects],
        tau2_method=tau2_method, knha=True)
    ratio = measure in RATIO
    def disp(v):
        return math.exp(v) if (ratio and v is not None) else v
    pool.extra = {
        "measure": measure,
        "ratio": ratio,
        "display": {
            "estimate": disp(pool.estimate),
            "ci_low": disp(pool.ci_low), "ci_high": disp(pool.ci_high),
            "pi_low": disp(pool.pi_low), "pi_high": disp(pool.pi_high),
        },
        "per_study": [
            {"name": e.name,
             "est": disp(e.yi),
             "lo": disp(e.yi - common.Z975 * math.sqrt(e.vi)),
             "hi": disp(e.yi + common.Z975 * math.sqrt(e.vi)),
             "weight": pool.weights[i],
             "raw": e.raw}
            for i, e in enumerate(effects)
        ],
    }
    return pool, effects, skipped
