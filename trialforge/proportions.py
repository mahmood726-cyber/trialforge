"""metaforge.proportions — single-arm proportion / prevalence meta-analysis.

Transformations:
  logit  : pool log-odds, back-transform with the inverse logit
  DAS    : Freeman-Tukey double arcsine (variance-stabilising; handles 0/1
           proportions). Back-transformed with the harmonic-mean-of-n
           correction (Miller 1978) to reduce bias.
  raw    : pool proportions directly (only for quick checks)

Random effects throughout (Paule-Mandel tau^2 by default).
"""
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Optional
from . import common


@dataclass
class PropStudy:
    name: str
    e: int
    n: int
    yi: float
    vi: float


def _logit_effect(name, e, n):
    # Haldane-Anscombe 0.5 correction for 0/everything cells
    a = e + 0.5
    b = n - e + 0.5
    yi = math.log(a / b)
    vi = 1 / a + 1 / b
    return PropStudy(name, e, n, yi, vi)


def _das_effect(name, e, n):
    # Freeman-Tukey double arcsine transform
    yi = 0.5 * (math.asin(math.sqrt(e / (n + 1))) + math.asin(math.sqrt((e + 1) / (n + 1))))
    vi = 1.0 / (4 * n + 2)
    return PropStudy(name, e, n, yi, vi)


def _raw_effect(name, e, n):
    p = e / n
    vi = max(p * (1 - p) / n, 1e-9)
    return PropStudy(name, e, n, p, vi)


def _das_backtransform(t, n_harm):
    """Back-transform a pooled double-arcsine value using the harmonic mean
    of sample sizes (Miller 1978)."""
    # p = 0.5 * (1 - sign(cos2t) * sqrt(1 - (sin2t + (sin2t - 1/sin2t)/n)^2))
    # Use the standard Miller inverse:
    try:
        val = 0.5 * (1 - math.copysign(1, math.cos(2 * t)) *
                     math.sqrt(max(0.0, 1 - (math.sin(2 * t) +
                               (math.sin(2 * t) - 1 / math.sin(2 * t)) / n_harm) ** 2)))
    except (ValueError, ZeroDivisionError):
        val = math.sin(t) ** 2
    return min(1.0, max(0.0, val))


def analyze(studies, method="DAS", tau2_method="PM"):
    """Pool single-arm proportions. `studies` = list of {name, e, n}."""
    builders = {"logit": _logit_effect, "DAS": _das_effect, "raw": _raw_effect}
    build = builders.get(method, _das_effect)
    ps = []
    skipped = []
    for i, s in enumerate(studies):
        name = s.get("name") or f"Study {i+1}"
        e, n = s.get("e"), s.get("n")
        if e is None or n is None or n <= 0 or e < 0 or e > n:
            skipped.append(name)
            continue
        ps.append(build(name, int(e), int(n)))
    if not ps:
        return None, ps, skipped

    pool = common.pool_inverse_variance([p.yi for p in ps], [p.vi for p in ps],
                                        tau2_method=tau2_method, knha=True)

    def back(v):
        if v is None:
            return None
        if method == "logit":
            return 1 / (1 + math.exp(-v))
        if method == "raw":
            return min(1.0, max(0.0, v))
        # DAS
        n_harm = len(ps) / sum(1.0 / p.n for p in ps)
        return _das_backtransform(v, n_harm)

    pool.extra = {
        "method": method,
        "display": {
            "estimate": back(pool.estimate),
            "ci_low": back(pool.ci_low), "ci_high": back(pool.ci_high),
            "pi_low": back(pool.pi_low), "pi_high": back(pool.pi_high),
        },
        "per_study": [
            {"name": p.name, "est": p.e / p.n,
             "lo": _wilson(p.e, p.n)[0], "hi": _wilson(p.e, p.n)[1],
             "weight": pool.weights[i], "raw": {"e": p.e, "n": p.n}}
            for i, p in enumerate(ps)
        ],
    }
    return pool, ps, skipped


def _wilson(e, n):
    """Wilson score interval for a single study's display CI."""
    if n == 0:
        return (0.0, 1.0)
    p = e / n
    z = common.Z975
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = (z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))
