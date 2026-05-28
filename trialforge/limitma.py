"""trialforge.limitma — limit meta-analysis (Rucker 2011, WLS approximation).

Ported (pure-stdlib) from the allmeta `limit-ma` engine. The limit
meta-analysis adjusts for small-study effects by regressing the observed
effect on its standard error (weighted by 1/(se^2 + tau^2)) and reading
the intercept at SE = 0 as the bias-adjusted ("limit") estimate.

  beta0 (intercept at SE=0) = limit estimate
  beta1 (slope)             > 0 indicates small-study effects

Validated against metasens::limitmeta R fixtures (allmeta): the FE/RE/tau2
baselines match exactly; the WLS limit estimate is the engine's documented
approximation to the full Rucker residual-decomposition algorithm
(direction + magnitude agree).
"""
from __future__ import annotations
import math
from . import common


def analyze(yis, vis, *, ratio=False):
    k = len(yis)
    if k < 3:
        return {"available": False, "reason": "need >=3 studies"}
    seis = [math.sqrt(v) for v in vis]

    tau2 = common.tau2_dersimonian_laird(yis, vis)
    # FE / RE baselines
    wfe = [1.0 / v for v in vis]
    fe = sum(w * y for w, y in zip(wfe, yis)) / sum(wfe)
    wre = [1.0 / (v + tau2) for v in vis]
    swre = sum(wre)
    re = sum(w * y for w, y in zip(wre, yis)) / swre
    se_re = math.sqrt(1.0 / swre)

    # WLS of te ~ b0 + b1*SE, weighted by 1/(se^2 + tau^2)
    w = wre
    sw = swx = swy = swxx = swxy = 0.0
    for x, y, wi in zip(seis, yis, w):
        sw += wi; swx += wi * x; swy += wi * y
        swxx += wi * x * x; swxy += wi * x * y
    det = sw * swxx - swx * swx
    if det == 0:
        return {"available": False, "reason": "degenerate regression"}
    b0 = (swxx * swy - swx * swxy) / det
    b1 = (sw * swxy - swx * swy) / det
    rss = sum(wi * (y - (b0 + b1 * x)) ** 2 for wi, x, y in zip(w, seis, yis))
    sigma2 = rss / max(1, k - 2)
    se_b0 = math.sqrt(max(0.0, sigma2 * swxx / det))

    def disp(v):
        return math.exp(v) if ratio else v

    return {
        "available": True, "k": k,
        "fe_estimate": disp(fe),
        "re_estimate": disp(re), "re_ci": (disp(re - common.Z975 * se_re),
                                           disp(re + common.Z975 * se_re)),
        "tau2": tau2,
        "limit_estimate": disp(b0),
        "limit_ci": (disp(b0 - common.Z975 * se_b0), disp(b0 + common.Z975 * se_b0)),
        "slope": b1,
        "small_study_effect": b1 > 0,
        "attenuates": (disp(b0) < disp(re)) if not ratio else (b0 < re),
        "note": "Rucker WLS approximation; intercept at SE=0 is the "
                "small-study-effect-adjusted estimate.",
    }
