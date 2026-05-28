"""trialforge.copas — Copas selection-model publication-bias sensitivity.

Ported (pure-stdlib re-implementation) from the allmeta `copas` engine.
The Copas model posits that a study is published with probability
    p_i = Phi(a - gamma / se_i)
so less-precise (small) studies are preferentially published when their
effect is large. We profile the adjusted pooled estimate over a grid of
the selection-strength parameter gamma (equivalently, an assumed
proportion of unpublished studies), using the Heckman-Tobit (HT)
inverse-probability reweighting approximation:
    w_i = (1/se_i^2) / p_i,   pooled = sum(w_i y_i)/sum(w_i),  se = sqrt(1/sum w_i)

This is the fast HT approximation (not the full bivariate MLE of
metasens::copas), used for a *sensitivity* read on how far the pooled
estimate could move under increasing assumed publication bias.

advanced-stats.md: Copas needs k>=15 for stable MLE; the HT sensitivity
profile is informative at smaller k but should be read as directional.
"""
from __future__ import annotations
import math
from . import common


def _ht_pool(yis, seis, gamma, p_unobs):
    """HT-reweighted pool at a given selection strength."""
    if gamma <= 0:
        w = [1.0 / (s * s) for s in seis]
        sw = sum(w)
        mu = sum(wi * y for wi, y in zip(w, yis)) / sw
        return mu, math.sqrt(1.0 / sw), [1.0] * len(seis)
    # binary search intercept a so mean(p_i) = 1 - p_unobs
    lo, hi = -8.0, 8.0
    target = 1.0 - p_unobs
    for _ in range(60):
        a = 0.5 * (lo + hi)
        mean_p = sum(common.norm_cdf(a - gamma / s) for s in seis) / len(seis)
        if mean_p < target:
            lo = a
        else:
            hi = a
    a = 0.5 * (lo + hi)
    p_sel = [max(0.01, common.norm_cdf(a - gamma / s)) for s in seis]
    w = [(1.0 / (s * s)) / p for s, p in zip(seis, p_sel)]
    sw = sum(w)
    mu = sum(wi * y for wi, y in zip(w, yis)) / sw
    return mu, math.sqrt(1.0 / sw), p_sel


def analyze(yis, vis, *, ratio=False, p_unobs_grid=None, gamma=1.0):
    """Copas sensitivity profile.

    Returns unadjusted FE/RE pools and adjusted estimates across a grid of
    assumed proportions of unpublished studies (default 5%..50%).
    """
    k = len(yis)
    if k < 3:
        return {"available": False, "reason": "need >=3 studies"}
    seis = [math.sqrt(v) for v in vis]

    # True fixed-effect pool (tau^2 = 0) is the Copas baseline.
    _w = [1.0 / v for v in vis]
    fe_estimate_raw = sum(wi * y for wi, y in zip(_w, yis)) / sum(_w)
    re = common.pool_inverse_variance(yis, vis, tau2_method="DL")

    if p_unobs_grid is None:
        p_unobs_grid = [0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50]

    def disp(v):
        return math.exp(v) if ratio else v

    profile = []
    for pu in p_unobs_grid:
        mu, se, _ = _ht_pool(yis, seis, gamma if pu > 0 else 0.0, pu)
        profile.append({
            "p_unpublished": pu,
            "estimate": disp(mu),
            "ci_low": disp(mu - common.Z975 * se),
            "ci_high": disp(mu + common.Z975 * se),
        })
    # worst case = largest assumed unpublished fraction
    worst = profile[-1]
    # slope of estimate vs assumed unpublished fraction (sensitivity)
    xs = [p["p_unpublished"] for p in profile]
    ys = [(math.log(p["estimate"]) if ratio else p["estimate"]) for p in profile]
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    denom = sum((x - mx) ** 2 for x in xs)
    slope = (sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / denom) if denom else 0.0

    return {
        "available": True, "k": k,
        "fe_estimate": disp(fe_estimate_raw),
        "re_estimate": disp(re.estimate),
        "re_ci": (disp(re.ci_low), disp(re.ci_high)),
        "profile": profile,
        "worst_case": worst,
        "sensitivity_slope": slope,
        "attenuates": (worst["estimate"] < re.estimate) if not ratio
                      else (worst["estimate"] < disp(re.estimate)),
        "note": "HT approximation; read as directional sensitivity (full MLE "
                "needs k>=15).",
    }
