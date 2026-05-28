"""metaforge.common — distributions, heterogeneity estimators, pooling.

Pure Python (standard library only). Optional numpy acceleration is used
by the NMA module but is NOT required here. These primitives are shared by
every analysis type (pairwise, proportions, NMA, dose-response).

Methodology standard (consistent across metaforge):
  * tau^2: Paule-Mandel (default) or DerSimonian-Laird or REML
  * Knapp-Hartung (HKSJ) confidence interval with a floor at 1
  * Prediction interval via t_(k-1) per Cochrane Handbook v6.5 (k>=3)
"""
from __future__ import annotations
import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

Z975 = 1.959964  # qnorm(0.975)


# ---------------------------------------------------------------------------
# Distributions
# ---------------------------------------------------------------------------
def norm_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def norm_sf(z: float) -> float:
    return 1.0 - norm_cdf(z)


def norm_ppf(p: float) -> float:
    """Inverse normal CDF (Acklam). Good to ~1e-9 in the central region."""
    if p <= 0.0:
        return -math.inf
    if p >= 1.0:
        return math.inf
    a = [-3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2,
         1.383577518672690e2, -3.066479806614716e1, 2.506628277459239e0]
    b = [-5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2,
         6.680131188771972e1, -1.328068155288572e1]
    c = [-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838e0,
         -2.549732539343734e0, 4.374664141464968e0, 2.938163982698783e0]
    d = [7.784695709041462e-3, 3.224671290700398e-1, 2.445134137142996e0,
         3.754408661907416e0]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


_T_TABLE = {
    1: 12.7062, 2: 4.30265, 3: 3.18245, 4: 2.77645, 5: 2.57058, 6: 2.44691,
    7: 2.36462, 8: 2.306, 9: 2.26216, 10: 2.22814, 11: 2.20099, 12: 2.1788,
    13: 2.16037, 14: 2.14479, 15: 2.13145, 16: 2.11991, 17: 2.10982,
    18: 2.10092, 19: 2.09302, 20: 2.08596, 21: 2.07961, 22: 2.07387,
    23: 2.06866, 24: 2.0639, 25: 2.05954, 26: 2.05553, 27: 2.05183,
    28: 2.04841, 29: 2.04523, 30: 2.04227, 40: 2.02108, 50: 2.00856,
    60: 2.0003, 80: 1.99006, 100: 1.98397, 120: 1.97993,
}


def t_ppf975(df: int) -> float:
    """Two-sided 97.5% Student-t quantile."""
    if df <= 0:
        return float("nan")
    if df in _T_TABLE:
        return _T_TABLE[df]
    if df > 120:
        return Z975 + (Z975**3 + Z975) / (4 * df)  # near-normal
    z = Z975
    g1 = (z**3 + z) / 4
    g2 = (5*z**5 + 16*z**3 + 3*z) / 96
    g3 = (3*z**7 + 19*z**5 + 17*z**3 - 15*z) / 384
    g4 = (79*z**9 + 776*z**7 + 1482*z**5 - 1920*z**3 - 945*z) / 92160
    return z + g1/df + g2/df**2 + g3/df**3 + g4/df**4


def _gammaincc(s: float, x: float) -> float:
    if x <= 0 or s <= 0:
        return 1.0
    if x < s + 1:
        ap, summ, term = s, 1.0 / s, 1.0 / s
        for _ in range(300):
            ap += 1
            term *= x / ap
            summ += term
            if abs(term) < abs(summ) * 1e-14:
                break
        return 1.0 - summ * math.exp(-x + s * math.log(x) - math.lgamma(s))
    b, c, d, h = x + 1 - s, 1e300, 0.0, 0.0
    d = 1.0 / b
    h = d
    for i in range(1, 300):
        an = -i * (i - s)
        b += 2
        d = an * d + b
        if abs(d) < 1e-30:
            d = 1e-30
        c = b + an / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delt = d * c
        h *= delt
        if abs(delt - 1.0) < 1e-14:
            break
    return h * math.exp(-x + s * math.log(x) - math.lgamma(s))


def chi2_sf(x: float, df: int) -> float:
    """Upper-tail chi-square survival function."""
    if df <= 0 or x < 0:
        return float("nan")
    return _gammaincc(df / 2.0, x / 2.0)


# ---------------------------------------------------------------------------
# Heterogeneity (tau^2)
# ---------------------------------------------------------------------------
def tau2_paule_mandel(yis: Sequence[float], vis: Sequence[float]) -> float:
    """Paule-Mandel via bisection (monotone, always converges)."""
    k = len(yis)
    if k < 2:
        return 0.0

    def genQ(t2: float) -> float:
        ws = [1.0 / (v + t2) for v in vis]
        sw = sum(ws)
        mu = sum(w * y for w, y in zip(ws, yis)) / sw
        return sum(w * (y - mu) ** 2 for w, y in zip(ws, yis))

    target = k - 1
    if genQ(0.0) <= target:
        return 0.0
    hi = 1.0
    for _ in range(200):
        if genQ(hi) < target:
            break
        hi *= 2.0
    else:
        return hi
    lo = 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        q = genQ(mid)
        if abs(q - target) < 1e-12:
            return mid
        if q > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def tau2_dersimonian_laird(yis: Sequence[float], vis: Sequence[float]) -> float:
    k = len(yis)
    if k < 2:
        return 0.0
    w = [1.0 / v for v in vis]
    sw = sum(w)
    mu = sum(wi * y for wi, y in zip(w, yis)) / sw
    Q = sum(wi * (y - mu) ** 2 for wi, y in zip(w, yis))
    sw2 = sum(wi * wi for wi in w)
    c = sw - sw2 / sw
    if c <= 0:
        return 0.0
    return max(0.0, (Q - (k - 1)) / c)


def _reml_loglik(t2: float, yis, vis) -> float:
    """Restricted log-likelihood (up to a constant) for the intercept-only
    random-effects model at a given tau^2."""
    w = [1.0 / (v + t2) for v in vis]
    sw = sum(w)
    mu = sum(wi * y for wi, y in zip(w, yis)) / sw
    ll = -0.5 * sum(math.log(v + t2) for v in vis)
    ll += -0.5 * math.log(sw)
    ll += -0.5 * sum(wi * (y - mu) ** 2 for wi, y in zip(w, yis))
    return ll


def tau2_reml(yis: Sequence[float], vis: Sequence[float]) -> float:
    """REML by directly maximising the restricted log-likelihood over
    tau^2 >= 0 (golden-section search). Robust and formula-free; the
    boundary tau^2 = 0 is handled explicitly."""
    k = len(yis)
    if k < 2:
        return 0.0
    # Upper bound: a few times the total observed variance spread.
    spread = max(yis) - min(yis)
    hi = max(1e-6, (spread ** 2) + max(vis)) * 4 + 1.0
    # If the likelihood is decreasing from 0, the MLE is at the boundary.
    if _reml_loglik(1e-9, yis, vis) >= _reml_loglik(1e-6, yis, vis):
        # check a slightly larger step to avoid numerical noise at 0
        if _reml_loglik(0.0 + 1e-9, yis, vis) >= _reml_loglik(hi * 1e-3, yis, vis):
            return 0.0
    gr = (math.sqrt(5) - 1) / 2
    a, b = 0.0, hi
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc = _reml_loglik(c, yis, vis)
    fd = _reml_loglik(d, yis, vis)
    for _ in range(200):
        if fc < fd:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = _reml_loglik(d, yis, vis)
        else:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = _reml_loglik(c, yis, vis)
        if abs(b - a) < 1e-10:
            break
    best = 0.5 * (a + b)
    # Compare with the boundary; REML can peak at 0.
    return best if _reml_loglik(best, yis, vis) > _reml_loglik(0.0 + 1e-12, yis, vis) else 0.0


TAU2_METHODS = {
    "PM": tau2_paule_mandel,
    "DL": tau2_dersimonian_laird,
    "REML": tau2_reml,
}


# ---------------------------------------------------------------------------
# Generic inverse-variance random-effects pool
# ---------------------------------------------------------------------------
@dataclass
class Pool:
    k: int
    estimate: float          # on the analysis (input) scale
    se: float
    ci_low: float
    ci_high: float
    pi_low: Optional[float]
    pi_high: Optional[float]
    tau2: float
    i2: float
    h2: float
    Q: float
    Q_df: int
    Q_p: float
    z: float
    p: float
    weights: list
    hksj_floor_applied: bool
    tau2_method: str
    extra: dict = field(default_factory=dict)


def pool_inverse_variance(yis: Sequence[float], vis: Sequence[float],
                          tau2_method: str = "PM",
                          knha: bool = True) -> Optional[Pool]:
    """Random-effects inverse-variance pool on the analysis scale.

    yis/vis are effect sizes and variances already on the scale to be
    pooled (e.g. log OR, MD, logit p, Fisher z). Back-transformation is
    the caller's responsibility.
    """
    k = len(yis)
    if k == 0:
        return None
    if k == 1:
        se = math.sqrt(vis[0])
        return Pool(k=1, estimate=yis[0], se=se,
                    ci_low=yis[0] - Z975 * se, ci_high=yis[0] + Z975 * se,
                    pi_low=None, pi_high=None, tau2=0.0, i2=0.0, h2=1.0,
                    Q=0.0, Q_df=0, Q_p=float("nan"),
                    z=yis[0] / se if se else float("nan"),
                    p=2 * norm_sf(abs(yis[0] / se)) if se else float("nan"),
                    weights=[1.0], hksj_floor_applied=False,
                    tau2_method=tau2_method)

    est_fn = TAU2_METHODS.get(tau2_method, tau2_paule_mandel)
    tau2 = est_fn(yis, vis)
    w = [1.0 / (v + tau2) for v in vis]
    sw = sum(w)
    mu = sum(wi * y for wi, y in zip(w, yis)) / sw
    var_mu = 1.0 / sw

    # Q, I^2, H^2 (fixed-effect weights)
    wf = [1.0 / v for v in vis]
    muf = sum(wi * y for wi, y in zip(wf, yis)) / sum(wf)
    Q = sum(wi * (y - muf) ** 2 for wi, y in zip(wf, yis))
    Q_df = k - 1
    Q_p = chi2_sf(Q, Q_df)
    i2 = max(0.0, (Q - Q_df) / Q) * 100 if Q > 0 else 0.0
    h2 = max(1.0, Q / Q_df) if Q_df > 0 else 1.0

    floor_applied = False
    if knha:
        scale_raw = sum(wi * (y - mu) ** 2 for wi, y in zip(w, yis)) / (k - 1)
        scale = max(1.0, scale_raw)
        floor_applied = scale_raw < 1.0
        se = math.sqrt(scale * var_mu)
        crit = t_ppf975(k - 1)
    else:
        se = math.sqrt(var_mu)
        crit = Z975

    ci_low, ci_high = mu - crit * se, mu + crit * se
    if k >= 3:
        pi_se = math.sqrt(tau2 + se ** 2)
        pcrit = t_ppf975(k - 1)
        pi_low, pi_high = mu - pcrit * pi_se, mu + pcrit * pi_se
    else:
        pi_low = pi_high = None

    z = mu / se if se > 0 else float("nan")
    p = 2 * norm_sf(abs(z)) if math.isfinite(z) else float("nan")
    weights = [100.0 * wi / sw for wi in w]

    return Pool(k=k, estimate=mu, se=se, ci_low=ci_low, ci_high=ci_high,
                pi_low=pi_low, pi_high=pi_high, tau2=tau2, i2=i2, h2=h2,
                Q=Q, Q_df=Q_df, Q_p=Q_p, z=z, p=p, weights=weights,
                hksj_floor_applied=floor_applied, tau2_method=tau2_method)
