"""trialforge.dta — diagnostic test accuracy meta-analysis.

Pure-stdlib re-implementation of the allmeta `hsroc` engine: a bivariate
DerSimonian-Laird approximation that pools logit(sensitivity) and
logit(false-positive rate) separately, with an empirical between-study
correlation, then summarises as a pooled Se/Sp operating point and an
SROC curve.

Per study from a 2x2 table (TP, FP, FN, TN):
  Se   = TP/(TP+FN);     logit(Se),   var = 1/TP + 1/FN
  FPR  = FP/(FP+TN);     logit(FPR),  var = 1/FP + 1/TN
  (continuity correction +0.5 to all cells only when any cell is 0)

Parameterisation matches mada::reitsma (logit FPR, NOT logit Spec).
The DL approximation differs from full bivariate REML — read the pooled
point as an approximation; rho is constrained to [-0.95, 0.95]
(advanced-stats.md). A Spearman threshold-effect check is reported: strong
negative logit(Se)-vs-logit(FPR) correlation favours reporting the SROC
curve over a single pooled point.
"""
from __future__ import annotations
import math
from . import common


def _logit(p):
    return math.log(p / (1 - p))


def _invlogit(x):
    return 1.0 / (1.0 + math.exp(-x))


def _dl_pool(ys, vs):
    w = [1.0 / v for v in vs]
    sw = sum(w)
    mu_fe = sum(wi * y for wi, y in zip(w, ys)) / sw
    Q = sum(wi * (y - mu_fe) ** 2 for wi, y in zip(w, ys))
    df = len(ys) - 1
    sw2 = sum(wi * wi for wi in w)
    tau2 = max(0.0, (Q - df) / (sw - sw2 / sw)) if (sw - sw2 / sw) > 0 else 0.0
    w_re = [1.0 / (v + tau2) for v in vs]
    sw_re = sum(w_re)
    mu = sum(wi * y for wi, y in zip(w_re, ys)) / sw_re
    return {"mu": mu, "se": math.sqrt(1.0 / sw_re), "tau2": tau2}


def _spearman(x, y):
    n = len(x)
    if n < 3:
        return 0.0
    rx = _ranks(x)
    ry = _ranks(y)
    mx = sum(rx) / n
    my = sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else 0.0


def _ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    for rank, i in enumerate(order, 1):
        r[i] = rank
    return r


def analyze(studies):
    """studies: list of {name, tp, fp, fn, tn}. Returns pooled Se/Sp + SROC."""
    logit_se, v_se, logit_fpr, v_fpr = [], [], [], []
    rows = []
    for i, s in enumerate(studies):
        tp, fp, fn, tn = s.get("tp"), s.get("fp"), s.get("fn"), s.get("tn")
        if None in (tp, fp, fn, tn) or (tp + fn) == 0 or (fp + tn) == 0:
            continue
        cc = 0.5 if 0 in (tp, fp, fn, tn) else 0.0
        a, b, c, d = tp + cc, fp + cc, fn + cc, tn + cc
        se = a / (a + c)
        fpr = b / (b + d)
        lse, lfpr = _logit(se), _logit(fpr)
        logit_se.append(lse); v_se.append(1 / a + 1 / c)
        logit_fpr.append(lfpr); v_fpr.append(1 / b + 1 / d)
        rows.append({"name": s.get("name", f"Study {i+1}"),
                     "se": tp / (tp + fn), "sp": tn / (tn + fp),
                     "fpr": fp / (fp + tn)})
    k = len(rows)
    if k < 2:
        return {"available": False, "reason": "need >=2 informative 2x2 tables"}

    pse = _dl_pool(logit_se, v_se)
    pfpr = _dl_pool(logit_fpr, v_fpr)

    # empirical between-study correlation of the logit random effects
    rho = max(-0.95, min(0.95, _correlation(logit_se, logit_fpr)))
    thresh = _spearman(logit_se, logit_fpr)

    mu1, mu2 = pse["mu"], pfpr["mu"]
    se_pool = _invlogit(mu1)
    sp_pool = 1 - _invlogit(mu2)

    # CI on the pooled point (delta method via logit SE)
    def ci_invlogit(mu, se):
        return _invlogit(mu - common.Z975 * se), _invlogit(mu + common.Z975 * se)
    se_lo, se_hi = ci_invlogit(mu1, pse["se"])
    fpr_lo, fpr_hi = ci_invlogit(mu2, pfpr["se"])

    # SROC curve: Moses-Littenberg-style line in logit space using the ratio
    # of the random-effects SDs as the slope (HSROC approximation).
    b = math.sqrt(pse["tau2"] / pfpr["tau2"]) if pfpr["tau2"] > 0 else 1.0
    curve = []
    for i in range(41):
        fpr = 0.01 + 0.98 * i / 40
        lfpr = _logit(fpr)
        lse = mu1 + b * (lfpr - mu2)
        curve.append({"fpr": fpr, "sensitivity": _invlogit(lse), "specificity": 1 - fpr})

    # diagnostic odds ratio
    dor = math.exp(mu1 - mu2)
    return {
        "available": True, "k": k,
        "sensitivity": se_pool, "sensitivity_ci": (se_lo, se_hi),
        "specificity": sp_pool, "specificity_ci": (1 - fpr_hi, 1 - fpr_lo),
        "mu_logit_se": mu1, "mu_logit_fpr": mu2,
        "tau2_se": pse["tau2"], "tau2_fpr": pfpr["tau2"], "rho": rho,
        "dor": dor,
        "threshold_corr": thresh,
        "threshold_effect": thresh < -0.6,
        "sroc": curve,
        "per_study": rows,
        "note": "Bivariate DL approximation (not full REML). Strong negative "
                "threshold correlation favours the SROC curve over a single point.",
    }


def _correlation(x, y):
    n = len(x)
    mx = sum(x) / n
    my = sum(y) / n
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    den = math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
    return num / den if den else 0.0
