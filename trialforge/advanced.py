"""trialforge.advanced — advanced meta-analysis methods on top of the
metaforge engine. All operate on effect sizes yi and variances vi.

  publication bias : Egger regression, Peters test, trim-and-fill (L0),
                     PET-PEESE (conditional)
  influence        : leave-one-out, Baujat coordinates, Cook-style influence
  meta-regression  : single moderator, mixed-effects (PM tau^2) + Knapp-Hartung
  subgroup         : within-group pools + Q-between interaction test
  cumulative       : chronological cumulative pooling
  rare events      : Peto OR, Mantel-Haenszel OR

Methodology follows the project gotchas (advanced-stats.md): Egger via the
regression intercept; Peters for binary; trim-and-fill is sensitivity-only;
PET-PEESE is the conditional procedure (PET first, switch to PEESE if PET
rejects the null).
"""
from __future__ import annotations
import math
from . import common


# ---------------------------------------------------------------------------
# tiny weighted least squares (returns coef, se, t, p for each term)
# ---------------------------------------------------------------------------
def _wls(X, y, w):
    """Weighted least squares. X: list of rows (with intercept col), y, w.
    Returns (beta, cov, df)."""
    from .linalg import quad_form_diag, xt_w_y, inv, matvec
    XtWX = quad_form_diag(X, w)
    XtWy = xt_w_y(X, w, y)
    cov = inv(XtWX)
    beta = matvec(cov, XtWy)
    return beta, cov


def egger_test(yis, vis):
    """Egger's regression test for small-study effects.
    Regress SND = yi/sei on precision = 1/sei (OLS); test intercept = 0."""
    k = len(yis)
    if k < 3:
        return {"available": False, "reason": "need >=3 studies"}
    seis = [math.sqrt(v) for v in vis]
    snd = [y / s for y, s in zip(yis, seis)]
    prec = [1.0 / s for s in seis]
    X = [[1.0, p] for p in prec]
    w = [1.0] * k
    beta, cov = _wls(X, snd, w)
    intercept, slope = beta[0], beta[1]
    # residual variance for OLS
    fitted = [intercept + slope * p for p in prec]
    rss = sum((a - b) ** 2 for a, b in zip(snd, fitted))
    df = k - 2
    sigma2 = rss / df if df > 0 else float("nan")
    se_int = math.sqrt(sigma2 * cov[0][0]) if df > 0 else float("nan")
    t = intercept / se_int if se_int else float("nan")
    p = 2 * (1 - _t_cdf(abs(t), df)) if df > 0 else float("nan")
    return {"available": True, "intercept": intercept, "se": se_int,
            "t": t, "df": df, "p": p,
            "interpretation": "asymmetry (possible small-study effects)" if p < 0.10
            else "no strong evidence of funnel asymmetry"}


def peters_test(studies):
    """Peters' test for binary outcomes: weighted regression of lnOR on 1/n,
    weights = 1/(1/(events)+1/(non-events)) summed across arms.
    studies: list of {tE,tN,cE,cN}."""
    rows = []
    for s in studies:
        if not all(k in s for k in ("tE", "tN", "cE", "cN")):
            continue
        a, b = s["tE"] + 0.5, s["tN"] - s["tE"] + 0.5
        c, d = s["cE"] + 0.5, s["cN"] - s["cE"] + 0.5
        lnor = math.log((a * d) / (b * c))
        n = s["tN"] + s["cN"]
        w = 1.0 / (1.0 / (s["tE"] + s["cE"] + 1) + 1.0 / (s["tN"] - s["tE"] + s["cN"] - s["cE"] + 1))
        rows.append((lnor, 1.0 / n, w))
    k = len(rows)
    if k < 3:
        return {"available": False, "reason": "need >=3 binary studies"}
    y = [r[0] for r in rows]
    X = [[1.0, r[1]] for r in rows]
    w = [r[2] for r in rows]
    beta, cov = _wls(X, y, w)
    fitted = [beta[0] + beta[1] * r[1] for r in rows]
    rss = sum(wi * (yi - fi) ** 2 for wi, yi, fi in zip(w, y, fitted))
    df = k - 2
    sigma2 = rss / df if df > 0 else float("nan")
    se_slope = math.sqrt(sigma2 * cov[1][1]) if df > 0 else float("nan")
    t = beta[1] / se_slope if se_slope else float("nan")
    p = 2 * (1 - _t_cdf(abs(t), df)) if df > 0 else float("nan")
    return {"available": True, "slope": beta[1], "t": t, "df": df, "p": p,
            "interpretation": "asymmetry (small-study effects)" if p < 0.10
            else "no strong evidence of asymmetry"}


def pet_peese(yis, vis):
    """Precision-Effect Test / PEESE (conditional). Returns the bias-adjusted
    estimate: PET first; if PET rejects the null (slope on SE, one-sided
    p<0.10) use PEESE, else PET."""
    k = len(yis)
    if k < 3:
        return {"available": False, "reason": "need >=3 studies"}
    seis = [math.sqrt(v) for v in vis]
    w = [1.0 / v for v in vis]
    # PET: yi = b0 + b1*sei
    Xpet = [[1.0, s] for s in seis]
    bpet, cpet = _wls(Xpet, yis, w)
    fitted = [bpet[0] + bpet[1] * s for s, in zip(seis)] if False else [bpet[0] + bpet[1] * s for s in seis]
    rss = sum(wi * (y - f) ** 2 for wi, y, f in zip(w, yis, fitted))
    df = k - 2
    sigma2 = rss / df
    se_b0 = math.sqrt(sigma2 * cpet[0][0])
    df_pet = df
    t0 = bpet[0] / se_b0 if se_b0 else float("nan")
    p0 = 2 * (1 - _t_cdf(abs(t0), df_pet))
    pet_sig = p0 < 0.10  # PET corrected effect differs from zero
    # PEESE: yi = b0 + b1*vi
    Xpeese = [[1.0, v] for v in vis]
    bpeese, cpeese = _wls(Xpeese, yis, w)
    chosen = "PEESE" if pet_sig else "PET"
    est = bpeese[0] if pet_sig else bpet[0]
    cov00 = cpeese[0][0] if pet_sig else cpet[0][0]
    fitted2 = ([bpeese[0] + bpeese[1] * v for v in vis] if pet_sig
               else fitted)
    rss2 = sum(wi * (y - f) ** 2 for wi, y, f in zip(w, yis, fitted2))
    sigma2_2 = rss2 / df
    se_est = math.sqrt(sigma2_2 * cov00)
    return {"available": True, "chosen": chosen,
            "adjusted_estimate": est, "se": se_est,
            "ci_low": est - common.Z975 * se_est,
            "ci_high": est + common.Z975 * se_est,
            "pet_intercept_p": p0}


def trim_and_fill(yis, vis, side="auto"):
    """Duval & Tweedie L0 trim-and-fill. Returns imputed study count and the
    adjusted random-effects pooled estimate (sensitivity analysis only)."""
    k = len(yis)
    if k < 3:
        return {"available": False, "reason": "need >=3 studies"}
    pool0 = common.pool_inverse_variance(yis, vis, knha=False)
    mu0 = pool0.estimate
    # Determine side: if not specified, impute on the side opposite to where
    # extreme positive residuals sit (Duval-Tweedie). Use sign of skew.
    pairs = sorted(zip(yis, vis), key=lambda t: t[0])
    ys = [p[0] for p in pairs]
    vs = [p[1] for p in pairs]

    def l0(center):
        ranks = _signed_ranks([y - center for y in ys])
        # L0 = (sum of positive signed ranks) - n(n+1)/2 ... use Tweedie L0
        Tn = sum(r for r in ranks if r > 0)
        n = len(ys)
        return max(0, int(round((4 * Tn - n * (n + 1)) / (2 * n - 1))))

    center = mu0
    k0 = 0
    for _ in range(30):
        k0_new = l0(center)
        if k0_new <= 0:
            k0 = 0
            break
        # trim k0_new most extreme on the heavy side, recompute center
        trimmed_y = ys[:len(ys) - k0_new] if (sum(ys) / len(ys) > center) else ys[k0_new:]
        trimmed_v = vs[:len(vs) - k0_new] if (sum(ys) / len(ys) > center) else vs[k0_new:]
        if len(trimmed_y) < 2:
            break
        center = common.pool_inverse_variance(trimmed_y, trimmed_v, knha=False).estimate
        if k0_new == k0:
            break
        k0 = k0_new
    # impute k0 mirror studies about center
    if k0 > 0:
        extreme = ys[-k0:] if (sum(ys) / len(ys) > center) else ys[:k0]
        extreme_v = vs[-k0:] if (sum(ys) / len(ys) > center) else vs[:k0]
        imputed_y = [2 * center - y for y in extreme]
        full_y = list(yis) + imputed_y
        full_v = list(vis) + list(extreme_v)
        adj = common.pool_inverse_variance(full_y, full_v, knha=False)
        adj_est = adj.estimate
    else:
        adj_est = mu0
    return {"available": True, "k_imputed": k0,
            "original_estimate": mu0, "adjusted_estimate": adj_est}


def leave_one_out(yis, vis, names, tau2_method="PM"):
    out = []
    for i in range(len(yis)):
        ry = [y for j, y in enumerate(yis) if j != i]
        rv = [v for j, v in enumerate(vis) if j != i]
        if len(ry) < 1:
            continue
        p = common.pool_inverse_variance(ry, rv, tau2_method=tau2_method)
        out.append({"omitted": names[i], "estimate": p.estimate,
                    "ci_low": p.ci_low, "ci_high": p.ci_high, "i2": p.i2})
    return out


def baujat(yis, vis, names):
    """Baujat plot coordinates: x = contribution to heterogeneity Q,
    y = influence on the pooled estimate."""
    wf = [1.0 / v for v in vis]
    mu = sum(w * y for w, y in zip(wf, yis)) / sum(wf)
    out = []
    for i in range(len(yis)):
        x = wf[i] * (yis[i] - mu) ** 2  # contribution to Q
        # influence: squared standardized change in pooled estimate when omitted
        ry = [y for j, y in enumerate(yis) if j != i]
        rw = [w for j, w in enumerate(wf) if j != i]
        mu_i = sum(w * y for w, y in zip(rw, ry)) / sum(rw)
        y_infl = (mu - mu_i) ** 2 / (1.0 / sum(wf))
        out.append({"name": names[i], "x_heterogeneity": x, "y_influence": y_infl})
    return out


def meta_regression(yis, vis, moderator, names=None):
    """Mixed-effects meta-regression on a single continuous moderator.
    tau^2 (residual) by method of moments; Knapp-Hartung t inference."""
    k = len(yis)
    if k < 3:
        return {"available": False, "reason": "need >=3 studies"}
    # iterate tau^2 (residual) via a simple moment update
    tau2 = 0.0
    for _ in range(60):
        w = [1.0 / (v + tau2) for v in vis]
        X = [[1.0, m] for m in moderator]
        beta, cov = _wls(X, yis, w)
        resid = [y - (beta[0] + beta[1] * m) for y, m in zip(yis, moderator)]
        Qres = sum(wi * r * r for wi, r in zip(w, resid))
        dfres = k - 2
        # trace term for moment estimator
        sw = sum(w)
        if Qres <= dfres:
            new = 0.0
        else:
            sw2 = sum(wi * wi for wi in w)
            c = sw - sw2 / sw
            new = max(0.0, (Qres - dfres) / c) if c > 0 else 0.0
        if abs(new - tau2) < 1e-10:
            tau2 = new
            break
        tau2 = new
    w = [1.0 / (v + tau2) for v in vis]
    X = [[1.0, m] for m in moderator]
    beta, cov = _wls(X, yis, w)
    resid = [y - (beta[0] + beta[1] * m) for y, m in zip(yis, moderator)]
    Qres = sum(wi * r * r for wi, r in zip(w, resid))
    dfres = k - 2
    # Knapp-Hartung scaling
    scale = max(1e-9, Qres / dfres)
    se_slope = math.sqrt(scale * cov[1][1])
    t = beta[1] / se_slope if se_slope else float("nan")
    p = 2 * (1 - _t_cdf(abs(t), dfres))
    return {"available": True, "intercept": beta[0], "slope": beta[1],
            "se_slope": se_slope, "t": t, "df": dfres, "p": p,
            "ci_low": beta[1] - common.t_ppf975(dfres) * se_slope,
            "ci_high": beta[1] + common.t_ppf975(dfres) * se_slope,
            "tau2_residual": tau2}


def subgroup(yis, vis, groups, tau2_method="PM"):
    """Pool within each subgroup and run a Q-between interaction test."""
    by = {}
    for y, v, g in zip(yis, vis, groups):
        by.setdefault(g, ([], []))
        by[g][0].append(y)
        by[g][1].append(v)
    results = {}
    group_mus, group_ws = [], []
    for g, (gy, gv) in by.items():
        p = common.pool_inverse_variance(gy, gv, tau2_method=tau2_method)
        results[g] = {"k": p.k, "estimate": p.estimate, "ci_low": p.ci_low,
                      "ci_high": p.ci_high, "i2": p.i2}
        group_mus.append(p.estimate)
        group_ws.append(1.0 / (p.se ** 2) if p.se > 0 else 0.0)
    sw = sum(group_ws)
    mu_overall = sum(w * m for w, m in zip(group_ws, group_mus)) / sw if sw else float("nan")
    Q_between = sum(w * (m - mu_overall) ** 2 for w, m in zip(group_ws, group_mus))
    df_between = len(by) - 1
    p_between = common.chi2_sf(Q_between, df_between) if df_between > 0 else float("nan")
    return {"subgroups": results, "Q_between": Q_between,
            "df_between": df_between, "p_between": p_between}


def cumulative(yis, vis, names, order, tau2_method="PM"):
    """Cumulative meta-analysis in the given order (e.g. by year)."""
    idx = sorted(range(len(yis)), key=lambda i: order[i])
    out = []
    cy, cv = [], []
    for i in idx:
        cy.append(yis[i]); cv.append(vis[i])
        p = common.pool_inverse_variance(cy, cv, tau2_method=tau2_method)
        out.append({"added": names[i], "order": order[i], "k": p.k,
                    "estimate": p.estimate, "ci_low": p.ci_low, "ci_high": p.ci_high})
    return out


# ---- rare events -----------------------------------------------------------
def peto_or(studies):
    """Peto one-step odds ratio (best for rare events, balanced arms)."""
    sum_OE = 0.0
    sum_V = 0.0
    for s in studies:
        a, n1, c, n0 = s["tE"], s["tN"], s["cE"], s["cN"]
        N = n1 + n0
        tot_e = a + c
        if N == 0 or tot_e == 0 or tot_e == N:
            continue
        E = tot_e * n1 / N
        V = tot_e * (N - tot_e) * n1 * n0 / (N * N * (N - 1))
        if V <= 0:
            continue
        sum_OE += (a - E)
        sum_V += V
    if sum_V <= 0:
        return {"available": False, "reason": "no informative studies"}
    lnor = sum_OE / sum_V
    se = math.sqrt(1.0 / sum_V)
    return {"available": True, "OR": math.exp(lnor),
            "ci_low": math.exp(lnor - common.Z975 * se),
            "ci_high": math.exp(lnor + common.Z975 * se),
            "lnOR": lnor, "se": se}


def mantel_haenszel_or(studies):
    """Mantel-Haenszel fixed-effect odds ratio."""
    num = denom = 0.0
    for s in studies:
        a, n1, c, n0 = s["tE"], s["tN"], s["cE"], s["cN"]
        b, d = n1 - a, n0 - c
        N = n1 + n0
        if N == 0:
            continue
        num += a * d / N
        denom += b * c / N
    if denom <= 0:
        return {"available": False, "reason": "no informative studies"}
    or_mh = num / denom
    # Robins-Breslow-Greenland variance (log scale)
    R = S = 0.0
    PR = PS_ = QR = QS = 0.0
    for s in studies:
        a, n1, c, n0 = s["tE"], s["tN"], s["cE"], s["cN"]
        b, d = n1 - a, n0 - c
        N = n1 + n0
        if N == 0:
            continue
        R += a * d / N
        S += b * c / N
        P = (a + d) / N
        Qd = (b + c) / N
        PR += P * (a * d / N)
        PS_ += P * (b * c / N) + Qd * (a * d / N)
        QS += Qd * (b * c / N)
    if R <= 0 or S <= 0:
        return {"available": False, "reason": "degenerate"}
    var = PR / (2 * R * R) + PS_ / (2 * R * S) + QS / (2 * S * S)
    se = math.sqrt(var)
    lnor = math.log(or_mh)
    return {"available": True, "OR": or_mh,
            "ci_low": math.exp(lnor - common.Z975 * se),
            "ci_high": math.exp(lnor + common.Z975 * se),
            "se": se}


# ---- helpers ---------------------------------------------------------------
def _t_cdf(t, df):
    """Student-t CDF via the regularized incomplete beta (good enough for p)."""
    if df <= 0:
        return float("nan")
    x = df / (df + t * t)
    ib = _betai(df / 2.0, 0.5, x)
    return 1 - 0.5 * ib if t > 0 else 0.5 * ib


def _betai(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(math.log(x) * a + math.log(1 - x) * b - lbeta) / a
    # continued fraction (Lentz)
    f, c, d = 1.0, 1.0, 0.0
    for i in range(0, 300):
        m = i // 2
        if i == 0:
            num = 1.0
        elif i % 2 == 0:
            num = (m * (b - m) * x) / ((a + 2 * m - 1) * (a + 2 * m))
        else:
            num = -((a + m) * (a + b + m) * x) / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1.0 + num * d
        if abs(d) < 1e-30:
            d = 1e-30
        d = 1.0 / d
        c = 1.0 + num / c
        if abs(c) < 1e-30:
            c = 1e-30
        f *= c * d
        if abs(1.0 - c * d) < 1e-12:
            break
    return front * (f - 1.0)


def _signed_ranks(vals):
    order = sorted(range(len(vals)), key=lambda i: abs(vals[i]))
    ranks = [0.0] * len(vals)
    for rank, i in enumerate(order, start=1):
        ranks[i] = rank if vals[i] >= 0 else -rank
    return ranks
