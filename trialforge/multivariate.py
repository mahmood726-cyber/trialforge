"""trialforge.multivariate — bivariate (two-outcome) meta-analysis with
borrowing of strength (Riley 2007).

Jointly pools two correlated outcomes. Each study contributes (y1, y2) with
a known within-study covariance built from (se1, se2, within-study rho_w),
plus an estimated between-study covariance (tau1^2, tau2^2, rho_b). Because
the outcomes are correlated, the pooled estimate of each outcome BORROWS
information from the other — tightening its CI, especially when some studies
report only one outcome.

Fit by REML over (tau1^2, tau2^2, rho_b) (Nelder-Mead); given the variance
components the pooled means are GLS closed-form. Reports each outcome's
univariate vs bivariate SE so the strength borrowed is visible.

Studies missing one outcome (y2 = null) are allowed — that is exactly where
borrowing helps most.
"""
from __future__ import annotations
import math
from . import common, linalg


def _inv2(M):
    a, b, c, d = M[0][0], M[0][1], M[1][0], M[1][1]
    det = a * d - b * c
    return [[d / det, -b / det], [-c / det, a / det]], det


def _gls(studies, t1, t2, rb):
    """Return pooled (mu1, mu2), covariance 2x2, and REML loglik."""
    between = [[t1 * t1, rb * t1 * t2], [rb * t1 * t2, t2 * t2]]
    A = [[0.0, 0.0], [0.0, 0.0]]   # sum Sigma_i^-1
    b = [0.0, 0.0]                 # sum Sigma_i^-1 y_i
    logdet = 0.0
    sum_yty = 0.0                  # sum y_i' Sigma_i^-1 y_i
    for s in studies:
        has2 = s.get("y2") is not None
        if has2:
            se1, se2 = s["se1"], s["se2"]
            rw = s.get("rho_w", 0.0)
            W = [[se1 * se1, rw * se1 * se2], [rw * se1 * se2, se2 * se2]]
            Sig = [[W[0][0] + between[0][0], W[0][1] + between[0][1]],
                   [W[1][0] + between[1][0], W[1][1] + between[1][1]]]
            Sinv, det = _inv2(Sig)
            y = [s["y1"], s["y2"]]
        else:
            # outcome 1 only -> 1x1 problem embedded
            v = s["se1"] ** 2 + between[0][0]
            Sinv = [[1.0 / v, 0.0], [0.0, 0.0]]
            det = v
            y = [s["y1"], 0.0]
        for i in range(2):
            for j in range(2):
                A[i][j] += Sinv[i][j]
            b[i] += sum(Sinv[i][k] * y[k] for k in range(2))
        logdet += math.log(abs(det))
        sum_yty += sum(y[i] * sum(Sinv[i][j] * y[j] for j in range(2)) for i in range(2))

    Ainv, detA = _inv2(A) if A[1][1] != 0 else (None, None)
    if Ainv is None:
        mu1 = b[0] / A[0][0]
        # REML with residual quad = sum_yty - mu1*b0 ; |A| = A00
        quad = sum_yty - mu1 * b[0]
        ll = -0.5 * logdet - 0.5 * quad - 0.5 * math.log(abs(A[0][0]))
        return (mu1, None), [[1.0 / A[0][0], 0], [0, 0]], ll
    mu = [Ainv[0][0] * b[0] + Ainv[0][1] * b[1],
          Ainv[1][0] * b[0] + Ainv[1][1] * b[1]]
    # residual quadratic form = sum_yty - mu' b  (since mu = A^-1 b)
    quad = sum_yty - (mu[0] * b[0] + mu[1] * b[1])
    ll = -0.5 * logdet - 0.5 * quad - 0.5 * math.log(abs(detA))
    return (mu[0], mu[1]), Ainv, ll


def analyze(studies, *, ratio=False):
    """studies: [{name, y1, se1, y2, se2, rho_w}] (y2/se2 may be null)."""
    n2 = sum(1 for s in studies if s.get("y2") is not None)
    if len(studies) < 3 or n2 < 2:
        return {"available": False, "reason": "need >=3 studies and >=2 with both outcomes"}

    # REML optimisation over (t1, t2, rb) via Nelder-Mead on (log t1, log t2, atanh rb)
    def negll(p):
        t1 = math.exp(p[0]); t2 = math.exp(p[1]); rb = math.tanh(p[2])
        rb = max(-0.95, min(0.95, rb))
        try:
            _, _, ll = _gls(studies, t1, t2, rb)
        except ZeroDivisionError:
            return 1e18
        return -ll

    p = _nelder_mead(negll, [math.log(0.1), math.log(0.1), 0.0])
    t1 = math.exp(p[0]); t2 = math.exp(p[1]); rb = max(-0.95, min(0.95, math.tanh(p[2])))
    (mu1, mu2), cov, _ = _gls(studies, t1, t2, rb)

    # univariate SE for outcome 1 (no borrowing) for comparison
    y1s = [s["y1"] for s in studies]
    v1s = [s["se1"] ** 2 for s in studies]
    uni1 = common.pool_inverse_variance(y1s, v1s, tau2_method="DL")
    s2_studies = [s for s in studies if s.get("y2") is not None]
    uni2 = common.pool_inverse_variance([s["y2"] for s in s2_studies],
                                        [s["se2"] ** 2 for s in s2_studies],
                                        tau2_method="DL")

    def disp(v):
        return math.exp(v) if (ratio and v is not None) else v

    se1_bi = math.sqrt(cov[0][0])
    se2_bi = math.sqrt(cov[1][1]) if mu2 is not None else None
    z = common.Z975
    return {
        "available": True, "k": len(studies), "k_both": n2,
        "tau1": t1, "tau2": t2, "rho_between": rb,
        "outcome1": {
            "bivariate": disp(mu1),
            "ci_low": disp(mu1 - z * se1_bi), "ci_high": disp(mu1 + z * se1_bi),
            "se_bivariate": se1_bi, "se_univariate": uni1.se,
            "borrowed_precision_pct": 100 * (1 - se1_bi / uni1.se) if uni1.se else 0.0,
        },
        "outcome2": ({
            "bivariate": disp(mu2),
            "ci_low": disp(mu2 - z * se2_bi), "ci_high": disp(mu2 + z * se2_bi),
            "se_bivariate": se2_bi, "se_univariate": uni2.se,
            "borrowed_precision_pct": 100 * (1 - se2_bi / uni2.se) if uni2.se else 0.0,
        } if mu2 is not None else None),
        "note": "Bivariate REML (Riley 2007); each outcome borrows strength "
                "from the other via the between-study correlation.",
    }


def _nelder_mead(f, x0, iters=300, step=0.5):
    n = len(x0)
    simplex = [list(x0)]
    for i in range(n):
        x = list(x0); x[i] += step
        simplex.append(x)
    fv = [f(x) for x in simplex]
    for _ in range(iters):
        order = sorted(range(n + 1), key=lambda i: fv[i])
        simplex = [simplex[i] for i in order]; fv = [fv[i] for i in order]
        cent = [sum(simplex[i][j] for i in range(n)) / n for j in range(n)]
        xr = [cent[j] + (cent[j] - simplex[-1][j]) for j in range(n)]
        fr = f(xr)
        if fv[0] <= fr < fv[-2]:
            simplex[-1], fv[-1] = xr, fr
        elif fr < fv[0]:
            xe = [cent[j] + 2 * (cent[j] - simplex[-1][j]) for j in range(n)]
            fe = f(xe)
            simplex[-1], fv[-1] = (xe, fe) if fe < fr else (xr, fr)
        else:
            xc = [cent[j] + 0.5 * (simplex[-1][j] - cent[j]) for j in range(n)]
            fc = f(xc)
            if fc < fv[-1]:
                simplex[-1], fv[-1] = xc, fc
            else:
                for i in range(1, n + 1):
                    simplex[i] = [simplex[0][j] + 0.5 * (simplex[i][j] - simplex[0][j])
                                  for j in range(n)]
                    fv[i] = f(simplex[i])
        if abs(fv[-1] - fv[0]) < 1e-9:
            break
    return simplex[0]
