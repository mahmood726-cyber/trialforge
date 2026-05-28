"""trialforge.glmm — rare-events binomial/hypergeometric-normal GLMM.

Pure-stdlib re-implementation of the allmeta `rare-events-glmm` engine
(Stijnen 2010). Models sparse 2x2 data with the conditional
(non-central hypergeometric) likelihood per study and a normal
random-effect on the log odds ratio, integrated by 10-point Gauss-Hermite
quadrature. Avoids the +0.5 continuity correction that biases the OR
toward the null when events are rare, and handles zero cells natively.

theta_i ~ N(theta, tau^2);  a_i | m1_i ~ NoncentralHypergeometric(psi=e^{theta_i})
Maximise sum_i log integral over theta_i  w.r.t. (theta, tau^2).

Implements the same conditional formulation as
metafor::rma.glmm(measure="OR", model="CM.AL", method="ML"). (No R-parity
fixture is bundled; the test checks agreement with the Peto OR on rare
balanced data rather than an exact metafor reference.)
"""
from __future__ import annotations
import math
from . import common

# 10-point Gauss-Hermite nodes/weights (for integral f(x) e^{-x^2} dx)
_GH_X = [-3.4361591188377376, -2.5327316742327897, -1.7566836492998817,
         -1.0366108297895136, -0.34290132722370461, 0.34290132722370461,
         1.0366108297895136, 1.7566836492998817, 2.5327316742327897,
         3.4361591188377376]
_GH_W = [7.6404328552326206e-6, 0.0013436457467812327, 0.033874394455481063,
         0.24013861108231469, 0.6108626337353258, 0.6108626337353258,
         0.24013861108231469, 0.033874394455481063, 0.0013436457467812327,
         7.6404328552326206e-6]
_INV_SQRT_PI = 1.0 / math.sqrt(math.pi)


def _logchoose(n, k):
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _log_nchg_pmf(a, n1, n0, m1, log_psi):
    """log P(X=a) for Fisher non-central hypergeometric with log odds ratio
    log_psi, drawing m1 'events' split across arms of size n1 and n0."""
    lo = max(0, m1 - n0)
    hi = min(n1, m1)
    if a < lo or a > hi:
        return -math.inf
    # log numerator terms over support for normalisation (log-sum-exp)
    terms = []
    for j in range(lo, hi + 1):
        terms.append(_logchoose(n1, j) + _logchoose(n0, m1 - j) + j * log_psi)
    mx = max(terms)
    log_denom = mx + math.log(sum(math.exp(t - mx) for t in terms))
    log_num = _logchoose(n1, a) + _logchoose(n0, m1 - a) + a * log_psi
    return log_num - log_denom


def _study_loglik(study, theta, tau):
    a = study["tE"]; n1 = study["tN"]; c = study["cE"]; n0 = study["cN"]
    m1 = a + c
    if m1 == 0 or m1 == n1 + n0:
        return 0.0  # non-informative study (conditional likelihood is degenerate)
    if tau <= 0:
        return _log_nchg_pmf(a, n1, n0, m1, theta)
    # integrate over theta_i = theta + sqrt(2)*tau*x
    vals = []
    for x, w in zip(_GH_X, _GH_W):
        lp = _log_nchg_pmf(a, n1, n0, m1, theta + math.sqrt(2.0) * tau * x)
        vals.append(math.log(w) + lp)
    mx = max(vals)
    if mx == -math.inf:
        return -math.inf
    return math.log(_INV_SQRT_PI) + mx + math.log(sum(math.exp(v - mx) for v in vals))


def _total_loglik(studies, theta, tau):
    return sum(_study_loglik(s, theta, tau) for s in studies)


def analyze(studies):
    """studies: list of {name, tE, tN, cE, cN}. Returns pooled OR + tau^2."""
    informative = [s for s in studies
                   if (s["tE"] + s["cE"]) not in (0, s["tN"] + s["cN"])]
    if len(informative) < 2:
        return {"available": False, "reason": "need >=2 informative studies "
                "(at least one event and one non-event across arms)"}

    # 2-D maximisation over (theta, tau) by coordinate ascent + golden section.
    def f(theta, tau):
        return _total_loglik(informative, theta, max(0.0, tau))

    theta, tau = 0.0, 0.1

    def golden(g, a, b):
        gr = (math.sqrt(5) - 1) / 2
        c = b - gr * (b - a); d = a + gr * (b - a)
        fc = g(c) ; fd = g(d)
        for _ in range(60):
            if fc < fd:
                a, c, fc = c, d, fd; d = a + gr * (b - a); fd = g(d)
            else:
                b, d, fd = d, c, fc; c = b - gr * (b - a); fc = g(c)
            if abs(b - a) < 1e-7:
                break
        return 0.5 * (a + b)

    for _ in range(40):
        new_theta = golden(lambda t: f(t, tau), -5.0, 5.0)
        new_tau = golden(lambda u: f(new_theta, u), 0.0, 3.0)
        if abs(new_theta - theta) < 1e-6 and abs(new_tau - tau) < 1e-6:
            theta, tau = new_theta, new_tau
            break
        theta, tau = new_theta, new_tau

    # SE of theta from the curvature of the PROFILE log-likelihood: at each
    # perturbed theta, tau is re-optimised (so the SE accounts for the
    # theta-tau covariance and is not anti-conservative).
    def profile(th):
        return f(th, golden(lambda u: f(th, u), 0.0, 3.0))
    h = 1e-2
    ll0 = profile(theta)
    llp = profile(theta + h)
    llm = profile(theta - h)
    curv = (llp - 2 * ll0 + llm) / (h * h)
    se = math.sqrt(-1.0 / curv) if curv < 0 else float("nan")

    return {
        "available": True, "k": len(informative),
        "logOR": theta, "OR": math.exp(theta),
        "se": se,
        "ci_low": math.exp(theta - common.Z975 * se) if se == se else None,
        "ci_high": math.exp(theta + common.Z975 * se) if se == se else None,
        "tau2": tau * tau, "tau": tau,
        "note": "Conditional hypergeometric-normal GLMM (Stijnen 2010); "
                "no continuity correction; same conditional formulation as "
                "metafor::rma.glmm CM.AL (not exactly parity-tested).",
    }
