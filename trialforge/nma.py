"""metaforge.nma — frequentist contrast-based random-effects network MA.

Model (consistency assumption):
  * Treatments indexed 0..T-1; treatment 0 (or a chosen reference) is the
    baseline. Basic parameters d = (d_1..d_{T-1}) are effects vs reference.
  * Each study contributes contrasts of its non-baseline arms vs its own
    baseline arm. Within-study covariance: Var(contrast_j)=v0+vj,
    Cov(contrast_j,contrast_k)=v0 (shared baseline arm).
  * Random effects: total covariance Sigma_i = S_i + tau^2 * R_i, with
    R_i = I on the diagonal and 1/2 off-diagonal (homogeneous heterogeneity
    variance, shared-baseline correlation 1/2).
  * GLS: d_hat = (sum X_i' Sigma_i^-1 X_i)^-1 (sum X_i' Sigma_i^-1 y_i).
  * tau^2 by the network generalisation of Paule-Mandel: solve
    Q_gen(tau^2) = (#contrasts - #basic params) by bisection (monotone).

Outputs: all pairwise relative effects with 95% CI, SUCRA ranking via
Monte-Carlo over the estimated multivariate-normal of d, and the network's
generalised Q heterogeneity/inconsistency statistic.

Binary arm variance uses log-odds (so the network effect measure is the
log OR). Continuous arms can be supplied as arm means+sd+n (mean-difference
network) via measure="MD".
"""
from __future__ import annotations
import math
import random
from dataclasses import dataclass, field
from typing import Optional
from . import common, linalg


@dataclass
class _StudyContrasts:
    name: str
    base: int            # baseline treatment index
    arms: list           # other treatment indices (in order)
    y: list              # contrast effects vs baseline
    S: list              # within-study covariance matrix (len(arms) x len(arms))


def _arm_logodds_var(e, n):
    a = e + 0.5
    b = n - e + 0.5
    return math.log(a / b), 1 / a + 1 / b


def _arm_mean_var(m, sd, n):
    return m, (sd ** 2) / n


def _build_study(study, tindex, measure):
    """Build contrasts (vs the study's first arm) + within-study covariance."""
    arms = study["arms"]
    if len(arms) < 2:
        return None
    theta, v = [], []
    for arm in arms:
        if measure == "OR":
            th, vv = _arm_logodds_var(arm["e"], arm["n"])
        elif measure == "MD":
            th, vv = _arm_mean_var(arm["m"], arm["sd"], arm["n"])
        else:
            return None
        theta.append(th)
        v.append(vv)
    base_t = tindex[arms[0]["t"]]
    other_t = [tindex[a["t"]] for a in arms[1:]]
    y = [theta[j] - theta[0] for j in range(1, len(arms))]
    m = len(y)
    v0 = v[0]
    S = [[v0 + (v[i + 1] if i == j else 0.0) for j in range(m)] for i in range(m)]
    return _StudyContrasts(study.get("name", "study"), base_t, other_t, y, S), v0


def _design_row(base_t, arm_t, ref, T):
    """Row of the design matrix mapping a contrast (arm vs base) onto basic
    parameters d_1..d_{T-1} (vs ref). +1 at arm, -1 at base, ref column omitted."""
    row = [0.0] * (T - 1)
    cols = [t for t in range(T) if t != ref]
    pos = {t: i for i, t in enumerate(cols)}
    if arm_t != ref:
        row[pos[arm_t]] += 1.0
    if base_t != ref:
        row[pos[base_t]] -= 1.0
    return row


def _assemble(studies_c, ref, T, tau2):
    """Build GLS normal equations LHS (XtSX) and RHS (XtSy) and return also
    the per-study (X_i, Sigma_i^-1, y_i) for residual computation."""
    p = T - 1
    XtSX = [[0.0] * p for _ in range(p)]
    XtSy = [0.0] * p
    blocks = []
    for sc in studies_c:
        m = len(sc.y)
        R = [[1.0 if i == j else 0.5 for j in range(m)] for i in range(m)]
        Sigma = [[sc.S[i][j] + tau2 * R[i][j] for j in range(m)] for i in range(m)]
        Sinv = linalg.inv(Sigma)
        X = [_design_row(sc.base, sc.arms[i], ref, T) for i in range(m)]
        # accumulate
        Xt = linalg.transpose(X)
        XtSinv = linalg.matmul(Xt, Sinv)          # p x m
        XtSX_i = linalg.matmul(XtSinv, X)          # p x p
        XtSy_i = linalg.matvec(XtSinv, sc.y)       # p
        for a in range(p):
            XtSy[a] += XtSy_i[a]
            for b in range(p):
                XtSX[a][b] += XtSX_i[a][b]
        blocks.append((X, Sinv, sc.y))
    return XtSX, XtSy, blocks


def _gen_Q(blocks, d):
    """Generalised Q = sum residual' Sigma^-1 residual."""
    Q = 0.0
    for X, Sinv, y in blocks:
        fitted = linalg.matvec(X, d)
        r = [yi - fi for yi, fi in zip(y, fitted)]
        Sr = linalg.matvec(Sinv, r)
        Q += sum(ri * sri for ri, sri in zip(r, Sr))
    return Q


def _cholesky(A):
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                val = A[i][i] - s
                L[i][j] = math.sqrt(val) if val > 0 else 1e-9
            else:
                L[i][j] = (A[i][j] - s) / L[j][j]
    return L


def analyze(studies, measure="OR", reference=None, tau2_assumed=None,
            smaller_better=True, n_sims=4000, seed=20260528):
    """Run a frequentist random-effects NMA.

    studies: list of {name, arms:[{t, e, n} | {t, m, sd, n}, ...]}
    measure: 'OR' (binary) or 'MD' (continuous)
    """
    # Collect treatments
    treatments = []
    for s in studies:
        for a in s["arms"]:
            if a["t"] not in treatments:
                treatments.append(a["t"])
    treatments.sort()
    T = len(treatments)
    if T < 2:
        return None
    tindex = {t: i for i, t in enumerate(treatments)}
    ref = tindex[reference] if reference in tindex else 0

    studies_c = []
    for s in studies:
        built = _build_study(s, tindex, measure)
        if built is None:
            continue
        studies_c.append(built[0])
    if not studies_c:
        return None

    n_contrasts = sum(len(sc.y) for sc in studies_c)
    p = T - 1
    df = n_contrasts - p

    # --- estimate tau^2 by network Paule-Mandel (bisection) ---------------
    def fit(tau2):
        XtSX, XtSy, blocks = _assemble(studies_c, ref, T, tau2)
        cov = linalg.inv(XtSX)
        d = linalg.matvec(cov, XtSy)
        return d, cov, blocks

    if tau2_assumed is not None:
        tau2 = float(tau2_assumed)
    elif df <= 0:
        tau2 = 0.0
    else:
        d0, _, blocks0 = fit(0.0)
        if _gen_Q(blocks0, d0) <= df:
            tau2 = 0.0
        else:
            lo, hi = 0.0, 1.0
            for _ in range(200):
                d_h, _, blk = fit(hi)
                if _gen_Q(blk, d_h) < df:
                    break
                hi *= 2.0
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                d_m, _, blk = fit(mid)
                q = _gen_Q(blk, d_m)
                if abs(q - df) < 1e-9:
                    break
                if q > df:
                    lo = mid
                else:
                    hi = mid
            tau2 = 0.5 * (lo + hi)

    d, cov, blocks = fit(tau2)
    Q_total = _gen_Q(blocks, d)
    Q_df = max(0, df)
    Q_p = common.chi2_sf(Q_total, Q_df) if Q_df > 0 else float("nan")

    # Full effect vector incl. reference (=0)
    cols = [t for t in range(T) if t != ref]
    full = [0.0] * T
    for i, t in enumerate(cols):
        full[t] = d[i]

    ratio = measure == "OR"

    def disp(v):
        return math.exp(v) if ratio else v

    # Pairwise relative effects (treatment a vs treatment b) with SEs
    def rel_effect(a, b):
        # effect of a vs b = full[a]-full[b]; var via cov of basic params
        ea = 0.0 if a == ref else d[cols.index(a)]
        eb = 0.0 if b == ref else d[cols.index(b)]
        est = ea - eb
        ia = None if a == ref else cols.index(a)
        ib = None if b == ref else cols.index(b)
        var = 0.0
        if ia is not None:
            var += cov[ia][ia]
        if ib is not None:
            var += cov[ib][ib]
        if ia is not None and ib is not None:
            var -= 2 * cov[ia][ib]
        se = math.sqrt(var) if var > 0 else 0.0
        return est, se

    pairwise = []
    for ai in range(T):
        for bi in range(T):
            if ai == bi:
                continue
            est, se = rel_effect(ai, bi)
            pairwise.append({
                "a": treatments[ai], "b": treatments[bi],
                "estimate": disp(est),
                "ci_low": disp(est - common.Z975 * se),
                "ci_high": disp(est + common.Z975 * se),
                "se": se,
            })

    # --- SUCRA via Monte-Carlo over MVN(d, cov) --------------------------
    rng = random.Random(seed)
    L = _cholesky(cov) if p > 0 else [[]]
    rank_sum = [0.0] * T
    better_counts = [[0] * T for _ in range(T)]
    sims = max(500, n_sims)
    for _ in range(sims):
        z = [rng.gauss(0, 1) for _ in range(p)]
        draw = [sum(L[i][j] * z[j] for j in range(i + 1)) for i in range(p)]
        sample = [0.0] * T
        for i, t in enumerate(cols):
            sample[t] = d[i] + draw[i]
        # lower effect = better when smaller_better
        order = sorted(range(T), key=lambda t: sample[t], reverse=not smaller_better)
        for rank, t in enumerate(order):
            # rank 0 = best
            rank_sum[t] += rank
        for a in range(T):
            for b in range(T):
                if a == b:
                    continue
                a_better = (sample[a] < sample[b]) if smaller_better else (sample[a] > sample[b])
                if a_better:
                    better_counts[a][b] += 1
    mean_rank = [rank_sum[t] / sims for t in range(T)]  # 0=best
    sucra = [(T - 1 - mean_rank[t]) / (T - 1) * 100 if T > 1 else 100.0 for t in range(T)]

    ranking = sorted(
        [{"treatment": treatments[t], "sucra": sucra[t], "mean_rank": mean_rank[t] + 1}
         for t in range(T)],
        key=lambda r: -r["sucra"])

    return {
        "measure": measure, "ratio": ratio,
        "treatments": treatments, "reference": treatments[ref],
        "n_studies": len(studies_c), "n_contrasts": n_contrasts,
        "tau2": tau2, "tau": math.sqrt(tau2),
        "Q": Q_total, "Q_df": Q_df, "Q_p": Q_p,
        "smaller_better": smaller_better,
        "pairwise": pairwise,
        "ranking": ranking,
        "vs_reference": [
            {"treatment": treatments[t],
             "estimate": disp(full[t]),
             "ci_low": disp(full[t] - common.Z975 * (math.sqrt(cov[cols.index(t)][cols.index(t)]) if t != ref else 0.0)),
             "ci_high": disp(full[t] + common.Z975 * (math.sqrt(cov[cols.index(t)][cols.index(t)]) if t != ref else 0.0)),
             }
            for t in range(T) if t != ref
        ],
    }
