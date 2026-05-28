"""trialforge.survival — RMST (restricted mean survival time) meta-analysis.

Two entry points:
  * rmst_from_km(points, tau): area under a Kaplan-Meier step curve up to
    tau* (points = [(time, survival), ...], survival starting at 1.0).
  * analyze(studies, tau): random-effects pool of per-study RMST DIFFERENCES
    (intervention - comparator). Differences are pooled, never ratios
    (advanced-stats.md). tau* must be stated and common across studies.

Each study supplies either a precomputed rmst_diff + se (or ci), or two
arms' KM step coordinates (km_t / km_c) from which RMST and the difference
are computed; in the coordinate case a standard error must still be given
(reconstructing IPD-level SE from a published curve needs Guyot digitisation
+ numbers-at-risk and is out of scope — we never claim IPD accuracy).
"""
from __future__ import annotations
import math
from . import common


def rmst_from_km(points, tau):
    """Area under the KM step function up to tau. points sorted by time,
    survival is the step value AFTER each event time; S(0)=1."""
    pts = sorted(points, key=lambda p: p[0])
    area = 0.0
    prev_t = 0.0
    prev_s = 1.0
    for t, s in pts:
        t_clip = min(t, tau)
        if t_clip > prev_t:
            area += prev_s * (t_clip - prev_t)
            prev_t = t_clip
        prev_s = s
        if prev_t >= tau:
            break
    if prev_t < tau:
        area += prev_s * (tau - prev_t)
    return area


def analyze(studies, tau=None, tau2_method="PM"):
    """Pool RMST differences. Returns pooled difference (time units) + CI/PI."""
    yis, vis, rows = [], [], []
    skipped = []
    for i, s in enumerate(studies):
        name = s.get("name", f"Study {i+1}")
        if "rmst_diff" in s and ("se" in s or ("ci_low" in s and "ci_high" in s)):
            diff = s["rmst_diff"]
            if "se" in s:
                se = s["se"]
            else:
                se = (s["ci_high"] - s["ci_low"]) / (2 * common.Z975)
        elif "km_t" in s and "km_c" in s and tau is not None and "se" in s:
            r_t = rmst_from_km(s["km_t"], tau)
            r_c = rmst_from_km(s["km_c"], tau)
            diff = r_t - r_c
            se = s["se"]
        else:
            skipped.append(name)
            continue
        if se <= 0:
            skipped.append(name)
            continue
        yis.append(diff); vis.append(se * se)
        rows.append({"name": name, "diff": diff,
                     "lo": diff - common.Z975 * se, "hi": diff + common.Z975 * se})
    if len(yis) < 1:
        return {"available": False, "reason": "no usable RMST data", "skipped": skipped}
    pool = common.pool_inverse_variance(yis, vis, tau2_method=tau2_method)
    for i, r in enumerate(rows):
        r["weight"] = pool.weights[i]
    return {
        "available": True, "k": pool.k, "tau_star": tau,
        "rmst_difference": pool.estimate,
        "ci_low": pool.ci_low, "ci_high": pool.ci_high,
        "pi_low": pool.pi_low, "pi_high": pool.pi_high,
        "i2": pool.i2, "tau2": pool.tau2,
        "per_study": rows, "skipped": skipped,
        "note": "RMST differences pooled (not ratios); tau* must be common. "
                "KM-reconstructed inputs do not carry IPD-level precision.",
    }
