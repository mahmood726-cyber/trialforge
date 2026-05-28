"""metaforge.doseresponse — two-stage dose-response meta-analysis.

Stage 1 (per study): fit a linear trend of the (log) effect on dose
through the study's reference dose, by weighted least squares over the
dose-specific effect estimates. Returns a per-study slope (effect per
unit dose) and its variance.

Stage 2: pool the per-study slopes with a random-effects model
(Paule-Mandel tau^2 + Knapp-Hartung), then predict the dose-response
curve effect(dose) = slope * (dose - reference_dose).

Input per study, either:
  * binary arm data:  doses: [{dose, e, n}, ...]  (first/lowest = reference)
  * precomputed:      doses: [{dose, effect, ci_low, ci_high}, ...]
                      (effect of that dose vs the reference dose; the
                       reference row may be omitted or given as effect=1)

Note: this is the *approximate* two-stage trend (within-study correlation
between dose levels is not reconstructed — that needs cases + person-time
per level via the Greenland-Longnecker covariance). It is the standard,
robust first-line dose-response summary; the slope and pooled trend are
unbiased, the per-study SE is mildly conservative.
"""
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Optional
from . import common


@dataclass
class StudySlope:
    name: str
    slope: float
    var: float
    points: list  # [(dose, y, v)]


def _study_points(study, measure):
    doses = study.get("doses", [])
    if len(doses) < 2:
        return None
    # reference = lowest dose
    doses = sorted(doses, key=lambda d: d["dose"])
    ref_dose = doses[0]["dose"]
    pts = []
    if all("e" in d and "n" in d for d in doses):
        e0, n0 = doses[0]["e"], doses[0]["n"]
        a0, b0 = e0 + 0.5, n0 - e0 + 0.5
        for d in doses[1:]:
            e, n = d["e"], d["n"]
            a, b = e + 0.5, n - e + 0.5
            if measure == "OR":
                y = math.log((a * (n0 - e0 + 0.5)) / (b * (e0 + 0.5)))
                v = 1/a + 1/b + 1/a0 + 1/b0
            else:  # RR
                y = math.log((a / (n + 1)) / (a0 / (n0 + 1)))
                v = 1/a - 1/(n + 1) + 1/a0 - 1/(n0 + 1)
            pts.append((d["dose"] - ref_dose, y, max(v, 1e-9)))
    elif all("effect" in d and "ci_low" in d and "ci_high" in d for d in doses if d["dose"] != ref_dose):
        for d in doses:
            if d["dose"] == ref_dose:
                continue
            eff = d["effect"]
            if eff <= 0 or d["ci_low"] <= 0 or d["ci_high"] <= 0:
                continue
            y = math.log(eff)
            se = (math.log(d["ci_high"]) - math.log(d["ci_low"])) / (2 * common.Z975)
            pts.append((d["dose"] - ref_dose, y, max(se * se, 1e-9)))
    else:
        return None
    return pts


def _fit_slope(pts):
    """Weighted least squares slope through the origin (reference dose)."""
    sxx = sum((x ** 2) / v for x, y, v in pts)
    sxy = sum((x * y) / v for x, y, v in pts)
    if sxx <= 0:
        return None
    slope = sxy / sxx
    var = 1.0 / sxx
    return slope, var


def analyze(studies, measure="RR", tau2_method="PM", predict_doses=None):
    slopes = []
    skipped = []
    for i, s in enumerate(studies):
        name = s.get("name") or f"Study {i+1}"
        pts = _study_points(s, measure)
        if not pts:
            skipped.append(name)
            continue
        fit = _fit_slope(pts)
        if fit is None:
            skipped.append(name)
            continue
        slopes.append(StudySlope(name, fit[0], fit[1], pts))
    if not slopes:
        return None, slopes, skipped

    pool = common.pool_inverse_variance(
        [s.slope for s in slopes], [s.var for s in slopes],
        tau2_method=tau2_method, knha=True)

    ratio = measure in ("OR", "RR")
    # Predicted dose-response curve relative to reference dose (=0 increment)
    if predict_doses is None:
        # default grid: 0 .. max observed increment, 25 points
        max_inc = max((x for s in slopes for x, _, _ in s.points), default=1.0)
        predict_doses = [max_inc * i / 24 for i in range(25)]
    slope = pool.estimate
    se = pool.se
    crit = common.Z975
    curve = []
    for dose_inc in predict_doses:
        eff = slope * dose_inc
        lo = (slope - crit * se) * dose_inc if dose_inc >= 0 else (slope + crit * se) * dose_inc
        hi = (slope + crit * se) * dose_inc if dose_inc >= 0 else (slope - crit * se) * dose_inc
        curve.append({
            "dose_increment": dose_inc,
            "effect": math.exp(eff) if ratio else eff,
            "ci_low": math.exp(lo) if ratio else lo,
            "ci_high": math.exp(hi) if ratio else hi,
        })

    pool.extra = {
        "measure": measure, "ratio": ratio,
        "slope_per_unit": slope,
        "slope_ci": (pool.ci_low, pool.ci_high),
        "slope_display_per_unit": math.exp(slope) if ratio else slope,
        "curve": curve,
        "per_study": [
            {"name": s.name, "slope": s.slope,
             "slope_display": math.exp(s.slope) if ratio else s.slope,
             "weight": pool.weights[i], "n_levels": len(s.points) + 1}
            for i, s in enumerate(slopes)
        ],
    }
    return pool, slopes, skipped
