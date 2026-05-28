"""trialforge.evalue — E-value for unmeasured confounding (VanderWeele & Ding 2017).

The E-value is the minimum strength of association (on the risk-ratio
scale) that an unmeasured confounder would need with BOTH the exposure and
the outcome to fully explain away an observed association. Larger E-value =
more robust to confounding. Useful for pooled estimates from observational
studies.

For a risk ratio RR (>1; if RR<1 use 1/RR first):
    E = RR + sqrt(RR * (RR - 1))
The E-value for the CI limit closest to the null is computed the same way
(=1 if the CI already crosses the null).

Approximate conversions to RR (VanderWeele 2017) are applied when the
pooled estimate is an OR (rare-outcome: RR~=OR; common-outcome: a sqrt
approximation) or HR.
"""
from __future__ import annotations
import math


def _e(rr):
    if rr < 1:
        rr = 1.0 / rr
    if rr <= 1:
        return 1.0
    return rr + math.sqrt(rr * (rr - 1))


def _to_rr(estimate, measure, rare=True):
    m = measure.upper()
    if m == "RR":
        return estimate
    if m == "OR":
        if rare:
            return estimate          # rare outcome: OR ~= RR
        return math.sqrt(estimate)   # common-outcome sqrt(OR) approximation
    if m == "HR":
        if rare or estimate == 1:
            return estimate          # rare outcome: HR ~= RR
        # common-outcome HR->RR (VanderWeele-Ding 2017):
        return (1 - 0.5 ** math.sqrt(estimate)) / (1 - 0.5 ** math.sqrt(1.0 / estimate))
    return estimate


def analyze(point, ci_low, ci_high, measure="RR", rare_outcome=True):
    """Compute the E-value for a pooled point estimate and the CI limit
    nearest the null. `point`, `ci_low`, `ci_high` are on the natural
    (ratio) scale."""
    rr = _to_rr(point, measure, rare_outcome)
    # CI limit closest to the null (RR=1)
    if ci_low <= 1 <= ci_high:
        e_ci = 1.0
        limit_used = 1.0
    else:
        # both bounds same side; nearest-to-null bound
        if rr > 1:
            limit_used = ci_low
        else:
            limit_used = ci_high
        e_ci = _e(_to_rr(limit_used, measure, rare_outcome))
    return {
        "available": True,
        "measure": measure,
        "rr_used": rr,
        "evalue_point": _e(rr),
        "evalue_ci": e_ci,
        "ci_crosses_null": ci_low <= 1 <= ci_high,
        "interpretation": (
            f"An unmeasured confounder would need a risk-ratio association of "
            f">= {_e(rr):.2f} with both exposure and outcome to explain away the "
            f"point estimate"
            + (", but the confidence interval already includes the null."
               if ci_low <= 1 <= ci_high else
               f"; >= {e_ci:.2f} to shift the CI to the null.")),
        "note": "VanderWeele-Ding E-value; for OR/HR an approximate RR "
                "conversion is applied (set rare_outcome=false for common outcomes).",
    }
