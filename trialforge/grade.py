"""trialforge.grade — GRADE certainty rating + Summary-of-Findings row.

Automates the data-driven GRADE domains for a pairwise meta-analysis and
surfaces the judgment domains for the reviewer:

  * risk of bias    (judgment — user supplies)
  * inconsistency   (DATA: I^2 — >50% one level, >75% two)
  * indirectness    (judgment — user supplies)
  * imprecision     (DATA: CI crosses null, or optimal information size
                     not met) -> downgrade
  * publication bias(DATA: Egger p<0.10 with k>=10 -> suspected)

Starts at HIGH (RCT evidence) and downgrades. Also builds the SoF absolute
effect: given a comparator (baseline) risk, converts the pooled relative
effect to events per 1000 with a 95% CI.
"""
from __future__ import annotations
import math
from . import common, advanced

_LEVELS = ["high", "moderate", "low", "very low"]


def _level(steps):
    return _LEVELS[min(len(_LEVELS) - 1, max(0, steps))]


def _abs_per_1000(baseline_risk, rr_point, rr_lo, rr_hi):
    """Absolute risk per 1000 from a baseline risk and a risk ratio."""
    def risk(rr):
        return min(1.0, max(0.0, baseline_risk * rr))
    base = baseline_risk * 1000
    return {
        "comparator_per_1000": round(base),
        "intervention_per_1000": round(risk(rr_point) * 1000),
        "difference_per_1000": round((risk(rr_point) - baseline_risk) * 1000),
        "diff_ci_per_1000": (round((risk(rr_lo) - baseline_risk) * 1000),
                             round((risk(rr_hi) - baseline_risk) * 1000)),
    }


def rate(*, measure, estimate, ci_low, ci_high, k, n_total, i2,
         baseline_risk=None, egger_p=None,
         risk_of_bias="not assessed", indirectness="not assessed",
         design="randomised trials"):
    ratio = measure.upper() in ("OR", "RR", "HR")
    null = 1.0 if ratio else 0.0
    domains = {}
    steps = 0

    # inconsistency (data)
    if i2 is None:
        domains["inconsistency"] = "not estimable"
    elif i2 > 75:
        domains["inconsistency"] = "serious (-2)"; steps += 2
    elif i2 > 50:
        domains["inconsistency"] = "serious (-1)"; steps += 1
    else:
        domains["inconsistency"] = "not serious"

    # imprecision (data): CI crossing null OR few events/participants
    crosses = ci_low <= null <= ci_high
    small = (n_total is not None and n_total < 400)
    if crosses and small:
        domains["imprecision"] = "serious (-2)"; steps += 2
    elif crosses or small:
        domains["imprecision"] = "serious (-1)"; steps += 1
    else:
        domains["imprecision"] = "not serious"

    # publication bias (data)
    if egger_p is not None and k >= 10 and egger_p < 0.10:
        domains["publication_bias"] = "suspected (-1)"; steps += 1
    elif k < 10:
        domains["publication_bias"] = "undetected (too few studies to assess)"
    else:
        domains["publication_bias"] = "undetected"

    # judgment domains
    for name, val in (("risk_of_bias", risk_of_bias), ("indirectness", indirectness)):
        domains[name] = val
        if isinstance(val, str):
            if "serious (-2)" in val or val == "major concerns":
                steps += 2
            elif "serious (-1)" in val or val == "some concerns" or val == "serious":
                steps += 1

    certainty = _level(steps)
    sof = None
    if ratio and baseline_risk is not None and measure.upper() == "RR":
        sof = _abs_per_1000(baseline_risk, estimate, ci_low, ci_high)
    elif ratio and baseline_risk is not None and measure.upper() == "OR":
        # convert OR -> RR at the baseline risk (rare-outcome friendly)
        def or_to_rr(orv):
            return orv / (1 - baseline_risk + baseline_risk * orv)
        sof = _abs_per_1000(baseline_risk, or_to_rr(estimate),
                            or_to_rr(ci_low), or_to_rr(ci_high))

    return {
        "available": True,
        "design": design, "certainty": certainty, "downgrade_steps": steps,
        "domains": domains,
        "relative_effect": {"measure": measure, "estimate": estimate,
                            "ci_low": ci_low, "ci_high": ci_high},
        "n_studies": k, "n_participants": n_total,
        "absolute_effect": sof,
        "note": "Data-driven domains (inconsistency, imprecision, publication "
                "bias) are automated; risk of bias and indirectness need "
                "reviewer judgment. Observational designs should start at 'low'.",
    }
