"""trialforge.tsa — Trial Sequential Analysis (information-based).

Treats a cumulative meta-analysis like a single ongoing trial subject to
sequential testing, and asks whether the accrued *information* has reached
the Required Information Size (RIS) and whether the cumulative test
statistic has crossed an O'Brien-Fleming monitoring boundary that controls
the inflated type-I error from repeated looks.

Method (advanced-stats.md gotchas applied):
  * Accrued information  I_acc = 1 / var(pooled RE)
  * RIS (fixed)          I_req = (z_{alpha/2} + z_beta)^2 / delta^2
        where delta is the effect to detect on the analysis scale (e.g. a
        log OR). Defaults to the observed RE estimate (post-hoc) with a
        caveat if not supplied.
  * Heterogeneity / Diversity adjustment (D), NOT the cluster design
        effect:  D = 1 + tau^2 * (sum(1/v_i^2)/(sum(1/v_i))^2 * k - 1);
        adjusted RIS = I_req * D
  * Information fraction  t_k = I_acc / RIS_adj
  * O'Brien-Fleming boundary  z_k = z_{alpha/2} / sqrt(t_k)   (two-sided)
  * Non-binding futility only (we do not act on a futility boundary).

Conclusion: firm evidence requires |z_cum| >= z_k AND t_k <= 1 reached, or
t_k >= 1 (RIS met). Otherwise "more information needed".
"""
from __future__ import annotations
import math
from . import common


def analyze(yis, vis, *, alpha=0.05, power=0.90, delta=None, ratio=False):
    k = len(yis)
    if k < 2:
        return {"available": False, "reason": "need >=2 studies"}

    tau2 = common.tau2_dersimonian_laird(yis, vis)
    wre = [1.0 / (v + tau2) for v in vis]
    swre = sum(wre)
    mu = sum(w * y for w, y in zip(wre, yis)) / swre
    var_re = 1.0 / swre
    se_re = math.sqrt(var_re)
    z_cum = mu / se_re if se_re else 0.0

    z_a = common.norm_ppf(1 - alpha / 2)
    z_b = common.norm_ppf(power)

    delta_used = delta if delta is not None else mu
    if delta_used == 0:
        return {"available": False, "reason": "effect to detect (delta) is 0"}

    I_acc = 1.0 / var_re
    I_req_fixed = (z_a + z_b) ** 2 / (delta_used ** 2)

    # Diversity / heterogeneity design effect (advanced-stats.md form)
    s1 = sum(1.0 / v for v in vis)
    s2 = sum(1.0 / (v * v) for v in vis)
    D = 1.0 + tau2 * (s2 / (s1 * s1) * k - 1.0) if s1 > 0 else 1.0
    D = max(1.0, D)
    I_req_adj = I_req_fixed * D

    t_k = I_acc / I_req_adj if I_req_adj > 0 else float("inf")
    t_k_capped = min(1.0, t_k)
    z_boundary = z_a / math.sqrt(t_k_capped) if t_k_capped > 0 else float("inf")

    ris_met = t_k >= 1.0
    crossed = abs(z_cum) >= z_boundary
    if ris_met:
        conclusion = ("firm evidence: required information size reached and "
                      "the cumulative effect is conclusive"
                      if abs(z_cum) >= z_a else
                      "required information size reached; effect not significant")
    elif crossed:
        conclusion = ("firm evidence: cumulative z crossed the O'Brien-Fleming "
                      "monitoring boundary before reaching RIS")
    else:
        conclusion = ("inconclusive: more information needed (boundary not "
                      "crossed and RIS not reached)")

    def disp(v):
        return math.exp(v) if ratio else v

    return {
        "available": True, "k": k,
        "z_cumulative": z_cum,
        "z_boundary": z_boundary,
        "alpha": alpha, "power": power,
        "delta": disp(delta_used), "delta_assumed": delta is None,
        "tau2": tau2, "diversity_D": D,
        "information_accrued": I_acc,
        "RIS_fixed": I_req_fixed, "RIS_adjusted": I_req_adj,
        "information_fraction": t_k,
        "ris_met": ris_met, "boundary_crossed": crossed,
        "conclusion": conclusion,
        "note": ("delta defaulted to the observed pooled effect (post-hoc) — "
                 "supply a pre-specified delta for a proper a-priori TSA."
                 if delta is None else ""),
    }
