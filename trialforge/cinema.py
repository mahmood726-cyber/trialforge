"""trialforge.cinema — Confidence In Network Meta-Analysis (CINeMA-style).

A structured per-comparison confidence rating inspired by the CINeMA
framework (Nikolakopoulou 2020). Six domains:
  1. within-study bias     (judgment — user supplies, else 'not assessed')
  2. reporting bias        (judgment — user supplies, else 'not assessed')
  3. indirectness          (judgment — user supplies, else 'not assessed')
  4. imprecision           (DATA: CI width vs a clinical threshold)
  5. heterogeneity         (DATA: network tau / prediction interval)
  6. incoherence           (DATA: Bucher loop inconsistency, if available)

Each domain -> no / some / major concerns. Overall confidence starts at
'high' and is downgraded GRADE-style: one 'some' -> moderate; one 'major'
-> low; two or more 'major' -> very low.

This automates the data-driven domains (imprecision, heterogeneity,
incoherence) and surfaces the judgment domains for the reviewer rather
than guessing them — an honest partial CINeMA.
"""
from __future__ import annotations
import math

_LEVELS = ["high", "moderate", "low", "very low"]


def _downgrade(domains):
    some = sum(1 for d in domains.values() if d == "some concerns")
    major = sum(1 for d in domains.values() if d == "major concerns")
    steps = some * 1 + major * 2
    idx = min(len(_LEVELS) - 1, steps)
    return _LEVELS[idx]


def rate(nma_result, *, loops=None, judgments=None, ratio=True,
         imprecision_threshold=None, null=None):
    """nma_result: output of trialforge.nma.analyze.
    loops: output of trialforge.nodesplit.loop_inconsistency (optional).
    judgments: {comparison_key: {within_study, reporting, indirectness}}.
    """
    judgments = judgments or {}
    null = null if null is not None else (1.0 if ratio else 0.0)
    tau = nma_result.get("tau", 0.0)
    # loop inconsistency lookup by treatment pair
    incoh = {}
    if loops and loops.get("loops"):
        for lp in loops["loops"]:
            edge = lp["edge"]  # "B vs A"
            incoh[edge] = lp

    out = []
    for pw in nma_result["pairwise"]:
        if pw["a"] >= pw["b"]:
            continue
        key = f"{pw['a']} vs {pw['b']}"
        lo, hi = pw["ci_low"], pw["ci_high"]
        # imprecision: CI crosses null, or very wide
        if ratio:
            width_ratio = hi / lo if lo > 0 else float("inf")
            crosses = lo <= null <= hi
            wide = width_ratio > (imprecision_threshold or 3.0)
        else:
            crosses = lo <= null <= hi
            wide = (hi - lo) > (imprecision_threshold or 1.0)
        if crosses and wide:
            imprecision = "major concerns"
        elif crosses or wide:
            imprecision = "some concerns"
        else:
            imprecision = "no concerns"
        # heterogeneity from network tau
        if tau > 0.5:
            heterogeneity = "major concerns"
        elif tau > 0.2:
            heterogeneity = "some concerns"
        else:
            heterogeneity = "no concerns"
        # incoherence
        incoherence = "not assessed"
        for ek, lp in incoh.items():
            ta, tb = [t.strip() for t in ek.split(" vs ")]
            if {ta, tb} == {pw["a"], pw["b"]}:
                incoherence = ("major concerns" if lp["inconsistent"]
                               else "no concerns")
                break
        # judgments
        j = judgments.get(key, {})
        domains = {
            "within-study bias": j.get("within_study", "not assessed"),
            "reporting bias": j.get("reporting", "not assessed"),
            "indirectness": j.get("indirectness", "not assessed"),
            "imprecision": imprecision,
            "heterogeneity": heterogeneity,
            "incoherence": incoherence,
        }
        # only count assessed domains for downgrading
        assessed = {k: v for k, v in domains.items() if v != "not assessed"}
        confidence = _downgrade(assessed)
        out.append({
            "comparison": key,
            "estimate": pw["estimate"], "ci_low": lo, "ci_high": hi,
            "domains": domains,
            "confidence": confidence,
        })
    return {
        "available": True,
        "comparisons": out,
        "note": "Partial CINeMA: imprecision, heterogeneity and incoherence "
                "are computed from the data; within-study bias, reporting bias "
                "and indirectness require reviewer judgment (supply via "
                "'cinema_judgments' or they stay 'not assessed').",
    }
