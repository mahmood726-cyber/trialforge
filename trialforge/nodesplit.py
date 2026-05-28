"""trialforge.nodesplit — NMA inconsistency via Bucher closed-loop testing.

For every closed triangle of treatments that all have DIRECT head-to-head
evidence, compare the direct estimate of one edge with the indirect
estimate obtained from the other two edges (Bucher 1997). A large
inconsistency factor (IF) signals that direct and indirect evidence
disagree — the key validity check before trusting an NMA ranking.

Works on the same study format as trialforge.nma:
  studies: [{name, arms:[{t, e, n} | {t, m, sd, n}]}]
"""
from __future__ import annotations
import math
from itertools import combinations
from . import common, pairwise


def _direct_estimates(studies, measure):
    """Pool direct head-to-head evidence for each treatment pair.
    Returns {(a,b): {"est":logscale, "var":, "k":}} with a<b lexically."""
    pair_studies = {}
    for s in studies:
        arms = s["arms"]
        # all unordered arm pairs in this study contribute a contrast
        for i, j in combinations(range(len(arms)), 2):
            ai, aj = arms[i], arms[j]
            a_t, b_t = ai["t"], aj["t"]
            # orient so contrast is (b vs a) with a<b lexically
            if a_t == b_t:
                continue
            lo, hi = sorted([a_t, b_t])
            # build a 2-arm pairwise effect: treat hi as "intervention"
            int_arm, ctrl_arm = (aj, ai) if b_t == hi else (ai, aj)
            if measure == "OR":
                eff = pairwise.effect_OR(s.get("name", "s"),
                                         int_arm["e"], int_arm["n"],
                                         ctrl_arm["e"], ctrl_arm["n"])
            elif measure == "MD":
                eff = pairwise.effect_MD(s.get("name", "s"),
                                         int_arm["m"], int_arm["sd"], int_arm["n"],
                                         ctrl_arm["m"], ctrl_arm["sd"], ctrl_arm["n"])
            else:
                eff = None
            if eff is None:
                continue
            pair_studies.setdefault((lo, hi), []).append(eff)
    out = {}
    for pair, effs in pair_studies.items():
        pool = common.pool_inverse_variance([e.yi for e in effs],
                                            [e.vi for e in effs])
        out[pair] = {"est": pool.estimate, "var": pool.se ** 2, "k": pool.k}
    return out


def _edge(direct, a, b):
    """Return (est, var) for the contrast b vs a (lexical orientation handled)."""
    lo, hi = sorted([a, b])
    if (lo, hi) not in direct:
        return None
    e = direct[(lo, hi)]
    est, var = e["est"], e["var"]
    # direct stores (hi vs lo); flip if we asked for (a=hi, b=lo) i.e. lo vs hi
    if a == lo and b == hi:
        return est, var          # b(hi) vs a(lo) == stored
    return -est, var             # a(hi) vs b(lo): negate


def loop_inconsistency(studies, measure="OR"):
    """Return Bucher inconsistency factors for all closed triangles."""
    direct = _direct_estimates(studies, measure)
    treatments = sorted({t for s in studies for t in (a["t"] for a in s["arms"])})
    ratio = measure == "OR"
    loops = []
    for A, B, C in combinations(treatments, 3):
        # need direct evidence on all three edges
        if not all(tuple(sorted(p)) in direct for p in [(A, B), (A, C), (B, C)]):
            continue
        # direct B vs A
        dir_BA = _edge(direct, A, B)
        # indirect B vs A via C = (B vs C) - (A vs C) ... = (B vs C)+(C vs A)
        bc = _edge(direct, C, B)   # B vs C
        ca = _edge(direct, A, C)   # C vs A  -> need C vs A
        # _edge(direct, A, C) returns C vs A? it returns (b=C vs a=A) = C vs A
        if dir_BA is None or bc is None or ca is None:
            continue
        ind_est = bc[0] + ca[0]
        ind_var = bc[1] + ca[1]
        IF = dir_BA[0] - ind_est
        IF_var = dir_BA[1] + ind_var
        se = math.sqrt(IF_var) if IF_var > 0 else float("nan")
        z = IF / se if se else float("nan")
        p = 2 * common.norm_sf(abs(z)) if math.isfinite(z) else float("nan")
        loops.append({
            "loop": f"{A}-{B}-{C}",
            "edge": f"{B} vs {A}",
            "direct": math.exp(dir_BA[0]) if ratio else dir_BA[0],
            "indirect": math.exp(ind_est) if ratio else ind_est,
            "IF": math.exp(IF) if ratio else IF,
            "IF_raw": IF, "z": z, "p": p,
            "inconsistent": (p < 0.05) if math.isfinite(p) else False,
        })
    return {"n_loops": len(loops), "loops": loops,
            "any_inconsistent": any(l["inconsistent"] for l in loops),
            "n_direct_edges": len(direct)}
