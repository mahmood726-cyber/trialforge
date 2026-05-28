"""trialforge.cnma — additive component network meta-analysis.

Pure-stdlib re-implementation of the allmeta `component-nma` engine.
Multi-component interventions (e.g. "A+B", "A+C") are decomposed into
their components; the additive CNMA model assumes the effect of a
combination equals the sum of its component effects:
    d(combination vs reference) = sum_c (component_c effect)

We build a study-contrast design matrix on COMPONENTS (not whole
treatments) and fit by weighted least squares, exactly as for a standard
NMA but with the component design. This yields each component's
incremental effect and lets us predict any (even unobserved) combination.

Input studies: [{name, arms:[{components:[...], e, n} | {components:[...], m, sd, n}]}]
Each arm lists its components; the reference is the empty component set
(e.g. placebo) or a chosen baseline.
"""
from __future__ import annotations
import math
from . import common, linalg, pairwise


def _arm_effect(arm, measure):
    if measure == "OR":
        a = arm["e"] + 0.5
        b = arm["n"] - arm["e"] + 0.5
        return math.log(a / b), 1 / a + 1 / b
    elif measure == "MD":
        return arm["m"], (arm["sd"] ** 2) / arm["n"]
    return None


def analyze(studies, measure="OR"):
    # collect components
    components = []
    for s in studies:
        for arm in s["arms"]:
            for c in arm.get("components", []):
                if c not in components:
                    components.append(c)
    components.sort()
    C = len(components)
    if C < 1:
        return {"available": False, "reason": "no components found"}
    cindex = {c: i for i, c in enumerate(components)}

    # build contrasts vs each study's first arm
    X_rows, y_rows, var_rows = [], [], []
    for s in studies:
        arms = s["arms"]
        if len(arms) < 2:
            continue
        base = arms[0]
        be = _arm_effect(base, measure)
        if be is None:
            continue
        base_comps = set(base.get("components", []))
        for arm in arms[1:]:
            ae = _arm_effect(arm, measure)
            if ae is None:
                continue
            y = ae[0] - be[0]
            v = ae[1] + be[1]
            # design row: +1 for components added, -1 for components removed
            arm_comps = set(arm.get("components", []))
            row = [0.0] * C
            for c in arm_comps - base_comps:
                row[cindex[c]] += 1.0
            for c in base_comps - arm_comps:
                row[cindex[c]] -= 1.0
            if not any(row):
                continue
            X_rows.append(row)
            y_rows.append(y)
            var_rows.append(v)

    if len(X_rows) < C:
        return {"available": False,
                "reason": f"under-identified: {len(X_rows)} contrasts for {C} components"}

    w = [1.0 / v for v in var_rows]
    XtWX = linalg.quad_form_diag(X_rows, w)
    XtWy = linalg.xt_w_y(X_rows, w, y_rows)
    try:
        cov = linalg.inv(XtWX)
    except ValueError:
        return {"available": False, "reason": "component design is singular "
                "(components not separately identifiable from this network)"}
    beta = linalg.matvec(cov, XtWy)

    ratio = measure == "OR"
    def disp(v):
        return math.exp(v) if ratio else v

    comp_effects = []
    for i, c in enumerate(components):
        se = math.sqrt(cov[i][i]) if cov[i][i] > 0 else 0.0
        comp_effects.append({
            "component": c,
            "effect": disp(beta[i]),
            "ci_low": disp(beta[i] - common.Z975 * se),
            "ci_high": disp(beta[i] + common.Z975 * se),
            "se": se,
        })
    comp_effects.sort(key=lambda r: r["effect"])
    return {
        "available": True, "measure": measure, "ratio": ratio,
        "components": components,
        "n_contrasts": len(X_rows),
        "component_effects": comp_effects,
        "note": "Additive CNMA: combination effect = sum of component effects "
                "(assumes no component interaction).",
    }


def predict_combination(result, comps):
    """Predict the additive effect of a component combination."""
    if not result.get("available"):
        return None
    lookup = {c["component"]: c for c in result["component_effects"]}
    ratio = result["ratio"]
    total = 0.0
    for c in comps:
        if c not in lookup:
            return None
        eff = lookup[c]["effect"]
        total += math.log(eff) if ratio else eff
    return math.exp(total) if ratio else total
