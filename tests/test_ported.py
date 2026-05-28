"""Tests for the allmeta-ported modules: copas, dta (hsroc), cnma.
Reference values from allmeta R-parity fixtures (metasens / mada)."""
import math, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from trialforge import copas, dta, cnma, common


# ---- Copas (fixture from allmeta copas/tests/test_against_metasens.py) ----
_YI = [0.25, 0.18, 0.40, 0.30, 0.12, 0.55, 0.22, 0.32, 0.45, 0.28]
_SEI = [0.08, 0.10, 0.15, 0.09, 0.07, 0.20, 0.08, 0.12, 0.17, 0.10]
_VI = [s * s for s in _SEI]


def test_copas_fe_matches_metasens():
    """Unadjusted FE pool must match metasens te_fe = 0.246944262521 (tol 1e-4)."""
    r = copas.analyze(_YI, _VI)
    assert r["available"]
    assert math.isclose(r["fe_estimate"], 0.246944262521, abs_tol=1e-4)


def test_copas_attenuates_under_assumed_bias():
    """Worst-case (50% unpublished) should attenuate toward the null for this
    small-study-effect fixture, and the sensitivity slope should be negative."""
    r = copas.analyze(_YI, _VI)
    assert r["worst_case"]["estimate"] < r["re_estimate"]
    assert r["sensitivity_slope"] < 0


def test_copas_needs_three():
    r = copas.analyze([0.2, 0.3], [0.01, 0.01])
    assert not r["available"]


# ---- DTA / HSROC (fixture from allmeta hsroc/tests/test_against_mada.py) ----
_DTA_ROWS = [
    {"tp": 80, "fp": 40, "fn": 20, "tn": 160},
    {"tp": 95, "fp": 5, "fn": 5, "tn": 95},
    {"tp": 60, "fp": 30, "fn": 40, "tn": 170},
    {"tp": 88, "fp": 12, "fn": 12, "tn": 88},
    {"tp": 70, "fp": 50, "fn": 30, "tn": 150},
    {"tp": 92, "fp": 8, "fn": 8, "tn": 192},
    {"tp": 55, "fp": 45, "fn": 45, "tn": 155},
]


def test_dta_pooled_se_sp_sensible():
    r = dta.analyze(_DTA_ROWS)
    assert r["available"] and r["k"] == 7
    assert 0.5 < r["sensitivity"] < 0.95
    assert 0.5 < r["specificity"] < 0.95
    # mu logit(Se) > 0 (Se>0.5); mu logit(FPR) < 0 (FPR<0.5)
    assert r["mu_logit_se"] > 0 and r["mu_logit_fpr"] < 0
    assert r["dor"] > 1
    assert len(r["sroc"]) == 41


def test_dta_rho_constrained():
    r = dta.analyze(_DTA_ROWS)
    assert -0.95 <= r["rho"] <= 0.95


def test_dta_continuity_correction():
    rows = [{"tp": 0, "fp": 5, "fn": 10, "tn": 85},
            {"tp": 20, "fp": 4, "fn": 6, "tn": 90},
            {"tp": 30, "fp": 3, "fn": 5, "tn": 92}]
    r = dta.analyze(rows)
    assert r["available"]  # zero cell handled by +0.5


# ---- Component NMA ----------------------------------------------------------
def test_cnma_additive_recovery():
    """A network where combination effects are additive in components should
    recover each component's incremental effect."""
    # Reference = placebo (no components). Components A and B each lower odds.
    studies = [
        {"name": "S1", "arms": [
            {"components": [], "e": 120, "n": 500},
            {"components": ["A"], "e": 90, "n": 500}]},
        {"name": "S2", "arms": [
            {"components": [], "e": 115, "n": 480},
            {"components": ["B"], "e": 95, "n": 480}]},
        {"name": "S3", "arms": [
            {"components": ["A"], "e": 88, "n": 460},
            {"components": ["A", "B"], "e": 70, "n": 460}]},
    ]
    r = cnma.analyze(studies, "OR")
    assert r["available"]
    assert set(r["components"]) == {"A", "B"}
    # both components protective (OR < 1)
    for c in r["component_effects"]:
        assert c["effect"] < 1.2
    # predict A+B combination = product of ORs
    comb = cnma.predict_combination(r, ["A", "B"])
    assert comb is not None and comb < 1.0


def test_cnma_underidentified():
    studies = [{"name": "S1", "arms": [
        {"components": [], "e": 100, "n": 500},
        {"components": ["A", "B", "C"], "e": 80, "n": 500}]}]
    r = cnma.analyze(studies, "OR")
    assert not r["available"]  # 1 contrast, 3 components


if __name__ == "__main__":
    import traceback
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for fn in fns:
        try:
            fn(); print(f"  PASS {fn.__name__}"); passed += 1
        except Exception:
            print(f"  FAIL {fn.__name__}"); traceback.print_exc()
    print(f"\n{passed}/{len(fns)} passed")
    sys.exit(0 if passed == len(fns) else 1)
