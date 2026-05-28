"""Tests for glmm, multivariate, survival/RMST, grade (allmeta-ported batch 4)."""
import math, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from trialforge import glmm, multivariate, survival, grade, advanced


# ---- rare-events GLMM ----
def test_glmm_rare_events_protective():
    studies = [
        {"tE": 1, "tN": 1000, "cE": 5, "cN": 1000},
        {"tE": 0, "tN": 800,  "cE": 4, "cN": 820},
        {"tE": 2, "tN": 1200, "cE": 7, "cN": 1180},
        {"tE": 0, "tN": 600,  "cE": 3, "cN": 610},
        {"tE": 3, "tN": 1500, "cE": 9, "cN": 1490},
    ]
    r = glmm.analyze(studies)
    assert r["available"]
    assert r["OR"] < 1.0          # intervention protective
    # roughly agrees with Peto on rare balanced data
    peto = advanced.peto_or(studies)
    assert abs(r["OR"] - peto["OR"]) < 0.15


def test_glmm_handles_zero_cells():
    # studies with zero treatment-arm events must not crash / need no +0.5
    studies = [
        {"tE": 0, "tN": 500, "cE": 6, "cN": 500},
        {"tE": 0, "tN": 450, "cE": 5, "cN": 460},
        {"tE": 1, "tN": 520, "cE": 8, "cN": 510},
    ]
    r = glmm.analyze(studies)
    assert r["available"] and r["OR"] < 1.0


def test_glmm_needs_informative():
    studies = [{"tE": 0, "tN": 100, "cE": 0, "cN": 100}]
    assert not glmm.analyze(studies)["available"]


# ---- multivariate borrowing of strength ----
def test_multivariate_borrows_for_missing_outcome():
    studies = [
        {"y1": -1.2, "se1": 0.30, "y2": -0.9, "se2": 0.35, "rho_w": 0.6},
        {"y1": -1.0, "se1": 0.28, "y2": -0.8, "se2": 0.33, "rho_w": 0.6},
        {"y1": -1.5, "se1": 0.40, "y2": -1.1, "se2": 0.45, "rho_w": 0.6},
        {"y1": -0.8, "se1": 0.25, "y2": None, "se2": None, "rho_w": 0.6},
        {"y1": -1.3, "se1": 0.35, "y2": -1.0, "se2": 0.40, "rho_w": 0.6},
        {"y1": -1.1, "se1": 0.30, "y2": None, "se2": None, "rho_w": 0.6},
    ]
    r = multivariate.analyze(studies)
    assert r["available"]
    # outcome 2 (2 missing) should gain precision from outcome 1
    assert r["outcome2"]["se_bivariate"] <= r["outcome2"]["se_univariate"] + 1e-9
    assert r["outcome2"]["borrowed_precision_pct"] >= 0
    assert -0.95 <= r["rho_between"] <= 0.95


def test_multivariate_needs_enough_data():
    studies = [{"y1": 1, "se1": 0.2, "y2": None, "se2": None}]
    assert not multivariate.analyze(studies)["available"]


# ---- RMST / survival ----
def test_rmst_from_km_area():
    # survival 1.0 to t=10 then 0.5 to t=20; tau=20
    pts = [(10, 0.5), (20, 0.5)]
    area = survival.rmst_from_km(pts, 20)
    # 1.0*10 + 0.5*10 = 15
    assert abs(area - 15.0) < 1e-9


def test_rmst_pool():
    studies = [
        {"name": "A", "rmst_diff": 1.8, "se": 0.6},
        {"name": "B", "rmst_diff": 2.4, "se": 0.8},
        {"name": "C", "rmst_diff": 1.2, "se": 0.5},
        {"name": "D", "rmst_diff": 3.0, "ci_low": 1.4, "ci_high": 4.6},
    ]
    r = survival.analyze(studies, tau=24)
    assert r["available"] and r["k"] == 4
    assert 1.0 < r["rmst_difference"] < 3.0
    assert r["ci_low"] < r["rmst_difference"] < r["ci_high"]


# ---- GRADE ----
def test_grade_high_when_clean():
    r = grade.rate(measure="RR", estimate=0.75, ci_low=0.65, ci_high=0.86,
                   k=12, n_total=8000, i2=10, baseline_risk=0.20, egger_p=0.4,
                   risk_of_bias="not serious", indirectness="not serious")
    assert r["certainty"] == "high"
    assert r["absolute_effect"]["difference_per_1000"] < 0  # protective


def test_grade_downgrades():
    r = grade.rate(measure="RR", estimate=0.9, ci_low=0.6, ci_high=1.35,
                   k=4, n_total=200, i2=80, baseline_risk=0.10,
                   risk_of_bias="serious (-1)", indirectness="not serious")
    # high heterogeneity (-2) + imprecision (crosses null & small, -2) + RoB (-1)
    assert r["certainty"] == "very low"


def test_grade_sof_absolute_effect():
    r = grade.rate(measure="RR", estimate=0.5, ci_low=0.4, ci_high=0.65,
                   k=8, n_total=5000, i2=20, baseline_risk=0.20)
    ae = r["absolute_effect"]
    assert ae["comparator_per_1000"] == 200
    assert ae["intervention_per_1000"] == 100   # 0.5 * 200
    assert ae["difference_per_1000"] == -100


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
