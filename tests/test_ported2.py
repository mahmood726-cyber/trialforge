"""Tests for limit-MA, TSA, E-value (allmeta-ported batch 2).
Limit-MA reference from allmeta limit-ma/tests (metasens::limitmeta)."""
import math, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from trialforge import limitma, tsa, evalue

# limit-ma fixture (allmeta limit-tiny): te / se pairs
_TE = [0.55, 0.48, 0.70, 0.60, 0.30, 0.75, 0.52, 0.65, 0.45, 0.58]
_SE = [0.08, 0.10, 0.15, 0.09, 0.07, 0.20, 0.08, 0.12, 0.17, 0.10]
_VI = [s * s for s in _SE]


def test_limitma_baselines_match_metasens():
    r = limitma.analyze(_TE, _VI)
    assert r["available"]
    assert math.isclose(r["fe_estimate"], 0.5121282168817, abs_tol=1e-6)
    assert math.isclose(r["re_estimate"], 0.5300923657391, abs_tol=1e-6)
    assert math.isclose(r["tau2"], 0.007231869265627, abs_tol=1e-6)


def test_limitma_attenuates_and_within_tol():
    r = limitma.analyze(_TE, _VI)
    assert r["limit_estimate"] < r["re_estimate"]
    assert r["slope"] > 0  # small-study effect
    # JS WLS approximation within 0.15 of metasens te_adj=0.412
    assert math.isclose(r["limit_estimate"], 0.411998010092, abs_tol=0.15)


def test_limitma_needs_three():
    assert not limitma.analyze([0.2, 0.3], [0.01, 0.01])["available"]


def test_tsa_inconclusive_when_underpowered():
    # 3 small studies, modest effect -> RIS not reached
    yis = [0.1, 0.12, 0.08]
    vis = [0.05, 0.06, 0.05]
    r = tsa.analyze(yis, vis, delta=0.2)
    assert r["available"]
    assert r["information_fraction"] < 1.0
    assert "inconclusive" in r["conclusion"] or "more information" in r["conclusion"]


def test_tsa_firm_when_large_consistent():
    # many precise consistent studies -> RIS reached, significant
    yis = [0.5] * 8
    vis = [0.01] * 8
    r = tsa.analyze(yis, vis, delta=0.5)
    assert r["available"]
    assert r["z_cumulative"] > 3
    assert r["ris_met"] or r["boundary_crossed"]


def test_tsa_diversity_ge_one():
    yis = [0.3, 0.6, 0.1, 0.5]
    vis = [0.02, 0.03, 0.02, 0.04]
    r = tsa.analyze(yis, vis, delta=0.3)
    assert r["diversity_D"] >= 1.0


def test_evalue_basic():
    # RR 2.0 -> E-value = 2 + sqrt(2*1) = 3.414
    r = evalue.analyze(2.0, 1.5, 2.7, measure="RR")
    assert math.isclose(r["evalue_point"], 2 + math.sqrt(2), abs_tol=1e-6)
    assert r["evalue_ci"] > 1.0  # CI excludes null
    assert not r["ci_crosses_null"]


def test_evalue_protective_rr():
    # RR 0.5 -> treat as 1/0.5=2 -> same E-value
    r = evalue.analyze(0.5, 0.4, 0.65, measure="RR")
    assert math.isclose(r["evalue_point"], 2 + math.sqrt(2), abs_tol=1e-6)


def test_evalue_null_ci():
    r = evalue.analyze(1.2, 0.9, 1.6, measure="RR")
    assert r["ci_crosses_null"]
    assert r["evalue_ci"] == 1.0


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
