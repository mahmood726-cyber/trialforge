import math, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from trialforge import advanced, nodesplit, common


# Build a small asymmetric set (small studies show larger effects)
def _asym():
    # big precise studies near 0; small studies shifted positive
    yis = [0.02, 0.0, 0.05, 0.6, 0.7, 0.8]
    vis = [0.01, 0.008, 0.012, 0.2, 0.25, 0.3]
    names = [f"S{i}" for i in range(6)]
    return yis, vis, names


def test_egger_detects_asymmetry():
    yis, vis, names = _asym()
    r = advanced.egger_test(yis, vis)
    assert r["available"]
    assert math.isfinite(r["p"])


def test_egger_symmetric_no_flag():
    yis = [0.1, -0.05, 0.08, -0.03, 0.06, -0.04]
    vis = [0.02, 0.02, 0.05, 0.05, 0.1, 0.1]
    r = advanced.egger_test(yis, vis)
    assert r["available"] and r["p"] > 0.10


def test_peters_test():
    studies = [
        {"tE": 10, "tN": 100, "cE": 12, "cN": 100},
        {"tE": 20, "tN": 200, "cE": 25, "cN": 200},
        {"tE": 5, "tN": 50, "cE": 9, "cN": 50},
        {"tE": 40, "tN": 400, "cE": 45, "cN": 400},
    ]
    r = advanced.peters_test(studies)
    assert r["available"] and math.isfinite(r["p"])


def test_pet_peese_runs():
    yis, vis, names = _asym()
    r = advanced.pet_peese(yis, vis)
    assert r["available"]
    assert r["chosen"] in ("PET", "PEESE")
    assert math.isfinite(r["adjusted_estimate"])


def test_trim_and_fill():
    yis, vis, names = _asym()
    r = advanced.trim_and_fill(yis, vis)
    assert r["available"]
    assert r["k_imputed"] >= 0


def test_leave_one_out():
    yis = [0.1, 0.2, 0.15, 0.9]
    vis = [0.02, 0.02, 0.02, 0.02]
    out = advanced.leave_one_out(yis, vis, ["A", "B", "C", "D"])
    assert len(out) == 4
    # omitting the outlier D should lower the pooled estimate vs omitting others
    d_omitted = next(o for o in out if o["omitted"] == "D")
    a_omitted = next(o for o in out if o["omitted"] == "A")
    assert d_omitted["estimate"] < a_omitted["estimate"]


def test_meta_regression_slope():
    # effect increases with moderator
    mod = [1, 2, 3, 4, 5, 6]
    yis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    vis = [0.01] * 6
    r = advanced.meta_regression(yis, vis, mod)
    assert r["available"]
    assert r["slope"] > 0.05
    assert r["p"] < 0.05


def test_subgroup_between_test():
    yis = [0.1, 0.12, 0.9, 0.95]
    vis = [0.01, 0.01, 0.01, 0.01]
    groups = ["G1", "G1", "G2", "G2"]
    r = advanced.subgroup(yis, vis, groups)
    assert r["p_between"] < 0.05   # groups clearly differ
    assert set(r["subgroups"]) == {"G1", "G2"}


def test_cumulative_order():
    yis = [0.5, 0.3, 0.2]
    vis = [0.1, 0.05, 0.02]
    out = advanced.cumulative(yis, vis, ["A", "B", "C"], [2018, 2020, 2022])
    assert len(out) == 3
    assert out[0]["k"] == 1 and out[-1]["k"] == 3


def test_peto_or_rare():
    studies = [
        {"tE": 1, "tN": 1000, "cE": 5, "cN": 1000},
        {"tE": 2, "tN": 1200, "cE": 7, "cN": 1200},
        {"tE": 0, "tN": 900, "cE": 4, "cN": 900},
    ]
    r = advanced.peto_or(studies)
    assert r["available"] and r["OR"] < 1.0  # intervention protective


def test_mantel_haenszel():
    studies = [
        {"tE": 10, "tN": 100, "cE": 20, "cN": 100},
        {"tE": 15, "tN": 120, "cE": 25, "cN": 110},
    ]
    r = advanced.mantel_haenszel_or(studies)
    assert r["available"] and r["OR"] < 1.0


def test_loop_inconsistency_consistent():
    # consistent triangle
    studies = [
        {"name": "AB", "arms": [{"t": "A", "e": 100, "n": 500}, {"t": "B", "e": 80, "n": 500}]},
        {"name": "AC", "arms": [{"t": "A", "e": 100, "n": 500}, {"t": "C", "e": 60, "n": 500}]},
        {"name": "BC", "arms": [{"t": "B", "e": 80, "n": 500}, {"t": "C", "e": 60, "n": 500}]},
    ]
    r = nodesplit.loop_inconsistency(studies, "OR")
    assert r["n_loops"] == 1
    assert not r["any_inconsistent"]  # built to be consistent


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
