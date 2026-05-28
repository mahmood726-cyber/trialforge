"""Tests for p-curve, GOSH, CINeMA (allmeta-ported batch 3).
References from allmeta p-curve & gosh R-parity fixtures."""
import math, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from trialforge import pcurve, gosh, cinema, nma, nodesplit


# ---- p-curve (allmeta pcurve-tiny: 5 p-values) ----
def test_pcurve_matches_r():
    r = pcurve.analyze([0.003, 0.012, 0.021, 0.038, 0.047])
    assert r["available"] and r["n_significant"] == 5
    assert math.isclose(r["fisher_chisq"], 10.88867977905, abs_tol=1e-5)
    assert r["fisher_df"] == 10
    assert math.isclose(r["fisher_p"], 0.3662564540401, abs_tol=1e-6)
    assert math.isclose(r["prop_low"], 0.6, abs_tol=1e-9)
    assert math.isclose(r["z_flat"], 0.4472135955, abs_tol=1e-6)
    assert math.isclose(r["p_flat"], 0.6547208460186, abs_tol=1e-6)


def test_pcurve_evidential_value():
    # all very small p -> strong right skew -> significant Fisher
    r = pcurve.analyze([0.001, 0.002, 0.001, 0.004, 0.003])
    assert r["evidential_value"]


def test_pcurve_drops_nonsignificant():
    r = pcurve.analyze([0.01, 0.02, 0.5, 0.8])
    assert r["n_significant"] == 2


# ---- GOSH (allmeta gosh-tiny: 5 studies) ----
_Y = [0.22, 0.30, 0.27, 0.15, 0.45]
_V = [0.0064, 0.0196, 0.0081, 0.0049, 0.0324]


def test_gosh_full_estimate_matches_metafor():
    r = gosh.analyze(_Y, _V, method="FE")
    assert r["available"]
    assert r["n_subsets"] == 26          # all size>=2 subsets of 5
    assert math.isclose(r["full_estimate"], 0.2254227851625, abs_tol=1e-6)
    assert math.isclose(r["median_estimate"], 0.2217261558338, abs_tol=1e-4)
    assert r["min_estimate"] <= r["median_estimate"] <= r["max_estimate"]
    assert r["full_i2"] == 0.0


def test_gosh_samples_large_k():
    y = [0.1 * i for i in range(20)]
    v = [0.05] * 20
    r = gosh.analyze(y, v, method="FE", n_samples=1000)
    assert r["sampled"] and r["n_subsets"] == 1000


# ---- CINeMA ----
def _nma_res():
    studies = [
        {"name": "S1", "arms": [{"t": "Placebo", "e": 120, "n": 500}, {"t": "A", "e": 90, "n": 500}]},
        {"name": "S2", "arms": [{"t": "Placebo", "e": 110, "n": 480}, {"t": "B", "e": 70, "n": 470}]},
        {"name": "S3", "arms": [{"t": "A", "e": 85, "n": 460}, {"t": "B", "e": 72, "n": 455}]},
    ]
    return studies, nma.analyze(studies, measure="OR", reference="Placebo")


def test_cinema_rates_each_comparison():
    studies, res = _nma_res()
    loops = nodesplit.loop_inconsistency(studies, "OR")
    cm = cinema.rate(res, loops=loops, ratio=True)
    assert cm["available"]
    assert len(cm["comparisons"]) == 3   # A-B, A-Placebo, B-Placebo
    for c in cm["comparisons"]:
        assert c["confidence"] in ("high", "moderate", "low", "very low")
        # data-driven domains are assessed
        assert c["domains"]["imprecision"] != "not assessed"
        assert c["domains"]["heterogeneity"] != "not assessed"


def test_cinema_judgments_downgrade():
    studies, res = _nma_res()
    cm = cinema.rate(res, ratio=True,
                     judgments={"A vs Placebo": {"within_study": "major concerns"}})
    ap = next(c for c in cm["comparisons"] if c["comparison"] == "A vs Placebo")
    assert ap["domains"]["within-study bias"] == "major concerns"
    assert ap["confidence"] in ("low", "very low")


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
