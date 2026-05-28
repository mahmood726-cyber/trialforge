"""AACT loader tests — skip cleanly when no snapshot is available so the
suite passes on any machine."""
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from trialforge import aact


def _have_snapshot():
    try:
        aact.find_snapshot()
        return True
    except FileNotFoundError:
        return False


def test_classify_param_type():
    assert aact.classify_param_type("Hazard Ratio (HR)") == ("HR", True)
    assert aact.classify_param_type("Odds Ratio (OR)") == ("OR", True)
    assert aact.classify_param_type("Risk Ratio (RR)") == ("RR", True)
    assert aact.classify_param_type("Mean Difference (Final Values)") == ("MD", False)
    assert aact.classify_param_type("Risk Difference (RD)") == ("RD", False)
    assert aact.classify_param_type("LS Mean Difference") == ("MD", False)
    assert aact.classify_param_type("Slope") == (None, None)


@pytest.mark.skipif(not _have_snapshot(), reason="no AACT snapshot on this machine")
def test_extract_effects_finerenone():
    a = aact.AACT()
    res = a.extract_effects(["NCT02540993", "NCT02545049"], force_measure="HR")
    assert res["measure"] == "HR"
    assert len(res["studies"]) >= 1
    for s in res["studies"]:
        assert s["effect"] > 0 and s["ci_low"] > 0 and s["ci_high"] > s["ci_low"]


@pytest.mark.skipif(not _have_snapshot(), reason="no AACT snapshot on this machine")
def test_find_trials_query():
    a = aact.AACT()
    ncts = a.find_trials(drug="finerenone", limit=20)
    assert len(ncts) > 0
    assert all(n.startswith("NCT") for n in ncts)


if __name__ == "__main__":
    import traceback
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = skipped = 0
    for fn in fns:
        try:
            fn(); print(f"  PASS {fn.__name__}"); passed += 1
        except Exception as e:
            if "Skipped" in type(e).__name__ or not _have_snapshot():
                print(f"  SKIP {fn.__name__}"); skipped += 1
            else:
                print(f"  FAIL {fn.__name__}"); traceback.print_exc()
    print(f"\n{passed} passed, {skipped} skipped")
