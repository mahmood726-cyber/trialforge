import json, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from trialforge import checks


def _codes(issues):
    return {x["code"] for x in issues}


def test_clean_pairwise_passes():
    cfg = {"type": "pairwise", "measure": "OR", "studies": [
        {"name": "A", "tE": 30, "tN": 300, "cE": 40, "cN": 300},
        {"name": "B", "tE": 25, "tN": 280, "cE": 35, "cN": 290},
        {"name": "C", "tE": 20, "tN": 250, "cE": 30, "cN": 250},
    ]}
    issues, s = checks.check(cfg)
    assert s["ok"] and s["errors"] == 0


def test_impossible_count_is_error():
    cfg = {"type": "pairwise", "studies": [
        {"name": "Bad", "tE": 95, "tN": 22, "cE": 3, "cN": 79},
        {"name": "Ok", "tE": 10, "tN": 100, "cE": 12, "cN": 100},
    ]}
    issues, s = checks.check(cfg)
    assert not s["ok"]
    assert "impossible_count" in _codes(issues)


def test_inverted_ci_is_error():
    cfg = {"type": "pairwise", "measure": "HR", "studies": [
        {"name": "A", "effect": 0.8, "ci_low": 0.95, "ci_high": 0.7},
        {"name": "B", "effect": 0.9, "ci_low": 0.8, "ci_high": 1.0},
    ]}
    issues, s = checks.check(cfg)
    assert "inverted_ci" in _codes(issues)


def test_nonpositive_ratio_error():
    cfg = {"type": "pairwise", "measure": "OR", "studies": [
        {"name": "A", "effect": 0.0, "ci_low": -0.1, "ci_high": 0.5},
        {"name": "B", "effect": 0.8, "ci_low": 0.6, "ci_high": 1.0},
    ]}
    issues, s = checks.check(cfg)
    assert "nonpositive_ratio" in _codes(issues)


def test_pi_undefined_warning_small_k():
    cfg = {"type": "pairwise", "studies": [
        {"name": "A", "tE": 10, "tN": 100, "cE": 12, "cN": 100},
        {"name": "B", "tE": 9, "tN": 95, "cE": 13, "cN": 99},
    ]}
    issues, s = checks.check(cfg)
    assert "pi_undefined" in _codes(issues)
    assert s["ok"]  # warning, not error


def test_pubbias_low_power_warning():
    cfg = {"type": "pairwise", "advanced": ["egger"], "studies": [
        {"name": f"S{i}", "tE": 10, "tN": 100, "cE": 12, "cN": 100} for i in range(5)
    ]}
    issues, s = checks.check(cfg)
    assert "pubbias_low_power" in _codes(issues)


def test_duplicate_names_and_bad_nct():
    cfg = {"type": "pairwise", "studies": [
        {"name": "Dup", "nct": "NCT123", "tE": 10, "tN": 100, "cE": 12, "cN": 100},
        {"name": "Dup", "tE": 9, "tN": 95, "cE": 13, "cN": 99},
        {"name": "X", "tE": 8, "tN": 90, "cE": 11, "cN": 90},
    ]}
    issues, s = checks.check(cfg)
    assert "duplicate_name" in _codes(issues)
    assert "nct_format" in _codes(issues)


def test_no_studies_error():
    issues, s = checks.check({"type": "nma"})
    assert not s["ok"] and "no_studies" in _codes(issues)


def test_check_flag_blocks_bad_config(tmp_path):
    cfg = {"type": "pairwise", "studies": [
        {"name": "Bad", "tE": 95, "tN": 22, "cE": 3, "cN": 79},
        {"name": "Ok", "tE": 10, "tN": 100, "cE": 12, "cN": 100},
    ]}
    cp = tmp_path / "bad.json"
    cp.write_text(json.dumps(cfg), encoding="utf-8")
    r = subprocess.run([sys.executable, str(ROOT / "run.py"), str(cp), "--check"],
                       capture_output=True, text=True, timeout=60)
    assert r.returncode == 2  # errors -> exit 2
    assert "impossible_count" in (r.stdout + r.stderr)


def test_check_flag_passes_clean(tmp_path):
    cfg = {"type": "pairwise", "measure": "OR", "studies": [
        {"name": "A", "tE": 30, "tN": 300, "cE": 40, "cN": 300},
        {"name": "B", "tE": 25, "tN": 280, "cE": 35, "cN": 290},
        {"name": "C", "tE": 20, "tN": 250, "cE": 30, "cN": 250},
    ]}
    cp = tmp_path / "good.json"
    cp.write_text(json.dumps(cfg), encoding="utf-8")
    r = subprocess.run([sys.executable, str(ROOT / "run.py"), str(cp), "--check"],
                       capture_output=True, text=True, timeout=60)
    assert r.returncode == 0


def test_build_blocks_on_error(tmp_path):
    cfg = {"type": "pairwise", "studies": [
        {"name": "Bad", "tE": 95, "tN": 22, "cE": 3, "cN": 79},
        {"name": "Ok", "tE": 10, "tN": 100, "cE": 12, "cN": 100},
    ]}
    cp = tmp_path / "bad2.json"
    cp.write_text(json.dumps(cfg), encoding="utf-8")
    out = tmp_path / "x.html"
    r = subprocess.run([sys.executable, str(ROOT / "run.py"), str(cp), "--out", str(out)],
                       capture_output=True, text=True, timeout=60)
    assert r.returncode == 2  # build blocked by data error
    assert not out.exists()


if __name__ == "__main__":
    import traceback, tempfile
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    p = 0
    for fn in fns:
        try:
            import inspect
            if "tmp_path" in inspect.signature(fn).parameters:
                with tempfile.TemporaryDirectory() as td:
                    fn(Path(td))
            else:
                fn()
            print(f"  PASS {fn.__name__}"); p += 1
        except Exception:
            print(f"  FAIL {fn.__name__}"); traceback.print_exc()
    print(f"\n{p}/{len(fns)} passed")
