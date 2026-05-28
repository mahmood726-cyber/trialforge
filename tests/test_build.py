import subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent


def _build(name, tmp_path):
    out = tmp_path / (name + ".html")
    r = subprocess.run(
        [sys.executable, str(ROOT / "run.py"),
         str(ROOT / "configs" / (name + ".json")), "--out", str(out)],
        capture_output=True, text=True, timeout=180)
    assert r.returncode == 0, r.stderr
    html = out.read_text(encoding="utf-8")
    for bad in (">None<", "None%", ">nan<", "${"):
        assert bad not in html, f"{bad!r} leaked into {name}"
    assert "<svg" in html
    return html


def test_pairwise_advanced(tmp_path):
    h = _build("example_pairwise_advanced", tmp_path)
    assert "Advanced diagnostics" in h
    assert "Egger" in h
    assert "Leave-one-out" in h
    assert "Cumulative" in h
    assert "Meta-regression" in h
    assert "Copas" in h


def test_dta(tmp_path):
    h = _build("example_dta", tmp_path)
    assert "SROC" in h
    assert "Sensitivity" in h


def test_cnma(tmp_path):
    h = _build("example_cnma", tmp_path)
    assert "component" in h.lower()
    assert "CBT" in h


def test_nma_inconsistency(tmp_path):
    h = _build("example_nma_inconsistency", tmp_path)
    assert "SUCRA" in h
    assert "inconsistency" in h.lower()


def test_rareevents(tmp_path):
    h = _build("example_rareevents", tmp_path)
    assert "GLMM" in h


def test_multivariate(tmp_path):
    h = _build("example_multivariate", tmp_path)
    assert "borrow" in h.lower()


def test_rmst(tmp_path):
    h = _build("example_rmst", tmp_path)
    assert "RMST" in h


def test_grade_in_pairwise(tmp_path):
    h = _build("example_pairwise_advanced", tmp_path)
    assert "GRADE" in h
    assert "certainty" in h.lower()


def test_webr_panel_in_pairwise(tmp_path):
    # The report must be capsule-like: it carries an in-browser WebR/metafor
    # cross-validation panel that re-pools the SAME yi/sei.
    h = _build("example_sglt2_hf", tmp_path)
    assert "Verify in R (WebR)" in h
    assert 'id="tf-webr-data"' in h
    assert 'id="tf-webr-btn"' in h
    assert "webr.r-wasm.org/latest/webr.mjs" in h
    assert "metafor" in h
    # the embedded JSON must carry trialforge's own pooled estimate to diff against
    import re, json
    m = re.search(r'<script type="application/json" id="tf-webr-data">(.*?)</script>',
                  h, re.S)
    assert m, "WebR data tag missing"
    d = json.loads(m.group(1).replace("\\u003c", "<"))
    assert d["measure"] == "HR" and d["method"] == "PM" and d["ratio"] is True
    assert len(d["yi"]) == len(d["sei"]) == 5
    assert abs(d["tf"]["estimate"] - 0.771) < 0.01


def test_webr_only_for_pairwise(tmp_path):
    # DTA/RMST reports are not pairwise pools — no WebR panel there.
    h = _build("example_rmst", tmp_path)
    assert "tf-webr-data" not in h


if __name__ == "__main__":
    import tempfile, traceback
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for fn in fns:
        try:
            with tempfile.TemporaryDirectory() as td:
                fn(Path(td))
            print(f"  PASS {fn.__name__}"); passed += 1
        except Exception:
            print(f"  FAIL {fn.__name__}"); traceback.print_exc()
    print(f"\n{passed}/{len(fns)} passed")
    sys.exit(0 if passed == len(fns) else 1)
