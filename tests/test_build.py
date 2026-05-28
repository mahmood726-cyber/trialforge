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


def test_nma_inconsistency(tmp_path):
    h = _build("example_nma_inconsistency", tmp_path)
    assert "SUCRA" in h
    assert "inconsistency" in h.lower()


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
