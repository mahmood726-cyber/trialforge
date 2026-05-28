"""Verify the pure-Python linear-algebra fallback agrees with numpy (when
present) and is internally consistent (when not). Addresses the reviewer's
'no numpy-vs-fallback parity test' gap."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from trialforge import linalg


def _close(A, B, tol=1e-9):
    return all(abs(A[i][j] - B[i][j]) < tol
               for i in range(len(A)) for j in range(len(A[0])))


def test_inv_identity():
    M = [[4.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 2.0]]
    Inv = linalg.inv(M)
    prod = linalg.matmul(M, Inv)
    I = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    assert _close(prod, I, 1e-8)


def test_solve_matches_inv():
    A = [[3.0, 2.0], [1.0, 4.0]]
    b = [7.0, 9.0]
    x = linalg.solve(A, b)
    # 3x+2y=7, x+4y=9 -> x=1, y=2
    assert abs(x[0] - 1.0) < 1e-9 and abs(x[1] - 2.0) < 1e-9


def test_pure_python_fallback_matches_numpy_if_present():
    if not linalg.HAVE_NUMPY:
        return  # nothing to compare against; fallback exercised everywhere else
    import numpy as np
    M = [[5.0, 2.0, 1.0], [2.0, 6.0, 2.0], [1.0, 2.0, 4.0]]
    # force the pure-python path
    saved = linalg.HAVE_NUMPY
    try:
        linalg.HAVE_NUMPY = False
        py_inv = linalg.inv(M)
    finally:
        linalg.HAVE_NUMPY = saved
    np_inv = np.linalg.inv(np.array(M)).tolist()
    assert _close(py_inv, np_inv, 1e-8)


def test_quad_form_diag():
    X = [[1.0, 0.0], [1.0, 1.0], [1.0, 2.0]]
    w = [1.0, 2.0, 0.5]
    R = linalg.quad_form_diag(X, w)
    # X' diag(w) X computed by hand: [[3.5, 3.0],[3.0,4.0]]
    assert abs(R[0][0] - 3.5) < 1e-9
    assert abs(R[0][1] - 3.0) < 1e-9
    assert abs(R[1][1] - 4.0) < 1e-9


if __name__ == "__main__":
    import traceback
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    p = 0
    for fn in fns:
        try:
            fn(); print(f"  PASS {fn.__name__}"); p += 1
        except Exception:
            print(f"  FAIL {fn.__name__}"); traceback.print_exc()
    print(f"\n{p}/{len(fns)} passed")
