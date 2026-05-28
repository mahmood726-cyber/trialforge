"""metaforge.linalg — tiny linear algebra for network meta-analysis.

Uses numpy when available (fast, robust for large networks); otherwise a
pure-Python fallback (Gauss-Jordan inverse + matmul) that is exact for
the small symmetric positive-definite systems NMA produces (typically
5-30 treatments). No install required for the fallback path.
"""
from __future__ import annotations
from typing import List

try:
    import numpy as _np
    HAVE_NUMPY = True
except Exception:  # pragma: no cover
    _np = None
    HAVE_NUMPY = False

Matrix = List[List[float]]
Vector = List[float]


def matmul(A: Matrix, B: Matrix) -> Matrix:
    if HAVE_NUMPY:
        return _np.asarray(A).dot(_np.asarray(B)).tolist()
    n, m, p = len(A), len(B), len(B[0])
    out = [[0.0] * p for _ in range(n)]
    for i in range(n):
        Ai = A[i]
        outi = out[i]
        for k in range(m):
            a = Ai[k]
            if a == 0.0:
                continue
            Bk = B[k]
            for j in range(p):
                outi[j] += a * Bk[j]
    return out


def matvec(A: Matrix, x: Vector) -> Vector:
    if HAVE_NUMPY:
        return _np.asarray(A).dot(_np.asarray(x)).tolist()
    return [sum(a * xi for a, xi in zip(row, x)) for row in A]


def transpose(A: Matrix) -> Matrix:
    return [list(col) for col in zip(*A)]


def inv(A: Matrix) -> Matrix:
    """Inverse of a square matrix."""
    if HAVE_NUMPY:
        return _np.linalg.inv(_np.asarray(A)).tolist()
    n = len(A)
    # Augment [A | I] and Gauss-Jordan with partial pivoting.
    M = [list(A[i]) + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-15:
            raise ValueError("matrix is singular (network may be disconnected)")
        M[col], M[piv] = M[piv], M[col]
        pivval = M[col][col]
        M[col] = [v / pivval for v in M[col]]
        for r in range(n):
            if r == col:
                continue
            factor = M[r][col]
            if factor != 0.0:
                M[r] = [a - factor * b for a, b in zip(M[r], M[col])]
    return [row[n:] for row in M]


def solve(A: Matrix, b: Vector) -> Vector:
    """Solve A x = b."""
    if HAVE_NUMPY:
        return _np.linalg.solve(_np.asarray(A), _np.asarray(b)).tolist()
    return matvec(inv(A), b)


def quad_form_diag(X: Matrix, w: Vector) -> Matrix:
    """Compute X^T diag(w) X efficiently (weighted normal equations LHS)."""
    if HAVE_NUMPY:
        Xn = _np.asarray(X)
        return (Xn.T * _np.asarray(w)).dot(Xn).tolist()
    n, p = len(X), len(X[0])
    out = [[0.0] * p for _ in range(p)]
    for i in range(n):
        wi = w[i]
        Xi = X[i]
        for a in range(p):
            xa = Xi[a] * wi
            if xa == 0.0:
                continue
            outr = out[a]
            for b in range(p):
                outr[b] += xa * Xi[b]
    return out


def xt_w_y(X: Matrix, w: Vector, y: Vector) -> Vector:
    """Compute X^T diag(w) y (weighted normal equations RHS)."""
    if HAVE_NUMPY:
        Xn = _np.asarray(X)
        return (Xn.T * _np.asarray(w)).dot(_np.asarray(y)).tolist()
    p = len(X[0])
    out = [0.0] * p
    for i in range(len(X)):
        wy = w[i] * y[i]
        Xi = X[i]
        for a in range(p):
            out[a] += Xi[a] * wy
    return out
