"""Cross-validation against the E156 flagship 'SGLT2 inhibitors in heart
failure' living capsule (F:/E156/flagship/sglt2-hf-capsule.html).

The capsule pools the published HR + 95% CI per trial with Paule-Mandel
tau^2 and a Knapp-Hartung CI with floor max(1, Q_gen/df) — the same
methodology as trialforge. This test replicates the capsule's pool
arithmetic in Python and confirms trialforge reproduces its pooled HR.
"""
import math, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from trialforge import pairwise, common

# The capsule's PUBLISHED trial table (sglt2-hf-capsule.html, lines 322-326)
TRIALS = [
    ("DAPA-HF", 0.74, 0.65, 0.85),
    ("EMPEROR-Reduced", 0.75, 0.65, 0.86),
    ("EMPEROR-Preserved", 0.79, 0.69, 0.90),
    ("DELIVER", 0.82, 0.73, 0.92),
    ("SOLOIST-WHF", 0.67, 0.52, 0.85),
]
Z = 1.959964


def _capsule_pool():
    """Replicate the capsule pool(): PM tau^2, weighted mean of log-HR."""
    y = [math.log(hr) for _, hr, lo, hi in TRIALS]
    v = [((math.log(hi) - math.log(lo)) / (2 * Z)) ** 2 for _, hr, lo, hi in TRIALS]
    tau2 = common.tau2_paule_mandel(y, v)
    w = [1.0 / (vi + tau2) for vi in v]
    sw = sum(w)
    re = sum(wi * yi for wi, yi in zip(w, y)) / sw
    # I^2 at fixed-effect weights
    wf = [1.0 / vi for vi in v]
    mf = sum(wi * yi for wi, yi in zip(wf, y)) / sum(wf)
    Q = sum(wi * (yi - mf) ** 2 for wi, yi in zip(wf, y))
    df = len(y) - 1
    i2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    return math.exp(re), i2, tau2


def test_trialforge_reproduces_capsule_pool():
    cap_hr, cap_i2, cap_tau2 = _capsule_pool()
    studies = [{"name": n, "effect": hr, "ci_low": lo, "ci_high": hi}
               for n, hr, lo, hi in TRIALS]
    pool, eff, skip = pairwise.analyze(studies, "HR", tau2_method="PM")
    tf_hr = pool.extra["display"]["estimate"]
    # point estimate (PM) must match the capsule to high precision
    assert math.isclose(tf_hr, cap_hr, abs_tol=1e-4), f"{tf_hr} vs capsule {cap_hr}"
    assert math.isclose(pool.i2, cap_i2, abs_tol=1e-6)
    assert math.isclose(pool.tau2, cap_tau2, abs_tol=1e-8)
    # sanity: pooled HR is a meaningful benefit around 0.77
    assert 0.74 < tf_hr < 0.80


def test_flagship_direction_and_significance():
    studies = [{"name": n, "effect": hr, "ci_low": lo, "ci_high": hi}
               for n, hr, lo, hi in TRIALS]
    pool, _, _ = pairwise.analyze(studies, "HR", tau2_method="PM")
    d = pool.extra["display"]
    assert d["ci_high"] < 1.0   # CI excludes the null -> significant benefit


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
