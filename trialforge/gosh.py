"""trialforge.gosh — GOSH (Graphical Of Subsets Heterogeneity) diagnostic.

Ported (pure-stdlib) from the allmeta `gosh` engine. Pools every subset of
studies (size >= 2) and summarises the distribution of the resulting
estimates and I^2. Multimodality or wide spread reveals influential studies
or distinct study clusters driving heterogeneity.

  * full enumeration for k <= 15 (all 2^k - 1 - k subsets of size >= 2)
  * random subset sampling for k > 15 (advanced-stats.md rule)
  * FE (inverse-variance) or DL random-effects pooling per subset

Validated against metafor::gosh: 5-study fixture -> 26 subsets, full-k FE
estimate 0.22542, median subset estimate 0.22173.
"""
from __future__ import annotations
import math
import random
from itertools import combinations
from . import common


def _pool(idx, yis, vis, method):
    if method == "DL":
        ys = [yis[i] for i in idx]
        vs = [vis[i] for i in idx]
        tau2 = common.tau2_dersimonian_laird(ys, vs)
        w = [1.0 / (v + tau2) for v in vs]
    else:  # FE
        w = [1.0 / vis[i] for i in idx]
    sw = sum(w)
    est = sum(wi * yis[i] for wi, i in zip(w, idx)) / sw
    # I^2 (always FE-weight Q)
    wf = [1.0 / vis[i] for i in idx]
    muf = sum(wi * yis[i] for wi, i in zip(wf, idx)) / sum(wf)
    Q = sum(wi * (yis[i] - muf) ** 2 for wi, i in zip(wf, idx))
    dfree = len(idx) - 1
    i2 = max(0.0, (Q - dfree) / Q) * 100 if Q > 0 else 0.0
    return est, i2


def analyze(yis, vis, *, method="FE", ratio=False, max_full=15,
            n_samples=5000, seed=20260528):
    k = len(yis)
    if k < 3:
        return {"available": False, "reason": "need >=3 studies"}
    idx_all = list(range(k))
    subsets = []
    if k <= max_full:
        for size in range(2, k + 1):
            subsets.extend(combinations(idx_all, size))
        sampled = False
    else:
        rng = random.Random(seed)
        seen = set()
        while len(subsets) < n_samples:
            size = rng.randint(2, k)
            sub = tuple(sorted(rng.sample(idx_all, size)))
            if sub not in seen:
                seen.add(sub)
                subsets.append(sub)
        sampled = True

    ests, i2s = [], []
    for sub in subsets:
        e, i2 = _pool(sub, yis, vis, method)
        ests.append(e)
        i2s.append(i2)

    def disp(v):
        return math.exp(v) if ratio else v

    ests_sorted = sorted(ests)
    n = len(ests_sorted)

    def q(p):
        # type-7 linear interpolation (R quantile default)
        if n == 1:
            return ests_sorted[0]
        h = (n - 1) * p
        lo = int(math.floor(h))
        hi = min(n - 1, lo + 1)
        return ests_sorted[lo] + (h - lo) * (ests_sorted[hi] - ests_sorted[lo])

    full_est, full_i2 = _pool(idx_all, yis, vis, method)
    return {
        "available": True, "k": k, "method": method,
        "n_subsets": n, "sampled": sampled,
        "median_estimate": disp(q(0.5)),
        "q25_estimate": disp(q(0.25)), "q75_estimate": disp(q(0.75)),
        "min_estimate": disp(min(ests)), "max_estimate": disp(max(ests)),
        "median_i2": sorted(i2s)[len(i2s) // 2],
        "max_i2": max(i2s),
        "full_estimate": disp(full_est), "full_i2": full_i2,
        "note": "Wide spread or multimodality across subsets indicates "
                "influential studies or distinct clusters.",
    }
