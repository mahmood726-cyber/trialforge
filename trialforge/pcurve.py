"""trialforge.pcurve — p-curve analysis (Simonsohn 2014) for evidential value.

Ported (pure-stdlib) from the allmeta `p-curve` engine. Examines the
distribution of statistically significant p-values (< 0.05):
  * a right-skewed p-curve (many very-small p-values) indicates the studies
    contain evidential value;
  * a flat / left-skewed curve suggests no evidential value or p-hacking.

Tests implemented:
  * Fisher's combined test on pp = p/0.05:  T = -2 sum ln(pp) ~ chi2(2k);
    a small fisher_p => right-skew => evidential value.
  * Flatness (binomial) test: proportion of pp <= 0.5 vs the 0.5 expected
    under a flat curve.

Validated against the allmeta R fixture (pchisq/pnorm): 5 significant
p-values -> fisher_chisq 10.8887 (df 10), fisher_p 0.3663, prop_low 0.6,
z_flat 0.4472, p_flat 0.6547.
"""
from __future__ import annotations
import math
from . import common


def analyze(p_values, alpha=0.05):
    sig = [p for p in p_values if 0 < p < alpha]
    k = len(sig)
    if k < 2:
        return {"available": False, "reason": "need >=2 significant (p<0.05) studies",
                "n_significant": k}
    pp = [p / alpha for p in sig]
    # Fisher combined on pp
    T = -2.0 * sum(math.log(x) for x in pp)
    df = 2 * k
    fisher_p = common.chi2_sf(T, df)
    # Binomial flatness: proportion of pp <= 0.5
    n_low = sum(1 for x in pp if x <= 0.5)
    prop_low = n_low / k
    z_flat = (prop_low - 0.5) / math.sqrt(0.25 / k)
    p_flat = 2 * (1 - common.norm_cdf(abs(z_flat)))
    return {
        "available": True,
        "n_significant": k,
        "fisher_chisq": T, "fisher_df": df, "fisher_p": fisher_p,
        "prop_low": prop_low, "z_flat": z_flat, "p_flat": p_flat,
        "evidential_value": fisher_p < 0.05,
        "interpretation": (
            "right-skewed p-curve: the studies contain evidential value"
            if fisher_p < 0.05 else
            "p-curve not significantly right-skewed: limited evidence of "
            "evidential value (could reflect low power, null effects, or "
            "selective reporting)"),
        "note": "p-curve uses only statistically significant results; it "
                "assesses evidential value, not effect size.",
    }
