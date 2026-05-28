"""trialforge.checks — data-quality + config sanity checks.

Runs BEFORE the analysis (and is also available as a standalone `--check`
preflight). Surfaces problems that would otherwise produce silently-wrong
or misleading results: impossible event counts, inverted/again-null CIs,
non-positive values on ratio scales, too-few studies for the chosen
diagnostics, duplicate / malformed identifiers, etc.

Each issue is {level, code, msg} with level in {error, warning, info}.
Errors should block; warnings are surfaced in the report's "Data checks"
section and on the console.
"""
from __future__ import annotations
import re

_NCT_RE = re.compile(r"^NCT\d{7,8}$")
RATIO = {"OR", "RR", "HR"}


def _err(code, msg):
    return {"level": "error", "code": code, "msg": msg}


def _warn(code, msg):
    return {"level": "warning", "code": code, "msg": msg}


def _info(code, msg):
    return {"level": "info", "code": code, "msg": msg}


def _check_binary_row(name, t, out):
    tE, tN, cE, cN = t.get("tE"), t.get("tN"), t.get("cE"), t.get("cN")
    if None in (tE, tN, cE, cN):
        return
    for lbl, e, n in (("intervention", tE, tN), ("comparator", cE, cN)):
        if n <= 0:
            out.append(_err("zero_n", f"{name}: {lbl} arm has N={n} (must be >0)."))
        elif e < 0 or e > n:
            out.append(_err("impossible_count",
                            f"{name}: {lbl} events {e} outside [0, {n}] — "
                            "events cannot exceed sample size."))
    if 0 in (tE, cE, tN - (tE or 0), cN - (cE or 0)):
        out.append(_info("zero_cell",
                         f"{name}: a zero cell is present — a 0.5 continuity "
                         "correction is applied (or use type-aware rare-event "
                         "methods: glmm/peto/mh)."))


def _check_precomputed_row(name, t, ratio, out):
    eff, lo, hi = t.get("effect"), t.get("ci_low"), t.get("ci_high")
    if None in (eff, lo, hi):
        return
    if lo > hi:
        out.append(_err("inverted_ci",
                        f"{name}: ci_low ({lo}) > ci_high ({hi})."))
    if not (lo <= eff <= hi):
        out.append(_warn("estimate_outside_ci",
                         f"{name}: effect {eff} lies outside its CI [{lo}, {hi}]."))
    if ratio and min(eff, lo, hi) <= 0:
        out.append(_err("nonpositive_ratio",
                        f"{name}: ratio-scale value <= 0 (effect/CI must be >0 "
                        "for OR/RR/HR)."))


def check(cfg):
    """Return (issues, summary) for a config. Does not raise."""
    out = []
    typ = (cfg.get("type") or "pairwise").lower()
    measure = (cfg.get("measure") or ("RR" if typ == "doseresponse" else "OR")).upper()
    ratio = measure in RATIO
    studies = cfg.get("studies") or cfg.get("trials") or []

    if not studies:
        out.append(_err("no_studies", f"'{typ}' needs a non-empty 'studies' array."))
        return out, _summary(out)

    # duplicate / missing names + NCT format
    seen = {}
    for i, s in enumerate(studies):
        name = s.get("name") or s.get("nct") or f"(study {i+1})"
        seen[name] = seen.get(name, 0) + 1
        nct = s.get("nct")
        if nct and not _NCT_RE.match(str(nct)):
            out.append(_warn("nct_format",
                             f"{name}: '{nct}' is not a valid NCT id (NCT + 7-8 digits)."))
    for name, n in seen.items():
        if n > 1:
            out.append(_warn("duplicate_name", f"study name '{name}' appears {n} times."))

    # per-row data checks (pairwise-like)
    if typ in ("pairwise", "doseresponse", "cnma", "nma"):
        for i, s in enumerate(studies):
            name = s.get("name") or s.get("nct") or f"(study {i+1})"
            if typ in ("nma", "cnma"):
                for arm in s.get("arms", []):
                    if all(k in arm for k in ("e", "n")):
                        if arm["n"] <= 0 or arm["e"] < 0 or arm["e"] > arm["n"]:
                            out.append(_err("impossible_count",
                                            f"{name}: arm {arm.get('t','?')} has "
                                            f"events {arm['e']} / N {arm['n']} out of range."))
            else:
                _check_binary_row(name, s, out)
                _check_precomputed_row(name, s, ratio, out)

    # k-based methodological warnings
    k = len(studies)
    if typ in ("pairwise", "proportion", "doseresponse") and k < 2:
        out.append(_err("too_few", f"only {k} study — need >=2 to pool."))
    if typ == "pairwise":
        if k < 3:
            out.append(_warn("pi_undefined",
                             f"k={k}: the prediction interval is undefined (needs >=3)."))
        adv = set(cfg.get("advanced", []))
        if adv & {"egger", "peters", "petpeese", "copas", "limitma"} and k < 10:
            out.append(_warn("pubbias_low_power",
                             f"k={k}: publication-bias tests have low power and are "
                             "unreliable below ~10 studies (advanced-stats.md)."))
        if "copas" in adv and k < 15:
            out.append(_info("copas_k",
                             f"k={k}: Copas needs >=15 studies for stable MLE; the "
                             "HT sensitivity profile is directional only."))

    # GRADE needs a baseline risk for absolute effects
    if "grade" in set(cfg.get("advanced", [])) and cfg.get("baseline_risk") is None:
        out.append(_info("grade_no_baseline",
                         "advanced 'grade' without 'baseline_risk' — the GRADE "
                         "Summary-of-Findings absolute effects will be omitted."))

    return out, _summary(out)


def _summary(out):
    n_err = sum(1 for x in out if x["level"] == "error")
    n_warn = sum(1 for x in out if x["level"] == "warning")
    n_info = sum(1 for x in out if x["level"] == "info")
    return {"errors": n_err, "warnings": n_warn, "info": n_info, "ok": n_err == 0}
