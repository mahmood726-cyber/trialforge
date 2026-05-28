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


def _num(name, field, val, out):
    """Coerce a field to float; emit a bad_type error and return None on failure.
    Keeps check() to its 'does not raise' contract for sloppy JSON."""
    if val is None:
        return None
    if isinstance(val, bool):  # bool is an int subclass — reject explicitly
        out.append(_err("bad_type", f"{name}: '{field}' is a boolean, expected a number."))
        return None
    if isinstance(val, (int, float)):
        import math as _m
        if not _m.isfinite(val):
            out.append(_err("nonfinite_value", f"{name}: '{field}' is not finite ({val})."))
            return None
        return float(val)
    out.append(_err("bad_type", f"{name}: '{field}'={val!r} is not a number."))
    return None


def _check_binary_row(name, t, out):
    if not all(k in t for k in ("tE", "tN", "cE", "cN")):
        return
    tE = _num(name, "tE", t["tE"], out); tN = _num(name, "tN", t["tN"], out)
    cE = _num(name, "cE", t["cE"], out); cN = _num(name, "cN", t["cN"], out)
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
    if not all(k in t for k in ("effect", "ci_low", "ci_high")):
        return
    eff = _num(name, "effect", t["effect"], out)
    lo = _num(name, "ci_low", t["ci_low"], out)
    hi = _num(name, "ci_high", t["ci_high"], out)
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


def _check_proportion_row(name, t, out):
    e = _num(name, "e", t.get("e"), out)
    n = _num(name, "n", t.get("n"), out)
    if None in (e, n):
        return
    if n <= 0:
        out.append(_err("zero_n", f"{name}: n={n} (must be >0)."))
    elif e < 0 or e > n:
        out.append(_err("impossible_count",
                        f"{name}: events {e} outside [0, {n}]."))


def _check_dta_row(name, t, out):
    vals = {k: _num(name, k, t.get(k), out) for k in ("tp", "fp", "fn", "tn")}
    for k, v in vals.items():
        if v is not None and v < 0:
            out.append(_err("impossible_count", f"{name}: {k}={v} is negative."))
    if all(v is not None for v in vals.values()):
        if vals["tp"] + vals["fn"] == 0 or vals["fp"] + vals["tn"] == 0:
            out.append(_err("empty_margin",
                            f"{name}: a diseased/non-diseased margin is zero."))


def _check_se_field(name, t, key, out):
    se = _num(name, key, t.get(key), out)
    if se is not None and se <= 0:
        out.append(_err("nonpositive_se", f"{name}: {key}={se} must be > 0."))


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

    # per-row data checks, dispatched by analysis type
    for i, s in enumerate(studies):
        name = s.get("name") or s.get("nct") or f"(study {i+1})"
        if typ in ("nma", "cnma"):
            for arm in s.get("arms", []):
                if all(k in arm for k in ("e", "n")):
                    e = _num(name, "e", arm["e"], out)
                    n = _num(name, "n", arm["n"], out)
                    if None not in (e, n) and (n <= 0 or e < 0 or e > n):
                        out.append(_err("impossible_count",
                                        f"{name}: arm {arm.get('t','?')} has "
                                        f"events {e} / N {n} out of range."))
        elif typ == "proportion":
            _check_proportion_row(name, s, out)
        elif typ == "dta":
            _check_dta_row(name, s, out)
        elif typ == "multivariate":
            for key in ("se1", "se2", "se"):
                if key in s:
                    _check_se_field(name, s, key, out)
            for key in ("y1", "y2", "y", "r"):
                if key in s:
                    _num(name, key, s[key], out)
        elif typ in ("rmst", "survival"):
            for key in ("rmst_diff", "ci_low", "ci_high"):
                if key in s:
                    _num(name, key, s[key], out)
            if "se" in s and s["se"] is not None:
                _check_se_field(name, s, "se", out)
        else:  # pairwise, doseresponse, and anything precomputed
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
