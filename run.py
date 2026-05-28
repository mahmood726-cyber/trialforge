#!/usr/bin/env python
"""trialforge — AACT-powered advanced meta-analysis from one config.

USAGE
    python run.py configs/your_analysis.json
    python run.py configs/your_analysis.json --out output/x.html

The config selects the base analysis (`type`: pairwise | proportion | nma |
doseresponse), optionally pulls the data straight from an AACT snapshot
(`source: "aact"`), and runs any requested advanced diagnostics
(`advanced: [...]`).

Advanced options (effect-level analyses; binary 2x2 needed for peters/peto/mh):
    egger, peters, trimfill, petpeese, loo, baujat, cumulative,
    subgroup, metareg, peto, mh        (pairwise)
    loops                              (nma — Bucher inconsistency)

Everything runs offline and deterministically (the only optional external
input is a local AACT snapshot). Exit codes: 0 ok, 2 bad config.
"""
from __future__ import annotations
import argparse, json, sys, io
from pathlib import Path

if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from trialforge import (pairwise, proportions, nma, doseresponse, report,  # noqa: E402
                        advanced, nodesplit, tfreport, copas, dta, cnma,
                        limitma, tsa, evalue, pcurve, gosh, cinema,
                        glmm, multivariate, survival, grade, checks)


def die(msg):
    print(f"\n[CONFIG ERROR] {msg}\n", file=sys.stderr)
    print("See README.md or copy a file from configs/.", file=sys.stderr)
    sys.exit(2)


def load_from_aact(cfg):
    """Populate cfg['studies'] + cfg['measure'] from an AACT snapshot."""
    from trialforge import aact
    spec = cfg.get("aact", {})
    try:
        a = aact.AACT(root=spec.get("root"))
    except FileNotFoundError as e:
        die(str(e))
    ncts = spec.get("ncts")
    if not ncts:
        ncts = a.find_trials(drug=spec.get("drug"), condition=spec.get("condition"),
                             limit=spec.get("limit", 500))
        if not ncts:
            die("AACT query returned no trials (check drug/condition spelling; "
                "use specific drug names, not class names).")
    res = a.extract_effects(ncts, prefer_primary=spec.get("prefer_primary", True),
                            force_measure=spec.get("force_measure"))
    if not res["studies"]:
        die(f"AACT found {len(ncts)} trials but no poolable effect estimates "
            f"({res['report'].get('reason','')}). Try force_measure or supply data manually.")
    cfg["studies"] = res["studies"]
    cfg.setdefault("measure", res["measure"])
    return res["report"]


def run_advanced(cfg, effects, measure, ratio):
    """Compute requested advanced diagnostics. Returns dict for tfreport."""
    want = set(cfg.get("advanced", []))
    if not want:
        return {}
    yis = [e.yi for e in effects]
    vis = [e.vi for e in effects]
    names = [e.name for e in effects]
    studies = cfg["studies"]
    adv = {}
    if "egger" in want:
        adv["egger"] = advanced.egger_test(yis, vis)
    if "peters" in want:
        adv["peters"] = advanced.peters_test(studies)
    if "trimfill" in want:
        adv["trimfill"] = advanced.trim_and_fill(yis, vis)
    if "petpeese" in want:
        adv["petpeese"] = advanced.pet_peese(yis, vis)
    if "loo" in want:
        adv["loo"] = advanced.leave_one_out(yis, vis, names)
    if "baujat" in want:
        adv["baujat"] = advanced.baujat(yis, vis, names)
    if "cumulative" in want:
        order = [s.get("year", i) for i, s in enumerate(cfg["studies"])]
        adv["cumulative"] = advanced.cumulative(yis, vis, names, order)
    if "subgroup" in want:
        groups = [s.get("subgroup", "all") for s in cfg["studies"]]
        if len(set(groups)) > 1:
            adv["subgroup"] = advanced.subgroup(yis, vis, groups)
    if "metareg" in want:
        mod = [s.get("moderator") for s in cfg["studies"]]
        if all(m is not None for m in mod):
            adv["metareg"] = advanced.meta_regression(yis, vis, mod, names)
    if "peto" in want:
        adv["peto"] = advanced.peto_or(studies)
    if "mh" in want:
        adv["mh"] = advanced.mantel_haenszel_or(studies)
    if "copas" in want:
        adv["copas"] = copas.analyze(yis, vis, ratio=ratio)
    if "limitma" in want:
        adv["limitma"] = limitma.analyze(yis, vis, ratio=ratio)
    if "tsa" in want:
        adv["tsa"] = tsa.analyze(yis, vis, delta=cfg.get("tsa_delta"), ratio=ratio)
    if "evalue" in want:
        # need the pooled display estimate + CI; compute a quick RE pool
        from trialforge import common as _c
        pool = _c.pool_inverse_variance(yis, vis)
        import math as _m
        disp = (lambda v: _m.exp(v)) if ratio else (lambda v: v)
        adv["evalue"] = evalue.analyze(disp(pool.estimate), disp(pool.ci_low),
                                       disp(pool.ci_high), measure=measure,
                                       rare_outcome=cfg.get("rare_outcome", True))
    if "gosh" in want:
        adv["gosh"] = gosh.analyze(yis, vis, ratio=ratio)
    if "pcurve" in want:
        pv = [s.get("p_value") for s in cfg["studies"]]
        if all(p is not None for p in pv):
            adv["pcurve"] = pcurve.analyze(pv)
    if "glmm" in want:
        if all(all(kk in s for kk in ("tE", "tN", "cE", "cN")) for s in studies):
            adv["glmm"] = glmm.analyze(studies)
    if "grade" in want:
        from trialforge import common as _c
        pool = _c.pool_inverse_variance(yis, vis)
        import math as _m
        dsp = (lambda v: _m.exp(v)) if ratio else (lambda v: v)
        n_total = sum((s.get("tN", 0) + s.get("cN", 0)) for s in studies) or None
        eg = adv.get("egger", {}) or advanced.egger_test(yis, vis)
        adv["grade"] = grade.rate(
            measure=measure, estimate=dsp(pool.estimate),
            ci_low=dsp(pool.ci_low), ci_high=dsp(pool.ci_high),
            k=len(studies), n_total=n_total, i2=pool.i2,
            baseline_risk=cfg.get("baseline_risk"),
            egger_p=eg.get("p") if eg.get("available") else None,
            risk_of_bias=cfg.get("risk_of_bias", "not assessed"),
            indirectness=cfg.get("indirectness", "not assessed"))
    return adv


def main():
    ap = argparse.ArgumentParser(description="trialforge — AACT-powered advanced meta-analysis.")
    ap.add_argument("config")
    ap.add_argument("--out")
    ap.add_argument("--check", action="store_true",
                    help="validate the config + data and report issues without building")
    args = ap.parse_args()

    p = Path(args.config)
    if not p.exists():
        die(f"config not found: {p}")
    try:
        cfg = json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        die(f"config is not valid JSON: {e}")

    source_note = ""
    if (cfg.get("source") or "").lower() == "aact":
        rep = load_from_aact(cfg)
        src = (f"Data auto-extracted from {rep.get('source','AACT')} "
               f"({rep.get('n_with_effects','?')} of {rep.get('n_ncts_queried','?')} "
               f"trials had poolable effects; measure distribution "
               f"{rep.get('measure_distribution','')}").rstrip()
        skipped = rep.get("malformed_rows_skipped")
        if skipped:
            src += f"; {skipped} malformed table row(s) skipped during scan"
        source_note = src + "). Review every value against the source publication before use."

    # Pre-flight data-quality checks (run on the resolved studies).
    issues, csum = checks.check(cfg)
    if args.check or not csum["ok"]:
        for x in issues:
            stream = sys.stderr if x["level"] == "error" else sys.stdout
            print(f"  [{x['level'].upper()}] {x['code']}: {x['msg']}", file=stream)
        print(f"\n  checks: {csum['errors']} error(s), {csum['warnings']} warning(s), "
              f"{csum['info']} info.")
        if args.check:
            sys.exit(0 if csum["ok"] else 2)
        if not csum["ok"]:
            print("\n[BLOCKED] data-quality errors above must be fixed before building.",
                  file=sys.stderr)
            sys.exit(2)

    typ = (cfg.get("type") or "pairwise").lower()
    out = Path(args.out) if args.out else (HERE / "output" / f"{p.stem}.html")
    out.parent.mkdir(parents=True, exist_ok=True)
    tau2_method = cfg.get("tau2_method", "PM")

    # Every analysis type needs a studies array — fail closed with a clear
    # message (not a raw traceback) for all types, not just pairwise.
    if not cfg.get("studies") and not cfg.get("trials"):
        die(f"'{typ}' needs a 'studies' array in the config.")

    if typ == "pairwise":
        measure = (cfg.get("measure") or "OR").upper()
        trials = cfg.get("studies") or cfg.get("trials")
        if not trials:
            die("pairwise needs 'studies'.")
        pool, effects, skip = pairwise.analyze(trials, measure, tau2_method)
        if pool is None:
            die("no usable studies.")
        base = report.render_pairwise(cfg, pool, measure)
        adv = run_advanced(cfg, effects, measure, pool.extra["ratio"])
        adv_html = tfreport.advanced_section(adv, ratio=pool.extra["ratio"])
        html_doc = tfreport.splice(base, adv_html)
        summary = (f"{pool.k} studies · {measure} {pool.extra['display']['estimate']:.3f} "
                   f"· I2={pool.i2:.0f}% · advanced: {','.join(adv) or 'none'}")

    elif typ == "proportion":
        pool, ps, skip = proportions.analyze(cfg["studies"], cfg.get("method", "DAS"), tau2_method)
        if pool is None:
            die("no usable studies.")
        html_doc = report.render_proportion(cfg, pool)
        d = pool.extra["display"]
        summary = f"{pool.k} studies · proportion {d['estimate']*100:.1f}%"

    elif typ == "nma":
        res = nma.analyze(cfg["studies"], measure=(cfg.get("measure") or "OR").upper(),
                          reference=cfg.get("reference"),
                          smaller_better=cfg.get("smaller_better", True))
        if res is None:
            die("NMA failed (need >=2 connected treatments).")
        base = report.render_nma(cfg, res)
        adv = {}
        want_nma = set(cfg.get("advanced", []))
        if "loops" in want_nma:
            adv["loops"] = nodesplit.loop_inconsistency(
                cfg["studies"], (cfg.get("measure") or "OR").upper())
        if "cinema" in want_nma:
            adv["cinema"] = cinema.rate(
                res, loops=adv.get("loops"), ratio=res["ratio"],
                judgments=cfg.get("cinema_judgments"))
        adv_html = tfreport.advanced_section(adv, ratio=res["ratio"])
        html_doc = tfreport.splice(base, adv_html)
        best = res["ranking"][0]
        summary = (f"{res['n_studies']} studies · {len(res['treatments'])} treatments · "
                   f"top {best['treatment']} ({best['sucra']:.0f}%)"
                   + (f" · loops:{adv['loops']['n_loops']}" if 'loops' in adv else ""))

    elif typ == "doseresponse":
        pool, slopes, skip = doseresponse.analyze(
            cfg["studies"], measure=(cfg.get("measure") or "RR").upper(),
            tau2_method=tau2_method, predict_doses=cfg.get("predict_doses"))
        if pool is None:
            die("no usable studies (need >=2 dose levels each).")
        html_doc = report.render_doseresponse(cfg, pool)
        summary = f"{pool.k} studies · {pool.extra['measure']} per unit {pool.extra['slope_display_per_unit']:.3f}"

    elif typ == "dta":
        res = dta.analyze(cfg["studies"])
        if not res.get("available"):
            die(f"DTA failed: {res.get('reason','need 2x2 tables tp/fp/fn/tn')}.")
        html_doc = tfreport.render_dta(cfg, res)
        summary = (f"{res['k']} studies · Se {res['sensitivity']*100:.0f}% "
                   f"Sp {res['specificity']*100:.0f}% · DOR {res['dor']:.1f}")

    elif typ == "cnma":
        res = cnma.analyze(cfg["studies"], measure=(cfg.get("measure") or "OR").upper())
        if not res.get("available"):
            die(f"component-NMA failed: {res.get('reason','')}.")
        html_doc = tfreport.render_cnma(cfg, res)
        summary = f"{len(res['components'])} components · {res['n_contrasts']} contrasts"

    elif typ == "multivariate":
        res = multivariate.analyze(cfg["studies"],
                                   ratio=cfg.get("ratio", False))
        if not res.get("available"):
            die(f"multivariate failed: {res.get('reason','')}.")
        html_doc = tfreport.render_multivariate(cfg, res)
        o1 = res["outcome1"]
        summary = (f"{res['k']} studies · outcome1 {o1['bivariate']:.3f}; "
                   f"borrowed precision {o1['borrowed_precision_pct']:.0f}%")

    elif typ == "rmst":
        res = survival.analyze(cfg["studies"], tau=cfg.get("tau_star"),
                               tau2_method=tau2_method)
        if not res.get("available"):
            die(f"RMST failed: {res.get('reason','')}.")
        html_doc = tfreport.render_rmst(cfg, res)
        summary = (f"{res['k']} studies · RMST difference {res['rmst_difference']:.2f} "
                   f"(95% CI {res['ci_low']:.2f} to {res['ci_high']:.2f})")

    else:
        die(f"unknown type '{typ}'. Use pairwise, proportion, nma, "
            "doseresponse, dta, cnma, multivariate, rmst.")

    # Data-checks section in every report (warnings are advisory; errors
    # already blocked above).
    html_doc = tfreport.splice(html_doc, tfreport.checks_section(issues))

    # Ensure the AACT "verify against source" caveat is in the rendered HTML
    # for EVERY type (it previously only appeared when an advanced-diagnostics
    # section happened to be present).
    if source_note:
        aact_section = (
            '<section><div class="disclaimer">'
            '<strong>AACT-assisted extraction.</strong> ' + source_note
            + '</div></section>')
        html_doc = tfreport.splice(html_doc, aact_section)

    out.write_text(html_doc, encoding="utf-8")
    if csum["warnings"]:
        print(f"  ({csum['warnings']} data-quality warning(s) — see the report's "
              "'Data checks' section)")
    print(f"\nBuilt: {out}")
    print(f"  {summary}")
    if source_note:
        print(f"  (AACT-assisted extraction — review values before use)")
    print(f"\nOpen the file in any browser (works offline).")


if __name__ == "__main__":
    main()
