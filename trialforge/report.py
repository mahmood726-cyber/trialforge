"""metaforge.report — assemble a standalone HTML report for any analysis type."""
from __future__ import annotations
import html
import math
from . import plots

MEASURE_NAME = {
    "OR": "Odds ratio", "RR": "Risk ratio", "RD": "Risk difference",
    "HR": "Hazard ratio", "MD": "Mean difference",
    "SMD": "Standardised mean difference (Hedges g)",
}


def _h(s):
    return "&mdash;" if s is None else html.escape(str(s))


def _n(v, nd=2):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "&mdash;"
    return f"{v:.{nd}f}"


def _p(v):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "&mdash;"
    return f"{v:.3f}" if v >= 0.001 else f"{v:.2e}"


_CSS = """
:root{--ink:#1e293b;--muted:#64748b;--line:#e2e8f0;--accent:#2563eb;--bg:#f8fafc;}
*{box-sizing:border-box;} body{margin:0;background:var(--bg);color:var(--ink);
font-family:system-ui,-apple-system,sans-serif;line-height:1.55;}
header{background:#0f172a;color:#fff;padding:26px 24px;}
header h1{margin:0 0 6px;font-size:21px;} header .sub{color:#cbd5e1;font-size:13px;}
header .badge{display:inline-block;background:#2563eb;color:#fff;font-size:11px;
padding:2px 9px;border-radius:10px;margin-right:6px;text-transform:uppercase;letter-spacing:.04em;}
main{max-width:900px;margin:0 auto;padding:24px;}
section{background:#fff;border:1px solid var(--line);border-radius:8px;padding:20px 22px;margin-bottom:18px;}
h2{font-size:15px;margin:0 0 12px;color:#0f172a;}
.headline{font-size:16px;line-height:1.6;}
.kpis{display:flex;gap:12px;flex-wrap:wrap;margin-top:8px;}
.kpi{background:var(--bg);border:1px solid var(--line);border-radius:6px;padding:10px 14px;min-width:110px;}
.kpi .v{font-size:19px;font-weight:600;color:var(--accent);} .kpi .l{font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:.04em;}
table{width:100%;border-collapse:collapse;font-size:13px;} th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--line);}
th{background:var(--bg);font-size:12px;} td.num,th.num{text-align:right;font-variant-numeric:tabular-nums;}
.fig{overflow-x:auto;} .fig svg{max-width:100%;height:auto;border:1px solid var(--line);border-radius:6px;}
.disclaimer{background:#fffbeb;border:1px solid #fde68a;color:#92400e;border-radius:6px;padding:12px 16px;font-size:12px;}
footer{max-width:900px;margin:0 auto;padding:0 24px 32px;color:var(--muted);font-size:11px;}
code{background:var(--bg);padding:1px 5px;border-radius:3px;}
"""


def _shell(title, badge, sub, body):
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{_h(title)}</title><meta name="description" content="{_h(title)} — metaforge meta-analysis">
<style>{_CSS}</style></head><body>
<header><h1>{_h(title)}</h1><div class="sub"><span class="badge">{_h(badge)}</span>{sub}</div></header>
<main>{body}
<section><div class="disclaimer"><strong>Auto-generated draft.</strong> metaforge computes the
statistics from the data you supplied. Verify every number against source publications and confirm
the studies are clinically combinable before using or citing any result.</div></section>
</main>
<footer>Generated offline by <code>metaforge</code> · random-effects (Paule&ndash;Mandel &tau;&sup2;,
Knapp&ndash;Hartung CI, Cochrane v6.5 prediction interval)</footer></body></html>"""


def _het_table(pool):
    return f"""<table>
<tr><td>&tau;&sup2; (between-study variance)</td><td class="num">{_n(pool.tau2,4)}</td></tr>
<tr><td>I&sup2;</td><td class="num">{_n(pool.i2,1)}%</td></tr>
<tr><td>Cochran Q (df = {pool.Q_df})</td><td class="num">{_n(pool.Q,2)}, p = {_p(pool.Q_p)}</td></tr>
<tr><td>&tau;&sup2; estimator</td><td class="num">{_h(pool.tau2_method)}</td></tr>
<tr><td>Test of overall effect</td><td class="num">z = {_n(pool.z,2)}, p = {_n(pool.p,4)}</td></tr>
</table>"""


# ---------------------------------------------------------------------------
def render_pairwise(cfg, pool, measure):
    d = pool.extra["display"]
    ratio = pool.extra["ratio"]
    mname = MEASURE_NAME.get(measure, measure)
    null = 1.0 if ratio else 0.0
    crosses = d["ci_low"] <= null <= d["ci_high"]
    direction = ("no statistically significant difference" if crosses else
                 "a statistically significant benefit of the intervention"
                 if ((ratio and d["estimate"] < 1) or (not ratio and d["estimate"] < 0))
                 else "a statistically significant effect favouring the comparator")
    het = ("low" if pool.i2 < 30 else "moderate" if pool.i2 < 60 else
           "substantial" if pool.i2 < 75 else "considerable")
    rows = [{"name": s["name"], "est": s["est"], "lo": s["lo"], "hi": s["hi"],
             "weight": s["weight"]} for s in pool.extra["per_study"]]
    forest = plots.forest_svg(
        rows, {"est": d["estimate"], "lo": d["ci_low"], "hi": d["ci_high"]},
        ratio=ratio, title=f"{measure} (95% CI)",
        pi=(d["pi_low"], d["pi_high"]) if d["pi_low"] is not None else None,
        null_label_left=f"favours {cfg.get('intervention','intervention')}",
        null_label_right=f"favours {cfg.get('comparator','comparator')}")
    pi_txt = (f"{_n(d['pi_low'])} to {_n(d['pi_high'])}" if d["pi_low"] is not None
              else "undefined (k&lt;3)")
    body = f"""
<section><h2>Result</h2>
<p class="headline">Across {pool.k} studies the random-effects pool shows {direction}
({mname} {_n(d['estimate'])}, 95% CI {_n(d['ci_low'])} to {_n(d['ci_high'])}).
Heterogeneity is {het} (I&sup2; {_n(pool.i2,1)}%).</p>
<div class="kpis">
<div class="kpi"><div class="v">{_n(d['estimate'])}</div><div class="l">{mname}</div></div>
<div class="kpi"><div class="v">{_n(d['ci_low'])}&ndash;{_n(d['ci_high'])}</div><div class="l">95% CI</div></div>
<div class="kpi"><div class="v">{pool.k}</div><div class="l">Studies</div></div>
<div class="kpi"><div class="v">{_n(pool.i2,1)}%</div><div class="l">I&sup2;</div></div>
<div class="kpi"><div class="v">{pi_txt}</div><div class="l">95% prediction interval</div></div>
</div></section>
<section><h2>Forest plot</h2><div class="fig">{forest}</div></section>
<section><h2>Methods &amp; heterogeneity</h2>{_het_table(pool)}</section>"""
    sub = f"{_h(cfg.get('intervention',''))} vs {_h(cfg.get('comparator',''))} &middot; {mname}"
    return _shell(cfg.get("title", "Pairwise meta-analysis"), "pairwise", sub, body)


def render_proportion(cfg, pool):
    d = pool.extra["display"]
    method = pool.extra["method"]
    rows = [{"name": s["name"], "est": s["est"], "lo": s["lo"], "hi": s["hi"],
             "weight": s["weight"]} for s in pool.extra["per_study"]]
    forest = plots.forest_svg(
        rows, {"est": d["estimate"], "lo": d["ci_low"], "hi": d["ci_high"]},
        ratio=False, title="Proportion (95% CI)",
        pi=(d["pi_low"], d["pi_high"]) if d["pi_low"] is not None else None,
        null_label_left="lower", null_label_right="higher")
    pct = lambda v: "&mdash;" if v is None else f"{v*100:.1f}%"
    pi_txt = (f"{pct(d['pi_low'])} to {pct(d['pi_high'])}" if d["pi_low"] is not None
              else "undefined (k&lt;3)")
    body = f"""
<section><h2>Result</h2>
<p class="headline">The pooled proportion across {pool.k} studies is
<strong>{pct(d['estimate'])}</strong> (95% CI {pct(d['ci_low'])} to {pct(d['ci_high'])}),
by the {_h(method)} random-effects model. Heterogeneity I&sup2; {_n(pool.i2,1)}%.</p>
<div class="kpis">
<div class="kpi"><div class="v">{pct(d['estimate'])}</div><div class="l">Pooled proportion</div></div>
<div class="kpi"><div class="v">{pct(d['ci_low'])}&ndash;{pct(d['ci_high'])}</div><div class="l">95% CI</div></div>
<div class="kpi"><div class="v">{pool.k}</div><div class="l">Studies</div></div>
<div class="kpi"><div class="v">{_n(pool.i2,1)}%</div><div class="l">I&sup2;</div></div>
<div class="kpi"><div class="v">{pi_txt}</div><div class="l">95% prediction interval</div></div>
</div></section>
<section><h2>Forest plot</h2><div class="fig">{forest}</div></section>
<section><h2>Methods &amp; heterogeneity</h2>{_het_table(pool)}
<p style="font-size:12px;color:#64748b">Transformation: {_h(method)}
({'double-arcsine, back-transformed with the harmonic mean of n' if method=='DAS' else method}).</p></section>"""
    return _shell(cfg.get("title", "Proportion meta-analysis"), "proportion",
                  _h(cfg.get("outcome", "")), body)


def render_nma(cfg, res):
    ratio = res["ratio"]
    treatments = res["treatments"]
    # network geometry
    edge_counts, study_counts = {}, {}
    for s in cfg.get("studies", []):
        ts = [a["t"] for a in s["arms"]]
        for t in ts:
            study_counts[t] = study_counts.get(t, 0) + 1
        for i in range(len(ts)):
            for j in range(i + 1, len(ts)):
                key = tuple(sorted((ts[i], ts[j])))
                edge_counts[key] = edge_counts.get(key, 0) + 1
    net = plots.network_svg(treatments, edge_counts, study_counts)
    sucra = plots.sucra_svg(res["ranking"])
    # vs-reference forest
    rows = [{"name": r["treatment"], "est": r["estimate"], "lo": r["ci_low"],
             "hi": r["ci_high"], "weight": None} for r in res["vs_reference"]]
    forest = plots.forest_svg(
        rows, {"est": 1.0 if ratio else 0.0, "lo": 1.0 if ratio else 0.0,
               "hi": 1.0 if ratio else 0.0},
        ratio=ratio, title=f"vs {res['reference']} ({res['measure']}, 95% CI)",
        null_label_left=f"favours treatment", null_label_right=f"favours {res['reference']}") \
        if rows else "<p>(reference only)</p>"
    rank_rows = "\n".join(
        f"<tr><td>{_h(r['treatment'])}</td><td class='num'>{_n(r['sucra'],0)}%</td>"
        f"<td class='num'>{_n(r['mean_rank'],1)}</td></tr>" for r in res["ranking"])
    pw_rows = "\n".join(
        f"<tr><td>{_h(p['a'])} vs {_h(p['b'])}</td>"
        f"<td class='num'>{_n(p['estimate'])} ({_n(p['ci_low'])}, {_n(p['ci_high'])})</td></tr>"
        for p in res["pairwise"] if p["a"] < p["b"])
    body = f"""
<section><h2>Network</h2>
<p class="headline">{res['n_studies']} studies, {len(treatments)} treatments,
{res['n_contrasts']} contrasts. Reference: <strong>{_h(res['reference'])}</strong>.
Between-study SD &tau; = {_n(res['tau'],3)}; network heterogeneity/inconsistency
Q = {_n(res['Q'],2)} (df {res['Q_df']}, p = {_p(res['Q_p'])}).</p>
<div class="fig">{net}</div></section>
<section><h2>Treatment ranking (SUCRA)</h2>
<p style="font-size:12px;color:#64748b">{'Lower' if res['smaller_better'] else 'Higher'} effect is treated as better.</p>
<div class="fig">{sucra}</div>
<table><thead><tr><th>Treatment</th><th class="num">SUCRA</th><th class="num">Mean rank</th></tr></thead>
<tbody>{rank_rows}</tbody></table></section>
<section><h2>Relative effects vs {_h(res['reference'])}</h2><div class="fig">{forest}</div></section>
<section><h2>All pairwise relative effects</h2>
<table><thead><tr><th>Comparison</th><th class="num">{res['measure']} (95% CI)</th></tr></thead>
<tbody>{pw_rows}</tbody></table>
<p style="font-size:12px;color:#64748b">Frequentist contrast-based random-effects NMA
(consistency model). Formal node-splitting inconsistency testing is not included;
the network Q above is the global heterogeneity/inconsistency indicator.</p></section>"""
    return _shell(cfg.get("title", "Network meta-analysis"), "network MA",
                  f"{len(treatments)} treatments &middot; {res['measure']}", body)


def render_doseresponse(cfg, pool):
    e = pool.extra
    ratio = e["ratio"]
    curve = plots.dose_svg(
        e["curve"], ratio=ratio,
        xlabel=cfg.get("dose_label", "Dose increment"),
        ylabel=f"{e['measure']} vs reference dose")
    slope_disp = e["slope_display_per_unit"]
    unit = cfg.get("dose_unit", "unit")
    rows = "\n".join(
        f"<tr><td>{_h(s['name'])}</td><td class='num'>{_n(s['slope_display'],3)}</td>"
        f"<td class='num'>{s['n_levels']}</td><td class='num'>{_n(s['weight'],1)}%</td></tr>"
        for s in e["per_study"])
    body = f"""
<section><h2>Result</h2>
<p class="headline">Pooled dose-response trend: each additional {_h(unit)} of dose multiplies
the {e['measure']} by <strong>{_n(slope_disp,3)}</strong>
(95% CI {_n(math.exp(e['slope_ci'][0]) if ratio else e['slope_ci'][0],3)} to
{_n(math.exp(e['slope_ci'][1]) if ratio else e['slope_ci'][1],3)}) across {pool.k} studies.
Heterogeneity I&sup2; {_n(pool.i2,1)}%.</p>
<div class="kpis">
<div class="kpi"><div class="v">{_n(slope_disp,3)}</div><div class="l">{e['measure']} per {_h(unit)}</div></div>
<div class="kpi"><div class="v">{pool.k}</div><div class="l">Studies</div></div>
<div class="kpi"><div class="v">{_n(pool.i2,1)}%</div><div class="l">I&sup2;</div></div>
</div></section>
<section><h2>Dose-response curve</h2><div class="fig">{curve}</div></section>
<section><h2>Per-study slopes</h2>
<table><thead><tr><th>Study</th><th class="num">{e['measure']}/{_h(unit)}</th>
<th class="num">Dose levels</th><th class="num">Weight</th></tr></thead>
<tbody>{rows}</tbody></table></section>
<section><h2>Methods &amp; heterogeneity</h2>{_het_table(pool)}
<p style="font-size:12px;color:#64748b">Approximate two-stage linear trend
(within-study dose-level correlation not reconstructed; needs cases + person-time
per level for the Greenland-Longnecker covariance).</p></section>"""
    return _shell(cfg.get("title", "Dose-response meta-analysis"), "dose-response",
                  _h(cfg.get("outcome", "")), body)
