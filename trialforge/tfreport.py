"""trialforge.tfreport — render the 'Advanced diagnostics' HTML section and
splice it into a metaforge base report; plus standalone DTA and
component-NMA reports."""
from __future__ import annotations
import math
import json as _json
from . import plots, report as _mf_report


import html as _html


def _h(s):
    return "&mdash;" if s is None else _html.escape(str(s))


def _n(v, nd=3):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "&mdash;"
    return f"{v:.{nd}f}"


def _p(v):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "&mdash;"
    return f"{v:.3f}" if v >= 0.001 else f"{v:.2e}"


def advanced_section(adv: dict, *, ratio: bool, source_note: str = "") -> str:
    """adv: dict of method-name -> result dict (only computed ones present)."""
    blocks = []

    # --- publication bias ---
    pb = []
    if "egger" in adv and adv["egger"].get("available"):
        e = adv["egger"]
        pb.append(f"<tr><td>Egger regression intercept</td>"
                  f"<td class='num'>{_n(e['intercept'])} (t={_n(e['t'],2)}, df={e['df']}, p={_p(e['p'])})</td>"
                  f"<td>{e['interpretation']}</td></tr>")
    if "peters" in adv and adv["peters"].get("available"):
        e = adv["peters"]
        pb.append(f"<tr><td>Peters test (binary)</td>"
                  f"<td class='num'>slope p={_p(e['p'])}</td><td>{e['interpretation']}</td></tr>")
    if "trimfill" in adv and adv["trimfill"].get("available"):
        e = adv["trimfill"]
        disp = (lambda v: math.exp(v)) if ratio else (lambda v: v)
        pb.append(f"<tr><td>Trim-and-fill (sensitivity)</td>"
                  f"<td class='num'>{e['k_imputed']} study(ies) imputed</td>"
                  f"<td>adjusted estimate {_n(disp(e['adjusted_estimate']))} "
                  f"(was {_n(disp(e['original_estimate']))})</td></tr>")
    if "petpeese" in adv and adv["petpeese"].get("available"):
        e = adv["petpeese"]
        disp = (lambda v: math.exp(v)) if ratio else (lambda v: v)
        pb.append(f"<tr><td>PET-PEESE (chose {e['chosen']})</td>"
                  f"<td class='num'>{_n(disp(e['adjusted_estimate']))} "
                  f"({_n(disp(e['ci_low']))}, {_n(disp(e['ci_high']))})</td>"
                  f"<td>small-study-adjusted effect</td></tr>")
    if pb:
        blocks.append("<h3>Publication bias / small-study effects</h3>"
                      "<table><thead><tr><th>Test</th><th class='num'>Result</th>"
                      "<th>Interpretation</th></tr></thead><tbody>"
                      + "\n".join(pb) + "</tbody></table>")

    # --- Copas selection-model sensitivity ---
    if "copas" in adv and adv["copas"].get("available"):
        c = adv["copas"]
        prof = "\n".join(
            f"<tr><td class='num'>{int(p['p_unpublished']*100)}%</td>"
            f"<td class='num'>{_n(p['estimate'])} ({_n(p['ci_low'])}, {_n(p['ci_high'])})</td></tr>"
            for p in c["profile"])
        blocks.append("<h3>Copas selection-model sensitivity</h3>"
                      "<table><thead><tr><th class='num'>Assumed unpublished</th>"
                      "<th class='num'>Adjusted estimate (95% CI)</th></tr></thead>"
                      f"<tbody>{prof}</tbody></table>"
                      f"<p style='font-size:12px;color:#64748b'>Unadjusted RE "
                      f"{_n(c['re_estimate'])}; worst-case (50% unpublished) "
                      f"{_n(c['worst_case']['estimate'])}. {c['note']}</p>")

    # --- limit meta-analysis ---
    if "limitma" in adv and adv["limitma"].get("available"):
        L = adv["limitma"]
        blocks.append("<h3>Limit meta-analysis (small-study-effect adjusted)</h3>"
                      "<table>"
                      f"<tr><td>Random-effects estimate</td><td class='num'>{_n(L['re_estimate'])} "
                      f"({_n(L['re_ci'][0])}, {_n(L['re_ci'][1])})</td></tr>"
                      f"<tr><td>Limit (adjusted) estimate</td><td class='num'>{_n(L['limit_estimate'])} "
                      f"({_n(L['limit_ci'][0])}, {_n(L['limit_ci'][1])})</td></tr>"
                      f"<tr><td>Rucker regression slope</td><td class='num'>{_n(L['slope'],3)} "
                      f"({'small-study effect present' if L['small_study_effect'] else 'no small-study effect'})</td></tr>"
                      "</table>"
                      f"<p style='font-size:12px;color:#64748b'>{L['note']}</p>")

    # --- E-value ---
    if "evalue" in adv and adv["evalue"].get("available"):
        E = adv["evalue"]
        blocks.append("<h3>E-value (robustness to unmeasured confounding)</h3>"
                      "<table>"
                      f"<tr><td>E-value (point estimate)</td><td class='num'>{_n(E['evalue_point'],2)}</td></tr>"
                      f"<tr><td>E-value (CI limit nearest null)</td><td class='num'>{_n(E['evalue_ci'],2)}</td></tr>"
                      "</table>"
                      f"<p style='font-size:12px;color:#64748b'>{E['interpretation']} {E['note']}</p>")

    # --- trial sequential analysis ---
    if "tsa" in adv and adv["tsa"].get("available"):
        T = adv["tsa"]
        blocks.append("<h3>Trial sequential analysis</h3>"
                      "<table>"
                      f"<tr><td>Cumulative z</td><td class='num'>{_n(T['z_cumulative'],2)}</td></tr>"
                      f"<tr><td>O'Brien-Fleming boundary z</td><td class='num'>{_n(T['z_boundary'],2)}</td></tr>"
                      f"<tr><td>Information fraction (accrued / required)</td>"
                      f"<td class='num'>{_n(T['information_fraction']*100,1)}%</td></tr>"
                      f"<tr><td>Required information size (heterogeneity-adjusted)</td>"
                      f"<td class='num'>D&times; diversity = {_n(T['diversity_D'],2)}</td></tr>"
                      "</table>"
                      f"<p style='font-size:12px;color:#64748b'><strong>{T['conclusion']}.</strong> "
                      f"{T['note']}</p>")

    # --- influence ---
    if "loo" in adv and adv["loo"]:
        disp = (lambda v: math.exp(v)) if ratio else (lambda v: v)
        rows = "\n".join(
            f"<tr><td>{r['omitted']}</td><td class='num'>{_n(disp(r['estimate']))} "
            f"({_n(disp(r['ci_low']))}, {_n(disp(r['ci_high']))})</td>"
            f"<td class='num'>{_n(r['i2'],1)}%</td></tr>" for r in adv["loo"])
        blocks.append("<h3>Leave-one-out influence</h3>"
                      "<table><thead><tr><th>Study omitted</th>"
                      "<th class='num'>Pooled (95% CI)</th><th class='num'>I&sup2;</th></tr></thead>"
                      f"<tbody>{rows}</tbody></table>")
    if "baujat" in adv and adv["baujat"]:
        top = sorted(adv["baujat"], key=lambda r: -r["x_heterogeneity"])[:5]
        rows = "\n".join(
            f"<tr><td>{r['name']}</td><td class='num'>{_n(r['x_heterogeneity'],2)}</td>"
            f"<td class='num'>{_n(r['y_influence'],3)}</td></tr>" for r in top)
        blocks.append("<h3>Baujat (top heterogeneity contributors)</h3>"
                      "<table><thead><tr><th>Study</th>"
                      "<th class='num'>Q contribution</th><th class='num'>Influence</th></tr></thead>"
                      f"<tbody>{rows}</tbody></table>")

    # --- cumulative ---
    if "cumulative" in adv and adv["cumulative"]:
        disp = (lambda v: math.exp(v)) if ratio else (lambda v: v)
        rows = "\n".join(
            f"<tr><td>{r['order']}</td><td>+ {r['added']}</td><td class='num'>{r['k']}</td>"
            f"<td class='num'>{_n(disp(r['estimate']))} ({_n(disp(r['ci_low']))}, {_n(disp(r['ci_high']))})</td></tr>"
            for r in adv["cumulative"])
        blocks.append("<h3>Cumulative meta-analysis</h3>"
                      "<table><thead><tr><th>Order</th><th>Study added</th>"
                      "<th class='num'>k</th><th class='num'>Pooled (95% CI)</th></tr></thead>"
                      f"<tbody>{rows}</tbody></table>")

    # --- subgroup ---
    if "subgroup" in adv and adv["subgroup"]:
        sg = adv["subgroup"]
        disp = (lambda v: math.exp(v)) if ratio else (lambda v: v)
        rows = "\n".join(
            f"<tr><td>{g}</td><td class='num'>{r['k']}</td>"
            f"<td class='num'>{_n(disp(r['estimate']))} ({_n(disp(r['ci_low']))}, {_n(disp(r['ci_high']))})</td>"
            f"<td class='num'>{_n(r['i2'],1)}%</td></tr>" for g, r in sg["subgroups"].items())
        blocks.append("<h3>Subgroup analysis</h3>"
                      "<table><thead><tr><th>Subgroup</th><th class='num'>k</th>"
                      "<th class='num'>Pooled (95% CI)</th><th class='num'>I&sup2;</th></tr></thead>"
                      f"<tbody>{rows}</tbody></table>"
                      f"<p style='font-size:12px;color:#64748b'>Between-subgroup test: "
                      f"Q={_n(sg['Q_between'],2)} (df {sg['df_between']}), p={_p(sg['p_between'])} "
                      f"&mdash; {'significant interaction' if (sg['p_between']==sg['p_between'] and sg['p_between']<0.05) else 'no significant interaction'}.</p>")

    # --- meta-regression ---
    if "metareg" in adv and adv["metareg"].get("available"):
        m = adv["metareg"]
        blocks.append("<h3>Meta-regression (single moderator)</h3>"
                      "<table>"
                      f"<tr><td>Slope (per unit moderator, log scale if ratio)</td>"
                      f"<td class='num'>{_n(m['slope'])} (95% CI {_n(m['ci_low'])} to {_n(m['ci_high'])})</td></tr>"
                      f"<tr><td>Test of slope</td><td class='num'>t={_n(m['t'],2)} (df {m['df']}), p={_p(m['p'])}</td></tr>"
                      f"<tr><td>Residual &tau;&sup2;</td><td class='num'>{_n(m['tau2_residual'],4)}</td></tr>"
                      "</table>")

    # --- rare events ---
    re_rows = []
    if "peto" in adv and adv["peto"].get("available"):
        e = adv["peto"]
        re_rows.append(f"<tr><td>Peto one-step OR</td>"
                       f"<td class='num'>{_n(e['OR'])} ({_n(e['ci_low'])}, {_n(e['ci_high'])})</td></tr>")
    if "mh" in adv and adv["mh"].get("available"):
        e = adv["mh"]
        re_rows.append(f"<tr><td>Mantel-Haenszel OR (fixed)</td>"
                       f"<td class='num'>{_n(e['OR'])} ({_n(e['ci_low'])}, {_n(e['ci_high'])})</td></tr>")
    if "glmm" in adv and adv["glmm"].get("available"):
        e = adv["glmm"]
        re_rows.append(f"<tr><td>Binomial-normal GLMM (no continuity correction)</td>"
                       f"<td class='num'>{_n(e['OR'])} ({_n(e['ci_low'])}, {_n(e['ci_high'])}); "
                       f"&tau;&sup2;={_n(e['tau2'],3)}</td></tr>")
    if re_rows:
        blocks.append("<h3>Rare-event estimators</h3>"
                      "<table><thead><tr><th>Method</th><th class='num'>OR (95% CI)</th></tr></thead>"
                      f"<tbody>{''.join(re_rows)}</tbody></table>")

    # --- GRADE certainty + SoF ---
    if "grade" in adv and adv["grade"].get("available"):
        G = adv["grade"]
        drows = "\n".join(f"<tr><td>{_h(k.replace('_',' '))}</td><td>{_h(v)}</td></tr>"
                          for k, v in G["domains"].items())
        sof = ""
        if G["absolute_effect"]:
            ae = G["absolute_effect"]
            sof = ("<p style='font-size:13px'><strong>Absolute effect</strong> "
                   f"(per 1000): comparator {ae['comparator_per_1000']}, "
                   f"intervention {ae['intervention_per_1000']} "
                   f"(difference {ae['difference_per_1000']}, 95% CI "
                   f"{ae['diff_ci_per_1000'][0]} to {ae['diff_ci_per_1000'][1]}).</p>")
        blocks.append("<h3>GRADE certainty of evidence</h3>"
                      f"<p class='headline'>Certainty: <strong>{_h(G['certainty']).upper()}</strong> "
                      f"({G['n_studies']} studies, {_h(G['n_participants'])} participants; "
                      f"starting design: {_h(G['design'])}).</p>"
                      "<table><thead><tr><th>Domain</th><th>Rating</th></tr></thead>"
                      f"<tbody>{drows}</tbody></table>{sof}"
                      f"<p style='font-size:12px;color:#64748b'>{G['note']}</p>")

    # --- NMA inconsistency ---
    if "loops" in adv and adv["loops"]:
        li = adv["loops"]
        if li["n_loops"]:
            rows = "\n".join(
                f"<tr><td>{l['loop']}</td><td>{l['edge']}</td>"
                f"<td class='num'>{_n(l['direct'])}</td><td class='num'>{_n(l['indirect'])}</td>"
                f"<td class='num'>{_n(l['IF'])}</td><td class='num'>{_p(l['p'])}</td>"
                f"<td>{'&#9888; inconsistent' if l['inconsistent'] else 'consistent'}</td></tr>"
                for l in li["loops"])
            blocks.append("<h3>Network inconsistency (Bucher closed loops)</h3>"
                          "<table><thead><tr><th>Loop</th><th>Edge</th>"
                          "<th class='num'>Direct</th><th class='num'>Indirect</th>"
                          "<th class='num'>IF</th><th class='num'>p</th><th>Verdict</th></tr></thead>"
                          f"<tbody>{rows}</tbody></table>"
                          f"<p style='font-size:12px;color:#64748b'>{li['n_loops']} closed loop(s); "
                          f"{'at least one shows inconsistency' if li['any_inconsistent'] else 'no significant inconsistency detected'}.</p>")

    # --- GOSH subset diagnostic ---
    if "gosh" in adv and adv["gosh"].get("available"):
        g = adv["gosh"]
        blocks.append("<h3>GOSH (subset-heterogeneity diagnostic)</h3>"
                      "<table>"
                      f"<tr><td>Subsets pooled</td><td class='num'>{g['n_subsets']}"
                      f"{' (sampled)' if g['sampled'] else ''}</td></tr>"
                      f"<tr><td>Full-sample estimate</td><td class='num'>{_n(g['full_estimate'])} "
                      f"(I&sup2; {_n(g['full_i2'],1)}%)</td></tr>"
                      f"<tr><td>Subset estimates (median, IQR)</td><td class='num'>{_n(g['median_estimate'])} "
                      f"[{_n(g['q25_estimate'])}, {_n(g['q75_estimate'])}]</td></tr>"
                      f"<tr><td>Subset estimate range</td><td class='num'>{_n(g['min_estimate'])} to {_n(g['max_estimate'])}</td></tr>"
                      f"<tr><td>Median / max subset I&sup2;</td><td class='num'>{_n(g['median_i2'],1)}% / {_n(g['max_i2'],1)}%</td></tr>"
                      "</table>"
                      f"<p style='font-size:12px;color:#64748b'>{g['note']}</p>")

    # --- p-curve ---
    if "pcurve" in adv and adv["pcurve"].get("available"):
        pc = adv["pcurve"]
        blocks.append("<h3>p-curve (evidential value)</h3>"
                      "<table>"
                      f"<tr><td>Significant studies used</td><td class='num'>{pc['n_significant']}</td></tr>"
                      f"<tr><td>Fisher's test (right-skew)</td><td class='num'>&chi;&sup2;={_n(pc['fisher_chisq'],2)} "
                      f"(df {pc['fisher_df']}), p={_p(pc['fisher_p'])}</td></tr>"
                      f"<tr><td>Flatness test</td><td class='num'>{int(pc['prop_low']*100)}% of pp&le;0.5, "
                      f"z={_n(pc['z_flat'],2)}, p={_p(pc['p_flat'])}</td></tr>"
                      "</table>"
                      f"<p style='font-size:12px;color:#64748b'>{pc['interpretation']}. {pc['note']}</p>")

    # --- CINeMA confidence ---
    if "cinema" in adv and adv["cinema"].get("available"):
        cm = adv["cinema"]
        rows = "\n".join(
            f"<tr><td>{_h(c['comparison'])}</td>"
            f"<td class='num'>{_n(c['estimate'])} ({_n(c['ci_low'])}, {_n(c['ci_high'])})</td>"
            f"<td>{c['domains']['imprecision']}</td>"
            f"<td>{c['domains']['heterogeneity']}</td>"
            f"<td>{c['domains']['incoherence']}</td>"
            f"<td><strong>{c['confidence']}</strong></td></tr>" for c in cm["comparisons"])
        blocks.append("<h3>Confidence in NMA (CINeMA-style)</h3>"
                      "<table><thead><tr><th>Comparison</th><th class='num'>Effect (95% CI)</th>"
                      "<th>Imprecision</th><th>Heterogeneity</th><th>Incoherence</th>"
                      "<th>Confidence</th></tr></thead>"
                      f"<tbody>{rows}</tbody></table>"
                      f"<p style='font-size:12px;color:#64748b'>{cm['note']}</p>")

    if not blocks:
        return ""
    sn = f"<p style='font-size:12px;color:#64748b'>{source_note}</p>" if source_note else ""
    return ('<section><h2>Advanced diagnostics</h2>' + sn + "".join(blocks) +
            "</section>")


# --------------------------------------------------------------------------
# WebR (in-browser R + metafor) cross-validation panel.
# Like the E156 living capsule, the report can verify its own pooled estimate
# by running metafor::rma() against the SAME yi/sei in the reader's browser.
# Plain (non-f-string) JS + a JSON data tag avoids brace / </script> escaping
# pitfalls (see rules/lessons.md).
# --------------------------------------------------------------------------
_WEBR_JS = r"""
(function(){
  var el = document.getElementById('tf-webr-data');
  if(!el) return;
  var D = JSON.parse(el.textContent);
  var btn = document.getElementById('tf-webr-btn');
  var out = document.getElementById('tf-webr-out');
  function fmt(x,n){ if(x===null||x===undefined||isNaN(x)) return '—';
    return Number(x).toFixed(n===undefined?4:n); }
  btn.addEventListener('click', async function(){
    btn.disabled = true;
    out.textContent = 'Loading R + metafor in your browser (first run ~20–40s)…';
    try{
      var mod = await import('https://webr.r-wasm.org/latest/webr.mjs');
      var webR = new mod.WebR();
      await webR.init();
      out.textContent = 'Installing metafor…';
      await webR.installPackages(['metafor']);
      out.textContent = 'Running rma()…';
      var knha = D.knha ? ',test="knha"' : '';
      var transf = D.ratio ? 'transf=exp,' : '';
      var code = 'suppressMessages(library(metafor)); yi<-c('+D.yi.join(',')+
        '); sei<-c('+D.sei.join(',')+'); res<-rma(yi=yi,sei=sei,method="'+D.method+
        '"'+knha+'); pr<-predict(res,'+transf+'pi=TRUE); '+
        'c(pr$pred,pr$ci.lb,pr$ci.ub,pr$pi.lb,pr$pi.ub,res$I2)';
      var r = await webR.evalR(code);
      var js = await r.toJs();
      try{ webR.destroy(r); }catch(e){}
      var v = js.values || [];
      var pred=v[0], clb=v[1], cub=v[2], plb=v[3], pub=v[4], i2=v[5];
      var tf = D.tf;
      var d = Math.abs(pred - tf.estimate);
      var ok = d < 1e-3;
      var lab = D.measure + ', ' + D.method + (D.knha ? ', Knapp–Hartung' : '');
      out.innerHTML =
        '<div><b>R metafor ('+lab+'):</b> '+fmt(pred)+' ('+fmt(clb)+'–'+fmt(cub)+'), '+
        'PI '+fmt(plb)+'–'+fmt(pub)+', I²='+fmt(i2,0)+'%</div>'+
        '<div><b>trialforge:</b> '+fmt(tf.estimate)+' ('+fmt(tf.ci_low)+'–'+fmt(tf.ci_high)+'), '+
        'PI '+fmt(tf.pi_low)+'–'+fmt(tf.pi_high)+', I²='+fmt(tf.i2,0)+'%</div>'+
        '<div style="margin-top:6px">Δ pooled point estimate = <b>'+d.toExponential(1)+'</b> '+
        (ok ? '<span style="color:#16a34a;font-weight:700">✓ agreement</span>'
            : '<span style="color:#b45309;font-weight:700">differs — investigate</span>')+'</div>'+
        (D.floor ? '<div style="font-size:12px;color:#64748b;margin-top:4px">The 95% CIs '+
          'differ by design here: trialforge floors the Knapp–Hartung variance scale at 1 '+
          'under under-dispersion (Q&lt;k−1), which metafor does not. The pooled point '+
          'estimate is the scale-invariant cross-check.</div>' : '');
    }catch(e){
      out.innerHTML = '<span style="color:#b45309">WebR could not run here — it needs '+
        'network access and an https:// origin (a local file:// page is often blocked by the '+
        'browser, and the metafor WebAssembly build must be reachable): '+
        String((e&&e.message)||e)+'.</span> The offline trialforge result and the copyable '+
        'R script above remain valid.';
    }finally{ btn.disabled = false; }
  });
})();
"""


def webr_section(yi, sei, *, ratio, method, knha, measure,
                 tf_estimate, tf_ci_low, tf_ci_high,
                 tf_pi_low, tf_pi_high, tf_i2, floor_applied=False) -> str:
    """A 'Verify in R (WebR)' section: a copyable metafor script + a button
    that runs rma() on the SAME yi/sei in-browser and diffs R's pooled
    estimate against trialforge's. Degrades gracefully offline."""
    if not yi or len(yi) < 2:
        return ""
    data = {
        "yi": [round(float(v), 8) for v in yi],
        "sei": [round(float(s), 8) for s in sei],
        "ratio": bool(ratio), "method": method, "knha": bool(knha),
        "measure": measure, "floor": bool(floor_applied),
        "tf": {"estimate": tf_estimate, "ci_low": tf_ci_low, "ci_high": tf_ci_high,
               "pi_low": tf_pi_low, "pi_high": tf_pi_high, "i2": tf_i2},
    }
    data_json = _json.dumps(data).replace("<", "\\u003c")
    yi_str = ", ".join(f"{float(v):.6f}" for v in yi)
    sei_str = ", ".join(f"{float(s):.6f}" for s in sei)
    knha_arg = ', test="knha"' if knha else ""
    transf = "transf=exp, " if ratio else ""
    scale_note = " (back-transformed)" if ratio else ""
    rcode = (
        "library(metafor)\n"
        f"yi  <- c({yi_str})\n"
        f"sei <- c({sei_str})\n"
        f'res <- rma(yi = yi, sei = sei, method = "{method}"{knha_arg})\n'
        f"predict(res, {transf}pi = TRUE)  # pooled {measure}{scale_note}, 95% CI, 95% PI"
    )
    return (
        '<section><h2>Verify in R (WebR)</h2>'
        '<p style="font-size:13px;color:#334155">This report can check its own '
        'arithmetic. The script below pools the <em>same</em> per-study effects '
        '(<code>yi</code>) and standard errors (<code>sei</code>) with '
        '<code>metafor::rma()</code>. Run it in any R session, or press the button '
        'to run it right here in your browser via '
        '<a href="https://docs.r-wasm.org/webr/latest/" rel="noopener">WebR</a> '
        '(R + metafor compiled to WebAssembly) and compare against trialforge’s '
        'pooled estimate.</p>'
        f'<pre style="background:#0f172a;color:#e2e8f0;padding:12px 14px;border-radius:8px;'
        f'overflow-x:auto;font-size:12.5px;line-height:1.5">{_h(rcode)}</pre>'
        '<p style="margin:10px 0"><button id="tf-webr-btn" type="button" '
        'style="background:#2563eb;color:#fff;border:0;border-radius:7px;'
        'padding:9px 16px;font-size:14px;font-weight:600;cursor:pointer">'
        '▶ Run metafor in this browser</button></p>'
        '<div id="tf-webr-out" style="font-size:13.5px;line-height:1.6;'
        'background:#f8fafc;border:1px solid #e2e8f0;border-radius:8px;padding:12px 14px;'
        'min-height:1.6em;color:#334155">Press the button to load R + metafor and '
        'cross-validate the pooled estimate (needs network + an https:// origin).</div>'
        '<script type="application/json" id="tf-webr-data">' + data_json + '</script>'
        '<script>' + _WEBR_JS + '</script>'
        '</section>'
    )


def splice(base_html: str, advanced_html: str) -> str:
    """Insert the advanced section before the disclaimer/footer."""
    if not advanced_html:
        return base_html
    marker = '<section><div class="disclaimer">'
    if marker in base_html:
        return base_html.replace(marker, advanced_html + marker, 1)
    return base_html.replace("</main>", advanced_html + "</main>", 1)


def checks_section(issues):
    """Render a 'Data checks' section from checks.check() issues."""
    if not issues:
        return ('<section><h2>Data checks</h2>'
                '<p style="color:#16a34a;font-size:13px">No data-quality issues '
                'detected by the automated pre-flight checks.</p></section>')
    colour = {"error": "#b91c1c", "warning": "#92400e", "info": "#475569"}
    rows = "\n".join(
        f"<tr><td style='color:{colour.get(x['level'],'#475569')};font-weight:600'>"
        f"{x['level'].upper()}</td><td>{_h(x['code'])}</td><td>{_h(x['msg'])}</td></tr>"
        for x in issues)
    return ('<section><h2>Data checks</h2>'
            '<table><thead><tr><th>Level</th><th>Code</th><th>Detail</th></tr></thead>'
            f'<tbody>{rows}</tbody></table>'
            '<p style="font-size:12px;color:#64748b">Automated pre-flight checks. '
            'Errors block the analysis; warnings are advisory.</p></section>')


def _pct(v):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "&mdash;"
    return f"{v*100:.1f}%"


def render_dta(cfg, res):
    """Standalone diagnostic test accuracy report."""
    sroc = plots.sroc_svg(res["per_study"],
                          {"sensitivity": res["sensitivity"], "specificity": res["specificity"]},
                          res["sroc"])
    rows = "\n".join(
        f"<tr><td>{_h(p['name'])}</td><td class='num'>{_pct(p['se'])}</td>"
        f"<td class='num'>{_pct(p['sp'])}</td></tr>" for p in res["per_study"])
    te = ("Strong threshold effect detected (Se and FPR strongly correlated) "
          "&mdash; prefer the SROC curve over a single pooled point."
          if res["threshold_effect"] else
          "No strong threshold effect; the pooled operating point is reasonable.")
    body = f"""
<section><h2>Result</h2>
<p class="headline">Pooled sensitivity <strong>{_pct(res['sensitivity'])}</strong>
(95% CI {_pct(res['sensitivity_ci'][0])} to {_pct(res['sensitivity_ci'][1])}) and
pooled specificity <strong>{_pct(res['specificity'])}</strong>
(95% CI {_pct(res['specificity_ci'][0])} to {_pct(res['specificity_ci'][1])})
across {res['k']} studies. Diagnostic odds ratio {_n(res['dor'],1)}.</p>
<div class="kpis">
<div class="kpi"><div class="v">{_pct(res['sensitivity'])}</div><div class="l">Sensitivity</div></div>
<div class="kpi"><div class="v">{_pct(res['specificity'])}</div><div class="l">Specificity</div></div>
<div class="kpi"><div class="v">{_n(res['dor'],1)}</div><div class="l">Diagnostic OR</div></div>
<div class="kpi"><div class="v">{res['k']}</div><div class="l">Studies</div></div>
</div></section>
<section><h2>Summary ROC (SROC)</h2><div class="fig">{sroc}</div>
<p style="font-size:12px;color:#64748b">Red point = pooled operating point; blue line = SROC curve;
grey points = individual studies.</p></section>
<section><h2>Per-study estimates</h2>
<table><thead><tr><th>Study</th><th class="num">Sensitivity</th><th class="num">Specificity</th></tr></thead>
<tbody>{rows}</tbody></table></section>
<section><h2>Heterogeneity &amp; threshold</h2>
<table>
<tr><td>&tau;&sup2; logit(Se) / logit(FPR)</td><td class="num">{_n(res['tau2_se'],3)} / {_n(res['tau2_fpr'],3)}</td></tr>
<tr><td>Between-study correlation &rho;</td><td class="num">{_n(res['rho'],3)}</td></tr>
<tr><td>Threshold correlation (Spearman)</td><td class="num">{_n(res['threshold_corr'],3)}</td></tr>
</table>
<p style="font-size:12px;color:#64748b">{te} {res['note']}</p></section>"""
    return _mf_report._shell(cfg.get("title", "Diagnostic test accuracy"),
                             "diagnostic accuracy",
                             _h(cfg.get("outcome", "Sensitivity & specificity")), body)


def render_cnma(cfg, res):
    """Standalone additive component-NMA report."""
    ratio = res["ratio"]
    rows = [{"name": c["component"], "est": c["effect"], "lo": c["ci_low"],
             "hi": c["ci_high"], "weight": None} for c in res["component_effects"]]
    forest = plots.forest_svg(
        rows, {"est": 1.0 if ratio else 0.0, "lo": 1.0 if ratio else 0.0,
               "hi": 1.0 if ratio else 0.0},
        ratio=ratio, title=f"Component effect ({res['measure']}, 95% CI)",
        null_label_left="component lowers the outcome",
        null_label_right="component raises it")
    crows = "\n".join(
        f"<tr><td>{_h(c['component'])}</td>"
        f"<td class='num'>{_n(c['effect'])} ({_n(c['ci_low'])}, {_n(c['ci_high'])})</td></tr>"
        for c in res["component_effects"])
    body = f"""
<section><h2>Component effects</h2>
<p class="headline">Additive component NMA across {res['n_contrasts']} contrasts
decomposes the interventions into {len(res['components'])} components. Each
component's incremental {res['measure']} (vs not having it) is estimated below;
a combination's effect is the sum of its components (on the {'log ' if ratio else ''}scale).</p>
<div class="fig">{forest}</div></section>
<section><h2>Estimated component effects</h2>
<table><thead><tr><th>Component</th><th class="num">{res['measure']} (95% CI)</th></tr></thead>
<tbody>{crows}</tbody></table>
<p style="font-size:12px;color:#64748b">{res['note']}</p></section>"""
    return _mf_report._shell(cfg.get("title", "Component network meta-analysis"),
                             "component NMA",
                             f"{len(res['components'])} components &middot; {res['measure']}", body)


def render_multivariate(cfg, res):
    """Standalone bivariate multi-outcome report."""
    o1, o2 = res["outcome1"], res["outcome2"]
    def block(name, o):
        if o is None:
            return ""
        return (f"<tr><td>{_h(name)}</td>"
                f"<td class='num'>{_n(o['bivariate'])} ({_n(o['ci_low'])}, {_n(o['ci_high'])})</td>"
                f"<td class='num'>{_n(o['se_univariate'],4)}</td>"
                f"<td class='num'>{_n(o['se_bivariate'],4)}</td>"
                f"<td class='num'>{_n(o['borrowed_precision_pct'],1)}%</td></tr>")
    n1 = cfg.get("outcome1_name", "Outcome 1")
    n2 = cfg.get("outcome2_name", "Outcome 2")
    ratio = cfg.get("ratio", False)
    frows = [{"name": n1, "est": o1["bivariate"], "lo": o1["ci_low"],
              "hi": o1["ci_high"], "weight": None}]
    if o2 is not None:
        frows.append({"name": n2, "est": o2["bivariate"], "lo": o2["ci_low"],
                      "hi": o2["ci_high"], "weight": None})
    forest = plots.forest_svg(
        frows, {"est": 1.0 if ratio else 0.0, "lo": 1.0 if ratio else 0.0,
                "hi": 1.0 if ratio else 0.0},
        ratio=ratio, title="Pooled effect (95% CI)",
        null_label_left="lower", null_label_right="higher")
    body = f"""
<section><h2>Bivariate result (borrowing of strength)</h2>
<p class="headline">Jointly modelling two correlated outcomes
(between-study correlation &rho; = {_n(res['rho_between'],2)}) lets each
outcome borrow precision from the other. The "borrowed precision" column
shows how much each outcome's standard error shrank versus a univariate
analysis.</p>
<div class="fig">{forest}</div>
<table><thead><tr><th>Outcome</th><th class="num">Bivariate estimate (95% CI)</th>
<th class="num">SE univariate</th><th class="num">SE bivariate</th>
<th class="num">Precision gained</th></tr></thead>
<tbody>{block(n1, o1)}{block(n2, o2)}</tbody></table>
<p style="font-size:12px;color:#64748b">{res['note']} &tau;&#8321;={_n(res['tau1'],3)},
&tau;&#8322;={_n(res['tau2'],3)}. {res['k_both']} of {res['k']} studies reported both outcomes.</p>
</section>"""
    return _mf_report._shell(cfg.get("title", "Multivariate (multi-outcome) meta-analysis"),
                             "multivariate", f"2 outcomes &middot; &rho;={_n(res['rho_between'],2)}", body)


def render_rmst(cfg, res):
    """Standalone RMST / survival report."""
    rows = "\n".join(
        f"<tr><td>{_h(r['name'])}</td><td class='num'>{_n(r['diff'],2)} "
        f"({_n(r['lo'],2)}, {_n(r['hi'],2)})</td><td class='num'>{_n(r['weight'],1)}%</td></tr>"
        for r in res["per_study"])
    rowsf = [{"name": r["name"], "est": r["diff"], "lo": r["lo"], "hi": r["hi"],
              "weight": r["weight"]} for r in res["per_study"]]
    forest = plots.forest_svg(
        rowsf, {"est": res["rmst_difference"], "lo": res["ci_low"], "hi": res["ci_high"]},
        ratio=False, title="RMST difference (95% CI)",
        pi=(res["pi_low"], res["pi_high"]) if res["pi_low"] is not None else None,
        null_label_left="favours comparator", null_label_right="favours intervention")
    tau_txt = f"{res['tau_star']}" if res["tau_star"] is not None else "study-specified"
    body = f"""
<section><h2>Result</h2>
<p class="headline">Pooled RMST difference (intervention &minus; comparator) over
&tau;* = {tau_txt}: <strong>{_n(res['rmst_difference'],2)}</strong> time units
(95% CI {_n(res['ci_low'],2)} to {_n(res['ci_high'],2)}) across {res['k']} studies.
Heterogeneity I&sup2; {_n(res['i2'],1)}%.</p>
<div class="kpis">
<div class="kpi"><div class="v">{_n(res['rmst_difference'],2)}</div><div class="l">RMST difference</div></div>
<div class="kpi"><div class="v">{res['k']}</div><div class="l">Studies</div></div>
<div class="kpi"><div class="v">{_n(res['i2'],1)}%</div><div class="l">I&sup2;</div></div>
</div></section>
<section><h2>Forest plot</h2><div class="fig">{forest}</div></section>
<section><h2>Per-study RMST differences</h2>
<table><thead><tr><th>Study</th><th class="num">RMST diff (95% CI)</th><th class="num">Weight</th></tr></thead>
<tbody>{rows}</tbody></table>
<p style="font-size:12px;color:#64748b">{res['note']}</p></section>"""
    return _mf_report._shell(cfg.get("title", "RMST / survival meta-analysis"),
                             "RMST survival", f"&tau;* = {tau_txt}", body)
