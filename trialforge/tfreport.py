"""trialforge.tfreport — render the 'Advanced diagnostics' HTML section and
splice it into a metaforge base report."""
from __future__ import annotations
import math


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
    if re_rows:
        blocks.append("<h3>Rare-event estimators</h3>"
                      "<table><thead><tr><th>Method</th><th class='num'>OR (95% CI)</th></tr></thead>"
                      f"<tbody>{''.join(re_rows)}</tbody></table>")

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

    if not blocks:
        return ""
    sn = f"<p style='font-size:12px;color:#64748b'>{source_note}</p>" if source_note else ""
    return ('<section><h2>Advanced diagnostics</h2>' + sn + "".join(blocks) +
            "</section>")


def splice(base_html: str, advanced_html: str) -> str:
    """Insert the advanced section before the disclaimer/footer."""
    if not advanced_html:
        return base_html
    marker = '<section><div class="disclaimer">'
    if marker in base_html:
        return base_html.replace(marker, advanced_html + marker, 1)
    return base_html.replace("</main>", advanced_html + "</main>", 1)
