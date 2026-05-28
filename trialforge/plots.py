"""metaforge.plots — inline SVG (no JS, no libraries, prints + offline).

forest_svg     pairwise / proportion forest with pooled diamond + PI bar
network_svg    NMA network geometry (nodes sized by studies, edges by data)
sucra_svg      NMA SUCRA bar chart
dose_svg       dose-response curve with CI band
"""
from __future__ import annotations
import math


def _esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _fmt(v):
    return "—" if v is None else f"{v:.2f}"


def forest_svg(rows, pooled, *, ratio, null_label_left="favours intervention",
               null_label_right="favours comparator", title="Effect (95% CI)",
               pi=None):
    """rows: [{name, est, lo, hi, weight}]; pooled: {est, lo, hi}."""
    k = len(rows)
    row_h = 26
    top = 64
    width = 760
    x0, x1 = 300, 560
    height = top + (k + 2) * row_h + 96
    null = 1.0 if ratio else 0.0

    pts = []
    for r in rows:
        pts += [r["lo"], r["hi"]]
    pts += [pooled["lo"], pooled["hi"]]
    if pi:
        pts += [pi[0], pi[1]]
    pts = [p for p in pts if p is not None and math.isfinite(p) and (p > 0 if ratio else True)]
    if ratio:
        lo_b, hi_b = max(min(pts) * 0.8, 1e-3), max(pts) * 1.2
        tlo, thi = math.log(lo_b), math.log(hi_b)
        xp = lambda v: x0 + (math.log(max(v, 1e-6)) - tlo) / (thi - tlo) * (x1 - x0)
    else:
        span = (max(pts) - min(pts)) or 1.0
        lo_b, hi_b = min(pts) - span * 0.1, max(pts) + span * 0.1
        xp = lambda v: x0 + (v - lo_b) / (hi_b - lo_b) * (x1 - x0)

    s = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
         f'font-family="system-ui,sans-serif" font-size="12">',
         f'<rect width="{width}" height="{height}" fill="#fff"/>']
    s.append(f'<text x="16" y="30" font-weight="600">Study</text>')
    s.append(f'<text x="200" y="30" font-weight="600" text-anchor="end">Weight</text>')
    s.append(f'<text x="{(x0+x1)//2}" y="30" font-weight="600" text-anchor="middle">{_esc(title)}</text>')
    s.append(f'<line x1="16" y1="40" x2="{width-16}" y2="40" stroke="#cbd5e1"/>')
    xn = xp(null)
    pbot = top + (k + 1) * row_h
    s.append(f'<line x1="{xn:.1f}" y1="{top-6}" x2="{xn:.1f}" y2="{pbot}" stroke="#94a3b8" stroke-dasharray="4 3"/>')
    weights = [r.get("weight") for r in rows if r.get("weight") is not None]
    maxw = max(weights) if weights else 1.0
    for i, r in enumerate(rows):
        y = top + i * row_h + row_h / 2
        s.append(f'<text x="16" y="{y+4:.0f}">{_esc(r["name"])[:32]}</text>')
        if r.get("weight") is not None:
            s.append(f'<text x="200" y="{y+4:.0f}" text-anchor="end" fill="#475569">{r["weight"]:.1f}%</text>')
        xl, xh, xc = xp(r["lo"]), xp(r["hi"]), xp(r["est"])
        xl = max(x0 - 4, min(x1 + 4, xl)); xh = max(x0 - 4, min(x1 + 4, xh))
        s.append(f'<line x1="{xl:.1f}" y1="{y:.1f}" x2="{xh:.1f}" y2="{y:.1f}" stroke="#334155"/>')
        wv = r.get("weight")
        bs = 4 + 8 * math.sqrt((wv if wv is not None else maxw) / maxw)
        s.append(f'<rect x="{xc-bs/2:.1f}" y="{y-bs/2:.1f}" width="{bs:.1f}" height="{bs:.1f}" fill="#2563eb"/>')
        s.append(f'<text x="{x1+12}" y="{y+4:.0f}">{_fmt(r["est"])} ({_fmt(r["lo"])}, {_fmt(r["hi"])})</text>')
    yd = top + k * row_h + row_h / 2 + 4
    xl, xr, xc = xp(pooled["lo"]), xp(pooled["hi"]), xp(pooled["est"])
    xl = max(x0 - 4, xl); xr = min(x1 + 4, xr)
    s.append(f'<polygon points="{xl:.1f},{yd:.1f} {xc:.1f},{yd-7:.1f} {xr:.1f},{yd:.1f} {xc:.1f},{yd+7:.1f}" fill="#dc2626"/>')
    s.append(f'<text x="16" y="{yd+4:.0f}" font-weight="600">Pooled</text>')
    s.append(f'<text x="{x1+12}" y="{yd+4:.0f}" font-weight="600" fill="#b91c1c">{_fmt(pooled["est"])} ({_fmt(pooled["lo"])}, {_fmt(pooled["hi"])})</text>')
    if pi and pi[0] is not None:
        yp = yd + row_h
        xpl, xph = max(x0 - 4, xp(pi[0])), min(x1 + 4, xp(pi[1]))
        s.append(f'<line x1="{xpl:.1f}" y1="{yp:.1f}" x2="{xph:.1f}" y2="{yp:.1f}" stroke="#dc2626" stroke-width="2" stroke-dasharray="2 2"/>')
        s.append(f'<text x="16" y="{yp+4:.0f}" fill="#b91c1c">95% prediction interval</text>')
        s.append(f'<text x="{x1+12}" y="{yp+4:.0f}" fill="#b91c1c">{_fmt(pi[0])} to {_fmt(pi[1])}</text>')
    ay = pbot + 16
    s.append(f'<line x1="{x0}" y1="{ay}" x2="{x1}" y2="{ay}" stroke="#334155"/>')
    if ratio:
        for tick in (0.1, 0.2, 0.5, 1, 2, 5, 10):
            if lo_b <= tick <= hi_b:
                xt = xp(tick)
                s.append(f'<line x1="{xt:.1f}" y1="{ay}" x2="{xt:.1f}" y2="{ay+4}" stroke="#334155"/>')
                s.append(f'<text x="{xt:.1f}" y="{ay+15}" text-anchor="middle" fill="#475569">{tick:g}</text>')
    s.append(f'<text x="{xn-8:.1f}" y="{ay+34}" text-anchor="end" fill="#16a34a">&#9664; {_esc(null_label_left)}</text>')
    s.append(f'<text x="{xn+8:.1f}" y="{ay+34}" fill="#b91c1c">{_esc(null_label_right)} &#9654;</text>')
    s.append("</svg>")
    return "\n".join(s)


def network_svg(treatments, edge_counts, study_counts):
    """Circular network layout. edge_counts: {(a,b): n}; study_counts: {t: n}."""
    T = len(treatments)
    size = 420
    cx, cy, R = size / 2, size / 2, size / 2 - 70
    pos = {}
    for i, t in enumerate(treatments):
        ang = 2 * math.pi * i / T - math.pi / 2
        pos[t] = (cx + R * math.cos(ang), cy + R * math.sin(ang))
    maxe = max(edge_counts.values()) if edge_counts else 1
    maxs = max(study_counts.values()) if study_counts else 1
    s = [f'<svg viewBox="0 0 {size} {size}" xmlns="http://www.w3.org/2000/svg" '
         f'font-family="system-ui,sans-serif" font-size="12">',
         f'<rect width="{size}" height="{size}" fill="#fff"/>']
    for (a, b), n in edge_counts.items():
        x1, y1 = pos[a]; x2, y2 = pos[b]
        w = 1 + 5 * n / maxe
        s.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#94a3b8" stroke-width="{w:.1f}"/>')
    for t in treatments:
        x, y = pos[t]
        r = 8 + 14 * (study_counts.get(t, 1) / maxs)
        s.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" fill="#2563eb" opacity="0.85"/>')
        s.append(f'<text x="{x:.1f}" y="{y-r-4:.1f}" text-anchor="middle" font-weight="600">{_esc(t)}</text>')
    s.append("</svg>")
    return "\n".join(s)


def sucra_svg(ranking):
    """ranking: [{treatment, sucra}] sorted desc."""
    n = len(ranking)
    row_h = 30
    width, height = 560, 50 + n * row_h
    x0, x1 = 160, 500
    s = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
         f'font-family="system-ui,sans-serif" font-size="12">',
         f'<rect width="{width}" height="{height}" fill="#fff"/>',
         f'<text x="16" y="28" font-weight="600">Treatment</text>',
         f'<text x="{(x0+x1)//2}" y="28" font-weight="600" text-anchor="middle">SUCRA (% — higher = ranked better)</text>']
    for i, r in enumerate(ranking):
        y = 44 + i * row_h
        s.append(f'<text x="16" y="{y+14}">{_esc(r["treatment"])[:20]}</text>')
        w = (x1 - x0) * r["sucra"] / 100
        s.append(f'<rect x="{x0}" y="{y}" width="{x1-x0}" height="18" fill="#e2e8f0"/>')
        s.append(f'<rect x="{x0}" y="{y}" width="{w:.1f}" height="18" fill="#2563eb"/>')
        s.append(f'<text x="{x0+w+6:.1f}" y="{y+14}" fill="#1e293b">{r["sucra"]:.0f}%</text>')
    s.append("</svg>")
    return "\n".join(s)


def dose_svg(curve, *, ratio, observed=None, xlabel="Dose increment", ylabel="Effect"):
    """curve: [{dose_increment, effect, ci_low, ci_high}]."""
    width, height = 640, 380
    x0, y0, x1, y1 = 70, 30, 600, 320
    xs = [c["dose_increment"] for c in curve]
    ys = [c["effect"] for c in curve] + [c["ci_low"] for c in curve] + [c["ci_high"] for c in curve]
    xmin, xmax = min(xs), max(xs) or 1
    ymin, ymax = min(ys), max(ys)
    if ratio:
        ymin = min(ymin, 1.0); ymax = max(ymax, 1.0)
    pad = (ymax - ymin) * 0.1 or 0.1
    ymin -= pad; ymax += pad
    xp = lambda v: x0 + (v - xmin) / ((xmax - xmin) or 1) * (x1 - x0)
    yp = lambda v: y1 - (v - ymin) / ((ymax - ymin) or 1) * (y1 - y0)
    s = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
         f'font-family="system-ui,sans-serif" font-size="12">',
         f'<rect width="{width}" height="{height}" fill="#fff"/>']
    # axes
    s.append(f'<line x1="{x0}" y1="{y1}" x2="{x1}" y2="{y1}" stroke="#334155"/>')
    s.append(f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y1}" stroke="#334155"/>')
    # null line
    if ratio:
        yn = yp(1.0)
        s.append(f'<line x1="{x0}" y1="{yn:.1f}" x2="{x1}" y2="{yn:.1f}" stroke="#94a3b8" stroke-dasharray="4 3"/>')
    # CI band
    top_pts = " ".join(f"{xp(c['dose_increment']):.1f},{yp(c['ci_high']):.1f}" for c in curve)
    bot_pts = " ".join(f"{xp(c['dose_increment']):.1f},{yp(c['ci_low']):.1f}" for c in reversed(curve))
    s.append(f'<polygon points="{top_pts} {bot_pts}" fill="#2563eb" opacity="0.15"/>')
    # central line
    line_pts = " ".join(f"{xp(c['dose_increment']):.1f},{yp(c['effect']):.1f}" for c in curve)
    s.append(f'<polyline points="{line_pts}" fill="none" stroke="#2563eb" stroke-width="2"/>')
    # ticks
    for i in range(5):
        xv = xmin + (xmax - xmin) * i / 4
        s.append(f'<text x="{xp(xv):.1f}" y="{y1+16}" text-anchor="middle" fill="#475569">{xv:g}</text>')
        yv = ymin + (ymax - ymin) * i / 4
        s.append(f'<text x="{x0-8}" y="{yp(yv)+4:.1f}" text-anchor="end" fill="#475569">{yv:.2f}</text>')
    s.append(f'<text x="{(x0+x1)//2}" y="{height-6}" text-anchor="middle" font-weight="600">{_esc(xlabel)}</text>')
    s.append(f'<text x="16" y="{(y0+y1)//2}" transform="rotate(-90 16 {(y0+y1)//2})" text-anchor="middle" font-weight="600">{_esc(ylabel)}</text>')
    s.append("</svg>")
    return "\n".join(s)
