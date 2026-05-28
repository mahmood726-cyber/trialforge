# AGENTS.md — instructions for AI CLIs (Gemini CLI, Claude Code, etc.)

trialforge runs **advanced meta-analyses** and can pull trial data straight
from a local **AACT (ClinicalTrials.gov) snapshot**. Your job is cheap:
write a small JSON config, run one command. Do **not** write statistics
code — the engine (`trialforge/`) is finished and tested.

## Workflow

1. Copy a `configs/example_*.json` to `configs/<name>.json`.
2. Fill in the user's data — or set `"source":"aact"` to auto-extract it.
3. Run: `python run.py configs/<name>.json`
4. Open `output/<name>.html`.

`run.py` is deterministic and makes no API calls (the only optional input
is a local AACT snapshot), so you can build many analyses without using quota.
On `[CONFIG ERROR]`, fix only the field named.

## Base analysis (`type`)
`pairwise` | `proportion` | `nma` | `doseresponse` (metaforge engine)
plus three advanced types ported from allmeta:
  * `dta`  — diagnostic test accuracy; studies are 2x2 tables
             `{name, tp, fp, fn, tn}`; output = pooled Se/Sp + SROC curve.
  * `cnma` — additive component NMA; each arm lists `components:[...]`
             (the reference arm has `components:[]`); output = per-component
             effects + combination prediction.
See the per-type field reference in `configs/example_*.json`.

## Two ways to get the data

### A. Inline (you provide the numbers)
```json
{ "type":"pairwise", "measure":"OR", "title":"...",
  "studies":[ {"name":"T1","year":2020,"tE":50,"tN":500,"cE":70,"cN":500}, ... ] }
```

### B. From AACT (auto-extract from a local ClinicalTrials.gov snapshot)
```json
{ "source":"aact", "type":"pairwise", "title":"...",
  "aact": { "ncts":["NCT02540993","NCT02545049"], "force_measure":"HR" },
  "advanced":["egger","loo","cumulative"] }
```
or query instead of listing NCTs:
```json
{ "source":"aact", "type":"pairwise", "title":"...",
  "aact": { "drug":"semaglutide", "condition":"diabetes", "force_measure":"HR", "limit":50 } }
```
- The snapshot is found via the `TRIALFORGE_AACT` env var or `aact.root`.
- `force_measure` (HR/OR/RR/RD/MD) aligns trials on one effect type.
- **AACT extraction is assisted, not authoritative**: outcomes are not
  harmonised across trials. Always tell the user to review the extracted
  values against the source publications before trusting the pool.

## Advanced diagnostics (`advanced: [...]`)

For `pairwise`:
`egger`, `peters` (binary), `trimfill`, `petpeese`, `copas`, `loo`,
`baujat`, `cumulative` (uses each study's `year`), `subgroup` (uses each
study's `subgroup`), `metareg` (uses each study's numeric `moderator`),
`peto` / `mh` (binary rare-events). `copas` profiles the pooled estimate
under increasing assumed publication bias (selection-model sensitivity).

For `nma`:
`loops` (Bucher closed-loop inconsistency — direct vs indirect).

Add only the diagnostics that fit the data: `peters`/`peto`/`mh` need
binary `tE/tN/cE/cN`; `cumulative` needs `year`; `subgroup` needs
`subgroup`; `metareg` needs `moderator`.

## Rules
- Use real numbers / real NCT ids. Never invent trial data.
- pairwise/proportion/doseresponse need >=2 studies; nma needs >=2 treatments.
- The report already carries a "verify against source" disclaimer.

## Don't
- Don't `pip install` (pure stdlib Python 3.8+; numpy used only if present).
- Don't edit `trialforge/` or `run.py`.
