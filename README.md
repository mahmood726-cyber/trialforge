# trialforge

An **AACT-powered, advanced** meta-analysis engine. It does everything
[metaforge](https://github.com/mahmood726-cyber/metaforge) does (pairwise,
proportions, network MA, dose-response) **plus**:

- **Builds datasets directly from a local ClinicalTrials.gov / AACT snapshot**
  — give it NCT ids or a drug + condition and it extracts the effect
  estimates for you.
- A full **advanced-diagnostics suite**: publication-bias tests (Egger,
  Peters, trim-and-fill, PET-PEESE, **Copas** selection-model sensitivity),
  influence (leave-one-out, Baujat), meta-regression, subgroup with
  interaction test, cumulative meta-analysis, rare-event estimators (Peto,
  Mantel-Haenszel), and **NMA inconsistency** via Bucher closed loops.
- Three further advanced methods ported from the `allmeta` suite:
  **diagnostic test accuracy** (`dta` — bivariate Se/Sp pooling + SROC
  curve), **additive component NMA** (`cnma` — decompose multi-component
  interventions), and the **Copas** selection model above. The DTA engine
  is checked against the allmeta `mada::reitsma` fixture; for Copas and
  limit-MA the **unadjusted FE/RE baselines** match the `metasens` fixtures
  to 1e-6, while the bias-**adjusted** estimates are documented HT/WLS
  approximations checked for direction (full MLE is not reproduced).

Same low-token, offline-first design as the sibling kits: the engine is
pre-built, a CLI only writes a small config, and `run.py` is deterministic
with no API calls.

> The four-tool series, increasing power:
> 1. [meta-starter-kit](https://github.com/mahmood726-cyber/meta-starter-kit) — single-file static pairwise
> 2. [rapidmeta-kit](https://github.com/mahmood726-cyber/rapidmeta-kit) — full interactive RapidMeta workbench
> 3. [metaforge](https://github.com/mahmood726-cyber/metaforge) — multi-method engine (pairwise + NMA + dose-response + proportions)
> 4. **trialforge** (this) — AACT data ingestion + advanced diagnostics on top of the metaforge engine

---

## ⬇️ Download (one click)

**[Download as ZIP »](https://github.com/mahmood726-cyber/trialforge/archive/refs/heads/main.zip)**
 · or the [Releases page](https://github.com/mahmood726-cyber/trialforge/releases).

## ▶️ Use it in 3 steps

1. **Install Python** (once): <https://www.python.org/downloads/> (Windows:
   tick "Add Python to PATH").
2. **Run an example.**
   - **Windows:** double-click **`RUN_EXAMPLES.bat`**.
   - **Mac/Linux:** `bash run_examples.sh`.
   - **One at a time:** `python run.py configs/example_pairwise_advanced.json`
3. **Open** the file it prints in `output/` — a full report with forest plot
   *and* an Advanced diagnostics section, working offline.

## Pre-flight data checks

Validate a config and its data before (or instead of) building:

```
python run.py configs/my_review.json --check
```
It reports **errors** (impossible event counts like events > sample size,
inverted or non-positive confidence intervals, no studies) and **warnings**
(too few studies for the prediction interval or publication-bias tests,
duplicate study names, malformed NCT ids). Exit 0 = clean, exit 2 = errors.

A normal build runs the same checks automatically: data **errors block the
build**; warnings are surfaced in the report's "Data checks" section and on
the console. This stops the classic failure mode where a typo'd 2×2 table
silently produces a confident but wrong pooled estimate.

## The AACT superpower

If you have a local [AACT static copy](https://aact.ctti-clinicaltrials.org/)
(the pipe-delimited `.txt` tables), trialforge can build the meta-analysis
for you:

```json
{ "source":"aact", "type":"pairwise",
  "title":"Finerenone cardiorenal HR (from AACT)",
  "aact": { "ncts":["NCT02540993","NCT02545049","NCT04435626"], "force_measure":"HR" },
  "advanced":["egger","loo","cumulative"] }
```
Point trialforge at the snapshot with the `TRIALFORGE_AACT` environment
variable (or `"aact":{"root":"..."}`). It streams only the tables it needs,
filtered by your NCT set, so even multi-GB files are processed in one pass.

You can also **query** instead of listing NCTs:
```json
"aact": { "drug":"semaglutide", "condition":"diabetes", "force_measure":"HR", "limit":50 }
```

> **AACT extraction is assisted, not authoritative.** ClinicalTrials.gov
> outcomes are not harmonised across trials, so the auto-extracted effects
> may mix outcomes. trialforge picks the primary-outcome effect of the
> requested measure where it can, but **you must review every extracted
> value against the source publication** before trusting the pool. The
> report says so prominently.

## Advanced diagnostics

Add a `"advanced": [...]` list to a `pairwise` config:

| Option | What it does |
|---|---|
| `egger` | Egger regression test for funnel asymmetry |
| `peters` | Peters test (better for binary outcomes) |
| `trimfill` | Duval–Tweedie trim-and-fill (sensitivity) |
| `petpeese` | PET-PEESE conditional small-study-effect adjustment |
| `copas` | Copas selection-model publication-bias sensitivity profile |
| `limitma` | Rücker limit meta-analysis (small-study-effect adjusted estimate) |
| `tsa` | trial sequential analysis (RIS + O'Brien-Fleming boundary) |
| `evalue` | E-value: robustness of the pooled estimate to unmeasured confounding |
| `gosh` | GOSH subset-heterogeneity diagnostic (median/IQR across all subsets) |
| `pcurve` | p-curve evidential-value test (needs a `p_value` per study) |
| `glmm` | rare-events binomial/hypergeometric-normal GLMM (no continuity correction) |
| `grade` | GRADE certainty rating + Summary-of-Findings (add `baseline_risk`) |
| `loo` | leave-one-out influence table |
| `baujat` | top heterogeneity contributors |
| `cumulative` | cumulative meta-analysis by `year` |
| `subgroup` | within-`subgroup` pools + interaction test |
| `metareg` | meta-regression on a numeric `moderator` |
| `peto` / `mh` | rare-event Peto / Mantel-Haenszel OR |

For an `nma` config, add `"advanced":["loops","cinema"]` to run Bucher
closed-loop inconsistency (direct vs indirect evidence) and a CINeMA-style
confidence-in-the-evidence rating per comparison.

### Diagnostic test accuracy (`type: "dta"`) and component NMA (`type: "cnma"`)

```json
{ "type":"dta", "title":"...",
  "studies":[ {"name":"Study 1","tp":80,"fp":40,"fn":20,"tn":160}, ... ] }
```
gives pooled sensitivity/specificity, a diagnostic odds ratio, an SROC
curve, and a threshold-effect check.

```json
{ "type":"cnma", "measure":"OR", "title":"...",
  "studies":[ {"name":"S1","arms":[
      {"components":[],"e":120,"n":500},
      {"components":["CBT"],"e":90,"n":500}]}, ... ] }
```
estimates each component's incremental effect under the additive model and
predicts any (even unobserved) combination.

## Worked real-world example (cross-validated)

`configs/example_sglt2_hf.json` reproduces the 5-trial **SGLT2 inhibitors
in heart failure** meta-analysis (DAPA-HF, EMPEROR-Reduced,
EMPEROR-Preserved, DELIVER, SOLOIST-WHF) from the E156 flagship living
capsule. trialforge pools the published hazard ratios to **HR 0.77** with
the same Paule-Mandel + Knapp-Hartung methodology, and
`tests/test_sglt2_flagship.py` asserts trialforge's pooled HR / I² / τ²
match the capsule's own engine to 1e-4 / machine precision — a real
cross-engine validation, not just a synthetic fixture.

## Methods

All pooling uses the metaforge core (random effects, Paule–Mandel τ²,
Knapp–Hartung CI with floor, Cochrane v6.5 prediction interval). The
advanced methods follow standard references (Egger 1997, Peters 2006,
Duval–Tweedie 2000, Stanley–Doucouliagos PET-PEESE, Bucher 1997). Pure
standard-library Python; numpy accelerates NMA if already installed.

## What NOT to edit

Only create/edit files in `configs/`. `trialforge/` (engine + AACT loader +
advanced suite) and `run.py` are finished. `output/` is regenerable.

## Honesty

Every report is an **auto-generated draft**, and AACT extraction is
**assisted** — verify inputs against source publications, confirm clinical
combinability and (for NMA) transitivity, and treat publication-bias tests
as low-powered when there are few studies. Trim-and-fill is a sensitivity
analysis, not a correction.
