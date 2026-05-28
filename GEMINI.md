# GEMINI.md

Gemini CLI: read **[AGENTS.md](AGENTS.md)** — the full guide.

Summary (save tokens):

> trialforge runs advanced meta-analyses and can auto-extract trial data
> from a local AACT (ClinicalTrials.gov) snapshot. To make one: copy a
> `configs/example_*.json`, set the study data (or `"source":"aact"` with
> NCTs / a drug+condition query), list any `advanced` diagnostics, and run
> `python run.py configs/<file>.json`. The report lands in `output/`. Do
> **not** edit `trialforge/` or `run.py`. AACT extraction is assisted —
> always tell the user to verify extracted numbers against source papers.
> The build is deterministic with no API calls, conserving quota.
