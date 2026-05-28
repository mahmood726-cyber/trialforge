"""trialforge.aact — build meta-analysis datasets directly from an AACT
(ClinicalTrials.gov) static snapshot.

AACT ships pipe-delimited flat files. This module locates a snapshot, then
scans only the tables it needs, filtering by a target set of NCT ids (set
membership, O(1)), so even multi-GB tables are processed in a single
streamed pass without loading them into memory.

Primary extraction path = `outcome_analyses` (pre-computed effect + 95% CI
per trial), which maps cleanly onto a metaforge pairwise/NMA dataset.
Secondary path = `outcome_counts` (per-arm event counts) for 2x2 tables.

Snapshot discovery order (never hardcode one drive):
  1. env var TRIALFORGE_AACT
  2. explicit `root=` argument
  3. candidate roots on F:/D:/C: (configurable list)
Fails closed with a clear message if no snapshot with the expected tables
is found (lessons.md: validate snapshot, don't assume a drive).
"""
from __future__ import annotations
import os
from pathlib import Path
from typing import Iterable, Optional

# AACT param_type -> (metaforge measure, is_ratio)
_RATIO_KEYS = [
    ("hazard", "HR"), ("odds", "OR"), ("risk ratio", "RR"),
    ("rate ratio", "RR"), ("relative risk", "RR"),
    ("geometric mean ratio", "RR"), ("gmt ratio", "RR"), ("gmr", "RR"),
]
_DIFF_KEYS = [
    ("risk difference", "RD"), ("difference in proportion", "RD"),
    ("mean difference", "MD"), ("ls mean", "MD"), ("lsmean", "MD"),
    ("least square", "MD"), ("least squares", "MD"),
    ("adjusted mean difference", "MD"), ("difference in percentage", "RD"),
]

CANDIDATE_ROOTS = [
    r"F:\AACT-storage\AACT", r"D:\AACT-storage\AACT", r"C:\AACT-storage\AACT",
    r"F:\AACT", r"D:\AACT", r"C:\AACT",
]
REQUIRED = ("outcome_analyses.txt", "outcomes.txt", "interventions.txt",
            "conditions.txt")


def classify_param_type(pt: str):
    """Return (measure, is_ratio) or (None, None) if not poolable."""
    if not pt:
        return None, None
    low = pt.lower()
    for key, meas in _RATIO_KEYS:
        if key in low:
            return meas, True
    for key, meas in _DIFF_KEYS:
        if key in low:
            return meas, False
    return None, None


def find_snapshot(root: Optional[str] = None) -> Path:
    cands = []
    if os.environ.get("TRIALFORGE_AACT"):
        cands.append(os.environ["TRIALFORGE_AACT"])
    if root:
        cands.append(root)
    cands += CANDIDATE_ROOTS
    for c in cands:
        base = Path(c)
        if not base.exists():
            continue
        # snapshot may be base itself or base/<date>/
        candidates = [base]
        if base.is_dir():
            candidates += sorted([p for p in base.iterdir() if p.is_dir()], reverse=True)
        for cand in candidates:
            if all((cand / t).exists() for t in REQUIRED):
                return cand
    raise FileNotFoundError(
        "No AACT snapshot found. Set the TRIALFORGE_AACT environment variable "
        "to a folder containing outcome_analyses.txt etc., or pass root=. "
        f"Looked in: {cands}")


def _scan(path: Path, ncol: int, nct_field: int, ncts: set, stats: dict = None):
    """Yield split rows of a pipe-delimited AACT table where field[nct_field]
    is in `ncts` and the row has the expected column count (guards against
    rows broken by embedded newlines in text fields). If `stats` is given,
    counts column-count-mismatched rows that mention an NCT into
    stats['skipped'] so silent data loss is visible to the caller."""
    with path.open("r", encoding="utf-8", errors="replace") as f:
        f.readline()  # header
        for line in f:
            # cheap pre-filter: must contain an NCT id at all
            if "NCT" not in line:
                continue
            fields = line.rstrip("\n").split("|")
            if len(fields) != ncol:
                if stats is not None:
                    stats["skipped"] = stats.get("skipped", 0) + 1
                continue
            if fields[nct_field] in ncts:
                yield fields


def _ncol(path: Path) -> int:
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return len(f.readline().rstrip("\n").split("|"))


class AACT:
    def __init__(self, root: Optional[str] = None):
        self.dir = find_snapshot(root)

    def _t(self, name):
        return self.dir / f"{name}.txt"

    # ---- discovery ------------------------------------------------------
    def find_trials(self, drug: Optional[str] = None,
                    condition: Optional[str] = None,
                    limit: int = 500) -> list:
        """Return NCT ids whose interventions match `drug` (substring, any
        drug/biological type) AND whose conditions match `condition`."""
        drug_ncts = None
        if drug:
            drug_ncts = set()
            dl = drug.lower()
            p = self._t("interventions")
            ncol = _ncol(p)
            with p.open("r", encoding="utf-8", errors="replace") as f:
                f.readline()
                for line in f:
                    if "NCT" not in line:
                        continue
                    fields = line.rstrip("\n").split("|")
                    if len(fields) != ncol:
                        continue
                    # id|nct_id|intervention_type|name|description
                    if dl in fields[3].lower():
                        drug_ncts.add(fields[1])
        cond_ncts = None
        if condition:
            cond_ncts = set()
            cl = condition.lower()
            p = self._t("conditions")
            ncol = _ncol(p)
            with p.open("r", encoding="utf-8", errors="replace") as f:
                f.readline()
                for line in f:
                    if "NCT" not in line:
                        continue
                    fields = line.rstrip("\n").split("|")
                    if len(fields) != ncol:
                        continue
                    # id|nct_id|name|downcase_name
                    if cl in fields[2].lower():
                        cond_ncts.add(fields[1])
        if drug_ncts is not None and cond_ncts is not None:
            ncts = drug_ncts & cond_ncts
        else:
            ncts = drug_ncts or cond_ncts or set()
        return sorted(ncts)[:limit]

    # ---- outcome metadata ----------------------------------------------
    def _primary_outcome_ids(self, ncts: set) -> set:
        """outcome_id values whose outcome_type is Primary, for target NCTs."""
        p = self._t("outcomes")
        ncol = _ncol(p)
        prim = set()
        for fields in _scan(p, ncol, 1, ncts):
            # id|nct_id|outcome_type|...
            if fields[2].strip().lower() == "primary":
                prim.add(fields[0])
        return prim

    # ---- extraction: precomputed effects -------------------------------
    def extract_effects(self, ncts: Iterable[str], prefer_primary: bool = True,
                        force_measure: Optional[str] = None) -> dict:
        """Pull one effect estimate per NCT from outcome_analyses.

        Returns {"measure": <dominant>, "studies": [...], "report": {...}}.
        Each study: {name, nct, effect, ci_low, ci_high}. When CIs are
        missing the row is skipped (we need them to weight the study).
        """
        ncts = set(ncts)
        if not ncts:
            return {"measure": None, "studies": [], "report": {"reason": "no NCTs"}}
        primary = self._primary_outcome_ids(ncts) if prefer_primary else set()

        p = self._t("outcome_analyses")
        ncol = _ncol(p)
        # header indices
        # id|nct_id|outcome_id|...|param_type(5)|param_value(6)|...|ci_lower(13)|ci_upper(14)
        scan_stats = {"skipped": 0}
        rows_by_nct = {}
        for fields in _scan(p, ncol, 1, ncts, stats=scan_stats):
            nct = fields[1]
            outcome_id = fields[2]
            param_type = fields[5]
            param_value = fields[6]
            ci_lo = fields[13]
            ci_hi = fields[14]
            measure, is_ratio = classify_param_type(param_type)
            if measure is None:
                continue
            if force_measure and measure != force_measure:
                continue
            try:
                est = float(param_value)
                lo = float(ci_lo)
                hi = float(ci_hi)
            except (ValueError, TypeError):
                continue
            if is_ratio and (est <= 0 or lo <= 0 or hi <= 0):
                continue
            if hi <= lo:
                continue
            is_primary = outcome_id in primary
            rec = {"nct": nct, "measure": measure, "is_ratio": is_ratio,
                   "effect": est, "ci_low": lo, "ci_high": hi,
                   "primary": is_primary}
            rows_by_nct.setdefault(nct, []).append(rec)

        # choose one analysis per NCT. Pick the dominant measure preferring
        # PRIMARY-outcome analyses (align on the primary effect type across
        # trials, not whatever outcome is most numerous).
        from collections import Counter
        measure_counter = Counter()
        primary_counter = Counter()
        for recs in rows_by_nct.values():
            for r in recs:
                measure_counter[r["measure"]] += 1
                if r["primary"]:
                    primary_counter[r["measure"]] += 1
        if not measure_counter:
            return {"measure": None, "studies": [],
                    "report": {"reason": "no poolable effect estimates found",
                               "n_ncts": len(ncts)}}
        if force_measure:
            dominant = force_measure
        elif primary_counter:
            dominant = primary_counter.most_common(1)[0][0]
        else:
            dominant = measure_counter.most_common(1)[0][0]

        studies = []
        for nct, recs in rows_by_nct.items():
            pool = [r for r in recs if r["measure"] == dominant]
            if not pool:
                continue
            pool.sort(key=lambda r: (not r["primary"],))  # primary first
            r = pool[0]
            studies.append({"name": nct, "nct": nct,
                            "effect": r["effect"], "ci_low": r["ci_low"],
                            "ci_high": r["ci_high"]})
        studies.sort(key=lambda s: s["nct"])
        return {
            "measure": dominant,
            "studies": studies,
            "report": {
                "n_ncts_queried": len(ncts),
                "n_with_effects": len(studies),
                "measure_distribution": dict(measure_counter),
                "malformed_rows_skipped": scan_stats["skipped"],
                "source": "AACT outcome_analyses",
                "snapshot": str(self.dir),
            },
        }
