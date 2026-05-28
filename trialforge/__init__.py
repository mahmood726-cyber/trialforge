"""trialforge — AACT-powered, advanced meta-analysis engine.

Builds on the bundled metaforge engine (pairwise / proportion / NMA /
dose-response) and adds:
  * AACT (ClinicalTrials.gov) data ingestion        (trialforge.aact)
  * advanced methods: publication-bias suite, influence diagnostics,
    meta-regression, subgroup, cumulative, rare-events  (trialforge.advanced)
  * NMA inconsistency via Bucher closed loops          (trialforge.nodesplit)
  * Copas selection-model publication-bias sensitivity (trialforge.copas)
  * diagnostic test accuracy (bivariate DL + SROC)     (trialforge.dta)
  * additive component network meta-analysis           (trialforge.cnma)
  * limit meta-analysis (Rucker small-study adjustment)(trialforge.limitma)
  * trial sequential analysis (RIS + OBF boundary)     (trialforge.tsa)
  * E-value for unmeasured confounding                 (trialforge.evalue)
  * p-curve evidential-value analysis                  (trialforge.pcurve)
  * GOSH subset-heterogeneity diagnostic               (trialforge.gosh)
  * CINeMA-style confidence-in-NMA rating              (trialforge.cinema)
  * rare-events binomial/hypergeometric-normal GLMM    (trialforge.glmm)
  * bivariate multi-outcome borrowing of strength      (trialforge.multivariate)
  * RMST / survival meta-analysis                      (trialforge.survival)
  * GRADE certainty + Summary-of-Findings automation   (trialforge.grade)

The statistics modules (common, pairwise, proportions, nma, doseresponse,
plots, report) are the metaforge engine, bundled so the kit is
self-contained and offline. The advanced methods are re-implemented
(pure stdlib) from the allmeta suite.
"""
from . import common, pairwise, proportions, nma, doseresponse  # noqa: F401
from . import advanced, nodesplit, tfreport, copas, dta, cnma  # noqa: F401
from . import limitma, tsa, evalue, pcurve, gosh, cinema  # noqa: F401
from . import glmm, multivariate, survival, grade, checks  # noqa: F401

__version__ = "1.6.0"
