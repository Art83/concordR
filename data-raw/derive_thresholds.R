# =============================================================================
# derive_thresholds.R
# Defines atlas_thresholds, the numeric cut-offs used by the triage() gates.
# =============================================================================
# Sourced by prep_atlas_data.R after all atlas objects exist, so the objects
# concordance_atlas, tissue_atlas, hpa_tissue, and localisation_flags are in
# scope here and can be interrogated when deriving a value.
#
# =============================================================================

atlas_thresholds <- list(
  protein_tau_specific = 0.85,   # gate input, IHC panel tau, antimode of a bimodal distribution
  ihc_min_detected = 3L,         # below this many IHC-detected tissues, specificity is sparse or na
  bulk_breadth_ihc_uninformative = 0.5, # bulk-detected share of covered tissues at or above which a sparse-IHC gene reads na
  tau_specific = 0.85,           # evidence only, RNA log-mean tau
  protein_breadth_broad = 0.75,  # evidence only, coverage-confounded
  cross_platform_agree_max = 0.30,
  n_platforms_min = 2L
)

# ---- Sanity ------------------------------------------------------------------
stopifnot(
  is.numeric(atlas_thresholds$protein_tau_specific),
  atlas_thresholds$protein_tau_specific > 0,
  atlas_thresholds$protein_tau_specific < 1,
  atlas_thresholds$ihc_min_detected >= 2L,
  atlas_thresholds$bulk_breadth_ihc_uninformative > 0,
  atlas_thresholds$bulk_breadth_ihc_uninformative <= 1,
  is.numeric(atlas_thresholds$tau_specific),
  atlas_thresholds$tau_specific > 0, atlas_thresholds$tau_specific < 1,
  atlas_thresholds$protein_breadth_broad > 0,
  atlas_thresholds$protein_breadth_broad <= 1,
  atlas_thresholds$cross_platform_agree_max >= 0,
  atlas_thresholds$cross_platform_agree_max <= 1,
  atlas_thresholds$n_platforms_min >= 1L
)

cat("\n==== atlas_thresholds ====\n")
str(atlas_thresholds)
