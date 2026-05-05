# =============================================================================
# prep_atlas_data.R
# Assembles the concordance_atlas internal dataset for concordR
# from the existing pipeline outputs.
#
# Run once, then: usethis::use_data(concordance_atlas, internal = TRUE)
#
# Inputs (adjust paths if your project root differs):
#   - results/ml/gene_discordance_profiles.csv
#   - results/ml/df_ml_v4_*_mincells50.parquet   (picks the latest match)
#   - results/ml/four_way_decomposition.csv      (optional: mechanistic tier)
#
# Output bundled columns:
#   gene_symbol, ensembl_id, gene_class, protein_confidence, detection_rate,
#   suppression_rate, mean_rna_rank,
#   n_tissues_concordant, n_tissues_suppressed, n_tissues,
#   mechanism_tier,
#   is_secreted, is_membrane, is_nuclear, is_cytoplasmic, is_mitochondrial,
#   is_er_golgi, is_cytoskeleton, is_projection, is_multilocal,
#   tau, n_tissues_for_tau, specificity_class
# =============================================================================

library(data.table)
library(arrow)

# ---- Paths -----------------------------------------------------------------
PROJECT_ROOT  <- "D:/Bioinformatics_projects/hpa_to_tabula"
ATLAS_CSV     <- file.path(PROJECT_ROOT, "results/ml/gene_discordance_profiles.csv")
TIER_CSV      <- file.path(PROJECT_ROOT, "results/ml/four_way_decomposition.csv")
PARQUET_DIR   <- file.path(PROJECT_ROOT, "results/ml")

# Auto-pick the latest matching parquet so the script doesn't go stale on rerun
parquet_candidates <- list.files(
  PARQUET_DIR,
  pattern    = "^df_ml_v4_.*_mincells50\\.parquet$",
  full.names = TRUE
)
stopifnot(length(parquet_candidates) >= 1)
PARQUET_PATH <- parquet_candidates[which.max(file.info(parquet_candidates)$mtime)]
cat("Using parquet:", basename(PARQUET_PATH), "\n")

# ---- 1. Gene discordance profiles (core atlas) -----------------------------
atlas <- fread(ATLAS_CSV)

required <- c("gene_symbol", "ensembl_id", "gene_class",
              "protein_confidence", "detection_rate", "mean_rna_rank")
missing_required <- setdiff(required, names(atlas))
if (length(missing_required) > 0) {
  stop("Missing required columns in atlas CSV: ",
       paste(missing_required, collapse = ", "))
}

# Optional columns we want if present (some older atlas builds may not have them)
optional_atlas_cols <- c("suppression_rate",
                         "n_tissues_concordant",
                         "n_tissues_suppressed")
missing_optional <- setdiff(optional_atlas_cols, names(atlas))
if (length(missing_optional) > 0) {
  warning("Atlas CSV is missing optional columns: ",
          paste(missing_optional, collapse = ", "),
          ". They will be computed from the parquet if possible, or left NA.")
}

# Deduplicate by gene symbol (keep first)
atlas <- atlas[!duplicated(atlas[["gene_symbol"]])]
cat("Atlas rows after dedup:", nrow(atlas), "\n")

# ---- 2. Load parquet (robust to tibble/data.frame/data.table) --------------
# arrow::read_parquet returns a tibble; force to plain data.frame then data.table
df_raw <- as.data.frame(arrow::read_parquet(PARQUET_PATH))
df     <- as.data.table(df_raw)
rm(df_raw)
cat("Parquet rows:", nrow(df), " cols:", ncol(df), "\n")
stopifnot("data.table" %in% class(df))

# Sanity: required columns for downstream steps
stopifnot(all(c("gene_symbol", "ts_tissue", "rna_rank") %in% names(df)))

# ---- 3. Gene features (subcellular localisation flags) ---------------------
loc_cols <- c("is_secreted", "is_membrane", "is_nuclear",
              "is_cytoplasmic", "is_mitochondrial", "is_er_golgi",
              "is_cytoskeleton", "is_projection", "is_multilocal")
loc_cols <- intersect(loc_cols, names(df))
cat("Localisation columns found:", paste(loc_cols, collapse = ", "), "\n")

# First non-NA value per gene for each localisation flag
gene_feat <- df[,
  lapply(.SD, function(x) {
    v <- x[!is.na(x)]
    if (length(v) == 0) NA_integer_ else as.integer(v[1])
  }),
  by    = gene_symbol,
  .SDcols = loc_cols
]

# ---- 4. Tissue count per gene (from observations) --------------------------
n_tissues <- df[,
  list(n_tissues = uniqueN(ts_tissue)),
  by = gene_symbol
]

# ---- 5. Tissue specificity (Yanai tau) on rna_rank -------------------------
# Aggregate to one value per (gene, tissue), averaging across cell types.
# Using rna_rank (not rna_mean) because rna_rank was computed per-source in
# script 02/04 and is therefore comparable across Tabula and Brain atlases.
gene_tissue <- df[,
  list(rna_tissue = mean(rna_rank, na.rm = TRUE)),
  by = list(gene_symbol, ts_tissue)
]

tau_dt <- gene_tissue[,
  {
    n <- .N
    if (n < 10) {
      list(tau = NA_real_, n_tissues_for_tau = n)
    } else {
      mx <- max(rna_tissue, na.rm = TRUE)
      if (!is.finite(mx) || mx == 0) {
        list(tau = NA_real_, n_tissues_for_tau = n)
      } else {
        x <- rna_tissue / mx
        list(tau               = sum(1 - x, na.rm = TRUE) / (n - 1),
             n_tissues_for_tau = n)
      }
    }
  },
  by = gene_symbol
]

tau_dt[, specificity_class := fcase(
  is.na(tau),  "insufficient_data",
  tau >= 0.55, "tissue_specific",
  tau >= 0.40, "tissue_enriched",
  tau >= 0.20, "mixed",
  default    = "broadly_expressed"
)]

cat("tau computed for ", sum(!is.na(tau_dt$tau)), " / ", nrow(tau_dt),
    " genes\n", sep = "")

# ---- 6. Mechanistic tier (optional) ----------------------------------------
if (file.exists(TIER_CSV)) {
  tiers <- fread(TIER_CSV)
  if (all(c("ensembl_id", "mechanism_tier") %in% names(tiers))) {
    tiers <- unique(tiers[, list(ensembl_id, mechanism_tier)])
    atlas <- merge(atlas, tiers, by = "ensembl_id", all.x = TRUE)
    cat("Mechanism tier merged for ",
        sum(!is.na(atlas$mechanism_tier)), " genes\n", sep = "")
  } else {
    warning("Tier CSV present but lacks expected columns ",
            "(ensembl_id, mechanism_tier); skipping.")
  }
} else {
  cat("No mechanism tier CSV found; skipping.\n")
}

# ---- 6b. Gene × tissue collapsed (cell-types collapsed within tissue) ----
# For each (gene, tissue): take the maximum RNA rank across cell types
# (the question is "is this gene present in this tissue at all?"), the
# fraction of cell types in the tissue where protein is detected, and
# the maximum protein score across cell types.
tissue_atlas <- df[, list(
  rna_rank_max         = max(rna_rank,     na.rm = TRUE),
  rna_mean_tissue      = mean(rna_mean,    na.rm = TRUE),
  detect_fraction      = mean(y_bin,       na.rm = TRUE),
  protein_max          = max(protein_score, na.rm = TRUE),
  n_cell_types         = uniqueN(ts_cell_type),
  is_brain             = max(is_brain,     na.rm = TRUE)
), by = list(gene_symbol, ts_tissue)]

# Per-(gene, tissue) class label.
# Thresholds mirror the global atlas call but use the collapsed
# tissue-level statistics. `detect_fraction` is the share of cell
# types in the tissue with protein detection (not raw detection
# probability). The "tissue gate" used by triage() asks: is this gene
# expressed at protein level in this tissue at all?
tissue_atlas[, tissue_class := data.table::fcase(
  rna_mean_tissue == 0 & detect_fraction == 0,            "low_expression",
  rna_rank_max <  0.20 & detect_fraction <  0.10,         "low_expression",
  rna_rank_max >= 0.50 & detect_fraction >= 0.50,         "concordant_in_tissue",
  rna_rank_max >= 0.50 & detect_fraction <  0.20,         "suppressed_in_tissue",
  default = "variable_in_tissue"
)]

# Sanity check before saving
cat("tissue_atlas: ", nrow(tissue_atlas), " (gene, tissue) rows\n", sep = "")
cat("  unique genes:   ", uniqueN(tissue_atlas$gene_symbol), "\n", sep = "")
cat("  unique tissues: ", uniqueN(tissue_atlas$ts_tissue),   "\n", sep = "")
cat("  size (MB):      ",
    format(object.size(tissue_atlas), units = "Mb"), "\n", sep = "")
print(table(tissue_atlas$tissue_class, useNA = "ifany"))



# ---- 7. Merge everything on gene_symbol ------------------------------------
concordance_atlas <- atlas
concordance_atlas <- merge(concordance_atlas, gene_feat, by = "gene_symbol", all.x = TRUE)
concordance_atlas <- merge(concordance_atlas, n_tissues, by = "gene_symbol", all.x = TRUE)
concordance_atlas <- merge(concordance_atlas, tau_dt,    by = "gene_symbol", all.x = TRUE)

# ---- 8. Fill NAs for binary localisation flags -----------------------------
for (col in loc_cols) {
  if (col %in% names(concordance_atlas)) {
    concordance_atlas[[col]][is.na(concordance_atlas[[col]])] <- 0L
    concordance_atlas[[col]] <- as.integer(concordance_atlas[[col]])
  }
}

# ---- 9. Final column ordering (optional, for readability) ------------------
preferred_order <- c(
  "gene_symbol", "ensembl_id",
  "gene_class",  "protein_confidence", "detection_rate", "suppression_rate",
  "mean_rna_rank",
  "n_tissues", "n_tissues_concordant", "n_tissues_suppressed",
  "mechanism_tier",
  loc_cols,
  "tau", "n_tissues_for_tau", "specificity_class"
)
preferred_order <- intersect(preferred_order, names(concordance_atlas))
remaining       <- setdiff(names(concordance_atlas), preferred_order)
setcolorder(concordance_atlas, c(preferred_order, remaining))

# ---- 10. Convert to data.frame for R package -------------------------------
concordance_atlas <- as.data.frame(concordance_atlas)

# ---- 11. Report -------------------------------------------------------------
cat("\n==== concordance_atlas ====\n")
cat("Dimensions: ", nrow(concordance_atlas), " genes x ",
    ncol(concordance_atlas), " columns\n", sep = "")
cat("Columns: ", paste(names(concordance_atlas), collapse = ", "), "\n", sep = "")

cat("\nClass distribution:\n")
print(table(concordance_atlas$gene_class, useNA = "ifany"))

if ("specificity_class" %in% names(concordance_atlas)) {
  cat("\nSpecificity class distribution:\n")
  print(table(concordance_atlas$specificity_class, useNA = "ifany"))
}

if ("mechanism_tier" %in% names(concordance_atlas)) {
  cat("\nMechanism tier distribution:\n")
  print(table(concordance_atlas$mechanism_tier, useNA = "ifany"))
}

