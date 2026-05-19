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
library(dplyr)
library(tidyr)

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

# ---- 6b. Gene x tissue collapsed (cell-types collapsed within tissue) ------
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
), by = list(gene_symbol, ensembl_id, ts_tissue)]

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


# ===========================================================================
# External platform data: PaxDB, GTEx protein, GTEx RNA
# ===========================================================================
# NOTE: tissue_atlas is already in memory from section 6b above.
#       Previous versions had a `load("R/sysdata.rda")` here which
#       silently overwrote the freshly-built tissue_atlas. Removed.


# --- PaxDB (merged via ENSP -> ENSG, not symbol) ---------------------------
#
# PaxDb indexes proteins by Ensembl protein ID (ENSP). Previous versions
# merged by gene symbol via an intermediate map, which silently dropped
# genes where PaxDb's symbol annotation diverged from HGNC (e.g. LILRA4:
# ENSP00000291759 in PaxDb, ENSG00000239961 in the atlas).
#
# Fix: read the ENSP from PaxDb, map to ENSG via biomaRt, merge on
# ensembl_id. The ENSP is retained as paxdb_ensp for traceability.

PAXDB_TISSUE_MAP <- list(
  "cerebral-cortex" = "Cerebral cortex",
  "heart"           = "Heart",
  "kidney"          = "Kidney",
  "liver"           = "Liver",
  "lung"            = "Lung",
  "lymph_node"      = "Lymph Node",
  "pancreas"        = "Pancreas",
  "skin"            = "Skin",
  "spleen"          = "Spleen",
  "stomach"         = "Stomach",
  "testis"          = "Testis"
)

PAXDB_DIR <- file.path(PROJECT_ROOT, "external_data/paxdb")

paxdb_raw <- lapply(names(PAXDB_TISSUE_MAP), function(name) {
  df <- read.csv(file.path(PAXDB_DIR, paste0(name, ".txt")),
                 skip = 11, sep = "\t", header = FALSE,
                 stringsAsFactors = FALSE, comment.char = "")
  # PaxDb TSV: col1 = internal_id, col2 = string_external_id, col3 = abundance
  # string_external_id format: "9606.ENSP00000XXXXXX"
  colnames(df) <- c("paxdb_internal_id", "paxdb_string_id", "abundance")
  df$ts_tissue <- PAXDB_TISSUE_MAP[[name]]
  df
})
paxdb_raw <- do.call(rbind, paxdb_raw)

# Extract bare ENSP from "9606.ENSP..." format
paxdb_raw$ensp <- sub("^9606\\.", "", paxdb_raw$paxdb_string_id)

cat(sprintf("PaxDb raw: %d rows, %d unique ENSP\n",
            nrow(paxdb_raw), length(unique(paxdb_raw$ensp))))

# --- ENSP -> ENSG mapping via biomaRt (cached to avoid repeated queries) ---
ENSP_MAP_CACHE <- file.path(PROJECT_ROOT, "external_data/ensp_to_ensg_biomart.txt")

if (file.exists(ENSP_MAP_CACHE)) {
  cat("Loading cached ENSP -> ENSG map from:", ENSP_MAP_CACHE, "\n")
  ensp_map <- read.csv(ENSP_MAP_CACHE,sep="\t", stringsAsFactors = FALSE)
}

colnames(ensp_map) <- c("ensembl_id", "ensp", "hgnc_symbol")

# Deduplicate: keep one ENSG per ENSP (multiple isoforms -> same gene)
ensp_map <- ensp_map[!duplicated(ensp_map$ensp), ]

cat(sprintf("ENSP map: %d ENSP -> %d unique ENSG\n",
            nrow(ensp_map), length(unique(ensp_map$ensembl_id))))

# Join PaxDb to ENSG, then collapse isoforms: max abundance per (gene, tissue)
paxdb <- paxdb_raw %>%
  inner_join(ensp_map[, c("ensp", "ensembl_id")], by = "ensp") %>%
  group_by(ensembl_id, ts_tissue) %>%
  summarise(
    abundance  = max(abundance, na.rm = TRUE),
    paxdb_ensp = paste(unique(ensp), collapse = ";"),  # traceability
    .groups    = "drop"
  )

# Diagnostic: how many atlas genes now have PaxDb coverage?
n_atlas_with_paxdb <- sum(unique(tissue_atlas$ensembl_id) %in% paxdb$ensembl_id)
cat(sprintf("PaxDb after ENSP mapping: %d (gene, tissue) rows, %d atlas genes covered\n",
            nrow(paxdb), n_atlas_with_paxdb))


# --- GTEx shared columns ---
cols <- c("Brain.Cortex", "Breast", "Heart.Atrial", "Heart.Ventricle",
          "Liver", "Lung", "Minor.Salivary", "Muscle.Skeletal", "Pancreas",
          "Prostate", "Skin.Unexpo", "Small.Intestine", "Spleen", "Stomach",
          "Testis", "Uterus")

rename_vec <- c(
  "Brain.Cortex"     = "Cerebral cortex",
  "Breast"           = "Mammary",
  "Minor.Salivary"   = "Salivary Gland",
  "Muscle.Skeletal"  = "Muscle",
  "Skin.Unexpo"      = "Skin",
  "Small.Intestine"  = "Small Intestine"
)

EXTERNAL_DIR <- file.path(PROJECT_ROOT, "external_data")

# =========================================================
# GTEx Medians (already numeric, no flags)
# =========================================================
gtex_medians <- read.csv(file.path(EXTERNAL_DIR, "gtex_protein_tissue_median.csv"))

gtex_medians <- gtex_medians %>%
  dplyr::select(all_of(c("gene.id", cols))) %>%
  rowwise() %>%
  mutate(Heart = median(c(Heart.Atrial, Heart.Ventricle), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(-Heart.Ventricle, -Heart.Atrial) %>%
  rename(ensembl_id = gene.id) %>%
  tidyr::pivot_longer(-ensembl_id, names_to = "ts_tissue", values_to = "abundance_gtex") %>%
  mutate(ts_tissue = recode(ts_tissue, !!!rename_vec))


# =========================================================
# GTEx Protein TS scores
# =========================================================
gtex_ts_raw <- read.csv(file.path(EXTERNAL_DIR, "gtex_protein_ts_score.csv"))

gtex_ts <- gtex_ts_raw[, c("ensembl_id", cols)] %>%
  tidyr::pivot_longer(-ensembl_id, names_to = "ts_tissue", values_to = "raw") %>%
  mutate(
    ts_protein_gtex = as.numeric(sapply(strsplit(raw, ";", fixed = TRUE), `[`, 1)),
    flag = sapply(strsplit(raw, ";", fixed = TRUE), `[`, 2)
  ) %>%
  filter(is.na(flag) | flag != "NA_one_rep_in_raw") %>%
  dplyr::select(-raw, -flag) %>%
  mutate(ts_tissue = recode(ts_tissue, !!!rename_vec,
                            "Heart.Atrial" = "Heart", "Heart.Ventricle" = "Heart")) %>%
  group_by(ensembl_id, ts_tissue) %>%
  summarise(ts_protein_gtex = median(ts_protein_gtex, na.rm = TRUE), .groups = "drop")


# =========================================================
# GTEx RNA TS scores
# =========================================================
gtex_rna_raw <- read.csv(file.path(EXTERNAL_DIR, "gtex_rna_ts_score.csv"))

gtex_rna_ts <- gtex_rna_raw[, c("ensembl_id", cols)] %>%
  tidyr::pivot_longer(-ensembl_id, names_to = "ts_tissue", values_to = "raw") %>%
  mutate(
    ts_rna_gtex = as.numeric(sapply(strsplit(raw, "[; ]"), `[`, 1)),
    flag = sapply(strsplit(raw, "[; ]"), `[`, 2)
  ) %>%
  filter(is.na(flag) | flag != "NA_raw_tpm_less_1") %>%
  filter(is.na(flag) | flag != "all_tissues_tpm_less_1") %>%
  dplyr::select(-raw, -flag) %>%
  mutate(ts_tissue = recode(ts_tissue, !!!rename_vec,
                            "Heart.Atrial" = "Heart", "Heart.Ventricle" = "Heart")) %>%
  group_by(ensembl_id, ts_tissue) %>%
  summarise(ts_rna_gtex = median(ts_rna_gtex, na.rm = TRUE), .groups = "drop")


# ===========================================================================
# Augment tissue_atlas with external platforms
# ===========================================================================
# NOTE: tissue_atlas was built in section 6b above.
#       No load("R/sysdata.rda") — that was a bug in previous versions
#       which silently overwrote the freshly-built tissue_atlas.

tissue_atlas <- tissue_atlas %>%
  left_join(gtex_medians, by = c("ensembl_id", "ts_tissue")) %>%
  left_join(gtex_ts,      by = c("ensembl_id", "ts_tissue")) %>%
  left_join(gtex_rna_ts,  by = c("ensembl_id", "ts_tissue")) %>%
  left_join(paxdb,         by = c("ensembl_id", "ts_tissue")) %>%
  mutate(
    ms_detected    = !is.na(ts_protein_gtex),
    paxdb_detected = !is.na(abundance),
    ihc_detected   = protein_max >= 1,
    protein_evidence_sources = ms_detected + paxdb_detected + ihc_detected,
    ihc_ms_agreement = case_when(
      ihc_detected & ms_detected  ~ "both_detected",
      ihc_detected & !ms_detected ~ "ihc_only",
      !ihc_detected & ms_detected ~ "ms_only",
      TRUE                        ~ "both_absent"
    )
  ) %>%
  mutate(
    # GTEx internal concordance
    gtex_concordance = case_when(
      is.na(ts_rna_gtex) | is.na(ts_protein_gtex) ~ NA_character_,
      ts_rna_gtex > 0 & ts_protein_gtex > 0       ~ "gtex_concordant",
      ts_rna_gtex > 0 & ts_protein_gtex <= 0      ~ "gtex_suppressed",
      ts_rna_gtex <= 0 & ts_protein_gtex > 0      ~ "gtex_stabilised",
      TRUE                                         ~ "gtex_low_both"
    ),
    # HPA vs GTEx agreement on direction
    platforms_agree = case_when(
      is.na(gtex_concordance) ~ NA,
      tissue_class == "suppressed_in_tissue" &
        gtex_concordance == "gtex_suppressed"       ~ TRUE,
      tissue_class == "concordant_in_tissue" &
        gtex_concordance == "gtex_concordant"        ~ TRUE,
      tissue_class %in% c("suppressed_in_tissue", "concordant_in_tissue") ~ FALSE,
      TRUE                                           ~ NA  # variable/low: can't compare
    ),
    # Protein detection consensus (all three platforms)
    protein_detection_consensus = case_when(
      protein_evidence_sources >= 2 ~ "multi_source",
      protein_evidence_sources == 1 ~ "single_source",
      TRUE                          ~ "no_detection"
    )
  )


# ===========================================================================
# Gene-level summaries from external platforms
# ===========================================================================
gtex_gene_summary <- tissue_atlas %>%
  group_by(gene_symbol) %>%
  summarise(
    gtex_n_tissues_detected    = sum(ms_detected, na.rm = TRUE),
    paxdb_n_tissues_detected   = sum(paxdb_detected, na.rm = TRUE),
    max_protein_evidence       = max(protein_evidence_sources, na.rm = TRUE),
    cross_platform_agree_rate  = mean(platforms_agree, na.rm = TRUE),
    gtex_suppression_rate      = mean(gtex_concordance == "gtex_suppressed", na.rm = TRUE),
    .groups = "drop"
  )

concordance_atlas <- concordance_atlas %>%
  left_join(gtex_gene_summary, by = "gene_symbol")


# ===========================================================================
# Brain-specific summary
# ===========================================================================
brain_summary <- tissue_atlas %>%
  filter(is_brain == 1, rna_rank_max > 0.70) %>%
  group_by(gene_symbol, ensembl_id) %>%
  summarise(
    brain_suppression_rate    = mean(tissue_class == "suppressed_in_tissue"),
    brain_n_tissues           = n(),
    brain_ms_detected         = any(ms_detected, na.rm = TRUE),
    brain_platforms_agree     = mean(platforms_agree, na.rm = TRUE),
    brain_protein_consensus   = max(protein_evidence_sources, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(brain_class = case_when(
    brain_suppression_rate > 0.80 ~ "brain_suppressed",
    brain_suppression_rate < 0.20 ~ "brain_concordant",
    TRUE                          ~ "brain_variable"
  ))


# ===========================================================================
# Save sysdata
# ===========================================================================
concordance_atlas <- as.data.frame(concordance_atlas)
tissue_atlas      <- as.data.frame(tissue_atlas)
brain_summary     <- as.data.frame(brain_summary)

cat("\n==== Final dimensions ====\n")
cat("  concordance_atlas:", nrow(concordance_atlas), "genes x",
    ncol(concordance_atlas), "cols\n")
cat("  tissue_atlas:     ", nrow(tissue_atlas), "(gene, tissue) rows x",
    ncol(tissue_atlas), "cols\n")
cat("  brain_summary:    ", nrow(brain_summary), "genes x",
    ncol(brain_summary), "cols\n")

# Sanity: confirm PaxDb coverage improved
cat("\n==== PaxDb coverage check ====\n")
cat("  tissue_atlas rows with paxdb_detected = TRUE:",
    sum(tissue_atlas$paxdb_detected, na.rm = TRUE), "\n")
cat("  concordance_atlas genes with paxdb_n_tissues > 0:",
    sum(concordance_atlas$paxdb_n_tissues_detected > 0, na.rm = TRUE), "\n")

usethis::use_data(concordance_atlas, tissue_atlas, brain_summary,
                  internal = TRUE, overwrite = TRUE)
