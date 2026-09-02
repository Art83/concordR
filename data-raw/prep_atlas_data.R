# =============================================================================
# prep_atlas_data.R
# Assembles the concordR internal datasets from the existing pipeline outputs.
# =============================================================================

library(data.table)
library(arrow)
library(dplyr)
library(tidyr)

# ---- Paths -----------------------------------------------------------------
# These are all outputs from companion paper on rna to protein concordance
PROJECT_ROOT  <- "D:/Bioinformatics_projects/hpa_to_tabula"
ATLAS_CSV <- file.path(PROJECT_ROOT, "results/ml/gene_discordance_profiles.csv")
TIER_CSV <- file.path(PROJECT_ROOT, "results/tables/three_way_decomposition.csv")
PARQUET_DIR <- file.path(PROJECT_ROOT, "results/ml")
EXTERNAL_DIR <- file.path(PROJECT_ROOT, "external_data")
LOC_CSV <- "D:/Bioinformatics_projects/localisation_uniprot.csv"

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
if (length(missing_required) > 0)
  stop("Missing required columns in atlas CSV: ",
       paste(missing_required, collapse = ", "))

optional_atlas_cols <- c("suppression_rate",
                         "n_concordant_pos",
                         "n_suppressed")
missing_optional <- setdiff(optional_atlas_cols, names(atlas))
if (length(missing_optional) > 0)
  warning("Atlas CSV is missing optional columns: ",
          paste(missing_optional, collapse = ", "), ".")

atlas <- atlas[!duplicated(atlas[["gene_symbol"]])]
cat("Atlas rows after dedup:", nrow(atlas), "\n")

# ---- 2. Load parquet -------------------------------------------------------
df_raw <- as.data.frame(arrow::read_parquet(PARQUET_PATH))
df <- as.data.table(df_raw)
rm(df_raw)
cat("Parquet rows:", nrow(df), " cols:", ncol(df), "\n")
stopifnot("data.table" %in% class(df))
stopifnot(all(c("gene_symbol", "ts_tissue", "rna_rank") %in% names(df)))

# ---- 3. Tissue count per gene ----------------------------------------------
# By ensembl_id: a few symbols map to two IDs and must not be pooled.
n_tissues <- df[, list(n_tissues = uniqueN(ts_tissue)), by = ensembl_id]

# ---- 4. Mechanistic tier (optional) ----------------------------------------
if (file.exists(TIER_CSV)) {
  tiers <- fread(TIER_CSV)
  if (all(c("ensembl_id", "mechanism_tier") %in% names(tiers))) {
    tiers <- unique(tiers[, list(ensembl_id, mechanism_tier)])
    atlas <- merge(atlas, tiers, by = "ensembl_id", all.x = TRUE)
    cat("Mechanism tier merged for ",
        sum(!is.na(atlas$mechanism_tier)), " genes\n", sep = "")
  } else {
    warning("Tier CSV present but lacks expected columns; skipping.")
  }
} else {
  cat("No mechanism tier CSV found; skipping.\n")
}

# ---- 5. Gene x tissue collapsed (concordance-gate machinery) ----------------
tissue_atlas <- df[, list(
  rna_rank_max = max(rna_rank,     na.rm = TRUE),
  rna_mean_tissue = mean(rna_mean,    na.rm = TRUE),
  detect_fraction = mean(y_bin,       na.rm = TRUE),
  protein_max = max(protein_score, na.rm = TRUE),
  n_cell_types = uniqueN(ts_cell_type),
  is_brain = max(is_brain,     na.rm = TRUE)
), by = list(gene_symbol, ensembl_id, ts_tissue)]

tissue_atlas[, tissue_class := data.table::fcase(
  rna_mean_tissue == 0 & detect_fraction == 0, "low_expression",
  rna_rank_max < 0.20 & detect_fraction < 0.10, "low_expression",
  rna_rank_max >= 0.50 & detect_fraction >= 0.50, "concordant_in_tissue",
  rna_rank_max >= 0.50 & detect_fraction <  0.20, "suppressed_in_tissue",
  default = "variable_in_tissue"
)]

cat("tissue_atlas: ", nrow(tissue_atlas), " (gene, tissue) rows\n", sep = "")
print(table(tissue_atlas$tissue_class, useNA = "ifany"))

# ---- 6. Merge gene-level annotations ---------------------------------------
concordance_atlas <- atlas
if ("n_tissues" %in% names(concordance_atlas))
  concordance_atlas[["n_tissues"]] <- NULL
concordance_atlas <- merge(concordance_atlas, n_tissues, by = "ensembl_id", all.x = TRUE)

# ---- 7. Column ordering ----------------------------------------------------
preferred_order <- c(
  "gene_symbol", "ensembl_id",
  "gene_class",  "protein_confidence", "detection_rate", "suppression_rate",
  "mean_rna_rank",
  "n_tissues", "n_concordant_pos", "n_suppressed",
  "mechanism_tier"
)
preferred_order <- intersect(preferred_order, names(concordance_atlas))
remaining  <- setdiff(names(concordance_atlas), preferred_order)
setcolorder(concordance_atlas, c(preferred_order, remaining))
concordance_atlas <- as.data.frame(concordance_atlas)

cat("\n==== concordance_atlas ====\n")
cat("Dimensions: ", nrow(concordance_atlas), " genes x ",
    ncol(concordance_atlas), " columns\n", sep = "")
print(table(concordance_atlas$gene_class, useNA = "ifany"))

# ===========================================================================
# External platform data: PaxDB, GTEx protein TS, GTEx RNA TS, GTEx medians
# ===========================================================================
# PaxDb is built once (atlas tissue names) and reused: directly for the
# tissue_atlas concordance machinery, and remapped to HPA names for the
# hpa_tissue presence vote.

PAXDB_TISSUE_MAP <- list(
  "cerebral-cortex" = "Cerebral cortex",
  "heart" = "Heart",
  "kidney" = "Kidney",
  "liver" = "Liver",
  "lung" = "Lung",
  "lymph_node" = "Lymph Node",
  "pancreas" = "Pancreas",
  "skin" = "Skin",
  "spleen" = "Spleen",
  "stomach" = "Stomach",
  "testis" = "Testis",
  "colon" = "Colon" 
)

PAXDB_DIR <- file.path(EXTERNAL_DIR, "paxdb")

paxdb_raw <- lapply(names(PAXDB_TISSUE_MAP), function(name) {
  df0 <- read.csv(file.path(PAXDB_DIR, paste0(name, ".txt")),
                  skip = 11, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE, comment.char = "")
  colnames(df0) <- c("paxdb_symbol", "paxdb_string_id", "abundance")
  df0$ts_tissue <- PAXDB_TISSUE_MAP[[name]]
  df0
})
paxdb_raw <- do.call(rbind, paxdb_raw)
paxdb_raw$ensp <- sub("^9606\\.", "", paxdb_raw$paxdb_string_id)
cat(sprintf("PaxDb raw: %d rows, %d unique ENSP\n",
            nrow(paxdb_raw), length(unique(paxdb_raw$ensp))))

ENSP_MAP_CACHE <- file.path(EXTERNAL_DIR, "ensp_to_ensg_biomart.txt")
if (!file.exists(ENSP_MAP_CACHE))
  stop("ENSP -> ENSG cache not found: ", ENSP_MAP_CACHE, call. = FALSE)
ensp_map <- read.csv(ENSP_MAP_CACHE, sep = "\t", stringsAsFactors = FALSE)
colnames(ensp_map) <- c("ensembl_id", "ensp", "hgnc_symbol")
ensp_map <- ensp_map[!duplicated(ensp_map$ensp), ]

paxdb <- paxdb_raw %>%
  inner_join(ensp_map[, c("ensp", "ensembl_id")], by = "ensp") %>%
  group_by(ensembl_id, ts_tissue) %>%
  summarise(abundance = max(abundance, na.rm = TRUE),
            paxdb_ensp = paste(unique(ensp), collapse = ";"),
            .groups = "drop")

n_atlas_with_paxdb <- sum(unique(tissue_atlas$ensembl_id) %in% paxdb$ensembl_id)
cat(sprintf("PaxDb: %d (gene, tissue) rows, %d atlas genes covered\n",
            nrow(paxdb), n_atlas_with_paxdb))

# --- GTEx protein TS and RNA TS (direction, for tissue_atlas concordance) ---
# Every GTEx column that maps to one of the 24 atlas tissues. 

GTEX_TO_ATLAS <- c(
  "Brain.Cerebellum" = "Cerebellum",
  "Brain.Cortex" = "Cerebral cortex",
  "Heart.Atrial" = "Heart",
  "Heart.Ventricle" = "Heart",
  "Breast" = "Mammary",
  "Liver" = "Liver",
  "Lung" = "Lung",
  "Minor.Salivary" = "Salivary Gland",
  "Muscle.Skeletal" = "Muscle",
  "Pancreas" = "Pancreas",
  "Prostate" = "Prostate",
  "Skin.Unexpo" = "Skin",
  "Skin.SunExpo" = "Skin",
  "Small.Intestine" = "Small Intestine",
  "Spleen" = "Spleen",
  "Stomach" = "Stomach",
  "Testis" = "Testis",
  "Uterus" = "Uterus"
)

gtex_ts <- read.csv(file.path(EXTERNAL_DIR, "gtex_protein_ts_score.csv")) %>%
  dplyr::select(ensembl_id, any_of(names(GTEX_TO_ATLAS))) %>%
  tidyr::pivot_longer(-ensembl_id, names_to = "gtex_col", values_to = "raw") %>%
  mutate(
    ts_protein_gtex = as.numeric(sapply(strsplit(raw, ";", fixed = TRUE), `[`, 1)),
    flag = sapply(strsplit(raw, ";", fixed = TRUE), `[`, 2),
    ts_tissue = GTEX_TO_ATLAS[gtex_col]
  ) %>%
  filter(is.na(flag) | flag != "NA_one_rep_in_raw") %>%
  group_by(ensembl_id, ts_tissue) %>%
  summarise(ts_protein_gtex = median(ts_protein_gtex, na.rm = TRUE), .groups = "drop")

gtex_rna_ts <- read.csv(file.path(EXTERNAL_DIR, "gtex_rna_ts_score.csv")) %>%
  dplyr::select(ensembl_id, any_of(names(GTEX_TO_ATLAS))) %>%
  tidyr::pivot_longer(-ensembl_id, names_to = "gtex_col", values_to = "raw") %>%
  mutate(
    ts_rna_gtex = as.numeric(sapply(strsplit(raw, "[; ]"), `[`, 1)),
    flag = sapply(strsplit(raw, "[; ]"), `[`, 2),
    ts_tissue = GTEX_TO_ATLAS[gtex_col]
  ) %>%
  filter(is.na(flag) | !flag %in% c("NA_raw_tpm_less_1", "all_tissues_tpm_less_1")) %>%
  group_by(ensembl_id, ts_tissue) %>%
  summarise(ts_rna_gtex = median(ts_rna_gtex, na.rm = TRUE), .groups = "drop")

# --- GTEx medians (presence only, for the hpa_tissue vote) ------------------
# Read once. Non-NA = quantified in the tissue. Not joined onto tissue_atlas.
gtex_medians_raw <- read.csv(file.path(EXTERNAL_DIR, "gtex_protein_tissue_median.csv"),
                             stringsAsFactors = FALSE)

# ===========================================================================
# tissue_atlas
# ===========================================================================
tissue_atlas <- tissue_atlas %>%
  left_join(gtex_ts, by = c("ensembl_id", "ts_tissue")) %>%
  left_join(gtex_rna_ts, by = c("ensembl_id", "ts_tissue")) %>%
  left_join(paxdb, by = c("ensembl_id", "ts_tissue")) %>%
  mutate(
    ms_detected = !is.na(ts_protein_gtex),
    paxdb_detected = !is.na(abundance),
    ihc_detected = protein_max >= 1,
    protein_evidence_sources = ms_detected + paxdb_detected + ihc_detected,
    gtex_concordance = case_when(
      is.na(ts_rna_gtex) | is.na(ts_protein_gtex) ~ NA_character_,
      ts_rna_gtex > 0 & ts_protein_gtex > 0  ~ "gtex_concordant",
      ts_rna_gtex > 0 & ts_protein_gtex <= 0 ~ "gtex_suppressed",
      ts_rna_gtex <= 0 & ts_protein_gtex > 0 ~ "gtex_stabilised",
      TRUE ~ "gtex_low_both"
    ),
    platforms_agree = case_when(
      is.na(gtex_concordance) ~ NA,
      tissue_class == "suppressed_in_tissue" &
        gtex_concordance == "gtex_suppressed"  ~ TRUE,
      tissue_class == "concordant_in_tissue" &
        gtex_concordance == "gtex_concordant" ~ TRUE,
      tissue_class %in% c("suppressed_in_tissue", "concordant_in_tissue") ~ FALSE,
      TRUE ~ NA
    )
  )

# ===========================================================================
# Tissue gate: three-source union presence vote
# ===========================================================================
# Built over the union of HPA, GTEx, and PaxDb: every gene any source sees,
# across every tissue any source covers. Genes HPA never stained still carry
# GTEx and PaxDb evidence.
HPA_IHC_PATH <- file.path(EXTERNAL_DIR, "normal_ihc_data.tsv")
stopifnot(file.exists(HPA_IHC_PATH))

hpa_raw <- fread(HPA_IHC_PATH)
setnames(hpa_raw, c("Gene", "Gene name", "Tissue", "IHC tissue name",
                    "Cell type", "Level", "Reliability"),
         c("ensembl_id", "gene_symbol", "tissue", "ihc_tissue_name",
           "cell_type", "level", "reliability"))
hpa_raw <- hpa_raw[reliability != "Uncertain"]
hpa_raw[, level_num := fcase(
  level == "High", 3L, level == "Medium", 2L, level == "Low", 1L, default = 0L)]
# Reliability tiers after the Uncertain filter, for reporting under details.
hpa_raw[, reliability_num := fcase(
  reliability == "Enhanced", 3L, reliability == "Supported", 2L,
  reliability == "Approved", 1L, default = 0L)]

# HPA IHC summary per (gene, tissue). hpa_call is present for every HPA pair.
hpa_ihc <- hpa_raw[,
  list(
    ihc_max_level = max(level_num, na.rm = TRUE),
    ihc_detected = any(level_num > 0, na.rm = TRUE),
    ihc_n_cell_types = .N,
    ihc_n_detected = sum(level_num > 0, na.rm = TRUE),
    ihc_reliability_max = max(reliability_num, na.rm = TRUE)
  ),
  by = list(ensembl_id, gene_symbol, tissue)
]
hpa_ihc[, hpa_call := ifelse(ihc_detected, "detected", "not_detected")]
hpa_ihc[, ihc_level := fcase(
  ihc_max_level == 3L, "High", ihc_max_level == 2L, "Medium",
  ihc_max_level == 1L, "Low", default = "Not detected")]
hpa_ihc <- as.data.frame(hpa_ihc)
cat("HPA IHC:", length(unique(hpa_ihc$ensembl_id)), "genes,",
    length(unique(hpa_ihc$tissue)), "tissues\n")

# --- GTEx presence in HPA tissue names (in-universe absence votes) ----------
GTEX_TO_HPA <- c(
  "Adrenal.Gland" = "Adrenal gland",  "Brain.Cerebellum" = "Cerebellum",
  "Brain.Cortex" = "Cerebral cortex","Breast" = "Breast",
  "Colon.Sigmoid" = "Colon", "Colon.Transverse" = "Colon",
  "Esophagus.Mucosa" = "Esophagus", "Esophagus.Muscle" = "Esophagus",
  "Heart.Atrial" = "Heart muscle", "Heart.Ventricle" = "Heart muscle",
  "Liver" = "Liver", "Lung" = "Lung",
  "Minor.Salivary" = "Salivary gland", "Muscle.Skeletal" = "Skeletal muscle",
  "Ovary" = "Ovary", "Pancreas" = "Pancreas",
  "Pituitary" = "Pituitary gland","Prostate" = "Prostate",
  "Skin.Unexpo" = "Skin", "Small.Intestine"  = "Small intestine",
  "Spleen" = "Spleen", "Stomach" = "Stomach",
  "Testis" = "Testis", "Thyroid" = "Thyroid gland",
  "Vagina" = "Vagina"
)

gtex_presence <- gtex_medians_raw %>%
  rename(ensembl_id = gene.id) %>%
  dplyr::select(ensembl_id, any_of(names(GTEX_TO_HPA))) %>%
  tidyr::pivot_longer(-ensembl_id, names_to = "gtex_col", values_to = "val") %>%
  mutate(tissue = GTEX_TO_HPA[gtex_col]) %>%
  group_by(ensembl_id, tissue) %>%
  summarise(val = { v <- val[!is.na(val)]; if (length(v)) median(v) else NA_real_ },
            .groups = "drop") %>%
  group_by(ensembl_id) %>%
  mutate(gtex_in_universe = any(!is.na(val))) %>%
  ungroup() %>%
  mutate(gtex_call = case_when(
    !is.na(val) ~ "detected",
    is.na(val) & gtex_in_universe ~ "not_detected",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(gtex_call)) %>%
  dplyr::select(ensembl_id, tissue, gtex_call)

# --- PaxDb presence in HPA tissue names (in-universe absence votes) ---------
PAXDB_TO_HPA <- c(
  "Cerebral cortex" = "Cerebral cortex", "Heart" = "Heart muscle",
  "Kidney" = "Kidney", "Liver" = "Liver",
  "Lung" = "Lung", "Lymph Node" = "Lymph node",
  "Pancreas" = "Pancreas", "Skin" = "Skin",
  "Spleen" = "Spleen", "Stomach" = "Stomach",
  "Testis" = "Testis", "Colon"   = "Colon"
)
paxdb_universe <- unique(paxdb$ensembl_id)
paxdb_hpa_tissues <- unique(unname(PAXDB_TO_HPA))

paxdb_presence <- tidyr::expand_grid(ensembl_id = paxdb_universe,
                                     tissue = paxdb_hpa_tissues) %>%
  left_join(paxdb %>%
              filter(ts_tissue %in% names(PAXDB_TO_HPA)) %>%
              transmute(ensembl_id, tissue = PAXDB_TO_HPA[ts_tissue], seen = TRUE),
            by = c("ensembl_id", "tissue")) %>%
  mutate(paxdb_call = ifelse(is.na(seen), "not_detected", "detected")) %>%
  dplyr::select(ensembl_id, tissue, paxdb_call)

# --- Union the three sources on (ensembl_id, tissue) ------------------------
# A row exists wherever any source has a call. HPA, GTEx, or PaxDb vote
# detected or not_detected; an absent source is NA and casts no vote. Presence
# and absence are asserted only with two or more voters, so a single hit is
# detected_by_one and a single platform cannot assert absence.
is_det <- function(x) !is.na(x) & x == "detected"
is_abs <- function(x) !is.na(x) & x == "not_detected"

hpa_tissue <- hpa_ihc %>%
  dplyr::full_join(gtex_presence,  by = c("ensembl_id", "tissue")) %>%
  dplyr::full_join(paxdb_presence, by = c("ensembl_id", "tissue")) %>%
  mutate(
    n_detected = is_det(hpa_call) + is_det(gtex_call) + is_det(paxdb_call),
    n_absent = is_abs(hpa_call) + is_abs(gtex_call) + is_abs(paxdb_call),
    n_voters = n_detected + n_absent,
    seen_by = sub(";$", "", paste0(
      ifelse(is_det(hpa_call),   "HPA;",   ""),
      ifelse(is_det(gtex_call),  "GTEx;",  ""),
      ifelse(is_det(paxdb_call), "PaxDb;", ""))),
    verdict = case_when(
      n_voters < 2L & n_detected == 1L ~ "detected_by_one",
      n_voters < 2L ~ "insufficient_coverage",
      n_detected > 0 & n_absent == 0 ~ "detected",
      n_absent   > 0 & n_detected == 0 ~ "concordant_absent",
      TRUE ~ "discordant"),
    any_detected = n_detected > 0,
    confident_detected = (!is.na(ihc_max_level) & ihc_max_level >= 2L) |
                          (n_detected >= 2)
  )

# Backfill gene_symbol for genes outside HPA (GTEx/PaxDb only) from the atlas
# and the PaxDb biomaRt map; leave NA where no source names the gene.
sym_sources <- rbind(
  data.frame(ensembl_id = hpa_ihc$ensembl_id, gene_symbol = hpa_ihc$gene_symbol,  stringsAsFactors = FALSE),
  data.frame(ensembl_id = atlas$ensembl_id,   gene_symbol = atlas$gene_symbol,    stringsAsFactors = FALSE),
  data.frame(ensembl_id = ensp_map$ensembl_id, gene_symbol = ensp_map$hgnc_symbol, stringsAsFactors = FALSE)
)
sym_sources <- sym_sources[!is.na(sym_sources$gene_symbol) & sym_sources$gene_symbol != "", ]
sym_sources <- sym_sources[!duplicated(sym_sources$ensembl_id), ]
need <- is.na(hpa_tissue$gene_symbol)
hpa_tissue$gene_symbol[need] <-
  sym_sources$gene_symbol[match(hpa_tissue$ensembl_id[need], sym_sources$ensembl_id)]

cat("\n==== tissue gate (three-source union) ====\n")
cat("genes:          ", length(unique(hpa_tissue$ensembl_id)), "\n")
cat("tissues:        ", length(unique(hpa_tissue$tissue)), "\n")
cat("rows:           ", nrow(hpa_tissue), "\n")
cat("rows w/ symbol: ", sum(!is.na(hpa_tissue$gene_symbol)), "\n")
print(table(hpa_tissue$verdict, useNA = "ifany"))

# ===========================================================================
# Gene-level tissue summaries (protein_tissues from HPA, rna_tissues from atlas)
# ===========================================================================
protein_tissue_summary <- as.data.table(hpa_tissue)[any_detected == TRUE,
  list(protein_tissues = paste(sort(unique(tissue)), collapse = ", "),
       n_tissues_protein = uniqueN(tissue)),
  by = ensembl_id]
protein_tissue_summary <- as.data.frame(protein_tissue_summary)

rna_tissue_summary <- as.data.table(tissue_atlas)[rna_rank_max >= 0.5,
  list(rna_tissues = paste(sort(unique(ts_tissue)), collapse = ", "),
       n_tissues_rna = uniqueN(ts_tissue)),
  by = ensembl_id]
rna_tissue_summary <- as.data.frame(rna_tissue_summary)

concordance_atlas <- concordance_atlas %>%
  left_join(protein_tissue_summary, by = "ensembl_id") %>%
  left_join(rna_tissue_summary, by = "ensembl_id")

# ===========================================================================
# Specificity inputs: tau (RNA magnitude) and protein_breadth (any platform)
# ===========================================================================
# tau over log1p(rna_mean_tissue) across the tissues in
# which the gene was profiled. The gate asks whether
# production concentrates anywhere.
# protein_breadth is the fraction of profiled tissues with protein evidence
# from any of IHC, GTEx MS, or PaxDb. The denominator is the number of tissues
# the gene was profiled in, not the atlas total.
# Grouped by ensembl_id: a few symbols map to two IDs and must not be pooled.
TAU_MIN_TISSUES <- 10L

specificity_inputs <- as.data.table(tissue_atlas)[, {
  x <- log1p(rna_mean_tissue)
  x[!is.finite(x)] <- NA_real_
  n <- sum(!is.na(x))
  mx <- if (n > 0) max(x, na.rm = TRUE) else NA_real_
  tau <- if (n >= TAU_MIN_TISSUES && !is.na(mx) && mx > 0)
    sum(1 - x[!is.na(x)] / mx) / (n - 1) else NA_real_
  list(
    tau = tau,
    protein_breadth = mean(protein_evidence_sources > 0, na.rm = TRUE),
    n_tissues_profiled = .N
  )
}, by = ensembl_id]
specificity_inputs <- as.data.frame(specificity_inputs)

concordance_atlas <- concordance_atlas %>%
  left_join(specificity_inputs, by = "ensembl_id")

cat("\n==== specificity inputs ====\n")
cat("tau:             ", sum(!is.na(concordance_atlas$tau)), " genes\n", sep = "")
cat("protein_breadth: ", sum(!is.na(concordance_atlas$protein_breadth)), " genes\n", sep = "")
print(summary(concordance_atlas$tau))
print(summary(concordance_atlas$protein_breadth))

# ===========================================================================
# Protein specificity: tau over the HPA IHC panel (the gate input)
# ===========================================================================
# The specificity gate is protein-level, so it reads IHC and nothing else.
# RNA tau (above) and protein_breadth are shipped as evidence only.
# Per (gene, tissue) value: ihc_max_level, the highest staining level in
# any cell type of that tissue, over every tissue HPA stained (the
# 59-tissue panel in hpa_tissue, not the 24-tissue atlas). Level alone.
# Weighting by fraction of cell types detected was tried and rejected: HPA
# scores many tissues on one cell type, so any High there scores the
# maximum, thousands of genes tie, and the tiebreak decides. Tau over that
# level per ensembl_id. NA when the gene has no IHC anywhere or fewer than
# PROTEIN_TAU_MIN_TISSUES stained tissues.
# primary_tissue_ihc is the set of tissues tied at the gene's maximum
# level, HPA names, alphabetical, comma-joined. A single winner is
# undefined for four genes in five (median six tissues share the maximum),
# so the set is the only honest representation. NA when no tissue scores
# above zero.
# Known source gaps (HPA annotation, not our filter): ALB, APOB, CP are
# "Not detected" in hepatocytes at Enhanced or Supported reliability. These
# genes read specificity "na", which is the honest protein-level answer.
PROTEIN_TAU_MIN_TISSUES <- 10L

ihc_panel <- as.data.table(hpa_tissue)[!is.na(ihc_max_level),
  list(ensembl_id, tissue, ihc_max_level, ihc_n_detected, ihc_n_cell_types)]

# Bulk (GTEx MS, PaxDb) detection breadth over the same panel. Read by the
# specificity gate only to decide whether IHC is informative for a gene
# with few IHC hits: a one-hit IHC gene that bulk sees in most covered
# tissues is an antibody failure, one that bulk also sees in one tissue is
# a thin but consistent restriction. Not a specificity measure.
bulk_panel <- as.data.table(hpa_tissue)[, {
  gt <- !is.na(gtex_call); px <- !is.na(paxdb_call)
  cov <- gt | px
  det <- (gt & gtex_call == "detected") | (px & paxdb_call == "detected")
  list(n_tissues_bulk_covered = sum(cov), n_tissues_bulk_detected = sum(det))
}, by = ensembl_id]

protein_specificity <- ihc_panel[, {
  s <- as.numeric(ihc_max_level)
  n <- .N
  mx <- max(s)
  ptau <- if (n >= PROTEIN_TAU_MIN_TISSUES && mx > 0)
    sum(1 - s / mx) / (n - 1) else NA_real_
  top <- if (mx > 0) sort(tissue[s == mx]) else character(0)
  list(
    protein_tau = ptau,
    primary_tissue_ihc = if (length(top)) paste(top, collapse = ", ") else NA_character_,
    n_primary_tissues = length(top),
    n_tissues_ihc = n,
    n_tissues_ihc_detected = sum(s > 0)
  )
}, by = ensembl_id]
protein_specificity <- as.data.frame(protein_specificity)
protein_specificity <- merge(protein_specificity, as.data.frame(bulk_panel),
                             by = "ensembl_id", all.x = TRUE)

concordance_atlas <- concordance_atlas %>%
  left_join(protein_specificity, by = "ensembl_id")

cat("\n==== protein specificity (IHC panel) ====\n")
cat("protein_tau:            ", sum(!is.na(concordance_atlas$protein_tau)), " genes\n", sep = "")
cat("stained, never detected: ", sum(!is.na(concordance_atlas$n_tissues_ihc) &
                                     concordance_atlas$n_tissues_ihc_detected == 0), " genes\n", sep = "")
cat("detected in one tissue: ", sum(concordance_atlas$n_tissues_ihc_detected == 1, na.rm = TRUE), " genes\n", sep = "")
one_hit <- !is.na(concordance_atlas$n_tissues_ihc_detected) & concordance_atlas$n_tissues_ihc_detected == 1
bulk_frac <- concordance_atlas$n_tissues_bulk_detected / concordance_atlas$n_tissues_bulk_covered
cat("  of which bulk-detected in >= half of covered tissues: ",
    sum(one_hit & !is.na(bulk_frac) & bulk_frac >= 0.5), "\n", sep = "")
cat("  of which bulk-detected in exactly one tissue:         ",
    sum(one_hit & concordance_atlas$n_tissues_bulk_detected %in% 1L), "\n", sep = "")
print(summary(concordance_atlas$protein_tau))

# ===========================================================================
# localisation_flags: single UniProt source for the biofluid plausibility gate
# ===========================================================================
localisation_flags <- read.csv(LOC_CSV, stringsAsFactors = FALSE)

# Source spells the id column ensemble_id. Probably my own typo from earlier
names(localisation_flags)[names(localisation_flags) == "ensemble_id"] <- "ensembl_id"

loc_required <- c("ensembl_id", "gene_symbol",
                  "is_secreted", "is_membrane", "is_nuclear", "is_cytoplasmic",
                  "is_mitochondrial", "is_er_golgi", "is_cytoskeleton",
                  "is_projection", "is_multilocal")
loc_missing_cols <- setdiff(loc_required, names(localisation_flags))
if (length(loc_missing_cols) > 0)
  stop("localisation_flags missing columns: ",
       paste(loc_missing_cols, collapse = ", "), call. = FALSE)

for (col in grep("^is_", names(localisation_flags), value = TRUE))
  localisation_flags[[col]] <- as.integer(localisation_flags[[col]])

# Coverage check.
# Coverage check on ensembl_id, which is what .join_localisation() matches
# first. Symbol is the fallback at run time, so a symbol-only miss here is
# not a miss in the package.
atlas_no_loc <- setdiff(concordance_atlas$ensembl_id, localisation_flags$ensembl_id)
sym_rescued <- toupper(concordance_atlas$gene_symbol[concordance_atlas$ensembl_id %in% atlas_no_loc]) %in%
               toupper(localisation_flags$gene_symbol)
if (length(atlas_no_loc) > 0)
  warning(length(atlas_no_loc), " atlas genes have no ensembl_id match in ",
          "localisation_flags (", sum(sym_rescued), " recoverable by symbol, ",
          sum(!sym_rescued), " will read as localisation unknown): ",
          paste(head(concordance_atlas$gene_symbol[concordance_atlas$ensembl_id %in% atlas_no_loc][!sym_rescued], 20),
                collapse = ", "),
          if (sum(!sym_rescued) > 20) " ..." else "")

# ===========================================================================
# Save
# ===========================================================================
concordance_atlas <- as.data.frame(concordance_atlas)
tissue_atlas <- as.data.frame(tissue_atlas)
hpa_tissue <- as.data.frame(hpa_tissue)
localisation_flags <- as.data.frame(localisation_flags)

cat("\n==== Final dimensions ====\n")
cat("concordance_atlas: ", nrow(concordance_atlas), "x", ncol(concordance_atlas), "\n")
cat("tissue_atlas:      ", nrow(tissue_atlas), "x", ncol(tissue_atlas), "\n")
cat("hpa_tissue:        ", nrow(hpa_tissue), "x", ncol(hpa_tissue), "\n")
cat("localisation_flags:", nrow(localisation_flags), "x", ncol(localisation_flags), "\n")

# ---- Thresholds --------------------------------------------------------------
# Defined in data-raw/derive_thresholds.R so the characterisation work that
# justifies them has one home. Sourced here so the values ship with the data.
source("data-raw/derive_thresholds.R", local = TRUE)
stopifnot(is.list(atlas_thresholds))

usethis::use_data(concordance_atlas, tissue_atlas, hpa_tissue, localisation_flags,
                  atlas_thresholds,
                  internal = TRUE, overwrite = TRUE)
