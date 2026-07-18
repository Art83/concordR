# =============================================================================
# prep_atlas_data.R
# Assembles the concordR internal datasets from the existing pipeline outputs.
#
# Run once, then: usethis::use_data(..., internal = TRUE)
#
# Three objects, non-overlapping jobs:
#   concordance_atlas  gene-level, concordance gate only (gene_class etc.)
#   tissue_atlas       gene x tissue, concordance-gate machinery only
#   hpa_tissue         gene x tissue, three-source UNION (HPA + GTEx + PaxDb),
#                      the tissue gate source plus primary_tissue/breadth
#   localisation_flags single UniProt localisation source for the compartment gate
#
# CHANGES vs the previous build:
#   - removed is_* / gene_feat from the gene atlas (localisation is one source)
#   - removed tau and specificity_class (RNA signal, superseded by protein breadth)
#   - removed the unused gtex_gene_summary join on the gene atlas
#   - dropped the unused gtex_medians join on tissue_atlas
#   - widened the GTEx/PaxDb -> HPA tissue maps in the hpa_tissue block
#   - folded the three-source presence vote into hpa_tissue (no separate object)
#   - load + validate localisation_flags, assert atlas coverage
# =============================================================================

library(data.table)
library(arrow)
library(dplyr)
library(tidyr)

# ---- Paths -----------------------------------------------------------------
PROJECT_ROOT  <- "D:/Bioinformatics_projects/hpa_to_tabula"
ATLAS_CSV     <- file.path(PROJECT_ROOT, "results/ml/gene_discordance_profiles.csv")
TIER_CSV      <- file.path(PROJECT_ROOT, "results/tables/three_way_decomposition.csv")
PARQUET_DIR   <- file.path(PROJECT_ROOT, "results/ml")
EXTERNAL_DIR  <- file.path(PROJECT_ROOT, "external_data")
LOC_CSV       <- "D:/Bioinformatics_projects/localisation_uniprot.csv"

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
                         "n_tissues_concordant",
                         "n_tissues_suppressed")
missing_optional <- setdiff(optional_atlas_cols, names(atlas))
if (length(missing_optional) > 0)
  warning("Atlas CSV is missing optional columns: ",
          paste(missing_optional, collapse = ", "), ".")

atlas <- atlas[!duplicated(atlas[["gene_symbol"]])]
cat("Atlas rows after dedup:", nrow(atlas), "\n")

# ---- 2. Load parquet -------------------------------------------------------
df_raw <- as.data.frame(arrow::read_parquet(PARQUET_PATH))
df     <- as.data.table(df_raw)
rm(df_raw)
cat("Parquet rows:", nrow(df), " cols:", ncol(df), "\n")
stopifnot("data.table" %in% class(df))
stopifnot(all(c("gene_symbol", "ts_tissue", "rna_rank") %in% names(df)))

# ---- 3. Tissue count per gene ----------------------------------------------
# (was section 4; localisation gene_feat block removed, is_* no longer shipped)
n_tissues <- df[, list(n_tissues = uniqueN(ts_tissue)), by = gene_symbol]

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
  rna_rank_max         = max(rna_rank,     na.rm = TRUE),
  rna_mean_tissue      = mean(rna_mean,    na.rm = TRUE),
  detect_fraction      = mean(y_bin,       na.rm = TRUE),
  protein_max          = max(protein_score, na.rm = TRUE),
  n_cell_types         = uniqueN(ts_cell_type),
  is_brain             = max(is_brain,     na.rm = TRUE)
), by = list(gene_symbol, ensembl_id, ts_tissue)]

tissue_atlas[, tissue_class := data.table::fcase(
  rna_mean_tissue == 0 & detect_fraction == 0,            "low_expression",
  rna_rank_max <  0.20 & detect_fraction <  0.10,         "low_expression",
  rna_rank_max >= 0.50 & detect_fraction >= 0.50,         "concordant_in_tissue",
  rna_rank_max >= 0.50 & detect_fraction <  0.20,         "suppressed_in_tissue",
  default = "variable_in_tissue"
)]

cat("tissue_atlas: ", nrow(tissue_atlas), " (gene, tissue) rows\n", sep = "")
print(table(tissue_atlas$tissue_class, useNA = "ifany"))

# ---- 6. Merge gene-level annotations ---------------------------------------
# (gene_feat and tau merges removed; only n_tissues remains)
concordance_atlas <- atlas
if ("n_tissues" %in% names(concordance_atlas))
  concordance_atlas[["n_tissues"]] <- NULL
concordance_atlas <- merge(concordance_atlas, n_tissues, by = "gene_symbol", all.x = TRUE)

# ---- 7. Column ordering ----------------------------------------------------
preferred_order <- c(
  "gene_symbol", "ensembl_id",
  "gene_class",  "protein_confidence", "detection_rate", "suppression_rate",
  "mean_rna_rank",
  "n_tissues", "n_tissues_concordant", "n_tissues_suppressed",
  "mechanism_tier"
)
preferred_order <- intersect(preferred_order, names(concordance_atlas))
remaining       <- setdiff(names(concordance_atlas), preferred_order)
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
  "heart"           = "Heart",
  "kidney"          = "Kidney",
  "liver"           = "Liver",
  "lung"            = "Lung",
  "lymph_node"      = "Lymph Node",
  "pancreas"        = "Pancreas",
  "skin"            = "Skin",
  "spleen"          = "Spleen",
  "stomach"         = "Stomach",
  "testis"          = "Testis",
  "colon"           = "Colon"     # requires external_data/paxdb/colon.txt; drop if absent
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
  summarise(abundance  = max(abundance, na.rm = TRUE),
            paxdb_ensp = paste(unique(ensp), collapse = ";"),
            .groups    = "drop")

n_atlas_with_paxdb <- sum(unique(tissue_atlas$ensembl_id) %in% paxdb$ensembl_id)
cat(sprintf("PaxDb: %d (gene, tissue) rows, %d atlas genes covered\n",
            nrow(paxdb), n_atlas_with_paxdb))

# --- GTEx protein TS and RNA TS (direction, for tissue_atlas concordance) ---
# Every GTEx column that maps to one of the 24 atlas tissues, not a hand-picked
# subset. Split sub-regions fold together in the group_by below, recovering
# Cerebellum and SunExpo skin that the old list dropped.
GTEX_TO_ATLAS <- c(
  "Brain.Cerebellum" = "Cerebellum",
  "Brain.Cortex"     = "Cerebral cortex",
  "Heart.Atrial"     = "Heart",
  "Heart.Ventricle"  = "Heart",
  "Breast"           = "Mammary",
  "Liver"            = "Liver",
  "Lung"             = "Lung",
  "Minor.Salivary"   = "Salivary Gland",
  "Muscle.Skeletal"  = "Muscle",
  "Pancreas"         = "Pancreas",
  "Prostate"         = "Prostate",
  "Skin.Unexpo"      = "Skin",
  "Skin.SunExpo"     = "Skin",
  "Small.Intestine"  = "Small Intestine",
  "Spleen"           = "Spleen",
  "Stomach"          = "Stomach",
  "Testis"           = "Testis",
  "Uterus"           = "Uterus"
)

gtex_ts <- read.csv(file.path(EXTERNAL_DIR, "gtex_protein_ts_score.csv")) %>%
  dplyr::select(ensembl_id, any_of(names(GTEX_TO_ATLAS))) %>%
  tidyr::pivot_longer(-ensembl_id, names_to = "gtex_col", values_to = "raw") %>%
  mutate(
    ts_protein_gtex = as.numeric(sapply(strsplit(raw, ";", fixed = TRUE), `[`, 1)),
    flag      = sapply(strsplit(raw, ";", fixed = TRUE), `[`, 2),
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
    flag      = sapply(strsplit(raw, "[; ]"), `[`, 2),
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
# tissue_atlas: augment with direction (concordance-gate machinery only)
# ===========================================================================
tissue_atlas <- tissue_atlas %>%
  left_join(gtex_ts,     by = c("ensembl_id", "ts_tissue")) %>%
  left_join(gtex_rna_ts, by = c("ensembl_id", "ts_tissue")) %>%
  left_join(paxdb,       by = c("ensembl_id", "ts_tissue")) %>%
  mutate(
    ms_detected    = !is.na(ts_protein_gtex),
    paxdb_detected = !is.na(abundance),
    ihc_detected   = protein_max >= 1,
    protein_evidence_sources = ms_detected + paxdb_detected + ihc_detected,
    gtex_concordance = case_when(
      is.na(ts_rna_gtex) | is.na(ts_protein_gtex) ~ NA_character_,
      ts_rna_gtex > 0 & ts_protein_gtex > 0       ~ "gtex_concordant",
      ts_rna_gtex > 0 & ts_protein_gtex <= 0      ~ "gtex_suppressed",
      ts_rna_gtex <= 0 & ts_protein_gtex > 0      ~ "gtex_stabilised",
      TRUE                                         ~ "gtex_low_both"
    ),
    platforms_agree = case_when(
      is.na(gtex_concordance) ~ NA,
      tissue_class == "suppressed_in_tissue" &
        gtex_concordance == "gtex_suppressed"       ~ TRUE,
      tissue_class == "concordant_in_tissue" &
        gtex_concordance == "gtex_concordant"        ~ TRUE,
      tissue_class %in% c("suppressed_in_tissue", "concordant_in_tissue") ~ FALSE,
      TRUE                                           ~ NA
    )
  )

# ===========================================================================
# Tissue gate: three-source UNION presence vote   (object name: hpa_tissue)
# ===========================================================================
# Built over the union of HPA, GTEx, and PaxDb: every gene any source sees,
# across every tissue any source covers. Genes HPA never stained still carry
# GTEx and PaxDb evidence, which is the whole point of decoupling. The object
# keeps the name hpa_tissue so load_atlas() and triage() need no change; the
# ihc_* columns are simply NA for genes outside the HPA panel.
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

# HPA IHC summary per (gene, tissue). hpa_call is present for every HPA pair.
hpa_ihc <- hpa_raw[,
  list(
    ihc_max_level    = max(level_num, na.rm = TRUE),
    ihc_detected     = any(level_num > 0, na.rm = TRUE),
    ihc_n_cell_types = .N,
    ihc_n_detected   = sum(level_num > 0, na.rm = TRUE)
  ),
  by = list(ensembl_id, gene_symbol, tissue)
]
hpa_ihc[, hpa_call  := ifelse(ihc_detected, "detected", "not_detected")]
hpa_ihc[, ihc_level := fcase(
  ihc_max_level == 3L, "High", ihc_max_level == 2L, "Medium",
  ihc_max_level == 1L, "Low", default = "Not detected")]
hpa_ihc <- as.data.frame(hpa_ihc)
cat("HPA IHC:", length(unique(hpa_ihc$ensembl_id)), "genes,",
    length(unique(hpa_ihc$tissue)), "tissues\n")

# --- GTEx presence in HPA tissue names (in-universe absence votes) ----------
GTEX_TO_HPA <- c(
  "Adrenal.Gland"    = "Adrenal gland",  "Brain.Cerebellum" = "Cerebellum",
  "Brain.Cortex"     = "Cerebral cortex","Breast"           = "Breast",
  "Colon.Sigmoid"    = "Colon",          "Colon.Transverse" = "Colon",
  "Esophagus.Mucosa" = "Esophagus",      "Esophagus.Muscle" = "Esophagus",
  "Heart.Atrial"     = "Heart muscle",   "Heart.Ventricle"  = "Heart muscle",
  "Liver"            = "Liver",          "Lung"             = "Lung",
  "Minor.Salivary"   = "Salivary gland", "Muscle.Skeletal"  = "Skeletal muscle",
  "Ovary"            = "Ovary",          "Pancreas"         = "Pancreas",
  "Pituitary"        = "Pituitary gland","Prostate"         = "Prostate",
  "Skin.Unexpo"      = "Skin",           "Small.Intestine"  = "Small intestine",
  "Spleen"           = "Spleen",         "Stomach"          = "Stomach",
  "Testis"           = "Testis",         "Thyroid"          = "Thyroid gland",
  "Vagina"           = "Vagina"
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
    !is.na(val)                   ~ "detected",
    is.na(val) & gtex_in_universe ~ "not_detected",
    TRUE                          ~ NA_character_)) %>%
  filter(!is.na(gtex_call)) %>%
  dplyr::select(ensembl_id, tissue, gtex_call)

# --- PaxDb presence in HPA tissue names (in-universe absence votes) ---------
PAXDB_TO_HPA <- c(
  "Cerebral cortex" = "Cerebral cortex", "Heart"   = "Heart muscle",
  "Kidney"          = "Kidney",          "Liver"   = "Liver",
  "Lung"            = "Lung",            "Lymph Node" = "Lymph node",
  "Pancreas"        = "Pancreas",        "Skin"    = "Skin",
  "Spleen"          = "Spleen",          "Stomach" = "Stomach",
  "Testis"          = "Testis",          "Colon"   = "Colon"
)
paxdb_universe    <- unique(paxdb$ensembl_id)
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
# and absence are asserted only with two or more voters, so a lone hit is
# detected_by_one and a single platform cannot assert absence.
is_det <- function(x) !is.na(x) & x == "detected"
is_abs <- function(x) !is.na(x) & x == "not_detected"

hpa_tissue <- hpa_ihc %>%
  dplyr::full_join(gtex_presence,  by = c("ensembl_id", "tissue")) %>%
  dplyr::full_join(paxdb_presence, by = c("ensembl_id", "tissue")) %>%
  mutate(
    n_detected = is_det(hpa_call) + is_det(gtex_call) + is_det(paxdb_call),
    n_absent   = is_abs(hpa_call) + is_abs(gtex_call) + is_abs(paxdb_call),
    n_voters   = n_detected + n_absent,
    seen_by    = sub(";$", "", paste0(
      ifelse(is_det(hpa_call),   "HPA;",   ""),
      ifelse(is_det(gtex_call),  "GTEx;",  ""),
      ifelse(is_det(paxdb_call), "PaxDb;", ""))),
    verdict = case_when(
      n_voters < 2L & n_detected == 1L ~ "detected_by_one",
      n_voters < 2L                    ~ "insufficient_coverage",
      n_detected > 0 & n_absent == 0   ~ "detected",
      n_absent   > 0 & n_detected == 0 ~ "concordant_absent",
      TRUE                             ~ "discordant"),
    any_detected       = n_detected > 0,
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
cat("  genes:          ", length(unique(hpa_tissue$ensembl_id)), "\n")
cat("  tissues:        ", length(unique(hpa_tissue$tissue)), "\n")
cat("  rows:           ", nrow(hpa_tissue), "\n")
cat("  rows w/ symbol: ", sum(!is.na(hpa_tissue$gene_symbol)), "\n")
print(table(hpa_tissue$verdict, useNA = "ifany"))

# ===========================================================================
# Gene-level tissue summaries (protein_tissues from HPA, rna_tissues from atlas)
# ===========================================================================
protein_tissue_summary <- as.data.table(hpa_tissue)[any_detected == TRUE,
  list(protein_tissues   = paste(sort(unique(tissue)), collapse = ", "),
       n_tissues_protein = uniqueN(tissue)),
  by = gene_symbol]
protein_tissue_summary <- as.data.frame(protein_tissue_summary)

rna_tissue_summary <- as.data.table(tissue_atlas)[rna_rank_max >= 0.5,
  list(rna_tissues   = paste(sort(unique(ts_tissue)), collapse = ", "),
       n_tissues_rna = uniqueN(ts_tissue)),
  by = gene_symbol]
rna_tissue_summary <- as.data.frame(rna_tissue_summary)

concordance_atlas <- concordance_atlas %>%
  left_join(protein_tissue_summary, by = "gene_symbol") %>%
  left_join(rna_tissue_summary,     by = "gene_symbol")

# ===========================================================================
# localisation_flags: single UniProt source for the compartment gate
# ===========================================================================
localisation_flags <- read.csv(LOC_CSV, stringsAsFactors = FALSE)

# Source spells the id column ensemble_id; normalise.
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

# Coverage check. With one localisation source and no fallback, any atlas gene
# absent here reads as localisation unknown at run time. Report the gap.
# Flip to stop() if you want to enforce full UniProt coverage of the atlas.
atlas_no_loc <- setdiff(toupper(concordance_atlas$gene_symbol),
                        toupper(localisation_flags$gene_symbol))
if (length(atlas_no_loc) > 0)
  warning(length(atlas_no_loc),
          " atlas genes are absent from localisation_flags and will read as ",
          "localisation unknown: ",
          paste(head(atlas_no_loc, 20), collapse = ", "),
          if (length(atlas_no_loc) > 20) " ..." else "")

# ===========================================================================
# Save
# ===========================================================================
concordance_atlas  <- as.data.frame(concordance_atlas)
tissue_atlas       <- as.data.frame(tissue_atlas)
hpa_tissue         <- as.data.frame(hpa_tissue)
localisation_flags <- as.data.frame(localisation_flags)

cat("\n==== Final dimensions ====\n")
cat("  concordance_atlas: ", nrow(concordance_atlas), "x", ncol(concordance_atlas), "\n")
cat("  tissue_atlas:      ", nrow(tissue_atlas), "x", ncol(tissue_atlas), "\n")
cat("  hpa_tissue:        ", nrow(hpa_tissue), "x", ncol(hpa_tissue), "\n")
cat("  localisation_flags:", nrow(localisation_flags), "x", ncol(localisation_flags), "\n")

usethis::use_data(concordance_atlas, tissue_atlas, hpa_tissue, localisation_flags,
                  internal = TRUE, overwrite = TRUE)
