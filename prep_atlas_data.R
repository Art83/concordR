# =============================================================================
# prep_atlas_data.R
# Assembles the concordance_atlas internal dataset for concordR
# from your existing pipeline outputs.
#
# Run once, then: usethis::use_data(concordance_atlas, internal = TRUE)
# =============================================================================

library(data.table)

# --- 1. Gene discordance profiles (core atlas) ---
atlas <- fread("D:/Bioinformatics_projects/hpa_to_tabula/results/ml/gene_discordance_profiles.csv")

# Verify expected columns exist
stopifnot(all(c("gene_symbol", "gene_class", "protein_confidence",
                "detection_rate", "mean_rna_rank") %in% names(atlas)))

# Deduplicate (keep first per gene)
atlas <- atlas[!duplicated(gene_symbol)]

# --- 2. Gene features (subcellular localisation flags) ---
# These come from master matrix — extract per-gene first-row
df <- arrow::read_parquet("D:/Bioinformatics_projects/hpa_to_tabula/results/ml/df_ml_v4_2026-04-08_mincells50.parquet")

loc_cols <- c("is_secreted", "is_membrane", "is_nuclear",
              "is_cytoplasmic", "is_mitochondrial", "is_er_golgi",
              "is_cytoskeleton", "is_multilocal")
# Keep only columns that exist
loc_cols <- intersect(loc_cols, names(df))

gene_feat <- as.data.table(df)[,
  lapply(.SD, first),
  by = gene_symbol,
  .SDcols = loc_cols
]

# --- 3. Tissue count ---
# Number of tissues with observations per gene
n_tissues <- as.data.table(df)[,
  .(n_tissues = uniqueN(ts_tissue)),
  by = gene_symbol
]

# --- 4. Mechanistic tier (if available) ---
# Expected columns: gene_symbol, mechanistic_tier
tier_file <- "D:/Bioinformatics_projects/hpa_to_tabula/results/ml/four_way_decomposition.csv"

tiers <- fread(tier_file)
atlas <- merge(atlas, tiers[, .(ensembl_id, mechanism_tier)],
                 by = "ensembl_id", all.x = TRUE)


# --- 5. Merge ---
concordance_atlas <- merge(atlas, gene_feat, by = "gene_symbol", all.x = TRUE)


# --- 6. Fill NAs for binary flags ---
for (col in loc_cols) {
  if (col %in% names(concordance_atlas)) {
    concordance_atlas[[col]][is.na(concordance_atlas[[col]])] <- 0L
  }
}

# --- 7. Convert to data.frame (not data.table) for R package ---
concordance_atlas <- as.data.frame(concordance_atlas)

cat("Atlas dimensions:", nrow(concordance_atlas), "genes x",
    ncol(concordance_atlas), "columns\n")
cat("Columns:", paste(names(concordance_atlas), collapse = ", "), "\n")
cat("Class distribution:\n")
print(table(concordance_atlas$gene_class))



