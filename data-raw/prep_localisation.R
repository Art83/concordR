# =============================================================================
# prep_localisation.R
# Builds localisation_uniprot.csv from a UniProt export and an Ensembl to
# UniProt map.
# =============================================================================
#
# Note: One UniProt entry maps to several Ensembl IDs whenever a protein is
# encoded by a gene family (histone clusters, CT47A, GAGE, SPANX, SSX,
# defensins) or sits on an alternative haplotype (MHC region). The join
# below produces one row per (Entry, Ensembl ID) and keeps all of them.
# An earlier version deduplicated on Entry after the join, which kept one
# Ensembl ID per protein and left 758 atlas genes without a match.
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

UNIPROT_TSV <- "D:/Bioinformatics_projects/hpa_to_tabula/external_data/uniprotkb_AND_model_organism_9606_2026_05_26.tsv"
MAP_TSV     <- "D:/Bioinformatics_projects/hpa_to_tabula/external_data/ensembl_uniprot_map.tsv"
LOC_CSV     <- "D:/Bioinformatics_projects/localisation_uniprot.csv"

uniprot <- read.delim(UNIPROT_TSV, stringsAsFactors = FALSE, check.names = TRUE)
map     <- read.delim(MAP_TSV,     stringsAsFactors = FALSE, check.names = TRUE)

stopifnot(all(c("Entry", "Gene.Names", "Subcellular.location..CC.") %in% names(uniprot)))
colnames(map)[2] <- "Entry"
stopifnot(all(c("Gene.stable.ID", "Entry", "HGNC.symbol") %in% names(map)))

flag <- function(x, pattern) as.integer(grepl(pattern, x, ignore.case = TRUE))

has_reviewed <- "Reviewed" %in% names(uniprot)
if (!has_reviewed)
  cat("NOTE: export has no Reviewed column. Symbol pass disabled: it would admit\n",
      "      TrEMBL fragments. Add 'Reviewed' to the UniProt export columns.\n")

uniprot_loc <- uniprot %>%
  dplyr::rename(location = Subcellular.location..CC.) %>%
  mutate(gene_primary = sub("[ ;].*$", "", trimws(Gene.Names)),
         reviewed = if (has_reviewed) tolower(Reviewed) == "reviewed" else NA) %>%
  dplyr::select(Entry, gene_primary, location, reviewed)

# Pass 1: by accession (the map's UniProt column is Swiss-Prot, so this is
# reviewed only). Pass 2: entries with no accession match, by primary gene
# name against the map's HGNC symbol, restricted to reviewed entries when
# the export says which are. BioMart and UniProt do not always agree on
# which accession represents a gene, and the second pass recovers reviewed
# entries the first misses for that reason.
by_acc <- uniprot_loc %>%
  inner_join(map, by = "Entry", relationship = "many-to-many")
by_sym <- uniprot_loc %>%
  filter(has_reviewed, !Entry %in% by_acc$Entry, nzchar(gene_primary),
         reviewed %in% TRUE) %>%
  inner_join(map %>% dplyr::select(-Entry) %>% distinct(),
             by = c("gene_primary" = "HGNC.symbol"),
             relationship = "many-to-many") %>%
  mutate(HGNC.symbol = gene_primary)
cat("matched by accession:", length(unique(by_acc$Entry)),
    " by symbol:", length(unique(by_sym$Entry)), "\n")

uniprot_loc <- bind_rows(by_acc, by_sym) %>%
  mutate(
    is_secreted = flag(location, "Secreted"),
    is_membrane = flag(location, "Cell membrane|Single-pass|Multi-pass|Cell surface|GPI-anchor"),
    is_nuclear = flag(location, "Nucleus|Nucleoplasm|Nucleolus|Nuclei"),
    is_cytoplasmic = flag(location, "Cytoplasm|Cytosol|Lysosome|Endosome|Peroxisome|Cytoplasmic vesicle|Sarcoplasmic reticulum"),
    is_mitochondrial = flag(location, "Mitochondri"),
    is_er_golgi = flag(location, "Endoplasmic reticulum|Golgi"),
    is_cytoskeleton  = flag(location, "Cytoskeleton|Microtubule"),
    is_projection = flag(location, paste(
      "Cell projection|Cilium|Axon|growth cone|Ranvier|dendrite|synaptic",
      "synapse|neurite|dendritic spine|postsynaptic|presynaptic",
      "nerve terminal|neuronal projection", sep = "|")),
    is_multilocal = as.integer(is_secreted + is_membrane + is_nuclear +
                                  is_cytoplasmic + is_mitochondrial +
                                  is_er_golgi + is_cytoskeleton +
                                  is_projection > 1L),
    has_location = as.integer(!is.na(location) & nzchar(location))
  ) %>%
  dplyr::rename(gene_symbol = HGNC.symbol, ensembl_id = Gene.stable.ID) %>%
  filter(!is.na(ensembl_id), nzchar(ensembl_id))

# One row per Ensembl ID. Where an ID maps to several entries (rare), OR
# the flags: a gene whose products are anywhere in a set of locations is
# reported as being in all of them.
flag_cols <- grep("^is_|^has_", names(uniprot_loc), value = TRUE)
localisation <- uniprot_loc %>%
  group_by(ensembl_id) %>%
  summarise(
    gene_symbol = {
      s <- gene_symbol[!is.na(gene_symbol) & nzchar(gene_symbol)]
      if (length(s)) s[1L] else NA_character_
    },
    across(all_of(flag_cols), ~ as.integer(any(.x == 1L))),
    .groups = "drop"
  )

cat("UniProt entries:      ", length(unique(uniprot$Entry)), "\n")
cat("rows after join:      ", nrow(uniprot_loc), "\n")
cat("Ensembl IDs written:  ", nrow(localisation), "\n")
cat("IDs with no location: ", sum(localisation$has_location == 0L), "\n")
cat("multi-ID entries:     ",
    sum(table(uniprot_loc$Entry) > 1L), " entries covering ",
    sum(table(uniprot_loc$Entry)[table(uniprot_loc$Entry) > 1L]), " IDs\n", sep = "")

write.csv(localisation, LOC_CSV, row.names = FALSE)
cat("written:", LOC_CSV, "\n")
