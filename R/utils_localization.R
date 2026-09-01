# =============================================================================
# Internal: subcellular localisation logic and tissue-name helpers
# =============================================================================
# Localisation flags come from a single UniProt source (localisation_flags),
# joined onto every queried gene in triage(). The tissue gate reads hpa_tissue
# (HPA tissue names). The concordance machinery reads tissue_atlas (atlas
# tissue names). The two namespaces are bridged by .hpa_to_atlas_tissue().


#' Sample-type to plausible localisation flags
#' @keywords internal
#' @noRd
.PLAUSIBLE_FLAGS <- list(
  plasma = c("is_secreted", "is_membrane"),
  serum  = c("is_secreted", "is_membrane"),
  csf = c("is_secreted", "is_membrane"),
  urine = c("is_secreted", "is_membrane"),
  cell_lysate            = NULL,
  tissue                 = NULL,
  nuclear_fraction       = c("is_nuclear"),
  mitochondrial_fraction = c("is_mitochondrial"),
  membrane_fraction      = c("is_membrane", "is_er_golgi")
)

#' Substrate-to-tissue mapping (HPA tissue namespace, since the tissue gate
#' reads hpa_tissue). NULL means systemic, no tissue is implied.
#' @keywords internal
#' @noRd
.SAMPLE_TISSUE_MAP <- list(
  plasma = NULL,
  serum = NULL,
  csf = NULL,
  urine = "Kidney",
  cell_lysate = NULL,
  tissue = NULL,
  nuclear_fraction = NULL,
  mitochondrial_fraction = NULL,
  membrane_fraction = NULL
)

#' Sample types where projection proteins can be released by damage and
#' detected in biofluids (the NEFL-in-plasma pattern).
#' @keywords internal
#' @noRd
.DAMAGE_RELEASE_TYPES <- c("plasma", "serum", "csf")

#' Brain regions. Identical strings in the atlas and HPA namespaces, so this
#' list serves both the gate and the bridge.
#' @keywords internal
#' @noRd
.BRAIN_REGIONS <- c("Caudate", "Cerebellum", "Cerebral cortex", "Hippocampus")

#' All recognised sample types
#' @keywords internal
#' @noRd
.SAMPLE_TYPES <- names(.PLAUSIBLE_FLAGS)

#' All localisation flag columns
#' @keywords internal
#' @noRd
.LOC_COLS <- c("is_secreted", "is_membrane", "is_nuclear",
               "is_cytoplasmic", "is_mitochondrial", "is_er_golgi",
               "is_cytoskeleton", "is_projection", "is_multilocal")


# -----------------------------------------------------------------------------
# Tissue-name bridge: HPA namespace (the gate) -> atlas namespace (tissue_atlas)
# -----------------------------------------------------------------------------
# Three labels are genuine renames between the two namespaces. The remainder
# differ at most by case, so a case-insensitive match against the atlas labels
# resolves them. An HPA tissue with no atlas counterpart returns NA, which
# simply means the atlas-side checks (IHC rescue, suppressed cross-platform)
# do not fire for it. Brain regions match outright, so CSF and plasma brain
# work needs no bridge.
.RENAME_HPA_TO_ATLAS <- c(
  "heart muscle"    = "Heart",
  "skeletal muscle" = "Muscle",
  "breast" = "Mammary"
)

#' Translate HPA tissue labels to atlas (tissue_atlas) labels.
#'
#' @param hpa_tissues Character vector of HPA-namespace tissue names.
#' @param tissue_atlas tissue_atlas data frame (for its ts_tissue vocabulary).
#' @return Character vector of atlas labels, NA where no counterpart exists.
#' @keywords internal
#' @noRd
.hpa_to_atlas_tissue <- function(hpa_tissues, tissue_atlas) {
  atlas_tis <- unique(tissue_atlas$ts_tissue)
  vapply(hpa_tissues, function(t) {
    if (is.na(t)) return(NA_character_)
    key <- tolower(t)
    if (key %in% names(.RENAME_HPA_TO_ATLAS))
      return(unname(.RENAME_HPA_TO_ATLAS[[key]]))
    hit <- atlas_tis[tolower(atlas_tis) == key]
    if (length(hit) == 1L) return(hit)
    NA_character_
  }, character(1), USE.NAMES = FALSE)
}


#' Validate and resolve a user-supplied tissue name.
#'
#' Validated against the gate vocabulary (hpa_tissue tissue names). "Brain"
#' is rejected with a directed message listing the four regions.
#'
#' @param context_tissue User-supplied character scalar.
#' @param valid_tissues Character vector of permitted tissue names
#'   (\code{unique(hpa_tissue$tissue)}).
#' @return Canonical tissue string.
#' @keywords internal
#' @noRd
.validate_context_tissue <- function(context_tissue, valid_tissues) {

  if (!is.character(context_tissue) || length(context_tissue) != 1L)
    stop("`context_tissue` must be a single character string.", call. = FALSE)

  if (tolower(context_tissue) == "brain")
    stop("The atlas resolves brain into four regions; please specify one ",
         "directly: ", paste(.BRAIN_REGIONS, collapse = ", "), ". ",
         "Aggregating them is left to the user.", call. = FALSE)

  hit <- which(tolower(valid_tissues) == tolower(context_tissue))
  if (length(hit) == 1L) return(valid_tissues[hit])

  stop("`context_tissue` value '", context_tissue, "' not found. ",
       "Available tissues: ", paste(sort(valid_tissues), collapse = ", "), ".",
       call. = FALSE)
}


#' Vectorised biofluid plausibility across a data frame of genes.
#'
#' @param atlas_subset Data frame with is_* columns (one row per gene).
#' @param sample_type Character scalar.
#' @return Logical vector. NA where the gene has no localisation data, so the
#'   gate cannot evaluate and does not reject.
#' @keywords internal
#' @noRd
.check_plausibility <- function(atlas_subset, sample_type) {
  if (!sample_type %in% .SAMPLE_TYPES)
    stop("Unknown sample type: ", sample_type)

  plausible_cols <- .PLAUSIBLE_FLAGS[[sample_type]]
  if (is.null(plausible_cols)) return(rep(TRUE, nrow(atlas_subset)))

  available <- intersect(plausible_cols, names(atlas_subset))
  if (length(available) == 0L) return(rep(NA, nrow(atlas_subset)))

  # Distinguish "has localisation data, none plausible for the biofluid" (a
  # real implausible call) from "no localisation data" (NA, not evaluable).
  all_loc_cols <- grep("^is_", names(atlas_subset), value = TRUE)

  loc_data <- as.matrix(atlas_subset[, all_loc_cols, drop = FALSE])
  plaus_data <- as.matrix(atlas_subset[, available,    drop = FALSE])
  storage.mode(loc_data) <- "numeric"
  storage.mode(plaus_data) <- "numeric"

  has_any_positive <- rowSums(loc_data   == 1, na.rm = TRUE) > 0
  has_plausible <- rowSums(plaus_data == 1, na.rm = TRUE) > 0

  ifelse(!has_any_positive, NA, ifelse(has_plausible, TRUE, FALSE))
}


#' Vectorised location labels across a data frame.
#' @keywords internal
#' @noRd
.location_labels <- function(atlas_subset) {
  n <- nrow(atlas_subset)
  if (n == 0L) return(character(0))

  mapping <- c(
    is_secreted = "Secreted",
    is_membrane = "Membrane",
    is_nuclear = "Nuclear",
    is_cytoplasmic = "Cytoplasmic",
    is_mitochondrial = "Mitochondrial",
    is_er_golgi = "ER/Golgi",
    is_cytoskeleton = "Cytoskeletal",
    is_projection = "Cell projection"
  )
  cols_to_use <- intersect(names(mapping), colnames(atlas_subset))
  if (length(cols_to_use) == 0L) return(rep(NA_character_, n))

  mat <- as.matrix(atlas_subset[, cols_to_use, drop = FALSE])
  apply(mat, 1, function(row) {
    if (all(is.na(row))) return(NA_character_)
    found <- mapping[cols_to_use[which(row == 1)]]
    if (length(found) == 0L) return(NA_character_)
    paste(found, collapse = ", ")
  })
}


#' Does the sample type support damage-release detection of intracellular
#' proteins (e.g. NEFL in plasma after axonal injury)?
#' @keywords internal
#' @noRd
.is_damage_release_context <- function(sample_type) {
  tolower(sample_type) %in% .DAMAGE_RELEASE_TYPES
}
