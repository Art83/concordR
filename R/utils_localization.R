# =============================================================================
# Internal: subcellular localisation logic (binary-flag version)
# =============================================================================
# Uses the binary is_* columns from the atlas (UniProt-derived by default;
# HPA IF fallback) instead of parsing semicolon-separated location strings.


#' Sample-type to plausible localisation flags
#'
#' Maps each sample type to the is_* columns that are biologically
#' plausible for detecting a protein in that compartment.
#' NULL means all localisations are plausible (no filtering).
#'
#' @keywords internal
#' @noRd
.PLAUSIBLE_FLAGS <- list(
  plasma                 = c("is_secreted", "is_membrane"),
  serum                  = c("is_secreted", "is_membrane"),
  csf                    = c("is_secreted", "is_membrane"),
  urine                  = c("is_secreted", "is_membrane"),
  cell_lysate            = NULL,
  tissue                 = NULL,
  nuclear_fraction       = c("is_nuclear"),
  mitochondrial_fraction = c("is_mitochondrial"),
  membrane_fraction      = c("is_membrane", "is_er_golgi")
)

#' Substrate-to-tissue mapping
#'
#' For biofluids with a clear anatomical source, maps to the tissues
#' that the package checks automatically. NULL means systemic — no
#' tissue is implied and gene-level evidence is used.
#' @keywords internal
#' @noRd
.SAMPLE_TISSUE_MAP <- list(
  plasma                 = NULL,
  serum                  = NULL,
  csf                    = c("Caudate", "Cerebellum", "Cerebral cortex", "Hippocampus"),
  urine                  = "Kidney",
  cell_lysate            = NULL,
  tissue                 = NULL,
  nuclear_fraction       = NULL,
  mitochondrial_fraction = NULL,
  membrane_fraction      = NULL
)


#' Sample types where projection proteins (axonal/dendritic) can be
#' released by damage and detected via regulated-release or biomarker
#' pathways (NEFL-in-plasma pattern). Used to avoid over-flagging
#' neurodegeneration damage markers as "implausible".
#' @keywords internal
#' @noRd
.DAMAGE_RELEASE_TYPES <- c("plasma", "serum", "csf")


#' Brain regions in the source atlas
#'
#' The source data resolves brain into four regions rather than a single
#' "Brain" tissue. The package treats these as distinct tissues and
#' refuses to aggregate them, since the underlying data does not.
#' @keywords internal
#' @noRd
.BRAIN_REGIONS <- c("Caudate", "Cerebellum", "Cerebral cortex", "Hippocampus")


#' Validate and resolve a user-supplied tissue name against the atlas.
#'
#' Returns the canonical \code{ts_tissue} value if a unique
#' case-insensitive match is found. If the user passes "Brain" or
#' "brain", emits a directed error listing the four regions. If no
#' match, lists available tissues.
#'
#' @param context_tissue User-supplied character.
#' @param tissue_atlas Tissue-resolved atlas (with \code{ts_tissue}).
#' @return Canonical tissue string.
#' @keywords internal
#' @noRd
.validate_context_tissue <- function(context_tissue, tissue_atlas) {

  if (!is.character(context_tissue) || length(context_tissue) != 1L)
    stop("`context_tissue` must be a single character string.",
         call. = FALSE)

  available <- unique(tissue_atlas$ts_tissue)

  # Special-case: discourage aggregating brain regions
  if (tolower(context_tissue) == "brain")
    stop("The atlas resolves brain into four regions; please specify one ",
         "directly: ", paste(.BRAIN_REGIONS, collapse = ", "), ". ",
         "Aggregating them is left to the user.",
         call. = FALSE)

  hit <- which(tolower(available) == tolower(context_tissue))
  if (length(hit) == 1L) return(available[hit])

  stop("`context_tissue` value '", context_tissue, "' not found in atlas. ",
       "Available tissues: ", paste(available, collapse = ", "), ".",
       call. = FALSE)
}

#' All recognised sample types
#' @keywords internal
#' @noRd
.SAMPLE_TYPES <- names(.PLAUSIBLE_FLAGS)

#' All localisation flag columns (must include is_projection)
#' @keywords internal
#' @noRd
.LOC_COLS <- c("is_secreted", "is_membrane", "is_nuclear",
               "is_cytoplasmic", "is_mitochondrial", "is_er_golgi",
               "is_cytoskeleton", "is_projection", "is_multilocal")


#' Vectorised plausibility check across a data frame of genes
#'
#' @param atlas_subset Data frame with is_* columns (one row per gene).
#' @param sample_type Character scalar.
#' @return Logical vector, one per row. NA when no localisation data.
#' @keywords internal
#' @noRd
.check_plausibility <- function(atlas_subset, sample_type) {
  if (!sample_type %in% .SAMPLE_TYPES)
    stop("Unknown sample type: ", sample_type)
  
  plausible_cols <- .PLAUSIBLE_FLAGS[[sample_type]]
  if (is.null(plausible_cols)) return(rep(TRUE, nrow(atlas_subset)))
  
  available <- intersect(plausible_cols, names(atlas_subset))
  if (length(available) == 0L) return(rep(NA, nrow(atlas_subset)))
  
  # Distinguish "has localisation data, none in this compartment" (legitimate
  # compartment_implausible) from "no localisation data recorded" (NA- gate
  # cannot evaluate). Without this distinction, a gene with all zero/NA
  # localisation flags fails compartment in every sample_type, treating
  # absence of evidence as evidence of absence.
  all_loc_cols <- grep("^is_", names(atlas_subset), value = TRUE)
  
  loc_data    <- as.matrix(atlas_subset[, all_loc_cols, drop = FALSE])
  plaus_data  <- as.matrix(atlas_subset[, available,    drop = FALSE])
  storage.mode(loc_data)   <- "numeric"
  storage.mode(plaus_data) <- "numeric"
  
  has_any_positive <- rowSums(loc_data   == 1, na.rm = TRUE) > 0
  has_plausible    <- rowSums(plaus_data == 1, na.rm = TRUE) > 0
  
  ifelse(!has_any_positive, NA,
         ifelse(has_plausible,     TRUE,
                FALSE))
}


#' Vectorised location labels across a data frame
#'
#' Handles empty inputs and rows where all localisation columns are NA
#' (e.g. genes not found in the atlas with missing = "keep").
#'
#' @param atlas_subset Data frame with is_* columns.
#' @return Character vector.
#' @keywords internal
#' @noRd
.location_labels <- function(atlas_subset) {
  n <- nrow(atlas_subset)
  if (n == 0L) return(character(0))
  
  mapping <- c(
    is_secreted      = "Secreted",
    is_membrane      = "Membrane",
    is_nuclear       = "Nuclear",
    is_cytoplasmic   = "Cytoplasmic",
    is_mitochondrial = "Mitochondrial",
    is_er_golgi      = "ER/Golgi",
    is_cytoskeleton  = "Cytoskeletal",
    is_projection   = "Cell projection"
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


#' Summarise compartment composition of a gene set
#'
#' @param atlas_subset Data frame with is_* columns (one row per gene).
#' @return Named numeric vector of fractions.
#' @keywords internal
#' @noRd
.compartment_summary <- function(atlas_subset) {
  cols <- c(is_secreted = "secreted", is_membrane = "membrane",
            is_nuclear = "nuclear", is_cytoplasmic = "cytoplasmic",
            is_mitochondrial = "mitochondrial", is_er_golgi = "er_golgi",
            is_cytoskeleton = "cytoskeletal", is_projection = "projection")
  
  n <- nrow(atlas_subset)
  if (n == 0L) return(setNames(rep(NA_real_, length(cols)), unname(cols)))
  
  result <- vapply(names(cols), function(col) {
    if (!col %in% names(atlas_subset)) return(NA_real_)
    sum(atlas_subset[[col]] == 1, na.rm = TRUE) / n
  }, numeric(1))
  names(result) <- unname(cols)
  result
}


#' Does the sample type support "damage-release" detection of intracellular
#' proteins (e.g. NEFL in plasma after axonal injury)?
#' @keywords internal
#' @noRd
.is_damage_release_context <- function(sample_type) {
  tolower(sample_type) %in% .DAMAGE_RELEASE_TYPES
}
