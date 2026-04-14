# =============================================================================
# Internal: subcellular localisation logic (binary-flag version)
# =============================================================================
# Uses the binary is_* columns from the atlas (matching the Streamlit app
# schema) instead of parsing semicolon-separated location strings.


#' Sample-type to plausible localisation flags
#'
#' Maps each sample type to the is_* columns that are biologically
#' plausible for detecting a protein in that compartment.
#' NULL means all localisations are plausible (no filtering).
#'
#' @keywords internal
#' @noRd
.PLAUSIBLE_FLAGS <- list(
  plasma   = c("is_secreted", "is_membrane"),
  serum    = c("is_secreted", "is_membrane"),
  csf      = c("is_secreted", "is_membrane"),
  urine    = c("is_secreted", "is_membrane"),
  cell_lysate          = NULL,
  tissue               = NULL,
  nuclear_fraction     = c("is_nuclear"),
  mitochondrial_fraction = c("is_mitochondrial"),
  membrane_fraction    = c("is_membrane", "is_er_golgi")
)

#' All recognised sample types
#' @keywords internal
#' @noRd
.SAMPLE_TYPES <- names(.PLAUSIBLE_FLAGS)

#' All localisation flag columns
#' @keywords internal
#' @noRd
.LOC_COLS <- c("is_secreted", "is_membrane", "is_nuclear",
               "is_cytoplasmic", "is_mitochondrial", "is_er_golgi",
               "is_cytoskeleton", "is_multilocal")


#' Vectorised plausibility check across a data frame of genes
#'
#' @param atlas_subset Data frame with is_* columns (one row per gene).
#' @param sample_type Character scalar.
#' @return Logical vector, one per row.
#' @keywords internal
#' @noRd
.check_plausibility <- function(atlas_subset, sample_type) {
  sample_type <- tolower(sample_type)
  if (!sample_type %in% .SAMPLE_TYPES)
    stop("Unknown sample_type '", sample_type, "'. Available: ",
         paste(.SAMPLE_TYPES, collapse = ", "), ".")

  plausible_cols <- .PLAUSIBLE_FLAGS[[sample_type]]
  if (is.null(plausible_cols)) return(rep(TRUE, nrow(atlas_subset)))

  available <- intersect(plausible_cols, names(atlas_subset))
  if (length(available) == 0L) return(rep(NA, nrow(atlas_subset)))

  apply(atlas_subset[, available, drop = FALSE], 1, function(row) {
    if (all(is.na(row))) return(NA)
    any(row == 1, na.rm = TRUE)
  })
}


#' Get human-readable localisation label from binary flags
#'
#' @param gene_row Named vector or single-row data frame with is_* columns.
#' @return Character scalar (comma-separated labels).
#' @keywords internal
#' @noRd
.location_label <- function(gene_row) {
  mapping <- c(
    is_secreted      = "Secreted",
    is_membrane      = "Membrane",
    is_nuclear       = "Nuclear",
    is_cytoplasmic   = "Cytoplasmic",
    is_mitochondrial = "Mitochondrial",
    is_er_golgi      = "ER/Golgi",
    is_cytoskeleton  = "Cytoskeletal"
  )
  present <- names(mapping)[vapply(names(mapping), function(col) {
    v <- gene_row[[col]]
    !is.null(v) && length(v) == 1L && !is.na(v) && v == 1
  }, logical(1))]

  if (length(present) == 0L) return(NA_character_)
  paste(mapping[present], collapse = ", ")
}


#' Vectorised location labels across a data frame
#'
#' @param atlas_subset Data frame with is_* columns.
#' @return Character vector.
#' @keywords internal
#' @noRd
.location_labels <- function(atlas_subset) {
  apply(atlas_subset, 1, .location_label)
}


#' Summarise compartment composition of a gene set
#'
#' @param atlas_subset Data frame with is_* columns (one row per gene).
#' @return Named numeric vector of fractions.
#' @keywords internal
#' @noRd
.compartment_summary <- function(atlas_subset) {
  n <- nrow(atlas_subset)
  if (n == 0L) {
    return(setNames(rep(NA_real_, 7),
                    c("secreted", "membrane", "nuclear", "cytoplasmic",
                      "mitochondrial", "er_golgi", "cytoskeletal")))
  }

  cols <- c(is_secreted = "secreted", is_membrane = "membrane",
            is_nuclear = "nuclear", is_cytoplasmic = "cytoplasmic",
            is_mitochondrial = "mitochondrial", is_er_golgi = "er_golgi",
            is_cytoskeleton = "cytoskeletal")

  result <- vapply(names(cols), function(col) {
    if (!col %in% names(atlas_subset)) return(NA_real_)
    sum(atlas_subset[[col]] == 1, na.rm = TRUE) / n
  }, numeric(1))
  names(result) <- unname(cols)
  result
}
