# =============================================================================
# Internal: shared constants and helpers for the triage() gates
# =============================================================================
# Every gate in triage() is a pure function of the atlas data and the sample
# type. Nothing here reads another gate's output. Thresholds come from the
# shipped atlas_thresholds object through .thresholds().
# the idea is to let a user decide how to combine multiple gates without imposing a funnel.

#' Sample types and their route rule. For biofluids the route gate reports
#' how a protein reaches the fluid, from localisation alone: secreted,
#' membrane (shed), or intracellular (damage release). For subcellular
#' fractions it reports whether the protein's localisation matches the
#' fraction. NULL means the sample type has no route rule ("na").
#' @keywords internal
#' @noRd
.ROUTE_RULES <- list(
  plasma = "biofluid",
  serum = "biofluid",
  csf = "biofluid",
  urine = "biofluid",
  cell_lysate = NULL,
  tissue = NULL,
  nuclear_fraction = c("is_nuclear"),
  mitochondrial_fraction = c("is_mitochondrial"),
  membrane_fraction = c("is_membrane", "is_er_golgi")
)

#' All recognised sample types.
#' @keywords internal
#' @noRd
.SAMPLE_TYPES <- names(.ROUTE_RULES)

#' Sample type to implied tissue, in the hpa_tissue namespace. NULL means
#' systemic: no tissue is implied and the tissue gate reports "systemic_na"
#' unless context_tissue is supplied. CSF is deliberately NULL: roughly 80% of
#' its protein derives from plasma, so automatic brain attribution is unsound, 
#' unlike what a lot of researchers in NDD space seem to believe.
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

#' Biofluids proper. In these, cellular proteomics under-detects secreted and
#' projection proteins that have left the cell, so a suppressed call from the
#' atlas is reported as unconfirmed rather than as suppressed.
#' @keywords internal
#' @noRd
.BIOFLUID_TYPES <- c("plasma", "serum", "csf", "urine")

#' Brain regions as named in hpa_tissue. Used for the directed error when a
#' user passes "Brain", and by audit_call() to expand a brain context to all
#' regions when deriving accept_tissues.
#' @keywords internal
#' @noRd
.BRAIN_REGIONS <- c("Caudate", "Cerebellum", "Cerebral cortex", "Hippocampus",
                    "Hypothalamus", "Choroid plexus", "Dorsal raphe",
                    "Substantia nigra")

#' Localisation flag columns shipped in localisation_flags.
#' @keywords internal
#' @noRd
.LOC_COLS <- c("is_secreted", "is_membrane", "is_nuclear",
               "is_cytoplasmic", "is_mitochondrial", "is_er_golgi",
               "is_cytoskeleton", "is_projection", "is_multilocal")


# -----------------------------------------------------------------------------
# Thresholds
# -----------------------------------------------------------------------------

#' Shipped thresholds (data-raw/derive_thresholds.R). Fails if the data
#' object is missing, since a silent default would hide a stale sysdata.rda.
#' @keywords internal
#' @noRd
.thresholds <- function() {
  if (is.null(.atlas_env$thresholds)) {
    ns <- asNamespace("concordR")
    if (!exists("atlas_thresholds", envir = ns, inherits = FALSE))
      stop("atlas_thresholds not found in package data. Rebuild sysdata.rda ",
           "with data-raw/prep_atlas_data.R and restart R.", call. = FALSE)
    .atlas_env$thresholds <- get("atlas_thresholds", envir = ns)
  }
  .atlas_env$thresholds
}


# -----------------------------------------------------------------------------
# Input validation
# -----------------------------------------------------------------------------

#' Validate sample_type and return it lower-cased.
#' @keywords internal
#' @noRd
.validate_sample_type <- function(sample_type) {
  if (missing(sample_type) || is.null(sample_type) ||
      !is.character(sample_type) || length(sample_type) != 1L)
    stop("`sample_type` is required (single character). Available: ",
         paste(.SAMPLE_TYPES, collapse = ", "), ".", call. = FALSE)
  sample_type <- tolower(trimws(sample_type))
  if (!sample_type %in% .SAMPLE_TYPES)
    stop("Unknown `sample_type`: '", sample_type, "'. Available: ",
         paste(.SAMPLE_TYPES, collapse = ", "), ".", call. = FALSE)
  sample_type
}

#' Validate a user-supplied tissue name against a character vocabulary.
#'
#' @param context_tissue Character scalar.
#' @param valid_tissues Character vector of permitted names, typically
#'   \code{unique(hpa_tissue$tissue)}.
#' @return Canonical tissue string.
#' @keywords internal
#' @noRd
.validate_context_tissue <- function(context_tissue, valid_tissues) {

  if (!is.character(context_tissue) || length(context_tissue) != 1L)
    stop("`context_tissue` must be a single character string.", call. = FALSE)
  if (!is.character(valid_tissues))
    stop("`valid_tissues` must be a character vector.", call. = FALSE)

  key <- tolower(trimws(context_tissue))

  if (key == "brain")
    stop("The atlas resolves brain into regions; specify one directly, e.g. ",
         paste(.BRAIN_REGIONS, collapse = ", "), ". Aggregating them is ",
         "left to the user.", call. = FALSE)

  hit <- which(tolower(valid_tissues) == key)
  if (length(hit) >= 1L) return(valid_tissues[hit[1L]])

  stop("`context_tissue` value '", context_tissue, "' not found. ",
       "Available tissues: ", paste(sort(valid_tissues), collapse = ", "), ".",
       call. = FALSE)
}


# -----------------------------------------------------------------------------
# Tissue-name bridge: hpa_tissue names to tissue_atlas names
# -----------------------------------------------------------------------------
# The tissue gate reads hpa_tissue. The concordance gate's cross-platform
# clause reads platforms_agree from tissue_atlas, whose tissue names differ for
# three labels and by case elsewhere. NA where no counterpart exists, in which
# case the clause does not fire.
.RENAME_HPA_TO_ATLAS <- c(
  "heart muscle" = "Heart",
  "skeletal muscle" = "Muscle",
  "breast" = "Mammary"
)

#' @keywords internal
#' @noRd
.to_atlas_tissue <- function(hpa_tissue_name, tissue_atlas) {
  if (is.null(hpa_tissue_name) || is.na(hpa_tissue_name)) return(NA_character_)
  key <- tolower(hpa_tissue_name)
  if (key %in% names(.RENAME_HPA_TO_ATLAS))
    return(unname(.RENAME_HPA_TO_ATLAS[[key]]))
  atlas_names <- unique(tissue_atlas$ts_tissue)
  hit <- atlas_names[tolower(atlas_names) == key]
  if (length(hit) >= 1L) hit[1L] else NA_character_
}


# -----------------------------------------------------------------------------
# Localisation join
# -----------------------------------------------------------------------------

#' Attach localisation flags to a gene data frame.
#'
#' Joins on ensembl_id where the row has one, then falls back to a
#' case-insensitive symbol match for rows without an ensembl_id or whose
#' ensembl_id is absent from localisation_flags. Rows with no match keep NA
#' in every is_* column, which the route gate reads as "unknown".
#'
#' @param df Data frame with gene_symbol and, optionally, ensembl_id.
#' @return \code{df} with the is_* columns appended (existing is_* columns are
#' replaced).
#' @keywords internal
#' @noRd
.join_localisation <- function(df) {
  loc <- .ensure_localisation_flags()
  n <- nrow(df)

  df[.LOC_COLS] <- NULL
  if (is.null(loc) || n == 0L) {
    for (col in .LOC_COLS) df[[col]] <- rep(NA_integer_, n)
    return(df)
  }

  idx <- rep(NA_integer_, n)
  if ("ensembl_id" %in% names(df)) {
    has_id <- !is.na(df$ensembl_id) & nzchar(df$ensembl_id)
    idx[has_id] <- match(df$ensembl_id[has_id], loc$ensembl_id)
  }
  need <- is.na(idx)
  if (any(need)) {
    idx[need] <- match(toupper(df$gene_symbol[need]), toupper(loc$gene_symbol))
  }

  for (col in .LOC_COLS) {
    df[[col]] <- if (col %in% names(loc)) as.integer(loc[[col]][idx])
                 else rep(NA_integer_, n)
  }
  df
}


# -----------------------------------------------------------------------------
# Route gate
# -----------------------------------------------------------------------------

#' Route by which a protein reaches the sample, from localisation alone.
#'
#' Reads no other gate. Presence in a tissue is a separate fact and the
#' combiners handle the conjunction.
#'
#' @param df Data frame with the is_* columns (one row per gene).
#' @param sample_type Validated sample type.
#' @return Character vector. For biofluids: \code{"secreted"},
#'   \code{"membrane"}, \code{"intracellular"}, \code{"unknown"}. For
#'   subcellular fractions: \code{"matched"}, \code{"unmatched"},
#'   \code{"unknown"}. \code{"na"} for sample types with no route rule.
#' @keywords internal
#' @noRd
.gate_route <- function(df, sample_type) {
  n <- nrow(df)
  rule <- .ROUTE_RULES[[sample_type]]
  if (is.null(rule)) return(rep("na", n))
  if (n == 0L) return(character(0))

  loc_cols <- intersect(.LOC_COLS, names(df))
  if (length(loc_cols) == 0L) return(rep("unknown", n))
  loc_mat <- as.matrix(df[, loc_cols, drop = FALSE])
  storage.mode(loc_mat) <- "numeric"
  has_any <- rowSums(loc_mat == 1, na.rm = TRUE) > 0
  flag <- function(col) if (col %in% names(df)) !is.na(df[[col]]) & df[[col]] == 1
                        else rep(FALSE, n)

  out <- rep("unknown", n)
  if (identical(rule, "biofluid")) {
    out[has_any] <- "intracellular"
    out[flag("is_membrane")] <- "membrane"
    out[flag("is_secreted")] <- "secreted"
  } else {
    matched <- Reduce(`|`, lapply(rule, flag), init = rep(FALSE, n))
    out[has_any & !matched] <- "unmatched"
    out[matched] <- "matched"
  }
  out
}

#' Human-readable localisation label, for details = TRUE.
#' @keywords internal
#' @noRd
.location_labels <- function(df) {
  n <- nrow(df)
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
  cols <- intersect(names(mapping), names(df))
  if (length(cols) == 0L) return(rep(NA_character_, n))

  mat <- as.matrix(df[, cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  apply(mat, 1L, function(row) {
    hit <- which(!is.na(row) & row == 1)
    if (length(hit) == 0L) return(NA_character_)
    paste(mapping[cols[hit]], collapse = ", ")
  })
}
