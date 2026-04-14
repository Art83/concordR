# =============================================================================
# flag_genes — pre-analysis sanity check
# =============================================================================

#' Flag genes by concordance reliability and compartment plausibility
#'
#' The first function to run on any proteomics gene list. For each gene,
#' reports whether RNA-based evidence for it is trustworthy, whether its
#' subcellular localisation is consistent with the measurement compartment,
#' and whether its tissue expression is broad enough to undermine
#' tissue-of-origin claims.
#'
#' @section Traffic-light logic:
#' Each gene receives a \code{flag} combining concordance and localisation:
#' \describe{
#'   \item{\code{"reliable"}}{Concordant gene with plausible localisation.
#'     RNA evidence is trustworthy for this gene in this context.}
#'   \item{\code{"caution"}}{Variable gene, or concordant with borderline
#'     localisation. Interpret with care.}
#'   \item{\code{"unreliable"}}{Suppressed gene and/or implausible
#'     localisation. RNA evidence is uninformative for this gene.
#'     Protein-level measurements are required.}
#'   \item{\code{NA}}{Gene not in atlas.}
#' }
#'
#' @param genes Character vector of HGNC gene symbols.
#' @param sample_type Character. The biological compartment where proteins
#'   were measured. One of \code{"plasma"}, \code{"serum"}, \code{"csf"},
#'   \code{"urine"}, \code{"cell_lysate"}, \code{"tissue"},
#'   \code{"nuclear_fraction"}, \code{"mitochondrial_fraction"},
#'   \code{"membrane_fraction"}. Determines which localisations are
#'   considered plausible.
#' @param max_breadth Integer. Genes detected in more than this many
#'   tissues (out of 24) are flagged as \code{broad_expression = TRUE},
#'   indicating they are poor candidates for tissue-specific claims.
#'   Default \code{10}.
#'
#' @return A data frame with one row per gene, containing:
#' \describe{
#'   \item{\code{gene_symbol}}{Gene symbol.}
#'   \item{\code{gene_class}}{From the atlas.}
#'   \item{\code{protein_confidence}}{From the atlas.}
#'   \item{\code{detection_rate}}{From the atlas.}
#'   \item{\code{localisation}}{Human-readable localisation label.}
#'   \item{\code{loc_plausible}}{Logical. Is the localisation consistent
#'     with \code{sample_type}?}
#'   \item{\code{broad_expression}}{Logical. Detected in > \code{max_breadth}
#'     tissues.}
#'   \item{\code{flag}}{Traffic light: \code{"reliable"}, \code{"caution"},
#'     or \code{"unreliable"}.}
#' }
#'
#' @examples
#' \dontrun{
#' # Plasma proteomics panel
#' flag_genes(
#'   genes       = c("NEFL", "GFAP", "TREM2", "NRGN", "HDAC1", "TP53"),
#'   sample_type = "plasma"
#' )
#'
#' # CSF biomarker candidates
#' flag_genes(
#'   genes       = c("NEFL", "VILIP1", "NRGN", "CST3"),
#'   sample_type = "csf"
#' )
#' }
#'
#' @seealso \code{\link{query_atlas}}, \code{\link{triage}}
#' @export
flag_genes <- function(genes,
                       sample_type,
                       max_breadth = 10L) {

  if (!is.character(genes) || length(genes) == 0L)
    stop("`genes` must be a non-empty character vector.")

  sample_type <- tolower(sample_type)
  if (!sample_type %in% .SAMPLE_TYPES)
    stop("Unknown `sample_type`: '", sample_type, "'. Available: ",
         paste(.SAMPLE_TYPES, collapse = ", "), ".")

  # --- Query atlas ---
  q <- query_atlas(genes, missing = "keep")

  # --- Localisation plausibility (using binary flags) ---
  q$loc_plausible <- .check_plausibility(q, sample_type)

  # --- Human-readable localisation ---
  q$localisation <- .location_labels(q)

  # --- Tissue breadth (protein-based) ---
  q$broad_expression <- !is.na(q$n_tissues) & q$n_tissues > max_breadth

  # --- Traffic light ---
  q$flag <- .assign_flag(q$gene_class, q$loc_plausible)

  # --- Select output columns ---
  out_cols <- c("gene_symbol", "gene_class", "protein_confidence",
                "detection_rate", "localisation", "loc_plausible",
                "broad_expression", "flag")
  out_cols <- intersect(out_cols, names(q))
  q[, out_cols, drop = FALSE]
}


#' Assign traffic-light flag from gene class and localisation
#' @keywords internal
#' @noRd
.assign_flag <- function(gene_class, loc_plausible) {
  mapply(function(gc, lp) {
    if (is.na(gc)) return(NA_character_)

    # Suppressed genes are unreliable regardless of localisation
    if (gc == "consistently_suppressed") return("unreliable")

    # Implausible localisation makes any gene unreliable in that context
    if (!is.na(lp) && !lp) return("unreliable")

    # Concordant + plausible (or unscored) = reliable
    if (gc == "consistently_concordant") {
      if (is.na(lp) || lp) return("reliable")
      return("caution")
    }

    # Variable genes = caution
    if (gc == "variable") return("caution")

    # Low expression
    if (gc == "low_expression") return("caution")

    NA_character_
  }, gene_class, loc_plausible, USE.NAMES = FALSE)
}
