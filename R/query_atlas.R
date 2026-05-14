# =============================================================================
# query_atlas for gene-level lookup
# =============================================================================

#' Look up genes in the RNA--protein concordance atlas
#'
#' The foundational accessor. Returns concordance class, protein confidence,
#' detection rate, mechanistic tier, and subcellular localisation flags
#' for one or more genes. All other concordR functions call this internally.
#'
#' @param genes Character vector of HGNC gene symbols.
#' @param atlas Optional. A gene-level atlas data frame (as returned by
#'   \code{\link{load_atlas}("gene")}). If \code{NULL}, the atlas is loaded
#'   automatically from package data.
#' @param missing Character. How to handle genes not found in the atlas.
#'   \code{"keep"} (default) returns a row with \code{gene_class = NA};
#'   \code{"drop"} silently removes them; \code{"warn"} keeps them and
#'   prints a message.
#'
#' @return A data frame with one row per queried gene (in input order),
#'   containing all atlas columns plus a logical \code{found} column.
#'   Genes not in the atlas have \code{found = FALSE} and \code{NA} for
#'   all atlas fields.
#'
#' @examples
#' \dontrun{
#' # Single gene
#' query_atlas("TREM2")
#'
#' # Panel of NDD-relevant genes
#' query_atlas(c("NEFL", "GFAP", "TREM2", "NRGN", "VILIP1"))
#'
#' # Drop unknowns
#' query_atlas(c("NEFL", "FAKEGENE"), missing = "drop")
#' }
#'
#' @export
query_atlas <- function(genes,
                        atlas   = NULL,
                        missing = c("keep", "drop", "warn")) {

  missing <- match.arg(missing)

  if (!is.character(genes) || length(genes) == 0L)
    stop("`genes` must be a non-empty character vector of gene symbols.")

  genes <- unique(trimws(genes))

  if (is.null(atlas)) atlas <- .ensure_atlas()

  # Match genes (case-insensitive)
  idx <- match(toupper(genes), toupper(atlas$gene_symbol))
  found <- !is.na(idx)

  if (missing == "warn" && any(!found))
    message(sum(!found), " gene(s) not found in atlas: ",
            paste(head(genes[!found], 10), collapse = ", "),
            if (sum(!found) > 10) "..." else "")

  # Build result from matched rows
  result <- atlas[idx[found], , drop = FALSE]
  result$found <- TRUE

  if (any(!found) && missing != "drop") {
    missing_df <- data.frame(
      gene_symbol = genes[!found],
      found       = FALSE,
      stringsAsFactors = FALSE
    )
    # Fill atlas columns with NA
    for (col in setdiff(names(atlas), "gene_symbol")) {
      missing_df[[col]] <- NA
    }
    missing_df <- missing_df[, names(result), drop = FALSE]

    # Interleave to preserve input order
    out <- rbind(result, missing_df)
    out <- out[order(match(toupper(out$gene_symbol), toupper(genes))), ,
               drop = FALSE]
    result <- out
  }

  rownames(result) <- NULL
  result
}


#' Quick summary of atlas coverage for a gene list
#'
#' @param genes Character vector of gene symbols.
#' @return A named list with coverage counts.
#' @keywords internal
#' @noRd
.atlas_summary <- function(genes) {
  q <- query_atlas(genes, missing = "drop")
  n <- nrow(q)
  classes <- table(factor(q$gene_class,
                          levels = c("consistently_concordant", "variable",
                                     "consistently_suppressed",
                                     "low_expression")))
  list(
    n_queried     = length(genes),
    n_found       = n,
    n_concordant  = as.integer(classes["consistently_concordant"]),
    n_variable    = as.integer(classes["variable"]),
    n_suppressed  = as.integer(classes["consistently_suppressed"]),
    n_low         = as.integer(classes["low_expression"]),
    pct_concordant = if (n > 0) classes["consistently_concordant"] / n * 100
                     else NA_real_,
    pct_suppressed = if (n > 0) classes["consistently_suppressed"] / n * 100
                     else NA_real_
  )
}


#' Look up genes in the tissue-resolved atlas
#'
#' @param genes Character vector of HGNC gene symbols.
#' @param tissues Character vector of tissue names, or NULL for all.
#' @return Data frame filtered to matching gene-tissue pairs.
#' @export
query_tissue <- function(genes, tissues = NULL) {
  ta <- .ensure_tissue_atlas()
  idx <- toupper(ta$gene_symbol) %in% toupper(genes)
  if (!is.null(tissues)) {
    idx <- idx & ta$ts_tissue %in% tissues
  }
  ta[idx, , drop = FALSE]
}