# =============================================================================
# query_atlas for gene-level lookup
# =============================================================================

#' Look up genes in the RNA to protein concordance atlas
#'
#' Returns the gene-level atlas row for one or more genes: concordance
#' class, protein confidence, suppression rate, mechanism tier, and the
#' specificity inputs. Localisation flags are not in this table; see
#' \code{\link{triage}}.
#'
#' @param genes Character vector of HGNC gene symbols. Matching is
#'   case-insensitive and repeats are collapsed, so \code{c("NEFL", "nefl")}
#'   returns one row.
#' @param atlas Optional. A gene-level atlas data frame (as returned by
#'   \code{\link{load_atlas}("gene")}). If \code{NULL}, the atlas is loaded
#'   from package data.
#' @param missing Character. How to handle genes not found in the atlas.
#'   \code{"keep"} (default) returns a row with \code{found = FALSE} and
#'   \code{NA} in every atlas column; \code{"drop"} removes them;
#'   \code{"warn"} keeps them and prints a message.
#'
#' @return A data frame with one row per distinct queried gene, in input
#'   order, containing all atlas columns plus a logical \code{found}. For
#'   found genes \code{gene_symbol} carries the atlas spelling; for missing
#'   genes it carries the input as typed.
#'
#' @examples
#' \dontrun{
#' query_atlas("TREM2")
#' query_atlas(c("NEFL", "GFAP", "TREM2", "NRGN", "VILIP1"))
#' query_atlas(c("NEFL", "FAKEGENE"), missing = "drop")
#' }
#'
#' @export
query_atlas <- function(genes,
                        atlas   = NULL,
                        missing = c("keep", "drop", "warn")) {

  missing <- match.arg(missing)

  if (!is.character(genes) || length(genes) == 0L)
    stop("`genes` must be a non-empty character vector of gene symbols.",
         call. = FALSE)

  genes <- trimws(genes)
  genes <- genes[!duplicated(toupper(genes))]

  if (is.null(atlas)) atlas <- .ensure_atlas()

  idx <- match(toupper(genes), toupper(atlas$gene_symbol))
  found <- !is.na(idx)

  if (missing == "warn" && any(!found))
    message(sum(!found), " gene(s) not found in atlas: ",
            paste(head(genes[!found], 10), collapse = ", "),
            if (sum(!found) > 10) "..." else "")

  if (missing == "drop") {
    genes <- genes[found]
    idx <- idx[found]
    found <- found[found]
  }

  # Indexing with NA yields an all-NA row of the right types, in input order.
  result <- atlas[idx, , drop = FALSE]
  result$gene_symbol[!found] <- genes[!found]
  result$found <- found
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
    n_queried = length(genes),
    n_found = n,
    n_concordant = as.integer(classes["consistently_concordant"]),
    n_variable = as.integer(classes["variable"]),
    n_suppressed = as.integer(classes["consistently_suppressed"]),
    n_low = as.integer(classes["low_expression"]),
    pct_concordant = if (n > 0) classes["consistently_concordant"] / n * 100
                     else NA_real_,
    pct_suppressed = if (n > 0) classes["consistently_suppressed"] / n * 100
                     else NA_real_
  )
}


#' Look up genes in the tissue-resolved atlas
#'
#' Returns every (gene, tissue) row of the tissue-resolved atlas for the
#' requested genes, optionally restricted to named tissues. Gene and tissue
#' matching is case-insensitive.
#'
#' @param genes Character vector of HGNC gene symbols.
#' @param tissues Character vector of atlas tissue names, or \code{NULL}
#'   (default) for all 24. Unknown names are an error that lists the
#'   vocabulary.
#' @param missing Character. \code{"warn"} (default) prints a message naming
#'   genes with no rows; \code{"silent"} does not.
#'
#' @return Data frame of matching rows, ordered by input gene then tissue.
#'   Zero rows if nothing matches.
#'
#' @examples
#' \dontrun{
#' query_tissue("GFAP")
#' query_tissue(c("NEFL", "GFAP"), tissues = c("Cerebral cortex", "Liver"))
#' }
#'
#' @export
query_tissue <- function(genes, tissues = NULL, missing = c("warn", "silent")) {
  missing <- match.arg(missing)
  if (!is.character(genes) || length(genes) == 0L)
    stop("`genes` must be a non-empty character vector of gene symbols.",
         call. = FALSE)
  ta <- .ensure_tissue_atlas()

  genes <- trimws(genes)
  genes <- genes[!duplicated(toupper(genes))]
  gene_key <- toupper(genes)

  keep <- toupper(ta$gene_symbol) %in% gene_key

  if (!is.null(tissues)) {
    if (!is.character(tissues) || length(tissues) == 0L)
      stop("`tissues` must be a character vector or NULL.", call. = FALSE)
    vocab <- unique(ta$ts_tissue)
    hit <- match(tolower(trimws(tissues)), tolower(vocab))
    if (anyNA(hit))
      stop("Unknown tissue(s): ",
           paste(tissues[is.na(hit)], collapse = ", "),
           ". Available: ", paste(sort(vocab), collapse = ", "), ".",
           call. = FALSE)
    keep <- keep & ta$ts_tissue %in% vocab[hit]
  }

  out <- ta[keep, , drop = FALSE]

  if (missing == "warn") {
    absent <- genes[!gene_key %in% toupper(out$gene_symbol)]
    if (length(absent) > 0L)
      message(length(absent), " gene(s) with no rows: ",
              paste(head(absent, 10), collapse = ", "),
              if (length(absent) > 10) "..." else "")
  }

  out <- out[order(match(toupper(out$gene_symbol), gene_key), out$ts_tissue), ,
             drop = FALSE]
  rownames(out) <- NULL
  out
}


#' Look up protein presence votes across the HPA tissue panel
#'
#' Returns the three-source presence vote (HPA IHC, GTEx MS, PaxDb) for
#' the requested genes across the 59 HPA tissues. This is the object the
#' tissue gate of \code{\link{triage}} reads, exposed so a user can see
#' every tissue at once rather than one context at a time. Gene and tissue
#' matching is case-insensitive. Tissue names are HPA names.
#'
#' @param genes Character vector of HGNC gene symbols.
#' @param tissues Character vector of HPA tissue names, or \code{NULL}
#'   (default) for all 59. Unknown names are an error that lists the
#'   vocabulary.
#' @param missing Character. \code{"warn"} (default) prints a message naming
#'   genes with no rows; \code{"silent"} does not.
#'
#' @return Data frame with one row per (gene, tissue), ordered by input
#'   gene then tissue. Columns: \code{gene_symbol}, \code{ensembl_id},
#'   \code{tissue}, \code{verdict} (the vote), \code{n_detected},
#'   \code{n_absent}, \code{n_voters}, \code{seen_by}, and the per-source
#'   calls \code{ihc_level}, \code{gtex_call}, \code{paxdb_call}.
#'
#' @examples
#' \dontrun{
#' query_presence("NEFL")
#' query_presence(c("NEFL", "GFAP"),
#'                tissues = c("Cerebral cortex", "Hippocampus", "Liver"))
#' }
#'
#' @seealso \code{\link{triage}}, \code{\link{query_tissue}}
#' @export
query_presence <- function(genes, tissues = NULL, missing = c("warn", "silent")) {
  missing <- match.arg(missing)
  if (!is.character(genes) || length(genes) == 0L)
    stop("`genes` must be a non-empty character vector of gene symbols.",
         call. = FALSE)
  hpa <- .ensure_hpa_tissue()
  if (is.null(hpa))
    stop("hpa_tissue not found in package data.", call. = FALSE)

  genes <- trimws(genes)
  genes <- genes[!duplicated(toupper(genes))]
  gene_key <- toupper(genes)

  keep <- !is.na(hpa$gene_symbol) & toupper(hpa$gene_symbol) %in% gene_key

  if (!is.null(tissues)) {
    tissues <- .validate_context_tissue(tissues, unique(hpa$tissue))
    keep <- keep & hpa$tissue %in% tissues
  }

  cols <- intersect(c("gene_symbol", "ensembl_id", "tissue", "verdict",
                      "n_detected", "n_absent", "n_voters", "seen_by",
                      "ihc_level", "gtex_call", "paxdb_call"), names(hpa))
  out <- hpa[keep, cols, drop = FALSE]

  if (missing == "warn") {
    absent <- genes[!gene_key %in% toupper(out$gene_symbol)]
    if (length(absent) > 0L)
      message(length(absent), " gene(s) with no rows: ",
              paste(head(absent, 10), collapse = ", "),
              if (length(absent) > 10) "..." else "")
  }

  out <- out[order(match(toupper(out$gene_symbol), gene_key), out$tissue), ,
             drop = FALSE]
  rownames(out) <- NULL
  out
}
