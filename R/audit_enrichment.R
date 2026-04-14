# =============================================================================
# audit_enrichment — concordance-corrected enrichment QC
# =============================================================================

#' Audit enrichment results against the concordance atlas
#'
#' Takes any standard enrichment result — from clusterProfiler, enrichR,
#' fgsea, g:Profiler, or a simple data frame — and recomputes significance
#' after removing or downweighting genes with poor RNA--protein concordance.
#' Where the standard and corrected results agree, the finding is robust.
#' Where they diverge, the standard analysis was driven by discordant genes.
#'
#' @section What this function does:
#' For each enriched pathway:
#' \enumerate{
#'   \item Identifies which genes in the pathway overlap the query list.
#'   \item Looks up each overlap gene's concordance class and localisation.
#'   \item Computes: fraction concordant, fraction suppressed, fraction
#'     localisation-plausible.
#'   \item Reruns the hypergeometric test after removing suppressed-class
#'     genes from both the overlap and the background (if
#'     \code{method = "remove"}) or downweighting them (if
#'     \code{method = "weight"}).
#'   \item Reports both original and corrected p-values.
#' }
#'
#' @section Input format:
#' Requires at minimum a data frame with columns mappable to a pathway
#' identifier, a gene list per pathway, and a p-value. Column names are
#' auto-detected from common enrichment tool outputs.
#'
#' @param enrichment_result Data frame of enrichment results.
#' @param gene_list Character vector. The original query gene list.
#' @param gene_sets Named list of character vectors (pathway -> genes).
#'   Required if \code{enrichment_result} lacks a gene overlap column.
#' @param universe Character vector. Background gene universe. Required
#'   for recomputing the hypergeometric test. If \code{NULL}, only
#'   annotation is performed.
#' @param sample_type Character or \code{NULL}. If supplied, adds
#'   localisation plausibility scoring per pathway.
#' @param method Character. \code{"remove"} (default) excludes suppressed
#'   genes; \code{"weight"} downweights by protein confidence.
#' @param alpha Numeric. Significance threshold. Default \code{0.05}.
#'
#' @return The original data frame with additional columns:
#'   \code{n_overlap}, \code{n_concordant}, \code{n_suppressed},
#'   \code{pct_concordant}, \code{pct_suppressed},
#'   \code{pct_loc_plausible} (if \code{sample_type} supplied),
#'   \code{p_corrected}, \code{padj_corrected}, \code{status_change},
#'   \code{concordance_verdict}.
#'
#' @examples
#' \dontrun{
#' audited <- audit_enrichment(
#'   enrichment_result = as.data.frame(ego),
#'   gene_list         = my_genes,
#'   universe          = all_genes,
#'   sample_type       = "plasma"
#' )
#' # Which pathways collapse?
#' audited[audited$status_change == "lost", c("Description", "pct_suppressed")]
#' }
#'
#' @seealso \code{\link{plot_audit}}, \code{\link{flag_genes}}
#' @export
audit_enrichment <- function(enrichment_result,
                             gene_list,
                             gene_sets   = NULL,
                             universe    = NULL,
                             sample_type = NULL,
                             method      = c("remove", "weight"),
                             alpha       = 0.05) {

  method <- match.arg(method)

  if (!is.data.frame(enrichment_result) || nrow(enrichment_result) == 0L)
    stop("`enrichment_result` must be a non-empty data frame.")
  if (!is.character(gene_list) || length(gene_list) == 0L)
    stop("`gene_list` must be a non-empty character vector.")

  atlas <- .ensure_atlas()

  # --- Detect columns ---
  pathway_col <- .detect_col(enrichment_result,
    c("ID", "term", "pathway", "Description", "set", "Term_name",
      "term_name", "source_term", "GO.ID", "category"),
    "pathway/term")

  pval_col <- .detect_col(enrichment_result,
    c("pvalue", "p.adjust", "padj", "p_value", "FDR", "qvalue",
      "P.value", "nominal_p"),
    "p-value", required = FALSE)

  gene_col <- .detect_col(enrichment_result,
    c("geneID", "core_enrichment", "overlap_genes", "Genes", "genes",
      "leadingEdge", "gene_ids"),
    "overlap genes", required = FALSE)

  # --- Extract per-pathway gene lists ---
  pathways <- as.character(enrichment_result[[pathway_col]])

  if (!is.null(gene_sets)) {
    pw_genes <- lapply(pathways, function(pw) {
      if (pw %in% names(gene_sets)) gene_sets[[pw]] else character(0)
    })
  } else if (!is.null(gene_col)) {
    pw_genes <- lapply(enrichment_result[[gene_col]], .parse_gene_string)
  } else {
    stop("Cannot determine pathway genes. Supply `gene_sets` or ensure ",
         "the enrichment result has a gene column.")
  }

  # --- Annotate each pathway ---
  annotations <- lapply(seq_along(pathways), function(i) {
    overlap <- intersect(gene_list, pw_genes[[i]])
    .annotate_pathway(overlap, atlas, sample_type)
  })
  ann_df <- do.call(rbind, annotations)

  # --- Corrected p-values ---
  if (!is.null(universe) && !is.null(gene_sets)) {
    p_corr <- vapply(seq_along(pathways), function(i) {
      .corrected_pvalue(gene_list, pw_genes[[i]], universe, atlas, method)
    }, numeric(1))

    ann_df$p_corrected    <- p_corr
    ann_df$padj_corrected <- p.adjust(p_corr, method = "BH")

    if (!is.null(pval_col)) {
      orig_sig <- enrichment_result[[pval_col]] <= alpha
      corr_sig <- ann_df$padj_corrected <= alpha
      ann_df$status_change <- ifelse(
        orig_sig & corr_sig, "robust",
        ifelse(orig_sig & !corr_sig, "lost",
               ifelse(!orig_sig & corr_sig, "gained", "ns")))
    }
  } else if (!is.null(universe) && is.null(gene_sets)) {
    message("Corrected p-values require `gene_sets` (full pathway memberships). ",
            "The `geneID` column only contains overlap genes. ",
            "Skipping p-value correction; annotation columns still added.")
  }

  # --- Concordance verdict ---
  ann_df$concordance_verdict <- .verdict(
    ann_df$pct_concordant,
    ann_df$pct_suppressed,
    if ("status_change" %in% names(ann_df)) ann_df$status_change else NULL
  )

  cbind(enrichment_result, ann_df)
}


# =============================================================================
# Internal helpers
# =============================================================================

#' @keywords internal
#' @noRd
.detect_col <- function(df, candidates, desc, required = TRUE) {
  for (col in candidates) {
    if (col %in% names(df)) return(col)
  }
  if (required)
    stop("Could not detect ", desc, " column. Expected one of: ",
         paste(candidates, collapse = ", "), ".")
  NULL
}

#' @keywords internal
#' @noRd
.parse_gene_string <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character(0))
  if (is.list(x)) return(unlist(x))
  trimws(strsplit(as.character(x), "[/,;]")[[1]])
}

#' Annotate a single pathway's overlap genes against atlas
#' @keywords internal
#' @noRd
.annotate_pathway <- function(overlap_genes, atlas, sample_type = NULL) {
  n <- length(overlap_genes)
  if (n == 0L) {
    out <- data.frame(
      n_overlap = 0L, n_concordant = 0L, n_suppressed = 0L,
      n_variable = 0L, pct_concordant = NA_real_,
      pct_suppressed = NA_real_, stringsAsFactors = FALSE)
    if (!is.null(sample_type)) out$pct_loc_plausible <- NA_real_
    return(out)
  }

  idx <- match(toupper(overlap_genes), toupper(atlas$gene_symbol))
  matched <- atlas[idx[!is.na(idx)], , drop = FALSE]
  n_found <- nrow(matched)
  classes <- matched$gene_class

  n_conc <- sum(classes == "consistently_concordant", na.rm = TRUE)
  n_supp <- sum(classes == "consistently_suppressed", na.rm = TRUE)
  n_var  <- sum(classes == "variable", na.rm = TRUE)

  out <- data.frame(
    n_overlap      = n,
    n_concordant   = n_conc,
    n_suppressed   = n_supp,
    n_variable     = n_var,
    pct_concordant = if (n_found > 0) n_conc / n_found * 100 else NA_real_,
    pct_suppressed = if (n_found > 0) n_supp / n_found * 100 else NA_real_,
    stringsAsFactors = FALSE)

  if (!is.null(sample_type)) {
    plausible <- .check_plausibility(matched, sample_type)
    out$pct_loc_plausible <- if (length(plausible) > 0)
      sum(plausible, na.rm = TRUE) / sum(!is.na(plausible)) * 100
    else NA_real_
  }

  out
}

#' Recompute hypergeometric p-value after concordance correction
#' @keywords internal
#' @noRd
.corrected_pvalue <- function(gene_list, pathway_genes, universe,
                              atlas, method) {

  if (method == "remove") {
    suppressed <- atlas$gene_symbol[
      !is.na(atlas$gene_class) &
        atlas$gene_class == "consistently_suppressed"]
    suppressed <- toupper(suppressed)

    gene_list_c     <- setdiff(toupper(gene_list), suppressed)
    pathway_genes_c <- setdiff(toupper(pathway_genes), suppressed)
    universe_c      <- setdiff(toupper(universe), suppressed)

    overlap <- length(intersect(gene_list_c, pathway_genes_c))
    n_list  <- length(intersect(gene_list_c, universe_c))
    n_pw    <- length(intersect(pathway_genes_c, universe_c))
    n_univ  <- length(universe_c)

    if (n_list == 0 || n_pw == 0 || n_univ == 0) return(1.0)
    phyper(overlap - 1, n_pw, n_univ - n_pw, n_list, lower.tail = FALSE)

  } else {
    # Weighted: use protein_confidence as weight
    overlap_genes <- intersect(toupper(gene_list), toupper(pathway_genes))
    overlap_genes <- intersect(overlap_genes, toupper(universe))
    if (length(overlap_genes) == 0L) return(1.0)

    idx <- match(overlap_genes, toupper(atlas$gene_symbol))
    weights <- ifelse(is.na(idx), 0.5, atlas$protein_confidence[idx])
    weights[is.na(weights)] <- 0.5

    effective_overlap <- sum(weights)
    n_list <- length(intersect(toupper(gene_list), toupper(universe)))
    n_pw   <- length(intersect(toupper(pathway_genes), toupper(universe)))
    n_univ <- length(universe)

    phyper(round(effective_overlap) - 1, n_pw, n_univ - n_pw, n_list,
           lower.tail = FALSE)
  }
}

#' @keywords internal
#' @noRd
.verdict <- function(pct_concordant, pct_suppressed, status_change = NULL) {
  mapply(function(pc, ps, sc) {
    if (is.na(pc)) return(NA_character_)
    if (pc >= 60 && (is.null(sc) || sc %in% c("robust", "ns")))
      return("trustworthy")
    if (pc < 30 || ps > 50)
      return("unreliable")
    "questionable"
  },
  pct_concordant,
  pct_suppressed,
  if (is.null(status_change)) rep(list(NULL), length(pct_concordant))
  else status_change,
  USE.NAMES = FALSE)
}
