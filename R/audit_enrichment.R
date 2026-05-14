# =============================================================================
# audit_enrichment — concordance-corrected enrichment QC
# =============================================================================

#' Audit enrichment results against the concordance atlas
#'
#' Takes any standard enrichment result — from clusterProfiler, enrichR,
#' fgsea, g:Profiler, or a simple data frame — and recomputes significance
#' after removing genes with poor RNA--protein concordance. Also surfaces
#' the fraction of pathway-driving genes that are broadly expressed across
#' tissues, and (when \code{context_tissue} is supplied) the fraction
#' detected at protein level in the claimed tissue.
#'
#' @section What this function does:
#' For each enriched pathway:
#' \enumerate{
#'   \item Identifies which genes in the pathway overlap the query list.
#'   \item Looks up each overlap gene's concordance class and localisation.
#'   \item Computes: fraction concordant, fraction suppressed, fraction
#'     localisation-plausible, fraction broadly expressed, and (if
#'     \code{context_tissue} supplied) fraction detected at protein level
#'     in the claimed tissue.
#'   \item Reruns the hypergeometric test after removing suppressed-class
#'     genes from both the overlap and the background.
#'   \item Reports both original and corrected p-values.
#' }
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
#' @param context_tissue Character or \code{NULL}. The tissue claimed as
#'   the source of the enrichment (e.g. \code{"Cerebellum"}). When
#'   supplied, adds \code{pct_detected_in_context_tissue} column and
#'   downgrades the verdict to \code{"tissue_unsupported"} for pathways
#'   where the majority of overlap genes are not detected at protein
#'   level in the claimed tissue. Brain is resolved as four regions —
#'   pass one explicitly. See \code{\link{triage}} for details.
#' @param method Character. \code{"remove"} (default) excludes suppressed
#'   genes when recomputing the hypergeometric test.
#' @param alpha Numeric. Significance threshold. Default \code{0.05}.
#'
#' @return The original data frame with additional columns:
#'   \code{n_overlap}, \code{n_concordant}, \code{n_suppressed},
#'   \code{n_variable}, \code{pct_concordant}, \code{pct_suppressed},
#'   \code{pct_broadly_expressed}, \code{pct_loc_plausible} (if
#'   \code{sample_type} supplied), \code{pct_detected_in_context_tissue}
#'   (if \code{context_tissue} supplied), \code{p_corrected},
#'   \code{padj_corrected}, \code{status_change}, \code{concordance_verdict}.
#'
#' @examples
#' \dontrun{
#' audited <- audit_enrichment(
#'   enrichment_result = as.data.frame(ego),
#'   gene_list         = my_genes,
#'   gene_sets         = my_gene_sets,
#'   universe          = all_genes,
#'   sample_type       = "csf",
#'   context_tissue    = "Cerebellum"
#' )
#' audited[audited$concordance_verdict == "tissue_unsupported",
#'         c("Description", "pct_detected_in_context_tissue")]
#' }
#'
#' @seealso \code{\link{plot_audit}}, \code{\link{flag_genes}},
#'   \code{\link{triage}}
#' @export
audit_enrichment <- function(enrichment_result,
                             gene_list,
                             gene_sets      = NULL,
                             universe       = NULL,
                             sample_type    = NULL,
                             context_tissue = NULL,
                             method         = c("remove"),
                             alpha          = 0.05) {

  method <- match.arg(method)

  if (!is.data.frame(enrichment_result) || nrow(enrichment_result) == 0L)
    stop("`enrichment_result` must be a non-empty data frame.")
  if (!is.character(gene_list) || length(gene_list) == 0L)
    stop("`gene_list` must be a non-empty character vector.")

  atlas <- .ensure_atlas()

  # --- Resolve tissue context (lazy: only loads tissue atlas if asked) ---
  tissue_atlas    <- NULL
  resolved_tissue <- NULL
  if (!is.null(context_tissue)) {
    tissue_atlas    <- .ensure_tissue_atlas()
    resolved_tissue <- .validate_context_tissue(context_tissue, tissue_atlas)
  }

  # --- Detect columns ---
  pathway_col <- .detect_col(enrichment_result,
                             c("ID", "term", "pathway", "Description", "set", "Term_name",
                               "term_name", "source_term", "GO.ID", "category"),
                             "pathway/term")

  # Detect raw and adjusted p-value columns separately. We use the
  # adjusted column for status comparison so that we are comparing like
  # to like (BH-adjusted original vs BH-adjusted corrected). Falling
  # back to raw if no adjusted column is present.
  padj_col <- .detect_col(enrichment_result,
                          c("p.adjust", "padj", "p.adj", "p_adj", "padjust",
                            "FDR", "qvalue", "q.value", "q_value"),
                          "adjusted p-value", required = FALSE)

  praw_col <- .detect_col(enrichment_result,
                          c("pvalue", "p_value", "p.value", "P.value",
                            "Pvalue", "PValue", "nominal_p"),
                          "raw p-value", required = FALSE)

  # Column to use for the original-vs-corrected status comparison.
  pval_col <- if (!is.null(padj_col)) padj_col else praw_col

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
    .annotate_pathway(overlap, atlas, sample_type,
                      tissue_atlas, resolved_tissue)
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
    message("Corrected p-values require `gene_sets` (full pathway ",
            "memberships). The `geneID` column only contains overlap ",
            "genes. Skipping p-value correction; annotation columns ",
            "still added.")
  }

  # --- Concordance verdict ---
  ann_df$concordance_verdict <- .verdict(
    ann_df$pct_concordant,
    ann_df$pct_suppressed,
    if ("pct_broadly_expressed" %in% names(ann_df))
      ann_df$pct_broadly_expressed else rep(NA_real_, nrow(ann_df)),
    if ("status_change" %in% names(ann_df)) ann_df$status_change else NULL,
    if ("pct_detected_in_context_tissue" %in% names(ann_df))
      ann_df$pct_detected_in_context_tissue else rep(NA_real_, nrow(ann_df)),
    if ("mean_cross_platform_agree" %in% names(ann_df))
      ann_df$mean_cross_platform_agree else NULL
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
#'
#' If \code{tissue_atlas} and \code{resolved_tissue} are supplied, the
#' returned row includes \code{pct_detected_in_context_tissue}: among
#' overlap genes for which we have tissue-resolved data, the fraction
#' detected at protein level in the named tissue.
#'
#' @keywords internal
#' @noRd
.annotate_pathway <- function(overlap_genes, atlas, sample_type = NULL,
                              tissue_atlas = NULL, resolved_tissue = NULL) {
  n <- length(overlap_genes)

  empty_row <- function() {
    out <- data.frame(
      n_overlap = 0L, n_concordant = 0L, n_suppressed = 0L,
      n_variable = 0L,
      pct_concordant = NA_real_, pct_suppressed = NA_real_,
      pct_broadly_expressed = NA_real_,
      mean_protein_evidence = NA_real_,
      mean_cross_platform_agree = NA_real_,
      stringsAsFactors = FALSE)
    if (!is.null(sample_type))    out$pct_loc_plausible              <- NA_real_
    if (!is.null(resolved_tissue)) out$pct_detected_in_context_tissue <- NA_real_
    out
  }

  if (n == 0L) return(empty_row())

  idx <- match(toupper(overlap_genes), toupper(atlas$gene_symbol))
  matched <- atlas[idx[!is.na(idx)], , drop = FALSE]
  n_found <- nrow(matched)
  if (n_found == 0L) return(empty_row())

  classes <- matched$gene_class
  n_conc <- sum(classes == "consistently_concordant", na.rm = TRUE)
  n_supp <- sum(classes == "consistently_suppressed", na.rm = TRUE)
  n_var  <- sum(classes == "variable", na.rm = TRUE)

  # Tissue-specificity: fraction broadly-expressed
  if ("specificity_class" %in% names(matched)) {
    n_broad <- sum(matched$specificity_class == "broadly_expressed",
                   na.rm = TRUE)
    n_spec_known <- sum(!is.na(matched$specificity_class))
    pct_broad <- if (n_spec_known > 0) n_broad / n_spec_known * 100
    else NA_real_
  } else if ("tau" %in% names(matched)) {
    n_broad <- sum(matched$tau < 0.20, na.rm = TRUE)
    n_tau_known <- sum(!is.na(matched$tau))
    pct_broad <- if (n_tau_known > 0) n_broad / n_tau_known * 100
    else NA_real_
  } else {
    pct_broad <- NA_real_
  }
  
  out <- data.frame(
    n_overlap             = n,
    n_concordant          = n_conc,
    n_suppressed          = n_supp,
    n_variable            = n_var,
    pct_concordant        = n_conc / n_found * 100,
    pct_suppressed        = n_supp / n_found * 100,
    pct_broadly_expressed = pct_broad,
    stringsAsFactors = FALSE)
  
  # Multi-platform protein evidence
  if ("mean_protein_evidence" %in% names(matched)) {
    out$mean_protein_evidence <- mean(matched$mean_protein_evidence, na.rm = TRUE)
  }
  if ("cross_platform_agree_rate" %in% names(matched)) {
    out$mean_cross_platform_agree <- mean(matched$cross_platform_agree_rate, na.rm = TRUE)
  }

  if (!is.null(sample_type)) {
    plausible <- .check_plausibility(matched, sample_type)
    out$pct_loc_plausible <- if (any(!is.na(plausible)))
      sum(plausible, na.rm = TRUE) / sum(!is.na(plausible)) * 100
    else NA_real_
  }

  # --- Tissue-of-origin coherence ---
  if (!is.null(resolved_tissue) && !is.null(tissue_atlas)) {
    out$pct_detected_in_context_tissue <-
      .pct_detected_in_tissue(overlap_genes, tissue_atlas, resolved_tissue)
  }

  out
}


#' Fraction of overlap genes detected in the claimed tissue.
#'
#' "Detected" means the per-(gene, tissue) row exists in the tissue
#' atlas with \code{tissue_class == "concordant_in_tissue"} or
#' \code{tissue_class == "variable_in_tissue"} with non-zero
#' \code{detect_fraction}. "Not detected" means \code{low_expression},
#' \code{suppressed_in_tissue}, or \code{variable_in_tissue} with zero
#' detection. Genes absent from the tissue atlas are excluded from the
#' denominator.
#'
#' @keywords internal
#' @noRd
.pct_detected_in_tissue <- function(overlap_genes, tissue_atlas, ts_tissue) {
  ta_sub <- tissue_atlas[tissue_atlas$ts_tissue == ts_tissue, , drop = FALSE]
  idx <- match(toupper(overlap_genes), toupper(ta_sub$gene_symbol))
  matched <- !is.na(idx)
  if (!any(matched)) return(NA_real_)
  
  cls <- ta_sub$tissue_class[idx[matched]]
  det <- ta_sub$detect_fraction[idx[matched]]
  ms  <- if ("ms_detected" %in% names(ta_sub))
    ta_sub$ms_detected[idx[matched]] else rep(FALSE, sum(matched))
  
  # Detected by IHC OR by MS
  detected <- (cls == "concordant_in_tissue") |
    (cls == "variable_in_tissue" & !is.na(det) & det > 0) |
    (!is.na(ms) & ms)
  sum(detected, na.rm = TRUE) / sum(matched) * 100
}


#' Recompute hypergeometric p-value after concordance correction
#' @keywords internal
#' @noRd
.corrected_pvalue <- function(gene_list, pathway_genes, universe,
                              atlas, method) {

  gene_list_u     <- toupper(gene_list)
  pathway_genes_u <- toupper(pathway_genes)
  universe_u      <- toupper(universe)

  suppressed <- toupper(atlas$gene_symbol[
    !is.na(atlas$gene_class) &
      atlas$gene_class == "consistently_suppressed"])

  gene_list_c     <- setdiff(gene_list_u, suppressed)
  pathway_genes_c <- setdiff(pathway_genes_u, suppressed)
  universe_c      <- setdiff(universe_u, suppressed)

  overlap <- length(intersect(gene_list_c, pathway_genes_c))
  n_list  <- length(intersect(gene_list_c, universe_c))
  n_pw    <- length(intersect(pathway_genes_c, universe_c))
  n_univ  <- length(universe_c)

  if (n_list == 0L || n_pw == 0L || n_univ == 0L) return(1.0)
  phyper(overlap - 1L, n_pw, n_univ - n_pw, n_list, lower.tail = FALSE)
}


#' Concordance verdict factoring tissue specificity and tissue context
#'
#' Tissue-of-origin coherence (when supplied) takes priority: a pathway
#' where most overlap genes are not detected at protein level in the
#' claimed tissue is flagged as \code{tissue_unsupported} regardless
#' of the concordance picture. This is the case where statistical
#' enrichment exists but biological grounding does not.
#' @keywords internal
#' @noRd
.verdict <- function(pct_concordant, pct_suppressed,
                     pct_broadly_expressed, status_change = NULL,
                     pct_detected_in_context_tissue = NULL,
                     mean_cross_platform_agree = NULL) {
  n <- length(pct_concordant)
  out <- character(n)
  for (i in seq_len(n)) {
    pc <- pct_concordant[i]
    ps <- pct_suppressed[i]
    pb <- pct_broadly_expressed[i]
    sc <- if (is.null(status_change)) NA_character_ else status_change[i]
    pt <- if (is.null(pct_detected_in_context_tissue)) NA_real_
    else pct_detected_in_context_tissue[i]
    mcpa <- if (is.null(mean_cross_platform_agree)) NA_real_
    else mean_cross_platform_agree[i]
    
    if (is.na(pc)) { out[i] <- NA_character_; next }
    
    # --- Tissue context overrides (only when supplied) ---
    if (!is.na(pt) && pt < 50) {
      out[i] <- "tissue_unsupported"; next
    }
    
    # --- High suppression fraction ---
    if (pc < 30 || (!is.na(ps) && ps > 50)) {
      # Cross-platform disagreement weakens the suppression calls
      if (!is.na(mcpa) && mcpa < 0.3) {
        out[i] <- "questionable"; next
      }
      out[i] <- "unreliable"; next
    }
    
    # --- Concordance OK but broadly expressed ---
    if (!is.na(pb) && pb > 80) {
      out[i] <- "concordant_but_broad"; next
    }
    
    # --- Strong concordance, robust correction ---
    if (pc >= 60 && (is.na(sc) || sc %in% c("robust", "ns"))) {
      # High cross-platform agreement strengthens the call
      if (!is.na(mcpa) && mcpa > 0.7) {
        out[i] <- "trustworthy"; next
      }
      out[i] <- "trustworthy"; next
    }
    
    out[i] <- "questionable"
  }
  out
}