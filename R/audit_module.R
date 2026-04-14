# =============================================================================
# audit_module — WGCNA module quality check
# =============================================================================

#' Assess WGCNA module coherence against the concordance atlas
#'
#' Takes a WGCNA module gene list and evaluates whether the co-expression
#' pattern at RNA level is likely to reflect co-regulation at protein level.
#' A module where 60\% of genes are suppressed-class means the RNA
#' correlation structure does not predict protein abundance — interpreting
#' that module as a "biological process" is misleading.
#'
#' @section Module quality dimensions:
#' \describe{
#'   \item{Concordance composition}{Fraction of module genes in each
#'     concordance class.}
#'   \item{Compartment heterogeneity}{How many distinct subcellular
#'     compartments are represented.}
#'   \item{Mechanistic coherence}{Among suppressed genes, what fraction
#'     share the same mechanistic tier.}
#' }
#'
#' @param genes Character vector of gene symbols in the module.
#' @param module_name Character. Optional label. Default \code{NULL}.
#' @param sample_type Character or \code{NULL}. If supplied, adds
#'   localisation plausibility.
#'
#' @return A list with class \code{"concordR_module_audit"} containing:
#'   \code{module_name}, \code{n_genes}, \code{n_found},
#'   \code{class_composition}, \code{compartment_composition},
#'   \code{n_compartments}, \code{mechanistic_tiers},
#'   \code{quality_score}, \code{verdict}, \code{gene_table}.
#'
#' @examples
#' \dontrun{
#' result <- audit_module(c("NEFL", "GFAP", "TREM2", "HDAC1", "TP53"),
#'                        module_name = "blue")
#' print(result)
#' }
#'
#' @seealso \code{\link{query_atlas}}, \code{\link{audit_enrichment}}
#' @export
audit_module <- function(genes,
                         module_name = NULL,
                         sample_type = NULL) {

  if (!is.character(genes) || length(genes) == 0L)
    stop("`genes` must be a non-empty character vector.")

  atlas <- .ensure_atlas()

  # --- Query atlas ---
  q <- query_atlas(genes, missing = "keep")
  found <- q[q$found, , drop = FALSE]
  n_found <- nrow(found)

  # --- Concordance class composition ---
  class_levels <- c("consistently_concordant", "variable",
                    "consistently_suppressed", "low_expression")
  if (n_found > 0) {
    class_tab <- table(factor(found$gene_class, levels = class_levels))
    class_frac <- as.numeric(class_tab) / n_found
    names(class_frac) <- class_levels
  } else {
    class_frac <- setNames(rep(NA_real_, 4), class_levels)
  }

  # --- Compartment heterogeneity ---
  comp <- .compartment_summary(found)
  n_compartments <- sum(comp > 0, na.rm = TRUE)

  # --- Mechanistic tiers (suppressed genes only) ---
  tier_levels <- c("mRNA_decay", "translational_block",
                   "post_translational", "unexplained")
  suppressed <- found[!is.na(found$gene_class) &
                        found$gene_class == "consistently_suppressed", ,
                      drop = FALSE]
  if (nrow(suppressed) > 0 && "mechanistic_tier" %in% names(suppressed)) {
    tier_tab <- table(factor(suppressed$mechanistic_tier,
                             levels = tier_levels))
    tier_frac <- as.numeric(tier_tab) / nrow(suppressed)
    names(tier_frac) <- tier_levels
  } else {
    tier_frac <- setNames(rep(NA_real_, 4), tier_levels)
  }

  # --- Quality score ---
  pct_conc <- class_frac["consistently_concordant"]
  pct_supp <- class_frac["consistently_suppressed"]
  compartment_focus <- if (n_compartments > 0) 1 / n_compartments else 0

  quality <- if (is.na(pct_conc)) NA_real_
             else pct_conc * 0.6 + (1 - pct_supp) * 0.2 +
                  compartment_focus * 0.2

  verdict <- if (is.na(quality)) NA_character_
             else if (quality >= 0.5) "interpretable"
             else if (quality >= 0.3) "mixed"
             else "misleading"

  result <- list(
    module_name             = module_name %||% "unnamed",
    n_genes                 = length(genes),
    n_found                 = n_found,
    class_composition       = class_frac,
    compartment_composition = comp,
    n_compartments          = n_compartments,
    mechanistic_tiers       = tier_frac,
    quality_score           = quality,
    verdict                 = verdict,
    gene_table              = q
  )

  class(result) <- "concordR_module_audit"
  result
}


#' @export
print.concordR_module_audit <- function(x, ...) {
  cat("Module:", x$module_name, "\n")
  cat("Genes:", x$n_genes, "(", x$n_found, "in atlas)\n")
  cat("\nConcordance composition:\n")
  labels <- c(consistently_concordant = "concordant",
              variable = "variable",
              consistently_suppressed = "suppressed",
              low_expression = "low expression")
  for (nm in names(x$class_composition)) {
    lab <- if (nm %in% names(labels)) labels[nm] else nm
    cat(sprintf("  %-16s %5.1f%%\n", lab, x$class_composition[nm] * 100))
  }
  cat("\nCompartments represented:", x$n_compartments, "\n")
  cat("Quality score:", sprintf("%.2f", x$quality_score), "\n")
  cat("Verdict:", toupper(x$verdict), "\n")
  invisible(x)
}


#' Null-coalescing operator
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
