# =============================================================================
# load_atlas — session-cached atlas loader
# =============================================================================

# Module-level cache (package environment)
.atlas_env <- new.env(parent = emptyenv())


#' Load the RNA--protein concordance atlas
#'
#' Loads the precomputed atlas into a session cache. Subsequent calls
#' within the same R session return the cached object without reloading.
#' All other concordR functions call this internally when needed, so
#' explicit use is optional. Both atlases ship as internal package data
#' (\code{R/sysdata.rda}) and require no network access.
#'
#' @section Atlas contents (gene-level):
#' One row per gene (up to 11,154 genes) with columns including
#' \code{gene_symbol}, \code{gene_class}, \code{protein_confidence},
#' \code{detection_rate}, \code{mean_rna_rank}, \code{n_tissues},
#' \code{mechanistic_tier}, \code{tau}, \code{specificity_class}, and
#' the \code{is_*} subcellular localisation flags. See package
#' documentation for full column descriptions.
#'
#' @section Atlas contents (tissue-resolved):
#' One row per (gene, tissue) pair with collapsed cell-type statistics:
#' \describe{
#'   \item{\code{gene_symbol}}{HGNC symbol.}
#'   \item{\code{ts_tissue}}{Tissue label. Brain is resolved into four
#'     regions (Caudate, Cerebellum, Cerebral cortex, Hippocampus); the
#'     package does not aggregate them.}
#'   \item{\code{rna_rank_max}}{Max RNA rank across cell types in the
#'     tissue.}
#'   \item{\code{rna_mean_tissue}}{Mean RNA expression across cell types.}
#'   \item{\code{detect_fraction}}{Fraction of cell types in the tissue
#'     where protein was detected.}
#'   \item{\code{protein_max}}{Max protein score across cell types.}
#'   \item{\code{n_cell_types}}{Number of cell types in the tissue.}
#'   \item{\code{is_brain}}{Convenience flag, 1 for the four brain regions.}
#'   \item{\code{tissue_class}}{Per-tissue concordance class:
#'     \code{"low_expression"}, \code{"suppressed_in_tissue"},
#'     \code{"variable_in_tissue"}, \code{"concordant_in_tissue"}.}
#' }
#'
#' @param component Character. \code{"gene"} (default) loads the gene-level
#'   summary atlas. \code{"tissue"} loads the tissue-resolved atlas.
#'   \code{"both"} loads both.
#' @param reload Logical. Force reload from package data. Default \code{FALSE}.
#'
#' @return For \code{"gene"}: a data frame (one row per gene). For
#'   \code{"tissue"}: a data frame (one row per gene x tissue). For
#'   \code{"both"}: a named list with elements \code{gene} and \code{tissue}.
#'
#' @examples
#' \dontrun{
#' atlas <- load_atlas()
#' head(atlas)
#'
#' # Tissue-resolved atlas
#' tissue <- load_atlas("tissue")
#' head(tissue[tissue$gene_symbol == "NEFL", ])
#' }
#'
#' @export
load_atlas <- function(component = c("gene", "tissue", "both"),
                       reload    = FALSE) {
  component <- match.arg(component)

  if (component %in% c("gene", "both")) {
    if (isTRUE(reload) || is.null(.atlas_env$gene)) {
      if (exists("concordance_atlas", envir = asNamespace("concordR"),
                 inherits = FALSE)) {
        .atlas_env$gene <- get("concordance_atlas",
                               envir = asNamespace("concordR"))
      } else {
        stop("Gene-level atlas not found in package data. ",
             "Reinstall the package.", call. = FALSE)
      }
    }
  }

  if (component %in% c("tissue", "both")) {
    if (isTRUE(reload) || is.null(.atlas_env$tissue)) {
      if (exists("tissue_atlas", envir = asNamespace("concordR"),
                 inherits = FALSE)) {
        .atlas_env$tissue <- get("tissue_atlas",
                                 envir = asNamespace("concordR"))
      } else {
        stop("Tissue-resolved atlas not found in package data. ",
             "Reinstall the package; if you have just rebuilt sysdata.rda, ",
             "restart R.", call. = FALSE)
      }
    }
  }

  switch(component,
    gene   = .atlas_env$gene,
    tissue = .atlas_env$tissue,
    both   = list(gene = .atlas_env$gene, tissue = .atlas_env$tissue)
  )
}


#' Internal helper: ensure gene-level atlas is loaded
#' @keywords internal
#' @noRd
.ensure_atlas <- function() {
  if (is.null(.atlas_env$gene)) load_atlas("gene")
  .atlas_env$gene
}


#' Internal helper: ensure tissue-resolved atlas is loaded
#' @keywords internal
#' @noRd
.ensure_tissue_atlas <- function() {
  if (is.null(.atlas_env$tissue)) load_atlas("tissue")
  # Coerce to plain data.frame: tissue_atlas may have been built with
  # data.table, whose `[i, , drop = FALSE]` semantics differ from
  # data.frame and produce wrong subsets in downstream code.
  if (inherits(.atlas_env$tissue, "data.table"))
    .atlas_env$tissue <- as.data.frame(.atlas_env$tissue)
  .atlas_env$tissue
}
