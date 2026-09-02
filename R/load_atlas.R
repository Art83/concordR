# =============================================================================
# load_atlas, cached atlas loader
# =============================================================================

# Module-level cache (package environment)
.atlas_env <- new.env(parent = emptyenv())


#' Load the RNA-protein concordance atlas
#'
#' Loads the precomputed atlas into a session cache. Subsequent calls
#' within the same R session return the cached object without reloading.
#' All other concordR functions call this internally when needed, so
#' explicit use is optional. Both atlases ship as internal package data
#' (\code{R/sysdata.rda}) and require no network access.
#'
#' @section Atlas contents (gene-level):
#' One row per gene (11,151 genes) with columns including
#' \code{gene_symbol}, \code{ensembl_id}, \code{gene_class},
#' \code{protein_confidence}, \code{detection_rate},
#' \code{suppression_rate}, \code{mean_rna_rank}, \code{n_tissues},
#' \code{mechanism_tier}, \code{tau} (RNA, evidence only),
#' \code{protein_breadth} (evidence only), \code{protein_tau} (IHC panel,
#' the specificity gate input), \code{primary_tissue_ihc},
#' \code{n_tissues_ihc}, and \code{n_tissues_ihc_detected}. Localisation
#' flags are not in this table; \code{\link{triage}} joins them at run time
#' from a separate UniProt-derived object.
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
#'   \code{"hpa"} loads the three-source protein presence vote (HPA IHC,
#'   GTEx MS, PaxDb) over 59 HPA tissues, which the tissue gate reads.
#'   \code{"both"} loads all three.
#' @param reload Logical. Force reload from package data. Default \code{FALSE}.
#'
#' @return For \code{"gene"}: a data frame (one row per gene). For
#'   \code{"tissue"}: a data frame (one row per gene x tissue). For
#'   \code{"hpa"}: a data frame (one row per gene x tissue, full panel).
#'   For \code{"both"}: a named list with elements \code{gene},
#'   \code{tissue}, and \code{hpa}.
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
load_atlas <- function(component = c("gene", "tissue", "hpa", "both"),
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

  if (component %in% c("hpa", "both")) {
    if (isTRUE(reload) || is.null(.atlas_env$hpa)) {
      if (exists("hpa_tissue", envir = asNamespace("concordR"),
                 inherits = FALSE)) {
        .atlas_env$hpa <- get("hpa_tissue",
                              envir = asNamespace("concordR"))
        if (inherits(.atlas_env$hpa, "data.table"))
          .atlas_env$hpa <- as.data.frame(.atlas_env$hpa)
      } else {
        .atlas_env$hpa <- NULL
      }
    }
  }

  switch(component,
    gene   = .atlas_env$gene,
    tissue = .atlas_env$tissue,
    hpa = .atlas_env$hpa,
    both = list(gene = .atlas_env$gene, tissue = .atlas_env$tissue,
                  hpa = .atlas_env$hpa)
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
  if (inherits(.atlas_env$tissue, "data.table"))
    .atlas_env$tissue <- as.data.frame(.atlas_env$tissue)
  .atlas_env$tissue
}


#' Internal helper: ensure HPA full tissue panel is loaded
#' @keywords internal
#' @noRd
.ensure_hpa_tissue <- function() {
  if (is.null(.atlas_env$hpa)) load_atlas("hpa")
  if (inherits(.atlas_env$hpa, "data.table"))
    .atlas_env$hpa <- as.data.frame(.atlas_env$hpa)
  .atlas_env$hpa
}


#' Internal helper: ensure localisation flags are loaded
#' @keywords internal
#' @noRd
.ensure_localisation_flags <- function() {
  if (is.null(.atlas_env$loc)) {
    if (exists("localisation_flags", envir = asNamespace("concordR"),
               inherits = FALSE)) {
      .atlas_env$loc <- get("localisation_flags",
                            envir = asNamespace("concordR"))
      if (inherits(.atlas_env$loc, "data.table"))
        .atlas_env$loc <- as.data.frame(.atlas_env$loc)
    } else {
      .atlas_env$loc <- NULL
    }
  }
  .atlas_env$loc
}
