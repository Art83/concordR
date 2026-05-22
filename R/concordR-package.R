#' concordR: Concordance-Aware Quality Control for Proteomic Interpretation
#'
#' An atlas of RNA-protein concordance across human tissues, packaged as
#' a practical QC layer for biofluid proteomic studies. Five functions,
#' one job: given a gene list from a proteomics experiment, tell the user
#' what they can and cannot conclude from it.
#'
#' The core workflow: \code{\link{triage}()} evaluates drug-target
#' defensibility through four gates (compartment plausibility,
#' tissue-of-origin coherence, RNA--protein concordance, and tissue
#' specificity). \code{\link{query_atlas}()} and
#' \code{\link{query_tissue}()} provide direct atlas lookups.
#' \code{\link{plot_gene}()} visualises tissue-resolved concordance
#' profiles. \code{\link{load_atlas}()} manages session-cached data
#' loading.
#'
#' Both the gene-level and tissue-resolved atlases ship as internal
#' package data and require no network access.
#'
#' @name concordR-package
#' @keywords internal
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom graphics barplot par abline legend text axis mtext points segments
#' @importFrom grDevices rgb
"_PACKAGE"
