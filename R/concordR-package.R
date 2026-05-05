#' concordR: Concordance-Aware Quality Control for Proteomic Interpretation
#'
#' An atlas of RNA--protein concordance across human tissues, packaged as
#' a practical QC layer for proteomic studies. Seven functions, one job:
#' given a gene list from a proteomics experiment, tell the user what they
#' can and cannot conclude from it.
#'
#' @section Core functions:
#'
#' **Understand your gene list:**
#' \itemize{
#'   \item \code{\link{query_atlas}} — look up concordance class, suppression
#'     rate, mechanistic tier, and localisation for any gene.
#'   \item \code{\link{flag_genes}} — pre-analysis sanity check: classify
#'     genes by concordance reliability and compartment plausibility
#'     using binary localisation flags (is_secreted, is_membrane, etc.).
#' }
#'
#' **Audit your results:**
#' \itemize{
#'   \item \code{\link{audit_enrichment}} — recompute enrichment with
#'     concordance weights; show what survives correction.
#' }
#'
#' **Interpret your candidates:**
#' \itemize{
#'   \item \code{\link{triage}} — drug-target defensibility audit:
#'     concordance, compartment plausibility, tissue-of-origin coherence,
#'     and tissue specificity gates applied per gene.
#' }
#'
#' **Visualise:**
#' \itemize{
#'   \item \code{\link{plot_gene}} — tissue-resolved concordance profile.
#'   \item \code{\link{plot_audit}} — standard vs corrected enrichment.
#' }
#'
#' @section Atlas data:
#' Both the gene-level and tissue-resolved atlases ship as internal
#' package data (\code{R/sysdata.rda}) and are loaded once per session
#' via \code{\link{load_atlas}}. No network access is required.
#'
#' @name concordR-package
#' @keywords internal
#' @importFrom stats setNames median quantile p.adjust phyper sd cor
#' @importFrom utils head
#' @importFrom graphics barplot par abline legend text axis mtext rect
#' @importFrom grDevices colorRampPalette rgb
"_PACKAGE"
