#' concordR: Concordance-Aware Quality Control for Proteomic Interpretation
#'
#' An atlas of RNA--protein concordance across human tissues, packaged as
#' a practical QC layer for proteomic studies. Eight functions, one job:
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
#'   \item \code{\link{audit_module}} — assess WGCNA module coherence
#'     against the atlas.
#' }
#'
#' **Interpret your candidates:**
#' \itemize{
#'   \item \code{\link{triage}} — classify genes as biomarker-plausible,
#'     target-RNA-supported, or artifact-likely given a claimed role and
#'     measurement context.
#' }
#'
#' **Visualise:**
#' \itemize{
#'   \item \code{\link{plot_gene}} — tissue-resolved concordance profile.
#'   \item \code{\link{plot_audit}} — standard vs corrected enrichment.
#' }
#'
#' @section Atlas data:
#' The atlas is loaded once via \code{\link{load_atlas}} and cached for the
#' session. It ships as package data for genes with complete annotation;
#' tissue-resolved profiles are loaded on first use from Zenodo.
#'
#' @docType package
#' @name concordR-package
#' @keywords internal
#' @importFrom stats setNames median quantile p.adjust phyper sd cor
#' @importFrom utils head
#' @importFrom graphics barplot par abline legend text axis mtext rect
#' @importFrom grDevices colorRampPalette rgb
"_PACKAGE"
