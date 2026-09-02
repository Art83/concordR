#' concordR: Protein-Level Audit of Gene Lists
#'
#' Audits any gene list at the protein level against an integrated atlas of
#' RNA to protein concordance. Given proteins nominated as drug targets or
#' mechanisms from a biofluid or tissue sample, \code{\link{triage}()} scores
#' each on four independent gates and returns one column per gate:
#' concordance, tissue presence, biofluid plausibility, and specificity.
#' \code{\link{audit_call}()} and \code{\link{funnel_class}()} combine the
#' gates. \code{\link{query_atlas}()} and \code{\link{query_tissue}()} are
#' direct lookups, \code{\link{plot_gene}()} draws a tissue profile, and
#' \code{\link{load_atlas}()} manages the session cache.
#'
#' All atlas objects and thresholds ship as internal package data and need
#' no network access.
#'
#' @name concordR-package
#' @keywords internal
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom graphics barplot par abline legend text axis mtext points segments
#' @importFrom grDevices rgb
"_PACKAGE"
