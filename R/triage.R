# =============================================================================
# triage — biomarker / target / enrichment hit classifier
# =============================================================================

#' Triage candidate genes by claimed role and measurement context
#'
#' For each gene, checks whether the claimed role — biomarker, therapeutic
#' target, or enrichment hit — is consistent with the concordance atlas.
#' A nuclear protein claimed as a plasma biomarker is flagged. A
#' suppressed-class gene claimed as a drug target from RNA evidence is
#' flagged.
#'
#' @section Verdict logic:
#'
#' **For \code{claim = "biomarker"}:**
#' \describe{
#'   \item{\code{"plausible"}}{Localisation consistent with sample type.
#'     Concordant or variable class.}
#'   \item{\code{"questionable"}}{Plausible localisation but suppressed
#'     class, or borderline. Requires protein-level validation.}
#'   \item{\code{"implausible"}}{Intracellular protein in extracellular
#'     fluid. Presence likely reflects cell death.}
#' }
#'
#' **For \code{claim = "target"}:**
#' \describe{
#'   \item{\code{"RNA-supported"}}{Concordant class. RNA reliably predicts
#'     protein for this gene.}
#'   \item{\code{"RNA-unreliable"}}{Suppressed class. Target may be valid
#'     but RNA-only evidence is not defensible.}
#'   \item{\code{"caution"}}{Variable class. Tissue-context-dependent.}
#' }
#'
#' **For \code{claim = "enrichment_hit"}:**
#' \describe{
#'   \item{\code{"supported"}}{Concordant and localisation-plausible.}
#'   \item{\code{"unreliable"}}{Suppressed or localisation-implausible.}
#' }
#'
#' @param genes Character vector of HGNC gene symbols.
#' @param claim Character. One of \code{"biomarker"}, \code{"target"}, or
#'   \code{"enrichment_hit"}.
#' @param sample_type Character. Required for \code{"biomarker"} and
#'   \code{"enrichment_hit"} claims. See \code{\link{flag_genes}}.
#' @param max_breadth Integer. Tissue breadth threshold. Default \code{10}.
#'
#' @return A data frame with atlas annotations plus \code{claim},
#'   \code{verdict}, and \code{reason} columns.
#'
#' @examples
#' \dontrun{
#' triage(c("NEFL", "GFAP", "NRGN", "HDAC1"), "biomarker", "plasma")
#' triage(c("AR", "ESR1", "CD274", "TP53"), "target")
#' }
#'
#' @seealso \code{\link{flag_genes}}, \code{\link{audit_enrichment}}
#' @export
triage <- function(genes,
                   claim,
                   sample_type = NULL,
                   max_breadth = 10L) {

  claim <- match.arg(claim, c("biomarker", "target", "enrichment_hit"))

  if (claim %in% c("biomarker", "enrichment_hit") && is.null(sample_type))
    stop("`sample_type` is required for '", claim, "' claims.")

  # --- Get atlas annotations ---
  q <- query_atlas(genes, missing = "keep")

  # --- Localisation plausibility ---
  if (!is.null(sample_type)) {
    q$loc_plausible <- .check_plausibility(q, sample_type)
  } else {
    q$loc_plausible <- NA
  }

  # --- Localisation label ---
  q$localisation <- .location_labels(q)

  # --- Tissue breadth ---
  q$broad_expression <- !is.na(q$n_tissues) & q$n_tissues > max_breadth

  # --- Biofluid interpretation (for biomarker claims) ---
  q$biofluid_note <- .biofluid_note(q)

  # --- Apply verdict logic ---
  q$claim <- claim
  verdicts <- mapply(.triage_one,
                     q$gene_class, q$loc_plausible, q$broad_expression,
                     MoreArgs = list(claim = claim),
                     SIMPLIFY = FALSE)

  q$verdict <- vapply(verdicts, `[[`, character(1), "verdict")
  q$reason  <- vapply(verdicts, `[[`, character(1), "reason")

  # --- Select output columns ---
  out_cols <- c("gene_symbol", "gene_class", "protein_confidence",
                "detection_rate", "localisation", "loc_plausible",
                "broad_expression", "biofluid_note",
                "claim", "verdict", "reason")
  out_cols <- intersect(out_cols, names(q))
  q[, out_cols, drop = FALSE]
}


#' Generate biofluid interpretation note (matches Streamlit logic)
#' @keywords internal
#' @noRd
.biofluid_note <- function(atlas_df) {
  vapply(seq_len(nrow(atlas_df)), function(i) {
    row <- atlas_df[i, ]
    if (is.na(row$found) || !row$found) return(NA_character_)

    is_sec  <- isTRUE(row$is_secreted == 1)
    is_mem  <- isTRUE(row$is_membrane == 1)
    is_nuc  <- isTRUE(row$is_nuclear == 1)
    is_cyto <- isTRUE(row$is_cytoplasmic == 1)
    is_mito <- isTRUE(row$is_mitochondrial == 1)
    conf    <- row$protein_confidence

    if (is_sec) return("Expected in biofluids (secreted)")
    if (is_mem) return("May appear via shedding/exosomes")
    if ((is_nuc || is_cyto) && !is.na(conf) && conf > 0.7)
      return("Biofluid elevation suggests tissue damage")
    if (is_mito && !is.na(conf) && conf > 0.7)
      return("Biofluid elevation suggests mitochondrial damage")
    "Context-dependent"
  }, character(1))
}


#' Triage logic for a single gene
#' @keywords internal
#' @noRd
.triage_one <- function(gc, lp, broad, claim) {
  if (is.na(gc))
    return(list(verdict = NA_character_, reason = "Gene not in atlas."))

  switch(claim,
    biomarker      = .triage_biomarker(gc, lp, broad),
    target         = .triage_target(gc),
    enrichment_hit = .triage_enrichment(gc, lp)
  )
}


#' @keywords internal
#' @noRd
.triage_biomarker <- function(gc, lp, broad) {
  if (!is.na(lp) && !lp) {
    reason <- "Intracellular protein in extracellular compartment."
    if (broad)
      reason <- paste(reason, "Expressed broadly; tissue-specific",
                      "provenance unsupported.")
    return(list(verdict = "implausible", reason = reason))
  }
  if (gc == "consistently_suppressed")
    return(list(verdict = "questionable",
                reason  = paste("Suppressed: RNA poorly predicts protein.",
                                "Protein-level validation required.")))
  if (broad && gc %in% c("consistently_concordant", "variable"))
    return(list(verdict = "questionable",
                reason  = paste("Expressed in many tissues;",
                                "tissue-specific provenance is weak.")))
  if (gc == "consistently_concordant")
    return(list(verdict = "plausible",
                reason  = "Concordant, localisation-consistent."))
  if (gc == "variable")
    return(list(verdict = "questionable",
                reason  = "Variable concordance; context-dependent."))
  list(verdict = "questionable", reason = "Low expression class.")
}


#' @keywords internal
#' @noRd
.triage_target <- function(gc) {
  switch(gc,
    consistently_concordant = list(
      verdict = "RNA-supported",
      reason  = paste("Concordant. RNA reliably predicts protein.",
                      "RNA-based target evidence is defensible.")),
    consistently_suppressed = list(
      verdict = "RNA-unreliable",
      reason  = paste("Suppressed. RNA does not predict protein.",
                      "Requires protein-level evidence.")),
    variable = list(
      verdict = "caution",
      reason  = paste("Variable. RNA-protein relationship is",
                      "tissue-dependent.")),
    list(verdict = "caution", reason = "Low expression class.")
  )
}


#' @keywords internal
#' @noRd
.triage_enrichment <- function(gc, lp) {
  if (gc == "consistently_suppressed" || (!is.na(lp) && !lp))
    return(list(verdict = "unreliable",
                reason  = "Suppressed and/or localisation-implausible."))
  if (gc == "consistently_concordant" && (is.na(lp) || lp))
    return(list(verdict = "supported",
                reason  = "Concordant, localisation-consistent."))
  list(verdict = "supported",
       reason  = "Variable class; context-dependent.")
}
