# =============================================================================
# flag_genes — pre-analysis sanity check
# =============================================================================

#' Flag genes by concordance reliability and compartment plausibility
#'
#' The first function to run on any proteomics gene list. For each gene,
#' reports whether RNA-based evidence for it is trustworthy, whether its
#' subcellular localisation is consistent with the measurement compartment,
#' and whether its tissue expression is broad enough to undermine
#' tissue-of-origin claims.
#'
#' @section Traffic-light logic:
#' Each gene receives a \code{flag} combining concordance and localisation:
#' \describe{
#'   \item{\code{"reliable"}}{Concordant gene with plausible localisation.
#'     RNA evidence is trustworthy for this gene in this context.}
#'   \item{\code{"caution"}}{Variable gene, or concordant with borderline
#'     localisation. Interpret with care.}
#'   \item{\code{"unreliable"}}{Suppressed gene and/or implausible
#'     localisation. RNA evidence is uninformative for this gene.
#'     Protein-level measurements are required.}
#'   \item{\code{NA}}{Gene not in atlas.}
#' }
#'
#' Projection proteins (axonal/dendritic) in biofluid contexts (plasma,
#' serum, CSF) receive a \code{"caution"} rather than \code{"unreliable"}
#' flag because damage-release biology (e.g. NEFL in plasma after axonal
#' injury) is a legitimate detection pathway.
#'
#' @param genes Character vector of HGNC gene symbols.
#' @param sample_type Character. The biological compartment where proteins
#'   were measured. One of \code{"plasma"}, \code{"serum"}, \code{"csf"},
#'   \code{"urine"}, \code{"cell_lysate"}, \code{"tissue"},
#'   \code{"nuclear_fraction"}, \code{"mitochondrial_fraction"},
#'   \code{"membrane_fraction"}. Determines which localisations are
#'   considered plausible.
#'
#' @return A data frame with one row per gene, containing:
#' \describe{
#'   \item{\code{gene_symbol}}{Gene symbol.}
#'   \item{\code{gene_class}}{From the atlas.}
#'   \item{\code{protein_confidence}}{From the atlas.}
#'   \item{\code{detection_rate}}{From the atlas.}
#'   \item{\code{mechanism_tier}}{From the atlas (NA for non-suppressed).}
#'   \item{\code{localisation}}{Human-readable localisation label.}
#'   \item{\code{loc_plausible}}{Logical. Is the localisation consistent
#'     with \code{sample_type}?}
#'   \item{\code{tau}}{Yanai tissue-specificity index.}
#'   \item{\code{specificity_class}}{Categorical specificity band:
#'     tissue_specific / tissue_enriched / mixed / broadly_expressed.}
#'   \item{\code{flag}}{Traffic light: \code{"reliable"}, \code{"caution"},
#'     or \code{"unreliable"}.}
#'   \item{\code{caveat}}{Non-NA where additional interpretive context
#'     applies (projection proteins, broadly expressed genes).}
#' }
#'
#' @examples
#' \dontrun{
#' # Plasma proteomics panel
#' flag_genes(
#'   genes       = c("NEFL", "GFAP", "TREM2", "NRGN", "HDAC1", "TP53"),
#'   sample_type = "plasma"
#' )
#'
#' # CSF biomarker candidates
#' flag_genes(
#'   genes       = c("NEFL", "VILIP1", "NRGN", "CST3"),
#'   sample_type = "csf"
#' )
#' }
#'
#' @seealso \code{\link{query_atlas}}, \code{\link{triage}}
#' @export
flag_genes <- function(genes, sample_type) {
  
  if (!is.character(genes) || length(genes) == 0L)
    stop("`genes` must be a non-empty character vector.")
  
  sample_type <- tolower(sample_type)
  if (!sample_type %in% .SAMPLE_TYPES)
    stop("Unknown `sample_type`: '", sample_type, "'. Available: ",
         paste(.SAMPLE_TYPES, collapse = ", "), ".")
  
  # --- Query atlas ---
  q <- query_atlas(genes, missing = "keep")
  
  # --- Localisation plausibility (found rows only) ---
  q$loc_plausible <- NA
  found_idx <- which(q$found)
  if (length(found_idx) > 0L) {
    q$loc_plausible[found_idx] <-
      .check_plausibility(q[found_idx, , drop = FALSE], sample_type)
  }
  
  # --- Human-readable localisation ---
  q$localisation <- NA_character_
  if (length(found_idx) > 0L) {
    q$localisation[found_idx] <-
      .location_labels(q[found_idx, , drop = FALSE])
  }
  
  # --- Projection status (for damage-release caveat logic) ---
  proj <- if ("is_projection" %in% names(q))
    !is.na(q$is_projection) & q$is_projection == 1
  else rep(FALSE, nrow(q))
  
  # --- Traffic light ---
  damage_context <- .is_damage_release_context(sample_type)
  q$flag <- .assign_flag(q$gene_class, q$loc_plausible, proj, damage_context)
  
  # --- Caveats ---
  q$caveat <- .build_caveat(q, proj, damage_context)
  
  # --- Select + order output columns ---
  candidate_cols <- c("gene_symbol", "gene_class", "protein_confidence",
                      "detection_rate", "mechanism_tier",
                      "localisation", "loc_plausible",
                      "tau", "specificity_class",
                      "flag", "caveat")
  out_cols <- intersect(candidate_cols, names(q))
  q[, out_cols, drop = FALSE]
}


#' Assign traffic-light flag from gene class, localisation, and context
#'
#' Projection proteins in damage-release contexts (plasma/serum/CSF) are
#' downgraded to "caution" rather than "unreliable" — NEFL-in-plasma is
#' a legitimate biomarker pattern.
#' @keywords internal
#' @noRd
.assign_flag <- function(gene_class, loc_plausible, is_projection,
                         damage_context) {
  mapply(function(gc, lp, proj) {
    if (is.na(gc)) return(NA_character_)
    
    if (gc == "consistently_suppressed") {
      if (proj) return("caution")
      return("unreliable")
    }
    
    # Implausible localisation in sample compartment
    if (!is.na(lp) && !lp) {
      # Damage-release context salvages projection proteins
      if (proj && damage_context) return("caution")
      return("unreliable")
    }
    
    if (gc == "consistently_concordant") {
      if (is.na(lp) || lp) return("reliable")
      return("caution")
    }
    
    if (gc == "variable")       return("caution")
    if (gc == "low_expression") return("caution")
    
    NA_character_
  }, gene_class, loc_plausible, is_projection, USE.NAMES = FALSE)
}


#' Build per-gene caveat text combining multiple interpretive flags
#' @keywords internal
#' @noRd
.build_caveat <- function(q, proj, damage_context) {
  n <- nrow(q)
  caveats <- rep(NA_character_, n)
  
  for (i in seq_len(n)) {
    msgs <- character(0)
    
    # Projection + suppressed: assay limitation
    if (!is.na(q$gene_class[i]) &&
        q$gene_class[i] == "consistently_suppressed" && proj[i]) {
      msgs <- c(msgs,
                "Localises to cell projections/axons — IHC may under-detect.")
    }
    
    # Projection + damage-release context: biomarker-via-damage
    if (proj[i] && damage_context &&
        !is.na(q$loc_plausible[i]) && !q$loc_plausible[i]) {
      msgs <- c(msgs,
                "Projection protein detectable via damage-release in biofluids.")
    }
    
    # Broadly expressed: tissue-of-origin claims unsupported
    if ("specificity_class" %in% names(q) &&
        !is.na(q$specificity_class[i]) &&
        q$specificity_class[i] == "broadly_expressed") {
      msgs <- c(msgs,
                "Broadly expressed (tau low); tissue-of-origin claims unsupported.")
    }
    
    if (length(msgs) > 0L) caveats[i] <- paste(msgs, collapse = " ")
  }
  
  caveats
}
