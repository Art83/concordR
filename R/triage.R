# =============================================================================
# triage — drug target defensibility
# =============================================================================

#' Triage candidate genes for drug target defensibility
#'
#' Evaluates whether a gene's nomination as a drug target is defensible
#' given the concordance atlas. Drug-target claims demand a coherent
#' mechanistic chain: RNA must reliably predict protein, the protein
#' must be plausibly detected in the sample compartment that nominated
#' it, and (where supplied) the gene must be expressed at protein level
#' in the claimed tissue. Failure of any link breaks the chain.
#'
#' @section Why target-only:
#' Earlier versions of \code{triage()} offered separate \code{biomarker}
#' and \code{enrichment_hit} modes. Biomarker triage is now handled by
#' \code{\link{flag_genes}}. Pathway-level questions are handled by
#' \code{\link{audit_enrichment}}. \code{triage()} is reserved for the
#' strictest claim — that a gene should be drugged.
#'
#' @section Gates (in order of severity):
#' \enumerate{
#'   \item \strong{Compartment plausibility}: subcellular localisation
#'     inconsistent with sample compartment. A damage-release rescue
#'     applies in biofluid contexts (\code{plasma}, \code{serum},
#'     \code{csf}) for projection-localised proteins (\code{is_projection
#'     = 1}), because HPA microscopy systematically under-detects axonal
#'     and dendritic proteins.
#'   \item \strong{Tissue-of-origin coherence}: when \code{context_tissue}
#'     is supplied, gene must be detected at protein level in that tissue.
#'     Three failure modes are flagged: \code{low_expression} (no RNA, no
#'     protein), \code{suppressed_in_tissue} (RNA present but no protein),
#'     and \code{variable_in_tissue} with zero protein detection across
#'     all cell types in the tissue. Brain is resolved as four regions —
#'     pass one explicitly (e.g. \code{"Cerebellum"}), not \code{"Brain"}.
#'   \item \strong{Concordance class}: suppressed class — RNA does not
#'     predict protein.
#'   \item \strong{Tissue specificity}: concordant and plausible, but
#'     broadly expressed across tissues — tissue-of-origin claims (e.g.
#'     cell-type-specific drug targeting) are unsupported.
#' }
#'
#' Multiple gates may fail simultaneously. The reported \code{verdict}
#' is the most severe failed gate; the \code{reason} string accumulates
#' all triggered gates so the full picture is visible. The
#' \code{compartment_rescue} column flags genes whose compartment gate
#' was rescued by the projection-protein exception.
#'
#' @section Verdicts:
#' \describe{
#'   \item{\code{RNA-supported}}{All gates pass. RNA-based target
#'     evidence is defensible.}
#'   \item{\code{RNA-supported_but_broad}}{Concordant and plausible, but
#'     broadly expressed. Use for non-tissue-specific target claims only.}
#'   \item{\code{RNA-unreliable}}{Suppressed class. Target may be valid
#'     but RNA-based evidence is not defensible; protein-level evidence
#'     required.}
#'   \item{\code{tissue_origin_mismatch}}{Gene not expressed at protein
#'     level in claimed tissue.}
#'   \item{\code{compartment_implausible}}{Subcellular localisation
#'     inconsistent with sample compartment, with no rescue available.
#'     The proteomic observation that nominated this target is suspect.}
#'   \item{\code{caution}}{Variable class or low expression — sparse or
#'     tissue-context-dependent evidence.}
#'   \item{\code{NA}}{Gene not in atlas.}
#' }
#'
#' @param genes Character vector of HGNC gene symbols.
#' @param sample_type Character. The biological compartment where proteins
#'   were measured. See \code{\link{flag_genes}} for the full list.
#' @param context_tissue Character or \code{NULL}. The tissue claimed as
#'   the source of the proteomic signal (e.g. \code{"Cerebellum"}). When
#'   supplied, enables the tissue-of-origin gate using the tissue-resolved
#'   atlas. The atlas resolves brain into four regions (Caudate,
#'   Cerebellum, Cerebral cortex, Hippocampus); pass one explicitly.
#'   Genes absent from the tissue-resolved atlas are not gated.
#'
#' @return A data frame with atlas annotations plus
#'   \code{compartment_rescue}, \code{context_tissue_class},
#'   \code{context_tissue_detected}, \code{verdict} and \code{reason}.
#'   The \code{context_tissue_*} columns are \code{NA} when
#'   \code{context_tissue} is not supplied or the gene is absent from
#'   the tissue-resolved atlas.
#'
#' @examples
#' \dontrun{
#' # Drug target audit on a plasma proteomics nomination list
#' triage(c("AR", "ESR1", "CD274", "TP53", "HDAC1", "NEFL", "GFAP"),
#'        sample_type = "plasma")
#'
#' # CSF target nominations claimed to originate in the cerebellum
#' triage(c("NEFL", "GFAP", "SPANXA1"),
#'        sample_type    = "csf",
#'        context_tissue = "Cerebellum")
#' }
#'
#' @seealso \code{\link{flag_genes}} for biomarker-style compartment
#'   checks; \code{\link{audit_enrichment}} for pathway-level audits.
#' @export
triage <- function(genes, sample_type, context_tissue = NULL) {

  if (!is.character(genes) || length(genes) == 0L)
    stop("`genes` must be a non-empty character vector.")

  if (missing(sample_type) || is.null(sample_type) ||
      !is.character(sample_type) || length(sample_type) != 1L)
    stop("`sample_type` is required (single character). Available: ",
         paste(.SAMPLE_TYPES, collapse = ", "), ".")

  sample_type <- tolower(sample_type)
  if (!sample_type %in% .SAMPLE_TYPES)
    stop("Unknown `sample_type`: '", sample_type, "'. Available: ",
         paste(.SAMPLE_TYPES, collapse = ", "), ".")

  # --- Atlas annotations ---
  q <- query_atlas(genes, missing = "keep")
  found_idx <- which(q$found)

  # --- Compartment plausibility (found rows only) ---
  q$loc_plausible <- NA
  if (length(found_idx) > 0L) {
    q$loc_plausible[found_idx] <-
      .check_plausibility(q[found_idx, , drop = FALSE], sample_type)
  }

  # --- Atlas blind-spot rescue ---
  # Cellular proteomics atlases systematically under-detect proteins
  # whose primary location is outside the cell of origin: axonal/dendritic
  # projection proteins (NEFL, MAP2, SYN1) and highly-secreted proteins
  # (ALB, APOE, immunoglobulins). For these genes in biofluid contexts,
  # both the compartment gate and the concordance gate are unreliable —
  # the atlas is blind to the relevant compartment. We rescue both gates
  # uniformly. `compartment_rescue` flags the rescue for transparency.
  q$compartment_rescue <- FALSE
  q$atlas_blind_spot   <- FALSE
  is_biofluid <- sample_type %in% c("plasma", "serum", "csf")
  if (is_biofluid && length(found_idx) > 0L) {
    is_proj <- "is_projection" %in% names(q) &
               !is.na(q$is_projection) & q$is_projection == 1
    is_sec  <- "is_secreted"   %in% names(q) &
               !is.na(q$is_secreted)   & q$is_secreted   == 1
    blind_mask <- (is_proj | is_sec)
    blind_mask[is.na(blind_mask)] <- FALSE
    q$atlas_blind_spot[blind_mask] <- TRUE

    # Compartment rescue: only flag when compartment gate would have failed
    rescue_mask <- blind_mask & !is.na(q$loc_plausible) & !q$loc_plausible
    if (any(rescue_mask)) {
      q$compartment_rescue[rescue_mask] <- TRUE
      q$loc_plausible[rescue_mask]      <- TRUE
    }
  }

  # --- Localisation label ---
  q$localisation <- NA_character_
  if (length(found_idx) > 0L) {
    q$localisation[found_idx] <-
      .location_labels(q[found_idx, , drop = FALSE])
  }

  # --- Biofluid interpretation note ---
  q$biofluid_note <- .biofluid_note(q)

  # --- Broadly-expressed flag ---
  broad <- if ("specificity_class" %in% names(q)) {
    !is.na(q$specificity_class) & q$specificity_class == "broadly_expressed"
  } else if ("tau" %in% names(q)) {
    !is.na(q$tau) & q$tau < 0.20
  } else {
    rep(FALSE, nrow(q))
  }

  # --- Tissue-of-origin gate ---
  q$context_tissue_class    <- NA_character_
  q$context_tissue_detected <- NA
  tissue_mismatch           <- rep(FALSE, nrow(q))
  resolved_tissue           <- NA_character_

  if (!is.null(context_tissue)) {
    ta <- .ensure_tissue_atlas()
    resolved_tissue <- .validate_context_tissue(context_tissue, ta)

    ta_sub <- ta[ta$ts_tissue == resolved_tissue, , drop = FALSE]
    idx <- match(toupper(q$gene_symbol), toupper(ta_sub$gene_symbol))
    matched <- !is.na(idx)

    q$context_tissue_class[matched]    <- ta_sub$tissue_class[idx[matched]]
    q$context_tissue_detected[matched] <- ta_sub$detect_fraction[idx[matched]] > 0

    # Mismatch criteria:
    #   - low_expression:        no RNA, no protein in tissue
    #   - suppressed_in_tissue:  RNA there, no protein
    #   - variable_in_tissue with detect_fraction == 0:
    #     trace RNA in some cell types, no protein anywhere
    hard_fail <- !is.na(q$context_tissue_class) &
                  q$context_tissue_class %in%
                    c("low_expression", "suppressed_in_tissue")
    soft_fail <- !is.na(q$context_tissue_class) &
                  q$context_tissue_class == "variable_in_tissue" &
                  !is.na(q$context_tissue_detected) &
                  !q$context_tissue_detected
    tissue_mismatch <- hard_fail | soft_fail
  }

  # --- Apply target verdicts ---
  verdicts <- mapply(.triage_target,
                     gc                    = q$gene_class,
                     lp                    = q$loc_plausible,
                     rescued               = q$compartment_rescue,
                     atlas_blind           = q$atlas_blind_spot,
                     broad                 = broad,
                     tissue_mismatch       = tissue_mismatch,
                     tissue_class_in_ctx   = q$context_tissue_class,
                     MoreArgs = list(context_tissue_name = resolved_tissue),
                     SIMPLIFY = FALSE)

  q$verdict <- vapply(verdicts, `[[`, character(1), "verdict")
  q$reason  <- vapply(verdicts, `[[`, character(1), "reason")

  # --- Output columns ---
  candidate_cols <- c("gene_symbol", "gene_class", "protein_confidence",
                      "detection_rate", "mechanism_tier",
                      "localisation", "loc_plausible",
                      "compartment_rescue", "atlas_blind_spot",
                      "tau", "specificity_class",
                      "context_tissue_class", "context_tissue_detected",
                      "biofluid_note", "verdict", "reason")
  out_cols <- intersect(candidate_cols, names(q))
  q[, out_cols, drop = FALSE]
}


#' Generate biofluid interpretation note (target context)
#' @keywords internal
#' @noRd
.biofluid_note <- function(atlas_df) {
  n <- nrow(atlas_df)
  if (n == 0L) return(character(0))

  vapply(seq_len(n), function(i) {
    row <- atlas_df[i, ]
    if (is.na(row$found) || !row$found) return(NA_character_)

    is_sec  <- isTRUE(row$is_secreted == 1)
    is_mem  <- isTRUE(row$is_membrane == 1)
    is_nuc  <- isTRUE(row$is_nuclear == 1)
    is_cyto <- isTRUE(row$is_cytoplasmic == 1)
    is_mito <- isTRUE(row$is_mitochondrial == 1)
    is_proj <- "is_projection" %in% names(row) &&
      isTRUE(row$is_projection == 1)
    conf    <- row$protein_confidence

    if (is_sec)  return("Expected in biofluids (secreted)")
    if (is_mem)  return("May appear via shedding/exosomes")
    if (is_proj) return("May appear via damage-release (axonal/dendritic injury)")
    if ((is_nuc || is_cyto) && !is.na(conf) && conf > 0.7)
      return("Biofluid elevation suggests tissue damage")
    if (is_mito && !is.na(conf) && conf > 0.7)
      return("Biofluid elevation suggests mitochondrial damage")
    "Context-dependent"
  }, character(1))
}


#' Triage logic for a single gene (target-only)
#'
#' Gates apply in order of severity. The most severe failed gate
#' determines the verdict; the reason string accumulates all triggered
#' gates so the user sees the full picture.
#'
#' @keywords internal
#' @noRd
.triage_target <- function(gc, lp, rescued, atlas_blind, broad,
                           tissue_mismatch, tissue_class_in_ctx,
                           context_tissue_name) {

  # Helper: turn a tissue_class string into a human-readable phrase.
  tissue_failure_text <- function() {
    where <- if (!is.na(context_tissue_name))
      paste0(" in ", context_tissue_name) else ""
    switch(as.character(tissue_class_in_ctx),
      "low_expression"        = paste0("Gene not expressed", where,
                                       " (no RNA, no protein detected)."),
      "suppressed_in_tissue"  = paste0("Gene transcribed but not translated",
                                       where, " (RNA present, no protein",
                                       " detected across cell types)."),
      "variable_in_tissue"    = paste0("Trace RNA in some cell types but",
                                       " no protein detected", where, "."),
      paste0("Gene not detected at protein level", where, ".")
    )
  }

  # Not in atlas
  if (is.na(gc))
    return(list(verdict = NA_character_, reason = "Gene not in atlas."))

  # --- Gate 1: compartment plausibility (most severe) ---
  if (!is.na(lp) && !lp) {
    r <- "Subcellular localisation inconsistent with sample compartment."
    extras <- character(0)
    if (gc == "consistently_suppressed")
      extras <- c(extras, "Also suppressed class; RNA evidence unreliable.")
    if (broad)
      extras <- c(extras, "Also broadly expressed; tissue-of-origin unsupported.")
    if (tissue_mismatch)
      extras <- c(extras, paste("Also:", tolower(tissue_failure_text())))
    return(list(
      verdict = "compartment_implausible",
      reason  = paste(c(r, extras), collapse = " ")
    ))
  }

  # --- Gate 2: tissue-of-origin mismatch ---
  if (tissue_mismatch) {
    r <- tissue_failure_text()
    extras <- character(0)
    if (gc == "consistently_suppressed")
      extras <- c(extras, "Also suppressed class; RNA evidence unreliable.")
    if (broad)
      extras <- c(extras, "Also broadly expressed.")
    return(list(
      verdict = "tissue_origin_mismatch",
      reason  = paste(c(r, extras), collapse = " ")
    ))
  }

  # If the compartment gate was only passed via the projection-protein
  # rescue, prepend a transparency note.
  rescue_note <- if (isTRUE(rescued))
    paste("Compartment gate rescued: projection-localised protein in",
          "biofluid (atlas microscopy under-detects axonal/dendritic",
          "compartments). Verdict reflects concordance and tissue gates only.")
  else
    NULL

  # --- Gate 3: concordance class ---
  if (gc == "consistently_suppressed") {
    if (isTRUE(atlas_blind)) {
      # Atlas is blind to this gene's compartment in biofluid context —
      # the suppressed call is an artifact of measurement, not biology.
      return(list(
        verdict = "RNA-supported",
        reason  = paste(c(rescue_note,
                          "Suppressed concordance call is an atlas artifact:",
                          "cellular proteomics under-detects secreted and",
                          "projection proteins that have left the cell of",
                          "origin. RNA-protein concordance cannot be",
                          "evaluated from this atlas for this gene in",
                          "biofluid context."),
                        collapse = " ")
      ))
    }
    return(list(
      verdict = "RNA-unreliable",
      reason  = paste(c(rescue_note,
                        "Suppressed. RNA does not predict protein.",
                        "RNA-based target evidence is not defensible;",
                        "protein-level evidence required."),
                      collapse = " ")
    ))
  }

  if (gc == "variable") {
    r <- "Variable concordance; RNA-protein relationship is tissue-dependent."
    if (broad)
      r <- paste(r, "Additionally broadly expressed; tissue-of-origin unsupported.")
    return(list(
      verdict = "caution",
      reason  = paste(c(rescue_note, r), collapse = " ")
    ))
  }

  if (gc == "low_expression")
    return(list(
      verdict = "caution",
      reason  = paste(c(rescue_note,
                        "Low expression class — sparse atlas data;",
                        "insufficient evidence."),
                      collapse = " ")
    ))

  # --- Gate 4: tissue specificity (only matters when concordant) ---
  if (gc == "consistently_concordant") {
    if (broad)
      return(list(
        verdict = "RNA-supported_but_broad",
        reason  = paste(c(rescue_note,
                          "Concordant and compartment-plausible, but",
                          "broadly expressed across tissues. Tissue-of-origin",
                          "claims (e.g. cell-type-specific drug targeting)",
                          "are unsupported. Use for non-tissue-specific",
                          "target claims only."),
                        collapse = " ")
      ))
    return(list(
      verdict = "RNA-supported",
      reason  = paste(c(rescue_note,
                        "Concordant and compartment-plausible.",
                        "RNA-based target evidence is defensible."),
                      collapse = " ")
    ))
  }

  # Fallback — unrecognised gene_class value
  list(verdict = "caution",
       reason  = paste(c(rescue_note, "Unrecognised gene class."),
                       collapse = " "))
}
