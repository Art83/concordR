# =============================================================================
# triage: drug target defensibility
# =============================================================================

#' Triage candidate genes for drug target defensibility
#'
#' Evaluates whether a gene's nomination as a drug target is defensible
#' given the concordance atlas. Drug-target claims demand a coherent
#' mechanistic chain: RNA must reliably predict protein, the protein
#' must be plausibly detected in the sample compartment that nominated
#' it, and the gene must be expressed at protein level in a tissue
#' relevant to the claim. Failure of any link breaks the chain.
#'
#' @section Automatic tissue inference:
#' For substrate-specific biofluids, the package automatically infers
#' which tissues are relevant: CSF implies the four brain regions
#' (Caudate, Cerebellum, Cerebral cortex, Hippocampus), urine implies
#' Kidney. Plasma and serum are systemic, no tissue is implied. When
#' \code{context_tissue} is supplied, it overrides the automatic
#' inference and narrows to the specified region.
#'
#' @section Gates (in order of severity):
#' \enumerate{
#'   \item \strong{Compartment plausibility}: subcellular localisation
#'     inconsistent with sample compartment. A damage-release rescue
#'     applies in biofluid contexts for projection-localised proteins,
#'     and a secretion rescue applies when the protein is secreted FROM
#'     a tissue relevant to the claim (not merely secreted anywhere).
#'   \item \strong{Tissue-of-origin coherence}: gene must be detected
#'     at protein level (IHC, MS, or PaxDB) in at least one tissue
#'     implied by the substrate. Fires automatically for CSF and urine;
#'     always fires when \code{context_tissue} is supplied.
#'   \item \strong{Concordance class}: suppressed class-RNA does not
#'     predict protein. Cross-platform evidence (HPA vs GTEx) is
#'     evaluated using data from the relevant tissues only.
#'   \item \strong{Tissue specificity}: concordant and plausible, but
#'     broadly expressed-tissue-of-origin claims unsupported.
#' }
#'
#' @section Verdicts:
#' \describe{
#'   \item{\code{RNA-supported}}{All gates pass.}
#'   \item{\code{RNA-supported_but_broad}}{Concordant and plausible, but
#'     broadly expressed.}
#'   \item{\code{RNA-unreliable}}{Suppressed class, confirmed
#'     cross-platform. RNA evidence not defensible.}
#'   \item{\code{tissue_origin_mismatch}}{Gene not expressed at protein
#'     level in tissue(s) relevant to the claim.}
#'   \item{\code{compartment_implausible}}{Subcellular localisation
#'     inconsistent with sample compartment.}
#'   \item{\code{caution}}{Variable class, low expression, or suppressed
#'     but cross-platform disagreement.}
#'   \item{\code{NA}}{Gene not in atlas.}
#' }
#'
#' @param genes Character vector of HGNC gene symbols.
#' @param sample_type Character. The biological compartment where proteins
#'   were measured. See \code{\link{flag_genes}} for the full list.
#' @param context_tissue Character or \code{NULL}. Overrides automatic
#'   tissue inference. When supplied, narrows the tissue gate to this
#'   specific tissue. Brain is resolved as four regions, pass one
#'   explicitly (e.g. \code{"Cerebellum"}), not \code{"Brain"}.
#'
#' @return A data frame with atlas annotations plus verdict columns.
#'
#' @examples
#' \dontrun{
#' # CSF: automatically checks brain tissues
#' triage(c("NEFL", "GFAP", "SPANXA1"), sample_type = "csf")
#'
#' # CSF narrowed to one region
#' triage(c("NEFL", "GFAP", "SPANXA1"),
#'        sample_type = "csf", context_tissue = "Cerebellum")
#'
#' # Plasma: systemic, no automatic tissue inference
#' triage(c("AR", "ESR1", "CD274", "TP53"), sample_type = "plasma")
#' }
#'
#' @seealso \code{\link{flag_genes}}, \code{\link{audit_enrichment}}
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
  
  # --- Resolve effective tissues ---
  # context_tissue overrides; otherwise infer from sample_type
  implied_tissues <- .SAMPLE_TISSUE_MAP[[sample_type]]
  is_biofluid     <- sample_type %in% c("plasma", "serum", "csf")
  is_systemic     <- is.null(implied_tissues) && is.null(context_tissue)
  
  if (!is.null(context_tissue)) {
    ta <- .ensure_tissue_atlas()
    resolved_tissue  <- .validate_context_tissue(context_tissue, ta)
    effective_tissues <- resolved_tissue
  } else if (!is.null(implied_tissues)) {
    ta <- .ensure_tissue_atlas()
    effective_tissues <- implied_tissues
  } else {
    ta <- NULL
    effective_tissues <- NULL
  }
  
  # --- Compute tissue-level evidence for effective tissues ---
  q$tissue_detected_any        <- NA
  q$tissue_best_class          <- NA_character_
  q$tissue_protein_evidence    <- NA_integer_
  q$tissue_cross_platform      <- NA_real_
  q$tissue_ms_detected         <- NA
  q$primary_tissue             <- NA_character_
  q$n_tissues_detected_relevant <- NA_integer_
  tissue_mismatch              <- rep(FALSE, nrow(q))
  
  if (!is.null(effective_tissues) && !is.null(ta)) {
    tissue_info <- .compute_tissue_evidence(
      q$gene_symbol, ta, effective_tissues
    )
    # merge by gene_symbol
    idx <- match(toupper(q$gene_symbol), toupper(tissue_info$gene_symbol))
    matched <- !is.na(idx)
    
    q$tissue_detected_any[matched]         <- tissue_info$detected_any[idx[matched]]
    q$tissue_best_class[matched]           <- tissue_info$best_class[idx[matched]]
    q$tissue_protein_evidence[matched]     <- tissue_info$max_protein_evidence[idx[matched]]
    q$tissue_cross_platform[matched]       <- tissue_info$cross_platform_agree[idx[matched]]
    q$tissue_ms_detected[matched]          <- tissue_info$ms_detected_any[idx[matched]]
    q$n_tissues_detected_relevant[matched] <- tissue_info$n_detected[idx[matched]]
    
    # Tissue mismatch: gene is in tissue atlas for these tissues but
    # not detected by ANY platform (IHC, MS, PaxDB)
    not_detected <- matched & !is.na(q$tissue_detected_any) & !q$tissue_detected_any
    absent       <- q$found & is.na(q$tissue_detected_any)
    tissue_mismatch <- not_detected | absent
  }
  
  # --- Primary tissue (where is this gene most expressed?) ---
  # Computed from full tissue atlas, not just effective tissues
  if (!is.null(ta) && length(found_idx) > 0L) {
    primary <- .compute_primary_tissue(q$gene_symbol[found_idx], ta)
    pidx <- match(toupper(q$gene_symbol), toupper(primary$gene_symbol))
    pmatched <- !is.na(pidx)
    q$primary_tissue[pmatched] <- primary$primary_tissue[pidx[pmatched]]
  }
  
  # --- Atlas blind-spot rescue (tissue-aware) ---
  q$compartment_rescue <- FALSE
  q$atlas_blind_spot   <- FALSE
  if (is_biofluid && length(found_idx) > 0L) {
    is_proj <- "is_projection" %in% names(q) &
      !is.na(q$is_projection) & q$is_projection == 1
    is_sec  <- "is_secreted"   %in% names(q) &
      !is.na(q$is_secreted)   & q$is_secreted   == 1
    
    if (is_systemic) {
      blind_mask <- (is_proj | is_sec)
    } else {
      # Projection proteins: IHC under-detects, so any tissue-level
      # evidence (IHC, MS, PaxDB) in relevant tissues justifies rescue
      proj_expressed <- is_proj & !is.na(q$tissue_detected_any) & q$tissue_detected_any
      
      # Secreted proteins: bulk MS/PaxDB can't distinguish real tissue
      # expression from blood contamination. Only IHC cell-type detection
      # (detect_fraction > 0) confirms the protein is made in this tissue.
      sec_expressed <- is_sec & .check_ihc_in_tissues(
        q$gene_symbol, ta, effective_tissues
      )
      
      blind_mask <- proj_expressed | sec_expressed
    }
    blind_mask[is.na(blind_mask)] <- FALSE
    q$atlas_blind_spot[blind_mask] <- TRUE
    
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
  
  # --- Apply target verdicts ---
  verdicts <- mapply(.triage_target,
                     gc                    = q$gene_class,
                     lp                    = q$loc_plausible,
                     rescued               = q$compartment_rescue,
                     atlas_blind           = q$atlas_blind_spot,
                     broad                 = broad,
                     tissue_mismatch       = tissue_mismatch,
                     tissue_best_class     = q$tissue_best_class,
                     cross_platform_agree  = q$tissue_cross_platform,
                     protein_evidence      = q$tissue_protein_evidence,
                     primary_tissue        = q$primary_tissue,
                     specificity_class     = if ("specificity_class" %in% names(q))
                       q$specificity_class
                     else rep(NA_character_, nrow(q)),
                     MoreArgs = list(
                       effective_tissue_names = paste(effective_tissues, collapse = ", "),
                       is_systemic            = is_systemic
                     ),
                     SIMPLIFY = FALSE)
  
  q$verdict <- vapply(verdicts, `[[`, character(1), "verdict")
  q$reason  <- vapply(verdicts, `[[`, character(1), "reason")
  
  # --- Output columns ---
  candidate_cols <- c("gene_symbol", "gene_class", "protein_confidence",
                      "detection_rate", "mechanism_tier",
                      "localisation", "loc_plausible",
                      "compartment_rescue", "atlas_blind_spot",
                      "tau", "specificity_class", "primary_tissue",
                      "tissue_best_class", "tissue_detected_any",
                      "tissue_protein_evidence", "tissue_cross_platform",
                      "n_tissues_detected_relevant",
                      "biofluid_note", "verdict", "reason")
  out_cols <- intersect(candidate_cols, names(q))
  q[, out_cols, drop = FALSE]
}


# =============================================================================
# Tissue evidence helpers
# =============================================================================

#' Compute per-gene tissue evidence for a set of target tissues
#'
#' Multi-evidence detection: a gene is "detected" in a tissue if ANY
#' platform sees it (IHC detect_fraction > 0, GTEx MS, or PaxDB).
#'
#' @param gene_symbols Character vector.
#' @param tissue_atlas Tissue-resolved atlas data frame.
#' @param tissues Character vector of tissue names to check.
#' @return Data frame with one row per gene.
#' @keywords internal
#' @noRd
.compute_tissue_evidence <- function(gene_symbols, tissue_atlas, tissues) {
  
  ta_sub <- tissue_atlas[tissue_atlas$ts_tissue %in% tissues &
                           toupper(tissue_atlas$gene_symbol) %in% toupper(gene_symbols), ,
                         drop = FALSE]
  
  if (nrow(ta_sub) == 0L) {
    return(data.frame(
      gene_symbol          = gene_symbols,
      detected_any         = NA,
      best_class           = NA_character_,
      max_protein_evidence = NA_integer_,
      cross_platform_agree = NA_real_,
      ms_detected_any      = NA,
      n_detected           = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }
  
  # Multi-evidence detection per row
  ihc_det <- !is.na(ta_sub$detect_fraction) & ta_sub$detect_fraction > 0
  ms_det  <- if ("ms_detected" %in% names(ta_sub))
    !is.na(ta_sub$ms_detected) & ta_sub$ms_detected
  else rep(FALSE, nrow(ta_sub))
  pax_det <- if ("paxdb_detected" %in% names(ta_sub))
    !is.na(ta_sub$paxdb_detected) & ta_sub$paxdb_detected
  else rep(FALSE, nrow(ta_sub))
  ta_sub$any_detected <- ihc_det | ms_det | pax_det
  ta_sub$ms_det       <- ms_det
  
  if (!"protein_evidence_sources" %in% names(ta_sub)) {
    ta_sub$protein_evidence_sources <- as.integer(ihc_det) +
      as.integer(ms_det) +
      as.integer(pax_det)
  }
  
  # Class severity
  class_rank <- c("concordant_in_tissue" = 1L, "variable_in_tissue" = 2L,
                  "suppressed_in_tissue" = 3L, "low_expression" = 4L)
  ta_sub$class_rank <- class_rank[ta_sub$tissue_class]
  
  # Aggregate per gene
  key <- toupper(ta_sub$gene_symbol)
  agg <- data.frame(
    gene_symbol_upper    = tapply(key, key, `[`, 1),
    detected_any         = as.logical(tapply(ta_sub$any_detected, key, any, na.rm = TRUE)),
    max_protein_evidence = as.integer(tapply(ta_sub$protein_evidence_sources, key, max, na.rm = TRUE)),
    ms_detected_any      = as.logical(tapply(ta_sub$ms_det, key, any, na.rm = TRUE)),
    n_detected           = as.integer(tapply(ta_sub$any_detected, key, sum, na.rm = TRUE)),
    best_rank            = tapply(ta_sub$class_rank, key, min, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Cross-platform agreement
  if ("platforms_agree" %in% names(ta_sub)) {
    agg$cross_platform_agree <- tapply(ta_sub$platforms_agree, key, mean, na.rm = TRUE)
  } else {
    agg$cross_platform_agree <- NA_real_
  }
  
  # Map rank back to class name
  rank_to_class <- setNames(names(class_rank), as.character(class_rank))
  agg$best_class <- rank_to_class[as.character(agg$best_rank)]
  agg$best_class[!is.finite(agg$best_rank)] <- NA_character_
  
  # Fix -Inf from max on empty
  agg$max_protein_evidence[!is.finite(agg$max_protein_evidence)] <- NA_integer_
  
  # Map back to original gene symbols
  gene_key <- toupper(gene_symbols)
  idx <- match(gene_key, agg$gene_symbol_upper)
  
  data.frame(
    gene_symbol          = gene_symbols,
    detected_any         = agg$detected_any[idx],
    best_class           = agg$best_class[idx],
    max_protein_evidence = agg$max_protein_evidence[idx],
    cross_platform_agree = agg$cross_platform_agree[idx],
    ms_detected_any      = agg$ms_detected_any[idx],
    n_detected           = agg$n_detected[idx],
    stringsAsFactors = FALSE
  )
}


#' Find the primary tissue of expression for each gene
#'
#' Returns the tissue with highest rna_rank_max for each gene.
#' Used to report WHERE a tissue-specific gene is expressed.
#'
#' @keywords internal
#' @noRd
.compute_primary_tissue <- function(gene_symbols, tissue_atlas) {
  ta_sub <- tissue_atlas[toupper(tissue_atlas$gene_symbol) %in% toupper(gene_symbols), ,
                         drop = FALSE]
  
  if (nrow(ta_sub) == 0L) {
    return(data.frame(gene_symbol = gene_symbols,
                      primary_tissue = NA_character_,
                      stringsAsFactors = FALSE))
  }
  
  ta_sub <- ta_sub[order(-ta_sub$protein_max, -ta_sub$detect_fraction,
                         -ta_sub$rna_rank_max), ]
  ta_sub <- ta_sub[!duplicated(toupper(ta_sub$gene_symbol)), ]
  
  idx <- match(toupper(gene_symbols), toupper(ta_sub$gene_symbol))
  data.frame(
    gene_symbol    = gene_symbols,
    primary_tissue = ta_sub$ts_tissue[idx],
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# Biofluid interpretation note
# =============================================================================

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


# =============================================================================
# Per-gene verdict logic
# =============================================================================

#' Triage logic for a single gene (target-only, tissue-aware)
#'
#' @keywords internal
#' @noRd
.triage_target <- function(gc, lp, rescued, atlas_blind, broad,
                           tissue_mismatch, tissue_best_class,
                           cross_platform_agree, protein_evidence,
                           primary_tissue, specificity_class,
                           effective_tissue_names, is_systemic) {
  
  # Helper: describe why the tissue gate failed
  tissue_failure_text <- function() {
    where <- if (nzchar(effective_tissue_names))
      paste0(" in ", effective_tissue_names) else ""
    
    if (is.na(tissue_best_class)) {
      return(paste0("Gene absent from tissue atlas", where, "."))
    }
    
    # Add primary tissue info when available and informative
    primary_note <- ""
    if (!is.na(primary_tissue) && !is.na(specificity_class) &&
        specificity_class %in% c("tissue_specific", "tissue_enriched")) {
      primary_note <- paste0(" Primary expression tissue: ", primary_tissue, ".")
    }
    
    msg <- switch(as.character(tissue_best_class),
                  "low_expression"        = paste0("Gene not expressed", where,
                                                   " (no RNA, no protein detected)."),
                  "suppressed_in_tissue"  = paste0("Gene transcribed but not translated",
                                                   where, " (RNA present, no protein",
                                                   " detected)."),
                  "variable_in_tissue"    = paste0("Trace RNA in some cell types but",
                                                   " no protein detected", where, "."),
                  paste0("Gene not detected at protein level", where, ".")
    )
    paste0(msg, primary_note)
  }
  
  # Not in atlas
  if (is.na(gc))
    return(list(verdict = NA_character_, reason = "Gene not in atlas."))
  
  # --- Gate 1: compartment plausibility ---
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
  
  # --- Gate 2: tissue-of-origin ---
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
  
  # Transparency note for rescued compartment gate
  rescue_note <- if (isTRUE(rescued))
    paste("Compartment gate rescued: projection-localised protein in",
          "biofluid (atlas microscopy under-detects axonal/dendritic",
          "compartments). Verdict reflects concordance and tissue gates only.")
  else
    NULL
  
  # --- Gate 3: concordance class ---
  if (gc == "consistently_suppressed") {
    if (isTRUE(atlas_blind)) {
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
    
    # Cross-platform check (tissue-filtered): GTEx disagrees AND
    # multiple platforms detect protein in the relevant tissue
    if (!is.na(cross_platform_agree) && cross_platform_agree < 0.3 &&
        !is.na(protein_evidence) && protein_evidence >= 2) {
      return(list(
        verdict = "caution",
        reason  = paste(c(rescue_note,
                          "Classified as suppressed by HPA IHC, but GTEx MS",
                          "does not confirm suppression in relevant tissue(s)",
                          "and protein is detected by multiple platforms.",
                          "IHC may undercount this protein;",
                          "interpret RNA-protein discordance with caution."),
                        collapse = " ")
      ))
    }
    
    return(list(
      verdict = "RNA-unreliable",
      reason  = paste(c(rescue_note,
                        "Suppressed. RNA does not predict protein.",
                        "RNA-based target evidence is not defensible;",
                        "protein-level evidence required.",
                        if (!is.na(cross_platform_agree) && cross_platform_agree > 0.7)
                          "Cross-platform agreement confirms suppression."
                        else NULL),
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
  
  if (gc == "low_expression") {
    # Check if multi-platform evidence contradicts the HPA low-expression
    # call in relevant tissues
    if (!is.na(protein_evidence) && protein_evidence >= 2) {
      return(list(
        verdict = "caution",
        reason  = paste(c(rescue_note,
                          "Low expression in HPA cell-type atlas, but protein",
                          "detected by multiple platforms in relevant tissue(s).",
                          "HPA classification may reflect coverage gap."),
                        collapse = " ")
      ))
    }
    return(list(
      verdict = "caution",
      reason  = paste(c(rescue_note,
                        "Low expression class-sparse atlas data;",
                        "insufficient evidence."),
                      collapse = " ")
    ))
  }
  
  # --- Gate 4: tissue specificity ---
  if (gc == "consistently_concordant") {
    if (broad) {
      broad_reason <- if (!is_systemic)
        paste("Concordant and compartment-plausible, but broadly expressed",
              "across tissues. CSF/biofluid detection likely reflects plasma",
              "filtration rather than tissue-specific origin.")
      else
        paste("Concordant and compartment-plausible, but broadly expressed",
              "across tissues. Tissue-of-origin claims unsupported.")
      return(list(
        verdict = "RNA-supported_but_broad",
        reason  = paste(c(rescue_note, broad_reason), collapse = " ")
      ))
    }
    return(list(
      verdict = "RNA-supported",
      reason  = paste(c(rescue_note,
                        "Concordant and compartment-plausible.",
                        "RNA-based target evidence is defensible."),
                      collapse = " ")
    ))
  }
  
  # Fallback
  list(verdict = "caution",
       reason  = paste(c(rescue_note, "Unrecognised gene class."),
                       collapse = " "))
}


#' Check whether genes have IHC protein detection in any target tissue
#'
#' For secreted proteins, bulk tissue MS/PaxDB cannot distinguish
#' genuine tissue expression from blood contamination. Only cell-type
#' resolved IHC (detect_fraction > 0) confirms the protein is
#' produced locally.
#'
#' @keywords internal
#' @noRd
.check_ihc_in_tissues <- function(gene_symbols, tissue_atlas, tissues) {
  ta_sub <- tissue_atlas[tissue_atlas$ts_tissue %in% tissues &
                           toupper(tissue_atlas$gene_symbol) %in% toupper(gene_symbols), ,
                         drop = FALSE]
  
  if (nrow(ta_sub) == 0L) return(rep(FALSE, length(gene_symbols)))
  
  key <- toupper(ta_sub$gene_symbol)
  any_ihc <- tapply(ta_sub$detect_fraction > 0, key, any, na.rm = TRUE)
  
  idx <- match(toupper(gene_symbols), names(any_ihc))
  result <- any_ihc[idx]
  result[is.na(result)] <- FALSE
  as.logical(result)
}