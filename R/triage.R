# =============================================================================
# triage: independent gate audit for biofluid drug-target nominations
# =============================================================================
# Gates, decoupled by data source, each scored on its own footprint:
#   biofluid_plausibility  localisation_flags (UniProt), proteome-wide
#   tissue                 hpa_tissue (HPA + GTEx + PaxDb union), wide
#   specificity            HPA protein breadth, wide
#   concordance            concordance_atlas (the 11k atlas), narrow
# No rollup verdict is returned; the columns are independent and the caller
# combines them. tissue_atlas supplies only the cell-type IHC used by the
# biofluid rescue and the suppressed cross-platform check.

#' Audit candidate genes for biofluid drug-target defensibility
#'
#' Reports, per gene, what each gate independently found. There is no single
#' verdict: each gate is its own controlled-vocabulary column, scored on its
#' own data source and NA only where that source is silent. The caller filters
#' and combines the columns, widest gate first, then narrowing into the atlas.
#'
#' @section Automatic tissue inference:
#' Plasma, serum, and CSF are systemic and receive no automatic tissue gate.
#' Although CSF is anatomically proximal to brain, roughly 80% of its protein
#' derives from plasma, so automatic brain attribution is unsound. Urine
#' implies Kidney. A supplied \code{context_tissue} overrides inference and
#' narrows to that region. Brain must be given as one of four regions.
#'
#' @section Gate columns:
#' \describe{
#'   \item{\code{concordance}}{Atlas only. concordant, variable, suppressed,
#'     low_expression, not_in_atlas.}
#'   \item{\code{biofluid_plausibility}}{UniProt localisation, raw call.
#'     plausible, implausible, no_data.}
#'   \item{\code{biofluid_rescue}}{Logical. TRUE where a projection or secreted
#'     protein would be rescued in a biofluid (microscopy under-detects
#'     proteins that have left the cell of origin). Reported beside
#'     plausibility, not folded into it.}
#'   \item{\code{tissue}}{Union gate. detected, detected_1platform, discordant,
#'     mismatch, coverage_limited, not_in_panel, systemic_na.}
#'   \item{\code{specificity}}{HPA breadth. specific, broad, na.}
#'   \item{\code{protein_tissues}}{Tissues of detection, from the union gate,
#'     populated off-atlas.}
#'   \item{\code{rna_tissues}}{Atlas-bound, NA off-atlas.}
#'   \item{\code{reason}}{Prose summary for reading, not for filtering.}
#' }
#'
#' @param genes Character vector of HGNC gene symbols.
#' @param sample_type Character. The biofluid or fraction sampled.
#' @param context_tissue Character or \code{NULL}. Overrides automatic tissue
#'   inference. Brain must be one region (e.g. \code{"Cerebellum"}).
#'
#' @return A data frame, one row per gene, with the independent gate columns
#'   above plus a prose \code{reason}. No rollup verdict.
#' @seealso \code{\link{query_atlas}}
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

  # --- Concordance lookup (gene_class etc.) ---
  q <- query_atlas(genes, missing = "keep")

  # --- Localisation: single UniProt source for every gene ---
  q <- .attach_localisation(q, .ensure_localisation_flags())
  q$loc_plausible <- .check_plausibility(q, sample_type)
  q$localisation  <- .location_labels(q)

  # --- Resolve effective tissues (HPA namespace) ---
  hpa <- .ensure_hpa_tissue()
  ta  <- .ensure_tissue_atlas()
  is_biofluid <- sample_type %in% c("plasma", "serum", "csf")

  implied_tissues <- .SAMPLE_TISSUE_MAP[[sample_type]]
  if (!is.null(context_tissue)) {
    effective_tissues <- .validate_context_tissue(context_tissue,
                                                  unique(hpa$tissue))
  } else {
    effective_tissues <- implied_tissues
  }
  is_systemic       <- is.null(effective_tissues)
  tissue_gate_active <- !is.null(effective_tissues)

  # Atlas-namespace equivalents for the tissue_atlas-side checks
  atlas_effective <- if (tissue_gate_active)
    stats::na.omit(.hpa_to_atlas_tissue(effective_tissues, ta)) else character(0)

  # --- Tissue gate from hpa_tissue verdict ---
  q$tissue_status           <- NA_character_
  q$tissue_detected_any     <- FALSE
  q$n_relevant_insufficient <- NA_integer_
  q$insufficient_tissues    <- NA_character_
  if (tissue_gate_active) {
    ts  <- .tissue_status_from_hpa(q$gene_symbol, hpa, effective_tissues)
    idx <- match(toupper(q$gene_symbol), toupper(ts$gene_symbol))
    q$tissue_status           <- ts$tissue_status[idx]
    q$tissue_detected_any     <- ts$tissue_detected_any[idx]
    q$n_relevant_insufficient <- ts$n_relevant_insufficient[idx]
    q$insufficient_tissues    <- ts$insufficient_tissues[idx]
  }

  # --- Primary tissue and protein breadth (full HPA panel, all genes) ---
  primary <- .compute_primary_tissue(q$gene_symbol, hpa)
  q$primary_tissue <- primary$primary_tissue[
    match(toupper(q$gene_symbol), toupper(primary$gene_symbol))]

  pb <- .protein_tissue_breadth_hpa(q$gene_symbol, hpa)
  q$protein_breadth <- pb$protein_breadth[
    match(toupper(q$gene_symbol), toupper(pb$gene_symbol))]
  q$broad <- !is.na(q$protein_breadth) & q$protein_breadth >= .BROAD_BREADTH

  # --- Biofluid rescue (projection/secreted blind spot) ---
  # Surfaced as its own boolean. Raw localisation plausibility is left intact
  # so biofluid_plausibility reports the localisation call alone.
  q$biofluid_rescue  <- FALSE
  q$atlas_blind_spot <- FALSE
  if (is_biofluid) {
    is_proj <- !is.na(q$is_projection) & q$is_projection == 1
    is_sec  <- !is.na(q$is_secreted)   & q$is_secreted   == 1
    if (is_systemic) {
      blind <- is_proj | is_sec
    } else {
      proj_expressed <- is_proj & q$tissue_detected_any
      sec_expressed  <- is_sec &
        .check_ihc_in_tissues(q$gene_symbol, ta, atlas_effective)
      blind <- proj_expressed | sec_expressed
    }
    blind[is.na(blind)] <- FALSE
    q$atlas_blind_spot[blind] <- TRUE
    # Rescue only flags genes whose raw localisation looked implausible.
    rescue <- blind & !is.na(q$loc_plausible) & !q$loc_plausible
    q$biofluid_rescue[rescue] <- TRUE
  }

  # --- Suppressed cross-platform evidence (tissue-filtered) ---
  xp <- .suppressed_xplatform(q$gene_symbol, ta, atlas_effective)
  q$cross_platform_agree <- xp$cross_platform_agree[
    match(toupper(q$gene_symbol), toupper(xp$gene_symbol))]
  q$protein_evidence <- xp$protein_evidence[
    match(toupper(q$gene_symbol), toupper(xp$gene_symbol))]

  # --- Biofluid interpretation note ---
  q$biofluid_note <- .biofluid_note(q)

  # --- Gate columns: each scored on its own data footprint, no rollup -------
  # Concordance (atlas only, not_in_atlas where the gene has no class).
  q$concordance <- c(
    consistently_concordant = "concordant",
    variable                = "variable",
    consistently_suppressed = "suppressed",
    low_expression          = "low_expression"
  )[q$gene_class]
  q$concordance[is.na(q$concordance)] <- "not_in_atlas"

  # Biofluid plausibility (UniProt localisation only, raw call). The rescue is
  # reported separately in biofluid_rescue, not folded in here.
  q$biofluid_plausibility <- ifelse(
    is.na(q$loc_plausible), "no_data",
    ifelse(q$loc_plausible, "plausible", "implausible"))

  # Tissue (union gate only). systemic_na means no tissue-of-origin gate ran.
  .TMAP <- c(detected = "detected", detected_weak = "detected_1platform",
             discordant = "discordant", absent = "mismatch",
             coverage_unknown = "coverage_limited", not_in_panel = "not_in_panel")
  q$tissue <- ifelse(is.na(q$tissue_status), "systemic_na", .TMAP[q$tissue_status])

  # Specificity (HPA breadth only).
  q$specificity <- ifelse(is.na(q$protein_breadth), "na",
                          ifelse(q$broad, "broad", "specific"))

  # protein_tissues from the union gate (populated off-atlas); rna_tissues
  # stays atlas-bound and NA off-atlas, since there is no RNA to list.
  q$protein_tissues <- .protein_tissues_from_gate(q$gene_symbol, hpa)

  # Suppressed-class rescue flag, retained for the prose reason only.
  supp_rescued <- q$concordance == "suppressed" &
    (q$atlas_blind_spot %in% TRUE |
       (!is.na(q$cross_platform_agree) & q$cross_platform_agree < 0.3 &
          !is.na(q$protein_evidence)   & q$protein_evidence >= 2))

  # --- Prose reason (for humans, not for filtering) -------------------------
  where <- if (nzchar(paste(effective_tissues, collapse = "")))
    paste0(" in ", paste(effective_tissues, collapse = ", ")) else ""
  q$reason <- .reason_from_axes(q$concordance, q$biofluid_plausibility,
                                q$biofluid_rescue, q$tissue, q$specificity,
                                supp_rescued, where)

  # --- Output: independent gate columns, no verdict -------------------------
  candidate_cols <- c("gene_symbol",
                      "concordance", "biofluid_plausibility", "biofluid_rescue",
                      "tissue", "specificity",
                      "primary_tissue", "protein_breadth", "tissue_detected_any",
                      "biofluid_note", "protein_confidence", "detection_rate",
                      "localisation", "protein_tissues", "rna_tissues", "reason")
  q[, intersect(candidate_cols, names(q)), drop = FALSE]
}


#' Protein breadth above which a gene is treated as broadly expressed.
#' @keywords internal
#' @noRd
.BROAD_BREADTH <- 0.65


# =============================================================================
# Data assembly helpers
# =============================================================================

#' Attach localisation flags from the single UniProt source onto every gene.
#' @keywords internal
#' @noRd
.attach_localisation <- function(q, loc) {
  for (col in .LOC_COLS) q[[col]] <- NA_integer_
  if (is.null(loc) || nrow(loc) == 0L) return(q)
  idx     <- match(toupper(q$gene_symbol), toupper(loc$gene_symbol))
  matched <- !is.na(idx)
  for (col in intersect(.LOC_COLS, names(loc)))
    q[[col]][matched] <- as.integer(loc[[col]][idx[matched]])
  q
}


#' Gene-level tissue status over the relevant tissues, from hpa_tissue.
#'
#' Presence in any relevant tissue wins. Corroborated absence (verdict
#' \code{concordant_absent}) only registers when nothing was detected. Single
#' platform (\code{detected_by_one}, \code{insufficient_coverage}) and absent
#' panel coverage are tracked so the verdict can report them without rejecting.
#'
#' @return Data frame: gene_symbol, tissue_status, tissue_detected_any,
#'   n_relevant_insufficient, insufficient_tissues.
#' @keywords internal
#' @noRd
.tissue_status_from_hpa <- function(gene_symbols, hpa, tissues) {
  n   <- length(gene_symbols)
  out <- data.frame(
    gene_symbol             = gene_symbols,
    tissue_status           = NA_character_,
    tissue_detected_any     = FALSE,
    n_relevant_insufficient = 0L,
    insufficient_tissues    = NA_character_,
    stringsAsFactors        = FALSE
  )
  sub <- hpa[hpa$tissue %in% tissues &
               toupper(hpa$gene_symbol) %in% toupper(gene_symbols), ,
             drop = FALSE]
  key <- toupper(gene_symbols)

  for (i in seq_len(n)) {
    rows <- sub[toupper(sub$gene_symbol) == key[i], , drop = FALSE]
    rows <- rows[!is.na(rows$verdict), , drop = FALSE]   # ignore unscored cells
    if (nrow(rows) == 0L) {
      out$tissue_status[i] <- "not_in_panel"
      next
    }
    v      <- rows$verdict
    insuff <- rows$tissue[v == "insufficient_coverage"]
    out$n_relevant_insufficient[i] <- length(insuff)
    if (length(insuff))
      out$insufficient_tissues[i] <- paste(sort(unique(insuff)), collapse = ", ")

    out$tissue_status[i] <-
      if (any(v == "detected", na.rm = TRUE))               "detected"
      else if (any(v == "detected_by_one", na.rm = TRUE))    "detected_weak"
      else if (any(v == "discordant", na.rm = TRUE))         "discordant"
      else if (any(v == "concordant_absent", na.rm = TRUE))  "absent"
      else                                                   "coverage_unknown"

    out$tissue_detected_any[i] <-
      out$tissue_status[i] %in% c("detected", "detected_weak", "discordant")
  }
  out
}


#' Primary tissue of expression, from the union gate.
#'
#' Among the tissues where the gene is detected, narrows to the leaders by
#' n_detected, then ihc_max_level, then ihc_n_detected, using only keys that
#' have a value. Returns the single leading tissue; "multiple" when two or more
#' are tied on every available key, which is an honest tie at vote resolution
#' rather than a measured abundance tie; NA when the gene is detected nowhere.
#' @keywords internal
#' @noRd
.compute_primary_tissue <- function(gene_symbols, hpa) {
  out <- rep(NA_character_, length(gene_symbols))
  sub <- hpa[toupper(hpa$gene_symbol) %in% toupper(gene_symbols) &
               hpa$any_detected %in% TRUE, , drop = FALSE]
  if (nrow(sub) == 0L)
    return(data.frame(gene_symbol = gene_symbols, primary_tissue = out,
                      stringsAsFactors = FALSE))

  parts <- split(sub, toupper(sub$gene_symbol))
  for (i in seq_along(gene_symbols)) {
    rows <- parts[[toupper(gene_symbols[i])]]
    if (is.null(rows) || nrow(rows) == 0L) next     # detected nowhere -> NA

    cand <- rows
    for (k in c("n_detected", "ihc_max_level", "ihc_n_detected")) {
      v <- cand[[k]]
      if (all(is.na(v))) next                       # key has no value, skip it
      m <- max(v, na.rm = TRUE)
      cand <- cand[!is.na(v) & v == m, , drop = FALSE]
    }
    out[i] <- if (nrow(cand) == 1L) cand$tissue[1L] else "multiple"
  }
  data.frame(gene_symbol = gene_symbols, primary_tissue = out,
             stringsAsFactors = FALSE)
}


#' Protein tissue breadth, full HPA panel. Fraction of tissues with confident
#' detection (Medium/High IHC or two-platform agreement).
#' @keywords internal
#' @noRd
.protein_tissue_breadth_hpa <- function(gene_symbols, hpa) {
  n_total <- length(unique(hpa$tissue))
  sub <- hpa[toupper(hpa$gene_symbol) %in% toupper(gene_symbols), , drop = FALSE]
  if (nrow(sub) == 0L)
    return(data.frame(gene_symbol = gene_symbols,
                      protein_breadth = NA_real_, stringsAsFactors = FALSE))

  det_col <- if ("confident_detected" %in% names(sub))
    "confident_detected" else "any_detected"
  key   <- toupper(sub$gene_symbol)
  n_det <- tapply(sub[[det_col]], key, sum, na.rm = TRUE)
  idx   <- match(toupper(gene_symbols), names(n_det))
  data.frame(gene_symbol = gene_symbols,
             protein_breadth = as.numeric(n_det[idx]) / n_total,
             stringsAsFactors = FALSE)
}


#' Cell-type IHC detection in the relevant atlas tissues (secreted rescue).
#' @keywords internal
#' @noRd
.check_ihc_in_tissues <- function(gene_symbols, tissue_atlas, atlas_tissues) {
  if (length(atlas_tissues) == 0L) return(rep(FALSE, length(gene_symbols)))
  sub <- tissue_atlas[tissue_atlas$ts_tissue %in% atlas_tissues &
                        toupper(tissue_atlas$gene_symbol) %in% toupper(gene_symbols), ,
                      drop = FALSE]
  if (nrow(sub) == 0L) return(rep(FALSE, length(gene_symbols)))
  key     <- toupper(sub$gene_symbol)
  any_ihc <- tapply(!is.na(sub$detect_fraction) & sub$detect_fraction > 0,
                    key, any, na.rm = TRUE)
  res <- as.logical(any_ihc[match(toupper(gene_symbols), names(any_ihc))])
  res[is.na(res)] <- FALSE
  res
}


#' Suppressed cross-platform evidence in the relevant atlas tissues.
#' @keywords internal
#' @noRd
.suppressed_xplatform <- function(gene_symbols, tissue_atlas, atlas_tissues) {
  base <- data.frame(gene_symbol = gene_symbols,
                     cross_platform_agree = NA_real_,
                     protein_evidence = NA_integer_, stringsAsFactors = FALSE)
  if (length(atlas_tissues) == 0L) return(base)
  sub <- tissue_atlas[tissue_atlas$ts_tissue %in% atlas_tissues &
                        toupper(tissue_atlas$gene_symbol) %in% toupper(gene_symbols), ,
                      drop = FALSE]
  if (nrow(sub) == 0L) return(base)
  key   <- toupper(sub$gene_symbol)
  agree <- tapply(sub$platforms_agree, key,
                  function(x) { m <- mean(x, na.rm = TRUE)
                               if (is.nan(m)) NA_real_ else m })
  pev   <- tapply(sub$protein_evidence_sources, key,
                  function(x) { m <- max(x, na.rm = TRUE)
                               if (is.finite(m)) as.integer(m) else NA_integer_ })
  gkey <- toupper(gene_symbols)
  base$cross_platform_agree <- as.numeric(agree[match(gkey, names(agree))])
  base$protein_evidence     <- as.integer(pev[match(gkey, names(pev))])
  base
}


# =============================================================================
# Biofluid interpretation note
# =============================================================================

#' @keywords internal
#' @noRd
.biofluid_note <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(character(0))
  vapply(seq_len(n), function(i) {
    row <- df[i, ]
    has_loc <- any(!is.na(row[intersect(.LOC_COLS, names(row))]))
    if (!has_loc) return(NA_character_)

    is_sec  <- isTRUE(row$is_secreted == 1)
    is_mem  <- isTRUE(row$is_membrane == 1)
    is_nuc  <- isTRUE(row$is_nuclear == 1)
    is_cyto <- isTRUE(row$is_cytoplasmic == 1)
    is_mito <- isTRUE(row$is_mitochondrial == 1)
    is_proj <- isTRUE(row$is_projection == 1)
    conf    <- row$protein_confidence

    if (is_sec)  return("Expected in biofluids (secreted)")
    if (is_mem)  return("May appear via shedding or exosomes")
    if (is_proj) return("May appear via damage-release (axonal or dendritic injury)")
    if ((is_nuc || is_cyto) && !is.na(conf) && conf > 0.7)
      return("Biofluid elevation suggests tissue damage")
    if (is_mito && !is.na(conf) && conf > 0.7)
      return("Biofluid elevation suggests mitochondrial damage")
    "Context-dependent"
  }, character(1))
}


# =============================================================================
# protein_tissues from the union gate, and the prose reason
# =============================================================================

#' Tissues where the protein is detected, from the union gate. Populated for
#' any detected gene, including those outside the concordance atlas.
#' @keywords internal
#' @noRd
.protein_tissues_from_gate <- function(gene_symbols, hpa) {
  sub <- hpa[toupper(hpa$gene_symbol) %in% toupper(gene_symbols) &
               hpa$any_detected %in% TRUE, , drop = FALSE]
  if (nrow(sub) == 0L) return(rep(NA_character_, length(gene_symbols)))
  ag <- tapply(sub$tissue, toupper(sub$gene_symbol),
               function(x) paste(sort(unique(x)), collapse = ", "))
  unname(ag[match(toupper(gene_symbols), names(ag))])
}


#' Build a short, axis-consistent reason string. Biofluid plausibility is the
#' raw localisation call; the rescue is reported as its own clause.
#' @keywords internal
#' @noRd
.reason_from_axes <- function(concordance, plausibility, rescue, tissue,
                              specificity, supp_rescued, where) {
  conc <- c(
    concordant     = "RNA-protein concordant",
    variable       = "variable concordance",
    suppressed     = "suppressed, RNA does not predict protein",
    low_expression = "low expression",
    not_in_atlas   = "concordance not assessable, gene outside the atlas"
  )[concordance]
  plaus <- c(
    plausible   = "localisation plausible for the biofluid",
    implausible = "localisation inconsistent with the biofluid",
    no_data     = "no localisation data"
  )[plausibility]
  tis <- c(
    detected           = paste0("protein detected", where),
    detected_1platform = paste0("protein detected", where, " by one platform only"),
    discordant         = paste0("platforms disagree on protein presence", where),
    mismatch           = paste0("protein not detected", where, ", platforms agree on absence"),
    coverage_limited   = paste0("protein presence", where, " not assessable, thin coverage"),
    not_in_panel       = paste0("gene absent from the protein panel", where),
    systemic_na        = "systemic sample, no tissue-of-origin gate"
  )[tissue]
  spec <- c(
    specific = "tissue-specific at protein level",
    broad    = "broadly expressed at protein level",
    na       = "protein breadth unknown"
  )[specificity]

  base <- paste(conc, plaus, tis, spec, sep = "; ")
  base <- ifelse(rescue %in% TRUE,
                 paste0(base, "; rescued, projection or secreted protein in biofluid"),
                 base)
  ifelse(supp_rescued %in% TRUE,
         paste0(base, "; suppression not corroborated cross-platform."),
         paste0(base, "."))
}
