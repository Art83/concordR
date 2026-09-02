# =============================================================================
# triage: independent gate columns for drug target assessment
# =============================================================================

#' Triage candidate genes for drug target assessment
#'
#' Scores each gene on four independent gates and returns one column per
#' gate. No verdict is produced here. It's up to user to work with the gates and
#' work out a decision making strategy, if verdict is needed there's audit function.
#' Combine the gates with
#' \code{\link{audit_call}} for an accept/caution/decline call, or with
#' \code{\link{funnel_class}} for a first-failure partition in a stated
#' order. Each gate is a function of the atlas data and the sample type
#' only. No gate reads another gate's result.
#'
#' @section Gates:
#' \describe{
#'   \item{\code{concordance}}{Gene-level RNA to protein concordance class
#'     from the atlas. \code{"concordant"}, \code{"variable"},
#'     \code{"suppressed"}, \code{"suppressed_unconfirmed"},
#'     \code{"low_expression"}, or \code{"not_in_atlas"}.
#'     \code{"suppressed_unconfirmed"} means the atlas calls the gene
#'     suppressed but that call cannot be trusted, either because the
#'     protein is secreted or projection-localised and the sample is a
#'     biofluid (cellular proteomics under-detects proteins that have left
#'     the cell), or because the other platforms in the presence vote
#'     disagree with the IHC call in the context tissue.}
#'   \item{\code{tissue}}{Protein presence in the context tissue from the
#'     three-source vote (HPA IHC, GTEx MS, PaxDb). \code{"detected"} (two
#'     or more voters, all present), \code{"detected_1platform"} (one voter,
#'     present), \code{"discordant"} (voters disagree), \code{"mismatch"}
#'     (two or more voters, all absent), \code{"insufficient_coverage"}
#'     (one voter, absent), \code{"not_in_panel"} (no row), or
#'     \code{"systemic_na"} (no tissue implied and none supplied).}
#'   \item{\code{biofluid_plausibility}}{Whether the protein's subcellular
#'     localisation (UniProt) is consistent with the sample type.
#'     \code{"plausible"}, \code{"implausible"}, \code{"unknown"} (no
#'     localisation data), or \code{"na"} (sample type has no rule).}
#'   \item{\code{biofluid_rescue}}{Logical. \code{TRUE} when plausibility is
#'     \code{"implausible"} but the protein is projection-localised, the
#'     sample is plasma, serum, or CSF, and either no context tissue was
#'     given or the tissue gate finds the protein present there. This is
#'     the damage-release route (the NEFL pattern). \code{NA} when
#'     plausibility is \code{"na"}.}
#'   \item{\code{specificity}}{Protein-level tissue specificity from the HPA
#'     IHC panel. \code{"specific"} when protein tau (highest IHC level in
#'     any cell type, per tissue, across every tissue HPA stained) is at or
#'     above the shipped threshold, \code{"broad"} below
#'     it, \code{"sparse"} when the gene is IHC-detected in too few tissues
#'     for tau to describe a profile and the bulk platforms agree it is
#'     restricted, \code{"na"} when the gene has no reliable IHC anywhere
#'     or when its few IHC hits are contradicted by broad bulk detection
#'     (the antibody is uninformative for this gene). RNA is not used in
#'     this gate.}
#'   \item{\code{primary_tissue}}{For specific genes, the set of HPA
#'     tissues tied at the gene's highest IHC level, alphabetical and
#'     comma-joined. A single winner is undefined for most genes (median
#'     six tissues share the maximum), so the set is reported.
#'     \code{"multiple"} for broad genes. \code{NA} for sparse and na.}
#' }
#'
#' @section Tissue inference:
#' Plasma, serum, and CSF are systemic. Although CSF is anatomically
#' proximal to brain, roughly 80\% of its protein content derives from
#' plasma, so no brain region is implied automatically. Urine implies
#' Kidney. \code{context_tissue} overrides any implication and is validated
#' against the tissue vocabulary of the presence vote. Brain is resolved
#' into regions; pass one explicitly (\code{"Cerebral cortex"}), not
#' \code{"Brain"}.
#'
#' @param genes Character vector of HGNC gene symbols.
#' @param sample_type Character. Where the proteins were measured. One of
#'   \code{"plasma"}, \code{"serum"}, \code{"csf"}, \code{"urine"},
#'   \code{"tissue"}, \code{"cell_lysate"}, \code{"nuclear_fraction"},
#'   \code{"mitochondrial_fraction"}, \code{"membrane_fraction"}.
#' @param context_tissue Character or \code{NULL}. Tissue the claim is
#'   about. Overrides automatic inference.
#' @param details Logical. When \code{TRUE}, append the evidence behind each
#'   gate (\code{ensembl_id}, \code{gene_class}, \code{protein_confidence},
#'   \code{suppression_rate}, \code{tau}, \code{protein_breadth},
#'   \code{localisation}, \code{n_platforms_context},
#'   \code{mechanism_tier}). Default \code{FALSE}.
#'
#' @return A data frame with one row per input gene, in input order. Seven
#'   columns by default: \code{gene_symbol}, \code{concordance},
#'   \code{tissue}, \code{biofluid_plausibility}, \code{biofluid_rescue},
#'   \code{specificity}, \code{primary_tissue}. The resolved sample type and
#'   context tissue are stored in the \code{"triage_context"} attribute.
#'
#' @examples
#' \dontrun{
#' # CSF with no tissue implied
#' triage(c("NEFL", "GFAP", "SPANXA1"), sample_type = "csf")
#'
#' # CSF against a named brain region
#' triage(c("NEFL", "GFAP", "SPANXA1"), sample_type = "csf",
#'        context_tissue = "Cerebral cortex")
#'
#' # Evidence behind the gates
#' triage("NEFL", sample_type = "csf", details = TRUE)
#' }
#'
#' @seealso \code{\link{audit_call}}, \code{\link{funnel_class}},
#'   \code{\link{query_atlas}}
#' @export
triage <- function(genes, sample_type, context_tissue = NULL,
                   details = FALSE) {

  if (!is.character(genes) || length(genes) == 0L)
    stop("`genes` must be a non-empty character vector.", call. = FALSE)
  if (!is.logical(details) || length(details) != 1L || is.na(details))
    stop("`details` must be TRUE or FALSE.", call. = FALSE)
  sample_type <- .validate_sample_type(sample_type)

  thr <- .thresholds()
  hpa <- .ensure_hpa_tissue()
  if (is.null(hpa))
    stop("hpa_tissue not found in package data. sysdata.rda problem", call. = FALSE)
  ta <- .ensure_tissue_atlas()

  # --- Atlas annotations and localisation flags ---
  # query_atlas() deduplicates; re-expand so the output has one row per input
  # gene, in input order, including repeats and case variants.
  q <- query_atlas(genes, missing = "keep")
  q <- q[match(toupper(trimws(genes)), toupper(q$gene_symbol)), , drop = FALSE]
  q$gene_symbol <- ifelse(q$found %in% TRUE, q$gene_symbol, genes)
  rownames(q) <- NULL
  q <- .join_localisation(q)
  n <- nrow(q)

  # --- Context tissue (hpa_tissue namespace) ---
  target <- if (!is.null(context_tissue))
    .validate_context_tissue(context_tissue, unique(hpa$tissue))
  else
    .SAMPLE_TISSUE_MAP[[sample_type]]

  # --- Gate: tissue ---
  tis <- .gate_tissue(q, hpa, target)

  # --- Cross-platform agreement in the context tissue ---
  cross_agree <- .cross_platform_agree(q, ta, target)

  # --- Gate: concordance ---
  concordance <- .gate_concordance(q, sample_type, tis$n_platforms,
                                   cross_agree, thr)

  # --- Gate: biofluid plausibility and rescue ---
  plaus <- .gate_plausibility(q, sample_type)
  rescue <- .gate_rescue(q, plaus, tis$call, sample_type)

  # --- Gate: specificity and primary tissue ---
  spec <- .gate_specificity(q, thr)
  primary <- .primary_tissue(q, spec)

  # --- Assemble ---
  out <- data.frame(
    gene_symbol = q$gene_symbol,
    concordance = concordance,
    tissue = tis$call,
    biofluid_plausibility = plaus,
    biofluid_rescue = rescue,
    specificity = spec,
    primary_tissue = primary,
    stringsAsFactors = FALSE
  )

  if (details) {
    pick <- function(col) if (col %in% names(q)) q[[col]] else rep(NA, n)
    out$ensembl_id <- pick("ensembl_id")
    out$gene_class <- pick("gene_class")
    out$protein_confidence <- pick("protein_confidence")
    out$suppression_rate <- pick("suppression_rate")
    out$protein_tau <- pick("protein_tau")
    out$n_tissues_ihc <- pick("n_tissues_ihc")
    out$n_tissues_ihc_detected <- pick("n_tissues_ihc_detected")
    out$n_primary_tissues <- pick("n_primary_tissues")
    out$n_tissues_bulk_detected <- pick("n_tissues_bulk_detected")
    out$n_tissues_bulk_covered <- pick("n_tissues_bulk_covered")
    out$tau <- pick("tau")
    out$protein_breadth <- pick("protein_breadth")
    out$localisation <- .location_labels(q)
    out$n_platforms_context <- tis$n_platforms
    out$mechanism_tier <- pick("mechanism_tier")
  }

  rownames(out) <- NULL
  attr(out, "triage_context") <- list(sample_type    = sample_type,
                                      context_tissue = target)
  out
}


# ---internal function-helpers and data structures---

# =============================================================================
# Gate: tissue
# =============================================================================

#' Map hpa_tissue vote verdicts onto the gate vocabulary.
#' @keywords internal
#' @noRd
.TISSUE_VERDICT_MAP <- c(
  detected = "detected",
  detected_by_one = "detected_1platform",
  discordant = "discordant",
  concordant_absent = "mismatch",
  insufficient_coverage = "insufficient_coverage"
)

#' Tissue gate values that count as the protein being present.
#' @keywords internal
#' @noRd
.TISSUE_PRESENT <- c("detected", "detected_1platform", "discordant")

#' Look up each gene's presence vote in one tissue.
#'
#' Matches on ensembl_id where the gene has one, then falls back to a
#' case-insensitive symbol match. The fallback is what lets genes outside
#' the concordance atlas still be scored on tissue, since hpa_tissue covers
#' more genes than the atlas does.
#'
#' @return List with \code{call} (character, gate vocabulary) and
#'   \code{n_platforms} (integer, number of platforms detecting the protein
#'   in the tissue, NA when systemic or not in panel).
#' @keywords internal
#' @noRd
.gate_tissue <- function(q, hpa, target) {
  n <- nrow(q)
  if (is.null(target))
    return(list(call = rep("systemic_na", n), n_platforms = rep(NA_integer_, n)))

  sub <- hpa[hpa$tissue == target, , drop = FALSE]

  idx <- rep(NA_integer_, n)
  if ("ensembl_id" %in% names(q)) {
    has_id <- !is.na(q$ensembl_id) & nzchar(q$ensembl_id)
    idx[has_id] <- match(q$ensembl_id[has_id], sub$ensembl_id)
  }
  need <- is.na(idx)
  if (any(need))
    idx[need] <- match(toupper(q$gene_symbol[need]), toupper(sub$gene_symbol))

  call <- unname(.TISSUE_VERDICT_MAP[sub$verdict[idx]])
  call[is.na(idx)] <- "not_in_panel"
  call[is.na(call)] <- "not_in_panel"

  n_platforms <- as.integer(sub$n_detected[idx])

  list(call = call, n_platforms = n_platforms)
}


# =============================================================================
# Gate: concordance
# =============================================================================

#' Mean platforms_agree for each gene in the context tissue (tissue_atlas).
#' NA when no context tissue, no atlas counterpart, or no row.
#' @keywords internal
#' @noRd
.cross_platform_agree <- function(q, ta, target) {
  n <- nrow(q)
  atlas_tissue <- .to_atlas_tissue(target, ta)
  if (is.na(atlas_tissue) || !"platforms_agree" %in% names(ta) ||
      !"ensembl_id" %in% names(q))
    return(rep(NA_real_, n))

  sub <- ta[ta$ts_tissue == atlas_tissue, , drop = FALSE]
  idx <- match(q$ensembl_id, sub$ensembl_id)
  as.numeric(sub$platforms_agree[idx])
}

#' @keywords internal
#' @noRd
.gate_concordance <- function(q, sample_type, n_platforms, cross_agree, thr) {
  n <- nrow(q)
  gc <- if ("gene_class" %in% names(q)) as.character(q$gene_class)
        else rep(NA_character_, n)

  out <- rep("not_in_atlas", n)
  out[!is.na(gc) & gc == "consistently_concordant"] <- "concordant"
  out[!is.na(gc) & gc == "variable"]                <- "variable"
  out[!is.na(gc) & gc == "low_expression"]          <- "low_expression"
  out[!is.na(gc) & gc == "consistently_suppressed"] <- "suppressed"

  suppressed <- out == "suppressed"
  if (!any(suppressed)) return(out)

  # Cellular proteomics blind spot in biofluids
  flag <- function(col) if (col %in% names(q)) !is.na(q[[col]]) & q[[col]] == 1
                        else rep(FALSE, n)
  blind <- sample_type %in% .BIOFLUID_TYPES &
           (flag("is_secreted") | flag("is_projection"))

  # GTEx disagrees in the context tissue and several platforms
  # detect the protein there
  contradicted <- !is.na(cross_agree) & cross_agree < thr$cross_platform_agree_max &
                  !is.na(n_platforms) & n_platforms >= thr$n_platforms_min

  out[suppressed & (blind | contradicted)] <- "suppressed_unconfirmed"
  out
}


# =============================================================================
# Gate: biofluid rescue
# =============================================================================

#' @keywords internal
#' @noRd
.gate_rescue <- function(q, plaus, tissue_call, sample_type) {
  n <- nrow(q)
  out <- rep(FALSE, n)
  out[plaus == "na"] <- NA

  if (!sample_type %in% .RESCUE_TYPES) return(out)

  is_proj <- if ("is_projection" %in% names(q))
    !is.na(q$is_projection) & q$is_projection == 1 else rep(FALSE, n)
  present <- tissue_call == "systemic_na" | tissue_call %in% .TISSUE_PRESENT

  out[plaus == "implausible" & is_proj & present] <- TRUE
  out
}


# =============================================================================
# Gate: specificity and primary tissue (protein-level, IHC panel)
# =============================================================================

#' @keywords internal
#' @noRd
.gate_specificity <- function(q, thr) {
  n <- nrow(q)
  col <- function(nm) if (nm %in% names(q)) as.numeric(q[[nm]]) else rep(NA_real_, n)
  ptau     <- col("protein_tau")
  n_det    <- col("n_tissues_ihc_detected")
  bulk_det <- col("n_tissues_bulk_detected")
  bulk_cov <- col("n_tissues_bulk_covered")

  out <- rep("na", n)
  out[!is.na(ptau) & ptau <  thr$protein_tau_specific] <- "broad"
  out[!is.na(ptau) & ptau >= thr$protein_tau_specific] <- "specific"

  # Too few IHC detections for tau to describe a profile. The bulk platforms
  # decide whether IHC is informative for the gene at all: broad bulk
  # detection means the antibody is missing a protein that is present
  # widely, so the call is na. Otherwise the restriction is consistent
  # across platforms but thin, and the call is sparse.
  sparse <- !is.na(ptau) & !is.na(n_det) & n_det < thr$ihc_min_detected
  bulk_frac <- ifelse(!is.na(bulk_cov) & bulk_cov > 0, bulk_det / bulk_cov, NA_real_)
  ihc_uninformative <- sparse & !is.na(bulk_frac) &
                       bulk_frac >= thr$bulk_breadth_ihc_uninformative
  out[sparse] <- "sparse"
  out[ihc_uninformative] <- "na"
  out
}

#' Tissues tied at the gene's maximum IHC level (HPA names, comma-joined),
#' computed at prep. "multiple" for broad genes, NA for sparse and na.
#' @keywords internal
#' @noRd
.primary_tissue <- function(q, spec) {
  n <- nrow(q)
  out <- rep(NA_character_, n)
  out[spec == "broad"] <- "multiple"
  want <- spec == "specific"
  if (any(want) && "primary_tissue_ihc" %in% names(q))
    out[want] <- as.character(q$primary_tissue_ihc[want])
  out
}
