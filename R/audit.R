# =============================================================================
# audit_call and funnel_class: combiners over triage() output
# =============================================================================
# triage() returns independent gate columns and no verdict. These two
# functions are the only places gates are combined. Both are pure functions
# of the triage() data frame and their arguments.


#' Accept, caution, or decline from triage gates
#'
#' Combines the gate columns returned by \code{\link{triage}} into one call
#' per gene.
#'
#' @section Rules:
#' A gene is \strong{declined} if any of the following holds: the tissue gate
#' is \code{"mismatch"}; plausibility is \code{"implausible"} without
#' \code{biofluid_rescue}; concordance is \code{"suppressed"}; the gene is
#' specific but none of its \code{primary_tissue} set is among
#' \code{accept_tissues} when a context tissue was given; or the gene is
#' absent from every layer
#' (\code{"not_in_atlas"}, \code{"not_in_panel"}, and specificity
#' \code{"na"}).
#'
#' A gene is \strong{accepted} if it is not declined, concordance is
#' \code{"concordant"}, plausibility is \code{"plausible"} or \code{"na"} or
#' rescued, specificity is \code{"specific"}, and the tissue gate is either
#' \code{"systemic_na"} or \code{"detected"} with an accepted tissue in the
#' \code{primary_tissue} set.
#'
#' Everything else is \strong{caution}. This includes
#' \code{"suppressed_unconfirmed"}, \code{"variable"},
#' \code{"low_expression"}, \code{"unknown"} plausibility, and
#' \code{"broad"} specificity.
#'
#' @param df Data frame returned by \code{\link{triage}}.
#' @param accept_tissues Character vector of HPA tissue names. A specific
#'   gene passes the tissue check if any tissue in its \code{primary_tissue}
#'   set is in this vector. Default \code{NULL} derives it from the context
#'   tissue stored in \code{df}, expanded to all brain regions when the
#'   context is one of them. Ignored when no context tissue was given.
#'
#' @return A factor with levels \code{"accept"}, \code{"caution"},
#'   \code{"decline"}, one element per row of \code{df}.
#'
#' @examples
#' \dontrun{
#' res <- triage(c("NEFL", "GFAP", "SPANXA1"), sample_type = "csf",
#'               context_tissue = "Cerebral cortex")
#' res$call <- audit_call(res)
#' table(res$call)
#' }
#'
#' @seealso \code{\link{triage}}, \code{\link{funnel_class}}
#' @export
audit_call <- function(df, accept_tissues = NULL) {
  .check_triage_df(df)
  accept_tissues <- .resolve_accept_tissues(df, accept_tissues)

  rescued <- df$biofluid_rescue %in% TRUE
  plaus_ok <- df$biofluid_plausibility %in% c("plausible", "na") | rescued
  systemic <- df$tissue == "systemic_na"

  primary_ok <- .primary_in(df$primary_tissue, accept_tissues)

  wrong_tissue <- !systemic & df$specificity == "specific" &
                  !is.na(df$primary_tissue) & df$primary_tissue != "multiple" &
                  !primary_ok

  no_support <- df$concordance == "not_in_atlas" &
                df$tissue == "not_in_panel" &
                df$specificity == "na"

  decline <- df$tissue == "mismatch" |
             (df$biofluid_plausibility == "implausible" & !rescued) |
             df$concordance == "suppressed" |
             wrong_tissue |
             no_support

  accept <- !decline &
            df$concordance == "concordant" &
            plaus_ok &
            df$specificity == "specific" &
            (systemic | (df$tissue == "detected" & primary_ok))

  factor(ifelse(decline, "decline", ifelse(accept, "accept", "caution")),
         levels = c("accept", "caution", "decline"))
}


#' First failing gate in a stated order
#'
#' Partitions genes by the first gate they fail when the gates are read in
#' \code{order}. A gene that fails no gate is \code{"pass"}, and the set of
#' \code{"pass"} genes is identical to \code{audit_call() == "accept"} for
#' the same \code{accept_tissues}. Per-gate pass conditions are the accept
#' conditions from \code{\link{audit_call}}, so the partition is a
#' decomposition of the accept rule, not a separate rule set.
#'
#' @param df Data frame returned by \code{\link{triage}}.
#' @param order Character vector naming the gates in funnel order. Any
#'   permutation of \code{"tissue"}, \code{"biofluid_plausibility"},
#'   \code{"concordance"}, \code{"specificity"}.
#' @param accept_tissues As in \code{\link{audit_call}}.
#' @param table Logical. When \code{TRUE}, return the survivor table (one
#'   row per stage with genes entering, failing, and surviving) instead of
#'   the per-gene vector. Default \code{FALSE}.
#'
#' @return A factor with levels \code{order} followed by \code{"pass"}, one
#'   element per row of \code{df}. With \code{table = TRUE}, a data frame
#'   with columns \code{stage}, \code{n_entering}, \code{n_failed},
#'   \code{n_surviving}.
#'
#' @examples
#' \dontrun{
#' res <- triage(genes, sample_type = "csf", context_tissue = "Cerebral cortex")
#' funnel_class(res)
#' funnel_class(res, table = TRUE)
#' }
#'
#' @seealso \code{\link{triage}}, \code{\link{audit_call}}
#' @export
funnel_class <- function(df,
                         order = c("tissue", "biofluid_plausibility",
                                   "concordance", "specificity"),
                         accept_tissues = NULL,
                         table = FALSE) {
  .check_triage_df(df)
  gates <- c("tissue", "biofluid_plausibility", "concordance", "specificity")
  if (!is.character(order) || length(order) != 4L || !setequal(order, gates))
    stop("`order` must be a permutation of: ", paste(gates, collapse = ", "),
         ".", call. = FALSE)
  accept_tissues <- .resolve_accept_tissues(df, accept_tissues)

  pass <- .gate_pass_matrix(df, accept_tissues)[, order, drop = FALSE]

  n <- nrow(df)
  first_fail <- rep("pass", n)
  for (g in rev(order)) first_fail[!pass[, g]] <- g
  out <- factor(first_fail, levels = c(order, "pass"))

  if (!table) return(out)

  entering <- n
  rows <- lapply(order, function(g) {
    failed <- sum(out == g)
    r <- data.frame(stage = g, n_entering = entering, n_failed = failed,
                    n_surviving = entering - failed, stringsAsFactors = FALSE)
    entering <<- entering - failed
    r
  })
  do.call(rbind, rows)
}


# =============================================================================
# Internal helpers
# =============================================================================

#' @keywords internal
#' @noRd
.TRIAGE_COLS <- c("gene_symbol", "concordance", "tissue",
                  "biofluid_plausibility", "biofluid_rescue",
                  "specificity", "primary_tissue")

#' @keywords internal
#' @noRd
.check_triage_df <- function(df) {
  if (!is.data.frame(df))
    stop("`df` must be a data frame returned by triage().", call. = FALSE)
  missing_cols <- setdiff(.TRIAGE_COLS, names(df))
  if (length(missing_cols) > 0L)
    stop("`df` is missing triage() columns: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  invisible(TRUE)
}

#' Derive accept_tissues from the triage context when not supplied.
#' primary_tissue and context_tissue share the hpa_tissue namespace, so no
#' name bridge is needed. Returns NULL when no context tissue applies,
#' meaning primary tissue is not checked.
#' @keywords internal
#' @noRd
.resolve_accept_tissues <- function(df, accept_tissues) {
  if (!is.null(accept_tissues)) {
    if (!is.character(accept_tissues) || length(accept_tissues) == 0L)
      stop("`accept_tissues` must be a character vector.", call. = FALSE)
    return(accept_tissues)
  }
  ctx <- attr(df, "triage_context")
  if (is.null(ctx) || is.null(ctx$context_tissue)) return(NULL)
  if (ctx$context_tissue %in% .BRAIN_REGIONS) return(.BRAIN_REGIONS)
  ctx$context_tissue
}

#' Whether any tissue in a comma-joined primary_tissue set is accepted.
#' NA and "multiple" give FALSE.
#' @keywords internal
#' @noRd
.primary_in <- function(primary_tissue, accept_tissues) {
  if (is.null(accept_tissues)) return(rep(TRUE, length(primary_tissue)))
  vapply(strsplit(as.character(primary_tissue), ", ", fixed = TRUE),
         function(x) !anyNA(x) && any(x %in% accept_tissues), logical(1))
}

#' Per-gate pass conditions as a logical matrix (rows genes, cols gates).
#' These are the accept conditions of audit_call() split by gate.
#' @keywords internal
#' @noRd
.gate_pass_matrix <- function(df, accept_tissues) {
  rescued  <- df$biofluid_rescue %in% TRUE
  systemic <- df$tissue == "systemic_na"
  primary_ok <- .primary_in(df$primary_tissue, accept_tissues)

  cbind(
    tissue = systemic | df$tissue == "detected",
    biofluid_plausibility = df$biofluid_plausibility %in% c("plausible", "na") |
                            rescued,
    concordance = df$concordance == "concordant",
    specificity = df$specificity == "specific" & (systemic | primary_ok)
  )
}
