# =============================================================================
# load_atlas — session-cached atlas loader
# =============================================================================

# Module-level cache (package environment)
.atlas_env <- new.env(parent = emptyenv())


#' Load the RNA--protein concordance atlas
#'
#' Loads the precomputed gene-level atlas into a session cache. Subsequent
#' calls within the same R session return the cached object without
#' reloading. All other concordR functions call this internally when
#' needed, so explicit use is optional.
#'
#' @section Atlas contents:
#' The atlas is a data frame with one row per gene (up to 11,154 genes)
#' and the following columns:
#' \describe{
#'   \item{\code{gene_symbol}}{HGNC symbol.}
#'   \item{\code{gene_class}}{Character: \code{"consistently_concordant"},
#'     \code{"variable"}, \code{"consistently_suppressed"},
#'     \code{"low_expression"}.}
#'   \item{\code{protein_confidence}}{Numeric [0, 1]. Equal to
#'     \code{1 - suppression_rate}. The core metric.}
#'   \item{\code{detection_rate}}{Numeric [0, 1]. Fraction of observations
#'     where protein was detected across all tissues.}
#'   \item{\code{mean_rna_rank}}{Numeric [0, 1]. Mean RNA expression
#'     percentile across all tissues.}
#'   \item{\code{n_tissues}}{Integer. Number of tissues with data for
#'     this gene in the atlas.}
#'   \item{\code{mechanistic_tier}}{Character: \code{"mRNA_decay"},
#'     \code{"translational_block"}, \code{"post_translational"},
#'     \code{"unexplained"}, \code{NA} (for non-suppressed genes).}
#'   \item{\code{is_secreted}}{Binary (0/1). HPA: secreted protein.}
#'   \item{\code{is_membrane}}{Binary (0/1). HPA: plasma membrane.}
#'   \item{\code{is_nuclear}}{Binary (0/1). HPA: nuclear localisation.}
#'   \item{\code{is_cytoplasmic}}{Binary (0/1). HPA: cytoplasmic.}
#'   \item{\code{is_mitochondrial}}{Binary (0/1). HPA: mitochondrial.}
#'   \item{\code{is_er_golgi}}{Binary (0/1). HPA: ER/Golgi.}
#'   \item{\code{is_cytoskeleton}}{Binary (0/1). HPA: cytoskeletal.}
#'   \item{\code{is_multilocal}}{Binary (0/1). Multiple localisations.}
#' }
#'
#' @section Tissue-resolved profiles:
#' An optional extended atlas with per-gene, per-tissue concordance data
#' is available via \code{component = "tissue"}. This is a larger object
#' (~50 MB) downloaded from Zenodo on first use.
#'
#' @param component Character. \code{"gene"} (default) loads the gene-level
#'   summary atlas. \code{"tissue"} loads the tissue-resolved concordance
#'   matrix. \code{"both"} loads both.
#' @param reload Logical. Force reload even if cached. Default \code{FALSE}.
#'
#' @return For \code{"gene"}: a data frame (one row per gene). For
#'   \code{"tissue"}: a data frame (one row per gene x tissue). For
#'   \code{"both"}: a named list with elements \code{gene} and \code{tissue}.
#'
#' @examples
#' \dontrun{
#' atlas <- load_atlas()
#' head(atlas)
#'
#' # Tissue-resolved profiles
#' tissue_atlas <- load_atlas("tissue")
#' }
#'
#' @export
load_atlas <- function(component = c("gene", "tissue", "both"),
                       reload    = FALSE) {
  component <- match.arg(component)

  # --- Gene-level atlas (ships with package or from sysdata.rda) ---
  if (component %in% c("gene", "both")) {
    if (isTRUE(reload) || is.null(.atlas_env$gene)) {
      if (exists("concordance_atlas", envir = asNamespace("concordR"),
                 inherits = FALSE)) {
        .atlas_env$gene <- get("concordance_atlas",
                               envir = asNamespace("concordR"))
      } else {
        stop("Gene-level atlas not found. ",
             "Ensure the package was installed with data intact.")
      }
    }
  }

  # --- Tissue-resolved atlas (from Zenodo, lazy-loaded) ---
  if (component %in% c("tissue", "both")) {
    if (isTRUE(reload) || is.null(.atlas_env$tissue)) {
      cache_dir <- .concordR_cache_dir()
      tissue_file <- file.path(cache_dir, "concordR_tissue_atlas.rds")

      if (!file.exists(tissue_file)) {
        message("Downloading tissue-resolved atlas from Zenodo ",
                "(one-time, ~50 MB)...")
        .download_tissue_atlas(tissue_file)
      }
      .atlas_env$tissue <- readRDS(tissue_file)
    }
  }

  switch(component,
    gene   = .atlas_env$gene,
    tissue = .atlas_env$tissue,
    both   = list(gene = .atlas_env$gene, tissue = .atlas_env$tissue)
  )
}


#' Cache directory for concordR
#' @keywords internal
#' @noRd
.concordR_cache_dir <- function() {
  tryCatch(
    tools::R_user_dir("concordR", which = "cache"),
    error = function(e) file.path(path.expand("~"), ".cache", "concordR")
  )
}


#' Download tissue atlas from Zenodo
#' @keywords internal
#' @noRd
.download_tissue_atlas <- function(dest_file) {
  # TODO: update record ID and filename once Zenodo deposit is finalised
  record_id <- "XXXXXXX"
  filename  <- "concordR_tissue_atlas.rds"
  url <- paste0("https://zenodo.org/api/records/", record_id,
                "/files/", filename, "/content")

  dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = 300L)

  code <- tryCatch(
    utils::download.file(url, dest_file, mode = "wb", quiet = FALSE),
    error = function(e) NA_integer_
  )

  if (is.na(code) || code != 0L || !file.exists(dest_file))
    stop("Download failed. Check network connectivity.\nURL: ", url)

  invisible(dest_file)
}


#' Internal helper: ensure gene-level atlas is loaded
#' @keywords internal
#' @noRd
.ensure_atlas <- function() {
  if (is.null(.atlas_env$gene)) load_atlas("gene")
  .atlas_env$gene
}
