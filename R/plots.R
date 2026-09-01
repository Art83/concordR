# =============================================================================
# plot_gene, tissue-resolved concordance profile
# =============================================================================

#' Plot tissue-resolved concordance profile for a gene
#'
#' Displays a gene's concordance landscape across tissues: RNA expression
#' rank, protein detection rate, and concordance status per tissue.
#'
#' @param gene Character. A single HGNC gene symbol.
#' @param atlas_tissue Optional. Tissue-resolved atlas data frame (from
#'   \code{load_atlas("tissue")}). If \code{NULL}, attempts to load
#'   automatically.
#' @param sort_by Character. \code{"detection"} (default), \code{"rna"},
#'   or \code{"name"}.
#' @param show_class Logical. Show the gene's overall class in the title.
#'   Default \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link[graphics]{barplot}}.
#'
#' @return Invisible data frame of tissue-level data for the gene.
#'
#' @examples
#' \dontrun{
#' plot_gene("TREM2")
#' plot_gene("NEFL", sort_by = "rna")
#' }
#'
#' @seealso \code{\link{query_atlas}}, \code{\link{query_tissue}}
#' @export
plot_gene <- function(gene,
                      atlas_tissue = NULL,
                      sort_by      = c("detection", "rna", "name"),
                      show_class   = TRUE,
                      ...) {

  sort_by <- match.arg(sort_by)

  if (length(gene) != 1L)
    stop("`gene` must be a single gene symbol.")

  gene_info <- query_atlas(gene, missing = "warn")
  if (!gene_info$found[1])
    stop("Gene '", gene, "' not found in atlas.")

  if (is.null(atlas_tissue)) {
    atlas_tissue <- tryCatch(load_atlas("tissue"),
                             error = function(e) NULL)
    if (is.null(atlas_tissue))
      stop("Tissue-resolved atlas not available. ",
           "Supply via `atlas_tissue` or run load_atlas('tissue').")
  }

  # Filter to gene - match on gene_symbol column
  sym_col <- if ("gene_symbol" %in% names(atlas_tissue)) "gene_symbol"
             else if ("gene" %in% names(atlas_tissue)) "gene"
             else stop("Tissue atlas has no gene_symbol/gene column.")

  gene_data <- atlas_tissue[toupper(atlas_tissue[[sym_col]]) == toupper(gene), ,
                             drop = FALSE]
  if (nrow(gene_data) == 0L)
    stop("No tissue data for '", gene, "'.")

  # Detect tissue and value columns
  tissue_col <- if ("ts_tissue" %in% names(gene_data)) "ts_tissue"
                else if ("tissue" %in% names(gene_data)) "tissue"
                else stop("No tissue column found.")

  det_col <- if ("detect_fraction" %in% names(gene_data)) "detect_fraction"
             else if ("detection_rate" %in% names(gene_data)) "detection_rate"
             else if ("y_bin" %in% names(gene_data)) "y_bin"
             else NULL

  rna_col <- if ("rna_rank_max" %in% names(gene_data)) "rna_rank_max"
             else if ("mean_rna_rank" %in% names(gene_data)) "mean_rna_rank"
             else if ("rna_rank" %in% names(gene_data)) "rna_rank"
             else NULL

  # Sort
  gene_data <- switch(sort_by,
    detection = if (!is.null(det_col))
      gene_data[order(gene_data[[det_col]], decreasing = TRUE), ] else gene_data,
    rna  = if (!is.null(rna_col))
      gene_data[order(gene_data[[rna_col]], decreasing = TRUE), ] else gene_data,
    name = gene_data[order(gene_data[[tissue_col]]), ]
  )

  # --- Plot ---
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(10, 4, 3, 2), las = 2)
  
  # Bars = RNA expression (always visible); colour = protein detection status
  if (!is.null(rna_col)) {
    bar_vals <- gene_data[[rna_col]]
  } else if (!is.null(det_col)) {
    bar_vals <- gene_data[[det_col]]
  } else {
    stop("No rna_rank or detection_rate column in tissue atlas.")
  }
  
  # Colour by protein detection
  if (!is.null(det_col)) {
    det_vals <- gene_data[[det_col]]
    cols <- ifelse(det_vals >= 0.5, "#2E86AB",
                   ifelse(det_vals > 0,    "#F6AE2D",
                          "#E8432A"))
  } else {
    cols <- rep("#2E86AB", length(bar_vals))
  }
  
  bp <- barplot(bar_vals,
                names.arg = gene_data[[tissue_col]],
                col = cols,
                border = NA,
                ylim = c(0, max(bar_vals, na.rm = TRUE) * 1.1),
                ylab = "RNA expression rank",
                main = if (show_class)
                  paste0(gene, " (", gene_info$gene_class[1], ")")
                else gene,
                cex.names = 0.7,
                ...)
  
  # Overlay protein detection as points
  if (!is.null(det_col) && !is.null(rna_col)) {
    points(bp, gene_data[[det_col]], pch = 18, col = "#1B1B1B", cex = 1.2)
  }
  
  # Legend below x-axis labels
  par(xpd = TRUE)
  legend_y <- -0.35 * max(bar_vals, na.rm = TRUE)
  legend(x = mean(range(bp)), y = legend_y,
         xjust = 0.5,
         title = "Bar = RNA rank, coloured by protein detection:",
         legend = c(">= 50% cell types", "1-50% cell types",
                    "Not detected", "Protein detection rate"),
         fill = c("#2E86AB", "#F6AE2D", "#E8432A", NA),
         pch = c(NA, NA, NA, 18),
         col = c(NA, NA, NA, "#1B1B1B"),
         border = NA, bty = "n", cex = 0.70, horiz = TRUE)
  par(xpd = FALSE)

  invisible(gene_data)
}
