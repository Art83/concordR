# =============================================================================
# plot_gene — tissue-resolved concordance profile
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
#' @seealso \code{\link{query_atlas}}, \code{\link{plot_audit}}
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

  # Filter to gene — match on gene_symbol column
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
  par(mar = c(8, 4, 3, 4), las = 2)

  # Use detection rate as bar heights if available, else RNA rank
  if (!is.null(det_col)) {
    bar_vals <- gene_data[[det_col]]
    ylab <- "Protein detection rate"
  } else if (!is.null(rna_col)) {
    bar_vals <- gene_data[[rna_col]]
    ylab <- "RNA expression rank"
  } else {
    stop("No detection_rate or rna_rank column in tissue atlas.")
  }

  # Colour by detection level
  cols <- ifelse(bar_vals >= 0.5, "#2E86AB",
          ifelse(bar_vals > 0, "#F6AE2D", "#E8432A"))

  bp <- barplot(bar_vals,
                names.arg = gene_data[[tissue_col]],
                col       = cols,
                border    = NA,
                ylim      = c(0, 1),
                ylab      = ylab,
                main      = if (show_class)
                  paste0(gene, " (", gene_info$gene_class[1], ")")
                else gene,
                cex.names = 0.7,
                ...)

  # Overlay RNA rank as points if we plotted detection as bars
  if (!is.null(det_col) && !is.null(rna_col)) {
    points(bp, gene_data[[rna_col]], pch = 18, col = "#1B1B1B", cex = 1.2)
    legend("topright",
           legend = c("Detection >= 50%", "Detection 1-50%",
                      "Not detected", "RNA rank"),
           fill   = c("#2E86AB", "#F6AE2D", "#E8432A", NA),
           pch    = c(NA, NA, NA, 18),
           col    = c(NA, NA, NA, "#1B1B1B"),
           border = NA, bty = "n", cex = 0.8)
  }

  invisible(gene_data)
}


# =============================================================================
# plot_audit — standard vs corrected enrichment comparison
# =============================================================================

#' Plot standard vs concordance-corrected enrichment
#'
#' Visualises the output of \code{\link{audit_enrichment}} as a paired
#' comparison: each pathway is shown with its original and corrected
#' p-value, connected by a segment coloured by whether the pathway
#' survived correction.
#'
#' @param audit_result Data frame returned by \code{\link{audit_enrichment}}.
#'   Must contain \code{padj_corrected}.
#' @param top_n Integer. Max pathways to display. Default \code{20}.
#' @param alpha Numeric. Significance line. Default \code{0.05}.
#' @param label_col Character or \code{NULL}. Pathway label column.
#'   Auto-detected if \code{NULL}.
#'
#' @return Invisible. Called for side effect (plot).
#'
#' @seealso \code{\link{audit_enrichment}}
#' @export
plot_audit <- function(audit_result,
                       top_n     = 20L,
                       alpha     = 0.05,
                       label_col = NULL) {

  if (!"padj_corrected" %in% names(audit_result))
    stop("audit_result must contain `padj_corrected`. ",
         "Re-run audit_enrichment() with `universe` supplied.")

  if (is.null(label_col)) {
    label_col <- .detect_col(audit_result,
      c("ID", "Description", "term", "pathway", "set", "Term_name"),
      "pathway label")
  }

  pval_col <- .detect_col(audit_result,
    c("p.adjust", "padj", "pvalue", "p_value", "FDR", "qvalue"),
    "original p-value")

  audit_result <- audit_result[order(audit_result[[pval_col]]), ]
  audit_result <- head(audit_result, top_n)
  n <- nrow(audit_result)

  if (n == 0L) { message("No pathways to plot."); return(invisible(NULL)) }

  labels <- as.character(audit_result[[label_col]])
  labels <- ifelse(nchar(labels) > 45,
                   paste0(substr(labels, 1, 42), "..."), labels)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(5, 18, 3, 2), las = 1)

  y_pos <- seq_len(n)

  orig_log <- -log10(audit_result[[pval_col]])
  corr_log <- -log10(audit_result$padj_corrected)
  xlim <- c(0, max(c(orig_log, corr_log), na.rm = TRUE) * 1.1)

  plot(NULL, xlim = xlim, ylim = c(0.5, n + 0.5),
       xlab = expression(-log[10](p[adj])),
       ylab = "", yaxt = "n",
       main = "Standard vs Concordance-Corrected Enrichment")

  axis(2, at = y_pos, labels = rev(labels), tick = FALSE, cex.axis = 0.7)
  abline(v = -log10(alpha), lty = 2, col = "grey50")

  status <- rev(audit_result$status_change)
  seg_col <- ifelse(status == "lost", "#E8432A",
             ifelse(status == "robust", "#2E86AB", "grey70"))

  segments(x0 = rev(orig_log), x1 = rev(corr_log),
           y0 = y_pos, y1 = y_pos, col = seg_col, lwd = 2)
  points(rev(orig_log), y_pos, pch = 16, col = "grey30", cex = 1.0)
  points(rev(corr_log), y_pos, pch = 17, col = seg_col, cex = 1.0)

  legend("bottomright",
         legend = c("Original", "Corrected", "Lost", "Robust"),
         pch    = c(16, 17, NA, NA),
         lty    = c(NA, NA, 1, 1),
         lwd    = c(NA, NA, 2, 2),
         col    = c("grey30", "grey30", "#E8432A", "#2E86AB"),
         bty    = "n", cex = 0.8)

  invisible(audit_result)
}
