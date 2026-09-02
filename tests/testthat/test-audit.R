# Tests for audit_call() and funnel_class().

# A synthetic triage() frame lets the rules be tested without the atlas.
mk <- function(...) {
  df <- data.frame(..., stringsAsFactors = FALSE)
  attr(df, "triage_context") <- list(sample_type = "csf",
                                     context_tissue = "Cerebral cortex")
  df
}

full <- function(gene, concordance = "concordant", tissue = "detected",
                 plaus = "plausible", rescue = FALSE, spec = "specific",
                 primary = "Cerebral cortex") {
  mk(gene_symbol = gene, concordance = concordance, tissue = tissue,
     biofluid_plausibility = plaus, biofluid_rescue = rescue,
     specificity = spec, primary_tissue = primary)
}

# -----------------------------------------------------------------------------
# audit_call
# -----------------------------------------------------------------------------

test_that("audit_call accepts a gene that passes every gate", {
  expect_equal(as.character(audit_call(full("A"))), "accept")
})

test_that("audit_call declines on each decline rule", {
  expect_equal(as.character(audit_call(full("A", tissue = "mismatch"))), "decline")
  expect_equal(as.character(audit_call(full("A", plaus = "implausible"))), "decline")
  expect_equal(as.character(audit_call(full("A", concordance = "suppressed"))), "decline")
  expect_equal(as.character(audit_call(full("A", primary = "Liver"))), "decline")
  expect_equal(as.character(audit_call(full("A", concordance = "not_in_atlas",
                                            tissue = "not_in_panel", spec = "na",
                                            primary = NA))), "decline")
})

test_that("rescue converts implausible from decline to not-decline", {
  res <- audit_call(full("A", plaus = "implausible", rescue = TRUE))
  expect_equal(as.character(res), "accept")
})

test_that("plausibility na passes for sample types with no rule", {
  df <- full("A", plaus = "na", rescue = NA)
  attr(df, "triage_context")$sample_type <- "tissue"
  expect_equal(as.character(audit_call(df)), "accept")
})

test_that("caution for the intermediate values", {
  expect_equal(as.character(audit_call(full("A", concordance = "suppressed_unconfirmed"))), "caution")
  expect_equal(as.character(audit_call(full("A", concordance = "variable"))), "caution")
  expect_equal(as.character(audit_call(full("A", concordance = "low_expression"))), "caution")
  expect_equal(as.character(audit_call(full("A", plaus = "unknown"))), "caution")
  expect_equal(as.character(audit_call(full("A", spec = "broad", primary = "multiple"))), "caution")
  expect_equal(as.character(audit_call(full("A", spec = "sparse", primary = NA))), "caution")
  expect_equal(as.character(audit_call(full("A", spec = "na", primary = NA))), "caution")
  expect_equal(as.character(audit_call(full("A", tissue = "detected_1platform"))), "caution")
  expect_equal(as.character(audit_call(full("A", tissue = "discordant"))), "caution")
})

test_that("brain context accepts any brain region as primary tissue", {
  expect_equal(as.character(audit_call(full("A", primary = "Hippocampus"))), "accept")
  expect_equal(as.character(audit_call(full("A", primary = "Substantia nigra"))), "accept")
})

test_that("primary_tissue is a set: any accepted member passes", {
  expect_equal(as.character(audit_call(full("A", primary = "Cerebellum, Placenta"))), "accept")
  expect_equal(as.character(audit_call(full("A", primary = "Liver, Placenta, Testis"))), "decline")
  expect_equal(as.character(audit_call(full("A", primary = "Liver", tissue = "detected"),
                                       accept_tissues = c("Liver", "Kidney"))), "accept")
})

test_that("systemic context does not check primary tissue", {
  df <- full("A", tissue = "systemic_na", primary = "Liver")
  attr(df, "triage_context")$context_tissue <- NULL
  expect_equal(as.character(audit_call(df)), "accept")
})

test_that("accept_tissues overrides the derived set", {
  df <- full("A", primary = "Liver")
  expect_equal(as.character(audit_call(df, accept_tissues = "Liver")), "accept")
  expect_error(audit_call(df, accept_tissues = 1), "character")
})

test_that("audit_call returns a factor with fixed levels", {
  res <- audit_call(full("A"))
  expect_s3_class(res, "factor")
  expect_equal(levels(res), c("accept", "caution", "decline"))
})

test_that("audit_call rejects frames that are not from triage", {
  expect_error(audit_call(data.frame(x = 1)), "missing")
  expect_error(audit_call(list()), "data frame")
})

# -----------------------------------------------------------------------------
# funnel_class
# -----------------------------------------------------------------------------

test_that("funnel_class reports the first failing gate in the stated order", {
  df <- full("A", tissue = "mismatch", concordance = "suppressed")
  expect_equal(as.character(funnel_class(df)), "tissue")
  expect_equal(as.character(funnel_class(df,
    order = c("concordance", "tissue", "biofluid_plausibility", "specificity"))),
    "concordance")
})

test_that("funnel_class levels are the order followed by pass", {
  ord <- c("specificity", "concordance", "tissue", "biofluid_plausibility")
  res <- funnel_class(full("A"), order = ord)
  expect_equal(levels(res), c(ord, "pass"))
  expect_equal(as.character(res), "pass")
})

test_that("funnel_class validates order", {
  expect_error(funnel_class(full("A"), order = c("tissue", "concordance")),
               "permutation")
  expect_error(funnel_class(full("A"), order = c("a", "b", "c", "d")),
               "permutation")
})

test_that("funnel_class table has one row per stage and conserves counts", {
  df <- rbind(full("A"),
              full("B", tissue = "mismatch"),
              full("C", plaus = "implausible"),
              full("D", concordance = "variable"),
              full("E", spec = "broad", primary = "multiple"))
  attr(df, "triage_context") <- attr(full("A"), "triage_context")
  tab <- funnel_class(df, table = TRUE)
  expect_equal(nrow(tab), 4L)
  expect_equal(tab$stage, c("tissue", "biofluid_plausibility",
                            "concordance", "specificity"))
  expect_equal(tab$n_entering[1], 5L)
  expect_equal(tab$n_failed, c(1L, 1L, 1L, 1L))
  expect_equal(tab$n_surviving[4], 1L)
  expect_equal(tab$n_entering[-1], tab$n_surviving[-4])
})

# -----------------------------------------------------------------------------
# Invariant: pass in the funnel is accept in audit_call
# -----------------------------------------------------------------------------

test_that("funnel pass equals audit accept on real triage output", {
  genes <- c("NEFL", "GFAP", "SPANXA1", "APOE", "APP", "MAPT", "SNCA",
             "TREM2", "ALB", "GAPDH", "UMOD", "NOT_A_GENE_XYZ")
  for (ctx in list(NULL, "Cerebral cortex", "Hippocampus", "Liver")) {
    res <- triage(genes, sample_type = "csf", context_tissue = ctx)
    a <- audit_call(res) == "accept"
    p <- funnel_class(res) == "pass"
    expect_equal(p, a, info = paste("context:", if (is.null(ctx)) "NULL" else ctx))
    for (ord in list(c("concordance", "specificity", "tissue", "biofluid_plausibility"),
                     c("specificity", "biofluid_plausibility", "concordance", "tissue"))) {
      expect_equal(funnel_class(res, order = ord) == "pass", a)
    }
  }
})
