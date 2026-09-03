# Tests for triage(). Expected values were read from the shipped
# sysdata.rda on 2026-09-02 (HPA IHC, GTEx MS, PaxDb, UniProt, Tabula
# Sapiens v2, Human Brain Cell Atlas v1). A rebuild of the data can change
# them; if a test here fails after a rebuild, check the value against the
# new data before changing the test.

CORE <- c("gene_symbol", "concordance", "tissue", "route",
          "specificity", "primary_tissue")

# -----------------------------------------------------------------------------
# Shape
# -----------------------------------------------------------------------------

test_that("triage returns the six core columns in input order", {
  res <- triage(c("NEFL", "GFAP", "SPANXA1"), sample_type = "csf")
  expect_s3_class(res, "data.frame")
  expect_named(res, CORE)
  expect_equal(res$gene_symbol, c("NEFL", "GFAP", "SPANXA1"))
  expect_equal(nrow(res), 3L)
})

test_that("details = TRUE appends evidence columns after the core", {
  res <- triage("NEFL", sample_type = "csf", details = TRUE)
  expect_equal(names(res)[seq_along(CORE)], CORE)
  expect_true(all(c("protein_tau", "n_tissues_ihc", "tau", "protein_breadth",
                    "localisation", "n_platforms_context") %in% names(res)))
})

test_that("triage_context attribute records the resolved inputs", {
  res <- triage("NEFL", sample_type = "csf", context_tissue = "cerebral cortex")
  ctx <- attr(res, "triage_context")
  expect_equal(ctx$sample_type, "csf")
  expect_equal(ctx$context_tissue, "Cerebral cortex")
})

test_that("gate columns only ever take vocabulary values", {
  res <- triage(c("NEFL", "GFAP", "SPANXA1", "NOT_A_GENE_XYZ"),
                sample_type = "csf", context_tissue = "Cerebral cortex")
  expect_true(all(res$concordance %in%
    c("concordant", "variable", "suppressed", "suppressed_unconfirmed",
      "low_expression", "not_in_atlas")))
  expect_true(all(res$tissue %in%
    c("systemic_na", "detected", "detected_1platform", "discordant",
      "mismatch", "insufficient_coverage", "not_in_panel")))
  expect_true(all(res$route %in%
    c("secreted", "membrane", "intracellular", "unknown", "na")))
  expect_true(all(res$specificity %in% c("specific", "broad", "sparse", "na")))
  expect_false(anyNA(res$concordance))
  expect_false(anyNA(res$tissue))
  expect_false(anyNA(res$route))
  expect_false(anyNA(res$specificity))
})

# -----------------------------------------------------------------------------
# Known genes, CSF, no context tissue
# -----------------------------------------------------------------------------

test_that("NEFL in CSF: intracellular, brain-specific", {
  res <- triage("NEFL", sample_type = "csf")
  expect_equal(res$concordance, "suppressed_unconfirmed")
  expect_equal(res$tissue, "systemic_na")
  expect_equal(res$route, "intracellular")
  expect_equal(res$specificity, "specific")
  expect_equal(res$primary_tissue, "Cerebral cortex, Hippocampus")
})

test_that("GFAP in CSF: intracellular, specific", {
  res <- triage("GFAP", sample_type = "csf")
  expect_equal(res$concordance, "variable")
  expect_equal(res$route, "intracellular")
  expect_equal(res$specificity, "specific")
  # HPA scores GFAP High in cerebellum and placenta. The set reports both.
  expect_equal(res$primary_tissue, "Cerebellum, Placenta")
})

test_that("secreted and membrane routes are read from the flags", {
  res <- triage(c("ALB", "APP"), sample_type = "csf")
  expect_equal(res$route[1], "secreted")
  expect_true(res$route[2] %in% c("secreted", "membrane"))
})

test_that("SPANXA1 in CSF: nuclear via the shared UniProt entry, one IHC hit reads sparse", {
  res <- triage("SPANXA1", sample_type = "csf")
  expect_equal(res$concordance, "low_expression")
  # SPANXA1 and SPANXA2 encode one protein and share a UniProt entry. The
  # localisation build expands the entry to both Ensembl IDs.
  expect_equal(res$route, "intracellular")
  # One IHC detection (Testis), bulk detects it in 1 of 12 covered tissues:
  # restriction is consistent across platforms but the evidence is thin.
  expect_equal(res$specificity, "sparse")
  expect_equal(res$primary_tissue, "Testis")
})

test_that("one IHC hit contradicted by broad bulk detection reads na", {
  # NTN1, SMOC1, SFRP1: one IHC hit each (Testis), bulk-detected in 23 or
  # 24 of 24 covered tissues. IHC is uninformative for these genes.
  res <- triage(c("NTN1", "SMOC1", "SFRP1"), sample_type = "csf")
  expect_equal(res$specificity, rep("na", 3))
  expect_true(all(is.na(res$primary_tissue)))
  d <- triage("NTN1", sample_type = "csf", details = TRUE)
  expect_equal(d$n_tissues_ihc_detected, 1)
  expect_true(d$n_tissues_bulk_detected / d$n_tissues_bulk_covered >= 0.5)
})

# -----------------------------------------------------------------------------
# Context tissue
# -----------------------------------------------------------------------------

test_that("context_tissue switches the tissue gate on", {
  res <- triage(c("NEFL", "GFAP", "SPANXA1"), sample_type = "csf",
                context_tissue = "Cerebral cortex")
  expect_equal(res$tissue, c("detected", "detected", "mismatch"))
  # Other gates unchanged by the context
  base <- triage(c("NEFL", "GFAP", "SPANXA1"), sample_type = "csf")
  expect_equal(res$concordance, base$concordance)
  expect_equal(res$route, base$route)
  expect_equal(res$specificity, base$specificity)
  expect_equal(res$primary_tissue, base$primary_tissue)
})

test_that("context_tissue is case-insensitive and canonicalised", {
  a <- triage("NEFL", sample_type = "csf", context_tissue = "cerebral cortex")
  b <- triage("NEFL", sample_type = "csf", context_tissue = "Cerebral cortex")
  expect_equal(a$tissue, b$tissue)
})

test_that("Brain is rejected with a directed message", {
  expect_error(triage("NEFL", sample_type = "csf", context_tissue = "Brain"),
               "regions")
})

test_that("several context tissues give long output with a context_tissue column", {
  regions <- c("Cerebral cortex", "Hippocampus", "Liver")
  res <- triage(c("NEFL", "GFAP"), sample_type = "csf", context_tissue = regions)
  expect_equal(nrow(res), 6L)
  expect_equal(names(res)[1:2], c("gene_symbol", "context_tissue"))
  expect_equal(res$context_tissue, rep(regions, each = 2))
  expect_equal(res$gene_symbol, rep(c("NEFL", "GFAP"), 3))
  # Context-independent gates repeat across a gene's rows
  for (g in c("NEFL", "GFAP")) {
    sub <- res[res$gene_symbol == g, ]
    expect_length(unique(sub$route), 1L)
    expect_length(unique(sub$specificity), 1L)
    expect_length(unique(sub$primary_tissue), 1L)
  }
  # Each block equals the single-context call
  for (r in regions) {
    single <- triage(c("NEFL", "GFAP"), sample_type = "csf", context_tissue = r)
    block <- res[res$context_tissue == r, setdiff(names(res), "context_tissue")]
    rownames(block) <- NULL
    expect_equal(block, single, ignore_attr = TRUE)
  }
  expect_equal(attr(res, "triage_context")$context_tissue, regions)
})

test_that("duplicate context tissues collapse and one bad name fails the call", {
  res <- triage("NEFL", sample_type = "csf",
                context_tissue = c("cerebral cortex", "Cerebral cortex"))
  expect_equal(nrow(res), 1L)
  expect_false("context_tissue" %in% names(res))
  expect_error(triage("NEFL", sample_type = "csf",
                      context_tissue = c("Cerebral cortex", "Mars")), "Mars")
})

test_that("unknown context_tissue errors and lists the vocabulary", {
  expect_error(triage("NEFL", sample_type = "csf", context_tissue = "Mars"),
               "not found")
  expect_error(triage("NEFL", sample_type = "csf", context_tissue = "Mars"),
               "Cerebellum")
})

test_that("urine implies Kidney", {
  res <- triage("UMOD", sample_type = "urine")
  expect_equal(attr(res, "triage_context")$context_tissue, "Kidney")
  expect_true(res$tissue != "systemic_na")
  # Two IHC hits (Kidney, Liver), bulk in 11 of 24: below the detection
  # floor and not contradicted, so sparse. The set is reported so the
  # caution can be read; Kidney is the only tissue at the maximum level.
  expect_equal(res$specificity, "sparse")
  expect_equal(res$primary_tissue, "Kidney")
})

# -----------------------------------------------------------------------------
# Sample types
# -----------------------------------------------------------------------------

test_that("tissue and cell_lysate have no route rule", {
  res <- triage("NEFL", sample_type = "tissue")
  expect_equal(res$route, "na")
  res <- triage("NEFL", sample_type = "cell_lysate")
  expect_equal(res$route, "na")
})

test_that("subcellular fractions report matched or unmatched", {
  # NEFL is cytoskeletal, not nuclear
  res <- triage("NEFL", sample_type = "nuclear_fraction")
  expect_equal(res$route, "unmatched")
  # HDAC1 is nuclear
  res <- triage("HDAC1", sample_type = "nuclear_fraction")
  expect_equal(res$route, "matched")
  # Unknown localisation stays unknown in any fraction. Readthrough loci
  # have no UniProt entry.
  res <- triage("NME1-NME2", sample_type = "nuclear_fraction")
  expect_equal(res$route, "unknown")
})

test_that("sample_type is validated", {
  expect_error(triage("NEFL"), "sample_type")
  expect_error(triage("NEFL", sample_type = "saliva"), "Unknown")
  expect_error(triage("NEFL", sample_type = c("csf", "plasma")))
})

test_that("genes and details are validated", {
  expect_error(triage(character(0), sample_type = "csf"), "non-empty")
  expect_error(triage(123, sample_type = "csf"), "character")
  expect_error(triage("NEFL", sample_type = "csf", details = NA), "details")
})

# -----------------------------------------------------------------------------
# Genes outside the atlas
# -----------------------------------------------------------------------------

test_that("a gene absent from every layer reads not_in_atlas, not_in_panel, na", {
  res <- triage("NOT_A_GENE_XYZ", sample_type = "csf",
                context_tissue = "Cerebral cortex")
  expect_equal(res$concordance, "not_in_atlas")
  expect_equal(res$tissue, "not_in_panel")
  expect_equal(res$route, "unknown")
  expect_equal(res$specificity, "na")
  expect_true(is.na(res$primary_tissue))
})

test_that("a gene outside the atlas but in the presence vote still gets a tissue call", {
  # hpa_tissue covers more genes than concordance_atlas. Pick one from the
  # difference at test time rather than hard-coding a symbol.
  hpa <- load_atlas("hpa")
  atlas <- load_atlas("gene")
  extra <- setdiff(unique(hpa$gene_symbol[hpa$tissue == "Liver" &
                                            hpa$verdict == "detected"]),
                   atlas$gene_symbol)
  extra <- extra[!is.na(extra) & nzchar(extra)]
  skip_if(length(extra) == 0L, "no gene outside the atlas with a Liver vote")
  res <- triage(extra[1L], sample_type = "plasma", context_tissue = "Liver")
  expect_equal(res$concordance, "not_in_atlas")
  expect_equal(res$tissue, "detected")
})

test_that("input order and duplicates are preserved", {
  res <- triage(c("GFAP", "NEFL", "GFAP"), sample_type = "csf")
  expect_equal(res$gene_symbol, c("GFAP", "NEFL", "GFAP"))
  expect_equal(res$concordance[1], res$concordance[3])
})
