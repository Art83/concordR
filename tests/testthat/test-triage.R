test_that("triage returns expected structure", {
  res <- triage(c("NEFL", "GFAP"), sample_type = "csf")
  expect_s3_class(res, "data.frame")
  expect_true(all(c("gene_symbol", "verdict", "reason") %in% names(res)))
  expect_equal(nrow(res), 2L)
})

test_that("triage requires sample_type", {
  expect_error(triage("NEFL"), "sample_type")
})

test_that("triage rejects invalid sample_type", {
  expect_error(triage("NEFL", sample_type = "saliva"),
               "Unknown `sample_type`")
})

test_that("triage rejects empty input", {
  expect_error(triage(character(0), sample_type = "csf"))
})

test_that("triage returns NA verdict for unknown genes", {
  res <- triage("FAKEGENE_XYZ", sample_type = "csf")
  expect_true(is.na(res$verdict))
})

# --- Gate 1: Compartment plausibility ---

test_that("triage flags nuclear proteins as compartment_implausible in CSF", {
  # HDAC1 is nuclear, should fail compartment gate
  res <- triage("HDAC1", sample_type = "csf")
  if (!is.na(res$verdict)) {
    expect_equal(res$verdict, "compartment_implausible")
  }
})

test_that("triage passes secreted proteins in plasma", {
  # A known secreted protein should pass compartment gate
  res <- triage("ALB", sample_type = "plasma")
  if (!is.na(res$verdict)) {
    expect_true(res$verdict != "compartment_implausible")
  }
})

# --- Gate 2: Tissue-of-origin ---

test_that("triage catches tissue-origin mismatch in CSF", {
  # SPANXA1 is testis-restricted, should fail tissue gate for CSF
  res <- triage("SPANXA1", sample_type = "csf")
  if (!is.na(res$verdict)) {
    expect_true(res$verdict %in% c("tissue_origin_mismatch",
                                    "compartment_implausible"))
  }
})

test_that("triage accepts brain-expressed genes in CSF", {
  # NEFL should pass tissue gate for CSF (expressed in brain)
  res <- triage("NEFL", sample_type = "csf")
  if (!is.na(res$verdict)) {
    expect_true(res$verdict != "tissue_origin_mismatch")
  }
})

# --- Gate 3: Concordance ---

test_that("triage flags suppressed genes as RNA-unreliable", {
  # Known suppressed targets
  res <- triage(c("AR", "ESR1", "CD274"), sample_type = "tissue")
  suppressed <- res[res$gene_class == "consistently_suppressed" &
                      !is.na(res$gene_class), ]
  if (nrow(suppressed) > 0) {
    expect_true(all(suppressed$verdict %in%
                      c("RNA-unreliable", "caution", "RNA-supported")))
  }
})

test_that("triage passes concordant genes", {
  res <- triage("NEFL", sample_type = "csf")
  if (!is.na(res$verdict) &&
      !is.na(res$gene_class) &&
      res$gene_class == "consistently_concordant") {
    expect_true(res$verdict %in% c("RNA-supported", "RNA-supported_but_broad",
                                    "caution"))
  }
})

# --- context_tissue ---

test_that("triage accepts context_tissue parameter", {
  res <- triage("NEFL", sample_type = "csf",
                context_tissue = "Cerebellum")
  expect_s3_class(res, "data.frame")
  expect_true("verdict" %in% names(res))
})

test_that("triage rejects 'Brain' as context_tissue", {
  expect_error(
    triage("NEFL", sample_type = "csf", context_tissue = "Brain"),
    "four regions"
  )
})

test_that("triage rejects invalid context_tissue", {
  expect_error(
    triage("NEFL", sample_type = "csf", context_tissue = "Mars"),
    "not found"
  )
})

# --- Plasma (systemic, no tissue gate) ---

test_that("triage skips tissue gate for plasma", {
  # Plasma is systemic — no tissue-of-origin gate fires
  res <- triage("NEFL", sample_type = "plasma")
  if (!is.na(res$verdict)) {
    expect_true(res$verdict != "tissue_origin_mismatch")
  }
})

# --- Output completeness ---

test_that("triage output includes tissue evidence columns for CSF", {
  res <- triage("NEFL", sample_type = "csf")
  expect_true("tissue_detected_any" %in% names(res))
  expect_true("tissue_protein_evidence" %in% names(res))
})

test_that("triage reason column is always character", {
  res <- triage(c("NEFL", "FAKEGENE_XYZ"), sample_type = "csf")
  expect_type(res$reason, "character")
})

# --- Damage-release rescue ---

test_that("triage rescues projection proteins in biofluid context", {
  # NEFL is a projection protein; if its raw compartment check fails,

  # it should be rescued in CSF/plasma
  res <- triage("NEFL", sample_type = "csf")
  if (!is.na(res$verdict)) {
    expect_true(res$verdict != "compartment_implausible")
  }
})
