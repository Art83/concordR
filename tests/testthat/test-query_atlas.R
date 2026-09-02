test_that("query_atlas returns expected structure", {
  res <- query_atlas("NEFL")
  expect_s3_class(res, "data.frame")
  expect_true("gene_symbol" %in% names(res))
  expect_true("gene_class" %in% names(res))
  expect_true("found" %in% names(res))
  expect_equal(nrow(res), 1L)
  expect_true(res$found)
})

test_that("query_atlas handles multiple genes", {
  genes <- c("NEFL", "GFAP", "TREM2")
  res <- query_atlas(genes)
  expect_equal(nrow(res), 3L)
  expect_true(all(res$found))
})

test_that("query_atlas handles missing genes with keep", {
  res <- query_atlas(c("NEFL", "FAKEGENE_XYZ"), missing = "keep")
  expect_equal(nrow(res), 2L)
  expect_true(res$found[toupper(res$gene_symbol) == "NEFL"])
  expect_false(res$found[toupper(res$gene_symbol) == "FAKEGENE_XYZ"])
})

test_that("query_atlas handles missing genes with drop", {
  res <- query_atlas(c("NEFL", "FAKEGENE_XYZ"), missing = "drop")
  expect_equal(nrow(res), 1L)
  expect_equal(toupper(res$gene_symbol), "NEFL")
})

test_that("query_atlas handles missing genes with warn", {
  expect_message(
    query_atlas(c("NEFL", "FAKEGENE_XYZ"), missing = "warn"),
    "not found"
  )
})

test_that("query_atlas rejects empty input", {
  expect_error(query_atlas(character(0)))
  expect_error(query_atlas(123))
})

test_that("query_atlas preserves input order", {
  genes <- c("GFAP", "NEFL", "TREM2")
  res <- query_atlas(genes)
  expect_equal(toupper(res$gene_symbol), toupper(genes))
})

test_that("query_atlas is case-insensitive", {
  res1 <- query_atlas("NEFL")
  res2 <- query_atlas("nefl")
  expect_equal(res1$gene_class, res2$gene_class)
})

test_that("query_atlas collapses case variants and repeats to one row", {
  res <- query_atlas(c("NEFL", "nefl", " NEFL ", "GFAP", "GFAP"))
  expect_equal(nrow(res), 2L)
  expect_equal(res$gene_symbol, c("NEFL", "GFAP"))
})

test_that("query_atlas missing rows keep the input spelling and NA typed columns", {
  res <- query_atlas(c("GFAP", "FakeGene_xyz", "NEFL"))
  expect_equal(res$gene_symbol, c("GFAP", "FakeGene_xyz", "NEFL"))
  expect_equal(res$found, c(TRUE, FALSE, TRUE))
  expect_true(is.na(res$gene_class[2]))
  expect_type(res$suppression_rate, "double")
})

test_that("query_tissue returns tissue-level data", {
  res <- query_tissue("NEFL")
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0L)
  expect_true("ts_tissue" %in% names(res))
  expect_true("gene_symbol" %in% names(res))
})

test_that("query_tissue filters by tissue", {
  res <- query_tissue("NEFL", tissues = "Cerebellum")
  expect_true(all(res$ts_tissue == "Cerebellum"))
})

test_that("query_tissue tissue matching is case-insensitive", {
  a <- query_tissue("NEFL", tissues = "cerebellum")
  b <- query_tissue("NEFL", tissues = "Cerebellum")
  expect_equal(nrow(a), nrow(b))
  expect_equal(nrow(a), 1L)
})

test_that("query_tissue rejects unknown tissues with the vocabulary", {
  expect_error(query_tissue("NEFL", tissues = "Mars"), "Unknown tissue")
  expect_error(query_tissue("NEFL", tissues = "Mars"), "Cerebellum")
})

test_that("query_tissue orders by input gene then tissue and reports absent genes", {
  expect_message(res <- query_tissue(c("GFAP", "FAKEGENE_XYZ", "NEFL")), "no rows")
  expect_equal(unique(res$gene_symbol), c("GFAP", "NEFL"))
  g <- res[res$gene_symbol == "GFAP", ]
  expect_equal(g$ts_tissue, sort(g$ts_tissue))
  expect_silent(query_tissue(c("GFAP", "FAKEGENE_XYZ"), missing = "silent"))
})
