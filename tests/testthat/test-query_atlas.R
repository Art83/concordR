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
