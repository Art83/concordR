test_that("load_atlas gene returns a data frame", {
  atlas <- load_atlas("gene")
  expect_s3_class(atlas, "data.frame")
  expect_true(nrow(atlas) > 0L)
  expect_true("gene_symbol" %in% names(atlas))
  expect_true("gene_class" %in% names(atlas))
})

test_that("load_atlas tissue returns a data frame", {
  ta <- load_atlas("tissue")
  expect_s3_class(ta, "data.frame")
  expect_true(nrow(ta) > 0L)
  expect_true("ts_tissue" %in% names(ta))
  expect_true("gene_symbol" %in% names(ta))
})

test_that("load_atlas both returns a named list", {
  both <- load_atlas("both")
  expect_type(both, "list")
  expect_named(both, c("gene", "tissue"))
  expect_s3_class(both$gene, "data.frame")
  expect_s3_class(both$tissue, "data.frame")
})

test_that("load_atlas caches across calls", {
  a1 <- load_atlas("gene")
  a2 <- load_atlas("gene")
  # Same object in memory (cached)
  expect_identical(a1, a2)
})

test_that("load_atlas reload forces fresh load", {
  a1 <- load_atlas("gene")
  a2 <- load_atlas("gene", reload = TRUE)
  expect_equal(nrow(a1), nrow(a2))
})

test_that("load_atlas rejects invalid component", {
  expect_error(load_atlas("brain"))
  expect_error(load_atlas("invalid"))
})

test_that("plot_gene rejects multiple genes", {
  expect_error(plot_gene(c("NEFL", "GFAP")), "single gene")
})

test_that("plot_gene rejects unknown genes", {
  expect_error(plot_gene("FAKEGENE_XYZ"), "not found")
})

test_that("plot_gene returns data invisibly", {
  res <- plot_gene("NEFL")
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0L)
})
