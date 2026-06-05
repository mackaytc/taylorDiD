# dCDH smoke tests. DIDmultiplegtDYN 2.3.x needs the polars package at runtime,
# so these skip cleanly where it is unavailable.

panel <- readRDS(system.file("extdata", "synthetic_panel.rds",
                             package = "taylorDiD"))

test_that("did.simple.dyn returns a static result (binary treatment)", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  res <- did.simple.dyn(panel, "y", "treat", "year", "id", verbose = FALSE)
  expect_s3_class(res, "taylorDiD_did")
  expect_true(is.finite(res$estimate))
})

test_that("did.event.study.dyn returns the shared event-study object", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  es <- did.event.study.dyn(panel, "y", "treat", "year", "id",
                            n.pre = 3, n.post = 3, verbose = FALSE)
  expect_s3_class(es, "taylorDiD_es")
  expect_equal(nrow(es$coefficients), 6)
  expect_true(is.finite(es$avg.post.effect$estimate))
})

test_that("did.event.study.dyn handles a non-binary treatment", {
  skip_if_not_installed("DIDmultiplegtDYN")
  skip_if_not_installed("polars")
  es <- did.event.study.dyn(panel, "y", "d", "year", "id",
                            n.pre = 3, n.post = 3, normalized = TRUE,
                            verbose = FALSE)
  expect_s3_class(es, "taylorDiD_es")
})
