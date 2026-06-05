# did2s is experimental; these confirm the generalized data prep and that the
# wrapper returns the shared result object.

panel <- readRDS(system.file("extdata", "synthetic_panel.rds",
                             package = "taylorDiD"))

test_that("prepare.did2s.data generalizes columns and builds event dummies", {
  prepared <- prepare.did2s.data(
    data.table::as.data.table(panel), "y", "g", "year", "id",
    pre.window = -4:-1, post.window = 0:3)

  expect_identical(prepared$status, "ok")
  expect_true(all(c("treat.did2s", "rel.year.did2s") %in% names(prepared$data)))
  # The omitted reference period (-1) is never a regressor.
  expect_false(-1L %in% prepared$event.map$event.time)
  # One 0/1 dummy column per mapped event time.
  expect_true(all(prepared$event.map$term %in% names(prepared$data)))
})

test_that("did2s_event_study returns the shared event-study object", {
  skip_if_not_installed("did2s")
  es <- did2s_event_study(panel, "y", "g", "year", "id",
                          pre.window = -4:-1, post.window = 0:3,
                          verbose = FALSE)
  expect_s3_class(es, "taylorDiD_es")
  expect_identical(es$estimator, "did2s")
  expect_true(is.finite(es$pretrends.test$p.value))
})
