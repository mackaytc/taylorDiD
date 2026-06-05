# BJS imputation smoke tests on the bundled synthetic panel.

panel <- readRDS(system.file("extdata", "synthetic_panel.rds",
                             package = "taylorDiD"))

test_that("did.simple returns a static result object", {
  res <- did.simple(panel, "y", "g", "year", "id", verbose = FALSE)
  expect_s3_class(res, "taylorDiD_did")
  expect_true(is.finite(res$estimate))
  expect_true(is.finite(res$pre.mean.treated))
})

test_that("did.event.study returns the shared event-study object", {
  es <- did.event.study(panel, "y", "g", "year", "id",
                        pre.window = -4:-1, post.window = 0:3, verbose = FALSE)
  expect_s3_class(es, "taylorDiD_es")
  expect_true(all(c("coefficients", "pretrends.test", "avg.post.effect",
                    "pre.mean.treated", "n.obs") %in% names(es)))
  expect_true(is.finite(es$pretrends.test$p.value))
  expect_true(is.finite(es$avg.post.effect$estimate))

  # The simulated post-treatment effect ramps positive.
  coefs <- es$coefficients
  post <- coefs$estimate[coefs$event.time >= 0]
  expect_gt(mean(post), 0)
})

test_that("did.event.study with trends.lin runs and flags itself", {
  es <- did.event.study(panel, "y", "g", "year", "id",
                        pre.window = -4:-1, post.window = 0:3,
                        trends.lin = TRUE, verbose = FALSE)
  expect_s3_class(es, "taylorDiD_es")
  expect_true(isTRUE(es$trends.lin))
})

test_that("BJS treats negative never-treated codes like Inf", {
  inf_panel <- data.table::copy(data.table::as.data.table(panel))
  neg_panel <- data.table::copy(inf_panel)
  neg_panel[is.infinite(g), g := -1]   # never-treated as a negative code

  es_inf <- did.event.study(inf_panel, "y", "g", "year", "id",
                            pre.window = -4:-1, post.window = 0:3,
                            verbose = FALSE)
  es_neg <- did.event.study(neg_panel, "y", "g", "year", "id",
                            pre.window = -4:-1, post.window = 0:3,
                            verbose = FALSE)

  expect_equal(es_neg$coefficients$estimate, es_inf$coefficients$estimate)
  expect_equal(es_neg$avg.post.effect$estimate, es_inf$avg.post.effect$estimate)
  expect_equal(es_neg$n.treated, es_inf$n.treated)
})
