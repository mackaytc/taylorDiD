# Joint pre-trends F-test and the pre-treatment mean.

test_that("joint_pretrends_F diagonal matches the sum of squared t-stats", {
  est <- c(0.1, -0.05, 0.2)
  se <- c(0.08, 0.07, 0.09)
  res <- joint_pretrends_F(est = est, se = se)

  wald <- sum((est / se)^2)
  expect_equal(res$f.stat, wald / 3)
  expect_equal(res$df, 3)
  expect_equal(res$p.value,
               pf(wald / 3, df1 = 3, df2 = Inf, lower.tail = FALSE))
})

test_that("joint_pretrends_F full vcov equals diagonal when vcov is diagonal", {
  est <- c(0.1, -0.05, 0.2)
  se <- c(0.08, 0.07, 0.09)
  expect_equal(joint_pretrends_F(est, vcov = diag(se^2))$f.stat,
               joint_pretrends_F(est, se = se)$f.stat)
})

test_that("joint_pretrends_F returns NA with no pre-periods", {
  res <- joint_pretrends_F(est = numeric(0), se = numeric(0), n_pre = 0)
  expect_true(is.na(res$f.stat))
  expect_true(is.na(res$p.value))
})

test_that("joint_pretrends_F errors when neither se nor vcov is supplied", {
  expect_error(joint_pretrends_F(est = c(0.1, 0.2)))
})

test_that("pretreat_mean averages treated pre-period outcomes only", {
  panel <- data.frame(
    y = c(1, 2, 3, 4, 5, 6),
    g = c(2002, 2002, 2002, Inf, Inf, Inf),
    t = c(2000, 2001, 2002, 2000, 2001, 2002)
  )
  # Treated unit's pre-periods are t < g = {2000, 2001}; never-treated excluded.
  expect_equal(pretreat_mean(panel, "y", "g", "t"), mean(c(1, 2)))
})
