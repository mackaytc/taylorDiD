# Never-treated coding is the one place a silent mistake breaks everything, so
# its edge cases are tested directly.

test_that("is_never_treated flags NA, Inf, and non-positive timing", {
  g <- c(2005, NA, Inf, 0, -1, 2010)
  expect_equal(is_never_treated(g), c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE))
})

test_that("is_never_treated keeps genuine cohorts", {
  expect_false(any(is_never_treated(c(1990, 2001, 2020))))
})

test_that("code_never_treated recodes only never-treated units", {
  g <- c(2005, NA, Inf, 0, -2, 2010)
  expect_equal(code_never_treated(g, "zero"), c(2005, 0, 0, 0, 0, 2010))
  expect_equal(code_never_treated(g, "Inf"), c(2005, Inf, Inf, Inf, Inf, 2010))

  na_out <- code_never_treated(g, "NA")
  expect_equal(which(is.na(na_out)), c(2L, 3L, 4L, 5L))
  expect_equal(na_out[c(1, 6)], c(2005, 2010))
})

test_that("code_never_treated default target is zero", {
  expect_equal(code_never_treated(c(2005, Inf)), c(2005, 0))
})

test_that("code_never_treated round-trips through is_never_treated", {
  g <- c(2005, NA, Inf, 0, 2010)
  for (tgt in c("zero", "NA", "Inf")) {
    expect_equal(is_never_treated(code_never_treated(g, tgt)),
                 is_never_treated(g))
  }
})

test_that("code_never_treated rejects an unknown target", {
  expect_error(code_never_treated(c(2005, Inf), target = "bogus"))
})
