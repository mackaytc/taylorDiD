# did2s (Gardner two-stage) tests, on par with the BJS and dCDH suites.

panel <- readRDS(system.file("extdata", "synthetic_panel.rds",
                             package = "taylorDiD"))

test_that("prepare.did2s.data generalizes columns and builds event dummies", {
  prepared <- prepare.did2s.data(
    data.table::as.data.table(panel), "y", "g", "year", "id",
    pre.window = -4:-1, post.window = 0:3)

  expect_identical(prepared$status, "ok")
  expect_true(all(c("treat.did2s", "rel.year.did2s") %in% names(prepared$data)))
  expect_false(-1L %in% prepared$event.map$event.time)   # reference omitted
  expect_true(all(prepared$event.map$term %in% names(prepared$data)))
})

test_that("did2s_event_study returns the shared event-study object", {
  es <- did2s_event_study(panel, "y", "g", "year", "id",
                          pre.window = -4:-1, post.window = 0:3, verbose = FALSE)
  expect_s3_class(es, "taylorDiD_es")
  expect_identical(es$estimator, "did2s")
  expect_true(all(c("coefficients", "pretrends.test", "avg.post.effect",
                    "pre.mean.treated", "n.obs") %in% names(es)))
  # 3 pre (-4, -3, -2; -1 is the reference) + 4 post (0:3) = 7 coefficients.
  expect_equal(nrow(es$coefficients), 7)
  expect_equal(es$pretrends.test$df, 3)
  expect_equal(es$avg.post.effect$window, c(0, 3))
  expect_true(is.finite(es$pretrends.test$p.value))
  expect_true(is.finite(es$avg.post.effect$estimate))

  post <- es$coefficients$estimate[es$coefficients$event.time >= 0]
  expect_gt(mean(post), 0)
})

test_that("did2s_event_study matches an independent did2s::did2s reconstruction", {
  es <- did2s_event_study(panel, "y", "g", "year", "id",
                          pre.window = -4:-1, post.window = 0:3, verbose = FALSE)

  # Rebuild the estimation sample and event dummies from scratch, then call
  # did2s::did2s() directly -- the wrapper must reproduce this exactly.
  d <- data.table::copy(data.table::as.data.table(panel))
  d[, rel := ifelse(is.finite(g), year - g, NA_integer_)]
  d[, treat.did2s := as.integer(is.finite(g) & year >= g)]
  keep.events <- sort(unique(c(-4:-1, 0:3)))
  d <- d[is.infinite(g) | (!is.na(rel) & rel %in% keep.events)]
  events <- setdiff(keep.events, -1L)
  term.of <- function(e) if (e < 0) paste0("event.m", abs(e)) else
    paste0("event.p", e)
  for (e in events) d[, (term.of(e)) := as.integer(!is.na(rel) & rel == e)]
  terms <- vapply(events, term.of, character(1))

  m <- did2s::did2s(d, yname = "y", first_stage = ~ 0 | id + year,
    second_stage = stats::as.formula(paste("~", paste(terms, collapse = " + "))),
    treatment = "treat.did2s", cluster_var = "id", verbose = FALSE)
  ct <- as.data.frame(fixest::coeftable(m))
  ref <- data.table::data.table(term = rownames(ct), est = ct[, 1], se = ct[, 2])
  got <- merge(es$coefficients[, .(term, estimate, std.error)], ref, by = "term")

  expect_equal(got$estimate, got$est, tolerance = 1e-8)
  expect_equal(got$std.error, got$se, tolerance = 1e-8)
  expect_equal(nrow(got), length(events))
})

test_that("did2s_event_study is invariant to column names", {
  es <- did2s_event_study(panel, "y", "g", "year", "id",
                          pre.window = -4:-1, post.window = 0:3, verbose = FALSE)

  p2 <- data.table::copy(data.table::as.data.table(panel))
  data.table::setnames(p2, c("id", "year", "g"), c("area_uid", "yr", "cohort"))
  es2 <- did2s_event_study(p2, "y", "cohort", "yr", "area_uid",
                           pre.window = -4:-1, post.window = 0:3, verbose = FALSE)

  expect_equal(es$coefficients$estimate, es2$coefficients$estimate)
  expect_equal(es$pretrends.test$p.value, es2$pretrends.test$p.value)
})

test_that("did2s_event_study passes cluster.var to the post-average", {
  pnl <- data.table::copy(data.table::as.data.table(panel))
  pnl[, grp := id %% 5L]   # coarser clustering than the unit id

  es_id <- did2s_event_study(pnl, "y", "g", "year", "id",
                             pre.window = -4:-1, post.window = 0:3,
                             cluster.var = "id", verbose = FALSE)
  es_grp <- did2s_event_study(pnl, "y", "g", "year", "id",
                              pre.window = -4:-1, post.window = 0:3,
                              cluster.var = "grp", verbose = FALSE)

  # Different clustering changes the post-average standard error.
  expect_false(isTRUE(all.equal(es_id$avg.post.effect$std.error,
                                es_grp$avg.post.effect$std.error)))
})

test_that("did2s preparation fails cleanly without sufficient support", {
  # All units adopt in the same late period with no never-treated units, so the
  # support window leaves no post-treatment observations.
  pnl <- data.table::CJ(id = 1:30, year = 2000:2010)
  pnl[, g := 2009]
  pnl[, y := rnorm(.N)]

  prepared <- prepare.did2s.data(pnl, "y", "g", "year", "id",
                                 pre.window = -4:-1, post.window = 0:3)
  expect_match(prepared$status, "^failed")
  expect_error(
    did2s_event_study(pnl, "y", "g", "year", "id",
                      pre.window = -4:-1, post.window = 0:3, verbose = FALSE),
    "preparation failed"
  )
})
