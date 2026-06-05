################################################################################
#
# FILE: bjs.R
#
# OVERVIEW: BJS imputation estimator (Borusyak, Jaravel, and Spiess) wrappers
# around didimputation::did_imputation(). Provides a static ATT (did.simple) and
# an event study (did.event.study) with the optional unit-specific linear-trend
# pre-residualization. Both delegate inference to the shared helpers in
# inference.R / event-study.R.
#
# AUTHOR: Taylor Mackay (tmackay@fullerton.edu)
#
################################################################################

#' Static DiD via the BJS imputation estimator
#'
#' @param data Data frame or data.table.
#' @param yname Outcome variable name (string).
#' @param gname Treatment-timing variable name; never-treated coded as
#'   `NA`/`Inf`/non-positive (see [is_never_treated()]).
#' @param tname Time-period variable name.
#' @param idname Unit identifier.
#' @param cluster.var Clustering variable (default: `idname`).
#' @param verbose If `TRUE` (default), print a console summary.
#' @return A `taylorDiD_did` result object (list) with `estimate`, `std.error`,
#'   `p.value`, `pre.mean.treated`, `n.obs`, `n.groups`, `n.treated`.
#' @examples
#' panel <- readRDS(system.file("extdata", "synthetic_panel.rds",
#'                              package = "taylorDiD"))
#' \donttest{
#' did.simple(panel, "y", "g", "year", "id")
#' }
#' @export
did.simple <- function(data, yname, gname, tname, idname, cluster.var = NULL,
                       verbose = TRUE) {
  if (is.null(cluster.var)) cluster.var <- idname
  dt <- as.data.table(data)

  unit.treated <- dt[, .(is.treated = any(!is_never_treated(get(gname)))),
                     by = c(idname)]
  n.groups <- nrow(unit.treated)
  n.treated <- sum(unit.treated$is.treated)

  fit <- didimputation::did_imputation(
    data = dt, yname = yname, gname = gname, tname = tname,
    idname = idname, cluster_var = cluster.var
  )

  est <- fit$estimate
  se <- fit$std.error
  pre.mean <- pretreat_mean(dt, yname, gname, tname)

  result <- new_did_result(
    estimate = est, std.error = se, p.value = 2 * pnorm(-abs(est / se)),
    pre.mean.treated = pre.mean, n.obs = nrow(dt), n.groups = n.groups,
    n.treated = n.treated, estimator = "bjs", yname = yname
  )
  if (verbose) print(result)
  invisible(result)
}

#' Event study via the BJS imputation estimator
#'
#' Estimates coefficient-by-event-time effects, a joint pre-trends F-test, and a
#' post-treatment average effect (auxiliary static DiD on the never-/not-yet-/
#' treated-in-window subset, per house style). With `trends.lin = TRUE` the
#' outcome is first pre-residualized on unit-specific linear time trends fit from
#' pre-treatment observations.
#'
#' @inheritParams did.simple
#' @param pre.window Pre-treatment event times (e.g. `-4:-1`).
#' @param post.window Post-treatment event times (e.g. `0:3`).
#' @param trends.lin If `TRUE`, pre-residualize the outcome on unit-specific
#'   pre-treatment linear trends before estimation (default `FALSE`).
#' @return A `taylorDiD_es` result object (list) with `coefficients`,
#'   `pre.mean.treated`, `pretrends.test`, `avg.post.effect`, `n.obs`,
#'   `n.groups`, `n.treated`.
#' @examples
#' panel <- readRDS(system.file("extdata", "synthetic_panel.rds",
#'                              package = "taylorDiD"))
#' \donttest{
#' es <- did.event.study(panel, "y", "g", "year", "id",
#'                       pre.window = -4:-1, post.window = 0:3)
#' es$avg.post.effect
#' }
#' @export
did.event.study <- function(data, yname, gname, tname, idname,
                            pre.window, post.window, cluster.var = NULL,
                            trends.lin = FALSE, verbose = TRUE) {
  if (is.null(cluster.var)) cluster.var <- idname
  dt <- as.data.table(copy(data))

  # Optionally pre-residualize the outcome on unit-specific pre-treatment linear
  # trends -- a robustness check that absorbs differential trends.
  if (trends.lin) {
    dt[, `:=`(.g.temp = get(gname), .t.temp = get(tname),
              .y.temp = get(yname), .id.temp = get(idname))]
    dt[, .is.pre := is_never_treated(.g.temp) | .t.temp < .g.temp]

    # Fit y = a + b * t on each unit's pre-treatment observations.
    trend.fits <- dt[.is.pre == TRUE, .(
      intercept = if (.N >= 2) coef(lm(.y.temp ~ .t.temp))[1]
                  else mean(.y.temp, na.rm = TRUE),
      slope     = if (.N >= 2) coef(lm(.y.temp ~ .t.temp))[2] else 0
    ), by = .id.temp]

    dt <- merge(dt, trend.fits, by = ".id.temp", all.x = TRUE)
    dt[is.na(intercept), intercept := mean(dt$.y.temp, na.rm = TRUE)]
    dt[is.na(slope), slope := 0]
    dt[, .y.detrended := .y.temp - (intercept + slope * .t.temp)]
    dt[[yname]] <- dt$.y.detrended

    dt[, c(".g.temp", ".t.temp", ".y.temp", ".id.temp", ".is.pre",
           "intercept", "slope", ".y.detrended") := NULL]
  }

  unit.treated <- dt[, .(is.treated = any(!is_never_treated(get(gname)))),
                     by = c(idname)]
  n.groups <- nrow(unit.treated)
  n.treated <- sum(unit.treated$is.treated)

  raw <- didimputation::did_imputation(
    data = dt, yname = yname, gname = gname, tname = tname, idname = idname,
    horizon = post.window, pretrends = pre.window, cluster_var = cluster.var
  )

  coefs <- as.data.table(raw)
  coefs[, event.time := as.numeric(term)]
  coefs[, p.value := 2 * pnorm(-abs(estimate / std.error))]
  setorder(coefs, event.time)

  pre.mean <- pretreat_mean(dt, yname, gname, tname)

  pre.coefs <- coefs[event.time < 0]
  pretrends <- joint_pretrends_F(est = pre.coefs$estimate,
                                 se = pre.coefs$std.error,
                                 n_pre = nrow(pre.coefs))

  avg.post <- post_avg_did(dt, yname, gname, tname, idname,
                           post.window = post.window, cluster.var = cluster.var,
                           pre.mean.treated = pre.mean, estimator = "bjs")

  result <- new_es_result(
    coefficients = coefs, pre.mean.treated = pre.mean,
    pretrends.test = pretrends, avg.post.effect = avg.post,
    n.obs = nrow(dt), n.groups = n.groups, n.treated = n.treated,
    estimator = "bjs", yname = yname, trends.lin = trends.lin
  )
  if (verbose) print(result)
  invisible(result)
}

################################################################################
# End of File
################################################################################
