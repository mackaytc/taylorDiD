################################################################################
#
# FILE: event-study.R
#
# OVERVIEW: The shared result objects that every estimator returns, their print
# methods (which hold the console summary so the estimator bodies stay quiet),
# and post_avg_did() -- the house-rule auxiliary static DiD used to summarize the
# post-treatment average effect.
#
# Two S3 classes:
#   taylorDiD_did -- a static (single-coefficient) estimate.
#   taylorDiD_es  -- an event study: coef-by-coef estimates plus the
#                    post-treatment average, joint pre-trends test, pre-treatment
#                    mean, and N.
#
# AUTHOR: Taylor Mackay (tmackay@fullerton.edu)
#
################################################################################

# Pretty estimator label for printing.
estimator_label <- function(estimator) {
  switch(estimator,
    bjs   = "BJS imputation",
    dcdh  = "de Chaisemartin & D'Haultfoeuille",
    did2s = "Gardner two-stage (did2s)",
    estimator)
}

#' Construct a static DiD result object
#'
#' Bundles a single-coefficient DiD estimate into a `taylorDiD_did` object.
#' Extra estimator-specific fields (e.g. `raw.result`) pass through via `...`.
#'
#' @param estimate,std.error,p.value The point estimate and its inference.
#' @param pre.mean.treated Pre-treatment outcome mean for treated units.
#' @param n.obs,n.groups,n.treated Sample sizes.
#' @param estimator Estimator tag (`"bjs"`, `"dcdh"`, or `"did2s"`).
#' @param ... Additional fields stored on the result.
#' @return An object of class `taylorDiD_did` (a named list).
#' @keywords internal
#' @export
new_did_result <- function(estimate, std.error, p.value, pre.mean.treated,
                           n.obs, n.groups, n.treated, estimator, ...) {
  structure(
    c(list(
      estimate = estimate,
      std.error = std.error,
      p.value = p.value,
      pre.mean.treated = pre.mean.treated,
      n.obs = n.obs,
      n.groups = n.groups,
      n.treated = n.treated,
      estimator = estimator
    ), list(...)),
    class = "taylorDiD_did"
  )
}

#' Construct an event-study result object
#'
#' Bundles the pieces of an event study into a `taylorDiD_es` object consumed by
#' [plot_event_study()] and [tex_event_study_table()]. The field layout matches
#' the lists the original helpers returned, so existing `$` access keeps working;
#' this only adds an S3 class and an `estimator` tag.
#'
#' @param coefficients data.table of coefficient-by-event-time estimates.
#' @param pre.mean.treated Pre-treatment outcome mean for treated units.
#' @param pretrends.test List from [joint_pretrends_F()].
#' @param avg.post.effect List describing the post-treatment average effect.
#' @param n.obs,n.groups,n.treated Sample sizes.
#' @param estimator Estimator tag (`"bjs"`, `"dcdh"`, or `"did2s"`).
#' @param ... Additional fields stored on the result.
#' @return An object of class `taylorDiD_es` (a named list).
#' @keywords internal
#' @export
new_es_result <- function(coefficients, pre.mean.treated, pretrends.test,
                          avg.post.effect, n.obs, n.groups, n.treated,
                          estimator, ...) {
  structure(
    c(list(
      coefficients = coefficients,
      pre.mean.treated = pre.mean.treated,
      pretrends.test = pretrends.test,
      avg.post.effect = avg.post.effect,
      n.obs = n.obs,
      n.groups = n.groups,
      n.treated = n.treated,
      estimator = estimator
    ), list(...)),
    class = "taylorDiD_es"
  )
}

#' @export
print.taylorDiD_did <- function(x, digits = 4, ...) {
  cat(sprintf("\n=== Static DiD Estimate (%s) ===\n", estimator_label(x$estimator)))
  if (!is.null(x$yname)) cat(sprintf("Outcome: %s\n", x$yname))
  if (!is.null(x$treatment.note)) cat(x$treatment.note, "\n", sep = "")
  cat(sprintf("Total groups: %d | Treated groups: %d\n", x$n.groups, x$n.treated))
  cat(sprintf("Estimate: %s (se %s, p %s)\n",
              fmt_num(x$estimate, digits), fmt_num(x$std.error, digits),
              fmt_num(x$p.value, digits)))
  cat(sprintf("Pre-treatment mean (treated): %s\n",
              fmt_num(x$pre.mean.treated, digits)))
  invisible(x)
}

#' @export
print.taylorDiD_es <- function(x, digits = 4, ...) {
  cat(sprintf("\n=== Event Study Estimates (%s) ===\n",
              estimator_label(x$estimator)))
  if (!is.null(x$yname)) cat(sprintf("Outcome: %s\n", x$yname))
  if (!is.null(x$treatment.note)) cat(x$treatment.note, "\n", sep = "")
  cat(sprintf("Total groups: %d | Treated groups: %d\n", x$n.groups, x$n.treated))
  cat(sprintf("Pre-treatment mean (treated): %s\n\n",
              fmt_num(x$pre.mean.treated, digits)))

  cat("Coefficients:\n")
  print(x$coefficients[, .(event.time, estimate, std.error, p.value,
                           conf.low, conf.high)], digits = digits)

  pt <- x$pretrends.test
  cat(sprintf("\nJoint pre-trends test: F = %s (df = %d), p = %s\n",
              fmt_num(pt$f.stat, digits), pt$df, fmt_num(pt$p.value, digits)))

  ap <- x$avg.post.effect
  cat(sprintf("Avg post effect (t = %d to %d): %s (se %s, p %s)\n",
              ap$window[1], ap$window[2], fmt_num(ap$estimate, digits),
              fmt_num(ap$std.error, digits), fmt_num(ap$p.value, digits)))
  if (!is.null(ap$pct.change) && is.finite(ap$pct.change)) {
    cat(sprintf("Percent change vs pre-treatment mean: %.1f%%\n", ap$pct.change))
  }
  invisible(x)
}

#' Post-treatment average effect via an auxiliary static DiD
#'
#' Implements the house rule for the summary post-treatment effect: rather than
#' averaging the event-study coefficients, fit a static DiD on the subset of
#' never-treated, not-yet-treated, and treated-within-window observations
#' (`event_time %in% 0:max(post.window)`). Supported for the imputation (`"bjs"`)
#' and two-stage (`"did2s"`) estimators.
#'
#' @param data Data frame or data.table (already detrended if applicable).
#' @param yname,gname,tname,idname Column names (strings).
#' @param post.window Integer vector of post-treatment event times; only
#'   `max(post.window)` is used to bound the window.
#' @param cluster.var Clustering variable (default: `idname`).
#' @param pre.mean.treated Pre-treatment mean, used to express the percent change.
#' @param estimator One of `"bjs"` or `"did2s"`.
#' @return List with `estimate`, `std.error`, `p.value`, `pct.change`, `window`.
#' @keywords internal
#' @export
post_avg_did <- function(data, yname, gname, tname, idname, post.window,
                         cluster.var = idname, pre.mean.treated = NA_real_,
                         estimator = c("bjs", "did2s")) {
  estimator <- match.arg(estimator)
  dt <- as.data.table(data)

  g <- dt[[gname]]
  evt <- dt[[tname]] - g
  max.post <- max(post.window)

  # Subset: never-treated OR not-yet-treated OR treated within the post window.
  keep <- is_never_treated(g) | dt[[tname]] < g |
    (!is.na(evt) & evt >= 0 & evt <= max.post)
  sub <- dt[keep]

  if (estimator == "bjs") {
    fit <- didimputation::did_imputation(
      data = sub, yname = yname, gname = gname, tname = tname,
      idname = idname, cluster_var = cluster.var
    )
    est <- fit$estimate
    se <- fit$std.error
  } else {
    # did2s static DiD: one absorbing treated-in-window indicator.
    sub[, .treat.static := as.integer(
      !is_never_treated(get(gname)) &
        (get(tname) - get(gname)) %in% 0:max.post
    )]
    first.stage <- as.formula(paste0("~ 0 | ", idname, " + ", tname))
    fit <- did2s::did2s(
      data = sub, yname = yname, first_stage = first.stage,
      second_stage = ~ .treat.static, treatment = ".treat.static",
      cluster_var = idname, verbose = FALSE
    )
    est <- unname(coef(fit)[1])
    se <- sqrt(diag(vcov(fit)))[1]
  }

  list(
    estimate = est,
    std.error = se,
    p.value = 2 * pnorm(-abs(est / se)),
    pct.change = est / pre.mean.treated * 100,
    window = c(0, max.post)
  )
}

################################################################################
# End of File
################################################################################
