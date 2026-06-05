################################################################################
#
# FILE: dcdh.R
#
# OVERVIEW: de Chaisemartin and D'Haultfoeuille estimator wrappers around
# DIDmultiplegtDYN::did_multiplegt_dyn(). Provides a static effect
# (did.simple.dyn) and an event study (did.event.study.dyn) that handle binary,
# discrete, and continuous treatments. Treatment enters as a level (dname), not
# a cohort, so switchers are units whose treatment changes from its first-period
# level, and pre-treatment periods precede that first change.
#
# NOTE: the post-treatment average here is the mean of the dynamic effects (with
# a conservative pooled SE), because the dCDH estimator has no clean static-DiD
# analog; its natural aggregate is the mean of the period-specific effects.
#
# AUTHOR: Taylor Mackay (tmackay@fullerton.edu)
#
################################################################################

# DIDmultiplegtDYN evaluates polars expressions that reference the `pl` object on
# the search path, but it does not import polars itself -- so attach polars on
# demand (with a clear message if it is missing).
ensure_polars <- function() {
  if (!requireNamespace("polars", quietly = TRUE)) {
    stop("The 'polars' package is required by DIDmultiplegtDYN. Install it with:\n",
         "  install.packages('polars', repos = 'https://rpolars.r-universe.dev')",
         call. = FALSE)
  }
  if (!"package:polars" %in% search()) attachNamespace("polars")
  invisible(TRUE)
}

# Human-readable description of the treatment specification, for printing.
dcdh_treatment_note <- function(dname, continuous, normalized, trends.lin) {
  treat.type <- if (normalized) {
    " (normalized: effect per unit)"
  } else if (!is.null(continuous)) {
    sprintf(" (continuous, order %d)", continuous)
  } else {
    " (binary)"
  }
  note <- sprintf("Treatment: %s%s", dname, treat.type)
  if (trends.lin) note <- paste0(note, "\nGroup-specific linear trends: YES")
  note
}

# Pull the Placebos (pre) and Effects (post) matrices out of a did_multiplegt_dyn
# result into one tidy coefficient table. Effect_1 is event time 0, and
# Placebo_i is event time -i.
extract_dcdh_coefficients <- function(raw) {
  rows <- list()

  placebos <- raw$results$Placebos
  if (!is.null(placebos) && nrow(placebos) > 0) {
    for (i in seq_len(nrow(placebos))) {
      rows[[length(rows) + 1]] <- data.table(
        term = as.character(-i), event.time = -i,
        estimate = placebos[i, "Estimate"], std.error = placebos[i, "SE"],
        conf.low = placebos[i, "LB CI"], conf.high = placebos[i, "UB CI"]
      )
    }
  }

  effects <- raw$results$Effects
  if (!is.null(effects) && nrow(effects) > 0) {
    for (i in seq_len(nrow(effects))) {
      rows[[length(rows) + 1]] <- data.table(
        term = as.character(i - 1), event.time = i - 1,
        estimate = effects[i, "Estimate"], std.error = effects[i, "SE"],
        conf.low = effects[i, "LB CI"], conf.high = effects[i, "UB CI"]
      )
    }
  }

  coefs <- rbindlist(rows)
  setorder(coefs, event.time)
  coefs[]
}

# dCDH post-treatment average: mean of the dynamic effects, with a conservative
# pooled standard error (root-mean-square of the period SEs).
dcdh_post_avg <- function(post.coefs, pre.mean.treated, n.post) {
  if (nrow(post.coefs) == 0) {
    return(list(estimate = NA_real_, std.error = NA_real_, p.value = NA_real_,
                pct.change = NA_real_, window = c(0, n.post - 1)))
  }
  est <- mean(post.coefs$estimate)
  se <- sqrt(mean(post.coefs$std.error^2))
  list(
    estimate = est, std.error = se, p.value = 2 * pnorm(-abs(est / se)),
    pct.change = est / pre.mean.treated * 100, window = c(0, n.post - 1)
  )
}

#' Static DiD via the de Chaisemartin & D'Haultfoeuille estimator
#'
#' @param data Data frame or data.table.
#' @param yname Outcome variable name (string).
#' @param dname Treatment variable name (binary, discrete, or continuous level).
#' @param tname Time-period variable name.
#' @param idname Unit identifier.
#' @param cluster.var Clustering variable (default: `idname`).
#' @param continuous Polynomial order if treatment is continuous (`NULL` for
#'   binary/discrete).
#' @param normalized If `TRUE`, estimate the effect of a one-unit treatment
#'   increase.
#' @param trends.lin If `TRUE`, include group-specific linear time trends.
#' @param verbose If `TRUE` (default), print a console summary.
#' @return A `taylorDiD_did` result object; also carries `continuous`,
#'   `normalized`, `trends.lin`, and the underlying `raw.result`.
#' @examples
#' panel <- readRDS(system.file("extdata", "synthetic_panel.rds",
#'                              package = "taylorDiD"))
#' \donttest{
#' did.simple.dyn(panel, "y", "d", "year", "id")
#' }
#' @export
did.simple.dyn <- function(data, yname, dname, tname, idname,
                           cluster.var = NULL, continuous = NULL,
                           normalized = FALSE, trends.lin = FALSE,
                           verbose = TRUE) {
  if (is.null(cluster.var)) cluster.var <- idname
  ensure_polars()
  dt <- as.data.table(data)

  changes <- dcdh_first_change(dt, dname, tname, idname)
  n.groups <- nrow(changes)
  n.treated <- sum(!is.na(changes$first.change))

  raw <- DIDmultiplegtDYN::did_multiplegt_dyn(
    df = dt, outcome = yname, group = idname, time = tname, treatment = dname,
    effects = 1, placebo = 0, cluster = cluster.var, continuous = continuous,
    normalized = normalized, trends_lin = trends.lin, graph_off = TRUE
  )

  effects.mat <- raw$results$Effects
  est <- effects.mat[1, "Estimate"]
  se <- effects.mat[1, "SE"]
  pre.mean <- pretreat_mean_dyn(dt, yname, dname, tname, idname)

  result <- new_did_result(
    estimate = est, std.error = se, p.value = 2 * pnorm(-abs(est / se)),
    pre.mean.treated = pre.mean, n.obs = nrow(dt), n.groups = n.groups,
    n.treated = n.treated, estimator = "dcdh", yname = yname,
    continuous = continuous, normalized = normalized, trends.lin = trends.lin,
    raw.result = raw,
    treatment.note = dcdh_treatment_note(dname, continuous, normalized,
                                         trends.lin)
  )
  if (verbose) print(result)
  invisible(result)
}

#' Event study via the de Chaisemartin & D'Haultfoeuille estimator
#'
#' Estimates dynamic (event-study) effects with `n.pre` placebos and `n.post`
#' effects, a joint pre-trends F-test, and a post-treatment average effect (the
#' mean of the dynamic effects, since the dCDH estimator has no clean static-DiD
#' analog).
#'
#' @inheritParams did.simple.dyn
#' @param n.pre Number of pre-treatment periods (placebos).
#' @param n.post Number of post-treatment periods (dynamic effects).
#' @return A `taylorDiD_es` result object; also carries `continuous`,
#'   `normalized`, `trends.lin`, and the underlying `raw.result`.
#' @examples
#' panel <- readRDS(system.file("extdata", "synthetic_panel.rds",
#'                              package = "taylorDiD"))
#' \donttest{
#' es <- did.event.study.dyn(panel, "y", "d", "year", "id",
#'                           n.pre = 3, n.post = 3)
#' es$coefficients
#' }
#' @export
did.event.study.dyn <- function(data, yname, dname, tname, idname,
                                n.pre, n.post, cluster.var = NULL,
                                continuous = NULL, normalized = FALSE,
                                trends.lin = FALSE, verbose = TRUE) {
  if (is.null(cluster.var)) cluster.var <- idname
  ensure_polars()
  dt <- as.data.table(data)

  changes <- dcdh_first_change(dt, dname, tname, idname)
  n.groups <- nrow(changes)
  n.treated <- sum(!is.na(changes$first.change))

  raw <- DIDmultiplegtDYN::did_multiplegt_dyn(
    df = dt, outcome = yname, group = idname, time = tname, treatment = dname,
    effects = n.post, placebo = n.pre, cluster = cluster.var,
    continuous = continuous, normalized = normalized, trends_lin = trends.lin,
    graph_off = TRUE
  )

  coefs <- extract_dcdh_coefficients(raw)
  coefs[, p.value := 2 * pnorm(-abs(estimate / std.error))]

  pre.mean <- pretreat_mean_dyn(dt, yname, dname, tname, idname)

  pre.coefs <- coefs[event.time < 0]
  pretrends <- joint_pretrends_F(est = pre.coefs$estimate,
                                 se = pre.coefs$std.error,
                                 n_pre = nrow(pre.coefs))

  avg.post <- dcdh_post_avg(coefs[event.time >= 0], pre.mean, n.post)

  result <- new_es_result(
    coefficients = coefs, pre.mean.treated = pre.mean,
    pretrends.test = pretrends, avg.post.effect = avg.post,
    n.obs = nrow(dt), n.groups = n.groups, n.treated = n.treated,
    estimator = "dcdh", yname = yname, continuous = continuous,
    normalized = normalized, trends.lin = trends.lin, raw.result = raw,
    treatment.note = dcdh_treatment_note(dname, continuous, normalized,
                                         trends.lin)
  )
  if (verbose) print(result)
  invisible(result)
}

################################################################################
# End of File
################################################################################
