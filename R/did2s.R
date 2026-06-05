################################################################################
#
# FILE: did2s.R
#
# OVERVIEW: Gardner (2021) two-stage estimator via did2s::did2s(). A top-level
# did2s_event_study() prepares the panel, fits the two-stage model, extracts the
# event-study coefficients, and returns the shared result object. Outcome, time,
# unit, and treatment timing enter as the yname / tname / idname / gname
# parameters.
#
# The internal helpers below are not exported.
#
# AUTHOR: Taylor Mackay (tmackay@fullerton.edu)
#
################################################################################

# Event-time -> second-stage term name, e.g. -2 -> "event.m2", 3 -> "event.p3".
# @noRd
event.term.name <- function(event.time) {
  if (event.time < 0) paste0("event.m", abs(event.time))
  else paste0("event.p", event.time)
}

# Build the did2s estimation panel: drop units treated too early to have full
# pre-window support, find the not-yet-treated support window, truncate to it
# (recoding units treated after the window as never-treated), then add the
# treatment indicator, relative-event-time, and one 0/1 dummy per event time.
# Works on internal canonical columns .g/.t/.id derived from gname/tname/idname;
# never-treated units are identified by is_never_treated().
# @noRd
prepare.did2s.data <- function(panel.dt, yname, gname, tname, idname,
                               pre.window, post.window) {
  dt <- as.data.table(copy(panel.dt))
  dt[, `:=`(.g = as.numeric(get(gname)), .t = get(tname), .id = get(idname))]

  earliest.allowed.g <- min(dt$.t) + abs(min(pre.window))

  early.units <- dt[!is_never_treated(.g) & .g < earliest.allowed.g, unique(.id)]
  if (length(early.units) > 0) dt <- dt[!.id %in% early.units]

  support.before <- dt[, .(
    n.areas = uniqueN(.id),
    n.not.yet.treated.or.inf = uniqueN(.id[is_never_treated(.g) | .g > .t]),
    n.already.treated = uniqueN(.id[!is_never_treated(.g) & .g <= .t])
  ), by = .t][order(.t)]

  support.years <- support.before[n.not.yet.treated.or.inf > 0, .t]
  if (length(support.years) == 0) {
    return(list(
      data = data.table(), support = support.before, event.map = data.table(),
      earliest.allowed.g = earliest.allowed.g,
      n.early.units.excluded = length(early.units),
      latest.support.year = NA_integer_,
      n.after.truncation.units.set.inf = 0L,
      status = "failed_no_not_yet_treated_support"
    ))
  }

  latest.support.year <- max(support.years)
  dt <- dt[.t <= latest.support.year]
  after.truncation.units <- dt[!is_never_treated(.g) &
                                 .g > latest.support.year, unique(.id)]
  if (length(after.truncation.units) > 0) dt[.id %in% after.truncation.units,
                                             .g := Inf]

  dt[, treat.did2s := as.integer(!is_never_treated(.g) & .t >= .g)]
  dt[, rel.year.did2s := NA_integer_]
  dt[!is_never_treated(.g), rel.year.did2s := as.integer(.t - .g)]

  keep.events <- sort(unique(c(pre.window, post.window)))
  dt <- dt[is_never_treated(.g) |
             (!is.na(rel.year.did2s) & rel.year.did2s %in% keep.events)]
  dt <- dt[!is.na(get(yname))]

  # Event dummies for every kept event time except the omitted reference (-1).
  event.map <- data.table(event.time = setdiff(keep.events, -1L))
  event.map[, term := vapply(event.time, event.term.name, character(1))]
  for (k in seq_len(nrow(event.map))) {
    dt[, (event.map$term[k]) := as.integer(!is.na(rel.year.did2s) &
          rel.year.did2s == event.map$event.time[k])]
  }

  support.after <- dt[, .(
    n.areas = uniqueN(.id),
    n.not.yet.treated.or.inf = uniqueN(.id[is_never_treated(.g) | .g > .t]),
    n.already.treated = uniqueN(.id[!is_never_treated(.g) & .g <= .t]),
    n.post.treated = uniqueN(.id[treat.did2s == 1])
  ), by = .t][order(.t)]

  setnames(support.before, ".t", tname)
  setnames(support.after, ".t", tname)
  dt[, c(".g", ".t", ".id") := NULL]

  # did2s needs both treated and control observations and at least one active
  # event dummy; without them the second stage has no usable variation.
  has.controls <- any(dt$treat.did2s == 0)
  has.treated <- any(dt$treat.did2s == 1)
  dummy.active <- any(vapply(event.map$term,
                             function(col) sum(dt[[col]]) > 0, logical(1)))
  status <- if (nrow(dt) > 0 && has.controls && has.treated && dummy.active) {
    "ok"
  } else {
    "failed_insufficient_support"
  }

  list(
    data = dt, support = support.after, event.map = event.map,
    earliest.allowed.g = earliest.allowed.g,
    n.early.units.excluded = length(early.units),
    latest.support.year = latest.support.year,
    n.after.truncation.units.set.inf = length(after.truncation.units),
    status = status
  )
}

# Tidy the fixest second-stage coefficient table into the shared coefficient
# layout, matching rows to event times via event.map.
# @noRd
extract.did2s.coefficients <- function(model, event.map) {
  ct <- as.data.table(fixest::coeftable(model), keep.rownames = "term")
  est.col <- names(ct)[grepl("^Estimate$", names(ct), ignore.case = TRUE)][1]
  se.col  <- names(ct)[grepl("Std\\.? Error", names(ct), ignore.case = TRUE)][1]
  p.col   <- names(ct)[grepl("Pr\\(", names(ct))][1]
  if (is.na(p.col)) p.col <- names(ct)[grepl("p", names(ct),
                                             ignore.case = TRUE)][1]

  out <- merge(event.map, ct, by = "term", all.x = TRUE)
  out[, `:=`(
    estimate = get(est.col),
    std.error = get(se.col),
    p.value = if (!is.na(p.col)) get(p.col)
              else 2 * pnorm(-abs(get(est.col) / get(se.col)))
  )]
  out[, conf.low := estimate - qnorm(0.975) * std.error]
  out[, conf.high := estimate + qnorm(0.975) * std.error]
  out[, .(term, event.time, estimate, std.error, p.value, conf.low, conf.high)]
}

#' Event study via the Gardner two-stage estimator
#'
#' Runs the full did2s flow: prepare the estimation panel, fit `did2s::did2s()`
#' with unit and time fixed effects in the first stage and event-time dummies in
#' the second, extract the coefficients, and summarize with the house-style joint
#' pre-trends F-test and the auxiliary-static-DiD post-treatment average.
#'
#' @details
#' The joint pre-trends test uses the full coefficient covariance from the fixest
#' fit (via [joint_pretrends_F()]); with `df2 = Inf` this F p-value equals the
#' chi-square Wald p-value. The post-treatment average is the auxiliary static
#' DiD on the never-/not-yet-/treated-in-window subset (see [post_avg_did()]), as
#' in the BJS estimator. Never-treated units are coded via [is_never_treated()]
#' (`NA`/`Inf`/non-positive). The omitted event-study reference is event time -1.
#' The estimation sample keeps never-treated units and treated observations whose
#' event time lies in `pre.window` or `post.window`; treated observations outside
#' that window are dropped so they do not enter the reference category.
#'
#' @param data Data frame or data.table.
#' @param yname Outcome variable name (string).
#' @param gname Treatment-timing (cohort) variable; never-treated coded as
#'   `NA`/`Inf`/non-positive.
#' @param tname Time-period variable name.
#' @param idname Unit identifier.
#' @param pre.window Pre-treatment event times (e.g. `-4:-1`).
#' @param post.window Post-treatment event times (e.g. `0:3`).
#' @param cluster.var Clustering variable (default: `idname`).
#' @param verbose If `TRUE` (default), print a console summary.
#' @return A `taylorDiD_es` result object; also carries the `prepared` panel
#'   list and the underlying fixest model in `raw.result`.
#'
#' @examples
#' panel <- readRDS(system.file("extdata", "synthetic_panel.rds",
#'                              package = "taylorDiD"))
#' \donttest{
#' es <- did2s_event_study(panel, "y", "g", "year", "id",
#'                         pre.window = -4:-1, post.window = 0:3)
#' es$coefficients
#' }
#' @export
did2s_event_study <- function(data, yname, gname, tname, idname,
                              pre.window, post.window, cluster.var = NULL,
                              verbose = TRUE) {
  if (is.null(cluster.var)) cluster.var <- idname
  dt <- as.data.table(data)

  prepared <- prepare.did2s.data(dt, yname, gname, tname, idname,
                                 pre.window, post.window)
  if (!identical(prepared$status, "ok")) {
    stop("did2s data preparation failed: ", prepared$status)
  }

  event.terms <- prepared$event.map$term
  first.stage <- as.formula(paste0("~ 0 | ", idname, " + ", tname))
  second.stage <- as.formula(paste0("~ ", paste(event.terms, collapse = " + ")))

  model <- did2s::did2s(
    data = prepared$data, yname = yname, first_stage = first.stage,
    second_stage = second.stage, treatment = "treat.did2s",
    cluster_var = cluster.var, verbose = FALSE
  )

  coefs <- extract.did2s.coefficients(model, prepared$event.map)
  setorder(coefs, event.time)

  # Treated units' pre-treatment observations carry rel.year.did2s < 0.
  pre.mean <- mean(
    prepared$data[!is.na(rel.year.did2s) & rel.year.did2s < 0][[yname]],
    na.rm = TRUE
  )

  pre.terms <- prepared$event.map[event.time < 0, term]
  pre.idx <- match(pre.terms, names(coef(model)))
  pre.idx <- pre.idx[!is.na(pre.idx)]
  pretrends <- joint_pretrends_F(
    est = coef(model)[pre.idx],
    vcov = vcov(model)[pre.idx, pre.idx, drop = FALSE],
    n_pre = length(pre.idx)
  )

  avg.post <- post_avg_did(dt, yname, gname, tname, idname,
                           post.window = post.window, cluster.var = cluster.var,
                           pre.mean.treated = pre.mean, estimator = "did2s")

  n.groups <- uniqueN(dt[[idname]])
  n.treated <- uniqueN(dt[!is_never_treated(get(gname)), get(idname)])

  result <- new_es_result(
    coefficients = coefs, pre.mean.treated = pre.mean,
    pretrends.test = pretrends, avg.post.effect = avg.post,
    n.obs = nrow(prepared$data), n.groups = n.groups, n.treated = n.treated,
    estimator = "did2s", yname = yname, prepared = prepared, raw.result = model
  )
  if (verbose) print(result)
  invisible(result)
}

################################################################################
# End of File
################################################################################
