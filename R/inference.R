################################################################################
#
# FILE: inference.R
#
# OVERVIEW: Shared inference building blocks used by every estimator:
#   - never-treated coding (one place, since packages disagree on the sentinel)
#   - the joint pre-trends F-test (Wald / n_pre, p-value from an F with df2 = Inf)
#   - the pre-treatment outcome mean for treated units
# Coding never-treated wrong fails silently downstream, so it lives here and is
# unit-tested.
#
# AUTHOR: Taylor Mackay (tmackay@fullerton.edu)
#
################################################################################

#' Identify never-treated units from a treatment-timing variable
#'
#' Treats `NA`, `Inf`, and non-positive timing values as never-treated. Packages
#' disagree on the never-treated sentinel (`didimputation` accepts 0/NA/Inf, the
#' did2s helpers use `Inf`), so this predicate is the single source of truth.
#'
#' @param g Numeric vector of treatment-timing (cohort) values.
#' @return Logical vector, `TRUE` where `g` denotes a never-treated unit.
#' @examples
#' is_never_treated(c(2005, NA, Inf, 0, -1, 2010))
#' @export
is_never_treated <- function(g) {
  is.na(g) | is.infinite(g) | g <= 0
}

#' Normalize never-treated coding to a chosen sentinel
#'
#' Recodes never-treated values (see [is_never_treated()]) to the sentinel a
#' given package expects, leaving genuine treatment-timing values untouched.
#'
#' @param g Numeric vector of treatment-timing (cohort) values.
#' @param target One of `"zero"`, `"NA"`, or `"Inf"` -- the sentinel to assign
#'   never-treated units.
#' @return Numeric vector with never-treated units recoded to `target`.
#' @examples
#' code_never_treated(c(2005, NA, Inf, 0), target = "Inf")
#' code_never_treated(c(2005, NA, Inf, 0), target = "zero")
#' @export
code_never_treated <- function(g, target = c("zero", "NA", "Inf")) {
  target <- match.arg(target)
  sentinel <- switch(target, zero = 0, "NA" = NA_real_, "Inf" = Inf)
  g <- as.numeric(g)
  g[is_never_treated(g)] <- sentinel
  g
}

#' Joint pre-trends F-test
#'
#' Reports the F-statistic version of the joint test that all pre-treatment
#' coefficients are zero: the Wald statistic divided by the number of
#' pre-periods, with a p-value from `pf(F, df1 = n_pre, df2 = Inf)`. When a full
#' covariance matrix `vcov` is supplied the Wald statistic uses it; otherwise it
#' falls back to the diagonal approximation `sum((est / se)^2)` (which is what
#' the imputation and DIDmultiplegtDYN outputs support, since they expose only
#' point estimates and standard errors).
#'
#' @param est Numeric vector of pre-treatment coefficient estimates.
#' @param se Numeric vector of standard errors (diagonal approximation). Ignored
#'   when `vcov` is supplied.
#' @param vcov Optional full covariance matrix of the pre-treatment estimates.
#' @param n_pre Number of pre-treatment periods (defaults to `length(est)`).
#' @return List with `f.stat`, `df`, and `p.value`.
#' @examples
#' joint_pretrends_F(est = c(0.1, -0.05, 0.02), se = c(0.08, 0.07, 0.09))
#' @export
joint_pretrends_F <- function(est, se = NULL, vcov = NULL,
                              n_pre = length(est)) {
  if (n_pre <= 0 || length(est) == 0) {
    return(list(f.stat = NA_real_, df = n_pre, p.value = NA_real_))
  }

  wald <- if (!is.null(vcov)) {
    as.numeric(t(est) %*% qr.solve(vcov, est))
  } else {
    if (is.null(se)) stop("Supply either `se` or `vcov`.")
    sum((est / se)^2)
  }

  f.stat <- wald / n_pre
  list(
    f.stat = f.stat,
    df = n_pre,
    p.value = pf(f.stat, df1 = n_pre, df2 = Inf, lower.tail = FALSE)
  )
}

#' Pre-treatment outcome mean for treated units (cohort timing)
#'
#' Averages the outcome over treated units' pre-treatment observations, where a
#' unit is treated when its timing value is a genuine cohort (see
#' [is_never_treated()]) and an observation is pre-treatment when `t < g`. This
#' is the convention used by the BJS and did2s estimators, which carry cohort
#' timing in `gname`.
#'
#' @param data Data frame or data.table.
#' @param yname Outcome variable name (string).
#' @param gname Treatment-timing (cohort) variable name (string).
#' @param tname Time-period variable name (string).
#' @return Numeric scalar, the pre-treatment outcome mean for treated units.
#' @examples
#' panel <- data.frame(
#'   y = c(1, 2, 3, 4, 5, 6),
#'   g = c(2002, 2002, 2002, Inf, Inf, Inf),
#'   t = c(2000, 2001, 2002, 2000, 2001, 2002)
#' )
#' pretreat_mean(panel, "y", "g", "t")
#' @export
pretreat_mean <- function(data, yname, gname, tname) {
  dt <- as.data.table(data)
  g <- dt[[gname]]
  treated <- !is_never_treated(g)
  pre <- dt[[tname]] < g
  mean(dt[[yname]][treated & pre], na.rm = TRUE)
}

# Per-unit first treatment change for the dCDH (level-treatment) case: the
# earliest period a unit's treatment differs from its first-period level, or NA
# if it never changes. dCDH identifies effects from these switchers, so a unit
# with a constant treatment (including a constant nonzero level) is not treated.
dcdh_first_change <- function(data, dname, tname, idname) {
  dt <- as.data.table(data)
  dt[, {
    o <- order(get(tname))
    d <- get(dname)[o]
    tt <- get(tname)[o]
    changed <- which(d != d[1])
    list(first.change = if (length(changed)) as.numeric(tt[changed[1]])
                        else NA_real_)
  }, by = c(idname)]
}

# Pre-treatment outcome mean for the dCDH case: the mean outcome over switchers'
# observations that precede their first treatment change.
pretreat_mean_dyn <- function(data, yname, dname, tname, idname) {
  dt <- as.data.table(data)
  changes <- dcdh_first_change(dt, dname, tname, idname)
  dt <- merge(dt, changes, by = idname, all.x = TRUE)
  pre <- !is.na(dt$first.change) & dt[[tname]] < dt$first.change
  mean(dt[[yname]][pre], na.rm = TRUE)
}

################################################################################
# End of File
################################################################################
