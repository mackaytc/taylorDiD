################################################################################
#
# FILE: taylorDiD-package.R
#
# OVERVIEW: Package-level documentation and namespace imports. Also declares the
# data.table column symbols that are referenced via non-standard evaluation, so
# that R CMD check does not flag them as undefined globals.
#
# AUTHOR: Taylor Mackay (tmackay@fullerton.edu)
#
################################################################################

#' @keywords internal
"_PACKAGE"

#' @import data.table
#' @import ggplot2
#' @importFrom stats pf pnorm qnorm lm coef vcov nobs as.formula
NULL

# Column names used unquoted inside dt[...] (data.table NSE) and ggplot aes().
# "." is data.table's .() helper, which R CMD check otherwise flags.
utils::globalVariables(c(
  ".",
  ".g", ".t", ".id",
  ".g.temp", ".t.temp", ".y.temp", ".id.temp", ".is.pre", ".y.detrended",
  ".first.treat", ".treat.static", "intercept", "slope",
  "event.time", "term", "estimate", "std.error", "p.value",
  "conf.low", "conf.high", "is.treated", "first.treat.time",
  "treat.did2s", "rel.year.did2s",
  "n.areas", "n.not.yet.treated.or.inf", "n.already.treated", "n.post.treated"
))

################################################################################
# End of File
################################################################################
