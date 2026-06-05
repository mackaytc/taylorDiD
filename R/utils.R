################################################################################
#
# FILE: utils.R
#
# OVERVIEW: Small formatting helpers shared by the table and plot builders:
# significance stars (plain and LaTeX) and fixed-decimal number formatting.
#
# AUTHOR: Taylor Mackay (tmackay@fullerton.edu)
#
################################################################################

#' Significance stars from p-values
#'
#' Maps p-values to the house-style significance markers: `***` for p < 0.01,
#' `**` for p < 0.05, `*` for p < 0.10, and `""` otherwise. Vectorized; `NA`
#' p-values return `""`.
#'
#' @param p Numeric vector of p-values.
#' @return Character vector of stars, the same length as `p`.
#' @examples
#' add_stars(c(0.001, 0.03, 0.2, NA))
#' @export
add_stars <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.01, "***",
      ifelse(p < 0.05, "**",
        ifelse(p < 0.10, "*", ""))))
}

# LaTeX significance stars, e.g. "$^{***}$"; "" when not significant.
stars_latex <- function(p) {
  plain <- add_stars(p)
  ifelse(plain == "", "", paste0("$^{", plain, "}$"))
}

# Format a number to a fixed number of decimal places.
fmt_num <- function(x, digits = 4) {
  sprintf(paste0("%.", digits, "f"), x)
}

################################################################################
# End of File
################################################################################
