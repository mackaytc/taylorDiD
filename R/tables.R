################################################################################
#
# FILE: tables.R
#
# OVERVIEW: LaTeX table builders. tex_event_study_table() is the two-panel
# event-study table (Panel A = coefficients by event time; Panel B = post
# average effect, pre-trends test, pre-treatment mean, N), extracted from the
# old per-estimator plot functions so both estimators share one builder.
# demographic.table.output() is a general lm()-to-LaTeX helper.
#
# AUTHOR: Taylor Mackay (tmackay@fullerton.edu)
#
################################################################################

#' Two-panel event-study LaTeX table
#'
#' Builds the house-style table: Panel A reports event-study coefficients (with
#' standard errors and significance stars) by event time across one or more
#' models; Panel B reports the post-treatment average effect, the joint
#' pre-trends test p-value, the pre-treatment outcome mean, and N.
#'
#' @param results A single `taylorDiD_es` object or a list of them (one per
#'   column).
#' @param model.names Optional column titles; defaults to `Model 1`, `Model 2`,
#'   ...
#' @param decimal.places Coefficient/SE precision (default 4).
#' @param effect.label Row label for the Panel B average effect (default
#'   `"Treatment Effect"`).
#' @param output.file Optional path; if given, the LaTeX is written there.
#' @return Character vector of LaTeX lines (invisibly returned semantics: also
#'   returned so it can be `cat()`-ed).
#' @examples
#' es <- new_es_result(
#'   coefficients = data.table::data.table(
#'     event.time = -2:2, term = as.character(-2:2),
#'     estimate = c(0.02, -0.01, 0.30, 0.42, 0.55), std.error = rep(0.1, 5),
#'     p.value = c(0.8, 0.9, 0.01, 0.001, 0.001),
#'     conf.low = c(-0.18, -0.21, 0.10, 0.22, 0.35),
#'     conf.high = c(0.22, 0.19, 0.50, 0.62, 0.75)),
#'   pre.mean.treated = 2.0,
#'   pretrends.test = list(f.stat = 0.3, df = 2, p.value = 0.74),
#'   avg.post.effect = list(estimate = 0.42, std.error = 0.1, p.value = 0.001,
#'                          pct.change = 21, window = c(0, 2)),
#'   n.obs = 500L, n.groups = 50L, n.treated = 20L, estimator = "bjs")
#' cat(tex_event_study_table(es), sep = "\n")
#' @export
tex_event_study_table <- function(results, model.names = NULL,
                                  decimal.places = 4,
                                  effect.label = "Treatment Effect",
                                  output.file = NULL) {
  # Accept a single result or a list of results.
  if (!is.null(results$coefficients)) results <- list(results)
  n.models <- length(results)

  if (is.null(model.names)) {
    model.names <- paste0("Model ", seq_len(n.models))
  } else if (length(model.names) != n.models) {
    stop("Length of model.names must match the number of models.")
  }

  fmt <- function(x) fmt_num(x, decimal.places)
  all.event.times <- sort(unique(unlist(lapply(results, function(r)
    r$coefficients$event.time))))
  n.cols <- 1 + n.models
  col.spec <- paste0("l", paste(rep("c", n.models), collapse = ""))

  tex <- c(
    sprintf("\\begin{tabular}{%s}", col.spec),
    "\\hline\\hline",
    paste0(paste0(" & (", seq_len(n.models), ")", collapse = ""), " \\\\"),
    paste0(" & ", paste(model.names, collapse = " & "), " \\\\"),
    "\\hline",
    sprintf("\\multicolumn{%d}{l}{\\textbf{Panel A: Event Study Estimates}} \\\\",
            n.cols),
    "\\hline"
  )

  # Panel A: one coefficient row and one SE row per event time.
  for (et in all.event.times) {
    coef.vals <- vapply(results, function(r) {
      row <- r$coefficients[r$coefficients$event.time == et, ]
      if (nrow(row) == 0) return("")
      paste0(fmt(row$estimate), stars_latex(row$p.value))
    }, character(1))
    se.vals <- vapply(results, function(r) {
      row <- r$coefficients[r$coefficients$event.time == et, ]
      if (nrow(row) == 0) return("")
      sprintf("(%s)", fmt(row$std.error))
    }, character(1))
    tex <- c(tex,
      sprintf("$t = %d$ & %s \\\\", et, paste(coef.vals, collapse = " & ")),
      sprintf(" & %s \\\\", paste(se.vals, collapse = " & ")))
  }

  tex <- c(tex, "\\hline",
    sprintf(paste0("\\multicolumn{%d}{l}{\\textbf{Panel B: ",
                   "Post-Treatment Average Effect}} \\\\"), n.cols),
    "\\hline")

  # Panel B: average effect, pre-trends test, pre-treatment mean, N.
  avg.est <- vapply(results, function(r)
    paste0(fmt(r$avg.post.effect$estimate),
           stars_latex(r$avg.post.effect$p.value)), character(1))
  avg.se <- vapply(results, function(r)
    sprintf("(%s)", fmt(r$avg.post.effect$std.error)), character(1))
  f.pvals <- vapply(results, function(r)
    fmt_num(r$pretrends.test$p.value, 3), character(1))
  pre.means <- vapply(results, function(r) fmt(r$pre.mean.treated), character(1))
  n.vals <- vapply(results, function(r)
    format(r$n.obs, big.mark = ","), character(1))

  tex <- c(tex,
    sprintf("%s & %s \\\\", effect.label, paste(avg.est, collapse = " & ")),
    sprintf(" & %s \\\\", paste(avg.se, collapse = " & ")),
    "\\hline",
    sprintf("Pre-Treatment Joint Sig. Test P-Value & %s \\\\",
            paste(f.pvals, collapse = " & ")),
    sprintf("Pre-Treatment Outcome Mean & %s \\\\",
            paste(pre.means, collapse = " & ")),
    sprintf("\\textit{N} & %s \\\\", paste(n.vals, collapse = " & ")),
    "\\hline\\hline", "\\end{tabular}")

  if (!is.null(output.file)) writeLines(tex, output.file)
  tex
}

#' LaTeX regression table from a list of lm() models
#'
#' General-purpose coefficient table: one column per model, coefficient over
#' standard error with significance stars, and an N row. Variable display names
#' and ordering follow `var.names` when supplied.
#'
#' @param models List of `lm()` (or compatible) model objects.
#' @param model.names Optional column titles (`NULL` = no title row).
#' @param var.names Optional named vector mapping variable names to display
#'   labels (its order sets the row order).
#' @param decimal.places Coefficient/SE precision (default 3).
#' @param output.file Optional path; if given, the LaTeX is written there.
#' @return The LaTeX string, invisibly.
#' @examples
#' m1 <- lm(mpg ~ wt, data = mtcars)
#' m2 <- lm(mpg ~ wt + hp, data = mtcars)
#' cat(demographic.table.output(list(m1, m2),
#'                              model.names = c("(a)", "(b)")))
#' @export
demographic.table.output <- function(models, model.names = NULL,
                                     var.names = NULL, decimal.places = 3,
                                     output.file = NULL) {
  n.models <- length(models)
  model.numbers <- paste0("(", seq_len(n.models), ")")

  model.results <- lapply(models, function(m) {
    coef.table <- summary(m)$coefficients
    data.frame(
      var = rownames(coef.table),
      coef = coef.table[, "Estimate"],
      se = coef.table[, "Std. Error"],
      pval = coef.table[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  })

  all.vars.raw <- unique(unlist(lapply(model.results, function(x) x$var)))
  if (!is.null(var.names)) {
    all.vars <- c(
      names(var.names)[names(var.names) %in% all.vars.raw],
      all.vars.raw[!all.vars.raw %in% names(var.names)]
    )
  } else {
    all.vars <- all.vars.raw
  }

  n.obs <- vapply(models, nobs, numeric(1))

  format.coef <- function(coef, pval) {
    if (is.na(coef)) return("")
    paste0(fmt_num(coef, decimal.places), add_stars(pval))
  }
  format.se <- function(se) {
    if (is.na(se)) return("")
    sprintf("(%s)", fmt_num(se, decimal.places))
  }

  col.spec <- paste0("l", paste(rep("c", n.models), collapse = ""))
  lines <- c(paste0("\\begin{tabular}{", col.spec, "}"), "\\hline\\hline")
  lines <- c(lines, paste0(paste(c("", model.numbers), collapse = " & "), " \\\\"))
  if (!is.null(model.names)) {
    lines <- c(lines,
               paste0(paste(c("", model.names), collapse = " & "), " \\\\"))
  }
  lines <- c(lines, "\\hline")

  for (v in all.vars) {
    display.name <- if (!is.null(var.names) && v %in% names(var.names)) {
      var.names[v]
    } else {
      v
    }
    coef.vals <- vapply(seq_len(n.models), function(i) {
      idx <- which(model.results[[i]]$var == v)
      if (length(idx) == 0) return("")
      format.coef(model.results[[i]]$coef[idx], model.results[[i]]$pval[idx])
    }, character(1))
    se.vals <- vapply(seq_len(n.models), function(i) {
      idx <- which(model.results[[i]]$var == v)
      if (length(idx) == 0) return("")
      format.se(model.results[[i]]$se[idx])
    }, character(1))
    lines <- c(lines,
      paste0(paste(c(display.name, coef.vals), collapse = " & "), " \\\\"),
      paste0(paste(c("", se.vals), collapse = " & "), " \\\\"))
  }

  lines <- c(lines, "\\hline",
    paste0(paste(c("\\textit{N}", format(n.obs, big.mark = ",")),
                 collapse = " & "), " \\\\"),
    "\\hline\\hline", "\\end{tabular}")

  tex.output <- paste(lines, collapse = "\n")
  if (!is.null(output.file)) writeLines(tex.output, output.file)
  invisible(tex.output)
}

################################################################################
# End of File
################################################################################
