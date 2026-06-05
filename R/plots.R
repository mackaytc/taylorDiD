################################################################################
#
# FILE: plots.R
#
# OVERVIEW: Event-study plotting. plot_event_study() is the shared ggplot with
# the house defaults (point estimates, 95% CI errorbars, dashed zero line, serif
# font). The two retained wrappers -- did.event.study.plot() (BJS-style, with a
# diagnostic caption) and did.event.study.plot.dyn() (dCDH-style, with y-axis
# titles) -- preserve the original call signatures and return both the plot(s)
# and the two-panel LaTeX table.
#
# AUTHOR: Taylor Mackay (tmackay@fullerton.edu)
#
################################################################################

#' Event-study plot with house defaults
#'
#' Point estimates with 95% confidence-interval errorbars, a dashed zero line,
#' and a serif theme, with integer breaks on the event-time axis.
#'
#' @param result A `taylorDiD_es` object (uses its `coefficients` table, which
#'   must carry `event.time`, `estimate`, `conf.low`, `conf.high`).
#' @param y.title Optional y-axis title (default none).
#' @param x.title X-axis title (default `"Event Time"`; pass `NULL` to omit).
#' @param title Optional plot title.
#' @param caption Optional caption (used for diagnostic annotations).
#' @param point.fill Fill color for the point estimates (default `"coral3"`).
#' @return A `ggplot` object.
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
#' plot_event_study(es, y.title = "Effect on outcome")
#' @export
plot_event_study <- function(result, y.title = NULL, x.title = "Event Time",
                             title = NULL, caption = NULL,
                             point.fill = "coral3") {
  coefs <- as.data.table(result$coefficients)
  x.breaks <- seq(min(coefs$event.time), max(coefs$event.time), by = 1)

  ggplot(coefs, aes(x = .data$event.time, y = .data$estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = .data$conf.low, ymax = .data$conf.high),
                  width = 0.2, color = "gray40", linewidth = 0.625) +
    geom_point(shape = 21, fill = point.fill, color = "black",
               stroke = 0.75, size = 3) +
    scale_x_continuous(breaks = x.breaks) +
    labs(x = x.title, y = y.title, title = title, caption = caption) +
    theme_minimal() +
    theme(
      text = element_text(family = "serif", size = 14),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks.x = element_line(),
      axis.text = element_text(size = 14),
      plot.title = element_text(size = 14, face = "italic", hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0.5, face = "italic",
                                  lineheight = 1.2)
    )
}

# Diagnostic two-line caption used by the BJS-style plot.
es_caption <- function(result) {
  pt <- result$pretrends.test
  ap <- result$avg.post.effect
  line1 <- sprintf("Groups: %d total, %d treated | Pre-trends F-stat: %.2f (p = %.3f)",
                   result$n.groups, result$n.treated, pt$f.stat, pt$p.value)
  line2 <- sprintf(paste0("Pre-treatment mean: %.3f | Avg post effect: %.3f ",
                          "(p = %.4f) | Pct change: %.1f%%"),
                   result$pre.mean.treated, ap$estimate, ap$p.value,
                   ap$pct.change)
  paste(line1, line2, sep = "\n")
}

#' Event-study plot and LaTeX table (BJS-style)
#'
#' Produces a diagnostic-captioned event-study plot per model and the two-panel
#' LaTeX table. Retained for backward compatibility with existing project code;
#' delegates to [plot_event_study()] and [tex_event_study_table()].
#'
#' @param es.results A single `taylorDiD_es` object or a list of them.
#' @param model.names Optional column/title names (one per model).
#' @param decimal.places Table precision (default 4).
#' @param output.file Optional path for the LaTeX table.
#' @return List with `plots` (a `ggplot` or list of them) and `table` (LaTeX
#'   lines).
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
#' out <- did.event.study.plot(es, model.names = "Outcome")
#' @export
did.event.study.plot <- function(es.results, model.names = NULL,
                                 decimal.places = 4, output.file = NULL) {
  if (!is.null(es.results$coefficients)) es.results <- list(es.results)
  n.models <- length(es.results)
  if (is.null(model.names)) model.names <- paste0("Model ", seq_len(n.models))

  plots <- lapply(seq_len(n.models), function(m) {
    plot_event_study(es.results[[m]],
                     title = paste("Event Study:", model.names[m]),
                     caption = es_caption(es.results[[m]]))
  })
  if (n.models == 1) plots <- plots[[1]]

  table <- tex_event_study_table(es.results, model.names = model.names,
                                 decimal.places = decimal.places,
                                 output.file = output.file)
  list(plots = plots, table = table)
}

#' Event-study plot and LaTeX table (dCDH-style)
#'
#' Produces a clean event-study plot per model (no caption, with per-model y-axis
#' titles) and the two-panel LaTeX table. Retained for backward compatibility;
#' delegates to [plot_event_study()] and [tex_event_study_table()].
#'
#' @inheritParams did.event.study.plot
#' @param y.titles Optional y-axis titles (one per model; default: `model.names`).
#' @return List with `plots` and `table`.
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
#'   n.obs = 500L, n.groups = 50L, n.treated = 20L, estimator = "dcdh")
#' out <- did.event.study.plot.dyn(es, model.names = "Outcome")
#' @export
did.event.study.plot.dyn <- function(es.results, model.names = NULL,
                                     y.titles = NULL, decimal.places = 4,
                                     output.file = NULL) {
  if (!is.null(es.results$coefficients)) es.results <- list(es.results)
  n.models <- length(es.results)
  if (is.null(model.names)) model.names <- paste0("Model ", seq_len(n.models))
  if (is.null(y.titles)) y.titles <- model.names

  plots <- lapply(seq_len(n.models), function(m) {
    plot_event_study(es.results[[m]], y.title = y.titles[m], x.title = NULL)
  })
  if (n.models == 1) plots <- plots[[1]]

  table <- tex_event_study_table(es.results, model.names = model.names,
                                 decimal.places = decimal.places,
                                 output.file = output.file)
  list(plots = plots, table = table)
}

################################################################################
# End of File
################################################################################
