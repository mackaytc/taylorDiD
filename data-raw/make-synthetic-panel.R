################################################################################
#
# FILE: make-synthetic-panel.R
#
# OVERVIEW: Generates the deterministic synthetic staggered-adoption panel that
# ships in inst/extdata for examples, tests, and the vignette. One row is a
# unit-year. Treated units adopt at a cohort year (g); never-treated units carry
# g = Inf. Outcome y has unit and year components plus a dynamic post-treatment
# effect that ramps with event time.
#
# INPUTS:  N/A (fully simulated)
# OUTPUTS:
#  - inst/extdata/synthetic_panel.rds
#  - inst/extdata/synthetic_panel.csv
#
# AUTHOR: Taylor Mackay (tmackay@fullerton.edu)
#
################################################################################

library(data.table)

out.dir <- "C:/Users/macka/OneDrive/Documents/github-repos/taylorDiD/inst/extdata"

set.seed(20260605)

n.units <- 60
years <- 2000:2012

# Cohort assignment: three treated cohorts plus never-treated (Inf), weighted so
# roughly a third of units never adopt. Cohorts sit far enough inside the panel
# to support a -4:-1 pre-window and a 0:3 post-window.
cohort.pool <- c(2004, 2006, 2008, Inf, Inf)
unit.cohort <- sample(cohort.pool, n.units, replace = TRUE)
unit.fe <- rnorm(n.units, mean = 0, sd = 1)
unit.intensity <- sample(1:3, n.units, replace = TRUE)

panel <- CJ(id = seq_len(n.units), year = years)
panel[, g := unit.cohort[id]]
panel[, event.time := ifelse(is.finite(g), year - g, NA_integer_)]
panel[, treat := as.integer(is.finite(g) & year >= g)]

# d is a non-binary treatment intensity (0 pre, a unit-specific level post),
# used to exercise the dCDH continuous/non-binary path; treat is the binary
# absorbing version.
panel[, d := ifelse(treat == 1L, unit.intensity[id], 0L)]

# Dynamic effect ramps with event time among post-treatment observations.
panel[, tau := ifelse(!is.na(event.time) & event.time >= 0,
                      0.5 * (event.time + 1), 0)]
panel[, y := unit.fe[id] + 0.1 * (year - min(years)) + tau +
        rnorm(.N, mean = 0, sd = 0.5)]

panel <- panel[, .(id, year, g, treat, d, y)]
setorder(panel, id, year)

saveRDS(panel, file.path(out.dir, "synthetic_panel.rds"))
fwrite(panel, file.path(out.dir, "synthetic_panel.csv"))

################################################################################
# End of File
################################################################################
