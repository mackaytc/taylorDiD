# taylorDiD

Modern difference-in-differences and event-study helpers in one installable
package, providing a single consistent implementation across projects.

It wraps three estimators behind a common interface and a shared result object:

| Estimator | Backend | Status |
|---|---|---|
| BJS imputation (Borusyak, Jaravel, Spiess) | `didimputation` | stable |
| dCDH (de Chaisemartin & D'Haultfoeuille) | `DIDmultiplegtDYN` | stable |
| did2s (Gardner two-stage) | `did2s` | stable |

All estimators share the same house-style conventions: standard errors clustered
at the unit level, a joint pre-trends F-test, a post-treatment average effect,
the pre-treatment outcome mean, a single never-treated normalizer, and matching
event-study plot and two-panel LaTeX-table builders.

## Installation

```r
# install.packages("remotes")
remotes::install_github("mackaytc/taylorDiD")
```

The `Remotes:` field pulls `didimputation` from GitHub automatically. The dCDH
backend (`DIDmultiplegtDYN`) additionally needs the **polars** package, which is
distributed from r-universe rather than CRAN:

```r
install.packages("polars", repos = "https://rpolars.r-universe.dev")
```

To pin a version for a specific paper, install a tagged release, e.g.
`remotes::install_github("mackaytc/taylorDiD@v0.1.0")`.

## Quick start

```r
library(taylorDiD)

panel <- readRDS(system.file("extdata", "synthetic_panel.rds",
                             package = "taylorDiD"))

# BJS imputation event study
bjs <- did.event.study(panel, yname = "y", gname = "g", tname = "year",
                       idname = "id", pre.window = -4:-1, post.window = 0:3)

# dCDH event study (binary, discrete, or continuous treatment)
dcdh <- did.event.study.dyn(panel, yname = "y", dname = "treat", tname = "year",
                            idname = "id", n.pre = 3, n.post = 3)

# did2s (Gardner two-stage) event study
es <- did2s_event_study(panel, yname = "y", gname = "g", tname = "year",
                        idname = "id", pre.window = -4:-1, post.window = 0:3)

# Plot and LaTeX table from any result
plot_event_study(bjs, y.title = "Effect on y")
cat(tex_event_study_table(bjs, model.names = "BJS"), sep = "\n")
```

See `vignette("walkthrough", package = "taylorDiD")` for the full tour.

## Conventions

- **Never-treated coding** (`is_never_treated()`, `code_never_treated()`): `NA`,
  `Inf`, and non-positive timing are all never-treated. One normalizer, so the
  estimators never disagree silently.
- **Joint pre-trends** (`joint_pretrends_F()`): F = Wald / n_pre, p-value from an
  F with `df2 = Inf`. Full covariance where available (did2s), diagonal
  approximation otherwise (BJS/dCDH).
- **Post-treatment average**: an auxiliary static DiD on the
  never-/not-yet-/treated-in-window subset (BJS, did2s). dCDH uses the mean of
  the dynamic effects, since it has no clean static-DiD analog.
- **Result object** (`taylorDiD_es`): bundles coefficients, pre-trends test,
  average effect, pre-treatment mean, and N; consumed by the plot/table builders.

## License

MIT © Taylor Mackay
