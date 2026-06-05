# taylorDiD

Modern difference-in-differences and event-study helpers, consolidated into one
installable package so every project loads the same code instead of hand-copying
helper files between `code/utility-code/` folders.

It wraps three estimators behind a common interface and a shared result object:

| Estimator | Backend | Status |
|---|---|---|
| BJS imputation (Borusyak, Jaravel, Spiess) | `didimputation` | stable |
| dCDH (de Chaisemartin & D'Haultfoeuille) | `DIDmultiplegtDYN` | stable |
| did2s (Gardner two-stage) | `did2s` | **experimental** |

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

# did2s (experimental)
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
  never-/not-yet-/treated-in-window subset (BJS, did2s). dCDH instead uses the
  mean of the dynamic effects, since it has no clean static-DiD analog -- a
  deliberate, documented exception.
- **Result object** (`taylorDiD_es`): bundles coefficients, pre-trends test,
  average effect, pre-treatment mean, and N; consumed by the plot/table builders.

## did2s status (experimental)

did2s is scaffolded this release: the ported helpers are generalized off their
original hardcoded column names and unit-tested, and `did2s_event_study()` runs
the full flow, but the module has not been validated against a reference
`did2s::did2s()` run. The GAP LIST lives in `?did2s_event_study`:

1. Validate against a known `did2s::did2s()` run on a real panel.
2. Confirm the auxiliary-static-DiD post average (the equal-weighted
   `linear.combo` is retained internally for A/B comparison).
3. Confirm the generalized column handling on a second dataset.
4. Confirm the full-covariance F-test convention for pre-trends.
5. Decide whether the project-specific `split.source.ids` helper belongs here.
6. Write a dedicated did2s walkthrough.

Prefer the BJS or dCDH estimators for finalized results.

## License

MIT © Taylor Mackay
