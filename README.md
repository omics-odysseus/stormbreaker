# STORMBREAKER

STORMBREAKER is an R package for building and comparing modelling routes for
transcriptomic prediction. A route pairs a compatible predictor with optional
ordered feature-selection steps; it can also use a predefined gene set directly.
Compare candidate routes on shared folds, choose a route, then refit and save it
for use in new cohorts. The package supports binary classification and
right-censored survival prediction, with additional multiclass interfaces.

The prediction menus contain 92 binary and 88 survival configurations. Inspect
available methods and engines with `models("binary")` or `models("survival")`.

## Installation

R 4.3 or later is required:

```r
# install.packages("remotes")
remotes::install_github("omics-odysseus/stormbreaker",
  dependencies = NA, build_vignettes = FALSE)
```

Required dependencies are installed with the package. Install optional engines
for the methods you intend to use; the method menus report their status.

## Quick start

Run both examples with the bundled synthetic data. Install `glmnet` and
`ggplot2` if needed. The script reads the input files from the installed package.

```r
source(system.file("examples", "quick-start.R",
  package = "stormbreaker", mustWork = TRUE))

quick_start$binary$cv_auc[c("route", "estimate")]
quick_start$binary$selected
quick_start$binary$test_auc[c("metric", "estimate")]
quick_start$survival$metrics[c("metric", "estimate")]
```

The [binary example](inst/examples/scenarios/07-route-workflow.R) follows
the same objects from comparison through reuse:

1. **Build:** define the direct and Lasso routes. Set the Lasso selection
   rule to `lambda.min` for this example.
2. **Compare and choose:** fit both routes on the same three development
   folds, including gene selection inside each training split. Choose by ROC
   AUC, keeping the first listed route if the estimates tie.
3. **Refit:** fit the chosen route to the development data. Use the separate
   synthetic test cohort for evaluation.
4. **Reuse:** save and reload the fit, then supply IDs and predictor columns
   to predict. View the test ROC plot in `quick_start$binary$figures$ROC`.

The [survival example](inst/examples/scenarios/04-survival.R) fits a Cox
model, predicts in a separate synthetic cohort without using its outcomes
as input, and reports concordance, time-dependent AUC and Brier score.
See the [example-data notes](inst/extdata/demo-data.txt) for the input files
and columns.

For your own data, see the guide to
[input format](inst/guides/user-guide.md#input-format) and
[ordered pipelines](inst/guides/user-guide.md#selected-models-and-ordered-pipelines).

## Study examples

- [HNSC survival prediction](inst/guides/paper-cases.md#hnsc-survival-prediction)
- [Tumour–normal tissue classification](inst/guides/paper-cases.md#tumour-normal-tissue-classification)
- [Multicancer immunotherapy-response prediction](inst/guides/paper-cases.md#multicancer-immunotherapy-response-prediction)

The cases explain the candidate sets, selection rules and final routes in the
paper, with code to fit and apply those routes to prepared inputs.

## Documentation

- [User guide](inst/guides/user-guide.md): complete workflows, advanced use,
  runnable examples and parameter reference.
- [Parameter reference](inst/guides/user-guide.md#parameter-reference): arguments,
  defaults and returned objects.
- R help: `?bi_run`, `?surv_run`, `?finish_run` and `?multi_run`.
- Installed guide: build with `build_vignettes = TRUE` and open
  `vignette("stormbreaker-short-scenarios")`.
- [Changelog](NEWS.md): current and historical changes.

## Citation and license

Run `citation("stormbreaker")` for the software citation.
MIT; see [LICENSE.md](LICENSE.md) and [third-party notices](inst/COPYRIGHTS).
