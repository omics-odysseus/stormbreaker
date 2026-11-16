# STORMBREAKER

STORMBREAKER is an R package for transcriptomic prediction. It connects feature
selection, shared-fold model comparison, final fitting, external evaluation and
interpretation for binary classification and right-censored survival, with
additional multiclass interfaces.

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

Compare two models on the bundled synthetic development data, then evaluate the
selected model on a separate test cohort. Each row is one sample: `patient_id`
is the unique ID, `outcome` is 0/1 (1 is positive), and `gene_a`, `gene_b`
and `gene_c` are numeric predictors. The ROC plot uses `ggplot2`
(`install.packages("ggplot2")` if needed).

```r
library(stormbreaker)

read_demo <- function(name) {
  read.csv(system.file("extdata", name, package = "stormbreaker",
    mustWork = TRUE), colClasses = c(patient_id = "character"))
}
development <- read_demo("demo-binary-development.csv")
test <- read_demo("demo-binary-external-a.csv")
genes <- c("gene_a", "gene_b", "gene_c")
head(development[c("patient_id", "outcome", genes)], 3)

run <- bi_run(development, outcome = "outcome", id = "patient_id",
  predictors = genes, models = c("Logistic", "NaiveBayes"),
  folds = 3, seed = 42, plots = FALSE)
cv_auc <- result_table(run$comparison, metric = "roc_auc")
route_ids <- vapply(run$pipelines, function(route) route$id, character(1))
cv_auc$model <- names(route_ids)[match(cv_auc$pipeline_id, route_ids)]
cv_auc[c("model", "estimate")]

# Select using development-fold AUC only; the test cohort is not used here.
selected <- cv_auc$model[which.max(cv_auc$estimate)]
finished <- finish_run(run, selected, external = list(Test = test), plots = FALSE)
test_eval <- finished$evaluations$Test
result_table(test_eval, metric = "roc_auc")[c("metric", "estimate")]
print(plot_roc(test_eval))

model_file <- tempfile(fileext = ".rds")
saveRDS(finished$fit, model_file)
test_predictors <- test[c("patient_id", genes)]
prediction <- predict_final(readRDS(model_file), test_predictors, "patient_id")
head(prediction$predictions[c("patient_id", "probability")], 3)
```

You should see two development AUCs, one test AUC, a ROC curve, and a preview
of the predictions. The reloaded model needs only IDs and predictor columns,
not outcomes. To use your own data, replace the two CSVs and specify your
predictor columns in `genes`. Keep the test data out of model selection.

The full menu and ordered feature-selection pipelines are covered in the
[user guide](inst/guides/user-guide.md#selected-models-and-ordered-pipelines).

## Documentation

- [User guide](inst/guides/user-guide.md): complete workflows, advanced use,
  runnable examples and parameter reference.
- R help: `?bi_run`, `?surv_run`, `?finish_run` and `?multi_run`.
- Installed guide (only when built with `build_vignettes = TRUE`):
  `vignette("stormbreaker-short-scenarios")`.
- [Changelog](NEWS.md): current and historical changes.

## Citation and license

Run `citation("stormbreaker")` for the software citation.
MIT; see [LICENSE.md](LICENSE.md) and [third-party notices](inst/COPYRIGHTS).
