# STORMBREAKER user guide

STORMBREAKER connects transcriptomic model comparison, feature selection, final
fitting, external evaluation and prediction explanation. The main workflows
cover binary classification and right-censored survival; multiclass analysis
uses its own task and comparison interfaces. This guide describes version
**0.3.47.1**.

## Contents

1. [Installation and data](#installation-and-data)
2. [Model comparison](#model-comparison)
3. [Final fitting and external evaluation](#final-fitting-and-external-evaluation)
4. [Diagnostics and explanations](#diagnostics-and-explanations)
5. [Saving and reporting](#saving-and-reporting)
6. [Advanced workflows](#advanced-workflows) ([clinical binary workflow](#clinical-binary-workflow))
7. [Multiclass and independent analysis](#multiclass-and-independent-analysis)
8. [Runnable examples](#runnable-examples)
9. [Parameter reference](#parameter-reference)

## Installation and data

### Installation

R 4.3 or later is required. Install from GitHub:

```r
# install.packages("remotes")
remotes::install_github("omics-odysseus/stormbreaker",
  dependencies = NA, build_vignettes = FALSE)
library(stormbreaker)
```

Required dependencies are installed with the package. Browse optional engines
with `models("binary")`, `models("survival")` or `models("multiclass")`.
The lower-level `sb_model_catalog("binary")` reports `model_id`, `r_packages`
and `dependency_status`; it accepts binary and survival tasks.

Install the engines for the methods you intend to use. Typical optional
packages include `glmnet` for regularized regression, `randomForest` or `ranger`
for forests, `mda` for FDA/PDA/MDA, and `ggplot2` for plots. The six short examples
also use `patchwork` and `glmnet`. Cohort synthesis uses `metafor`.

Bioconductor engines such as `mixOmics`, `cancerclass` and `switchBox` use
[BiocManager](https://www.bioconductor.org/install/). See the
[CatBoost R instructions](https://catboost.ai/docs/en/installation/r-installation-binary-installation)
for that engine. The historical `fastAdaboost` package is available in its
[CRAN archive](https://cran.r-project.org/src/contrib/Archive/fastAdaboost/).
Kernel-survival configurations require `survivalmodels` and `distr6`; consult
the [distr6 distribution](https://mlr-org.r-universe.dev/distr6) when it is not
available from the selected CRAN repository.

For a supplied local source archive, use
`install.packages(file.choose(), repos = NULL, type = "source")`; use
`type = "win.binary"` for a compatible Windows binary. Check the archive's
version and R requirements. After updating a loaded package, restart R and
check `packageVersion("stormbreaker")` and `find.package("stormbreaker")`.
Use `lib = your_project_library` for an isolated installation.
Use `sb_start()` to open and copy local workflow examples.

### Input format

Each row is one observation, with a unique nonmissing identifier and numeric
predictor columns. Binary outcomes are 0/1, with 1 positive. Survival inputs
contain follow-up time and an event indicator: 1 for an event, 0 for censoring.
Keep distinct endpoints, such as overall and disease-free survival, separate.

Use **months** consistently in these examples. Convert days with
`time / 30.4375` for average months, or years with `time * 12`. If a study uses
another conversion convention, apply it consistently to all cohorts and
prediction horizons. Convert the values before setting `time_unit`.

Transpose a gene-by-sample expression matrix before constructing a task.
Align expression and clinical data by identifier. `sb_read_cohort()` preserves
textual IDs such as `001` and matches a separate clinical table by ID:

```r
cohort <- sb_read_cohort("expression.csv", "clinical.csv")
genes <- c("EGFR", "TP53")
task <- bi_task(cohort, "patient_id", "outcome", predictors = genes)
```

The CSV example requires your own files. Prepare missing values, feature
mapping and expression scales for the study; apply data-dependent preprocessing
within the relevant analysis folds. Supply every required external feature on
the same scale and with the same definition used in development.
For repeated observations, use unique row IDs and declare the patient/group
column through the task's `group_col` argument.

The main examples below use packaged synthetic data:

```r
development <- bi_example(n = 300)
validation <- bi_example(n = 100, seed = 2)
survival_data <- surv_example(n = 300)
```

The [example-data notes](../extdata/demo-data.txt) describe the CSV columns,
generation seeds and units. Declare candidate genes, primary metrics, folds,
thresholds, horizons and the model-selection rule before examining test outcomes.

## Model comparison

For a first run, use the [two-route workflow](#runnable-examples) below. It
compares direct Logistic with Lasso followed by Logistic on the same folds. The
examples below without an explicit model list compare the full menu.

### Binary classification

```r
binary <- bi_run(
  development, outcome = "outcome", id = "patient_id",
  decision_thresholds = seq(0.05, 0.95, by = 0.05)
)
result_table(binary$comparison, metric = "roc_auc")
binary$plots$performance
names(binary$pipelines)
```

The default compares the complete binary menu and supported two-stage
pipelines, with built-in model settings and shared five-fold cross-validation.
This can involve many routes. Pass `predictors = genes` when the data contain
other descriptive columns that are not candidate predictors.

### Survival prediction

```r
survival_run <- surv_run(
  survival_data, time = "time", event = "event", id = "patient_id"
)
result_table(survival_run$comparison, metric = "harrell_c")
survival_run$plots$performance
```

The default includes the full survival menu, two-stage combinations and built-in
three-stage pipelines. Add `prediction_horizon = 6` for a six-month probability
target, and `decision_thresholds` for decision curves. Choose a horizon
supported by follow-up in the evaluation cohorts. Risk scores and absolute
probabilities are different prediction outputs.

### Selected models and ordered pipelines

```r
selected <- bi_run(
  development, outcome = "outcome", id = "patient_id",
  models = c("Logistic", "RF", "SVM")
)

routes <- list(
  Direct = bi_pipeline("Logistic"),
  Selected = bi_pipeline(c("Lasso", "Logistic"))
)
custom <- bi_run(
  development, outcome = "outcome", id = "patient_id", pipelines = routes
)
```

Earlier binary stages extract features and the last predicts. A three-stage
route therefore needs two supported extractors followed by a predictor.
Survival uses `surv_pipeline()` and ordered `surv_selector()` steps. With an
explicit model pool, set `stages = c(1, 2, 3)` to request supported depths;
explicit model lists otherwise default to single-stage comparisons.

Selection and method-specific tuning use the relevant training partitions.
An empty selected set remains a failed-route record. Inspect `run$issues` and
`result_table(run$comparison, "runs")`; retain failures when comparing methods.

Model-specific options belong in the pipeline specification. For example:

```r
lasso_route <- sb_pipeline("lasso_binary_nonzero", "gbm_binary",
  selector_options = list(lambda_rule = "lambda.min"))
step_route <- sb_pipeline("stepglm_binary_aic", "gaussian_naive_bayes",
  selector_options = list(direction = "forward"))
forest_route <- sb_pipeline("rf_binary_importance", "logistic",
  selector_options = list(threshold = 0))
```

The documented defaults include `lambda.1se` for Lasso selection, both
directions for stepwise AIC, 100 iterations and 0.1 shrinkage for glmBoost
selection, and 500 trees with permutation importance above zero for the forest
selector. Consult `?sb_pipeline` and the menus for supported combinations and
parameters. In `sb_compare_models()`, omitting `models` uses its older small
starting set; explicitly pass `models = NULL` for that interface's full menu.

## Final fitting and external evaluation

Choose a successful route using the prespecified development or validation
rule. Use its exact name from `names(binary$pipelines)`:

```r
finished <- finish_run(
  binary, "Lasso + Logistic",
  external = list(Validation = validation)
)
result_table(finished$evaluations$Validation, metric = "roc_auc")
head(finished$predictions$Validation)
```

`finish_run()` retains the comparison and reuses an existing final fit or fits
the chosen pipeline on all development data. `$fit` is the saved model,
`$genes` contains stage-wise features, and `$predictions` and `$evaluations`
contain results by cohort.

```r
views <- clinical_plots(
  finished, detail = "full",
  thresholds = c(0.1, 0.3, 0.5, 0.7, 0.9), confusion_threshold = 0.5
)
views$plots$ROC
views$plots$Calibration
views$plots[["Decision curve"]]
views$plots[["Precision-recall"]]
views$plots[["Probability distribution"]]
views$plots[["Threshold performance"]]
views$plots[["Confusion matrix"]]
```

Declare thresholds during modeling for decision curves. For the confusion
matrix, use a threshold chosen in advance; the example uses 0.5. Check
`views$status` for figure availability and reasons.

```r
survival_validation <- surv_example(n = 100, seed = 2)
survival_finished <- finish_run(
  survival_run, "Cox", external = list(Validation = survival_validation)
)
survival_finished$plots$external$Validation
```

Risk-group cutoffs come from development. Horizon-specific calibration needs
saved survival probabilities; survival decision curves also need declared
thresholds. `survival_probability` is the probability of remaining alive;
decision curves use its event-risk complement.

For separate prediction and outcome handling, use an existing task and route:

```r
fit <- fit_final(task, bi_pipeline("Logistic"), capture_diagnostics = TRUE)
prediction <- predict_final(fit, new_samples, "patient_id")
endpoint <- bi_endpoint(new_outcomes, "patient_id", "outcome")
evaluation <- eval_external(prediction, endpoint)
```

Prepare `task`, `new_samples` and `new_outcomes` as above. Predict using IDs
and the columns in `fit$required_features`, then supply outcomes separately.
Use matching IDs for predictions and endpoints.

The corresponding task-based survival workflow declares the same time unit
and horizon during development and evaluation:

```r
survival_task <- surv_task(survival_data, "patient_id", "time", "event",
  c("gene_a", "gene_b", "gene_c"), "month", prediction_horizon = 6)
survival_fit <- fit_final(survival_task, surv_pipeline("Cox"))
survival_prediction <- predict_final(survival_fit,
  survival_validation[c("patient_id", "gene_a", "gene_b", "gene_c")],
  "patient_id")
survival_endpoint <- surv_endpoint(survival_validation,
  "patient_id", "time", "event", "month", prediction_horizon = 6)
survival_evaluation <- eval_external(survival_prediction, survival_endpoint)
```

To compare several locked models on **one** cohort, pass their prediction
objects and its endpoint to `compare_external(predictions, endpoint)`.

### Reading metrics

`result_table()` reads saved metrics, predictions or run records. It also accepts
whole repeated-CV, cohort-comparison and stacking results; cohort roles,
repetitions, pipeline identities, horizons and units remain separate.
`sb_failure_table()` extracts recorded failures from these objects.
`sb_export_failures(x, "new-failures.csv")` writes their saved cohort, pipeline,
stage, fold, metric and reason fields. Review error text before sharing it.

Read `estimate` alongside `status`, `run_status`, `scope` and `aggregation`.
Use pipeline and patient columns together to identify rows across models.
External evaluation tables include outcomes. Set `show_unavailable = TRUE` to
show unavailable labels; missing and nonfinite estimates remain unavailable.

Compare routes using the saved assessment-fold summaries. Binary ROC AUC and
Harrell C aggregate within-fold pairs, and average precision uses its fold-aware
definition. Treat variation across folds as descriptive.

Match the saved target and unit when comparing `time_auc`, `uno_c` and
`ipcw_brier`. The weighted IPCW Brier estimate can exceed one. Check
`not_requested` for successful horizon-free results and `not_saved` for metrics
absent from the stored results. Plots read the saved estimates.

## Diagnostics and explanations

### Native and stage-wise diagnostics

```r
finished$genes
diagnostics <- plots(finished$fit)
names(diagnostics$plots)
diagnostics$unavailable
diagnostics$plots$selector_1_cv
diagnostics$plots$selector_1_coefficient_path
```

Capture the native diagnostics you plan to inspect when fitting. Use
`plots(fit)` to find the available diagnostics, including the first selector's
CV and coefficient path in the Lasso example. Supply query data for diagnostics
that use it: `plots(fit, new_data = query)`.

To inspect a specific stage or learner, use `native_diag(fit, complete = TRUE)`,
`path_diag()`, `forest_diag()` or another supported extractor, then `plot_diag()`.

### SHAP and coefficients

```r
shap <- model_shap(
  finished$fit, validation, "patient_id", background_data = development
)
shap_plots <- list(
  importance = plot_shap_imp(shap),
  distribution = plot_shap(shap, feature_values = validation, id_col = "patient_id"),
  waterfall = plot_shap_one(shap)
)
shap_plots$importance
shap_plots$distribution
shap_plots$waterfall
```

`model_shap()` computes approximate permutation contributions using a declared
background. Binary explanations use positive-class probability; survival
explanations use the saved risk score. `additive_shap()` instead gives exact
contributions for supported linear models, on log-odds or linear log-risk scales.
Check `score_scale` when interpreting either output. Use `plot_shap_dep()` for
a feature dependence plot and reuse the computed SHAP object across views.

`model_coefs(fit, reference_data = reference_features)` retains original-scale
coefficients, zeros and the intercept. `plot_coefs(..., scale = "reference_sd")`
rescales coefficients for display only. `nomogram()` and `plot_nomogram()`
represent supported original linear models and their reference ranges.

### Survival curves

Use the horizon-aware linear Cox fit above and the original development task.
`add_surv_curve()` checks identifiers, time unit and the saved Breslow or Efron
baseline point, then supplies the complete baseline using the same coefficients.

```r
curve_fit <- add_surv_curve(survival_fit, survival_task)
curves <- surv_curves(curve_fit, survival_validation, "patient_id")
plot_surv(curves, rows = 1)
```

The baseline corresponds to a linear predictor of zero. The plotted curves
are model predictions within observed development follow-up.

### Other saved-result views

| Question | Functions | Required material |
|---|---|---|
| Discrimination and precision-recall | `plot_roc()`, `plot_pr()` | Predictions and outcomes |
| Calibration and clinical utility | `plot_calib()`, `plot_dca()` | Saved coordinates, thresholds and any survival horizon |
| Model/cohort performance | `plot_perf()`, `plot_heatmap()`, `plot_rank()` | Compatible metrics and declared ranking cohorts |
| Fixed-threshold errors | `plot_cm()`, `plot_thresholds()` | Predictions and declared thresholds |
| Survival groups and timed ROC | `surv_risk()`, `surv_roc()` | Development risk reference, follow-up and requested times |
| Feature stability | `plot_freq()`, `plot_jaccard()`, `plot_upset()` | Saved resampling selections |
| External uncertainty | `plot_ext_ci()`, `plot_ext_diff()` | Saved bootstrap results |
| Cohort synthesis | `plot_forest()` | Saved meta-analysis results |

Figures use concise titles, axes, labels and legends. Detailed methods and
interpretation belong in captions; returned notes, captions and status tables
retain the relevant information. Keep apparent training, development CV and
external evaluation clearly identified. Missing values stay missing rather
than becoming zero.

## Saving and reporting

Algorithm figures are collected by `finish_run()`. Attach separately created
clinical and SHAP figures before exporting:

```r
finished$plots$clinical <- views$plots
finished$plots$shap <- shap_plots
finished$explanation <- shap
output_dir <- "prediction-results"
export_run(finished, output_dir)
saveRDS(finished$fit, file.path(output_dir, "model.rds"))

fit <- readRDS(file.path(output_dir, "model.rds"))
prediction <- predict_final(fit, validation, "patient_id")
saved_diagnostics <- plots(fit)
```

Use a new output directory. Exported CSVs include model comparisons, stage-wise
genes, cohort predictions, evaluations and issues. `figures.csv` maps plot names
to PDF files. The model retains captured diagnostic materials;
`explanation.rds` retains a supplied SHAP result. Save the whole result with `saveRDS(finished, ...)`
when the comparison and evaluations are needed later.

Keep the recorded engine versions when moving models between environments.
For example, newly saved XGBoost states use JSON, while older saved states may
need the engine version that wrote their original format. Saved models,
diagnostics and predictions can contain sample-level information; inspect them
before sharing.

### Reports and result browsing

```r
clinical_report(views, file.path(output_dir, "clinical-report.html"))
sb_render_report(list(Development = binary$comparison,
  Validation = finished$evaluations$Validation),
  file.path(output_dir, "results-report.html"))
panel <- plot_panel(views, c("ROC", "Calibration"))
ggplot2::ggsave(file.path(output_dir, "clinical-panel.pdf"), panel,
  width = 12, height = 5)
```

`clinical_report(..., complete = TRUE)` includes the complete supplied clinical
view. Add native diagnostics, selector evidence, SHAP and nomograms through
`clinical_plots()` when needed. `sb_render_model_report()` writes a model report
from a saved fit, reference features and query samples; `sb_render_visual_report()`
writes an offline gallery of supplied plots. Report paths must be new.

`sb_view_results()` opens a model-by-cohort page with metric selection, search,
cohort visibility and CSV/SVG export. Use `open = FALSE` and a new `path` for
noninteractive output. Sort or import saved summaries to browse results.
Validation averages use compatible displayed validation cohorts only. CSV
retains filtered routes; SVG and printing follow the displayed page.

## Advanced workflows

### Clinical binary workflow

The synthetic example compares six model families with shared five-fold CV
on the development partition. It declares four fixed configurations and a
small grid for LightGBM and Extra Trees. The selected configuration is locked
on development CV before predicting the 30% holdout.

```r
library(stormbreaker)
source(system.file("examples", "advanced-examples.R",
  package = "stormbreaker", mustWork = TRUE))
example <- run_advanced_example("clinical-binary")
example$split$counts
example$candidates[c("candidate", "status", "roc_auc")]
example$metrics

# Optional: run and save the figures, tables and fitted model in a new directory.
# report <- stormbreaker_clinical_example("clinical-results",
#   figure_format = "png", font_family = "Times New Roman")
```

Use `sb_split_binary(task, train_prop = .7, seed = 82L)` with your own
`bi_task()`. Use unique observation IDs and set `group_col` for patients with
repeated observations. The split keeps each group together and both outcomes
represented; the realized proportions can differ from 70:30 for groups of
different sizes or small strata. Fit and select on `split$development`, then
pass the unchanged internal holdout to
`finish_run(..., external = list(Holdout = split$assessment))`. The assessment
plots use the holdout label.

Specify `selector = "rf_binary_importance"` and
`selector_options = list(top_n = 5L, spearman_cutoff = .70)` for optional RF
ranking and redundancy pruning inside each training fold and on final
development data. Set the pruning cutoff and probability threshold from your
analysis plan; the example uses 0.70 and 0.5, respectively. Pruning can retain
fewer than five predictors. Report the patient-bootstrap AUC interval with the
ROC curve. For repeated observations, define the evaluation unit and use
group-aware intervals.

The example exports an outcome bar chart, RF importance, a training
correlation heatmap, the existing clinical curves and confusion matrix, and
five SHAP views. Its SHAP reference comes from development predictors;
the heatmap shows signed contributions to the fitted predictions. Review
calibration with its curves, Brier score and calibration summaries.

#### Compare more binary methods

Inspect every binary menu entry and build one configuration per compatible
predictor under the same upstream selector, using the registered learner
settings. Start with a small named subset to limit runtime. Change `selector`
for a separate run with a different upstream method.

```r
d <- bi_example(120L, seed = 42L)
task <- bi_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
split <- sb_split_binary(task)
plan <- sb_binary_workflow_plan()    # Inspect only; no training or files
table(plan$candidates$status)
# Specify models to limit the batch; omit models to use the full plan.
plan <- sb_binary_workflow_plan(models = c("Logistic", "Lasso", "NaiveBayes"))
run <- bi_run(split$development, pipelines = plan$pipelines,
  folds = 5L, seed = 11L, plots = FALSE)
run$comparison$comparison
run$issues

# Choose using development evidence only; "Logistic" is illustrative.
finished <- finish_run(run, "Logistic", external = list(Holdout = split$assessment),
  plots = FALSE, capture_diagnostics = TRUE)
finished$metrics$Holdout

# Later: draw only what is needed, from saved evaluation or diagnostic data.
finished$plots <- list(Holdout = list(ROC = plot_roc(finished$evaluations$Holdout)))
sb_run_figure_keys(finished)
# export_run(finished, "selected-figures", figure_keys = "/Holdout/ROC",
#   include_data = FALSE, figure_format = "png")
# export_run(finished, "tables-only", figure_keys = character(0))
```

In this package version, the binary menu contains 92 prediction configurations
and 12 selector-only entries. RF is registered for 90 prediction configurations.
Pre and Cancerclass are limited to routes without an upstream selector; inspect
them with `selector = "none"` or explicitly include them as distinct routes using
`include_unfiltered = TRUE`. Review optional package locations and incompatible
pairs before fitting, then check the recorded fit failures. The batch uses shared
folds, recalculates screening on each training fold and saves out-of-fold
predictions for each candidate. Refit the chosen candidate for the holdout, then
request supported diagnostics or SHAP for that fit. Use `export_run()` to write
figures and tables.

### Task-based repeated validation

The following small example is independent of the full-menu examples above.

```r
adv_development <- bi_example(n = 120, seed = 42)
adv_genes <- c("gene_a", "gene_b", "gene_c")
adv_task <- bi_task(adv_development, "patient_id", "outcome", adv_genes)
adv_routes <- list(
  Direct = bi_pipeline("Logistic"),
  Selected = bi_pipeline(c("Lasso", "Logistic"))
)
adv_run <- bi_run(adv_task, pipelines = adv_routes,
  folds = 5, seed = 42, plots = FALSE)
result_table(adv_run$comparison, "runs")

repeated <- sb_compare_repeated_cv(adv_task, adv_routes,
  folds = 5, seeds = c(11, 12), workers = 1)
repeated$comparison
```

Each seed defines one repetition; pipelines share splits within a repetition.
Use the saved results and out-of-fold predictions to assess variation across
development splits.

`workers = 1` is serial. `workers = 2` runs two repetitions in separate local
processes and retains input order. Allow memory for two data copies and active
fits. Built-in parallel execution keeps engines single-threaded; use serial
execution when the learner requires the caller's RNG stream.

### Multiple cohorts and explicit selection

```r
validation_a <- bi_example(n = 80, seed = 43)
validation_b <- bi_example(n = 70, seed = 44)
validation_a$patient_id <- paste0("val_a_", validation_a$patient_id)
validation_b$patient_id <- paste0("val_b_", validation_b$patient_id)
cohort_result <- sb_compare_cohorts(
  adv_task,
  cohorts = list(
    Validation_A = bi_task(validation_a, "patient_id", "outcome", adv_genes),
    Validation_B = bi_task(validation_b, "patient_id", "outcome", adv_genes)
  ),
  roles = c(Validation_A = "validation", Validation_B = "validation"),
  pipelines = adv_routes, folds = 5, seed = 42, view = FALSE
)
cohort_result$cohorts$Validation_A$comparison
cohort_result$cohorts$Validation_B$comparison
```

Supply cohorts with disjoint patients and groups. For survival, use matching
time units and horizons. Fit routes on development data, retain its out-of-fold
results, and apply each fitted route to the other cohorts.

Label each cohort `"validation"`, `"test"` or `"external"` (the default).
Use `sb_select_model(cohort_result)` to select by the declared metric from the
validation cohorts. It weights eligible cohorts equally and resolves ties in
declaration order. Keep fitting confined to development data and reserve a
separate test cohort for evaluation. Review `$model_selection$ranking` on the
returned model.

Alternatively choose a route explicitly. The successful `Direct` example is
evaluated on a separately generated test cohort:

```r
test_data <- bi_example(n = 90, seed = 45)
test_data$patient_id <- paste0("test_", test_data$patient_id)
adv_finished <- finish_run(adv_run, "Direct",
  external = list(Test = test_data))
result_table(adv_finished$evaluations$Test, metric = "roc_auc")
```

Use the same export and reload steps as in [Saving and reporting](#saving-and-reporting).
Review failed routes in the comparison results before choosing a route.

### Stacking and held-out cohorts

Whole results from `sb_run_stack_cv()` and `sb_run_internal_external_cv()`
use the same `result_table()`, plotting, report and failure-export interfaces.
Stack views identify the member and stack roles from outer CV.
Use internal-external CV for method development and a separate final test set
for independent evaluation. Context columns such as `repeat_id`,
`heldout_cohort_id`, `model_role`, `validation_scheme` and `view_id` keep saved
entries distinct by repetition, cohort and survival target. A `view_id`
identifies a results view.

### Computation reuse

`compare_cv()` defaults to `reuse_selectors = TRUE`. Eligible identical
fold-local selections, supported CV/final fits and some survival-curve
calculations can share computation. Matching includes ordered features,
samples, outcomes, method settings, seeds and relevant RNG/diagnostic context.
Outer analysis folds, inner analysis folds and final development data remain
separate. Ordered chains can share identical prefixes; later stages still run.

The cache retains successful fits without warnings. Pipeline identities,
configured training budgets, predictions and diagnostics remain distinct.
Each R process maintains its own in-memory context for the current scope.
CV and final learner entries share a **256 MiB** cache budget with
least-recently-used eviction.

```r
unshared <- compare_cv(adv_task, adv_routes, folds = 5, reuse_selectors = FALSE)
measured <- sb_with_selector_reuse({
  cv_results <- lapply(adv_routes, function(p) run_cv(adv_task, p, folds = 5))
  final_fits <- lapply(adv_routes, function(p) fit_final(adv_task, p))
  list(cv = cv_results, fits = final_fits)
}, diagnostics = TRUE)
measured$reuse
```

For an execution comparison, check features/order, splits, predictions,
metrics, diagnostics and random behavior. Runtime gains depend on genuinely
shared work, not reduced candidate counts or training budgets.

### Uncertainty, stability and cohort synthesis

`sb_bootstrap_external()` quantifies external-patient sampling uncertainty
conditional on the locked model and predictions. `sb_bootstrap_calibration()`
uses the declared task and pipeline for calibration assessment. Bootstrap
draws, effective repetitions, failures and prediction horizons belong in the
analysis report. Development fitting and model-selection uncertainty are not
included in a conditional external bootstrap interval.

The installed visualization example demonstrates bootstrap intervals, paired
differences, five cohort-synthesis views, selection stability and SHAP:

```r
source(system.file("examples", "advanced-examples.R",
  package = "stormbreaker", mustWork = TRUE))
example <- stormbreaker_visualization_example()
plot_ext_ci(example$bootstraps, "roc_auc")
plot_ext_diff(example$paired, "roc_auc")
plot_forest(example$meta$Discrimination)
plot_forest(example$meta$Calibration_slope)
plot_forest(example$meta$Mean_calibration)
plot_forest(example$meta$Discrimination_difference)
plot_forest(example$meta$Prediction_error_difference)
```

Its 40 bootstrap draws are for a small synthetic demonstration. Interval plots
read saved results. Bootstrap metric names are used exactly as stored, for
example `ipcw_brier_score`, versus `ipcw_brier` in result views. Cohort intervals,
pooled intervals and between-cohort prediction intervals have different meanings.
Use feature-selection frequency and overlap to inspect saved resampling
selections. Choose `source = "glmnet_active_set"` for saved active-set records
in supported selection-frequency views.

## Multiclass and independent analysis

### Multiclass tasks

```r
multi_data <- iris
multi_data$patient_id <- paste0("sample_", seq_len(nrow(multi_data)))
multi_genes <- names(iris)[1:4]
multi_task_data <- multi_task(multi_data, "patient_id", "Species",
  predictors = multi_genes, classes = levels(iris$Species))
multi_result <- multi_run(multi_data, "Species", predictors = multi_genes,
  id = "patient_id", models = "RF", folds = 3)
multi_comparison <- sb_compare_multiclass(multi_task_data,
  c("lasso", "elastic_net", "random_forest"), folds = 5, inner_folds = 5)
multi_comparison$comparison
multi_comparison$results$lasso$metrics$per_class
plot_multi_cm(multi_comparison$results$lasso$metrics)
```

Use `models("multiclass")` and `sb_method_catalog()` for supported methods and
dependencies. Multiclass interfaces are not a Cartesian extension of the binary
selector menu. Class labels must match the data; small-class and discriminant
model assumptions still apply. Supported saved fits use
`sb_predict_multiclass(fit, new_patients)` for per-class probabilities.

### Independent method analysis

`analyze_methods()` fits each selected method independently to the same
candidate features using all supplied samples. To assess a complete
selection-and-prediction pipeline on held-out data, use route comparison
and run selection inside each training fold. Use binary, survival or
multiclass task objects as appropriate.

```r
analysis <- analyze_methods(adv_task,
  methods = c("LASSO", "Elastic Net", "random_forest"))
analysis$status
selected_genes(analysis, "LASSO")
sb_method_features(analysis, "LASSO")
sb_method_features(analysis, "random_forest")
analysis$results$lasso$plots$coefficient_path
analysis$results$lasso$plots$cv_curve
```

`selected_genes()` returns a defined selected set; `sb_method_features()` retains
full feature evidence. For multiclass results, use `class = "setosa"` or another
actual class to inspect class-specific coefficients. Use the forest
importance table to inspect feature rankings.

| Method | Retained evidence / selection rule |
|---|---|
| Boruta | Importance, hit history; confirmed features after RoughFix |
| varSelRF | Subset size and OOB error; smallest subset within one standard error |
| VSURF | Importance and stage sizes; prediction-stage subset |
| Lasso / Elastic Net | Coefficient path and CV; nonzero coefficients at selected lambda |
| MCP / SCAD | Coefficient path; minimum-BIC nonzero features |
| ABESS | Support size, GIC and selected subset |
| Random forest | Importance and tree evidence; rankings |
| mRMR, JMI, CMIM, DISR, JMIM | Mutual-information selection scores and selected features |

In independent analysis, Lasso uses `lambda.1se`; Elastic Net uses `alpha = 0.5`.
Their internal CV defaults are five folds and `inner_seed = 1701`. Method
failures and plotting failures are recorded separately. Native evidence can
include class statistics, tree splits, selection histories and support criteria;
inspect the returned method-specific tables and `plot_status`.

Two menu configurations have input-scale assumptions worth checking:
NearestCentroid uses class means, Euclidean distance and equal priors without
automatic standardization; its softmax distance score is not separately
calibrated risk. BernoulliNB uses `x > 0`, Laplace smoothing 1 and training
priors. Positive-only expression data can lose quantitative information under
that fixed presence/absence rule.

## Runnable examples

Start with the [Quick start](https://github.com/omics-odysseus/stormbreaker#quick-start)
to run the binary route workflow and the Cox survival example together on
packaged synthetic data.

For the three applications in the manuscript, see the
[study examples](paper-cases.md). Use the scripts below for individual tasks.

| Script | Workflow | Main outputs |
|---|---|---|
| `01-binary.R` | 120-row synthetic binary development cohort | CV metrics, ROC, calibration and decision curves |
| `02-compare.R` | Specified routes on shared folds | Comparison table, performance plot and ROC |
| `03-external.R` | Development and two validation cohorts | One final model, predictions, metrics and ROC |
| `04-survival.R` | Development and validation, with month units | Risk, horizon probabilities, concordance, timed AUC, Brier, calibration and DCA |
| `05-diagnostics.R` | Lasso followed by Logistic | Selected genes, native diagnostics and coefficients |
| `06-shap.R` | Saved model and eight query samples | SHAP values, importance, beeswarm, dependence and waterfall |
| `07-route-workflow.R` | Direct Logistic versus Lasso followed by Logistic; shared comparison, final fitting and reuse | Route AUCs, chosen route, test ROC and predictions after reloading |

Each script returns `example_result` and `figures`. Replace the filename to run
another workflow with the packaged synthetic CSVs. Install the required engines
before running a script; export the returned figures separately as needed.

Source the installed `examples/advanced-examples.R`, list the demonstrations
with `advanced_examples()`, and run one with
`run_advanced_example("multiclass")` or another listed name. Multi-cohort,
custom-model and paper-figure examples use the same entry point.
The `user-models` entry is a file-based template: replace its two input paths
with your own gene list and development CSV before running it.

For a file-based develop/finalize workflow, the same script provides
`workflow_demo()`, `workflow_develop()`
and `workflow_finalize()`. The latter requires an explicit selected candidate
and selection reason, with external predictor and outcome files supplied
separately. The corresponding `development.rds` and final model can be reopened
for later prediction. `stormbreaker_visualization_example()` runs the gallery
above. Call `regenerate_demo_data(directory)` to create synthetic CSVs in a
new directory.

Short names call the corresponding public functions.
Use `stormbreaker::fit_final()` when another package has an overlapping name.
The [workflow name map](../guides/preferred-name-map.csv) and
[diagnostic name map](../guides/specialty-name-map.csv) list all aliases.
If vignettes were built during installation, `vignette("stormbreaker-short-scenarios")`
renders this same English guide. The GitHub installation command above skips
the build; use the Markdown guide directly or install with `build_vignettes = TRUE`.
Use `citation("stormbreaker")` for the software citation.

### Where older tutorials moved

| Old vignette name | Section in this guide |
| --- | --- |
| `stormbreaker-calibration-dca` | [Reading metrics](#reading-metrics) and [Diagnostics and explanations](#diagnostics-and-explanations) |
| `stormbreaker-external-result-views` | [Final fitting and external evaluation](#final-fitting-and-external-evaluation) |
| `stormbreaker-html-report` | [Reports and result browsing](#reports-and-result-browsing) |
| `stormbreaker-result-views` | [Other saved-result views](#other-saved-result-views) |
| `stormbreaker-visualization` | [Diagnostics and explanations](#diagnostics-and-explanations) |
| `stormbreaker-workflow` | [Model comparison](#model-comparison) and [Object handoffs](#object-handoffs) |

## Parameter reference

The following reference collects argument types, defaults and object handoffs.
Use `?function_name` for the full help of an individual function.

### Object handoffs

```text
data -> bi_run / surv_run -> run$comparison -> result_table
                         -> finish_run -> finished$fit -> predict_final
                                                       -> plots
                                                       -> model_shap
finished$evaluations -> clinical_plots
finished$plots -> export_run
finished$fit -> saveRDS -> readRDS -> prediction and diagnostics
```

`run` holds the comparison, `finished` holds the selected pipeline's results,
and `finished$fit` is the model for prediction and explanation. `finish_run()`
uses a name from `names(run$pipelines)`, distinct from the internal
`pipeline_id` in metric tables.

### Data, examples and menus

| Function | Parameter | Type and default | Meaning / return |
|---|---|---|---|
| `bi_example()` | `n`, `seed` | Integers; `60`, `20260818` | At least 10 synthetic samples; returns ID, outcome and three genes |
| `surv_example()` | `n`, `seed` | Integers; `60`, `20260812` | At least 10 synthetic samples; returns ID, time, event and three genes |
| `models()` | `task` | String; `"binary"` | Also `"survival"` or `"multiclass"`; returns names, configurations and dependency status |

Example seeds are positive integers. The binary menu additionally reports
`can_predict`, `can_extract_features` and `stage_selector`, distinguishing
prediction methods from selection-only methods.

See [Multiclass and independent analysis](#multiclass-and-independent-analysis) for multiclass workflows.

### `bi_run()` and `surv_run()`

#### Shared arguments

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `data` | Data frame or corresponding task object | Required | Development observations; an existing task supplies its saved columns and settings |
| `predictors` | Character vector or `NULL` | `NULL` | All columns other than ID and endpoints; specify explicitly when other descriptive fields exist |
| `id` | Column-name string or `NULL` | `NULL` | Uses the supplied ID column or creates `.sb_row_id` |
| `models` | Character vector or `NULL` | `NULL` | Full menu by default; specified model names default to single-model comparisons |
| `stages` | Integer vector | Conditional, below | Single-stage, two-stage and supported three-stage pipelines |
| `selector` | Selector-chain object or `NULL` | `NULL` | Survival accepts `surv_selector()`; new binary code uses ordered names in `bi_pipeline()` |
| `folds` | Integer | `5` | Outer cross-validation folds |
| `seed` | Positive integer | `1` | Outer resampling seed |
| `plots` | `NULL`, `TRUE` or `FALSE` | `NULL` | Controls overview figures and additional final fits, below |
| `pipelines` | Uniquely named list or `NULL` | `NULL` | Supply complete routes here instead of `models` and `selector`. Use `stages` for menu-generated routes. |
| `decision_thresholds` | Numeric probability vector or `NULL` | `NULL` | Declared thresholds for saved decision curves; none are requested by default |

`stages` defaults to `c(1, 2)` for full binary runs and `c(1, 2, 3)` for full
survival runs, whose three-stage options use built-in presets. When `models` is
specified, `stages` defaults to `1`.

`plots = NULL` draws comparison overviews and additionally fits a final model
only when exactly one pipeline is requested. `plots = TRUE` fits and draws every
candidate. `plots = FALSE` skips plots and additional final fits. After a batch
comparison, use `finish_run()` for the selected pipeline. Detailed all-pipeline
plots add runtime, memory and output size.

#### Endpoint arguments

| Function | Argument | Type and default | Meaning |
|---|---|---|---|
| `bi_run()` | `outcome` | Column-name string; `NULL` | Supply the 0/1 outcome column for a data-frame input |
| `surv_run()` | `time`, `event` | Column-name strings; `NULL` | Supply follow-up time and 0/1 event columns |
| `surv_run()` | `time_unit` | String; `"month"` | Unit of the supplied follow-up times and prediction horizon |
| `surv_run()` | `prediction_horizon` | Positive number or `NULL`; `NULL` | Risk ranking by default; `6` requests a six-month probability target |

Set `time_unit` to the unit of the supplied values. For fixed-horizon evaluation,
supply model probabilities and cohort follow-up for that horizon.

#### Returned fields

`$comparison` contains CV results for `result_table()`; `$pipelines` contains
requested pipelines; `$plots` contains overview figures; `$details` contains
any additional final fits and figures; `$issues` records problems; and
`$development_task` retains the development data needed by `finish_run()`.

### Pipeline constructors

- `bi_pipeline(model = "RF", selector = NULL)`: `model` is a string or ordered
  character vector. In `c("Lasso", "Logistic")`, earlier stages select features
  and the last predicts. `selector` is a compatibility argument; use ordered
  model names in new binary code.
- `surv_pipeline(model = "Cox", selector = NULL)`: one final model string and an
  optional survival selector specification.
- `surv_selector(methods)`: required ordered character vector of selector names,
  for example `c("Lasso", "StepCox")`.

These functions create pipeline specifications. Put them in a named list and
pass it as `pipelines` to the appropriate run function for fitting.

### `finish_run()`

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `run` | Binary/survival comparison result | Required | Completed development comparison |
| `pipeline` | String | Required | Exact selected name from `names(run$pipelines)` |
| `external` | Uniquely named list of data frames | `list()` | Cohorts with consistent ID, predictor and endpoint columns |
| `plots` | Logical | `TRUE` | Algorithm and external figures; captures diagnostic materials when a new fit is required |
| `explain` | External cohort name or `NULL` | `NULL` | Optional cohort to explain using development background data |

`finish_run()` keeps the original comparison and reuses an available final fit.
Otherwise it fits the selected pipeline on all development data. Capture the
diagnostic materials you plan to use when fitting the model.

Returns `$fit`, `$genes`, `$predictions`, `$evaluations`, `$metrics`, `$plots`,
`$explanation` and `$issues`. Use `$evaluated_predictions` for the predictions
used in evaluation and `$predictions` for the original model predictions.
Survival `$risk_groups` includes the development cutoff, group counts and
available log-rank results.

### `result_table()`

| Argument | Type and default | Meaning |
|---|---|---|
| `x` | Required CV comparison/result or external evaluation | Usually `run$comparison` or `finished$evaluations$Validation` |
| `component` | String; `"metrics"` | Also `"predictions"` or `"runs"` |
| `metric` | String or `NULL`; `NULL` | Explicitly required for metrics; leave `NULL` for the other components |
| `cohort_label` | String or `NULL`; `NULL` | Optional display label |

Binary metrics include `roc_auc`, `average_precision`, `brier_score` and
`log_loss`. Survival metrics include `harrell_c` and, when requested and
evaluable, `time_auc`, `uno_c` and `ipcw_brier`. The returned data frame reports
metric `status`, overall `run_status`, and the saved `aggregation` convention.

### `clinical_plots()`

| Argument | Type and default | Meaning |
|---|---|---|
| `x` | Required saved evaluation/result | A run displays development CV; a completed run displays external evaluations |
| `detail` | String; `"brief"` | Brief requests ROC, calibration and DCA; `"full"` adds PR and probability distributions |
| `thresholds` | Numeric vector in [0,1] or `NULL`; `NULL` | Adds binary threshold-metric plots |
| `confusion_threshold` | One probability or `NULL`; `NULL` | Adds a confusion matrix; must occur in `thresholds` |
| `native` | Native-diagnostic result or named list; `list()` | Already extracted model diagnostics |
| `selectors` | Named selector-diagnostic list; `list()` | Already extracted selector diagnostics |
| `shap` | SHAP result or `NULL`; `NULL` | Adds saved importance and beeswarm views |
| `feature_values` | Data frame or `NULL`; `NULL` | Values for SHAP colors/dependence |
| `id_col` | Column-name string or `NULL`; `NULL` | Explicit ID matching for the feature values |
| `dependence_features` | Character vector or `NULL`; `NULL` | Named SHAP features; requires SHAP, values and ID |
| `nomogram` | Nomogram result or `NULL`; `NULL` | Already constructed exact-model nomogram |

Returns `$plots`, `$status`, `$diagnostics`, `$captions` and grouping information.
`names(views$plots)` lists figures, and `views$status` records availability.
Threshold metrics and confusion matrices use saved binary predictions. Decision
curves must be calculated during modeling with `decision_thresholds`. Survival
results use saved horizon-specific calibration and decision curves.

### `plots()`

`plots(fit, new_data = NULL, id_col = NULL, background_data = NULL)` takes a
successful final fit. Optional `new_data` supplies query rows for local/model
diagnostics. Supplying `background_data` additionally requests SHAP and requires
explanation rows and their ID column. It returns `$plots`, `$pages`,
`$unavailable`, `$native` and `$shap`. `plots(fit)` reads saved diagnostics;
`plots(fit, new_data = query)` supplies samples for query-dependent views.

By default, `finish_run()` places algorithm figures in `$plots$algorithm` for
`export_run()`. Saving `$fit` also retains its captured diagnostic materials.
See [Diagnostics and explanations](#diagnostics-and-explanations).

### SHAP

| `model_shap()` argument | Type and default | Meaning |
|---|---|---|
| `fit` | Required final fit | Saved binary or survival model |
| `new_data` | Required data frame | Samples to explain, with IDs and required features |
| `id_col` | Required string | ID column |
| `background_data` | Required data frame | Reference features, normally from development |
| `permutations` | Positive integer; `32` | Random feature orders |
| `background_n` | Positive integer; `20` | Maximum sampled background rows |
| `seed` | Positive integer; `1` | Permutation and background-sampling seed |

Returns `$contributions`, `$reference_value`, `$scores` and method metadata.
This is approximate permutation SHAP: binary explanations use positive-class
probability and survival explanations use the saved risk score.

- `plot_shap_imp(x, top_n = 20, groupby = NULL)`: required SHAP result, positive
  integer display limit and optional per-sample grouping vector. Named groups
  are aligned by ID. Returns mean-absolute-contribution plots.
- `plot_shap(x, top_n = 20, feature_values = NULL, id_col = NULL)`: required SHAP
  result and display limit; optional values plus ID enable feature-value colors.
- `plot_shap_one(x, row = 1, top_n = 10)`: required SHAP result, one-based sample
  index and positive integer feature limit. Returns a waterfall, including for
  a single feature or a single observation; undisplayed terms are summed.

These return ggplot objects from the supplied SHAP result. Use
`sb_plot_shap_heatmap(shap, top_n = 20)` for signed contributions by observation.

### Prediction, export and reload

`predict_final(fit, new_data, id_col)` requires a saved final fit, a data frame
with IDs and required predictors, and the ID-column string. Predictions are in
`$predictions`: probabilities for binary models; risk scores and supported
horizon predictions for survival models.

`export_run(x, directory)` requires a completed run and a **new** directory
string. It invisibly returns the absolute directory path and writes:

- `model-comparison.csv`, `genes-by-stage.csv` and `cohorts.csv`;
- cohort predictions, evaluated predictions and available metric tables;
- `figure-*.pdf` by default, or PNG with `figure_format = "png", dpi = 300L`;
  the `figures.csv` name index and `figure-notes.csv`;
- issue records and, when supplied, `explanation.rds`.

`x$plots` supplies the figures to export; add separately created clinical and
SHAP figures there first. `font_family = "Times New Roman"` is optional and
requires the font on the export machine. Explicitly mapped fonts or embedded
labels in custom objects are not rewritten.
Save the model with `saveRDS(x$fit, ...)`, or the whole
result with `saveRDS(x, ...)`. Saved objects may contain sample-level data.

Select three pipelines from the completed `binary` comparison and export their
models and figures to new directories:

```r
chosen <- c("Logistic", "RF", "SVM")
for (i in seq_along(chosen)) {
  item <- finish_run(binary, chosen[i])
  destination <- file.path("algorithm-results", sprintf("model-%03d", i))
  export_run(item, destination)
  saveRDS(item$fit, file.path(destination, "model.rds"))
  writeLines(chosen[i], file.path(destination, "model-name.txt"))
}
```

Choose the needed subset of `names(binary$pipelines)`. Full combinations can
contain thousands of pipelines, so choose the export scope to match the study.

### Troubleshooting

| Question | Where to look |
|---|---|
| Missing model metrics | `result_table(run$comparison, "runs")` and `run$issues` |
| Unknown pipeline name | `names(run$pipelines)` |
| Missing decision curve | Original `decision_thresholds`; survival horizon, probability support and follow-up |
| Missing algorithm figure | `plots(fit)$unavailable`, captured evidence and explicit query rows |
| Separately drawn figure not exported | Whether it was added to `finished$plots` |
| Existing output directory | Select a new directory |
| Full R help | `?bi_run`, `?surv_run`, `?finish_run`, `?clinical_plots`, `?plots`, `?model_shap` |
