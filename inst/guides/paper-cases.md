# Study examples

Prepare the inputs below, then fit and apply the selected route for your case.
To build and compare your own candidate routes, use the
[README quick start](../../README.md#quick-start),
[runnable example](../examples/scenarios/07-route-workflow.R) or the
[user guide](user-guide.md).

## Prepared inputs

Provide one row per sample, a unique `analysis_id` for each row, a `patient_id`
to group repeated or paired records in development, and numeric expression
columns named by the case's gene vector. Survival data also require `time`
(days in this example) and `event` (1 = event, 0 = censored); binary data require
`outcome` (1 = positive, 0 = negative). Apply the same feature definitions and
scales in development and each comparison cohort. Each `*_comparison`
data frame contains one prepared cohort. Use Supplementary Tables S1–S3
for the study's configurations, gene lists, cohort sources and preprocessing.

Run this helper once. It fits the selected route, saves and reloads the model,
then predicts using IDs and retained features. Replace `tempfile()` with a new
file path to keep the fitted model. Use that saved fit with `predict_final()`
for additional cohorts.

```r
library(stormbreaker)

fit_and_reuse <- function(task, route, new_data, model_file) {
  fit <- fit_final(task, route)
  if (fit$status != "success") stop("The chosen route did not fit.")
  saveRDS(fit, model_file)
  reloaded <- readRDS(model_file)
  columns <- c("analysis_id", reloaded$required_features)
  predict_final(reloaded, new_data[columns], "analysis_id")
}
```

The study used the same external comparison cohorts to select configurations
and report performance.

## HNSC survival prediction

- **Prepare:** Use the 518 TCGA-HNSC training patients in
  `survival_development`, the 118-gene panel in `survival_genes`, and one
  comparison cohort in `survival_comparison`. The study applied the model to
  429 patients across nine comparison cohorts. Supply survival time in days
  and event status; install `glmnet` and `pamr` for this route.
- **Compare and choose:** Of the 1,718 candidate routes, 1,706 trained
  successfully. The study selected from routes evaluated across the
  comparison cohorts by equally weighted mean external C-index.
- **Selected route:** Ridge-Cox weight-based selection followed by prediction
  analysis for microarrays (PAMR), study configuration `aligned_1112`. The
  selector uses five-fold inner cross-validation and ranks genes by absolute
  coefficient multiplied by training-set standard deviation. The study fit
  retained 59 genes.

**Fit, apply and inspect:** View the C-index and patient-level predictions.
See Figures 3, 4 and 7 for the study results.

```r
survival_route <- surv_pipeline("pamr", selector = surv_selector("Ridge"))
survival_task <- surv_task(survival_development, "analysis_id", "time", "event",
  survival_genes, time_unit = "day", group_col = "patient_id")
survival_prediction <- fit_and_reuse(survival_task, survival_route,
  survival_comparison, tempfile("hnsc-survival-", fileext = ".rds"))
survival_endpoint <- surv_endpoint(survival_comparison, "analysis_id",
  "time", "event", time_unit = "day")
survival_evaluation <- eval_external(survival_prediction, survival_endpoint)
result_table(survival_evaluation, metric = "harrell_c")
head(survival_prediction$predictions)
```

## Tumour-normal tissue classification

- **Prepare:** Put the 564 TCGA-HNSC training specimens (520 tumour, 44 normal)
  in `tissue_development`, the 118-gene panel in `tissue_genes`, and one
  comparison cohort in `tissue_comparison`. Encode tumour as `1` and
  normal or source-labelled non-malignant tissue as `0`. Give paired
  specimens distinct sample IDs and the same patient ID. The study applied the
  model to 656 specimens across nine comparison cohorts. For the two
  input-column exceptions, use the preparation rules in Supplementary Table S3.
  Install `praznik` and `xgboost` for this route.
- **Compare and choose:** The study considered 2,340 candidates: 92 with no
  additional selection, 92 with predefined gene sets, and 2,156 with a
  learned selection stage. It selected by equally weighted mean external
  ROC AUC.
- **Selected route:** Minimum-redundancy maximum-relevance (mRMR) selection
  followed by linear XGBoost, study configuration `binary_1353`. The
  predictor uses 100 boosting rounds and a learning rate of 0.1. The study
  fit retained 89 genes.

**Fit, apply and inspect:** View ROC AUC and specimen-level predictions.
See Figure 5 for the study results.

```r
tissue_route <- bi_pipeline(c("mRMR", "xgboost_linear"))
tissue_task <- bi_task(tissue_development, "analysis_id", "outcome",
  tissue_genes, group_col = "patient_id")
tissue_prediction <- fit_and_reuse(tissue_task, tissue_route, tissue_comparison,
  tempfile("tumour-normal-", fileext = ".rds"))
tissue_endpoint <- bi_endpoint(tissue_comparison, "analysis_id", "outcome")
tissue_evaluation <- eval_external(tissue_prediction, tissue_endpoint)
result_table(tissue_evaluation, metric = "roc_auc")
head(tissue_prediction$predictions)
```

## Multicancer immunotherapy-response prediction

- **Prepare:** Put the 697 training patients in `response_development`.
  The study split 872 patients from four source cohorts, stratifying by cohort
  and response, and used the other 175 for internal ranking. Put one of the ten
  external comparison cohorts (450 patients in total) in `response_comparison`.
  Encode response as `1` and nonresponse as `0`. Use the same 41 measurable
  genes from the 118-gene panel in `response_genes` and transform each row
  to within-sample midrank percentiles `(rank - 0.5) / p`, where `p` is the
  number of genes ranked. Install `randomForest` and `mclust` for this route.
- **Compare and choose:** The study used the same 2,340 binary candidates and
  selected by equally weighted mean external ROC AUC.
- **Selected route:** Random-forest screening followed by `mclustDA`, study
  configuration `binary_1585`. BIC selects among one to three mixture
  components. The study fit retained 33 genes.

**Fit, apply and inspect:** View ROC AUC and patient-level predictions.
See Figure 6 for the study results.

```r
response_route <- bi_pipeline(c("RF", "mclustda_bic"))
response_task <- bi_task(response_development, "analysis_id", "outcome",
  response_genes, group_col = "patient_id")
response_prediction <- fit_and_reuse(response_task, response_route,
  response_comparison, tempfile("immunotherapy-", fileext = ".rds"))
response_endpoint <- bi_endpoint(response_comparison, "analysis_id", "outcome")
response_evaluation <- eval_external(response_prediction, response_endpoint)
result_table(response_evaluation, metric = "roc_auc")
head(response_prediction$predictions)
```
