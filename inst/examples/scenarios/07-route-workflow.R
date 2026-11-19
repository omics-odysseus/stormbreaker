## ---- scenario-route-workflow ----
library(stormbreaker)

read_demo <- function(name) {
  read.csv(system.file("extdata", name, package = "stormbreaker",
    mustWork = TRUE), colClasses = c(patient_id = "character"))
}
development <- read_demo("demo-binary-development.csv")
genes <- c("gene_a", "gene_b", "gene_c")

# 1. Define two routes with the same predictor.
routes <- list(
  Direct = bi_pipeline("Logistic"),
  Lasso = sb_pipeline("lasso_binary_nonzero", "logistic",
    selector_options = list(lambda_rule = "lambda.min"))
)

# 2. Compare on shared development folds; ties keep the first listed route.
run <- bi_run(development, outcome = "outcome", id = "patient_id",
  predictors = genes, pipelines = routes, folds = 3, seed = 42,
  plots = FALSE)
route_ids <- vapply(routes, \(route) route$id, character(1))
cv_auc <- result_table(run$comparison, metric = "roc_auc")
cv_auc <- cv_auc[match(route_ids, cv_auc$pipeline_id), , drop = FALSE]
cv_auc$route <- names(route_ids)
selected <- cv_auc$route[which.max(cv_auc$estimate)]

# 3. Refit the chosen route and evaluate it in a separate synthetic cohort.
test <- read_demo("demo-binary-external-a.csv")
finished <- finish_run(run, selected, external = list(Test = test), plots = FALSE)
test_auc <- result_table(finished$evaluations$Test, metric = "roc_auc")
figures <- list(ROC = plot_roc(finished$evaluations$Test))

# 4. Reload the fitted model and predict without test outcomes.
model_file <- tempfile(fileext = ".rds")
saveRDS(finished$fit, model_file)
reloaded <- readRDS(model_file)
test_predictors <- test[c("patient_id", genes)]
prediction <- predict_final(reloaded, test_predictors, "patient_id")

example_result <- list(routes = routes, run = run, cv_auc = cv_auc,
  selected = selected, finished = finished, test_auc = test_auc,
  reloaded = reloaded, test_predictors = test_predictors,
  prediction = prediction, figures = figures)

if (interactive()) print(figures$ROC)
