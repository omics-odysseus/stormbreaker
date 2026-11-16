# ---- shap-prediction-reuse ----
local({
test_that("SHAP validates the saved fit once without changing predictions", {
  d <- sb_example_binary(60, 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, bi_pipeline("Logistic"))
  validate <- sb_validate_final_fit_for_prediction
  calls <- 0L
  local_mocked_bindings(sb_validate_final_fit_for_prediction = function(fit) {
    calls <<- calls + 1L
    validate(fit)
  })

  result <- sb_model_shap(fit, d[1:2, ], "patient_id", d,
    permutations = 2, background_n = 3, seed = 12)

  expect_identical(calls, 1L)
  expected <- sb_predict_final(fit, d[1:2, ], "patient_id")$predictions$probability
  expect_identical(result$scores$score, expected)
  expect_equal(unname(result$reference_value + rowSums(result$contributions)), expected)
  expect_identical(result, model_shap(fit, d[1:2, ], "patient_id", d,
    permutations = 2, background_n = 3, seed = 12))
})

test_that("tree SHAP does not validate a second time through the generic wrapper", {
  skip_if_not_installed("rpart")
  d <- sb_example_binary(60, 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "decision_tree_binary"))
  expect_identical(fit$status, "success")
  validate <- sb_validate_final_fit_for_prediction
  calls <- 0L
  local_mocked_bindings(sb_validate_final_fit_for_prediction = function(fit) {
    calls <<- calls + 1L
    validate(fit)
  })

  result <- sb_tree_shap(fit, d[1:2, ], "patient_id", d,
    permutations = 2, background_n = 3)

  expect_identical(calls, 1L)
  expect_s3_class(result, "stormbreaker_tree_shap")
  expect_equal(unname(result$reference_value + rowSums(result$contributions)), result$scores$score)
})

test_that("the SHAP inner predictor still rejects invalid probability values", {
  d <- sb_example_binary(60, 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, bi_pipeline("Logistic"))
  local_mocked_bindings(sb_predict_final_state = function(state, x) {
    list(probability = rep(1.1, nrow(x)))
  })

  error <- tryCatch(sb_model_shap(fit, d[1:2, ], "patient_id", d,
    permutations = 2, background_n = 3), error = identity)

  expect_s3_class(error, "stormbreaker_invalid_final_probability")
})

test_that("survival SHAP explains the saved risk and preserves patient order", {
  d <- sb_example_survival(60, 42)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b"), "month")
  fit <- sb_fit_final(task, surv_pipeline("Cox"))
  query <- d[c(4, 1, 3), ]

  result <- sb_model_shap(fit, query, "patient_id", d,
    permutations = 2, background_n = 3)

  expect_identical(result$scores$patient_id, query$patient_id)
  expect_identical(result$scores$score,
    sb_predict_final(fit, query, "patient_id")$predictions$risk)
  expect_equal(unname(result$reference_value + rowSums(result$contributions)), result$scores$score)
  expect_identical(result$score_scale, "saved_risk_score")
})

test_that("survival SHAP still rejects non-finite inner predictions", {
  d <- sb_example_survival(60, 42)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b"), "month")
  fit <- sb_fit_final(task, surv_pipeline("Cox"))
  local_mocked_bindings(sb_predict_final_state = function(state, x) {
    list(risk = rep(Inf, nrow(x)))
  })

  error <- tryCatch(sb_model_shap(fit, d[1:2, ], "patient_id", d,
    permutations = 2, background_n = 3), error = identity)

  expect_s3_class(error, "stormbreaker_invalid_final_risk")
})
})

# ---- logrank-status ----
local({
test_that("log-rank status distinguishes a successful test from unavailable groups", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  d <- sb_example_survival(60, 42)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b"), "month")
  fit <- sb_fit_final(task, surv_pipeline("Cox"))
  reference <- surv_risk_reference(fit, task)

  result <- surv_risk_plots(reference, task)

  native <- survival::survdiff(survival::Surv(time, event) ~ group, result$patients)
  expect_identical(result$logrank_status, "success")
  expect_identical(result$logrank_error, NA_character_)
  expect_equal(result$logrank_p, stats::pchisq(native$chisq, 1, lower.tail = FALSE))
  reference$cutoff <- Inf
  single <- surv_risk_plots(reference, task)
  expect_identical(single$status, "single_group")
  expect_identical(single$logrank_status, "single_group")
  expect_identical(single$logrank_p, NA_real_)
})

test_that("log-rank reports no events without inventing a P value", {
  d <- data.frame(time = 1:6, event = 0L,
    group = factor(rep(c("Low", "High"), each = 3)))

  result <- sb_risk_logrank(d)

  expect_identical(result$status, "no_events")
  expect_identical(result$p, NA_real_)
  expect_identical(result$error, NA_character_)
})

test_that("a log-rank error remains visible without discarding usable panels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  d <- sb_example_survival(60, 42)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b"), "month")
  fit <- sb_fit_final(task, surv_pipeline("Cox"))
  reference <- surv_risk_reference(fit, task)
  local_mocked_bindings(survdiff = function(...) stop("SYNTHETIC_LOGRANK_FAILURE"),
    .package = "survival")
  warnings <- character()

  result <- withCallingHandlers(surv_risk_plots(reference, task), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })

  expect_identical(result$logrank_status, "failed")
  expect_identical(result$logrank_error, "SYNTHETIC_LOGRANK_FAILURE")
  expect_identical(result$logrank_p, NA_real_)
  expect_match(warnings, "SYNTHETIC_LOGRANK_FAILURE")
  expect_s3_class(result$plot, "patchwork")
  expect_equal(nrow(result$patients), nrow(d))
})

test_that("completion retains and exports a failed log-rank calculation", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  d <- sb_example_survival(60, 42)
  run <- surv_run(d, "time", "event", c("gene_a", "gene_b"),
    "patient_id", models = "Cox", folds = 3, plots = FALSE)
  external <- d
  external$patient_id <- paste0("EXT_", external$patient_id)
  local_mocked_bindings(survdiff = function(...) stop("SYNTHETIC_LOGRANK_FAILURE"),
    .package = "survival")
  warnings <- character()

  result <- withCallingHandlers(sb_complete_run(run, "Cox", list(Validation = external)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  expect_identical(result$evaluations$Validation$status, "success")
  expect_identical(result$risk_groups$Validation$logrank_status, "failed")
  expect_match(result$issues$risk$Validation, "SYNTHETIC_LOGRANK_FAILURE")
  expect_match(warnings, "SYNTHETIC_LOGRANK_FAILURE")
  expect_s3_class(result$plots$external$Validation, "patchwork")
  output <- file.path(withr::local_tempdir(), "export")
  sb_export_run(result, output)
  exported <- read.csv(file.path(output, "risk-groups.csv"))
  expect_identical(exported$logrank_status, "failed")
  expect_identical(exported$logrank_error, "SYNTHETIC_LOGRANK_FAILURE")
  expect_match(paste(read.csv(file.path(output, "issues.csv"))$message, collapse = " "),
    "SYNTHETIC_LOGRANK_FAILURE")
})

test_that("completion does not claim a log-rank test when plots were not requested", {
  d <- sb_example_survival(60, 42)
  run <- surv_run(d, "time", "event", c("gene_a", "gene_b"),
    "patient_id", models = "Cox", folds = 3, plots = FALSE)
  external <- d
  external$patient_id <- paste0("EXT_", external$patient_id)
  local_mocked_bindings(survdiff = function(...) stop("Unexpected log-rank test"),
    .package = "survival")

  result <- sb_complete_run(run, "Cox", list(Validation = external), plots = FALSE)

  expect_identical(result$evaluations$Validation$status, "success")
  expect_identical(result$risk_groups$Validation$logrank_status, "not_computed")
  expect_identical(result$risk_groups$Validation$logrank_p, NA_real_)
  expect_identical(result$risk_groups$Validation$logrank_error, NA_character_)
})
})

# ---- short-scenarios ----
local({
run_short_scenario <- function(file) {
  environment <- new.env(parent = globalenv())
  sys.source(system.file("examples", "scenarios", file,
    package = "stormbreaker", mustWork = TRUE), envir = environment)
  environment$example_result
}

test_that("binary scenario connects saved CV tables and three clinical plots", {
  skip_if_not_installed("ggplot2")
  x <- run_short_scenario("01-binary.R")
  expect_equal(nrow(x$input), 120L)
  expect_equal(nrow(x$metrics), 1L)
  expect_identical(x$metrics$status, "success")
  expect_setequal(names(x$figures), c("ROC", "Calibration", "DCA"))
  expect_gt(nrow(x$figures$ROC$data), 0L)
  expect_identical(x$metrics$estimate, x$comparison$results[[1]]$metrics$roc_auc$estimate)
})

test_that("comparison scenario returns both declared candidates", {
  skip_if_not_installed("ggplot2")
  x <- run_short_scenario("02-compare.R")
  expect_equal(nrow(x$metrics), 2L)
  expect_identical(x$metrics$status, rep("success", 2))
  expect_length(x$comparison$results, 2L)
  expect_s3_class(x$figures$Performance, "ggplot")
  expect_gt(nrow(x$figures$ROC$data), 0L)
})

test_that("external scenario keeps both evaluations and the original predictions", {
  skip_if_not_installed("ggplot2")
  x <- run_short_scenario("03-external.R")
  expect_identical(x$finished$fit$status, "success")
  expect_equal(nrow(x$finished$predictions$Validation_A), 80L)
  expect_equal(nrow(x$finished$predictions$Validation_B), 70L)
  expect_identical(x$metrics$status, rep("success", 2))
  expect_setequal(x$metrics$cohort_label, c("Validation_A", "Validation_B"))
  expect_length(x$finished$direction_correction, 0L)
  expect_null(x$finished$evaluations$Validation_A$direction_correction)
  expect_null(x$finished$evaluations$Validation_B$direction_correction)
  original_a <- predict_final(x$finished$fit, x$external$Validation_A,
    "patient_id")$predictions
  original_b <- predict_final(x$finished$fit, x$external$Validation_B,
    "patient_id")$predictions
  expect_identical(x$finished$predictions$Validation_A, original_a)
  expect_identical(x$finished$predictions$Validation_B, original_b)
  expect_identical(x$finished$predictions$Validation_A$probability,
    x$finished$evaluated_predictions$Validation_A$probability)
  expect_identical(x$finished$predictions$Validation_B$probability,
    x$finished$evaluated_predictions$Validation_B$probability)
  expect_s3_class(x$figures$ROC, "ggplot")
})

test_that("survival scenario keeps horizon units and complete external rows", {
  skip_if_not_installed("ggplot2")
  x <- run_short_scenario("04-survival.R")
  p <- x$prediction$predictions
  expect_identical(x$fit$status, "success")
  expect_identical(x$evaluation$status, "success")
  expect_equal(nrow(p), 100L)
  expect_equal(unique(p$prediction_horizon), 6)
  expect_identical(unique(p$time_unit), "month")
  expect_gte(min(p$survival_probability), 0)
  expect_lte(max(p$survival_probability), 1)
  expect_setequal(x$horizon_metrics$metric, c("time_auc", "ipcw_brier"))
  expect_equal(unique(x$horizon_metrics$prediction_horizon), 6)
  expect_null(x$evaluation$direction_correction)
  expect_gt(nrow(x$figures$Calibration$data), 0L)
})

test_that("selector scenario retains native CV and coefficient evidence", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("glmnet")
  x <- run_short_scenario("05-diagnostics.R")
  expect_identical(x$fit$status, "success")
  expect_setequal(x$genes$gene, x$fit$required_features)
  expect_contains(c("gene_a", "gene_b", "gene_c"), x$genes$gene)
  expect_gt(nrow(x$genes), 0L)
  expect_contains(names(x$figures), c("selector_1_cv", "selector_1_coefficient_path", "Coefficients"))
  expect_s3_class(x$figures$selector_1_cv, "ggplot")
  expect_length(x$diagnostics$unavailable, 0L)
})

test_that("SHAP scenario reconstructs probabilities from its declared background", {
  skip_if_not_installed("ggplot2")
  x <- run_short_scenario("06-shap.R")
  expect_identical(dim(x$shap$contributions), c(8L, 3L))
  expect_identical(rownames(x$shap$contributions), x$query$patient_id)
  expect_identical(x$shap$score_scale, "class_1_probability")
  expect_equal(unname(x$shap$reference_value + rowSums(x$shap$contributions)),
    x$shap$scores$score, tolerance = 1e-10)
  expect_setequal(names(x$figures), c("Importance", "Beeswarm", "Dependence", "Waterfall"))
  expect_s3_class(x$figures$Waterfall, "ggplot")
})

test_that("tutorial code can be extracted and replayed independently", {
  skip_if_not_installed("knitr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("glmnet")
  input <- test_path("..", "..", "vignettes", "stormbreaker-short-scenarios.Rmd")
  if (!file.exists(input)) {
    input <- system.file("doc", "stormbreaker-short-scenarios.Rmd",
      package = "stormbreaker", mustWork = TRUE)
  }
  directory <- withr::local_tempdir()
  output <- file.path(directory, "tutorial.R")
  knitr::purl(input, output = output, quiet = TRUE)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    paste0("root <- ", paste(deparse(find.package("stormbreaker")), collapse = "")),
    "if (file.exists(file.path(root, 'Meta', 'package.rds'))) {",
    "  library(stormbreaker, lib.loc = dirname(root))",
    "} else {",
    "  pkgload::load_all(root, helpers = FALSE)",
    "}",
    "set.seed(42)",
    "grDevices::pdf('tutorial.pdf')",
    "source('tutorial.R', echo = FALSE, print.eval = TRUE)",
    "grDevices::dev.off()",
    "saveRDS(list(metrics = metrics, contributions = shap$contributions,",
    "  reference_value = shap$reference_value, scores = shap$scores$score), 'result.rds')"
  ), file.path(directory, "replay.R"), useBytes = TRUE)
  withr::local_dir(directory)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "replay.R"),
    stdout = "replay.log", stderr = "replay.log")
  expect_identical(status, 0L,
    info = paste(tail(readLines("replay.log", warn = FALSE), 20), collapse = "\n"))
  if (status != 0L) return(invisible(NULL))
  result <- readRDS("result.rds")
  expect_s3_class(result$metrics, "data.frame")
  expect_identical(dim(result$contributions), c(8L, 3L))
  expect_equal(
    unname(result$reference_value + rowSums(result$contributions)),
    result$scores, tolerance = 1e-10)
})
})

# ---- read-cohort ----
local({
csv_fixture <- function(lines) {
  path <- tempfile(fileext = ".csv")
  writeLines(lines, path)
  path
}

test_that("CSV cohort joins by preserved IDs, not row position", {
  expression <- csv_fixture(c("patient_id,EGFR,TP53", "001,1.2,4", "002,2.3,5"))
  clinical <- csv_fixture(c("patient_id,outcome", "002,0", "001,1"))
  d <- sb_read_cohort(expression, clinical)
  expect_identical(d$patient_id, c("001", "002"))
  expect_equal(d$outcome, c(1, 0))
  expect_equal(d$EGFR, c(1.2, 2.3))
  expect_named(sb_read_cohort(expression), c("patient_id", "EGFR", "TP53"))
  expect_s3_class(sb_binary_task(d, "patient_id", "outcome", c("EGFR", "TP53")),
    "stormbreaker_binary_task")
})

test_that("unmatched samples identify both missing and extra clinical IDs", {
  expression <- csv_fixture(c("patient_id,EGFR", "A,1", "B,2"))
  clinical <- csv_fixture(c("patient_id,outcome", "B,0", "C,1"))
  expect_error(sb_read_cohort(expression, clinical), "A.*C",
    class = "stormbreaker_unmatched_cohort_ids")
})

test_that("CSV input rejects duplicate IDs, columns and overlapping values", {
  expression <- csv_fixture(c("patient_id,EGFR", "A,1", "A,2"))
  expect_error(sb_read_cohort(expression), "A", class = "stormbreaker_duplicate_id")
  duplicate <- csv_fixture(c("patient_id,EGFR,EGFR", "A,1,2"))
  expect_error(sb_read_cohort(duplicate), "EGFR", class = "stormbreaker_duplicate_column_names")
  expression <- csv_fixture(c("patient_id,EGFR", "A,1", "B,2"))
  clinical <- csv_fixture(c("patient_id,EGFR", "A,0", "B,1"))
  expect_error(sb_read_cohort(expression, clinical), "EGFR",
    class = "stormbreaker_overlapping_cohort_columns")
  no_id <- csv_fixture(c("sample,EGFR", "A,1", "B,2"))
  expect_error(sb_read_cohort(no_id), "patient_id", class = "stormbreaker_missing_endpoint_columns")
})

test_that("missing measurements and character outcomes are not silently repaired", {
  expression <- csv_fixture(c("patient_id,EGFR", "A,NA", "B,2"))
  clinical <- csv_fixture(c("patient_id,outcome", "A,control", "B,cancer"))
  d <- sb_read_cohort(expression, clinical)
  expect_identical(is.na(d$EGFR), c(TRUE, FALSE))
  expect_identical(d$outcome, c("control", "cancer"))
  expect_error(sb_binary_task(d, "patient_id", "outcome", "EGFR"),
    "0/1", class = "stormbreaker_bad_binary_outcome")
})
})

# ---- example-data ----
local({
test_that("the synthetic survival example is valid and reproducible", {
  first <- sb_example_survival(n = 50, seed = 42)
  second <- sb_example_survival(n = 50, seed = 42)

  expect_identical(first, second)
  expect_identical(
    names(first),
    c("patient_id", "time", "event", "gene_a", "gene_b", "gene_c")
  )
  expect_identical(nrow(first), 50L)
  expect_identical(sum(first$event == 1L), 25L)
  expect_identical(sum(first$event == 0L), 25L)
  expect_identical(anyNA(first), FALSE)
  expect_identical(all(first$time > 0), TRUE)

  task <- sb_surv_task(
    first, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  expect_s3_class(task, "stormbreaker_surv_task")
})

test_that("synthetic data arguments are checked", {
  n_error <- capture_sb_error(sb_example_survival(n = 9, seed = 1))
  expect_s3_class(n_error, "stormbreaker_bad_integer")

  seed_error <- capture_sb_error(sb_example_survival(n = 20, seed = 0))
  expect_s3_class(seed_error, "stormbreaker_bad_integer")
})

test_that("the synthetic binary example is valid and reproducible", {
  first <- sb_example_binary(n = 50, seed = 42)
  second <- sb_example_binary(n = 50, seed = 42)

  expect_identical(first, second)
  expect_identical(
    names(first),
    c("patient_id", "outcome", "gene_a", "gene_b", "gene_c")
  )
  expect_identical(nrow(first), 50L)
  expect_identical(sum(first$outcome == 1L), 25L)
  expect_identical(sum(first$outcome == 0L), 25L)
  expect_identical(anyNA(first), FALSE)

  task <- sb_binary_task(
    first, "patient_id", "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  expect_s3_class(task, "stormbreaker_binary_task")
})

test_that("synthetic binary data arguments are checked", {
  n_error <- capture_sb_error(sb_example_binary(n = 9, seed = 1))
  expect_s3_class(n_error, "stormbreaker_bad_integer")

  seed_error <- capture_sb_error(sb_example_binary(n = 20, seed = 0))
  expect_s3_class(seed_error, "stormbreaker_bad_integer")
})
})
