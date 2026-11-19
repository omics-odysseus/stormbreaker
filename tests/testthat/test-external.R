# ---- external-evaluation ----
local({
# ---- external-evaluation ----

# Original-direction metric tests explicitly opt out of direction correction.
evaluate_external_original <- function(prediction, endpoint) {
  sb_evaluate_external(prediction, endpoint, reverse = FALSE)
}

make_sb_binary_final_prediction <- function(
    patient_ids,
    probability,
    pipeline_id = "none__logistic",
    learner = "logistic") {
  n <- length(patient_ids)
  result <- list(
    status = "success",
    task_type = "binary",
    pipeline_id = pipeline_id,
    learner = learner,
    id_col = "sample_id",
    predictions = data.frame(
      patient_id = patient_ids,
      pipeline_id = rep(pipeline_id, n),
      learner = rep(learner, n),
      positive_class = rep.int(1L, n),
      probability = probability,
      status = rep("success", n),
      error_class = rep(NA_character_, n),
      error_message = rep(NA_character_, n),
      stringsAsFactors = FALSE
    ),
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )
  class(result) <- "stormbreaker_final_prediction"
  result
}

make_sb_surv_final_prediction <- function(
    patient_ids,
    risk,
    survival_probability = NULL,
    prediction_horizon = NULL,
    time_unit = NULL,
    pipeline_id = "none__coxph",
    learner = "coxph",
    decision_thresholds = NULL) {
  n <- length(patient_ids)
  rows <- data.frame(
    patient_id = patient_ids,
    pipeline_id = rep(pipeline_id, n),
    learner = rep(learner, n),
    risk = risk,
    stringsAsFactors = FALSE
  )
  if (!is.null(prediction_horizon)) {
    rows$prediction_horizon <- rep(prediction_horizon, n)
    rows$time_unit <- rep(time_unit, n)
    rows$survival_probability <- survival_probability
  }
  rows$status <- rep("success", n)
  rows$error_class <- rep(NA_character_, n)
  rows$error_message <- rep(NA_character_, n)
  result <- list(
    status = "success",
    task_type = "survival",
    pipeline_id = pipeline_id,
    learner = learner,
    id_col = "sample_id",
    predictions = rows,
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )
  if (!is.null(decision_thresholds)) {
    result$decision_thresholds <- decision_thresholds
  }
  class(result) <- "stormbreaker_final_prediction"
  result
}

expect_external_failure <- function(result, stage, class_name) {
  expect_s3_class(result, "stormbreaker_external_evaluation")
  expect_identical(result$status, "failed")
  expect_identical(result$error_stage, stage)
  expect_identical(result$error_class, class_name)
  expect_identical(nrow(result$scored_predictions), 0L)
  expect_identical(result$metrics, list())
}

test_that("external evaluation keeps prediction-first arguments and explicit direction control", {
  expect_identical(
    names(formals(sb_evaluate_external)),
    c("prediction", "endpoint", "reverse")
  )
})

test_that("binary outcomes align by ID into endpoint order", {
  prediction <- make_sb_binary_final_prediction(
    c("p3", "p1", "p2"),
    c(0.9, 0.1, 0.8)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(
      subject = c("p1", "p2", "p3"),
      truth = c(0L, 1L, 1L)
    ),
    "subject",
    "truth"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_s3_class(result, "stormbreaker_external_evaluation")
  expect_identical(result$status, "success")
  expect_identical(result$scored_predictions$patient_id, c("p1", "p2", "p3"))
  expect_identical(result$scored_predictions$outcome, c(0L, 1L, 1L))
  expect_identical(result$scored_predictions$probability, c(0.1, 0.8, 0.9))
  expect_identical(result$patient_match$status, "success")
  expect_identical(result$patient_match$exact_set, TRUE)
  expect_identical(result$patient_match$output_order, "external_endpoint")
  expect_identical(
    names(result$scored_predictions),
    c("patient_id", "outcome", "probability", "status")
  )
})

test_that("binary external mean calibration is the observed-to-expected ratio", {
  prediction <- make_sb_binary_final_prediction(
    c("p3", "p1", "p2"),
    c(0.7, 0.1, 0.4)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(
      subject = c("p1", "p2", "p3"),
      truth = c(0L, 1L, 1L)
    ),
    "subject",
    "truth"
  )

  result <- evaluate_external_original(prediction, endpoint)
  calibration <- result$metrics$mean_calibration

  expect_identical(calibration$status, "success")
  expect_identical(calibration$ratio_direction, "observed_over_expected")
  expect_identical(calibration$observed_events, 2L)
  expect_equal(calibration$expected_events, 1.2, tolerance = 1e-15)
  expect_equal(calibration$observed_event_probability, 2 / 3)
  expect_equal(calibration$mean_predicted_event_probability, 0.4)
  expect_equal(calibration$mean_calibration_ratio, 2 / 1.2)
  expect_identical(calibration$target, 1)
  expect_match(
    paste(capture.output(print(result)), collapse = "\n"),
    "External mean calibration O:E:",
    fixed = TRUE
  )
})

test_that("binary mean calibration does not invent an O:E denominator", {
  prediction <- make_sb_binary_final_prediction(
    c("p1", "p2"),
    c(0, 0)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(id = c("p1", "p2"), outcome = c(0L, 1L)),
    "id",
    "outcome"
  )

  calibration <- evaluate_external_original(prediction, endpoint)$metrics$
    mean_calibration

  expect_identical(calibration$status, "not_evaluable")
  expect_true(is.na(calibration$mean_calibration_ratio))
  expect_identical(
    calibration$error_class,
    "stormbreaker_binary_mean_calibration_zero_expected"
  )
})

test_that("binary mean calibration does not label an overflowed O:E successful", {
  calibration <- sb_external_binary_mean_calibration_metric(data.frame(
    outcome = 1L,
    probability = 1e-309
  ))

  expect_identical(calibration$status, "not_evaluable")
  expect_identical(is.na(calibration$estimate), TRUE)
  expect_identical(is.na(calibration$mean_calibration_ratio), TRUE)
  expect_identical(calibration$expected_events, 1e-309)
  expect_identical(
    calibration$error_class,
    "stormbreaker_binary_mean_calibration_nonfinite_ratio"
  )
})

test_that("survival outcomes align by ID rather than row position", {
  prediction <- make_sb_surv_final_prediction(
    c("p2", "p3", "p1"),
    c(2, 3, 1)
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      subject = c("p1", "p2", "p3"),
      follow_up = c(1, 2, 3),
      event = c(1L, 1L, 0L)
    ),
    "subject",
    "follow_up",
    "event",
    "month"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(result$scored_predictions$patient_id, c("p1", "p2", "p3"))
  expect_identical(result$scored_predictions$risk, c(1, 2, 3))
  expect_identical(result$scored_predictions$time, c(1, 2, 3))
  expect_identical(result$scored_predictions$event, c(1L, 1L, 0L))
})

test_that("external evaluation rejects any nonidentical patient set", {
  prediction <- make_sb_binary_final_prediction(
    c("p1", "p2", "p3"),
    c(0.1, 0.5, 0.9)
  )
  endpoint_inputs <- list(
    data.frame(id = c("p1", "p2"), y = c(0L, 1L)),
    data.frame(id = c("p1", "p2", "p3", "p4"), y = c(0L, 1L, 0L, 1L)),
    data.frame(id = c("p1", "p2", "p4"), y = c(0L, 1L, 0L))
  )

  for (data in endpoint_inputs) {
    endpoint <- sb_binary_endpoint(data, "id", "y")
    result <- evaluate_external_original(prediction, endpoint)
    expect_external_failure(
      result,
      "patient_match",
      "stormbreaker_external_patient_mismatch"
    )
    expect_identical(result$patient_match$status, "failed")
    expect_identical(result$patient_match$exact_set, FALSE)
  }
})

test_that("invalid prediction identity and provenance fail before scoring", {
  endpoint <- sb_binary_endpoint(
    data.frame(id = c("p1", "p2"), y = c(0L, 1L)),
    "id",
    "y"
  )
  base <- make_sb_binary_final_prediction(c("p1", "p2"), c(0.2, 0.8))
  duplicate_id <- base
  duplicate_id$predictions$patient_id[2L] <- "p1"
  blank_id <- base
  blank_id$predictions$patient_id[1L] <- ""
  wrong_class <- base
  wrong_class$predictions$positive_class[1L] <- 0L
  wrong_pipeline <- base
  wrong_pipeline$predictions$pipeline_id[1L] <- "other"
  wrong_learner <- base
  wrong_learner$predictions$learner[1L] <- "other"
  failed_row <- base
  failed_row$predictions$status[1L] <- "failed"
  invalid <- list(
    duplicate_id,
    blank_id,
    wrong_class,
    wrong_pipeline,
    wrong_learner,
    failed_row
  )

  for (prediction in invalid) {
    class(prediction) <- "stormbreaker_final_prediction"
    result <- evaluate_external_original(prediction, endpoint)
    expect_external_failure(
      result,
      "prediction_contract",
      "stormbreaker_bad_external_prediction"
    )
  }

  failed <- base
  failed$status <- "failed"
  result <- evaluate_external_original(failed, endpoint)
  expect_external_failure(
    result,
    "prediction_contract",
    "stormbreaker_failed_external_prediction"
  )
})

test_that("endpoint mutation and task-type mismatch remain explicit", {
  binary_prediction <- make_sb_binary_final_prediction(
    c("p1", "p2"),
    c(0.2, 0.8)
  )
  binary_endpoint <- sb_binary_endpoint(
    data.frame(id = c("p1", "p2"), y = c(0L, 1L)),
    "id",
    "y"
  )
  binary_endpoint$n_positive <- 2L
  mutated <- evaluate_external_original(binary_prediction, binary_endpoint)
  expect_external_failure(
    mutated,
    "endpoint_contract",
    "stormbreaker_mutated_endpoint"
  )

  survival_endpoint <- sb_surv_endpoint(
    data.frame(
      id = c("p1", "p2"),
      time = c(1, 2),
      event = c(1L, 0L)
    ),
    "id",
    "time",
    "event",
    "month"
  )
  mismatch <- evaluate_external_original(binary_prediction, survival_endpoint)
  expect_external_failure(
    mismatch,
    "endpoint_contract",
    "stormbreaker_external_task_mismatch"
  )

  modeling_data <- data.frame(
    id = c("p1", "p2"),
    y = c(0L, 1L),
    x = c(1, 2)
  )
  modeling_task <- sb_binary_task(
    modeling_data,
    "id",
    "y",
    "x"
  )
  task_result <- evaluate_external_original(binary_prediction, modeling_task)
  expect_external_failure(
    task_result,
    "endpoint_contract",
    "stormbreaker_bad_endpoint"
  )
})

test_that("an existing prediction is evaluated without reprediction or a fit", {
  development <- sb_example_binary(n = 30L, seed = 701L)
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  external <- sb_example_binary(n = 12L, seed = 702L)
  external$patient_id <- paste0("external_", seq_len(nrow(external)))
  prediction <- sb_predict_final(
    fit,
    external[c("patient_id", "gene_a", "gene_b")],
    "patient_id"
  )
  endpoint <- sb_binary_endpoint(
    external,
    "patient_id",
    "outcome"
  )
  rm(fit)
  testthat::local_mocked_bindings(
    sb_predict_final = function(...) stop("must not repredict")
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(result$status, "success")
  expect_identical(
    result$evaluation_boundary$model_access,
    "none"
  )
  expect_identical(
    result$evaluation_boundary$outcome_use,
    "metric_evaluation_only"
  )
  expect_identical(
    result$scored_predictions$probability,
    prediction$predictions$probability
  )
  expect_null(result$fit)
  expect_null(result$model)
  expect_null(result[["task"]])
})

test_that("saved survival predictions integrate with endpoint-only evaluation", {
  development <- sb_example_survival(n = 30L, seed = 703L)
  task <- sb_surv_task(
    development,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 6
  )
  fit <- sb_fit_final(task, sb_pipeline())
  external <- sb_example_survival(n = 12L, seed = 704L)
  external$patient_id <- paste0("external_", seq_len(nrow(external)))
  prediction <- sb_predict_final(
    fit,
    external[c("patient_id", "gene_a", "gene_b")],
    "patient_id"
  )
  endpoint <- sb_surv_endpoint(
    external,
    "patient_id",
    "time",
    "event",
    "month"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_false(identical(result$status, "failed"))
  expect_identical(result$patient_match$exact_set, TRUE)
  expect_identical(
    result$scored_predictions$risk,
    prediction$predictions$risk
  )
  expect_identical(
    result$scored_predictions$survival_probability,
    prediction$predictions$survival_probability
  )
  expect_identical(
    unique(result$scored_predictions$prediction_horizon),
    6
  )
  calibration <- result$metrics$calibration
  direct_calibration <- sb_survival_calibration_metric(
    patient_id = result$scored_predictions$patient_id,
    time = result$scored_predictions$time,
    event = result$scored_predictions$event,
    survival_probability = result$scored_predictions$survival_probability,
    prediction_status = result$scored_predictions$status,
    prediction_horizon = 6,
    time_unit = "month",
    scope = "patient_level_external_cohort_at_locked_horizon"
  )
  expect_identical(calibration$status, "success")
  expect_identical(
    calibration$scope,
    "patient_level_external_cohort_at_locked_horizon"
  )
  expect_identical(calibration$knots, direct_calibration$knots)
  expect_identical(calibration$curve, direct_calibration$curve)
  expect_identical(
    calibration$patient_calibration,
    direct_calibration$patient_calibration
  )
  expect_identical(
    calibration[c("ici", "e50", "e90")],
    direct_calibration[c("ici", "e50", "e90")]
  )
  weak_calibration <- result$metrics$weak_calibration
  direct_weak_calibration <- sb_survival_weak_calibration_metric(
    patient_id = result$scored_predictions$patient_id,
    time = result$scored_predictions$time,
    event = result$scored_predictions$event,
    survival_probability = result$scored_predictions$survival_probability,
    prediction_status = result$scored_predictions$status,
    prediction_horizon = 6,
    time_unit = "month",
    scope = "patient_level_external_cohort_at_locked_horizon"
  )
  expect_identical(weak_calibration$status, "success")
  expect_identical(
    weak_calibration$patient_predictions$patient_id,
    result$scored_predictions$patient_id
  )
  expect_identical(weak_calibration$per_fold, data.frame(
    fold = integer(),
    status = character(),
    n = integer(),
    n_scored = integer(),
    n_events = integer(),
    n_censored = integer(),
    n_boundary_probability = integer(),
    n_followed_to_horizon = integer(),
    error_class = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  ))
  expect_equal(
    weak_calibration$mean_calibration_ratio,
    direct_weak_calibration$mean_calibration_ratio,
    tolerance = 0
  )
  expect_equal(
    weak_calibration$calibration_slope,
    direct_weak_calibration$calibration_slope,
    tolerance = 0
  )

  shuffled_external <- external[rev(seq_len(nrow(external))), , drop = FALSE]
  shuffled <- evaluate_external_original(
    prediction,
    sb_surv_endpoint(
      shuffled_external,
      "patient_id",
      "time",
      "event",
      "month"
    )
  )
  expect_identical(
    shuffled$metrics$calibration$patient_calibration$patient_id,
    shuffled_external$patient_id
  )
  expect_equal(
    shuffled$metrics$calibration[c("ici", "e50", "e90")],
    calibration[c("ici", "e50", "e90")],
    tolerance = 1e-14
  )
  expect_identical(shuffled$metrics$calibration$knots, calibration$knots)
  expect_equal(
    shuffled$metrics$calibration$curve,
    calibration$curve,
    tolerance = 1e-14
  )
  expect_identical(
    shuffled$metrics$weak_calibration$patient_predictions$patient_id,
    shuffled_external$patient_id
  )
  expect_equal(
    shuffled$metrics$weak_calibration[c(
      "mean_calibration_ratio", "calibration_slope"
    )],
    weak_calibration[c("mean_calibration_ratio", "calibration_slope")],
    tolerance = 1e-13
  )
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(printed, "External survival calibration ICI:", fixed = TRUE)
  expect_match(printed, "External mean calibration O:E:", fixed = TRUE)
  expect_match(printed, "External survival calibration slope:", fixed = TRUE)
})

test_that("external survival weak calibration is optional and isolated", {
  ids <- paste0("weak_external_", seq_len(12L))
  prediction <- make_sb_surv_final_prediction(
    ids,
    risk = seq(1.1, -0.7, length.out = 12L),
    survival_probability = exp(-exp(seq(-1.2, 0.9, length.out = 12L))),
    prediction_horizon = 6,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
      event = rep(c(1L, 0L), 6L)
    ),
    "id",
    "time",
    "event",
    "month"
  )
  baseline <- evaluate_external_original(prediction, endpoint)
  testthat::local_mocked_bindings(
    sb_survival_weak_calibration_metric = function(...) {
      condition <- simpleError("external weak calibration exploded")
      class(condition) <- c(
        "test_external_weak_calibration_error",
        class(condition)
      )
      stop(condition)
    }
  )

  observed <- evaluate_external_original(prediction, endpoint)

  expect_identical(observed$status, baseline$status)
  expect_identical(observed$scored_predictions, baseline$scored_predictions)
  expect_identical(
    observed$metrics[c(
      "harrell_c", "ipcw_brier", "time_auc", "uno_c", "calibration"
    )],
    baseline$metrics[c(
      "harrell_c", "ipcw_brier", "time_auc", "uno_c", "calibration"
    )]
  )
  expect_identical(
    observed$metrics$weak_calibration$status,
    "metric_failed"
  )
  expect_identical(
    observed$metrics$weak_calibration$error_class,
    "test_external_weak_calibration_error"
  )
})

test_that("external survival calibration is optional and isolated", {
  development <- sb_example_survival(n = 30L, seed = 703L)
  task <- sb_surv_task(
    development,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 6
  )
  fit <- sb_fit_final(task, sb_pipeline())
  external <- sb_example_survival(n = 12L, seed = 704L)
  external$patient_id <- paste0("external_", seq_len(nrow(external)))
  prediction <- sb_predict_final(
    fit,
    external[c("patient_id", "gene_a", "gene_b")],
    "patient_id"
  )
  endpoint <- sb_surv_endpoint(
    external,
    "patient_id",
    "time",
    "event",
    "month"
  )
  baseline <- evaluate_external_original(prediction, endpoint)
  testthat::local_mocked_bindings(
    sb_survival_calibration_metric = function(...) {
      condition <- simpleError("external calibration exploded")
      class(condition) <- c("test_external_calibration_error", class(condition))
      stop(condition)
    }
  )

  observed <- evaluate_external_original(prediction, endpoint)

  expect_identical(observed$status, baseline$status)
  expect_identical(observed$scored_predictions, baseline$scored_predictions)
  expect_identical(
    observed$metrics[c("harrell_c", "ipcw_brier", "time_auc", "uno_c")],
    baseline$metrics[c("harrell_c", "ipcw_brier", "time_auc", "uno_c")]
  )
  expect_identical(observed$metrics$calibration$status, "metric_failed")
  expect_identical(
    observed$metrics$calibration$error_class,
    "test_external_calibration_error"
  )
  expect_identical(
    observed$metrics$calibration$error_message,
    "external calibration exploded"
  )
  expect_match(
    paste(capture.output(print(observed)), collapse = "\n"),
    "External survival calibration: unavailable",
    fixed = TRUE
  )
})

test_that("boundary survival probability only disables external calibration", {
  ids <- paste0("p", seq_len(12L))
  prediction <- make_sb_surv_final_prediction(
    ids,
    risk = seq(1.2, 0.1, length.out = 12L),
    survival_probability = c(0, seq(0.15, 0.85, length.out = 11L)),
    prediction_horizon = 3,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 2, 2.5, 4, 5, 6, 1.5, 3.5, 7, 8, 9, 10),
      event = c(1L, 1L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L, 0L)
    ),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 3
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(result$status, "success")
  expect_identical(result$metrics$harrell_c$status, "success")
  expect_identical(result$metrics$ipcw_brier$status, "success")
  expect_identical(result$metrics$time_auc$status, "success")
  expect_identical(result$metrics$uno_c$status, "success")
  expect_identical(result$metrics$calibration$status, "incomplete")
  expect_identical(result$metrics$weak_calibration$status, "incomplete")
  expect_identical(
    result$metrics$weak_calibration$mean_calibration$status,
    "success"
  )
  expect_identical(
    result$metrics$weak_calibration$slope_fit$status,
    "not_evaluable"
  )
  expect_identical(
    result$metrics$calibration$error_class,
    "stormbreaker_survival_calibration_boundary_probability"
  )
  expect_identical(
    result$scored_predictions$survival_probability,
    prediction$predictions$survival_probability
  )
})

test_that("external binary metrics match hand calculations with ties", {
  prediction <- make_sb_binary_final_prediction(
    paste0("p", 1:4),
    c(0.9, 0.9, 0.2, 0.1)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(
      id = paste0("p", 1:4),
      y = c(1L, 0L, 1L, 0L)
    ),
    "id",
    "y"
  )
  result <- evaluate_external_original(prediction, endpoint)
  expected_brier <- mean((c(1, 0, 1, 0) - c(0.9, 0.9, 0.2, 0.1))^2)
  expected_log_loss <- mean(c(
    -log(0.9),
    -log1p(-0.9),
    -log(0.2),
    -log1p(-0.1)
  ))

  expect_equal(result$metrics$roc_auc$estimate, 0.625)
  expect_identical(result$metrics$roc_auc$comparable_pairs, 4L)
  expect_identical(result$metrics$roc_auc$correctly_ordered_pairs, 2)
  expect_identical(result$metrics$roc_auc$tied_pairs, 1)
  expect_equal(result$metrics$average_precision$estimate, 7 / 12)
  expect_equal(result$metrics$brier_score$estimate, expected_brier)
  expect_equal(result$metrics$log_loss$estimate, expected_log_loss)
  expect_identical(result$metrics$roc_auc$scope, "pooled_external_cohort")
  expect_identical(
    result$metrics$brier_score$scope,
    "patient_level_external_cohort"
  )
  expect_identical("fold" %in% names(result$metrics$roc_auc), FALSE)
  expect_identical("per_fold" %in% names(result$metrics$roc_auc), FALSE)
  expect_identical("reference" %in% names(result$metrics$brier_score), FALSE)
})

test_that("external calibration is descriptive and never changes probability", {
  probability <- c(0.2, 0.3, 0.7, 0.8)
  prediction <- make_sb_binary_final_prediction(
    paste0("p", 1:4),
    probability
  )
  endpoint <- sb_binary_endpoint(
    data.frame(
      id = paste0("p", 1:4),
      y = c(0L, 1L, 0L, 1L)
    ),
    "id",
    "y"
  )

  result <- evaluate_external_original(prediction, endpoint)
  calibration <- result$metrics$calibration

  expect_identical(calibration$status, "success")
  expect_true(is.finite(calibration$calibration_intercept))
  expect_true(is.finite(calibration$calibration_slope))
  expect_identical(
    result$scored_predictions$probability,
    probability
  )
  expect_null(result$recalibrated_probability)
})

test_that("constant external probabilities retain a typed calibration failure", {
  ids <- paste0("p", seq_len(8L))
  probability <- rep(0.5, length(ids))
  prediction <- make_sb_binary_final_prediction(ids, probability)
  endpoint <- sb_binary_endpoint(
    data.frame(id = ids, y = rep(c(0L, 1L), 4L)),
    "id",
    "y"
  )

  result <- evaluate_external_original(prediction, endpoint)
  calibration <- result$metrics$calibration

  expect_identical(calibration$status, "incomplete")
  expect_identical(calibration$intercept_fit$status, "success")
  expect_identical(calibration$slope_fit$status, "failed")
  expect_identical(
    calibration$slope_fit$error_class,
    "stormbreaker_binary_calibration_nonfinite_coefficients"
  )
  expect_identical(
    calibration$error_class,
    "stormbreaker_binary_calibration_nonfinite_coefficients"
  )
  expect_identical(result$scored_predictions$probability, probability)
})

test_that("external binary moderate calibration uses one aligned locked fit", {
  ids <- paste0("p", seq_len(24L))
  probability <- seq(0.03, 0.97, length.out = length(ids))
  outcome <- rep(c(0L, 1L), length.out = length(ids))
  prediction <- make_sb_binary_final_prediction(
    rev(ids),
    rev(probability)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(id = ids, y = outcome),
    "id",
    "y"
  )
  original_metric <- sb_binary_moderate_calibration_metric
  calls <- 0L
  observed_scored <- NULL
  testthat::local_mocked_bindings(
    sb_binary_moderate_calibration_metric = function(scored, scope) {
      calls <<- calls + 1L
      observed_scored <<- scored
      original_metric(scored, scope)
    }
  )

  result <- evaluate_external_original(prediction, endpoint)
  direct <- original_metric(
    result$scored_predictions,
    "patient_level_external_cohort"
  )

  expect_identical(calls, 1L)
  expect_identical(observed_scored, result$scored_predictions)
  expect_identical(result$metrics$moderate_calibration$status, "success")
  expect_identical(
    result$metrics$moderate_calibration$patient_calibration$patient_id,
    ids
  )
  expect_identical(result$metrics$moderate_calibration$knots, direct$knots)
  expect_identical(result$metrics$moderate_calibration$curve, direct$curve)
  expect_equal(
    result$metrics$moderate_calibration[c("ici", "e50", "e90")],
    direct[c("ici", "e50", "e90")],
    tolerance = 0
  )
  expect_identical(result$scored_predictions$probability, probability)
  expect_match(
    paste(capture.output(print(result)), collapse = "\n"),
    "External binary moderate calibration ICI:",
    fixed = TRUE
  )
})

test_that("external binary moderate calibration is endpoint-order invariant", {
  ids <- paste0("p", seq_len(24L))
  probability <- seq(0.04, 0.96, length.out = length(ids))
  outcome <- rep(c(0L, 1L), length.out = length(ids))
  prediction <- make_sb_binary_final_prediction(ids, probability)
  first <- evaluate_external_original(
    prediction,
    sb_binary_endpoint(data.frame(id = ids, y = outcome), "id", "y")
  )
  order <- c(seq(2L, length(ids), by = 2L), seq(1L, length(ids), by = 2L))
  second <- evaluate_external_original(
    prediction,
    sb_binary_endpoint(
      data.frame(id = ids[order], y = outcome[order]),
      "id",
      "y"
    )
  )

  expect_identical(first$metrics$moderate_calibration$status, "success")
  expect_identical(second$metrics$moderate_calibration$status, "success")
  expect_equal(
    second$metrics$moderate_calibration[c("ici", "e50", "e90")],
    first$metrics$moderate_calibration[c("ici", "e50", "e90")],
    tolerance = 1e-13
  )
  expect_equal(
    second$metrics$moderate_calibration$curve,
    first$metrics$moderate_calibration$curve,
    tolerance = 1e-13
  )
})

test_that("an external binary moderate-calibration error is isolated", {
  ids <- paste0("p", seq_len(24L))
  probability <- seq(0.03, 0.97, length.out = length(ids))
  endpoint <- sb_binary_endpoint(
    data.frame(id = ids, y = rep(c(0L, 1L), 12L)),
    "id",
    "y"
  )
  prediction <- make_sb_binary_final_prediction(ids, probability)
  control_result <- evaluate_external_original(prediction, endpoint)
  calls <- 0L
  testthat::local_mocked_bindings(
    sb_binary_moderate_calibration_metric = function(...) {
      calls <<- calls + 1L
      condition <- simpleError("external binary moderate calibration exploded")
      class(condition) <- c(
        "test_external_binary_moderate_calibration_error",
        class(condition)
      )
      stop(condition)
    }
  )

  observed <- evaluate_external_original(prediction, endpoint)

  expect_identical(calls, 1L)
  expect_identical(observed$status, control_result$status)
  expect_identical(
    observed$scored_predictions,
    control_result$scored_predictions
  )
  expect_identical(
    observed$metrics[names(observed$metrics) != "moderate_calibration"],
    control_result$metrics[
      names(control_result$metrics) != "moderate_calibration"
    ]
  )
  expect_identical(
    observed$metrics$moderate_calibration$status,
    "metric_failed"
  )
  expect_identical(
    observed$metrics$moderate_calibration$error_class,
    "test_external_binary_moderate_calibration_error"
  )
  expect_match(
    paste(capture.output(print(observed)), collapse = "\n"),
    "External binary moderate calibration: unavailable",
    fixed = TRUE
  )
})

test_that("wrong boundary probabilities retain infinite unclipped log loss", {
  prediction <- make_sb_binary_final_prediction(
    c("p1", "p2"),
    c(0, 1)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(id = c("p1", "p2"), y = c(1L, 0L)),
    "id",
    "y"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(result$metrics$log_loss$status, "success")
  expect_identical(result$metrics$log_loss$estimate, Inf)
  expect_identical(result$metrics$log_loss$n_infinite_loss, 2L)
  expect_identical(result$metrics$log_loss$probability_clipping, "none")
  expect_identical(result$metrics$calibration$status, "incomplete")
  expect_identical(
    result$metrics$calibration$input_status,
    "not_evaluable"
  )
  expect_identical(
    result$metrics$calibration$error_class,
    "stormbreaker_binary_calibration_boundary_probability"
  )
  expect_identical(result$metrics$moderate_calibration$status, "incomplete")
  expect_identical(
    result$metrics$moderate_calibration$n_boundary_probability,
    2L
  )
  expect_identical(
    result$metrics$moderate_calibration$probability_clipping,
    "none"
  )
})

test_that("single-patient and single-class binary evidence stays partial", {
  cases <- list(
    list(id = "p1", y = 1L, p = 0.7),
    list(id = c("p1", "p2"), y = c(1L, 1L), p = c(0.7, 0.8))
  )

  for (case in cases) {
    prediction <- make_sb_binary_final_prediction(case$id, case$p)
    endpoint <- sb_binary_endpoint(
      data.frame(id = case$id, y = case$y),
      "id",
      "y"
    )
    result <- evaluate_external_original(prediction, endpoint)

    expect_identical(result$status, "incomplete")
    expect_identical(result$metrics$roc_auc$status, "not_evaluable")
    expect_identical(
      result$metrics$average_precision$status,
      "not_evaluable"
    )
    expect_identical(result$metrics$brier_score$status, "success")
    expect_identical(result$metrics$log_loss$status, "success")
    expect_identical(result$metrics$calibration$status, "incomplete")
    expect_identical(
      result$metrics$moderate_calibration$status,
      "incomplete"
    )
    expect_identical(
      result$metrics$moderate_calibration$input_status,
      "not_evaluable"
    )
    expect_identical(is.nan(result$metrics$brier_score$estimate), FALSE)
  }
})

test_that("unrelated endpoint input columns cannot change binary scoring", {
  prediction <- make_sb_binary_final_prediction(
    c("p1", "p2", "p3", "p4"),
    c(0.1, 0.4, 0.6, 0.9)
  )
  first_data <- data.frame(
    id = c("p1", "p2", "p3", "p4"),
    y = c(0L, 0L, 1L, 1L),
    gene = c(1, 2, 3, 4)
  )
  second_data <- first_data[c("gene", "y", "id")]
  second_data$gene <- c(Inf, NA, -100, 100)
  first <- evaluate_external_original(
    prediction,
    sb_binary_endpoint(first_data, "id", "y")
  )
  second <- evaluate_external_original(
    prediction,
    sb_binary_endpoint(second_data, "id", "y")
  )

  expect_identical(first$scored_predictions, second$scored_predictions)
  expect_identical(first$metrics, second$metrics)

  changed <- first_data
  changed$y <- 1L - changed$y
  changed_result <- evaluate_external_original(
    prediction,
    sb_binary_endpoint(changed, "id", "y")
  )
  expect_identical(
    changed_result$scored_predictions$probability,
    first$scored_predictions$probability
  )
  expect_false(identical(changed_result$metrics, first$metrics))
})

test_that("external Harrell C uses the complete supplied cohort", {
  ids <- paste0("p", 1:5)
  risk <- c(1.2, -0.5, 0.7, 2.1, 0.1)
  time <- c(1, 4, 3, 2, 5)
  event <- c(1L, 0L, 1L, 1L, 0L)
  prediction <- make_sb_surv_final_prediction(ids, risk)
  endpoint <- sb_surv_endpoint(
    data.frame(id = ids, time = time, event = event),
    "id",
    "time",
    "event",
    "month"
  )
  expected <- survival::concordance(
    survival::Surv(time, event) ~ risk,
    reverse = TRUE,
    timewt = "n"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(result$status, "success")
  expect_equal(result$metrics$harrell_c$estimate, expected$concordance)
  expect_identical(result$metrics$harrell_c$scope, "external_cohort")
  expect_identical("ipcw_brier" %in% names(result$metrics), FALSE)
  expect_identical("time_auc" %in% names(result$metrics), FALSE)
  expect_identical("uno_c" %in% names(result$metrics), FALSE)
  expect_false("weak_calibration" %in% names(result$metrics))
})

test_that("external Harrell C preserves its count-ratio rounding", {
  scored <- data.frame(
    time = c(1, 7, 5),
    event = c(1L, 1L, 1L),
    risk = c(-0.8, 0.6, -0.8)
  )
  official <- survival::concordance(
    survival::Surv(time, event) ~ risk,
    data = scored,
    reverse = TRUE,
    timewt = "n"
  )
  expected <- unname(
    (official$count[["concordant"]] + 0.5 * official$count[["tied.x"]]) /
      sum(official$count[c("concordant", "discordant", "tied.x")])
  )

  observed <- sb_external_harrell_c_metric(scored)

  expect_identical(observed$estimate, expected)
  expect_false(identical(observed$estimate, official$concordance))
})

test_that("external Harrell C maps all tie counts without swapping", {
  count_fields <- c(
    "concordant", "discordant", "tied_risk", "tied_time",
    "tied_both", "comparable_pairs"
  )
  cases <- list(
    list(risk = c(0.9, 0.8, 0.2), expected = c(2, 0, 0, 1, 0, 2)),
    list(risk = c(0.5, 0.5, 0.2), expected = c(2, 0, 0, 0, 1, 2))
  )

  for (case in cases) {
    observed <- sb_external_harrell_c_metric(data.frame(
      time = c(1, 1, 2),
      event = c(1L, 1L, 0L),
      risk = case$risk
    ))

    expect_equal(
      unname(unlist(observed[count_fields], use.names = FALSE)),
      case$expected,
      tolerance = 0
    )
  }
})

test_that("external IPCW Brier follows the declared patient-level formula", {
  ids <- c("a", "b", "c", "d", "e")
  probability <- c(0.2, 0.4, 0.6, 0.8, 0.3)
  prediction <- make_sb_surv_final_prediction(
    ids,
    risk = c(2, 0, 1, -1, 1.5),
    survival_probability = probability,
    prediction_horizon = 3,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 2, 3, 4, 2),
      event = c(1L, 0L, 0L, 1L, 1L)
    ),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 3
  )

  testthat::local_mocked_bindings(
    sb_fit_event_km = function(...) stop("event KM must not run")
  )
  result <- evaluate_external_original(prediction, endpoint)
  brier <- result$metrics$ipcw_brier
  time_auc <- result$metrics$time_auc
  uno_c <- result$metrics$uno_c

  expected_weighted_loss <-
    0.2^2 +
    0.3^2 +
    (1 - 0.6)^2 * 1.5 +
    (1 - 0.8)^2 * 1.5
  expect_identical(brier$status, "success")
  expect_equal(brier$estimate, expected_weighted_loss / 5)
  expect_equal(brier$weighted_loss_sum, expected_weighted_loss)
  expect_equal(brier$weight_sum, 5)
  expect_identical(brier$n_scored, 4L)
  expect_identical(brier$n_censored_at_horizon, 1L)
  expect_identical(brier$censoring_model, "external_reverse_kaplan_meier")
  expect_identical(
    brier$censoring_assumption,
    "marginal_noninformative_censoring"
  )
  expect_identical("reference" %in% names(brier), FALSE)
  expect_identical("skill_score" %in% names(brier), FALSE)
  expect_identical("per_fold" %in% names(brier), FALSE)
  expect_identical(time_auc$status, "success")
  expect_equal(time_auc$estimate, 1, tolerance = 0)
  expect_identical(time_auc$n_cases, 2L)
  expect_identical(time_auc$n_controls, 1L)
  expect_identical(time_auc$n_excluded_censored, 2L)
  expect_equal(time_auc$case_weight_sum, 2, tolerance = 0)
  expect_equal(time_auc$weighted_pair_credit, 2, tolerance = 0)
  expect_equal(time_auc$pair_denominator, 2, tolerance = 0)
  expect_identical(time_auc$marker, "one_minus_survival_probability")
  expect_identical(uno_c$status, "success")
  expect_equal(uno_c$estimate, 1, tolerance = 0)
  expect_identical(uno_c$n_eligible_events, 2L)
  expect_identical(uno_c$n_comparable_pairs, 7L)
  expect_identical(uno_c$n_tied_time_comparable_pairs, 1L)
  expect_equal(uno_c$weighted_pair_credit, 7, tolerance = 0)
  expect_equal(uno_c$weighted_pair_denominator, 7, tolerance = 0)
  expect_identical(uno_c$marker, "risk")
})

test_that("external time AUC reuses one curve without Brier status coupling", {
  ids <- c("case", "control", "censored")
  prediction <- make_sb_surv_final_prediction(
    ids,
    risk = c(2, 0, 1),
    survival_probability = c(0.2, 0.8, 0.5),
    prediction_horizon = 2,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 3, 1.5),
      event = c(1L, 1L, 0L)
    ),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 2
  )
  original_fit <- sb_fit_censoring_km
  original_brier <- sb_external_ipcw_brier_metric
  censoring_calls <- 0L
  testthat::local_mocked_bindings(
    sb_fit_censoring_km = function(...) {
      censoring_calls <<- censoring_calls + 1L
      original_fit(...)
    },
    sb_external_ipcw_brier_metric = function(...) {
      metric <- original_brier(...)
      metric$status <- "incomplete"
      metric$estimate <- NA_real_
      metric
    }
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(censoring_calls, 1L)
  expect_identical(result$status, "incomplete")
  expect_identical(result$metrics$ipcw_brier$status, "incomplete")
  expect_identical(result$metrics$time_auc$status, "success")
  expect_equal(result$metrics$time_auc$estimate, 1, tolerance = 0)
  expect_identical(result$metrics$uno_c$status, "success")
  expect_equal(result$metrics$uno_c$estimate, 1, tolerance = 0)
})

test_that("external Uno C uses the aligned external censoring distribution", {
  ids <- paste0("external_g_", 1:5)
  time <- c(1, 2, 3, 4, 5)
  event <- c(1L, 0L, 1L, 0L, 0L)
  risk <- c(0, 1, 3, 2, -1)
  prediction <- make_sb_surv_final_prediction(
    ids,
    risk,
    survival_probability = c(0.8, 0.7, 0.4, 0.3, 0.2),
    prediction_horizon = 5,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(id = ids, time = time, event = event),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 5
  )

  result <- evaluate_external_original(prediction, endpoint)
  scored <- result$scored_predictions
  external_curve <- sb_fit_censoring_km(time, event)
  development_curve <- sb_fit_censoring_km(1:5, rep.int(1L, 5L))
  expected <- sb_external_uno_c_metric(
    scored,
    list(time_unit = "month"),
    external_curve
  )
  counterfactual <- sb_external_uno_c_metric(
    scored,
    list(time_unit = "month"),
    development_curve
  )

  expect_identical(result$metrics$uno_c$status, "success")
  expect_equal(result$metrics$uno_c$estimate, expected$estimate, tolerance = 0)
  expect_false(isTRUE(all.equal(expected$estimate, counterfactual$estimate)))
  expect_equal(expected$estimate, 4.55555555555556 / 7.55555555555556)
  expect_equal(counterfactual$estimate, 0.5)
})

test_that("an Uno exception remains inside the external metric", {
  ids <- paste0("uno_error_", 1:6)
  prediction <- make_sb_surv_final_prediction(
    ids,
    risk = c(3, 2, 1, 0, -1, -2),
    survival_probability = c(0.2, 0.3, 0.4, 0.6, 0.7, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = 1:6,
      event = c(1L, 1L, 0L, 1L, 0L, 1L)
    ),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 3
  )
  baseline <- evaluate_external_original(prediction, endpoint)
  testthat::local_mocked_bindings(
    sb_uno_c_stratum = function(...) {
      condition <- simpleError("injected Uno failure")
      class(condition) <- c("test_uno_error", class(condition))
      stop(condition)
    }
  )

  observed <- evaluate_external_original(prediction, endpoint)

  expect_identical(observed$status, "incomplete")
  expect_identical(observed$scored_predictions, baseline$scored_predictions)
  expect_identical(
    observed$metrics[names(observed$metrics) != "uno_c"],
    baseline$metrics[names(baseline$metrics) != "uno_c"]
  )
  expect_identical(observed$metrics$uno_c$status, "metric_failed")
  expect_identical(observed$metrics$uno_c$error_class, "test_uno_error")
  expect_identical(
    observed$metrics$uno_c$error_message,
    "injected Uno failure"
  )
})

test_that("all-event IPCW remains valid beyond the largest observed time", {
  ids <- paste0("p", 1:4)
  probability <- c(0.1, 0.2, 0.3, 0.4)
  prediction <- make_sb_surv_final_prediction(
    ids,
    risk = c(4, 3, 2, 1),
    survival_probability = probability,
    prediction_horizon = 5,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(id = ids, time = 1:4, event = 1L),
    "id",
    "time",
    "event",
    "month"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(result$status, "incomplete")
  expect_identical(result$metrics$ipcw_brier$status, "success")
  expect_identical(result$metrics$time_auc$status, "not_evaluable")
  expect_identical(
    result$metrics$time_auc$error_class,
    "stormbreaker_no_time_auc_control"
  )
  expect_equal(
    result$metrics$ipcw_brier$estimate,
    mean(probability^2)
  )
  expect_true(all(
    result$metrics$ipcw_brier$censoring_curve$survival_after == 1
  ))
})

test_that("survival endpoint information propagates without fabricated metrics", {
  one_prediction <- make_sb_surv_final_prediction("p1", 1)
  one_endpoint <- sb_surv_endpoint(
    data.frame(id = "p1", time = 2, event = 1L),
    "id",
    "time",
    "event",
    "month"
  )
  one <- evaluate_external_original(one_prediction, one_endpoint)
  expect_identical(one$status, "incomplete")
  expect_identical(one$metrics$harrell_c$status, "not_evaluable")
  expect_identical(is.nan(one$metrics$harrell_c$estimate), FALSE)

  after_ids <- c("a", "b", "c")
  after_prediction <- make_sb_surv_final_prediction(
    after_ids,
    risk = c(1, 2, 3),
    survival_probability = c(0.8, 0.7, 0.6),
    prediction_horizon = 3,
    time_unit = "month"
  )
  after_endpoint <- sb_surv_endpoint(
    data.frame(id = after_ids, time = c(4, 5, 6), event = 0L),
    "id",
    "time",
    "event",
    "month"
  )
  after <- evaluate_external_original(after_prediction, after_endpoint)
  expect_identical(after$status, "incomplete")
  expect_identical(after$metrics$harrell_c$status, "not_evaluable")
  expect_identical(after$metrics$ipcw_brier$status, "success")
  expect_equal(
    after$metrics$ipcw_brier$estimate,
    mean((1 - c(0.8, 0.7, 0.6))^2)
  )

  before_prediction <- make_sb_surv_final_prediction(
    after_ids,
    risk = c(1, 2, 3),
    survival_probability = c(0.8, 0.7, 0.6),
    prediction_horizon = 4,
    time_unit = "month"
  )
  before_endpoint <- sb_surv_endpoint(
    data.frame(id = after_ids, time = c(1, 2, 3), event = 0L),
    "id",
    "time",
    "event",
    "month"
  )
  before <- evaluate_external_original(before_prediction, before_endpoint)
  expect_identical(before$status, "incomplete")
  expect_identical(before$metrics$ipcw_brier$status, "incomplete")
  expect_identical(before$metrics$ipcw_brier$n_scored, 0L)
  expect_identical(
    before$metrics$ipcw_brier$error_class,
    "stormbreaker_no_evaluable_horizon_status"
  )
  expect_identical(is.nan(before$metrics$ipcw_brier$estimate), FALSE)
})

test_that("survival horizon and time-unit states are never silently changed", {
  ids <- c("p1", "p2", "p3")
  risk_only <- make_sb_surv_final_prediction(ids, c(3, 2, 1))
  horizon <- make_sb_surv_final_prediction(
    ids,
    c(3, 2, 1),
    survival_probability = c(0.2, 0.5, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )
  data <- data.frame(
    id = ids,
    time = c(1, 3, 4),
    event = c(1L, 0L, 1L)
  )
  no_horizon <- sb_surv_endpoint(
    data, "id", "time", "event", "month"
  )
  declared <- sb_surv_endpoint(
    data, "id", "time", "event", "month", 3
  )

  risk_result <- evaluate_external_original(risk_only, no_horizon)
  expect_identical(risk_result$status, "success")
  expect_identical("ipcw_brier" %in% names(risk_result$metrics), FALSE)
  expect_identical("time_auc" %in% names(risk_result$metrics), FALSE)
  expect_identical("uno_c" %in% names(risk_result$metrics), FALSE)
  expect_false("weak_calibration" %in% names(risk_result$metrics))

  risk_mismatch <- evaluate_external_original(risk_only, declared)
  expect_external_failure(
    risk_mismatch,
    "horizon_semantics",
    "stormbreaker_external_horizon_mismatch"
  )

  inherited <- evaluate_external_original(horizon, no_horizon)
  expect_identical(inherited$status, "success")
  expect_identical(
    unique(inherited$scored_predictions$prediction_horizon),
    3
  )

  matched <- evaluate_external_original(horizon, declared)
  expect_identical(matched$status, "success")

  different_horizon <- sb_surv_endpoint(
    data, "id", "time", "event", "month", 2
  )
  horizon_mismatch <- evaluate_external_original(horizon, different_horizon)
  expect_external_failure(
    horizon_mismatch,
    "horizon_semantics",
    "stormbreaker_external_horizon_mismatch"
  )

  different_unit <- sb_surv_endpoint(
    data, "id", "time", "event", "day"
  )
  unit_mismatch <- evaluate_external_original(horizon, different_unit)
  expect_external_failure(
    unit_mismatch,
    "horizon_semantics",
    "stormbreaker_external_time_unit_mismatch"
  )
})

test_that("partial or row-inconsistent survival probability semantics fail", {
  ids <- c("p1", "p2", "p3")
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 4, 5),
      event = c(1L, 0L, 1L)
    ),
    "id",
    "time",
    "event",
    "month"
  )
  base <- make_sb_surv_final_prediction(
    ids,
    c(3, 2, 1),
    survival_probability = c(0.2, 0.5, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )
  partial <- base
  partial$predictions$survival_probability <- NULL
  different_horizon <- base
  different_horizon$predictions$prediction_horizon[2L] <- 4
  different_unit <- base
  different_unit$predictions$time_unit[2L] <- "day"

  for (prediction in list(partial, different_horizon, different_unit)) {
    result <- evaluate_external_original(prediction, endpoint)
    expect_external_failure(
      result,
      "prediction_contract",
      "stormbreaker_bad_external_prediction"
    )
  }
})

test_that("a required zero censoring denominator stays incomplete", {
  ids <- c("p1", "p2", "p3")
  prediction <- make_sb_surv_final_prediction(
    ids,
    c(3, 2, 1),
    survival_probability = c(0.2, 0.5, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 4, 5),
      event = c(1L, 0L, 1L)
    ),
    "id",
    "time",
    "event",
    "month"
  )
  testthat::local_mocked_bindings(
    sb_censoring_survival_at = function(...) 0
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(result$status, "incomplete")
  expect_identical(result$metrics$ipcw_brier$status, "incomplete")
  expect_identical(
    result$metrics$ipcw_brier$error_class,
    "stormbreaker_zero_external_censoring_survival"
  )
})

test_that("unrelated survival endpoint input columns never enter scoring", {
  ids <- c("p1", "p2", "p3")
  prediction <- make_sb_surv_final_prediction(
    ids,
    c(3, 2, 1),
    survival_probability = c(0.2, 0.5, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )
  first_data <- data.frame(
    id = ids,
    time = c(1, 4, 5),
    event = c(1L, 0L, 1L),
    gene = c(1, 2, 3)
  )
  second_data <- first_data[c("gene", "event", "time", "id")]
  second_data$gene <- c(NA, Inf, -Inf)

  first <- evaluate_external_original(
    prediction,
    sb_surv_endpoint(
      first_data, "id", "time", "event", "month"
    )
  )
  second <- evaluate_external_original(
    prediction,
    sb_surv_endpoint(
      second_data, "id", "time", "event", "month"
    )
  )

  expect_identical(first$scored_predictions, second$scored_predictions)
  expect_identical(first$metrics, second$metrics)
})

test_that("horizon probabilities may vary by patient but remain bounded", {
  ids <- c("p1", "p2", "p3")
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 4, 5),
      event = c(1L, 0L, 1L)
    ),
    "id",
    "time",
    "event",
    "month"
  )
  valid <- make_sb_surv_final_prediction(
    ids,
    c(3, 2, 1),
    survival_probability = c(0.2, 0.5, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )
  result <- evaluate_external_original(valid, endpoint)
  expect_identical(result$status, "success")
  expect_identical(
    result$scored_predictions$survival_probability,
    c(0.2, 0.5, 0.8)
  )

  invalid <- valid
  invalid$predictions$survival_probability[2L] <- 2
  failure <- evaluate_external_original(invalid, endpoint)
  expect_external_failure(
    failure,
    "prediction_contract",
    "stormbreaker_bad_external_prediction"
  )
})

test_that("external evaluation survives prediction serialization", {
  binary_prediction <- make_sb_binary_final_prediction(
    c("b1", "b2"),
    c(0.2, 0.8)
  )
  binary_endpoint <- sb_binary_endpoint(
    data.frame(id = c("b1", "b2"), y = c(0L, 1L)),
    "id",
    "y"
  )
  survival_prediction <- make_sb_surv_final_prediction(
    c("s1", "s2", "s3"),
    c(3, 2, 1),
    survival_probability = c(0.2, 0.5, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )
  survival_endpoint <- sb_surv_endpoint(
    data.frame(
      id = c("s1", "s2", "s3"),
      time = c(1, 4, 5),
      event = c(1L, 0L, 1L)
    ),
    "id",
    "time",
    "event",
    "month"
  )
  risk_prediction <- make_sb_surv_final_prediction(
    c("s1", "s2", "s3"),
    c(3, 2, 1)
  )

  binary_before <- evaluate_external_original(
    binary_prediction,
    binary_endpoint
  )
  binary_after <- evaluate_external_original(
    unserialize(serialize(binary_prediction, NULL)),
    binary_endpoint
  )
  survival_before <- evaluate_external_original(
    survival_prediction,
    survival_endpoint
  )
  survival_after <- evaluate_external_original(
    unserialize(serialize(survival_prediction, NULL)),
    survival_endpoint
  )
  risk_before <- evaluate_external_original(risk_prediction, survival_endpoint)
  risk_after <- evaluate_external_original(
    unserialize(serialize(risk_prediction, NULL)),
    survival_endpoint
  )

  expect_identical(binary_after, binary_before)
  expect_identical(survival_after, survival_before)
  expect_identical(risk_after, risk_before)
})

test_that("binary external decision curve uses development-declared thresholds", {
  development <- make_sb_binary_fixture(n = 60L, seed = 20260910L)
  external <- make_sb_binary_fixture(n = 30L, seed = 20260911L)
  external$patient_id <- paste0("external_", seq_len(nrow(external)))
  predictors <- c("gene_a", "gene_b", "gene_c")
  thresholds <- c(0.20, 0.40, 0.60)
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    predictors,
    decision_thresholds = thresholds
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  prediction <- sb_predict_final(
    fit,
    external[c("patient_id", predictors)],
    "patient_id"
  )
  endpoint <- sb_binary_endpoint(
    external, "patient_id", "outcome"
  )

  result <- evaluate_external_original(prediction, endpoint)
  metric <- result$metrics$decision_curve

  expect_identical(result$status, "success")
  expect_identical(metric$status, "success")
  expect_identical(metric$scope, "patient_level_external_cohort")
  expect_identical(metric$curve$threshold_probability, thresholds)
  for (index in seq_along(thresholds)) {
    threshold <- thresholds[[index]]
    selected <- result$scored_predictions$probability >= threshold
    true_positive <- as.integer(sum(
      selected & result$scored_predictions$outcome == 1L
    ))
    false_positive <- as.integer(sum(
      selected & result$scored_predictions$outcome == 0L
    ))
    expected <- true_positive / nrow(external) -
      false_positive / nrow(external) * threshold / (1 - threshold)

    expect_identical(metric$curve$true_positive[[index]], true_positive)
    expect_identical(metric$curve$false_positive[[index]], false_positive)
    expect_equal(
      metric$curve$model_net_benefit[[index]],
      expected,
      tolerance = 1e-15
    )
  }

  reordered <- external[rev(seq_len(nrow(external))), , drop = FALSE]
  reordered_result <- evaluate_external_original(
    prediction,
    sb_binary_endpoint(reordered, "patient_id", "outcome")
  )
  expect_identical(
    reordered_result$metrics$decision_curve$curve,
    metric$curve
  )
})

test_that("binary external evaluation omits undeclared decision curves", {
  prediction <- make_sb_binary_final_prediction(
    c("p1", "p2", "p3", "p4"),
    c(0.1, 0.4, 0.6, 0.9)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(
      id = c("p1", "p2", "p3", "p4"),
      y = c(0L, 0L, 1L, 1L)
    ),
    "id",
    "y"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_false("decision_curve" %in% names(result$metrics))
})

test_that("invalid prediction-side decision thresholds fail before scoring", {
  prediction <- make_sb_binary_final_prediction(
    c("p1", "p2"),
    c(0.2, 0.8)
  )
  prediction$decision_thresholds <- c(0.2, 1)
  endpoint <- sb_binary_endpoint(
    data.frame(id = c("p1", "p2"), y = c(0L, 1L)),
    "id",
    "y"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_external_failure(
    result,
    "prediction_contract",
    "stormbreaker_bad_decision_thresholds"
  )
})

test_that("external evaluation printing reports metrics without proof claims", {
  prediction <- make_sb_binary_final_prediction(
    c("p1", "p2"),
    c(0.2, 0.8)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(id = c("p1", "p2"), y = c(0L, 1L)),
    "id",
    "y"
  )
  result <- evaluate_external_original(prediction, endpoint)

  output <- paste(capture.output(print(result)), collapse = "\n")

  expect_match(output, "precomputed-prediction external evaluation")
  expect_match(output, "External ROC AUC")
  expect_match(output, "Patients: 2")
  expect_no_match(output, "transportability.*proven", ignore.case = TRUE)
  expect_no_match(output, "clinical utility.*proven", ignore.case = TRUE)
})

test_that("external survival DCA uses fixed predictions after ID alignment", {
  endpoint_ids <- paste0("p", 1:6)
  prediction_order <- c(6L, 2L, 5L, 1L, 4L, 3L)
  time <- 1:6
  event <- c(1L, 0L, 1L, 0L, 1L, 0L)
  event_probability <- c(0.80, 0.70, 0.60, 0.40, 0.55, 0.10)
  prediction <- make_sb_surv_final_prediction(
    endpoint_ids[prediction_order],
    risk = seq_along(prediction_order),
    survival_probability = (1 - event_probability)[prediction_order],
    prediction_horizon = 4,
    time_unit = "month",
    decision_thresholds = 0.50
  )
  endpoint <- sb_surv_endpoint(
    data.frame(id = endpoint_ids, time = time, event = event),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 4
  )

  result <- evaluate_external_original(prediction, endpoint)
  metric <- result$metrics$decision_curve
  row <- metric$curve[1L, ]

  expect_identical(result$status, "success")
  expect_identical(metric$status, "success")
  expect_identical(metric$scope, "patient_level_external_cohort")
  expect_identical(result$scored_predictions$patient_id, endpoint_ids)
  expect_equal(
    1 - result$scored_predictions$survival_probability,
    event_probability
  )
  expect_identical(row$n_selected, 4L)
  expect_equal(row$selected_survival_probability, 3 / 8, tolerance = 0)
  expect_equal(row$model_net_benefit, 1 / 6, tolerance = 1e-15)
  expect_equal(row$treat_all_net_benefit, -1 / 4, tolerance = 0)
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(
    printed,
    "External survival decision curve: 1 declared threshold"
  )
})

test_that("external survival DCA retains endpoint-specific estimability", {
  ids <- c("p1", "p2")
  for (event in list(c(1L, 1L), c(0L, 0L))) {
    prediction <- make_sb_surv_final_prediction(
      ids,
      risk = c(1, 2),
      survival_probability = c(0.40, 0.40),
      prediction_horizon = 2,
      time_unit = "month",
      decision_thresholds = 0.50
    )
    endpoint <- sb_surv_endpoint(
      data.frame(id = ids, time = c(1, 3), event = event),
      "id",
      "time",
      "event",
      "month",
      prediction_horizon = 2
    )
    result <- evaluate_external_original(prediction, endpoint)

    expect_identical(result$metrics$decision_curve$status, "success")
    expect_true(is.finite(
      result$metrics$decision_curve$curve$model_net_benefit
    ))
  }

  prediction <- make_sb_surv_final_prediction(
    c("p1", "p2", "p3"),
    risk = c(1, 2, 3),
    survival_probability = c(0.20, 0.20, 0.90),
    prediction_horizon = 3,
    time_unit = "month",
    decision_thresholds = 0.50
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = c("p1", "p2", "p3"),
      time = c(1, 2, 4),
      event = c(1L, 0L, 0L)
    ),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 3
  )
  incomplete <- evaluate_external_original(prediction, endpoint)

  expect_identical(incomplete$status, "incomplete")
  expect_identical(
    incomplete$metrics$decision_curve$curve$horizon_support,
    "insufficient_follow_up"
  )
})

test_that("external survival thresholds are optional and strictly validated", {
  ids <- c("p1", "p2", "p3")
  plain <- make_sb_surv_final_prediction(
    ids,
    risk = c(1, 2, 3),
    survival_probability = c(0.8, 0.5, 0.2),
    prediction_horizon = 2,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 2, 3),
      event = c(1L, 0L, 1L)
    ),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 2
  )
  plain_result <- evaluate_external_original(plain, endpoint)
  expect_false("decision_curve" %in% names(plain_result$metrics))

  invalid <- plain
  invalid$decision_thresholds <- c(0.20, 1)
  invalid_result <- evaluate_external_original(invalid, endpoint)
  expect_external_failure(
    invalid_result,
    "prediction_contract",
    "stormbreaker_bad_decision_thresholds"
  )

  risk_only <- make_sb_surv_final_prediction(
    ids,
    risk = c(1, 2, 3),
    decision_thresholds = 0.50
  )
  risk_endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 2, 3),
      event = c(1L, 0L, 1L)
    ),
    "id",
    "time",
    "event",
    "month"
  )
  risk_result <- evaluate_external_original(risk_only, risk_endpoint)
  expect_external_failure(
    risk_result,
    "prediction_contract",
    "stormbreaker_bad_external_prediction"
  )
})

# ---- external-performance ----

test_that("external performance preserves metric values and cohort identity", {
  skip_if_not_installed("ggplot2")
  for (type in c("binary", "survival")) {
    x <- make_sb_curve_view_fixture(type, comparison = TRUE)$external
    metrics <- if (type == "binary") c("roc_auc", "average_precision", "brier_score", "log_loss") else "harrell_c"
    for (metric in metrics) {
      p <- sb_plot_performance(x, metric, cohort_label = "Independent cohort")
      expect_identical(p$data$estimate, sb_result_table(x, metric = metric)$estimate)
      expect_identical(levels(p$data$pipeline_id), names(x$evaluations))
      expect_identical(p$labels$title, "External performance")
      expect_match(p$labels$subtitle, "Independent cohort", fixed = TRUE)
      expect_match(p$labels$subtitle, x$evaluations[[1]]$metrics[[metric]]$aggregation, fixed = TRUE)
      expect_length(p$layers, 1L)
      expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
      expect_equal(ggplot2::ggplot_build(p)$data[[1]]$y, p$data$estimate)
    }
  }
})

test_that("external failures and nonfinite values are retained but not drawn", {
  skip_if_not_installed("ggplot2")
  a <- make_sb_curve_view_fixture("binary")$external
  b <- a
  b$pipeline_id <- "failed"
  b$status <- "failed"
  b$metrics <- list()
  c <- a
  c$pipeline_id <- "infinite"
  c$metrics$log_loss$estimate <- Inf
  x <- structure(list(evaluations = list(b, a, c)), class = "stormbreaker_external_comparison")
  p <- sb_plot_performance(x, "log_loss")
  expect_identical(p$data$estimate, c(NA_real_, a$metrics$log_loss$estimate, Inf))
  expect_identical(nrow(p$layers[[1]]$data), 1L)
  expect_identical(levels(p$data$pipeline_id), c("failed", a$pipeline_id, "infinite"))
  expect_match(p$labels$caption, "2 pipeline", fixed = TRUE)
  x$evaluations <- list(b, c)
  empty <- sb_plot_performance(x, "log_loss", show_unavailable = TRUE)
  expect_identical(nrow(ggplot2::ggplot_build(empty)$data[[1]]), 0L)
  expect_identical(nrow(empty$data), 2L)
})

test_that("external performance never refits or evaluates and CV remains unchanged", {
  skip_if_not_installed("ggplot2")
  fixture <- make_sb_curve_view_fixture("survival")
  x <- fixture$external
  before <- x
  local_mocked_bindings(
    sb_fit_final = function(...) stop("must not fit"),
    sb_evaluate_external = function(...) stop("must not evaluate"),
    sb_external_harrell_c_metric = function(...) stop("must not compute")
  )
  p <- sb_plot_performance(x, "harrell_c")
  expect_identical(p$data$estimate, x$metrics$harrell_c$estimate)
  expect_identical(x, before)
  cv <- sb_plot_performance(fixture$cv, "harrell_c")
  expect_identical(cv$labels$title, "CV performance")
  expect_identical(names(sb_result_table(fixture$cv, metric = "harrell_c")), c("pipeline_id", "task_type", "metric", "estimate", "status", "run_status", "scope", "aggregation"))
})
})

# ---- external-inputs ----
local({
# ---- external-endpoint ----

test_that("binary endpoints retain only patient ID and outcome", {
  expect_identical(
    names(formals(sb_binary_endpoint)),
    c("data", "id_col", "outcome_col")
  )
  data <- data.frame(
    patient_id = c("e1", "e2", "e3"),
    outcome = c(1L, 1L, 0L),
    gene_a = c(10, 20, 30),
    note = c("a", "b", "c"),
    stringsAsFactors = FALSE
  )

  endpoint <- sb_binary_endpoint(data, "patient_id", "outcome")

  expect_s3_class(endpoint, "stormbreaker_binary_endpoint")
  expect_identical(names(endpoint$data), c("patient_id", "outcome"))
  expect_identical(endpoint$data, data[c("patient_id", "outcome")])
  expect_identical(endpoint$n, 3L)
  expect_identical(endpoint$n_positive, 2L)
  expect_identical(endpoint$n_negative, 1L)
  expect_identical(endpoint$positive_class, 1L)
  expect_identical(endpoint$negative_class, 0L)
  expect_identical(endpoint$independent_unit, "patient")
})

test_that("binary endpoints record one patient and one observed class", {
  one_positive <- data.frame(patient_id = "e1", outcome = 1L)
  one_negative_class <- data.frame(
    patient_id = c("e1", "e2"),
    outcome = c(0L, 0L)
  )

  positive <- sb_binary_endpoint(
    one_positive,
    "patient_id",
    "outcome"
  )
  negative <- sb_binary_endpoint(
    one_negative_class,
    "patient_id",
    "outcome"
  )

  expect_identical(positive$n, 1L)
  expect_identical(positive$n_positive, 1L)
  expect_identical(positive$n_negative, 0L)
  expect_identical(negative$n_positive, 0L)
  expect_identical(negative$n_negative, 2L)
})

test_that("logical binary endpoints retain explicit zero-one semantics", {
  data <- data.frame(
    patient_id = c("e1", "e2"),
    outcome = c(TRUE, FALSE)
  )

  endpoint <- sb_binary_endpoint(data, "patient_id", "outcome")

  expect_identical(endpoint$data$outcome, c(TRUE, FALSE))
  expect_identical(endpoint$n_positive, 1L)
  expect_identical(endpoint$n_negative, 1L)
})

test_that("survival endpoints retain only ID time and event semantics", {
  expect_identical(
    names(formals(sb_surv_endpoint)),
    c(
      "data", "id_col", "time_col", "event_col", "time_unit",
      "prediction_horizon"
    )
  )
  data <- data.frame(
    patient_id = c("e1", "e2", "e3"),
    follow_up = c(2, 5, 7),
    died = c(1L, 0L, 1L),
    gene_a = c(4, 5, 6),
    stringsAsFactors = FALSE
  )

  endpoint <- sb_surv_endpoint(
    data,
    "patient_id",
    "follow_up",
    "died",
    "month",
    prediction_horizon = 6
  )

  expect_s3_class(endpoint, "stormbreaker_surv_endpoint")
  expect_identical(
    names(endpoint$data),
    c("patient_id", "follow_up", "died")
  )
  expect_identical(
    endpoint$data,
    data[c("patient_id", "follow_up", "died")]
  )
  expect_identical(endpoint$n, 3L)
  expect_identical(endpoint$n_events, 2L)
  expect_identical(endpoint$n_censored, 1L)
  expect_identical(endpoint$time_unit, "month")
  expect_identical(endpoint$prediction_horizon, 6)
  expect_identical(endpoint$independent_unit, "patient")
})

test_that("survival endpoints allow one row all events and all censoring", {
  one <- data.frame(patient_id = "e1", time = 2, event = 1L)
  all_events <- data.frame(
    patient_id = paste0("e", 1:3),
    time = c(1, 2, 3),
    event = 1L
  )
  all_censored <- data.frame(
    patient_id = paste0("c", 1:3),
    time = c(4, 5, 6),
    event = 0L
  )

  one_endpoint <- sb_surv_endpoint(
    one, "patient_id", "time", "event", "day"
  )
  event_endpoint <- sb_surv_endpoint(
    all_events, "patient_id", "time", "event", "month"
  )
  censored_endpoint <- sb_surv_endpoint(
    all_censored, "patient_id", "time", "event", "year"
  )

  expect_identical(one_endpoint$n, 1L)
  expect_identical(event_endpoint$n_events, 3L)
  expect_identical(event_endpoint$n_censored, 0L)
  expect_identical(censored_endpoint$n_events, 0L)
  expect_identical(censored_endpoint$n_censored, 3L)
})

test_that("endpoint constructors require one nonambiguous endpoint row", {
  empty_binary <- data.frame(
    patient_id = character(),
    outcome = integer()
  )
  empty_error <- capture_sb_error(
    sb_binary_endpoint(empty_binary, "patient_id", "outcome")
  )
  expect_s3_class(empty_error, "stormbreaker_empty_endpoint")

  binary <- data.frame(patient_id = "e1", outcome = 1L)
  overlap_error <- capture_sb_error(
    sb_binary_endpoint(binary, "patient_id", "patient_id")
  )
  expect_s3_class(
    overlap_error,
    "stormbreaker_overlapping_endpoint_roles"
  )

  duplicate <- data.frame(
    patient_id = c("e1", "e2"),
    outcome = c(0L, 1L),
    second = c(1L, 0L),
    check.names = FALSE
  )
  names(duplicate)[3L] <- "outcome"
  duplicate_error <- capture_sb_error(
    sb_binary_endpoint(duplicate, "patient_id", "outcome")
  )
  expect_s3_class(
    duplicate_error,
    "stormbreaker_duplicate_endpoint_column"
  )

  survival <- data.frame(patient_id = "e1", time = 1, event = 1L)
  survival_role_error <- capture_sb_error(
    sb_surv_endpoint(
      survival, "patient_id", "time", "time", "month"
    )
  )
  expect_s3_class(
    survival_role_error,
    "stormbreaker_overlapping_endpoint_roles"
  )
})

test_that("endpoint patient identifiers remain complete and unique", {
  duplicate <- data.frame(
    patient_id = c("e1", "e1"),
    outcome = c(0L, 1L)
  )
  duplicate_error <- capture_sb_error(
    sb_binary_endpoint(duplicate, "patient_id", "outcome")
  )
  expect_s3_class(duplicate_error, "stormbreaker_duplicate_id")

  blank <- data.frame(patient_id = " ", outcome = 0L)
  blank_error <- capture_sb_error(
    sb_binary_endpoint(blank, "patient_id", "outcome")
  )
  expect_s3_class(blank_error, "stormbreaker_bad_id")
})

test_that("endpoint values and time declarations are explicit", {
  binary <- data.frame(patient_id = "e1", outcome = 2L)
  binary_error <- capture_sb_error(
    sb_binary_endpoint(binary, "patient_id", "outcome")
  )
  expect_s3_class(binary_error, "stormbreaker_bad_binary_outcome")

  survival <- data.frame(patient_id = "e1", time = 1, event = 1L)
  bad_time <- survival
  bad_time$time <- -1
  time_error <- capture_sb_error(
    sb_surv_endpoint(
      bad_time, "patient_id", "time", "event", "month"
    )
  )
  expect_s3_class(time_error, "stormbreaker_bad_time")

  bad_event <- survival
  bad_event$event <- 2L
  event_error <- capture_sb_error(
    sb_surv_endpoint(
      bad_event, "patient_id", "time", "event", "month"
    )
  )
  expect_s3_class(event_error, "stormbreaker_bad_event")

  unit_error <- capture_sb_error(
    sb_surv_endpoint(
      survival, "patient_id", "time", "event", "months"
    )
  )
  expect_s3_class(unit_error, "stormbreaker_bad_time_unit")

  for (horizon in list(0, -1, NA_real_, Inf, "1", c(1, 2))) {
    horizon_error <- capture_sb_error(
      sb_surv_endpoint(
        survival,
        "patient_id",
        "time",
        "event",
        "month",
        prediction_horizon = horizon
      )
    )
    expect_s3_class(
      horizon_error,
      "stormbreaker_bad_prediction_horizon"
    )
  }
})

test_that("endpoint contract validation detects changed endpoint metadata", {
  binary <- sb_binary_endpoint(
    data.frame(patient_id = c("e1", "e2"), outcome = c(0L, 1L)),
    "patient_id",
    "outcome"
  )
  binary$n_positive <- 2L
  binary_error <- capture_sb_error(
    sb_validate_binary_endpoint_contract(binary)
  )
  expect_s3_class(binary_error, "stormbreaker_mutated_endpoint")

  survival <- sb_surv_endpoint(
    data.frame(
      patient_id = c("e1", "e2"),
      time = c(1, 2),
      event = c(1L, 0L)
    ),
    "patient_id",
    "time",
    "event",
    "month"
  )
  survival$n_events <- 2L
  survival_error <- capture_sb_error(
    sb_validate_surv_endpoint_contract(survival)
  )
  expect_s3_class(survival_error, "stormbreaker_mutated_endpoint")
})

test_that("endpoint print methods describe outcome-only objects", {
  binary <- sb_binary_endpoint(
    data.frame(patient_id = "e1", outcome = 1L),
    "patient_id",
    "outcome"
  )
  survival <- sb_surv_endpoint(
    data.frame(patient_id = "e1", time = 1, event = 1L),
    "patient_id",
    "time",
    "event",
    "month",
    prediction_horizon = 2
  )

  binary_output <- paste(capture.output(print(binary)), collapse = "\n")
  survival_output <- paste(
    capture.output(print(survival)),
    collapse = "\n"
  )

  expect_match(binary_output, "outcome-only evaluation endpoint")
  expect_match(binary_output, "Patients: 1")
  expect_match(survival_output, "outcome-only evaluation endpoint")
  expect_match(survival_output, "Prediction horizon: 2 \\[ month \\]")
  expect_no_match(binary_output, "fitted")
  expect_no_match(survival_output, "validated")
})

# ---- external-result-table ----

test_that("external metric tables copy all supported saved scalar values", {
  for (type in c("binary", "survival")) {
    x <- make_sb_curve_view_fixture(type, comparison = TRUE)$external
    before <- x
    metrics <- if (type == "binary") c("roc_auc", "average_precision", "brier_score", "log_loss") else "harrell_c"
    for (metric in metrics) {
      tab <- sb_result_table(x, metric = metric, cohort_label = "Independent cohort")
      expect_identical(tab$pipeline_id, names(x$evaluations))
      expect_identical(tab$cohort_label, rep("Independent cohort", 2L))
      for (i in seq_along(x$evaluations)) {
        saved <- x$evaluations[[i]]
        expect_identical(tab$estimate[i], saved$metrics[[metric]]$estimate)
        expect_identical(tab$status[i], saved$metrics[[metric]]$status)
        expect_identical(tab$run_status[i], saved$status)
        expect_identical(tab$scope[i], saved$metrics[[metric]]$scope)
        expect_identical(tab$aggregation[i], saved$metrics[[metric]]$aggregation)
        expect_identical(tab$task_type[i], type)
      }
      single <- sb_result_table(x$evaluations[[1]], metric = metric)
      expect_identical(single$cohort_label, NA_character_)
      expect_identical(single$estimate, tab$estimate[1])
    }
    expect_identical(x, before)
  }
})

test_that("external predictions and evaluation runs preserve their distinct meaning", {
  for (type in c("binary", "survival")) {
    x <- make_sb_curve_view_fixture(type, comparison = TRUE)$external
    for (i in seq_along(x$evaluations)) {
      saved <- x$evaluations[[i]]
      tab <- sb_result_table(saved, "predictions", cohort_label = "Test cohort")
      expect_identical(tab[names(saved$scored_predictions)], saved$scored_predictions)
      expect_identical(tab$task_type, rep(type, nrow(tab)))
      expect_identical(tab$cohort_label, rep("Test cohort", nrow(tab)))
      runs <- sb_result_table(saved, "runs")
      expect_identical(nrow(runs), 1L)
      expect_identical(runs$status, saved$status)
      expect_identical(runs$patient_match_status, saved$patient_match$status)
      expect_identical(runs$error_stage, saved$error_stage)
      expect_identical(runs$error_class, saved$error_class)
      expect_identical(runs$error_message, saved$error_message)
      expect_identical(runs$n, saved$n)
      expect_identical("fold" %in% names(runs), FALSE)
    }
    expected <- do.call(rbind, lapply(unname(x$evaluations), `[[`, "scored_predictions"))
    rownames(expected) <- NULL
    tab <- sb_result_table(x, "predictions")
    expect_identical(tab[names(expected)], expected)
    expect_identical(sb_result_table(x, "runs")$pipeline_id, names(x$evaluations))
  }
})

test_that("real external evaluation failure remains named without invented metrics", {
  d <- sb_example_binary(120, seed = 1)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  ext <- sb_example_binary(100, seed = 2)
  pred <- sb_predict_final(fit, ext[c("patient_id", "gene_a", "gene_b")], id_col = "patient_id")
  ep <- sb_binary_endpoint(ext, "patient_id", "outcome")
  good <- sb_evaluate_external(pred, ep)
  pred$pipeline_id <- "failed_id_match"
  pred$predictions$pipeline_id <- "failed_id_match"
  pred$predictions$patient_id[1] <- "unknown_patient"
  failed <- sb_evaluate_external(pred, ep)
  expect_identical(failed$status, "failed")
  x <- structure(list(evaluations = list(failed, good)), class = "stormbreaker_external_comparison")
  tab <- sb_result_table(x, metric = "roc_auc")
  expect_identical(tab$pipeline_id, c("failed_id_match", good$pipeline_id))
  expect_identical(tab$estimate, c(NA_real_, good$metrics$roc_auc$estimate))
  expect_identical(tab$status, c("not_saved", good$metrics$roc_auc$status))
  expect_identical(tab$run_status, c("failed", good$status))
  expect_identical(tab$scope[1], NA_character_)
  expect_identical(tab$aggregation[1], NA_character_)
  expect_identical(tab$error_stage[1], failed$error_stage)
  expect_identical(tab$error_message[1], failed$error_message)
  expect_identical(nrow(sb_result_table(failed, "predictions")), 0L)
  expect_identical(nrow(sb_result_table(x, "predictions")), nrow(good$scored_predictions))
})

test_that("external display is read-only and rejects unsupported requests", {
  x <- make_sb_curve_view_fixture("binary")$external
  before <- x
  local_mocked_bindings(
    sb_fit_final = function(...) stop("must not fit"),
    sb_evaluate_external = function(...) stop("must not evaluate"),
    sb_external_binary_auc_metric = function(...) stop("must not compute")
  )
  expect_identical(sb_result_table(x, metric = "roc_auc")$estimate, x$metrics$roc_auc$estimate)
  expect_identical(x, before)
  expect_error(sb_result_table(x, metric = "harrell_c"), class = "stormbreaker_result_metric_not_saved")
  expect_error(sb_result_table(x, metric = "calibration"), class = "stormbreaker_bad_result_metric")
  expect_error(sb_result_table(x, "runs", metric = "roc_auc"), class = "stormbreaker_unused_result_metric")
  for (label in list("", "  ", NA_character_, c("A", "B"), 1)) {
    expect_error(sb_result_table(x, metric = "roc_auc", cohort_label = label), class = "stormbreaker_bad_cohort_label")
  }
  expect_error(sb_result_table(list(x), metric = "roc_auc"), class = "stormbreaker_bad_result_view")
  expect_error(sb_result_table(structure(list(evaluations = list()), class = "stormbreaker_external_comparison"), metric = "roc_auc"), class = "stormbreaker_bad_result_view")
})

# ---- prediction-variation ----

variation_comparison <- function() {
  make_task <- function(prefix, seed) {
    d <- sb_example_binary(60, seed)
    d$patient_id <- paste0(prefix, d$patient_id)
    sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  }
  sb_compare_cohorts(make_task("dev", 1),
    list(Validation = make_task("val", 2), Test = make_task("test", 3)),
    roles = c(Validation = "validation", Test = "test"),
    models = c("logistic", "gaussian_naive_bayes"), folds = 3,
    view = FALSE, progress = FALSE)
}

test_that("variation distinguishes exact ties, small differences and singleton batches", {
  rows <- data.frame(status = "success", risk = c(1, 1))
  expect_identical(sb_prediction_variation(rows, "survival")$status, "constant")
  rows$risk[2] <- 1 + 1e-12
  expect_identical(sb_prediction_variation(rows, "survival")$status, "variable")
  expect_identical(sb_prediction_variation(rows[1, ], "survival")$status,
    "insufficient_patients")
  rows$risk[2] <- NA_real_
  expect_identical(sb_prediction_variation(rows, "survival")$status, "incomplete")
})

test_that("different constants in different CV folds do not count as discrimination", {
  rows <- data.frame(status = "success", probability = rep(c(.2, .4, .6), each = 3),
    fold = rep(1:3, each = 3))
  diagnostic <- sb_prediction_variation(rows, "binary")
  expect_identical(diagnostic$status, "constant")
  expect_identical(diagnostic$constant_folds, 1:3)
  rows$probability[1] <- .3
  expect_identical(sb_prediction_variation(rows, "binary")$status, "partly_constant")
})

test_that("constant validation predictions cannot beat a worse-than-chance model", {
  compared <- variation_comparison()
  ids <- unname(vapply(compared$model_pipelines, function(x) x$id, ""))
  compared$cohorts$Validation$comparison$roc_auc <- c(.5, .4)
  compared$cohorts$Validation$evaluations[[ids[1]]]$scored_predictions$probability <- .5
  selected <- sb_select_model(compared)
  expect_identical(selected$pipeline$id, ids[2])
  ranking <- selected$model_selection$ranking
  expect_identical(ranking$eligible[match(ids[1], ranking$pipeline_id)], FALSE)
  expect_match(ranking$exclusion_reason[match(ids[1], ranking$pipeline_id)],
    "constant_validation_predictions:Validation", fixed = TRUE)
  expect_equal(compared$cohorts$Validation$comparison$roc_auc[1], .5)
})

test_that("constant development folds are ineligible without reading test predictions", {
  compared <- variation_comparison()
  ids <- unname(vapply(compared$model_pipelines, function(x) x$id, ""))
  compared$cohorts$Validation$comparison$roc_auc <- c(.9, .6)
  rows <- compared$development$results[[ids[1]]]$predictions
  rows$probability <- rows$fold / 10
  compared$development$results[[ids[1]]]$predictions <- rows
  selected <- sb_select_model(compared)
  expect_identical(selected$pipeline$id, ids[2])
  compared$cohorts$Test$evaluations[[ids[2]]]$scored_predictions$probability <- .5
  expect_identical(sb_select_model(compared)$pipeline$id, ids[2])
  compared$cohorts$Validation$evaluations[[ids[2]]]$scored_predictions$probability <- .5
  error <- tryCatch(sb_select_model(compared), error = identity)
  expect_s3_class(error, "stormbreaker_no_eligible_model")
  expect_match(conditionMessage(error), "constant", fixed = TRUE)
})

test_that("constant final predictions stay available and explicitly labelled", {
  compared <- variation_comparison()
  fit <- compared$fits[[1]]
  d <- sb_example_binary(10, 4)
  d$gene_a <- 0
  d$gene_b <- 0
  prediction <- sb_predict_final(fit, d, "patient_id")
  expect_identical(prediction$status, "success")
  expect_identical(prediction$prediction_variation$status, "constant")
  endpoint <- sb_binary_endpoint(d, "patient_id", "outcome")
  external <- sb_compare_external(list(prediction), endpoint)
  expect_identical(external$comparison$prediction_variation, "constant")
  expect_identical(external$comparison$roc_auc_status, "success")
  expect_equal(external$comparison$roc_auc, .5)
})
})

# ---- external-direction ----
local({
# ---- external-reverse ----

reverse_fixture <- function(probability, outcome = c(0L, 0L, 1L, 1L)) {
  ids <- paste0("p", seq_along(outcome))
  prediction <- structure(list(
    status = "success", task_type = "binary", pipeline_id = "test",
    learner = "logistic", id_col = "id",
    predictions = data.frame(
      patient_id = ids, pipeline_id = "test", learner = "logistic",
      positive_class = 1L, probability = probability, status = "success"
    )
  ), class = "stormbreaker_final_prediction")
  endpoint <- sb_binary_endpoint(
    data.frame(id = rev(ids), outcome = rev(outcome)), "id", "outcome"
  )
  list(prediction = prediction, endpoint = endpoint)
}

test_that("reverse flips only low AUC and retains aligned originals", {
  x <- reverse_fixture(c(.9, .8, .2, .1))
  before <- x$prediction
  original <- sb_evaluate_external(x$prediction, x$endpoint, reverse = FALSE)
  off <- sb_evaluate_external(x$prediction, x$endpoint, reverse = FALSE)
  on <- sb_evaluate_external(x$prediction, x$endpoint, reverse = TRUE)
  expect_identical(sb_evaluate_external(x$prediction, x$endpoint), off)
  expect_identical(eval_external(x$prediction, x$endpoint), off)
  expect_identical(formals(sb_evaluate_external)$reverse, FALSE)
  expect_identical(formals(eval_external)$reverse, FALSE)
  expect_null(off$direction_correction)
  expect_equal(off$metrics$roc_auc$estimate, 0)
  expect_identical(off, original)
  expect_identical(x$prediction, before)
  expect_identical(on$original_metrics, original$metrics)
  expect_identical(on$original_scored_predictions, original$scored_predictions)
  expect_equal(on$direction_correction$auc_original, 0)
  expect_equal(on$direction_correction$auc_processed, 1)
  expect_identical(on$direction_correction$reversed, TRUE)
  expect_equal(on$scored_predictions$probability,
               1 - original$scored_predictions$probability)
  expect_equal(on$metrics$brier_score$estimate,
    mean((on$scored_predictions$outcome - on$scored_predictions$probability)^2))
  expect_identical(on$evaluation_boundary$outcome_use,
                  "metrics_and_post_hoc_direction_selection")
})

test_that("good, tied and unavailable AUC do not reverse", {
  for (p in list(c(.1, .2, .8, .9), rep(.5, 4), c(.2, .8, .2, .8))) {
    x <- reverse_fixture(p)
    result <- sb_evaluate_external(x$prediction, x$endpoint, reverse = TRUE)
    expect_identical(result$direction_correction$reversed, FALSE)
    expect_equal(result$scored_predictions$probability,
                 result$original_scored_predictions$probability)
  }
  x <- reverse_fixture(c(.1, .2, .3, .4), rep(1L, 4))
  result <- sb_evaluate_external(x$prediction, x$endpoint, reverse = TRUE)
  expect_identical(result$direction_correction$status, "not_evaluable")
  expect_identical(result$direction_correction$reversed, FALSE)
})

test_that("reverse validates its flag and preserves prediction failures", {
  x <- reverse_fixture(c(.1, .2, .8, .9))
  for (flag in list(NA, 1, "TRUE", logical(), c(TRUE, FALSE))) {
    error <- tryCatch(
      sb_evaluate_external(x$prediction, x$endpoint, reverse = flag),
      error = identity
    )
    expect_s3_class(error, "stormbreaker_bad_reverse")
  }
  result <- sb_evaluate_external(list(), x$endpoint, reverse = TRUE)
  expect_identical(result$status, "failed")
})

test_that("tied scores retain half credit after reversal", {
  x <- reverse_fixture(c(.8, .8, .8, .2))
  result <- sb_evaluate_external(x$prediction, x$endpoint, reverse = TRUE)
  expect_equal(result$direction_correction$auc_original, .25)
  expect_equal(result$direction_correction$auc_processed, .75)
})

test_that("omitted reverse preserves survival evaluation", {
  x <- reverse_fixture(c(.8, .7, .3, .2))
  p <- x$prediction
  p$task_type <- "survival"
  p$predictions$risk <- c(.8, .7, .3, .2)
  e <- sb_surv_endpoint(data.frame(id = paste0("p", 1:4),
    time = 1:4, event = c(1, 0, 1, 0)), "id", "time", "event", "day")
  expect_identical(sb_evaluate_external(p, e),
                   sb_evaluate_external(p, e, reverse = FALSE))
  expect_s3_class(tryCatch(sb_evaluate_external(p, e, reverse = TRUE),
    error = identity), "stormbreaker_bad_reverse")
})

# ---- short-names ----

test_that("short names are exact aliases with unchanged arguments", {
  pairs <- c(
    models = "sb_models", bi_task = "sb_binary_task", surv_task = "sb_surv_task",
    multi_task = "sb_multiclass_task", bi_endpoint = "sb_binary_endpoint",
    surv_endpoint = "sb_surv_endpoint", run_cv = "sb_run_cv", compare_cv = "sb_compare_cv",
    fit_final = "sb_fit_final", predict_final = "sb_predict_final",
    eval_external = "sb_evaluate_external", compare_external = "sb_compare_external",
    finish_run = "sb_complete_run", clinical_plots = "sb_clinician_views", plots = "sb_plots",
    result_table = "sb_result_table", selected_genes = "sb_selected_features",
    export_run = "sb_export_run", clinical_report = "sb_render_clinician_report",
    model_shap = "sb_model_shap", model_coefs = "sb_model_coefficients",
    bi_example = "sb_example_binary", surv_example = "sb_example_survival",
    plot_roc = "sb_plot_roc", plot_pr = "sb_plot_pr", plot_calib = "sb_plot_calibration",
    plot_dca = "sb_plot_dca", plot_perf = "sb_plot_performance",
    plot_heatmap = "sb_plot_performance_heatmap", plot_confusion = "sb_plot_confusion",
    plot_prob = "sb_plot_probability_distribution", plot_survival = "sb_plot_survival_curves",
    plot_nomogram = "sb_plot_nomogram", plot_coefs = "sb_plot_coefficients",
    plot_shap = "sb_plot_shap_beeswarm", plot_shap_importance = "sb_plot_shap_importance",
    plot_shap_dependence = "sb_plot_shap_dependence", plot_shap_waterfall = "sb_plot_shap_waterfall",
    plot_rank = "sb_plot_model_ranking", plot_cohorts = "sb_plot_cohort_metrics",
    plot_panel = "sb_plot_clinician_panel", plot_pages = "sb_plot_pages",
    plot_forest = "sb_plot_meta_forest", plot_selector = "sb_plot_selector_diagnostic",
    plot_algorithm = "sb_plot_algorithm_diagnostic"
  )
  exports <- getNamespaceExports("stormbreaker")
  expect_length(pairs, 45L)
  expect_contains(exports, names(pairs))
  expect_contains(exports, unname(pairs))
  for (short in names(pairs)) {
    old <- getExportedValue("stormbreaker", pairs[[short]])
    new <- getExportedValue("stormbreaker", short)
    expect_identical(new, old, info = short)
    expect_identical(formals(new), formals(old), info = short)
    expect_lt(nchar(short), nchar(pairs[[short]]))
  }
})

test_that("short workflow preserves predictions metrics and clinical coordinates", {
  skip_if_not_installed("ggplot2")
  development <- bi_example(120, seed = 42)
  validation <- bi_example(80, seed = 17)
  task <- bi_task(development, "patient_id", "outcome", c("gene_a", "gene_b"))
  expect_identical(task, sb_binary_task(development, "patient_id", "outcome", c("gene_a", "gene_b")))
  run <- bi_run(task, models = "Logistic", folds = 3, plots = FALSE)
  old <- sb_complete_run(run, "Logistic", external = list(validation = validation), plots = FALSE)
  new <- finish_run(run, "Logistic", external = list(validation = validation), plots = FALSE)
  expect_identical(new$predictions, old$predictions)
  expect_identical(new$metrics, old$metrics)
  expect_identical(result_table(new$evaluations$validation, metric = "roc_auc"),
    sb_result_table(old$evaluations$validation, metric = "roc_auc"))
  expect_null(new$direction_correction$validation)
  expect_identical(clinical_plots(new)$plots$ROC$data, sb_clinician_views(old)$plots$ROC$data)
  expect_identical(plot_roc(new)$data, sb_plot_roc(old)$data)
})

test_that("short survival workflow preserves raw risk and direction defaults", {
  data <- surv_example(80, seed = 42)
  genes <- c("gene_a", "gene_b")
  old_task <- sb_surv_task(data, "patient_id", "time", "event", genes, "month")
  new_task <- surv_task(data, "patient_id", "time", "event", genes, "month")
  expect_identical(new_task, old_task)
  route <- surv_pipeline("Cox")
  old_fit <- sb_fit_final(old_task, route)
  new_fit <- fit_final(new_task, route)
  expect_identical(new_fit$status, "success")
  query <- data[c("patient_id", genes)]
  old_prediction <- sb_predict_final(old_fit, query, "patient_id")
  new_prediction <- predict_final(new_fit, query, "patient_id")
  expect_identical(new_prediction$predictions, old_prediction$predictions)
  endpoint <- surv_endpoint(data, "patient_id", "time", "event", "month")
  old_evaluation <- sb_evaluate_external(old_prediction, endpoint)
  new_evaluation <- eval_external(new_prediction, endpoint)
  expect_identical(new_evaluation$metrics, old_evaluation$metrics)
  expect_null(new_evaluation$direction_correction)
})

test_that("short task names preserve validation errors", {
  data <- sb_example_binary(30, seed = 42)
  data$outcome[1] <- 2
  old <- tryCatch(sb_binary_task(data, "patient_id", "outcome", c("gene_a", "gene_b")), error = identity)
  new <- tryCatch(bi_task(data, "patient_id", "outcome", c("gene_a", "gene_b")), error = identity)
  expect_s3_class(old, "error")
  expect_s3_class(new, "error")
  expect_identical(class(new), class(old))
  expect_identical(conditionMessage(new), conditionMessage(old))
})

# ---- specialty-short-names ----

test_that("every specialty short name preserves the existing function and arguments", {
  pairs <- read.csv(system.file("guides", "specialty-name-map.csv",
    package = "stormbreaker", mustWork = TRUE), stringsAsFactors = FALSE)
  exports <- getNamespaceExports("stormbreaker")
  expect_gt(nrow(pairs), 0L)
  expect_contains(exports, pairs$short)
  expect_contains(exports, pairs$original)
  for (i in seq_len(nrow(pairs))) {
    short <- getExportedValue("stormbreaker", pairs$short[i])
    old <- getExportedValue("stormbreaker", pairs$original[i])
    expect_identical(short, old, info = pairs$short[i])
    expect_identical(formals(short), formals(old), info = pairs$short[i])
    expect_lt(nchar(pairs$short[i]), nchar(pairs$original[i]))
  }
})

test_that("all diagnostic and drawing exports have a preferred short name", {
  preferred <- read.csv(system.file("guides", "preferred-name-map.csv",
    package = "stormbreaker", mustWork = TRUE), stringsAsFactors = FALSE)
  exports <- getNamespaceExports("stormbreaker")
  scoped <- grep("^sb_plot_|_diagnostics$", exports, value = TRUE)
  expect_contains(preferred$original, scoped)
  expect_identical(anyDuplicated(preferred$short), 0L)
  expect_identical(plot_diag, plot_algorithm)
  expect_identical(plot_shap_imp, plot_shap_importance)
  expect_identical(plot_shap_dep, plot_shap_dependence)
  expect_identical(plot_shap_one, plot_shap_waterfall)
})

test_that("short Lasso interfaces preserve native tables and coordinates", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ggplot2")
  data <- bi_example(80, seed = 42)
  task <- bi_task(data, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  fit <- fit_final(task, bi_pipeline("Lasso"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  cv <- path_diag(fit, "cv_curve")
  path <- path_diag(fit, "coefficient_path")
  expect_identical(cv, sb_algorithm_diagnostics(fit, "cv_curve"))
  expect_identical(path, sb_algorithm_diagnostics(fit, "coefficient_path"))
  expect_gt(nrow(path$table), 0L)
  expect_identical(plot_diag(cv)$data, sb_plot_algorithm_diagnostic(cv)$data)
  expect_identical(native_diag(fit, complete = TRUE)$status,
    sb_native_diagnostics(fit, complete = TRUE)$status)
})

test_that("short forest interface preserves original importance", {
  skip_if_not_installed("randomForest")
  skip_if_not_installed("ggplot2")
  data <- bi_example(80, seed = 42)
  task <- bi_task(data, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  fit <- fit_final(task, bi_pipeline("RF"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  diagnostic <- forest_diag(fit, "importance")
  expect_identical(diagnostic, sb_forest_diagnostics(fit, "importance"))
  expect_equal(nrow(diagnostic$table), length(fit$required_features))
  expect_identical(plot_diag(diagnostic)$data, sb_plot_algorithm_diagnostic(diagnostic)$data)
})

test_that("short PAMR interface preserves both centroid stages", {
  skip_if_not_installed("pamr")
  skip_if_not_installed("ggplot2")
  data <- bi_example(80, seed = 42)
  task <- bi_task(data, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  fit <- fit_final(task, bi_pipeline("pamr"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  diagnostic <- pamr_diag(fit)
  expect_identical(diagnostic, sb_pamr_diagnostics(fit))
  expect_setequal(diagnostic$table$stage, c("Before", "After"))
  expect_equal(nrow(diagnostic$table), length(fit$required_features) * 4L)
  expect_identical(plot_diag(diagnostic)$data, sb_plot_algorithm_diagnostic(diagnostic)$data)
})

test_that("short diagnostic interfaces preserve missing-capture failures", {
  data <- bi_example(60, seed = 42)
  task <- bi_task(data, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- fit_final(task, bi_pipeline("Logistic"))
  old <- tryCatch(sb_algorithm_diagnostics(fit), error = identity)
  new <- tryCatch(path_diag(fit), error = identity)
  expect_s3_class(new, "stormbreaker_missing_diagnostic_artifacts")
  expect_identical(class(new), class(old))
  expect_identical(conditionMessage(new), conditionMessage(old))
  expect_null(fit$diagnostic_artifacts)
})

# ---- stable-defaults ----

test_that("new defaults keep historical pipeline settings readable", {
  nn <- sb_pipeline(learner = "neuralnet_multilayer_binary")
  expect_identical(nn$learner_params$hidden, c(3L, 2L))
  expect_identical(nn$learner_params$error_function, "sse")
  expect_identical(nn$learner_params$stepmax, 500000L)
  old <- nn
  old$learner_params$hidden <- c(5L, 3L)
  old$learner_params$error_function <- "ce"
  old$learner_params$stepmax <- 100000L
  expect_identical(sb_validate_pipeline_contract(old, "binary"), old)
  fda <- sb_pipeline(learner = "fda_binary")
  expect_identical(fda$learner_params$degree, 1L)
  fda$learner_params$degree <- 2L
  expect_identical(sb_validate_pipeline_contract(fda, "binary"), fda)
  lda <- sb_pipeline(learner = "loclda_binary")
  expect_identical(lda$learner_params$k_rule, "dimension_aware")
  lda$learner_params$k_rule <- NULL
  expect_identical(sb_validate_pipeline_contract(lda, "binary"), lda)
  svm <- sb_pipeline(learner = "survivalsvm_rbf")
  expect_identical(svm$learner_params$kernel_parameter_rule, "analysis_median_squared_distance")
  svm$learner_params$kernel_parameter_rule <- "native_inverse_analysis_feature_count"
  expect_identical(sb_validate_pipeline_contract(svm, "survival"), svm)
  svm$learner_params$opt_meth <- "unregistered"
  expect_s3_class(tryCatch(sb_validate_pipeline_contract(svm, "survival"), error = identity),
    "stormbreaker_unresolved_pipeline")
})

test_that("local LDA neighbourhood accounts for dimension and class count", {
  p <- sb_classical_binary_params("loclda_binary")
  expect_equal(sb_loclda_neighbors(p, 3L, 100L, 2L), 30L)
  expect_equal(sb_loclda_neighbors(p, 60L, 400L, 2L), 124L)
  expect_equal(sb_loclda_neighbors(p, 60L, 400L, 3L), 126L)
  expect_equal(sb_loclda_neighbors(p, 60L, 80L, 2L), 80L)
  p$k_rule <- NULL
  expect_equal(sb_loclda_neighbors(p, 60L, 400L, 2L), 30L)
})

# ---- unified-reverse-workflow ----

merge_reverse_fixture <- function() {
  ids <- paste0("p", 1:8)
  data <- data.frame(id = ids, gene = seq_along(ids), y = rep(0:1, each = 4))
  prediction <- structure(list(
    status = "success", task_type = "binary", pipeline_id = "merge-test",
    learner = "logistic", id_col = "id",
    predictions = data.frame(patient_id = ids, pipeline_id = "merge-test",
      learner = "logistic", positive_class = 1L,
      probability = c(.9, .8, .7, .6, .4, .3, .2, .1), status = "success")
  ), class = "stormbreaker_final_prediction")
  list(data = data, prediction = prediction,
    endpoint = sb_binary_endpoint(data[8:1, ], "id", "y"))
}

test_that("merged batch comparison retains raw direction and variation diagnostics", {
  x <- merge_reverse_fixture()
  point <- sb_evaluate_external(x$prediction, x$endpoint)
  batch <- sb_compare_external(list(x$prediction), x$endpoint)
  expect_identical(formals(sb_evaluate_external)$reverse, FALSE)
  expect_equal(point$metrics$roc_auc$estimate, 0)
  expect_null(point$direction_correction)
  corrected <- sb_evaluate_external(x$prediction, x$endpoint, reverse = TRUE)
  expect_equal(corrected$metrics$roc_auc$estimate, 1)
  expect_equal(batch$comparison$roc_auc, 0)
  expect_identical(batch$comparison$prediction_variation,
    sb_prediction_variation(x$prediction$predictions, "binary")$status)
  expect_null(batch$evaluations[[1]]$direction_correction)
  x$prediction$predictions$probability[] <- .5
  tied <- sb_compare_external(list(x$prediction), x$endpoint)
  expect_identical(tied$comparison$prediction_variation, "constant")
  expect_equal(tied$comparison$roc_auc, .5)
})

test_that("completion and export keep raw and processed values explicitly paired", {
  x <- merge_reverse_fixture()
  run <- list(
    task_type = "binary", pipelines = list(chosen = sb_pipeline(learner = "logistic")),
    development_task = sb_binary_task(x$data, "id", "y", "gene"),
    comparison = list(comparison = data.frame(model = "chosen")),
    details = list(chosen = list(fit = list(status = "success", selection = list(state = NULL))))
  )
  local_mocked_bindings(
    sb_predict_final = function(...) x$prediction,
    sb_fit_final = function(...) stop("Unexpected refit"),
    sb_compare_cv = function(...) stop("Unexpected model search"))
  original_run <- run
  result <- sb_complete_run(run, "chosen", list(validation = x$data[8:1, ]),
    plots = FALSE, reverse = TRUE)
  expect_identical(run, original_run)
  expect_identical(result$predictions$validation, x$prediction$predictions)
  expect_identical(result$direction_correction$validation$reversed, TRUE)
  evaluated <- result$evaluated_predictions$validation
  expect_equal(result$metrics$validation$roc_auc$estimate, 1)
  expect_equal(evaluated$probability, 1 - evaluated$probability_original)
  expect_identical(evaluated$patient_id, rev(x$prediction$predictions$patient_id))
  raw <- sb_complete_run(run, "chosen", list(validation = x$data), plots = FALSE, reverse = FALSE)
  expect_identical(sb_complete_run(run, "chosen", list(validation = x$data),
    plots = FALSE), raw)
  expect_identical(finish_run(run, "chosen", list(validation = x$data),
    plots = FALSE), raw)
  expect_identical(finish_run(run, "chosen", list(validation = x$data[8:1, ]),
    plots = FALSE, reverse = TRUE), result)
  expect_equal(raw$metrics$validation$roc_auc$estimate, 0)
  expect_null(raw$direction_correction$validation)
  path <- tempfile("merged-export-")
  sb_export_run(result, path)
  correction <- read.csv(file.path(path, "cohort-001-direction-correction.csv"))
  scored <- read.csv(file.path(path, "cohort-001-evaluated-predictions.csv"))
  saved <- read.csv(file.path(path, "cohort-001-predictions.csv"))
  original_metrics <- read.csv(file.path(path, "cohort-001-original-metrics.csv"))
  expect_identical(correction$reversed, TRUE)
  expect_equal(correction$auc_original, 0)
  expect_equal(correction$auc_processed, 1)
  expect_equal(scored$probability, evaluated$probability)
  expect_equal(saved$probability, x$prediction$predictions$probability)
  expect_equal(original_metrics$estimate[original_metrics$metric == "roc_auc"], 0)
})

test_that("reverse exports still expose the expanded APIs and complete menus", {
  expected <- c("bi_run", "surv_run", "multi_run", "sb_binary_menu_pipelines",
    "sb_survival_menu_pipelines", "sb_selected_features", "sb_algorithm_diagnostics",
    "sb_compare_multiclass", "sb_model_shap", "sb_complete_run", "sb_export_run")
  expect_identical(all(expected %in% getNamespaceExports("stormbreaker")), TRUE)
  expect_identical(formals(sb_complete_run)$reverse, FALSE)
  expect_identical(formals(finish_run)$reverse, FALSE)
})
})

# ---- external-comparison ----
local({
# ---- bootstrap-external-comparison ----

make_sb_paired_binary_prediction <- function(patient_ids, probability,
                                             pipeline_id) {
  n <- length(patient_ids)
  result <- list(
    status = "success",
    task_type = "binary",
    pipeline_id = pipeline_id,
    learner = "logistic",
    id_col = "patient_id",
    predictions = data.frame(
      patient_id = patient_ids,
      pipeline_id = rep(pipeline_id, n),
      learner = rep("logistic", n),
      positive_class = rep.int(1L, n),
      probability = probability,
      status = rep("success", n),
      error_class = rep(NA_character_, n),
      error_message = rep(NA_character_, n),
      stringsAsFactors = FALSE
    ),
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )
  class(result) <- "stormbreaker_final_prediction"
  result
}

make_sb_paired_survival_prediction <- function(
    patient_ids,
    risk,
    survival_probability,
    pipeline_id) {
  n <- length(patient_ids)
  result <- list(
    status = "success",
    task_type = "survival",
    pipeline_id = pipeline_id,
    learner = "coxph",
    id_col = "patient_id",
    predictions = data.frame(
      patient_id = patient_ids,
      pipeline_id = rep(pipeline_id, n),
      learner = rep("coxph", n),
      risk = risk,
      prediction_horizon = rep(30, n),
      time_unit = rep("month", n),
      survival_probability = survival_probability,
      status = rep("success", n),
      error_class = rep(NA_character_, n),
      error_message = rep(NA_character_, n),
      stringsAsFactors = FALSE
    ),
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )
  class(result) <- "stormbreaker_final_prediction"
  result
}

test_that("binary model differences use the same patient bootstrap draws", {
  ids <- sprintf("p%02d", seq_len(40L))
  outcome <- rep(c(0L, 1L), 20L)
  reference_probability <- seq(0.08, 0.92, length.out = length(ids))
  candidate_probability <- pmin(
    0.97,
    pmax(0.03, reference_probability + rep(c(-0.04, 0.06), 20L))
  )
  predictions <- list(
    make_sb_paired_binary_prediction(
      ids, reference_probability, "reference_binary"
    ),
    make_sb_paired_binary_prediction(
      ids, candidate_probability, "candidate_binary"
    ),
    make_sb_paired_binary_prediction(
      ids, rev(candidate_probability), "candidate_binary_2"
    )
  )
  endpoint <- sb_binary_endpoint(
    data.frame(patient_id = ids, outcome = outcome),
    "patient_id",
    "outcome"
  )
  comparison <- sb_compare_external(predictions, endpoint)

  set.seed(9400L)
  caller_state <- .Random.seed
  observed <- sb_bootstrap_external_comparison(
    comparison,
    reference_pipeline_id = "reference_binary",
    reps = 30L,
    confidence_level = 0.90,
    seed = 9401L
  )
  expect_identical(.Random.seed, caller_state)

  reference <- sb_bootstrap_external(
    comparison$evaluations$reference_binary,
    reps = 30L,
    confidence_level = 0.90,
    seed = 9401L
  )
  candidate <- sb_bootstrap_external(
    comparison$evaluations$candidate_binary,
    reps = 30L,
    confidence_level = 0.90,
    seed = 9401L
  )
  reference_rows <- reference$replicates[
    reference$replicates$metric %in% observed$metrics,
    ,
    drop = FALSE
  ]
  candidate_rows <- candidate$replicates[
    candidate$replicates$metric %in% observed$metrics,
    ,
    drop = FALSE
  ]
  expected_key <- paste(reference_rows$replicate, reference_rows$metric)
  candidate_rows <- candidate_rows[match(
    expected_key,
    paste(candidate_rows$replicate, candidate_rows$metric)
  ), , drop = FALSE]
  difference_rows <- observed$differences[
    observed$differences$candidate_pipeline_id == "candidate_binary",
    ,
    drop = FALSE
  ]
  difference_rows <- difference_rows[match(
    expected_key,
    paste(difference_rows$replicate, difference_rows$metric)
  ), , drop = FALSE]

  expect_s3_class(
    observed,
    "stormbreaker_external_bootstrap_comparison"
  )
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$reference_pipeline_id, "reference_binary")
  expect_identical(observed$prediction_target, list(
    target_status = "complete",
    has_horizon = FALSE,
    prediction_horizon = NULL,
    time_unit = NULL,
    horizon_source = "not_applicable"
  ))
  expect_identical(
    observed$metrics,
    c("roc_auc", "average_precision", "brier_score", "log_loss")
  )
  expect_identical(observed$resample_audit, reference$resample_audit)
  expect_equal(
    difference_rows$reference_estimate,
    reference_rows$estimate,
    tolerance = 0
  )
  expect_equal(
    difference_rows$candidate_estimate,
    candidate_rows$estimate,
    tolerance = 0
  )
  expect_equal(
    difference_rows$estimate,
    candidate_rows$estimate - reference_rows$estimate,
    tolerance = 0
  )
  brier_differences <- difference_rows$estimate[
    difference_rows$metric == "brier_score"
  ]
  brier_interval <- observed$intervals[
    observed$intervals$candidate_pipeline_id == "candidate_binary" &
      observed$intervals$metric == "brier_score",
    ,
    drop = FALSE
  ]
  expect_equal(
    unname(c(brier_interval$lower, brier_interval$upper)),
    unname(stats::quantile(
      brier_differences,
      probs = c(0.05, 0.95),
      names = FALSE,
      type = 7
    )),
    tolerance = 0
  )
  expect_identical(
    brier_interval$metric_direction,
    "lower_is_better"
  )
  second_candidate <- observed$differences[
    observed$differences$candidate_pipeline_id == "candidate_binary_2",
    ,
    drop = FALSE
  ]
  second_candidate <- second_candidate[match(
    expected_key,
    paste(second_candidate$replicate, second_candidate$metric)
  ), , drop = FALSE]
  expect_equal(
    second_candidate$reference_estimate,
    difference_rows$reference_estimate,
    tolerance = 0
  )
  expect_identical(
    observed$candidate_status$status,
    c("success", "success")
  )
  expect_identical(
    any(grepl("calibration|decision|rank|winner", observed$metrics)),
    FALSE
  )
})

test_that("survival horizon differences preserve paired censoring draws", {
  ids <- sprintf("s%02d", seq_len(60L))
  time <- as.numeric(seq_len(60L))
  event <- rep.int(1L, length(ids))
  reference_risk <- seq(-1.2, 1.2, length.out = length(ids))
  candidate_risk <- reference_risk + sin(seq_along(ids)) / 8
  predictions <- list(
    make_sb_paired_survival_prediction(
      ids,
      reference_risk,
      stats::plogis(-reference_risk),
      "reference_survival"
    ),
    make_sb_paired_survival_prediction(
      ids,
      candidate_risk,
      stats::plogis(-candidate_risk),
      "candidate_survival"
    )
  )
  endpoint <- sb_surv_endpoint(
    data.frame(patient_id = ids, time = time, event = event),
    "patient_id",
    "time",
    "event",
    "month",
    prediction_horizon = 30
  )
  comparison <- sb_compare_external(predictions, endpoint)

  observed <- sb_bootstrap_external_comparison(
    comparison,
    reference_pipeline_id = "reference_survival",
    reps = 20L,
    confidence_level = 0.90,
    seed = 9501L
  )
  reference <- sb_bootstrap_external(
    comparison$evaluations$reference_survival,
    reps = 20L,
    confidence_level = 0.90,
    seed = 9501L
  )
  candidate <- sb_bootstrap_external(
    comparison$evaluations$candidate_survival,
    reps = 20L,
    confidence_level = 0.90,
    seed = 9501L
  )
  reference_rows <- reference$replicates[
    reference$replicates$metric %in% observed$metrics,
    ,
    drop = FALSE
  ]
  candidate_rows <- candidate$replicates[
    candidate$replicates$metric %in% observed$metrics,
    ,
    drop = FALSE
  ]
  key <- paste(reference_rows$replicate, reference_rows$metric)
  candidate_rows <- candidate_rows[match(
    key,
    paste(candidate_rows$replicate, candidate_rows$metric)
  ), , drop = FALSE]
  difference_rows <- observed$differences[match(
    key,
    paste(observed$differences$replicate, observed$differences$metric)
  ), , drop = FALSE]

  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "survival")
  expect_identical(observed$prediction_target, list(
    target_status = "complete",
    has_horizon = TRUE,
    prediction_horizon = 30,
    time_unit = "month",
    horizon_source = "endpoint"
  ))
  expect_identical(
    observed$metrics,
    c("harrell_c", "ipcw_brier_score", "time_auc", "uno_c")
  )
  expect_identical(observed$resample_audit, reference$resample_audit)
  expect_equal(
    difference_rows$estimate,
    candidate_rows$estimate - reference_rows$estimate,
    tolerance = 0
  )
  expect_identical(
    observed$intervals$metric_direction,
    c(
      "higher_is_better", "lower_is_better",
      "higher_is_better", "higher_is_better"
    )
  )
})

test_that("a failed candidate cannot erase valid paired comparisons", {
  ids <- sprintf("f%02d", seq_len(40L))
  outcome <- rep(c(0L, 1L), 20L)
  base_probability <- seq(0.10, 0.90, length.out = length(ids))
  predictions <- list(
    make_sb_paired_binary_prediction(ids, base_probability, "reference"),
    make_sb_paired_binary_prediction(
      ids,
      rev(base_probability),
      "valid_candidate"
    ),
    make_sb_paired_binary_prediction(
      ids,
      pmin(base_probability + 0.02, 0.99),
      "failed_candidate"
    )
  )
  predictions[[3L]]$predictions$patient_id[[1L]] <- "wrong_patient"
  endpoint <- sb_binary_endpoint(
    data.frame(patient_id = ids, outcome = outcome),
    "patient_id",
    "outcome"
  )
  comparison <- sb_compare_external(predictions, endpoint)

  observed <- sb_bootstrap_external_comparison(
    comparison,
    reference_pipeline_id = "reference",
    reps = 12L,
    confidence_level = 0.90,
    seed = 9601L
  )

  expect_identical(comparison$status, "incomplete")
  expect_identical(observed$status, "incomplete")
  expect_identical(
    observed$candidate_status$candidate_pipeline_id,
    c("valid_candidate", "failed_candidate")
  )
  expect_identical(
    observed$candidate_status$status,
    c("success", "failed")
  )
  expect_identical(
    observed$candidate_status$error_class,
    c(NA_character_, "stormbreaker_failed_external_evaluation")
  )
  expect_identical(
    unique(observed$intervals$candidate_pipeline_id),
    "valid_candidate"
  )
  expect_identical(
    unique(observed$differences$candidate_pipeline_id),
    "valid_candidate"
  )
  expect_identical(
    "failed_candidate" %in% observed$model_replicates$pipeline_id,
    FALSE
  )
  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "No ranking, winner, or p value", fixed = TRUE)
})

# ---- collect-external-comparisons ----

make_sb_binary_external_comparison <- function(
  prefix,
  reverse = FALSE,
  second_signature_id = "published_b",
  fail_second = FALSE
) {
  index <- seq_len(24L)
  data <- data.frame(
    patient_id = paste0(prefix, "_", index),
    outcome = rep(c(0L, 1L), 12L),
    gene_a = seq(-1, 1, length.out = 24L),
    gene_b = cos(index),
    stringsAsFactors = FALSE
  )
  signatures <- list(
    sb_linear_signature(
      "published_a", "binary", c(gene_a = 0.8, gene_b = -0.2),
      intercept = -0.1
    ),
    sb_linear_signature(
      second_signature_id, "binary", c(gene_a = -0.3, gene_b = 0.6),
      intercept = 0.2
    )
  )
  if (reverse) {
    signatures <- rev(signatures)
  }
  predictions <- lapply(signatures, function(signature) {
    sb_predict_signature(
      signature,
      data[c("patient_id", "gene_a", "gene_b")],
      "patient_id"
    )
  })
  if (fail_second) {
    predictions[[2L]]$predictions$patient_id[[1L]] <- "wrong_patient"
  }
  sb_compare_external(
    predictions,
    sb_binary_endpoint(data, "patient_id", "outcome")
  )
}

make_sb_survival_external_comparison <- function(
  prefix,
  time_unit,
  horizon = NULL,
  endpoint_horizon = TRUE,
  fail_prediction = FALSE
) {
  index <- seq_len(24L)
  data <- data.frame(
    patient_id = paste0(prefix, "_", index),
    time = index,
    event = rep(c(1L, 0L, 1L), 8L),
    gene_a = seq(-1, 1, length.out = 24L),
    stringsAsFactors = FALSE
  )
  signature_args <- list(
    signature_id = "published_survival",
    task_type = "survival",
    coefficients = c(gene_a = 0.4),
    intercept = 0
  )
  if (!is.null(horizon)) {
    signature_args$prediction_horizon <- horizon
    signature_args$time_unit <- time_unit
    signature_args$baseline_survival <- 0.8
  }
  signature <- do.call(sb_linear_signature, signature_args)
  predictor_data <- if (fail_prediction) {
    data["patient_id"]
  } else {
    data[c("patient_id", "gene_a")]
  }
  prediction <- sb_predict_signature(
    signature,
    predictor_data,
    "patient_id"
  )
  endpoint_args <- list(
    data = data,
    id_col = "patient_id",
    time_col = "time",
    event_col = "event",
    time_unit = time_unit
  )
  if (!is.null(horizon) && endpoint_horizon) {
    endpoint_args$prediction_horizon <- horizon
  }
  sb_compare_external(
    list(prediction),
    do.call(
      sb_surv_endpoint,
      endpoint_args
    )
  )
}

test_that("survival collections retain the common evaluation target", {
  horizon_collection <- sb_collect_external_comparisons(list(
    cohort_a = make_sb_survival_external_comparison("horizon_a", "month", 6),
    cohort_b = make_sb_survival_external_comparison("horizon_b", "month", 6)
  ))
  risk_collection <- sb_collect_external_comparisons(list(
    cohort_a = make_sb_survival_external_comparison("risk_a", "month"),
    cohort_b = make_sb_survival_external_comparison("risk_b", "month")
  ))
  inherited_collection <- sb_collect_external_comparisons(list(
    cohort_a = make_sb_survival_external_comparison(
      "inherited_a", "month", 6, endpoint_horizon = FALSE
    ),
    cohort_b = make_sb_survival_external_comparison(
      "inherited_b", "month", 6, endpoint_horizon = FALSE
    )
  ))
  failed_inherited <- make_sb_survival_external_comparison(
    "failed_inherited",
    "month",
    6,
    endpoint_horizon = FALSE,
    fail_prediction = TRUE
  )
  mixed_status_collection <- sb_collect_external_comparisons(list(
    successful = inherited_collection$comparisons$cohort_a,
    failed = failed_inherited
  ))

  expect_identical(horizon_collection$status, "success")
  expect_identical(horizon_collection$task_type, "survival")
  expect_identical(horizon_collection$cohorts$time_unit, c("month", "month"))
  expect_identical(horizon_collection$cohorts$prediction_horizon, c(6, 6))
  expect_identical(
    horizon_collection$results$cohort_id,
    c("cohort_a", "cohort_b")
  )
  expect_identical(risk_collection$status, "success")
  expect_identical(
    "prediction_horizon" %in% names(risk_collection$cohorts),
    FALSE
  )
  expect_identical(
    inherited_collection$cohorts$prediction_horizon,
    c(6, 6)
  )
  expect_identical(
    c("ipcw_brier", "time_auc", "uno_c") %in%
      names(inherited_collection$results),
    rep(TRUE, 3L)
  )
  expect_identical(mixed_status_collection$status, "incomplete")
  expect_identical(
    mixed_status_collection$cohorts$prediction_horizon,
    c(6, 6)
  )
})

test_that("external comparisons retain results in cohort order", {
  cohort_a <- make_sb_binary_external_comparison("a")
  cohort_b <- make_sb_binary_external_comparison("b", reverse = TRUE)

  observed <- sb_collect_external_comparisons(list(
    cohort_a = cohort_a,
    cohort_b = cohort_b
  ))

  expect_s3_class(observed, "stormbreaker_external_collection")
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expected_ids <- c("signature__published_a", "signature__published_b")
  expect_identical(observed$pipeline_ids, expected_ids)
  expect_identical(observed$cohorts$cohort_id, c("cohort_a", "cohort_b"))
  expect_identical(observed$cohorts$n_patients, c(24L, 24L))
  expect_identical(
    observed$results$cohort_id,
    rep(c("cohort_a", "cohort_b"), each = 2L)
  )
  expect_identical(
    observed$results$pipeline_id,
    rep(expected_ids, 2L)
  )
  expect_equal(
    observed$results$roc_auc[1:2],
    cohort_a$comparison$roc_auc,
    tolerance = 0
  )
  expect_equal(
    observed$results$roc_auc[3:4],
    rev(cohort_b$comparison$roc_auc),
    tolerance = 0
  )
  expect_identical(observed$comparisons, list(
    cohort_a = cohort_a,
    cohort_b = cohort_b
  ))

  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "No metric is pooled", fixed = TRUE)
  expect_match(printed, "no pipeline is ranked or selected", fixed = TRUE)
})

test_that("incomplete cohort comparisons remain visible", {
  cohort_a <- make_sb_binary_external_comparison("complete")
  cohort_b <- make_sb_binary_external_comparison(
    "incomplete",
    fail_second = TRUE
  )

  observed <- sb_collect_external_comparisons(list(
    complete = cohort_a,
    incomplete = cohort_b
  ))

  expect_identical(observed$status, "incomplete")
  expect_identical(observed$cohorts$status, c("success", "incomplete"))
  expect_identical(
    observed$results$status,
    c("success", "success", "success", "failed")
  )
  expect_identical(
    observed$results$error_class[[4L]],
    "stormbreaker_external_patient_mismatch"
  )
})

test_that("external collection rejects ambiguous or incompatible inputs", {
  binary <- make_sb_binary_external_comparison("binary")
  other_pipelines <- make_sb_binary_external_comparison(
    "other",
    second_signature_id = "different_signature"
  )
  month_6 <- make_sb_survival_external_comparison("month_6", "month", 6)
  month_12 <- make_sb_survival_external_comparison("month_12", "month", 12)
  year_6 <- make_sb_survival_external_comparison("year_6", "year", 6)
  risk_only <- make_sb_survival_external_comparison("risk_only", "month")
  inherited_6 <- make_sb_survival_external_comparison(
    "inherited_6", "month", 6, endpoint_horizon = FALSE
  )
  inherited_12 <- make_sb_survival_external_comparison(
    "inherited_12", "month", 12, endpoint_horizon = FALSE
  )
  failed_inherited_6 <- make_sb_survival_external_comparison(
    "failed_inherited_6",
    "month",
    6,
    endpoint_horizon = FALSE,
    fail_prediction = TRUE
  )
  failed_inherited_12 <- make_sb_survival_external_comparison(
    "failed_inherited_12",
    "month",
    12,
    endpoint_horizon = FALSE,
    fail_prediction = TRUE
  )
  unknown_target <- failed_inherited_6
  unknown_target$prediction_target <- list(
    target_status = "unknown",
    has_horizon = NA,
    prediction_horizon = NULL,
    time_unit = "month",
    horizon_source = "unknown"
  )

  empty <- capture_sb_error(sb_collect_external_comparisons(list()))
  unnamed <- capture_sb_error(sb_collect_external_comparisons(list(binary)))
  direct <- capture_sb_error(sb_collect_external_comparisons(binary))
  duplicate_names <- capture_sb_error(sb_collect_external_comparisons(
    structure(list(binary, binary), names = c("same", "same"))
  ))
  task_mismatch <- capture_sb_error(sb_collect_external_comparisons(list(
    cohort_a = binary,
    cohort_b = month_6
  )))
  pipeline_mismatch <- capture_sb_error(sb_collect_external_comparisons(list(
    cohort_a = binary,
    cohort_b = other_pipelines
  )))
  horizon_mismatch <- capture_sb_error(sb_collect_external_comparisons(list(
    cohort_a = month_6,
    cohort_b = month_12
  )))
  horizon_presence_mismatch <- capture_sb_error(
    sb_collect_external_comparisons(list(
      cohort_a = month_6,
      cohort_b = risk_only
    ))
  )
  unit_mismatch <- capture_sb_error(sb_collect_external_comparisons(list(
    cohort_a = month_6,
    cohort_b = year_6
  )))
  inherited_horizon_mismatch <- capture_sb_error(
    sb_collect_external_comparisons(list(
      cohort_a = inherited_6,
      cohort_b = inherited_12
    ))
  )
  failed_horizon_mismatch <- capture_sb_error(
    sb_collect_external_comparisons(list(
      cohort_a = failed_inherited_6,
      cohort_b = failed_inherited_12
    ))
  )
  unknown_target_error <- capture_sb_error(
    sb_collect_external_comparisons(list(cohort_a = unknown_target))
  )

  expect_s3_class(empty, "stormbreaker_bad_external_collection")
  expect_s3_class(unnamed, "stormbreaker_bad_external_collection")
  expect_s3_class(direct, "stormbreaker_bad_external_collection")
  expect_s3_class(duplicate_names, "stormbreaker_bad_external_collection")
  expect_s3_class(task_mismatch, "stormbreaker_external_task_mismatch")
  expect_s3_class(
    pipeline_mismatch,
    "stormbreaker_external_pipeline_mismatch"
  )
  expect_s3_class(
    horizon_mismatch,
    "stormbreaker_external_target_mismatch"
  )
  expect_s3_class(
    horizon_presence_mismatch,
    "stormbreaker_external_target_mismatch"
  )
  expect_s3_class(unit_mismatch, "stormbreaker_external_target_mismatch")
  expect_s3_class(
    inherited_horizon_mismatch,
    "stormbreaker_external_target_mismatch"
  )
  expect_s3_class(
    failed_horizon_mismatch,
    "stormbreaker_external_target_mismatch"
  )
  expect_s3_class(
    unknown_target_error,
    "stormbreaker_external_target_mismatch"
  )
})

# ---- compare-external ----

test_that("binary locked predictions share one external endpoint", {
  development <- make_sb_binary_fixture(n = 60L, seed = 8101L)
  development_task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )
  fits <- lapply(pipelines, function(pipeline) {
    sb_fit_final(development_task, pipeline)
  })

  external <- make_sb_binary_fixture(n = 36L, seed = 8102L)
  external$patient_id <- paste0("external_", external$patient_id)
  predictors <- external[c(
    "patient_id", "gene_a", "gene_b", "gene_c"
  )]
  predictions <- lapply(fits, function(fit) {
    sb_predict_final(fit, predictors, "patient_id")
  })
  endpoint <- sb_binary_endpoint(
    external,
    "patient_id",
    "outcome"
  )

  observed <- sb_compare_external(predictions, endpoint)
  direct <- lapply(predictions, sb_evaluate_external, endpoint = endpoint, reverse = FALSE)
  pipeline_ids <- vapply(predictions, `[[`, character(1L), "pipeline_id")
  names(direct) <- pipeline_ids

  expect_s3_class(observed, "stormbreaker_external_comparison")
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$n_patients, 36L)
  expect_identical(observed$evaluations, direct)
  expect_identical(observed$comparison$pipeline_id, pipeline_ids)
  expect_equal(
    observed$comparison$roc_auc,
    unname(vapply(
      direct, function(x) x$metrics$roc_auc$estimate, numeric(1L)
    )),
    tolerance = 0
  )
  expect_equal(
    observed$comparison$brier_score,
    unname(vapply(
      direct, function(x) x$metrics$brier_score$estimate, numeric(1L)
    )),
    tolerance = 0
  )
  expect_equal(
    observed$comparison$mean_calibration_ratio,
    unname(vapply(
      direct,
      function(x) x$metrics$mean_calibration$mean_calibration_ratio,
      numeric(1L)
    )),
    tolerance = 0
  )
  expect_identical(
    unname(vapply(
      observed$evaluations,
      function(x) x$patient_match$output_order,
      character(1L)
    )),
    rep("external_endpoint", 2L)
  )
  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "No pipeline is ranked or selected", fixed = TRUE)
})

test_that("survival horizon comparisons copy existing external metrics", {
  development <- make_sb_fixture(n = 72L, seed = 8201L)
  development_task <- sb_surv_task(
    development,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 6
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    )
  )
  fits <- lapply(pipelines, function(pipeline) {
    sb_fit_final(development_task, pipeline)
  })

  external <- make_sb_fixture(n = 60L, seed = 8202L)
  external$patient_id <- paste0("external_", external$patient_id)
  predictors <- external[c(
    "patient_id", "gene_a", "gene_b", "gene_c"
  )]
  predictions <- lapply(fits, function(fit) {
    sb_predict_final(fit, predictors, "patient_id")
  })
  endpoint <- sb_surv_endpoint(
    external,
    "patient_id",
    "time",
    "event",
    "month",
    prediction_horizon = 6
  )

  observed <- sb_compare_external(predictions, endpoint)
  direct <- lapply(predictions, sb_evaluate_external, endpoint = endpoint, reverse = FALSE)
  pipeline_ids <- vapply(predictions, `[[`, character(1L), "pipeline_id")
  names(direct) <- pipeline_ids

  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "survival")
  expect_identical(observed$endpoint$prediction_horizon, 6)
  expect_identical(observed$evaluations, direct)
  expect_equal(
    observed$comparison$harrell_c,
    unname(vapply(
      direct, function(x) x$metrics$harrell_c$estimate, numeric(1L)
    )),
    tolerance = 0
  )
  expect_equal(
    observed$comparison$ipcw_brier,
    unname(vapply(
      direct, function(x) x$metrics$ipcw_brier$estimate, numeric(1L)
    )),
    tolerance = 0
  )
  expect_equal(
    observed$comparison$time_auc,
    unname(vapply(
      direct, function(x) x$metrics$time_auc$estimate, numeric(1L)
    )),
    tolerance = 0
  )
  expect_equal(
    observed$comparison$uno_c,
    unname(vapply(
      direct, function(x) x$metrics$uno_c$estimate, numeric(1L)
    )),
    tolerance = 0
  )
})

test_that("external comparison retains a horizon inherited from predictions", {
  index <- seq_len(30L)
  data <- data.frame(
    patient_id = paste0("inherited_", index),
    time = index,
    event = rep(c(1L, 0L, 1L), 10L),
    gene_a = seq(-1, 1, length.out = 30L),
    stringsAsFactors = FALSE
  )
  signatures <- list(
    sb_linear_signature(
      "inherited_a", "survival", c(gene_a = 0.3), intercept = 0,
      prediction_horizon = 6, time_unit = "month",
      baseline_survival = 0.8
    ),
    sb_linear_signature(
      "inherited_b", "survival", c(gene_a = -0.2), intercept = 0,
      prediction_horizon = 6, time_unit = "month",
      baseline_survival = 0.75
    )
  )
  predictions <- lapply(signatures, function(signature) {
    sb_predict_signature(
      signature,
      data[c("patient_id", "gene_a")],
      "patient_id"
    )
  })
  endpoint <- sb_surv_endpoint(
    data, "patient_id", "time", "event", "month"
  )

  observed <- sb_compare_external(predictions, endpoint)
  direct <- lapply(predictions, sb_evaluate_external, endpoint = endpoint, reverse = FALSE)

  expect_identical(observed$prediction_target$has_horizon, TRUE)
  expect_identical(observed$prediction_target$prediction_horizon, 6)
  expect_identical(observed$prediction_target$time_unit, "month")
  expect_identical(observed$prediction_target$horizon_source, "prediction")
  expect_identical("prediction_horizon" %in% names(observed$endpoint), FALSE)
  expect_identical(
    c("ipcw_brier", "time_auc", "uno_c") %in% names(observed$comparison),
    rep(TRUE, 3L)
  )
  expect_equal(
    observed$comparison$ipcw_brier,
    vapply(direct, function(x) x$metrics$ipcw_brier$estimate, numeric(1L)),
    tolerance = 0
  )
})

test_that("comparison rejects unlike inherited survival targets", {
  index <- seq_len(24L)
  data <- data.frame(
    patient_id = paste0("mixed_", index),
    time = index,
    event = rep(c(1L, 0L, 1L), 8L),
    gene_a = seq(-1, 1, length.out = 24L),
    stringsAsFactors = FALSE
  )
  make_prediction <- function(signature_id, horizon = NULL) {
    arguments <- list(
      signature_id = signature_id,
      task_type = "survival",
      coefficients = c(gene_a = 0.3),
      intercept = 0
    )
    if (!is.null(horizon)) {
      arguments$prediction_horizon <- horizon
      arguments$time_unit <- "month"
      arguments$baseline_survival <- 0.8
    }
    sb_predict_signature(
      do.call(sb_linear_signature, arguments),
      data[c("patient_id", "gene_a")],
      "patient_id"
    )
  }
  endpoint <- sb_surv_endpoint(
    data, "patient_id", "time", "event", "month"
  )
  different_horizons <- capture_sb_error(sb_compare_external(
    list(make_prediction("h6", 6), make_prediction("h12", 12)),
    endpoint
  ))
  mixed_presence <- capture_sb_error(sb_compare_external(
    list(make_prediction("risk"), make_prediction("horizon", 6)),
    endpoint
  ))

  expect_s3_class(
    different_horizons,
    "stormbreaker_external_target_mismatch"
  )
  expect_s3_class(mixed_presence, "stormbreaker_external_target_mismatch")
})

test_that("failed predictions retain recoverable survival targets", {
  index <- seq_len(20L)
  data <- data.frame(
    patient_id = paste0("failed_", index),
    time = index,
    event = rep(c(1L, 0L), 10L),
    gene_a = seq(-1, 1, length.out = 20L),
    stringsAsFactors = FALSE
  )
  signature <- sb_linear_signature(
    "failed_horizon", "survival", c(gene_a = 0.4), intercept = 0,
    prediction_horizon = 6, time_unit = "month",
    baseline_survival = 0.8
  )
  failed <- sb_predict_signature(
    signature,
    data["patient_id"],
    "patient_id"
  )
  endpoint <- sb_surv_endpoint(
    data, "patient_id", "time", "event", "month"
  )

  observed <- sb_compare_external(list(failed), endpoint)

  expect_identical(failed$status, "failed")
  expect_identical(observed$status, "incomplete")
  expect_identical(observed$prediction_target$target_status, "complete")
  expect_identical(observed$prediction_target$has_horizon, TRUE)
  expect_identical(observed$prediction_target$prediction_horizon, 6)
  expect_identical(observed$prediction_target$time_unit, "month")

  unknown <- failed
  unknown$signature <- NULL
  unknown$predictions <- unknown$predictions[0, , drop = FALSE]
  unknown_comparison <- sb_compare_external(list(unknown), endpoint)
  expect_identical(
    unknown_comparison$prediction_target$target_status,
    "unknown"
  )
  expect_identical(
    is.na(unknown_comparison$prediction_target$has_horizon),
    TRUE
  )
})

test_that("a failed prediction with unknown task makes the target partial", {
  index <- seq_len(20L)
  data <- data.frame(
    patient_id = paste0("partial_", index),
    time = index,
    event = rep(c(1L, 0L), 10L),
    gene_a = seq(-1, 1, length.out = 20L),
    stringsAsFactors = FALSE
  )
  known_signature <- sb_linear_signature(
    "known_horizon", "survival", c(gene_a = 0.4), intercept = 0,
    prediction_horizon = 6, time_unit = "month",
    baseline_survival = 0.8
  )
  known <- sb_predict_signature(
    known_signature,
    data[c("patient_id", "gene_a")],
    "patient_id"
  )
  unknown <- sb_predict_signature(
    list(
      signature_id = "unknown_task",
      task_type = "not_a_task",
      coefficients = c(gene_a = 0.4),
      intercept = 0
    ),
    data[c("patient_id", "gene_a")],
    "patient_id"
  )
  endpoint <- sb_surv_endpoint(
    data, "patient_id", "time", "event", "month"
  )

  observed <- sb_compare_external(list(known, unknown), endpoint)
  collection_error <- capture_sb_error(
    sb_collect_external_comparisons(list(cohort_a = observed))
  )

  expect_identical(unknown$status, "failed")
  expect_identical(is.na(unknown$task_type), TRUE)
  expect_identical(observed$status, "incomplete")
  expect_identical(observed$prediction_target$target_status, "partial")
  expect_s3_class(
    collection_error,
    "stormbreaker_external_target_mismatch"
  )
})

test_that("one invalid external prediction remains visible without fallback", {
  development <- make_sb_binary_fixture(n = 60L, seed = 8301L)
  development_task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )
  fits <- lapply(pipelines, function(pipeline) {
    sb_fit_final(development_task, pipeline)
  })

  external <- make_sb_binary_fixture(n = 30L, seed = 8302L)
  external$patient_id <- paste0("external_", external$patient_id)
  predictors <- external[c(
    "patient_id", "gene_a", "gene_b", "gene_c"
  )]
  predictions <- lapply(fits, function(fit) {
    sb_predict_final(fit, predictors, "patient_id")
  })
  predictions[[2L]]$predictions$patient_id[[1L]] <- "wrong_patient"
  endpoint <- sb_binary_endpoint(
    external,
    "patient_id",
    "outcome"
  )

  observed <- sb_compare_external(predictions, endpoint)

  expect_identical(observed$status, "incomplete")
  expect_identical(observed$comparison$status, c("success", "failed"))
  expect_identical(
    observed$comparison$error_stage,
    c(NA_character_, "patient_match")
  )
  expect_identical(
    observed$comparison$error_class,
    c(NA_character_, "stormbreaker_external_patient_mismatch")
  )
  expect_identical(is.na(observed$comparison$roc_auc), c(FALSE, TRUE))
  expect_identical(observed$evaluations[[1L]]$status, "success")
  expect_identical(observed$evaluations[[2L]]$status, "failed")
  expect_identical(
    observed$evaluations[[2L]]$error_stage,
    "patient_match"
  )
})

test_that("external comparison requires distinct final predictions", {
  development <- make_sb_binary_fixture(n = 30L, seed = 8401L)
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  prediction <- sb_predict_final(
    fit,
    development[c("patient_id", "gene_a", "gene_b", "gene_c")],
    "patient_id"
  )
  endpoint <- sb_binary_endpoint(development, "patient_id", "outcome")

  empty <- capture_sb_error(sb_compare_external(list(), endpoint))
  duplicate <- capture_sb_error(sb_compare_external(
    list(prediction, prediction),
    endpoint
  ))

  expect_s3_class(empty, "stormbreaker_bad_prediction_collection")
  expect_s3_class(duplicate, "stormbreaker_duplicate_pipeline_id")
})
})
