# ---- final-fit ----
local({
final_fit_contains_function <- function(value) {
  if (is.function(value) || is.environment(value) ||
    typeof(value) == "externalptr") {
    return(TRUE)
  }
  if (!is.list(value)) {
    return(FALSE)
  }
  any(vapply(unclass(value), final_fit_contains_function, logical(1)))
}

final_state_oracle <- function(state, x) {
  switch(state$learner,
    coxph = sb_predict_cox_learner_state(state, x),
    glmnet_cox = sb_predict_cox_learner_state(state, x),
    coxboost = sb_predict_cox_learner_state(state, x),
    plsrcox = sb_predict_cox_learner_state(state, x),
    superpc = sb_predict_cox_learner_state(state, x),
    survivalsvm = sb_predict_cox_learner_state(state, x),
    logistic = sb_predict_logistic_learner_state(state, x),
    glmnet_logistic = sb_predict_logistic_learner_state(state, x),
    gaussian_naive_bayes = sb_predict_naive_bayes_learner_state(state, x),
    lightgbm_binary = sb_predict_lightgbm_learner_state(state, x),
    catboost_binary = sb_predict_catboost_learner_state(state, x),
    catboost_cox = sb_predict_catboost_cox_learner_state(state, x),
    gbm_cox = sb_predict_gbm_cox_learner_state(state, x),
    ranger_binary = sb_predict_ranger_binary_learner_state(state, x),
    adaboost_binary =
      sb_predict_classical_binary_boosting_learner_state(state, x),
    logitboost_binary =
      sb_predict_classical_binary_boosting_learner_state(state, x),
    ranger_survival = sb_predict_ranger_survival_learner_state(state, x)
  )
}

expect_final_prediction_failure <- function(result, error_class, n_rows) {
  expect_s3_class(result, "stormbreaker_final_prediction")
  expect_identical(result$status, "failed")
  expect_identical(result$error_class, error_class)
  expect_identical(nrow(result$predictions), as.integer(n_rows))
  if (n_rows > 0L) {
    expect_true(all(result$predictions$status == "failed"))
    expect_true(all(result$predictions$error_class == error_class))
  }
}

test_that("ordinary Cox final fit uses every development patient", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_b", "gene_a"),
    "month"
  )
  pipeline <- sb_pipeline()
  fit <- sb_fit_final(task, pipeline)

  expect_s3_class(fit, "stormbreaker_final_fit")
  expect_identical(fit$status, "success")
  expect_identical(fit$task_type, "survival")
  expect_identical(fit$pipeline, pipeline)
  expect_identical(
    fit$development$patient_ids,
    as.character(d$patient_id)
  )
  expect_identical(fit$development$n, 60L)
  expect_identical(fit$development$independent_unit, "patient")
  expect_identical(
    fit$development$candidate_features,
    c("gene_b", "gene_a")
  )
  expect_identical(fit$development$endpoint$n_events, task$n_events)
  expect_identical(fit$development$endpoint$n_censored, task$n_censored)
  expect_identical(fit$selection$status, "success")
  expect_identical(
    fit$selection$selected_features,
    c("gene_b", "gene_a")
  )
  expect_s3_class(fit$selection$state, "stormbreaker_selector_fit")
  expect_identical(fit$tuning$status, "not_required")
  expect_identical(fit$learner$status, "success")
  expect_s3_class(fit$learner_state, "stormbreaker_learner_state")
  expect_identical(fit$required_features, c("gene_b", "gene_a"))
  expect_identical(fit$prediction$risk_formula, "x_beta")
  expect_identical(fit$prediction$risk_reference, "zero")
  expect_true(all(c(
    "package", "package_version", "r_version", "survival_version"
  ) %in% names(fit$provenance)))
  expect_false(any(c(
    "metrics", "folds", "comparison", "selection_stability",
    "predictions", "runs", "task", "data"
  ) %in% names(fit)))
  expect_false("data" %in% names(fit$development))
  expect_false(final_fit_contains_function(fit))

  printed <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(printed, "final development fit", fixed = TRUE)
  expect_match(printed, pipeline$id, fixed = TRUE)
  expect_match(printed, "survival", fixed = TRUE)
  expect_match(printed, "Required features: 2", fixed = TRUE)
  expect_false(grepl("AUC|Brier|concordance|validation performance", printed))
})

test_that("Top-K Cox final selection is fitted on all development rows", {
  d <- make_sb_fixture(n = 70L, seed = 20260820L)
  predictors <- c("gene_c", "gene_a", "gene_b")
  task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month"
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "coxph",
    selector_n_features = 2L
  )
  fit <- sb_fit_final(task, pipeline)
  direct <- sb_fit_univ_cox_top_k(
    d[predictors],
    d$time,
    as.integer(d$event),
    n_features = 2L
  )

  expect_identical(fit$status, "success")
  expect_identical(
    fit$selection$selected_features,
    direct$selected_features
  )
  expect_identical(fit$selection$scores, direct$scores)
  expect_identical(fit$selection$state, direct)
  expect_identical(fit$required_features, direct$selected_features)

  new_data <- d[
    1:4,
    c("patient_id", fit$required_features),
    drop = FALSE
  ]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  expect_identical(prediction$status, "success")
  expect_identical(nrow(prediction$predictions), 4L)
})

test_that("glmnet Cox final tuning covers the complete development cohort", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60L, seed = 20260821L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month"
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 1701L
  )
  fit <- sb_fit_final(task, pipeline)
  direct <- sb_prepare_runtime_learner(
    sb_resolve_pipeline_runtime(pipeline)$learner,
    pipeline,
    d[predictors],
    list(
      task_type = "survival",
      time = d$time,
      event = as.integer(d$event)
    ),
    d$patient_id
  )$tuning

  expect_identical(fit$status, "success")
  expect_identical(fit$tuning$status, "success")
  expect_identical(fit$tuning$lambda, direct$lambda)
  expect_identical(fit$tuning$curve, direct$curve)
  expect_identical(fit$tuning$fold_id, direct$fold_id)
  expect_identical(names(fit$tuning$fold_id), as.character(d$patient_id))
  expect_identical(length(fit$tuning$fold_id), nrow(d))
  expect_true(fit$tuning$lambda %in% fit$tuning$curve$lambda)
  expect_identical(
    fit$provenance$glmnet_version,
    as.character(utils::packageVersion("glmnet"))
  )
})

test_that("core binary learners produce portable final states", {
  d <- make_sb_binary_fixture(n = 80L, seed = 20260822L)
  predictors <- c("gene_b", "gene_a", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(
      learner = "gaussian_naive_bayes",
      naive_bayes_var_smoothing = 1e-8
    )
  )

  for (pipeline in pipelines) {
    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success", info = pipeline$id)
    expect_identical(fit$task_type, "binary", info = pipeline$id)
    expect_identical(fit$pipeline, pipeline, info = pipeline$id)
    expect_identical(fit$selection$status, "success", info = pipeline$id)
    expect_identical(fit$tuning$status, "not_required", info = pipeline$id)
    expect_identical(fit$learner$status, "success", info = pipeline$id)
    expect_s3_class(
      fit$learner_state,
      "stormbreaker_learner_state"
    )
    expect_identical(fit$required_features, predictors, info = pipeline$id)
    expect_identical(fit$prediction$positive_class, 1L, info = pipeline$id)
    expect_false(final_fit_contains_function(fit), info = pipeline$id)
  }
})

test_that("glmnet logistic final tuning uses all development patients", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 80L, seed = 20260823L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0.25,
    inner_folds = 4L,
    inner_seed = 1702L,
    glmnet_standardize = FALSE,
    glmnet_lambda_rule = "lambda.1se",
    glmnet_tuning_measure = "auc"
  )
  fit <- sb_fit_final(task, pipeline)
  direct <- sb_prepare_runtime_learner(
    sb_resolve_pipeline_runtime(pipeline)$learner,
    pipeline,
    d[predictors],
    list(task_type = "binary", outcome = as.integer(d$outcome)),
    d$patient_id
  )$tuning

  expect_identical(fit$status, "success")
  expect_identical(fit$tuning$status, "success")
  expect_identical(fit$tuning$lambda, direct$lambda)
  expect_identical(fit$tuning$curve, direct$curve)
  expect_identical(fit$tuning$fold_id, direct$fold_id)
  expect_identical(fit$tuning$standardize, FALSE)
  expect_identical(fit$tuning$selection_rule, "lambda.1se")
  expect_identical(fit$tuning$measure, "auc")
  expect_identical(fit$learner_state$params$standardize, FALSE)
  expect_identical(names(fit$tuning$fold_id), as.character(d$patient_id))
  expect_identical(length(fit$tuning$fold_id), nrow(d))
  expect_true(fit$tuning$lambda %in% fit$tuning$curve$lambda)
})

test_that("fixed tree learners produce portable final states", {
  d <- make_sb_binary_fixture(n = 80L, seed = 20260824L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)

  if (sb_lightgbm_available()) {
    pipeline <- sb_pipeline(
      learner = "lightgbm_binary",
      lightgbm_nrounds = 5L,
      lightgbm_num_leaves = 7L,
      lightgbm_min_data_in_leaf = 3L,
      lightgbm_seed = 101L
    )
    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    expect_type(fit$learner_state$model_string, "character")
    expect_false("booster" %in% names(fit$learner_state))
    expect_false(final_fit_contains_function(fit))
    expect_identical(
      fit$provenance$lightgbm_version,
      as.character(utils::packageVersion("lightgbm"))
    )
  }

  if (sb_catboost_available()) {
    pipeline <- sb_pipeline(
      learner = "catboost_binary",
      catboost_iterations = 5L,
      catboost_depth = 3L,
      catboost_seed = 101L
    )
    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    expect_type(fit$learner_state$model_raw, "raw")
    expect_false("booster" %in% names(fit$learner_state))
    expect_false(final_fit_contains_function(fit))
    expect_identical(
      fit$provenance$catboost_version,
      as.character(utils::packageVersion("catboost"))
    )
  }
})

test_that("final fitting requires one explicit task-matched pipeline", {
  d <- make_sb_fixture(n = 40L, seed = 20260825L)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month"
  )

  expect_s3_class(
    capture_sb_error(sb_fit_final(task, NULL)),
    "stormbreaker_bad_pipeline"
  )
  expect_s3_class(
    capture_sb_error(sb_fit_final(task, list(sb_pipeline()))),
    "stormbreaker_bad_pipeline"
  )
  expect_s3_class(
    capture_sb_error(sb_fit_final(
      task,
      structure(list(), class = "stormbreaker_cv_result")
    )),
    "stormbreaker_bad_pipeline"
  )
  expect_s3_class(
    capture_sb_error(sb_fit_final(
      task,
      sb_pipeline(learner = "logistic")
    )),
    "stormbreaker_pipeline_task_mismatch"
  )
})

test_that("final-fit failures retain completed earlier stages", {
  selector_data <- make_sb_fixture(n = 30L, seed = 20260826L)
  selector_data$gene_a <- 1
  selector_data$gene_b <- 1
  selector_task <- sb_surv_task(
    selector_data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month"
  )
  selector_fit <- sb_fit_final(
    selector_task,
    sb_pipeline(
      selector = "univ_cox_top_k",
      learner = "coxph",
      selector_n_features = 1L
    )
  )
  expect_identical(selector_fit$status, "failed")
  expect_identical(selector_fit$error_stage, "selection")
  expect_identical(selector_fit$selection$status, "failed")
  expect_identical(selector_fit$learner$status, "not_fitted")
  expect_null(selector_fit$learner_state)

  tuning_data <- make_sb_binary_fixture(n = 10L, seed = 20260827L)
  tuning_task <- sb_binary_task(
    tuning_data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  tuning_fit <- sb_fit_final(
    tuning_task,
    sb_pipeline(
      learner = "glmnet_logistic",
      inner_folds = 6L,
      inner_seed = 19L
    )
  )
  expect_identical(tuning_fit$status, "failed")
  expect_identical(tuning_fit$error_stage, "tuning")
  expect_identical(tuning_fit$selection$status, "success")
  expect_identical(tuning_fit$tuning$status, "failed")
  expect_identical(tuning_fit$learner$status, "not_fitted")
  expect_null(tuning_fit$learner_state)

  missing_pipeline <- sb_pipeline(
    learner = "lightgbm_binary",
    lightgbm_nrounds = 3L,
    lightgbm_num_leaves = 5L,
    lightgbm_min_data_in_leaf = 2L
  )
  testthat::local_mocked_bindings(
    sb_lightgbm_available = function() FALSE
  )
  missing_fit <- sb_fit_final(tuning_task, missing_pipeline)
  expect_identical(missing_fit$status, "failed")
  expect_identical(missing_fit$error_stage, "learner_fit")
  expect_identical(missing_fit$selection$status, "success")
  expect_identical(missing_fit$learner$status, "failed")
  expect_identical(
    missing_fit$error_class,
    "stormbreaker_missing_lightgbm"
  )
})

test_that("learner and state failures remain distinct", {
  d <- make_sb_fixture(n = 12L, seed = 20260828L)
  withr::with_seed(20260828L, {
    for (index in seq_len(12L)) {
      d[[sprintf("x_%02d", index)]] <- stats::rnorm(nrow(d))
    }
  })
  predictors <- sprintf("x_%02d", seq_len(12L))
  task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month"
  )
  learner_fit <- sb_fit_final(task, sb_pipeline())

  expect_identical(learner_fit$status, "failed")
  expect_identical(learner_fit$error_stage, "learner_fit")
  expect_identical(learner_fit$selection$status, "success")
  expect_identical(learner_fit$learner$status, "failed")
  expect_null(learner_fit$learner_state)

  ordinary_task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month"
  )
  pipeline <- sb_pipeline()
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  runtime$learner$extract_state <- function(model, context = list()) {
    sb_abort(
      "Synthetic state extraction failure.",
      "stormbreaker_test_state_extraction_failure"
    )
  }
  testthat::local_mocked_bindings(
    sb_resolve_pipeline_runtime = function(pipeline) runtime
  )
  state_fit <- sb_fit_final(ordinary_task, pipeline)

  expect_identical(state_fit$status, "failed")
  expect_identical(state_fit$error_stage, "state_extraction")
  expect_identical(state_fit$selection$status, "success")
  expect_identical(state_fit$learner$status, "success")
  expect_identical(state_fit$required_features, c("gene_a", "gene_b"))
  expect_null(state_fit$learner_state)
  expect_identical(
    state_fit$error_class,
    "stormbreaker_test_state_extraction_failure"
  )
})

test_that("final prediction has an endpoint-blind three-argument interface", {
  expect_identical(
    names(formals(sb_predict_final)),
    c("fit", "new_data", "id_col")
  )
})

test_that("survival final prediction uses only saved state", {
  d <- make_sb_fixture(n = 60L, seed = 20260829L)
  predictors <- c("gene_b", "gene_a")
  rows <- c(9L, 2L, 14L, 5L)

  risk_task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month"
  )
  pipelines <- list(sb_pipeline())
  if (sb_glmnet_available()) {
    pipelines[[length(pipelines) + 1L]] <- sb_pipeline(
      learner = "glmnet_cox",
      glmnet_alpha = 0.5,
      inner_folds = 3L,
      inner_seed = 2026L
    )
  }

  for (pipeline in pipelines) {
    fit <- sb_fit_final(risk_task, pipeline)
    new_data <- d[rows, c("patient_id", rev(predictors)), drop = FALSE]
    prediction <- sb_predict_final(fit, new_data, "patient_id")
    oracle <- final_state_oracle(
      fit$learner_state,
      new_data[fit$required_features]
    )

    expect_identical(prediction$status, "success", info = pipeline$id)
    expect_identical(
      prediction$predictions$patient_id,
      as.character(d$patient_id[rows]),
      info = pipeline$id
    )
    expect_true(all(
      prediction$predictions$pipeline_id == pipeline$id
    ), info = pipeline$id)
    expect_equal(
      prediction$predictions$risk,
      oracle$risk,
      tolerance = 0,
      info = pipeline$id
    )
    expect_false(any(c("time", "event", "outcome") %in%
      names(prediction$predictions)))

    restored <- unserialize(serialize(fit, NULL))
    restored_prediction <- sb_predict_final(
      restored,
      new_data,
      "patient_id"
    )
    expect_identical(restored_prediction, prediction, info = pipeline$id)
  }

  horizon_task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 12
  )
  horizon_fit <- sb_fit_final(horizon_task, sb_pipeline())
  horizon_data <- d[rows, c("patient_id", predictors), drop = FALSE]
  horizon_prediction <- sb_predict_final(
    horizon_fit,
    horizon_data,
    "patient_id"
  )
  risk <- unname(as.vector(
    as.matrix(horizon_data[predictors]) %*%
      horizon_fit$learner_state$coefficient
  ))
  probability <- exp(
    -horizon_fit$learner_state$baseline_hazard$
      cumulative_baseline_hazard * exp(risk)
  )

  expect_identical(horizon_prediction$status, "success")
  expect_equal(horizon_prediction$predictions$risk, risk, tolerance = 1e-12)
  expect_equal(
    horizon_prediction$predictions$survival_probability,
    probability,
    tolerance = 1e-12
  )
  expect_true(all(
    horizon_prediction$predictions$prediction_horizon == 12
  ))
  expect_true(all(horizon_prediction$predictions$time_unit == "month"))
})

test_that("declared survival outputs remain bound to their saved state", {
  d <- make_sb_fixture(n = 60L, seed = 20260907L)
  predictors <- c("gene_a", "gene_b")
  horizon_task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 12
  )
  horizon_fit <- sb_fit_final(horizon_task, sb_pipeline())
  new_data <- d[1:4, c("patient_id", predictors)]

  missing_baseline <- horizon_fit
  missing_baseline$learner_state$baseline_hazard <- NULL
  expect_final_prediction_failure(
    sb_predict_final(missing_baseline, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    4L
  )

  missing_unit <- horizon_fit
  missing_unit$learner_state$baseline_hazard$time_unit <- NULL
  expect_final_prediction_failure(
    sb_predict_final(missing_unit, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    4L
  )

  wrong_horizon <- horizon_fit
  wrong_horizon$learner_state$baseline_hazard$prediction_horizon <- 6
  expect_final_prediction_failure(
    sb_predict_final(wrong_horizon, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    4L
  )

  for (bad_hazard in list(-1, NA_real_, Inf)) {
    invalid_hazard <- horizon_fit
    invalid_hazard$learner_state$baseline_hazard$
      cumulative_baseline_hazard <- bad_hazard
    expect_final_prediction_failure(
      sb_predict_final(invalid_hazard, new_data, "patient_id"),
      "stormbreaker_final_state_semantics_mismatch",
      4L
    )
  }

  risk_task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month"
  )
  risk_fit <- sb_fit_final(risk_task, sb_pipeline())
  risk_fit$learner_state$baseline_hazard <-
    horizon_fit$learner_state$baseline_hazard
  expect_final_prediction_failure(
    sb_predict_final(risk_fit, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    4L
  )
})

test_that("every binary final learner predicts class-one probability", {
  d <- make_sb_binary_fixture(n = 80L, seed = 20260830L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(
      learner = "gaussian_naive_bayes",
      naive_bayes_var_smoothing = 1e-8
    )
  )
  if (sb_glmnet_available()) {
    pipelines[[length(pipelines) + 1L]] <- sb_pipeline(
      learner = "glmnet_logistic",
      glmnet_alpha = 0.5,
      inner_folds = 4L,
      inner_seed = 303L
    )
  }
  if (sb_lightgbm_available()) {
    pipelines[[length(pipelines) + 1L]] <- sb_pipeline(
      learner = "lightgbm_binary",
      lightgbm_nrounds = 5L,
      lightgbm_num_leaves = 7L,
      lightgbm_min_data_in_leaf = 3L,
      lightgbm_seed = 304L
    )
  }
  if (sb_catboost_available()) {
    pipelines[[length(pipelines) + 1L]] <- sb_pipeline(
      learner = "catboost_binary",
      catboost_iterations = 5L,
      catboost_depth = 3L,
      catboost_seed = 305L
    )
  }

  rows <- c(11L, 3L, 8L, 1L, 17L)
  new_data <- d[rows, c("patient_id", rev(predictors)), drop = FALSE]
  for (pipeline in pipelines) {
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(fit, new_data, "patient_id")
    oracle <- final_state_oracle(
      fit$learner_state,
      new_data[fit$required_features]
    )

    expect_identical(prediction$status, "success", info = pipeline$id)
    expect_identical(
      prediction$predictions$patient_id,
      as.character(d$patient_id[rows]),
      info = pipeline$id
    )
    expect_true(all(
      prediction$predictions$pipeline_id == pipeline$id
    ), info = pipeline$id)
    expect_equal(
      prediction$predictions$probability,
      oracle$probability,
      tolerance = 0,
      info = pipeline$id
    )
    expect_true(all(prediction$predictions$status == "success"))
    expect_false(any(c("outcome", "time", "event") %in%
      names(prediction$predictions)))

    restored <- unserialize(serialize(fit, NULL))
    restored_prediction <- sb_predict_final(
      restored,
      new_data,
      "patient_id"
    )
    expect_identical(restored_prediction, prediction, info = pipeline$id)
  }
})

test_that("binary decision thresholds propagate without changing probabilities", {
  d <- make_sb_binary_fixture(n = 80L, seed = 20260908L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  thresholds <- c(0.10, 0.25, 0.50)
  plain_task <- sb_binary_task(
    d, "patient_id", "outcome", predictors
  )
  threshold_task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    predictors,
    decision_thresholds = thresholds
  )
  pipeline <- sb_pipeline(learner = "logistic")
  plain_fit <- sb_fit_final(plain_task, pipeline)
  threshold_fit <- sb_fit_final(threshold_task, pipeline)
  new_data <- d[1:8, c("patient_id", predictors), drop = FALSE]
  plain_prediction <- sb_predict_final(
    plain_fit, new_data, "patient_id"
  )
  threshold_prediction <- sb_predict_final(
    threshold_fit, new_data, "patient_id"
  )

  expect_identical(
    threshold_fit$development$endpoint$decision_thresholds,
    thresholds
  )
  expect_identical(
    threshold_fit$prediction$decision_thresholds,
    thresholds
  )
  expect_identical(
    threshold_prediction$decision_thresholds,
    thresholds
  )
  expect_equal(
    threshold_prediction$predictions$probability,
    plain_prediction$predictions$probability,
    tolerance = 0
  )
  expect_false("decision_thresholds" %in% names(plain_prediction))
  expect_false(any(c(
    "predicted_class", "decision", "action", "treat", "selected"
  ) %in% names(threshold_prediction$predictions)))
})

test_that("binary final-fit threshold copies must remain valid and consistent", {
  d <- make_sb_binary_fixture(n = 60L, seed = 20260909L)
  predictors <- c("gene_a", "gene_b")
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    predictors,
    decision_thresholds = c(0.2, 0.5)
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  new_data <- d[1:5, c("patient_id", predictors), drop = FALSE]

  missing_copy <- fit
  missing_copy$development$endpoint$decision_thresholds <- NULL
  expect_final_prediction_failure(
    sb_predict_final(missing_copy, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    5L
  )

  inconsistent_copy <- fit
  inconsistent_copy$development$endpoint$decision_thresholds <- c(0.2, 0.6)
  expect_final_prediction_failure(
    sb_predict_final(inconsistent_copy, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    5L
  )

  invalid_copy <- fit
  invalid_copy$prediction$decision_thresholds <- c(0.2, 1)
  expect_final_prediction_failure(
    sb_predict_final(invalid_copy, new_data, "patient_id"),
    "stormbreaker_bad_decision_thresholds",
    5L
  )
})

test_that("irrelevant endpoint columns are never inspected", {
  binary <- make_sb_binary_fixture(n = 60L, seed = 20260831L)
  binary_task <- sb_binary_task(
    binary,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  binary_fit <- sb_fit_final(
    binary_task,
    sb_pipeline(learner = "logistic")
  )
  binary_new <- binary[1:8, c("patient_id", "gene_a", "gene_b")]
  binary_new$outcome <- binary$outcome[1:8]
  binary_new$time <- seq_len(8L)
  binary_new$event <- rep(c(0L, 1L), 4L)
  original <- sb_predict_final(binary_fit, binary_new, "patient_id")
  changed <- binary_new
  changed$outcome <- c(NA, Inf, "changed", rep("irrelevant", 5L))
  changed$time <- list(NULL, NaN, Inf, -Inf, "x", NA, list(1), factor("a"))
  changed$event <- rev(changed$outcome)
  observed <- sb_predict_final(binary_fit, changed, "patient_id")

  expect_identical(observed, original)

  survival <- make_sb_fixture(n = 60L, seed = 20260901L)
  survival_task <- sb_surv_task(
    survival,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 12
  )
  survival_fit <- sb_fit_final(survival_task, sb_pipeline())
  survival_new <- survival[1:8, c(
    "patient_id", "gene_a", "gene_b", "time", "event"
  )]
  first <- sb_predict_final(survival_fit, survival_new, "patient_id")
  survival_new$time <- rep(NA_character_, 8L)
  survival_new$event <- rep(Inf, 8L)
  second <- sb_predict_final(survival_fit, survival_new, "patient_id")

  expect_identical(second, first)
})

test_that("prediction restores fitted feature order and keeps row order", {
  d <- make_sb_binary_fixture(n = 60L, seed = 20260902L)
  names(d)[names(d) == "gene_a"] <- "gene A"
  names(d)[names(d) == "gene_b"] <- "gene-B"
  d[["gene-B"]] <- as.integer(seq_len(nrow(d)) %% 7L)
  predictors <- c("gene A", "gene-B")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  rows <- c(12L, 1L, 7L, 4L)
  new_data <- d[rows, c("gene-B", "patient_id", "gene A"), drop = FALSE]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  oracle <- sb_predict_logistic_learner_state(
    fit$learner_state,
    new_data[predictors]
  )

  expect_identical(prediction$status, "success")
  expect_identical(
    prediction$predictions$patient_id,
    as.character(d$patient_id[rows])
  )
  expect_equal(
    prediction$predictions$probability,
    oracle$probability,
    tolerance = 0
  )
})

test_that("invalid prediction inputs fail visibly without fallback", {
  d <- make_sb_binary_fixture(n = 60L, seed = 20260903L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  valid <- d[1:6, c("patient_id", "gene_a", "gene_b")]

  missing <- valid[c("patient_id", "gene_a")]
  expect_final_prediction_failure(
    sb_predict_final(fit, missing, "patient_id"),
    "stormbreaker_missing_final_features",
    6L
  )

  duplicate_feature <- valid
  duplicate_feature$copy <- duplicate_feature$gene_a
  names(duplicate_feature)[names(duplicate_feature) == "copy"] <- "gene_a"
  expect_final_prediction_failure(
    sb_predict_final(fit, duplicate_feature, "patient_id"),
    "stormbreaker_duplicate_final_feature_columns",
    6L
  )

  duplicate_id_column <- valid
  duplicate_id_column$copy <- duplicate_id_column$patient_id
  names(duplicate_id_column)[names(duplicate_id_column) == "copy"] <-
    "patient_id"
  expect_final_prediction_failure(
    sb_predict_final(fit, duplicate_id_column, "patient_id"),
    "stormbreaker_duplicate_prediction_id_column",
    0L
  )

  duplicate_ids <- valid
  duplicate_ids$patient_id[2L] <- duplicate_ids$patient_id[1L]
  expect_final_prediction_failure(
    sb_predict_final(fit, duplicate_ids, "patient_id"),
    "stormbreaker_duplicate_prediction_id",
    0L
  )

  for (bad_id in list(NA_character_, "")) {
    bad <- valid
    bad$patient_id[2L] <- bad_id
    expect_final_prediction_failure(
      sb_predict_final(fit, bad, "patient_id"),
      "stormbreaker_bad_prediction_id",
      0L
    )
  }

  nonnumeric <- valid
  nonnumeric$gene_a <- as.character(nonnumeric$gene_a)
  expect_final_prediction_failure(
    sb_predict_final(fit, nonnumeric, "patient_id"),
    "stormbreaker_bad_adapter_predictor_type",
    6L
  )

  for (bad_value in list(NA_real_, Inf)) {
    invalid <- valid
    invalid$gene_a[1L] <- bad_value
    expect_final_prediction_failure(
      sb_predict_final(fit, invalid, "patient_id"),
      "stormbreaker_bad_adapter_predictor_values",
      6L
    )
  }
})

test_that("failed fits and malformed states remain failed predictions", {
  d <- make_sb_binary_fixture(n = 10L, seed = 20260904L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  failed_fit <- sb_fit_final(
    task,
    sb_pipeline(
      learner = "glmnet_logistic",
      inner_folds = 6L,
      inner_seed = 29L
    )
  )
  new_data <- d[1:4, c("patient_id", "gene_a", "gene_b")]
  expect_final_prediction_failure(
    sb_predict_final(failed_fit, new_data, "patient_id"),
    "stormbreaker_failed_final_fit",
    4L
  )

  good_fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  bad_state <- good_fit
  bad_state$learner_state$coefficient[1L] <- NA_real_
  expect_final_prediction_failure(
    sb_predict_final(bad_state, new_data, "patient_id"),
    "stormbreaker_bad_learner_state",
    4L
  )

  testthat::local_mocked_bindings(
    sb_predict_final_state = function(state, x) {
      list(probability = rep.int(2, nrow(x)))
    }
  )
  invalid_output <- sb_predict_final(good_fit, new_data, "patient_id")
  expect_final_prediction_failure(
    invalid_output,
    "stormbreaker_invalid_final_probability",
    4L
  )
})

test_that("a corrupt saved tree fails at state loading", {
  skip_if_not_installed("lightgbm")

  d <- make_sb_binary_fixture(n = 60L, seed = 20260905L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  fit <- sb_fit_final(
    task,
    sb_pipeline(
      learner = "lightgbm_binary",
      lightgbm_nrounds = 3L,
      lightgbm_num_leaves = 5L,
      lightgbm_min_data_in_leaf = 2L
    )
  )
  fit$learner_state$model_string <- "not a LightGBM model"
  new_data <- d[1:4, c("patient_id", "gene_a", "gene_b")]
  prediction <- sb_predict_final(fit, new_data, "patient_id")

  expect_final_prediction_failure(
    prediction,
    "stormbreaker_invalid_lightgbm_state",
    4L
  )
})

test_that("final prediction printing reports no validation performance", {
  d <- make_sb_binary_fixture(n = 40L, seed = 20260906L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  prediction <- sb_predict_final(
    fit,
    d[1:5, c("patient_id", "gene_a", "gene_b")],
    "patient_id"
  )
  printed <- paste(capture.output(print(prediction)), collapse = "\n")

  expect_match(printed, "final state-only prediction", fixed = TRUE)
  expect_match(printed, fit$pipeline$id, fixed = TRUE)
  expect_match(printed, "Patients: 5", fixed = TRUE)
  expect_false(grepl("AUC|Brier|concordance|validation performance", printed))
})

test_that("survival final prediction carries declared event-risk thresholds", {
  d <- make_sb_fixture(n = 60L, seed = 20260910L)
  predictors <- c("gene_a", "gene_b")
  thresholds <- c(0.20, 0.50)
  plain_task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month",
    prediction_horizon = 12
  )
  threshold_task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month",
    prediction_horizon = 12,
    decision_thresholds = thresholds
  )
  plain_fit <- sb_fit_final(plain_task, sb_pipeline())
  threshold_fit <- sb_fit_final(threshold_task, sb_pipeline())
  new_data <- d[1:8, c("patient_id", predictors), drop = FALSE]
  plain_prediction <- sb_predict_final(
    plain_fit,
    new_data,
    "patient_id"
  )
  threshold_prediction <- sb_predict_final(
    threshold_fit,
    new_data,
    "patient_id"
  )

  expect_identical(
    threshold_fit$development$endpoint$decision_thresholds,
    thresholds
  )
  expect_identical(threshold_fit$prediction$decision_thresholds, thresholds)
  expect_identical(threshold_prediction$decision_thresholds, thresholds)
  expect_equal(
    threshold_prediction$predictions$risk,
    plain_prediction$predictions$risk,
    tolerance = 0
  )
  expect_equal(
    threshold_prediction$predictions$survival_probability,
    plain_prediction$predictions$survival_probability,
    tolerance = 0
  )
  expect_false("decision_thresholds" %in% names(plain_prediction))
  expect_false(any(c(
    "predicted_class", "decision", "action", "treat", "selected"
  ) %in% names(threshold_prediction$predictions)))
})

test_that("survival final-fit threshold copies must remain consistent", {
  d <- make_sb_fixture(n = 60L, seed = 20260911L)
  predictors <- c("gene_a", "gene_b")
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 12,
    decision_thresholds = c(0.20, 0.50)
  )
  fit <- sb_fit_final(task, sb_pipeline())
  new_data <- d[1:5, c("patient_id", predictors), drop = FALSE]

  missing_copy <- fit
  missing_copy$development$endpoint$decision_thresholds <- NULL
  expect_final_prediction_failure(
    sb_predict_final(missing_copy, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    5L
  )

  inconsistent_copy <- fit
  inconsistent_copy$development$endpoint$decision_thresholds <- c(0.20, 0.60)
  expect_final_prediction_failure(
    sb_predict_final(inconsistent_copy, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    5L
  )

  invalid_copy <- fit
  invalid_copy$prediction$decision_thresholds <- c(0.20, 1)
  expect_final_prediction_failure(
    sb_predict_final(invalid_copy, new_data, "patient_id"),
    "stormbreaker_bad_decision_thresholds",
    5L
  )
})
})

# ---- pipeline-fit-core ----
local({
test_that("shared selector fitting uses the declared endpoint", {
  survival_data <- make_sb_fixture(n = 60L, seed = 20260819L)
  survival_x <- survival_data[c("gene_a", "gene_b", "gene_c")]
  survival_endpoint <- list(
    task_type = "survival",
    time = survival_data$time,
    event = as.integer(survival_data$event)
  )
  survival_selector <- sb_univ_cox_top_k_selector_adapter(2L)
  survival_fit <- sb_fit_runtime_selector(
    survival_selector,
    survival_x,
    survival_endpoint
  )
  direct_survival <- survival_selector$fit(
    survival_x,
    survival_endpoint$time,
    survival_endpoint$event
  )

  expect_identical(
    survival_fit$selected_features,
    direct_survival$selected_features
  )
  expect_identical(survival_fit$scores, direct_survival$scores)

  binary_data <- sb_example_binary(n = 60L, seed = 20260819L)
  binary_x <- binary_data[c("gene_a", "gene_b", "gene_c")]
  binary_endpoint <- list(
    task_type = "binary",
    outcome = as.integer(binary_data$outcome)
  )
  binary_selector <- sb_univ_binary_auc_top_k_selector_adapter(2L)
  binary_fit <- sb_fit_runtime_selector(
    binary_selector,
    binary_x,
    binary_endpoint
  )
  direct_binary <- binary_selector$fit(
    binary_x,
    binary_endpoint$outcome
  )

  expect_identical(
    binary_fit$selected_features,
    direct_binary$selected_features
  )
  expect_identical(binary_fit$scores, direct_binary$scores)

  none_fit <- sb_fit_runtime_selector(
    sb_none_selector_adapter(),
    binary_x,
    binary_endpoint
  )
  expect_identical(none_fit$selected_features, names(binary_x))
})

test_that("shared learner preparation reruns glmnet inner folds on its input", {
  skip_if_not_installed("glmnet")

  survival_data <- make_sb_fixture(n = 60L, seed = 20260819L)
  survival_x <- survival_data[c("gene_a", "gene_b", "gene_c")]
  survival_ids <- as.character(survival_data$patient_id)
  survival_pipeline <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 101L
  )
  survival_runtime <- sb_resolve_pipeline_runtime(survival_pipeline)
  survival_prepared <- sb_prepare_runtime_learner(
    survival_runtime$learner,
    survival_pipeline,
    survival_x,
    list(
      task_type = "survival",
      time = survival_data$time,
      event = as.integer(survival_data$event)
    ),
    survival_ids
  )

  expect_identical(survival_prepared$learner$id, "glmnet_cox")
  expect_identical(
    names(survival_prepared$tuning$fold_id),
    survival_ids
  )
  expect_length(survival_prepared$tuning$fold_id, nrow(survival_data))

  binary_data <- sb_example_binary(n = 60L, seed = 20260819L)
  binary_x <- binary_data[c("gene_a", "gene_b", "gene_c")]
  binary_ids <- as.character(binary_data$patient_id)
  binary_pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 101L
  )
  binary_runtime <- sb_resolve_pipeline_runtime(binary_pipeline)
  binary_prepared <- sb_prepare_runtime_learner(
    binary_runtime$learner,
    binary_pipeline,
    binary_x,
    list(
      task_type = "binary",
      outcome = as.integer(binary_data$outcome)
    ),
    binary_ids
  )

  expect_identical(binary_prepared$learner$id, "glmnet_logistic")
  expect_identical(names(binary_prepared$tuning$fold_id), binary_ids)
  expect_length(binary_prepared$tuning$fold_id, nrow(binary_data))
})

test_that("shared learner fitting dispatches survival and binary endpoints", {
  survival_data <- make_sb_fixture(n = 40L, seed = 20260819L)
  survival_x <- survival_data[c("gene_a", "gene_b")]
  survival_model <- sb_fit_runtime_learner(
    sb_coxph_learner_adapter(),
    survival_x,
    list(
      task_type = "survival",
      time = survival_data$time,
      event = as.integer(survival_data$event)
    )
  )
  expect_s3_class(survival_model, "coxph")

  binary_data <- sb_example_binary(n = 60L, seed = 20260819L)
  binary_x <- binary_data[c("gene_a", "gene_b")]
  binary_model <- sb_fit_runtime_learner(
    sb_logistic_learner_adapter(),
    binary_x,
    list(
      task_type = "binary",
      outcome = as.integer(binary_data$outcome)
    )
  )
  expect_s3_class(binary_model, "glm")
})
})

# ---- learner-reuse-execution ----
local({
# ---- cv-prediction-reuse ----

test_that("survival CV uses each bagged-kNN assessment curve once", {
  skip_if_not_installed("bnnSurvival")
  skip_if_not_installed("survivalmodels")
  skip_if_not_installed("distr6")
  d <- sb_example_survival(90L, 82L)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month",
    prediction_horizon = as.numeric(median(d$time)))
  bagged <- sb_predict_bagged_knn_survival_curves
  akritas <- sb_akritas_curves
  calls <- c(bagged = 0L, akritas = 0L)
  local_mocked_bindings(
    sb_predict_bagged_knn_survival_curves = function(...) {
      calls[["bagged"]] <<- calls[["bagged"]] + 1L
      bagged(...)
    },
    sb_akritas_curves = function(...) {
      calls[["akritas"]] <<- calls[["akritas"]] + 1L
      akritas(...)
    })
  for (id in c("bagged_knn_survival", "akritas_survival")) {
    result <- sb_run_cv(task, sb_pipeline("none", id), folds = 3L, seed = 82L)
    expect_identical(result$complete_oof, TRUE)
    expect_identical(result$complete_survival_probability_oof, TRUE)
  }
  expect_identical(calls, c(bagged = 3L, akritas = 6L))
})

test_that("shared assessment predictions preserve complete CV and horizon-only failures", {
  skip_if_not_installed("bnnSurvival")
  skip_if_not_installed("survivalmodels")
  skip_if_not_installed("distr6")
  withr::local_seed(17L)
  d <- sb_example_survival(90L, 82L)
  # One unique maximum makes exactly its assessment fold unable to support h.
  d$time[1L] <- max(d$time) + 1
  features <- c("gene_a", "gene_b", "gene_c")
  bagged <- sb_bagged_knn_survival_learner_adapter
  akritas <- sb_akritas_learner_adapter
  plain <- function(task, pipeline) {
    with_mocked_bindings(sb_run_cv(task, pipeline, folds = 3L, seed = 82L),
      sb_bagged_knn_survival_learner_adapter = function(...) {
        a <- bagged(...); a$predict_details <- NULL; a
      },
      sb_akritas_learner_adapter = function(...) {
        a <- akritas(...); a$predict_details <- NULL; a
      })
  }
  before <- .Random.seed
  for (h in c(median(d$time), max(d$time))) {
    task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
      prediction_horizon = as.numeric(h))
    for (id in c("bagged_knn_survival", "akritas_survival")) {
      pipelines <- list(sb_pipeline("none", id),
        sb_pipeline("univ_cox_top_k", id, selector_n_features = 2L),
        sb_pipeline("fixed_features", id, selector_features = features[1:2], feature_set_id = "two"))
      for (p in pipelines) {
        expected <- plain(task, p)
        actual <- sb_run_cv(task, p, folds = 3L, seed = 82L)
        expect_identical(actual, expected)
        expect_identical(actual$complete_oof, TRUE)
        expect_identical(actual$complete_survival_probability_oof, h < max(d$time))
        if (h == max(d$time)) {
          unsupported <- vapply(actual$baseline_hazards, function(item) {
            identical(item$error_class, "stormbreaker_prediction_horizon_unsupported")
          }, logical(1))
          expect_identical(sum(unsupported), 1L)
        }
        expect_identical(.Random.seed, before)
      }
    }
  }
  d$gene_b <- d$gene_a
  singular <- sb_surv_task(d, "patient_id", "time", "event", features[1:2], "month",
    prediction_horizon = as.numeric(median(d$time)))
  p <- sb_pipeline("none", "bagged_knn_survival")
  expected <- plain(singular, p)
  actual <- sb_run_cv(singular, p, folds = 3L, seed = 82L)
  expect_identical(actual, expected)
  expect_identical(actual$complete_oof, FALSE)
})

test_that("bagged-kNN shared curves preserve direct horizon queries", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(72L, 42L)
  d$time <- rep(seq_len(12L), length.out = nrow(d))
  for (p in c(1L, 3L)) {
    x <- d[c("gene_a", "gene_b", "gene_c")][seq_len(p)]
    adapter <- sb_bagged_knn_survival_learner_adapter("bagged_knn_survival")
    model <- adapter$fit(x, d$time, d$event)
    query <- x[c(1L, 9L, 30L), , drop = FALSE]
    details <- adapter$predict_details(model, query)
    expect_identical(details$risk, adapter$predict(model, query))
    for (h in c(0.5, 1, 1.5, 6, 6.75, 12)) {
      expected <- adapter$predict_survival(model, x, d$time, d$event, query, h)
      actual <- adapter$predict_survival(model, x, d$time, d$event, query, h,
        prediction_details = details)
      expect_identical(actual, expected)
    }
  }
})

test_that("nonfatal curve warnings are preserved and prevent prediction reuse", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(72L, 42L)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month",
    prediction_horizon = as.numeric(median(d$time)))
  original <- sb_predict_bagged_knn_survival_curves
  factory <- sb_bagged_knn_survival_learner_adapter
  local_mocked_bindings(sb_predict_bagged_knn_survival_curves = function(...) {
    warning("prediction warning", call. = FALSE)
    original(...)
  })
  run <- function(shared) {
    warnings <- character()
    value <- withCallingHandlers(with_mocked_bindings(
      sb_run_cv(task, sb_pipeline("none", "bagged_knn_survival"), folds = 3L, seed = 82L),
      sb_bagged_knn_survival_learner_adapter = function(...) {
        a <- factory(...)
        if (!shared) a$predict_details <- NULL
        a
      }), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(value = value, warnings = warnings)
  }
  actual <- run(TRUE)
  expected <- run(FALSE)
  expect_identical(actual, expected)
  expect_identical(actual$warnings, rep("prediction warning", 6L))
})

test_that("risk-only CV does not request shared prediction details", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(72L, 42L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month")
  factory <- sb_bagged_knn_survival_learner_adapter
  local_mocked_bindings(sb_bagged_knn_survival_learner_adapter = function(...) {
    a <- factory(...)
    a$predict_details <- function(...) stop("Risk-only CV must not request details")
    a
  })
  result <- sb_run_cv(task, sb_pipeline("none", "bagged_knn_survival"), folds = 3L, seed = 82L)
  expect_identical(result$complete_oof, TRUE)
})

test_that("native bagged-kNN prediction errors are not converted to successful risks", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(72L, 42L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month",
    prediction_horizon = as.numeric(median(d$time)))
  p <- sb_pipeline("none", "bagged_knn_survival")
  factory <- sb_bagged_knn_survival_learner_adapter
  local_mocked_bindings(sb_predict_bagged_knn_engine = function(...) {
    warning("native engine warning", call. = FALSE)
    stop("native engine prediction failed")
  })
  actual <- sb_run_cv(task, p, folds = 3L, seed = 82L)
  expected <- with_mocked_bindings(sb_run_cv(task, p, folds = 3L, seed = 82L),
    sb_bagged_knn_survival_learner_adapter = function(...) {
      a <- factory(...); a$predict_details <- NULL; a
    })
  expect_identical(actual, expected)
  expect_identical(actual$complete_oof, FALSE)
  expect_identical(actual$runs$error_class, rep("simpleError", 3L))
  expect_identical(actual$runs$status, rep("failed", 3L))
  expect_true(all(grepl("native engine prediction failed", actual$runs$error_message, fixed = TRUE)))
})

test_that("ordinary bagged-kNN warnings preserve valid shared and plain predictions", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(72L, 42L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month",
    prediction_horizon = as.numeric(median(d$time)))
  p <- sb_pipeline("none", "bagged_knn_survival")
  native <- sb_predict_bagged_knn_engine
  factory <- sb_bagged_knn_survival_learner_adapter
  baseline <- sb_run_cv(task, p, folds = 3L, seed = 82L)
  local_mocked_bindings(sb_predict_bagged_knn_engine = function(...) {
    warning("native engine warning", call. = FALSE)
    native(...)
  })
  run <- function(shared) {
    messages <- character()
    value <- withCallingHandlers({
      if (shared) sb_run_cv(task, p, folds = 3L, seed = 82L) else
        with_mocked_bindings(sb_run_cv(task, p, folds = 3L, seed = 82L),
          sb_bagged_knn_survival_learner_adapter = function(...) {
            a <- factory(...); a$predict_details <- NULL; a
          })
    }, warning = function(w) {
      messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    list(value = value, messages = messages)
  }
  shared <- run(TRUE)
  plain <- run(FALSE)
  expect_identical(shared$value, plain$value)
  expect_identical(shared$value, baseline)
  expect_true(length(shared$messages) > 0L)
  expect_true(length(plain$messages) > 0L)
  expect_true(all(shared$messages == "native engine warning"))
  expect_true(all(plain$messages == "native engine warning"))
  expect_identical(shared$value$complete_oof, TRUE)
})

# ---- final-cv-learner-reuse ----

make_final_cv_reuse_task <- function() {
  withr::with_seed(8901L, {
    x <- as.data.frame(matrix(stats::rnorm(90L * 12L), nrow = 90L))
    names(x) <- sprintf("gene_%02d", seq_len(12L))
    lp <- x[[1]] - x[[2]] + 0.6 * x[[3]]
    et <- stats::rexp(90L, rate = exp(lp) / 15)
    ct <- stats::rexp(90L, rate = 1 / 25)
    d <- data.frame(patient_id = paste0("patient_", seq_len(90L)),
      time = pmax(pmin(et, ct), 0.05), event = as.integer(et <= ct), x)
    sb_surv_task(d, "patient_id", "time", "event", names(x), "month",
      prediction_horizon = 7)
  })
}

test_that("identical final CV learners share complete fitting and diagnostics", {
  skip_if_not_installed("CoxBoost")
  skip_if_not_installed("plsRcox")
  skip_if_not_installed("penAFT")
  task <- make_final_cv_reuse_task()
  native <- sb_fit_runtime_learner
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  withr::local_seed(38L)
  rng <- .Random.seed
  for (id in c("coxboost_cv", "plsrcox_cv", "penaft_gehan_elastic_net")) {
    args <- if (id == "penaft_gehan_elastic_net") list() else
      list(learner_tuning_folds = 3L, learner_seed = if (id == "coxboost_cv") 8902L else 8903L)
    pp <- list(do.call(sb_pipeline, c(list("none", id), args)),
      do.call(sb_pipeline, c(list("fixed_features", id, selector_features = task$predictors,
        feature_set_id = "all"), args)))
    run <- function(enabled) sb_with_selector_reuse(lapply(pp, function(p) {
      sb_fit_final(task, p, capture_diagnostics = TRUE)
    }), enabled = enabled, diagnostics = TRUE)
    calls <- 0L
    plain <- run(FALSE)
    expect_identical(calls, 2L, info = id)
    calls <- 0L
    shared <- run(TRUE)
    expect_identical(calls, 1L, info = id)
    expect_identical(shared$reuse$final_learner_hits, 1L, info = id)
    for (i in 1:2) {
      a <- plain$result[[i]]; b <- shared$result[[i]]
      expect_identical(b$status, "success", info = id)
      expect_equal(b, a, info = id)
      expect_identical(b$learner_state, a$learner_state, info = id)
      expect_identical(b$tuning, a$tuning, info = id)
      expect_identical(sb_predict_final(unserialize(serialize(b, NULL)),
        task$data[c("patient_id", task$predictors)], "patient_id"),
        sb_predict_final(a, task$data[c("patient_id", task$predictors)], "patient_id"), info = id)
    }
    expect_identical(.Random.seed, rng, info = id)
  }
})

test_that("stored component-search warnings prevent final model reuse", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(...) {
    calls <<- calls + 1L
    structure(list(value = 1L), stormbreaker_fit_warnings = "late-component warning")
  })
  learner <- sb_plsrcox_cv_learner_adapter()
  x <- data.frame(a = 1:6, b = 6:1)
  endpoint <- list(task_type = "survival", time = 1:6, event = rep(0:1, 3))
  z <- sb_with_selector_reuse({
    lapply(1:2, function(i) sb_fit_final_runtime_learner(
      learner, x, endpoint, letters[1:6], "none"))
  }, diagnostics = TRUE)
  expect_identical(calls, 2L)
  expect_identical(z$reuse$final_learner_hits, 0L)
  expect_identical(z$reuse$final_learner_entries, 0L)
  expect_identical(vapply(z$result, function(m) attr(m, "stormbreaker_fit_warnings"), ""),
    rep("late-component warning", 2L))
})

test_that("final CV learner reuse distinguishes seeds and tuning folds", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(...) {
    calls <<- calls + 1L
    list(value = calls)
  })
  x <- data.frame(a = 1:6, b = 6:1)
  e <- list(task_type = "survival", time = 1:6, event = rep(0:1, 3))
  for (l in list(sb_coxboost_cv_learner_adapter(), sb_plsrcox_cv_learner_adapter(),
    sb_rule_aft_learner_adapter("penaft_gehan_elastic_net"))) {
    calls <- 0L
    run <- function(learner = l) sb_fit_final_runtime_learner(learner, x, e, letters[1:6], "none")
    z <- sb_with_selector_reuse({
      first <- run()
      expect_identical(run(), first)
      changed <- l; changed$params$seed <- changed$params$seed + 1L; run(changed)
      changed <- l
      if (l$id == "penaft_gehan_elastic_net") changed$params$internal_folds <- 3L else
        changed$params$inner_folds <- 3L
      run(changed)
    }, diagnostics = TRUE)
    expect_identical(calls, 3L)
    expect_identical(z$reuse$final_learner_hits, 1L)
    expect_identical(z$reuse$final_learner_entries, 3L)
  }
})

test_that("native PLS search warnings remain on every unshared final route", {
  skip_if_not_installed("plsRcox")
  task <- make_final_cv_reuse_task()
  message <- "Loglik converged before variable  3 ; coefficient may be infinite. "
  cv_calls <- 0L
  local_mocked_bindings(
    cv.plsRcox = function(...) {
      cv_calls <<- cv_calls + 1L
      warning(message)
      list(lambda.min5 = 2L, folds = split(seq_len(90L), rep(1:3, each = 30L)),
        cv.error5 = seq(0, 1, length.out = 11L), cv.se5 = rep(0.01, 11L),
        completed.cv5 = matrix(TRUE, 11L, 3L))
    },
    plsRcox = function(Xplan, nt, ...) {
      list(Coeffs = matrix(0.1, nrow = ncol(Xplan), ncol = 1L,
        dimnames = list(names(Xplan), NULL)), computed_nt = nt)
    }, .package = "plsRcox")
  pp <- list(sb_pipeline("none", "plsrcox_cv", learner_tuning_folds = 3L),
    sb_pipeline("fixed_features", "plsrcox_cv", learner_tuning_folds = 3L,
      selector_features = task$predictors, feature_set_id = "all"))
  run <- function(enabled) sb_with_selector_reuse(
    lapply(pp, function(p) sb_fit_final(task, p)), enabled = enabled, diagnostics = TRUE)
  plain <- run(FALSE)
  cv_calls <- 0L
  shared <- run(TRUE)
  expect_identical(cv_calls, 2L)
  expect_identical(shared$reuse$final_learner_entries, 0L)
  expect_identical(shared$reuse$final_learner_hits, 0L)
  expect_equal(shared$result, plain$result)
  for (fit in shared$result) {
    expect_identical(fit$status, "success")
    expect_identical(fit$tuning$warnings, paste0("Learner refit: ", message))
  }
})

test_that("final tuned learners preserve absent RNG and later diagnostic requests", {
  skip_if_not_installed("CoxBoost")
  skip_if_not_installed("plsRcox")
  skip_if_not_installed("penAFT")
  task <- make_final_cv_reuse_task()
  withr::local_seed(39L)
  rm(".Random.seed", envir = .GlobalEnv)
  for (id in c("coxboost_cv", "plsrcox_cv", "penaft_gehan_elastic_net")) {
    args <- if (id == "penaft_gehan_elastic_net") list() else
      list(learner_tuning_folds = 3L, learner_seed = if (id == "coxboost_cv") 8902L else 8903L)
    p <- do.call(sb_pipeline, c(list("none", id), args))
    q <- do.call(sb_pipeline, c(list("fixed_features", id,
      selector_features = task$predictors, feature_set_id = "all"), args))
    expected <- sb_fit_final(task, q, capture_diagnostics = TRUE)
    actual <- sb_with_selector_reuse({
      first <- sb_fit_final(task, p)
      expect_null(first$diagnostic_artifacts)
      sb_fit_final(task, q, capture_diagnostics = TRUE)
    }, diagnostics = TRUE)
    expect_identical(actual$result$status, "success", info = id)
    expect_identical(actual$reuse$final_learner_hits, 1L, info = id)
    expect_equal(actual$result, expected, info = id)
    expect_identical(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), FALSE,
      info = id)
  }
})

# ---- final-learned-selector-reuse ----

test_that("learned and declared selectors share identical final SVM inputs", {
  skip_if_not_installed("survivalsvm")
  d <- sb_example_survival(72L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
    prediction_horizon = 7)
  p <- sb_pipeline("univ_cox_top_k", "survivalsvm", selector_n_features = 2L)
  first <- sb_fit_final(task, p)
  expect_identical(first$status, "success")
  q <- sb_pipeline("fixed_features", "survivalsvm",
    selector_features = first$required_features, feature_set_id = "same_selected")
  native <- sb_fit_survivalsvm
  calls <- 0L
  local_mocked_bindings(sb_fit_survivalsvm = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  run <- function(enabled) sb_with_selector_reuse(
    lapply(list(p, q), function(pipeline) sb_fit_final(task, pipeline,
      capture_diagnostics = TRUE)), enabled = enabled, diagnostics = TRUE)
  withr::local_seed(42L)
  before <- .Random.seed
  plain <- run(FALSE)
  expect_identical(calls, 2L)
  calls <- 0L
  shared <- run(TRUE)
  expect_identical(calls, 1L)
  expect_identical(shared$reuse$final_learner_hits, 1L)
  expect_identical(.Random.seed, before)
  for (i in 1:2) {
    expect_identical(shared$result[[i]]$status, "success")
    expect_equal(shared$result[[i]], plain$result[[i]])
    expect_identical(shared$result[[i]]$learner_state, plain$result[[i]]$learner_state)
    expect_identical(sb_predict_final(unserialize(serialize(shared$result[[i]], NULL)),
      d[c("patient_id", features)], "patient_id"),
      sb_predict_final(plain$result[[i]], d[c("patient_id", features)], "patient_id"))
  }
})

# ---- learner-reuse-analysis-cv ----

test_that("analysis-local searches reuse the complete fitted model", {
  for (package in c("SLOPE", "randomGLM", "cancerclass", "BART", "CoxBoost",
    "plsRcox", "superpc", "gbm", "pre", "penAFT")) skip_if_not_installed(package)
  # Load before both arms: an installed-package build-version warning is not
  # a fit warning and is emitted only on the first namespace load.
  if (!"package:gbm" %in% search()) {
    suppressWarnings(suppressPackageStartupMessages(library("gbm", character.only = TRUE)))
    withr::defer(detach("package:gbm", unload = FALSE))
  }
  features <- c("gene_a", "gene_b", "gene_c")
  calls <- 0L
  native <- sb_fit_runtime_learner
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L
    native(learner, x, endpoint)
  })
  binary_ids <- c("random_glm_binary", "cancerclass_binary",
    "pre_binary_rule_ensemble")
  survival_ids <- c("bart_survival", "bart_aft_survival", "coxboost_cv",
    "plsrcox_cv", "superpc_cv", "gbm_cox_cv", "pre_survival_rule_ensemble",
    "penaft_gehan_elastic_net")
  for (id in c(survival_ids, binary_ids)) {
    features <- c("gene_a", "gene_b", "gene_c")
    if (id %in% binary_ids) {
      d <- sb_example_binary(120L, 1208L)
      task <- sb_binary_task(d, "patient_id", "outcome", features)
    } else {
      d <- sb_example_survival(180L, 11032L)
      task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
        prediction_horizon = 7)
    }
    pipelines <- list(sb_pipeline("none", id), sb_pipeline("fixed_features", id,
      selector_features = features, feature_set_id = "all"))
    for (final in c(FALSE, TRUE)) {
      run <- function(enabled) withr::with_seed(101L, {
        warnings <- character()
        value <- withCallingHandlers({
          if (final) sb_with_selector_reuse(lapply(pipelines,
            function(p) sb_fit_final(task, p, capture_diagnostics = TRUE)), enabled = enabled)
          else sb_compare_cv(task, pipelines, folds = 3L, seed = 11L, reuse_selectors = enabled)
        }, warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
        list(value = value, warnings = warnings, rng = .Random.seed)
      })
      calls <- 0L; plain <- run(FALSE); plain_calls <- calls
      calls <- 0L; shared <- run(TRUE)
      expect_equal(shared, plain, info = id)
      expect_identical(plain_calls, if (final) 2L else 6L, info = id)
      # This SuperPC fixture has one genuine native nonconvergence fold.
      # Both routes must retry/report that failed fit, sharing only successes.
      expected_calls <- if (final) 1L else if (id == "superpc_cv") 4L else 3L
      expect_identical(calls, expected_calls, info = id)
      if (final) {
        expect_identical(vapply(shared$value, `[[`, "", "status"), rep("success", 2L), info = id)
        for (i in 1:2) expect_identical(
          sb_predict_final(unserialize(serialize(shared$value[[i]], NULL)),
            d[c("patient_id", features)], "patient_id"),
          sb_predict_final(plain$value[[i]], d[c("patient_id", features)], "patient_id"), info = id)
      } else {
        expect_identical(unique(shared$value$comparison$status),
          if (id == "superpc_cv") "failed" else "success", info = id)
        if (id == "superpc_cv") expect_identical(
          shared$value$results[[1L]]$runs$error_class[[1L]],
          "stormbreaker_superpc_cv_warning")
      }
    }
  }
})

# ---- learner-reuse-cv ----

test_that("identical learner fits share within, never across, analysis sets", {
  skip_if_not_installed("survivalsvm")
  skip_if_not_installed("mboost")
  skip_if_not_installed("plsRglm")
  native <- sb_fit_runtime_learner
  calls <- list()
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls[[length(calls) + 1L]] <<- rownames(x)
    native(learner, x, endpoint)
  })
  features <- c("gene_a", "gene_b", "gene_c")
  binary_ids <- c("logistic", "gaussian_naive_bayes", sb_specialized_binary_linear_learner_ids())
  for (id in c(binary_ids, "coxph", sb_survivalsvm_learner_ids())) {
    if (id %in% binary_ids) {
      d <- sb_example_binary(120L, 1208L)
      task <- sb_binary_task(d, "patient_id", "outcome", features)
    } else {
      d <- sb_example_survival(120L, 11032L)
      task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
        prediction_horizon = 7)
    }
    pipelines <- list(sb_pipeline("none", id), sb_pipeline("fixed_features", id,
      selector_features = features, feature_set_id = "all"))
    calls <- list()
    plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 11L, reuse_selectors = FALSE)
    expect_length(calls, 6L)
    calls <- list()
    shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 11L)
    expect_length(calls, 3L)
    expect_equal(shared, plain)
    expect_identical(unique(shared$comparison$status), "success")
    expected <- lapply(shared$results[[1]]$folds, function(f) rownames(task$data)[f$analysis_rows])
    expect_identical(calls, expected)
    run <- function(enabled) sb_with_selector_reuse(
      lapply(pipelines, function(p) sb_fit_final(task, p, capture_diagnostics = TRUE)),
      enabled = enabled)
    calls <- list(); plain_final <- run(FALSE)
    expect_length(calls, 2L)
    calls <- list(); shared_final <- run(TRUE)
    expect_length(calls, 1L)
    expect_equal(shared_final, plain_final)
    for (i in seq_along(pipelines)) expect_identical(
      sb_predict_final(unserialize(serialize(shared_final[[i]], NULL)),
        d[c("patient_id", features)], "patient_id"),
      sb_predict_final(plain_final[[i]], d[c("patient_id", features)], "patient_id"))
  }
})

test_that("CV reuse separates full input identities and retries warning or failed fits", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L
    if (identical(learner$params$mode, "warning")) warning("native warning")
    if (identical(learner$params$mode, "error")) stop("native failure")
    list(coefficient = c(a = 1, b = 2))
  })
  learner <- sb_logistic_learner_adapter()
  x <- data.frame(a = 1:6, b = 6:1)
  ids <- paste0("p", 1:6)
  endpoint <- list(task_type = "binary", outcome = rep(0:1, 3))
  result <- sb_with_selector_reuse({
    fit <- function(l = learner, z = x, e = endpoint, p = ids) {
      sb_fit_reusable_runtime_learner(l, z, e, p, "none")
    }
    first <- fit(); changed <- first; changed$coefficient[] <- 0
    expect_identical(fit(), first)
    fit(z = x[6:1, ]); fit(z = x[c("b", "a")]); fit(p = rev(ids))
    changed_endpoint <- endpoint; changed_endpoint$outcome <- 1L - endpoint$outcome
    fit(e = changed_endpoint)
    changed_learner <- learner; changed_learner$params <- list(mode = "different")
    fit(l = changed_learner)
    for (i in 1:2) {
      changed_learner$params$mode <- "warning"
      expect_warning(fit(l = changed_learner), "native warning")
      changed_learner$params$mode <- "error"
      expect_error(fit(l = changed_learner), "native failure")
    }
    sb_fit_final_runtime_learner(learner, x, endpoint, ids, "none")
  }, diagnostics = TRUE)
  expect_identical(calls, 11L)
  expect_identical(result$reuse$cv_learner_hits, 1L)
  expect_equal(result$reuse$cv_learner_entries, 6L)
  expect_identical(result$reuse$final_learner_misses, 1L)
  expect_null(.sb_selector_reuse_context$active)
})

# ---- learner-reuse-mutable ----

test_that("shared boosting fits do not alias mutable native model environments", {
  skip_if_not_installed("mboost")
  d <- sb_example_survival(80L, 11032L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "survival", time = d$time, event = d$event)
  learner <- sb_mboost_cox_learner_adapter("glmboost_cox")
  result <- sb_with_selector_reuse({
    fit <- function() sb_fit_reusable_runtime_learner(learner, x, endpoint, d$patient_id, "none")
    first <- fit()
    original <- learner$predict(first, x)
    mboost::mstop(first$engine) <- 1L
    second <- fit()
    expect_identical(as.integer(mboost::mstop(second$engine)), 100L)
    expect_identical(learner$predict(second, x), original)
    mboost::mstop(second$engine) <- 2L
    third <- fit()
    expect_identical(as.integer(mboost::mstop(third$engine)), 100L)
    expect_identical(learner$predict(third, x), original)
  }, diagnostics = TRUE)
  expect_identical(result$reuse$cv_learner_hits, 2L)
})

test_that("shared XGBoost fits own independent native boosters", {
  skip_if_not_installed("xgboost")
  d <- sb_example_binary(80L, 1208L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "binary", outcome = d$outcome)
  learner <- sb_xgboost_learner_adapter("xgboost_binary")
  for (phase in c("cv", "final")) {
    result <- sb_with_selector_reuse({
      fit <- function() sb_fit_reusable_runtime_learner(learner, x, endpoint,
        d$patient_id, "none", phase = phase)
      first <- fit()
      original <- learner$predict(first, x)
      xgboost::xgb.attr(first$booster, "consumer_note") <- "first"
      second <- fit()
      expect_null(xgboost::xgb.attr(second$booster, "consumer_note"))
      expect_identical(learner$predict(second, x), original)
      xgboost::xgb.attr(second$booster, "consumer_note") <- "second"
      third <- fit()
      expect_null(xgboost::xgb.attr(third$booster, "consumer_note"))
      expect_identical(learner$predict(third, x), original)
      expect_identical(xgboost::xgb.attr(first$booster, "consumer_note"), "first")
    }, diagnostics = TRUE)
    expect_identical(result$reuse[[paste0(phase, "_learner_hits")]], 2L)
  }
})

test_that("shared LightGBM fits isolate native iteration changes", {
  skip_if_not_installed("lightgbm")
  d <- sb_example_binary(100L, 1208L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "binary", outcome = d$outcome)
  learner <- sb_lightgbm_binary_learner_adapter()
  out <- sb_with_selector_reuse({
    fit <- function() sb_fit_reusable_runtime_learner(learner, x, endpoint,
      d$patient_id, "none")
    first <- fit()
    original <- learner$predict(first, x)
    original_text <- first$booster$save_model_to_string()
    first$booster$rollback_one_iter()
    expect_false(identical(first$booster$save_model_to_string(), original_text))
    second <- fit()
    expect_identical(second$booster$save_model_to_string(), original_text)
    expect_identical(learner$predict(second, x), original)
    # Loaded boosters need not retain rollback history. Use the already changed
    # native booster to verify stale source-precision diagnostics are not reused.
    second$booster <- first$booster
    expect_false(identical(second$booster$save_model_to_string(), original_text))
    changed <- second
    changed$reused_native_importance <- NULL
    expect_identical(sb_extract_lightgbm_tree_model(second),
      sb_extract_lightgbm_tree_model(changed))
    third <- fit()
    expect_identical(third$booster$save_model_to_string(), original_text)
    expect_identical(learner$predict(third, x), original)
  }, diagnostics = TRUE)
  expect_identical(out$reuse$cv_learner_hits, 2L)
})

# ---- learner-reuse-retention ----

test_that("final fits reuse a nonconsecutive identical analysis input", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L
    list(features = names(x), means = colMeans(x))
  })
  learner <- list(id = "logistic", params = list(), fit = function(...) NULL)
  x <- data.frame(a = 1:6, b = 6:1)
  endpoint <- list(task_type = "binary", outcome = rep(0:1, 3L))
  ids <- paste0("p", 1:6)
  out <- sb_with_selector_reuse({
    first <- sb_fit_final_runtime_learner(learner, x, endpoint, ids, "none")
    different <- sb_fit_final_runtime_learner(learner, x["a"], endpoint, ids, "fixed_features")
    again <- sb_fit_final_runtime_learner(learner, x, endpoint, ids, "none")
    expect_identical(again, first)
    expect_identical(different$features, "a")
  }, diagnostics = TRUE)
  expect_identical(calls, 2L)
  expect_identical(out$reuse$final_learner_hits, 1L)
  expect_identical(out$reuse$final_learner_entries, 2L)
})

test_that("retained learner payloads obey a shared CV and final memory budget", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L
    list(value = sum(x), payload = raw(20000L))
  })
  learner <- list(id = "logistic", params = list(), fit = function(...) NULL)
  endpoint <- list(task_type = "binary", outcome = rep(0:1, 3L))
  ids <- paste0("p", 1:6)
  out <- sb_with_selector_reuse({
    context <- .sb_selector_reuse_context$active
    context$learner_byte_limit <- 70000
    fit <- function(v, phase = "final") sb_fit_reusable_runtime_learner(
      learner, data.frame(a = seq_len(6L) + v), endpoint, ids, "none", phase)
    a <- fit(0); fit(10, "cv"); expect_identical(fit(0), a)
    fit(20) # A was recently touched: evict the older CV entry, not A.
    expect_identical(fit(0), a)
    expect_identical(calls, 3L)
    fit(10, "cv")
    expect_identical(calls, 4L)
    expect_lte(context$learner_bytes, context$learner_byte_limit)
    expect_gt(context$learner_evictions, 0L)
  }, diagnostics = TRUE)
  expect_lte(out$reuse$learner_bytes, 70000)
  expect_null(.sb_selector_reuse_context$active)
})

test_that("oversized learner states are returned but not retained", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(...) {
    calls <<- calls + 1L
    list(payload = raw(20000L))
  })
  learner <- list(id = "logistic", params = list(), fit = function(...) NULL)
  out <- sb_with_selector_reuse({
    .sb_selector_reuse_context$active$learner_byte_limit <- 1000
    fit <- function() sb_fit_final_runtime_learner(learner, data.frame(a = 1:6),
      list(outcome = rep(0:1, 3L)), paste0("p", 1:6), "none")
    expect_identical(fit(), fit())
  }, diagnostics = TRUE)
  expect_identical(calls, 2L)
  expect_identical(out$reuse$final_learner_entries, 0L)
  expect_identical(out$reuse$learner_oversized, 2L)
})

test_that("eviction preserves native final model predictions", {
  d <- sb_example_binary(100L, 1208L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  a <- sb_pipeline("none", "logistic")
  b <- sb_pipeline("fixed_features", "logistic",
    selector_features = features[1:2], feature_set_id = "two")
  plain <- lapply(list(a, b, a), function(p) sb_fit_final(task, p))
  shared <- sb_with_selector_reuse({
    first <- sb_fit_final(task, a)
    context <- .sb_selector_reuse_context$active
    context$learner_byte_limit <- context$learner_bytes
    list(first, sb_fit_final(task, b), sb_fit_final(task, a))
  }, diagnostics = TRUE)
  expect_gt(shared$reuse$learner_evictions, 0L)
  expect_lte(shared$reuse$learner_bytes, shared$reuse$learner_byte_limit)
  expect_equal(shared$result, plain)
  for (i in 1:3) expect_identical(
    sb_predict_final(unserialize(serialize(shared$result[[i]], NULL)),
      d[c("patient_id", features)], "patient_id"),
    sb_predict_final(plain[[i]], d[c("patient_id", features)], "patient_id"))
})

# ---- learner-reuse-tuning ----

test_that("matching glmnet preparation shares native tuning and refits", {
  skip_if_not_installed("glmnet")
  preparations <- 0L
  # Observe the actual native lambda searches, not only the wrapper invocation.
  native_binary <- sb_select_glmnet_logistic_lambda
  native_survival <- sb_select_glmnet_cox_lambda
  local_mocked_bindings(
    sb_select_glmnet_logistic_lambda = function(...) {
      preparations <<- preparations + 1L; native_binary(...)
    },
    sb_select_glmnet_cox_lambda = function(...) {
      preparations <<- preparations + 1L; native_survival(...)
    })
  features <- c("gene_a", "gene_b", "gene_c")
  for (type in c("binary", "survival")) {
    if (type == "binary") {
      d <- sb_example_binary(120L, 1208L)
      task <- sb_binary_task(d, "patient_id", "outcome", features)
      id <- "glmnet_logistic"
    } else {
      d <- sb_example_survival(120L, 11032L)
      task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
        prediction_horizon = 7)
      id <- "glmnet_cox"
    }
    pipelines <- list(sb_pipeline("none", id, inner_folds = 3L, inner_seed = 101L),
      sb_pipeline("fixed_features", id, selector_features = features,
        feature_set_id = "all", inner_folds = 3L, inner_seed = 101L))
    preparations <- 0L
    plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 11L, reuse_selectors = FALSE)
    expect_identical(preparations, 6L)
    preparations <- 0L
    shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 11L)
    expect_identical(preparations, 3L)
    expect_identical(shared, plain)
    expect_identical(unique(shared$comparison$status), "success")
    run <- function(enabled) sb_with_selector_reuse(
      lapply(pipelines, function(p) sb_fit_final(task, p)), enabled = enabled,
      diagnostics = TRUE)
    preparations <- 0L; plain_final <- run(FALSE)
    expect_identical(preparations, 2L)
    preparations <- 0L; shared_final <- run(TRUE)
    expect_identical(preparations, 1L)
    expect_equal(shared_final$result, plain_final$result)
    expect_identical(shared_final$reuse$final_learner_hits, 1L)
    for (i in seq_along(pipelines)) expect_identical(
      sb_predict_final(unserialize(serialize(shared_final$result[[i]], NULL)),
        d[c("patient_id", features)], "patient_id"),
      sb_predict_final(plain_final$result[[i]], d[c("patient_id", features)], "patient_id"))
  }
})

test_that("preparation identity includes upstream data, selection and grouping", {
  skip_if_not_installed("glmnet")
  d <- sb_example_binary(120L, 1208L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "binary", outcome = d$outcome)
  p <- sb_pipeline("univ_binary_auc_top_k", "glmnet_logistic",
    selector_n_features = 2L, inner_folds = 3L, inner_seed = 101L)
  runtime <- sb_resolve_pipeline_runtime(p)
  selected <- sb_fit_runtime_selector_uncached(runtime$selector, x, endpoint)
  selected_x <- runtime$selector$apply(selected, x)
  sc <- list(selector = runtime$selector, fitted_selector = selected, x = x)
  native <- sb_select_glmnet_logistic_lambda
  calls <- 0L
  local_mocked_bindings(sb_select_glmnet_logistic_lambda = function(...) {
    calls <<- calls + 1L; native(...)
  })
  result <- sb_with_selector_reuse({
    prepare <- function(pipeline = p, context = sc, groups = d$patient_id) {
      sb_prepare_runtime_learner(sb_resolve_pipeline_runtime(pipeline)$learner,
        pipeline, selected_x, endpoint, d$patient_id, context, groups)
    }
    a <- prepare(); b <- prepare()
    expect_equal(a$tuning, b$tuning)
    fixed <- sb_pipeline("fixed_features", "glmnet_logistic",
      selector_features = names(selected_x), feature_set_id = "same_outer_genes",
      inner_folds = 3L, inner_seed = 101L)
    prepare(fixed)
    changed <- sc
    omitted <- setdiff(names(x), names(selected_x))
    changed$x[[omitted]] <- -changed$x[[omitted]]
    prepare(context = changed)
    prepare(groups = rev(d$patient_id))
  }, diagnostics = TRUE)
  expect_identical(calls, 4L)
  expect_identical(result$reuse$preparation_hits, 1L)
})

test_that("preparation warnings and failures are never hidden by reuse", {
  calls <- 0L
  mode <- "warning"
  local_mocked_bindings(sb_prepare_runtime_learner_uncached = function(...) {
    calls <<- calls + 1L
    if (mode == "error") stop("native tuning failure")
    list(learner = list(), tuning = list(warnings = "native retained warning"))
  })
  p <- sb_pipeline("none", "glmnet_logistic", inner_folds = 3L, inner_seed = 1L)
  learner <- sb_resolve_pipeline_runtime(p)$learner
  result <- sb_with_selector_reuse({
    prepare <- function() sb_prepare_runtime_learner(learner, p,
      data.frame(a = 1:6, b = 6:1), list(task_type = "binary", outcome = rep(0:1, 3)),
      as.character(1:6))
    for (i in 1:2) expect_identical(prepare()$tuning$warnings, "native retained warning")
    mode <- "error"
    for (i in 1:2) expect_error(prepare(), "native tuning failure")
  }, diagnostics = TRUE)
  expect_identical(calls, 4L)
  expect_identical(result$reuse$preparation_hits, 0L)
})

# ---- native-hotspot-reuse ----

test_that("Lasso selection reuses the full analysis path already fitted by CV", {
  skip_if_not_installed("glmnet")
  d <- sb_example_survival(120L, 11032L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  native <- glmnet::glmnet
  calls <- 0L
  local_mocked_bindings(glmnet = function(...) {
    calls <<- calls + 1L
    native(...)
  }, .package = "glmnet")
  actual <- sb_fit_lasso_nonzero(x, d$time, d$event, 3L, 17L,
    capture_history = TRUE)
  expect_identical(calls, 4L)
  model <- sb_glmnet_cox_learner_adapter(1, actual$lambda)$fit(x, d$time, d$event)
  coefficient <- drop(as.matrix(stats::coef(model, s = actual$lambda)))
  expect_identical(actual$scores$coefficient, unname(coefficient))
  expect_identical(actual$coefficient_path, list(lambda = model$lambda, beta = model$beta))
  expect_identical(actual$selected_features, names(coefficient)[coefficient != 0])
})

test_that("RSF selector avoids proximity without changing native selection", {
  skip_if_not_installed("randomForestSRC")
  d <- sb_example_survival(90L, 11032L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  native <- randomForestSRC::rfsrc
  requested <- list()
  local_mocked_bindings(rfsrc = function(...) {
    args <- list(...)
    requested[[length(requested) + 1L]] <<- args$proximity
    native(...)
  }, .package = "randomForestSRC")
  withr::local_seed(112L)
  before <- .Random.seed
  actual <- sb_fit_rsf_minimal_depth(x, d$time, d$event, 5L, 17L)
  after <- .Random.seed
  expect_identical(requested[[1L]], FALSE)
  original <- sb_fit_rsf_minimal_depth
  body(original) <- parse(text = sub("proximity = FALSE", "proximity = TRUE",
    paste(deparse(body(original)), collapse = "\n"), fixed = TRUE))[[1L]]
  assign(".Random.seed", before, .GlobalEnv)
  expected <- original(x, d$time, d$event, 5L, 17L)
  expect_identical(actual, expected)
  expect_identical(.Random.seed, after)
})

test_that("RSF native depth traversal can use private reference containers", {
  skip_if_not_installed("randomForestSRC")
  d <- sb_example_survival(90L, 11032L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  data <- sb_additional_survival_formula_data(x, d$time, d$event,
    random_forest_src = TRUE)
  engine <- withr::with_options(list(rf.cores = 1L, mc.cores = 1L),
    randomForestSRC::rfsrc(data$formula, data$data, ntree = 30L,
      nodesize = 5L, splitrule = "logrank", importance = TRUE, seed = 17L))
  ns <- asNamespace("randomForestSRC")
  parser <- get("rfsrcParseTree", ns)
  maximum <- randomForestSRC::max.subtree
  before <- .Random.seed
  expected <- withr::with_options(list(rf.cores = 1L, mc.cores = 1L),
    randomForestSRC::var.select(engine, conservative = "high", verbose = FALSE))
  after <- .Random.seed
  assign(".Random.seed", before, .GlobalEnv)
  actual <- withr::with_options(list(rf.cores = 1L, mc.cores = 1L),
    sb_rsf_minimal_depth_selection(engine))
  expect_identical(actual, expected)
  expect_identical(.Random.seed, after)
  expect_identical(get("rfsrcParseTree", ns), parser)
  expect_identical(randomForestSRC::max.subtree, maximum)
})
})

# ---- learner-reuse-native ----
local({
# ---- additional-survival-curve-reuse ----

test_that("saved BART and cforest states evaluate a native curve once", {
  skip_if_not_installed("BART")
  skip_if_not_installed("partykit")
  skip_if_not_installed("coin")
  d <- sb_example_survival(60L, 992L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  query <- x[c(3L, 12L, 27L, 59L), , drop = FALSE]
  bart <- sb_bart_survival_curves
  party <- sb_partykit_survival_distributions
  calls <- 0L
  local_mocked_bindings(
    sb_bart_survival_curves = function(...) { calls <<- calls + 1L; bart(...) },
    sb_partykit_survival_distributions = function(...) { calls <<- calls + 1L; party(...) })
  withr::local_seed(93L)
  for (id in c("bart_survival", "cforest_survival")) {
    a <- sb_additional_survival_learner_adapter(id)
    model <- a$fit(x, d$time, d$event)
    for (h in c(0.001, as.numeric(median(d$time)), as.numeric(max(d$time)))) {
      state <- unserialize(serialize(a$extract_state(model, list(x = x, time = d$time,
        event = d$event, prediction_horizon = h, time_unit = "month")), NULL))
      expected <- list(risk = a$predict(model, query), prediction_horizon = h,
        time_unit = "month", survival_probability = a$predict_survival(
          model, x, d$time, d$event, query, h)$probability)
      calls <- 0L; rng <- .Random.seed
      actual <- a$predict_state(state, query)
      expect_identical(actual, expected, info = id)
      expect_identical(calls, 1L, info = id)
      expect_identical(.Random.seed, rng, info = id)
    }
  }
})

test_that("curve-capable adapters expose exact CV details only for pure prediction", {
  skip_if_not_installed("BART")
  skip_if_not_installed("partykit")
  skip_if_not_installed("coin")
  d <- sb_example_survival(60L, 992L)
  x <- d[c("gene_a", "gene_b")]
  for (id in c("bart_survival", "cforest_survival")) {
    a <- sb_additional_survival_learner_adapter(id)
    expect_identical(is.function(a$predict_details), TRUE, info = id)
    if (!is.function(a$predict_details)) next
    model <- a$fit(x, d$time, d$event)
    details <- a$predict_details(model, x[1:3, , drop = FALSE])
    expect_identical(details$risk, a$predict(model, x[1:3, , drop = FALSE]))
    for (h in c(0.001, as.numeric(median(d$time)))) {
      expected <- a$predict_survival(model, x, d$time, d$event, x[1:3, , drop = FALSE], h)
      actual <- a$predict_survival(model, x, d$time, d$event, x[1:3, , drop = FALSE], h,
        prediction_details = details)
      expect_identical(actual, expected, info = id)
    }
  }
  for (id in c("rsf_survival", "ctree_survival", "aorsf_survival", "rpart_survival",
    "survreg_weibull")) {
    expect_null(sb_additional_survival_learner_adapter(id)$predict_details)
  }
})

test_that("shared native curves preserve CV and horizon-only failure boundaries", {
  skip_if_not_installed("BART")
  skip_if_not_installed("partykit")
  skip_if_not_installed("coin")
  d <- sb_example_survival(72L, 992L)
  d$time[1L] <- max(d$time) + 1
  factory <- sb_additional_survival_learner_adapter
  bart <- sb_bart_survival_curves
  party <- sb_partykit_survival_distributions
  calls <- 0L
  local_mocked_bindings(
    sb_bart_survival_curves = function(...) { calls <<- calls + 1L; bart(...) },
    sb_partykit_survival_distributions = function(...) { calls <<- calls + 1L; party(...) })
  withr::local_seed(94L)
  before <- .Random.seed
  for (id in c("bart_survival", "cforest_survival")) {
    for (h in c(as.numeric(median(d$time)), as.numeric(max(d$time)))) {
      task <- sb_surv_task(d, "patient_id", "time", "event",
        c("gene_a", "gene_b", "gene_c"), "month", prediction_horizon = h)
      p <- sb_pipeline("univ_cox_top_k", id, selector_n_features = 2L)
      calls <- 0L
      actual <- sb_run_cv(task, p, folds = 3L, seed = 82L)
      expect_identical(calls, 6L, info = id)
      calls <- 0L
      expected <- with_mocked_bindings(sb_run_cv(task, p, folds = 3L, seed = 82L),
        sb_additional_survival_learner_adapter = function(...) {
          a <- factory(...); a$predict_details <- NULL; a
        })
      expect_identical(calls, if (h < max(d$time)) 9L else 8L, info = id)
      expect_identical(actual, expected, info = id)
      expect_identical(actual$complete_oof, TRUE, info = id)
      expect_identical(actual$complete_survival_probability_oof, h < max(d$time), info = id)
      expect_identical(.Random.seed, before, info = id)
    }
  }
})

test_that("public final predictions replay shared distributions without outcomes", {
  skip_if_not_installed("BART")
  skip_if_not_installed("partykit")
  skip_if_not_installed("coin")
  d <- sb_example_survival(60L, 992L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
    prediction_horizon = as.numeric(median(d$time)))
  for (id in c("bart_survival", "ctree_survival", "cforest_survival")) {
    fit <- sb_fit_final(task, sb_pipeline("fixed_features", id,
      selector_features = features[1:2], feature_set_id = "two"))
    expect_identical(fit$status, "success", info = id)
    blind <- d[c("patient_id", rev(features))]
    actual <- sb_predict_final(unserialize(serialize(fit, NULL)), blind, "patient_id")
    expect_identical(actual$status, "success", info = id)
    expect_identical(actual, sb_predict_final(fit, blind, "patient_id"), info = id)
    expect_identical(actual$predictions$patient_id, d$patient_id, info = id)
  }
})

test_that("a warned native distribution is not reused during state prediction", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("coin")
  d <- sb_example_survival(60L, 992L)
  x <- d[c("gene_a", "gene_b")]
  a <- sb_additional_survival_learner_adapter("cforest_survival")
  model <- a$fit(x, d$time, d$event)
  state <- a$extract_state(model, list(x = x, time = d$time, event = d$event,
    prediction_horizon = 5, time_unit = "month"))
  expected <- a$predict_state(state, x[1:3, , drop = FALSE])
  original <- sb_partykit_survival_distributions
  calls <- 0L
  local_mocked_bindings(sb_partykit_survival_distributions = function(...) {
    calls <<- calls + 1L
    warning("native prediction warning", call. = FALSE)
    original(...)
  })
  warnings <- character()
  actual <- withCallingHandlers(a$predict_state(state, x[1:3, , drop = FALSE]),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  expect_identical(actual, expected)
  expect_identical(calls, 2L)
  expect_identical(warnings, rep("native prediction warning", 2L))
})

test_that("risk-only state calls once and native errors remain failures", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("coin")
  d <- sb_example_survival(60L, 992L)
  x <- d[c("gene_a", "gene_b")]
  a <- sb_additional_survival_learner_adapter("ctree_survival")
  model <- a$fit(x, d$time, d$event)
  state <- a$extract_state(model)
  original <- sb_partykit_survival_distributions
  calls <- 0L
  local_mocked_bindings(sb_partykit_survival_distributions = function(...) {
    calls <<- calls + 1L; original(...)
  })
  expected <- a$predict(model, x)
  calls <- 0L
  withr::local_seed(95L)
  rm(".Random.seed", envir = .GlobalEnv)
  expect_identical(a$predict_state(state, x), list(risk = expected))
  expect_identical(calls, 1L)
  expect_identical(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), FALSE)
  local_mocked_bindings(sb_partykit_survival_distributions = function(...) {
    sb_abort("native distribution failed", "stormbreaker_test_prediction_failure")
  })
  error <- capture_sb_error(a$predict_state(state, x))
  expect_s3_class(error, "stormbreaker_test_prediction_failure")
})

# ---- cforest-prediction-batches ----

test_that("survival forest predictions agree across assessment batch sizes", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("coin")
  withr::local_seed(123)
  x <- data.frame(a = rnorm(66), b = rnorm(66), c = rnorm(66))
  time <- rexp(66, exp(x$a) / 10)
  event <- rep(c(1L, 1L, 0L), 22)
  adapter <- sb_resolve_pipeline_runtime(sb_pipeline("none", "cforest_survival"))$learner
  model <- adapter$fit(x, time, event)
  large <- x[rep(seq_len(66), length.out = 518), , drop = FALSE]
  actual <- sb_partykit_survival_distributions(model, large)
  expected <- lapply(seq_len(nrow(large)), function(i) {
    predict(model$engine, newdata = large[i, , drop = FALSE], type = "prob")[[1L]]
  })
  expect_length(actual, 518L)
  expect_equal(lapply(actual, `[[`, "surv"), lapply(expected, `[[`, "surv"))
  expect_equal(lapply(actual, `[[`, "time"), lapply(expected, `[[`, "time"))
  small <- sb_partykit_survival_distributions(model, x[1:10, ])
  expect_equal(small, predict(model$engine, newdata = x[1:10, ], type = "prob"))
})

# ---- final-ensemble-reuse ----

test_that("identical final GBM and RLT routes share fitting without changing outputs", {
  skip_if_not_installed("gbm")
  skip_if_not_installed("RLT")
  d <- sb_example_survival(90L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
    prediction_horizon = 7)
  gbm_fit <- sb_fit_gbm_cox_cv
  rlt_fit <- sb_fit_rlt_learner
  calls <- 0L
  local_mocked_bindings(
    sb_fit_gbm_cox_cv = function(...) { calls <<- calls + 1L; gbm_fit(...) },
    sb_fit_rlt_learner = function(...) { calls <<- calls + 1L; rlt_fit(...) })
  withr::local_seed(35L)
  before <- .Random.seed
  for (id in c("gbm_cox_cv", "rlt_reinforcement_survival")) {
    args <- if (id == "gbm_cox_cv") list(learner_tuning_folds = 3L,
      gbm_cv_max_trees = 30L) else list()
    pipelines <- list(do.call(sb_pipeline, c(list("none", id), args)),
      do.call(sb_pipeline, c(list("fixed_features", id,
        selector_features = features, feature_set_id = "all"), args)))
    run <- function(enabled) sb_with_selector_reuse(
      lapply(pipelines, function(p) sb_fit_final(task, p, capture_diagnostics = TRUE)),
      enabled = enabled, diagnostics = TRUE)
    calls <- 0L
    plain <- run(FALSE)
    expect_identical(calls, 2L, info = id)
    calls <- 0L
    shared <- run(TRUE)
    expect_identical(calls, 1L, info = id)
    expect_identical(shared$reuse$final_learner_hits, 1L, info = id)
    for (i in seq_along(pipelines)) {
      expect_identical(shared$result[[i]]$status, "success", info = id)
      expect_equal(shared$result[[i]], plain$result[[i]], info = id)
      expect_identical(shared$result[[i]]$learner_state,
        plain$result[[i]]$learner_state, info = id)
      expect_identical(sb_predict_final(unserialize(serialize(shared$result[[i]], NULL)),
        d[c("patient_id", features)], "patient_id"),
        sb_predict_final(plain$result[[i]], d[c("patient_id", features)], "patient_id"), info = id)
    }
    expect_identical(.Random.seed, before, info = id)
  }
  expect_null(.sb_selector_reuse_context$active)
})

test_that("final ensemble reuse keeps tuning settings and failures distinct", {
  x <- data.frame(a = 1:6, b = 6:1)
  e <- list(task_type = "survival", time = 1:6, event = rep(0:1, 3))
  ids <- letters[1:6]
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, ...) {
    calls <<- calls + 1L
    if (identical(learner$params$mode, "warning")) warning("native warning")
    if (identical(learner$params$mode, "error")) sb_abort("native failure", "native_failure")
    list(call = calls)
  })
  for (learner in list(sb_gbm_cox_cv_learner_adapter(),
    sb_rlt_learner_adapter("rlt_reinforcement_survival"))) {
    calls <- 0L
    run <- function(l = learner, selector = "none") {
      sb_fit_final_runtime_learner(l, x, e, ids, selector)
    }
    z <- sb_with_selector_reuse({
      first <- run()
      expect_identical(run(selector = "fixed_features"), first)
      altered <- learner; altered$params$seed <- 9L; run(altered)
      altered <- learner
      if (learner$id == "gbm_cox_cv") altered$params$inner_folds <- 3L else
        altered$params$num_trees <- 20L
      run(altered)
      run(selector = "univ_cox_top_k")
      run(selector = "univ_cox_top_k")
      altered <- learner; altered$params$mode <- "warning"
      for (i in 1:2) expect_warning(run(altered), "native warning")
      altered$params$mode <- "error"
      for (i in 1:2) expect_error(run(altered), class = "native_failure")
    }, diagnostics = TRUE)
    expect_identical(calls, 7L)
    expect_identical(z$reuse$final_learner_hits, 3L)
    expect_identical(z$reuse$final_learner_entries, 3L)
  }
})

test_that("GBM CV is not reused across assessment folds or horizon extraction", {
  skip_if_not_installed("gbm")
  d <- sb_example_survival(180L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- function(h) sb_surv_task(d, "patient_id", "time", "event", features,
    "month", prediction_horizon = h)
  p <- sb_pipeline("none", "gbm_cox_cv", learner_tuning_folds = 3L, gbm_cv_max_trees = 30L)
  native <- sb_fit_gbm_cox_cv
  calls <- 0L
  local_mocked_bindings(sb_fit_gbm_cox_cv = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  plain <- lapply(c(3, 9), function(h) sb_fit_final(task(h), p))
  calls <- 0L
  shared <- sb_with_selector_reuse({
    fits <- lapply(c(3, 9), function(h) sb_fit_final(task(h), p))
    cv <- sb_run_cv(task(3), p, folds = 3L, seed = 27L, progress = FALSE)
    list(fits = fits, cv = cv)
  }, diagnostics = TRUE)
  expect_identical(calls, 4L)
  expect_identical(shared$reuse$final_learner_hits, 1L)
  expect_equal(shared$result$fits, plain)
  expect_identical(shared$result$cv$complete_oof, TRUE)
  expect_identical(shared$result$cv$complete_survival_probability_oof, TRUE)
  expect_identical(shared$result$cv,
    sb_run_cv(task(3), p, folds = 3L, seed = 27L, progress = FALSE))
  a <- sb_predict_final(shared$result$fits[[1]], d[c("patient_id", features)], "patient_id")
  b <- sb_predict_final(shared$result$fits[[2]], d[c("patient_id", features)], "patient_id")
  expect_identical(a$predictions$risk, b$predictions$risk)
  expect_identical(identical(a$predictions$survival_probability, b$predictions$survival_probability), FALSE)
})

test_that("cohort comparisons reuse identical CV and final GBM fits separately", {
  skip_if_not_installed("gbm")
  features <- c("gene_a", "gene_b", "gene_c")
  make_task <- function(seed, prefix) {
    d <- sb_example_survival(180L, seed)
    d$patient_id <- paste0(prefix, d$patient_id)
    sb_surv_task(d, "patient_id", "time", "event", features, "month")
  }
  development <- make_task(11032L, "dev")
  external <- make_task(11033L, "external")
  pipelines <- list(sb_pipeline("none", "gbm_cox_cv", learner_tuning_folds = 3L,
      gbm_cv_max_trees = 30L), sb_pipeline("fixed_features", "gbm_cox_cv",
      selector_features = features, feature_set_id = "all", learner_tuning_folds = 3L,
      gbm_cv_max_trees = 30L))
  native <- sb_fit_gbm_cox_cv
  calls <- 0L
  local_mocked_bindings(sb_fit_gbm_cox_cv = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  shared <- sb_compare_cohorts(development, list(External = external),
    pipelines = pipelines, folds = 3L, seed = 27L, progress = FALSE)
  expect_identical(calls, 4L)
  expect_identical(unname(vapply(shared$fits, `[[`, "", "status")), rep("success", 2L))
  calls <- 0L
  local_mocked_bindings(sb_fit_final_runtime_learner = function(learner, x, endpoint, ...) {
    sb_fit_runtime_learner(learner, x, endpoint)
  })
  plain <- sb_compare_cohorts(development, list(External = external),
    pipelines = pipelines, folds = 3L, seed = 27L, progress = FALSE)
  expect_identical(calls, 5L)
  expect_equal(shared$fits, plain$fits)
  expect_identical(shared$predictions, plain$predictions)
  expect_equal(shared$cohorts, plain$cohorts)
})

test_that("seeded final ensembles do not create caller RNG state", {
  skip_if_not_installed("gbm")
  skip_if_not_installed("RLT")
  d <- sb_example_survival(90L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  withr::local_seed(35L)
  rm(".Random.seed", envir = .GlobalEnv)
  for (p in list(sb_pipeline("none", "gbm_cox_cv", learner_tuning_folds = 3L,
      gbm_cv_max_trees = 30L), sb_pipeline("none", "rlt_reinforcement_survival"))) {
    fits <- sb_with_selector_reuse(list(sb_fit_final(task, p), sb_fit_final(task, p)),
      diagnostics = TRUE)
    expect_identical(fits$result[[1]]$status, "success")
    expect_equal(fits$result[[1]], fits$result[[2]])
    expect_identical(fits$reuse$final_learner_hits, 1L)
    expect_identical(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), FALSE)
  }
})

# ---- final-survivalsvm-reuse ----

test_that("identical none and fixed-feature final SVM fits share only native fitting", {
  skip_if_not_installed("survivalsvm")
  d <- sb_example_survival(72L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
    prediction_horizon = 7)
  pipelines <- list(sb_pipeline("none", "survivalsvm"),
    sb_pipeline("fixed_features", "survivalsvm",
      selector_features = features, feature_set_id = "declared"))
  native <- sb_fit_survivalsvm
  calls <- 0L
  local_mocked_bindings(sb_fit_survivalsvm = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  run <- function(enabled) sb_with_selector_reuse(
    lapply(pipelines, function(p) sb_fit_final(task, p, capture_diagnostics = TRUE)),
    enabled = enabled, diagnostics = TRUE)
  plain <- run(FALSE)
  expect_identical(calls, 2L)
  calls <- 0L
  shared <- run(TRUE)
  expect_identical(calls, 1L)
  expect_identical(shared$reuse$final_learner_hits, 1L)
  for (i in seq_along(pipelines)) {
    a <- plain$result[[i]]; b <- shared$result[[i]]
    expect_identical(b$status, "success")
    expect_identical(b$pipeline, pipelines[[i]])
    expect_equal(b, a)
    expect_identical(sb_predict_final(unserialize(serialize(b, NULL)),
      d[c("patient_id", features)], "patient_id"),
      sb_predict_final(unserialize(serialize(a, NULL)),
        d[c("patient_id", features)], "patient_id"))
  }
  expect_null(.sb_selector_reuse_context$active)
})

test_that("all fixed survival-SVM variants replay exactly without consuming RNG", {
  skip_if_not_installed("survivalsvm")
  d <- sb_example_survival(72L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
    prediction_horizon = 7)
  native <- sb_fit_survivalsvm
  calls <- 0L
  local_mocked_bindings(sb_fit_survivalsvm = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  withr::local_seed(21L)
  before <- .Random.seed
  for (id in sb_survivalsvm_learner_ids()) {
    p <- sb_pipeline("none", id)
    q <- sb_pipeline("fixed_features", id, selector_features = features,
      feature_set_id = "all")
    plain <- sb_fit_final(task, q, capture_diagnostics = TRUE)
    count <- calls
    shared <- sb_with_selector_reuse({
      first <- sb_fit_final(task, p)
      second <- sb_fit_final(task, q, capture_diagnostics = TRUE)
      expect_null(first$diagnostic_artifacts)
      second
    }, diagnostics = TRUE)
    expect_identical(calls - count, 1L, info = id)
    expect_identical(plain$status, "success", info = id)
    expect_identical(shared$result$status, "success", info = id)
    expect_equal(shared$result, plain, info = id)
    expect_identical(sb_predict_final(unserialize(serialize(shared$result, NULL)),
      d[c("patient_id", features)], "patient_id"),
      sb_predict_final(plain, d[c("patient_id", features)], "patient_id"), info = id)
    expect_identical(.Random.seed, before, info = id)
  }
  rm(".Random.seed", envir = .GlobalEnv)
  sb_with_selector_reuse({
    sb_fit_final(task, sb_pipeline("none", "survivalsvm"))
    sb_fit_final(task, sb_pipeline("none", "survivalsvm"))
  })
  expect_identical(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), FALSE)
})

test_that("final learner reuse distinguishes complete data and retains distinct fits", {
  learner <- sb_survivalsvm_learner_adapter()
  x <- data.frame(a = 1:6, b = 6:1)
  endpoint <- list(task_type = "survival", time = 1:6, event = rep(0:1, 3),
    group_ids = NULL)
  ids <- paste0("patient", 1:6)
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(...) {
    calls <<- calls + 1L
    list(value = calls)
  })
  run <- function(l = learner, data = x, e = endpoint, patients = ids,
                  selector = "none") {
    sb_fit_final_runtime_learner(l, data, e, patients, selector)
  }
  out <- sb_with_selector_reuse({
    a <- run()
    expect_identical(run(selector = "fixed_features"), a)
    run(data = x[c("b", "a")])
    run(data = x[6:1, ])
    changed <- x; changed$a[1] <- 99
    run(data = changed)
    changed <- x; rownames(changed) <- ids
    run(data = changed)
    e <- endpoint; e$event[1] <- 1L
    run(e = e)
    e <- endpoint; e$time[1] <- 1.5
    run(e = e)
    e <- endpoint; e$group_ids <- rep(1:3, 2)
    run(e = e)
    run(patients = paste0("other", ids))
    l <- learner; l$params$gamma_mu <- 2
    run(l = l)
    l <- learner; body(l$fit) <- quote(stop("changed implementation"))
    run(l = l)
    run()
  }, diagnostics = TRUE)
  expect_identical(calls, 11L)
  expect_identical(out$reuse$final_learner_hits, 2L)
  expect_identical(out$reuse$final_learner_entries, 11L)
  expect_identical(out$reuse$hits, 0L)
})

test_that("warnings and failures are retried while successful eligible fits reuse", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, ...) {
    calls <<- calls + 1L
    if (learner$params$mode == "warning") warning("native warning")
    if (learner$params$mode == "error") stop("native failure")
    list(value = calls)
  })
  l <- sb_survivalsvm_learner_adapter()
  x <- data.frame(a = 1:4); e <- list(task_type = "survival", time = 1:4, event = c(1, 0, 1, 0))
  ids <- letters[1:4]
  run <- function(selector = "none") sb_fit_final_runtime_learner(l, x, e, ids, selector)
  z <- sb_with_selector_reuse({
    l$params$mode <- "warning"
    for (i in 1:2) expect_warning(run(), "native warning")
    l$params$mode <- "error"
    for (i in 1:2) expect_error(run(), "native failure")
    l$params$mode <- "ordinary"
    l$id <- "gbm_cox_cv"
    for (i in 1:2) run("univ_cox_top_k")
    l$id <- "coxboost"
    for (i in 1:2) run()
  }, diagnostics = TRUE)
  expect_identical(calls, 6L)
  expect_identical(z$reuse$final_learner_entries, 2L)
  expect_identical(z$reuse$final_learner_hits, 2L)
})

test_that("final learner scopes are isolated and CV still fits every assessment fold", {
  skip_if_not_installed("survivalsvm")
  d <- sb_example_survival(90L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  p <- sb_pipeline("none", "survivalsvm")
  native <- sb_fit_survivalsvm
  calls <- 0L
  local_mocked_bindings(sb_fit_survivalsvm = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  out <- sb_with_selector_reuse({
    first <- sb_fit_final(task, p)
    sb_with_selector_reuse(sb_fit_final(task, p), enabled = FALSE)
    sb_with_selector_reuse(sb_fit_final(task, p))
    expect_equal(sb_fit_final(task, p), first)
    cv <- sb_run_cv(task, p, folds = 3L, seed = 21L, progress = FALSE)
    cv
  }, diagnostics = TRUE)
  expect_identical(calls, 6L)
  expect_identical(out$reuse$final_learner_hits, 1L)
  expect_null(.sb_selector_reuse_context$active)
  expect_error(sb_with_selector_reuse({sb_fit_final(task, p); stop("scope exit")}), "scope exit")
  expect_null(.sb_selector_reuse_context$active)
})

test_that("horizon extraction stays route-local after a shared final SVM fit", {
  skip_if_not_installed("survivalsvm")
  d <- sb_example_survival(72L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  make_task <- function(h) sb_surv_task(d, "patient_id", "time", "event", features,
    "month", prediction_horizon = h)
  p <- sb_pipeline("none", "survivalsvm")
  plain <- lapply(c(3, 9), function(h) sb_fit_final(make_task(h), p))
  shared <- sb_with_selector_reuse(
    lapply(c(3, 9), function(h) sb_fit_final(make_task(h), p)), diagnostics = TRUE)
  expect_identical(shared$reuse$final_learner_hits, 1L)
  expect_equal(shared$result, plain)
  a <- sb_predict_final(shared$result[[1]], d[c("patient_id", features)], "patient_id")
  b <- sb_predict_final(shared$result[[2]], d[c("patient_id", features)], "patient_id")
  expect_identical(a$predictions$risk, b$predictions$risk)
  expect_identical(identical(a$predictions$survival_probability, b$predictions$survival_probability), FALSE)
})

test_that("cohort comparisons share matching CV and final fits in separate scopes", {
  skip_if_not_installed("survivalsvm")
  features <- c("gene_a", "gene_b", "gene_c")
  make_task <- function(seed, prefix) {
    d <- sb_example_survival(72L, seed)
    d$patient_id <- paste0(prefix, d$patient_id)
    sb_surv_task(d, "patient_id", "time", "event", features, "month")
  }
  development <- make_task(11032L, "dev")
  external <- make_task(11033L, "external")
  pipelines <- list(sb_pipeline("none", "survivalsvm"),
    sb_pipeline("fixed_features", "survivalsvm",
      selector_features = features, feature_set_id = "all"))
  native <- sb_fit_survivalsvm
  calls <- 0L
  local_mocked_bindings(sb_fit_survivalsvm = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  shared <- sb_compare_cohorts(development, list(External = external),
    pipelines = pipelines, folds = 3L, seed = 21L, progress = FALSE)
  expect_identical(calls, 4L)
  calls <- 0L
  local_mocked_bindings(sb_fit_final_runtime_learner = function(learner, x, endpoint, ...) {
    sb_fit_runtime_learner(learner, x, endpoint)
  })
  plain <- sb_compare_cohorts(development, list(External = external),
    pipelines = pipelines, folds = 3L, seed = 21L, progress = FALSE)
  expect_identical(calls, 5L)
  expect_equal(shared$fits, plain$fits)
  expect_identical(shared$predictions, plain$predictions)
  expect_equal(shared$cohorts, plain$cohorts)
})

# ---- learner-reuse-classical ----

test_that("classical binary final fits retain native results under reuse", {
  d <- sb_example_binary(120L, 1208L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  native <- sb_fit_runtime_learner
  calls <- 0L
  native_calibration <- sb_calibrate_svm_probability
  calibrations <- 0L
  local_mocked_bindings(sb_calibrate_svm_probability = function(...) {
    calibrations <<- calibrations + 1L
    native_calibration(...)
  })
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L
    native(learner, x, endpoint)
  })
  for (id in sb_classical_binary_learner_ids()) {
    pipelines <- list(sb_pipeline("none", id), sb_pipeline("fixed_features", id,
      selector_features = features, feature_set_id = "all"))
    run <- function(enabled) withr::with_seed(101L, {
      value <- sb_with_selector_reuse(lapply(pipelines, function(p) {
        sb_fit_final(task, p, capture_diagnostics = TRUE)
      }), enabled = enabled, diagnostics = TRUE)
      list(value = value, rng = .Random.seed)
    })
    calls <- 0L; calibrations <- 0L; plain <- run(FALSE)
    plain_calls <- calls
    plain_calibrations <- calibrations
    calls <- 0L; calibrations <- 0L; shared <- run(TRUE)
    comparable <- function(value) {
      if (id %in% c("c50_binary", "c50_rules_binary")) {
        clean_clock <- function(engine) {
          for (field in intersect(c("names", "output"), names(engine))) {
            engine[[field]] <- gsub("(?m)^\\| on [^\\n]+$", "| on <native clock>",
              engine[[field]], perl = TRUE)
            engine[[field]] <- gsub("(?m)^(C5\\.0 \\[Release [^\\]]+\\][ \\t]+)[^\\n]+$",
              "\\1<native clock>", engine[[field]], perl = TRUE)
          }
          engine
        }
        value <- lapply(value, function(fit) {
          fit$learner_state$engine_raw <- clean_clock(unserialize(fit$learner_state$engine_raw))
          fit$diagnostic_artifacts$learner$engine <- clean_clock(fit$diagnostic_artifacts$learner$engine)
          fit
        })
      }
      # A serialization round trip changes environment reference encoding,
      # not the restored model; keep full engine and exact prediction checks.
      if (id %in% c("gamboost_binary", "blackboost_binary",
          "gaussian_process_binary", "neuralnet_multilayer_binary",
          "ctree_binary", "cforest_binary", "loclda_binary")) {
        value <- lapply(value, function(fit) {
          if (is.raw(fit$learner_state$engine_raw)) {
            fit$learner_state$engine_raw <- unserialize(fit$learner_state$engine_raw)
          }
          fit
        })
      }
      value
    }
    expect_equal(comparable(shared$value$result), comparable(plain$value$result), info = id)
    expect_identical(shared$rng, plain$rng, info = id)
    if (all(vapply(plain$value$result, function(f) f$status == "success", logical(1)))) {
      expect_identical(shared$value$reuse$final_learner_hits, 1L, info = id)
      expect_identical(calls, plain_calls - 1L, info = id)
      if (id %in% sb_svm_probability_learner_ids()) {
        expect_identical(plain_calibrations, 2L, info = id)
        expect_identical(calibrations, 1L, info = id)
      }
      for (i in seq_along(pipelines)) expect_identical(
        sb_predict_final(unserialize(serialize(shared$value$result[[i]], NULL)),
          d[c("patient_id", features)], "patient_id"),
        sb_predict_final(plain$value$result[[i]], d[c("patient_id", features)], "patient_id"),
        info = id)
    } else {
      expect_identical(shared$value$reuse$final_learner_hits, 0L, info = id)
      expect_identical(calls, plain_calls, info = id)
    }
  }
})

test_that("classical binary CV reuse preserves held-out predictions and failures", {
  d <- sb_example_binary(120L, 1208L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  native <- sb_fit_runtime_learner
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L
    native(learner, x, endpoint)
  })
  for (id in sb_classical_binary_learner_ids()) {
    pipelines <- list(sb_pipeline("none", id), sb_pipeline("fixed_features", id,
      selector_features = features, feature_set_id = "all"))
    run <- function(enabled) withr::with_seed(101L, {
      warnings <- character()
      value <- withCallingHandlers(sb_compare_cv(task, pipelines, folds = 3L,
        seed = 11L, reuse_selectors = enabled), warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      list(value = value, warnings = warnings, rng = .Random.seed)
    })
    calls <- 0L; plain <- run(FALSE); plain_calls <- calls
    calls <- 0L; shared <- run(TRUE)
    expect_equal(shared$value, plain$value, info = id)
    for (i in seq_along(pipelines)) expect_identical(
      shared$value$results[[i]]$predictions, plain$value$results[[i]]$predictions,
      info = id)
    expect_identical(shared$warnings, plain$warnings, info = id)
    expect_identical(shared$rng, plain$rng, info = id)
    expect_lte(calls, plain_calls)
    if (all(shared$value$comparison$status == "success")) {
      expect_identical(calls, 3L, info = id)
    }
  }
})

test_that("SVM fit reuse includes its inner calibration coefficients", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(100L, 1208L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "binary", outcome = d$outcome)
  adapter <- sb_classical_binary_learner_adapter("svm_radial_binary")
  first <- sb_svm_calibrated_runtime_learner(adapter, list(A = -1, B = .1, positive_class = "1"))
  second <- sb_svm_calibrated_runtime_learner(adapter, list(A = -2, B = .2, positive_class = "1"))
  out <- sb_with_selector_reuse({
    a <- sb_fit_reusable_runtime_learner(first, x, endpoint, d$patient_id, "none")
    b <- sb_fit_reusable_runtime_learner(second, x, endpoint, d$patient_id, "none")
    expect_false(identical(first$predict(a, x), second$predict(b, x)))
    sb_fit_reusable_runtime_learner(first, x, endpoint, d$patient_id, "none")
  }, diagnostics = TRUE)
  expect_identical(out$reuse$cv_learner_misses, 2L)
  expect_identical(out$reuse$cv_learner_hits, 1L)
})

# ---- learner-reuse-parametric ----

test_that("AdaBoost reuse preserves native tree CV and caller random state", {
  skip_if_not_installed("fastAdaboost")
  d <- sb_example_binary(80L, 1208L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "binary", outcome = d$outcome)
  learner <- sb_resolve_pipeline_runtime(sb_pipeline("none", "adaboost_binary"))$learner
  for (phase in c("cv", "final")) {
    native <- withr::with_seed(55L, {
      model <- sb_fit_runtime_learner(learner, x, endpoint)
      list(model = model, rng = .Random.seed)
    })
    out <- sb_with_selector_reuse({
      run <- function() withr::with_seed(55L, {
        model <- sb_fit_reusable_runtime_learner(learner, x, endpoint,
          d$patient_id, "none", phase = phase)
        list(model = model, rng = .Random.seed)
      })
      first <- run(); second <- run()
      expect_equal(first$model, native$model)
      expect_equal(second$model, native$model)
      expect_identical(first$rng, native$rng)
      expect_identical(second$rng, native$rng)
    }, diagnostics = TRUE)
    expect_identical(out$reuse[[paste0(phase, "_learner_hits")]], 1L)
    expect_identical(out$reuse[[paste0(phase, "_learner_misses")]], 1L)
  }
})

check_parametric_learner_reuse <- function(akritas_only = FALSE) {
  dependencies <- if (akritas_only) c("survivalmodels", "distr6") else
    c("flexsurv", "timereg", "ahaz", "coxphf", "ncvreg", "abess")
  for (p in dependencies) {
    skip_if_not_installed(p)
  }
  remaining_binary <- c("catboost_binary", "rlt_reinforcement_binary",
    "adaboost_binary", "logitboost_binary", "wsrf_binary")
  remaining_survival <- c("catboost_cox", "catboost_aft_survival",
    "gbm_cox", "coxboost", "plsrcox", "superpc",
    "mtlr_survival", "pamr_survival",
    "bagged_knn_survival", "bujar_linear_boost",
    "rlt_reinforcement_survival", "rpart_survival", "rsf_survival",
    "ctree_survival", "cforest_survival", "aorsf_survival")
  extra_binary <- c(sb_omics_specialist_binary_learner_ids(), sb_mclust_da_learner_ids(),
    sb_mixomics_da_learner_ids(), sb_hdda_binary_learner_ids())
  ids <- c(remaining_binary, remaining_survival,
    sb_lightgbm_binary_learner_ids(), sb_xgboost_learner_ids(), extra_binary,
    sb_ranger_learner_ids(), sb_survreg_learner_ids(),
    sb_flexsurv_parametric_learner_ids(),
    sb_flexsurvspline_learner_ids(), sb_timereg_survival_learner_ids(),
    sb_ahaz_survival_learner_ids(), sb_coxphf_learner_ids(),
    sb_ncvreg_learner_ids(), sb_abess_learner_ids())
  if (akritas_only) ids <- "akritas_survival"
  features <- c("gene_a", "gene_b", "gene_c")
  calls <- 0L
  native <- sb_fit_runtime_learner
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L
    native(learner, x, endpoint)
  })
  for (id in ids) {
    features <- c("gene_a", "gene_b", "gene_c")
    binary <- id %in% c(remaining_binary, sb_lightgbm_binary_learner_ids(), sb_xgboost_binary_learner_ids(),
      extra_binary, sb_ranger_binary_learner_ids(),
      sb_ncvreg_binary_learner_ids(), sb_abess_binary_learner_ids())
    if (binary) {
      d <- sb_example_binary(180L, 1208L)
      if (id %in% extra_binary) {
        d$gene_d <- withr::with_seed(702L, stats::rnorm(nrow(d)))
        features <- c(features, "gene_d")
      }
      task <- sb_binary_task(d, "patient_id", "outcome", features)
    } else {
      d <- sb_example_survival(180L, 11032L)
      if (id == "aalen_additive_survival") d <- withr::with_seed(9431L, {
        x <- matrix(stats::runif(180L * 3L), ncol = 3L)
        event_time <- stats::rexp(180L, 0.08 + drop(x %*% c(0.08, 0.05, 0.03)))
        censor_time <- stats::rexp(180L, 0.025)
        data.frame(patient_id = paste0("A", seq_len(180L)),
          time = pmax(pmin(event_time, censor_time), 0.01),
          event = as.integer(event_time <= censor_time),
          gene_a = x[, 1L], gene_b = x[, 2L], gene_c = x[, 3L])
      })
      task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
        prediction_horizon = 7)
    }
    pipelines <- list(sb_pipeline("none", id), sb_pipeline("fixed_features", id,
      selector_features = features, feature_set_id = "all"))
    run <- function(enabled, final) withr::with_seed(101L, {
      warnings <- character()
      value <- withCallingHandlers({
        if (final) sb_with_selector_reuse(lapply(pipelines,
          function(p) sb_fit_final(task, p, capture_diagnostics = TRUE)), enabled = enabled)
        else sb_compare_cv(task, pipelines, folds = 3L, seed = 11L, reuse_selectors = enabled)
      }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      })
      list(value = value, warnings = warnings, rng = .Random.seed)
    })
    for (final in c(FALSE, TRUE)) {
      calls <- 0L; plain <- run(FALSE, final); nplain <- calls
      calls <- 0L; shared <- run(TRUE, final)
      comparable <- function(value) {
        if (final && id %in% c("aorsf_survival", "ctree_survival",
            "cforest_survival")) value <- lapply(value, function(fit) {
          # Environment reference encoding can differ; compare the restored
          # native model and, below, exact portable prediction frames.
          fit$learner_state$engine_raw <- unserialize(fit$learner_state$engine_raw)
          fit
        })
        if (id %in% c("catboost_binary", "catboost_cox", "catboost_aft_survival")) {
          clean_cat <- function(x) {
            if (inherits(x, "catboost.Model")) {
              return(list(fields = x[setdiff(names(x), "cpp_obj")],
                parameters = catboost::catboost.get_model_params(x)))
            }
            if (is.list(x)) {
              y <- lapply(seq_along(x), function(i) {
                if (!is.null(names(x)) && names(x)[[i]] == "model_raw" && is.raw(x[[i]])) {
                  clean_cat(unserialize(x[[i]]))
                } else clean_cat(x[[i]])
              })
              attributes(y) <- attributes(x)
              return(y)
            }
            x
          }
          # Native CatBoost model identities are volatile between two fits;
          # compare parameters/audits and exact OOF/portable predictions instead.
          value <- clean_cat(value)
        }
        if (final && id == "rsf_survival") value <- lapply(value, function(fit) {
          clean_rsf <- function(engine) {
            engine[c("ctime.internal", "ctime.external")] <- NULL
            engine$forest$sampsize <- list(body = body(engine$forest$sampsize),
              formals = formals(engine$forest$sampsize))
            engine
          }
          fit$learner_state$engine_raw <- clean_rsf(unserialize(fit$learner_state$engine_raw))
          fit$diagnostic_artifacts$learner$engine <- clean_rsf(fit$diagnostic_artifacts$learner$engine)
          fit
        })
        if (final && id %in% sb_lightgbm_binary_learner_ids()) {
          value <- lapply(value, function(fit) {
            booster <- fit$diagnostic_artifacts$learner$booster
            fields <- intersect(c("best_iter", "best_score", "params", "record_evals",
              "data_processor", "raw"), names(booster))
            fit$diagnostic_artifacts$learner$booster <- list(
              model = booster$save_model_to_string(),
              metadata = setNames(lapply(fields, function(name) booster[[name]]), fields))
            # Internal source-precision audit storage is not a separate model
            # result; the actual public importance tables remain compared below.
            fit$diagnostic_artifacts$learner$reused_native_importance <- NULL
            fit
          })
        }
        if (final && id %in% sb_xgboost_learner_ids()) {
          value <- lapply(value, function(fit) {
            booster <- fit$diagnostic_artifacts$learner$booster
            # Native pointers differ even between two non-shared fits.
            fit$diagnostic_artifacts$learner$booster <- list(
              model = xgboost::xgb.save.raw(booster, raw_format = "json"),
              metadata = if (is.list(booster)) booster[setdiff(names(booster), "handle")] else NULL)
            fit
          })
        }
        value
      }
      expect_equal(comparable(shared$value), comparable(plain$value), info = id)
      expect_identical(shared$warnings, plain$warnings, info = id)
      expect_identical(shared$rng, plain$rng, info = id)
      expect_identical(calls, if (id %in% c("rsf_survival", "adaboost_binary")) nplain else if (final) 1L else 3L, info = id)
      expect_identical(nplain, if (final) 2L else 6L, info = id)
      if (final) {
        expect_identical(vapply(shared$value, `[[`, "", "status"), rep("success", 2L), info = id)
        for (i in 1:2) expect_identical(
          sb_predict_final(unserialize(serialize(shared$value[[i]], NULL)),
            d[c("patient_id", features)], "patient_id"),
          sb_predict_final(plain$value[[i]], d[c("patient_id", features)], "patient_id"), info = id)
      } else {
        expect_identical(unique(shared$value$comparison$status), "success", info = id)
        for (i in 1:2) expect_identical(shared$value$results[[i]]$predictions,
          plain$value$results[[i]]$predictions, info = id)
      }
    }
  }
}

test_that("parametric, penalized and discriminant learners reuse identical analysis fits", {
  check_parametric_learner_reuse()
})

test_that("Akritas reuses identical analysis fits when its dependencies are available", {
  check_parametric_learner_reuse(akritas_only = TRUE)
})

test_that("failed Akritas fits are never reused as successful models", {
  d <- sb_example_survival(48L, 42L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month")
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(...) {
    calls <<- calls + 1L
    stop("synthetic learner failure")
  })
  for (phase in c("cv", "final")) {
    calls <- 0L
    results <- sb_with_selector_reuse(lapply(1:2, function(i) {
      p <- sb_pipeline("none", "akritas_survival")
      if (phase == "final") sb_fit_final(task, p) else
        sb_run_cv(task, p, folds = 3L, seed = 82L)
    }))
    expect_identical(calls, if (phase == "final") 2L else 6L)
    for (result in results) {
      if (phase == "final") expect_identical(result$status, "failed") else {
        expect_identical(result$complete_oof, FALSE)
        expect_true(all(result$runs$status == "failed"))
      }
    }
  }
})

# ---- learner-reuse-smooth ----

test_that("smooth and boosting families preserve CV and final predictions under reuse", {
  for (package in c("mboost", "mgcv", "partykit")) skip_if_not_installed(package)
  ids <- c(sb_mboost_cox_learner_ids(), sb_mboost_aft_learner_ids(),
    sb_mboost_ranking_learner_ids(), sb_mgcv_gam_learner_ids())
  features <- c("gene_a", "gene_b", "gene_c")
  native <- sb_fit_runtime_learner
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L; native(learner, x, endpoint)
  })
  for (id in ids) {
    if (id == "gam_binary_reml") {
      d <- sb_example_binary(120L, 1208L)
      task <- sb_binary_task(d, "patient_id", "outcome", features)
    } else {
      d <- sb_example_survival(120L, 11032L)
      task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
        prediction_horizon = 7)
    }
    pipelines <- list(sb_pipeline("none", id), sb_pipeline("fixed_features", id,
      selector_features = features, feature_set_id = "all"))
    run <- function(enabled, final) withr::with_seed(101L, {
      warnings <- character()
      value <- withCallingHandlers({
        if (final) sb_with_selector_reuse(lapply(pipelines,
          function(p) sb_fit_final(task, p, capture_diagnostics = TRUE)), enabled = enabled)
        else sb_compare_cv(task, pipelines, folds = 3L, seed = 11L, reuse_selectors = enabled)
      }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      })
      list(value = value, warnings = warnings, rng = .Random.seed)
    })
    for (final in c(FALSE, TRUE)) {
      calls <- 0L; plain <- run(FALSE, final); nplain <- calls
      calls <- 0L; shared <- run(TRUE, final)
      # Serializing an isolated mboost environment can change reference encoding,
      # not the model. Compare restored engines rather than their raw bytes.
      comparable <- function(value) {
        if (final) value <- lapply(value, function(fit) {
          if (is.raw(fit$learner_state$engine_raw)) {
            fit$learner_state$engine_raw <- unserialize(fit$learner_state$engine_raw)
          }
          fit
        })
        value
      }
      expect_equal(comparable(shared$value), comparable(plain$value), info = id)
      expect_identical(shared$warnings, plain$warnings, info = id)
      expect_identical(shared$rng, plain$rng, info = id)
      expected_calls <- if (id == "blackboost_cox") nplain else if (final) 1L else 3L
      expect_identical(calls, expected_calls, info = id)
      expect_identical(nplain, if (final) 2L else 6L, info = id)
      if (final) {
        expect_identical(vapply(shared$value, `[[`, "", "status"), rep("success", 2L), info = id)
        for (i in 1:2) expect_identical(
          sb_predict_final(unserialize(serialize(shared$value[[i]], NULL)),
            d[c("patient_id", features)], "patient_id"),
          sb_predict_final(plain$value[[i]], d[c("patient_id", features)], "patient_id"), info = id)
      } else {
        expect_identical(unique(shared$value$comparison$status), "success", info = id)
        for (i in 1:2) expect_identical(shared$value$results[[i]]$predictions,
          plain$value$results[[i]]$predictions, info = id)
      }
    }
  }
})

test_that("native RNG-dependent fits reuse only the same starting random state", {
  skip_if_not_installed("mboost")
  skip_if_not_installed("partykit")
  d <- sb_example_survival(80L, 11032L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "survival", time = d$time, event = d$event)
  learner <- sb_mboost_cox_learner_adapter("blackboost_cox")
  native <- withr::with_seed(55L, {
    model <- sb_fit_runtime_learner(learner, x, endpoint)
    list(model = model, rng = .Random.seed)
  })
  out <- sb_with_selector_reuse({
    run <- function() withr::with_seed(55L, {
      model <- sb_fit_reusable_runtime_learner(learner, x, endpoint, d$patient_id, "none")
      list(model = model, rng = .Random.seed)
    })
    first <- run(); second <- run()
    expect_equal(first$model, native$model)
    expect_equal(second$model, native$model)
    expect_identical(first$rng, native$rng)
    expect_identical(second$rng, native$rng)
  }, diagnostics = TRUE)
  expect_identical(out$reuse$cv_learner_hits, 1L)
  expect_identical(out$reuse$cv_learner_misses, 1L)
})

# ---- rsf-reuse-payload ----

test_that("RSF sampling callback does not retain the native training frame", {
  skip_if_not_installed("randomForestSRC")
  d <- sb_example_survival(80L, 11032L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_additional_survival_learner_adapter("rsf_survival")
  model <- adapter$fit(x, d$time, d$event)
  expect_false(exists("nativeOutput", environment(model$engine$forest$sampsize), inherits = FALSE))
  expect_identical(model$engine$forest$sampsize(80L), 80L * 0.632)
  prepared <- sb_additional_survival_formula_data(x, d$time, d$event,
    random_forest_src = TRUE)
  native <- withr::with_options(list(rf.cores = 1L, mc.cores = 1L),
    randomForestSRC::rfsrc(prepared$formula,
      data = prepared$data, ntree = 1000L, nodesize = 5L, splitrule = "logrank",
      importance = FALSE, proximity = FALSE, forest = TRUE, seed = 1L))
  reference <- model
  reference$engine <- native
  expect_identical(withr::with_seed(22L, adapter$predict(model, x)),
    withr::with_seed(22L, adapter$predict(reference, x)))
  expect_identical(withr::with_seed(23L,
    adapter$predict_survival(model, x, d$time, d$event, x, 7)),
    withr::with_seed(23L,
      adapter$predict_survival(reference, x, d$time, d$event, x, 7)))
})

# ---- survival-prediction-reuse ----

test_that("bagged kNN inverts each bootstrap covariance once per prediction batch", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(72L, 42L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_bagged_knn_survival_learner_adapter("bagged_knn_survival")
  model <- adapter$fit(x, d$time, d$event)
  native_solve <- base::solve.default
  count <- 0L
  local_mocked_bindings(solve.default = function(...) {
    count <<- count + 1L
    native_solve(...)
  }, .package = "base")
  curves <- sb_predict_bagged_knn_survival_curves(model, x[1:4, ])
  expect_identical(count, model$params$num_base_learners)
  expect_identical(nrow(curves$survival), 4L)
})

test_that("Akritas multivariate prediction avoids rebuilding an R6 distribution", {
  skip_if_not_installed("distr6"); skip_if_not_installed("survivalmodels")
  d <- sb_example_survival(72L, 42L); x <- d[c("gene_a", "gene_b", "gene_c")]
  model <- sb_fit_akritas_survival(x, d$time, d$event, "akritas_survival", sb_akritas_params())
  old <- sb_akritas_distribution
  calls <- 0L
  local_mocked_bindings(sb_akritas_distribution = function(basis) {
    calls <<- calls + 1L
    old(basis)
  })
  sb_akritas_curves(model, x[1:4, ])
  expect_identical(calls, 0L)
})

test_that("Akritas state prediction evaluates its risk grid and horizon in one call", {
  skip_if_not_installed("distr6"); skip_if_not_installed("survivalmodels")
  d <- sb_example_survival(72L, 42L); x <- d[c("gene_a", "gene_b")]
  model <- sb_fit_akritas_survival(x, d$time, d$event, "akritas_survival", sb_akritas_params())
  h <- as.numeric(median(d$time))
  state <- sb_extract_akritas_state(model, list(time = d$time, prediction_horizon = h, time_unit = "month"))
  original <- sb_akritas_curves
  calls <- 0L
  local_mocked_bindings(sb_akritas_curves = function(...) {
    calls <<- calls + 1L
    original(...)
  })
  result <- sb_predict_akritas_state(state, x[1:4, ])
  expect_identical(calls, 1L)
  expect_identical(result$survival_probability, as.numeric(original(state, x[1:4, ], h)))
  expect_identical(result$risk, sb_predict_akritas_risk(model, x[1:4, ]))
})

test_that("bagged kNN preserves native curves, ties, one-feature inputs and failures", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(72L, 42L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene A", "gene-B", "gene C")
  x[2L, ] <- x[1L, ]
  adapter <- sb_bagged_knn_survival_learner_adapter("bagged_knn_survival")
  native <- get("predict", asNamespace("bnnSurvival"))
  old_method <- methods::getMethod("predict", "bnnSurvivalEnsemble", where = asNamespace("bnnSurvival"))
  withr::local_seed(19)
  withr::local_options(mc.cores = 7L)
  before <- .Random.seed
  for (p in c(1L, 3L)) {
    train <- x[seq_len(p)]
    model <- adapter$fit(train, d$time, d$event)
    for (rows in list(1L, c(1L, 2L, 1L, 4L, 8L))) {
      query <- train[rows, , drop = FALSE]
      internal <- sb_bagged_knn_internal_data(query)$data
      expected <- sb_with_bagged_knn_single_core(native(model$engine, internal))
      observed <- sb_predict_bagged_knn_survival_curves(model, query)
      expect_identical(observed$survival, as.matrix(bnnSurvival::predictions(expected)))
      expect_identical(observed$time, as.numeric(bnnSurvival::timepoints(expected)))
      state <- adapter$extract_state(model, list(time = d$time,
        prediction_horizon = as.numeric(median(d$time)), time_unit = "month"))
      replay <- adapter$predict_state(unserialize(serialize(state, NULL)), query)
      expect_identical(replay$risk, observed$risk)
      expect_identical(.Random.seed, before)
      expect_identical(getOption("mc.cores"), 7L)
    }
  }
  expect_identical(methods::getMethod("predict", "bnnSurvivalEnsemble", where = asNamespace("bnnSurvival")), old_method)
  singular <- x[c(1L, 1L)]
  names(singular) <- c("a", "duplicate")
  model <- adapter$fit(singular, d$time, d$event)
  query <- singular[1:2, ]; internal <- sb_bagged_knn_internal_data(query)$data
  expected <- tryCatch(sb_with_bagged_knn_single_core(native(model$engine, internal)), error = identity)
  observed <- tryCatch(sb_predict_bagged_knn_survival_curves(model, query), error = identity)
  expect_s3_class(expected, "error")
  expect_identical(class(observed), class(expected))
  expect_identical(conditionMessage(observed), conditionMessage(expected))
})

test_that("Akritas native CDF and batched times retain native boundary behavior", {
  skip_if_not_installed("distr6"); skip_if_not_installed("survivalmodels")
  withr::local_seed(9281)
  for (p in c(1L, 2L, 20L)) {
    basis <- matrix(rnorm(60L * p), ncol = p)
    basis[2L, ] <- basis[1L, ]
    query <- rbind(basis[1:3, , drop = FALSE], rep(-1e6, p), rep(1e6, p))
    expected <- as.numeric(sb_akritas_distribution(basis)$cdf(data = query))
    expect_identical(sb_akritas_empirical_cdf(basis, query), expected)
  }
  d <- sb_example_survival(72L, 42L)
  d$time <- rep(seq_len(12L), length.out = nrow(d))
  x <- d[c("gene_a", "gene_b", "gene_c")]
  before <- .Random.seed
  for (p in c(1L, 3L)) {
    train <- x[seq_len(p)]
    model <- sb_fit_akritas_survival(train, d$time, d$event, "akritas_survival", sb_akritas_params())
    native <- utils::getFromNamespace("C_Akritas", "survivalmodels")
    distribution <- sb_akritas_distribution(unname(as.matrix(train)))
    query <- train[c(1, 2, 7), , drop = FALSE]
    q <- as.numeric(distribution$cdf(data = unname(as.matrix(query))))
    for (h in c(min(d$time) / 2, min(d$time), median(d$time), median(d$time) + 0.25, max(d$time))) {
      state <- sb_extract_akritas_state(model, list(time = d$time,
        prediction_horizon = as.numeric(h), time_unit = "month"))
      expected <- as.numeric(native(model$truth, as.numeric(h), model$grid,
        model$cdf_train, q, model$params$lambda))
      replay <- sb_predict_akritas_state(unserialize(serialize(state, NULL)), query)
      expect_identical(replay$survival_probability, expected)
      expect_identical(replay$risk, sb_predict_akritas_risk(model, query))
    }
  }
  expect_identical(.Random.seed, before)
})
})

# ---- ordinary-native-warning-policy ----
local({
for (id in c("adaboost_binary", "logitboost_binary",
  "abess_binary_gic", "abess_cox_gic")) {
  test_that(paste(id, "preserves native results with ordinary warnings"), {
    package <- switch(id, adaboost_binary = "fastAdaboost",
      logitboost_binary = "caTools", "abess")
    skip_if_not_installed(package)
    survival_task <- identical(id, "abess_cox_gic")
    d <- if (survival_task) make_sb_fixture(n = 48L, seed = 8612L) else {
      make_sb_binary_fixture(n = 48L, seed = 8611L)
    }
    x <- d[c("gene_a", "gene_b", "gene_c")]
    adapter <- if (startsWith(id, "abess")) sb_abess_learner_adapter(id) else {
      sb_classical_binary_boosting_learner_adapter(id, 5L)
    }
    fit <- function() {
      if (survival_task) adapter$fit(x, d$time, d$event) else adapter$fit(x, d$outcome)
    }
    reference <- fit()
    expected <- adapter$predict(reference, x)
    if (startsWith(id, "abess")) {
      native_fit <- abess::abess
      local_mocked_bindings(abess = function(...) {
        warning("informational fit warning")
        native_fit(...)
      }, .package = "abess")
    } else if (id == "adaboost_binary") {
      native_fit <- sb_fastadaboost_fit
      local_mocked_bindings(sb_fastadaboost_fit = function(...) {
        warning("informational fit warning")
        native_fit(...)
      })
    } else {
      native_fit <- sb_catools_logitboost_fit
      local_mocked_bindings(sb_catools_logitboost_fit = function(...) {
        warning("informational fit warning")
        native_fit(...)
      })
    }
    native_predict <- stats::predict
    local_mocked_bindings(predict = function(...) {
      warning("informational prediction warning")
      native_predict(...)
    }, .package = "stats")
    seen <- character()
    observed <- withCallingHandlers({
      model <- fit()
      adapter$predict(model, x)
    }, warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    expect_identical(observed, expected)
    expect_identical(all(c("informational fit warning", "informational prediction warning") %in% seen), TRUE)
  })
}

test_that("ABESS coercion and changed-model warnings remain failures", {
  messages <- c(
    "c.max should be an integer. It is coerced to as.integer(c.max).",
    "x should be a *numeric* matrix/data.frame! The factor value are coerced to as.numeric(x).",
    "y is a binary variable and is not match to family = 'multinomial' or 'ordinal'. We change to family = 'binomial'",
    "The number of the category of y is relative large compare to nvars. The numerical result might be unstable."
  )
  for (message in messages) {
    error <- tryCatch(sb_abess_stop_on_warning({warning(message); 42},
      "ABESS fit", "stormbreaker_abess_fit_warning"), error = identity)
    expect_s3_class(error, "stormbreaker_abess_fit_warning")
  }
})
})
