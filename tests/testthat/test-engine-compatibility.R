# ---- xgboost ----
local({
# ---- xgboost-aft-variants ----

test_that("XGBoost AFT boosters have one explicit native specification", {
  learners <- c(
    "xgboost_aft",
    "xgboost_aft_dart",
    "xgboost_aft_linear"
  )
  expect_identical(sb_xgboost_aft_learner_ids(), learners)

  for (learner in learners) {
    params <- sb_xgboost_params(learner)
    expect_identical(params$nrounds, 100L)
    expect_identical(params$nthread, 1L)
    expect_identical(params$seed, 1L)
    expect_identical(params$objective, "survival:aft")
    expect_identical(params$aft_loss_distribution, "normal")
    expect_identical(params$aft_loss_distribution_scale, 1)
  }
  expect_identical(
    sb_xgboost_booster("xgboost_aft"),
    "gbtree"
  )
  expect_identical(
    sb_xgboost_booster("xgboost_aft_dart"),
    "dart"
  )
  expect_identical(
    sb_xgboost_booster("xgboost_aft_linear"),
    "gblinear"
  )

  registry <- sb_pipeline_registry()
  selectors <- expected_fixed_learner_selectors("survival")
  for (learner in learners) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(unique(rows$task_type), "survival")
    expect_setequal(rows$selector, selectors)
    expect_true(all(grepl("xgboost", rows$optional_dependency, fixed = TRUE)))
  }
})

test_that("XGBoost AFT uses right-censoring bounds and declared controls", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  data <- make_sb_fixture(n = 56L, seed = 8831L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  original_train <- xgboost::xgb.train
  observed <- NULL

  testthat::with_mocked_bindings(
    xgb.train = function(params = list(), data, nrounds, ...) {
      observed <<- list(
        params = params,
        lower = xgboost::getinfo(data, "label_lower_bound"),
        upper = xgboost::getinfo(data, "label_upper_bound")
      )
      original_train(
        params = params,
        data = data,
        nrounds = nrounds,
        ...
      )
    },
    {
      params <- sb_xgboost_params("xgboost_aft")
      params$nrounds <- 4L
      sb_xgboost_learner_adapter(
        "xgboost_aft",
        params
      )$fit(x, data$time, data$event)
    },
    .package = "xgboost"
  )

  expect_equal(observed$lower, data$time, tolerance = 1e-5)
  expect_equal(
    observed$upper[data$event == 1L],
    data$time[data$event == 1L],
    tolerance = 1e-5
  )
  expect_true(all(is.infinite(observed$upper[data$event == 0L])))
  expect_identical(observed$params$objective, "survival:aft")
  expect_identical(observed$params$aft_loss_distribution, "normal")
  expect_identical(observed$params$aft_loss_distribution_scale, 1)
})

test_that("XGBoost AFT risk, native survival, and raw state replay exactly", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  data <- make_sb_fixture(n = 72L, seed = 8832L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene one", "2gene", "gene-three")
  new_x <- x[c(4L, 29L, 61L), , drop = FALSE]
  withr::local_seed(8833L)

  for (learner in sb_xgboost_aft_learner_ids()) {
    params <- sb_xgboost_params(learner)
    params$nrounds <- 12L
    adapter <- sb_xgboost_learner_adapter(learner, params)
    rng_before <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, rng_before)
    repeated <- adapter$fit(x, data$time, data$event)
    expect_identical(
      xgboost::xgb.save.raw(model$booster),
      xgboost::xgb.save.raw(repeated$booster)
    )

    native_input <- xgboost::xgb.DMatrix(as.matrix(new_x), nthread = 1L)
    native_margin <- as.numeric(stats::predict(
      model$booster,
      native_input,
      outputmargin = TRUE
    ))
    native_time <- as.numeric(stats::predict(
      model$booster,
      native_input
    ))
    risk <- adapter$predict(model, new_x)
    probability <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      new_x,
      prediction_horizon = 7
    )
    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = 7,
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, new_x)

    expect_equal(native_time, exp(native_margin), tolerance = 1e-6)
    expect_equal(risk, -native_margin, tolerance = 0)
    expect_equal(
      probability$probability,
      stats::pnorm(native_margin - log(7)),
      tolerance = 0
    )
    expect_identical(probability$source, "xgboost_aft_normal_native")
    expect_identical(probability$ties, "not_applicable")
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      probability$probability,
      tolerance = 0
    )
    expect_identical(state$risk_reference, "xgboost_aft_normal_native")
    expect_identical(
      state$risk_formula,
      "negative_xgboost_aft_log_time_margin"
    )
    expect_type(state$model_raw, "raw")
    expect_false("booster" %in% names(state))
    expect_false(any(c("time", "event", "outcome") %in% names(state)))
  }
})

test_that("XGBoost AFT supports one feature and monotone horizon survival", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  data <- make_sb_fixture(n = 64L, seed = 8834L)
  x <- data["gene_a"]
  names(x) <- "1 gene-A"
  params <- sb_xgboost_params("xgboost_aft")
  params$nrounds <- 10L
  adapter <- sb_xgboost_learner_adapter(
    "xgboost_aft",
    params
  )
  model <- adapter$fit(x, data$time, data$event)
  horizons <- c(4, 7, 10)
  probability <- vapply(
    horizons,
    function(horizon) {
      adapter$predict_survival(
        model,
        x,
        data$time,
        data$event,
        x[1:6, , drop = FALSE],
        prediction_horizon = horizon
      )$probability
    },
    numeric(6L)
  )

  expect_true(all(probability >= 0 & probability <= 1))
  expect_true(all(probability[, 1L] >= probability[, 2L]))
  expect_true(all(probability[, 2L] >= probability[, 3L]))
})

test_that("XGBoost AFT boosters execute through CV and final prediction", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  data <- make_sb_fixture(n = 72L, seed = 8835L)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  for (learner in sb_xgboost_aft_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 2L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 8836L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:5, c("patient_id", "gene_c", "gene_a", "gene_b")],
      "patient_id"
    )

    expect_identical(pipeline$risk_type, "aft_log_time")
    expect_identical(
      pipeline$risk_reference,
      "xgboost_aft_normal_native"
    )
    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(
      comparison$complete_survival_probability_oof,
      TRUE
    )
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(
      prediction$predictions$survival_probability >= 0 &
        prediction$predictions$survival_probability <= 1
    ))
  }
})

# ---- xgboost-booster-variants ----

test_that("XGBoost DART and linear boosters are distinct fixed learners", {
  binary_ids <- c(
    "xgboost_binary_dart",
    "xgboost_binary_linear"
  )
  cox_ids <- c(
    "xgboost_cox_dart",
    "xgboost_cox_linear"
  )

  expect_setequal(sb_xgboost_binary_learner_ids(), c(
    "xgboost_binary", binary_ids
  ))
  expect_setequal(sb_xgboost_cox_learner_ids(), c(
    "xgboost_cox", cox_ids
  ))

  dart <- sb_xgboost_params("xgboost_binary_dart")
  expect_identical(dart, list(
    nrounds = 100L,
    booster = "dart",
    eta = 0.1,
    max_depth = 3L,
    min_child_weight = 1,
    subsample = 1,
    colsample_bytree = 1,
    sample_type = "uniform",
    normalize_type = "tree",
    rate_drop = 0.1,
    skip_drop = 0,
    one_drop = FALSE,
    nthread = 1L,
    seed = 1L,
    objective = "binary:logistic"
  ))
  expect_identical(
    sb_xgboost_params("xgboost_cox_dart")$objective,
    "survival:cox"
  )

  linear <- sb_xgboost_params("xgboost_binary_linear")
  expect_identical(linear, list(
    nrounds = 100L,
    booster = "gblinear",
    eta = 0.1,
    lambda = 1,
    alpha = 0,
    updater = "coord_descent",
    feature_selector = "cyclic",
    nthread = 1L,
    seed = 1L,
    objective = "binary:logistic"
  ))
  expect_identical(
    sb_xgboost_params("xgboost_cox_linear")$objective,
    "survival:cox"
  )

  registry <- sb_pipeline_registry()
  expected_selectors <- list(
    binary = expected_fixed_learner_selectors("binary"),
    survival = expected_fixed_learner_selectors("survival")
  )
  for (learner in c(binary_ids, cox_ids)) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    task <- if (learner %in% binary_ids) "binary" else "survival"
    expect_identical(unique(rows$task_type), task)
    expect_setequal(rows$selector, expected_selectors[[task]])
    expect_true(all(grepl("xgboost", rows$optional_dependency, fixed = TRUE)))
  }
})

test_that("XGBoost booster variants preserve native predictions and state", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  binary <- make_sb_binary_fixture(n = 64L, seed = 8711L)
  binary_x <- binary[c("gene_a", "gene_b", "gene_c")]
  names(binary_x) <- c("gene A", "gene-B", "gene_c")
  binary_new <- binary_x[c(4L, 27L, 51L), , drop = FALSE]
  withr::local_seed(8710L)

  for (learner in c(
    "xgboost_binary_dart",
    "xgboost_binary_linear"
  )) {
    params <- sb_xgboost_params(learner)
    params$nrounds <- 12L
    adapter <- sb_xgboost_learner_adapter(learner, params)
    rng_before <- .Random.seed
    model <- adapter$fit(binary_x, binary$outcome)
    expect_identical(.Random.seed, rng_before)
    repeated <- adapter$fit(binary_x, binary$outcome)
    expect_identical(
      xgboost::xgb.save.raw(model$booster),
      xgboost::xgb.save.raw(repeated$booster)
    )

    probability <- adapter$predict(model, binary_new)
    native <- as.numeric(stats::predict(
      model$booster,
      as.matrix(binary_new)
    ))
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, binary_new)$probability
    config <- xgboost::xgb.config(model$booster)

    expect_equal(probability, native, tolerance = 0)
    expect_equal(replay, probability, tolerance = 0)
    expect_true(all(probability >= 0 & probability <= 1))
    if (is.list(config)) {
      expect_identical(config$learner$learner_train_param$booster, params$booster)
    } else {
      expect_true(grepl(
        paste0('"booster":"', params$booster, '"'),
        config,
        fixed = TRUE
      ))
    }
    expect_type(state$model_raw, "raw")
    expect_false("booster" %in% names(state))
  }

  survival <- make_sb_fixture(n = 70L, seed = 8712L)
  survival_x <- survival[c("gene_a", "gene_b", "gene_c")]
  survival_new <- survival_x[c(6L, 31L, 58L), , drop = FALSE]
  for (learner in c(
    "xgboost_cox_dart",
    "xgboost_cox_linear"
  )) {
    params <- sb_xgboost_params(learner)
    params$nrounds <- 12L
    adapter <- sb_xgboost_learner_adapter(learner, params)
    model <- adapter$fit(survival_x, survival$time, survival$event)
    risk <- adapter$predict(model, survival_new)
    native <- as.numeric(stats::predict(
      model$booster,
      as.matrix(survival_new),
      outputmargin = TRUE
    ))
    probability <- adapter$predict_survival(
      model,
      survival_x,
      survival$time,
      survival$event,
      survival_new,
      prediction_horizon = 7
    )
    state <- adapter$extract_state(model, list(
      x = survival_x,
      time = survival$time,
      event = survival$event,
      prediction_horizon = 7,
      time_unit = "month"
    ))
    replay <- adapter$predict_state(
      unserialize(serialize(state, NULL)),
      survival_new
    )

    expect_equal(risk, native, tolerance = 0)
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      probability$probability,
      tolerance = 0
    )
    expect_true(all(
      probability$probability >= 0 & probability$probability <= 1
    ))
    expect_identical(probability$source, "xgboost_cox_breslow_calibration")
  }
})

test_that("XGBoost forwards the declared DART and linear controls", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  data <- make_sb_binary_fixture(n = 40L, seed = 8713L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  original_train <- xgboost::xgb.train
  observed <- list()

  testthat::with_mocked_bindings(
    xgb.train = function(params = list(), data, nrounds, ...) {
      observed[[length(observed) + 1L]] <<- params
      original_train(
        params = params,
        data = data,
        nrounds = nrounds,
        ...
      )
    },
    {
      for (learner in c(
        "xgboost_binary_dart",
        "xgboost_binary_linear"
      )) {
        params <- sb_xgboost_params(learner)
        params$nrounds <- 3L
        sb_xgboost_learner_adapter(learner, params)$fit(
          x,
          data$outcome
        )
      }
    },
    .package = "xgboost"
  )

  expect_identical(observed[[1L]], list(
    objective = "binary:logistic",
    booster = "dart",
    eta = 0.1,
    max_depth = 3L,
    min_child_weight = 1,
    subsample = 1,
    colsample_bytree = 1,
    nthread = 1L,
    sample_type = "uniform",
    normalize_type = "tree",
    rate_drop = 0.1,
    skip_drop = 0,
    one_drop = FALSE
  ))
  expect_identical(observed[[2L]], list(
    objective = "binary:logistic",
    booster = "gblinear",
    eta = 0.1,
    lambda = 1,
    alpha = 0,
    updater = "coord_descent",
    feature_selector = "cyclic",
    nthread = 1L
  ))
})

test_that("XGBoost booster variants execute through CV and final prediction", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  binary <- make_sb_binary_fixture(n = 60L, seed = 8721L)
  binary_task <- sb_binary_task(
    binary,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  for (learner in c(
    "xgboost_binary_dart",
    "xgboost_binary_linear"
  )) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 2L
    )
    comparison <- sb_run_cv(binary_task, pipeline, folds = 3L, seed = 8722L)
    fit <- sb_fit_final(binary_task, pipeline)
    prediction <- sb_predict_final(
      fit,
      binary[1:5, c("patient_id", "gene_c", "gene_a", "gene_b")],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$probability, 5L)
  }

  survival <- make_sb_fixture(n = 66L, seed = 8723L)
  survival_task <- sb_surv_task(
    survival,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  for (learner in c(
    "xgboost_cox_dart",
    "xgboost_cox_linear"
  )) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 2L
    )
    comparison <- sb_run_cv(
      survival_task,
      pipeline,
      folds = 3L,
      seed = 8724L
    )
    fit <- sb_fit_final(survival_task, pipeline)
    prediction <- sb_predict_final(
      fit,
      survival[1:5, c("patient_id", "gene_b", "gene_c", "gene_a")],
      "patient_id"
    )

    expect_identical(pipeline$risk_type, "cox_raw_margin")
    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(
      comparison$complete_survival_probability_oof,
      TRUE
    )
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
  }
})

# ---- xgboost-fixed-adapters ----

test_that("fixed XGBoost specifications are explicit", {
  binary <- sb_xgboost_params("xgboost_binary")
  survival <- sb_xgboost_params("xgboost_cox")

  expect_identical(binary$nrounds, 100L)
  expect_identical(binary$eta, 0.1)
  expect_identical(binary$max_depth, 3L)
  expect_identical(binary$min_child_weight, 1)
  expect_identical(binary$subsample, 1)
  expect_identical(binary$colsample_bytree, 1)
  expect_identical(binary$nthread, 1L)
  expect_identical(binary$seed, 1L)
  expect_identical(binary$objective, "binary:logistic")
  expect_identical(survival$objective, "survival:cox")
})

test_that("fixed XGBoost uses one thread for matrix construction", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  data <- make_sb_binary_fixture(n = 36L, seed = 7900L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  params <- sb_xgboost_params("xgboost_binary")
  params$nrounds <- 5L
  adapter <- sb_xgboost_learner_adapter(
    "xgboost_binary",
    params
  )
  original_dmatrix <- xgboost::xgb.DMatrix
  observed_threads <- integer()
  observed_matrix_input <- logical()

  testthat::with_mocked_bindings(
    xgb.DMatrix = function(data, nthread = NULL, ...) {
      observed_threads <<- c(
        observed_threads,
        if (is.null(nthread)) NA_integer_ else as.integer(nthread)
      )
      observed_matrix_input <<- c(observed_matrix_input, is.matrix(data))
      original_dmatrix(
        data = data,
        nthread = nthread,
        ...
      )
    },
    {
      model <- adapter$fit(x, data$outcome)
      adapter$predict(model, x[1:3, , drop = FALSE])
    },
    .package = "xgboost"
  )

  expect_identical(
    observed_threads[observed_matrix_input],
    c(1L, 1L)
  )
})

test_that("fixed XGBoost adapters preserve native predictions and state", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  binary <- make_sb_binary_fixture(n = 54L, seed = 7901L)
  binary_x <- binary[c("gene_a", "gene_b", "gene_c")]
  binary_params <- sb_xgboost_params("xgboost_binary")
  binary_params$nrounds <- 12L
  binary_adapter <- sb_xgboost_learner_adapter(
    "xgboost_binary",
    binary_params
  )
  binary_model <- binary_adapter$fit(binary_x, binary$outcome)
  binary_new <- binary_x[c(4L, 21L, 43L), , drop = FALSE]
  binary_probability <- binary_adapter$predict(binary_model, binary_new)
  binary_state <- binary_adapter$extract_state(binary_model)
  binary_replay <- binary_adapter$predict_state(
    unserialize(serialize(binary_state, NULL)),
    binary_new
  )

  expect_equal(
    binary_probability,
    as.numeric(stats::predict(binary_model$booster, as.matrix(binary_new))),
    tolerance = 0
  )
  expect_equal(binary_replay$probability, binary_probability, tolerance = 0)
  expect_type(binary_state$model_raw, "raw")
  expect_false("booster" %in% names(binary_state))

  survival <- make_sb_fixture(n = 60L, seed = 7902L)
  survival_x <- survival[c("gene_a", "gene_b", "gene_c")]
  survival_params <- sb_xgboost_params("xgboost_cox")
  survival_params$nrounds <- 12L
  survival_adapter <- sb_xgboost_learner_adapter(
    "xgboost_cox",
    survival_params
  )
  survival_model <- survival_adapter$fit(
    survival_x,
    survival$time,
    survival$event
  )
  survival_new <- survival_x[c(3L, 28L, 51L), , drop = FALSE]
  risk <- survival_adapter$predict(survival_model, survival_new)
  probability <- survival_adapter$predict_survival(
    survival_model,
    survival_x,
    survival$time,
    survival$event,
    survival_new,
    prediction_horizon = 7
  )
  survival_state <- survival_adapter$extract_state(survival_model, list(
    x = survival_x,
    time = survival$time,
    event = survival$event,
    prediction_horizon = 7,
    time_unit = "month"
  ))
  survival_replay <- survival_adapter$predict_state(
    unserialize(serialize(survival_state, NULL)),
    survival_new
  )

  native_margin <- as.numeric(stats::predict(
    survival_model$booster,
    as.matrix(survival_new),
    outputmargin = TRUE
  ))
  expect_equal(risk, native_margin, tolerance = 0)
  expect_equal(survival_replay$risk, risk, tolerance = 0)
  expect_equal(
    survival_replay$survival_probability,
    probability$probability,
    tolerance = 0
  )
  expect_identical(probability$source, "xgboost_cox_breslow_calibration")
  expect_identical(probability$ties, "breslow")
  expect_true(all(probability$probability >= 0 & probability$probability <= 1))
})

test_that("XGBoost routes run binary and survival workflows", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  binary <- make_sb_binary_fixture(n = 60L, seed = 7911L)
  binary_task <- sb_binary_task(
    binary,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  binary_pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "xgboost_binary",
    selector_n_features = 2L
  )
  binary_result <- sb_run_cv(
    binary_task,
    binary_pipeline,
    folds = 3L,
    seed = 7912L
  )
  binary_fit <- sb_fit_final(binary_task, binary_pipeline)
  binary_prediction <- sb_predict_final(
    binary_fit,
    binary[c(5L, 24L, 47L), c("patient_id", "gene_c", "gene_a", "gene_b")],
    "patient_id"
  )

  expect_identical(binary_result$status, "success")
  expect_true(binary_result$complete_oof)
  expect_identical(binary_fit$status, "success")
  expect_identical(binary_prediction$status, "success")
  expect_true(all(binary_prediction$predictions$probability >= 0 &
    binary_prediction$predictions$probability <= 1))

  survival <- make_sb_fixture(n = 60L, seed = 7913L)
  survival_task <- sb_surv_task(
    survival,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  survival_pipeline <- sb_pipeline(learner = "xgboost_cox")
  survival_result <- sb_run_cv(
    survival_task,
    survival_pipeline,
    folds = 3L,
    seed = 7914L
  )
  survival_fit <- sb_fit_final(survival_task, survival_pipeline)
  survival_prediction <- sb_predict_final(
    survival_fit,
    survival[c(7L, 30L, 55L), c("patient_id", "gene_b", "gene_c", "gene_a")],
    "patient_id"
  )

  expect_identical(survival_pipeline$risk_type, "cox_raw_margin")
  expect_identical(survival_pipeline$risk_reference, "xgboost_native")
  expect_identical(survival_result$status, "success")
  expect_true(survival_result$complete_oof)
  expect_true(survival_result$complete_survival_probability_oof)
  expect_identical(survival_fit$status, "success")
  expect_identical(survival_prediction$status, "success")
})
})

# ---- lightgbm ----
local({
# ---- lightgbm-binary-fixed-adapter ----

test_that("fixed LightGBM model text reconstructs class-one probabilities", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  x <- data.frame(
    gene_a = c(-2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
    gene_b = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    unused = 1,
    check.names = FALSE
  )
  outcome <- c(0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 1L)
  new_x <- data.frame(
    gene_a = c(-1.5, 0.75, 3.25),
    gene_b = c(1, 0, 1),
    unused = 1,
    check.names = FALSE
  )
  adapter <- sb_lightgbm_binary_learner_adapter(
    nrounds = 12L,
    learning_rate = 0.1,
    num_leaves = 7L,
    min_data_in_leaf = 2L,
    seed = 101L
  )
  model <- adapter$fit(x, outcome)
  observed <- adapter$predict(model, new_x)
  audit <- adapter$extract_tree_model(model)
  engine_params <- model$booster$params
  expected_engine_params <- list(
    learning_rate = 0.1,
    num_leaves = 7L,
    min_data_in_leaf = 2L,
    seed = 101L,
    objective = "binary",
    boosting = "gbdt",
    device_type = "cpu",
    num_threads = 1L,
    deterministic = TRUE,
    force_col_wise = TRUE,
    feature_fraction = 1,
    bagging_fraction = 1,
    bagging_freq = 0L
  )
  restored <- lightgbm::lgb.load(model_str = audit$model_string)
  rebuilt <- as.numeric(stats::predict(
    restored,
    as.matrix(new_x),
    type = "response",
    num_iteration = audit$nrounds,
    params = list(num_threads = 1L)
  ))

  expect_s3_class(model, "stormbreaker_lightgbm_binary")
  expect_identical(
    engine_params[names(expected_engine_params)],
    expected_engine_params
  )
  expect_identical(engine_params$num_iterations, 12L)
  expect_identical(as.integer(model$booster$num_iter()), 12L)
  expect_identical(as.integer(model$booster$num_trees()), 12L)
  expect_identical(model$booster$best_iter, -1L)
  expect_identical(model$booster$record_evals, list())
  expect_equal(observed, rebuilt, tolerance = 0)
  expect_true(all(is.finite(observed)))
  expect_true(all(observed >= 0 & observed <= 1))
  expect_identical(audit$nrounds, 12L)
  expect_identical(audit$features, names(x))
  expect_identical(audit$params, adapter$params)
  expect_type(audit$model_string, "character")
  expect_length(audit$model_string, 1L)
  expect_named(
    audit$feature_importance,
    c("feature", "gain", "cover", "frequency"),
    ignore.order = FALSE
  )
  expect_identical(audit$feature_importance$feature, names(x))
  expect_identical(
    audit$feature_importance$frequency[audit$feature_importance$feature ==
      "unused"],
    0L
  )
})

test_that("fixed LightGBM predictions do not depend on assessment batching", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  d <- make_sb_binary_fixture(n = 50L, seed = 2801L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_lightgbm_binary_learner_adapter(
    nrounds = 10L,
    learning_rate = 0.1,
    num_leaves = 7L,
    min_data_in_leaf = 3L,
    seed = 2901L
  )
  model <- adapter$fit(x[1:40, , drop = FALSE], d$outcome[1:40])
  new_x <- x[41:50, , drop = FALSE]

  batch <- adapter$predict(model, new_x)
  one_at_a_time <- vapply(seq_len(nrow(new_x)), function(row) {
    adapter$predict(model, new_x[row, , drop = FALSE])
  }, numeric(1))
  order <- c(10L, 2L, 8L, 4L, 6L, 1L, 9L, 3L, 7L, 5L)
  reordered <- adapter$predict(model, new_x[order, , drop = FALSE])
  restored_order <- reordered[match(seq_along(order), order)]

  expect_equal(batch, one_at_a_time, tolerance = 0)
  expect_equal(batch, restored_order, tolerance = 0)
})

test_that("fixed LightGBM is deterministic and locks predictor schema", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  d <- make_sb_binary_fixture(n = 40L, seed = 3001L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_lightgbm_binary_learner_adapter(
    nrounds = 8L,
    learning_rate = 0.1,
    num_leaves = 7L,
    min_data_in_leaf = 3L,
    seed = 3101L
  )
  first <- adapter$fit(x, d$outcome)
  second <- adapter$fit(x, d$outcome)

  expect_identical(
    adapter$extract_tree_model(first)$model_string,
    adapter$extract_tree_model(second)$model_string
  )
  expect_equal(
    adapter$predict(first, x),
    adapter$predict(second, x),
    tolerance = 0
  )

  schema_error <- capture_sb_error(
    adapter$predict(first, x[c("gene_b", "gene_a", "gene_c")])
  )
  expect_s3_class(schema_error, "stormbreaker_learner_schema_mismatch")
})

test_that("fixed LightGBM dependency and class failures remain explicit", {
  x <- data.frame(
    gene_a = seq_len(8),
    gene_b = rev(seq_len(8))
  )
  adapter <- sb_lightgbm_binary_learner_adapter(
    nrounds = 5L,
    learning_rate = 0.1,
    num_leaves = 3L,
    min_data_in_leaf = 1L,
    seed = 1L
  )

  local_mocked_bindings(sb_lightgbm_available = function() FALSE)
  missing_error <- capture_sb_error(
    adapter$fit(x, rep(c(0L, 1L), each = 4L))
  )
  expect_s3_class(missing_error, "stormbreaker_missing_lightgbm")

  local_mocked_bindings(sb_lightgbm_available = function() TRUE)
  class_error <- capture_sb_error(adapter$fit(x, rep(0L, 8L)))
  expect_s3_class(class_error, "stormbreaker_lightgbm_single_class")
})

# ---- lightgbm-bioconductor ----

test_that("LightGBM logloss aliases survive Bioconductor list formatting", {
  skip_if_not_installed("lightgbm")
  skip_if_not_installed("BiocGenerics")
  requireNamespace("BiocGenerics", quietly = TRUE)
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  methods <- c("lightgbm", "lightgbm_dart", "lightgbm_goss", "lightgbm_rf")
  for (method in methods) {
    fit <- sb_fit_multiclass(task, method, 3, 1)
    x <- as.matrix(d[task$predictors])
    colnames(x) <- fit$feature_codes
    params <- fit$parameters[setdiff(names(fit$parameters), c("nrounds", "metric"))]
    # The objective's native default is multi_logloss, independently of aliases.
    reference <- lightgbm::lgb.train(params,
      lightgbm::lgb.Dataset(x, label = as.integer(factor(d$Species, levels = task$classes)) - 1L),
      nrounds = fit$parameters$nrounds, verbose = -1L)
    expect_match(fit$engine, "\\[metric: multi_logloss\\]")
    expect_equal(unname(sb_predict_multiclass(fit, d)), unname(predict(reference, x)), tolerance = 0)
  }
  x <- d[1:100, 1:4]
  outcome <- rep(0:1, each = 50)
  for (id in sb_lightgbm_binary_learner_ids()) {
    params <- sb_lightgbm_binary_params(learner = id)
    fit <- sb_fit_lightgbm_binary(x, outcome, params, id)
    reference <- lightgbm::lgb.train(params[setdiff(names(params), c("nrounds", "metric"))],
      lightgbm::lgb.Dataset(as.matrix(x), label = outcome), nrounds = params$nrounds, verbose = -1L)
    expect_match(fit$booster$save_model_to_string(), "\\[metric: binary_logloss\\]")
    expect_equal(predict(fit$booster, as.matrix(x)), predict(reference, as.matrix(x)), tolerance = 0)
  }
})

# ---- lightgbm-training-regimes ----

test_that("LightGBM DART, GOSS, and RF are distinct training regimes", {
  learners <- c(
    "lightgbm_binary_dart",
    "lightgbm_binary_goss",
    "lightgbm_binary_rf"
  )
  expect_setequal(
    sb_lightgbm_binary_learner_ids(),
    c("lightgbm_binary", learners)
  )

  dart <- sb_lightgbm_binary_params(
    learner = "lightgbm_binary_dart"
  )
  expect_identical(dart$boosting, "dart")
  expect_identical(dart$drop_rate, 0.1)
  expect_identical(dart$max_drop, 50L)
  expect_identical(dart$skip_drop, 0)
  expect_identical(dart$uniform_drop, FALSE)
  expect_identical(dart$xgboost_dart_mode, FALSE)

  goss <- sb_lightgbm_binary_params(
    learner = "lightgbm_binary_goss"
  )
  expect_identical(goss$boosting, "gbdt")
  expect_identical(goss$data_sample_strategy, "goss")
  expect_identical(goss$top_rate, 0.2)
  expect_identical(goss$other_rate, 0.1)
  expect_identical(goss$bagging_fraction, 1)
  expect_identical(goss$bagging_freq, 0L)

  forest <- sb_lightgbm_binary_params(
    learner = "lightgbm_binary_rf"
  )
  expect_identical(forest$boosting, "rf")
  expect_identical(forest$feature_fraction, 0.8)
  expect_identical(forest$bagging_fraction, 0.8)
  expect_identical(forest$bagging_freq, 1L)

  registry <- sb_pipeline_registry()
  expected_selectors <- expected_fixed_learner_selectors("binary")
  for (learner in learners) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(unique(rows$task_type), "binary")
    expect_setequal(rows$selector, expected_selectors)
    expect_true(all(grepl(
      "lightgbm", rows$optional_dependency,
      fixed = TRUE
    )))
    expect_identical(unique(rows$tuning_semantics), "none")

    pipeline <- sb_pipeline(
      learner = learner,
      lightgbm_nrounds = 12L,
      lightgbm_learning_rate = 0.05,
      lightgbm_num_leaves = 7L,
      lightgbm_min_data_in_leaf = 3L,
      lightgbm_seed = 29L
    )
    expect_identical(pipeline$learner, learner)
    expect_identical(pipeline$learner_params$nrounds, 12L)
    expect_identical(pipeline$learner_params$seed, 29L)
    expect_match(pipeline$id, paste0("none__", learner), fixed = TRUE)
  }
})

test_that("new LightGBM regimes preserve native probabilities and state", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  data <- make_sb_binary_fixture(n = 72L, seed = 8911L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene A", "gene-B", "gene_c")
  assessment <- x[c(5L, 21L, 46L, 67L), , drop = FALSE]
  learners <- c(
    "lightgbm_binary_dart",
    "lightgbm_binary_goss",
    "lightgbm_binary_rf"
  )
  withr::local_seed(8912L)

  for (learner in learners) {
    adapter <- sb_lightgbm_binary_learner_adapter(
      nrounds = 12L,
      learning_rate = 0.1,
      num_leaves = 7L,
      min_data_in_leaf = 3L,
      seed = 1L,
      learner = learner
    )
    rng_before <- .Random.seed
    model <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, rng_before)
    repeated <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, rng_before)
    expect_identical(
      model$booster$save_model_to_string(num_iteration = 12L),
      repeated$booster$save_model_to_string(num_iteration = 12L)
    )

    observed <- adapter$predict(model, assessment)
    native <- as.numeric(stats::predict(
      model$booster,
      as.matrix(assessment),
      type = "response",
      num_iteration = 12L,
      params = list(num_threads = 1L)
    ))
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, assessment)$probability

    expect_equal(observed, native, tolerance = 0)
    expect_equal(replay, observed, tolerance = 0)
    expect_identical(state$learner, learner)
    expect_type(state$model_string, "character")
    expect_identical("booster" %in% names(state), FALSE)
    expect_true(all(observed >= 0 & observed <= 1))
  }
})

test_that("new LightGBM regimes execute through fold-local CV and final fit", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  data <- make_sb_binary_fixture(n = 72L, seed = 8921L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  learners <- c(
    "lightgbm_binary_dart",
    "lightgbm_binary_goss",
    "lightgbm_binary_rf"
  )

  for (learner in learners) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 2L,
      lightgbm_nrounds = 12L,
      lightgbm_num_leaves = 7L,
      lightgbm_min_data_in_leaf = 3L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 8922L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:6, c("patient_id", "gene_c", "gene_a", "gene_b")],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(sum(comparison$runs$status == "success"), 3L)
    expect_identical(
      comparison$provenance$lightgbm_version,
      as.character(utils::packageVersion("lightgbm"))
    )
    expect_identical(fit$status, "success")
    expect_identical(fit$learner_state$learner, learner)
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$probability, 6L)
  }
})
})

# ---- fda-stable-probability ----
local({
test_that("FDA probabilities match native posteriors on ordinary observations", {
  skip_if_not_installed("mda")
  d <- sb_example_binary(80, seed = 9301)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  for (id in c("fda_binary", "pda_binary")) {
    adapter <- sb_classical_binary_learner_adapter(id)
    fit <- adapter$fit(x, d$outcome)
    expected <- unname(predict(fit$engine, newdata = x, type = "posterior")[, "1"])
    expect_equal(adapter$predict(fit, x), expected, tolerance = 1e-12)
  }
})

test_that("FDA computes finite posteriors when native exponentials underflow", {
  skip_if_not_installed("mda")
  d <- sb_example_binary(80, seed = 9301)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  extreme <- x[1:3, , drop = FALSE] * 1000
  for (id in c("fda_binary", "pda_binary")) {
    adapter <- sb_classical_binary_learner_adapter(id)
    fit <- adapter$fit(x, d$outcome)
    native <- predict(fit$engine, newdata = extreme, type = "posterior")
    distances <- predict(fit$engine, newdata = extreme, type = "distances")
    expect_true(any(!is.finite(native)))
    expected <- plogis(-0.5 * (distances[, "1"] - distances[, "0"]) +
      log(fit$engine$prior["1"] / fit$engine$prior["0"]))
    probability <- adapter$predict(fit, extreme)
    expect_equal(probability, unname(expected), tolerance = 1e-12)
    expect_true(all(is.finite(probability) & probability >= 0 & probability <= 1))
    restored <- unserialize(serialize(adapter$extract_state(fit), NULL))
    expect_equal(adapter$predict_state(restored, extreme)$probability, probability)
    expect_equal(adapter$predict(fit, extreme[3:1, , drop = FALSE]), rev(probability))
  }
})
})

# ---- selector-reuse-core ----
local({
# ---- selector-budget-reuse ----

test_that("explicit Boruta budgets retain defaults and survive contracts", {
  for (id in sb_boruta_selector_ids()) {
    learner <- if (grepl("survival", id)) "coxph" else "logistic"
    old <- sb_pipeline(id, learner)
    expect_identical(old$id, paste(id, learner, sep = "__"))
    p <- sb_pipeline(id, learner, selector_options = list(num_trees = 50, maximum_runs = 20))
    expect_identical(p$selector_params$num_trees, 50L)
    expect_identical(p$selector_params$maximum_runs, 20L)
    expect_identical(sb_validate_pipeline_contract(p, p$task_type), p)
    expect_false(identical(p$id, old$id))
  }
  expect_error(sb_pipeline("boruta_survival_roughfix", "coxph",
    selector_options = list(num_threads = 8)), class = "stormbreaker_bad_selector_options")
})
test_that("Boruta budget reaches the native call", {
  skip_if_not_installed("Boruta")
  seen <- NULL
  local_mocked_bindings(sb_call_boruta = function(x, response, params) {
    seen <<- params
    stop("observed native request")
  })
  d <- sb_example_survival(60, 11)
  p <- sb_pipeline("boruta_survival_roughfix", "coxph",
    selector_options = list(num_trees = 50L, maximum_runs = 20L))
  a <- sb_resolve_selector_adapter(p$selector, p$selector_params)
  expect_error(a$fit(d[c("gene_a", "gene_b")], d$time, d$event), "observed native request")
  expect_identical(seen$num_trees, 50L)
  expect_identical(seen$maximum_runs, 20L)
  expect_identical(seen$num_threads, 1L)
  captured <- sb_captured_selector_adapter(p$selector, p$selector_params)
  expect_error(captured$fit(d[c("gene_a", "gene_b")], d$time, d$event), "observed native request")
  expect_identical(seen$num_trees, 50L)
  default_captured <- sb_captured_selector_adapter(p$selector, NULL)
  expect_error(default_captured$fit(d[c("gene_a", "gene_b")], d$time, d$event), "observed native request")
  expect_identical(seen$num_trees, 500L)
})

test_that("retired PST selectors are rejected rather than substituted", {
  for (task in c("binary", "survival")) {
    id <- if (task == "binary") "vsolassobag_binary_pst" else "vsolassobag_cox_pst"
    learner <- if (task == "binary") "logistic" else "coxph"
    expect_false(id %in% sb_pipeline_registry(task_type = task)$selector)
    expect_false(id %in% unname(sb_doctor_selector_menu(task)))
    expect_error(sb_pipeline(id, learner), class = "stormbreaker_unsupported_selector")
    expect_error(sb_doctor_selector("VSOLassoBag_PST", task), "Unknown selector")
  }
})

# ---- selector-reuse-chain-rfe ----

test_that("SVM-RFE shares the same learned prefix while retaining inner refits", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(120L, 1208L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  pipelines <- lapply(c("logistic", "gaussian_naive_bayes"), function(learner) {
    sb_new_chain_pipeline(sb_pipeline("none", learner), list(
      sb_pipeline("univ_binary_auc_top_k", learner, selector_n_features = 2L),
      sb_pipeline("svm_rfe_binary", learner, selector_options = list(inner_folds = 3L))))
  })
  native <- sb_fit_svm_rfe
  calls <- 0L
  local_mocked_bindings(sb_fit_svm_rfe = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L, reuse_selectors = FALSE)
  expect_identical(calls, 6L)
  calls <- 0L
  shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L)
  expect_identical(calls, 3L)
  expect_identical(shared, plain)
  expect_identical(unique(shared$comparison$status), "success")

  # Changing an upstream input must not hit merely because selected columns match.
  base <- sb_pipeline("none", "logistic")
  chain <- sb_new_chain_pipeline(base, list(
    sb_pipeline("fixed_features", "logistic", selector_features = features[1:2],
      feature_set_id = "first_two"),
    sb_pipeline("svm_rfe_binary", "logistic", selector_options = list(inner_folds = 3L))))
  adapter <- sb_resolve_pipeline_runtime(chain)$selector
  endpoint <- list(task_type = "binary", outcome = d$outcome, group_ids = d$patient_id)
  x <- d[features]; changed <- x; changed$gene_c <- -changed$gene_c
  calls <- 0L
  sb_with_selector_reuse({
    a <- adapter$fit(x, endpoint)
    expect_identical(adapter$fit(x, endpoint), a)
    adapter$fit(changed, endpoint)
  })
  expect_identical(calls, 2L)
})

# ---- selector-reuse-doctor ----

test_that("doctor comparison reuses selectors without additional user settings", {
  d <- sb_example_survival(120, 37)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month")
  pipelines <- list(
    cox = surv_pipeline("Cox", surv_selector("StepCox")),
    chain = surv_pipeline("Cox", surv_selector(c("StepCox", "StepCox"))))
  original <- sb_fit_runtime_selector_uncached
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    if (identical(selector$id, "stepcox")) calls <<- calls + 1L
    original(selector, x, endpoint)
  })
  plain <- sb_compare_cv(task, unname(pipelines), folds = 3, seed = 19,
    reuse_selectors = FALSE, progress = FALSE)
  plain_calls <- calls
  calls <- 0L
  run <- surv_run(task, pipelines = pipelines, folds = 3, seed = 19,
    plots = FALSE)
  expect_identical(run$comparison, plain)
  expect_identical(run$comparison$comparison$status, rep("success", 2))
  expect_lt(calls, plain_calls)
  expect_length(run$details, 0L)
})

# ---- selector-reuse-inner ----

test_that("sharing respects outer and inner glmnet analysis folds", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("Boruta")
  skip_if_not_installed("ranger")
  calls <- list()
  local_mocked_bindings(sb_fit_runtime_selector_uncached=function(selector,x,endpoint) {
    calls[[length(calls)+1L]] <<- list(rows=rownames(x),endpoint=endpoint)
    sb_new_forest_selector_fit(selector$id,names(x),names(x)[1:2],
      sb_forest_selector_params(selector$id),
      data.frame(feature=names(x),selected=names(x)%in%names(x)[1:2]))
  })
  d <- sb_example_binary(n=120,seed=1208)
  task <- sb_binary_task(d,"patient_id","outcome",c("gene_a","gene_b","gene_c"))
  pipelines <- lapply(c(1,.5),function(alpha)sb_pipeline(
    "boruta_binary_roughfix","glmnet_logistic",glmnet_alpha=alpha,
    inner_folds=3L,inner_seed=101L))
  plain <- sb_compare_cv(task,pipelines,folds=3L,seed=1208L,reuse_selectors=FALSE)
  plain_calls <- calls; calls <- list()
  shared <- sb_compare_cv(task,pipelines,folds=3L,seed=1208L)
  expect_true(all(shared$comparison$status=="success"))
  expect_equal(shared,plain,tolerance=1e-12)
  expect_length(plain_calls,24L)
  expect_length(calls,12L)
  # No full-development selector fit is allowed during CV, and inner fits
  # must have smaller analysis sets than the outer fits they belong to.
  sizes <- vapply(calls,function(x)length(x$rows),integer(1))
  expect_equal(sum(sizes==80L),3L)
  expected <- c(
    lapply(shared$results[[1]]$folds,function(f)rownames(task$data)[f$analysis_rows]),
    unlist(lapply(shared$results[[1]]$tuning,function(t)
      lapply(t$inner_selections,function(s)s$analysis_patient_ids)),recursive=FALSE))
  expect_length(expected,12L)
  # Stratifying each class can yield unequal inner-fold sizes. Compare the
  # exact declared analysis identities instead of assuming balanced sizes.
  expect_true(all(vapply(expected,function(ids)
    sum(vapply(calls,function(call)identical(call$rows,ids),logical(1)))==1L,
    logical(1))))
  expect_true(all(sizes<nrow(task$data)))
  expect_null(.sb_selector_reuse_context$active)
})

test_that("seeded selector reuse leaves caller random state unchanged", {
  withr::local_seed(891L)
  local_mocked_bindings(sb_fit_runtime_selector_uncached=function(selector,x,endpoint) {
    withr::with_seed(1L,list(selected_features=names(x),value=runif(3)))
  })
  selector <- sb_resolve_selector_adapter("boruta_binary_roughfix",list())
  x <- data.frame(a=1:4); endpoint <- list(task_type="binary",outcome=c(0,1,0,1))
  before <- .Random.seed
  out <- sb_with_selector_reuse({
    first <- sb_fit_runtime_selector(selector,x,endpoint)
    expect_identical(.Random.seed,before)
    second <- sb_fit_runtime_selector(selector,x,endpoint)
    expect_identical(first,second)
    expect_identical(.Random.seed,before)
  },diagnostics=TRUE)
  expect_identical(out$reuse$hits,1L)
  expect_identical(.Random.seed,before)
})

# ---- selector-reuse-lasso-prefix ----

test_that("Lasso selection and its StepCox prefix share actual fitted states", {
  skip_if_not_installed("glmnet")
  withr::local_seed(241L)
  x <- data.frame(a = rnorm(120), b = rnorm(120), noise = rnorm(120))
  rownames(x) <- paste0("patient_", seq_len(nrow(x)))
  time <- rexp(120, exp(1.5 * x$a - x$b))
  censor <- rexp(120, 0.12)
  e <- list(task_type = "survival", time = pmin(time, censor),
    event = as.integer(time <= censor))
  first <- sb_lasso_nonzero_selector_adapter(3L, 17L)
  chain <- sb_lasso_then_stepcox_selector_adapter(3L, 17L, "forward")
  before <- .Random.seed
  expected <- sb_with_selector_reuse({
    list(first = sb_fit_runtime_selector(first, x, e),
      chain = sb_fit_runtime_selector(chain, x, e))
  }, enabled = FALSE)
  actual <- sb_with_selector_reuse({
    list(first = sb_fit_runtime_selector(first, x, e),
      chain = sb_fit_runtime_selector(chain, x, e))
  }, diagnostics = TRUE)
  expect_identical(actual$result, expected)
  expect_identical(actual$result$first, actual$result$chain$stages$first)
  expect_identical(actual$reuse$hits, 1L)
  expect_identical(actual$reuse$entries, 2L)
  expect_identical(.Random.seed, before)

  settings <- sb_with_selector_reuse({
    sb_fit_runtime_selector(first, x, e)
    sb_fit_runtime_selector(sb_lasso_nonzero_selector_adapter(3L, 18L), x, e)
    sb_fit_runtime_selector(sb_lasso_nonzero_selector_adapter(4L, 17L), x, e)
  }, diagnostics = TRUE)
  expect_identical(settings$reuse$hits, 0L)
  expect_identical(settings$reuse$misses, 3L)
})

test_that("Lasso reuse preserves outer selection and complete predictions", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("survivalsvm")
  withr::local_seed(241L)
  d <- data.frame(patient_id = paste0("p", 1:120), a = rnorm(120), b = rnorm(120))
  time <- rexp(120, exp(1.5 * d$a - d$b))
  censor <- rexp(120, 0.12)
  d$time <- pmin(time, censor)
  d$event <- as.integer(time <= censor)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("a", "b"), "day")
  pipelines <- lapply(c(1, 2), function(gamma)
    sb_pipeline("lasso_nonzero", "survivalsvm",
      selector_inner_folds = 3L, selector_seed = 17L,
      survivalsvm_gamma_mu = gamma))
  original <- sb_fit_runtime_selector_uncached
  calls <- list()
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls[[length(calls) + 1L]] <<- rownames(x)
    original(selector, x, endpoint)
  })
  plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 81L, reuse_selectors = FALSE)
  plain_calls <- length(calls)
  calls <- list()
  shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 81L)
  expect_true(all(shared$comparison$status == "success"))
  expect_identical(shared, plain)
  expect_lt(length(calls), plain_calls)
  expected <- lapply(shared$results[[1]]$folds,
    function(f) rownames(task$data)[f$analysis_rows])
  expect_true(all(vapply(expected, function(ids)
    sum(vapply(calls, function(rows) identical(rows, ids), logical(1))) == 1L,
    logical(1))))
  expect_true(all(lengths(calls) < nrow(d)))
})

# ---- selector-reuse-models ----

test_that("model selectors reuse native fits without changing results or RNG", {
  for (package in c("glmnet", "mboost", "randomForest", "gbm", "e1071",
                    "xgboost", "plsRglm", "abess")) {
    skip_if_not_installed(package)
    requireNamespace(package, quietly = TRUE)
  }
  native <- sb_fit_runtime_selector_uncached
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    native(selector, x, endpoint)
  })
  d <- sb_example_binary(160L, 1208L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "binary", outcome = d$outcome,
    group_ids = d$patient_id)
  withr::local_seed(37L)
  before <- .Random.seed
  for (id in c(sb_binary_model_selector_ids(), sb_abess_selector_ids())) {
    if (id == "abess_cox_gic_nonzero") {
      s <- sb_example_survival(160L, 11032L)
      endpoint <- list(task_type = "survival", time = s$time, event = s$event,
        group_ids = s$patient_id)
      x <- s[c("gene_a", "gene_b", "gene_c")]
    }
    adapter <- if (id %in% sb_binary_model_selector_ids()) {
      sb_binary_model_selector_adapter(id, sb_binary_model_selector_options(id))
    } else sb_abess_selector_adapter(id)
    plain <- native(adapter, x, endpoint)
    calls <- 0L
    shared <- sb_with_selector_reuse({
      a <- sb_fit_runtime_selector(adapter, x, endpoint)
      b <- sb_fit_runtime_selector(adapter, x, endpoint)
      list(a, b)
    }, diagnostics = TRUE)
    expect_identical(shared$result[[1]], plain, info = id)
    expect_identical(shared$result[[2]], plain, info = id)
    expect_identical(calls, 1L, info = id)
    expect_identical(shared$reuse$hits, 1L, info = id)
    expect_identical(.Random.seed, before, info = id)
  }
})

test_that("shared selectors keep diagnostic history requests separate", {
  for (package in c("glmnet", "randomForest", "abess")) skip_if_not_installed(package)
  d <- sb_example_binary(160L, 1208L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "binary", outcome = d$outcome)
  for (id in c("lasso_binary_nonzero", "rf_binary_importance", "abess_binary_gic_nonzero")) {
    adapters <- lapply(c(FALSE, TRUE), function(history) {
      if (id == "abess_binary_gic_nonzero") sb_abess_selector_adapter(id, history)
      else sb_binary_model_selector_adapter(id, sb_binary_model_selector_options(id), history)
    })
    plain <- lapply(adapters, sb_fit_runtime_selector_uncached, x = x, endpoint = endpoint)
    shared <- sb_with_selector_reuse({
      lapply(c(adapters, adapters), sb_fit_runtime_selector, x = x, endpoint = endpoint)
    }, diagnostics = TRUE)
    expect_identical(shared$result, c(plain, plain), info = id)
    expect_identical(shared$reuse$misses, 2L, info = id)
    expect_identical(shared$reuse$hits, 2L, info = id)
  }
})

test_that("ABESS shares distinct outer and inner fits without sharing across folds", {
  skip_if_not_installed("abess")
  skip_if_not_installed("glmnet")
  native <- sb_fit_runtime_selector_uncached
  calls <- list()
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls[[length(calls) + 1L]] <<- rownames(x)
    native(selector, x, endpoint)
  })
  for (type in c("binary", "survival")) {
    features <- c("gene_a", "gene_b", "gene_c")
    if (type == "binary") {
      d <- sb_example_binary(160L, 1208L)
      task <- sb_binary_task(d, "patient_id", "outcome", features)
      selector <- "abess_binary_gic_nonzero"; learner <- "glmnet_logistic"
    } else {
      d <- sb_example_survival(160L, 11032L)
      task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
      selector <- "abess_cox_gic_nonzero"; learner <- "glmnet_cox"
    }
    pipelines <- lapply(c(1, .5), function(alpha) sb_pipeline(selector, learner,
      glmnet_alpha = alpha, inner_folds = 3L, inner_seed = 101L))
    calls <- list()
    plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L, reuse_selectors = FALSE)
    expect_length(calls, 24L)
    calls <- list()
    shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L)
    expect_identical(shared, plain)
    expect_identical(unique(shared$comparison$status), "success")
    expect_length(calls, 12L)
    expected <- c(lapply(shared$results[[1]]$folds, function(f) rownames(task$data)[f$analysis_rows]),
      unlist(lapply(shared$results[[1]]$tuning, function(t)
        lapply(t$inner_selections, function(s) s$analysis_patient_ids)), recursive = FALSE))
    expect_identical(vapply(expected, function(ids) sum(vapply(calls,
      function(rows) identical(rows, ids), logical(1))), integer(1)), rep(1L, 12L))
  }
})

test_that("binary model selectors preserve complete comparisons with sharing", {
  for (package in c("glmnet", "mboost", "randomForest", "gbm", "e1071",
                    "xgboost", "plsRglm")) skip_if_not_installed(package)
  d <- sb_example_binary(160L, 1208L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  native <- sb_fit_runtime_selector_uncached
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    native(selector, x, endpoint)
  })
  for (id in sb_binary_model_selector_ids()) {
    pipelines <- lapply(c("logistic", "gaussian_naive_bayes"),
      function(learner) sb_pipeline(id, learner))
    calls <- 0L
    plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L, reuse_selectors = FALSE)
    expect_identical(calls, 6L, info = id)
    calls <- 0L
    shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L)
    expect_identical(shared, plain, info = id)
    failed <- sum(vapply(shared$results[[1]]$selections,
      function(s) !identical(s$status, "success"), logical(1)))
    expect_identical(calls, 3L + failed, info = id)
    # Empty native selections are retried for each consumer, not cached or filled.
    if (failed) expect_identical(unique(shared$comparison$status), "failed", info = id)
  }
})

# ---- selector-reuse ----

test_that("reuse matches complete inputs, parameters and RNG kind", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector,x,endpoint) {
    calls <<- calls + 1L
    list(selected_features=names(x), value=calls)
  })
  selector <- sb_resolve_selector_adapter("boruta_binary_roughfix",list())
  x <- data.frame(a=1:6,b=6:1,row.names=paste0("p",1:6))
  endpoint <- list(task_type="binary",outcome=rep(0:1,3))
  original_kind <- RNGkind()
  withr::defer(do.call(RNGkind,as.list(original_kind)))
  result <- sb_with_selector_reuse({
    first <- sb_fit_runtime_selector(selector,x,endpoint)
    expect_identical(sb_fit_runtime_selector(selector,x,endpoint),first)
    changed <- x; changed$a[1] <- 99L
    sb_fit_runtime_selector(selector,changed,endpoint)
    sb_fit_runtime_selector(selector,x[6:1,],endpoint)
    sb_fit_runtime_selector(selector,x[c("b","a")],endpoint)
    changed <- x; rownames(changed) <- paste0("q",1:6)
    sb_fit_runtime_selector(selector,changed,endpoint)
    changed_endpoint <- endpoint; changed_endpoint$outcome[1] <- 1L
    sb_fit_runtime_selector(selector,x,changed_endpoint)
    changed_selector <- selector; changed_selector$params <- list(seed=2L)
    sb_fit_runtime_selector(changed_selector,x,endpoint)
    RNGkind(normal.kind=if(original_kind[2]=="Inversion")"Box-Muller" else "Inversion")
    sb_fit_runtime_selector(selector,x,endpoint)
    do.call(RNGkind,as.list(original_kind))
    first
  },diagnostics=TRUE)
  expect_identical(calls,8L)
  expect_identical(result$reuse$hits,1L)
  expect_null(.sb_selector_reuse_context$active)
  sb_fit_runtime_selector(selector,x,endpoint)
  expect_identical(calls,9L)
})

test_that("warnings, errors and unoptimized selectors are executed again", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached=function(selector,x,endpoint) {
    calls <<- calls + 1L
    if(selector$params$mode=="warning")warning("same warning")
    if(selector$params$mode=="error")stop("same error")
    list(selected_features=names(x))
  })
  selector <- sb_resolve_selector_adapter("boruta_binary_roughfix",list())
  x <- data.frame(a=1:4); endpoint <- list(task_type="binary",outcome=c(0,1,0,1))
  sb_with_selector_reuse({
    selector$params <- list(mode="warning")
    for(i in 1:2)expect_warning(sb_fit_runtime_selector(selector,x,endpoint),"same warning")
    selector$params <- list(mode="error")
    for(i in 1:2)expect_error(sb_fit_runtime_selector(selector,x,endpoint),"same error")
    selector$id <- "none"; selector$params <- list(mode="none")
    for(i in 1:2)sb_fit_runtime_selector(selector,x,endpoint)
  })
  expect_identical(calls,6L)
  expect_error(sb_with_selector_reuse(stop("cleanup")),"cleanup")
  expect_null(.sb_selector_reuse_context$active)
})

test_that("disabled and nested contexts restore the outer context", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached=function(selector,x,endpoint) {
    calls <<- calls + 1L; list(selected_features=names(x))
  })
  selector <- sb_resolve_selector_adapter("boruta_binary_roughfix",list())
  x <- data.frame(a=1:4); endpoint <- list(task_type="binary",outcome=c(0,1,0,1))
  out <- sb_with_selector_reuse({
    sb_fit_runtime_selector(selector,x,endpoint)
    sb_with_selector_reuse({
      sb_fit_runtime_selector(selector,x,endpoint)
      sb_fit_runtime_selector(selector,x,endpoint)
    },enabled=FALSE)
    sb_with_selector_reuse({
      sb_fit_runtime_selector(selector,x,endpoint)
      sb_fit_runtime_selector(selector,x,endpoint)
    })
    sb_fit_runtime_selector(selector,x,endpoint)
  },diagnostics=TRUE)
  expect_identical(calls,4L)
  expect_identical(out$reuse$hits,1L)
  expect_null(.sb_selector_reuse_context$active)
  expect_error(sb_with_selector_reuse(1,enabled=NA),class="stormbreaker_bad_selector_reuse")
})
})

# ---- selector-reuse-native ----
local({
# ---- selector-reuse-boosting ----

test_that("survival boosting selectors share only identical successful fits", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    list(selected_features = names(x), params = selector$params)
  })
  x <- data.frame(a = 1:6, b = 6:1)
  endpoint <- list(task_type = "survival", time = 1:6, event = rep(0:1, 3))
  out <- sb_with_selector_reuse({
    for (id in sb_survival_boosting_selector_ids()) {
      p <- sb_pipeline(id, "coxph")
      s <- sb_resolve_selector_adapter(id, p$selector_params)
      first <- sb_fit_runtime_selector(s, x, endpoint)
      expect_identical(sb_fit_runtime_selector(s, x, endpoint), first)
      s$params$seed <- 9L
      sb_fit_runtime_selector(s, x, endpoint)
      sb_fit_runtime_selector(s, x[1:4, ], list(task_type = "survival", time = 1:4, event = rep(0:1, 2)))
    }
  }, diagnostics = TRUE)
  expect_identical(calls, 15L)
  expect_identical(out$reuse$hits, 5L)
  expect_null(.sb_selector_reuse_context$active)
})

test_that("boosting warnings and failures are still evaluated on every request", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    if (selector$params$mode == "warning") warning("original warning")
    if (selector$params$mode == "error") sb_abort("original error", "stormbreaker_probe_failure")
    list(selected_features = names(x))
  })
  s <- sb_resolve_selector_adapter("gbm_cox_importance", list())
  sb_with_selector_reuse({
    s$params <- list(mode = "warning")
    for (i in 1:2) expect_warning(sb_fit_runtime_selector(s, data.frame(x = 1:4), list()), "original warning")
    s$params <- list(mode = "error")
    for (i in 1:2) expect_error(sb_fit_runtime_selector(s, data.frame(x = 1:4), list()), class = "stormbreaker_probe_failure")
  })
  expect_identical(calls, 4L)
})

test_that("native boosting selection, final predictions and RNG survive reuse", {
  for (pkg in c("gbm", "mboost", "superpc", "xgboost")) skip_if_not_installed(pkg)
  d <- sb_example_survival(150, 82)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
  pipelines <- unlist(lapply(sb_survival_boosting_selector_ids(), function(id) {
    lapply(c("coxph", "survreg_weibull"), function(learner) sb_pipeline(id, learner))
  }), recursive = FALSE)
  run <- function(enabled) sb_with_selector_reuse(lapply(pipelines, function(p) {
    f <- sb_fit_final(task, p, capture_diagnostics = TRUE)
    expect_identical(f$status, "success")
    prediction <- sb_predict_final(unserialize(serialize(f, NULL)), d[c("patient_id", task$predictors)], "patient_id")
    list(selection = f$selection, prediction = prediction)
  }), enabled = enabled, diagnostics = TRUE)
  withr::local_seed(847L)
  before <- .Random.seed
  plain <- run(FALSE)
  expect_identical(.Random.seed, before)
  shared <- run(TRUE)
  expect_identical(.Random.seed, before)
  expect_identical(shared$result, plain$result)
  expect_identical(shared$reuse$hits, 5L)
})

test_that("boosting reuse retains separate outer and glmnet inner analysis folds", {
  skip_if_not_installed("mboost")
  skip_if_not_installed("glmnet")
  d <- sb_example_survival(150, 82)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
  calls <- list()
  original <- sb_fit_runtime_selector_uncached
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls[[length(calls) + 1L]] <<- rownames(x)
    original(selector, x, endpoint)
  })
  pipelines <- lapply(c(0, .5), function(alpha) sb_pipeline("gamboost_cox_nonzero", "glmnet_cox",
    glmnet_alpha = alpha, inner_folds = 3L, inner_seed = 7L))
  plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 82L, reuse_selectors = FALSE)
  plain_calls <- calls; calls <- list()
  shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 82L)
  expect_identical(shared, plain)
  expect_length(plain_calls, 24L)
  expect_length(calls, 12L)
  expected <- c(lapply(shared$results[[1L]]$folds, function(f) rownames(task$data)[f$analysis_rows]),
    unlist(lapply(shared$results[[1L]]$tuning, function(t) lapply(t$inner_selections, function(s) s$analysis_patient_ids)), recursive = FALSE))
  expect_true(all(lengths(calls) < nrow(task$data)))
  expect_true(all(vapply(expected, function(ids) sum(vapply(calls, identical, logical(1), y = ids)) == 1L, logical(1))))
})

# ---- selector-reuse-coxboost-rsf ----

test_that("native CoxBoost and RSF selections replay without RNG side effects", {
  skip_if_not_installed("CoxBoost")
  skip_if_not_installed("randomForestSRC")
  withr::local_seed(8601L)
  x <- as.data.frame(matrix(rnorm(90 * 6), nrow = 90))
  names(x) <- paste0("gene", 1:6)
  rownames(x) <- paste0("patient", 1:90)
  time <- rexp(90, exp(1.2 * x$gene1 - 0.9 * x$gene2 + 0.6 * x$gene3) / 15)
  censor <- rexp(90, 1 / 20)
  endpoint <- list(task_type = "survival", time = pmax(pmin(time, censor), 0.05),
    event = as.integer(time <= censor))
  for (selector in list(sb_coxboost_nonzero_selector_adapter(3L, 1L),
    sb_rsf_minimal_depth_selector_adapter(5L, 1L))) {
    before <- .Random.seed
    plain <- sb_with_selector_reuse(sb_fit_runtime_selector(selector, x, endpoint),
      enabled = FALSE)
    out <- sb_with_selector_reuse({
      first <- sb_fit_runtime_selector(selector, x, endpoint)
      second <- sb_fit_runtime_selector(selector, x, endpoint)
      expect_identical(first, second)
      first
    }, diagnostics = TRUE)
    expect_s3_class(plain, "stormbreaker_selector_fit")
    expect_identical(out$result, plain)
    expect_identical(out$reuse$hits, 1L)
    expect_identical(out$reuse$misses, 1L)
    expect_identical(.Random.seed, before)
  }
})

test_that("CoxBoost and RSF reuse distinguishes settings and never stores failures", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    if (identical(endpoint$mode, "warning")) warning("native warning")
    if (identical(endpoint$mode, "error")) stop("native error")
    list(selected_features = names(x))
  })
  x <- data.frame(a = 1:12, b = 12:1, row.names = paste0("p", 1:12))
  endpoint <- list(task_type = "survival", time = as.numeric(1:12), event = rep(0:1, 6))
  for (selector in list(sb_coxboost_nonzero_selector_adapter(3L, 1L),
    sb_rsf_minimal_depth_selector_adapter(5L, 1L))) {
    out <- sb_with_selector_reuse({
      sb_fit_runtime_selector(selector, x, endpoint)
      sb_fit_runtime_selector(selector, x, endpoint)
      changed <- selector
      changed$params$seed <- 2L
      sb_fit_runtime_selector(changed, x, endpoint)
      changed <- selector
      if (selector$id == "coxboost_nonzero") changed$params$inner_folds <- 4L
      else changed$params$nodesize <- 6L
      sb_fit_runtime_selector(changed, x, endpoint)
      changed_endpoint <- endpoint
      changed_endpoint$time[1] <- 0.5
      sb_fit_runtime_selector(selector, x, changed_endpoint)
      sb_fit_runtime_selector(selector, x[12:1, ], endpoint)
      for (i in 1:2) {
        expect_warning(sb_fit_runtime_selector(selector, x,
          c(endpoint, list(mode = "warning"))), "native warning")
        expect_error(sb_fit_runtime_selector(selector, x,
          c(endpoint, list(mode = "error"))), "native error")
      }
    }, diagnostics = TRUE)
    expect_identical(out$reuse$hits, 1L)
    expect_identical(out$reuse$entries, 5L)
    expect_identical(out$reuse$misses, 9L)
  }
  expect_identical(calls, 18L)
})

# ---- selector-reuse-information ----

test_that("native information selectors share matching outer and inner analysis fits", {
  skip_if_not_installed("praznik")
  skip_if_not_installed("glmnet")
  native <- sb_fit_runtime_selector_uncached
  calls <- list()
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls[[length(calls) + 1L]] <<- rownames(x)
    native(selector, x, endpoint)
  })
  withr::local_seed(1208L)
  d <- sb_example_binary(120L, 1208L)
  for (i in seq_len(21L)) d[[paste0("noise", i)]] <- rnorm(nrow(d))
  features <- c("gene_a", "gene_b", "gene_c", paste0("noise", seq_len(21L)))
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  before <- .Random.seed
  for (selector in c("mrmr_binary_nlogn", sb_information_binary_selector_ids())) {
    pipelines <- lapply(c(1, .5), function(alpha) sb_pipeline(selector,
      "glmnet_logistic", glmnet_alpha = alpha, inner_folds = 3L, inner_seed = 101L))
    calls <- list()
    plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L,
      reuse_selectors = FALSE)
    expect_length(calls, 24L)
    calls <- list()
    shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L)
    expect_identical(shared, plain)
    expect_identical(unique(shared$comparison$status), "success")
    expect_length(calls, 12L)
    expected <- c(lapply(shared$results[[1]]$folds,
      function(f) rownames(task$data)[f$analysis_rows]),
      unlist(lapply(shared$results[[1]]$tuning, function(t)
        lapply(t$inner_selections, function(s) s$analysis_patient_ids)), recursive = FALSE))
    expect_identical(vapply(expected, function(ids) sum(vapply(calls,
      function(rows) identical(rows, ids), logical(1))), integer(1)), rep(1L, 12L))
    expect_identical(.Random.seed, before)
  }
})

# ---- selector-reuse-omics ----

test_that("seeded omics and PAMR selectors share complete native selections", {
  for (package in c("SIS", "VSURF", "lpc", "stabs", "pamr")) {
    skip_if_not_installed(package)
    requireNamespace(package, quietly = TRUE)
  }
  fixture <- withr::with_seed(9201L, {
    n <- 96L; p <- 60L
    y <- rep(0:1, length.out = n)
    x <- as.data.frame(matrix(rnorm(n * p), n, p))
    x[[1]] <- x[[1]] + 2 * y
    x[[2]] <- x[[2]] - 1.5 * y
    t <- rexp(n, exp(x[[1]] - x[[2]]))
    c <- rexp(n, .3)
    list(x = x, y = y, time = pmin(t, c) + .01, event = as.integer(t <= c))
  })
  native <- sb_fit_runtime_selector_uncached
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    native(selector, x, endpoint)
  })
  ids <- c("sis_binary_nlogn", "vsurf_binary_prediction", "lpc_binary_nlogn",
    "lpc_cox_nlogn", "stabsel_binary_lasso", "stabsel_cox_lasso", sb_pamr_selector_ids())
  for (id in ids) {
    adapter <- if (id %in% sb_pamr_selector_ids()) {
      sb_pamr_selector_adapter(id, sb_pamr_selector_params(id))
    } else sb_high_dimensional_selector_adapter(id)
    endpoint <- if (grepl("cox|survival", id)) {
      list(task_type = "survival", time = fixture$time, event = fixture$event)
    } else list(task_type = "binary", outcome = fixture$y)
    withr::with_seed(37L, {
      before <- .Random.seed
      plain <- native(adapter, fixture$x, endpoint)
      calls <- 0L
      shared <- sb_with_selector_reuse({
        list(sb_fit_runtime_selector(adapter, fixture$x, endpoint),
          sb_fit_runtime_selector(adapter, fixture$x, endpoint))
      }, diagnostics = TRUE)
      expect_identical(shared$result, list(plain, plain), info = id)
      expect_identical(calls, 1L, info = id)
      expect_identical(shared$reuse$hits, 1L, info = id)
      expect_identical(.Random.seed, before, info = id)
    })
  }
})

# ---- selector-reuse-stepcox ----

test_that("real StepCox fits share only identical analysis and direction", {
  withr::local_seed(241L)
  x <- data.frame(a = rnorm(120), b = rnorm(120), c = rnorm(120))
  rownames(x) <- paste0("patient_", seq_len(nrow(x)))
  endpoint <- list(task_type = "survival",
    time = rexp(120, exp(1.5 * x$a - x$b)), event = rep(1L, 120))
  before <- .Random.seed
  out <- sb_with_selector_reuse({
    for (direction in c("forward", "backward", "both")) {
      selector <- sb_stepcox_selector_adapter(direction)
      first <- sb_fit_runtime_selector(selector, x, endpoint)
      expect_identical(sb_fit_runtime_selector(selector, x, endpoint), first)
      expect_identical(first, sb_fit_runtime_selector_uncached(selector, x, endpoint))
    }
    sb_fit_runtime_selector(selector, x[120:1, ],
      list(task_type = "survival", time = rev(endpoint$time), event = endpoint$event))
  }, diagnostics = TRUE)
  expect_identical(out$reuse$hits, 3L)
  expect_identical(out$reuse$misses, 4L)
  expect_identical(.Random.seed, before)
  expect_null(.sb_selector_reuse_context$active)
})

test_that("StepCox warning and error executions are not reused", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    if (selector$params$direction == "forward") warning("native warning")
    else stop("native failure")
    list(selected_features = names(x))
  })
  out <- sb_with_selector_reuse({
    for (i in 1:2) {
      expect_warning(sb_fit_runtime_selector(sb_stepcox_selector_adapter("forward"),
        data.frame(a = 1:4), list(task_type = "survival")), "native warning")
      expect_error(sb_fit_runtime_selector(sb_stepcox_selector_adapter("both"),
        data.frame(a = 1:4), list(task_type = "survival")), "native failure")
    }
  }, diagnostics = TRUE)
  expect_identical(calls, 4L)
  expect_identical(out$reuse$misses, 4L)
  expect_identical(out$reuse$entries, 0L)
})

test_that("StepCox reuse preserves complete nested Cox comparisons", {
  skip_if_not_installed("glmnet")
  withr::local_seed(241L)
  d <- data.frame(patient_id = paste0("p", 1:120), a = rnorm(120), b = rnorm(120))
  d$time <- rexp(120, exp(1.5 * d$a - d$b))
  censor <- rexp(120, 0.12)
  d$event <- as.integer(d$time <= censor)
  d$time <- pmin(d$time, censor)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("a", "b"), "day")
  pipelines <- lapply(c(0, 1), function(alpha) sb_pipeline("stepcox", "glmnet_cox",
    stepcox_direction = "forward", glmnet_alpha = alpha, inner_folds = 3L,
    inner_seed = 17L))
  original <- sb_fit_runtime_selector_uncached
  calls <- list()
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls[[length(calls) + 1L]] <<- rownames(x)
    original(selector, x, endpoint)
  })
  plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 81L, reuse_selectors = FALSE)
  plain_calls <- length(calls)
  calls <- list()
  shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 81L)
  expect_true(all(shared$comparison$status == "success"))
  expect_identical(shared, plain)
  expect_lt(length(calls), plain_calls)
  expect_true(all(lengths(calls) < nrow(d)))
  expected <- c(lapply(shared$results[[1]]$folds,
    function(f) rownames(task$data)[f$analysis_rows]),
    unlist(lapply(shared$results[[1]]$tuning, function(t)
      lapply(t$inner_selections, function(s) s$analysis_patient_ids)), recursive = FALSE))
  expect_true(all(vapply(expected, function(ids)
    sum(vapply(calls, function(rows) identical(rows, ids), logical(1))) == 1L,
    logical(1))))
})

# ---- selector-reuse-univ ----

test_that("native univariable selectors share only matching outer and inner fits", {
  skip_if_not_installed("glmnet")
  native <- sb_fit_runtime_selector_uncached
  calls <- list()
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls[[length(calls) + 1L]] <<- rownames(x)
    native(selector, x, endpoint)
  })
  for (type in c("binary", "survival")) {
    features <- c("gene_a", "gene_b", "gene_c")
    if (type == "binary") {
      d <- sb_example_binary(120L, 1208L)
      task <- sb_binary_task(d, "patient_id", "outcome", features)
      selector <- "univ_binary_auc_top_k"; learner <- "glmnet_logistic"
    } else {
      d <- sb_example_survival(120L, 11032L)
      task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
      selector <- "univ_cox_top_k"; learner <- "glmnet_cox"
    }
    pipelines <- lapply(c(1, .5), function(alpha) sb_pipeline(selector, learner,
      selector_n_features = 2L, glmnet_alpha = alpha, inner_folds = 3L, inner_seed = 101L))
    withr::local_seed(23L)
    before <- .Random.seed
    calls <- list()
    plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L, reuse_selectors = FALSE)
    expect_length(calls, 24L)
    calls <- list()
    shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L)
    expect_identical(shared, plain)
    expect_identical(unique(shared$comparison$status), "success")
    expect_length(calls, 12L)
    expected <- c(lapply(shared$results[[1]]$folds, function(f) rownames(task$data)[f$analysis_rows]),
      unlist(lapply(shared$results[[1]]$tuning, function(t)
        lapply(t$inner_selections, function(s) s$analysis_patient_ids)), recursive = FALSE))
    expect_identical(length(expected), 12L)
    expect_identical(vapply(expected, function(ids) sum(vapply(calls,
      function(rows) identical(rows, ids), logical(1))), integer(1)), rep(1L, 12L))
    expect_identical(.Random.seed, before)
  }
})

# ---- selector-reuse-weights ----

weight_reuse_pipeline <- function(id, learner = "coxph") {
  if (id %in% c("lasso_nonzero", "enet_cox_nonzero")) {
    sb_pipeline(id, learner, selector_inner_folds = 10L, selector_seed = 1L)
  } else sb_pipeline(id, learner)
}

test_that("weight and Elastic Net selectors reuse only the same analysis inputs", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    list(selected_features = names(x), params = selector$params)
  })
  x <- data.frame(a = 1:6, b = 6:1, row.names = paste0("p", 1:6))
  endpoint <- list(task_type = "survival", time = 1:6, event = rep(0:1, 3),
    group_ids = c(1, 1, 2, 2, 3, 3))
  ids <- c(sb_survival_weight_selector_ids(), "enet_cox_nonzero")
  out <- sb_with_selector_reuse({
    for (id in ids) {
      p <- weight_reuse_pipeline(id)
      s <- sb_resolve_selector_adapter(id, p$selector_params)
      first <- sb_fit_runtime_selector(s, x, endpoint)
      expect_identical(sb_fit_runtime_selector(s, x, endpoint), first)
      changed <- endpoint; changed$group_ids <- 1:6
      sb_fit_runtime_selector(s, x, changed)
      changed <- endpoint; changed$event[1] <- 1L
      sb_fit_runtime_selector(s, x, changed)
      sb_fit_runtime_selector(s, x[c("b", "a")], endpoint)
      changed_s <- s; changed_s$params$probe_parameter <- 2L
      sb_fit_runtime_selector(changed_s, x, endpoint)
    }
  }, diagnostics = TRUE)
  expect_identical(calls, 5L * length(ids))
  expect_identical(out$reuse$hits, length(ids))
  expect_null(.sb_selector_reuse_context$active)
})

test_that("weight and Elastic Net warning or failed fits are not reused", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    if (selector$params$probe == "warning") warning("native warning")
    if (selector$params$probe == "error") sb_abort("native error", "stormbreaker_probe_failure")
    list(selected_features = names(x))
  })
  sb_with_selector_reuse({
    for (id in c(sb_survival_weight_selector_ids(), "enet_cox_nonzero")) {
      p <- weight_reuse_pipeline(id)
      s <- sb_resolve_selector_adapter(id, p$selector_params)
      s$params$probe <- "warning"
      for (i in 1:2) expect_warning(sb_fit_runtime_selector(s, data.frame(x = 1:4), list()), "native warning")
      s$params$probe <- "error"
      for (i in 1:2) expect_error(sb_fit_runtime_selector(s, data.frame(x = 1:4), list()), class = "stormbreaker_probe_failure")
    }
  })
  expect_identical(calls, 16L)
})

test_that("sharing keeps captured coefficient paths and loadings distinct", {
  for (pkg in c("glmnet", "plsRcox")) skip_if_not_installed(pkg)
  d <- sb_example_survival(150, 82)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "survival", time = d$time, event = d$event)
  for (id in c("ridge_cox_weights", "plsrcox_weights", "lasso_nonzero", "enet_cox_nonzero")) {
    p <- weight_reuse_pipeline(id)
    ordinary <- sb_resolve_selector_adapter(id, p$selector_params)
    captured <- sb_captured_selector_adapter(id, p$selector_params)
    run <- function(enabled) sb_with_selector_reuse({
      a <- sb_fit_runtime_selector(ordinary, x, endpoint)
      b <- sb_fit_runtime_selector(captured, x, endpoint)
      c <- sb_fit_runtime_selector(captured, x, endpoint)
      list(a, b, c)
    }, enabled = enabled, diagnostics = TRUE)
    plain <- run(FALSE); shared <- run(TRUE)
    expect_identical(shared$result, plain$result)
    expect_identical(shared$reuse$hits, 1L)
    history <- if (id == "plsrcox_weights") shared$result[[2]]$evidence$loadings else
      if (id == "ridge_cox_weights") shared$result[[2]]$evidence$coefficient_path else shared$result[[2]]$coefficient_path
    expect_false(is.null(history))
  }
})

test_that("CoxBoost reuse also distinguishes requested diagnostic history", {
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    list(history = environment(selector$fit)$capture_history)
  })
  out <- sb_with_selector_reuse({
    x <- data.frame(a = 1:6, b = 6:1)
    ordinary <- sb_coxboost_nonzero_selector_adapter(3L, 1L)
    captured <- sb_coxboost_nonzero_selector_adapter(3L, 1L, TRUE)
    list(sb_fit_runtime_selector(ordinary, x, list()),
      sb_fit_runtime_selector(captured, x, list()),
      sb_fit_runtime_selector(captured, x, list()))
  }, diagnostics = TRUE)
  expect_identical(out$result, list(list(history = FALSE), list(history = TRUE), list(history = TRUE)))
  expect_identical(out$reuse$hits, 1L)
})

test_that("native weight reuse preserves RNG and outcome-blind final predictions", {
  for (pkg in c("glmnet", "plsRcox", "survivalsvm")) skip_if_not_installed(pkg)
  d <- sb_example_survival(150, 82)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
  pipelines <- unlist(lapply(c(sb_survival_weight_selector_ids(), "enet_cox_nonzero"), function(id)
    lapply(c("coxph", "survreg_weibull"), function(learner) weight_reuse_pipeline(id, learner))), recursive = FALSE)
  run <- function(enabled) sb_with_selector_reuse(lapply(pipelines, function(p) {
    f <- sb_fit_final(task, p, capture_diagnostics = TRUE)
    expect_identical(f$status, "success")
    pred <- sb_predict_final(unserialize(serialize(f, NULL)), d[c("patient_id", task$predictors)], "patient_id")
    expect_identical(pred$status, "success")
    list(selection = f$selection, prediction = pred)
  }), enabled = enabled, diagnostics = TRUE)
  withr::local_seed(847L)
  before <- .Random.seed
  plain <- run(FALSE); expect_identical(.Random.seed, before)
  shared <- run(TRUE); expect_identical(.Random.seed, before)
  expect_identical(shared$result, plain$result)
  expect_identical(shared$reuse$hits, 4L)
})

test_that("weight reuse does not substitute outer selection for glmnet inner selection", {
  for (pkg in c("glmnet", "plsRcox")) skip_if_not_installed(pkg)
  d <- sb_example_survival(150, 82)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
  calls <- list()
  original <- sb_fit_runtime_selector_uncached
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls[[length(calls) + 1L]] <<- rownames(x)
    original(selector, x, endpoint)
  })
  pipelines <- lapply(c(0, .5), function(alpha) sb_pipeline("plsrcox_weights", "glmnet_cox",
    glmnet_alpha = alpha, inner_folds = 3L, inner_seed = 7L))
  plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 82L, reuse_selectors = FALSE)
  plain_calls <- calls; calls <- list()
  shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 82L)
  expect_identical(shared, plain)
  expect_length(plain_calls, 24L)
  expect_length(calls, 12L)
  expected <- c(lapply(shared$results[[1]]$folds, function(f) rownames(task$data)[f$analysis_rows]),
    unlist(lapply(shared$results[[1]]$tuning, function(t) lapply(t$inner_selections, function(s) s$analysis_patient_ids)), recursive = FALSE))
  expect_true(all(lengths(calls) < nrow(task$data)))
  expect_true(all(vapply(expected, function(ids) sum(vapply(calls, identical, logical(1), y = ids)) == 1L, logical(1))))
})
})

# ---- user-models ----
local({
test_that("the menu includes all direct learners and three glmnet settings", {
  registry <- sb_pipeline_registry()
  for (type in c("binary", "survival")) {
    menu <- sb_model_catalog(type)
    expected <- unique(registry$learner[registry$task_type == type &
      registry$selector == "none" & registry$admission_status == "admitted"])
    expect_setequal(unique(menu$learner), expected)
    expect_setequal(menu$glmnet_alpha[!is.na(menu$glmnet_alpha)], c(0, .5, 1))
    pipelines <- sb_model_pipelines(type)
    expect_identical(names(pipelines), menu$model_id)
    expect_true(all(vapply(pipelines, function(p) p$selector == "none", logical(1))))
    expect_length(unique(vapply(pipelines, function(p) p$id, character(1))), nrow(menu))
  }
  expect_error(sb_model_catalog(NULL), "explicitly")
  expect_error(sb_model_pipelines("binary", "coxph"), "model IDs")
  expect_error(sb_model_pipelines("survival", character()), "model IDs")
  expect_error(sb_model_pipelines("binary", c("logistic", "logistic")), "distinct")
})

test_that("users choose task and genes and can request specific models", {
  genes <- c("gene_c", "gene_a")
  for (type in c("binary", "survival")) {
    d <- if (type == "binary") sb_example_binary(60, seed = 231) else
      sb_example_survival(60, seed = 231)
    task <- if (type == "binary") sb_binary_task(d, "patient_id", "outcome", genes) else
      sb_surv_task(d, "patient_id", "time", "event", genes, "month")
    models <- if (type == "binary") c("gaussian_naive_bayes", "logistic") else
      c("survreg_weibull", "coxph")
    result <- sb_compare_models(task, models = models, folds = 3L, seed = 91L, view = FALSE)
    expect_identical(result$input_genes, genes)
    expect_identical(result$requested_models, models)
    expect_identical(result$comparison$model_id, models)
    expect_true(all(result$comparison$complete_oof))
    expect_true(all(vapply(result$results, function(r)
      identical(r$provenance$fold_assignment, result$fold_assignment), logical(1))))
    expect_gte(result$elapsed_seconds, 0)
  }
})

test_that("all-model dispatch keeps missing dependencies and is independent of display", {
  task <- sb_binary_task(sb_example_binary(60, seed = 231), "patient_id", "outcome", "gene_a")
  menu <- sb_model_catalog("binary")
  menu$dependency_status <- "missing"
  dispatched <- list()
  local_mocked_bindings(
    sb_model_catalog = function(...) menu,
    sb_compare_cv = function(task, pipelines, ...) {
      dispatched[[length(dispatched) + 1L]] <<- pipelines
      list(comparison = data.frame(pipeline_id = vapply(pipelines, function(p) p$id, character(1))))
    },
    sb_attach_result_page = function(x, ...) x,
    sb_native_diagnostics = function(...) stop("Comparison must not extract diagnostics")
  )
  quiet <- sb_compare_models(task, models = NULL, view = FALSE)
  visible <- sb_compare_models(task, models = NULL, view = TRUE)
  expect_identical(quiet$requested_models, menu$model_id)
  expect_identical(visible$requested_models, menu$model_id)
  expect_identical(dispatched[[1L]], dispatched[[2L]])
  selected <- rev(menu$model_id[c(1L, 3L)])
  explicit <- sb_compare_models(task, models = selected, view = FALSE)
  expect_identical(explicit$requested_models, selected)
  expect_identical(names(dispatched[[3L]]), selected)
})
})
