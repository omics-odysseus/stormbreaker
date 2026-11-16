# ---- model-catboost ----
local({
# ---- catboost-aft-survival-fixed-adapter ----

test_that("fixed CatBoost AFT reproduces native risk and Normal survival", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  withr::local_seed(8801L)
  x <- data.frame(
    gene_a = stats::rnorm(48L),
    gene_b = stats::rnorm(48L),
    gene_c = stats::rnorm(48L),
    check.names = FALSE
  )
  time <- exp(1 - 0.7 * x$gene_a + 0.25 * x$gene_b +
    stats::rnorm(48L, sd = 0.35))
  event <- rep(c(1L, 1L, 0L), length.out = nrow(x))
  new_x <- x[c(4L, 19L, 37L), , drop = FALSE]
  adapter <- sb_catboost_aft_survival_learner_adapter(
    iterations = 10L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 2.5,
    seed = 8802L
  )
  model <- adapter$fit(x, time, event)

  native_pool <- catboost::catboost.load_pool(
    data = as.matrix(x),
    label = cbind(time, ifelse(event == 1L, time, -1))
  )
  native <- catboost::catboost.train(
    learn_pool = native_pool,
    test_pool = NULL,
    params = adapter$params
  )
  native_margin <- as.numeric(catboost::catboost.predict(
    native,
    catboost::catboost.load_pool(as.matrix(new_x)),
    prediction_type = "RawFormulaVal",
    ntree_start = 0L,
    ntree_end = 10L,
    thread_count = 1L
  ))
  risk <- adapter$predict(model, new_x)
  survival <- adapter$predict_survival(
    model,
    x,
    time,
    event,
    new_x,
    prediction_horizon = 3
  )
  audit <- adapter$extract_tree_model(model)
  state <- unserialize(serialize(adapter$extract_state(model, list(
    x = x,
    time = time,
    event = event,
    prediction_horizon = 3,
    time_unit = "month"
  )), NULL))
  replay <- adapter$predict_state(state, new_x)

  expect_s3_class(model, "stormbreaker_catboost_aft_survival")
  expect_identical(adapter$params$loss_function,
    "SurvivalAft:dist=Normal;scale=1"
  )
  expect_equal(risk, -native_margin, tolerance = 0)
  expect_equal(
    survival$probability,
    stats::pnorm((native_margin - log(3)) / 1),
    tolerance = 0
  )
  expect_identical(survival$source, "catboost_aft_normal_native")
  expect_identical(survival$ties, "not_applicable")
  expect_identical(audit$tree_count, 10L)
  expect_identical(audit$features, names(x))
  expect_equal(replay$risk, risk, tolerance = 0)
  expect_equal(
    replay$survival_probability,
    survival$probability,
    tolerance = 0
  )
  expect_identical(replay$prediction_horizon, 3)
  expect_identical(replay$time_unit, "month")
})

test_that("fixed CatBoost AFT keeps dependency and schema failures explicit", {
  d <- make_sb_fixture(n = 24L, seed = 8811L)
  x <- d[c("gene_a", "gene_b")]
  adapter <- sb_catboost_aft_survival_learner_adapter(
    iterations = 5L,
    learning_rate = 0.1,
    depth = 2L,
    l2_leaf_reg = 3,
    seed = 8812L
  )

  local_mocked_bindings(sb_catboost_available = function() FALSE)
  missing <- capture_sb_error(adapter$fit(x, d$time, d$event))
  expect_s3_class(missing, "stormbreaker_missing_catboost")

  local_mocked_bindings(sb_catboost_available = function() TRUE)
  nonpositive <- d$time
  nonpositive[[1L]] <- 0
  bad_time <- capture_sb_error(adapter$fit(x, nonpositive, d$event))
  no_event <- capture_sb_error(adapter$fit(x, d$time, rep(0L, nrow(d))))
  expect_s3_class(bad_time, "stormbreaker_catboost_aft_nonpositive_time")
  expect_s3_class(no_event, "stormbreaker_catboost_aft_no_event")
})

test_that("CatBoost AFT runs audited survival CV and final prediction", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_fixture(n = 60L, seed = 8821L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )
  pipeline <- sb_pipeline(
    learner = "catboost_aft_survival",
    catboost_iterations = 6L,
    catboost_learning_rate = 0.1,
    catboost_depth = 3L,
    catboost_l2_leaf_reg = 3,
    catboost_seed = 8822L
  )
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 8823L)

  expect_identical(pipeline$task_type, "survival")
  expect_identical(pipeline$risk_type, "aft_log_time")
  expect_identical(pipeline$risk_reference, "catboost_aft_normal_native")
  expect_identical(
    pipeline$risk_formula,
    "negative_catboost_aft_log_time_margin"
  )
  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(result$complete_survival_probability_oof, TRUE)
  expect_length(result$tree_models, 3L)
  expect_identical(
    vapply(result$tree_models, `[[`, "", "status"),
    rep("success", 3L)
  )
  expect_identical(
    vapply(result$baseline_hazards, `[[`, "", "status"),
    rep("not_applicable", 3L)
  )
  expect_identical(
    vapply(result$baseline_hazards, `[[`, "", "source"),
    rep("catboost_aft_normal_native", 3L)
  )
  expect_identical(
    result$provenance$catboost_version,
    as.character(utils::packageVersion("catboost"))
  )

  fit <- sb_fit_final(task, pipeline)
  rows <- c(5L, 21L, 49L)
  new_data <- d[rows, c("patient_id", rev(predictors)), drop = FALSE]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  direct <- sb_predict_catboost_aft_survival_learner_state(
    fit$learner_state,
    new_data[fit$required_features]
  )

  expect_identical(fit$status, "success")
  expect_identical(
    fit$learner_state$learner,
    "catboost_aft_survival"
  )
  expect_identical(prediction$status, "success")
  expect_equal(prediction$predictions$risk, direct$risk, tolerance = 0)
  expect_equal(
    prediction$predictions$survival_probability,
    direct$survival_probability,
    tolerance = 0
  )
})

# ---- catboost-binary-fixed-adapter ----

test_that("fixed CatBoost raw model reconstructs class-one probabilities", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

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
  adapter <- sb_catboost_binary_learner_adapter(
    iterations = 12L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 2.5,
    seed = 101L
  )
  model <- adapter$fit(x, outcome)
  observed <- adapter$predict(model, new_x)
  audit <- adapter$extract_tree_model(model)
  engine_params <- catboost::catboost.get_model_params(model$booster)
  expected_flat_params <- list(
    iterations = 12L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 2.5,
    random_seed = 101L,
    loss_function = "Logloss",
    eval_metric = "Logloss",
    task_type = "CPU",
    thread_count = 1L,
    bootstrap_type = "Bayesian",
    bagging_temperature = 1L,
    random_strength = 1L,
    rsm = 1L,
    grow_policy = "SymmetricTree",
    boosting_type = "Plain",
    use_best_model = FALSE,
    boost_from_average = TRUE,
    allow_writing_files = FALSE,
    save_snapshot = FALSE,
    logging_level = "Silent"
  )
  restored <- unserialize(audit$model_raw)
  rebuilt_pool <- catboost::catboost.load_pool(as.matrix(new_x))
  rebuilt <- as.numeric(catboost::catboost.predict(
    restored,
    rebuilt_pool,
    prediction_type = "Probability",
    thread_count = 1L
  ))

  expect_s3_class(model, "stormbreaker_catboost_binary")
  expect_identical(
    engine_params$flat_params[names(expected_flat_params)],
    expected_flat_params
  )
  expect_identical(engine_params$boosting_options$iterations, 12L)
  expect_identical(engine_params$boosting_options$od_config$type, "None")
  expect_identical(model$booster$tree_count, 12L)
  expect_equal(observed, rebuilt, tolerance = 0)
  expect_true(all(is.finite(observed)))
  expect_true(all(observed >= 0 & observed <= 1))
  expect_identical(audit$iterations, 12L)
  expect_identical(audit$tree_count, 12L)
  expect_identical(audit$features, names(x))
  expect_identical(audit$params, adapter$params)
  expect_identical(audit$model_serialization, "R_xdr_v3_catboost.Model")
  expect_type(audit$model_raw, "raw")
  expect_gt(length(audit$model_raw), 0L)
  expect_named(
    audit$feature_importance,
    c("feature", "importance"),
    ignore.order = FALSE
  )
  expect_identical(audit$feature_importance$feature, names(x))
  expect_identical(
    audit$feature_importance$importance[
      audit$feature_importance$feature == "unused"
    ],
    0
  )
})

test_that("fixed CatBoost predictions do not depend on assessment batching", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_binary_fixture(n = 50L, seed = 4201L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_catboost_binary_learner_adapter(
    iterations = 10L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 3,
    seed = 4301L
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

test_that("fixed CatBoost is deterministic and locks predictor schema", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_binary_fixture(n = 40L, seed = 4401L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_catboost_binary_learner_adapter(
    iterations = 8L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 3,
    seed = 4501L
  )
  first <- adapter$fit(x, d$outcome)
  second <- adapter$fit(x, d$outcome)

  expect_equal(
    adapter$predict(first, x),
    adapter$predict(second, x),
    tolerance = 0
  )
  expect_equal(
    adapter$extract_tree_model(first)$feature_importance,
    adapter$extract_tree_model(second)$feature_importance,
    tolerance = 0
  )

  schema_error <- capture_sb_error(
    adapter$predict(first, x[c("gene_b", "gene_a", "gene_c")])
  )
  expect_s3_class(schema_error, "stormbreaker_learner_schema_mismatch")
})

test_that("fixed CatBoost dependency and class failures remain explicit", {
  x <- data.frame(
    gene_a = seq_len(8),
    gene_b = rev(seq_len(8))
  )
  adapter <- sb_catboost_binary_learner_adapter(
    iterations = 5L,
    learning_rate = 0.1,
    depth = 2L,
    l2_leaf_reg = 3,
    seed = 1L
  )

  local_mocked_bindings(sb_catboost_available = function() FALSE)
  missing_error <- capture_sb_error(
    adapter$fit(x, rep(c(0L, 1L), each = 4L))
  )
  expect_s3_class(missing_error, "stormbreaker_missing_catboost")

  local_mocked_bindings(sb_catboost_available = function() TRUE)
  class_error <- capture_sb_error(adapter$fit(x, rep(0L, 8L)))
  expect_s3_class(class_error, "stormbreaker_catboost_single_class")
})

# ---- catboost-cox-fixed-adapter ----

test_that("fixed CatBoost Cox state reconstructs risk and horizon probability", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_fixture(n = 48L, seed = 7501L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(3L, 17L, 29L), , drop = FALSE]
  adapter <- sb_catboost_cox_learner_adapter(
    iterations = 12L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 2.5,
    seed = 7502L
  )
  model <- adapter$fit(x, d$time, d$event)
  observed_risk <- adapter$predict(model, new_x)
  probability <- adapter$predict_survival(
    model,
    x,
    d$time,
    d$event,
    new_x,
    prediction_horizon = 7
  )
  audit <- adapter$extract_tree_model(model)

  restored <- unserialize(audit$model_raw)
  rebuilt_pool <- catboost::catboost.load_pool(as.matrix(new_x))
  rebuilt_risk <- as.numeric(catboost::catboost.predict(
    restored,
    rebuilt_pool,
    prediction_type = "RawFormulaVal",
    ntree_start = 0L,
    ntree_end = 12L,
    thread_count = 1L
  ))
  analysis_risk <- adapter$predict(model, x)
  expected_hazard <- sb_breslow_baseline_hazard(
    d$time,
    d$event,
    analysis_risk,
    prediction_horizon = 7
  )
  expected_probability <- sb_survival_probability_from_hazard(
    observed_risk,
    expected_hazard
  )
  state <- adapter$extract_state(model, list(
    x = x,
    time = d$time,
    event = d$event,
    prediction_horizon = 7,
    time_unit = "month"
  ))
  state_prediction <- adapter$predict_state(state, new_x)

  expect_s3_class(model, "stormbreaker_catboost_cox")
  expect_equal(observed_risk, rebuilt_risk, tolerance = 0)
  expect_identical(model$booster$tree_count, 12L)
  expect_identical(adapter$params$loss_function, "Cox")
  expect_identical(adapter$params$eval_metric, "Cox")
  expect_identical(probability$prediction_horizon, 7)
  expect_identical(probability$cumulative_baseline_hazard, expected_hazard)
  expect_identical(probability$source, "catboost_cox_breslow_raw_margin")
  expect_identical(probability$ties, "breslow")
  expect_equal(probability$probability, expected_probability, tolerance = 0)
  expect_identical(audit$features, names(x))
  expect_identical(audit$params, adapter$params)
  expect_identical(audit$model_serialization, "R_xdr_v3_catboost.Model")
  expect_type(audit$model_raw, "raw")
  expect_gt(length(audit$model_raw), 0L)
  expect_identical(state$learner, "catboost_cox")
  expect_false("booster" %in% names(state))
  expect_equal(state_prediction$risk, observed_risk, tolerance = 0)
  expect_equal(
    state_prediction$survival_probability,
    expected_probability,
    tolerance = 0
  )
  expect_identical(state_prediction$prediction_horizon, 7)
  expect_identical(state_prediction$time_unit, "month")
})

test_that("fixed CatBoost Cox is deterministic and batching independent", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_fixture(n = 50L, seed = 7511L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_catboost_cox_learner_adapter(
    iterations = 10L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 3,
    seed = 7512L
  )
  first <- adapter$fit(x[1:40, , drop = FALSE], d$time[1:40], d$event[1:40])
  second <- adapter$fit(x[1:40, , drop = FALSE], d$time[1:40], d$event[1:40])
  new_x <- x[41:50, , drop = FALSE]

  batch <- adapter$predict(first, new_x)
  repeated <- adapter$predict(second, new_x)
  one_at_a_time <- vapply(seq_len(nrow(new_x)), function(row) {
    adapter$predict(first, new_x[row, , drop = FALSE])
  }, numeric(1))

  expect_equal(batch, repeated, tolerance = 0)
  expect_equal(batch, one_at_a_time, tolerance = 0)
  expect_error(
    adapter$predict(first, new_x[c("gene_b", "gene_a", "gene_c")]),
    class = "stormbreaker_learner_schema_mismatch"
  )
})

test_that("fixed CatBoost Cox failures remain explicit", {
  d <- make_sb_fixture(n = 24L, seed = 7521L)
  x <- d[c("gene_a", "gene_b")]
  adapter <- sb_catboost_cox_learner_adapter(
    iterations = 5L,
    learning_rate = 0.1,
    depth = 2L,
    l2_leaf_reg = 3,
    seed = 7522L
  )

  local_mocked_bindings(sb_catboost_available = function() FALSE)
  expect_error(
    adapter$fit(x, d$time, d$event),
    class = "stormbreaker_missing_catboost"
  )

  local_mocked_bindings(sb_catboost_available = function() TRUE)
  nonpositive <- d$time
  nonpositive[[1L]] <- 0
  expect_error(
    adapter$fit(x, nonpositive, d$event),
    class = "stormbreaker_catboost_cox_nonpositive_time"
  )
  expect_error(
    adapter$fit(x, d$time, rep(0L, nrow(d))),
    class = "stormbreaker_catboost_cox_no_event"
  )
})

test_that("CatBoost Cox runs audited survival CV and final prediction", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_fixture(n = 72L, seed = 7531L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )
  pipeline <- sb_pipeline(
    learner = "catboost_cox",
    catboost_iterations = 8L,
    catboost_learning_rate = 0.1,
    catboost_depth = 3L,
    catboost_l2_leaf_reg = 3,
    catboost_seed = 7532L
  )
  screened <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "catboost_cox",
    selector_n_features = 2L,
    catboost_iterations = 8L,
    catboost_depth = 3L,
    catboost_seed = 7533L
  )
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 7534L)
  screened_result <- sb_run_cv(task, screened, folds = 3L, seed = 7534L)

  expect_identical(pipeline$task_type, "survival")
  expect_identical(pipeline$risk_type, "cox_raw_margin")
  expect_identical(pipeline$risk_reference, "catboost_native")
  expect_identical(pipeline$risk_formula, "catboost_raw_formula_value")
  expect_identical(result$status, "success")
  expect_true(result$complete_oof)
  expect_true(result$complete_survival_probability_oof)
  expect_true(all(is.finite(result$predictions$risk)))
  expect_true(all(is.finite(result$predictions$survival_probability)))
  expect_identical(unique(result$predictions$learner), "catboost_cox")
  expect_identical(result$metrics$harrell_c$status, "success")
  expect_identical(result$metrics$ipcw_brier$status, "success")
  expect_length(result$tree_models, 3L)
  expect_true(all(vapply(
    result$tree_models,
    function(audit) identical(audit$status, "success") &&
      is.raw(audit$model_raw) && length(audit$model_raw) > 0L,
    logical(1)
  )))
  expect_true(all(vapply(
    result$baseline_hazards,
    function(audit) identical(
      audit$source,
      "catboost_cox_breslow_raw_margin"
    ),
    logical(1)
  )))
  expect_identical(
    result$provenance$catboost_version,
    as.character(utils::packageVersion("catboost"))
  )
  expect_identical(screened_result$status, "success")
  expect_length(screened_result$selections, 3L)
  expect_true(all(vapply(
    screened_result$selections,
    function(selection) length(selection$selected_features) == 2L,
    logical(1)
  )))

  fit <- sb_fit_final(task, pipeline)
  rows <- c(5L, 21L, 49L)
  new_data <- d[rows, c("patient_id", rev(predictors)), drop = FALSE]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  selected_x <- new_data[fit$required_features]
  direct <- sb_predict_catboost_cox_learner_state(
    fit$learner_state,
    selected_x
  )

  expect_identical(fit$status, "success")
  expect_s3_class(fit$learner_state, "stormbreaker_learner_state")
  expect_identical(fit$learner_state$learner, "catboost_cox")
  expect_type(fit$learner_state$model_raw, "raw")
  expect_identical(prediction$status, "success")
  expect_equal(prediction$predictions$risk, direct$risk, tolerance = 0)
  expect_equal(
    prediction$predictions$survival_probability,
    direct$survival_probability,
    tolerance = 0
  )
  expect_identical(
    fit$provenance$catboost_version,
    as.character(utils::packageVersion("catboost"))
  )
})
})

# ---- model-gbm ----
local({
# ---- algorithm-gbm-xgboost ----

test_that("all nine XGBoost variants preserve native scores and zero-use features", {
  skip_if_not_installed("xgboost")
  for (id in sb_xgboost_learner_ids()) {
    binary <- id %in% sb_xgboost_binary_learner_ids()
    d <- if (binary) make_sb_binary_fixture(n = 100L, seed = 108L) else make_sb_fixture(n = 100L, seed = 8612L)
    d$constant <- 0
    features <- c("gene_a", "gene_b", "gene_c", "constant")
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    original <- fit$learner_state
    q <- sb_boosting_diagnostics(fit, "importance")
    native_model <- sb_load_xgboost_state_model(fit$learner_state, id)$booster
    native <- as.data.frame(xgboost::xgb.importance(feature_names = features, model = native_model))
    expect_identical(q$native_table, native)
    expect_identical(q$table$feature, features)
    score <- if (grepl("linear", id, fixed = TRUE)) "Weight" else "Gain"
    expect_equal(q$table$value[match(native$Feature, features)], native[[score]], tolerance = 0)
    expect_equal(q$table$value[q$table$feature == "constant"], 0)
    expect_identical(fit$learner_state, original)
    err <- tryCatch(sb_boosting_diagnostics(fit, "training_loss"), error = identity)
    expect_s3_class(err, "stormbreaker_unavailable_diagnostic")
    if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(
      ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
  }
})

test_that("GBM final refits preserve the selected count and native predictions", {
  skip_if_not_installed("gbm")
  for (signal in c(FALSE, TRUE)) {
    d <- make_sb_fixture(n = 100L, seed = 8612L)
    if (signal) {
      d$time <- exp(2 - d$gene_a)
      d$event[] <- 1L
    }
    x <- d[c("gene_a", "gene_b", "gene_c")]
    params <- sb_gbm_cox_cv_params(inner_folds = 3L, maximum_trees = 30L)
    model <- sb_fit_gbm_cox_cv(x, d$time, d$event, params)
    count <- model$n_trees
    expect_equal(count, which.min(model$tuning$cv_error))
    expect_equal(count == 1L, !signal)
    formula_data <- sb_additional_survival_formula_data(x, d$time, d$event)
    native <- withr::with_seed(params$seed, gbm::gbm(
      formula = formula_data$formula, data = formula_data$data,
      distribution = params$distribution, n.trees = count,
      interaction.depth = params$interaction_depth, n.minobsinnode = params$min_obs_node,
      shrinkage = params$shrinkage, bag.fraction = params$bag_fraction,
      cv.folds = if (count == 1L) 0L else params$inner_folds,
      n.cores = 1L, keep.data = FALSE, verbose = FALSE))
    expect_equal(sb_predict_gbm_cox(model, x),
      unname(as.numeric(stats::predict(native, newdata = x, n.trees = count, type = "link"))), tolerance = 0)
  }
})

test_that("GBM paths distinguish native final fits and original CV including a one-tree optimum", {
  skip_if_not_installed("gbm")
  for (id in sb_gbm_diagnostic_ids()) {
    binary <- id == "gbm_binary"
    d <- if (binary) make_sb_binary_fixture(n = 100L, seed = 108L) else make_sb_fixture(n = 100L, seed = 8612L)
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    pipeline <- if (id == "gbm_cox_cv") sb_pipeline(learner = id, gbm_cv_max_trees = 30L, learner_tuning_folds = 3L) else sb_pipeline(learner = id)
    fit <- sb_fit_final(task, pipeline, capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    original <- fit$learner_state
    engine <- if (binary) fit$diagnostic_artifacts$learner$engine else fit$diagnostic_artifacts$learner$booster
    q <- sb_boosting_diagnostics(fit)
    expect_equal(q$table$loss, engine$train.error, tolerance = 0)
    expect_equal(nrow(q$table), engine$n.trees)
    importance <- sb_boosting_diagnostics(fit, "importance")
    expect_equal(importance$table$value, unname(gbm::relative.influence(engine,
      n.trees = engine$n.trees, scale. = FALSE, sort. = FALSE)))
    if (id == "gbm_cox_cv") {
      cv <- sb_boosting_diagnostics(fit, "cv_curve")
      expect_equal(cv$table$loss, fit$diagnostic_artifacts$learner$tuning$cv_error)
      expect_equal(cv$selected_iteration, which.min(cv$table$loss))
      expect_equal(cv$selected_iteration, 1L)
      expect_equal(nrow(cv$table), 30L)
      expect_equal(nrow(q$table), 1L)
      expect_null(engine$cv.error)
      if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(
        ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(cv)), "ggplot_built")
    }
    expect_identical(fit$learner_state, original)
    if (requireNamespace("ggplot2", quietly = TRUE)) for (z in list(q, importance)) expect_s3_class(
      ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    fit$diagnostic_artifacts <- NULL
    expect_identical(sb_boosting_diagnostics(fit, "importance"), importance)
    if (binary) {
      err <- tryCatch(sb_boosting_diagnostics(fit), error = identity)
      expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
    } else expect_identical(sb_boosting_diagnostics(fit), q)
  }
})

# ---- gbm-cox-fixed-adapter ----

test_that("fixed GBM-Cox state reconstructs risk and horizon probability", {
  skip_if_not_installed("gbm", minimum_version = "2.2.0")

  d <- make_sb_fixture(n = 54L, seed = 7701L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(4L, 19L, 37L), , drop = FALSE]
  adapter <- sb_gbm_cox_learner_adapter(
    n_trees = 12L,
    interaction_depth = 2L,
    min_obs_node = 3L,
    shrinkage = 0.05,
    bag_fraction = 1,
    seed = 7702L
  )
  model <- adapter$fit(x, d$time, d$event)
  observed_risk <- adapter$predict(model, new_x)
  probability <- adapter$predict_survival(
    model,
    x,
    d$time,
    d$event,
    new_x,
    prediction_horizon = 7
  )
  audit <- adapter$extract_tree_model(model)

  restored <- unserialize(audit$model_raw)
  restored_risk <- as.numeric(stats::predict(
    restored,
    newdata = as.matrix(new_x),
    n.trees = 12L,
    type = "link"
  ))
  analysis_risk <- adapter$predict(model, x)
  expected_hazard <- sb_breslow_baseline_hazard(
    d$time,
    d$event,
    analysis_risk,
    prediction_horizon = 7
  )
  expected_probability <- sb_survival_probability_from_hazard(
    observed_risk,
    expected_hazard
  )
  state <- adapter$extract_state(model, list(
    x = x,
    time = d$time,
    event = d$event,
    prediction_horizon = 7,
    time_unit = "month"
  ))
  state_prediction <- adapter$predict_state(state, new_x)

  expect_s3_class(model, "stormbreaker_gbm_cox")
  expect_identical(as.integer(model$booster$n.trees), 12L)
  expect_null(model$booster$data)
  expect_equal(observed_risk, restored_risk, tolerance = 0)
  expect_identical(probability$prediction_horizon, 7)
  expect_identical(probability$cumulative_baseline_hazard, expected_hazard)
  expect_identical(
    probability$source,
    "gbm_cox_breslow_log_hazard"
  )
  expect_identical(probability$ties, "breslow")
  expect_equal(probability$probability, expected_probability, tolerance = 0)
  expect_identical(audit$n_trees, 12L)
  expect_identical(audit$features, names(x))
  expect_identical(audit$params, adapter$params)
  expect_identical(audit$importance_type, "relative_influence")
  expect_identical(audit$model_serialization, "R_xdr_v3_gbm")
  expect_type(audit$model_raw, "raw")
  expect_gt(length(audit$model_raw), 0L)
  expect_identical(audit$feature_importance$feature, names(x))
  expect_true(all(is.finite(audit$feature_importance$importance)))
  expect_identical(state$learner, "gbm_cox")
  expect_false("booster" %in% names(state))
  expect_equal(state_prediction$risk, observed_risk, tolerance = 0)
  expect_equal(
    state_prediction$survival_probability,
    expected_probability,
    tolerance = 0
  )
  expect_identical(state_prediction$prediction_horizon, 7)
  expect_identical(state_prediction$time_unit, "month")
})

test_that("fixed GBM-Cox settings are explicit and reproducible", {
  skip_if_not_installed("gbm", minimum_version = "2.2.0")

  pipeline <- sb_pipeline(learner = "gbm_cox")
  expect_identical(
    pipeline$id,
    paste0(
      "none__gbm_cox__trees100__depth1__minnode10",
      "__shrink0.1__bag1__seed1"
    )
  )
  expect_identical(pipeline$task_type, "survival")
  expect_identical(pipeline$risk_type, "cox_log_hazard")
  expect_identical(pipeline$risk_reference, "gbm_native")
  expect_identical(pipeline$risk_formula, "gbm_link_prediction")
  expect_identical(pipeline$learner_params$distribution, "coxph")
  expect_identical(pipeline$learner_params$keep_data, FALSE)

  d <- make_sb_fixture(n = 50L, seed = 7711L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_gbm_cox_learner_adapter(
    n_trees = 10L,
    interaction_depth = 1L,
    min_obs_node = 3L,
    shrinkage = 0.1,
    bag_fraction = 0.8,
    seed = 7712L
  )
  first <- adapter$fit(x[1:40, , drop = FALSE], d$time[1:40], d$event[1:40])
  second <- adapter$fit(x[1:40, , drop = FALSE], d$time[1:40], d$event[1:40])
  new_x <- x[41:50, , drop = FALSE]

  expect_equal(
    adapter$predict(first, new_x),
    adapter$predict(second, new_x),
    tolerance = 0
  )
  expect_error(
    adapter$predict(first, new_x[c("gene_b", "gene_a", "gene_c")]),
    class = "stormbreaker_learner_schema_mismatch"
  )

  expect_error(
    sb_pipeline(learner = "gbm_cox", gbm_shrinkage = 0),
    class = "stormbreaker_bad_gbm_shrinkage"
  )
  expect_error(
    sb_pipeline(learner = "gbm_cox", gbm_bag_fraction = 1.1),
    class = "stormbreaker_bad_gbm_bag_fraction"
  )
  expect_error(
    sb_pipeline(learner = "coxph", gbm_n_trees = 10L),
    class = "stormbreaker_irrelevant_gbm_parameter"
  )
})

test_that("fixed GBM-Cox failures remain explicit", {
  d <- make_sb_fixture(n = 24L, seed = 7721L)
  x <- d[c("gene_a", "gene_b")]
  adapter <- sb_gbm_cox_learner_adapter(
    n_trees = 5L,
    min_obs_node = 3L,
    seed = 7722L
  )

  local_mocked_bindings(sb_gbm_available = function() FALSE)
  expect_error(
    adapter$fit(x, d$time, d$event),
    class = "stormbreaker_missing_gbm"
  )
  expect_error(
    adapter$fit(x, d$time, rep(0L, nrow(d))),
    class = "stormbreaker_gbm_cox_no_event"
  )
})

test_that("GBM-Cox runs selected survival CV and final prediction", {
  skip_if_not_installed("gbm", minimum_version = "2.2.0")

  d <- make_sb_fixture(n = 72L, seed = 7731L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "gbm_cox",
    selector_n_features = 2L,
    gbm_n_trees = 10L,
    gbm_min_obs_node = 3L,
    gbm_seed = 7732L
  )
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 7733L)

  expect_identical(result$status, "success")
  expect_true(result$complete_oof)
  expect_true(result$complete_survival_probability_oof)
  expect_true(all(is.finite(result$predictions$risk)))
  expect_true(all(is.finite(result$predictions$survival_probability)))
  expect_identical(unique(result$predictions$learner), "gbm_cox")
  expect_length(result$selections, 3L)
  expect_true(all(vapply(
    result$selections,
    function(selection) length(selection$selected_features) == 2L,
    logical(1)
  )))
  expect_length(result$tree_models, 3L)
  expect_true(all(vapply(
    result$tree_models,
    function(audit) identical(audit$status, "success") &&
      is.raw(audit$model_raw) && length(audit$model_raw) > 0L,
    logical(1)
  )))
  expect_true(all(vapply(
    result$baseline_hazards,
    function(audit) identical(
      audit$source,
      "gbm_cox_breslow_log_hazard"
    ),
    logical(1)
  )))
  expect_identical(
    result$provenance$gbm_version,
    as.character(utils::packageVersion("gbm"))
  )

  fit <- sb_fit_final(task, pipeline)
  rows <- c(5L, 21L, 49L)
  new_data <- d[rows, c("patient_id", rev(predictors)), drop = FALSE]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  direct <- sb_predict_gbm_cox_learner_state(
    fit$learner_state,
    new_data[fit$required_features]
  )

  expect_identical(fit$status, "success")
  expect_identical(fit$learner_state$learner, "gbm_cox")
  expect_type(fit$learner_state$model_raw, "raw")
  expect_identical(prediction$status, "success")
  expect_equal(prediction$predictions$risk, direct$risk, tolerance = 0)
  expect_equal(
    prediction$predictions$survival_probability,
    direct$survival_probability,
    tolerance = 0
  )
  expect_identical(
    fit$provenance$gbm_version,
    as.character(utils::packageVersion("gbm"))
  )
})

# ---- gbm-cv-reuse ----

test_that("GBM CV reuse preserves the original two-call native fit", {
  skip_if_not_installed("gbm")
  d <- sb_example_survival(90, seed = 61)
  x <- d[c("gene_a", "gene_b")]
  for (steps in c(1L, 2L, 100L)) {
    p <- sb_gbm_cox_cv_params(inner_folds = 3L, maximum_trees = steps)
    p$shrinkage <- 0.1
    data <- sb_additional_survival_formula_data(x, d$time, d$event)
    native <- function(nt, folds = p$inner_folds) {
      withr::with_seed(p$seed, gbm::gbm(data$formula, data = data$data,
        distribution = p$distribution, n.trees = nt,
        interaction.depth = p$interaction_depth, n.minobsinnode = p$min_obs_node,
        shrinkage = p$shrinkage, bag.fraction = p$bag_fraction,
        cv.folds = folds, n.cores = 1L, keep.data = FALSE, verbose = FALSE))
    }
    withr::local_seed(91)
    before <- .Random.seed
    invisible(capture.output({
      search <- if (steps == 1L) NULL else native(steps)
      best <- if (steps == 1L) 1L else which.min(search$cv.error)
      expected <- native(best, if (best == 1L) 0L else p$inner_folds)
      actual <- sb_fit_gbm_cox_cv(x, d$time, d$event, p)
    }))
    expect_identical(.Random.seed, before)
    expect_identical(actual$n_trees, best)
    expect_identical(actual$tuning$cv_error, as.numeric(search$cv.error))
    expect_identical(actual$booster$trees, expected$trees)
    expect_identical(actual$booster$cv.error, expected$cv.error)
    expect_identical(actual$booster$cv.fitted, expected$cv.fitted)
    expect_identical(actual$booster$train.error, expected$train.error)
    expect_identical(actual$booster$oobag.improve, expected$oobag.improve)
    expect_identical(sb_predict_gbm_cox(actual, x),
      as.numeric(predict(expected, data$data, n.trees = best, type = "link")))
    restored <- unserialize(serialize(actual, NULL))
    expect_identical(sb_predict_gbm_cox(restored, x),
      sb_predict_gbm_cox(actual, x))
  }
})

test_that("warned GBM searches are not reused", {
  skip_if_not_installed("gbm")
  calls <- 0L
  local_mocked_bindings(gbm = function(...) {
    calls <<- calls + 1L
    warning("native search warning")
    structure(list(n.trees = 2L, cv.error = c(2, 1),
      call = quote(gbm()), m = quote(gbm())), class = "gbm")
  }, .package = "gbm")
  d <- sb_example_survival(60, seed = 1)
  seen <- character()
  withCallingHandlers(sb_fit_gbm_cox_cv(d[c("gene_a", "gene_b")], d$time,
    d$event, sb_gbm_cox_cv_params(inner_folds = 3L, maximum_trees = 2L)),
    warning = function(w) { seen <<- c(seen, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_identical(calls, 2L)
  expect_identical(seen, rep("native search warning", 2L))
})

test_that("an unremarkable upper-bound GBM optimum needs only one native call", {
  skip_if_not_installed("gbm")
  calls <- 0L
  local_mocked_bindings(gbm = function(...) {
    calls <<- calls + 1L
    structure(list(n.trees = 2L, cv.error = c(2, 1),
      call = quote(gbm()), m = quote(gbm())), class = "gbm")
  }, .package = "gbm")
  d <- sb_example_survival(60, seed = 1)
  fit <- sb_fit_gbm_cox_cv(d[c("gene_a", "gene_b")], d$time, d$event,
    sb_gbm_cox_cv_params(inner_folds = 3L, maximum_trees = 2L))
  expect_identical(calls, 1L)
  expect_identical(fit$n_trees, 2L)
})

# ---- gbm-single-candidate ----

test_that("a single allowed GBM tree is fitted without inventing a CV search", {
  skip_if_not_installed("gbm")
  d <- make_sb_fixture(n = 100L, seed = 8612L)
  features <- c("gene_a", "gene_b", "gene_c")
  x <- d[features]
  params <- sb_gbm_cox_cv_params(inner_folds = 3L, maximum_trees = 1L)
  model <- sb_fit_gbm_cox_cv(x, d$time, d$event, params)
  expect_identical(model$n_trees, 1L)
  expect_identical(model$tuning$cv_error, numeric())
  expect_identical(model$tuning$search_status, "single_candidate_no_search")
  expect_null(model$booster$cv.error)
  formula_data <- sb_additional_survival_formula_data(x, d$time, d$event)
  native <- withr::with_seed(params$seed, gbm::gbm(
    formula = formula_data$formula, data = formula_data$data,
    distribution = params$distribution, n.trees = 1L,
    interaction.depth = params$interaction_depth, n.minobsinnode = params$min_obs_node,
    shrinkage = params$shrinkage, bag.fraction = params$bag_fraction,
    cv.folds = 0L, n.cores = 1L, keep.data = FALSE, verbose = FALSE))
  expect_equal(sb_predict_gbm_cox(model, x),
    unname(as.numeric(stats::predict(native, newdata = x, n.trees = 1L, type = "link"))), tolerance = 0)
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  fit <- sb_fit_final(task, sb_pipeline(learner = "gbm_cox_cv", gbm_cv_max_trees = 1L,
    learner_tuning_folds = 3L), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  expect_identical(fit$learner_state$tuning$search_status, "single_candidate_no_search")
  err <- tryCatch(sb_boosting_diagnostics(fit, "cv_curve"), error = identity)
  expect_s3_class(err, "stormbreaker_unavailable_diagnostic")
  loss <- sb_boosting_diagnostics(fit, "training_loss")
  expect_equal(loss$table$loss, fit$diagnostic_artifacts$learner$booster$train.error)
  expect_equal(nrow(loss$table), 1L)
  expect_identical(serialize(fit, NULL, version = 2L), before)
  pipeline <- fit$pipeline
  compact <- sb_fit_final(task, pipeline)
  expect_identical(compact$status, "success")
  expect_identical(compact$learner_state, fit$learner_state)
  cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 7733L)
  expect_identical(cv$status, "success")
  expect_true(cv$complete_oof)
  expect_true(all(is.finite(cv$predictions$risk)))
})
})

# ---- model-lightgbm ----
local({
test_that("all LightGBM variants keep native gain and count actual rather than configured trees", {
  skip_if_not_installed("lightgbm")
  skip_if_not_installed("jsonlite")
  d <- make_sb_binary_fixture(n = 100L, seed = 108L)
  d$constant <- 0
  features <- c("gene_a", "gene_b", "gene_c", "constant")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  for (id in sb_lightgbm_binary_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    original <- fit$learner_state
    importance <- sb_boosting_diagnostics(fit, "importance")
    expect_identical(importance$native_table, original$feature_importance)
    expect_equal(importance$table$value, original$feature_importance$gain)
    expect_equal(importance$table$value[importance$table$feature == "constant"], 0)
    sizes <- sb_boosting_diagnostics(fit, "tree_size")
    # Native model text independently exposes a num_leaves line for each tree.
    lines <- strsplit(original$model_string, "\n", fixed = TRUE)[[1]]
    native_leaves <- as.integer(sub("num_leaves=", "", grep("^num_leaves=", lines, value = TRUE), fixed = TRUE))
    expect_equal(sizes$table$terminal_nodes, native_leaves)
    expect_equal(sizes$table$native_tree_index, seq_along(native_leaves) - 1L)
    expect_equal(sizes$requested_iterations, original$nrounds)
    if (id == "lightgbm_binary_goss") expect_lt(nrow(sizes$table), sizes$requested_iterations)
    expect_identical(fit$learner_state, original)
    err <- tryCatch(sb_boosting_diagnostics(fit, "training_loss"), error = identity)
    expect_s3_class(err, "stormbreaker_unavailable_diagnostic")
    if (requireNamespace("ggplot2", quietly = TRUE)) for (q in list(importance, sizes)) expect_s3_class(
      ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
  }
})

test_that("native empty-pair failures are explicit and preserve their original cause", {
  skip_if_not_installed("catboost")
  d <- make_sb_binary_fixture(n = 100L, seed = 108L)
  fit <- sb_fit_final(
    sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b")),
    sb_pipeline(learner = "catboost_binary", catboost_depth = 1L)
  )
  expect_identical(fit$status, "success")
  original <- fit$learner_state
  native <- tryCatch(catboost::catboost.get_feature_importance(
    unserialize(original$model_raw), type = "Interaction", thread_count = 1L
  ), error = identity)
  observed <- tryCatch(sb_boosting_diagnostics(fit, "interaction"), error = identity)
  if (inherits(native, "error")) {
    expect_s3_class(observed, "stormbreaker_unavailable_diagnostic")
    expect_identical(conditionMessage(observed$parent), conditionMessage(native))
  } else {
    # Future CatBoost releases may fix the native empty-result wrapper.
    expect_s3_class(observed, "stormbreaker_algorithm_diagnostic")
    expect_equal(nrow(observed$table), 0L)
  }
  expect_identical(fit$learner_state, original)
})

test_that("all CatBoost task types map native interaction indices and limit only display", {
  skip_if_not_installed("catboost")
  for (id in sb_catboost_diagnostic_ids()) {
    binary <- id == "catboost_binary"
    d <- if (binary) make_sb_binary_fixture(n = 100L, seed = 108L) else make_sb_fixture(n = 100L, seed = 8612L)
    features <- c("gene_c", "gene_a", "gene_b")
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    original <- fit$learner_state
    importance <- sb_boosting_diagnostics(fit, "importance")
    expect_identical(importance$native_table, original$feature_importance)
    q <- sb_boosting_diagnostics(fit, "interaction", max_pairs = 1L)
    native <- catboost::catboost.get_feature_importance(unserialize(original$model_raw), type = "Interaction", thread_count = 1L)
    expect_identical(q$native_table, native)
    expect_equal(q$table$feature_1, features[native[, 1] + 1L])
    expect_equal(q$table$feature_2, features[native[, 2] + 1L])
    expect_equal(q$table$score, unname(native[, 3]))
    expect_identical(fit$learner_state, original)
    expect_equal(nrow(q$table), nrow(native))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- sb_plot_algorithm_diagnostic(q)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_equal(nrow(p$data), 1L)
      expect_equal(p$data$score, max(native[, 3]))
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(importance)), "ggplot_built")
      empty <- q
      empty$table <- empty$table[0, ]
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(empty)), "ggplot_built")
    }
    err <- tryCatch(sb_boosting_diagnostics(fit, "interaction", max_pairs = 0L), error = identity)
    expect_s3_class(err, "error")
  }
})

test_that("one-feature CatBoost fits have no invented interaction pairs", {
  skip_if_not_installed("catboost")
  for (id in sb_catboost_diagnostic_ids()) {
  binary <- id == "catboost_binary"
  d <- if (binary) make_sb_binary_fixture(n = 100L, seed = 108L) else make_sb_fixture(n = 100L, seed = 8612L)
  task <- if (binary) sb_binary_task(d, "patient_id", "outcome", "gene_a") else sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month")
  fit <- sb_fit_final(task, sb_pipeline(learner = id))
  expect_identical(fit$status, "success")
  original <- fit$learner_state
  q <- sb_boosting_diagnostics(fit, "interaction")
  expect_equal(nrow(q$table), 0L)
  expect_equal(dim(q$native_table), c(0L, 3L))
  expect_identical(fit$learner_state, original)
  expect_match(q$interpretation, "extractor was not called", fixed = TRUE)
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(
    ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
  }
})
})

# ---- model-mboost ----
local({
# ---- algorithm-boosting ----

test_that("all twenty mboost options retain original loss, scale, stopping and components", {
  skip_if_not_installed("mboost")
  skip_if_not_installed("partykit")
  ids <- sb_mboost_diagnostic_ids()
  expect_length(ids, 20L)
  for (id in ids) {
    binary <- grepl("binary", id, fixed = TRUE)
    d <- if (binary) make_sb_binary_fixture(n = 100L, seed = 108L) else make_sb_fixture(n = 100L, seed = 8612L)
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    expect_null(fit$diagnostic_artifacts)
    original <- serialize(fit, NULL)
    engine <- unserialize(fit$learner_state$engine_raw)
    q <- sb_boosting_diagnostics(fit)
    expect_equal(q$table$loss, unname(mboost::risk(engine)), tolerance = 0)
    expect_equal(q$table$iteration, seq.int(0L, mboost::mstop(engine)))
    expect_equal(q$stopping_iteration, mboost::mstop(engine))
    if (grepl("unoc", id, fixed = TRUE)) expect_equal(sum(q$table$loss < 0), nrow(q$table))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    }
    if (!startsWith(id, "blackboost_")) {
      components <- sb_boosting_diagnostics(fit, "component_usage")
      native_selected <- as.integer(mboost::selected(engine))
      expect_identical(components$selections$component_id, native_selected)
      expect_equal(components$table$updates, tabulate(native_selected, nrow(components$table)))
      expect_equal(sum(components$table$updates), mboost::mstop(engine))
      expect_setequal(setdiff(components$table$feature, "(Intercept)"), features)
      expect_equal(components$selections$feature,
        components$table$feature[components$selections$component_id])
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        p <- sb_plot_algorithm_diagnostic(components)
        expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
        expect_equal(nrow(p$data), nrow(components$table))
      }
    } else {
      err <- tryCatch(sb_boosting_diagnostics(fit, "component_usage"), error = identity)
      expect_s3_class(err, "stormbreaker_unavailable_diagnostic")
    }
    expect_identical(serialize(fit, NULL), original)
    expect_equal(length(grepRaw("B001", serialize(q, NULL), fixed = TRUE)), 0L)
  }
})

# ---- mboost-aft-fixed ----

make_mboost_aft_fixture <- function(n = 84L, p = 6L, seed = 9811L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n)
    signal <- -0.8 * x[, 1L]
    if (p >= 2L) {
      signal <- signal + 0.5 * x[, 2L]
    }
    event_time <- exp(2.2 + signal + stats::rnorm(n, sd = 0.65))
    censor_time <- exp(stats::rnorm(n, mean = 2.55, sd = 0.8))
    feature_names <- c(
      "gene one",
      if (p >= 2L) "2gene" else NULL,
      if (p >= 3L) paste0("gene-", seq_len(p - 2L)) else NULL
    )
    colnames(x) <- feature_names
    data.frame(
      patient_id = sprintf("MBA%03d", seq_len(n)),
      time = pmin(event_time, censor_time),
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    )
  })
}

mboost_aft_expected_specifications <- function() {
  data.frame(
    learner = c(
      "gamboost_weibull_aft",
      "gamboost_loglogistic_aft",
      "gamboost_lognormal_aft",
      "glmboost_weibull_aft",
      "glmboost_loglogistic_aft",
      "glmboost_lognormal_aft",
      "blackboost_weibull_aft",
      "blackboost_loglogistic_aft",
      "blackboost_lognormal_aft"
    ),
    algorithm = rep(c("gamboost", "glmboost", "blackboost"), each = 3L),
    family = rep(c("Weibull", "Loglog", "Lognormal"), times = 3L),
    distribution = rep(
      c("weibull", "loglogistic", "lognormal"),
      times = 3L
    ),
    stringsAsFactors = FALSE
  )
}

mboost_aft_native_survival <- function(distribution, link, scale, horizon) {
  transformed <- (log(horizon) - link) / scale
  switch(distribution,
    weibull = exp(-exp(transformed)),
    loglogistic = stats::plogis(-transformed),
    lognormal = stats::pnorm(transformed, lower.tail = FALSE)
  )
}

test_that("mboost AFT declares nine distribution-by-base-learner methods", {
  expected <- mboost_aft_expected_specifications()
  expect_identical(sb_mboost_aft_learner_ids(), expected$learner)

  for (i in seq_len(nrow(expected))) {
    learner <- expected$learner[[i]]
    params <- sb_mboost_aft_params(learner)
    expect_identical(params$algorithm, expected$algorithm[[i]])
    expect_identical(params$family, expected$family[[i]])
    expect_identical(params$distribution, expected$distribution[[i]])
    expect_identical(params$mstop, 100L)
    expect_identical(params$nu, 0.1)
    expect_identical(params$risk, "inbag")
    expect_identical(params$nuisance_range, c(0.1, 100))
    if (identical(params$algorithm, "blackboost")) {
      expect_identical(params$seed, 1L)
    } else {
      expect_false("seed" %in% names(params))
    }
  }

  registry <- sb_pipeline_registry()
  rows <- registry[
    registry$learner %in% expected$learner, ,
    drop = FALSE
  ]
  expect_identical(nrow(rows), 9L * length(expected_fixed_learner_selectors("survival")))
  expect_identical(anyDuplicated(rows[c("selector", "learner")]), 0L)
  expect_setequal(
    unique(rows$selector),
    expected_fixed_learner_selectors("survival")
  )
  expect_true(all(grepl("mboost", rows$optional_dependency, fixed = TRUE)))
  black_rows <- rows[grepl("^blackboost", rows$learner), , drop = FALSE]
  expect_true(all(grepl("partykit", black_rows$optional_dependency, fixed = TRUE)))
})

test_that("mboost AFT risk, direct survival, RNG, and state match native fits", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  data <- make_mboost_aft_fixture(seed = 9812L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[c(3L, 27L, 68L), , drop = FALSE]
  horizon <- as.numeric(stats::quantile(data$time, 0.45))

  for (learner in sb_mboost_aft_learner_ids()) {
    adapter <- sb_mboost_aft_learner_adapter(learner)
    set.seed(9813L)
    caller_rng <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, caller_rng)

    engine_assessment <- assessment
    names(engine_assessment) <- model$engine_features
    native_link <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_assessment,
      type = "link"
    ))
    native_survival <- mboost_aft_native_survival(
      model$params$distribution,
      native_link,
      model$nuisance_scale,
      horizon
    )
    risk <- adapter$predict(model, assessment)
    survival <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      assessment,
      horizon
    )
    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = horizon,
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, assessment)

    expect_s3_class(model$engine, model$params$algorithm)
    expect_identical(as.integer(mboost::mstop(model$engine)), 100L)
    expect_equal(
      model$nuisance_scale,
      as.numeric(mboost::nuisance(model$engine)),
      tolerance = 0
    )
    expect_equal(risk, -native_link, tolerance = 0)
    expect_equal(survival$probability, native_survival, tolerance = 0)
    expect_identical(survival$ties, "not_applicable")
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      survival$probability,
      tolerance = 0
    )
    expect_type(state$engine_raw, "raw")
    expect_identical(
      state$training_data_retention,
      "required_by_mboost_prediction_state"
    )
    expect_false(any(c("time", "event", "outcome") %in% names(state)))

    if (identical(model$params$algorithm, "blackboost")) {
      repeated <- adapter$fit(x, data$time, data$event)
      expect_identical(
        serialize(model$engine, NULL),
        serialize(repeated$engine, NULL)
      )
    }
  }
})

test_that("mboost AFT supports one feature, p greater than n, and horizons", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  cases <- list(
    one_feature = make_mboost_aft_fixture(n = 54L, p = 1L, seed = 9814L),
    high_dimensional = make_mboost_aft_fixture(n = 44L, p = 52L, seed = 9815L)
  )
  for (data in cases) {
    features <- setdiff(names(data), c("patient_id", "time", "event"))
    x <- data[features]
    horizons <- as.numeric(stats::quantile(data$time, c(0.3, 0.7)))
    for (learner in sb_mboost_aft_learner_ids()) {
      adapter <- sb_mboost_aft_learner_adapter(learner)
      model <- adapter$fit(x, data$time, data$event)
      early <- adapter$predict_survival(
        model, x, data$time, data$event,
        x[1:5, , drop = FALSE], horizons[[1L]]
      )$probability
      late <- adapter$predict_survival(
        model, x, data$time, data$event,
        x[1:5, , drop = FALSE], horizons[[2L]]
      )$probability

      expect_true(all(is.finite(early)))
      expect_true(all(is.finite(late)))
      expect_true(all(early >= 0 & early <= 1))
      expect_true(all(late >= 0 & late <= 1))
      expect_true(all(early >= late))
    }
  }
})

test_that("mboost AFT methods execute through CV and final prediction", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  data <- make_mboost_aft_fixture(n = 72L, seed = 9816L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    features,
    "month",
    prediction_horizon = as.numeric(stats::quantile(data$time, 0.45))
  )
  new_data <- data[1:6, c("patient_id", rev(features)), drop = FALSE]

  for (learner in sb_mboost_aft_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 3L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 9817L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      new_data,
      "patient_id"
    )

    expect_identical(pipeline$risk_type, "aft_log_time")
    expect_identical(comparison$status, "success")
    expect_true(comparison$complete_oof)
    expect_true(comparison$complete_survival_probability_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})

# ---- mboost-cox-fixed-adapters ----

test_that("mboost Cox learners keep three distinct native algorithms", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  expect_identical(
    sb_mboost_cox_learner_ids(),
    c(
      "gamboost_cox",
      "glmboost_cox",
      "blackboost_cox"
    )
  )

  d <- make_sb_fixture(n = 60L, seed = 7961L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  expected_class <- c(
    gamboost_cox = "gamboost",
    glmboost_cox = "glmboost",
    blackboost_cox = "blackboost"
  )

  for (learner in sb_mboost_cox_learner_ids()) {
    params <- sb_mboost_cox_params(learner)
    adapter <- sb_mboost_cox_learner_adapter(learner, params)
    model <- adapter$fit(x, d$time, d$event)

    expect_s3_class(model$engine, expected_class[[learner]])
    expect_identical(as.integer(mboost::mstop(model$engine)), 100L)
    expect_identical(model$params$mstop, 100L)
    expect_identical(model$params$nu, 0.1)
    expect_identical(model$params$family, "CoxPH")
  }
})

test_that("mboost Cox adapters return native risk, survival, and replayable state", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  d <- make_sb_fixture(n = 60L, seed = 7962L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(3L, 22L, 51L), , drop = FALSE]

  for (learner in sb_mboost_cox_learner_ids()) {
    adapter <- sb_mboost_cox_learner_adapter(learner)
    model <- adapter$fit(x, d$time, d$event)
    risk <- adapter$predict(model, new_x)
    native_risk <- as.numeric(stats::predict(
      model$engine,
      newdata = new_x,
      type = "link"
    ))
    probability <- adapter$predict_survival(
      model,
      x,
      d$time,
      d$event,
      new_x,
      prediction_horizon = 7
    )
    state <- adapter$extract_state(model, list(
      x = x,
      time = d$time,
      event = d$event,
      prediction_horizon = 7,
      time_unit = "month"
    ))
    replay <- adapter$predict_state(
      unserialize(serialize(state, NULL)),
      new_x
    )

    expect_equal(risk, native_risk, tolerance = 1e-14)
    expect_true(all(is.finite(risk)))
    expect_true(all(probability$probability >= 0 &
      probability$probability <= 1))
    expect_equal(replay$risk, risk, tolerance = 1e-14)
    expect_equal(
      replay$survival_probability,
      probability$probability,
      tolerance = 1e-14
    )
    expect_identical(
      state$training_data_retention,
      "required_by_mboost_prediction_state"
    )
  }
})

test_that("mboost Cox routes complete CV and final outcome-blind prediction", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  d <- make_sb_fixture(n = 60L, seed = 7963L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  new_data <- d[
    c(5L, 27L, 54L),
    c("patient_id", "gene_c", "gene_a", "gene_b")
  ]

  for (learner in sb_mboost_cox_learner_ids()) {
    pipeline <- sb_pipeline(learner = learner)
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 7964L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(fit, new_data, "patient_id")

    expect_identical(result$status, "success")
    expect_true(result$complete_oof)
    expect_true(result$complete_survival_probability_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})

# ---- mboost-ranking-fixed ----

make_mboost_ranking_fixture <- function(n = 84L, p = 6L, seed = 9871L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n)
    signal <- -0.9 * x[, 1L]
    if (p >= 2L) {
      signal <- signal + 0.6 * x[, 2L]
    }
    if (p >= 3L) {
      signal <- signal - 0.4 * x[, 3L]
    }
    event_time <- exp(2.2 + signal + stats::rnorm(n, sd = 0.65))
    censor_time <- exp(stats::rnorm(n, mean = 2.65, sd = 0.8))
    feature_names <- c(
      "gene one",
      if (p >= 2L) "2gene" else NULL,
      if (p >= 3L) paste0("gene-", seq_len(p - 2L)) else NULL
    )
    colnames(x) <- feature_names
    data.frame(
      patient_id = sprintf("MBR%03d", seq_len(n)),
      time = pmin(event_time, censor_time),
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    )
  })
}

mboost_ranking_expected_specifications <- function() {
  data.frame(
    learner = c(
      "gamboost_gehan_aft",
      "gamboost_unoc",
      "glmboost_gehan_aft",
      "glmboost_unoc",
      "blackboost_gehan_aft",
      "blackboost_unoc"
    ),
    algorithm = rep(c("gamboost", "glmboost", "blackboost"), each = 2L),
    family = rep(c("Gehan", "Cindex"), times = 3L),
    objective = rep(c("gehan_rank_aft", "uno_ipcw_concordance"), times = 3L),
    stringsAsFactors = FALSE
  )
}

mboost_ranking_manual_residual_survival <- function(
  analysis_time, analysis_event, analysis_link, assessment_link, horizon
) {
  fit <- survival::survfit(
    survival::Surv(log(analysis_time) - analysis_link, analysis_event) ~ 1
  )
  index <- findInterval(log(horizon) - assessment_link, fit$time)
  probability <- rep(1, length(index))
  use <- index > 0L
  probability[use] <- fit$surv[index[use]]
  unname(probability)
}

test_that("installed mboost ranking losses use longer-survival link direction", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")

  time <- c(1, 2, 3, 4)
  response <- survival::Surv(time, rep(1L, length(time)))
  weights <- rep(1, length(time))
  expect_lt(
    mboost::Gehan()@risk(response, log(time), weights),
    mboost::Gehan()@risk(response, -log(time), weights)
  )
  expect_lt(
    mboost::Cindex(sigma = 0.1, ipcw = 1)@risk(response, time, weights),
    mboost::Cindex(sigma = 0.1, ipcw = 1)@risk(response, -time, weights)
  )
})

test_that("mboost ranking declares six objective-by-base-learner methods", {
  expected <- mboost_ranking_expected_specifications()
  expect_identical(sb_mboost_ranking_learner_ids(), expected$learner)

  for (i in seq_len(nrow(expected))) {
    learner <- expected$learner[[i]]
    params <- sb_mboost_ranking_params(learner)
    expect_identical(params$algorithm, expected$algorithm[[i]])
    expect_identical(params$family, expected$family[[i]])
    expect_identical(params$objective, expected$objective[[i]])
    expect_identical(params$mstop, 100L)
    expect_identical(params$nu, 0.1)
    expect_identical(params$risk, "inbag")
    if (identical(params$family, "Cindex")) {
      expect_identical(params$sigma, 0.1)
      expect_identical(params$ipcw, 1)
      expect_identical(params$calibration_ties, "breslow")
    }
    if (identical(params$algorithm, "blackboost")) {
      expect_identical(params$seed, 1L)
    } else {
      expect_false("seed" %in% names(params))
    }
  }

  registry <- sb_pipeline_registry()
  rows <- registry[
    registry$learner %in% expected$learner, ,
    drop = FALSE
  ]
  expect_identical(nrow(rows), 6L * length(expected_fixed_learner_selectors("survival")))
  expect_identical(anyDuplicated(rows[c("selector", "learner")]), 0L)
  expect_setequal(
    unique(rows$selector),
    expected_fixed_learner_selectors("survival")
  )
  expect_true(all(grepl("mboost", rows$optional_dependency, fixed = TRUE)))
  black_rows <- rows[grepl("^blackboost", rows$learner), , drop = FALSE]
  expect_true(all(grepl("partykit", black_rows$optional_dependency, fixed = TRUE)))
})

test_that("mboost ranking risk and horizon survival match declared methods", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  data <- make_mboost_ranking_fixture(seed = 9872L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[c(4L, 29L, 71L), , drop = FALSE]
  horizon <- as.numeric(stats::quantile(data$time, 0.45))

  for (learner in sb_mboost_ranking_learner_ids()) {
    adapter <- sb_mboost_ranking_learner_adapter(learner)
    set.seed(9873L)
    caller_rng <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, caller_rng)

    engine_analysis <- x
    engine_assessment <- assessment
    names(engine_analysis) <- model$engine_features
    names(engine_assessment) <- model$engine_features
    native_analysis_link <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_analysis,
      type = "link"
    ))
    native_assessment_link <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_assessment,
      type = "link"
    ))
    risk <- adapter$predict(model, assessment)
    survival <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      assessment,
      horizon
    )
    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = horizon,
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, assessment)

    expect_s3_class(model$engine, model$params$algorithm)
    expect_identical(as.integer(mboost::mstop(model$engine)), 100L)
    expect_equal(risk, -native_assessment_link, tolerance = 0)
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      survival$probability,
      tolerance = 0
    )
    expect_type(state$engine_raw, "raw")
    expect_identical(
      state$training_data_retention,
      "required_by_mboost_prediction_state"
    )
    expect_false(any(c("time", "event", "outcome") %in% names(state)))

    if (identical(model$params$family, "Gehan")) {
      expected_probability <- mboost_ranking_manual_residual_survival(
        data$time,
        data$event,
        native_analysis_link,
        native_assessment_link,
        horizon
      )
      expect_equal(survival$probability, expected_probability, tolerance = 0)
      expect_identical(
        survival$source,
        sb_mboost_ranking_direct_survival_source(learner)
      )
      expect_null(state$baseline_hazard)
      expect_true(is.list(state$residual_survival))
    } else {
      calibration <- survival::coxph(
        survival::Surv(data$time, data$event) ~ raw_risk,
        data = data.frame(
          time = data$time,
          event = data$event,
          raw_risk = -native_analysis_link
        ),
        ties = "breslow",
        singular.ok = FALSE
      )
      coefficient <- unname(as.numeric(stats::coef(calibration)))
      calibrated_analysis <- coefficient * -native_analysis_link
      calibrated_assessment <- coefficient * -native_assessment_link
      baseline <- sb_breslow_baseline_hazard(
        data$time,
        data$event,
        calibrated_analysis,
        horizon
      )
      expected_probability <- sb_survival_probability_from_hazard(
        calibrated_assessment,
        baseline
      )
      expect_equal(
        survival$calibration_coefficient,
        coefficient,
        tolerance = 0
      )
      expect_equal(survival$probability, expected_probability, tolerance = 0)
      expect_null(state$residual_survival)
      expect_equal(state$calibration_coefficient, coefficient, tolerance = 0)
      expect_true(is.list(state$baseline_hazard))
    }

    if (identical(model$params$algorithm, "blackboost")) {
      repeated <- adapter$fit(x, data$time, data$event)
      expect_identical(
        serialize(model$engine, NULL),
        serialize(repeated$engine, NULL)
      )
    }
  }
})

test_that("mboost ranking supports one feature and p greater than n", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  cases <- list(
    one_feature = make_mboost_ranking_fixture(n = 54L, p = 1L, seed = 9874L),
    high_dimensional = make_mboost_ranking_fixture(
      n = 44L,
      p = 52L,
      seed = 9875L
    )
  )
  for (data in cases) {
    features <- setdiff(names(data), c("patient_id", "time", "event"))
    x <- data[features]
    horizons <- as.numeric(stats::quantile(data$time, c(0.3, 0.7)))
    for (learner in sb_mboost_ranking_learner_ids()) {
      adapter <- sb_mboost_ranking_learner_adapter(learner)
      model <- adapter$fit(x, data$time, data$event)
      risk <- adapter$predict(model, x[1:5, , drop = FALSE])
      early <- adapter$predict_survival(
        model, x, data$time, data$event,
        x[1:5, , drop = FALSE], horizons[[1L]]
      )$probability
      late <- adapter$predict_survival(
        model, x, data$time, data$event,
        x[1:5, , drop = FALSE], horizons[[2L]]
      )$probability

      expect_true(all(is.finite(risk)))
      expect_true(all(is.finite(early)))
      expect_true(all(is.finite(late)))
      expect_true(all(early >= 0 & early <= 1))
      expect_true(all(late >= 0 & late <= 1))
      expect_true(all(early >= late))
    }
  }
})

test_that("mboost ranking methods execute through CV and final prediction", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  data <- make_mboost_ranking_fixture(n = 72L, seed = 9876L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    features,
    "month",
    prediction_horizon = as.numeric(stats::quantile(data$time, 0.45))
  )
  new_data <- data[1:6, c("patient_id", rev(features)), drop = FALSE]

  for (learner in sb_mboost_ranking_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 3L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 9877L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      new_data,
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_true(comparison$complete_oof)
    expect_true(comparison$complete_survival_probability_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})
})
