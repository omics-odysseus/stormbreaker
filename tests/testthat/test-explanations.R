# ---- additive-shap ----
local({
test_that("logistic additive SHAP exactly reconstructs class-one log odds", {
  development <- make_sb_binary_fixture(n = 80L, seed = 9601L)
  predictors <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    predictors
  )
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 2L
  )
  fit <- sb_fit_final(task, pipeline)

  new_data <- make_sb_binary_fixture(n = 8L, seed = 9602L)
  new_data$patient_id <- paste0("new_", new_data$patient_id)
  background <- development[1:30, c(
    "outcome", rev(fit$required_features)
  )]
  observed <- sb_additive_shap(
    fit,
    new_data,
    "patient_id",
    background
  )

  state <- fit$learner_state
  x <- as.matrix(new_data[fit$required_features])
  background_x <- as.matrix(background[fit$required_features])
  new_terms <- sweep(x, 2L, state$coefficient, "*")
  background_terms <- sweep(
    background_x, 2L, state$coefficient, "*"
  )
  reference_terms <- colMeans(background_terms)
  expected_contributions <- sweep(
    new_terms, 2L, reference_terms, "-"
  )
  expected_reference <- state$intercept + sum(reference_terms)
  expected_score <- expected_reference + rowSums(expected_contributions)
  direct_prediction <- sb_predict_final(
    fit,
    new_data,
    "patient_id"
  )

  expect_s3_class(observed, "stormbreaker_additive_shap")
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$pipeline_id, pipeline$id)
  expect_identical(observed$learner, "logistic")
  expect_identical(
    observed$method,
    "exact_interventional_shap_for_additive_score"
  )
  expect_identical(observed$score_scale, "class_1_log_odds")
  expect_identical(observed$background$n, 30L)
  expect_identical(
    observed$background$feature_summary$feature,
    fit$required_features
  )
  expect_equal(
    observed$background$feature_summary$mean_feature_value,
    unname(colMeans(background_x)),
    tolerance = 0
  )
  expect_equal(
    observed$background$feature_summary$mean_score_term,
    unname(reference_terms),
    tolerance = 0
  )
  expect_equal(observed$reference_value, expected_reference, tolerance = 0)
  expect_equal(
    unname(observed$contributions),
    unname(expected_contributions),
    tolerance = 0
  )
  expect_identical(
    rownames(observed$contributions),
    as.character(new_data$patient_id)
  )
  expect_identical(colnames(observed$contributions), fit$required_features)
  expect_equal(observed$scores$score, expected_score, tolerance = 1e-15)
  expect_equal(
    stats::plogis(observed$scores$score),
    direct_prediction$predictions$probability,
    tolerance = 1e-15
  )
  expect_identical(observed$prediction, direct_prediction)
  expect_false(any(c("outcome", "time", "event") %in% names(observed)))

  changed_new <- new_data
  changed_new$outcome <- 1L - changed_new$outcome
  changed_background <- background
  changed_background$outcome <- 1L - changed_background$outcome
  expect_identical(
    sb_additive_shap(
      fit,
      changed_new,
      "patient_id",
      changed_background
    ),
    observed
  )
})

test_that("Cox additive SHAP stays on the saved risk-score scale", {
  development <- make_sb_fixture(n = 72L, seed = 9701L)
  task <- sb_surv_task(
    development,
    "patient_id",
    "time",
    "event",
    c("gene_b", "gene_a", "gene_c"),
    "month",
    prediction_horizon = 6
  )
  fit <- sb_fit_final(task, sb_pipeline())
  new_data <- make_sb_fixture(n = 7L, seed = 9702L)
  new_data$patient_id <- paste0("new_", new_data$patient_id)
  background <- development[1:24, rev(fit$required_features)]

  observed <- sb_additive_shap(
    fit,
    new_data,
    "patient_id",
    background
  )
  state <- fit$learner_state
  x <- as.matrix(new_data[fit$required_features])
  background_x <- as.matrix(background[fit$required_features])
  expected_reference <- sum(colMeans(background_x) * state$coefficient)
  expected_contributions <- sweep(
    sweep(x, 2L, state$coefficient, "*"),
    2L,
    colMeans(background_x) * state$coefficient,
    "-"
  )

  expect_identical(observed$task_type, "survival")
  expect_identical(
    observed$score_scale,
    "zero_reference_log_hazard_ratio"
  )
  expect_equal(observed$reference_value, expected_reference, tolerance = 0)
  expect_equal(
    unname(observed$contributions),
    unname(expected_contributions),
    tolerance = 1e-15
  )
  expect_equal(
    observed$scores$score,
    observed$prediction$predictions$risk,
    tolerance = 1e-15
  )
  expect_identical(
    "survival_probability" %in% names(observed$prediction$predictions),
    TRUE
  )
  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "not probability-scale effects", fixed = TRUE)
  expect_match(printed, "No feature is ranked", fixed = TRUE)
})

test_that("Gaussian Naive Bayes SHAP reconstructs posterior log odds", {
  development <- make_sb_binary_fixture(n = 84L, seed = 9801L)
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipeline <- sb_pipeline(
    learner = "gaussian_naive_bayes",
    naive_bayes_var_smoothing = 1e-8
  )
  fit <- sb_fit_final(task, pipeline)
  new_data <- make_sb_binary_fixture(n = 9L, seed = 9802L)
  new_data$patient_id <- paste0("new_", new_data$patient_id)
  background <- development[1:28, fit$required_features]

  observed <- sb_additive_shap(
    fit,
    new_data,
    "patient_id",
    background
  )
  model <- structure(
    fit$learner_state$payload,
    class = "stormbreaker_gaussian_naive_bayes"
  )
  term_matrix <- function(data) {
    x <- as.matrix(data[model$features])
    class_term <- function(class) {
      centered <- sweep(x, 2L, model$means[class, ], "-")
      sweep(
        centered^2,
        2L,
        model$variances[class, ],
        "/"
      ) + matrix(
        log(2 * pi * model$variances[class, ]),
        nrow = nrow(x),
        ncol = ncol(x),
        byrow = TRUE
      )
    }
    -0.5 * class_term("1") + 0.5 * class_term("0")
  }
  new_terms <- term_matrix(new_data)
  background_terms <- term_matrix(background)
  reference_terms <- colMeans(background_terms)
  expected_reference <- log(model$priors[["1"]]) -
    log(model$priors[["0"]]) + sum(reference_terms)
  expected_contributions <- sweep(
    new_terms, 2L, reference_terms, "-"
  )

  expect_identical(
    observed$score_scale,
    "class_1_log_posterior_odds"
  )
  expect_equal(observed$reference_value, expected_reference, tolerance = 1e-15)
  expect_equal(
    unname(observed$contributions),
    unname(expected_contributions),
    tolerance = 1e-15
  )
  expect_equal(
    stats::plogis(observed$scores$score),
    observed$prediction$predictions$probability,
    tolerance = 1e-14
  )
})

test_that("glmnet logistic uses the same exact additive explanation", {
  skip_if_not_installed("glmnet")

  development <- make_sb_binary_fixture(n = 72L, seed = 9851L)
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 9852L
  )
  fit <- sb_fit_final(task, pipeline)
  observed <- sb_additive_shap(
    fit,
    development[1:5, ],
    "patient_id",
    development[6:30, ]
  )

  expect_identical(observed$learner, "glmnet_logistic")
  expect_equal(
    stats::plogis(observed$scores$score),
    observed$prediction$predictions$probability,
    tolerance = 1e-14
  )
})

test_that("additive SHAP rejects invalid references and unsupported learners", {
  development <- make_sb_binary_fixture(n = 60L, seed = 9901L)
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  new_data <- development[1:4, ]
  missing_background <- development[5:20, c("gene_a", "gene_b")]
  duplicated_background <- development[5:20, fit$required_features]
  duplicated_background[[ncol(duplicated_background) + 1L]] <-
    duplicated_background[[fit$required_features[[1L]]]]
  names(duplicated_background)[[ncol(duplicated_background)]] <-
    fit$required_features[[1L]]

  expect_s3_class(
    capture_sb_error(sb_additive_shap(
      fit,
      new_data,
      "patient_id",
      missing_background
    )),
    "stormbreaker_missing_final_features"
  )
  expect_s3_class(
    capture_sb_error(sb_additive_shap(
      fit,
      new_data,
      "patient_id",
      duplicated_background
    )),
    "stormbreaker_duplicate_final_feature_columns"
  )
  expect_s3_class(
    capture_sb_error(sb_additive_shap(
      fit,
      new_data,
      "patient_id",
      development[FALSE, fit$required_features]
    )),
    "stormbreaker_bad_adapter_predictors"
  )

  if (requireNamespace("ranger", quietly = TRUE)) {
    ranger_fit <- sb_fit_final(
      task,
      sb_pipeline(
        learner = "ranger_binary",
        ranger_num_trees = 20L,
        ranger_seed = 9902L
      )
    )
    expect_s3_class(
      capture_sb_error(sb_additive_shap(
        ranger_fit,
        new_data,
        "patient_id",
        development[5:20, ]
      )),
      "stormbreaker_unsupported_additive_shap"
    )
  }
})
})

# ---- tree-shap ----
local({
test_that("tree SHAP preserves prediction additivity and seeds", {
  d <- sb_example_binary(90, 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  for (model in c("RF", "GBM", "XGBoost")) {
    fit <- sb_fit_final(task, bi_pipeline(model))
    expect_identical(fit$status, "success")
    s <- sb_tree_shap(fit, d[1:2, ], "patient_id", d, permutations = 3, background_n = 4)
    expect_s3_class(s, "stormbreaker_tree_shap")
    expect_equal(unname(s$reference_value + rowSums(s$contributions)), s$scores$score)
    expect_equal(s, sb_tree_shap(fit, d[1:2, ], "patient_id", d, permutations = 3, background_n = 4))
    expect_identical(s$method, "approximate_permutation_interventional_shap")
    expect_s3_class(sb_plot_shap_importance(s), "ggplot")
    expect_s3_class(sb_plot_shap_beeswarm(s), "ggplot")
    expect_s3_class(sb_plot_shap_waterfall(s), "ggplot")
    expect_s3_class(sb_plot_shap_dependence(s, "gene_a", d[1:2, ], "patient_id"), "ggplot")
  }
})

test_that("one-feature SHAP equals prediction minus empirical reference", {
  d <- sb_example_binary(70, 12)
  t <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  f <- sb_fit_final(t, bi_pipeline("RF"))
  s <- sb_tree_shap(f, d[1:2, ], "patient_id", d[3:7, ], permutations = 1, background_n = 5)
  expected <- sb_predict_final(f, d[3:7, c("patient_id", "gene_a")], "patient_id")$predictions$probability
  expect_equal(s$reference_value, mean(expected))
  expect_equal(as.numeric(s$contributions), s$scores$score - mean(expected))
})

test_that("survival tree SHAP explains saved risk not event probability", {
  d <- sb_example_survival(70, 12)
  t <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month")
  f <- sb_fit_final(t, surv_pipeline("RF"))
  s <- sb_tree_shap(f, d[1:2, ], "patient_id", d, permutations = 2, background_n = 3)
  expect_identical(s$score_scale, "saved_risk_score")
  expect_equal(unname(rowSums(s$contributions) + s$reference_value), s$scores$score)
})
})

# ---- model-shap ----
local({
test_that("LogitBoost and RDA SHAP use original probability predictions", {
  render <- function(p) ggplot2::ggplotGrob(p + ggplot2::theme(text = ggplot2::element_text(family = "sans")))
  d <- sb_example_binary(100, 23)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  for (id in c("logitboost_binary", "rda_binary")) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    s <- sb_model_shap(fit, d[1:3, ], "patient_id", d, permutations = 3, background_n = 4)
    expect_s3_class(s, "stormbreaker_model_shap")
    expect_identical(s$score_scale, "class_1_probability")
    expect_equal(unname(s$reference_value + rowSums(s$contributions)), s$scores$score)
    expect_equal(s$scores$score, sb_predict_final(fit, d[1:3, ], "patient_id")$predictions$probability)
    expect_equal(s, sb_model_shap(fit, d[1:3, ], "patient_id", d, permutations = 3, background_n = 4))
    expect_s3_class(render(sb_plot_shap_beeswarm(s, feature_values = d[1:3, ], id_col = "patient_id")), "gtable")
    expect_s3_class(render(sb_plot_shap_importance(s, groupby = c("A", "B", "A"))), "gtable")
    expect_s3_class(render(sb_plot_shap_waterfall(s)), "gtable")
  }
})
})

# ---- model-coefficients ----
local({
# ---- cox-forest ----

test_that("Cox forest matches native Cox and retains unsuccessful genes", {
  d <- sb_example_survival(100, 42)
  d$constant <- 1
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "constant"), "month")
  result <- sb_cox_forest(task)
  native <- survival::coxph(survival::Surv(time, event) ~ gene_a, data = d)
  row <- result$table[1, ]
  expect_equal(row$hr, unname(exp(stats::coef(native))))
  expect_equal(c(row$lower, row$upper), unname(as.numeric(exp(stats::confint(native)))))
  expect_equal(result$table$p_adjusted, stats::p.adjust(result$table$p_value, "BH"))
  expect_identical(result$table$status[3], "failed")
  expect_identical(result$table$plot_available[3], FALSE)
  expect_s3_class(result$plot, "ggplot")
  expect_equal(nrow(ggplot2::ggplot_build(result$plot)$data[[3]]), 2)
})

# ---- model-coefficients ----

test_that("coefficients retain prediction units and optional reference scaling", {
  d <- make_sb_binary_fixture(n = 60L, seed = 108L)
  fit <- sb_fit_final(sb_binary_task(
    d, "patient_id", "outcome",
    c("gene_a", "gene_b")
  ), sb_pipeline(learner = "logistic"))
  x <- sb_model_coefficients(fit, d)
  expect_equal(x$table$coefficient, fit$learner_state$coefficient)
  expect_equal(x$table$reference_sd, unname(vapply(d[fit$required_features], sd, numeric(1))))
  expect_equal(x$table$scaled_coefficient, x$table$coefficient * x$table$reference_sd)
  expect_equal(x$intercept, fit$learner_state$intercept)
  expect_identical(x$score_scale, "class_1_log_odds")
  expect_identical(x$reference_n, nrow(d))
  expect_equal(sb_model_coefficients(fit)$table$coefficient, x$table$coefficient)
  skip_if_not_installed("ggplot2")
  p <- sb_plot_coefficients(x, scale = "reference_sd")
  expect_s3_class(p, "ggplot")
  expect_equal(sort(p$data$value), sort(x$table$scaled_coefficient))
})
})

# ---- model-curves ----
local({
# ---- model-survival-curves ----

test_that("fixed Efron baseline agrees with native Cox and original predictions", {
  d <- make_sb_fixture(n = 72L, seed = 9701L)
  d$time <- ceiling(d$time)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b"), "month",
    prediction_horizon = 6
  )
  fit <- sb_fit_final(task, sb_pipeline())
  augmented <- sb_add_survival_curve(fit, task)
  expect_identical(augmented$learner_state, fit$learner_state)
  expect_identical(
    sb_predict_final(augmented, d, "patient_id"),
    sb_predict_final(fit, d, "patient_id")
  )
  native <- survival::coxph(survival::Surv(time, event) ~ gene_a + gene_b, d)
  base <- survival::basehaz(native, centered = FALSE)
  at <- augmented$survival_curve$table
  expect_equal(at$cumulative_hazard,
    c(0, base$hazard)[findInterval(at$time, base$time) + 1L],
    tolerance = 1e-7
  )
  curves <- sb_survival_curves(augmented, d[1:3, ], "patient_id")
  expect_equal(dim(curves$survival), c(nrow(at), 3L))
  expect_equal(curves$survival[1, ], rep(1, 3))
  expect_lte(max(apply(curves$survival, 2, diff)), 0)
  expect_equal(
    sb_survival_curves(augmented, d[1, ], "patient_id", times = 6)$survival,
    curves$survival[at$time == 6, 1, drop = FALSE]
  )
  skip_if_not_installed("ggplot2")
  expect_s3_class(sb_plot_survival_curves(curves), "ggplot")
  expect_s3_class(sb_plot_survival_curves(curves, rows = integer()), "ggplot")
  sparse <- sb_survival_curves(augmented, d[1:3, ], "patient_id", times = 6)
  expect_equal(sb_plot_survival_curves(sparse)$data, sb_plot_survival_curves(curves)$data)
})

test_that("Breslow recovery equals an independent risk-set calculation", {
  d <- make_sb_fixture(n = 55L, seed = 9715L)
  d$time <- seq_len(nrow(d))
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 20
  )
  fit <- sb_fit_final(task, sb_pipeline())
  # With unique event times the native Efron and Breslow estimates coincide.
  fit$learner_state$baseline_hazard$ties <- "breslow"
  augmented <- sb_add_survival_curve(fit, task)
  risk <- drop(as.matrix(d[fit$required_features]) %*% fit$learner_state$coefficient)
  t <- d$time[d$event == 1]
  jumps <- vapply(t, function(point) 1 / sum(exp(risk[d$time >= point])), numeric(1))
  expected <- vapply(
    augmented$survival_curve$table$time,
    function(point) sum(jumps[t <= point]), numeric(1)
  )
  expect_equal(augmented$survival_curve$table$cumulative_hazard, expected, tolerance = 1e-12)
})

test_that("wrong development and extrapolation are refused", {
  d <- make_sb_fixture(n = 60L, seed = 105L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 6
  )
  fit <- sb_fit_final(task, sb_pipeline())
  other <- d
  other$patient_id <- paste0("other", other$patient_id)
  wrong <- sb_surv_task(other, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 6
  )
  expect_error(sb_add_survival_curve(fit, wrong), class = "stormbreaker_curve_development_mismatch")
  augmented <- sb_add_survival_curve(fit, task)
  expect_error(sb_survival_curves(augmented, d, "patient_id", max(d$time) + 1),
    class = "stormbreaker_curve_time_unsupported"
  )
  changed <- d
  changed$gene_a <- changed$gene_a + 10
  changed_task <- sb_surv_task(changed, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 6
  )
  expect_error(sb_add_survival_curve(fit, changed_task), class = "stormbreaker_curve_baseline_mismatch")
})

# ---- stable-survival ----

stable_survival_fixture <- function() {
  d <- sb_example_survival(100, seed = 20260819)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month", prediction_horizon = 6)
  fit <- sb_fit_final(task, sb_pipeline())
  external <- sb_example_survival(60, seed = 8921)
  external$patient_id <- paste0("external_", seq_len(nrow(external)))
  external$gene_a <- seq(10, 15, length.out = nrow(external)) /
    fit$learner_state$coefficient[1]
  external$gene_b <- external$gene_c <- 0
  new_data <- external[c("patient_id", "gene_a", "gene_b", "gene_c")]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  endpoint <- sb_surv_endpoint(external, "patient_id", "time", "event",
    "month", prediction_horizon = 6)
  list(fit = fit, external = external, new_data = new_data,
    prediction = prediction, endpoint = endpoint)
}

test_that("saturated Cox predictions preserve their model-defined finite scale", {
  x <- stable_survival_fixture()
  rows <- x$prediction$predictions
  expect_identical(x$prediction$status, "success")
  expect_true(all(rows$survival_probability == 0))
  expect_length(rows$log_cumulative_hazard, nrow(rows))
  expect_equal(rows$log_cumulative_hazard,
    log(x$fit$learner_state$baseline_hazard$cumulative_baseline_hazard) + rows$risk)
  fit <- x$fit
  d <- sb_example_survival(40, seed = 99)
  ordinary <- sb_predict_final(fit, d, "patient_id")
  expect_false("log_cumulative_hazard" %in% names(ordinary$predictions))
  expect_identical(fit, x$fit)
  # Genuine baseline-zero probabilities must not acquire a fabricated scale.
  fit$learner_state$baseline_hazard$cumulative_baseline_hazard <- 0
  expect_null(sb_final_cox_log_cumulative_hazard(fit,
    list(risk = c(1, 2), survival_probability = c(1, 1))))
  fit$learner_state$baseline_hazard$cumulative_baseline_hazard <- 0.5
  expect_equal(sb_final_cox_log_cumulative_hazard(fit,
    list(risk = c(-50, -60), survival_probability = c(1, 1))),
    log(0.5) + c(-50, -60))
  fit$learner_state$baseline_hazard$cumulative_baseline_hazard <- 1e-320
  expect_null(sb_final_cox_log_cumulative_hazard(fit,
    list(risk = c(710, 711), survival_probability = c(0, 0))))
  expect_error(sb_validate_stable_survival_scale(c(0, 1), c(0.3, 0.3)), "disagree")
  expect_error(sb_validate_stable_survival_scale(c(0, Inf), c(0.3, 0.3)), "finite")
})

test_that("external calibration and AUC use the scale, not clipped probabilities", {
  x <- stable_survival_fixture()
  evaluation <- sb_evaluate_external(x$prediction, x$endpoint)
  rows <- evaluation$scored_predictions
  expect_equal(rows$log_cumulative_hazard, x$prediction$predictions$log_cumulative_hazard)
  manual <- survival::coxph(survival::Surv(pmin(time, 6), event == 1 & time <= 6) ~
    log_cumulative_hazard, data = rows, ties = "efron")
  expect_equal(evaluation$metrics$weak_calibration$calibration_slope,
    unname(stats::coef(manual)), tolerance = 1e-10)
  expect_equal(evaluation$metrics$calibration$patient_calibration$transformed_prediction,
    rows$log_cumulative_hazard)
  expect_identical(evaluation$metrics$time_auc$marker, "log_cumulative_hazard")
  expect_identical(evaluation$metrics$calibration$status, "success")
  old <- x$prediction
  old$predictions$log_cumulative_hazard <- NULL
  legacy <- sb_evaluate_external(old, x$endpoint)
  expect_identical(legacy$metrics$weak_calibration$slope_fit$error_class,
    "stormbreaker_survival_weak_calibration_boundary_probability")
  expect_identical(legacy$metrics$calibration$error_class,
    "stormbreaker_survival_calibration_boundary_probability")
  expect_identical(evaluation$metrics$ipcw_brier, legacy$metrics$ipcw_brier)
  expect_identical(evaluation$metrics$harrell_c, legacy$metrics$harrell_c)
  reordered <- x$prediction
  reordered$predictions <- reordered$predictions[rev(seq_len(nrow(rows))), ]
  expect_equal(sb_evaluate_external(reordered, x$endpoint)$scored_predictions,
    evaluation$scored_predictions)
  malformed <- x$prediction
  malformed$predictions$log_cumulative_hazard[1] <- Inf
  expect_identical(sb_evaluate_external(malformed, x$endpoint)$status, "failed")
  contract <- sb_external_bootstrap_contract(evaluation)
  expect_equal(contract$scored$log_cumulative_hazard, rows$log_cumulative_hazard)
  boot <- sb_bootstrap_external(evaluation, reps = 12L, seed = 71L)
  expect_false(identical(boot$status, "failed"))
  expect_true(all(boot$intervals$n_finite_successful[
    boot$intervals$metric %in% c("time_auc", "calibration_slope", "calibration_ici")] > 0L))
})
})

# ---- nomogram ----
local({
# ---- nomogram-aft ----

test_that("AFT nomograms reconstruct native scores and distribution survival", {
  skip_if_not_installed("survival")
  d <- sb_example_survival(160L, seed = 8612L)
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month", prediction_horizon = 7)
  for (id in sb_survreg_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id),
      capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    n <- sb_nomogram(fit, d)
    v <- sb_nomogram_values(n, d)
    engine <- fit$diagnostic_artifacts$learner$engine
    eta <- drop(cbind(1, as.matrix(d[features])) %*% stats::coef(engine))
    expect_equal(v$score, -eta, tolerance = 1e-12)
    expect_equal(n$intercept, -unname(stats::coef(engine)[1]), tolerance = 0)
    expect_equal(n$features$coefficient, -unname(stats::coef(engine)[-1]), tolerance = 0)
    native_distribution <- survival::survreg.distributions[[engine$dist]]
    if (!is.null(native_distribution$dist)) native_distribution <-
      survival::survreg.distributions[[native_distribution$dist]]
    expected <- native_distribution$density((log(7) - eta) / engine$scale)[, 2]
    expect_equal(v$probability, unname(expected), tolerance = 1e-12)
    expect_identical(v$probability_type, "survival_at_saved_horizon")
    expect_null(n$baseline)
    expect_identical(n$score_scale, "negative_log_time_location")
    probabilities <- c(.001, .1, .5, .9, .999)
    expect_equal(sb_nomogram_aft_probability(n$aft, sb_nomogram_aft_score(n$aft, probabilities)),
      probabilities, tolerance = 1e-12)
    p <- sb_plot_nomogram(n)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    expect_identical(tail(attr(p, "stormbreaker_nomogram_axes")$name, 1), "Survival at 7 month")
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_identical(sb_nomogram(compact, d), n)
    no_horizon <- sb_fit_final(sb_surv_task(d, "patient_id", "time", "event", features, "month"),
      sb_pipeline(learner = id))
    plain <- sb_nomogram(no_horizon, d)
    expect_null(sb_nomogram_values(plain, d)$probability)
    expect_equal(sb_nomogram_values(plain, d)$score, -eta, tolerance = 1e-12)
    expect_identical(tail(attr(sb_plot_nomogram(plain), "stormbreaker_nomogram_axes")$name, 1),
      "Negative log-time location")
    constant <- d
    for (feature in features) constant[[feature]] <- d[[feature]][1]
    flat <- sb_nomogram(fit, constant)
    expect_equal(flat$total_points, 0)
    expect_equal(sb_nomogram_values(flat, constant)$probability, rep(v$probability[1], nrow(d)), tolerance = 1e-12)
    expect_s3_class(ggplot2::ggplot_build(sb_plot_nomogram(flat)), "ggplot_built")
    outside <- d[1, , drop = FALSE]
    outside$gene_a <- max(d$gene_a) + 10
    vv <- sb_nomogram_values(n, outside)
    expect_identical(vv$outside_display_range, TRUE)
    expect_equal(vv$score, -unname(drop(cbind(1, as.matrix(outside[features])) %*% stats::coef(engine))), tolerance = 1e-12)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})

# ---- nomogram-readability ----

test_that("nomograms keep instructions outside the plotted figure", {
  d <- sb_example_binary(80, seed = 11)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  nomogram <- sb_nomogram(fit, d)

  plot <- sb_plot_nomogram(nomogram)

  expect_null(plot$labels$caption)
  expect_null(plot$labels$subtitle)
  expect_identical(plot$labels$title, "Nomogram")
  expect_match(paste(attr(plot, "stormbreaker_notes"), collapse = " "),
    "Sum feature points", fixed = TRUE)
  expect_equal(sb_nomogram_values(nomogram, d)$probability,
    sb_predict_final(fit, d, "patient_id")$predictions$probability)
})

# ---- nomogram ----

test_that("nomogram points exactly reconstruct the fitted logistic predictor", {
  d <- make_sb_binary_fixture(n = 60L, seed = 108L)
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b")),
    sb_pipeline(learner = "logistic"))
  original <- fit
  n <- sb_nomogram(fit, d)
  v <- sb_nomogram_values(n, d)
  expected <- drop(as.matrix(d[fit$learner_state$features]) %*% fit$learner_state$coefficient) + fit$learner_state$intercept
  expect_equal(v$score, expected, tolerance = 1e-12)
  expect_equal(v$probability, plogis(expected), tolerance = 1e-12)
  expect_identical(fit, original)
  expect_equal(sum(v$outside_display_range), 0)
  outside <- d[1, ]
  outside$gene_a <- max(d$gene_a) + 100
  vv <- sb_nomogram_values(n, outside)
  expect_identical(vv$outside_display_range, TRUE)
  expect_equal(vv$score, sum(unlist(outside[fit$learner_state$features]) * fit$learner_state$coefficient) + fit$learner_state$intercept)
  skip_if_not_installed("ggplot2")
  expect_s3_class(ggplot2::ggplot_build(sb_plot_nomogram(n)), "ggplot_built")
})

test_that("nomogram display preserves zero terms and spaces short axes", {
  d <- make_sb_binary_fixture(n = 60L, seed = 108L)
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b")),
    sb_pipeline(learner = "logistic"))
  # Synthetic coefficient configurations exercise exact point conversion.
  fit$learner_state$coefficient[] <- c(0, -.00001)
  n <- sb_nomogram(fit, d)
  p <- sb_plot_nomogram(n)
  expect_identical(attr(p, "stormbreaker_nomogram_omitted_zero_terms"), "gene_a")
  expect_equal(nrow(n$features), 2)
  expect_equal(nrow(attr(sb_plot_nomogram(n, show_zero = TRUE), "stormbreaker_nomogram_axes")), 5)
  expect_equal(sb_nomogram_values(n, d)$score,
    fit$learner_state$intercept + d$gene_b * -.00001, tolerance = 1e-12)
  fit$learner_state$coefficient[] <- 0
  z <- sb_nomogram(fit, d)
  expect_equal(z$total_points, 0)
  expect_equal(sb_nomogram_values(z, d)$probability,
    rep(plogis(fit$learner_state$intercept), nrow(d)))
  expect_s3_class(ggplot2::ggplot_build(sb_plot_nomogram(z)), "ggplot_built")
  constant <- d
  constant$gene_a <- 5
  constant$gene_b <- 10
  c <- sb_nomogram(fit, constant)
  expect_equal(sb_nomogram_values(c, constant)$total_points, rep(0, nrow(d)))
})

test_that("penalized nomogram reproduces the original portable predictor", {
  skip_if_not_installed("glmnet")
  d <- make_sb_binary_fixture(n = 90L, seed = 108L)
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b")),
    sb_pipeline(learner = "glmnet_logistic", inner_folds = 3L, inner_seed = 5L))
  expect_identical(fit$status, "success")
  n <- sb_nomogram(fit, d)
  v <- sb_nomogram_values(n, d)
  prediction <- sb_predict_final(fit, d, "patient_id")
  expect_identical(prediction$status, "success")
  expect_equal(v$probability, prediction$predictions$probability, tolerance = 1e-12)
})

test_that("nomogram preserves the original Cox horizon baseline", {
  d <- make_sb_fixture(n = 72L, seed = 9701L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month", prediction_horizon = 6)
  fit <- sb_fit_final(task, sb_pipeline())
  n <- sb_nomogram(fit, d)
  v <- sb_nomogram_values(n, d)
  expected <- drop(as.matrix(d[fit$learner_state$features]) %*% fit$learner_state$coefficient)
  expect_equal(v$score, expected, tolerance = 1e-12)
  expect_equal(v$probability, exp(-fit$learner_state$baseline_hazard$cumulative_baseline_hazard * exp(expected)), tolerance = 1e-12)
  expect_identical(n$baseline, fit$learner_state$baseline_hazard)
  expect_identical(v$probability_type, "survival_at_saved_horizon")
})
})

# ---- plot-shap ----
local({
# ---- plot-shap-single-feature ----

test_that("single-feature waterfalls preserve labels and reconstruct saved scores", {
  skip_if_not_installed("ggplot2")
  for (n in c(1L, 3L)) {
    values <- c(0.4, -0.2, 0)[seq_len(n)]
    contributions <- matrix(values, ncol = 1L,
      dimnames = list(paste0("private-", seq_len(n)), "gene A"))
    x <- structure(list(status = "success", pipeline_id = "single-feature",
      score_scale = "class_1_probability", reference_value = 0.3,
      contributions = contributions,
      scores = data.frame(patient_id = rownames(contributions), score = 0.3 + values),
      background = list(n = 10L)), class = "stormbreaker_model_shap")
    file <- tempfile(fileext = ".rds")
    saveRDS(x, file)
    restored <- readRDS(file)
    for (row in seq_len(n)) for (top_n in c(1L, 10L)) {
      p <- plot_shap_one(restored, row = row, top_n = top_n)
      expect_identical(p$data$feature, c("Background reference", "gene A", "Final model score"))
      expect_identical(p$data$kind, c("reference", "feature", "prediction"))
      expect_equal(p$data$contribution, c(0.3, values[row], x$scores$score[row]))
      expect_equal(tail(p$data$end, 1L), x$scores$score[row])
      expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
    }
    expect_identical(restored, x)
  }
})

# ---- plot-shap ----

local({


test_that("SHAP importance is mean absolute saved contribution on the score scale", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_shap_plot_fixture()
  before <- x
  p <- sb_plot_shap_importance(x)
  expect_equal(p$data$mean_absolute_shap, c(0.75, 0.35, 0.15))
  expect_identical(p$data$feature, c("a", "b", "c"))
  expect_identical(p$labels$x, "Mean |SHAP| (log odds)")
  expect_identical(x, before)
})

test_that("groupby reuses contributions and aligns named groups by ID", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_shap_plot_fixture()
  before <- x
  p <- sb_plot_shap_importance(x, groupby = c("Normal", "Disease"), top_n = 2)
  expect_equal(p$data$mean_absolute_shap, c(1, .3, .5, .4))
  expect_identical(p$data$feature, rep(c("a", "b"), 2))
  expect_equal(p$data$n_observations, rep(1, 4))
  groups <- setNames(c("Disease", "Normal"), rev(rownames(x$contributions)))
  expect_identical(sb_plot_shap_importance(x, top_n = 2, groupby = groups)$data, p$data)
  q <- sb_plot_shap_importance(x, groupby = c("All", "All"))
  expect_equal(q$data$mean_absolute_shap, sb_plot_shap_importance(x)$data$mean_absolute_shap)
  one <- sb_plot_shap_importance(x, top_n = 1, groupby = factor(c("Normal", "Disease"),
    levels = c("Disease", "Normal", "Unused")))
  expect_equal(one$data$mean_absolute_shap, c(.5, 1))
  expect_identical(levels(one$data$group), c("Disease", "Normal"))
  built <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(built$data[[1]]$y)), 4)
  expect_equal(built$data[[1]]$y, built$data[[2]]$y)
  expect_null(p$labels$subtitle)
  expect_identical(x, before)
  expect_error(sb_plot_shap_importance(x, groupby = "Normal"), class = "stormbreaker_bad_visualization_argument")
  expect_error(sb_plot_shap_importance(x, groupby = c(NA, "Normal")), class = "stormbreaker_bad_visualization_argument")
  expect_error(sb_plot_shap_importance(x, groupby = c("", "Normal")), class = "stormbreaker_bad_visualization_argument")
  expect_error(sb_plot_shap_importance(x, groupby = c(wrong = "Normal", other = "Disease")),
    class = "stormbreaker_bad_visualization_argument")
})

test_that("waterfall retains baseline and remainder and reconstructs the full score", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_shap_plot_fixture()
  p <- sb_plot_shap_waterfall(x, row = 1, top_n = 1)
  expect_identical(p$data$kind, c("reference", "feature", "remainder", "prediction"))
  expect_equal(p$data$contribution, c(-0.4, 1, -0.1, 0.5))
  expect_equal(p$data$end, c(-0.4, 0.6, 0.5, 0.5))
  expect_equal(tail(p$data$end, 1), x$scores$score[1])
  expect_match(p$labels$subtitle, "Observation 1", fixed = TRUE)
  expect_identical(any(grepl("patient-secret", unlist(p$data))), FALSE)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that("beeswarm is deterministic and feature values align by explicit IDs", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_shap_plot_fixture()
  set.seed(29)
  before <- .Random.seed
  p <- sb_plot_shap_beeswarm(x)
  q <- sb_plot_shap_beeswarm(x)
  expect_identical(p$data, q$data)
  expect_identical(.Random.seed, before)
  expect_equal(nrow(p$data), 6L)
  feature_values <- data.frame(id = rev(rownames(x$contributions)), a = c(20, 10), b = c(40, 30), c = c(60, 50))
  coloured <- sb_plot_shap_beeswarm(x, feature_values = feature_values, id_col = "id")
  expect_equal(coloured$data$feature_value, c(10, 20, 30, 40, 50, 60))
  expect_identical(any(grepl("patient-secret", unlist(coloured$data))), FALSE)
  expect_error(sb_plot_shap_beeswarm(x, feature_values = feature_values),
    class = "stormbreaker_bad_visualization_argument")
})

test_that("bad attribution states are not silently relabelled as valid explanations", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_shap_plot_fixture()
  expect_error(sb_plot_shap_waterfall(x, row = 3), class = "stormbreaker_bad_visualization_argument")
  x$scores$score[1] <- 100
  expect_error(sb_plot_shap_waterfall(x), class = "stormbreaker_bad_shap_visualization")
})

})

# ---- shap-dependence ----

test_that("dependence aligns IDs and preserves all contributions", {
  skip_if_not_installed("ggplot2")
  m <- matrix(c(1, 2, -1, -2), 2, dimnames = list(c("secret1", "secret2"), c("a", "b")))
  x <- structure(list(status = "success", pipeline_id = "logistic", score_scale = "log_odds",
    reference_value = 0, contributions = m, scores = data.frame(patient_id = rownames(m), score = rowSums(m)),
    background = list(n = 10L)), class = "stormbreaker_additive_shap")
  values <- data.frame(id = c("secret2", "secret1"), a = c(20, 10), b = c(4, 3))
  before <- x
  p <- sb_plot_shap_dependence(x, "a", values, "id", colour_by = "b")
  expect_equal(p$data$feature_value, c(10, 20))
  expect_equal(p$data$contribution, c(1, 2))
  expect_equal(p$data$colour_value, c(3, 4))
  expect_identical(x, before)
  expect_identical(any(grepl("secret", unlist(p$data))), FALSE)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

# ---- shap-visual-style ----

local({


test_that("importance carries numeric labels without altering the ranked values", {
  x <- make_sb_shap_plot_fixture()
  p <- sb_plot_shap_importance(x)
  expect_identical(p$data$value_label, c("0.75", "0.35", "0.15"))
  expect_equal(p$data$mean_absolute_shap, unname(colMeans(abs(x$contributions))))
})

test_that("waterfall has a readable separate value column and unchanged endpoints", {
  x <- make_sb_shap_plot_fixture()
  p <- sb_plot_shap_waterfall(x, top_n = 1)
  expect_identical(p$data$value_label, c("-0.4", "+1", "-0.1", "0.5"))
  expect_equal(length(unique(p$data$label_x)), 1L)
  expect_gt(p$data$label_x[1], max(p$data$end))
  expect_equal(tail(p$data$end, 1), x$scores$score[1])
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that("zero and tiny contributions remain finite and reconstruct exactly", {
  x <- make_sb_shap_plot_fixture()
  x$contributions[,] <- 0
  x$reference_value <- 0
  x$scores$score[] <- 0
  p <- sb_plot_shap_waterfall(x)
  expect_identical(all(is.finite(p$data$label_x)), TRUE)
  expect_identical(p$data$value_label, c("0", "+0", "+0", "+0", "0"))
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  x$contributions[1,1:2] <- c(1e-9, -1e-9)
  q <- sb_plot_shap_waterfall(x)
  expect_identical(q$data$value_label[2:3], c("+1e-09", "-1e-09"))
  expect_equal(tail(q$data$end, 1), 0)
})
})
})

test_that("SHAP heatmap keeps every observation and shows signed saved contributions", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_shap_plot_fixture()
  p <- sb_plot_shap_heatmap(x, top_n = 2L)
  expect_identical(p$data$observation, c(1L, 2L, 1L, 2L))
  expect_identical(as.character(p$data$feature), c("a", "a", "b", "b"))
  expect_equal(p$data$contribution, c(1, -.5, -.3, .4))
  expect_false(any(grepl("patient-secret", unlist(p$data))))
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  one <- x
  one$contributions <- one$contributions[1L, "a", drop = FALSE] * 0
  one$scores <- data.frame(patient_id = rownames(one$contributions), score = one$reference_value)
  expect_equal(sb_plot_shap_heatmap(one)$data$contribution, 0)
})
