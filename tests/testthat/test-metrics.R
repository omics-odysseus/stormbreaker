# ---- metrics-binary ----
local({
# ---- metrics-binary-decision-curve ----

test_that("binary decision curve matches hand-counted net benefit", {
  scored <- data.frame(
    outcome = c(1L, 1L, 0L, 0L),
    probability = c(0.90, 0.50, 0.70, 0.20),
    status = rep("success", 4L)
  )

  metric <- stormbreaker:::sb_binary_decision_curve_metric(
    scored,
    c(0.30, 0.50, 0.80),
    "hand_calculation"
  )

  expect_identical(metric$status, "success")
  expect_identical(
    metric$decision_rule,
    "probability_greater_than_or_equal_to_threshold"
  )
  expect_identical(metric$curve$true_positive, c(2L, 2L, 1L))
  expect_identical(metric$curve$false_positive, c(1L, 1L, 0L))
  expect_identical(metric$curve$true_negative, c(1L, 1L, 2L))
  expect_identical(metric$curve$false_negative, c(0L, 0L, 1L))
  expect_identical(metric$curve$n_selected, c(3L, 3L, 1L))
  expect_equal(metric$curve$sensitivity, c(1, 1, 0.5), tolerance = 0)
  expect_equal(metric$curve$specificity, c(0.5, 0.5, 1), tolerance = 0)
  expect_equal(
    metric$curve$positive_predictive_value,
    c(2 / 3, 2 / 3, 1),
    tolerance = 0
  )
  expect_equal(
    metric$curve$negative_predictive_value,
    c(1, 1, 2 / 3),
    tolerance = 0
  )
  expect_equal(
    metric$curve$model_net_benefit,
    c(0.5 - 0.25 * 0.3 / 0.7, 0.5 - 0.25, 0.25),
    tolerance = 1e-15
  )
  expect_equal(
    metric$curve$treat_all_net_benefit,
    0.5 - 0.5 * c(0.3 / 0.7, 1, 4),
    tolerance = 1e-15
  )
  expect_identical(
    metric$curve$treat_none_net_benefit,
    c(0, 0, 0)
  )
})

test_that("binary decision curve is invariant to patient order", {
  scored <- data.frame(
    outcome = c(1L, 0L, 1L, 0L),
    probability = c(0.8, 0.6, 0.4, 0.1),
    status = rep("success", 4L)
  )

  direct <- stormbreaker:::sb_binary_decision_curve_metric(
    scored,
    c(0.2, 0.5),
    "test"
  )
  reversed <- stormbreaker:::sb_binary_decision_curve_metric(
    scored[4:1, ],
    c(0.2, 0.5),
    "test"
  )

  expect_identical(reversed$curve, direct$curve)
})

test_that("binary decision curve never promotes a partial patient subset", {
  scored <- data.frame(
    outcome = c(1L, 0L, 1L),
    probability = c(0.8, NA_real_, 0.2),
    status = c("success", "failed", "success")
  )

  metric <- stormbreaker:::sb_binary_decision_curve_metric(
    scored,
    c(0.2, 0.5),
    "test"
  )

  expect_identical(metric$status, "incomplete")
  expect_identical(metric$n, 3L)
  expect_identical(metric$n_scored, 0L)
  expect_true(all(is.na(metric$curve$model_net_benefit)))
  expect_match(metric$error_message, "incomplete or invalid")
})

test_that("binary decision curve handles boundaries and one class", {
  scored <- data.frame(
    outcome = c(0L, 0L, 0L),
    probability = c(0, 0.5, 1),
    status = rep("success", 3L)
  )

  metric <- stormbreaker:::sb_binary_decision_curve_metric(
    scored,
    c(0.25, 0.75),
    "test"
  )

  expect_identical(metric$status, "success")
  expect_identical(metric$n_positive, 0L)
  expect_identical(metric$n_negative, 3L)
  expect_identical(metric$curve$true_positive, c(0L, 0L))
  expect_identical(metric$curve$false_positive, c(2L, 1L))
  expect_true(all(is.na(metric$curve$sensitivity)))
  expect_equal(metric$curve$specificity, c(1 / 3, 2 / 3), tolerance = 0)
  expect_identical(metric$curve$positive_predictive_value, c(0, 0))
  expect_identical(metric$curve$negative_predictive_value, c(1, 1))
  expect_equal(
    metric$curve$treat_all_net_benefit,
    -c(0.25 / 0.75, 0.75 / 0.25),
    tolerance = 1e-15
  )

  scored$outcome <- rep.int(1L, nrow(scored))
  positive_only <- stormbreaker:::sb_binary_decision_curve_metric(
    scored,
    c(0.25, 0.75),
    "test"
  )
  expect_equal(
    positive_only$curve$sensitivity,
    c(2 / 3, 1 / 3),
    tolerance = 0
  )
  expect_true(all(is.na(positive_only$curve$specificity)))
  expect_identical(positive_only$curve$positive_predictive_value, c(1, 1))
  expect_identical(positive_only$curve$negative_predictive_value, c(0, 0))
})

test_that("undefined predictive-value denominators stay explicit", {
  scored <- data.frame(
    outcome = c(1L, 0L),
    probability = c(0.9, 0.8),
    status = rep("success", 2L)
  )

  metric <- stormbreaker:::sb_binary_decision_curve_metric(
    scored,
    c(0.5, 0.95),
    "test"
  )

  expect_identical(metric$curve$sensitivity, c(1, 0))
  expect_identical(metric$curve$specificity, c(0, 1))
  expect_equal(
    metric$curve$positive_predictive_value,
    c(0.5, NA_real_),
    tolerance = 0
  )
  expect_equal(
    metric$curve$negative_predictive_value,
    c(NA_real_, 0.5),
    tolerance = 0
  )
})

test_that("binary decision curve exposes reconstructible curve columns", {
  scored <- data.frame(
    outcome = c(1L, 0L),
    probability = c(0.9, 0.1),
    status = rep("success", 2L)
  )

  metric <- stormbreaker:::sb_binary_decision_curve_metric(
    scored,
    0.5,
    "test"
  )

  expect_named(
    metric$curve,
    c(
      "threshold_probability",
      "threshold_odds",
      "n",
      "n_positive",
      "n_negative",
      "n_selected",
      "true_positive",
      "false_positive",
      "true_negative",
      "false_negative",
      "sensitivity",
      "specificity",
      "positive_predictive_value",
      "negative_predictive_value",
      "model_net_benefit",
      "treat_all_net_benefit",
      "treat_none_net_benefit"
    )
  )
  expect_identical(metric$n_thresholds, 1L)
  expect_identical(metric$threshold_source, "study_declared_binary_task")
})

# ---- metrics-binary-moderate-calibration ----

make_binary_moderate_scored <- function(n = 24L) {
  probability <- seq(0.03, 0.97, length.out = n)
  data.frame(
    patient_id = paste0("patient_", seq_len(n)),
    fold = rep(seq_len(4L), length.out = n),
    outcome = rep(c(0L, 1L), length.out = n),
    probability = probability,
    status = rep("success", n),
    stringsAsFactors = FALSE
  )
}

manual_binary_rcs3_basis <- function(x, knots) {
  d <- function(value, knot) {
    (pmax(value - knot, 0)^3 - pmax(value - knots[[3L]], 0)^3) /
      (knots[[3L]] - knot)
  }
  cbind(
    rcs_linear = x,
    rcs_nonlinear = d(x, knots[[1L]]) - d(x, knots[[2L]])
  )
}

test_that("binary moderate calibration agrees with an independent RCS oracle", {
  scored <- make_binary_moderate_scored()
  result <- sb_binary_moderate_calibration_metric(
    scored,
    "pooled_patient_level_out_of_fold"
  )

  logit_probability <- stats::qlogis(scored$probability)
  knots <- unname(stats::quantile(
    logit_probability,
    c(0.10, 0.50, 0.90),
    names = FALSE,
    type = 7
  ))
  manual_basis <- manual_binary_rcs3_basis(logit_probability, knots)
  manual_design <- cbind("(Intercept)" = 1, manual_basis)
  manual_fit <- stats::glm.fit(
    x = manual_design,
    y = scored$outcome,
    family = stats::binomial()
  )
  manual_probability <- stats::plogis(drop(
    manual_design %*% manual_fit$coefficients
  ))
  manual_error <- abs(manual_probability - scored$probability)

  expect_identical(result$status, "success")
  expect_identical(result$input_status, "success")
  expect_identical(result$metric, "binary_calibration_curve")
  expect_identical(result$calibration_level, "moderate")
  expect_identical(
    result$method,
    "logistic_three_knot_restricted_natural_cubic_spline"
  )
  expect_identical(result$prediction_transform, "logit")
  expect_identical(result$probability_clipping, "none")
  expect_identical(result$spline$n_knots, 3L)
  expect_equal(result$knots$logit_probability, knots, tolerance = 1e-14)
  expect_equal(result$spline$boundary_knots, knots[c(1L, 3L)])
  expect_equal(result$spline$interior_knots, knots[[2L]])
  expect_identical(
    result$patient_calibration$patient_id,
    scored$patient_id
  )
  expect_identical(
    result$patient_calibration$row_instance,
    seq_len(nrow(scored))
  )
  expect_equal(
    result$patient_calibration$calibrated_probability,
    manual_probability,
    tolerance = 1e-10
  )
  expect_equal(result$ici, mean(manual_error), tolerance = 1e-10)
  expect_equal(
    result$e50,
    unname(stats::quantile(manual_error, 0.50, names = FALSE, type = 7)),
    tolerance = 1e-10
  )
  expect_equal(
    result$e90,
    unname(stats::quantile(manual_error, 0.90, names = FALSE, type = 7)),
    tolerance = 1e-10
  )
  expect_identical(result$fit$terms, result$fit$coefficients$term)
  expect_identical(result$fit$rank, 3L)
  expect_identical(result$fit$converged, TRUE)
  expect_identical(result$fit$boundary, FALSE)
  expect_identical(result$fit$warnings, character())
  expect_identical(result$fit$separation_status, "not_detected")
  expect_identical(result$fit$separation_detected, FALSE)
  expect_identical(result$per_fold$status, rep("success", 4L))
})

test_that("saved audit alone reconstructs the binary calibration evidence", {
  scored <- make_binary_moderate_scored()
  result <- sb_binary_moderate_calibration_metric(
    scored,
    "pooled_patient_level_out_of_fold"
  )
  saved_knots <- result$knots$logit_probability
  saved_basis <- splines::ns(
    result$patient_calibration$logit_probability,
    knots = result$spline$interior_knots,
    Boundary.knots = result$spline$boundary_knots,
    intercept = result$spline$intercept
  )
  colnames(saved_basis) <- c("calibration_ns1", "calibration_ns2")
  saved_design <- cbind("(Intercept)" = 1, saved_basis)
  saved_coefficients <- stats::setNames(
    result$fit$coefficients$estimate,
    result$fit$coefficients$term
  )[result$fit$terms]
  reconstructed <- stats::plogis(drop(saved_design %*% saved_coefficients))
  reconstructed_error <- abs(reconstructed - scored$probability)
  expected_curve <- data.frame(
    probability = sort(unique(scored$probability)),
    stringsAsFactors = FALSE
  )
  first <- match(expected_curve$probability, scored$probability)
  expected_curve$logit_probability <- stats::qlogis(
    expected_curve$probability
  )
  expected_curve$calibrated_probability <- reconstructed[first]
  expected_curve$n_row_instances <- as.integer(tabulate(
    match(scored$probability, expected_curve$probability),
    nbins = nrow(expected_curve)
  ))

  expect_equal(saved_knots, c(
    result$spline$boundary_knots[[1L]],
    result$spline$interior_knots[[1L]],
    result$spline$boundary_knots[[2L]]
  ))
  expect_equal(
    reconstructed,
    result$patient_calibration$calibrated_probability,
    tolerance = 1e-14
  )
  expect_equal(result$curve, expected_curve, tolerance = 1e-14)
  expect_equal(result$ici, mean(reconstructed_error), tolerance = 1e-14)
  expect_equal(
    result$e50,
    unname(stats::quantile(
      reconstructed_error,
      0.50,
      names = FALSE,
      type = 7
    )),
    tolerance = 1e-14
  )
  expect_equal(
    result$e90,
    unname(stats::quantile(
      reconstructed_error,
      0.90,
      names = FALSE,
      type = 7
    )),
    tolerance = 1e-14
  )
})

test_that("binary moderate calibration preserves ties rows and pooled meaning", {
  scored <- make_binary_moderate_scored()
  scored$probability[c(3L, 4L)] <- scored$probability[[2L]]
  scored$patient_id[c(7L, 8L)] <- "resampled_patient"
  result <- sb_binary_moderate_calibration_metric(
    scored,
    "patient_level_external_bootstrap_sample"
  )

  expect_identical(result$status, "success")
  expect_identical(nrow(result$patient_calibration), nrow(scored))
  expect_identical(
    sum(result$patient_calibration$patient_id == "resampled_patient"),
    2L
  )
  expect_identical(
    sum(result$curve$n_row_instances),
    as.integer(nrow(scored))
  )
  expect_identical(max(result$curve$n_row_instances), 3L)

  one_class_fold <- make_binary_moderate_scored()
  one_class_fold$outcome[one_class_fold$fold == 1L] <- 0L
  pooled <- sb_binary_moderate_calibration_metric(
    one_class_fold,
    "pooled_patient_level_out_of_fold"
  )
  expect_identical(pooled$status, "success")
  expect_identical(pooled$per_fold$status, rep("success", 4L))
  expect_identical(pooled$per_fold$n_positive[[1L]], 0L)
})

test_that("binary moderate calibration is invariant to patient row order", {
  scored <- make_binary_moderate_scored()
  original <- sb_binary_moderate_calibration_metric(
    scored,
    "patient_level_external_cohort"
  )
  order <- c(seq(2L, nrow(scored), by = 2L), seq(1L, nrow(scored), by = 2L))
  shuffled_scored <- scored[order, , drop = FALSE]
  shuffled <- sb_binary_moderate_calibration_metric(
    shuffled_scored,
    "patient_level_external_cohort"
  )

  expect_identical(shuffled$status, "success")
  expect_equal(shuffled[c("ici", "e50", "e90")], original[c(
    "ici", "e50", "e90"
  )], tolerance = 1e-13)
  expect_equal(shuffled$knots, original$knots, tolerance = 1e-13)
  expect_equal(shuffled$curve, original$curve, tolerance = 1e-13)
  aligned <- shuffled$patient_calibration[
    match(original$patient_calibration$patient_id,
      shuffled$patient_calibration$patient_id
    ),
    ,
    drop = FALSE
  ]
  expect_equal(
    aligned$calibrated_probability,
    original$patient_calibration$calibrated_probability,
    tolerance = 1e-13
  )
})

test_that("binary moderate calibration retains concrete unavailable states", {
  base <- make_binary_moderate_scored()

  failed <- base
  failed$status[[1L]] <- "prediction_failed"
  failed_result <- sb_binary_moderate_calibration_metric(failed, "test")
  expect_identical(failed_result$status, "incomplete")
  expect_identical(failed_result$input_status, "prediction_failed")

  one_class <- base
  one_class$outcome <- 0L
  class_result <- sb_binary_moderate_calibration_metric(one_class, "test")
  expect_identical(class_result$status, "incomplete")
  expect_identical(class_result$input_status, "not_evaluable")
  expect_match(class_result$error_message, "Both outcome classes", fixed = TRUE)

  boundary <- base
  boundary$probability[[1L]] <- 0
  boundary_result <- sb_binary_moderate_calibration_metric(boundary, "test")
  expect_identical(boundary_result$status, "incomplete")
  expect_identical(boundary_result$input_status, "not_evaluable")
  expect_identical(boundary_result$n_boundary_probability, 1L)
  expect_identical(boundary_result$probability_clipping, "none")

  tied <- base
  tied$probability <- 0.5
  tied_result <- sb_binary_moderate_calibration_metric(tied, "test")
  expect_identical(tied_result$status, "incomplete")
  expect_identical(tied_result$input_status, "not_evaluable")
  expect_match(tied_result$error_message, "strictly increasing", fixed = TRUE)
})

test_that("extreme finite probabilities are neither clipped nor rejected", {
  probability <- c(
    1e-300, 1e-20, 1e-8, 1e-4, 0.01, 0.05, 0.15, 0.30,
    0.50, 0.70, 0.85, 0.95, 0.99, 1 - 1e-4, 1 - 1e-8,
    1 - .Machine$double.eps
  )
  scored <- data.frame(
    patient_id = paste0("extreme_", seq_along(probability)),
    outcome = rep(c(0L, 1L), length.out = length(probability)),
    probability = probability,
    status = rep("success", length(probability)),
    stringsAsFactors = FALSE
  )
  result <- sb_binary_moderate_calibration_metric(scored, "test")

  expect_identical(result$status, "success")
  expect_identical(result$n_boundary_probability, 0L)
  expect_identical(
    result$patient_calibration$probability,
    probability
  )
  expect_equal(
    result$patient_calibration$logit_probability,
    stats::qlogis(probability),
    tolerance = 0
  )
  expect_identical(
    all(is.finite(result$patient_calibration$logit_probability)),
    TRUE
  )
  expect_identical(result$probability_clipping, "none")
  expect_identical(all(is.finite(c(result$ici, result$e50, result$e90))), TRUE)
})

test_that("probability warnings retain an estimable calibration fit", {
  scored <- make_binary_moderate_scored()
  expected <- sb_binary_moderate_calibration_metric(scored, "test")
  native_fit <- stats::glm.fit
  probability_warning <-
    "glm.fit: fitted probabilities numerically 0 or 1 occurred"
  local_mocked_bindings(
    glm.fit = function(...) {
      fit <- native_fit(...)
      warning(probability_warning, call. = FALSE)
      fit
    },
    .package = "stats"
  )

  result <- sb_binary_moderate_calibration_metric(scored, "test")

  expect_identical(result$status, "success")
  expect_identical(result$fit$status, "success")
  expect_identical(result$fit$warnings, probability_warning)
  expect_identical(result$fit$converged, TRUE)
  expect_identical(result$fit$rank, 3L)
  expect_identical(result$fit$separation_status, "not_detected")
  expect_identical(result$fit$separation_detected, FALSE)
  expect_equal(result$fit$coefficients, expected$fit$coefficients, tolerance = 0)
  expect_equal(result$patient_calibration, expected$patient_calibration,
    tolerance = 0
  )
  expect_equal(result[c("ici", "e50", "e90")],
    expected[c("ici", "e50", "e90")], tolerance = 0
  )
})

test_that("calibration warnings do not rescue invalid native fits", {
  scored <- make_binary_moderate_scored()
  native_fit <- stats::glm.fit
  probability_warning <-
    "glm.fit: fitted probabilities numerically 0 or 1 occurred"
  cases <- list(
    list(warning = "unrelated native warning", changes = list(),
      error = "glm_warning"),
    list(warning = "glm.fit: algorithm did not converge", changes = list(),
      error = "glm_warning"),
    list(warning = probability_warning, changes = list(converged = FALSE),
      error = "nonconvergence"),
    list(warning = probability_warning, changes = list(boundary = TRUE),
      error = "boundary_fit"),
    list(warning = probability_warning, changes = list(rank = 2L),
      error = "singular_fit")
  )
  for (case in cases) {
    result <- with_mocked_bindings(
      sb_binary_moderate_calibration_metric(scored, "test"),
      glm.fit = function(...) {
        fit <- native_fit(...)
        fit[names(case$changes)] <- case$changes
        warning(case$warning, call. = FALSE)
        fit
      },
      .package = "stats"
    )
    expect_identical(result$status, "fit_failed")
    expect_identical(result$fit$status, "failed")
    expect_identical(result$fit$warnings, case$warning)
    expect_identical(result$error_class, paste0(
      "stormbreaker_binary_moderate_calibration_", case$error
    ))
    expect_identical(is.na(result$ici), TRUE)
  }
})

test_that("a separated spline fit remains visible without fallback", {
  probability <- seq(0.03, 0.97, length.out = 40L)
  scored <- data.frame(
    patient_id = paste0("separated_", seq_along(probability)),
    outcome = as.integer(probability >= 0.5),
    probability = probability,
    status = rep("success", length(probability)),
    stringsAsFactors = FALSE
  )
  result <- sb_binary_moderate_calibration_metric(scored, "test")

  expect_identical(result$status, "fit_failed")
  expect_identical(result$input_status, "success")
  expect_identical(result$fit$status, "failed")
  expect_gt(length(result$fit$warnings), 0L)
  expect_true(is.na(result$ici))
  expect_identical(result$method, "logistic_three_knot_restricted_natural_cubic_spline")
})

test_that("silent finite IRLS output cannot hide nonlinear separation", {
  probability <- c(
    0.013529309201985598,
    0.232556497720070199,
    0.292985442220233372,
    0.376448372132144882,
    0.875341706429608180,
    0.887932648812420622
  )
  scored <- data.frame(
    patient_id = paste0("silent_separation_", seq_along(probability)),
    outcome = c(1L, 0L, 0L, 0L, 1L, 1L),
    probability = probability,
    status = rep("success", length(probability)),
    stringsAsFactors = FALSE
  )

  result <- sb_binary_moderate_calibration_metric(scored, "test")

  expect_identical(result$status, "fit_failed")
  expect_identical(result$input_status, "success")
  expect_identical(result$fit$status, "failed")
  expect_identical(result$fit$converged, TRUE)
  expect_identical(result$fit$boundary, FALSE)
  expect_identical(result$fit$warnings, character())
  expect_identical(result$fit$separation_status, "detected")
  expect_identical(result$fit$separation_detected, TRUE)
  expect_identical(result$fit$separation_solver_status, 0L)
  expect_gt(
    result$fit$separation_objective,
    result$fit$separation_tolerance
  )
  expect_gte(
    result$fit$separation_minimum_signed_margin,
    -result$fit$separation_tolerance
  )
  expect_gt(
    result$fit$separation_maximum_signed_margin,
    result$fit$separation_tolerance
  )
  expect_identical(
    result$error_class,
    "stormbreaker_binary_moderate_calibration_separation"
  )
  expect_identical(all(is.finite(result$fit$coefficients$estimate)), TRUE)
  expect_identical(is.na(result$ici), TRUE)
})

test_that("the signed-margin program detects quasi-complete separation", {
  probability <- c(0.05, 0.20, 0.50, 0.50, 0.80, 0.95)
  outcome <- c(0L, 0L, 0L, 1L, 1L, 1L)
  logit_probability <- stats::qlogis(probability)
  knots <- sb_binary_moderate_calibration_knots(
    logit_probability
  )$logit_probability
  basis <- sb_binary_moderate_calibration_basis(logit_probability, knots)
  design <- cbind("(Intercept)" = 1, basis)

  result <- sb_binary_moderate_calibration_separation(design, outcome)

  expect_identical(result$status, "detected")
  expect_identical(result$detected, TRUE)
  expect_identical(result$solver_status, 0L)
  expect_gt(result$objective, result$tolerance)
  expect_lte(abs(result$minimum_signed_margin), result$tolerance)
  expect_gt(result$maximum_signed_margin, result$tolerance)
})

# ---- metrics-binary ----

test_that("fold-stratified ROC AUC aggregates only within-model pairs", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:8),
    fold = rep(1:2, each = 4),
    outcome = c(1L, 1L, 0L, 0L, 1L, 1L, 0L, 0L),
    probability = c(0.9, 0.5, 0.4, 0.1, 0.6, 0.2, 0.6, 0.3),
    status = "success",
    stringsAsFactors = FALSE
  )

  metric <- sb_binary_auc_metric(predictions)

  expect_identical(metric$status, "success")
  expect_identical(metric$scope, "within_outer_assessment_fold")
  expect_identical(metric$comparable_pairs, 8)
  expect_identical(metric$correctly_ordered_pairs, 5)
  expect_identical(metric$tied_pairs, 1)
  expect_equal(metric$estimate, 5.5 / 8, tolerance = 0)
  expect_equal(metric$per_fold$estimate, c(1, 0.375), tolerance = 0)
})

test_that("failed binary predictions cannot produce a headline AUC", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:4),
    fold = rep(1:2, each = 2),
    outcome = c(1L, 0L, 1L, 0L),
    probability = c(0.8, 0.2, NA_real_, NA_real_),
    status = c("success", "success", "failed", "failed"),
    stringsAsFactors = FALSE
  )

  metric <- sb_binary_auc_metric(predictions)

  expect_identical(metric$status, "incomplete")
  expect_identical(is.na(metric$estimate), TRUE)
  expect_identical(
    metric$per_fold$status,
    c("success", "prediction_failed")
  )
})

test_that("binary AUC stratum exposes exact pair accounting", {
  audit <- sb_binary_auc_stratum(
    outcome = c(1L, 1L, 0L, 0L),
    probability = c(0.6, 0.2, 0.6, 0.3)
  )

  expect_identical(audit$status, "success")
  expect_identical(audit$n, 4L)
  expect_identical(audit$n_positive, 2L)
  expect_identical(audit$n_negative, 2L)
  expect_equal(audit$correctly_ordered_pairs, 1, tolerance = 0)
  expect_equal(audit$tied_pairs, 1, tolerance = 0)
  expect_equal(audit$comparable_pairs, 4, tolerance = 0)
  expect_equal(audit$estimate, 1.5 / 4, tolerance = 0)
})

test_that("binary AUC stratum is explicitly not evaluable with one class", {
  audit <- sb_binary_auc_stratum(
    outcome = c(0L, 0L),
    probability = c(0.1, 0.2)
  )

  expect_identical(audit$status, "not_evaluable")
  expect_identical(is.na(audit$estimate), TRUE)
  expect_equal(audit$comparable_pairs, 0, tolerance = 0)
})

test_that("average precision is tie-invariant and positive-weighted by fold", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:10),
    fold = c(rep(1L, 4L), rep(2L, 6L)),
    outcome = c(1L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L),
    probability = c(0.9, 0.8, 0.8, 0.1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4),
    status = "success",
    stringsAsFactors = FALSE
  )

  metric <- sb_binary_average_precision_metric(predictions)

  expect_identical(metric$metric, "average_precision")
  expect_identical(metric$scope, "within_outer_assessment_fold")
  expect_identical(
    metric$aggregation,
    "positive_weighted_mean_of_fold_average_precision"
  )
  expect_identical(metric$interpolation, "none")
  expect_identical(metric$status, "success")
  expect_identical(metric$n, 10L)
  expect_identical(metric$n_scored, 10L)
  expect_identical(metric$n_positive, 5L)
  expect_identical(metric$n_negative, 5L)
  expect_equal(
    metric$per_fold$estimate,
    c(5 / 6, 8 / 15),
    tolerance = 1e-15
  )
  expect_equal(
    metric$per_fold$positive_weighted_sum,
    c(5 / 3, 8 / 5),
    tolerance = 1e-15
  )
  expect_equal(metric$positive_weighted_sum, 49 / 15, tolerance = 1e-15)
  expect_equal(metric$estimate, 49 / 75, tolerance = 1e-15)
  expect_false(isTRUE(all.equal(
    metric$estimate,
    mean(metric$per_fold$estimate),
    tolerance = 0
  )))
  expect_length(metric$curves, 2L)
  expect_equal(metric$curves[[1L]]$threshold, c(0.9, 0.8, 0.1))
  expect_identical(metric$curves[[1L]]$n_at_threshold, c(1L, 2L, 1L))
  expect_identical(metric$curves[[1L]]$positive_at_threshold, c(1L, 1L, 0L))
  expect_identical(metric$curves[[1L]]$negative_at_threshold, c(0L, 1L, 1L))
  expect_equal(metric$curves[[1L]]$precision, c(1, 2 / 3, 1 / 2))
  expect_equal(metric$curves[[1L]]$recall, c(1 / 2, 1, 1))
  expect_equal(metric$curves[[1L]]$recall_increment, c(1 / 2, 1 / 2, 0))
  expect_equal(metric$curves[[1L]]$contribution, c(1 / 2, 1 / 3, 0))

  permuted <- predictions[c(4L, 3L, 2L, 1L, 10L, 6L, 9L, 5L, 8L, 7L), ]
  permuted_metric <- sb_binary_average_precision_metric(permuted)
  expect_identical(permuted_metric$per_fold, metric$per_fold)
  expect_identical(permuted_metric$curves, metric$curves)
  expect_identical(permuted_metric$estimate, metric$estimate)
})

test_that("average precision never ranks probabilities from different folds", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:4),
    fold = rep(1:2, each = 2L),
    outcome = rep(c(1L, 0L), 2L),
    probability = c(0.6, 0.5, 0.2, 0.1),
    status = "success",
    stringsAsFactors = FALSE
  )

  metric <- sb_binary_average_precision_metric(predictions)

  expect_identical(metric$status, "success")
  expect_equal(metric$per_fold$estimate, c(1, 1), tolerance = 0)
  expect_equal(metric$estimate, 1, tolerance = 0)

  pooled_curve <- sb_binary_precision_recall_curve(
    predictions$outcome,
    predictions$probability,
    fold_id = 1L
  )
  pooled_average_precision <- sum(pooled_curve$contribution)
  expect_equal(pooled_average_precision, 5 / 6, tolerance = 1e-15)
  expect_false(isTRUE(all.equal(
    metric$estimate,
    pooled_average_precision,
    tolerance = 0
  )))
})

test_that("incomplete folds cannot produce headline average precision", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:8),
    fold = rep(1:2, each = 4L),
    outcome = c(1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L),
    probability = c(0.9, 0.8, 0.4, 0.1, NA, NA, NA, NA),
    status = c(rep("success", 4L), rep("failed", 4L)),
    stringsAsFactors = FALSE
  )

  metric <- sb_binary_average_precision_metric(predictions)

  expect_identical(metric$status, "incomplete")
  expect_identical(is.na(metric$estimate), TRUE)
  expect_identical(metric$n, 8L)
  expect_identical(metric$n_scored, 4L)
  expect_identical(metric$n_evaluable_folds, 1L)
  expect_identical(
    metric$per_fold$status,
    c("success", "prediction_failed")
  )
  expect_equal(metric$positive_weighted_sum, 5 / 3, tolerance = 1e-15)
  expect_identical(nrow(metric$curves[[2L]]), 0L)

  out_of_range <- predictions
  out_of_range$status <- "success"
  out_of_range$probability[5:8] <- c(1.1, 0.7, 0.4, 0.2)
  invalid_metric <- sb_binary_average_precision_metric(out_of_range)
  expect_identical(invalid_metric$status, "incomplete")
  expect_identical(
    invalid_metric$per_fold$status,
    c("success", "prediction_failed")
  )

  one_class <- predictions
  one_class$status <- "success"
  one_class$probability[5:8] <- c(0.9, 0.7, 0.4, 0.2)
  one_class$outcome[5:8] <- 1L
  one_class_metric <- sb_binary_average_precision_metric(one_class)
  expect_identical(one_class_metric$status, "incomplete")
  expect_identical(
    one_class_metric$per_fold$status,
    c("success", "not_evaluable")
  )
})

test_that("OOF Brier score is patient-weighted across unequal folds", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:5),
    fold = c(1L, 1L, 2L, 2L, 2L),
    outcome = c(1L, 0L, 1L, 0L, 1L),
    probability = c(0.8, 0.3, 0.6, 0.4, 0.9),
    status = "success",
    stringsAsFactors = FALSE
  )
  squared_error <- (predictions$outcome - predictions$probability)^2

  metric <- sb_binary_brier_metric(predictions)

  expect_identical(metric$metric, "brier_score")
  expect_identical(metric$scope, "patient_level_out_of_fold")
  expect_identical(
    metric$aggregation,
    "summed_squared_error_over_patients"
  )
  expect_identical(metric$status, "success")
  expect_identical(metric$n, 5L)
  expect_identical(metric$n_scored, 5L)
  expect_identical(metric$n_folds, 2L)
  expect_identical(metric$n_evaluable_folds, 2L)
  expect_equal(metric$squared_error_sum, sum(squared_error), tolerance = 0)
  expect_equal(metric$estimate, mean(squared_error), tolerance = 0)
  expect_identical(metric$per_fold$n, c(2L, 3L))
  expect_equal(
    metric$per_fold$squared_error_sum,
    c(sum(squared_error[1:2]), sum(squared_error[3:5])),
    tolerance = 0
  )
  expect_equal(
    metric$per_fold$estimate,
    c(mean(squared_error[1:2]), mean(squared_error[3:5])),
    tolerance = 0
  )
  expect_false(isTRUE(all.equal(
    metric$estimate,
    mean(metric$per_fold$estimate),
    tolerance = 0
  )))
})

test_that("incomplete OOF probabilities cannot produce a headline Brier score", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:4),
    fold = rep(1:2, each = 2),
    outcome = c(1L, 0L, 1L, 0L),
    probability = c(0.8, 0.2, NA_real_, NA_real_),
    status = c("success", "success", "failed", "failed"),
    stringsAsFactors = FALSE
  )

  metric <- sb_binary_brier_metric(predictions)

  expect_identical(metric$status, "incomplete")
  expect_identical(is.na(metric$estimate), TRUE)
  expect_identical(metric$n, 4L)
  expect_identical(metric$n_scored, 2L)
  expect_equal(metric$squared_error_sum, 0.08, tolerance = 1e-15)
  expect_identical(
    metric$per_fold$status,
    c("success", "prediction_failed")
  )

  out_of_range <- predictions
  out_of_range$status <- "success"
  out_of_range$probability[3:4] <- c(1.1, 0.4)
  invalid_metric <- sb_binary_brier_metric(out_of_range)
  expect_identical(invalid_metric$status, "incomplete")
  expect_identical(is.na(invalid_metric$estimate), TRUE)
  expect_identical(
    invalid_metric$per_fold$status,
    c("success", "prediction_failed")
  )
})

test_that("fold-local prevalence reference yields patient-weighted Brier skill", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:5),
    fold = c(1L, 1L, 2L, 2L, 2L),
    outcome = c(1L, 0L, 1L, 0L, 1L),
    probability = c(0.8, 0.3, 0.6, 0.4, 0.9),
    status = "success",
    stringsAsFactors = FALSE
  )
  reference_input <- data.frame(
    fold = 1:2,
    analysis_prevalence = c(0.25, 0.75)
  )
  reference_probability <- c(0.25, 0.25, 0.75, 0.75, 0.75)
  model_loss <- (predictions$outcome - predictions$probability)^2
  reference_loss <- (predictions$outcome - reference_probability)^2

  metric <- sb_binary_brier_metric(predictions, reference_input)
  reference <- metric$reference
  skill <- metric$skill_score

  expect_identical(reference$type, "outer_analysis_prevalence")
  expect_identical(reference$source, "outer_analysis_outcome_only")
  expect_identical(reference$status, "success")
  expect_identical(reference$n, 5L)
  expect_identical(reference$n_scored, 5L)
  expect_equal(
    reference$squared_error_sum,
    sum(reference_loss),
    tolerance = 0
  )
  expect_equal(reference$estimate, mean(reference_loss), tolerance = 0)
  expect_identical(
    reference$per_fold$analysis_prevalence,
    c(0.25, 0.75)
  )
  expect_equal(
    reference$per_fold$squared_error_sum,
    c(sum(reference_loss[1:2]), sum(reference_loss[3:5])),
    tolerance = 0
  )
  expect_identical(skill$metric, "brier_skill_score")
  expect_identical(skill$status, "success")
  expect_equal(
    skill$estimate,
    1 - sum(model_loss) / sum(reference_loss),
    tolerance = 0
  )
})

test_that("Brier skill is unavailable when model OOF is incomplete", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:4),
    fold = rep(1:2, each = 2),
    outcome = c(1L, 0L, 1L, 0L),
    probability = c(0.8, 0.2, NA_real_, NA_real_),
    status = c("success", "success", "failed", "failed"),
    stringsAsFactors = FALSE
  )
  reference_input <- data.frame(
    fold = 1:2,
    analysis_prevalence = c(0.4, 0.6)
  )

  metric <- sb_binary_brier_metric(predictions, reference_input)

  expect_identical(metric$status, "incomplete")
  expect_identical(metric$reference$status, "success")
  expect_identical(metric$reference$n_scored, 4L)
  expect_identical(metric$skill_score$status, "incomplete")
  expect_identical(is.na(metric$skill_score$estimate), TRUE)
})

test_that("OOF log loss is patient-weighted with a fold-local reference", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:5),
    fold = c(1L, 1L, 2L, 2L, 2L),
    outcome = c(1L, 0L, 1L, 0L, 1L),
    probability = c(0.8, 0.25, 0.6, 0.4, 0.9),
    status = "success",
    stringsAsFactors = FALSE
  )
  reference_input <- data.frame(
    fold = 1:2,
    analysis_prevalence = c(0.5, 0.4)
  )
  case_loss <- c(
    -log(0.8), -log(0.75), -log(0.6), -log(0.6), -log(0.9)
  )
  reference_case_loss <- c(
    -log(0.5), -log(0.5), -log(0.4), -log(0.6), -log(0.4)
  )

  metric <- sb_binary_log_loss_metric(predictions, reference_input)

  expect_identical(metric$metric, "log_loss")
  expect_identical(metric$scope, "patient_level_out_of_fold")
  expect_identical(metric$aggregation, "summed_case_loss_over_patients")
  expect_identical(metric$probability_clipping, "none")
  expect_identical(metric$status, "success")
  expect_identical(metric$n, 5L)
  expect_identical(metric$n_scored, 5L)
  expect_identical(metric$n_finite_loss, 5L)
  expect_identical(metric$n_infinite_loss, 0L)
  expect_equal(metric$loss_sum, sum(case_loss), tolerance = 1e-15)
  expect_equal(metric$estimate, mean(case_loss), tolerance = 1e-15)
  expect_identical(metric$per_fold$n, c(2L, 3L))
  expect_equal(
    metric$per_fold$loss_sum,
    c(sum(case_loss[1:2]), sum(case_loss[3:5])),
    tolerance = 1e-15
  )
  expect_equal(
    metric$per_fold$estimate,
    c(mean(case_loss[1:2]), mean(case_loss[3:5])),
    tolerance = 1e-15
  )
  expect_false(isTRUE(all.equal(
    metric$estimate,
    mean(metric$per_fold$estimate),
    tolerance = 0
  )))
  expect_identical(metric$reference$status, "success")
  expect_equal(
    metric$reference$loss_sum,
    sum(reference_case_loss),
    tolerance = 1e-15
  )
  expect_equal(
    metric$reference$estimate,
    mean(reference_case_loss),
    tolerance = 1e-15
  )
  expect_equal(
    metric$reference$per_fold$analysis_prevalence,
    c(0.5, 0.4),
    tolerance = 0
  )
  expect_identical(metric$skill_score$status, "success")
  expect_equal(
    metric$skill_score$estimate,
    1 - sum(case_loss) / sum(reference_case_loss),
    tolerance = 1e-15
  )
})

test_that("log loss does not clip exact probability boundaries", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:4),
    fold = 1L,
    outcome = c(1L, 0L, 1L, 0L),
    probability = c(1, 0, 0, 1),
    status = "success",
    stringsAsFactors = FALSE
  )
  reference_input <- data.frame(
    fold = 1L,
    analysis_prevalence = 0.5
  )

  metric <- sb_binary_log_loss_metric(predictions, reference_input)

  expect_identical(metric$status, "success")
  expect_identical(metric$n_finite_loss, 2L)
  expect_identical(metric$n_infinite_loss, 2L)
  expect_identical(metric$per_fold$n_finite_loss, 2L)
  expect_identical(metric$per_fold$n_infinite_loss, 2L)
  expect_identical(metric$loss_sum, Inf)
  expect_identical(metric$estimate, Inf)
  expect_identical(metric$reference$status, "success")
  expect_identical(is.finite(metric$reference$estimate), TRUE)
  expect_identical(metric$skill_score$status, "success")
  expect_identical(metric$skill_score$estimate, -Inf)
})

test_that("failed log-loss folds retain evidence without a headline", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:4),
    fold = rep(1:2, each = 2L),
    outcome = rep(c(1L, 0L), 2L),
    probability = c(0.8, 0.2, NA_real_, NA_real_),
    status = c("success", "success", "failed", "failed"),
    stringsAsFactors = FALSE
  )
  reference_input <- data.frame(
    fold = 1:2,
    analysis_prevalence = c(0.5, 0.5)
  )

  metric <- sb_binary_log_loss_metric(predictions, reference_input)

  expect_identical(metric$status, "incomplete")
  expect_identical(is.na(metric$estimate), TRUE)
  expect_identical(metric$n_scored, 2L)
  expect_equal(metric$loss_sum, -2 * log(0.8), tolerance = 1e-15)
  expect_identical(
    metric$per_fold$status,
    c("success", "prediction_failed")
  )
  expect_identical(metric$reference$status, "success")
  expect_identical(metric$skill_score$status, "incomplete")
  expect_identical(is.na(metric$skill_score$estimate), TRUE)

  out_of_range <- predictions
  out_of_range$status <- "success"
  out_of_range$probability[3:4] <- c(1.1, 0.2)
  invalid_metric <- sb_binary_log_loss_metric(out_of_range, reference_input)
  expect_identical(invalid_metric$status, "incomplete")
  expect_identical(
    invalid_metric$per_fold$status,
    c("success", "prediction_failed")
  )

  non_finite <- out_of_range
  non_finite$probability[3:4] <- c(Inf, 0.2)
  non_finite_metric <- sb_binary_log_loss_metric(
    non_finite,
    reference_input
  )
  expect_identical(non_finite_metric$status, "incomplete")
  expect_identical(
    non_finite_metric$per_fold$status,
    c("success", "prediction_failed")
  )

  complete_predictions <- predictions
  complete_predictions$status <- "success"
  complete_predictions$probability[3:4] <- c(0.7, 0.3)
  missing_reference <- sb_binary_log_loss_metric(
    complete_predictions,
    reference_input[1L, , drop = FALSE]
  )
  expect_identical(missing_reference$status, "success")
  expect_identical(missing_reference$reference$status, "incomplete")
  expect_identical(missing_reference$skill_score$status, "incomplete")

  no_reference <- sb_binary_log_loss_metric(complete_predictions)
  expect_identical(no_reference$status, "success")
  expect_identical(no_reference$reference$status, "not_requested")
  expect_identical(no_reference$skill_score$status, "incomplete")
})

test_that("OOF calibration intercept and slope match direct logistic fits", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:20),
    fold = rep(1:5, each = 4L),
    outcome = c(
      1L, 1L, 0L, 1L,
      0L, 0L, 1L, 1L,
      1L, 0L, 1L, 0L,
      0L, 1L, 1L, 0L,
      1L, 0L, 1L, 0L
    ),
    probability = seq(0.08, 0.92, length.out = 20L),
    status = "success",
    stringsAsFactors = FALSE
  )
  logit_probability <- stats::qlogis(predictions$probability)
  direct_intercept <- stats::glm(
    predictions$outcome ~ 1,
    family = stats::binomial(),
    offset = logit_probability
  )
  direct_slope <- stats::glm(
    predictions$outcome ~ logit_probability,
    family = stats::binomial()
  )

  metric <- sb_binary_calibration_metric(predictions)

  expect_identical(metric$metric, "calibration_intercept_slope")
  expect_identical(metric$scope, "pooled_patient_level_out_of_fold")
  expect_identical(metric$purpose, "descriptive_evaluation_only")
  expect_identical(metric$probability_clipping, "none")
  expect_identical(metric$status, "success")
  expect_identical(metric$input_status, "success")
  expect_identical(metric$n, 20L)
  expect_identical(metric$n_scored, 20L)
  expect_identical(metric$n_positive, 11L)
  expect_identical(metric$n_negative, 9L)
  expect_identical(metric$n_boundary_probability, 0L)
  expect_identical(metric$intercept_target, 0)
  expect_identical(metric$slope_target, 1)
  expect_equal(
    metric$calibration_intercept,
    unname(stats::coef(direct_intercept)[[1L]]),
    tolerance = 1e-12
  )
  expect_equal(
    metric$calibration_slope,
    unname(stats::coef(direct_slope)[[2L]]),
    tolerance = 1e-12
  )
  expect_identical(metric$intercept_fit$status, "success")
  expect_identical(metric$slope_fit$status, "success")
  expect_identical(metric$intercept_fit$converged, TRUE)
  expect_identical(metric$slope_fit$converged, TRUE)
  expect_equal(
    metric$slope_fit$coefficients$estimate,
    unname(stats::coef(direct_slope)),
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(
    metric$calibration_intercept,
    metric$slope_fit$coefficients$estimate[[1L]],
    tolerance = 1e-6
  )))
  expect_identical(metric$per_fold$n, rep(4L, 5L))
  expect_identical(metric$per_fold$n_positive, c(3L, rep(2L, 4L)))
  expect_identical(metric$per_fold$n_negative, c(1L, rep(2L, 4L)))

  permuted <- predictions[c(20:1), ]
  permuted_metric <- sb_binary_calibration_metric(permuted)
  expect_equal(
    permuted_metric$calibration_intercept,
    metric$calibration_intercept,
    tolerance = 1e-12
  )
  expect_equal(
    permuted_metric$calibration_slope,
    metric$calibration_slope,
    tolerance = 1e-12
  )
  relabeled <- predictions
  relabeled$fold <- 6L - relabeled$fold
  relabeled_metric <- sb_binary_calibration_metric(relabeled)
  expect_equal(
    relabeled_metric$calibration_intercept,
    metric$calibration_intercept,
    tolerance = 1e-12
  )
  expect_equal(
    relabeled_metric$calibration_slope,
    metric$calibration_slope,
    tolerance = 1e-12
  )
})

test_that("empirically ideal probability groups hit calibration targets", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:30),
    fold = rep(1:5, each = 6L),
    outcome = c(
      rep(1L, 2L), rep(0L, 8L),
      rep(1L, 5L), rep(0L, 5L),
      rep(1L, 8L), rep(0L, 2L)
    ),
    probability = rep(c(0.2, 0.5, 0.8), each = 10L),
    status = "success",
    stringsAsFactors = FALSE
  )

  metric <- sb_binary_calibration_metric(predictions)

  expect_identical(metric$status, "success")
  expect_equal(metric$calibration_intercept, 0, tolerance = 1e-10)
  expect_equal(metric$calibration_slope, 1, tolerance = 1e-10)

  logit_probability <- stats::qlogis(predictions$probability)
  recalibration_intercept <-
    metric$slope_fit$coefficients$estimate[[1L]]
  recalibrated <- stats::plogis(
    recalibration_intercept + metric$calibration_slope * logit_probability
  )
  expect_equal(
    sum(predictions$outcome - recalibrated),
    0,
    tolerance = 1e-10
  )
  expect_equal(
    sum((predictions$outcome - recalibrated) * logit_probability),
    0,
    tolerance = 1e-10
  )
})

test_that("OOF calibration keeps non-estimable inputs explicit", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:8),
    fold = rep(1:2, each = 4L),
    outcome = rep(c(0L, 1L), 4L),
    probability = seq(0.15, 0.85, length.out = 8L),
    status = "success",
    stringsAsFactors = FALSE
  )

  failed <- predictions
  failed$status[5:8] <- "failed"
  failed$probability[5:8] <- NA_real_
  failed_metric <- sb_binary_calibration_metric(failed)
  expect_identical(failed_metric$status, "incomplete")
  expect_identical(failed_metric$input_status, "prediction_failed")
  expect_identical(failed_metric$n_scored, 4L)
  expect_identical(
    failed_metric$per_fold$status,
    c("success", "prediction_failed")
  )
  expect_identical(failed_metric$intercept_fit$status, "not_fitted")
  expect_identical(failed_metric$slope_fit$status, "not_fitted")
  expect_identical(is.na(failed_metric$calibration_intercept), TRUE)
  expect_identical(is.na(failed_metric$calibration_slope), TRUE)

  boundary <- predictions
  boundary$probability[[1L]] <- 0
  boundary_metric <- sb_binary_calibration_metric(boundary)
  expect_identical(boundary_metric$status, "incomplete")
  expect_identical(boundary_metric$input_status, "not_evaluable")
  expect_identical(boundary_metric$n_scored, 8L)
  expect_identical(boundary_metric$n_boundary_probability, 1L)
  expect_identical(
    boundary_metric$error_class,
    "stormbreaker_binary_calibration_boundary_probability"
  )
  expect_match(
    boundary_metric$error_message,
    "without probability clipping",
    fixed = TRUE
  )

  out_of_range <- predictions
  out_of_range$probability[[1L]] <- 1.1
  invalid_metric <- sb_binary_calibration_metric(out_of_range)
  expect_identical(invalid_metric$status, "incomplete")
  expect_identical(invalid_metric$input_status, "prediction_failed")

  one_class <- predictions
  one_class$outcome <- 0L
  one_class_metric <- sb_binary_calibration_metric(one_class)
  expect_identical(one_class_metric$status, "incomplete")
  expect_identical(one_class_metric$input_status, "not_evaluable")
  expect_match(
    one_class_metric$error_message,
    "Both outcome classes",
    fixed = TRUE
  )

  constant_probability <- predictions
  constant_probability$probability <- 0.5
  constant_metric <- sb_binary_calibration_metric(constant_probability)
  expect_identical(constant_metric$status, "incomplete")
  expect_identical(constant_metric$input_status, "success")
  expect_identical(constant_metric$intercept_fit$status, "success")
  expect_identical(constant_metric$slope_fit$status, "failed")
  expect_identical(
    constant_metric$slope_fit$error_class,
    "stormbreaker_binary_calibration_nonfinite_coefficients"
  )
  expect_identical(
    constant_metric$error_class,
    "stormbreaker_binary_calibration_nonfinite_coefficients"
  )
  expect_match(
    constant_metric$slope_fit$error_message,
    "non-finite coefficients",
    fixed = TRUE
  )

  separated <- predictions
  separated$outcome <- c(rep(0L, 4L), rep(1L, 4L))
  separated$probability <- c(rep(0.2, 4L), rep(0.8, 4L))
  separated_metric <- sb_binary_calibration_metric(separated)
  expect_identical(separated_metric$status, "incomplete")
  expect_identical(separated_metric$input_status, "success")
  expect_identical(separated_metric$intercept_fit$status, "success")
  expect_identical(separated_metric$slope_fit$status, "failed")
  expect_identical(
    separated_metric$slope_fit$error_class,
    "stormbreaker_calibration_separation"
  )
  expect_match(
    separated_metric$slope_fit$error_message,
    "no finite calibration slope",
    fixed = TRUE
  )
  expect_identical(is.na(separated_metric$calibration_slope), TRUE)

  reversed_separation <- separated
  reversed_separation$probability <- 1 - reversed_separation$probability
  reversed_metric <- sb_binary_calibration_metric(reversed_separation)
  expect_identical(reversed_metric$status, "incomplete")
  expect_identical(
    reversed_metric$slope_fit$error_class,
    "stormbreaker_calibration_separation"
  )

  quasi_separation <- separated
  quasi_separation$probability <- c(0.2, 0.3, 0.4, 0.5, 0.5, 0.6, 0.7, 0.8)
  quasi_metric <- sb_binary_calibration_metric(quasi_separation)
  expect_identical(quasi_metric$status, "incomplete")
  expect_identical(
    quasi_metric$slope_fit$error_class,
    "stormbreaker_calibration_separation"
  )
})
})

# ---- metrics-survival-calibration ----
local({
# ---- metrics-survival-brier ----

test_that("event KM is right-continuous and ignores censor-only jumps", {
  time <- c(1, 2, 2, 3, 4)
  event <- c(0L, 1L, 0L, 0L, 1L)
  curve <- sb_fit_event_km(time, event)
  independent <- survival::survfit(
    survival::Surv(time, event) ~ 1,
    se.fit = FALSE
  )

  expect_identical(
    names(curve),
    c(
      "time", "n_risk", "n_event", "n_censoring",
      "survival_before", "survival_after"
    )
  )
  expect_identical(curve$n_risk, c(5L, 4L, 2L, 1L))
  expect_identical(curve$n_event, c(0L, 1L, 0L, 1L))
  expect_identical(curve$n_censoring, c(1L, 1L, 1L, 0L))
  expect_equal(curve$survival_before, c(1, 1, 0.75, 0.75))
  expect_equal(curve$survival_after, c(1, 0.75, 0.75, 0))
  expect_equal(
    curve$survival_after,
    summary(independent, times = curve$time, extend = TRUE)$surv,
    tolerance = 1e-15
  )
  expect_equal(sb_event_survival_at(curve, 1), 1)
  expect_equal(sb_event_survival_at(curve, 2), 0.75)
  expect_equal(sb_event_survival_at(curve, 3), 0.75)
  expect_equal(sb_event_survival_at(curve, 4), 0)
})

test_that("IPCW Brier reduces to patient-level squared error without censoring", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:5),
    fold = c(1L, 1L, 2L, 2L, 2L),
    time = c(1, 3, 1, 3, 4),
    event = rep(1L, 5L),
    survival_probability = c(0.8, 0.4, 0.7, 0.5, 0.2),
    survival_probability_status = rep("success", 5L),
    status = rep("success", 5L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(
    list(fold = 1L, time = c(1, 2, 3, 4), event = rep(1L, 4L)),
    list(fold = 2L, time = c(1, 2, 3, 4), event = rep(1L, 4L))
  )

  observed <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 2
  )
  expected_loss <- c(0.8^2, (1 - 0.4)^2, 0.7^2,
                     (1 - 0.5)^2, (1 - 0.2)^2)

  expect_identical(observed$status, "success")
  expect_identical(observed$n, 5L)
  expect_identical(observed$n_scored, 5L)
  expect_equal(observed$weight_sum, 5, tolerance = 0)
  expect_equal(
    observed$weighted_loss_sum,
    sum(expected_loss),
    tolerance = 1e-15
  )
  expect_equal(
    observed$estimate,
    sum(expected_loss) / 5,
    tolerance = 1e-15
  )
  expect_equal(
    observed$per_fold$estimate,
    c(sum(expected_loss[1:2]) / 2, sum(expected_loss[3:5]) / 3),
    tolerance = 0
  )
  expect_gt(
    abs(observed$estimate - mean(observed$per_fold$estimate)),
    0
  )
  expect_identical(observed$reference$status, "success")
  expect_equal(
    observed$reference$per_fold$reference_survival_horizon,
    c(0.5, 0.5),
    tolerance = 0
  )
  expect_equal(observed$reference$weighted_loss_sum, 1.25, tolerance = 0)
  expect_equal(observed$reference$estimate, 0.25, tolerance = 0)
  expect_identical(observed$skill_score$status, "success")
  expect_equal(
    observed$skill_score$estimate,
    1 - sum(expected_loss) / 1.25,
    tolerance = 1e-15
  )
})

test_that("survival Brier reference and skill aggregate patient losses", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:5),
    fold = c(1L, 1L, 2L, 2L, 2L),
    time = c(1, 3, 1, 3, 4),
    event = rep(1L, 5L),
    survival_probability = c(0.8, 0.4, 0.7, 0.5, 0.2),
    survival_probability_status = rep("success", 5L),
    status = rep("success", 5L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(
    list(fold = 1L, time = c(1, 2, 3, 4), event = rep(1L, 4L)),
    list(
      fold = 2L,
      time = c(1, 2, 3, 4),
      event = c(1L, 0L, 1L, 1L)
    )
  )
  model_loss <- c(0.8^2, (1 - 0.4)^2, 0.7^2,
                  (1 - 0.5)^2, (1 - 0.2)^2)
  reference_loss <- c(0.5^2, (1 - 0.5)^2, 0.75^2,
                      (1 - 0.75)^2, (1 - 0.75)^2)

  observed <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 2
  )
  reference <- observed$reference

  expect_identical(reference$type, "outer_analysis_kaplan_meier")
  expect_identical(reference$source, "outer_analysis_event_only")
  expect_identical(reference$status, "success")
  expect_identical(reference$n, 5L)
  expect_identical(reference$n_scored, 5L)
  expect_equal(
    reference$per_fold$reference_survival_horizon,
    c(0.5, 0.75),
    tolerance = 0
  )
  expect_equal(
    reference$per_fold$weighted_loss_sum,
    c(sum(reference_loss[1:2]), sum(reference_loss[3:5])),
    tolerance = 0
  )
  expect_equal(reference$weighted_loss_sum, sum(reference_loss), tolerance = 0)
  expect_equal(reference$estimate, sum(reference_loss) / 5, tolerance = 0)
  expect_gt(
    abs(reference$estimate - mean(reference$per_fold$estimate)),
    0
  )
  expect_identical(observed$skill_score$metric, "ipcw_brier_skill_score")
  expect_identical(observed$skill_score$status, "success")
  expect_equal(
    observed$skill_score$estimate,
    1 - sum(model_loss) / sum(reference_loss),
    tolerance = 1e-15
  )
})

test_that("IPCW Brier uses analysis-fold G(T-) and G(t)", {
  predictions <- data.frame(
    patient_id = c("event", "survivor", "censored"),
    fold = rep(1L, 3L),
    time = c(2, 4, 2.5),
    event = c(1L, 1L, 0L),
    survival_probability = c(0.6, 0.7, 0.5),
    survival_probability_status = rep("success", 3L),
    status = rep("success", 3L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(list(
    fold = 1L,
    time = c(1, 2, 2, 3, 4),
    event = c(0L, 1L, 0L, 0L, 1L)
  ))

  observed <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 3
  )
  curve <- observed$censoring_curves[[1L]]
  expected_event_loss <- (1 / 0.8) * 0.6^2
  expected_survivor_loss <- (1 / (0.8 * 2 / 3)) * (1 - 0.7)^2

  expect_identical(
    names(curve),
    c(
      "time", "n_risk", "n_event", "n_censoring",
      "n_risk_for_censoring", "survival_before", "survival_after"
    )
  )
  expect_equal(curve$survival_before, c(1, 0.8, 0.8 * 2 / 3, 0.8 / 3))
  expect_equal(curve$survival_after, c(0.8, 0.8 * 2 / 3, 0.8 / 3, 0.8 / 3))
  expect_identical(observed$per_fold$n_events_by_horizon, 1L)
  expect_identical(observed$per_fold$n_surviving_past_horizon, 1L)
  expect_identical(observed$per_fold$n_censored_by_horizon, 1L)
  expect_equal(
    observed$per_fold$censoring_survival_horizon,
    0.8 * 2 / 3
  )
  expect_equal(observed$weight_sum, 1 / 0.8 + 1 / (0.8 * 2 / 3))
  expect_equal(
    observed$weighted_loss_sum,
    expected_event_loss + expected_survivor_loss
  )
  expect_equal(
    observed$estimate,
    (expected_event_loss + expected_survivor_loss) / 3
  )
})

test_that("reverse KM orders tied events before censoring", {
  analysis_time <- c(1, 2, 2, 3, 4)
  analysis_event <- c(0L, 1L, 0L, 0L, 1L)

  curve <- sb_fit_censoring_km(analysis_time, analysis_event)
  censor_rows <- curve$n_censoring > 0L
  shifted_time <- analysis_time + ifelse(analysis_event == 0L, 0.25, 0)
  independent <- survival::survfit(
    survival::Surv(shifted_time, 1L - analysis_event) ~ 1,
    se.fit = FALSE
  )
  independent_after <- summary(
    independent,
    times = c(1.25, 2.25, 3.25),
    extend = TRUE
  )$surv

  expect_identical(curve$n_event, c(0L, 1L, 0L, 1L))
  expect_identical(curve$n_risk_for_censoring, c(5L, 3L, 2L, 0L))
  expect_equal(
    curve$survival_after[censor_rows],
    independent_after,
    tolerance = 1e-15
  )
  expect_equal(
    sb_censoring_survival_at(curve, 2, left_limit = TRUE),
    0.8,
    tolerance = 1e-15
  )
  expect_equal(
    sb_censoring_survival_at(curve, 3, left_limit = TRUE),
    0.8 * 2 / 3,
    tolerance = 1e-15
  )
})

test_that("censoring at the horizon has a known survivor outcome", {
  predictions <- data.frame(
    patient_id = c("event_at", "censored_at", "survivor"),
    fold = rep(1L, 3L),
    time = c(3, 3, 4),
    event = c(1L, 0L, 1L),
    survival_probability = c(0.6, 0.7, 0.8),
    survival_probability_status = rep("success", 3L),
    status = rep("success", 3L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(list(
    fold = 1L,
    time = c(1, 2, 2, 3, 4),
    event = c(0L, 1L, 0L, 0L, 1L)
  ))

  observed <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 3
  )
  censoring_survival <- 0.8 * 2 / 3
  expected_loss <- (
    0.6^2 + (1 - 0.7)^2 + (1 - 0.8)^2
  ) / censoring_survival

  expect_identical(observed$status, "success")
  expect_identical(observed$n_scored, 3L)
  expect_equal(observed$weight_sum, 3 / censoring_survival)
  expect_equal(observed$weighted_loss_sum, expected_loss)
  expect_equal(observed$estimate, expected_loss / 3)
})

test_that("left-continuous G retains support at a censoring horizon", {
  predictions <- data.frame(
    patient_id = c("P1", "P2"),
    fold = c(1L, 1L),
    time = c(1, 3),
    event = c(1L, 1L),
    survival_probability = c(0.6, 0.7),
    survival_probability_status = c("success", "success"),
    status = c("success", "success"),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(list(
    fold = 1L,
    time = c(1, 2, 3),
    event = c(1L, 0L, 0L)
  ))

  observed <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 3
  )

  expect_identical(observed$status, "success")
  expect_equal(observed$per_fold$censoring_survival_horizon, 0.5)
  expect_equal(observed$weight_sum, 3)
  expect_equal(observed$weighted_loss_sum, 0.6^2 + 2 * 0.7^2)
  expect_equal(observed$estimate, (0.6^2 + 2 * 0.7^2) / 2)
  expect_length(observed$censoring_curves, 1L)
})

test_that("invalid censoring survival withholds model, reference, and skill", {
  predictions <- data.frame(
    patient_id = c("P1", "P2"),
    fold = rep(1L, 2L),
    time = c(1, 3),
    event = c(1L, 1L),
    survival_probability = c(0.6, 0.7),
    survival_probability_status = rep("success", 2L),
    status = rep("success", 2L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(list(
    fold = 1L,
    time = c(1, 2),
    event = c(1L, 0L)
  ))
  forced_survival <- 0
  local_mocked_bindings(
    sb_fit_censoring_km = function(...) {
      data.frame(
        time = c(1, 2),
        n_risk = c(2L, 1L),
        n_event = c(1L, 0L),
        n_censoring = c(0L, 1L),
        n_risk_for_censoring = c(1L, 1L),
        survival_before = c(1, forced_survival),
        survival_after = c(forced_survival, forced_survival)
      )
    }
  )

  zero <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 2
  )
  forced_survival <- Inf
  non_finite <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 2
  )

  for (observed in list(zero, non_finite)) {
    expect_identical(observed$status, "incomplete")
    expect_identical(observed$per_fold$status, "censoring_failed")
    expect_identical(observed$reference$status, "incomplete")
    expect_identical(
      observed$reference$per_fold$status,
      "censoring_failed"
    )
    expect_identical(observed$skill_score$status, "incomplete")
    expect_identical(is.na(observed$skill_score$estimate), TRUE)
  }
})

test_that("failed probability folds retain partial IPCW audit evidence", {
  predictions <- data.frame(
    patient_id = c("P1", "P2", "P3", "P4"),
    fold = c(1L, 1L, 2L, 2L),
    time = c(1, 3, 1, 3),
    event = rep(1L, 4L),
    survival_probability = c(0.6, 0.7, NA, NA),
    survival_probability_status = c(
      "success", "success", "failed", "failed"
    ),
    status = rep("success", 4L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(
    list(fold = 1L, time = c(1, 2, 3), event = c(1L, 0L, 1L)),
    list(fold = 2L, time = c(1, 2, 3), event = c(1L, 0L, 1L))
  )

  observed <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 2
  )

  expect_identical(observed$status, "incomplete")
  expect_true(is.na(observed$estimate))
  expect_identical(
    observed$per_fold$status,
    c("success", "prediction_failed")
  )
  expect_gt(observed$weighted_loss_sum, 0)
  expect_identical(observed$reference$status, "success")
  expect_identical(observed$reference$n_scored, 4L)
  expect_identical(observed$skill_score$status, "incomplete")
  expect_identical(is.na(observed$skill_score$estimate), TRUE)
})

test_that("boundary reference probabilities are valid but zero loss has no skill", {
  survivor_predictions <- data.frame(
    patient_id = c("P1", "P2"),
    fold = rep(1L, 2L),
    time = c(3, 4),
    event = c(1L, 0L),
    survival_probability = c(0.8, 0.7),
    survival_probability_status = rep("success", 2L),
    status = rep("success", 2L),
    stringsAsFactors = FALSE
  )
  survivor_analysis <- list(list(
    fold = 1L,
    time = c(3, 4),
    event = c(1L, 0L)
  ))

  reference_one <- sb_survival_ipcw_brier_metric(
    survivor_predictions,
    survivor_analysis,
    prediction_horizon = 2
  )

  expect_identical(reference_one$reference$status, "success")
  expect_identical(
    reference_one$reference$per_fold$reference_survival_horizon,
    1
  )
  expect_identical(reference_one$reference$weighted_loss_sum, 0)
  expect_identical(reference_one$reference$estimate, 0)
  expect_identical(reference_one$skill_score$status, "incomplete")
  expect_identical(is.na(reference_one$skill_score$estimate), TRUE)

  event_predictions <- survivor_predictions[1L, , drop = FALSE]
  event_predictions$time <- 3
  event_predictions$event <- 1L
  event_predictions$survival_probability <- 0.5
  event_analysis <- list(list(
    fold = 1L,
    time = c(1, 2),
    event = c(1L, 1L)
  ))
  reference_zero <- sb_survival_ipcw_brier_metric(
    event_predictions,
    event_analysis,
    prediction_horizon = 2
  )

  expect_identical(reference_zero$reference$status, "success")
  expect_identical(
    reference_zero$reference$per_fold$reference_survival_horizon,
    0
  )
  expect_identical(reference_zero$reference$estimate, 1)
  expect_identical(reference_zero$skill_score$status, "success")
  expect_identical(reference_zero$skill_score$estimate, 0.75)
})

# ---- metrics-survival-calibration ----

make_sb_survival_calibration_fixture <- function(
    n = 160L,
    seed = 20260820L,
    horizon = 2) {
  withr::with_seed(seed, {
    prognostic_index <- stats::rnorm(n)
    event_time <- stats::rexp(
      n,
      rate = 0.15 * exp(0.45 * prognostic_index)
    )
    censoring_time <- stats::rexp(n, rate = 0.05)
    time <- pmin(event_time, censoring_time)
    event <- as.integer(event_time <= censoring_time)
    survival_probability <- exp(
      -0.15 * horizon * exp(0.45 * prognostic_index)
    )
    data.frame(
      patient_id = sprintf("C%03d", seq_len(n)),
      time = time,
      event = event,
      survival_probability = survival_probability,
      status = rep("success", n),
      fold = rep(seq_len(5L), length.out = n),
      stringsAsFactors = FALSE
    )
  })
}

manual_three_knot_rcs_basis <- function(x, knots) {
  knot1 <- knots[[1L]]
  knot2 <- knots[[2L]]
  knot3 <- knots[[3L]]
  scale <- (knot3 - knot1)^(2 / 3)
  nonlinear <- pmax((x - knot1) / scale, 0)^3 +
    (
      (knot2 - knot1) * pmax((x - knot3) / scale, 0)^3 -
        (knot3 - knot1) * pmax((x - knot2) / scale, 0)^3
    ) / (knot3 - knot2)
  cbind(rcs_linear = x, rcs_nonlinear = nonlinear)
}

manual_survival_calibration_oracle <- function(data, horizon) {
  x <- log(-log(data$survival_probability))
  knots <- unname(stats::quantile(
    x,
    probs = c(0.10, 0.50, 0.90),
    names = FALSE,
    type = 7
  ))
  basis <- manual_three_knot_rcs_basis(x, knots)
  fit_data <- data.frame(
    time = data$time,
    event = data$event,
    rcs_linear = basis[, 1L],
    rcs_nonlinear = basis[, 2L]
  )
  fit <- survival::coxph(
    survival::Surv(time, event) ~ rcs_linear + rcs_nonlinear,
    data = fit_data,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
  new_data <- fit_data
  new_data$time <- rep(horizon, nrow(new_data))
  new_data$event <- rep(1L, nrow(new_data))
  cumulative_hazard_predict <- unname(stats::predict(
    fit,
    newdata = new_data,
    type = "expected"
  ))

  baseline <- survival::basehaz(fit, centered = TRUE)
  eligible <- baseline$time <= horizon
  baseline_horizon <- if (any(eligible)) {
    unname(utils::tail(baseline$hazard[eligible], 1L))
  } else {
    0
  }
  centered_basis <- sweep(basis, 2L, fit$means, FUN = "-")
  centered_lp <- drop(centered_basis %*% stats::coef(fit))
  cumulative_hazard_centered <- if (baseline_horizon == 0) {
    rep(0, nrow(data))
  } else {
    exp(log(baseline_horizon) + centered_lp)
  }
  calibrated_survival <- exp(-cumulative_hazard_predict)
  absolute_error <- abs(
    calibrated_survival - data$survival_probability
  )
  list(
    knots = knots,
    fit = fit,
    cumulative_hazard_predict = cumulative_hazard_predict,
    cumulative_hazard_centered = cumulative_hazard_centered,
    calibrated_survival_probability = calibrated_survival,
    ici = mean(absolute_error),
    e50 = unname(stats::quantile(absolute_error, 0.50, type = 7)),
    e90 = unname(stats::quantile(absolute_error, 0.90, type = 7))
  )
}

call_sb_survival_calibration <- function(data, horizon = 2, fold = data$fold) {
  sb_survival_calibration_metric(
    patient_id = data$patient_id,
    time = data$time,
    event = data$event,
    survival_probability = data$survival_probability,
    prediction_status = data$status,
    prediction_horizon = horizon,
    time_unit = "year",
    scope = "pooled_patient_level_out_of_fold",
    fold = fold
  )
}

test_that("survival calibration matches an independent three-knot RCS oracle", {
  data <- make_sb_survival_calibration_fixture()
  expected <- manual_survival_calibration_oracle(data, horizon = 2)
  observed <- call_sb_survival_calibration(data)

  expect_identical(observed$metric, "survival_calibration_curve")
  expect_identical(observed$calibration_level, "moderate")
  expect_identical(
    observed$method,
    "cox_ph_three_knot_restricted_natural_cubic_spline"
  )
  expect_identical(observed$scope, "pooled_patient_level_out_of_fold")
  expect_identical(observed$purpose, "descriptive_evaluation_only")
  expect_identical(observed$event_probability_transform, "complementary_log_log")
  expect_identical(observed$probability_clipping, "none")
  expect_identical(observed$ties, "efron")
  expect_identical(observed$baseline_reference, "coxph_mean_basis_centered")
  expect_identical(observed$status, "success")
  expect_identical(observed$input_status, "success")
  expect_equal(observed$knots$quantile_probability, c(0.10, 0.50, 0.90))
  expect_equal(observed$knots$transformed_value, expected$knots, tolerance = 0)
  expect_equal(
    expected$cumulative_hazard_predict,
    expected$cumulative_hazard_centered,
    tolerance = 1e-11
  )
  expect_equal(
    observed$patient_calibration$calibrated_survival_probability,
    expected$calibrated_survival_probability,
    tolerance = 1e-10
  )
  expect_equal(observed$ici, expected$ici, tolerance = 1e-11)
  expect_equal(observed$e50, expected$e50, tolerance = 1e-11)
  expect_equal(observed$e90, expected$e90, tolerance = 1e-11)
  expect_equal(
    observed$fit$means,
    unname(observed$fit$means),
    ignore_attr = TRUE
  )
  expect_identical(
    observed$patient_calibration$row_instance,
    seq_len(nrow(data))
  )
  expect_identical(
    observed$patient_calibration$patient_id,
    data$patient_id
  )
  expect_identical(observed$n_folds, 5L)
  expect_identical(observed$n_evaluable_folds, 5L)
})

test_that("survival calibration keeps extreme interior probabilities finite", {
  data <- make_sb_survival_calibration_fixture(n = 180L, seed = 17L)
  data$survival_probability[1:3] <- c(
    1e-20,
    1e-300,
    1 - .Machine$double.eps
  )
  observed <- call_sb_survival_calibration(data)

  expect_identical(observed$status, "success")
  expect_identical(observed$n_boundary_probability, 0L)
  expect_true(all(is.finite(
    observed$patient_calibration$transformed_prediction[1:3]
  )))
  expect_equal(
    observed$patient_calibration$transformed_prediction[1:3],
    log(-log(data$survival_probability[1:3])),
    tolerance = 0
  )
  expect_equal(
    observed$patient_calibration$predicted_event_probability[1:3],
    -expm1(log(data$survival_probability[1:3])),
    tolerance = 0
  )
  expect_true(all(is.finite(
    observed$patient_calibration$absolute_error[1:3]
  )))
})

test_that("survival calibration exposes exact boundaries without clipping", {
  data <- make_sb_survival_calibration_fixture()
  data$survival_probability[1:2] <- c(0, 1)
  observed <- call_sb_survival_calibration(data)

  expect_identical(observed$status, "incomplete")
  expect_identical(observed$input_status, "not_evaluable")
  expect_identical(
    observed$error_class,
    "stormbreaker_survival_calibration_boundary_probability"
  )
  expect_identical(observed$n_boundary_probability, 2L)
  expect_identical(nrow(observed$patient_calibration), nrow(data))
  expect_equal(
    observed$patient_calibration$survival_probability[1:2],
    c(0, 1),
    tolerance = 0
  )
  expect_true(all(is.na(c(observed$ici, observed$e50, observed$e90))))
})

test_that("survival calibration does not move duplicated knots", {
  data <- make_sb_survival_calibration_fixture()
  data$survival_probability <- rep(0.8, nrow(data))
  data$survival_probability[seq_len(5L)] <- 0.7
  observed <- call_sb_survival_calibration(data)

  expect_identical(observed$status, "incomplete")
  expect_identical(observed$input_status, "not_evaluable")
  expect_identical(
    observed$error_class,
    "stormbreaker_survival_calibration_duplicate_knots"
  )
  expect_identical(observed$n_unique_survival_probability, 2L)
  expect_identical(nrow(observed$knots), 3L)
  expect_true(anyDuplicated(observed$knots$transformed_value) > 0L)
})

test_that("survival calibration distinguishes horizon support cases", {
  data <- make_sb_survival_calibration_fixture(n = 180L, seed = 29L)
  after_all_followup <- call_sb_survival_calibration(
    data,
    horizon = max(data$time) + 1
  )
  expect_identical(after_all_followup$status, "incomplete")
  expect_identical(
    after_all_followup$error_class,
    "stormbreaker_survival_calibration_horizon_unsupported"
  )
  expect_identical(after_all_followup$n_followed_to_horizon, 0L)

  early <- data
  early$time <- early$time + 1
  before_first_event <- call_sb_survival_calibration(early, horizon = 0.5)
  expect_identical(before_first_event$status, "success")
  expect_equal(
    before_first_event$fit$baseline_cumulative_hazard_horizon,
    0,
    tolerance = 0
  )
  expect_true(all(
    before_first_event$patient_calibration$calibrated_survival_probability == 1
  ))

  after_last_event <- data
  after_last_event$event[after_last_event$time > 3] <- 0L
  after_last_event$time[[1L]] <- max(after_last_event$time) + 5
  after_last_event$event[[1L]] <- 0L
  supported <- call_sb_survival_calibration(after_last_event, horizon = 4)
  expect_identical(supported$status, "success")
  expect_gt(supported$n_followed_to_horizon, 0L)
  expect_true(max(after_last_event$time[after_last_event$event == 1L]) < 4)
})

test_that("survival calibration keeps complete follow-up after the horizon", {
  data <- make_sb_survival_calibration_fixture(n = 180L, seed = 31L)
  horizon <- 1
  changed <- data
  after <- changed$time > horizon
  changed$event[after] <- 1L - changed$event[after]
  if (sum(changed$event) == 0L) {
    changed$event[which(after)[1L]] <- 1L
  }
  original <- call_sb_survival_calibration(data, horizon = horizon)
  modified <- call_sb_survival_calibration(changed, horizon = horizon)

  expect_identical(original$status, "success")
  expect_identical(modified$status, "success")
  expect_false(isTRUE(all.equal(
    original$patient_calibration$calibrated_survival_probability,
    modified$patient_calibration$calibrated_survival_probability,
    tolerance = 1e-12
  )))
})

test_that("survival calibration retains bootstrap row instances and time zero", {
  data <- make_sb_survival_calibration_fixture(n = 180L, seed = 41L)
  data$time[1:2] <- 0
  data$event[1:2] <- c(1L, 0L)
  sampled <- data[c(1:170, 10, 10, 20, 20, 30, 30, 40, 40, 50, 50), ]
  rownames(sampled) <- NULL
  observed <- call_sb_survival_calibration(sampled)

  expect_identical(observed$status, "success")
  expect_true(anyDuplicated(sampled$patient_id) > 0L)
  expect_identical(
    observed$patient_calibration$row_instance,
    seq_len(nrow(sampled))
  )
  expect_identical(
    observed$patient_calibration$patient_id,
    sampled$patient_id
  )
  expect_identical(nrow(observed$patient_calibration), nrow(sampled))
})

test_that("survival calibration is row-order invariant after row alignment", {
  data <- make_sb_survival_calibration_fixture(n = 180L, seed = 43L)
  order <- c(seq(2L, nrow(data), by = 2L), seq(1L, nrow(data), by = 2L))
  observed <- call_sb_survival_calibration(data)
  shuffled <- call_sb_survival_calibration(data[order, ])
  aligned <- match(
    observed$patient_calibration$patient_id,
    shuffled$patient_calibration$patient_id
  )

  expect_identical(observed$status, "success")
  expect_identical(shuffled$status, "success")
  expect_equal(observed$ici, shuffled$ici, tolerance = 1e-12)
  expect_equal(observed$e50, shuffled$e50, tolerance = 1e-12)
  expect_equal(observed$e90, shuffled$e90, tolerance = 1e-12)
  expect_equal(
    observed$patient_calibration$calibrated_survival_probability,
    shuffled$patient_calibration$calibrated_survival_probability[aligned],
    tolerance = 1e-11
  )
})

test_that("survival calibration fit warnings and errors remain metric-local", {
  data <- make_sb_survival_calibration_fixture()
  original <- sb_survival_calibration_fit_cox
  testthat::local_mocked_bindings(
    sb_survival_calibration_fit_cox = function(fit_data) {
      warning("forced calibration warning")
      original(fit_data)
    }
  )
  warned <- call_sb_survival_calibration(data)
  expect_identical(warned$status, "fit_failed")
  expect_match(warned$fit$warnings, "forced calibration warning")
  expect_true(all(is.na(c(warned$ici, warned$e50, warned$e90))))
})

test_that("survival calibration error audits preserve the original condition", {
  error <- structure(
    list(message = "forced optional calibration error", call = NULL),
    class = c("forced_calibration_error", "error", "condition")
  )
  observed <- sb_survival_calibration_error_metric(
    error = error,
    prediction_horizon = 2,
    time_unit = "year",
    scope = "patient_level_external_cohort_at_locked_horizon",
    n = 25L
  )

  expect_identical(observed$status, "metric_failed")
  expect_identical(observed$input_status, "unexpected_error")
  expect_identical(observed$error_class, "forced_calibration_error")
  expect_identical(
    observed$error_message,
    "forced optional calibration error"
  )
  expect_identical(observed$n, 25L)
})

# ---- metrics-survival-decision-curve ----

make_sb_survival_dca_predictions <- function(
    time,
    event,
    event_probability,
    status = rep("success", length(time)),
    probability_status = rep("success", length(time))) {
  data.frame(
    patient_id = paste0("s", seq_along(time)),
    time = time,
    event = as.integer(event),
    survival_probability = 1 - event_probability,
    status = status,
    survival_probability_status = probability_status,
    stringsAsFactors = FALSE
  )
}

test_that("survival DCA reconstructs the published Kaplan-Meier formula", {
  predictions <- make_sb_survival_dca_predictions(
    time = 1:6,
    event = c(1L, 0L, 1L, 0L, 1L, 0L),
    event_probability = c(0.80, 0.70, 0.60, 0.40, 0.55, 0.10)
  )

  metric <- sb_survival_decision_curve_metric(
    predictions,
    decision_thresholds = 0.50,
    prediction_horizon = 4,
    time_unit = "month",
    scope = "test"
  )
  row <- metric$curve[1L, ]

  expect_identical(metric$status, "success")
  expect_identical(row$n_selected, 4L)
  expect_equal(row$selected_fraction, 4 / 6, tolerance = 0)
  expect_equal(row$selected_survival_probability, 3 / 8, tolerance = 0)
  expect_equal(row$selected_event_probability, 5 / 8, tolerance = 0)
  expect_equal(row$estimated_true_positive_rate, 5 / 12, tolerance = 1e-15)
  expect_equal(row$estimated_false_positive_rate, 1 / 4, tolerance = 0)
  expect_equal(row$model_net_benefit, 1 / 6, tolerance = 1e-15)
  expect_equal(metric$overall_survival_probability, 5 / 8, tolerance = 0)
  expect_equal(row$treat_all_net_benefit, -1 / 4, tolerance = 0)
  expect_identical(row$treat_none_net_benefit, 0)
})

test_that("survival DCA applies exact action and empty-selection boundaries", {
  predictions <- make_sb_survival_dca_predictions(
    time = c(1, 2, 3, 4),
    event = c(1L, 0L, 1L, 0L),
    event_probability = c(0.80, 0.50, 0.30, 0.10)
  )
  metric <- sb_survival_decision_curve_metric(
    predictions,
    decision_thresholds = c(0.05, 0.50, 0.90),
    prediction_horizon = 3,
    time_unit = "month",
    scope = "test"
  )

  all_row <- metric$curve[metric$curve$threshold_probability == 0.05, ]
  tied_row <- metric$curve[metric$curve$threshold_probability == 0.50, ]
  empty_row <- metric$curve[metric$curve$threshold_probability == 0.90, ]
  expect_identical(all_row$n_selected, 4L)
  expect_equal(
    all_row$model_net_benefit,
    all_row$treat_all_net_benefit,
    tolerance = 0
  )
  expect_identical(tied_row$n_selected, 2L)
  expect_identical(empty_row$n_selected, 0L)
  expect_identical(empty_row$horizon_support, "empty_selection")
  expect_identical(empty_row$selected_survival_probability, 1)
  expect_identical(empty_row$selected_event_probability, 0)
  expect_identical(empty_row$model_net_benefit, 0)
})

test_that("event-censor ties process the event before censoring", {
  predictions <- make_sb_survival_dca_predictions(
    time = c(1, 1),
    event = c(1L, 0L),
    event_probability = c(0.80, 0.80)
  )
  metric <- sb_survival_decision_curve_metric(
    predictions,
    decision_thresholds = 0.50,
    prediction_horizon = 1,
    time_unit = "month",
    scope = "test"
  )

  expect_identical(metric$status, "success")
  expect_equal(
    metric$curve$selected_survival_probability,
    0.5,
    tolerance = 0
  )
  expect_equal(metric$overall_survival_probability, 0.5, tolerance = 0)
})

test_that("zero survival is estimable beyond the last selected follow-up", {
  predictions <- make_sb_survival_dca_predictions(
    time = c(1, 2, 4),
    event = c(1L, 1L, 0L),
    event_probability = c(0.80, 0.80, 0.10)
  )
  metric <- sb_survival_decision_curve_metric(
    predictions,
    decision_thresholds = 0.50,
    prediction_horizon = 3,
    time_unit = "month",
    scope = "test"
  )

  expect_identical(metric$status, "success")
  expect_identical(
    metric$curve$horizon_support,
    "survival_zero_before_horizon"
  )
  expect_identical(metric$curve$selected_survival_probability, 0)
  expect_equal(metric$curve$model_net_benefit, 2 / 3, tolerance = 0)
  expect_equal(
    metric$curve$treat_all_net_benefit,
    1 / 3,
    tolerance = 1e-15
  )
})

test_that("positive KM survival is not extrapolated past observed support", {
  predictions <- make_sb_survival_dca_predictions(
    time = c(1, 2, 4),
    event = c(1L, 0L, 0L),
    event_probability = c(0.80, 0.80, 0.10)
  )
  metric <- sb_survival_decision_curve_metric(
    predictions,
    decision_thresholds = 0.50,
    prediction_horizon = 3,
    time_unit = "month",
    scope = "test"
  )

  expect_identical(metric$status, "incomplete")
  expect_identical(metric$n_evaluable_thresholds, 0L)
  expect_identical(metric$curve$status, "not_evaluable")
  expect_identical(metric$curve$horizon_support, "insufficient_follow_up")
  expect_true(is.na(metric$curve$model_net_benefit))
  expect_true(is.finite(metric$curve$treat_all_net_benefit))
  expect_identical(
    metric$curve$error_class,
    "stormbreaker_insufficient_survival_dca_follow_up"
  )
})

test_that("the treat-all reference uses the same zero-versus-positive rule", {
  zero <- make_sb_survival_dca_predictions(
    time = c(1, 2),
    event = c(1L, 1L),
    event_probability = c(0.10, 0.10)
  )
  positive <- make_sb_survival_dca_predictions(
    time = c(1, 2),
    event = c(1L, 0L),
    event_probability = c(0.10, 0.10)
  )
  zero_metric <- sb_survival_decision_curve_metric(
    zero, 0.50, 3, "month", "test"
  )
  positive_metric <- sb_survival_decision_curve_metric(
    positive, 0.50, 3, "month", "test"
  )

  expect_identical(zero_metric$overall_horizon_support,
                   "survival_zero_before_horizon")
  expect_identical(zero_metric$status, "success")
  expect_identical(zero_metric$curve$model_net_benefit, 0)
  expect_identical(zero_metric$curve$treat_all_net_benefit, 1)
  expect_identical(positive_metric$overall_horizon_support,
                   "insufficient_follow_up")
  expect_identical(positive_metric$status, "incomplete")
  expect_true(is.na(positive_metric$curve$treat_all_net_benefit))
})

test_that("invalid or failed horizon predictions leave DCA incomplete", {
  predictions <- make_sb_survival_dca_predictions(
    time = c(1, 2, 3),
    event = c(1L, 0L, 1L),
    event_probability = c(0, 0.5, 1),
    probability_status = c("success", "failed", "success")
  )
  metric <- sb_survival_decision_curve_metric(
    predictions, 0.50, 2, "month", "test"
  )

  expect_identical(metric$status, "incomplete")
  expect_identical(metric$curve$status, "prediction_failed")
  expect_true(is.na(metric$curve$model_net_benefit))
  expect_identical(
    metric$curve$error_class,
    "stormbreaker_invalid_survival_prediction"
  )
})

# ---- metrics-survival-weak-calibration ----

make_sb_survival_weak_calibration_fixture <- function() {
  transformed_prediction <- c(0.8, -0.4, 1.2, -0.9, 0.3, 1.0, -1.1, 0.1)
  data.frame(
    patient_id = paste0("weak_", seq_along(transformed_prediction)),
    fold = rep(1:4, each = 2L),
    time = c(0.7, 1.3, 2.1, 3.4, 4.0, 5.2, 6.1, 7.3),
    event = c(1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L),
    survival_probability = exp(-exp(transformed_prediction)),
    stringsAsFactors = FALSE
  )
}

manual_sb_event_probability_at <- function(time, event, horizon) {
  event_times <- sort(unique(time[time <= horizon & event == 1L]))
  survival_probability <- 1
  for (event_time in event_times) {
    n_risk <- sum(time >= event_time)
    n_event <- sum(time == event_time & event == 1L)
    survival_probability <- survival_probability * (1 - n_event / n_risk)
  }
  1 - survival_probability
}

manual_sb_no_tie_cox_slope <- function(time, event, predictor, horizon) {
  administrative_time <- pmin(time, horizon)
  event_by_horizon <- as.integer(event == 1L & time <= horizon)
  event_rows <- which(event_by_horizon == 1L)
  partial_log_likelihood <- function(coefficient) {
    sum(vapply(event_rows, function(row) {
      risk_set <- administrative_time >= administrative_time[[row]]
      coefficient * predictor[[row]] - log(sum(exp(
        coefficient * predictor[risk_set]
      )))
    }, numeric(1)))
  }
  stats::optimize(
    function(coefficient) -partial_log_likelihood(coefficient),
    interval = c(-10, 10),
    tol = 1e-12
  )$minimum
}

manual_sb_efron_cox_slope <- function(time, event, predictor, horizon) {
  administrative_time <- pmin(time, horizon)
  event_by_horizon <- as.integer(event == 1L & time <= horizon)
  event_times <- sort(unique(administrative_time[event_by_horizon == 1L]))
  partial_log_likelihood <- function(coefficient) {
    total <- 0
    for (event_time in event_times) {
      event_rows <- which(
        administrative_time == event_time & event_by_horizon == 1L
      )
      risk_rows <- which(administrative_time >= event_time)
      n_tied_events <- length(event_rows)
      risk_sum <- sum(exp(coefficient * predictor[risk_rows]))
      tied_event_sum <- sum(exp(coefficient * predictor[event_rows]))
      denominators <- vapply(
        0:(n_tied_events - 1L),
        function(tie_index) {
          risk_sum - tie_index / n_tied_events * tied_event_sum
        },
        numeric(1)
      )
      total <- total + coefficient * sum(predictor[event_rows]) -
        sum(log(denominators))
    }
    total
  }
  stats::optimize(
    function(coefficient) -partial_log_likelihood(coefficient),
    interval = c(-10, 10),
    tol = 1e-12
  )$minimum
}

call_sb_survival_weak_calibration <- function(
    data,
    horizon = 4,
    fold = data$fold,
    prediction_status = rep("success", nrow(data))) {
  sb_survival_weak_calibration_metric(
    patient_id = data$patient_id,
    time = data$time,
    event = data$event,
    survival_probability = data$survival_probability,
    prediction_status = prediction_status,
    prediction_horizon = horizon,
    time_unit = "month",
    scope = "test_fixed_horizon",
    fold = fold
  )
}

test_that("survival weak calibration matches independent scalar oracles", {
  data <- make_sb_survival_weak_calibration_fixture()
  horizon <- 4

  observed <- call_sb_survival_weak_calibration(data, horizon)
  transformed_prediction <- log(-log(data$survival_probability))
  expected_event_probability <- manual_sb_event_probability_at(
    data$time,
    data$event,
    horizon
  )
  expected_mean_prediction <- mean(-expm1(log(
    data$survival_probability
  )))
  expected_slope <- manual_sb_no_tie_cox_slope(
    data$time,
    data$event,
    transformed_prediction,
    horizon
  )

  expect_identical(observed$status, "success")
  expect_identical(observed$input_status, "success")
  expect_identical(observed$mean_calibration$status, "success")
  expect_identical(observed$slope_fit$status, "success")
  expect_equal(
    observed$mean_calibration$observed_event_probability,
    expected_event_probability,
    tolerance = 1e-12
  )
  expect_equal(
    observed$mean_calibration$mean_predicted_event_probability,
    expected_mean_prediction,
    tolerance = 1e-15
  )
  expect_equal(
    observed$mean_calibration_ratio,
    expected_event_probability / expected_mean_prediction,
    tolerance = 1e-12
  )
  expect_equal(observed$calibration_slope, expected_slope, tolerance = 5e-6)
  expect_identical(observed$mean_calibration$ratio_direction, "observed_over_expected")
  expect_identical(observed$mean_calibration_target, 1)
  expect_identical(observed$slope_target, 1)
  expect_identical(observed$slope_fit$ties, "efron")
  expect_identical(
    observed$slope_fit$administrative_censoring_at_horizon,
    TRUE
  )
})

test_that("survival weak calibration uses the Efron tied-event likelihood", {
  transformed_prediction <- c(-1.1, 0.2, 0.8, -0.4, 1.5, 0.5)
  data <- data.frame(
    patient_id = paste0("tie_", seq_along(transformed_prediction)),
    fold = rep(1L, length(transformed_prediction)),
    time = c(1, 1, 2, 3, 4, 5),
    event = c(1L, 1L, 1L, 0L, 1L, 0L),
    survival_probability = exp(-exp(transformed_prediction)),
    stringsAsFactors = FALSE
  )

  observed <- call_sb_survival_weak_calibration(data, horizon = 5)
  expected_efron <- manual_sb_efron_cox_slope(
    data$time,
    data$event,
    transformed_prediction,
    5
  )
  breslow <- unname(stats::coef(survival::coxph(
    survival::Surv(time, event) ~ transformed_prediction,
    data = data,
    ties = "breslow"
  ))[[1L]])

  expect_identical(observed$slope_fit$status, "success")
  expect_identical(observed$slope_fit$ties, "efron")
  expect_equal(observed$calibration_slope, expected_efron, tolerance = 5e-6)
  expect_gt(abs(observed$calibration_slope - breslow), 0.05)
})

test_that("survival weak calibration preserves row evidence and row invariance", {
  data <- make_sb_survival_weak_calibration_fixture()
  original <- call_sb_survival_weak_calibration(data)
  order <- c(8, 2, 6, 1, 7, 3, 5, 4)
  shuffled <- call_sb_survival_weak_calibration(data[order, ])

  expect_equal(
    shuffled$mean_calibration_ratio,
    original$mean_calibration_ratio,
    tolerance = 1e-12
  )
  expect_equal(
    shuffled$calibration_slope,
    original$calibration_slope,
    tolerance = 1e-12
  )
  expect_identical(original$patient_predictions$row_instance, seq_len(nrow(data)))
  expect_identical(original$patient_predictions$patient_id, data$patient_id)
  expect_equal(
    original$patient_predictions$predicted_event_probability,
    -expm1(log(data$survival_probability)),
    tolerance = 0
  )
  expect_equal(
    original$patient_predictions$transformed_prediction,
    log(-log(data$survival_probability)),
    tolerance = 0
  )
  expect_identical(
    original$patient_predictions$administrative_time,
    pmin(data$time, 4)
  )
  expect_identical(
    original$patient_predictions$event_by_horizon,
    as.integer(data$event == 1L & data$time <= 4)
  )

  sampled <- data[c(1, 1, 3, 4, 6, 6, 8, 2), ]
  repeated <- call_sb_survival_weak_calibration(sampled)
  expect_identical(repeated$patient_predictions$row_instance, seq_len(nrow(sampled)))
  expect_identical(repeated$patient_predictions$patient_id, sampled$patient_id)
})

test_that("survival weak calibration ignores outcomes strictly after horizon", {
  data <- make_sb_survival_weak_calibration_fixture()
  changed <- data
  later <- changed$time > 4
  changed$time[later] <- changed$time[later] + 100
  changed$event[later] <- 1L - changed$event[later]

  original <- call_sb_survival_weak_calibration(data)
  modified <- call_sb_survival_weak_calibration(changed)

  expect_equal(
    modified$mean_calibration_ratio,
    original$mean_calibration_ratio,
    tolerance = 0
  )
  expect_equal(
    modified$calibration_slope,
    original$calibration_slope,
    tolerance = 0
  )
})

test_that("survival weak calibration keeps exact-horizon event semantics", {
  data <- data.frame(
    patient_id = paste0("horizon_", 1:7),
    fold = rep(1L, 7L),
    time = c(1, 2, 4, 4, 5, 6, 7),
    event = c(1L, 0L, 1L, 0L, 1L, 0L, 1L),
    survival_probability = exp(-exp(c(-0.5, 0.7, -0.1, 1.1, 0.4, -0.8, 0.2))),
    stringsAsFactors = FALSE
  )

  observed <- call_sb_survival_weak_calibration(data, horizon = 4)

  expect_identical(observed$patient_predictions$event_by_horizon, c(
    1L, 0L, 1L, 0L, 0L, 0L, 0L
  ))
  expect_identical(observed$patient_predictions$administrative_time, c(
    1, 2, 4, 4, 4, 4, 4
  ))
  expect_equal(
    observed$mean_calibration$observed_event_probability,
    manual_sb_event_probability_at(data$time, data$event, 4),
    tolerance = 0
  )
})

test_that("mean calibration survives boundary and constant predictions", {
  data <- make_sb_survival_weak_calibration_fixture()
  boundary <- data
  boundary$survival_probability[1:2] <- c(0, 1)
  boundary_result <- call_sb_survival_weak_calibration(boundary)

  expect_identical(boundary_result$mean_calibration$status, "success")
  expect_true(is.finite(boundary_result$mean_calibration_ratio))
  expect_identical(boundary_result$slope_fit$status, "not_evaluable")
  expect_identical(
    boundary_result$slope_fit$error_class,
    "stormbreaker_survival_weak_calibration_boundary_probability"
  )
  expect_identical(is.na(boundary_result$calibration_slope), TRUE)
  expect_identical(boundary_result$status, "incomplete")
  expect_identical(boundary_result$probability_clipping, "none")

  constant <- data
  constant$survival_probability <- rep(0.7, nrow(constant))
  constant_result <- call_sb_survival_weak_calibration(constant)

  expect_identical(constant_result$mean_calibration$status, "success")
  expect_true(is.finite(constant_result$mean_calibration_ratio))
  expect_identical(constant_result$slope_fit$status, "not_evaluable")
  expect_identical(
    constant_result$slope_fit$error_class,
    "stormbreaker_survival_weak_calibration_constant_prediction"
  )
})

test_that("events after horizon cannot rescue the calibration slope", {
  data <- data.frame(
    patient_id = paste0("later_", 1:6),
    fold = rep(1L, 6L),
    time = c(5, 6, 7, 8, 9, 10),
    event = c(1L, 0L, 1L, 0L, 1L, 0L),
    survival_probability = exp(-exp(c(-1, -0.5, 0, 0.5, 1, 1.5))),
    stringsAsFactors = FALSE
  )

  observed <- call_sb_survival_weak_calibration(data, horizon = 4)

  expect_identical(observed$mean_calibration$status, "success")
  expect_identical(observed$mean_calibration$observed_event_probability, 0)
  expect_identical(observed$mean_calibration_ratio, 0)
  expect_identical(observed$slope_fit$status, "not_evaluable")
  expect_identical(
    observed$slope_fit$error_class,
    "stormbreaker_survival_weak_calibration_no_event"
  )
  expect_identical(is.na(observed$calibration_slope), TRUE)
})

test_that("mean horizon support cannot erase an independently estimable slope", {
  transformed_prediction <- c(-0.8, 0.4, 1.1, -0.2, 0.7, -1.0, 0.1, 1.4)
  data <- data.frame(
    patient_id = paste0("short_", seq_along(transformed_prediction)),
    fold = rep(1L, length(transformed_prediction)),
    time = c(0.5, 0.9, 1.4, 2.0, 2.6, 3.1, 3.7, 4.4),
    event = c(1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L),
    survival_probability = exp(-exp(transformed_prediction)),
    stringsAsFactors = FALSE
  )

  observed <- call_sb_survival_weak_calibration(data, horizon = 10)
  expected_slope <- manual_sb_no_tie_cox_slope(
    data$time,
    data$event,
    transformed_prediction,
    10
  )

  expect_identical(observed$mean_calibration$status, "not_evaluable")
  expect_identical(
    observed$mean_calibration$error_class,
    "stormbreaker_survival_weak_calibration_horizon_unsupported"
  )
  expect_identical(is.na(observed$mean_calibration_ratio), TRUE)
  expect_identical(observed$slope_fit$status, "success")
  expect_equal(observed$calibration_slope, expected_slope, tolerance = 5e-6)
  expect_identical(observed$status, "incomplete")

  absorbing <- data.frame(
    patient_id = paste0("absorbing_", 1:4),
    fold = rep(1L, 4L),
    time = c(1, 2, 3, 4),
    event = rep(1L, 4L),
    survival_probability = exp(-exp(c(-1, -0.2, 0.4, 1.1))),
    stringsAsFactors = FALSE
  )
  absorbing_result <- call_sb_survival_weak_calibration(
    absorbing,
    horizon = 10
  )
  expect_identical(absorbing_result$mean_calibration$status, "success")
  expect_identical(
    absorbing_result$mean_calibration$horizon_support,
    "survival_zero_before_horizon"
  )
  expect_identical(
    absorbing_result$mean_calibration$observed_event_probability,
    1
  )
})

test_that("survival weak calibration fit conditions stay component-local", {
  data <- make_sb_survival_weak_calibration_fixture()
  baseline <- call_sb_survival_weak_calibration(data)
  original_fit <- sb_survival_weak_calibration_fit_slope

  testthat::local_mocked_bindings(
    sb_survival_weak_calibration_fit_slope = function(fit_data) {
      warning("secondary slope warning")
      original_fit(fit_data)
    }
  )
  warned <- call_sb_survival_weak_calibration(data)

  expect_identical(warned$mean_calibration, baseline$mean_calibration)
  expect_identical(warned$slope_fit$status, "failed")
  expect_identical(
    warned$slope_fit$error_class,
    "stormbreaker_survival_weak_calibration_cox_warning"
  )
  expect_match(warned$slope_fit$error_message, "secondary slope warning", fixed = TRUE)

  testthat::local_mocked_bindings(
    sb_survival_weak_calibration_fit_slope = function(fit_data) {
      condition <- simpleError("secondary slope error")
      class(condition) <- c("test_weak_slope_error", class(condition))
      stop(condition)
    }
  )
  failed <- call_sb_survival_weak_calibration(data)

  expect_identical(failed$mean_calibration, baseline$mean_calibration)
  expect_identical(failed$slope_fit$status, "failed")
  expect_identical(failed$slope_fit$error_class, "test_weak_slope_error")
  expect_identical(failed$slope_fit$error_message, "secondary slope error")
})

test_that("survival weak calibration exposes input and unexpected failures", {
  data <- make_sb_survival_weak_calibration_fixture()
  failed_status <- rep("success", nrow(data))
  failed_status[[2L]] <- "failed"
  failed <- call_sb_survival_weak_calibration(
    data,
    prediction_status = failed_status
  )

  expect_identical(failed$status, "incomplete")
  expect_identical(failed$input_status, "prediction_failed")
  expect_identical(failed$mean_calibration$status, "not_evaluated")
  expect_identical(failed$slope_fit$status, "not_fitted")

  condition <- simpleError("outer weak calibration error")
  class(condition) <- c("test_weak_calibration_error", class(condition))
  outer <- sb_survival_weak_calibration_error_metric(
    error = condition,
    prediction_horizon = 4,
    time_unit = "month",
    scope = "test_fixed_horizon",
    n = nrow(data),
    fold = data$fold
  )

  expect_identical(outer$status, "metric_failed")
  expect_identical(outer$input_status, "unexpected_error")
  expect_identical(outer$error_class, "test_weak_calibration_error")
  expect_identical(outer$error_message, "outer weak calibration error")
  expect_identical(outer$n_folds, 4L)
})
})

# ---- metrics-survival-discrimination ----
local({
# ---- metrics-concordance ----

test_that("Harrell C matches a three-patient hand calculation", {
  predictions <- data.frame(
    fold = rep(1L, 3L),
    time = c(2, 5, 6),
    event = c(1L, 1L, 0L),
    risk = c(0.8, 0.3, 0.5),
    status = rep("success", 3L)
  )

  metric <- sb_harrell_c_metric(predictions)

  required_metric_fields <- c(
    "metric", "scope", "aggregation", "time_weight",
    "risk_direction", "status", "estimate", "n_folds",
    "n_evaluable_folds", "concordant", "discordant",
    "tied_risk", "tied_time", "tied_both", "comparable_pairs",
    "per_fold"
  )
  expect_length(setdiff(required_metric_fields, names(metric)), 0L)
  expect_identical(metric$metric, "harrell_c")
  expect_identical(metric$scope, "within_outer_assessment_fold")
  expect_identical(metric$aggregation, "summed_comparable_pairs")
  expect_identical(metric$time_weight, "n")
  expect_identical(
    metric$risk_direction,
    "higher_is_higher_hazard"
  )
  expect_identical(metric$status, "success")
  expect_equal(metric$estimate, 2 / 3, tolerance = 1e-15)
  expect_identical(metric$n_folds, 1L)
  expect_identical(metric$n_evaluable_folds, 1L)
  expect_equal(metric$concordant, 2, tolerance = 0)
  expect_equal(metric$discordant, 1, tolerance = 0)
  expect_equal(metric$tied_risk, 0, tolerance = 0)
  expect_equal(metric$comparable_pairs, 3, tolerance = 0)
  required_fold_fields <- c(
    "fold", "status", "n", "n_events", "estimate",
    "concordant", "discordant", "tied_risk", "tied_time",
    "tied_both", "comparable_pairs", "error_class",
    "error_message"
  )
  expect_length(setdiff(required_fold_fields, names(metric$per_fold)), 0L)
})

test_that("Harrell C stratum exposes survival concordance counts", {
  audit <- sb_harrell_c_stratum(
    time = c(2, 5, 6),
    event = c(1L, 1L, 0L),
    risk = c(0.8, 0.3, 0.5)
  )

  expect_identical(audit$status, "success")
  expect_equal(audit$estimate, 2 / 3, tolerance = 1e-15)
  expect_equal(audit$concordant, 2, tolerance = 0)
  expect_equal(audit$discordant, 1, tolerance = 0)
  expect_equal(audit$tied_risk, 0, tolerance = 0)
  expect_equal(audit$comparable_pairs, 3, tolerance = 0)
})

test_that("Harrell C stratum keeps tied-risk half credit", {
  audit <- sb_harrell_c_stratum(
    time = c(1, 2),
    event = c(1L, 0L),
    risk = c(0.5, 0.5)
  )

  expect_identical(audit$status, "success")
  expect_equal(audit$estimate, 0.5, tolerance = 0)
  expect_equal(audit$tied_risk, 1, tolerance = 0)
})

test_that("Harrell C stratum and fold wrapper distinguish all tie counts", {
  count_fields <- c(
    "concordant", "discordant", "tied_risk", "tied_time",
    "tied_both", "comparable_pairs"
  )
  cases <- list(
    list(risk = c(0.9, 0.8, 0.2), expected = c(2, 0, 0, 1, 0, 2)),
    list(risk = c(0.5, 0.5, 0.2), expected = c(2, 0, 0, 0, 1, 2))
  )

  for (case in cases) {
    audit <- sb_harrell_c_stratum(
      time = c(1, 1, 2),
      event = c(1L, 1L, 0L),
      risk = case$risk
    )
    predictions <- data.frame(
      fold = rep(1L, 3L),
      time = c(1, 1, 2),
      event = c(1L, 1L, 0L),
      risk = case$risk,
      status = rep("success", 3L)
    )
    metric <- sb_harrell_c_metric(predictions)

    expect_equal(
      unname(unlist(audit[count_fields], use.names = FALSE)),
      case$expected,
      tolerance = 0
    )
    expect_equal(
      unname(unlist(metric$per_fold[1L, count_fields], use.names = FALSE)),
      case$expected,
      tolerance = 0
    )
  }
})

test_that("Harrell C compares only patients predicted by the same fold model", {
  predictions <- data.frame(
    fold = rep(1:2, each = 3L),
    time = c(2, 5, 6, 1, 4, 7),
    event = c(1L, 1L, 0L, 1L, 0L, 1L),
    risk = c(0.8, 0.3, 0.5, 0.9, 0.4, 0.2),
    status = rep("success", 6L)
  )

  observed <- sb_harrell_c_metric(predictions)
  shifted <- predictions
  shifted$risk[shifted$fold == 2L] <-
    shifted$risk[shifted$fold == 2L] + 1000
  shifted_observed <- sb_harrell_c_metric(shifted)
  reversed <- predictions
  reversed$risk <- -reversed$risk
  reversed_observed <- sb_harrell_c_metric(reversed)
  official <- survival::concordance(
    survival::Surv(time, event) ~ risk + survival::strata(fold),
    data = predictions,
    reverse = TRUE,
    timewt = "n"
  )

  expect_equal(observed$estimate, 4 / 5, tolerance = 1e-15)
  expect_equal(shifted_observed$estimate, observed$estimate, tolerance = 0)
  expect_equal(reversed_observed$estimate, 1 - observed$estimate,
               tolerance = 1e-15)
  expect_equal(observed$estimate, unname(official$concordance),
               tolerance = 1e-15)
  expect_equal(observed$per_fold$estimate, c(2 / 3, 1),
               tolerance = 1e-15)
  expect_equal(observed$per_fold$comparable_pairs, c(3, 2),
               tolerance = 0)
  expect_equal(
    observed$estimate,
    sum(observed$per_fold$concordant) /
      sum(observed$per_fold$comparable_pairs),
    tolerance = 1e-15
  )
})

test_that("risk-score ties receive half credit", {
  predictions <- data.frame(
    fold = c(1L, 1L),
    time = c(1, 2),
    event = c(1L, 0L),
    risk = c(0.5, 0.5),
    status = c("success", "success")
  )

  metric <- sb_harrell_c_metric(predictions)

  expect_identical(metric$status, "success")
  expect_equal(metric$estimate, 0.5, tolerance = 0)
  expect_equal(metric$tied_risk, 1, tolerance = 0)
  expect_equal(metric$comparable_pairs, 1, tolerance = 0)
})

test_that("a fold without comparable pairs is visibly not evaluable", {
  predictions <- data.frame(
    fold = rep(1L, 3L),
    time = c(5, 1, 2),
    event = c(1L, 0L, 0L),
    risk = c(0.8, 0.3, 0.5),
    status = rep("success", 3L)
  )

  metric <- sb_harrell_c_metric(predictions)

  expect_identical(metric$status, "incomplete")
  expect_identical(metric$n_evaluable_folds, 0L)
  expect_identical(is.na(metric$estimate), TRUE)
  expect_identical(metric$per_fold$status, "not_evaluable")
  expect_match(metric$per_fold$error_message, "No comparable patient pairs")
})

test_that("failed OOF predictions cannot yield a complete concordance", {
  predictions <- data.frame(
    fold = rep(1:2, each = 3L),
    time = c(2, 5, 6, 1, 4, 7),
    event = c(1L, 1L, 0L, 1L, 0L, 1L),
    risk = c(0.8, 0.3, 0.5, NA, NA, NA),
    status = c(rep("success", 3L), rep("failed", 3L))
  )

  metric <- sb_harrell_c_metric(predictions)

  expect_identical(metric$status, "incomplete")
  expect_identical(is.na(metric$estimate), TRUE)
  expect_identical(metric$n_folds, 2L)
  expect_identical(metric$n_evaluable_folds, 1L)
  expect_identical(
    metric$per_fold$status,
    c("success", "prediction_failed")
  )
  expect_match(
    metric$per_fold$error_message[2],
    "failed or non-finite"
  )
})

# ---- metrics-survival-time-auc ----

make_sb_time_auc_curve <- function(time, event) {
  sb_fit_censoring_km(time, event)
}

make_sb_time_auc_stratum <- function(time, event, marker, horizon, curve,
                                     prediction_ok = NULL, fold = 1L) {
  if (is.null(prediction_ok)) {
    prediction_ok <- rep(TRUE, length(time))
  }
  sb_time_auc_stratum(
    time = time,
    event = event,
    survival_probability = 1 - marker,
    prediction_ok = prediction_ok,
    censoring_curve = curve,
    prediction_horizon = horizon,
    fold = fold
  )
}

test_that("time AUC matches a literal unequal-IPCW pair oracle", {
  horizon <- 3
  time <- c(1.5, 2.5, 4, 5)
  event <- c(1L, 1L, 0L, 1L)
  marker <- c(0.9, 0.5, 0.5, 0.2)
  curve <- make_sb_time_auc_curve(
    c(1, 2, 3, 4, 5),
    c(0L, 0L, 1L, 1L, 1L)
  )
  cases <- event == 1L & time <= horizon
  controls <- time > horizon
  case_weights <- 1 / vapply(
    time[cases],
    function(value) sb_censoring_survival_at(
      curve,
      value,
      left_limit = TRUE
    ),
    numeric(1)
  )
  pair_credit <- outer(
    marker[cases],
    marker[controls],
    function(case, control) {
      (case > control) + 0.5 * (case == control)
    }
  )
  expected_concordant <- sum(case_weights * rowSums(
    outer(marker[cases], marker[controls], `>`)
  ))
  expected_tied <- sum(case_weights * rowSums(
    outer(marker[cases], marker[controls], `==`)
  ))
  expected_numerator <- sum(case_weights * rowSums(pair_credit))
  expected_denominator <- sum(case_weights) * sum(controls)

  observed <- make_sb_time_auc_stratum(
    time,
    event,
    marker,
    horizon,
    curve
  )

  expect_identical(observed$status, "success")
  expect_identical(observed$n_cases, 2L)
  expect_identical(observed$n_controls, 2L)
  expect_identical(observed$n_excluded_censored, 0L)
  expect_equal(observed$case_weight_sum, sum(case_weights), tolerance = 1e-15)
  expect_equal(
    observed$weighted_concordant_pairs,
    expected_concordant,
    tolerance = 1e-15
  )
  expect_equal(
    observed$weighted_tied_pairs,
    expected_tied,
    tolerance = 1e-15
  )
  expect_equal(
    observed$weighted_pair_credit,
    expected_numerator,
    tolerance = 1e-15
  )
  expect_equal(
    observed$pair_denominator,
    expected_denominator,
    tolerance = 1e-15
  )
  expect_equal(
    observed$estimate,
    expected_numerator / expected_denominator,
    tolerance = 1e-15
  )
})

test_that("time AUC fixes the discrete horizon case and control boundary", {
  horizon <- 2
  time <- c(2, 2, 3, 1.5)
  event <- c(1L, 0L, 1L, 0L)
  marker <- c(0.8, 0.9, 0.2, 0.1)
  curve <- make_sb_time_auc_curve(1:4, rep(1L, 4L))

  observed <- make_sb_time_auc_stratum(
    time,
    event,
    marker,
    horizon,
    curve
  )

  expect_identical(observed$status, "success")
  expect_identical(observed$n_cases, 1L)
  expect_identical(observed$n_controls, 1L)
  expect_identical(observed$n_excluded_censored, 2L)
  expect_equal(observed$estimate, 1, tolerance = 0)

  before_only <- make_sb_time_auc_stratum(
    time,
    event,
    marker,
    horizon - .Machine$double.eps^0.5,
    curve
  )
  expect_identical(before_only$status, "not_evaluable")
  expect_identical(before_only$n_cases, 0L)
})

test_that("time AUC never requires the cancelled control G at the horizon", {
  horizon <- 2
  curve <- make_sb_time_auc_curve(
    c(1, 2, 2),
    c(1L, 1L, 0L)
  )
  expect_equal(
    sb_censoring_survival_at(curve, horizon, left_limit = TRUE),
    1,
    tolerance = 0
  )
  expect_equal(
    sb_censoring_survival_at(curve, horizon, left_limit = FALSE),
    0,
    tolerance = 0
  )

  observed <- make_sb_time_auc_stratum(
    time = c(2, 3),
    event = c(1L, 1L),
    marker = c(0.8, 0.2),
    horizon = horizon,
    curve = curve
  )

  expect_identical(observed$status, "success")
  expect_equal(observed$case_weight_sum, 1, tolerance = 0)
  expect_equal(observed$pair_denominator, 1, tolerance = 0)
  expect_equal(observed$estimate, 1, tolerance = 0)
})

test_that("time AUC exposes target-specific unevaluable and failure states", {
  curve <- make_sb_time_auc_curve(1:4, rep(1L, 4L))
  no_case <- make_sb_time_auc_stratum(
    c(3, 4), c(1L, 0L), c(0.8, 0.2), 2, curve
  )
  no_control <- make_sb_time_auc_stratum(
    c(1, 2), c(1L, 0L), c(0.8, 0.2), 2, curve
  )
  bad_prediction <- make_sb_time_auc_stratum(
    c(1, 3), c(1L, 1L), c(0.8, 0.2), 2, curve,
    prediction_ok = c(TRUE, FALSE)
  )
  zero_curve <- make_sb_time_auc_curve(c(1, 1), c(1L, 0L))
  zero_case_support <- make_sb_time_auc_stratum(
    c(2, 3), c(1L, 1L), c(0.8, 0.2), 2, zero_curve
  )
  empty_curve <- make_sb_time_auc_stratum(
    c(1, 3), c(1L, 1L), c(0.8, 0.2), 2,
    sb_empty_censoring_curve()
  )

  expect_identical(no_case$status, "not_evaluable")
  expect_identical(no_case$error_class, "stormbreaker_no_time_auc_case")
  expect_identical(no_control$status, "not_evaluable")
  expect_identical(no_control$error_class, "stormbreaker_no_time_auc_control")
  expect_identical(bad_prediction$status, "prediction_failed")
  expect_identical(
    bad_prediction$error_class,
    "stormbreaker_invalid_survival_prediction"
  )
  expect_identical(zero_case_support$status, "censoring_failed")
  expect_identical(
    zero_case_support$error_class,
    "stormbreaker_zero_censoring_survival"
  )
  expect_identical(empty_curve$status, "censoring_failed")
  expect_identical(
    empty_curve$error_class,
    "stormbreaker_missing_censoring_curve"
  )
  expect_true(all(is.na(c(
    no_case$estimate,
    no_control$estimate,
    bad_prediction$estimate,
    zero_case_support$estimate,
    empty_curve$estimate
  ))))
})

test_that("time AUC is order invariant and obeys strict ranking transforms", {
  time <- c(1, 1.5, 3, 4)
  event <- rep(1L, 4L)
  marker <- c(0.8, 0.6, 0.7, 0.5)
  curve <- make_sb_time_auc_curve(1:4, rep(1L, 4L))
  observed <- make_sb_time_auc_stratum(time, event, marker, 2, curve)
  shuffled <- make_sb_time_auc_stratum(
    time[c(4, 2, 1, 3)],
    event[c(4, 2, 1, 3)],
    marker[c(4, 2, 1, 3)],
    2,
    curve
  )
  transformed <- make_sb_time_auc_stratum(
    time,
    event,
    marker^2,
    2,
    curve
  )
  reversed <- make_sb_time_auc_stratum(
    time,
    event,
    1 - marker,
    2,
    curve
  )

  expect_equal(observed$estimate, 0.75, tolerance = 0)
  expect_equal(shuffled$estimate, observed$estimate, tolerance = 0)
  expect_equal(transformed$estimate, observed$estimate, tolerance = 0)
  expect_equal(reversed$estimate, 1 - observed$estimate, tolerance = 0)
})

test_that("CV time AUC sums only within-fold weighted pair quantities", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:5),
    fold = c(1L, 1L, 2L, 2L, 2L),
    time = c(1, 3, 1, 3, 4),
    event = rep(1L, 5L),
    survival_probability = 1 - c(0.2, 0.1, 0.9, 0.8, 0.7),
    survival_probability_status = rep("success", 5L),
    status = rep("success", 5L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(
    list(fold = 1L, time = 1:4, event = rep(1L, 4L)),
    list(fold = 2L, time = 1:4, event = rep(1L, 4L))
  )
  curves <- lapply(analysis_inputs, function(input) {
    make_sb_time_auc_curve(input$time, input$event)
  })
  names(curves) <- c("1", "2")

  observed <- sb_survival_time_auc_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 2,
    censoring_curves = curves
  )

  expect_identical(observed$status, "success")
  expect_equal(observed$estimate, 1, tolerance = 0)
  expect_equal(observed$weighted_pair_credit, 3, tolerance = 0)
  expect_equal(observed$pair_denominator, 3, tolerance = 0)
  expect_identical(observed$n_folds, 2L)
  expect_identical(observed$n_evaluable_folds, 2L)
  expect_equal(observed$per_fold$pair_denominator, c(1, 2), tolerance = 0)

  marker <- 1 - predictions$survival_probability
  pooled_credit <- outer(
    marker[predictions$time <= 2],
    marker[predictions$time > 2],
    function(case, control) (case > control) + 0.5 * (case == control)
  )
  expect_equal(mean(pooled_credit), 4 / 6, tolerance = 0)
  expect_gt(abs(observed$estimate - mean(pooled_credit)), 0)
  expect_equal(
    observed$pair_denominator,
    sum(observed$per_fold$pair_denominator),
    tolerance = 0
  )
})

test_that("CV time AUC requires exact analysis and censoring fold mapping", {
  predictions <- data.frame(
    patient_id = c("case", "control"),
    fold = rep(1L, 2L),
    time = c(1, 3),
    event = rep(1L, 2L),
    survival_probability = c(0.2, 0.8),
    survival_probability_status = rep("success", 2L),
    status = rep("success", 2L),
    stringsAsFactors = FALSE
  )
  input <- list(fold = 1L, time = 1:3, event = rep(1L, 3L))
  curve <- make_sb_time_auc_curve(input$time, input$event)
  duplicated_input <- sb_survival_time_auc_metric(
    predictions,
    list(input, input),
    2,
    list(`1` = curve)
  )
  missing_curve <- sb_survival_time_auc_metric(
    predictions,
    list(input),
    2,
    list(`2` = curve)
  )
  duplicate_curve <- sb_survival_time_auc_metric(
    predictions,
    list(input),
    2,
    structure(list(curve, curve), names = c("1", "1"))
  )

  expect_identical(duplicated_input$status, "incomplete")
  expect_identical(
    duplicated_input$per_fold$error_class,
    "stormbreaker_bad_time_auc_analysis_input"
  )
  expect_identical(missing_curve$status, "incomplete")
  expect_identical(
    missing_curve$per_fold$error_class,
    "stormbreaker_missing_censoring_curve"
  )
  expect_identical(duplicate_curve$status, "incomplete")
  expect_identical(
    duplicate_curve$per_fold$error_class,
    "stormbreaker_bad_time_auc_censoring_curve"
  )
})

test_that("an incomplete time AUC retains successful fold evidence", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:4),
    fold = c(1L, 1L, 2L, 2L),
    time = c(1, 3, 1, 1.5),
    event = rep(1L, 4L),
    survival_probability = c(0.2, 0.8, 0.2, 0.8),
    survival_probability_status = rep("success", 4L),
    status = rep("success", 4L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(
    list(fold = 1L, time = 1:3, event = rep(1L, 3L)),
    list(fold = 2L, time = 1:3, event = rep(1L, 3L))
  )
  curves <- lapply(analysis_inputs, function(input) {
    make_sb_time_auc_curve(input$time, input$event)
  })
  names(curves) <- c("1", "2")

  observed <- sb_survival_time_auc_metric(
    predictions,
    analysis_inputs,
    2,
    curves
  )

  expect_identical(observed$status, "incomplete")
  expect_true(is.na(observed$estimate))
  expect_identical(observed$n_evaluable_folds, 1L)
  expect_equal(observed$weighted_pair_credit, 1, tolerance = 0)
  expect_equal(observed$pair_denominator, 1, tolerance = 0)
  expect_identical(observed$per_fold$status, c("success", "not_evaluable"))
})

# ---- metrics-survival-uno-c ----

uno_test_curve <- function(time = c(1, 2, 3, 4, 5),
                           survival_after = c(1, 0.8, 0.8, 0.4, 0.4)) {
  data.frame(
    time = time,
    n_risk = rep.int(1L, length(time)),
    n_event = rep.int(0L, length(time)),
    n_censoring = rep.int(0L, length(time)),
    n_risk_for_censoring = rep.int(1L, length(time)),
    survival_before = c(1, head(survival_after, -1L)),
    survival_after = survival_after,
    stringsAsFactors = FALSE
  )
}

literal_uno_pairs <- function(time, event, risk, curve, horizon) {
  numerator <- denominator <- concordant <- tied_risk <- 0
  comparable_pairs <- tied_time_pairs <- 0L
  eligible <- which(event == 1L & time < horizon)
  for (index in eligible) {
    g_left <- if (any(curve$time < time[[index]])) {
      utils::tail(curve$survival_after[curve$time < time[[index]]], 1L)
    } else {
      1
    }
    weight <- 1 / g_left^2
    comparators <- which(
      time > time[[index]] |
        (time == time[[index]] & event == 0L)
    )
    comparable_pairs <- comparable_pairs + length(comparators)
    tied_time_pairs <- tied_time_pairs + sum(
      time[comparators] == time[[index]]
    )
    concordant <- concordant + weight * sum(
      risk[[index]] > risk[comparators]
    )
    tied_risk <- tied_risk + weight * sum(
      risk[[index]] == risk[comparators]
    )
    denominator <- denominator + weight * length(comparators)
  }
  numerator <- concordant + 0.5 * tied_risk
  list(
    estimate = numerator / denominator,
    numerator = numerator,
    denominator = denominator,
    concordant = concordant,
    tied_risk = tied_risk,
    comparable_pairs = as.integer(comparable_pairs),
    tied_time_pairs = as.integer(tied_time_pairs)
  )
}

test_that("Uno C uses squared left-limit IPCW and literal pair credit", {
  curve <- uno_test_curve()
  time <- c(1.5, 2.5, 3.5, 3.5, 5.5)
  event <- c(1L, 1L, 1L, 0L, 0L)
  risk <- c(3, 1, 2, 2, 0)
  expected <- literal_uno_pairs(time, event, risk, curve, 5)

  result <- sb_uno_c_stratum(
    time,
    event,
    risk,
    rep(TRUE, length(time)),
    curve,
    prediction_horizon = 5,
    fold = 1L
  )

  expect_identical(result$status, "success")
  expect_equal(result$weighted_concordant_pairs, 7.125)
  expect_equal(result$weighted_tied_risk_pairs, 1.5625)
  expect_equal(result$weighted_pair_credit, 7.90625)
  expect_equal(result$weighted_pair_denominator, 11.8125)
  expect_equal(result$estimate, expected$estimate)
  expect_identical(result$n_comparable_pairs, 9L)
  expect_identical(result$n_tied_time_comparable_pairs, 1L)
  expect_equal(
    result$estimate,
    result$weighted_pair_credit / result$weighted_pair_denominator
  )
  expect_false(isTRUE(all.equal(result$estimate, 7.25 / 10.25)))
})

test_that("Uno C locks strict horizon and observed-time tie conventions", {
  curve <- uno_test_curve(c(0.5, 1.5, 2.5), c(1, 1, 1))
  boundary <- sb_uno_c_stratum(
    time = c(1.999999, 2, 3),
    event = c(1L, 1L, 0L),
    risk = c(2, 100, 0),
    prediction_ok = rep(TRUE, 3L),
    censoring_curve = curve,
    prediction_horizon = 2,
    fold = 1L
  )
  expect_identical(boundary$status, "success")
  expect_identical(boundary$n_eligible_events, 1L)
  expect_identical(boundary$n_events_at_horizon_excluded, 1L)
  expect_equal(boundary$estimate, 0.5)

  tied_time <- sb_uno_c_stratum(
    time = c(1, 1, 1, 2),
    event = c(1L, 0L, 1L, 0L),
    risk = c(3, 1, 100, 0),
    prediction_ok = rep(TRUE, 4L),
    censoring_curve = curve,
    prediction_horizon = 2.5,
    fold = 1L
  )
  expect_identical(tied_time$status, "success")
  expect_identical(tied_time$n_comparable_pairs, 4L)
  expect_identical(tied_time$n_tied_time_comparable_pairs, 2L)
  expect_equal(tied_time$estimate, 1)
})

test_that("Uno C gives half credit only to exactly equal risks", {
  result <- sb_uno_c_stratum(
    time = c(1, 2, 3),
    event = c(1L, 0L, 0L),
    risk = c(1, 1, 1 + 1e-12),
    prediction_ok = rep(TRUE, 3L),
    censoring_curve = uno_test_curve(c(0.5, 2.5), c(1, 1)),
    prediction_horizon = 4,
    fold = 1L
  )
  expect_identical(result$status, "success")
  expect_equal(result$weighted_tied_risk_pairs, 1)
  expect_equal(result$estimate, 0.25)
})

test_that("Uno C reports unavailable inputs without fallback", {
  curve <- uno_test_curve(c(0.5, 1.5), c(1, 1))
  no_case <- sb_uno_c_stratum(
    c(1, 2), c(0L, 0L), c(2, 1), c(TRUE, TRUE), curve, 3, 1L
  )
  expect_identical(no_case$status, "not_evaluable")
  expect_identical(no_case$error_class, "stormbreaker_no_uno_c_case")
  expect_true(is.na(no_case$estimate))

  no_pair <- sb_uno_c_stratum(
    c(2, 1), c(1L, 0L), c(2, 1), c(TRUE, TRUE), curve, 3, 1L
  )
  expect_identical(no_pair$status, "not_evaluable")
  expect_identical(
    no_pair$error_class,
    "stormbreaker_no_uno_c_comparable_pair"
  )

  bad_risk <- sb_uno_c_stratum(
    c(1, 2), c(1L, 0L), c(NA_real_, 1), c(TRUE, TRUE), curve, 3, 1L
  )
  expect_identical(bad_risk$status, "prediction_failed")
  expect_identical(
    bad_risk$error_class,
    "stormbreaker_invalid_uno_c_risk"
  )

  bad_status <- sb_uno_c_stratum(
    c(1, 2), c(1L, 0L), c(2, 1), c(TRUE, FALSE), curve, 3, 1L
  )
  expect_identical(bad_status$status, "prediction_failed")

  zero_horizon_support <- sb_uno_c_stratum(
    c(0.5, 3),
    c(1L, 0L),
    c(2, 1),
    c(TRUE, TRUE),
    uno_test_curve(c(1, 2), c(0.5, 0)),
    3,
    1L
  )
  expect_identical(zero_horizon_support$status, "censoring_failed")
  expect_identical(
    zero_horizon_support$error_class,
    "stormbreaker_zero_uno_c_horizon_censoring_survival"
  )

  missing_curve <- sb_uno_c_stratum(
    c(1, 2),
    c(1L, 0L),
    c(2, 1),
    c(TRUE, TRUE),
    sb_empty_censoring_curve(),
    3,
    1L
  )
  expect_identical(missing_curve$status, "censoring_failed")
  expect_identical(
    missing_curve$error_class,
    "stormbreaker_missing_uno_c_censoring_curve"
  )
})

test_that("Uno C is invariant to row order and strict risk transformations", {
  curve <- uno_test_curve()
  time <- c(1.5, 2.5, 3.5, 5.5)
  event <- c(1L, 1L, 1L, 0L)
  risk <- c(3, 1, 2, 0)
  base <- sb_uno_c_stratum(
    time, event, risk, rep(TRUE, 4L), curve, 5, 1L
  )
  order <- c(4, 2, 1, 3)
  reordered <- sb_uno_c_stratum(
    time[order], event[order], risk[order], rep(TRUE, 4L), curve, 5, 1L
  )
  transformed <- sb_uno_c_stratum(
    time, event, exp(risk), rep(TRUE, 4L), curve, 5, 1L
  )
  reversed <- sb_uno_c_stratum(
    time, event, -risk, rep(TRUE, 4L), curve, 5, 1L
  )

  expect_equal(reordered$estimate, base$estimate)
  expect_equal(transformed$estimate, base$estimate)
  expect_equal(reversed$estimate, 1 - base$estimate)
})

test_that("Uno C aggregates weighted pair totals without cross-fold ranking", {
  curve <- uno_test_curve(c(0.25, 6), c(1, 1))
  predictions <- data.frame(
    patient_id = paste0("p", 1:6),
    fold = c(1L, 1L, 2L, 2L, 2L, 2L),
    time = c(0.5, 4, 2, 3, 4, 5),
    event = c(1L, 0L, 1L, 0L, 0L, 0L),
    risk = c(10, 0, 0, 10, 20, 30),
    status = "success",
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(
    list(fold = 1L, time = c(1, 2), event = c(1L, 1L)),
    list(fold = 2L, time = c(1, 2), event = c(1L, 1L))
  )
  curves <- list(`1` = curve, `2` = curve)
  result <- sb_survival_uno_c_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 6,
    censoring_curves = curves
  )

  expect_identical(result$status, "success")
  expect_equal(result$weighted_pair_credit, 1)
  expect_equal(result$weighted_pair_denominator, 4)
  expect_equal(result$estimate, 0.25)
  expect_false(isTRUE(all.equal(result$estimate, mean(c(1, 0)))))

  pooled <- literal_uno_pairs(
    predictions$time,
    predictions$event,
    predictions$risk,
    curve,
    6
  )
  expect_false(isTRUE(all.equal(result$estimate, pooled$estimate)))

  shifted <- predictions
  shifted$risk[shifted$fold == 2L] <-
    shifted$risk[shifted$fold == 2L] + 1000
  shifted_result <- sb_survival_uno_c_metric(
    shifted,
    analysis_inputs,
    prediction_horizon = 6,
    censoring_curves = curves
  )
  expect_equal(shifted_result$estimate, result$estimate)
})

test_that("Uno C requires exactly one analysis input and curve per fold", {
  predictions <- data.frame(
    patient_id = c("p1", "p2"),
    fold = c(1L, 1L),
    time = c(1, 2),
    event = c(1L, 0L),
    risk = c(2, 1),
    status = "success",
    stringsAsFactors = FALSE
  )
  curve <- uno_test_curve(c(0.5, 3), c(1, 1))

  missing_input <- sb_survival_uno_c_metric(
    predictions, list(), 3, list(`1` = curve)
  )
  expect_identical(
    missing_input$per_fold$error_class,
    "stormbreaker_bad_uno_c_analysis_input"
  )

  duplicate_input <- sb_survival_uno_c_metric(
    predictions,
    list(
      list(fold = 1L, time = 1, event = 1L),
      list(fold = 1L, time = 2, event = 0L)
    ),
    3,
    list(`1` = curve)
  )
  expect_identical(
    duplicate_input$per_fold$error_class,
    "stormbreaker_bad_uno_c_analysis_input"
  )

  missing_curve <- sb_survival_uno_c_metric(
    predictions,
    list(list(fold = 1L, time = c(1, 2), event = c(1L, 0L))),
    3,
    list()
  )
  expect_identical(
    missing_curve$per_fold$error_class,
    "stormbreaker_missing_uno_c_censoring_curve"
  )

  duplicate_curves <- list(curve, curve)
  names(duplicate_curves) <- c("1", "1")
  duplicate_curve <- sb_survival_uno_c_metric(
    predictions,
    list(list(fold = 1L, time = c(1, 2), event = c(1L, 0L))),
    3,
    duplicate_curves
  )
  expect_identical(
    duplicate_curve$per_fold$error_class,
    "stormbreaker_bad_uno_c_censoring_curve"
  )
})

test_that("unexpected Uno errors become fold-specific metric evidence", {
  predictions <- data.frame(
    patient_id = c("p1", "p2"),
    fold = c(1L, 1L),
    time = c(1, 2),
    event = c(1L, 0L),
    risk = c(2, 1),
    status = "success",
    stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    sb_uno_c_stratum = function(...) {
      condition <- simpleError("injected fold Uno failure")
      class(condition) <- c("test_fold_uno_error", class(condition))
      stop(condition)
    }
  )

  result <- sb_survival_uno_c_metric(
    predictions,
    list(list(fold = 1L, time = c(1, 2), event = c(1L, 0L))),
    3,
    list(`1` = uno_test_curve(c(0.5, 3), c(1, 1)))
  )

  expect_identical(result$status, "incomplete")
  expect_true(is.na(result$estimate))
  expect_identical(result$per_fold$status, "metric_failed")
  expect_identical(result$per_fold$error_class, "test_fold_uno_error")
  expect_identical(
    result$per_fold$error_message,
    "injected fold Uno failure"
  )
})
})
