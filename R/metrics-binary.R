# ---- metrics-binary ----

# ---- metrics-binary-decision-curve ----

sb_binary_decision_curve_metric <- function(
    predictions,
    decision_thresholds,
    scope) {
  decision_thresholds <- sb_validate_decision_thresholds(
    decision_thresholds
  )
  n <- as.integer(nrow(predictions))
  outcome <- predictions$outcome
  probability <- predictions$probability
  prediction_status <- predictions$status
  n_positive <- as.integer(sum(outcome == 1L, na.rm = TRUE))
  n_negative <- as.integer(sum(outcome == 0L, na.rm = TRUE))
  threshold_odds <- decision_thresholds / (1 - decision_thresholds)

  curve <- data.frame(
    threshold_probability = decision_thresholds,
    threshold_odds = threshold_odds,
    n = rep.int(n, length(decision_thresholds)),
    n_positive = rep.int(n_positive, length(decision_thresholds)),
    n_negative = rep.int(n_negative, length(decision_thresholds)),
    n_selected = rep.int(NA_integer_, length(decision_thresholds)),
    true_positive = rep.int(NA_integer_, length(decision_thresholds)),
    false_positive = rep.int(NA_integer_, length(decision_thresholds)),
    true_negative = rep.int(NA_integer_, length(decision_thresholds)),
    false_negative = rep.int(NA_integer_, length(decision_thresholds)),
    sensitivity = rep.int(NA_real_, length(decision_thresholds)),
    specificity = rep.int(NA_real_, length(decision_thresholds)),
    positive_predictive_value = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    negative_predictive_value = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    model_net_benefit = rep.int(NA_real_, length(decision_thresholds)),
    treat_all_net_benefit = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    treat_none_net_benefit = rep.int(NA_real_, length(decision_thresholds)),
    stringsAsFactors = FALSE
  )
  result <- list(
    metric = "decision_curve",
    task_type = "binary",
    scope = scope,
    purpose = "descriptive_clinical_utility_evaluation",
    threshold_source = "study_declared_binary_task",
    decision_rule = "probability_greater_than_or_equal_to_threshold",
    net_benefit_formula = paste0(
      "true_positive / n - false_positive / n * threshold / ",
      "(1 - threshold)"
    ),
    reference_strategies = c("treat_all", "treat_none"),
    status = "incomplete",
    n = n,
    n_scored = 0L,
    n_positive = n_positive,
    n_negative = n_negative,
    n_thresholds = as.integer(length(decision_thresholds)),
    curve = curve,
    error_class = "stormbreaker_incomplete_decision_curve_input",
    error_message = paste0(
      "Binary outcomes, probabilities, or prediction statuses are ",
      "incomplete or invalid."
    )
  )

  complete <- n >= 1L &&
    (is.numeric(outcome) || is.logical(outcome)) &&
    length(outcome) == n && !anyNA(outcome) &&
    (!is.numeric(outcome) || all(is.finite(outcome))) &&
    all(outcome %in% c(0, 1)) &&
    is.numeric(probability) && length(probability) == n &&
    !anyNA(probability) && all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1) &&
    is.character(prediction_status) && length(prediction_status) == n &&
    !anyNA(prediction_status) && all(prediction_status == "success")
  if (!complete) {
    return(result)
  }

  outcome <- as.integer(outcome)
  true_positive <- vapply(decision_thresholds, function(threshold) {
    as.integer(sum(probability >= threshold & outcome == 1L))
  }, integer(1))
  false_positive <- vapply(decision_thresholds, function(threshold) {
    as.integer(sum(probability >= threshold & outcome == 0L))
  }, integer(1))
  true_negative <- vapply(decision_thresholds, function(threshold) {
    as.integer(sum(probability < threshold & outcome == 0L))
  }, integer(1))
  false_negative <- vapply(decision_thresholds, function(threshold) {
    as.integer(sum(probability < threshold & outcome == 1L))
  }, integer(1))
  n_selected <- true_positive + false_positive
  n_not_selected <- true_negative + false_negative

  result$status <- "success"
  result$n_scored <- n
  result$curve$n_selected <- n_selected
  result$curve$true_positive <- true_positive
  result$curve$false_positive <- false_positive
  result$curve$true_negative <- true_negative
  result$curve$false_negative <- false_negative
  result$curve$sensitivity <- if (n_positive > 0L) {
    true_positive / n_positive
  } else {
    rep.int(NA_real_, length(decision_thresholds))
  }
  result$curve$specificity <- if (n_negative > 0L) {
    true_negative / n_negative
  } else {
    rep.int(NA_real_, length(decision_thresholds))
  }
  has_selected <- n_selected > 0L
  has_not_selected <- n_not_selected > 0L
  result$curve$positive_predictive_value[has_selected] <-
    true_positive[has_selected] / n_selected[has_selected]
  result$curve$negative_predictive_value[has_not_selected] <-
    true_negative[has_not_selected] / n_not_selected[has_not_selected]
  result$curve$model_net_benefit <- true_positive / n -
    false_positive / n * threshold_odds
  result$curve$treat_all_net_benefit <- n_positive / n -
    n_negative / n * threshold_odds
  result$curve$treat_none_net_benefit <- rep.int(
    0,
    length(decision_thresholds)
  )
  result$error_class <- NA_character_
  result$error_message <- NA_character_
  result
}

# ---- metrics-binary-moderate-calibration ----

sb_empty_binary_moderate_calibration_knots <- function() {
  data.frame(
    role = character(),
    quantile_probability = numeric(),
    logit_probability = numeric(),
    probability = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_empty_binary_moderate_calibration_patients <- function() {
  data.frame(
    row_instance = integer(),
    patient_id = character(),
    fold = integer(),
    outcome = integer(),
    probability = numeric(),
    logit_probability = numeric(),
    calibrated_probability = numeric(),
    absolute_error = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_empty_binary_moderate_calibration_curve <- function() {
  data.frame(
    probability = numeric(),
    logit_probability = numeric(),
    calibrated_probability = numeric(),
    n_row_instances = integer(),
    stringsAsFactors = FALSE
  )
}

sb_empty_binary_moderate_calibration_folds <- function() {
  data.frame(
    fold = integer(),
    status = character(),
    n = integer(),
    n_scored = integer(),
    n_positive = integer(),
    n_negative = integer(),
    n_boundary_probability = integer(),
    error_class = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
}

sb_binary_moderate_calibration_empty_fit <- function() {
  terms <- c("(Intercept)", "calibration_ns1", "calibration_ns2")
  coefficients <- data.frame(
    term = terms,
    estimate = rep(NA_real_, length(terms)),
    stringsAsFactors = FALSE
  )
  names(coefficients$estimate) <- terms
  list(
    status = "not_fitted",
    method = "binomial_glm_three_knot_restricted_natural_cubic_spline",
    family = "binomial",
    link = "logit",
    basis = "splines_ns_intercept_false",
    terms = terms,
    coefficients = coefficients,
    rank = NA_integer_,
    converged = NA,
    boundary = NA,
    iterations = NA_integer_,
    deviance = NA_real_,
    null_deviance = NA_real_,
    separation_method = "bounded_linear_program_signed_margin",
    separation_status = "not_checked",
    separation_detected = NA,
    separation_solver_status = NA_integer_,
    separation_objective = NA_real_,
    separation_tolerance = sqrt(.Machine$double.eps),
    separation_minimum_signed_margin = NA_real_,
    separation_maximum_signed_margin = NA_real_,
    warnings = character(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_binary_moderate_calibration_empty_result <- function(scope, n = 0L) {
  safe_scope <- if (is.character(scope) && length(scope) >= 1L) {
    as.character(scope[[1L]])
  } else {
    NA_character_
  }
  list(
    metric = "binary_calibration_curve",
    calibration_level = "moderate",
    estimand = "smoothed_observed_class_1_probability",
    method = "logistic_three_knot_restricted_natural_cubic_spline",
    scope = safe_scope,
    purpose = "descriptive_evaluation_only",
    prediction_scale = "class_1_probability",
    prediction_transform = "logit",
    probability_clipping = "none",
    positive_class = 1L,
    spline = list(
      type = "restricted_natural_cubic_spline",
      n_knots = 3L,
      knot_probabilities = c(0.10, 0.50, 0.90),
      quantile_type = 7L,
      outcome_tuning = "none",
      basis = "splines_ns",
      intercept = FALSE,
      boundary_knots = c(NA_real_, NA_real_),
      interior_knots = NA_real_
    ),
    status = "incomplete",
    input_status = "not_checked",
    ici = NA_real_,
    e50 = NA_real_,
    e90 = NA_real_,
    n = as.integer(n),
    n_scored = 0L,
    n_positive = 0L,
    n_negative = 0L,
    n_boundary_probability = 0L,
    n_unique_probability = 0L,
    n_folds = 0L,
    n_evaluable_folds = 0L,
    knots = sb_empty_binary_moderate_calibration_knots(),
    fit = sb_binary_moderate_calibration_empty_fit(),
    curve = sb_empty_binary_moderate_calibration_curve(),
    patient_calibration = sb_empty_binary_moderate_calibration_patients(),
    per_fold = sb_empty_binary_moderate_calibration_folds(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_binary_moderate_calibration_fail <- function(
    result,
    status,
    input_status,
    error_class,
    error_message) {
  result$status <- status
  result$input_status <- input_status
  result$error_class <- error_class
  result$error_message <- error_message
  result
}

sb_binary_moderate_calibration_error_metric <- function(error, scope, n = 0L) {
  result <- sb_binary_moderate_calibration_empty_result(scope, n)
  result$status <- "metric_failed"
  result$input_status <- "unexpected_error"
  result$error_class <- class(error)[[1L]]
  result$error_message <- conditionMessage(error)
  result
}

sb_binary_moderate_calibration_fold_audit <- function(patients, status) {
  if (all(is.na(patients$fold))) {
    return(sb_empty_binary_moderate_calibration_folds())
  }
  fold_ids <- sort(unique(patients$fold))
  rows <- lapply(fold_ids, function(fold_id) {
    index <- patients$fold == fold_id
    successful <- all(status[index] == "success")
    data.frame(
      fold = as.integer(fold_id),
      status = if (successful) "success" else "prediction_failed",
      n = as.integer(sum(index)),
      n_scored = if (successful) as.integer(sum(index)) else 0L,
      n_positive = as.integer(sum(patients$outcome[index] == 1L)),
      n_negative = as.integer(sum(patients$outcome[index] == 0L)),
      n_boundary_probability = as.integer(sum(
        patients$probability[index] %in% c(0, 1),
        na.rm = TRUE
      )),
      error_class = if (successful) {
        NA_character_
      } else {
        "stormbreaker_invalid_binary_moderate_calibration_prediction"
      },
      error_message = if (successful) {
        NA_character_
      } else {
        "One or more saved predictions in this fold failed."
      },
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

sb_binary_moderate_calibration_knots <- function(logit_probability) {
  values <- unname(stats::quantile(
    logit_probability,
    probs = c(0.10, 0.50, 0.90),
    names = FALSE,
    type = 7
  ))
  data.frame(
    role = c("lower_boundary", "interior", "upper_boundary"),
    quantile_probability = c(0.10, 0.50, 0.90),
    logit_probability = values,
    probability = stats::plogis(values),
    stringsAsFactors = FALSE
  )
}

sb_binary_moderate_calibration_basis <- function(logit_probability, knots) {
  basis <- splines::ns(
    logit_probability,
    knots = knots[[2L]],
    Boundary.knots = knots[c(1L, 3L)],
    intercept = FALSE
  )
  colnames(basis) <- c("calibration_ns1", "calibration_ns2")
  basis
}

sb_binary_moderate_calibration_separation <- function(design, outcome) {
  column_scale <- apply(abs(design), 2L, max)
  normalized_design <- sweep(design, 2L, column_scale, "/")
  outcome_sign <- ifelse(outcome == 1L, 1, -1)
  signed_design <- sweep(normalized_design, 1L, outcome_sign, "*")
  n_terms <- ncol(signed_design)
  objective <- c(colSums(signed_design), -colSums(signed_design))
  constraints <- rbind(
    cbind(signed_design, -signed_design),
    rep(1, 2L * n_terms)
  )
  solver <- tryCatch(
    lpSolve::lp(
      direction = "max",
      objective.in = objective,
      const.mat = constraints,
      const.dir = c(rep(">=", nrow(signed_design)), "<="),
      const.rhs = c(rep(0, nrow(signed_design)), 1)
    ),
    error = function(error) error
  )
  tolerance <- sqrt(.Machine$double.eps)
  if (inherits(solver, "condition")) {
    return(list(
      status = "failed",
      detected = NA,
      solver_status = NA_integer_,
      objective = NA_real_,
      tolerance = tolerance,
      minimum_signed_margin = NA_real_,
      maximum_signed_margin = NA_real_,
      error_class = class(solver)[[1L]],
      error_message = conditionMessage(solver)
    ))
  }
  if (!identical(as.integer(solver$status), 0L)) {
    return(list(
      status = "failed",
      detected = NA,
      solver_status = as.integer(solver$status),
      objective = unname(solver$objval),
      tolerance = tolerance,
      minimum_signed_margin = NA_real_,
      maximum_signed_margin = NA_real_,
      error_class =
        "stormbreaker_binary_moderate_calibration_separation_check_failed",
      error_message = paste0(
        "The bounded separation linear program returned solver status ",
        as.integer(solver$status), "."
      )
    ))
  }

  positive <- solver$solution[seq_len(n_terms)]
  negative <- solver$solution[n_terms + seq_len(n_terms)]
  separator <- positive - negative
  signed_margin <- drop(signed_design %*% separator)
  objective_value <- unname(solver$objval)
  minimum_margin <- unname(min(signed_margin))
  maximum_margin <- unname(max(signed_margin))
  feasible <- minimum_margin >= -tolerance
  detected <- objective_value > tolerance && feasible &&
    maximum_margin > tolerance
  if (objective_value > tolerance && !feasible) {
    return(list(
      status = "failed",
      detected = NA,
      solver_status = as.integer(solver$status),
      objective = objective_value,
      tolerance = tolerance,
      minimum_signed_margin = minimum_margin,
      maximum_signed_margin = maximum_margin,
      error_class =
        "stormbreaker_binary_moderate_calibration_separation_check_failed",
      error_message = paste0(
        "The bounded separation linear program returned a solution outside ",
        "its signed-margin tolerance."
      )
    ))
  }
  list(
    status = if (detected) "detected" else "not_detected",
    detected = detected,
    solver_status = as.integer(solver$status),
    objective = objective_value,
    tolerance = tolerance,
    minimum_signed_margin = minimum_margin,
    maximum_signed_margin = maximum_margin,
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_binary_moderate_calibration_curve_table <- function(patients) {
  probability <- sort(unique(patients$probability))
  first <- match(probability, patients$probability)
  data.frame(
    probability = probability,
    logit_probability = patients$logit_probability[first],
    calibrated_probability = patients$calibrated_probability[first],
    n_row_instances = as.integer(tabulate(
      match(patients$probability, probability),
      nbins = length(probability)
    )),
    stringsAsFactors = FALSE
  )
}

sb_binary_moderate_calibration_metric <- function(scored, scope) {
  n <- if (is.data.frame(scored)) nrow(scored) else 0L
  result <- sb_binary_moderate_calibration_empty_result(scope, n)
  required <- c("patient_id", "outcome", "probability", "status")
  required_once <- is.data.frame(scored) && n >= 1L &&
    all(vapply(required, function(name) {
      sum(names(scored) == name) == 1L
    }, logical(1)))
  fold_count <- if (is.data.frame(scored)) sum(names(scored) == "fold") else 0L
  valid_scope <- is.character(scope) && length(scope) == 1L &&
    !is.na(scope) && nzchar(scope)
  if (!required_once || fold_count > 1L || !valid_scope) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "input_failed",
      "stormbreaker_bad_binary_moderate_calibration_input",
      paste0(
        "Binary moderate calibration requires one non-empty scored table, ",
        "one valid scope, required columns, and at most one fold column."
      )
    ))
  }

  patient_id <- scored$patient_id
  outcome <- scored$outcome
  probability <- scored$probability
  prediction_status <- scored$status
  fold <- if (fold_count == 1L) scored$fold else rep(NA_integer_, n)
  valid_patient_id <- is.character(patient_id) && length(patient_id) == n &&
    !anyNA(patient_id) && all(nzchar(trimws(patient_id)))
  valid_outcome <- (is.numeric(outcome) || is.logical(outcome)) &&
    length(outcome) == n && !anyNA(outcome) &&
    (!is.numeric(outcome) || all(is.finite(outcome))) &&
    all(outcome %in% c(0, 1))
  valid_probability_type <- is.numeric(probability) &&
    length(probability) == n
  valid_status <- is.character(prediction_status) &&
    length(prediction_status) == n && !anyNA(prediction_status)
  valid_fold <- fold_count == 0L || (
    is.numeric(fold) && length(fold) == n && !anyNA(fold) &&
      all(is.finite(fold)) && all(fold == trunc(fold))
  )
  if (!valid_patient_id || !valid_outcome || !valid_probability_type ||
      !valid_status || !valid_fold) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "input_failed",
      "stormbreaker_bad_binary_moderate_calibration_input",
      paste0(
        "Binary moderate calibration requires valid aligned patient IDs, ",
        "0/1 outcomes, probabilities, statuses, and optional integer folds."
      )
    ))
  }

  valid_probability <- !anyNA(probability) && all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  boundary <- !is.na(probability) & is.finite(probability) &
    probability %in% c(0, 1)
  interior <- !is.na(probability) & is.finite(probability) &
    probability > 0 & probability < 1
  logit_probability <- rep(NA_real_, n)
  logit_probability[interior] <- stats::qlogis(probability[interior])
  patients <- data.frame(
    row_instance = seq_len(n),
    patient_id = patient_id,
    fold = if (fold_count == 1L) as.integer(fold) else rep(NA_integer_, n),
    outcome = as.integer(outcome),
    probability = unname(probability),
    logit_probability = logit_probability,
    calibrated_probability = rep(NA_real_, n),
    absolute_error = rep(NA_real_, n),
    stringsAsFactors = FALSE
  )
  result$patient_calibration <- patients
  result$n_positive <- as.integer(sum(patients$outcome == 1L))
  result$n_negative <- as.integer(sum(patients$outcome == 0L))
  result$n_boundary_probability <- as.integer(sum(boundary))
  result$n_unique_probability <- as.integer(length(unique(
    probability[!is.na(probability) & is.finite(probability)]
  )))
  result$per_fold <- sb_binary_moderate_calibration_fold_audit(
    patients,
    prediction_status
  )
  result$n_folds <- as.integer(nrow(result$per_fold))
  result$n_evaluable_folds <- as.integer(sum(
    result$per_fold$status == "success"
  ))

  if (any(prediction_status != "success")) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "prediction_failed",
      "stormbreaker_invalid_binary_moderate_calibration_prediction",
      "Complete successful probabilities are required for calibration."
    ))
  }
  if (!valid_probability) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "input_failed",
      "stormbreaker_bad_binary_moderate_calibration_input",
      "Successful binary predictions must be finite probabilities in [0, 1]."
    ))
  }
  result$n_scored <- as.integer(n)
  if (result$n_positive == 0L || result$n_negative == 0L) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_binary_moderate_calibration_one_class",
      "Both outcome classes are required for binary moderate calibration."
    ))
  }
  if (result$n_boundary_probability > 0L) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_binary_moderate_calibration_boundary_probability",
      paste0(
        "Exact probabilities zero or one cannot enter a finite-logit ",
        "calibration fit without clipping, which is not performed."
      )
    ))
  }
  if (any(!is.finite(logit_probability))) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_binary_moderate_calibration_nonfinite_logit",
      "The unmodified probabilities did not produce finite logits."
    ))
  }

  knots <- sb_binary_moderate_calibration_knots(logit_probability)
  result$knots <- knots
  knot_values <- knots$logit_probability
  result$spline$boundary_knots <- unname(knot_values[c(1L, 3L)])
  result$spline$interior_knots <- unname(knot_values[[2L]])
  if (anyNA(knot_values) || any(!is.finite(knot_values)) ||
      any(diff(knot_values) <= 0)) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_binary_moderate_calibration_duplicate_knots",
      paste0(
        "The fixed 10th, 50th, and 90th percentile logit knots are not ",
        "finite and strictly increasing."
      )
    ))
  }

  basis <- tryCatch(
    sb_binary_moderate_calibration_basis(logit_probability, knot_values),
    error = function(error) error
  )
  if (inherits(basis, "condition")) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      class(basis)[[1L]],
      conditionMessage(basis)
    ))
  }
  design <- cbind("(Intercept)" = 1, basis)
  captured_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm.fit(
        x = design,
        y = patients$outcome,
        family = stats::binomial(),
        control = stats::glm.control()
      ),
      warning = function(warning) {
        captured_warnings <<- c(captured_warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  result$fit$warnings <- unique(captured_warnings)
  if (inherits(fit, "condition")) {
    result$fit$status <- "failed"
    result$fit$error_class <- class(fit)[[1L]]
    result$fit$error_message <- conditionMessage(fit)
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }

  terms <- result$fit$terms
  coefficient_values <- unname(fit$coefficients[terms])
  names(coefficient_values) <- terms
  result$fit$coefficients$estimate <- coefficient_values
  result$fit$rank <- as.integer(fit$rank)
  result$fit$converged <- isTRUE(fit$converged)
  result$fit$boundary <- isTRUE(fit$boundary)
  result$fit$iterations <- as.integer(fit$iter)
  result$fit$deviance <- unname(fit$deviance)
  result$fit$null_deviance <- unname(fit$null.deviance)
  # Saturated fitted probabilities alone do not establish separation. Keep
  # this known warning and still require every fit and separation check below.
  unexpected_warnings <- setdiff(
    captured_warnings,
    "glm.fit: fitted probabilities numerically 0 or 1 occurred"
  )
  if (length(unexpected_warnings) > 0L) {
    result$fit$status <- "failed"
    result$fit$error_class <-
      "stormbreaker_binary_moderate_calibration_glm_warning"
    result$fit$error_message <- paste(unique(captured_warnings), collapse = " ")
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }
  if (!result$fit$converged) {
    result$fit$status <- "failed"
    result$fit$error_class <-
      "stormbreaker_binary_moderate_calibration_nonconvergence"
    result$fit$error_message <- "The binary calibration GLM did not converge."
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }
  if (result$fit$boundary) {
    result$fit$status <- "failed"
    result$fit$error_class <-
      "stormbreaker_binary_moderate_calibration_boundary_fit"
    result$fit$error_message <- "The binary calibration GLM stopped at a boundary."
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }
  valid_coefficients <- identical(result$fit$rank, 3L) &&
    length(coefficient_values) == 3L && !anyNA(coefficient_values) &&
    all(is.finite(coefficient_values))
  if (!valid_coefficients) {
    result$fit$status <- "failed"
    result$fit$error_class <-
      "stormbreaker_binary_moderate_calibration_singular_fit"
    result$fit$error_message <- paste0(
      "The binary calibration GLM returned missing, singular, or non-finite ",
      "basis coefficients."
    )
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }

  # IRLS can report convergence with finite coefficients even when separation
  # makes the logistic maximum-likelihood estimate nonexistent. This post-fit
  # linear program checks that data-dependent failure without using a
  # coefficient cutoff or replacing the requested calibration model.
  separation <- sb_binary_moderate_calibration_separation(
    design,
    patients$outcome
  )
  result$fit$separation_status <- separation$status
  result$fit$separation_detected <- separation$detected
  result$fit$separation_solver_status <- separation$solver_status
  result$fit$separation_objective <- separation$objective
  result$fit$separation_tolerance <- separation$tolerance
  result$fit$separation_minimum_signed_margin <-
    separation$minimum_signed_margin
  result$fit$separation_maximum_signed_margin <-
    separation$maximum_signed_margin
  if (identical(separation$status, "failed")) {
    result$fit$status <- "failed"
    result$fit$error_class <- separation$error_class
    result$fit$error_message <- separation$error_message
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }
  if (isTRUE(separation$detected)) {
    result$fit$status <- "failed"
    result$fit$error_class <-
      "stormbreaker_binary_moderate_calibration_separation"
    result$fit$error_message <- paste0(
      "The binary calibration design has complete or quasi-complete outcome ",
      "separation, so no finite logistic maximum-likelihood estimate is ",
      "reported."
    )
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }

  calibrated_probability <- stats::plogis(drop(
    design %*% coefficient_values
  ))
  valid_calibrated <- !anyNA(calibrated_probability) &&
    all(is.finite(calibrated_probability)) &&
    all(calibrated_probability >= 0 & calibrated_probability <= 1)
  if (!valid_calibrated) {
    result$fit$status <- "failed"
    result$fit$error_class <-
      "stormbreaker_binary_moderate_calibration_invalid_probability"
    result$fit$error_message <- paste0(
      "The binary calibration GLM produced invalid probabilities."
    )
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }

  absolute_error <- abs(calibrated_probability - probability)
  result$patient_calibration$calibrated_probability <- calibrated_probability
  result$patient_calibration$absolute_error <- absolute_error
  result$curve <- sb_binary_moderate_calibration_curve_table(
    result$patient_calibration
  )
  result$ici <- unname(mean(absolute_error))
  result$e50 <- unname(stats::quantile(
    absolute_error,
    probs = 0.50,
    names = FALSE,
    type = 7
  ))
  result$e90 <- unname(stats::quantile(
    absolute_error,
    probs = 0.90,
    names = FALSE,
    type = 7
  ))
  result$fit$status <- "success"
  result$status <- "success"
  result$input_status <- "success"
  result
}

# ---- metrics-binary ----

sb_binary_auc_stratum <- function(outcome, probability) {
  outcome <- as.integer(outcome)
  n_positive <- as.integer(sum(outcome == 1L))
  n_negative <- as.integer(sum(outcome == 0L))
  comparable_pairs <- unname(n_positive * n_negative)
  result <- list(
    status = "not_evaluable",
    estimate = NA_real_,
    n = as.integer(length(outcome)),
    n_positive = n_positive,
    n_negative = n_negative,
    correctly_ordered_pairs = 0,
    tied_pairs = 0,
    comparable_pairs = comparable_pairs
  )
  if (n_positive == 0L || n_negative == 0L) {
    return(result)
  }

  ranks <- rank(probability, ties.method = "average")
  pair_credit <- sum(ranks[outcome == 1L]) -
    n_positive * (n_positive + 1) / 2
  ordered <- order(probability)
  runs <- rle(probability[ordered])
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  ordered_outcome <- outcome[ordered]
  tied_pairs <- sum(vapply(seq_along(starts), function(index) {
    values <- ordered_outcome[starts[[index]]:ends[[index]]]
    sum(values == 1L) * sum(values == 0L)
  }, numeric(1)))

  result$status <- "success"
  result$estimate <- unname(pair_credit / comparable_pairs)
  result$correctly_ordered_pairs <- unname(pair_credit - 0.5 * tied_pairs)
  result$tied_pairs <- unname(tied_pairs)
  result
}

sb_binary_auc_metric <- function(predictions) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  per_fold <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      c("outcome", "probability", "status"),
      drop = FALSE
    ]
    audit <- data.frame(
      fold = fold_id,
      status = "prediction_failed",
      n = nrow(fold_predictions),
      n_positive = as.integer(sum(fold_predictions$outcome == 1L)),
      n_negative = as.integer(sum(fold_predictions$outcome == 0L)),
      estimate = NA_real_,
      correctly_ordered_pairs = 0,
      tied_pairs = 0,
      comparable_pairs = 0,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    prediction_failed <- anyNA(fold_predictions$status) ||
      any(fold_predictions$status != "success") ||
      anyNA(fold_predictions$probability) ||
      any(!is.finite(fold_predictions$probability))
    if (prediction_failed) {
      audit$error_message <- paste0(
        "OOF probabilities for this fold are failed or non-finite."
      )
      return(audit)
    }
    stratum <- sb_binary_auc_stratum(
      outcome = fold_predictions$outcome,
      probability = fold_predictions$probability
    )
    audit$n <- stratum$n
    audit$n_positive <- stratum$n_positive
    audit$n_negative <- stratum$n_negative
    audit$correctly_ordered_pairs <- stratum$correctly_ordered_pairs
    audit$tied_pairs <- stratum$tied_pairs
    audit$comparable_pairs <- stratum$comparable_pairs
    if (!identical(stratum$status, "success")) {
      audit$status <- "not_evaluable"
      audit$error_message <- "Both outcome classes are required in this fold."
      return(audit)
    }

    audit$status <- "success"
    audit$estimate <- stratum$estimate
    audit
  })
  per_fold <- do.call(rbind, per_fold)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  totals <- c(
    correctly_ordered_pairs = 0,
    tied_pairs = 0,
    comparable_pairs = 0
  )
  if (any(evaluable)) {
    totals <- colSums(per_fold[
      evaluable,
      names(totals),
      drop = FALSE
    ])
  }
  complete <- length(fold_ids) > 0L && all(evaluable)
  estimate <- if (complete) {
    unname(
      (totals[["correctly_ordered_pairs"]] +
         0.5 * totals[["tied_pairs"]]) /
        totals[["comparable_pairs"]]
    )
  } else {
    NA_real_
  }

  list(
    metric = "roc_auc",
    scope = "within_outer_assessment_fold",
    aggregation = "summed_positive_negative_pairs",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = if (complete) "success" else "incomplete",
    estimate = estimate,
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    correctly_ordered_pairs = unname(
      totals[["correctly_ordered_pairs"]]
    ),
    tied_pairs = unname(totals[["tied_pairs"]]),
    comparable_pairs = unname(totals[["comparable_pairs"]]),
    per_fold = per_fold
  )
}

sb_empty_binary_precision_recall_curve <- function() {
  data.frame(
    fold = integer(),
    threshold = numeric(),
    n_at_threshold = integer(),
    positive_at_threshold = integer(),
    negative_at_threshold = integer(),
    true_positive = integer(),
    false_positive = integer(),
    precision = numeric(),
    recall = numeric(),
    recall_increment = numeric(),
    contribution = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_binary_precision_recall_curve <- function(
    outcome,
    probability,
    fold_id) {
  ordered <- order(probability, decreasing = TRUE)
  ordered_probability <- probability[ordered]
  ordered_outcome <- as.integer(outcome[ordered])
  score_runs <- rle(ordered_probability)
  ends <- cumsum(score_runs$lengths)
  starts <- ends - score_runs$lengths + 1L
  positive_at_threshold <- vapply(seq_along(starts), function(index) {
    sum(ordered_outcome[starts[[index]]:ends[[index]]] == 1L)
  }, integer(1))
  negative_at_threshold <- score_runs$lengths - positive_at_threshold
  true_positive <- cumsum(positive_at_threshold)
  false_positive <- cumsum(negative_at_threshold)
  n_positive <- sum(ordered_outcome == 1L)
  precision <- true_positive / (true_positive + false_positive)
  recall <- true_positive / n_positive
  recall_increment <- c(recall[[1L]], diff(recall))

  data.frame(
    fold = rep(as.integer(fold_id), length(score_runs$values)),
    threshold = as.numeric(score_runs$values),
    n_at_threshold = as.integer(score_runs$lengths),
    positive_at_threshold = as.integer(positive_at_threshold),
    negative_at_threshold = as.integer(negative_at_threshold),
    true_positive = as.integer(true_positive),
    false_positive = as.integer(false_positive),
    precision = unname(precision),
    recall = unname(recall),
    recall_increment = unname(recall_increment),
    contribution = unname(recall_increment * precision),
    stringsAsFactors = FALSE
  )
}

sb_binary_average_precision_metric <- function(predictions) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  fold_results <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      c("outcome", "probability", "status"),
      drop = FALSE
    ]
    audit <- data.frame(
      fold = fold_id,
      status = "prediction_failed",
      n = as.integer(nrow(fold_predictions)),
      n_scored = 0L,
      n_positive = as.integer(sum(fold_predictions$outcome == 1L)),
      n_negative = as.integer(sum(fold_predictions$outcome == 0L)),
      n_thresholds = 0L,
      positive_weighted_sum = 0,
      estimate = NA_real_,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    curve <- sb_empty_binary_precision_recall_curve()
    prediction_failed <- anyNA(fold_predictions$status) ||
      any(fold_predictions$status != "success") ||
      anyNA(fold_predictions$outcome) ||
      any(!fold_predictions$outcome %in% c(0L, 1L)) ||
      anyNA(fold_predictions$probability) ||
      any(!is.finite(fold_predictions$probability)) ||
      any(
        fold_predictions$probability < 0 |
          fold_predictions$probability > 1
      )
    if (prediction_failed) {
      audit$error_message <- paste0(
        "OOF outcomes or probabilities for this fold are failed or invalid."
      )
      return(list(audit = audit, curve = curve))
    }
    if (audit$n_positive == 0L || audit$n_negative == 0L) {
      audit$status <- "not_evaluable"
      audit$error_message <- "Both outcome classes are required in this fold."
      return(list(audit = audit, curve = curve))
    }

    curve <- sb_binary_precision_recall_curve(
      as.integer(fold_predictions$outcome),
      fold_predictions$probability,
      fold_id
    )
    positive_weighted_sum <- sum(
      curve$positive_at_threshold * curve$precision
    )
    audit$status <- "success"
    audit$n_scored <- audit$n
    audit$n_thresholds <- as.integer(nrow(curve))
    audit$positive_weighted_sum <- unname(positive_weighted_sum)
    audit$estimate <- unname(positive_weighted_sum / audit$n_positive)
    list(audit = audit, curve = curve)
  })
  per_fold <- do.call(rbind, lapply(fold_results, `[[`, "audit"))
  rownames(per_fold) <- NULL
  curves <- lapply(fold_results, `[[`, "curve")
  names(curves) <- as.character(fold_ids)
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  positive_weighted_sum <- unname(sum(
    per_fold$positive_weighted_sum[evaluable]
  ))
  n_positive <- as.integer(sum(per_fold$n_positive))

  list(
    metric = "average_precision",
    scope = "within_outer_assessment_fold",
    aggregation = "positive_weighted_mean_of_fold_average_precision",
    definition = "sum((recall_i - recall_previous) * precision_i)",
    interpolation = "none",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = if (complete) "success" else "incomplete",
    estimate = if (complete) {
      unname(positive_weighted_sum / n_positive)
    } else {
      NA_real_
    },
    n = as.integer(sum(per_fold$n)),
    n_scored = as.integer(sum(per_fold$n_scored[evaluable])),
    n_positive = n_positive,
    n_negative = as.integer(sum(per_fold$n_negative)),
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    positive_weighted_sum = positive_weighted_sum,
    per_fold = per_fold,
    curves = curves
  )
}

sb_binary_brier_reference <- function(predictions, fold_ids, reference_input) {
  n <- as.integer(nrow(predictions))
  if (is.null(reference_input)) {
    return(list(
      type = "outer_analysis_prevalence",
      source = "outer_analysis_outcome_only",
      status = "not_requested",
      estimate = NA_real_,
      n = n,
      n_scored = 0L,
      n_folds = as.integer(length(fold_ids)),
      n_evaluable_folds = 0L,
      squared_error_sum = 0,
      per_fold = NULL
    ))
  }

  valid_schema <- is.data.frame(reference_input) &&
    all(c("fold", "analysis_prevalence") %in% names(reference_input))
  per_fold <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      c("outcome"),
      drop = FALSE
    ]
    audit <- data.frame(
      fold = fold_id,
      status = "reference_failed",
      analysis_prevalence = NA_real_,
      n = as.integer(nrow(fold_predictions)),
      n_scored = 0L,
      squared_error_sum = 0,
      estimate = NA_real_,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    matched <- if (valid_schema) {
      which(as.integer(reference_input$fold) == fold_id)
    } else {
      integer()
    }
    if (length(matched) != 1L) {
      audit$error_message <- paste0(
        "Exactly one outer-analysis prevalence is required for this fold."
      )
      return(audit)
    }
    probability <- reference_input$analysis_prevalence[[matched]]
    valid_probability <- is.numeric(probability) && length(probability) == 1L &&
      !is.na(probability) && is.finite(probability) &&
      probability >= 0 && probability <= 1
    valid_outcome <- !anyNA(fold_predictions$outcome) &&
      all(fold_predictions$outcome %in% c(0L, 1L))
    if (!valid_probability || !valid_outcome) {
      audit$error_message <- paste0(
        "The analysis prevalence or assessment outcomes are invalid."
      )
      return(audit)
    }

    squared_error <- (
      as.integer(fold_predictions$outcome) - probability
    )^2
    audit$status <- "success"
    audit$analysis_prevalence <- probability
    audit$n_scored <- audit$n
    audit$squared_error_sum <- unname(sum(squared_error))
    audit$estimate <- unname(mean(squared_error))
    audit
  })
  per_fold <- do.call(rbind, per_fold)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  n_scored <- as.integer(sum(per_fold$n_scored[evaluable]))
  if (complete) {
    matched <- match(
      as.integer(predictions$fold),
      as.integer(reference_input$fold)
    )
    probability <- reference_input$analysis_prevalence[matched]
    squared_error <- (
      as.integer(predictions$outcome) - probability
    )^2
    squared_error_sum <- unname(sum(squared_error))
    estimate <- unname(mean(squared_error))
  } else {
    squared_error_sum <- unname(sum(
      per_fold$squared_error_sum[evaluable]
    ))
    estimate <- NA_real_
  }

  list(
    type = "outer_analysis_prevalence",
    source = "outer_analysis_outcome_only",
    status = if (complete) "success" else "incomplete",
    estimate = estimate,
    n = n,
    n_scored = n_scored,
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    squared_error_sum = squared_error_sum,
    per_fold = per_fold
  )
}

sb_binary_brier_metric <- function(predictions, reference_input = NULL) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  per_fold <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      c("outcome", "probability", "status"),
      drop = FALSE
    ]
    audit <- data.frame(
      fold = fold_id,
      status = "prediction_failed",
      n = as.integer(nrow(fold_predictions)),
      n_scored = 0L,
      squared_error_sum = 0,
      estimate = NA_real_,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    prediction_failed <- anyNA(fold_predictions$status) ||
      any(fold_predictions$status != "success") ||
      anyNA(fold_predictions$outcome) ||
      any(!fold_predictions$outcome %in% c(0L, 1L)) ||
      anyNA(fold_predictions$probability) ||
      any(!is.finite(fold_predictions$probability)) ||
      any(
        fold_predictions$probability < 0 |
          fold_predictions$probability > 1
      )
    if (prediction_failed) {
      audit$error_message <- paste0(
        "OOF outcomes or probabilities for this fold are failed or invalid."
      )
      return(audit)
    }

    squared_error <- (
      as.integer(fold_predictions$outcome) - fold_predictions$probability
    )^2
    audit$status <- "success"
    audit$n_scored <- audit$n
    audit$squared_error_sum <- unname(sum(squared_error))
    audit$estimate <- unname(mean(squared_error))
    audit
  })
  per_fold <- do.call(rbind, per_fold)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  n <- as.integer(sum(per_fold$n))
  n_scored <- as.integer(sum(per_fold$n_scored[evaluable]))
  if (complete) {
    squared_error <- (
      as.integer(predictions$outcome) - predictions$probability
    )^2
    squared_error_sum <- unname(sum(squared_error))
    estimate <- unname(mean(squared_error))
  } else {
    squared_error_sum <- unname(sum(
      per_fold$squared_error_sum[evaluable]
    ))
    estimate <- NA_real_
  }
  reference <- sb_binary_brier_reference(
    predictions,
    fold_ids,
    reference_input
  )
  skill_complete <- complete && identical(reference$status, "success") &&
    is.finite(reference$squared_error_sum) &&
    reference$squared_error_sum > 0
  skill_score <- list(
    metric = "brier_skill_score",
    reference = reference$type,
    formula = "1 - model_squared_error / reference_squared_error",
    status = if (skill_complete) "success" else "incomplete",
    estimate = if (skill_complete) {
      unname(1 - squared_error_sum / reference$squared_error_sum)
    } else {
      NA_real_
    }
  )

  list(
    metric = "brier_score",
    scope = "patient_level_out_of_fold",
    aggregation = "summed_squared_error_over_patients",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = if (complete) "success" else "incomplete",
    estimate = estimate,
    n = n,
    n_scored = n_scored,
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    squared_error_sum = squared_error_sum,
    per_fold = per_fold,
    reference = reference,
    skill_score = skill_score
  )
}

sb_binary_case_log_loss <- function(outcome, probability) {
  outcome <- as.integer(outcome)
  positive <- outcome == 1L
  case_loss <- numeric(length(outcome))
  case_loss[positive] <- -log(probability[positive])
  case_loss[!positive] <- -log1p(-probability[!positive])
  unname(case_loss)
}

sb_binary_log_loss_reference <- function(predictions, fold_ids, reference_input) {
  n <- as.integer(nrow(predictions))
  if (is.null(reference_input)) {
    return(list(
      type = "outer_analysis_prevalence",
      source = "outer_analysis_outcome_only",
      probability_clipping = "none",
      status = "not_requested",
      estimate = NA_real_,
      n = n,
      n_scored = 0L,
      n_finite_loss = 0L,
      n_infinite_loss = 0L,
      n_folds = as.integer(length(fold_ids)),
      n_evaluable_folds = 0L,
      loss_sum = 0,
      per_fold = NULL
    ))
  }

  valid_schema <- is.data.frame(reference_input) &&
    all(c("fold", "analysis_prevalence") %in% names(reference_input))
  per_fold <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      "outcome",
      drop = FALSE
    ]
    audit <- data.frame(
      fold = fold_id,
      status = "reference_failed",
      analysis_prevalence = NA_real_,
      n = as.integer(nrow(fold_predictions)),
      n_scored = 0L,
      n_finite_loss = 0L,
      n_infinite_loss = 0L,
      loss_sum = 0,
      estimate = NA_real_,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    matched <- if (valid_schema) {
      which(as.integer(reference_input$fold) == fold_id)
    } else {
      integer()
    }
    if (length(matched) != 1L) {
      audit$error_message <- paste0(
        "Exactly one outer-analysis prevalence is required for this fold."
      )
      return(audit)
    }
    probability <- reference_input$analysis_prevalence[[matched]]
    valid_probability <- is.numeric(probability) && length(probability) == 1L &&
      !is.na(probability) && is.finite(probability) &&
      probability >= 0 && probability <= 1
    valid_outcome <- !anyNA(fold_predictions$outcome) &&
      all(fold_predictions$outcome %in% c(0L, 1L))
    if (!valid_probability || !valid_outcome) {
      audit$error_message <- paste0(
        "The analysis prevalence or assessment outcomes are invalid."
      )
      return(audit)
    }

    case_loss <- sb_binary_case_log_loss(
      fold_predictions$outcome,
      rep(probability, nrow(fold_predictions))
    )
    audit$status <- "success"
    audit$analysis_prevalence <- probability
    audit$n_scored <- audit$n
    audit$n_finite_loss <- as.integer(sum(is.finite(case_loss)))
    audit$n_infinite_loss <- as.integer(sum(is.infinite(case_loss)))
    audit$loss_sum <- unname(sum(case_loss))
    audit$estimate <- unname(mean(case_loss))
    audit
  })
  per_fold <- do.call(rbind, per_fold)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  n_scored <- as.integer(sum(per_fold$n_scored[evaluable]))
  loss_sum <- unname(sum(per_fold$loss_sum[evaluable]))

  list(
    type = "outer_analysis_prevalence",
    source = "outer_analysis_outcome_only",
    probability_clipping = "none",
    status = if (complete) "success" else "incomplete",
    estimate = if (complete) unname(loss_sum / n) else NA_real_,
    n = n,
    n_scored = n_scored,
    n_finite_loss = as.integer(sum(per_fold$n_finite_loss[evaluable])),
    n_infinite_loss = as.integer(sum(per_fold$n_infinite_loss[evaluable])),
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    loss_sum = loss_sum,
    per_fold = per_fold
  )
}

sb_binary_log_loss_metric <- function(predictions, reference_input = NULL) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  per_fold <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      c("outcome", "probability", "status"),
      drop = FALSE
    ]
    audit <- data.frame(
      fold = fold_id,
      status = "prediction_failed",
      n = as.integer(nrow(fold_predictions)),
      n_scored = 0L,
      n_finite_loss = 0L,
      n_infinite_loss = 0L,
      loss_sum = 0,
      estimate = NA_real_,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    prediction_failed <- anyNA(fold_predictions$status) ||
      any(fold_predictions$status != "success") ||
      anyNA(fold_predictions$outcome) ||
      any(!fold_predictions$outcome %in% c(0L, 1L)) ||
      anyNA(fold_predictions$probability) ||
      any(!is.finite(fold_predictions$probability)) ||
      any(
        fold_predictions$probability < 0 |
          fold_predictions$probability > 1
      )
    if (prediction_failed) {
      audit$error_message <- paste0(
        "OOF outcomes or probabilities for this fold are failed or invalid."
      )
      return(audit)
    }

    case_loss <- sb_binary_case_log_loss(
      fold_predictions$outcome,
      fold_predictions$probability
    )
    audit$status <- "success"
    audit$n_scored <- audit$n
    audit$n_finite_loss <- as.integer(sum(is.finite(case_loss)))
    audit$n_infinite_loss <- as.integer(sum(is.infinite(case_loss)))
    audit$loss_sum <- unname(sum(case_loss))
    audit$estimate <- unname(mean(case_loss))
    audit
  })
  per_fold <- do.call(rbind, per_fold)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  n <- as.integer(sum(per_fold$n))
  n_scored <- as.integer(sum(per_fold$n_scored[evaluable]))
  loss_sum <- unname(sum(per_fold$loss_sum[evaluable]))
  reference <- sb_binary_log_loss_reference(
    predictions,
    fold_ids,
    reference_input
  )
  skill_complete <- complete && identical(reference$status, "success") &&
    is.finite(reference$loss_sum) && reference$loss_sum > 0
  skill_score <- list(
    metric = "log_loss_skill_score",
    reference = reference$type,
    formula = "1 - model_log_loss / reference_log_loss",
    status = if (skill_complete) "success" else "incomplete",
    estimate = if (skill_complete) {
      unname(1 - loss_sum / reference$loss_sum)
    } else {
      NA_real_
    }
  )

  list(
    metric = "log_loss",
    scope = "patient_level_out_of_fold",
    aggregation = "summed_case_loss_over_patients",
    definition = "negative_binary_log_likelihood",
    probability_clipping = "none",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    metric_direction = "lower_is_better",
    status = if (complete) "success" else "incomplete",
    estimate = if (complete) unname(loss_sum / n) else NA_real_,
    n = n,
    n_scored = n_scored,
    n_finite_loss = as.integer(sum(per_fold$n_finite_loss[evaluable])),
    n_infinite_loss = as.integer(sum(per_fold$n_infinite_loss[evaluable])),
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    loss_sum = loss_sum,
    per_fold = per_fold,
    reference = reference,
    skill_score = skill_score
  )
}

sb_binary_empty_calibration_fit <- function(fit_type, terms) {
  list(
    fit_type = fit_type,
    status = "not_fitted",
    coefficients = data.frame(
      term = terms,
      estimate = rep(NA_real_, length(terms)),
      stringsAsFactors = FALSE
    ),
    converged = NA,
    iterations = 0L,
    boundary = NA,
    deviance = NA_real_,
    warnings = character(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_binary_calibration_has_separation <- function(outcome, logit_probability) {
  if (length(unique(logit_probability)) < 2L) {
    return(FALSE)
  }
  positive <- logit_probability[outcome == 1L]
  negative <- logit_probability[outcome == 0L]
  min(positive) >= max(negative) || max(positive) <= min(negative)
}

sb_binary_fit_calibration <- function(outcome, logit_probability, fit_type) {
  intercept_only <- identical(fit_type, "calibration_in_the_large")
  terms <- if (intercept_only) {
    "(Intercept)"
  } else {
    c("(Intercept)", "logit_probability")
  }
  audit <- sb_binary_empty_calibration_fit(fit_type, terms)
  if (!intercept_only && sb_binary_calibration_has_separation(
    outcome,
    logit_probability
  )) {
    audit$status <- "failed"
    audit$error_class <- "stormbreaker_calibration_separation"
    audit$error_message <- paste0(
      "Complete or quasi-complete separation means that no finite ",
      "calibration slope maximum-likelihood estimate exists."
    )
    return(audit)
  }
  design <- if (intercept_only) {
    matrix(
      1,
      nrow = length(outcome),
      ncol = 1L,
      dimnames = list(NULL, "(Intercept)")
    )
  } else {
    cbind(
      "(Intercept)" = 1,
      logit_probability = logit_probability
    )
  }
  offset <- if (intercept_only) {
    logit_probability
  } else {
    rep(0, length(outcome))
  }
  captured_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm.fit(
        x = design,
        y = as.integer(outcome),
        family = stats::binomial(),
        offset = offset,
        control = stats::glm.control()
      ),
      warning = function(warning) {
        captured_warnings <<- c(
          captured_warnings,
          conditionMessage(warning)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  audit$warnings <- unique(captured_warnings)
  if (inherits(fit, "condition")) {
    audit$status <- "failed"
    audit$error_class <- class(fit)[[1L]]
    audit$error_message <- conditionMessage(fit)
    return(audit)
  }

  audit$coefficients$estimate <- unname(fit$coefficients[terms])
  audit$converged <- isTRUE(fit$converged)
  audit$iterations <- as.integer(fit$iter)
  audit$boundary <- isTRUE(fit$boundary)
  audit$deviance <- unname(fit$deviance)
  if (!audit$converged) {
    audit$status <- "failed"
    audit$error_class <-
      "stormbreaker_binary_calibration_nonconvergence"
    audit$error_message <- paste0(
      "The ", fit_type, " calibration fit did not converge."
    )
    return(audit)
  }
  if (any(!is.finite(audit$coefficients$estimate))) {
    audit$status <- "failed"
    audit$error_class <-
      "stormbreaker_binary_calibration_nonfinite_coefficients"
    audit$error_message <- paste0(
      "The ", fit_type,
      " calibration fit returned non-finite coefficients."
    )
    return(audit)
  }
  audit$status <- "success"
  audit
}

sb_binary_calibration_metric <- function(predictions) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  per_fold <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      c("outcome", "probability", "status"),
      drop = FALSE
    ]
    prediction_failed <- anyNA(fold_predictions$status) ||
      any(fold_predictions$status != "success") ||
      anyNA(fold_predictions$outcome) ||
      any(!fold_predictions$outcome %in% c(0L, 1L)) ||
      anyNA(fold_predictions$probability) ||
      any(!is.finite(fold_predictions$probability)) ||
      any(
        fold_predictions$probability < 0 |
          fold_predictions$probability > 1
      )
    data.frame(
      fold = fold_id,
      status = if (prediction_failed) {
        "prediction_failed"
      } else {
        "success"
      },
      n = as.integer(nrow(fold_predictions)),
      n_scored = if (prediction_failed) {
        0L
      } else {
        as.integer(nrow(fold_predictions))
      },
      n_positive = as.integer(sum(
        fold_predictions$outcome == 1L,
        na.rm = TRUE
      )),
      n_negative = as.integer(sum(
        fold_predictions$outcome == 0L,
        na.rm = TRUE
      )),
      n_boundary_probability = if (prediction_failed) {
        0L
      } else {
        as.integer(sum(fold_predictions$probability %in% c(0, 1)))
      },
      error_message = if (prediction_failed) {
        "OOF outcomes or probabilities for this fold are failed or invalid."
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  })
  per_fold <- do.call(rbind, per_fold)
  rownames(per_fold) <- NULL
  evaluable_folds <- per_fold$status == "success"
  complete_predictions <- length(fold_ids) > 0L && all(evaluable_folds)
  n <- as.integer(sum(per_fold$n))
  n_scored <- as.integer(sum(per_fold$n_scored[evaluable_folds]))
  n_positive <- as.integer(sum(per_fold$n_positive))
  n_negative <- as.integer(sum(per_fold$n_negative))
  n_boundary_probability <- as.integer(sum(
    per_fold$n_boundary_probability[evaluable_folds]
  ))
  intercept_fit <- sb_binary_empty_calibration_fit(
    "calibration_in_the_large",
    "(Intercept)"
  )
  slope_fit <- sb_binary_empty_calibration_fit(
    "logistic_recalibration",
    c("(Intercept)", "logit_probability")
  )
  result <- list(
    metric = "calibration_intercept_slope",
    scope = "pooled_patient_level_out_of_fold",
    purpose = "descriptive_evaluation_only",
    probability_clipping = "none",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = "incomplete",
    input_status = if (complete_predictions) {
      "success"
    } else {
      "prediction_failed"
    },
    calibration_intercept = NA_real_,
    calibration_slope = NA_real_,
    intercept_target = 0,
    slope_target = 1,
    n = n,
    n_scored = n_scored,
    n_positive = n_positive,
    n_negative = n_negative,
    n_boundary_probability = n_boundary_probability,
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable_folds)),
    per_fold = per_fold,
    intercept_fit = intercept_fit,
    slope_fit = slope_fit,
    error_class = NA_character_,
    error_message = NA_character_
  )
  if (!complete_predictions) {
    result$error_message <- paste0(
      "Complete successful OOF probabilities are required for calibration."
    )
    return(result)
  }
  if (n_positive == 0L || n_negative == 0L) {
    result$input_status <- "not_evaluable"
    result$error_message <- paste0(
      "Both outcome classes are required for calibration."
    )
    return(result)
  }
  if (n_boundary_probability > 0L) {
    result$input_status <- "not_evaluable"
    result$error_class <-
      "stormbreaker_binary_calibration_boundary_probability"
    result$error_message <- paste0(
      "Exact boundary probabilities cannot enter a finite-logit calibration ",
      "fit without probability clipping, which is not performed."
    )
    return(result)
  }

  logit_probability <- stats::qlogis(predictions$probability)
  intercept_fit <- sb_binary_fit_calibration(
    predictions$outcome,
    logit_probability,
    "calibration_in_the_large"
  )
  slope_fit <- sb_binary_fit_calibration(
    predictions$outcome,
    logit_probability,
    "logistic_recalibration"
  )
  result$intercept_fit <- intercept_fit
  result$slope_fit <- slope_fit
  if (identical(intercept_fit$status, "success")) {
    result$calibration_intercept <-
      intercept_fit$coefficients$estimate[[1L]]
  }
  if (identical(slope_fit$status, "success")) {
    result$calibration_slope <- slope_fit$coefficients$estimate[[2L]]
  }
  fits_complete <- identical(intercept_fit$status, "success") &&
    identical(slope_fit$status, "success")
  if (fits_complete) {
    result$status <- "success"
    return(result)
  }
  failed_fit <- if (!identical(intercept_fit$status, "success")) {
    intercept_fit
  } else {
    slope_fit
  }
  result$error_class <- failed_fit$error_class
  result$error_message <- paste(
    stats::na.omit(c(intercept_fit$error_message, slope_fit$error_message)),
    collapse = " "
  )
  result
}
