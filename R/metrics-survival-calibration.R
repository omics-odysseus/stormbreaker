# ---- metrics-survival-calibration ----

# ---- metrics-survival-brier ----

sb_fit_censoring_km <- function(time, event) {
  time <- sb_validate_endpoint_vector(time, "time", length(time))
  event <- sb_validate_endpoint_vector(event, "event", length(time))
  if (any(time < 0) || !all(event %in% c(0, 1))) {
    sb_abort(
      "The censoring estimate requires non-negative time and 0/1 event.",
      "stormbreaker_bad_censoring_input"
    )
  }

  times <- sort(unique(time))
  n_risk <- vapply(times, function(value) sum(time >= value), integer(1))
  n_event <- vapply(
    times,
    function(value) sum(time == value & event == 1L),
    integer(1)
  )
  n_censoring <- vapply(
    times,
    function(value) sum(time == value & event == 0L),
    integer(1)
  )
  n_risk_for_censoring <- n_risk - n_event
  survival_before <- numeric(length(times))
  survival_after <- numeric(length(times))
  current <- 1
  for (index in seq_along(times)) {
    survival_before[index] <- current
    if (n_censoring[index] > 0L) {
      current <- current * (
        1 - n_censoring[index] / n_risk_for_censoring[index]
      )
    }
    survival_after[index] <- current
  }

  data.frame(
    time = times,
    n_risk = n_risk,
    n_event = n_event,
    n_censoring = n_censoring,
    n_risk_for_censoring = n_risk_for_censoring,
    survival_before = survival_before,
    survival_after = survival_after,
    stringsAsFactors = FALSE
  )
}

sb_censoring_survival_at <- function(curve, time, left_limit = FALSE) {
  eligible <- if (left_limit) curve$time < time else curve$time <= time
  if (!any(eligible)) {
    return(1)
  }
  unname(utils::tail(curve$survival_after[eligible], 1L))
}

sb_empty_censoring_curve <- function() {
  data.frame(
    time = numeric(),
    n_risk = integer(),
    n_event = integer(),
    n_censoring = integer(),
    n_risk_for_censoring = integer(),
    survival_before = numeric(),
    survival_after = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_fit_event_km <- function(time, event) {
  time <- sb_validate_endpoint_vector(time, "time", length(time))
  event <- sb_validate_endpoint_vector(event, "event", length(time))
  if (any(time < 0) || !all(event %in% c(0, 1))) {
    sb_abort(
      "The survival reference requires non-negative time and 0/1 event.",
      "stormbreaker_bad_survival_reference_input"
    )
  }

  times <- sort(unique(time))
  n_risk <- vapply(times, function(value) sum(time >= value), integer(1))
  n_event <- vapply(
    times,
    function(value) sum(time == value & event == 1L),
    integer(1)
  )
  n_censoring <- vapply(
    times,
    function(value) sum(time == value & event == 0L),
    integer(1)
  )
  survival_before <- numeric(length(times))
  survival_after <- numeric(length(times))
  current <- 1
  for (index in seq_along(times)) {
    survival_before[index] <- current
    if (n_event[index] > 0L) {
      current <- current * (1 - n_event[index] / n_risk[index])
    }
    survival_after[index] <- current
  }

  data.frame(
    time = times,
    n_risk = n_risk,
    n_event = n_event,
    n_censoring = n_censoring,
    survival_before = survival_before,
    survival_after = survival_after,
    stringsAsFactors = FALSE
  )
}

sb_event_survival_at <- function(curve, time) {
  eligible <- curve$time <= time
  if (!any(eligible)) {
    return(1)
  }
  unname(utils::tail(curve$survival_after[eligible], 1L))
}

sb_empty_event_curve <- function() {
  data.frame(
    time = numeric(),
    n_risk = integer(),
    n_event = integer(),
    n_censoring = integer(),
    survival_before = numeric(),
    survival_after = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_survival_ipcw_brier_metric <- function(predictions, analysis_inputs,
                                          prediction_horizon) {
  prediction_horizon <- as.numeric(prediction_horizon)
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  fold_results <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      ,
      drop = FALSE
    ]
    audit <- data.frame(
      fold = fold_id,
      status = "analysis_input_failed",
      prediction_horizon = prediction_horizon,
      n = as.integer(nrow(fold_predictions)),
      n_scored = 0L,
      n_events_by_horizon = as.integer(sum(
        fold_predictions$time <= prediction_horizon &
          fold_predictions$event == 1L
      )),
      n_surviving_past_horizon = as.integer(sum(
        fold_predictions$time > prediction_horizon
      )),
      n_censored_at_horizon = as.integer(sum(
        fold_predictions$time == prediction_horizon &
          fold_predictions$event == 0L
      )),
      n_censored_by_horizon = as.integer(sum(
        fold_predictions$time <= prediction_horizon &
          fold_predictions$event == 0L
      )),
      censoring_survival_horizon = NA_real_,
      weight_sum = 0,
      weighted_loss_sum = 0,
      estimate = NA_real_,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    reference_audit <- data.frame(
      fold = fold_id,
      status = "analysis_input_failed",
      prediction_horizon = prediction_horizon,
      reference_survival_horizon = NA_real_,
      n = as.integer(nrow(fold_predictions)),
      n_scored = 0L,
      censoring_survival_horizon = NA_real_,
      weight_sum = 0,
      weighted_loss_sum = 0,
      estimate = NA_real_,
      source = "outer_analysis_event_only",
      ties = "events_before_censoring",
      horizon_evaluation = "right_continuous_after_events",
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    curve <- sb_empty_censoring_curve()
    reference_curve <- sb_empty_event_curve()
    input_folds <- vapply(
      analysis_inputs,
      function(input) as.integer(input$fold),
      integer(1)
    )
    matched <- which(input_folds == fold_id)
    if (length(matched) != 1L) {
      audit$error_class <- "stormbreaker_bad_ipcw_analysis_input"
      audit$error_message <- paste0(
        "Exactly one outer-analysis censoring input is required for this fold."
      )
      reference_audit$error_class <- audit$error_class
      reference_audit$error_message <- audit$error_message
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }

    analysis <- analysis_inputs[[matched]]
    curve_result <- tryCatch(
      sb_fit_censoring_km(analysis$time, analysis$event),
      error = function(error) error
    )
    if (inherits(curve_result, "error")) {
      audit$error_class <- class(curve_result)[1]
      audit$error_message <- conditionMessage(curve_result)
      reference_audit$error_class <- audit$error_class
      reference_audit$error_message <- audit$error_message
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }
    curve <- curve_result
    reference_curve_result <- tryCatch(
      sb_fit_event_km(analysis$time, analysis$event),
      error = function(error) error
    )
    if (inherits(reference_curve_result, "error")) {
      reference_audit$status <- "reference_failed"
      reference_audit$error_class <- class(reference_curve_result)[1]
      reference_audit$error_message <- conditionMessage(
        reference_curve_result
      )
    } else {
      reference_curve <- reference_curve_result
    }
    if (prediction_horizon > max(analysis$time)) {
      audit$status <- "horizon_unsupported"
      audit$error_class <- "stormbreaker_prediction_horizon_unsupported"
      audit$error_message <- paste0(
        "The prediction horizon exceeds this outer analysis fold's ",
        "observed follow-up."
      )
      reference_audit$status <- audit$status
      reference_audit$error_class <- audit$error_class
      reference_audit$error_message <- audit$error_message
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }

    censoring_survival_horizon <- sb_censoring_survival_at(
      curve,
      prediction_horizon,
      left_limit = TRUE
    )
    audit$censoring_survival_horizon <- censoring_survival_horizon
    reference_audit$censoring_survival_horizon <-
      censoring_survival_horizon
    if (!is.finite(censoring_survival_horizon) ||
        censoring_survival_horizon <= 0) {
      audit$status <- "censoring_failed"
      audit$error_class <- "stormbreaker_zero_censoring_survival"
      audit$error_message <- paste0(
        "The outer-analysis censoring survival is zero at the prediction ",
        "horizon."
      )
      reference_audit$status <- audit$status
      reference_audit$error_class <- audit$error_class
      reference_audit$error_message <- audit$error_message
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }

    event_by_horizon <- fold_predictions$time <= prediction_horizon &
      fold_predictions$event == 1L
    surviving_at_horizon <- fold_predictions$time > prediction_horizon |
      (
        fold_predictions$time == prediction_horizon &
          fold_predictions$event == 0L
      )
    event_survival <- vapply(
      fold_predictions$time[event_by_horizon],
      function(value) sb_censoring_survival_at(
        curve,
        value,
        left_limit = TRUE
      ),
      numeric(1)
    )
    if (any(!is.finite(event_survival)) || any(event_survival <= 0)) {
      audit$status <- "censoring_failed"
      audit$error_class <- "stormbreaker_zero_censoring_survival"
      audit$error_message <- paste0(
        "The outer-analysis censoring survival is zero before an observed ",
        "assessment event."
      )
      reference_audit$status <- audit$status
      reference_audit$error_class <- audit$error_class
      reference_audit$error_message <- audit$error_message
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }

    weights <- numeric(nrow(fold_predictions))
    weights[event_by_horizon] <- 1 / event_survival
    weights[surviving_at_horizon] <- 1 / censoring_survival_horizon
    if (anyNA(weights) || any(!is.finite(weights))) {
      audit$status <- "censoring_failed"
      audit$error_class <- "stormbreaker_invalid_ipcw_weight"
      audit$error_message <- "The IPCW weights are missing or non-finite."
      reference_audit$status <- audit$status
      reference_audit$error_class <- audit$error_class
      reference_audit$error_message <- audit$error_message
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }

    n_scored <- as.integer(sum(event_by_horizon | surviving_at_horizon))
    weight_sum <- unname(sum(weights))
    if (!inherits(reference_curve_result, "error")) {
      reference_survival <- sb_event_survival_at(
        reference_curve,
        prediction_horizon
      )
      valid_reference <- length(reference_survival) == 1L &&
        !is.na(reference_survival) &&
        is.finite(reference_survival) &&
        reference_survival >= 0 && reference_survival <= 1
      if (valid_reference) {
        reference_losses <- numeric(nrow(fold_predictions))
        reference_losses[event_by_horizon] <- reference_survival^2
        reference_losses[surviving_at_horizon] <- (1 - reference_survival)^2
        reference_weighted_loss <- weights * reference_losses
        if (anyNA(reference_weighted_loss) ||
            any(!is.finite(reference_weighted_loss))) {
          reference_audit$status <- "reference_failed"
          reference_audit$error_class <-
            "stormbreaker_invalid_reference_weighted_loss"
          reference_audit$error_message <- paste0(
            "The reference IPCW weighted losses are missing or non-finite."
          )
        } else {
          reference_audit$status <- "success"
          reference_audit$reference_survival_horizon <- reference_survival
          reference_audit$n_scored <- n_scored
          reference_audit$weight_sum <- weight_sum
          reference_audit$weighted_loss_sum <- unname(sum(
            reference_weighted_loss
          ))
          reference_audit$estimate <- unname(
            reference_audit$weighted_loss_sum / reference_audit$n
          )
        }
      } else {
        reference_audit$status <- "reference_failed"
        reference_audit$error_class <-
          "stormbreaker_invalid_reference_survival"
        reference_audit$error_message <- paste0(
          "The outer-analysis reference survival is invalid."
        )
      }
    }

    prediction_failed <-
      anyNA(fold_predictions$status) ||
      any(fold_predictions$status != "success") ||
      anyNA(fold_predictions$survival_probability_status) ||
      any(fold_predictions$survival_probability_status != "success") ||
      anyNA(fold_predictions$survival_probability) ||
      any(!is.finite(fold_predictions$survival_probability)) ||
      any(
        fold_predictions$survival_probability < 0 |
          fold_predictions$survival_probability > 1
      )
    if (prediction_failed) {
      audit$status <- "prediction_failed"
      audit$error_class <- "stormbreaker_invalid_survival_prediction"
      audit$error_message <- paste0(
        "OOF survival probabilities for this fold are failed or invalid."
      )
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }

    losses <- numeric(nrow(fold_predictions))
    losses[event_by_horizon] <-
      fold_predictions$survival_probability[event_by_horizon]^2
    losses[surviving_at_horizon] <- (
      1 - fold_predictions$survival_probability[surviving_at_horizon]
    )^2
    weighted_loss <- weights * losses
    if (anyNA(weighted_loss) || any(!is.finite(weighted_loss))) {
      audit$status <- "prediction_failed"
      audit$error_class <- "stormbreaker_invalid_ipcw_weighted_loss"
      audit$error_message <- paste0(
        "The model IPCW weighted losses are missing or non-finite."
      )
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }

    audit$status <- "success"
    audit$n_scored <- n_scored
    audit$weight_sum <- weight_sum
    audit$weighted_loss_sum <- unname(sum(weighted_loss))
    audit$estimate <- unname(audit$weighted_loss_sum / audit$n)
    list(
      audit = audit,
      curve = curve,
      reference_audit = reference_audit,
      reference_curve = reference_curve
    )
  })

  per_fold <- do.call(rbind, lapply(fold_results, `[[`, "audit"))
  rownames(per_fold) <- NULL
  censoring_curves <- lapply(fold_results, `[[`, "curve")
  names(censoring_curves) <- as.character(fold_ids)
  reference_per_fold <- do.call(
    rbind,
    lapply(fold_results, `[[`, "reference_audit")
  )
  rownames(reference_per_fold) <- NULL
  reference_curves <- lapply(fold_results, `[[`, "reference_curve")
  names(reference_curves) <- as.character(fold_ids)
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  n <- as.integer(sum(per_fold$n))
  weighted_loss_sum <- unname(sum(per_fold$weighted_loss_sum[evaluable]))
  reference_evaluable <- reference_per_fold$status == "success"
  reference_complete <- length(fold_ids) > 0L && all(reference_evaluable)
  reference_weighted_loss_sum <- unname(sum(
    reference_per_fold$weighted_loss_sum[reference_evaluable]
  ))
  reference <- list(
    metric = "ipcw_brier_reference",
    type = "outer_analysis_kaplan_meier",
    source = "outer_analysis_event_only",
    scope = "patient_level_out_of_fold_at_declared_horizon",
    censoring_model = "outer_analysis_reverse_kaplan_meier",
    aggregation = "summed_weighted_loss_over_all_assessment_patients",
    prediction_horizon = prediction_horizon,
    ties = "events_before_censoring",
    horizon_evaluation = "right_continuous_after_events",
    status = if (reference_complete) "success" else "incomplete",
    estimate = if (reference_complete) {
      reference_weighted_loss_sum / n
    } else {
      NA_real_
    },
    n = n,
    n_scored = as.integer(sum(
      reference_per_fold$n_scored[reference_evaluable]
    )),
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(reference_evaluable)),
    weight_sum = unname(sum(
      reference_per_fold$weight_sum[reference_evaluable]
    )),
    weighted_loss_sum = reference_weighted_loss_sum,
    per_fold = reference_per_fold,
    survival_curves = reference_curves
  )
  skill_complete <- complete && reference_complete &&
    is.finite(reference_weighted_loss_sum) &&
    reference_weighted_loss_sum > 0
  skill_score <- list(
    metric = "ipcw_brier_skill_score",
    reference = reference$type,
    formula = paste0(
      "1 - model_weighted_loss_sum / reference_weighted_loss_sum"
    ),
    status = if (skill_complete) "success" else "incomplete",
    estimate = if (skill_complete) {
      unname(1 - weighted_loss_sum / reference_weighted_loss_sum)
    } else {
      NA_real_
    }
  )

  list(
    metric = "ipcw_brier_score",
    scope = "patient_level_out_of_fold_at_declared_horizon",
    censoring_model = "outer_analysis_reverse_kaplan_meier",
    censoring_ties = "events_before_censoring",
    censoring_evaluation = "left_continuous_before_censoring_at_time",
    aggregation = "summed_weighted_loss_over_all_assessment_patients",
    prediction_horizon = prediction_horizon,
    status = if (complete) "success" else "incomplete",
    estimate = if (complete) weighted_loss_sum / n else NA_real_,
    n = n,
    n_scored = as.integer(sum(per_fold$n_scored[evaluable])),
    n_events_by_horizon = as.integer(sum(per_fold$n_events_by_horizon)),
    n_surviving_past_horizon = as.integer(sum(
      per_fold$n_surviving_past_horizon
    )),
    n_censored_at_horizon = as.integer(sum(
      per_fold$n_censored_at_horizon
    )),
    n_censored_by_horizon = as.integer(sum(
      per_fold$n_censored_by_horizon
    )),
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    weight_sum = unname(sum(per_fold$weight_sum[evaluable])),
    weighted_loss_sum = weighted_loss_sum,
    per_fold = per_fold,
    censoring_curves = censoring_curves,
    reference = reference,
    skill_score = skill_score
  )
}

# ---- metrics-survival-calibration ----

sb_empty_survival_calibration_knots <- function() {
  data.frame(
    role = character(),
    quantile_probability = numeric(),
    transformed_value = numeric(),
    survival_probability = numeric(),
    event_probability = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_empty_survival_calibration_patients <- function() {
  data.frame(
    row_instance = integer(),
    patient_id = character(),
    fold = integer(),
    time = numeric(),
    event = integer(),
    survival_probability = numeric(),
    predicted_event_probability = numeric(),
    transformed_prediction = numeric(),
    calibrated_survival_probability = numeric(),
    calibrated_event_probability = numeric(),
    absolute_error = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_empty_survival_calibration_curve <- function() {
  data.frame(
    survival_probability = numeric(),
    predicted_event_probability = numeric(),
    transformed_prediction = numeric(),
    calibrated_survival_probability = numeric(),
    calibrated_event_probability = numeric(),
    n_row_instances = integer(),
    stringsAsFactors = FALSE
  )
}

sb_empty_survival_calibration_folds <- function() {
  data.frame(
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
  )
}

sb_survival_calibration_empty_fit <- function() {
  list(
    status = "not_fitted",
    method = "cox_ph_three_knot_restricted_natural_cubic_spline",
    formula = paste0(
      "Surv(time, event) ~ calibration_ns1 + calibration_ns2"
    ),
    ties = "efron",
    baseline_reference = "coxph_mean_basis_centered",
    coefficients = data.frame(
      term = c("calibration_ns1", "calibration_ns2"),
      estimate = c(NA_real_, NA_real_),
      stringsAsFactors = FALSE
    ),
    means = c(calibration_ns1 = NA_real_, calibration_ns2 = NA_real_),
    baseline_cumulative_hazard_horizon = NA_real_,
    iterations = NA_integer_,
    log_likelihood = c(initial = NA_real_, final = NA_real_),
    warnings = character(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_survival_calibration_empty_result <- function(
    prediction_horizon,
    time_unit,
    scope,
    n = 0L) {
  safe_horizon <- if (is.numeric(prediction_horizon) &&
      length(prediction_horizon) >= 1L) {
    as.numeric(prediction_horizon[[1L]])
  } else {
    NA_real_
  }
  safe_unit <- if (is.character(time_unit) && length(time_unit) >= 1L) {
    as.character(time_unit[[1L]])
  } else {
    NA_character_
  }
  safe_scope <- if (is.character(scope) && length(scope) >= 1L) {
    as.character(scope[[1L]])
  } else {
    NA_character_
  }
  list(
    metric = "survival_calibration_curve",
    calibration_level = "moderate",
    estimand = "fixed_horizon_smoothed_calibration",
    method = "cox_ph_three_knot_restricted_natural_cubic_spline",
    scope = safe_scope,
    purpose = "descriptive_evaluation_only",
    prediction_horizon = safe_horizon,
    time_unit = safe_unit,
    prediction_scale = "survival_probability",
    event_probability_transform = "complementary_log_log",
    transform_formula = "log(-log(S(h)))",
    probability_clipping = "none",
    spline = list(
      type = "restricted_natural_cubic_spline",
      n_knots = 3L,
      knot_probabilities = c(0.10, 0.50, 0.90),
      quantile_type = 7L,
      outcome_tuning = "none"
    ),
    censoring_handling = "cox_partial_likelihood_complete_follow_up",
    censoring_assumption = "independent_noninformative_censoring",
    proportional_hazards_assumption = TRUE,
    administrative_censoring_at_horizon = FALSE,
    ties = "efron",
    baseline_reference = "coxph_mean_basis_centered",
    status = "incomplete",
    input_status = "not_checked",
    ici = NA_real_,
    e50 = NA_real_,
    e90 = NA_real_,
    n = as.integer(n),
    n_scored = 0L,
    n_events = 0L,
    n_censored = 0L,
    n_boundary_probability = 0L,
    n_followed_to_horizon = 0L,
    n_unique_survival_probability = 0L,
    n_folds = 0L,
    n_evaluable_folds = 0L,
    knots = sb_empty_survival_calibration_knots(),
    fit = sb_survival_calibration_empty_fit(),
    curve = sb_empty_survival_calibration_curve(),
    patient_calibration = sb_empty_survival_calibration_patients(),
    per_fold = sb_empty_survival_calibration_folds(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_survival_calibration_fail <- sb_binary_moderate_calibration_fail

sb_survival_calibration_error_metric <- function(
    error,
    prediction_horizon,
    time_unit,
    scope,
    n,
    fold = NULL) {
  result <- sb_survival_calibration_empty_result(
    prediction_horizon,
    time_unit,
    scope,
    n
  )
  if (!is.null(fold) && is.numeric(fold)) {
    fold_ids <- sort(unique(as.integer(fold[is.finite(fold)])))
    result$n_folds <- as.integer(length(fold_ids))
  }
  result$status <- "metric_failed"
  result$input_status <- "unexpected_error"
  result$error_class <- class(error)[[1L]]
  result$error_message <- conditionMessage(error)
  result
}

sb_survival_calibration_patient_table <- function(
    patient_id,
    time,
    event,
    survival_probability,
    fold,
    log_cumulative_hazard = NULL) {
  n <- length(time)
  safe_patient_id <- if (is.character(patient_id) &&
      length(patient_id) == n) {
    patient_id
  } else {
    rep(NA_character_, n)
  }
  safe_event <- if ((is.numeric(event) || is.logical(event)) &&
      length(event) == n) {
    as.integer(event)
  } else {
    rep(NA_integer_, n)
  }
  safe_probability <- if (is.numeric(survival_probability) &&
      length(survival_probability) == n) {
    as.numeric(survival_probability)
  } else {
    rep(NA_real_, n)
  }
  safe_fold <- if (is.null(fold)) {
    rep(NA_integer_, n)
  } else if (is.numeric(fold) && length(fold) == n) {
    as.integer(fold)
  } else {
    rep(NA_integer_, n)
  }
  valid_probability <- !is.na(safe_probability) &
    is.finite(safe_probability) &
    safe_probability >= 0 & safe_probability <= 1
  event_probability <- rep(NA_real_, n)
  transformed <- rep(NA_real_, n)
  event_probability[valid_probability] <- -expm1(
    log(safe_probability[valid_probability])
  )
  transformed[valid_probability] <- log(
    -log(safe_probability[valid_probability])
  )
  if (!is.null(log_cumulative_hazard)) {
    transformed <- sb_validate_stable_survival_scale(
      log_cumulative_hazard, safe_probability)
  }
  data.frame(
    row_instance = seq_len(n),
    patient_id = safe_patient_id,
    fold = safe_fold,
    time = if (is.numeric(time) && length(time) == n) {
      as.numeric(time)
    } else {
      rep(NA_real_, n)
    },
    event = safe_event,
    survival_probability = safe_probability,
    predicted_event_probability = event_probability,
    transformed_prediction = transformed,
    calibrated_survival_probability = rep(NA_real_, n),
    calibrated_event_probability = rep(NA_real_, n),
    absolute_error = rep(NA_real_, n),
    stringsAsFactors = FALSE
  )
}

sb_survival_calibration_fold_audit <- function(
    patients,
    prediction_status,
    prediction_horizon,
    fold_supplied) {
  if (!fold_supplied) {
    return(sb_empty_survival_calibration_folds())
  }
  fold_ids <- sort(unique(patients$fold))
  rows <- lapply(fold_ids, function(fold_id) {
    index <- patients$fold == fold_id
    valid_status <- is.character(prediction_status) &&
      length(prediction_status) == nrow(patients) &&
      !anyNA(prediction_status[index]) &&
      all(prediction_status[index] == "success")
    valid_rows <- valid_status &&
      !anyNA(patients$time[index]) &&
      all(is.finite(patients$time[index])) &&
      all(patients$time[index] >= 0) &&
      !anyNA(patients$event[index]) &&
      all(patients$event[index] %in% c(0L, 1L)) &&
      !anyNA(patients$survival_probability[index]) &&
      all(is.finite(patients$survival_probability[index])) &&
      all(
        patients$survival_probability[index] >= 0 &
          patients$survival_probability[index] <= 1
      )
    data.frame(
      fold = as.integer(fold_id),
      status = if (valid_rows) "success" else "prediction_failed",
      n = as.integer(sum(index)),
      n_scored = if (valid_rows) as.integer(sum(index)) else 0L,
      n_events = as.integer(sum(patients$event[index] == 1L, na.rm = TRUE)),
      n_censored = as.integer(sum(
        patients$event[index] == 0L,
        na.rm = TRUE
      )),
      n_boundary_probability = as.integer(sum(
        patients$survival_probability[index] %in% c(0, 1),
        na.rm = TRUE
      )),
      n_followed_to_horizon = as.integer(sum(
        patients$time[index] >= prediction_horizon,
        na.rm = TRUE
      )),
      error_class = if (valid_rows) {
        NA_character_
      } else {
        "stormbreaker_invalid_survival_calibration_prediction"
      },
      error_message = if (valid_rows) {
        NA_character_
      } else {
        "OOF survival endpoint or probability rows are failed or invalid."
      },
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

sb_survival_calibration_knots <- function(transformed_prediction) {
  values <- unname(stats::quantile(
    transformed_prediction,
    probs = c(0.10, 0.50, 0.90),
    names = FALSE,
    type = 7
  ))
  cumulative_hazard <- exp(values)
  data.frame(
    role = c("lower_boundary", "interior", "upper_boundary"),
    quantile_probability = c(0.10, 0.50, 0.90),
    transformed_value = values,
    survival_probability = exp(-cumulative_hazard),
    event_probability = -expm1(-cumulative_hazard),
    stringsAsFactors = FALSE
  )
}

sb_survival_calibration_basis <- function(transformed_prediction, knots) {
  basis <- splines::ns(
    transformed_prediction,
    knots = knots[[2L]],
    Boundary.knots = knots[c(1L, 3L)],
    intercept = FALSE
  )
  colnames(basis) <- c("calibration_ns1", "calibration_ns2")
  basis
}

sb_survival_calibration_fit_cox <- function(fit_data) {
  survival::coxph(
    survival::Surv(time, event) ~ calibration_ns1 + calibration_ns2,
    data = fit_data,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
}

sb_survival_probability_from_log_cumulative_hazard <- function(value) {
  exp(-exp(value))
}

sb_survival_calibration_curve_table <- function(patients) {
  survival_values <- sort(
    unique(patients$survival_probability),
    decreasing = TRUE
  )
  group <- match(patients$survival_probability, survival_values)
  first <- match(survival_values, patients$survival_probability)
  data.frame(
    survival_probability = survival_values,
    predicted_event_probability =
      patients$predicted_event_probability[first],
    transformed_prediction = patients$transformed_prediction[first],
    calibrated_survival_probability =
      patients$calibrated_survival_probability[first],
    calibrated_event_probability =
      patients$calibrated_event_probability[first],
    n_row_instances = as.integer(tabulate(
      group,
      nbins = length(survival_values)
    )),
    stringsAsFactors = FALSE
  )
}

sb_survival_calibration_metric <- function(
    patient_id,
    time,
    event,
    survival_probability,
    prediction_status,
    prediction_horizon,
    time_unit,
    scope,
    fold = NULL,
    log_cumulative_hazard = NULL) {
  n <- length(time)
  result <- sb_survival_calibration_empty_result(
    prediction_horizon,
    time_unit,
    scope,
    n
  )
  patients <- sb_survival_calibration_patient_table(
    patient_id,
    time,
    event,
    survival_probability,
    fold,
    log_cumulative_hazard
  )
  result$patient_calibration <- patients
  if (!is.null(log_cumulative_hazard)) result$transform_formula <- "saved_log_cumulative_hazard"

  valid_horizon <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
    prediction_horizon > 0
  valid_unit <- is.character(time_unit) && length(time_unit) == 1L &&
    !is.na(time_unit) && time_unit %in% c("day", "month", "year")
  valid_scope <- is.character(scope) && length(scope) == 1L &&
    !is.na(scope) && nzchar(scope)
  valid_patient_id <- is.character(patient_id) && length(patient_id) == n &&
    !anyNA(patient_id) && all(nzchar(trimws(patient_id)))
  valid_endpoint <- is.numeric(time) && length(time) == n &&
    (is.numeric(event) || is.logical(event)) && length(event) == n &&
    !anyNA(time) && !anyNA(event) &&
    all(is.finite(time)) && all(time >= 0) &&
    (!is.numeric(event) || all(is.finite(event))) &&
    all(event %in% c(0, 1))
  valid_probability <- is.numeric(survival_probability) &&
    length(survival_probability) == n &&
    !anyNA(survival_probability) &&
    all(is.finite(survival_probability)) &&
    all(survival_probability >= 0 & survival_probability <= 1)
  valid_status <- is.character(prediction_status) &&
    length(prediction_status) == n && !anyNA(prediction_status)
  fold_supplied <- !is.null(fold)
  valid_fold <- !fold_supplied || (
    is.numeric(fold) && length(fold) == n && !anyNA(fold) &&
      all(is.finite(fold)) && all(fold == trunc(fold))
  )

  if (valid_horizon) {
    result$n_followed_to_horizon <- as.integer(sum(
      patients$time >= prediction_horizon,
      na.rm = TRUE
    ))
  }
  if (valid_endpoint) {
    result$n_events <- as.integer(sum(as.integer(event) == 1L))
    result$n_censored <- as.integer(sum(as.integer(event) == 0L))
  }
  if (valid_probability) {
    result$n_boundary_probability <- as.integer(sum(
      survival_probability %in% c(0, 1)
    ))
    result$n_unique_survival_probability <- as.integer(length(unique(
      survival_probability
    )))
  }
  result$per_fold <- if (valid_fold && valid_horizon) {
    sb_survival_calibration_fold_audit(
      patients,
      prediction_status,
      prediction_horizon,
      fold_supplied
    )
  } else {
    sb_empty_survival_calibration_folds()
  }
  result$n_folds <- as.integer(nrow(result$per_fold))
  result$n_evaluable_folds <- as.integer(sum(
    result$per_fold$status == "success"
  ))

  if (!valid_horizon || !valid_unit || !valid_scope ||
      !valid_patient_id || !valid_endpoint || !valid_probability ||
      !valid_status || !valid_fold) {
    return(sb_survival_calibration_fail(
      result,
      "incomplete",
      "input_failed",
      "stormbreaker_bad_survival_calibration_input",
      paste0(
        "Survival calibration requires valid aligned endpoint, probability, ",
        "horizon, unit, scope, status, and optional fold vectors."
      )
    ))
  }
  if (any(prediction_status != "success")) {
    return(sb_survival_calibration_fail(
      result,
      "incomplete",
      "prediction_failed",
      "stormbreaker_invalid_survival_calibration_prediction",
      "Complete successful survival probabilities are required for calibration."
    ))
  }
  result$n_scored <- as.integer(n)
  if (result$n_events == 0L) {
    return(sb_survival_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_survival_calibration_no_event",
      "At least one observed event is required for survival calibration."
    ))
  }
  if (result$n_followed_to_horizon == 0L) {
    return(sb_survival_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_survival_calibration_horizon_unsupported",
      paste0(
        "The prediction horizon exceeds all observed endpoint follow-up."
      )
    ))
  }
  if (result$n_boundary_probability > 0L && is.null(log_cumulative_hazard)) {
    return(sb_survival_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_survival_calibration_boundary_probability",
      paste0(
        "Exact survival probabilities zero or one cannot enter a finite ",
        "complementary-log-log calibration fit without clipping, which is ",
        "not performed."
      )
    ))
  }

  transformed_prediction <- patients$transformed_prediction
  knots <- sb_survival_calibration_knots(transformed_prediction)
  result$knots <- knots
  knot_values <- knots$transformed_value
  if (anyNA(knot_values) || any(!is.finite(knot_values)) ||
      any(diff(knot_values) <= 0)) {
    return(sb_survival_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_survival_calibration_duplicate_knots",
      paste0(
        "The fixed 10th, 50th, and 90th percentile spline knots are not ",
        "finite and strictly increasing."
      )
    ))
  }

  basis_attempt <- tryCatch(
    sb_survival_calibration_basis(transformed_prediction, knot_values),
    error = function(error) error
  )
  if (inherits(basis_attempt, "condition")) {
    return(sb_survival_calibration_fail(
      result,
      "fit_failed",
      "success",
      class(basis_attempt)[[1L]],
      conditionMessage(basis_attempt)
    ))
  }
  basis <- basis_attempt
  fit_data <- data.frame(
    time = as.numeric(time),
    event = as.integer(event),
    calibration_ns1 = basis[, 1L],
    calibration_ns2 = basis[, 2L]
  )
  captured_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      sb_survival_calibration_fit_cox(fit_data),
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
  result$fit$warnings <- unique(captured_warnings)
  if (inherits(fit, "condition")) {
    result$fit$status <- "failed"
    result$fit$error_class <- class(fit)[[1L]]
    result$fit$error_message <- conditionMessage(fit)
    result$status <- "fit_failed"
    result$input_status <- "success"
    result$error_class <- result$fit$error_class
    result$error_message <- result$fit$error_message
    return(result)
  }

  coefficients <- stats::coef(fit)
  means <- fit$means
  result$fit$coefficients$estimate <- unname(coefficients[
    c("calibration_ns1", "calibration_ns2")
  ])
  result$fit$means <- unname(means[
    c("calibration_ns1", "calibration_ns2")
  ])
  names(result$fit$means) <- c("calibration_ns1", "calibration_ns2")
  result$fit$iterations <- as.integer(fit$iter[[1L]])
  if (length(fit$loglik) >= 2L) {
    result$fit$log_likelihood <- c(
      initial = unname(fit$loglik[[1L]]),
      final = unname(fit$loglik[[2L]])
    )
  }
  if (length(captured_warnings) > 0L) {
    result$fit$status <- "failed"
    result$fit$error_class <- "stormbreaker_survival_calibration_cox_warning"
    result$fit$error_message <- paste(unique(captured_warnings), collapse = " ")
    result$status <- "fit_failed"
    result$input_status <- "success"
    result$error_class <- result$fit$error_class
    result$error_message <- result$fit$error_message
    return(result)
  }
  valid_fit_parameters <- length(coefficients) == 2L &&
    all(c("calibration_ns1", "calibration_ns2") %in% names(coefficients)) &&
    !anyNA(coefficients) && all(is.finite(coefficients)) &&
    length(means) == 2L &&
    all(c("calibration_ns1", "calibration_ns2") %in% names(means)) &&
    !anyNA(means) && all(is.finite(means))
  if (!valid_fit_parameters) {
    result$fit$status <- "failed"
    result$fit$error_class <- "stormbreaker_survival_calibration_singular_fit"
    result$fit$error_message <- paste0(
      "The calibration Cox fit returned missing, singular, or non-finite ",
      "basis parameters."
    )
    result$status <- "fit_failed"
    result$input_status <- "success"
    result$error_class <- result$fit$error_class
    result$error_message <- result$fit$error_message
    return(result)
  }

  baseline <- tryCatch(
    withCallingHandlers(
      survival::basehaz(fit, centered = TRUE),
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
  result$fit$warnings <- unique(captured_warnings)
  if (inherits(baseline, "condition") || length(captured_warnings) > 0L) {
    result$fit$status <- "failed"
    result$fit$error_class <- if (inherits(baseline, "condition")) {
      class(baseline)[[1L]]
    } else {
      "stormbreaker_survival_calibration_baseline_warning"
    }
    result$fit$error_message <- if (inherits(baseline, "condition")) {
      conditionMessage(baseline)
    } else {
      paste(unique(captured_warnings), collapse = " ")
    }
    result$status <- "fit_failed"
    result$input_status <- "success"
    result$error_class <- result$fit$error_class
    result$error_message <- result$fit$error_message
    return(result)
  }
  valid_baseline <- is.data.frame(baseline) &&
    all(c("hazard", "time") %in% names(baseline)) &&
    !anyNA(baseline$hazard) && !anyNA(baseline$time) &&
    all(is.finite(baseline$hazard)) && all(baseline$hazard >= 0) &&
    all(is.finite(baseline$time))
  if (!valid_baseline) {
    result$fit$status <- "failed"
    result$fit$error_class <- "stormbreaker_survival_calibration_bad_baseline"
    result$fit$error_message <- paste0(
      "The centered calibration baseline cumulative hazard is invalid."
    )
    result$status <- "fit_failed"
    result$input_status <- "success"
    result$error_class <- result$fit$error_class
    result$error_message <- result$fit$error_message
    return(result)
  }
  eligible <- baseline$time <= prediction_horizon
  baseline_horizon <- if (any(eligible)) {
    unname(utils::tail(baseline$hazard[eligible], 1L))
  } else {
    0
  }
  result$fit$baseline_cumulative_hazard_horizon <- baseline_horizon

  centered_basis <- sweep(
    basis,
    2L,
    means[c("calibration_ns1", "calibration_ns2")],
    FUN = "-"
  )
  centered_lp <- drop(centered_basis %*% coefficients[
    c("calibration_ns1", "calibration_ns2")
  ])
  if (baseline_horizon == 0) {
    calibrated_survival <- rep(1, n)
    calibrated_event <- rep(0, n)
  } else {
    log_cumulative_hazard <- log(baseline_horizon) + centered_lp
    calibrated_survival <-
      sb_survival_probability_from_log_cumulative_hazard(
        log_cumulative_hazard
      )
    calibrated_event <- -expm1(-exp(log_cumulative_hazard))
  }
  valid_calibrated <- !anyNA(calibrated_survival) &&
    !anyNA(calibrated_event) &&
    all(is.finite(calibrated_survival)) &&
    all(is.finite(calibrated_event)) &&
    all(calibrated_survival >= 0 & calibrated_survival <= 1) &&
    all(calibrated_event >= 0 & calibrated_event <= 1)
  if (!valid_calibrated) {
    result$fit$status <- "failed"
    result$fit$error_class <-
      "stormbreaker_survival_calibration_invalid_probability"
    result$fit$error_message <- paste0(
      "The calibration fit produced invalid horizon probabilities."
    )
    result$status <- "fit_failed"
    result$input_status <- "success"
    result$error_class <- result$fit$error_class
    result$error_message <- result$fit$error_message
    return(result)
  }

  absolute_error <- abs(calibrated_survival - survival_probability)
  result$patient_calibration$calibrated_survival_probability <-
    calibrated_survival
  result$patient_calibration$calibrated_event_probability <- calibrated_event
  result$patient_calibration$absolute_error <- absolute_error
  result$curve <- sb_survival_calibration_curve_table(
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

# ---- metrics-survival-decision-curve ----

sb_survival_dca_horizon_km <- function(time, event, prediction_horizon) {
  if (length(time) == 0L) {
    return(list(
      status = "success",
      horizon_support = "empty_selection",
      survival_probability = 1,
      event_probability = 0,
      max_follow_up = NA_real_,
      n_followed_to_horizon = 0L
    ))
  }

  curve <- sb_fit_event_km(time, event)
  survival_probability <- sb_event_survival_at(
    curve,
    prediction_horizon
  )
  max_follow_up <- max(time)
  observed_to_horizon <- max_follow_up >= prediction_horizon
  zero_before_horizon <- !observed_to_horizon &&
    isTRUE(survival_probability == 0)

  list(
    status = if (observed_to_horizon || zero_before_horizon) {
      "success"
    } else {
      "not_evaluable"
    },
    horizon_support = if (observed_to_horizon) {
      "observed_to_horizon"
    } else if (zero_before_horizon) {
      "survival_zero_before_horizon"
    } else {
      "insufficient_follow_up"
    },
    survival_probability = survival_probability,
    event_probability = 1 - survival_probability,
    max_follow_up = unname(max_follow_up),
    n_followed_to_horizon = as.integer(sum(time >= prediction_horizon))
  )
}

sb_survival_dca_empty_curve <- function(decision_thresholds, n) {
  threshold_odds <- decision_thresholds / (1 - decision_thresholds)
  data.frame(
    threshold_probability = decision_thresholds,
    threshold_odds = threshold_odds,
    n = rep.int(as.integer(n), length(decision_thresholds)),
    n_selected = rep.int(NA_integer_, length(decision_thresholds)),
    selected_fraction = rep.int(NA_real_, length(decision_thresholds)),
    n_selected_events_by_horizon = rep.int(
      NA_integer_,
      length(decision_thresholds)
    ),
    n_selected_followed_to_horizon = rep.int(
      NA_integer_,
      length(decision_thresholds)
    ),
    max_selected_follow_up = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    horizon_support = rep.int(NA_character_, length(decision_thresholds)),
    selected_survival_probability = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    selected_event_probability = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    estimated_true_positive_rate = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    estimated_false_positive_rate = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    model_net_benefit = rep.int(NA_real_, length(decision_thresholds)),
    treat_all_net_benefit = rep.int(NA_real_, length(decision_thresholds)),
    treat_none_net_benefit = rep.int(0, length(decision_thresholds)),
    status = rep.int("input_failed", length(decision_thresholds)),
    error_class = rep.int(NA_character_, length(decision_thresholds)),
    error_message = rep.int(NA_character_, length(decision_thresholds)),
    stringsAsFactors = FALSE
  )
}

sb_survival_dca_valid_input <- function(predictions, prediction_horizon,
                                        time_unit, scope) {
  required <- c("time", "event", "survival_probability", "status")
  valid_frame <- is.data.frame(predictions) && nrow(predictions) > 0L &&
    all(vapply(required, function(name) {
      sum(names(predictions) == name) == 1L
    }, logical(1)))
  if (!valid_frame) {
    return(FALSE)
  }

  n <- nrow(predictions)
  time <- predictions$time
  event <- predictions$event
  probability <- predictions$survival_probability
  status <- predictions$status
  valid_probability_status <- TRUE
  if ("survival_probability_status" %in% names(predictions)) {
    valid_probability_status <-
      sum(names(predictions) == "survival_probability_status") == 1L &&
      is.character(predictions$survival_probability_status) &&
      length(predictions$survival_probability_status) == n &&
      !anyNA(predictions$survival_probability_status) &&
      all(predictions$survival_probability_status == "success")
  }

  is.numeric(time) && length(time) == n && !anyNA(time) &&
    all(is.finite(time)) && all(time >= 0) &&
    (is.numeric(event) || is.logical(event)) && length(event) == n &&
    !anyNA(event) && (!is.numeric(event) || all(is.finite(event))) &&
    all(event %in% c(0, 1)) &&
    is.numeric(probability) && length(probability) == n &&
    !anyNA(probability) && all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1) &&
    is.character(status) && length(status) == n && !anyNA(status) &&
    all(status == "success") && valid_probability_status &&
    is.numeric(prediction_horizon) && length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
    prediction_horizon > 0 &&
    is.character(time_unit) && length(time_unit) == 1L &&
    !is.na(time_unit) && time_unit %in% c("day", "month", "year") &&
    is.character(scope) && length(scope) == 1L && !is.na(scope) &&
    nzchar(scope)
}

sb_survival_decision_curve_metric <- function(
    predictions,
    decision_thresholds,
    prediction_horizon,
    time_unit,
    scope) {
  decision_thresholds <- sb_validate_decision_thresholds(
    decision_thresholds
  )
  n <- if (is.data.frame(predictions)) nrow(predictions) else 0L
  curve <- sb_survival_dca_empty_curve(decision_thresholds, n)
  result <- list(
    metric = "decision_curve",
    task_type = "survival",
    scope = scope,
    purpose = "descriptive_clinical_utility_evaluation",
    prediction_target = "event_by_declared_horizon",
    threshold_source = "study_declared_survival_task",
    decision_rule = paste0(
      "one_minus_survival_probability_greater_than_or_equal_to_threshold"
    ),
    event_probability_estimator =
      "kaplan_meier_within_threshold_selected_group",
    net_benefit_formula = paste0(
      "selected_fraction * event_probability - selected_fraction * ",
      "survival_probability * threshold / (1 - threshold)"
    ),
    reference_strategies = c("treat_all", "treat_none"),
    censoring_assumption = paste0(
      "independent_noninformative_censoring_within_relevant_risk_groups"
    ),
    censoring_predictor_assumption =
      "censoring_independent_of_model_predictors",
    survival_ties = "events_before_censoring",
    prediction_horizon = as.numeric(prediction_horizon)[[1L]],
    time_unit = as.character(time_unit)[[1L]],
    status = "incomplete",
    n = as.integer(n),
    n_scored = 0L,
    n_thresholds = as.integer(length(decision_thresholds)),
    n_evaluable_thresholds = 0L,
    overall_horizon_support = NA_character_,
    overall_max_follow_up = NA_real_,
    overall_survival_probability = NA_real_,
    overall_event_probability = NA_real_,
    curve = curve,
    error_class = "stormbreaker_invalid_survival_prediction",
    error_message = paste0(
      "Survival endpoints, locked-horizon probabilities, or prediction ",
      "statuses are incomplete or invalid."
    )
  )

  if (!sb_survival_dca_valid_input(
    predictions,
    prediction_horizon,
    time_unit,
    scope
  )) {
    result$curve$status <- "prediction_failed"
    result$curve$error_class <- result$error_class
    result$curve$error_message <- result$error_message
    return(result)
  }

  time <- predictions$time
  event <- as.integer(predictions$event)
  event_probability <- 1 - predictions$survival_probability
  overall <- sb_survival_dca_horizon_km(
    time,
    event,
    prediction_horizon
  )
  result$n_scored <- as.integer(n)
  result$overall_horizon_support <- overall$horizon_support
  result$overall_max_follow_up <- overall$max_follow_up
  result$overall_survival_probability <- overall$survival_probability
  result$overall_event_probability <- overall$event_probability

  rows <- lapply(seq_along(decision_thresholds), function(index) {
    threshold <- decision_thresholds[[index]]
    selected <- event_probability >= threshold
    selected_n <- as.integer(sum(selected))
    selected_fraction <- selected_n / n
    selected_km <- sb_survival_dca_horizon_km(
      time[selected],
      event[selected],
      prediction_horizon
    )
    true_positive_rate <- if (identical(selected_km$status, "success")) {
      selected_fraction * selected_km$event_probability
    } else {
      NA_real_
    }
    false_positive_rate <- if (identical(selected_km$status, "success")) {
      selected_fraction * selected_km$survival_probability
    } else {
      NA_real_
    }
    model_net_benefit <- if (identical(selected_km$status, "success")) {
      true_positive_rate - false_positive_rate *
        curve$threshold_odds[[index]]
    } else {
      NA_real_
    }
    treat_all_net_benefit <- if (identical(overall$status, "success")) {
      overall$event_probability - overall$survival_probability *
        curve$threshold_odds[[index]]
    } else {
      NA_real_
    }
    row_success <- identical(selected_km$status, "success") &&
      identical(overall$status, "success")
    failed_support <- if (!identical(selected_km$status, "success")) {
      "The threshold-selected group has insufficient follow-up at the horizon."
    } else {
      "The complete evaluation sample has insufficient follow-up at the horizon."
    }

    data.frame(
      threshold_probability = threshold,
      threshold_odds = curve$threshold_odds[[index]],
      n = as.integer(n),
      n_selected = selected_n,
      selected_fraction = selected_fraction,
      n_selected_events_by_horizon = as.integer(sum(
        selected & time <= prediction_horizon & event == 1L
      )),
      n_selected_followed_to_horizon =
        selected_km$n_followed_to_horizon,
      max_selected_follow_up = selected_km$max_follow_up,
      horizon_support = selected_km$horizon_support,
      selected_survival_probability =
        selected_km$survival_probability,
      selected_event_probability = selected_km$event_probability,
      estimated_true_positive_rate = true_positive_rate,
      estimated_false_positive_rate = false_positive_rate,
      model_net_benefit = model_net_benefit,
      treat_all_net_benefit = treat_all_net_benefit,
      treat_none_net_benefit = 0,
      status = if (row_success) "success" else "not_evaluable",
      error_class = if (row_success) {
        NA_character_
      } else {
        "stormbreaker_insufficient_survival_dca_follow_up"
      },
      error_message = if (row_success) NA_character_ else failed_support,
      stringsAsFactors = FALSE
    )
  })
  result$curve <- do.call(rbind, rows)
  rownames(result$curve) <- NULL
  result$n_evaluable_thresholds <- as.integer(sum(
    result$curve$status == "success"
  ))
  complete <- result$n_evaluable_thresholds == result$n_thresholds
  result$status <- if (complete) "success" else "incomplete"
  result$error_class <- if (complete) {
    NA_character_
  } else {
    "stormbreaker_insufficient_survival_dca_follow_up"
  }
  result$error_message <- if (complete) {
    NA_character_
  } else {
    "At least one declared survival decision threshold is not evaluable."
  }
  result
}

# ---- metrics-survival-weak-calibration ----

sb_empty_survival_weak_calibration_patients <- function() {
  data.frame(
    row_instance = integer(),
    patient_id = character(),
    fold = integer(),
    time = numeric(),
    event = integer(),
    survival_probability = numeric(),
    predicted_event_probability = numeric(),
    transformed_prediction = numeric(),
    administrative_time = numeric(),
    event_by_horizon = integer(),
    stringsAsFactors = FALSE
  )
}

sb_empty_survival_weak_calibration_folds <- function() {
  sb_empty_survival_calibration_folds()
}

sb_survival_weak_calibration_empty_mean <- function() {
  list(
    status = "not_evaluated",
    method = "kaplan_meier_observed_over_mean_prediction",
    ratio_direction = "observed_over_expected",
    horizon_support = "not_checked",
    observed_event_probability = NA_real_,
    mean_predicted_event_probability = NA_real_,
    mean_calibration_ratio = NA_real_,
    n_followed_to_horizon = 0L,
    maximum_follow_up = NA_real_,
    km_survival_probability = NA_real_,
    event_curve = sb_empty_event_curve(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_survival_weak_calibration_empty_slope <- function() {
  list(
    status = "not_fitted",
    method = "secondary_cox_proportional_hazards",
    formula = paste0(
      "Surv(administrative_time, event_by_horizon) ~ ",
      "transformed_prediction"
    ),
    ties = "efron",
    administrative_censoring_at_horizon = TRUE,
    proportional_hazards_assumption = TRUE,
    censoring_assumption = "independent_noninformative_censoring",
    coefficient = NA_real_,
    n_events_by_horizon = 0L,
    n_unique_transformed_prediction = 0L,
    iterations = NA_integer_,
    log_likelihood = c(initial = NA_real_, final = NA_real_),
    warnings = character(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_survival_weak_calibration_empty_result <- function(
    prediction_horizon,
    time_unit,
    scope,
    n = 0L) {
  safe_horizon <- if (is.numeric(prediction_horizon) &&
      length(prediction_horizon) >= 1L) {
    as.numeric(prediction_horizon[[1L]])
  } else {
    NA_real_
  }
  safe_unit <- if (is.character(time_unit) && length(time_unit) >= 1L) {
    as.character(time_unit[[1L]])
  } else {
    NA_character_
  }
  safe_scope <- if (is.character(scope) && length(scope) >= 1L) {
    as.character(scope[[1L]])
  } else {
    NA_character_
  }
  list(
    metric = "survival_weak_calibration",
    calibration_level = "weak",
    estimand = "fixed_horizon_mean_calibration_and_slope",
    scope = safe_scope,
    purpose = "descriptive_evaluation_only",
    prediction_horizon = safe_horizon,
    time_unit = safe_unit,
    prediction_scale = "survival_probability",
    event_probability_transform = "complementary_log_log",
    transform_formula = "log(-log(S(h)))",
    probability_clipping = "none",
    status = "incomplete",
    input_status = "not_checked",
    mean_calibration_ratio = NA_real_,
    calibration_slope = NA_real_,
    mean_calibration_target = 1,
    slope_target = 1,
    n = as.integer(n),
    n_scored = 0L,
    n_events = 0L,
    n_censored = 0L,
    n_events_by_horizon = 0L,
    n_boundary_probability = 0L,
    n_followed_to_horizon = 0L,
    n_folds = 0L,
    n_evaluable_folds = 0L,
    mean_calibration = sb_survival_weak_calibration_empty_mean(),
    slope_fit = sb_survival_weak_calibration_empty_slope(),
    patient_predictions = sb_empty_survival_weak_calibration_patients(),
    per_fold = sb_empty_survival_weak_calibration_folds(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_survival_weak_calibration_error_metric <- function(
    error,
    prediction_horizon,
    time_unit,
    scope,
    n,
    fold = NULL) {
  result <- sb_survival_weak_calibration_empty_result(
    prediction_horizon,
    time_unit,
    scope,
    n
  )
  if (!is.null(fold) && is.numeric(fold)) {
    fold_ids <- sort(unique(as.integer(fold[is.finite(fold)])))
    result$n_folds <- as.integer(length(fold_ids))
  }
  result$status <- "metric_failed"
  result$input_status <- "unexpected_error"
  result$error_class <- class(error)[[1L]]
  result$error_message <- conditionMessage(error)
  result
}

sb_survival_weak_calibration_patient_table <- function(
    patient_id,
    time,
    event,
    survival_probability,
    fold,
    prediction_horizon,
    log_cumulative_hazard = NULL) {
  patients <- sb_survival_calibration_patient_table(
    patient_id,
    time,
    event,
    survival_probability,
    fold,
    log_cumulative_hazard
  )
  n <- nrow(patients)
  valid_horizon <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon)
  valid_time <- is.numeric(time) && length(time) == n
  valid_event <- (is.numeric(event) || is.logical(event)) &&
    length(event) == n
  patients$administrative_time <- if (valid_horizon && valid_time) {
    pmin(as.numeric(time), as.numeric(prediction_horizon))
  } else {
    rep(NA_real_, n)
  }
  patients$event_by_horizon <- if (
      valid_horizon && valid_time && valid_event) {
    as.integer(as.integer(event) == 1L & time <= prediction_horizon)
  } else {
    rep(NA_integer_, n)
  }
  patients[c(
    "row_instance",
    "patient_id",
    "fold",
    "time",
    "event",
    "survival_probability",
    "predicted_event_probability",
    "transformed_prediction",
    "administrative_time",
    "event_by_horizon"
  )]
}

sb_survival_weak_calibration_fit_slope <- function(fit_data) {
  survival::coxph(
    survival::Surv(administrative_time, event_by_horizon) ~
      transformed_prediction,
    data = fit_data,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
}

sb_survival_weak_calibration_component_error <- function(
    component,
    status,
    error_class,
    error_message) {
  component$status <- status
  component$error_class <- error_class
  component$error_message <- error_message
  component
}

sb_survival_weak_calibration_finish <- function(result) {
  if (identical(result$mean_calibration$status, "success") &&
      identical(result$slope_fit$status, "success")) {
    result$status <- "success"
    return(result)
  }
  failed_component <- if (!identical(
      result$mean_calibration$status,
      "success"
    )) {
    result$mean_calibration
  } else {
    result$slope_fit
  }
  result$status <- "incomplete"
  result$error_class <- failed_component$error_class
  result$error_message <- failed_component$error_message
  result
}

sb_survival_weak_calibration_metric <- function(
    patient_id,
    time,
    event,
    survival_probability,
    prediction_status,
    prediction_horizon,
    time_unit,
    scope,
    fold = NULL,
    log_cumulative_hazard = NULL) {
  n <- length(time)
  result <- sb_survival_weak_calibration_empty_result(
    prediction_horizon,
    time_unit,
    scope,
    n
  )
  patients <- sb_survival_weak_calibration_patient_table(
    patient_id,
    time,
    event,
    survival_probability,
    fold,
    prediction_horizon,
    log_cumulative_hazard
  )
  result$patient_predictions <- patients
  if (!is.null(log_cumulative_hazard)) result$transform_formula <- "saved_log_cumulative_hazard"

  valid_horizon <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
    prediction_horizon > 0
  valid_unit <- is.character(time_unit) && length(time_unit) == 1L &&
    !is.na(time_unit) && time_unit %in% c("day", "month", "year")
  valid_scope <- is.character(scope) && length(scope) == 1L &&
    !is.na(scope) && nzchar(scope)
  valid_patient_id <- n > 0L && is.character(patient_id) &&
    length(patient_id) == n && !anyNA(patient_id) &&
    all(nzchar(trimws(patient_id)))
  valid_endpoint <- is.numeric(time) && length(time) == n &&
    (is.numeric(event) || is.logical(event)) && length(event) == n &&
    !anyNA(time) && !anyNA(event) &&
    all(is.finite(time)) && all(time >= 0) &&
    (!is.numeric(event) || all(is.finite(event))) &&
    all(event %in% c(0, 1))
  valid_probability <- is.numeric(survival_probability) &&
    length(survival_probability) == n &&
    !anyNA(survival_probability) &&
    all(is.finite(survival_probability)) &&
    all(survival_probability >= 0 & survival_probability <= 1)
  valid_status <- is.character(prediction_status) &&
    length(prediction_status) == n && !anyNA(prediction_status)
  fold_supplied <- !is.null(fold)
  valid_fold <- !fold_supplied || (
    is.numeric(fold) && length(fold) == n && !anyNA(fold) &&
      all(is.finite(fold)) && all(fold == trunc(fold))
  )

  if (valid_horizon && valid_endpoint) {
    result$n_followed_to_horizon <- as.integer(sum(time >= prediction_horizon))
    result$n_events_by_horizon <- as.integer(sum(
      event == 1L & time <= prediction_horizon
    ))
  }
  if (valid_endpoint) {
    result$n_events <- as.integer(sum(event == 1L))
    result$n_censored <- as.integer(sum(event == 0L))
  }
  if (valid_probability) {
    result$n_boundary_probability <- as.integer(sum(
      survival_probability %in% c(0, 1)
    ))
  }
  result$per_fold <- if (valid_fold && valid_horizon) {
    sb_survival_calibration_fold_audit(
      patients,
      prediction_status,
      prediction_horizon,
      fold_supplied
    )
  } else {
    sb_empty_survival_weak_calibration_folds()
  }
  result$n_folds <- as.integer(nrow(result$per_fold))
  result$n_evaluable_folds <- as.integer(sum(
    result$per_fold$status == "success"
  ))

  if (!valid_horizon || !valid_unit || !valid_scope ||
      !valid_patient_id || !valid_endpoint || !valid_probability ||
      !valid_status || !valid_fold) {
    result$input_status <- "input_failed"
    result$error_class <- "stormbreaker_bad_survival_weak_calibration_input"
    result$error_message <- paste0(
      "Survival weak calibration requires valid aligned endpoint, ",
      "probability, horizon, unit, scope, status, and optional fold vectors."
    )
    return(result)
  }
  if (any(prediction_status != "success")) {
    result$input_status <- "prediction_failed"
    result$error_class <-
      "stormbreaker_invalid_survival_weak_calibration_prediction"
    result$error_message <- paste0(
      "Complete successful survival probabilities are required for weak ",
      "calibration."
    )
    return(result)
  }

  result$input_status <- "success"
  result$n_scored <- as.integer(n)

  event_curve <- sb_fit_event_km(time, event)
  maximum_follow_up <- max(time)
  terminal_survival <- utils::tail(event_curve$survival_after, 1L)
  horizon_support <- if (maximum_follow_up >= prediction_horizon) {
    "observed_to_horizon"
  } else if (terminal_survival == 0) {
    "survival_zero_before_horizon"
  } else {
    "unsupported"
  }
  result$mean_calibration$event_curve <- event_curve
  result$mean_calibration$maximum_follow_up <- maximum_follow_up
  result$mean_calibration$n_followed_to_horizon <-
    result$n_followed_to_horizon
  result$mean_calibration$horizon_support <- horizon_support
  if (horizon_support == "unsupported") {
    result$mean_calibration <- sb_survival_weak_calibration_component_error(
      result$mean_calibration,
      "not_evaluable",
      "stormbreaker_survival_weak_calibration_horizon_unsupported",
      "The prediction horizon exceeds all observed endpoint follow-up."
    )
  } else {
    km_survival <- sb_event_survival_at(event_curve, prediction_horizon)
    observed_event_probability <- -expm1(log(km_survival))
    mean_predicted_event_probability <- mean(
      patients$predicted_event_probability
    )
    result$mean_calibration$km_survival_probability <- km_survival
    result$mean_calibration$observed_event_probability <-
      observed_event_probability
    result$mean_calibration$mean_predicted_event_probability <-
      mean_predicted_event_probability
    if (mean_predicted_event_probability == 0) {
      result$mean_calibration <- sb_survival_weak_calibration_component_error(
        result$mean_calibration,
        "not_evaluable",
        "stormbreaker_survival_weak_calibration_zero_mean_prediction",
        "The mean predicted event probability is zero."
      )
    } else {
      ratio <- observed_event_probability / mean_predicted_event_probability
      if (is.finite(ratio)) {
        result$mean_calibration$status <- "success"
        result$mean_calibration$mean_calibration_ratio <- unname(ratio)
        result$mean_calibration_ratio <- unname(ratio)
      } else {
        result$mean_calibration <- sb_survival_weak_calibration_component_error(
          result$mean_calibration,
          "not_evaluable",
          "stormbreaker_survival_weak_calibration_nonfinite_ratio",
          "The observed-to-expected calibration ratio is not finite."
        )
      }
    }
  }

  transformed_prediction <- patients$transformed_prediction
  result$slope_fit$n_events_by_horizon <- result$n_events_by_horizon
  result$slope_fit$n_unique_transformed_prediction <- as.integer(length(
    unique(transformed_prediction[is.finite(transformed_prediction)])
  ))
  if (result$n_boundary_probability > 0L && is.null(log_cumulative_hazard)) {
    result$slope_fit <- sb_survival_weak_calibration_component_error(
      result$slope_fit,
      "not_evaluable",
      "stormbreaker_survival_weak_calibration_boundary_probability",
      paste0(
        "Exact survival probabilities zero or one do not have a finite ",
        "complementary-log-log value; probabilities are not clipped."
      )
    )
    return(sb_survival_weak_calibration_finish(result))
  }
  if (result$n_events_by_horizon == 0L) {
    result$slope_fit <- sb_survival_weak_calibration_component_error(
      result$slope_fit,
      "not_evaluable",
      "stormbreaker_survival_weak_calibration_no_event",
      "At least one observed event by the prediction horizon is required."
    )
    return(sb_survival_weak_calibration_finish(result))
  }
  if (result$slope_fit$n_unique_transformed_prediction < 2L) {
    result$slope_fit <- sb_survival_weak_calibration_component_error(
      result$slope_fit,
      "not_evaluable",
      "stormbreaker_survival_weak_calibration_constant_prediction",
      "The transformed survival prediction is constant."
    )
    return(sb_survival_weak_calibration_finish(result))
  }

  fit_data <- data.frame(
    administrative_time = patients$administrative_time,
    event_by_horizon = patients$event_by_horizon,
    transformed_prediction = transformed_prediction
  )
  captured_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      sb_survival_weak_calibration_fit_slope(fit_data),
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
  result$slope_fit$warnings <- unique(captured_warnings)
  if (inherits(fit, "condition")) {
    result$slope_fit <- sb_survival_weak_calibration_component_error(
      result$slope_fit,
      "failed",
      class(fit)[[1L]],
      conditionMessage(fit)
    )
    return(sb_survival_weak_calibration_finish(result))
  }
  if (length(captured_warnings) > 0L) {
    result$slope_fit <- sb_survival_weak_calibration_component_error(
      result$slope_fit,
      "failed",
      "stormbreaker_survival_weak_calibration_cox_warning",
      paste(unique(captured_warnings), collapse = " ")
    )
    return(sb_survival_weak_calibration_finish(result))
  }

  coefficients <- stats::coef(fit)
  coefficient <- if (
      "transformed_prediction" %in% names(coefficients)) {
    unname(coefficients[["transformed_prediction"]])
  } else {
    NA_real_
  }
  result$slope_fit$iterations <- if (!is.null(fit$iter) &&
      length(fit$iter) >= 1L) {
    as.integer(fit$iter[[1L]])
  } else {
    NA_integer_
  }
  if (!is.null(fit$loglik) && length(fit$loglik) >= 2L) {
    result$slope_fit$log_likelihood <- c(
      initial = unname(fit$loglik[[1L]]),
      final = unname(fit$loglik[[2L]])
    )
  }
  if (length(coefficient) != 1L || is.na(coefficient) ||
      !is.finite(coefficient)) {
    result$slope_fit <- sb_survival_weak_calibration_component_error(
      result$slope_fit,
      "failed",
      "stormbreaker_survival_weak_calibration_singular_fit",
      "The secondary Cox fit did not return a finite named coefficient."
    )
    return(sb_survival_weak_calibration_finish(result))
  }

  result$slope_fit$status <- "success"
  result$slope_fit$coefficient <- coefficient
  result$calibration_slope <- coefficient
  sb_survival_weak_calibration_finish(result)
}
