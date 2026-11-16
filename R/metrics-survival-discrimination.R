# ---- metrics-survival-discrimination ----

# ---- metrics-concordance ----

sb_harrell_c_stratum <- function(time, event, risk) {
  data <- data.frame(
    time = time,
    event = event,
    risk = risk
  )
  result <- list(
    status = "not_evaluable",
    estimate = NA_real_,
    n = as.integer(length(time)),
    n_events = as.integer(sum(event)),
    concordant = 0,
    discordant = 0,
    tied_risk = 0,
    tied_time = 0,
    tied_both = 0,
    comparable_pairs = 0,
    error_class = NA_character_,
    error_message = "No comparable patient pairs."
  )
  concordance <- tryCatch(
    survival::concordance(
      survival::Surv(time, event) ~ risk,
      data = data,
      reverse = TRUE,
      timewt = "n"
    ),
    error = function(error) error
  )
  if (inherits(concordance, "error")) {
    result$status <- "metric_failed"
    result$error_class <- class(concordance)[1L]
    result$error_message <- conditionMessage(concordance)
    return(result)
  }

  counts <- concordance$count
  result$concordant <- unname(counts[["concordant"]])
  result$discordant <- unname(counts[["discordant"]])
  result$tied_risk <- unname(counts[["tied.x"]])
  result$tied_time <- unname(counts[["tied.y"]])
  result$tied_both <- unname(counts[["tied.xy"]])
  result$comparable_pairs <- unname(
    result$concordant + result$discordant + result$tied_risk
  )
  if (!is.finite(result$comparable_pairs) ||
      result$comparable_pairs <= 0) {
    return(result)
  }

  result$status <- "success"
  result$estimate <- unname(concordance$concordance)
  result$error_message <- NA_character_
  result
}

sb_harrell_c_metric <- function(predictions) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  per_fold <- lapply(
    fold_ids,
    function(fold_id) {
      fold_predictions <- predictions[
        predictions$fold == fold_id,
        c("time", "event", "risk", "status"),
        drop = FALSE
      ]
      audit <- data.frame(
        fold = fold_id,
        status = "prediction_failed",
        n = nrow(fold_predictions),
        n_events = as.integer(sum(fold_predictions$event)),
        estimate = NA_real_,
        concordant = 0,
        discordant = 0,
        tied_risk = 0,
        tied_time = 0,
        tied_both = 0,
        comparable_pairs = 0,
        error_class = NA_character_,
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
      prediction_failed <- anyNA(fold_predictions$status) ||
        any(fold_predictions$status != "success") ||
        anyNA(fold_predictions$risk) ||
        any(!is.finite(fold_predictions$risk))
      if (prediction_failed) {
        audit$error_message <- paste0(
          "OOF predictions for this fold are failed or non-finite."
        )
        return(audit)
      }

      stratum <- sb_harrell_c_stratum(
        time = fold_predictions$time,
        event = fold_predictions$event,
        risk = fold_predictions$risk
      )
      audit$n <- stratum$n
      audit$n_events <- stratum$n_events
      audit$concordant <- stratum$concordant
      audit$discordant <- stratum$discordant
      audit$tied_risk <- stratum$tied_risk
      audit$tied_time <- stratum$tied_time
      audit$tied_both <- stratum$tied_both
      audit$comparable_pairs <- stratum$comparable_pairs
      if (identical(stratum$status, "metric_failed")) {
        audit$status <- "metric_failed"
        audit$error_class <- stratum$error_class
        audit$error_message <- stratum$error_message
        return(audit)
      }
      if (identical(stratum$status, "not_evaluable")) {
        audit$status <- "not_evaluable"
        audit$error_message <- "No comparable patient pairs in this fold."
        return(audit)
      }

      audit$status <- "success"
      audit$estimate <- stratum$estimate
      audit
    }
  )
  per_fold <- do.call(rbind, per_fold)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  count_columns <- c(
    "concordant", "discordant", "tied_risk", "tied_time",
    "tied_both", "comparable_pairs"
  )
  totals <- stats::setNames(numeric(length(count_columns)), count_columns)
  if (any(evaluable)) {
    totals <- colSums(per_fold[evaluable, count_columns, drop = FALSE])
  }
  complete <- length(fold_ids) > 0L && all(evaluable)
  estimate <- if (complete) {
    unname(
      (totals["concordant"] + 0.5 * totals["tied_risk"]) /
        totals["comparable_pairs"]
    )
  } else {
    NA_real_
  }

  list(
    metric = "harrell_c",
    scope = "within_outer_assessment_fold",
    aggregation = "summed_comparable_pairs",
    time_weight = "n",
    risk_direction = "higher_is_higher_hazard",
    status = if (complete) "success" else "incomplete",
    estimate = estimate,
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    concordant = unname(totals["concordant"]),
    discordant = unname(totals["discordant"]),
    tied_risk = unname(totals["tied_risk"]),
    tied_time = unname(totals["tied_time"]),
    tied_both = unname(totals["tied_both"]),
    comparable_pairs = unname(totals["comparable_pairs"]),
    per_fold = per_fold
  )
}

# ---- metrics-survival-time-auc ----

sb_time_auc_audit <- function(fold, prediction_horizon, time, event) {
  valid_endpoint <- is.numeric(time) &&
    (is.numeric(event) || is.logical(event)) &&
    length(event) == length(time) &&
    !anyNA(time) && !anyNA(event) &&
    all(is.finite(time)) && all(time >= 0) &&
    (!is.numeric(event) || all(is.finite(event))) &&
    all(event %in% c(0, 1))
  valid_horizon <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
    prediction_horizon > 0
  if (valid_endpoint && valid_horizon) {
    event <- as.integer(event)
    cases <- event == 1L & time <= prediction_horizon
    controls <- time > prediction_horizon
    excluded_censored <- event == 0L & time <= prediction_horizon
  } else {
    cases <- controls <- excluded_censored <- logical(length(time))
  }
  data.frame(
    fold = as.integer(fold),
    status = "input_failed",
    prediction_horizon = as.numeric(prediction_horizon)[[1L]],
    n = as.integer(length(time)),
    n_cases = as.integer(sum(cases)),
    n_controls = as.integer(sum(controls)),
    n_excluded_censored = as.integer(sum(excluded_censored)),
    case_weight_sum = 0,
    weighted_concordant_pairs = 0,
    weighted_tied_pairs = 0,
    weighted_pair_credit = 0,
    pair_denominator = 0,
    estimate = NA_real_,
    error_class = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

sb_time_auc_fail <- function(audit, status, error_class, error_message) {
  audit$status <- status
  audit$error_class <- error_class
  audit$error_message <- error_message
  audit
}

sb_time_auc_stratum <- function(time, event, survival_probability,
                                 prediction_ok, censoring_curve,
                                 prediction_horizon, fold = NA_integer_,
                                 log_cumulative_hazard = NULL) {
  audit <- sb_time_auc_audit(fold, prediction_horizon, time, event)
  valid_endpoint <- is.numeric(time) &&
    (is.numeric(event) || is.logical(event)) &&
    length(event) == length(time) &&
    !anyNA(time) && !anyNA(event) &&
    all(is.finite(time)) && all(time >= 0) &&
    (!is.numeric(event) || all(is.finite(event))) &&
    all(event %in% c(0, 1))
  valid_horizon <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
    prediction_horizon > 0
  if (!valid_endpoint || !valid_horizon) {
    return(sb_time_auc_fail(
      audit,
      "input_failed",
      "stormbreaker_bad_time_auc_input",
      "Time AUC requires a valid survival endpoint and prediction horizon."
    ))
  }

  valid_prediction <- is.numeric(survival_probability) &&
    length(survival_probability) == length(time) &&
    !anyNA(survival_probability) &&
    all(is.finite(survival_probability)) &&
    all(survival_probability >= 0 & survival_probability <= 1) &&
    is.logical(prediction_ok) && length(prediction_ok) == length(time) &&
    !anyNA(prediction_ok) && all(prediction_ok)
  if (!valid_prediction) {
    return(sb_time_auc_fail(
      audit,
      "prediction_failed",
      "stormbreaker_invalid_survival_prediction",
      paste0(
        "Locked-horizon survival probabilities are failed, missing, or ",
        "invalid."
      )
    ))
  }
  if (audit$n_cases == 0L) {
    return(sb_time_auc_fail(
      audit,
      "not_evaluable",
      "stormbreaker_no_time_auc_case",
      "No observed event occurred at or before the prediction horizon."
    ))
  }
  if (audit$n_controls == 0L) {
    return(sb_time_auc_fail(
      audit,
      "not_evaluable",
      "stormbreaker_no_time_auc_control",
      "No patient has observed follow-up beyond the prediction horizon."
    ))
  }

  valid_curve <- is.data.frame(censoring_curve) &&
    nrow(censoring_curve) > 0L &&
    all(c("time", "survival_after") %in% names(censoring_curve))
  if (!valid_curve) {
    return(sb_time_auc_fail(
      audit,
      "censoring_failed",
      "stormbreaker_missing_censoring_curve",
      "The required marginal censoring curve is unavailable."
    ))
  }

  event <- as.integer(event)
  cases <- event == 1L & time <= prediction_horizon
  controls <- time > prediction_horizon
  case_survival <- vapply(
    time[cases],
    function(value) {
      sb_censoring_survival_at(
        censoring_curve,
        value,
        left_limit = TRUE
      )
    },
    numeric(1)
  )
  if (anyNA(case_survival) || any(!is.finite(case_survival)) ||
      any(case_survival <= 0)) {
    return(sb_time_auc_fail(
      audit,
      "censoring_failed",
      "stormbreaker_zero_censoring_survival",
      paste0(
        "The marginal censoring survival is non-positive or non-finite ",
        "before an observed case event."
      )
    ))
  }
  case_weights <- 1 / case_survival
  if (anyNA(case_weights) || any(!is.finite(case_weights))) {
    return(sb_time_auc_fail(
      audit,
      "censoring_failed",
      "stormbreaker_invalid_time_auc_weight",
      "A required time-AUC case weight is missing or non-finite."
    ))
  }

  marker <- if (is.null(log_cumulative_hazard)) {
    1 - survival_probability
  } else sb_validate_stable_survival_scale(log_cumulative_hazard, survival_probability)
  case_marker <- marker[cases]
  control_marker <- marker[controls]
  marker_values <- sort(unique(c(case_marker, control_marker)))
  case_groups <- match(case_marker, marker_values)
  control_groups <- match(control_marker, marker_values)
  control_counts <- tabulate(
    control_groups,
    nbins = length(marker_values)
  )
  case_weight_by_group <- vapply(
    seq_along(marker_values),
    function(group) sum(case_weights[case_groups == group]),
    numeric(1)
  )
  controls_below <- cumsum(control_counts) - control_counts
  weighted_concordant <- unname(sum(
    case_weight_by_group * controls_below
  ))
  weighted_tied <- unname(sum(case_weight_by_group * control_counts))
  weighted_pair_credit <- unname(
    weighted_concordant + 0.5 * weighted_tied
  )
  case_weight_sum <- unname(sum(case_weights))
  pair_denominator <- unname(case_weight_sum * audit$n_controls)
  valid_totals <- all(is.finite(c(
    weighted_concordant,
    weighted_tied,
    weighted_pair_credit,
    case_weight_sum,
    pair_denominator
  ))) && pair_denominator > 0
  if (!valid_totals) {
    return(sb_time_auc_fail(
      audit,
      "censoring_failed",
      "stormbreaker_invalid_time_auc_weight",
      "The weighted time-AUC pair totals are non-finite or non-positive."
    ))
  }

  audit$status <- "success"
  audit$case_weight_sum <- case_weight_sum
  audit$weighted_concordant_pairs <- weighted_concordant
  audit$weighted_tied_pairs <- weighted_tied
  audit$weighted_pair_credit <- weighted_pair_credit
  audit$pair_denominator <- pair_denominator
  audit$estimate <- unname(weighted_pair_credit / pair_denominator)
  audit
}

sb_survival_time_auc_metric <- function(predictions, analysis_inputs,
                                         prediction_horizon,
                                         censoring_curves) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  input_folds <- vapply(analysis_inputs, function(input) {
    fold <- input[["fold"]]
    if (is.numeric(fold) && length(fold) == 1L && !is.na(fold) &&
        is.finite(fold)) {
      as.integer(fold)
    } else {
      NA_integer_
    }
  }, integer(1))
  curve_names <- names(censoring_curves)
  if (is.null(curve_names)) {
    curve_names <- rep("", length(censoring_curves))
  }

  fold_results <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      ,
      drop = FALSE
    ]
    prediction_ok <- !is.na(fold_predictions$status) &
      fold_predictions$status == "success" &
      !is.na(fold_predictions$survival_probability_status) &
      fold_predictions$survival_probability_status == "success"
    empty_audit <- function() {
      sb_time_auc_stratum(
        fold_predictions$time,
        fold_predictions$event,
        fold_predictions$survival_probability,
        prediction_ok,
        sb_empty_censoring_curve(),
        prediction_horizon,
        fold_id
      )
    }

    matched_input <- which(input_folds == fold_id)
    if (length(matched_input) != 1L) {
      return(sb_time_auc_fail(
        empty_audit(),
        "analysis_input_failed",
        "stormbreaker_bad_time_auc_analysis_input",
        paste0(
          "Exactly one outer-analysis censoring input is required for this ",
          "fold."
        )
      ))
    }
    matched_curve <- which(curve_names == as.character(fold_id))
    if (length(matched_curve) == 0L) {
      return(sb_time_auc_fail(
        empty_audit(),
        "censoring_failed",
        "stormbreaker_missing_censoring_curve",
        "The matching outer-analysis censoring curve is unavailable."
      ))
    }
    if (length(matched_curve) != 1L) {
      return(sb_time_auc_fail(
        empty_audit(),
        "censoring_failed",
        "stormbreaker_bad_time_auc_censoring_curve",
        "Exactly one matching outer-analysis censoring curve is required."
      ))
    }

    sb_time_auc_stratum(
      fold_predictions$time,
      fold_predictions$event,
      fold_predictions$survival_probability,
      prediction_ok,
      censoring_curves[[matched_curve]],
      prediction_horizon,
      fold_id
    )
  })
  per_fold <- do.call(rbind, fold_results)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  weighted_concordant <- unname(sum(
    per_fold$weighted_concordant_pairs[evaluable]
  ))
  weighted_tied <- unname(sum(per_fold$weighted_tied_pairs[evaluable]))
  weighted_pair_credit <- unname(sum(
    per_fold$weighted_pair_credit[evaluable]
  ))
  pair_denominator <- unname(sum(per_fold$pair_denominator[evaluable]))

  list(
    metric = "time_auc",
    estimand = "cumulative_dynamic_auc",
    scope = "fold_stratified_out_of_fold_at_declared_horizon",
    marker = "one_minus_survival_probability",
    marker_direction = "higher_is_higher_event_risk",
    case_definition = "observed_event_time_less_than_or_equal_to_horizon",
    control_definition = "observed_follow_up_strictly_greater_than_horizon",
    excluded = "censored_at_or_before_horizon",
    censoring_model = "outer_analysis_reverse_kaplan_meier",
    censoring_ties = "events_before_censoring",
    case_weight = "one_over_censoring_survival_at_event_left_limit",
    control_weight = "constant_at_horizon_cancelled_not_evaluated",
    marker_ties = "exact_half_credit",
    aggregation = "summed_within_fold_weighted_pair_credit",
    prediction_horizon = as.numeric(prediction_horizon),
    status = if (complete) "success" else "incomplete",
    estimate = if (complete) {
      unname(weighted_pair_credit / pair_denominator)
    } else {
      NA_real_
    },
    n = as.integer(sum(per_fold$n)),
    n_cases = as.integer(sum(per_fold$n_cases)),
    n_controls = as.integer(sum(per_fold$n_controls)),
    n_excluded_censored = as.integer(sum(
      per_fold$n_excluded_censored
    )),
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    case_weight_sum = unname(sum(per_fold$case_weight_sum[evaluable])),
    weighted_concordant_pairs = weighted_concordant,
    weighted_tied_pairs = weighted_tied,
    weighted_pair_credit = weighted_pair_credit,
    pair_denominator = pair_denominator,
    per_fold = per_fold
  )
}

# ---- metrics-survival-uno-c ----

sb_uno_c_audit <- function(fold, prediction_horizon, time, event) {
  valid_endpoint <- is.numeric(time) &&
    (is.numeric(event) || is.logical(event)) &&
    length(event) == length(time) &&
    !anyNA(time) && !anyNA(event) &&
    all(is.finite(time)) && all(time >= 0) &&
    (!is.numeric(event) || all(is.finite(event))) &&
    all(event %in% c(0, 1))
  valid_horizon <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
    prediction_horizon > 0
  if (valid_endpoint && valid_horizon) {
    event <- as.integer(event)
    eligible <- event == 1L & time < prediction_horizon
    excluded_at_horizon <- event == 1L & time == prediction_horizon
  } else {
    eligible <- excluded_at_horizon <- logical(length(time))
  }

  data.frame(
    fold = as.integer(fold),
    status = "input_failed",
    prediction_horizon = if (valid_horizon) {
      as.numeric(prediction_horizon)
    } else {
      NA_real_
    },
    n = as.integer(length(time)),
    n_eligible_events = as.integer(sum(eligible)),
    n_events_at_horizon_excluded = as.integer(sum(excluded_at_horizon)),
    n_comparable_pairs = 0L,
    n_tied_time_comparable_pairs = 0L,
    censoring_survival_horizon = NA_real_,
    weighted_concordant_pairs = 0,
    weighted_tied_risk_pairs = 0,
    weighted_pair_credit = 0,
    weighted_pair_denominator = 0,
    estimate = NA_real_,
    error_class = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

sb_uno_c_fail <- sb_time_auc_fail

sb_uno_c_stratum <- function(time, event, risk, prediction_ok,
                              censoring_curve, prediction_horizon,
                              fold = NA_integer_) {
  audit <- sb_uno_c_audit(fold, prediction_horizon, time, event)
  valid_endpoint <- is.numeric(time) &&
    (is.numeric(event) || is.logical(event)) &&
    length(event) == length(time) &&
    !anyNA(time) && !anyNA(event) &&
    all(is.finite(time)) && all(time >= 0) &&
    (!is.numeric(event) || all(is.finite(event))) &&
    all(event %in% c(0, 1))
  valid_horizon <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
    prediction_horizon > 0
  if (!valid_endpoint || !valid_horizon) {
    return(sb_uno_c_fail(
      audit,
      "input_failed",
      "stormbreaker_bad_uno_c_input",
      "Uno C requires a valid survival endpoint and prediction horizon."
    ))
  }

  valid_prediction <- is.numeric(risk) && length(risk) == length(time) &&
    !anyNA(risk) && all(is.finite(risk)) &&
    is.logical(prediction_ok) && length(prediction_ok) == length(time) &&
    !anyNA(prediction_ok) && all(prediction_ok)
  if (!valid_prediction) {
    return(sb_uno_c_fail(
      audit,
      "prediction_failed",
      "stormbreaker_invalid_uno_c_risk",
      "Uno C requires one finite successful risk score per patient."
    ))
  }
  if (audit$n_eligible_events == 0L) {
    return(sb_uno_c_fail(
      audit,
      "not_evaluable",
      "stormbreaker_no_uno_c_case",
      "No observed event occurred strictly before the prediction horizon."
    ))
  }

  valid_curve <- is.data.frame(censoring_curve) &&
    nrow(censoring_curve) > 0L &&
    all(c("time", "survival_after") %in% names(censoring_curve))
  if (!valid_curve) {
    return(sb_uno_c_fail(
      audit,
      "censoring_failed",
      "stormbreaker_missing_uno_c_censoring_curve",
      "The required marginal censoring curve is unavailable."
    ))
  }

  horizon_survival <- sb_censoring_survival_at(
    censoring_curve,
    prediction_horizon,
    left_limit = TRUE
  )
  audit$censoring_survival_horizon <- horizon_survival
  if (length(horizon_survival) != 1L || is.na(horizon_survival) ||
      !is.finite(horizon_survival) || horizon_survival <= 0) {
    return(sb_uno_c_fail(
      audit,
      "censoring_failed",
      "stormbreaker_zero_uno_c_horizon_censoring_survival",
      paste0(
        "The marginal censoring survival immediately before the prediction ",
        "horizon is non-positive or non-finite."
      )
    ))
  }

  event <- as.integer(event)
  eligible <- which(event == 1L & time < prediction_horizon)
  event_survival <- vapply(
    time[eligible],
    function(value) {
      sb_censoring_survival_at(
        censoring_curve,
        value,
        left_limit = TRUE
      )
    },
    numeric(1)
  )
  if (anyNA(event_survival) || any(!is.finite(event_survival)) ||
      any(event_survival <= 0)) {
    return(sb_uno_c_fail(
      audit,
      "censoring_failed",
      "stormbreaker_zero_uno_c_event_censoring_survival",
      paste0(
        "The marginal censoring survival immediately before an eligible ",
        "event is non-positive or non-finite."
      )
    ))
  }
  case_weights <- 1 / event_survival^2
  if (anyNA(case_weights) || any(!is.finite(case_weights))) {
    return(sb_uno_c_fail(
      audit,
      "censoring_failed",
      "stormbreaker_invalid_uno_c_weight",
      "A required squared Uno IPCW weight is missing or non-finite."
    ))
  }

  weighted_concordant <- 0
  weighted_tied_risk <- 0
  weighted_denominator <- 0
  n_comparable <- 0L
  n_tied_time <- 0L
  for (position in seq_along(eligible)) {
    index <- eligible[[position]]
    comparators <- which(
      time > time[[index]] |
        (time == time[[index]] & event == 0L)
    )
    if (length(comparators) == 0L) {
      next
    }
    weight <- case_weights[[position]]
    n_comparable <- n_comparable + length(comparators)
    n_tied_time <- n_tied_time + sum(
      time[comparators] == time[[index]]
    )
    weighted_concordant <- weighted_concordant + weight * sum(
      risk[[index]] > risk[comparators]
    )
    weighted_tied_risk <- weighted_tied_risk + weight * sum(
      risk[[index]] == risk[comparators]
    )
    weighted_denominator <- weighted_denominator +
      weight * length(comparators)
  }

  audit$n_comparable_pairs <- as.integer(n_comparable)
  audit$n_tied_time_comparable_pairs <- as.integer(n_tied_time)
  if (n_comparable == 0L) {
    return(sb_uno_c_fail(
      audit,
      "not_evaluable",
      "stormbreaker_no_uno_c_comparable_pair",
      "No patient pair is comparable for Uno C before the horizon."
    ))
  }

  weighted_credit <- weighted_concordant + 0.5 * weighted_tied_risk
  valid_totals <- all(is.finite(c(
    weighted_concordant,
    weighted_tied_risk,
    weighted_credit,
    weighted_denominator
  ))) && weighted_denominator > 0
  if (!valid_totals) {
    return(sb_uno_c_fail(
      audit,
      "censoring_failed",
      "stormbreaker_invalid_uno_c_weight",
      "The weighted Uno C pair totals are non-finite or non-positive."
    ))
  }

  audit$status <- "success"
  audit$weighted_concordant_pairs <- unname(weighted_concordant)
  audit$weighted_tied_risk_pairs <- unname(weighted_tied_risk)
  audit$weighted_pair_credit <- unname(weighted_credit)
  audit$weighted_pair_denominator <- unname(weighted_denominator)
  audit$estimate <- unname(weighted_credit / weighted_denominator)
  audit
}

sb_survival_uno_c_metric <- function(predictions, analysis_inputs,
                                      prediction_horizon,
                                      censoring_curves) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  input_folds <- vapply(analysis_inputs, function(input) {
    fold <- input[["fold"]]
    if (is.numeric(fold) && length(fold) == 1L && !is.na(fold) &&
        is.finite(fold)) {
      as.integer(fold)
    } else {
      NA_integer_
    }
  }, integer(1))
  curve_names <- names(censoring_curves)
  if (is.null(curve_names)) {
    curve_names <- rep("", length(censoring_curves))
  }

  fold_results <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      ,
      drop = FALSE
    ]
    prediction_ok <- !is.na(fold_predictions$status) &
      fold_predictions$status == "success"
    empty_audit <- function() {
      sb_uno_c_audit(
        fold_id,
        prediction_horizon,
        fold_predictions$time,
        fold_predictions$event
      )
    }

    matched_input <- which(input_folds == fold_id)
    if (length(matched_input) != 1L) {
      return(sb_uno_c_fail(
        empty_audit(),
        "analysis_input_failed",
        "stormbreaker_bad_uno_c_analysis_input",
        paste0(
          "Exactly one outer-analysis censoring input is required for this ",
          "fold."
        )
      ))
    }
    matched_curve <- which(curve_names == as.character(fold_id))
    if (length(matched_curve) == 0L) {
      return(sb_uno_c_fail(
        empty_audit(),
        "censoring_failed",
        "stormbreaker_missing_uno_c_censoring_curve",
        "The matching outer-analysis censoring curve is unavailable."
      ))
    }
    if (length(matched_curve) != 1L) {
      return(sb_uno_c_fail(
        empty_audit(),
        "censoring_failed",
        "stormbreaker_bad_uno_c_censoring_curve",
        "Exactly one matching outer-analysis censoring curve is required."
      ))
    }

    result <- tryCatch(
      sb_uno_c_stratum(
        fold_predictions$time,
        fold_predictions$event,
        fold_predictions$risk,
        prediction_ok,
        censoring_curves[[matched_curve]],
        prediction_horizon,
        fold_id
      ),
      error = function(error) error
    )
    if (inherits(result, "condition")) {
      return(sb_uno_c_fail(
        empty_audit(),
        "metric_failed",
        class(result)[[1L]],
        conditionMessage(result)
      ))
    }
    result
  })

  per_fold <- do.call(rbind, fold_results)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  weighted_concordant <- unname(sum(
    per_fold$weighted_concordant_pairs[evaluable]
  ))
  weighted_tied_risk <- unname(sum(
    per_fold$weighted_tied_risk_pairs[evaluable]
  ))
  weighted_credit <- unname(sum(
    per_fold$weighted_pair_credit[evaluable]
  ))
  weighted_denominator <- unname(sum(
    per_fold$weighted_pair_denominator[evaluable]
  ))

  list(
    metric = "uno_c",
    estimand = "truncated_concordance_probability",
    scope = "fold_stratified_out_of_fold_at_declared_horizon",
    marker = "risk",
    marker_direction = "higher_is_higher_event_risk",
    case_definition = "observed_event_time_strictly_less_than_horizon",
    comparator_definition = paste0(
      "later_observed_time_or_censoring_at_same_observed_time"
    ),
    censoring_model = "outer_analysis_reverse_kaplan_meier",
    censoring_assumption = "marginal_random_censoring",
    censoring_ties = "events_before_censoring",
    case_weight = "one_over_squared_censoring_survival_at_event_left_limit",
    time_ties = "event_censor_comparable_event_event_not_comparable",
    marker_ties = "exact_half_credit",
    aggregation = "summed_within_fold_weighted_pair_credit",
    prediction_horizon = as.numeric(prediction_horizon),
    status = if (complete) "success" else "incomplete",
    estimate = if (complete && weighted_denominator > 0) {
      unname(weighted_credit / weighted_denominator)
    } else {
      NA_real_
    },
    n = as.integer(sum(per_fold$n)),
    n_eligible_events = as.integer(sum(per_fold$n_eligible_events)),
    n_events_at_horizon_excluded = as.integer(sum(
      per_fold$n_events_at_horizon_excluded
    )),
    n_comparable_pairs = as.integer(sum(per_fold$n_comparable_pairs)),
    n_tied_time_comparable_pairs = as.integer(sum(
      per_fold$n_tied_time_comparable_pairs
    )),
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    weighted_concordant_pairs = weighted_concordant,
    weighted_tied_risk_pairs = weighted_tied_risk,
    weighted_pair_credit = weighted_credit,
    weighted_pair_denominator = weighted_denominator,
    per_fold = per_fold
  )
}

# ---- external-metrics ----

sb_external_binary_auc_metric <- function(scored) {
  stratum <- sb_binary_auc_stratum(
    outcome = scored$outcome,
    probability = scored$probability
  )
  result <- list(
    metric = "roc_auc",
    scope = "pooled_external_cohort",
    aggregation = "all_external_positive_negative_pairs",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = "not_evaluable",
    estimate = NA_real_,
    n = stratum$n,
    n_positive = stratum$n_positive,
    n_negative = stratum$n_negative,
    correctly_ordered_pairs = 0,
    tied_pairs = 0,
    comparable_pairs = stratum$comparable_pairs,
    error_class = NA_character_,
    error_message = "Both outcome classes are required."
  )
  if (!identical(stratum$status, "success")) {
    return(result)
  }

  result$status <- "success"
  result$estimate <- stratum$estimate
  result$correctly_ordered_pairs <- stratum$correctly_ordered_pairs
  result$tied_pairs <- stratum$tied_pairs
  result$error_message <- NA_character_
  result
}

sb_external_binary_average_precision_metric <- function(scored) {
  outcome <- as.integer(scored$outcome)
  probability <- scored$probability
  n_positive <- as.integer(sum(outcome == 1L))
  n_negative <- as.integer(sum(outcome == 0L))
  empty_curve <- sb_empty_binary_precision_recall_curve()
  empty_curve$fold <- NULL
  result <- list(
    metric = "average_precision",
    scope = "pooled_external_cohort",
    aggregation = "positive_weighted_precision",
    definition = "sum((recall_i - recall_previous) * precision_i)",
    interpolation = "none",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = "not_evaluable",
    estimate = NA_real_,
    n = as.integer(length(outcome)),
    n_positive = n_positive,
    n_negative = n_negative,
    n_thresholds = 0L,
    positive_weighted_sum = 0,
    curve = empty_curve,
    error_class = NA_character_,
    error_message = "Both outcome classes are required."
  )
  if (n_positive == 0L || n_negative == 0L) {
    return(result)
  }
  curve <- sb_binary_precision_recall_curve(outcome, probability, 1L)
  curve$fold <- NULL
  positive_weighted_sum <- sum(
    curve$positive_at_threshold * curve$precision
  )
  result$status <- "success"
  result$estimate <- unname(positive_weighted_sum / n_positive)
  result$n_thresholds <- as.integer(nrow(curve))
  result$positive_weighted_sum <- unname(positive_weighted_sum)
  result$curve <- curve
  result$error_message <- NA_character_
  result
}

sb_external_binary_brier_metric <- function(scored) {
  squared_error <- (as.integer(scored$outcome) - scored$probability)^2
  list(
    metric = "brier_score",
    scope = "patient_level_external_cohort",
    aggregation = "mean_squared_error_over_external_patients",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    metric_direction = "lower_is_better",
    status = "success",
    estimate = unname(mean(squared_error)),
    n = as.integer(nrow(scored)),
    n_scored = as.integer(nrow(scored)),
    squared_error_sum = unname(sum(squared_error))
  )
}

sb_external_binary_log_loss_metric <- function(scored) {
  case_loss <- sb_binary_case_log_loss(
    scored$outcome,
    scored$probability
  )
  list(
    metric = "log_loss",
    scope = "patient_level_external_cohort",
    aggregation = "mean_case_loss_over_external_patients",
    definition = "negative_binary_log_likelihood",
    probability_clipping = "none",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    metric_direction = "lower_is_better",
    status = "success",
    estimate = unname(mean(case_loss)),
    n = as.integer(nrow(scored)),
    n_scored = as.integer(nrow(scored)),
    n_finite_loss = as.integer(sum(is.finite(case_loss))),
    n_infinite_loss = as.integer(sum(is.infinite(case_loss))),
    loss_sum = unname(sum(case_loss))
  )
}

sb_external_binary_mean_calibration_metric <- function(scored) {
  outcome <- as.integer(scored$outcome)
  probability <- scored$probability
  observed_events <- as.integer(sum(outcome == 1L))
  expected_events <- sum(probability)
  result <- list(
    metric = "mean_calibration_ratio",
    calibration_level = "mean",
    estimand = "binary_observed_over_expected_event_ratio",
    scope = "pooled_external_cohort",
    purpose = "descriptive_evaluation_only",
    ratio_direction = "observed_over_expected",
    probability_clipping = "none",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = "not_evaluable",
    estimate = NA_real_,
    mean_calibration_ratio = NA_real_,
    target = 1,
    n = as.integer(nrow(scored)),
    observed_events = observed_events,
    expected_events = unname(expected_events),
    observed_event_probability = unname(mean(outcome)),
    mean_predicted_event_probability = unname(mean(probability)),
    error_class = NA_character_,
    error_message = NA_character_
  )
  if (expected_events == 0) {
    result$error_class <-
      "stormbreaker_binary_mean_calibration_zero_expected"
    result$error_message <- paste0(
      "The total expected number of events is zero, so the observed-to-",
      "expected ratio is undefined."
    )
    return(result)
  }

  ratio <- observed_events / expected_events
  if (!is.finite(ratio)) {
    result$error_class <-
      "stormbreaker_binary_mean_calibration_nonfinite_ratio"
    result$error_message <- paste0(
      "The observed-to-expected ratio is not finite on the available ",
      "numeric scale."
    )
    return(result)
  }
  result$status <- "success"
  result$estimate <- unname(ratio)
  result$mean_calibration_ratio <- unname(ratio)
  result
}

sb_external_binary_calibration_metric <- function(scored) {
  outcome <- as.integer(scored$outcome)
  probability <- scored$probability
  n_positive <- as.integer(sum(outcome == 1L))
  n_negative <- as.integer(sum(outcome == 0L))
  n_boundary <- as.integer(sum(probability %in% c(0, 1)))
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
    scope = "pooled_external_cohort",
    purpose = "descriptive_evaluation_only",
    probability_clipping = "none",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = "incomplete",
    input_status = "success",
    calibration_intercept = NA_real_,
    calibration_slope = NA_real_,
    intercept_target = 0,
    slope_target = 1,
    n = as.integer(nrow(scored)),
    n_positive = n_positive,
    n_negative = n_negative,
    n_boundary_probability = n_boundary,
    intercept_fit = intercept_fit,
    slope_fit = slope_fit,
    error_class = NA_character_,
    error_message = NA_character_
  )
  if (n_positive == 0L || n_negative == 0L) {
    result$input_status <- "not_evaluable"
    result$error_message <- "Both outcome classes are required for calibration."
    return(result)
  }
  if (n_boundary > 0L) {
    result$input_status <- "not_evaluable"
    result$error_class <-
      "stormbreaker_binary_calibration_boundary_probability"
    result$error_message <- paste0(
      "Exact boundary probabilities cannot enter finite-logit calibration ",
      "without clipping, which is not performed."
    )
    return(result)
  }

  logit_probability <- stats::qlogis(probability)
  intercept_fit <- sb_binary_fit_calibration(
    outcome,
    logit_probability,
    "calibration_in_the_large"
  )
  slope_fit <- sb_binary_fit_calibration(
    outcome,
    logit_probability,
    "logistic_recalibration"
  )
  result$intercept_fit <- intercept_fit
  result$slope_fit <- slope_fit
  if (identical(intercept_fit$status, "success")) {
    result$calibration_intercept <- intercept_fit$coefficients$estimate[[1L]]
  }
  if (identical(slope_fit$status, "success")) {
    result$calibration_slope <- slope_fit$coefficients$estimate[[2L]]
  }
  if (identical(intercept_fit$status, "success") &&
      identical(slope_fit$status, "success")) {
    result$status <- "success"
  } else {
    failed_fit <- if (!identical(intercept_fit$status, "success")) {
      intercept_fit
    } else {
      slope_fit
    }
    result$error_class <- failed_fit$error_class
    result$error_message <- paste(
      stats::na.omit(c(
        intercept_fit$error_message,
        slope_fit$error_message
      )),
      collapse = " "
    )
  }
  result
}

sb_external_harrell_c_metric <- function(scored) {
  result <- list(
    metric = "harrell_c",
    scope = "external_cohort",
    aggregation = "all_external_comparable_pairs",
    time_weight = "n",
    risk_direction = "higher_is_higher_hazard",
    status = "not_evaluable",
    estimate = NA_real_,
    n = as.integer(nrow(scored)),
    n_events = as.integer(sum(scored$event == 1L)),
    concordant = 0,
    discordant = 0,
    tied_risk = 0,
    tied_time = 0,
    tied_both = 0,
    comparable_pairs = 0,
    error_class = NA_character_,
    error_message = "No comparable patient pairs in the external endpoint."
  )
  if (nrow(scored) < 2L || sum(scored$event == 1L) == 0L) {
    return(result)
  }
  stratum <- sb_harrell_c_stratum(
    time = scored$time,
    event = scored$event,
    risk = scored$risk
  )
  result$concordant <- stratum$concordant
  result$discordant <- stratum$discordant
  result$tied_risk <- stratum$tied_risk
  result$tied_time <- stratum$tied_time
  result$tied_both <- stratum$tied_both
  result$comparable_pairs <- stratum$comparable_pairs
  if (identical(stratum$status, "metric_failed")) {
    result$status <- "metric_failed"
    result$error_class <- stratum$error_class
    result$error_message <- stratum$error_message
    return(result)
  }
  if (identical(stratum$status, "not_evaluable")) {
    return(result)
  }
  result$status <- "success"
  result$estimate <- unname(
    (result$concordant + 0.5 * result$tied_risk) /
      result$comparable_pairs
  )
  result$error_message <- NA_character_
  result
}

sb_external_ipcw_brier_metric <- function(scored, endpoint) {
  horizon <- unique(as.numeric(scored$prediction_horizon))
  time <- scored$time
  event <- as.integer(scored$event)
  probability <- scored$survival_probability
  n <- as.integer(nrow(scored))
  curve_result <- tryCatch(
    sb_fit_censoring_km(time, event),
    error = function(error) error
  )
  if (inherits(curve_result, "condition")) {
    return(list(
      metric = "ipcw_brier_score",
      scope = "patient_level_external_cohort_at_locked_horizon",
      censoring_model = "external_reverse_kaplan_meier",
      censoring_assumption = "marginal_noninformative_censoring",
      status = "incomplete",
      estimate = NA_real_,
      prediction_horizon = horizon,
      time_unit = endpoint$time_unit,
      n = n,
      n_scored = 0L,
      error_class = class(curve_result)[[1L]],
      error_message = conditionMessage(curve_result),
      censoring_curve = sb_empty_censoring_curve()
    ))
  }
  curve <- curve_result
  event_by_horizon <- time <= horizon & event == 1L
  surviving_at_horizon <- time > horizon |
    (time == horizon & event == 0L)
  n_scored <- as.integer(sum(event_by_horizon | surviving_at_horizon))
  horizon_survival <- sb_censoring_survival_at(
    curve,
    horizon,
    left_limit = TRUE
  )
  result <- list(
    metric = "ipcw_brier_score",
    scope = "patient_level_external_cohort_at_locked_horizon",
    censoring_model = "external_reverse_kaplan_meier",
    censoring_assumption = "marginal_noninformative_censoring",
    censoring_ties = "events_before_censoring",
    censoring_evaluation = "left_continuous_before_censoring_at_time",
    aggregation = "summed_weighted_loss_over_all_external_patients",
    denominator = "all_external_patients_n",
    prediction_horizon = horizon,
    time_unit = endpoint$time_unit,
    status = "incomplete",
    estimate = NA_real_,
    n = n,
    n_scored = n_scored,
    n_events_by_horizon = as.integer(sum(event_by_horizon)),
    n_surviving_at_horizon = as.integer(sum(surviving_at_horizon)),
    n_censored_before_horizon = as.integer(sum(
      time < horizon & event == 0L
    )),
    n_censored_at_horizon = as.integer(sum(
      time == horizon & event == 0L
    )),
    censoring_survival_horizon = unname(horizon_survival),
    weight_sum = 0,
    weighted_loss_sum = 0,
    censoring_curve = curve,
    error_class = NA_character_,
    error_message = NA_character_
  )
  if (n_scored == 0L) {
    result$error_class <- "stormbreaker_no_evaluable_horizon_status"
    result$error_message <- paste0(
      "No external patient has observed event status at the locked horizon."
    )
    return(result)
  }

  event_survival <- vapply(
    time[event_by_horizon],
    function(value) sb_censoring_survival_at(
      curve,
      value,
      left_limit = TRUE
    ),
    numeric(1)
  )
  valid_event_survival <- all(is.finite(event_survival)) &&
    all(event_survival > 0)
  valid_horizon_survival <- !any(surviving_at_horizon) ||
    (is.finite(horizon_survival) && horizon_survival > 0)
  if (!valid_event_survival || !valid_horizon_survival) {
    result$error_class <- "stormbreaker_zero_external_censoring_survival"
    result$error_message <- paste0(
      "A reverse-Kaplan-Meier denominator required for external IPCW is ",
      "non-positive or non-finite."
    )
    return(result)
  }

  weights <- numeric(n)
  weights[event_by_horizon] <- 1 / event_survival
  if (any(surviving_at_horizon)) {
    weights[surviving_at_horizon] <- 1 / horizon_survival
  }
  losses <- numeric(n)
  losses[event_by_horizon] <- probability[event_by_horizon]^2
  losses[surviving_at_horizon] <-
    (1 - probability[surviving_at_horizon])^2
  weighted_loss <- weights * losses
  if (anyNA(weighted_loss) || any(!is.finite(weighted_loss))) {
    result$error_class <- "stormbreaker_invalid_external_ipcw_loss"
    result$error_message <- "External IPCW losses are missing or non-finite."
    return(result)
  }

  result$status <- "success"
  result$weight_sum <- unname(sum(weights))
  result$weighted_loss_sum <- unname(sum(weighted_loss))
  result$estimate <- unname(result$weighted_loss_sum / n)
  result
}

sb_external_time_auc_metric <- function(scored, endpoint, censoring_curve) {
  audit <- sb_time_auc_stratum(
    time = scored$time,
    event = scored$event,
    survival_probability = scored$survival_probability,
    log_cumulative_hazard = scored$log_cumulative_hazard,
    prediction_ok = !is.na(scored$status) & scored$status == "success",
    censoring_curve = censoring_curve,
    prediction_horizon = unique(as.numeric(scored$prediction_horizon)),
    fold = NA_integer_
  )
  list(
    metric = "time_auc",
    estimand = "cumulative_dynamic_auc",
    scope = "external_cohort_at_locked_horizon",
    marker = if (is.null(scored$log_cumulative_hazard)) {
      "one_minus_survival_probability"
    } else "log_cumulative_hazard",
    marker_direction = "higher_is_higher_event_risk",
    case_definition = "observed_event_time_less_than_or_equal_to_horizon",
    control_definition = "observed_follow_up_strictly_greater_than_horizon",
    excluded = "censored_at_or_before_horizon",
    censoring_model = "external_reverse_kaplan_meier",
    censoring_assumption = "marginal_noninformative_censoring",
    censoring_ties = "events_before_censoring",
    case_weight = "one_over_censoring_survival_at_event_left_limit",
    control_weight = "constant_at_horizon_cancelled_not_evaluated",
    marker_ties = "exact_half_credit",
    aggregation = "all_external_case_control_pairs",
    prediction_horizon = audit$prediction_horizon,
    time_unit = endpoint$time_unit,
    status = audit$status,
    estimate = audit$estimate,
    n = audit$n,
    n_cases = audit$n_cases,
    n_controls = audit$n_controls,
    n_excluded_censored = audit$n_excluded_censored,
    case_weight_sum = audit$case_weight_sum,
    weighted_concordant_pairs = audit$weighted_concordant_pairs,
    weighted_tied_pairs = audit$weighted_tied_pairs,
    weighted_pair_credit = audit$weighted_pair_credit,
    pair_denominator = audit$pair_denominator,
    error_class = audit$error_class,
    error_message = audit$error_message
  )
}

sb_external_uno_c_metric <- function(scored, endpoint, censoring_curve) {
  horizon <- unique(as.numeric(scored$prediction_horizon))
  audit <- tryCatch(
    sb_uno_c_stratum(
      time = scored$time,
      event = scored$event,
      risk = scored$risk,
      prediction_ok = !is.na(scored$status) & scored$status == "success",
      censoring_curve = censoring_curve,
      prediction_horizon = horizon,
      fold = NA_integer_
    ),
    error = function(error) error
  )
  if (inherits(audit, "condition")) {
    audit <- sb_uno_c_fail(
      sb_uno_c_audit(NA_integer_, horizon, scored$time, scored$event),
      "metric_failed",
      class(audit)[[1L]],
      conditionMessage(audit)
    )
  }

  list(
    metric = "uno_c",
    estimand = "truncated_concordance_probability",
    scope = "external_cohort_at_locked_horizon",
    marker = "risk",
    marker_direction = "higher_is_higher_event_risk",
    case_definition = "observed_event_time_strictly_less_than_horizon",
    comparator_definition = paste0(
      "later_observed_time_or_censoring_at_same_observed_time"
    ),
    censoring_model = "external_reverse_kaplan_meier",
    censoring_assumption = "marginal_random_censoring",
    censoring_ties = "events_before_censoring",
    case_weight = "one_over_squared_censoring_survival_at_event_left_limit",
    time_ties = "event_censor_comparable_event_event_not_comparable",
    marker_ties = "exact_half_credit",
    aggregation = "all_external_comparable_pairs",
    prediction_horizon = audit$prediction_horizon,
    time_unit = endpoint$time_unit,
    status = audit$status,
    estimate = audit$estimate,
    n = audit$n,
    n_eligible_events = audit$n_eligible_events,
    n_events_at_horizon_excluded = audit$n_events_at_horizon_excluded,
    n_comparable_pairs = audit$n_comparable_pairs,
    n_tied_time_comparable_pairs = audit$n_tied_time_comparable_pairs,
    censoring_survival_horizon = audit$censoring_survival_horizon,
    weighted_concordant_pairs = audit$weighted_concordant_pairs,
    weighted_tied_risk_pairs = audit$weighted_tied_risk_pairs,
    weighted_pair_credit = audit$weighted_pair_credit,
    weighted_pair_denominator = audit$weighted_pair_denominator,
    error_class = audit$error_class,
    error_message = audit$error_message
  )
}
