# ---- model-ahaz-bujar ----

# ---- ahaz-adapters ----

sb_ahaz_survival_learner_ids <- function() {
  c(
    "ahaz_lin_ying",
    "ahaz_lasso_bic",
    "ahaz_sscad_bic"
  )
}

sb_ahaz_available <- function() {
  requireNamespace("ahaz", quietly = TRUE)
}

sb_require_ahaz_dependency <- function() {
  if (!sb_ahaz_available()) {
    sb_abort(
      "Package `ahaz` is required for additive-hazards learners.",
      "stormbreaker_missing_ahaz"
    )
  }
  invisible(TRUE)
}

sb_ahaz_learner_params <- function(learner) {
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_ahaz_survival_learner_ids()) {
    sb_abort(
      "The requested additive-hazards learner is not supported.",
      "stormbreaker_unsupported_learner"
    )
  }
  if (identical(learner, "ahaz_lin_ying")) {
    return(list(
      estimator = "lin_ying",
      robust = FALSE
    ))
  }
  if (identical(learner, "ahaz_lasso_bic")) {
    return(list(
      estimator = "penalized",
      penalty = "lasso",
      alpha = 1,
      standardize = TRUE,
      tuning = "bic",
      nlambda = 100L,
      dfmax_rule = "min_30_p"
    ))
  }
  list(
    estimator = "penalized",
    penalty = "sscad",
    a = 3.7,
    nsteps = 1L,
    standardize = TRUE,
    tuning = "bic",
    nlambda = 100L,
    dfmax_rule = "min_30_p"
  )
}

sb_ahaz_parameter_semantics <- function(learner) {
  params <- sb_ahaz_learner_params(learner)
  if (identical(params$estimator, "lin_ying")) {
    return(paste0(
      "unpenalized Lin-Ying constant-coefficient additive-hazards fit ",
      "with model-based variance"
    ))
  }
  penalty <- if (identical(params$penalty, "lasso")) {
    "Lasso alpha 1"
  } else {
    "one-step SCAD with a=3.7"
  }
  paste0(
    penalty, ", standardization, 100-lambda path, and support cap ",
    "min(30,p) fixed; lambda selected by BIC inside each analysis fold"
  )
}

sb_ahaz_stop_on_warning <- function(code, stage) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      sb_abort(
        paste0(
          "The ahaz ", stage, " raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_ahaz_warning",
        parent = warning
      )
    }
  )
}

sb_ahaz_wrap_error <- function(code, stage, subclass) {
  tryCatch(
    force(code),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "The ahaz ", stage, " failed: ", conditionMessage(error)
        ),
        subclass,
        parent = error
      )
    }
  )
}

sb_ahaz_engine_x <- function(x) {
  engine_x <- as.matrix(x)
  rownames(engine_x) <- NULL
  colnames(engine_x) <- paste0("x", seq_len(ncol(engine_x)))
  engine_x
}

sb_fit_ahaz_engine <- function(surv, engine_x, learner, params) {
  if (identical(params$estimator, "lin_ying")) {
    return(ahaz::ahaz(
      surv,
      engine_x,
      robust = params$robust
    ))
  }
  penalty <- if (identical(params$penalty, "lasso")) {
    ahaz::lasso.control(alpha = params$alpha)
  } else {
    ahaz::sscad.control(a = params$a, nsteps = params$nsteps)
  }
  ahaz::tune.ahazpen(
    surv,
    engine_x,
    standardize = params$standardize,
    penalty = penalty,
    tune = ahaz::bic.control(),
    dfmax = min(30L, ncol(engine_x)),
    nlambda = params$nlambda
  )
}

sb_ahaz_fit_audit <- function(engine, learner, coefficient) {
  if (identical(learner, "ahaz_lin_ying")) {
    return(list(
      selected_lambda = NA_real_,
      path_length = NA_integer_,
      selected_df = as.integer(sum(coefficient != 0))
    ))
  }
  selected_lambda <- as.numeric(engine$lambda.min)
  path <- as.numeric(engine$lambda)
  valid <- length(selected_lambda) == 1L && !is.na(selected_lambda) &&
    is.finite(selected_lambda) && selected_lambda > 0 &&
    length(path) > 0L && !anyNA(path) && all(is.finite(path)) &&
    all(path > 0) && is.numeric(engine$df) &&
    length(engine$df) == length(path)
  if (!valid) {
    sb_abort(
      paste0("Learner `", learner, "` returned an invalid BIC path."),
      "stormbreaker_invalid_ahaz_fit"
    )
  }
  selected_index <- which.min(abs(path - selected_lambda))
  list(
    selected_lambda = selected_lambda,
    path_length = as.integer(length(path)),
    selected_df = as.integer(engine$df[[selected_index]])
  )
}

sb_validate_ahaz_model <- function(model) {
  valid <- inherits(model, "stormbreaker_ahaz_survival_model") &&
    is.list(model) &&
    is.character(model$learner) && length(model$learner) == 1L &&
    !is.na(model$learner) &&
    model$learner %in% sb_ahaz_survival_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.numeric(model$coefficient) &&
    length(model$coefficient) == length(model$features) &&
    !anyNA(model$coefficient) && all(is.finite(model$coefficient)) &&
    identical(model$params, sb_ahaz_learner_params(model$learner)) &&
    is.list(model$audit)
  if (!valid) {
    sb_abort(
      "The supplied model is not a fitted ahaz survival learner.",
      "stormbreaker_bad_ahaz_fit"
    )
  }
  model
}

sb_ahaz_fit_times <- function(time, event, seed = 10101L) {
  event_times <- time[event == 1L]
  tied_times <- sort(unique(event_times[duplicated(event_times)]))
  fitted_time <- time
  audit <- list(
    method = "seeded_tied_event_offsets_fit_only",
    seed = seed,
    tied_event_groups = length(tied_times),
    tied_events = 0L,
    changed_events = 0L,
    max_abs_shift = 0
  )
  if (!length(tied_times)) return(list(time = time, audit = audit))
  distinct <- sort(unique(c(0, time)))
  withr::with_seed(seed, {
    for (value in tied_times) {
      ii <- which(time == value & event == 1L)
      previous <- distinct[match(value, distinct) - 1L]
      width <- min(value * 1e-8, (value - previous) / 4)
      # Keep deaths at/before their recorded time, so censoring at that time
      # remains after them. Random ranks avoid collisions from random draws.
      offsets <- sample.int(length(ii)) - 1L
      fitted_time[ii] <- value - width * offsets / length(ii)
      if (anyDuplicated(fitted_time[ii]) || any(fitted_time[ii] <= previous)) {
        sb_abort(
          "Tied event times cannot be separated at this numeric precision.",
          "stormbreaker_ahaz_time_precision"
        )
      }
      audit$tied_events <- audit$tied_events + length(ii)
    }
  })
  audit$changed_events <- sum(fitted_time != time)
  audit$max_abs_shift <- max(abs(fitted_time - time))
  list(time = fitted_time, audit = audit)
}

sb_fit_ahaz_survival <- function(x, time, event, learner, params) {
  inputs <- sb_validate_classical_survival_inputs(x, time, event, learner)
  sb_require_ahaz_dependency()
  engine_x <- sb_ahaz_engine_x(inputs$x)
  fit_times <- sb_ahaz_fit_times(inputs$time, inputs$event)
  surv <- survival::Surv(fit_times$time, inputs$event)
  engine <- sb_ahaz_wrap_error(
    sb_ahaz_stop_on_warning(
      sb_fit_ahaz_engine(surv, engine_x, learner, params),
      paste0("fit for `", learner, "`")
    ),
    paste0("fit for `", learner, "`"),
    "stormbreaker_ahaz_fit_failed"
  )
  coefficient <- unname(as.numeric(stats::coef(engine)))
  if (length(coefficient) != ncol(inputs$x) || anyNA(coefficient) ||
    any(!is.finite(coefficient))) {
    sb_abort(
      paste0(
        "Learner `", learner,
        "` did not return one finite coefficient per predictor."
      ),
      "stormbreaker_invalid_ahaz_fit"
    )
  }
  model <- structure(
    list(
      learner = learner,
      features = names(inputs$x),
      params = params,
      coefficient = coefficient,
      audit = c(
        sb_ahaz_fit_audit(engine, learner, coefficient),
        list(time_ties = fit_times$audit)
      )
    ),
    class = "stormbreaker_ahaz_survival_model"
  )
  sb_predict_ahaz_survival(model, inputs$x)
  if (!identical(learner, "ahaz_lin_ying")) model$diagnostic_path <- list(
    beta = engine$fit$beta, lambda = engine$lambda, df = engine$df,
    criterion = engine$tunem, selected_lambda = engine$lambda.min,
    criterion_type = engine$tune$type)
  model
}

sb_validate_ahaz_model_x <- function(model, x) {
  model <- sb_validate_ahaz_model(model)
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted ahaz learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  list(model = model, x = x)
}

sb_predict_ahaz_survival <- function(model, x) {
  validated <- sb_validate_ahaz_model_x(model, x)
  risk <- unname(as.vector(
    as.matrix(validated$x) %*% validated$model$coefficient
  ))
  if (length(risk) != nrow(validated$x) || anyNA(risk) ||
    any(!is.finite(risk))) {
    sb_abort(
      "The ahaz learner produced incomplete or non-finite risk scores.",
      "stormbreaker_invalid_ahaz_prediction"
    )
  }
  risk
}

sb_extract_ahaz_coefficients <- function(model) {
  model <- sb_validate_ahaz_model(model)
  data.frame(
    feature = model$features,
    coefficient = model$coefficient,
    nonzero = model$coefficient != 0,
    direction = ifelse(
      model$coefficient > 0,
      "higher_additive_hazard",
      ifelse(
        model$coefficient < 0,
        "lower_additive_hazard",
        "zero"
      )
    ),
    stringsAsFactors = FALSE
  )
}

sb_ahaz_survival_source <- function(learner) {
  paste0(learner, "_analysis_cox_breslow_calibration")
}

sb_fit_ahaz_calibration <- function(time, event, risk, prediction_horizon) {
  calibration <- sb_ahaz_wrap_error(
    sb_ahaz_stop_on_warning(
      survival::coxph(
        survival::Surv(time, event) ~ risk,
        data = data.frame(time = time, event = event, risk = risk),
        ties = "breslow",
        singular.ok = FALSE,
        model = FALSE,
        x = FALSE,
        y = FALSE
      ),
      "analysis-fold Cox/Breslow calibration"
    ),
    "analysis-fold Cox/Breslow calibration",
    "stormbreaker_ahaz_calibration_failed"
  )
  coefficient <- unname(as.numeric(stats::coef(calibration)))
  if (length(coefficient) != 1L || is.na(coefficient) ||
    !is.finite(coefficient)) {
    sb_abort(
      "The ahaz analysis-fold calibration coefficient is invalid.",
      "stormbreaker_ahaz_calibration_failed"
    )
  }
  if (coefficient < 0) {
    sb_abort(
      paste0(
        "The ahaz analysis-fold calibration reversed the native ",
        "higher-risk direction."
      ),
      "stormbreaker_ahaz_calibration_direction"
    )
  }
  calibrated_risk <- coefficient * risk
  if (anyNA(calibrated_risk) || any(!is.finite(calibrated_risk))) {
    sb_abort(
      "The ahaz calibrated analysis-fold risk is invalid.",
      "stormbreaker_ahaz_calibration_failed"
    )
  }
  list(
    coefficient = coefficient,
    cumulative_baseline_hazard = sb_breslow_baseline_hazard(
      time,
      event,
      calibrated_risk,
      prediction_horizon
    )
  )
}

sb_predict_ahaz_survival_probability <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  model <- sb_validate_ahaz_model(model)
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted ahaz learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_risk <- sb_predict_ahaz_survival(model, inputs$x)
  assessment_risk <- sb_predict_ahaz_survival(model, x_assessment)
  calibration <- sb_fit_ahaz_calibration(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  result <- sb_survival_prediction_result(
    prediction_horizon,
    calibration$cumulative_baseline_hazard,
    source = sb_ahaz_survival_source(model$learner),
    ties = "breslow",
    risk = calibration$coefficient * assessment_risk
  )
  result$calibration_coefficient <- calibration$coefficient
  result
}

sb_extract_ahaz_survival_state <- function(model, context = list()) {
  model <- sb_validate_ahaz_model(model)
  baseline_hazard <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_ahaz_survival_probability(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    baseline_hazard <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      cumulative_baseline_hazard =
        probability$cumulative_baseline_hazard,
      baseline_survival = probability$baseline_survival,
      source = probability$source,
      ties = probability$ties,
      calibration_coefficient = probability$calibration_coefficient
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      coefficient = model$coefficient,
      params = model$params,
      audit = model$audit,
      risk_reference = "ahaz_native",
      risk_formula = "x_beta",
      baseline_hazard = baseline_hazard
    ),
    class = c(
      "stormbreaker_ahaz_survival_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_ahaz_survival_state <- function(state, x) {
  learner <- if (is.list(state)) state$learner else NULL
  if (!inherits(state, "stormbreaker_ahaz_survival_state") ||
    !is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_ahaz_survival_learner_ids()) {
    sb_bad_learner_state("The saved ahaz learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  coefficient <- state$coefficient
  valid_coefficient <- is.numeric(coefficient) &&
    length(coefficient) == ncol(x) && !anyNA(coefficient) &&
    all(is.finite(coefficient))
  if (!valid_coefficient) {
    sb_bad_learner_state("The saved ahaz coefficients are invalid.")
  }
  risk <- unname(as.vector(as.matrix(x) %*% coefficient))
  if (anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "The saved ahaz state produced invalid risk predictions.",
      "stormbreaker_invalid_ahaz_prediction"
    )
  }
  result <- list(risk = risk)
  baseline <- state$baseline_hazard
  if (!is.null(baseline)) {
    calibration <- baseline$calibration_coefficient
    valid <- is.list(baseline) &&
      is.numeric(baseline$prediction_horizon) &&
      length(baseline$prediction_horizon) == 1L &&
      !is.na(baseline$prediction_horizon) &&
      is.finite(baseline$prediction_horizon) &&
      baseline$prediction_horizon > 0 &&
      is.numeric(calibration) && length(calibration) == 1L &&
      !is.na(calibration) && is.finite(calibration) && calibration >= 0 &&
      identical(baseline$source, sb_ahaz_survival_source(learner)) &&
      identical(baseline$ties, "breslow")
    if (!valid) {
      sb_bad_learner_state("The saved ahaz horizon state is invalid.")
    }
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      calibration * risk,
      baseline$cumulative_baseline_hazard
    )
  }
  result
}

sb_ahaz_survival_learner_adapter <- function(
  learner, params = sb_ahaz_learner_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_ahaz_survival(x, time, event, learner, params)
    },
    predict = sb_predict_ahaz_survival,
    predict_survival = sb_predict_ahaz_survival_probability,
    extract_coefficients = sb_extract_ahaz_coefficients,
    extract_state = sb_extract_ahaz_survival_state,
    predict_state = sb_predict_ahaz_survival_state
  )
}

# ---- algorithm-ahaz ----

#' @rdname sb_flexsurv_diagnostics
#' @export
sb_ahaz_diagnostics <- function(fit, type = "coefficients") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% sb_ahaz_survival_learner_ids()) sb_abort(
    "This diagnostic requires an ahaz learner.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("coefficients", "bic_path", "coefficient_path")) sb_abort(
    "Choose coefficients, bic_path or coefficient_path.", "stormbreaker_unsupported_diagnostic")
  if (!is.numeric(s$coefficient) || length(s$coefficient) != length(s$features) || any(!is.finite(s$coefficient))) sb_abort(
    "Saved additive-hazard coefficients are invalid.", "stormbreaker_bad_learner_state")
  path <- NULL
  if (type == "coefficients") {
    tab <- data.frame(feature = s$features, coefficient = s$coefficient)
  } else {
    if (s$learner == "ahaz_lin_ying") sb_abort(
      "Unpenalized Lin-Ying has no penalized BIC path.", "stormbreaker_unsupported_diagnostic")
    path <- fit$diagnostic_artifacts$learner$diagnostic_path
    if (is.null(path)) sb_abort("Original ahaz paths were not captured.", "stormbreaker_missing_diagnostic_artifacts")
    n <- length(path$lambda)
    if (!n || !identical(path$criterion_type, "BIC") || any(!is.finite(path$lambda)) || any(path$lambda <= 0) ||
        length(path$criterion) != n || length(path$df) != n ||
        any(!is.finite(path$criterion)) || any(!is.finite(path$df)) ||
        !identical(dim(path$beta), c(as.integer(length(s$features)), as.integer(n)))) sb_abort(
      "Original additive-hazard path is invalid.", "stormbreaker_invalid_diagnostic")
    selected <- which(path$lambda == path$selected_lambda)
    if (length(selected) != 1L || !isTRUE(all.equal(as.numeric(path$beta[, selected]), s$coefficient, tolerance = 1e-12))) sb_abort(
      "Original selected path does not match the saved model.", "stormbreaker_invalid_diagnostic")
    tab <- if (type == "bic_path") data.frame(lambda = path$lambda, criterion = path$criterion,
      df = path$df, selected = seq_len(n) == selected) else data.frame(
        feature = rep(s$features, each = n), lambda = rep(path$lambda, length(s$features)),
        coefficient = as.numeric(t(as.matrix(path$beta))))
  }
  structure(list(type = paste0("ahaz_", type), table = tab, native_path = path,
    features = s$features, selected_lambda = s$audit$selected_lambda,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "original ahaz additive-hazard coefficients and saved native BIC path",
    interpretation = paste("Original input-unit additive-hazard coefficients, not log hazard ratios or causal effects.",
      "Native BIC is the original selection criterion, not held-out performance; marked lambda is the original choice.",
      "Full path and zero coefficients retained, without refitting or uncertainty reconstruction.",
      "Any separate saved Cox-Breslow probability calibration is not represented by these raw coefficients.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_ahaz_diagnostic <- function(x) {
  feature <- coefficient <- lambda <- criterion <- NULL
  d <- x$table
  if (x$type == "ahaz_coefficients") {
    d$feature <- factor(d$feature, levels = rev(x$features))
    p <- ggplot2::ggplot(d, ggplot2::aes(coefficient, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = .4) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::labs(title = "Additive-hazard coefficients", x = "Coefficient per input unit", y = NULL)
  } else if (x$type == "ahaz_bic_path") {
    p <- ggplot2::ggplot(d, ggplot2::aes(lambda, criterion)) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::geom_point(data = d[d$selected, , drop = FALSE], colour = "#C96942", size = 2.5) +
      ggplot2::scale_x_log10() + ggplot2::labs(title = "Additive-hazard BIC path", x = "Lambda", y = "Native BIC criterion")
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(lambda, coefficient, colour = feature)) +
      ggplot2::geom_line() + ggplot2::geom_vline(xintercept = x$selected_lambda, colour = "grey60", linetype = 2) +
      ggplot2::scale_x_log10() + ggplot2::labs(title = "Additive-hazard coefficient path", x = "Lambda", y = "Coefficient", colour = NULL)
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = if (x$type == "ahaz_coefficients") max(5, 2 + .3 * length(x$features)) else 5)
  p
}

# ---- algorithm-bujar ----

#' @rdname sb_flexsurv_diagnostics
#' @export
sb_bujar_diagnostics <- function(fit, type = "coefficients") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!identical(s$learner, sb_bujar_linear_learner_ids())) sb_abort(
    "This diagnostic requires the Bujar linear learner.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("coefficients", "residual_survival")) sb_abort(
    "Choose coefficients or residual_survival.", "stormbreaker_unsupported_diagnostic")
  if (type == "coefficients") {
    if (length(s$coefficient) != length(s$features) || any(!is.finite(s$coefficient)) ||
        length(s$intercept) != 1L || !is.finite(s$intercept)) sb_abort(
      "Original Bujar coefficients are invalid.", "stormbreaker_bad_learner_state")
    tab <- data.frame(feature = s$features, coefficient = s$coefficient)
  } else {
    if (is.null(s$survival_probability)) sb_abort(
      "No residual distribution was saved: the original fit did not request a prediction horizon.",
      "stormbreaker_unavailable_diagnostic")
    probability <- sb_validate_bujar_probability_state(s)
    tab <- data.frame(residual = probability$residual_time, survival = probability$residual_survival)
  }
  structure(list(type = paste0("bujar_", type), table = tab,
    features = s$features, intercept = s$intercept, learner = s$learner,
    pipeline_id = fit$pipeline$id,
    source = "original Bujar compact coefficients and saved analysis residual Kaplan-Meier mapping",
    interpretation = paste("Coefficients act on predicted log time, not log hazard ratios; risk is negative predicted log time.",
      "The intercept is retained separately. Residual event steps are training-derived, not held-out calibration or patient survival.",
      "No refitting, smoothing, invented iteration history, censored-only steps or uncertainty bands.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_bujar_diagnostic <- function(x) {
  feature <- coefficient <- residual <- survival <- NULL
  d <- x$table
  if (x$type == "bujar_coefficients") {
    d$feature <- factor(d$feature, levels = rev(x$features))
    p <- ggplot2::ggplot(d, ggplot2::aes(coefficient, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = .4) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::labs(title = "Buckley-James coefficients", x = "Log-time coefficient", y = NULL)
  } else p <- ggplot2::ggplot(d, ggplot2::aes(residual, survival)) +
    ggplot2::geom_step(direction = "hv", colour = "#2878A0") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = "Residual survival distribution", x = "Log-time residual", y = "Residual survival")
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = if (x$type == "bujar_coefficients") max(5, 2 + .3 * length(x$features)) else 5)
  p
}

# ---- bujar-linear-adapter ----

sb_bujar_linear_learner_ids <- function() {
  "bujar_linear_boost"
}

sb_bujar_available <- function() {
  requireNamespace("bujar", quietly = TRUE)
}

sb_require_bujar_dependency <- function() {
  if (!sb_bujar_available()) {
    sb_abort(
      "Package `bujar` is required for Buckley-James linear boosting.",
      "stormbreaker_missing_bujar"
    )
  }
  invisible(TRUE)
}

sb_bujar_linear_params <- function(
  learner = sb_bujar_linear_learner_ids()
) {
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !identical(learner, sb_bujar_linear_learner_ids())) {
    sb_abort(
      "The requested Buckley-James learner is not supported.",
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    response_scale = "log_time",
    learner = "linear.regression",
    degree = 1,
    center = TRUE,
    mimpu = NULL,
    iter_bj = 20L,
    max_cycle = 10L,
    nu = 0.1,
    mstop = 50L,
    twin = FALSE,
    mstop2 = 100L,
    tuning = FALSE,
    cv = FALSE,
    vimpint = FALSE,
    rel_inf = FALSE,
    tolerance = 1e-4,
    cores = 1L,
    rng = 1L,
    seed = 1L,
    trace = FALSE
  )
}

sb_bujar_linear_parameter_semantics <- function(
  learner = sb_bujar_linear_learner_ids()
) {
  params <- sb_bujar_linear_params(learner)
  paste0(
    "Buckley-James log-time iterative imputation with componentwise ",
    "linear L2 boosting; ", params$iter_bj, " BJ iterations, ",
    params$mstop, " boosting steps, shrinkage ", params$nu,
    ", tolerance 1e-4, seed 1, and one-core execution fixed before ",
    "resampling"
  )
}

sb_bujar_stop_on_warning <- function(code, stage) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      sb_abort(
        paste0(
          "The Bujar ", stage, " raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_bujar_warning",
        parent = warning
      )
    }
  )
}

sb_bujar_wrap_error <- function(code, stage, subclass) {
  tryCatch(
    force(code),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "The Bujar ", stage, " failed: ", conditionMessage(error)
        ),
        subclass,
        parent = error
      )
    }
  )
}

sb_bujar_engine_x <- sb_ahaz_engine_x

sb_validate_bujar_linear_model <- function(model) {
  valid <- inherits(model, "stormbreaker_bujar_linear_model") &&
    is.list(model) &&
    identical(model$learner, sb_bujar_linear_learner_ids()) &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.numeric(model$intercept) && length(model$intercept) == 1L &&
    !is.na(model$intercept) && is.finite(model$intercept) &&
    is.numeric(model$coefficient) &&
    length(model$coefficient) == length(model$features) &&
    !anyNA(model$coefficient) && all(is.finite(model$coefficient)) &&
    identical(model$params, sb_bujar_linear_params(model$learner))
  if (!valid) {
    sb_abort(
      "The supplied model is not a fitted Bujar linear learner.",
      "stormbreaker_bad_bujar_fit"
    )
  }
  model
}

sb_fit_bujar_linear <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_bujar_dependency()
  engine_x <- sb_bujar_engine_x(inputs$x)
  engine <- sb_bujar_wrap_error(
    withr::with_seed(
      params$seed,
      bujar::bujar(
        y = log(inputs$time),
        cens = inputs$event,
        x = engine_x,
        degree = params$degree,
        learner = params$learner,
        center = params$center,
        mimpu = params$mimpu,
        iter.bj = params$iter_bj,
        max.cycle = params$max_cycle,
        nu = params$nu,
        mstop = params$mstop,
        twin = params$twin,
        mstop2 = params$mstop2,
        tuning = params$tuning,
        cv = params$cv,
        vimpint = params$vimpint,
        rel.inf = params$rel_inf,
        tol = params$tolerance,
        n.cores = params$cores,
        rng = params$rng,
        trace = params$trace
      )
    ),
    "linear-boosting fit",
    "stormbreaker_bujar_fit_failed"
  )
  converged <- is.list(engine) &&
    is.numeric(engine$contype) && length(engine$contype) == 1L &&
    !is.na(engine$contype) && engine$contype == 0 &&
    is.logical(engine$nonconv) && length(engine$nonconv) == 1L &&
    !is.na(engine$nonconv) && !engine$nonconv
  if (!converged) {
    sb_abort(
      paste0(
        "The Bujar linear learner did not converge (contype=",
        paste(engine$contype, collapse = ","), "; nonconv=",
        paste(engine$nonconv, collapse = ","), ")."
      ),
      "stormbreaker_bujar_nonconvergence"
    )
  }
  coefficient <- unname(as.numeric(engine$coef.bj))
  if (length(coefficient) != ncol(engine_x) + 1L ||
    anyNA(coefficient) || any(!is.finite(coefficient))) {
    sb_abort(
      "The Bujar linear learner did not return portable coefficients.",
      "stormbreaker_invalid_bujar_fit"
    )
  }
  intercept <- coefficient[[1L]]
  coefficient <- coefficient[-1L]
  portable <- unname(intercept + as.vector(engine_x %*% coefficient))
  native <- sb_bujar_wrap_error(
    as.numeric(stats::predict(engine, newx = engine_x)),
    "linear-boosting prediction",
    "stormbreaker_bujar_prediction_failed"
  )
  if (!isTRUE(all.equal(
    portable,
    unname(native),
    tolerance = 0,
    check.attributes = FALSE
  ))) {
    sb_abort(
      paste0(
        "Portable Bujar coefficients did not reproduce the native ",
        "training prediction."
      ),
      "stormbreaker_invalid_bujar_fit"
    )
  }
  structure(
    list(
      learner = learner,
      features = names(inputs$x),
      params = params,
      intercept = intercept,
      coefficient = coefficient
    ),
    class = "stormbreaker_bujar_linear_model"
  )
}

sb_validate_bujar_model_x <- function(model, x) {
  model <- sb_validate_bujar_linear_model(model)
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted Bujar learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  list(model = model, x = x)
}

sb_bujar_predicted_log_time <- function(model, x) {
  validated <- sb_validate_bujar_model_x(model, x)
  link <- unname(
    validated$model$intercept +
      as.vector(as.matrix(validated$x) %*% validated$model$coefficient)
  )
  if (length(link) != nrow(validated$x) || anyNA(link) ||
    any(!is.finite(link))) {
    sb_abort(
      "The Bujar linear learner produced invalid predicted log times.",
      "stormbreaker_invalid_bujar_prediction"
    )
  }
  link
}

sb_predict_bujar_linear <- function(model, x) {
  -sb_bujar_predicted_log_time(model, x)
}

sb_bujar_residual_curve <- function(time, event, predicted_log_time) {
  residual <- log(time) - predicted_log_time
  curve <- sb_bujar_wrap_error(
    sb_bujar_stop_on_warning(
      survival::survfit(survival::Surv(residual, event) ~ 1),
      "analysis-fold residual Kaplan-Meier calibration"
    ),
    "analysis-fold residual Kaplan-Meier calibration",
    "stormbreaker_bujar_residual_calibration_failed"
  )
  event_steps <- curve$n.event > 0
  curve_time <- unname(as.numeric(curve$time[event_steps]))
  curve_survival <- unname(as.numeric(curve$surv[event_steps]))
  valid <- length(curve_time) > 0L &&
    length(curve_time) == length(curve_survival) &&
    !anyNA(curve_time) && all(is.finite(curve_time)) &&
    !is.unsorted(curve_time) && !anyNA(curve_survival) &&
    all(is.finite(curve_survival)) &&
    all(curve_survival >= 0 & curve_survival <= 1)
  if (!valid) {
    sb_abort(
      "The Bujar residual Kaplan-Meier curve is invalid.",
      "stormbreaker_bujar_residual_calibration_failed"
    )
  }
  list(time = curve_time, survival = curve_survival)
}

sb_bujar_residual_probability <- function(
  curve_time, curve_survival, predicted_log_time, prediction_horizon
) {
  threshold <- log(prediction_horizon) - predicted_log_time
  index <- findInterval(threshold, curve_time)
  probability <- rep(1, length(index))
  use <- index > 0L
  probability[use] <- curve_survival[index[use]]
  if (anyNA(probability) || any(!is.finite(probability)) ||
    any(probability < 0 | probability > 1)) {
    sb_abort(
      "The Bujar residual calibration produced invalid survival values.",
      "stormbreaker_invalid_bujar_prediction"
    )
  }
  unname(probability)
}

sb_bujar_direct_survival_source <- function() {
  "bujar_linear_boost_analysis_residual_kaplan_meier"
}

sb_predict_bujar_linear_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  model <- sb_validate_bujar_linear_model(model)
  inputs <- sb_validate_additional_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted Bujar learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_link <- sb_bujar_predicted_log_time(model, inputs$x)
  assessment_link <- sb_bujar_predicted_log_time(model, x_assessment)
  curve <- sb_bujar_residual_curve(
    inputs$time,
    inputs$event,
    analysis_link
  )
  list(
    prediction_horizon = prediction_horizon,
    source = sb_bujar_direct_survival_source(),
    ties = "not_applicable",
    probability = sb_bujar_residual_probability(
      curve$time,
      curve$survival,
      assessment_link,
      prediction_horizon
    ),
    residual_time = curve$time,
    residual_survival = curve$survival
  )
}

sb_extract_bujar_linear_state <- function(model, context = list()) {
  model <- sb_validate_bujar_linear_model(model)
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_bujar_linear_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    survival_probability <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      source = probability$source,
      ties = probability$ties,
      residual_time = probability$residual_time,
      residual_survival = probability$residual_survival
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      intercept = model$intercept,
      coefficient = model$coefficient,
      params = model$params,
      risk_reference = "bujar_native",
      risk_formula = "negative_predicted_log_time",
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_bujar_linear_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_validate_bujar_probability_state <- function(state) {
  probability <- state$survival_probability
  valid <- is.list(probability) &&
    is.numeric(probability$prediction_horizon) &&
    length(probability$prediction_horizon) == 1L &&
    !is.na(probability$prediction_horizon) &&
    is.finite(probability$prediction_horizon) &&
    probability$prediction_horizon > 0 &&
    identical(probability$source, sb_bujar_direct_survival_source()) &&
    identical(probability$ties, "not_applicable") &&
    is.numeric(probability$residual_time) &&
    length(probability$residual_time) > 0L &&
    length(probability$residual_time) ==
      length(probability$residual_survival) &&
    !anyNA(probability$residual_time) &&
    all(is.finite(probability$residual_time)) &&
    !is.unsorted(probability$residual_time) &&
    is.numeric(probability$residual_survival) &&
    !anyNA(probability$residual_survival) &&
    all(is.finite(probability$residual_survival)) &&
    all(probability$residual_survival >= 0 &
      probability$residual_survival <= 1)
  if (!valid) {
    sb_bad_learner_state("The saved Bujar horizon state is invalid.")
  }
  probability
}

sb_predict_bujar_linear_state <- function(state, x) {
  learner <- if (is.list(state)) state$learner else NULL
  if (!inherits(state, "stormbreaker_bujar_linear_state") ||
    !identical(learner, sb_bujar_linear_learner_ids())) {
    sb_bad_learner_state("The saved Bujar learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  valid <- is.numeric(state$intercept) && length(state$intercept) == 1L &&
    !is.na(state$intercept) && is.finite(state$intercept) &&
    is.numeric(state$coefficient) &&
    length(state$coefficient) == ncol(x) &&
    !anyNA(state$coefficient) && all(is.finite(state$coefficient))
  if (!valid) {
    sb_bad_learner_state("The saved Bujar coefficients are invalid.")
  }
  link <- unname(
    state$intercept + as.vector(as.matrix(x) %*% state$coefficient)
  )
  if (anyNA(link) || any(!is.finite(link))) {
    sb_abort(
      "The saved Bujar state produced invalid predictions.",
      "stormbreaker_invalid_bujar_prediction"
    )
  }
  result <- list(risk = -link)
  if (!is.null(state$survival_probability)) {
    probability <- sb_validate_bujar_probability_state(state)
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- sb_bujar_residual_probability(
      probability$residual_time,
      probability$residual_survival,
      link,
      probability$prediction_horizon
    )
  }
  result
}

sb_bujar_linear_learner_adapter <- function(
  learner = sb_bujar_linear_learner_ids(),
  params = sb_bujar_linear_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_bujar_linear(x, time, event, learner, params)
    },
    predict = sb_predict_bujar_linear,
    predict_survival = sb_predict_bujar_linear_survival,
    extract_state = sb_extract_bujar_linear_state,
    predict_state = sb_predict_bujar_linear_state
  )
}

# ---- model-timereg ----

# ---- algorithm-timereg ----

#' @rdname sb_flexsurv_diagnostics
#' @export
sb_timereg_diagnostics <- function(fit, type = "effects") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% sb_timereg_survival_learner_ids()) sb_abort(
    "This diagnostic requires an Aalen or proportional-odds timereg fit.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) || !type %in% c("effects", "baseline")) sb_abort(
    "Choose effects or baseline.", "stormbreaker_unsupported_diagnostic")
  aalen <- s$learner == "aalen_additive_survival"
  zero <- as.data.frame(as.list(stats::setNames(rep(0, length(s$features)), s$features)), check.names = FALSE)
  sb_restore_timereg_survival_model(s, zero)
  if (type == "effects" && aalen) {
    tab <- data.frame(time = rep(s$time, length(s$features)),
      feature = rep(s$features, each = length(s$time)),
      value = as.numeric(s$cumulative_coefficient))
    scale <- "Cumulative coefficient"
  } else if (type == "effects") {
    tab <- data.frame(feature = s$features, value = s$coefficient)
    scale <- "Log cumulative-odds coefficient"
  } else {
    tab <- data.frame(time = s$time, value = if (aalen) s$cumulative_intercept else s$baseline_cumulative_odds)
    scale <- if (aalen) "Cumulative intercept" else "Baseline cumulative odds"
  }
  structure(list(type = paste0("timereg_", if (aalen) "aalen_" else "odds_", type),
    table = tab, value_scale = scale, time_unit = fit$development$endpoint$time_unit,
    features = s$features, time = s$time, learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "original compact timereg coefficient and baseline processes",
    interpretation = paste(if (aalen) "Aalen effects are cumulative coefficients, not instantaneous hazard effects or hazard ratios." else
      "Proportional-odds coefficients act on cumulative event odds, not hazards.",
      "All saved times and coefficients retained. Steps are right-continuous; no smoothing, differentiation, clipping or monotonicity correction.",
      "Baseline refers to all predictors equal to zero, not an average or representative patient.",
      "Feature panels use independent vertical scales. No unavailable variance, confidence bands or tests reconstructed; conditional model parameters are not causal effects.",
      if (aalen && type == "baseline" && (any(tab$value < 0) || any(diff(tab$value) < 0)))
        "This unconstrained cumulative-intercept estimate includes negative values or decreases; do not interpret it as a valid monotone patient cumulative-risk curve." else "")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_timereg_diagnostic <- function(x) {
  time <- value <- feature <- NULL
  d <- x$table
  if (x$type == "timereg_odds_effects") {
    d$feature <- factor(d$feature, levels = rev(x$features))
    p <- ggplot2::ggplot(d, ggplot2::aes(value, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = .4) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::labs(title = "Proportional-odds coefficients", x = x$value_scale, y = NULL)
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(time, value)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = .4) +
      ggplot2::geom_step(direction = "hv", colour = "#2878A0", linewidth = .6) +
      ggplot2::labs(title = if (x$type == "timereg_aalen_effects") "Aalen cumulative effects" else x$value_scale,
        x = paste0("Time (", x$time_unit, ")"), y = x$value_scale)
    if (x$type == "timereg_aalen_effects") {
      p$data$feature <- factor(p$data$feature, levels = x$features)
      p <- p + ggplot2::facet_wrap(~feature, ncol = 1, scales = "free_y")
    }
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = if (x$type == "timereg_aalen_effects") max(5, 2.5 * length(x$features)) else 5)
  p
}

# ---- timereg-survival-adapters ----

sb_timereg_survival_learner_ids <- function() {
  c("aalen_additive_survival", "prop_odds_survival")
}

sb_timereg_survival_params <- function(learner) {
  switch(learner,
    aalen_additive_survival = list(
      coefficient_form = "time_varying",
      robust = 0L,
      n_sim = 0L,
      weighted_test = 0L,
      covariance = 0L,
      resample_iid = 0L,
      residuals = 0L,
      seed = 1L,
      risk_time_rule = "last_analysis_event"
    ),
    prop_odds_survival = list(
      iterations = 20L,
      n_sim = 0L,
      weighted_test = 0L,
      profile = 1L,
      sym = 0L,
      baseline_variance = 1L,
      seed = 1L
    ),
    sb_abort(
      paste0("Unknown timereg survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_timereg_available <- function() {
  requireNamespace("timereg", quietly = TRUE)
}

sb_require_timereg <- function(learner) {
  if (!sb_timereg_available()) {
    sb_abort(
      paste0("Package `timereg` is required for learner `", learner, "`."),
      "stormbreaker_missing_timereg_dependency"
    )
  }
  invisible(TRUE)
}

sb_timereg_formula_data <- function(x, time, event, learner) {
  analysis <- x
  time_name <- ".stormbreaker_time"
  while (time_name %in% names(analysis)) {
    time_name <- paste0(time_name, "_")
  }
  event_name <- ".stormbreaker_event"
  while (event_name %in% c(names(analysis), time_name)) {
    event_name <- paste0(event_name, "_")
  }
  analysis[[time_name]] <- time
  analysis[[event_name]] <- event

  response <- if (identical(learner, "aalen_additive_survival")) {
    paste0(
      "survival::Surv(", sb_quote_name(time_name), ", ",
      sb_quote_name(event_name), ")"
    )
  } else {
    paste0(
      "timereg::Event(", sb_quote_name(time_name), ", ",
      sb_quote_name(event_name), ")"
    )
  }
  rhs <- paste(vapply(names(x), sb_quote_name, ""), collapse = " + ")
  list(
    data = analysis,
    formula = stats::as.formula(
      paste(response, "~", rhs),
      env = baseenv()
    )
  )
}

sb_new_timereg_survival_model <- function(learner, features, params, time,
                                           cumulative_intercept = NULL,
                                           cumulative_coefficient = NULL,
                                           risk_reference_time = NULL,
                                           baseline_cumulative_odds = NULL,
                                           coefficient = NULL) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      time = time,
      cumulative_intercept = cumulative_intercept,
      cumulative_coefficient = cumulative_coefficient,
      risk_reference_time = risk_reference_time,
      baseline_cumulative_odds = baseline_cumulative_odds,
      coefficient = coefficient
    ),
    class = "stormbreaker_timereg_survival_fixed"
  )
}

sb_fit_timereg_survival <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_timereg(learner)
  formula_data <- sb_timereg_formula_data(
    inputs$x,
    inputs$time,
    inputs$event,
    learner
  )

  engine <- if (identical(learner, "aalen_additive_survival")) {
    withr::with_seed(
      params$seed,
      timereg::aalen(
        formula = formula_data$formula,
        data = formula_data$data,
        robust = params$robust,
        residuals = params$residuals,
        n.sim = params$n_sim,
        weighted.test = params$weighted_test,
        covariance = params$covariance,
        resample.iid = params$resample_iid,
        silent = 1L
      )
    )
  } else {
    withr::with_seed(
      params$seed,
      timereg::prop.odds(
        formula = formula_data$formula,
        data = formula_data$data,
        Nit = params$iterations,
        detail = 0L,
        n.sim = params$n_sim,
        weighted.test = params$weighted_test,
        profile = params$profile,
        sym = params$sym,
        baselinevar = params$baseline_variance
      )
    )
  }

  cumulative <- as.matrix(engine$cum)
  valid_cumulative <- is.numeric(cumulative) && nrow(cumulative) > 0L &&
    !anyNA(cumulative) && all(is.finite(cumulative)) &&
    all(diff(cumulative[, 1L]) >= 0)
  if (!valid_cumulative) {
    sb_abort(
      paste0("Learner `", learner, "` did not return a finite time process."),
      "stormbreaker_incomplete_timereg_survival_fit"
    )
  }

  if (identical(learner, "aalen_additive_survival")) {
    valid <- ncol(cumulative) == ncol(inputs$x) + 2L
    if (!valid) {
      sb_abort(
        "Aalen fitting did not retain one cumulative coefficient per feature.",
        "stormbreaker_incomplete_timereg_survival_fit"
      )
    }
    coefficient <- cumulative[, -(1:2), drop = FALSE]
    colnames(coefficient) <- names(inputs$x)
    model <- sb_new_timereg_survival_model(
      learner = learner,
      features = names(inputs$x),
      params = params,
      time = unname(cumulative[, 1L]),
      cumulative_intercept = unname(cumulative[, 2L]),
      cumulative_coefficient = unname(coefficient),
      risk_reference_time = unname(cumulative[nrow(cumulative), 1L])
    )
  } else {
    gamma <- as.matrix(engine$gamma)
    valid <- ncol(cumulative) == 2L && nrow(gamma) == ncol(inputs$x) &&
      ncol(gamma) >= 1L && is.numeric(gamma) && !anyNA(gamma[, 1L]) &&
      all(is.finite(gamma[, 1L]))
    if (!valid) {
      sb_abort(
        paste0(
          "Proportional-odds fitting did not retain one finite coefficient ",
          "per feature."
        ),
        "stormbreaker_incomplete_timereg_survival_fit"
      )
    }
    model <- sb_new_timereg_survival_model(
      learner = learner,
      features = names(inputs$x),
      params = params,
      time = unname(cumulative[, 1L]),
      baseline_cumulative_odds = unname(cumulative[, 2L]),
      coefficient = unname(gamma[, 1L])
    )
  }

  sb_predict_timereg_survival(model, inputs$x)
  model
}

sb_validate_timereg_survival_model_x <- function(model, x) {
  valid <- inherits(model, "stormbreaker_timereg_survival_fixed") &&
    is.list(model) && model$learner %in% sb_timereg_survival_learner_ids()
  if (!isTRUE(valid)) {
    sb_abort(
      "The supplied model is not a fitted timereg survival learner.",
      "stormbreaker_bad_timereg_survival_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted timereg learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  x
}

sb_timereg_process_index <- function(time, horizon) {
  eligible <- which(time <= horizon)
  if (length(eligible) == 0L) 0L else max(eligible)
}

sb_predict_timereg_survival <- function(model, x) {
  x <- sb_validate_timereg_survival_model_x(model, x)
  risk <- if (identical(model$learner, "aalen_additive_survival")) {
    row <- length(model$time)
    model$cumulative_intercept[[row]] +
      as.vector(as.matrix(x) %*% model$cumulative_coefficient[row, ])
  } else {
    as.vector(as.matrix(x) %*% model$coefficient)
  }
  risk <- unname(as.numeric(risk))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid risk scores."),
      "stormbreaker_invalid_timereg_survival_prediction"
    )
  }
  risk
}

sb_timereg_direct_survival_source <- function(learner) {
  switch(learner,
    aalen_additive_survival =
      "timereg_aalen_native_additive_survival_direct",
    prop_odds_survival =
      "timereg_proportional_odds_survival_direct",
    NULL
  )
}

sb_timereg_direct_survival_probability <- function(model, x, horizon) {
  x <- sb_validate_timereg_survival_model_x(model, x)
  valid_horizon <- is.numeric(horizon) && length(horizon) == 1L &&
    !is.na(horizon) && is.finite(horizon) && horizon > 0
  if (!valid_horizon) {
    sb_abort(
      "`horizon` must be one finite positive number.",
      "stormbreaker_bad_prediction_horizon"
    )
  }
  index <- sb_timereg_process_index(model$time, horizon)

  probability <- if (identical(
    model$learner,
    "aalen_additive_survival"
  )) {
    cumulative_hazard <- if (index == 0L) {
      rep(0, nrow(x))
    } else {
      model$cumulative_intercept[[index]] +
        as.vector(
          as.matrix(x) %*% model$cumulative_coefficient[index, ]
        )
    }
    exp(-cumulative_hazard)
  } else {
    baseline_odds <- if (index == 0L) {
      0
    } else {
      model$baseline_cumulative_odds[[index]]
    }
    if (!is.finite(baseline_odds) || baseline_odds < 0) {
      rep(NA_real_, nrow(x))
    } else if (baseline_odds == 0) {
      rep(1, nrow(x))
    } else {
      stats::plogis(-(log(baseline_odds) +
        sb_predict_timereg_survival(model, x)))
    }
  }

  probability <- unname(as.numeric(probability))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      paste0(
        "Learner `", model$learner,
        "` produced an invalid native horizon survival probability; ",
        "the estimate was not clipped or recalibrated."
      ),
      "stormbreaker_invalid_timereg_survival_prediction"
    )
  }
  probability
}

sb_predict_timereg_survival_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  inputs <- sb_validate_additional_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted timereg learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  list(
    prediction_horizon = prediction_horizon,
    source = sb_timereg_direct_survival_source(model$learner),
    ties = "not_applicable",
    probability = sb_timereg_direct_survival_probability(
      model,
      x_assessment,
      prediction_horizon
    )
  )
}

sb_extract_timereg_survival_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_timereg_survival_fixed")) {
    sb_bad_learner_state("The fitted timereg survival learner is invalid.")
  }
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_timereg_survival_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x,
      context$prediction_horizon
    )
    survival_probability <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      source = probability$source,
      ties = probability$ties
    )
  }

  state <- unclass(model)
  state$survival_probability <- survival_probability
  state$risk_reference <- if (identical(
    model$learner,
    "aalen_additive_survival"
  )) {
    "timereg_native_last_analysis_event_time"
  } else {
    "zero"
  }
  state$risk_formula <- if (identical(
    model$learner,
    "aalen_additive_survival"
  )) {
    "cumulative_hazard_at_last_analysis_event"
  } else {
    "x_beta"
  }
  structure(
    state,
    class = c(
      "stormbreaker_timereg_survival_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_timereg_survival_model <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_timereg_survival_learner_ids()) {
    sb_bad_learner_state("The saved timereg survival state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  valid_time <- is.numeric(state$time) && length(state$time) > 0L &&
    !anyNA(state$time) && all(is.finite(state$time)) &&
    all(diff(state$time) >= 0)
  if (!valid_time) {
    sb_bad_learner_state("The saved timereg time process is invalid.")
  }

  if (identical(learner, "aalen_additive_survival")) {
    valid <- is.numeric(state$cumulative_intercept) &&
      length(state$cumulative_intercept) == length(state$time) &&
      !anyNA(state$cumulative_intercept) &&
      all(is.finite(state$cumulative_intercept)) &&
      is.matrix(state$cumulative_coefficient) &&
      nrow(state$cumulative_coefficient) == length(state$time) &&
      ncol(state$cumulative_coefficient) == ncol(x) &&
      !anyNA(state$cumulative_coefficient) &&
      all(is.finite(state$cumulative_coefficient)) &&
      is.numeric(state$risk_reference_time) &&
      length(state$risk_reference_time) == 1L &&
      identical(
        state$risk_reference_time,
        state$time[[length(state$time)]]
      )
  } else {
    valid <- is.numeric(state$baseline_cumulative_odds) &&
      length(state$baseline_cumulative_odds) == length(state$time) &&
      !anyNA(state$baseline_cumulative_odds) &&
      all(is.finite(state$baseline_cumulative_odds)) &&
      is.numeric(state$coefficient) && length(state$coefficient) == ncol(x) &&
      !anyNA(state$coefficient) && all(is.finite(state$coefficient))
  }
  if (!isTRUE(valid)) {
    sb_bad_learner_state("The saved timereg prediction parameters are invalid.")
  }

  model <- sb_new_timereg_survival_model(
    learner = learner,
    features = state$features,
    params = state$params,
    time = state$time,
    cumulative_intercept = state$cumulative_intercept,
    cumulative_coefficient = state$cumulative_coefficient,
    risk_reference_time = state$risk_reference_time,
    baseline_cumulative_odds = state$baseline_cumulative_odds,
    coefficient = state$coefficient
  )
  list(model = model, x = x)
}

sb_predict_timereg_survival_state <- function(state, x) {
  restored <- sb_restore_timereg_survival_model(state, x)
  result <- list(
    risk = sb_predict_timereg_survival(restored$model, restored$x)
  )
  probability <- state$survival_probability
  if (!is.null(probability)) {
    valid <- is.list(probability) &&
      is.numeric(probability$prediction_horizon) &&
      length(probability$prediction_horizon) == 1L &&
      !is.na(probability$prediction_horizon) &&
      is.finite(probability$prediction_horizon) &&
      probability$prediction_horizon > 0 &&
      is.character(probability$time_unit) &&
      length(probability$time_unit) == 1L && !is.na(probability$time_unit) &&
      nzchar(probability$time_unit) &&
      identical(
        probability$source,
        sb_timereg_direct_survival_source(state$learner)
      ) && identical(probability$ties, "not_applicable")
    if (!valid) {
      sb_bad_learner_state(
        "The saved timereg survival-probability state is invalid."
      )
    }
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- sb_timereg_direct_survival_probability(
      restored$model,
      restored$x,
      probability$prediction_horizon
    )
  }
  result
}

sb_timereg_survival_learner_adapter <- function(
  learner,
  params = sb_timereg_survival_params(learner)
) {
  if (!learner %in% sb_timereg_survival_learner_ids()) {
    sb_abort(
      paste0("Unknown timereg survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_timereg_survival(x, time, event, learner, params)
    },
    predict = sb_predict_timereg_survival,
    predict_survival = sb_predict_timereg_survival_survival,
    extract_state = sb_extract_timereg_survival_state,
    predict_state = sb_predict_timereg_survival_state
  )
}

# ---- model-mtlr ----

# ---- algorithm-mtlr ----

#' @rdname sb_flexsurv_diagnostics
#' @export
sb_mtlr_diagnostics <- function(fit, type = "feature_weights") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (fit$pipeline$learner != "mtlr_survival") sb_abort(
    "This diagnostic requires an MTLR final fit.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("feature_weights", "biases")) sb_abort(
    "Choose feature_weights or biases.", "stormbreaker_unsupported_diagnostic")
  s <- fit$learner_state
  sb_restore_mtlr_survival_model(s)
  e <- s$engine
  w <- e$weight_matrix
  times <- e$time_points
  if (!is.matrix(w) || nrow(w) != length(times) ||
      !identical(colnames(w), c("Bias", s$internal_features)) || any(!is.finite(w)) ||
      !length(times) || any(!is.finite(times)) || any(diff(times) <= 0)) sb_abort(
    "Original MTLR weight/time mapping is invalid.", "stormbreaker_bad_learner_state")
  center <- rep(0, length(s$features))
  sd <- rep(1, length(s$features))
  if (!is.null(e$scale)) {
    center <- unname(e$scale$center)
    sd <- unname(e$scale$sd)
    if (length(center) != length(s$features) || length(sd) != length(s$features) ||
        any(!is.finite(center)) || any(!is.finite(sd)) || any(sd < 0)) sb_abort(
      "Original MTLR normalization is invalid.", "stormbreaker_bad_learner_state")
  }
  slopes <- sweep(w[, -1, drop = FALSE], 2, ifelse(sd == 0, 1, sd), "/")
  slopes[, sd == 0] <- 0
  raw <- cbind(w[, 1] - drop(slopes %*% center), slopes)
  dimnames(raw) <- dimnames(w)
  tab <- if (type == "feature_weights") data.frame(
    time_index = rep(seq_along(times), length(s$features)),
    time = rep(times, length(s$features)), feature = rep(s$features, each = length(times)),
    weight = as.numeric(w[, -1, drop = FALSE])) else data.frame(
      time_index = seq_along(times), time = times, bias = as.numeric(w[, 1]))
  structure(list(type = paste0("mtlr_", type), table = tab, native_weights = w,
    original_unit_weights = raw, time_points = times,
    normalization = data.frame(feature = s$features, center = center, sd = sd),
    requested_intervals = s$requested_intervals, effective_time_points = s$effective_time_points,
    time_unit = fit$development$endpoint$time_unit, learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "original compact MTLR effective time grid, weights and normalization",
    interpretation = paste("Native normalized-input feature weights and biases are separate; all effective time points retained.",
      "Columns represent ordered native time points, not equal-length time intervals.",
      "MTLR couples these potentials to produce a survival distribution; individual weights are not hazard ratios, independent logits or SHAP.",
      "Exact input-unit equivalents are retained separately; no refitting, optimization history or uncertainty is added.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_mtlr_diagnostic <- function(x) {
  time_index <- feature <- weight <- time <- bias <- NULL
  d <- x$table
  if (x$type == "mtlr_feature_weights") {
    d$feature <- factor(d$feature, levels = rev(x$normalization$feature))
    ticks <- unique(round(seq(1, length(x$time_points), length.out = min(8, length(x$time_points)))))
    bound <- max(abs(d$weight))
    if (bound == 0) bound <- 1
    p <- ggplot2::ggplot(d, ggplot2::aes(time_index, feature, fill = weight)) +
      ggplot2::geom_tile(colour = "white", linewidth = .2) +
      ggplot2::scale_x_continuous(breaks = ticks, labels = format(signif(x$time_points[ticks], 3), trim = TRUE), expand = c(0, 0)) +
      ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#C96942", limits = c(-bound, bound),
        guide = ggplot2::guide_colourbar(display = "rectangles")) +
      ggplot2::labs(title = "MTLR feature weights", x = paste0("Ordered time points (", x$time_unit, ")"), y = NULL, fill = "Weight")
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(time, bias)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = .4) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::labs(title = "MTLR biases", x = paste0("Time (", x$time_unit, ")"), y = "Native bias")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = if (x$type == "mtlr_feature_weights") max(5, 2 + .3 * nrow(x$normalization)) else 5)
  p
}

# ---- mtlr-survival-adapter ----

sb_mtlr_survival_learner_ids <- function() {
  "mtlr_survival"
}

sb_mtlr_survival_params <- function(learner = sb_mtlr_survival_learner_ids()) {
  if (!identical(learner, sb_mtlr_survival_learner_ids())) {
    sb_abort(
      paste0("Unknown MTLR survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    interval_rule = "ceiling_sqrt_analysis_n",
    normalize = TRUE,
    C1 = 1,
    train_biases = TRUE,
    train_uncensored = TRUE,
    threshold = 1e-5,
    maximum_iterations = 5000L,
    lower_weight_bound = -15,
    upper_weight_bound = 15
  )
}

sb_mtlr_survival_parameter_semantics <- function() {
  paste0(
    "multi-task logistic regression with ceiling-sqrt analysis-fold time ",
    "intervals, normalization, C1 1, and bounded optimization fixed before ",
    "resampling"
  )
}

sb_mtlr_survival_available <- function() {
  requireNamespace("MTLR", quietly = TRUE)
}

sb_require_mtlr_survival_dependency <- function() {
  if (!sb_mtlr_survival_available()) {
    sb_abort(
      "Package `MTLR` is required for learner `mtlr_survival`.",
      "stormbreaker_missing_mtlr_survival_dependency"
    )
  }
  invisible(TRUE)
}

sb_mtlr_internal_features <- function(p) {
  paste0(".stormbreaker_mtlr_x", seq_len(p))
}

sb_mtlr_formula_data <- function(x, time, event) {
  internal_features <- sb_mtlr_internal_features(ncol(x))
  analysis <- x
  names(analysis) <- internal_features
  time_name <- ".stormbreaker_mtlr_time"
  event_name <- ".stormbreaker_mtlr_event"
  analysis[[time_name]] <- time
  analysis[[event_name]] <- event
  response <- paste0(
    "survival::Surv(", time_name, ", ", event_name, ")"
  )
  list(
    data = analysis,
    formula = stats::reformulate(
      internal_features,
      response = response,
      env = baseenv()
    ),
    internal_features = internal_features
  )
}

sb_mtlr_prediction_engine <- function(engine) {
  engine[c("x", "y", "response", "Call")] <- NULL
  attr(engine$Terms, ".Environment") <- baseenv()
  engine
}

sb_new_mtlr_survival_model <- function(
  features, internal_features, params, requested_intervals,
  effective_time_points, engine
) {
  structure(
    list(
      learner = sb_mtlr_survival_learner_ids(),
      features = features,
      internal_features = internal_features,
      params = params,
      requested_intervals = requested_intervals,
      effective_time_points = effective_time_points,
      engine = engine
    ),
    class = "stormbreaker_mtlr_survival_model"
  )
}

sb_fit_mtlr_survival <- function(x, time, event, learner, params) {
  inputs <- sb_validate_classical_survival_inputs(x, time, event, learner)
  sb_require_mtlr_survival_dependency()
  formula_data <- sb_mtlr_formula_data(
    inputs$x,
    inputs$time,
    inputs$event
  )
  intervals <- as.integer(ceiling(sqrt(nrow(inputs$x))))
  engine <- MTLR::mtlr(
    formula = formula_data$formula,
    data = formula_data$data,
    nintervals = intervals,
    normalize = params$normalize,
    C1 = params$C1,
    train_biases = params$train_biases,
    train_uncensored = params$train_uncensored,
    threshold = params$threshold,
    maxit = params$maximum_iterations,
    lower = params$lower_weight_bound,
    upper = params$upper_weight_bound
  )
  engine <- sb_mtlr_prediction_engine(engine)
  valid <- inherits(engine, "mtlr") &&
    is.matrix(engine$weight_matrix) &&
    nrow(engine$weight_matrix) == length(engine$time_points) &&
    ncol(engine$weight_matrix) == length(formula_data$internal_features) + 1L &&
    !anyNA(engine$weight_matrix) && all(is.finite(engine$weight_matrix)) &&
    is.numeric(engine$time_points) && length(engine$time_points) >= 1L &&
    !anyNA(engine$time_points) && all(is.finite(engine$time_points))
  if (!valid) {
    sb_abort(
      "MTLR did not return finite weights and a survival time grid.",
      "stormbreaker_incomplete_mtlr_survival_fit"
    )
  }
  sb_new_mtlr_survival_model(
    names(inputs$x),
    formula_data$internal_features,
    params,
    intervals,
    length(engine$time_points),
    engine
  )
}

sb_mtlr_prediction_data <- function(model, x) {
  if (!inherits(model, "stormbreaker_mtlr_survival_model")) {
    sb_abort(
      "The supplied model is not a fitted MTLR survival learner.",
      "stormbreaker_bad_mtlr_survival_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted MTLR learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  names(x) <- model$internal_features
  x
}

sb_predict_mtlr_survival_risk <- function(model, x) {
  sb_require_mtlr_survival_dependency()
  prediction_data <- sb_mtlr_prediction_data(model, x)
  mean_time <- as.numeric(stats::predict(
    model$engine,
    newdata = prediction_data,
    type = "mean_time"
  ))
  risk <- -mean_time
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "MTLR produced invalid native mean-survival risk scores.",
      "stormbreaker_invalid_mtlr_survival_prediction"
    )
  }
  unname(risk)
}

sb_mtlr_survival_probability <- function(model, x, prediction_horizon) {
  sb_require_mtlr_survival_dependency()
  prediction_data <- sb_mtlr_prediction_data(model, x)
  curves <- stats::predict(
    model$engine,
    newdata = prediction_data,
    type = "prob_times",
    times = prediction_horizon
  )
  probability <- if (is.data.frame(curves) && nrow(curves) >= 1L) {
    unname(as.numeric(curves[nrow(curves), -1L, drop = TRUE]))
  } else {
    numeric()
  }
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      "MTLR produced invalid native survival probabilities.",
      "stormbreaker_invalid_mtlr_survival_probability"
    )
  }
  probability
}

sb_mtlr_direct_survival_source <- function() {
  "MTLR_native_multitask_survival_curve"
}

sb_predict_mtlr_survival_probability <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted MTLR learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  list(
    prediction_horizon = prediction_horizon,
    source = sb_mtlr_direct_survival_source(),
    ties = "not_applicable",
    probability = sb_mtlr_survival_probability(
      model,
      x_assessment,
      prediction_horizon
    )
  )
}

sb_extract_mtlr_survival_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_mtlr_survival_model")) {
    sb_bad_learner_state("The fitted MTLR survival learner state is invalid.")
  }
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    horizon <- sb_validate_prediction_horizon(
      context$prediction_horizon,
      context$time
    )
    survival_probability <- list(
      prediction_horizon = horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      source = sb_mtlr_direct_survival_source(),
      ties = "not_applicable"
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      internal_features = model$internal_features,
      params = model$params,
      requested_intervals = model$requested_intervals,
      effective_time_points = model$effective_time_points,
      engine = model$engine,
      training_data_retention = "none",
      risk_reference = "MTLR_native_time_grid",
      risk_formula = "negative_native_mean_predicted_survival_time",
      baseline_hazard = NULL,
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_mtlr_survival_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_mtlr_survival_model <- function(state) {
  valid <- inherits(state, "stormbreaker_mtlr_survival_state") &&
    identical(state$learner, sb_mtlr_survival_learner_ids()) &&
    is.character(state$internal_features) &&
    length(state$internal_features) == length(state$features) &&
    is.list(state$params) && inherits(state$engine, "mtlr") &&
    identical(state$training_data_retention, "none")
  if (!isTRUE(valid)) {
    sb_bad_learner_state("The saved MTLR survival learner state is invalid.")
  }
  sb_new_mtlr_survival_model(
    state$features,
    state$internal_features,
    state$params,
    state$requested_intervals,
    state$effective_time_points,
    state$engine
  )
}

sb_predict_mtlr_survival_state <- function(state, x) {
  x <- sb_validate_learner_state_x(
    state,
    x,
    sb_mtlr_survival_learner_ids()
  )
  model <- sb_restore_mtlr_survival_model(state)
  result <- list(risk = sb_predict_mtlr_survival_risk(model, x))
  if (!is.null(state$survival_probability)) {
    probability <- state$survival_probability
    valid <- is.numeric(probability$prediction_horizon) &&
      length(probability$prediction_horizon) == 1L &&
      !is.na(probability$prediction_horizon) &&
      is.finite(probability$prediction_horizon) &&
      probability$prediction_horizon > 0 &&
      is.character(probability$time_unit) &&
      length(probability$time_unit) == 1L &&
      !is.na(probability$time_unit) && nzchar(probability$time_unit) &&
      identical(probability$source, sb_mtlr_direct_survival_source()) &&
      identical(probability$ties, "not_applicable")
    if (!isTRUE(valid)) {
      sb_bad_learner_state(
        "The saved MTLR direct survival-probability state is invalid."
      )
    }
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- sb_mtlr_survival_probability(
      model,
      x,
      probability$prediction_horizon
    )
  }
  result
}

sb_mtlr_survival_learner_adapter <- function(
  learner = sb_mtlr_survival_learner_ids(),
  params = sb_mtlr_survival_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_mtlr_survival(x, time, event, learner, params)
    },
    predict = sb_predict_mtlr_survival_risk,
    predict_survival = sb_predict_mtlr_survival_probability,
    extract_state = sb_extract_mtlr_survival_state,
    predict_state = sb_predict_mtlr_survival_state
  )
}
