#' Survival-model diagnostics
#' @name sb_flexsurv_diagnostics
#' @description Inspect fitted survival-model parameters, effects and native tuning diagnostics.
#' @rawRd \arguments{
#' \item{fit}{
#' \strong{For \code{sb_flexsurv_diagnostics()}:} A successful fixed flexsurv or flexsurvspline final fit.
#' 
#' \strong{For \code{sb_aft_diagnostics()}:} A successful Weibull, exponential, log-normal or log-logistic
#' survreg final fit.
#' 
#' \strong{For \code{sb_ahaz_diagnostics()}:} A successful Lin-Ying, LASSO or SSCAD ahaz final fit.
#' 
#' \strong{For \code{sb_bujar_diagnostics()}:} A successful Bujar linear final fit.
#' 
#' \strong{For \code{sb_penaft_diagnostics()}:} A successful penAFT final fit.
#' 
#' \strong{For \code{sb_timereg_diagnostics()}:} A successful timereg Aalen or proportional-odds final fit.
#' 
#' \strong{For \code{sb_mtlr_diagnostics()}:} A successful MTLR final fit.
#' }
#' \item{type}{
#' \strong{For \code{sb_flexsurv_diagnostics()}:} location_coefficients, or spline_baseline for spline models.
#' 
#' \strong{For \code{sb_aft_diagnostics()}:} log_time_coefficients or time_ratios.
#' 
#' \strong{For \code{sb_ahaz_diagnostics()}:} coefficients, bic_path or coefficient_path. Paths require
#' capture_diagnostics = TRUE during the original penalized final fit.
#' 
#' \strong{For \code{sb_bujar_diagnostics()}:} coefficients or residual_survival. The latter requires an
#' originally requested prediction horizon and its saved residual curve.
#' 
#' \strong{For \code{sb_penaft_diagnostics()}:} coefficients, cv_curve or coefficient_path. Paths require
#' capture_diagnostics = TRUE during the original final fit.
#' 
#' \strong{For \code{sb_timereg_diagnostics()}:} effects or baseline.
#' 
#' \strong{For \code{sb_mtlr_diagnostics()}:} feature_weights or biases.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_flexsurv_diagnostics()}}{
#' Inspect original flexible parametric survival parameters
#' 
#' \strong{Returns:} Original regression coefficients on their actual parameter scale,
#' or the saved spline's survival curve with all input predictors set to zero.
#' This mathematical baseline is not an average or representative patient and
#' zero can be outside the training predictor range. The curve is confined to
#' saved boundary knots; internal knot positions are shown without refitting.
#' Compact parameters suffice. No hazard-ratio conversion or uncertainty is
#' invented. Distribution parameters and full-precision knots remain separate.
#' }
#' 
#' \subsection{\code{sb_aft_diagnostics()}}{
#' Inspect original accelerated failure-time coefficients
#' 
#' \strong{Returns:} An algorithm diagnostic with original feature coefficients and
#' exp(coefficient) time ratios per input unit. Intercept and distribution scale
#' are retained separately. A ratio above one multiplies conditional survival
#' time quantiles upward; it is not a hazard ratio or causal treatment effect.
#' Only point estimates are provided; missing covariance is not reconstructed.
#' Compact saved parameters suffice and no model is refitted.
#' }
#' 
#' \subsection{\code{sb_ahaz_diagnostics()}}{
#' Inspect original constant-coefficient additive-hazard models
#' 
#' \strong{Returns:} Original input-unit additive-hazard coefficients or native penalized
#' coefficient/BIC paths. All path points, zero coefficients and original
#' selected lambda are retained. A separately saved Cox-Breslow calibration
#' is not the raw additive-hazard equation and is not plotted here. No missing
#' path, uncertainty, baseline or likelihood is reconstructed.
#' }
#' 
#' \subsection{\code{sb_bujar_diagnostics()}}{
#' Inspect saved Buckley-James linear boosting effects
#' 
#' \strong{Returns:} Original log-time coefficients, or the saved event-step residual
#' Kaplan-Meier distribution used by the model's probability mapping.
#' The intercept is retained separately. Coefficients are not log hazard
#' ratios; the package risk score is negative predicted log time. The residual
#' distribution is training-derived, not held-out calibration or a patient
#' survival curve. No variance, censored-only steps or training iteration
#' history is reconstructed from compact state.
#' }
#' 
#' \subsection{\code{sb_penaft_diagnostics()}}{
#' Inspect original Gehan elastic-net AFT coefficients and tuning
#' 
#' \strong{Returns:} Original input-unit log-time coefficients and centered intercept,
#' or the complete original path and linear-predictor CV selection score.
#' Standardized native path coefficients are converted using saved means/SDs;
#' the native path remains available separately. The CV score is not AUC or
#' held-out performance. No standard errors, surrogate fit or missing path
#' is reconstructed. Package risk is negative centered log-time prediction.
#' }
#' 
#' \subsection{\code{sb_timereg_diagnostics()}}{
#' Inspect original additive-hazard and proportional-odds processes
#' 
#' \strong{Returns:} A diagnostic preserving every saved time and cumulative coefficient
#' for Aalen, or original log-odds coefficients for proportional odds. Baselines
#' retain the original cumulative intercept or cumulative odds, respectively.
#' Steps are not smoothed or differentiated; negative/nonmonotone estimates
#' are retained. No confidence bands or resampling tests are reconstructed.
#' The zero-predictor baseline is not a representative-patient prediction.
#' }
#' 
#' \subsection{\code{sb_mtlr_diagnostics()}}{
#' Inspect the original multi-task logistic survival weights
#' 
#' \strong{Returns:} An algorithm diagnostic preserving every effective native time point
#' and original weight. Feature weights act on the native normalized inputs;
#' biases are displayed separately. These coupled interval potentials are not
#' independent logistic models, hazard ratios or feature-attribution scores.
#' Exact original-input-unit weights are retained separately for inspection,
#' with centering absorbed into each bias and native zero-SD behavior retained.
#' Neither an optimization trace nor uncertainty is reconstructed.
#' }
#' }
NULL

# ---- model-flexsurv ----

# ---- algorithm-flexsurv ----

#' @rdname sb_flexsurv_diagnostics
#' @export
sb_flexsurv_diagnostics <- function(fit, type = "location_coefficients") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  spline <- s$learner %in% sb_flexsurvspline_learner_ids()
  if (!spline && !s$learner %in% sb_flexsurv_parametric_learner_ids()) sb_abort(
    "This diagnostic requires a fixed flexsurv or flexsurvspline model.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("location_coefficients", "spline_baseline") || (type == "spline_baseline" && !spline)) sb_abort(
    "Choose location_coefficients, or spline_baseline for a spline model.", "stormbreaker_unsupported_diagnostic")
  if (spline) sb_restore_flexsurvspline_model(s, s$learner) else sb_restore_flexsurv_parametric_model(s, s$learner)
  scale <- if (spline) switch(s$params$scale, hazard = "Log cumulative hazard",
    odds = "Log cumulative odds", normal = "Probit event probability") else
    if (s$params$distribution == "gengamma") "Log-time location (mu)" else "Log rate"
  if (type == "location_coefficients") {
    tab <- data.frame(feature = s$features, coefficient = s$coefficient)
    interpretation <- paste("Original location-regression coefficients per input unit on the stated native scale.",
      "All fitted features retained; ancillary parameters and intercept/baseline coefficients remain separate.",
      "Not a universal hazard ratio, SHAP importance, confidence interval or causal effect.")
  } else {
    times <- sort(unique(c(exp(seq(min(s$knots), max(s$knots), length.out = 201)), exp(s$knots))))
    probability <- flexsurv::psurvspline(times, gamma = s$gamma, knots = s$knots,
      scale = s$params$scale, timescale = s$params$timescale, spline = s$params$spline, lower.tail = FALSE)
    tab <- data.frame(time = times, survival = as.numeric(probability))
    interpretation <- paste("Original spline baseline with every input predictor set to zero, not a typical patient or cohort average.",
      "Zero may lie outside the training predictor range. Time is restricted to the saved boundary knots.",
      "Vertical markers are original internal knots. No uncertainty or clinical validation is implied.")
  }
  structure(list(type = paste0("flexsurv_", type), table = tab, parameter_scale = scale,
    intercept = s$intercept, ancillary = s$ancillary, gamma = s$gamma, knots = s$knots,
    params = s$params, learner = s$learner, pipeline_id = fit$pipeline$id,
    time_unit = fit$development$endpoint$time_unit,
    source = "original compact flexsurv prediction parameters", interpretation = interpretation),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_flexsurv_diagnostic <- function(x) {
  feature <- coefficient <- time <- survival <- NULL
  d <- x$table
  if (x$type == "flexsurv_location_coefficients") {
    d$feature <- factor(d$feature, levels = rev(d$feature))
    p <- ggplot2::ggplot(d, ggplot2::aes(coefficient, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey65", linewidth = .4) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::labs(title = "Location coefficients", x = paste0(x$parameter_scale, " coefficient"), y = NULL)
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(time, survival)) +
      ggplot2::geom_vline(xintercept = exp(x$knots[-c(1, length(x$knots))]), colour = "grey65", linetype = 2) +
      ggplot2::geom_line(colour = "#2878A0", linewidth = .7) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = "Spline baseline",
        x = if (length(x$time_unit)) paste0("Time (", x$time_unit, ")") else "Time (model units)", y = "Survival probability")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 10, height = 5)
  p
}

# ---- flexsurv-numerics ----

# Numerical coordinates change neither the declared model nor the feature set.
# Compact prediction parameters are mapped back to original input units.
sb_flexsurv_numerical_settings <- function(x) {
  center <- colMeans(x)
  scale <- vapply(x, stats::sd, numeric(1))
  if (any(!is.finite(scale) | scale <= 0)) {
    sb_abort("Flexible survival fitting requires positive finite training predictor scales.",
      "stormbreaker_flexsurv_invalid_scale")
  }
  list(center = center, scale = scale, maxit = 1000L,
    standardization = "training_center_sd", initialization = "native")
}

sb_flexsurv_numerical_predictors <- function(x, settings) {
  as.data.frame(sweep(sweep(as.matrix(x), 2L, settings$center, "-"),
    2L, settings$scale, "/"), check.names = FALSE)
}

sb_flexsurvspline_numerical_fit <- function(formula, data, k, scale,
                                           timescale, spline, hessian, maxit) {
  baseline <- flexsurv::flexsurvspline(stats::update(formula, . ~ 1),
    data = data, k = k, scale = scale, timescale = timescale, spline = spline,
    hessian = TRUE, control = list(maxit = maxit))
  if (!identical(baseline$opt$convergence, 0L) || any(!is.finite(stats::coef(baseline)))) {
    sb_abort("The training-fold intercept-only spline initialization did not converge.",
      "stormbreaker_flexsurvspline_nonconvergence")
  }
  mm <- stats::model.matrix(formula, stats::model.frame(formula, data = data))
  init <- c(unname(stats::coef(baseline)), rep(0, ncol(mm) - 1L))
  fit <- flexsurv::flexsurvspline(formula, data = data, k = k, scale = scale,
    timescale = timescale, spline = spline, hessian = hessian, inits = init,
    control = list(maxit = maxit))
  if (!isTRUE(all.equal(baseline$knots, fit$knots))) {
    sb_abort("Spline initialization changed the training-fold knots.",
      "stormbreaker_flexsurvspline_fit_failure")
  }
  fit
}

# ---- flexsurv-parametric-fixed-adapters ----

sb_flexsurv_parametric_learner_ids <- function() {
  c(
    "flexsurv_gamma",
    "flexsurv_gompertz",
    "flexsurv_gengamma"
  )
}

sb_flexsurv_parametric_params <- function(learner) {
  switch(learner,
    flexsurv_gamma = list(
      distribution = "gamma",
      parameterization = "shape_rate",
      location_parameter = "rate",
      risk_transform = "log_location"
    ),
    flexsurv_gompertz = list(
      distribution = "gompertz",
      parameterization = "shape_rate",
      location_parameter = "rate",
      risk_transform = "log_location"
    ),
    flexsurv_gengamma = list(
      distribution = "gengamma",
      parameterization = "stable_mu_sigma_Q",
      location_parameter = "mu",
      risk_transform = "negative_location"
    ),
    sb_abort(
      paste0("Unknown fixed flexsurv learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_flexsurv_parametric_parameter_semantics <- function(learner) {
  params <- sb_flexsurv_parametric_params(learner)
  switch(params$distribution,
    gamma = paste0(
      "Gamma shape-rate distribution and native rate regression fixed; ",
      "parameters estimated in each analysis fold"
    ),
    gompertz = paste0(
      "Gompertz shape-rate proportional-hazards distribution fixed; ",
      "parameters estimated in each analysis fold"
    ),
    gengamma = paste0(
      "stable generalized-gamma mu-sigma-Q parameterization fixed; ",
      "parameters estimated in each analysis fold"
    )
  )
}

sb_flexsurv_parametric_dependency <- function(learner) {
  if (!learner %in% sb_flexsurv_parametric_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed flexsurv learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    package = "flexsurv",
    packages = "flexsurv",
    available = requireNamespace("flexsurv", quietly = TRUE)
  )
}

sb_require_flexsurv_parametric_dependency <- function(learner) {
  dependency <- sb_flexsurv_parametric_dependency(learner)
  if (!dependency$available) {
    sb_abort(
      paste0("Package `flexsurv` is required for learner `", learner, "`."),
      "stormbreaker_missing_flexsurv_dependency"
    )
  }
  invisible(dependency)
}

sb_flexsurv_parametric_ancillary_names <- function(distribution) {
  switch(distribution,
    gamma = "shape",
    gompertz = "shape",
    gengamma = c("sigma", "Q")
  )
}

sb_new_flexsurv_parametric_model <- function(
  learner, features, engine_features, params, prediction, engine = NULL
) {
  structure(
    list(
      learner = learner,
      features = features,
      engine_features = engine_features,
      params = params,
      prediction = prediction,
      engine = engine
    ),
    class = "stormbreaker_flexsurv_parametric_fixed"
  )
}

sb_fit_flexsurv_parametric <- function(
  x, time, event, learner, params
) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_flexsurv_parametric_dependency(learner)

  numerical <- sb_flexsurv_numerical_settings(inputs$x)
  engine_x <- sb_flexsurv_numerical_predictors(inputs$x, numerical)
  engine_features <- sprintf(".stormbreaker_x%04d", seq_len(ncol(engine_x)))
  names(engine_x) <- engine_features
  formula_data <- sb_additional_survival_formula_data(
    engine_x,
    inputs$time,
    inputs$event
  )
  native_arguments <- list(
    formula = formula_data$formula,
    data = formula_data$data,
    dist = params$distribution,
    hessian = TRUE,
    control = list(maxit = numerical$maxit)
  )
  if (identical(params$distribution, "gompertz")) {
    # Use native normalized-time starting values and optimizer units, while
    # keeping the actual fit and all returned parameters in the input units.
    # Otherwise a long censored time can round its initial CDF to exactly one.
    time_scale <- mean(inputs$time)
    native_arguments$inits <- c(0.001 / time_scale, 1 / time_scale,
                               rep(0, ncol(engine_x)))
    native_arguments$control <- list(
      parscale = c(1 / time_scale, rep(1, ncol(engine_x) + 1L)),
      reltol = 1e-10,
      maxit = numerical$maxit
    )
  }
  critical_warnings <- character()
  engine <- tryCatch(
    withCallingHandlers(
      do.call(flexsurv::flexsurvreg, native_arguments),
      warning = function(warning) {
        message <- conditionMessage(warning)
        critical <- grepl(
          paste0(
            "Hessian not positive definite|Hessian is not finite|",
            "not converged to the maximum likelihood"
          ),
          message,
          ignore.case = TRUE
        )
        if (critical) {
          critical_warnings <<- c(critical_warnings, message)
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(error) {
      sb_abort(
        paste0(
          "Fixed flexsurv ", params$distribution, " fitting failed: ",
          conditionMessage(error)
        ),
        "stormbreaker_flexsurv_fit_failure",
        parent = error
      )
    }
  )
  if (length(critical_warnings) > 0L) {
    sb_abort(
      paste0(
        "Fixed flexsurv ", params$distribution,
        " reported an unreliable likelihood fit: ",
        paste(unique(critical_warnings), collapse = " | ")
      ),
      "stormbreaker_flexsurv_nonconvergence"
    )
  }
  convergence <- if (is.list(engine$opt)) engine$opt$convergence else NULL
  covariance_valid <- is.matrix(engine$cov) &&
    !anyNA(engine$cov) && all(is.finite(engine$cov))
  if (length(convergence) != 1L || is.na(convergence) || convergence != 0) {
    sb_abort(
      paste0(
        "Fixed flexsurv ", params$distribution,
        " did not converge (optimizer code ",
        if (length(convergence) == 1L && !is.na(convergence)) {
          convergence
        } else {
          "unavailable"
        },
        ")."
      ),
      "stormbreaker_flexsurv_nonconvergence"
    )
  }
  if (!covariance_valid) {
    sb_abort(
      paste0(
        "Fixed flexsurv ", params$distribution,
        " did not return a finite Hessian covariance matrix."
      ),
      "stormbreaker_flexsurv_nonconvergence"
    )
  }

  coefficients <- stats::coef(engine)
  ancillary_names <- sb_flexsurv_parametric_ancillary_names(
    params$distribution
  )
  ancillary <- unname(engine$res[ancillary_names, "est"])
  names(ancillary) <- ancillary_names
  valid <- is.numeric(coefficients) &&
    params$location_parameter %in% names(coefficients) &&
    all(engine_features %in% names(coefficients)) &&
    !anyNA(coefficients) && all(is.finite(coefficients)) &&
    is.numeric(ancillary) && length(ancillary) == length(ancillary_names) &&
    !anyNA(ancillary) && all(is.finite(ancillary))
  if (valid && params$distribution %in% c("gamma", "gengamma")) {
    valid <- ancillary[[1L]] > 0
  }
  if (!valid) {
    sb_abort(
      paste0(
        "Fixed flexsurv ", params$distribution,
        " did not return finite prediction parameters."
      ),
      "stormbreaker_incomplete_flexsurv_fit"
    )
  }
  prediction <- list(
    intercept = unname(coefficients[[params$location_parameter]]) -
      sum(unname(coefficients[engine_features]) / numerical$scale * numerical$center),
    coefficient = unname(coefficients[engine_features] / numerical$scale),
    ancillary = ancillary
  )
  model <- sb_new_flexsurv_parametric_model(
    learner,
    names(inputs$x),
    engine_features,
    params,
    prediction,
    engine
  )
  model$numerical_fitting <- numerical
  sb_predict_flexsurv_parametric(model, inputs$x)
  model
}

sb_validate_flexsurv_parametric_model_x <- function(model, x) {
  if (!inherits(model, "stormbreaker_flexsurv_parametric_fixed") ||
    !model$learner %in% sb_flexsurv_parametric_learner_ids()) {
    sb_abort(
      "The supplied model is not a fitted flexsurv parametric learner.",
      "stormbreaker_bad_flexsurv_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the flexsurv learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  sb_require_flexsurv_parametric_dependency(model$learner)
  x
}

sb_flexsurv_parametric_location <- function(model, x) {
  model$prediction$intercept +
    as.vector(as.matrix(x) %*% model$prediction$coefficient)
}

sb_predict_flexsurv_parametric <- function(model, x) {
  x <- sb_validate_flexsurv_parametric_model_x(model, x)
  location <- sb_flexsurv_parametric_location(model, x)
  risk <- if (identical(model$params$risk_transform, "log_location")) {
    location
  } else {
    -location
  }
  risk <- unname(as.numeric(risk))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid risk scores."),
      "stormbreaker_invalid_flexsurv_prediction"
    )
  }
  risk
}

sb_flexsurv_parametric_direct_probability <- function(model, x, horizon) {
  x <- sb_validate_flexsurv_parametric_model_x(model, x)
  location <- sb_flexsurv_parametric_location(model, x)
  ancillary <- model$prediction$ancillary
  probability <- switch(model$params$distribution,
    gamma = stats::pgamma(
      horizon,
      shape = ancillary[["shape"]],
      rate = exp(location),
      lower.tail = FALSE
    ),
    gompertz = flexsurv::pgompertz(
      horizon,
      shape = ancillary[["shape"]],
      rate = exp(location),
      lower.tail = FALSE
    ),
    gengamma = flexsurv::pgengamma(
      horizon,
      mu = location,
      sigma = ancillary[["sigma"]],
      Q = ancillary[["Q"]],
      lower.tail = FALSE
    )
  )
  probability <- unname(as.numeric(probability))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      paste0(
        "Learner `", model$learner,
        "` produced invalid horizon survival probabilities."
      ),
      "stormbreaker_invalid_flexsurv_prediction"
    )
  }
  probability
}

sb_flexsurv_parametric_direct_survival_source <- function(learner) {
  params <- sb_flexsurv_parametric_params(learner)
  paste0("flexsurv_", params$distribution, "_distribution_direct")
}

sb_predict_flexsurv_parametric_survival <- function(
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
      "Analysis columns do not match the flexsurv learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  list(
    prediction_horizon = prediction_horizon,
    source = sb_flexsurv_parametric_direct_survival_source(model$learner),
    ties = "not_applicable",
    probability = sb_flexsurv_parametric_direct_probability(
      model,
      x_assessment,
      prediction_horizon
    )
  )
}

sb_flexsurv_parametric_risk_reference <- function(learner) {
  params <- sb_flexsurv_parametric_params(learner)
  paste0(
    "flexsurv_", params$distribution, "_", params$parameterization,
    "_native"
  )
}

sb_flexsurv_parametric_risk_formula <- function(learner) {
  params <- sb_flexsurv_parametric_params(learner)
  if (identical(params$risk_transform, "log_location")) {
    "log_rate_location"
  } else {
    "negative_mu_log_time_location"
  }
}

sb_extract_flexsurv_parametric_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_flexsurv_parametric_fixed")) {
    sb_bad_learner_state("The fitted flexsurv learner state is invalid.")
  }
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_flexsurv_parametric_survival(
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
      ties = probability$ties
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      intercept = model$prediction$intercept,
      coefficient = model$prediction$coefficient,
      ancillary = model$prediction$ancillary,
      numerical_fitting = model$numerical_fitting,
      risk_reference = sb_flexsurv_parametric_risk_reference(model$learner),
      risk_formula = sb_flexsurv_parametric_risk_formula(model$learner),
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_flexsurv_parametric_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_flexsurv_parametric_model <- function(state, learner) {
  params <- sb_flexsurv_parametric_params(learner)
  ancillary_names <- sb_flexsurv_parametric_ancillary_names(
    params$distribution
  )
  valid <- sb_pipeline_parameters_equal(state$params, params) &&
    is.numeric(state$intercept) && length(state$intercept) == 1L &&
    !is.na(state$intercept) && is.finite(state$intercept) &&
    is.numeric(state$coefficient) &&
    length(state$coefficient) == length(state$features) &&
    !anyNA(state$coefficient) && all(is.finite(state$coefficient)) &&
    is.numeric(state$ancillary) &&
    identical(names(state$ancillary), ancillary_names) &&
    !anyNA(state$ancillary) && all(is.finite(state$ancillary))
  if (valid && params$distribution %in% c("gamma", "gengamma")) {
    valid <- state$ancillary[[1L]] > 0
  }
  if (!valid) {
    sb_bad_learner_state("The saved flexsurv prediction state is invalid.")
  }
  sb_new_flexsurv_parametric_model(
    learner,
    state$features,
    sprintf(".stormbreaker_x%04d", seq_along(state$features)),
    params,
    list(
      intercept = state$intercept,
      coefficient = state$coefficient,
      ancillary = state$ancillary
    )
  )
}

sb_predict_flexsurv_parametric_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_flexsurv_parametric_learner_ids()) {
    sb_bad_learner_state("The saved flexsurv learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  sb_require_flexsurv_parametric_dependency(learner)
  model <- sb_restore_flexsurv_parametric_model(state, learner)
  result <- list(risk = sb_predict_flexsurv_parametric(model, x))
  if (!is.null(state$survival_probability)) {
    probability <- state$survival_probability
    valid <- is.list(probability) &&
      is.numeric(probability$prediction_horizon) &&
      length(probability$prediction_horizon) == 1L &&
      !is.na(probability$prediction_horizon) &&
      is.finite(probability$prediction_horizon) &&
      probability$prediction_horizon > 0 &&
      is.character(probability$time_unit) &&
      length(probability$time_unit) == 1L &&
      !is.na(probability$time_unit) && nzchar(probability$time_unit) &&
      identical(
        probability$source,
        sb_flexsurv_parametric_direct_survival_source(learner)
      ) &&
      identical(probability$ties, "not_applicable")
    if (!valid) {
      sb_bad_learner_state(
        "The saved flexsurv survival-probability state is invalid."
      )
    }
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <-
      sb_flexsurv_parametric_direct_probability(
        model,
        x,
        probability$prediction_horizon
      )
  }
  result
}

sb_flexsurv_parametric_learner_adapter <- function(
  learner,
  params = sb_flexsurv_parametric_params(learner)
) {
  if (!learner %in% sb_flexsurv_parametric_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed flexsurv learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_flexsurv_parametric(x, time, event, learner, params)
    },
    predict = sb_predict_flexsurv_parametric,
    predict_survival = sb_predict_flexsurv_parametric_survival,
    extract_state = sb_extract_flexsurv_parametric_state,
    predict_state = sb_predict_flexsurv_parametric_state
  )
}

# ---- model-flexspline ----

sb_flexsurvspline_learner_ids <- function() {
  c(
    "flexsurvspline_hazard_k1",
    "flexsurvspline_odds_k1",
    "flexsurvspline_normal_k1"
  )
}

sb_flexsurvspline_params <- function(learner) {
  scale <- switch(learner,
    flexsurvspline_hazard_k1 = "hazard",
    flexsurvspline_odds_k1 = "odds",
    flexsurvspline_normal_k1 = "normal",
    sb_abort(
      paste0("Unknown fixed flexsurvspline learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
  list(
    scale = scale,
    k = 1L,
    timescale = "log",
    spline = "rp",
    knot_rule = "analysis_event_log_time_quantiles"
  )
}

sb_flexsurvspline_parameter_semantics <- function(learner) {
  params <- sb_flexsurvspline_params(learner)
  scale <- switch(params$scale,
    hazard = "log cumulative-hazard",
    odds = "log cumulative-odds",
    normal = "normal survival-transform"
  )
  paste0(
    "Royston-Parmar ", scale,
    " scale with one internal analysis-event log-time knot; ",
    "spline coefficients and knots estimated in each analysis fold"
  )
}

sb_flexsurvspline_dependency <- function(learner) {
  if (!learner %in% sb_flexsurvspline_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed flexsurvspline learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    package = "flexsurv",
    packages = "flexsurv",
    available = requireNamespace("flexsurv", quietly = TRUE)
  )
}

sb_require_flexsurvspline_dependency <- function(learner) {
  dependency <- sb_flexsurvspline_dependency(learner)
  if (!dependency$available) {
    sb_abort(
      paste0("Package `flexsurv` is required for learner `", learner, "`."),
      "stormbreaker_missing_flexsurvspline_dependency"
    )
  }
  invisible(dependency)
}

sb_flexsurvspline_gamma_names <- function(params) {
  paste0("gamma", 0:(params$k + 1L))
}

sb_new_flexsurvspline_model <- function(
  learner, features, engine_features, params, prediction, engine = NULL
) {
  structure(
    list(
      learner = learner,
      features = features,
      engine_features = engine_features,
      params = params,
      prediction = prediction,
      engine = engine
    ),
    class = "stormbreaker_flexsurvspline_fixed"
  )
}

sb_fit_flexsurvspline <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_flexsurvspline_dependency(learner)

  numerical <- sb_flexsurv_numerical_settings(inputs$x)
  engine_x <- sb_flexsurv_numerical_predictors(inputs$x, numerical)
  engine_features <- sprintf(".stormbreaker_x%04d", seq_len(ncol(engine_x)))
  names(engine_x) <- engine_features
  formula_data <- sb_additional_survival_formula_data(
    engine_x,
    inputs$time,
    inputs$event
  )
  critical_warnings <- character()
  engine <- tryCatch(
    withCallingHandlers(
      sb_flexsurvspline_numerical_fit(
        formula = formula_data$formula,
        data = formula_data$data,
        k = params$k,
        scale = params$scale,
        timescale = params$timescale,
        spline = params$spline,
        hessian = TRUE,
        maxit = numerical$maxit
      ),
      warning = function(warning) {
        message <- conditionMessage(warning)
        critical <- grepl(
          paste0(
            "Hessian not positive definite|Hessian is not finite|",
            "not converged to the maximum likelihood"
          ),
          message,
          ignore.case = TRUE
        )
        if (critical) {
          critical_warnings <<- c(critical_warnings, message)
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(error) {
      sb_abort(
        paste0(
          "Fixed flexsurvspline ", params$scale, "-scale fitting failed: ",
          conditionMessage(error)
        ),
        "stormbreaker_flexsurvspline_fit_failure",
        parent = error
      )
    }
  )
  if (length(critical_warnings) > 0L) {
    sb_abort(
      paste0(
        "Fixed flexsurvspline ", params$scale,
        "-scale fit reported an unreliable likelihood fit: ",
        paste(unique(critical_warnings), collapse = " | ")
      ),
      "stormbreaker_flexsurvspline_nonconvergence"
    )
  }
  convergence <- if (is.list(engine$opt)) engine$opt$convergence else NULL
  covariance_valid <- is.matrix(engine$cov) &&
    !anyNA(engine$cov) && all(is.finite(engine$cov))
  if (length(convergence) != 1L || is.na(convergence) || convergence != 0) {
    sb_abort(
      paste0(
        "Fixed flexsurvspline ", params$scale,
        "-scale fit did not converge (optimizer code ",
        if (length(convergence) == 1L && !is.na(convergence)) {
          convergence
        } else {
          "unavailable"
        },
        ")."
      ),
      "stormbreaker_flexsurvspline_nonconvergence"
    )
  }
  if (!covariance_valid) {
    sb_abort(
      paste0(
        "Fixed flexsurvspline ", params$scale,
        "-scale fit did not return a finite Hessian covariance matrix."
      ),
      "stormbreaker_flexsurvspline_nonconvergence"
    )
  }

  coefficients <- stats::coef(engine)
  gamma_names <- sb_flexsurvspline_gamma_names(params)
  gamma <- unname(coefficients[gamma_names])
  names(gamma) <- gamma_names
  coefficient <- unname(coefficients[engine_features])
  knots <- unname(engine$knots)
  valid <- is.numeric(gamma) && length(gamma) == length(gamma_names) &&
    !anyNA(gamma) && all(is.finite(gamma)) &&
    is.numeric(coefficient) && length(coefficient) == length(engine_features) &&
    !anyNA(coefficient) && all(is.finite(coefficient)) &&
    is.numeric(knots) && length(knots) == params$k + 2L &&
    !anyNA(knots) && all(is.finite(knots)) &&
    length(engine$k) == 1L && isTRUE(engine$k == params$k) &&
    identical(engine$scale, params$scale) &&
    identical(engine$timescale, params$timescale) &&
    identical(engine$spline, params$spline)
  if (!valid) {
    sb_abort(
      paste0(
        "Fixed flexsurvspline ", params$scale,
        "-scale fit did not return finite prediction parameters."
      ),
      "stormbreaker_incomplete_flexsurvspline_fit"
    )
  }
  gamma[[1L]] <- gamma[[1L]] - sum(coefficient / numerical$scale * numerical$center)
  coefficient <- unname(coefficient / numerical$scale)
  numerical$initialization <- "training_intercept_only_spline"
  prediction <- list(
    gamma = gamma,
    coefficient = coefficient,
    knots = knots
  )
  model <- sb_new_flexsurvspline_model(
    learner,
    names(inputs$x),
    engine_features,
    params,
    prediction,
    engine
  )
  model$numerical_fitting <- numerical
  sb_predict_flexsurvspline(model, inputs$x)
  model
}

sb_validate_flexsurvspline_model_x <- function(model, x) {
  if (!inherits(model, "stormbreaker_flexsurvspline_fixed") ||
    !model$learner %in% sb_flexsurvspline_learner_ids()) {
    sb_abort(
      "The supplied model is not a fitted flexsurvspline learner.",
      "stormbreaker_bad_flexsurvspline_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the flexsurvspline learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  sb_require_flexsurvspline_dependency(model$learner)
  x
}

sb_flexsurvspline_location <- function(model, x) {
  model$prediction$gamma[[1L]] +
    as.vector(as.matrix(x) %*% model$prediction$coefficient)
}

sb_predict_flexsurvspline <- function(model, x) {
  x <- sb_validate_flexsurvspline_model_x(model, x)
  risk <- unname(as.numeric(sb_flexsurvspline_location(model, x)))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid risk scores."),
      "stormbreaker_invalid_flexsurvspline_prediction"
    )
  }
  risk
}

sb_flexsurvspline_direct_probability <- function(model, x, horizon) {
  x <- sb_validate_flexsurvspline_model_x(model, x)
  gamma <- model$prediction$gamma
  gamma_matrix <- matrix(
    rep(gamma, each = nrow(x)),
    nrow = nrow(x),
    ncol = length(gamma)
  )
  gamma_matrix[, 1L] <- sb_flexsurvspline_location(model, x)
  probability <- flexsurv::psurvspline(
    q = rep(horizon, nrow(x)),
    gamma = gamma_matrix,
    knots = model$prediction$knots,
    scale = model$params$scale,
    timescale = model$params$timescale,
    spline = model$params$spline,
    lower.tail = FALSE
  )
  probability <- unname(as.numeric(probability))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      paste0(
        "Learner `", model$learner,
        "` produced invalid horizon survival probabilities."
      ),
      "stormbreaker_invalid_flexsurvspline_prediction"
    )
  }
  probability
}

sb_flexsurvspline_direct_survival_source <- function(learner) {
  params <- sb_flexsurvspline_params(learner)
  paste0(
    "flexsurvspline_", params$scale, "_k", params$k, "_",
    params$spline, "_distribution_direct"
  )
}

sb_predict_flexsurvspline_survival <- function(
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
      "Analysis columns do not match the flexsurvspline learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  list(
    prediction_horizon = prediction_horizon,
    source = sb_flexsurvspline_direct_survival_source(model$learner),
    ties = "not_applicable",
    probability = sb_flexsurvspline_direct_probability(
      model,
      x_assessment,
      prediction_horizon
    )
  )
}

sb_flexsurvspline_risk_reference <- function(learner) {
  params <- sb_flexsurvspline_params(learner)
  paste0(
    "flexsurvspline_", params$scale, "_k", params$k, "_",
    params$spline, "_native_location"
  )
}

sb_flexsurvspline_risk_formula <- function(learner) {
  sb_flexsurvspline_params(learner)
  "gamma0_plus_x_beta_native_location"
}

sb_extract_flexsurvspline_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_flexsurvspline_fixed")) {
    sb_bad_learner_state("The fitted flexsurvspline learner state is invalid.")
  }
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_flexsurvspline_survival(
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
      ties = probability$ties
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      gamma = model$prediction$gamma,
      coefficient = model$prediction$coefficient,
      knots = model$prediction$knots,
      numerical_fitting = model$numerical_fitting,
      risk_reference = sb_flexsurvspline_risk_reference(model$learner),
      risk_formula = sb_flexsurvspline_risk_formula(model$learner),
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_flexsurvspline_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_flexsurvspline_model <- function(state, learner) {
  params <- sb_flexsurvspline_params(learner)
  gamma_names <- sb_flexsurvspline_gamma_names(params)
  valid <- sb_pipeline_parameters_equal(state$params, params) &&
    is.numeric(state$gamma) &&
    identical(names(state$gamma), gamma_names) &&
    !anyNA(state$gamma) && all(is.finite(state$gamma)) &&
    is.numeric(state$coefficient) &&
    length(state$coefficient) == length(state$features) &&
    !anyNA(state$coefficient) && all(is.finite(state$coefficient)) &&
    is.numeric(state$knots) && length(state$knots) == params$k + 2L &&
    !anyNA(state$knots) && all(is.finite(state$knots))
  if (!valid) {
    sb_bad_learner_state(
      "The saved flexsurvspline prediction state is invalid."
    )
  }
  sb_new_flexsurvspline_model(
    learner,
    state$features,
    sprintf(".stormbreaker_x%04d", seq_along(state$features)),
    params,
    list(
      gamma = state$gamma,
      coefficient = state$coefficient,
      knots = state$knots
    )
  )
}

sb_predict_flexsurvspline_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_flexsurvspline_learner_ids()) {
    sb_bad_learner_state("The saved flexsurvspline learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  sb_require_flexsurvspline_dependency(learner)
  model <- sb_restore_flexsurvspline_model(state, learner)
  result <- list(risk = sb_predict_flexsurvspline(model, x))
  if (!is.null(state$survival_probability)) {
    probability <- state$survival_probability
    valid <- is.list(probability) &&
      is.numeric(probability$prediction_horizon) &&
      length(probability$prediction_horizon) == 1L &&
      !is.na(probability$prediction_horizon) &&
      is.finite(probability$prediction_horizon) &&
      probability$prediction_horizon > 0 &&
      is.character(probability$time_unit) &&
      length(probability$time_unit) == 1L &&
      !is.na(probability$time_unit) && nzchar(probability$time_unit) &&
      identical(
        probability$source,
        sb_flexsurvspline_direct_survival_source(learner)
      ) &&
      identical(probability$ties, "not_applicable")
    if (!valid) {
      sb_bad_learner_state(
        "The saved flexsurvspline survival-probability state is invalid."
      )
    }
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- sb_flexsurvspline_direct_probability(
      model,
      x,
      probability$prediction_horizon
    )
  }
  result
}

sb_flexsurvspline_learner_adapter <- function(
  learner,
  params = sb_flexsurvspline_params(learner)
) {
  if (!learner %in% sb_flexsurvspline_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed flexsurvspline learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_flexsurvspline(x, time, event, learner, params)
    },
    predict = sb_predict_flexsurvspline,
    predict_survival = sb_predict_flexsurvspline_survival,
    extract_state = sb_extract_flexsurvspline_state,
    predict_state = sb_predict_flexsurvspline_state
  )
}

# ---- model-rule-aft ----

# ---- algorithm-aft ----

#' @rdname sb_flexsurv_diagnostics
#' @export
sb_aft_diagnostics <- function(fit, type = "time_ratios") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  state <- fit$learner_state
  if (!state$learner %in% sb_survreg_learner_ids()) sb_abort(
    "This diagnostic requires one of the four survreg AFT models.",
    "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("log_time_coefficients", "time_ratios")) sb_abort(
    "Choose log_time_coefficients or time_ratios.", "stormbreaker_unsupported_diagnostic")
  sb_restore_additional_survival_model(state, state$learner)
  tab <- data.frame(feature = state$features, coefficient = state$coefficient,
    time_ratio = exp(state$coefficient))
  structure(list(type = paste0("aft_", type), table = tab,
    intercept = state$intercept, distribution_scale = state$scale,
    distribution = state$params$distribution, learner = state$learner,
    pipeline_id = fit$pipeline$id,
    source = "original saved survreg log-time location coefficients",
    interpretation = paste("Point estimates per original input unit, with all fitted features retained.",
      "exp(coefficient) multiplies conditional time quantiles; it is not a hazard ratio or a causal effect.",
      "The intercept and distribution scale are retained separately, not ranked as feature effects.",
      "No uncertainty is fabricated; ratios outside floating-point range remain in the table and are omitted from the ratio plot.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_aft_diagnostic <- function(x) {
  feature <- value <- NULL
  tab <- x$table
  ratio <- x$type == "aft_time_ratios"
  tab$value <- if (ratio) tab$time_ratio else tab$coefficient
  tab$value[!is.finite(tab$value) | (ratio & tab$value <= 0)] <- NA_real_
  if (!any(is.finite(tab$value))) sb_abort("No finite values are available on this axis; use log_time_coefficients.",
    "stormbreaker_unavailable_diagnostic")
  tab$feature <- factor(tab$feature, levels = rev(tab$feature))
  p <- ggplot2::ggplot(tab, ggplot2::aes(value, feature)) +
    ggplot2::geom_vline(xintercept = if (ratio) 1 else 0, colour = "grey65", linewidth = .4) +
    ggplot2::geom_point(colour = "#2878A0", size = 2.5, na.rm = TRUE) +
    ggplot2::labs(title = if (ratio) "AFT time ratios" else "AFT log-time coefficients",
      x = if (ratio) "Time ratio per input unit" else "Log-time coefficient per input unit", y = NULL) +
    sb_viz_theme()
  if (ratio) p <- p + ggplot2::scale_x_log10()
  attr(p, "stormbreaker_suggested_size") <- c(width = 9, height = max(5, 2 + .3 * nrow(tab)))
  p
}

# ---- algorithm-pre-rules ----

#' @rdname sb_gam_effects
#' @export
sb_rule_ensemble <- function(fit, type = "rule_coefficients") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!fit$pipeline$learner %in% sb_pre_rule_ensemble_learner_ids()) {
    sb_abort("This diagnostic requires an original PRE rule ensemble.",
      "stormbreaker_unsupported_diagnostic")
  }
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("rule_coefficients", "rule_cv_curve")) {
    sb_abort("Choose rule_coefficients or rule_cv_curve.",
      "stormbreaker_unsupported_diagnostic")
  }
  sb_require_rule_aft_dependency(fit$pipeline$learner)
  state <- fit$learner_state
  engine <- unserialize(state$engine_raw)
  if (!inherits(engine, "pre") || !inherits(engine$glmnet.fit, "cv.glmnet")) {
    sb_abort("The saved PRE engine is incomplete.", "stormbreaker_invalid_diagnostic")
  }
  rule <- state$params$penalty_rule
  if (!identical(rule, "lambda.1se")) {
    sb_abort("The saved PRE selection rule is unsupported.", "stormbreaker_invalid_diagnostic")
  }
  cv <- engine$glmnet.fit
  selected <- cv[[rule]]
  if (!is.numeric(selected) || length(selected) != 1L || !is.finite(selected) ||
      selected <= 0 || !selected %in% cv$lambda) {
    sb_abort("The original selected PRE lambda is invalid.", "stormbreaker_invalid_diagnostic")
  }
  if (type == "rule_coefficients") {
    tab <- stats::coef(engine, penalty.par.val = rule)
    if (!is.data.frame(tab) || !all(c("rule", "coefficient", "description") %in% names(tab)) ||
        anyNA(tab) || !is.numeric(tab$coefficient) || any(!is.finite(tab$coefficient)) ||
        anyDuplicated(tab$rule)) {
      sb_abort("Original PRE term definitions or coefficients are invalid.",
        "stormbreaker_invalid_diagnostic")
    }
    tab$term_type <- ifelse(tab$rule == "(Intercept)", "intercept",
      ifelse(tab$rule %in% engine$rules$rule, "rule", "linear"))
    bound_index <- match(tab$rule, engine$wins_points$varname)
    tab$winsor_lower <- engine$wins_points$lb[bound_index]
    tab$winsor_upper <- engine$wins_points$ub[bound_index]
    tab$selected <- tab$coefficient != 0
    rownames(tab) <- NULL
  } else {
    values <- list(lambda = cv$lambda, mean_value = cv$cvm,
      standard_error = cv$cvsd, nonzero = cv$nzero)
    if (!length(cv$lambda) || !all(vapply(values, function(z) {
      is.numeric(z) && length(z) == length(cv$lambda) && all(is.finite(z))
    }, logical(1))) || any(cv$lambda <= 0) || anyDuplicated(cv$lambda) ||
        any(cv$cvsd < 0) || any(cv$nzero < 0) || any(cv$nzero != floor(cv$nzero))) {
      sb_abort("The saved PRE coefficient-validation path is invalid.",
        "stormbreaker_invalid_diagnostic")
    }
    tab <- as.data.frame(values)
    tab$selected <- tab$lambda == selected
  }
  structure(list(type = type, table = tab, learner = fit$pipeline$learner,
    pipeline_id = fit$pipeline$id, selected_lambda = selected, selection_rule = rule,
    measure = cv$name, scale = if (engine$family == "binomial") "Log odds" else "Log hazard",
    source = "original saved PRE engine; no refitting",
    validation_scope = "coefficient-only CV after rule induction on complete development data"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_pre_diagnostic <- function(x) {
  if (x$type == "rule_cv_curve") {
    lambda <- mean_value <- standard_error <- NULL
    return(ggplot2::ggplot(x$table, ggplot2::aes(lambda, mean_value)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_value - standard_error,
        ymax = mean_value + standard_error), width = 0, colour = "#A6A6A6") +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::geom_point(size = 1, colour = "#2878A0") +
      ggplot2::geom_vline(xintercept = x$selected_lambda, linetype = 2, colour = "#D77935") +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = "PRE coefficient-only validation",
        subtitle = paste("Original selection:", x$selection_rule),
        x = "Lambda (log scale)", y = x$measure,
        caption = "Rules learned before internal coefficient CV; not validation of complete rule learning.\nBars: mean +/- one standard error, not a confidence interval. Dashed line: selected lambda.") +
      sb_viz_theme())
  }
  tab <- x$table[x$table$selected & x$table$term_type != "intercept", , drop = FALSE]
  tab$rule <- factor(tab$rule, levels = rev(tab$rule[order(abs(tab$coefficient), decreasing = TRUE)]))
  coefficient <- rule <- term_type <- NULL
  intercept <- x$table$coefficient[x$table$term_type == "intercept"]
  subtitle <- if (length(intercept)) paste("Original intercept:", format(signif(intercept, 4))) else
    "Cox ensemble: no intercept"
  p <- ggplot2::ggplot(tab, ggplot2::aes(coefficient, rule)) +
    ggplot2::labs(title = "PRE nonzero terms", subtitle = subtitle,
      x = paste(x$scale, "coefficient (native term units)"), y = "Native term ID",
      caption = "Complete term definitions and zero coefficients are retained in the diagnostic table.\nRule coefficients and linear slopes have different units; this is not feature importance.") +
    sb_viz_theme()
  if (nrow(tab)) {
    p <- p + ggplot2::geom_vline(xintercept = 0, colour = "#9A9A9A") +
      ggplot2::geom_col(ggplot2::aes(fill = term_type), width = .65) +
      ggplot2::scale_fill_manual(values = c(rule = "#2878A0", linear = "#D77935"), name = "Term")
  } else {
    p <- p + ggplot2::annotate("text", x = 0, y = 1,
      label = "No nonzero rule or linear terms at the original selected lambda.", size = 3.2) +
      ggplot2::scale_x_continuous(limits = c(-1, 1)) +
      ggplot2::scale_y_continuous(limits = c(0, 2)) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  }
  attr(p, "stormbreaker_suggested_size") <- c(width = 10, height = max(5, 2.5 + .28 * nrow(tab)))
  p
}

# ---- rule-aft-adapters ----

sb_rule_aft_learner_ids <- function() {
  c(
    "bart_aft_survival",
    "pre_binary_rule_ensemble",
    "pre_survival_rule_ensemble",
    "penaft_gehan_elastic_net"
  )
}

sb_pre_rule_ensemble_learner_ids <- function() {
  c("pre_binary_rule_ensemble", "pre_survival_rule_ensemble")
}

sb_rule_aft_survival_learner_ids <- function() {
  c(
    "bart_aft_survival",
    "pre_survival_rule_ensemble",
    "penaft_gehan_elastic_net"
  )
}

sb_rule_aft_params <- function(learner) {
  switch(learner,
    bart_aft_survival = list(
      type = "abart",
      error_distribution = "normal",
      num_trees = 50L,
      num_cutpoints = 100L,
      time_grid_size = 50L,
      posterior_draws = 250L,
      burn_in = 100L,
      keep_every = 1L,
      prior_k = 2,
      tree_power = 2,
      tree_base = 0.95,
      sparse = FALSE,
      seed = 1L,
      cores = 1L
    ),
    pre_binary_rule_ensemble = list(
      family = "binomial",
      ensemble_type = "both",
      sample_fraction = 0.5,
      maximum_depth = 3L,
      learning_rate = 0.01,
      mtry_rule = "ceiling_sqrt_p",
      num_trees = 500L,
      normalize = TRUE,
      standardize = FALSE,
      internal_folds = 5L,
      penalty_rule = "lambda.1se",
      seed = 1L
    ),
    pre_survival_rule_ensemble = list(
      family = "cox",
      ensemble_type = "both",
      sample_fraction = 0.5,
      maximum_depth = 3L,
      learning_rate = 0.01,
      mtry_rule = "ceiling_sqrt_p",
      num_trees = 500L,
      normalize = TRUE,
      standardize = FALSE,
      internal_folds = 5L,
      penalty_rule = "lambda.1se",
      seed = 1L
    ),
    penaft_gehan_elastic_net = list(
      penalty = "EN",
      alpha = 0.5,
      num_lambda = 50L,
      minimum_lambda_ratio = 0.1,
      standardize = TRUE,
      internal_folds = 5L,
      absolute_tolerance = 1e-8,
      relative_tolerance = 2.5e-4,
      maximum_iterations = 10000L,
      selection_rule = "minimum_cv_linear_predictor_score",
      seed = 1L
    ),
    sb_abort(
      paste0("Unknown rule/AFT learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_rule_aft_parameter_semantics <- function(learner) {
  params <- sb_rule_aft_params(learner)
  switch(learner,
    bart_aft_survival = paste0(
      "normal-error AFT BART with ", params$num_trees,
      " trees, ", params$posterior_draws, " posterior draws, ",
      params$burn_in, " burn-in draws, seed ", params$seed,
      ", and one-core execution fixed before resampling"
    ),
    pre_binary_rule_ensemble = paste0(
      "binary prediction-rule ensemble with linear terms, ",
      params$num_trees, " trees, depth ", params$maximum_depth,
      ", and ceiling-square-root mtry; the rule library is learned from the ",
      "complete analysis fold before seeded ", params$internal_folds,
      "-fold coefficient-only lambda.1se CV on that fixed library"
    ),
    pre_survival_rule_ensemble = paste0(
      "Cox prediction-rule ensemble with linear terms, ",
      params$num_trees, " trees, depth ", params$maximum_depth,
      ", and ceiling-square-root mtry; the rule library is learned from the ",
      "complete analysis fold before seeded ", params$internal_folds,
      "-fold coefficient-only lambda.1se CV on that fixed library"
    ),
    penaft_gehan_elastic_net = paste0(
      "semiparametric Gehan AFT elastic net with alpha ", params$alpha,
      ", ", params$num_lambda, " lambdas, and seeded event-stratified ",
      params$internal_folds,
      "-fold linear-predictor-score selection inside each analysis fold"
    )
  )
}

sb_rule_aft_dependency <- function(learner) {
  packages <- switch(learner,
    bart_aft_survival = "BART",
    pre_binary_rule_ensemble = "pre",
    pre_survival_rule_ensemble = c("pre", "mboost"),
    penaft_gehan_elastic_net = "penAFT",
    sb_abort(
      paste0("Unknown rule/AFT learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
  installed <- stats::setNames(
    vapply(packages, requireNamespace, logical(1), quietly = TRUE),
    packages
  )
  list(
    package = packages[[1L]],
    packages = packages,
    installed = installed,
    available = all(installed)
  )
}

sb_require_rule_aft_dependency <- function(learner) {
  dependency <- sb_rule_aft_dependency(learner)
  if (!isTRUE(dependency$available)) {
    missing <- dependency$packages[!dependency$installed]
    sb_abort(
      paste0(
        "Package", if (length(missing) > 1L) "s" else "", " `",
        paste(missing, collapse = "`, `"), "` ",
        if (length(missing) > 1L) "are" else "is",
        " required for learner `", learner, "`."
      ),
      "stormbreaker_missing_rule_aft_dependency"
    )
  }
  invisible(dependency)
}

sb_rule_aft_stop_on_warning <- function(
  code, learner, stage, allow_partial_glmnet_path = FALSE
) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      message <- conditionMessage(warning)
      partial_glmnet_path <- isTRUE(allow_partial_glmnet_path) &&
        grepl(
          paste0(
            "^from glmnet C\\+\\+ code \\(error code -[0-9]+\\); ",
            ".*solutions for larger (values of lambda|lambdas) returned$"
          ),
          message
        )
      if (partial_glmnet_path) {
        invokeRestart("muffleWarning")
      }
      sb_abort(
        paste0(
          "Learner `", learner, "` raised a ", stage, " warning: ",
          message
        ),
        "stormbreaker_rule_aft_warning",
        parent = warning
      )
    }
  )
}

sb_new_rule_aft_model <- function(learner, features, params, ...) {
  structure(
    c(
      list(learner = learner, features = features, params = params),
      list(...)
    ),
    class = "stormbreaker_rule_aft_model"
  )
}

sb_validate_rule_aft_model_x <- function(model, x) {
  valid <- inherits(model, "stormbreaker_rule_aft_model") &&
    is.list(model) && model$learner %in% sb_rule_aft_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features)
  if (!isTRUE(valid)) {
    sb_abort(
      "The supplied rule/AFT model is incomplete or invalid.",
      "stormbreaker_bad_rule_aft_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted rule/AFT learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  x
}

sb_compact_pre_engine <- function(engine) {
  engine[c("call", "weights", "data", "modmat")] <- NULL
  environment(engine$formula) <- baseenv()
  engine
}

sb_fit_pre_rule_ensemble <- function(x, endpoint, learner, params) {
  if (identical(learner, "pre_binary_rule_ensemble")) {
    inputs <- sb_validate_classical_binary_inputs(x, endpoint, learner)
    formula_data <- sb_classical_binary_formula_data(
      inputs$x,
      inputs$outcome
    )
  } else {
    inputs <- sb_validate_additional_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      learner
    )
    formula_data <- sb_additional_survival_formula_data(
      inputs$x,
      inputs$time,
      inputs$event
    )
  }
  sb_require_rule_aft_dependency(learner)
  effective_mtry <- as.integer(ceiling(sqrt(ncol(inputs$x))))
  engine <- NULL
  invisible(utils::capture.output(
    engine <- withr::with_seed(
      params$seed,
      sb_rule_aft_stop_on_warning(
        pre::pre(
          formula = formula_data$formula,
          data = formula_data$data,
          family = params$family,
          type = params$ensemble_type,
          sampfrac = params$sample_fraction,
          maxdepth = params$maximum_depth,
          learnrate = params$learning_rate,
          mtry = effective_mtry,
          ntrees = params$num_trees,
          normalize = params$normalize,
          standardize = params$standardize,
          nfolds = params$internal_folds,
          tree.unbiased = TRUE,
          removecomplements = TRUE,
          removeduplicates = TRUE,
          verbose = FALSE
        ),
        learner,
        "fitting",
        allow_partial_glmnet_path = TRUE
      )
    ),
    type = "output"
  ))
  if (!inherits(engine, "pre")) {
    sb_abort(
      "The PRE fit did not return a prediction-rule ensemble.",
      "stormbreaker_incomplete_pre_fit"
    )
  }
  sb_new_rule_aft_model(
    learner,
    names(inputs$x),
    params,
    engine = sb_compact_pre_engine(engine),
    effective_mtry = effective_mtry
  )
}

sb_fit_bart_aft <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_rule_aft_dependency(learner)
  engine <- NULL
  invisible(utils::capture.output(
    engine <- withr::with_seed(
      params$seed,
      sb_rule_aft_stop_on_warning(
        BART::abart(
          x.train = as.matrix(inputs$x),
          times = inputs$time,
          delta = inputs$event,
          K = params$time_grid_size,
          type = params$type,
          ntype = 1L,
          sparse = params$sparse,
          k = params$prior_k,
          power = params$tree_power,
          base = params$tree_base,
          ntree = params$num_trees,
          numcut = params$num_cutpoints,
          ndpost = params$posterior_draws,
          nskip = params$burn_in,
          keepevery = params$keep_every,
          printevery = .Machine$integer.max,
          mc.cores = params$cores,
          seed = params$seed
        ),
        learner,
        "fitting"
      )
    ),
    type = "output"
  ))
  post_sigma <- utils::tail(
    as.numeric(engine$sigma),
    params$posterior_draws
  )
  kept_columns <- as.integer(engine$rm.const)
  if (length(kept_columns) == 0L) {
    kept_columns <- seq_len(ncol(inputs$x))
  }
  valid <- inherits(engine, "abart") && is.list(engine$treedraws) &&
    is.numeric(engine$offset) && length(engine$offset) == 1L &&
    is.finite(engine$offset) &&
    length(post_sigma) == params$posterior_draws &&
    !anyNA(post_sigma) && all(is.finite(post_sigma)) &&
    all(post_sigma > 0) && length(kept_columns) > 0L &&
    !anyNA(kept_columns) && !anyDuplicated(kept_columns) &&
    all(kept_columns >= 1L & kept_columns <= ncol(inputs$x))
  if (!isTRUE(valid)) {
    sb_abort(
      "AFT-BART did not return a portable posterior prediction state.",
      "stormbreaker_incomplete_bart_aft_fit"
    )
  }
  compact_engine <- structure(
    list(
      treedraws = engine$treedraws,
      offset = unname(as.numeric(engine$offset)),
      sigma = unname(post_sigma),
      kept_columns = kept_columns
    ),
    class = "stormbreaker_bart_aft_engine"
  )
  sb_new_rule_aft_model(
    learner,
    names(inputs$x),
    params,
    engine = compact_engine
  )
}

sb_fit_penaft <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_rule_aft_dependency(learner)
  fold_id <- sb_make_stratified_fold_id(
    inputs$event,
    params$internal_folds,
    params$seed,
    stratum_label = "event"
  )
  cv_index <- lapply(seq_len(params$internal_folds), function(fold) {
    which(fold_id == fold)
  })
  fit <- NULL
  invisible(utils::capture.output(
    fit <- withr::with_seed(
      params$seed,
      sb_rule_aft_stop_on_warning(
        penAFT::penAFT.cv(
          X = as.matrix(inputs$x),
          logY = log(inputs$time),
          delta = inputs$event,
          nlambda = params$num_lambda,
          lambda.ratio.min = params$minimum_lambda_ratio,
          penalty = params$penalty,
          alpha = params$alpha,
          tol.abs = params$absolute_tolerance,
          tol.rel = params$relative_tolerance,
          standardize = params$standardize,
          cv.index = cv_index,
          admm.max.iter = params$maximum_iterations,
          quiet = TRUE
        ),
        learner,
        "fitting"
      )
    ),
    type = "output"
  ))
  cv_error <- as.numeric(fit$cv.err.linPred)
  selected_index <- if (length(cv_error) > 0L &&
    !anyNA(cv_error) && all(is.finite(cv_error))) {
    min(which(cv_error == min(cv_error)))
  } else {
    integer()
  }
  full <- fit$full.fit
  valid <- inherits(fit, "penAFT.cv") && length(selected_index) == 1L &&
    is.numeric(full$lambda) && length(full$lambda) == length(cv_error) &&
    !anyNA(full$lambda) && all(is.finite(full$lambda)) &&
    nrow(full$beta) == ncol(inputs$x) &&
    ncol(full$beta) == length(full$lambda)
  if (!isTRUE(valid)) {
    sb_abort(
      "penAFT did not return a complete cross-validated coefficient path.",
      "stormbreaker_incomplete_penaft_fit"
    )
  }
  coefficient <- as.numeric(as.matrix(full$beta[, selected_index]))
  if (isTRUE(full$standardize)) {
    coefficient <- coefficient / as.numeric(full$X.sd)
    intercept <- -sum(as.numeric(full$X.mean) * coefficient)
  } else {
    intercept <- 0
  }
  if (length(coefficient) != ncol(inputs$x) || anyNA(coefficient) ||
    any(!is.finite(coefficient)) || !is.finite(intercept)) {
    sb_abort(
      "penAFT did not return portable selected coefficients.",
      "stormbreaker_incomplete_penaft_fit"
    )
  }
  sb_new_rule_aft_model(
    learner,
    names(inputs$x),
    params,
    intercept = unname(intercept),
    coefficient = unname(coefficient),
    selected_lambda = unname(as.numeric(full$lambda[[selected_index]])),
    selected_cv_error = unname(cv_error[[selected_index]]),
    diagnostic_path = list(full_fit = full, cv_error = cv_error,
      selected_index = selected_index)
  )
}

sb_predict_pre_rule_ensemble <- function(model, x) {
  x <- sb_validate_rule_aft_model_x(model, x)
  sb_require_rule_aft_dependency(model$learner)
  prediction <- NULL
  prediction_type <- if (identical(
    model$learner,
    "pre_binary_rule_ensemble"
  )) {
    "response"
  } else {
    "link"
  }
  invisible(utils::capture.output(
    prediction <- sb_rule_aft_stop_on_warning(
      stats::predict(
        model$engine,
        newdata = x,
        type = prediction_type,
        penalty.par.val = model$params$penalty_rule
      ),
      model$learner,
      "prediction"
    ),
    type = "output"
  ))
  prediction <- unname(as.numeric(prediction))
  valid <- length(prediction) == nrow(x) && !anyNA(prediction) &&
    all(is.finite(prediction))
  if (identical(model$learner, "pre_binary_rule_ensemble")) {
    valid <- valid && all(prediction >= 0 & prediction <= 1)
  }
  if (!isTRUE(valid)) {
    sb_abort(
      "PRE produced incomplete or invalid predictions.",
      "stormbreaker_invalid_pre_prediction"
    )
  }
  prediction
}

sb_bart_aft_log_time_draws <- function(model, x) {
  x <- sb_validate_rule_aft_model_x(model, x)
  sb_require_rule_aft_dependency(model$learner)
  engine <- model$engine
  matrix_x <- as.matrix(x)[, engine$kept_columns, drop = FALSE]
  draws <- NULL
  invisible(utils::capture.output(
    draws <- sb_rule_aft_stop_on_warning(
      BART::pwbart(
        matrix_x,
        engine$treedraws,
        mu = engine$offset,
        mc.cores = model$params$cores,
        dodraws = TRUE
      ),
      model$learner,
      "prediction"
    ),
    type = "output"
  ))
  draws <- as.matrix(draws)
  valid <- nrow(draws) == length(engine$sigma) &&
    ncol(draws) == nrow(x) && !anyNA(draws) && all(is.finite(draws))
  if (!isTRUE(valid)) {
    sb_abort(
      "AFT-BART produced incomplete posterior log-time draws.",
      "stormbreaker_invalid_bart_aft_prediction"
    )
  }
  draws
}

sb_penaft_log_time_link <- function(model, x) {
  x <- sb_validate_rule_aft_model_x(model, x)
  link <- unname(
    model$intercept + as.vector(as.matrix(x) %*% model$coefficient)
  )
  if (length(link) != nrow(x) || anyNA(link) || any(!is.finite(link))) {
    sb_abort(
      "penAFT produced incomplete predicted log times.",
      "stormbreaker_invalid_penaft_prediction"
    )
  }
  link
}

sb_predict_rule_aft <- function(model, x) {
  learner <- model$learner
  if (learner %in% sb_pre_rule_ensemble_learner_ids()) {
    return(sb_predict_pre_rule_ensemble(model, x))
  }
  if (identical(learner, "bart_aft_survival")) {
    return(-colMeans(sb_bart_aft_log_time_draws(model, x)))
  }
  if (identical(learner, "penaft_gehan_elastic_net")) {
    return(-sb_penaft_log_time_link(model, x))
  }
  sb_abort(
    "The supplied rule/AFT model is not supported.",
    "stormbreaker_bad_rule_aft_fit"
  )
}

sb_bart_aft_survival_probability <- function(model, x, prediction_horizon) {
  draws <- sb_bart_aft_log_time_draws(model, x)
  standardized <- (log(prediction_horizon) - draws) / model$engine$sigma
  probability <- unname(colMeans(stats::pnorm(
    standardized,
    lower.tail = FALSE
  )))
  if (length(probability) != ncol(draws) || anyNA(probability) ||
    any(!is.finite(probability)) || any(probability < 0 | probability > 1)) {
    sb_abort(
      "AFT-BART produced invalid horizon survival probabilities.",
      "stormbreaker_invalid_bart_aft_prediction"
    )
  }
  probability
}

sb_rule_aft_direct_survival_source <- function(learner) {
  switch(learner,
    bart_aft_survival = "BART_abart_normal_posterior_survival_direct",
    penaft_gehan_elastic_net =
      "penAFT_gehan_analysis_residual_kaplan_meier",
    NULL
  )
}

sb_predict_rule_aft_survival <- function(
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
      "Analysis columns do not match the fitted rule/AFT learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  if (identical(model$learner, "bart_aft_survival")) {
    return(list(
      prediction_horizon = prediction_horizon,
      source = sb_rule_aft_direct_survival_source(model$learner),
      ties = "not_applicable",
      probability = sb_bart_aft_survival_probability(
        model,
        x_assessment,
        prediction_horizon
      )
    ))
  }
  if (identical(model$learner, "pre_survival_rule_ensemble")) {
    analysis_risk <- sb_predict_pre_rule_ensemble(model, inputs$x)
    assessment_risk <- sb_predict_pre_rule_ensemble(model, x_assessment)
    baseline <- sb_breslow_baseline_hazard(
      inputs$time,
      inputs$event,
      analysis_risk,
      prediction_horizon
    )
    return(sb_survival_prediction_result(
      prediction_horizon,
      baseline,
      source = "pre_rule_ensemble_analysis_breslow",
      ties = "breslow",
      risk = assessment_risk
    ))
  }
  analysis_link <- sb_penaft_log_time_link(model, inputs$x)
  assessment_link <- sb_penaft_log_time_link(model, x_assessment)
  curve <- sb_fit_mboost_gehan_residual_survival(
    inputs$time,
    inputs$event,
    analysis_link
  )
  list(
    prediction_horizon = prediction_horizon,
    source = sb_rule_aft_direct_survival_source(model$learner),
    ties = "not_applicable",
    probability = sb_predict_mboost_gehan_residual_survival(
      curve,
      assessment_link,
      prediction_horizon
    ),
    residual_survival = curve
  )
}

sb_extract_rule_aft_coefficients <- function(model) {
  if (!identical(model$learner, "penaft_gehan_elastic_net")) {
    sb_abort(
      "Only the penAFT learner exposes portable linear coefficients.",
      "stormbreaker_rule_aft_coefficients_not_available"
    )
  }
  values <- model$coefficient
  data.frame(
    feature = model$features,
    coefficient = values,
    nonzero = values != 0,
    direction = ifelse(
      values > 0,
      "longer_predicted_survival",
      ifelse(values < 0, "shorter_predicted_survival", "zero")
    ),
    stringsAsFactors = FALSE
  )
}

sb_extract_rule_aft_state <- function(model, context = list()) {
  learner <- model$learner
  is_pre <- learner %in% sb_pre_rule_ensemble_learner_ids()
  state <- list(
    learner = learner,
    features = model$features,
    params = model$params,
    baseline_hazard = NULL,
    survival_probability = NULL
  )
  if (is_pre || identical(learner, "bart_aft_survival")) {
    state$engine_raw <- serialize(
      model$engine,
      connection = NULL,
      ascii = FALSE,
      xdr = TRUE,
      version = 3L
    )
  }
  if (is_pre) {
    state$effective_mtry <- model$effective_mtry
  }
  if (identical(learner, "penaft_gehan_elastic_net")) {
    state$intercept <- model$intercept
    state$coefficient <- model$coefficient
    state$selected_lambda <- model$selected_lambda
    state$selected_cv_error <- model$selected_cv_error
  }
  if (identical(learner, "pre_binary_rule_ensemble")) {
    state$prediction_reference <- "pre_native_rule_ensemble"
    state$prediction_formula <- "positive_class_probability"
  } else {
    risk <- switch(learner,
      bart_aft_survival = list(
        reference = "BART_abart_normal_native",
        formula = "negative_posterior_mean_log_time"
      ),
      pre_survival_rule_ensemble = list(
        reference = "pre_native_rule_ensemble",
        formula = "pre_cox_link_prediction"
      ),
      penaft_gehan_elastic_net = list(
        reference = "penAFT_native",
        formula = "negative_centered_elastic_net_log_time_prediction"
      )
    )
    state$risk_reference <- risk$reference
    state$risk_formula <- risk$formula
  }
  if (learner %in% sb_rule_aft_survival_learner_ids() &&
    "prediction_horizon" %in% names(context)) {
    probability <- sb_predict_rule_aft_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    time_unit <- if ("time_unit" %in% names(context)) {
      context$time_unit
    } else {
      NA_character_
    }
    if (identical(learner, "pre_survival_rule_ensemble")) {
      state$baseline_hazard <- list(
        prediction_horizon = probability$prediction_horizon,
        time_unit = time_unit,
        cumulative_baseline_hazard = probability$cumulative_baseline_hazard,
        baseline_survival = probability$baseline_survival,
        source = probability$source,
        ties = probability$ties
      )
    } else {
      state$survival_probability <- list(
        prediction_horizon = probability$prediction_horizon,
        time_unit = time_unit,
        source = probability$source,
        ties = probability$ties
      )
      if (identical(learner, "penaft_gehan_elastic_net")) {
        state$survival_probability$residual_survival <-
          probability$residual_survival
      }
    }
  }
  structure(
    state,
    class = c("stormbreaker_rule_aft_state", "stormbreaker_learner_state")
  )
}

sb_restore_rule_aft_model <- function(state) {
  learner <- state$learner
  if (learner %in% sb_pre_rule_ensemble_learner_ids() ||
    identical(learner, "bart_aft_survival")) {
    if (!is.raw(state$engine_raw) || length(state$engine_raw) == 0L) {
      sb_bad_learner_state("The saved rule/AFT engine is invalid.")
    }
    engine <- tryCatch(
      unserialize(state$engine_raw),
      error = function(error) {
        sb_abort(
          paste0(
            "The saved `", learner, "` engine could not be loaded: ",
            conditionMessage(error)
          ),
          "stormbreaker_invalid_rule_aft_state",
          parent = error
        )
      }
    )
    valid_engine <- if (learner %in% sb_pre_rule_ensemble_learner_ids()) {
      inherits(engine, "pre")
    } else {
      inherits(engine, "stormbreaker_bart_aft_engine")
    }
    if (!isTRUE(valid_engine)) {
      sb_bad_learner_state("The saved rule/AFT engine is invalid.")
    }
    return(sb_new_rule_aft_model(
      learner,
      state$features,
      state$params,
      engine = engine,
      effective_mtry = state$effective_mtry
    ))
  }
  valid <- is.numeric(state$intercept) && length(state$intercept) == 1L &&
    is.finite(state$intercept) && is.numeric(state$coefficient) &&
    length(state$coefficient) == length(state$features) &&
    !anyNA(state$coefficient) && all(is.finite(state$coefficient))
  if (!isTRUE(valid)) {
    sb_bad_learner_state("The saved penAFT coefficients are invalid.")
  }
  sb_new_rule_aft_model(
    learner,
    state$features,
    state$params,
    intercept = state$intercept,
    coefficient = state$coefficient,
    selected_lambda = state$selected_lambda,
    selected_cv_error = state$selected_cv_error
  )
}

sb_validate_rule_aft_horizon_state <- function(state, learner) {
  probability <- if (identical(learner, "pre_survival_rule_ensemble")) {
    state$baseline_hazard
  } else {
    state$survival_probability
  }
  valid <- is.list(probability) &&
    is.numeric(probability$prediction_horizon) &&
    length(probability$prediction_horizon) == 1L &&
    is.finite(probability$prediction_horizon) &&
    probability$prediction_horizon > 0 &&
    is.character(probability$time_unit) &&
    length(probability$time_unit) == 1L &&
    !is.na(probability$time_unit) && nzchar(probability$time_unit)
  if (!isTRUE(valid)) {
    sb_bad_learner_state("The saved rule/AFT horizon state is invalid.")
  }
  probability
}

sb_predict_rule_aft_state <- function(state, x) {
  learner <- if (is.list(state)) state$learner else NULL
  if (!inherits(state, "stormbreaker_rule_aft_state") ||
    !is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_rule_aft_learner_ids()) {
    sb_bad_learner_state("The saved rule/AFT learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_restore_rule_aft_model(state)
  prediction <- sb_predict_rule_aft(model, x)
  if (identical(learner, "pre_binary_rule_ensemble")) {
    return(list(probability = prediction))
  }
  result <- list(risk = prediction)
  if (!is.null(state$baseline_hazard)) {
    baseline <- sb_validate_rule_aft_horizon_state(state, learner)
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      prediction,
      baseline$cumulative_baseline_hazard
    )
  } else if (!is.null(state$survival_probability)) {
    probability <- sb_validate_rule_aft_horizon_state(state, learner)
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    if (identical(learner, "bart_aft_survival")) {
      result$survival_probability <- sb_bart_aft_survival_probability(
        model,
        x,
        probability$prediction_horizon
      )
    } else {
      result$survival_probability <-
        sb_predict_mboost_gehan_residual_survival(
          probability$residual_survival,
          -prediction,
          probability$prediction_horizon
        )
    }
  }
  result
}

sb_rule_aft_learner_adapter <- function(
  learner,
  params = sb_rule_aft_params(learner)
) {
  if (!learner %in% sb_rule_aft_learner_ids()) {
    sb_abort(
      paste0("Unknown rule/AFT learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  fit <- if (identical(learner, "pre_binary_rule_ensemble")) {
    function(x, outcome) {
      sb_fit_pre_rule_ensemble(x, outcome, learner, params)
    }
  } else if (identical(learner, "pre_survival_rule_ensemble")) {
    function(x, time, event) {
      sb_fit_pre_rule_ensemble(
        x,
        list(time = time, event = event),
        learner,
        params
      )
    }
  } else if (identical(learner, "bart_aft_survival")) {
    function(x, time, event) {
      sb_fit_bart_aft(x, time, event, learner, params)
    }
  } else {
    function(x, time, event) {
      sb_fit_penaft(x, time, event, learner, params)
    }
  }
  adapter <- list(
    id = learner,
    params = params,
    fit = fit,
    predict = sb_predict_rule_aft,
    extract_state = sb_extract_rule_aft_state,
    predict_state = sb_predict_rule_aft_state
  )
  if (learner %in% sb_rule_aft_survival_learner_ids()) {
    adapter$predict_survival <- sb_predict_rule_aft_survival
  }
  if (identical(learner, "penaft_gehan_elastic_net")) {
    adapter$extract_coefficients <- sb_extract_rule_aft_coefficients
  }
  adapter
}
