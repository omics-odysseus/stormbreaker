# ---- specialized-binary-linear-adapters ----

sb_specialized_binary_linear_learner_ids <- function() {
  c(
    "stepglm_binary",
    "glmboost_binary",
    "plsrglm_binary"
  )
}

sb_validate_stepglm_direction <- function(direction) {
  if (!is.character(direction) || length(direction) != 1L ||
    is.na(direction) || !direction %in% c("forward", "backward", "both")) {
    sb_abort(
      "`stepglm_direction` must be one of `forward`, `backward`, or `both`.",
      "stormbreaker_bad_stepglm_direction"
    )
  }
  unname(as.character(direction))
}

sb_specialized_binary_linear_params <- function(
  learner,
  stepglm_direction = NULL
) {
  if (!learner %in% sb_specialized_binary_linear_learner_ids()) {
    sb_abort(
      paste0("Unknown specialized binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  if (identical(learner, "stepglm_binary")) {
    if (is.null(stepglm_direction)) {
      stepglm_direction <- "forward"
    }
    return(list(
      direction = sb_validate_stepglm_direction(stepglm_direction),
      criterion = "AIC",
      k = 2,
      max_steps = 1000L
    ))
  }
  if (!is.null(stepglm_direction)) {
    sb_abort(
      paste0(
        "`stepglm_direction` applies only when ",
        "`learner = 'stepglm_binary'`."
      ),
      "stormbreaker_irrelevant_stepglm_parameter"
    )
  }
  if (identical(learner, "glmboost_binary")) {
    return(list(
      iterations = 100L,
      shrinkage = 0.1,
      center = FALSE,
      family = "binomial_logit"
    ))
  }
  list(
    components = 2L,
    scale_x = TRUE,
    model = "pls-glm-logistic"
  )
}

sb_specialized_binary_linear_dependency <- function(learner) {
  package <- switch(learner,
    stepglm_binary = NA_character_,
    glmboost_binary = "mboost",
    plsrglm_binary = "plsRglm",
    sb_abort(
      paste0("Unknown specialized binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
  list(
    package = package,
    available = is.na(package) || requireNamespace(package, quietly = TRUE)
  )
}

sb_require_specialized_binary_linear_dependency <- function(learner) {
  dependency <- sb_specialized_binary_linear_dependency(learner)
  if (!dependency$available) {
    sb_abort(
      paste0(
        "Package `", dependency$package, "` is required for learner `",
        learner, "`."
      ),
      "stormbreaker_missing_specialized_binary_dependency"
    )
  }
  invisible(dependency)
}

sb_specialized_binary_stop_on_warning <- function(code, learner) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      sb_abort(
        paste0(
          "Learner `", learner, "` raised a fitting warning: ",
          conditionMessage(warning)
        ),
        paste0(
          "stormbreaker_",
          sub("_binary$", "", learner),
          "_warning"
        ),
        parent = warning
      )
    }
  )
}

sb_new_specialized_binary_linear_model <- function(
  learner, features, params, intercept, coefficient, center, effective
) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      intercept = intercept,
      coefficient = coefficient,
      center = center,
      effective = effective
    ),
    class = "stormbreaker_specialized_binary_linear"
  )
}

sb_fit_stepglm_binary <- function(inputs, params, capture_history = FALSE) {
  analysis <- inputs$x
  outcome_name <- ".stormbreaker_outcome"
  while (outcome_name %in% names(analysis)) {
    outcome_name <- paste0(outcome_name, "_")
  }
  analysis[[outcome_name]] <- inputs$outcome
  rhs <- paste(vapply(names(inputs$x), sb_quote_name, ""), collapse = " + ")
  full_formula <- stats::as.formula(
    paste(sb_quote_name(outcome_name), "~", rhs),
    env = environment()
  )
  null_formula <- stats::as.formula(
    paste(sb_quote_name(outcome_name), "~ 1"),
    env = environment()
  )
  start_formula <- if (identical(params$direction, "forward")) {
    null_formula
  } else {
    full_formula
  }
  model <- sb_specialized_binary_stop_on_warning(
    {
      start <- stats::glm(
        start_formula,
        data = analysis,
        family = stats::binomial(link = "logit"),
        model = FALSE,
        x = FALSE,
        y = FALSE
      )
      stats::step(
        start,
        scope = list(lower = null_formula, upper = full_formula),
        direction = params$direction,
        trace = 0,
        keep = NULL,
        steps = params$max_steps,
        k = params$k
      )
    },
    "stepglm_binary"
  )
  fitted <- stats::coef(model)
  selected <- setdiff(names(fitted), "(Intercept)")
  valid <- isTRUE(model$converged) &&
    identical(names(fitted)[[1L]], "(Intercept)") &&
    all(selected %in% names(inputs$x)) &&
    identical(as.integer(model$rank), as.integer(length(fitted))) &&
    !anyNA(fitted) && all(is.finite(fitted))
  if (!valid) {
    sb_abort(
      paste0(
        "Fixed Stepwise Logistic did not return a converged, full-rank fit ",
        "with finite coefficients."
      ),
      "stormbreaker_incomplete_stepglm_fit"
    )
  }
  coefficient <- stats::setNames(rep(0, ncol(inputs$x)), names(inputs$x))
  coefficient[selected] <- fitted[selected]
  list(
    intercept = unname(fitted[["(Intercept)"]]),
    coefficient = coefficient,
    center = stats::setNames(rep(0, ncol(inputs$x)), names(inputs$x)),
    effective = list(
      selected_features = selected,
      direction = params$direction,
      aic = unname(stats::AIC(model)),
      selection_path = if (capture_history) model$anova else NULL
    )
  )
}

sb_fit_glmboost_binary <- function(inputs, params, capture_history = FALSE) {
  model <- sb_specialized_binary_stop_on_warning(
    mboost::glmboost(
      x = as.matrix(inputs$x),
      y = factor(inputs$outcome, levels = c(0L, 1L)),
      family = mboost::Binomial(type = "glm", link = "logit"),
      control = mboost::boost_control(
        mstop = params$iterations,
        nu = params$shrinkage,
        risk = "inbag",
        stopintern = FALSE,
        trace = FALSE
      ),
      center = params$center
    ),
    "glmboost_binary"
  )
  fitted <- stats::coef(model)
  offset <- attr(fitted, "offset", exact = TRUE)
  selected <- names(fitted)
  center <- model$center
  valid <- is.numeric(fitted) && all(selected %in% names(inputs$x)) &&
    !anyNA(fitted) && all(is.finite(fitted)) &&
    is.numeric(offset) && length(offset) == 1L &&
    !is.na(offset) && is.finite(offset) &&
    is.numeric(center) && length(center) == ncol(inputs$x) &&
    !anyNA(center) && all(is.finite(center)) &&
    identical(as.integer(mboost::mstop(model)), params$iterations)
  if (!valid) {
    sb_abort(
      "Fixed glmBoost did not return finite linear prediction parameters.",
      "stormbreaker_incomplete_glmboost_fit"
    )
  }
  coefficient <- stats::setNames(rep(0, ncol(inputs$x)), names(inputs$x))
  coefficient[selected] <- fitted
  center <- stats::setNames(as.numeric(center), names(inputs$x))
  list(
    intercept = unname(offset),
    coefficient = coefficient,
    center = center,
    effective = list(
      selected_features = selected,
      iterations = params$iterations,
      training_loss = if (capture_history) data.frame(step = seq_along(mboost::risk(model)) - 1L,
        loss = as.numeric(mboost::risk(model))) else NULL
    )
  )
}

sb_fit_plsrglm_binary <- function(inputs, params) {
  if ("Intercept" %in% names(inputs$x)) {
    sb_abort(
      "Predictor name `Intercept` is reserved by the plsRglm engine.",
      "stormbreaker_reserved_plsrglm_intercept_name"
    )
  }
  model <- sb_specialized_binary_stop_on_warning(
    plsRglm::PLS_glm(
      dataY = inputs$outcome,
      dataX = inputs$x,
      nt = params$components,
      modele = params$model,
      typeVC = "none",
      scaleX = params$scale_x,
      pvals.expli = FALSE,
      MClassed = FALSE,
      sparse = FALSE,
      naive = FALSE,
      verbose = FALSE
    ),
    "plsrglm_binary"
  )
  fitted <- if (is.list(model)) drop(model$Coeffs) else NULL
  expected_terms <- c("Intercept", names(inputs$x))
  effective_components <- if (is.list(model)) model$computed_nt else NULL
  valid <- is.numeric(fitted) && identical(names(fitted), expected_terms) &&
    !anyNA(fitted) && all(is.finite(fitted)) &&
    is.numeric(effective_components) && length(effective_components) == 1L &&
    !is.na(effective_components) && is.finite(effective_components) &&
    effective_components >= 1L &&
    effective_components <= min(params$components, ncol(inputs$x)) &&
    is.list(model$FinalModel) && isTRUE(model$FinalModel$converged) &&
    !isTRUE(model$FinalModel$boundary)
  if (!valid) {
    sb_abort(
      "Fixed plsRglm did not return finite logistic prediction parameters.",
      "stormbreaker_incomplete_plsrglm_fit"
    )
  }
  list(
    intercept = unname(fitted[["Intercept"]]),
    coefficient = fitted[names(inputs$x)],
    center = stats::setNames(rep(0, ncol(inputs$x)), names(inputs$x)),
    effective = list(components = as.integer(effective_components)),
    diagnostic_components = list(loadings = model$pp, scores = model$tt,
      loading_definition = "Native plsRglm X reconstruction loadings")
  )
}

sb_fit_specialized_binary_linear <- function(
  x, outcome, learner, params
) {
  inputs <- sb_validate_classical_binary_inputs(x, outcome, learner)
  sb_validate_binary_predictor_names(names(inputs$x))
  sb_require_specialized_binary_linear_dependency(learner)
  fitted <- switch(learner,
    stepglm_binary = sb_fit_stepglm_binary(inputs, params),
    glmboost_binary = sb_fit_glmboost_binary(inputs, params),
    plsrglm_binary = sb_fit_plsrglm_binary(inputs, params)
  )
  model <- sb_new_specialized_binary_linear_model(
    learner = learner,
    features = names(inputs$x),
    params = params,
    intercept = fitted$intercept,
    coefficient = fitted$coefficient,
    center = fitted$center,
    effective = fitted$effective
  )
  if (!is.null(fitted$diagnostic_components)) model$diagnostic_components <- fitted$diagnostic_components
  model
}

sb_predict_specialized_binary_linear <- function(model, x) {
  if (!inherits(model, "stormbreaker_specialized_binary_linear")) {
    sb_abort(
      "The supplied model is not a fitted specialized binary learner.",
      "stormbreaker_bad_specialized_binary_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the specialized binary learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  matrix_x <- as.matrix(x)
  centered <- sweep(matrix_x, 2L, model$center, "-")
  probability <- stats::plogis(
    model$intercept + as.vector(centered %*% model$coefficient)
  )
  sb_class_one_probability(probability, nrow(x), model$learner)
}

sb_extract_specialized_binary_linear_coefficients <- function(model) {
  original_intercept <- model$intercept -
    sum(model$center * model$coefficient)
  data.frame(
    term = c("(Intercept)", model$features),
    coefficient = unname(c(original_intercept, model$coefficient)),
    stringsAsFactors = FALSE
  )
}

sb_extract_specialized_binary_linear_state <- function(model) {
  if (!inherits(model, "stormbreaker_specialized_binary_linear")) {
    sb_bad_learner_state(
      "The fitted specialized binary learner state is invalid."
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      intercept = model$intercept,
      coefficient = model$coefficient,
      center = model$center,
      effective = model$effective
    ),
    class = c(
      "stormbreaker_specialized_binary_linear_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_specialized_binary_linear_model <- function(state, learner) {
  valid <- is.list(state) && identical(state$learner, learner) &&
    is.character(state$features) && length(state$features) > 0L &&
    !anyNA(state$features) && !anyDuplicated(state$features) &&
    is.list(state$params) &&
    is.numeric(state$intercept) && length(state$intercept) == 1L &&
    !is.na(state$intercept) && is.finite(state$intercept) &&
    is.numeric(state$coefficient) &&
    identical(names(state$coefficient), state$features) &&
    !anyNA(state$coefficient) && all(is.finite(state$coefficient)) &&
    is.numeric(state$center) &&
    identical(names(state$center), state$features) &&
    !anyNA(state$center) && all(is.finite(state$center)) &&
    is.list(state$effective)
  if (!valid) {
    sb_bad_learner_state("The saved specialized binary learner state is invalid.")
  }
  sb_new_specialized_binary_linear_model(
    learner = learner,
    features = state$features,
    params = state$params,
    intercept = state$intercept,
    coefficient = state$coefficient,
    center = state$center,
    effective = state$effective
  )
}

sb_predict_specialized_binary_linear_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_specialized_binary_linear_learner_ids()) {
    sb_bad_learner_state("The saved specialized binary learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_restore_specialized_binary_linear_model(state, learner)
  list(probability = sb_predict_specialized_binary_linear(model, x))
}

sb_specialized_binary_linear_learner_adapter <- function(
  learner,
  params = sb_specialized_binary_linear_params(learner)
) {
  if (!learner %in% sb_specialized_binary_linear_learner_ids()) {
    sb_abort(
      paste0("Unknown specialized binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_specialized_binary_linear(x, outcome, learner, params)
    },
    predict = sb_predict_specialized_binary_linear,
    extract_coefficients = sb_extract_specialized_binary_linear_coefficients,
    extract_state = function(model, context = list()) {
      sb_extract_specialized_binary_linear_state(model)
    },
    predict_state = sb_predict_specialized_binary_linear_state
  )
}

# ---- model-logistic ----

# ---- extended-logistic-parameters ----

sb_extended_logistic_learner_ids <- function() c(sb_specialized_binary_linear_learner_ids(),
  "brglm2_logistic", "bayesglm_logistic_cauchy",
  "liblinear_l2_logistic_binary", "liblinear_l1_logistic_binary",
  "robust_logistic_mqle")

sb_extended_logistic_parameters <- function(state) {
  if (identical(state$learner, "brglm2_logistic") && sb_is_constant_sum_state(state)) {
    return(sb_constant_sum_linear_parameters(state))
  }
  if (state$learner %in% c("liblinear_l2_logistic_binary",
                          "liblinear_l1_logistic_binary")) {
    e <- unserialize(state$engine_raw)
    m <- e$model
    features <- state$features
    expected_type <- if (state$learner == "liblinear_l2_logistic_binary") 0L else 6L
    if (!inherits(m, "LiblineaR") || !identical(as.integer(m$Type), expected_type) ||
        !identical(sort(as.integer(m$ClassNames)), c(0L, 1L)) ||
        !is.matrix(m$W) || !identical(dim(m$W), c(1L, length(features) + 1L)) ||
        any(!is.finite(m$W)) || !identical(names(e$center), features) ||
        !identical(names(e$scale), features) || any(!is.finite(e$center)) ||
        any(!is.finite(e$scale)) || any(e$scale <= 0) ||
        length(m$Bias) != 1L || !is.finite(m$Bias) || m$Bias < 0) sb_abort(
      "Original LIBLINEAR logistic weights, class order or scaling are invalid.",
      "stormbreaker_bad_learner_state")
    orientation <- if (m$ClassNames[1L] == 1L) 1 else -1
    w <- as.numeric(m$W[1L, seq_along(features)])
    return(list(coefficient = stats::setNames(orientation * w / e$scale, features),
      intercept = unname(orientation * (m$W[1L, length(features) + 1L] * m$Bias -
        sum(w * e$center / e$scale)))))
  }
  if (state$learner %in% sb_specialized_binary_linear_learner_ids()) {
    model <- sb_restore_specialized_binary_linear_model(state, state$learner)
    return(list(coefficient = model$coefficient,
      intercept = model$intercept - sum(model$center * model$coefficient)))
  }
  if (!state$learner %in% c("brglm2_logistic", "bayesglm_logistic_cauchy",
                          "robust_logistic_mqle")) sb_abort(
    "This learner has no extended logistic parameter adapter.", "stormbreaker_unsupported_coefficients")
  p <- if (state$learner == "robust_logistic_mqle") unserialize(state$engine_raw) else state$engine
  if (!is.numeric(p$coefficient) || !identical(names(p$coefficient), state$features) ||
      any(!is.finite(p$coefficient)) || !is.numeric(p$intercept) || length(p$intercept) != 1L ||
      !is.finite(p$intercept)) sb_abort("Saved logistic parameters are invalid.", "stormbreaker_bad_learner_state")
  list(coefficient = p$coefficient, intercept = p$intercept)
}

sb_linear_coefficients_diagnostic <- function(fit) {
  parameters <- sb_model_coefficients(fit)
  structure(list(type = "linear_coefficients", table = parameters$table,
    parameters = parameters, learner = parameters$learner, pipeline_id = parameters$pipeline_id,
    source = "original saved linear predictor in input units",
    interpretation = paste("Signed coefficients on", parameters$score_scale,
      "; original intercept retained separately. Not SHAP, causal effects or confidence intervals.")),
    class = "stormbreaker_algorithm_diagnostic")
}

# ---- logistic-adapter ----

sb_validate_binary_predictor_names <- function(predictors) {
  if ("(Intercept)" %in% predictors) {
    sb_abort(
      paste0(
        "Predictor name `(Intercept)` is reserved for the model intercept ",
        "in binary pipelines."
      ),
      "stormbreaker_reserved_intercept_name"
    )
  }
  invisible(predictors)
}

sb_validate_binary_outcome_vector <- function(outcome, n) {
  valid <- (is.numeric(outcome) || is.logical(outcome)) &&
    length(outcome) == n && !anyNA(outcome) &&
    (!is.numeric(outcome) || all(is.finite(outcome))) &&
    all(outcome %in% c(0, 1))
  if (!valid) {
    sb_abort(
      "`outcome` must contain one finite 0/1 value per analysis patient.",
      "stormbreaker_bad_adapter_binary_outcome"
    )
  }
  as.integer(outcome)
}

sb_fit_logistic <- function(x, outcome) {
  x <- sb_validate_predictor_frame(x, "x")
  sb_validate_binary_predictor_names(names(x))
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "A logistic analysis fold must contain both outcome classes.",
      "stormbreaker_logistic_single_class"
    )
  }

  outcome_name <- ".stormbreaker_outcome"
  while (outcome_name %in% names(x)) {
    outcome_name <- paste0(outcome_name, "_")
  }
  analysis <- x
  analysis[[outcome_name]] <- outcome
  rhs <- paste(vapply(names(x), sb_quote_name, ""), collapse = " + ")
  formula <- stats::as.formula(
    paste(sb_quote_name(outcome_name), "~", rhs),
    env = parent.frame()
  )
  fit_warnings <- character()
  model <- withCallingHandlers(
    stats::glm(
      formula,
      data = analysis,
      family = stats::binomial(link = "logit"),
      x = TRUE,
      y = TRUE,
      model = TRUE
    ),
    warning = function(warning) {
      if (identical(
        conditionMessage(warning),
        "glm.fit: fitted probabilities numerically 0 or 1 occurred"
      )) {
        fit_warnings <<- c(fit_warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
      sb_abort(
        paste0(
          "The logistic learner raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_logistic_warning",
        parent = warning
      )
    }
  )
  coefficients <- stats::coef(model)
  fitted_matrix <- unname(model$x[, -1L, drop = FALSE])
  expected_matrix <- unname(as.matrix(x))
  storage.mode(expected_matrix) <- "double"
  valid_fit <- isTRUE(model$converged) && NROW(model$y) == nrow(analysis) &&
    length(coefficients) == ncol(x) + 1L &&
    model$rank == ncol(x) + 1L &&
    identical(fitted_matrix, expected_matrix) &&
    !anyNA(coefficients) && all(is.finite(coefficients))
  if (!valid_fit) {
    sb_abort(
      paste0(
        "The logistic learner did not produce a converged, full-rank fit ",
        "with one finite coefficient per declared term."
      ),
      "stormbreaker_invalid_logistic_fit"
    )
  }
  if (length(fit_warnings) > 0L) {
    # A saturated probability can coexist with a finite MLE. Reuse the
    # existing post-fit separation calculation rather than ignoring a warning
    # that can also accompany an unbounded fit.
    separation <- sb_binary_moderate_calibration_separation(model$x, outcome)
    if (identical(separation$status, "failed")) {
      sb_abort(
        paste0(
          "The logistic separation calculation failed: ",
          separation$error_message
        ),
        "stormbreaker_logistic_separation_check_failed"
      )
    }
    if (isTRUE(separation$detected)) {
      sb_abort(
        paste0(
          "The logistic design has complete or quasi-complete separation; ",
          "no finite maximum-likelihood fit is reported. Use fewer predictors ",
          "or explicitly choose a regularized learner."
        ),
        "stormbreaker_logistic_separation"
      )
    }
  }
  attr(model, "stormbreaker_features") <- names(x)
  attr(model, "stormbreaker_fit_warnings") <- fit_warnings
  model
}

sb_predict_logistic <- function(model, x) {
  x <- sb_validate_predictor_frame(x, "x")
  expected_features <- attr(model, "stormbreaker_features", exact = TRUE)
  if (!identical(names(x), expected_features)) {
    sb_abort(
      "Prediction columns do not match the fitted logistic learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  probability <- unname(stats::predict(model, newdata = x, type = "response"))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      "The logistic learner produced invalid class-one probabilities.",
      "stormbreaker_invalid_logistic_prediction"
    )
  }
  probability
}

sb_logistic_learner_adapter <- function() {
  list(
    id = "logistic",
    fit = sb_fit_logistic,
    predict = sb_predict_logistic,
    extract_coefficients = function(model) {
      coefficients <- stats::coef(model)
      features <- attr(model, "stormbreaker_features", exact = TRUE)
      data.frame(
        term = c("(Intercept)", features),
        coefficient = unname(coefficients),
        stringsAsFactors = FALSE
      )
    },
    extract_state = function(model, context = list()) {
      coefficients <- stats::coef(model)
      features <- attr(model, "stormbreaker_features", exact = TRUE)
      sb_extract_logistic_learner_state(
        learner = "logistic",
        coefficients = data.frame(
          term = c("(Intercept)", features),
          coefficient = unname(coefficients),
          stringsAsFactors = FALSE
        )
      )
    },
    predict_state = sb_predict_logistic_learner_state
  )
}

# ---- model-glmnet ----

# ---- glmnet-active-set-stability ----

sb_glmnet_active_set_stability <- function(pipeline, coefficients,
                                            task_predictors) {
  is_glmnet <- pipeline$learner %in% c("glmnet_cox", "glmnet_logistic")
  alpha <- if (is_glmnet) {
    as.numeric(pipeline$learner_params$alpha)
  } else {
    NA_real_
  }
  applicable <- is_glmnet && alpha > 0
  fixed_features <- identical(pipeline$selector, "fixed_features")
  candidate_features <- if (fixed_features) {
    pipeline$selector_params$features
  } else {
    task_predictors
  }
  active_sets <- lapply(coefficients, function(audit) {
    success <- identical(audit$status, "success")
    list(
      fold = as.integer(audit$fold),
      status = if (success) "success" else audit$status,
      selected_features = if (success) {
        as.character(audit$nonzero_features)
      } else {
        character()
      },
      error_class = NA_character_,
      error_message = NA_character_
    )
  })
  summary <- sb_selection_stability(
    if (applicable) "glmnet_active_set" else "none",
    active_sets,
    candidate_features
  )

  all_folds_fitted <- summary$n_folds >= 2L &&
    summary$n_successful_folds == summary$n_folds
  empty_union <- any(summary$pairwise$status == "not_evaluable")
  if (applicable && identical(summary$status, "incomplete") &&
      all_folds_fitted && empty_union) {
    summary$status <- "not_evaluable"
    summary$error_message <- paste0(
      "Glmnet active-set stability is not evaluable because at least one ",
      "outer-fold pair has an empty feature union."
    )
  }

  list(
    learner = pipeline$learner,
    alpha = alpha,
    status = summary$status,
    active_set_definition = if (applicable) {
      "exact_nonzero_predictor_coefficients_at_selected_lambda"
    } else {
      NA_character_
    },
    scope = "outer_analysis_fold_glmnet_active_sets",
    aggregation = summary$aggregation,
    frequency_denominator = "successful_glmnet_coefficient_folds",
    n_candidate_features = summary$n_candidate_features,
    n_folds = summary$n_folds,
    n_successful_folds = summary$n_successful_folds,
    n_pairs = summary$n_pairs,
    n_evaluable_pairs = summary$n_evaluable_pairs,
    mean_pairwise_jaccard = summary$mean_pairwise_jaccard,
    feature_frequency = summary$feature_frequency,
    per_fold = summary$per_fold,
    pairwise = summary$pairwise,
    error_message = summary$error_message
  )
}

# ---- glmnet-cox-adapter ----

sb_glmnet_available <- function() {
  requireNamespace("glmnet", quietly = TRUE)
}

sb_glmnet_lambda_in_path <- function(lambda, path) {
  limits <- range(path)
  # Screening can change column order and round an otherwise identical path
  # endpoint. Accept roundoff only; retain the requested lambda/native path.
  tolerance <- 16 * .Machine$double.eps
  (limits[[1L]] - lambda) / limits[[1L]] <= tolerance &&
    (lambda - limits[[2L]]) / limits[[2L]] <= tolerance
}

sb_predict_glmnet_cox <- function(model, x) {
  x <- sb_validate_predictor_frame(x, "x")
  expected_features <- attr(
    model,
    "stormbreaker_features",
    exact = TRUE
  )
  if (!identical(names(x), expected_features)) {
    sb_abort(
      paste0(
        "Prediction columns do not match the fitted glmnet Cox ",
        "learner schema."
      ),
      "stormbreaker_learner_schema_mismatch"
    )
  }
  lambda <- attr(model, "stormbreaker_lambda", exact = TRUE)
  risk <- drop(stats::predict(
    model,
    newx = as.matrix(x),
    s = lambda,
    type = "link"
  ))
  if (length(risk) != nrow(x) || anyNA(risk) ||
      any(!is.finite(risk))) {
    sb_abort(
      paste0(
        "The glmnet Cox learner produced incomplete or non-finite ",
        "risk predictions."
      ),
      "stormbreaker_invalid_glmnet_prediction"
    )
  }
  unname(risk)
}

sb_predict_glmnet_cox_survival <- function(model, x_analysis, time, event,
                                           x_assessment,
                                           prediction_horizon) {
  x_analysis <- sb_validate_predictor_frame(x_analysis, "x_analysis")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x_analysis))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x_analysis))
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    time
  )
  analysis_risk <- sb_predict_glmnet_cox(model, x_analysis)
  assessment_risk <- sb_predict_glmnet_cox(model, x_assessment)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    time,
    as.integer(event),
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = "glmnet_fixed_beta_breslow",
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_glmnet_cox_learner_adapter <- function(alpha, lambda,
                                          standardize = TRUE) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha < 0 || alpha > 1) {
    sb_abort(
      "`alpha` must be one finite number between 0 and 1.",
      "stormbreaker_bad_glmnet_alpha"
    )
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) ||
      !is.finite(lambda) || lambda <= 0) {
    sb_abort(
      "`lambda` must be one finite positive number.",
      "stormbreaker_bad_glmnet_lambda"
    )
  }
  if (!is.logical(standardize) || length(standardize) != 1L ||
      is.na(standardize)) {
    sb_abort(
      "`standardize` must be `TRUE` or `FALSE`.",
      "stormbreaker_bad_glmnet_standardize"
    )
  }

  params <- list(
    alpha = as.numeric(alpha),
    lambda = as.numeric(lambda),
    standardize = standardize
  )
  list(
    id = "glmnet_cox",
    params = params,
    fit = function(x, time, event) {
      if (!sb_glmnet_available()) {
        sb_abort(
          paste0(
            "Package `glmnet` is required for the glmnet Cox learner."
          ),
          "stormbreaker_missing_glmnet"
        )
      }
      x <- sb_validate_predictor_frame(x, "x")
      time <- sb_validate_endpoint_vector(time, "time", nrow(x))
      event <- sb_validate_endpoint_vector(event, "event", nrow(x))
      if (any(time <= 0)) {
        sb_abort(
          paste0(
            "The glmnet Cox learner requires strictly positive ",
            "survival times."
          ),
          "stormbreaker_glmnet_nonpositive_time"
        )
      }
      if (!all(event %in% c(0, 1))) {
        sb_abort(
          "`event` must contain only 0/1 values.",
          "stormbreaker_bad_adapter_event"
        )
      }

      path_warnings <- character()
      model <- withCallingHandlers(
        glmnet::glmnet(
          x = as.matrix(x),
          y = survival::Surv(time, event),
          family = "cox",
          cox.ties = "breslow",
          alpha = params$alpha,
          standardize = params$standardize
        ),
        warning = function(warning) {
          if (sb_is_glmnet_tail_convergence_warning(warning)) {
            path_warnings <<- c(path_warnings, conditionMessage(warning))
            invokeRestart("muffleWarning")
          }
          sb_abort(
            paste0(
              "The glmnet Cox learner raised a warning: ",
              conditionMessage(warning)
            ),
            "stormbreaker_glmnet_fit_warning",
            parent = warning
          )
        }
      )
      lambda_min <- min(model$lambda)
      lambda_max <- max(model$lambda)
      if (!sb_glmnet_lambda_in_path(params$lambda, model$lambda)) {
        sb_abort(
          paste0(
            "The requested `lambda` is outside the fitted glmnet Cox ",
            "path [", format(lambda_min), ", ", format(lambda_max), "]."
          ),
          "stormbreaker_glmnet_lambda_outside_path"
        )
      }
      attr(model, "stormbreaker_features") <- names(x)
      attr(model, "stormbreaker_lambda") <- params$lambda
      attr(model, "stormbreaker_fit_warnings") <- path_warnings
      model
    },
    predict = sb_predict_glmnet_cox,
    predict_survival = sb_predict_glmnet_cox_survival,
    extract_coefficients = function(model) {
      features <- attr(model, "stormbreaker_features", exact = TRUE)
      lambda <- attr(model, "stormbreaker_lambda", exact = TRUE)
      coefficient_matrix <- as.matrix(stats::coef(model, s = lambda))
      coefficient <- unname(coefficient_matrix[features, 1L])
      nonzero <- coefficient != 0
      data.frame(
        feature = features,
        coefficient = coefficient,
        nonzero = nonzero,
        direction = ifelse(
          coefficient > 0,
          "higher_hazard",
          ifelse(coefficient < 0, "lower_hazard", "zero")
        ),
        stringsAsFactors = FALSE
      )
    },
    extract_state = function(model, context = list()) {
      features <- attr(model, "stormbreaker_features", exact = TRUE)
      lambda <- attr(model, "stormbreaker_lambda", exact = TRUE)
      coefficient_matrix <- as.matrix(stats::coef(model, s = lambda))
      coefficient <- unname(coefficient_matrix[features, 1L])
      coefficients <- data.frame(
        feature = features,
        coefficient = coefficient,
        stringsAsFactors = FALSE
      )
      sb_extract_cox_learner_state(
        model,
        learner = "glmnet_cox",
        coefficients = coefficients,
        predict_survival = sb_predict_glmnet_cox_survival,
        context = context,
        params = params
      )
    },
    predict_state = sb_predict_cox_learner_state
  )
}

# ---- glmnet-logistic-adapter ----

sb_glmnet_logistic_learner_adapter <- function(alpha, lambda,
                                               standardize = TRUE) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
    !is.finite(alpha) || alpha < 0 || alpha > 1) {
    sb_abort(
      "`alpha` must be one finite number between 0 and 1.",
      "stormbreaker_bad_glmnet_alpha"
    )
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) ||
    !is.finite(lambda) || lambda <= 0) {
    sb_abort(
      "`lambda` must be one finite positive number.",
      "stormbreaker_bad_glmnet_lambda"
    )
  }
  if (!is.logical(standardize) || length(standardize) != 1L ||
    is.na(standardize)) {
    sb_abort(
      "`standardize` must be `TRUE` or `FALSE`.",
      "stormbreaker_bad_glmnet_standardize"
    )
  }

  params <- list(
    alpha = as.numeric(alpha),
    lambda = as.numeric(lambda),
    standardize = standardize
  )
  list(
    id = "glmnet_logistic",
    params = params,
    fit = function(x, outcome) {
      if (!sb_glmnet_available()) {
        sb_abort(
          "Package `glmnet` is required for the glmnet logistic learner.",
          "stormbreaker_missing_glmnet"
        )
      }
      x <- sb_validate_predictor_frame(x, "x")
      sb_validate_binary_predictor_names(names(x))
      outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
      if (!all(c(0L, 1L) %in% outcome)) {
        sb_abort(
          "A glmnet logistic analysis fold must contain both classes.",
          "stormbreaker_logistic_single_class"
        )
      }

      path_warnings <- character()
      model <- withCallingHandlers(
        glmnet::glmnet(
          x = as.matrix(x),
          y = outcome,
          family = "binomial",
          alpha = params$alpha,
          standardize = params$standardize
        ),
        warning = function(warning) {
          if (sb_is_glmnet_tail_convergence_warning(warning)) {
            path_warnings <<- c(path_warnings, conditionMessage(warning))
            invokeRestart("muffleWarning")
          }
          sb_abort(
            paste0(
              "The glmnet logistic learner raised a warning: ",
              conditionMessage(warning)
            ),
            "stormbreaker_glmnet_logistic_fit_warning",
            parent = warning
          )
        }
      )
      lambda_min <- min(model$lambda)
      lambda_max <- max(model$lambda)
      if (!sb_glmnet_lambda_in_path(params$lambda, model$lambda)) {
        sb_abort(
          paste0(
            "The requested `lambda` is outside the fitted glmnet logistic ",
            "path [", format(lambda_min), ", ", format(lambda_max), "]."
          ),
          "stormbreaker_glmnet_lambda_outside_path"
        )
      }
      attr(model, "stormbreaker_features") <- names(x)
      attr(model, "stormbreaker_lambda") <- params$lambda
      attr(model, "stormbreaker_fit_warnings") <- path_warnings
      model
    },
    predict = function(model, x) {
      x <- sb_validate_predictor_frame(x, "x")
      expected_features <- attr(
        model,
        "stormbreaker_features",
        exact = TRUE
      )
      if (!identical(names(x), expected_features)) {
        sb_abort(
          paste0(
            "Prediction columns do not match the fitted glmnet logistic ",
            "learner schema."
          ),
          "stormbreaker_learner_schema_mismatch"
        )
      }
      lambda <- attr(model, "stormbreaker_lambda", exact = TRUE)
      probability <- drop(stats::predict(
        model,
        newx = as.matrix(x),
        s = lambda,
        type = "response"
      ))
      valid <- length(probability) == nrow(x) && !anyNA(probability) &&
        all(is.finite(probability)) &&
        all(probability >= 0 & probability <= 1)
      if (!valid) {
        sb_abort(
          paste0(
            "The glmnet logistic learner produced invalid class-one ",
            "probabilities."
          ),
          "stormbreaker_invalid_glmnet_logistic_prediction"
        )
      }
      unname(probability)
    },
    extract_coefficients = function(model) {
      features <- attr(model, "stormbreaker_features", exact = TRUE)
      lambda <- attr(model, "stormbreaker_lambda", exact = TRUE)
      coefficient_matrix <- as.matrix(stats::coef(model, s = lambda))
      terms <- c("(Intercept)", features)
      coefficient <- unname(coefficient_matrix[terms, 1L])
      data.frame(
        term = terms,
        coefficient = coefficient,
        nonzero = coefficient != 0,
        stringsAsFactors = FALSE
      )
    },
    extract_state = function(model, context = list()) {
      features <- attr(model, "stormbreaker_features", exact = TRUE)
      lambda <- attr(model, "stormbreaker_lambda", exact = TRUE)
      coefficient_matrix <- as.matrix(stats::coef(model, s = lambda))
      terms <- c("(Intercept)", features)
      sb_extract_logistic_learner_state(
        learner = "glmnet_logistic",
        coefficients = data.frame(
          term = terms,
          coefficient = unname(coefficient_matrix[terms, 1L]),
          stringsAsFactors = FALSE
        ),
        params = params
      )
    },
    predict_state = sb_predict_logistic_learner_state
  )
}

# ---- model-sparse ----

# ---- abess-adapters ----

sb_abess_binary_learner_ids <- function() {
  "abess_binary_gic"
}

sb_abess_survival_learner_ids <- function() {
  "abess_cox_gic"
}

sb_abess_learner_ids <- function() {
  c(sb_abess_binary_learner_ids(), sb_abess_survival_learner_ids())
}

sb_abess_selector_ids <- function() {
  c("abess_binary_gic_nonzero", "abess_cox_gic_nonzero")
}

sb_abess_available <- function() {
  requireNamespace("abess", quietly = TRUE)
}

sb_abess_dependency <- function() {
  list(
    package = "abess",
    available = sb_abess_available()
  )
}

sb_abess_learner_params <- function(learner) {
  if (!learner %in% sb_abess_learner_ids()) {
    sb_abort(
      paste0("Unknown ABESS learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  family <- if (learner %in% sb_abess_binary_learner_ids()) {
    "binomial"
  } else {
    "cox"
  }
  list(
    family = family,
    tune_path = "sequence",
    tune_type = "gic",
    minimum_support = 0L,
    maximum_support = 10L,
    normalize = if (identical(family, "binomial")) 2L else 3L,
    lambda = 0,
    c_max = 2L,
    splicing_type = 2L,
    max_splicing_iter = 20L,
    newton = "exact",
    max_newton_iter = 10L,
    early_stop = FALSE,
    warm_start = TRUE,
    important_search = "all_predictors",
    num_threads = 1L,
    seed = 1L
  )
}

sb_abess_selector_params <- function(selector) {
  learner <- switch(selector,
    abess_binary_gic_nonzero = "abess_binary_gic",
    abess_cox_gic_nonzero = "abess_cox_gic",
    sb_abort(
      paste0("Unknown ABESS selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  )
  params <- sb_abess_learner_params(learner)
  params$minimum_support <- 2L
  params
}

sb_abess_stop_on_warning <- function(code, label, subclass) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      message <- conditionMessage(warning)
      changes_model <- grepl(
        "coerc|change to family|numerical result might be unstable",
        message,
        ignore.case = TRUE
      )
      if (!changes_model) return(invisible(NULL))
      sb_abort(
        paste0(label, " raised a warning: ", message),
        subclass,
        parent = warning
      )
    }
  )
}

sb_abess_support_path <- function(params, n_predictors, label) {
  maximum <- min(params$maximum_support, n_predictors)
  if (maximum < params$minimum_support) {
    sb_abort(
      paste0(
        label,
        " requires at least ",
        params$minimum_support,
        " predictors."
      ),
      "stormbreaker_abess_too_few_predictors"
    )
  }
  seq.int(params$minimum_support, maximum)
}

sb_fit_abess_engine <- function(x, endpoint, params, label) {
  if (!sb_abess_available()) {
    sb_abort(
      "Package `abess` is required for the requested ABESS component.",
      "stormbreaker_missing_abess"
    )
  }
  matrix_x <- as.matrix(x)
  storage.mode(matrix_x) <- "double"
  if (identical(params$family, "cox") && ncol(matrix_x) < 2L) {
    sb_abort(
      paste0(
        label,
        " requires at least 2 predictors because the native ABESS Cox ",
        "path does not preserve a one-column matrix."
      ),
      "stormbreaker_abess_cox_too_few_predictors"
    )
  }
  support_size <- sb_abess_support_path(params, ncol(matrix_x), label)
  important_search <- ncol(matrix_x)
  fit_intercept <- identical(params$family, "binomial")

  engine <- sb_abess_stop_on_warning(
    withr::with_seed(params$seed, {
      if (identical(params$family, "binomial")) {
        abess::abess(
          x = matrix_x,
          y = endpoint$outcome,
          family = "binomial",
          tune.path = params$tune_path,
          tune.type = params$tune_type,
          normalize = params$normalize,
          fit.intercept = fit_intercept,
          c.max = params$c_max,
          support.size = support_size,
          lambda = params$lambda,
          splicing.type = params$splicing_type,
          max.splicing.iter = params$max_splicing_iter,
          important.search = important_search,
          warm.start = params$warm_start,
          newton = params$newton,
          max.newton.iter = params$max_newton_iter,
          early.stop = params$early_stop,
          num.threads = params$num_threads,
          seed = params$seed
        )
      } else {
        abess::abess(
          x = matrix_x,
          y = survival::Surv(endpoint$time, endpoint$event),
          family = "cox",
          tune.path = params$tune_path,
          tune.type = params$tune_type,
          normalize = params$normalize,
          fit.intercept = fit_intercept,
          c.max = params$c_max,
          support.size = support_size,
          lambda = params$lambda,
          splicing.type = params$splicing_type,
          max.splicing.iter = params$max_splicing_iter,
          important.search = important_search,
          warm.start = params$warm_start,
          newton = params$newton,
          max.newton.iter = params$max_newton_iter,
          early.stop = params$early_stop,
          num.threads = params$num_threads,
          seed = params$seed
        )
      }
    }),
    paste0(label, " fit"),
    "stormbreaker_abess_fit_warning"
  )
  valid <- inherits(engine, "abess") &&
    is.numeric(engine$best.size) && length(engine$best.size) == 1L &&
    !is.na(engine$best.size) && engine$best.size %in% support_size
  if (!valid) {
    sb_abort(
      paste0(label, " did not return one valid GIC-selected support size."),
      "stormbreaker_invalid_abess_fit"
    )
  }
  engine$call <- NULL
  engine
}

sb_abess_coefficient_matrix <- function(engine) {
  coefficients <- as.matrix(stats::coef(
    engine,
    support.size = engine$best.size
  ))
  if (ncol(coefficients) != 1L || is.null(rownames(coefficients)) ||
    anyNA(coefficients) || any(!is.finite(coefficients))) {
    sb_abort(
      "The fitted ABESS model did not return portable coefficients.",
      "stormbreaker_invalid_abess_coefficients"
    )
  }
  coefficients
}

sb_new_abess_model <- function(learner, features, params, engine) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      selected_support_size = as.integer(engine$best.size),
      support_path = as.integer(engine$support.size),
      criterion_path = data.frame(
        support_size = as.integer(engine$support.size),
        gic = as.numeric(engine$tune.value),
        deviance = as.numeric(engine$dev),
        stringsAsFactors = FALSE
      ),
      engine = engine
    ),
    class = "stormbreaker_abess_model"
  )
}

sb_validate_abess_model_x <- function(model, x) {
  if (!inherits(model, "stormbreaker_abess_model") ||
    !model$learner %in% sb_abess_learner_ids()) {
    sb_abort(
      "The supplied model is not a fitted ABESS learner.",
      "stormbreaker_bad_abess_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted ABESS learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  x
}

sb_fit_abess_learner <- function(x, endpoint, learner, params) {
  if (learner %in% sb_abess_binary_learner_ids()) {
    inputs <- sb_validate_classical_binary_inputs(
      x,
      endpoint$outcome,
      learner
    )
    endpoint <- list(outcome = inputs$outcome)
  } else {
    inputs <- sb_validate_classical_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      learner
    )
    endpoint <- list(time = inputs$time, event = inputs$event)
  }
  engine <- sb_fit_abess_engine(inputs$x, endpoint, params, learner)
  sb_new_abess_model(learner, names(inputs$x), params, engine)
}

sb_predict_abess_learner <- function(model, x) {
  x <- sb_validate_abess_model_x(model, x)
  type <- if (model$learner %in% sb_abess_binary_learner_ids()) {
    "response"
  } else {
    "link"
  }
  prediction <- sb_abess_stop_on_warning(
    stats::predict(
      model$engine,
      as.matrix(x),
      type = type,
      support.size = model$selected_support_size
    ),
    paste0(model$learner, " prediction"),
    "stormbreaker_abess_prediction_warning"
  )
  prediction <- unname(as.numeric(prediction))
  valid <- length(prediction) == nrow(x) && !anyNA(prediction) &&
    all(is.finite(prediction))
  if (model$learner %in% sb_abess_binary_learner_ids()) {
    valid <- valid && all(prediction >= 0 & prediction <= 1)
  }
  if (!valid) {
    sb_abort(
      "The fitted ABESS learner produced invalid predictions.",
      "stormbreaker_invalid_abess_prediction"
    )
  }
  prediction
}

sb_extract_abess_coefficients <- function(model) {
  if (!inherits(model, "stormbreaker_abess_model")) {
    sb_abort(
      "The supplied model is not a fitted ABESS learner.",
      "stormbreaker_bad_abess_fit"
    )
  }
  coefficients <- sb_abess_coefficient_matrix(model$engine)
  expected_rows <- c("(intercept)", model$features)
  if (!identical(rownames(coefficients), expected_rows)) {
    sb_abort(
      "The fitted ABESS coefficient names do not match the learner schema.",
      "stormbreaker_invalid_abess_coefficients"
    )
  }
  values <- unname(as.numeric(coefficients[, 1L]))
  if (model$learner %in% sb_abess_binary_learner_ids()) {
    return(data.frame(
      term = c("(Intercept)", model$features),
      coefficient = values,
      nonzero = values != 0,
      stringsAsFactors = FALSE
    ))
  }
  predictor_values <- values[-1L]
  data.frame(
    feature = model$features,
    coefficient = predictor_values,
    nonzero = predictor_values != 0,
    direction = ifelse(
      predictor_values > 0,
      "higher_hazard",
      ifelse(predictor_values < 0, "lower_hazard", "zero")
    ),
    stringsAsFactors = FALSE
  )
}

sb_predict_abess_cox_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_risk <- sb_predict_abess_learner(model, inputs$x)
  assessment_risk <- sb_predict_abess_learner(model, x_assessment)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = "abess_cox_gic_breslow_native_link",
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_extract_abess_learner_state <- function(model, context = list()) {
  coefficients <- sb_extract_abess_coefficients(model)
  state_params <- c(
    model$params,
    list(selected_support_size = model$selected_support_size)
  )
  if (model$learner %in% sb_abess_binary_learner_ids()) {
    return(sb_extract_logistic_learner_state(
      model$learner,
      coefficients,
      params = state_params
    ))
  }
  sb_extract_cox_learner_state(
    model,
    model$learner,
    coefficients,
    sb_predict_abess_cox_survival,
    context,
    params = state_params,
    audit = list(
      support_path = model$support_path,
      criterion_path = model$criterion_path,
      support_selection = "analysis_fold_gic"
    )
  )
}

sb_predict_abess_learner_state <- function(state, x) {
  if (state$learner %in% sb_abess_binary_learner_ids()) {
    return(sb_predict_logistic_learner_state(state, x))
  }
  sb_predict_cox_learner_state(state, x)
}

sb_abess_learner_adapter <- function(
  learner,
  params = sb_abess_learner_params(learner)
) {
  if (!learner %in% sb_abess_learner_ids()) {
    sb_abort(
      paste0("Unknown ABESS learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  adapter <- list(
    id = learner,
    params = params,
    fit = if (learner %in% sb_abess_binary_learner_ids()) {
      function(x, outcome) {
        sb_fit_abess_learner(
          x,
          list(outcome = outcome),
          learner,
          params
        )
      }
    } else {
      function(x, time, event) {
        sb_fit_abess_learner(
          x,
          list(time = time, event = event),
          learner,
          params
        )
      }
    },
    predict = sb_predict_abess_learner,
    extract_coefficients = sb_extract_abess_coefficients,
    extract_state = sb_extract_abess_learner_state,
    predict_state = sb_predict_abess_learner_state
  )
  if (learner %in% sb_abess_survival_learner_ids()) {
    adapter$predict_survival <- sb_predict_abess_cox_survival
  }
  adapter
}

sb_fit_abess_selector <- function(x, endpoint, selector, capture_history = FALSE) {
  params <- sb_abess_selector_params(selector)
  learner <- if (identical(selector, "abess_binary_gic_nonzero")) {
    "abess_binary_gic"
  } else {
    "abess_cox_gic"
  }
  if (learner %in% sb_abess_binary_learner_ids()) {
    inputs <- sb_validate_classical_binary_inputs(
      x,
      endpoint$outcome,
      selector
    )
    endpoint <- list(outcome = inputs$outcome)
  } else {
    inputs <- sb_validate_classical_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      selector
    )
    endpoint <- list(time = inputs$time, event = inputs$event)
  }
  engine <- sb_fit_abess_engine(inputs$x, endpoint, params, selector)
  coefficients <- sb_abess_coefficient_matrix(engine)
  if (!identical(rownames(coefficients), c("(intercept)", names(inputs$x)))) {
    sb_abort(
      "The ABESS selector coefficient names do not match its input schema.",
      "stormbreaker_invalid_abess_selector"
    )
  }
  values <- unname(as.numeric(coefficients[-1L, 1L]))
  scores <- data.frame(
    feature = names(inputs$x),
    coefficient = values,
    absolute_coefficient = abs(values),
    selected = values != 0,
    stringsAsFactors = FALSE
  )
  selected_features <- scores$feature[scores$selected]
  if (length(selected_features) < 2L) {
    sb_abort(
      "The ABESS GIC selector retained fewer than two predictors.",
      "stormbreaker_abess_selector_too_few_features"
    )
  }
  structure(
    list(
      id = selector,
      source_features = names(inputs$x),
      selected_features = selected_features,
      selected_support_size = as.integer(engine$best.size),
      params = params,
      scores = scores,
      coefficient_path = if (capture_history) list(beta = engine$beta,
        support_size = engine$support.size) else NULL,
      criterion_path = data.frame(
        support_size = as.integer(engine$support.size),
        gic = as.numeric(engine$tune.value),
        deviance = as.numeric(engine$dev),
        stringsAsFactors = FALSE
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_apply_abess_selector <- function(fitted, x) {
  valid <- inherits(fitted, "stormbreaker_selector_fit") &&
    fitted$id %in% sb_abess_selector_ids() &&
    is.character(fitted$source_features) &&
    is.character(fitted$selected_features) &&
    length(fitted$selected_features) >= 2L &&
    !anyNA(fitted$selected_features) &&
    !anyDuplicated(fitted$selected_features) &&
    all(fitted$selected_features %in% fitted$source_features)
  if (!valid) {
    sb_abort(
      "The fitted selector is not a valid ABESS GIC selector state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      "Predictor columns do not match the fitted ABESS selector schema.",
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_abess_selector_adapter <- function(selector, capture_history = FALSE) {
  if (!selector %in% sb_abess_selector_ids()) {
    sb_abort(
      paste0("Unknown ABESS selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  }
  list(
    id = selector,
    params = list(),
    capture_history = capture_history,
    needs_outcome = TRUE,
    fit = if (identical(selector, "abess_binary_gic_nonzero")) {
      function(x, outcome) {
        sb_fit_abess_selector(x, list(outcome = outcome), selector, capture_history)
      }
    } else {
      function(x, time, event) {
        sb_fit_abess_selector(
          x,
          list(time = time, event = event),
          selector, capture_history
        )
      }
    },
    apply = sb_apply_abess_selector
  )
}

# ---- algorithm-ncvreg-path ----

sb_ncvreg_path_diagnostic <- function(fit, type) {
  model <- fit$diagnostic_artifacts$learner
  path <- model$diagnostic_path
  if (!is.list(path) || !is.numeric(path$lambda) || !length(path$lambda) ||
      any(!is.finite(path$lambda)) || any(path$lambda <= 0) ||
      !is.numeric(path$bic) || length(path$bic) != length(path$lambda) ||
      any(!is.finite(path$bic))) sb_abort(
    "The original MCP/SCAD path is missing or invalid. Older fits may lack path capture.",
    "stormbreaker_invalid_diagnostic")
  beta <- as.matrix(path$beta)
  expected <- if (model$learner %in% sb_ncvreg_binary_learner_ids()) c("(Intercept)", model$features) else model$features
  if (!is.numeric(beta) || ncol(beta) != length(path$lambda) ||
      !identical(rownames(beta), expected) || any(!is.finite(beta))) sb_abort(
    "Native MCP/SCAD coefficient path does not match the fitted features.",
    "stormbreaker_invalid_diagnostic")
  if (type == "bic_path") {
    tab <- data.frame(lambda = path$lambda, bic = path$bic)
  } else {
    tab <- data.frame(lambda = rep(path$lambda, each = nrow(beta)),
      feature = rep(rownames(beta), ncol(beta)), coefficient = as.vector(beta))
  }
  structure(list(type = type, table = tab, learner = fit$pipeline$learner,
    pipeline_id = fit$pipeline$id, engine_class = path$engine_class,
    selected_lambda = model$effective$selected_lambda,
    source = "native original-fit MCP/SCAD path; intercept explicitly labelled when present"),
    class = "stormbreaker_algorithm_diagnostic")
}

# ---- ncvreg-adapters ----

sb_ncvreg_binary_learner_ids <- function() {
  c("ncvreg_mcp_binary_bic", "ncvreg_scad_binary_bic")
}

sb_ncvreg_survival_learner_ids <- function() {
  c("ncvreg_mcp_cox_bic", "ncvreg_scad_cox_bic")
}

sb_ncvreg_learner_ids <- function() {
  c(sb_ncvreg_binary_learner_ids(), sb_ncvreg_survival_learner_ids())
}

sb_ncvreg_available <- function() {
  requireNamespace("ncvreg", quietly = TRUE)
}

sb_ncvreg_dependency <- function() {
  list(package = "ncvreg", available = sb_ncvreg_available())
}

sb_ncvreg_learner_params <- function(learner) {
  if (!learner %in% sb_ncvreg_learner_ids()) {
    sb_abort(
      paste0("Unknown ncvreg learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  penalty <- if (grepl("_mcp_", learner, fixed = TRUE)) {
    "MCP"
  } else {
    "SCAD"
  }
  list(
    family = if (learner %in% sb_ncvreg_binary_learner_ids()) {
      "binomial"
    } else {
      "cox"
    },
    penalty = penalty,
    gamma = if (identical(penalty, "MCP")) 3 else 3.7,
    alpha = 1,
    nlambda = 100L,
    lambda_min_ratio_n_gt_p = 0.001,
    lambda_min_ratio_n_le_p = 0.05,
    eps = 1e-4,
    max_iter = 10000L,
    convex_diagnostics = TRUE,
    standardization = "native_internal",
    selection = "analysis_fold_bic",
    bic_tie_rule = "first_minimum_on_descending_lambda_path",
    warning_policy = "fail_fit",
    return_x = FALSE
  )
}

sb_ncvreg_parameter_semantics <- function(learner) {
  params <- sb_ncvreg_learner_params(learner)
  paste0(
    params$penalty,
    " gamma ",
    sb_format_roundtrip_double(params$gamma),
    ", alpha 1, native standardization, 100-lambda path, and documented ",
    "lambda-min ratio rule fixed; lambda selected by BIC inside each ",
    "analysis fold"
  )
}

sb_require_ncvreg_dependency <- function() {
  if (!sb_ncvreg_available()) {
    sb_abort(
      "Package `ncvreg` is required for MCP and SCAD learners.",
      "stormbreaker_missing_ncvreg"
    )
  }
  invisible(TRUE)
}

sb_ncvreg_stop_on_warning <- function(code, learner, stage) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      sb_abort(
        paste0(
          "Learner `", learner, "` raised a ", stage, " warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_ncvreg_fit_warning",
        parent = warning
      )
    }
  )
}

sb_ncvreg_lambda_min_ratio <- function(params, n, p) {
  if (n > p) {
    params$lambda_min_ratio_n_gt_p
  } else {
    params$lambda_min_ratio_n_le_p
  }
}

sb_ncvreg_selected_coefficients <- function(engine, features, learner,
                                             selected_lambda) {
  coefficient <- drop(stats::coef(engine, lambda = selected_lambda))
  expected <- if (learner %in% sb_ncvreg_binary_learner_ids()) {
    c("(Intercept)", features)
  } else {
    features
  }
  if (is.numeric(coefficient) && length(coefficient) == length(expected) &&
    is.null(names(coefficient))) {
    names(coefficient) <- expected
  }
  valid <- is.numeric(coefficient) &&
    identical(names(coefficient), expected) &&
    !anyNA(coefficient) && all(is.finite(coefficient))
  if (!valid) {
    sb_abort(
      paste0(
        "Learner `", learner,
        "` did not return finite coefficients in the fitted feature order."
      ),
      "stormbreaker_invalid_ncvreg_coefficients"
    )
  }
  coefficient
}

sb_fit_ncvreg_engine <- function(inputs, learner, params) {
  sb_require_ncvreg_dependency()
  matrix_x <- as.matrix(inputs$x)
  storage.mode(matrix_x) <- "double"
  lambda_min <- sb_ncvreg_lambda_min_ratio(
    params,
    nrow(matrix_x),
    ncol(matrix_x)
  )
  dfmax <- ncol(matrix_x) + if (identical(params$family, "binomial")) 1L else 0L

  engine <- sb_ncvreg_stop_on_warning(
    if (identical(params$family, "binomial")) {
      ncvreg::ncvreg(
        X = matrix_x,
        y = as.numeric(inputs$outcome),
        family = "binomial",
        penalty = params$penalty,
        gamma = params$gamma,
        alpha = params$alpha,
        lambda.min = lambda_min,
        nlambda = params$nlambda,
        eps = params$eps,
        max.iter = params$max_iter,
        convex = params$convex_diagnostics,
        dfmax = dfmax,
        penalty.factor = rep(1, ncol(matrix_x)),
        warn = TRUE,
        returnX = params$return_x
      )
    } else {
      response <- cbind(
        time = as.numeric(inputs$time),
        event = as.numeric(inputs$event)
      )
      ncvreg::ncvsurv(
        X = matrix_x,
        y = response,
        penalty = params$penalty,
        gamma = params$gamma,
        alpha = params$alpha,
        lambda.min = lambda_min,
        nlambda = params$nlambda,
        eps = params$eps,
        max.iter = params$max_iter,
        convex = params$convex_diagnostics,
        dfmax = dfmax,
        penalty.factor = rep(1, ncol(matrix_x)),
        warn = TRUE,
        returnX = params$return_x
      )
    },
    learner,
    "fitting"
  )
  expected_class <- if (identical(params$family, "binomial")) {
    "ncvreg"
  } else {
    "ncvsurv"
  }
  if (!inherits(engine, expected_class)) {
    sb_abort(
      paste0("Learner `", learner, "` returned an invalid native fit."),
      "stormbreaker_invalid_ncvreg_fit"
    )
  }

  bic <- sb_ncvreg_stop_on_warning(
    stats::BIC(engine),
    learner,
    "BIC calculation"
  )
  valid_bic <- is.numeric(bic) && length(bic) == length(engine$lambda) &&
    length(bic) > 0L && !anyNA(bic) && all(is.finite(bic)) &&
    is.numeric(engine$lambda) && all(is.finite(engine$lambda)) &&
    all(engine$lambda > 0)
  if (!valid_bic) {
    sb_abort(
      paste0("Learner `", learner, "` returned an invalid BIC path."),
      "stormbreaker_invalid_ncvreg_bic"
    )
  }
  selected_index <- which.min(bic)
  selected_lambda <- unname(engine$lambda[[selected_index]])
  coefficient <- sb_ncvreg_selected_coefficients(
    engine,
    names(inputs$x),
    learner,
    selected_lambda
  )
  convex_min <- engine$convex.min
  if (!is.numeric(convex_min) || length(convex_min) != 1L ||
    is.na(convex_min) || !is.finite(convex_min)) {
    convex_min <- NA_integer_
  } else {
    convex_min <- as.integer(convex_min)
  }

  list(
    coefficient = coefficient,
    diagnostic_path = list(lambda = engine$lambda, beta = engine$beta, bic = bic,
      engine_class = class(engine)),
    effective = list(
      selected_index = as.integer(selected_index),
      selected_lambda = selected_lambda,
      selected_bic = unname(bic[[selected_index]]),
      path_length = as.integer(length(engine$lambda)),
      lambda_min_ratio = lambda_min,
      convex_min = convex_min
    )
  )
}

sb_new_ncvreg_model <- function(learner, features, params, fitted) {
  binary <- learner %in% sb_ncvreg_binary_learner_ids()
  coefficient <- fitted$coefficient
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      intercept = if (binary) unname(coefficient[[1L]]) else 0,
      coefficient = stats::setNames(
        unname(if (binary) coefficient[-1L] else coefficient),
        features
      ),
      effective = fitted$effective,
      diagnostic_path = fitted$diagnostic_path
    ),
    class = "stormbreaker_ncvreg_model"
  )
}

sb_fit_ncvreg_learner <- function(x, endpoint, learner, params) {
  inputs <- if (learner %in% sb_ncvreg_binary_learner_ids()) {
    sb_validate_classical_binary_inputs(
      x,
      endpoint$outcome,
      learner
    )
  } else {
    sb_validate_classical_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      learner
    )
  }
  fitted <- sb_fit_ncvreg_engine(inputs, learner, params)
  sb_new_ncvreg_model(learner, names(inputs$x), params, fitted)
}

sb_validate_ncvreg_model_x <- function(model, x) {
  valid <- inherits(model, "stormbreaker_ncvreg_model") &&
    model$learner %in% sb_ncvreg_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.numeric(model$coefficient) &&
    identical(names(model$coefficient), model$features) &&
    !anyNA(model$coefficient) && all(is.finite(model$coefficient))
  if (!valid) {
    sb_abort(
      "The supplied model is not a valid fitted ncvreg learner.",
      "stormbreaker_bad_ncvreg_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted ncvreg learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  x
}

sb_predict_ncvreg_learner <- function(model, x) {
  x <- sb_validate_ncvreg_model_x(model, x)
  score <- model$intercept + drop(as.matrix(x) %*% model$coefficient)
  prediction <- if (model$learner %in% sb_ncvreg_binary_learner_ids()) {
    stats::plogis(score)
  } else {
    score
  }
  valid <- length(prediction) == nrow(x) && !anyNA(prediction) &&
    all(is.finite(prediction))
  if (model$learner %in% sb_ncvreg_binary_learner_ids()) {
    valid <- valid && all(prediction >= 0 & prediction <= 1)
  }
  if (!valid) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid predictions."),
      "stormbreaker_invalid_ncvreg_prediction"
    )
  }
  unname(prediction)
}

sb_extract_ncvreg_coefficients <- function(model) {
  sb_validate_ncvreg_model_x(
    model,
    as.data.frame(
      matrix(
        0,
        nrow = 1L,
        ncol = length(model$features),
        dimnames = list(NULL, model$features)
      ),
      check.names = FALSE
    )
  )
  if (model$learner %in% sb_ncvreg_binary_learner_ids()) {
    values <- c(model$intercept, unname(model$coefficient))
    return(data.frame(
      term = c("(Intercept)", model$features),
      coefficient = values,
      nonzero = values != 0,
      stringsAsFactors = FALSE
    ))
  }
  values <- unname(model$coefficient)
  data.frame(
    feature = model$features,
    coefficient = values,
    nonzero = values != 0,
    direction = ifelse(
      values > 0,
      "higher_hazard",
      ifelse(values < 0, "lower_hazard", "zero")
    ),
    stringsAsFactors = FALSE
  )
}

sb_predict_ncvreg_cox_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  if (!model$learner %in% sb_ncvreg_survival_learner_ids()) {
    sb_abort(
      "A binary ncvreg learner cannot predict survival probability.",
      "stormbreaker_bad_ncvreg_fit"
    )
  }
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_risk <- sb_predict_ncvreg_learner(model, inputs$x)
  assessment_risk <- sb_predict_ncvreg_learner(model, x_assessment)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = paste0(model$learner, "_breslow_native_link"),
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_extract_ncvreg_state <- function(model, context = list()) {
  coefficients <- sb_extract_ncvreg_coefficients(model)
  state_params <- c(model$params, model$effective)
  if (model$learner %in% sb_ncvreg_binary_learner_ids()) {
    return(sb_extract_logistic_learner_state(
      model$learner,
      coefficients,
      params = state_params
    ))
  }
  sb_extract_cox_learner_state(
    model,
    model$learner,
    coefficients,
    sb_predict_ncvreg_cox_survival,
    context,
    params = state_params,
    audit = list(lambda_selection = "analysis_fold_bic")
  )
}

sb_predict_ncvreg_state <- function(state, x) {
  if (state$learner %in% sb_ncvreg_binary_learner_ids()) {
    return(sb_predict_logistic_learner_state(state, x))
  }
  sb_predict_cox_learner_state(state, x)
}

sb_ncvreg_learner_adapter <- function(
  learner,
  params = sb_ncvreg_learner_params(learner)
) {
  if (!learner %in% sb_ncvreg_learner_ids()) {
    sb_abort(
      paste0("Unknown ncvreg learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  adapter <- list(
    id = learner,
    params = params,
    fit = if (learner %in% sb_ncvreg_binary_learner_ids()) {
      function(x, outcome) {
        sb_fit_ncvreg_learner(
          x,
          list(outcome = outcome),
          learner,
          params
        )
      }
    } else {
      function(x, time, event) {
        sb_fit_ncvreg_learner(
          x,
          list(time = time, event = event),
          learner,
          params
        )
      }
    },
    predict = sb_predict_ncvreg_learner,
    extract_coefficients = sb_extract_ncvreg_coefficients,
    extract_state = sb_extract_ncvreg_state,
    predict_state = sb_predict_ncvreg_state
  )
  if (learner %in% sb_ncvreg_survival_learner_ids()) {
    adapter$predict_survival <- sb_predict_ncvreg_cox_survival
  }
  adapter
}

# ---- model-slope-randomglm ----

# ---- algorithm-randomglm ----

#' @rdname sb_model_coefficients
#' @export
sb_randomglm_diagnostics <- function(fit, type = "inclusion_frequency") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (fit$pipeline$learner != "random_glm_binary") sb_abort(
    "This diagnostic requires randomGLM.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("inclusion_frequency", "bag_coefficients", "bag_size")) sb_abort(
    "Choose inclusion_frequency, bag_coefficients or bag_size.", "stormbreaker_unsupported_diagnostic")
  state <- fit$learner_state
  n <- length(state$intercepts)
  features <- state$features
  if (!n || length(state$bag_features) != n || length(state$bag_coefficients) != n) sb_abort(
    "Saved randomGLM ensemble is incomplete.", "stormbreaker_bad_learner_state")
  coefficient <- matrix(0, nrow = n, ncol = length(features))
  included <- matrix(FALSE, nrow = n, ncol = length(features))
  for (i in seq_len(n)) {
    index <- state$bag_features[[i]]
    beta <- state$bag_coefficients[[i]]
    if (length(index) != length(beta) || anyNA(index) || anyDuplicated(index) ||
        any(index < 1L | index > length(features)) || any(!is.finite(beta))) sb_abort(
      "Saved randomGLM bag terms are invalid.", "stormbreaker_bad_learner_state")
    coefficient[i, index] <- beta
    included[i, index] <- TRUE
  }
  if (type == "inclusion_frequency") {
    tab <- data.frame(feature = features, bags = colSums(included),
      fraction = colMeans(included), total_bags = n)
    interpretation <- "Feature inclusion in original ensemble bags; not CV selection stability, uncertainty, causal effect or permutation importance."
  } else if (type == "bag_coefficients") {
    tab <- data.frame(bag = rep(seq_len(n), length(features)), feature = rep(features, each = n),
      coefficient = as.vector(coefficient), included = as.vector(included))
    interpretation <- paste("Original per-bag log-odds coefficients in input units; omitted features have zero contribution and remain flagged.",
      "Across-bag spread is not a confidence interval. The ensemble averages probabilities, not these coefficients.")
  } else {
    tab <- data.frame(bag = seq_len(n), selected_features = rowSums(included))
    interpretation <- "Number of selected features in each original ensemble bag, including intercept-only bags; not learning progress or patient counts."
  }
  structure(list(type = paste0("randomglm_", type), table = tab, features = features,
    intercepts = state$intercepts, total_bags = n, learner = fit$pipeline$learner,
    pipeline_id = fit$pipeline$id, source = "original saved randomGLM ensemble members",
    interpretation = interpretation), class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_randomglm_diagnostic <- function(x) {
  feature <- fraction <- coefficient <- included <- selected_features <- NULL
  tab <- x$table
  if (x$type != "randomglm_bag_size") tab$feature <- factor(tab$feature, levels = rev(x$features))
  if (x$type == "randomglm_inclusion_frequency") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(fraction, feature)) +
      ggplot2::geom_col(fill = "#2878A0", width = .65) +
      ggplot2::scale_x_continuous(limits = c(0, 1), labels = function(z) paste0(round(100 * z), "%")) +
      ggplot2::labs(title = "randomGLM feature inclusion", x = "Ensemble bags including feature", y = NULL)
  } else if (x$type == "randomglm_bag_coefficients") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(coefficient, feature, colour = included)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey75", linewidth = .4) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0, height = .18, seed = 1),
        size = 1.3, alpha = .65) +
      ggplot2::scale_colour_manual(values = c("FALSE" = "#9AA6AF", "TRUE" = "#2878A0"),
        breaks = c("TRUE", "FALSE"), labels = c("Included", "Omitted")) +
      ggplot2::labs(title = "randomGLM bag coefficients", x = "Bag coefficient (log odds per input unit)",
        y = NULL, colour = NULL)
  } else {
    p <- ggplot2::ggplot(tab, ggplot2::aes(selected_features)) +
      ggplot2::geom_bar(fill = "#2878A0", width = .65) +
      ggplot2::scale_x_continuous(breaks = sort(unique(tab$selected_features))) +
      ggplot2::labs(title = "randomGLM bag sizes", x = "Selected predictors", y = "Number of bags")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = if (x$type == "randomglm_bag_size") 5 else max(5, 2 + .3 * length(x$features)))
  p
}

# ---- algorithm-slope ----

#' @rdname sb_model_coefficients
#' @export
sb_slope_diagnostics <- function(fit) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (fit$pipeline$learner != "slope_binary") sb_abort(
    "This diagnostic requires SLOPE.", "stormbreaker_unsupported_diagnostic")
  path <- fit$diagnostic_artifacts$learner$diagnostic_path
  if (is.null(path)) sb_abort("The original SLOPE CV summary was not captured.",
    "stormbreaker_missing_diagnostic_artifacts")
  tab <- path$summary
  columns <- c("q", "gamma", "alpha", "measure", "mean", "se", "lo", "hi")
  if (!is.data.frame(tab) || !all(columns %in% names(tab)) || !nrow(tab) ||
      any(!is.finite(as.matrix(tab[setdiff(columns, "measure")]))) || any(tab$alpha <= 0) ||
      any(tab$se < 0) || any(tab$lo > tab$mean | tab$hi < tab$mean) ||
      !identical(unique(tab$measure), "deviance")) sb_abort(
    "The original SLOPE deviance path is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
  optimum <- path$optima[1L, , drop = FALSE]
  tab$selected <- tab$q == optimum$q & tab$gamma == optimum$gamma & tab$alpha == optimum$alpha
  if (sum(tab$selected) != 1L || !isTRUE(all.equal(tab$mean[tab$selected], fit$learner_state$tuning$mean))) sb_abort(
    "The native SLOPE optimum does not match the saved prediction model.", "stormbreaker_invalid_diagnostic")
  structure(list(type = "slope_cv_curve", table = tab, native_table = path$summary,
    selected = optimum, score_rows = path$score_rows,
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    source = "original final-fit SLOPE cross-validation summary",
    interpretation = paste("Original mean CV deviance and native t-based 95% bounds; lower is better.",
      "Orange point identifies the original chosen q/gamma/alpha, not a new selection.",
      "Native early-stopped paths retained; not external-validation uncertainty or a coefficient path.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_slope_diagnostic <- function(x) {
  alpha <- mean <- lo <- hi <- NULL
  tab <- x$table
  p <- ggplot2::ggplot(tab, ggplot2::aes(alpha, mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "#2878A0", alpha = .16) +
    ggplot2::geom_line(colour = "#2878A0") +
    ggplot2::geom_point(data = tab[tab$selected, , drop = FALSE], colour = "#D77935", size = 2.5) +
    ggplot2::facet_wrap(~q + gamma, labeller = ggplot2::label_both) +
    ggplot2::scale_x_log10() + ggplot2::labs(title = "SLOPE cross-validation", x = "Alpha (log scale)",
      y = "CV deviance") + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = max(5, 2 + 2.5 * ceiling(nrow(unique(tab[c("q", "gamma")])) / 3)))
  p
}

# ---- slope-randomglm-binary-adapters ----


sb_slope_binary_learner_ids <- function() character()

sb_random_glm_binary_learner_ids <- function() {
  "random_glm_binary"
}

sb_slope_randomglm_binary_learner_ids <- function() {
  c(sb_slope_binary_learner_ids(), sb_random_glm_binary_learner_ids())
}

sb_slope_randomglm_binary_params <- function(learner) {
  switch(learner,
    slope_binary = list(
      q = c(0.05, 0.1, 0.2),
      gamma = 0,
      inner_folds = 5L,
      repeats = 1L,
      measure = "deviance",
      center = "mean",
      scale = "l2",
      path_length = 30L,
      threads = 1L,
      seed = 1L
    ),
    random_glm_binary = list(
      n_bags = 100L,
      n_candidate_covariates = 3L,
      max_interaction_order = 1L,
      replace = TRUE,
      threads = 1L,
      seed = 1L
    ),
    sb_abort(
      paste0("Unknown SLOPE/randomGLM binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_slope_randomglm_binary_dependency <- function(learner) {
  package <- switch(learner,
    slope_binary = "SLOPE",
    random_glm_binary = "randomGLM",
    sb_abort(
      paste0("Unknown SLOPE/randomGLM binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
  list(
    package = package,
    available = requireNamespace(package, quietly = TRUE)
  )
}

sb_require_slope_randomglm_binary_dependency <- function(learner) {
  dependency <- sb_slope_randomglm_binary_dependency(learner)
  if (!dependency$available) {
    sb_abort(
      paste0(
        "Package `", dependency$package, "` is required for learner `",
        learner, "`."
      ),
      if (identical(learner, "slope_binary")) {
        "stormbreaker_missing_slope"
      } else {
        "stormbreaker_missing_randomglm"
      }
    )
  }
  invisible(TRUE)
}

sb_slope_randomglm_binary_parameter_semantics <- function(learner) {
  params <- sb_slope_randomglm_binary_params(learner)
  switch(learner,
    slope_binary = paste0(
      "SLOPE binomial path with q in {",
      paste(params$q, collapse = ","), "}, gamma=", params$gamma,
      ", mean centering, L2 scaling, and ", params$path_length,
      " alpha values; q and alpha selected by native ",
      params$inner_folds, "-fold deviance CV inside each analysis fold with ",
      "seed ", params$seed
    ),
    random_glm_binary = paste0(
      "randomGLM ensemble with ", params$n_bags,
      " bootstrap bags, the native random-subspace rule, at most ",
      params$n_candidate_covariates,
      " forward-selection candidates per bag, no interactions, serial ",
      "execution, and seed ", params$seed, " fixed before resampling"
    )
  )
}

sb_new_slope_binary_model <- function(features, params, intercept, coefficient,
                                      tuning) {
  structure(
    list(
      learner = "slope_binary",
      features = features,
      params = params,
      intercept = intercept,
      coefficient = coefficient,
      tuning = tuning
    ),
    class = "stormbreaker_slope_binary_model"
  )
}

sb_fit_slope_binary <- function(inputs, params) {
  sb_reject_retired_learner("slope_binary")
}


sb_compact_random_glm_binary <- function(fit, features, params) {
  intercepts <- unname(as.numeric(fit$interceptOfForwardRegression))
  bag_coefficients <- unname(fit$coefOfForwardRegression)
  bag_features <- lapply(fit$featuresInForwardRegression, function(indices) {
    as.integer(indices)
  })
  valid <- length(intercepts) == params$n_bags &&
    length(bag_coefficients) == params$n_bags &&
    length(bag_features) == params$n_bags &&
    !anyNA(intercepts) && all(is.finite(intercepts))
  if (valid) {
    valid <- all(vapply(seq_len(params$n_bags), function(index) {
      coefficients <- unname(as.numeric(bag_coefficients[[index]]))
      indices <- bag_features[[index]]
      length(coefficients) == length(indices) &&
        !anyNA(coefficients) && all(is.finite(coefficients)) &&
        !anyNA(indices) && all(indices >= 1L & indices <= length(features))
    }, logical(1)))
  }
  if (!isTRUE(valid)) {
    sb_abort(
      "randomGLM returned an invalid bagged logistic ensemble.",
      "stormbreaker_invalid_randomglm_fit"
    )
  }
  structure(
    list(
      learner = "random_glm_binary",
      features = features,
      params = params,
      intercepts = intercepts,
      bag_features = bag_features,
      bag_coefficients = lapply(bag_coefficients, function(x) {
        unname(as.numeric(x))
      })
    ),
    class = "stormbreaker_random_glm_binary_model"
  )
}

sb_fit_random_glm_binary <- function(inputs, params) {
  fit <- withr::with_seed(
    params$seed,
    randomGLM::randomGLM(
      x = as.matrix(inputs$x),
      y = factor(inputs$outcome, levels = c(0L, 1L)),
      type = "binary",
      maxInteractionOrder = params$max_interaction_order,
      nBags = params$n_bags,
      replace = params$replace,
      nCandidateCovariates = min(
        params$n_candidate_covariates,
        ncol(inputs$x)
      ),
      keepModels = FALSE,
      randomSeed = params$seed,
      nThreads = params$threads,
      verbose = 0
    )
  )
  sb_compact_random_glm_binary(fit, names(inputs$x), params)
}

sb_fit_slope_randomglm_binary <- function(x, outcome, learner, params) {
  sb_reject_retired_learner(learner)
  sb_require_slope_randomglm_binary_dependency(learner)
  inputs <- sb_validate_classical_binary_inputs(x, outcome, learner)
  if (identical(learner, "slope_binary")) {
    return(sb_fit_slope_binary(inputs, params))
  }
  sb_fit_random_glm_binary(inputs, params)
}

sb_validate_slope_randomglm_binary_x <- function(model, x) {
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      paste0("Prediction columns do not match learner `", model$learner, "`."),
      "stormbreaker_learner_schema_mismatch"
    )
  }
  x
}

sb_predict_slope_randomglm_binary <- function(model, x) {
  x <- sb_validate_slope_randomglm_binary_x(model, x)
  matrix_x <- as.matrix(x)
  probability <- if (inherits(model, "stormbreaker_slope_binary_model")) {
    stats::plogis(model$intercept + matrix_x %*% model$coefficient)
  } else if (inherits(model, "stormbreaker_random_glm_binary_model")) {
    bag_probabilities <- lapply(seq_along(model$intercepts), function(index) {
      selected <- model$bag_features[[index]]
      linear_predictor <- rep(model$intercepts[[index]], nrow(matrix_x))
      if (length(selected) > 0L) {
        linear_predictor <- linear_predictor + as.numeric(
          matrix_x[, selected, drop = FALSE] %*%
            model$bag_coefficients[[index]]
        )
      }
      stats::plogis(linear_predictor)
    })
    Reduce(`+`, bag_probabilities) / length(bag_probabilities)
  } else {
    sb_abort(
      "The supplied SLOPE/randomGLM binary model is invalid.",
      "stormbreaker_bad_slope_randomglm_fit"
    )
  }
  sb_class_one_probability(probability, nrow(x), model$learner)
}

sb_extract_slope_binary_coefficients <- function(model) {
  data.frame(
    term = c("(Intercept)", model$features),
    coefficient = c(model$intercept, unname(model$coefficient)),
    nonzero = c(model$intercept, unname(model$coefficient)) != 0,
    stringsAsFactors = FALSE
  )
}

sb_extract_slope_randomglm_binary_state <- function(model, context = list()) {
  state <- unclass(model)
  state$diagnostic_path <- NULL
  class(state) <- c(
    "stormbreaker_slope_randomglm_binary_state",
    "stormbreaker_learner_state"
  )
  state
}

sb_restore_slope_randomglm_binary_model <- function(state) {
  if (identical(state$learner, "slope_binary")) {
    return(sb_new_slope_binary_model(
      state$features,
      state$params,
      state$intercept,
      state$coefficient,
      state$tuning
    ))
  }
  if (identical(state$learner, "random_glm_binary")) {
    model <- state
    class(model) <- "stormbreaker_random_glm_binary_model"
    return(model)
  }
  sb_bad_learner_state("The saved SLOPE/randomGLM binary state is invalid.")
}

sb_predict_slope_randomglm_binary_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_slope_randomglm_binary_learner_ids()) {
    sb_bad_learner_state("The saved SLOPE/randomGLM binary state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_restore_slope_randomglm_binary_model(state)
  list(probability = sb_predict_slope_randomglm_binary(model, x))
}

sb_slope_randomglm_binary_learner_adapter <- function(
  learner,
  params = sb_slope_randomglm_binary_params(learner)
) {
  adapter <- list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_slope_randomglm_binary(x, outcome, learner, params)
    },
    predict = sb_predict_slope_randomglm_binary,
    extract_state = sb_extract_slope_randomglm_binary_state,
    predict_state = sb_predict_slope_randomglm_binary_state
  )
  if (identical(learner, "slope_binary")) {
    adapter$extract_coefficients <- sb_extract_slope_binary_coefficients
  }
  adapter
}

# ---- diagnostics-linear ----

# ---- algorithm-discriminant ----

sb_discriminant_parameters <- function(state) {
  id <- state$learner
  if (!id %in% c("lda_binary", "qda_binary", "rda_binary")) sb_abort(
    "This learner has no supported discriminant parameter adapter.", "stormbreaker_unsupported_diagnostic")
  if (sb_is_constant_sum_state(state)) return(sb_constant_sum_linear_parameters(state))
  if (!is.raw(state$engine_raw) || !length(state$engine_raw)) sb_abort(
    "Saved discriminant parameters are missing.", "stormbreaker_bad_learner_state")
  e <- unserialize(state$engine_raw)
  expected <- switch(id, lda_binary = "lda", qda_binary = "qda", rda_binary = "rda")
  features <- state$features
  p <- length(features)
  classes <- if (id == "rda_binary") e$classes else e$lev
  means <- if (id == "rda_binary") t(e$means) else e$means
  if (!inherits(e, expected) || !identical(classes, c("0", "1")) ||
      !is.matrix(means) || !identical(dim(means), c(2L, as.integer(p))) ||
      !identical(colnames(means), features) || any(!is.finite(means)) ||
      !identical(names(e$prior), classes) || any(!is.finite(e$prior)) || any(e$prior <= 0)) sb_abort(
    "Original discriminant class means or priors are invalid.", "stormbreaker_bad_learner_state")
  prior <- e$prior
  covariance <- NULL
  if (id == "lda_binary") {
    if (!is.matrix(e$scaling) || nrow(e$scaling) != p || !ncol(e$scaling) ||
        !identical(rownames(e$scaling), features) || any(!is.finite(e$scaling))) sb_abort(
      "Original LDA scaling is invalid.", "stormbreaker_bad_learner_state")
    center <- colSums(prior * means)
    projected <- sweep(means, 2L, center) %*% e$scaling
    beta <- drop(e$scaling %*% (projected[2L, ] - projected[1L, ]))
    intercept <- log(prior[2L] / prior[1L]) - sum(center * beta) -
      0.5 * (sum(projected[2L, ]^2) - sum(projected[1L, ]^2))
    quadratic <- matrix(0, p, p)
  } else {
    if (id == "qda_binary") {
      precision <- lapply(seq_len(2L), function(i) tcrossprod(matrix(e$scaling[, , i], p, p)))
      covariance <- lapply(precision, solve)
      logdet <- e$ldet
    } else {
      gamma <- e$regularization[1L]
      lambda <- e$regularization[2L]
      if (length(c(gamma, lambda)) != 2L || any(!is.finite(c(gamma, lambda))) ||
          any(c(gamma, lambda) < 0 | c(gamma, lambda) > 1)) sb_abort(
        "Original RDA regularization is invalid.", "stormbreaker_bad_learner_state")
      covariance <- lapply(seq_len(2L), function(i) {
        pooled <- (1 - lambda) * matrix(e$covariances[, , i], p, p) + lambda * e$covpooled
        (1 - gamma) * pooled + gamma * mean(diag(pooled)) * diag(p)
      })
      precision <- lapply(covariance, solve)
      logdet <- vapply(covariance, function(z) as.numeric(determinant(z, logarithm = TRUE)$modulus), 0)
    }
    quadratic <- 0.5 * (precision[[1L]] - precision[[2L]])
    beta <- drop(precision[[2L]] %*% means[2L, ] - precision[[1L]] %*% means[1L, ])
    intercept <- log(prior[2L] / prior[1L]) - 0.5 * (logdet[2L] - logdet[1L] +
      sum(means[2L, ] * (precision[[2L]] %*% means[2L, ])) -
      sum(means[1L, ] * (precision[[1L]] %*% means[1L, ])))
    covariance <- lapply(covariance, function(z) {
      dimnames(z) <- list(features, features)
      z
    })
    names(covariance) <- classes
  }
  if (any(!is.finite(c(beta, intercept, quadratic)))) sb_abort(
    "Original discriminant equation is not finite.", "stormbreaker_bad_learner_state")
  names(beta) <- features
  dimnames(quadratic) <- list(features, features)
  list(coefficient = beta, intercept = unname(intercept), quadratic = quadratic,
    means = means, priors = prior, covariance = covariance)
}

#' @rdname sb_pls_diagnostics
#' @export
sb_discriminant_diagnostics <- function(fit, type = "class_means") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("class_means", "class_correlation")) sb_abort(
    "Choose class_means or class_correlation.", "stormbreaker_unsupported_diagnostic")
  parameters <- sb_discriminant_parameters(fit$learner_state)
  features <- fit$learner_state$features
  if (type == "class_means") {
    tab <- data.frame(class = rep(c("0", "1"), length(features)),
      feature = rep(features, each = 2L), value = as.vector(parameters$means))
    interpretation <- "Original fitted class means in input units; not uncertainty intervals or variable importance."
  } else {
    if (is.null(parameters$covariance)) sb_abort(
      "Rank-reduced LDA scaling does not retain a recoverable full within-class covariance.",
      "stormbreaker_unavailable_diagnostic")
    tab <- do.call(rbind, lapply(c("0", "1"), function(class) {
      correlation <- stats::cov2cor(parameters$covariance[[class]])
      data.frame(class = class, feature = rep(features, length(features)),
        feature2 = rep(features, each = length(features)), value = as.vector(correlation))
    }))
    interpretation <- paste("Correlations implied by the original class-specific prediction covariance.",
      "RDA includes its original pooling and spherical shrinkage; not raw-data correlation, importance or SHAP interaction.")
  }
  structure(list(type = paste0("discriminant_", type), table = tab, features = features,
    parameters = parameters, learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    source = "original saved discriminant model parameters", interpretation = interpretation),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_discriminant_diagnostic <- function(x) {
  feature <- feature2 <- class <- value <- NULL
  tab <- x$table
  tab$feature <- factor(tab$feature, levels = x$features)
  if (x$type == "discriminant_class_means") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(value, class, colour = class)) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::facet_wrap(~feature, scales = "free_x", ncol = min(3L, length(x$features))) +
      ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#D77935"), guide = "none") +
      ggplot2::labs(title = "Fitted class means", x = "Mean (original units)", y = "Encoded class") + sb_viz_theme()
    size <- c(width = 11, height = max(4, 2 + 2 * ceiling(length(x$features) / 3)))
  } else {
    tab$feature2 <- factor(tab$feature2, levels = rev(x$features))
    p <- ggplot2::ggplot(tab, ggplot2::aes(feature, feature2, fill = value)) +
      ggplot2::geom_tile() + ggplot2::facet_wrap(~class, nrow = 1L, labeller = ggplot2::label_both) +
      ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#D77935", limits = c(-1, 1),
        guide = ggplot2::guide_colourbar(display = "rectangles")) +
      ggplot2::coord_fixed() + ggplot2::labs(title = "Fitted within-class correlations", x = NULL, y = NULL,
        fill = "Correlation") + sb_viz_theme() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    size <- c(width = max(11, 4 + .6 * length(x$features)), height = max(5, 2 + .3 * length(x$features)))
  }
  attr(p, "stormbreaker_suggested_size") <- size
  p
}

# ---- algorithm-fda ----

#' @rdname sb_pls_diagnostics
#' @export
sb_fda_diagnostics <- function(fit, type = "centers", newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% c("fda_binary", "pda_binary", "mda_binary") ||
      !is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("centers", "query_variates", "mixture_weights")) {
    sb_abort("Choose centers, query_variates or mixture_weights for FDA/PDA/MDA.",
      "stormbreaker_unsupported_diagnostic")
  }
  if (!requireNamespace("mda", quietly = TRUE)) sb_abort(
    "These diagnostics require mda.", "stormbreaker_unavailable_diagnostic")
  e <- unserialize(s$engine_raw)
  mixture <- s$learner == "mda_binary"
  if (!inherits(e, if (mixture) "mda" else "fda") || !is.matrix(e$means) ||
      any(!is.finite(e$means))) sb_abort(
    "Original discriminant centers are invalid.", "stormbreaker_invalid_diagnostic")
  mapping <- if (mixture) do.call(rbind, lapply(names(e$assign.theta), function(cls) {
    data.frame(center = e$assign.theta[[cls]], class = cls,
      component = seq_along(e$assign.theta[[cls]]))
  })) else data.frame(center = seq_len(nrow(e$means)), class = rownames(e$means), component = 1L)
  mapping <- mapping[order(mapping$center), ]
  if (!identical(as.integer(mapping$center), seq_len(nrow(e$means)))) sb_abort(
    "Original subclass assignments are invalid.", "stormbreaker_invalid_diagnostic")
  if (type == "centers") {
    tab <- do.call(rbind, lapply(seq_len(ncol(e$means)), function(j) {
      data.frame(mapping, axis = j, coordinate = e$means[, j])
    }))
  } else if (type == "mixture_weights") {
    if (!mixture) sb_abort("Only MDA has fitted mixture subclasses.", "stormbreaker_unavailable_diagnostic")
    tab <- mapping
    tab$weight <- vapply(seq_len(nrow(tab)), function(i) {
      unname(e$sub.prior[[tab$class[i]]][tab$component[i]])
    }, numeric(1))
  } else {
    if (is.null(newdata)) sb_abort("Supply explicit query rows for discriminant variates.", "stormbreaker_unavailable_diagnostic")
    x <- sb_select_final_feature_frame(newdata, s$features, "Query data", "query fitted features")
    if (!nrow(x)) sb_abort("Queries must contain at least one row.", "stormbreaker_invalid_diagnostic")
    m <- sb_predict_fda_native(e, x, type = "variates", dimension = e$dimension)
    if (!is.matrix(m) || nrow(m) != nrow(x) || any(!is.finite(m))) sb_abort(
      "Native query variates are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(query = rep(seq_len(nrow(m)), ncol(m)),
      axis = rep(seq_len(ncol(m)), each = nrow(m)), coordinate = as.numeric(m))
  }
  structure(list(type = paste0("fda_", type), table = tab,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original mda-package centers, subclass priors and native variates",
    interpretation = "Discriminant axes are not input features. MDA components are latent model structure, not clinical subgroups. Axis signs are basis conventions. No refitting or new projection."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_fda_diagnostic <- function(x) {
  coordinate <- class <- component <- weight <- query <- NULL
  if (x$type == "fda_centers") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(coordinate, factor(component), colour = class)) +
      ggplot2::geom_point(size = 2.5) + ggplot2::facet_wrap(~axis, labeller = ggplot2::label_both) +
      ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#B86236")) +
      ggplot2::labs(title = "Discriminant centers", x = "Discriminant coordinate", y = "Component", colour = "Class")
  } else if (x$type == "fda_mixture_weights") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(factor(component), weight)) +
      ggplot2::geom_col(fill = "#2878A0", width = .6) +
      ggplot2::facet_wrap(~class, labeller = ggplot2::label_both) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = "Within-class mixture weights", x = "Component", y = "Weight")
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(coordinate, factor(query))) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::facet_wrap(~axis, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "Discriminant query coordinates", x = "Discriminant coordinate", y = "Query")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}

# ---- algorithm-linear-margin ----

sb_linear_margin_diagnostic_ids <- function() {
  c("liblinear_l2_l2loss_svc_binary", "liblinear_l2_l1loss_svc_binary",
    "liblinear_l1_l2loss_svc_binary", "dwd_binary")
}

sb_linear_margin_diagnostic <- function(fit, type) {
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("coefficients", "decision_distribution", "support_counts")) sb_abort(
    "Choose coefficients or decision_distribution for this linear margin model.", "stormbreaker_unsupported_diagnostic")
  if (type == "support_counts") sb_abort(
    "This engine does not retain an e1071-style support-vector set; counts are unavailable.",
    "stormbreaker_unavailable_diagnostic")
  state <- fit$learner_state
  features <- state$features
  id <- fit$pipeline$learner
  if (id == "dwd_binary") {
    e <- state$engine
    beta <- e$coefficient
    intercept <- e$intercept
    native <- list(coefficient = beta, intercept = intercept)
    orientation <- 1
  } else {
    e <- unserialize(state$engine_raw)
    m <- e$model
    if (!inherits(m, "LiblineaR") || !m$Type %in% c(1L, 3L, 5L) ||
        !setequal(m$ClassNames, c(0L, 1L)) || !is.matrix(m$W) ||
        !identical(dim(m$W), c(1L, length(features) + 1L)) ||
        !identical(names(e$center), features) || !identical(names(e$scale), features) ||
        any(!is.finite(e$scale)) || any(e$scale <= 0) || any(!is.finite(e$center)) ||
        length(m$Bias) != 1L || !is.finite(m$Bias) || m$Bias < 0) sb_abort(
      "Original LIBLINEAR weights, bias or scaling are invalid.", "stormbreaker_invalid_diagnostic")
    orientation <- if (m$ClassNames[1L] == 1L) 1 else -1
    w <- as.numeric(m$W[1L, seq_along(features)])
    beta <- stats::setNames(orientation * w / e$scale, features)
    intercept <- orientation * (m$W[1L, length(features) + 1L] * m$Bias - sum(w * e$center / e$scale))
    native <- list(W = m$W, Bias = m$Bias, ClassNames = m$ClassNames,
      center = e$center, scale = e$scale)
  }
  if (!is.numeric(beta) || !identical(names(beta), features) || any(!is.finite(beta)) ||
      length(intercept) != 1L || !is.finite(intercept)) sb_abort(
    "Original linear decision function is invalid.", "stormbreaker_invalid_diagnostic")
  if (type == "coefficients") {
    tab <- data.frame(feature = features, coefficient = as.numeric(beta))
    interpretation <- "Exact original-unit decision weights toward class 1; intercept returned separately.\nNot log-odds, odds ratios, confidence intervals or feature importance; no refit."
  } else {
    a <- fit$diagnostic_artifacts
    if (is.null(a$predictors)) sb_abort("Development predictors were not captured; enable capture_diagnostics at final fitting.",
      "stormbreaker_missing_diagnostic_artifacts")
    if (!identical(names(a$predictors), features) || length(a$endpoint$outcome) != nrow(a$predictors) ||
        anyNA(a$endpoint$outcome) || !all(a$endpoint$outcome %in% c(0, 1))) sb_abort(
      "Captured predictor order or outcomes are invalid.", "stormbreaker_invalid_diagnostic")
    decision <- as.numeric(data.matrix(a$predictors) %*% beta + intercept)
    if (any(!is.finite(decision))) sb_abort("Original linear decision scores are nonfinite.",
      "stormbreaker_invalid_diagnostic")
    tab <- data.frame(observation = seq_along(decision), class = as.character(a$endpoint$outcome), decision = decision)
    native <- NULL
    mapping <- if (id == "dwd_binary") "plogis(2 * score)" else "plogis(score)"
    interpretation <- paste0("Native decision values for development samples, oriented toward class 1.\n",
      "Zero marks the class boundary. The package transforms the decision values with the fixed function ",
      mapping, ".")
  }
  structure(list(type = paste0("svm_", type), table = tab, native_table = native,
    intercept = unname(intercept), orientation = orientation, learner = id, pipeline_id = fit$pipeline$id,
    sensitive_training_data = type == "decision_distribution", interpretation = interpretation,
    source = "original saved linear margin model; no fitting"), class = "stormbreaker_algorithm_diagnostic")
}

# ---- algorithm-shrinkage-discriminant ----

#' @rdname sb_pls_diagnostics
#' @export
sb_shrinkage_discriminant_diagnostics <- function(fit, type = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  is_sda <- identical(s$learner, "sda_binary")
  if (!s$learner %in% c("sda_binary", "slda_binary")) sb_abort(
    "This diagnostic requires SDA or SLDA.", "stormbreaker_unsupported_diagnostic")
  choices <- if (is_sda) c("coefficients", "regularization") else c("loadings", "class_means")
  if (is.null(type)) type <- choices[1L]
  if (!is.character(type) || length(type) != 1L || is.na(type) || !type %in% choices) sb_abort(
    paste("Choose", paste(choices, collapse = " or ")), "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  intercept <- NULL
  if (is_sda) {
    if (!inherits(e, "sda") || !identical(colnames(e$beta), s$features) ||
        !identical(rownames(e$beta), c("0", "1")) || any(!is.finite(e$beta))) sb_abort(
      "Original SDA coefficients are invalid.", "stormbreaker_invalid_diagnostic")
    if (type == "coefficients") {
      tab <- data.frame(feature = s$features, coefficient = as.numeric(e$beta["1", ] - e$beta["0", ]))
      intercept <- unname(e$alpha["1"] - e$alpha["0"])
    } else {
      tab <- data.frame(parameter = names(e$regularization), value = as.numeric(e$regularization))
    }
  } else {
    if (!inherits(e, "slda") || !is.matrix(e$scores) || nrow(e$scores) != length(s$features) ||
        any(!is.finite(e$scores)) || !identical(e$mylda$lev, c("0", "1"))) sb_abort(
      "Original SLDA subspace is invalid.", "stormbreaker_invalid_diagnostic")
    if (type == "loadings") {
      tab <- data.frame(feature = rep(s$features, ncol(e$scores)),
        axis = rep(seq_len(ncol(e$scores)), each = nrow(e$scores)), loading = as.numeric(e$scores))
    } else {
      m <- e$mylda$means
      tab <- data.frame(class = rep(rownames(m), ncol(m)),
        axis = rep(seq_len(ncol(m)), each = nrow(m)), mean = as.numeric(m))
    }
  }
  structure(list(type = paste0("shrinkage_discriminant_", type), table = tab,
    intercept = intercept, learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original SDA alpha/beta/shrinkage or SLDA scores/mylda",
    interpretation = "SDA coefficients describe pre-rounding class-1 log odds; native probabilities remain unchanged. SLDA loading signs are basis conventions, not importance or clinical effects. No refitting."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_shrinkage_discriminant_diagnostic <- function(x) {
  feature <- coefficient <- parameter <- value <- axis <- loading <- class <- mean <- NULL
  if (x$type == "shrinkage_discriminant_coefficients") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(coefficient, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey65", linewidth = .3) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::labs(title = "SDA class-1 coefficients", x = "Log-odds coefficient", y = NULL)
  } else if (x$type == "shrinkage_discriminant_regularization") {
    labels <- c(lambda = "Correlation", lambda.var = "Variance", lambda.freqs = "Class frequency")
    p <- ggplot2::ggplot(x$table, ggplot2::aes(value, parameter)) +
      ggplot2::geom_col(fill = "#2878A0", width = .6) +
      ggplot2::scale_y_discrete(labels = labels) + ggplot2::scale_x_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = "SDA shrinkage", x = "Shrinkage intensity", y = NULL)
  } else if (x$type == "shrinkage_discriminant_loadings") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(loading, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey65", linewidth = .3) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::facet_wrap(~axis, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "SLDA subspace loadings", x = "Loading", y = NULL)
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(mean, class)) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::facet_wrap(~axis, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "SLDA class centers", x = "Subspace coordinate", y = "Class")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 10, height = 4)
  p
}

# ---- fda-probability ----

sb_predict_fda_native <- function(engine, x, type, dimension = ncol(engine$means)) {
  n <- nrow(x)
  # mda::polybasis drops the row dimension for singleton interaction queries.
  # Its polynomial evaluation is row-local; duplicate only to retain dimensions.
  singleton <- n == 1L && inherits(engine$fit, "polyreg") &&
    engine$fit$degree > 1L && !isTRUE(engine$fit$monomial) && ncol(x) > 1L
  if (singleton) x <- x[c(1L, 1L), , drop = FALSE]
  native <- get("predict.fda", asNamespace("mda"))
  result <- native(engine, newdata = x, type = type, dimension = dimension)
  result[seq_len(n), , drop = FALSE]
}

sb_predict_fda_probability <- function(engine, x) {
  distances <- sb_predict_fda_native(engine, x, type = "distances")
  # The native posterior divides exponentiated class weights. Their log odds
  # give the same binary posterior without a zero/zero underflow.
  log_odds <- -0.5 * (distances[, "1"] - distances[, "0"]) +
    log(engine$prior[["1"]]) - log(engine$prior[["0"]])
  stats::plogis(log_odds)
}
