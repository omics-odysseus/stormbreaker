# ---- classical-survival-fixed-adapters ----

sb_survivalsvm_linear_learner_ids <- function() {
  c(
    "survivalsvm",
    "survivalsvm_vanbelle1",
    "survivalsvm_vanbelle2",
    "survivalsvm_hybrid"
  )
}

sb_survivalsvm_additive_learner_ids <- function() {
  c(
    "survivalsvm_additive",
    "survivalsvm_vanbelle1_additive",
    "survivalsvm_vanbelle2_additive",
    "survivalsvm_hybrid_additive"
  )
}

sb_survivalsvm_rbf_learner_ids <- function() {
  c(
    "survivalsvm_rbf",
    "survivalsvm_vanbelle1_rbf",
    "survivalsvm_vanbelle2_rbf",
    "survivalsvm_hybrid_rbf"
  )
}

sb_survivalsvm_kernel_learner_ids <- function() {
  c(
    sb_survivalsvm_additive_learner_ids(),
    sb_survivalsvm_rbf_learner_ids()
  )
}

sb_survivalsvm_learner_ids <- function() {
  c(
    sb_survivalsvm_linear_learner_ids(),
    sb_survivalsvm_kernel_learner_ids()
  )
}

sb_survivalsvm_type <- function(learner) {
  type <- switch(learner,
    survivalsvm = "regression",
    survivalsvm_vanbelle1 = "vanbelle1",
    survivalsvm_vanbelle2 = "vanbelle2",
    survivalsvm_hybrid = "hybrid",
    survivalsvm_additive = "regression",
    survivalsvm_vanbelle1_additive = "vanbelle1",
    survivalsvm_vanbelle2_additive = "vanbelle2",
    survivalsvm_hybrid_additive = "hybrid",
    survivalsvm_rbf = "regression",
    survivalsvm_vanbelle1_rbf = "vanbelle1",
    survivalsvm_vanbelle2_rbf = "vanbelle2",
    survivalsvm_hybrid_rbf = "hybrid",
    NULL
  )
  if (is.null(type)) {
    sb_abort(
      paste0("Unknown fixed survival-SVM learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  type
}

sb_survivalsvm_kernel <- function(learner) {
  if (learner %in% sb_survivalsvm_linear_learner_ids()) {
    return("lin_kernel")
  }
  if (learner %in% sb_survivalsvm_additive_learner_ids()) {
    return("add_kernel")
  }
  if (learner %in% sb_survivalsvm_rbf_learner_ids()) {
    return("rbf_kernel")
  }
  sb_abort(
    paste0("Unknown fixed survival-SVM learner `", learner, "`."),
    "stormbreaker_unsupported_learner"
  )
}

sb_classical_survival_learner_ids <- function() {
  c(
    "coxboost",
    "plsrcox",
    "superpc",
    sb_survivalsvm_linear_learner_ids()
  )
}

sb_analysis_tuned_survival_learner_ids <- function() {
  c("coxboost_cv", "plsrcox_cv", "superpc_cv", "gbm_cox_cv")
}

sb_classical_survival_linear_learner_ids <- function() {
  c(
    sb_classical_survival_learner_ids(),
    setdiff(sb_analysis_tuned_survival_learner_ids(), "gbm_cox_cv")
  )
}

sb_linear_survival_learner_ids <- function() {
  c(
    "coxph",
    sb_coxphf_learner_ids(),
    "glmnet_cox",
    sb_classical_survival_linear_learner_ids(),
    sb_abess_survival_learner_ids(),
    sb_ncvreg_survival_learner_ids()
  )
}

sb_survival_learner_ids <- function() {
  c(
    sb_linear_survival_learner_ids(),
    sb_survivalsvm_kernel_learner_ids(),
    "catboost_cox",
    "catboost_aft_survival",
    "gbm_cox",
    "gbm_cox_cv",
    sb_mboost_cox_learner_ids(),
    sb_mboost_aft_learner_ids(),
    sb_mboost_ranking_learner_ids(),
    sb_ranger_survival_learner_ids(),
    sb_rlt_survival_learner_ids(),
    sb_bagged_knn_survival_learner_ids(), sb_akritas_learner_ids(),
    sb_ahaz_survival_learner_ids(),
    sb_bujar_linear_learner_ids(),
    sb_rule_aft_survival_learner_ids(),
    sb_xgboost_survival_learner_ids(),
    sb_flexsurv_parametric_learner_ids(),
    sb_flexsurvspline_learner_ids(),
    sb_mtlr_survival_learner_ids(),
    sb_pamr_survival_learner_ids(),
    sb_additional_survival_learner_ids(),
    sb_timereg_survival_learner_ids(),
    sb_mgcv_gam_survival_learner_ids()
  )
}

sb_coxboost_available <- function() {
  requireNamespace("CoxBoost", quietly = TRUE)
}

sb_plsrcox_available <- function() {
  requireNamespace("plsRcox", quietly = TRUE)
}

sb_superpc_available <- function() {
  requireNamespace("superpc", quietly = TRUE)
}

sb_survivalsvm_available <- function() {
  requireNamespace("survivalsvm", quietly = TRUE)
}

sb_validate_classical_survival_inputs <- function(x, time, event, learner) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time < 0)) {
    sb_abort(
      "`time` must contain non-negative values.",
      "stormbreaker_bad_adapter_time"
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  if (!any(event == 1L)) {
    sb_abort(
      paste0("A fixed ", learner, " analysis fold must contain an event."),
      "stormbreaker_classical_survival_no_event"
    )
  }
  list(x = x, time = as.numeric(time), event = as.integer(event))
}

sb_validate_positive_number <- function(value, argument, subclass) {
  valid <- is.numeric(value) && length(value) == 1L &&
    !is.na(value) && is.finite(value) && value > 0
  if (!valid) {
    sb_abort(
      paste0("`", argument, "` must be one finite positive number."),
      subclass
    )
  }
  as.numeric(value)
}

sb_classical_survival_stop_on_warning <- function(code, label, subclass) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      sb_abort(
        paste0(
          "The ", label, " raised a warning: ",
          conditionMessage(warning)
        ),
        subclass,
        parent = warning
      )
    }
  )
}

sb_coxboost_params <- function(stepno = 100L, penalty = NULL) {
  if (!is.null(penalty)) {
    penalty <- sb_validate_positive_number(
      penalty,
      "coxboost_penalty",
      "stormbreaker_bad_coxboost_penalty"
    )
  }
  list(
    stepno = sb_check_positive_integer(
      stepno,
      "coxboost_stepno",
      minimum = 1L
    ),
    penalty = penalty,
    penalty_rule = if (is.null(penalty)) {
      "9_times_analysis_event_count"
    } else {
      "fixed_numeric"
    },
    standardize = TRUE,
    criterion = "pscore"
  )
}

sb_coxboost_cv_params <- function(inner_folds = 10L, seed = 1L,
                                  maximum_steps = 500L,
                                  start_penalty = 500) {
  list(
    inner_folds = sb_check_positive_integer(
      inner_folds,
      "learner_tuning_folds",
      minimum = 3L
    ),
    seed = sb_check_positive_integer(seed, "learner_seed", minimum = 1L),
    maximum_steps = sb_check_positive_integer(
      maximum_steps,
      "maximum_steps",
      minimum = 1L
    ),
    start_penalty = sb_validate_positive_number(
      start_penalty,
      "start_penalty",
      "stormbreaker_bad_coxboost_penalty"
    ),
    standardize = TRUE,
    tuning_type = "verweij"
  )
}

sb_plsrcox_params <- function(components = 2L) {
  list(
    components = sb_check_positive_integer(
      components,
      "plsrcox_components",
      minimum = 1L
    ),
    scale_x = TRUE,
    scale_y = FALSE,
    sparse_stop = FALSE
  )
}

sb_plsrcox_cv_params <- function(inner_folds = 5L, seed = 1L,
                                 maximum_components = 10L) {
  list(
    inner_folds = sb_check_positive_integer(
      inner_folds,
      "learner_tuning_folds",
      minimum = 3L
    ),
    seed = sb_check_positive_integer(seed, "learner_seed", minimum = 1L),
    maximum_components = sb_check_positive_integer(
      maximum_components,
      "maximum_components",
      minimum = 1L
    ),
    selection_rule = "minimum_cv_error",
    scale_x = TRUE,
    scale_y = FALSE,
    sparse_stop = FALSE
  )
}

sb_superpc_params <- function(threshold = 1, components = 1L) {
  valid_threshold <- is.numeric(threshold) && length(threshold) == 1L &&
    !is.na(threshold) && is.finite(threshold) && threshold >= 0
  if (!valid_threshold) {
    sb_abort(
      "`superpc_threshold` must be one finite non-negative number.",
      "stormbreaker_bad_superpc_threshold"
    )
  }
  components <- sb_check_positive_integer(
    components,
    "superpc_components",
    minimum = 1L
  )
  if (components > 3L) {
    sb_abort(
      "`superpc_components` must be an integer from 1 through 3.",
      "stormbreaker_bad_superpc_components"
    )
  }
  list(
    threshold = as.numeric(threshold),
    components = components,
    s0_perc = 0.5
  )
}

sb_superpc_cv_params <- function(inner_folds = 10L, seed = 1L) {
  list(
    inner_folds = sb_check_positive_integer(
      inner_folds,
      "learner_tuning_folds",
      minimum = 3L
    ),
    seed = sb_check_positive_integer(seed, "learner_seed", minimum = 1L),
    threshold_count = 20L,
    tuning_components = 3L,
    prediction_components = 1L,
    minimum_features = 2L,
    s0_perc = 0.5
  )
}

sb_survivalsvm_params <- function(
  gamma_mu = 1,
  learner = "survivalsvm"
) {
  type <- sb_survivalsvm_type(learner)
  params <- list(
    gamma_mu = sb_validate_positive_number(
      gamma_mu,
      "survivalsvm_gamma_mu",
      "stormbreaker_bad_survivalsvm_gamma_mu"
    ),
    type = type,
    kernel = sb_survivalsvm_kernel(learner),
    opt_meth = "quadprog",
    standardize = TRUE,
    calibration_ties = "breslow"
  )
  if (identical(params$kernel, "rbf_kernel")) {
    params["kernel_pars"] <- list(NULL)
    params$kernel_parameter_rule <-
      "analysis_median_squared_distance"
  }
  if (!identical(type, "regression")) {
    params <- append(params, list(diff_meth = "makediff3"), after = 2L)
  }
  params
}

sb_new_classical_survival_model <- function(
  learner, features, coefficient, params, audit, engine,
  baseline_source
) {
  valid <- learner %in% sb_classical_survival_linear_learner_ids() &&
    is.character(features) && length(features) > 0L &&
    !anyNA(features) && !anyDuplicated(features) &&
    is.numeric(coefficient) && length(coefficient) == length(features) &&
    !anyNA(coefficient) && all(is.finite(coefficient))
  if (!valid) {
    sb_abort(
      paste0("Fixed ", learner, " did not yield portable coefficients."),
      "stormbreaker_invalid_classical_survival_fit"
    )
  }
  structure(
    list(
      learner = learner,
      features = features,
      coefficient = unname(coefficient),
      params = params,
      audit = audit,
      baseline_source = baseline_source,
      engine = engine
    ),
    class = c(
      paste0("stormbreaker_", learner),
      "stormbreaker_classical_survival_fixed"
    )
  )
}

sb_predict_classical_survival <- function(model, x) {
  if (!inherits(model, "stormbreaker_classical_survival_fixed")) {
    sb_abort(
      "The supplied model is not a fitted classical survival learner.",
      "stormbreaker_bad_classical_survival_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  risk <- unname(as.vector(as.matrix(x) %*% model$coefficient))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "The fitted survival learner produced invalid risk predictions.",
      "stormbreaker_invalid_classical_survival_prediction"
    )
  }
  risk
}

sb_extract_classical_survival_coefficients <- function(model) {
  if (!inherits(model, "stormbreaker_classical_survival_fixed")) {
    sb_abort(
      "The supplied model is not a fitted classical survival learner.",
      "stormbreaker_bad_classical_survival_fit"
    )
  }
  data.frame(
    feature = model$features,
    coefficient = model$coefficient,
    nonzero = model$coefficient != 0,
    direction = ifelse(
      model$coefficient > 0,
      "higher_hazard",
      ifelse(model$coefficient < 0, "lower_hazard", "zero")
    ),
    stringsAsFactors = FALSE
  )
}

sb_predict_classical_survival_survival <- function(
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
  analysis_risk <- sb_predict_classical_survival(model, inputs$x)
  assessment_risk <- sb_predict_classical_survival(
    model,
    x_assessment
  )
  cumulative_hazard <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = model$baseline_source,
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_fit_coxboost_cv_native <- function(x, time, event, params,
                                      label = "CoxBoost CV") {
  inputs <- sb_validate_classical_survival_inputs(x, time, event, label)
  if (!sb_coxboost_available()) {
    sb_abort(
      "Package `CoxBoost` is required for the CV-tuned CoxBoost learner.",
      "stormbreaker_missing_coxboost"
    )
  }
  fold_id <- sb_make_stratified_fold_id(
    inputs$event,
    params$inner_folds,
    params$seed,
    stratum_label = "event"
  )
  folds <- lapply(seq_len(params$inner_folds), function(fold) {
    which(fold_id == fold)
  })
  x_matrix <- as.matrix(inputs$x)
  penalty_search <- sb_classical_survival_stop_on_warning(
    withr::with_seed(
      params$seed,
      CoxBoost::optimCoxBoostPenalty(
        time = inputs$time,
        status = inputs$event,
        x = x_matrix,
        start.penalty = params$start_penalty,
        parallel = FALSE,
        trace = FALSE
      )
    ),
    paste0(label, " penalty search"),
    "stormbreaker_coxboost_cv_warning"
  )
  penalty <- as.numeric(penalty_search$penalty)
  if (length(penalty) != 1L || !is.finite(penalty) || penalty <= 0) {
    sb_abort(
      "CV-tuned CoxBoost returned an invalid penalty.",
      "stormbreaker_invalid_coxboost_cv"
    )
  }
  cross_validation <- sb_classical_survival_stop_on_warning(
    withr::with_seed(
      params$seed,
      CoxBoost::cv.CoxBoost(
        time = inputs$time,
        status = inputs$event,
        x = x_matrix,
        maxstepno = params$maximum_steps,
        K = params$inner_folds,
        type = params$tuning_type,
        penalty = penalty,
        folds = folds,
        parallel = FALSE,
        multicore = FALSE,
        trace = FALSE
      )
    ),
    paste0(label, " step search"),
    "stormbreaker_coxboost_cv_warning"
  )
  optimal_step <- as.integer(cross_validation$optimal.step)
  if (length(optimal_step) != 1L || is.na(optimal_step) ||
    optimal_step < 1L || optimal_step > params$maximum_steps) {
    sb_abort(
      "CV-tuned CoxBoost returned an invalid optimal step.",
      "stormbreaker_invalid_coxboost_cv"
    )
  }
  engine <- sb_classical_survival_stop_on_warning(
    CoxBoost::CoxBoost(
      time = inputs$time,
      status = inputs$event,
      x = x_matrix,
      stepno = optimal_step,
      penalty = penalty,
      standardize = params$standardize,
      trace = FALSE
    ),
    paste0(label, " final fit"),
    "stormbreaker_coxboost_cv_warning"
  )
  native_coefficient <- as.numeric(stats::coef(
    engine,
    at.step = optimal_step,
    scaled = FALSE
  ))
  scale <- as.numeric(engine$sdx)
  valid <- length(native_coefficient) == ncol(inputs$x) &&
    length(scale) == ncol(inputs$x) && !anyNA(native_coefficient) &&
    all(is.finite(native_coefficient)) && !anyNA(scale) &&
    all(is.finite(scale)) && all(scale > 0)
  if (!valid) {
    sb_abort(
      "CV-tuned CoxBoost returned invalid coefficients.",
      "stormbreaker_invalid_coxboost_cv"
    )
  }
  names(native_coefficient) <- names(inputs$x)
  list(
    inputs = inputs,
    engine = engine,
    native_coefficient = native_coefficient,
    coefficient = native_coefficient / scale,
    penalty = penalty,
    optimal_step = optimal_step,
    fold_id = fold_id,
    curve = data.frame(
      step = seq_along(cross_validation$mean.logplik) - 1L,
      mean_log_likelihood = as.numeric(cross_validation$mean.logplik),
      standard_error = as.numeric(cross_validation$se.logplik),
      stringsAsFactors = FALSE
    )
  )
}

sb_fit_coxboost_cv <- function(x, time, event, params) {
  native <- sb_fit_coxboost_cv_native(x, time, event, params)
  sb_new_classical_survival_model(
    learner = "coxboost_cv",
    features = names(native$inputs$x),
    coefficient = native$coefficient,
    params = params,
    audit = list(
      effective_penalty = native$penalty,
      optimal_step = native$optimal_step,
      fold_id = native$fold_id,
      tuning_curve = native$curve,
      standardize = params$standardize,
      tuning_type = params$tuning_type
    ),
    engine = native$engine,
    baseline_source = "coxboost_cv_breslow_log_hazard"
  )
}

sb_fit_coxboost <- function(x, time, event, params) {
  inputs <- sb_validate_classical_survival_inputs(
    x,
    time,
    event,
    "CoxBoost"
  )
  if (!sb_coxboost_available()) {
    sb_abort(
      "Package `CoxBoost` is required for the fixed CoxBoost learner.",
      "stormbreaker_missing_coxboost"
    )
  }
  effective_penalty <- if (is.null(params$penalty)) {
    9 * sum(inputs$event == 1L)
  } else {
    params$penalty
  }
  engine <- sb_classical_survival_stop_on_warning(
    CoxBoost::CoxBoost(
      time = inputs$time,
      status = inputs$event,
      x = as.matrix(inputs$x),
      standardize = params$standardize,
      stepno = params$stepno,
      penalty = effective_penalty,
      criterion = params$criterion,
      trace = FALSE
    ),
    "fixed CoxBoost learner",
    "stormbreaker_coxboost_fit_warning"
  )
  internal <- as.numeric(stats::coef(
    engine,
    at.step = params$stepno,
    scaled = FALSE
  ))
  scale <- as.numeric(engine$sdx)
  coefficient <- internal / scale
  sb_new_classical_survival_model(
    learner = "coxboost",
    features = names(inputs$x),
    coefficient = coefficient,
    params = params,
    audit = list(
      effective_penalty = as.numeric(effective_penalty),
      analysis_events = as.integer(sum(inputs$event == 1L)),
      fitted_steps = as.integer(engine$stepno),
      standardize = isTRUE(engine$standardize),
      criterion = params$criterion
    ),
    engine = engine,
    baseline_source = "coxboost_breslow_log_hazard"
  )
}

sb_fit_plsrcox_cv <- function(x, time, event, params) {
  inputs <- sb_validate_classical_survival_inputs(x, time, event, "plsRcox CV")
  if (!sb_plsrcox_available()) {
    sb_abort(
      "Package `plsRcox` is required for the CV-tuned plsRcox learner.",
      "stormbreaker_missing_plsrcox"
    )
  }
  single_predictor <- ncol(inputs$x) == 1L
  # Small late-component effects can trigger the native infinite-coefficient
  # warning when likelihood convergence precedes score convergence. Tighten
  # likelihood accuracy; keep the native Cholesky and infinite-coefficient
  # warning thresholds. toler.inf otherwise changes implicitly with eps.
  cox_control <- survival::coxph.control(
    eps = 1e-11,
    iter.max = 50L,
    toler.inf = survival::coxph.control()$toler.inf
  )
  search_warnings <- list()
  warning_indices <- list()
  abort_search_warning <- function(warning) {
    sb_abort(
      paste0("The plsRcox component search raised a warning: ",
             conditionMessage(warning)),
      "stormbreaker_plsrcox_cv_warning",
      parent = warning
    )
  }
  cross_validation <- if (single_predictor) {
    # There is exactly one feasible component, so the native component-count
    # search has no choice to make. cv.plsRcox drops one-column folds to a
    # vector; retain its native final learner at the unique component count.
    list(
      lambda.min5 = 1L,
      folds = list(),
      cv.error5 = NA_real_,
      cv.se5 = NA_real_
    )
  } else {
    withCallingHandlers(
      withr::with_seed(
        params$seed,
        plsRcox::cv.plsRcox(
          list(x = inputs$x, time = inputs$time, status = inputs$event),
          nfold = params$inner_folds,
          nt = params$maximum_components,
          plot.it = FALSE,
          folddetails = TRUE,
          verbose = FALSE,
          control = cox_control
        )
      ),
      warning = function(warning) {
        pattern <- paste0(
          "^Loglik converged before variable[[:space:]]+",
          "([0-9,[:space:]]+); coefficient may be infinite\\.[[:space:]]*$"
        )
        message <- conditionMessage(warning)
        matched <- regmatches(message, regexec(pattern, message))[[1L]]
        if (length(matched) != 2L) abort_search_warning(warning)
        indices <- as.integer(strsplit(
          gsub("[[:space:]]", "", matched[[2L]]), ",", fixed = TRUE
        )[[1L]])
        if (anyNA(indices) || any(indices < 1L |
                                 indices > params$maximum_components)) {
          abort_search_warning(warning)
        }
        search_warnings[[length(search_warnings) + 1L]] <<- warning
        warning_indices[[length(warning_indices) + 1L]] <<- indices
        invokeRestart("muffleWarning")
      }
    )
  }
  selected_components <- as.integer(cross_validation$lambda.min5)
  if (length(selected_components) != 1L || is.na(selected_components) ||
    selected_components < 1L ||
    selected_components > params$maximum_components) {
    sb_abort(
      "CV-tuned plsRcox returned an invalid component count.",
      "stormbreaker_invalid_plsrcox_cv"
    )
  }
  if (length(search_warnings) > 0L) {
    # Native component k fits Cox models with at most k columns. A warning
    # naming only later columns cannot come from the selected prefix.
    affected <- which(vapply(warning_indices, function(indices) {
      any(indices <= selected_components)
    }, logical(1)))
    if (length(affected) > 0L) {
      abort_search_warning(search_warnings[[affected[[1L]]]])
    }
    completed <- cross_validation$completed.cv5
    valid <- is.matrix(completed) &&
      nrow(completed) >= selected_components + 1L &&
      ncol(completed) == params$inner_folds &&
      !anyNA(completed[selected_components + 1L, ]) &&
      all(completed[selected_components + 1L, ])
    if (!valid) {
      sb_abort(
        "CV-tuned plsRcox did not complete every fold at the selected component count.",
        "stormbreaker_invalid_plsrcox_cv"
      )
    }
  }
  search_messages <- vapply(search_warnings, conditionMessage, character(1))
  engine <- sb_classical_survival_stop_on_warning(
    plsRcox::plsRcox(
      inputs$x,
      time = inputs$time,
      event = inputs$event,
      nt = selected_components,
      control = cox_control,
      scaleX = params$scale_x,
      scaleY = params$scale_y,
      sparseStop = params$sparse_stop,
      allres = FALSE,
      verbose = FALSE
    ),
    "CV-tuned plsRcox final fit",
    "stormbreaker_plsrcox_cv_warning"
  )
  coefficient_matrix <- engine$Coeffs
  coefficient <- if (is.matrix(coefficient_matrix) &&
    ncol(coefficient_matrix) >= 1L) {
    as.numeric(coefficient_matrix[, 1L])
  } else {
    numeric()
  }
  if (!identical(rownames(coefficient_matrix), names(inputs$x)) ||
    !identical(as.integer(engine$computed_nt), selected_components)) {
    coefficient <- numeric()
  }
  model <- sb_new_classical_survival_model(
    learner = "plsrcox_cv",
    features = names(inputs$x),
    coefficient = coefficient,
    params = params,
    audit = list(
      selected_components = selected_components,
      maximum_components = params$maximum_components,
      component_selection = if (single_predictor) {
        "single_predictor_has_one_feasible_component"
      } else {
        params$selection_rule
      },
      folds = cross_validation$folds,
      cv_error = as.numeric(cross_validation$cv.error5),
      cv_standard_error = as.numeric(cross_validation$cv.se5),
      component_search_warnings = search_messages,
      cox_control = cox_control
    ),
    engine = engine,
    baseline_source = "plsrcox_cv_breslow_log_hazard"
  )
  if (length(search_messages) > 0L) {
    attr(model, "stormbreaker_fit_warnings") <- search_messages
  }
  model
}

sb_fit_plsrcox <- function(x, time, event, params) {
  inputs <- sb_validate_classical_survival_inputs(
    x,
    time,
    event,
    "plsRcox"
  )
  if (!sb_plsrcox_available()) {
    sb_abort(
      "Package `plsRcox` is required for the fixed plsRcox learner.",
      "stormbreaker_missing_plsrcox"
    )
  }
  if (params$components > ncol(inputs$x)) {
    sb_abort(
      "`plsrcox_components` cannot exceed the fitted predictor count.",
      "stormbreaker_plsrcox_too_many_components"
    )
  }
  engine <- sb_classical_survival_stop_on_warning(
    plsRcox::plsRcox(
      inputs$x,
      time = inputs$time,
      event = inputs$event,
      nt = params$components,
      scaleX = params$scale_x,
      scaleY = params$scale_y,
      sparseStop = params$sparse_stop,
      allres = FALSE,
      verbose = FALSE
    ),
    "fixed plsRcox learner",
    "stormbreaker_plsrcox_fit_warning"
  )
  coefficient_matrix <- engine$Coeffs
  coefficient <- if (is.matrix(coefficient_matrix) &&
    ncol(coefficient_matrix) >= 1L) {
    as.numeric(coefficient_matrix[, 1L])
  } else {
    numeric()
  }
  coefficient_names <- rownames(coefficient_matrix)
  if (!is.null(coefficient_names) &&
    !identical(coefficient_names, names(inputs$x))) {
    coefficient <- numeric()
  }
  if (!identical(as.integer(engine$computed_nt), params$components)) {
    coefficient <- numeric()
  }
  sb_new_classical_survival_model(
    learner = "plsrcox",
    features = names(inputs$x),
    coefficient = coefficient,
    params = params,
    audit = list(
      computed_components = as.integer(engine$computed_nt),
      scale_x = params$scale_x,
      scale_y = params$scale_y,
      sparse_stop = params$sparse_stop
    ),
    engine = engine,
    baseline_source = "plsrcox_breslow_log_hazard"
  )
}

sb_superpc_portable_coefficient <- function(predicted, features,
                                            components) {
  selected <- as.logical(predicted$which.features)
  u <- predicted$u
  if (is.null(dim(u))) {
    u <- matrix(u, ncol = 1L)
  }
  d <- as.numeric(predicted$d)
  outcome_coefficient <- as.numeric(predicted$coef)
  valid <- length(selected) == length(features) &&
    nrow(u) == sum(selected) && ncol(u) == components &&
    length(d) == components && length(outcome_coefficient) == components
  coefficient <- numeric()
  if (valid) {
    scaling <- colSums(sweep(abs(u), 2L, d, "/"))
    selected_coefficient <- as.vector(
      u %*% (outcome_coefficient / (scaling * d))
    )
    coefficient <- numeric(length(features))
    coefficient[selected] <- selected_coefficient
  }
  list(coefficient = coefficient, selected = selected)
}

sb_fit_superpc_cv <- function(x, time, event, params) {
  inputs <- sb_validate_classical_survival_inputs(x, time, event, "SuperPC CV")
  if (!sb_superpc_available()) {
    sb_abort(
      "Package `superpc` is required for the CV-tuned SuperPC learner.",
      "stormbreaker_missing_superpc"
    )
  }
  if (ncol(inputs$x) < params$minimum_features) {
    sb_abort(
      paste0(
        "SuperPC CV threshold tuning requires at least ",
        params$minimum_features, " predictors but received ", ncol(inputs$x),
        "."
      ),
      "stormbreaker_superpc_cv_too_few_features"
    )
  }
  superpc_data <- list(
    x = t(as.matrix(inputs$x)),
    y = inputs$time,
    censoring.status = inputs$event,
    featurenames = names(inputs$x),
    data.names = as.character(seq_len(nrow(inputs$x)))
  )
  trained <- sb_classical_survival_stop_on_warning(
    withr::with_seed(
      params$seed,
      superpc::superpc.train(
        data = superpc_data,
        type = "survival",
        s0.perc = params$s0_perc
      )
    ),
    "SuperPC training call",
    "stormbreaker_superpc_cv_warning"
  )
  cross_validation <- sb_classical_survival_stop_on_warning(
    withr::with_seed(
      params$seed,
      superpc::superpc.cv(
        fit = trained,
        data = superpc_data,
        n.threshold = params$threshold_count,
        n.fold = params$inner_folds,
        n.components = params$tuning_components,
        min.features = params$minimum_features,
        max.features = nrow(superpc_data$x),
        compute.fullcv = TRUE,
        compute.preval = TRUE
      )
    ),
    "SuperPC threshold search",
    "stormbreaker_superpc_cv_warning"
  )
  best_index <- which.max(cross_validation$scor[1L, ])
  threshold <- as.numeric(cross_validation$thresholds[[best_index]])
  if (length(threshold) != 1L || !is.finite(threshold) || threshold < 0) {
    sb_abort(
      "CV-tuned SuperPC returned an invalid threshold.",
      "stormbreaker_invalid_superpc_cv"
    )
  }
  predicted <- sb_classical_survival_stop_on_warning(
    superpc::superpc.predict(
      object = trained,
      data = superpc_data,
      newdata = superpc_data,
      threshold = threshold,
      n.components = params$prediction_components,
      prediction.type = "continuous"
    ),
    "CV-tuned SuperPC final prediction call",
    "stormbreaker_superpc_cv_warning"
  )
  portable <- sb_superpc_portable_coefficient(
    predicted,
    names(inputs$x),
    params$prediction_components
  )
  sb_new_classical_survival_model(
    learner = "superpc_cv",
    features = names(inputs$x),
    coefficient = portable$coefficient,
    params = params,
    audit = list(
      selected_features = names(inputs$x)[portable$selected],
      selected_threshold = threshold,
      selected_threshold_index = as.integer(best_index),
      threshold_count = params$threshold_count,
      tuning_components = params$tuning_components,
      prediction_components = params$prediction_components,
      requested_folds = params$inner_folds,
      effective_folds = as.integer(cross_validation$n.fold),
      cv_scores = cross_validation$scor
    ),
    engine = list(
      trained = trained,
      cross_validation = cross_validation,
      prediction = predicted
    ),
    baseline_source = "superpc_cv_breslow_log_hazard"
  )
}

sb_superpc_single_feature_prediction <- function(trained, data, selected) {
  # The sole principal component is the centered selected predictor, up to
  # sign. Preserve matrix shape and use SuperPC's own outcome fit unchanged.
  x <- data$x[selected, , drop = FALSE]
  centered <- sweep(x, 1L, rowMeans(x), "-")
  decomposition <- svd(centered, nu = 1L, nv = 0L)
  u <- decomposition$u[, 1L, drop = FALSE]
  d <- decomposition$d[1L]
  scaling <- sum(abs(u) / d)
  score <- (t(centered) %*% u) / (scaling * d)
  coefficient <- superpc::superpc.fit.to.outcome(
    trained, data, score, print = FALSE
  )$results$coef
  out <- scale(score, center = FALSE, scale = sign(coefficient))
  list(
    u = u, d = d, coef = coefficient, which.features = selected,
    v.pred = out,
    v.pred.1df = apply(
      scale(out, center = FALSE, scale = 1 / abs(coefficient)), 1L, sum
    ),
    n.components = 1L, prediction.type = "continuous"
  )
}

sb_fit_superpc <- function(x, time, event, params) {
  inputs <- sb_validate_classical_survival_inputs(
    x,
    time,
    event,
    "SuperPC"
  )
  if (!sb_superpc_available()) {
    sb_abort(
      "Package `superpc` is required for the fixed SuperPC learner.",
      "stormbreaker_missing_superpc"
    )
  }
  if (ncol(inputs$x) < 2L) {
    sb_abort(
      paste0(
        "Fixed SuperPC training requires at least 2 predictors because ",
        "the native survival-score path does not preserve a one-row ",
        "feature matrix."
      ),
      "stormbreaker_superpc_too_few_features"
    )
  }
  superpc_data <- list(
    x = t(as.matrix(inputs$x)),
    y = inputs$time,
    censoring.status = inputs$event,
    featurenames = names(inputs$x),
    data.names = as.character(seq_len(nrow(inputs$x)))
  )
  trained <- sb_classical_survival_stop_on_warning(
    superpc::superpc.train(
      data = superpc_data,
      type = "survival",
      s0.perc = params$s0_perc
    ),
    "fixed SuperPC training call",
    "stormbreaker_superpc_train_warning"
  )
  predicted <- tryCatch(
    sb_classical_survival_stop_on_warning(
      {
        selected <- abs(trained$feature.scores) >= params$threshold
        if (sum(selected) == 1L && params$components == 1L) {
          sb_superpc_single_feature_prediction(trained, superpc_data, selected)
        } else {
          superpc::superpc.predict(
            object = trained,
            data = superpc_data,
            newdata = superpc_data,
            threshold = params$threshold,
            n.components = params$components,
            prediction.type = "continuous"
          )
        }
      },
      "fixed SuperPC prediction call",
      "stormbreaker_superpc_predict_warning"
    ),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) stop(error)
      selected_count <- sum(abs(trained$feature.scores) >= params$threshold)
      sb_abort(
        paste0(
          "Native SuperPC prediction failed after selecting ", selected_count,
          " of ", ncol(inputs$x), " predictors (threshold=", params$threshold,
          ", components=", params$components, "): ", conditionMessage(error)
        ),
        "stormbreaker_superpc_predict_error",
        parent = error
      )
    }
  )
  portable <- sb_superpc_portable_coefficient(
    predicted,
    names(inputs$x),
    params$components
  )
  sb_new_classical_survival_model(
    learner = "superpc",
    features = names(inputs$x),
    coefficient = portable$coefficient,
    params = params,
    audit = list(
      selected_features = names(inputs$x)[portable$selected],
      threshold = params$threshold,
      components = params$components,
      s0_perc = params$s0_perc
    ),
    engine = list(trained = trained, prediction = predicted),
    baseline_source = "superpc_breslow_log_hazard"
  )
}

sb_fit_survivalsvm_calibration <- function(time, event, raw_risk, ties) {
  calibration_data <- data.frame(
    time = time,
    event = event,
    raw_risk = raw_risk
  )
  sb_classical_survival_stop_on_warning(
    survival::coxph(
      survival::Surv(time, event) ~ raw_risk,
      data = calibration_data,
      ties = ties,
      singular.ok = FALSE
    ),
    "survival-SVM Cox calibration",
    "stormbreaker_survivalsvm_calibration_warning"
  )
}

sb_survivalsvm_native_components <- function(engine, type, feature_count) {
  fit <- engine$model.fit
  basis <- matrix(numeric(), nrow = 0L, ncol = 0L)
  row_weight <- numeric()
  intercept <- numeric()
  valid <- FALSE
  if (identical(type, "regression")) {
    basis <- as.matrix(fit$SV)
    row_weight <- as.numeric(fit$Beta)
    intercept <- as.numeric(fit$b0)
    valid <- ncol(basis) == feature_count &&
      nrow(basis) == length(row_weight) && length(intercept) == 1L &&
      is.finite(intercept) && !anyNA(basis) && all(is.finite(basis)) &&
      !anyNA(row_weight) && all(is.finite(row_weight))
  } else {
    basis <- as.matrix(fit$Xtrain)
    differences <- as.matrix(fit$DifMat$Mat)
    alpha <- as.numeric(fit$Alpha)
    valid <- ncol(basis) == feature_count &&
      nrow(basis) == ncol(differences) &&
      nrow(differences) == length(alpha) &&
      !anyNA(basis) && all(is.finite(basis)) &&
      !anyNA(differences) && all(is.finite(differences)) &&
      !anyNA(alpha) && all(is.finite(alpha))
    if (valid) {
      row_weight <- as.vector(crossprod(differences, alpha))
      intercept <- 0
      if (identical(type, "hybrid")) {
        beta <- as.numeric(fit$Beta)
        beta_star <- as.numeric(fit$Betarstar)
        delta <- as.numeric(fit$Delta)
        intercept <- as.numeric(fit$b0)
        valid <- length(beta) == nrow(basis) &&
          length(beta_star) == nrow(basis) &&
          length(delta) == nrow(basis) && length(intercept) == 1L &&
          is.finite(intercept) && !anyNA(beta) && all(is.finite(beta)) &&
          !anyNA(beta_star) && all(is.finite(beta_star)) &&
          !anyNA(delta) && all(is.finite(delta))
        if (valid) {
          row_weight <- row_weight + beta - delta * beta_star
        }
      }
    }
  }
  valid <- valid && nrow(basis) > 0L && length(row_weight) == nrow(basis) &&
    !anyNA(row_weight) && all(is.finite(row_weight))
  if (!valid) {
    sb_abort(
      "Fixed survival-SVM returned invalid native prediction components.",
      "stormbreaker_invalid_survivalsvm_fit"
    )
  }
  list(basis = basis, row_weight = row_weight, intercept = intercept)
}

sb_survivalsvm_linear_components <- function(engine, type, feature_count) {
  native <- sb_survivalsvm_native_components(engine, type, feature_count)
  weight <- as.vector(crossprod(native$basis, native$row_weight))
  if (length(weight) != feature_count || anyNA(weight) ||
    any(!is.finite(weight))) {
    sb_abort(
      "Fixed survival-SVM returned an invalid linear model.",
      "stormbreaker_invalid_survivalsvm_fit"
    )
  }
  list(weight = weight, intercept = native$intercept)
}

sb_survivalsvm_additive_kernel_matrix <- function(basis, x) {
  basis <- as.matrix(basis)
  x <- as.matrix(x)
  ranges <- apply(basis, 2L, max) - apply(basis, 2L, min)
  if (any(ranges == 0)) {
    sb_abort(
      "The fitted additive survival-SVM basis contains a constant feature.",
      "stormbreaker_invalid_survivalsvm_fit"
    )
  }
  transposed_x <- t(x)
  rows <- lapply(seq_len(nrow(basis)), function(index) {
    contribution <- t(
      (ranges - abs(transposed_x - basis[index, ])) / ranges
    )
    rowSums(contribution)
  })
  do.call(rbind, rows)
}

sb_survivalsvm_rbf_kernel_matrix <- function(basis, x, kernel_pars) {
  basis <- as.matrix(basis)
  x <- as.matrix(x)
  valid <- is.numeric(kernel_pars) && length(kernel_pars) > 0L &&
    !anyNA(kernel_pars) && all(is.finite(kernel_pars)) &&
    all(kernel_pars > 0)
  if (!valid) {
    sb_abort(
      "The fitted RBF survival-SVM has invalid kernel parameters.",
      "stormbreaker_invalid_survivalsvm_fit"
    )
  }
  basis_norm <- crossprod(
    t(rowSums(basis^2)),
    t(rep(1, nrow(x)))
  )
  x_norm <- crossprod(
    t(rowSums(x^2)),
    t(rep(1, nrow(basis)))
  )
  squared_distance <- basis_norm + t(x_norm) -
    2 * tcrossprod(basis, x)
  exp(-squared_distance / (2 * kernel_pars[[1L]]))
}

sb_survivalsvm_kernel_matrix <- function(basis, x, kernel, kernel_pars) {
  switch(kernel,
    add_kernel = sb_survivalsvm_additive_kernel_matrix(basis, x),
    rbf_kernel = sb_survivalsvm_rbf_kernel_matrix(
      basis,
      x,
      kernel_pars
    ),
    sb_abort(
      paste0("Unsupported portable survival-SVM kernel `", kernel, "`."),
      "stormbreaker_invalid_survivalsvm_fit"
    )
  )
}

sb_survivalsvm_kernel_reference <- function(kernel) {
  switch(kernel,
    add_kernel = "survivalsvm_additive_kernel_native",
    rbf_kernel = "survivalsvm_rbf_kernel_native",
    sb_abort(
      paste0("Unsupported portable survival-SVM kernel `", kernel, "`."),
      "stormbreaker_invalid_survivalsvm_fit"
    )
  )
}

sb_survivalsvm_kernel_model_class <- function(kernel) {
  switch(kernel,
    add_kernel = "stormbreaker_survivalsvm_additive",
    rbf_kernel = "stormbreaker_survivalsvm_rbf",
    sb_abort(
      paste0("Unsupported portable survival-SVM kernel `", kernel, "`."),
      "stormbreaker_invalid_survivalsvm_fit"
    )
  )
}

sb_survivalsvm_kernel_state_class <- function(kernel) {
  switch(kernel,
    add_kernel = "stormbreaker_survivalsvm_additive_learner_state",
    rbf_kernel = "stormbreaker_survivalsvm_rbf_learner_state",
    sb_abort(
      paste0("Unsupported portable survival-SVM kernel `", kernel, "`."),
      "stormbreaker_invalid_survivalsvm_fit"
    )
  )
}

sb_new_survivalsvm_kernel_model <- function(
  learner, features, basis, row_weight, intercept,
  calibration_coefficient, kernel_pars, params, audit, engine,
  baseline_source
) {
  kernel <- sb_survivalsvm_kernel(learner)
  valid_kernel_pars <- if (identical(kernel, "rbf_kernel")) {
    is.numeric(kernel_pars) && length(kernel_pars) > 0L &&
      !anyNA(kernel_pars) && all(is.finite(kernel_pars)) &&
      all(kernel_pars > 0)
  } else {
    is.null(kernel_pars)
  }
  valid <- learner %in% sb_survivalsvm_kernel_learner_ids() &&
    identical(params$kernel, kernel) && valid_kernel_pars &&
    is.character(features) && length(features) > 0L &&
    !anyNA(features) && !anyDuplicated(features) &&
    is.matrix(basis) && is.numeric(basis) && ncol(basis) == length(features) &&
    nrow(basis) > 0L && !anyNA(basis) && all(is.finite(basis)) &&
    is.numeric(row_weight) && length(row_weight) == nrow(basis) &&
    !anyNA(row_weight) && all(is.finite(row_weight)) &&
    is.numeric(intercept) && length(intercept) == 1L &&
    !is.na(intercept) && is.finite(intercept) &&
    is.numeric(calibration_coefficient) &&
    length(calibration_coefficient) == 1L &&
    !is.na(calibration_coefficient) && is.finite(calibration_coefficient)
  if (!valid) {
    sb_abort(
      paste0("Fixed ", learner, " did not yield a portable kernel model."),
      "stormbreaker_invalid_survivalsvm_fit"
    )
  }
  structure(
    list(
      learner = learner,
      features = features,
      basis = basis,
      row_weight = unname(row_weight),
      intercept = unname(intercept),
      calibration_coefficient = unname(calibration_coefficient),
      kernel_pars = kernel_pars,
      params = params,
      audit = audit,
      baseline_source = baseline_source,
      engine = engine
    ),
    class = c(
      paste0("stormbreaker_", learner),
      sb_survivalsvm_kernel_model_class(kernel),
      "stormbreaker_survivalsvm_kernel_fixed"
    )
  )
}

sb_predict_survivalsvm_kernel <- function(model, x) {
  if (!inherits(model, "stormbreaker_survivalsvm_kernel_fixed")) {
    sb_abort(
      "The supplied model is not a fitted kernel survival-SVM learner.",
      "stormbreaker_bad_survivalsvm_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted survival-SVM schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  scaled_x <- sweep(
    sweep(as.matrix(x), 2L, model$audit$center, "-"),
    2L,
    model$audit$scale,
    "/"
  )
  kernel <- sb_survivalsvm_kernel_matrix(
    model$basis,
    scaled_x,
    model$params$kernel,
    model$kernel_pars
  )
  native_score <- as.vector(crossprod(model$row_weight, kernel)) +
    model$intercept
  risk <- -model$calibration_coefficient * native_score
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "The fitted kernel survival-SVM produced invalid risk predictions.",
      "stormbreaker_invalid_survivalsvm_prediction"
    )
  }
  unname(risk)
}

sb_predict_survivalsvm_kernel_survival <- function(
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
  analysis_risk <- sb_predict_survivalsvm_kernel(model, inputs$x)
  assessment_risk <- sb_predict_survivalsvm_kernel(
    model,
    x_assessment
  )
  # A common log-risk offset cancels in H0(t) * exp(risk). Store the
  # baseline at the training mean to avoid overflow from native intercepts.
  risk_reference <- mean(analysis_risk)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    analysis_risk - risk_reference,
    prediction_horizon
  )
  result <- sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = model$baseline_source,
    ties = "breslow",
    risk = assessment_risk - risk_reference
  )
  result$risk_reference <- risk_reference
  result
}

sb_extract_survivalsvm_kernel_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_survivalsvm_kernel_fixed")) {
    sb_bad_learner_state(
      "The fitted kernel survival-SVM learner state is invalid."
    )
  }
  kernel <- sb_survivalsvm_kernel(model$learner)
  basis <- model$basis
  rownames(basis) <- NULL
  baseline_hazard <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_survivalsvm_kernel_survival(
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
      cumulative_baseline_hazard = probability$cumulative_baseline_hazard,
      risk_reference = probability$risk_reference,
      baseline_survival = probability$baseline_survival,
      source = probability$source,
      ties = probability$ties
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      audit = model$audit,
      center = model$audit$center,
      scale = model$audit$scale,
      basis = basis,
      row_weight = model$row_weight,
      intercept = model$intercept,
      calibration_coefficient = model$calibration_coefficient,
      kernel_pars = model$kernel_pars,
      risk_reference = sb_survivalsvm_kernel_reference(kernel),
      risk_formula =
        "negative_native_kernel_score_times_analysis_cox_coefficient",
      training_data_retention =
        "required_scaled_kernel_basis_without_outcomes",
      baseline_hazard = baseline_hazard
    ),
    class = c(
      sb_survivalsvm_kernel_state_class(kernel),
      "stormbreaker_survivalsvm_kernel_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_survivalsvm_kernel_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_survivalsvm_kernel_learner_ids()) {
    sb_bad_learner_state("The saved kernel survival-SVM state is invalid.")
  }
  kernel_type <- sb_survivalsvm_kernel(learner)
  x <- sb_validate_learner_state_x(state, x, learner)
  basis <- state[["basis"]]
  row_weight <- state[["row_weight"]]
  intercept <- state[["intercept"]]
  calibration_coefficient <- state[["calibration_coefficient"]]
  kernel_pars <- state[["kernel_pars"]]
  center <- state[["center"]]
  scale <- state[["scale"]]
  valid <- is.matrix(basis) && is.numeric(basis) &&
    identical(dim(basis), c(length(row_weight), ncol(x))) &&
    nrow(basis) > 0L && !anyNA(basis) && all(is.finite(basis)) &&
    is.numeric(row_weight) && !anyNA(row_weight) &&
    all(is.finite(row_weight)) &&
    is.numeric(intercept) && length(intercept) == 1L &&
    !is.na(intercept) && is.finite(intercept) &&
    is.numeric(calibration_coefficient) &&
    length(calibration_coefficient) == 1L &&
    !is.na(calibration_coefficient) && is.finite(calibration_coefficient) &&
    is.numeric(center) && length(center) == ncol(x) &&
    !anyNA(center) && all(is.finite(center)) &&
    is.numeric(scale) && length(scale) == ncol(x) &&
    !anyNA(scale) && all(is.finite(scale)) && all(scale > 0)
  valid_kernel_pars <- if (identical(kernel_type, "rbf_kernel")) {
    is.numeric(kernel_pars) && length(kernel_pars) > 0L &&
      !anyNA(kernel_pars) && all(is.finite(kernel_pars)) &&
      all(kernel_pars > 0)
  } else {
    is.null(kernel_pars)
  }
  valid <- valid && valid_kernel_pars
  if (!valid) {
    sb_bad_learner_state("The saved kernel survival-SVM state is invalid.")
  }
  scaled_x <- sweep(
    sweep(as.matrix(x), 2L, center, "-"),
    2L,
    scale,
    "/"
  )
  kernel <- sb_survivalsvm_kernel_matrix(
    basis,
    scaled_x,
    kernel_type,
    kernel_pars
  )
  native_score <- as.vector(crossprod(row_weight, kernel)) + intercept
  risk <- unname(-calibration_coefficient * native_score)
  if (anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "The saved kernel survival-SVM state produced invalid risk scores.",
      "stormbreaker_invalid_survivalsvm_prediction"
    )
  }
  result <- list(risk = risk)
  baseline <- state[["baseline_hazard"]]
  if (!is.null(baseline)) {
    valid_baseline <- is.list(baseline) &&
      is.numeric(baseline[["prediction_horizon"]]) &&
      length(baseline[["prediction_horizon"]]) == 1L &&
      !is.na(baseline[["prediction_horizon"]]) &&
      is.finite(baseline[["prediction_horizon"]]) &&
      baseline[["prediction_horizon"]] > 0
    if (!valid_baseline) {
      sb_bad_learner_state(
        "The saved kernel survival-SVM baseline-hazard state is invalid."
      )
    }
    result$prediction_horizon <- baseline[["prediction_horizon"]]
    result$time_unit <- baseline[["time_unit"]]
    # Legacy states used the uncentered baseline (reference zero).
    risk_reference <- baseline[["risk_reference"]]
    if (is.null(risk_reference)) {
      risk_reference <- 0
    }
    if (!is.numeric(risk_reference) || length(risk_reference) != 1L ||
      is.na(risk_reference) || !is.finite(risk_reference)) {
      sb_bad_learner_state(
        "The saved kernel survival-SVM baseline risk reference is invalid."
      )
    }
    result$survival_probability <- sb_survival_probability_from_hazard(
      risk - risk_reference,
      baseline[["cumulative_baseline_hazard"]]
    )
  }
  result
}

sb_fit_survivalsvm <- function(
  x,
  time,
  event,
  params,
  learner = "survivalsvm"
) {
  inputs <- sb_validate_classical_survival_inputs(
    x,
    time,
    event,
    "survival-SVM"
  )
  if (!sb_survivalsvm_available()) {
    sb_abort(
      "Package `survivalsvm` is required for the fixed survival-SVM learner.",
      "stormbreaker_missing_survivalsvm"
    )
  }
  minimum_features <- if (
    identical(params$kernel, "add_kernel") ||
      (identical(params$kernel, "rbf_kernel") &&
        !identical(params$type, "regression"))
  ) {
    2L
  } else {
    1L
  }
  if (ncol(inputs$x) < minimum_features) {
    sb_abort(
      paste0(
        "Native ", params$type, " survival-SVM with ", params$kernel,
        " requires at least ", minimum_features, " predictors but received ",
        ncol(inputs$x), "."
      ),
      "stormbreaker_survivalsvm_too_few_features"
    )
  }
  center <- vapply(inputs$x, mean, numeric(1))
  scale <- vapply(inputs$x, stats::sd, numeric(1))
  scale[!is.finite(scale) | scale == 0] <- 1
  scaled_x <- as.data.frame(
    sweep(sweep(as.matrix(inputs$x), 2L, center, "-"), 2L, scale, "/"),
    optional = TRUE,
    stringsAsFactors = FALSE
  )
  names(scaled_x) <- names(inputs$x)
  time_name <- ".stormbreaker_time"
  while (time_name %in% names(scaled_x)) {
    time_name <- paste0(time_name, "_")
  }
  event_name <- ".stormbreaker_event"
  while (event_name %in% c(names(scaled_x), time_name)) {
    event_name <- paste0(event_name, "_")
  }
  analysis <- scaled_x
  analysis[[time_name]] <- inputs$time
  analysis[[event_name]] <- inputs$event
  response <- paste0(
    "survival::Surv(",
    sb_quote_name(time_name), ", ",
    sb_quote_name(event_name), ")"
  )
  rhs <- paste(vapply(names(inputs$x), sb_quote_name, ""), collapse = " + ")
  formula <- stats::as.formula(
    paste(response, "~", rhs),
    env = parent.frame()
  )
  if ("diff_meth" %in% names(params)) {
    # survivalsvm 0.0.6 resolves this exported function in the caller frame.
    makediff3 <- survivalsvm::makediff3
  }
  kernel_parameters <- params$kernel_pars
  if (identical(params$kernel_parameter_rule, "analysis_median_squared_distance")) {
    kernel_parameters <- stats::median(as.numeric(stats::dist(scaled_x))^2)
  }
  engine <- sb_classical_survival_stop_on_warning(
    survivalsvm::survivalsvm(
      formula = formula,
      data = analysis,
      type = params$type,
      diff.meth = if ("diff_meth" %in% names(params)) {
        params$diff_meth
      } else {
        NULL
      },
      gamma.mu = if (identical(params$type, "hybrid")) {
        rep(params$gamma_mu, 2L)
      } else {
        params$gamma_mu
      },
      opt.meth = params$opt_meth,
      kernel = params$kernel,
      kernel.pars = kernel_parameters
    ),
    "fixed survival-SVM learner",
    "stormbreaker_survivalsvm_fit_warning"
  )
  kernel_pars <- NULL
  if (params$kernel %in% c("add_kernel", "rbf_kernel")) {
    components <- sb_survivalsvm_native_components(
      engine,
      params$type,
      ncol(inputs$x)
    )
    if (identical(params$kernel, "rbf_kernel")) {
      kernel_pars <- as.numeric(engine$model.fit$Kernel$Kernpar)
    }
    kernel <- sb_survivalsvm_kernel_matrix(
      components$basis,
      scaled_x,
      params$kernel,
      kernel_pars
    )
    native_score <- as.vector(crossprod(components$row_weight, kernel)) +
      components$intercept
    raw_risk <- -native_score
  } else {
    components <- sb_survivalsvm_linear_components(
      engine,
      params$type,
      ncol(inputs$x)
    )
    scaled_weight <- components$weight
    raw_risk <- -(
      as.vector(as.matrix(scaled_x) %*% scaled_weight) + components$intercept
    )
  }
  calibration <- sb_fit_survivalsvm_calibration(
    inputs$time,
    inputs$event,
    raw_risk,
    params$calibration_ties
  )
  calibration_coefficient <- as.numeric(stats::coef(calibration))
  audit <- list(
    calibration_coefficient = calibration_coefficient,
    raw_score = switch(params$type,
      regression = "negative_predicted_survival_time",
      vanbelle1 = "negative_native_vanbelle1_ranking_score",
      vanbelle2 = "negative_native_vanbelle2_ranking_score",
      hybrid = "negative_native_hybrid_score"
    ),
    center = center,
    scale = scale,
    type = params$type,
    diff_meth = if ("diff_meth" %in% names(params)) {
      params$diff_meth
    } else {
      NULL
    },
    kernel = params$kernel,
    kernel_pars = kernel_pars,
    opt_meth = params$opt_meth
  )
  baseline_source <- paste0(
    learner,
    "_breslow_calibrated_log_hazard"
  )
  if (params$kernel %in% c("add_kernel", "rbf_kernel")) {
    return(sb_new_survivalsvm_kernel_model(
      learner = learner,
      features = names(inputs$x),
      basis = components$basis,
      row_weight = components$row_weight,
      intercept = components$intercept,
      calibration_coefficient = calibration_coefficient,
      kernel_pars = kernel_pars,
      params = params,
      audit = audit,
      engine = engine,
      baseline_source = baseline_source
    ))
  }
  coefficient <- -calibration_coefficient * scaled_weight / scale
  sb_new_classical_survival_model(
    learner = learner,
    features = names(inputs$x),
    coefficient = coefficient,
    params = params,
    audit = audit,
    engine = engine,
    baseline_source = baseline_source
  )
}

sb_classical_survival_adapter <- function(id, params, fit_function) {
  list(
    id = id,
    params = params,
    fit = function(x, time, event) {
      fit_function(x, time, event, params)
    },
    predict = sb_predict_classical_survival,
    predict_survival = sb_predict_classical_survival_survival,
    extract_coefficients = sb_extract_classical_survival_coefficients,
    extract_state = function(model, context = list()) {
      sb_extract_cox_learner_state(
        model,
        learner = id,
        coefficients = sb_extract_classical_survival_coefficients(model),
        predict_survival = sb_predict_classical_survival_survival,
        context = context,
        params = model$params,
        audit = model$audit
      )
    },
    predict_state = sb_predict_cox_learner_state
  )
}

sb_coxboost_learner_adapter <- function(stepno = 100L,
                                              penalty = NULL) {
  sb_classical_survival_adapter(
    "coxboost",
    sb_coxboost_params(stepno, penalty),
    sb_fit_coxboost
  )
}

sb_coxboost_cv_learner_adapter <- function(inner_folds = 10L,
                                           seed = 1L,
                                           maximum_steps = 500L,
                                           start_penalty = 500) {
  sb_classical_survival_adapter(
    "coxboost_cv",
    sb_coxboost_cv_params(
      inner_folds = inner_folds,
      seed = seed,
      maximum_steps = maximum_steps,
      start_penalty = start_penalty
    ),
    sb_fit_coxboost_cv
  )
}

sb_plsrcox_learner_adapter <- function(components = 2L) {
  sb_classical_survival_adapter(
    "plsrcox",
    sb_plsrcox_params(components),
    sb_fit_plsrcox
  )
}

sb_plsrcox_cv_learner_adapter <- function(inner_folds = 5L,
                                          seed = 1L,
                                          maximum_components = 10L) {
  sb_classical_survival_adapter(
    "plsrcox_cv",
    sb_plsrcox_cv_params(
      inner_folds = inner_folds,
      seed = seed,
      maximum_components = maximum_components
    ),
    sb_fit_plsrcox_cv
  )
}

sb_superpc_learner_adapter <- function(threshold = 1,
                                             components = 1L) {
  sb_classical_survival_adapter(
    "superpc",
    sb_superpc_params(threshold, components),
    sb_fit_superpc
  )
}

sb_superpc_cv_learner_adapter <- function(inner_folds = 10L,
                                          seed = 1L) {
  sb_classical_survival_adapter(
    "superpc_cv",
    sb_superpc_cv_params(inner_folds = inner_folds, seed = seed),
    sb_fit_superpc_cv
  )
}

sb_survivalsvm_learner_adapter <- function(
  gamma_mu = 1,
  learner = "survivalsvm",
  params = sb_survivalsvm_params(gamma_mu, learner)
) {
  if (learner %in% sb_survivalsvm_kernel_learner_ids()) {
    return(list(
      id = learner,
      params = params,
      fit = function(x, time, event) {
        sb_fit_survivalsvm(x, time, event, params, learner)
      },
      predict = sb_predict_survivalsvm_kernel,
      predict_survival = sb_predict_survivalsvm_kernel_survival,
      extract_state = sb_extract_survivalsvm_kernel_state,
      predict_state = sb_predict_survivalsvm_kernel_state
    ))
  }
  sb_classical_survival_adapter(
    learner,
    params,
    function(x, time, event, params) {
      sb_fit_survivalsvm(x, time, event, params, learner)
    }
  )
}

# ---- additional-survival-fixed-adapters ----

sb_survreg_learner_ids <- function() {
  c(
    "survreg_weibull",
    "survreg_exponential",
    "survreg_lognormal",
    "survreg_loglogistic"
  )
}

sb_additional_survival_learner_ids <- function() {
  c(
    sb_survreg_learner_ids(),
    "rpart_survival",
    "rsf_survival",
    "ctree_survival",
    "cforest_survival",
    "bart_survival",
    "aorsf_survival"
  )
}

sb_additional_survival_params <- function(learner) {
  switch(learner,
    survreg_weibull = list(
      distribution = "weibull",
      scale = "estimated"
    ),
    survreg_exponential = list(
      distribution = "exponential",
      scale = "fixed_1"
    ),
    survreg_lognormal = list(
      distribution = "lognormal",
      scale = "estimated"
    ),
    survreg_loglogistic = list(
      distribution = "loglogistic",
      scale = "estimated"
    ),
    rpart_survival = list(
      method = "exp",
      cp = 0.01,
      minsplit = 20L,
      maxdepth = 30L,
      xval = 0L
    ),
    rsf_survival = list(
      num_trees = 1000L,
      mtry_rule = "randomForestSRC_default",
      min_node_size = 5L,
      splitrule = "logrank",
      seed = 1L
    ),
    ctree_survival = list(
      mincriterion = 0.95,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = Inf
    ),
    cforest_survival = list(
      num_trees = 500L,
      mtry_rule = "ceiling_sqrt_p",
      replace = TRUE,
      sample_fraction = 1,
      mincriterion = 0,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = Inf,
      seed = 1L
    ),
    bart_survival = list(
      type = "pbart",
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
    aorsf_survival = list(
      num_trees = 500L,
      num_split = 5L,
      num_retry = 3L,
      mtry_rule = "ceiling_sqrt_p",
      sample_with_replacement = TRUE,
      sample_fraction = 0.632,
      control_method = "glm",
      control_scale_x = TRUE,
      control_ties = "efron",
      control_max_iter = 1L,
      seed = 1L,
      threads = 1L,
      importance = "none",
      attach_data = FALSE
    ),
    sb_abort(
      paste0("Unknown fixed survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_additional_survival_dependency <- function(learner) {
  packages <- switch(learner,
    survreg_weibull = "survival",
    survreg_exponential = "survival",
    survreg_lognormal = "survival",
    survreg_loglogistic = "survival",
    rpart_survival = "rpart",
    rsf_survival = "randomForestSRC",
    ctree_survival = c("partykit", "coin"),
    cforest_survival = c("partykit", "coin"),
    bart_survival = "BART",
    aorsf_survival = "aorsf"
  )
  list(
    package = packages[[1L]],
    packages = packages,
    available = all(vapply(
      packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    ))
  )
}

sb_require_additional_survival_dependency <- function(learner) {
  dependency <- sb_additional_survival_dependency(learner)
  if (!dependency$available) {
    missing <- dependency$packages[!vapply(
      dependency$packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )]
    sb_abort(
      paste0(
        "Package", if (length(missing) > 1L) "s" else "", " `",
        paste(missing, collapse = "`, `"), "` ",
        if (length(missing) > 1L) "are" else "is",
        " required for learner `", learner, "`."
      ),
      "stormbreaker_missing_additional_survival_dependency"
    )
  }
  invisible(dependency)
}

sb_validate_additional_survival_inputs <- function(x, time, event, learner) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time <= 0)) {
    sb_abort(
      paste0("Learner `", learner, "` requires strictly positive times."),
      "stormbreaker_additional_survival_nonpositive_time"
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  if (!any(event == 1L)) {
    sb_abort(
      paste0("Learner `", learner, "` requires at least one event."),
      "stormbreaker_additional_survival_no_event"
    )
  }
  list(x = x, time = time, event = as.integer(event))
}

sb_additional_survival_formula_data <- function(x, time, event,
                                                random_forest_src = FALSE) {
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
  rhs <- paste(vapply(names(x), sb_quote_name, ""), collapse = " + ")
  response <- if (random_forest_src) {
    paste0(
      "Surv(", sb_quote_name(time_name), ", ",
      sb_quote_name(event_name), ")"
    )
  } else {
    paste0(
      "survival::Surv(", sb_quote_name(time_name), ", ",
      sb_quote_name(event_name), ")"
    )
  }
  formula <- stats::as.formula(
    paste(response, "~", rhs),
    env = if (random_forest_src) asNamespace("survival") else baseenv()
  )
  list(data = analysis, formula = formula)
}

sb_new_additional_survival_model <- function(
  learner, features, params, engine = NULL, prediction = NULL,
  effective_mtry = NULL
) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      engine = engine,
      prediction = prediction,
      effective_mtry = effective_mtry
    ),
    class = "stormbreaker_additional_survival_fixed"
  )
}

sb_bart_survival_prediction_engine <- function(engine) {
  keep <- c(
    "treedraws", "rm.const", "binaryOffset", "offset", "type", "times", "K"
  )
  engine[setdiff(names(engine), keep)] <- NULL
  engine
}

sb_fit_bart_survival <- function(x, time, event, params) {
  engine <- NULL
  invisible(utils::capture.output(
    engine <- withr::with_seed(
      params$seed,
      withCallingHandlers(
        BART::surv.bart(
          x.train = as.matrix(x),
          times = time,
          delta = event,
          K = params$time_grid_size,
          sparse = params$sparse,
          type = params$type,
          k = params$prior_k,
          power = params$tree_power,
          base = params$tree_base,
          ntree = params$num_trees,
          numcut = params$num_cutpoints,
          ndpost = params$posterior_draws,
          nskip = params$burn_in,
          keepevery = params$keep_every,
          printevery = .Machine$integer.max,
          seed = params$seed,
          mc.cores = params$cores
        ),
        warning = function(warning) {
          sb_abort(
            paste0(
              "Fixed survival BART raised a fitting warning: ",
              conditionMessage(warning)
            ),
            "stormbreaker_bart_survival_fit_warning",
            parent = warning
          )
        }
      )
    ),
    type = "output"
  ))
  valid <- inherits(engine, "survbart") && !is.null(engine$treedraws) &&
    is.numeric(engine$times) && length(engine$times) > 0L &&
    !anyNA(engine$times) && all(is.finite(engine$times)) &&
    all(engine$times > 0) && identical(as.integer(engine$K),
      as.integer(length(engine$times)))
  if (!valid) {
    sb_abort(
      "Fixed survival BART did not return posterior trees and a time grid.",
      "stormbreaker_incomplete_bart_survival_fit"
    )
  }
  sb_bart_survival_prediction_engine(engine)
}

sb_fit_survreg <- function(formula_data, learner, params) {
  arguments <- list(
    formula = formula_data$formula,
    data = formula_data$data,
    dist = params$distribution,
    model = FALSE,
    x = FALSE,
    y = FALSE
  )
  if (!identical(params$distribution, "exponential")) {
    arguments$scale <- 0
  }
  withCallingHandlers(
    do.call(survival::survreg, arguments),
    warning = function(warning) {
      if (identical(
        conditionMessage(warning),
        "Ran out of iterations and did not converge"
      )) {
        sb_abort(
          paste0(
            "Fixed ", params$distribution, " SurvReg did not converge: ",
            conditionMessage(warning)
          ),
          "stormbreaker_survreg_nonconvergence",
          parent = warning
        )
      }
    }
  )
}

# The native default callback closes over rfsrc's full training frame. Keep the
# same default SWOR sample-size rule in the package namespace instead.
sb_rsf_native_sample_size <- function(x) x * 0.632

sb_fit_additional_survival <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_additional_survival_dependency(learner)
  formula_x <- inputs$x
  if (identical(learner, "aorsf_survival")) {
    names(formula_x) <- make.names(names(formula_x), unique = TRUE)
  }
  formula_data <- sb_additional_survival_formula_data(
    formula_x,
    inputs$time,
    inputs$event,
    random_forest_src = identical(learner, "rsf_survival")
  )

  engine <- if (learner %in% sb_survreg_learner_ids()) {
    sb_fit_survreg(formula_data, learner, params)
  } else {
    switch(learner,
    rpart_survival = rpart::rpart(
      formula_data$formula,
      data = formula_data$data,
      method = params$method,
      model = FALSE,
      x = FALSE,
      y = FALSE,
      control = rpart::rpart.control(
        cp = params$cp,
        minsplit = params$minsplit,
        maxdepth = params$maxdepth,
        xval = params$xval
      )
    ),
    rsf_survival = withr::with_options(
      list(rf.cores = 1L, mc.cores = 1L),
      randomForestSRC::rfsrc(
        formula_data$formula,
        data = formula_data$data,
        ntree = params$num_trees,
        nodesize = params$min_node_size,
        sampsize = sb_rsf_native_sample_size,
        splitrule = params$splitrule,
        importance = FALSE,
        proximity = FALSE,
        forest = TRUE,
        seed = params$seed
      )
    ),
    ctree_survival = partykit::ctree(
      formula_data$formula,
      data = formula_data$data,
      control = partykit::ctree_control(
        mincriterion = params$mincriterion,
        minsplit = params$minsplit,
        minbucket = params$minbucket,
        maxdepth = params$maxdepth
      )
    ),
    cforest_survival = withr::with_seed(
      params$seed,
      partykit::cforest(
        formula_data$formula,
        data = formula_data$data,
        ntree = params$num_trees,
        mtry = max(1L, ceiling(sqrt(ncol(inputs$x)))),
        perturb = list(
          replace = params$replace,
          fraction = params$sample_fraction
        ),
        control = partykit::ctree_control(
          teststat = "quad",
          testtype = "Univ",
          mincriterion = params$mincriterion,
          saveinfo = FALSE,
          minsplit = params$minsplit,
          minbucket = params$minbucket,
          maxdepth = params$maxdepth
        )
      )
    ),
    bart_survival = sb_fit_bart_survival(
      inputs$x,
      inputs$time,
      inputs$event,
      params
    ),
    aorsf_survival = {
      tree_seeds <- withr::with_seed(
        params$seed,
        sample.int(.Machine$integer.max, params$num_trees)
      )
      aorsf::orsf(
        data = formula_data$data,
        formula = formula_data$formula,
        control = aorsf::orsf_control_survival(
          method = params$control_method,
          scale_x = params$control_scale_x,
          ties = params$control_ties,
          max_iter = params$control_max_iter
        ),
        n_tree = params$num_trees,
        n_split = params$num_split,
        n_retry = params$num_retry,
        n_thread = params$threads,
        mtry = max(1L, ceiling(sqrt(ncol(inputs$x)))),
        sample_with_replacement = params$sample_with_replacement,
        sample_fraction = params$sample_fraction,
        importance = params$importance,
        tree_seeds = tree_seeds,
        attach_data = params$attach_data,
        na_action = "fail",
        verbose_progress = FALSE
      )
    }
    )
  }

  prediction <- NULL
  effective_mtry <- NULL
  if (learner %in% sb_survreg_learner_ids()) {
    coefficients <- stats::coef(engine)
    coefficient_names <- colnames(stats::model.matrix(
      formula_data$formula, data = formula_data$data
    ))
    valid <- is.numeric(coefficients) &&
      identical(names(coefficients), coefficient_names) &&
      !anyNA(coefficients) && all(is.finite(coefficients)) &&
      is.numeric(engine$scale) && length(engine$scale) == 1L &&
      !is.na(engine$scale) && is.finite(engine$scale) && engine$scale > 0
    if (!valid) {
      sb_abort(
        paste0(
          "Fixed ", params$distribution,
          " SurvReg did not return finite prediction parameters."
        ),
        "stormbreaker_incomplete_survreg_fit"
      )
    }
    prediction <- list(
      intercept = unname(coefficients[[1L]]),
      coefficient = unname(coefficients[-1L]),
      scale = unname(engine$scale)
    )
  } else if (identical(learner, "rpart_survival")) {
    if (!inherits(engine, "rpart")) {
      sb_abort(
        "Fixed survival Rpart did not return a fitted tree.",
        "stormbreaker_incomplete_rpart_survival_fit"
      )
    }
  } else if (identical(learner, "rsf_survival")) {
    if (!inherits(engine, "rfsrc") || is.null(engine$forest)) {
      sb_abort(
        "Fixed RSF did not return a prediction forest.",
        "stormbreaker_incomplete_rsf_fit"
      )
    }
    effective_mtry <- as.integer(engine$mtry)
  } else if (identical(learner, "ctree_survival")) {
    if (!inherits(engine, "constparty")) {
      sb_abort(
        "Fixed survival Ctree did not return a fitted tree.",
        "stormbreaker_incomplete_ctree_survival_fit"
      )
    }
  } else if (identical(learner, "bart_survival")) {
    if (!inherits(engine, "survbart") || is.null(engine$treedraws)) {
      sb_abort(
        "Fixed survival BART did not return a prediction engine.",
        "stormbreaker_incomplete_bart_survival_fit"
      )
    }
  } else if (identical(learner, "aorsf_survival")) {
    if (!inherits(engine, "ObliqueForestSurvival")) {
      sb_abort(
        "Fixed aorsf survival did not return an oblique survival forest.",
        "stormbreaker_incomplete_aorsf_survival_fit"
      )
    }
    effective_mtry <- max(1L, ceiling(sqrt(ncol(inputs$x))))
  } else if (!inherits(engine, "cforest")) {
    sb_abort(
      "Fixed survival CForest did not return a fitted forest.",
      "stormbreaker_incomplete_cforest_survival_fit"
    )
  } else {
    effective_mtry <- max(1L, ceiling(sqrt(ncol(inputs$x))))
  }

  model <- sb_new_additional_survival_model(
    learner,
    names(inputs$x),
    params,
    engine,
    prediction,
    effective_mtry
  )
  risk <- sb_predict_additional_survival(model, inputs$x)
  if (length(risk) != nrow(inputs$x)) {
    sb_abort(
      paste0("Learner `", learner, "` did not produce one fitted risk per row."),
      "stormbreaker_invalid_additional_survival_prediction"
    )
  }
  model
}

sb_validate_additional_survival_model_x <- function(model, x) {
  if (!inherits(model, "stormbreaker_additional_survival_fixed") ||
    !model$learner %in% sb_additional_survival_learner_ids()) {
    sb_abort(
      "The supplied model is not a fitted additional survival learner.",
      "stormbreaker_bad_additional_survival_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  sb_require_additional_survival_dependency(model$learner)
  x
}

sb_survfit_last_cumulative_hazard <- function(fit) {
  value <- fit$cumhaz
  if (is.null(value) || length(value) == 0L) {
    return(0)
  }
  unname(as.numeric(value[[length(value)]]))
}

sb_partykit_survival_distributions <- function(model, x) {
  batch_size <- if (inherits(model$engine, "cforest")) {
    nrow(model$engine$data)
  } else {
    nrow(x)
  }
  # Bound native forest-weight allocation when assessment exceeds training size.
  # Each observation uses the same fitted forest and native probability method.
  if (nrow(x) > batch_size) {
    batches <- split(seq_len(nrow(x)), ceiling(seq_len(nrow(x)) / batch_size))
    distributions <- unlist(lapply(batches, function(ii) {
      stats::predict(model$engine, newdata = x[ii, , drop = FALSE], type = "prob")
    }), recursive = FALSE, use.names = FALSE)
  } else {
    distributions <- stats::predict(model$engine, newdata = x, type = "prob")
  }
  valid <- is.list(distributions) && length(distributions) == nrow(x) &&
    all(vapply(distributions, inherits, logical(1), what = "survfit"))
  if (!valid) {
    sb_abort(
      "The partykit survival learner produced invalid distributions.",
      "stormbreaker_invalid_additional_survival_prediction"
    )
  }
  distributions
}

sb_rsf_survival_curves <- function(model, x) {
  prediction <- withr::with_options(
    list(rf.cores = 1L, mc.cores = 1L),
    stats::predict(model$engine, newdata = x)
  )
  risk <- unname(as.numeric(prediction$predicted))
  survival <- as.matrix(prediction$survival)
  time <- as.numeric(prediction$time.interest)
  valid <- length(risk) == nrow(x) && !anyNA(risk) && all(is.finite(risk)) &&
    nrow(survival) == nrow(x) && ncol(survival) == length(time) &&
    length(time) > 0L && !anyNA(time) && all(is.finite(time)) &&
    !anyNA(survival) && all(is.finite(survival)) &&
    all(survival >= 0 & survival <= 1)
  if (!valid) {
    sb_abort(
      "Fixed RSF produced invalid mortality or survival predictions.",
      "stormbreaker_invalid_additional_survival_prediction"
    )
  }
  list(risk = risk, time = time, survival = survival)
}

sb_bart_survival_curves <- function(model, x) {
  x <- sb_validate_additional_survival_model_x(model, x)
  grid <- unname(as.numeric(model$engine$times))
  n <- nrow(x)
  k <- length(grid)
  matrix_x <- as.matrix(x)
  augmented <- cbind(
    t = rep(grid, times = n),
    matrix_x[rep(seq_len(n), each = k), , drop = FALSE]
  )
  removed <- model$engine$rm.const
  if (length(removed) > 0L) {
    augmented <- augmented[, removed, drop = FALSE]
  }
  prediction <- NULL
  predict_method <- get(
    "predict.survbart",
    envir = asNamespace("BART"),
    inherits = FALSE
  )
  invisible(utils::capture.output(
    prediction <- predict_method(
      model$engine,
      augmented,
      mc.cores = 1L,
      openmp = FALSE
    ),
    type = "output"
  ))
  survival <- matrix(
    unname(as.numeric(prediction$surv.test.mean)),
    nrow = n,
    ncol = k,
    byrow = TRUE
  )
  valid <- length(grid) > 0L && !anyNA(grid) && all(is.finite(grid)) &&
    all(grid > 0) && !is.unsorted(grid) &&
    nrow(survival) == n && ncol(survival) == k &&
    !anyNA(survival) && all(is.finite(survival)) &&
    all(survival >= 0 & survival <= 1)
  if (!valid) {
    sb_abort(
      "Fixed survival BART produced invalid posterior survival curves.",
      "stormbreaker_invalid_additional_survival_prediction"
    )
  }
  list(time = grid, survival = survival)
}

sb_bart_survival_negative_rmst <- function(curves) {
  interval <- diff(c(0, curves$time))
  previous <- cbind(
    1,
    curves$survival[, -ncol(curves$survival), drop = FALSE]
  )
  -as.vector(previous %*% interval)
}

sb_additional_survival_prediction_details <- function(model, x) {
  x <- sb_validate_additional_survival_model_x(model, x)
  details <- if (identical(model$learner, "bart_survival")) {
    list(curves = sb_bart_survival_curves(model, x))
  } else {
    list(distributions = sb_partykit_survival_distributions(model, x))
  }
  details$risk <- sb_predict_additional_survival(model, x, details)
  details
}

sb_predict_additional_survival <- function(model, x, prediction_details = NULL) {
  x <- sb_validate_additional_survival_model_x(model, x)
  learner <- model$learner
  risk <- if (learner %in% sb_survreg_learner_ids()) {
    location <- model$prediction$intercept +
      as.vector(as.matrix(x) %*% model$prediction$coefficient)
    -location
  } else {
    switch(learner,
    rpart_survival = {
      rate <- as.numeric(stats::predict(model$engine, newdata = x, type = "vector"))
      if (anyNA(rate) || any(!is.finite(rate)) || any(rate <= 0)) {
        sb_abort(
          "Fixed survival Rpart produced invalid relative rates.",
          "stormbreaker_invalid_additional_survival_prediction"
        )
      }
      log(rate)
    },
    rsf_survival = sb_rsf_survival_curves(model, x)$risk,
    ctree_survival = vapply(
      sb_partykit_survival_distributions(model, x),
      sb_survfit_last_cumulative_hazard,
      numeric(1)
    ),
    cforest_survival = vapply(
      if (is.null(prediction_details)) sb_partykit_survival_distributions(model, x) else
        prediction_details$distributions,
      sb_survfit_last_cumulative_hazard,
      numeric(1)
    ),
    bart_survival = sb_bart_survival_negative_rmst(
      if (is.null(prediction_details)) sb_bart_survival_curves(model, x) else
        prediction_details$curves
    ),
    aorsf_survival = {
      names(x) <- model$engine$get_names_x()
      stats::predict(
        model$engine,
        new_data = x,
        pred_type = "mort",
        n_thread = model$params$threads
      )
    }
    )
  }
  risk <- unname(as.numeric(risk))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      paste0("Learner `", learner, "` produced invalid risk predictions."),
      "stormbreaker_invalid_additional_survival_prediction"
    )
  }
  risk
}

sb_survival_curve_probability_at_horizon <- function(time, survival,
                                                     horizon) {
  eligible <- which(time <= horizon)
  if (length(eligible) == 0L) {
    return(1)
  }
  unname(as.numeric(survival[[max(eligible)]]))
}

sb_additional_direct_survival_probability <- function(model, x, horizon,
                                                      prediction_details = NULL) {
  x <- sb_validate_additional_survival_model_x(model, x)
  learner <- model$learner
  probability <- if (learner %in% sb_survreg_learner_ids()) {
    location <- model$prediction$intercept +
      as.vector(as.matrix(x) %*% model$prediction$coefficient)
    transformed <- (log(horizon) - location) / model$prediction$scale
    switch(model$params$distribution,
      weibull = exp(-exp(transformed)),
      exponential = exp(-exp(transformed)),
      lognormal = stats::pnorm(transformed, lower.tail = FALSE),
      loglogistic = stats::plogis(transformed, lower.tail = FALSE)
    )
  } else {
    switch(learner,
    rsf_survival = {
      curves <- sb_rsf_survival_curves(model, x)
      eligible <- which(curves$time <= horizon)
      if (length(eligible) == 0L) {
        rep(1, nrow(x))
      } else {
        unname(curves$survival[, max(eligible)])
      }
    },
    ctree_survival = vapply(
      sb_partykit_survival_distributions(model, x),
      function(distribution) {
        sb_survival_curve_probability_at_horizon(
          distribution$time,
          distribution$surv,
          horizon
        )
      },
      numeric(1)
    ),
    bart_survival = {
      curves <- if (is.null(prediction_details)) sb_bart_survival_curves(model, x) else
        prediction_details$curves
      eligible <- which(curves$time <= horizon)
      if (length(eligible) == 0L) {
        rep(1, nrow(x))
      } else {
        unname(curves$survival[, max(eligible)])
      }
    },
    cforest_survival = vapply(
      if (is.null(prediction_details)) sb_partykit_survival_distributions(model, x) else
        prediction_details$distributions,
      function(distribution) {
        sb_survival_curve_probability_at_horizon(
          distribution$time,
          distribution$surv,
          horizon
        )
      },
      numeric(1)
    ),
    aorsf_survival = {
      names(x) <- model$engine$get_names_x()
      stats::predict(
        model$engine,
        new_data = x,
        pred_type = "surv",
        pred_horizon = horizon,
        n_thread = model$params$threads
      )
    },
    sb_abort(
      "This learner requires Breslow calibration rather than direct survival.",
      "stormbreaker_bad_additional_survival_fit"
    )
    )
  }
  probability <- unname(as.numeric(probability))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      paste0("Learner `", learner, "` produced invalid horizon survival."),
      "stormbreaker_invalid_additional_survival_prediction"
    )
  }
  probability
}

sb_additional_direct_survival_source <- function(learner) {
  if (learner %in% sb_survreg_learner_ids()) {
    params <- sb_additional_survival_params(learner)
    return(paste0(
      "survreg_", params$distribution, "_distribution_direct"
    ))
  }
  switch(learner,
    rsf_survival = "randomForestSRC_survival_forest_direct",
    ctree_survival = "partykit_ctree_distribution_direct",
    cforest_survival = "partykit_cforest_distribution_direct",
    bart_survival = "BART_discrete_time_posterior_survival_direct",
    aorsf_survival = "aorsf_oblique_survival_forest_direct",
    NULL
  )
}

sb_predict_additional_survival_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon,
  prediction_details = NULL
) {
  inputs <- sb_validate_additional_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  if (identical(model$learner, "rpart_survival")) {
    analysis_risk <- sb_predict_additional_survival(model, inputs$x)
    assessment_risk <- sb_predict_additional_survival(
      model,
      x_assessment
    )
    baseline <- sb_breslow_baseline_hazard(
      inputs$time,
      inputs$event,
      analysis_risk,
      prediction_horizon
    )
    return(sb_survival_prediction_result(
      prediction_horizon,
      baseline,
      source = "rpart_survival_breslow_calibration",
      ties = "breslow",
      risk = assessment_risk
    ))
  }
  list(
    prediction_horizon = prediction_horizon,
    source = sb_additional_direct_survival_source(model$learner),
    ties = "not_applicable",
    probability = sb_additional_direct_survival_probability(
      model,
      x_assessment,
      prediction_horizon,
      prediction_details = prediction_details
    )
  )
}

sb_additional_survival_prediction_engine <- function(learner, engine) {
  if (identical(learner, "rpart_survival")) {
    engine[intersect(c("model", "x", "y", "where"), names(engine))] <- NULL
    return(engine)
  }
  if (identical(learner, "bart_survival")) {
    return(sb_bart_survival_prediction_engine(engine))
  }
  if (identical(learner, "aorsf_survival")) {
    return(engine)
  }
  if (!identical(learner, "rsf_survival")) {
    return(engine)
  }
  drop_fields <- c(
    "call", "yvar", "xvar", "subj", "proximity", "forest.wt",
    "case.depth", "distance", "membership", "tdc.membership", "inbag",
    "var.used", "imputed.indv", "imputed.data", "split.depth",
    "node.stats", "ensemble", "holdout.array", "empr.risk",
    "oob.empr.risk", "chf", "chf.oob", "predicted", "predicted.oob",
    "hazard", "hazard.oob", "survival", "survival.oob", "cif", "cif.oob",
    "err.rate"
  )
  engine[intersect(drop_fields, names(engine))] <- NULL
  engine
}

sb_extract_additional_survival_learner_state <- function(
  model, context = list()
) {
  if (!inherits(model, "stormbreaker_additional_survival_fixed")) {
    sb_bad_learner_state(
      "The fitted additional survival learner state is invalid."
    )
  }
  baseline_hazard <- NULL
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_additional_survival_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    if (identical(model$learner, "rpart_survival")) {
      baseline_hazard <- list(
        prediction_horizon = probability$prediction_horizon,
        time_unit = if ("time_unit" %in% names(context)) {
          context$time_unit
        } else {
          NA_character_
        },
        cumulative_baseline_hazard = probability$cumulative_baseline_hazard,
        baseline_survival = probability$baseline_survival,
        source = probability$source,
        ties = probability$ties
      )
    } else {
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
  }

  engine_raw <- NULL
  is_survreg <- model$learner %in% sb_survreg_learner_ids()
  if (!is_survreg) {
    engine <- sb_additional_survival_prediction_engine(
      model$learner,
      model$engine
    )
    engine_raw <- serialize(
      engine,
      connection = NULL,
      ascii = FALSE,
      xdr = TRUE,
      version = 3L
    )
  }
  state <- list(
    learner = model$learner,
    features = model$features,
    params = model$params,
    effective_mtry = model$effective_mtry,
    risk_reference = if (is_survreg) {
      paste0(model$params$distribution, "_aft_native")
    } else switch(model$learner,
      rpart_survival = "rpart_native",
      rsf_survival = "randomForestSRC_native",
      ctree_survival = "partykit_native",
      cforest_survival = "partykit_native",
      bart_survival = "BART_discrete_time_posterior",
      aorsf_survival = "aorsf_native"
    ),
    risk_formula = if (is_survreg) {
      "negative_log_time_location"
    } else switch(model$learner,
      rpart_survival = "log_relative_event_rate",
      rsf_survival = "mortality_sum_cumulative_hazard",
      ctree_survival = "terminal_distribution_cumulative_hazard",
      cforest_survival = "ensemble_distribution_cumulative_hazard",
      bart_survival = "negative_step_RMST_over_analysis_time_grid",
      aorsf_survival = "native_mortality_over_event_times"
    ),
    baseline_hazard = baseline_hazard,
    survival_probability = survival_probability
  )
  if (is_survreg) {
    state$intercept <- model$prediction$intercept
    state$coefficient <- model$prediction$coefficient
    state$scale <- model$prediction$scale
  } else {
    state$engine_raw <- engine_raw
  }
  if (model$learner %in% c(
    "ctree_survival",
    "cforest_survival"
  )) {
    state$training_data_retention <-
      "required_by_partykit_conditional_distribution_state"
  } else if (identical(model$learner, "rsf_survival")) {
    state$training_data_retention <-
      "event_time_metadata_required_by_randomForestSRC_forest"
  }
  structure(
    state,
    class = c(
      "stormbreaker_additional_survival_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_additional_survival_model <- function(state, learner) {
  if (learner %in% sb_survreg_learner_ids()) {
    valid <- is.numeric(state$intercept) && length(state$intercept) == 1L &&
      !is.na(state$intercept) && is.finite(state$intercept) &&
      is.numeric(state$coefficient) &&
      length(state$coefficient) == length(state$features) &&
      !anyNA(state$coefficient) && all(is.finite(state$coefficient)) &&
      is.numeric(state$scale) && length(state$scale) == 1L &&
      !is.na(state$scale) && is.finite(state$scale) && state$scale > 0
    if (!valid) {
      sb_bad_learner_state("The saved fixed SurvReg state is invalid.")
    }
    return(sb_new_additional_survival_model(
      learner,
      state$features,
      state$params,
      prediction = list(
        intercept = state$intercept,
        coefficient = state$coefficient,
        scale = state$scale
      )
    ))
  }
  if (!is.raw(state$engine_raw) || length(state$engine_raw) == 0L) {
    sb_bad_learner_state("The saved survival learner engine is invalid.")
  }
  engine <- tryCatch(
    unserialize(state$engine_raw),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved `", learner, "` engine could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_additional_survival_state",
        parent = error
      )
    }
  )
  valid_engine <- switch(learner,
    rpart_survival = inherits(engine, "rpart"),
    rsf_survival = inherits(engine, "rfsrc"),
    ctree_survival = inherits(engine, "constparty"),
    cforest_survival = inherits(engine, "cforest"),
    bart_survival = inherits(engine, "survbart"),
    aorsf_survival = inherits(engine, "ObliqueForestSurvival")
  )
  if (!valid_engine) {
    sb_bad_learner_state("The saved survival learner engine is invalid.")
  }
  sb_new_additional_survival_model(
    learner,
    state$features,
    state$params,
    engine = engine,
    effective_mtry = state$effective_mtry
  )
}

sb_validate_additional_survival_probability_state <- function(state, learner) {
  probability <- state$survival_probability
  expected_source <- sb_additional_direct_survival_source(learner)
  valid <- is.list(probability) &&
    is.numeric(probability$prediction_horizon) &&
    length(probability$prediction_horizon) == 1L &&
    !is.na(probability$prediction_horizon) &&
    is.finite(probability$prediction_horizon) &&
    probability$prediction_horizon > 0 &&
    is.character(probability$time_unit) &&
    length(probability$time_unit) == 1L &&
    !is.na(probability$time_unit) && nzchar(probability$time_unit) &&
    identical(probability$source, expected_source) &&
    identical(probability$ties, "not_applicable")
  if (!valid) {
    sb_bad_learner_state(
      "The saved direct survival-probability state is invalid."
    )
  }
  probability
}

sb_predict_additional_survival_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_additional_survival_learner_ids()) {
    sb_bad_learner_state(
      "The saved additional survival learner state is invalid."
    )
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  sb_require_additional_survival_dependency(learner)
  model <- sb_restore_additional_survival_model(state, learner)
  prediction_details <- NULL
  if (!is.null(state$survival_probability) && learner %in%
      c("bart_survival", "cforest_survival")) {
    warned <- FALSE
    prediction_details <- withCallingHandlers(
      sb_additional_survival_prediction_details(model, x),
      warning = function(w) warned <<- TRUE)
    risk <- prediction_details$risk
    if (warned) prediction_details <- NULL
  } else {
    risk <- sb_predict_additional_survival(model, x)
  }
  result <- list(risk = risk)
  if (!is.null(state$baseline_hazard)) {
    baseline <- state$baseline_hazard
    valid <- identical(learner, "rpart_survival") &&
      is.list(baseline) &&
      is.numeric(baseline$prediction_horizon) &&
      length(baseline$prediction_horizon) == 1L &&
      !is.na(baseline$prediction_horizon) &&
      is.finite(baseline$prediction_horizon) &&
      baseline$prediction_horizon > 0 &&
      is.character(baseline$time_unit) && length(baseline$time_unit) == 1L &&
      !is.na(baseline$time_unit) && nzchar(baseline$time_unit)
    if (!valid) {
      sb_bad_learner_state(
        "The saved Rpart survival baseline-hazard state is invalid."
      )
    }
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      risk,
      baseline$cumulative_baseline_hazard
    )
  } else if (!is.null(state$survival_probability)) {
    probability <- sb_validate_additional_survival_probability_state(
      state,
      learner
    )
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- sb_additional_direct_survival_probability(
      model,
      x,
      probability$prediction_horizon,
      prediction_details = prediction_details
    )
  }
  result
}

sb_additional_survival_learner_adapter <- function(
  learner,
  params = sb_additional_survival_params(learner)
) {
  if (!learner %in% sb_additional_survival_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_additional_survival(x, time, event, learner, params)
    },
    predict = sb_predict_additional_survival,
    predict_details = if (learner %in%
      c("bart_survival", "cforest_survival")) {
      sb_additional_survival_prediction_details
    } else NULL,
    predict_survival = sb_predict_additional_survival_survival,
    extract_state = function(model, context = list()) {
      sb_extract_additional_survival_learner_state(model, context)
    },
    predict_state = sb_predict_additional_survival_learner_state
  )
}

# ---- model-cox ----

# ---- cox-adapter ----

sb_fit_coxph <- function(x, time, event) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time < 0)) {
    sb_abort(
      "`time` must contain non-negative values.",
      "stormbreaker_bad_adapter_time"
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  if (ncol(x) >= nrow(x)) {
    sb_abort(
      paste0(
        "The unpenalized Cox learner cannot safely fit ", ncol(x),
        " predictors with ", nrow(x),
        " patients; use a fold-local selector or glmnet Cox."
      ),
      "stormbreaker_high_dimensional_cox"
    )
  }

  time_name <- ".stormbreaker_time"
  while (time_name %in% names(x)) {
    time_name <- paste0(time_name, "_")
  }
  event_name <- ".stormbreaker_event"
  while (event_name %in% c(names(x), time_name)) {
    event_name <- paste0(event_name, "_")
  }
  analysis <- x
  analysis[[time_name]] <- time
  analysis[[event_name]] <- event
  rhs <- paste(vapply(names(x), sb_quote_name, ""), collapse = " + ")
  response <- paste0(
    "survival::Surv(",
    sb_quote_name(time_name), ", ",
    sb_quote_name(event_name), ")"
  )
  formula <- stats::as.formula(
    paste(response, "~", rhs),
    env = parent.frame()
  )

  model <- withCallingHandlers(
    survival::coxph(
      formula,
      data = analysis,
      ties = "efron",
      singular.ok = FALSE,
      x = TRUE,
      y = TRUE,
      model = TRUE
    ),
    warning = function(warning) {
      sb_abort(
        paste0("The Cox learner raised a warning: ",
               conditionMessage(warning)),
        "stormbreaker_cox_warning",
        parent = warning
      )
    }
  )

  if (NROW(model$y) != nrow(analysis)) {
    sb_abort(
      "The Cox learner dropped analysis patients.",
      "stormbreaker_rows_dropped"
    )
  }
  coefficients <- stats::coef(model)
  if (length(coefficients) < 1L || anyNA(coefficients) ||
      any(!is.finite(coefficients))) {
    sb_abort(
      "The Cox learner produced missing or non-finite coefficients.",
      "stormbreaker_invalid_cox_fit"
    )
  }

  model
}

sb_predict_coxph <- function(model, x) {
  x <- sb_validate_predictor_frame(x, "x")
  expected_features <- attr(model, "stormbreaker_features", exact = TRUE)
  if (!identical(names(x), expected_features)) {
    sb_abort(
      "Prediction columns do not match the fitted Cox learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  risk <- unname(stats::predict(
    model,
    newdata = x,
    type = "lp",
    reference = "zero"
  ))
  if (length(risk) != nrow(x) || anyNA(risk) ||
      any(!is.finite(risk))) {
    sb_abort(
      "The Cox learner produced incomplete or non-finite risk predictions.",
      "stormbreaker_invalid_cox_prediction"
    )
  }
  risk
}

sb_extract_coxph_coefficients <- function(model) {
  features <- attr(model, "stormbreaker_features", exact = TRUE)
  coefficient <- unname(stats::coef(model))
  data.frame(
    feature = features,
    coefficient = coefficient,
    nonzero = coefficient != 0,
    direction = ifelse(
      coefficient > 0,
      "higher_hazard",
      ifelse(coefficient < 0, "lower_hazard", "zero")
    ),
    stringsAsFactors = FALSE
  )
}

sb_predict_coxph_survival <- function(model, x_analysis, time, event,
                                      x_assessment, prediction_horizon) {
  x_analysis <- sb_validate_predictor_frame(x_analysis, "x_analysis")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x_analysis))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x_analysis))
  expected_features <- attr(model, "stormbreaker_features", exact = TRUE)
  if (!identical(names(x_analysis), expected_features)) {
    sb_abort(
      "Analysis columns do not match the fitted Cox learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    time
  )
  risk <- sb_predict_coxph(model, x_assessment)
  baseline <- withCallingHandlers(
    survival::basehaz(model, centered = FALSE),
    warning = function(warning) {
      sb_abort(
        paste0(
          "The Cox baseline-hazard calculation raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_baseline_hazard_warning",
        parent = warning
      )
    }
  )
  eligible <- baseline$time <= prediction_horizon
  cumulative_hazard <- if (any(eligible)) {
    unname(utils::tail(baseline$hazard[eligible], 1L))
  } else {
    0
  }
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = "coxph_basehaz_zero_reference",
    ties = model$method,
    risk = risk
  )
}

sb_coxph_learner_adapter <- function() {
  list(
    id = "coxph",
    fit = function(x, time, event) {
      model <- sb_fit_coxph(x, time, event)
      attr(model, "stormbreaker_features") <- names(x)
      model
    },
    predict = sb_predict_coxph,
    predict_survival = sb_predict_coxph_survival,
    extract_coefficients = sb_extract_coxph_coefficients,
    extract_state = function(model, context = list()) {
      sb_extract_cox_learner_state(
        model,
        learner = "coxph",
        coefficients = sb_extract_coxph_coefficients(model),
        predict_survival = sb_predict_coxph_survival,
        context = context
      )
    },
    predict_state = sb_predict_cox_learner_state
  )
}

# ---- coxphf-adapter ----

sb_coxphf_learner_ids <- function() {
  "coxphf_firth"
}

sb_coxphf_available <- function() {
  requireNamespace("coxphf", quietly = TRUE)
}

sb_require_coxphf_dependency <- function() {
  if (!sb_coxphf_available()) {
    sb_abort(
      "Package `coxphf` is required for Firth-corrected Cox regression.",
      "stormbreaker_missing_coxphf"
    )
  }
  invisible(TRUE)
}

sb_coxphf_learner_params <- function(
  learner = sb_coxphf_learner_ids()
) {
  if (!learner %in% sb_coxphf_learner_ids()) {
    sb_abort(
      paste0("Unknown Firth Cox learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    firth = TRUE,
    penalty = 0.5,
    profile_likelihood = FALSE,
    max_iterations = 100L,
    max_step_halvings = 5L,
    convergence_epsilon = 1e-6,
    gradient_convergence = 1e-4,
    maximum_standardized_step = 0.5,
    ties = "breslow"
  )
}

sb_coxphf_parameter_semantics <- function(learner) {
  params <- sb_coxphf_learner_params(learner)
  paste0(
    "Firth bias-reduced Cox regression with penalty ", params$penalty,
    ", Breslow ties, Wald summaries without profile-likelihood fitting, and ",
    params$max_iterations, " maximum iterations fixed before resampling"
  )
}

sb_coxphf_stop_on_warning <- function(code, stage) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      sb_abort(
        paste0(
          "The Firth Cox learner raised a ", stage, " warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_coxphf_warning",
        parent = warning
      )
    }
  )
}

sb_coxphf_internal_names <- function(p) {
  paste0(".stormbreaker_coxphf_x", seq_len(p))
}

sb_coxphf_formula_data <- function(inputs) {
  analysis <- inputs$x
  internal_features <- sb_coxphf_internal_names(ncol(analysis))
  names(analysis) <- internal_features
  analysis$.stormbreaker_coxphf_time <- inputs$time
  analysis$.stormbreaker_coxphf_event <- inputs$event
  rownames(analysis) <- as.character(seq_len(nrow(analysis)))
  response <- paste0(
    "Surv(.stormbreaker_coxphf_time, .stormbreaker_coxphf_event)"
  )
  formula <- stats::as.formula(
    paste(response, "~", paste(internal_features, collapse = " + ")),
    env = asNamespace("coxphf")
  )
  list(
    data = analysis,
    formula = formula,
    internal_features = internal_features
  )
}

sb_validate_coxphf_model <- function(model) {
  valid <- inherits(model, "stormbreaker_coxphf_model") && is.list(model) &&
    identical(model$learner, sb_coxphf_learner_ids()) &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.character(model$internal_features) &&
    length(model$internal_features) == length(model$features) &&
    is.numeric(model$coefficient) &&
    length(model$coefficient) == length(model$features) &&
    !anyNA(model$coefficient) && all(is.finite(model$coefficient)) &&
    is.list(model$params)
  if (!isTRUE(valid)) {
    sb_abort(
      "The supplied model is not a valid fitted Firth Cox learner.",
      "stormbreaker_bad_coxphf_fit"
    )
  }
  model
}

sb_fit_coxphf_learner <- function(x, time, event, learner, params) {
  sb_require_coxphf_dependency()
  inputs <- sb_validate_classical_survival_inputs(
    x,
    time,
    event,
    learner
  )
  prepared <- sb_coxphf_formula_data(inputs)
  engine <- sb_coxphf_stop_on_warning(
    coxphf::coxphf(
      formula = prepared$formula,
      data = prepared$data,
      pl = params$profile_likelihood,
      maxit = params$max_iterations,
      maxhs = params$max_step_halvings,
      epsilon = params$convergence_epsilon,
      gconv = params$gradient_convergence,
      maxstep = params$maximum_standardized_step,
      firth = params$firth,
      penalty = params$penalty
    ),
    "fitting"
  )
  coefficient <- unname(engine$coefficients)
  valid <- inherits(engine, "coxphf") &&
    identical(names(engine$coefficients), prepared$internal_features) &&
    length(coefficient) == ncol(inputs$x) && !anyNA(coefficient) &&
    all(is.finite(coefficient)) && identical(engine$method.ties, params$ties)
  if (!isTRUE(valid)) {
    sb_abort(
      "The Firth Cox learner returned invalid coefficients.",
      "stormbreaker_invalid_coxphf_fit"
    )
  }
  sb_validate_coxphf_model(structure(
    list(
      learner = learner,
      features = names(inputs$x),
      internal_features = prepared$internal_features,
      coefficient = coefficient,
      params = params,
      iterations = as.integer(engine$iter)
    ),
    class = "stormbreaker_coxphf_model"
  ))
}

sb_validate_coxphf_model_x <- function(model, x) {
  model <- sb_validate_coxphf_model(model)
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted Firth Cox learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  list(model = model, x = x)
}

sb_predict_coxphf_learner <- function(model, x) {
  validated <- sb_validate_coxphf_model_x(model, x)
  risk <- unname(as.vector(
    as.matrix(validated$x) %*% validated$model$coefficient
  ))
  if (length(risk) != nrow(validated$x) || anyNA(risk) ||
    any(!is.finite(risk))) {
    sb_abort(
      "The Firth Cox learner produced invalid risk predictions.",
      "stormbreaker_invalid_coxphf_prediction"
    )
  }
  risk
}

sb_extract_coxphf_coefficients <- function(model) {
  model <- sb_validate_coxphf_model(model)
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

sb_predict_coxphf_survival <- function(
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
  analysis_risk <- sb_predict_coxphf_learner(model, inputs$x)
  assessment_risk <- sb_predict_coxphf_learner(model, x_assessment)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = "coxphf_firth_breslow",
    ties = model$params$ties,
    risk = assessment_risk
  )
}

sb_extract_coxphf_state <- function(model, context = list()) {
  model <- sb_validate_coxphf_model(model)
  sb_extract_cox_learner_state(
    model,
    learner = model$learner,
    coefficients = sb_extract_coxphf_coefficients(model),
    predict_survival = sb_predict_coxphf_survival,
    context = context,
    params = model$params,
    audit = list(
      native_linear_predictors_used = FALSE,
      fitting_iterations = model$iterations,
      estimation = "Firth_penalized_partial_likelihood"
    )
  )
}

sb_coxphf_learner_adapter <- function(
  learner,
  params = sb_coxphf_learner_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_coxphf_learner(x, time, event, learner, params)
    },
    predict = sb_predict_coxphf_learner,
    predict_survival = sb_predict_coxphf_survival,
    extract_coefficients = sb_extract_coxphf_coefficients,
    extract_state = sb_extract_coxphf_state,
    predict_state = sb_predict_cox_learner_state
  )
}

# ---- diagnostics-survival ----

# ---- algorithm-coxboost ----

sb_coxboost_diagnostic <- function(fit, type) {
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("coefficient_path", "nonzero_path", "cv_curve")) {
    sb_abort("Choose coefficient_path, nonzero_path or cv_curve for CoxBoost.",
      "stormbreaker_unavailable_diagnostic")
  }
  id <- fit$pipeline$learner
  native <- NULL
  if (type == "cv_curve") {
    if (id != "coxboost_cv") sb_abort("Fixed CoxBoost has no CV search.",
      "stormbreaker_unavailable_diagnostic")
    audit <- fit$learner_state$audit
    tab <- audit$tuning_curve
    stop <- audit$optimal_step
    if (!is.data.frame(tab) || !all(c("step", "mean_log_likelihood", "standard_error") %in% names(tab)) ||
        !nrow(tab) || !all(vapply(tab, function(z) is.numeric(z) && all(is.finite(z)), logical(1))) ||
        any(tab$standard_error < 0) || !identical(as.integer(tab$step), seq_len(nrow(tab)) - 1L) ||
        length(stop) != 1L || !stop %in% tab$step) sb_abort(
      "Original CoxBoost CV curve is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    native <- tab
    measure <- "Cross-validated partial log-likelihood"
    interpretation <- paste0("Original ", audit$tuning_type,
      " CV search; ribbon is +/- one standard error, not a confidence interval.\n",
      "Penalty was chosen before this step-count search; not validation of the complete pipeline.")
  } else {
    model <- fit$diagnostic_artifacts$learner
    if (is.null(model)) sb_abort(
      "Original CoxBoost path was not captured; enable capture_diagnostics at final fitting.",
      "stormbreaker_missing_diagnostic_artifacts")
    engine <- model$engine
    if (!inherits(engine, "CoxBoost") || !identical(model$features, engine$xnames)) sb_abort(
      "Captured CoxBoost engine or feature ordering is invalid.", "stormbreaker_invalid_diagnostic")
    stop <- engine$stepno
    native <- as.matrix(engine$coefficients)
    scale <- engine$sdx
    if (length(stop) != 1L || !is.finite(stop) || stop < 1 || stop != floor(stop) ||
        !identical(dim(native), c(as.integer(stop + 1L), length(model$features))) ||
        any(!is.finite(native)) || length(scale) != ncol(native) || any(!is.finite(scale)) || any(scale <= 0)) sb_abort(
      "Original CoxBoost coefficient path is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    # Use the same conversion as the original package fitting adapter, not a
    # newly fitted model or a different coefficient convention from coef().
    beta <- sweep(native, 2L, scale, "/")
    if (!isTRUE(all.equal(as.numeric(beta[nrow(beta), ]),
      as.numeric(fit$learner_state$coefficient), tolerance = 1e-12))) sb_abort(
      "Path endpoint disagrees with the saved prediction coefficients.", "stormbreaker_invalid_diagnostic")
    tab <- if (type == "coefficient_path") {
      data.frame(step = rep(seq.int(0L, stop), times = ncol(beta)),
        feature = rep(model$features, each = nrow(beta)), coefficient = as.vector(beta))
    } else {
      data.frame(step = seq.int(0L, stop), nonzero = rowSums(beta != 0))
    }
    measure <- if (type == "coefficient_path") "Original-input CoxBoost coefficient path" else "CoxBoost nonzero predictor count"
    interpretation <- "Original final-fit path including step zero; no additional training or change in stopping.\nCoefficients follow the saved prediction scale, not hazard ratios or uncertainty intervals."
  }
  structure(list(type = paste0("coxboost_", type), table = tab, native_table = native,
    stopping_iteration = as.integer(stop), learner = id, pipeline_id = fit$pipeline$id,
    measure = measure, interpretation = interpretation,
    source = "original CoxBoost artifacts; no refitting"), class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_coxboost_diagnostic <- function(x) {
  step <- feature <- coefficient <- nonzero <- mean_log_likelihood <- standard_error <- NULL
  tab <- x$table
  if (x$type == "coxboost_cv_curve") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(step, mean_log_likelihood)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_log_likelihood - standard_error,
        ymax = mean_log_likelihood + standard_error), fill = "#2878A0", alpha = .15) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::labs(y = "Partial log-likelihood (higher is better)")
  } else if (x$type == "coxboost_coefficient_path") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(step, coefficient, colour = feature)) +
      ggplot2::geom_line() + ggplot2::labs(y = "Coefficient on saved prediction scale", colour = "Feature")
  } else {
    p <- ggplot2::ggplot(tab, ggplot2::aes(step, nonzero)) +
      ggplot2::geom_step(colour = "#2878A0") +
      ggplot2::scale_y_continuous(breaks = function(limits) unique(round(pretty(limits))), limits = c(0, NA)) +
      ggplot2::labs(y = "Number of nonzero predictors")
  }
  p + ggplot2::geom_vline(xintercept = x$stopping_iteration, colour = "#D77935", linetype = 2) +
    ggplot2::labs(title = x$measure, subtitle = paste0(x$learner, "; original stopping step: ", x$stopping_iteration),
      x = "Boosting step", caption = x$interpretation) + sb_viz_theme()
}

# ---- algorithm-penaft ----

#' @rdname sb_flexsurv_diagnostics
#' @export
sb_penaft_diagnostics <- function(fit, type = "coefficients") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!identical(s$learner, "penaft_gehan_elastic_net")) sb_abort(
    "This diagnostic requires penAFT.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("coefficients", "cv_curve", "coefficient_path")) sb_abort(
    "Choose coefficients, cv_curve or coefficient_path.", "stormbreaker_unsupported_diagnostic")
  path <- NULL
  intercept_path <- NULL
  if (type == "coefficients") tab <- data.frame(feature = s$features, coefficient = s$coefficient) else {
    path <- fit$diagnostic_artifacts$learner$diagnostic_path
    if (is.null(path)) sb_abort("Original penAFT paths were not captured.", "stormbreaker_missing_diagnostic_artifacts")
    f <- path$full_fit
    beta <- as.matrix(f$beta)
    n <- length(f$lambda)
    j <- path$selected_index
    if (!n || !identical(dim(beta), c(as.integer(length(s$features)), as.integer(n))) ||
        length(path$cv_error) != n || any(!is.finite(path$cv_error)) ||
        any(!is.finite(beta)) || any(!is.finite(f$lambda)) || any(f$lambda <= 0) ||
        length(j) != 1L || !is.finite(j) || j < 1 || j > n) sb_abort(
      "Original penAFT path is invalid.", "stormbreaker_invalid_diagnostic")
    intercept_path <- rep(0, n)
    if (isTRUE(f$standardize)) {
      if (length(f$X.sd) != nrow(beta) || length(f$X.mean) != nrow(beta) ||
          any(!is.finite(f$X.sd)) || any(f$X.sd <= 0) || any(!is.finite(f$X.mean))) sb_abort(
        "Original penAFT normalization is invalid.", "stormbreaker_invalid_diagnostic")
      beta <- sweep(beta, 1, f$X.sd, "/")
      intercept_path <- -as.numeric(crossprod(f$X.mean, beta))
    }
    if (!isTRUE(all.equal(as.numeric(beta[, j]), s$coefficient, tolerance = 1e-12)) ||
        !isTRUE(all.equal(intercept_path[j], s$intercept, tolerance = 1e-12)) ||
        f$lambda[j] != s$selected_lambda || path$cv_error[j] != s$selected_cv_error) sb_abort(
      "Original selected penAFT path does not match the saved prediction state.", "stormbreaker_invalid_diagnostic")
    tab <- if (type == "cv_curve") data.frame(lambda = f$lambda, score = path$cv_error,
      selected = seq_len(n) == j) else data.frame(feature = rep(s$features, each = n),
        lambda = rep(f$lambda, length(s$features)), coefficient = as.numeric(t(beta)))
  }
  structure(list(type = paste0("penaft_", type), table = tab, native_path = path,
    intercept = s$intercept, intercept_path = intercept_path, features = s$features,
    selected_lambda = s$selected_lambda, learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "original penAFT Gehan elastic-net model and captured linear-predictor CV path",
    interpretation = paste("Original input-unit centered log-time coefficients, not log hazard ratios; risk has the opposite sign.",
      "Centering intercepts are retained separately. Native standardized coefficients and all path points are preserved.",
      "The marked lambda is the original selection; CV score is not AUC or untouched-test performance. No invented uncertainty.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_penaft_diagnostic <- function(x) {
  feature <- coefficient <- lambda <- score <- NULL
  d <- x$table
  if (x$type == "penaft_coefficients") {
    d$feature <- factor(d$feature, levels = rev(x$features))
    p <- ggplot2::ggplot(d, ggplot2::aes(coefficient, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = .4) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::labs(title = "Gehan elastic-net coefficients", x = "Log-time coefficient", y = NULL)
  } else if (x$type == "penaft_cv_curve") p <- ggplot2::ggplot(d, ggplot2::aes(lambda, score)) +
    ggplot2::geom_line(colour = "#2878A0") +
    ggplot2::geom_point(data = d[d$selected, , drop = FALSE], colour = "#C96942", size = 2.5) +
    ggplot2::scale_x_log10() + ggplot2::labs(title = "Gehan elastic-net selection", x = "Lambda", y = "CV linear-predictor score") else
      p <- ggplot2::ggplot(d, ggplot2::aes(lambda, coefficient, colour = feature)) +
        ggplot2::geom_line() + ggplot2::geom_vline(xintercept = x$selected_lambda, linetype = 2, colour = "grey60") +
        ggplot2::scale_x_log10() + ggplot2::labs(title = "Gehan elastic-net coefficient path", x = "Lambda", y = "Log-time coefficient", colour = NULL)
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = if (x$type == "penaft_coefficients") max(5, 2 + .3 * length(x$features)) else 5)
  p
}
