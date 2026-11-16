#' Nonlinear and boosting diagnostics
#' @name sb_gam_effects
#' @description Read native BART, GAM, MARS, neural, boosting and rule-model diagnostics.
#' @rawRd \arguments{
#' \item{fit}{
#' \strong{For \code{sb_gam_effects()}:} A final mgcv GAM fit with capture_diagnostics = TRUE.
#' 
#' \strong{For \code{sb_bart_posterior()}:} Successful binary, discrete-time survival or AFT BART final fit.
#' 
#' \strong{For \code{sb_earth_diagnostics()}:} A successful earth binary final fit.
#' 
#' \strong{For \code{sb_neural_diagnostics()}:} A successful nnet or multilayer neuralnet final fit.
#' 
#' \strong{For \code{sb_boosting_diagnostics()}:} A final supported mboost, GBM, XGBoost, LightGBM, CatBoost or
#' CoxBoost, AdaBoost, LogitBoost or hinge-boosting model. Saved state is
#' sufficient except binary GBM training loss, CoxBoost paths and hinge
#' iteration diagnostics require diagnostic capture.
#' 
#' \strong{For \code{sb_rule_ensemble()}:} A final binary or survival PRE model. Its saved prediction state
#' already contains the native engine; extra diagnostic capture is not needed.
#' }
#' \item{features}{
#' Predictor names; NULL displays every fitted predictor.
#' }
#' \item{n}{
#' Grid points spanning each predictor's observed training range.
#' }
#' \item{level}{
#' Pointwise normal-approximation interval level.
#' }
#' \item{newdata}{
#' Explicit query data. No training queries are implicit.
#' }
#' \item{type}{
#' \strong{For \code{sb_earth_diagnostics()}:} basis_coefficients or pruning_gcv.
#' 
#' \strong{For \code{sb_neural_diagnostics()}:} topology or weights.
#' 
#' \strong{For \code{sb_boosting_diagnostics()}:} training_loss or component_usage for mboost; training_loss,
#' cv_curve or importance for GBM where recorded; importance for XGBoost;
#' importance/tree_size for LightGBM and importance/interaction for CatBoost.
#' For CoxBoost use coefficient_path, nonzero_path or cv_curve.
#' AdaBoost/LogitBoost support split_usage; AdaBoost also supports tree_weights
#' and tree_size, while LogitBoost supports stump_sequence.
#' 
#' \strong{For \code{sb_rule_ensemble()}:} Either rule_coefficients or rule_cv_curve.
#' }
#' \item{max_pairs}{
#' Maximum native CatBoost interaction pairs displayed (20 by
#' default). Only the plot is limited; the complete reported table is retained.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_gam_effects()}}{
#' Recover original GAM smooth effects and conditional uncertainty
#' 
#' \strong{Returns:} An algorithm diagnostic with original centered smooth contributions,
#' conditional standard errors, and observed predictor values without IDs.
#' Values can still be sensitive. Other predictors are held at training
#' medians; this adapter supports only the package's additive univariate GAMs.
#' These are link-scale associations, not probabilities or causal effects.
#' }
#' 
#' \subsection{\code{sb_bart_posterior()}}{
#' Inspect original BART posterior predictions
#' 
#' \strong{Returns:} All retained posterior draws and their mean and pointwise 95 percent
#' equal-tail credible limits (empirical quantiles, type 7). Binary values are
#' class-1 probabilities; survival values are survival probabilities at native
#' time points; AFT values are latent log-time locations, not realized times.
#' These describe posterior parameter uncertainty, not prediction intervals
#' for individual outcomes, validation performance or simultaneous bands.
#' }
#' 
#' \subsection{\code{sb_earth_diagnostics()}}{
#' Inspect the original MARS logistic basis and pruning criterion
#' 
#' \strong{Returns:} An algorithm diagnostic retaining original selected logistic basis
#' coefficients, or the native basis-pruning GCV sequence. Logistic parameters
#' come from glm.coefficients, not the separate least-squares coefficients.
#' Native basis labels are retained for display; full-precision cuts and
#' directions remain available separately and must not be reconstructed from
#' rounded labels. The GCV sequence needs original captured artifacts and is
#' the basis-selection criterion, not held-out binomial deviance, calibration
#' or the final logistic model's cross-validation performance. Term counts
#' include the intercept. No refitting, surrogate model or uncertainty is added.
#' }
#' 
#' \subsection{\code{sb_neural_diagnostics()}}{
#' Inspect original neural-network connections
#' 
#' \strong{Returns:} An algorithm diagnostic from compact original connection weights.
#' Bias connections and all hidden units are retained. Topology edges indicate
#' connectivity, not importance or causality. Weight matrices use the model's
#' original standardized inputs; weights are not log odds, SHAP or input
#' importance, and are not comparable as feature effects across layers.
#' Bias nodes are repeated visually per receiving layer for readability and
#' represent the same constant value one. No new fit or optimization-history
#' reconstruction is performed. The saved \code{parameters} permit native forward
#' calculation with logistic activations at each hidden and output layer.
#' nnet's native activation is exactly zero below -15 and one above 15;
#' its saturation thresholds are retained rather than replaced by plogis.
#' }
#' 
#' \subsection{\code{sb_boosting_diagnostics()}}{
#' Inspect saved boosting loss and component selection
#' 
#' \strong{Returns:} An algorithm diagnostic. training_loss retains the native objective
#' at iterations zero through the original stopping iteration. Loss scale and
#' sign are unchanged, including the negative concordance objective. Values
#' from different families or samples are not directly comparable. This is
#' in-bag training loss, not an external validation or early-stopping curve.
#' 
#' component_usage retains every original component, including unused terms
#' and an intercept when present, and a per-iteration selection table. Counts
#' describe how often a component was updated, not its coefficient magnitude,
#' predictive importance or across-resample stability. Internal feature aliases
#' are mapped back to their saved input names. BlackBoost updates a single
#' tree base learner; it is not given a misleading variable-selection plot.
#' Neither diagnostic changes stopping, fits, coefficients or predictions.
#' 
#' GBM keeps final-fit training errors separate from its original CV search.
#' XGBoost returns native gain/cover/frequency shares for tree boosters, or
#' signed coefficients for linear boosters (bias excluded). AFT coefficients
#' are on the log-time scale, not the package's reversed risk scale. Missing
#' zero-use feature rows are restored as zeros, with the native table retained.
#' Original pipelines without evaluation logs do not get invented curves.
#' LightGBM retains unscaled gain/cover/split counts and counts actual saved
#' trees, which can be fewer than configured rounds. CatBoost interaction is
#' a model-internal strength score, not SHAP interaction or a causal effect.
#' Native zero-based feature indices are mapped to the saved input names;
#' unreported pairs are not assigned fabricated zero scores.
#' A single input has no distinct feature pairs and returns an empty table
#' without calling the native pair extractor. Native extraction failures for
#' other models remain explicit unavailable diagnostics with their cause.
#' CoxBoost supports coefficient_path and nonzero_path with diagnostic
#' capture, and cv_curve from saved tuning for its CV variant without capture.
#' Coefficients retain the original prediction adapter's scale conversion;
#' the CV curve is conditional on the previously selected penalty.
#' AdaBoost and LogitBoost expose split_usage without capture. AdaBoost also
#' exposes tree_weights and tree_size for actual fitted trees, including
#' native convergence before the configured maximum. LogitBoost exposes its
#' stump_sequence with exact native thresholds and signed log-odds votes.
#' These are model-specific split rules, not clinical decision thresholds.
#' Hinge boosting supports training_loss, component_usage and coefficient_path
#' from newly captured original iterations. Native losses start after update
#' one; no initial loss is invented. Coefficient paths include the initial
#' margin equation and its centering-adjusted intercept, not probability or
#' log-odds coefficients. Historical compact fits cannot supply these paths.
#' }
#' 
#' \subsection{\code{sb_rule_ensemble()}}{
#' Inspect an original PRE rule ensemble and its coefficient tuning
#' 
#' \strong{Returns:} An algorithm diagnostic. For coefficients, table contains every
#' original term, its exact native description, coefficient and term type,
#' including zero coefficients and the intercept where present. Linear terms
#' retain native winsorization bounds and coefficient rescaling. Their native
#' interval description denotes clipping, not an interval-indicator rule.
#' Coefficient size is
#' not comparable feature importance. The plot excludes zero coefficients and
#' the intercept, but the complete table remains unchanged.
#' 
#' For tuning, table contains the saved coefficient-only cross-validation
#' sequence, with its original standard errors and selected lambda. Rules were
#' induced on the complete development data before this internal validation;
#' this is not validation of the complete rule-learning pipeline. No model is
#' refitted. Rule descriptions are returned as text, not evaluated.
#' }
#' }
NULL

# ---- model-gam ----

# ---- algorithm-gam-effects ----

#' @rdname sb_gam_effects
#' @export
sb_gam_effects <- function(fit, features = NULL, n = 100L, level = .95) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!fit$pipeline$learner %in% sb_mgcv_gam_learner_ids()) sb_abort(
    "Smooth-effect diagnostics require an original mgcv GAM learner.",
    "stormbreaker_unsupported_diagnostic")
  artifact <- fit$diagnostic_artifacts
  if (is.null(artifact)) sb_abort(
    "Original GAM and reference predictors were not captured; enable capture_diagnostics at final fitting.",
    "stormbreaker_missing_diagnostic_artifacts")
  model <- artifact$learner
  if (!inherits(model, "stormbreaker_mgcv_gam")) sb_abort(
    "Captured GAM model is invalid.", "stormbreaker_invalid_diagnostic")
  n <- sb_viz_positive_integer(n, "n")
  if (n < 2 || !is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) sb_abort(
    "Use n >= 2 and an interval level strictly between zero and one.",
    "stormbreaker_bad_visualization_argument")
  if (is.null(features)) features <- model$features
  if (!is.character(features) || !length(features) || anyNA(features) ||
      anyDuplicated(features) || !all(features %in% model$features)) sb_abort(
    "Specify distinct fitted predictor names.", "stormbreaker_bad_visualization_argument")
  ref <- sb_select_final_feature_frame(artifact$predictors, model$features,
    "Captured development data", "GAM reference")
  engine <- sb_validate_mgcv_engine(model$engine, model$learner, model$engine_features)
  medians <- vapply(ref, stats::median, numeric(1))
  tab <- support <- list()
  z <- stats::qnorm((1 + level) / 2)
  for (feature in features) {
    i <- match(feature, model$features)
    grid <- as.data.frame(as.list(medians), check.names = FALSE)[rep(1L, n), , drop = FALSE]
    grid[[feature]] <- seq(min(ref[[feature]]), max(ref[[feature]]), length.out = n)
    values <- grid[[feature]]
    names(grid) <- model$engine_features
    term <- engine$smooth[[i]]$label
    predicted <- mgcv::predict.gam(engine, newdata = grid, type = "terms", se.fit = TRUE,
      terms = term, unconditional = FALSE)
    effect <- as.numeric(predicted$fit[, term])
    se <- as.numeric(predicted$se.fit[, term])
    if (length(effect) != n || length(se) != n || any(!is.finite(effect)) ||
        any(!is.finite(se)) || any(se < 0)) sb_abort(
      "The original GAM could not return finite smooth contributions and uncertainty.",
      "stormbreaker_invalid_diagnostic")
    tab[[feature]] <- data.frame(feature = feature, value = values, effect = effect,
      standard_error = se, lower = effect - z * se, upper = effect + z * se)
    support[[feature]] <- data.frame(feature = feature, value = ref[[feature]])
  }
  structure(list(type = "gam_effects", table = do.call(rbind, tab),
    support = do.call(rbind, support), level = level, learner = model$learner,
    pipeline_id = fit$pipeline$id, score_scale = if (fit$task_type == "binary")
      "Centered smooth contribution to log odds" else "Centered smooth contribution to log hazard",
    source = "original mgcv predict(type='terms'); conditional on fitted smoothing parameters"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_gam_effects <- function(x) {
  value <- effect <- lower <- upper <- feature <- NULL
  ggplot2::ggplot(x$table, ggplot2::aes(value, effect)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#B9D7E5", alpha = .6) +
    ggplot2::geom_hline(yintercept = 0, colour = "#909090", linetype = 2) +
    ggplot2::geom_line(colour = "#2878A0") +
    ggplot2::geom_rug(data = x$support, ggplot2::aes(x = value), inherit.aes = FALSE,
      sides = "b", alpha = .25) +
    ggplot2::facet_wrap(~feature, scales = "free_x", ncol = 2) +
    ggplot2::labs(title = "GAM smooth effects", x = "Predictor value (original units)",
      y = x$score_scale,
      caption = paste0(round(x$level * 100, 1), "% pointwise normal-approximation bands, conditional on fitted smoothing parameters.\nRugs: all captured training values. Centered terms exclude the intercept; not probabilities or causal effects.")) +
    sb_viz_theme()
}

# ---- mgcv-gam-adapters ----


sb_mgcv_gam_binary_learner_ids <- function() character()

sb_mgcv_gam_survival_learner_ids <- function() character()

sb_mgcv_gam_learner_ids <- function() {
  c(
    sb_mgcv_gam_binary_learner_ids(),
    sb_mgcv_gam_survival_learner_ids()
  )
}

sb_mgcv_gam_params <- function(learner) {
  if (!learner %in% sb_mgcv_gam_learner_ids()) {
    sb_abort(
      paste0("Unknown mgcv GAM learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  common <- list(
    basis = "cs",
    basis_dimension = 5L,
    method = "REML",
    threads = 1L
  )
  if (learner %in% sb_mgcv_gam_binary_learner_ids()) {
    return(c(common, list(family = "binomial", link = "logit")))
  }
  c(common, list(family = "Cox PH", link = "identity", ties = "Peto"))
}

sb_mgcv_available <- function() {
  requireNamespace("mgcv", quietly = TRUE)
}

sb_require_mgcv <- function() {
  if (!sb_mgcv_available()) {
    sb_abort(
      "Package `mgcv` is required for GAM learners.",
      "stormbreaker_missing_mgcv"
    )
  }
  invisible(TRUE)
}

sb_mgcv_internal_feature_names <- function(n_features) {
  paste0(".stormbreaker_mgcv_x", seq_len(n_features))
}

sb_mgcv_formula_data <- function(x, response, params) {
  engine_features <- sb_mgcv_internal_feature_names(ncol(x))
  response_name <- ".stormbreaker_mgcv_response"
  analysis <- x
  names(analysis) <- engine_features
  analysis[[response_name]] <- response
  smooth_terms <- paste0(
    "s(", engine_features,
    ", bs = '", params$basis,
    "', k = ", params$basis_dimension, ")"
  )
  formula <- stats::as.formula(
    paste(response_name, "~", paste(smooth_terms, collapse = " + ")),
    env = asNamespace("mgcv")
  )
  list(
    data = analysis,
    formula = formula,
    engine_features = engine_features,
    response_name = response_name
  )
}

sb_mgcv_fit_call <- function(formula_data, endpoint, learner, params) {
  sb_reject_retired_learner(learner)
}


sb_validate_mgcv_engine <- function(engine, learner, engine_features,
                                    require_baseline = TRUE) {
  params <- sb_mgcv_gam_params(learner)
  coefficients <- tryCatch(
    stats::coef(engine),
    error = function(error) NULL
  )
  smooth_terms <- if (inherits(engine, "gam") && is.list(engine$smooth)) {
    vapply(engine$smooth, function(smooth) {
      as.character(smooth$term)[[1L]]
    }, "")
  } else {
    character()
  }
  smooth_dimensions <- if (inherits(engine, "gam") &&
      is.list(engine$smooth)) {
    vapply(engine$smooth, function(smooth) {
      as.integer(smooth$bs.dim)[[1L]]
    }, integer(1))
  } else {
    integer()
  }
  smooth_ok <- inherits(engine, "gam") &&
    identical(unname(smooth_terms), engine_features) &&
    identical(
      unname(smooth_dimensions),
      rep(params$basis_dimension, length(engine_features))
    )
  convergence_ok <- if (learner %in% sb_mgcv_gam_binary_learner_ids()) {
    isTRUE(engine$converged)
  } else {
    is.null(engine$converged) || isTRUE(engine$converged)
  }
  common_ok <- smooth_ok && convergence_ok &&
    identical(engine$method, params$method) &&
    is.list(engine$outer.info) &&
    identical(engine$outer.info$conv, "full convergence") &&
    is.numeric(coefficients) && length(coefficients) > 0L &&
    !anyNA(coefficients) && all(is.finite(coefficients)) &&
    is.numeric(engine$sp) && length(engine$sp) == length(engine_features) &&
    !anyNA(engine$sp) && all(is.finite(engine$sp)) && all(engine$sp >= 0)
  family_ok <- identical(engine$family$family, params$family) &&
    identical(engine$family$link, params$link)
  baseline_ok <- TRUE
  if (learner %in% sb_mgcv_gam_survival_learner_ids() && require_baseline) {
    baseline <- engine$family$data
    baseline_ok <- is.list(baseline) &&
      is.numeric(baseline$tr) && is.numeric(baseline$h) &&
      length(baseline$tr) > 0L &&
      length(baseline$tr) == length(baseline$h) &&
      !anyNA(baseline$tr) && all(is.finite(baseline$tr)) &&
      all(baseline$tr > 0) &&
      !anyNA(baseline$h) && all(is.finite(baseline$h)) &&
      all(baseline$h >= 0)
  }
  if (!isTRUE(common_ok && family_ok && baseline_ok)) {
    sb_abort(
      paste0("Learner `", learner, "` returned an incomplete mgcv GAM."),
      "stormbreaker_incomplete_mgcv_fit"
    )
  }
  engine
}

sb_mgcv_prediction_engine <- function(engine, learner,
                                      retain_cox_baseline = TRUE) {
  patient_fields <- c(
    "model", "y", "z", "weights", "prior.weights", "working.weights",
    "fitted.values", "residuals", "linear.predictors", "offset", "hat",
    "call"
  )
  engine[intersect(patient_fields, names(engine))] <- NULL
  if (learner %in% sb_mgcv_gam_survival_learner_ids()) {
    if (isTRUE(retain_cox_baseline)) {
      engine$family$data <- engine$family$data[c("tr", "h")]
    } else {
      engine$family$data <- NULL
    }
  }
  engine
}

sb_new_mgcv_gam_model <- function(learner, features, engine_features,
                                  response_name, params, engine) {
  structure(
    list(
      learner = learner,
      features = features,
      engine_features = engine_features,
      response_name = response_name,
      params = params,
      engine = engine
    ),
    class = "stormbreaker_mgcv_gam"
  )
}

sb_validate_mgcv_binary_inputs <- function(x, outcome, learner) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      paste0("A `", learner, "` analysis fold must contain both classes."),
      "stormbreaker_mgcv_binary_single_class"
    )
  }
  list(x = x, outcome = outcome)
}

sb_validate_mgcv_survival_inputs <- function(x, time, event, learner) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time <= 0)) {
    sb_abort(
      paste0("Learner `", learner, "` requires strictly positive times."),
      "stormbreaker_mgcv_nonpositive_time"
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
      "stormbreaker_mgcv_no_event"
    )
  }
  list(x = x, time = time, event = as.integer(event))
}

sb_fit_mgcv_gam <- function(x, endpoint, learner, params) {
  sb_reject_retired_learner(learner)
}


sb_validate_mgcv_model_x <- function(model, x) {
  valid_model <- inherits(model, "stormbreaker_mgcv_gam") &&
    is.character(model$learner) && length(model$learner) == 1L &&
    model$learner %in% sb_mgcv_gam_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.character(model$engine_features) &&
    length(model$engine_features) == length(model$features) &&
    is.character(model$response_name) && length(model$response_name) == 1L &&
    is.list(model$params)
  if (!valid_model) {
    sb_abort(
      "The fitted mgcv GAM is incomplete or invalid.",
      "stormbreaker_bad_mgcv_model"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted mgcv GAM.",
      "stormbreaker_mgcv_schema_mismatch"
    )
  }
  x
}

sb_mgcv_predict_call <- function(engine, newdata, type, learner) {
  tryCatch(
    stats::predict(engine, newdata = newdata, type = type),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "Learner `", learner, "` could not predict from its mgcv GAM: ",
          conditionMessage(error)
        ),
        "stormbreaker_mgcv_prediction_error",
        parent = error
      )
    }
  )
}

sb_predict_mgcv_gam <- function(model, x) {
  sb_require_mgcv()
  x <- sb_validate_mgcv_model_x(model, x)
  sb_validate_mgcv_engine(
    model$engine,
    model$learner,
    model$engine_features,
    require_baseline = FALSE
  )
  newdata <- x
  names(newdata) <- model$engine_features
  type <- if (model$learner %in% sb_mgcv_gam_binary_learner_ids()) {
    "response"
  } else {
    "link"
  }
  values <- unname(sb_mgcv_predict_call(
    model$engine,
    newdata,
    type,
    model$learner
  ))
  valid <- length(values) == nrow(x) && !anyNA(values) &&
    all(is.finite(values))
  if (model$learner %in% sb_mgcv_gam_binary_learner_ids()) {
    valid <- valid && all(values >= 0 & values <= 1)
  }
  if (!valid) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid predictions."),
      if (model$learner %in% sb_mgcv_gam_binary_learner_ids()) {
        "stormbreaker_invalid_binary_prediction"
      } else {
        "stormbreaker_invalid_survival_prediction"
      }
    )
  }
  values
}

sb_mgcv_cox_baseline_hazard <- function(model, prediction_horizon) {
  baseline <- model$engine$family$data
  valid <- is.list(baseline) && is.numeric(baseline$tr) &&
    is.numeric(baseline$h) && length(baseline$tr) > 0L &&
    length(baseline$tr) == length(baseline$h) &&
    !anyNA(baseline$tr) && all(is.finite(baseline$tr)) &&
    !anyNA(baseline$h) && all(is.finite(baseline$h)) &&
    all(baseline$h >= 0)
  if (!valid) {
    sb_abort(
      "The fitted mgcv Cox baseline hazard is incomplete or invalid.",
      "stormbreaker_invalid_baseline_hazard"
    )
  }
  ordering <- order(baseline$tr, seq_along(baseline$tr))
  times <- baseline$tr[ordering]
  hazards <- baseline$h[ordering]
  index <- findInterval(prediction_horizon, times)
  cumulative_hazard <- if (index == 0L) 0 else hazards[[index]]
  if (!is.finite(cumulative_hazard) || cumulative_hazard < 0) {
    sb_abort(
      "The fitted mgcv Cox baseline hazard is incomplete or invalid.",
      "stormbreaker_invalid_baseline_hazard"
    )
  }
  unname(cumulative_hazard)
}

sb_predict_mgcv_gam_survival <- function(model, x_analysis, time, event,
                                         x_assessment,
                                         prediction_horizon) {
  if (!inherits(model, "stormbreaker_mgcv_gam") ||
      !model$learner %in% sb_mgcv_gam_survival_learner_ids()) {
    sb_abort(
      "The fitted mgcv Cox GAM is incomplete or invalid.",
      "stormbreaker_bad_mgcv_model"
    )
  }
  inputs <- sb_validate_mgcv_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted mgcv Cox GAM.",
      "stormbreaker_mgcv_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  risk <- sb_predict_mgcv_gam(model, x_assessment)
  baseline_hazard <- sb_mgcv_cox_baseline_hazard(model, prediction_horizon)
  sb_survival_prediction_result(
    prediction_horizon,
    baseline_hazard,
    source = "mgcv_cox_ph_native_baseline",
    ties = "Peto",
    risk = risk
  )
}

sb_extract_mgcv_gam_learner_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_mgcv_gam")) {
    sb_bad_learner_state("The fitted mgcv GAM learner state is invalid.")
  }
  baseline_hazard <- NULL
  if (model$learner %in% sb_mgcv_gam_survival_learner_ids() &&
      "prediction_horizon" %in% names(context)) {
    prediction <- sb_predict_mgcv_gam_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    baseline_hazard <- list(
      prediction_horizon = prediction$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      cumulative_baseline_hazard = prediction$cumulative_baseline_hazard,
      baseline_survival = prediction$baseline_survival,
      source = prediction$source,
      ties = prediction$ties
    )
  }
  engine <- sb_mgcv_prediction_engine(
    model$engine,
    model$learner,
    retain_cox_baseline = FALSE
  )
  structure(
    list(
      learner = model$learner,
      features = model$features,
      engine_features = model$engine_features,
      response_name = model$response_name,
      params = model$params,
      engine = engine,
      baseline_hazard = baseline_hazard
    ),
    class = c(
      "stormbreaker_mgcv_gam_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_mgcv_gam_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
      !learner %in% sb_mgcv_gam_learner_ids()) {
    sb_bad_learner_state("The saved mgcv GAM learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_new_mgcv_gam_model(
    learner,
    state$features,
    state$engine_features,
    state$response_name,
    state$params,
    state$engine
  )
  values <- sb_predict_mgcv_gam(model, x)
  if (learner %in% sb_mgcv_gam_binary_learner_ids()) {
    return(list(probability = values))
  }
  result <- list(risk = values)
  baseline <- state$baseline_hazard
  if (!is.null(baseline)) {
    valid_baseline <- is.list(baseline) &&
      is.numeric(baseline$prediction_horizon) &&
      length(baseline$prediction_horizon) == 1L &&
      !is.na(baseline$prediction_horizon) &&
      is.finite(baseline$prediction_horizon) &&
      baseline$prediction_horizon > 0 &&
      is.character(baseline$time_unit) && length(baseline$time_unit) == 1L &&
      !is.na(baseline$time_unit) && nzchar(baseline$time_unit)
    if (!valid_baseline) {
      sb_bad_learner_state(
        "The saved mgcv Cox baseline-hazard state is invalid."
      )
    }
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      values,
      baseline$cumulative_baseline_hazard
    )
  }
  result
}

sb_mgcv_gam_learner_adapter <- function(learner,
                                        params = sb_mgcv_gam_params(learner)) {
  if (!learner %in% sb_mgcv_gam_learner_ids()) {
    sb_abort(
      paste0("Unknown mgcv GAM learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  adapter <- list(
    id = learner,
    params = params,
    fit = if (learner %in% sb_mgcv_gam_binary_learner_ids()) {
      function(x, outcome) {
        sb_fit_mgcv_gam(x, list(outcome = outcome), learner, params)
      }
    } else {
      function(x, time, event) {
        sb_fit_mgcv_gam(
          x,
          list(time = time, event = event),
          learner,
          params
        )
      }
    },
    predict = sb_predict_mgcv_gam,
    extract_state = sb_extract_mgcv_gam_learner_state,
    predict_state = sb_predict_mgcv_gam_learner_state
  )
  if (learner %in% sb_mgcv_gam_survival_learner_ids()) {
    adapter$predict_survival <- sb_predict_mgcv_gam_survival
  }
  adapter
}

# ---- diagnostics-smooth ----

# ---- algorithm-bart ----

#' @rdname sb_gam_effects
#' @export
sb_bart_posterior <- function(fit, newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  id <- s$learner
  if (!id %in% c("bart_binary", "bart_survival", "bart_aft_survival")) sb_abort(
    "This posterior view requires BART.", "stormbreaker_unsupported_diagnostic")
  if (is.null(newdata)) sb_abort("Supply explicit queries for posterior BART predictions.",
    "stormbreaker_unavailable_diagnostic")
  x <- sb_select_final_feature_frame(newdata, s$features, "Query data", "query fitted features")
  if (!nrow(x)) sb_abort("Queries must contain at least one row.", "stormbreaker_invalid_diagnostic")
  e <- if (is.null(s$engine_raw)) s$engine else unserialize(s$engine_raw)
  n <- nrow(x)
  time <- NA_real_
  if (id == "bart_aft_survival") {
    draws <- sb_bart_aft_log_time_draws(sb_restore_rule_aft_model(s), x)
    scale <- "Latent log-time location"
  } else {
    matrix_x <- as.matrix(x)
    if (id == "bart_survival") {
      time <- as.numeric(e$times)
      matrix_x <- cbind(t = rep(time, n), matrix_x[rep(seq_len(n), each = length(time)), , drop = FALSE])
      if (length(e$rm.const)) matrix_x <- matrix_x[, e$rm.const, drop = FALSE]
    } else matrix_x <- sb_bart_prediction_matrix(e, matrix_x)
    fn <- get(if (id == "bart_binary") "predict.pbart" else "predict.survbart", asNamespace("BART"))
    invisible(utils::capture.output(p <- fn(e, matrix_x, mc.cores = 1L, openmp = FALSE)))
    draws <- if (id == "bart_binary") p$prob.test else p$surv.test
    scale <- if (id == "bart_binary") "Class-1 probability" else "Survival probability"
  }
  if (!is.matrix(draws) || ncol(draws) != n * length(time) || nrow(draws) < 2L ||
      any(!is.finite(draws))) sb_abort("Original posterior draws are incomplete.", "stormbreaker_invalid_diagnostic")
  interval <- apply(draws, 2L, stats::quantile, probs = c(.025, .975), names = FALSE, type = 7)
  tab <- data.frame(query = rep(seq_len(n), each = length(time)), time = rep(time, n),
    mean = colMeans(draws), lower = interval[1L, ], upper = interval[2L, ])
  structure(list(type = "bart_posterior", table = tab, draws = draws, score_scale = scale,
    learner = id, pipeline_id = fit$pipeline$id, source = "Original saved BART posterior tree draws; native prediction",
    interpretation = paste("Pointwise 95% equal-tail posterior credible limits, not simultaneous bands or validation.",
      "AFT values are latent log-time locations, not realized outcome prediction intervals.",
      "No new MCMC, training, surrogate model or patient prediction card.")), class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_bart_posterior <- function(x) {
  query <- time <- mean <- lower <- upper <- NULL
  if (x$learner == "bart_survival") {
    ribbon <- do.call(rbind, lapply(split(x$table, x$table$query), function(d) {
      z <- d[rep(seq_len(nrow(d)), each = 2L)[- (2L * nrow(d))], ]
      z$time <- c(d$time[1L], rep(d$time[-1L], each = 2L))
      z
    }))
    p <- ggplot2::ggplot(x$table, ggplot2::aes(time, mean)) +
      ggplot2::geom_ribbon(data = ribbon, ggplot2::aes(ymin = lower, ymax = upper), fill = "#2878A0", alpha = .18) +
      ggplot2::geom_step(colour = "#2878A0") + ggplot2::facet_wrap(~query, labeller = ggplot2::label_both) +
      ggplot2::labs(x = "Time")
  } else p <- ggplot2::ggplot(x$table, ggplot2::aes(factor(query), mean, ymin = lower, ymax = upper)) +
    ggplot2::geom_pointrange(colour = "#2878A0") + ggplot2::labs(x = "Query")
  p + ggplot2::labs(title = "BART posterior predictions", y = x$score_scale) + sb_viz_theme()
}

# ---- algorithm-earth ----

#' @rdname sb_gam_effects
#' @export
sb_earth_diagnostics <- function(fit, type = "basis_coefficients") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (fit$pipeline$learner != "earth_binary") sb_abort(
    "This diagnostic requires the earth binary learner.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("basis_coefficients", "pruning_gcv")) sb_abort(
    "Choose basis_coefficients or pruning_gcv.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(fit$learner_state$engine_raw)
  if (type == "basis_coefficients") {
    beta <- e$glm.coefficients
    if (!is.matrix(beta) || ncol(beta) != 1L || nrow(beta) != length(e$selected.terms) ||
        is.null(rownames(beta)) || any(!is.finite(beta))) sb_abort(
      "Original MARS logistic basis coefficients are invalid.", "stormbreaker_bad_learner_state")
    tab <- data.frame(term = rownames(beta), coefficient = as.numeric(beta),
      basis_index = e$selected.terms, intercept = e$selected.terms == 1L)
    interpretation <- paste("Original selected logistic basis coefficients, including the intercept; not the least-squares coefficients.",
      "Basis labels may round knot locations; exact cuts/directions are retained separately.",
      "Not raw-feature effects, SHAP or uncertainty.")
  } else {
    native <- fit$diagnostic_artifacts$learner$engine
    if (is.null(native$gcv.per.subset)) sb_abort("Original MARS pruning diagnostics were not captured.",
      "stormbreaker_missing_diagnostic_artifacts")
    gcv <- native$gcv.per.subset
    if (!is.numeric(gcv) || !length(gcv) || !any(is.finite(gcv))) sb_abort(
      "The native pruning sequence has no finite criterion values.", "stormbreaker_unavailable_diagnostic")
    tab <- data.frame(terms = seq_along(gcv), gcv = as.numeric(gcv),
      selected = seq_along(gcv) == length(e$selected.terms))
    interpretation <- paste("Native basis-pruning GCV; smaller is better. Term count includes the intercept.",
      "Orange point marks the original selected basis size, not a new optimum.",
      "This is not final logistic-model held-out deviance or a confidence interval; nonfinite criteria remain in the table and leave plot gaps.")
  }
  structure(list(type = paste0("earth_", type), table = tab,
    directions = e$dirs, cuts = e$cuts, selected_terms = e$selected.terms,
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    source = "original saved earth logistic basis and native pruning criterion", interpretation = interpretation),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_additive_shap_earth_terms <- function(state, x) {
  e <- unserialize(state$engine_raw)
  directions <- e$dirs[e$selected.terms, , drop = FALSE]
  if (any(rowSums(directions != 0) > 1L)) sb_abort_unsupported_additive_shap(state$learner)
  basis <- get("model.matrix.earth", asNamespace("earth"))(e, x)
  beta <- as.numeric(e$glm.coefficients)
  contributions <- matrix(0, nrow(x), ncol(x), dimnames = list(NULL, state$features))
  constant <- 0
  for (j in seq_len(ncol(basis))) {
    feature <- which(directions[j, ] != 0)
    if (length(feature)) contributions[, feature] <- contributions[, feature] + basis[, j] * beta[j] else
      constant <- constant + beta[j]
  }
  list(constant = constant, terms = contributions, score_scale = "class_1_log_odds")
}

sb_plot_earth_diagnostic <- function(x) {
  term <- coefficient <- terms <- gcv <- NULL
  tab <- x$table
  if (x$type == "earth_basis_coefficients") {
    tab$term <- factor(tab$term, levels = rev(tab$term))
    p <- ggplot2::ggplot(tab, ggplot2::aes(coefficient, term)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = .4) +
      ggplot2::geom_col(fill = "#2878A0", width = .65) +
      ggplot2::labs(title = "MARS logistic basis coefficients", x = "Logistic coefficient per basis unit", y = NULL)
  } else {
    tab$gcv[!is.finite(tab$gcv)] <- NA_real_
    p <- ggplot2::ggplot(tab, ggplot2::aes(terms, gcv)) +
      ggplot2::geom_line(colour = "#2878A0", na.rm = TRUE) +
      ggplot2::geom_point(colour = "#2878A0", na.rm = TRUE) +
      ggplot2::geom_point(data = tab[tab$selected, , drop = FALSE], colour = "#D77935", size = 2.5, na.rm = TRUE) +
      ggplot2::labs(title = "MARS basis pruning", x = "Number of basis terms", y = "GCV")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = if (x$type == "earth_basis_coefficients") max(5, 2 + .3 * nrow(tab)) else 5)
  p
}

# ---- algorithm-neural ----

sb_neural_parameters <- function(state) {
  if (!state$learner %in% c("nnet_binary", "neuralnet_multilayer_binary")) sb_abort(
    "This diagnostic requires a supported neural network.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(state$engine_raw)
  m <- e$model
  if (state$learner == "nnet_binary") {
    sizes <- m$n
    weights <- list()
    offset <- c(0L, cumsum(sizes))
    for (layer in 1:2) {
      w <- matrix(0, sizes[layer] + 1L, sizes[layer + 1L])
      for (j in seq_len(ncol(w))) {
        destination <- offset[layer + 1L] + j
        at <- seq.int(m$nconn[destination + 1L] + 1L, m$nconn[destination + 2L])
        sources <- m$conn[at]
        row <- ifelse(sources == 0, 1L, sources - offset[layer] + 1L)
        if (length(at) != nrow(w) || anyDuplicated(row) || any(row < 1L | row > nrow(w))) sb_abort(
          "Saved nnet connections do not match the supported feed-forward topology.", "stormbreaker_invalid_diagnostic")
        w[row, j] <- m$wts[at]
      }
      weights[[layer]] <- w
    }
  } else {
    weights <- m$weights[[1L]]
    sizes <- c(nrow(weights[[1L]]) - 1L, vapply(weights, ncol, integer(1)))
  }
  if (sizes[1L] != length(state$features) || utils::tail(sizes, 1L) != 1L ||
      any(!is.finite(unlist(weights))) || any(vapply(seq_along(weights), function(i)
        nrow(weights[[i]]) != sizes[i] + 1L, logical(1)))) sb_abort(
    "Saved neural weight dimensions or values are invalid.", "stormbreaker_bad_learner_state")
  list(weights = weights, sizes = sizes, center = e$center, scale = e$scale,
    activation = if (state$learner == "nnet_binary") "nnet_saturated_logistic" else "logistic",
    saturation = if (state$learner == "nnet_binary") c(lower = -15, upper = 15) else NULL)
}

#' @rdname sb_gam_effects
#' @export
sb_neural_diagnostics <- function(fit, type = "topology") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("topology", "weights")) sb_abort(
    "Choose topology or weights.", "stormbreaker_unsupported_diagnostic")
  p <- sb_neural_parameters(fit$learner_state)
  layers <- length(p$weights)
  labels <- c(list(fit$learner_state$features), lapply(seq_len(layers), function(i)
    if (i == layers) "Class 1" else paste0("H", i, ".", seq_len(p$sizes[i + 1L]))))
  nodes <- do.call(rbind, lapply(seq_along(labels), function(i) {
    n <- length(labels[[i]])
    data.frame(id = paste0("L", i - 1L, "_", seq_len(n)), label = labels[[i]], layer = i - 1L,
      position = rev(seq_len(n)) - (n + 1) / 2, bias = FALSE)
  }))
  bias <- data.frame(id = paste0("B", seq_len(layers)), label = "Bias", layer = seq_len(layers) - 1L,
    position = max(p$sizes) / 2 + 1, bias = TRUE)
  nodes <- rbind(nodes, bias)
  edges <- do.call(rbind, lapply(seq_len(layers), function(i) {
    w <- p$weights[[i]]
    data.frame(stage = paste0("Layer ", i - 1L, " -> ", i),
      source = rep(c(paste0("B", i), paste0("L", i - 1L, "_", seq_len(nrow(w) - 1L))), ncol(w)),
      target = rep(paste0("L", i, "_", seq_len(ncol(w))), each = nrow(w)),
      source_label = rep(c("Bias", labels[[i]]), ncol(w)),
      target_label = rep(labels[[i + 1L]], each = nrow(w)), weight = as.vector(w))
  }))
  structure(list(type = paste0("neural_", type), table = edges, nodes = nodes, parameters = p,
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    source = "original saved neural-network weights, connections and input scaling",
    interpretation = paste("Original feed-forward logistic network, retaining every connection and bias.",
      "Input coordinates use the fitted training centering/scaling; weights are not feature importance, log odds or SHAP.",
      "Repeated bias nodes represent constant one. No optimization history or new fit.",
      if (p$activation == "nnet_saturated_logistic") "Native nnet activations saturate to 0 below -15 and 1 above 15." else
        "Original neuralnet logistic activation retained.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_neural_diagnostic <- function(x) {
  layer <- position <- label <- xend <- yend <- source_label <- target_label <- weight <- NULL
  start <- end <- NULL
  nodes <- x$nodes
  tab <- x$table
  if (x$type == "neural_topology") {
    from <- match(tab$source, nodes$id)
    to <- match(tab$target, nodes$id)
    tab$layer <- nodes$layer[from]
    tab$position <- nodes$position[from]
    tab$xend <- nodes$layer[to]
    tab$yend <- nodes$position[to]
    tab$start <- tab$layer + .12
    tab$end <- tab$xend - .12
    p <- ggplot2::ggplot(tab, ggplot2::aes(layer, position)) +
      ggplot2::geom_segment(ggplot2::aes(xend = start, yend = position), colour = "#9AA6AF", alpha = .55, linewidth = .35) +
      ggplot2::geom_segment(ggplot2::aes(x = start, xend = end, yend = yend), colour = "#9AA6AF", alpha = .55, linewidth = .35) +
      ggplot2::geom_segment(ggplot2::aes(x = end, y = yend, xend = xend, yend = yend), colour = "#9AA6AF", alpha = .55, linewidth = .35) +
      ggplot2::geom_point(data = nodes, shape = 21, fill = "white", colour = "#2878A0", size = 4) +
      ggplot2::geom_label(data = nodes, ggplot2::aes(label = label), nudge_y = .3,
        size = 3, linewidth = 0, fill = "white", label.padding = grid::unit(.08, "lines")) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = .35)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = .6)) +
      ggplot2::labs(title = "Neural-network structure") + sb_viz_theme() +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  } else {
    tab$source_label <- factor(tab$source_label, levels = rev(unique(tab$source_label)))
    tab$target_label <- factor(tab$target_label, levels = unique(tab$target_label))
    p <- ggplot2::ggplot(tab, ggplot2::aes(target_label, source_label, fill = weight)) +
      ggplot2::geom_tile() + ggplot2::facet_wrap(~stage, scales = "free") +
      sb_neural_weight_scale() +
      ggplot2::labs(title = "Neural connection weights", x = "Receiving unit", y = "Sending unit", fill = "Weight") +
      sb_viz_theme() + ggplot2::theme(panel.grid = ggplot2::element_blank())
  }
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = max(6, 2 + .5 * max(x$parameters$sizes)))
  if (x$type != "neural_topology") attr(p, "stormbreaker_notes") <- c(attr(p, "stormbreaker_notes"),
    "Weight colors use a symmetric pseudo-log scale (sigma = 1); legend values are original weights. No clipping or rescaling of model coefficients.")
  p
}

sb_neural_weight_scale <- function() {
  ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#D77935",
    transform = scales::pseudo_log_trans(sigma = 1),
    breaks = function(limits) {
      maximum <- max(abs(limits))
      if (maximum == 0) return(0)
      if (maximum <= 10) return(pretty(limits, n = 5))
      powers <- 10^seq(max(0, floor(log10(maximum)) - 2), floor(log10(maximum)))
      c(-rev(powers), 0, powers)
    }, guide = ggplot2::guide_colourbar(display = "rectangles"))
}
