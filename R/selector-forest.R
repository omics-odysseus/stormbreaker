# ---- selector-forest ----

# ---- selector-forest-wrappers ----

sb_boruta_selector_ids <- function() {
  c("boruta_binary_roughfix", "boruta_survival_roughfix")
}

sb_varselrf_selector_ids <- function() {
  "varselrf_binary_oob_1se"
}

sb_forest_selector_ids <- function() {
  c(sb_boruta_selector_ids(), sb_varselrf_selector_ids())
}

sb_selector_budget_options <- function(id, options = NULL) {
  defaults <- sb_forest_selector_params(id)[c("num_trees", "maximum_runs", "seed")]
  if (is.null(options) || identical(options, list())) return(list())
  if (!is.list(options) || is.null(names(options)) || anyDuplicated(names(options)) ||
      anyNA(names(options)) || any(!names(options) %in% names(defaults))) {
    sb_abort(paste0("selector_options for ", id, " accepts only ",
      paste(names(defaults), collapse = ", "), "."), "stormbreaker_bad_selector_options")
  }
  for (name in names(options)) {
    options[[name]] <- sb_check_positive_integer(options[[name]], name,
      minimum = if (name == "maximum_runs") 11L else 1L)
  }
  # Canonical order and omission of explicit defaults preserve existing IDs.
  options <- options[intersect(names(defaults), names(options))]
  options[!vapply(names(options), function(n) identical(options[[n]], defaults[[n]]), logical(1))]
}

sb_forest_selector_params <- function(selector) {
  selector <- sub("_multiclass_", "_binary_", selector, fixed = TRUE)
  if (selector %in% sb_boruta_selector_ids()) {
    return(list(
      method = "Boruta",
      importance = "ranger_normalized_permutation",
      p_value = 0.01,
      multiple_comparison_adjustment = TRUE,
      maximum_runs = 100L,
      num_trees = 500L,
      num_threads = 1L,
      tentative_rule = "TentativeRoughFix",
      seed = 1L
    ))
  }
  if (selector %in% sb_varselrf_selector_ids()) {
    return(list(
      method = "varSelRF",
      importance = "randomForest_mean_decrease_accuracy",
      c_sd = 1,
      mtry_factor = 1,
      initial_trees = 5000L,
      iteration_trees = 2000L,
      drop_fraction = 0.2,
      whole_range = TRUE,
      recompute_importance = FALSE,
      seed = 1L
    ))
  }
  sb_abort(
    paste0("Unknown forest selector `", selector, "`."),
    "stormbreaker_unsupported_selector"
  )
}

sb_forest_selector_dependency <- function(selector) {
  selector <- sub("_multiclass_", "_binary_", selector, fixed = TRUE)
  if (selector %in% sb_boruta_selector_ids()) {
    return(list(
      package = "Boruta",
      engine = "ranger",
      available = requireNamespace("Boruta", quietly = TRUE) &&
        requireNamespace("ranger", quietly = TRUE)
    ))
  }
  if (selector %in% sb_varselrf_selector_ids()) {
    return(list(
      package = "varSelRF",
      engine = "randomForest",
      available = requireNamespace("varSelRF", quietly = TRUE) &&
        requireNamespace("randomForest", quietly = TRUE)
    ))
  }
  sb_abort(
    paste0("Unknown forest selector `", selector, "`."),
    "stormbreaker_unsupported_selector"
  )
}

sb_forest_selector_provenance <- function(selector) {
  packages <- if (selector %in% sb_boruta_selector_ids()) {
    c("Boruta", "ranger")
  } else if (selector %in% sb_varselrf_selector_ids()) {
    c("varSelRF", "randomForest")
  } else {
    sb_abort(
      paste0("Unknown forest selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  }
  fields <- paste0(tolower(packages), "_version")
  versions <- lapply(packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
      as.character(utils::packageVersion(package))
    } else {
      NA_character_
    }
  })
  stats::setNames(versions, fields)
}

sb_require_forest_selector_dependency <- function(selector) {
  dependency <- sb_forest_selector_dependency(selector)
  if (!dependency$available) {
    sb_abort(
      paste0(
        "Packages `", dependency$package, "` and `", dependency$engine,
        "` are required for selector `", selector, "`."
      ),
      "stormbreaker_missing_forest_selector_dependency"
    )
  }
  invisible(dependency)
}

sb_forest_selector_stop_on_warning <- function(code, selector) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      sb_abort(
        paste0(
          "Selector `", selector, "` raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_forest_selector_warning",
        parent = warning
      )
    }
  )
}

sb_validate_forest_selector_inputs <- function(x, endpoint, selector) {
  if (selector %in% c("boruta_multiclass_roughfix", "varselrf_multiclass_oob_1se")) {
    return(sb_multiclass_selector_inputs(x, endpoint$outcome))
  }
  if (identical(selector, "boruta_survival_roughfix")) {
    inputs <- sb_validate_classical_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      selector
    )
    return(list(
      x = inputs$x,
      endpoint = list(time = inputs$time, event = inputs$event)
    ))
  }
  inputs <- sb_validate_classical_binary_inputs(
    x,
    endpoint$outcome,
    selector
  )
  list(x = inputs$x, endpoint = list(outcome = inputs$outcome))
}

sb_new_forest_selector_fit <- function(selector, source_features,
                                        selected_features, params, scores,
                                        selection_history = NULL,
                                        audit = list()) {
  valid <- is.character(source_features) && length(source_features) > 0L &&
    !anyNA(source_features) && !anyDuplicated(source_features) &&
    is.character(selected_features) &&
    (length(selected_features) > 0L || identical(selector, "boruta_multiclass_roughfix")) &&
    !anyNA(selected_features) && !anyDuplicated(selected_features) &&
    all(selected_features %in% source_features) && is.data.frame(scores) &&
    identical(scores$feature, source_features) &&
    identical(scores$selected, source_features %in% selected_features)
  if (!valid) {
    sb_abort(
      "The native forest selector did not return a valid feature subset.",
      "stormbreaker_invalid_forest_selector"
    )
  }
  structure(
    list(
      id = selector,
      source_features = source_features,
      selected_features = selected_features,
      params = params,
      scores = scores,
      selection_history = selection_history,
      audit = audit
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_call_boruta <- function(x, response, params) {
  Boruta::Boruta(
    x,
    response,
    getImp = Boruta::getImpRfZ,
    pValue = params$p_value,
    mcAdj = params$multiple_comparison_adjustment,
    maxRuns = params$maximum_runs,
    doTrace = 0L,
    holdHistory = TRUE,
    num.trees = params$num_trees,
    num.threads = params$num_threads
  )
}

sb_fit_boruta_selector <- function(x, endpoint, selector, capture_history = FALSE,
                                  options = list()) {
  sb_require_forest_selector_dependency(selector)
  params <- utils::modifyList(sb_forest_selector_params(selector), options)
  inputs <- sb_validate_forest_selector_inputs(x, endpoint, selector)
  source_features <- names(inputs$x)
  response <- if (identical(selector, "boruta_survival_roughfix")) {
    survival::Surv(inputs$endpoint$time, inputs$endpoint$event)
  } else if (identical(selector, "boruta_multiclass_roughfix")) {
    inputs$endpoint$outcome
  } else {
    factor(inputs$endpoint$outcome, levels = c(0L, 1L))
  }

  fitted <- sb_forest_selector_stop_on_warning(
    withr::with_seed(
      params$seed,
      sb_call_boruta(inputs$x, response, params)
    ),
    selector
  )
  runs <- nrow(fitted$ImpHistory)
  initial_decision <- fitted$finalDecision
  rough_fixed <- any(fitted$finalDecision == "Tentative")
  if (rough_fixed) {
    fitted <- sb_forest_selector_stop_on_warning(
      Boruta::TentativeRoughFix(fitted),
      selector
    )
  }
  statistics <- Boruta::attStats(fitted)
  valid_statistics <- is.data.frame(statistics) &&
    identical(rownames(statistics), source_features) &&
    all(c(
      "meanImp", "medianImp", "minImp", "maxImp", "normHits", "decision"
    ) %in% names(statistics)) &&
    !anyNA(statistics$decision) &&
    all(as.character(statistics$decision) %in% c("Confirmed", "Rejected"))
  if (!valid_statistics) {
    sb_abort(
      "Boruta did not return complete rough-fixed feature decisions.",
      "stormbreaker_invalid_boruta_selector"
    )
  }
  native_selected <- Boruta::getSelectedAttributes(
    fitted,
    withTentative = FALSE
  )
  selected_features <- source_features[source_features %in% native_selected]
  if (length(selected_features) < 1L && !identical(selector, "boruta_multiclass_roughfix")) {
    sb_abort(
      "Boruta did not confirm any predictor after its rough fix.",
      "stormbreaker_boruta_no_selected_features"
    )
  }
  scores <- data.frame(
    feature = source_features,
    mean_importance = as.numeric(statistics$meanImp),
    median_importance = as.numeric(statistics$medianImp),
    minimum_importance = as.numeric(statistics$minImp),
    maximum_importance = as.numeric(statistics$maxImp),
    normalized_hits = as.numeric(statistics$normHits),
    decision = as.character(statistics$decision),
    selected = source_features %in% selected_features,
    stringsAsFactors = FALSE
  )
  if (anyNA(scores[-c(1L, 7L, 8L)]) ||
    any(!is.finite(as.matrix(scores[c(
      "mean_importance", "median_importance", "minimum_importance",
      "maximum_importance", "normalized_hits"
    )])))) {
    sb_abort(
      "Boruta returned non-finite feature-importance statistics.",
      "stormbreaker_invalid_boruta_selector"
    )
  }
  decision_counts <- table(factor(
    scores$decision,
    levels = c("Confirmed", "Rejected")
  ))
  result <- sb_new_forest_selector_fit(
    selector,
    source_features,
    selected_features,
    params,
    scores,
    audit = list(
      completed_runs = as.integer(runs),
      tentative_rough_fix_applied = rough_fixed,
      confirmed = as.integer(decision_counts[["Confirmed"]]),
      rejected = as.integer(decision_counts[["Rejected"]])
    )
  )
  if (capture_history) result$diagnostic_history <- list(
    importance = fitted$ImpHistory,
    decision_before_rough_fix = initial_decision,
    decision_after_rough_fix = fitted$finalDecision)
  result
}

sb_fit_varselrf_selector <- function(x, outcome, selector) {
  sb_require_forest_selector_dependency(selector)
  params <- sb_forest_selector_params(selector)
  inputs <- sb_validate_forest_selector_inputs(
    x,
    list(outcome = outcome),
    selector
  )
  source_features <- names(inputs$x)
  if (length(source_features) < 2L) {
    sb_abort(
      "varSelRF requires at least two predictors.",
      "stormbreaker_varselrf_too_few_predictors"
    )
  }
  fitted <- sb_forest_selector_stop_on_warning(
    withr::with_seed(params$seed, varSelRF::varSelRF(
      xdata = inputs$x,
      Class = if (identical(selector, "varselrf_multiclass_oob_1se"))
        inputs$endpoint$outcome else factor(inputs$endpoint$outcome, levels = c(0L, 1L)),
      c.sd = params$c_sd,
      mtryFactor = params$mtry_factor,
      ntree = params$initial_trees,
      ntreeIterat = params$iteration_trees,
      vars.drop.num = NULL,
      vars.drop.frac = params$drop_fraction,
      whole.range = params$whole_range,
      recompute.var.imp = params$recompute_importance,
      verbose = FALSE,
      returnFirstForest = FALSE,
      keep.forest = FALSE
    )),
    selector
  )
  selected_features <- fitted$selected.vars
  valid_selected <- is.character(selected_features) &&
    length(selected_features) >= 2L && !anyNA(selected_features) &&
    !anyDuplicated(selected_features) &&
    all(selected_features %in% source_features)
  importance <- as.matrix(fitted$initialImportances)
  valid_importance <- ncol(importance) == 1L &&
    identical(rownames(importance), source_features) &&
    !anyNA(importance) && all(is.finite(importance))
  history <- fitted$selec.history
  required_history <- c("Number.Variables", "OOB", "sd.OOB")
  valid_history <- is.data.frame(history) && nrow(history) > 0L &&
    all(required_history %in% names(history)) &&
    all(is.finite(as.matrix(history[required_history])))
  if (!valid_selected || !valid_importance || !valid_history) {
    sb_abort(
      "varSelRF did not return a valid OOB-selected feature subset.",
      "stormbreaker_invalid_varselrf_selector"
    )
  }
  scores <- data.frame(
    feature = source_features,
    mean_decrease_accuracy = as.numeric(importance[, 1L]),
    selected = source_features %in% selected_features,
    stringsAsFactors = FALSE
  )
  selection_history <- data.frame(
    number_of_variables = as.integer(history$Number.Variables),
    oob_error = as.numeric(history$OOB),
    oob_standard_error = as.numeric(history$sd.OOB),
    stringsAsFactors = FALSE
  )
  sb_new_forest_selector_fit(
    selector,
    source_features,
    selected_features,
    params,
    scores,
    selection_history = selection_history,
    audit = list(
      selection_rule = "smallest_subset_within_one_sd_of_minimum_oob_error",
      best_model_n_features = as.integer(fitted$best.model.nvars)
    )
  )
}

sb_apply_forest_selector <- function(fitted, x) {
  valid <- inherits(fitted, "stormbreaker_selector_fit") &&
    fitted$id %in% sb_forest_selector_ids() &&
    is.character(fitted$source_features) &&
    is.character(fitted$selected_features) &&
    length(fitted$selected_features) > 0L &&
    !anyNA(fitted$selected_features) &&
    !anyDuplicated(fitted$selected_features) &&
    all(fitted$selected_features %in% fitted$source_features)
  if (!valid) {
    sb_abort(
      "The fitted selector is not a valid forest-selector state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      "Predictor columns do not match the fitted forest selector schema.",
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_forest_selector_adapter <- function(selector, capture_history = FALSE, options = list()) {
  if (is.null(options)) options <- list()
  if (!selector %in% sb_forest_selector_ids()) {
    sb_abort(
      paste0("Unknown forest selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  }
  fit <- if (identical(selector, "boruta_survival_roughfix")) {
    function(x, time, event) {
      sb_fit_boruta_selector(x, list(time = time, event = event), selector, capture_history, options)
    }
  } else if (identical(selector, "boruta_binary_roughfix")) {
    function(x, outcome) {
      sb_fit_boruta_selector(x, list(outcome = outcome), selector, capture_history, options)
    }
  } else {
    function(x, outcome) {
      sb_fit_varselrf_selector(x, outcome, selector)
    }
  }
  list(
    id = selector,
    params = c(options, if (capture_history) list(capture_history = TRUE) else list()),
    needs_outcome = TRUE,
    fit = fit,
    apply = sb_apply_forest_selector
  )
}

# ---- selector-survival-weights ----

sb_survival_weight_selector_ids <- function() c("ridge_cox_weights", "plsrcox_weights", "svm_survival_weights")

sb_survival_weight_selector_params <- function(id) {
  model <- switch(id,
    ridge_cox_weights = list(alpha = 0, inner_folds = 5L, seed = 1L, lambda_rule = "lambda.min"),
    plsrcox_weights = sb_plsrcox_params(),
    svm_survival_weights = sb_survivalsvm_params())
  list(model = model, fraction = .5, minimum = 2L, rule = "absolute_coefficient_times_training_sd")
}

sb_survival_weight_selector_adapter <- function(id, params, capture_history = FALSE) {
  list(id = id, params = c(params, if (capture_history) list(capture_history = TRUE)), needs_outcome = TRUE,
    fit = function(x, time, event, group_ids = NULL) {
      x <- sb_validate_predictor_frame(x, "x")
      evidence <- list()
      if (id == "ridge_cox_weights") {
        tuned <- sb_fit_lasso_nonzero(x, time, event, params$model$inner_folds,
          params$model$seed, group_ids, alpha = 0, id = id, capture_history = capture_history)
        beta <- tuned$scores$coefficient[match(names(x), tuned$scores$feature)]
        evidence <- list(lambda = tuned$lambda, tuning_curve = tuned$tuning_curve,
          fold_id = tuned$fold_id, group_ids = group_ids, coefficient_path = tuned$coefficient_path)
      } else if (id == "plsrcox_weights") {
        controls <- params$model
        controls$components <- min(controls$components, ncol(x))
        model <- sb_fit_plsrcox(x, time, event, controls)
        beta <- model$coefficient
        evidence$components <- controls$components
        if (capture_history) evidence$loadings <- model$engine$pp
      } else {
        model <- sb_fit_survivalsvm(x, time, event, params$model)
        beta <- model$coefficient
      }
      sd <- vapply(x, stats::sd, numeric(1))
      score <- abs(beta * sd)
      if (length(score) != ncol(x) || any(!is.finite(score))) stop("Native screening weights are unavailable.")
      eligible <- which(score > 0)
      if (!length(eligible)) stop("No nonzero screening weights; no fallback genes were added.")
      n_keep <- min(length(eligible), max(params$minimum, ceiling(ncol(x) * params$fraction)))
      chosen <- eligible[order(-score[eligible], eligible)][seq_len(n_keep)]
      selected <- seq_len(ncol(x)) %in% chosen
      structure(list(id = id, params = params, source_features = names(x),
        selected_features = names(x)[selected],
        scores = data.frame(feature = names(x), coefficient = as.numeric(beta),
          training_sd = unname(sd), score = unname(score), selected = selected),
        evidence = evidence), class = "stormbreaker_selector_fit")
    },
    apply = function(fitted, x) {
      x <- sb_validate_predictor_frame(x, "x")
      if (!identical(names(x), fitted$source_features)) stop("Predictors differ from the fitted screening model.")
      x[fitted$selected_features]
    })
}

# ---- selector-reuse ----

.sb_selector_reuse_context <- new.env(parent = emptyenv())
.sb_selector_reuse_context$active <- NULL

#' @rdname sb_compare_cv
#' @export
sb_with_selector_reuse <- function(expr, enabled = TRUE, diagnostics = FALSE) {
  if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled) ||
      !is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) {
    sb_abort("enabled and diagnostics must be single non-missing logical values.",
      "stormbreaker_bad_selector_reuse")
  }
  previous <- .sb_selector_reuse_context$active
  on.exit(.sb_selector_reuse_context$active <- previous, add = TRUE)
  context <- new.env(parent = emptyenv())
  context$entries <- new.env(parent = emptyenv())
  context$hits <- 0L
  context$misses <- 0L
  context$n_entries <- 0L
  context$final_learners <- new.env(parent = emptyenv())
  context$final_learner_hits <- 0L
  context$final_learner_misses <- 0L
  context$cv_learners <- new.env(parent = emptyenv())
  context$cv_learner_hits <- 0L
  context$cv_learner_misses <- 0L
  context$learner_byte_limit <- 256 * 1024^2
  context$learner_bytes <- 0
  context$learner_clock <- 0
  context$learner_evictions <- 0L
  context$learner_oversized <- 0L
  context$preparations <- new.env(parent = emptyenv())
  context$preparation_hits <- 0L
  context$preparation_misses <- 0L
  .sb_selector_reuse_context$active <- if (enabled) context else NULL
  value <- eval(substitute(expr), envir = parent.frame())
  if (diagnostics) {
    return(list(result = value, reuse = list(enabled = enabled,
      hits = context$hits, misses = context$misses, entries = context$n_entries,
      final_learner_hits = context$final_learner_hits,
      final_learner_misses = context$final_learner_misses,
      final_learner_entries = sum(lengths(as.list(context$final_learners))),
      cv_learner_hits = context$cv_learner_hits,
      cv_learner_misses = context$cv_learner_misses,
      cv_learner_entries = sum(lengths(as.list(context$cv_learners))),
      learner_bytes = context$learner_bytes,
      learner_byte_limit = context$learner_byte_limit,
      learner_evictions = context$learner_evictions,
      learner_oversized = context$learner_oversized,
      preparation_hits = context$preparation_hits,
      preparation_misses = context$preparation_misses)))
  }
  value
}

# Internal handoff for bounded Boruta workers. The controller owns persistence;
# these objects contain analysis predictors/endpoints and are not portable models.
# Consumers still go through the ordinary exact-input/signature matcher.
sb_boruta_worker_checkpoint <- function(id) {
  context <- .sb_selector_reuse_context$active
  if (!id %in% c("boruta_binary_roughfix", "boruta_survival_roughfix") ||
      is.null(context)) {
    stop("Boruta worker handoff requires an active reuse context and Boruta selector.")
  }
  list(id = id, entries = context$entries[[id]])
}

sb_restore_boruta_worker_checkpoint <- function(checkpoint) {
  context <- .sb_selector_reuse_context$active
  if (is.null(context) || !checkpoint$id %in%
      c("boruta_binary_roughfix", "boruta_survival_roughfix")) {
    stop("Boruta worker handoff requires an active reuse context and Boruta selector.")
  }
  previous <- length(context$entries[[checkpoint$id]])
  context$entries[[checkpoint$id]] <- checkpoint$entries
  context$n_entries <- context$n_entries - previous + length(checkpoint$entries)
  invisible(NULL)
}

sb_fit_final_runtime_learner <- function(learner, x, endpoint, patient_ids,
                                       selector_id) {
  sb_fit_reusable_runtime_learner(learner, x, endpoint, patient_ids,
    selector_id, phase = "final")
}

sb_fit_reusable_runtime_learner <- function(learner, x, endpoint, patient_ids,
                                           selector_id, phase = "cv") {
  context <- .sb_selector_reuse_context$active
  eligible <- c(sb_survivalsvm_learner_ids(), "logistic", "coxph", "gaussian_naive_bayes",
    "glmnet_logistic", "glmnet_cox", sb_specialized_binary_linear_learner_ids(),
    sb_classical_binary_learner_ids(), sb_mboost_cox_learner_ids(),
    sb_mboost_aft_learner_ids(), sb_mboost_ranking_learner_ids(),
    sb_mgcv_gam_learner_ids(), sb_flexsurv_parametric_learner_ids(),
    sb_flexsurvspline_learner_ids(), sb_timereg_survival_learner_ids(),
    sb_ahaz_survival_learner_ids(), sb_coxphf_learner_ids(),
    sb_ncvreg_learner_ids(), sb_abess_learner_ids(),
    sb_omics_specialist_binary_learner_ids(), sb_mclust_da_learner_ids(),
    sb_mixomics_da_learner_ids(), sb_hdda_binary_learner_ids(),
    sb_ranger_learner_ids(), sb_survreg_learner_ids(),
    sb_xgboost_learner_ids(), sb_lightgbm_binary_learner_ids(),
    "catboost_binary", "catboost_cox", "catboost_aft_survival",
    "gbm_cox", "coxboost", "plsrcox", "superpc",
    sb_mtlr_survival_learner_ids(), sb_pamr_survival_learner_ids(),
    sb_akritas_learner_ids(), sb_bagged_knn_survival_learner_ids(),
    sb_bujar_linear_learner_ids(), sb_rlt_learner_ids(), sb_wsrf_binary_learner_ids(),
    "adaboost_binary", "logitboost_binary", "rpart_survival",
    "rsf_survival", "ctree_survival", "cforest_survival",
    "aorsf_survival", "slope_binary", "random_glm_binary",
    "cancerclass_binary", "bart_survival", "bart_aft_survival",
    "coxboost_cv", "plsrcox_cv", "superpc_cv", "gbm_cox_cv",
    "pre_binary_rule_ensemble", "pre_survival_rule_ensemble",
    "penaft_gehan_elastic_net")
  if (is.null(context) || !learner$id %in% eligible) {
    return(sb_fit_runtime_learner(learner, x, endpoint))
  }
  # These native fits consume caller RNG (including AdaBoost's rpart CV).
  # Matching its starting state and replaying the transition preserves both
  # the model/diagnostics and the random stream of subsequent computations.
  rng_sensitive <- sb_learner_uses_caller_rng(learner$id)
  rng_input <- if (rng_sensitive && exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", .GlobalEnv, inherits = FALSE)
  } else NULL
  signature <- list(id = learner$id, params = learner$params,
    fit_context = learner$fit_context,
    rng_input = rng_input,
    fit_body = body(learner$fit), fit_formals = formals(learner$fit),
    rng_kind = RNGkind())
  store <- context[[paste0(phase, "_learners")]]
  hits <- paste0(phase, "_learner_hits")
  misses <- paste0(phase, "_learner_misses")
  entries <- store[[learner$id]]
  for (index in seq_along(entries)) {
    entry <- entries[[index]]
    if (identical(signature, entry$signature) &&
      identical(x, entry$x) && identical(endpoint, entry$endpoint) &&
      identical(patient_ids, entry$patient_ids)) {
      context[[hits]] <- context[[hits]] + 1L
      context$learner_clock <- context$learner_clock + 1
      entries[[index]]$last_used <- context$learner_clock
      store[[learner$id]] <- entries
      if (rng_sensitive) assign(".Random.seed", entry$rng_after, .GlobalEnv)
      return(sb_restore_reused_learner_model(unserialize(entry$model), learner$id))
    }
  }
  context[[misses]] <- context[[misses]] + 1L
  warned <- FALSE
  model <- withCallingHandlers(sb_fit_runtime_learner(learner, x, endpoint),
    warning = function(w) warned <<- TRUE)
  if (!warned && length(attr(model, "stormbreaker_fit_warnings", exact = TRUE)) == 0L &&
      (!rng_sensitive || !is.null(rng_input))) {
    entry <- list(signature = signature,
      x = x, endpoint = endpoint, patient_ids = patient_ids,
      model = serialize(sb_store_reused_learner_model(model, learner$id), NULL),
      rng_after = if (rng_sensitive) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL)
    sb_retain_shared_learner(context, phase, learner$id, entry)
  }
  model
}

sb_retain_shared_learner <- function(context, phase, id, entry) {
  context$learner_clock <- context$learner_clock + 1
  entry$last_used <- context$learner_clock
  entry$bytes <- 0
  entry$bytes <- as.numeric(utils::object.size(entry))
  if (entry$bytes > context$learner_byte_limit) {
    context$learner_oversized <- context$learner_oversized + 1L
    return(invisible(NULL))
  }
  while (context$learner_bytes + entry$bytes > context$learner_byte_limit) {
    oldest <- Inf
    for (name in c("cv_learners", "final_learners")) {
      store <- context[[name]]
      for (key in ls(store, all.names = TRUE)) {
        entries <- store[[key]]
        for (i in seq_along(entries)) if (entries[[i]]$last_used < oldest) {
          oldest <- entries[[i]]$last_used
          victim <- list(store = name, key = key, index = i)
        }
      }
    }
    store <- context[[victim$store]]
    entries <- store[[victim$key]]
    context$learner_bytes <- context$learner_bytes - entries[[victim$index]]$bytes
    store[[victim$key]] <- entries[-victim$index]
    context$learner_evictions <- context$learner_evictions + 1L
  }
  store <- context[[paste0(phase, "_learners")]]
  store[[id]] <- c(store[[id]], list(entry))
  context$learner_bytes <- context$learner_bytes + entry$bytes
  invisible(NULL)
}

sb_store_reused_learner_model <- function(model, learner) {
  if (learner %in% sb_lightgbm_binary_learner_ids()) {
    text <- model$booster$save_model_to_string()
    fields <- intersect(c("best_iter", "best_score", "params", "record_evals",
      "data_processor", "raw"), names(model$booster))
    metadata <- stats::setNames(lapply(fields, function(name) model$booster[[name]]), fields)
    # LightGBM's text format rounds split gains, although predictions are exact.
    # Retain the original diagnostic values alongside the independent booster.
    model$reused_native_importance <- list(model_string = text,
      observed = as.data.frame(lightgbm::lgb.importance(model$booster, percentage = FALSE),
        stringsAsFactors = FALSE))
    model$booster <- NULL
    return(list(format = "lightgbm", value = model, text = text, metadata = metadata))
  }
  if (learner %in% sb_xgboost_learner_ids()) {
    raw <- xgboost::xgb.save.raw(model$booster, raw_format = "json")
    # Older XGBoost keeps R-side call/iteration metadata alongside the pointer.
    metadata <- if (is.list(model$booster)) {
      model$booster[setdiff(names(model$booster), "handle")]
    } else NULL
    model$booster <- NULL
    return(list(format = "xgboost", value = model, raw = raw, metadata = metadata))
  }
  if (learner %in% c(sb_mboost_cox_learner_ids(),
      sb_mboost_aft_learner_ids(), sb_mboost_ranking_learner_ids(),
      "gamboost_binary", "blackboost_binary", "catboost_binary",
      "catboost_cox", "catboost_aft_survival")) {
    # mstop<- mutates native environments: every consumer needs an isolated fit.
    return(list(format = "serialized", value = serialize(model, NULL)))
  }
  list(format = "value", value = model)
}

sb_restore_reused_learner_model <- function(stored, learner) {
  if (stored$format == "serialized") return(unserialize(stored$value))
  if (stored$format == "lightgbm") {
    model <- stored$value
    model$booster <- lightgbm::lgb.load(model_str = stored$text)
    for (name in names(stored$metadata)) model$booster[[name]] <- stored$metadata[[name]]
    return(model)
  }
  if (stored$format == "xgboost") {
    model <- stored$value
    model$booster <- sb_load_xgboost_state_model(list(model_raw = stored$raw,
      params = model$params, features = model$features), learner)$booster
    if (is.list(model$booster) && !is.null(stored$metadata)) {
      model$booster[names(stored$metadata)] <- stored$metadata
    }
    return(model)
  }
  stored$value
}

sb_fit_runtime_selector <- function(selector, x, endpoint) {
  context <- .sb_selector_reuse_context$active
  if (is.null(context) ||
      !selector$id %in% c(sb_forest_selector_ids(), "stepcox", "lasso_nonzero",
        "univ_cox_top_k", "univ_binary_auc_top_k",
        sb_high_dimensional_selector_ids(), sb_pamr_selector_ids(),
        sb_binary_model_selector_ids(), sb_abess_selector_ids(),
        "enet_cox_nonzero", sb_survival_weight_selector_ids(),
        "coxboost_nonzero", "rsf_minimal_depth", sb_survival_boosting_selector_ids())) {
    return(sb_fit_runtime_selector_uncached(selector, x, endpoint))
  }
  # Source-reference environments change identity after worker serialization;
  # executable body/formals, rather than debugger metadata, identify the fit.
  fit_definition <- utils::removeSource(selector$fit)
  signature <- list(id = selector$id, adapter_params = selector$params,
    capture_history = selector$capture_history,
    upstream_context = selector$upstream_context,
    method_params = if (selector$id %in% sb_forest_selector_ids()) {
      sb_forest_selector_params(selector$id)
    } else NULL,
    needs_outcome = selector$needs_outcome,
    fit_body = body(fit_definition), fit_formals = formals(fit_definition),
    rng_kind = RNGkind())
  entries <- context$entries[[selector$id]]
  for (entry in entries) {
    if (identical(signature, entry$signature) && identical(x, entry$x) &&
        identical(endpoint, entry$endpoint)) {
      context$hits <- context$hits + 1L
      return(entry$value)
    }
  }
  context$misses <- context$misses + 1L
  warned <- FALSE
  value <- withCallingHandlers(
    sb_fit_runtime_selector_uncached(selector, x, endpoint),
    warning = function(w) warned <<- TRUE
  )
  if (!warned) {
    entries[[length(entries) + 1L]] <- list(signature = signature,
      x = x, endpoint = endpoint, value = value)
    context$entries[[selector$id]] <- entries
    context$n_entries <- context$n_entries + 1L
  }
  value
}
