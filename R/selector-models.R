# ---- selector-binary-models ----

sb_binary_model_selector_ids <- function() {
  c("lasso_binary_nonzero", "stepglm_binary_aic",
    "glmboost_binary_nonzero", "rf_binary_importance",
    "enet_binary_nonzero", "gbm_binary_importance", "svm_rfe_binary", "gamboost_binary_nonzero", "xgboost_binary_importance", "plsrglm_binary_weights")
}

sb_binary_model_selector_options <- function(id, options = NULL) {
  defaults <- switch(id,
    plsrglm_binary_weights = list(),
    xgboost_binary_importance = list(seed = 1L),
    gamboost_binary_nonzero = list(seed = 1L),
    svm_rfe_binary = list(inner_folds = 5L, seed = 1L),
    lasso_binary_nonzero = list(inner_folds = 5L, seed = 1L, lambda_rule = "lambda.1se"),
    enet_binary_nonzero = list(inner_folds = 5L, seed = 1L, lambda_rule = "lambda.1se", alpha = .5),
    gbm_binary_importance = list(seed = 1L),
    stepglm_binary_aic = list(direction = "both"),
    glmboost_binary_nonzero = list(iterations = 100L, shrinkage = .1),
    rf_binary_importance = list(trees = 500L, seed = 1L, importance = "permutation", threshold = 0))
  if (is.null(defaults)) sb_abort("Unknown binary model selector.", "stormbreaker_unsupported_selector")
  optional <- if (id == "rf_binary_importance") c("top_n", "spearman_cutoff") else character()
  if (!is.null(options) && (!is.list(options) || (length(options) &&
      (is.null(names(options)) || anyNA(names(options)) || anyDuplicated(names(options)) ||
       !all(names(options) %in% c(names(defaults), optional)))))) sb_abort(
    "Unknown or duplicated selector_options; see sb_pipeline help.", "stormbreaker_bad_selector_options")
  defaults[intersect(names(defaults), names(options))] <-
    options[intersect(names(defaults), names(options))]
  for (n in intersect(names(defaults), c("inner_folds", "seed", "iterations", "trees"))) {
    defaults[[n]] <- sb_check_positive_integer(defaults[[n]], n, minimum = if (n == "inner_folds") 3L else 1L)
  }
  if (id %in% c("lasso_binary_nonzero", "enet_binary_nonzero")) defaults$lambda_rule <- sb_validate_glmnet_lambda_rule(defaults$lambda_rule)
  if (id == "enet_binary_nonzero" && (!is.numeric(defaults$alpha) || length(defaults$alpha) != 1L ||
      !is.finite(defaults$alpha) || defaults$alpha <= 0 || defaults$alpha > 1)) {
    sb_abort("Enet selection requires alpha in (0,1].", "stormbreaker_bad_selector_options")
  }
  if (id == "stepglm_binary_aic") defaults$direction <- sb_validate_stepglm_direction(defaults$direction)
  if (id == "glmboost_binary_nonzero" && (!is.numeric(defaults$shrinkage) ||
      length(defaults$shrinkage) != 1L || !is.finite(defaults$shrinkage) ||
      defaults$shrinkage <= 0 || defaults$shrinkage > 1)) sb_abort(
    "shrinkage must be in (0,1].", "stormbreaker_bad_selector_options")
  if (id == "rf_binary_importance") {
    if (!is.character(defaults$importance) || length(defaults$importance) != 1L ||
        is.na(defaults$importance) || !defaults$importance %in% c("permutation", "gini") ||
        !is.numeric(defaults$threshold) || length(defaults$threshold) != 1L ||
        !is.finite(defaults$threshold)) sb_abort(
      "RF needs permutation/gini importance and a finite threshold.", "stormbreaker_bad_selector_options")
    if (!is.null(options$top_n)) {
      if (defaults$threshold != 0) sb_abort(
        "Use top_n or a nondefault importance threshold, not both.",
        "stormbreaker_bad_selector_options")
      defaults$top_n <- sb_check_positive_integer(options$top_n, "top_n")
    }
    cutoff <- options$spearman_cutoff
    if (!is.null(cutoff)) {
      if (is.null(defaults$top_n) || !is.numeric(cutoff) || length(cutoff) != 1L ||
          !is.finite(cutoff) || cutoff <= 0 || cutoff > 1) {
        sb_abort("spearman_cutoff requires top_n and a value in (0, 1].",
          "stormbreaker_bad_selector_options")
      }
      defaults$spearman_cutoff <- as.numeric(cutoff)
    }
  }
  defaults
}

sb_rf_ranked_screen <- function(x, scores, top_n, spearman_cutoff = NULL) {
  top <- names(x)[utils::head(order(-scores, seq_along(scores)), top_n)]
  varying <- vapply(x[top], function(v) any(v != v[1L]), logical(1))
  retained <- top[varying]
  correlation <- NULL
  if (!is.null(spearman_cutoff) && length(retained)) {
    correlation <- stats::cor(x[retained], method = "spearman")
    kept <- integer()
    for (i in seq_along(retained)) {
      if (!any(abs(correlation[i, kept]) > spearman_cutoff)) kept <- c(kept, i)
    }
    retained <- retained[kept]
  }
  list(top_features = top, constant_features = top[!varying],
    correlation = correlation, retained_features = retained)
}

sb_fit_binary_model_selector <- function(x, outcome, id, params, group_ids = NULL, capture_history = FALSE) {
  if (id == "svm_rfe_binary") {
    fit <- sb_fit_svm_rfe(x, outcome, group_ids, params$inner_folds, params$seed)
    fit$params <- params
    return(fit)
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (length(outcome) != nrow(x) || anyNA(outcome) || !all(outcome %in% c(0, 1)) ||
      length(unique(outcome)) != 2L) sb_abort("Binary selection needs both 0/1 classes.", "stormbreaker_bad_binary_outcome")
  original_names <- names(x)
  work <- x
  names(work) <- paste0("feature_", seq_len(ncol(x)))
  evidence <- list()
  if (id %in% c("lasso_binary_nonzero", "enet_binary_nonzero")) {
    alpha <- if (id == "enet_binary_nonzero") params$alpha else 1
    tuned <- sb_select_glmnet_logistic_lambda(work, outcome,
      inner_folds = params$inner_folds, seed = params$seed,
      selection_rule = params$lambda_rule, group_ids = group_ids, alpha = alpha,
      capture_path = capture_history)
    model <- glmnet::glmnet(as.matrix(work), outcome, family = "binomial", alpha = alpha,
      lambda = tuned$lambda, standardize = TRUE)
    score <- as.numeric(stats::coef(model))[-1L]
    selected <- score != 0
    evidence <- list(lambda = tuned$lambda, tuning_curve = tuned$curve,
      fold_id = tuned$fold_id, group_ids = group_ids)
    if (capture_history) {
      evidence$coefficient_path <- tuned$coefficient_path
      rownames(evidence$coefficient_path$beta) <- original_names
    }
  } else if (id == "plsrglm_binary_weights") {
    controls <- sb_specialized_binary_linear_params("plsrglm_binary")
    controls$components <- min(controls$components, ncol(work))
    model <- sb_fit_specialized_binary_linear(work, outcome, "plsrglm_binary", controls)
    score <- abs(model$coefficient * vapply(work, stats::sd, numeric(1)))
    eligible <- which(score > 0 & is.finite(score))
    n_keep <- min(length(eligible), max(2L, ceiling(ncol(work) / 2)))
    chosen <- eligible[order(-score[eligible], eligible)][seq_len(n_keep)]
    selected <- seq_len(ncol(work)) %in% chosen
    evidence <- list(rule = "absolute_coefficient_times_training_sd_top_half_minimum_two",
      components = model$effective$components)
    if (capture_history) evidence$loadings <- model$diagnostic_components$loadings
  } else if (id == "xgboost_binary_importance") {
    native_params <- sb_xgboost_params("xgboost_binary")
    native_params$seed <- params$seed
    model <- sb_fit_xgboost(work, outcome, "xgboost_binary", native_params)
    score <- sb_xgboost_selector_gain(model)
    selected <- score > 0
    evidence <- list(rule = "positive_native_gain", xgboost_version = as.character(utils::packageVersion("xgboost")))
  } else if (id == "gamboost_binary_nonzero") {
    native_params <- sb_classical_binary_params("gamboost_binary")
    native_params$seed <- params$seed
    model <- sb_fit_classical_binary(work, outcome, "gamboost_binary", native_params)$engine
    score <- sb_mboost_spline_selector_scores(model, names(work))
    selected <- score != 0
    loss <- as.numeric(mboost::risk(model))
    evidence <- list(rule = "active_bbs_spline_terms", score_definition = "spline_coefficient_L2_norm",
      training_loss = data.frame(step = seq_along(loss) - 1L, loss = loss))
  } else if (id == "gbm_binary_importance") {
    if (!requireNamespace("gbm", quietly = TRUE)) sb_abort("Install gbm for GBM selection.", "stormbreaker_missing_dependency")
    model <- withr::with_seed(params$seed, gbm::gbm.fit(x = as.matrix(work), y = outcome, distribution = "bernoulli",
      n.trees = 100L, interaction.depth = 1L, n.minobsinnode = 10L,
      shrinkage = .1, bag.fraction = 1, nTrain = nrow(work), keep.data = FALSE,
      verbose = FALSE, var.names = names(work)))
    importance <- summary(model, n.trees = 100L, plotit = FALSE)
    score <- importance$rel.inf[match(names(work), importance$var)]
    selected <- score > 0
    evidence <- list(rule = "positive_relative_influence", trees = 100L,
      training_loss = data.frame(trees = seq_along(model$train.error), loss = model$train.error))
  } else if (id %in% c("stepglm_binary_aic", "glmboost_binary_nonzero")) {
    if (id == "stepglm_binary_aic") {
      native <- sb_fit_stepglm_binary(list(x = work, outcome = outcome),
        list(direction = params$direction, k = 2, max_steps = 1000L), capture_history)
    } else {
      if (!requireNamespace("mboost", quietly = TRUE)) sb_abort("Install mboost for glmBoost selection.", "stormbreaker_missing_dependency")
      native <- sb_fit_glmboost_binary(list(x = work, outcome = outcome),
        list(iterations = params$iterations, shrinkage = params$shrinkage, center = FALSE), capture_history)
    }
    score <- as.numeric(native$coefficient[names(work)])
    selected <- if (id == "stepglm_binary_aic") names(work) %in% native$effective$selected_features else score != 0
    evidence <- native$effective
    evidence$selected_features <- original_names[match(evidence$selected_features, names(work))]
  } else {
    if (!requireNamespace("randomForest", quietly = TRUE)) sb_abort("Install randomForest for RF selection.", "stormbreaker_missing_dependency")
    model <- withr::with_seed(params$seed, randomForest::randomForest(x = work, y = factor(outcome, levels = c(0, 1)),
      ntree = params$trees, importance = TRUE, keep.forest = FALSE))
    importance <- randomForest::importance(model, type = if (params$importance == "permutation") 1 else 2, scale = FALSE)
    score <- as.numeric(importance[names(work), 1L])
    selected <- score > params$threshold
    if (capture_history) evidence$oob_error <- data.frame(
      trees = seq_len(nrow(model$err.rate)), model$err.rate, check.names = FALSE)
  }
  if (length(score) != ncol(x) || any(!is.finite(score))) sb_abort(
    "Selector returned invalid feature scores.", "stormbreaker_invalid_selector_scores")
  if (identical(id, "rf_binary_importance") && !is.null(params$top_n)) {
    evidence$rf_screen <- sb_rf_ranked_screen(x, score, params$top_n,
      params$spearman_cutoff)
    selected <- original_names %in% evidence$rf_screen$retained_features
  }
  if (!any(selected)) sb_abort("No features met the declared selection rule; no fallback genes were added.", "stormbreaker_empty_selection")
  structure(list(id = id, source_features = original_names,
    selected_features = original_names[selected], params = params,
    scores = data.frame(feature = original_names, score = score, selected = selected),
    evidence = evidence), class = "stormbreaker_selector_fit")
}

sb_binary_model_selector_adapter <- function(id, params, capture_history = FALSE) {
  list(id = id, params = params, capture_history = capture_history, needs_outcome = TRUE,
    fit = function(x, outcome, group_ids = NULL) {
      if (capture_history) sb_fit_binary_model_selector(x, outcome, id, params, group_ids, TRUE)
      else sb_fit_binary_model_selector(x, outcome, id, params, group_ids)
    },
    apply = function(fitted, x) {
      x <- sb_validate_predictor_frame(x, "x")
      if (!inherits(fitted, "stormbreaker_selector_fit") || fitted$id != id ||
          !identical(names(x), fitted$source_features)) sb_abort(
        "Predictor columns do not match the fitted selector.", "stormbreaker_selector_schema_mismatch")
      x[fitted$selected_features]
    })
}

# ---- selector-survival ----

# ---- selector-mime-survival ----

sb_mime_rsf_selector_ids <- function() {
  c("rsf_minimal_depth", "rsf_then_stepcox")
}

sb_mime_tuned_selector_ids <- function() {
  c(
    "coxboost_nonzero",
    "lasso_nonzero",
    "enet_cox_nonzero",
    "coxboost_then_stepcox",
    "lasso_then_stepcox"
  )
}

sb_mime_composite_selector_ids <- function() {
  c(
    "rsf_then_stepcox",
    "coxboost_then_stepcox",
    "lasso_then_stepcox"
  )
}

sb_mime_survival_selector_ids <- function() {
  c(
    "rsf_minimal_depth",
    "coxboost_nonzero",
    "lasso_nonzero",
    "enet_cox_nonzero",
    sb_mime_composite_selector_ids()
  )
}

sb_validate_mime_survival_selector_inputs <- function(x, time, event,
                                                       label) {
  inputs <- sb_validate_classical_survival_inputs(x, time, event, label)
  if (any(inputs$time <= 0)) {
    sb_abort(
      paste0(label, " selection requires strictly positive survival times."),
      "stormbreaker_mime_selector_nonpositive_time"
    )
  }
  inputs
}

sb_validate_mime_selector_features <- function(selected_features,
                                                source_features,
                                                label) {
  valid <- is.character(selected_features) &&
    length(selected_features) > 0L &&
    !anyNA(selected_features) &&
    anyDuplicated(selected_features) == 0L &&
    all(selected_features %in% source_features)
  if (!valid) {
    sb_abort(
      paste0(
        label, " did not return a nonempty set of valid predictors (returned ",
        length(selected_features), " of ", length(source_features),
        " candidate predictors)."
      ),
      "stormbreaker_mime_selector_too_few_features"
    )
  }
  unname(selected_features)
}

sb_apply_mime_survival_selector <- function(fitted, x) {
  valid_fit <- inherits(fitted, "stormbreaker_selector_fit") &&
    fitted$id %in% c(
      "rsf_minimal_depth",
      "coxboost_nonzero",
      "lasso_nonzero",
      "enet_cox_nonzero",
      "rsf_then_stepcox",
      "coxboost_then_stepcox",
      "lasso_then_stepcox"
    ) &&
    is.character(fitted$source_features) &&
    is.character(fitted$selected_features) &&
    length(fitted$selected_features) > 0L &&
    anyDuplicated(fitted$selected_features) == 0L &&
    all(fitted$selected_features %in% fitted$source_features)
  if (!valid_fit) {
    sb_abort(
      "The fitted selector is not a valid MIME survival selector state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      paste0(
        "Predictor columns must exactly match the fitted selector schema: ",
        paste(fitted$source_features, collapse = ", "),
        "."
      ),
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_rsf_parse_tree_in_place <- function(recursiveObject, ...) {
  # Keep the native traversal, arithmetic and diagnostics. Reference containers
  # avoid copying the accumulated arrays at every recursive node visit.
  state <- list2env(recursiveObject, parent = emptyenv())
  state$subtree <- list2env(state$subtree, parent = emptyenv())
  native <- get("rfsrcParseTree", asNamespace("randomForestSRC"))
  result <- native(state, ...)
  result$subtree <- as.list(result$subtree)
  as.list(result)
}

sb_rsf_minimal_depth_selection <- function(engine) {
  # Rebind only this invocation; never alter the dependency's namespace.
  native <- asNamespace("randomForestSRC")
  scope <- new.env(parent = native)
  select <- get("var.select.rfsrc", native)
  maximum <- randomForestSRC::max.subtree
  environment(select) <- environment(maximum) <- scope
  scope$max.subtree <- maximum
  scope$rfsrcParseTree <- sb_rsf_parse_tree_in_place
  select(object = engine, conservative = "high", verbose = FALSE)
}

sb_fit_rsf_minimal_depth <- function(x, time, event, nodesize, seed) {
  inputs <- sb_validate_mime_survival_selector_inputs(
    x,
    time,
    event,
    "RSF minimal-depth"
  )
  if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
    sb_abort(
      "Package `randomForestSRC` is required for RSF feature selection.",
      "stormbreaker_missing_randomforestsrc"
    )
  }
  formula_data <- sb_additional_survival_formula_data(
    inputs$x,
    inputs$time,
    inputs$event,
    random_forest_src = TRUE
  )
  engine <- sb_classical_survival_stop_on_warning(
    withr::with_options(
      list(rf.cores = 1L, mc.cores = 1L),
      randomForestSRC::rfsrc(
        formula_data$formula,
        data = formula_data$data,
        ntree = 1000L,
        nodesize = nodesize,
        splitrule = "logrank",
        importance = TRUE,
        proximity = FALSE,
        forest = TRUE,
        seed = seed
      )
    ),
    "RSF minimal-depth selector fit",
    "stormbreaker_rsf_selector_warning"
  )
  selection <- sb_classical_survival_stop_on_warning(
    withr::with_options(
      list(rf.cores = 1L, mc.cores = 1L),
      sb_rsf_minimal_depth_selection(engine)
    ),
    "RSF minimal-depth variable selection",
    "stormbreaker_rsf_selector_warning"
  )
  selected_features <- sb_validate_mime_selector_features(
    selection$topvars,
    names(inputs$x),
    "RSF minimal-depth"
  )
  native_scores <- selection$varselect
  positions <- match(names(inputs$x), rownames(native_scores))
  if (anyNA(positions)) {
    sb_abort(
      "RSF minimal-depth selection returned an incomplete score table.",
      "stormbreaker_invalid_rsf_selector"
    )
  }
  scores <- data.frame(
    feature = names(inputs$x),
    minimal_depth = as.numeric(native_scores$depth[positions]),
    importance = as.numeric(native_scores$vimp[positions]),
    selected = names(inputs$x) %in% selected_features,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      id = "rsf_minimal_depth",
      source_features = names(inputs$x),
      selected_features = selected_features,
      scores = scores,
      conservative = "high",
      num_trees = 1000L,
      nodesize = nodesize,
      splitrule = "logrank",
      seed = seed,
      model_size = as.integer(length(selected_features)),
      error_rate = as.numeric(selection$err.rate)
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_rsf_minimal_depth_selector_adapter <- function(nodesize = 5L,
                                                   seed = 1L) {
  nodesize <- sb_check_positive_integer(
    nodesize,
    "selector_rsf_nodesize",
    minimum = 1L
  )
  seed <- sb_check_positive_integer(seed, "selector_seed", minimum = 1L)
  list(
    id = "rsf_minimal_depth",
    params = list(nodesize = nodesize, seed = seed),
    needs_outcome = TRUE,
    fit = function(x, time, event) {
      sb_fit_rsf_minimal_depth(x, time, event, nodesize, seed)
    },
    apply = sb_apply_mime_survival_selector
  )
}

sb_fit_coxboost_nonzero <- function(x, time, event, inner_folds, seed, capture_history = FALSE) {
  params <- sb_coxboost_cv_params(inner_folds = inner_folds, seed = seed)
  native <- sb_fit_coxboost_cv_native(
    x,
    time,
    event,
    params,
    label = "CoxBoost selector"
  )
  inputs <- native$inputs
  coefficient <- native$native_coefficient
  selected_features <- sb_validate_mime_selector_features(
    names(inputs$x)[coefficient != 0],
    names(inputs$x),
    "CoxBoost"
  )
  scores <- data.frame(
    feature = names(inputs$x),
    coefficient = unname(coefficient),
    selected = coefficient != 0,
    stringsAsFactors = FALSE
  )
  structure(
    list(
      id = "coxboost_nonzero",
      source_features = names(inputs$x),
      selected_features = selected_features,
      scores = scores,
      penalty = native$penalty,
      start_penalty = params$start_penalty,
      optimal_step = native$optimal_step,
      maximum_steps = params$maximum_steps,
      inner_folds = inner_folds,
      fold_id = native$fold_id,
      seed = seed,
      tuning_type = params$tuning_type,
      coefficient_path = if (capture_history) sweep(as.matrix(native$engine$coefficients),
        2, native$engine$sdx, "/") else NULL,
      tuning_curve = native$curve
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_coxboost_nonzero_selector_adapter <- function(inner_folds = 10L,
                                                  seed = 1L, capture_history = FALSE) {
  inner_folds <- sb_check_positive_integer(
    inner_folds,
    "selector_inner_folds",
    minimum = 3L
  )
  seed <- sb_check_positive_integer(seed, "selector_seed", minimum = 1L)
  list(
    id = "coxboost_nonzero",
    params = c(list(inner_folds = inner_folds, seed = seed),
      if (capture_history) list(capture_history = TRUE)),
    needs_outcome = TRUE,
    fit = function(x, time, event) {
      sb_fit_coxboost_nonzero(x, time, event, inner_folds, seed, capture_history)
    },
    apply = sb_apply_mime_survival_selector
  )
}

sb_fit_lasso_nonzero <- function(x, time, event, inner_folds, seed, group_ids = NULL,
                                  alpha = 1, id = "lasso_nonzero", capture_history = FALSE) {
  inputs <- sb_validate_mime_survival_selector_inputs(
    x,
    time,
    event,
    "Lasso"
  )
  tuning <- sb_select_glmnet_cox_lambda(
    inputs$x,
    inputs$time,
    inputs$event,
    inner_folds = inner_folds,
    seed = seed,
    alpha = alpha,
    standardize = TRUE,
    selection_rule = "lambda.min",
    measure = "deviance",
    group_ids = group_ids,
    retain_analysis_fit = TRUE
  )
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = alpha,
    lambda = tuning$lambda,
    standardize = TRUE
  )
  # cv.glmnet has already fitted this exact full-analysis lambda path.
  engine <- tuning$analysis_fit
  attr(engine, "stormbreaker_features") <- names(inputs$x)
  attr(engine, "stormbreaker_lambda") <- tuning$lambda
  scores <- learner$extract_coefficients(engine)
  selected_features <- sb_validate_mime_selector_features(
    scores$feature[scores$nonzero],
    names(inputs$x),
    "Lasso"
  )

  structure(
    list(
      id = id,
      source_features = names(inputs$x),
      selected_features = selected_features,
      scores = data.frame(
        feature = scores$feature,
        coefficient = scores$coefficient,
        selected = scores$nonzero,
        stringsAsFactors = FALSE
      ),
      lambda = tuning$lambda,
      alpha = alpha,
      standardize = TRUE,
      selection_rule = "lambda.min",
      measure = "deviance",
      inner_folds = inner_folds,
      fold_id = tuning$fold_id,
      seed = seed,
      tuning_curve = tuning$curve,
      coefficient_path = if (capture_history) list(lambda = engine$lambda,
        beta = engine$beta) else NULL,
      tuning_warnings = tuning$warnings
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_lasso_nonzero_selector_adapter <- function(inner_folds = 10L,
                                               seed = 1L, id = "lasso_nonzero", capture_history = FALSE) {
  inner_folds <- sb_check_positive_integer(
    inner_folds,
    "selector_inner_folds",
    minimum = 3L
  )
  seed <- sb_check_positive_integer(seed, "selector_seed", minimum = 1L)
  list(
    id = id,
    params = c(list(inner_folds = inner_folds, seed = seed),
      if (capture_history) list(capture_history = TRUE)),
    needs_outcome = TRUE,
    fit = function(x, time, event, group_ids = NULL) {
      sb_fit_lasso_nonzero(x, time, event, inner_folds, seed, group_ids,
        alpha = if (id == "enet_cox_nonzero") .5 else 1, id = id,
        capture_history = capture_history)
    },
    apply = sb_apply_mime_survival_selector
  )
}

sb_sequential_survival_selector_adapter <- function(first, second, id) {
  valid_adapter <- function(adapter) {
    is.list(adapter) && is.character(adapter$id) &&
      length(adapter$id) == 1L && is.list(adapter$params) &&
      isTRUE(adapter$needs_outcome) && is.function(adapter$fit) &&
      is.function(adapter$apply)
  }
  valid_id <- is.character(id) && length(id) == 1L && !is.na(id) &&
    id %in% c(
      "rsf_then_stepcox",
      "coxboost_then_stepcox",
      "lasso_then_stepcox"
    )
  if (!valid_adapter(first) || !valid_adapter(second) || !valid_id) {
    sb_abort(
      "Sequential survival selection requires two valid selector adapters.",
      "stormbreaker_bad_sequential_selector"
    )
  }

  list(
    id = id,
    params = list(
      first = list(id = first$id, params = first$params),
      second = list(id = second$id, params = second$params)
    ),
    needs_outcome = TRUE,
    fit = function(x, time, event) {
      inputs <- sb_validate_mime_survival_selector_inputs(
        x,
        time,
        event,
        "Sequential survival"
      )
      endpoint <- list(task_type = "survival", time = inputs$time,
        event = inputs$event)
      fitted_first <- sb_fit_runtime_selector(
        first,
        inputs$x,
        endpoint
      )
      selected_first <- first$apply(fitted_first, inputs$x)
      fitted_second <- sb_fit_runtime_selector(
        second,
        selected_first,
        endpoint
      )
      selected_features <- fitted_second$selected_features
      valid_selection <- is.character(selected_features) &&
        length(selected_features) > 0L && !anyNA(selected_features) &&
        anyDuplicated(selected_features) == 0L &&
        all(selected_features %in% fitted_first$selected_features)
      if (!valid_selection) {
        sb_abort(
          "Sequential survival selection returned invalid final features.",
          "stormbreaker_invalid_sequential_selector"
        )
      }
      scores <- data.frame(
        feature = names(inputs$x),
        first_stage_selected = names(inputs$x) %in%
          fitted_first$selected_features,
        selected = names(inputs$x) %in% selected_features,
        stringsAsFactors = FALSE
      )
      structure(
        list(
          id = id,
          source_features = names(inputs$x),
          selected_features = selected_features,
          scores = scores,
          stages = list(first = fitted_first, second = fitted_second)
        ),
        class = "stormbreaker_selector_fit"
      )
    },
    apply = sb_apply_mime_survival_selector
  )
}

sb_rsf_then_stepcox_selector_adapter <- function(nodesize = 5L,
                                                  seed = 1L,
                                                  direction) {
  adapter <- sb_sequential_survival_selector_adapter(
    sb_rsf_minimal_depth_selector_adapter(nodesize, seed),
    sb_stepcox_selector_adapter(direction),
    "rsf_then_stepcox"
  )
  adapter$params <- list(
    nodesize = sb_check_positive_integer(
      nodesize,
      "selector_rsf_nodesize",
      minimum = 1L
    ),
    seed = sb_check_positive_integer(
      seed,
      "selector_seed",
      minimum = 1L
    ),
    direction = sb_validate_stepcox_direction(direction)
  )
  adapter
}

sb_coxboost_then_stepcox_selector_adapter <- function(inner_folds = 10L,
                                                       seed = 1L,
                                                       direction) {
  adapter <- sb_sequential_survival_selector_adapter(
    sb_coxboost_nonzero_selector_adapter(inner_folds, seed),
    sb_stepcox_selector_adapter(direction),
    "coxboost_then_stepcox"
  )
  adapter$params <- list(
    inner_folds = sb_check_positive_integer(
      inner_folds,
      "selector_inner_folds",
      minimum = 3L
    ),
    seed = sb_check_positive_integer(
      seed,
      "selector_seed",
      minimum = 1L
    ),
    direction = sb_validate_stepcox_direction(direction)
  )
  adapter
}

sb_lasso_then_stepcox_selector_adapter <- function(inner_folds = 10L,
                                                    seed = 1L,
                                                    direction) {
  adapter <- sb_sequential_survival_selector_adapter(
    sb_lasso_nonzero_selector_adapter(inner_folds, seed),
    sb_stepcox_selector_adapter(direction),
    "lasso_then_stepcox"
  )
  adapter$params <- list(
    inner_folds = sb_check_positive_integer(
      inner_folds,
      "selector_inner_folds",
      minimum = 3L
    ),
    seed = sb_check_positive_integer(
      seed,
      "selector_seed",
      minimum = 1L
    ),
    direction = sb_validate_stepcox_direction(direction)
  )
  adapter
}

# ---- selector-stepcox ----

sb_validate_stepcox_direction <- function(direction) {
  valid <- is.character(direction) && length(direction) == 1L &&
    !is.na(direction) && direction %in% c("forward", "backward", "both")
  if (!valid) {
    sb_abort(
      paste0(
        "`stepcox_direction` must be one of `forward`, `backward`, or `both`."
      ),
      "stormbreaker_bad_stepcox_direction"
    )
  }
  unname(as.character(direction))
}

sb_stepcox_selection_path <- function(model, internal_features,
                                      source_features, final_aic) {
  if (is.null(model$anova)) {
    return(data.frame(
      Step = "",
      AIC = final_aic,
      stringsAsFactors = FALSE
    ))
  }
  path <- as.data.frame(model$anova, stringsAsFactors = FALSE)
  actions <- as.character(path$Step)
  actions[is.na(actions)] <- ""
  for (index in seq_along(actions)) {
    action <- trimws(actions[[index]])
    if (!nzchar(action)) {
      actions[[index]] <- ""
      next
    }
    operation <- substr(action, 1L, 1L)
    internal <- trimws(substring(action, 2L))
    position <- match(internal, internal_features)
    if (!is.na(position) && operation %in% c("+", "-")) {
      actions[[index]] <- paste(operation, source_features[[position]])
    }
  }
  path$Step <- actions
  rownames(path) <- NULL
  path
}

sb_fit_stepcox <- function(x, time, event, direction) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  direction <- sb_validate_stepcox_direction(direction)
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

  source_features <- names(x)
  internal_features <- paste0(
    ".stormbreaker_stepcox_x",
    seq_along(source_features)
  )
  analysis <- x
  names(analysis) <- internal_features
  analysis$.stormbreaker_stepcox_time <- time
  analysis$.stormbreaker_stepcox_event <- event
  response <- paste0(
    "survival::Surv(",
    ".stormbreaker_stepcox_time, .stormbreaker_stepcox_event)"
  )
  full_formula <- stats::reformulate(
    internal_features,
    response = response,
    env = environment()
  )
  null_formula <- stats::as.formula(
    paste(response, "~ 1"),
    env = environment()
  )
  start_formula <- if (identical(direction, "forward")) {
    null_formula
  } else {
    full_formula
  }
  # Tighten likelihood convergence without changing the native warning or
  # Cholesky thresholds. Candidate-model warnings still stop selection.
  cox_control <- survival::coxph.control(
    eps = 1e-11,
    iter.max = 50L,
    toler.inf = survival::coxph.control()$toler.inf
  )

  model <- withCallingHandlers(
    {
      start <- survival::coxph(
        start_formula,
        data = analysis,
        ties = "efron",
        singular.ok = FALSE,
        x = TRUE,
        y = TRUE,
        model = TRUE,
        control = cox_control
      )
      stats::step(
        start,
        scope = list(lower = null_formula, upper = full_formula),
        direction = direction,
        trace = 0,
        k = 2
      )
    },
    warning = function(warning) {
      sb_abort(
        paste0(
          "The StepCox selector raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_stepcox_warning",
        parent = warning
      )
    }
  )
  if (NROW(model$y) != nrow(analysis)) {
    sb_abort(
      "The StepCox selector dropped analysis patients.",
      "stormbreaker_selector_rows_dropped"
    )
  }

  coefficients <- stats::coef(model)
  if (length(coefficients) < 1L) {
    sb_abort(
      paste0(
        "StepCox selected no predictor with direction `", direction, "`."
      ),
      "stormbreaker_stepcox_empty_selection"
    )
  }
  selected_positions <- match(names(coefficients), internal_features)
  valid_coefficients <- !anyNA(selected_positions) &&
    !anyNA(coefficients) && all(is.finite(coefficients))
  if (!valid_coefficients) {
    sb_abort(
      "The StepCox selector produced an invalid fitted model.",
      "stormbreaker_invalid_stepcox_fit"
    )
  }
  selected_features <- source_features[selected_positions]
  final_aic <- unname(stats::extractAIC(model, k = 2)[[2L]])
  if (!is.finite(final_aic)) {
    sb_abort(
      "The StepCox selector produced a non-finite AIC.",
      "stormbreaker_invalid_stepcox_fit"
    )
  }

  scores <- data.frame(
    feature = source_features,
    original_order = as.integer(seq_along(source_features)),
    selected = source_features %in% selected_features,
    coefficient = NA_real_,
    stringsAsFactors = FALSE
  )
  scores$coefficient[selected_positions] <- unname(coefficients)

  structure(
    list(
      id = "stepcox",
      direction = direction,
      criterion = "AIC",
      penalty = 2,
      start_model = if (identical(direction, "forward")) "null" else "full",
      ties = "efron",
      cox_control = cox_control,
      source_features = source_features,
      selected_features = selected_features,
      scores = scores,
      final_aic = final_aic,
      selection_path = sb_stepcox_selection_path(
        model,
        internal_features,
        source_features,
        final_aic
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_apply_stepcox <- function(fitted, x) {
  valid_fit <- inherits(fitted, "stormbreaker_selector_fit") &&
    identical(fitted$id, "stepcox") &&
    is.character(fitted$direction) && length(fitted$direction) == 1L &&
    fitted$direction %in% c("forward", "backward", "both") &&
    is.character(fitted$source_features) &&
    is.character(fitted$selected_features) &&
    length(fitted$selected_features) > 0L &&
    anyDuplicated(fitted$selected_features) == 0L &&
    all(fitted$selected_features %in% fitted$source_features)
  if (!valid_fit) {
    sb_abort(
      "The fitted selector is not a valid StepCox state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      paste0(
        "Predictor columns must exactly match the fitted selector schema: ",
        paste(fitted$source_features, collapse = ", "),
        "."
      ),
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_stepcox_selector_adapter <- function(direction) {
  direction <- sb_validate_stepcox_direction(direction)
  list(
    id = "stepcox",
    params = list(direction = direction),
    needs_outcome = TRUE,
    fit = function(x, time, event) {
      sb_fit_stepcox(x, time, event, direction)
    },
    apply = sb_apply_stepcox
  )
}

# ---- selector-survival-boosting ----

sb_survival_boosting_selector_ids <- function() c("gbm_cox_importance", "glmboost_cox_nonzero", "superpc_survival_scores", "gamboost_cox_nonzero", "xgboost_cox_importance")

sb_mboost_spline_selector_scores <- function(model, features) {
  beta <- stats::coef(model)
  index <- match(sub("^bbs\\(([^,)]+).*", "\\1", names(beta)), features)
  keep <- !is.na(index)
  score <- numeric(length(features))
  score[index[keep]] <- vapply(beta[keep], function(b) sqrt(sum(b^2)), numeric(1))
  score
}

sb_survival_boosting_selector_params <- function(id) {
  if (id == "xgboost_cox_importance") return(sb_xgboost_params("xgboost_cox"))
  if (id == "superpc_survival_scores") return(sb_superpc_params())
  if (id == "gamboost_cox_nonzero") return(sb_mboost_cox_params("gamboost_cox"))
  if (id == "gbm_cox_importance") sb_gbm_cox_params() else
    sb_mboost_cox_params("glmboost_cox")
}

sb_survival_boosting_selector_adapter <- function(id, params) {
  list(id = id, params = params, needs_outcome = TRUE,
    fit = function(x, time, event) {
      x <- sb_validate_predictor_frame(x, "x")
      work <- x
      names(work) <- paste0("feature_", seq_len(ncol(x)))
      if (id == "superpc_survival_scores") {
        inputs <- sb_validate_classical_survival_inputs(work, time, event, "SuperPC selector")
        if (!requireNamespace("superpc", quietly = TRUE)) stop("Install superpc for SuperPC selection.")
        if (ncol(work) < 2L) stop("SuperPC selection needs at least two input features.")
        model <- superpc::superpc.train(list(x = t(as.matrix(inputs$x)),
          y = inputs$time, censoring.status = inputs$event, featurenames = names(work)),
          type = "survival", s0.perc = params$s0_perc)
        score <- as.numeric(model$feature.scores)
        selected <- abs(score) >= params$threshold
        history <- numeric()
      } else if (id == "xgboost_cox_importance") {
        model <- sb_fit_xgboost(work, list(time = time, event = event), "xgboost_cox", params)
        score <- sb_xgboost_selector_gain(model)
        selected <- score > 0
        history <- numeric()
      } else if (id == "gbm_cox_importance") {
        model <- sb_fit_gbm_cox(work, time, event, params)$booster
        importance <- summary(model, n.trees = params$n_trees, plotit = FALSE)
        score <- importance$rel.inf[match(names(work), importance$var)]
        selected <- score > 0
        history <- model$train.error
      } else {
        spline <- id == "gamboost_cox_nonzero"
        model <- sb_fit_mboost_cox(work, time, event,
          if (spline) "gamboost_cox" else "glmboost_cox", params)$engine
        if (spline) score <- sb_mboost_spline_selector_scores(model, names(work)) else {
          beta <- stats::coef(model, off2int = FALSE)
          score <- numeric(ncol(work))
          index <- match(names(beta), names(work))
          keep <- !is.na(index)
          score[index[keep]] <- as.numeric(beta[keep])
        }
        selected <- score != 0
        history <- mboost::risk(model)
      }
      if (length(score) != ncol(x) || any(!is.finite(score))) stop("Invalid native selector scores.")
      if (!any(selected)) stop("No features met the selection rule; no fallback genes were added.")
      structure(list(id = id, params = params, source_features = names(x),
        selected_features = names(x)[selected],
        scores = data.frame(feature = names(x), score = score, selected = selected),
        evidence = list(training_loss = data.frame(
          step = seq_along(history) - as.integer(id %in% c("glmboost_cox_nonzero", "gamboost_cox_nonzero")), loss = as.numeric(history)))),
        class = "stormbreaker_selector_fit")
    },
    apply = function(fitted, x) {
      x <- sb_validate_predictor_frame(x, "x")
      if (!identical(names(x), fitted$source_features)) stop("Predictor columns differ from the fitted selector.")
      x[fitted$selected_features]
    })
}

# ---- selector-univ-cox-top-k ----

sb_score_univariable_cox_feature <- function(feature, x, survival_response,
                                              weights, control,
                                              original_order) {
  empty <- data.frame(
    feature = feature,
    original_order = as.integer(original_order),
    coefficient = NA_real_,
    standard_error = NA_real_,
    wald_z = NA_real_,
    p_value = NA_real_,
    status = "failed",
    error_class = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
  tryCatch(
    {
      values <- x[[feature]]
      n_above_minimum <- sum(values > min(values))
      n_below_maximum <- sum(values < max(values))
      if (n_above_minimum < 2L || n_below_maximum < 2L) {
        sb_abort(
          paste0(
            "Univariable Cox scoring cannot use `", feature,
            "`: fewer than two patients lie away from one expression ",
            "extreme, so a single patient could define the feature score."
          ),
          "stormbreaker_selector_insufficient_feature_support"
        )
      }
      model <- withCallingHandlers(
        survival::coxph.fit(
          x = matrix(values, ncol = 1L),
          y = survival_response,
          strata = NULL,
          offset = NULL,
          init = NULL,
          control = control,
          weights = weights,
          method = "efron",
          rownames = NULL,
          resid = FALSE,
          nocenter = c(-1, 0, 1)
        ),
        warning = function(warning) {
          sb_abort(
            paste0(
              "Univariable Cox scoring raised a warning for `",
              feature, "`: ", conditionMessage(warning)
            ),
            "stormbreaker_selector_cox_warning",
            parent = warning
          )
        }
      )
      coefficient <- unname(model$coefficients[[1L]])
      standard_error <- unname(sqrt(model$var[[1L]]))
      if (!is.finite(coefficient) || !is.finite(standard_error) ||
          standard_error <= 0) {
        sb_abort(
          paste0(
            "Univariable Cox scoring produced invalid statistics for `",
            feature, "`."
          ),
          "stormbreaker_invalid_selector_score"
        )
      }
      wald_z <- coefficient / standard_error
      p_value <- 2 * stats::pnorm(abs(wald_z), lower.tail = FALSE)
      if (!is.finite(wald_z) || !is.finite(p_value)) {
        sb_abort(
          paste0(
            "Univariable Cox scoring produced a non-finite Wald score for `",
            feature, "`."
          ),
          "stormbreaker_invalid_selector_score"
        )
      }
      empty$coefficient <- coefficient
      empty$standard_error <- standard_error
      empty$wald_z <- wald_z
      empty$p_value <- p_value
      empty$status <- "success"
      empty
    },
    error = function(error) {
      empty$error_class <- class(error)[1]
      empty$error_message <- conditionMessage(error)
      empty
    }
  )
}

sb_fit_univ_cox_top_k <- function(x, time, event, n_features) {
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
  if (n_features > ncol(x)) {
    sb_abort(
      paste0(
        "The selector requested ", n_features,
        " features, but the task has only ", ncol(x),
        " candidate features."
      ),
      "stormbreaker_insufficient_candidate_features"
    )
  }

  survival_response <- survival::Surv(time, event)
  weights <- rep.int(1, nrow(x))
  control <- survival::coxph.control()

  scores <- do.call(
    rbind,
    lapply(seq_along(x), function(index) {
      sb_score_univariable_cox_feature(
        names(x)[index],
        x,
        survival_response,
        weights,
        control,
        index
      )
    })
  )
  rownames(scores) <- NULL
  scores$rank <- NA_integer_
  scores$selected <- FALSE
  valid <- which(scores$status == "success" & is.finite(scores$p_value))
  ordered_valid <- valid[order(
    scores$p_value[valid],
    scores$original_order[valid]
  )]
  if (length(ordered_valid) < n_features) {
    sb_abort(
      paste0(
        "The selector requested ", n_features,
        " features, but only ", length(ordered_valid),
        " features had valid univariable Cox scores."
      ),
      "stormbreaker_insufficient_valid_features"
    )
  }
  scores$rank[ordered_valid] <- seq_along(ordered_valid)
  selected_rows <- ordered_valid[seq_len(n_features)]
  scores$selected[selected_rows] <- TRUE
  selected_features <- scores$feature[selected_rows]

  structure(
    list(
      id = "univ_cox_top_k",
      n_features = n_features,
      source_features = names(x),
      selected_features = selected_features,
      scores = scores
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_apply_univ_cox_top_k <- function(fitted, x) {
  if (!inherits(fitted, "stormbreaker_selector_fit") ||
      !identical(fitted$id, "univ_cox_top_k") ||
      !is.character(fitted$source_features) ||
      !is.character(fitted$selected_features) ||
      length(fitted$selected_features) != fitted$n_features) {
    sb_abort(
      "The fitted selector is not a valid univariable Cox Top-K state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      paste0(
        "Predictor columns must exactly match the fitted selector schema: ",
        paste(fitted$source_features, collapse = ", "),
        "."
      ),
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_univ_cox_top_k_selector_adapter <- function(n_features) {
  n_features <- sb_check_positive_integer(
    n_features,
    "n_features",
    minimum = 1L
  )
  list(
    id = "univ_cox_top_k",
    params = list(n_features = n_features),
    needs_outcome = TRUE,
    fit = function(x, time, event) {
      sb_fit_univ_cox_top_k(x, time, event, n_features)
    },
    apply = sb_apply_univ_cox_top_k
  )
}
