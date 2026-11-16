# ---- classical-binary-boosting-fixed-adapters ----

sb_classical_binary_boosting_learner_ids <- function() {
  c("adaboost_binary", "logitboost_binary")
}

sb_fastadaboost_available <- function() {
  requireNamespace("fastAdaboost", quietly = TRUE)
}

sb_catools_available <- function() {
  requireNamespace("caTools", quietly = TRUE)
}

sb_classical_binary_boosting_dependency <- function(learner) {
  switch(learner,
    adaboost_binary = list(
      field = "fastadaboost_version",
      package = "fastAdaboost",
      available = sb_fastadaboost_available()
    ),
    logitboost_binary = list(
      field = "catools_version",
      package = "caTools",
      available = sb_catools_available()
    ),
    NULL
  )
}

sb_classical_binary_boosting_params <- function(learner, iterations = NULL) {
  if (!learner %in% sb_classical_binary_boosting_learner_ids()) {
    sb_abort(
      "The requested classical binary boosting learner is not supported.",
      "stormbreaker_bad_classical_binary_boosting_learner"
    )
  }
  if (is.null(iterations)) {
    iterations <- if (identical(learner, "adaboost_binary")) 50L else 51L
  }
  iterations <- sb_check_positive_integer(
    iterations,
    if (identical(learner, "adaboost_binary")) {
      "adaboost_iterations"
    } else {
      "logitboost_iterations"
    },
    minimum = 1L
  )
  if (identical(learner, "adaboost_binary")) {
    return(list(
      iterations = iterations,
      method = "AdaBoost.M1",
      base_learner = "rpart_classification_tree",
      rpart_cp = 0,
      iteration_semantics = "maximum_until_native_convergence",
      tuning = "none"
    ))
  }
  list(
    iterations = iterations,
    method = "LogitBoost",
    base_learner = "decision_stump",
    iteration_semantics = "exact",
    tuning = "none"
  )
}



sb_fastadaboost_fit <- function(formula, data, iterations) {
  fastAdaboost::adaboost(formula, data, nIter = iterations)
}

sb_catools_logitboost_fit <- function(x, outcome, iterations) {
  caTools::LogitBoost(x, outcome, nIter = iterations)
}

sb_validate_classical_binary_boosting_engine <- function(engine, learner,
                                                         iterations) {
  valid <- if (identical(learner, "adaboost_binary")) {
    inherits(engine, "adaboost") &&
      length(engine$trees) >= 1L && length(engine$trees) <= iterations &&
      length(engine$weights) == length(engine$trees) &&
      all(vapply(engine$trees, inherits, logical(1), "rpart")) &&
      all(is.finite(engine$weights)) &&
      identical(unname(engine$classnames), c("0", "1"))
  } else {
    inherits(engine, "LogitBoost") && is.matrix(engine$Stump) &&
      nrow(engine$Stump) == iterations && ncol(engine$Stump) == 3L &&
      !anyNA(engine$Stump) && all(is.finite(engine$Stump)) &&
      identical(as.character(engine$lablist), c("0", "1"))
  }
  if (!valid) {
    sb_abort(
      "The classical binary boosting fit is incomplete or invalid.",
      "stormbreaker_invalid_classical_binary_boosting_fit"
    )
  }
  engine
}

sb_fit_classical_binary_boosting <- function(x, outcome, learner,
                                                   params) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "A classical binary boosting analysis fold must contain both classes.",
      "stormbreaker_classical_binary_boosting_single_class"
    )
  }
  params <- sb_classical_binary_boosting_params(learner, params$iterations)

  if (identical(learner, "adaboost_binary")) {
    if (!sb_fastadaboost_available()) {
      sb_abort(
        "Package `fastAdaboost` is required for fixed AdaBoost.M1.",
        "stormbreaker_missing_fastadaboost"
      )
    }
    response <- ".stormbreaker_outcome"
    while (response %in% names(x)) {
      response <- paste0(".", response)
    }
    training <- data.frame(
      factor(outcome, levels = c(0L, 1L)),
      x,
      check.names = FALSE
    )
    names(training)[[1L]] <- response
    formula <- stats::as.formula(
      paste0("`", response, "` ~ ."),
      env = baseenv()
    )
    engine <- sb_fastadaboost_fit(formula, training, params$iterations)
  } else {
    if (!sb_catools_available()) {
      sb_abort(
        "Package `caTools` is required for fixed LogitBoost.",
        "stormbreaker_missing_catools"
      )
    }
    engine <- sb_catools_logitboost_fit(
      as.matrix(x),
      outcome,
      params$iterations
    )
  }
  engine <- sb_validate_classical_binary_boosting_engine(
    engine,
    learner,
    params$iterations
  )

  structure(
    list(
      learner = learner,
      features = names(x),
      params = params,
      iterations = params$iterations,
      fitted_iterations = if (identical(
        learner,
        "adaboost_binary"
      )) {
        as.integer(length(engine$trees))
      } else {
        as.integer(nrow(engine$Stump))
      },
      engine = engine
    ),
    class = "stormbreaker_classical_binary_boosting_fixed"
  )
}

sb_predict_classical_binary_boosting <- function(model, x) {
  if (!inherits(model, "stormbreaker_classical_binary_boosting_fixed") ||
    !model$learner %in% sb_classical_binary_boosting_learner_ids()) {
    sb_abort(
      "The supplied model is not a fitted classical binary boosting learner.",
      "stormbreaker_bad_classical_binary_boosting_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the binary boosting schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }

  if (identical(model$learner, "adaboost_binary")) {
    if (!sb_fastadaboost_available()) {
      sb_abort(
        "Package `fastAdaboost` is required for fixed AdaBoost.M1 prediction.",
        "stormbreaker_missing_fastadaboost"
      )
    }
    prediction <- stats::predict(model$engine, newdata = x)
    positive_column <- match("1", unname(model$engine$classnames))
    probability <- prediction$prob[, positive_column]
  } else {
    if (!sb_catools_available()) {
      sb_abort(
        "Package `caTools` is required for fixed LogitBoost prediction.",
        "stormbreaker_missing_catools"
      )
    }
    prediction <- stats::predict(
      model$engine,
      as.matrix(x),
      type = "raw",
      nIter = model$iterations
    )
    positive_column <- match("1", as.character(model$engine$lablist))
    probability <- prediction[, positive_column]
  }
  probability <- unname(as.numeric(probability))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      "Classical binary boosting produced invalid class-one probabilities.",
      "stormbreaker_invalid_classical_binary_boosting_prediction"
    )
  }
  probability
}

sb_classical_binary_boosting_feature_usage <- function(model) {
  split_features <- if (identical(model$learner, "adaboost_binary")) {
    unlist(lapply(model$engine$trees, function(tree) {
      variables <- as.character(tree$frame$var)
      variables[variables != "<leaf>"]
    }), use.names = FALSE)
  } else {
    indices <- as.integer(model$engine$Stump[, "feature"])
    if (any(indices < 1L | indices > length(model$features))) {
      sb_abort(
        "LogitBoost returned a stump outside the fitted feature schema.",
        "stormbreaker_invalid_classical_binary_boosting_fit"
      )
    }
    model$features[indices]
  }
  if (length(setdiff(split_features, model$features)) > 0L) {
    sb_abort(
      "Binary boosting returned a split outside the fitted feature schema.",
      "stormbreaker_invalid_classical_binary_boosting_fit"
    )
  }
  as.numeric(tabulate(
    match(split_features, model$features),
    nbins = length(model$features)
  ))
}

sb_extract_classical_binary_boosting_model <- function(model) {
  if (!inherits(model, "stormbreaker_classical_binary_boosting_fixed")) {
    sb_abort(
      "The supplied model is not a fitted classical binary boosting learner.",
      "stormbreaker_bad_classical_binary_boosting_fit"
    )
  }
  list(
    iterations = model$iterations,
    fitted_iterations = model$fitted_iterations,
    features = model$features,
    params = model$params,
    importance_type = "split_count",
    model_serialization = paste0(
      "R_xdr_v3_",
      if (identical(model$learner, "adaboost_binary")) {
        "adaboost"
      } else {
        "LogitBoost"
      }
    ),
    model_raw = serialize(
      model$engine,
      connection = NULL,
      ascii = FALSE,
      xdr = TRUE,
      version = 3L
    ),
    feature_importance = data.frame(
      feature = model$features,
      importance = sb_classical_binary_boosting_feature_usage(model),
      stringsAsFactors = FALSE
    )
  )
}

sb_classical_binary_boosting_learner_adapter <- function(learner,
                                                         iterations = NULL) {
  params <- sb_classical_binary_boosting_params(learner, iterations)
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_classical_binary_boosting(x, outcome, learner, params)
    },
    predict = sb_predict_classical_binary_boosting,
    extract_tree_model = sb_extract_classical_binary_boosting_model,
    extract_state = function(model, context = list()) {
      sb_extract_classical_binary_boosting_learner_state(model)
    },
    predict_state = sb_predict_classical_binary_boosting_learner_state
  )
}

# ---- diagnostics-boosting ----

# ---- algorithm-classical-boosting ----

sb_classical_boosting_diagnostic <- function(fit, type) {
  id <- fit$pipeline$learner
  state <- fit$learner_state
  ada <- id == "adaboost_binary"
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("split_usage", if (ada) c("tree_weights", "tree_size") else "stump_sequence")) {
    sb_abort("Use split_usage for either learner; tree_weights/tree_size for AdaBoost or stump_sequence for LogitBoost. Original evaluation histories are not retained.",
      "stormbreaker_unavailable_diagnostic")
  }
  engine <- sb_validate_classical_binary_boosting_engine(unserialize(state$model_raw), id, state$iterations)
  n <- if (ada) length(engine$trees) else nrow(engine$Stump)
  native <- NULL
  if (type == "split_usage") {
    tab <- state$feature_importance
    if (!is.data.frame(tab) || !identical(tab$feature, state$features) ||
        !is.numeric(tab$importance) || any(!is.finite(tab$importance)) ||
        any(tab$importance < 0) || any(tab$importance != floor(tab$importance))) sb_abort(
      "Original boosting split counts are invalid.", "stormbreaker_invalid_diagnostic")
    native <- tab
    tab <- data.frame(feature = tab$feature, count = tab$importance)
    measure <- "Original split usage"
    interpretation <- "Unweighted primary split count; not gain, permutation importance, effect direction or selection stability."
  } else if (type == "tree_weights") {
    native <- engine$weights
    tab <- data.frame(iteration = seq_len(n), weight = as.numeric(native))
    measure <- "Original AdaBoost tree weights"
    interpretation <- "Saved ensemble voting weights, not tree accuracy, probabilities or a training-loss curve."
  } else if (type == "tree_size") {
    tab <- data.frame(tree = seq_len(n), terminal_nodes = vapply(engine$trees,
      function(tree) sum(tree$frame$var == "<leaf>"), integer(1)))
    measure <- "Original AdaBoost tree-size distribution"
    interpretation <- "Every original fitted tree counted; one tree is not the complete ensemble."
  } else {
    native <- engine$Stump
    index <- as.integer(native[, 1])
    if (any(index < 1L) || any(index > length(state$features)) ||
        any(index != native[, 1]) || any(!native[, 3] %in% c(-1, 1))) sb_abort(
      "Native LogitBoost stump indices or directions are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(iteration = seq_len(n), feature = state$features[index],
      native_feature_index = index, threshold = as.numeric(native[, 2]),
      contribution_at_or_below = as.numeric(native[, 3]),
      contribution_above = -as.numeric(native[, 3]))
    measure <- "Original LogitBoost stump sequence"
    interpretation <- "Colour: contribution to class-1 log-odds when value <= the saved threshold; above it, the sign reverses.\nExact feature-specific thresholds remain in the table; these are model splits, not clinical cutoffs."
  }
  structure(list(type = paste0("classical_boosting_", type), table = tab, native_table = native,
    learner = id, pipeline_id = fit$pipeline$id, fitted_iterations = n,
    requested_iterations = state$iterations, features = state$features,
    measure = measure, interpretation = interpretation,
    source = "original saved boosting model; no refit or patient prediction"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_classical_boosting_diagnostic <- function(x) {
  iteration <- weight <- feature <- count <- direction <- NULL
  tab <- x$table
  if (x$type == "classical_boosting_tree_size") {
    z <- x
    z$type <- "forest_tree_size"
    p <- sb_plot_forest_diagnostic(z) +
      ggplot2::labs(x = "Terminal nodes per original AdaBoost tree")
  } else if (x$type == "classical_boosting_tree_weights") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(iteration, weight)) +
      ggplot2::geom_segment(ggplot2::aes(xend = iteration, yend = 0), colour = "#2878A0") +
      ggplot2::geom_point(colour = "#2878A0") +
      ggplot2::labs(x = "Original tree iteration", y = "Native voting weight") + sb_viz_theme()
  } else if (x$type == "classical_boosting_split_usage") {
    tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$count)])
    p <- ggplot2::ggplot(tab, ggplot2::aes(count, feature)) +
      ggplot2::geom_col(fill = "#2878A0", width = .65) +
      ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::labs(x = "Primary split count", y = NULL) + sb_viz_theme()
  } else {
    tab$feature <- factor(tab$feature, levels = rev(x$features))
    tab$direction <- factor(tab$contribution_at_or_below, levels = c(-1, 1), labels = c("-1", "+1"))
    p <- ggplot2::ggplot(tab, ggplot2::aes(iteration, feature, colour = direction)) +
      ggplot2::geom_point(size = 2) + ggplot2::scale_y_discrete(drop = FALSE) +
      ggplot2::scale_colour_manual(values = c("-1" = "#2878A0", "+1" = "#D77935"), drop = FALSE) +
      ggplot2::labs(x = "Original stump iteration", y = NULL, colour = "At/below threshold") + sb_viz_theme()
  }
  p + ggplot2::labs(title = x$measure,
    subtitle = paste0(x$learner, "; fitted iterations: ", x$fitted_iterations,
      if (x$learner == "adaboost_binary") "; configured maximum: " else "; configured iterations: ",
      x$requested_iterations), caption = x$interpretation)
}

# ---- algorithm-hinge-boosting ----

sb_hinge_boosting_diagnostic <- function(fit, type) {
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("training_loss", "component_usage", "coefficient_path")) sb_abort(
    "Hinge boosting supports training_loss, component_usage and coefficient_path.",
    "stormbreaker_unavailable_diagnostic")
  model <- fit$diagnostic_artifacts$learner
  path <- model$diagnostic_path
  if (is.null(path)) sb_abort(
    "Original hinge-boosting iterations were not captured; enable capture_diagnostics on a new final fit. Historical compact fits contain only the final risk.",
    "stormbreaker_missing_diagnostic_artifacts")
  n <- model$params$iterations
  features <- model$features
  if (!identical(path$family, "hinge") || !identical(path$learner, "ls") ||
      length(path$risk) != n || length(path$selected) != n || length(path$updates) != n ||
      any(!is.finite(path$risk)) || any(!is.finite(path$updates)) ||
      any(!is.finite(path$selected)) || any(path$selected != floor(path$selected)) ||
      any(path$selected < 1L) || any(path$selected > length(features)) ||
      length(path$meanx) != length(features) || any(!is.finite(path$meanx)) ||
      length(path$offset) != 1L || !is.finite(path$offset) ||
      length(path$learning_rate) != 1L || !is.finite(path$learning_rate) ||
      !is.null(path$intercept_updates)) sb_abort(
    "Captured native hinge-boosting path is invalid or uses unsupported intercept updates.",
    "stormbreaker_invalid_diagnostic")
  selections <- data.frame(iteration = seq_len(n), feature = features[path$selected],
    native_feature_index = path$selected, native_update = path$updates)
  if (type == "training_loss") {
    tab <- data.frame(iteration = seq_len(n), loss = as.numeric(path$risk))
    measure <- "Native mean hinge loss after each update"
  } else if (type == "component_usage") {
    tab <- data.frame(feature = features, updates = tabulate(path$selected, nbins = length(features)))
    measure <- "Native hinge-boosting component updates"
  } else {
    beta <- matrix(0, nrow = n + 1L, ncol = length(features))
    for (i in seq_len(n)) {
      beta[i + 1L, ] <- beta[i, ]
      j <- path$selected[i]
      beta[i + 1L, j] <- beta[i + 1L, j] + path$learning_rate * path$updates[i]
    }
    intercept <- as.numeric(path$offset - beta %*% path$meanx)
    if (!isTRUE(all.equal(as.numeric(beta[n + 1L, ]), as.numeric(model$engine$coefficient), tolerance = 1e-12)) ||
        !isTRUE(all.equal(intercept[n + 1L], model$engine$intercept, tolerance = 1e-12))) sb_abort(
      "Hinge path endpoint disagrees with the original margin equation.", "stormbreaker_invalid_diagnostic")
    beta <- cbind(intercept, beta)
    tab <- data.frame(iteration = rep(seq.int(0L, n), times = ncol(beta)),
      feature = rep(c("(Intercept)", features), each = nrow(beta)), coefficient = as.vector(beta))
    measure <- "Original hinge-margin coefficient path"
  }
  structure(list(type = paste0("hinge_", type), table = tab, selections = selections,
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id, measure = measure,
    source = "captured original hinge fit; no additional training",
    interpretation = "Native hinge-margin scale, not probability or log-odds. Original stopping retained.\nTraining diagnostics, not external validation, uncertainty intervals or causal effects."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_hinge_boosting_diagnostic <- function(x) {
  iteration <- loss <- feature <- coefficient <- NULL
  if (x$type == "hinge_component_usage") {
    z <- x
    z$type <- "boosting_component_usage"
    return(sb_plot_boosting_diagnostic(z) + ggplot2::labs(title = x$measure))
  }
  p <- if (x$type == "hinge_training_loss") {
    ggplot2::ggplot(x$table, ggplot2::aes(iteration, loss)) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::labs(y = "Mean native hinge loss", x = "Completed boosting update (starts at 1)")
  } else {
    ggplot2::ggplot(x$table, ggplot2::aes(iteration, coefficient, colour = feature)) +
      ggplot2::geom_line() + ggplot2::labs(y = "Coefficient in original margin equation",
        x = "Boosting update (zero = initial equation)", colour = "Term")
  }
  p + ggplot2::labs(title = x$measure, subtitle = x$learner, caption = x$interpretation) + sb_viz_theme()
}
