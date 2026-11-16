# ---- model-akritas ----

# ---- akritas-survival-adapter ----

sb_akritas_learner_ids <- function() "akritas_survival"

sb_akritas_params <- function() list(lambda = 0.5)

sb_require_akritas_dependencies <- function() {
  for (package in c("survivalmodels", "distr6")) if (!requireNamespace(package, quietly = TRUE))
    sb_abort(paste0("Install `", package, "` for native Akritas prediction."),
      "stormbreaker_missing_dependency")
}

sb_akritas_distribution <- function(basis) {
  if (ncol(basis) == 1L) distr6::Empirical$new(basis) else distr6::EmpiricalMV$new(basis)
}

sb_akritas_empirical_cdf <- function(basis, x) {
  if (ncol(basis) == 1L) {
    return(as.numeric(sb_akritas_distribution(basis)$cdf(data = x)))
  }
  # Same C routine as distr6::EmpiricalMV$.cdf, without reconstructing R6
  # support/parameter objects for already validated numeric matrices.
  native <- utils::getFromNamespace("C_EmpiricalMVCdf", "distr6")
  as.numeric(native(x = x, data = basis))
}

sb_fit_akritas_survival <- function(x, time, event, learner, params) {
  sb_require_akritas_dependencies()
  inputs <- sb_validate_classical_survival_inputs(x, time, event, learner)
  basis <- unname(as.matrix(inputs$x))
  # The native C engine expects decreasing follow-up time. Keep the empirical
  # distribution basis in input order, and reorder only the risk-set inputs.
  ord <- order(inputs$time, decreasing = TRUE)
  list(learner = learner, features = names(inputs$x), params = params,
    basis = basis, truth = unname(cbind(inputs$time, inputs$event)[ord, , drop = FALSE]),
    cdf_train = sb_akritas_empirical_cdf(basis, basis)[ord],
    grid = sort(unique(inputs$time)))
}

sb_akritas_curves <- function(model, x, times = model$grid) {
  sb_require_akritas_dependencies()
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) sb_abort(
    "Prediction columns differ from the fitted Akritas predictors.", "stormbreaker_learner_schema_mismatch")
  # survivalmodels 0.1.191 retains this native engine but removed the exported
  # R wrapper. Use the installed engine, not a replacement nearest-neighbour fit.
  native <- utils::getFromNamespace("C_Akritas", "survivalmodels")
  cdf <- sb_akritas_empirical_cdf(model$basis, unname(as.matrix(x)))
  distinct_cdf <- unique(cdf)
  # The native query input is this scalar CDF, not the predictor row. Compute
  # identical inputs once within this call, then restore patient order exactly.
  curve <- native(model$truth, as.numeric(times), model$grid,
    model$cdf_train, distinct_cdf,
    model$params$lambda)
  curve <- curve[match(cdf, distinct_cdf), , drop = FALSE]
  if (!is.matrix(curve) || !identical(dim(curve), c(nrow(x), length(times))) ||
      any(!is.finite(curve)) || any(curve < 0 | curve > 1))
    sb_abort("Native Akritas returned invalid survival curves.", "stormbreaker_invalid_akritas_prediction")
  curve
}

sb_predict_akritas_risk <- function(model, x) {
  curve <- sb_akritas_curves(model, x)
  sb_akritas_risk_from_curves(curve, model$grid)
}

sb_akritas_risk_from_curves <- function(curve, grid) {
  # Integral of the right-continuous survival step function up to max analysis
  # follow-up, not the historical wrapper's incomplete event-mass mean.
  before <- cbind(1, curve)[, seq_along(grid), drop = FALSE]
  -as.numeric(before %*% diff(c(0, grid)))
}

sb_akritas_direct_source <- function() "native_akritas_conditional_product_limit"

sb_akritas_survival_at_horizon <- function(model, x_analysis, time, event,
                                         x_assessment, prediction_horizon, ...) {
  horizon <- sb_validate_prediction_horizon(prediction_horizon, time)
  list(prediction_horizon = horizon, source = sb_akritas_direct_source(),
    ties = "native_grouped_event_times",
    probability = as.numeric(sb_akritas_curves(model, x_assessment, horizon)))
}

sb_extract_akritas_state <- function(model, context = list()) {
  state <- model
  state$training_data_retention <- "predictor_basis_and_paired_outcomes_required_for_conditional_product_limit"
  state$survival_probability <- if (!is.null(context$prediction_horizon)) list(
    prediction_horizon = sb_validate_prediction_horizon(context$prediction_horizon, context$time),
    time_unit = context$time_unit, source = sb_akritas_direct_source(), ties = "native_grouped_event_times") else NULL
  class(state) <- c("stormbreaker_akritas_state", "stormbreaker_learner_state")
  state
}

sb_predict_akritas_state <- function(state, x) {
  x <- sb_validate_learner_state_x(state, x, state$learner)
  times <- sort(unique(c(state$grid, state$survival_probability$prediction_horizon)))
  curves <- sb_akritas_curves(state, x, times)
  result <- list(risk = sb_akritas_risk_from_curves(
    curves[, match(state$grid, times), drop = FALSE], state$grid))
  if (!is.null(state$survival_probability)) {
    p <- state$survival_probability
    result$survival_probability <- as.numeric(curves[, match(p$prediction_horizon, times)])
    result$prediction_horizon <- p$prediction_horizon
    result$time_unit <- p$time_unit
  }
  result
}

sb_akritas_learner_adapter <- function(learner, params = sb_akritas_params()) {
  list(id = learner, params = params,
    fit = function(x, time, event) sb_fit_akritas_survival(x, time, event, learner, params),
    predict = sb_predict_akritas_risk, predict_survival = sb_akritas_survival_at_horizon,
    extract_state = sb_extract_akritas_state, predict_state = sb_predict_akritas_state)
}

# ---- algorithm-akritas ----

# Internal native diagnostics used by the common plotting registry.
# Neighborhood counts use the saved empirical-CDF uniform kernel.
# Curves use the original conditional product-limit engine without refitting.
sb_akritas_diagnostics <- function(fit, type = "neighborhood", newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!identical(s$learner, "akritas_survival") || !is.character(type) ||
      length(type) != 1L || is.na(type) ||
      !type %in% c("neighborhood", "survival_curves")) {
    sb_abort("Choose neighborhood or survival_curves for an Akritas fit.",
      "stormbreaker_unsupported_diagnostic")
  }
  if (is.null(newdata)) {
    sb_abort("Supply explicit queries for Akritas diagnostics.",
      "stormbreaker_unavailable_diagnostic")
  }
  x <- sb_select_final_feature_frame(newdata, s$features, "Query data", "query fitted features")
  if (!nrow(x)) {
    sb_abort("Queries must contain at least one row.", "stormbreaker_invalid_diagnostic")
  }
  if (type == "neighborhood") {
    cdf <- sb_akritas_empirical_cdf(s$basis, unname(as.matrix(x)))
    counts <- lapply(cdf, function(value) {
      inside <- abs(value - s$cdf_train) <= s$params$lambda
      c(neighbors = sum(inside), events = sum(s$truth[inside, 2L]))
    })
    counts <- do.call(rbind, counts)
    tab <- data.frame(query = seq_len(nrow(x)), empirical_cdf = cdf,
      neighbors = as.integer(counts[, "neighbors"]), events = as.integer(counts[, "events"]))
    tab$censored <- tab$neighbors - tab$events
  } else {
    times <- sort(unique(c(0, s$grid)))
    curves <- sb_akritas_curves(s, x, times)
    tab <- data.frame(query = rep(seq_len(nrow(x)), each = length(times)),
      time = rep(times, nrow(x)), survival = as.numeric(t(curves)))
  }
  structure(list(type = paste0("akritas_", type), table = tab,
    learner = s$learner, pipeline_id = fit$pipeline$id, lambda = s$params$lambda,
    time_unit = fit$development$endpoint$time_unit,
    source = "Saved empirical-CDF uniform kernel and native Akritas conditional product-limit engine",
    interpretation = "Neighborhood counts describe training support. Curves are original model predictions for the supplied queries; no refitting or bandwidth search.",
    sensitive_training_data = TRUE), class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_akritas_diagnostic <- function(x) {
  empirical_cdf <- neighbors <- time <- survival <- query <- NULL
  if (x$type == "akritas_neighborhood") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(empirical_cdf, neighbors)) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.2) +
      ggplot2::labs(title = "Kernel neighborhood", x = "Query empirical CDF",
        y = "Training neighbors")
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(time, survival, group = query))
    if (length(unique(x$table$query)) <= 10L) {
      p <- p + ggplot2::geom_step(ggplot2::aes(colour = factor(query)), linewidth = 0.7) +
        ggplot2::labs(colour = "Query")
    } else {
      p <- p + ggplot2::geom_step(colour = "#2878A0", alpha = 0.25, linewidth = 0.4)
    }
    time_label <- if (is.null(x$time_unit)) "Time" else paste0("Time (", x$time_unit, ")")
    p <- p + ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = "Conditional survival", x = time_label,
        y = "Survival probability")
  }
  p + sb_viz_theme()
}

# ---- diagnostics-kernels ----

# ---- algorithm-gaussian-process ----

#' @rdname sb_knn_neighborhood
#' @export
sb_gaussian_process_diagnostics <- function(fit, type = "basis_weights", newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (s$learner != "gaussian_process_binary" || !is.character(type) ||
      length(type) != 1L || is.na(type) || !type %in% c("basis_weights", "query_contributions")) {
    sb_abort("Choose basis_weights or query_contributions for a Gaussian-process classifier.",
      "stormbreaker_unsupported_diagnostic")
  }
  if (!requireNamespace("kernlab", quietly = TRUE)) sb_abort(
    "Gaussian-process diagnostics require kernlab.", "stormbreaker_unavailable_diagnostic")
  e <- unserialize(s$engine_raw)
  if (!methods::is(e, "gausspr") || e@type != "classification" ||
      !identical(e@lev, c("0", "1")) || length(e@alpha) != 1L ||
      !identical(colnames(e@xmatrix), s$features)) {
    sb_abort("Original binary Gaussian-process model is invalid.", "stormbreaker_invalid_diagnostic")
  }
  indices <- e@alphaindex[[1L]]
  weights <- -as.numeric(e@alpha[[1L]])
  if (length(indices) != length(weights) || any(!is.finite(weights))) {
    sb_abort("Original Gaussian-process basis weights are invalid.", "stormbreaker_invalid_diagnostic")
  }
  summary <- NULL
  if (type == "basis_weights") {
    tab <- data.frame(basis = seq_along(weights), reference_row = indices, weight = weights)
  } else {
    if (is.null(newdata)) sb_abort("Supply explicit queries for kernel contributions.", "stormbreaker_unavailable_diagnostic")
    x <- as.matrix(sb_select_final_feature_frame(newdata, s$features, "Query data", "query fitted features"))
    if (!nrow(x)) sb_abort("Queries must contain at least one row.", "stormbreaker_invalid_diagnostic")
    native <- kernlab::predict(e, x, type = "probabilities")
    if (is.list(e@scaling)) {
      selected <- e@scaling$scaled
      x[, selected] <- scale(x[, selected, drop = FALSE],
        center = e@scaling$x.scale[["scaled:center"]], scale = e@scaling$x.scale[["scaled:scale"]])
    }
    k <- as.matrix(kernlab::kernelMatrix(e@kernelf, x, e@xmatrix[indices, , drop = FALSE]))
    terms <- sweep(k, 2L, weights, "*")
    tab <- data.frame(query = rep(seq_len(nrow(x)), each = length(weights)),
      basis = rep(seq_along(weights), nrow(x)), reference_row = rep(indices, nrow(x)),
      similarity = as.numeric(t(k)), contribution = as.numeric(t(terms)))
    summary <- data.frame(query = seq_len(nrow(x)), latent_score = rowSums(terms),
      probability = as.numeric(native[, "1"]))
  }
  structure(list(type = paste0("gaussian_process_", type), table = tab, summary = summary,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original kernlab gausspr basis coefficients and saved kernel",
    interpretation = "Class-1-oriented reference-row kernel terms, not feature effects or SHAP. No uncertainty band, refit or calibration. Reference indices may be sensitive."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_gaussian_process_diagnostic <- function(x) {
  basis <- weight <- contribution <- NULL
  if (x$type == "gaussian_process_basis_weights") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(basis, weight)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey65", linewidth = .3) +
      ggplot2::geom_point(colour = "#2878A0", size = 1.3) +
      ggplot2::labs(title = "Gaussian-process basis weights", x = "Native basis", y = "Class-1 weight")
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(basis, contribution)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey65", linewidth = .3) +
      ggplot2::geom_point(colour = "#2878A0", size = 1.3) +
      ggplot2::facet_wrap(~query, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "Kernel contributions", x = "Native basis", y = "Class-1 latent contribution")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}

# ---- algorithm-lssvm ----

sb_lssvm_diagnostic <- function(fit, type) {
  if (!is.character(type) || length(type) != 1L || is.na(type)) sb_abort(
    "Choose decision_distribution for LS-SVM.", "stormbreaker_unsupported_diagnostic")
  if (type != "decision_distribution") sb_abort(
    "Radial LS-SVM retains a dense kernel basis, not a sparse support-vector set or input-feature coefficients.",
    "stormbreaker_unavailable_diagnostic")
  a <- fit$diagnostic_artifacts
  if (is.null(a$predictors)) sb_abort("Development predictors were not captured; enable capture_diagnostics at final fitting.",
    "stormbreaker_missing_diagnostic_artifacts")
  e <- unserialize(fit$learner_state$engine_raw)
  features <- fit$learner_state$features
  if (!identical(names(a$predictors), features) || !identical(e$class_labels, c("0", "1")) ||
      !identical(names(e$center), features) || !identical(names(e$scale), features) ||
      any(!is.finite(e$scale)) || any(e$scale <= 0) || !is.matrix(e$alpha) || ncol(e$alpha) != 2L ||
      !is.matrix(e$scaled_training) || nrow(e$scaled_training) != nrow(e$alpha) ||
      ncol(e$scaled_training) != length(features) || length(e$intercept) != 2L ||
      length(e$sigma) != 1L || !is.finite(e$sigma) || e$sigma <= 0) sb_abort(
    "Original LS-SVM kernel basis or class mapping is invalid.", "stormbreaker_invalid_diagnostic")
  z <- sweep(sweep(data.matrix(a$predictors), 2L, e$center, "-"), 2L, e$scale, "/")
  # Same saved kernel expansion, not a refit or inverse probability transform.
  distance <- outer(rowSums(z^2), rowSums(e$scaled_training^2), "+") - 2 * tcrossprod(z, e$scaled_training)
  scores <- sweep(exp(-e$sigma * pmax(distance, 0)) %*% e$alpha, 2L, e$intercept, "+")
  outcome <- a$endpoint$outcome
  if (any(!is.finite(scores)) || length(outcome) != nrow(scores) || anyNA(outcome) ||
      !all(outcome %in% c(0, 1))) sb_abort("Original LS-SVM scores or captured outcomes are invalid.",
    "stormbreaker_invalid_diagnostic")
  tab <- data.frame(observation = seq_len(nrow(scores)), class = as.character(outcome),
    decision = as.numeric(scores[, 2L] - scores[, 1L]))
  structure(list(type = "svm_decision_distribution", table = tab, native_table = unname(scores),
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id, sensitive_training_data = TRUE,
    interpretation = "Original development class-score contrast (class 1 minus class 0); not held-out performance.\nZero is the class boundary. The package uses plogis(contrast), not a fitted calibration; scores are not probabilities.",
    source = "original saved LS-SVM kernel expansion; no fitting"), class = "stormbreaker_algorithm_diagnostic")
}

# ---- algorithm-survival-svm ----

#' @rdname sb_knn_neighborhood
#' @export
sb_survival_svm_diagnostics <- function(fit, type = "score_distribution") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!fit$pipeline$learner %in% sb_survivalsvm_learner_ids()) sb_abort(
    "This diagnostic requires a survival-SVM learner.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("score_distribution", "basis_weights", "coefficients")) sb_abort(
    "Choose score_distribution, basis_weights or coefficients.", "stormbreaker_unsupported_diagnostic")
  a <- fit$diagnostic_artifacts
  model <- a$learner
  if (is.null(model)) sb_abort("Original survival-SVM artifacts were not captured; enable capture_diagnostics at final fitting.",
    "stormbreaker_missing_diagnostic_artifacts")
  features <- model$features
  formulation <- model$params$type
  kernel <- model$params$kernel
  native <- sb_survivalsvm_native_components(model$engine, formulation, length(features))
  intercept <- NULL
  if (type == "basis_weights") {
    tab <- data.frame(basis_row = seq_along(native$row_weight), value = native$row_weight)
    label <- if (formulation == "regression") "Native support-basis weights" else "Native combined kernel-row weights"
    interpretation <- "Signed weights of the original kernel expansion; not feature or patient importance.\nRanking/hybrid basis rows are not an e1071-style sparse support-vector set. Zero weights are retained."
  } else if (type == "coefficients") {
    if (kernel != "lin_kernel") sb_abort("This nonlinear kernel has no single exact input-feature coefficient vector.",
      "stormbreaker_unavailable_diagnostic")
    w <- as.numeric(crossprod(native$basis, native$row_weight))
    beta <- w / model$audit$scale
    intercept <- native$intercept - sum(w * model$audit$center / model$audit$scale)
    tab <- data.frame(feature = features, value = as.numeric(beta))
    label <- "Native survival-SVM original-unit weights"
    interpretation <- "Original native score equation before Cox calibration; intercept is returned separately.\nNot calibrated log-hazard coefficients, hazard ratios, uncertainty or feature importance."
  } else {
    if (!identical(names(a$predictors), features)) sb_abort("Captured predictor order does not match survival-SVM features.",
      "stormbreaker_invalid_diagnostic")
    if (!requireNamespace("survivalsvm", quietly = TRUE)) sb_abort("Package survivalsvm is required for native score evaluation.",
      "stormbreaker_missing_dependency")
    x <- as.data.frame(sweep(sweep(as.matrix(a$predictors), 2L, model$audit$center, "-"),
      2L, model$audit$scale, "/"))
    score <- as.numeric(stats::predict(model$engine, newdata = x)$predicted)
    if (length(score) != nrow(x)) sb_abort("Native survival-SVM score length is invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(observation = seq_along(score), value = score)
    label <- if (formulation == "regression") "Native predicted-time score" else "Native ranking / hybrid score"
    interpretation <- "Original development-score distribution before Cox calibration; not held-out discrimination.\nScores are not survival probabilities or calibrated log hazards. Censoring is not treated as a binary class."
  }
  if (any(!is.finite(tab$value)) || (!is.null(intercept) && !is.finite(intercept))) sb_abort(
    "Original survival-SVM diagnostic values are nonfinite.", "stormbreaker_invalid_diagnostic")
  structure(list(type = paste0("survival_svm_", type), table = tab, intercept = unname(intercept),
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    formulation = formulation, kernel = kernel, label = label,
    calibration_coefficient = model$audit$calibration_coefficient,
    sensitive_training_data = type != "coefficients", interpretation = interpretation,
    source = "original captured survival-SVM; no fitting"), class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_survival_svm_diagnostic <- function(x) {
  feature <- value <- NULL
  tab <- x$table
  if (x$type == "survival_svm_coefficients") {
    tab$feature <- factor(tab$feature, levels = rev(tab$feature))
    p <- ggplot2::ggplot(tab, ggplot2::aes(value, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
      ggplot2::geom_point(colour = "#2878A0") + ggplot2::labs(y = NULL)
  } else {
    p <- ggplot2::ggplot(tab, ggplot2::aes(value)) +
      ggplot2::stat_ecdf(geom = "step", pad = FALSE, colour = "#2878A0") +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(y = "Cumulative fraction")
    if (length(unique(tab$value)) == 1L) p <- p +
      ggplot2::geom_point(data = data.frame(value = tab$value[1L]), ggplot2::aes(y = 1), colour = "#2878A0")
  }
  p + ggplot2::labs(title = x$label, subtitle = x$learner, x = x$label,
    caption = x$interpretation) + sb_viz_theme()
}

# ---- algorithm-svm ----

#' @rdname sb_knn_neighborhood
#' @export
sb_svm_diagnostics <- function(fit, type = "support_counts") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (fit$pipeline$learner == "lssvm_radial_binary") return(sb_lssvm_diagnostic(fit, type))
  if (fit$pipeline$learner %in% sb_linear_margin_diagnostic_ids()) return(sb_linear_margin_diagnostic(fit, type))
  ids <- c("svm_linear_binary", "svm_radial_binary",
    "svm_radial_weighted_binary", "svm_polynomial_binary", "svm_sigmoid_binary")
  id <- fit$pipeline$learner
  if (!id %in% ids) sb_abort("This diagnostic requires a supported binary LIBSVM learner.",
    "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("support_counts", "decision_distribution", "coefficients")) sb_abort(
    "Choose support_counts, decision_distribution or coefficients.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(fit$learner_state$engine_raw)
  features <- fit$learner_state$features
  if (!inherits(e, "svm") || e$nclasses != 2L || !setequal(e$levels, c("0", "1")) ||
      !setequal(e$labels, 1:2) || length(e$nSV) != 2L || any(!is.finite(e$nSV)) ||
      any(e$nSV < 0) || any(e$nSV != floor(e$nSV)) || sum(e$nSV) != nrow(e$SV) ||
      !identical(colnames(e$SV), features)) sb_abort(
    "Original LIBSVM structure or class mapping is invalid.", "stormbreaker_invalid_diagnostic")
  classes <- e$levels[e$labels]
  orientation <- if (classes[1L] == "1") 1 else -1
  sensitive <- FALSE
  intercept <- NULL
  if (type == "support_counts") {
    tab <- data.frame(class = classes, support_vectors = as.integer(e$nSV))
    native <- list(nSV = e$nSV, labels = e$labels, levels = e$levels)
    interpretation <- "Original support-vector counts by encoded class; not class prevalence or feature importance."
  } else if (type == "coefficients") {
    if (id != "svm_linear_binary" || e$kernel != 0) sb_abort(
      "Nonlinear kernels have no single exact input-feature coefficient vector.", "stormbreaker_unavailable_diagnostic")
    w <- as.numeric(crossprod(e$coefs, e$SV))
    center <- rep(0, length(features))
    scale <- rep(1, length(features))
    if (any(e$scaled)) {
      center[e$scaled] <- e$x.scale[["scaled:center"]]
      scale[e$scaled] <- e$x.scale[["scaled:scale"]]
    }
    if (length(w) != length(features) || any(!is.finite(w)) || any(!is.finite(scale)) ||
        any(scale <= 0) || any(!is.finite(center)) || length(e$rho) != 1L || !is.finite(e$rho)) sb_abort(
      "Original linear LIBSVM scaling or coefficients are invalid.", "stormbreaker_invalid_diagnostic")
    beta <- orientation * w / scale
    intercept <- orientation * (-e$rho - sum(w * center / scale))
    tab <- data.frame(feature = features, coefficient = beta)
    native <- list(coefs = e$coefs, rho = e$rho, scaled_weights = w, scale = scale, center = center)
    interpretation <- "Exact decision-score weights in original predictor units, oriented toward class 1.\nIntercept is returned separately; not log-odds, odds ratios, uncertainty or feature importance."
  } else {
    a <- fit$diagnostic_artifacts
    if (is.null(a$predictors)) sb_abort("Development predictors were not captured; enable capture_diagnostics at final fitting.",
      "stormbreaker_missing_diagnostic_artifacts")
    if (!identical(names(a$predictors), features)) sb_abort(
      "Captured predictor order does not match the original SVM.", "stormbreaker_invalid_diagnostic")
    if (!requireNamespace("e1071", quietly = TRUE)) sb_abort("Package e1071 is required to evaluate native decision values.",
      "stormbreaker_missing_dependency")
    pred <- stats::predict(e, newdata = data.matrix(a$predictors), decision.values = TRUE, probability = FALSE)
    raw <- attr(pred, "decision.values")
    outcome <- a$endpoint$outcome
    if (!is.matrix(raw) || ncol(raw) != 1L || nrow(raw) != nrow(a$predictors) || any(!is.finite(raw)) ||
        length(outcome) != nrow(raw) || anyNA(outcome) || !all(outcome %in% c(0, 1))) sb_abort(
      "Native decision values or captured outcomes are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(observation = seq_len(nrow(raw)), class = as.character(outcome),
      decision = orientation * as.numeric(raw))
    native <- unname(raw)
    sensitive <- TRUE
    interpretation <- "Original development-sample decision scores, oriented toward class 1; not held-out performance.\nZero: native classification boundary, not necessarily calibrated probability 0.5. Scores are not probabilities."
  }
  structure(list(type = paste0("svm_", type), table = tab, native_table = native,
    intercept = intercept, orientation = orientation, learner = id, pipeline_id = fit$pipeline$id,
    sensitive_training_data = sensitive, interpretation = interpretation,
    source = "original binary LIBSVM model; no fitting"), class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_svm_diagnostic <- function(x) {
  class <- support_vectors <- feature <- coefficient <- decision <- NULL
  tab <- x$table
  if (x$type == "svm_support_counts") {
    tab$class <- factor(tab$class, levels = c("0", "1"))
    p <- ggplot2::ggplot(tab, ggplot2::aes(class, support_vectors)) +
      ggplot2::geom_col(fill = "#2878A0", width = .5) +
      ggplot2::scale_y_continuous(breaks = function(limits) unique(round(pretty(limits)))) +
      ggplot2::labs(title = "Support-vector counts", x = "Encoded class", y = "Number of support vectors")
  } else if (x$type == "svm_coefficients") {
    tab$feature <- factor(tab$feature, levels = rev(tab$feature))
    p <- ggplot2::ggplot(tab, ggplot2::aes(coefficient, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
      ggplot2::geom_point(colour = "#2878A0") +
      ggplot2::labs(title = "Original-unit decision weights", x = "Decision-score coefficient", y = NULL)
  } else {
    p <- ggplot2::ggplot(tab, ggplot2::aes(decision, colour = class)) +
      ggplot2::stat_ecdf(geom = "step", pad = FALSE) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey50", linetype = 2) +
      ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#D77935")) +
      ggplot2::labs(title = "Development decision-score distributions", x = "Decision score toward class 1",
        y = "Cumulative fraction within class", colour = "Observed class")
  }
  p + ggplot2::labs(subtitle = x$learner, caption = x$interpretation) + sb_viz_theme()
}

# ---- svm-probability-calibration ----

# The sigmoid optimizer below adapts LIBSVM's sigmoid_train algorithm.
# Copyright (c) 2000-2023 Chih-Chung Chang and Chih-Jen Lin.
# See inst/COPYRIGHTS for the complete BSD-3-Clause notice and disclaimer.

sb_svm_probability_learner_ids <- function() {
  c("svm_linear_binary", "svm_radial_binary",
    "svm_radial_weighted_binary", "svm_polynomial_binary",
    "svm_sigmoid_binary")
}

sb_make_svm_calibration_folds <- function(outcome, seed, group_ids = NULL) {
  n <- length(outcome)
  class_order <- unlist(lapply(unique(outcome), function(value) {
    which(outcome == value)
  }), use.names = FALSE)
  # LIBSVM first groups by class and then shuffles with one uniform per row.
  # Keep that ordering, including the subproblem training order, for singleton
  # groups. R's sample.int() does not reproduce the native shuffle.
  shuffled <- withr::with_seed(seed, {
    perm <- seq_len(n)
    for (i in seq_len(n)) {
      j <- i + floor(stats::runif(1L) * (n - i + 1L))
      tmp <- perm[i]
      perm[i] <- perm[j]
      perm[j] <- tmp
    }
    class_order[perm]
  })
  if (!is.null(group_ids)) {
    group_ids <- sb_validate_group_ids(group_ids, n)
  }
  grouped <- !is.null(group_ids) && anyDuplicated(group_ids) > 0L
  if (grouped) {
    fold_id <- sb_make_group_stratified_fold_id(
      outcome, group_ids, 5L, seed, "stormbreaker_infeasible_inner_folds",
      require_both_classes = FALSE
    )
  } else {
    fold_id <- integer(n)
    for (fold in seq_len(5L)) {
      begin <- floor((fold - 1L) * n / 5L)
      end <- floor(fold * n / 5L)
      if (end > begin) fold_id[shuffled[seq.int(begin + 1L, end)]] <- fold
    }
  }
  list(fold_id = fold_id, order = shuffled, class_order = class_order,
       source = if (grouped) "group_stratified" else "libsvm_row_shuffle")
}

sb_fit_libsvm_sigmoid <- function(decision, positive) {
  prior1 <- sum(positive)
  prior0 <- length(positive) - prior1
  target <- ifelse(positive, (prior1 + 1) / (prior1 + 2), 1 / (prior0 + 2))
  loss <- function(A, B) {
    z <- decision * A + B
    sum(ifelse(z >= 0, target * z, (target - 1) * z) +
          log1p(exp(-abs(z))))
  }
  A <- 0
  B <- log((prior0 + 1) / (prior1 + 1))
  value <- loss(A, B)
  termination <- "max_iterations"
  for (iteration in seq_len(100L)) {
    z <- decision * A + B
    small <- exp(-abs(z))
    p <- ifelse(z >= 0, small / (1 + small), 1 / (1 + small))
    q <- ifelse(z >= 0, 1 / (1 + small), small / (1 + small))
    d2 <- p * q
    h11 <- 1e-12 + sum(decision^2 * d2)
    h22 <- 1e-12 + sum(d2)
    h21 <- sum(decision * d2)
    g1 <- sum(decision * (target - p))
    g2 <- sum(target - p)
    if (abs(g1) < 1e-5 && abs(g2) < 1e-5) {
      termination <- "gradient_tolerance"
      break
    }
    determinant <- h11 * h22 - h21^2
    dA <- -(h22 * g1 - h21 * g2) / determinant
    dB <- -(-h21 * g1 + h11 * g2) / determinant
    directional <- g1 * dA + g2 * dB
    step <- 1
    while (step >= 1e-10) {
      newA <- A + step * dA
      newB <- B + step * dB
      newvalue <- loss(newA, newB)
      if (newvalue < value + 0.0001 * step * directional) {
        A <- newA
        B <- newB
        value <- newvalue
        break
      }
      step <- step / 2
    }
    if (step < 1e-10) {
      termination <- "line_search_failed"
      break
    }
  }
  list(A = A, B = B, termination = termination, iterations = iteration,
       negative_log_likelihood = value,
       positive_target = (prior1 + 1) / (prior1 + 2),
       negative_target = 1 / (prior0 + 2))
}

sb_calibrate_svm_probability <- function(x, outcome, learner, params,
                                          patient_ids, group_ids = NULL,
                                          selector = NULL) {
  inputs <- sb_validate_classical_binary_inputs(x, outcome, learner)
  sb_require_classical_binary_dependency(learner)
  x <- inputs$x
  outcome <- inputs$outcome
  patient_ids <- as.character(patient_ids)
  if (length(patient_ids) != nrow(x) || anyNA(patient_ids) ||
      any(!nzchar(trimws(patient_ids))) || anyDuplicated(patient_ids)) {
    sb_abort("SVM calibration requires a unique analysis ID for each row.",
             "stormbreaker_bad_runtime_patient_ids")
  }
  if (!is.null(group_ids)) {
    group_ids <- sb_validate_group_ids(group_ids, nrow(x))
  }
  folds <- sb_make_svm_calibration_folds(outcome, params$seed, group_ids)
  positive <- outcome[[1L]]
  decision <- rep(NA_real_, nrow(x))
  records <- vector("list", 5L)
  warnings <- character()
  fit_params <- params
  fit_params$probability <- FALSE
  for (fold in seq_len(5L)) {
    train <- folds$order[folds$fold_id[folds$order] != fold]
    assess <- folds$order[folds$fold_id[folds$order] == fold]
    # Native binary subproblems use +1/-1, with +1 first after class grouping.
    train <- c(train[outcome[train] == positive],
               train[outcome[train] != positive])
    features <- colnames(x)
    center <- scale <- class_weights <- NULL
    fit_warnings <- character()
    if (!length(assess)) {
      fit_status <- "empty_assessment"
    } else if (length(unique(outcome[train])) < 2L) {
      # Preserve LIBSVM's defined constant margin for single-class subproblems.
      decision[assess] <- if (!length(train)) 0 else {
        if (outcome[train[[1L]]] == positive) 1 else -1
      }
      features <- character()
      fit_status <- "single_class_constant"
    } else {
      train_x <- x[train, , drop = FALSE]
      assess_x <- x[assess, , drop = FALSE]
      model <- withCallingHandlers({
        if (!is.null(selector)) {
          selected <- sb_fit_runtime_selector(selector, train_x,
              list(task_type = "binary", outcome = outcome[train], group_ids = group_ids[train]))
          train_x <- selector$apply(selected, train_x)
          assess_x <- selector$apply(selected, assess_x)
        }
        sb_fit_classical_binary(train_x, outcome[train], learner,
                                       fit_params)
      }, warning = function(warning) {
        fit_warnings <<- c(fit_warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      })
      features <- model$features
      predicted <- stats::predict(model$engine, as.matrix(assess_x),
                                   decision.values = TRUE)
      native_positive <- as.character(model$engine$levels[
        model$engine$labels[[1L]]])
      direction <- if (native_positive == as.character(positive)) 1 else -1
      decision[assess] <- direction * as.numeric(attr(predicted,
                                                      "decision.values"))
      center <- model$engine$x.scale[["scaled:center"]]
      scale <- model$engine$x.scale[["scaled:scale"]]
      if (identical(params$class_weights, "inverse_frequency_normalized")) {
        counts <- table(factor(outcome[train], levels = c(0L, 1L)))
        class_weights <- stats::setNames(length(train) / (2 * as.numeric(counts)),
                                         names(counts))
      }
      fit_status <- "success"
    }
    records[[fold]] <- list(
      fold = fold, analysis_patient_ids = patient_ids[train],
      assessment_patient_ids = patient_ids[assess], selected_features = features,
      status = fit_status, scale_center = center, scale_scale = scale,
      class_weights = class_weights, warnings = fit_warnings
    )
    warnings <- c(warnings, if (length(fit_warnings)) {
      paste0("Calibration fold ", fold, ": ", fit_warnings)
    })
  }
  if (any(!is.finite(decision))) {
    sb_abort("SVM calibration did not produce finite out-of-fold margins.",
             "stormbreaker_svm_calibration_prediction")
  }
  # Retain LIBSVM's class-grouped order for its optimizer, not assessment order.
  sigmoid <- sb_fit_libsvm_sigmoid(decision[folds$class_order],
                                   outcome[folds$class_order] == positive)
  if (!identical(sigmoid$termination, "gradient_tolerance")) {
    warnings <- c(warnings, paste0("LIBSVM sigmoid: ", sigmoid$termination))
  }
  list(
    lambda = NA_real_, alpha = NA_real_, selection_rule = "libsvm_platt",
    measure = "smoothed_binomial_likelihood", inner_folds = 5L,
    seed = params$seed, fold_id = stats::setNames(folds$fold_id, patient_ids),
    alignment = "decision_value",
    inner_selection_method = if (is.null(selector)) "none" else {
      "svm_fold_local_probability_calibration"
    },
    inner_selections = records, curve = NULL, warnings = warnings,
    standardize = params$standardize,
    calibration = c(list(
      method = "libsvm_platt_oof", positive_class = as.character(positive),
      fold_source = folds$source,
      group_ids = stats::setNames(if (is.null(group_ids)) patient_ids else {
        group_ids
      }, patient_ids),
      out_of_fold = data.frame(patient_id = patient_ids, outcome = outcome,
        fold = folds$fold_id, decision = decision, stringsAsFactors = FALSE)
    ), sigmoid)
  )
}

sb_svm_calibrated_runtime_learner <- function(learner, calibration) {
  id <- learner$id
  params <- learner$params
  learner$fit_context <- calibration
  learner$fit <- function(x, outcome) {
    fit_params <- params
    fit_params$probability <- FALSE
    model <- sb_fit_classical_binary(x, outcome, id, fit_params)
    positive <- as.character(model$engine$levels[model$engine$labels[[1L]]])
    model$engine$probA <- calibration$A
    model$engine$probB <- if (positive == calibration$positive_class) {
      calibration$B
    } else -calibration$B
    model$engine$compprob <- TRUE
    model$params <- params
    model
  }
  learner
}
