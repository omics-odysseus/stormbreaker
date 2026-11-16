# ---- multiclass-workflow ----

# ---- multiclass-evaluation ----

sb_validate_multiclass_probabilities <- function(probability, classes) {
  if (!is.matrix(probability) || !is.numeric(probability) || !nrow(probability) ||
      !identical(colnames(probability), classes) || any(!is.finite(probability)) ||
      any(probability < 0 | probability > 1) || any(abs(rowSums(probability) - 1) > 1e-6)) sb_abort(
    "Supply finite class probabilities summing to one, in declared class order.",
    "stormbreaker_bad_multiclass_probabilities")
  invisible(probability)
}

#' @rdname multi_run
#' @export
sb_multiclass_metrics <- function(truth, probability) {
  classes <- colnames(probability)
  if (is.null(classes) || length(classes) < 3L || anyNA(classes) ||
      anyDuplicated(classes) || any(!nzchar(classes))) sb_abort(
    "Probability columns must identify at least three distinct classes.",
    "stormbreaker_bad_multiclass_probabilities")
  sb_validate_multiclass_probabilities(probability, classes)
  if (length(truth) != nrow(probability) || anyNA(truth) ||
      any(!as.character(truth) %in% classes)) sb_abort(
    "Truth must contain known classes, aligned with probability rows.",
    "stormbreaker_bad_multiclass_outcome")
  truth <- factor(truth, levels = classes)
  predicted <- factor(classes[max.col(probability, ties.method = "first")], levels = classes)
  counts <- table(truth = truth, predicted = predicted)
  tp <- diag(counts)
  actual <- rowSums(counts)
  guessed <- colSums(counts)
  per_class <- data.frame(class = classes, n = as.integer(actual),
    precision = ifelse(guessed > 0, tp / guessed, NA_real_),
    recall = ifelse(actual > 0, tp / actual, NA_real_),
    f1 = ifelse(actual > 0, 2 * tp / (actual + guessed), NA_real_),
    auc = vapply(seq_along(classes), function(k)
      sb_binary_auc_stratum(as.integer(truth == classes[k]), probability[, k])$estimate,
      numeric(1)), row.names = NULL)
  target <- outer(as.character(truth), classes, `==`)
  list(overall = data.frame(accuracy = mean(truth == predicted),
    balanced_accuracy = mean(per_class$recall), macro_f1 = mean(per_class$f1),
    macro_auc = mean(per_class$auc),
    log_loss = -mean(log(pmax(probability[cbind(seq_along(truth), as.integer(truth))], 1e-15))),
    brier = mean(rowSums((probability - target)^2))),
    per_class = per_class, confusion = counts, predicted = predicted)
}

#' @rdname multi_run
#' @export
sb_compare_multiclass <- function(task, methods = c("lasso", "random_forest", "svm_radial"), folds = 5L, seed = 1L,
                                   inner_folds = 5L, inner_seed = 1701L) {
  task <- sb_validate_multiclass_task(task)
  methods <- sb_method_names(methods)
  catalog <- sb_method_catalog("multiclass")
  catalog$supported[catalog$method %in% sb_multiclass_standalone_selectors()] <- FALSE
  if (anyDuplicated(methods) || any(!methods %in% catalog$method[catalog$supported])) sb_abort(
    "Choose distinct supported multiclass methods.", "stormbreaker_bad_method")
  assignment <- sb_multiclass_folds(task, folds, seed)
  results <- lapply(methods, function(method) {
    probability <- matrix(NA_real_, task$n, length(task$classes),
      dimnames = list(as.character(task$data[[task$id_col]]), task$classes))
    fold_results <- vector("list", folds)
    for (fold in seq_len(folds)) {
      train <- which(assignment != fold)
      test <- which(assignment == fold)
      training <- sb_multiclass_task(task$data[train, , drop = FALSE], task$id_col,
        task$outcome_col, task$predictors, task$classes, task$group_col)
      fold_results[[fold]] <- tryCatch({
        fit <- sb_fit_multiclass(training, method, inner_folds, inner_seed,
          seed = seed + fold)
        probability[test, ] <- sb_predict_multiclass(fit, task$data[test, , drop = FALSE])
        list(status = "success", fold = fold, training_rows = train, assessment_rows = test,
          inner_fold_id = fit$fold_id, parameters = fit$parameters,
          warnings = fit$warnings,
          features = sb_multiclass_evidence(fit))
      }, error = function(e) list(status = "failed", fold = fold,
        error_message = conditionMessage(e)))
    }
    failed <- vapply(fold_results, function(x) x$status == "failed", logical(1))
    list(method = method, status = if (any(failed)) "failed" else "success",
      error_message = if (any(failed)) paste(vapply(fold_results[failed],
        function(x) paste0("Fold ", x$fold, ": ", x$error_message), character(1)), collapse = "; ") else "",
      probability = probability, folds = fold_results,
      metrics = if (any(failed)) NULL else sb_multiclass_metrics(task$data[[task$outcome_col]], probability))
  })
  names(results) <- methods
  comparison <- do.call(rbind, lapply(results, function(result) {
    row <- data.frame(method = result$method, status = result$status, reason = result$error_message)
    metrics <- if (result$status == "success") result$metrics$overall else data.frame(
      accuracy = NA_real_, balanced_accuracy = NA_real_, macro_f1 = NA_real_,
      macro_auc = NA_real_, log_loss = NA_real_, brier = NA_real_)
    cbind(row, metrics)
  }))
  rownames(comparison) <- NULL
  structure(list(comparison = comparison, results = results, assignment = assignment,
    ids = as.character(task$data[[task$id_col]]), truth = task$data[[task$outcome_col]],
    classes = task$classes, predictors = task$predictors, seed = seed,
    interpretation = "Shared outer-fold evaluation; multinomial tuning within training folds"),
    class = "stormbreaker_multiclass_comparison")
}

#' @rdname multi_run
#' @export
sb_plot_multiclass_confusion <- function(x) {
  if (!is.list(x) || !inherits(x$confusion, "table")) sb_abort(
    "Supply sb_multiclass_metrics() output.", "stormbreaker_bad_visualization_argument")
  truth <- predicted <- Freq <- NULL
  ggplot2::ggplot(as.data.frame(x$confusion), ggplot2::aes(predicted, truth, fill = Freq)) +
    ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = Freq)) +
    ggplot2::scale_fill_gradient(low = "#eff6fa", high = "#73aac5", guide = "none") +
    ggplot2::coord_equal() + sb_viz_theme() +
    ggplot2::labs(title = "Confusion matrix", x = "Predicted class", y = "True class")
}

# ---- multiclass-information ----

sb_multiclass_information_catalog <- function() {
  data.frame(method = c("mrmr", "jmi", "cmim", "disr", "jmim"),
    label = c("mRMR", "JMI", "CMIM", "DISR", "JMIM"), supported = TRUE,
    engine = "praznik", feature_result = "selection",
    selection_rule = "Native sequential information ranking; at most min(p, max(2, floor(n/log(n)))) genes")
}

sb_multiclass_standalone_selectors <- function() {
  c("boruta", "varselrf", "vsurf", sb_multiclass_information_catalog()$method,
    sb_multiclass_omics_catalog()$method)
}

sb_analyze_multiclass_information <- function(task, spec) {
  if (!requireNamespace("praznik", quietly = TRUE)) sb_abort(
    "Install praznik for this selector.", "stormbreaker_missing_dependency")
  x <- task$data[task$predictors]
  x[] <- lapply(x, as.numeric)
  y <- factor(task$data[[task$outcome_col]], levels = task$classes)
  size <- sb_high_dimensional_nlogn_size(nrow(x), ncol(x))
  args <- list(X = x, Y = y, k = size, threads = 1L)
  if (spec$method == "mrmr") args$positive <- FALSE
  native <- do.call(getExportedValue("praznik", toupper(spec$method)), args)
  selected <- as.integer(native$selection)
  if (length(selected) > size || anyNA(selected) || anyDuplicated(selected) ||
      any(!selected %in% seq_along(x)) || length(native$score) != length(selected) ||
      any(!is.finite(native$score))) sb_abort(
    "The information selector returned an invalid feature subset.",
    "stormbreaker_invalid_information_selector")
  tab <- data.frame(feature = names(x), score = NA_real_,
    rank = match(seq_along(x), selected), selected = seq_along(x) %in% selected)
  tab$score[selected] <- as.numeric(native$score)
  params <- list(screen_size = size, screen_size_rule = spec$selection_rule,
    threads = 1L, native_method = toupper(spec$method), selected_size = length(selected),
    discretization = "native_equal_width_bins_for_numeric_predictors")
  if (spec$method == "mrmr") params$positive <- FALSE
  fit <- structure(list(id = paste0(spec$method, "_multiclass_nlogn"),
    source_features = names(x), selected_features = names(x)[selected],
    scores = tab, params = params, classes = task$classes), class = "stormbreaker_selector_fit")
  feature <- score <- NULL
  shown <- tab[selected, ]
  shown$feature <- factor(shown$feature, levels = rev(shown$feature))
  plot <- ggplot2::ggplot(shown, ggplot2::aes(score, feature)) +
    ggplot2::geom_point(colour = "#2878A0") + ggplot2::theme_minimal() +
    ggplot2::labs(title = spec$label, x = "Selection score", y = NULL)
  list(status = "success", method = spec$method, label = spec$label,
    feature_result = "selection", selection_rule = spec$selection_rule,
    features = tab, class_features = NULL, evidence = tab, parameters = params,
    fit = fit, diagnostics = list(selection_scores = tab), plots = list(selection_scores = plot),
    notes = list(interpretation = paste("Full-input feature selection, not held-out evaluation.",
      "Scores are native sequential selection criteria, not patient risks or marginal gene effects.")),
    plot_status = data.frame(view = "selection_scores", status = "available", reason = ""))
}

# ---- multiclass-models ----

sb_fit_multiclass <- function(task, method, inner_folds, inner_seed,
                              alpha = 0.5, lambda_rule = "lambda.1se", seed = 1L) {
  catalog <- sb_method_catalog("multiclass")
  if (!method %in% catalog$method[catalog$supported] ||
      method %in% sb_multiclass_standalone_selectors()) sb_abort(
    paste(method, "is not implemented for multiclass tasks."),
    "stormbreaker_unsupported_method_task")
  x <- task$data[task$predictors]
  codes <- paste0("class", seq_along(task$classes))
  y <- factor(codes[match(as.character(task$data[[task$outcome_col]]), task$classes)], levels = codes)
  fit <- list(method = method, classes = task$classes, codes = codes,
    features = task$predictors, seed = seed, fold_id = NULL, lambda = NULL,
    warnings = character())
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  if (method %in% c("pre", "slope")) fit$fold_id <- sb_multiclass_folds(task, inner_folds, inner_seed)
  if (method == "abess") return(sb_fit_multiclass_abess(fit, x, y))
  if (method %in% sb_multiclass_native_catalog()$method) return(sb_fit_multiclass_native(fit, x, y))
  if (method %in% sb_multiclass_classical_catalog()$method) return(
    sb_fit_multiclass_classical(fit, x, y))
  if (method == "random_forest") {
    if (!requireNamespace("randomForest", quietly = TRUE)) sb_abort(
      "Install randomForest for this method.", "stormbreaker_missing_dependency")
    fit$engine <- withCallingHandlers(withr::with_seed(seed,
      randomForest::randomForest(x, y, ntree = 500L, importance = TRUE)),
      warning = record_warning)
    fit$parameters <- list(ntree = 500L, mtry = fit$engine$mtry, seed = seed)
  } else {
    if (!requireNamespace("glmnet", quietly = TRUE)) sb_abort(
      "Install glmnet for this method.", "stormbreaker_missing_dependency")
    if (!is.character(lambda_rule) || length(lambda_rule) != 1L || is.na(lambda_rule) ||
        !lambda_rule %in% c("lambda.1se", "lambda.min")) sb_abort(
      "Choose lambda.1se or lambda.min.", "stormbreaker_bad_method")
    inner_folds <- sb_check_positive_integer(inner_folds, "inner_folds", minimum = 3L)
    alpha <- switch(method, lasso = 1, ridge = 0, elastic_net = alpha)
    fit$fold_id <- sb_multiclass_folds(task, inner_folds, inner_seed)
    fit$engine <- withCallingHandlers(glmnet::cv.glmnet(as.matrix(x), y, family = "multinomial",
      alpha = alpha, foldid = fit$fold_id,
      type.measure = "deviance", type.multinomial = "ungrouped", standardize = TRUE),
      warning = record_warning)
    fit$lambda <- fit$engine[[lambda_rule]]
    fit$parameters <- list(alpha = alpha,
      lambda_rule = lambda_rule, lambda = fit$lambda, standardize = TRUE,
      type_multinomial = "ungrouped", inner_folds = inner_folds, inner_seed = inner_seed)
  }
  structure(fit, class = "stormbreaker_multiclass_fit")
}

#' @rdname multi_run
#' @export
sb_predict_multiclass <- function(fit, newdata) {
  if (!inherits(fit, "stormbreaker_multiclass_fit")) sb_abort(
    "Supply a saved multiclass fit.", "stormbreaker_bad_multiclass_fit")
  if (!is.data.frame(newdata) || !all(fit$features %in% names(newdata))) sb_abort(
    "New data must contain all fitted predictors.", "stormbreaker_missing_predictors")
  x <- sb_validate_predictor_frame(newdata[fit$features], "newdata")
  if (fit$method == "abess") {
    if (!requireNamespace("abess", quietly = TRUE)) sb_abort(
      "Install abess to predict with this saved fit.", "stormbreaker_missing_dependency")
    probability <- stats::predict(fit$engine, as.matrix(x),
      support.size = fit$engine$best.size, type = "response")[[1]]
  } else if (fit$method %in% sb_multiclass_native_catalog()$method) {
    probability <- sb_predict_multiclass_native(fit, x)
  } else if (fit$method %in% sb_multiclass_classical_catalog()$method) {
    probability <- sb_predict_multiclass_classical(fit, x)
  } else if (fit$method == "random_forest") {
    if (!requireNamespace("randomForest", quietly = TRUE)) sb_abort(
      "Install randomForest to predict with this saved fit.", "stormbreaker_missing_dependency")
    probability <- stats::predict(fit$engine, x, type = "prob")[, fit$codes, drop = FALSE]
  } else {
    if (!requireNamespace("glmnet", quietly = TRUE)) sb_abort(
      "Install glmnet to predict with this saved fit.", "stormbreaker_missing_dependency")
    array <- stats::predict(fit$engine, as.matrix(x), s = fit$lambda, type = "response")
    probability <- matrix(array[, fit$codes, 1L], nrow = nrow(x), ncol = length(fit$codes))
  }
  colnames(probability) <- fit$classes
  rownames(probability) <- rownames(newdata)
  sb_validate_multiclass_probabilities(probability, fit$classes)
  probability
}

sb_multiclass_evidence <- function(fit) {
  if (fit$method %in% sb_multiclass_native_catalog()$method) return(sb_multiclass_native_evidence(fit))
  if (fit$method == "abess") return(sb_multiclass_abess_evidence(fit))
  if (fit$method %in% sb_multiclass_classical_catalog()$method) return(
    sb_multiclass_classical_evidence(fit))
  if (fit$method == "random_forest") {
    native <- fit$engine$importance
    features <- data.frame(feature = fit$features,
      importance = unname(native[fit$features, "MeanDecreaseGini"]), selected = NA)
    features$rank <- rank(-features$importance, ties.method = "min")
    by_class <- do.call(rbind, lapply(seq_along(fit$classes), function(k) data.frame(
      feature = fit$features, class = fit$classes[k],
      importance = unname(native[fit$features, fit$codes[k]]), selected = NA)))
    by_class$rank <- stats::ave(-by_class$importance, by_class$class,
      FUN = function(x) rank(x, ties.method = "min"))
    return(list(features = features, by_class = by_class))
  }
  beta <- stats::coef(fit$engine, s = fit$lambda)
  by_class <- do.call(rbind, lapply(seq_along(fit$classes), function(k) {
    b <- as.matrix(beta[[fit$codes[k]]])[fit$features, 1L]
    data.frame(feature = fit$features, class = fit$classes[k],
      coefficient = unname(b), selected = unname(b != 0))
  }))
  features <- data.frame(feature = fit$features, selected = vapply(fit$features,
    function(gene) any(by_class$selected[by_class$feature == gene]), logical(1)))
  if (fit$method == "ridge") {
    features$selected <- NA
    by_class$selected <- NA
  }
  list(features = features, by_class = by_class)
}

sb_multiclass_plots <- function(fit, evidence, new_data = NULL) {
  if (fit$method %in% sb_multiclass_native_catalog()$method) return(sb_multiclass_native_plots(fit, new_data))
  if (fit$method == "abess") return(sb_multiclass_abess_plots(fit, evidence))
  if (fit$method %in% sb_multiclass_classical_catalog()$method) return(
    sb_multiclass_classical_plots(fit, evidence))
  feature <- importance <- class <- iteration <- error <- lambda <- coefficient <- NULL
  mean_deviance <- standard_error <- NULL
  if (fit$method == "random_forest") {
    tab <- evidence$features
    tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$importance)])
    importance_plot <- ggplot2::ggplot(tab, ggplot2::aes(importance, feature)) +
      ggplot2::geom_col(fill = "#2878A0") + sb_viz_theme() +
      ggplot2::labs(title = "Random forest: Importance", x = "Mean decrease in Gini", y = NULL)
    errors <- fit$engine$err.rate
    history <- data.frame(iteration = rep(seq_len(nrow(errors)), ncol(errors)),
      series = rep(seq_len(ncol(errors)), each = nrow(errors)), error = as.vector(errors))
    history$class <- factor(history$series, levels = seq_len(ncol(errors)),
      labels = c("Overall", paste0("Class: ", fit$classes)))
    oob <- ggplot2::ggplot(history, ggplot2::aes(iteration, error, colour = class)) +
      ggplot2::geom_line() + sb_viz_theme() +
      ggplot2::labs(title = "Random forest: OOB error", x = "Trees", y = "Error", colour = NULL)
    class_data <- evidence$by_class
    class_data$class <- factor(class_data$class, levels = fit$classes)
    class_plot <- ggplot2::ggplot(class_data,
      ggplot2::aes(importance, feature)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::facet_wrap(~class) + sb_viz_theme() +
      ggplot2::labs(title = "Random forest: Class importance",
        x = "OOB permutation accuracy decrease", y = NULL)
    return(list(plots = list(importance = importance_plot, oob_error = oob,
      class_importance = class_plot), tables = list(importance = tab,
        oob_error = history, class_importance = evidence$by_class)))
  }
  engine <- fit$engine$glmnet.fit
  path <- do.call(rbind, lapply(seq_along(fit$classes), function(k) {
    b <- as.matrix(engine$beta[[fit$codes[k]]])
    data.frame(feature = rep(rownames(b), ncol(b)),
      lambda = rep(engine$lambda, each = nrow(b)), coefficient = as.vector(b),
      class = fit$classes[k])
  }))
  path$class <- factor(path$class, levels = fit$classes)
  label <- switch(fit$method, lasso = "LASSO", ridge = "Ridge", elastic_net = "Elastic Net")
  path_plot <- ggplot2::ggplot(path, ggplot2::aes(lambda, coefficient, colour = feature)) +
    ggplot2::geom_line() + ggplot2::scale_x_log10() + ggplot2::facet_wrap(~class) +
    sb_viz_theme() + ggplot2::labs(title = paste0(label, ": Class coefficient paths"),
      x = "Lambda (log scale)", y = "Coefficient", colour = "Feature")
  cv <- data.frame(lambda = fit$engine$lambda, mean_deviance = fit$engine$cvm,
    standard_error = fit$engine$cvsd)
  curve <- ggplot2::ggplot(cv, ggplot2::aes(lambda, mean_deviance)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_deviance - standard_error,
      ymax = mean_deviance + standard_error), width = 0, colour = "grey65") +
    ggplot2::geom_line(colour = "#2878A0") + ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept = fit$lambda, linetype = "dashed") +
    sb_viz_theme() + ggplot2::labs(title = paste0(label, ": Cross-validation"),
      x = "Lambda (log scale)", y = "Multinomial deviance")
  list(plots = list(coefficient_path = path_plot, cv_curve = curve),
    tables = list(coefficient_path = path, cv_curve = cv))
}

sb_analyze_multiclass_method <- function(task, spec, inner_folds, inner_seed,
                                         elastic_net_alpha, lambda_rule) {
  if (spec$method == "lpc") return(sb_analyze_multiclass_lpc(task, spec))
  if (spec$method %in% sb_multiclass_information_catalog()$method) return(
    sb_analyze_multiclass_information(task, spec))
  if (spec$method %in% c("boruta", "varselrf", "vsurf")) return(
    sb_analyze_multiclass_selector(task, spec))
  fit <- sb_fit_multiclass(task, spec$method, inner_folds, inner_seed,
    elastic_net_alpha, lambda_rule)
  evidence <- sb_multiclass_evidence(fit)
  drawn <- tryCatch(sb_multiclass_plots(fit, evidence, task$data[task$predictors]), error = identity)
  failed <- inherits(drawn, "error")
  figure_notes <- list()
  if (!failed) {
    drawn$plots <- lapply(drawn$plots, sb_compact_plot)
    figure_notes <- lapply(drawn$plots, attr, which = "stormbreaker_notes")
  }
  list(status = "success", method = spec$method, label = spec$label,
    feature_result = spec$feature_result, selection_rule = spec$selection_rule,
    features = evidence$features, class_features = evidence$by_class,
    evidence = evidence, parameters = fit$parameters, fit = fit,
    plots = if (failed) list() else drawn$plots,
    diagnostics = if (failed) list() else drawn$tables,
    notes = c(figure_notes, list(interpretation = if (fit$method %in% sb_multiclass_native_catalog()$method ||
        fit$method == "ridge") paste("Full-input exploratory fit, not held-out evaluation.",
      "No gene-selection rule is defined. SVM uses native one-versus-one probability coupling.",
      "Support-vector counts are model structure, not accuracy or gene importance.",
      "Only diagnostics implemented for this engine are returned; no surrogate model is fitted.") else
      if (fit$method %in% sb_multiclass_classical_catalog()$method) paste(
      "Full-input exploration, not held-out evaluation. Class means and SD",
      "are feature distributions, not importance, confidence intervals or selected genes.",
      "Tree importance is native split importance; used_for_split separately",
      "identifies actual primary split variables. No fallback model is fitted.") else paste(
      "Full-input exploration, not held-out evaluation. Class coefficients are",
      "the native multinomial parameterization, not binary odds ratios.",
      "Overall selection is the union of nonzero coefficients across classes.",
      "Forest importance is ranking only; OOB class importance is unscaled",
      "permutation accuracy decrease, not Gini importance. Native OOB sampling",
      "is row-based, not patient-group-isolated; use grouped outer evaluation",
      "for repeated patients."))),
    plot_status = data.frame(view = if (failed) "multiclass_plots" else names(drawn$plots),
      status = if (failed) "error" else rep("available", length(drawn$plots)),
      reason = if (failed) conditionMessage(drawn) else rep("", length(drawn$plots))))
}

# ---- multiclass-native ----

sb_multiclass_native_catalog <- function() {
  basic <- data.frame(method = c("svm_linear", "svm_radial", "svm_radial_weighted",
    "svm_polynomial", "svm_sigmoid", "kernel_naive_bayes", "rda", "nnet", "knn"),
    label = c("Linear SVM", "Radial SVM", "Weighted radial SVM", "Polynomial SVM",
      "Sigmoid SVM", "Kernel naive Bayes", "Regularized discriminant analysis",
      "Neural network", "K-nearest neighbors"), supported = TRUE,
    engine = c(rep("e1071_svm", 5), "klaR_NaiveBayes", "klaR_rda", "nnet_softmax", "knn_votes"),
    feature_result = "model_features", selection_rule = "No standalone gene-selection rule")
  rbind(basic, sb_multiclass_tree_catalog(), sb_multiclass_discriminant_catalog(),
    sb_multiclass_boosting_catalog(), sb_multiclass_specialist_catalog(),
    sb_multiclass_multinomial_catalog())
}

sb_fit_multiclass_native <- function(fit, x, y) {
  if (fit$method %in% sb_multiclass_multinomial_catalog()$method) return(sb_fit_multiclass_multinomial(fit, x, y))
  if (fit$method %in% sb_multiclass_specialist_catalog()$method) return(sb_fit_multiclass_specialist(fit, x, y))
  if (fit$method %in% sb_multiclass_boosting_catalog()$method) return(sb_fit_multiclass_boosting(fit, x, y))
  if (fit$method %in% sb_multiclass_discriminant_catalog()$method) return(
    sb_fit_multiclass_discriminant(fit, x, y))
  if (fit$method %in% sb_multiclass_tree_catalog()$method) return(sb_fit_multiclass_tree(fit, x, y))
  method <- fit$method
  package <- if (startsWith(method, "svm_")) "e1071" else
    switch(method, kernel_naive_bayes = , rda = "klaR", nnet = "nnet", knn = NA_character_)
  if (!is.na(package) && !requireNamespace(package, quietly = TRUE)) sb_abort(
    paste("Install", package, "for this method."), "stormbreaker_missing_dependency")
  fit$package <- package
  params <- sb_classical_binary_params(paste0(method, "_binary"))
  fit$parameters <- params
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  if (startsWith(method, "svm_")) {
    weights <- if (method == "svm_radial_weighted") {
      counts <- table(y)
      stats::setNames(length(y) / (length(counts) * as.numeric(counts)), names(counts))
    } else NULL
    fit$engine <- withCallingHandlers(withr::with_seed(fit$seed,
      e1071::svm(x, y, type = "C-classification", kernel = params$kernel,
        cost = params$cost, gamma = 1 / ncol(x), degree = 3L, coef0 = 0,
        scale = TRUE, probability = TRUE, class.weights = weights)), warning = record_warning)
    fit$parameters$gamma <- 1 / ncol(x)
    fit$parameters$class_weights <- weights
    fit$parameters$multiclass_strategy <- "native_one_vs_one_probability_coupling"
    fit$parameters$probability_calibration <- "libsvm_internal_row_based"
  } else if (method == "kernel_naive_bayes") {
    fit$engine <- withCallingHandlers(klaR::NaiveBayes(x, y, usekernel = TRUE,
      fL = params$fL, adjust = params$adjust), warning = record_warning)
  } else if (method == "rda") {
    fit$engine <- withCallingHandlers(klaR::rda(x, y, gamma = params$gamma,
      lambda = params$lambda, estimate.error = FALSE), warning = record_warning)
  } else {
    fit$center <- colMeans(x)
    fit$scale <- vapply(x, stats::sd, numeric(1))
    fit$scale[fit$scale == 0] <- 1
    z <- sweep(sweep(as.matrix(x), 2L, fit$center, "-"), 2L, fit$scale, "/")
    if (method == "nnet") {
      target <- outer(as.character(y), fit$codes, `==`) * 1
      colnames(target) <- fit$codes
      fit$engine <- withCallingHandlers(withr::with_seed(fit$seed,
        nnet::nnet(z, target, size = params$size, decay = params$decay,
          maxit = params$maxit, rang = params$rang, softmax = TRUE,
          MaxNWts = params$max_weights, trace = FALSE)), warning = record_warning)
      fit$parameters$entropy <- NULL
      fit$parameters$softmax <- TRUE
      fit$parameters$convergence <- fit$engine$convergence
      if (fit$engine$convergence != 0L) fit$warnings <- c(fit$warnings,
        "nnet reached its iteration limit; convergence was not established.")
    } else {
      if (nrow(x) < params$k) sb_abort("KNN requires at least k training rows.",
        "stormbreaker_bad_multiclass_fit")
      fit$engine <- list(training = z, outcome = y, k = params$k)
    }
  }
  fit$parameters$seed <- fit$seed
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_predict_multiclass_native <- function(fit, x) {
  if (fit$method %in% sb_multiclass_multinomial_catalog()$method) return(sb_predict_multiclass_multinomial(fit, x))
  if (fit$method %in% sb_multiclass_specialist_catalog()$method) return(sb_predict_multiclass_specialist(fit, x))
  if (fit$method %in% sb_multiclass_boosting_catalog()$method) return(sb_predict_multiclass_boosting(fit, x))
  if (fit$method %in% sb_multiclass_discriminant_catalog()$method) return(
    sb_predict_multiclass_discriminant(fit, x))
  if (fit$method %in% sb_multiclass_tree_catalog()$method) return(sb_predict_multiclass_tree(fit, x))
  if (!is.na(fit$package) && !requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "to predict with this saved fit."), "stormbreaker_missing_dependency")
  if (startsWith(fit$method, "svm_")) return(
    attr(stats::predict(fit$engine, x, probability = TRUE), "probabilities")[, fit$codes, drop = FALSE])
  if (fit$method == "rda") {
    predict_method <- get("predict.rda", envir = asNamespace("klaR"), inherits = FALSE)
    return(predict_method(fit$engine, newdata = x)$posterior[, fit$codes, drop = FALSE])
  }
  if (fit$method == "kernel_naive_bayes") return(
    stats::predict(fit$engine, x)$posterior[, fit$codes, drop = FALSE])
  z <- sweep(sweep(as.matrix(x), 2L, fit$center, "-"), 2L, fit$scale, "/")
  if (fit$method == "nnet") return(stats::predict(fit$engine, z, type = "raw")[, fit$codes, drop = FALSE])
  t(vapply(seq_len(nrow(z)), function(i) {
    distance <- rowSums(sweep(fit$engine$training, 2L, z[i, ], "-")^2)
    index <- order(distance, seq_along(distance))[seq_len(fit$engine$k)]
    tabulate(as.integer(fit$engine$outcome[index]), nbins = length(fit$codes)) / fit$engine$k
  }, numeric(length(fit$codes))))
}

sb_multiclass_native_evidence <- function(fit) {
  if (fit$method %in% sb_multiclass_boosting_catalog()$method) return(sb_multiclass_boosting_evidence(fit))
  if (fit$method %in% sb_multiclass_tree_catalog()$method) return(sb_multiclass_tree_evidence(fit))
  list(features = data.frame(feature = fit$features, selected = NA), by_class = NULL)
}

sb_multiclass_native_plots <- function(fit, new_data = NULL) {
  if (fit$method == "loclda") return(sb_multiclass_loclda_plots(fit, new_data))
  if (fit$method %in% c("knn", "kernel_weighted_knn")) return(sb_multiclass_neighbor_plots(fit, new_data))
  if (fit$method == "pamr") return(sb_multiclass_pamr_plots(fit))
  if (fit$method == "gam") return(sb_multiclass_gam_plots(fit))
  if (fit$method == "pre") return(sb_multiclass_pre_plots(fit))
  if (fit$method == "gaussian_process") return(sb_multiclass_gp_plots(fit))
  if (fit$method %in% c("cforest", "bagging", "rlt", "aorsf", "wsrf")) return(sb_multiclass_forest_structure_plots(fit))
  if (fit$method %in% c("edda", "mclustda", "hdda", "hda")) return(sb_multiclass_specialist_plots(fit))
  if (fit$method %in% c("liblinear_l2_logistic", "liblinear_l1_logistic", "xgboost_linear")) return(sb_multiclass_linear_weight_plots(fit))
  if (fit$method %in% c("fda", "mda", "pda", "slda", "sda")) return(sb_multiclass_projection_plots(fit))
  if (fit$method %in% c("kernel_naive_bayes", "rda")) return(sb_multiclass_distribution_plots(fit))
  if (fit$method %in% c("logistic", "stepglm", "brglm2")) return(sb_multiclass_logistic_plots(fit))
  if (fit$method == "nnet") return(sb_multiclass_neural_plots(fit))
  if (fit$method %in% sb_multiclass_boosting_catalog()$method) return(sb_multiclass_boosting_plots(fit))
  if (fit$method %in% sb_multiclass_tree_catalog()$method) return(sb_multiclass_tree_plots(fit))
  if (startsWith(fit$method, "svm_")) {
    codes <- fit$engine$levels[fit$engine$labels]
    tab <- data.frame(class = fit$classes[match(codes, fit$codes)], count = fit$engine$nSV)
    tab$class <- factor(tab$class, levels = fit$classes)
    class <- count <- NULL
    plot <- ggplot2::ggplot(tab, ggplot2::aes(class, count)) +
      ggplot2::geom_col(fill = "#2878A0") + sb_viz_theme() +
      ggplot2::labs(title = "SVM: Support vectors", x = NULL, y = "Support vectors")
    decisions <- fit$engine$decision.values
    if (is.null(decisions)) return(list(plots = list(support_vectors = plot), tables = list(support_vectors = tab)))
    pairs <- colnames(decisions)
    for (i in order(nchar(fit$codes), decreasing = TRUE))
      pairs <- gsub(fit$codes[i], fit$classes[i], pairs, fixed = TRUE)
    scores <- data.frame(sample = rep(seq_len(nrow(decisions)), ncol(decisions)),
      pair = rep(pairs, each = nrow(decisions)), decision = as.vector(decisions))
    decision <- NULL
    distribution <- ggplot2::ggplot(scores, ggplot2::aes(decision)) +
      ggplot2::stat_ecdf(geom = "step", colour = "#2878A0") +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
      ggplot2::facet_wrap(~pair) + sb_viz_theme() +
      ggplot2::labs(title = "SVM: Pairwise decision scores", x = "Decision score", y = "Cumulative fraction",
        caption = "Original training-row one-vs-one decision scores, not probabilities or held-out performance. Every training row is retained for every native class pair; zero is the pairwise decision boundary.")
    return(list(plots = list(support_vectors = plot, decision_distribution = distribution),
      tables = list(support_vectors = tab, decision_distribution = scores)))
  }
  list(plots = list(), tables = list())
}

# ---- multiclass-selection ----

sb_multiclass_selector_inputs <- function(x, y) {
  x <- sb_validate_predictor_frame(x, "x")
  if (!is.factor(y) || nlevels(y) < 3L || length(y) != nrow(x) ||
      anyNA(y) || any(table(y) == 0L)) sb_abort(
    "Supply a factor with at least three observed classes.", "stormbreaker_bad_multiclass_outcome")
  list(x = x, endpoint = list(outcome = y))
}

sb_analyze_multiclass_selector <- function(task, spec) {
  x <- task$data[task$predictors]
  y <- factor(task$data[[task$outcome_col]], levels = task$classes)
  method <- spec$method
  fit <- switch(method,
    boruta = sb_fit_boruta_selector(x, list(outcome = y),
      "boruta_multiclass_roughfix", capture_history = TRUE),
    varselrf = sb_fit_varselrf_selector(x, y, "varselrf_multiclass_oob_1se"),
    vsurf = sb_fit_vsurf_classification_selector(x, y, "vsurf_multiclass_prediction"))
  fit$classes <- task$classes
  types <- switch(method, boruta = c("boruta_summary", "boruta_hits", "boruta_history"),
    varselrf = "selection_history", vsurf = c("vsurf_importance", "vsurf_stages"))
  diagnostics <- plots <- notes <- list()
  status <- lapply(types, function(type) {
    error <- tryCatch({
      diagnostic <- sb_selector_diagnostics(fit, type)
      plot <- sb_plot_selector_diagnostic(diagnostic)
      notes[[type]] <<- list(labels = plot$labels, source = diagnostic$source)
      diagnostics[[type]] <<- diagnostic
      plots[[type]] <<- sb_method_plot(plot, type, spec$label)
      NULL
    }, error = identity)
    data.frame(view = type, status = if (is.null(error)) "available" else "error",
      reason = if (is.null(error)) "" else conditionMessage(error))
  })
  notes$interpretation <- paste("Full-input multiclass feature selection, not held-out evaluation.",
    "These are selectors, not probability-producing classifiers.",
    "Native forest OOB resampling is row-based, not patient-group-isolated.")
  list(status = "success", method = method, label = spec$label,
    feature_result = "selection", selection_rule = spec$selection_rule,
    features = fit$scores, class_features = NULL, evidence = fit$scores,
    parameters = fit$params, fit = fit, diagnostics = diagnostics, plots = plots,
    notes = notes, plot_status = do.call(rbind, status))
}

sb_fit_multiclass_abess <- function(fit, x, y) {
  if (!requireNamespace("abess", quietly = TRUE)) sb_abort(
    "Install abess for this method.", "stormbreaker_missing_dependency")
  support <- 0:min(10L, ncol(x))
  args <- list(x = as.matrix(x), y = y, family = "multinomial",
      tune.path = "sequence", tune.type = "gic", support.size = support,
      normalize = 2L, fit.intercept = TRUE, lambda = 0, c.max = 2L,
      splicing.type = 2L, max.splicing.iter = 20L, newton = "approx",
      max.newton.iter = 10L, early.stop = FALSE, warm.start = TRUE,
      important.search = ncol(x), num.threads = 1L, seed = fit$seed)
  run <- function(args) {
    warnings <- character()
    engine <- withCallingHandlers(withr::with_seed(args$seed,
      do.call(abess::abess, args)), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(engine = engine, warnings = warnings)
  }
  environment(run) <- baseenv()
  if (.Platform$OS.type == "windows") {
    # Native multinomial ABESS can terminate a mixed-engine Windows session.
    # A fresh process retains the exact algorithm while containing native faults.
    worker <- parallel::makePSOCKcluster(1L)
    on.exit(parallel::stopCluster(worker), add = TRUE)
    result <- parallel::clusterCall(worker, function(run, args, paths) {
      .libPaths(paths)
      run(args)
    }, run, args, .libPaths())[[1L]]
  } else result <- run(args)
  fit$engine <- result$engine
  fit$warnings <- c(fit$warnings, result$warnings)
  fit$engine$call <- NULL
  fit$parameters <- list(family = "multinomial", tune_type = "gic", support_size = support,
    selected_support_size = fit$engine$best.size, normalize = 2L, lambda = 0,
    newton = "approx", max_newton_iter = 10L,
    reference_class = utils::tail(fit$classes, 1L), seed = fit$seed)
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_multiclass_abess_evidence <- function(fit) {
  beta <- as.matrix(stats::coef(fit$engine, support.size = fit$engine$best.size)[[1]])
  # ABESS drops the last response column as reference. Its 0.4.11 coefficient
  # column labels instead drop the first; use the declared factor-level order.
  beta <- cbind(beta[fit$features, , drop = FALSE], 0)
  by_class <- do.call(rbind, lapply(seq_along(fit$classes), function(k)
    data.frame(feature = fit$features, class = fit$classes[k],
      coefficient = unname(beta[, k]), selected = unname(beta[, k] != 0),
      reference = k == length(fit$classes))))
  features <- data.frame(feature = fit$features, selected = rowSums(beta != 0) > 0,
    row.names = NULL)
  list(features = features, by_class = by_class)
}

sb_multiclass_abess_plots <- function(fit, evidence) {
  support_size <- gic <- coefficient <- feature <- NULL
  criterion <- data.frame(support_size = fit$engine$support.size, gic = fit$engine$tune.value)
  curve <- ggplot2::ggplot(criterion, ggplot2::aes(support_size, gic)) +
    ggplot2::geom_line(colour = "#2878A0") + ggplot2::geom_point(colour = "#2878A0") +
    ggplot2::geom_vline(xintercept = fit$engine$best.size, linetype = "dashed") +
    sb_viz_theme() + ggplot2::labs(title = "ABESS: Subset selection",
      x = "Features", y = "GIC")
  tab <- evidence$by_class
  tab$class <- factor(tab$class, levels = fit$classes)
  coefficients <- ggplot2::ggplot(tab, ggplot2::aes(coefficient, feature)) +
    ggplot2::geom_col(fill = "#2878A0") + ggplot2::facet_wrap(~class) +
    sb_viz_theme() + ggplot2::labs(title = "ABESS: Class coefficients",
      x = "Coefficient", y = NULL)
  list(plots = list(support_criterion = curve, class_coefficients = coefficients),
    tables = list(support_criterion = criterion, class_coefficients = tab))
}
