# ---- multiclass-classical ----

# ---- multiclass-classical ----

sb_multiclass_classical_catalog <- function() {
  data.frame(method = c("lda", "qda", "gaussian_naive_bayes", "decision_tree"),
    label = c("LDA", "QDA", "Gaussian naive Bayes", "Decision tree"),
    supported = TRUE, engine = c("MASS_lda", "MASS_qda", "e1071_naiveBayes", "rpart_class"),
    feature_result = c(rep("model_features", 3), "ranking"),
    selection_rule = "No standalone gene-selection rule")
}

sb_multiclass_classical_package <- function(method) {
  switch(method, lda = , qda = "MASS", gaussian_naive_bayes = "e1071", decision_tree = "rpart")
}

sb_fit_multiclass_classical <- function(fit, x, y) {
  package <- sb_multiclass_classical_package(fit$method)
  if (!requireNamespace(package, quietly = TRUE)) sb_abort(
    paste("Install", package, "for this method."), "stormbreaker_missing_dependency")
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  if (fit$method %in% c("lda", "qda")) {
    fit$engine <- withCallingHandlers(
      if (fit$method == "lda") MASS::lda(x, y, method = "moment") else
        MASS::qda(x, y, method = "moment"), warning = record_warning)
    fit$parameters <- list(method = "moment", prior = fit$engine$prior)
  } else if (fit$method == "gaussian_naive_bayes") {
    fit$engine <- withCallingHandlers(e1071::naiveBayes(x, y, laplace = 0),
      warning = record_warning)
    fit$parameters <- list(laplace = 0, prediction_eps = 0, prediction_threshold = 0.001)
  } else {
    fit$feature_codes <- paste0("feature", seq_along(fit$features))
    names(x) <- fit$feature_codes
    x$.outcome <- y
    fit$engine <- withCallingHandlers(rpart::rpart(.outcome ~ ., data = x,
      method = "class", parms = list(split = "gini"),
      control = rpart::rpart.control(cp = .01, minsplit = 20L, minbucket = 7L,
        maxdepth = 30L, maxsurrogate = 0L, usesurrogate = 0L, xval = 0L)),
      warning = record_warning)
    attr(fit$engine$terms, ".Environment") <- baseenv()
    fit$parameters <- list(cp = .01, minsplit = 20L, minbucket = 7L, maxdepth = 30L,
      split = "gini", maxsurrogate = 0L, usesurrogate = 0L, xval = 0L)
  }
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_predict_multiclass_classical <- function(fit, x) {
  package <- sb_multiclass_classical_package(fit$method)
  if (!requireNamespace(package, quietly = TRUE)) sb_abort(
    paste("Install", package, "to predict with this saved fit."), "stormbreaker_missing_dependency")
  if (fit$method %in% c("lda", "qda")) return(
    stats::predict(fit$engine, x)$posterior[, fit$codes, drop = FALSE])
  if (fit$method == "gaussian_naive_bayes") return(
    stats::predict(fit$engine, x, type = "raw", eps = fit$parameters$prediction_eps,
      threshold = fit$parameters$prediction_threshold)[, fit$codes, drop = FALSE])
  names(x) <- fit$feature_codes
  stats::predict(fit$engine, x, type = "prob")[, fit$codes, drop = FALSE]
}

sb_multiclass_classical_evidence <- function(fit) {
  features <- data.frame(feature = fit$features, selected = NA)
  if (fit$method == "decision_tree") {
    importance <- fit$engine$variable.importance
    features$importance <- 0
    present <- match(fit$feature_codes, names(importance))
    features$importance[!is.na(present)] <- unname(importance[present[!is.na(present)]])
    features$rank <- rank(-features$importance, ties.method = "min")
    features$used_for_split <- fit$feature_codes %in% fit$engine$frame$var
    return(list(features = features, by_class = NULL))
  }
  by_class <- do.call(rbind, lapply(seq_along(fit$classes), function(k) {
    tab <- data.frame(feature = fit$features, class = fit$classes[k], selected = NA)
    if (fit$method == "gaussian_naive_bayes") {
      tab$mean <- vapply(fit$features, function(gene) fit$engine$tables[[gene]][fit$codes[k], 1], numeric(1))
      tab$sd <- vapply(fit$features, function(gene) fit$engine$tables[[gene]][fit$codes[k], 2], numeric(1))
    } else tab$mean <- unname(fit$engine$means[fit$codes[k], fit$features])
    tab
  }))
  list(features = features, by_class = by_class)
}

sb_multiclass_classical_plots <- function(fit, evidence) {
  feature <- importance <- class <- mean <- sd <- NULL
  label <- sb_multiclass_classical_catalog()$label[
    match(fit$method, sb_multiclass_classical_catalog()$method)]
  if (fit$method == "decision_tree") {
    tab <- evidence$features
    tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$importance)])
    plot <- ggplot2::ggplot(tab, ggplot2::aes(importance, feature)) +
      ggplot2::geom_col(fill = "#2878A0") + sb_viz_theme() +
      ggplot2::labs(title = "Decision tree: Importance", x = "Native split importance", y = NULL)
    tree <- sb_multiclass_tree_structure(fit)
    return(list(plots = list(importance = plot, tree_structure = sb_plot_rpart_tree(tree)),
      tables = list(importance = tab, tree_structure = tree$table)))
  }
  tab <- evidence$by_class
  tab$class <- factor(tab$class, levels = fit$classes)
  plot <- ggplot2::ggplot(tab, ggplot2::aes(mean, feature)) + ggplot2::geom_point(colour = "#2878A0") +
    ggplot2::facet_wrap(~class) + sb_viz_theme() +
    ggplot2::labs(title = paste0(label, ": Class means"), x = "Mean", y = NULL)
  if (fit$method == "gaussian_naive_bayes") plot <- plot +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = mean - sd, xmax = mean + sd),
      orientation = "y", width = .2, colour = "#2878A0") +
    ggplot2::labs(title = "Gaussian naive Bayes: Class distributions", x = "Mean \u00b1 SD")
  list(plots = list(class_statistics = plot), tables = list(class_statistics = tab))
}

# ---- multiclass-distribution-plots ----

sb_multiclass_distribution_plots <- function(fit) {
  m <- fit$engine
  value <- density <- class <- feature <- NULL
  if (fit$method == "kernel_naive_bayes") {
    d <- do.call(rbind, lapply(fit$features, function(feature) {
      do.call(rbind, lapply(seq_along(fit$codes), function(k) {
        native <- m$tables[[feature]][[fit$codes[k]]]
        data.frame(feature = feature, class = fit$classes[k],
          value = native$x, density = native$y)
      }))
    }))
    p <- ggplot2::ggplot(d, ggplot2::aes(value, density, colour = class)) +
      ggplot2::geom_line() + ggplot2::facet_wrap(~feature, scales = "free", ncol = 2) +
      sb_viz_theme() + ggplot2::labs(title = "Class-conditional densities", x = "Feature value",
        y = "Density", colour = NULL,
        caption = "Saved kernel density grids used by the naive Bayes model; no new smoothing or density estimation is performed.")
    return(list(plots = list(class_densities = p), tables = list(class_densities = d)))
  }
  d <- data.frame(feature = rep(rownames(m$means), ncol(m$means)),
    class = rep(fit$classes[match(colnames(m$means), fit$codes)], each = nrow(m$means)),
    value = as.vector(m$means))
  p <- ggplot2::ggplot(d, ggplot2::aes(value, feature)) + ggplot2::geom_point(colour = "#2878A0") +
    ggplot2::facet_wrap(~class) + sb_viz_theme() +
    ggplot2::labs(title = "Class means", x = "Mean", y = NULL,
      caption = "Original fitted class means in input units, not variable importance or uncertainty intervals.")
  correlations <- do.call(rbind, lapply(seq_along(m$classes), function(k) {
    covariance <- (1 - m$regularization[2]) * m$covariances[, , k] + m$regularization[2] * m$covpooled
    covariance <- (1 - m$regularization[1]) * covariance +
      m$regularization[1] * mean(diag(covariance)) * diag(nrow(covariance))
    correlation <- stats::cov2cor(covariance)
    data.frame(feature = rep(rownames(correlation), ncol(correlation)),
      other = rep(colnames(correlation), each = nrow(correlation)),
      class = fit$classes[match(m$classes[k], fit$codes)], correlation = as.vector(correlation))
  }))
  other <- correlation <- NULL
  heat <- ggplot2::ggplot(correlations, ggplot2::aes(other, feature, fill = correlation)) +
    ggplot2::geom_tile() + ggplot2::facet_wrap(~class) +
    ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#C77335", limits = c(-1, 1)) +
    sb_viz_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(title = "Regularized class correlations", x = NULL, y = NULL, fill = "Correlation",
      caption = "Within-class correlations from the covariance matrices used for prediction, after both native RDA regularizations; not pooled empirical correlations.")
  list(plots = list(class_means = p, class_correlations = heat),
    tables = list(class_means = d, class_correlations = correlations))
}

# ---- multiclass-gp-plots ----

sb_multiclass_gp_plots <- function(fit) {
  m <- fit$engine
  pairs <- utils::combn(seq_along(m@lev), 2L)
  classes <- fit$classes[match(m@lev, fit$codes)]
  d <- do.call(rbind, lapply(seq_along(m@alpha), function(k) {
    data.frame(pair = paste(classes[pairs[, k]], collapse = " vs "),
      basis = seq_along(m@alpha[[k]]), reference_row = m@alphaindex[[k]],
      weight = as.numeric(m@alpha[[k]]))
  }))
  basis <- weight <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(basis, weight)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey65") +
    ggplot2::geom_point(colour = "#2878A0", size = 1.2) +
    ggplot2::facet_wrap(~pair) + sb_viz_theme() +
    ggplot2::labs(title = "Pairwise kernel weights", x = "Native basis", y = "Kernel weight",
      caption = "All original kernlab Gaussian-process one-vs-one kernel coefficients. Positive latent scores favor the first named class. These weights belong to reference observations, not genes; they are not feature importance or SHAP. No posterior uncertainty band is invented, and no observations or fitted coefficients are removed.")
  list(plots = list(kernel_weights = p), tables = list(kernel_weights = d))
}

# ---- multiclass-neighbor-plots ----

sb_multiclass_neighbor_plots <- function(fit, new_data) {
  if (is.null(new_data)) stop("Supply new_data for prediction-neighborhood plots.", call. = FALSE)
  x <- sb_validate_predictor_frame(new_data[fit$features], "new_data")
  m <- fit$engine
  if (fit$method == "knn") {
    z <- sweep(sweep(as.matrix(x), 2, fit$center, "-"), 2, fit$scale, "/")
    rows <- lapply(seq_len(nrow(z)), function(i) {
      distance <- rowSums(sweep(m$training, 2, z[i, ], "-")^2)
      index <- order(distance, seq_along(distance))[seq_len(m$k)]
      data.frame(query = i, rank = seq_along(index), reference_row = index,
        distance = sqrt(distance[index]), class = as.character(m$outcome[index]),
        weight = 1 / m$k)
    })
    tab <- do.call(rbind, rows)
  } else {
    names(x) <- fit$feature_codes
    p <- fit$parameters
    native <- kknn::kknn(stats::as.formula(".outcome ~ .", env = baseenv()),
      train = m$training, test = x, k = p$k, distance = p$distance,
      kernel = p$kernel, scale = p$scale)
    tab <- data.frame(query = rep(seq_len(nrow(x)), each = ncol(native$C)),
      rank = rep(seq_len(ncol(native$C)), nrow(x)), reference_row = as.integer(t(native$C)),
      distance = as.numeric(t(native$D)), class = as.character(t(native$CL)),
      weight = as.numeric(t(native$W / rowSums(native$W))))
  }
  tab$class <- fit$classes[match(tab$class, fit$codes)]
  rank <- distance <- class <- weight <- NULL
  p <- ggplot2::ggplot(tab, ggplot2::aes(rank, distance, colour = class, size = weight)) +
    ggplot2::geom_point(alpha = .8) + ggplot2::facet_wrap(~query) +
    ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
    sb_viz_theme() + ggplot2::labs(title = "Prediction neighborhoods", x = "Neighbor rank",
      y = "Distance", colour = "Class", size = "Voting weight",
      caption = "Original neighbor selection, distances and normalized class-voting weights for the supplied queries. Training queries retain self-neighbors; these are explanations, not held-out or leave-one-out performance. All queries and selected neighbors are retained across pages. Reference-row identities and classes may be sensitive.")
  list(plots = list(neighborhoods = p), tables = list(neighborhoods = tab))
}

# ---- multiclass-neural-plots ----

sb_multiclass_neural_plots <- function(fit) {
  m <- fit$engine
  labels <- c("Bias", fit$features, paste("Hidden", seq_len(m$n[2L])), fit$classes)
  destination <- rep(seq_len(length(m$nconn) - 1L) - 1L, diff(m$nconn))
  tab <- data.frame(source = labels[m$conn + 1L],
    destination = labels[destination + 1L], weight = m$wts,
    layer = ifelse(destination <= sum(m$n[1:2]), "Input to hidden", "Hidden to output"))
  tab$source <- factor(tab$source, levels = unique(labels))
  tab$destination <- factor(tab$destination, levels = unique(labels))
  source <- weight <- NULL
  plot <- ggplot2::ggplot(tab, ggplot2::aes(destination, source, fill = weight)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~layer, scales = "free") +
    sb_neural_weight_scale() +
    sb_viz_theme() + ggplot2::labs(title = "Neural network weights", x = NULL,
      y = NULL, fill = "Weight", caption = paste(
        "Original fitted connections, including biases. Inputs were standardized;",
        "outputs use softmax. Connection weights are not gene importance or causal effects.",
        "Weight colors use symmetric pseudo-log scaling (sigma = 1), with original weights on the legend; no coefficient clipping."))
  list(plots = list(connection_weights = plot), tables = list(connection_weights = tab))
}

# ---- multiclass-discriminants ----

# ---- multiclass-discriminants ----

sb_multiclass_discriminant_catalog <- function() {
  data.frame(method = c("fda", "mda", "pda", "slda", "sda", "loclda", "pamr",
    "gaussian_process", "kernel_weighted_knn", "liblinear_l2_logistic", "liblinear_l1_logistic"),
    label = c("Flexible discriminant analysis", "Mixture discriminant analysis",
      "Penalized discriminant analysis", "Stabilized LDA", "Shrinkage discriminant analysis",
      "Local LDA", "Nearest shrunken centroids", "Gaussian process", "Kernel-weighted KNN",
      "Liblinear L2 logistic", "Liblinear L1 logistic"), supported = TRUE,
    engine = c(rep("mda", 3), "ipred", "sda", "klaR", "pamr", "kernlab", "kknn", rep("LiblineaR", 2)),
    feature_result = "model_features", selection_rule = "No standalone gene-selection rule")
}

sb_fit_multiclass_discriminant <- function(fit, x, y) {
  method <- fit$method
  catalog <- sb_multiclass_discriminant_catalog()
  fit$package <- catalog$engine[match(method, catalog$method)]
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "for this method."), "stormbreaker_missing_dependency")
  params <- sb_classical_binary_params(paste0(method, "_binary"))
  params$seed <- fit$seed
  fit$feature_codes <- paste0("feature", seq_along(fit$features))
  names(x) <- fit$feature_codes
  matrix_x <- as.matrix(x)
  d <- x
  d$.outcome <- y
  formula <- stats::as.formula(".outcome ~ .", env = baseenv())
  priors <- as.numeric(table(y)) / length(y)
  if (method %in% c("fda", "pda")) params$dimension <- min(nlevels(y) - 1L, ncol(x))
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  fit$engine <- withCallingHandlers(withr::with_seed(fit$seed, {
    switch(method,
      fda = mda::fda(formula, data = d, dimension = params$dimension, method = mda::polyreg,
        degree = params$degree, monomial = params$monomial, keep.fitted = FALSE),
      pda = mda::fda(formula, data = d, dimension = params$dimension, method = mda::gen.ridge,
        lambda = params$lambda, keep.fitted = FALSE),
      mda = mda::mda(formula, data = d, subclasses = params$subclasses_per_class,
        dimension = params$dimension, iter = params$iterations, method = mda::polyreg,
        degree = params$degree, keep.fitted = FALSE, trace = FALSE),
      slda = ipred::slda(formula, data = d, q = params$q, prior = priors),
      sda = sda::sda(matrix_x, y, diagonal = FALSE, verbose = FALSE),
      loclda = klaR::loclda(matrix_x, y,
        k = sb_loclda_neighbors(params, ncol(x), nrow(x), nlevels(y)),
        weighted.apriori = params$weighted_apriori),
      pamr = {
        native <- NULL
        invisible(utils::capture.output(native <- pamr::pamr.train(
          list(x = t(matrix_x), y = y, geneid = fit$feature_codes, genenames = fit$feature_codes),
          threshold = params$threshold, scale.sd = params$scale_sd,
          offset.percent = params$offset_percent, prior = priors,
          remove.zeros = params$remove_zero_thresholds, sign.contrast = "both")))
        native
      },
      gaussian_process = kernlab::gausspr(matrix_x, y, type = "classification",
        kernel = params$kernel, kpar = "automatic", scaled = TRUE,
        tol = params$tolerance, cross = 0L, fit = FALSE),
      kernel_weighted_knn = list(training = d),
      {
        fit$center <- colMeans(matrix_x)
        fit$scale <- apply(matrix_x, 2L, stats::sd)
        fit$scale[fit$scale == 0] <- 1
        z <- sweep(sweep(matrix_x, 2L, fit$center, "-"), 2L, fit$scale, "/")
        params$score_mapping <- "native_multiclass_OVR_probabilities"
        LiblineaR::LiblineaR(z, y, type = params$solver_type, cost = params$cost,
          epsilon = params$tolerance, bias = params$bias, verbose = FALSE)
      })
  }), warning = record_warning)
  fit$parameters <- params
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_predict_multiclass_discriminant <- function(fit, x) {
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "to predict with this saved fit."), "stormbreaker_missing_dependency")
  names(x) <- fit$feature_codes
  method <- fit$method
  p <- if (method %in% c("fda", "mda", "pda")) {
    # mda's quadratic polybasis drops dimensions for one query row.
    query <- if (method == "fda" && nrow(x) == 1L) x[c(1L, 1L), , drop = FALSE] else x
    native <- stats::predict(fit$engine, newdata = query, type = "posterior")
    native[seq_len(nrow(x)), , drop = FALSE]
  } else if (method == "sda") {
    stats::predict(fit$engine, as.matrix(x), verbose = FALSE)$posterior
  } else if (method == "slda") {
    stats::predict(fit$engine, x)$posterior
  } else if (method == "loclda") {
    native <- stats::predict(fit$engine, as.matrix(x))$posterior
    colnames(native) <- fit$engine$lev
    native
  } else if (method == "pamr") {
    pamr::pamr.predict(fit$engine, t(as.matrix(x)), threshold = fit$parameters$threshold,
      type = "posterior")
  } else if (method == "gaussian_process") {
    kernlab::predict(fit$engine, as.matrix(x), type = "probabilities")
  } else if (method == "kernel_weighted_knn") {
    kknn::kknn(stats::as.formula(".outcome ~ .", env = baseenv()),
      train = fit$engine$training, test = x, k = fit$parameters$k,
      distance = fit$parameters$distance, kernel = fit$parameters$kernel,
      scale = fit$parameters$scale)$prob
  } else {
    z <- sweep(sweep(as.matrix(x), 2L, fit$center, "-"), 2L, fit$scale, "/")
    stats::predict(fit$engine, z, proba = TRUE)$probabilities
  }
  p[, fit$codes, drop = FALSE]
}

# ---- multiclass-loclda-plots ----

sb_multiclass_loclda_plots <- function(fit, new_data) {
  if (is.null(new_data)) stop("Supply new_data for local-LDA neighborhoods.", call. = FALSE)
  x <- as.matrix(sb_validate_predictor_frame(new_data[fit$features], "new_data"))
  tab <- sb_loclda_neighbor_table(fit$engine, x)
  tab$class <- fit$classes[match(tab$class, fit$codes)]
  distance <- weight <- class <- used <- NULL
  p <- ggplot2::ggplot(tab, ggplot2::aes(distance, weight, colour = class, shape = used)) +
    ggplot2::geom_point(size = 1.8, alpha = .8) + ggplot2::facet_wrap(~query) +
    ggplot2::scale_shape_manual(values = c("FALSE" = 4, "TRUE" = 16),
      labels = c("FALSE" = "Excluded", "TRUE" = "Used")) + sb_viz_theme() +
    ggplot2::labs(title = "Local-LDA neighborhoods", x = "Distance", y = "Local weight",
      colour = "Class", shape = NULL,
      caption = "Original distance-radius neighborhoods, including distance ties. Classes represented by at most one neighbor are excluded by the native local-LDA rule, but remain visible. All supplied queries are retained across pages. Self-neighbors remain for training queries. Local fitting weights are not probabilities, feature importance or validation performance; reference rows and classes may be sensitive.")
  list(plots = list(neighborhoods = p), tables = list(neighborhoods = tab))
}

# ---- multiclass-linear ----

# ---- multiclass-linear-weight-plots ----

sb_multiclass_linear_weight_plots <- function(fit) {
  if (fit$method == "xgboost_linear") {
    booster <- sb_load_xgboost_state_model(list(model_raw = fit$engine,
      params = fit$parameters, features = fit$features), "xgboost_binary_linear")$booster
    model <- jsonlite::fromJSON(rawToChar(xgboost::xgb.save.raw(booster, raw_format = "json")))
    weights <- model$learner$gradient_booster$model$weights
    coefficients <- matrix(weights, nrow = length(fit$classes))
    classes <- fit$classes
    features <- c(fit$features, "Bias")
    note <- "Original full-precision gblinear class-margin weights and learned biases from the saved model. The global base score is separate. These are not tree gain or reference-class logistic coefficients."
  } else {
    coefficients <- fit$engine$W
    features <- colnames(coefficients)
    at <- match(features, fit$feature_codes)
    features[!is.na(at)] <- fit$features[at[!is.na(at)]]
    classes <- fit$classes[match(as.character(fit$engine$ClassNames), fit$codes)]
    note <- "Original LIBLINEAR one-versus-rest weights on standardized inputs, including the bias weight. Not multinomial reference-class coefficients; class probabilities are subsequently normalized by the native predictor."
  }
  d <- data.frame(feature = rep(features, each = nrow(coefficients)),
    class = rep(classes, ncol(coefficients)), weight = as.vector(coefficients))
  d$feature <- factor(d$feature, levels = rev(features))
  feature <- weight <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(weight, feature)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
    ggplot2::geom_point(colour = "#2878A0") + ggplot2::facet_wrap(~class) +
    sb_viz_theme() + ggplot2::labs(title = "Class-specific linear weights", x = "Weight",
      y = NULL, caption = note)
  list(plots = list(linear_weights = p), tables = list(linear_weights = d))
}

# ---- multiclass-logistic-plots ----

sb_multiclass_logistic_plots <- function(fit) {
  coefficients <- stats::coef(fit$engine)
  terms <- colnames(coefficients)
  at <- match(terms, fit$feature_codes)
  terms[!is.na(at)] <- fit$features[at[!is.na(at)]]
  classes <- fit$classes[match(rownames(coefficients), fit$codes)]
  d <- data.frame(feature = rep(terms, each = nrow(coefficients)),
    contrast = rep(paste(classes, "vs", fit$classes[1L]), ncol(coefficients)),
    coefficient = as.vector(coefficients))
  d$feature <- factor(d$feature, levels = rev(terms))
  coefficient <- feature <- NULL
  plot <- ggplot2::ggplot(d, ggplot2::aes(coefficient, feature)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
    ggplot2::geom_point(colour = "#2878A0") +
    ggplot2::facet_wrap(~contrast) + sb_viz_theme() +
    ggplot2::labs(title = "Multinomial coefficients", x = "Log-odds coefficient", y = NULL,
      caption = "Original fitted coefficients, including the intercept, relative to the reference class. No refitting or coefficient thresholding is performed.")
  result <- list(plots = list(coefficients = plot), tables = list(coefficients = d))
  if (fit$method == "stepglm") {
    path <- fit$engine$anova
    path$step <- seq_len(nrow(path)) - 1L
    step <- AIC <- NULL
    result$plots$selection_path <- ggplot2::ggplot(path, ggplot2::aes(step, AIC)) +
      ggplot2::geom_line(colour = "#2878A0") + ggplot2::geom_point(colour = "#2878A0") +
      sb_viz_theme() + ggplot2::labs(title = "Stepwise selection", x = "Step", y = "AIC",
        caption = "Actual accepted stepAIC steps, including the initial model; not a cross-validation performance curve.")
    result$tables$selection_path <- path
  }
  result
}

# ---- multiclass-multinomial ----

sb_multiclass_multinomial_catalog <- function() {
  data.frame(method = c("logistic", "stepglm", "brglm2", "gam", "pre", "slope"),
    label = c("Multinomial logistic", "Stepwise multinomial logistic", "Bias-reduced multinomial",
      "Multinomial GAM", "Prediction rule ensemble", "Multinomial SLOPE"),
    supported = c(rep(TRUE, 5L), FALSE),
    engine = c("nnet", "nnet", "brglm2", "mgcv", "pre", "SLOPE"),
    feature_result = "model_features", selection_rule = c(rep(
      "Predictive model; no standalone gene-selection interface", 5L),
      "Experimental adapter; multiclass validation incomplete; not available"))
}

sb_fit_multiclass_multinomial <- function(fit, x, y) {
  method <- fit$method
  catalog <- sb_multiclass_multinomial_catalog()
  fit$package <- catalog$engine[match(method, catalog$method)]
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "for this method."), "stormbreaker_missing_dependency")
  fit$feature_codes <- paste0("feature", seq_along(fit$features))
  names(x) <- fit$feature_codes
  d <- x
  d$.outcome <- y
  formula <- stats::as.formula(".outcome ~ .", env = environment())
  params <- list()
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  fit$engine <- withCallingHandlers(withr::with_seed(fit$seed, {
    if (method %in% c("logistic", "stepglm")) {
      params <- list(maxit = 200L, decay = 0, reference_class = fit$classes[1L])
      native <- nnet::multinom(formula, data = d, trace = FALSE, maxit = params$maxit,
        MaxNWts = 1000000L, decay = 0, model = FALSE)
      if (method == "stepglm") {
        params$selection <- "both_direction_AIC"
        native <- MASS::stepAIC(native, direction = "both", trace = FALSE)
      }
      native
    } else if (method == "brglm2") {
      params <- list(type = "AS_mean", reference_class = fit$classes[1L])
      brglm2::brmultinom(formula, data = d, ref = 1L, type = "AS_mean")
    } else if (method == "gam") {
      params <- list(family = "multinom", basis = "cs", basis_dimension = 5L, method = "REML")
      d$.outcome <- as.integer(y) - 1L
      rhs <- paste0("s(", names(x), ",bs='cs',k=5)", collapse = " + ")
      formulas <- lapply(seq_len(nlevels(y) - 1L), function(k)
        stats::as.formula(paste(if (k == 1L) ".outcome" else "", "~", rhs),
          env = asNamespace("mgcv")))
      mgcv::gam(formulas, data = d, family = mgcv::multinom(K = nlevels(y) - 1L),
        method = "REML", control = mgcv::gam.control(nthreads = 1L))
    } else if (method == "pre") {
      params <- sb_rule_aft_params("pre_binary_rule_ensemble")
      params$family <- "multinomial"
      params$internal_folds <- length(unique(fit$fold_id))
      pre::pre(formula, data = d, family = "multinomial", type = params$ensemble_type,
        sampfrac = params$sample_fraction, maxdepth = params$maximum_depth,
        learnrate = params$learning_rate, mtry = max(1L, ceiling(sqrt(ncol(x)))),
        ntrees = params$num_trees, normalize = TRUE, standardize = FALSE,
        nfolds = params$internal_folds, foldid = fit$fold_id, tree.unbiased = TRUE,
        removecomplements = TRUE, removeduplicates = TRUE, verbose = FALSE)
    } else {
      params <- sb_slope_randomglm_binary_params("slope_binary")
      params$family <- "multinomial"
      # The native gamma=0 relaxed multinomial refit can silently produce
      # NaN intercepts. Use ordinary penalized SLOPE, not relaxed refitting.
      params$gamma <- 1
      params$internal_resampling <- "task_stratified_or_grouped_folds"
      params$inner_folds <- length(unique(fit$fold_id))
      params$probability_mapping <- "multinomial_inverse_link_last_class_reference"
      native <- sb_multiclass_slope_cv(as.matrix(x), y, fit$fold_id, params)
      params$optimum <- native$optima
      fit$tuning <- native$summary
      if (any(!is.finite(unlist(native$model$intercepts))) ||
          any(!is.finite(as.matrix(native$model$coefficients[[1L]])))) sb_abort(
        "SLOPE returned non-finite multinomial coefficients.", "stormbreaker_bad_multiclass_fit")
      native$model
    }
  }), warning = record_warning)
  fit$parameters <- params
  fit$parameters$seed <- fit$seed
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_predict_multiclass_multinomial <- function(fit, x) {
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "to predict with this saved fit."), "stormbreaker_missing_dependency")
  names(x) <- fit$feature_codes
  if (fit$method == "slope") {
    # SLOPE stores K-1 logits against the last class. Its native response
    # method omits that class and uses sum(1 + exp(eta)) as denominator.
    eta <- stats::predict(fit$engine, as.matrix(x), type = "link", simplify = FALSE)
    if (!identical(dim(eta), c(nrow(x), length(fit$codes) - 1L, 1L))) sb_abort(
      "Expected one selected SLOPE model with K-1 class logits.", "stormbreaker_bad_multiclass_fit")
    eta <- cbind(matrix(eta, nrow(x), length(fit$codes) - 1L), 0)
    p <- exp(eta - apply(eta, 1L, max))
    p <- p / rowSums(p)
    colnames(p) <- fit$engine$class_names
    return(p[, fit$codes, drop = FALSE])
  }
  p <- switch(fit$method,
    logistic = , stepglm = stats::predict(fit$engine, x, type = "probs"),
    brglm2 = stats::predict(fit$engine, x, type = "probs"),
    gam = stats::predict(fit$engine, x, type = "response"),
    pre = stats::predict(fit$engine, x, type = "response", penalty.par.val = "lambda.1se"))
  if (fit$method == "gam") return(p)
  if (is.null(dim(p)) && nrow(x) == 1L) p <- matrix(p, nrow = 1L, dimnames = list(NULL, names(p)))
  p[, fit$codes, drop = FALSE]
}

sb_multiclass_slope_cv <- function(x, y, folds, params) {
  # The installed cvSLOPE multinomial path gives NaN fold scores even when
  # SLOPE fits on those training rows are finite. Evaluate native fits here.
  train <- function(rows, q, alpha = "path") SLOPE::SLOPE(
    x[rows, , drop = FALSE], y[rows], family = "multinomial", q = q,
    alpha = alpha, gamma = params$gamma, center = params$center, scale = params$scale,
    path_length = params$path_length, threads = 1L,
    tol_dev_change = 0, tol_dev_ratio = 1)
  summaries <- lapply(params$q, function(q) {
    path <- train(seq_len(nrow(x)), q)
    alpha <- path$alpha
    loss <- matrix(NA_real_, nrow(x), length(alpha))
    for (fold in sort(unique(folds))) {
      test <- which(folds == fold)
      model <- train(which(folds != fold), q, alpha)
      eta <- stats::predict(model, x[test, , drop = FALSE], type = "link", simplify = FALSE)
      if (!identical(dim(eta), c(length(test), nlevels(y) - 1L, length(alpha))) ||
          any(!is.finite(eta))) sb_abort(
        "SLOPE returned an incomplete or non-finite inner-CV path.", "stormbreaker_bad_multiclass_fit")
      for (j in seq_along(alpha)) {
        logits <- cbind(matrix(eta[, , j], length(test), nlevels(y) - 1L), 0)
        maximum <- apply(logits, 1L, max)
        normalizer <- maximum + log(rowSums(exp(logits - maximum)))
        loss[test, j] <- 2 * (normalizer - logits[cbind(seq_along(test), as.integer(y[test]))])
      }
    }
    data.frame(q = q, gamma = params$gamma, alpha = alpha, measure = "deviance",
      mean = colMeans(loss))
  })
  summary <- do.call(rbind, summaries)
  optimum <- summary[which.min(summary$mean), , drop = FALSE]
  list(model = train(seq_len(nrow(x)), optimum$q, optimum$alpha),
    optima = optimum, summary = summary)
}
