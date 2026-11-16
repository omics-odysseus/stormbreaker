# ---- multiclass-omics ----

# ---- multiclass-omics ----

sb_multiclass_omics_catalog <- function() {
  data.frame(method = "lpc", label = "LPC", supported = TRUE, engine = "lpc",
    feature_result = "selection",
    selection_rule = "Top min(p, max(2, floor(n/log(n)))) absolute native LPC scores")
}

sb_analyze_multiclass_lpc <- function(task, spec) {
  if (!requireNamespace("lpc", quietly = TRUE)) sb_abort(
    "Install lpc for this selector.", "stormbreaker_missing_dependency")
  x <- task$data[task$predictors]
  # Native automatic threshold selection averages ranks 1:50 unconditionally.
  if (ncol(x) < 50L) sb_abort(
    "Native LPC automatic threshold selection requires at least 50 candidate genes.",
    "stormbreaker_lpc_too_few_predictors")
  y <- as.integer(factor(task$data[[task$outcome_col]], levels = task$classes))
  native <- NULL
  invisible(utils::capture.output(native <- withr::with_seed(1L,
    lpc::LPC(t(as.matrix(x)), y, type = "multiclass", soft.thresh = NULL))))
  scores <- as.numeric(native$lpcscores)
  if (length(scores) != ncol(x) || any(!is.finite(scores)) ||
      length(native$soft.thresh) != 1L || !is.finite(native$soft.thresh)) sb_abort(
    "LPC did not return finite gene scores and a selected threshold.",
    "stormbreaker_invalid_lpc_selector")
  size <- sb_high_dimensional_nlogn_size(nrow(x), ncol(x))
  ranking <- order(-abs(scores), seq_along(scores))
  selected <- utils::head(ranking, size)
  tab <- data.frame(feature = names(x), lpc_score = scores,
    base_score = as.numeric(native$tscores), rank = match(seq_along(x), ranking),
    selected = seq_along(x) %in% selected)
  params <- list(native_type = "multiclass", seed = 1L, screen_size = size,
    soft_threshold = native$soft.thresh, internal_resampling = "native_row_based_threshold_search")
  fit <- structure(list(id = "lpc_multiclass_nlogn", source_features = names(x),
    selected_features = names(x)[selected], scores = tab, params = params,
    classes = task$classes), class = "stormbreaker_selector_fit")
  feature <- lpc_score <- NULL
  shown <- tab[selected, ]
  shown$feature <- factor(shown$feature, levels = rev(shown$feature))
  plot <- ggplot2::ggplot(shown, ggplot2::aes(lpc_score, feature)) +
    ggplot2::geom_point(colour = "#2878A0") + ggplot2::theme_minimal() +
    ggplot2::labs(title = "LPC", x = "LPC score", y = NULL)
  list(status = "success", method = spec$method, label = spec$label,
    feature_result = "selection", selection_rule = spec$selection_rule,
    features = tab, class_features = NULL, evidence = tab, parameters = params,
    fit = fit, diagnostics = list(selection_scores = tab), plots = list(selection_scores = plot),
    notes = list(interpretation = paste("Full-input exploratory feature screening.",
      "Native threshold search is row-based, not patient-group-isolated.",
      "LPC scores do not give patient risk probabilities.")),
    plot_status = data.frame(view = "selection_scores", status = "available", reason = ""))
}

# ---- multiclass-pamr-plots ----

sb_multiclass_pamr_plots <- function(fit) {
  m <- fit$engine
  centers <- pamr::pamr.predict(m, newx = NULL,
    threshold = fit$parameters$threshold, type = "centroid")
  contrast <- (centers - m$centroid.overall) / m$sd
  d <- data.frame(feature = rep(fit$features[match(rownames(contrast), fit$feature_codes)], ncol(contrast)),
    class = rep(fit$classes[match(colnames(contrast), fit$codes)], each = nrow(contrast)),
    contrast = as.vector(contrast))
  d$feature <- factor(d$feature, levels = rev(fit$features))
  feature <- contrast <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(contrast, feature)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey65") +
    ggplot2::geom_col(fill = "#2878A0") + ggplot2::facet_wrap(~class) +
    sb_viz_theme() + ggplot2::labs(title = "Shrunken centroids",
      x = "Standardized centroid difference", y = NULL,
      caption = paste0("Native PAMR centroid contrasts at fitted threshold ", fit$parameters$threshold,
        "; centered on the overall centroid and divided by the native within-class scale, as in pamr.plotcen. All features, including zero contrasts, are retained. This is not a cross-validation curve; the configured fit uses one fixed threshold."))
  list(plots = list(shrunken_centroids = p), tables = list(shrunken_centroids = d))
}

# ---- multiclass-projection-plots ----

sb_multiclass_projection_plots <- function(fit) {
  m <- fit$engine
  result <- list(plots = list(), tables = list())
  add_matrix <- function(values, rows, columns, key, title, axis, note) {
    d <- data.frame(label = rep(rows, ncol(values)), panel = rep(columns, each = nrow(values)),
      value = as.vector(values))
    d$label <- factor(d$label, levels = rev(unique(rows)))
    value <- label <- NULL
    p <- ggplot2::ggplot(d, ggplot2::aes(value, label)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey75") +
      ggplot2::geom_point(colour = "#2878A0") + ggplot2::facet_wrap(~panel) +
      sb_viz_theme() + ggplot2::labs(title = title, x = axis, y = NULL, caption = note)
    result$plots[[key]] <<- p
    result$tables[[key]] <<- d
  }
  if (fit$method %in% c("fda", "mda", "pda")) {
    rows <- rownames(m$means)
    if (fit$method == "mda") {
      for (k in seq_along(m$assign.theta)) {
        at <- m$assign.theta[[k]]
        rows[at] <- paste(fit$classes[match(names(m$assign.theta)[k], fit$codes)],
          "subclass", seq_along(at))
      }
    } else rows <- fit$classes[match(rows, fit$codes)]
    add_matrix(m$means, rows, paste("Axis", seq_len(ncol(m$means))), "projected_centroids",
      "Discriminant centroids", "Discriminant coordinate",
      "Saved class or subclass centroids in the fitted discriminant space; these are not original gene expression means.")
    if (fit$method == "mda") {
      priors <- do.call(rbind, lapply(seq_along(m$sub.prior), function(k) {
        data.frame(class = fit$classes[match(names(m$sub.prior)[k], fit$codes)],
          subclass = names(m$sub.prior[[k]]), probability = as.numeric(m$sub.prior[[k]]))
      }))
      subclass <- probability <- NULL
      result$plots$subclass_weights <- ggplot2::ggplot(priors, ggplot2::aes(subclass, probability)) +
        ggplot2::geom_col(fill = "#2878A0") + ggplot2::facet_wrap(~class) + sb_viz_theme() +
        ggplot2::labs(title = "Subclass weights", x = NULL, y = "Within-class probability",
          caption = "Saved mixture probabilities within each outcome class, not overall class prevalence.")
      result$tables$subclass_weights <- priors
    }
  } else if (fit$method == "slda") {
    add_matrix(m$scores, fit$features, paste("Axis", seq_len(ncol(m$scores))), "projection_loadings",
      "Projection loadings", "Loading", "Original stabilized-LDA projection matrix; signs depend on axis orientation and are not causal effects.")
    add_matrix(m$mylda$means, fit$classes[match(rownames(m$mylda$means), fit$codes)],
      paste("Axis", seq_len(ncol(m$mylda$means))), "projected_centroids", "Projected class means",
      "Projected coordinate", "Fitted class means after the native stabilizing projection, not original gene expression means.")
  } else {
    add_matrix(t(m$beta), fit$features[match(colnames(m$beta), fit$feature_codes)],
      fit$classes[match(rownames(m$beta), fit$codes)], "discriminant_coefficients",
      "Discriminant coefficients", "Coefficient", "Original class-specific discriminant-score coefficients, not standalone gene importance.")
    d <- data.frame(parameter = names(m$regularization), intensity = unname(m$regularization))
    parameter <- intensity <- NULL
    result$plots$shrinkage <- ggplot2::ggplot(d, ggplot2::aes(intensity, parameter)) +
      ggplot2::geom_point(colour = "#2878A0") + sb_viz_theme() +
      ggplot2::labs(title = "Shrinkage intensities", x = "Intensity", y = NULL,
        caption = "Native shrinkage estimates for correlation, variance and class frequencies; not a tuning-search history.")
    result$tables$shrinkage <- d
  }
  result
}

# ---- multiclass-specialists ----

# ---- multiclass-gam-plots ----

sb_multiclass_gam_plots <- function(fit) {
  m <- fit$engine
  lpi <- attr(m$formula, "lpi")
  curves <- rugs <- list()
  for (j in seq_along(m$smooth)) {
    smooth <- m$smooth[[j]]
    term <- smooth$term
    index <- smooth$first.para:smooth$last.para
    lp <- which(vapply(lpi, function(z) all(index %in% z), logical(1)))
    values <- m$model[[term]]
    grid <- stats::setNames(data.frame(seq(min(values), max(values), length.out = 100)), term)
    X <- mgcv::PredictMat(smooth, grid)
    estimate <- drop(X %*% stats::coef(m)[index])
    se <- sqrt(rowSums((X %*% m$Vp[index, index, drop = FALSE]) * X))
    feature <- fit$features[match(term, fit$feature_codes)]
    curves[[j]] <- data.frame(feature = feature, class = fit$classes[lp + 1L],
      value = grid[[1]], estimate = estimate, se = se,
      lower = estimate - stats::qnorm(.975) * se, upper = estimate + stats::qnorm(.975) * se)
    rugs[[j]] <- data.frame(feature = feature, class = fit$classes[lp + 1L], value = values)
  }
  d <- do.call(rbind, curves)
  rug <- do.call(rbind, rugs)
  value <- estimate <- lower <- upper <- NULL
  plots <- list()
  for (cl in unique(d$class)) {
    plots[[paste0("smooths_", cl)]] <- ggplot2::ggplot(d[d$class == cl, ], ggplot2::aes(value, estimate)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#DDEBF3") +
      ggplot2::geom_hline(yintercept = 0, colour = "grey65", linetype = 2) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::geom_rug(data = rug[rug$class == cl, ], ggplot2::aes(x = value),
        inherit.aes = FALSE, sides = "b", alpha = .25) +
      ggplot2::facet_wrap(~feature, ncol = 1, scales = "free") + sb_viz_theme() +
      ggplot2::labs(title = paste(cl, "vs", fit$classes[1], "smooth effects"),
        x = "Feature value", y = "Partial log-odds",
        caption = paste0("Native centered GAM smooth terms for ", cl, " versus reference ", fit$classes[1],
          ". Pointwise approximate 95% Bayesian intervals use the fitted coefficient covariance conditional on smoothing parameters; not simultaneous bands or predicted probabilities. Rug marks retain every training value. No refitting, extrapolation or causal interpretation."))
  }
  list(plots = plots, tables = list(smooths = d, training_rug = rug))
}

# ---- multiclass-pre-plots ----

sb_multiclass_pre_plots <- function(fit) {
  m <- fit$engine
  native <- stats::coef(m, penalty.par.val = "lambda.1se")
  tab <- do.call(rbind, lapply(seq_along(fit$codes), function(k) {
    data.frame(feature = native$rule, description = native$description,
      class = fit$classes[k], coefficient = native[[paste0("coefficient.", fit$codes[k])]])
  }))
  for (i in order(nchar(fit$feature_codes), decreasing = TRUE))
    tab$description <- gsub(fit$feature_codes[i], fit$features[i], tab$description, fixed = TRUE)
  tab$term_type <- ifelse(tab$feature == "(Intercept)", "intercept",
    ifelse(tab$feature %in% m$rules$rule, "rule", "linear"))
  cv <- m$glmnet.fit
  tuning <- data.frame(lambda = cv$lambda, mean_value = cv$cvm, standard_error = cv$cvsd)
  plots <- list(coefficient_cv = sb_plot_pre_diagnostic(list(type = "rule_cv_curve",
    table = tuning, selected_lambda = cv$lambda.1se, selection_rule = "lambda.1se", measure = cv$name)))
  coefficient <- feature <- term_type <- NULL
  for (cl in fit$classes) {
    d <- tab[tab$class == cl & tab$coefficient != 0 & tab$term_type != "intercept", ]
    if (!nrow(d)) next
    d$feature <- factor(d$feature, levels = d$feature[order(d$coefficient)])
    plots[[paste0("nonzero_terms_", cl)]] <- ggplot2::ggplot(d, ggplot2::aes(coefficient, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey65") +
      ggplot2::geom_col(ggplot2::aes(fill = term_type)) +
      ggplot2::scale_fill_manual(values = c(rule = "#2878A0", linear = "#D77935"), name = "Term") +
      sb_viz_theme() + ggplot2::labs(title = paste(cl, "nonzero terms"),
        x = "Class-score coefficient", y = "Native term ID",
        caption = "Nonzero non-intercept terms at the saved lambda.1se. Complete class coefficients, intercepts, zero terms and rule definitions remain in the diagnostic table. Native multinomial scores use softmax, not one-reference-class log odds. Rule and linear terms have different units; coefficient magnitude is not feature importance. No refitting.")
  }
  list(plots = plots, tables = list(coefficients = tab, coefficient_cv = tuning))
}

# ---- multiclass-specialist-plots ----

sb_multiclass_specialist_plots <- function(fit) {
  m <- fit$engine
  feature <- value <- component <- class <- NULL
  if (fit$method %in% c("edda", "mclustda")) {
    means <- do.call(rbind, lapply(seq_along(m$models), function(k) {
      a <- m$models[[k]]$parameters$mean
      data.frame(feature = rep(fit$features[match(rownames(a), fit$feature_codes)], ncol(a)),
        component = rep(seq_len(ncol(a)), each = nrow(a)),
        class = fit$classes[match(names(m$models)[k], fit$codes)], value = as.vector(a))
    }))
    p <- ggplot2::ggplot(means, ggplot2::aes(value, feature, colour = factor(component))) +
      ggplot2::geom_point() + ggplot2::facet_wrap(~class) + sb_viz_theme() +
      ggplot2::labs(title = "Mixture component means", x = "Mean", y = NULL, colour = "Component",
        caption = "Original selected Gaussian mixture component means in input units. Component numbers are local to each class, not matched biological subtypes.")
    result <- list(plots = list(component_means = p), tables = list(component_means = means))
    bic <- do.call(rbind, lapply(seq_along(m$models), function(k) {
      a <- m$models[[k]]$BIC
      if (!is.matrix(a)) return(NULL)
      data.frame(components = rep(rownames(a), ncol(a)), covariance = rep(colnames(a), each = nrow(a)),
        class = fit$classes[match(names(m$models)[k], fit$codes)], BIC = as.vector(a))
    }))
    if (!is.null(bic)) {
      components <- covariance <- BIC <- NULL
      result$plots$bic_search <- ggplot2::ggplot(bic, ggplot2::aes(components, covariance,
        fill = ifelse(is.finite(BIC), BIC, NA_real_))) +
        ggplot2::geom_tile() + ggplot2::facet_wrap(~class) + sb_viz_theme() +
        ggplot2::scale_fill_gradient(low = "#EFF4F8", high = "#2878A0", na.value = "grey85") +
        ggplot2::labs(title = "Mixture model selection", x = "Components", y = "Covariance model", fill = "BIC",
          caption = "Saved candidate BIC values, with larger values preferred by mclust. Grey cells have no finite native estimate; the full table retains those unavailable candidates. This is model selection, not validation accuracy.")
      result$tables$bic_search <- bic
    }
    return(result)
  }
  if (fit$method == "hdda") {
    values <- as.matrix(m$ev)
    d <- data.frame(class = rep(fit$classes[match(m$kname, fit$codes)], ncol(values)),
      axis = rep(seq_len(ncol(values)), each = nrow(values)), value = as.vector(values))
    d$retained <- d$axis <= rep(as.numeric(m$d), ncol(values))
    axis <- retained <- NULL
    p <- ggplot2::ggplot(d, ggplot2::aes(axis, value)) +
      ggplot2::geom_line(colour = "grey60") + ggplot2::geom_point(ggplot2::aes(shape = retained), colour = "#2878A0") +
      ggplot2::facet_wrap(~class) + sb_viz_theme() +
      ggplot2::labs(title = "Class eigenvalue spectra", x = "Axis", y = "Eigenvalue", shape = "Retained",
        caption = "All saved class eigenvalues and the native retained dimensions, on the fitted scaled-input representation. No new PCA, threshold selection or truncation is performed.")
    return(list(plots = list(eigenvalue_spectrum = p), tables = list(eigenvalue_spectrum = d)))
  }
  a <- m$hda.loadings
  d <- data.frame(feature = rep(fit$features[match(rownames(a), fit$feature_codes)], ncol(a)),
    axis = rep(seq_len(ncol(a)), each = nrow(a)), value = as.vector(a))
  d$retained <- d$axis <= m$reduced.dimension
  axis <- retained <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(factor(axis), feature, fill = value)) +
    ggplot2::geom_tile() + ggplot2::facet_wrap(~retained, scales = "free_x", labeller = ggplot2::as_labeller(
      c(`TRUE` = "Used for classification", `FALSE` = "Remaining axes"))) +
    ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#C77335") + sb_viz_theme() +
    ggplot2::labs(title = "HDA projection loadings", x = "Axis", y = NULL, fill = "Loading",
      caption = "All original projection loadings, separating the configured classification dimensions from the remaining saved axes. Axis signs are arbitrary; these are not gene effect sizes.")
  list(plots = list(projection_loadings = p), tables = list(projection_loadings = d))
}

# ---- multiclass-specialists ----

sb_multiclass_specialist_catalog <- function() {
  data.frame(method = c("edda", "mclustda", "hdda", "hda", "rlt", "aorsf", "wsrf"),
    label = c("EDDA", "Mixture discriminant analysis (mclust)", "High-dimensional discriminant analysis",
      "Heteroscedastic discriminant analysis", "Reinforcement-learning forest",
      "Oblique random forest", "Weighted-subspace forest"), supported = TRUE,
    engine = c("mclust", "mclust", "HDclassif", "hda", "RLT", "aorsf", "wsrf"),
    feature_result = "model_features", selection_rule = "No standalone gene-selection rule")
}

sb_fit_multiclass_specialist <- function(fit, x, y) {
  method <- fit$method
  catalog <- sb_multiclass_specialist_catalog()
  fit$package <- catalog$engine[match(method, catalog$method)]
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "for this method."), "stormbreaker_missing_dependency")
  fit$feature_codes <- paste0("feature", seq_along(fit$features))
  names(x) <- fit$feature_codes
  matrix_x <- as.matrix(x)
  params <- list()
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  fit$engine <- withCallingHandlers(withr::with_seed(fit$seed, {
    if (method %in% c("edda", "mclustda")) {
      params <- sb_mclust_da_params(paste0(method, "_binary_bic"))
      params$active_covariance_models <- sb_mclust_da_active_models(ncol(x), params)
      if (method == "mclustda") params$initialization$seed <- fit$seed
      sb_mclust_da_fit_call(matrix_x, y, params, method)
    } else if (method == "hdda") {
      params <- sb_hdda_binary_params()
      params$probability_mapping <- "native_multiclass_posterior"
      HDclassif::hdda(matrix_x, cls = y, model = params$model, graph = FALSE,
        d_select = params$dimension_selection, threshold = params$threshold, show = FALSE,
        scaling = params$scaling, noise.ctrl = params$noise_control)
    } else if (method == "hda") {
      params <- list(reduced_dimension = min(nlevels(y) - 1L, ncol(x)),
        outer_iterations = 7L, inner_iterations = 10L, classification_rule = "naive_bayes")
      native <- NULL
      invisible(utils::capture.output(native <- hda::hda(matrix_x, grouping = y,
        newdim = params$reduced_dimension, crule = TRUE, reg.lamb = NULL, reg.gamm = NULL,
        initial.loadings = NULL, noutit = params$outer_iterations, ninit = params$inner_iterations,
        verbose = FALSE)))
      native
    } else if (method == "rlt") {
      params <- sb_rlt_learner_params("rlt_reinforcement_binary")
      RLT::RLT(x, y, model = "classification", print.summary = 0, use.cores = 1,
        ntrees = params$num_trees, nmin = params$native_nmin, reinforcement = TRUE,
        combsplit = 1L, embed.nmin = params$embedded_stop_size, importance = FALSE, track.obs = FALSE)
    } else if (method == "aorsf") {
      params <- sb_classical_binary_params("aorsf_binary")
      d <- x
      d$.outcome <- y
      aorsf::orsf(data = d, formula = stats::as.formula(".outcome ~ .", env = baseenv()),
        control = aorsf::orsf_control_classification(method = params$control_method,
          scale_x = params$control_scale_x, max_iter = params$control_max_iter),
        n_tree = params$num_trees, n_split = params$num_split, n_retry = params$num_retry,
        n_thread = 1L, mtry = max(1L, ceiling(sqrt(ncol(x)))),
        sample_with_replacement = params$sample_with_replacement, sample_fraction = params$sample_fraction,
        importance = "none", tree_seeds = sample.int(.Machine$integer.max, params$num_trees),
        attach_data = FALSE, na_action = "fail", verbose_progress = FALSE)
    } else {
      params <- sb_wsrf_binary_params()
      wsrf::wsrf(x, y, ntree = params$num_trees, mtry = min(ncol(x), floor(log2(ncol(x)) + 1)),
        weights = TRUE, parallel = FALSE, importance = FALSE, nodesize = params$leaf_node_size)
    }
  }), warning = record_warning)
  params$seed <- fit$seed
  fit$parameters <- params
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_predict_multiclass_specialist <- function(fit, x) {
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "to predict with this saved fit."), "stormbreaker_missing_dependency")
  names(x) <- fit$feature_codes
  p <- switch(fit$method,
    edda = , mclustda = stats::predict(fit$engine, newdata = as.matrix(x))$z,
    hdda = stats::predict(fit$engine, data = as.matrix(x))$posterior,
    hda = stats::predict(fit$engine, newdata = as.matrix(x), task = "c")$posteriors,
    rlt = {
      native <- stats::predict(fit$engine, x)$ProbPrediction
      colnames(native) <- fit$engine$ylevels
      native
    },
    aorsf = stats::predict(fit$engine, new_data = x, pred_type = "prob", n_thread = 1L),
    wsrf = stats::predict(fit$engine, newdata = x, type = "prob")$prob)
  p[, fit$codes, drop = FALSE]
}
