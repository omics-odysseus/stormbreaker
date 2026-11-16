#' Neighborhood and kernel diagnostics
#' @name sb_knn_neighborhood
#' @description Inspect neighborhoods, kernels, class distributions and mixture components of fitted models.
#' @rawRd \arguments{
#' \item{fit}{
#' \strong{For \code{sb_knn_neighborhood()}:} A successful ordinary or kernel-weighted binary KNN final fit.
#' 
#' \strong{For \code{sb_bagged_knn_diagnostics()}:} Successful bagged survival KNN final fit.
#' 
#' \strong{For \code{sb_loclda_neighborhood()}:} Successful local-LDA final fit.
#' 
#' \strong{For \code{sb_lvq_diagnostics()}:} A successful LVQ1, OLVQ1, LVQ2.1 or LVQ3 final fit.
#' 
#' \strong{For \code{sb_svm_diagnostics()}:} A successful final binary e1071 SVM, LIBLINEAR SVC, LS-SVM or DWD fit.
#' 
#' \strong{For \code{sb_survival_svm_diagnostics()}:} A final survival-SVM fit with capture_diagnostics = TRUE.
#' 
#' \strong{For \code{sb_gaussian_process_diagnostics()}:} Successful Gaussian-process binary final fit.
#' 
#' \strong{For \code{sb_naive_bayes_diagnostics()}:} A final Gaussian or kernel Naive Bayes model.
#' 
#' \strong{For \code{sb_mclust_da_diagnostics()}:} Successful EDDA or MclustDA final fit.
#' }
#' \item{newdata}{
#' \strong{For \code{sb_knn_neighborhood()}:} Explicit query rows with the original fitted features.
#' NULL reports unavailable rather than silently querying the training set.
#' 
#' \strong{For \code{sb_bagged_knn_diagnostics()}:} Explicit query rows required for bag_curves.
#' 
#' \strong{For \code{sb_loclda_neighborhood()}:} Explicit query rows; no implicit training queries.
#' 
#' \strong{For \code{sb_gaussian_process_diagnostics()}:} Explicit queries for query_contributions.
#' }
#' \item{type}{
#' \strong{For \code{sb_bagged_knn_diagnostics()}:} bootstrap or bag_curves.
#' 
#' \strong{For \code{sb_lvq_diagnostics()}:} prototypes or prototype_distances.
#' 
#' \strong{For \code{sb_svm_diagnostics()}:} support_counts, decision_distribution or coefficients (linear only).
#' 
#' \strong{For \code{sb_survival_svm_diagnostics()}:} score_distribution, basis_weights or coefficients (linear kernel).
#' 
#' \strong{For \code{sb_gaussian_process_diagnostics()}:} basis_weights or query_contributions.
#' 
#' \strong{For \code{sb_naive_bayes_diagnostics()}:} conditional_density or log_likelihood_ratio (Gaussian only).
#' 
#' \strong{For \code{sb_mclust_da_diagnostics()}:} centers, correlations, or mixture_weights.
#' }
#' \item{features}{
#' Predictor names to display; NULL retains all fitted predictors.
#' }
#' \item{grid_points}{
#' Gaussian grid size per feature. Kernel models retain their
#' original density grids and ignore this argument.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_knn_neighborhood()}}{
#' Explain original nearest-neighbor predictions
#' 
#' \strong{Returns:} A diagnostic with one row per query and neighbor, original training
#' row index, distance, class, native weight and normalized voting weight.
#' Prediction probabilities are returned separately. Ordinary KNN uses its
#' original standardization, Euclidean distance and training-row tie rule.
#' Weighted KNN uses the original kknn distance/kernel/scaling settings.
#' Training-row indices and classes are sensitive. Query rows that duplicate
#' training rows retain self-neighbors, exactly as prediction does; this is
#' not leave-one-out validation. No tuning or patient card is generated.
#' }
#' 
#' \subsection{\code{sb_bagged_knn_diagnostics()}}{
#' Inspect the original bagged survival KNN ensemble
#' 
#' \strong{Returns:} Original bootstrap multiplicities and per-bag coverage, or native
#' per-bag survival curves and their ensemble mean. Bag-to-bag variation is
#' not a confidence interval. No refitting or resampling occurs. Bootstrap
#' row indices are sensitive; repeated rows remain repeated in prediction.
#' }
#' 
#' \subsection{\code{sb_loclda_neighborhood()}}{
#' Inspect the original local-LDA prediction neighborhood
#' 
#' \strong{Returns:} All neighbors within the native kth-distance radius, including ties,
#' their native weights, classes and whether the native class-count rule uses
#' them. Native posterior probabilities are also returned. Reference indices
#' and classes may be sensitive. No tuning or separate classifier is fitted.
#' }
#' 
#' \subsection{\code{sb_lvq_diagnostics()}}{
#' Inspect the original learned LVQ prototypes
#' 
#' \strong{Returns:} An algorithm diagnostic using compact saved prediction state.
#' Prototype profiles retain the model's training-centered/scaled coordinates;
#' the table additionally contains their exact original-unit coordinates.
#' Distances are squared Euclidean distances in that same prediction metric,
#' without clustering, projection, refitting or restandardizing prototypes.
#' Native prototype order and encoded classes are retained. Class counts are
#' prototype counts, not sample counts or estimated prevalence. Learned
#' prototypes are not individual patients or class means. No training-history
#' or calibration claim is made: this learner maps the difference between
#' nearest-class squared distances through a fixed logistic function.
#' }
#' 
#' \subsection{\code{sb_svm_diagnostics()}}{
#' Inspect original binary SVM and DWD model structure
#' 
#' \strong{Returns:} An algorithm diagnostic. Counts and linear weights use portable state.
#' Decision distributions require capture_diagnostics = TRUE and evaluate the
#' original model on captured development predictors without fitting. Scores
#' are oriented toward encoded class 1, not calibrated probabilities or
#' geometric distances. Zero is the native classification boundary, not
#' necessarily the calibrated probability 0.5 boundary. Row-level scores remain
#' sensitive even though patient identifiers are omitted. Linear coefficients
#' and intercept are converted back to original predictor units; they are not
#' odds ratios or feature importance. No tuning or support-vector search is run.
#' }
#' 
#' \subsection{\code{sb_survival_svm_diagnostics()}}{
#' Inspect native survival-SVM structure without refitting
#' 
#' \strong{Returns:} An algorithm diagnostic. All twelve formulation/kernel combinations
#' retain their native score scale, before the package's Cox calibration.
#' Regression scores follow native predicted-time units; ranking/hybrid scores
#' are not probabilities or log hazards. Basis weights refer to kernel rows,
#' not features or patient importance. Linear coefficients include a separately
#' returned intercept and reproduce native scores in original input units.
#' Score and basis-row tables omit identifiers but remain sensitive training
#' data. No fitting, tuning, calibration or inverse probability calculation.
#' }
#' 
#' \subsection{\code{sb_gaussian_process_diagnostics()}}{
#' Inspect native Gaussian-process classifier kernel terms
#' 
#' \strong{Returns:} Original basis weights, or per-query kernel terms and native
#' probabilities. Terms are oriented toward class 1 and sum to its native
#' latent score; they are reference-row contributions, not feature SHAP.
#' No posterior interval is generated: this classifier does not save the
#' regression posterior covariance required by the native variance method.
#' Reference indices can be sensitive. No refitting or calibration occurs.
#' }
#' 
#' \subsection{\code{sb_naive_bayes_diagnostics()}}{
#' Inspect fitted Naive Bayes class-conditional distributions
#' 
#' \strong{Returns:} An algorithm diagnostic using compact saved state, without capture or
#' refitting. Gaussian curves use the original smoothed variances and a
#' model-implied grid spanning both class means plus/minus four fitted SDs.
#' Kernel curves retain every native density-grid point, with no new density
#' estimate or bandwidth choice. Density is not class probability or importance;
#' axes differ across features because their units differ. Gaussian likelihood
#' ratios are log f(x|class 1) minus log f(x|class 0), excluding prior log odds.
#' Prior log odds plus all feature contributions yield the original posterior
#' log odds under the fitted conditional-independence model.
#' The result's \code{source} and \code{interpretation} fields provide separate figure
#' legend text; plotting does not put these explanations inside the panels.
#' }
#' 
#' \subsection{\code{sb_mclust_da_diagnostics()}}{
#' Inspect original model-based discriminant components
#' 
#' \strong{Returns:} Original component means, within-component correlations, or
#' within-class mixing weights. Components are latent model structure, not
#' clinical subgroups. Covariances are original fitted parameters, not
#' empirical sample correlations. BIC paths absent from compact fits are not
#' reconstructed. No mixture refitting or component reselection occurs.
#' }
#' }
#' @rawRd \examples{
#' # sb_naive_bayes_diagnostics
#' d <- sb_example_binary(n = 60L, seed = 108L)
#' task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#' fit <- sb_fit_final(task, sb_pipeline(learner = "gaussian_naive_bayes"))
#' diagnostic <- sb_naive_bayes_diagnostics(fit)
#' diagnostic$interpretation
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   sb_plot_algorithm_diagnostic(diagnostic)
#' }
#' }
NULL

# ---- model-naive-bayes ----

# ---- algorithm-naive-bayes ----

#' @rdname sb_knn_neighborhood
#' @export
sb_naive_bayes_diagnostics <- function(fit, type = "conditional_density", features = NULL, grid_points = 101L) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  id <- fit$pipeline$learner
  if (!id %in% c("gaussian_naive_bayes", "kernel_naive_bayes_binary")) sb_abort(
    "This diagnostic requires Gaussian or kernel Naive Bayes.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("conditional_density", "log_likelihood_ratio")) sb_abort(
    "Choose conditional_density or log_likelihood_ratio.", "stormbreaker_unsupported_diagnostic")
  all_features <- fit$learner_state$features
  if (is.null(features)) features <- all_features
  if (!is.character(features) || !length(features) || anyNA(features) || anyDuplicated(features) ||
      !all(features %in% all_features)) sb_abort("Choose unique fitted predictor names.", "stormbreaker_invalid_diagnostic")
  prior_log_odds <- NULL
  if (id == "gaussian_naive_bayes") {
    points <- sb_check_positive_integer(grid_points, "grid_points", minimum = 2L)
    m <- fit$learner_state$payload
    prior_log_odds <- log(m$priors[["1"]]) - log(m$priors[["0"]])
    tables <- lapply(features, function(feature) {
      mu <- m$means[c("0", "1"), feature]
      sd <- sqrt(m$variances[c("0", "1"), feature])
      if (any(!is.finite(mu)) || any(!is.finite(sd)) || any(sd <= 0)) sb_abort(
        "Original Gaussian parameters are invalid.", "stormbreaker_invalid_diagnostic")
      x <- seq(min(mu - 4 * sd), max(mu + 4 * sd), length.out = points)
      if (type == "log_likelihood_ratio") return(data.frame(feature = feature, x = x,
        value = stats::dnorm(x, mu[2L], sd[2L], log = TRUE) - stats::dnorm(x, mu[1L], sd[1L], log = TRUE)))
      rbind(data.frame(feature = feature, class = "0", x = x, value = stats::dnorm(x, mu[1L], sd[1L])),
        data.frame(feature = feature, class = "1", x = x, value = stats::dnorm(x, mu[2L], sd[2L])))
    })
    native <- list(means = m$means, variances = m$variances, priors = m$priors, epsilon = m$epsilon)
    source <- "original Gaussian means and smoothed variances; model-implied +/-4 SD grid"
  } else {
    if (type == "log_likelihood_ratio") sb_abort(
      "Kernel likelihood contributions require native prediction-floor semantics; this view exposes only the original densities.",
      "stormbreaker_unavailable_diagnostic")
    e <- unserialize(fit$learner_state$engine_raw)
    if (!inherits(e, "NaiveBayes") || !isTRUE(e$usekernel) || !identical(e$levels, c("0", "1"))) sb_abort(
      "Original kernel Naive Bayes model is invalid.", "stormbreaker_invalid_diagnostic")
    tables <- lapply(features, function(feature) do.call(rbind, lapply(c("0", "1"), function(class) {
      density <- e$tables[[feature]][[class]]
      if (!inherits(density, "density") || length(density$x) != length(density$y) ||
          length(density$x) < 2L || any(diff(density$x) <= 0)) sb_abort(
        "Original kernel density grid is invalid.", "stormbreaker_invalid_diagnostic")
      data.frame(feature = feature, class = class, x = density$x, value = density$y)
    })))
    native <- e$tables[features]
    source <- "original saved class-conditional kernel-density grids; no bandwidth change or refitting"
  }
  tab <- do.call(rbind, tables)
  rownames(tab) <- NULL
  if (any(!is.finite(tab$x)) || any(!is.finite(tab$value)) ||
      (type == "conditional_density" && any(tab$value < 0))) sb_abort(
    "Class-conditional curve contains invalid values.", "stormbreaker_invalid_diagnostic")
  structure(list(type = paste0("naive_bayes_", type), table = tab, native_table = native,
    features = features, prior_log_odds = prior_log_odds, learner = id, pipeline_id = fit$pipeline$id,
    source = source, interpretation = if (type == "conditional_density")
      "Fitted class-conditional densities, not posterior class probabilities or variable importance.\nFeature axes use their own units; do not compare density peak heights across features." else
      "Original log likelihood contribution toward class 1; excludes prior log odds and other features.\nNot SHAP, a causal effect, or a validated clinical threshold; original variance smoothing retained."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_naive_bayes_diagnostic <- function(x) {
  feature <- value <- class <- NULL
  tab <- x$table
  tab$feature <- factor(tab$feature, levels = x$features)
  if (x$type == "naive_bayes_conditional_density") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(x, value, colour = class)) +
      ggplot2::geom_line() + ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#D77935")) +
      ggplot2::scale_y_continuous(labels = function(z) sprintf("%.1e", z)) +
      ggplot2::labs(title = "Fitted class-conditional distributions", y = "Density", colour = "Encoded class")
  } else {
    p <- ggplot2::ggplot(tab, ggplot2::aes(x, value)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey65", linetype = 2) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::labs(title = "Gaussian Naive Bayes log likelihood contributions", y = "Log likelihood ratio: class 1 / class 0")
  }
  p + ggplot2::facet_wrap(~feature,
    scales = if (x$type == "naive_bayes_conditional_density") "free" else "free_x",
    ncol = min(3L, length(x$features))) +
    ggplot2::labs(x = "Feature value (original units)") + sb_viz_theme()
}

# ---- centroid-bernoulli ----

sb_predict_centroid_binary <- function(engine, x) {
  # Equivalent to (squared distance to class 0 - distance to class 1)/2,
  # avoiding cancellation of the common x^2 term.
  c0 <- engine$centroids[1L, ]
  c1 <- engine$centroids[2L, ]
  stats::plogis(as.numeric(sweep(x, 2L, (c0 + c1) / 2, "-") %*% (c1 - c0)))
}

sb_fit_bernoulli_binary <- function(x, outcome, params) {
  binary <- x > params$binarize
  counts <- c(sum(outcome == 0L), sum(outcome == 1L))
  probability <- rbind(
    (colSums(binary[outcome == 0L, , drop = FALSE]) + params$alpha) /
      (counts[1L] + 2 * params$alpha),
    (colSums(binary[outcome == 1L, , drop = FALSE]) + params$alpha) /
      (counts[2L] + 2 * params$alpha))
  list(feature_probability = probability, prior = counts / sum(counts))
}

sb_predict_bernoulli_binary <- function(engine, x, params) {
  binary <- 1 * (x > params$binarize)
  p <- engine$feature_probability
  joint <- binary %*% t(log(p)) + (1 - binary) %*% t(log1p(-p))
  joint <- sweep(joint, 2L, log(engine$prior), "+")
  stats::plogis(joint[, 2L] - joint[, 1L])
}

# ---- gaussian-naive-bayes-adapter ----

sb_validate_naive_bayes_var_smoothing <- function(var_smoothing) {
  valid <- is.numeric(var_smoothing) && length(var_smoothing) == 1L &&
    !is.na(var_smoothing) && is.finite(var_smoothing) &&
    var_smoothing >= 0
  if (!valid) {
    sb_abort(
      "`var_smoothing` must be one finite, non-negative number.",
      "stormbreaker_bad_naive_bayes_var_smoothing"
    )
  }
  as.numeric(var_smoothing)
}

sb_population_column_variance <- function(x) {
  centered <- sweep(x, 2L, colMeans(x), "-")
  colMeans(centered^2)
}

sb_fit_gaussian_naive_bayes <- function(x, outcome, var_smoothing) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "A Gaussian Naive Bayes analysis fold must contain both classes.",
      "stormbreaker_naive_bayes_single_class"
    )
  }
  var_smoothing <- sb_validate_naive_bayes_var_smoothing(var_smoothing)
  matrix_x <- as.matrix(x)
  global_variance <- sb_population_column_variance(matrix_x)
  epsilon <- unname(var_smoothing * max(global_variance))
  classes <- c(0L, 1L)
  counts <- vapply(classes, function(class) {
    sum(outcome == class)
  }, integer(1))
  means <- do.call(rbind, lapply(classes, function(class) {
    colMeans(matrix_x[outcome == class, , drop = FALSE])
  }))
  variances <- do.call(rbind, lapply(seq_along(classes), function(index) {
    class_x <- matrix_x[outcome == classes[[index]], , drop = FALSE]
    sb_population_column_variance(class_x) + epsilon
  }))
  rownames(means) <- rownames(variances) <- as.character(classes)
  colnames(means) <- colnames(variances) <- names(x)
  valid_variance <- is.finite(epsilon) &&
    all(is.finite(variances)) && all(variances > 0)
  if (!valid_variance) {
    sb_abort(
      paste0(
        "Gaussian Naive Bayes did not produce one positive, finite ",
        "variance per class and predictor."
      ),
      "stormbreaker_invalid_naive_bayes_variance"
    )
  }
  priors <- counts / length(outcome)
  names(counts) <- names(priors) <- as.character(classes)

  structure(
    list(
      features = names(x),
      classes = classes,
      var_smoothing = var_smoothing,
      epsilon = epsilon,
      counts = counts,
      priors = priors,
      means = means,
      variances = variances
    ),
    class = "stormbreaker_gaussian_naive_bayes"
  )
}

sb_predict_gaussian_naive_bayes <- function(model, x) {
  if (!inherits(model, "stormbreaker_gaussian_naive_bayes")) {
    sb_abort(
      "The supplied model is not a fitted Gaussian Naive Bayes learner.",
      "stormbreaker_bad_naive_bayes_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the Gaussian Naive Bayes schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  matrix_x <- as.matrix(x)
  class_names <- as.character(model$classes)
  log_joint <- vapply(class_names, function(class) {
    centered <- sweep(matrix_x, 2L, model$means[class, ], "-")
    scaled_squared_distance <- sweep(
      centered^2,
      2L,
      model$variances[class, ],
      "/"
    )
    log(model$priors[[class]]) - 0.5 * (
      sum(log(2 * pi * model$variances[class, ])) +
        rowSums(scaled_squared_distance)
    )
  }, numeric(nrow(matrix_x)))
  log_joint <- matrix(
    log_joint,
    nrow = nrow(matrix_x),
    ncol = length(class_names),
    dimnames = list(NULL, class_names)
  )
  probability <- unname(stats::plogis(
    log_joint[, "1"] - log_joint[, "0"]
  ))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      "Gaussian Naive Bayes produced invalid class-one probabilities.",
      "stormbreaker_invalid_naive_bayes_prediction"
    )
  }
  probability
}

sb_extract_gaussian_naive_bayes_parameters <- function(model) {
  if (!inherits(model, "stormbreaker_gaussian_naive_bayes")) {
    sb_abort(
      "The supplied model is not a fitted Gaussian Naive Bayes learner.",
      "stormbreaker_bad_naive_bayes_fit"
    )
  }
  values <- do.call(rbind, lapply(model$classes, function(class) {
    class_name <- as.character(class)
    data.frame(
      feature = model$features,
      class = as.integer(class),
      mean = unname(model$means[class_name, ]),
      variance = unname(model$variances[class_name, ]),
      standard_deviation = unname(sqrt(model$variances[class_name, ])),
      stringsAsFactors = FALSE
    )
  }))
  rownames(values) <- NULL
  list(
    var_smoothing = model$var_smoothing,
    epsilon = model$epsilon,
    priors = data.frame(
      class = model$classes,
      n = as.integer(model$counts[as.character(model$classes)]),
      prior = unname(model$priors[as.character(model$classes)]),
      stringsAsFactors = FALSE
    ),
    values = values
  )
}

sb_gaussian_naive_bayes_learner_adapter <- function(var_smoothing = 1e-9) {
  var_smoothing <- sb_validate_naive_bayes_var_smoothing(var_smoothing)
  params <- list(
    var_smoothing = var_smoothing,
    variance_estimator = "population",
    prior = "empirical"
  )
  list(
    id = "gaussian_naive_bayes",
    params = params,
    fit = function(x, outcome) {
      sb_fit_gaussian_naive_bayes(x, outcome, var_smoothing)
    },
    predict = sb_predict_gaussian_naive_bayes,
    extract_parameters = sb_extract_gaussian_naive_bayes_parameters,
    extract_state = function(model, context = list()) {
      sb_extract_naive_bayes_learner_state(model)
    },
    predict_state = sb_predict_naive_bayes_learner_state
  )
}

# ---- model-bagged-knn ----

# ---- algorithm-bagged-knn ----

#' @rdname sb_knn_neighborhood
#' @export
sb_bagged_knn_diagnostics <- function(fit, type = "bootstrap", newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (s$learner != "bagged_knn_survival" || !is.character(type) ||
      length(type) != 1L || is.na(type) || !type %in% c("bootstrap", "bag_curves")) sb_abort(
    "Choose bootstrap or bag_curves for bagged survival KNN.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  bags <- methods::slot(e, "base_learners")
  train <- methods::slot(e, "train_data")
  if (type == "bootstrap") {
    tab <- do.call(rbind, lapply(seq_along(bags), function(i) {
      index <- methods::slot(bags[[i]], "bootstrap_sample")
      data.frame(bag = i, training_row = seq_len(nrow(train)),
        multiplicity = tabulate(index, nbins = nrow(train)))
    }))
    summary <- do.call(rbind, lapply(split(tab, tab$bag), function(z)
      data.frame(bag = z$bag[1L], sampled_rows = sum(z$multiplicity),
        unique_rows = sum(z$multiplicity > 0), out_of_bag_rows = sum(z$multiplicity == 0))))
  } else {
    if (is.null(newdata)) sb_abort("Supply explicit queries for bag-specific survival curves.",
      "stormbreaker_unavailable_diagnostic")
    x <- sb_select_final_feature_frame(newdata, s$features, "Query data", "query fitted features")
    if (!nrow(x)) sb_abort("Query data must contain at least one row.", "stormbreaker_invalid_diagnostic")
    x <- data.matrix(sb_bagged_knn_internal_data(x)$data)
    time <- methods::slot(e, "timepoints")
    native_predict <- get("predict", asNamespace("bnnSurvival"))
    values <- lapply(bags, function(b) native_predict(b, train, x, time,
      methods::slot(e, "metric"), methods::slot(e, "weighting_function"), methods::slot(e, "k")))
    if (any(vapply(values, function(v) !is.matrix(v) ||
        !identical(dim(v), c(nrow(x), length(time))) || any(!is.finite(v)) ||
        any(v < 0 | v > 1), logical(1)))) sb_abort(
      "Native bag-specific survival curves are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- do.call(rbind, lapply(seq_along(values), function(i)
      data.frame(bag = i, query = rep(seq_len(nrow(x)), each = length(time)),
        time = rep(time, nrow(x)), survival = as.numeric(t(values[[i]])))))
    summary <- stats::aggregate(survival ~ query + time, tab, mean)
  }
  structure(list(type = paste0("bagged_knn_", type), table = tab, summary = summary,
    learner = s$learner, pipeline_id = fit$pipeline$id, sensitive_training_data = type == "bootstrap",
    source = "Original bnnSurvival bootstrap indices and native base-learner predictions",
    interpretation = paste("All original bags retained. No new bootstrap or fitted model.",
      "Grey curves show individual bags; blue is the original ensemble mean, not a confidence band.",
      "Bootstrap coverage is not OOB predictive performance. Explicit queries are not validation.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_bagged_knn_diagnostic <- function(x) {
  bag <- unique_rows <- time <- survival <- query <- NULL
  if (x$type == "bagged_knn_bootstrap") {
    p <- ggplot2::ggplot(x$summary, ggplot2::aes(bag, unique_rows)) +
      ggplot2::geom_point(colour = "#2878A0") +
      ggplot2::labs(title = "Bootstrap coverage", x = "Bag", y = "Unique training rows")
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(time, survival, group = bag)) +
      ggplot2::geom_step(colour = "#A5ADB4", alpha = 0.35, linewidth = 0.3) +
      ggplot2::geom_step(data = x$summary, ggplot2::aes(group = query), colour = "#2878A0", linewidth = 0.9) +
      ggplot2::facet_wrap(~query, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "Bag-specific survival predictions", x = "Time", y = "Survival probability")
  }
  p + sb_viz_theme()
}

# ---- bagged-knn-survival-adapter ----

sb_bagged_knn_survival_learner_ids <- function() {
  "bagged_knn_survival"
}

sb_bagged_knn_survival_available <- function() {
  requireNamespace("bnnSurvival", quietly = TRUE)
}

sb_require_bagged_knn_survival_dependency <- function() {
  if (!sb_bagged_knn_survival_available()) {
    sb_abort(
      paste0(
        "Package `bnnSurvival` is required for bagged k-nearest-neighbor ",
        "survival prediction."
      ),
      "stormbreaker_missing_bnn_survival"
    )
  }
  invisible(TRUE)
}

sb_bagged_knn_survival_params <- function(
  learner = sb_bagged_knn_survival_learner_ids()
) {
  if (!learner %in% sb_bagged_knn_survival_learner_ids()) {
    sb_abort(
      paste0("Unknown bagged kNN survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    num_base_learners = 50L,
    neighbor_rule = "max(1,floor(n_analysis/10))",
    feature_rule = "all_selected_features",
    metric = "mahalanobis",
    distance_weighting = "constant_one",
    replace = TRUE,
    sample_fraction = 1,
    cores = 1L,
    seed = 1L
  )
}

sb_bagged_knn_survival_parameter_semantics <- function(learner) {
  params <- sb_bagged_knn_survival_params(learner)
  paste0(
    "bagged k-nearest-neighbor survival probability with ",
    params$num_base_learners, " bootstrap learners, k=max(1,floor(n/10)), ",
    "all selected features, Mahalanobis distance, constant neighbor weights, ",
    "one-core execution, and seed ", params$seed,
    " fixed before resampling"
  )
}



sb_with_bagged_knn_single_core <- function(code) {
  previous <- getOption("mc.cores")
  on.exit(options(mc.cores = previous), add = TRUE)
  options(mc.cores = 1L)
  force(code)
}

sb_bagged_knn_constant_weight <- function(distance) {
  distance * 0 + 1
}

sb_bagged_knn_internal_names <- function(p) {
  paste0(".stormbreaker_bknn_x", seq_len(p))
}

sb_bagged_knn_internal_data <- function(x) {
  internal <- x
  internal_features <- sb_bagged_knn_internal_names(ncol(x))
  names(internal) <- internal_features
  uses_dummy <- ncol(x) == 1L
  if (uses_dummy) {
    internal$.stormbreaker_bknn_zero <- 0
  }
  list(
    data = internal,
    internal_features = names(internal),
    uses_zero_variance_bridge = uses_dummy
  )
}

sb_bagged_knn_formula_data <- function(inputs) {
  prepared <- sb_bagged_knn_internal_data(inputs$x)
  analysis <- prepared$data
  rownames(analysis) <- NULL
  analysis$.stormbreaker_bknn_time <- inputs$time
  analysis$.stormbreaker_bknn_event <- inputs$event
  response <- paste0(
    "Surv(.stormbreaker_bknn_time, .stormbreaker_bknn_event)"
  )
  formula <- stats::as.formula(
    paste(response, "~", paste(prepared$internal_features, collapse = " + ")),
    env = asNamespace("survival")
  )
  c(prepared, list(analysis = analysis, formula = formula))
}

sb_validate_bagged_knn_engine <- function(engine, internal_features,
                                          effective_k, num_base_learners) {
  valid <- methods::is(engine, "bnnSurvivalEnsemble") &&
    identical(
      as.integer(methods::slot(engine, "num_base_learners")),
      as.integer(num_base_learners)
    ) &&
    identical(as.integer(methods::slot(engine, "k")), as.integer(effective_k)) &&
    identical(
      as.integer(methods::slot(engine, "num_features_per_base_learner")),
      as.integer(length(internal_features))
    ) &&
    identical(
      colnames(methods::slot(engine, "train_data"))[-c(1L, 2L)],
      internal_features
    )
  timepoints <- if (isTRUE(valid)) {
    as.numeric(methods::slot(engine, "timepoints"))
  } else {
    numeric()
  }
  valid <- isTRUE(valid) && length(timepoints) > 0L && !anyNA(timepoints) &&
    all(is.finite(timepoints)) && all(timepoints >= 0) &&
    all(diff(timepoints) > 0)
  if (!isTRUE(valid)) {
    sb_abort(
      "The bagged kNN survival learner returned an invalid native fit.",
      "stormbreaker_invalid_bagged_knn_survival_fit"
    )
  }
  engine
}

sb_validate_bagged_knn_model <- function(model) {
  valid <- inherits(model, "stormbreaker_bagged_knn_survival_model") &&
    is.list(model) &&
    identical(model$learner, sb_bagged_knn_survival_learner_ids()) &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.character(model$internal_features) &&
    length(model$internal_features) %in% c(
      length(model$features),
      length(model$features) + 1L
    ) &&
    is.logical(model$uses_zero_variance_bridge) &&
    length(model$uses_zero_variance_bridge) == 1L &&
    is.numeric(model$effective_k) && length(model$effective_k) == 1L &&
    !is.na(model$effective_k) && model$effective_k >= 1L &&
    is.list(model$params)
  if (!isTRUE(valid)) {
    sb_abort(
      "The supplied model is not a valid fitted bagged kNN survival learner.",
      "stormbreaker_bad_bagged_knn_survival_fit"
    )
  }
  sb_validate_bagged_knn_engine(
    model$engine,
    model$internal_features,
    model$effective_k,
    model$params$num_base_learners
  )
  model
}

sb_fit_bagged_knn_survival <- function(x, time, event, learner, params) {
  sb_require_bagged_knn_survival_dependency()
  inputs <- sb_validate_classical_survival_inputs(x, time, event, learner)
  prepared <- sb_bagged_knn_formula_data(inputs)
  effective_k <- max(1L, floor(nrow(inputs$x) / 10L))
  engine <- sb_with_bagged_knn_single_core(
    withr::with_seed(
      params$seed,
      bnnSurvival::bnnSurvival(
        formula = prepared$formula,
        data = prepared$analysis,
        k = effective_k,
        num_base_learners = params$num_base_learners,
        num_features_per_base_learner = length(prepared$internal_features),
        metric = params$metric,
        weighting_function = sb_bagged_knn_constant_weight,
        replace = params$replace,
        sample_fraction = params$sample_fraction
      )
    )
  )
  model <- structure(
    list(
      learner = learner,
      features = names(inputs$x),
      internal_features = prepared$internal_features,
      uses_zero_variance_bridge = prepared$uses_zero_variance_bridge,
      effective_k = as.integer(effective_k),
      params = params,
      engine = engine
    ),
    class = "stormbreaker_bagged_knn_survival_model"
  )
  sb_validate_bagged_knn_model(model)
}

sb_bagged_knn_prediction_data <- function(model, x) {
  model <- sb_validate_bagged_knn_model(model)
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      paste0(
        "Prediction columns do not match the fitted bagged kNN survival ",
        "learner schema."
      ),
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prepared <- sb_bagged_knn_internal_data(x)
  if (!identical(prepared$internal_features, model$internal_features) ||
    !identical(
      prepared$uses_zero_variance_bridge,
      model$uses_zero_variance_bridge
    )) {
    sb_abort(
      "The internal bagged kNN survival predictor schema is invalid.",
      "stormbreaker_bad_bagged_knn_survival_fit"
    )
  }
  list(model = model, x = x, internal = prepared$data)
}

sb_bagged_knn_distance_function <- function(k) {
  inverse <- product <- transpose_product <- diagonal <- NULL
  bound_scale <- train_magnitude <- NULL
  function(x, center, cov, ...) {
    x <- as.matrix(x)
    if (is.null(inverse)) {
      inverse <<- solve(cov, ...)
      product <<- x %*% inverse
      transpose_product <<- x %*% t(inverse)
      diagonal <<- rowSums(product * x)
      train_magnitude <<- max(abs(x))
      bound_scale <<- 128 * .Machine$double.eps * ncol(x)^3 *
        max(abs(inverse))
    }
    # Expanded quadratic products screen neighbors only. Retain both cross
    # terms: a numerically inverted covariance need not be exactly symmetric.
    values <- diagonal - as.vector(product %*% center) -
      as.vector(transpose_product %*% center) +
      sum((center %*% inverse) * center)
    # Cancellation can change near-ties. This conservative absolute rounding
    # envelope covers both the expanded and centered quadratic evaluations;
    # the normal-minimum term also covers underflow. Recompute every possible
    # top-k neighbor with the original centered native expression below.
    bound <- bound_scale * (train_magnitude + max(abs(center)))^2 +
      128 * .Machine$double.xmin * ncol(x)^3
    exact <- function(rows = NULL) {
      xx <- if (is.null(rows)) x else x[rows, , drop = FALSE]
      stats::mahalanobis(xx, center, inverse, inverted = TRUE)
    }
    if (any(!is.finite(values)) || !is.finite(bound) || k >= nrow(x)) {
      return(exact())
    }
    index <- order(values)
    selected <- which(values <= values[index[k]] + 2 * bound)
    values[selected] <- exact(selected)
    stats::setNames(values, rownames(x))
  }
}

sb_predict_bagged_knn_engine <- function(object, test_data) {
  ns <- asNamespace("bnnSurvival")
  # The native 0.1.5 S4 methods store their implementation as a .local closure.
  # Use local copies of those closures, not namespace/S4-table replacement or
  # a second implementation of native sorting, KM and ensemble aggregation.
  base_predict <- body(methods::getMethod("predict",
    "bnnSurvivalBaseLearner", where = ns))[[2L]][[3L]]
  predict_base <- function(object, train_data, test_data, timepoints,
                           metric, weighting_function, k) {
    # A fresh closure belongs only to this bootstrap learner and this call.
    distance <- sb_bagged_knn_distance_function(k)
    local_base <- base_predict
    environment(local_base) <- list2env(list(mahalanobis = distance),
      parent = environment(base_predict))
    local_base(object, train_data, test_data, timepoints,
      metric, weighting_function, k)
  }
  native_ensemble <- body(methods::getMethod("predict",
    "bnnSurvivalEnsemble", where = ns))[[2L]][[3L]]
  local_ensemble <- native_ensemble
  environment(local_ensemble) <- list2env(list(predict = predict_base),
    parent = environment(native_ensemble))
  local_ensemble(object, test_data)
}

sb_predict_bagged_knn_survival_curves <- function(model, x) {
  prepared <- sb_bagged_knn_prediction_data(model, x)
  prediction_object <- sb_with_bagged_knn_single_core(
    sb_predict_bagged_knn_engine(prepared$model$engine, prepared$internal)
  )
  survival <- as.matrix(bnnSurvival::predictions(prediction_object))
  timepoints <- as.numeric(bnnSurvival::timepoints(prediction_object))
  valid <- nrow(survival) == nrow(prepared$x) &&
    ncol(survival) == length(timepoints) && length(timepoints) > 0L &&
    !anyNA(timepoints) && all(is.finite(timepoints)) &&
    all(timepoints >= 0) && all(diff(timepoints) > 0) &&
    !anyNA(survival) && all(is.finite(survival)) &&
    all(survival >= 0 & survival <= 1) &&
    all(apply(survival, 1L, function(values) {
      all(diff(values) <= 1e-10)
    }))
  if (!isTRUE(valid)) {
    sb_abort(
      "The bagged kNN learner produced invalid survival curves.",
      "stormbreaker_invalid_bagged_knn_survival_prediction"
    )
  }
  interval <- diff(c(0, timepoints))
  previous <- if (ncol(survival) == 1L) {
    matrix(1, nrow = nrow(survival), ncol = 1L)
  } else {
    cbind(1, survival[, -ncol(survival), drop = FALSE])
  }
  risk <- -as.vector(previous %*% interval)
  if (length(risk) != nrow(prepared$x) || anyNA(risk) ||
    any(!is.finite(risk))) {
    sb_abort(
      "The bagged kNN learner produced invalid restricted-mean risk scores.",
      "stormbreaker_invalid_bagged_knn_survival_prediction"
    )
  }
  list(
    risk = unname(risk),
    time = timepoints,
    survival = survival
  )
}

sb_predict_bagged_knn_survival_risk <- function(model, x) {
  sb_predict_bagged_knn_survival_curves(model, x)$risk
}

sb_bagged_knn_direct_survival_source <- function() {
  "bnn_survival_bagged_knn_direct"
}

sb_predict_bagged_knn_survival_probability <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon,
  prediction_details = NULL
) {
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      paste0(
        "Analysis columns do not match the fitted bagged kNN survival ",
        "learner schema."
      ),
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  curves <- if (is.null(prediction_details)) {
    sb_predict_bagged_knn_survival_curves(model, x_assessment)
  } else {
    prediction_details
  }
  eligible <- which(curves$time <= prediction_horizon)
  probability <- if (length(eligible) == 0L) {
    rep(1, length(curves$risk))
  } else {
    unname(curves$survival[, max(eligible)])
  }
  list(
    prediction_horizon = prediction_horizon,
    source = sb_bagged_knn_direct_survival_source(),
    ties = "not_applicable",
    probability = probability
  )
}

sb_extract_bagged_knn_survival_state <- function(model, context = list()) {
  model <- sb_validate_bagged_knn_model(model)
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    horizon <- sb_validate_prediction_horizon(
      context$prediction_horizon,
      context$time
    )
    survival_probability <- list(
      prediction_horizon = horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      source = sb_bagged_knn_direct_survival_source(),
      ties = "not_applicable"
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      internal_features = model$internal_features,
      uses_zero_variance_bridge = model$uses_zero_variance_bridge,
      effective_k = model$effective_k,
      params = model$params,
      engine_serialization = "R_xdr_v3_bnnSurvival_prediction_engine",
      engine_raw = serialize(
        model$engine,
        connection = NULL,
        ascii = FALSE,
        xdr = TRUE,
        version = 3L
      ),
      training_data_retention =
        "required_by_bagged_knn_prediction_state",
      risk_reference = "bnnSurvival_native_training_time_grid",
      risk_formula = "negative_step_RMST_over_native_training_time_grid",
      baseline_hazard = NULL,
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_bagged_knn_survival_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_bagged_knn_survival_model <- function(state) {
  if (!is.raw(state$engine_raw) || length(state$engine_raw) == 0L) {
    sb_bad_learner_state(
      "The saved bagged kNN survival prediction engine is invalid."
    )
  }
  engine <- tryCatch(
    unserialize(state$engine_raw),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved bagged kNN survival engine could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_bagged_knn_survival_state",
        parent = error
      )
    }
  )
  model <- structure(
    list(
      learner = state$learner,
      features = state$features,
      internal_features = state$internal_features,
      uses_zero_variance_bridge = state$uses_zero_variance_bridge,
      effective_k = state$effective_k,
      params = state$params,
      engine = engine
    ),
    class = "stormbreaker_bagged_knn_survival_model"
  )
  sb_validate_bagged_knn_model(model)
}

sb_validate_bagged_knn_probability_state <- function(state) {
  probability <- state$survival_probability
  valid <- is.list(probability) &&
    is.numeric(probability$prediction_horizon) &&
    length(probability$prediction_horizon) == 1L &&
    !is.na(probability$prediction_horizon) &&
    is.finite(probability$prediction_horizon) &&
    probability$prediction_horizon > 0 &&
    is.character(probability$time_unit) &&
    length(probability$time_unit) == 1L &&
    !is.na(probability$time_unit) && nzchar(probability$time_unit) &&
    identical(
      probability$source,
      sb_bagged_knn_direct_survival_source()
    ) && identical(probability$ties, "not_applicable")
  if (!isTRUE(valid)) {
    sb_bad_learner_state(
      "The saved bagged kNN direct survival-probability state is invalid."
    )
  }
  probability
}

sb_predict_bagged_knn_survival_state <- function(state, x) {
  valid <- inherits(state, "stormbreaker_bagged_knn_survival_state") &&
    is.list(state) &&
    identical(state$learner, sb_bagged_knn_survival_learner_ids()) &&
    identical(
      state$training_data_retention,
      "required_by_bagged_knn_prediction_state"
    )
  if (!isTRUE(valid)) {
    sb_bad_learner_state("The saved bagged kNN survival state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, state$learner)
  sb_require_bagged_knn_survival_dependency()
  model <- sb_restore_bagged_knn_survival_model(state)
  curves <- sb_predict_bagged_knn_survival_curves(model, x)
  result <- list(risk = curves$risk)
  if (!is.null(state$survival_probability)) {
    probability <- sb_validate_bagged_knn_probability_state(state)
    eligible <- which(curves$time <= probability$prediction_horizon)
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- if (length(eligible) == 0L) {
      rep(1, length(curves$risk))
    } else {
      unname(curves$survival[, max(eligible)])
    }
  }
  result
}

sb_bagged_knn_survival_learner_adapter <- function(
  learner,
  params = sb_bagged_knn_survival_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_bagged_knn_survival(x, time, event, learner, params)
    },
    predict = sb_predict_bagged_knn_survival_risk,
    predict_details = sb_predict_bagged_knn_survival_curves,
    predict_survival = sb_predict_bagged_knn_survival_probability,
    extract_state = sb_extract_bagged_knn_survival_state,
    predict_state = sb_predict_bagged_knn_survival_state
  )
}

# ---- model-mclust ----

# ---- algorithm-mclust-da ----

#' @rdname sb_knn_neighborhood
#' @export
sb_mclust_da_diagnostics <- function(fit, type = "centers") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% sb_mclust_da_learner_ids() || !is.character(type) ||
      length(type) != 1L || is.na(type) || !type %in% c("centers", "correlations", "mixture_weights")) sb_abort(
    "Choose centers, correlations or mixture_weights for EDDA/MclustDA.", "stormbreaker_unsupported_diagnostic")
  e <- s$engine
  sb_validate_mclust_da_prediction_engine(e, s$params)
  p <- length(s$features)
  tab <- do.call(rbind, lapply(names(e$models), function(cls) {
    model <- e$models[[cls]]
    par <- model$parameters
    g <- length(par$pro)
    means <- matrix(par$mean, nrow = p, ncol = g)
    if (type == "centers") return(data.frame(class = cls,
      component = rep(seq_len(g), each = p), feature = rep(s$features, g), mean = as.numeric(means)))
    if (type == "mixture_weights") return(data.frame(class = cls, component = seq_len(g),
      weight = as.numeric(par$pro), class_prior = unname(e$prop[cls]), model = model$modelName))
    covariance <- par$variance$sigma
    if (p == 1L) covariance <- array(rep(par$variance$sigmasq, length.out = g), c(1L, 1L, g))
    if (!is.array(covariance) || !identical(as.integer(dim(covariance)), as.integer(c(p, p, g))) ||
        any(!is.finite(covariance))) sb_abort("Original component covariances are invalid.", "stormbreaker_invalid_diagnostic")
    do.call(rbind, lapply(seq_len(g), function(j) {
      m <- matrix(covariance[, , j], p, p)
      if (any(diag(m) <= 0)) sb_abort("Original component variances must be positive.", "stormbreaker_invalid_diagnostic")
      correlation <- stats::cov2cor(m)
      data.frame(class = cls, component = j, feature = rep(s$features, p),
        feature2 = rep(s$features, each = p), correlation = as.numeric(correlation))
    }))
  }))
  structure(list(type = paste0("mclust_da_", type), table = tab,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original saved mclust discriminant mixture parameters",
    interpretation = "Latent components are not clinical subgroups. Correlations are fitted component correlations, not raw sample correlations. Mixing weights are conditional on class. No refitting or BIC reconstruction."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_mclust_da_diagnostic <- function(x) {
  mean <- feature <- feature2 <- correlation <- component <- weight <- NULL
  if (x$type == "mclust_da_centers") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(mean, feature)) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::facet_wrap(~class + component, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "Mixture component centers", x = "Mean", y = NULL)
  } else if (x$type == "mclust_da_correlations") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(feature2, feature, fill = correlation)) +
      ggplot2::geom_tile(colour = "white", linewidth = .2) +
      ggplot2::facet_wrap(~class + component, labeller = ggplot2::label_both) +
      ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#B86236", limits = c(-1, 1),
        guide = ggplot2::guide_colourbar(display = "rectangles")) +
      ggplot2::labs(title = "Fitted component correlations", x = NULL, y = NULL, fill = "Correlation")
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(factor(component), weight)) +
      ggplot2::geom_col(fill = "#2878A0", width = .6) + ggplot2::facet_wrap(~class, labeller = ggplot2::label_both) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = "Within-class mixing weights", x = "Component", y = "Weight")
  }
  p <- p + sb_viz_theme()
  if (x$type == "mclust_da_correlations") p <- p + ggplot2::coord_fixed() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank())
  attr(p, "stormbreaker_suggested_size") <- c(width = 12, height = 6)
  p
}

# ---- mclust-da-adapters ----

sb_mclust_da_learner_ids <- function() {
  c("edda_binary_bic", "mclustda_binary_bic")
}

sb_mclust_da_multivariate_models <- function() {
  c(
    "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE",
    "VEE", "EVE", "VVE", "EEV", "VEV", "EVV", "VVV"
  )
}

sb_mclust_da_params <- function(learner) {
  if (!learner %in% sb_mclust_da_learner_ids()) {
    sb_abort(
      paste0("Unknown mclust discriminant learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    model_type = if (identical(learner, "edda_binary_bic")) {
      "EDDA"
    } else {
      "MclustDA"
    },
    components = if (identical(learner, "edda_binary_bic")) 1L else 1:3,
    univariate_covariance_models = c("E", "V"),
    multivariate_covariance_models = sb_mclust_da_multivariate_models(),
    criterion = "BIC",
    class_prior = "analysis_fold_empirical",
    parameter_prior = "none",
    scale = "original",
    em_control = list(
      eps = .Machine$double.eps,
      tol = c(1e-5, sqrt(.Machine$double.eps)),
      itmax = rep.int(.Machine$integer.max, 2L),
      equalPro = FALSE
    ),
    initialization = if (identical(learner, "edda_binary_bic")) {
      list(method = "known_class_labels")
    } else {
      list(
        method = "native_hierarchical",
        subset_max = 2000L,
        hc_model_name = "VVV",
        hc_use = "SVD",
        seed = 1L
      )
    },
    native_candidate_warnings = FALSE,
    all_model_failure = "typed"
  )
}

sb_mclust_available <- function() {
  requireNamespace("mclust", quietly = TRUE)
}

sb_require_mclust <- function() {
  if (!sb_mclust_available()) {
    sb_abort(
      "Package `mclust` is required for model-based discriminant learners.",
      "stormbreaker_missing_mclust"
    )
  }
  invisible(TRUE)
}

sb_mclust_da_active_models <- function(n_features, params) {
  if (identical(n_features, 1L)) {
    params$univariate_covariance_models
  } else {
    params$multivariate_covariance_models
  }
}

sb_validate_mclust_da_inputs <- function(x, outcome, learner) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      paste0("A `", learner, "` analysis fold must contain both classes."),
      "stormbreaker_mclust_da_single_class"
    )
  }
  list(x = x, outcome = outcome)
}

sb_mclust_da_fit_call <- function(matrix_x, grouping, params, learner) {
  # MclustDA resolves native fitting helpers in its caller environment.
  call_environment <- new.env(parent = asNamespace("mclust"))
  call_environment$data <- matrix_x
  call_environment$grouping <- grouping
  call_environment$components <- params$components
  call_environment$covariance_models <- params$active_covariance_models
  call_environment$model_type <- params$model_type
  call_environment$control <- do.call(mclust::emControl, params$em_control)
  call_environment$native_candidate_warnings <-
    params$native_candidate_warnings
  if (identical(params$model_type, "MclustDA")) {
    initialization <- params$initialization
    fixed_options <- list(
      subset = initialization$subset_max,
      hcModelName = initialization$hc_model_name,
      hcUse = initialization$hc_use
    )
    call_environment$mclust.options <- function(option) {
      if (option %in% names(fixed_options)) {
        return(fixed_options[[option]])
      }
      mclust::mclust.options(option)
    }
    local_mclust_bic <- mclust::mclustBIC
    environment(local_mclust_bic) <- call_environment
    call_environment$mclustBIC <- local_mclust_bic
  }

  tryCatch(
    withCallingHandlers(
      {
        fit_call <- quote(MclustDA(
          data = data,
          class = grouping,
          G = components,
          modelNames = covariance_models,
          modelType = model_type,
          prior = NULL,
          control = control,
          warn = native_candidate_warnings,
          verbose = FALSE
        ))
        if (identical(params$model_type, "MclustDA")) {
          withr::with_seed(
            params$initialization$seed,
            eval(fit_call, envir = call_environment)
          )
        } else {
          eval(fit_call, envir = call_environment)
        }
      },
      warning = function(warning) {
        sb_abort(
          paste0(
            "Learner `", learner, "` raised a fitting warning: ",
            conditionMessage(warning)
          ),
          "stormbreaker_mclust_da_fit_warning",
          parent = warning
        )
      }
    ),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "Learner `", learner, "` could not fit its mclust model: ",
          conditionMessage(error)
        ),
        "stormbreaker_mclust_da_fit_error",
        parent = error
      )
    }
  )
}

sb_validate_mclust_da_search_grid <- function(engine, params) {
  covariance_models <- if (!is.null(params$active_covariance_models)) {
    params$active_covariance_models
  } else {
    sb_mclust_da_active_models(as.integer(engine$d), params)
  }
  if (identical(params$model_type, "EDDA")) {
    small_class_models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
    native_models <- unique(unlist(lapply(engine$models, function(model) {
      if (model$n <= engine$d) {
        intersect(covariance_models, small_class_models)
      } else {
        covariance_models
      }
    }), use.names = FALSE))
    valid <- length(native_models) == length(covariance_models) &&
      setequal(native_models, covariance_models)
  } else {
    expected_components <- as.character(params$components)
    valid <- inherits(engine, "MclustDA") && is.list(engine$models) &&
      length(engine$models) == 2L
    if (valid) {
      valid <- all(vapply(engine$models, function(model) {
        bic <- model$BIC
        is.matrix(bic) &&
          identical(rownames(bic), expected_components) &&
          identical(colnames(bic), covariance_models) &&
          any(is.finite(bic))
      }, logical(1)))
    }
  }
  if (!valid) {
    sb_abort(
      paste0(
        "mclust reduced or did not return the predeclared BIC search grid; ",
        "the learner will not silently fit a narrower specification."
      ),
      "stormbreaker_mclust_da_search_reduction"
    )
  }
  invisible(TRUE)
}

sb_validate_mclust_da_engine <- function(engine, params) {
  expected_classes <- c("0", "1")
  valid <- inherits(engine, "MclustDA") &&
    identical(engine$type, params$model_type) &&
    is.numeric(engine$prop) && length(engine$prop) == 2L &&
    identical(names(engine$prop), expected_classes) &&
    !anyNA(engine$prop) && all(is.finite(engine$prop)) &&
    all(engine$prop > 0) && abs(sum(engine$prop) - 1) < 1e-12 &&
    is.list(engine$models) &&
    identical(names(engine$models), expected_classes) &&
    is.numeric(engine$bic) && length(engine$bic) == 1L &&
    !is.na(engine$bic) && is.finite(engine$bic)
  if (valid) {
    allowed_names <- c(
      params$active_covariance_models,
      "X", "XII", "XXI", "XXX"
    )
    valid <- all(vapply(engine$models, function(model) {
      is.list(model) &&
        is.character(model$modelName) && length(model$modelName) == 1L &&
        !is.na(model$modelName) && model$modelName %in% allowed_names &&
        is.numeric(model$G) && length(model$G) == 1L &&
        !is.na(model$G) && model$G %in% params$components &&
        is.list(model$parameters) &&
        is.numeric(model$parameters$pro) &&
        length(model$parameters$pro) == model$G &&
        !anyNA(model$parameters$pro) &&
        all(is.finite(model$parameters$pro)) &&
        all(model$parameters$pro > 0)
    }, logical(1)))
  }
  if (valid && identical(params$model_type, "EDDA")) {
    valid <- length(unique(vapply(
      engine$models,
      `[[`,
      character(1),
      "modelName"
    ))) == 1L
  }
  if (!valid) {
    sb_abort(
      "The fitted mclust discriminant model is incomplete or invalid.",
      "stormbreaker_incomplete_mclust_da_fit"
    )
  }
  sb_validate_mclust_da_search_grid(engine, params)
  engine
}

sb_mclust_da_audit <- function(engine, params) {
  selected <- data.frame(
    class = names(engine$models),
    model_name = vapply(engine$models, `[[`, character(1), "modelName"),
    components = as.integer(vapply(engine$models, `[[`, numeric(1), "G")),
    stringsAsFactors = FALSE
  )
  candidate_summary <- if (identical(params$model_type, "MclustDA")) {
    data.frame(
      class = names(engine$models),
      requested = vapply(
        engine$models,
        function(model) length(model$BIC),
        integer(1)
      ),
      estimable = vapply(
        engine$models,
        function(model) sum(is.finite(model$BIC)),
        integer(1)
      ),
      initialization_subset_size = vapply(
        engine$models,
        function(model) {
          initialization <- attr(model$BIC, "initialization", exact = TRUE)
          length(initialization$subset)
        },
        integer(1)
      ),
      stringsAsFactors = FALSE
    )
  } else {
    classes <- names(engine$models)
    data.frame(
      class = classes,
      requested = rep.int(
        length(params$active_covariance_models),
        length(classes)
      ),
      estimable = rep.int(NA_integer_, length(classes)),
      initialization_subset_size = rep.int(NA_integer_, length(classes)),
      stringsAsFactors = FALSE
    )
  }
  list(
    class_prior = engine$prop,
    selected = selected,
    candidate_summary = candidate_summary
  )
}

sb_mclust_da_prediction_engine <- function(engine) {
  models <- lapply(engine$models, function(model) {
    model[c("modelName", "parameters")]
  })
  names(models) <- names(engine$models)
  structure(
    list(type = engine$type, prop = engine$prop, models = models),
    class = "MclustDA"
  )
}

sb_validate_mclust_da_prediction_engine <- function(engine, params) {
  expected_classes <- c("0", "1")
  valid <- inherits(engine, "MclustDA") &&
    identical(engine$type, params$model_type) &&
    is.numeric(engine$prop) && length(engine$prop) == 2L &&
    identical(names(engine$prop), expected_classes) &&
    !anyNA(engine$prop) && all(is.finite(engine$prop)) &&
    all(engine$prop > 0) && abs(sum(engine$prop) - 1) < 1e-12 &&
    is.list(engine$models) &&
    identical(names(engine$models), expected_classes) &&
    all(vapply(engine$models, function(model) {
      is.list(model) &&
        is.character(model$modelName) && length(model$modelName) == 1L &&
        !is.na(model$modelName) && is.list(model$parameters)
    }, logical(1)))
  if (!valid) {
    sb_abort(
      "The saved mclust discriminant engine is incomplete or invalid.",
      "stormbreaker_bad_mclust_da_model"
    )
  }
  engine
}

sb_new_mclust_da_model <- function(learner, features, params, engine, audit) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      engine = engine,
      audit = audit
    ),
    class = "stormbreaker_mclust_da"
  )
}

sb_fit_mclust_da <- function(x, outcome, learner, params) {
  sb_require_mclust()
  inputs <- sb_validate_mclust_da_inputs(x, outcome, learner)
  active_params <- params
  active_params$active_covariance_models <- sb_mclust_da_active_models(
    ncol(inputs$x),
    params
  )
  engine <- sb_mclust_da_fit_call(
    data.matrix(inputs$x),
    factor(inputs$outcome, levels = c(0L, 1L)),
    active_params,
    learner
  )
  engine <- sb_validate_mclust_da_engine(engine, active_params)
  audit <- sb_mclust_da_audit(engine, active_params)
  engine <- sb_mclust_da_prediction_engine(engine)
  sb_validate_mclust_da_prediction_engine(engine, active_params)
  sb_new_mclust_da_model(
    learner,
    names(inputs$x),
    active_params,
    engine,
    audit
  )
}

sb_validate_mclust_da_model_x <- function(model, x) {
  valid <- inherits(model, "stormbreaker_mclust_da") &&
    is.character(model$learner) && length(model$learner) == 1L &&
    model$learner %in% sb_mclust_da_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.list(model$params) && is.list(model$audit)
  if (!valid) {
    sb_abort(
      "The fitted mclust discriminant model is incomplete or invalid.",
      "stormbreaker_bad_mclust_da_model"
    )
  }
  sb_validate_mclust_da_prediction_engine(model$engine, model$params)
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted mclust model.",
      "stormbreaker_mclust_da_schema_mismatch"
    )
  }
  x
}

sb_predict_mclust_da <- function(model, x) {
  sb_require_mclust()
  x <- sb_validate_mclust_da_model_x(model, x)
  prediction <- tryCatch(
    withCallingHandlers(
      stats::predict(model$engine, newdata = data.matrix(x)),
      warning = function(warning) {
        sb_abort(
          paste0(
            "Learner `", model$learner, "` raised a prediction warning: ",
            conditionMessage(warning)
          ),
          "stormbreaker_mclust_da_prediction_warning",
          parent = warning
        )
      }
    ),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "Learner `", model$learner, "` could not produce probabilities: ",
          conditionMessage(error)
        ),
        "stormbreaker_mclust_da_prediction_error",
        parent = error
      )
    }
  )
  probability <- if (is.list(prediction) && is.matrix(prediction$z) &&
      "1" %in% colnames(prediction$z)) {
    unname(prediction$z[, "1"])
  } else {
    numeric()
  }
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) && all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid probabilities."),
      "stormbreaker_invalid_binary_prediction"
    )
  }
  probability
}

sb_extract_mclust_da_learner_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_mclust_da")) {
    sb_bad_learner_state("The fitted mclust discriminant state is invalid.")
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      engine = model$engine,
      audit = model$audit
    ),
    class = c(
      "stormbreaker_mclust_da_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_mclust_da_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
      !learner %in% sb_mclust_da_learner_ids()) {
    sb_bad_learner_state("The saved mclust discriminant state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_new_mclust_da_model(
    learner,
    state$features,
    state$params,
    state$engine,
    state$audit
  )
  list(probability = sb_predict_mclust_da(model, x))
}

sb_mclust_da_learner_adapter <- function(
    learner,
    params = sb_mclust_da_params(learner)) {
  if (!learner %in% sb_mclust_da_learner_ids()) {
    sb_abort(
      paste0("Unknown mclust discriminant learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_mclust_da(x, outcome, learner, params)
    },
    predict = sb_predict_mclust_da,
    extract_state = sb_extract_mclust_da_learner_state,
    predict_state = sb_predict_mclust_da_learner_state
  )
}

# ---- diagnostics-neighbors ----

# ---- algorithm-knn ----

#' @rdname sb_knn_neighborhood
#' @export
sb_knn_neighborhood <- function(fit, newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% c("knn_binary", "kernel_weighted_knn_binary")) sb_abort(
    "This neighborhood view requires binary KNN.", "stormbreaker_unsupported_diagnostic")
  if (is.null(newdata)) sb_abort("Supply explicit query data for a nearest-neighbor explanation.",
    "stormbreaker_unavailable_diagnostic")
  x <- sb_select_final_feature_frame(newdata, s$features, "Query data", "query fitted features")
  if (!nrow(x)) sb_abort("Query data must contain at least one row.", "stormbreaker_invalid_diagnostic")
  e <- unserialize(s$engine_raw)
  n <- nrow(x)
  k <- e$k
  if (s$learner == "knn_binary") {
    z <- sweep(sweep(as.matrix(x), 2L, e$center, "-"), 2L, e$scale, "/")
    native <- lapply(seq_len(n), function(i) {
      distance <- rowSums(sweep(e$training, 2L, z[i, ], "-")^2)
      index <- order(distance, seq_along(distance))[seq_len(k)]
      list(index = index, distance = sqrt(distance[index]))
    })
    index <- do.call(rbind, lapply(native, `[[`, "index"))
    distance <- do.call(rbind, lapply(native, `[[`, "distance"))
    weight <- matrix(1 / k, n, k)
    classes <- matrix(as.character(e$outcome[index]), n, k)
    probability <- rowMeans(classes == "1")
    metric <- "Euclidean distance in original training-standardized feature space"
  } else {
    training <- e$training
    outcome_name <- ".stormbreaker_outcome"
    while (outcome_name %in% names(training)) outcome_name <- paste0(outcome_name, "_")
    training[[outcome_name]] <- e$outcome
    rhs <- paste(vapply(names(x), sb_quote_name, ""), collapse = " + ")
    formula <- stats::as.formula(paste(sb_quote_name(outcome_name), "~", rhs), env = baseenv())
    native <- kknn::kknn(formula, train = training, test = x, k = k,
      distance = e$distance, kernel = e$kernel, scale = e$scale)
    index <- native$C
    distance <- native$D
    weight <- native$W
    classes <- native$CL
    probability <- native$prob[, "1"]
    metric <- "Native kknn distance under original distance/kernel/scaling settings"
  }
  if (any(vapply(list(index, distance, weight, classes), function(z)
      !is.matrix(z) || !identical(dim(z), c(n, as.integer(k))), logical(1))) ||
      any(!is.finite(index)) || any(index != floor(index)) || any(index < 1 | index > nrow(e$training)) ||
      any(!is.finite(distance)) || any(distance < 0) || any(!is.finite(weight)) ||
      any(weight < 0) || any(rowSums(weight) <= 0) || any(!classes %in% c("0", "1"))) sb_abort(
    "Original nearest-neighbor result is invalid.", "stormbreaker_invalid_diagnostic")
  normalized <- weight / rowSums(weight)
  tab <- data.frame(query = rep(seq_len(n), each = k), rank = rep(seq_len(k), n),
    training_row = as.integer(t(index)), distance = as.numeric(t(distance)),
    class = as.character(t(classes)), native_weight = as.numeric(t(weight)),
    voting_weight = as.numeric(t(normalized)))
  structure(list(type = "knn_neighborhood", table = tab,
    predictions = data.frame(query = seq_len(n), probability = as.numeric(probability)),
    settings = s$params,
    learner = s$learner, pipeline_id = fit$pipeline$id, sensitive_training_data = TRUE,
    source = metric, interpretation = paste("Explicit query rows; original neighbor selection and voting.",
      "Native training row indices and classes are sensitive. Self-neighbors are not removed.",
      "Not validation, a tuning curve, causal effects or a patient prediction card.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_knn_neighborhood <- function(x) {
  distance <- rank <- class <- voting_weight <- query <- NULL
  p <- ggplot2::ggplot(x$table, ggplot2::aes(rank, distance, colour = class, size = voting_weight)) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::facet_wrap(~query, labeller = ggplot2::label_both) +
    ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#B86236")) +
    ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
    ggplot2::labs(title = "Prediction neighborhoods", x = "Neighbor rank", y = "Distance",
      colour = "Class", size = "Voting weight") + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}

# ---- algorithm-loclda ----

#' @rdname sb_knn_neighborhood
#' @export
sb_loclda_neighborhood <- function(fit, newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (s$learner != "loclda_binary") sb_abort(
    "This diagnostic requires local LDA.", "stormbreaker_unsupported_diagnostic")
  if (is.null(newdata)) sb_abort("Supply explicit query rows for local-LDA neighborhoods.", "stormbreaker_unavailable_diagnostic")
  if (!requireNamespace("klaR", quietly = TRUE)) sb_abort("Local-LDA diagnostics require klaR.", "stormbreaker_unavailable_diagnostic")
  e <- unserialize(s$engine_raw)
  x <- as.matrix(sb_select_final_feature_frame(newdata, s$features, "Query data", "query fitted features"))
  if (!nrow(x) || !inherits(e, "loclda") || !identical(colnames(e$learn), s$features)) sb_abort(
    "Original local-LDA data or query is invalid.", "stormbreaker_invalid_diagnostic")
  tab <- sb_loclda_neighbor_table(e, x)
  probability <- sb_predict_loclda_binary(e, x)
  structure(list(type = "loclda_neighborhood", table = tab,
    predictions = data.frame(query = seq_len(nrow(x)), probability = as.numeric(probability)),
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original learn/grouping/k/weight.func and native predict.loclda",
    interpretation = "Distance-radius neighbors include ties. Classes with at most one neighbor are excluded by the native prediction rule. Weights are local prediction weights, not feature importance. Reference indices and labels may be sensitive. No tuning or extra model fit."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_loclda_neighbor_table <- function(e, x) {
  do.call(rbind, lapply(seq_len(nrow(x)), function(i) {
    distances <- sqrt(rowSums(sweep(e$learn, 2L, x[i, ], "-")^2))
    radius <- sort(distances)[e$k]
    if (!is.finite(radius) || radius <= 0) sb_abort(
      "Native local-LDA neighborhood radius is zero or invalid.", "stormbreaker_unavailable_diagnostic")
    selected <- which(distances <= radius)
    classes <- as.character(e$grouping[selected])
    counts <- table(factor(classes, levels = e$lev))
    weights <- vapply(distances[selected] / radius, e$weight.func, numeric(1))
    used <- counts[classes] > 1L
    data.frame(query = i, reference_row = selected, class = classes,
      distance = distances[selected], radius = radius, weight = weights,
      used = as.logical(used), effective_weight = weights * as.logical(used))
  }))
}

sb_plot_loclda_neighborhood <- function(x) {
  distance <- weight <- class <- used <- NULL
  p <- ggplot2::ggplot(x$table, ggplot2::aes(distance, weight, colour = class, shape = used)) +
    ggplot2::geom_point(size = 2, alpha = .8) +
    ggplot2::facet_wrap(~query, labeller = ggplot2::label_both) +
    ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#B86236")) +
    ggplot2::scale_shape_manual(values = c("FALSE" = 4, "TRUE" = 16), labels = c("FALSE" = "Excluded", "TRUE" = "Used")) +
    ggplot2::labs(title = "Local-LDA neighborhood", x = "Distance", y = "Local weight", colour = "Class", shape = NULL) + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}

# ---- algorithm-lvq ----

#' @rdname sb_knn_neighborhood
#' @export
sb_lvq_diagnostics <- function(fit, type = "prototypes") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!fit$pipeline$learner %in% sb_lvq_binary_learner_ids()) sb_abort(
    "This diagnostic requires an LVQ learner.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("prototypes", "prototype_distances")) sb_abort(
    "Choose prototypes or prototype_distances.", "stormbreaker_unsupported_diagnostic")
  state <- fit$learner_state
  e <- unserialize(state$engine_raw)
  z <- e$prototypes
  features <- state$features
  if (!is.matrix(z) || nrow(z) < 2L || ncol(z) != length(features) ||
      any(!is.finite(z)) || length(e$prototype_class) != nrow(z) ||
      !identical(sort(unique(e$prototype_class)), c("0", "1")) ||
      length(e$center) != ncol(z) || length(e$scale) != ncol(z) ||
      any(!is.finite(c(e$center, e$scale))) || any(e$scale <= 0)) sb_abort(
    "Saved LVQ codebook or training scaling is invalid.", "stormbreaker_bad_learner_state")
  labels <- paste0("P", seq_len(nrow(z)), " (", e$prototype_class, ")")
  original <- sweep(sweep(z, 2L, e$scale, "*"), 2L, e$center, "+")
  if (type == "prototypes") {
    tab <- data.frame(prototype = rep(labels, ncol(z)),
      prototype_index = rep(seq_len(nrow(z)), ncol(z)),
      class = rep(e$prototype_class, ncol(z)),
      feature = rep(features, each = nrow(z)), value = as.vector(z),
      original_value = as.vector(original))
    interpretation <- paste("Original learned prototypes in the fitted training-standardized coordinates;",
      "original-unit coordinates remain in the table. Parentheses identify encoded classes.",
      "Not individual patients, class means, feature importance or SHAP.")
  } else {
    distances <- as.matrix(stats::dist(z))^2
    tab <- data.frame(prototype = rep(labels, nrow(z)), prototype2 = rep(labels, each = nrow(z)),
      class = rep(e$prototype_class, nrow(z)), class2 = rep(e$prototype_class, each = nrow(z)),
      value = as.vector(distances))
    interpretation <- paste("Squared Euclidean distances between original learned prototypes,",
      "using the model's training scaling. Parentheses identify encoded classes.",
      "Not patient distances, calibrated risk, a class separation test or validation performance.")
  }
  structure(list(type = paste0("lvq_", type), table = tab, features = features, prototypes = labels,
    prototype_counts = data.frame(class = c("0", "1"), count = as.integer(table(factor(e$prototype_class, levels = c("0", "1"))))),
    center = stats::setNames(e$center, features), scale = stats::setNames(e$scale, features),
    iterations = e$iterations, learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    source = "original saved LVQ codebook and training scaling", interpretation = interpretation),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_lvq_diagnostic <- function(x) {
  feature <- prototype <- prototype2 <- value <- NULL
  tab <- x$table
  tab$prototype <- factor(tab$prototype, levels = rev(x$prototypes))
  if (x$type == "lvq_prototypes") {
    tab$feature <- factor(tab$feature, levels = x$features)
    p <- ggplot2::ggplot(tab, ggplot2::aes(feature, prototype, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#D77935",
        guide = ggplot2::guide_colourbar(display = "rectangles")) +
      ggplot2::labs(title = "Learned LVQ prototypes", x = NULL, y = "Prototype (class)",
        fill = "Standardized\ncoordinate")
    size <- c(width = max(11, 3 + .45 * length(x$features)), height = max(5, 2 + .35 * length(x$prototypes)))
  } else {
    tab$prototype2 <- factor(tab$prototype2, levels = x$prototypes)
    p <- ggplot2::ggplot(tab, ggplot2::aes(prototype2, prototype, fill = value)) +
      ggplot2::geom_tile(colour = "#D9DFE3", linewidth = 0.3) + ggplot2::coord_fixed() +
      ggplot2::scale_fill_gradient(low = "white", high = "#2878A0",
        guide = ggplot2::guide_colourbar(display = "rectangles")) +
      ggplot2::labs(title = "LVQ prototype distances", x = "Prototype (class)", y = "Prototype (class)",
        fill = "Squared\ndistance")
    size <- c(width = max(8, 3 + .4 * length(x$prototypes)), height = max(6, 2 + .4 * length(x$prototypes)))
  }
  p <- p + sb_viz_theme() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  attr(p, "stormbreaker_suggested_size") <- size
  p
}

# ---- loclda-neighbors ----

sb_loclda_neighbors <- function(params, p, n, classes) {
  k <- params$k_max
  if (identical(params$k_rule, "dimension_aware")) {
    k <- max(k, 2L * (p + classes))
  }
  min(n, k)
}
