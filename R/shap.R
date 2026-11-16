#' SHAP explanations
#' @name sb_model_shap
#' @description Compute supported saved-model SHAP explanations and plot their contributions.
#' @rawRd \arguments{
#' \item{fit}{
#' \strong{For \code{sb_model_shap()}:} Successful saved binary or survival final fit, including non-tree models.
#' 
#' \strong{For \code{sb_additive_shap()}:} A successful additive-model \code{stormbreaker_final_fit} returned by
#' \code{\link[=sb_fit_final]{sb_fit_final()}}.
#' 
#' \strong{For \code{sb_tree_shap()}:} Successful saved binary or survival tree-model final fit.
#' }
#' \item{new_data}{
#' \strong{For \code{sb_model_shap()}, \code{sb_tree_shap()}:} Observations to explain, including ID and predictor columns.
#' 
#' \strong{For \code{sb_additive_shap()}:} A non-empty data frame containing one row per explained
#' patient, the requested ID column, and every fitted numeric feature. Other
#' columns are ignored.
#' }
#' \item{id_col}{
#' \strong{For \code{sb_model_shap()}, \code{sb_tree_shap()}:} ID column in new_data.
#' 
#' \strong{For \code{sb_additive_shap()}:} Name of the unique patient identifier in \code{new_data}.
#' 
#' \strong{For \code{sb_plot_shap_importance()}:} ID-column name, required with feature_values.
#' 
#' \strong{For \code{sb_plot_shap_dependence()}:} ID column, matching explained observations exactly.
#' }
#' \item{background_data}{
#' \strong{For \code{sb_model_shap()}, \code{sb_tree_shap()}:} Development/reference predictors, not outcome columns.
#' 
#' \strong{For \code{sb_additive_shap()}:} A non-empty data frame containing every fitted numeric
#' feature. Rows define the empirical interventional reference distribution.
#' An ID and outcomes are neither required nor read; unrelated columns are
#' ignored.
#' }
#' \item{permutations}{
#' Number of random feature orders, default 32.
#' }
#' \item{background_n}{
#' Maximum background rows, default 20 (seeded sampling).
#' }
#' \item{seed}{
#' Seed for orders and background sampling.
#' }
#' \item{x}{
#' \strong{For \code{sb_plot_shap_importance()}:} A successful \code{\link[=sb_additive_shap]{sb_additive_shap()}} or \code{\link[=sb_tree_shap]{sb_tree_shap()}} or \code{\link[=sb_model_shap]{sb_model_shap()}} result.
#' The saved exact/approximate method and score scale are preserved.
#' 
#' \strong{For \code{sb_plot_shap_dependence()}:} A successful \code{\link[=sb_additive_shap]{sb_additive_shap()}} or \code{\link[=sb_tree_shap]{sb_tree_shap()}} result.
#' 
#' \strong{For \code{sb_plot_shap_heatmap()}:} A successful exact or approximate
#' saved SHAP result, as for importance.
#' }
#' \item{top_n}{
#' Maximum features shown, selected by mean absolute contribution
#' for importance/beeswarm/heatmap and by absolute individual contribution for the
#' waterfall. Ties retain feature input order.
#' }
#' \item{groupby}{
#' Optional group vector in explained-observation order. A named
#' vector is instead aligned by observation IDs. Groups must be nonmissing.
#' Uses the same saved model, SHAP values and background for every group;
#' top_n is ranked globally, then mean absolute SHAP is computed within groups.
#' }
#' \item{feature_values}{
#' \strong{For \code{sb_plot_shap_importance()}:} Optional data frame containing the explained feature
#' values and an explicit ID column. Rows may be reordered but IDs must match
#' the explained observations exactly; only fitted feature columns are read.
#' 
#' \strong{For \code{sb_plot_shap_dependence()}:} Feature data with explicit observation IDs.
#' }
#' \item{row}{
#' One-based explained-observation index for the waterfall.
#' }
#' \item{feature}{
#' One explained feature name.
#' }
#' \item{colour_by}{
#' Optional other explained feature used only for colouring.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_model_shap()}}{
#' Approximate SHAP for a saved binary or survival model
#' 
#' \strong{Returns:} A stormbreaker_model_shap result accepted by the SHAP plotting functions.
#' 
#' \strong{Details:} Uses the same permutation interventional estimator as sb_tree_shap,
#' with the original model's prediction interface. Supports non-tree models,
#' including LogitBoost and RDA. Binary values explain class-1 probability;
#' survival values explain the saved risk score. No model is refitted.
#' Background data should come from the development cohort. This is an
#' approximation, not exact TreeSHAP; more permutations reduce Monte Carlo noise.
#' }
#' 
#' \subsection{\code{sb_additive_shap()}}{
#' Also supports saved additive log-odds predictors from stepwise GLM, GLMBoost,
#' plsRglm, bias-reduced logistic and Cauchy-prior logistic regression, retaining
#' fitted centering and original parameter estimates.
#' LDA uses its exact original posterior log-odds equation, not a refitted GLM.
#' SLOPE retains its original selected penalized logistic coefficients.
#' Degree-one MARS uses its exact nonlinear additive hinge-basis contributions.
#' 
#' \strong{Returns:} A \code{stormbreaker_additive_shap} object containing the model and score
#' semantics, background feature and score-term means, scalar reference
#' value, a patient-by-feature contribution matrix, reconstructed additive
#' scores, and the unchanged state-only final prediction.
#' 
#' \strong{Details:} \code{sb_additive_shap()} explains a successful final fit whose saved prediction
#' score is additive over features. It currently supports ordinary, glmnet,
#' MCP/SCAD, and ABESS best-subset logistic regression, every saved linear
#' survival learner, and Gaussian Naive Bayes. It reads only the requested ID
#' and fitted numeric
#' features from
#' \code{new_data}; outcomes and unrelated columns are ignored. \code{background_data}
#' likewise supplies only a caller-chosen empirical feature reference and need
#' not contain an ID or outcome.
#' 
#' For an additive score \eqn{f(x) = c + \sum_j g_j(x_j)}, the returned value for
#' feature \eqn{j} is
#' \eqn{g_j(x_j) - mean_b(g_j(x_{bj}))}, where the mean is over the supplied
#' background rows. The reference value is
#' \eqn{c + \sum_j mean_b(g_j(x_{bj}))}. Their sum reconstructs the model score.
#' These are exact empirical interventional SHAP values for this additive score;
#' they are not a coefficient-times-input shortcut mislabeled as SHAP.
#' 
#' The score is class-1 log odds for logistic regression, class-1 log posterior
#' odds for Gaussian Naive Bayes, and the saved zero-reference \code{x_beta}
#' log-hazard ratio for a supported linear survival learner. Contributions are
#' not changes in
#' probability, survival probability, or absolute risk. The returned final
#' prediction keeps those task-specific outputs separate.
#' 
#' The background choice is part of the scientific question and changes both
#' the reference value and contributions. It should be declared and reported;
#' the function never defaults to the explained patients or reads their
#' outcomes. This interventional explanation does not preserve correlations or
#' estimate conditional SHAP values. It is descriptive, not causal, and does
#' not rank features or measure explanation stability.
#' 
#' Tree learners are deliberately unsupported here because their prediction
#' functions are not additive feature scores. The function does not silently
#' replace exact model-specific SHAP with a permutation approximation.
#' }
#' 
#' \subsection{\code{sb_tree_shap()}}{
#' Approximate interventional SHAP for saved tree-model predictions
#' 
#' \strong{Returns:} A stormbreaker_tree_shap result accepted by existing SHAP plots.
#' 
#' \strong{Details:} This is Monte Carlo permutation interventional SHAP, NOT exact
#' TreeSHAP. Every path telescopes from the same empirical background mean
#' to the explained prediction, so additivity alone does not prove estimator
#' accuracy. Binary explanations use class-1 probability; survival uses the
#' saved risk score, not survival probability or causal effects. Runtime grows
#' with features, patients, permutations and background size. No refitting.
#' }
#' 
#' \subsection{\code{sb_plot_shap_importance()}}{
#' Plot existing additive-SHAP explanations
#' 
#' \strong{Returns:} A ggplot object containing plotted contributions and descriptive
#' feature labels, without the original patient identifiers. These graphics
#' can still disclose individual predictions or feature values; review them
#' before sharing. Optional ggplot2 is required.
#' 
#' \strong{Details:} Plots preserve the saved score scale and empirical background
#' reference. Additive explanations are not probability effects; tree SHAP
#' can explicitly explain class-1 probability on its declared score scale.
#' Importance is mean absolute contribution over the explained observations,
#' not a causal importance estimate. The beeswarm uses deterministic binned
#' offsets and does not change the random state. Without supplied feature
#' values all points have one colour, rather than invented high/low values.
#' Optional values are aligned by explicit IDs and colour-normalized within
#' each feature; constant features map to the colour midpoint.
#' 
#' Waterfalls retain a baseline, the displayed feature terms, one sum of any
#' remaining terms, and the full reconstructed score. The displayed terms
#' therefore add up even when top_n omits features. Observation numbers, not
#' patient IDs, label individual plots. No model is fitted or explained again.
#' Importance bars display numeric mean absolute contributions. Waterfalls
#' use directional arrows (red increases, blue decreases), neutral endpoint
#' markers and a separate numeric column so narrow terms remain readable.
#' Beeswarm colour instead encodes within-feature values, not SHAP signs.
#' }
#' 
#' \subsection{\code{sb_plot_shap_dependence()}}{
#' Relate saved SHAP contributions to observed feature values
#' 
#' \strong{Returns:} A ggplot object containing all observations without original IDs.
#' Colouring by another variable does not compute SHAP interaction values.
#' No threshold, smoothing, causal effect or new explanation is inferred.
#' }
#'
#' \subsection{\code{sb_plot_shap_heatmap()}}{
#' Show signed saved SHAP contributions for each observation and selected feature.
#' The horizontal axis follows the saved observation order without printing IDs;
#' all supplied observations remain. No model fitting or new SHAP computation.
#' \strong{Returns:} A ggplot using the saved score scale.
#' }
#' }
#' @rawRd \references{
#' Lundberg SM, Lee S-I (2017). A unified approach to interpreting
#' model predictions. \emph{Advances in Neural Information Processing Systems 30},
#' 4765-4774. \url{https://arxiv.org/abs/1705.07874}.
#' }
#' @rawRd \examples{
#' # sb_model_shap
#' d <- sb_example_binary(60, seed = 42)
#' task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#' fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
#' s <- sb_model_shap(fit, d[1:2, ], "patient_id", d,
#'   permutations = 2, background_n = 3)
#' s$reference_value + rowSums(s$contributions)
#' 
#' # sb_additive_shap
#' development <- sb_example_binary(n = 40, seed = 1)
#' task <- sb_binary_task(
#'   development,
#'   "patient_id",
#'   "outcome",
#'   c("gene_a", "gene_b")
#' )
#' fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
#' sb_additive_shap(
#'   fit,
#'   development[1:3, ],
#'   "patient_id",
#'   background_data = development[4:20, ]
#' )
#' 
#' # sb_tree_shap
#' if (requireNamespace("rpart", quietly = TRUE)) {
#'   d <- sb_example_binary(60, seed = 42)
#'   task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#'   fit <- sb_fit_final(task, sb_pipeline(learner = "decision_tree_binary"))
#'   sb_tree_shap(fit, d[1:2, ], "patient_id", d,
#'     permutations = 2, background_n = 3)
#' }
#' }
NULL

# ---- additive-shap ----

sb_additive_shap_linear_terms <- function(state, x, learner) {
  if (learner %in% c("lda_binary", sb_extended_logistic_learner_ids())) {
    p <- if (learner == "lda_binary") sb_discriminant_parameters(state) else sb_extended_logistic_parameters(state)
    return(list(constant = p$intercept, terms = sweep(as.matrix(x), 2L, p$coefficient, "*"),
      score_scale = "class_1_log_odds"))
  }
  if (learner %in% c(
    "logistic", "glmnet_logistic", "slope_binary", sb_abess_binary_learner_ids(),
    sb_ncvreg_binary_learner_ids()
  )) {
    return(list(
      constant = state$intercept,
      terms = sweep(as.matrix(x), 2L, state$coefficient, "*"),
      score_scale = "class_1_log_odds"
    ))
  }
  list(
    constant = 0,
    terms = sweep(as.matrix(x), 2L, state$coefficient, "*"),
    score_scale = "zero_reference_log_hazard_ratio"
  )
}

sb_additive_shap_naive_bayes_terms <- function(state, x) {
  model <- structure(
    state$payload,
    class = "stormbreaker_gaussian_naive_bayes"
  )
  matrix_x <- as.matrix(x)
  class_terms <- function(class) {
    centered <- sweep(matrix_x, 2L, model$means[class, ], "-")
    scaled_distance <- sweep(
      centered^2,
      2L,
      model$variances[class, ],
      "/"
    )
    log_variance <- matrix(
      log(2 * pi * model$variances[class, ]),
      nrow = nrow(matrix_x),
      ncol = ncol(matrix_x),
      byrow = TRUE
    )
    -0.5 * (log_variance + scaled_distance)
  }
  terms <- class_terms("1") - class_terms("0")
  dimnames(terms) <- list(NULL, names(x))
  list(
    constant = log(model$priors[["1"]]) - log(model$priors[["0"]]),
    terms = terms,
    score_scale = "class_1_log_posterior_odds"
  )
}

sb_abort_unsupported_additive_shap <- function(learner) {
  sb_abort(
    paste0(
      "Exact additive SHAP is not implemented for learner `", learner,
      "`; STORMBREAKER will not substitute a model-agnostic approximation."
    ),
    "stormbreaker_unsupported_additive_shap"
  )
}

sb_additive_shap_terms <- function(state, x) {
  learner <- state$learner
  if (learner == "earth_binary") return(sb_additive_shap_earth_terms(state, x))
  if (learner %in% c(
    "logistic", "glmnet_logistic", "lda_binary", "slope_binary", sb_abess_binary_learner_ids(),
    sb_ncvreg_binary_learner_ids(),
    sb_linear_survival_learner_ids(), sb_extended_logistic_learner_ids()
  )) {
    return(sb_additive_shap_linear_terms(state, x, learner))
  }
  if (identical(learner, "gaussian_naive_bayes")) {
    return(sb_additive_shap_naive_bayes_terms(state, x))
  }
  sb_abort_unsupported_additive_shap(learner)
}

#' @rdname sb_model_shap
#' @export
sb_additive_shap <- function(fit, new_data, id_col, background_data) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  learner <- fit$learner_state$learner
  supported <- learner %in% c(
    "earth_binary",
    "logistic", "glmnet_logistic", "lda_binary", "slope_binary", sb_abess_binary_learner_ids(),
    sb_ncvreg_binary_learner_ids(),
    sb_linear_survival_learner_ids(), sb_extended_logistic_learner_ids(),
    "gaussian_naive_bayes"
  )
  if (!supported) {
    sb_abort_unsupported_additive_shap(learner)
  }

  patient_ids <- sb_extract_final_prediction_ids(new_data, id_col)
  x <- sb_select_final_prediction_features(
    new_data,
    id_col,
    fit$required_features
  )
  background_x <- sb_select_final_feature_frame(
    background_data,
    fit$required_features,
    "Background data",
    "background_data fitted features"
  )

  prediction <- sb_predict_final(fit, new_data, id_col)
  if (!identical(prediction$status, "success")) {
    sb_abort(
      paste0(
        "The final prediction needed for additive SHAP failed: ",
        prediction$error_message
      ),
      "stormbreaker_additive_shap_prediction_failure"
    )
  }

  new_components <- sb_additive_shap_terms(fit$learner_state, x)
  background_components <- sb_additive_shap_terms(
    fit$learner_state,
    background_x
  )
  reference_terms <- colMeans(background_components$terms)
  reference_value <- background_components$constant + sum(reference_terms)
  contributions <- sweep(
    new_components$terms,
    2L,
    reference_terms,
    "-"
  )
  finite <- is.finite(reference_value) &&
    all(is.finite(reference_terms)) && all(is.finite(contributions))
  if (!finite) {
    sb_abort(
      "The additive SHAP calculation produced non-finite score values.",
      "stormbreaker_nonfinite_additive_shap"
    )
  }
  rownames(contributions) <- patient_ids
  colnames(contributions) <- fit$required_features
  scores <- data.frame(
    patient_id = patient_ids,
    score = unname(reference_value + rowSums(contributions)),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      status = "success",
      task_type = fit$task_type,
      pipeline_id = fit$pipeline$id,
      learner = learner,
      method = "exact_interventional_shap_for_additive_score",
      score_scale = new_components$score_scale,
      background = list(
        n = as.integer(nrow(background_x)),
        feature_summary = data.frame(
          feature = fit$required_features,
          mean_feature_value = unname(colMeans(background_x)),
          mean_score_term = unname(reference_terms),
          stringsAsFactors = FALSE
        )
      ),
      reference_value = unname(reference_value),
      contributions = contributions,
      scores = scores,
      prediction = prediction
    ),
    class = "stormbreaker_additive_shap"
  )
}

#' @export
print.stormbreaker_additive_shap <- function(x, ...) {
  cat("<stormbreaker_additive_shap:", x$status, ">\n")
  cat("exact interventional SHAP for an additive model score\n")
  cat(
    "Patients:", nrow(x$contributions),
    "| Features:", ncol(x$contributions),
    "| Background rows:", x$background$n, "\n"
  )
  cat("Pipeline:", x$pipeline_id, "| Learner:", x$learner, "\n")
  cat("Score scale:", x$score_scale, "\n")
  cat("Contributions are not probability-scale effects.\n")
  cat("No feature is ranked and no model output is changed.\n")
  invisible(x)
}

# ---- tree-shap ----

sb_tree_shap_ids <- function() c("random_forest_binary", "gbm_binary",
  "gbm_cox", "gbm_cox_cv", "rsf_survival", "decision_tree_binary",
  "rpart_survival", "cforest_binary", "cforest_survival",
  sb_ranger_learner_ids(), sb_xgboost_learner_ids(),
  sb_lightgbm_binary_learner_ids(), "catboost_binary", "catboost_cox")

#' @rdname sb_model_shap
#' @export
sb_tree_shap <- function(fit, new_data, id_col, background_data,
                         permutations = 32L, background_n = 20L, seed = 1L) {
  fit <- sb_validate_final_fit_for_prediction(sb_upgrade_legacy_names(fit))
  if (!fit$pipeline$learner %in% sb_tree_shap_ids()) stop("This learner is not in the tree-SHAP menu.")
  result <- sb_model_shap_validated(fit, new_data, id_col, background_data,
    permutations, background_n, seed)
  class(result) <- "stormbreaker_tree_shap"
  result
}

#' @rdname sb_model_shap
#' @export
sb_model_shap <- function(fit, new_data, id_col, background_data,
                          permutations = 32L, background_n = 20L, seed = 1L) {
  fit <- sb_validate_final_fit_for_prediction(sb_upgrade_legacy_names(fit))
  sb_model_shap_validated(fit, new_data, id_col, background_data,
    permutations, background_n, seed)
}

sb_model_shap_validated <- function(fit, new_data, id_col, background_data,
                                    permutations, background_n, seed) {
  if (!fit$task_type %in% c("binary", "survival")) stop("Supply a binary or survival fit.")
  permutations <- sb_check_positive_integer(permutations, "permutations")
  background_n <- sb_check_positive_integer(background_n, "background_n")
  seed <- sb_check_positive_integer(seed, "seed")
  ids <- sb_extract_final_prediction_ids(new_data, id_col)
  features <- fit$required_features
  x <- sb_select_final_prediction_features(new_data, id_col, features)
  b <- sb_select_final_feature_frame(background_data, features,
    "Background data", "background_data")
  withr::local_seed(seed)
  if (nrow(b) > background_n) b <- b[sample.int(nrow(b), background_n), , drop = FALSE]
  score <- function(d) {
    values <- sb_predict_final_state(fit$learner_state, d)
    values <- sb_validate_final_prediction_values(fit, values, nrow(d))
    if (fit$task_type == "binary") values$probability else values$risk
  }
  reference <- mean(score(b))
  target <- score(x)
  if (length(target) != nrow(x) || any(!is.finite(target)) || !is.finite(reference)) stop("Saved prediction score is unavailable.")
  phi <- matrix(0, nrow(x), ncol(x), dimnames = list(ids, features))
  for (i in seq_len(nrow(x))) for (iteration in seq_len(permutations)) {
    order <- sample.int(ncol(x))
    hybrid <- b
    previous <- reference
    for (j in order) {
      hybrid[[j]] <- x[[j]][i]
      current <- mean(score(hybrid))
      phi[i, j] <- phi[i, j] + (current - previous) / permutations
      previous <- current
    }
  }
  structure(list(status = "success", task_type = fit$task_type,
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    method = "approximate_permutation_interventional_shap",
    score_scale = if (fit$task_type == "binary") "class_1_probability" else "saved_risk_score",
    reference_value = reference, contributions = phi,
    scores = data.frame(patient_id = ids, score = target),
    background = list(n = nrow(b)), permutations = permutations, seed = seed),
    class = "stormbreaker_model_shap")
}
