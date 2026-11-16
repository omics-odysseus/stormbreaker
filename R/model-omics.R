#' Projection and discriminant diagnostics
#' @name sb_pls_diagnostics
#' @description Inspect projection components, discriminant geometry, gene-pair rules and class centroids.
#' @rawRd \arguments{
#' \item{fit}{
#' \strong{For \code{sb_pls_diagnostics()}:} A final plsRglm, PLS-DA, sparse PLS-DA or plsRcox model fitted with
#' capture_diagnostics = TRUE for components. CV curves use saved tuning and
#' do not require capture. Historical compact binary fits lack components.
#' 
#' \strong{For \code{sb_superpc_diagnostics()}:} A final SuperPC model fitted with capture_diagnostics = TRUE.
#' 
#' \strong{For \code{sb_cancerclass_diagnostics()}:} Successful Cancerclass final fit.
#' 
#' \strong{For \code{sb_high_dimensional_diagnostics()}:} Successful HDA or HDDA final fit.
#' 
#' \strong{For \code{sb_ktsp_diagnostics()}:} A successful k-TSP binary final fit.
#' 
#' \strong{For \code{sb_pamr_diagnostics()}:} A successful binary or survival PAMR final fit.
#' 
#' \strong{For \code{sb_fda_diagnostics()}:} Successful FDA, PDA or MDA final fit.
#' 
#' \strong{For \code{sb_discriminant_diagnostics()}:} A successful LDA, QDA or fixed-regularization klaR RDA final fit.
#' 
#' \strong{For \code{sb_shrinkage_discriminant_diagnostics()}:} Successful SDA or SLDA final fit.
#' }
#' \item{type}{
#' \strong{For \code{sb_pls_diagnostics()}:} loadings, score_distribution or cv_curve (CV plsRcox only).
#' 
#' \strong{For \code{sb_superpc_diagnostics()}:} feature_scores, loadings or cv_curve (CV variant only).
#' 
#' \strong{For \code{sb_cancerclass_diagnostics()}:} centroids or score_mapping.
#' 
#' \strong{For \code{sb_high_dimensional_diagnostics()}:} HDA: loadings or densities. HDDA: loadings, spectrum or centers.
#' NULL chooses loadings.
#' 
#' \strong{For \code{sb_fda_diagnostics()}:} centers, query_variates, or mixture_weights (MDA only).
#' 
#' \strong{For \code{sb_discriminant_diagnostics()}:} class_means or class_correlation (QDA/RDA only).
#' 
#' \strong{For \code{sb_shrinkage_discriminant_diagnostics()}:} coefficients or regularization for SDA; loadings or
#' class_means for SLDA. NULL chooses the first applicable view.
#' }
#' \item{score_range}{
#' Finite increasing native-score plotting range, default
#' c(-1, 1). This is a display range, not an observed training-score range.
#' }
#' \item{newdata}{
#' Explicit query rows for query_variates.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_pls_diagnostics()}}{
#' Inspect native PLS component structure
#' 
#' \strong{Returns:} An algorithm diagnostic with all original components and features.
#' mixOmics returns its native X component weights; plsRglm/plsRcox return X
#' reconstruction loadings (pp). These definitions are not interchangeable.
#' Scores describe the development samples, not held-out discrimination.
#' Score tables omit patient identifiers but remain sensitive row-level data.
#' No new decomposition, rotation, sign alignment or fitting is performed.
#' CV criterion 5 is integrated AUC (iAUC_sh), higher is better, despite the
#' historical field name cv_error. Native cv.se5 is retained as native spread,
#' not called a conventional standard error or confidence interval. The native
#' implementation uses fold standard deviation divided by the fold count.
#' }
#' 
#' \subsection{\code{sb_superpc_diagnostics()}}{
#' Inspect original supervised principal-component diagnostics
#' 
#' \strong{Returns:} An algorithm diagnostic retaining native feature scores, selected
#' singular-vector loadings, or the full component-by-threshold CV table.
#' Feature scores are training associations, not model coefficients. Loadings
#' retain the native SVD sign and are not hazard ratios or variable importance.
#' CV scores are geometric means of held-out likelihood-ratio statistics;
#' bounds use plus/minus one standard error on the log scale, not confidence
#' intervals. The original pipeline selects thresholds using one component,
#' even though native CV evaluates more components. No refitting is performed.
#' Compare thresholds within a component count; likelihood-ratio statistics
#' with different degrees of freedom are not a model ranking across counts.
#' }
#' 
#' \subsection{\code{sb_cancerclass_diagnostics()}}{
#' Inspect original Cancerclass centroids and probability mapping
#' 
#' \strong{Returns:} Original centroids or the saved inner-OOF-fitted logistic mapping.
#' The mapping is not a calibration-reliability curve or a validation result.
#' Original OOF scores are not present in compact fits and are not recreated.
#' No new calibration fit or surrogate model is introduced.
#' }
#' 
#' \subsection{\code{sb_high_dimensional_diagnostics()}}{
#' Inspect original HDA and HDDA discriminant subspaces
#' 
#' \strong{Returns:} Original retained subspace loadings, HDA Gaussian coordinate
#' densities, or HDDA eigenvalues, retained dimensions and original-unit
#' centers. HDDA loadings and eigenvalues use the fitted standardized space.
#' Spectra are not a newly computed PCA or explained-variance analysis.
#' Axis signs are basis conventions, not feature importance or effects.
#' }
#' 
#' \subsection{\code{sb_ktsp_diagnostics()}}{
#' Inspect the original k-top-scoring-pairs decision rules
#' 
#' \strong{Returns:} An algorithm diagnostic preserving all selected ordered feature pairs
#' and their original training scores. A first-greater-than-second comparison
#' contributes one class-one vote, equality contributes one half and a smaller
#' first value contributes zero. The package averages pair votes; this fraction
#' is not independently calibrated probability. Native pair scores are training
#' selection evidence, not validation accuracy, uncertainty or vote weights.
#' Compact saved state suffices; no pairs are retrained or reordered.
#' }
#' 
#' \subsection{\code{sb_pamr_diagnostics()}}{
#' Inspect original and shrunken PAMR centroids
#' 
#' \strong{Returns:} All feature/class centroids before and after the model's original
#' fixed threshold. Raw input-unit centroids and deviations from the overall
#' centroid divided by the saved PAMR SD are retained. This is the native
#' PAMR scaling, not a new feature-wise z-score. Survival classes are native
#' censoring-aware early/late groups, not external patient subgroups or
#' horizon-specific survival probabilities. No threshold tuning is performed.
#' }
#' 
#' \subsection{\code{sb_fda_diagnostics()}}{
#' Inspect original flexible, penalized and mixture discriminant geometry
#' 
#' \strong{Returns:} Original class/subclass centers, within-class mixture weights,
#' or native query discriminant coordinates. Axes are discriminant variates,
#' not original features. Mixture subclasses are fitted latent components,
#' not clinical subgroups. No refitting, new projection or covariance estimate.
#' }
#' 
#' \subsection{\code{sb_discriminant_diagnostics()}}{
#' Inspect original LDA, QDA and regularized discriminant parameters
#' 
#' \strong{Returns:} An algorithm diagnostic from compact saved model parameters, without
#' capture, refitting, covariance estimation or surrogate models. Class means
#' use original feature units; they are not confidence intervals. Correlations
#' describe the fitted within-class covariance used for prediction. QDA uses
#' its native precision factors; RDA applies both native regularizations.
#' LDA's rank-reduced discriminant scaling cannot identify the full covariance,
#' so this view is explicitly unavailable for LDA. Its original linear log-odds
#' equation is available through \code{\link[=sb_model_coefficients]{sb_model_coefficients()}} and \code{\link[=sb_nomogram]{sb_nomogram()}}.
#' The result retains the exact equation as \code{parameters}: intercept plus
#' x times coefficient plus x' quadratic x. Off-diagonal matrix entries occur
#' twice in that quadratic form. QDA/RDA are not treated as additive models.
#' Source and interpretation belong in a separate figure legend, not panels.
#' }
#' 
#' \subsection{\code{sb_shrinkage_discriminant_diagnostics()}}{
#' Inspect saved shrinkage and subspace discriminant models
#' 
#' \strong{Returns:} Native shrinkage parameters, class-1 log-score contrast, subspace
#' loadings, or class means in that subspace. SDA's native posterior method
#' rounds probabilities; coefficients describe its pre-rounding log odds.
#' SLDA axes are fitted subspace coordinates, not original features.
#' No covariance reconstruction, feature selection or refitting occurs.
#' }
#' }
NULL

# ---- omics-specialist-binary-adapters ----

sb_omics_specialist_binary_learner_ids <- function() {
  c(
    "brglm2_logistic",
    "bayesglm_logistic_cauchy",
    "bst_hinge_binary",
    "hda_binary",
    "dwd_binary",
    "ktsp_binary"
  )
}

sb_omics_specialist_binary_params <- function(learner) {
  switch(learner,
    brglm2_logistic = list(
      type = "AS_mean",
      link = "logit"
    ),
    bayesglm_logistic_cauchy = list(
      link = "logit",
      prior_mean = 0,
      prior_scale = 2.5,
      prior_df = 1,
      intercept_prior_mean = 0,
      intercept_prior_scale = 10,
      intercept_prior_df = 1,
      scaled = TRUE,
      maximum_iterations = 100L
    ),
    bst_hinge_binary = list(
      family = "hinge",
      cost = 0.5,
      base_learner = "ls",
      iterations = 50L,
      learning_rate = 0.1,
      center = TRUE,
      fit_intercept = FALSE,
      probability_mapping = "binary_softmax_native_margin"
    ),
    hda_binary = list(
      reduced_dimension = 1L,
      classification_rule = "naive_bayes",
      regularization_lambda = NULL,
      regularization_gamma = NULL,
      initial_loadings = "identity",
      outer_iterations = 7L,
      inner_iterations = 10L
    ),
    dwd_binary = list(
      exponent = 1,
      penalty_rule = "DWDLargeR_penaltyParameter",
      remove_zero_norm_features = TRUE,
      scale_features = TRUE,
      tolerance = 1e-5,
      maximum_iterations = 2000L,
      probability_mapping = "binary_softmax_native_margin"
    ),
    ktsp_binary = list(
      maximum_pairs = 10L,
      feature_filter = "wilcoxon",
      feature_cap = 100L,
      disjoint_pairs = TRUE,
      handle_ties = FALSE,
      probability_mapping = "positive_pair_vote_fraction"
    ),
    sb_abort(
      paste0("Unknown omics-specialist binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_omics_specialist_binary_dependency <- function(learner) {
  package <- switch(learner,
    brglm2_logistic = "brglm2",
    bayesglm_logistic_cauchy = "arm",
    bst_hinge_binary = "bst",
    hda_binary = "hda",
    dwd_binary = "DWDLargeR",
    ktsp_binary = "switchBox",
    sb_abort(
      paste0("Unknown omics-specialist binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
  list(
    package = package,
    packages = package,
    available = requireNamespace(package, quietly = TRUE)
  )
}

sb_require_omics_specialist_binary_dependency <- function(learner) {
  dependency <- sb_omics_specialist_binary_dependency(learner)
  if (!dependency$available) {
    sb_abort(
      paste0(
        "Package `", dependency$package, "` is required for learner `",
        learner, "`."
      ),
      "stormbreaker_missing_omics_specialist_binary_dependency"
    )
  }
  invisible(dependency)
}

sb_omics_specialist_binary_parameter_semantics <- function(learner) {
  switch(learner,
    brglm2_logistic = paste0(
      "mean bias-reduced logistic adjusted-score fit with logit link fixed ",
      "before resampling"
    ),
    bayesglm_logistic_cauchy = paste0(
      "Bayesian logistic posterior-mode fit with scaled Cauchy coefficient ",
      "priors and a Cauchy intercept prior fixed before resampling"
    ),
    bst_hinge_binary = paste0(
      "50-step componentwise linear hinge-loss boosting with shrinkage 0.1, ",
      "centering, and a fixed softmax mapping of the native margin"
    ),
    hda_binary = paste0(
      "unregularized one-dimensional heteroscedastic discriminant ",
      "transformation with identity initialization, nonsingular class ",
      "covariances, and a native Gaussian naive-Bayes rule"
    ),
    dwd_binary = paste0(
      "distance-weighted discrimination exponent 1 with native ",
      "analysis-fold penalty rule, scaling, and all-zero feature removal"
    ),
    ktsp_binary = paste0(
      "disjoint K-TSP with analysis-fold Wilcoxon filtering, at most ten ",
      "pairs, and native K selection"
    )
  )
}

sb_new_omics_specialist_binary_model <- function(
  learner, features, params, engine
) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      engine = engine
    ),
    class = "stormbreaker_omics_specialist_binary_model"
  )
}

sb_fit_brglm2_logistic <- function(x, outcome, params) {
  design <- cbind("(Intercept)" = 1, data.matrix(x))
  fit <- brglm2::brglmFit(
    x = design,
    y = outcome,
    family = stats::binomial(link = params$link),
    control = brglm2::brglmControl(type = params$type),
    intercept = TRUE,
    singular.ok = FALSE
  )
  coefficients <- unname(fit$coefficients)
  valid <- isTRUE(fit$converged) && identical(fit$type, params$type) &&
    identical(as.integer(fit$rank), as.integer(ncol(design))) &&
    length(coefficients) == ncol(design) &&
    !anyNA(coefficients) && all(is.finite(coefficients))
  if (!valid) {
    sb_abort(
      paste0(
        "Mean bias-reduced logistic regression did not return a converged, ",
        "full-rank finite fit."
      ),
      "stormbreaker_incomplete_brglm2_logistic_fit"
    )
  }
  list(
    intercept = coefficients[[1L]],
    coefficient = stats::setNames(coefficients[-1L], names(x)),
    iterations = as.integer(fit$iter),
    type = fit$type
  )
}

sb_fit_bayesglm_logistic_cauchy <- function(x, outcome, params) {
  design <- cbind("(Intercept)" = 1, data.matrix(x))
  fit <- suppressWarnings(arm::bayesglm.fit(
    x = design,
    y = outcome,
    family = stats::binomial(link = params$link),
    control = list(maxit = params$maximum_iterations),
    intercept = TRUE,
    prior.mean = params$prior_mean,
    prior.scale = params$prior_scale,
    prior.df = params$prior_df,
    prior.mean.for.intercept = params$intercept_prior_mean,
    prior.scale.for.intercept = params$intercept_prior_scale,
    prior.df.for.intercept = params$intercept_prior_df,
    scaled = params$scaled,
    Warning = FALSE
  ))
  coefficients <- unname(fit$coefficients)
  valid <- isTRUE(fit$converged) &&
    identical(as.integer(fit$rank), as.integer(ncol(design))) &&
    length(coefficients) == ncol(design) &&
    !anyNA(coefficients) && all(is.finite(coefficients))
  if (!valid) {
    sb_abort(
      "Bayesian Cauchy logistic regression did not return a finite fit.",
      "stormbreaker_incomplete_bayesglm_logistic_fit"
    )
  }
  list(
    intercept = coefficients[[1L]],
    coefficient = stats::setNames(coefficients[-1L], names(x)),
    iterations = as.integer(fit$iter)
  )
}

sb_fit_bst_hinge_binary <- function(x, outcome, params) {
  matrix_x <- data.matrix(x)
  labels <- ifelse(outcome == 1L, 1, -1)
  fit <- bst::bst(
    matrix_x,
    labels,
    cost = params$cost,
    family = params$family,
    learner = params$base_learner,
    ctrl = bst::bst_control(
      mstop = params$iterations,
      nu = params$learning_rate,
      center = params$center,
      intercept = params$fit_intercept,
      trace = FALSE
    )
  )
  coefficient <- numeric(ncol(matrix_x))
  for (iteration in seq_len(params$iterations)) {
    feature <- fit$ensemble[[iteration]]
    coefficient[[feature]] <- coefficient[[feature]] +
      params$learning_rate * fit$coef[[iteration]]
  }
  intercept <- fit$offset - sum(coefficient * fit$meanx)
  valid <- inherits(fit, "bst") && identical(fit$family, params$family) &&
    identical(fit$learner, params$base_learner) &&
    length(coefficient) == ncol(matrix_x) &&
    !anyNA(coefficient) && all(is.finite(coefficient)) &&
    length(intercept) == 1L && is.finite(intercept)
  if (!valid) {
    sb_abort(
      "Hinge-loss boosting did not return a finite linear margin.",
      "stormbreaker_incomplete_bst_hinge_fit"
    )
  }
  list(
    intercept = unname(intercept),
    coefficient = stats::setNames(coefficient, names(x)),
    selected_features = names(x)[coefficient != 0],
    final_training_risk = unname(fit$risk[[params$iterations]]),
    diagnostic_path = list(risk = fit$risk, selected = fit$ensemble,
      updates = fit$coef, meanx = fit$meanx, offset = fit$offset,
      learning_rate = params$learning_rate, family = fit$family,
      learner = fit$learner, intercept_updates = fit$int)
  )
}

sb_hda_binary_prediction_engine <- function(engine) {
  keep <- c("hda.loadings", "reduced.dimension", "naivebayes")
  lean <- engine[keep]
  class(lean) <- class(engine)
  lean
}

sb_predict_hda_binary <- function(engine, matrix_x) {
  prediction <- suppressWarnings(stats::predict(
    engine,
    newdata = matrix_x,
    task = "c"
  ))
  positive <- match("1", engine$naivebayes$levels)
  posterior <- prediction$posteriors
  if (is.na(positive) || !is.matrix(posterior) ||
    nrow(posterior) != nrow(matrix_x) || ncol(posterior) < positive) {
    sb_abort(
      "HDA did not return a class-one posterior.",
      "stormbreaker_invalid_hda_binary_prediction"
    )
  }
  unname(posterior[, positive])
}

sb_fit_hda_binary <- function(x, outcome, params) {
  engine <- NULL
  invisible(utils::capture.output(
    engine <- hda::hda(
      data.matrix(x),
      grouping = factor(outcome, levels = c(0L, 1L)),
      newdim = params$reduced_dimension,
      crule = TRUE,
      reg.lamb = params$regularization_lambda,
      reg.gamm = params$regularization_gamma,
      initial.loadings = NULL,
      noutit = params$outer_iterations,
      ninit = params$inner_iterations,
      verbose = FALSE
    ),
    type = "output"
  ))
  engine <- sb_hda_binary_prediction_engine(engine)
  valid <- inherits(engine, "hda") &&
    is.matrix(engine$hda.loadings) &&
    nrow(engine$hda.loadings) == ncol(x) &&
    identical(engine$reduced.dimension, params$reduced_dimension) &&
    inherits(engine$naivebayes, "naiveBayes") &&
    identical(engine$naivebayes$levels, c("0", "1"))
  if (!valid) {
    sb_abort(
      "HDA did not return a two-class discriminant prediction engine.",
      "stormbreaker_incomplete_hda_binary_fit"
    )
  }
  sb_predict_hda_binary(engine, data.matrix(x))
  engine
}

sb_dwd_native_feature_map <- function(matrix_x, params) {
  feature_norm <- sqrt(colSums(matrix_x * matrix_x))
  kept <- seq_len(ncol(matrix_x))
  if (isTRUE(params$remove_zero_norm_features) && any(feature_norm == 0)) {
    kept <- which(feature_norm > 0)
    feature_norm <- c(feature_norm[kept], 0)
  }

  multiplier <- rep(1, length(feature_norm))
  dimension <- length(feature_norm)
  n <- nrow(matrix_x)
  if (isTRUE(params$scale_features) && dimension > 0.5 * n &&
    max(feature_norm) > 2 * min(feature_norm)) {
    denominator <- if (dimension > 3 * n) {
      sqrt(feature_norm)
    } else {
      feature_norm
    }
    multiplier <- 1 / pmax(1, denominator)
  }
  list(
    kept = kept,
    multiplier = multiplier
  )
}

sb_fit_dwd_binary <- function(x, outcome, params) {
  matrix_x <- data.matrix(x)
  feature_map <- sb_dwd_native_feature_map(matrix_x, params)
  labels <- ifelse(outcome == 1L, 1, -1)
  sparse_x <- methods::as(t(matrix_x), "dgCMatrix")
  fit <- NULL
  penalty <- NULL
  withr::with_preserve_seed(
    invisible(utils::capture.output(
      {
        penalty <- DWDLargeR::penaltyParameter(
          sparse_x,
          labels,
          expon = params$exponent
        )
        fit <- DWDLargeR::genDWD(
          sparse_x,
          labels,
          C = penalty,
          expon = params$exponent,
          tol = params$tolerance,
          maxIter = params$maximum_iterations,
          printDetails = 0,
          rmzeroFea = as.integer(params$remove_zero_norm_features),
          scaleFea = as.integer(params$scale_features)
        )
      },
      type = "output"
    ))
  )
  native_coefficient <- if (is.list(fit)) as.numeric(fit$w) else numeric()
  intercept <- if (is.list(fit)) as.numeric(fit$beta) else numeric()
  valid <- is.numeric(penalty) && length(penalty) == 1L &&
    !is.na(penalty) && is.finite(penalty) && penalty > 0 &&
    length(intercept) == 1L && !is.na(intercept) && is.finite(intercept) &&
    length(native_coefficient) == length(feature_map$multiplier) &&
    !anyNA(native_coefficient) && all(is.finite(native_coefficient))
  if (!valid) {
    sb_abort(
      "Distance-weighted discrimination did not return a finite linear fit.",
      "stormbreaker_incomplete_dwd_binary_fit"
    )
  }
  coefficient <- numeric(ncol(matrix_x))
  if (length(feature_map$kept) > 0L) {
    kept_index <- seq_along(feature_map$kept)
    coefficient[feature_map$kept] <-
      native_coefficient[kept_index] * feature_map$multiplier[kept_index]
  }
  list(
    intercept = unname(intercept),
    coefficient = stats::setNames(coefficient, names(x)),
    penalty = unname(as.numeric(penalty)),
    iterations = if (is.list(fit$info)) as.integer(fit$info$iter) else NA_integer_,
    training_error = if (is.list(fit$info)) {
      as.numeric(fit$info$trainerr)
    } else {
      NA_real_
    }
  )
}

sb_fit_ktsp_binary <- function(x, outcome, params) {
  if (ncol(x) < 4L) {
    sb_abort(
      "K-TSP requires at least four features in the current analysis fold.",
      "stormbreaker_ktsp_too_few_features"
    )
  }
  maximum_pairs <- min(params$maximum_pairs, floor(ncol(x) / 2L))
  classifier <- switchBox::SWAP.Train.KTSP(
    inputMat = t(data.matrix(x)),
    phenoGroup = factor(outcome, levels = c(0L, 1L)),
    krange = seq_len(maximum_pairs),
    FilterFunc = switchBox::SWAP.Filter.Wilcoxon,
    handleTies = params$handle_ties,
    disjoint = params$disjoint_pairs,
    featureNo = min(params$feature_cap, ncol(x)),
    verbose = FALSE
  )
  pairs <- classifier$TSPs
  valid <- is.matrix(pairs) && ncol(pairs) == 2L && nrow(pairs) >= 1L &&
    nrow(pairs) <= maximum_pairs && all(pairs %in% names(x)) &&
    identical(as.character(classifier$labels), c("0", "1")) &&
    length(classifier$tieVote) == nrow(pairs) &&
    all(as.character(classifier$tieVote) == "both")
  if (!valid) {
    sb_abort(
      "K-TSP did not return an oriented set of class-one feature pairs.",
      "stormbreaker_incomplete_ktsp_binary_fit"
    )
  }
  list(
    pairs = unname(pairs),
    pair_scores = unname(as.numeric(classifier$score)),
    selected_pairs = as.integer(nrow(pairs)),
    labels = c("0", "1")
  )
}

sb_fit_omics_specialist_binary <- function(x, outcome, learner, params) {
  inputs <- sb_validate_classical_binary_inputs(x, outcome, learner)
  sb_require_omics_specialist_binary_dependency(learner)
  engine <- switch(learner,
    brglm2_logistic = sb_fit_brglm2_logistic(
      inputs$x,
      inputs$outcome,
      params
    ),
    bayesglm_logistic_cauchy =
      sb_fit_bayesglm_logistic_cauchy(
        inputs$x,
        inputs$outcome,
        params
      ),
    bst_hinge_binary = sb_fit_bst_hinge_binary(
      inputs$x,
      inputs$outcome,
      params
    ),
    hda_binary = sb_fit_hda_binary(
      inputs$x,
      inputs$outcome,
      params
    ),
    dwd_binary = sb_fit_dwd_binary(
      inputs$x,
      inputs$outcome,
      params
    ),
    ktsp_binary = sb_fit_ktsp_binary(
      inputs$x,
      inputs$outcome,
      params
    )
  )
  diagnostic_path <- NULL
  if (learner == "bst_hinge_binary") {
    diagnostic_path <- engine$diagnostic_path
    engine$diagnostic_path <- NULL
  }
  model <- sb_new_omics_specialist_binary_model(
    learner,
    names(inputs$x),
    params,
    engine
  )
  if (!is.null(diagnostic_path)) model$diagnostic_path <- diagnostic_path
  model
}

sb_predict_omics_specialist_binary <- function(model, x) {
  if (!inherits(model, "stormbreaker_omics_specialist_binary_model")) {
    sb_abort(
      "The supplied model is not a fitted omics-specialist binary learner.",
      "stormbreaker_bad_omics_specialist_binary_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      paste0(
        "Prediction columns do not match the omics-specialist binary ",
        "learner schema."
      ),
      "stormbreaker_learner_schema_mismatch"
    )
  }
  matrix_x <- data.matrix(x)
  probability <- switch(model$learner,
    brglm2_logistic = stats::plogis(
      model$engine$intercept +
        as.vector(matrix_x %*% model$engine$coefficient)
    ),
    bayesglm_logistic_cauchy = stats::plogis(
      model$engine$intercept +
        as.vector(matrix_x %*% model$engine$coefficient)
    ),
    bst_hinge_binary = {
      margin <- model$engine$intercept +
        as.vector(matrix_x %*% model$engine$coefficient)
      stats::plogis(2 * margin)
    },
    hda_binary = sb_predict_hda_binary(
      model$engine,
      matrix_x
    ),
    dwd_binary = {
      margin <- model$engine$intercept +
        as.vector(matrix_x %*% model$engine$coefficient)
      stats::plogis(2 * margin)
    },
    ktsp_binary = {
      first <- match(model$engine$pairs[, 1L], model$features)
      second <- match(model$engine$pairs[, 2L], model$features)
      first_value <- matrix_x[, first, drop = FALSE]
      second_value <- matrix_x[, second, drop = FALSE]
      rowMeans((first_value > second_value) +
        0.5 * (first_value == second_value))
    }
  )
  sb_class_one_probability(probability, nrow(x), model$learner)
}

sb_extract_omics_specialist_binary_coefficients <- function(model) {
  if (!model$learner %in% c(
    "brglm2_logistic",
    "bayesglm_logistic_cauchy",
    "bst_hinge_binary",
    "dwd_binary"
  )) {
    sb_abort(
      "This omics-specialist learner has no linear coefficient vector.",
      "stormbreaker_omics_specialist_coefficients_not_applicable"
    )
  }
  data.frame(
    term = c("(Intercept)", model$features),
    coefficient = unname(c(
      model$engine$intercept,
      model$engine$coefficient
    )),
    stringsAsFactors = FALSE
  )
}

sb_extract_omics_specialist_binary_state <- function(model) {
  if (!inherits(model, "stormbreaker_omics_specialist_binary_model")) {
    sb_bad_learner_state(
      "The fitted omics-specialist binary learner state is invalid."
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      engine = model$engine,
      training_data_retention = "none"
    ),
    class = c(
      "stormbreaker_omics_specialist_binary_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_omics_specialist_binary_state <- function(state, x) {
  valid <- inherits(state, "stormbreaker_omics_specialist_binary_state") &&
    is.list(state) &&
    state$learner %in% sb_omics_specialist_binary_learner_ids() &&
    is.list(state$params) && is.list(state$engine) &&
    identical(state$training_data_retention, "none")
  if (!isTRUE(valid)) {
    sb_bad_learner_state(
      "The saved omics-specialist binary learner state is invalid."
    )
  }
  x <- sb_validate_learner_state_x(state, x, state$learner)
  model <- sb_new_omics_specialist_binary_model(
    state$learner,
    state$features,
    state$params,
    state$engine
  )
  list(probability = sb_predict_omics_specialist_binary(model, x))
}

sb_omics_specialist_binary_learner_adapter <- function(
  learner,
  params = sb_omics_specialist_binary_params(learner)
) {
  adapter <- list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_omics_specialist_binary(x, outcome, learner, params)
    },
    predict = sb_predict_omics_specialist_binary,
    extract_state = function(model, context = list()) {
      sb_extract_omics_specialist_binary_state(model)
    },
    predict_state = sb_predict_omics_specialist_binary_state
  )
  if (learner %in% c(
    "brglm2_logistic",
    "bayesglm_logistic_cauchy",
    "bst_hinge_binary",
    "dwd_binary"
  )) {
    adapter$extract_coefficients <-
      sb_extract_omics_specialist_binary_coefficients
  }
  adapter
}

# ---- high-dimensional-omics-adapters ----

sb_high_dimensional_binary_selector_ids <- function() {
  c(
    "sis_binary_nlogn",
    "vsurf_binary_prediction",
    "mrmr_binary_nlogn",
    sb_information_binary_selector_ids(),
    "lpc_binary_nlogn",
    "stabsel_binary_lasso"
  )
}

sb_information_binary_selector_ids <- function() {
  c(
    "jmi_binary_nlogn",
    "cmim_binary_nlogn",
    "disr_binary_nlogn",
    "jmim_binary_nlogn"
  )
}

sb_high_dimensional_survival_selector_ids <- function() {
  c("lpc_cox_nlogn", "stabsel_cox_lasso")
}

sb_lpc_selector_ids <- function() {
  c("lpc_binary_nlogn", "lpc_cox_nlogn")
}

sb_high_dimensional_selector_ids <- function() {
  c(
    sb_high_dimensional_binary_selector_ids(),
    sb_high_dimensional_survival_selector_ids()
  )
}

sb_hdda_binary_learner_ids <- function() {
  "hdda_binary_cattell"
}

sb_high_dimensional_method_ids <- function() {
  c(sb_high_dimensional_selector_ids(), sb_hdda_binary_learner_ids())
}

sb_high_dimensional_method_packages <- function(method) {
  switch(method,
    sis_binary_nlogn = "SIS",
    vsurf_binary_prediction = "VSURF",
    mrmr_binary_nlogn = "praznik",
    jmi_binary_nlogn = "praznik",
    cmim_binary_nlogn = "praznik",
    disr_binary_nlogn = "praznik",
    jmim_binary_nlogn = "praznik",
    lpc_binary_nlogn = "lpc",
    lpc_cox_nlogn = "lpc",
    stabsel_binary_lasso = c("stabs", "glmnet"),
    stabsel_cox_lasso = c("stabs", "glmnet"),
    hdda_binary_cattell = "HDclassif",
    sb_abort(
      paste0("Unknown high-dimensional method `", method, "`."),
      "stormbreaker_unsupported_high_dimensional_method"
    )
  )
}

sb_require_high_dimensional_method <- function(method) {
  packages <- sb_high_dimensional_method_packages(method)
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0L) {
    sb_abort(
      paste0(
        "Package", if (length(missing) == 1L) " `" else "s `",
        paste(missing, collapse = "`, `"), "` ",
        if (length(missing) == 1L) "is" else "are",
        " required for method `", method, "`."
      ),
      "stormbreaker_missing_high_dimensional_dependency"
    )
  }
  invisible(packages)
}

sb_high_dimensional_method_provenance <- function(method) {
  packages <- sb_high_dimensional_method_packages(method)
  fields <- c(
    SIS = "sis_version",
    VSURF = "vsurf_version",
    praznik = "praznik_version",
    lpc = "lpc_version",
    stabs = "stabs_version",
    glmnet = "glmnet_version",
    HDclassif = "hdclassif_version"
  )
  versions <- lapply(packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
      as.character(utils::packageVersion(package))
    } else {
      NA_character_
    }
  })
  stats::setNames(versions, unname(fields[packages]))
}

sb_high_dimensional_native_call <- function(code, method, stage) {
  tryCatch(
    withCallingHandlers(
      force(code),
      warning = function(warning) {
        sb_abort(
          paste0(
            "Method `", method, "` raised a ", stage, " warning: ",
            conditionMessage(warning)
          ),
          "stormbreaker_high_dimensional_native_warning",
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
          "Method `", method, "` failed during ", stage, ": ",
          conditionMessage(error)
        ),
        "stormbreaker_high_dimensional_native_error",
        parent = error
      )
    }
  )
}

sb_validate_high_dimensional_binary_inputs <- function(x, outcome, method) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      paste0("A `", method, "` analysis fold must contain both classes."),
      "stormbreaker_high_dimensional_single_class"
    )
  }
  list(x = x, outcome = outcome)
}

sb_sis_binary_selector_params <- function() {
  list(
    family = "binomial",
    penalty = "lasso",
    tuning = "bic",
    selection_output = "sis.ix0_before_penalized_fit",
    screen_size_rule = "min(p, max(2, floor(n/log(n))))",
    iterative = FALSE,
    standardize = TRUE,
    parallel = FALSE,
    seed = 1L
  )
}

sb_high_dimensional_nlogn_size <- function(n, p) {
  as.integer(min(p, max(2L, floor(n / log(n)))))
}

sb_sis_binary_screen_size <- function(n, p) {
  sb_high_dimensional_nlogn_size(n, p)
}

sb_fit_sis_binary_selector <- function(x, outcome) {
  method <- "sis_binary_nlogn"
  sb_require_high_dimensional_method(method)
  inputs <- sb_validate_high_dimensional_binary_inputs(x, outcome, method)
  if (ncol(inputs$x) < 2L) {
    sb_abort(
      "SIS screening requires at least two candidate predictors.",
      "stormbreaker_sis_too_few_predictors"
    )
  }

  params <- sb_sis_binary_selector_params()
  screen_size <- sb_sis_binary_screen_size(nrow(inputs$x), ncol(inputs$x))
  native <- sb_high_dimensional_native_call(
    withr::with_seed(params$seed, SIS::SIS(
      x = data.matrix(inputs$x),
      y = inputs$outcome,
      family = params$family,
      penalty = params$penalty,
      tune = params$tuning,
      type.measure = "deviance",
      nsis = screen_size,
      iter = params$iterative,
      varISIS = "vanilla",
      perm = FALSE,
      greedy = FALSE,
      seed = params$seed,
      standardize = params$standardize,
      boot_ci = FALSE,
      parallel = params$parallel
    )),
    method,
    "fit"
  )

  selected_index <- as.integer(native$sis.ix0)
  ranking <- if (is.list(native$ix_list) && length(native$ix_list) >= 1L) {
    as.integer(native$ix_list[[1L]])
  } else {
    integer()
  }
  p <- ncol(inputs$x)
  valid <- length(selected_index) == screen_size &&
    !anyNA(selected_index) && !anyDuplicated(selected_index) &&
    all(selected_index >= 1L & selected_index <= p) &&
    length(ranking) == p && !anyNA(ranking) &&
    !anyDuplicated(ranking) && setequal(ranking, seq_len(p))
  if (!valid) {
    sb_abort(
      "SIS did not return a complete marginal screening ranking.",
      "stormbreaker_invalid_sis_selector"
    )
  }

  source_features <- names(inputs$x)
  selected_features <- source_features[selected_index]
  active_params <- params
  active_params$screen_size <- screen_size
  structure(
    list(
      id = method,
      source_features = source_features,
      selected_features = selected_features,
      params = active_params,
      scores = data.frame(
        feature = source_features,
        rank = as.integer(match(seq_len(p), ranking)),
        selected = source_features %in% selected_features,
        stringsAsFactors = FALSE
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_vsurf_binary_selector_params <- function() {
  list(
    method = "VSURF",
    selected_subset = "prediction",
    random_forest_engine = "randomForest",
    mtry_rule = "max(floor(p/3), 1)",
    threshold_trees = 500L,
    threshold_forests = 20L,
    interpretation_trees = 100L,
    interpretation_forests = 10L,
    prediction_trees = 100L,
    prediction_forests = 10L,
    nmin = 1L,
    nsd = 1,
    nmj = 1,
    parallel = FALSE,
    seed = 1L
  )
}

sb_fit_vsurf_binary_selector <- function(x, outcome) {
  method <- "vsurf_binary_prediction"
  inputs <- sb_validate_high_dimensional_binary_inputs(x, outcome, method)
  sb_fit_vsurf_classification_selector(inputs$x, factor(inputs$outcome, levels = c(0L, 1L)),
    "vsurf_binary_prediction")
}

sb_fit_vsurf_classification_selector <- function(x, response, method) {
  sb_require_high_dimensional_method("vsurf_binary_prediction")
  inputs <- list(x = x)
  if (ncol(inputs$x) < 2L) {
    sb_abort(
      "VSURF selection requires at least two candidate predictors.",
      "stormbreaker_vsurf_too_few_predictors"
    )
  }

  params <- sb_vsurf_binary_selector_params()
  mtry <- max(floor(ncol(inputs$x) / 3), 1L)
  native <- sb_high_dimensional_native_call(
    withr::with_seed(params$seed, VSURF::VSURF(
      x = data.matrix(inputs$x),
      y = response,
      mtry = mtry,
      ntree.thres = params$threshold_trees,
      nfor.thres = params$threshold_forests,
      nmin = params$nmin,
      ntree.interp = params$interpretation_trees,
      nfor.interp = params$interpretation_forests,
      nsd = params$nsd,
      ntree.pred = params$prediction_trees,
      nfor.pred = params$prediction_forests,
      nmj = params$nmj,
      RFimplem = params$random_forest_engine,
      parallel = params$parallel,
      verbose = FALSE
    )),
    method,
    "fit"
  )

  selected_index <- as.integer(native$varselect.pred)
  ranking <- as.integer(native$imp.mean.dec.ind)
  importance <- as.numeric(native$imp.mean.dec)
  importance_sd <- as.numeric(native$imp.sd.dec)
  p <- ncol(inputs$x)
  valid <- length(selected_index) > 0L && !anyNA(selected_index) &&
    !anyDuplicated(selected_index) &&
    all(selected_index >= 1L & selected_index <= p) &&
    length(ranking) == p && !anyNA(ranking) &&
    !anyDuplicated(ranking) && setequal(ranking, seq_len(p)) &&
    length(importance) == p && length(importance_sd) == p &&
    !anyNA(importance) && !anyNA(importance_sd) &&
    all(is.finite(importance)) && all(is.finite(importance_sd))
  if (!valid) {
    sb_abort(
      "VSURF did not return a valid prediction-stage feature subset.",
      "stormbreaker_invalid_vsurf_selector"
    )
  }

  source_features <- names(inputs$x)
  selected_features <- source_features[selected_index]
  rank_by_feature <- integer(p)
  mean_by_feature <- numeric(p)
  sd_by_feature <- numeric(p)
  rank_by_feature[ranking] <- seq_len(p)
  mean_by_feature[ranking] <- importance
  sd_by_feature[ranking] <- importance_sd
  active_params <- params
  active_params$mtry <- as.integer(mtry)
  structure(
    list(
      id = method,
      source_features = source_features,
      selected_features = selected_features,
      params = active_params,
      scores = data.frame(
        feature = source_features,
        importance = mean_by_feature,
        importance_sd = sd_by_feature,
        rank = as.integer(rank_by_feature),
        selected = source_features %in% selected_features,
        stringsAsFactors = FALSE
      ),
      selection_history = data.frame(
        stage = c("threshold", "interpretation", "prediction"),
        n_features = as.integer(c(
          length(native$varselect.thres),
          length(native$varselect.interp),
          length(native$varselect.pred)
        )),
        stringsAsFactors = FALSE
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_mrmr_binary_selector_params <- function() {
  list(
    method = "minimum_redundancy_maximum_relevance",
    screen_size_rule = "min(p, max(2, floor(n/log(n))))",
    positive_scores_only = FALSE,
    threads = 1L
  )
}

sb_fit_mrmr_binary_selector <- function(x, outcome) {
  method <- "mrmr_binary_nlogn"
  sb_require_high_dimensional_method(method)
  inputs <- sb_validate_high_dimensional_binary_inputs(x, outcome, method)
  if (ncol(inputs$x) < 2L) {
    sb_abort(
      "mRMR screening requires at least two candidate predictors.",
      "stormbreaker_mrmr_too_few_predictors"
    )
  }

  params <- sb_mrmr_binary_selector_params()
  screen_size <- sb_high_dimensional_nlogn_size(
    nrow(inputs$x),
    ncol(inputs$x)
  )
  native <- sb_high_dimensional_native_call(
    praznik::MRMR(
      X = inputs$x,
      Y = factor(inputs$outcome, levels = c(0L, 1L)),
      k = screen_size,
      positive = params$positive_scores_only,
      threads = params$threads
    ),
    method,
    "fit"
  )

  selected_index <- as.integer(native$selection)
  selected_score <- as.numeric(native$score)
  p <- ncol(inputs$x)
  valid <- length(selected_index) == screen_size &&
    !anyNA(selected_index) && !anyDuplicated(selected_index) &&
    all(selected_index >= 1L & selected_index <= p) &&
    length(selected_score) == screen_size && !anyNA(selected_score) &&
    all(is.finite(selected_score))
  if (!valid) {
    sb_abort(
      "mRMR did not return a valid ordered feature subset.",
      "stormbreaker_invalid_mrmr_selector"
    )
  }

  source_features <- names(inputs$x)
  selected_features <- source_features[selected_index]
  score_by_feature <- rep(NA_real_, p)
  score_by_feature[selected_index] <- selected_score
  active_params <- params
  active_params$screen_size <- screen_size
  structure(
    list(
      id = method,
      source_features = source_features,
      selected_features = selected_features,
      params = active_params,
      scores = data.frame(
        feature = source_features,
        mrmr_score = score_by_feature,
        rank = as.integer(match(seq_len(p), selected_index)),
        selected = source_features %in% selected_features,
        stringsAsFactors = FALSE
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_information_binary_selector_params <- function(selector) {
  native_methods <- c(
    jmi_binary_nlogn = "JMI",
    cmim_binary_nlogn = "CMIM",
    disr_binary_nlogn = "DISR",
    jmim_binary_nlogn = "JMIM"
  )
  criteria <- c(
    jmi_binary_nlogn = "joint_mutual_information",
    cmim_binary_nlogn = "conditional_mutual_information_maximization",
    disr_binary_nlogn = "double_input_symmetrical_relevance",
    jmim_binary_nlogn = "minimal_joint_mutual_information_maximization"
  )
  if (!selector %in% names(native_methods)) {
    sb_abort(
      paste0("Unknown information selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  }
  list(
    method = unname(criteria[[selector]]),
    native_method = unname(native_methods[[selector]]),
    screen_size_rule = "min(p, max(2, floor(n/log(n))))",
    threads = 1L
  )
}

sb_run_information_binary_selector <- function(selector, x, outcome,
                                               screen_size, threads) {
  switch(selector,
    jmi_binary_nlogn = praznik::JMI(
      X = x,
      Y = outcome,
      k = screen_size,
      threads = threads
    ),
    cmim_binary_nlogn = praznik::CMIM(
      X = x,
      Y = outcome,
      k = screen_size,
      threads = threads
    ),
    disr_binary_nlogn = praznik::DISR(
      X = x,
      Y = outcome,
      k = screen_size,
      threads = threads
    ),
    jmim_binary_nlogn = praznik::JMIM(
      X = x,
      Y = outcome,
      k = screen_size,
      threads = threads
    )
  )
}

sb_fit_information_binary_selector <- function(x, outcome, selector) {
  params <- sb_information_binary_selector_params(selector)
  sb_require_high_dimensional_method(selector)
  inputs <- sb_validate_high_dimensional_binary_inputs(x, outcome, selector)
  if (ncol(inputs$x) < 2L) {
    sb_abort(
      paste0(params$native_method, " screening requires at least two candidate predictors."),
      "stormbreaker_information_selector_too_few_predictors"
    )
  }

  screen_size <- sb_high_dimensional_nlogn_size(
    nrow(inputs$x),
    ncol(inputs$x)
  )
  native <- sb_high_dimensional_native_call(
    sb_run_information_binary_selector(
      selector,
      inputs$x,
      factor(inputs$outcome, levels = c(0L, 1L)),
      screen_size,
      params$threads
    ),
    selector,
    "fit"
  )

  selected_index <- as.integer(native$selection)
  selected_score <- as.numeric(native$score)
  p <- ncol(inputs$x)
  valid <- length(selected_index) == screen_size &&
    !anyNA(selected_index) && !anyDuplicated(selected_index) &&
    all(selected_index >= 1L & selected_index <= p) &&
    length(selected_score) == screen_size && !anyNA(selected_score) &&
    all(is.finite(selected_score))
  if (!valid) {
    sb_abort(
      paste0(
        params$native_method,
        " did not return a valid ordered feature subset."
      ),
      "stormbreaker_invalid_information_selector"
    )
  }

  source_features <- names(inputs$x)
  selected_features <- source_features[selected_index]
  score_by_feature <- rep(NA_real_, p)
  score_by_feature[selected_index] <- selected_score
  active_params <- params
  active_params$screen_size <- screen_size
  structure(
    list(
      id = selector,
      source_features = source_features,
      selected_features = selected_features,
      params = active_params,
      scores = data.frame(
        feature = source_features,
        information_score = score_by_feature,
        rank = as.integer(match(seq_len(p), selected_index)),
        selected = source_features %in% selected_features,
        stringsAsFactors = FALSE
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_lpc_selector_params <- function(selector) {
  if (!selector %in% sb_lpc_selector_ids()) {
    sb_abort(
      paste0("Unknown LPC selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  }
  list(
    method = "lassoed_principal_components",
    native_type = if (identical(selector, "lpc_binary_nlogn")) {
      "two class"
    } else {
      "survival"
    },
    ranking = "decreasing_absolute_lpc_score",
    screen_size_rule = "min(p, max(2, floor(n/log(n))))",
    soft_threshold = "native_10_half_split_predictive_advantage",
    seed = 1L
  )
}

sb_fit_lpc_selector <- function(x, endpoint, selector) {
  sb_require_high_dimensional_method(selector)
  params <- sb_lpc_selector_params(selector)
  inputs <- if (identical(selector, "lpc_binary_nlogn")) {
    binary <- sb_validate_high_dimensional_binary_inputs(
      x,
      endpoint$outcome,
      selector
    )
    list(
      x = binary$x,
      response = binary$outcome + 1L,
      censoring_status = NULL
    )
  } else {
    survival <- sb_validate_classical_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      selector
    )
    list(
      x = survival$x,
      response = survival$time,
      censoring_status = survival$event
    )
  }
  screen_size <- sb_high_dimensional_nlogn_size(
    nrow(inputs$x),
    ncol(inputs$x)
  )
  native <- sb_high_dimensional_native_call(
    {
      fitted <- NULL
      invisible(utils::capture.output(
        fitted <- withr::with_seed(
          params$seed,
          lpc::LPC(
            x = t(data.matrix(inputs$x)),
            y = inputs$response,
            type = params$native_type,
            soft.thresh = NULL,
            censoring.status = inputs$censoring_status
          )
        ),
        type = "output"
      ))
      fitted
    },
    selector,
    "fit"
  )

  lpc_score <- as.numeric(native$lpcscores)
  base_score <- as.numeric(native$tscores)
  soft_threshold <- as.numeric(native$soft.thresh)
  p <- ncol(inputs$x)
  valid <- length(lpc_score) == p && !anyNA(lpc_score) &&
    all(is.finite(lpc_score)) &&
    length(base_score) == p && !anyNA(base_score) &&
    all(is.finite(base_score)) &&
    length(soft_threshold) == 1L && !is.na(soft_threshold) &&
    is.finite(soft_threshold) && soft_threshold >= 0
  if (!valid) {
    sb_abort(
      if (p < 50L && length(soft_threshold) == 0L) {
        paste0("Native LPC automatic threshold search requires at least 50 predictors; ",
               "the supplied ", p, " predictors produced an empty threshold. ",
               "No padding or substitute threshold was used.")
      } else {
        "LPC did not return finite gene scores and one soft threshold."
      },
      if (p < 50L && length(soft_threshold) == 0L) {
        c("stormbreaker_lpc_dimension_unsupported", "stormbreaker_invalid_lpc_selector")
      } else "stormbreaker_invalid_lpc_selector"
    )
  }

  source_features <- names(inputs$x)
  ranking <- order(-abs(lpc_score), seq_len(p), method = "radix")
  selected_index <- ranking[seq_len(screen_size)]
  selected_features <- source_features[selected_index]
  active_params <- params
  active_params$screen_size <- screen_size
  active_params$realized_soft_threshold <- soft_threshold
  structure(
    list(
      id = selector,
      source_features = source_features,
      selected_features = selected_features,
      params = active_params,
      scores = data.frame(
        feature = source_features,
        lpc_score = lpc_score,
        base_score = base_score,
        rank = as.integer(match(seq_len(p), ranking)),
        selected = source_features %in% selected_features,
        stringsAsFactors = FALSE
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_fit_lpc_binary_selector <- function(x, outcome) {
  sb_fit_lpc_selector(
    x,
    list(outcome = outcome),
    "lpc_binary_nlogn"
  )
}

sb_fit_lpc_cox_selector <- function(x, time, event) {
  sb_fit_lpc_selector(
    x,
    list(time = time, event = event),
    "lpc_cox_nlogn"
  )
}

sb_stabsel_selector_params <- function(selector) {
  if (!selector %in% c("stabsel_binary_lasso", "stabsel_cox_lasso")) {
    sb_abort(
      paste0("Unknown stability selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  }
  params <- list(
    base_selector = "glmnet_lasso",
    family = if (identical(selector, "stabsel_binary_lasso")) {
      "binomial"
    } else {
      "cox"
    },
    subsample_q_rule = paste0(
      "min(max(2, floor(n/log(n))), floor(p/2)) with p >= 4"
    ),
    base_selection_type = "conservative",
    cutoff = 0.75,
    complementary_pairs = 50L,
    sampling_type = "SS",
    assumption = "unimodal",
    standardize = TRUE,
    parallel = FALSE,
    seed = 1L
  )
  if (identical(selector, "stabsel_cox_lasso")) {
    params$cox_ties <- "breslow"
  }
  params
}

sb_stabsel_subsample_q <- function(n, p) {
  as.integer(min(max(2L, floor(n / log(n))), floor(p / 2)))
}

sb_fit_stabsel_selector <- function(x, endpoint, selector) {
  sb_require_high_dimensional_method(selector)
  inputs <- if (identical(selector, "stabsel_binary_lasso")) {
    binary <- sb_validate_high_dimensional_binary_inputs(
      x,
      endpoint$outcome,
      selector
    )
    list(
      x = binary$x,
      response = binary$outcome
    )
  } else {
    survival <- sb_validate_classical_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      selector
    )
    list(
      x = survival$x,
      response = survival::Surv(survival$time, survival$event)
    )
  }
  if (ncol(inputs$x) < 4L) {
    sb_abort(
      "Stability selection requires at least four candidate predictors.",
      "stormbreaker_stabsel_too_few_predictors"
    )
  }

  params <- sb_stabsel_selector_params(selector)
  subsample_q <- sb_stabsel_subsample_q(
    nrow(inputs$x),
    ncol(inputs$x)
  )
  fitfun_args <- list(
    type = params$base_selection_type,
    family = params$family,
    standardize = params$standardize
  )
  if (identical(params$family, "cox")) {
    fitfun_args$cox.ties <- params$cox_ties
  }
  native <- sb_high_dimensional_native_call(
    withr::with_seed(params$seed, stabs::stabsel(
      x = data.matrix(inputs$x),
      y = inputs$response,
      fitfun = stabs::glmnet.lasso,
      args.fitfun = fitfun_args,
      cutoff = params$cutoff,
      q = subsample_q,
      B = params$complementary_pairs,
      sampling.type = params$sampling_type,
      assumption = params$assumption,
      papply = parallel::mclapply,
      mc.cores = 1L,
      verbose = FALSE
    )),
    selector,
    "fit"
  )

  selected_index <- unname(as.integer(stabs::selected(native)))
  source_features <- names(inputs$x)
  probabilities <- native$max[source_features]
  p <- ncol(inputs$x)
  valid <- length(probabilities) == p && !anyNA(probabilities) &&
    all(is.finite(probabilities)) &&
    all(probabilities >= 0 & probabilities <= 1) &&
    !anyNA(selected_index) && !anyDuplicated(selected_index) &&
    all(selected_index >= 1L & selected_index <= p) &&
    is.numeric(native$PFER) && length(native$PFER) == 1L &&
    is.finite(native$PFER) && native$PFER >= 0
  if (!valid) {
    sb_abort(
      "Stability selection did not return valid selection probabilities.",
      "stormbreaker_invalid_stabsel_selector"
    )
  }
  if (length(selected_index) == 0L) {
    sb_abort(
      "Stability selection did not retain any predictors at cutoff 0.75.",
      "stormbreaker_empty_stabsel_selector"
    )
  }

  selected_features <- source_features[selected_index]
  active_params <- params
  active_params$subsample_q <- subsample_q
  active_params$realized_pfer <- unname(as.numeric(native$PFER))
  structure(
    list(
      id = selector,
      source_features = source_features,
      selected_features = selected_features,
      params = active_params,
      scores = data.frame(
        feature = source_features,
        selection_probability = unname(as.numeric(probabilities)),
        rank = as.integer(rank(-probabilities, ties.method = "min")),
        selected = source_features %in% selected_features,
        stringsAsFactors = FALSE
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_fit_stabsel_binary_selector <- function(x, outcome) {
  sb_fit_stabsel_selector(
    x,
    list(outcome = outcome),
    "stabsel_binary_lasso"
  )
}

sb_fit_stabsel_cox_selector <- function(x, time, event) {
  sb_fit_stabsel_selector(
    x,
    list(time = time, event = event),
    "stabsel_cox_lasso"
  )
}

sb_apply_high_dimensional_selector <- function(fitted, x) {
  valid <- inherits(fitted, "stormbreaker_selector_fit") &&
    fitted$id %in% sb_high_dimensional_selector_ids() &&
    is.character(fitted$source_features) &&
    is.character(fitted$selected_features) &&
    length(fitted$selected_features) > 0L &&
    !anyNA(fitted$selected_features) &&
    !anyDuplicated(fitted$selected_features) &&
    all(fitted$selected_features %in% fitted$source_features)
  if (!valid) {
    sb_abort(
      "The fitted high-dimensional selector state is invalid.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      "Predictor columns do not match the fitted selector schema.",
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_high_dimensional_selector_adapter <- function(selector) {
  if (!selector %in% sb_high_dimensional_selector_ids()) {
    sb_abort(
      paste0("Unknown high-dimensional selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  }
  fit <- if (selector %in% sb_information_binary_selector_ids()) {
    function(x, outcome) {
      sb_fit_information_binary_selector(x, outcome, selector)
    }
  } else {
    switch(selector,
      sis_binary_nlogn = sb_fit_sis_binary_selector,
      vsurf_binary_prediction = sb_fit_vsurf_binary_selector,
      mrmr_binary_nlogn = sb_fit_mrmr_binary_selector,
      lpc_binary_nlogn = sb_fit_lpc_binary_selector,
      lpc_cox_nlogn = sb_fit_lpc_cox_selector,
      stabsel_binary_lasso = sb_fit_stabsel_binary_selector,
      stabsel_cox_lasso = sb_fit_stabsel_cox_selector
    )
  }
  list(
    id = selector,
    params = list(),
    needs_outcome = TRUE,
    fit = fit,
    apply = sb_apply_high_dimensional_selector
  )
}

sb_hdda_binary_params <- function(learner = "hdda_binary_cattell") {
  if (!learner %in% sb_hdda_binary_learner_ids()) {
    sb_abort(
      paste0("Unknown HDDA binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    model = "AkjBkQkDk",
    dimension_selection = "Cattell",
    threshold = 0.2,
    scaling = TRUE,
    noise_control = 1e-8,
    probability_mapping = "native_class_one_posterior"
  )
}

sb_hdda_binary_parameter_semantics <- function(
  learner = "hdda_binary_cattell"
) {
  sb_hdda_binary_params(learner)
  paste0(
    "HDDA AkjBkQkDk class-specific subspaces, Cattell dimension rule, ",
    "threshold 0.2, analysis-fold scaling, and native posterior probabilities"
  )
}

sb_fit_hdda_binary <- function(x, outcome, learner, params) {
  sb_require_high_dimensional_method(learner)
  inputs <- sb_validate_high_dimensional_binary_inputs(x, outcome, learner)
  engine <- sb_high_dimensional_native_call(
    HDclassif::hdda(
      data = data.matrix(inputs$x),
      cls = factor(inputs$outcome, levels = c(0L, 1L)),
      model = params$model,
      graph = FALSE,
      d_select = params$dimension_selection,
      threshold = params$threshold,
      show = FALSE,
      scaling = params$scaling,
      noise.ctrl = params$noise_control
    ),
    learner,
    "fit"
  )
  valid <- inherits(engine, "hdc") && identical(engine$K, 2L) &&
    identical(engine$kname, c("0", "1"))
  if (!valid) {
    sb_abort(
      "HDDA did not return a fitted two-class discriminant model.",
      "stormbreaker_invalid_hdda_fit"
    )
  }
  structure(
    list(
      learner = learner,
      features = names(inputs$x),
      params = params,
      engine = engine
    ),
    class = "stormbreaker_hdda_binary"
  )
}

sb_predict_hdda_binary <- function(model, x) {
  if (!inherits(model, "stormbreaker_hdda_binary") ||
    !model$learner %in% sb_hdda_binary_learner_ids()) {
    sb_abort(
      "The fitted HDDA binary model is invalid.",
      "stormbreaker_bad_hdda_model"
    )
  }
  sb_require_high_dimensional_method(model$learner)
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted HDDA model.",
      "stormbreaker_hdda_schema_mismatch"
    )
  }
  prediction <- sb_high_dimensional_native_call(
    stats::predict(model$engine, data = data.matrix(x)),
    model$learner,
    "prediction"
  )
  posterior <- prediction$posterior
  probability <- if (is.matrix(posterior) && "1" %in% colnames(posterior)) {
    unname(as.numeric(posterior[, "1"]))
  } else {
    numeric()
  }
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      "HDDA produced invalid class-one posterior probabilities.",
      "stormbreaker_invalid_binary_prediction"
    )
  }
  probability
}

sb_extract_hdda_binary_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_hdda_binary")) {
    sb_bad_learner_state("The fitted HDDA binary state is invalid.")
  }
  structure(
    unclass(model),
    class = c(
      "stormbreaker_hdda_binary_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_hdda_binary_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_hdda_binary_learner_ids()) {
    sb_bad_learner_state("The saved HDDA binary state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- structure(unclass(state), class = "stormbreaker_hdda_binary")
  list(probability = sb_predict_hdda_binary(model, x))
}

sb_hdda_binary_learner_adapter <- function(
  learner,
  params = sb_hdda_binary_params(learner)
) {
  if (!learner %in% sb_hdda_binary_learner_ids()) {
    sb_abort(
      paste0("Unknown HDDA binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_hdda_binary(x, outcome, learner, params)
    },
    predict = sb_predict_hdda_binary,
    extract_state = sb_extract_hdda_binary_state,
    predict_state = sb_predict_hdda_binary_state
  )
}

# ---- cancerclass-binary-adapter ----

sb_cancerclass_available <- function() {
  requireNamespace("cancerclass", quietly = TRUE) &&
    requireNamespace("Biobase", quietly = TRUE)
}

sb_cancerclass_binary_params <- function() {
  list(
    feature_method = "welch.test",
    hparam = 0.75,
    distance = "cor",
    feature_count = "all",
    score = "z",
    calibration = "oof_logistic",
    inner_folds = 5L,
    inner_seed = 1L
  )
}

sb_cancerclass_expression_set <- function(x, outcome) {
  sample_names <- paste0("sample_", seq_len(nrow(x)))
  expression <- t(data.matrix(x))
  colnames(expression) <- sample_names
  phenotype <- data.frame(
    class = factor(outcome, levels = c(0L, 1L)),
    row.names = sample_names
  )
  Biobase::ExpressionSet(
    assayData = expression,
    phenoData = Biobase::AnnotatedDataFrame(phenotype)
  )
}

sb_fit_cancerclass_native_predictor <- function(x, outcome, params) {
  search_before <- search()
  withr::local_package("Biobase", quietly = TRUE)
  attached_dependencies <- setdiff(
    search(),
    c(search_before, "package:Biobase")
  )
  withr::defer({
    for (package in attached_dependencies) {
      if (package %in% search()) {
        base::detach(package, character.only = TRUE)
      }
    }
  }, priority = "last")
  fitted <- suppressMessages(cancerclass::fit(
    sb_cancerclass_expression_set(x, outcome),
    class = "class",
    method = params$feature_method,
    hparam = params$hparam
  ))
  predictor <- if (methods::is(fitted, "predictor")) {
    methods::slot(fitted, "predictor")
  } else {
    NULL
  }
  valid <- is.matrix(predictor) && nrow(predictor) == ncol(x) &&
    setequal(rownames(predictor), names(x)) &&
    all(c(params$feature_method, "0", "1") %in% colnames(predictor))
  if (!valid) {
    sb_abort(
      "`cancerclass` did not return one usable centroid row per feature.",
      "stormbreaker_invalid_cancerclass_fit"
    )
  }
  predictor
}

sb_cancerclass_score_from_predictor <- function(
    predictor, x, features, distance) {
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), features)) {
    sb_abort(
      "Prediction columns do not match the fixed cancerclass schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  data <- t(data.matrix(x))
  data <- data[rownames(predictor), , drop = FALSE]
  distance_zero <- cancerclass::get.d(
    data,
    predictor[, "0"],
    method = distance
  )
  distance_one <- cancerclass::get.d(
    data,
    predictor[, "1"],
    method = distance
  )
  class_distance <- cancerclass::get.d(
    predictor[, "0"],
    predictor[, "1"],
    method = distance
  )
  score <- (distance_zero - distance_one) / class_distance
  valid <- length(class_distance) == 1L && !is.na(class_distance) &&
    is.finite(class_distance) && class_distance > 0 &&
    length(score) == nrow(x) && !anyNA(score) && all(is.finite(score))
  if (!valid) {
    sb_abort(
      "`cancerclass` produced an invalid native correlation-distance score.",
      "stormbreaker_invalid_cancerclass_score"
    )
  }
  unname(score)
}

sb_fit_cancerclass_binary <- function(x, outcome, params) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "A fixed cancerclass analysis fold must contain both classes.",
      "stormbreaker_cancerclass_single_class"
    )
  }
  if (!sb_cancerclass_available()) {
    sb_abort(
      paste0(
        "Packages `cancerclass` and `Biobase` are required for the fixed ",
        "cancerclass learner."
      ),
      "stormbreaker_missing_cancerclass"
    )
  }

  fold_id <- sb_make_stratified_fold_id(
    outcome,
    folds = params$inner_folds,
    seed = params$inner_seed,
    stratum_label = "class"
  )
  oof_score <- numeric(nrow(x))
  for (fold in seq_len(params$inner_folds)) {
    analysis_rows <- which(fold_id != fold)
    assessment_rows <- which(fold_id == fold)
    predictor <- sb_fit_cancerclass_native_predictor(
      x[analysis_rows, , drop = FALSE],
      outcome[analysis_rows],
      params
    )
    oof_score[assessment_rows] <- sb_cancerclass_score_from_predictor(
      predictor,
      x[assessment_rows, , drop = FALSE],
      names(x),
      params$distance
    )
  }

  calibration_warnings <- character()
  calibration <- withCallingHandlers(
    stats::glm.fit(
      x = cbind("(Intercept)" = 1, native_z_score = oof_score),
      y = outcome,
      family = stats::binomial()
    ),
    warning = function(warning) {
      calibration_warnings <<- c(
        calibration_warnings,
        conditionMessage(warning)
      )
      invokeRestart("muffleWarning")
    }
  )
  coefficients <- unname(as.numeric(calibration$coefficients))
  valid_calibration <- isTRUE(calibration$converged) &&
    length(coefficients) == 2L && !anyNA(coefficients) &&
    all(is.finite(coefficients))
  if (!valid_calibration) {
    sb_abort(
      "The inner out-of-fold cancerclass probability calibration failed.",
      "stormbreaker_cancerclass_calibration_failure"
    )
  }
  names(coefficients) <- c("intercept", "native_z_score")
  final_predictor <- sb_fit_cancerclass_native_predictor(x, outcome, params)

  structure(
    list(
      features = names(x),
      params = params,
      predictor = final_predictor,
      calibration_coefficients = coefficients,
      calibration_fold_id = fold_id,
      calibration_oof_score = oof_score,
      calibration_warnings = unique(calibration_warnings)
    ),
    class = "stormbreaker_cancerclass_binary"
  )
}

sb_cancerclass_raw_score <- function(model, x) {
  if (!inherits(model, "stormbreaker_cancerclass_binary")) {
    sb_abort(
      "The supplied model is not a fitted fixed cancerclass learner.",
      "stormbreaker_bad_cancerclass_fit"
    )
  }
  if (!sb_cancerclass_available()) {
    sb_abort(
      paste0(
        "Packages `cancerclass` and `Biobase` are required for fixed ",
        "cancerclass prediction."
      ),
      "stormbreaker_missing_cancerclass"
    )
  }
  sb_cancerclass_score_from_predictor(
    model$predictor,
    x,
    model$features,
    model$params$distance
  )
}

sb_predict_cancerclass_binary <- function(model, x) {
  score <- sb_cancerclass_raw_score(model, x)
  probability <- stats::plogis(
    model$calibration_coefficients[[1L]] +
      model$calibration_coefficients[[2L]] * score
  )
  if (length(probability) != length(score) || anyNA(probability) ||
    any(!is.finite(probability)) ||
    any(probability < 0 | probability > 1)) {
    sb_abort(
      "The fixed cancerclass learner produced invalid probabilities.",
      "stormbreaker_invalid_cancerclass_probability"
    )
  }
  unname(probability)
}

sb_extract_cancerclass_calibration_coefficients <- function(model) {
  data.frame(
    term = c("(Intercept)", "native_z_score"),
    coefficient = unname(model$calibration_coefficients),
    stringsAsFactors = FALSE
  )
}

sb_extract_cancerclass_binary_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_cancerclass_binary")) {
    sb_bad_learner_state("The fitted cancerclass learner state is invalid.")
  }
  structure(
    list(
      learner = "cancerclass_binary",
      features = model$features,
      params = model$params,
      predictor = model$predictor,
      calibration_coefficients = model$calibration_coefficients,
      score_reference = "cancerclass_native_welch_correlation_z",
      probability_formula = "plogis(intercept + slope * native_z_score)",
      training_data_retention = "none"
    ),
    class = c(
      "stormbreaker_cancerclass_binary_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_cancerclass_binary_learner_state <- function(state, x) {
  x <- sb_validate_learner_state_x(state, x, "cancerclass_binary")
  valid <- is.matrix(state[["predictor"]]) &&
    is.numeric(state[["calibration_coefficients"]]) &&
    length(state[["calibration_coefficients"]]) == 2L &&
    !anyNA(state[["calibration_coefficients"]]) &&
    all(is.finite(state[["calibration_coefficients"]]))
  if (!valid) {
    sb_bad_learner_state("The saved cancerclass learner state is invalid.")
  }
  model <- structure(
    list(
      features = state[["features"]],
      params = state[["params"]],
      predictor = state[["predictor"]],
      calibration_coefficients = state[["calibration_coefficients"]]
    ),
    class = "stormbreaker_cancerclass_binary"
  )
  list(probability = sb_predict_cancerclass_binary(model, x))
}

sb_cancerclass_binary_learner_adapter <- function() {
  params <- sb_cancerclass_binary_params()
  list(
    id = "cancerclass_binary",
    params = params,
    fit = function(x, outcome) {
      sb_fit_cancerclass_binary(x, outcome, params)
    },
    predict = sb_predict_cancerclass_binary,
    extract_coefficients = sb_extract_cancerclass_calibration_coefficients,
    extract_state = sb_extract_cancerclass_binary_state,
    predict_state = sb_predict_cancerclass_binary_learner_state
  )
}

# ---- diagnostics-omics ----

# ---- algorithm-cancerclass ----

#' @rdname sb_pls_diagnostics
#' @export
sb_cancerclass_diagnostics <- function(fit, type = "centroids", score_range = c(-1, 1)) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (s$learner != "cancerclass_binary" || !is.character(type) ||
      length(type) != 1L || is.na(type) || !type %in% c("centroids", "score_mapping")) sb_abort(
    "Choose centroids or score_mapping for Cancerclass.", "stormbreaker_unsupported_diagnostic")
  if (!is.matrix(s$predictor) || !setequal(rownames(s$predictor), s$features) ||
      !all(c("0", "1") %in% colnames(s$predictor))) sb_abort(
    "Original Cancerclass centroids are invalid.", "stormbreaker_invalid_diagnostic")
  if (type == "centroids") {
    m <- s$predictor[s$features, c("0", "1"), drop = FALSE]
    tab <- data.frame(feature = rep(s$features, 2L), class = rep(c("0", "1"), each = length(s$features)), mean = as.numeric(m))
  } else {
    if (!is.numeric(score_range) || length(score_range) != 2L ||
        any(!is.finite(score_range)) || score_range[2L] <= score_range[1L]) sb_abort(
      "score_range must contain two finite increasing values.", "stormbreaker_invalid_diagnostic")
    score <- seq(score_range[1L], score_range[2L], length.out = 201L)
    b <- s$calibration_coefficients
    tab <- data.frame(score = score, probability = stats::plogis(b[[1L]] + b[[2L]] * score))
  }
  structure(list(type = paste0("cancerclass_", type), table = tab,
    calibration_coefficients = s$calibration_coefficients,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original saved Cancerclass centroids and existing OOF-fitted score mapping",
    interpretation = "The score-mapping curve is the existing probability formula, not a calibration-reliability assessment. Its domain is an explicit display range, not observed support. No OOF reconstruction, refit or surrogate."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_cancerclass_diagnostic <- function(x) {
  mean <- feature <- class <- score <- probability <- NULL
  if (x$type == "cancerclass_centroids") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(mean, feature, colour = class)) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#B86236")) +
      ggplot2::labs(title = "Cancerclass centroids", x = "Centroid value", y = NULL, colour = "Class")
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(score, probability)) +
      ggplot2::geom_line(colour = "#2878A0", linewidth = .7) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = "Cancerclass probability mapping", x = "Native z score", y = "Class-1 probability")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 10, height = 4)
  p
}

# ---- algorithm-high-dimensional ----

#' @rdname sb_pls_diagnostics
#' @export
sb_high_dimensional_diagnostics <- function(fit, type = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% c("hda_binary", "hdda_binary_cattell")) sb_abort(
    "These diagnostics require HDA or HDDA.", "stormbreaker_unsupported_diagnostic")
  hda <- s$learner == "hda_binary"
  choices <- if (hda) c("loadings", "densities") else c("loadings", "spectrum", "centers")
  if (is.null(type)) type <- "loadings"
  if (!is.character(type) || length(type) != 1L || is.na(type) || !type %in% choices) sb_abort(
    paste("Choose", paste(choices, collapse = ", ")), "stormbreaker_unsupported_diagnostic")
  e <- s$engine
  p <- length(s$features)
  parameters <- NULL
  if (hda) {
    if (!inherits(e, "hda") || !identical(rownames(e$hda.loadings), s$features)) sb_abort(
      "Original HDA loadings are invalid.", "stormbreaker_invalid_diagnostic")
    k <- e$reduced.dimension
    if (type == "loadings") {
      m <- e$hda.loadings[, seq_len(k), drop = FALSE]
      tab <- data.frame(feature = rep(s$features, k), axis = rep(seq_len(k), each = p), loading = as.numeric(m), class = "Shared")
    } else {
      parameters <- do.call(rbind, lapply(seq_len(k), function(j) {
        m <- e$naivebayes$tables[[j]]
        data.frame(axis = j, class = rownames(m), mean = m[, 1L], sd = m[, 2L])
      }))
      tab <- do.call(rbind, lapply(seq_len(k), function(j) {
        z <- parameters[parameters$axis == j, ]
        grid <- seq(min(z$mean - 4 * z$sd), max(z$mean + 4 * z$sd), length.out = 201L)
        do.call(rbind, lapply(seq_len(nrow(z)), function(i) data.frame(axis = j, class = z$class[i],
          coordinate = grid, density = stats::dnorm(grid, z$mean[i], z$sd[i]))))
      }))
    }
  } else {
    if (!inherits(e, "hdc") || e$model != "AKJBKQKDK" || !identical(e$kname, c("0", "1"))) sb_abort(
      "Original HDDA model is invalid.", "stormbreaker_invalid_diagnostic")
    tab <- do.call(rbind, lapply(seq_along(e$kname), function(i) {
      k <- as.integer(e$d[i])
      if (type == "loadings") return(data.frame(class = e$kname[i],
        axis = rep(seq_len(k), each = p), feature = rep(s$features, k), loading = as.numeric(e$Q[[i]])))
      if (type == "spectrum") return(data.frame(class = e$kname[i], axis = seq_len(p),
        eigenvalue = as.numeric(e$ev[i, ]), retained = seq_len(p) <= k,
        model_variance = c(as.numeric(e$a[i, seq_len(k)]), rep(as.numeric(e$b[i]), p - k))))
      mu <- as.numeric(e$mu[i, ])
      if (length(e$scaling)) mu <- mu * e$scaling$sd + e$scaling$mu
      data.frame(class = e$kname[i], feature = s$features, mean = as.numeric(mu))
    }))
  }
  structure(list(type = paste0("high_dimensional_", type), table = tab, parameters = parameters,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original HDA loadings/naive-Bayes parameters or HDDA Q/ev/d/a/b/mu",
    interpretation = "Only retained native loadings are displayed. HDA densities use reduced coordinates; HDDA axes use standardized features and centers use original units. Signs are basis conventions, not effects or importance. No PCA recomputation or refitting."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_high_dimensional_diagnostic <- function(x) {
  loading <- feature <- coordinate <- density <- class <- axis <- eigenvalue <- retained <- mean <- NULL
  if (x$type == "high_dimensional_loadings") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(loading, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey65", linewidth = .3) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::facet_wrap(~class + axis, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "Retained discriminant loadings", x = "Loading", y = NULL)
  } else if (x$type == "high_dimensional_densities") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(coordinate, density, colour = class)) +
      ggplot2::geom_line(linewidth = .7) + ggplot2::facet_wrap(~axis, labeller = ggplot2::label_both) +
      ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#B86236")) +
      ggplot2::labs(title = "HDA reduced-space class densities", x = "Reduced coordinate", y = "Density", colour = "Class")
  } else if (x$type == "high_dimensional_spectrum") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(axis, eigenvalue)) +
      ggplot2::geom_line(colour = "grey60", linewidth = .4) + ggplot2::geom_point(ggplot2::aes(colour = retained), size = 2.5) +
      ggplot2::scale_colour_manual(values = c("FALSE" = "grey65", "TRUE" = "#2878A0"), labels = c("FALSE" = "Not retained", "TRUE" = "Retained")) +
      ggplot2::scale_x_continuous(breaks = unique(x$table$axis)) +
      ggplot2::facet_wrap(~class, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "HDDA class spectra", x = "Axis", y = "Eigenvalue", colour = NULL)
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(mean, feature)) + ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::facet_wrap(~class, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "HDDA class centers", x = "Mean in original units", y = NULL)
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}

# ---- algorithm-ktsp ----

#' @rdname sb_pls_diagnostics
#' @export
sb_ktsp_diagnostics <- function(fit) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (fit$pipeline$learner != "ktsp_binary") sb_abort(
    "This diagnostic requires the k-TSP binary learner.", "stormbreaker_unsupported_diagnostic")
  e <- fit$learner_state$engine
  if (!is.matrix(e$pairs) || ncol(e$pairs) != 2L || !nrow(e$pairs) ||
      !all(e$pairs %in% fit$learner_state$features) ||
      length(e$pair_scores) != nrow(e$pairs) || any(!is.finite(e$pair_scores)) ||
      !identical(e$labels, c("0", "1"))) sb_abort(
    "Original k-TSP pair mapping or scores are invalid.", "stormbreaker_bad_learner_state")
  tab <- data.frame(pair = seq_len(nrow(e$pairs)), first = e$pairs[, 1], second = e$pairs[, 2],
    native_score = e$pair_scores, class_1_rule = paste(e$pairs[, 1], ">", e$pairs[, 2]))
  structure(list(type = "ktsp_pair_scores", table = tab,
    vote_rule = c(greater = 1, equal = .5, smaller = 0), vote_aggregation = "unweighted_mean",
    unused_features = setdiff(fit$learner_state$features, as.vector(e$pairs)),
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    source = "original saved switchBox ordered pairs and pair selection scores",
    interpretation = paste("Each displayed comparison contributes one class-one vote when true, one half for equality, zero otherwise.",
      "Original pair order retained. Final output is the unweighted mean of votes; native scores are not voting weights.",
      "Scores describe training pair selection, not validation accuracy or uncertainty; vote fractions are not independently calibrated probabilities.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_ktsp_diagnostic <- function(x) {
  rule <- native_score <- NULL
  d <- x$table
  labels <- paste0("Pair ", d$pair, ": ", d$class_1_rule)
  d$rule <- factor(labels, levels = rev(labels))
  p <- ggplot2::ggplot(d, ggplot2::aes(native_score, rule)) +
    ggplot2::geom_point(size = 2.5, colour = "#2878A0") +
    ggplot2::labs(title = "k-TSP pair rules", x = "Native pair selection score", y = NULL) +
    sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = max(4.5, 2 + .4 * nrow(d)))
  p
}

# ---- algorithm-superpc ----

#' @rdname sb_pls_diagnostics
#' @export
sb_superpc_diagnostics <- function(fit, type = "feature_scores") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  id <- fit$pipeline$learner
  if (!id %in% c("superpc", "superpc_cv")) sb_abort(
    "This diagnostic requires an original SuperPC learner.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("feature_scores", "loadings", "cv_curve")) sb_abort(
    "Choose feature_scores, loadings or cv_curve.", "stormbreaker_unsupported_diagnostic")
  model <- fit$diagnostic_artifacts$learner
  if (is.null(model)) sb_abort("Original SuperPC objects were not captured; enable capture_diagnostics on final fitting.",
    "stormbreaker_missing_diagnostic_artifacts")
  engine <- model$engine
  features <- model$features
  selected <- engine$prediction$which.features
  threshold <- if (id == "superpc_cv") model$audit$selected_threshold else model$params$threshold
  if (!is.logical(selected) || length(selected) != length(features) || anyNA(selected)) sb_abort(
    "Original SuperPC feature selection is invalid.", "stormbreaker_invalid_diagnostic")
  folds <- NULL
  if (type == "feature_scores") {
    native <- engine$trained$feature.scores
    if (!is.numeric(native) || length(native) != length(features) || any(!is.finite(native)) ||
        !identical(names(native), features) || !identical(unname(abs(native) >= threshold), unname(selected))) sb_abort(
      "Native SuperPC scores disagree with the original selection.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(feature = features, score = as.numeric(native), selected = unname(selected))
    measure <- "Original SuperPC feature-screening scores"
    interpretation <- "Signed native survival-screening associations, not fitted coefficients or causal effects.\nDashed lines: original absolute-score threshold. All candidate features retained."
  } else if (type == "loadings") {
    native <- engine$prediction$u
    mat <- if (is.null(dim(native))) matrix(native, ncol = 1L) else native
    if (!is.matrix(mat) || nrow(mat) != sum(selected) ||
        ncol(mat) != engine$prediction$n.components || any(!is.finite(mat))) sb_abort(
      "Original SuperPC singular vectors are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(feature = rep(features[selected], times = ncol(mat)),
      component = rep(seq_len(ncol(mat)), each = nrow(mat)), loading = as.vector(mat))
    measure <- "Original SuperPC singular-vector loadings"
    interpretation <- "Selected features only; native SVD signs retained. Not coefficients, hazard ratios or importance.\nExcluded input features remain in the feature-screening table, not fabricated loading rows."
  } else {
    if (id != "superpc_cv") sb_abort("Fixed SuperPC has no CV threshold search.",
      "stormbreaker_unavailable_diagnostic")
    cv <- engine$cross_validation
    native <- list(thresholds = cv$thresholds, scores = cv$scor, lower = cv$scor.lower, upper = cv$scor.upper)
    if (!is.numeric(cv$thresholds) || !length(cv$thresholds) || any(!is.finite(cv$thresholds)) ||
        !is.matrix(cv$scor) || !nrow(cv$scor) || ncol(cv$scor) != length(cv$thresholds) ||
        !identical(dim(cv$scor), dim(cv$scor.lower)) || !identical(dim(cv$scor), dim(cv$scor.upper))) sb_abort(
      "Original SuperPC CV dimensions are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(threshold = rep(cv$thresholds, each = nrow(cv$scor)),
      component = rep(seq_len(nrow(cv$scor)), times = ncol(cv$scor)),
      score = as.vector(cv$scor), lower = as.vector(cv$scor.lower), upper = as.vector(cv$scor.upper))
    tab$selected <- tab$component == model$params$prediction_components &
      rep(seq_along(cv$thresholds) == model$audit$selected_threshold_index, each = nrow(cv$scor))
    if (sum(tab$selected) != 1L || !is.finite(tab$score[tab$selected]) ||
        !identical(as.numeric(tab$threshold[tab$selected]), as.numeric(threshold)) ||
        any(vapply(tab[c("score", "lower", "upper")], function(z) any(is.infinite(z)), logical(1)))) sb_abort(
      "Original SuperPC CV scores or selected threshold are invalid.", "stormbreaker_invalid_diagnostic")
    folds <- c(requested = model$audit$requested_folds, effective = cv$n.fold)
    measure <- "Original SuperPC threshold cross-validation"
    interpretation <- paste0("Geometric mean held-out likelihood-ratio statistic; bands are +/- one log-scale SE, not CI.\n",
      "Compare thresholds within a component count, not across counts with different degrees of freedom.\n",
      "Original one-component threshold selection retained; folds requested/effective: ", folds[1], "/", folds[2], ".")
  }
  structure(list(type = paste0("superpc_", type), table = tab, native_table = native,
    learner = id, pipeline_id = fit$pipeline$id, threshold = threshold, folds = folds,
    measure = measure, interpretation = interpretation, source = "original captured SuperPC fit; no refitting"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_superpc_diagnostic <- function(x) {
  feature <- score <- selected <- loading <- component <- threshold <- lower <- upper <- NULL
  tab <- x$table
  if (x$type == "superpc_feature_scores") {
    tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$score)])
    p <- ggplot2::ggplot(tab, ggplot2::aes(score, feature, fill = selected)) +
      ggplot2::geom_col(width = .65) +
      ggplot2::geom_vline(xintercept = c(-x$threshold, x$threshold), linetype = 2, colour = "#D77935") +
      ggplot2::scale_fill_manual(values = c("FALSE" = "#B8B8B8", "TRUE" = "#2878A0")) +
      ggplot2::labs(x = "Native screening score", y = NULL, fill = "Selected")
  } else if (x$type == "superpc_loadings") {
    tab$component <- factor(tab$component)
    p <- ggplot2::ggplot(tab, ggplot2::aes(loading, feature, fill = component)) +
      ggplot2::geom_col(position = "dodge", width = .65) +
      ggplot2::labs(x = "Native singular-vector loading", y = NULL, fill = "Component")
  } else {
    selected_component <- as.character(tab$component[tab$selected])
    tab$component <- factor(tab$component)
    tab$panel <- ifelse(as.character(tab$component) == selected_component,
      paste("Component", tab$component, "(used)"), paste("Component", tab$component))
    p <- ggplot2::ggplot(tab, ggplot2::aes(threshold, score)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#2878A0", alpha = .12, colour = NA, na.rm = TRUE) +
      ggplot2::geom_line(colour = "#2878A0", na.rm = TRUE) +
      ggplot2::geom_point(data = tab[tab$selected, ], shape = 21, colour = "black", size = 3) +
      ggplot2::geom_vline(data = tab[tab$selected, ], ggplot2::aes(xintercept = threshold), colour = "#D77935", linetype = 2) +
      ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
      ggplot2::labs(x = "Feature-screening threshold", y = "CV likelihood-ratio statistic")
  }
  p + ggplot2::labs(title = x$measure, subtitle = x$learner, caption = x$interpretation) + sb_viz_theme()
}

# ---- model-mixomics ----

sb_mixomics_da_learner_ids <- function() {
  c("plsda_binary", "splsda_binary")
}

sb_mixomics_da_params <- function(learner) {
  if (!learner %in% sb_mixomics_da_learner_ids()) {
    sb_abort(
      paste0("Unknown mixOmics discriminant learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    method = if (identical(learner, "plsda_binary")) {
      "PLS-DA"
    } else {
      "sPLS-DA"
    },
    requested_components = 2L,
    component_rule = "min_requested_features_analysis_rows_minus_one",
    scale = TRUE,
    tolerance = 1e-6,
    maximum_iterations = 100L,
    near_zero_variance_filter = FALSE,
    logratio = "none",
    keep_x_rule = if (identical(learner, "plsda_binary")) {
      "all_features"
    } else {
      "min_10_and_p_per_component"
    },
    probability_mapping = "binary_softmax_of_native_dummy_scores"
  )
}

sb_require_mixomics <- function() {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    sb_abort(
      "Package `mixOmics` is required for PLS-DA learners.",
      "stormbreaker_missing_mixomics"
    )
  }
  invisible(TRUE)
}

sb_mixomics_da_active_spec <- function(x, params) {
  components <- min(
    params$requested_components,
    ncol(x),
    nrow(x) - 1L
  )
  list(
    components = as.integer(components),
    keep_x = if (identical(params$method, "sPLS-DA")) {
      rep.int(as.integer(min(10L, ncol(x))), components)
    } else {
      NULL
    }
  )
}

sb_mixomics_da_native_fit <- function(x, grouping, params, active, learner) {
  fit <- function() {
    arguments <- list(
      X = data.matrix(x),
      Y = grouping,
      ncomp = active$components,
      scale = params$scale,
      tol = params$tolerance,
      max.iter = params$maximum_iterations,
      near.zero.var = params$near_zero_variance_filter,
      logratio = params$logratio,
      all.outputs = TRUE
    )
    if (identical(params$method, "sPLS-DA")) {
      arguments$keepX <- active$keep_x
      return(do.call(mixOmics::splsda, arguments))
    }
    do.call(mixOmics::plsda, arguments)
  }

  tryCatch(
    withCallingHandlers(
      fit(),
      warning = function(warning) {
        sb_abort(
          paste0(
            "Learner `", learner, "` raised a fitting warning: ",
            conditionMessage(warning)
          ),
          "stormbreaker_mixomics_da_fit_warning",
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
          "Learner `", learner, "` could not fit its mixOmics model: ",
          conditionMessage(error)
        ),
        "stormbreaker_mixomics_da_fit_error",
        parent = error
      )
    }
  )
}

sb_mixomics_da_projection <- function(native, features, components, learner) {
  x <- native$X
  y <- native$ind.mat
  variates <- native$variates$X
  loadings <- native$loadings$X
  index <- seq_len(components)
  p_matrix <- crossprod(x, variates)
  c_matrix <- crossprod(y, variates)
  coefficients <- loadings[, index, drop = FALSE] %*%
    solve(
      t(p_matrix[, index, drop = FALSE]) %*%
        loadings[, index, drop = FALSE]
    ) %*%
    t(c_matrix)[index, , drop = FALSE]
  rownames(coefficients) <- features

  engine <- list(
    coefficients = coefficients,
    x_center = unname(attr(x, "scaled:center")),
    x_scale = unname(attr(x, "scaled:scale")),
    y_center = unname(attr(y, "scaled:center")),
    y_scale = unname(attr(y, "scaled:scale"))
  )
  valid <- identical(dim(engine$coefficients), c(length(features), 2L)) &&
    all(vapply(engine, is.numeric, logical(1))) &&
    all(vapply(engine, function(value) {
      !anyNA(value) && all(is.finite(value))
    }, logical(1))) &&
    length(engine$x_center) == length(features) &&
    length(engine$x_scale) == length(features) &&
    length(engine$y_center) == 2L && length(engine$y_scale) == 2L &&
    all(engine$x_scale > 0) && all(engine$y_scale > 0)
  if (!valid) {
    sb_abort(
      paste0("Learner `", learner, "` returned an incomplete projection."),
      "stormbreaker_incomplete_mixomics_da_fit"
    )
  }
  engine
}

sb_fit_mixomics_da <- function(x, outcome, learner, params) {
  sb_require_mixomics()
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      paste0("A `", learner, "` analysis fold must contain both classes."),
      "stormbreaker_mixomics_da_single_class"
    )
  }

  active <- sb_mixomics_da_active_spec(x, params)
  native <- sb_mixomics_da_native_fit(
    x,
    factor(outcome, levels = c(0L, 1L)),
    params,
    active,
    learner
  )
  engine <- tryCatch(
    sb_mixomics_da_projection(native, names(x), active$components, learner),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "Learner `", learner, "` could not form its prediction state: ",
          conditionMessage(error)
        ),
        "stormbreaker_mixomics_da_projection_error",
        parent = error
      )
    }
  )
  component_features <- if (identical(params$method, "sPLS-DA")) {
    lapply(seq_len(active$components), function(component) {
      names(x)[native$loadings$X[, component] != 0]
    })
  } else {
    list()
  }
  active_params <- params
  active_params$effective_components <- active$components
  active_params$active_keep_x <- active$keep_x

  structure(
    list(
      learner = learner,
      features = names(x),
      params = active_params,
      engine = engine,
      diagnostic_components = list(loadings = native$loadings$X,
        scores = native$variates$X, loading_definition = "Native mixOmics X component weights"),
      audit = list(
        effective_components = active$components,
        active_keep_x = active$keep_x,
        native_iterations = as.integer(native$iter),
        component_features = component_features,
        probability_mapping = params$probability_mapping
      )
    ),
    class = "stormbreaker_mixomics_da"
  )
}

sb_predict_mixomics_da <- function(model, x) {
  if (!inherits(model, "stormbreaker_mixomics_da")) {
    sb_abort(
      "The fitted mixOmics discriminant model is invalid.",
      "stormbreaker_bad_mixomics_da_model"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted mixOmics model.",
      "stormbreaker_mixomics_da_schema_mismatch"
    )
  }

  matrix_x <- data.matrix(x)
  matrix_x <- sweep(matrix_x, 2L, model$engine$x_center, "-")
  matrix_x <- sweep(matrix_x, 2L, model$engine$x_scale, "/")
  scores <- matrix_x %*% model$engine$coefficients
  scores <- sweep(scores, 2L, model$engine$y_scale, "*")
  scores <- sweep(scores, 2L, model$engine$y_center, "+")
  probability <- stats::plogis(scores[, "1"] - scores[, "0"])
  probability <- unname(as.numeric(probability))
  if (length(probability) != nrow(x) || anyNA(probability) ||
      !all(is.finite(probability)) ||
      !all(probability >= 0 & probability <= 1)) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid probabilities."),
      "stormbreaker_invalid_mixomics_da_prediction"
    )
  }
  probability
}

sb_extract_mixomics_da_learner_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_mixomics_da")) {
    sb_bad_learner_state("The fitted mixOmics discriminant state is invalid.")
  }
  portable <- unclass(model)
  portable$diagnostic_components <- NULL
  structure(
    portable,
    class = c(
      "stormbreaker_mixomics_da_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_mixomics_da_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
      !learner %in% sb_mixomics_da_learner_ids()) {
    sb_bad_learner_state("The saved mixOmics discriminant state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- structure(unclass(state), class = "stormbreaker_mixomics_da")
  list(probability = sb_predict_mixomics_da(model, x))
}

sb_mixomics_da_learner_adapter <- function(
    learner,
    params = sb_mixomics_da_params(learner)) {
  if (!learner %in% sb_mixomics_da_learner_ids()) {
    sb_abort(
      paste0("Unknown mixOmics discriminant learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_mixomics_da(x, outcome, learner, params)
    },
    predict = sb_predict_mixomics_da,
    extract_state = sb_extract_mixomics_da_learner_state,
    predict_state = sb_predict_mixomics_da_learner_state
  )
}

# ---- model-pamr ----

# ---- algorithm-pamr ----

#' @rdname sb_pls_diagnostics
#' @export
sb_pamr_diagnostics <- function(fit) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% c("pamr_binary", "pamr_survival")) sb_abort(
    "This diagnostic requires PAMR.", "stormbreaker_unsupported_diagnostic")
  sb_require_pamr_survival_dependency()
  e <- if (s$learner == "pamr_survival") s$engine else unserialize(s$engine_raw)
  p <- length(s$features)
  original <- e$centroids
  if (!is.matrix(original) || nrow(original) != p || ncol(original) != 2L ||
      !identical(rownames(original), s$features) ||
      !identical(colnames(original), if (s$learner == "pamr_survival") c("1", "2") else c("0", "1")) ||
      any(!is.finite(original)) || length(e$centroid.overall) != p ||
      length(e$sd) != p || any(!is.finite(e$sd)) || any(e$sd <= 0)) sb_abort(
    "Saved PAMR centroids or scaling are invalid.", "stormbreaker_bad_learner_state")
  dummy <- matrix(0, nrow = p, ncol = 1L)
  shrunken <- pamr::pamr.predict(e, newx = dummy, threshold = s$params$threshold, type = "centroid")
  if (!identical(dim(shrunken), dim(original)) || any(!is.finite(shrunken))) sb_abort(
    "Native shrunken centroids are invalid.", "stormbreaker_invalid_diagnostic")
  classes <- colnames(original)
  labels <- if (s$learner == "pamr_survival") c("Early", "Late") else classes
  make <- function(z, stage) data.frame(feature = rep(s$features, 2L),
    class = rep(classes, each = p), column = rep(paste(stage, labels), each = p),
    stage = stage, original_value = as.numeric(z),
    value = as.numeric((z - e$centroid.overall) / e$sd))
  tab <- rbind(make(original, "Before"), make(shrunken, "After"))
  structure(list(type = "pamr_centroids", table = tab, features = s$features,
    columns = unique(tab$column), threshold = s$params$threshold,
    overall = e$centroid.overall, sd = e$sd, learner = s$learner,
    pipeline_id = fit$pipeline$id, source = "saved native PAMR centroids and pamr.predict centroid transformation",
    interpretation = paste("All features retained before and after the original fixed threshold; no retuning.",
      "Colour is deviation from overall centroid divided by saved PAMR SD; original-unit values remain in the table.",
      "Survival columns denote native early/late survival groups, not patient subgroup analysis or survival probabilities.",
      "No training history, confidence intervals or predictive performance is inferred.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_pamr_diagnostic <- function(x) {
  column <- feature <- value <- NULL
  d <- x$table
  d$column <- factor(d$column, levels = x$columns)
  d$feature <- factor(d$feature, levels = rev(x$features))
  p <- ggplot2::ggplot(d, ggplot2::aes(column, feature, fill = value)) +
    ggplot2::geom_tile(colour = "#D9DFE3", linewidth = .2) +
    ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#C96942", midpoint = 0,
      guide = ggplot2::guide_colourbar(display = "rectangles")) +
    ggplot2::labs(title = "PAMR centroid shrinkage", x = NULL, y = NULL, fill = "PAMR-scaled\ndeviation") + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = max(5, 2 + .3 * length(x$features)))
  p
}

# ---- pamr-survival-adapter ----

sb_pamr_survival_learner_ids <- function() {
  "pamr_survival"
}

sb_pamr_survival_params <- function(
  learner = sb_pamr_survival_learner_ids()
) {
  if (!identical(learner, sb_pamr_survival_learner_ids())) {
    sb_abort(
      paste0("Unknown PAMR survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    threshold = 1,
    scale_sd = TRUE,
    offset_percent = 50,
    remove_zero_thresholds = FALSE,
    prior = "native_equal_survival_groups",
    survival_groups = 2L,
    ties = "breslow"
  )
}

sb_pamr_survival_parameter_semantics <- function() {
  paste0(
    "two native censoring-aware survival groups, nearest shrunken-centroid ",
    "threshold 1, variance scaling, equal group priors, and offset 50 fixed ",
    "before resampling"
  )
}

sb_pamr_survival_available <- function() {
  requireNamespace("pamr", quietly = TRUE)
}

sb_require_pamr_survival_dependency <- function() {
  if (!sb_pamr_survival_available()) {
    sb_abort(
      "Package `pamr` is required for learner `pamr_survival`.",
      "stormbreaker_missing_pamr_survival_dependency"
    )
  }
  invisible(TRUE)
}

sb_new_pamr_survival_model <- function(features, params, audit, engine) {
  structure(
    list(
      learner = sb_pamr_survival_learner_ids(),
      features = features,
      params = params,
      audit = audit,
      engine = engine
    ),
    class = "stormbreaker_pamr_survival_model"
  )
}

sb_fit_pamr_survival <- function(x, time, event, learner, params) {
  inputs <- sb_validate_classical_survival_inputs(x, time, event, learner)
  sb_require_pamr_survival_dependency()
  if (ncol(inputs$x) < 2L) {
    sb_abort(
      paste0(
        "PAMR survival training requires at least 2 predictors because ",
        "the native subset path does not preserve a one-row feature matrix."
      ),
      "stormbreaker_pamr_survival_too_few_features"
    )
  }
  data <- list(
    x = t(data.matrix(inputs$x)),
    survival.time = inputs$time,
    censoring.status = inputs$event,
    geneid = names(inputs$x),
    genenames = names(inputs$x)
  )
  engine <- NULL
  invisible(utils::capture.output(
    engine <- pamr::pamr.train(
      data,
      threshold = params$threshold,
      scale.sd = params$scale_sd,
      offset.percent = params$offset_percent,
      prior = NULL,
      remove.zeros = params$remove_zero_thresholds,
      sign.contrast = "both",
      ngroup.survival = params$survival_groups
    ),
    type = "output"
  ))
  cutoff <- unname(as.numeric(engine$cutoffs.survival))
  valid <- inherits(engine, "pamrtrained") &&
    identical(engine$problem.type, "surv.km") &&
    length(cutoff) == 1L && is.finite(cutoff) && cutoff > 0 &&
    identical(as.integer(engine$ngroup.survival), params$survival_groups)
  if (!valid) {
    sb_abort(
      "PAMR did not return a two-group survival centroid fit.",
      "stormbreaker_incomplete_pamr_survival_fit"
    )
  }
  model <- sb_new_pamr_survival_model(
    names(inputs$x),
    params,
    list(
      survival_cutoff = cutoff,
      early_survival_class = "1",
      late_survival_class = "2"
    ),
    sb_pamr_binary_prediction_engine(engine)
  )
  sb_predict_pamr_survival_risk(model, inputs$x)
  model
}

sb_validate_pamr_survival_model_x <- function(model, x) {
  if (!inherits(model, "stormbreaker_pamr_survival_model")) {
    sb_abort(
      "The supplied model is not a fitted PAMR survival learner.",
      "stormbreaker_bad_pamr_survival_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted PAMR survival schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  list(model = model, x = x)
}

sb_predict_pamr_survival_risk <- function(model, x) {
  sb_require_pamr_survival_dependency()
  validated <- sb_validate_pamr_survival_model_x(model, x)
  posterior <- pamr::pamr.predict(
    validated$model$engine,
    newx = t(data.matrix(validated$x)),
    threshold = validated$model$params$threshold,
    type = "posterior"
  )
  risk <- if (is.matrix(posterior) && "1" %in% colnames(posterior)) {
    unname(posterior[, "1"])
  } else {
    numeric()
  }
  if (length(risk) != nrow(validated$x) || anyNA(risk) ||
    any(!is.finite(risk)) || any(risk < 0 | risk > 1)) {
    sb_abort(
      "PAMR produced invalid early-survival posterior risk scores.",
      "stormbreaker_invalid_pamr_survival_prediction"
    )
  }
  risk
}

sb_pamr_survival_source <- function() {
  "pamr_survival_analysis_cox_breslow_calibration"
}

sb_fit_pamr_survival_calibration <- function(
  time, event, risk, prediction_horizon
) {
  calibration <- sb_classical_survival_stop_on_warning(
    survival::coxph(
      survival::Surv(time, event) ~ risk,
      data = data.frame(time = time, event = event, risk = risk),
      ties = "breslow",
      singular.ok = FALSE,
      model = FALSE,
      x = FALSE,
      y = FALSE
    ),
    "PAMR survival Cox calibration",
    "stormbreaker_pamr_survival_calibration_warning"
  )
  coefficient <- unname(as.numeric(stats::coef(calibration)))
  if (length(coefficient) != 1L || is.na(coefficient) ||
    !is.finite(coefficient) || coefficient < 0) {
    sb_abort(
      "PAMR survival risk could not be calibrated in its native direction.",
      "stormbreaker_invalid_pamr_survival_calibration"
    )
  }
  list(
    coefficient = coefficient,
    cumulative_baseline_hazard = sb_breslow_baseline_hazard(
      time,
      event,
      coefficient * risk,
      prediction_horizon
    )
  )
}

sb_predict_pamr_survival_probability <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted PAMR survival schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_risk <- sb_predict_pamr_survival_risk(model, inputs$x)
  assessment_risk <- sb_predict_pamr_survival_risk(model, x_assessment)
  calibration <- sb_fit_pamr_survival_calibration(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  result <- sb_survival_prediction_result(
    prediction_horizon,
    calibration$cumulative_baseline_hazard,
    source = sb_pamr_survival_source(),
    ties = model$params$ties,
    risk = calibration$coefficient * assessment_risk
  )
  result$calibration_coefficient <- calibration$coefficient
  result
}

sb_extract_pamr_survival_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_pamr_survival_model")) {
    sb_bad_learner_state("The fitted PAMR survival learner state is invalid.")
  }
  baseline_hazard <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_pamr_survival_probability(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    baseline_hazard <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      cumulative_baseline_hazard = probability$cumulative_baseline_hazard,
      baseline_survival = probability$baseline_survival,
      source = probability$source,
      ties = probability$ties,
      calibration_coefficient = probability$calibration_coefficient
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      audit = model$audit,
      engine = model$engine,
      training_data_retention = "none",
      risk_reference = "pamr_native_survival_groups",
      risk_formula = "native_posterior_for_early_survival_group",
      baseline_hazard = baseline_hazard
    ),
    class = c(
      "stormbreaker_pamr_survival_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_pamr_survival_state <- function(state, x) {
  if (!inherits(state, "stormbreaker_pamr_survival_state") ||
    !identical(state$learner, sb_pamr_survival_learner_ids()) ||
    !inherits(state$engine, "pamrtrained") ||
    !identical(state$training_data_retention, "none")) {
    sb_bad_learner_state("The saved PAMR survival learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, state$learner)
  model <- sb_new_pamr_survival_model(
    state$features,
    state$params,
    state$audit,
    state$engine
  )
  risk <- sb_predict_pamr_survival_risk(model, x)
  result <- list(risk = risk)
  baseline <- state$baseline_hazard
  if (!is.null(baseline)) {
    coefficient <- baseline$calibration_coefficient
    valid <- is.numeric(coefficient) && length(coefficient) == 1L &&
      is.finite(coefficient) && coefficient >= 0 &&
      identical(baseline$source, sb_pamr_survival_source()) &&
      identical(baseline$ties, state$params$ties)
    if (!valid) {
      sb_bad_learner_state("The saved PAMR survival horizon state is invalid.")
    }
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      coefficient * risk,
      baseline$cumulative_baseline_hazard
    )
  }
  result
}

sb_pamr_survival_learner_adapter <- function(
  learner = sb_pamr_survival_learner_ids(),
  params = sb_pamr_survival_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_pamr_survival(x, time, event, learner, params)
    },
    predict = sb_predict_pamr_survival_risk,
    predict_survival = sb_predict_pamr_survival_probability,
    extract_state = sb_extract_pamr_survival_state,
    predict_state = sb_predict_pamr_survival_state
  )
}

# ---- diagnostics-projections ----

# ---- algorithm-pls-cv ----

sb_pls_cv_diagnostic <- function(fit) {
  if (fit$pipeline$learner != "plsrcox_cv") sb_abort(
    "This PLS learner has no component-count CV search.", "stormbreaker_unavailable_diagnostic")
  audit <- fit$learner_state$audit
  if (identical(audit$component_selection, "single_predictor_has_one_feasible_component")) sb_abort(
    "One predictor has one feasible component; no CV search was performed.", "stormbreaker_unavailable_diagnostic")
  value <- audit$cv_error
  spread <- audit$cv_standard_error
  maximum <- audit$maximum_components
  selected <- audit$selected_components
  if (!is.numeric(maximum) || length(maximum) != 1L || !is.finite(maximum) || maximum < 1L ||
      maximum != floor(maximum) || !is.numeric(value) || length(value) != maximum + 1L ||
      !is.numeric(spread) || length(spread) != length(value) || any(is.infinite(value)) ||
      any(is.infinite(spread)) || any(spread < 0, na.rm = TRUE) ||
      !is.numeric(selected) || length(selected) != 1L || !is.finite(selected) ||
      selected < 1L || selected > maximum || selected != floor(selected) ||
      !is.finite(value[selected + 1L])) sb_abort(
    "Saved plsRcox CV criterion or component count is invalid.", "stormbreaker_invalid_diagnostic")
  tab <- data.frame(component = seq.int(0L, maximum), iauc = value,
    native_spread = spread, selected = seq.int(0L, maximum) == selected)
  structure(list(type = "pls_cv_curve", table = tab,
    native_table = list(cv_error = value, cv_standard_error = spread),
    selected_components = as.integer(selected), folds = length(audit$folds),
    search_warnings = audit$component_search_warnings,
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    interpretation = paste0("Original component search, including zero-component baseline; not whole-pipeline validation.\n",
      "Bars: saved native cv.se5 spread (fold SD / fold count), not conventional SE or confidence intervals.\n",
      "Missing values remain missing; fold-level completion cannot be recovered from historical summaries."),
    source = "saved native plsRcox criterion 5 (iAUC_sh); no fitting"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_pls_cv_diagnostic <- function(x) {
  component <- iauc <- native_spread <- NULL
  tab <- x$table
  p <- ggplot2::ggplot(tab, ggplot2::aes(component, iauc)) +
    ggplot2::geom_line(colour = "#2878A0", na.rm = TRUE) +
    ggplot2::geom_point(colour = "#2878A0", na.rm = TRUE) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = iauc - native_spread,
      ymax = iauc + native_spread), width = .15, colour = "#2878A0", na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = x$selected_components, colour = "#D77935", linetype = 2) +
    ggplot2::scale_x_continuous(breaks = function(limits) unique(round(pretty(limits)))) +
    ggplot2::labs(title = "PLS-Cox component-count selection",
      subtitle = paste0(x$folds, " original folds; selected components: ", x$selected_components,
        if (length(x$search_warnings)) paste0("; retained search warnings: ", length(x$search_warnings)) else ""),
      x = "Number of components", y = "Integrated AUC (iAUC_sh; higher is better)",
      caption = x$interpretation) + sb_viz_theme()
  p
}

# ---- algorithm-pls ----

#' @rdname sb_pls_diagnostics
#' @export
sb_pls_diagnostics <- function(fit, type = "loadings") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  id <- fit$pipeline$learner
  if (!id %in% c("plsrglm_binary", "plsda_binary", "splsda_binary",
      "plsrcox", "plsrcox_cv")) sb_abort("This diagnostic requires a supported native PLS learner.",
    "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("loadings", "score_distribution", "cv_curve")) sb_abort(
    "Choose loadings, score_distribution or cv_curve.", "stormbreaker_unsupported_diagnostic")
  if (type == "cv_curve") return(sb_pls_cv_diagnostic(fit))
  model <- fit$diagnostic_artifacts$learner
  if (is.null(model)) sb_abort("Original PLS components were not captured; enable capture_diagnostics at final fitting.",
    "stormbreaker_missing_diagnostic_artifacts")
  native <- if (id %in% c("plsrcox", "plsrcox_cv")) list(
    loadings = model$engine$pp, scores = model$engine$tt,
    loading_definition = "Native plsRcox X reconstruction loadings") else model$diagnostic_components
  if (is.null(native)) sb_abort("This historical model did not retain native PLS components.",
    "stormbreaker_missing_diagnostic_artifacts")
  mat <- if (type == "loadings") native$loadings else native$scores
  expected <- if (type == "loadings") length(model$features) else nrow(fit$diagnostic_artifacts$predictors)
  if (!is.matrix(mat) || !is.numeric(mat) || nrow(mat) != expected || ncol(mat) < 1L ||
      any(!is.finite(mat))) sb_abort("Native PLS component matrix is invalid.", "stormbreaker_invalid_diagnostic")
  if (type == "loadings" && !identical(rownames(mat), model$features)) sb_abort(
    "Native PLS loading order does not match the fitted features.", "stormbreaker_invalid_diagnostic")
  tab <- if (type == "loadings") data.frame(feature = rep(model$features, times = ncol(mat)),
    component = rep(seq_len(ncol(mat)), each = nrow(mat)), loading = as.vector(mat)) else
    data.frame(observation = rep(seq_len(nrow(mat)), times = ncol(mat)),
      component = rep(seq_len(ncol(mat)), each = nrow(mat)), score = as.vector(mat))
  structure(list(type = paste0("pls_", type), table = tab,
    native_table = if (type == "loadings") mat else unname(mat),
    learner = id, pipeline_id = fit$pipeline$id, features = model$features,
    loading_definition = native$loading_definition,
    sensitive_training_data = type == "score_distribution",
    interpretation = "Original development-fit components; native signs and component order retained.\nNot feature importance, hazard ratios, probability effects or evidence of held-out class separation.",
    source = "original captured PLS components; no fitting or decomposition"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_pls_diagnostic <- function(x) {
  if (x$type == "pls_cv_curve") return(sb_plot_pls_cv_diagnostic(x))
  feature <- component <- loading <- score <- NULL
  tab <- x$table
  tab$component <- factor(tab$component, levels = sort(unique(tab$component)))
  if (x$type == "pls_loadings") {
    tab$feature <- factor(tab$feature, levels = rev(x$features))
    extent <- max(abs(tab$loading))
    if (extent == 0) extent <- 1
    p <- ggplot2::ggplot(tab, ggplot2::aes(component, feature, fill = loading)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#D77935",
        midpoint = 0, limits = c(-extent, extent),
        guide = ggplot2::guide_colourbar(display = "rectangles")) +
      ggplot2::labs(title = x$loading_definition, x = "Original component", y = NULL, fill = "Loading")
  } else {
    p <- ggplot2::ggplot(tab, ggplot2::aes(component, score)) +
      ggplot2::geom_boxplot(fill = "#B5D3E1", outlier.shape = NA, width = .45) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = .12, height = 0, seed = 1L),
        colour = "#2878A0", size = .8, alpha = .5) +
      ggplot2::labs(title = "Development-sample component scores",
        x = "Original component", y = "Native component score")
  }
  caption <- x$interpretation
  if (x$type == "pls_score_distribution") caption <- paste(caption,
    "Box: median and interquartile range; whiskers: observations within 1.5 IQR.\nAll observations shown; boxes and whiskers are not confidence intervals.", sep = "\n")
  p + ggplot2::labs(subtitle = x$learner, caption = caption) + sb_viz_theme()
}
