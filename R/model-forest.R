#' Tree and forest diagnostics
#' @name sb_forest_diagnostics
#' @description Inspect native tree structures, forest importance and saved ensemble diagnostics.
#' @rawRd \arguments{
#' \item{fit}{
#' \strong{For \code{sb_forest_diagnostics()}:} A final ranger, randomForest/RRF binary or randomForestSRC survival
#' model. Importance and tree sizes use saved prediction state. OOB error
#' requires capture_diagnostics = TRUE for randomForest/randomForestSRC.
#' 
#' \strong{For \code{sb_aorsf_diagnostics()}:} A successful aorsf binary or survival final fit.
#' 
#' \strong{For \code{sb_bagging_diagnostics()}:} A successful ipred bagging final fit.
#' 
#' \strong{For \code{sb_c50_diagnostics()}:} Successful C5.0 final fit with original diagnostic capture.
#' 
#' \strong{For \code{sb_cforest_diagnostics()}:} A successful binary or survival cforest final fit.
#' 
#' \strong{For \code{sb_classification_tree_diagnostics()}:} Successful tree_binary final fit.
#' 
#' \strong{For \code{sb_ctree_structure()}:} A successful binary or survival ctree final fit.
#' 
#' \strong{For \code{sb_evtree_structure()}:} Successful evolutionary tree final fit.
#' 
#' \strong{For \code{sb_ferns_diagnostics()}:} Successful random-fern final fit.
#' 
#' \strong{For \code{sb_rlt_diagnostics()}:} A successful binary or survival RLT final fit.
#' 
#' \strong{For \code{sb_rotation_forest_diagnostics()}:} Successful rotation-forest final fit.
#' 
#' \strong{For \code{sb_tree_structure()}:} A final rpart classification or survival tree, fitted with
#' capture_diagnostics = TRUE. Forests and boosted ensembles are not trees.
#' 
#' \strong{For \code{sb_wsrf_diagnostics()}:} A successful weighted-subspace forest final fit.
#' }
#' \item{type}{
#' \strong{For \code{sb_forest_diagnostics()}:} One of importance, tree_size or oob_error.
#' 
#' \strong{For \code{sb_aorsf_diagnostics()}:} tree_size or oob_score.
#' 
#' \strong{For \code{sb_bagging_diagnostics()}:} split_usage, tree_size or importance.
#' 
#' \strong{For \code{sb_c50_diagnostics()}:} usage, splits or size.
#' 
#' \strong{For \code{sb_cforest_diagnostics()}, \code{sb_rlt_diagnostics()}:} split_usage or tree_size.
#' 
#' \strong{For \code{sb_classification_tree_diagnostics()}:} structure or leaf_probabilities.
#' 
#' \strong{For \code{sb_ferns_diagnostics()}:} split_usage or thresholds.
#' 
#' \strong{For \code{sb_rotation_forest_diagnostics()}:} loadings, tree_size or tree_predictions.
#' 
#' \strong{For \code{sb_wsrf_diagnostics()}:} importance, tree_nodes or tree_oob. Importance and OOB require
#' original capture_diagnostics; compact OOB placeholders are never used.
#' }
#' \item{newdata}{
#' Explicit query rows required for tree_predictions.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_forest_diagnostics()}}{
#' Inspect native forest importance, tree sizes or saved OOB errors
#' 
#' \strong{Returns:} An algorithm diagnostic with the full native table and source labels.
#' Ranger importance is its saved impurity measure under the declared split
#' rule; randomForest importance is mean decrease in Gini. Neither is
#' permutation importance, a causal effect or directly comparable across
#' engines. All features, including zero/negative values, remain in the table.
#' 
#' Tree sizes count original terminal nodes, without extracting a surrogate
#' tree. OOB values are development-forest diagnostics, not validation of
#' preprocessing/feature selection or untouched-cohort performance. Missing
#' OOB values remain explicit. randomForestSRC may save only its final OOB
#' error: this produces one point, not an invented learning curve. Ranger OOB
#' computation was disabled in these pipelines and is not reconstructed.
#' No fitting, permutation, parameter selection or prediction occurs.
#' }
#' 
#' \subsection{\code{sb_aorsf_diagnostics()}}{
#' Inspect saved accelerated oblique forest diagnostics
#' 
#' \strong{Returns:} Original terminal-node counts, or saved OOB evaluation values on
#' the native evaluation schedule. Binary scores are AUC-ROC and survival
#' scores are Harrell's C-index, not error rates. A single final evaluation
#' stays a single point. Importance was disabled and is not reconstructed.
#' OOB evaluation does not validate prior preprocessing or feature selection.
#' Serialized R6 state is read into a separate object, without updating the fit.
#' }
#' 
#' \subsection{\code{sb_bagging_diagnostics()}}{
#' Inspect original bagged classification trees
#' 
#' \strong{Returns:} Original primary split counts, node/leaf counts, or complete
#' tree-by-feature native rpart importance. Absent feature importance entries
#' in a tree are zero, including stumps. The importance plot averages original
#' rpart scores across all bags; these may credit surrogate splits and are not
#' permutation importance or validation performance. OOB was disabled.
#' }
#' 
#' \subsection{\code{sb_c50_diagnostics()}}{
#' Inspect captured C5.0 tree or rule ensembles
#' 
#' \strong{Returns:} Native C5imp percentages or original trial sizes (leaves for trees,
#' rules for rule ensembles). Native unparseable importance is unavailable,
#' not replaced with fabricated zeros. Compact prediction state deliberately
#' omits these summaries, so original capture is required. No refitting occurs.
#' }
#' 
#' \subsection{\code{sb_cforest_diagnostics()}}{
#' Inspect original conditional inference forest structure
#' 
#' \strong{Returns:} All-feature split counts or all-tree node/leaf counts. Native
#' partynode child structure and original model-frame column mapping are
#' used; response columns are not mistaken for predictors. Counts describe
#' structure, not conditional permutation importance, statistical significance
#' or held-out performance. No new predictions, permutations or refits occur.
#' }
#' 
#' \subsection{\code{sb_classification_tree_diagnostics()}}{
#' Inspect the original tree-package classifier
#' 
#' \strong{Returns:} All native nodes and split text, or terminal node class probabilities
#' and original sample counts. Probabilities are training-leaf proportions,
#' not held-out calibration. Native split strings retain original precision;
#' the complete frame is returned, without surrogate fitting or pruning.
#' }
#' 
#' \subsection{\code{sb_ctree_structure()}}{
#' Inspect an original conditional inference tree
#' 
#' \strong{Returns:} An algorithm diagnostic containing every original node, branch rule,
#' training-node count, exact native split and available native node information.
#' Counts describe training membership, not validation performance. No new
#' hypothesis tests, probabilities or survival estimates are computed. The
#' compact saved engine suffices; no refitting or pruning is performed.
#' }
#' 
#' \subsection{\code{sb_evtree_structure()}}{
#' Inspect the original evolutionary classification tree
#' 
#' \strong{Returns:} All original nodes, split rules and retained training memberships.
#' The compact predictor schema maps native split indices to feature names;
#' fitted membership columns are never used as predictor labels. No new
#' evolutionary search, pruning, hypothesis test or surrogate tree is fitted.
#' }
#' 
#' \subsection{\code{sb_ferns_diagnostics()}}{
#' Inspect original random-fern splits
#' 
#' \strong{Returns:} Original split records and counts, without refitting. Split usage
#' describes randomly generated structure, not predictive importance.
#' Threshold plots currently require real-valued predictors.
#' }
#' 
#' \subsection{\code{sb_rlt_diagnostics()}}{
#' Inspect original reinforcement-learning forest structure
#' 
#' \strong{Returns:} Original axis-aligned split counts for every feature, or node and
#' terminal-node counts for every tree. Split counts are structural usage,
#' not permutation importance or predictive effects. The original pipeline
#' disabled importance computation; missing OOB or reinforcement histories
#' are not inferred. No model is retrained or replaced by a surrogate tree.
#' }
#' 
#' \subsection{\code{sb_rotation_forest_diagnostics()}}{
#' Inspect original rotation-forest geometry and predictions
#' 
#' \strong{Returns:} All original rotation matrices, tree sizes or per-tree probabilities.
#' Rotated axes are not the original features despite native reused names.
#' Loading signs are basis conventions, not effect directions or importance.
#' Per-tree spread is not a confidence interval. No new PCA, fitting or pruning.
#' }
#' 
#' \subsection{\code{sb_tree_structure()}}{
#' Inspect a fitted single rpart tree without replacing it
#' 
#' \strong{Returns:} An algorithm diagnostic with every node, parent, native split rule,
#' node count, node prediction, exact native split matrix and complexity table.
#' Binary node probabilities are training-node proportions, not validation
#' estimates. Survival node values are relative rates, not survival probabilities.
#' }
#' 
#' \subsection{\code{sb_wsrf_diagnostics()}}{
#' Inspect original weighted-subspace forest diagnostics
#' 
#' \strong{Returns:} Native mean-decrease-information-gain-ratio importance, total node
#' counts per tree, or individual-tree OOB errors. Individual-tree OOB is
#' not cumulative ensemble OOB or external performance. Artificial constant
#' bridge features, if present, remain separately identified in the table.
#' No permutation importance, learning curve or missing history is invented.
#' }
#' }
NULL

# ---- model-ranger ----

# ---- algorithm-forest ----

#' @rdname sb_forest_diagnostics
#' @export
sb_forest_diagnostics <- function(fit, type = "importance") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  learner <- fit$pipeline$learner
  ranger <- learner %in% sb_ranger_learner_ids()
  rf <- learner %in% c("random_forest_binary", "rrf_binary")
  src <- identical(learner, "rsf_survival")
  if (!ranger && !rf && !src) sb_abort("This forest engine is not yet supported.",
    "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("importance", "tree_size", "oob_error")) sb_abort(
    "Choose importance, tree_size or oob_error.", "stormbreaker_unsupported_diagnostic")
  state <- fit$learner_state
  if (type == "oob_error") {
    if (ranger) sb_abort(
      "Ranger OOB computation was disabled; no original OOB error sequence exists.",
      "stormbreaker_unavailable_diagnostic")
    if (is.null(fit$diagnostic_artifacts)) sb_abort(
      "Original OOB results require capture_diagnostics when fitting; they are not in prediction state.",
      "stormbreaker_missing_diagnostic_artifacts")
    engine <- fit$diagnostic_artifacts$learner$engine
    if (rf) {
      errors <- engine$err.rate
      if (!is.matrix(errors) || !identical(colnames(errors), c("OOB", "0", "1"))) sb_abort(
        "Native OOB error columns are not the expected overall/0/1 mapping.",
        "stormbreaker_invalid_diagnostic")
      tab <- data.frame(trees = rep(seq_len(nrow(errors)), ncol(errors)),
        series = rep(c("Overall", "Class 0", "Class 1"), each = nrow(errors)),
        error = as.vector(errors))
      measure <- "OOB misclassification rate"
    } else {
      errors <- engine$err.rate
      if (!is.numeric(errors) || !is.null(dim(errors))) sb_abort(
        "Expected a single-outcome native survival OOB error sequence.",
        "stormbreaker_invalid_diagnostic")
      tab <- data.frame(trees = seq_along(errors), series = "Overall", error = errors)
      measure <- "Native OOB error (1 - concordance)"
    }
    if (!nrow(tab) || max(tab$trees) != engine$ntree ||
        any(is.infinite(tab$error)) || any(tab$error < 0 | tab$error > 1, na.rm = TRUE)) sb_abort(
      "Native forest OOB errors are incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    tab$available <- is.finite(tab$error)
    source <- "original captured forest OOB results"
  } else if (type == "importance") {
    if (src) sb_abort(
      "This randomForestSRC fit disabled importance computation; none is invented.",
      "stormbreaker_unavailable_diagnostic")
    if (ranger) {
      tab <- state$feature_importance
      measure <- paste0("Native impurity importance (", state$params$splitrule, " splits)")
    } else {
      engine <- unserialize(state$engine_raw)
      importance <- engine$importance
      if (!is.matrix(importance) || !"MeanDecreaseGini" %in% colnames(importance)) sb_abort(
        "Native mean-decrease-Gini importance is missing.", "stormbreaker_invalid_diagnostic")
      tab <- data.frame(feature = rownames(importance), importance = importance[, "MeanDecreaseGini"])
      rownames(tab) <- NULL
      measure <- "Mean decrease in Gini"
    }
    if (!is.data.frame(tab) || !all(c("feature", "importance") %in% names(tab)) ||
        !is.numeric(tab$importance) || any(!is.finite(tab$importance)) ||
        anyDuplicated(tab$feature) || !setequal(tab$feature, state$features)) sb_abort(
      "Native forest importance is missing or invalid.", "stormbreaker_invalid_diagnostic")
    source <- "original forest importance retained in prediction state"
  } else {
    if (ranger) {
      forest <- unserialize(state$model_raw)
      leaves <- vapply(forest$child.nodeIDs, function(children) {
        if (length(children) != 2L || length(children[[1]]) != length(children[[2]])) sb_abort(
          "Native ranger child-node arrays are invalid.", "stormbreaker_invalid_diagnostic")
        sum(children[[1]] == 0 & children[[2]] == 0)
      }, integer(1))
      trees <- state$num_trees
    } else {
      engine <- unserialize(state$engine_raw)
      leaves <- if (rf) colSums(engine$forest$nodestatus == -1L) else engine$leaf.count
      trees <- engine$ntree
    }
    if (!is.numeric(leaves) || length(leaves) != trees || !length(leaves) ||
        any(!is.finite(leaves)) || any(leaves < 1) || any(leaves != floor(leaves))) sb_abort(
      "Native terminal-node counts are incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(tree = seq_along(leaves), terminal_nodes = as.integer(leaves))
    measure <- "Terminal nodes per original tree"
    source <- "original forest structure retained in prediction state"
  }
  structure(list(type = paste0("forest_", type), table = tab, measure = measure,
    learner = learner, pipeline_id = fit$pipeline$id, source = source),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_forest_diagnostic <- function(x) {
  importance <- feature <- terminal_nodes <- trees <- error <- series <- NULL
  if (x$type == "forest_importance") {
    tab <- x$table
    tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$importance)])
    p <- ggplot2::ggplot(tab, ggplot2::aes(importance, feature)) +
      ggplot2::geom_col(fill = "#2878A0", width = .65) +
      ggplot2::geom_vline(xintercept = 0, colour = "#909090") +
      ggplot2::labs(title = "Forest importance", subtitle = x$learner,
        x = x$measure, y = NULL,
        caption = "Native split-based importance, not permutation importance or a causal effect.\nCan favour variables with more possible splits; correlated features can share importance.") + sb_viz_theme()
    attr(p, "stormbreaker_suggested_size") <- c(width = 10, height = max(5, 2.5 + .28 * nrow(tab)))
    return(p)
  }
  if (x$type == "forest_tree_size") {
    counts <- as.data.frame(table(x$table$terminal_nodes), stringsAsFactors = FALSE)
    names(counts) <- c("terminal_nodes", "trees")
    counts$terminal_nodes <- as.integer(counts$terminal_nodes)
    return(ggplot2::ggplot(counts, ggplot2::aes(terminal_nodes, trees)) +
      ggplot2::geom_col(fill = "#2878A0", width = .8) +
      ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::scale_y_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::labs(title = "Forest tree-size distribution", subtitle = x$learner,
        x = x$measure, y = "Number of trees",
        caption = "Every fitted tree counted once. Original terminal-node counts; no pruning or surrogate.\nStructural complexity, not an estimate of predictive performance.") + sb_viz_theme())
  }
  tab <- x$table[x$table$available, , drop = FALSE]
  p <- ggplot2::ggplot(tab, ggplot2::aes(trees, error, colour = series)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_manual(values = c(Overall = "#2878A0", "Class 0" = "#909090", "Class 1" = "#D77935")) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_x_continuous(limits = c(0, max(x$table$trees)),
      breaks = function(z) unique(round(pretty(z)))) +
    ggplot2::labs(title = "Forest OOB results", subtitle = x$learner,
      x = "Number of accumulated trees", y = x$measure, colour = "Outcome",
      caption = "Only native recorded values shown; missing tree counts remain in the diagnostic table.\nOOB does not validate prior preprocessing/selection and is not external-cohort performance.") + sb_viz_theme()
  if (nrow(tab) > 1L) {
    # Keep missing tree positions in the line layer, so gaps are not interpolated.
    p <- p + ggplot2::geom_line(data = x$table, na.rm = TRUE, linewidth = .4)
  }
  if (!nrow(tab)) {
    p <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
      label = "No finite native OOB error was recorded.") + ggplot2::theme_void() +
      ggplot2::labs(title = "Forest OOB results", subtitle = x$learner)
  }
  p
}

# ---- ranger-binary-fixed-adapter ----

sb_ranger_binary_params <- function(
    num_trees = 500L,
    mtry = NULL,
    min_node_size = 10L,
    seed = 1L,
    learner = "ranger_binary") {
  if (!learner %in% sb_ranger_binary_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed ranger binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  mtry <- sb_validate_ranger_mtry(mtry)
  split_spec <- sb_ranger_splitrule_spec(learner)
  params <- list(
    num_trees = sb_check_positive_integer(
      num_trees,
      "ranger_num_trees",
      minimum = 1L
    ),
    mtry = mtry,
    mtry_rule = if (is.null(mtry)) "floor_sqrt_p" else "fixed",
    min_node_size = sb_check_positive_integer(
      min_node_size,
      "ranger_min_node_size",
      minimum = 1L
    ),
    seed = sb_check_positive_integer(
      seed,
      "ranger_seed",
      minimum = 1L
    ),
    probability = TRUE,
    splitrule = split_spec$splitrule,
    importance = "impurity",
    replace = split_spec$replace,
    sample_fraction = 1,
    max_depth = 0L,
    num_threads = 1L,
    oob_error = FALSE,
    write_forest = TRUE,
    na_action = "na.fail"
  )
  for (name in setdiff(names(split_spec), c("splitrule", "replace"))) {
    params[[name]] <- split_spec[[name]]
  }
  params
}

sb_fit_ranger_binary <- function(
  x, outcome, params, learner = "ranger_binary") {
  if (!learner %in% sb_ranger_binary_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed ranger binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "A fixed ranger binary analysis fold must contain both classes.",
      "stormbreaker_ranger_binary_single_class"
    )
  }
  if (!sb_ranger_available()) {
    sb_abort(
      "Package `ranger` is required for the fixed binary forest learner.",
      "stormbreaker_missing_ranger"
    )
  }

  effective_mtry <- if (is.null(params$mtry)) {
    max(1L, floor(sqrt(ncol(x))))
  } else {
    params$mtry
  }
  if (effective_mtry > ncol(x)) {
    sb_abort(
      "`ranger_mtry` cannot exceed the number of fitted features.",
      "stormbreaker_ranger_mtry_exceeds_features"
    )
  }

  fit_arguments <- list(
    x = as.matrix(x),
    y = factor(as.character(outcome), levels = c("0", "1")),
    num.trees = params$num_trees,
    mtry = effective_mtry,
    probability = params$probability,
    importance = params$importance,
    write.forest = params$write_forest,
    min.node.size = params$min_node_size,
    max.depth = params$max_depth,
    replace = params$replace,
    sample.fraction = params$sample_fraction,
    splitrule = params$splitrule,
    oob.error = params$oob_error,
    num.threads = params$num_threads,
    seed = params$seed,
    na.action = params$na_action,
    verbose = FALSE
  )
  if ("num_random_splits" %in% names(params)) {
    fit_arguments$num.random.splits <- params$num_random_splits
  }
  fitted <- withCallingHandlers(
    do.call(ranger::ranger, fit_arguments),
    warning = function(warning) {
      message <- conditionMessage(warning)
      if (!startsWith(message, "Unused arguments:") && !identical(message,
        "Argument 'num.random.splits' ignored if splitrule is not 'extratrees'.")) {
        return(invisible(NULL))
      }
      sb_abort(
        paste0(
          "The fixed ranger binary learner raised a warning: ",
          message
        ),
        "stormbreaker_ranger_binary_warning",
        parent = warning
      )
    }
  )
  valid_forest <- inherits(fitted$forest, "ranger.forest") &&
    identical(as.integer(fitted$num.trees), params$num_trees) &&
    identical(fitted$treetype, "Probability estimation") &&
    identical(fitted$forest$levels, c("0", "1"))
  if (!valid_forest) {
    sb_abort(
      "Fixed ranger binary did not return the declared probability forest.",
      "stormbreaker_incomplete_ranger_fit"
    )
  }
  importance <- fitted$variable.importance
  valid_importance <- is.numeric(importance) &&
    length(importance) == ncol(x) && !anyNA(importance) &&
    all(is.finite(importance)) &&
    setequal(names(importance), names(x))
  if (!valid_importance) {
    sb_abort(
      "Fixed ranger binary did not return one finite importance per feature.",
      "stormbreaker_invalid_ranger_importance"
    )
  }

  structure(
    list(
      learner = learner,
      features = names(x),
      params = params,
      num_trees = params$num_trees,
      effective_mtry = as.integer(effective_mtry),
      forest = fitted$forest,
      feature_importance = unname(importance[names(x)])
    ),
    class = "stormbreaker_ranger_binary"
  )
}

sb_predict_ranger_binary <- function(model, x) {
  if (!inherits(model, "stormbreaker_ranger_binary")) {
    sb_abort(
      "The supplied model is not a fitted fixed ranger binary learner.",
      "stormbreaker_bad_ranger_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed ranger binary schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (!sb_ranger_available()) {
    sb_abort(
      "Package `ranger` is required for fixed binary forest prediction.",
      "stormbreaker_missing_ranger"
    )
  }

  predictions <- stats::predict(
    model$forest,
    data = x,
    num.trees = model$num_trees,
    type = "response",
    seed = 0L,
    num.threads = 1L,
    verbose = FALSE
  )$predictions
  valid_matrix <- is.matrix(predictions) && nrow(predictions) == nrow(x) &&
    identical(colnames(predictions), c("0", "1")) &&
    !anyNA(predictions) && all(is.finite(predictions))
  if (!valid_matrix) {
    sb_abort(
      "Fixed ranger binary produced invalid class probability predictions.",
      "stormbreaker_invalid_ranger_prediction"
    )
  }
  probability <- unname(as.numeric(predictions[, "1"]))
  if (any(probability < 0 | probability > 1)) {
    sb_abort(
      "Fixed ranger binary produced invalid class-one probabilities.",
      "stormbreaker_invalid_ranger_prediction"
    )
  }
  probability
}

sb_extract_ranger_binary_tree_model <- function(model) {
  if (!inherits(model, "stormbreaker_ranger_binary")) {
    sb_abort(
      "The supplied model is not a fitted fixed ranger binary learner.",
      "stormbreaker_bad_ranger_fit"
    )
  }
  list(
    num_trees = model$num_trees,
    effective_mtry = model$effective_mtry,
    features = model$features,
    params = model$params,
    importance_type = model$params$importance,
    model_serialization = "R_xdr_v3_ranger.forest",
    model_raw = serialize(
      model$forest,
      connection = NULL,
      ascii = FALSE,
      xdr = TRUE,
      version = 3L
    ),
    feature_importance = data.frame(
      feature = model$features,
      importance = model$feature_importance,
      stringsAsFactors = FALSE
    )
  )
}

sb_ranger_binary_learner_adapter <- function(
    num_trees = 500L,
    mtry = NULL,
    min_node_size = 10L,
    seed = 1L,
    learner = "ranger_binary") {
  params <- sb_ranger_binary_params(
    num_trees = num_trees,
    mtry = mtry,
    min_node_size = min_node_size,
    seed = seed,
    learner = learner
  )
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_ranger_binary(x, outcome, params, learner)
    },
    predict = sb_predict_ranger_binary,
    extract_tree_model = sb_extract_ranger_binary_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_ranger_binary_learner_state(model)
    },
    predict_state = sb_predict_ranger_binary_learner_state
  )
}

# ---- ranger-survival-fixed-adapter ----

sb_ranger_binary_learner_ids <- function() {
  c(
    "ranger_binary",
    "ranger_extratrees_binary",
    "ranger_hellinger_binary"
  )
}

sb_ranger_survival_learner_ids <- function() {
  c(
    "ranger_survival",
    "ranger_extratrees_survival",
    "ranger_cindex_survival",
    "ranger_maxstat_survival"
  )
}

sb_ranger_learner_ids <- function() {
  c(sb_ranger_binary_learner_ids(), sb_ranger_survival_learner_ids())
}

sb_ranger_splitrule_spec <- function(learner) {
  switch(learner,
    ranger_binary = list(splitrule = "gini", replace = TRUE),
    ranger_extratrees_binary = list(
      splitrule = "extratrees",
      replace = FALSE,
      num_random_splits = 1L
    ),
    ranger_hellinger_binary = list(
      splitrule = "hellinger",
      replace = TRUE
    ),
    ranger_survival = list(splitrule = "logrank", replace = TRUE),
    ranger_extratrees_survival = list(
      splitrule = "extratrees",
      replace = FALSE,
      num_random_splits = 1L
    ),
    ranger_cindex_survival = list(splitrule = "C", replace = TRUE),
    ranger_maxstat_survival = list(
      splitrule = "maxstat",
      replace = TRUE,
      alpha = 0.5,
      minprop = 0.1
    ),
    sb_abort(
      paste0("Unknown fixed ranger learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_ranger_parameter_semantics <- function(learner) {
  switch(learner,
    ranger_binary = paste0(
      "Gini probability-forest trees, mtry rule, node size, bootstrap ",
      "sampling, and seed fixed before resampling"
    ),
    ranger_extratrees_binary = paste0(
      "Extremely randomized probability-forest split locations, one random ",
      "split per candidate variable, full no-replacement sampling, mtry rule, ",
      "node size, and seed fixed before resampling"
    ),
    ranger_hellinger_binary = paste0(
      "Hellinger-distance probability-forest trees, mtry rule, node size, ",
      "bootstrap sampling, and seed fixed before resampling"
    ),
    ranger_survival = paste0(
      "logrank survival-forest trees, mtry rule, node size, bootstrap ",
      "sampling, and seed fixed before resampling"
    ),
    ranger_extratrees_survival = paste0(
      "extremely randomized survival-forest split locations, one random ",
      "split per candidate variable, full no-replacement sampling, mtry rule, ",
      "node size, and seed fixed before resampling"
    ),
    ranger_cindex_survival = paste0(
      "C-index survival-forest split rule, mtry rule, node size, bootstrap ",
      "sampling, and seed fixed before resampling"
    ),
    ranger_maxstat_survival = paste0(
      "maximally selected rank-statistic survival-forest split rule, alpha ",
      "0.5, minprop 0.1, mtry rule, node size, bootstrap sampling, and seed ",
      "fixed before resampling"
    ),
    sb_abort(
      paste0("Unknown fixed ranger learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_ranger_available <- function() {
  requireNamespace("ranger", quietly = TRUE)
}

sb_validate_ranger_mtry <- function(mtry) {
  if (is.null(mtry)) {
    return(NULL)
  }
  sb_check_positive_integer(mtry, "ranger_mtry", minimum = 1L)
}

sb_ranger_survival_params <- function(
    num_trees = 500L,
    mtry = NULL,
    min_node_size = 3L,
    seed = 1L,
    learner = "ranger_survival") {
  if (!learner %in% sb_ranger_survival_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed ranger survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  mtry <- sb_validate_ranger_mtry(mtry)
  split_spec <- sb_ranger_splitrule_spec(learner)
  params <- list(
    num_trees = sb_check_positive_integer(
      num_trees,
      "ranger_num_trees",
      minimum = 1L
    ),
    mtry = mtry,
    mtry_rule = if (is.null(mtry)) "floor_sqrt_p" else "fixed",
    min_node_size = sb_check_positive_integer(
      min_node_size,
      "ranger_min_node_size",
      minimum = 1L
    ),
    seed = sb_check_positive_integer(
      seed,
      "ranger_seed",
      minimum = 1L
    ),
    splitrule = split_spec$splitrule,
    importance = "impurity",
    replace = split_spec$replace,
    sample_fraction = 1,
    max_depth = 0L,
    num_threads = 1L,
    oob_error = FALSE,
    write_forest = TRUE,
    na_action = "na.fail"
  )
  for (name in setdiff(names(split_spec), c("splitrule", "replace"))) {
    params[[name]] <- split_spec[[name]]
  }
  params
}

sb_fit_ranger_survival <- function(
  x, time, event, params, learner = "ranger_survival") {
  if (!learner %in% sb_ranger_survival_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed ranger survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
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
  if (!any(event == 1L)) {
    sb_abort(
      "A fixed ranger survival analysis fold must contain at least one event.",
      "stormbreaker_ranger_survival_no_event"
    )
  }
  if (!sb_ranger_available()) {
    sb_abort(
      "Package `ranger` is required for the fixed survival forest learner.",
      "stormbreaker_missing_ranger"
    )
  }

  effective_mtry <- if (is.null(params$mtry)) {
    max(1L, floor(sqrt(ncol(x))))
  } else {
    params$mtry
  }
  if (effective_mtry > ncol(x)) {
    sb_abort(
      "`ranger_mtry` cannot exceed the number of fitted features.",
      "stormbreaker_ranger_mtry_exceeds_features"
    )
  }

  if (identical(learner, "ranger_cindex_survival")) {
    threads <- getOption("stormbreaker.ranger_cindex_threads", params$num_threads)
    if (!is.numeric(threads) || length(threads) != 1L || is.na(threads) ||
      !threads %in% 1:4) {
      sb_abort("`stormbreaker.ranger_cindex_threads` must be an integer from 1 to 4.",
        "stormbreaker_bad_ranger_threads")
    }
    params$num_threads <- as.integer(threads)
  }
  fit_arguments <- list(
    x = as.matrix(x),
    y = survival::Surv(time, event),
    num.trees = params$num_trees,
    mtry = effective_mtry,
    importance = params$importance,
    write.forest = params$write_forest,
    min.node.size = params$min_node_size,
    max.depth = params$max_depth,
    replace = params$replace,
    sample.fraction = params$sample_fraction,
    splitrule = params$splitrule,
    oob.error = params$oob_error,
    num.threads = params$num_threads,
    seed = params$seed,
    na.action = params$na_action,
    verbose = FALSE
  )
  if ("num_random_splits" %in% names(params)) {
    fit_arguments$num.random.splits <- params$num_random_splits
  }
  if ("alpha" %in% names(params)) {
    fit_arguments$alpha <- params$alpha
  }
  if ("minprop" %in% names(params)) {
    fit_arguments$minprop <- params$minprop
  }
  fitted <- withCallingHandlers(
    do.call(ranger::ranger, fit_arguments),
    warning = function(warning) {
      message <- conditionMessage(warning)
      if (!startsWith(message, "Unused arguments:") && !identical(message,
        "Argument 'num.random.splits' ignored if splitrule is not 'extratrees'.")) {
        return(invisible(NULL))
      }
      sb_abort(
        paste0(
          "The fixed ranger survival learner raised a warning: ",
          message
        ),
        "stormbreaker_ranger_survival_warning",
        parent = warning
      )
    }
  )
  valid_forest <- inherits(fitted$forest, "ranger.forest") &&
    identical(as.integer(fitted$num.trees), params$num_trees) &&
    identical(fitted$treetype, "Survival")
  if (!valid_forest) {
    sb_abort(
      "Fixed ranger survival did not return the declared forest.",
      "stormbreaker_incomplete_ranger_fit"
    )
  }
  importance <- fitted$variable.importance
  valid_importance <- is.numeric(importance) &&
    length(importance) == ncol(x) && !anyNA(importance) &&
    all(is.finite(importance)) &&
    setequal(names(importance), names(x))
  if (!valid_importance) {
    sb_abort(
      "Fixed ranger survival did not return one finite importance per feature.",
      "stormbreaker_invalid_ranger_importance"
    )
  }
  importance <- unname(importance[names(x)])

  structure(
    list(
      learner = learner,
      features = names(x),
      params = params,
      num_trees = params$num_trees,
      effective_mtry = as.integer(effective_mtry),
      forest = fitted$forest,
      feature_importance = importance
    ),
    class = "stormbreaker_ranger_survival"
  )
}

sb_predict_ranger_survival_curves <- function(model, x) {
  if (!inherits(model, "stormbreaker_ranger_survival")) {
    sb_abort(
      "The supplied model is not a fitted fixed ranger survival learner.",
      "stormbreaker_bad_ranger_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed ranger survival schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (!sb_ranger_available()) {
    sb_abort(
      "Package `ranger` is required for fixed survival forest prediction.",
      "stormbreaker_missing_ranger"
    )
  }

  prediction <- stats::predict(
    model$forest,
    data = x,
    num.trees = model$num_trees,
    type = "response",
    seed = 0L,
    num.threads = 1L,
    verbose = FALSE
  )
  chf <- if (is.null(dim(prediction$chf))) {
    matrix(prediction$chf, nrow = nrow(x), byrow = TRUE)
  } else {
    as.matrix(prediction$chf)
  }
  survival_probability <- if (is.null(dim(prediction$survival))) {
    matrix(prediction$survival, nrow = nrow(x), byrow = TRUE)
  } else {
    as.matrix(prediction$survival)
  }
  times <- as.numeric(prediction$unique.death.times)
  valid <- nrow(chf) == nrow(x) && nrow(survival_probability) == nrow(x) &&
    ncol(chf) == length(times) &&
    identical(dim(chf), dim(survival_probability)) &&
    length(times) > 0L && !anyNA(times) && all(is.finite(times)) &&
    !anyNA(chf) && all(is.finite(chf)) &&
    !anyNA(survival_probability) && all(is.finite(survival_probability)) &&
    all(survival_probability >= 0 & survival_probability <= 1)
  if (!valid) {
    sb_abort(
      "Fixed ranger survival produced invalid prediction curves.",
      "stormbreaker_invalid_ranger_prediction"
    )
  }
  risk <- unname(rowSums(chf))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "Fixed ranger survival produced invalid mortality predictions.",
      "stormbreaker_invalid_ranger_prediction"
    )
  }
  list(
    risk = risk,
    time = times,
    survival = survival_probability
  )
}

sb_ranger_survival_probability_at_horizon <- function(curves,
                                                       prediction_horizon) {
  eligible <- which(curves$time <= prediction_horizon)
  probability <- if (length(eligible) == 0L) {
    rep(1, length(curves$risk))
  } else {
    unname(curves$survival[, max(eligible)])
  }
  if (length(probability) != length(curves$risk) || anyNA(probability) ||
      any(!is.finite(probability)) ||
      any(probability < 0 | probability > 1)) {
    sb_abort(
      "Fixed ranger survival produced invalid horizon probabilities.",
      "stormbreaker_invalid_ranger_prediction"
    )
  }
  probability
}

sb_predict_ranger_survival <- function(model, x) {
  sb_predict_ranger_survival_curves(model, x)$risk
}

sb_predict_ranger_survival_survival <- function(
    model, x_analysis, time, event, x_assessment, prediction_horizon) {
  x_analysis <- sb_validate_predictor_frame(x_analysis, "x_analysis")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x_analysis))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x_analysis))
  if (!identical(names(x_analysis), model$features)) {
    sb_abort(
      "Analysis columns do not match the fixed ranger survival schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
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
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    time
  )
  curves <- sb_predict_ranger_survival_curves(model, x_assessment)
  list(
    prediction_horizon = prediction_horizon,
    source = "ranger_survival_forest_direct",
    ties = "not_applicable",
    probability = sb_ranger_survival_probability_at_horizon(
      curves,
      prediction_horizon
    )
  )
}

sb_extract_ranger_survival_tree_model <- function(model) {
  if (!inherits(model, "stormbreaker_ranger_survival")) {
    sb_abort(
      "The supplied model is not a fitted fixed ranger survival learner.",
      "stormbreaker_bad_ranger_fit"
    )
  }
  list(
    num_trees = model$num_trees,
    effective_mtry = model$effective_mtry,
    features = model$features,
    params = model$params,
    importance_type = model$params$importance,
    model_serialization = "R_xdr_v3_ranger.forest",
    model_raw = serialize(
      model$forest,
      connection = NULL,
      ascii = FALSE,
      xdr = TRUE,
      version = 3L
    ),
    feature_importance = data.frame(
      feature = model$features,
      importance = model$feature_importance,
      stringsAsFactors = FALSE
    )
  )
}

sb_ranger_survival_learner_adapter <- function(
    num_trees = 500L,
    mtry = NULL,
    min_node_size = 3L,
    seed = 1L,
    learner = "ranger_survival") {
  params <- sb_ranger_survival_params(
    num_trees = num_trees,
    mtry = mtry,
    min_node_size = min_node_size,
    seed = seed,
    learner = learner
  )
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_ranger_survival(x, time, event, params, learner)
    },
    predict = sb_predict_ranger_survival,
    predict_survival = sb_predict_ranger_survival_survival,
    extract_tree_model = sb_extract_ranger_survival_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_ranger_survival_learner_state(model, context)
    },
    predict_state = sb_predict_ranger_survival_learner_state
  )
}

# ---- model-rlt ----

# ---- algorithm-rlt ----

#' @rdname sb_forest_diagnostics
#' @export
sb_rlt_diagnostics <- function(fit, type = "split_usage") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% sb_rlt_learner_ids()) sb_abort(
    "This diagnostic requires an RLT model.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("split_usage", "tree_size")) sb_abort(
    "Choose split_usage or tree_size.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  sb_validate_rlt_engine(e, s$learner, s$features, s$params)
  rows <- lapply(seq_along(e$FittedTrees), function(i) {
    z <- e$FittedTrees[[i]]
    if (!is.matrix(z) || !all(c("NodeType", "Node", "NextLeft", "NextRight", "NumOfComb", "SplitVar1") %in% rownames(z)) ||
        any(!z["NodeType", ] %in% c(1, 2))) sb_abort(
      "Saved native RLT tree structure is invalid.", "stormbreaker_invalid_diagnostic")
    internal <- z["NodeType", ] == 1
    v <- z["SplitVar1", internal]
    if (any(z["NumOfComb", internal] != 1) || any(!is.finite(v)) ||
        any(v < 1 | v > length(s$features) | v != floor(v))) sb_abort(
      "Expected original axis-aligned RLT split indices.", "stormbreaker_invalid_diagnostic")
    list(counts = tabulate(as.integer(v), nbins = length(s$features)),
      size = data.frame(tree = i, nodes = ncol(z), terminal_nodes = sum(!internal)))
  })
  tab <- if (type == "split_usage") data.frame(feature = s$features,
    splits = Reduce(`+`, lapply(rows, `[[`, "counts"))) else do.call(rbind, lapply(rows, `[[`, "size"))
  structure(list(type = paste0("rlt_", type), table = tab, learner = s$learner,
    pipeline_id = fit$pipeline$id, source = "original saved RLT node arrays",
    interpretation = paste("Every original tree and feature retained. Counts describe structural use, not variable importance or causal effects.",
      "Original importance computation was disabled; no OOB or reinforcement-training trajectory is reconstructed.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_rlt_diagnostic <- function(x) {
  feature <- splits <- terminal_nodes <- trees <- NULL
  d <- x$table
  if (x$type == "rlt_split_usage") {
    d$feature <- factor(d$feature, levels = d$feature[order(d$splits)])
    p <- ggplot2::ggplot(d, ggplot2::aes(splits, feature)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::labs(title = "RLT split-variable usage", x = "Number of original splits", y = NULL)
  } else {
    counts <- as.data.frame(table(d$terminal_nodes))
    names(counts) <- c("terminal_nodes", "trees")
    counts$terminal_nodes <- as.numeric(as.character(counts$terminal_nodes))
    p <- ggplot2::ggplot(counts, ggplot2::aes(terminal_nodes, trees)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::scale_y_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::labs(title = "RLT tree-size distribution", x = "Terminal nodes per tree", y = "Number of trees")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}

# ---- rlt-adapters ----

sb_rlt_binary_learner_ids <- function() {
  "rlt_reinforcement_binary"
}

sb_rlt_survival_learner_ids <- function() {
  "rlt_reinforcement_survival"
}

sb_rlt_learner_ids <- function() {
  c(sb_rlt_binary_learner_ids(), sb_rlt_survival_learner_ids())
}

sb_rlt_available <- function() {
  requireNamespace("RLT", quietly = TRUE)
}

sb_rlt_dependency <- function() {
  list(package = "RLT", available = sb_rlt_available())
}

sb_require_rlt_dependency <- function() {
  if (!sb_rlt_available()) {
    sb_abort(
      "Package `RLT` is required for reinforcement-learning forests.",
      "stormbreaker_missing_rlt"
    )
  }
  invisible(TRUE)
}

sb_rlt_learner_params <- function(learner) {
  if (!learner %in% sb_rlt_learner_ids()) {
    sb_abort(
      paste0("Unknown RLT learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    model = if (learner %in% sb_rlt_binary_learner_ids()) {
      "classification"
    } else {
      "survival"
    },
    num_trees = 100L,
    reinforcement = TRUE,
    combination_split_size = 1L,
    native_nmin = 5L,
    embedded_stop_size = 1L,
    cores = 1L,
    importance = FALSE,
    track_observations = FALSE,
    seed = 1L,
    split_rule = "axis_aligned",
    remaining_controls = "RLT_native_defaults"
  )
}

sb_rlt_parameter_semantics <- function(learner) {
  params <- sb_rlt_learner_params(learner)
  paste0(
    "reinforcement-learning forest with ", params$num_trees,
    " trees, native internal-node split threshold nmin ", params$native_nmin,
    ", axis-aligned splits, embedded-model stop size ",
    params$embedded_stop_size,
    ", one-core execution, and seed ", params$seed,
    " fixed; remaining split, resampling, and muting controls use native ",
    "RLT defaults"
  )
}

sb_rlt_stop_on_warning <- function(code, learner, stage) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      message <- conditionMessage(warning)
      resets_weights <- grepl("^(Subject|Variable) weights .*reset to equal weights$", message)
      mismatches_features <- identical(message,
        "test data variables names does not match training data...")
      if (!resets_weights && !mismatches_features) return(invisible(NULL))
      sb_abort(
        paste0(
          "Learner `", learner, "` raised a ", stage, " warning: ",
          message
        ),
        "stormbreaker_rlt_warning",
        parent = warning
      )
    }
  )
}

sb_validate_rlt_survival_native_domain <- function(inputs, learner) {
  if (any(inputs$time <= 0)) {
    sb_abort(
      "RLT survival fitting requires strictly positive follow-up times.",
      "stormbreaker_bad_rlt_time"
    )
  }
  distinct_event_times <- unique(inputs$time[inputs$event == 1L])
  if (length(distinct_event_times) < 2L) {
    # RLT 3.2.6 terminates the Windows R process in compiled code when its
    # survival time grid has one point. An R condition cannot catch that native
    # process exit, so this check belongs immediately before the unsafe call.
    sb_abort(
      paste0(
        "Learner `", learner, "` requires at least two distinct event times ",
        "because the native RLT survival engine cannot safely fit a one-point ",
        "event-time grid."
      ),
      "stormbreaker_rlt_insufficient_event_times"
    )
  }
  inputs
}

sb_validate_rlt_engine <- function(engine, learner, features, params) {
  expected_model <- params$model
  valid <- inherits(engine, "RLT") && is.list(engine) &&
    identical(engine$model, expected_model) &&
    is.numeric(engine$p) && length(engine$p) == 1L &&
    identical(as.integer(engine$p), as.integer(length(features))) &&
    identical(engine$variablenames, features) &&
    is.numeric(engine$ntrees) && length(engine$ntrees) == 1L &&
    identical(as.integer(engine$ntrees), params$num_trees) &&
    is.list(engine$FittedTrees) &&
    length(engine$FittedTrees) == params$num_trees &&
    is.numeric(engine$ncat) && length(engine$ncat) == length(features) &&
    is.numeric(engine$parameters.int) &&
    is.numeric(engine$parameters.double)
  if (valid && learner %in% sb_rlt_binary_learner_ids()) {
    valid <- identical(engine$ylevels, c("control", "case")) &&
      identical(as.integer(engine$nclass), 2L)
  }
  if (valid && learner %in% sb_rlt_survival_learner_ids()) {
    valid <- is.list(engine$FittedSurv) &&
      length(engine$FittedSurv) == params$num_trees &&
      is.numeric(engine$timepoints) && length(engine$timepoints) >= 2L &&
      !anyNA(engine$timepoints) && all(is.finite(engine$timepoints)) &&
      all(engine$timepoints > 0) && all(diff(engine$timepoints) > 0)
  }
  if (!isTRUE(valid)) {
    sb_abort(
      paste0("Learner `", learner, "` returned an invalid native RLT fit."),
      "stormbreaker_invalid_rlt_fit"
    )
  }
  engine
}

sb_rlt_prediction_engine <- function(engine, learner, features, params) {
  engine <- sb_validate_rlt_engine(engine, learner, features, params)
  keep <- c(
    "FittedTrees", "FittedSurv", "timepoints", "ylevels", "nclass",
    "model", "p", "variablenames", "ncat", "xlevels", "ntrees",
    "parameters.int", "parameters.double"
  )
  compact <- engine[intersect(keep, names(engine))]
  class(compact) <- "RLT"
  sb_validate_rlt_engine(compact, learner, features, params)
}

sb_fit_rlt_learner <- function(x, endpoint, learner, params) {
  sb_require_rlt_dependency()
  inputs <- if (learner %in% sb_rlt_binary_learner_ids()) {
    sb_validate_classical_binary_inputs(x, endpoint$outcome, learner)
  } else {
    sb_validate_classical_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      learner
    ) |>
      sb_validate_rlt_survival_native_domain(learner)
  }
  matrix_x <- as.matrix(inputs$x)
  storage.mode(matrix_x) <- "double"
  arguments <- list(
    x = matrix_x,
    model = params$model,
    print.summary = 0L,
    use.cores = params$cores,
    ntrees = params$num_trees,
    nmin = params$native_nmin,
    reinforcement = params$reinforcement,
    combsplit = params$combination_split_size,
    embed.n.th = params$embedded_stop_size,
    importance = params$importance,
    track.obs = params$track_observations
  )
  if (learner %in% sb_rlt_binary_learner_ids()) {
    arguments$y <- factor(
      inputs$outcome,
      levels = c(0L, 1L),
      labels = c("control", "case")
    )
  } else {
    arguments$y <- inputs$time
    arguments$censor <- inputs$event
  }
  engine <- sb_rlt_stop_on_warning(
    withr::with_seed(params$seed, do.call(RLT::RLT, arguments)),
    learner,
    "fitting"
  )
  engine <- sb_rlt_prediction_engine(
    engine,
    learner,
    names(inputs$x),
    params
  )
  structure(
    list(
      learner = learner,
      features = names(inputs$x),
      params = params,
      engine = engine
    ),
    class = "stormbreaker_rlt_model"
  )
}

sb_validate_rlt_model_x <- function(model, x) {
  valid <- inherits(model, "stormbreaker_rlt_model") && is.list(model) &&
    model$learner %in% sb_rlt_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.list(model$params)
  if (!isTRUE(valid)) {
    sb_abort(
      "The supplied model is not a valid fitted RLT learner.",
      "stormbreaker_bad_rlt_fit"
    )
  }
  sb_validate_rlt_engine(
    model$engine,
    model$learner,
    model$features,
    model$params
  )
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted RLT learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  x
}

sb_predict_rlt_native <- function(model, x) {
  x <- sb_validate_rlt_model_x(model, x)
  matrix_x <- as.matrix(x)
  storage.mode(matrix_x) <- "double"
  predict_method <- get(
    "predict.RLT",
    envir = asNamespace("RLT"),
    inherits = FALSE
  )
  prediction <- sb_rlt_stop_on_warning(
    predict_method(model$engine, testx = matrix_x),
    model$learner,
    "prediction"
  )
  list(x = x, prediction = prediction)
}

sb_predict_rlt_binary <- function(model, x) {
  native <- sb_predict_rlt_native(model, x)
  if (!model$learner %in% sb_rlt_binary_learner_ids()) {
    sb_abort(
      "A survival RLT model cannot produce binary probabilities.",
      "stormbreaker_bad_rlt_fit"
    )
  }
  probabilities <- as.matrix(native$prediction$ProbPrediction)
  valid <- nrow(probabilities) == nrow(native$x) &&
    ncol(probabilities) == 2L && !anyNA(probabilities) &&
    all(is.finite(probabilities)) &&
    all(probabilities >= 0 & probabilities <= 1) &&
    all(abs(rowSums(probabilities) - 1) <= 1e-10)
  if (!isTRUE(valid)) {
    sb_abort(
      "The RLT classification forest produced invalid class probabilities.",
      "stormbreaker_invalid_rlt_prediction"
    )
  }
  unname(probabilities[, 2L])
}

sb_predict_rlt_survival_curves <- function(model, x) {
  native <- sb_predict_rlt_native(model, x)
  if (!model$learner %in% sb_rlt_survival_learner_ids()) {
    sb_abort(
      "A binary RLT model cannot produce survival curves.",
      "stormbreaker_bad_rlt_fit"
    )
  }
  hazard <- as.matrix(native$prediction$SurvPred)
  timepoints <- as.numeric(native$prediction$timepoints)
  valid <- length(timepoints) >= 2L && !anyNA(timepoints) &&
    all(is.finite(timepoints)) && all(timepoints > 0) &&
    all(diff(timepoints) > 0) && nrow(hazard) == nrow(native$x) &&
    ncol(hazard) == length(timepoints) && !anyNA(hazard) &&
    all(is.finite(hazard)) && all(hazard >= 0 & hazard <= 1)
  if (!isTRUE(valid)) {
    sb_abort(
      "The RLT survival forest produced invalid discrete hazard increments.",
      "stormbreaker_invalid_rlt_prediction"
    )
  }
  survival <- matrix(
    NA_real_,
    nrow = nrow(hazard),
    ncol = ncol(hazard)
  )
  for (row in seq_len(nrow(hazard))) {
    survival[row, ] <- cumprod(1 - hazard[row, ])
  }
  valid_survival <- !anyNA(survival) && all(is.finite(survival)) &&
    all(survival >= 0 & survival <= 1) &&
    all(apply(survival, 1L, function(values) all(diff(values) <= 0)))
  if (!isTRUE(valid_survival)) {
    sb_abort(
      "The reconstructed RLT survival curves are invalid.",
      "stormbreaker_invalid_rlt_prediction"
    )
  }
  interval <- diff(c(0, timepoints))
  previous <- cbind(
    1,
    survival[, -ncol(survival), drop = FALSE]
  )
  risk <- -as.vector(previous %*% interval)
  if (length(risk) != nrow(native$x) || anyNA(risk) ||
    any(!is.finite(risk))) {
    sb_abort(
      "The RLT survival forest produced invalid restricted-mean risk scores.",
      "stormbreaker_invalid_rlt_prediction"
    )
  }
  list(
    risk = unname(risk),
    time = timepoints,
    hazard = hazard,
    survival = survival
  )
}

sb_predict_rlt_learner <- function(model, x) {
  if (model$learner %in% sb_rlt_binary_learner_ids()) {
    return(sb_predict_rlt_binary(model, x))
  }
  sb_predict_rlt_survival_curves(model, x)$risk
}

sb_rlt_direct_survival_source <- function() {
  "rlt_reconstructed_survival_from_native_discrete_hazard"
}

sb_predict_rlt_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  if (!model$learner %in% sb_rlt_survival_learner_ids()) {
    sb_abort(
      "A binary RLT learner cannot predict survival probability.",
      "stormbreaker_bad_rlt_fit"
    )
  }
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted RLT learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  curves <- sb_predict_rlt_survival_curves(model, x_assessment)
  eligible <- which(curves$time <= prediction_horizon)
  probability <- if (length(eligible) == 0L) {
    rep(1, length(curves$risk))
  } else {
    unname(curves$survival[, max(eligible)])
  }
  if (length(probability) != nrow(x_assessment) || anyNA(probability) ||
    any(!is.finite(probability)) || any(probability < 0 | probability > 1)) {
    sb_abort(
      "The RLT survival forest produced invalid horizon probabilities.",
      "stormbreaker_invalid_rlt_prediction"
    )
  }
  list(
    prediction_horizon = prediction_horizon,
    source = sb_rlt_direct_survival_source(),
    ties = "not_applicable",
    probability = probability
  )
}

sb_extract_rlt_state <- function(model, context = list()) {
  sb_validate_rlt_model_x(
    model,
    as.data.frame(
      matrix(
        0,
        nrow = 1L,
        ncol = length(model$features),
        dimnames = list(NULL, model$features)
      ),
      check.names = FALSE
    )
  )
  survival_probability <- NULL
  if (model$learner %in% sb_rlt_survival_learner_ids() &&
    "prediction_horizon" %in% names(context)) {
    probability <- sb_predict_rlt_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    survival_probability <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      source = probability$source,
      ties = probability$ties
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      engine_serialization = "R_xdr_v3_RLT_compact_prediction_engine",
      engine_raw = serialize(
        model$engine,
        connection = NULL,
        ascii = FALSE,
        xdr = TRUE,
        version = 3L
      ),
      risk_reference = if (
        model$learner %in% sb_rlt_survival_learner_ids()
      ) {
        "RLT_native_discrete_hazard"
      } else {
        NULL
      },
      risk_formula = if (
        model$learner %in% sb_rlt_survival_learner_ids()
      ) {
        "negative_step_RMST_over_native_event_grid"
      } else {
        NULL
      },
      baseline_hazard = NULL,
      survival_probability = survival_probability
    ),
    class = c("stormbreaker_rlt_learner_state", "stormbreaker_learner_state")
  )
}

sb_restore_rlt_model <- function(state, learner) {
  if (!is.raw(state$engine_raw) || length(state$engine_raw) == 0L) {
    sb_bad_learner_state("The saved RLT prediction engine is invalid.")
  }
  engine <- tryCatch(
    unserialize(state$engine_raw),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved RLT prediction engine could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_rlt_state",
        parent = error
      )
    }
  )
  engine <- sb_validate_rlt_engine(
    engine,
    learner,
    state$features,
    state$params
  )
  structure(
    list(
      learner = learner,
      features = state$features,
      params = state$params,
      engine = engine
    ),
    class = "stormbreaker_rlt_model"
  )
}

sb_validate_rlt_probability_state <- function(state) {
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
    identical(probability$source, sb_rlt_direct_survival_source()) &&
    identical(probability$ties, "not_applicable")
  if (!isTRUE(valid)) {
    sb_bad_learner_state(
      "The saved RLT direct survival-probability state is invalid."
    )
  }
  probability
}

sb_predict_rlt_state <- function(state, x) {
  learner <- if (is.list(state)) state$learner else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_rlt_learner_ids()) {
    sb_bad_learner_state("The saved RLT learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  sb_require_rlt_dependency()
  model <- sb_restore_rlt_model(state, learner)
  if (learner %in% sb_rlt_binary_learner_ids()) {
    return(list(probability = sb_predict_rlt_binary(model, x)))
  }
  curves <- sb_predict_rlt_survival_curves(model, x)
  result <- list(risk = curves$risk)
  if (!is.null(state$survival_probability)) {
    probability <- sb_validate_rlt_probability_state(state)
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

sb_rlt_learner_adapter <- function(
  learner,
  params = sb_rlt_learner_params(learner)
) {
  if (!learner %in% sb_rlt_learner_ids()) {
    sb_abort(
      paste0("Unknown RLT learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  adapter <- list(
    id = learner,
    params = params,
    fit = if (learner %in% sb_rlt_binary_learner_ids()) {
      function(x, outcome) {
        sb_fit_rlt_learner(x, list(outcome = outcome), learner, params)
      }
    } else {
      function(x, time, event) {
        sb_fit_rlt_learner(
          x,
          list(time = time, event = event),
          learner,
          params
        )
      }
    },
    predict = sb_predict_rlt_learner,
    extract_state = sb_extract_rlt_state,
    predict_state = sb_predict_rlt_state
  )
  if (learner %in% sb_rlt_survival_learner_ids()) {
    adapter$predict_survival <- sb_predict_rlt_survival
  }
  adapter
}

# ---- model-wsrf ----

# ---- algorithm-wsrf ----

#' @rdname sb_forest_diagnostics
#' @export
sb_wsrf_diagnostics <- function(fit, type = "importance") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!identical(s$learner, "wsrf_binary")) sb_abort(
    "This diagnostic requires WSRF.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("importance", "tree_nodes", "tree_oob")) sb_abort(
    "Choose importance, tree_nodes or tree_oob.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  n <- length(e$trees)
  m <- NULL
  if (type == "tree_nodes") {
    tab <- data.frame(tree = seq_len(n), nodes = lengths(e$trees))
  } else {
    m <- fit$diagnostic_artifacts$learner$diagnostic_metrics
    if (is.null(m)) sb_abort("Original WSRF metrics were not captured; compact placeholders are not diagnostic values.",
      "stormbreaker_missing_diagnostic_artifacts")
    if (type == "tree_oob") {
      if (!is.numeric(m$tree_oob) || length(m$tree_oob) != n || any(is.infinite(m$tree_oob)) ||
          any(m$tree_oob < 0 | m$tree_oob > 1, na.rm = TRUE)) sb_abort(
        "Original WSRF tree OOB values are invalid.", "stormbreaker_invalid_diagnostic")
      tab <- data.frame(tree = seq_len(n), error = m$tree_oob, available = is.finite(m$tree_oob))
    } else {
      z <- m$importance
      if (!is.matrix(z) || !"MeanDecreaseIGR" %in% colnames(z) ||
          !identical(rownames(z), s$internal_features) || any(!is.finite(z))) sb_abort(
        "Original WSRF information-gain-ratio importance is invalid.", "stormbreaker_invalid_diagnostic")
      bridge <- seq_len(nrow(z)) > length(s$features)
      tab <- data.frame(feature = c(s$features, if (any(bridge)) "Constant bridge"),
        internal_feature = s$internal_features, artificial_bridge = bridge,
        importance = unname(z[, "MeanDecreaseIGR"]))
    }
  }
  structure(list(type = paste0("wsrf_", type), table = tab, learner = s$learner,
    pipeline_id = fit$pipeline$id, forest_oob = m$forest_oob,
    source = "original WSRF node lists and separately captured pre-compaction metrics",
    interpretation = paste("Importance is native mean decrease in information gain ratio, not Gini or permutation importance.",
      "Node counts include internal and terminal nodes. OOB points are individual-tree errors, not a cumulative learning curve.",
      "Compact prediction placeholders are never reported as OOB results. No refitting or uncertainty reconstruction.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_wsrf_diagnostic <- function(x) {
  importance <- feature <- nodes <- trees <- tree <- error <- NULL
  d <- x$table
  if (x$type == "wsrf_importance") {
    d$feature <- factor(d$feature, levels = d$feature[order(d$importance)])
    p <- ggplot2::ggplot(d, ggplot2::aes(importance, feature)) +
      ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::labs(title = "Weighted-subspace forest importance", x = "Mean decrease in information gain ratio", y = NULL)
  } else if (x$type == "wsrf_tree_nodes") {
    counts <- as.data.frame(table(d$nodes))
    names(counts) <- c("nodes", "trees")
    counts$nodes <- as.numeric(as.character(counts$nodes))
    p <- ggplot2::ggplot(counts, ggplot2::aes(nodes, trees)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::labs(title = "Weighted-subspace forest sizes", x = "Total nodes per tree", y = "Number of trees")
  } else p <- ggplot2::ggplot(d[d$available, , drop = FALSE], ggplot2::aes(tree, error)) +
    ggplot2::geom_point(colour = "#2878A0", size = 1) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = "Individual-tree OOB errors", x = "Tree index", y = "OOB misclassification rate")
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}

# ---- wsrf-binary-adapter ----

sb_wsrf_binary_learner_ids <- function() {
  "wsrf_binary"
}

sb_wsrf_available <- function() {
  requireNamespace("wsrf", quietly = TRUE)
}

sb_require_wsrf_dependency <- function() {
  if (!sb_wsrf_available()) {
    sb_abort(
      paste0(
        "Package `wsrf` is required for the weighted-subspace binary forest."
      ),
      "stormbreaker_missing_wsrf"
    )
  }
  invisible(TRUE)
}

sb_wsrf_binary_params <- function(
  learner = sb_wsrf_binary_learner_ids()
) {
  if (!learner %in% sb_wsrf_binary_learner_ids()) {
    sb_abort(
      paste0("Unknown weighted-subspace forest learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    num_trees = 500L,
    mtry_rule = "floor(log2(p)+1)",
    weighted_subspace = TRUE,
    parallel = FALSE,
    importance = FALSE,
    leaf_node_size = 2L,
    seed = 1L,
    probability = "native_vote_proportion"
  )
}

sb_wsrf_binary_parameter_semantics <- function(learner) {
  params <- sb_wsrf_binary_params(learner)
  paste0(
    "weighted-subspace C4.5 classification forest with ",
    params$num_trees, " trees, mtry=floor(log2(p)+1), leaf node size ",
    params$leaf_node_size, ", native vote probabilities, serial execution, ",
    "and seed ", params$seed, " fixed before resampling"
  )
}



sb_wsrf_internal_names <- function(p) {
  paste0(".stormbreaker_wsrf_x", seq_len(p))
}

sb_wsrf_internal_data <- function(x) {
  internal <- x
  names(internal) <- sb_wsrf_internal_names(ncol(x))
  uses_dummy <- ncol(x) == 1L
  if (uses_dummy) {
    internal$.stormbreaker_wsrf_zero <- 0
  }
  list(
    data = internal,
    internal_features = names(internal),
    uses_zero_variance_bridge = uses_dummy
  )
}

sb_validate_wsrf_engine <- function(engine, internal_features, num_trees) {
  required <- c(
    "meta", "targetData", "trees", "treeOOBErrorRates", "OOBSets",
    "OOBPredictLabels", "treeIgrImportance"
  )
  valid <- inherits(engine, "wsrf") && is.list(engine) &&
    all(required %in% names(engine)) && is.list(engine$trees) &&
    length(engine$trees) == num_trees &&
    is.list(engine$meta) &&
    is.character(internal_features) && length(internal_features) > 0L
  if (!isTRUE(valid)) {
    sb_abort(
      "The weighted-subspace forest returned an invalid native fit.",
      "stormbreaker_invalid_wsrf_fit"
    )
  }
  engine
}

sb_wsrf_compact_prediction_engine <- function(engine, internal_features,
                                               num_trees) {
  engine <- sb_validate_wsrf_engine(engine, internal_features, num_trees)
  engine["targetData"] <- list(NULL)
  engine$treeOOBErrorRates <- rep(0, num_trees)
  engine$OOBSets <- replicate(num_trees, integer(), simplify = FALSE)
  engine$OOBPredictLabels <- replicate(
    num_trees,
    integer(),
    simplify = FALSE
  )
  engine$treeIgrImportance <- replicate(
    num_trees,
    numeric(),
    simplify = FALSE
  )
  removable <- c(
    "predicted", "oob.times", "confusion", "importance", "importanceSD",
    "RFOOBErrorRate", "strength", "correlation", "c_s2", "useweights",
    "mtry", "nodesize"
  )
  for (name in intersect(removable, names(engine))) {
    engine[name] <- list(NULL)
  }
  class(engine) <- "wsrf"
  sb_validate_wsrf_engine(engine, internal_features, num_trees)
}

sb_validate_wsrf_model <- function(model) {
  valid <- inherits(model, "stormbreaker_wsrf_binary_model") &&
    is.list(model) && identical(model$learner, sb_wsrf_binary_learner_ids()) &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.character(model$internal_features) &&
    length(model$internal_features) %in% c(
      length(model$features),
      length(model$features) + 1L
    ) &&
    is.logical(model$uses_zero_variance_bridge) &&
    length(model$uses_zero_variance_bridge) == 1L &&
    is.numeric(model$effective_mtry) && length(model$effective_mtry) == 1L &&
    !is.na(model$effective_mtry) && model$effective_mtry >= 1L &&
    is.list(model$params)
  if (!isTRUE(valid)) {
    sb_abort(
      "The supplied model is not a valid fitted weighted-subspace forest.",
      "stormbreaker_bad_wsrf_fit"
    )
  }
  sb_validate_wsrf_engine(
    model$engine,
    model$internal_features,
    model$params$num_trees
  )
  model
}

sb_fit_wsrf_binary <- function(x, outcome, learner, params) {
  sb_require_wsrf_dependency()
  inputs <- sb_validate_classical_binary_inputs(x, outcome, learner)
  prepared <- sb_wsrf_internal_data(inputs$x)
  effective_mtry <- as.integer(floor(log2(ncol(inputs$x)) + 1))
  grouping <- factor(
    inputs$outcome,
    levels = c(0L, 1L),
    labels = c("control", "case")
  )
  engine <- withr::with_seed(
    params$seed,
    wsrf::wsrf(
      x = prepared$data,
      y = grouping,
      mtry = effective_mtry,
      ntree = params$num_trees,
      weights = params$weighted_subspace,
      parallel = params$parallel,
      importance = params$importance,
      nodesize = params$leaf_node_size
    )
  )
  diagnostic_metrics <- list(importance = engine$importance,
    tree_oob = engine$treeOOBErrorRates, forest_oob = engine$RFOOBErrorRate)
  engine <- sb_wsrf_compact_prediction_engine(
    engine,
    prepared$internal_features,
    params$num_trees
  )
  model <- structure(
    list(
      learner = learner,
      features = names(inputs$x),
      internal_features = prepared$internal_features,
      uses_zero_variance_bridge = prepared$uses_zero_variance_bridge,
      effective_mtry = effective_mtry,
      diagnostic_metrics = diagnostic_metrics,
      params = params,
      engine = engine
    ),
    class = "stormbreaker_wsrf_binary_model"
  )
  sb_validate_wsrf_model(model)
}

sb_wsrf_prediction_data <- function(model, x) {
  model <- sb_validate_wsrf_model(model)
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      paste0(
        "Prediction columns do not match the fitted weighted-subspace forest ",
        "schema."
      ),
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prepared <- sb_wsrf_internal_data(x)
  if (!identical(prepared$internal_features, model$internal_features) ||
    !identical(
      prepared$uses_zero_variance_bridge,
      model$uses_zero_variance_bridge
    )) {
    sb_abort(
      "The internal weighted-subspace forest predictor schema is invalid.",
      "stormbreaker_bad_wsrf_fit"
    )
  }
  list(model = model, x = x, internal = prepared$data)
}

sb_predict_wsrf_binary <- function(model, x) {
  prepared <- sb_wsrf_prediction_data(model, x)
  prediction <- stats::predict(
    prepared$model$engine,
    newdata = prepared$internal,
    type = "prob"
  )
  probability_matrix <- prediction$prob
  valid <- is.matrix(probability_matrix) &&
    nrow(probability_matrix) == nrow(prepared$x) &&
    identical(colnames(probability_matrix), c("control", "case")) &&
    !anyNA(probability_matrix) && all(is.finite(probability_matrix)) &&
    all(probability_matrix >= 0 & probability_matrix <= 1) &&
    all(abs(rowSums(probability_matrix) - 1) <= 1e-10)
  if (!isTRUE(valid)) {
    sb_abort(
      "The weighted-subspace forest produced invalid class probabilities.",
      "stormbreaker_invalid_wsrf_prediction"
    )
  }
  unname(probability_matrix[, "case"])
}

sb_extract_wsrf_binary_state <- function(model, context = list()) {
  model <- sb_validate_wsrf_model(model)
  structure(
    list(
      learner = model$learner,
      features = model$features,
      internal_features = model$internal_features,
      uses_zero_variance_bridge = model$uses_zero_variance_bridge,
      effective_mtry = model$effective_mtry,
      params = model$params,
      positive_class = 1L,
      native_positive_label = "case",
      probability = "native_vote_proportion",
      engine_serialization = "R_xdr_v3_wsrf_compact_prediction_engine",
      engine_raw = serialize(
        model$engine,
        connection = NULL,
        ascii = FALSE,
        xdr = TRUE,
        version = 3L
      ),
      training_outcome_retained = FALSE,
      oob_patient_indices_retained = FALSE
    ),
    class = c(
      "stormbreaker_wsrf_binary_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_wsrf_binary_model <- function(state) {
  if (!is.raw(state$engine_raw) || length(state$engine_raw) == 0L) {
    sb_bad_learner_state(
      "The saved weighted-subspace forest prediction engine is invalid."
    )
  }
  engine <- tryCatch(
    unserialize(state$engine_raw),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved weighted-subspace forest engine could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_wsrf_state",
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
      effective_mtry = state$effective_mtry,
      params = state$params,
      engine = engine
    ),
    class = "stormbreaker_wsrf_binary_model"
  )
  sb_validate_wsrf_model(model)
}

sb_predict_wsrf_binary_state <- function(state, x) {
  valid <- inherits(state, "stormbreaker_wsrf_binary_state") &&
    is.list(state) && identical(state$learner, sb_wsrf_binary_learner_ids()) &&
    identical(state$positive_class, 1L) &&
    identical(state$native_positive_label, "case") &&
    identical(state$probability, "native_vote_proportion") &&
    identical(state$training_outcome_retained, FALSE) &&
    identical(state$oob_patient_indices_retained, FALSE)
  if (!isTRUE(valid)) {
    sb_bad_learner_state(
      "The saved weighted-subspace forest state is invalid."
    )
  }
  x <- sb_validate_learner_state_x(state, x, state$learner)
  sb_require_wsrf_dependency()
  model <- sb_restore_wsrf_binary_model(state)
  list(probability = sb_predict_wsrf_binary(model, x))
}

sb_wsrf_binary_learner_adapter <- function(
  learner,
  params = sb_wsrf_binary_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_wsrf_binary(x, outcome, learner, params)
    },
    predict = sb_predict_wsrf_binary,
    extract_state = sb_extract_wsrf_binary_state,
    predict_state = sb_predict_wsrf_binary_state
  )
}

# ---- diagnostics-trees ----

# ---- algorithm-c50 ----

#' @rdname sb_forest_diagnostics
#' @export
sb_c50_diagnostics <- function(fit,type="usage") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% c("c50_binary","c50_rules_binary") ||
      !is.character(type) || length(type)!=1L || is.na(type) || !type %in% c("usage","splits","size")) sb_abort(
    "Choose usage, splits or size for a C5.0 ensemble.","stormbreaker_unsupported_diagnostic")
  e <- fit$diagnostic_artifacts$learner$engine
  if(is.null(e)) sb_abort("Original C5.0 summaries were not captured.","stormbreaker_missing_diagnostic_artifacts")
  if(!inherits(e,"C5.0")) sb_abort("Captured C5.0 engine is invalid.","stormbreaker_invalid_diagnostic")
  if(type=="size") {
    if(length(e$size)!=e$trials[["Actual"]] || any(!is.finite(e$size)) || any(e$size<0)) sb_abort(
      "Native C5.0 sizes are invalid.","stormbreaker_invalid_diagnostic")
    tab <- data.frame(trial=seq_along(e$size),size=e$size)
  } else {
    native <- tryCatch(C50::C5imp(e,metric=type),error=identity)
    if(inherits(native,"error")) sb_abort(paste("Native C5.0 importance unavailable:",conditionMessage(native)),
      "stormbreaker_unavailable_diagnostic")
    if(is.data.frame(native) && "Overall" %in% names(native) && all(is.nan(native$Overall))) sb_abort(
      "Native C5.0 split percentages are undefined because no split importance was reported.","stormbreaker_unavailable_diagnostic")
    if(!is.data.frame(native) || !"Overall" %in% names(native) ||
        !setequal(rownames(native),s$features) || any(!is.finite(native$Overall))) sb_abort(
      "Native C5.0 feature summary is invalid.","stormbreaker_invalid_diagnostic")
    tab <- data.frame(feature=s$features,value=native[s$features,"Overall"])
  }
  structure(list(type=paste0("c50_",type),table=tab,learner=s$learner,pipeline_id=fit$pipeline$id,
    rules=isTRUE(e$rbm),source="Original captured C5.0 summaries and native C5imp",
    interpretation="Native usage and split percentages are different summaries, not permutation importance or causal effects. Trial sizes are not performance. No refitting or missing-value reconstruction."),
    class="stormbreaker_algorithm_diagnostic")
}

sb_plot_c50_diagnostic <- function(x) {
  feature <- value <- trial <- size <- NULL
  if(x$type=="c50_size") {
    p <- ggplot2::ggplot(x$table,ggplot2::aes(trial,size)) + ggplot2::geom_point(colour="#2878A0") +
      ggplot2::scale_x_continuous(breaks=x$table$trial) +
      ggplot2::scale_y_continuous(breaks=function(z) unique(round(pretty(z)))) +
      ggplot2::labs(title="C5.0 trial sizes",x="Trial",y=if(x$rules) "Number of rules" else "Number of leaves")
  } else {
    d <- x$table; d$feature <- factor(d$feature,levels=d$feature[order(d$value)])
    p <- ggplot2::ggplot(d,ggplot2::aes(value,feature)) + ggplot2::geom_col(fill="#2878A0") +
      ggplot2::labs(title="C5.0 native feature summary",x=if(x$type=="c50_usage") "Native usage (%)" else "Native splits (%)",y=NULL)
  }
  p + sb_viz_theme()
}

# ---- algorithm-classification-tree ----

#' @rdname sb_forest_diagnostics
#' @export
sb_classification_tree_diagnostics <- function(fit, type = "structure") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (s$learner != "tree_binary" || !is.character(type) || length(type) != 1L ||
      is.na(type) || !type %in% c("structure", "leaf_probabilities")) sb_abort(
    "Choose structure or leaf_probabilities for the tree-package classifier.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  f <- e$frame
  ids <- as.integer(rownames(f))
  if (!inherits(e, "tree") || !identical(attr(e,"ylevels"), c("0","1")) ||
      !is.matrix(f$yprob) || !identical(colnames(f$yprob), c("0","1"))) sb_abort(
    "Native classification tree or class probabilities are invalid.", "stormbreaker_invalid_diagnostic")
  internal <- sb_tree_binary_internal_names(length(s$features))
  vars <- as.character(f$var)
  terminal <- vars == "<leaf>"
  if (any(!vars[!terminal] %in% internal)) sb_abort("A native split feature is unknown.", "stormbreaker_invalid_diagnostic")
  feature <- s$features[match(vars,internal)]
  parent <- ifelse(ids == 1L, NA_integer_, ids %/% 2L)
  tab <- data.frame(node=ids,parent=parent,depth=floor(log2(ids)),terminal=terminal,
    incoming_rule="Root",split_feature=feature,n=f$n,probability=f$yprob[,"1"],
    x=NA_real_,y=-floor(log2(ids)))
  for (i in which(ids != 1L)) {
    j <- match(parent[i],ids)
    tab$incoming_rule[i] <- paste(feature[j], f$splits[j,if(ids[i] %% 2L == 0L) 1L else 2L])
  }
  position <- 0L
  layout <- function(id) {
    i <- match(id,ids)
    if (terminal[i]) { position <<- position+1L; tab$x[i] <<- position } else {
      layout(id*2L); layout(id*2L+1L)
      tab$x[i] <<- mean(range(tab$x[match(c(id*2L,id*2L+1L),ids)]))
    }
  }
  layout(1L)
  if (type == "leaf_probabilities") tab <- tab[tab$terminal, , drop=FALSE]
  structure(list(type=if(type=="structure") "ctree_structure" else "classification_tree_leaf_probabilities",
    table=tab,native_frame=f,learner=s$learner,pipeline_id=fit$pipeline$id,
    source="Original tree frame, native split strings and class probabilities",
    interpretation="Original training node counts and leaf probabilities; not held-out performance. No pruning or refitting."),
    class="stormbreaker_algorithm_diagnostic")
}

sb_plot_classification_tree_probabilities <- function(x) {
  node <- probability <- n <- NULL
  ggplot2::ggplot(x$table,ggplot2::aes(factor(node),probability,size=n)) +
    ggplot2::geom_point(colour="#2878A0") + ggplot2::scale_y_continuous(limits=c(0,1)) +
    ggplot2::labs(title="Terminal-node probabilities",x="Node",y="Class-1 probability",size="Training rows") + sb_viz_theme()
}

# ---- algorithm-ctree ----

#' @rdname sb_forest_diagnostics
#' @export
sb_ctree_structure <- function(fit) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!fit$pipeline$learner %in% c("ctree_binary", "ctree_survival")) sb_abort(
    "This diagnostic requires a single conditional inference tree.",
    "stormbreaker_unsupported_diagnostic")
  sb_party_tree_structure(fit)
}

#' @rdname sb_forest_diagnostics
#' @export
sb_evtree_structure <- function(fit) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (fit$pipeline$learner != "evtree_binary") sb_abort(
    "This diagnostic requires an evolutionary classification tree.", "stormbreaker_unsupported_diagnostic")
  sb_party_tree_structure(fit)
}

sb_party_tree_structure <- function(fit) {
  e <- unserialize(fit$learner_state$engine_raw)
  sb_party_native_structure(e, fit$pipeline$learner, fit$pipeline$id)
}

sb_party_native_structure <- function(e, learner, pipeline_id = NULL) {
  evolutionary <- learner %in% c("evtree_binary", "evtree_multiclass")
  if (!inherits(e, "constparty")) sb_abort("The original conditional tree is missing.",
    "stormbreaker_bad_learner_state")
  frame <- partykit::data_party(e)
  label_frame <- if (evolutionary) e$data else frame
  membership <- frame[["(fitted)"]]
  weights <- frame[["(weights)"]]
  if (is.null(membership) || is.null(weights) || length(membership) != length(weights)) sb_abort(
    "Original tree membership or weights are missing.", "stormbreaker_bad_learner_state")
  rows <- splits <- information <- list()
  leaf_position <- 0L
  visit <- function(node, parent = NA_integer_, depth = 0L, rule = "Root") {
    id <- partykit::id_node(node)
    key <- as.character(id)
    kids <- partykit::kids_node(node)
    split <- partykit::split_node(node)
    terminal <- !length(kids)
    rules <- if (terminal) NULL else partykit::character_split(split, data = label_frame, digits = 7)
    if (!terminal && length(rules$levels) != length(kids)) sb_abort(
      "Native split labels do not match child branches.", "stormbreaker_bad_learner_state")
    rows[[key]] <<- data.frame(node = id, parent = parent, depth = depth,
      terminal = terminal, incoming_rule = rule,
      split_feature = if (terminal) NA_character_ else rules$name,
      n = NA_integer_, weight = NA_real_, x = NA_real_, y = -depth)
    splits[key] <<- list(split)
    information[key] <<- list(partykit::info_node(node))
    if (terminal) {
      leaf_position <<- leaf_position + 1L
      descendants <- id
      position <- leaf_position
    } else {
      descendants <- unlist(lapply(seq_along(kids), function(j) {
        visit(kids[[j]], id, depth + 1L, paste(rules$name, rules$levels[j]))
      }), use.names = FALSE)
      positions <- vapply(kids, function(kid) rows[[as.character(partykit::id_node(kid))]]$x, 0)
      position <- mean(range(positions))
    }
    included <- membership %in% descendants
    rows[[key]]$n <<- sum(included)
    rows[[key]]$weight <<- sum(weights[included])
    rows[[key]]$x <<- position
    descendants
  }
  leaves <- visit(partykit::node_party(e))
  if (!all(membership %in% leaves)) sb_abort("Native membership contains unknown terminal nodes.",
    "stormbreaker_bad_learner_state")
  tab <- do.call(rbind, unname(rows))
  rownames(tab) <- NULL
  structure(list(type = "ctree_structure", table = tab, native_splits = splits,
    native_node_information = information, learner = learner,
    pipeline_id = pipeline_id,
    source = if (evolutionary) "original evolutionary tree and retained training memberships" else
      "all original conditional inference tree nodes and saved training memberships",
    interpretation = paste("Counts are original training-node memberships, not validation performance.",
      "Exact native split objects are retained separately from rounded display rules.",
      "Only originally retained node statistics are available; missing tests are not recomputed.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_ctree_structure <- function(x) {
  d <- x$table
  parent <- match(d$parent, d$node)
  edges <- data.frame(x = d$x[parent], y = d$y[parent] - .22,
    xend = d$x, yend = d$y + .22)
  edges <- edges[!is.na(parent), , drop = FALSE]
  d$label <- vapply(seq_len(nrow(d)), function(i) paste0("Node ", d$node[i],
    " | n = ", d$n[i], "\n",
    paste(strwrap(d$incoming_rule[i], width = 28), collapse = "\n")), "")
  y <- xend <- yend <- label <- terminal <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_segment(data = edges, ggplot2::aes(xend = xend, yend = yend), colour = "grey60") +
    ggplot2::geom_label(ggplot2::aes(label = label, fill = terminal), size = 3.1,
      linewidth = .2, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#EDF3F7", "TRUE" = "#D7EBE0")) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = .8)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = .5)) +
    ggplot2::labs(title = if (x$learner == "evtree_binary") "Evolutionary classification tree" else
      if (x$learner == "tree_binary") "Classification tree" else
      "Conditional inference tree", x = NULL, y = NULL) +
    sb_viz_theme() + ggplot2::theme(axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  attr(p, "stormbreaker_suggested_size") <- c(width = max(9, 2.3 * sum(d$terminal)),
    height = max(5, 1.6 * (max(d$depth) + 1)))
  p
}

# ---- algorithm-rpart-tree ----

#' @rdname sb_forest_diagnostics
#' @export
sb_tree_structure <- function(fit) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!fit$pipeline$learner %in% c("decision_tree_binary", "rpart_survival")) sb_abort(
    "This adapter requires a single rpart tree, not an ensemble or another tree implementation.",
    "stormbreaker_unsupported_diagnostic")
  if (is.null(fit$diagnostic_artifacts)) sb_abort(
    "The original tree was not captured; enable capture_diagnostics when fitting.",
    "stormbreaker_missing_diagnostic_artifacts")
  engine <- fit$diagnostic_artifacts$learner$engine
  if (!inherits(engine, "rpart")) sb_abort("The captured rpart engine is missing.",
    "stormbreaker_invalid_diagnostic")
  frame <- engine$frame
  nodes <- as.numeric(rownames(frame))
  if (!length(nodes) || any(!is.finite(nodes)) || anyDuplicated(nodes) ||
      nodes[1] != 1 || any(nodes < 1) || any(nodes != floor(nodes))) sb_abort(
    "Native tree node identifiers are invalid.", "stormbreaker_invalid_diagnostic")
  leaf <- frame$var == "<leaf>"
  if (any(!leaf & (!nodes * 2 %in% nodes | !(nodes * 2 + 1) %in% nodes))) sb_abort(
    "Native tree contains an incomplete branch.", "stormbreaker_invalid_diagnostic")
  tab <- data.frame(node = nodes, parent = ifelse(nodes == 1, NA_real_, floor(nodes / 2)),
    depth = floor(log2(nodes)), terminal = leaf, n = frame$n,
    split_feature = as.character(frame$var), incoming_rule = labels(engine, digits = 7, minlength = 0),
    native_value = frame$yval)
  if (identical(engine$method, "class")) {
    levels <- attr(engine, "ylevels")
    if (!identical(as.character(levels), c("0", "1")) || ncol(frame$yval2) != 6L) sb_abort(
      "The tree's binary class mapping is not the expected 0/1 mapping.", "stormbreaker_invalid_diagnostic")
    tab$prediction <- frame$yval2[, 5L]
    tab$class_0_n <- frame$yval2[, 2L]
    tab$class_1_n <- frame$yval2[, 3L]
    scale <- "Training-node P(class 1)"
  } else {
    tab$prediction <- frame$yval
    scale <- "Native relative event rate"
  }
  xpos <- rep(NA_real_, length(nodes))
  xpos[leaf] <- seq_len(sum(leaf))
  position <- function(node) {
    i <- match(node, nodes)
    if (is.finite(xpos[i])) return(xpos[i])
    xpos[i] <<- mean(c(position(2 * node), position(2 * node + 1)))
    xpos[i]
  }
  position(1)
  tab$x <- xpos
  tab$y <- -tab$depth
  structure(list(type = "tree_structure", table = tab,
    native_splits = engine$splits, native_categorical_splits = engine$csplit,
    complexity = engine$cptable, prediction_scale = scale,
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    source = "all nodes of original rpart fit; no pruning or refitting",
    suggested_width = max(9, 1.6 * sum(leaf)),
    suggested_height = max(5, 1.2 * (max(tab$depth) + 1))),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_rpart_tree <- function(x) {
  d <- x$table
  parent <- match(d$parent, d$node)
  edges <- data.frame(x = d$x[parent], y = d$y[parent], xend = d$x, yend = d$y)
  edges <- edges[!is.na(parent), , drop = FALSE]
  edges$y <- edges$y - .33
  edges$yend <- edges$yend + .33
  d$label <- vapply(seq_len(nrow(d)), function(i) {
    if ("display_label" %in% names(d)) return(d$display_label[i])
    rule <- paste(strwrap(d$incoming_rule[i], width = 24), collapse = "\n")
    paste0("Node ", d$node[i], " | n = ", d$n[i], "\n", rule, "\n",
      if (x$learner == "decision_tree_binary") "P(1) = " else "Relative rate = ",
      format(signif(d$prediction[i], 3), trim = TRUE))
  }, "")
  node_x <- node_y <- y <- xend <- yend <- label <- terminal <- NULL
  names(d)[names(d) == "x"] <- "node_x"
  names(d)[names(d) == "y"] <- "node_y"
  p <- ggplot2::ggplot(d, ggplot2::aes(node_x, node_y)) +
    ggplot2::geom_segment(data = edges, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE, colour = "#909090") +
    ggplot2::geom_label(ggplot2::aes(label = label, fill = terminal), size = 3,
      linewidth = .2, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#FFFFFF", "TRUE" = "#DDEBF3")) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = .75)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = .7)) +
    ggplot2::labs(title = "Single-tree decision structure",
      subtitle = x$prediction_scale,
      caption = "Every original node shown; a node's rule describes its incoming branch. Terminal nodes are shaded.\nDisplay thresholds are rounded; exact native splits remain in the diagnostic. No pruning or refitting.") +
    sb_viz_theme() + ggplot2::theme(axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank())
  attr(p, "stormbreaker_suggested_size") <- c(width = x$suggested_width, height = x$suggested_height)
  p
}

# ---- diagnostics-forests ----

# ---- algorithm-aorsf ----

#' @rdname sb_forest_diagnostics
#' @export
sb_aorsf_diagnostics <- function(fit, type = "tree_size") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% c("aorsf_binary", "aorsf_survival")) sb_abort(
    "This diagnostic requires an aorsf model.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("tree_size", "oob_score")) sb_abort(
    "Choose tree_size or oob_score.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  if (type == "tree_size") {
    leaves <- vapply(e$forest$leaf_summary, function(z) sum(z != 0), integer(1))
    if (length(leaves) != e$n_tree || anyNA(leaves) || any(leaves < 1)) sb_abort(
      "Original oblique forest leaf summaries are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(tree = seq_along(leaves), terminal_nodes = leaves)
    measure <- "Terminal nodes per tree"
  } else {
    measure <- e$eval_oobag$stat_type
    expected <- if (s$learner == "aorsf_binary") "AUC-ROC" else "Harrell's C-index"
    z <- e$eval_oobag$stat_values
    steps <- unique(c(seq.int(e$oobag_eval_every, e$n_tree, by = e$oobag_eval_every), e$n_tree))
    if (!identical(measure, expected) || !is.matrix(z) || ncol(z) != 1L ||
        nrow(z) != length(steps) || any(is.infinite(z)) || any(z < 0 | z > 1, na.rm = TRUE)) sb_abort(
      "Saved native OOB score or evaluation schedule is invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(trees = steps, score = as.numeric(z), available = is.finite(as.numeric(z)))
  }
  structure(list(type = paste0("aorsf_", type), table = tab, measure = measure,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "saved native aorsf leaf summaries and OOB evaluations",
    interpretation = paste("Counts use original native leaf summaries. OOB scores retain their original metric and evaluation schedule.",
      "AUC and concordance are scores, not error rates. A final-only evaluation remains one point, not an invented learning curve.",
      "OOB is not validation of prior preprocessing/selection or external performance; no importance, intervals or new predictions are computed.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_aorsf_diagnostic <- function(x) {
  terminal_nodes <- trees <- score <- NULL
  d <- x$table
  if (x$type == "aorsf_tree_size") {
    counts <- as.data.frame(table(d$terminal_nodes))
    names(counts) <- c("terminal_nodes", "trees")
    counts$terminal_nodes <- as.numeric(as.character(counts$terminal_nodes))
    p <- ggplot2::ggplot(counts, ggplot2::aes(terminal_nodes, trees)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::scale_y_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::labs(title = "Oblique forest tree sizes", x = x$measure, y = "Number of trees")
  } else p <- ggplot2::ggplot(d[d$available, , drop = FALSE], ggplot2::aes(trees, score)) +
    ggplot2::geom_point(colour = "#2878A0", size = 3) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_x_continuous(limits = c(0, max(d$trees))) +
    ggplot2::labs(title = "Recorded oblique forest OOB score", x = "Number of trees", y = x$measure)
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}

# ---- algorithm-bagging ----

#' @rdname sb_forest_diagnostics
#' @export
sb_bagging_diagnostics <- function(fit, type = "split_usage") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!identical(s$learner, "bagging_binary")) sb_abort(
    "This diagnostic requires ipred bagging.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("split_usage", "tree_size", "importance")) sb_abort(
    "Choose split_usage, tree_size or importance.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  if (!inherits(e, "classbagg") || length(e$mtrees) != s$params$num_trees) sb_abort(
    "Saved bagged trees are incomplete.", "stormbreaker_bad_learner_state")
  counts <- stats::setNames(integer(length(s$features)), s$features)
  rows <- lapply(seq_along(e$mtrees), function(i) {
    tree <- e$mtrees[[i]]$btree
    frame <- tree$frame
    if (!inherits(tree, "rpart") || !is.data.frame(frame) || !nrow(frame)) sb_abort(
      "An original bagged rpart tree is invalid.", "stormbreaker_invalid_diagnostic")
    vars <- frame$var[frame$var != "<leaf>"]
    if (any(!vars %in% s$features)) sb_abort("A bagged split has an unknown feature.", "stormbreaker_invalid_diagnostic")
    counts <<- counts + as.integer(table(factor(vars, levels = s$features)))
    values <- stats::setNames(numeric(length(s$features)), s$features)
    native <- tree$variable.importance
    if (length(native)) {
      if (any(!names(native) %in% s$features) || any(!is.finite(native))) sb_abort(
        "Native bagged-tree importance is invalid.", "stormbreaker_invalid_diagnostic")
      values[names(native)] <- native
    }
    list(size = data.frame(tree = i, nodes = nrow(frame), terminal_nodes = sum(frame$var == "<leaf>")),
      importance = data.frame(tree = i, feature = s$features, importance = unname(values)))
  })
  tab <- if (type == "split_usage") data.frame(feature = s$features, splits = unname(counts)) else
    do.call(rbind, lapply(rows, `[[`, if (type == "tree_size") "size" else "importance"))
  structure(list(type = paste0("bagging_", type), table = tab, learner = s$learner,
    pipeline_id = fit$pipeline$id, source = "original ipred bags and native rpart frames/importance",
    interpretation = paste("All original bags and features retained. Split usage counts primary splits only.",
      "Native rpart importance may credit surrogate splits; its across-bag mean is not permutation importance or a causal effect.",
      "OOB was disabled. No OOB history, new predictions, confidence intervals or surrogate model is reconstructed.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_bagging_diagnostic <- function(x) {
  feature <- value <- terminal_nodes <- trees <- NULL
  d <- x$table
  if (x$type != "bagging_tree_size") {
    if (x$type == "bagging_importance") {
      d <- stats::aggregate(importance ~ feature, data = d, FUN = mean)
      names(d)[2] <- "value"
    } else names(d)[2] <- "value"
    d$feature <- factor(d$feature, levels = d$feature[order(d$value)])
    p <- ggplot2::ggplot(d, ggplot2::aes(value, feature)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::labs(title = if (x$type == "bagging_importance") "Bagged-tree native importance" else "Bagged-tree split usage",
        x = if (x$type == "bagging_importance") "Mean native rpart importance" else "Number of primary splits", y = NULL)
  } else {
    counts <- as.data.frame(table(d$terminal_nodes))
    names(counts) <- c("terminal_nodes", "trees")
    counts$terminal_nodes <- as.numeric(as.character(counts$terminal_nodes))
    p <- ggplot2::ggplot(counts, ggplot2::aes(terminal_nodes, trees)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::scale_y_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::labs(title = "Bagged-tree sizes", x = "Terminal nodes per tree", y = "Number of trees")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}

# ---- algorithm-cforest ----

#' @rdname sb_forest_diagnostics
#' @export
sb_cforest_diagnostics <- function(fit, type = "split_usage") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% c("cforest_binary", "cforest_survival")) sb_abort(
    "This diagnostic requires cforest.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("split_usage", "tree_size")) sb_abort(
    "Choose split_usage or tree_size.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  if (!inherits(e, "cforest") || length(e$nodes) != s$params$num_trees) sb_abort(
    "Original conditional forest nodes are incomplete.", "stormbreaker_bad_learner_state")
  counts <- stats::setNames(integer(length(s$features)), s$features)
  visit <- function(node) {
    kids <- partykit::kids_node(node)
    if (!length(kids)) return(c(nodes = 1L, terminal_nodes = 1L))
    index <- partykit::varid_split(partykit::split_node(node))
    feature <- names(e$data)[index]
    if (length(feature) != 1L || is.na(feature) || !feature %in% s$features) sb_abort(
      "Original conditional forest split mapping is invalid.", "stormbreaker_invalid_diagnostic")
    counts[feature] <<- counts[feature] + 1L
    c(nodes = 1L, terminal_nodes = 0L) + Reduce(`+`, lapply(kids, visit))
  }
  sizes <- t(vapply(e$nodes, visit, c(nodes = 0L, terminal_nodes = 0L)))
  tab <- if (type == "split_usage") data.frame(feature = s$features, splits = unname(counts)) else
    data.frame(tree = seq_len(nrow(sizes)), nodes = sizes[, "nodes"], terminal_nodes = sizes[, "terminal_nodes"])
  structure(list(type = paste0("cforest_", type), table = tab, learner = s$learner,
    pipeline_id = fit$pipeline$id, source = "original saved partykit cforest nodes and model-frame mapping",
    interpretation = paste("Every original tree and feature is retained. Counts are structural use, not conditional permutation importance.",
      "No new significance tests, OOB evaluation, permutations, surrogate trees or predictions are computed.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_cforest_diagnostic <- function(x) {
  feature <- splits <- terminal_nodes <- trees <- NULL
  d <- x$table
  if (x$type == "cforest_split_usage") {
    d$feature <- factor(d$feature, levels = d$feature[order(d$splits)])
    p <- ggplot2::ggplot(d, ggplot2::aes(splits, feature)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::labs(title = "Conditional forest split usage", x = "Number of original splits", y = NULL)
  } else {
    counts <- as.data.frame(table(d$terminal_nodes))
    names(counts) <- c("terminal_nodes", "trees")
    counts$terminal_nodes <- as.numeric(as.character(counts$terminal_nodes))
    p <- ggplot2::ggplot(counts, ggplot2::aes(terminal_nodes, trees)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::scale_y_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::labs(title = "Conditional forest tree sizes", x = "Terminal nodes per tree", y = "Number of trees")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}

# ---- algorithm-ferns ----

#' @rdname sb_forest_diagnostics
#' @export
sb_ferns_diagnostics <- function(fit, type = "split_usage") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (s$learner != "random_ferns_binary" || length(type) != 1L ||
      !is.character(type) || is.na(type) || !type %in% c("split_usage", "thresholds")) {
    sb_abort("Choose split_usage or thresholds for random ferns.", "stormbreaker_unsupported_diagnostic")
  }
  e <- unserialize(s$engine_raw)
  depth <- unname(e$parameters["depth"])
  nferns <- unname(e$parameters["ferns"])
  idx <- e$model$splitAttIdxs
  features <- names(e$isStruct$predictorLevels)
  if (!inherits(e, "rFerns") || !identical(features, s$features) ||
      length(idx) != depth * nferns || anyNA(idx) ||
      any(idx < 0L | idx >= length(features))) {
    sb_abort("Original random-fern splits are invalid.", "stormbreaker_invalid_diagnostic")
  }
  splits <- data.frame(fern = rep(seq_len(nferns), each = depth),
    split = rep(seq_len(depth), nferns), feature = features[idx + 1L])
  if (type == "thresholds") {
    if (any(lengths(e$isStruct$predictorLevels) > 0L) ||
        any(e$isStruct$integerPredictors)) {
      sb_abort("Threshold plots require original real-valued predictors; integer and categorical encodings are not converted.",
        "stormbreaker_unavailable_diagnostic")
    }
    if (length(e$model$threReal) != nrow(splits) || any(!is.finite(e$model$threReal))) {
      sb_abort("Original real-valued thresholds are invalid.", "stormbreaker_invalid_diagnostic")
    }
    splits$threshold <- e$model$threReal
    tab <- splits
  } else {
    tab <- data.frame(feature = features, splits = tabulate(idx + 1L, nbins = length(features)))
  }
  structure(list(type = paste0("ferns_", type), table = tab, splits = splits,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original saved rFerns split indices and thresholds",
    interpretation = "Random split usage is not predictive importance. Thresholds describe fitted structure, not clinical cutoffs. No OOB performance is reconstructed."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_ferns_diagnostic <- function(x) {
  feature <- splits <- threshold <- NULL
  if (x$type == "ferns_split_usage") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(splits, feature)) +
      ggplot2::geom_col(fill = "#2878A0", width = .65) +
      ggplot2::labs(title = "Random-fern split usage", x = "Number of splits", y = NULL)
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(threshold, feature)) +
      ggplot2::geom_boxplot(fill = "#A9CEDF", width = .55, outlier.size = .8) +
      ggplot2::labs(title = "Random-fern split thresholds", x = "Original threshold", y = NULL)
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 10, height = max(4, length(unique(x$table$feature)) * .3 + 2))
  p
}

# ---- algorithm-rotation-forest ----

#' @rdname sb_forest_diagnostics
#' @export
sb_rotation_forest_diagnostics <- function(fit,type="loadings",newdata=NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if(s$learner!="rotation_forest_binary" || !is.character(type) || length(type)!=1L ||
      is.na(type) || !type %in% c("loadings","tree_size","tree_predictions")) sb_abort(
    "Choose loadings, tree_size or tree_predictions for rotation forest.","stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  if(!inherits(e,"rotationForest") || !identical(e$columnnames,s$features) || length(e$models)!=length(e$loadings)) sb_abort(
    "Original rotation forest is invalid.","stormbreaker_invalid_diagnostic")
  if(type=="loadings") {
    tab <- do.call(rbind,lapply(seq_along(e$loadings),function(i) {
      m <- as.matrix(e$loadings[[i]])
      if(!identical(dim(m),rep(length(s$features),2)) || any(!is.finite(m))) sb_abort(
        "An original rotation matrix is invalid.","stormbreaker_invalid_diagnostic")
      data.frame(tree=i,feature=rep(s$features,ncol(m)),axis=rep(seq_len(ncol(m)),each=nrow(m)),
        native_axis=rep(colnames(m),each=nrow(m)),loading=as.numeric(m))
    }))
  } else if(type=="tree_size") {
    tab <- do.call(rbind,lapply(seq_along(e$models),function(i) {
      f <- e$models[[i]]$frame
      data.frame(tree=i,nodes=nrow(f),terminal_nodes=sum(f$var=="<leaf>"))
    }))
  } else {
    if(is.null(newdata)) sb_abort("Supply explicit queries for rotation-forest tree predictions.","stormbreaker_unavailable_diagnostic")
    x <- sb_select_final_feature_frame(newdata,s$features,"Query data","query fitted features")
    if(!nrow(x)) sb_abort("Queries must contain at least one row.","stormbreaker_invalid_diagnostic")
    native_predict <- get("predict.rotationForest",asNamespace("rotationForest"))
    m <- native_predict(e,x,all=TRUE)
    if(!is.matrix(m) || nrow(m)!=nrow(x) || ncol(m)!=length(e$models) || any(!is.finite(m)) || any(m<0 | m>1)) sb_abort(
      "Native per-tree probabilities are invalid.","stormbreaker_invalid_diagnostic")
    tab <- data.frame(query=rep(seq_len(nrow(m)),ncol(m)),tree=rep(seq_len(ncol(m)),each=nrow(m)),probability=as.numeric(m))
  }
  structure(list(type=paste0("rotation_forest_",type),table=tab,learner=s$learner,pipeline_id=fit$pipeline$id,
    source="Original saved rotation matrices, rpart frames and native per-tree probabilities",
    interpretation="Rotated axes are not original features. Loading signs are basis conventions, not feature effects or importance. Tree variation is not a confidence interval. No PCA refit, pruning or validation."),
    class="stormbreaker_algorithm_diagnostic")
}

sb_plot_rotation_forest_diagnostic <- function(x) {
  axis <- feature <- loading <- tree <- terminal_nodes <- probability <- query <- NULL
  if(x$type=="rotation_forest_loadings") {
    p <- ggplot2::ggplot(x$table,ggplot2::aes(factor(axis),feature,fill=loading)) +
      ggplot2::geom_tile(colour="white",linewidth=.2) + ggplot2::facet_wrap(~tree,ncol=5,labeller=ggplot2::label_both) +
      ggplot2::scale_fill_gradient2(low="#2878A0",mid="white",high="#B86236",midpoint=0,
        guide=ggplot2::guide_colourbar(display="rectangles")) +
      ggplot2::labs(title="Original rotation matrices",x="Rotated axis",y=NULL,fill="Loading")
    size <- c(width=14,height=7)
  } else if(x$type=="rotation_forest_tree_size") {
    p <- ggplot2::ggplot(x$table,ggplot2::aes(tree,terminal_nodes)) + ggplot2::geom_point(colour="#2878A0") +
      ggplot2::scale_x_continuous(breaks=x$table$tree) +
      ggplot2::labs(title="Rotation-forest tree sizes",x="Tree",y="Terminal nodes")
    size <- c(width=11,height=5)
  } else {
    p <- ggplot2::ggplot(x$table,ggplot2::aes(tree,probability)) + ggplot2::geom_point(colour="#2878A0") +
      ggplot2::facet_wrap(~query,labeller=ggplot2::label_both) + ggplot2::scale_y_continuous(limits=c(0,1)) +
      ggplot2::scale_x_continuous(breaks=unique(x$table$tree)) +
      ggplot2::labs(title="Native per-tree predictions",x="Tree",y="Class-1 probability")
    size <- c(width=11,height=5)
  }
  p <- p + sb_viz_theme(); attr(p,"stormbreaker_suggested_size") <- size; p
}
