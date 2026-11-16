# ---- multiclass-boosting ----

sb_multiclass_boosting_catalog <- function() {
  data.frame(method = c("xgboost", "xgboost_dart", "xgboost_linear", "lightgbm",
    "lightgbm_dart", "lightgbm_goss", "lightgbm_rf", "catboost"),
    label = c("XGBoost", "XGBoost DART", "XGBoost linear", "LightGBM", "LightGBM DART",
      "LightGBM GOSS", "LightGBM RF", "CatBoost"), supported = TRUE,
    engine = c(rep("xgboost", 3), rep("lightgbm", 4), "catboost"),
    feature_result = c("ranking", "ranking", "model_features", rep("ranking", 5)),
    selection_rule = "No standalone gene-selection rule")
}

sb_fit_multiclass_boosting <- function(fit, x, y) {
  method <- fit$method
  catalog <- sb_multiclass_boosting_catalog()
  package <- catalog$engine[match(method, catalog$method)]
  if (!requireNamespace(package, quietly = TRUE)) sb_abort(
    paste("Install", package, "for this method."), "stormbreaker_missing_dependency")
  fit$package <- package
  fit$feature_codes <- paste0("feature", seq_along(fit$features))
  x <- as.matrix(x)
  colnames(x) <- fit$feature_codes
  target <- as.integer(y) - 1L
  fit$importance <- NULL
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  withCallingHandlers(withr::with_seed(fit$seed, {
    if (package == "xgboost") {
      id <- switch(method, xgboost = "xgboost_binary",
        xgboost_dart = "xgboost_binary_dart", xgboost_linear = "xgboost_binary_linear")
      params <- sb_xgboost_params(id)
      params$objective <- "multi:softprob"
      params$num_class <- length(fit$codes)
      params$eval_metric <- "mlogloss"
      params$seed <- fit$seed
      booster <- xgboost::xgb.train(params = params[setdiff(names(params), c("nrounds", "seed"))],
        data = xgboost::xgb.DMatrix(x, label = target), nrounds = params$nrounds, verbose = 0)
      fit$engine <- xgboost::xgb.save.raw(booster)
      if (method != "xgboost_linear") {
        importance <- xgboost::xgb.importance(feature_names = fit$feature_codes, model = booster)
        fit$importance <- stats::setNames(importance$Gain, importance$Feature)
        fit$importance_label <- "Gain share"
      }
    } else if (package == "lightgbm") {
      id <- switch(method, lightgbm = "lightgbm_binary", lightgbm_dart = "lightgbm_binary_dart",
        lightgbm_goss = "lightgbm_binary_goss", lightgbm_rf = "lightgbm_binary_rf")
      params <- sb_lightgbm_binary_params(learner = id, seed = fit$seed)
      params$objective <- "multiclass"
      params$num_class <- length(fit$codes)
      params$metric <- "multi_logloss"
      training_params <- params[setdiff(names(params), "nrounds")]
      # Native alias avoids BiocGenerics format.list truncating multi_logloss.
      training_params$metric <- "multiclass"
      booster <- lightgbm::lgb.train(params = training_params,
        data = lightgbm::lgb.Dataset(x, label = target), nrounds = params$nrounds, verbose = -1L)
      fit$engine <- booster$save_model_to_string()
      importance <- lightgbm::lgb.importance(booster, percentage = FALSE)
      fit$importance <- stats::setNames(importance$Gain, importance$Feature)
      fit$importance_label <- "Total split gain"
    } else {
      params <- sb_catboost_binary_params(seed = fit$seed)
      params$loss_function <- "MultiClass"
      params$eval_metric <- "MultiClass"
      params$classes_count <- length(fit$codes)
      params$boost_from_average <- FALSE
      booster <- catboost::catboost.train(catboost::catboost.load_pool(x, label = target),
        test_pool = NULL, params = params)
      fit$engine <- serialize(booster, NULL)
      fit$importance <- stats::setNames(catboost::catboost.get_feature_importance(booster,
        type = "PredictionValuesChange", thread_count = 1L), fit$feature_codes)
      fit$importance_label <- "Prediction value change"
    }
    fit$parameters <- params
  }), warning = record_warning)
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_predict_multiclass_boosting <- function(fit, x) {
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "to predict with this saved fit."), "stormbreaker_missing_dependency")
  x <- as.matrix(x)
  colnames(x) <- fit$feature_codes
  if (fit$package == "xgboost") {
    booster <- xgboost::xgb.load.raw(fit$engine)
    p <- stats::predict(booster, x, reshape = TRUE)
  } else if (fit$package == "lightgbm") {
    booster <- lightgbm::lgb.load(model_str = fit$engine)
    p <- stats::predict(booster, x)
  } else {
    p <- catboost::catboost.predict(unserialize(fit$engine), catboost::catboost.load_pool(x),
      prediction_type = "Probability", thread_count = 1L)
  }
  if (length(p) != nrow(x) * length(fit$codes)) sb_abort(
    "The native booster returned an unexpected probability shape.", "stormbreaker_bad_multiclass_probabilities")
  matrix(p, nrow = nrow(x), ncol = length(fit$codes))
}

sb_multiclass_boosting_evidence <- function(fit) {
  tab <- data.frame(feature = fit$features, selected = NA)
  if (!is.null(fit$importance)) {
    tab$importance <- 0
    index <- match(fit$feature_codes, names(fit$importance))
    tab$importance[!is.na(index)] <- unname(fit$importance[index[!is.na(index)]])
    tab$rank <- rank(-tab$importance, ties.method = "min")
  }
  list(features = tab, by_class = NULL)
}

sb_multiclass_boosting_plots <- function(fit) {
  if (is.null(fit$importance)) return(list(plots = list(), tables = list()))
  tab <- sb_multiclass_boosting_evidence(fit)$features
  tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$importance)])
  catalog <- sb_multiclass_boosting_catalog()
  label <- catalog$label[match(fit$method, catalog$method)]
  feature <- importance <- NULL
  plot <- ggplot2::ggplot(tab, ggplot2::aes(importance, feature)) +
    ggplot2::geom_col(fill = "#2878A0") + sb_viz_theme() +
    ggplot2::labs(title = paste0(label, ": Importance"), x = fit$importance_label, y = NULL)
  list(plots = list(importance = plot), tables = list(importance = tab))
}

# ---- multiclass-trees ----

# ---- multiclass-forest-structure-plots ----

sb_multiclass_forest_structure_plots <- function(fit) {
  m <- fit$engine
  counts <- stats::setNames(integer(length(fit$features)), fit$feature_codes)
  if (fit$method == "cforest") {
    visit <- function(node) {
      kids <- partykit::kids_node(node)
      if (!length(kids)) return(1L)
      name <- names(m$data)[partykit::varid_split(partykit::split_node(node))]
      counts[name] <<- counts[name] + 1L
      sum(vapply(kids, visit, 1L))
    }
    sizes <- vapply(m$nodes, visit, 1L)
  } else if (fit$method == "bagging") {
    sizes <- vapply(m$mtrees, function(tree) {
      vars <- as.character(tree$btree$frame$var)
      counts <<- counts + tabulate(match(vars[vars != "<leaf>"], fit$feature_codes), length(counts))
      sum(vars == "<leaf>")
    }, 1L)
  } else if (fit$method == "rlt") {
    sizes <- vapply(m$FittedTrees, function(tree) {
      internal <- tree["NodeType", ] == 1
      counts <<- counts + tabulate(as.integer(tree["SplitVar1", internal]), length(counts))
      sum(!internal)
    }, 1L)
  } else if (fit$method == "aorsf") {
    sizes <- vapply(m$forest$leaf_summary, function(x) sum(x != 0), 1L)
  } else sizes <- lengths(m$trees)
  d <- data.frame(tree = seq_along(sizes), nodes = as.integer(sizes))
  frequencies <- as.data.frame(table(d$nodes))
  names(frequencies) <- c("nodes", "trees")
  frequencies$nodes <- as.integer(as.character(frequencies$nodes))
  nodes <- trees <- NULL
  integer_breaks <- function(limits) {
    ticks <- pretty(limits, n = 5)
    ticks[ticks == floor(ticks) & ticks >= 0]
  }
  p <- ggplot2::ggplot(frequencies, ggplot2::aes(nodes, trees)) + ggplot2::geom_col(fill = "#2878A0") +
    ggplot2::scale_x_continuous(breaks = integer_breaks) +
    ggplot2::scale_y_continuous(breaks = integer_breaks) +
    sb_viz_theme() + ggplot2::labs(title = "Tree-size distribution",
      x = if (fit$method == "wsrf") "Total nodes per tree" else "Terminal nodes per tree", y = "Trees",
      caption = "Every original fitted tree is counted once; node counts describe structural complexity, not predictive performance. No pruning or surrogate tree is used.")
  result <- list(plots = list(tree_size = p), tables = list(tree_size = d))
  if (fit$method %in% c("cforest", "bagging", "rlt")) {
    d <- data.frame(feature = fit$features, splits = unname(counts))
    d$feature <- factor(d$feature, levels = d$feature[order(d$splits)])
    feature <- splits <- NULL
    result$plots$split_usage <- ggplot2::ggplot(d, ggplot2::aes(splits, feature)) +
      ggplot2::geom_col(fill = "#2878A0") + sb_viz_theme() +
      ggplot2::scale_x_continuous(breaks = integer_breaks) +
      ggplot2::labs(title = "Split usage", x = "Splits", y = NULL,
        caption = "Counts of original primary axis-aligned splits, including zero-use features; not impurity, permutation importance or causal effects.")
    result$tables$split_usage <- d
  }
  if (fit$method == "wsrf") {
    d <- data.frame(tree = seq_along(m$treeOOBErrorRates), error = as.numeric(m$treeOOBErrorRates))
    tree <- error <- NULL
    result$plots$individual_tree_oob <- ggplot2::ggplot(d, ggplot2::aes(tree, error)) +
      ggplot2::geom_point(colour = "#2878A0", size = 1) + sb_viz_theme() +
      ggplot2::labs(title = "Individual-tree OOB error", x = "Tree", y = "Error rate",
        caption = "Saved individual-tree out-of-bag errors, not a cumulative ensemble learning curve or held-out validation.")
    result$tables$individual_tree_oob <- d
  }
  result
}

# ---- multiclass-tree-native-plots ----

sb_multiclass_tree_native_plots <- function(fit) {
  forest <- fit$method == "rrf"
  native <- if (forest) fit$engine$importance else C50::C5imp(fit$engine, metric = "usage", pct = TRUE)
  d <- data.frame(feature = fit$features[match(rownames(native), fit$feature_codes)],
    importance = as.numeric(native[, 1L]))
  d$feature <- factor(d$feature, levels = d$feature[order(d$importance)])
  feature <- importance <- NULL
  plot <- ggplot2::ggplot(d, ggplot2::aes(importance, feature)) +
    ggplot2::geom_col(fill = "#2878A0") + sb_viz_theme() +
    ggplot2::labs(title = if (forest) "Regularized forest importance" else "C5.0 feature usage",
      x = if (forest) "Mean decrease in Gini" else "Usage (%)", y = NULL,
      caption = if (forest) "Original impurity importance, including zero-valued features; not permutation importance." else
        "Original C5.0 usage percentages; these are not mutually exclusive and need not sum to 100%.")
  result <- list(plots = list(importance = plot), tables = list(importance = d))
  if (forest) {
    errors <- fit$engine$err.rate
    names <- colnames(errors)
    at <- match(names, fit$codes)
    names[!is.na(at)] <- fit$classes[at[!is.na(at)]]
    d <- data.frame(trees = rep(seq_len(nrow(errors)), ncol(errors)),
      series = rep(names, each = nrow(errors)), error = as.vector(errors))
    trees <- error <- series <- NULL
    result$plots$oob_error <- ggplot2::ggplot(d, ggplot2::aes(trees, error, colour = series)) +
      ggplot2::geom_line() + sb_viz_theme() + ggplot2::labs(title = "Out-of-bag error",
        x = "Trees", y = "Error rate", colour = NULL,
        caption = "Original cumulative out-of-bag classification errors, overall and by class; not external validation.")
    result$tables$oob_error <- d
  }
  result
}

# ---- multiclass-tree-structure ----

sb_multiclass_tree_structure <- function(fit) {
  engine <- fit$engine
  frame <- engine$frame
  node <- as.numeric(rownames(frame))
  leaf <- frame$var == "<leaf>"
  rules <- if (inherits(engine, "tree")) labels(engine, pretty = FALSE) else
    labels(engine, digits = 5, minlength = 0)
  for (i in order(nchar(fit$feature_codes), decreasing = TRUE))
    rules <- gsub(fit$feature_codes[i], fit$features[i], rules, fixed = TRUE)
  codes <- attr(engine, "ylevels")
  predicted <- fit$classes[match(codes[frame$yval], fit$codes)]
  x <- rep(NA_real_, length(node))
  x[leaf] <- seq_len(sum(leaf))
  place <- function(n) {
    i <- match(n, node)
    if (is.na(x[i])) x[i] <<- mean(c(place(2 * n), place(2 * n + 1)))
    x[i]
  }
  place(1)
  tab <- data.frame(node = node, parent = ifelse(node == 1, NA_real_, floor(node / 2)),
    x = x, y = -floor(log2(node)), terminal = leaf, n = frame$n,
    incoming_rule = rules, prediction = predicted,
    display_label = paste0("n = ", frame$n, "\n", rules, "\nClass: ", predicted))
  list(table = tab, prediction_scale = "Training-node class",
    learner = if (inherits(engine, "tree")) "multiclass_tree" else "multiclass_rpart",
    suggested_width = max(11, 1.6 * sum(leaf)),
    suggested_height = max(7, 1.2 * (max(-tab$y) + 1)))
}

# ---- multiclass-trees ----

sb_multiclass_tree_catalog <- function() {
  data.frame(method = c("ranger", "ranger_extratrees", "c50", "c50_rules", "ctree",
    "cforest", "bagging", "rrf", "tree", "evtree"),
    label = c("Ranger forest", "Extremely randomized forest", "C5.0 trees", "C5.0 rules",
      "Conditional inference tree", "Conditional inference forest", "Bagged trees",
      "Regularized random forest", "Classification tree", "Evolutionary tree"),
    supported = TRUE, engine = c("ranger", "ranger", "C50", "C50", "partykit",
      "partykit", "ipred", "RRF", "tree", "evtree"),
    feature_result = c("ranking", "ranking", rep("model_features", 8)),
    selection_rule = "No standalone gene-selection rule")
}

sb_fit_multiclass_tree <- function(fit, x, y) {
  catalog <- sb_multiclass_tree_catalog()
  method <- fit$method
  package <- catalog$engine[match(method, catalog$method)]
  if (!requireNamespace(package, quietly = TRUE)) sb_abort(
    paste("Install", package, "for this method."), "stormbreaker_missing_dependency")
  fit$package <- package
  fit$feature_codes <- paste0("feature", seq_along(fit$features))
  names(x) <- fit$feature_codes
  d <- x
  d$.outcome <- y
  formula <- stats::as.formula(".outcome ~ .", env = baseenv())
  params <- if (method %in% c("ranger", "ranger_extratrees"))
    sb_ranger_binary_params(learner = paste0(method, "_binary")) else
    sb_classical_binary_params(paste0(method, "_binary"))
  params$seed <- fit$seed
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  fit$engine <- withCallingHandlers(withr::with_seed(fit$seed, {
    if (method %in% c("ranger", "ranger_extratrees")) {
      ranger::ranger(x = x, y = y, probability = TRUE, num.trees = params$num_trees,
        mtry = max(1L, floor(sqrt(ncol(x)))), min.node.size = params$min_node_size,
        splitrule = params$splitrule, importance = "impurity", replace = params$replace,
        sample.fraction = params$sample_fraction, num.threads = 1L,
        seed = fit$seed, oob.error = FALSE)
    } else if (method %in% c("c50", "c50_rules")) {
      C50::C5.0(x, y, trials = params$trials, rules = params$rules,
        control = C50::C5.0Control(winnow = FALSE, earlyStopping = FALSE, seed = fit$seed))
    } else if (method %in% c("ctree", "cforest")) {
      control <- partykit::ctree_control(teststat = params$teststat, testtype = params$testtype,
        mincriterion = params$mincriterion, minsplit = params$minsplit,
        minbucket = params$minbucket, maxdepth = params$maxdepth, saveinfo = params$saveinfo)
      if (method == "ctree") partykit::ctree(formula, data = d, control = control) else
        partykit::cforest(formula, data = d, ntree = params$num_trees,
          mtry = max(1L, ceiling(sqrt(ncol(x)))),
          perturb = list(replace = params$replace, fraction = params$sample_fraction),
          control = control, cores = 1L, trace = FALSE)
    } else if (method == "bagging") {
      ipred::bagging(formula, data = d, nbagg = params$num_trees, coob = FALSE,
        keepX = params$keep_x, control = rpart::rpart.control(cp = params$cp,
          minsplit = params$minsplit, minbucket = params$minbucket,
          maxdepth = params$maxdepth, xval = 0L))
    } else if (method == "rrf") {
      RRF::RRF(x, y, ntree = params$num_trees, mtry = max(1L, floor(sqrt(ncol(x)))),
        nodesize = params$min_node_size, replace = params$replace,
        coefReg = params$regularization_coefficient, flagReg = params$regularization_mode,
        importance = FALSE, keep.forest = TRUE)
    } else if (method == "tree") {
      tree::tree(formula, data = d, split = params$split,
        control = tree::tree.control(nrow(d), mincut = params$mincut,
          minsize = params$minsize, mindev = params$mindev),
        model = FALSE, x = FALSE, y = FALSE, wts = FALSE)
    } else {
      evtree::evtree(formula, data = d, control = evtree::evtree.control(
        minbucket = params$minbucket, minsplit = params$minsplit, maxdepth = params$maxdepth,
        niterations = params$iterations, ntrees = params$population_size,
        alpha = params$complexity_weight, seed = fit$seed))
    }
  }), warning = record_warning)
  fit$parameters <- params
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_predict_multiclass_tree <- function(fit, x) {
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "to predict with this saved fit."), "stormbreaker_missing_dependency")
  names(x) <- fit$feature_codes
  method <- fit$method
  probability <- if (method %in% c("ranger", "ranger_extratrees")) {
    stats::predict(fit$engine, data = x, num.threads = 1L)$predictions
  } else if (method == "cforest") {
    stats::predict(fit$engine, newdata = x, type = "prob", simplify = TRUE, OOB = FALSE)
  } else if (method == "tree") {
    stats::predict(fit$engine, newdata = x, type = "vector")
  } else stats::predict(fit$engine, newdata = x, type = "prob")
  if (is.list(probability)) probability <- do.call(rbind, probability)
  probability[, fit$codes, drop = FALSE]
}

sb_multiclass_tree_evidence <- function(fit) {
  tab <- data.frame(feature = fit$features, selected = NA)
  if (fit$method %in% c("ranger", "ranger_extratrees")) {
    tab$importance <- unname(fit$engine$variable.importance[fit$feature_codes])
    tab$rank <- rank(-tab$importance, ties.method = "min")
  }
  list(features = tab, by_class = NULL)
}

sb_multiclass_tree_plots <- function(fit) {
  if (fit$method == "tree") {
    structure <- sb_multiclass_tree_structure(fit)
    return(list(plots = list(tree_structure = sb_plot_rpart_tree(structure)),
      tables = list(tree_structure = structure$table)))
  }
  if (fit$method %in% c("ctree", "evtree")) {
    structure <- sb_party_native_structure(fit$engine, paste0(fit$method, "_multiclass"))
    for (i in order(nchar(fit$feature_codes), decreasing = TRUE)) {
      structure$table$incoming_rule <- gsub(fit$feature_codes[i], fit$features[i],
        structure$table$incoming_rule, fixed = TRUE)
      structure$table$split_feature <- gsub(fit$feature_codes[i], fit$features[i],
        structure$table$split_feature, fixed = TRUE)
    }
    p <- sb_plot_ctree_structure(structure) + ggplot2::labs(
      title = if (fit$method == "evtree") "Evolutionary tree" else "Conditional inference tree",
      caption = structure$interpretation)
    return(list(plots = list(tree_structure = p), tables = list(tree_structure = structure$table)))
  }
  if (fit$method %in% c("c50", "c50_rules", "rrf")) return(sb_multiclass_tree_native_plots(fit))
  if (!fit$method %in% c("ranger", "ranger_extratrees")) return(list(plots = list(), tables = list()))
  tab <- sb_multiclass_tree_evidence(fit)$features
  tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$importance)])
  feature <- importance <- NULL
  plot <- ggplot2::ggplot(tab, ggplot2::aes(importance, feature)) +
    ggplot2::geom_col(fill = "#2878A0") + sb_viz_theme() +
    ggplot2::labs(title = "Forest: Importance", x = "Impurity importance", y = NULL)
  list(plots = list(importance = plot), tables = list(importance = tab))
}
