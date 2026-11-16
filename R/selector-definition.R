# ---- selector-specification ----

# ---- doctor-selectors ----

sb_doctor_selector_menu <- function(task) {
  if (task == "binary") c(Lasso = "lasso_binary_nonzero", Stepglm = "stepglm_binary_aic",
    glmBoost = "glmboost_binary_nonzero", GamBoost = "gamboost_binary_nonzero", RF = "rf_binary_importance",
    Enet = "enet_binary_nonzero", PAMR = "pamr_binary_nonzero", PLSR = "plsrglm_binary_weights", GBM = "gbm_binary_importance", XGBoost = "xgboost_binary_importance", SVM_RFE = "svm_rfe_binary",
    Boruta = "boruta_binary_roughfix", varSelRF = "varselrf_binary_oob_1se",
    VSURF = "vsurf_binary_prediction", ABESS = "abess_binary_gic_nonzero",
    SIS = "sis_binary_nlogn", mRMR = "mrmr_binary_nlogn", JMI = "jmi_binary_nlogn",
    CMIM = "cmim_binary_nlogn", DISR = "disr_binary_nlogn", JMIM = "jmim_binary_nlogn",
    LPC = "lpc_binary_nlogn", Stabsel = "stabsel_binary_lasso") else
    c(Lasso = "lasso_nonzero", Enet = "enet_cox_nonzero", StepCox = "stepcox", CoxBoost = "coxboost_nonzero",
      RF = "rsf_minimal_depth", GBM = "gbm_cox_importance", glmBoost = "glmboost_cox_nonzero",
      GamBoost = "gamboost_cox_nonzero",
      XGBoost = "xgboost_cox_importance", PAMR = "pamr_survival_nonzero",
      SuperPC = "superpc_survival_scores", Ridge = "ridge_cox_weights",
      PLSR = "plsrcox_weights", SVM = "svm_survival_weights", Boruta = "boruta_survival_roughfix",
      ABESS = "abess_cox_gic_nonzero", LPC = "lpc_cox_nlogn", Stabsel = "stabsel_cox_lasso")
}

sb_doctor_selector <- function(methods, task) {
  menu <- sb_doctor_selector_menu(task)
  if (!is.character(methods) || !length(methods) || anyNA(methods)) stop("Supply ordered selector names.")
  methods[tolower(methods) == "elasticnet"] <- "Enet"
  if (task == "survival") methods[tolower(methods) == "mboost"] <- "glmBoost"
  if (task == "survival") methods[tolower(methods) == "rsf"] <- "RF"
  index <- match(tolower(methods), tolower(names(menu)))
  if (anyNA(index)) stop("Unknown selector. Choices: ", paste(names(menu), collapse = ", "))
  structure(list(task_type = task, methods = names(menu)[index], ids = unname(menu[index])),
    class = "stormbreaker_selector_chain_spec")
}

#' @rdname sb_pipeline
bi_selector <- function(methods) sb_doctor_selector(methods, "binary")

#' @rdname sb_pipeline
#' @export
surv_selector <- function(methods) sb_doctor_selector(methods, "survival")

sb_doctor_stage <- function(id, learner, base) {
  args <- list(selector = id, learner = learner)
  if (learner %in% c("glmnet_logistic", "glmnet_cox")) {
    args$inner_folds <- base$learner_params$inner_folds
    args$inner_seed <- base$learner_params$inner_seed
    args$glmnet_alpha <- base$learner_params$alpha
  }
  if (id == "stepcox") args$stepcox_direction <- "both"
  if (id %in% c("lasso_nonzero", "enet_cox_nonzero", "coxboost_nonzero")) {
    args$selector_inner_folds <- 5L
    args$selector_seed <- 1L
  }
  if (id == "rsf_minimal_depth") {
    args$selector_rsf_nodesize <- 5L
    args$selector_seed <- 1L
  }
  do.call(sb_pipeline, args)
}

sb_validate_chain_stages <- function(base, stages) {
  if (!is.list(stages) || !length(stages)) stop("A selector chain requires stages.")
  # Stage admission retains the restrictions on selection inside native CV.
  lapply(stages, function(s) {
    s <- sb_validate_pipeline_contract(s, base$task_type)
    if (!identical(s$learner, base$learner) || s$selector == "ordered_chain") stop("Invalid selector stage.")
    s
  })
}

sb_chain_pipeline_id <- function(base, stages) {
  paste(c("chain", vapply(stages, `[[`, character(1), "id"), base$id), collapse = "__then__")
}

sb_new_chain_pipeline <- function(base, stages) {
  if (!length(stages)) return(base)
  stages <- sb_validate_chain_stages(base, stages)
  result <- base
  result$selector <- "ordered_chain"
  result$selector_params <- list(stages = stages)
  result$chain_base <- base
  result$id <- sb_chain_pipeline_id(base, stages)
  result
}

sb_validate_chain_pipeline <- function(pipeline, task_type) {
  base <- sb_validate_pipeline_contract(pipeline$chain_base, task_type)
  stages <- pipeline$selector_params$stages
  stages <- sb_validate_chain_stages(base, stages)
  fields <- setdiff(sb_pipeline_execution_fields(task_type), c("id", "selector", "selector_params"))
  if (!identical(pipeline$id, sb_chain_pipeline_id(base, stages)) ||
    !all(vapply(fields, function(field) {
      if (field == "learner_params") {
        sb_pipeline_parameters_equal(pipeline[[field]], base[[field]])
      } else {
        identical(pipeline[[field]], base[[field]])
      }
    }, logical(1)))) {
    stop("Selector chain does not match its executable specification.")
  }
  pipeline$chain_base <- base
  pipeline$learner_params <- base$learner_params
  pipeline$selector_params$stages <- stages
  pipeline
}

sb_chain_selector_adapter <- function(params, capture_history = FALSE) {
  adapters <- lapply(params$stages, function(p) {
    if (capture_history) sb_captured_selector_adapter(p$selector, p$selector_params)
    else sb_resolve_selector_adapter(p$selector, p$selector_params)
  })
  list(id = "ordered_chain", params = params, needs_outcome = TRUE,
    fit = function(x, endpoint) {
      if (endpoint$task_type == "survival" && anyDuplicated(endpoint$group_ids) &&
          (any(vapply(params$stages, function(p) p$selector == "coxboost_nonzero", logical(1))) ||
           params$stages[[1L]]$learner %in% sb_analysis_tuned_survival_learner_ids())) {
        stop("This native survival CV engine cannot keep repeated patients together; choose a group-aware learner/selector.")
      }
      original <- names(x)
      source <- x
      states <- vector("list", length(adapters))
      for (i in seq_along(adapters)) {
        fit_stage <- function() {
          if (i == 1L || adapters[[i]]$id != "svm_rfe_binary") {
            return(sb_fit_runtime_selector(adapters[[i]], x, endpoint))
          }
          prefix <- sb_chain_selector_adapter(list(stages = params$stages[seq_len(i - 1L)]))
          prepare <- function(rows) {
            inner <- list(task_type = "binary", outcome = endpoint$outcome[rows],
              group_ids = endpoint$group_ids[rows])
            fitted <- prefix$fit(source[rows, , drop = FALSE], inner)
            list(train = prefix$apply(fitted, source[rows, , drop = FALSE]),
              test = prefix$apply(fitted, source[!rows, , drop = FALSE]))
          }
          controls <- adapters[[i]]$params
          stage <- adapters[[i]]
          # RFE's inner folds need the original input and the complete prefix,
          # not just the features selected on this outer analysis set.
          stage$upstream_context <- list(x = source,
            prefix = lapply(params$stages[seq_len(i - 1L)], function(p) {
              list(selector = p$selector, params = p$selector_params)
            }))
          stage$fit <- function(x, outcome, group_ids = NULL) {
            fitted <- sb_fit_svm_rfe(x, outcome, group_ids,
              controls$inner_folds, controls$seed, prepare_fold = prepare)
            fitted$params <- controls
            fitted
          }
          sb_fit_runtime_selector(stage, x, endpoint)
        }
        states[[i]] <- tryCatch(fit_stage(),
          error = function(e) sb_abort(paste0("Selector stage ", i, " (", adapters[[i]]$id, "): ", conditionMessage(e)),
            "stormbreaker_chain_stage_failed", parent = e))
        x <- adapters[[i]]$apply(states[[i]], x)
        if (!ncol(x)) sb_abort(paste("Selector stage", i, "retained no genes."), "stormbreaker_empty_selection")
      }
      structure(list(id = "ordered_chain", source_features = original,
        selected_features = names(x), params = params, stages = states,
        scores = data.frame(feature = original, selected = original %in% names(x))),
        class = "stormbreaker_selector_fit")
    },
    apply = function(fitted, x) {
      if (!identical(names(x), fitted$source_features)) stop("Predictors differ from the saved selector chain.")
      for (i in seq_along(adapters)) x <- adapters[[i]]$apply(fitted$stages[[i]], x)
      x
    })
}

# ---- selector-adapter ----

sb_fit_none_selector <- function(x) {
  x <- sb_validate_predictor_frame(x, "x")
  structure(
    list(
      id = "none",
      selected_features = names(x)
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_apply_none_selector <- function(fitted, x) {
  if (!inherits(fitted, "stormbreaker_selector_fit") ||
      !identical(fitted$id, "none") ||
      !is.character(fitted$selected_features) ||
      length(fitted$selected_features) < 1L) {
    sb_abort(
      "The fitted selector is not a valid `none` selector state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$selected_features)) {
    sb_abort(
      paste0(
        "Predictor columns must exactly match the fitted selector schema: ",
        paste(fitted$selected_features, collapse = ", "),
        "."
      ),
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x
}

sb_none_selector_adapter <- function() {
  list(
    id = "none",
    params = list(),
    needs_outcome = FALSE,
    fit = sb_fit_none_selector,
    apply = sb_apply_none_selector
  )
}

# ---- selector-fixed-features ----

sb_validate_fixed_features <- function(features) {
  valid <- is.character(features) && length(features) > 0L &&
    !anyNA(features) && all(nzchar(trimws(features))) &&
    anyDuplicated(features) == 0L
  if (!valid) {
    sb_abort(
      "`selector_features` must contain unique, non-empty feature names.",
      "stormbreaker_bad_fixed_features"
    )
  }
  features
}

sb_validate_feature_set_id <- function(feature_set_id) {
  valid <- is.character(feature_set_id) && length(feature_set_id) == 1L &&
    !is.na(feature_set_id) && nzchar(trimws(feature_set_id))
  if (!valid) {
    sb_abort(
      "`feature_set_id` must be one non-empty character value.",
      "stormbreaker_bad_feature_set_id"
    )
  }
  feature_set_id
}

sb_fit_fixed_features <- function(x, features, feature_set_id) {
  x <- sb_validate_predictor_frame(x, "x")
  features <- sb_validate_fixed_features(features)
  feature_set_id <- sb_validate_feature_set_id(feature_set_id)
  missing <- setdiff(features, names(x))
  if (length(missing) > 0L) {
    sb_abort(
      paste0(
        "The predeclared feature set is missing from the task predictors: ",
        paste(missing, collapse = ", "),
        "."
      ),
      "stormbreaker_missing_fixed_feature"
    )
  }
  structure(
    list(
      id = "fixed_features",
      feature_set_id = feature_set_id,
      source_features = names(x),
      selected_features = features
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_apply_fixed_features <- function(fitted, x) {
  valid_fit <- inherits(fitted, "stormbreaker_selector_fit") &&
    identical(fitted$id, "fixed_features") &&
    is.character(fitted$feature_set_id) &&
    length(fitted$feature_set_id) == 1L &&
    is.character(fitted$source_features) &&
    is.character(fitted$selected_features) &&
    length(fitted$selected_features) > 0L
  if (!valid_fit) {
    sb_abort(
      "The fitted selector is not a valid fixed-feature state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      paste0(
        "Predictor columns must exactly match the fitted selector schema: ",
        paste(fitted$source_features, collapse = ", "),
        "."
      ),
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_fixed_features_selector_adapter <- function(features, feature_set_id) {
  features <- sb_validate_fixed_features(features)
  feature_set_id <- sb_validate_feature_set_id(feature_set_id)
  list(
    id = "fixed_features",
    params = list(
      features = features,
      feature_set_id = feature_set_id,
      n_features = as.integer(length(features))
    ),
    needs_outcome = FALSE,
    fit = function(x) {
      sb_fit_fixed_features(x, features, feature_set_id)
    },
    apply = sb_apply_fixed_features
  )
}

sb_selector_is_learned <- function(selector) {
  selector %in% c(
    "ordered_chain",
    sb_survival_weight_selector_ids(),
    sb_survival_boosting_selector_ids(),
    sb_pamr_selector_ids(),
    sb_binary_model_selector_ids(),
    "univ_cox_top_k",
    "univ_binary_auc_top_k",
    "stepcox",
    "rsf_minimal_depth",
    "coxboost_nonzero",
    "lasso_nonzero",
    "enet_cox_nonzero",
    "rsf_then_stepcox",
    "coxboost_then_stepcox",
    "lasso_then_stepcox",
    sb_abess_selector_ids(),
    sb_forest_selector_ids(),
    sb_high_dimensional_selector_ids()
  )
}

# ---- selector-native ----

# ---- selector-pamr ----

sb_pamr_selector_ids <- function() c("pamr_binary_nonzero", "pamr_survival_nonzero")

sb_pamr_selector_params <- function(id) {
  if (id == "pamr_survival_nonzero") sb_pamr_survival_params() else
    sb_classical_binary_params("pamr_binary")
}

sb_pamr_selector_state <- function(engine, x, id, params) {
  selected_index <- pamr::pamr.predict(engine, t(as.matrix(x)),
    threshold = params$threshold, type = "nonzero")
  centroids <- pamr::pamr.predict(engine, t(as.matrix(x)),
    threshold = params$threshold, type = "centroid")
  score <- apply(abs((centroids - engine$centroid.overall) / engine$sd), 1L, max)
  if (length(score) != ncol(x) || any(!is.finite(score))) sb_abort(
    "PAMR returned invalid shrunken-centroid scores.", "stormbreaker_invalid_pamr_selector")
  if (!length(selected_index)) sb_abort(
    "PAMR selected no nonzero shrunken-centroid features; no genes were added.",
    "stormbreaker_empty_selection")
  selected <- seq_len(ncol(x)) %in% selected_index
  structure(list(id = id, params = params, source_features = names(x),
    selected_features = names(x)[selected],
    scores = data.frame(feature = names(x), score = unname(score), selected = selected),
    evidence = list(rule = "native_nonzero_shrunken_centroids",
      pamr_version = as.character(utils::packageVersion("pamr")))),
    class = "stormbreaker_selector_fit")
}

sb_pamr_selector_adapter <- function(id, params) {
  list(id = id, params = params, needs_outcome = TRUE,
    fit = if (id == "pamr_survival_nonzero") {
      function(x, time, event) {
        model <- sb_fit_pamr_survival(x, time, event, "pamr_survival", params)
        sb_pamr_selector_state(model$engine, x, id, params)
      }
    } else {
      function(x, outcome) {
        model <- sb_fit_classical_binary(x, outcome, "pamr_binary", params)
        sb_pamr_selector_state(model$engine, x, id, params)
      }
    },
    apply = function(fitted, x) {
      x <- sb_validate_predictor_frame(x, "x")
      if (!identical(names(x), fitted$source_features)) stop("Predictors differ from fitted PAMR selection.")
      x[fitted$selected_features]
    })
}

# ---- selector-svm-rfe ----

sb_svm_rfe_path <- function(x, outcome, assessment = NULL) {
  active <- names(x)
  subsets <- predictions <- vector("list", ncol(x))
  while (length(active)) {
    k <- length(active)
    train <- as.matrix(x[active])
    center <- colMeans(train)
    spread <- apply(train, 2L, stats::sd)
    spread[!is.finite(spread) | spread == 0] <- 1
    train <- sweep(sweep(train, 2L, center), 2L, spread, "/")
    model <- e1071::svm(train, factor(outcome, levels = c(0, 1)),
      kernel = "linear", cost = 1, scale = FALSE, probability = FALSE)
    subsets[[k]] <- active
    if (!is.null(assessment)) {
      test <- sweep(sweep(as.matrix(assessment[active]), 2L, center), 2L, spread, "/")
      predictions[[k]] <- as.integer(as.character(stats::predict(model, test)))
    }
    weights <- as.numeric(crossprod(model$coefs, model$SV))^2
    active <- active[-which.min(weights)]
  }
  list(subsets = subsets, predictions = predictions)
}

sb_fit_svm_rfe <- function(x, outcome, group_ids = NULL, folds = 5L, seed = 1L,
                           prepare_fold = NULL) {
  x <- sb_validate_predictor_frame(x, "x")
  if (length(outcome) != nrow(x) || anyNA(outcome) ||
      !all(outcome %in% c(0, 1)) || length(unique(outcome)) != 2L) {
    sb_abort("SVM-RFE requires both binary outcome classes.", "stormbreaker_bad_binary_outcome")
  }
  groups <- if (is.null(group_ids)) as.character(seq_len(nrow(x))) else
    sb_validate_group_ids(group_ids, nrow(x))
  fold_id <- sb_make_group_stratified_fold_id(outcome, groups, folds, seed,
    "stormbreaker_infeasible_inner_folds")
  correct <- vector("list", folds)
  fold_subsets <- vector("list", folds)
  for (fold in seq_len(folds)) {
    heldout <- fold_id == fold
    data <- if (is.null(prepare_fold)) list(train = x[!heldout, , drop = FALSE],
      test = x[heldout, , drop = FALSE]) else prepare_fold(!heldout)
    path <- sb_svm_rfe_path(data$train, outcome[!heldout], data$test)
    correct[[fold]] <- vapply(path$predictions,
      function(pred) sum(pred == outcome[heldout]), integer(1))
    fold_subsets[[fold]] <- path$subsets
  }
  candidates <- seq_len(min(ncol(x), lengths(correct)))
  correct <- do.call(rbind, lapply(correct, function(z) z[candidates]))
  accuracy <- colSums(correct) / nrow(x)
  best <- which.max(accuracy)
  final <- sb_svm_rfe_path(x, outcome)
  structure(list(id = "svm_rfe_binary", source_features = names(x),
    selected_features = final$subsets[[best]],
    evidence = list(curve = data.frame(n_features = candidates,
      accuracy = accuracy, error = 1 - accuracy, selected = candidates == best),
      fold_id = fold_id, group_ids = group_ids, fold_subsets = fold_subsets,
      subsets = final$subsets, correct = correct)), class = "stormbreaker_selector_fit")
}

# ---- selector-univ-binary-auc-top-k ----

sb_score_univariable_binary_auc_features <- function(x, outcome) {
  n_positive <- as.integer(sum(outcome == 1L))
  n_negative <- as.integer(sum(outcome == 0L))
  comparable_pairs <- as.numeric(n_positive) * as.numeric(n_negative)
  n_features <- ncol(x)
  auc <- rep(NA_real_, n_features)
  discrimination <- rep(NA_real_, n_features)
  direction <- rep(NA_character_, n_features)
  status <- rep("failed", n_features)
  error_class <- rep(NA_character_, n_features)
  error_message <- rep(NA_character_, n_features)

  for (index in seq_len(n_features)) {
    feature <- names(x)[index]
    values <- x[[index]]
    if (min(values) == max(values)) {
      error_class[index] <- "stormbreaker_selector_constant_feature"
      error_message[index] <- paste0(
        "Binary AUC scoring cannot rank constant predictor `",
        feature, "`."
      )
      next
    }

    ranks <- rank(values, ties.method = "average")
    pair_credit <- sum(ranks[outcome == 1L]) -
      n_positive * (n_positive + 1) / 2
    feature_auc <- unname(pair_credit / comparable_pairs)
    if (!is.finite(feature_auc) || feature_auc < 0 || feature_auc > 1) {
      error_class[index] <- "stormbreaker_invalid_selector_score"
      error_message[index] <- paste0(
        "Binary AUC scoring produced an invalid value for `",
        feature, "`."
      )
      next
    }

    auc[index] <- feature_auc
    discrimination[index] <- max(
      pair_credit,
      comparable_pairs - pair_credit
    ) / comparable_pairs
    direction[index] <- if (feature_auc > 0.5) {
      "higher_in_positive"
    } else if (feature_auc < 0.5) {
      "lower_in_positive"
    } else {
      "tied"
    }
    status[index] <- "success"
  }

  data.frame(
    feature = names(x),
    original_order = seq_len(n_features),
    auc = auc,
    discrimination = discrimination,
    direction = direction,
    n_positive = rep.int(n_positive, n_features),
    n_negative = rep.int(n_negative, n_features),
    comparable_pairs = rep.int(comparable_pairs, n_features),
    status = status,
    error_class = error_class,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

sb_fit_univ_binary_auc_top_k <- function(x, outcome, n_features) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "Binary AUC selection requires both outcome classes.",
      "stormbreaker_selector_single_class"
    )
  }
  n_features <- sb_check_positive_integer(
    n_features,
    "n_features",
    minimum = 1L
  )
  if (n_features > ncol(x)) {
    sb_abort(
      paste0(
        "The selector requested ", n_features,
        " features, but the task has only ", ncol(x),
        " candidate features."
      ),
      "stormbreaker_insufficient_candidate_features"
    )
  }

  scores <- sb_score_univariable_binary_auc_features(x, outcome)
  scores$rank <- NA_integer_
  scores$selected <- FALSE
  valid <- which(
    scores$status == "success" & is.finite(scores$discrimination)
  )
  ordered_valid <- valid[order(
    -scores$discrimination[valid],
    scores$original_order[valid]
  )]
  if (length(ordered_valid) < n_features) {
    sb_abort(
      paste0(
        "The selector requested ", n_features,
        " features, but only ", length(ordered_valid),
        " features had valid binary AUC scores."
      ),
      "stormbreaker_insufficient_valid_features"
    )
  }
  scores$rank[ordered_valid] <- seq_along(ordered_valid)
  selected_rows <- ordered_valid[seq_len(n_features)]
  scores$selected[selected_rows] <- TRUE
  selected_features <- scores$feature[selected_rows]

  structure(
    list(
      id = "univ_binary_auc_top_k",
      n_features = n_features,
      source_features = names(x),
      selected_features = selected_features,
      scores = scores
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_apply_univ_binary_auc_top_k <- function(fitted, x) {
  valid_fit <- inherits(fitted, "stormbreaker_selector_fit") &&
    identical(fitted$id, "univ_binary_auc_top_k") &&
    is.character(fitted$source_features) &&
    is.character(fitted$selected_features) &&
    length(fitted$selected_features) == fitted$n_features
  if (!valid_fit) {
    sb_abort(
      "The fitted selector is not a valid binary AUC Top-K state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      paste0(
        "Predictor columns must exactly match the fitted selector schema: ",
        paste(fitted$source_features, collapse = ", "),
        "."
      ),
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_univ_binary_auc_top_k_selector_adapter <- function(n_features) {
  n_features <- sb_check_positive_integer(
    n_features,
    "n_features",
    minimum = 1L
  )
  list(
    id = "univ_binary_auc_top_k",
    params = list(n_features = n_features),
    needs_outcome = TRUE,
    fit = function(x, outcome) {
      sb_fit_univ_binary_auc_top_k(x, outcome, n_features)
    },
    apply = sb_apply_univ_binary_auc_top_k
  )
}
