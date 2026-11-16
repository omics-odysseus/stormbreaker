# ---- results ----

# ---- composite-result-view ----

sb_is_composite_view <- function(x) {
  inherits(x, c("stormbreaker_repeated_cv_comparison",
                "stormbreaker_internal_external_cv", "stormbreaker_stack_cv_result"))
}

sb_view_bind <- function(rows) {
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    for (name in setdiff(columns, names(row))) row[[name]] <- rep(NA, nrow(row))
    row[columns]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

sb_view_context <- function() {
  list(repeat_id = NA_integer_, heldout_cohort_id = NA_character_,
       model_role = "pipeline", validation_scheme = "development_cv",
       view_id = NA_character_)
}

sb_view_leaves <- function(x, context = sb_view_context(), path = "result") {
  if (inherits(x, "stormbreaker_repeated_cv_comparison")) {
    if (!length(x$repetitions)) sb_abort("No saved repetitions.", "stormbreaker_bad_result_view")
    return(unlist(lapply(seq_along(x$repetitions), function(i) {
      ctx <- context
      ctx$repeat_id <- as.integer(i)
      ctx$validation_scheme <- "repeated_cv"
      sb_view_leaves(x$repetitions[[i]], ctx, paste0(path, "/repeat_", i))
    }), recursive = FALSE))
  }
  if (inherits(x, "stormbreaker_internal_external_cv")) {
    if (!length(x$folds)) sb_abort("No saved IECV folds.", "stormbreaker_bad_result_view")
    return(unlist(lapply(seq_along(x$folds), function(i) {
      fold <- x$folds[[i]]
      ctx <- context
      ctx$heldout_cohort_id <- fold$heldout_cohort_id
      ctx$validation_scheme <- "internal_external_cv"
      if (identical(x$strategy_type, "stack")) ctx$model_role <- "stack"
      sb_view_leaves(fold$comparison, ctx, paste0(path, "/cohort_", i))
    }), recursive = FALSE))
  }
  if (inherits(x, "stormbreaker_stack_cv_result")) {
    ctx <- context
    ctx$model_role <- "member"
    members <- sb_view_leaves(x$members, ctx, paste0(path, "/members"))
    # This local adapter only exposes the already saved stack outputs to readers.
    stack <- x[c("status", "task_type", "predictions", "metrics", "runs", "task")]
    stack$pipeline <- list(id = x$stack$id)
    # Stack descriptors are compact; retain unambiguous metadata saved per row.
    for (field in c("prediction_horizon", "time_unit")) {
      saved <- unique(stack$predictions[[field]])
      saved <- saved[!is.na(saved)]
      if (is.null(stack$task[[field]]) && length(saved) == 1L) {
        stack$task[[field]] <- saved[[1L]]
      }
    }
    for (component in c("predictions", "runs")) {
      if (!"pipeline_id" %in% names(stack[[component]])) {
        stack[[component]]$pipeline_id <- rep(x$stack$id, nrow(stack[[component]]))
      }
    }
    class(stack) <- if (identical(x$task_type, "binary")) {
      "stormbreaker_binary_cv_result"
    } else "stormbreaker_cv_result"
    ctx$model_role <- "stack"
    return(c(members, sb_view_leaves(stack, ctx, paste0(path, "/stack"))))
  }
  if (inherits(x, c("stormbreaker_cv_comparison", "stormbreaker_external_comparison"))) {
    values <- if (inherits(x, "stormbreaker_cv_comparison")) x$results else x$evaluations
    if (!length(values)) sb_abort("No saved comparison results.", "stormbreaker_bad_result_view")
    return(unlist(lapply(seq_along(values), function(i) {
      sb_view_leaves(values[[i]], context, paste0(path, "/pipeline_", i))
    }), recursive = FALSE))
  }
  if (!inherits(x, c("stormbreaker_binary_cv_result", "stormbreaker_cv_result",
                    "stormbreaker_external_evaluation"))) {
    sb_abort("Unsupported saved result in composite object.", "stormbreaker_bad_result_view")
  }
  if (inherits(x, "stormbreaker_external_evaluation") &&
      identical(context$validation_scheme, "development_cv")) {
    context$validation_scheme <- if (identical(x$assessment_role, "internal_holdout")) {
      "internal_holdout"
    } else "external_evaluation"
  }
  context$view_id <- path
  list(list(value = x, context = context))
}

sb_view_annotate <- function(tab, context) {
  for (name in names(context)) tab[[name]] <- rep(context[[name]], nrow(tab))
  tab
}

sb_view_boundary <- function(context) {
  switch(context$validation_scheme,
    internal_external_cv = "Development IECV (not final independent external validation)",
    internal_holdout = "Internal holdout",
    repeated_cv = "Development repeated CV", external_evaluation = "External evaluation",
    "Development CV")
}

sb_view_label <- function(context) {
  paste(c(if (!is.na(context$repeat_id)) paste("Repeat", context$repeat_id),
    if (!is.na(context$heldout_cohort_id)) paste("Held-out cohort", context$heldout_cohort_id),
    paste("Role:", context$model_role)), collapse = " | ")
}

sb_view_wrap_label <- function(labels) {
  vapply(strsplit(as.character(labels), "\n", fixed = TRUE), function(lines) {
    paste(unlist(lapply(lines, function(line) {
      if (!nzchar(line)) return("")
      starts <- seq.int(1L, nchar(line), by = 50L)
      substring(line, starts, pmin(starts + 49L, nchar(line)))
    }), use.names = FALSE), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

sb_composite_table <- function(x, component, metric, cohort_label) {
  context <- attr(x, "stormbreaker_view_context")
  attr(x, "stormbreaker_view_context") <- NULL
  leaves <- if (is.null(context)) sb_view_leaves(x) else list(list(value = x, context = context))
  sb_view_bind(lapply(leaves, function(leaf) {
    sb_view_annotate(sb_result_table(leaf$value, component, metric, cohort_label), leaf$context)
  }))
}

sb_composite_curves <- function(x, kind, cohort_label) {
  context <- attr(x, "stormbreaker_view_context")
  attr(x, "stormbreaker_view_context") <- NULL
  leaves <- if (is.null(context)) sb_view_leaves(x) else list(list(value = x, context = context))
  parts <- lapply(leaves, function(leaf) {
    view <- sb_curve_view_data(leaf$value, kind, cohort_label)
    for (name in c("data", "status")) {
      tab <- sb_view_annotate(view[[name]], leaf$context)
      tab$validation <- rep(leaf$context$validation_scheme, nrow(tab))
      target <- ifelse(tab$task_type == "binary", "class 1",
        ifelse(is.na(tab$prediction_horizon), "horizon not declared",
               paste(tab$prediction_horizon, tab$time_unit)))
      tab$panel <- paste(paste(tab$pipeline_id, target, sep = " | "),
        paste(c(cohort_label, sb_view_label(leaf$context)), collapse = " | "),
        sb_view_boundary(leaf$context),
        paste("curve:", tab$metric_status, "run:", tab$run_status),
        paste("Saved entry:", leaf$context$view_id), sep = "\n")
      view[[name]] <- tab
    }
    view
  })
  data <- sb_view_bind(lapply(parts, `[[`, "data"))
  status <- sb_view_bind(lapply(parts, `[[`, "status"))
  data$panel <- factor(data$panel, levels = unique(status$panel))
  data$line_group <- interaction(data$view_id, data$strategy, data$segment, drop = TRUE)
  list(data = data, status = status)
}

sb_composite_failures <- function(x, label) {
  context <- attr(x, "stormbreaker_view_context")
  attr(x, "stormbreaker_view_context") <- NULL
  leaves <- if (is.null(context)) sb_view_leaves(x) else list(list(value = x, context = context))
  sb_view_bind(lapply(leaves, function(leaf) {
    tab <- sb_view_annotate(sb_failure_rows(leaf$value, label), leaf$context)
    tab$validation <- rep(leaf$context$validation_scheme, nrow(tab))
    tab
  }))
}

# ---- failure-table ----

sb_failure_empty <- function() {
  data.frame(cohort_label=character(),validation=character(),task_type=character(),
    pipeline_id=character(),issue_type=character(),stage=character(),fold=integer(),
    metric=character(),status=character(),run_status=character(),
    source_error_stage=character(),error_class=character(),error_message=character(),
    estimate=numeric(),source_row=integer(),prediction_horizon=numeric(),time_unit=character(),
    stringsAsFactors=FALSE)
}

sb_failure_results <- function(x) {
  external <- inherits(x,c("stormbreaker_external_evaluation","stormbreaker_external_comparison"))
  is_cv <- function(value) inherits(value,c("stormbreaker_cv_result","stormbreaker_binary_cv_result"))
  results <- if (inherits(x,"stormbreaker_cv_comparison")) x$results else if (
    inherits(x,"stormbreaker_external_comparison")) x$evaluations else if (
      external || is_cv(x)) list(x) else NULL
  predicate <- if (external) function(value) inherits(value,"stormbreaker_external_evaluation") else is_cv
  if (!is.list(results) || !length(results) || !all(vapply(results,predicate,logical(1)))) {
    sb_abort("Supply supported CV or external results/comparisons.","stormbreaker_bad_failure_view")
  }
  list(results=results,external=external)
}

sb_failure_rows <- function(x,cohort_label) {
  if (sb_is_composite_view(x) || !is.null(attr(x, "stormbreaker_view_context"))) {
    return(sb_composite_failures(x, cohort_label))
  }
  parts <- sb_failure_results(x)
  rows <- list()
  value_or <- function(value,default) if(is.null(value) || length(value)!=1L) default else value
  problem <- function(status) !is.null(status) && length(status)==1L && !is.na(status) &&
    !status %in% c("success","not_requested","not_applicable")
  for (result in parts$results) {
    descriptor <- if(parts$external) result$endpoint else result$task
    type <- if(parts$external) result$task_type else if(inherits(result,"stormbreaker_binary_cv_result")) "binary" else "survival"
    pipeline <- if(parts$external) result$pipeline_id else result$pipeline$id
    add <- function(issue,stage,status,source,fold=NA_integer_,metric=NA_character_,
                    row=NA_integer_,estimate=NA_real_) {
      rows[[length(rows)+1L]] <<- data.frame(
        cohort_label=cohort_label,validation=if(parts$external) "external_evaluation" else "development_cv",
        task_type=type,pipeline_id=pipeline,issue_type=issue,stage=stage,fold=as.integer(fold),
        metric=metric,status=value_or(status,NA_character_),run_status=result$status,
        source_error_stage=value_or(source$error_stage,NA_character_),
        error_class=value_or(source$error_class,NA_character_),
        error_message=value_or(source$error_message,NA_character_),estimate=estimate,
        source_row=as.integer(row),
        prediction_horizon=value_or(source$prediction_horizon,value_or(descriptor$prediction_horizon,NA_real_)),
        time_unit=value_or(source$time_unit,value_or(descriptor$time_unit,NA_character_)),
        stringsAsFactors=FALSE)
    }
    if (problem(result$status)) {
      add("pipeline",if(parts$external) "external_evaluation" else "cv",result$status,result)
    }
    if (!parts$external && is.data.frame(result$runs)) {
      for(i in which(vapply(as.list(result$runs$status),problem,logical(1)))) {
        source <- as.list(result$runs[i,,drop=FALSE])
        add("run","cv_fold",source$status,source,fold=value_or(source$fold,NA_integer_))
      }
    }
    visit <- function(value,path) {
      if (!is.list(value) || is.data.frame(value)) return(invisible(NULL))
      if(problem(value$status)) {
        add("metric","metric",value$status,value,metric=path,
            estimate=value_or(value$estimate,NA_real_))
      } else if(identical(value$status,"success") && is.numeric(value$estimate) &&
                length(value$estimate)==1L && !is.finite(value$estimate)) {
        add("metric_nonfinite","metric",value$status,value,metric=path,estimate=value$estimate)
      }
      for(name in names(value)) {
        child <- value[[name]]
        if(is.data.frame(child) && "status" %in% names(child)) {
          # Only saved aggregate curve/fold diagnostic rows, never predictions.
          if(!name %in% c("curve","per_fold")) next
          for(i in which(vapply(as.list(child$status),problem,logical(1)))) {
            source <- as.list(child[i,,drop=FALSE])
            add(if(name=="curve") "curve_row" else "metric_fold",name,
                source$status,source,fold=value_or(source$fold,NA_integer_),metric=path,row=i)
          }
        } else if(is.list(child) && !is.data.frame(child)) {
          visit(child,paste(path,name,sep="."))
        }
      }
      invisible(NULL)
    }
    for(metric in names(result$metrics)) visit(result$metrics[[metric]],metric)
  }
  if(!length(rows)) return(sb_failure_empty())
  result <- do.call(rbind,rows); rownames(result) <- NULL
  result
}

#' @rdname sb_render_report
#' @export
sb_failure_table <- function(x,cohort_label=NULL) {
  label <- sb_result_cohort_label(cohort_label)
  supported <- inherits(x,c("stormbreaker_cv_result","stormbreaker_binary_cv_result",
    "stormbreaker_cv_comparison","stormbreaker_external_evaluation","stormbreaker_external_comparison"))
  if(supported || sb_is_composite_view(x)) return(sb_failure_rows(x,label))
  if(!is.null(cohort_label) || !is.list(x) || !length(x) || is.null(names(x)) ||
     anyNA(names(x)) || any(!nzchar(trimws(names(x)))) || anyDuplicated(names(x))) {
    sb_abort("Supply a supported result or a uniquely named non-empty result list.",
             "stormbreaker_bad_failure_view")
  }
  rows <- lapply(seq_along(x),function(i) sb_failure_rows(x[[i]],names(x)[i]))
  output <- sb_view_bind(rows); rownames(output) <- NULL
  output
}

sb_new_output_path <- function(file,extension) {
  if(!is.character(file) || length(file)!=1L || is.na(file) || !nzchar(trimws(file)) ||
     !endsWith(tolower(file),paste0(".",extension)) || !dir.exists(dirname(file))) {
    sb_abort(paste0("Supply a new .",extension," file in an existing directory."),
             "stormbreaker_bad_output_path")
  }
  if(file.exists(file)) sb_abort("Output already exists; choose a new filename.","stormbreaker_output_exists")
  file.path(normalizePath(dirname(file),winslash="/",mustWork=TRUE),basename(file))
}

#' @rdname sb_render_report
#' @export
sb_export_failures <- function(x,file,cohort_label=NULL) {
  file <- sb_new_output_path(file,"csv")
  tab <- sb_failure_table(x,cohort_label)
  file <- sb_new_output_path(file,"csv")
  utils::write.csv(tab,file,row.names=FALSE,na="NA",fileEncoding="UTF-8")
  invisible(normalizePath(file,winslash="/",mustWork=TRUE))
}

# ---- model-display ----

sb_display_learner <- function(learner, parameters = list(), model_id = NULL) {
  learner <- sb_canonical_learner_id(learner)

  if (is.null(learner) || !length(learner) || is.na(learner)) return("Model")
  if (learner %in% c("glmnet_logistic", "glmnet_cox")) {
    alpha <- parameters$alpha
    if (is.null(alpha)) alpha <- if (!is.null(model_id) && grepl("^ridge_", model_id)) 0 else
      if (!is.null(model_id) && grepl("^elastic_net_", model_id)) 0.5 else 1
    penalty <- if (identical(as.numeric(alpha), 0)) "Ridge" else
      if (identical(as.numeric(alpha), 1)) "Lasso" else "Elastic net"
    return(paste(penalty, if (learner == "glmnet_cox") "Cox" else "logistic"))
  }
  known <- c(
    logistic = "Logistic regression", coxph = "Cox regression",
    gaussian_naive_bayes = "Gaussian naive Bayes",
    ranger_binary = "Random forest (ranger)",
    ranger_survival = "Random survival forest (ranger)",
    random_forest_binary = "Random forest (randomForest)",
    rsf_survival = "Random survival forest (randomForestSRC)",
    lda_binary = "LDA", qda_binary = "QDA", fda_binary = "FDA",
    pda_binary = "Penalized FDA", mda_binary = "MDA",
    slda_binary = "Stabilized LDA", sda_binary = "Shrinkage DA",
    earth_binary = "MARS", decision_tree_binary = "Decision tree (rpart)",
    tree_binary = "Decision tree (tree)", rpart_survival = "Survival tree (rpart)",
    knn_binary = "k-nearest neighbours", kernel_weighted_knn_binary = "Weighted kNN",
    kernel_naive_bayes_binary = "Kernel naive Bayes",
    nearest_centroid_binary = "Nearest centroid",
    bernoulli_naive_bayes_binary = "Bernoulli naive Bayes",
    gaussian_process_binary = "Gaussian process",
    nnet_binary = "Neural network (nnet)",
    neuralnet_multilayer_binary = "Multilayer neural network",
    aorsf_binary = "Oblique random forest (aorsf)",
    aorsf_survival = "Oblique survival forest (aorsf)",
    brglm2_logistic = "Bias-reduced logistic",
    bayesglm_logistic_cauchy = "Bayesian logistic (Cauchy)",
    robust_logistic_mqle = "Robust logistic (Mqle)",
    coxphf_firth = "Firth Cox regression",
    hda_binary = "High-dimensional discriminant analysis",
    dwd_binary = "Distance-weighted discrimination",
    ktsp_binary = "k-top scoring pairs",
    rotation_forest_binary = "Rotation forest",
    rrf_binary = "Regularized random forest",
    random_ferns_binary = "Random ferns",
    c50_rules_binary = "C5.0 rules", c50_binary = "C5.0",
    bagging_binary = "Bagged trees",
    bagged_knn_survival = "Bagged survival kNN",
    akritas_survival = "Akritas",
    survreg_weibull = "Weibull AFT", survreg_exponential = "Exponential AFT",
    survreg_lognormal = "Log-normal AFT", survreg_loglogistic = "Log-logistic AFT",
    mtlr_survival = "MTLR", ctree_binary = "Conditional inference tree",
    ctree_survival = "Conditional inference survival tree",
    cforest_binary = "Conditional inference forest",
    cforest_survival = "Conditional inference survival forest"
  )
  if (learner %in% names(known)) return(unname(known[[learner]]))
  tokens <- strsplit(learner, "_", fixed = TRUE)[[1L]]
  tokens <- tokens[!tokens %in% c("fixed", "binary", "survival")]
  dictionary <- c(
    xgboost = "XGBoost", lightgbm = "LightGBM", catboost = "CatBoost",
    gbm = "GBM", cox = "Cox", aft = "AFT", dart = "DART", goss = "GOSS",
    rf = "RF", linear = "linear", svm = "SVM", lssvm = "LS-SVM",
    radial = "RBF", rbf = "RBF", polynomial = "polynomial", sigmoid = "sigmoid",
    survivalsvm = "Survival SVM", vanbelle1 = "Van Belle 1", vanbelle2 = "Van Belle 2",
    ranger = "ranger", extratrees = "Extra Trees", hellinger = "Hellinger",
    cindex = "C-index", maxstat = "maxstat", rlt = "RLT",
    adaboost = "AdaBoost", logitboost = "LogitBoost", bst = "BST",
    rda = "RDA", evtree = "Evolutionary tree", coxboost = "CoxBoost",
    plsrcox = "PLS-Cox", plsrglm = "PLS-GLM", superpc = "SuperPC",
    gamboost = "GamBoost", glmboost = "GLMBoost", blackboost = "BlackBoost",
    gehan = "Gehan", unoc = "Uno C", ncvreg = "ncvreg", mcp = "MCP",
    scad = "SCAD", bic = "BIC", gic = "GIC", abess = "ABESS",
    gam = "GAM", reml = "REML", bart = "BART", pamr = "PAM",
    flexsurv = "flexsurv", flexsurvspline = "flexsurvspline",
    liblinear = "LIBLINEAR", l1 = "L1", l2 = "L2", svc = "SVC",
    l2loss = "L2 loss", l1loss = "L1 loss", lvq1 = "LVQ1",
    olvq1 = "OLVQ1", lvq21 = "LVQ2.1", lvq3 = "LVQ3",
    stepglm = "Stepwise logistic", pre = "PRE", ahaz = "ahaz",
    aalen = "Aalen", bujar = "Bujar", penaft = "penAFT",
    cv = "CV-tuned", lasso = "Lasso", sscad = "S-SCAD",
    c50 = "C5.0"
  )
  matched <- match(tokens, names(dictionary))
  tokens[!is.na(matched)] <- unname(dictionary[matched[!is.na(matched)]])
  paste(tokens[nzchar(tokens)], collapse = " ")
}

sb_display_selector <- function(selector, parameters = list()) {
  if (identical(selector, "ordered_chain")) return(paste(vapply(parameters$stages,
    function(p) sb_display_selector(p$selector, p$selector_params), character(1)), collapse = " -> "))
  if (is.null(selector) || !length(selector) || is.na(selector)) return("Not recorded")
  known <- c(
    lasso_binary_nonzero = "LASSO", stepglm_binary_aic = "Stepglm",
    enet_binary_nonzero = "Enet", gbm_binary_importance = "GBM",
    glmboost_binary_nonzero = "glmBoost", rf_binary_importance = "RF importance",
    none = "None", univ_cox_top_k = "Univariate Cox", univ_binary_auc_top_k = "Univariate AUC",
    stepcox = "Stepwise Cox", abess_binary_gic_nonzero = "ABESS",
    abess_cox_gic_nonzero = "ABESS", boruta_binary_roughfix = "Boruta",
    boruta_survival_roughfix = "Boruta", varselrf_binary_oob_1se = "varSelRF",
    sis_binary_nlogn = "SIS", vsurf_binary_prediction = "VSURF",
    mrmr_binary_nlogn = "mRMR", jmi_binary_nlogn = "JMI", cmim_binary_nlogn = "CMIM",
    disr_binary_nlogn = "DISR", jmim_binary_nlogn = "JMIM",
    lpc_binary_nlogn = "LPC", lpc_cox_nlogn = "LPC",
    stabsel_binary_lasso = "Stability selection", stabsel_cox_lasso = "Stability selection",
    rsf_minimal_depth = "RSF minimal depth", coxboost_nonzero = "CoxBoost",
    lasso_nonzero = "Lasso", enet_cox_nonzero = "Enet", rsf_then_stepcox = "RSF + Stepwise Cox",
    coxboost_then_stepcox = "CoxBoost + Stepwise Cox",
    lasso_then_stepcox = "Lasso + Stepwise Cox", fixed_features = "Specified gene set"
  )
  label <- if (selector %in% names(known)) unname(known[[selector]]) else gsub("_", " ", selector)
  k <- if (!is.null(parameters$n_features)) parameters$n_features else parameters$k
  if (!is.null(k)) label <- paste0(label, " (", k, ")")
  if (!is.null(parameters$direction)) label <- paste0(label, " (", parameters$direction, ")")
  label
}

sb_display_parameters <- function(parameters) {
  if (!length(parameters)) return("")
  values <- vapply(parameters, function(value) {
    if (is.atomic(value) && length(value) <= 8L) paste(value, collapse = ", ") else
      paste0("[", length(value), " values]")
  }, character(1))
  paste(paste(names(values), values, sep = " = "), collapse = "; ")
}

# ---- result-table ----

#' @rdname sb_render_report
#' @export
sb_result_table <- function(x,
                            component = c("metrics", "predictions", "runs"),
                            metric = NULL, cohort_label = NULL) {
  component <- match.arg(component)
  label <- sb_result_cohort_label(cohort_label)
  if (sb_is_composite_view(x) || !is.null(attr(x, "stormbreaker_view_context"))) {
    return(sb_composite_table(x, component, metric, cohort_label))
  }
  if (inherits(x, c("stormbreaker_external_evaluation",
                    "stormbreaker_external_comparison"))) {
    return(sb_external_result_table(x, component, metric, cohort_label))
  }
  is_cv <- function(value) {
    inherits(value, c("stormbreaker_binary_cv_result", "stormbreaker_cv_result"))
  }
  if (inherits(x, "stormbreaker_cv_comparison")) {
    results <- x$results
  } else if (is_cv(x)) {
    results <- list(x)
  } else {
    sb_abort(
      "`x` must be a CV result/comparison or external evaluation/comparison.",
      "stormbreaker_bad_result_view"
    )
  }
  if (!is.list(results) || length(results) == 0L ||
      !all(vapply(results, is_cv, logical(1)))) {
    sb_abort(
      "`x` must contain at least one supported CV result.",
      "stormbreaker_bad_result_view"
    )
  }
  if (identical(component, "metrics")) {
    supported <- sb_scalar_view_metrics()
    if (!is.character(metric) || length(metric) != 1L ||
        is.na(metric) || !metric %in% supported) {
      sb_abort(
        paste0(
          "Unsupported `metric`; explicitly choose one of: ",
          paste(supported, collapse = ", "), "."
        ),
        "stormbreaker_bad_result_metric"
      )
    }
  } else if (!is.null(metric)) {
    sb_abort(
      "`metric` is only used for `component = \"metrics\"`; otherwise leave it NULL.",
      "stormbreaker_unused_result_metric"
    )
  }

  rows <- lapply(results, function(result) {
    task_type <- if (inherits(result, "stormbreaker_binary_cv_result")) {
      "binary"
    } else {
      "survival"
    }
    if (!identical(component, "metrics")) {
      tab <- result[[component]]
      tab$task_type <- rep(task_type, nrow(tab))
      return(tab)
    }
    if (metric %in% sb_horizon_view_metrics()) {
      return(sb_horizon_metric_row(result, metric, task_type))
    }
    value <- result$metrics[[metric]]
    if (is.null(value)) {
      sb_abort(
        paste0(
          "Metric `", metric, "` is not saved in pipeline `",
          result$pipeline$id, "`. No metric is recomputed by this function."
        ),
        "stormbreaker_result_metric_not_saved"
      )
    }
    data.frame(
      pipeline_id = result$pipeline$id,
      task_type = task_type,
      metric = metric,
      estimate = value$estimate,
      status = value$status,
      run_status = result$status,
      scope = value$scope,
      aggregation = value$aggregation,
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  if (!is.null(cohort_label)) output$cohort_label <- rep(label, nrow(output))
  output
}

# ---- select-model ----

#' @rdname sb_complete_run
#' @export
sb_select_model <- function(comparison) {
  if (!inherits(comparison, "stormbreaker_cohort_comparison")) {
    sb_abort("Supply the result of sb_compare_cohorts().",
      "stormbreaker_bad_cohort_comparison")
  }
  cohorts <- names(comparison$cohort_roles)[comparison$cohort_roles == "validation"]
  if (!length(cohorts)) {
    sb_abort("Model selection requires at least one cohort with role 'validation'.",
      "stormbreaker_no_validation_cohort")
  }
  metric <- if (identical(comparison$task_type, "binary")) "roc_auc" else "harrell_c"
  ids <- vapply(comparison$model_pipelines, function(x) x$id, character(1))
  fit_ids <- vapply(comparison$fits, function(x) x$pipeline$id, character(1))
  fit_index <- match(ids, fit_ids)
  cv <- comparison$development$comparison
  eligible <- !is.na(fit_index) & vapply(fit_index, function(i) {
    !is.na(i) && identical(comparison$fits[[i]]$status, "success")
  }, logical(1)) & cv$status[match(ids, cv$pipeline_id)] == "success"
  exclusion_reason <- rep("", length(ids))
  diagnostics <- list()
  check_variation <- function(results, field, cohort, scope) {
    states <- vapply(ids, function(id) {
      sb_prediction_variation(results[[id]][[field]], comparison$task_type)$status
    }, character(1))
    diagnostics[[length(diagnostics) + 1L]] <<- data.frame(
      pipeline_id = ids, cohort = cohort, prediction_variation = unname(states))
    invalid <- !states %in% c("variable", "partly_constant")
    reason <- ifelse(states == "constant", paste0("constant_", scope, "_predictions"),
      paste0("unavailable_", scope, "_prediction_evidence"))
    if (scope == "validation") reason <- paste0(reason, ":", cohort)
    add <- which(invalid)
    exclusion_reason[add] <<- ifelse(nzchar(exclusion_reason[add]),
      paste(exclusion_reason[add], reason[add], sep = ";"), reason[add])
    eligible <<- eligible & !invalid
  }
  check_variation(comparison$development$results, "predictions", "development_cv", "development")
  scores <- matrix(NA_real_, nrow = length(ids), ncol = length(cohorts),
    dimnames = list(ids, cohorts))
  for (cohort in cohorts) {
    tab <- comparison$cohorts[[cohort]]$comparison
    rows <- match(ids, tab$pipeline_id)
    scores[, cohort] <- tab[[metric]][rows]
    eligible <- eligible & tab$status[rows] == "success" &
      tab[[paste0(metric, "_status")]][rows] == "success" &
      is.finite(scores[, cohort])
    check_variation(comparison$cohorts[[cohort]]$evaluations, "scored_predictions", cohort, "validation")
  }
  eligible[is.na(eligible)] <- FALSE
  exclusion_reason[!eligible & !nzchar(exclusion_reason)] <- "incomplete_fit_or_metrics"
  means <- rowMeans(scores)
  means[!eligible] <- NA_real_
  ranking <- data.frame(pipeline_id = ids, eligible = eligible,
    validation_mean = unname(means), declared_order = seq_along(ids),
    exclusion_reason = exclusion_reason,
    stringsAsFactors = FALSE)
  if (!any(eligible)) {
    sb_abort(paste0("No eligible model has successful development CV, a final fit, ",
      "nonconstant prediction evidence and finite ", metric, " in every validation cohort. ",
      "Exclusions: ", paste(unique(exclusion_reason), collapse = "; "), "."),
      "stormbreaker_no_eligible_model")
  }
  best <- which.max(replace(means, !eligible, -Inf))
  ranking$selected <- seq_along(ids) == best
  model <- comparison$fits[[fit_index[best]]]
  model$model_selection <- list(metric = metric, cohorts = cohorts,
    rule = "maximum_equal_weight_validation_mean; ties_in_declaration_order",
    scores = scores,
    prediction_diagnostics = do.call(rbind, diagnostics),
    ranking = ranking[order(!eligible, -means, ranking$declared_order), , drop = FALSE])
  model
}
