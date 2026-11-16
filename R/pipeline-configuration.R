# ---- pipeline-configuration ----

# ---- pipeline-identity ----

sb_pipeline_model_id <- function(base_id, learner, learner_params, uses_inner_cv) {
  is_catboost <- learner %in% c(
    "catboost_binary",
    "catboost_cox",
    "catboost_aft_survival"
  )
  is_lightgbm <- learner %in% sb_lightgbm_binary_learner_ids()
  is_gbm_fixed <- identical(learner, "gbm_cox")
  is_analysis_tuned <- learner %in% sb_analysis_tuned_survival_learner_ids()
  is_coxboost <- identical(learner, "coxboost")
  is_plsrcox <- identical(learner, "plsrcox")
  is_superpc <- identical(learner, "superpc")
  is_survivalsvm <- learner %in% sb_survivalsvm_learner_ids()
  is_classical_binary_boosting <- learner %in%
    sb_classical_binary_boosting_learner_ids()
  is_ranger <- learner %in% sb_ranger_learner_ids()
  pipeline_id <- if (uses_inner_cv) {
    standardize_suffix <- if (isTRUE(learner_params$standardize)) {
      ""
    } else {
      "__standardize0"
    }
    lambda_rule_suffix <- if (identical(
      learner_params$selection_rule,
      "lambda.min"
    )) {
      ""
    } else {
      "__lambda1se"
    }
    measure_suffix <- if (identical(learner_params$measure, "deviance")) {
      ""
    } else {
      paste0("__measure", learner_params$measure)
    }
    paste0(
      base_id,
      "__alpha", format(
        learner_params$alpha,
        trim = TRUE,
        scientific = FALSE,
        digits = 15L
      ),
      standardize_suffix,
      lambda_rule_suffix,
      measure_suffix,
      "__inner", learner_params$inner_folds,
      "__seed", learner_params$inner_seed
    )
  } else if (identical(learner, "gaussian_naive_bayes")) {
    paste0(
      base_id,
      "__smooth",
      sb_format_roundtrip_double(learner_params$var_smoothing)
    )
  } else if (is_lightgbm) {
    paste0(
      base_id,
      "__rounds", learner_params$nrounds,
      "__lr", sb_format_roundtrip_double(
        learner_params$learning_rate
      ),
      "__leaves", learner_params$num_leaves,
      "__minleaf", learner_params$min_data_in_leaf,
      "__seed", learner_params$seed
    )
  } else if (is_catboost) {
    paste0(
      base_id,
      "__iter", learner_params$iterations,
      "__lr", sb_format_roundtrip_double(
        learner_params$learning_rate
      ),
      "__depth", learner_params$depth,
      "__l2reg", sb_format_roundtrip_double(
        learner_params$l2_leaf_reg
      ),
      "__seed", learner_params$random_seed
    )
  } else if (identical(learner, "gbm_cox_cv")) {
    paste0(
      base_id,
      "__inner", learner_params$inner_folds,
      "__seed", learner_params$seed,
      "__maxtrees", learner_params$maximum_trees
    )
  } else if (is_analysis_tuned) {
    paste0(
      base_id,
      "__inner", learner_params$inner_folds,
      "__seed", learner_params$seed
    )
  } else if (is_gbm_fixed) {
    paste0(
      base_id,
      "__trees", learner_params$n_trees,
      "__depth", learner_params$interaction_depth,
      "__minnode", learner_params$min_obs_node,
      "__shrink", sb_format_roundtrip_double(
        learner_params$shrinkage
      ),
      "__bag", sb_format_roundtrip_double(
        learner_params$bag_fraction
      ),
      "__seed", learner_params$seed
    )
  } else if (is_ranger) {
    paste0(
      base_id,
      "__trees", learner_params$num_trees,
      "__mtry", if (is.null(learner_params$mtry)) {
        "sqrt"
      } else {
        learner_params$mtry
      },
      "__minnode", learner_params$min_node_size,
      "__seed", learner_params$seed
    )
  } else if (is_coxboost) {
    paste0(
      base_id,
      "__steps", learner_params$stepno,
      "__penalty", if (is.null(learner_params$penalty)) {
        "event9"
      } else {
        sb_format_roundtrip_double(learner_params$penalty)
      }
    )
  } else if (is_plsrcox) {
    paste0(base_id, "__components", learner_params$components)
  } else if (is_superpc) {
    paste0(
      base_id,
      "__threshold", sb_format_roundtrip_double(learner_params$threshold),
      "__components", learner_params$components
    )
  } else if (is_survivalsvm) {
    paste0(
      base_id,
      "__gamma", sb_format_roundtrip_double(learner_params$gamma_mu)
    )
  } else if (is_classical_binary_boosting) {
    paste0(base_id, "__iterations", learner_params$iterations)
  } else if (identical(learner, "stepglm_binary")) {
    paste0(base_id, "__direction", learner_params$direction)
  } else {
    base_id
  }
  pipeline_id
}

# ---- pipeline-learner-config ----

sb_pipeline_tuning_params <- function(learner, pipeline_task_type, uses_inner_cv, options) {
  is_analysis_tuned <- learner %in% sb_analysis_tuned_survival_learner_ids()
  learner_params <- list()
  if (!uses_inner_cv) {
    if (!is.null(options$inner_folds) || !is.null(options$inner_seed)) {
      sb_abort(
        paste0(
          "`inner_folds` and `inner_seed` apply only when ",
          "a glmnet learner is requested."
        ),
        "stormbreaker_irrelevant_inner_parameter"
      )
    }
    if (!is.null(options$glmnet_alpha) || !is.null(options$glmnet_standardize) ||
      !is.null(options$glmnet_lambda_rule) || !is.null(options$glmnet_tuning_measure)) {
      sb_abort(
        paste0(
          "`glmnet_alpha`, `glmnet_standardize`, `glmnet_lambda_rule`, ",
          "and `glmnet_tuning_measure` ",
          "apply only when a glmnet learner is requested."
        ),
        "stormbreaker_irrelevant_glmnet_parameter"
      )
    }
  } else {
    if (is.null(options$inner_folds) || is.null(options$inner_seed)) {
      sb_abort(
        paste0(
          "`inner_folds` and `inner_seed` are required when ",
          "a glmnet learner is requested."
        ),
        "stormbreaker_missing_inner_parameter"
      )
    }
    if (is.null(options$glmnet_alpha)) {
      options$glmnet_alpha <- 1
    }
    if (!is.numeric(options$glmnet_alpha) || length(options$glmnet_alpha) != 1L ||
      is.na(options$glmnet_alpha) || !is.finite(options$glmnet_alpha) ||
      options$glmnet_alpha < 0 || options$glmnet_alpha > 1) {
      sb_abort(
        "`glmnet_alpha` must be one finite number between 0 and 1.",
        "stormbreaker_bad_glmnet_alpha"
      )
    }
    if (is.null(options$glmnet_standardize)) {
      options$glmnet_standardize <- TRUE
    }
    if (!is.logical(options$glmnet_standardize) ||
      length(options$glmnet_standardize) != 1L || is.na(options$glmnet_standardize)) {
      sb_abort(
        "`glmnet_standardize` must be `TRUE` or `FALSE`.",
        "stormbreaker_bad_glmnet_standardize"
      )
    }
    if (is.null(options$glmnet_lambda_rule)) {
      options$glmnet_lambda_rule <- "lambda.min"
    }
    options$glmnet_lambda_rule <- sb_validate_glmnet_lambda_rule(
      options$glmnet_lambda_rule,
      argument = "glmnet_lambda_rule"
    )
    if (is.null(options$glmnet_tuning_measure)) {
      options$glmnet_tuning_measure <- "deviance"
    }
    options$glmnet_tuning_measure <- sb_validate_glmnet_tuning_measure(
      options$glmnet_tuning_measure,
      pipeline_task_type,
      argument = "glmnet_tuning_measure"
    )
    learner_params <- list(
      alpha = as.numeric(options$glmnet_alpha),
      standardize = options$glmnet_standardize,
      inner_folds = sb_check_positive_integer(
        options$inner_folds,
        "inner_folds",
        minimum = 3L
      ),
      inner_seed = sb_check_positive_integer(
        options$inner_seed,
        "inner_seed",
        minimum = 1L
      ),
      selection_rule = options$glmnet_lambda_rule,
      measure = options$glmnet_tuning_measure
    )
  }
  has_analysis_tuning_control <- !is.null(options$learner_tuning_folds) ||
    !is.null(options$learner_seed) || !is.null(options$gbm_cv_max_trees)
  if (!is_analysis_tuned && has_analysis_tuning_control) {
    sb_abort(
      paste0(
        "`learner_tuning_folds`, `learner_seed`, and `gbm_cv_max_trees` ",
        "apply only to an analysis-fold tuned survival learner."
      ),
      "stormbreaker_irrelevant_learner_tuning_parameter"
    )
  }
  if (is_analysis_tuned) {
    if (is.null(options$learner_tuning_folds)) {
      options$learner_tuning_folds <- if (identical(learner, "plsrcox_cv")) {
        5L
      } else {
        10L
      }
    }
    if (is.null(options$learner_seed)) {
      options$learner_seed <- 1L
    }
    if (!identical(learner, "gbm_cox_cv") &&
      !is.null(options$gbm_cv_max_trees)) {
      sb_abort(
        "`gbm_cv_max_trees` applies only to `gbm_cox_cv`.",
        "stormbreaker_irrelevant_learner_tuning_parameter"
      )
    }
    learner_params <- switch(learner,
      coxboost_cv = sb_coxboost_cv_params(
        inner_folds = options$learner_tuning_folds,
        seed = options$learner_seed
      ),
      plsrcox_cv = sb_plsrcox_cv_params(
        inner_folds = options$learner_tuning_folds,
        seed = options$learner_seed
      ),
      superpc_cv = sb_superpc_cv_params(
        inner_folds = options$learner_tuning_folds,
        seed = options$learner_seed
      ),
      gbm_cox_cv = sb_gbm_cox_cv_params(
        inner_folds = options$learner_tuning_folds,
        seed = options$learner_seed,
        maximum_trees = if (is.null(options$gbm_cv_max_trees)) {
          10000L
        } else {
          options$gbm_cv_max_trees
        }
      )
    )
  }
  learner_params
}

sb_pipeline_tree_params <- function(learner, learner_params, options) {
  is_catboost <- learner %in% c(
    "catboost_binary",
    "catboost_cox",
    "catboost_aft_survival"
  )
  is_lightgbm <- learner %in% sb_lightgbm_binary_learner_ids()
  is_gbm_fixed <- identical(learner, "gbm_cox")
  is_ranger_survival <- learner %in% sb_ranger_survival_learner_ids()
  is_ranger <- learner %in% sb_ranger_learner_ids()
  if (identical(learner, "gaussian_naive_bayes")) {
    if (is.null(options$naive_bayes_var_smoothing)) {
      options$naive_bayes_var_smoothing <- 1e-9
    }
    valid_smoothing <- is.numeric(options$naive_bayes_var_smoothing) &&
      length(options$naive_bayes_var_smoothing) == 1L &&
      !is.na(options$naive_bayes_var_smoothing) &&
      is.finite(options$naive_bayes_var_smoothing) &&
      options$naive_bayes_var_smoothing >= 0
    if (!valid_smoothing) {
      sb_abort(
        paste0(
          "`naive_bayes_var_smoothing` must be one finite, ",
          "non-negative number."
        ),
        "stormbreaker_bad_naive_bayes_var_smoothing"
      )
    }
    learner_params <- list(
      var_smoothing = as.numeric(options$naive_bayes_var_smoothing),
      variance_estimator = "population",
      prior = "empirical"
    )
  } else if (!is.null(options$naive_bayes_var_smoothing)) {
    sb_abort(
      paste0(
        "`naive_bayes_var_smoothing` applies only when ",
        "Gaussian Naive Bayes is requested."
      ),
      "stormbreaker_irrelevant_naive_bayes_parameter"
    )
  }
  lightgbm_arguments <- list(
    nrounds = options$lightgbm_nrounds,
    learning_rate = options$lightgbm_learning_rate,
    num_leaves = options$lightgbm_num_leaves,
    min_data_in_leaf = options$lightgbm_min_data_in_leaf,
    seed = options$lightgbm_seed
  )
  if (is_lightgbm) {
    lightgbm_defaults <- list(
      nrounds = 100L,
      learning_rate = 0.1,
      num_leaves = 31L,
      min_data_in_leaf = 20L,
      seed = 1L
    )
    for (name in names(lightgbm_arguments)) {
      if (is.null(lightgbm_arguments[[name]])) {
        lightgbm_arguments[[name]] <- lightgbm_defaults[[name]]
      }
    }
    learner_params <- do.call(
      sb_lightgbm_binary_params,
      c(lightgbm_arguments, list(learner = learner))
    )
  } else if (any(!vapply(lightgbm_arguments, is.null, logical(1)))) {
    sb_abort(
      paste0(
        "LightGBM parameters apply only to fixed LightGBM learners."
      ),
      "stormbreaker_irrelevant_lightgbm_parameter"
    )
  }
  catboost_arguments <- list(
    iterations = options$catboost_iterations,
    learning_rate = options$catboost_learning_rate,
    depth = options$catboost_depth,
    l2_leaf_reg = options$catboost_l2_leaf_reg,
    seed = options$catboost_seed
  )
  if (is_catboost) {
    catboost_defaults <- list(
      iterations = 100L,
      learning_rate = 0.1,
      depth = 6L,
      l2_leaf_reg = 3,
      seed = 1L
    )
    for (name in names(catboost_arguments)) {
      if (is.null(catboost_arguments[[name]])) {
        catboost_arguments[[name]] <- catboost_defaults[[name]]
      }
    }
    parameter_constructor <- if (identical(learner, "catboost_cox")) {
      sb_catboost_cox_params
    } else if (identical(learner, "catboost_aft_survival")) {
      sb_catboost_aft_survival_params
    } else {
      sb_catboost_binary_params
    }
    learner_params <- do.call(parameter_constructor, catboost_arguments)
  } else if (any(!vapply(catboost_arguments, is.null, logical(1)))) {
    sb_abort(
      paste0(
        "CatBoost parameters apply only to a fixed CatBoost learner."
      ),
      "stormbreaker_irrelevant_catboost_parameter"
    )
  }
  gbm_arguments <- list(
    n_trees = options$gbm_n_trees,
    interaction_depth = options$gbm_interaction_depth,
    min_obs_node = options$gbm_min_obs_node,
    shrinkage = options$gbm_shrinkage,
    bag_fraction = options$gbm_bag_fraction,
    seed = options$gbm_seed
  )
  if (is_gbm_fixed) {
    gbm_defaults <- list(
      n_trees = 100L,
      interaction_depth = 1L,
      min_obs_node = 10L,
      shrinkage = 0.1,
      bag_fraction = 1,
      seed = 1L
    )
    for (name in names(gbm_arguments)) {
      if (is.null(gbm_arguments[[name]])) {
        gbm_arguments[[name]] <- gbm_defaults[[name]]
      }
    }
    learner_params <- do.call(sb_gbm_cox_params, gbm_arguments)
  } else if (any(!vapply(gbm_arguments, is.null, logical(1)))) {
    sb_abort(
      "GBM parameters apply only when `learner = 'gbm_cox'`.",
      "stormbreaker_irrelevant_gbm_parameter"
    )
  }
  ranger_arguments <- list(
    num_trees = options$ranger_num_trees,
    mtry = options$ranger_mtry,
    min_node_size = options$ranger_min_node_size,
    seed = options$ranger_seed
  )
  if (is_ranger) {
    ranger_defaults <- list(
      num_trees = 500L,
      min_node_size = if (is_ranger_survival) 3L else 10L,
      seed = 1L
    )
    for (name in names(ranger_defaults)) {
      if (is.null(ranger_arguments[[name]])) {
        ranger_arguments[[name]] <- ranger_defaults[[name]]
      }
    }
    parameter_constructor <- if (is_ranger_survival) {
      sb_ranger_survival_params
    } else {
      sb_ranger_binary_params
    }
    ranger_arguments$learner <- learner
    learner_params <- do.call(parameter_constructor, ranger_arguments)
  } else if (any(!vapply(ranger_arguments, is.null, logical(1)))) {
    sb_abort(
      paste0(
        "ranger parameters apply only to a fixed ranger learner."
      ),
      "stormbreaker_irrelevant_ranger_parameter"
    )
  }
  learner_params
}

sb_pipeline_fixed_params <- function(learner, learner_params, options) {
  is_coxboost <- identical(learner, "coxboost")
  is_plsrcox <- identical(learner, "plsrcox")
  is_superpc <- identical(learner, "superpc")
  is_survivalsvm <- learner %in% sb_survivalsvm_learner_ids()
  is_classical_binary_boosting <- learner %in%
    sb_classical_binary_boosting_learner_ids()
  is_classical_binary_fixed <- learner %in%
    sb_classical_binary_learner_ids()
  is_specialized_binary_linear <- learner %in%
    sb_specialized_binary_linear_learner_ids()
  is_omics_specialist_binary <- learner %in%
    sb_omics_specialist_binary_learner_ids()
  is_cancerclass <- identical(learner, "cancerclass_binary")
  is_xgboost <- learner %in% sb_xgboost_learner_ids()
  is_additional_survival <- learner %in%
    sb_additional_survival_learner_ids()
  is_flexsurv_parametric <- learner %in%
    sb_flexsurv_parametric_learner_ids()
  is_flexsurvspline <- learner %in% sb_flexsurvspline_learner_ids()
  is_mtlr_survival <- learner %in% sb_mtlr_survival_learner_ids()
  is_pamr_survival <- learner %in% sb_pamr_survival_learner_ids()
  is_timereg_survival <- learner %in% sb_timereg_survival_learner_ids()
  is_mboost_cox <- learner %in% sb_mboost_cox_learner_ids()
  is_mboost_aft <- learner %in% sb_mboost_aft_learner_ids()
  is_mboost_ranking <- learner %in% sb_mboost_ranking_learner_ids()
  is_abess <- learner %in% sb_abess_learner_ids()
  is_ncvreg <- learner %in% sb_ncvreg_learner_ids()
  is_rlt <- learner %in% sb_rlt_learner_ids()
  is_coxphf <- learner %in% sb_coxphf_learner_ids()
  is_bagged_knn_survival <- learner %in%
    sb_bagged_knn_survival_learner_ids()
  is_ahaz_survival <- learner %in% sb_ahaz_survival_learner_ids()
  is_bujar_linear <- learner %in% sb_bujar_linear_learner_ids()
  is_rule_aft <- learner %in% sb_rule_aft_learner_ids()
  is_wsrf <- learner %in% sb_wsrf_binary_learner_ids()
  is_slope_randomglm_binary <- learner %in%
    sb_slope_randomglm_binary_learner_ids()
  is_mgcv_gam <- learner %in% sb_mgcv_gam_learner_ids()
  is_mclust_da <- learner %in% sb_mclust_da_learner_ids()
  is_mixomics_da <- learner %in% sb_mixomics_da_learner_ids()
  is_hdda_binary <- learner %in% sb_hdda_binary_learner_ids()
  classical_arguments <- list(
    coxboost_stepno = options$coxboost_stepno,
    coxboost_penalty = options$coxboost_penalty,
    plsrcox_components = options$plsrcox_components,
    superpc_threshold = options$superpc_threshold,
    superpc_components = options$superpc_components,
    survivalsvm_gamma_mu = options$survivalsvm_gamma_mu
  )
  relevant_classical_arguments <- if (is_survivalsvm) {
    "survivalsvm_gamma_mu"
  } else {
    switch(learner,
      coxboost = c("coxboost_stepno", "coxboost_penalty"),
      plsrcox = "plsrcox_components",
      superpc = c("superpc_threshold", "superpc_components"),
      character()
    )
  }
  irrelevant_classical_arguments <- setdiff(
    names(classical_arguments),
    relevant_classical_arguments
  )
  has_irrelevant_classical_argument <- any(!vapply(
    classical_arguments[irrelevant_classical_arguments],
    is.null,
    logical(1)
  ))
  if (has_irrelevant_classical_argument) {
    sb_abort(
      paste0(
        "CoxBoost, plsRcox, SuperPC, and survival-SVM parameters apply ",
        "only to their matching learner."
      ),
      "stormbreaker_irrelevant_classical_survival_parameter"
    )
  }
  if (is_coxboost) {
    learner_params <- sb_coxboost_params(
      stepno = if (is.null(options$coxboost_stepno)) 100L else options$coxboost_stepno,
      penalty = options$coxboost_penalty
    )
  } else if (is_plsrcox) {
    learner_params <- sb_plsrcox_params(
      components = if (is.null(options$plsrcox_components)) {
        2L
      } else {
        options$plsrcox_components
      }
    )
  } else if (is_superpc) {
    learner_params <- sb_superpc_params(
      threshold = if (is.null(options$superpc_threshold)) 1 else options$superpc_threshold,
      components = if (is.null(options$superpc_components)) {
        1L
      } else {
        options$superpc_components
      }
    )
  } else if (is_survivalsvm) {
    learner_params <- sb_survivalsvm_params(
      gamma_mu = if (is.null(options$survivalsvm_gamma_mu)) {
        1
      } else {
        options$survivalsvm_gamma_mu
      },
      learner = learner
    )
  }

  classical_binary_boosting_arguments <- list(
    adaboost_iterations = options$adaboost_iterations,
    logitboost_iterations = options$logitboost_iterations
  )
  relevant_binary_boosting_argument <- switch(learner,
    adaboost_binary = "adaboost_iterations",
    logitboost_binary = "logitboost_iterations",
    character()
  )
  irrelevant_binary_boosting_arguments <- setdiff(
    names(classical_binary_boosting_arguments),
    relevant_binary_boosting_argument
  )
  if (any(!vapply(
    classical_binary_boosting_arguments[
      irrelevant_binary_boosting_arguments
    ],
    is.null,
    logical(1)
  ))) {
    sb_abort(
      paste0(
        "AdaBoost and LogitBoost iteration parameters apply only to their ",
        "matching learner."
      ),
      "stormbreaker_irrelevant_classical_binary_boosting_parameter"
    )
  }
  if (is_classical_binary_boosting) {
    learner_params <- sb_classical_binary_boosting_params(
      learner,
      classical_binary_boosting_arguments[[
        relevant_binary_boosting_argument
      ]]
    )
  }
  if (is_classical_binary_fixed) {
    learner_params <- sb_classical_binary_params(learner)
  }
  if (!identical(learner, "stepglm_binary") &&
    !is.null(options$stepglm_direction)) {
    sb_abort(
      paste0(
        "`stepglm_direction` applies only when ",
        "`learner = 'stepglm_binary'`."
      ),
      "stormbreaker_irrelevant_stepglm_parameter"
    )
  }
  if (is_specialized_binary_linear) {
    learner_params <- sb_specialized_binary_linear_params(
      learner,
      stepglm_direction = options$stepglm_direction
    )
  }
  if (is_omics_specialist_binary) {
    learner_params <- sb_omics_specialist_binary_params(learner)
  }
  if (is_cancerclass) {
    learner_params <- sb_cancerclass_binary_params()
  }
  if (is_xgboost) {
    learner_params <- sb_xgboost_params(learner)
  }
  if (is_additional_survival) {
    learner_params <- sb_additional_survival_params(learner)
  }
  if (is_flexsurv_parametric) {
    learner_params <- sb_flexsurv_parametric_params(learner)
  }
  if (is_flexsurvspline) {
    learner_params <- sb_flexsurvspline_params(learner)
  }
  if (is_mtlr_survival) {
    learner_params <- sb_mtlr_survival_params(learner)
  }
  if (is_pamr_survival) {
    learner_params <- sb_pamr_survival_params(learner)
  }
  if (is_timereg_survival) {
    learner_params <- sb_timereg_survival_params(learner)
  }
  if (is_mboost_cox) {
    learner_params <- sb_mboost_cox_params(learner)
  }
  if (is_mboost_aft) {
    learner_params <- sb_mboost_aft_params(learner)
  }
  if (is_mboost_ranking) {
    learner_params <- sb_mboost_ranking_params(learner)
  }
  if (is_abess) {
    learner_params <- sb_abess_learner_params(learner)
  }
  if (is_ncvreg) {
    learner_params <- sb_ncvreg_learner_params(learner)
  }
  if (is_rlt) {
    learner_params <- sb_rlt_learner_params(learner)
  }
  if (is_coxphf) {
    learner_params <- sb_coxphf_learner_params(learner)
  }
  if (learner %in% sb_akritas_learner_ids()) learner_params <- sb_akritas_params()
  if (is_bagged_knn_survival) {
    learner_params <- sb_bagged_knn_survival_params(learner)
  }
  if (is_ahaz_survival) {
    learner_params <- sb_ahaz_learner_params(learner)
  }
  if (is_bujar_linear) {
    learner_params <- sb_bujar_linear_params(learner)
  }
  if (is_rule_aft) {
    learner_params <- sb_rule_aft_params(learner)
  }
  if (is_wsrf) {
    learner_params <- sb_wsrf_binary_params(learner)
  }
  if (is_slope_randomglm_binary) {
    learner_params <- sb_slope_randomglm_binary_params(learner)
  }
  if (is_mgcv_gam) {
    learner_params <- sb_mgcv_gam_params(learner)
  }
  if (is_mclust_da) {
    learner_params <- sb_mclust_da_params(learner)
  }
  if (is_mixomics_da) {
    learner_params <- sb_mixomics_da_params(learner)
  }
  if (is_hdda_binary) {
    learner_params <- sb_hdda_binary_params(learner)
  }
  learner_params
}

# ---- pipeline-prediction-spec ----

sb_pipeline_prediction_spec <- function(pipeline_task_type, learner, learner_params) {
  is_gbm <- learner %in% c("gbm_cox", "gbm_cox_cv")
  is_survivalsvm_kernel <- learner %in%
    sb_survivalsvm_kernel_learner_ids()
  is_xgboost_cox <- learner %in% sb_xgboost_cox_learner_ids()
  is_xgboost_aft <- learner %in% sb_xgboost_aft_learner_ids()
  is_flexsurv_parametric <- learner %in%
    sb_flexsurv_parametric_learner_ids()
  is_flexsurvspline <- learner %in% sb_flexsurvspline_learner_ids()
  is_mtlr_survival <- learner %in% sb_mtlr_survival_learner_ids()
  is_pamr_survival <- learner %in% sb_pamr_survival_learner_ids()
  is_mboost_cox <- learner %in% sb_mboost_cox_learner_ids()
  is_mboost_aft <- learner %in% sb_mboost_aft_learner_ids()
  is_mboost_ranking <- learner %in% sb_mboost_ranking_learner_ids()
  is_bagged_knn_survival <- learner %in%
    sb_bagged_knn_survival_learner_ids()
  is_ahaz_survival <- learner %in% sb_ahaz_survival_learner_ids()
  is_bujar_linear <- learner %in% sb_bujar_linear_learner_ids()
  is_rule_aft <- learner %in% sb_rule_aft_learner_ids()
  is_ranger_survival <- learner %in% sb_ranger_survival_learner_ids()
  pipeline <- list()
  if (identical(pipeline_task_type, "binary")) {
    pipeline$prediction_type <- "probability"
    pipeline$positive_class <- 1L
    pipeline$prediction_direction <- "higher_is_more_likely_positive"
  } else if (is_survivalsvm_kernel) {
    pipeline$risk_type <- "calibrated_survival_svm_kernel_score"
    pipeline$risk_reference <- sb_survivalsvm_kernel_reference(
      learner_params$kernel
    )
    pipeline$risk_formula <-
      "negative_native_kernel_score_times_analysis_cox_coefficient"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "catboost_cox")) {
    pipeline$risk_type <- "cox_raw_margin"
    pipeline$risk_reference <- "catboost_native"
    pipeline$risk_formula <- "catboost_raw_formula_value"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "catboost_aft_survival")) {
    pipeline$risk_type <- "aft_log_time"
    pipeline$risk_reference <- "catboost_aft_normal_native"
    pipeline$risk_formula <- "negative_catboost_aft_log_time_margin"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_xgboost_cox) {
    pipeline$risk_type <- "cox_raw_margin"
    pipeline$risk_reference <- "xgboost_native"
    pipeline$risk_formula <- "xgboost_raw_margin"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (is_xgboost_aft) {
    pipeline$risk_type <- "aft_log_time"
    pipeline$risk_reference <- "xgboost_aft_normal_native"
    pipeline$risk_formula <- "negative_xgboost_aft_log_time_margin"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (is_gbm) {
    pipeline$risk_type <- "cox_log_hazard"
    pipeline$risk_reference <- "gbm_native"
    pipeline$risk_formula <- "gbm_link_prediction"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (is_mboost_cox) {
    pipeline$risk_type <- "cox_log_hazard"
    pipeline$risk_reference <- "mboost_native"
    pipeline$risk_formula <- "mboost_coxph_link_prediction"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (is_mboost_aft) {
    pipeline$risk_type <- "aft_log_time"
    pipeline$risk_reference <- paste0(
      "mboost_", learner_params$algorithm, "_",
      learner_params$distribution, "_aft_native"
    )
    pipeline$risk_formula <- "negative_mboost_aft_log_time_link"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_mboost_ranking) {
    pipeline$risk_type <- if (identical(learner_params$family, "Gehan")) {
      "gehan_aft_rank_score"
    } else {
      "uno_ipcw_concordance_score"
    }
    pipeline$risk_reference <- paste0(
      "mboost_", learner_params$algorithm, "_",
      tolower(learner_params$family), "_native"
    )
    pipeline$risk_formula <- "negative_mboost_ranking_link"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (identical(learner, "gam_cox_reml")) {
    pipeline$risk_type <- "cox_log_hazard"
    pipeline$risk_reference <- "mgcv_native"
    pipeline$risk_formula <- "mgcv_cox_ph_link_prediction"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (learner %in% sb_rlt_survival_learner_ids()) {
    pipeline$risk_type <- "reinforcement_survival_forest_rmst"
    pipeline$risk_reference <- "RLT_native_discrete_hazard"
    pipeline$risk_formula <- "negative_step_RMST_over_native_event_grid"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (learner %in% sb_akritas_learner_ids()) {
    pipeline$risk_type <- "akritas_conditional_survival_rmst"
    pipeline$risk_reference <- "maximum_analysis_followup"
    pipeline$risk_formula <- "negative_step_RMST_over_analysis_time_grid"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_bagged_knn_survival) {
    pipeline$risk_type <- "bagged_knn_survival_rmst"
    pipeline$risk_reference <- "bnnSurvival_native_training_time_grid"
    pipeline$risk_formula <-
      "negative_step_RMST_over_native_training_time_grid"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_ahaz_survival) {
    pipeline$risk_type <- "additive_hazard_linear_predictor"
    pipeline$risk_reference <- "ahaz_native"
    pipeline$risk_formula <- "x_beta"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (is_bujar_linear) {
    pipeline$risk_type <- "buckley_james_log_time"
    pipeline$risk_reference <- "bujar_native"
    pipeline$risk_formula <- "negative_predicted_log_time"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_rule_aft) {
    risk <- switch(learner,
      bart_aft_survival = list(
        type = "aft_log_time",
        reference = "BART_abart_normal_native",
        formula = "negative_posterior_mean_log_time",
        direction = "higher_is_shorter_predicted_survival"
      ),
      pre_survival_rule_ensemble = list(
        type = "cox_rule_ensemble_link",
        reference = "pre_native_rule_ensemble",
        formula = "pre_cox_link_prediction",
        direction = "higher_is_higher_hazard"
      ),
      penaft_gehan_elastic_net = list(
        type = "semiparametric_gehan_aft_log_time",
        reference = "penAFT_native",
        formula = "negative_centered_elastic_net_log_time_prediction",
        direction = "higher_is_shorter_predicted_survival"
      )
    )
    pipeline$risk_type <- risk$type
    pipeline$risk_reference <- risk$reference
    pipeline$risk_formula <- risk$formula
    pipeline$risk_direction <- risk$direction
  } else if (is_ranger_survival) {
    pipeline$risk_type <- "survival_forest_mortality"
    pipeline$risk_reference <- "ranger_native"
    pipeline$risk_formula <- "row_sums_cumulative_hazard"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "aalen_additive_survival")) {
    pipeline$risk_type <- "aalen_cumulative_hazard"
    pipeline$risk_reference <- "timereg_native_last_analysis_event_time"
    pipeline$risk_formula <- "cumulative_hazard_at_last_analysis_event"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "prop_odds_survival")) {
    pipeline$risk_type <- "proportional_odds_log_failure_odds"
    pipeline$risk_reference <- "zero"
    pipeline$risk_formula <- "x_beta"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (learner %in% sb_survreg_learner_ids()) {
    pipeline$risk_type <- "aft_log_time"
    pipeline$risk_reference <- paste0(
      learner_params$distribution,
      "_aft_native"
    )
    pipeline$risk_formula <- "negative_log_time_location"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (is_flexsurv_parametric) {
    pipeline$risk_type <- if (identical(
      learner_params$location_parameter,
      "mu"
    )) {
      "aft_log_time"
    } else if (identical(learner_params$distribution, "gompertz")) {
      "proportional_hazards_log_rate"
    } else {
      "accelerated_failure_time_log_rate"
    }
    pipeline$risk_reference <- sb_flexsurv_parametric_risk_reference(learner)
    pipeline$risk_formula <- sb_flexsurv_parametric_risk_formula(learner)
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_flexsurvspline) {
    pipeline$risk_type <- switch(learner_params$scale,
      hazard = "flexible_parametric_log_cumulative_hazard_location",
      odds = "flexible_parametric_log_cumulative_odds_location",
      normal = "flexible_parametric_normal_survival_location"
    )
    pipeline$risk_reference <- sb_flexsurvspline_risk_reference(learner)
    pipeline$risk_formula <- sb_flexsurvspline_risk_formula(learner)
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_mtlr_survival) {
    pipeline$risk_type <- "multitask_logistic_mean_survival"
    pipeline$risk_reference <- "MTLR_native_time_grid"
    pipeline$risk_formula <- "negative_native_mean_predicted_survival_time"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_pamr_survival) {
    pipeline$risk_type <- "nearest_shrunken_centroid_survival_group_posterior"
    pipeline$risk_reference <- "pamr_native_survival_groups"
    pipeline$risk_formula <- "native_posterior_for_early_survival_group"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "rpart_survival")) {
    pipeline$risk_type <- "survival_tree_log_rate"
    pipeline$risk_reference <- "rpart_native"
    pipeline$risk_formula <- "log_relative_event_rate"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "rsf_survival")) {
    pipeline$risk_type <- "survival_forest_mortality"
    pipeline$risk_reference <- "randomForestSRC_native"
    pipeline$risk_formula <- "mortality_sum_cumulative_hazard"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "ctree_survival")) {
    pipeline$risk_type <- "conditional_survival_cumulative_hazard"
    pipeline$risk_reference <- "partykit_native"
    pipeline$risk_formula <- "terminal_distribution_cumulative_hazard"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "cforest_survival")) {
    pipeline$risk_type <- "conditional_survival_cumulative_hazard"
    pipeline$risk_reference <- "partykit_native"
    pipeline$risk_formula <- "ensemble_distribution_cumulative_hazard"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "bart_survival")) {
    pipeline$risk_type <- "posterior_negative_restricted_mean_survival"
    pipeline$risk_reference <- "BART_analysis_fold_time_grid"
    pipeline$risk_formula <- "negative_step_RMST_from_posterior_mean_survival"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "aorsf_survival")) {
    pipeline$risk_type <- "oblique_survival_forest_mortality"
    pipeline$risk_reference <- "aorsf_native"
    pipeline$risk_formula <- "native_mortality_over_event_times"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else {
    pipeline$risk_type <- "linear_predictor"
    pipeline$risk_reference <- "zero"
    pipeline$risk_formula <- "x_beta"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  }
  pipeline
}

# ---- pipeline-selector-config ----

sb_pipeline_selector_config <- function(selector, admission, uses_inner_cv, options) {
  is_mime_selector <- selector %in% sb_mime_survival_selector_ids()
  is_mime_rsf_selector <- selector %in% sb_mime_rsf_selector_ids()
  is_mime_tuned_selector <- selector %in% sb_mime_tuned_selector_ids()
  is_mime_composite_selector <- selector %in%
    sb_mime_composite_selector_ids()
  requires_top_k <- identical(
    admission$selector_parameter[[1L]],
    "selector_n_features"
  )
  requires_fixed_features <- identical(
    admission$selector_parameter[[1L]],
    "selector_features_and_feature_set_id"
  )
  requires_stepcox <- identical(
    admission$selector_parameter[[1L]],
    "stepcox_direction"
  ) || is_mime_composite_selector
  if (!requires_top_k && !is.null(options$selector_n_features)) {
    sb_abort(
      "`selector_n_features` applies only to a Top-K selector.",
      "stormbreaker_irrelevant_selector_parameter"
    )
  }
  if (requires_top_k && is.null(options$selector_n_features)) {
    sb_abort(
      paste0(
        "`selector_n_features` is required when ",
        "a Top-K selector is requested."
      ),
      "stormbreaker_missing_selector_parameter"
    )
  }
  has_fixed_parameter <- !is.null(options$selector_features) ||
    !is.null(options$feature_set_id)
  if (!requires_fixed_features && has_fixed_parameter) {
    sb_abort(
      paste0(
        "`selector_features` and `feature_set_id` apply only when ",
        "`selector = 'fixed_features'`."
      ),
      "stormbreaker_irrelevant_selector_parameter"
    )
  }
  if (requires_fixed_features &&
    (is.null(options$selector_features) || is.null(options$feature_set_id))) {
    sb_abort(
      paste0(
        "`selector_features` and `feature_set_id` are both required when ",
        "`selector = 'fixed_features'`."
      ),
      "stormbreaker_missing_selector_parameter"
    )
  }
  if (!requires_stepcox && !is.null(options$stepcox_direction)) {
    sb_abort(
      paste0(
        "`stepcox_direction` applies only to `stepcox` or a sequential ",
        "MIME selector ending in StepCox."
      ),
      "stormbreaker_irrelevant_selector_parameter"
    )
  }
  if (requires_stepcox && is.null(options$stepcox_direction)) {
    sb_abort(
      "`stepcox_direction` is required when `selector = 'stepcox'`.",
      "stormbreaker_missing_selector_parameter"
    )
  }
  has_mime_selector_control <- !is.null(options$selector_inner_folds) ||
    !is.null(options$selector_seed) || !is.null(options$selector_rsf_nodesize)
  if (!is_mime_selector && has_mime_selector_control) {
    sb_abort(
      paste0(
        "`selector_inner_folds`, `selector_seed`, and ",
        "`selector_rsf_nodesize` apply only to MIME survival selectors."
      ),
      "stormbreaker_irrelevant_selector_parameter"
    )
  }
  if (is_mime_selector && is.null(options$selector_seed)) {
    sb_abort(
      "`selector_seed` is required for a MIME survival selector.",
      "stormbreaker_missing_selector_parameter"
    )
  }
  if (is_mime_rsf_selector && is.null(options$selector_rsf_nodesize)) {
    sb_abort(
      "`selector_rsf_nodesize` is required for an RSF selector.",
      "stormbreaker_missing_selector_parameter"
    )
  }
  if (!is_mime_rsf_selector && !is.null(options$selector_rsf_nodesize)) {
    sb_abort(
      "`selector_rsf_nodesize` applies only to an RSF selector.",
      "stormbreaker_irrelevant_selector_parameter"
    )
  }
  if (is_mime_tuned_selector && is.null(options$selector_inner_folds)) {
    sb_abort(
      paste0(
        "`selector_inner_folds` is required for a CoxBoost or Lasso ",
        "selector."
      ),
      "stormbreaker_missing_selector_parameter"
    )
  }
  if (!is_mime_tuned_selector && !is.null(options$selector_inner_folds)) {
    sb_abort(
      paste0(
        "`selector_inner_folds` applies only to a CoxBoost or Lasso ",
        "selector."
      ),
      "stormbreaker_irrelevant_selector_parameter"
    )
  }
  selector_params <- list()
  budget_options <- list()
  if (selector %in% sb_pamr_selector_ids()) selector_params <- sb_pamr_selector_params(selector)
  if (selector %in% sb_survival_weight_selector_ids()) selector_params <- sb_survival_weight_selector_params(selector)
  if (selector %in% sb_survival_boosting_selector_ids()) selector_params <- sb_survival_boosting_selector_params(selector)
  if (selector %in% sb_binary_model_selector_ids()) {
    selector_params <- sb_binary_model_selector_options(selector, options$selector_options)
  } else if (selector %in% sb_boruta_selector_ids()) {
    budget_options <- sb_selector_budget_options(selector, options$selector_options)
    selector_params <- utils::modifyList(selector_params, budget_options)
  } else if (!is.null(options$selector_options)) sb_abort(
    "selector_options applies only to binary model and Boruta selectors.", "stormbreaker_bad_selector_options")
  if (requires_top_k) {
    selector_params$n_features <- sb_check_positive_integer(
      options$selector_n_features,
      "selector_n_features",
      minimum = 1L
    )
  }
  if (requires_fixed_features) {
    selector_params <- list(
      features = sb_validate_fixed_features(options$selector_features),
      feature_set_id = sb_validate_feature_set_id(options$feature_set_id),
      n_features = as.integer(length(options$selector_features))
    )
  }
  if (requires_stepcox && !is_mime_selector) {
    selector_params <- list(
      direction = sb_validate_stepcox_direction(options$stepcox_direction)
    )
  }
  if (is_mime_rsf_selector) {
    selector_params <- list(
      nodesize = sb_check_positive_integer(
        options$selector_rsf_nodesize,
        "selector_rsf_nodesize",
        minimum = 1L
      ),
      seed = sb_check_positive_integer(
        options$selector_seed,
        "selector_seed",
        minimum = 1L
      )
    )
  } else if (is_mime_tuned_selector) {
    selector_params <- list(
      inner_folds = sb_check_positive_integer(
        options$selector_inner_folds,
        "selector_inner_folds",
        minimum = 3L
      ),
      seed = sb_check_positive_integer(
        options$selector_seed,
        "selector_seed",
        minimum = 1L
      )
    )
  }
  if (is_mime_composite_selector) {
    selector_params$direction <- sb_validate_stepcox_direction(
      options$stepcox_direction
    )
  }
  if (uses_inner_cv && "n_features" %in% names(selector_params) &&
    selector_params$n_features < 2L) {
    sb_abort(
      paste0(
        "A glmnet pipeline must retain at least 2 features; ",
        "glmnet cannot fit a path with only one variable."
      ),
      "stormbreaker_glmnet_too_few_selected_features"
    )
  }
  selector_suffix <- if (length(budget_options)) {
    paste0("__", paste(paste0(names(budget_options), unlist(budget_options)), collapse = "__"))
  } else if (selector %in% sb_binary_model_selector_ids()) {
    paste0("__", paste(vapply(names(selector_params), function(n) {
      v <- selector_params[[n]]
      paste0(n, if (is.numeric(v)) sb_format_roundtrip_double(as.numeric(v)) else v)
    }, character(1)), collapse = "__"))
  } else if (requires_top_k) {
    paste0("__k", selector_params$n_features)
  } else if (requires_fixed_features) {
    paste0("__set", selector_params$feature_set_id)
  } else if (is_mime_rsf_selector) {
    paste0(
      "__nodesize", selector_params$nodesize,
      "__seed", selector_params$seed,
      if (is_mime_composite_selector) {
        paste0("__direction", selector_params$direction)
      } else {
        ""
      }
    )
  } else if (is_mime_tuned_selector) {
    paste0(
      "__inner", selector_params$inner_folds,
      "__seed", selector_params$seed,
      if (is_mime_composite_selector) {
        paste0("__direction", selector_params$direction)
      } else {
        ""
      }
    )
  } else if (requires_stepcox) {
    paste0("__direction", selector_params$direction)
  } else {
    ""
  }
  list(params = selector_params, suffix = selector_suffix)
}

# ---- adapter-registry ----

sb_uniform_learner_adapter_families <- function() {
  list(
    list(sb_akritas_learner_ids(), sb_akritas_learner_adapter),
    list(sb_coxphf_learner_ids(), sb_coxphf_learner_adapter),
    list(sb_mboost_cox_learner_ids(), sb_mboost_cox_learner_adapter),
    list(sb_mboost_aft_learner_ids(), sb_mboost_aft_learner_adapter),
    list(
      sb_mboost_ranking_learner_ids(),
      sb_mboost_ranking_learner_adapter
    ),
    list(sb_xgboost_learner_ids(), sb_xgboost_learner_adapter),
    list(sb_mtlr_survival_learner_ids(), sb_mtlr_survival_learner_adapter),
    list(
      sb_pamr_survival_learner_ids(),
      sb_pamr_survival_learner_adapter
    ),
    list(
      sb_additional_survival_learner_ids(),
      sb_additional_survival_learner_adapter
    ),
    list(
      sb_flexsurv_parametric_learner_ids(),
      sb_flexsurv_parametric_learner_adapter
    ),
    list(
      sb_flexsurvspline_learner_ids(),
      sb_flexsurvspline_learner_adapter
    ),
    list(sb_timereg_survival_learner_ids(), sb_timereg_survival_learner_adapter),
    list(
      sb_classical_binary_learner_ids(),
      sb_classical_binary_learner_adapter
    ),
    list(
      sb_specialized_binary_linear_learner_ids(),
      sb_specialized_binary_linear_learner_adapter
    ),
    list(
      sb_omics_specialist_binary_learner_ids(),
      sb_omics_specialist_binary_learner_adapter
    ),
    list(sb_abess_learner_ids(), sb_abess_learner_adapter),
    list(sb_ncvreg_learner_ids(), sb_ncvreg_learner_adapter),
    list(sb_rlt_learner_ids(), sb_rlt_learner_adapter),
    list(
      sb_bagged_knn_survival_learner_ids(),
      sb_bagged_knn_survival_learner_adapter
    ),
    list(sb_ahaz_survival_learner_ids(), sb_ahaz_survival_learner_adapter),
    list(sb_bujar_linear_learner_ids(), sb_bujar_linear_learner_adapter),
    list(sb_rule_aft_learner_ids(), sb_rule_aft_learner_adapter),
    list(sb_wsrf_binary_learner_ids(), sb_wsrf_binary_learner_adapter),
    list(
      sb_slope_randomglm_binary_learner_ids(),
      sb_slope_randomglm_binary_learner_adapter
    ),
    list(sb_mclust_da_learner_ids(), sb_mclust_da_learner_adapter),
    list(sb_mixomics_da_learner_ids(), sb_mixomics_da_learner_adapter),
    list(sb_hdda_binary_learner_ids(), sb_hdda_binary_learner_adapter)
  )
}

sb_resolve_uniform_learner_adapter <- function(id, params) {
  for (family in sb_uniform_learner_adapter_families()) {
    if (id %in% family[[1L]]) {
      return(sb_constant_sum_learner_adapter(family[[2L]](id, params), family[[2L]]))
    }
  }
  NULL
}

sb_resolve_selector_adapter <- function(id, params) {
  if (id %in% sb_pamr_selector_ids()) return(sb_pamr_selector_adapter(id, params))
  if (id %in% sb_survival_weight_selector_ids()) return(sb_survival_weight_selector_adapter(id, params))
  if (id %in% sb_survival_boosting_selector_ids()) return(sb_survival_boosting_selector_adapter(id, params))
  if (identical(id, "ordered_chain")) return(sb_chain_selector_adapter(params))
  if (id %in% sb_binary_model_selector_ids()) return(sb_binary_model_selector_adapter(id, params))
  if (identical(id, "none")) {
    return(sb_none_selector_adapter())
  }
  if (identical(id, "univ_cox_top_k")) {
    return(sb_univ_cox_top_k_selector_adapter(params$n_features))
  }
  if (identical(id, "univ_binary_auc_top_k")) {
    return(sb_univ_binary_auc_top_k_selector_adapter(params$n_features))
  }
  if (identical(id, "stepcox")) {
    return(sb_stepcox_selector_adapter(params$direction))
  }
  if (identical(id, "rsf_minimal_depth")) {
    return(sb_rsf_minimal_depth_selector_adapter(
      params$nodesize,
      params$seed
    ))
  }
  if (identical(id, "coxboost_nonzero")) {
    return(sb_coxboost_nonzero_selector_adapter(
      params$inner_folds,
      params$seed
    ))
  }
  if (id %in% c("lasso_nonzero", "enet_cox_nonzero")) {
    return(sb_lasso_nonzero_selector_adapter(
      params$inner_folds,
      params$seed, id = id
    ))
  }
  if (identical(id, "rsf_then_stepcox")) {
    return(sb_rsf_then_stepcox_selector_adapter(
      params$nodesize,
      params$seed,
      params$direction
    ))
  }
  if (identical(id, "coxboost_then_stepcox")) {
    return(sb_coxboost_then_stepcox_selector_adapter(
      params$inner_folds,
      params$seed,
      params$direction
    ))
  }
  if (identical(id, "lasso_then_stepcox")) {
    return(sb_lasso_then_stepcox_selector_adapter(
      params$inner_folds,
      params$seed,
      params$direction
    ))
  }
  if (id %in% sb_abess_selector_ids()) {
    return(sb_abess_selector_adapter(id))
  }
  if (id %in% sb_forest_selector_ids()) {
    return(sb_forest_selector_adapter(id, options = params))
  }
  if (id %in% sb_high_dimensional_selector_ids()) {
    return(sb_high_dimensional_selector_adapter(id))
  }
  if (identical(id, "fixed_features")) {
    return(sb_fixed_features_selector_adapter(
      params$features,
      params$feature_set_id
    ))
  }
  sb_abort(
    paste0("No admitted selector adapter is registered for `", id, "`."),
    "stormbreaker_unresolved_selector_adapter"
  )
}

sb_resolve_learner_adapter <- function(id, params = list()) {
  id <- sb_canonical_learner_id(id)
  sb_reject_retired_learner(id)

  if (identical(id, "coxph")) {
    return(sb_coxph_learner_adapter())
  }
  if (identical(id, "glmnet_cox")) {
    return(list(
      id = "glmnet_cox",
      params = params,
      requires_tuning = TRUE
    ))
  }
  if (identical(id, "glmnet_logistic")) {
    return(list(
      id = "glmnet_logistic",
      params = params,
      requires_tuning = TRUE
    ))
  }
  if (identical(id, "logistic")) {
    return(sb_logistic_learner_adapter())
  }
  if (identical(id, "gaussian_naive_bayes")) {
    return(sb_gaussian_naive_bayes_learner_adapter(params$var_smoothing))
  }
  if (id %in% sb_lightgbm_binary_learner_ids()) {
    return(sb_lightgbm_binary_learner_adapter(
      nrounds = params$nrounds,
      learning_rate = params$learning_rate,
      num_leaves = params$num_leaves,
      min_data_in_leaf = params$min_data_in_leaf,
      seed = params$seed,
      learner = id
    ))
  }
  if (identical(id, "catboost_binary")) {
    return(sb_catboost_binary_learner_adapter(
      iterations = params$iterations,
      learning_rate = params$learning_rate,
      depth = params$depth,
      l2_leaf_reg = params$l2_leaf_reg,
      seed = params$random_seed
    ))
  }
  if (identical(id, "catboost_cox")) {
    return(sb_catboost_cox_learner_adapter(
      iterations = params$iterations,
      learning_rate = params$learning_rate,
      depth = params$depth,
      l2_leaf_reg = params$l2_leaf_reg,
      seed = params$random_seed
    ))
  }
  if (identical(id, "catboost_aft_survival")) {
    return(sb_catboost_aft_survival_learner_adapter(
      iterations = params$iterations,
      learning_rate = params$learning_rate,
      depth = params$depth,
      l2_leaf_reg = params$l2_leaf_reg,
      seed = params$random_seed
    ))
  }
  if (identical(id, "cancerclass_binary")) {
    return(sb_cancerclass_binary_learner_adapter())
  }
  if (identical(id, "gbm_cox")) {
    return(sb_gbm_cox_learner_adapter(
      n_trees = params$n_trees,
      interaction_depth = params$interaction_depth,
      min_obs_node = params$min_obs_node,
      shrinkage = params$shrinkage,
      bag_fraction = params$bag_fraction,
      seed = params$seed
    ))
  }
  if (identical(id, "gbm_cox_cv")) {
    return(sb_gbm_cox_cv_learner_adapter(
      inner_folds = params$inner_folds,
      seed = params$seed,
      maximum_trees = params$maximum_trees
    ))
  }
  if (identical(id, "coxboost")) {
    return(sb_coxboost_learner_adapter(
      stepno = params$stepno,
      penalty = params$penalty
    ))
  }
  if (identical(id, "coxboost_cv")) {
    return(sb_coxboost_cv_learner_adapter(
      inner_folds = params$inner_folds,
      seed = params$seed,
      maximum_steps = params$maximum_steps,
      start_penalty = params$start_penalty
    ))
  }
  if (identical(id, "plsrcox")) {
    return(sb_plsrcox_learner_adapter(
      components = params$components
    ))
  }
  if (identical(id, "plsrcox_cv")) {
    return(sb_plsrcox_cv_learner_adapter(
      inner_folds = params$inner_folds,
      seed = params$seed,
      maximum_components = params$maximum_components
    ))
  }
  if (identical(id, "superpc")) {
    return(sb_superpc_learner_adapter(
      threshold = params$threshold,
      components = params$components
    ))
  }
  if (id %in% sb_survivalsvm_learner_ids()) {
    return(sb_survivalsvm_learner_adapter(
      gamma_mu = params$gamma_mu,
      learner = id,
      params = params
    ))
  }
  if (id %in% sb_ranger_binary_learner_ids()) {
    return(sb_ranger_binary_learner_adapter(
      num_trees = params$num_trees,
      mtry = params$mtry,
      min_node_size = params$min_node_size,
      seed = params$seed,
      learner = id
    ))
  }
  if (id %in% sb_classical_binary_boosting_learner_ids()) {
    return(sb_classical_binary_boosting_learner_adapter(
      learner = id,
      iterations = params$iterations
    ))
  }
  if (identical(id, "superpc_cv")) {
    return(sb_superpc_cv_learner_adapter(
      inner_folds = params$inner_folds,
      seed = params$seed
    ))
  }
  if (id %in% sb_ranger_survival_learner_ids()) {
    return(sb_ranger_survival_learner_adapter(
      num_trees = params$num_trees,
      mtry = params$mtry,
      min_node_size = params$min_node_size,
      seed = params$seed,
      learner = id
    ))
  }
  adapter <- sb_resolve_uniform_learner_adapter(id, params)
  if (!is.null(adapter)) {
    return(adapter)
  }
  sb_abort(
    paste0("No admitted learner adapter is registered for `", id, "`."),
    "stormbreaker_unresolved_learner_adapter"
  )
}

sb_resolve_pipeline_runtime <- function(pipeline) {
  pipeline <- sb_upgrade_legacy_names(pipeline)

  if (!inherits(pipeline, "stormbreaker_pipeline")) {
    sb_abort(
      "`pipeline` must be created by `sb_pipeline()`.",
      "stormbreaker_bad_pipeline"
    )
  }
  runtime <- list(
    id = pipeline$id,
    selector = sb_resolve_selector_adapter(
      pipeline$selector,
      if (!is.null(pipeline$selector_params)) {
        pipeline$selector_params
      } else {
        list()
      }
    ),
    learner = sb_resolve_learner_adapter(
      pipeline$learner,
      if (!is.null(pipeline$learner_params)) {
        pipeline$learner_params
      } else {
        list()
      }
    )
  )
  class(runtime) <- "stormbreaker_pipeline_runtime"
  runtime
}

sb_validate_pipeline_runtime <- function(runtime, pipeline) {
  selector_params <- if (is.null(pipeline$selector_params)) {
    list()
  } else {
    pipeline$selector_params
  }
  learner_params <- if (is.null(pipeline$learner_params)) {
    list()
  } else {
    pipeline$learner_params
  }
  valid_runtime <- inherits(runtime, "stormbreaker_pipeline_runtime") &&
    identical(runtime$id, pipeline$id) &&
    is.list(runtime$selector) &&
    identical(runtime$selector$id, pipeline$selector) &&
    identical(runtime$selector$params, selector_params) &&
    is.logical(runtime$selector$needs_outcome) &&
    length(runtime$selector$needs_outcome) == 1L &&
    !is.na(runtime$selector$needs_outcome) &&
    is.function(runtime$selector$fit) &&
    is.function(runtime$selector$apply) &&
    is.list(runtime$learner) &&
    identical(runtime$learner$id, pipeline$learner)

  tuned_placeholder <- pipeline$learner %in% c(
    "glmnet_cox",
    "glmnet_logistic"
  )
  if (valid_runtime && tuned_placeholder) {
    valid_runtime <- identical(runtime$learner$params, learner_params) &&
      identical(runtime$learner$requires_tuning, TRUE)
  }
  if (valid_runtime && !tuned_placeholder) {
    required_functions <- c(
      "fit",
      "predict",
      "extract_state",
      "predict_state"
    )
    valid_runtime <- all(vapply(
      runtime$learner[required_functions],
      is.function,
      logical(1)
    ))
    if (valid_runtime && identical(pipeline$task_type, "survival")) {
      valid_runtime <- is.function(runtime$learner$predict_survival)
    }
    if (valid_runtime && (length(learner_params) > 0L ||
      "params" %in% names(runtime$learner))) {
      valid_runtime <- identical(runtime$learner$params, learner_params)
    }
  }
  if (!valid_runtime) {
    sb_abort(
      "The resolved pipeline runtime is incomplete or inconsistent.",
      "stormbreaker_bad_pipeline_runtime"
    )
  }
  runtime
}

# ---- learner-state ----

sb_bad_learner_state <- function(message) {
  sb_abort(message, "stormbreaker_bad_learner_state")
}

sb_validate_learner_state_x <- function(state, x, learner) {
  valid_state <- inherits(state, "stormbreaker_learner_state") &&
    is.list(state) && identical(state[["learner"]], learner) &&
    is.character(state[["features"]]) &&
    length(state[["features"]]) > 0L &&
    !anyNA(state[["features"]]) &&
    !anyDuplicated(state[["features"]])
  if (!valid_state) {
    sb_bad_learner_state("The saved learner state is incomplete or invalid.")
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), state[["features"]])) {
    sb_abort(
      "Prediction columns do not match the saved learner state.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  x
}

sb_extract_cox_learner_state <- function(model, learner, coefficients,
                                         predict_survival, context,
                                         params = list(), audit = list()) {
  valid_coefficients <- is.data.frame(coefficients) &&
    all(c("feature", "coefficient") %in% names(coefficients)) &&
    nrow(coefficients) > 0L &&
    is.character(coefficients[["feature"]]) &&
    !anyNA(coefficients[["feature"]]) &&
    !anyDuplicated(coefficients[["feature"]]) &&
    is.numeric(coefficients[["coefficient"]]) &&
    !anyNA(coefficients[["coefficient"]]) &&
    all(is.finite(coefficients[["coefficient"]]))
  if (!valid_coefficients) {
    sb_bad_learner_state(
      "The fitted Cox learner did not yield portable coefficients."
    )
  }

  baseline_hazard <- NULL
  prediction_horizon <- if ("prediction_horizon" %in% names(context)) {
    context[["prediction_horizon"]]
  } else {
    NULL
  }
  if (!is.null(prediction_horizon)) {
    x <- context[["x"]]
    probability_result <- predict_survival(
      model,
      x,
      context[["time"]],
      context[["event"]],
      x[1L, , drop = FALSE],
      prediction_horizon
    )
    baseline_hazard <- list(
      prediction_horizon = probability_result$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context[["time_unit"]]
      } else {
        NA_character_
      },
      cumulative_baseline_hazard =
        probability_result$cumulative_baseline_hazard,
      baseline_survival = probability_result$baseline_survival,
      source = probability_result$source,
      ties = probability_result$ties
    )
  }

  structure(
    list(
      learner = learner,
      features = coefficients[["feature"]],
      coefficient = unname(coefficients[["coefficient"]]),
      risk_reference = "zero",
      risk_formula = "x_beta",
      params = params,
      audit = audit,
      baseline_hazard = baseline_hazard
    ),
    class = c("stormbreaker_cox_learner_state", "stormbreaker_learner_state")
  )
}

sb_predict_cox_learner_state <- function(state, x) {
  valid_learner <- is.list(state) && state[["learner"]] %in%
    sb_linear_survival_learner_ids()
  if (!isTRUE(valid_learner)) {
    sb_bad_learner_state("The saved Cox learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, state[["learner"]])
  coefficient <- state[["coefficient"]]
  valid_coefficient <- is.numeric(coefficient) &&
    length(coefficient) == ncol(x) && !anyNA(coefficient) &&
    all(is.finite(coefficient))
  if (!valid_coefficient) {
    sb_bad_learner_state("The saved Cox coefficients are invalid.")
  }
  risk <- unname(as.vector(as.matrix(x) %*% coefficient))
  if (anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "The saved Cox state produced invalid risk predictions.",
      "stormbreaker_invalid_survival_prediction"
    )
  }
  result <- list(risk = risk)
  baseline <- state[["baseline_hazard"]]
  if (!is.null(baseline)) {
    valid_baseline <- is.list(baseline) &&
      is.numeric(baseline[["prediction_horizon"]]) &&
      length(baseline[["prediction_horizon"]]) == 1L &&
      !is.na(baseline[["prediction_horizon"]]) &&
      is.finite(baseline[["prediction_horizon"]]) &&
      baseline[["prediction_horizon"]] > 0
    if (!valid_baseline) {
      sb_bad_learner_state("The saved Cox baseline-hazard state is invalid.")
    }
    result$prediction_horizon <- baseline[["prediction_horizon"]]
    result$time_unit <- baseline[["time_unit"]]
    result$survival_probability <- sb_survival_probability_from_hazard(
      risk,
      baseline[["cumulative_baseline_hazard"]]
    )
  }
  result
}

sb_extract_logistic_learner_state <- function(learner, coefficients,
                                              params = list()) {
  valid_coefficients <- is.data.frame(coefficients) &&
    all(c("term", "coefficient") %in% names(coefficients)) &&
    nrow(coefficients) >= 2L &&
    identical(coefficients[["term"]][1L], "(Intercept)") &&
    is.character(coefficients[["term"]]) &&
    !anyNA(coefficients[["term"]]) &&
    !anyDuplicated(coefficients[["term"]]) &&
    is.numeric(coefficients[["coefficient"]]) &&
    !anyNA(coefficients[["coefficient"]]) &&
    all(is.finite(coefficients[["coefficient"]]))
  if (!valid_coefficients) {
    sb_bad_learner_state(
      "The fitted logistic learner did not yield portable coefficients."
    )
  }
  structure(
    list(
      learner = learner,
      features = coefficients[["term"]][-1L],
      intercept = unname(coefficients[["coefficient"]][1L]),
      coefficient = unname(coefficients[["coefficient"]][-1L]),
      positive_class = 1L,
      params = params
    ),
    class = c(
      "stormbreaker_logistic_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_logistic_learner_state <- function(state, x) {
  valid_learner <- is.list(state) && state[["learner"]] %in%
    c(
      "logistic", "glmnet_logistic", sb_abess_binary_learner_ids(),
      sb_ncvreg_binary_learner_ids()
    )
  if (!isTRUE(valid_learner)) {
    sb_bad_learner_state("The saved logistic learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, state[["learner"]])
  intercept <- state[["intercept"]]
  coefficient <- state[["coefficient"]]
  valid_coefficients <- is.numeric(intercept) && length(intercept) == 1L &&
    !is.na(intercept) && is.finite(intercept) &&
    is.numeric(coefficient) && length(coefficient) == ncol(x) &&
    !anyNA(coefficient) && all(is.finite(coefficient))
  if (!valid_coefficients) {
    sb_bad_learner_state("The saved logistic coefficients are invalid.")
  }
  probability <- unname(stats::plogis(
    intercept + as.vector(as.matrix(x) %*% coefficient)
  ))
  valid_probability <- length(probability) == nrow(x) &&
    !anyNA(probability) && all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid_probability) {
    sb_abort(
      "The saved logistic state produced invalid class-one probabilities.",
      "stormbreaker_invalid_logistic_prediction"
    )
  }
  list(probability = probability)
}

sb_extract_naive_bayes_learner_state <- function(model) {
  if (!inherits(model, "stormbreaker_gaussian_naive_bayes")) {
    sb_bad_learner_state("The fitted Gaussian Naive Bayes state is invalid.")
  }
  structure(
    list(
      learner = "gaussian_naive_bayes",
      features = model$features,
      payload = unclass(model)
    ),
    class = c(
      "stormbreaker_naive_bayes_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_naive_bayes_learner_state <- function(state, x) {
  x <- sb_validate_learner_state_x(state, x, "gaussian_naive_bayes")
  payload <- state[["payload"]]
  required <- c(
    "features", "classes", "var_smoothing", "epsilon", "counts",
    "priors", "means", "variances"
  )
  if (!is.list(payload) || length(setdiff(required, names(payload))) > 0L) {
    sb_bad_learner_state("The saved Gaussian Naive Bayes state is invalid.")
  }
  model <- structure(payload, class = "stormbreaker_gaussian_naive_bayes")
  list(probability = sb_predict_gaussian_naive_bayes(model, x))
}

sb_extract_lightgbm_learner_state <- function(model) {
  extracted <- sb_extract_lightgbm_tree_model(model)
  structure(
    c(list(learner = model$learner), extracted),
    class = c(
      "stormbreaker_lightgbm_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_lightgbm_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_lightgbm_binary_learner_ids()) {
    sb_bad_learner_state("The saved LightGBM learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  valid_state <- is.character(state[["model_string"]]) &&
    length(state[["model_string"]]) == 1L &&
    !is.na(state[["model_string"]]) &&
    nzchar(state[["model_string"]]) &&
    is.numeric(state[["nrounds"]]) && length(state[["nrounds"]]) == 1L &&
    !is.na(state[["nrounds"]]) && state[["nrounds"]] >= 1
  if (!valid_state) {
    sb_bad_learner_state("The saved LightGBM learner state is invalid.")
  }
  if (!sb_lightgbm_available()) {
    sb_abort(
      "Package `lightgbm` is required for fixed LightGBM prediction.",
      "stormbreaker_missing_lightgbm"
    )
  }
  booster <- tryCatch(
    lightgbm::lgb.load(model_str = state[["model_string"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved LightGBM state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_lightgbm_state",
        parent = error
      )
    }
  )
  model <- structure(
    list(
      learner = learner,
      features = state[["features"]],
      params = state[["params"]],
      nrounds = as.integer(state[["nrounds"]]),
      booster = booster
    ),
    class = "stormbreaker_lightgbm_binary"
  )
  list(probability = sb_predict_lightgbm_binary(model, x))
}

sb_extract_catboost_learner_state <- function(model) {
  extracted <- sb_extract_catboost_tree_model(model)
  structure(
    c(list(learner = "catboost_binary"), extracted),
    class = c(
      "stormbreaker_catboost_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_catboost_learner_state <- function(state, x) {
  x <- sb_validate_learner_state_x(state, x, "catboost_binary")
  valid_state <- is.raw(state[["model_raw"]]) &&
    length(state[["model_raw"]]) > 0L &&
    is.numeric(state[["iterations"]]) &&
    length(state[["iterations"]]) == 1L &&
    !is.na(state[["iterations"]]) && state[["iterations"]] >= 1
  if (!valid_state) {
    sb_bad_learner_state("The saved CatBoost learner state is invalid.")
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for fixed CatBoost prediction.",
      "stormbreaker_missing_catboost"
    )
  }
  booster <- tryCatch(
    unserialize(state[["model_raw"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved CatBoost state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_catboost_state",
        parent = error
      )
    }
  )
  model <- structure(
    list(
      features = state[["features"]],
      params = state[["params"]],
      iterations = as.integer(state[["iterations"]]),
      booster = booster
    ),
    class = "stormbreaker_catboost_binary"
  )
  list(probability = sb_predict_catboost_binary(model, x))
}

sb_extract_ranger_binary_learner_state <- function(model) {
  learner <- model[["learner"]]
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_ranger_binary_learner_ids()) {
    sb_bad_learner_state("The fitted ranger binary learner state is invalid.")
  }
  extracted <- sb_extract_ranger_binary_tree_model(model)
  structure(
    c(list(learner = learner), extracted),
    class = c(
      "stormbreaker_ranger_binary_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_ranger_binary_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_ranger_binary_learner_ids()) {
    sb_bad_learner_state("The saved ranger binary learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  valid_state <- is.raw(state[["model_raw"]]) &&
    length(state[["model_raw"]]) > 0L &&
    is.numeric(state[["num_trees"]]) &&
    length(state[["num_trees"]]) == 1L &&
    !is.na(state[["num_trees"]]) && state[["num_trees"]] >= 1 &&
    is.numeric(state[["effective_mtry"]]) &&
    length(state[["effective_mtry"]]) == 1L &&
    !is.na(state[["effective_mtry"]]) && state[["effective_mtry"]] >= 1 &&
    is.list(state[["params"]])
  if (!valid_state) {
    sb_bad_learner_state("The saved ranger binary learner state is invalid.")
  }
  if (!sb_ranger_available()) {
    sb_abort(
      "Package `ranger` is required for fixed binary forest prediction.",
      "stormbreaker_missing_ranger"
    )
  }
  forest <- tryCatch(
    unserialize(state[["model_raw"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved ranger binary state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_ranger_state",
        parent = error
      )
    }
  )
  if (!inherits(forest, "ranger.forest") ||
    !identical(forest$levels, c("0", "1"))) {
    sb_bad_learner_state("The saved ranger binary forest is invalid.")
  }
  model <- structure(
    list(
      learner = learner,
      features = state[["features"]],
      params = state[["params"]],
      num_trees = as.integer(state[["num_trees"]]),
      effective_mtry = as.integer(state[["effective_mtry"]]),
      forest = forest
    ),
    class = "stormbreaker_ranger_binary"
  )
  list(probability = sb_predict_ranger_binary(model, x))
}

sb_extract_classical_binary_boosting_learner_state <- function(model) {
  extracted <- sb_extract_classical_binary_boosting_model(model)
  structure(
    c(list(learner = model$learner), extracted),
    class = c(
      "stormbreaker_classical_binary_boosting_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_classical_binary_boosting_learner_state <- function(state, x) {
  learner <- state[["learner"]]
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_classical_binary_boosting_learner_ids()) {
    sb_bad_learner_state("The saved binary boosting learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  valid_state <- is.raw(state[["model_raw"]]) &&
    length(state[["model_raw"]]) > 0L &&
    is.numeric(state[["iterations"]]) &&
    length(state[["iterations"]]) == 1L &&
    !is.na(state[["iterations"]]) && state[["iterations"]] >= 1L &&
    is.numeric(state[["fitted_iterations"]]) &&
    length(state[["fitted_iterations"]]) == 1L &&
    !is.na(state[["fitted_iterations"]]) &&
    state[["fitted_iterations"]] >= 1L &&
    state[["fitted_iterations"]] <= state[["iterations"]] &&
    is.list(state[["params"]])
  if (!valid_state) {
    sb_bad_learner_state("The saved binary boosting learner state is invalid.")
  }
  dependency <- sb_classical_binary_boosting_dependency(learner)
  if (!dependency$available) {
    sb_abort(
      paste0(
        "Package `", dependency$package,
        "` is required for fixed binary boosting prediction."
      ),
      if (identical(learner, "adaboost_binary")) {
        "stormbreaker_missing_fastadaboost"
      } else {
        "stormbreaker_missing_catools"
      }
    )
  }
  engine <- tryCatch(
    unserialize(state[["model_raw"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved binary boosting state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_classical_binary_boosting_state",
        parent = error
      )
    }
  )
  engine <- sb_validate_classical_binary_boosting_engine(
    engine,
    learner,
    as.integer(state[["iterations"]])
  )
  fitted_iterations <- if (identical(learner, "adaboost_binary")) {
    length(engine$trees)
  } else {
    nrow(engine$Stump)
  }
  if (!identical(
    as.integer(state[["fitted_iterations"]]),
    as.integer(fitted_iterations)
  )) {
    sb_bad_learner_state("The saved binary boosting learner state is invalid.")
  }
  model <- structure(
    list(
      learner = learner,
      features = state[["features"]],
      params = state[["params"]],
      iterations = as.integer(state[["iterations"]]),
      fitted_iterations = as.integer(state[["fitted_iterations"]]),
      engine = engine
    ),
    class = "stormbreaker_classical_binary_boosting_fixed"
  )
  list(probability = sb_predict_classical_binary_boosting(model, x))
}

sb_extract_catboost_cox_learner_state <- function(model, context = list()) {
  extracted <- sb_extract_catboost_tree_model(model)
  baseline_hazard <- NULL
  prediction_horizon <- if ("prediction_horizon" %in% names(context)) {
    context[["prediction_horizon"]]
  } else {
    NULL
  }
  if (!is.null(prediction_horizon)) {
    probability_result <- sb_predict_catboost_cox_survival(
      model,
      context[["x"]],
      context[["time"]],
      context[["event"]],
      context[["x"]][1L, , drop = FALSE],
      prediction_horizon
    )
    baseline_hazard <- list(
      prediction_horizon = probability_result$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context[["time_unit"]]
      } else {
        NA_character_
      },
      cumulative_baseline_hazard =
        probability_result$cumulative_baseline_hazard,
      baseline_survival = probability_result$baseline_survival,
      source = probability_result$source,
      ties = probability_result$ties
    )
  }
  structure(
    c(
      list(
        learner = "catboost_cox",
        risk_reference = "catboost_native",
        risk_formula = "catboost_raw_formula_value"
      ),
      extracted,
      list(baseline_hazard = baseline_hazard)
    ),
    class = c(
      "stormbreaker_catboost_cox_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_catboost_cox_learner_state <- function(state, x) {
  x <- sb_validate_learner_state_x(state, x, "catboost_cox")
  valid_state <- is.raw(state[["model_raw"]]) &&
    length(state[["model_raw"]]) > 0L &&
    is.numeric(state[["iterations"]]) &&
    length(state[["iterations"]]) == 1L &&
    !is.na(state[["iterations"]]) && state[["iterations"]] >= 1
  if (!valid_state) {
    sb_bad_learner_state("The saved CatBoost Cox learner state is invalid.")
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for fixed CatBoost Cox prediction.",
      "stormbreaker_missing_catboost"
    )
  }
  booster <- tryCatch(
    unserialize(state[["model_raw"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved CatBoost Cox state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_catboost_state",
        parent = error
      )
    }
  )
  model <- structure(
    list(
      features = state[["features"]],
      params = state[["params"]],
      iterations = as.integer(state[["iterations"]]),
      booster = booster
    ),
    class = "stormbreaker_catboost_cox"
  )
  result <- list(risk = sb_predict_catboost_cox(model, x))
  baseline <- state[["baseline_hazard"]]
  if (!is.null(baseline)) {
    valid_baseline <- is.list(baseline) &&
      is.numeric(baseline[["prediction_horizon"]]) &&
      length(baseline[["prediction_horizon"]]) == 1L &&
      !is.na(baseline[["prediction_horizon"]]) &&
      is.finite(baseline[["prediction_horizon"]]) &&
      baseline[["prediction_horizon"]] > 0 &&
      is.character(baseline[["time_unit"]]) &&
      length(baseline[["time_unit"]]) == 1L &&
      !is.na(baseline[["time_unit"]]) && nzchar(baseline[["time_unit"]])
    if (!valid_baseline) {
      sb_bad_learner_state(
        "The saved CatBoost Cox baseline-hazard state is invalid."
      )
    }
    result$prediction_horizon <- baseline[["prediction_horizon"]]
    result$time_unit <- baseline[["time_unit"]]
    result$survival_probability <- sb_survival_probability_from_hazard(
      result$risk,
      baseline[["cumulative_baseline_hazard"]]
    )
  }
  result
}

sb_extract_gbm_cox_learner_state <- function(model, context = list()) {
  extracted <- sb_extract_gbm_cox_tree_model(model)
  learner <- if ("learner" %in% names(model)) {
    model[["learner"]]
  } else {
    "gbm_cox"
  }
  baseline_hazard <- NULL
  prediction_horizon <- if ("prediction_horizon" %in% names(context)) {
    context[["prediction_horizon"]]
  } else {
    NULL
  }
  if (!is.null(prediction_horizon)) {
    probability_result <- sb_predict_gbm_cox_survival(
      model,
      context[["x"]],
      context[["time"]],
      context[["event"]],
      context[["x"]][1L, , drop = FALSE],
      prediction_horizon
    )
    baseline_hazard <- list(
      prediction_horizon = probability_result$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context[["time_unit"]]
      } else {
        NA_character_
      },
      cumulative_baseline_hazard =
        probability_result$cumulative_baseline_hazard,
      baseline_survival = probability_result$baseline_survival,
      source = probability_result$source,
      ties = probability_result$ties
    )
  }
  structure(
    c(
      list(
        learner = learner,
        risk_reference = "gbm_native",
        risk_formula = "gbm_link_prediction"
      ),
      extracted,
      list(baseline_hazard = baseline_hazard)
    ),
    class = c(
      "stormbreaker_gbm_cox_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_gbm_cox_learner_state <- function(state, x) {
  learner <- state[["learner"]]
  if (!is.character(learner) || length(learner) != 1L ||
    !learner %in% c("gbm_cox", "gbm_cox_cv")) {
    sb_bad_learner_state("The saved GBM-Cox learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  valid_state <- is.raw(state[["model_raw"]]) &&
    length(state[["model_raw"]]) > 0L &&
    is.numeric(state[["n_trees"]]) &&
    length(state[["n_trees"]]) == 1L &&
    !is.na(state[["n_trees"]]) && state[["n_trees"]] >= 1 &&
    is.list(state[["params"]])
  if (!valid_state) {
    sb_bad_learner_state("The saved GBM-Cox learner state is invalid.")
  }
  if (!sb_gbm_available()) {
    sb_abort(
      "Package `gbm` is required for fixed GBM-Cox prediction.",
      "stormbreaker_missing_gbm"
    )
  }
  booster <- tryCatch(
    unserialize(state[["model_raw"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved GBM-Cox state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_gbm_state",
        parent = error
      )
    }
  )
  if (!inherits(booster, "gbm")) {
    sb_bad_learner_state("The saved GBM-Cox model is invalid.")
  }
  model <- structure(
    list(
      learner = learner,
      features = state[["features"]],
      params = state[["params"]],
      n_trees = as.integer(state[["n_trees"]]),
      booster = booster
    ),
    class = "stormbreaker_gbm_cox"
  )
  result <- list(risk = sb_predict_gbm_cox(model, x))
  baseline <- state[["baseline_hazard"]]
  if (!is.null(baseline)) {
    valid_baseline <- is.list(baseline) &&
      is.numeric(baseline[["prediction_horizon"]]) &&
      length(baseline[["prediction_horizon"]]) == 1L &&
      !is.na(baseline[["prediction_horizon"]]) &&
      is.finite(baseline[["prediction_horizon"]]) &&
      baseline[["prediction_horizon"]] > 0 &&
      is.character(baseline[["time_unit"]]) &&
      length(baseline[["time_unit"]]) == 1L &&
      !is.na(baseline[["time_unit"]]) && nzchar(baseline[["time_unit"]])
    if (!valid_baseline) {
      sb_bad_learner_state(
        "The saved GBM-Cox baseline-hazard state is invalid."
      )
    }
    result$prediction_horizon <- baseline[["prediction_horizon"]]
    result$time_unit <- baseline[["time_unit"]]
    result$survival_probability <- sb_survival_probability_from_hazard(
      result$risk,
      baseline[["cumulative_baseline_hazard"]]
    )
  }
  result
}

sb_extract_ranger_survival_learner_state <- function(model,
                                                     context = list()) {
  learner <- model[["learner"]]
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_ranger_survival_learner_ids()) {
    sb_bad_learner_state("The fitted ranger survival learner state is invalid.")
  }
  extracted <- sb_extract_ranger_survival_tree_model(model)
  survival_probability <- NULL
  prediction_horizon <- if ("prediction_horizon" %in% names(context)) {
    context[["prediction_horizon"]]
  } else {
    NULL
  }
  if (!is.null(prediction_horizon)) {
    probability_result <- sb_predict_ranger_survival_survival(
      model,
      context[["x"]],
      context[["time"]],
      context[["event"]],
      context[["x"]][1L, , drop = FALSE],
      prediction_horizon
    )
    survival_probability <- list(
      prediction_horizon = probability_result$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context[["time_unit"]]
      } else {
        NA_character_
      },
      source = probability_result$source,
      ties = probability_result$ties
    )
  }
  structure(
    c(
      list(
        learner = learner,
        risk_reference = "ranger_native",
        risk_formula = "row_sums_cumulative_hazard"
      ),
      extracted,
      list(
        baseline_hazard = NULL,
        survival_probability = survival_probability
      )
    ),
    class = c(
      "stormbreaker_ranger_survival_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_ranger_survival_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_ranger_survival_learner_ids()) {
    sb_bad_learner_state("The saved ranger survival learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  valid_state <- is.raw(state[["model_raw"]]) &&
    length(state[["model_raw"]]) > 0L &&
    is.numeric(state[["num_trees"]]) &&
    length(state[["num_trees"]]) == 1L &&
    !is.na(state[["num_trees"]]) && state[["num_trees"]] >= 1 &&
    is.numeric(state[["effective_mtry"]]) &&
    length(state[["effective_mtry"]]) == 1L &&
    !is.na(state[["effective_mtry"]]) && state[["effective_mtry"]] >= 1 &&
    is.list(state[["params"]])
  if (!valid_state) {
    sb_bad_learner_state("The saved ranger survival learner state is invalid.")
  }
  if (!sb_ranger_available()) {
    sb_abort(
      "Package `ranger` is required for fixed survival forest prediction.",
      "stormbreaker_missing_ranger"
    )
  }
  forest <- tryCatch(
    unserialize(state[["model_raw"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved ranger survival state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_ranger_state",
        parent = error
      )
    }
  )
  if (!inherits(forest, "ranger.forest")) {
    sb_bad_learner_state("The saved ranger survival forest is invalid.")
  }
  model <- structure(
    list(
      learner = learner,
      features = state[["features"]],
      params = state[["params"]],
      num_trees = as.integer(state[["num_trees"]]),
      effective_mtry = as.integer(state[["effective_mtry"]]),
      forest = forest
    ),
    class = "stormbreaker_ranger_survival"
  )
  curves <- sb_predict_ranger_survival_curves(model, x)
  result <- list(risk = curves$risk)
  probability_state <- state[["survival_probability"]]
  if (!is.null(probability_state)) {
    valid_probability_state <- is.list(probability_state) &&
      is.numeric(probability_state[["prediction_horizon"]]) &&
      length(probability_state[["prediction_horizon"]]) == 1L &&
      !is.na(probability_state[["prediction_horizon"]]) &&
      is.finite(probability_state[["prediction_horizon"]]) &&
      probability_state[["prediction_horizon"]] > 0 &&
      is.character(probability_state[["time_unit"]]) &&
      length(probability_state[["time_unit"]]) == 1L &&
      !is.na(probability_state[["time_unit"]]) &&
      nzchar(probability_state[["time_unit"]]) &&
      identical(
        probability_state[["source"]],
        "ranger_survival_forest_direct"
      ) &&
      identical(probability_state[["ties"]], "not_applicable")
    if (!valid_probability_state) {
      sb_bad_learner_state(
        "The saved ranger survival-probability state is invalid."
      )
    }
    result$prediction_horizon <-
      probability_state[["prediction_horizon"]]
    result$time_unit <- probability_state[["time_unit"]]
    result$survival_probability <-
      sb_ranger_survival_probability_at_horizon(
        curves,
        probability_state[["prediction_horizon"]]
      )
  }
  result
}

# ---- user-models ----

#' @rdname sb_models
#' @export
sb_model_catalog <- function(task_type) {
  if (!is.character(task_type) || length(task_type) != 1L || is.na(task_type) ||
      !task_type %in% c("binary", "survival")) {
    sb_abort("Choose task_type explicitly: 'binary' or 'survival'.",
      "stormbreaker_bad_task_type")
  }
  registry <- sb_pipeline_registry()
  rows <- registry[registry$task_type == task_type &
    registry$selector == "none" & registry$admission_status == "admitted", , drop = FALSE]
  rows <- rows[!duplicated(rows$learner), , drop = FALSE]
  rows$model_id <- rows$learner
  rows$glmnet_alpha <- ifelse(grepl("^glmnet_", rows$learner), 1, NA_real_)
  g <- rows[rows$learner == if (task_type == "binary") "glmnet_logistic" else "glmnet_cox", , drop = FALSE]
  ridge <- elastic <- g
  ridge$model_id <- paste0("ridge_", if (task_type == "binary") "logistic" else "cox")
  elastic$model_id <- paste0("elastic_net_", if (task_type == "binary") "logistic" else "cox")
  ridge$glmnet_alpha <- 0
  elastic$glmnet_alpha <- 0.5
  rows <- rbind(rows, ridge, elastic)
  rownames(rows) <- NULL
  rows <- rows[c("model_id", "learner", "glmnet_alpha",
    setdiff(names(rows), c("model_id", "learner", "glmnet_alpha")))]
  rows$model_label <- vapply(seq_len(nrow(rows)), function(i) {
    parameters <- if (is.na(rows$glmnet_alpha[i])) list() else list(alpha = rows$glmnet_alpha[i])
    sb_display_learner(rows$learner[i], parameters, rows$model_id[i])
  }, character(1))
  sb_add_dependency_status(rows)
}

#' @rdname sb_models
#' @export
sb_model_pipelines <- function(task_type, models = NULL,
                               inner_folds = 5L, inner_seed = 1701L) {
  menu <- sb_model_catalog(task_type)
  if (is.null(models)) models <- menu$model_id
  if (!is.character(models) || !length(models) || anyNA(models) ||
      anyDuplicated(models) || any(!models %in% menu$model_id)) {
    sb_abort("Supply distinct model IDs from sb_model_catalog() for this task.",
      "stormbreaker_bad_model_choice")
  }
  selected <- menu[match(models, menu$model_id), , drop = FALSE]
  pipelines <- lapply(seq_len(nrow(selected)), function(i) {
    args <- list(selector = "none", learner = selected$learner[i])
    if (!is.na(selected$glmnet_alpha[i])) {
      args$glmnet_alpha <- selected$glmnet_alpha[i]
      args$inner_folds <- inner_folds
      args$inner_seed <- inner_seed
    }
    do.call(sb_pipeline, args)
  })
  names(pipelines) <- models
  pipelines
}

#' @rdname sb_complete_run
#' @export
sb_compare_models <- function(task, models = "recommended", folds = 5L, seed = 1L,
                               inner_folds = 5L, inner_seed = 1701L,
                               view = FALSE, report_file = NULL) {
  sb_check_page_request(view, report_file)
  started <- proc.time()[["elapsed"]]
  type <- if (inherits(task, "stormbreaker_binary_task")) "binary" else
    if (inherits(task, "stormbreaker_surv_task")) "survival" else
      sb_abort("Create the task with sb_binary_task() or sb_surv_task().",
        "stormbreaker_bad_task")
  if (identical(models, "recommended")) models <- if (type == "binary")
    c("logistic", "glmnet_logistic", "random_forest_binary") else
    c("coxph", "glmnet_cox", "ranger_survival")
  pipelines <- sb_model_pipelines(type, models, inner_folds, inner_seed)
  compared <- sb_compare_cv(task, pipelines, folds = folds, seed = seed)
  ids <- vapply(pipelines, function(pipeline) pipeline$id, character(1))
  compared$comparison$model_id <- names(ids)[match(compared$comparison$pipeline_id, ids)]
  compared$model_pipelines <- pipelines
  compared$input_genes <- task$predictors
  compared$requested_models <- names(pipelines)
  compared$elapsed_seconds <- unname(proc.time()[["elapsed"]] - started)
  sb_attach_result_page(compared, view, report_file)
}
