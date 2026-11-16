# ---- component-fixtures ----

make_component_survival_task <- function() {
  set.seed(8901L)
  x <- as.data.frame(matrix(stats::rnorm(90L * 12L), nrow = 90L))
  names(x) <- sprintf("gene_%02d", seq_len(12L))
  lp <- x[[1]] - x[[2]] + .6 * x[[3]]
  et <- stats::rexp(90L, rate = exp(lp) / 15)
  ct <- stats::rexp(90L, rate = 1 / 25)
  d <- data.frame(id = seq_len(90L), time = pmax(pmin(et, ct), .05), event = as.integer(et <= ct), x)
  sb_surv_task(d, "id", "time", "event", names(x), "month")
}

# ---- curve-views ----

make_sb_curve_view_fixture <- function(task_type = "binary", comparison = FALSE) {
  thresholds <- c(0.1, 0.2, 0.3, 0.5, 0.7)
  predictors <- c("gene_a", "gene_b")
  if (identical(task_type, "binary")) {
    d <- sb_example_binary(120, seed = 1)
    task <- sb_binary_task(d, "patient_id", "outcome", predictors,
                           decision_thresholds = thresholds)
    pipelines <- list(sb_pipeline(learner = "logistic"))
    if (comparison) {
      pipelines[[2]] <- sb_pipeline(learner = "gaussian_naive_bayes")
    }
    external <- sb_example_binary(100, seed = 2)
    external$patient_id <- paste0("EXT_", external$patient_id)
    endpoint <- sb_binary_endpoint(external, "patient_id", "outcome")
  } else {
    d <- sb_example_survival(120, seed = 1)
    task <- sb_surv_task(d, "patient_id", "time", "event", predictors, "month",
                         prediction_horizon = 6, decision_thresholds = thresholds)
    pipelines <- list(sb_pipeline())
    if (comparison) {
      pipelines[[2]] <- sb_pipeline(selector = "univ_cox_top_k",
                                    selector_n_features = 1)
    }
    external <- sb_example_survival(100, seed = 2)
    external$patient_id <- paste0("EXT_", external$patient_id)
    endpoint <- sb_surv_endpoint(external, "patient_id", "time", "event", "month",
                                 prediction_horizon = 6)
  }
  predictions <- lapply(pipelines, function(pipeline) {
    fit <- sb_fit_final(task, pipeline)
    sb_predict_final(fit, external[c("patient_id", predictors)],
                     id_col = "patient_id")
  })
  if (comparison) {
    cv <- sb_compare_cv(task, pipelines, folds = 3, seed = 1)
    evaluated <- sb_compare_external(predictions, endpoint)
  } else {
    cv <- sb_run_cv(task, pipelines[[1]], folds = 3, seed = 1)
    evaluated <- sb_evaluate_external(predictions[[1]], endpoint)
  }
  list(cv = cv, external = evaluated, task = task)
}

# ---- fixtures ----

make_sb_fixture <- function(n = 60L, seed = 20260812L) {
  withr::with_seed(seed, {
    event <- sample(rep(c(0L, 1L), length.out = n))
    data.frame(
      patient_id = sprintf("P%03d", seq_len(n)),
      time = round(stats::rexp(n, rate = 0.08) + 0.1, 3),
      event = event,
      gene_a = stats::rnorm(n),
      gene_b = stats::rnorm(n),
      gene_c = stats::rnorm(n),
      check.names = FALSE
    )
  })
}

make_sb_binary_fixture <- function(n = 60L, seed = 20260818L) {
  withr::with_seed(seed, {
    outcome <- sample(rep(c(0L, 1L), length.out = n))
    data.frame(
      patient_id = sprintf("B%03d", seq_len(n)),
      outcome = outcome,
      gene_a = stats::rnorm(n, mean = 0.5 * outcome),
      gene_b = stats::rnorm(n),
      gene_c = stats::rnorm(n),
      check.names = FALSE
    )
  })
}

make_sb_high_dim_binary_fixture <- function(n = 60L, p = 80L,
                                            seed = 20260819L) {
  withr::with_seed(seed, {
    outcome <- sample(rep(c(0L, 1L), length.out = n))
    expression <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    expression[, 1L] <- expression[, 1L] + 0.8 * outcome
    expression[, 2L] <- expression[, 2L] - 0.5 * outcome
    colnames(expression) <- sprintf("gene_%03d", seq_len(p))
    data.frame(
      patient_id = sprintf("HD%03d", seq_len(n)),
      outcome = outcome,
      as.data.frame(expression, check.names = FALSE),
      check.names = FALSE
    )
  })
}

capture_sb_error <- function(code) {
  tryCatch(
    {
      force(code)
      NULL
    },
    stormbreaker_error = function(error) error
  )
}

population_variance <- function(x) {
  mean((x - mean(x))^2)
}

manual_gaussian_nb_probability <- function(x, priors, means, variances) {
  vapply(seq_len(nrow(x)), function(row) {
    log_joint <- vapply(c("0", "1"), function(class) {
      log(priors[[class]]) + sum(stats::dnorm(
        x[row, ],
        mean = means[class, ],
        sd = sqrt(variances[class, ]),
        log = TRUE
      ))
    }, numeric(1))
    unname(stats::plogis(log_joint[["1"]] - log_joint[["0"]]))
  }, numeric(1))
}

# ---- selector-routes ----

expected_fixed_learner_selectors <- function(task) {
  if (task == "binary") {
    c("none", "univ_binary_auc_top_k", "abess_binary_gic_nonzero",
      "pamr_binary_nonzero",
      "lasso_binary_nonzero", "stepglm_binary_aic", "glmboost_binary_nonzero",
      "rf_binary_importance", "enet_binary_nonzero", "gbm_binary_importance",
      "svm_rfe_binary", "gamboost_binary_nonzero", "xgboost_binary_importance",
      "plsrglm_binary_weights", "boruta_binary_roughfix", "varselrf_binary_oob_1se",
      "sis_binary_nlogn", "vsurf_binary_prediction", "mrmr_binary_nlogn",
      "jmi_binary_nlogn", "cmim_binary_nlogn", "disr_binary_nlogn",
      "jmim_binary_nlogn", "lpc_binary_nlogn", "stabsel_binary_lasso",
      "fixed_features")
  } else {
    c("none", "univ_cox_top_k", "stepcox", "abess_cox_gic_nonzero",
      "pamr_survival_nonzero",
      "ridge_cox_weights", "plsrcox_weights", "svm_survival_weights",
      "gbm_cox_importance", "glmboost_cox_nonzero", "superpc_survival_scores",
      "gamboost_cox_nonzero", "xgboost_cox_importance", "boruta_survival_roughfix",
      "lpc_cox_nlogn", "stabsel_cox_lasso", "lasso_nonzero", "enet_cox_nonzero",
      "coxboost_nonzero", "rsf_minimal_depth", "fixed_features")
  }
}

expect_route_dependencies <- function(rows, learner_packages) {
  selector_packages <- c(
    none = "", univ_binary_auc_top_k = "", univ_cox_top_k = "",
    stepcox = "", fixed_features = "", stepglm_binary_aic = "",
    abess_binary_gic_nonzero = "abess", abess_cox_gic_nonzero = "abess",
    pamr_binary_nonzero = "pamr", pamr_survival_nonzero = "pamr",
    lasso_binary_nonzero = "glmnet", enet_binary_nonzero = "glmnet",
    ridge_cox_weights = "glmnet", lasso_nonzero = "glmnet", enet_cox_nonzero = "glmnet",
    glmboost_binary_nonzero = "mboost", gamboost_binary_nonzero = "mboost",
    glmboost_cox_nonzero = "mboost", gamboost_cox_nonzero = "mboost",
    rf_binary_importance = "randomForest", gbm_binary_importance = "gbm",
    gbm_cox_importance = "gbm", svm_rfe_binary = "e1071",
    svm_survival_weights = "survivalsvm", plsrcox_weights = "plsRcox",
    plsrglm_binary_weights = "plsRglm", superpc_survival_scores = "superpc",
    xgboost_binary_importance = "xgboost", xgboost_cox_importance = "xgboost",
    boruta_binary_roughfix = "Boruta", boruta_survival_roughfix = "Boruta",
    varselrf_binary_oob_1se = "varSelRF", sis_binary_nlogn = "SIS",
    vsurf_binary_prediction = "VSURF", mrmr_binary_nlogn = "praznik",
    jmi_binary_nlogn = "praznik", cmim_binary_nlogn = "praznik",
    disr_binary_nlogn = "praznik", jmim_binary_nlogn = "praznik",
    lpc_binary_nlogn = "lpc", lpc_cox_nlogn = "lpc",
    stabsel_binary_lasso = "stabs + glmnet", stabsel_cox_lasso = "stabs + glmnet",
    coxboost_nonzero = "CoxBoost", rsf_minimal_depth = "randomForestSRC"
  )
  for (i in seq_len(nrow(rows))) {
    expected <- paste(selector_packages[[rows$selector[i]]],
      learner_packages[[rows$learner[i]]], sep = " + ")
    packages <- strsplit(expected, " + ", fixed = TRUE)[[1L]]
    expect_setequal(strsplit(rows$optional_dependency[i], " + ", fixed = TRUE)[[1L]],
      unique(packages[nzchar(packages)]))
  }
}

# ---- survival-svm-diagnostic ----

make_survival_svm_diagnostic_task <- function() {
  set.seed(11032L)
  x <- as.data.frame(matrix(rnorm(72L * 4L), nrow = 72L))
  names(x) <- c("gene A", "gene-B", "gene_c", "gene_d")
  lp <- .7 * x[[1]] - .5 * x[[2]] + .35 * x[[3]]^2 - .25 * x[[1]] * x[[4]]
  et <- rexp(72L, rate = .07 * exp(lp))
  ct <- rexp(72L, rate = .03)
  d <- data.frame(id = seq_len(72L), time = pmax(pmin(et, ct), .01), event = as.integer(et <= ct), x, check.names = FALSE)
  sb_surv_task(d, "id", "time", "event", names(x), "month")
}

# SHAP fixtures shared by the plotting scenarios.
make_sb_shap_plot_fixture <- function() {
  m <- rbind(c(a = 1, b = -0.3, c = 0.2), c(a = -0.5, b = 0.4, c = -0.1))
  rownames(m) <- c("patient-secret-1", "patient-secret-2")
  structure(list(status = "success", pipeline_id = "logistic", score_scale = "log_odds",
    reference_value = -0.4, contributions = m,
    scores = data.frame(patient_id = rownames(m), score = -0.4 + rowSums(m)),
    background = list(n = 10L)), class = "stormbreaker_additive_shap")
}
