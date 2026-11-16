# ---- pipeline-registry ----

sb_pipeline_registry_data <- function() {
  paired_rows <- function(task_type, learners, selector, dependencies,
                          tuning, base_semantics, selected_semantics) {
    data.frame(
      source = rep("stormbreaker_builtin", length(learners) * 2L),
      task_type = rep(task_type, length(learners) * 2L),
      selector = rep(c("none", selector), times = length(learners)),
      learner = rep(learners, each = 2L),
      admission_status = rep("admitted", length(learners) * 2L),
      selector_parameter = rep(
        c(NA, "selector_n_features"),
        times = length(learners)
      ),
      tuning_semantics = rep(tuning, each = 2L),
      optional_dependency = rep(dependencies, each = 2L),
      parameter_semantics = as.vector(rbind(
        base_semantics,
        selected_semantics
      )),
      stringsAsFactors = FALSE
    )
  }

  survival_learners <- c(
    "coxph", "glmnet_cox", "catboost_cox",
    "catboost_aft_survival", "gbm_cox",
    sb_mboost_cox_learner_ids(),
    sb_mboost_aft_learner_ids(),
    sb_mboost_ranking_learner_ids(),
    sb_ranger_survival_learner_ids(), sb_rlt_survival_learner_ids(),
    "coxboost", "plsrcox",
    "superpc", sb_survivalsvm_learner_ids(),
    sb_xgboost_cox_learner_ids(),
    sb_xgboost_aft_learner_ids(),
    sb_survreg_learner_ids(),
    sb_flexsurv_parametric_learner_ids(),
    sb_flexsurvspline_learner_ids(),
    sb_mtlr_survival_learner_ids(),
    sb_pamr_survival_learner_ids(), "rpart_survival",
    "rsf_survival",
    "ctree_survival", "cforest_survival", "bart_survival",
    "aorsf_survival", sb_timereg_survival_learner_ids(),
    sb_ncvreg_survival_learner_ids(),
    "abess_cox_gic",
    sb_coxphf_learner_ids(), sb_bagged_knn_survival_learner_ids(),
    sb_ahaz_survival_learner_ids(), sb_bujar_linear_learner_ids(), sb_akritas_learner_ids()
  )
  survivalsvm_base <- vapply(
    sb_survivalsvm_learner_ids(),
    function(learner) {
      kernel_type <- sb_survivalsvm_kernel(learner)
      kernel <- switch(kernel_type,
        lin_kernel = "linear",
        add_kernel = "additive",
        rbf_kernel = "radial-basis-function"
      )
      kernel_rule <- if (identical(kernel_type, "rbf_kernel")) {
        paste0(
          "; analysis-median-squared-distance kernel-parameter rule ",
          "fixed before resampling"
        )
      } else {
        ""
      }
      switch(sb_survivalsvm_type(learner),
        regression = paste0(
          kernel,
          " survival-SVM regression gamma fixed before resampling",
          kernel_rule
        ),
        vanbelle1 = paste0(
          kernel,
          " survival-SVM vanbelle1 ranking gamma and makediff3 rule fixed ",
          "before resampling",
          kernel_rule
        ),
        vanbelle2 = paste0(
          kernel,
          " survival-SVM vanbelle2 ranking gamma and makediff3 rule fixed ",
          "before resampling",
          kernel_rule
        ),
        hybrid = paste0(
          kernel,
          " survival-SVM hybrid regression/ranking gamma pair and ",
          "makediff3 rule fixed before resampling",
          kernel_rule
        )
      )
    },
    character(1)
  )
  survival_base <- c(
    "fixed ordinary Cox fit",
    paste0(
      "alpha, standardization setting, lambda rule, and tuning measure fixed; ",
      "lambda selected by inner CV"
    ),
    "all learner parameters fixed before resampling",
    paste0(
      "CatBoost SurvivalAft Normal loss, unit scale, tree controls, and seed ",
      "fixed before resampling"
    ),
    "all learner parameters fixed before resampling",
    paste0(
      "100 CoxPH boosting iterations, shrinkage 0.1, and spline base ",
      "learners fixed before resampling"
    ),
    paste0(
      "100 CoxPH boosting iterations, shrinkage 0.1, centering, and linear ",
      "base learners fixed before resampling"
    ),
    paste0(
      "100 CoxPH boosting iterations, shrinkage 0.1, and tree controls ",
      "fixed before resampling"
    ),
    vapply(
      sb_mboost_aft_learner_ids(),
      sb_mboost_aft_parameter_semantics,
      character(1)
    ),
    vapply(
      sb_mboost_ranking_learner_ids(),
      sb_mboost_ranking_parameter_semantics,
      character(1)
    ),
    vapply(
      sb_ranger_survival_learner_ids(),
      sb_ranger_parameter_semantics,
      character(1)
    ),
    vapply(
      sb_rlt_survival_learner_ids(),
      sb_rlt_parameter_semantics,
      character(1)
    ),
    "boosting steps and penalty rule fixed before resampling",
    "PLS component count fixed before resampling",
    "SuperPC threshold and component count fixed before resampling",
    survivalsvm_base,
    "XGBoost-Cox gbtree controls and seed fixed before resampling",
    paste0(
      "XGBoost-Cox DART tree and dropout controls and seed fixed before ",
      "resampling"
    ),
    paste0(
      "XGBoost-Cox coordinate-descent linear booster penalty, update order, ",
      "and seed fixed before resampling"
    ),
    paste0(
      "XGBoost-AFT normal-loss gbtree controls, scale, and seed fixed ",
      "before resampling"
    ),
    paste0(
      "XGBoost-AFT normal-loss DART tree/dropout controls, scale, and seed ",
      "fixed before resampling"
    ),
    paste0(
      "XGBoost-AFT normal-loss coordinate-descent linear booster penalty, ",
      "scale, update order, and seed fixed before resampling"
    ),
    "Weibull AFT distribution fixed; scale estimated in each analysis fold",
    "exponential AFT distribution and unit scale fixed before resampling",
    "log-normal AFT distribution fixed; scale estimated in each analysis fold",
    "log-logistic AFT distribution fixed; scale estimated in each analysis fold",
    vapply(
      sb_flexsurv_parametric_learner_ids(),
      sb_flexsurv_parametric_parameter_semantics,
      character(1)
    ),
    vapply(
      sb_flexsurvspline_learner_ids(),
      sb_flexsurvspline_parameter_semantics,
      character(1)
    ),
    sb_mtlr_survival_parameter_semantics(),
    sb_pamr_survival_parameter_semantics(),
    "survival-tree controls fixed before resampling",
    "random survival forest controls and seed fixed before resampling",
    "conditional-inference survival-tree controls fixed before resampling",
    "conditional-inference survival-forest controls fixed before resampling",
    paste0(
      "discrete-time probit BART trees, time grid, posterior draws, burn-in, ",
      "and seed fixed before resampling"
    ),
    paste0(
      "oblique survival-forest trees, split controls, mtry rule, sampling, ",
      "seed, and one-thread execution fixed before resampling"
    ),
    paste0(
      "Aalen additive-hazards time-varying coefficients and seeded native ",
      "tie handling fixed before resampling"
    ),
    paste0(
      "semiparametric proportional-odds iterations, profile fitting, ",
      "baseline variance, and seeded native tie handling fixed before ",
      "resampling"
    ),
    vapply(
      sb_ncvreg_survival_learner_ids(),
      sb_ncvreg_parameter_semantics,
      character(1)
    ),
    paste0(
      "ABESS sequence support path 0 through min(10,p), GIC criterion, ",
      "normalization, exact Newton updates, seed, and one-thread execution ",
      "fixed before resampling"
    ),
    sb_coxphf_parameter_semantics(sb_coxphf_learner_ids()),
    sb_bagged_knn_survival_parameter_semantics(
      sb_bagged_knn_survival_learner_ids()
    ),
    vapply(
      sb_ahaz_survival_learner_ids(),
      sb_ahaz_parameter_semantics,
      character(1)
    ),
    sb_bujar_linear_parameter_semantics(),
    "native Akritas conditional product-limit; analysis empirical joint CDF; fixed lambda 0.5; negative step-RMST through maximum analysis follow-up; no search; requires sensitive training predictor basis and paired outcomes"
  )
  survival_selected <- c(
    "fold-local Top-K then fixed ordinary Cox fit",
    paste0(
      "Top-K refit in each inner analysis fold; alpha and ",
      "standardization setting, lambda rule, and tuning measure fixed; ",
      "lambda selected by inner CV"
    ),
    paste0("fold-local Top-K; ", survival_base[-c(1L, 2L)])
  )
  survival <- paired_rows(
    "survival",
    survival_learners,
    "univ_cox_top_k",
    c(
      NA, "glmnet", "catboost", "catboost", "gbm", "mboost", "mboost",
      "mboost + partykit",
      vapply(
        sb_mboost_aft_learner_ids(),
        function(learner) {
          paste(sb_mboost_aft_dependency(learner)$packages, collapse = " + ")
        },
        character(1)
      ),
      vapply(
        sb_mboost_ranking_learner_ids(),
        function(learner) {
          paste(
            sb_mboost_ranking_dependency(learner)$packages,
            collapse = " + "
          )
        },
        character(1)
      ),
      rep("ranger", length(sb_ranger_survival_learner_ids())),
      rep("RLT", length(sb_rlt_survival_learner_ids())),
      "CoxBoost",
      "plsRcox", "superpc",
      rep("survivalsvm", length(sb_survivalsvm_learner_ids())),
      rep("xgboost", 6L),
      rep("survival", 4L),
      rep("flexsurv", length(sb_flexsurv_parametric_learner_ids())),
      rep("flexsurv", length(sb_flexsurvspline_learner_ids())),
      "MTLR",
      "pamr",
      "rpart", "randomForestSRC", "partykit + coin", "partykit + coin",
      "BART", "aorsf", rep("timereg", 2L),
      rep("ncvreg", length(sb_ncvreg_survival_learner_ids())),
      "abess", "coxphf", "bnnSurvival",
      rep("ahaz", length(sb_ahaz_survival_learner_ids())), "bujar", "survivalmodels + distr6"
    ),
    c(
      "none",
      "lambda_selected_by_inner_cv",
      rep("none", length(survival_learners) - 2L)
    ),
    survival_base,
    survival_selected
  )
  survival$tuning_semantics[
    survival$learner == "abess_cox_gic"
  ] <- "analysis_fold_gic"
  survival$tuning_semantics[
    survival$learner %in% sb_ncvreg_survival_learner_ids()
  ] <- "analysis_fold_bic"
  survival$tuning_semantics[
    survival$learner %in% c("ahaz_lasso_bic", "ahaz_sscad_bic")
  ] <- "analysis_fold_bic"
  survival <- rbind(
    survival,
    paired_rows(
      "survival",
      "bart_aft_survival",
      "univ_cox_top_k",
      "BART",
      "none",
      sb_rule_aft_parameter_semantics("bart_aft_survival"),
      paste0(
        "fold-local Top-K; ",
        sb_rule_aft_parameter_semantics("bart_aft_survival")
      )
    )
  )
  tuned_survival_map <- c(
    coxboost = "coxboost_cv",
    plsrcox = "plsrcox_cv",
    superpc = "superpc_cv",
    gbm_cox = "gbm_cox_cv"
  )
  tuned_survival <- do.call(rbind, lapply(
    names(tuned_survival_map),
    function(source_learner) {
      rows <- survival[survival$learner == source_learner, , drop = FALSE]
      rows$learner <- unname(tuned_survival_map[[source_learner]])
      rows$tuning_semantics <- "analysis_fold_internal_cv"
      rows$parameter_semantics <- ifelse(
        rows$selector == "none",
        "reference tuning repeated inside each analysis fold",
        paste0(
          "fold-local Top-K followed by reference tuning inside the same ",
          "analysis fold"
        )
      )
      rows
    }
  ))
  rownames(tuned_survival) <- NULL
  survival <- rbind(survival, tuned_survival)

  binary_learners <- c(
    "logistic", "glmnet_logistic", "gaussian_naive_bayes",
    sb_omics_specialist_binary_learner_ids(),
    sb_lightgbm_binary_learner_ids(),
    "catboost_binary", sb_ranger_binary_learner_ids(),
    sb_rlt_binary_learner_ids(),
    "adaboost_binary", "logitboost_binary",
    sb_classical_binary_learner_ids(),
    sb_specialized_binary_linear_learner_ids(),
    sb_xgboost_binary_learner_ids(), sb_ncvreg_binary_learner_ids(),
    "abess_binary_gic",
    sb_mclust_da_learner_ids(), sb_mixomics_da_learner_ids(),
    sb_hdda_binary_learner_ids(),
    sb_wsrf_binary_learner_ids(), sb_random_glm_binary_learner_ids()
  )
  binary_base <- c(
    "fixed ordinary logistic fit",
    paste0(
      "alpha, standardization setting, lambda rule, and tuning measure fixed; ",
      "lambda selected by inner CV"
    ),
    "variance smoothing fixed before resampling",
    vapply(
      sb_omics_specialist_binary_learner_ids(),
      sb_omics_specialist_binary_parameter_semantics,
      character(1)
    ),
    "LightGBM binary GBDT controls and seed fixed before resampling",
    paste0(
      "LightGBM binary DART tree and dropout controls and seed fixed before ",
      "resampling"
    ),
    paste0(
      "LightGBM binary GOSS top-gradient and random-small-gradient fractions ",
      "and seed fixed before resampling"
    ),
    paste0(
      "LightGBM binary random-forest row and feature fractions, frequency, ",
      "and seed fixed before resampling"
    ),
    "all learner parameters fixed before resampling",
    vapply(
      sb_ranger_binary_learner_ids(),
      sb_ranger_parameter_semantics,
      character(1)
    ),
    vapply(
      sb_rlt_binary_learner_ids(),
      sb_rlt_parameter_semantics,
      character(1)
    ),
    "AdaBoost.M1 maximum iteration count fixed before resampling",
    "LogitBoost iteration count fixed before resampling",
    "Euclidean nearest class mean, uniform priors, no scaling/shrinkage; uncalibrated squared-distance softmax",
    "Bernoulli NB with x > 0, Laplace alpha 1 and training-fold empirical priors",
    "empirical class priors fixed within each analysis fold",
    "empirical class priors fixed within each analysis fold",
    paste0(
      "kernel-density Naive Bayes priors, bandwidth adjustment, and Laplace ",
      "factor fixed before resampling"
    ),
    "RDA gamma and lambda covariance regularization fixed before resampling",
    "standardized Euclidean KNN with k=5",
    paste0(
      "single-hidden-layer neural-network size, decay, iteration limit, ",
      "scaling, and seed fixed before resampling"
    ),
    "classification tree controls fixed before resampling",
    "random-forest controls and seed fixed before resampling",
    "linear SVM cost and scaling fixed before resampling",
    "radial SVM cost, gamma rule, and scaling fixed before resampling",
    paste0(
      "radial SVM cost, gamma rule, scaling, and fold-local inverse-frequency ",
      "class weights fixed before resampling"
    ),
    "polynomial SVM cost, gamma rule, degree, and scaling fixed before resampling",
    "sigmoid SVM cost, gamma rule, offset, and scaling fixed before resampling",
    "GBM controls and seed fixed before resampling",
    paste0(
      "probit BART trees, posterior draws, burn-in, priors, and seed fixed ",
      "before resampling"
    ),
    "boosted C5.0 trials, pruning controls, and seed fixed before resampling",
    "kernel-weighted KNN k, distance, kernel, and scaling fixed before resampling",
    paste0(
      "Gaussian-process RBF kernel, automatic analysis-fold width rule, ",
      "scaling, tolerance, and seed fixed before resampling"
    ),
    paste0(
      "nearest shrunken-centroid threshold, variance scaling, offset, and ",
      "empirical analysis-fold priors fixed before resampling"
    ),
    paste0(
      "two-hidden-layer resilient-backprop architecture, stopping threshold, ",
      "scaling, and seed fixed before resampling"
    ),
    "conditional-inference classification-tree controls fixed before resampling",
    paste0(
      "conditional-inference classification-forest controls, mtry rule, ",
      "sampling, seed, and one-core execution fixed before resampling"
    ),
    "bootstrap Rpart bagging controls and seed fixed before resampling",
    paste0(
      "smooth additive Binomial boosting iterations, shrinkage, and spline ",
      "base learners fixed before resampling"
    ),
    paste0(
      "conditional-tree Binomial boosting iterations, shrinkage, and tree ",
      "controls fixed before resampling"
    ),
    paste0(
      "localized LDA neighbor cap and locally weighted empirical priors ",
      "fixed before resampling"
    ),
    paste0(
      "oblique classification-forest trees, split controls, mtry rule, ",
      "sampling, seed, and one-thread execution fixed before resampling"
    ),
    paste0(
      "Earth hinge-basis degree, basis cap, backward pruning rule, and ",
      "Binomial-logit family fixed before resampling"
    ),
    paste0(
      "flexible discriminant polynomial degree, basis convention, and one ",
      "discriminant dimension fixed before resampling"
    ),
    paste0(
      "mixture discriminant subclasses, iterations, seed, linear polyreg ",
      "map, and discriminant dimension fixed before resampling"
    ),
    vapply(
      sb_liblinear_binary_learner_ids(),
      sb_liblinear_binary_parameter_semantics,
      character(1)
    ),
    paste0(
      "evolutionary-tree population, iteration budget, depth controls, ",
      "complexity weight, and seed fixed before resampling"
    ),
    paste0(
      "shrinkage discriminant covariance, variance, and class-frequency ",
      "intensities estimated analytically inside each analysis fold"
    ),
    paste0(
      "rotation-forest tree count, variable-subset rule, tree controls, and ",
      "seed fixed before resampling"
    ),
    paste0(
      "regularized random-forest tree count, mtry rule, node size, ",
      "regularization coefficient, and seed fixed before resampling"
    ),
    paste0(
      "penalized discriminant ridge lambda 1, one discriminant dimension, ",
      "and empirical analysis-fold priors fixed before resampling"
    ),
    paste0(
      "stabilised LDA native PC_q score-dimension rule and empirical class ",
      "priors applied inside each analysis fold"
    ),
    paste0(
      "classic tree deviance split, minimum node sizes, and minimum deviance ",
      "reduction fixed before resampling"
    ),
    paste0(
      "C5.0 rule-set trials, pruning controls, and seed fixed before ",
      "resampling"
    ),
    paste0(
      "random-fern depth, ensemble size, seed, one-thread execution, and ",
      "fixed softmax mapping of native class scores"
    ),
    paste0(
      "robust Mqle logistic tuning constant, convergence controls, and ",
      "predictor-weighting rule fixed before resampling"
    ),
    paste0(
      "full radial least-squares SVM tau and scaling fixed; predictor-only ",
      "kernel width estimated by native sigest inside each analysis fold ",
      "with local seed 1; fixed softmax mapping of native decision scores"
    ),
    paste0(
      "LVQ1 equal-prior-target five-neighbor native initialization with the ",
      "native size rule and minimum two, analysis-fold scaling, 100 updates ",
      "per prototype, alpha 0.03, seed 1, and fixed nearest-prototype score mapping"
    ),
    paste0(
      "optimized LVQ1 equal-prior-target five-neighbor native initialization ",
      "with the native size rule and minimum two, analysis-fold scaling, 40 ",
      "updates per prototype, alpha 0.3, seed 1, and fixed nearest-prototype ",
      "score mapping"
    ),
    paste0(
      "LVQ2.1 equal-prior-target five-neighbor native initialization with the ",
      "native size rule and minimum two, analysis-fold scaling, 100 updates per ",
      "prototype, alpha 0.03, window 0.3, seed 1, and fixed nearest-prototype ",
      "score mapping"
    ),
    paste0(
      "LVQ3 equal-prior-target five-neighbor native initialization with the ",
      "native size rule and minimum two, analysis-fold scaling, 100 updates per ",
      "prototype, alpha 0.03, window 0.3, epsilon 0.1, seed 1, and fixed ",
      "nearest-prototype score mapping"
    ),
    "direction and AIC stopping rule fixed before resampling",
    paste0(
      "iteration count, shrinkage, family, and centering fixed before ",
      "resampling"
    ),
    "component cap, scaling, and logistic PLS model fixed before resampling",
    "XGBoost binary gbtree controls and seed fixed before resampling",
    paste0(
      "XGBoost binary DART tree and dropout controls and seed fixed before ",
      "resampling"
    ),
    paste0(
      "XGBoost binary coordinate-descent linear booster penalty, update ",
      "order, and seed fixed before resampling"
    ),
    vapply(
      sb_ncvreg_binary_learner_ids(),
      sb_ncvreg_parameter_semantics,
      character(1)
    ),
    paste0(
      "ABESS sequence support path 0 through min(10,p), GIC criterion, ",
      "normalization, exact Newton updates, seed, and one-thread execution ",
      "fixed before resampling"
    ),
    paste0(
      "one Gaussian component per class with a common covariance structure ",
      "selected by BIC inside each analysis fold; known class labels require ",
      "no random initialization"
    ),
    paste0(
      "one through three Gaussian components and covariance structure ",
      "selected separately by class using BIC inside each analysis fold; ",
      "hierarchical subset cap 2000, VVV/SVD initialization, and seed 1 fixed"
    ),
    paste0(
      "two-component cap, analysis-fold scaling, and fixed softmax mapping ",
      "of native PLS-DA dummy-response scores"
    ),
    paste0(
      "two-component cap, up to ten nonzero loadings per component, ",
      "analysis-fold scaling, and fixed softmax mapping of native sparse ",
      "PLS-DA dummy-response scores"
    ),
    sb_hdda_binary_parameter_semantics(sb_hdda_binary_learner_ids()),
    sb_wsrf_binary_parameter_semantics(sb_wsrf_binary_learner_ids()),
    sb_slope_randomglm_binary_parameter_semantics(
      sb_random_glm_binary_learner_ids()
    )
  )
  binary_selected <- c(
    "fold-local Top-K then fixed ordinary logistic fit",
    paste0(
      "Top-K refit in each inner analysis fold; alpha and ",
      "standardization setting, lambda rule, and tuning measure fixed; ",
      "lambda selected by inner CV"
    ),
    paste0("fold-local Top-K; ", binary_base[-c(1L, 2L)])
  )
  binary <- paired_rows(
    "binary",
    binary_learners,
    "univ_binary_auc_top_k",
    c(
      NA, "glmnet", NA,
      vapply(
        sb_omics_specialist_binary_learner_ids(),
        function(learner) {
          sb_omics_specialist_binary_dependency(learner)$package
        },
        character(1)
      ),
      rep("lightgbm", 4L), "catboost",
      rep("ranger", length(sb_ranger_binary_learner_ids())),
      rep("RLT", length(sb_rlt_binary_learner_ids())),
      "fastAdaboost", "caTools", NA, NA, "MASS", "MASS", "klaR", "klaR", NA,
      "nnet", "rpart", "randomForest", "e1071", "e1071", "e1071",
      "e1071", "e1071", "gbm", "BART", "C50", "kknn", "kernlab",
      "pamr", "neuralnet",
      "partykit", "partykit", "ipred", "mboost", "mboost", "klaR",
      "aorsf", "earth", "mda", "mda",
      rep("LiblineaR", length(sb_liblinear_binary_learner_ids())),
      "evtree", "sda",
      "rotationForest", "RRF", "mda", "ipred", "tree", "C50",
      "rFerns", "robustbase", "kernlab", rep("class", 4L),
      NA, "mboost", "plsRglm",
      rep("xgboost", 3L),
      rep("ncvreg", length(sb_ncvreg_binary_learner_ids())),
      "abess",
      "mclust", "mclust",
      "mixOmics", "mixOmics", "HDclassif", "wsrf", "randomGLM"
    ),
    c(
      "none",
      "lambda_selected_by_inner_cv",
      rep("none", length(binary_learners) - 2L)
    ),
    binary_base,
    binary_selected
  )
  binary$tuning_semantics[
    binary$learner == "abess_binary_gic"
  ] <- "analysis_fold_gic"
  binary$tuning_semantics[
    binary$learner %in% sb_mclust_da_learner_ids()
  ] <- "analysis_fold_bic"
  binary$tuning_semantics[
    binary$learner %in% sb_ncvreg_binary_learner_ids()
  ] <- "analysis_fold_bic"
  binary$tuning_semantics[
    binary$learner == "lssvm_radial_binary"
  ] <- "analysis_fold_predictor_only_sigest"
  binary$tuning_semantics[
    binary$learner %in% sb_svm_probability_learner_ids()
  ] <- "analysis_fold_oof_probability_calibration"
  registry <- rbind(survival, binary)
  stepcox <- registry[
    registry$task_type == "survival" & registry$selector == "none", ,
    drop = FALSE
  ]
  stepcox$selector <- "stepcox"
  stepcox$selector_parameter <- "stepcox_direction"
  stepcox$parameter_semantics <- ifelse(
    stepcox$tuning_semantics == "lambda_selected_by_inner_cv",
    paste0(
      "StepCox AIC selection refit in each inner analysis fold; ",
      "direction, alpha, standardization setting, lambda rule, and ",
      "tuning measure fixed; lambda selected by inner CV"
    ),
    paste0(
      "fold-local StepCox AIC selection with direction fixed before ",
      "resampling; ",
      stepcox$parameter_semantics
    )
  )
  registry <- rbind(registry, stepcox)

  combine_dependencies <- function(selector_dependency, learner_dependency) {
    if (is.na(selector_dependency) || !nzchar(selector_dependency)) {
      return(learner_dependency)
    }
    if (is.na(learner_dependency) || !nzchar(learner_dependency)) {
      return(selector_dependency)
    }
    dependencies <- unique(c(
      strsplit(selector_dependency, " + ", fixed = TRUE)[[1L]],
      strsplit(learner_dependency, " + ", fixed = TRUE)[[1L]]
    ))
    paste(dependencies, collapse = " + ")
  }
  abess_selector_rows <- function(task_type, selector, excluded_learners) {
    rows <- registry[
      registry$task_type == task_type & registry$selector == "none" &
        !registry$learner %in% excluded_learners, ,
      drop = FALSE
    ]
    rows$source <- "stormbreaker_abess_runtime"
    rows$selector <- selector
    rows$selector_parameter <- NA_character_
    rows$optional_dependency <- mapply(
      combine_dependencies,
      "abess",
      rows$optional_dependency,
      USE.NAMES = FALSE
    )
    rows$parameter_semantics <- paste0(
      "fold-local ABESS sequence path with GIC-selected nonzero support ",
      "from 2 through min(10,p); ",
      rows$parameter_semantics
    )
    rows
  }
  registry <- rbind(
    registry,
    abess_selector_rows(
      "binary",
      "abess_binary_gic_nonzero",
      "abess_binary_gic"
    ),
    abess_selector_rows(
      "survival",
      "abess_cox_gic_nonzero",
      c("abess_cox_gic", sb_analysis_tuned_survival_learner_ids())
    )
  )
  learned_selector_rows <- function(task_type, selector, selector_dependency,
                                    excluded_learners, selection_semantics,
                                    source) {
    rows <- registry[
      registry$task_type == task_type & registry$selector == "none" &
        !registry$learner %in% excluded_learners, ,
      drop = FALSE
    ]
    rows$source <- source
    rows$selector <- selector
    rows$selector_parameter <- NA_character_
    rows$optional_dependency <- mapply(
      combine_dependencies,
      selector_dependency,
      rows$optional_dependency,
      USE.NAMES = FALSE
    )
    rows$parameter_semantics <- paste0(
      selection_semantics,
      "; ",
      rows$parameter_semantics
    )
    rows
  }
  information_criteria <- c(
    jmi_binary_nlogn = "joint mutual-information (JMI)",
    cmim_binary_nlogn = "conditional mutual-information maximization (CMIM)",
    disr_binary_nlogn = "double-input symmetrical relevance (DISR)",
    jmim_binary_nlogn = "minimal joint mutual-information maximization (JMIM)"
  )
  model_selectors <- sb_binary_model_selector_ids()
  for (id in sb_pamr_selector_ids()) registry <- rbind(registry,
    learned_selector_rows(if (id == "pamr_survival_nonzero") "survival" else "binary",
      id, "pamr", sb_analysis_tuned_survival_learner_ids(),
      paste0("fold-local native PAMR nonzero shrunken-centroid features; fixed threshold 1; ",
        if (id == "pamr_survival_nonzero") "two native censoring-aware survival groups" else "empirical binary class priors",
        "; scale.sd TRUE, offset 50, no threshold search or feature padding"),
      "stormbreaker_pamr_selector_runtime"))
  for (i in seq_along(sb_survival_weight_selector_ids())) registry <- rbind(registry,
    learned_selector_rows("survival", sb_survival_weight_selector_ids()[i],
      c("glmnet", "plsRcox", "survivalsvm")[i], sb_analysis_tuned_survival_learner_ids(),
      paste0("fold-local ", c("Ridge Cox", "PLS-Cox", "linear survival-SVM")[i],
        " absolute coefficient times training SD; retain top half, minimum 2 positive weights when available"),
      "stormbreaker_survival_weight_selector_runtime"))
  for (i in seq_along(sb_survival_boosting_selector_ids())) registry <- rbind(registry,
    learned_selector_rows("survival", sb_survival_boosting_selector_ids()[i],
      c("gbm", "mboost", "superpc", "mboost", "xgboost")[i], sb_analysis_tuned_survival_learner_ids(),
      c("fold-local Cox GBM positive relative influence; 100 trees, depth 1, shrinkage 0.1",
        "fold-local linear Cox glmBoost nonzero coefficients; 100 steps, shrinkage 0.1",
        "fold-local SuperPC native survival scores; absolute threshold 1, s0 percentile 0.5",
        "fold-local Cox GamBoost active bbs spline terms; dfbase 4, 100 steps, shrinkage 0.1; score is spline coefficient L2 norm",
        "fold-local native XGBoost Cox positive gain; 100 rounds, eta 0.1, depth 3, full sampling, one thread, seed 1")[i],
      "stormbreaker_survival_boosting_selector_runtime"))
  for (i in seq_along(model_selectors)) registry <- rbind(registry,
    learned_selector_rows("binary", model_selectors[i],
      c("glmnet", NA, "mboost", "randomForest", "glmnet", "gbm", "e1071", "mboost", "xgboost", "plsRglm")[i], character(),
      c("fold-local binomial LASSO CV nonzero coefficients",
        "fold-local stepwise logistic AIC variables",
        "fold-local glmBoost nonzero coefficients",
        "fold-local RF importance strictly above declared threshold",
        "fold-local binomial elastic-net CV nonzero coefficients",
        "fold-local GBM positive relative influence",
        "linear SVM recursive elimination with grouped inner-CV feature count selection",
        "fold-local binomial GamBoost active bbs spline terms; dfbase 4, 100 steps, shrinkage 0.1; score is spline coefficient L2 norm",
        "fold-local native XGBoost binary positive gain; 100 rounds, eta 0.1, depth 3, full sampling, one thread, seed 1",
        "fold-local native logistic PLS; at most two scaled components; absolute coefficient times training SD; retain top half, minimum two positive weights when available; no search")[i],
      "stormbreaker_binary_model_selector_runtime"))
  information_selector_rows <- do.call(rbind, lapply(
    names(information_criteria),
    function(selector) {
      learned_selector_rows(
        "binary",
        selector,
        "praznik",
        character(),
        paste0(
          "fold-local ", information_criteria[[selector]],
          " selection using the native praznik ranking, ",
          "min(p, max(2, floor(n/log(n)))) retained predictors, all ",
          "finite native scores, and one thread"
        ),
        "stormbreaker_high_dimensional_runtime"
      )
    }
  ))
  registry <- rbind(
    registry,
    learned_selector_rows(
      "binary",
      "boruta_binary_roughfix",
      "Boruta",
      character(),
      paste0(
        "fold-local Boruta all-relevant selection using ranger normalized ",
        "permutation importance, multiplicity-adjusted p=0.01 decisions, ",
        "at most 100 runs, 500 trees, one thread, seed 1, and ",
        "TentativeRoughFix"
      ),
      "stormbreaker_forest_selector_runtime"
    ),
    learned_selector_rows(
      "binary",
      "varselrf_binary_oob_1se",
      "varSelRF",
      character(),
      paste0(
        "fold-local varSelRF backward elimination using fixed initial ",
        "MeanDecreaseAccuracy, the one-standard-error OOB subset rule, ",
        "5,000 initial trees, 2,000 iteration trees, 20% variable drops, ",
        "whole-range search, and seed 1"
      ),
      "stormbreaker_forest_selector_runtime"
    ),
    learned_selector_rows(
      "survival",
      "boruta_survival_roughfix",
      "Boruta",
      sb_analysis_tuned_survival_learner_ids(),
      paste0(
        "fold-local survival Boruta all-relevant selection using ranger ",
        "normalized permutation importance, multiplicity-adjusted p=0.01 ",
        "decisions, at most 100 runs, 500 trees, one thread, seed 1, and ",
        "TentativeRoughFix"
      ),
      "stormbreaker_forest_selector_runtime"
    ),
    learned_selector_rows(
      "binary",
      "sis_binary_nlogn",
      "SIS",
      character(),
      paste0(
        "fold-local binomial SIS marginal-utility ranking with ",
        "min(p, max(2, floor(n/log(n)))) retained predictors, native ",
        "standardization, non-iterative screening, and seed 1"
      ),
      "stormbreaker_high_dimensional_runtime"
    ),
    learned_selector_rows(
      "binary",
      "vsurf_binary_prediction",
      "VSURF",
      character(),
      paste0(
        "fold-local VSURF threshold, interpretation, and prediction stages ",
        "using its randomForest backend, native default forest counts, ",
        "serial execution, and seed 1"
      ),
      "stormbreaker_high_dimensional_runtime"
    ),
    learned_selector_rows(
      "binary",
      "mrmr_binary_nlogn",
      "praznik",
      character(),
      paste0(
        "fold-local minimum-redundancy maximum-relevance selection using ",
        "the native praznik ranking, min(p, max(2, floor(n/log(n)))) ",
        "retained predictors, all finite native scores, and one thread"
      ),
      "stormbreaker_high_dimensional_runtime"
    ),
    information_selector_rows,
    learned_selector_rows(
      "binary",
      "lpc_binary_nlogn",
      "lpc",
      character(),
      paste0(
        "fold-local Lassoed Principal Components using native two-class ",
        "t scores, the native 10 repeated half-split predictive-advantage ",
        "soft threshold, decreasing absolute LPC-score ranking, ",
        "min(p, max(2, floor(n/log(n)))) retained predictors, and seed 1"
      ),
      "stormbreaker_high_dimensional_runtime"
    ),
    learned_selector_rows(
      "survival",
      "lpc_cox_nlogn",
      "lpc",
      sb_analysis_tuned_survival_learner_ids(),
      paste0(
        "fold-local Lassoed Principal Components using native survival Cox ",
        "scores, the native 10 repeated half-split predictive-advantage ",
        "soft threshold, decreasing absolute LPC-score ranking, ",
        "min(p, max(2, floor(n/log(n)))) retained predictors, and seed 1"
      ),
      "stormbreaker_high_dimensional_runtime"
    ),
    learned_selector_rows(
      "binary",
      "stabsel_binary_lasso",
      "stabs + glmnet",
      "glmnet_logistic",
      paste0(
        "fold-local complementary-pairs stability selection with ",
        "conservative standardized binomial Lasso paths, 50 pairs, ",
        "min(max(2, floor(n/log(n))), floor(p/2)) variables per ",
        "subsample with at least four candidate predictors, ",
        "selection-probability cutoff 0.75, unimodality assumption, ",
        "serial execution, and seed 1"
      ),
      "stormbreaker_high_dimensional_runtime"
    ),
    learned_selector_rows(
      "survival",
      "stabsel_cox_lasso",
      "stabs + glmnet",
      c(sb_analysis_tuned_survival_learner_ids(), "glmnet_cox"),
      paste0(
        "fold-local complementary-pairs stability selection with ",
        "conservative standardized Cox Lasso paths with Breslow ties, ",
        "50 pairs, ",
        "min(max(2, floor(n/log(n))), floor(p/2)) variables per ",
        "subsample with at least four candidate predictors, ",
        "selection-probability cutoff 0.75, unimodality assumption, ",
        "serial execution, and seed 1"
      ),
      "stormbreaker_high_dimensional_runtime"
    )
  )
  mime_selector_rows <- function(selector, learners, selector_parameter,
                                 selector_dependency, selection_semantics) {
    rows <- registry[
      registry$task_type == "survival" &
        registry$selector == "none" &
        registry$learner %in% learners, ,
      drop = FALSE
    ]
    rows <- rows[match(learners, rows$learner), , drop = FALSE]
    rows$source <- "stormbreaker_mime117_runtime"
    rows$selector <- selector
    rows$selector_parameter <- selector_parameter
    rows$optional_dependency <- mapply(
      combine_dependencies,
      selector_dependency,
      rows$optional_dependency,
      USE.NAMES = FALSE
    )
    rows$parameter_semantics <- paste0(
      selection_semantics,
      "; ",
      rows$parameter_semantics
    )
    rows
  }

  registry <- rbind(
    registry,
    mime_selector_rows(
      "rsf_minimal_depth",
      c(
        "coxboost_cv", "glmnet_cox", "gbm_cox_cv",
        "plsrcox_cv", "superpc_cv", sb_survivalsvm_learner_ids()
      ),
      "selector_rsf_nodesize_and_seed",
      "randomForestSRC",
      paste0(
        "fold-local 1,000-tree RSF minimal-depth selection with ",
        "conservative='high'"
      )
    ),
    mime_selector_rows(
      "coxboost_nonzero",
      c(
        "glmnet_cox", "gbm_cox_cv", "plsrcox_cv",
        "superpc_cv", sb_survivalsvm_learner_ids()
      ),
      "selector_inner_folds_and_seed",
      "CoxBoost",
      paste0(
        "fold-local CoxBoost penalty and step tuning followed by nonzero ",
        "coefficient selection"
      )
    ),
    mime_selector_rows(
      "enet_cox_nonzero",
      unique(c(registry$learner[registry$selector == "lasso_nonzero"],
        "coxph", "glmnet_cox", "ranger_survival",
        "coxboost_cv", "gbm_cox_cv", "plsrcox_cv",
        "rsf_survival", "superpc_cv", sb_survivalsvm_learner_ids())),
      "selector_inner_folds_and_seed", "glmnet",
      "fold-local Cox Enet alpha=0.5 lambda.min nonzero coefficients"
    ),
    mime_selector_rows(
      "lasso_nonzero",
      c(
        "coxboost_cv", "gbm_cox_cv", "plsrcox_cv",
        "rsf_survival", "superpc_cv", sb_survivalsvm_learner_ids()
      ),
      "selector_inner_folds_and_seed",
      "glmnet",
      paste0(
        "fold-local Cox Lasso lambda.min tuning followed by nonzero ",
        "coefficient selection"
      )
    ),
    mime_selector_rows(
      "rsf_then_stepcox",
      "coxph",
      "selector_rsf_nodesize_seed_and_stepcox_direction",
      "randomForestSRC",
      paste0(
        "fold-local 1,000-tree RSF minimal-depth selection followed by ",
        "fold-local AIC StepCox"
      )
    ),
    mime_selector_rows(
      "coxboost_then_stepcox",
      "coxph",
      "selector_inner_folds_seed_and_stepcox_direction",
      "CoxBoost",
      paste0(
        "fold-local tuned CoxBoost nonzero selection followed by fold-local ",
        "AIC StepCox"
      )
    ),
    mime_selector_rows(
      "lasso_then_stepcox",
      "coxph",
      "selector_inner_folds_seed_and_stepcox_direction",
      "glmnet",
      paste0(
        "fold-local tuned Cox Lasso nonzero selection followed by fold-local ",
        "AIC StepCox"
      )
    )
  )
  for (selector in c("lasso_nonzero", "enet_cox_nonzero", "coxboost_nonzero", "rsf_minimal_depth")) {
    dependency <- switch(selector, lasso_nonzero = "glmnet", enet_cox_nonzero = "glmnet", coxboost_nonzero = "CoxBoost",
      rsf_minimal_depth = "randomForestSRC")
    extra <- learned_selector_rows("survival", selector, dependency,
      sb_analysis_tuned_survival_learner_ids(),
      "fold-local native feature selection; existing selector defaults retained",
      "stormbreaker_mime117_runtime")
    extra$selector_parameter <- if (selector == "rsf_minimal_depth")
      "selector_rsf_nodesize_and_seed" else "selector_inner_folds_and_seed"
    extra <- extra[!extra$learner %in% registry$learner[registry$selector == selector], , drop = FALSE]
    registry <- rbind(registry, extra)
  }
  restricted_internal_cv <- data.frame(
    source = rep("stormbreaker_builtin", 4L),
    task_type = c("binary", "survival", "survival", "binary"),
    selector = rep("none", 4L),
    learner = c(
      "pre_binary_rule_ensemble",
      "pre_survival_rule_ensemble",
      "penaft_gehan_elastic_net",
      "cancerclass_binary"
    ),
    admission_status = rep("admitted", 4L),
    selector_parameter = rep(NA_character_, 4L),
    tuning_semantics = c(
      "analysis_fold_native_coefficient_cv_after_rule_induction",
      "analysis_fold_native_coefficient_cv_after_rule_induction",
      "analysis_fold_internal_cv_fixed",
      "analysis_fold_oof_probability_calibration"
    ),
    optional_dependency = c(
      "pre", "pre + mboost", "penAFT",
      "cancerclass + Biobase"
    ),
    parameter_semantics = c(
      vapply(
        c(
          "pre_binary_rule_ensemble",
          "pre_survival_rule_ensemble",
          "penaft_gehan_elastic_net"
        ),
        sb_rule_aft_parameter_semantics,
        character(1)
      ),
      paste0(
        "native Welch ranking and correlation-distance z score over all ",
        "supplied features; five-fold OOF logistic calibration inside each ",
        "analysis fold"
      )
    ),
    stringsAsFactors = FALSE
  )
  registry <- rbind(registry, restricted_internal_cv)
  first_learner_row <- !duplicated(registry$learner)
  fixed <- registry[first_learner_row, , drop = FALSE]
  fixed$selector <- "fixed_features"
  fixed$selector_parameter <- "selector_features_and_feature_set_id"
  fixed$parameter_semantics <- paste0(
    "predeclared outcome-independent feature set; ",
    fixed$parameter_semantics
  )
  rownames(fixed) <- NULL
  rbind(registry, fixed)
}

sb_pipeline_registry_match <- function(selector, learner,
                                       registry = sb_pipeline_registry_data()) {
  registry[
    registry$selector == selector & registry$learner == learner, ,
    drop = FALSE
  ]
}

sb_pipeline_uses_inner_cv <- function(learner,
                                      registry = sb_pipeline_registry_data()) {
  semantics <- unique(registry$tuning_semantics[
    registry$learner == learner
  ])
  identical(semantics, "lambda_selected_by_inner_cv")
}

sb_pipeline_registry_learner_view <- function(registry) {
  learners <- registry[registry$selector == "none", , drop = FALSE]
  route_counts <- table(registry$learner)
  selectors <- vapply(
    learners$learner,
    function(learner) {
      paste(registry$selector[registry$learner == learner], collapse = ", ")
    },
    character(1)
  )
  data.frame(
    task_type = learners$task_type,
    learner = learners$learner,
    pipeline_count = unname(as.integer(route_counts[learners$learner])),
    selectors = selectors,
    tuning_semantics = learners$tuning_semantics,
    optional_dependency = learners$optional_dependency,
    parameter_semantics = learners$parameter_semantics,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' @rdname sb_models
#' @export
sb_pipeline_registry <- function(task_type = NULL,
                                 level = c("route", "learner"),
                                 check_installed = FALSE) {
  level <- match.arg(level)
  if (!is.null(task_type)) {
    task_type <- match.arg(task_type, c("survival", "binary"))
  }
  if (!is.logical(check_installed) || length(check_installed) != 1L ||
    is.na(check_installed)) {
    sb_abort(
      "`check_installed` must be `TRUE` or `FALSE`.",
      "stormbreaker_bad_registry_check_installed"
    )
  }
  registry <- sb_pipeline_registry_data()
  if (!is.null(task_type)) {
    registry <- registry[registry$task_type == task_type, , drop = FALSE]
    rownames(registry) <- NULL
  }
  if (identical(level, "learner")) {
    registry <- sb_pipeline_registry_learner_view(registry)
  }
  if (check_installed) {
    registry <- sb_add_dependency_status(registry)
  }
  registry
}

sb_pipeline_lint_report <- function(status, pipeline = NULL, spec = NULL,
                                    task_type = NULL, error = NULL) {
  valid_named_spec <- !inherits(spec, "stormbreaker_pipeline") &&
    is.list(spec) && (length(spec) == 0L ||
    (!is.null(names(spec)) && !any(is.na(names(spec))) &&
      !any(!nzchar(names(spec))) && anyDuplicated(names(spec)) == 0L))
  raw_value <- function(name, default = NA_character_) {
    if (!valid_named_spec) {
      return(NA_character_)
    }
    if (!name %in% names(spec)) {
      return(default)
    }
    value <- spec[[name]]
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
      return(NA_character_)
    }
    value
  }
  pipeline_value <- function(name) {
    if (!is.list(pipeline) || !name %in% names(pipeline)) {
      return(NA_character_)
    }
    value <- pipeline[[name]]
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
      return(NA_character_)
    }
    value
  }
  has_pipeline <- is.list(pipeline)
  pipeline_id <- pipeline_value("id")
  reported_task <- if (!is.null(task_type)) {
    task_type
  } else {
    pipeline_value("task_type")
  }
  selector <- if (has_pipeline) {
    pipeline_value("selector")
  } else {
    raw_value("selector", "none")
  }
  learner <- if (has_pipeline) {
    pipeline_value("learner")
  } else {
    raw_value("learner", "coxph")
  }
  issue_code <- if (is.null(error)) NA_character_ else class(error)[[1L]]
  message <- if (is.null(error)) NA_character_ else conditionMessage(error)

  data.frame(
    status = status,
    pipeline_id = pipeline_id,
    task_type = reported_task,
    selector = selector,
    learner = learner,
    issue_code = issue_code,
    message = message,
    stringsAsFactors = FALSE
  )
}

sb_is_expected_pipeline_lint_error <- function(error) {
  sb_is_expected_pipeline_constructor_error(error) ||
    inherits(error, c(
      "stormbreaker_bad_pipeline_spec",
      "stormbreaker_unknown_pipeline_argument",
      "stormbreaker_bad_pipeline_task_type",
      "stormbreaker_pipeline_task_mismatch",
      "stormbreaker_unresolved_pipeline"
    ))
}

sb_pipeline_from_lint_spec <- function(spec) {
  if (!is.list(spec)) {
    sb_abort(
      paste0(
        "`spec` must be one `stormbreaker_pipeline` or one named list ",
        "of arguments for `sb_pipeline()`."
      ),
      "stormbreaker_bad_pipeline_spec"
    )
  }
  if (length(spec) > 0L && (is.null(names(spec)) ||
    any(is.na(names(spec))) || any(!nzchar(names(spec))) ||
    anyDuplicated(names(spec)) > 0L)) {
    sb_abort(
      "A pipeline argument list must have unique, non-empty names.",
      "stormbreaker_bad_pipeline_spec"
    )
  }
  unknown <- setdiff(names(spec), names(formals(sb_pipeline)))
  if (length(unknown) > 0L) {
    sb_abort(
      paste0(
        "Unknown `sb_pipeline()` argument", if (length(unknown) > 1L) {
          "s"
        } else {
          ""
        }, ": ", paste0("`", unknown, "`", collapse = ", "), "."
      ),
      "stormbreaker_unknown_pipeline_argument"
    )
  }
  do.call(sb_pipeline, spec)
}

#' @rdname sb_pipeline
#' @export
sb_lint_pipeline <- function(spec, task_type = NULL) {
  if (!is.null(task_type) && (!is.character(task_type) ||
    length(task_type) != 1L || is.na(task_type) ||
    !task_type %in% c("survival", "binary"))) {
    sb_abort(
      "`task_type` must be `NULL`, `\"survival\"`, or `\"binary\"`.",
      "stormbreaker_bad_pipeline_lint_task"
    )
  }

  pipeline <- NULL
  result <- tryCatch(
    {
      pipeline <- if (inherits(spec, "stormbreaker_pipeline")) {
        spec
      } else {
        sb_pipeline_from_lint_spec(spec)
      }
      target_task <- if (is.null(task_type)) {
        pipeline[["task_type"]]
      } else {
        task_type
      }
      if (!is.character(target_task) || length(target_task) != 1L ||
        is.na(target_task) ||
        !target_task %in% c("survival", "binary")) {
        sb_abort(
          "The pipeline does not declare a supported task type.",
          "stormbreaker_bad_pipeline_task_type"
        )
      }
      pipeline <- sb_validate_pipeline_contract(pipeline, target_task)
      pipeline
    },
    stormbreaker_error = function(error) {
      if (!sb_is_expected_pipeline_lint_error(error)) {
        stop(error)
      }
      error
    }
  )

  if (inherits(result, "stormbreaker_error")) {
    return(sb_pipeline_lint_report(
      status = "rejected",
      pipeline = pipeline,
      spec = spec,
      task_type = task_type,
      error = result
    ))
  }
  runtime <- sb_resolve_pipeline_runtime(result)
  sb_validate_pipeline_runtime(runtime, result)
  sb_pipeline_lint_report(
    status = "admitted",
    pipeline = result,
    spec = spec,
    task_type = task_type
  )
}

# ---- candidate-menus ----

# ---- binary-menu-pipelines ----

sb_binary_menu_arguments <- function(row, selector_features, feature_set_id,
                                     selector_n_features, inner_folds,
                                     inner_seed) {
  args <- list(selector = row$selector, learner = row$learner)
  if (identical(row$selector_parameter, "selector_n_features")) {
    args$selector_n_features <- selector_n_features
  }
  if (identical(row$selector_parameter, "selector_features_and_feature_set_id")) {
    args$selector_features <- selector_features
    args$feature_set_id <- feature_set_id
  }
  if (identical(row$tuning_semantics, "lambda_selected_by_inner_cv")) {
    args$inner_folds <- inner_folds
    args$inner_seed <- inner_seed
  }
  if (identical(row$learner, "stepglm_binary")) {
    args$stepglm_direction <- "both"
  }
  args
}

#' @rdname sb_models
#' @export
sb_binary_workflow_plan <- function(models = NULL, selector = "rf_binary_importance",
                                    selector_options = NULL, selector_n_features = 5L,
                                    inner_folds = 5L, inner_seed = 1701L,
                                    include_unfiltered = FALSE, selector_features = NULL,
                                    feature_set_id = NULL) {
  registry <- sb_pipeline_registry("binary")
  registry <- registry[registry$admission_status == "admitted", , drop = FALSE]
  if (!is.character(selector) || length(selector) != 1L || is.na(selector)) {
    stop("Choose one binary selector ID or 'none'.", call. = FALSE)
  }
  aliases <- sb_doctor_selector_menu("binary")
  alias <- match(tolower(selector), tolower(names(aliases)))
  if (!is.na(alias)) selector <- unname(aliases[alias])
  if (!selector %in% registry$selector) {
    stop("Unknown binary selector: ", selector, ". See sb_pipeline_registry('binary').", call. = FALSE)
  }
  if (!is.logical(include_unfiltered) || length(include_unfiltered) != 1L || is.na(include_unfiltered)) {
    stop("include_unfiltered must be TRUE or FALSE.", call. = FALSE)
  }
  if (selector == "fixed_features") {
    selector_features <- sb_validate_fixed_features(selector_features)
    feature_set_id <- sb_validate_feature_set_id(feature_set_id)
  } else if (!is.null(selector_features) || !is.null(feature_set_id)) {
    stop("selector_features and feature_set_id require selector = 'fixed_features'.", call. = FALSE)
  }
  selector_n_features <- sb_check_positive_integer(selector_n_features, "selector_n_features", minimum = 2L)
  inner_folds <- sb_check_positive_integer(inner_folds, "inner_folds", minimum = 3L)
  inner_seed <- sb_check_positive_integer(inner_seed, "inner_seed")
  if (selector == "rf_binary_importance" && is.null(selector_options)) {
    selector_options <- list(top_n = 5L, spearman_cutoff = .70)
  }

  menu <- sb_models("binary")
  if (!is.null(models)) {
    rows <- sb_binary_method_ids(models, menu)
    if (anyDuplicated(rows)) stop("Choose each binary method only once.", call. = FALSE)
    menu <- menu[rows, , drop = FALSE]
  }
  chosen <- registry[registry$selector == selector, , drop = FALSE]
  route <- match(menu$learner, chosen$learner)
  if (include_unfiltered && selector != "none") {
    unfiltered <- registry[registry$selector == "none", , drop = FALSE]
    missing <- is.na(route) & menu$can_predict
    route[missing] <- nrow(chosen) + match(menu$learner[missing], unfiltered$learner)
    chosen <- rbind(chosen, unfiltered)
  }
  chosen <- sb_add_dependency_status(chosen)[route, , drop = FALSE]
  candidates <- data.frame(name = menu$name, model_id = menu$model_id, learner = menu$learner,
    glmnet_alpha = menu$glmnet_alpha, can_predict = menu$can_predict,
    stage_selector = menu$stage_selector, requested_selector = selector,
    effective_selector = chosen$selector, tuning_semantics = chosen$tuning_semantics,
    dependency_status = chosen$dependency_status, missing_r_packages = chosen$missing_r_packages,
    status = "ready", reason = NA_character_, pipeline_id = NA_character_)
  pipelines <- list()
  for (i in seq_len(nrow(candidates))) {
    row <- chosen[i, , drop = FALSE]
    if (!menu$can_predict[i]) {
      candidates$status[i] <- "selector_only"
      candidates$reason[i] <- "Use this method as an upstream selector, not as a final predictor."
    } else if (is.na(route[i])) {
      candidates$status[i] <- "selector_not_admitted"
      candidates$reason[i] <- "This selector-learner pair is not admitted."
    } else if (identical(row$dependency_status, "missing_r_packages")) {
      candidates$status[i] <- "missing_dependency"
      candidates$reason[i] <- paste("R packages not found:", row$missing_r_packages)
    } else {
      args <- sb_binary_menu_arguments(row, selector_features, feature_set_id,
        selector_n_features, inner_folds, inner_seed)
      if (identical(row$selector, selector)) args$selector_options <- selector_options
      if (!is.na(menu$glmnet_alpha[i])) args$glmnet_alpha <- menu$glmnet_alpha[i]
      pipeline <- tryCatch(do.call(sb_pipeline, args), error = function(error) {
        if (!sb_is_expected_pipeline_constructor_error(error)) stop(error)
        error
      })
      if (inherits(pipeline, "error")) {
        candidates$status[i] <- "invalid_specification"
        candidates$reason[i] <- conditionMessage(pipeline)
      } else {
        pipelines[[menu$name[i]]] <- pipeline
        candidates$pipeline_id[i] <- pipeline$id
        if (row$selector != selector) {
          candidates$reason[i] <- "Explicitly included without an upstream selector; compare as a distinct route."
        }
      }
    }
  }
  list(pipelines = pipelines, candidates = candidates)
}

#' @rdname sb_models
#' @export
sb_binary_menu_pipelines <- function(selector_features, feature_set_id,
                                      selector_n_features = 10L,
                                      inner_folds = 5L, inner_seed = 1L) {
  selector_features <- sb_validate_fixed_features(selector_features)
  feature_set_id <- sb_validate_feature_set_id(feature_set_id)
  selector_n_features <- sb_check_positive_integer(
    selector_n_features, "selector_n_features", minimum = 2L
  )
  inner_folds <- sb_check_positive_integer(
    inner_folds, "inner_folds", minimum = 3L
  )
  inner_seed <- sb_check_positive_integer(inner_seed, "inner_seed", minimum = 1L)
  registry <- sb_pipeline_registry(task_type = "binary")
  registry <- registry[order(
    registry$selector != "none", registry$selector, registry$learner
  ), , drop = FALSE]
  pipelines <- lapply(seq_len(nrow(registry)), function(i) {
    args <- sb_binary_menu_arguments(
      registry[i, , drop = FALSE], selector_features, feature_set_id,
      selector_n_features, inner_folds, inner_seed
    )
    do.call(sb_pipeline, args)
  })
  menu <- data.frame(
    config_id = sprintf("binary_%04d", seq_along(pipelines)),
    selector = registry$selector,
    learner = registry$learner,
    stringsAsFactors = FALSE
  )
  reference <- data.frame(
    official_method = c("nb", "svmRadialWeights", "rf", "kknn",
      "adaboost", "LogitBoost", "cancerclass"),
    learner = c("kernel_naive_bayes_binary", "svm_radial_weighted_binary",
      "random_forest_binary", "kernel_weighted_knn_binary",
      "adaboost_binary", "logitboost_binary", "cancerclass_binary"),
    relationship = "algorithm_family_only",
    stringsAsFactors = FALSE
  )
  keys <- paste(menu$selector, menu$learner)
  reference$config_id <- menu$config_id[match(paste("none", reference$learner), keys)]
  menu$mime_family <- reference$official_method[match(menu$config_id, reference$config_id)]
  menu$mime_family[is.na(menu$mime_family)] <- ""
  menu$comparison_role <- ifelse(nzchar(menu$mime_family),
    "MIME_family_reference", "STORMBREAKER_extension")
  menu$pipeline <- pipelines
  attr(menu, "mime7_family_mapping") <- reference
  attr(menu, "registry_pairs") <- nrow(registry)
  menu
}

# ---- candidate-registry ----

sb_candidate_catalog_text <- function(value) {
  if (length(value) == 0L || is.na(value)) {
    return("")
  }
  as.character(value)
}

sb_candidate_stepcox_direction <- function(label) {
  matched <- regexec(
    "^StepCox\\[(forward|backward|both)\\]$",
    label
  )
  parts <- regmatches(label, matched)[[1L]]
  if (length(parts) != 2L) {
    return(NULL)
  }
  parts[[2L]]
}

sb_candidate_stepglm_direction <- function(label) {
  matched <- regexec(
    "^Stepglm\\[(forward|backward|both)\\]$",
    label
  )
  parts <- regmatches(label, matched)[[1L]]
  if (length(parts) != 2L) {
    return(NULL)
  }
  parts[[2L]]
}

sb_candidate_glmnet_alpha <- function(label) {
  if (identical(label, "Ridge")) {
    return(0)
  }
  if (identical(label, "Lasso")) {
    return(1)
  }
  if (!grepl("^Enet\\[alpha=[0-9.]+\\]$", label)) {
    return(NULL)
  }
  alpha <- as.numeric(sub(
    "^Enet\\[alpha=([0-9.]+)\\]$",
    "\\1",
    label
  ))
  if (is.na(alpha) || !is.finite(alpha) || alpha < 0 || alpha > 1) {
    return(NULL)
  }
  alpha
}

sb_candidate_runtime_parameter <- function(name, value) {
  paste0(
    name,
    "=",
    format(value, trim = TRUE, scientific = FALSE)
  )
}

sb_candidate_source_spec_audit <- function(task_type, model_type,
                                           learner_label, selector) {
  if (identical(selector, "stepcox")) {
    return(list(
      status = "audited_source_identity_ambiguous",
      reference = paste0(
        "Mime1 0.0.0.9000 ML.Dev.Prog.Sig.R and the historical Fig.5A ",
        "reconstruction ML296_Fig5_Backend.R; audited 2026-08-17 and ",
        "2026-08-24."
      ),
      difference = paste0(
        "The shared catalog label does not identify one operational StepCox ",
        "specification. Mime1 applies a common full-training univariable-Cox ",
        "prefilter, starts every direction from a full Cox model, and omits a ",
        "two-stage result when fewer than two variables remain. The historical ",
        "Fig.5A reconstruction uses different forward-start, minimum-feature, ",
        "truncation, warning, and fallback behavior. STORMBREAKER has no ",
        "implicit common prefilter, refits inside analysis resampling, starts ",
        "forward from the null model and backward/both from the full model, ",
        "and reports failures instead of silently substituting a result."
      )
    ))
  }

  direct_survival_glmnet <- identical(task_type, "survival") &&
    identical(model_type, "single") && identical(selector, "none") &&
    !is.null(sb_candidate_glmnet_alpha(learner_label))
  if (direct_survival_glmnet) {
    difference <- if (identical(learner_label, "Ridge")) {
      paste0(
        "Official Mime1 fits alpha=0 but selects lambda.min from a separate ",
        "cv.glmnet call that omits alpha and therefore uses glmnet's default ",
        "alpha=1; the historical Fig.5A reconstruction cross-validates ",
        "alpha=0 with a sample-size-dependent fold count."
      )
    } else if (identical(learner_label, "Lasso")) {
      paste0(
        "Official Mime1 cross-validates alpha=1 with fixed 10 folds and ",
        "lambda.min; the historical Fig.5A reconstruction uses alpha=1 with ",
        "a sample-size-dependent fold count."
      )
    } else {
      paste0(
        "Official Mime1 cross-validates the label alpha with fixed 10 folds ",
        "and lambda.min; the historical Fig.5A reconstruction uses the same ",
        "label alpha with a sample-size-dependent fold count."
      )
    }
    difference <- paste0(
      difference,
      " Mime1's wrapper also applies a common full-training univariable-Cox ",
      "prefilter, intersects predictors across supplied cohorts, and ",
      "mean-imputes validation cohorts. The catalog supplies only ",
      "historical_fig5a_label_inferred and no DOI or figure, so it does not ",
      "establish which implementation generated the row. STORMBREAKER ",
      "cross-validates the declared alpha with declared inner folds, seed, ",
      "measure, lambda rule, and standardization inside each outer analysis ",
      "fold, without an implicit common prefilter or fallback."
    )
    return(list(
      status = "audited_source_identity_ambiguous",
      reference = paste0(
        "Mime1 0.0.0.9000 ML.Dev.Prog.Sig.R lines 764-786, 2013-2057, ",
        "2424-2447, and 2638-2697; historical Fig.5A reconstruction ",
        "ML296_Fig5_Backend.R lines 323-350. Catalog evidence is ",
        "historical_fig5a_label_inferred with no source DOI or figure, so ",
        "source identity is not uniquely established; audited 2026-08-24."
      ),
      difference = difference
    ))
  }

  direct_survival_fixed <- identical(task_type, "survival") &&
    identical(model_type, "single") && identical(selector, "none") &&
    learner_label %in% c(
      "survival-SVM", "CoxBoost", "GBM", "SuperPC", "plsRcox", "Ranger",
      "GamBoost", "GlmBoost", "BlackBoost"
    )
  if (direct_survival_fixed) {
    ambiguous <- learner_label %in% c(
      "survival-SVM", "CoxBoost", "GBM", "SuperPC"
    )
    difference <- switch(learner_label,
      CoxBoost = paste0(
        "Official Mime1 optimizes a penalty from start.penalty=500 and then ",
        "selects the boosting step by fixed 10-fold cross-validation with ",
        "maxstepno=500. The historical Fig.5A reconstruction instead uses ",
        "maxstepno=200 and a sample-size-dependent fold count. STORMBREAKER ",
        "uses a declared fixed step count and either a declared penalty or ",
        "its documented analysis-event-count penalty rule, without an ",
        "automatic penalty or step search."
      ),
      GBM = paste0(
        "Official Mime1 searches 10,000 trees by fixed 10-fold ",
        "cross-validation and refits at the selected tree count. The ",
        "historical Fig.5A reconstruction searches 2,000 trees with a ",
        "sample-size-dependent fold count and different shrinkage, without ",
        "Mime1's second fit. STORMBREAKER uses declared fixed tree, depth, ",
        "node-size, shrinkage, bag-fraction, and seed parameters with no ",
        "early stopping or automatic tree-count search."
      ),
      SuperPC = paste0(
        "Official Mime1 selects a threshold by fixed 10-fold superpc.cv, ",
        "retries indefinitely after an error, and supplies validation ",
        "outcomes in each superpc.predict test object. The historical Fig.5A ",
        "reconstruction uses a sample-size-dependent fold count and no retry. ",
        "STORMBREAKER uses a declared fixed threshold and component count, ",
        "does not receive validation outcomes for assessment prediction, and ",
        "does not retry or substitute a fallback result."
      ),
      `survival-SVM` = paste0(
        "Official Mime1 calls survivalsvm with gamma.mu=1 and the package ",
        "defaults, then reports the raw predicted score. The historical ",
        "Fig.5A reconstruction makes the regression, linear-kernel, and ",
        "quadprog choices explicit but retains raw package scores and its ",
        "optional outer fallback. STORMBREAKER uses an explicit linear ",
        "regression SVM, fold-local standardization, and analysis-fold Cox ",
        "calibration to a log-hazard score, with failures reported instead ",
        "of fallback."
      ),
      plsRcox = paste0(
        "The 101/117 implementations select the requested component count ",
        "with cv.plsRcox before fitting plsRcox; the historical Fig.5A ",
        "reconstruction exposes that learner only behind a selector and also ",
        "adjusts the candidate component limit to the fitted feature count. ",
        "STORMBREAKER uses a declared fixed component count and performs no ",
        "automatic component search. The compatibility label alone does not ",
        "identify which operational definition the direct row should reproduce."
      ),
      Ranger = paste0(
        "The cited Table 1 identifies the Ranger label and ranger package but ",
        "does not specify tree count, mtry, node size, sampling, split rule, ",
        "seed, or risk-score semantics. Mime1's RSF uses randomForestSRC and ",
        "is not ranger. STORMBREAKER uses its own declared fixed ranger ",
        "survival parameters and native survival-curve prediction; neither ",
        "source establishes that this is the catalog row's implementation."
      ),
      GamBoost = paste0(
        "Official Mime1 does not implement GamBoost. The historical Fig.5A ",
        "backend uses mboost::gamboost with CoxPH, mstop=100, and nu=0.1, ",
        "whereas a later local all-in-one reconstruction incorrectly dispatches ",
        "its GamBoost label to glmboost. STORMBREAKER uses the distinct native ",
        "gamboost spline learner with those fixed boosting controls."
      ),
      GlmBoost = paste0(
        "Official Mime1 does not implement survival GlmBoost. Local historical ",
        "reconstructions use mboost::glmboost with CoxPH, mstop=100, and nu=0.1. ",
        "STORMBREAKER uses that distinct native linear boosting family with ",
        "explicit predictors and no automatic iteration search."
      ),
      BlackBoost = paste0(
        "Official Mime1 and the historical Fig.5A backend do not implement ",
        "BlackBoost. A later local all-in-one reconstruction uses ",
        "mboost::blackboost with CoxPH, mstop=100, and nu=0.1. STORMBREAKER ",
        "uses the native tree-boosting learner with those fixed controls and ",
        "explicit default tree settings."
      )
    )
    if (ambiguous) {
      difference <- paste0(
        difference,
        " Mime1 also applies a common full-training univariable-Cox prefilter, ",
        "intersects predictors across supplied cohorts, and mean-imputes each ",
        "validation cohort. The reconstructed wrapper and STORMBREAKER do not ",
        "share that complete operational procedure."
      )
    }
    reference <- if (ambiguous) {
      paste0(
        "Mime1 0.0.0.9000 ML.Dev.Prog.Sig.R all-mode lines 1521-2008 and ",
        "single-mode lines 2474-2634; historical Fig.5A reconstruction ",
        "ML296_Fig5_Backend.R lines 279-305 and 366-409. The bundled row ",
        "records historical_fig5a_label_inferred but no source DOI or figure, ",
        "so it does not identify one operational implementation; audited ",
        "2026-08-24."
      )
    } else if (identical(learner_label, "plsRcox")) {
      paste0(
        "Catalog-construction audit: direct plsRcox is the sole 101/117 ",
        "compatibility addition absent from the reconstructed 296 single-model ",
        "set. Comparison references are the Mime paper DOI ",
        "10.1016/j.csbj.2024.06.035, Mime1 0.0.0.9000 ",
        "ML.Dev.Prog.Sig.R lines 1907-1923 and 2501-2524, and historical ",
        "ML296_Fig5_Backend.R lines 279-305 and 373-377. The bundled row ",
        "records published_101_117_compatibility_candidate but no source DOI ",
        "or figure, so its exact source identity remains unverified; audited ",
        "2026-08-24."
      )
    } else if (learner_label %in% c(
      "GamBoost", "GlmBoost", "BlackBoost"
    )) {
      paste0(
        "Comparison references: historical ML296_Fig5_Backend.R lines ",
        "234-248 and 386-393, plus the later local ml296_all_in_one.R ",
        "mboost dispatcher and fixed parameter table. Official Mime1 ",
        "0.0.0.9000 ML.Dev.Prog.Sig contains no mboost learner. The bundled ",
        "row supplies no source DOI or figure, so its exact source identity ",
        "remains unverified; audited 2026-08-26."
      )
    } else {
      paste0(
        "Comparison reference only: npj Precision Oncology 2026 Table 1 ",
        "lists Random Forest (Ranger implementation), abbreviation Ranger, ",
        "and package ranger. The bundled row records ",
        "table1_catalog_candidate_endpoint_unverified but no source DOI or ",
        "figure; Table 1 supplies no operational settings, and official Mime1 ",
        "uses RSF/randomForestSRC rather than ranger. Source identity remains ",
        "unverified; audited 2026-08-24."
      )
    }
    return(list(
      status = if (ambiguous) {
        "audited_source_identity_ambiguous"
      } else {
        "audited_source_identity_unverified"
      },
      reference = reference,
      difference = difference
    ))
  }

  direct_binary_linear <- identical(task_type, "binary") &&
    identical(model_type, "single") && identical(selector, "none") &&
    (identical(learner_label, "LogisticRegression") ||
      !is.null(sb_candidate_glmnet_alpha(learner_label)))
  if (direct_binary_linear) {
    difference <- if (identical(learner_label, "LogisticRegression")) {
      paste0(
        "The catalog generator contributes only a normalized ",
        "LogisticRegression label; it does not specify an engine, formula ",
        "options, preprocessing, resampling, convergence handling, or ",
        "probability convention. STORMBREAKER fits an unpenalized ",
        "binomial-logit GLM on each analysis fold, requires a converged ",
        "full-rank fit with finite coefficients, and reports any warning or ",
        "invalid fit as a failure."
      )
    } else {
      paste0(
        "The catalog generator contributes only the Ridge, Lasso, or Enet ",
        "label and alpha; it does not specify an engine, lambda selection, ",
        "inner folds, seed, tuning measure, lambda rule, standardization, ",
        "preprocessing, or failure handling. STORMBREAKER cross-validates ",
        "the declared alpha inside each outer analysis fold using declared ",
        "inner folds, seed, measure, lambda rule, and standardization; any ",
        "supervised selector is refitted in the corresponding inner analysis ",
        "folds, and failures are reported without fallback."
      )
    }
    difference <- paste0(
      difference,
      " The catalog does not establish that this operational implementation ",
      "generated the candidate."
    )
    return(list(
      status = "audited_source_identity_unverified",
      reference = paste0(
        "Catalog construction reference: build_classification_candidate_space() ",
        "lines 609-681 and summary lines 1035-1044 in the 2026-07-18 ",
        "union-space generator, plus its algorithm-label provenance CSV. The ",
        "rows record union_generated_candidate_endpoint_unverified, ",
        "source_mime_binary=no, and no source DOI or figure. Official Mime1's ",
        "category module exposes neither LogisticRegression nor glmnet; source ",
        "identity remains unverified; audited 2026-08-24."
      ),
      difference = difference
    ))
  }

  mime_reference_binary <- identical(task_type, "binary") &&
    identical(model_type, "single") && identical(selector, "none") &&
    learner_label %in% c("NaiveBayes", "AdaBoost", "LogitBoost")
  if (!mime_reference_binary) {
    return(list(
      status = "not_yet_audited",
      reference = "",
      difference = ""
    ))
  }

  difference <- switch(learner_label,
    NaiveBayes = paste0(
      "Mime1 calls caret's nb method with usekernel=TRUE and tunes 25 ",
      "fL/adjust combinations by 5-fold, 10-repeat ROC. STORMBREAKER uses ",
      "fixed Gaussian Naive Bayes with empirical priors and a declared ",
      "variance-smoothing value, not the kernel model or Mime1 grid."
    ),
    AdaBoost = paste0(
      "Mime1 tunes 10 iteration/method combinations, including AdaBoost.M1 ",
      "and Real AdaBoost, by 5-fold, 10-repeat ROC. STORMBREAKER uses fixed ",
      "AdaBoost.M1 only and does not search the Mime1 grid."
    ),
    LogitBoost = paste0(
      "Mime1 tunes ten iteration counts from 11 through 101 by 5-fold, ",
      "10-repeat ROC. STORMBREAKER uses one caller-declared fixed iteration ",
      "count and does not search the Mime1 grid."
    )
  )
  difference <- paste0(
    difference,
    " Mime1's surrounding wrapper also intersects predictors across supplied ",
    "training and validation cohorts and mean-imputes each cohort; ",
    "STORMBREAKER does neither implicitly. The catalog does not establish ",
    "that Mime1 generated this candidate."
  )
  list(
    status = "audited_source_identity_unverified",
    reference = paste0(
      "Comparison reference only: Mime1 0.0.0.9000 ",
      "ML.Dev.Pred.Category.Sig.R lines 48-115 and 265-379. The catalog ",
      "rows provide no source DOI or figure and record ",
      "union_generated_candidate_endpoint_unverified, so source identity ",
      "remains unverified; audited 2026-08-24."
    ),
    difference = difference
  )
}

sb_candidate_learner_map <- function(label, task_type) {
  alpha <- sb_candidate_glmnet_alpha(label)
  if (!is.null(alpha)) {
    return(list(
      learner = if (identical(task_type, "survival")) {
        "glmnet_cox"
      } else {
        "glmnet_logistic"
      },
      parameters = sb_candidate_runtime_parameter("glmnet_alpha", alpha),
      configuration_required = TRUE,
      unsupported_shape = FALSE,
      missing = character()
    ))
  }

  if (identical(task_type, "survival")) {
    direction <- sb_candidate_stepcox_direction(label)
    if (!is.null(direction)) {
      return(list(
        learner = NA_character_,
        parameters = character(),
        configuration_required = FALSE,
        unsupported_shape = TRUE,
        missing = paste0("learner:", label)
      ))
    }
    map <- c(
      Coxph = "coxph",
      CoxBoost = "coxboost",
      plsRcox = "plsrcox",
      SuperPC = "superpc",
      `survival-SVM` = "survivalsvm",
      GBM = "gbm_cox",
      GamBoost = "gamboost_cox",
      GlmBoost = "glmboost_cox",
      BlackBoost = "blackboost_cox",
      Ranger = "ranger_survival",
      XGBoost = "xgboost_cox",
      SurvReg = "survreg_weibull",
      Rpart = "rpart_survival",
      RSF = "rsf_survival",
      Ctree = "ctree_survival",
      CForest = "cforest_survival"
    )
    configuration_required <- !identical(label, "Coxph")
  } else {
    direction <- sb_candidate_stepglm_direction(label)
    if (!is.null(direction)) {
      return(list(
        learner = "stepglm_binary",
        parameters = sb_candidate_runtime_parameter(
          "stepglm_direction",
          direction
        ),
        configuration_required = TRUE,
        unsupported_shape = FALSE,
        missing = character()
      ))
    }
    map <- c(
      LogisticRegression = "logistic",
      NaiveBayes = "gaussian_naive_bayes",
      AdaBoost = "adaboost_binary",
      LogitBoost = "logitboost_binary",
      LDA = "lda_binary",
      QDA = "qda_binary",
      KNN = "knn_binary",
      DecisionTree = "decision_tree_binary",
      RF = "random_forest_binary",
      SVM = "svm_radial_binary",
      svmRadialWeights = "svm_radial_weighted_binary",
      GBM = "gbm_binary",
      glmBoost = "glmboost_binary",
      plsRglm = "plsrglm_binary",
      XGBoost = "xgboost_binary"
    )
    configuration_required <- !identical(label, "LogisticRegression")
  }

  learner <- unname(map[label])
  if (length(learner) == 0L || is.na(learner)) {
    return(list(
      learner = NA_character_,
      parameters = character(),
      configuration_required = FALSE,
      unsupported_shape = FALSE,
      missing = paste0("learner:", label)
    ))
  }
  list(
    learner = learner,
    parameters = character(),
    configuration_required = configuration_required,
    unsupported_shape = FALSE,
    missing = character()
  )
}

sb_candidate_selector_map <- function(label, task_type) {
  if (!nzchar(label)) {
    return(list(
      selector = "none",
      parameters = character(),
      missing = character()
    ))
  }
  direction <- if (identical(task_type, "survival")) {
    sb_candidate_stepcox_direction(label)
  } else {
    NULL
  }
  if (!is.null(direction)) {
    return(list(
      selector = "stepcox",
      parameters = sb_candidate_runtime_parameter(
        "stepcox_direction",
        direction
      ),
      missing = character()
    ))
  }
  list(
    selector = NA_character_,
    parameters = character(),
    missing = paste0("selector:", label)
  )
}

sb_candidate_mapping_result <- function(task_type, route_status,
                                        admission_status, reason,
                                        selector = NA_character_,
                                        learner = NA_character_,
                                        parameters = NA_character_,
                                        route_available = FALSE,
                                        configuration_status = "not_applicable",
                                        source_spec_status = "not_applicable",
                                        source_spec_reference = "",
                                        source_spec_difference = "",
                                        missing = NA_character_) {
  data.frame(
    runtime_task_type = task_type,
    runtime_selector = selector,
    runtime_learner = learner,
    runtime_parameters = parameters,
    runtime_route_status = route_status,
    runtime_route_available = route_available,
    runtime_configuration_status = configuration_status,
    runtime_source_spec_status = source_spec_status,
    runtime_source_spec_reference = source_spec_reference,
    runtime_source_spec_difference = source_spec_difference,
    missing_runtime_components = missing,
    runtime_admission_status = admission_status,
    runtime_mapping_reason = reason,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

sb_candidate_mapping_row <- function(row, route_keys) {
  task_type <- if (identical(
    sb_candidate_catalog_text(row$task_type),
    "survival_right_censored"
  )) {
    "survival"
  } else {
    "binary"
  }
  admission_status <- sb_candidate_catalog_text(row$admission_status)
  if (grepl("^blocked", admission_status)) {
    return(sb_candidate_mapping_result(
      task_type,
      "methodologically_blocked",
      "not_admitted_methodologically_blocked",
      paste0(
        "The source catalog blocks this candidate until its stated ",
        "methodological requirement is implemented."
      )
    ))
  }
  if (identical(
    sb_candidate_catalog_text(row$role_semantics_status),
    "stage_roles_unresolved"
  )) {
    return(sb_candidate_mapping_result(
      task_type,
      "stage_roles_unresolved",
      "not_admitted_stage_roles_unresolved",
      paste0(
        "The source catalog does not resolve stage roles; STORMBREAKER ",
        "does not infer selector or learner roles from stage order."
      )
    ))
  }

  model_type <- sb_candidate_catalog_text(row$model_type)
  learner_label <- sb_candidate_catalog_text(row$learner)
  if (identical(model_type, "single") &&
    identical(task_type, "survival")) {
    direction <- sb_candidate_stepcox_direction(learner_label)
  } else {
    direction <- NULL
  }
  if (!is.null(direction)) {
    selector <- list(
      selector = "stepcox",
      parameters = sb_candidate_runtime_parameter(
        "stepcox_direction",
        direction
      ),
      missing = character()
    )
    learner <- list(
      learner = "coxph",
      parameters = character(),
      configuration_required = FALSE,
      unsupported_shape = FALSE,
      missing = character()
    )
  } else if (identical(model_type, "single")) {
    selector <- sb_candidate_selector_map("", task_type)
    learner <- sb_candidate_learner_map(learner_label, task_type)
  } else if (identical(model_type, "selector_learner")) {
    selector <- sb_candidate_selector_map(
      sb_candidate_catalog_text(row$selector_1),
      task_type
    )
    learner <- sb_candidate_learner_map(learner_label, task_type)
  } else {
    return(sb_candidate_mapping_result(
      task_type,
      "unsupported_pipeline_shape",
      "not_admitted_runtime_route_unavailable",
      "The current runtime does not represent this catalog pipeline shape."
    ))
  }

  parameters <- c(selector$parameters, learner$parameters)
  parameters <- if (length(parameters) == 0L) {
    ""
  } else {
    paste(parameters, collapse = "; ")
  }
  if (learner$unsupported_shape) {
    return(sb_candidate_mapping_result(
      task_type = task_type,
      route_status = "unsupported_pipeline_shape",
      admission_status = "not_admitted_runtime_route_unavailable",
      reason = paste0(
        "The catalog learner itself performs StepCox selection; the current ",
        "one-selector/one-learner grammar does not silently collapse that ",
        "additional stage."
      ),
      selector = selector$selector,
      parameters = parameters,
      missing = learner$missing
    ))
  }

  missing <- c(selector$missing, learner$missing)
  if (length(missing) > 0L) {
    return(sb_candidate_mapping_result(
      task_type = task_type,
      route_status = "component_unavailable",
      admission_status = "not_admitted_runtime_route_unavailable",
      reason = paste0(
        "No current runtime adapter is assigned to: ",
        paste(missing, collapse = "; "),
        "."
      ),
      selector = selector$selector,
      learner = learner$learner,
      parameters = parameters,
      missing = paste(missing, collapse = "; ")
    ))
  }

  route_key <- paste(
    task_type,
    selector$selector,
    learner$learner,
    sep = "::"
  )
  route_available <- route_key %in% route_keys
  if (!route_available) {
    return(sb_candidate_mapping_result(
      task_type = task_type,
      route_status = "route_not_admitted",
      admission_status = "not_admitted_runtime_route_unavailable",
      reason = paste0(
        "Mapped components exist, but their task-aware selector-learner ",
        "route is not admitted by the current runtime."
      ),
      selector = selector$selector,
      learner = learner$learner,
      parameters = parameters
    ))
  }

  source_spec <- sb_candidate_source_spec_audit(
    task_type,
    model_type,
    learner_label,
    selector$selector
  )
  source_identity_ambiguous <- identical(
    source_spec$status,
    "audited_source_identity_ambiguous"
  )
  source_identity_unverified <- identical(
    source_spec$status,
    "audited_source_identity_unverified"
  )

  sb_candidate_mapping_result(
    task_type = task_type,
    route_status = "route_available",
    admission_status = if (source_identity_ambiguous) {
      "not_admitted_source_identity_ambiguous"
    } else if (source_identity_unverified) {
      "not_admitted_source_identity_unverified"
    } else {
      "not_admitted_source_to_spec_unverified"
    },
    reason = if (source_identity_ambiguous) {
      paste0(
        "A current runtime route matches the task and mapped stage families, ",
        "but the shared historical label does not identify one ",
        "operational source specification."
      )
    } else if (source_identity_unverified) {
      paste0(
        "A current runtime route matches the task and named learner family, ",
        "but the candidate row does not identify one operational source ",
        "specification; any documented implementation remains only a ",
        "comparison reference."
      )
    } else {
      paste0(
        "A current runtime route matches this candidate's task and mapped ",
        "stage families, but source-to-spec semantics and any required ",
        "configuration are not yet verified."
      )
    },
    selector = selector$selector,
    learner = learner$learner,
    parameters = parameters,
    route_available = TRUE,
    configuration_status = if (learner$configuration_required) {
      "additional_parameters_required"
    } else {
      "no_additional_runtime_parameters"
    },
    source_spec_status = source_spec$status,
    source_spec_reference = source_spec$reference,
    source_spec_difference = source_spec$difference
  )
}

sb_candidate_registry_data <- function() {
  catalog <- utils::read.csv(
    system.file(
      "extdata",
      "unified_pipeline_registry_1213_20260806.csv",
      package = "stormbreaker"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  routes <- sb_pipeline_registry_data()
  route_keys <- paste(
    routes$task_type,
    routes$selector,
    routes$learner,
    sep = "::"
  )
  mapping <- lapply(seq_len(nrow(catalog)), function(index) {
    sb_candidate_mapping_row(catalog[index, , drop = FALSE], route_keys)
  })
  mapping <- do.call(rbind, mapping)
  rownames(mapping) <- NULL
  result <- data.frame(
    catalog,
    mapping,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(result) <- NULL
  result
}

sb_candidate_route_statuses <- function() {
  c(
    "route_available",
    "component_unavailable",
    "route_not_admitted",
    "unsupported_pipeline_shape",
    "stage_roles_unresolved",
    "methodologically_blocked"
  )
}

#' @rdname sb_models
#' @export
sb_candidate_registry <- function(task_type = NULL, route_status = NULL) {
  if (!is.null(task_type) && (!is.character(task_type) ||
    length(task_type) != 1L || is.na(task_type) ||
    !task_type %in% c("survival", "binary"))) {
    sb_abort(
      "`task_type` must be `NULL`, `\"survival\"`, or `\"binary\"`.",
      "stormbreaker_bad_candidate_registry_task"
    )
  }
  if (!is.null(route_status) && (!is.character(route_status) ||
    length(route_status) != 1L || is.na(route_status) ||
    !route_status %in% sb_candidate_route_statuses())) {
    sb_abort(
      "`route_status` must be `NULL` or one documented route status.",
      "stormbreaker_bad_candidate_registry_route_status"
    )
  }

  registry <- sb_candidate_registry_data()
  if (!is.null(task_type)) {
    registry <- registry[
      registry$runtime_task_type == task_type, ,
      drop = FALSE
    ]
  }
  if (!is.null(route_status)) {
    registry <- registry[
      registry$runtime_route_status == route_status, ,
      drop = FALSE
    ]
  }
  rownames(registry) <- NULL
  registry
}

# ---- survival-menu-pipelines ----

sb_survival_menu_arguments <- function(row, selector_features, feature_set_id,
                                       selector_n_features) {
  args <- list(selector = row$selector, learner = row$learner)
  parameter <- row$selector_parameter
  if (identical(parameter, "selector_n_features")) {
    args$selector_n_features <- selector_n_features
  }
  if (parameter %in% c(
    "stepcox_direction", "selector_rsf_nodesize_seed_and_stepcox_direction",
    "selector_inner_folds_seed_and_stepcox_direction"
  )) {
    args$stepcox_direction <- "both"
  }
  if (parameter %in% c(
    "selector_rsf_nodesize_and_seed",
    "selector_rsf_nodesize_seed_and_stepcox_direction"
  )) {
    args$selector_rsf_nodesize <- 5L
    args$selector_seed <- 1L
  }
  if (parameter %in% c(
    "selector_inner_folds_and_seed",
    "selector_inner_folds_seed_and_stepcox_direction"
  )) {
    args$selector_inner_folds <- 10L
    args$selector_seed <- 1L
  }
  if (identical(parameter, "selector_features_and_feature_set_id")) {
    args$selector_features <- selector_features
    args$feature_set_id <- feature_set_id
  }
  if (identical(row$tuning_semantics, "lambda_selected_by_inner_cv")) {
    args$inner_folds <- 10L
    args$inner_seed <- 1L
    args$glmnet_standardize <- TRUE
    args$glmnet_lambda_rule <- "lambda.min"
    args$glmnet_tuning_measure <- "deviance"
  }
  if (identical(row$tuning_semantics, "analysis_fold_internal_cv")) {
    args$learner_tuning_folds <- if (row$learner == "plsrcox_cv") 5L else 10L
    args$learner_seed <- 1L
    if (row$learner == "gbm_cox_cv") {
      args$gbm_cv_max_trees <- 10000L
    }
  }
  args
}

#' @rdname sb_models
#' @export
sb_survival_menu_pipelines <- function(selector_features, feature_set_id,
                                       selector_n_features = 30L) {
  selector_features <- sb_validate_fixed_features(selector_features)
  feature_set_id <- sb_validate_feature_set_id(feature_set_id)
  selector_n_features <- sb_check_positive_integer(
    selector_n_features, "selector_n_features", minimum = 2L
  )
  registry <- sb_pipeline_registry(task_type = "survival")
  registry <- registry[order(
    registry$selector != "none", registry$selector, registry$learner
  ), , drop = FALSE]
  pipelines <- lapply(seq_len(nrow(registry)), function(i) {
    args <- sb_survival_menu_arguments(
      registry[i, , drop = FALSE], selector_features, feature_set_id,
      selector_n_features
    )
    do.call(sb_pipeline, args)
  })
  reference <- sb_mime117_pipelines(
    inner_folds = 10L, inner_seed = 1L,
    selector_inner_folds = 10L, selector_seed = 1L, rsf_nodesize = 5L,
    learner_tuning_folds = 10L, plsrcox_tuning_folds = 5L,
    learner_seed = 1L, gbm_cv_max_trees = 10000L
  )
  pipelines <- c(pipelines, reference$pipeline)
  ids <- vapply(pipelines, `[[`, "", "id")
  pipelines <- pipelines[!duplicated(ids)]
  ids <- vapply(pipelines, `[[`, "", "id")
  config_ids <- sprintf("aligned_%04d", seq_along(pipelines))
  mapping <- reference[setdiff(names(reference), "pipeline")]
  mapping$comparison_config_id <- config_ids[match(
    vapply(reference$pipeline, `[[`, "", "id"), ids
  )]
  menu <- data.frame(
    config_id = config_ids,
    selector = vapply(pipelines, `[[`, "", "selector"),
    learner = vapply(pipelines, `[[`, "", "learner")
  )
  menu$mime_labels <- vapply(config_ids, function(id) {
    paste(mapping$official_model_label[mapping$comparison_config_id == id],
      collapse = " | "
    )
  }, "")
  menu$comparison_role <- ifelse(
    nzchar(menu$mime_labels), "MIME_structure_matched", "STORMBREAKER_extension"
  )
  menu$pipeline <- pipelines
  attr(menu, "mime117_mapping") <- mapping
  attr(menu, "registry_pairs") <- nrow(registry)
  menu
}

# ---- mime117-pipelines ----

sb_mime117_alpha_text <- function(alpha) {
  format(alpha, nsmall = 1L, trim = TRUE, scientific = FALSE)
}

sb_mime117_component_label <- function(family, alpha = NA_real_,
                                        direction = NA_character_,
                                        single = FALSE) {
  if (identical(family, "Enet")) {
    return(paste0("Enet[\u03b1=", sb_mime117_alpha_text(alpha), "]"))
  }
  if (identical(family, "StepCox")) {
    return(paste0("StepCox[", direction, "]"))
  }
  if (identical(family, "survival-SVM") && single) {
    return("survival - SVM")
  }
  family
}

sb_mime117_variant_text <- function(family, alpha, direction) {
  if (identical(family, "Enet")) {
    return(paste0("alpha=", sb_mime117_alpha_text(alpha)))
  }
  if (identical(family, "StepCox")) {
    return(paste0("direction=", direction))
  }
  ""
}

sb_mime117_source_grid <- function() {
  alphas <- seq(0.1, 0.9, by = 0.1)
  directions <- c("both", "backward", "forward")
  rows <- list()

  add <- function(block, stage1_family, stage2_family = "",
                  alpha = NA_real_, direction = NA_character_) {
    has_stage2 <- nzchar(stage2_family)
    stage1_alpha <- if (identical(stage1_family, "Enet")) {
      alpha
    } else {
      NA_real_
    }
    stage2_alpha <- if (identical(stage2_family, "Enet")) {
      alpha
    } else {
      NA_real_
    }
    stage1_direction <- if (identical(stage1_family, "StepCox")) {
      direction
    } else {
      NA_character_
    }
    stage2_direction <- if (identical(stage2_family, "StepCox")) {
      direction
    } else {
      NA_character_
    }
    stage1_label <- sb_mime117_component_label(
      stage1_family,
      alpha = stage1_alpha,
      direction = stage1_direction,
      single = !has_stage2
    )
    label <- if (has_stage2) {
      paste(
        stage1_label,
        sb_mime117_component_label(
          stage2_family,
          alpha = stage2_alpha,
          direction = stage2_direction
        ),
        sep = " + "
      )
    } else {
      stage1_label
    }

    rows[[length(rows) + 1L]] <<- data.frame(
      config_id = "",
      official_order = 0L,
      block = block,
      pipeline_type = if (has_stage2) "two_stage" else "single",
      official_model_label = label,
      stage1_family = stage1_family,
      stage1_variant = sb_mime117_variant_text(
        stage1_family,
        stage1_alpha,
        stage1_direction
      ),
      stage2_family = stage2_family,
      stage2_variant = if (has_stage2) {
        sb_mime117_variant_text(
          stage2_family,
          stage2_alpha,
          stage2_direction
        )
      } else {
        ""
      },
      enet_alpha = if (identical(stage1_family, "Enet") ||
        identical(stage2_family, "Enet")) {
        alpha
      } else {
        NA_real_
      },
      stepcox_direction = if (identical(stage1_family, "StepCox") ||
        identical(stage2_family, "StepCox")) {
        direction
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  }

  add("RSF", "RSF")
  add("RSF", "RSF", "CoxBoost")
  for (alpha in alphas) {
    add("RSF", "RSF", "Enet", alpha = alpha)
  }
  for (stage2 in c("GBM", "Lasso", "plsRcox", "Ridge")) {
    add("RSF", "RSF", stage2)
  }
  for (direction in directions) {
    add("RSF", "RSF", "StepCox", direction = direction)
  }
  add("RSF", "RSF", "SuperPC")
  add("RSF", "RSF", "survival-SVM")

  for (alpha in alphas) {
    add("Enet", "Enet", alpha = alpha)
  }

  for (direction in directions) {
    add("StepCox", "StepCox", direction = direction)
  }
  for (direction in directions) {
    add("StepCox", "StepCox", "CoxBoost", direction = direction)
    for (alpha in alphas) {
      add(
        "StepCox",
        "StepCox",
        "Enet",
        alpha = alpha,
        direction = direction
      )
    }
    for (stage2 in c(
      "GBM", "Lasso", "plsRcox", "Ridge", "RSF", "SuperPC",
      "survival-SVM"
    )) {
      add("StepCox", "StepCox", stage2, direction = direction)
    }
  }

  add("CoxBoost", "CoxBoost")
  for (alpha in alphas) {
    add("CoxBoost", "CoxBoost", "Enet", alpha = alpha)
  }
  for (stage2 in c("GBM", "Lasso", "plsRcox", "Ridge")) {
    add("CoxBoost", "CoxBoost", stage2)
  }
  for (direction in directions) {
    add("CoxBoost", "CoxBoost", "StepCox", direction = direction)
  }
  add("CoxBoost", "CoxBoost", "SuperPC")
  add("CoxBoost", "CoxBoost", "survival-SVM")

  add("plsRcox", "plsRcox")
  add("SuperPC", "SuperPC")
  add("GBM", "GBM")
  add("survival-SVM", "survival-SVM")
  add("Ridge", "Ridge")

  add("Lasso", "Lasso")
  for (stage2 in c("CoxBoost", "GBM", "plsRcox", "RSF")) {
    add("Lasso", "Lasso", stage2)
  }
  for (direction in directions) {
    add("Lasso", "Lasso", "StepCox", direction = direction)
  }
  add("Lasso", "Lasso", "SuperPC")
  add("Lasso", "Lasso", "survival-SVM")

  grid <- do.call(rbind, rows)
  rownames(grid) <- NULL
  grid$official_order <- seq_len(nrow(grid))
  grid$config_id <- sprintf("MIME117-%03d", grid$official_order)
  grid
}

sb_mime117_learner <- function(family) {
  switch(
    family,
    RSF = "rsf_survival",
    Enet = "glmnet_cox",
    Lasso = "glmnet_cox",
    Ridge = "glmnet_cox",
    CoxBoost = "coxboost_cv",
    GBM = "gbm_cox_cv",
    plsRcox = "plsrcox_cv",
    SuperPC = "superpc_cv",
    `survival-SVM` = "survivalsvm",
    StepCox = "coxph"
  )
}

sb_mime117_selector <- function(stage1_family, final_family,
                                has_stage2) {
  if (!has_stage2) {
    if (identical(stage1_family, "StepCox")) "stepcox" else "none"
  } else if (identical(final_family, "StepCox")) {
    switch(
      stage1_family,
      RSF = "rsf_then_stepcox",
      CoxBoost = "coxboost_then_stepcox",
      Lasso = "lasso_then_stepcox"
    )
  } else {
    switch(
      stage1_family,
      RSF = "rsf_minimal_depth",
      StepCox = "stepcox",
      CoxBoost = "coxboost_nonzero",
      Lasso = "lasso_nonzero"
    )
  }
}

sb_mime117_pipeline_from_row <- function(row, inner_folds, inner_seed,
                                         selector_inner_folds,
                                         selector_seed, rsf_nodesize,
                                         learner_tuning_folds,
                                         plsrcox_tuning_folds,
                                         learner_seed,
                                         gbm_cv_max_trees) {
  has_stage2 <- identical(row$pipeline_type[[1L]], "two_stage")
  final_family <- if (has_stage2) {
    row$stage2_family[[1L]]
  } else {
    row$stage1_family[[1L]]
  }
  selector <- sb_mime117_selector(
    row$stage1_family[[1L]],
    final_family,
    has_stage2
  )
  learner <- sb_mime117_learner(final_family)
  args <- list(selector = selector, learner = learner)

  if (selector %in% c("stepcox", sb_mime_composite_selector_ids())) {
    args$stepcox_direction <- row$stepcox_direction[[1L]]
  }
  if (selector %in% sb_mime_rsf_selector_ids()) {
    args$selector_rsf_nodesize <- rsf_nodesize
    args$selector_seed <- selector_seed
  } else if (selector %in% sb_mime_tuned_selector_ids()) {
    args$selector_inner_folds <- selector_inner_folds
    args$selector_seed <- selector_seed
  }
  if (identical(learner, "glmnet_cox")) {
    alpha <- switch(
      final_family,
      Enet = row$enet_alpha[[1L]],
      Lasso = 1,
      Ridge = 0
    )
    args$inner_folds <- inner_folds
    args$inner_seed <- inner_seed
    args$glmnet_alpha <- alpha
    args$glmnet_standardize <- TRUE
    args$glmnet_lambda_rule <- "lambda.min"
    args$glmnet_tuning_measure <- "deviance"
  }
  if (learner %in% sb_analysis_tuned_survival_learner_ids()) {
    args$learner_tuning_folds <- if (identical(learner, "plsrcox_cv")) {
      plsrcox_tuning_folds
    } else {
      learner_tuning_folds
    }
    args$learner_seed <- learner_seed
    if (identical(learner, "gbm_cox_cv")) {
      args$gbm_cv_max_trees <- gbm_cv_max_trees
    }
  }
  do.call(sb_pipeline, args)
}

#' @rdname sb_models
#' @export
sb_mime117_pipelines <- function(inner_folds = 10L,
                                 inner_seed = 5201314L,
                                 selector_inner_folds = 10L,
                                 selector_seed = 5201314L,
                                 rsf_nodesize = 5L,
                                 learner_tuning_folds = 10L,
                                 plsrcox_tuning_folds = 5L,
                                 learner_seed = 5201314L,
                                 gbm_cv_max_trees = 10000L) {
  inner_folds <- sb_check_positive_integer(
    inner_folds,
    "inner_folds",
    minimum = 3L
  )
  inner_seed <- sb_check_positive_integer(inner_seed, "inner_seed")
  selector_inner_folds <- sb_check_positive_integer(
    selector_inner_folds,
    "selector_inner_folds",
    minimum = 3L
  )
  selector_seed <- sb_check_positive_integer(selector_seed, "selector_seed")
  rsf_nodesize <- sb_check_positive_integer(rsf_nodesize, "rsf_nodesize")
  learner_tuning_folds <- sb_check_positive_integer(
    learner_tuning_folds,
    "learner_tuning_folds",
    minimum = 3L
  )
  plsrcox_tuning_folds <- sb_check_positive_integer(
    plsrcox_tuning_folds,
    "plsrcox_tuning_folds",
    minimum = 3L
  )
  learner_seed <- sb_check_positive_integer(learner_seed, "learner_seed")
  gbm_cv_max_trees <- sb_check_positive_integer(
    gbm_cv_max_trees,
    "gbm_cv_max_trees"
  )

  suite <- sb_mime117_source_grid()
  pipelines <- lapply(seq_len(nrow(suite)), function(index) {
    sb_mime117_pipeline_from_row(
      suite[index, , drop = FALSE],
      inner_folds = inner_folds,
      inner_seed = inner_seed,
      selector_inner_folds = selector_inner_folds,
      selector_seed = selector_seed,
      rsf_nodesize = rsf_nodesize,
      learner_tuning_folds = learner_tuning_folds,
      plsrcox_tuning_folds = plsrcox_tuning_folds,
      learner_seed = learner_seed,
      gbm_cv_max_trees = gbm_cv_max_trees
    )
  })
  suite$operational_fidelity <- paste0(
    "fold_local_mime_training_reproduction_not_bug_for_bug"
  )
  suite[["pipeline"]] <- pipelines
  suite
}
