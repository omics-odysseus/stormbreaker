#' STORMBREAKER overview
#' @name stormbreaker-package
#' @description Build, compare and interpret transcriptomic prediction models.
#' @rawRd \arguments{
#' \item{view}{
#' Open the local page in a browser; defaults to interactive().
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{Package overview}{
#' STORMBREAKER stands for Systematic transcriptomic Risk Modeling and
#' Benchmarking for Reproducible, Evidence-Aware biomarker Evaluation in R.
#' It combines feature selection, model comparison, final fitting and external
#' evaluation for numeric transcriptomic predictors. The main workflows support
#' binary and right-censored survival outcomes; \code{\link[=multi_run]{multi_run()}} compares the
#' implemented multiclass learners.
#' 
#' \strong{Core workflow:} 
#' 
#' \itemize{
#' \item \code{\link[=bi_run]{bi_run()}} and \code{\link[=surv_run]{surv_run()}} compare candidate pipelines on shared folds.
#' Data contain sample IDs, outcomes and numeric predictors; survival data
#' also specify a time unit and, when needed, a prediction horizon.
#' \item \code{\link[=finish_run]{finish_run()}} fits a pipeline chosen from the development comparison and
#' evaluates named external cohorts. \code{\link[=predict_final]{predict_final()}} applies the saved model
#' to new samples using their IDs and required predictors.
#' \item \code{\link[=clinical_plots]{clinical_plots()}} draws evaluation results; \code{\link[=plots]{plots()}} draws algorithm and
#' feature-selection diagnostics; \code{\link[=model_shap]{model_shap()}} computes prediction explanations.
#' \item \code{\link[=export_run]{export_run()}} saves tables and figures. Use \code{\link[base:readRDS]{base::saveRDS()}} and
#' \code{\link[base:readRDS]{base::readRDS()}} to retain and reload the fitted model.
#' }
#' 
#' 
#' \strong{Evaluation:} 
#' 
#' Supervised feature selection and tuning use each resampling analysis set.
#' Declare patient groups for repeated observations; \code{\link[=sb_run_internal_external_cv]{sb_run_internal_external_cv()}}
#' holds out complete cohorts. Choose the final pipeline using development results.
#' 
#' \code{\link[=finish_run]{finish_run()}} and \code{\link[=sb_evaluate_external]{sb_evaluate_external()}} default to \code{reverse = FALSE},
#' keeping the development-defined direction fixed. For binary outcomes,
#' explicit \code{reverse = TRUE} uses external outcomes to flip probabilities when
#' AUC is below 0.5 and retains original predictions and correction records.
#' Survival evaluation requires \code{FALSE}.
#' 
#' Fold-to-fold variability describes resampling variation, not a confidence
#' interval. \code{\link[=sb_bootstrap_external]{sb_bootstrap_external()}} estimates external-patient sampling
#' uncertainty conditional on the saved model and predictions. Survival metrics
#' use the declared time unit, horizon and censoring assumptions.
#' Supply predictors on the study's preprocessing scale. Data-dependent
#' preprocessing belongs within the relevant analysis folds; glmnet's optional
#' internal standardization is fitted within each analysis set.
#' 
#' 
#' \strong{Formula and feature bridges:} 
#' 
#' \code{\link[=sb_linear_signature]{sb_linear_signature()}} records a supplied linear formula with fixed
#' coefficients for prediction. \code{\link[=sb_apply_feature_map]{sb_apply_feature_map()}} applies a one-to-one
#' feature mapping before prediction. Prepare annotation matches and any
#' expression aggregation or preprocessing before using these interfaces.
#' 
#' 
#' \strong{Scope:} 
#' 
#' Use \code{\link[=models]{models()}} to inspect the binary, survival and multiclass method menus
#' and their dependencies. Multiclass analysis uses \code{\link[=multi_run]{multi_run()}},
#' \code{\link[=analyze_methods]{analyze_methods()}} or \code{\link[=sb_compare_multiclass]{sb_compare_multiclass()}}; the binary and survival
#' selector-chain workflows have separate interfaces. Predictors are numeric,
#' and survival outcomes use right censoring rather than competing-risk coding.
#' 
#' }
#' 
#' \subsection{\code{sb_start()}}{
#' Open the local getting-started page
#' 
#' \strong{Returns:} Invisibly the path to a local HTML page with modeling and validation workflows and
#' copyable R examples. It does not read patient data, fit models, install
#' packages or upload anything. Examples are read from the installed package.
#' }
#' }
#' @rawRd \references{
#' McLernon DJ, Giardiello D, Van Calster B, et al. (2023).
#' Assessing performance and clinical usefulness in prediction models with
#' survival outcomes: practical guidance for Cox proportional hazards models.
#' \emph{Annals of Internal Medicine}, 176(1), 105-114.
#' \doi{10.7326/M22-0844}.
#' 
#' Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ (2011). On the C-statistics
#' for evaluating overall adequacy of risk prediction procedures with
#' censored survival data. \emph{Statistics in Medicine}, 30(10), 1105-1117.
#' \doi{10.1002/sim.4154}.
#' }
NULL

# ---- core-conditions ----

# ---- adapter-conditions ----

sb_validate_predictor_frame <- function(x, argument) {
  if (!is.data.frame(x) || nrow(x) < 1L || ncol(x) < 1L) {
    sb_abort(
      paste0("`", argument, "` must be a non-empty data frame."),
      "stormbreaker_bad_adapter_predictors"
    )
  }
  if (anyDuplicated(names(x)) || any(!nzchar(names(x)))) {
    sb_abort(
      paste0("`", argument, "` must have unique, non-empty column names."),
      "stormbreaker_bad_adapter_predictor_names"
    )
  }
  numeric_columns <- vapply(x, is.numeric, logical(1))
  if (!all(numeric_columns)) {
    sb_abort(
      paste0("`", argument, "` must contain only numeric columns."),
      "stormbreaker_bad_adapter_predictor_type"
    )
  }
  finite_columns <- vapply(
    x,
    function(column) !anyNA(column) && all(is.finite(column)),
    logical(1)
  )
  if (!all(finite_columns)) {
    sb_abort(
      paste0("`", argument, "` contains missing or non-finite values."),
      "stormbreaker_bad_adapter_predictor_values"
    )
  }
  x
}

sb_validate_endpoint_vector <- function(value, argument, n) {
  if (!is.numeric(value) || length(value) != n || anyNA(value) ||
      any(!is.finite(value))) {
    sb_abort(
      paste0(
        "`", argument,
        "` must contain one finite numeric value per analysis patient."
      ),
      "stormbreaker_bad_adapter_endpoint"
    )
  }
  value
}

# ---- comparison-progress ----

sb_check_progress <- function(progress) {
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    sb_abort("`progress` must be TRUE or FALSE.", "stormbreaker_bad_progress")
  }
}

sb_progress_message <- function(progress, ...) {
  if (isTRUE(progress)) {
    message("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  }
  invisible(NULL)
}

# ---- conditions ----

sb_abort <- function(message, subclass, parent = NULL) {
  condition <- structure(
    list(message = message, call = NULL, parent = parent),
    class = c(subclass, "stormbreaker_error", "error", "condition")
  )
  stop(condition)
}

sb_check_name <- function(value, argument) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    sb_abort(
      paste0("`", argument, "` must be one non-empty column name."),
      "stormbreaker_bad_column_argument"
    )
  }
}

sb_quote_name <- function(value) {
  if (grepl("`", value, fixed = TRUE)) {
    sb_abort(
      "Column names containing backticks are not supported.",
      "stormbreaker_unsupported_column_name"
    )
  }
  paste0("`", value, "`")
}

sb_check_positive_integer <- function(value, argument, minimum = 1L) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < minimum ||
      value > .Machine$integer.max || value != as.integer(value)) {
    sb_abort(
      paste0("`", argument, "` must be an integer >= ", minimum, "."),
      "stormbreaker_bad_integer"
    )
  }
  as.integer(value)
}

# ---- retired-learners ----

# Retirement concerns the binary/survival learners, not similarly named engines.
sb_retired_learner_ids <- function() {
  c("gam_cox_reml", "gam_binary_reml", "slope_binary")
}

sb_reject_retired_learner <- function(learner) {
  if (is.character(learner) && any(learner %in% sb_retired_learner_ids())) {
    sb_abort(
      paste0(
        "Learner `", paste(learner[learner %in% sb_retired_learner_ids()], collapse = ", "),
        "` was removed in 0.3.38 because of excessive runtime and computational cost."
      ),
      "stormbreaker_retired_learner"
    )
  }
  invisible(NULL)
}

# ---- runtime-dependencies ----

sb_dependency_components <- function(dependency) {
  if (is.na(dependency) || !nzchar(dependency)) {
    return(list(r_packages = character(), python_modules = character()))
  }
  parts <- trimws(strsplit(dependency, " + ", fixed = TRUE)[[1L]])
  python_start <- which(startsWith(parts, "Python "))
  if (length(python_start) == 0L) {
    return(list(r_packages = parts, python_modules = character()))
  }
  python_start <- python_start[[1L]]
  r_packages <- if (python_start > 1L) {
    parts[seq_len(python_start - 1L)]
  } else {
    character()
  }
  python_modules <- parts[seq.int(python_start, length(parts))]
  python_modules[[1L]] <- sub("^Python ", "", python_modules[[1L]])
  list(r_packages = r_packages, python_modules = python_modules)
}

sb_r_package_installed <- function(package) {
  length(find.package(package, quiet = TRUE)) == 1L
}

sb_dependency_status_row <- function(dependency) {
  components <- sb_dependency_components(dependency)
  missing <- components$r_packages[!vapply(
    components$r_packages,
    sb_r_package_installed,
    logical(1)
  )]
  collapse <- function(x) {
    if (length(x) == 0L) NA_character_ else paste(x, collapse = ", ")
  }
  status <- if (length(components$r_packages) == 0L &&
    length(components$python_modules) == 0L) {
    "no_optional_dependency"
  } else if (length(missing) > 0L) {
    "missing_r_packages"
  } else if (length(components$python_modules) > 0L) {
    "python_not_checked"
  } else {
    "r_packages_installed"
  }
  data.frame(
    dependency = dependency,
    r_packages = collapse(components$r_packages),
    python_modules = collapse(components$python_modules),
    missing_r_packages = collapse(missing),
    dependency_status = status,
    stringsAsFactors = FALSE
  )
}

sb_add_dependency_status <- function(registry) {
  dependencies <- unique(registry$optional_dependency)
  status <- do.call(
    rbind,
    lapply(dependencies, sb_dependency_status_row)
  )
  status <- status[
    match(registry$optional_dependency, status$dependency),
    c(
      "r_packages", "python_modules", "missing_r_packages",
      "dependency_status"
    ),
    drop = FALSE
  ]
  rownames(status) <- NULL
  cbind(registry, status)
}

# ---- legacy-name-compatibility ----

# Legacy identifiers remain accepted; newly constructed objects use canonical IDs.
sb_legacy_learner_names <- c(
  "aalen_additive_survival_fixed" = "aalen_additive_survival",
  "adaboost_binary_fixed" = "adaboost_binary",
  "ahaz_lin_ying_fixed" = "ahaz_lin_ying",
  "akritas_survival_fixed" = "akritas_survival",
  "aorsf_binary_fixed" = "aorsf_binary",
  "aorsf_survival_fixed" = "aorsf_survival",
  "bagged_knn_survival_fixed" = "bagged_knn_survival",
  "bagging_binary_fixed" = "bagging_binary",
  "bart_aft_survival_fixed" = "bart_aft_survival",
  "bart_binary_fixed" = "bart_binary",
  "bart_survival_fixed" = "bart_survival",
  "bayesglm_logistic_cauchy_fixed" = "bayesglm_logistic_cauchy",
  "bernoulli_naive_bayes_binary_fixed" = "bernoulli_naive_bayes_binary",
  "blackboost_binary_fixed" = "blackboost_binary",
  "blackboost_cox_fixed" = "blackboost_cox",
  "blackboost_gehan_aft_fixed" = "blackboost_gehan_aft",
  "blackboost_loglogistic_aft_fixed" = "blackboost_loglogistic_aft",
  "blackboost_lognormal_aft_fixed" = "blackboost_lognormal_aft",
  "blackboost_unoc_fixed" = "blackboost_unoc",
  "blackboost_weibull_aft_fixed" = "blackboost_weibull_aft",
  "brglm2_logistic_fixed" = "brglm2_logistic",
  "bst_hinge_binary_fixed" = "bst_hinge_binary",
  "bujar_linear_boost_fixed" = "bujar_linear_boost",
  "c50_binary_fixed" = "c50_binary",
  "c50_rules_binary_fixed" = "c50_rules_binary",
  "cancerclass_binary_fixed" = "cancerclass_binary",
  "catboost_aft_survival_fixed" = "catboost_aft_survival",
  "catboost_binary_fixed" = "catboost_binary",
  "catboost_cox_fixed" = "catboost_cox",
  "cforest_binary_fixed" = "cforest_binary",
  "cforest_survival_fixed" = "cforest_survival",
  "coxboost_fixed" = "coxboost",
  "coxphf_firth_fixed" = "coxphf_firth",
  "ctree_binary_fixed" = "ctree_binary",
  "ctree_survival_fixed" = "ctree_survival",
  "decision_tree_binary_fixed" = "decision_tree_binary",
  "dwd_binary_fixed" = "dwd_binary",
  "earth_binary_fixed" = "earth_binary",
  "evtree_binary_fixed" = "evtree_binary",
  "fda_binary_fixed" = "fda_binary",
  "flexsurv_gamma_fixed" = "flexsurv_gamma",
  "flexsurv_gengamma_fixed" = "flexsurv_gengamma",
  "flexsurv_gompertz_fixed" = "flexsurv_gompertz",
  "flexsurvspline_hazard_k1_fixed" = "flexsurvspline_hazard_k1",
  "flexsurvspline_normal_k1_fixed" = "flexsurvspline_normal_k1",
  "flexsurvspline_odds_k1_fixed" = "flexsurvspline_odds_k1",
  "gamboost_binary_fixed" = "gamboost_binary",
  "gamboost_cox_fixed" = "gamboost_cox",
  "gamboost_gehan_aft_fixed" = "gamboost_gehan_aft",
  "gamboost_loglogistic_aft_fixed" = "gamboost_loglogistic_aft",
  "gamboost_lognormal_aft_fixed" = "gamboost_lognormal_aft",
  "gamboost_unoc_fixed" = "gamboost_unoc",
  "gamboost_weibull_aft_fixed" = "gamboost_weibull_aft",
  "gaussian_process_binary_fixed" = "gaussian_process_binary",
  "gbm_binary_fixed" = "gbm_binary",
  "gbm_cox_fixed" = "gbm_cox",
  "glmboost_binary_fixed" = "glmboost_binary",
  "glmboost_cox_fixed" = "glmboost_cox",
  "glmboost_gehan_aft_fixed" = "glmboost_gehan_aft",
  "glmboost_loglogistic_aft_fixed" = "glmboost_loglogistic_aft",
  "glmboost_lognormal_aft_fixed" = "glmboost_lognormal_aft",
  "glmboost_unoc_fixed" = "glmboost_unoc",
  "glmboost_weibull_aft_fixed" = "glmboost_weibull_aft",
  "hda_binary_fixed" = "hda_binary",
  "kernel_naive_bayes_binary_fixed" = "kernel_naive_bayes_binary",
  "kernel_weighted_knn_binary_fixed" = "kernel_weighted_knn_binary",
  "knn_binary_fixed" = "knn_binary",
  "ktsp_binary_fixed" = "ktsp_binary",
  "liblinear_l1_l2loss_svc_binary_fixed" = "liblinear_l1_l2loss_svc_binary",
  "liblinear_l1_logistic_binary_fixed" = "liblinear_l1_logistic_binary",
  "liblinear_l2_l1loss_svc_binary_fixed" = "liblinear_l2_l1loss_svc_binary",
  "liblinear_l2_l2loss_svc_binary_fixed" = "liblinear_l2_l2loss_svc_binary",
  "liblinear_l2_logistic_binary_fixed" = "liblinear_l2_logistic_binary",
  "lightgbm_binary_dart_fixed" = "lightgbm_binary_dart",
  "lightgbm_binary_fixed" = "lightgbm_binary",
  "lightgbm_binary_goss_fixed" = "lightgbm_binary_goss",
  "lightgbm_binary_rf_fixed" = "lightgbm_binary_rf",
  "loclda_binary_fixed" = "loclda_binary",
  "logitboost_binary_fixed" = "logitboost_binary",
  "lssvm_radial_binary_fixed" = "lssvm_radial_binary",
  "lvq1_binary_fixed" = "lvq1_binary",
  "lvq21_binary_fixed" = "lvq21_binary",
  "lvq3_binary_fixed" = "lvq3_binary",
  "mda_binary_fixed" = "mda_binary",
  "mtlr_survival_fixed" = "mtlr_survival",
  "nearest_centroid_binary_fixed" = "nearest_centroid_binary",
  "neuralnet_multilayer_binary_fixed" = "neuralnet_multilayer_binary",
  "nnet_binary_fixed" = "nnet_binary",
  "olvq1_binary_fixed" = "olvq1_binary",
  "pamr_binary_fixed" = "pamr_binary",
  "pamr_survival_fixed" = "pamr_survival",
  "pda_binary_fixed" = "pda_binary",
  "plsda_binary_fixed" = "plsda_binary",
  "plsrcox_fixed" = "plsrcox",
  "plsrglm_binary_fixed" = "plsrglm_binary",
  "prop_odds_survival_fixed" = "prop_odds_survival",
  "random_ferns_binary_fixed" = "random_ferns_binary",
  "random_forest_binary_fixed" = "random_forest_binary",
  "ranger_binary_fixed" = "ranger_binary",
  "ranger_cindex_survival_fixed" = "ranger_cindex_survival",
  "ranger_extratrees_binary_fixed" = "ranger_extratrees_binary",
  "ranger_extratrees_survival_fixed" = "ranger_extratrees_survival",
  "ranger_hellinger_binary_fixed" = "ranger_hellinger_binary",
  "ranger_maxstat_survival_fixed" = "ranger_maxstat_survival",
  "ranger_survival_fixed" = "ranger_survival",
  "rda_binary_fixed" = "rda_binary",
  "rlt_reinforcement_binary_fixed" = "rlt_reinforcement_binary",
  "rlt_reinforcement_survival_fixed" = "rlt_reinforcement_survival",
  "robust_logistic_mqle_fixed" = "robust_logistic_mqle",
  "rotation_forest_binary_fixed" = "rotation_forest_binary",
  "rpart_survival_fixed" = "rpart_survival",
  "rrf_binary_fixed" = "rrf_binary",
  "rsf_survival_fixed" = "rsf_survival",
  "sda_binary_fixed" = "sda_binary",
  "slda_binary_fixed" = "slda_binary",
  "splsda_binary_fixed" = "splsda_binary",
  "stepglm_binary_fixed" = "stepglm_binary",
  "superpc_fixed" = "superpc",
  "survivalsvm_additive_fixed" = "survivalsvm_additive",
  "survivalsvm_fixed" = "survivalsvm",
  "survivalsvm_hybrid_additive_fixed" = "survivalsvm_hybrid_additive",
  "survivalsvm_hybrid_fixed" = "survivalsvm_hybrid",
  "survivalsvm_hybrid_rbf_fixed" = "survivalsvm_hybrid_rbf",
  "survivalsvm_rbf_fixed" = "survivalsvm_rbf",
  "survivalsvm_vanbelle1_additive_fixed" = "survivalsvm_vanbelle1_additive",
  "survivalsvm_vanbelle1_fixed" = "survivalsvm_vanbelle1",
  "survivalsvm_vanbelle1_rbf_fixed" = "survivalsvm_vanbelle1_rbf",
  "survivalsvm_vanbelle2_additive_fixed" = "survivalsvm_vanbelle2_additive",
  "survivalsvm_vanbelle2_fixed" = "survivalsvm_vanbelle2",
  "survivalsvm_vanbelle2_rbf_fixed" = "survivalsvm_vanbelle2_rbf",
  "survreg_exponential_fixed" = "survreg_exponential",
  "survreg_loglogistic_fixed" = "survreg_loglogistic",
  "survreg_lognormal_fixed" = "survreg_lognormal",
  "survreg_weibull_fixed" = "survreg_weibull",
  "svm_linear_binary_fixed" = "svm_linear_binary",
  "svm_polynomial_binary_fixed" = "svm_polynomial_binary",
  "svm_radial_binary_fixed" = "svm_radial_binary",
  "svm_radial_weighted_binary_fixed" = "svm_radial_weighted_binary",
  "svm_sigmoid_binary_fixed" = "svm_sigmoid_binary",
  "tree_binary_fixed" = "tree_binary",
  "wsrf_binary_fixed" = "wsrf_binary",
  "xgboost_aft_dart_fixed" = "xgboost_aft_dart",
  "xgboost_aft_fixed" = "xgboost_aft",
  "xgboost_aft_linear_fixed" = "xgboost_aft_linear",
  "xgboost_binary_dart_fixed" = "xgboost_binary_dart",
  "xgboost_binary_fixed" = "xgboost_binary",
  "xgboost_binary_linear_fixed" = "xgboost_binary_linear",
  "xgboost_cox_dart_fixed" = "xgboost_cox_dart",
  "xgboost_cox_fixed" = "xgboost_cox",
  "xgboost_cox_linear_fixed" = "xgboost_cox_linear"
)

sb_canonical_learner_id <- function(id) {
  if (!is.character(id)) return(id)
  at <- match(id, names(sb_legacy_learner_names))
  take <- !is.na(at)
  id[take] <- unname(sb_legacy_learner_names[at[take]])
  id
}

sb_canonical_name_text <- function(text) {
  if (!is.character(text) || !any(grepl("_fixed", text, fixed = TRUE))) return(text)
  for (old in names(sb_legacy_learner_names)) {
    text <- gsub(old, sb_legacy_learner_names[[old]], text, fixed = TRUE)
  }
  text
}

sb_upgrade_legacy_names <- function(object) {
  # Only software metadata is migrated; expression values, gene names, patient IDs,
  # fitted coefficients and native model objects are preserved.
  if (!is.list(object) || is.data.frame(object)) return(object)
  metadata <- c("id", "learner", "model_id", "pipeline_id", "pipeline_ids", "risk_reference",
                "risk_formula", "method", "algorithm", "source")
  children <- c("pipeline", "learner_state", "runtime", "specification",
                "provenance", "prediction_semantics", "learner", "prediction",
                "baseline", "baseline_hazard", "direct_probability", "stack", "chain_base", "inner")
  for (key in intersect(names(object), metadata)) {
    if (is.character(object[[key]])) object[[key]] <- sb_canonical_name_text(object[[key]])
  }
  for (key in intersect(names(object), children)) {
    object[[key]] <- sb_upgrade_legacy_names(object[[key]])
  }
  for (key in intersect(names(object), c("pipelines", "member_fits"))) {
    if (is.list(object[[key]]) && !is.data.frame(object[[key]])) {
      object[[key]] <- lapply(object[[key]], sb_upgrade_legacy_names)
    }
  }
  if (identical(object$selector, "ordered_chain") && is.list(object$selector_params$stages)) {
    object$selector_params$stages <- lapply(object$selector_params$stages, sb_upgrade_legacy_names)
  }
  object
}

sb_additional_survival_fixed_learner_adapter <- function(...) sb_additional_survival_learner_adapter(...)

sb_additional_survival_fixed_learner_ids <- function(...) sb_additional_survival_learner_ids(...)

sb_additional_survival_fixed_params <- function(...) sb_additional_survival_params(...)

sb_cancerclass_binary_fixed_learner_adapter <- function(...) sb_cancerclass_binary_learner_adapter(...)

sb_cancerclass_binary_fixed_params <- function(...) sb_cancerclass_binary_params(...)

sb_catboost_aft_survival_fixed_learner_adapter <- function(...) sb_catboost_aft_survival_learner_adapter(...)

sb_catboost_aft_survival_fixed_params <- function(...) sb_catboost_aft_survival_params(...)

sb_catboost_binary_fixed_learner_adapter <- function(...) sb_catboost_binary_learner_adapter(...)

sb_catboost_binary_fixed_params <- function(...) sb_catboost_binary_params(...)

sb_catboost_cox_fixed_learner_adapter <- function(...) sb_catboost_cox_learner_adapter(...)

sb_catboost_cox_fixed_params <- function(...) sb_catboost_cox_params(...)

sb_classical_binary_fixed_dependency <- function(...) sb_classical_binary_dependency(...)

sb_classical_binary_fixed_learner_adapter <- function(...) sb_classical_binary_learner_adapter(...)

sb_classical_binary_fixed_learner_ids <- function(...) sb_classical_binary_learner_ids(...)

sb_classical_binary_fixed_params <- function(...) sb_classical_binary_params(...)

sb_classical_survival_fixed_adapter <- function(...) sb_classical_survival_adapter(...)

sb_classical_survival_fixed_learner_ids <- function(...) sb_classical_survival_learner_ids(...)

sb_coxboost_fixed_learner_adapter <- function(...) sb_coxboost_learner_adapter(...)

sb_coxboost_fixed_params <- function(...) sb_coxboost_params(...)

sb_extract_additional_survival_fixed_learner_state <- function(...) sb_extract_additional_survival_learner_state(...)

sb_extract_classical_binary_fixed_learner_state <- function(...) sb_extract_classical_binary_learner_state(...)

sb_extract_xgboost_fixed_learner_state <- function(...) sb_extract_xgboost_learner_state(...)

sb_fit_additional_survival_fixed <- function(...) sb_fit_additional_survival(...)

sb_fit_bart_binary_fixed <- function(...) sb_fit_bart_binary(...)

sb_fit_bart_survival_fixed <- function(...) sb_fit_bart_survival(...)

sb_fit_bayesglm_logistic_cauchy_fixed <- function(...) sb_fit_bayesglm_logistic_cauchy(...)

sb_fit_brglm2_logistic_fixed <- function(...) sb_fit_brglm2_logistic(...)

sb_fit_bst_hinge_binary_fixed <- function(...) sb_fit_bst_hinge_binary(...)

sb_fit_c50_binary_fixed <- function(...) sb_fit_c50_binary(...)

sb_fit_cancerclass_binary_fixed <- function(...) sb_fit_cancerclass_binary(...)

sb_fit_catboost_aft_survival_fixed <- function(...) sb_fit_catboost_aft_survival(...)

sb_fit_catboost_binary_fixed <- function(...) sb_fit_catboost_binary(...)

sb_fit_catboost_cox_fixed <- function(...) sb_fit_catboost_cox(...)

sb_fit_classical_binary_boosting_fixed <- function(...) sb_fit_classical_binary_boosting(...)

sb_fit_classical_binary_fixed <- function(...) sb_fit_classical_binary(...)

sb_fit_coxboost_fixed <- function(...) sb_fit_coxboost(...)

sb_fit_dwd_binary_fixed <- function(...) sb_fit_dwd_binary(...)

sb_fit_earth_binary_fixed <- function(...) sb_fit_earth_binary(...)

sb_fit_flexsurv_parametric_fixed <- function(...) sb_fit_flexsurv_parametric(...)

sb_fit_flexsurvspline_fixed <- function(...) sb_fit_flexsurvspline(...)

sb_fit_gaussian_process_binary_fixed <- function(...) sb_fit_gaussian_process_binary(...)

sb_fit_gbm_cox_fixed <- function(...) sb_fit_gbm_cox(...)

sb_fit_glmboost_binary_fixed <- function(...) sb_fit_glmboost_binary(...)

sb_fit_hda_binary_fixed <- function(...) sb_fit_hda_binary(...)

sb_fit_kernel_weighted_knn_binary_fixed <- function(...) sb_fit_kernel_weighted_knn_binary(...)

sb_fit_ktsp_binary_fixed <- function(...) sb_fit_ktsp_binary(...)

sb_fit_liblinear_binary_fixed <- function(...) sb_fit_liblinear_binary(...)

sb_fit_lightgbm_binary_fixed <- function(...) sb_fit_lightgbm_binary(...)

sb_fit_lssvm_radial_binary_fixed <- function(...) sb_fit_lssvm_radial_binary(...)

sb_fit_lvq_binary_fixed <- function(...) sb_fit_lvq_binary(...)

sb_fit_mboost_aft_fixed <- function(...) sb_fit_mboost_aft(...)

sb_fit_mboost_cox_fixed <- function(...) sb_fit_mboost_cox(...)

sb_fit_mboost_ranking_fixed <- function(...) sb_fit_mboost_ranking(...)

sb_fit_neuralnet_multilayer_binary_fixed <- function(...) sb_fit_neuralnet_multilayer_binary(...)

sb_fit_pamr_binary_fixed <- function(...) sb_fit_pamr_binary(...)

sb_fit_pda_binary_fixed <- function(...) sb_fit_pda_binary(...)

sb_fit_plsrcox_fixed <- function(...) sb_fit_plsrcox(...)

sb_fit_plsrglm_binary_fixed <- function(...) sb_fit_plsrglm_binary(...)

sb_fit_random_ferns_binary_fixed <- function(...) sb_fit_random_ferns_binary(...)

sb_fit_ranger_binary_fixed <- function(...) sb_fit_ranger_binary(...)

sb_fit_ranger_survival_fixed <- function(...) sb_fit_ranger_survival(...)

sb_fit_robust_logistic_mqle_fixed <- function(...) sb_fit_robust_logistic_mqle(...)

sb_fit_rotation_forest_binary_fixed <- function(...) sb_fit_rotation_forest_binary(...)

sb_fit_stepglm_binary_fixed <- function(...) sb_fit_stepglm_binary(...)

sb_fit_superpc_fixed <- function(...) sb_fit_superpc(...)

sb_fit_survivalsvm_fixed <- function(...) sb_fit_survivalsvm(...)

sb_fit_survreg_fixed <- function(...) sb_fit_survreg(...)

sb_fit_timereg_survival_fixed <- function(...) sb_fit_timereg_survival(...)

sb_fit_xgboost_fixed <- function(...) sb_fit_xgboost(...)

sb_flexsurv_parametric_fixed_learner_adapter <- function(...) sb_flexsurv_parametric_learner_adapter(...)

sb_flexsurv_parametric_fixed_learner_ids <- function(...) sb_flexsurv_parametric_learner_ids(...)

sb_flexsurv_parametric_fixed_params <- function(...) sb_flexsurv_parametric_params(...)

sb_flexsurvspline_fixed_learner_adapter <- function(...) sb_flexsurvspline_learner_adapter(...)

sb_flexsurvspline_fixed_learner_ids <- function(...) sb_flexsurvspline_learner_ids(...)

sb_flexsurvspline_fixed_params <- function(...) sb_flexsurvspline_params(...)

sb_gbm_cox_fixed_learner_adapter <- function(...) sb_gbm_cox_learner_adapter(...)

sb_gbm_cox_fixed_params <- function(...) sb_gbm_cox_params(...)

sb_liblinear_binary_fixed_params <- function(...) sb_liblinear_binary_params(...)

sb_lightgbm_binary_fixed_learner_adapter <- function(...) sb_lightgbm_binary_learner_adapter(...)

sb_lightgbm_binary_fixed_params <- function(...) sb_lightgbm_binary_params(...)

sb_lvq_binary_fixed_learner_ids <- function(...) sb_lvq_binary_learner_ids(...)

sb_lvq_binary_fixed_params <- function(...) sb_lvq_binary_params(...)

sb_mboost_aft_fixed_learner_adapter <- function(...) sb_mboost_aft_learner_adapter(...)

sb_mboost_aft_fixed_learner_ids <- function(...) sb_mboost_aft_learner_ids(...)

sb_mboost_aft_fixed_params <- function(...) sb_mboost_aft_params(...)

sb_mboost_cox_fixed_learner_adapter <- function(...) sb_mboost_cox_learner_adapter(...)

sb_mboost_cox_fixed_learner_ids <- function(...) sb_mboost_cox_learner_ids(...)

sb_mboost_cox_fixed_params <- function(...) sb_mboost_cox_params(...)

sb_mboost_ranking_fixed_learner_adapter <- function(...) sb_mboost_ranking_learner_adapter(...)

sb_mboost_ranking_fixed_learner_ids <- function(...) sb_mboost_ranking_learner_ids(...)

sb_mboost_ranking_fixed_params <- function(...) sb_mboost_ranking_params(...)

sb_new_additional_survival_fixed_model <- function(...) sb_new_additional_survival_model(...)

sb_new_classical_binary_fixed_model <- function(...) sb_new_classical_binary_model(...)

sb_new_classical_survival_fixed_model <- function(...) sb_new_classical_survival_model(...)

sb_new_flexsurv_parametric_fixed_model <- function(...) sb_new_flexsurv_parametric_model(...)

sb_new_flexsurvspline_fixed_model <- function(...) sb_new_flexsurvspline_model(...)

sb_new_mboost_aft_fixed_model <- function(...) sb_new_mboost_aft_model(...)

sb_new_mboost_cox_fixed_model <- function(...) sb_new_mboost_cox_model(...)

sb_new_mboost_ranking_fixed_model <- function(...) sb_new_mboost_ranking_model(...)

sb_plsrcox_fixed_learner_adapter <- function(...) sb_plsrcox_learner_adapter(...)

sb_plsrcox_fixed_params <- function(...) sb_plsrcox_params(...)

sb_predict_additional_survival_fixed <- function(...) sb_predict_additional_survival(...)

sb_predict_additional_survival_fixed_learner_state <- function(...) sb_predict_additional_survival_learner_state(...)

sb_predict_additional_survival_fixed_survival <- function(...) sb_predict_additional_survival_survival(...)

sb_predict_bart_binary_fixed <- function(...) sb_predict_bart_binary(...)

sb_predict_c50_binary_fixed <- function(...) sb_predict_c50_binary(...)

sb_predict_cancerclass_binary_fixed <- function(...) sb_predict_cancerclass_binary(...)

sb_predict_catboost_aft_survival_fixed <- function(...) sb_predict_catboost_aft_survival(...)

sb_predict_catboost_aft_survival_fixed_survival <- function(...) sb_predict_catboost_aft_survival_survival(...)

sb_predict_catboost_binary_fixed <- function(...) sb_predict_catboost_binary(...)

sb_predict_catboost_cox_fixed <- function(...) sb_predict_catboost_cox(...)

sb_predict_catboost_cox_fixed_survival <- function(...) sb_predict_catboost_cox_survival(...)

sb_predict_cforest_binary_fixed <- function(...) sb_predict_cforest_binary(...)

sb_predict_classical_binary_boosting_fixed <- function(...) sb_predict_classical_binary_boosting(...)

sb_predict_classical_binary_fixed <- function(...) sb_predict_classical_binary(...)

sb_predict_classical_binary_fixed_learner_state <- function(...) sb_predict_classical_binary_learner_state(...)

sb_predict_classical_survival_fixed <- function(...) sb_predict_classical_survival(...)

sb_predict_classical_survival_fixed_survival <- function(...) sb_predict_classical_survival_survival(...)

sb_predict_earth_binary_fixed <- function(...) sb_predict_earth_binary(...)

sb_predict_flexsurv_parametric_fixed <- function(...) sb_predict_flexsurv_parametric(...)

sb_predict_flexsurv_parametric_fixed_survival <- function(...) sb_predict_flexsurv_parametric_survival(...)

sb_predict_flexsurvspline_fixed <- function(...) sb_predict_flexsurvspline(...)

sb_predict_flexsurvspline_fixed_survival <- function(...) sb_predict_flexsurvspline_survival(...)

sb_predict_gaussian_process_binary_fixed <- function(...) sb_predict_gaussian_process_binary(...)

sb_predict_gbm_cox_fixed <- function(...) sb_predict_gbm_cox(...)

sb_predict_gbm_cox_fixed_survival <- function(...) sb_predict_gbm_cox_survival(...)

sb_predict_hda_binary_fixed <- function(...) sb_predict_hda_binary(...)

sb_predict_kernel_weighted_knn_binary_fixed <- function(...) sb_predict_kernel_weighted_knn_binary(...)

sb_predict_klar_rda_binary_fixed <- function(...) sb_predict_klar_rda_binary(...)

sb_predict_knn_binary_fixed <- function(...) sb_predict_knn_binary(...)

sb_predict_liblinear_binary_fixed <- function(...) sb_predict_liblinear_binary(...)

sb_predict_lightgbm_binary_fixed <- function(...) sb_predict_lightgbm_binary(...)

sb_predict_loclda_binary_fixed <- function(...) sb_predict_loclda_binary(...)

sb_predict_lssvm_radial_binary_fixed <- function(...) sb_predict_lssvm_radial_binary(...)

sb_predict_lvq_binary_fixed <- function(...) sb_predict_lvq_binary(...)

sb_predict_mboost_aft_fixed <- function(...) sb_predict_mboost_aft(...)

sb_predict_mboost_cox_fixed <- function(...) sb_predict_mboost_cox(...)

sb_predict_mboost_ranking_fixed <- function(...) sb_predict_mboost_ranking(...)

sb_predict_neuralnet_multilayer_binary_fixed <- function(...) sb_predict_neuralnet_multilayer_binary(...)

sb_predict_pamr_binary_fixed <- function(...) sb_predict_pamr_binary(...)

sb_predict_random_ferns_binary_fixed <- function(...) sb_predict_random_ferns_binary(...)

sb_predict_ranger_binary_fixed <- function(...) sb_predict_ranger_binary(...)

sb_predict_ranger_survival_fixed <- function(...) sb_predict_ranger_survival(...)

sb_predict_ranger_survival_fixed_survival <- function(...) sb_predict_ranger_survival_survival(...)

sb_predict_robust_logistic_mqle_fixed <- function(...) sb_predict_robust_logistic_mqle(...)

sb_predict_sda_binary_fixed <- function(...) sb_predict_sda_binary(...)

sb_predict_slda_binary_fixed <- function(...) sb_predict_slda_binary(...)

sb_predict_survivalsvm_kernel_fixed <- function(...) sb_predict_survivalsvm_kernel(...)

sb_predict_timereg_survival_fixed <- function(...) sb_predict_timereg_survival(...)

sb_predict_timereg_survival_fixed_survival <- function(...) sb_predict_timereg_survival_survival(...)

sb_predict_tree_binary_fixed <- function(...) sb_predict_tree_binary(...)

sb_predict_xgboost_fixed_learner_state <- function(...) sb_predict_xgboost_learner_state(...)

sb_ranger_binary_fixed_learner_adapter <- function(...) sb_ranger_binary_learner_adapter(...)

sb_ranger_binary_fixed_params <- function(...) sb_ranger_binary_params(...)

sb_ranger_fixed_learner_ids <- function(...) sb_ranger_learner_ids(...)

sb_ranger_survival_fixed_learner_adapter <- function(...) sb_ranger_survival_learner_adapter(...)

sb_ranger_survival_fixed_params <- function(...) sb_ranger_survival_params(...)

sb_require_classical_binary_fixed_dependency <- function(...) sb_require_classical_binary_dependency(...)

sb_restore_additional_survival_fixed_model <- function(...) sb_restore_additional_survival_model(...)

sb_superpc_fixed_learner_adapter <- function(...) sb_superpc_learner_adapter(...)

sb_superpc_fixed_params <- function(...) sb_superpc_params(...)

sb_survivalsvm_fixed_learner_adapter <- function(...) sb_survivalsvm_learner_adapter(...)

sb_survivalsvm_fixed_params <- function(...) sb_survivalsvm_params(...)

sb_survreg_fixed_learner_ids <- function(...) sb_survreg_learner_ids(...)

sb_validate_classical_binary_fixed_inputs <- function(...) sb_validate_classical_binary_inputs(...)

sb_validate_xgboost_fixed_params <- function(...) sb_validate_xgboost_params(...)

sb_xgboost_fixed_booster <- function(...) sb_xgboost_booster(...)

sb_xgboost_fixed_learner_adapter <- function(...) sb_xgboost_learner_adapter(...)

sb_xgboost_fixed_learner_ids <- function(...) sb_xgboost_learner_ids(...)

sb_xgboost_fixed_params <- function(...) sb_xgboost_params(...)

# ---- stormbreaker-package ----

#' @rdname stormbreaker-package
"_PACKAGE"
