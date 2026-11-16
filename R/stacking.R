#' Stacked models
#' @name sb_stack
#' @description Specify, evaluate and fit predeclared stacks with nested resampling.
#' @rawRd \arguments{
#' \item{pipelines}{
#' An ordered list containing at least two distinct pipeline
#' specifications returned by \code{\link[=sb_pipeline]{sb_pipeline()}}.
#' }
#' \item{inner_folds}{
#' Number of event- or outcome-stratified folds used to
#' create member out-of-fold predictions and select the ridge penalty inside
#' the current analysis data. Must be at least three.
#' }
#' \item{inner_seed}{
#' Positive integer seed for the inner fold assignment.
#' }
#' \item{task}{
#' A task returned by \code{\link[=sb_surv_task]{sb_surv_task()}} or \code{\link[=sb_binary_task]{sb_binary_task()}}.
#' }
#' \item{stack}{
#' A stack specification returned by \code{\link[=sb_stack]{sb_stack()}}.
#' }
#' \item{folds}{
#' Either one outer-fold count or an exact patient-ID-named fold
#' assignment accepted by \code{\link[=sb_run_cv]{sb_run_cv()}}.
#' }
#' \item{seed}{
#' Positive integer used when \code{folds} is a scalar count.
#' }
#' \item{fit}{
#' A successful result from \code{\link[=sb_fit_stack_final]{sb_fit_stack_final()}}.
#' }
#' \item{new_data}{
#' A non-empty data frame containing the declared ID and every
#' numeric feature required by the saved member fits. Outcome columns are not
#' required and, when present, are ignored.
#' }
#' \item{id_col}{
#' Name of the unique patient identifier in \code{new_data}.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_stack()}}{
#' \code{sb_stack()} records an ordered set of existing STORMBREAKER pipelines for
#' stacked generalization. It does not fit models or inspect outcomes. Binary
#' stacks use a ridge logistic meta-learner; survival stacks use a ridge Cox
#' meta-learner. In both cases the ridge penalty is selected from member
#' out-of-fold predictions inside the current analysis data.
#' 
#' \strong{Returns:} A \code{stormbreaker_stack} specification. No model is fitted.
#' 
#' \strong{Details:} Members are never searched, dropped, replaced, or reordered. At least two
#' distinct pipelines with one common task type are required.
#' }
#' 
#' \subsection{\code{sb_run_stack_cv()}}{
#' \code{sb_run_stack_cv()} compares every stack member on one shared outer fold
#' assignment and evaluates their ridge stack on those same held-out patients.
#' Inside each outer analysis fold, members are rerun under an inner
#' event- or outcome-stratified assignment. Only those inner out-of-fold
#' predictions and analysis outcomes may train or tune the meta-learner. Outer
#' assessment outcomes are used only after stacked predictions are fixed.
#' 
#' \strong{Returns:} A \code{stormbreaker_stack_cv_result} with member results, stacked OOF
#' predictions and metrics, per-outer-fold inner OOF/meta audits, and one
#' descriptive member-plus-stack comparison table.
#' 
#' \strong{Details:} A member with incomplete inner or outer predictions causes the corresponding
#' stack fold to fail. The function never removes or replaces a member. The
#' returned comparison is descriptive and does not select a winner.
#' }
#' 
#' \subsection{\code{sb_fit_stack_final()}}{
#' \code{sb_fit_stack_final()} first creates member out-of-fold predictions within
#' the complete development task and fits the stack's fixed ridge meta-model to
#' those predictions. It then fits every predeclared member pipeline once on all
#' development patients. The external cohort is not read and no member is
#' selected, dropped, or replaced.
#' 
#' \strong{Returns:} A \code{stormbreaker_stack_final_fit} with stage-aware success or failure,
#' inner resampling audit, meta fit, ordered member fits, required feature
#' union, prediction semantics, and provenance.
#' 
#' \strong{Details:} The returned fit stores portable member fits and portable meta coefficients.
#' For a survival task with a declared horizon, its Breslow baseline is learned
#' from development OOF stacked risks rather than external outcomes.
#' }
#' 
#' \subsection{\code{sb_predict_stack_final()}}{
#' \code{sb_predict_stack_final()} asks every saved member fit for its outcome-blind
#' prediction, preserves the declared member order, and applies the saved ridge
#' meta coefficients. It does not refit a member, the meta-model, or a survival
#' baseline. A failed member prediction makes the complete stack prediction
#' fail; no member is removed or replaced.
#' 
#' \strong{Returns:} A \code{stormbreaker_final_prediction} with stacked patient predictions,
#' the aligned member-prediction matrix, and stage-aware failure information.
#' 
#' \strong{Details:} The returned object uses the same final-prediction contract as
#' \code{\link[=sb_predict_final]{sb_predict_final()}} and can be passed directly to \code{\link[=sb_evaluate_external]{sb_evaluate_external()}}.
#' }
#' }
#' @rawRd \examples{
#' # sb_stack
#' binary_stack <- sb_stack(list(
#'   sb_pipeline(learner = "logistic"),
#'   sb_pipeline(learner = "gaussian_naive_bayes")
#' ), inner_folds = 3, inner_seed = 2026)
#' binary_stack
#' 
#' # sb_run_stack_cv
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   d <- sb_example_binary(n = 48, seed = 1)
#'   task <- sb_binary_task(
#'     d, "patient_id", "outcome", c("gene_a", "gene_b")
#'   )
#'   stack <- sb_stack(list(
#'     sb_pipeline(learner = "logistic"),
#'     sb_pipeline(learner = "gaussian_naive_bayes")
#'   ), inner_folds = 3)
#'   sb_run_stack_cv(task, stack, folds = 3, seed = 2026)
#' }
#' 
#' # sb_fit_stack_final
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   d <- sb_example_binary(n = 48, seed = 1)
#'   task <- sb_binary_task(
#'     d, "patient_id", "outcome", c("gene_a", "gene_b")
#'   )
#'   stack <- sb_stack(list(
#'     sb_pipeline(learner = "logistic"),
#'     sb_pipeline(learner = "gaussian_naive_bayes")
#'   ), inner_folds = 3)
#'   sb_fit_stack_final(task, stack)
#' }
#' 
#' # sb_predict_stack_final
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   d <- sb_example_binary(n = 48, seed = 1)
#'   task <- sb_binary_task(
#'     d, "patient_id", "outcome", c("gene_a", "gene_b")
#'   )
#'   stack <- sb_stack(list(
#'     sb_pipeline(learner = "logistic"),
#'     sb_pipeline(learner = "gaussian_naive_bayes")
#'   ), inner_folds = 3)
#'   fit <- sb_fit_stack_final(task, stack)
#'   sb_predict_stack_final(
#'     fit,
#'     d[1:3, c("patient_id", "gene_a", "gene_b")],
#'     "patient_id"
#'   )
#' }
#' }
NULL

# ---- stacking ----

# ---- stacking-cv ----

sb_stack_outer_member_matrix <- function(comparison, patient_ids, stack) {
  task_type <- stack$task_type
  value_name <- if (identical(task_type, "binary")) {
    "probability"
  } else {
    "risk"
  }
  columns <- lapply(stack$pipeline_ids, function(pipeline_id) {
    result <- comparison$results[[pipeline_id]]
    predictions <- if (is.list(result)) result$predictions else NULL
    valid_frame <- is.data.frame(predictions) &&
      all(c("patient_id", "status", value_name) %in% names(predictions))
    if (!valid_frame) {
      sb_incomplete_stack_member_oof(pipeline_id)
    }
    prediction_ids <- as.character(predictions$patient_id)
    positions <- match(patient_ids, prediction_ids)
    status <- predictions$status[positions]
    value <- predictions[[value_name]][positions]
    valid <- !anyNA(positions) && !anyDuplicated(prediction_ids) &&
      is.character(status) && !anyNA(status) && all(status == "success") &&
      is.numeric(value) && !anyNA(value) && all(is.finite(value))
    if (identical(task_type, "binary")) {
      valid <- valid && all(value >= 0 & value <= 1)
    }
    if (!valid) {
      sb_incomplete_stack_member_oof(pipeline_id)
    }
    unname(value)
  })
  names(columns) <- stack$pipeline_ids
  data.frame(
    columns,
    row.names = patient_ids,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

sb_stack_fold_prediction_rows <- function(task, stack, split, values = NULL) {
  patient_ids <- as.character(task$data[[task$id_col]][split$assessment_rows])
  success <- !is.null(values)
  if (inherits(task, "stormbreaker_binary_task")) {
    return(data.frame(
      patient_id = patient_ids,
      fold = split$fold,
      outcome = as.integer(task$data[[task$outcome_col]][split$assessment_rows]),
      positive_class = task$positive_class,
      probability = if (success) values$probability else NA_real_,
      status = if (success) "success" else "failed",
      pipeline_id = stack$id,
      learner = stack$meta_learner,
      stringsAsFactors = FALSE
    ))
  }

  rows <- data.frame(
    patient_id = patient_ids,
    fold = split$fold,
    time = task$data[[task$time_col]][split$assessment_rows],
    event = as.integer(task$data[[task$event_col]][split$assessment_rows]),
    time_unit = task$time_unit,
    risk = if (success) values$risk else NA_real_,
    status = if (success) "success" else "failed",
    pipeline_id = stack$id,
    learner = stack$meta_learner,
    stringsAsFactors = FALSE
  )
  if ("prediction_horizon" %in% names(task)) {
    rows$prediction_horizon <- task$prediction_horizon
    rows$survival_probability <- if (success) {
      values$survival_probability
    } else {
      NA_real_
    }
    rows$survival_probability_status <- if (success) "success" else "failed"
  }
  rows
}

sb_run_one_stack_fold <- function(task, stack, outer_members, split,
                                  outer_seed) {
  inner_comparison <- NULL
  member_matrix <- NULL
  meta <- NULL
  tryCatch(
    {
      analysis_task <- sb_subset_task(task, split$analysis_rows)
      inner_comparison <- sb_compare_cv(
        analysis_task,
        stack$pipelines,
        folds = stack$inner_folds,
        seed = stack$inner_seed
      )
      member_matrix <- sb_stack_member_matrix(
        inner_comparison,
        analysis_task,
        stack
      )
      meta <- sb_fit_stack_meta(analysis_task, stack, member_matrix)
      assessment_matrix <- sb_stack_outer_member_matrix(
        outer_members,
        split$assessment_ids,
        stack
      )
      values <- sb_predict_stack_meta(meta$state, assessment_matrix)
      predictions <- sb_stack_fold_prediction_rows(
        task,
        stack,
        split,
        values
      )
      run <- data.frame(
        fold = split$fold,
        seed = outer_seed,
        stack_id = stack$id,
        meta_learner = stack$meta_learner,
        n_members = length(stack$pipelines),
        n_analysis = length(split$analysis_rows),
        n_assessment = length(split$assessment_rows),
        inner_folds = stack$inner_folds,
        inner_seed = stack$inner_seed,
        meta_lambda = meta$tuning$lambda,
        status = "success",
        error_class = NA_character_,
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
      if (inherits(task, "stormbreaker_binary_task")) {
        run$n_analysis_positive <- sum(
          task$data[[task$outcome_col]][split$analysis_rows]
        )
        run <- run[c(
          "fold", "seed", "stack_id", "meta_learner", "n_members",
          "n_analysis", "n_analysis_positive", "n_assessment",
          "inner_folds", "inner_seed", "meta_lambda", "status",
          "error_class", "error_message"
        )]
      }
      list(
        predictions = predictions,
        run = run,
        audit = list(
          fold = split$fold,
          status = "success",
          comparison = inner_comparison,
          member_matrix = member_matrix,
          meta = meta,
          error_class = NA_character_,
          error_message = NA_character_
        )
      )
    },
    error = function(error) {
      predictions <- sb_stack_fold_prediction_rows(task, stack, split)
      run <- data.frame(
        fold = split$fold,
        seed = outer_seed,
        stack_id = stack$id,
        meta_learner = stack$meta_learner,
        n_members = length(stack$pipelines),
        n_analysis = length(split$analysis_rows),
        n_assessment = length(split$assessment_rows),
        inner_folds = stack$inner_folds,
        inner_seed = stack$inner_seed,
        meta_lambda = NA_real_,
        status = "failed",
        error_class = class(error)[[1L]],
        error_message = conditionMessage(error),
        stringsAsFactors = FALSE
      )
      if (inherits(task, "stormbreaker_binary_task")) {
        run$n_analysis_positive <- sum(
          task$data[[task$outcome_col]][split$analysis_rows]
        )
        run <- run[c(
          "fold", "seed", "stack_id", "meta_learner", "n_members",
          "n_analysis", "n_analysis_positive", "n_assessment",
          "inner_folds", "inner_seed", "meta_lambda", "status",
          "error_class", "error_message"
        )]
      }
      list(
        predictions = predictions,
        run = run,
        audit = list(
          fold = split$fold,
          status = "failed",
          comparison = inner_comparison,
          member_matrix = member_matrix,
          meta = meta,
          error_class = class(error)[[1L]],
          error_message = conditionMessage(error)
        )
      )
    }
  )
}

sb_stack_binary_metrics <- function(task, predictions, runs) {
  brier_reference <- data.frame(
    fold = runs$fold,
    analysis_prevalence = runs$n_analysis_positive / runs$n_analysis
  )
  moderate_calibration <- tryCatch(
    sb_binary_moderate_calibration_metric(
      predictions,
      scope = "pooled_patient_level_out_of_fold"
    ),
    error = function(error) {
      sb_binary_moderate_calibration_error_metric(
        error,
        scope = "pooled_patient_level_out_of_fold",
        n = nrow(predictions)
      )
    }
  )
  metrics <- list(
    roc_auc = sb_binary_auc_metric(predictions),
    average_precision = sb_binary_average_precision_metric(predictions),
    brier_score = sb_binary_brier_metric(predictions, brier_reference),
    log_loss = sb_binary_log_loss_metric(predictions, brier_reference),
    calibration = sb_binary_calibration_metric(predictions),
    moderate_calibration = moderate_calibration
  )
  if ("decision_thresholds" %in% names(task)) {
    metrics$decision_curve <- sb_binary_decision_curve_metric(
      predictions,
      task$decision_thresholds,
      scope = "pooled_patient_level_out_of_fold"
    )
  }
  metrics
}

sb_stack_survival_metrics <- function(task, predictions, splits) {
  metrics <- list(harrell_c = sb_harrell_c_metric(predictions))
  if (!"prediction_horizon" %in% names(task)) {
    return(metrics)
  }
  analysis_inputs <- lapply(splits, function(split) {
    list(
      fold = split$fold,
      time = task$data[[task$time_col]][split$analysis_rows],
      event = as.integer(task$data[[task$event_col]][split$analysis_rows])
    )
  })
  horizon <- task$prediction_horizon
  metrics$ipcw_brier <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    horizon
  )
  metrics$time_auc <- sb_survival_time_auc_metric(
    predictions,
    analysis_inputs,
    horizon,
    metrics$ipcw_brier$censoring_curves
  )
  metrics$uno_c <- sb_survival_uno_c_metric(
    predictions,
    analysis_inputs,
    horizon,
    metrics$ipcw_brier$censoring_curves
  )
  calibration_status <- ifelse(
    predictions$status == "success" &
      predictions$survival_probability_status == "success",
    "success",
    "failed"
  )
  metrics$calibration <- tryCatch(
    sb_survival_calibration_metric(
      patient_id = predictions$patient_id,
      time = predictions$time,
      event = predictions$event,
      survival_probability = predictions$survival_probability,
      prediction_status = calibration_status,
      prediction_horizon = horizon,
      time_unit = task$time_unit,
      scope = "pooled_patient_level_out_of_fold",
      fold = predictions$fold
    ),
    error = function(error) {
      sb_survival_calibration_error_metric(
        error = error,
        prediction_horizon = horizon,
        time_unit = task$time_unit,
        scope = "pooled_patient_level_out_of_fold",
        n = nrow(predictions),
        fold = predictions$fold
      )
    }
  )
  metrics$weak_calibration <- tryCatch(
    sb_survival_weak_calibration_metric(
      patient_id = predictions$patient_id,
      time = predictions$time,
      event = predictions$event,
      survival_probability = predictions$survival_probability,
      prediction_status = calibration_status,
      prediction_horizon = horizon,
      time_unit = task$time_unit,
      scope = "pooled_patient_level_out_of_fold",
      fold = predictions$fold
    ),
    error = function(error) {
      sb_survival_weak_calibration_error_metric(
        error = error,
        prediction_horizon = horizon,
        time_unit = task$time_unit,
        scope = "pooled_patient_level_out_of_fold",
        n = nrow(predictions),
        fold = predictions$fold
      )
    }
  )
  if ("decision_thresholds" %in% names(task)) {
    metrics$decision_curve <- sb_survival_decision_curve_metric(
      predictions,
      task$decision_thresholds,
      horizon,
      task$time_unit,
      scope = "pooled_patient_level_out_of_fold"
    )
  }
  metrics
}

sb_stack_comparison_row <- function(result) {
  n_folds <- as.integer(nrow(result$runs))
  n_successful <- as.integer(sum(result$runs$status == "success"))
  prediction_variation <- sb_prediction_variation(
    result$predictions, result$task_type
  )$status
  if (identical(result$task_type, "binary")) {
    metrics <- result$metrics
    return(data.frame(
      pipeline_id = result$stack$id,
      selector = "stack_members",
      learner = result$stack$meta_learner,
      status = result$status,
      complete_oof = result$complete_oof,
      prediction_variation = prediction_variation,
      roc_auc_status = metrics$roc_auc$status,
      roc_auc = metrics$roc_auc$estimate,
      comparable_pairs = metrics$roc_auc$comparable_pairs,
      average_precision_status = metrics$average_precision$status,
      average_precision = metrics$average_precision$estimate,
      brier_status = metrics$brier_score$status,
      brier_score = metrics$brier_score$estimate,
      n_scored = metrics$brier_score$n_scored,
      reference_brier_status = metrics$brier_score$reference$status,
      reference_brier = metrics$brier_score$reference$estimate,
      brier_skill_status = metrics$brier_score$skill_score$status,
      brier_skill_score = metrics$brier_score$skill_score$estimate,
      log_loss_status = metrics$log_loss$status,
      log_loss = metrics$log_loss$estimate,
      reference_log_loss_status = metrics$log_loss$reference$status,
      reference_log_loss = metrics$log_loss$reference$estimate,
      log_loss_skill_status = metrics$log_loss$skill_score$status,
      log_loss_skill_score = metrics$log_loss$skill_score$estimate,
      calibration_status = metrics$calibration$status,
      calibration_intercept = metrics$calibration$calibration_intercept,
      calibration_slope = metrics$calibration$calibration_slope,
      moderate_calibration_status = metrics$moderate_calibration$status,
      moderate_calibration_ici = metrics$moderate_calibration$ici,
      moderate_calibration_e50 = metrics$moderate_calibration$e50,
      moderate_calibration_e90 = metrics$moderate_calibration$e90,
      selection_stability_status = "not_applicable",
      mean_pairwise_jaccard = NA_real_,
      n_successful_selection_folds = 0L,
      n_folds = n_folds,
      n_successful_folds = n_successful,
      glmnet_active_set_stability_status = "not_applicable",
      glmnet_active_set_mean_pairwise_jaccard = NA_real_,
      n_successful_glmnet_active_set_folds = 0L,
      stringsAsFactors = FALSE
    ))
  }

  metrics <- result$metrics
  row <- data.frame(
    pipeline_id = result$stack$id,
    selector = "stack_members",
    learner = result$stack$meta_learner,
    status = result$status,
    complete_oof = result$complete_oof,
    prediction_variation = prediction_variation,
    harrell_c_status = metrics$harrell_c$status,
    harrell_c = metrics$harrell_c$estimate,
    comparable_pairs = metrics$harrell_c$comparable_pairs,
    selection_stability_status = "not_applicable",
    mean_pairwise_jaccard = NA_real_,
    n_successful_selection_folds = 0L,
    n_folds = n_folds,
    n_successful_folds = n_successful,
    glmnet_active_set_stability_status = "not_applicable",
    glmnet_active_set_mean_pairwise_jaccard = NA_real_,
    n_successful_glmnet_active_set_folds = 0L,
    stringsAsFactors = FALSE
  )
  if (!"ipcw_brier" %in% names(metrics)) {
    return(row)
  }
  data.frame(
    pipeline_id = row$pipeline_id,
    selector = row$selector,
    learner = row$learner,
    status = row$status,
    complete_oof = row$complete_oof,
    prediction_variation = row$prediction_variation,
    complete_survival_probability_oof =
      result$complete_survival_probability_oof,
    prediction_horizon = metrics$ipcw_brier$prediction_horizon,
    harrell_c_status = row$harrell_c_status,
    harrell_c = row$harrell_c,
    comparable_pairs = row$comparable_pairs,
    ipcw_brier_status = metrics$ipcw_brier$status,
    ipcw_brier = metrics$ipcw_brier$estimate,
    reference_ipcw_brier_status = metrics$ipcw_brier$reference$status,
    reference_ipcw_brier = metrics$ipcw_brier$reference$estimate,
    ipcw_brier_skill_status = metrics$ipcw_brier$skill_score$status,
    ipcw_brier_skill_score = metrics$ipcw_brier$skill_score$estimate,
    time_auc_status = metrics$time_auc$status,
    time_auc = metrics$time_auc$estimate,
    time_auc_pair_denominator = metrics$time_auc$pair_denominator,
    uno_c_status = metrics$uno_c$status,
    uno_c = metrics$uno_c$estimate,
    uno_c_weighted_pair_denominator = metrics$uno_c$weighted_pair_denominator,
    calibration_status = metrics$calibration$status,
    calibration_ici = metrics$calibration$ici,
    calibration_e50 = metrics$calibration$e50,
    calibration_e90 = metrics$calibration$e90,
    mean_calibration_status = metrics$weak_calibration$mean_calibration$status,
    mean_calibration_ratio = metrics$weak_calibration$mean_calibration_ratio,
    calibration_slope_status = metrics$weak_calibration$slope_fit$status,
    calibration_slope = metrics$weak_calibration$calibration_slope,
    selection_stability_status = row$selection_stability_status,
    mean_pairwise_jaccard = row$mean_pairwise_jaccard,
    n_successful_selection_folds = row$n_successful_selection_folds,
    n_folds = row$n_folds,
    n_successful_folds = row$n_successful_folds,
    glmnet_active_set_stability_status =
      row$glmnet_active_set_stability_status,
    glmnet_active_set_mean_pairwise_jaccard =
      row$glmnet_active_set_mean_pairwise_jaccard,
    n_successful_glmnet_active_set_folds =
      row$n_successful_glmnet_active_set_folds,
    stringsAsFactors = FALSE
  )
}

#' @rdname sb_stack
#' @export
sb_run_stack_cv <- function(task, stack, folds = 5L, seed = 1L) {
  contract <- sb_stack_task_contract(task, stack)
  task <- contract$task
  stack <- contract$stack
  task_type <- contract$task_type
  split_object <- sb_make_outer_folds(task, folds, seed)
  outer_members <- sb_compare_cv(
    task,
    stack$pipelines,
    folds = folds,
    seed = seed
  )
  fold_results <- lapply(split_object$splits, function(split) {
    sb_run_one_stack_fold(
      task,
      stack,
      outer_members,
      split,
      split_object$seed
    )
  })
  predictions <- do.call(
    rbind,
    lapply(fold_results, function(result) result$predictions)
  )
  runs <- do.call(rbind, lapply(fold_results, function(result) result$run))
  inner_folds <- lapply(fold_results, function(result) result$audit)
  rownames(predictions) <- NULL
  rownames(runs) <- NULL
  patient_ids <- as.character(task$data[[task$id_col]])
  predictions <- predictions[
    match(patient_ids, predictions$patient_id),
    ,
    drop = FALSE
  ]
  rownames(predictions) <- NULL
  exactly_once <- nrow(predictions) == task$n &&
    !anyNA(predictions$patient_id) &&
    anyDuplicated(predictions$patient_id) == 0L &&
    setequal(predictions$patient_id, patient_ids)
  value_name <- if (identical(task_type, "binary")) "probability" else "risk"
  complete_oof <- exactly_once && all(predictions$status == "success") &&
    all(is.finite(predictions[[value_name]]))
  if (identical(task_type, "binary")) {
    complete_oof <- complete_oof &&
      all(predictions$probability >= 0 & predictions$probability <= 1)
    metrics <- sb_stack_binary_metrics(task, predictions, runs)
    requested_metrics_complete <- TRUE
    complete_probability <- NULL
  } else {
    complete_probability <- if ("prediction_horizon" %in% names(task)) {
      exactly_once &&
        all(predictions$survival_probability_status == "success") &&
        all(is.finite(predictions$survival_probability)) &&
        all(
          predictions$survival_probability >= 0 &
            predictions$survival_probability <= 1
        )
    } else {
      TRUE
    }
    metrics <- sb_stack_survival_metrics(
      task,
      predictions,
      split_object$splits
    )
    decision_curve_complete <-
      !("decision_thresholds" %in% names(task)) ||
        identical(metrics$decision_curve$status, "success")
    requested_metrics_complete <-
      !("prediction_horizon" %in% names(task)) ||
        (
          identical(metrics$ipcw_brier$status, "success") &&
            identical(metrics$time_auc$status, "success") &&
            identical(metrics$uno_c$status, "success") &&
            decision_curve_complete
        )
  }
  status <- if (all(runs$status == "success") && complete_oof &&
    (is.null(complete_probability) || complete_probability) &&
    requested_metrics_complete) {
    "success"
  } else {
    "failed"
  }
  fold_assignment <- data.frame(
    patient_id = patient_ids,
    fold = split_object$assignment,
    stringsAsFactors = FALSE
  )
  result <- list(
    status = status,
    task_type = task_type,
    complete_oof = complete_oof,
    predictions = predictions,
    metrics = metrics,
    runs = runs,
    folds = split_object$splits,
    inner_folds = inner_folds,
    stack = stack,
    members = outer_members,
    task = list(
      n = task$n,
      id_col = task$id_col,
      predictors = task$predictors,
      independent_unit = task$independent_unit
    ),
    provenance = list(
      seed = split_object$seed,
      folds = split_object$n_folds,
      fold_source = split_object$source,
      fold_assignment = fold_assignment,
      inner_folds = stack$inner_folds,
      inner_seed = stack$inner_seed,
      package = "stormbreaker",
      package_version = as.character(utils::packageVersion("stormbreaker")),
      r_version = R.version.string,
      glmnet_version = if (sb_glmnet_available()) {
        as.character(utils::packageVersion("glmnet"))
      } else {
        NA_character_
      }
    )
  )
  if (!is.null(complete_probability)) {
    result <- append(
      result,
      list(complete_survival_probability_oof = complete_probability),
      after = match("complete_oof", names(result))
    )
  }
  class(result) <- "stormbreaker_stack_cv_result"
  stack_row <- sb_stack_comparison_row(result)
  result$comparison <- rbind(
    outer_members$comparison,
    stack_row[names(outer_members$comparison)]
  )
  rownames(result$comparison) <- NULL
  if ("decision_thresholds" %in% names(task)) {
    result$decision_curves <- c(
      outer_members$decision_curves,
      stats::setNames(list(metrics$decision_curve), stack$id)
    )
  }
  result
}

#' @export
print.stormbreaker_stack_cv_result <- function(x, ...) {
  cat("<stormbreaker_stack_cv_result:", x$status, ">\n")
  cat(
    "Members:", length(x$stack$pipelines),
    "| Patients:", nrow(x$predictions),
    "| Outer folds:", nrow(x$runs), "\n"
  )
  cat("Complete stacked OOF:", x$complete_oof, "\n")
  if (identical(x$task_type, "binary")) {
    cat("ROC AUC:", format(round(x$metrics$roc_auc$estimate, 3)), "\n")
  } else {
    cat(
      "Harrell C:",
      format(round(x$metrics$harrell_c$estimate, 3)),
      "\n"
    )
  }
  cat("No member or stack is selected by this result.\n")
  invisible(x)
}

# ---- stacking-final ----

sb_stack_final_fit_failure <- function(result, stage, error) {
  result$status <- "failed"
  result$error_stage <- stage
  result$error_class <- class(error)[[1L]]
  result$error_message <- conditionMessage(error)
  class(result) <- "stormbreaker_stack_final_fit"
  result
}

sb_stack_final_prediction_semantics <- function(task, stack) {
  if (inherits(task, "stormbreaker_binary_task")) {
    semantics <- list(
      output = "probability",
      positive_class = 1L,
      direction = "higher_probability_of_class_1",
      meta_model = stack$meta_learner
    )
  } else {
    semantics <- list(
      output = "risk",
      risk_type = "stacked_ridge_cox_linear_predictor",
      risk_reference = "zero",
      risk_formula = "member_predictions_beta",
      risk_direction = "higher_is_higher_event_risk",
      meta_model = stack$meta_learner
    )
    if ("prediction_horizon" %in% names(task)) {
      semantics$output <- c("risk", "survival_probability")
      semantics$prediction_horizon <- task$prediction_horizon
      semantics$time_unit <- task$time_unit
    }
  }
  if ("decision_thresholds" %in% names(task)) {
    semantics$decision_thresholds <- task$decision_thresholds
  }
  semantics
}

sb_stack_member_fit_error <- function(pipeline_id, fit) {
  message <- paste0(
    "Stack member `", pipeline_id, "` failed during final fitting"
  )
  if (is.character(fit$error_message) && length(fit$error_message) == 1L &&
    !is.na(fit$error_message) && nzchar(fit$error_message)) {
    message <- paste0(message, ": ", fit$error_message)
  } else {
    message <- paste0(message, ".")
  }
  structure(
    list(message = message, call = NULL, parent = NULL),
    class = c(
      "stormbreaker_stack_member_fit_failure",
      "stormbreaker_error",
      "error",
      "condition"
    )
  )
}

#' @rdname sb_stack
#' @export
sb_fit_stack_final <- function(task, stack) {
  contract <- sb_stack_task_contract(task, stack)
  task <- contract$task
  stack <- contract$stack
  task_type <- contract$task_type
  result <- list(
    status = "not_fitted",
    task_type = task_type,
    stack = stack,
    pipeline = list(id = stack$id, learner = stack$meta_learner),
    development = sb_final_fit_development(task),
    inner = list(
      status = "not_fitted",
      folds = stack$inner_folds,
      seed = stack$inner_seed,
      fold_source = NA_character_,
      fold_assignment = NULL,
      member_comparison = NULL
    ),
    meta = NULL,
    member_fits = list(),
    required_features = character(),
    prediction = sb_stack_final_prediction_semantics(task, stack),
    provenance = list(
      package = "stormbreaker",
      package_version = as.character(utils::packageVersion("stormbreaker")),
      r_version = R.version.string,
      survival_version = as.character(utils::packageVersion("survival")),
      glmnet_version = if (sb_glmnet_available()) {
        as.character(utils::packageVersion("glmnet"))
      } else {
        NA_character_
      }
    ),
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )

  comparison_attempt <- sb_final_fit_attempt(sb_compare_cv(
    task,
    stack$pipelines,
    folds = stack$inner_folds,
    seed = stack$inner_seed
  ))
  if (!is.null(comparison_attempt$error)) {
    return(sb_stack_final_fit_failure(
      result,
      "inner_resampling",
      comparison_attempt$error
    ))
  }
  comparison <- comparison_attempt$value
  result$inner <- list(
    status = "success",
    folds = comparison$folds,
    seed = comparison$seed,
    fold_source = comparison$fold_source,
    fold_assignment = comparison$fold_assignment,
    member_comparison = comparison$comparison
  )

  matrix_attempt <- sb_final_fit_attempt(
    sb_stack_member_matrix(comparison, task, stack)
  )
  if (!is.null(matrix_attempt$error)) {
    result$inner$status <- "failed"
    return(sb_stack_final_fit_failure(
      result,
      "inner_oof",
      matrix_attempt$error
    ))
  }
  member_matrix <- matrix_attempt$value

  meta_attempt <- sb_final_fit_attempt(
    sb_fit_stack_meta(task, stack, member_matrix)
  )
  if (!is.null(meta_attempt$error)) {
    return(sb_stack_final_fit_failure(
      result,
      "meta_fit",
      meta_attempt$error
    ))
  }
  result$meta <- meta_attempt$value

  member_error <- NULL
  sb_with_selector_reuse(for (index in seq_along(stack$pipelines)) {
    pipeline_id <- stack$pipeline_ids[[index]]
    member_fit <- sb_fit_final(task, stack$pipelines[[index]])
    result$member_fits[[pipeline_id]] <- member_fit
    if (!identical(member_fit$status, "success")) {
      member_error <- sb_stack_member_fit_error(pipeline_id, member_fit)
      break
    }
  })
  if (!is.null(member_error)) {
    return(sb_stack_final_fit_failure(result, "member_fit", member_error))
  }

  used_features <- unique(unlist(lapply(
    result$member_fits,
    function(fit) fit$required_features
  )))
  result$required_features <- task$predictors[
    task$predictors %in% used_features
  ]
  result$status <- "success"
  class(result) <- "stormbreaker_stack_final_fit"
  result
}

#' @export
print.stormbreaker_stack_final_fit <- function(x, ...) {
  cat("<stormbreaker_stack_final_fit:", x$status, ">\n")
  cat(
    "Task:", x$task_type,
    "| Members:", length(x$stack$pipelines),
    "| Development patients:", x$development$n, "\n"
  )
  cat("Meta-learner:", x$stack$meta_learner, "\n")
  if (identical(x$status, "failed")) {
    cat("Failed stage:", x$error_stage, "|", x$error_class, "\n")
  }
  invisible(x)
}

sb_validate_stack_final_fit_for_prediction <- function(fit) {
  required <- c(
    "status", "task_type", "stack", "pipeline", "development", "inner",
    "meta", "member_fits", "required_features", "prediction", "provenance"
  )
  valid <- inherits(fit, "stormbreaker_stack_final_fit") && is.list(fit) &&
    length(setdiff(required, names(fit))) == 0L &&
    identical(fit$status, "success")
  if (!valid) {
    sb_abort(
      "`fit` must be a successful result from `sb_fit_stack_final()`.",
      "stormbreaker_bad_stack_final_fit"
    )
  }
  stack <- sb_validate_stack_contract(fit$stack, fit$task_type)
  valid_members <- identical(names(fit$member_fits), stack$pipeline_ids) &&
    length(fit$member_fits) == length(stack$pipelines)
  if (!valid_members) {
    sb_abort(
      "The saved stack member fits are incomplete or reordered.",
      "stormbreaker_bad_stack_final_fit"
    )
  }
  for (member_fit in fit$member_fits) {
    sb_validate_final_fit_for_prediction(member_fit)
  }
  valid_meta <- inherits(fit$meta, "stormbreaker_stack_meta_fit") &&
    is.list(fit$meta$state) &&
    identical(fit$meta$state$features, stack$pipeline_ids)
  valid_features <- is.character(fit$required_features) &&
    !anyNA(fit$required_features) && !anyDuplicated(fit$required_features)
  valid_pipeline <- is.list(fit$pipeline) &&
    identical(fit$pipeline$id, stack$id) &&
    identical(fit$pipeline$learner, stack$meta_learner)
  if (!valid_meta || !valid_features || !valid_pipeline) {
    sb_abort(
      "The saved stack meta state or prediction metadata are invalid.",
      "stormbreaker_bad_stack_final_fit"
    )
  }
  fit
}

sb_stack_member_prediction_error <- function(pipeline_id, prediction) {
  detail <- if (is.list(prediction) &&
    is.character(prediction$error_message) &&
    length(prediction$error_message) == 1L &&
    !is.na(prediction$error_message) && nzchar(prediction$error_message)) {
    paste0(": ", prediction$error_message)
  } else {
    "."
  }
  sb_abort(
    paste0(
      "Stack member `", pipeline_id,
      "` failed during outcome-blind prediction", detail
    ),
    "stormbreaker_stack_member_prediction_failure"
  )
}

sb_predict_stack_members <- function(fit, new_data, id_col, patient_ids) {
  columns <- vector("list", length(fit$stack$pipelines))
  names(columns) <- fit$stack$pipeline_ids
  for (index in seq_along(fit$stack$pipelines)) {
    pipeline_id <- fit$stack$pipeline_ids[[index]]
    prediction <- sb_predict_final(
      fit$member_fits[[pipeline_id]],
      new_data,
      id_col
    )
    if (!is.list(prediction) || !identical(prediction$status, "success")) {
      sb_stack_member_prediction_error(pipeline_id, prediction)
    }
    contract <- sb_external_prediction_contract(prediction)
    if (!identical(contract$task_type, fit$task_type) ||
      !setequal(contract$patient_ids, patient_ids)) {
      sb_abort(
        paste0(
          "Stack member `", pipeline_id,
          "` returned misaligned patient predictions."
        ),
        "stormbreaker_stack_member_prediction_failure"
      )
    }
    value <- if (identical(fit$task_type, "binary")) {
      contract$probability
    } else {
      contract$risk
    }
    columns[[pipeline_id]] <- unname(value[
      match(patient_ids, contract$patient_ids)
    ])
  }
  data.frame(
    columns,
    row.names = patient_ids,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' @rdname sb_stack
#' @export
sb_predict_stack_final <- function(fit, new_data, id_col) {
  id_attempt <- sb_final_fit_attempt(
    sb_extract_final_prediction_ids(new_data, id_col)
  )
  if (!is.null(id_attempt$error)) {
    return(sb_final_prediction_failure(
      fit,
      id_col,
      character(),
      "id",
      id_attempt$error
    ))
  }
  patient_ids <- id_attempt$value

  fit_attempt <- sb_final_fit_attempt(
    sb_validate_stack_final_fit_for_prediction(fit)
  )
  if (!is.null(fit_attempt$error)) {
    return(sb_final_prediction_failure(
      fit,
      id_col,
      patient_ids,
      "fit",
      fit_attempt$error
    ))
  }
  fit <- fit_attempt$value

  member_attempt <- sb_final_fit_attempt(
    sb_predict_stack_members(fit, new_data, id_col, patient_ids)
  )
  if (!is.null(member_attempt$error)) {
    return(sb_final_prediction_failure(
      fit,
      id_col,
      patient_ids,
      "member_prediction",
      member_attempt$error
    ))
  }
  member_matrix <- member_attempt$value

  meta_attempt <- sb_final_fit_attempt({
    values <- sb_predict_stack_meta(fit$meta$state, member_matrix)
    sb_validate_final_prediction_values(fit, values, nrow(member_matrix))
  })
  if (!is.null(meta_attempt$error)) {
    return(sb_final_prediction_failure(
      fit,
      id_col,
      patient_ids,
      "meta_prediction",
      meta_attempt$error
    ))
  }
  metadata <- sb_final_prediction_metadata(fit)
  result <- list(
    status = "success",
    task_type = metadata$task_type,
    pipeline_id = metadata$pipeline_id,
    learner = metadata$learner,
    id_col = id_col,
    predictions = sb_final_prediction_rows(
      fit,
      patient_ids,
      "success",
      NA_character_,
      NA_character_,
      meta_attempt$value
    ),
    member_predictions = member_matrix,
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )
  if (!is.null(metadata$decision_thresholds)) {
    result <- append(
      result,
      list(decision_thresholds = metadata$decision_thresholds),
      after = match("id_col", names(result))
    )
  }
  sb_as_final_prediction(result)
}

# ---- stacking-meta ----

sb_stack_task_contract <- function(task, stack) {
  stack <- sb_upgrade_legacy_names(stack)
  if (inherits(task, "stormbreaker_binary_task")) {
    task <- sb_validate_binary_task_contract(task)
    task_type <- "binary"
  } else if (inherits(task, "stormbreaker_surv_task")) {
    task <- sb_validate_task_contract(task)
    task_type <- "survival"
  } else {
    sb_abort(
      "`task` must be created by `sb_surv_task()` or `sb_binary_task()`.",
      "stormbreaker_bad_task"
    )
  }
  stack <- sb_validate_stack_contract(stack, task_type)
  list(task = task, stack = stack, task_type = task_type)
}

sb_incomplete_stack_member_oof <- function(pipeline_id) {
  sb_abort(
    paste0(
      "Stack member `", pipeline_id,
      "` does not have one complete finite out-of-fold prediction per patient."
    ),
    "stormbreaker_incomplete_stack_member_oof"
  )
}

sb_stack_member_matrix <- function(comparison, task, stack) {
  contract <- sb_stack_task_contract(task, stack)
  task <- contract$task
  stack <- contract$stack
  task_type <- contract$task_type
  valid_comparison <- inherits(
    comparison,
    "stormbreaker_cv_comparison"
  ) && is.list(comparison) && identical(comparison$task_type, task_type) &&
    is.list(comparison$results) &&
    setequal(names(comparison$results), stack$pipeline_ids)
  if (!valid_comparison) {
    sb_abort(
      "The member comparison does not match this task and stack.",
      "stormbreaker_bad_stack_member_comparison"
    )
  }

  patient_ids <- as.character(task$data[[task$id_col]])
  value_name <- if (identical(task_type, "binary")) {
    "probability"
  } else {
    "risk"
  }
  columns <- lapply(stack$pipeline_ids, function(pipeline_id) {
    result <- comparison$results[[pipeline_id]]
    predictions <- if (is.list(result)) result$predictions else NULL
    valid_frame <- isTRUE(result$complete_oof) &&
      is.data.frame(predictions) && nrow(predictions) == task$n &&
      all(c("patient_id", "status", value_name) %in% names(predictions))
    if (!valid_frame) {
      sb_incomplete_stack_member_oof(pipeline_id)
    }
    prediction_ids <- as.character(predictions$patient_id)
    valid_ids <- !anyNA(prediction_ids) &&
      !anyDuplicated(prediction_ids) && setequal(prediction_ids, patient_ids)
    status <- predictions$status
    value <- predictions[[value_name]]
    valid_values <- is.character(status) && length(status) == task$n &&
      !anyNA(status) && all(status == "success") &&
      is.numeric(value) && length(value) == task$n &&
      !anyNA(value) && all(is.finite(value))
    if (identical(task_type, "binary")) {
      valid_values <- valid_values && all(value >= 0 & value <= 1)
    }
    if (!valid_ids || !valid_values) {
      sb_incomplete_stack_member_oof(pipeline_id)
    }
    unname(value[match(patient_ids, prediction_ids)])
  })
  names(columns) <- stack$pipeline_ids
  data.frame(
    columns,
    row.names = patient_ids,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

sb_validate_stack_member_matrix <- function(member_matrix, task, stack) {
  patient_ids <- as.character(task$data[[task$id_col]])
  valid <- is.data.frame(member_matrix) && nrow(member_matrix) == task$n &&
    identical(names(member_matrix), stack$pipeline_ids) &&
    identical(rownames(member_matrix), patient_ids) &&
    all(vapply(member_matrix, is.numeric, logical(1))) &&
    all(vapply(
      member_matrix,
      function(column) !anyNA(column) && all(is.finite(column)),
      logical(1)
    ))
  if (valid && identical(stack$task_type, "binary")) {
    values <- as.matrix(member_matrix)
    valid <- all(values >= 0 & values <= 1)
  }
  if (!valid) {
    sb_abort(
      paste0(
        "`member_matrix` must contain one ordered finite column for every ",
        "stack member and one ordered row for every task patient."
      ),
      "stormbreaker_bad_stack_member_matrix"
    )
  }
  member_matrix
}

sb_fit_stack_meta <- function(task, stack, member_matrix) {
  contract <- sb_stack_task_contract(task, stack)
  task <- contract$task
  stack <- contract$stack
  task_type <- contract$task_type
  member_matrix <- sb_validate_stack_member_matrix(
    member_matrix,
    task,
    stack
  )
  patient_ids <- as.character(task$data[[task$id_col]])

  if (identical(task_type, "binary")) {
    outcome <- as.integer(task$data[[task$outcome_col]])
    tuned <- sb_select_glmnet_logistic_lambda(
      member_matrix,
      outcome,
      inner_folds = stack$inner_folds,
      seed = stack$inner_seed,
      alpha = stack$meta_alpha,
      standardize = stack$meta_standardize,
      selection_rule = stack$meta_selection_rule,
      measure = stack$meta_measure
    )
    names(tuned$fold_id) <- patient_ids
    adapter <- sb_glmnet_logistic_learner_adapter(
      alpha = tuned$alpha,
      lambda = tuned$lambda,
      standardize = tuned$standardize
    )
    model <- adapter$fit(member_matrix, outcome)
    coefficients <- adapter$extract_coefficients(model)
    state <- adapter$extract_state(model)
  } else {
    time <- task$data[[task$time_col]]
    event <- as.integer(task$data[[task$event_col]])
    tuned <- sb_select_glmnet_cox_lambda(
      member_matrix,
      time,
      event,
      inner_folds = stack$inner_folds,
      seed = stack$inner_seed,
      alpha = stack$meta_alpha,
      standardize = stack$meta_standardize,
      selection_rule = stack$meta_selection_rule,
      measure = stack$meta_measure
    )
    names(tuned$fold_id) <- patient_ids
    adapter <- sb_glmnet_cox_learner_adapter(
      alpha = tuned$alpha,
      lambda = tuned$lambda,
      standardize = tuned$standardize
    )
    model <- adapter$fit(member_matrix, time, event)
    coefficients <- adapter$extract_coefficients(model)
    context <- list(x = member_matrix, time = time, event = event)
    if ("prediction_horizon" %in% names(task)) {
      context$prediction_horizon <- task$prediction_horizon
      context$time_unit <- task$time_unit
    }
    state <- adapter$extract_state(model, context)
  }

  structure(
    list(
      task_type = task_type,
      meta_learner = stack$meta_learner,
      tuning = tuned,
      coefficients = coefficients,
      state = state
    ),
    class = "stormbreaker_stack_meta_fit"
  )
}

sb_predict_stack_meta <- function(meta_state, member_matrix) {
  member_matrix <- sb_validate_predictor_frame(
    member_matrix,
    "member_matrix"
  )
  sb_predict_final_state(meta_state, member_matrix)
}

# ---- stacking-spec ----

#' @rdname sb_stack
#' @export
sb_stack <- function(pipelines, inner_folds = 5L, inner_seed = 1L) {
  valid_collection <- is.list(pipelines) && length(pipelines) >= 2L &&
    all(vapply(
      pipelines,
      function(pipeline) inherits(pipeline, "stormbreaker_pipeline"),
      logical(1)
    ))
  if (!valid_collection) {
    sb_abort(
      paste0(
        "`pipelines` must contain at least two specifications returned by ",
        "`sb_pipeline()`."
      ),
      "stormbreaker_bad_stack_pipelines"
    )
  }

  task_types <- vapply(
    pipelines,
    function(pipeline) pipeline$task_type,
    character(1)
  )
  if (length(unique(task_types)) != 1L) {
    sb_abort(
      "Every stack member must have the same task type.",
      "stormbreaker_stack_task_mismatch"
    )
  }
  task_type <- task_types[[1L]]
  pipelines <- lapply(
    pipelines,
    sb_validate_pipeline_contract,
    task_type = task_type
  )
  pipeline_ids <- unname(vapply(
    pipelines,
    function(pipeline) pipeline$id,
    character(1)
  ))
  if (anyDuplicated(pipeline_ids) > 0L) {
    sb_abort(
      "Every stack member must have a distinct pipeline ID.",
      "stormbreaker_duplicate_stack_pipeline_id"
    )
  }

  inner_folds <- sb_check_positive_integer(
    inner_folds,
    "inner_folds",
    minimum = 3L
  )
  inner_seed <- sb_check_positive_integer(inner_seed, "inner_seed")
  meta_learner <- if (identical(task_type, "binary")) {
    "ridge_logistic"
  } else {
    "ridge_cox"
  }

  structure(
    list(
      id = paste0(
        "ridge_stack__",
        paste(pipeline_ids, collapse = "__")
      ),
      task_type = task_type,
      pipelines = pipelines,
      pipeline_ids = pipeline_ids,
      meta_learner = meta_learner,
      meta_alpha = 0,
      meta_standardize = TRUE,
      meta_selection_rule = "lambda.min",
      meta_measure = "deviance",
      inner_folds = inner_folds,
      inner_seed = inner_seed
    ),
    class = "stormbreaker_stack"
  )
}

sb_validate_stack_contract <- function(stack, task_type = NULL) {
  required <- c(
    "id", "task_type", "pipelines", "pipeline_ids", "meta_learner",
    "meta_alpha", "meta_standardize", "meta_selection_rule",
    "meta_measure", "inner_folds", "inner_seed"
  )
  if (!inherits(stack, "stormbreaker_stack") || !is.list(stack) ||
    length(setdiff(required, names(stack))) > 0L) {
    sb_abort(
      "`stack` must be created by `sb_stack()`.",
      "stormbreaker_bad_stack"
    )
  }

  if (!is.null(task_type) && !identical(stack$task_type, task_type)) {
    sb_abort(
      "The stack does not match the supplied task type.",
      "stormbreaker_stack_task_mismatch"
    )
  }
  stack
}

sb_subset_task <- function(task, rows) {
  valid_rows <- is.numeric(rows) && length(rows) >= 2L &&
    !anyNA(rows) && all(is.finite(rows)) &&
    all(rows == as.integer(rows)) && !anyDuplicated(rows)
  if (!valid_rows) {
    sb_abort(
      "Task row indices must be distinct finite whole numbers.",
      "stormbreaker_bad_task_rows"
    )
  }
  rows <- as.integer(rows)

  if (inherits(task, "stormbreaker_binary_task")) {
    task <- sb_validate_binary_task_contract(task)
    if (any(rows < 1L | rows > task$n)) {
      sb_abort(
        "Task row indices fall outside the available patients.",
        "stormbreaker_bad_task_rows"
      )
    }
    thresholds <- if ("decision_thresholds" %in% names(task)) {
      task$decision_thresholds
    } else {
      NULL
    }
    return(sb_binary_task(
      task$data[rows, , drop = FALSE],
      task$id_col,
      task$outcome_col,
      task$predictors,
      decision_thresholds = thresholds,
      group_col = task$group_col
    ))
  }

  if (inherits(task, "stormbreaker_surv_task")) {
    task <- sb_validate_task_contract(task)
    if (any(rows < 1L | rows > task$n)) {
      sb_abort(
        "Task row indices fall outside the available patients.",
        "stormbreaker_bad_task_rows"
      )
    }
    horizon <- if ("prediction_horizon" %in% names(task)) {
      task$prediction_horizon
    } else {
      NULL
    }
    thresholds <- if ("decision_thresholds" %in% names(task)) {
      task$decision_thresholds
    } else {
      NULL
    }
    return(sb_surv_task(
      task$data[rows, , drop = FALSE],
      task$id_col,
      task$time_col,
      task$event_col,
      task$predictors,
      task$time_unit,
      prediction_horizon = horizon,
      decision_thresholds = thresholds,
      group_col = task$group_col
    ))
  }

  sb_abort(
    "`task` must be created by `sb_surv_task()` or `sb_binary_task()`.",
    "stormbreaker_bad_task"
  )
}

#' @export
print.stormbreaker_stack <- function(x, ...) {
  cat("<stormbreaker_stack>\n")
  cat("Task:", x$task_type, "| Members:", length(x$pipelines), "\n")
  cat("Meta-learner:", x$meta_learner, "| Inner folds:", x$inner_folds, "\n")
  invisible(x)
}
