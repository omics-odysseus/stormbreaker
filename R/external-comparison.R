#' External model comparisons
#' @name sb_compare_external
#' @description Compare locked predictions on one external endpoint and collect compatible cohort results.
#' @rawRd \arguments{
#' \item{predictions}{
#' A non-empty list of distinct final-prediction objects
#' returned by \code{\link[=sb_predict_final]{sb_predict_final()}} or \code{\link[=sb_predict_signature]{sb_predict_signature()}}. Each object
#' keeps its own failure status, but every available pipeline ID must be
#' distinct.
#' }
#' \item{endpoint}{
#' One outcome-only object returned by \code{\link[=sb_binary_endpoint]{sb_binary_endpoint()}} or
#' \code{\link[=sb_surv_endpoint]{sb_surv_endpoint()}}. The complete patient-ID set is reused for every
#' prediction.
#' }
#' \item{comparisons}{
#' A non-empty named list of
#' \code{stormbreaker_external_comparison} objects returned by
#' \code{\link[=sb_compare_external]{sb_compare_external()}}. Names must be unique, non-blank cohort IDs. A
#' single comparison must still be wrapped, for example
#' \code{list(cohort_a = comparison)}.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_compare_external()}}{
#' \code{sb_compare_external()} evaluates a non-empty set of already computed final
#' predictions against the same outcome-only external endpoint. It calls
#' \code{\link[=sb_evaluate_external]{sb_evaluate_external()}} once per prediction, so exact patient-ID alignment,
#' horizon semantics, metric definitions, and failure behavior are unchanged.
#' The function has no access to fitted models or external predictor values and
#' does not refit, recalibrate, or replace a learner.
#' 
#' \strong{Returns:} A \code{stormbreaker_external_comparison} containing a descriptive
#' \code{comparison} table, every complete named external \code{evaluation}, the task
#' type, effective \code{prediction_target}, endpoint summary, and patient count.
#' When present, named complete \code{decision_curves} are copied from the
#' evaluations.
#' 
#' \strong{Details:} The returned table copies scalar metrics from each complete evaluation.
#' Every full \code{stormbreaker_external_evaluation} remains available in
#' \code{evaluations}; a failed or incomplete evaluation stays visible and does not
#' erase the other predictions. Optional decision curves are retained as named
#' complete objects rather than flattened or integrated.
#' 
#' For survival, the comparison retains one effective \code{prediction_target}. A
#' horizon declared by the endpoint is authoritative. If the endpoint omits a
#' horizon, a horizon already locked into the predictions is retained instead;
#' all otherwise evaluable survival predictions must then be either risk-only
#' or share that same time unit and horizon. This prevents unlike time targets
#' from being presented as one comparison while preserving the endpoint's
#' original declaration separately in \code{endpoint}. A failed prediction can
#' still contribute target metadata when its locked horizon and time unit are
#' recoverable. Otherwise the target is marked \code{partial} or \code{unknown} rather
#' than being mislabeled as risk-only.
#' 
#' Pass predictions fixed before inspecting the external outcomes. The result
#' compares those predictions on the supplied cohort. Use development data to
#' choose a pipeline and threshold; reserve final-test outcomes for evaluation.
#' For paired intervals, pass the result to
#' \code{\link[=sb_bootstrap_external_comparison]{sb_bootstrap_external_comparison()}}.
#' }
#' 
#' \subsection{\code{sb_collect_external_comparisons()}}{
#' \code{sb_collect_external_comparisons()} keeps several completed
#' \code{\link[=sb_compare_external]{sb_compare_external()}} results in one ordered object. The outer list names
#' are cohort IDs. Each source comparison and every failed model row are
#' retained; the long \code{results} table only adds \code{cohort_id} and normalizes
#' pipeline order to the first cohort.
#' 
#' \strong{Returns:} A \code{stormbreaker_external_collection} containing overall \code{status},
#' common \code{task_type}, ordered \code{pipeline_ids}, a one-row-per-cohort \code{cohorts}
#' table, a long \code{results} table by cohort, and the unchanged named
#' source \code{comparisons}.
#' 
#' \strong{Details:} Cohorts must contain the same task, pipeline IDs, and learner labels. For a
#' survival task they must also use the same time unit and either the same
#' effective prediction horizon or no horizon in every cohort. The effective
#' target is the endpoint-declared horizon, or the common locked prediction
#' horizon when the endpoint omitted it. These checks prevent unlike evaluation
#' targets from being placed in one table. A comparison whose target remained
#' unknown or only partially recoverable cannot be collected. The checks use
#' stored metadata. Interpret the results using the study's cohort roles and
#' independence, eligibility, outcome definitions, preprocessing and fitted
#' models or published formulas.
#' 
#' The collection presents observed metrics by cohort in the supplied pipeline
#' order. The original comparison objects retain predictions and uncertainty
#' estimates for each cohort.
#' }
#' }
#' @rawRd \examples{
#' # sb_compare_external
#' development <- sb_example_binary(n = 50, seed = 1)
#' task <- sb_binary_task(
#'   development,
#'   "patient_id",
#'   "outcome",
#'   c("gene_a", "gene_b", "gene_c")
#' )
#' pipelines <- list(
#'   sb_pipeline(learner = "logistic"),
#'   sb_pipeline(learner = "gaussian_naive_bayes")
#' )
#' fits <- lapply(pipelines, function(pipeline) {
#'   sb_fit_final(task, pipeline)
#' })
#' external <- sb_example_binary(n = 30, seed = 2)
#' external$patient_id <- paste0("external_", external$patient_id)
#' predictions <- lapply(fits, function(fit) {
#'   sb_predict_final(
#'     fit,
#'     external[c("patient_id", "gene_a", "gene_b", "gene_c")],
#'     "patient_id"
#'   )
#' })
#' endpoint <- sb_binary_endpoint(external, "patient_id", "outcome")
#' sb_compare_external(predictions, endpoint)
#' 
#' # sb_collect_external_comparisons
#' data <- data.frame(
#'   patient_id = paste0("p", 1:12),
#'   outcome = rep(c(0L, 1L), 6),
#'   gene_a = seq(-1, 1, length.out = 12)
#' )
#' signature <- sb_linear_signature(
#'   "published_example", "binary", c(gene_a = 0.5), intercept = 0
#' )
#' prediction <- sb_predict_signature(
#'   signature, data[c("patient_id", "gene_a")], "patient_id"
#' )
#' comparison <- sb_compare_external(
#'   list(prediction),
#'   sb_binary_endpoint(data, "patient_id", "outcome")
#' )
#' sb_collect_external_comparisons(list(cohort_a = comparison))
#' }
NULL

# ---- external-comparison ----

# ---- bootstrap-external-comparison ----

sb_external_paired_metric_spec <- function(task_type, has_horizon) {
  if (identical(task_type, "binary")) {
    return(data.frame(
      metric = c(
        "roc_auc", "average_precision", "brier_score", "log_loss"
      ),
      direction = c(
        "higher_is_better", "higher_is_better",
        "lower_is_better", "lower_is_better"
      ),
      stringsAsFactors = FALSE
    ))
  }
  metric <- "harrell_c"
  direction <- "higher_is_better"
  if (has_horizon) {
    metric <- c(metric, "ipcw_brier_score", "time_auc", "uno_c")
    direction <- c(
      direction,
      "lower_is_better", "higher_is_better", "higher_is_better"
    )
  }
  data.frame(metric = metric, direction = direction, stringsAsFactors = FALSE)
}

sb_external_paired_metric_rows <- function(scored, task_type, metric_spec) {
  rows <- sb_external_bootstrap_metric_rows(scored, task_type)
  positions <- match(metric_spec$metric, rows$metric)
  if (anyNA(positions)) {
    sb_abort(
      "A required paired-comparison metric is unavailable from the runtime.",
      "stormbreaker_missing_paired_metric"
    )
  }
  rows[positions, , drop = FALSE]
}

sb_external_paired_candidate_contract <- function(evaluation, reference) {
  candidate <- sb_external_bootstrap_contract(evaluation)
  same_design <- identical(candidate$task_type, reference$task_type) &&
    identical(candidate$n, reference$n) &&
    identical(candidate$has_horizon, reference$has_horizon) &&
    identical(candidate$prediction_horizon, reference$prediction_horizon) &&
    identical(candidate$time_unit, reference$time_unit)
  endpoint_columns <- if (identical(reference$task_type, "binary")) {
    c("patient_id", "outcome")
  } else {
    c("patient_id", "time", "event")
  }
  same_endpoint <- same_design && identical(
    candidate$scored[endpoint_columns],
    reference$scored[endpoint_columns]
  )
  if (!same_endpoint) {
    sb_abort(
      paste0(
        "Every paired external evaluation must contain the same patients ",
        "and endpoint rows in the same endpoint-defined order."
      ),
      "stormbreaker_unpaired_external_evaluation"
    )
  }
  candidate
}

sb_external_paired_empty_differences <- function(replicates = FALSE) {
  result <- data.frame(
    reference_pipeline_id = character(),
    candidate_pipeline_id = character(),
    contrast = character(),
    metric = character(),
    metric_direction = character(),
    reference_estimate = numeric(),
    candidate_estimate = numeric(),
    estimate = numeric(),
    status = character(),
    error_class = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
  if (replicates) {
    result <- data.frame(
      replicate = integer(),
      result,
      stringsAsFactors = FALSE
    )
  }
  result
}

sb_external_paired_empty_intervals <- function() {
  data.frame(
    reference_pipeline_id = character(),
    candidate_pipeline_id = character(),
    contrast = character(),
    metric = character(),
    metric_direction = character(),
    reference_point_estimate = numeric(),
    candidate_point_estimate = numeric(),
    point_estimate = numeric(),
    lower = numeric(),
    upper = numeric(),
    confidence_level = numeric(),
    status = character(),
    interval_available = logical(),
    replicate_completeness = character(),
    n_requested = integer(),
    n_finite_successful = integer(),
    n_unavailable = integer(),
    n_nonfinite = integer(),
    valid_fraction = numeric(),
    error_class = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
}

sb_external_paired_difference_rows <- function(
    source,
    reference_pipeline_id,
    candidate_pipeline_id,
    metric_spec) {
  candidate <- source[
    source$pipeline_id == candidate_pipeline_id,
    ,
    drop = FALSE
  ]
  reference <- source[
    source$pipeline_id == reference_pipeline_id,
    ,
    drop = FALSE
  ]
  has_replicate <- "replicate" %in% names(source)
  candidate_key <- if (has_replicate) {
    paste(candidate$replicate, candidate$metric, sep = "\r")
  } else {
    candidate$metric
  }
  reference_key <- if (has_replicate) {
    paste(reference$replicate, reference$metric, sep = "\r")
  } else {
    reference$metric
  }
  reference <- reference[match(candidate_key, reference_key), , drop = FALSE]

  available <- !is.na(candidate$status) & !is.na(reference$status) &
    candidate$status == "success" & reference$status == "success"
  difference <- candidate$estimate - reference$estimate
  difference[!available] <- NA_real_
  error_class <- rep(NA_character_, nrow(candidate))
  error_message <- rep(NA_character_, nrow(candidate))
  for (index in which(!available)) {
    candidate_failed <- !identical(candidate$status[[index]], "success")
    source_row <- if (candidate_failed) candidate else reference
    error_class[[index]] <- source_row$error_class[[index]]
    error_message[[index]] <- source_row$error_message[[index]]
    if (is.na(error_class[[index]]) || !nzchar(error_class[[index]])) {
      error_class[[index]] <- "stormbreaker_unavailable_paired_metric"
    }
    if (is.na(error_message[[index]]) || !nzchar(error_message[[index]])) {
      error_message[[index]] <- paste0(
        "Both models require an available `", candidate$metric[[index]],
        "` estimate in the same patient draw."
      )
    }
  }
  result <- data.frame(
    reference_pipeline_id = rep(reference_pipeline_id, nrow(candidate)),
    candidate_pipeline_id = rep(candidate_pipeline_id, nrow(candidate)),
    contrast = rep("candidate_minus_reference", nrow(candidate)),
    metric = candidate$metric,
    metric_direction = metric_spec$direction[match(
      candidate$metric,
      metric_spec$metric
    )],
    reference_estimate = reference$estimate,
    candidate_estimate = candidate$estimate,
    estimate = difference,
    status = ifelse(available, "success", "not_evaluable"),
    error_class = error_class,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
  if (has_replicate) {
    result <- data.frame(
      replicate = candidate$replicate,
      result,
      stringsAsFactors = FALSE
    )
  }
  result
}

sb_external_paired_interval_rows <- function(
    point_differences,
    differences,
    confidence_level,
    reps) {
  interval_input <- point_differences[c(
    "metric", "estimate", "status", "error_class", "error_message"
  )]
  replicate_input <- differences[c(
    "replicate", "metric", "estimate", "status",
    "error_class", "error_message"
  )]
  intervals <- sb_external_percentile_intervals(
    interval_input,
    replicate_input,
    confidence_level,
    reps
  )
  positions <- match(intervals$metric, point_differences$metric)
  data.frame(
    reference_pipeline_id = point_differences$reference_pipeline_id[positions],
    candidate_pipeline_id = point_differences$candidate_pipeline_id[positions],
    contrast = point_differences$contrast[positions],
    metric = intervals$metric,
    metric_direction = point_differences$metric_direction[positions],
    reference_point_estimate =
      point_differences$reference_estimate[positions],
    candidate_point_estimate =
      point_differences$candidate_estimate[positions],
    intervals[names(intervals) != "metric"],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

sb_external_paired_model_rows <- function(
    contract,
    metric_spec,
    pipeline_id,
    replicate = NULL,
    indices = NULL) {
  scored <- if (is.null(indices)) {
    contract$scored
  } else {
    sampled <- contract$scored[indices, , drop = FALSE]
    rownames(sampled) <- NULL
    sampled
  }
  rows <- sb_external_paired_metric_rows(
    scored,
    contract$task_type,
    metric_spec
  )
  result <- data.frame(
    pipeline_id = rep(pipeline_id, nrow(rows)),
    rows,
    stringsAsFactors = FALSE
  )
  if (!is.null(replicate)) {
    result <- data.frame(
      replicate = rep.int(as.integer(replicate), nrow(result)),
      result,
      stringsAsFactors = FALSE
    )
  }
  result
}

#' @rdname sb_bootstrap_external
#' @export
sb_bootstrap_external_comparison <- function(
    comparison,
    reference_pipeline_id,
    reps = 1000L,
    confidence_level = 0.95,
    seed = 1L) {
  valid_comparison <- inherits(
    comparison,
    "stormbreaker_external_comparison"
  ) && is.list(comparison) && is.list(comparison$evaluations) &&
    length(comparison$evaluations) >= 2L
  if (!valid_comparison) {
    sb_abort(
      paste0(
        "`comparison` must be returned by `sb_compare_external()` and ",
        "contain at least two predictions."
      ),
      "stormbreaker_bad_external_comparison"
    )
  }
  pipeline_ids <- vapply(comparison$evaluations, function(evaluation) {
    pipeline_id <- evaluation[["pipeline_id"]]
    if (sb_external_scalar_character(pipeline_id)) {
      pipeline_id
    } else {
      NA_character_
    }
  }, character(1L))
  valid_reference <- sb_external_scalar_character(reference_pipeline_id) &&
    reference_pipeline_id %in% pipeline_ids
  if (!valid_reference) {
    sb_abort(
      "`reference_pipeline_id` must identify one compared prediction.",
      "stormbreaker_bad_reference_pipeline"
    )
  }
  reps <- sb_check_positive_integer(reps, "reps", minimum = 2L)
  seed <- sb_check_positive_integer(seed, "seed")
  valid_confidence <- is.numeric(confidence_level) &&
    length(confidence_level) == 1L && !is.na(confidence_level) &&
    is.finite(confidence_level) && confidence_level > 0 &&
    confidence_level < 1
  if (!valid_confidence) {
    sb_abort(
      "`confidence_level` must be one finite number strictly between 0 and 1.",
      "stormbreaker_bad_confidence_level"
    )
  }
  confidence_level <- as.numeric(confidence_level)

  reference_position <- match(reference_pipeline_id, pipeline_ids)
  reference <- sb_external_bootstrap_contract(
    comparison$evaluations[[reference_position]]
  )
  metric_spec <- sb_external_paired_metric_spec(
    reference$task_type,
    reference$has_horizon
  )
  candidate_ids <- pipeline_ids[-reference_position]
  candidate_evaluations <- comparison$evaluations[-reference_position]
  candidate_attempts <- lapply(candidate_evaluations, function(evaluation) {
    sb_external_attempt(
      sb_external_paired_candidate_contract(evaluation, reference)
    )
  })
  candidate_status <- do.call(rbind, lapply(
    seq_along(candidate_ids),
    function(index) {
      error <- candidate_attempts[[index]]$error
      data.frame(
        candidate_pipeline_id = candidate_ids[[index]],
        status = if (is.null(error)) "success" else "failed",
        error_class = if (is.null(error)) {
          NA_character_
        } else {
          class(error)[[1L]]
        },
        error_message = if (is.null(error)) {
          NA_character_
        } else {
          conditionMessage(error)
        },
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(candidate_status) <- NULL
  eligible <- which(candidate_status$status == "success")
  contracts <- list(reference)
  names(contracts) <- reference_pipeline_id
  if (length(eligible) > 0L) {
    eligible_contracts <- lapply(
      candidate_attempts[eligible],
      `[[`,
      "value"
    )
    names(eligible_contracts) <- candidate_ids[eligible]
    contracts <- c(contracts, eligible_contracts)
  }

  model_point_estimates <- do.call(rbind, lapply(names(contracts), function(id) {
    sb_external_paired_model_rows(
      contracts[[id]],
      metric_spec,
      id
    )
  }))
  rownames(model_point_estimates) <- NULL
  index_plan <- withr::with_seed(seed, {
    lapply(seq_len(reps), function(replicate) {
      sample.int(reference$n, reference$n, replace = TRUE)
    })
  })
  sampled <- lapply(seq_len(reps), function(replicate) {
      indices <- index_plan[[replicate]]
      reference_sample <- reference$scored[indices, , drop = FALSE]
      rownames(reference_sample) <- NULL
      list(
        rows = do.call(rbind, lapply(names(contracts), function(id) {
          sb_external_paired_model_rows(
            contracts[[id]],
            metric_spec,
            id,
            replicate,
            indices
          )
        })),
        audit = sb_external_bootstrap_audit_row(
          reference_sample,
          reference$task_type,
          replicate
        )
      )
  })
  model_replicates <- do.call(rbind, lapply(sampled, `[[`, "rows"))
  resample_audit <- do.call(rbind, lapply(sampled, `[[`, "audit"))
  rownames(model_replicates) <- NULL
  rownames(resample_audit) <- NULL

  point_pieces <- lapply(candidate_ids[eligible], function(candidate_id) {
    sb_external_paired_difference_rows(
      model_point_estimates,
      reference_pipeline_id,
      candidate_id,
      metric_spec
    )
  })
  difference_pieces <- lapply(candidate_ids[eligible], function(candidate_id) {
    sb_external_paired_difference_rows(
      model_replicates,
      reference_pipeline_id,
      candidate_id,
      metric_spec
    )
  })
  if (length(point_pieces) == 0L) {
    point_differences <- sb_external_paired_empty_differences()
    differences <- sb_external_paired_empty_differences(replicates = TRUE)
    intervals <- sb_external_paired_empty_intervals()
  } else {
    point_differences <- do.call(rbind, point_pieces)
    differences <- do.call(rbind, difference_pieces)
    rownames(point_differences) <- NULL
    rownames(differences) <- NULL
    interval_pieces <- lapply(candidate_ids[eligible], function(candidate_id) {
      point_rows <- point_differences[
        point_differences$candidate_pipeline_id == candidate_id,
        ,
        drop = FALSE
      ]
      difference_rows <- differences[
        differences$candidate_pipeline_id == candidate_id,
        ,
        drop = FALSE
      ]
      sb_external_paired_interval_rows(
        point_rows,
        difference_rows,
        confidence_level,
        reps
      )
    })
    intervals <- do.call(rbind, interval_pieces)
    rownames(intervals) <- NULL
  }

  complete <- all(candidate_status$status == "success") &&
    nrow(intervals) > 0L && all(intervals$status == "success")
  result <- list(
    status = if (complete) "success" else "incomplete",
    task_type = reference$task_type,
    reference_pipeline_id = reference_pipeline_id,
    candidate_status = candidate_status,
    n = reference$n,
    prediction_target = reference$prediction_target,
    reps = reps,
    confidence_level = confidence_level,
    seed = seed,
    metrics = metric_spec$metric,
    pairing = list(
      resampling_unit = "patient",
      draw_sharing = "same_patient_indices_for_all_models",
      contrast = "candidate_minus_reference"
    ),
    interval_method = sb_external_bootstrap_interval_method(confidence_level),
    resampling = sb_external_bootstrap_resampling(
      reference$task_type,
      reference$has_horizon
    ),
    model_point_estimates = model_point_estimates,
    point_differences = point_differences,
    model_replicates = model_replicates,
    differences = differences,
    intervals = intervals,
    resample_audit = resample_audit,
    provenance = sb_external_bootstrap_provenance()
  )
  structure(result, class = "stormbreaker_external_bootstrap_comparison")
}

#' @export
print.stormbreaker_external_bootstrap_comparison <- function(x, ...) {
  cat("<stormbreaker_external_bootstrap_comparison:", x$status, ">\n")
  cat("paired fixed-prediction external patient bootstrap\n")
  cat("Reference:", x$reference_pipeline_id, "| Task:", x$task_type, "\n")
  cat(
    "Patients:", x$n,
    "| Replicates:", x$reps,
    "| Confidence:", format(x$confidence_level), "\n"
  )
  if (nrow(x$intervals) > 0L) {
    print(
      x$intervals[c(
        "candidate_pipeline_id", "metric", "point_estimate",
        "lower", "upper", "status"
      )],
      row.names = FALSE,
      digits = 3
    )
  }
  failed <- sum(x$candidate_status$status == "failed")
  if (failed > 0L) {
    cat("Failed candidates retained:", failed, "\n")
  }
  cat(
    paste0(
      "No ranking, winner, or p value is produced; the reference must be ",
      "declared before external outcomes are inspected.\n"
    )
  )
  invisible(x)
}

# ---- collect-external-comparisons ----

sb_external_collection_member <- function(comparison, cohort_id) {
  if (!inherits(comparison, "stormbreaker_external_comparison") ||
      !is.list(comparison)) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` is not a result from ",
        "`sb_compare_external()`."
      ),
      "stormbreaker_bad_external_collection"
    )
  }

  task_type <- comparison[["task_type"]]
  status <- comparison[["status"]]
  n_patients <- comparison[["n_patients"]]
  table <- comparison[["comparison"]]
  valid_core <- is.character(task_type) && length(task_type) == 1L &&
    !is.na(task_type) && task_type %in% c("binary", "survival") &&
    is.character(status) && length(status) == 1L && !is.na(status) &&
    status %in% c("success", "incomplete") &&
    is.numeric(n_patients) && length(n_patients) == 1L &&
    !is.na(n_patients) && is.finite(n_patients) && n_patients >= 1 &&
    n_patients == as.integer(n_patients) &&
    is.data.frame(table) && nrow(table) >= 1L
  if (!valid_core) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has incomplete external-comparison metadata."
      ),
      "stormbreaker_bad_external_collection"
    )
  }

  pipeline_ids <- table[["pipeline_id"]]
  learners <- table[["learner"]]
  valid_pipelines <- is.character(pipeline_ids) &&
    length(pipeline_ids) == nrow(table) && !anyNA(pipeline_ids) &&
    all(nzchar(trimws(pipeline_ids))) && anyDuplicated(pipeline_ids) == 0L &&
    is.character(learners) && length(learners) == nrow(table) &&
    !anyNA(learners) && all(nzchar(trimws(learners)))
  if (!valid_pipelines) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` does not retain unique pipeline IDs and learner labels."
      ),
      "stormbreaker_bad_external_collection"
    )
  }

  time_unit <- NULL
  prediction_horizon <- NULL
  if (identical(task_type, "survival")) {
    target <- comparison[["prediction_target"]]
    if (!is.list(target)) {
      sb_abort(
        paste0(
          "Cohort `", cohort_id,
          "` does not retain a usable survival evaluation target."
        ),
        "stormbreaker_bad_external_collection"
      )
    }
    target_status <- target[["target_status"]]
    valid_target_status <- is.character(target_status) &&
      length(target_status) == 1L && !is.na(target_status) &&
      target_status %in% c("complete", "partial", "unknown")
    if (!valid_target_status) {
      sb_abort(
        paste0(
          "Cohort `", cohort_id,
          "` does not retain a usable survival evaluation target."
        ),
        "stormbreaker_bad_external_collection"
      )
    }
    if (!identical(target_status, "complete")) {
      sb_abort(
        paste0(
          "Cohort `", cohort_id,
          "` has an unknown or partially recoverable survival target."
        ),
        "stormbreaker_external_target_mismatch"
      )
    }
    has_horizon <- target[["has_horizon"]]
    time_unit <- target[["time_unit"]]
    prediction_horizon <- target[["prediction_horizon"]]
    valid_has_horizon <- is.logical(has_horizon) &&
      length(has_horizon) == 1L && !is.na(has_horizon)
    valid_time_unit <- is.character(time_unit) && length(time_unit) == 1L &&
      !is.na(time_unit) && nzchar(trimws(time_unit))
    valid_horizon <- valid_has_horizon && (
      (!has_horizon && is.null(prediction_horizon)) ||
      (has_horizon && is.numeric(prediction_horizon) &&
        length(prediction_horizon) == 1L &&
        !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
        prediction_horizon > 0)
    )
    if (!valid_time_unit || !valid_horizon) {
      sb_abort(
        paste0(
          "Cohort `", cohort_id,
          "` does not retain a usable survival evaluation target."
        ),
        "stormbreaker_bad_external_collection"
      )
    }
  }

  list(
    status = status,
    task_type = task_type,
    n_patients = as.integer(n_patients),
    table = table,
    pipeline_ids = pipeline_ids,
    learners = learners,
    time_unit = time_unit,
    prediction_horizon = prediction_horizon
  )
}

#' @rdname sb_compare_external
#' @export
sb_collect_external_comparisons <- function(comparisons) {
  valid_collection <- is.list(comparisons) && length(comparisons) > 0L &&
    !inherits(comparisons, "stormbreaker_external_comparison") &&
    all(vapply(
      comparisons,
      function(comparison) {
        inherits(comparison, "stormbreaker_external_comparison")
      },
      logical(1L)
    ))
  if (!valid_collection) {
    sb_abort(
      paste0(
        "`comparisons` must be a non-empty named list returned by ",
        "`sb_compare_external()`; wrap one object as ",
        "`list(cohort_id = comparison)`."
      ),
      "stormbreaker_bad_external_collection"
    )
  }

  cohort_ids <- names(comparisons)
  valid_ids <- !is.null(cohort_ids) &&
    length(cohort_ids) == length(comparisons) &&
    !anyNA(cohort_ids) && all(nzchar(trimws(cohort_ids))) &&
    anyDuplicated(cohort_ids) == 0L
  if (!valid_ids) {
    sb_abort(
      "Outer list names must be unique, non-blank cohort IDs.",
      "stormbreaker_bad_external_collection"
    )
  }

  members <- lapply(seq_along(comparisons), function(index) {
    sb_external_collection_member(comparisons[[index]], cohort_ids[[index]])
  })
  task_types <- vapply(members, `[[`, character(1L), "task_type")
  if (length(unique(task_types)) != 1L) {
    sb_abort(
      "All collected external comparisons must use the same task type.",
      "stormbreaker_external_task_mismatch"
    )
  }

  reference_ids <- members[[1L]]$pipeline_ids
  reference_learners <- members[[1L]]$learners
  ordered_tables <- lapply(seq_along(members), function(index) {
    member <- members[[index]]
    same_ids <- length(member$pipeline_ids) == length(reference_ids) &&
      setequal(member$pipeline_ids, reference_ids)
    if (!same_ids) {
      sb_abort(
        paste0(
          "Cohort `", cohort_ids[[index]],
          "` does not contain the same pipeline ID set as the first cohort."
        ),
        "stormbreaker_external_pipeline_mismatch"
      )
    }
    positions <- match(reference_ids, member$pipeline_ids)
    if (!identical(member$learners[positions], reference_learners)) {
      sb_abort(
        paste0(
          "Cohort `", cohort_ids[[index]],
          "` does not retain the same learner labels as the first cohort."
        ),
        "stormbreaker_external_pipeline_mismatch"
      )
    }
    member$table[positions, , drop = FALSE]
  })

  task_type <- task_types[[1L]]
  if (identical(task_type, "survival")) {
    reference_unit <- members[[1L]]$time_unit
    reference_horizon <- members[[1L]]$prediction_horizon
    compatible_target <- vapply(members, function(member) {
      same_unit <- identical(member$time_unit, reference_unit)
      same_horizon_presence <- identical(
        is.null(member$prediction_horizon),
        is.null(reference_horizon)
      )
      same_horizon <- same_horizon_presence &&
        (is.null(reference_horizon) || identical(
          as.numeric(member$prediction_horizon),
          as.numeric(reference_horizon)
        ))
      same_unit && same_horizon
    }, logical(1L))
    if (!all(compatible_target)) {
      sb_abort(
        paste0(
          "Survival comparisons must use the same time unit and prediction ",
          "horizon."
        ),
        "stormbreaker_external_target_mismatch"
      )
    }
  }

  results <- lapply(seq_along(ordered_tables), function(index) {
    cbind(
      data.frame(
        cohort_id = rep(cohort_ids[[index]], nrow(ordered_tables[[index]])),
        stringsAsFactors = FALSE
      ),
      ordered_tables[[index]]
    )
  })
  results <- do.call(rbind, results)
  rownames(results) <- NULL

  cohorts <- data.frame(
    cohort_id = cohort_ids,
    status = vapply(members, `[[`, character(1L), "status"),
    n_patients = vapply(members, `[[`, integer(1L), "n_patients"),
    stringsAsFactors = FALSE
  )
  if (identical(task_type, "survival")) {
    cohorts$time_unit <- vapply(members, `[[`, character(1L), "time_unit")
    if (!is.null(members[[1L]]$prediction_horizon)) {
      cohorts$prediction_horizon <- vapply(
        members,
        `[[`,
        numeric(1L),
        "prediction_horizon"
      )
    }
  }

  complete <- all(cohorts$status == "success")
  structure(
    list(
      status = if (complete) "success" else "incomplete",
      task_type = task_type,
      pipeline_ids = reference_ids,
      cohorts = cohorts,
      results = results,
      comparisons = comparisons
    ),
    class = "stormbreaker_external_collection"
  )
}

#' @export
print.stormbreaker_external_collection <- function(x, ...) {
  cat("<stormbreaker_external_collection:", x$status, ">\n")
  cat(
    "Cohorts:", nrow(x$cohorts),
    "| Pipelines:", length(x$pipeline_ids),
    "| Task:", x$task_type, "\n"
  )
  print(x$cohorts, row.names = FALSE)
  cat(
    paste0(
      "No metric is pooled; no pipeline is ranked or selected from these ",
      "cohorts.\n"
    )
  )
  invisible(x)
}

# ---- compare-external ----

sb_external_comparison_value <- function(evaluation, path, default) {
  value <- evaluation
  for (name in path) {
    if (!is.list(value) || is.null(value[[name]])) {
      return(default)
    }
    value <- value[[name]]
  }
  if (length(value) != 1L) {
    return(default)
  }
  value[[1L]]
}

sb_external_binary_comparison_row <- function(evaluation) {
  value <- function(metric, field, default) {
    sb_external_comparison_value(
      evaluation,
      c("metrics", metric, field),
      default
    )
  }
  data.frame(
    pipeline_id = evaluation$pipeline_id,
    learner = evaluation$learner,
    status = evaluation$status,
    patient_match_status = evaluation$patient_match$status,
    roc_auc_status = value("roc_auc", "status", NA_character_),
    roc_auc = value("roc_auc", "estimate", NA_real_),
    comparable_pairs = value("roc_auc", "comparable_pairs", NA_real_),
    average_precision_status = value(
      "average_precision", "status", NA_character_
    ),
    average_precision = value(
      "average_precision", "estimate", NA_real_
    ),
    brier_status = value("brier_score", "status", NA_character_),
    brier_score = value("brier_score", "estimate", NA_real_),
    log_loss_status = value("log_loss", "status", NA_character_),
    log_loss = value("log_loss", "estimate", NA_real_),
    mean_calibration_status = value(
      "mean_calibration", "status", NA_character_
    ),
    mean_calibration_ratio = value(
      "mean_calibration", "mean_calibration_ratio", NA_real_
    ),
    calibration_status = value(
      "calibration", "status", NA_character_
    ),
    calibration_intercept = value(
      "calibration", "calibration_intercept", NA_real_
    ),
    calibration_slope = value(
      "calibration", "calibration_slope", NA_real_
    ),
    moderate_calibration_status = value(
      "moderate_calibration", "status", NA_character_
    ),
    moderate_calibration_ici = value(
      "moderate_calibration", "ici", NA_real_
    ),
    moderate_calibration_e50 = value(
      "moderate_calibration", "e50", NA_real_
    ),
    moderate_calibration_e90 = value(
      "moderate_calibration", "e90", NA_real_
    ),
    n = evaluation$n,
    error_stage = evaluation$error_stage,
    error_class = evaluation$error_class,
    error_message = evaluation$error_message,
    stringsAsFactors = FALSE
  )
}

sb_external_survival_comparison_row <- function(evaluation, has_horizon) {
  value <- function(metric, field, default) {
    sb_external_comparison_value(
      evaluation,
      c("metrics", metric, field),
      default
    )
  }
  row <- data.frame(
    pipeline_id = evaluation$pipeline_id,
    learner = evaluation$learner,
    status = evaluation$status,
    patient_match_status = evaluation$patient_match$status,
    harrell_c_status = value("harrell_c", "status", NA_character_),
    harrell_c = value("harrell_c", "estimate", NA_real_),
    comparable_pairs = value(
      "harrell_c", "comparable_pairs", NA_real_
    ),
    n = evaluation$n,
    error_stage = evaluation$error_stage,
    error_class = evaluation$error_class,
    error_message = evaluation$error_message,
    stringsAsFactors = FALSE
  )
  if (!has_horizon) {
    return(row)
  }

  weak_value <- function(component, field, default) {
    sb_external_comparison_value(
      evaluation,
      c("metrics", "weak_calibration", component, field),
      default
    )
  }
  data.frame(
    pipeline_id = row$pipeline_id,
    learner = row$learner,
    status = row$status,
    patient_match_status = row$patient_match_status,
    harrell_c_status = row$harrell_c_status,
    harrell_c = row$harrell_c,
    comparable_pairs = row$comparable_pairs,
    ipcw_brier_status = value(
      "ipcw_brier", "status", NA_character_
    ),
    ipcw_brier = value("ipcw_brier", "estimate", NA_real_),
    time_auc_status = value("time_auc", "status", NA_character_),
    time_auc = value("time_auc", "estimate", NA_real_),
    time_auc_pair_denominator = value(
      "time_auc", "pair_denominator", NA_real_
    ),
    uno_c_status = value("uno_c", "status", NA_character_),
    uno_c = value("uno_c", "estimate", NA_real_),
    uno_c_weighted_pair_denominator = value(
      "uno_c", "weighted_pair_denominator", NA_real_
    ),
    calibration_status = value(
      "calibration", "status", NA_character_
    ),
    calibration_ici = value("calibration", "ici", NA_real_),
    calibration_e50 = value("calibration", "e50", NA_real_),
    calibration_e90 = value("calibration", "e90", NA_real_),
    mean_calibration_status = weak_value(
      "mean_calibration", "status", NA_character_
    ),
    mean_calibration_ratio = sb_external_comparison_value(
      evaluation,
      c("metrics", "weak_calibration", "mean_calibration_ratio"),
      NA_real_
    ),
    calibration_slope_status = weak_value(
      "slope_fit", "status", NA_character_
    ),
    calibration_slope = sb_external_comparison_value(
      evaluation,
      c("metrics", "weak_calibration", "calibration_slope"),
      NA_real_
    ),
    n = row$n,
    error_stage = row$error_stage,
    error_class = row$error_class,
    error_message = row$error_message,
    stringsAsFactors = FALSE
  )
}

sb_external_prediction_target_state <- function(prediction, endpoint) {
  if (!is.list(prediction)) {
    return(list(known = FALSE))
  }
  task_type <- prediction[["task_type"]]
  if (identical(task_type, "binary")) {
    return(NULL)
  }
  if (!identical(task_type, "survival")) {
    return(list(known = FALSE))
  }

  contract_attempt <- sb_external_attempt(
    sb_external_prediction_contract(prediction)
  )
  if (is.null(contract_attempt$error)) {
    contract <- contract_attempt$value
    return(list(
      known = TRUE,
      has_horizon = contract$has_horizon,
      prediction_horizon = contract$prediction_horizon,
      time_unit = if (contract$has_horizon) {
        contract$time_unit
      } else {
        endpoint$time_unit
      }
    ))
  }

  signature <- prediction[["signature"]]
  if (inherits(signature, "stormbreaker_linear_signature")) {
    signature_attempt <- sb_external_attempt(
      sb_validate_linear_signature(signature)
    )
    if (is.null(signature_attempt$error) &&
        identical(signature_attempt$value$task_type, "survival")) {
      signature <- signature_attempt$value
      has_horizon <- "prediction_horizon" %in% names(signature)
      return(list(
        known = TRUE,
        has_horizon = has_horizon,
        prediction_horizon = if (has_horizon) {
          signature$prediction_horizon
        } else {
          NULL
        },
        time_unit = if (has_horizon) {
          signature$time_unit
        } else {
          endpoint$time_unit
        }
      ))
    }
  }

  rows <- prediction[["predictions"]]
  if (!is.data.frame(rows)) {
    return(list(known = FALSE))
  }
  horizon_fields <- c(
    "prediction_horizon", "time_unit", "survival_probability"
  )
  horizon_counts <- vapply(horizon_fields, function(name) {
    sum(names(rows) == name)
  }, integer(1L))
  if (all(horizon_counts == 0L) && sum(names(rows) == "risk") == 1L) {
    return(list(
      known = TRUE,
      has_horizon = FALSE,
      prediction_horizon = NULL,
      time_unit = endpoint$time_unit
    ))
  }
  if (!all(horizon_counts == 1L) || nrow(rows) == 0L) {
    return(list(known = FALSE))
  }

  horizon <- rows[["prediction_horizon"]]
  time_unit <- rows[["time_unit"]]
  valid_horizon <- is.numeric(horizon) && length(horizon) == nrow(rows) &&
    !anyNA(horizon) && all(is.finite(horizon)) && all(horizon > 0) &&
    length(unique(as.numeric(horizon))) == 1L
  valid_unit <- is.character(time_unit) && length(time_unit) == nrow(rows) &&
    !anyNA(time_unit) && length(unique(time_unit)) == 1L &&
    unique(time_unit) %in% c("day", "month", "year")
  if (!valid_horizon || !valid_unit) {
    return(list(known = FALSE))
  }
  list(
    known = TRUE,
    has_horizon = TRUE,
    prediction_horizon = as.numeric(horizon[[1L]]),
    time_unit = time_unit[[1L]]
  )
}

sb_external_comparison_target <- function(predictions, endpoint) {
  if (identical(endpoint$task_type, "binary")) {
    return(list(
      target_status = "complete",
      has_horizon = FALSE,
      prediction_horizon = NULL,
      time_unit = NULL,
      horizon_source = "not_applicable"
    ))
  }
  if (!is.null(endpoint$prediction_horizon)) {
    return(list(
      target_status = "complete",
      has_horizon = TRUE,
      prediction_horizon = endpoint$prediction_horizon,
      time_unit = endpoint$time_unit,
      horizon_source = "endpoint"
    ))
  }

  states <- lapply(predictions, function(prediction) {
    sb_external_prediction_target_state(prediction, endpoint)
  })
  states <- Filter(Negate(is.null), states)
  known <- Filter(function(state) isTRUE(state$known), states)
  if (length(known) == 0L) {
    return(list(
      target_status = "unknown",
      has_horizon = NA,
      prediction_horizon = NULL,
      time_unit = endpoint$time_unit,
      horizon_source = "unknown"
    ))
  }

  reference <- known[[1L]]
  same_target <- vapply(known, function(target) {
    same_presence <- identical(target$has_horizon, reference$has_horizon)
    same_presence && (!reference$has_horizon || (
      identical(target$time_unit, reference$time_unit) &&
        identical(
          as.numeric(target$prediction_horizon),
          as.numeric(reference$prediction_horizon)
        )
    ))
  }, logical(1L))
  if (!all(same_target)) {
    sb_abort(
      paste0(
        "Survival predictions compared without an endpoint horizon must all ",
        "be risk-only or share one locked time unit and horizon."
      ),
      "stormbreaker_external_target_mismatch"
    )
  }
  reference$known <- NULL
  reference$target_status <- if (length(known) == length(states)) {
    "complete"
  } else {
    "partial"
  }
  reference$horizon_source <- if (reference$has_horizon) {
    "prediction"
  } else {
    "none"
  }
  reference
}

#' @rdname sb_compare_external
#' @export
sb_compare_external <- function(predictions, endpoint) {
  valid_collection <- is.list(predictions) && length(predictions) > 0L &&
    all(vapply(
      predictions,
      function(prediction) {
        inherits(prediction, "stormbreaker_final_prediction")
      },
      logical(1L)
    ))
  if (!valid_collection) {
    sb_abort(
      paste0(
        "`predictions` must be a non-empty list returned by ",
        "`sb_predict_final()` or `sb_predict_signature()`."
      ),
      "stormbreaker_bad_prediction_collection"
    )
  }

  pipeline_ids <- vapply(predictions, function(prediction) {
    pipeline_id <- prediction[["pipeline_id"]]
    if (!is.character(pipeline_id) || length(pipeline_id) != 1L ||
        is.na(pipeline_id) || !nzchar(pipeline_id)) {
      return(NA_character_)
    }
    pipeline_id
  }, character(1L))
  if (anyNA(pipeline_ids)) {
    sb_abort(
      "Every final prediction must retain one non-empty pipeline ID.",
      "stormbreaker_bad_prediction_collection"
    )
  }
  if (anyDuplicated(pipeline_ids) > 0L) {
    sb_abort(
      "Every compared external prediction must have a distinct pipeline ID.",
      "stormbreaker_duplicate_pipeline_id"
    )
  }

  endpoint_contract <- sb_external_endpoint_contract(endpoint)
  prediction_target <- sb_external_comparison_target(
    predictions,
    endpoint_contract
  )
  evaluations <- lapply(predictions, function(prediction) {
    sb_evaluate_external(prediction, endpoint, reverse = FALSE)
  })
  names(evaluations) <- pipeline_ids
  has_horizon <- isTRUE(prediction_target$has_horizon)
  rows <- lapply(evaluations, function(evaluation) {
    if (identical(endpoint_contract$task_type, "binary")) {
      sb_external_binary_comparison_row(evaluation)
    } else {
      sb_external_survival_comparison_row(evaluation, has_horizon)
    }
  })
  comparison <- do.call(rbind, rows)
  rownames(comparison) <- NULL
  comparison$prediction_variation <- vapply(evaluations, function(evaluation) {
    sb_prediction_variation(evaluation$scored_predictions, endpoint_contract$task_type)$status
  }, character(1))

  complete <- all(vapply(
    evaluations,
    function(evaluation) identical(evaluation$status, "success"),
    logical(1L)
  ))
  result <- list(
    status = if (complete) "success" else "incomplete",
    task_type = endpoint_contract$task_type,
    comparison = comparison,
    evaluations = evaluations,
    n_patients = as.integer(length(endpoint_contract$patient_ids)),
    prediction_target = prediction_target,
    endpoint = endpoint_contract$summary
  )
  has_decision_curve <- any(vapply(
    evaluations,
    function(evaluation) {
      !is.null(evaluation$metrics[["decision_curve"]])
    },
    logical(1L)
  ))
  if (has_decision_curve) {
    decision_curves <- lapply(evaluations, function(evaluation) {
      evaluation$metrics[["decision_curve"]]
    })
    names(decision_curves) <- pipeline_ids
    result$decision_curves <- decision_curves
  }
  structure(result, class = "stormbreaker_external_comparison")
}

#' @export
print.stormbreaker_external_comparison <- function(x, ...) {
  cat("<stormbreaker_external_comparison:", x$status, ">\n")
  cat(
    "Pipelines:", nrow(x$comparison),
    "| Patients:", x$n_patients,
    "| Task:", x$task_type, "\n"
  )
  metric_columns <- if (identical(x$task_type, "binary")) {
    c("roc_auc", "average_precision", "brier_score", "log_loss")
  } else {
    columns <- "harrell_c"
    if ("ipcw_brier" %in% names(x$comparison)) {
      columns <- c(columns, "ipcw_brier", "time_auc", "uno_c")
    }
    columns
  }
  print(
    x$comparison[c("pipeline_id", "status", metric_columns)],
    row.names = FALSE,
    digits = 3
  )
  cat(
    paste0(
      "No pipeline is ranked or selected; external outcomes must not ",
      "choose the final model.\n"
    )
  )
  invisible(x)
}
