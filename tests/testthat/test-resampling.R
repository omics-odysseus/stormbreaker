# ---- resampling ----
local({
test_that("mixed-outcome patients stay together in stratified inner folds", {
  groups <- c(rep(paste0("paired", 1:7), each = 2),
              rep("large_negative", 9), paste0("negative", 1:6),
              rep("large_positive", 6), paste0("positive", 1:6))
  outcome <- c(rep(c(0L, 1L), 7), rep(0L, 15), rep(1L, 12))
  withr::local_seed(3749L)
  before <- .Random.seed
  observed <- sb_make_stratified_fold_id(
    outcome, 5L, 3748L, group_ids = groups
  )
  expect_identical(.Random.seed, before)
  expect_identical(observed, sb_make_stratified_fold_id(
    outcome, 5L, 3748L, group_ids = groups
  ))
  expect_identical(unname(vapply(split(observed, groups), function(x) {
    length(unique(x))
  }, integer(1))), rep(1L, length(unique(groups))))
  expect_equal(sort(unique(observed)), 1:5)
  expect_gt(min(table(observed, outcome)), 0)
  few_mixed_y <- c(rep(c(0L, 1L), 3L), rep(0L, 8L), rep(1L, 8L))
  few_mixed_groups <- c(rep(paste0("mixed", 1:3), each = 2L),
    rep(paste0("negative", 1:4), each = 2L),
    rep(paste0("positive", 1:4), each = 2L))
  few_mixed <- sb_make_stratified_fold_id(few_mixed_y, 5L, 3748L,
    group_ids = few_mixed_groups)
  expect_gt(min(table(few_mixed, few_mixed_y)), 0)
  expect_identical(unname(vapply(split(few_mixed, few_mixed_groups), function(x) {
    length(unique(x))
  }, integer(1))), rep(1L, 11L))
  unique_ids <- seq_along(outcome)
  expect_identical(
    sb_make_stratified_fold_id(outcome, 5L, 3748L, group_ids = unique_ids),
    sb_make_stratified_fold_id(outcome, 5L, 3748L)
  )
  error <- capture_sb_error(sb_make_stratified_fold_id(
    rep(c(0L, 1L), each = 10), 5L, 1L,
    group_ids = rep(c("only_negative", "only_positive"), each = 10)
  ))
  expect_s3_class(error, "stormbreaker_infeasible_inner_folds")
})

test_that("task grouping survives validation and protects preassigned folds", {
  d <- sb_example_binary(60L, seed = 3745L)
  d$group <- rep(paste0("patient", 1:30), each = 2)
  d$outcome <- rep(c(0L, 1L), 30)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"),
                         group_col = "group")
  expect_identical(sb_validate_binary_task_contract(task), task)
  subset <- sb_subset_task(task, 1:40)
  expect_identical(subset$group_col, "group")
  expect_identical(subset$data$group, d$group[1:40])
  folds <- sb_make_outer_folds(task, 5L, 3746L)
  expect_identical(folds$source, "generated_group_stratified")
  expect_identical(unname(vapply(split(folds$assignment, d$group), function(x) {
    length(unique(x))
  }, integer(1))), rep(1L, 30))
  supplied <- setNames(folds$assignment, d$patient_id)
  expect_identical(sb_make_outer_folds(task, rev(supplied), 1L)$assignment,
                   folds$assignment)
  supplied[2] <- supplied[1] %% 5L + 1L
  error <- capture_sb_error(sb_make_outer_folds(task, supplied, 1L))
  expect_s3_class(error, "stormbreaker_group_fold_overlap")
  d$group <- paste0("unique", seq_len(nrow(d)))
  unique_task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"),
                                group_col = "group")
  legacy_task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  expect_identical(sb_make_outer_folds(unique_task, 5L, 3746L),
                   sb_make_outer_folds(legacy_task, 5L, 3746L))
})

test_that("outer-fold dispatch preserves the scalar stratified path", {
  d <- make_sb_fixture(n = 30L, seed = 3701L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month"
  )
  expected <- sb_make_stratified_folds(task, folds = 3L, seed = 3702L)

  observed <- sb_make_outer_folds(task, folds = 3L, seed = 3702L)
  named_scalar <- sb_make_outer_folds(
    task,
    folds = c(requested = 3L),
    seed = 3702L
  )

  expect_identical(observed$assignment, expected$assignment)
  expect_identical(observed$splits, expected$splits)
  expect_identical(named_scalar$assignment, expected$assignment)
  expect_identical(named_scalar$splits, expected$splits)
  expect_identical(observed$source, "generated_stratified")
  expect_identical(observed$seed, 3702L)
  expect_identical(observed$n_folds, 3L)
})

test_that("preassigned outer folds align by patient ID rather than position", {
  d <- make_sb_fixture(n = 12L, seed = 3711L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month"
  )
  task_ids <- as.character(task$data[[task$id_col]])
  task_assignment <- rep(1:3, length.out = task$n)
  supplied_order <- c(12L, 2L, 9L, 1L, 7L, 4L, 10L, 3L, 6L, 8L, 5L, 11L)
  supplied <- setNames(
    task_assignment[supplied_order],
    task_ids[supplied_order]
  )

  observed <- sb_make_outer_folds(task, folds = supplied, seed = 999L)
  expected_assignment <- as.integer(supplied[match(task_ids, names(supplied))])

  expect_identical(observed$assignment, expected_assignment)
  expect_identical(observed$source, "preassigned_patient_id")
  expect_identical(observed$seed, NA_integer_)
  expect_identical(observed$n_folds, 3L)
  for (fold in seq_len(3L)) {
    expected_assessment <- which(expected_assignment == fold)
    expected_analysis <- which(expected_assignment != fold)
    expect_identical(observed$splits[[fold]]$fold, fold)
    expect_identical(
      observed$splits[[fold]]$assessment_rows,
      expected_assessment
    )
    expect_identical(observed$splits[[fold]]$analysis_rows, expected_analysis)
    expect_identical(
      observed$splits[[fold]]$assessment_ids,
      task_ids[expected_assessment]
    )
    expect_identical(
      observed$splits[[fold]]$analysis_ids,
      task_ids[expected_analysis]
    )
  }
})

test_that("ambiguous preassigned fold vectors fail before modeling", {
  d <- make_sb_binary_fixture(n = 12L, seed = 3721L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  ids <- as.character(task$data[[task$id_col]])
  valid <- setNames(rep(1:3, length.out = task$n), ids)
  duplicate_name <- valid
  names(duplicate_name)[[2L]] <- names(duplicate_name)[[1L]]
  blank_name <- valid
  names(blank_name)[[2L]] <- ""
  wrong_ids <- valid
  names(wrong_ids)[[2L]] <- "not_a_task_patient"
  missing_value <- valid
  missing_value[[2L]] <- NA_real_
  nonwhole_value <- valid
  nonwhole_value[[2L]] <- 1.5
  nonpositive_value <- valid
  nonpositive_value[[2L]] <- 0
  one_fold <- setNames(rep(1L, task$n), ids)
  nonconsecutive <- valid
  nonconsecutive[nonconsecutive == 3L] <- 4L
  character_value <- setNames(as.character(valid), ids)
  cases <- list(
    unnamed = list(value = unname(valid), pattern = "named"),
    duplicate = list(value = duplicate_name, pattern = "unique"),
    blank = list(value = blank_name, pattern = "non-empty"),
    wrong_ids = list(value = wrong_ids, pattern = "patient IDs"),
    missing = list(value = missing_value, pattern = "finite positive"),
    nonwhole = list(value = nonwhole_value, pattern = "finite positive"),
    nonpositive = list(value = nonpositive_value, pattern = "finite positive"),
    one_fold = list(value = one_fold, pattern = "at least two"),
    nonconsecutive = list(value = nonconsecutive, pattern = "1 through K"),
    character = list(value = character_value, pattern = "numeric")
  )

  for (case in cases) {
    error <- capture_sb_error(
      sb_make_outer_folds(task, folds = case$value, seed = 1L)
    )
    expect_s3_class(error, "stormbreaker_bad_fold_assignment")
    expect_match(conditionMessage(error), case$pattern)
  }
})

test_that("structurally valid outcome-segregated folds are not repaired", {
  d <- make_sb_binary_fixture(n = 20L, seed = 3731L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  ids <- as.character(task$data[[task$id_col]])
  supplied <- setNames(ifelse(d$outcome == 1L, 1L, 2L), ids)

  observed <- sb_make_outer_folds(task, folds = supplied, seed = 3732L)

  expect_identical(observed$assignment, as.integer(supplied))
  expect_identical(
    task$data[[task$outcome_col]][observed$splits[[1L]]$assessment_rows],
    rep(1L, task$n_positive)
  )
  expect_identical(
    task$data[[task$outcome_col]][observed$splits[[2L]]$assessment_rows],
    rep(0L, task$n_negative)
  )
})

test_that("preassigned fold construction does not consume the random stream", {
  d <- make_sb_fixture(n = 18L, seed = 3741L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month"
  )
  ids <- as.character(task$data[[task$id_col]])
  supplied <- setNames(rep(1:3, length.out = task$n), ids)
  set.seed(3742L)
  before <- .Random.seed

  first <- sb_make_outer_folds(task, folds = supplied, seed = 1L)
  after <- .Random.seed
  second <- sb_make_outer_folds(task, folds = supplied, seed = 999L)

  expect_identical(after, before)
  expect_identical(.Random.seed, before)
  expect_identical(first, second)
  expect_identical(first$seed, NA_integer_)
})
})

# ---- run-cv ----
local({
test_that("five-fold CV gives every patient one finite OOF risk", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  result <- sb_run_cv(task, sb_pipeline(), folds = 5, seed = 5201314)

  expect_s3_class(result, "stormbreaker_cv_result")
  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_equal(nrow(result$predictions), nrow(d))
  expect_identical(anyDuplicated(result$predictions$patient_id), 0L)
  expect_setequal(result$predictions$patient_id, d$patient_id)
  expect_equal(result$predictions$time, d$time)
  expect_identical(result$predictions$event, d$event)
  expect_identical(unique(result$predictions$time_unit), "month")
  expect_identical(all(is.finite(result$predictions$risk)), TRUE)
  expect_identical(unique(result$predictions$status), "success")
  expect_identical(nrow(result$runs), 5L)
  expect_identical(unique(result$runs$learner_used), "coxph")
  expect_length(result$coefficients, 5L)
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("success", 5L)
  )
  expect_named(result$metrics, "harrell_c")
  expect_identical(result$metrics$harrell_c$status, "success")
  expect_identical(result$metrics$harrell_c$n_folds, 5L)
  expect_identical(result$metrics$harrell_c$n_evaluable_folds, 5L)
  expect_identical(
    result$metrics$harrell_c$per_fold$fold,
    seq_len(5L)
  )
  expect_identical(
    result$metrics$harrell_c$per_fold$status,
    rep("success", 5L)
  )
  official <- survival::concordance(
    survival::Surv(time, event) ~ risk + survival::strata(fold),
    data = result$predictions,
    reverse = TRUE,
    timewt = "n"
  )
  expect_equal(
    result$metrics$harrell_c$estimate,
    unname(official$concordance),
    tolerance = 1e-15
  )
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(printed, "Fold-stratified Harrell C:")
  expect_match(printed, "Comparable pairs:")
  expect_identical(
    result$runs$n_analysis_events + result$runs$n_assessment_events,
    rep.int(sum(d$event), 5L)
  )
})

test_that("the public result schema includes audited metrics", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  result <- sb_run_cv(task, folds = 2, seed = 2)

  expect_named(
    result,
    c(
      "status", "complete_oof", "predictions", "metrics", "runs",
      "selection_stability", "coefficients", "folds", "task", "pipeline",
      "provenance", "glmnet_active_set_stability"
    ),
    ignore.order = FALSE
  )
  expect_identical(result$selection_stability$status, "not_applicable")
  expect_identical(result$selection_stability$n_folds, 2L)
  expect_false("weak_calibration" %in% names(result$metrics))
  expect_identical(result$selection_stability$n_successful_folds, 0L)
  expect_true(is.na(result$selection_stability$mean_pairwise_jaccard))
  expect_identical(
    result$glmnet_active_set_stability$status,
    "not_applicable"
  )
  expect_true(is.na(result$glmnet_active_set_stability$alpha))
  expect_named(
    result$predictions,
    c(
      "patient_id", "fold", "time", "event", "time_unit", "risk",
      "status", "pipeline_id", "learner"
    ),
    ignore.order = FALSE
  )
  expect_named(
    result$runs,
    c(
      "fold", "seed", "pipeline_id", "selector_requested",
      "learner_requested", "learner_used", "n_analysis",
      "n_analysis_events", "n_assessment", "n_assessment_events",
      "status", "error_class", "error_message"
    ),
    ignore.order = FALSE
  )
})

test_that("a declared horizon adds reconstructible Cox OOF survival probability", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    )
  )

  for (pipeline in pipelines) {
    result <- sb_run_cv(task, pipeline, folds = 5L, seed = 20260819L)

    expect_identical(result$status, "success")
    expect_identical(result$complete_oof, TRUE)
    expect_identical(result$complete_survival_probability_oof, TRUE)
    expect_named(
      result$metrics,
      c(
        "harrell_c", "ipcw_brier", "time_auc", "uno_c", "calibration",
        "weak_calibration"
      )
    )
    expect_identical(result$metrics$ipcw_brier$status, "success")
    expect_identical(result$metrics$ipcw_brier$prediction_horizon, 7)
    expect_true(is.finite(result$metrics$ipcw_brier$estimate))
    expect_identical(
      result$metrics$ipcw_brier$reference$status,
      "success"
    )
    expect_identical(
      result$metrics$ipcw_brier$skill_score$status,
      "success"
    )
    expect_true(is.finite(
      result$metrics$ipcw_brier$skill_score$estimate
    ))
    expect_identical(result$metrics$time_auc$status, "success")
    expect_identical(result$metrics$time_auc$prediction_horizon, 7)
    expect_true(is.finite(result$metrics$time_auc$estimate))
    expect_gt(result$metrics$time_auc$pair_denominator, 0)
    expect_identical(result$metrics$time_auc$n_evaluable_folds, 5L)
    expect_identical(result$metrics$uno_c$status, "success")
    expect_identical(result$metrics$uno_c$prediction_horizon, 7)
    expect_true(is.finite(result$metrics$uno_c$estimate))
    expect_gt(result$metrics$uno_c$weighted_pair_denominator, 0)
    expect_identical(result$metrics$uno_c$n_evaluable_folds, 5L)
    calibration <- result$metrics$calibration
    expect_identical(calibration$status, "success")
    expect_identical(
      calibration$scope,
      "pooled_patient_level_out_of_fold"
    )
    expect_identical(calibration$purpose, "descriptive_evaluation_only")
    expect_identical(calibration$prediction_horizon, 7)
    expect_identical(calibration$time_unit, "month")
    expect_identical(
      calibration$patient_calibration$patient_id,
      result$predictions$patient_id
    )
    direct_calibration <- sb_survival_calibration_metric(
      patient_id = result$predictions$patient_id,
      time = result$predictions$time,
      event = result$predictions$event,
      survival_probability = result$predictions$survival_probability,
      prediction_status = ifelse(
        result$predictions$status == "success" &
          result$predictions$survival_probability_status == "success",
        "success",
        "failed"
      ),
      prediction_horizon = result$task$prediction_horizon,
      time_unit = result$task$time_unit,
      scope = "pooled_patient_level_out_of_fold",
      fold = result$predictions$fold
    )
    expect_identical(calibration$knots, direct_calibration$knots)
    expect_identical(calibration$curve, direct_calibration$curve)
    expect_identical(
      calibration$patient_calibration$absolute_error,
      direct_calibration$patient_calibration$absolute_error
    )
    expect_identical(
      calibration[c("ici", "e50", "e90")],
      direct_calibration[c("ici", "e50", "e90")]
    )
    weak_calibration <- result$metrics$weak_calibration
    direct_weak_calibration <- sb_survival_weak_calibration_metric(
      patient_id = result$predictions$patient_id,
      time = result$predictions$time,
      event = result$predictions$event,
      survival_probability = result$predictions$survival_probability,
      prediction_status = ifelse(
        result$predictions$status == "success" &
          result$predictions$survival_probability_status == "success",
        "success",
        "failed"
      ),
      prediction_horizon = result$task$prediction_horizon,
      time_unit = result$task$time_unit,
      scope = "pooled_patient_level_out_of_fold",
      fold = result$predictions$fold
    )
    expect_identical(weak_calibration$status, "success")
    expect_identical(
      weak_calibration$patient_predictions$patient_id,
      result$predictions$patient_id
    )
    expect_identical(weak_calibration$per_fold$fold, seq_len(5L))
    expect_equal(
      weak_calibration$mean_calibration_ratio,
      direct_weak_calibration$mean_calibration_ratio,
      tolerance = 0
    )
    expect_equal(
      weak_calibration$calibration_slope,
      direct_weak_calibration$calibration_slope,
      tolerance = 0
    )
    expect_named(
      result$predictions,
      c(
        "patient_id", "fold", "time", "event", "time_unit", "risk",
        "prediction_horizon", "survival_probability",
        "survival_probability_status", "status", "pipeline_id", "learner"
      ),
      ignore.order = FALSE
    )
    expect_identical(
      unique(result$predictions$prediction_horizon),
      7
    )
    expect_identical(
      unique(result$predictions$survival_probability_status),
      "success"
    )
    expect_true(all(
      result$predictions$survival_probability >= 0 &
        result$predictions$survival_probability <= 1
    ))
    expect_length(result$baseline_hazards, 5L)

    for (fold in seq_along(result$folds)) {
      rows <- result$predictions$fold == fold
      coefficient_audit <- result$coefficients[[fold]]
      baseline_audit <- result$baseline_hazards[[fold]]
      x <- d[
        match(
          result$predictions$patient_id[rows],
          d$patient_id
        ),
        coefficient_audit$values$feature,
        drop = FALSE
      ]
      reconstructed_risk <- unname(as.vector(
        as.matrix(x) %*% coefficient_audit$values$coefficient
      ))
      reconstructed_probability <- exp(
        -baseline_audit$cumulative_baseline_hazard *
          exp(reconstructed_risk)
      )

      expect_identical(baseline_audit$status, "success")
      expect_identical(baseline_audit$ties, "efron")
      expect_identical(
        baseline_audit$source,
        "coxph_basehaz_zero_reference"
      )
      expect_equal(
        result$predictions$risk[rows],
        reconstructed_risk,
        tolerance = 1e-12
      )
      expect_equal(
        result$predictions$survival_probability[rows],
        reconstructed_probability,
        tolerance = 1e-12
      )
      reference_curve <-
        result$metrics$ipcw_brier$reference$survival_curves[[fold]]
      reference_per_fold <- result$metrics$ipcw_brier$reference$per_fold
      expect_equal(
        reference_per_fold$reference_survival_horizon[[fold]],
        sb_event_survival_at(reference_curve, 7),
        tolerance = 0
      )
    }
    printed <- paste(capture.output(print(result)), collapse = "\n")
    expect_match(printed, "IPCW Brier reference:", fixed = TRUE)
    expect_match(printed, "IPCW Brier skill score:", fixed = TRUE)
    expect_match(printed, "Time-dependent AUC:", fixed = TRUE)
    expect_match(printed, "Fold-stratified Uno C:", fixed = TRUE)
    expect_match(printed, "Survival calibration ICI:", fixed = TRUE)
    expect_match(printed, "Mean calibration O:E:", fixed = TRUE)
    expect_match(printed, "Survival calibration slope:", fixed = TRUE)
    expect_match(printed, "E50:", fixed = TRUE)
    expect_match(printed, "E90:", fixed = TRUE)
  }
})

test_that("survival weak calibration exceptions stay inside the optional metric", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  baseline <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)
  testthat::local_mocked_bindings(
    sb_survival_weak_calibration_metric = function(...) {
      condition <- simpleError("weak calibration exploded")
      class(condition) <- c("test_weak_calibration_error", class(condition))
      stop(condition)
    }
  )

  observed <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)

  expect_identical(observed$status, baseline$status)
  expect_identical(observed$predictions, baseline$predictions)
  expect_identical(observed$folds, baseline$folds)
  expect_identical(observed$coefficients, baseline$coefficients)
  expect_identical(observed$baseline_hazards, baseline$baseline_hazards)
  expect_identical(
    observed$metrics[c(
      "harrell_c", "ipcw_brier", "time_auc", "uno_c", "calibration"
    )],
    baseline$metrics[c(
      "harrell_c", "ipcw_brier", "time_auc", "uno_c", "calibration"
    )]
  )
  expect_identical(
    observed$metrics$weak_calibration$status,
    "metric_failed"
  )
  expect_identical(
    observed$metrics$weak_calibration$error_class,
    "test_weak_calibration_error"
  )
  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "Mean calibration O:E: unavailable", fixed = TRUE)
  expect_match(
    printed,
    "Survival calibration slope: unavailable",
    fixed = TRUE
  )
})

test_that("survival calibration incompleteness cannot alter the CV run", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  baseline <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)
  original_fit <- sb_fit_coxph
  learner_calls <- 0L
  testthat::local_mocked_bindings(
    sb_fit_coxph = function(...) {
      learner_calls <<- learner_calls + 1L
      original_fit(...)
    },
    sb_survival_calibration_metric = function(...) {
      metric <- sb_survival_calibration_empty_result(
        prediction_horizon = 7,
        time_unit = "month",
        scope = "pooled_patient_level_out_of_fold",
        n = nrow(d)
      )
      metric$input_status <- "not_evaluable"
      metric$error_class <- "test_calibration_unavailable"
      metric$error_message <- "deliberately unavailable"
      metric
    }
  )

  observed <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)

  expect_identical(learner_calls, 5L)
  expect_identical(observed$status, baseline$status)
  expect_identical(observed$predictions, baseline$predictions)
  expect_identical(observed$folds, baseline$folds)
  expect_identical(observed$coefficients, baseline$coefficients)
  expect_identical(
    observed$metrics[c("harrell_c", "ipcw_brier", "time_auc", "uno_c")],
    baseline$metrics[c("harrell_c", "ipcw_brier", "time_auc", "uno_c")]
  )
  expect_identical(observed$metrics$calibration$status, "incomplete")
})

test_that("survival calibration exceptions stay inside the optional metric", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  baseline <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)
  testthat::local_mocked_bindings(
    sb_survival_calibration_metric = function(...) {
      condition <- simpleError("calibration exploded")
      class(condition) <- c("test_calibration_error", class(condition))
      stop(condition)
    }
  )

  observed <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)

  expect_identical(observed$status, baseline$status)
  expect_identical(observed$predictions, baseline$predictions)
  expect_identical(observed$folds, baseline$folds)
  expect_identical(observed$coefficients, baseline$coefficients)
  expect_identical(
    observed$metrics[c("harrell_c", "ipcw_brier", "time_auc", "uno_c")],
    baseline$metrics[c("harrell_c", "ipcw_brier", "time_auc", "uno_c")]
  )
  expect_identical(observed$metrics$calibration$status, "metric_failed")
  expect_identical(
    observed$metrics$calibration$error_class,
    "test_calibration_error"
  )
  expect_identical(
    observed$metrics$calibration$error_message,
    "calibration exploded"
  )
  expect_match(
    paste(capture.output(print(observed)), collapse = "\n"),
    "Survival calibration: unavailable",
    fixed = TRUE
  )
})

test_that("horizon CV shares one censoring fit without Brier status coupling", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  original_fit <- sb_fit_censoring_km
  original_brier <- sb_survival_ipcw_brier_metric
  censoring_calls <- 0L
  testthat::local_mocked_bindings(
    sb_fit_censoring_km = function(...) {
      censoring_calls <<- censoring_calls + 1L
      original_fit(...)
    },
    sb_survival_ipcw_brier_metric = function(...) {
      metric <- original_brier(...)
      metric$status <- "incomplete"
      metric$estimate <- NA_real_
      metric
    }
  )

  result <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)

  expect_identical(censoring_calls, 5L)
  expect_identical(result$status, "failed")
  expect_identical(result$metrics$ipcw_brier$status, "incomplete")
  expect_identical(result$metrics$time_auc$status, "success")
  expect_true(is.finite(result$metrics$time_auc$estimate))
  expect_identical(result$metrics$uno_c$status, "success")
  expect_true(is.finite(result$metrics$uno_c$estimate))
})

test_that("glmnet Cox OOF survival probability uses saved Breslow baselines", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L
  )
  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 20260819L)

  expect_identical(result$status, "success")
  expect_identical(result$complete_survival_probability_oof, TRUE)
  for (fold in seq_along(result$folds)) {
    rows <- result$predictions$fold == fold
    coefficient_audit <- result$coefficients[[fold]]
    baseline_audit <- result$baseline_hazards[[fold]]
    x <- d[
      match(result$predictions$patient_id[rows], d$patient_id),
      coefficient_audit$values$feature,
      drop = FALSE
    ]
    reconstructed_risk <- unname(as.vector(
      as.matrix(x) %*% coefficient_audit$values$coefficient
    ))
    reconstructed_probability <- exp(
      -baseline_audit$cumulative_baseline_hazard *
        exp(reconstructed_risk)
    )

    expect_identical(baseline_audit$status, "success")
    expect_identical(baseline_audit$ties, "breslow")
    expect_identical(baseline_audit$source, "glmnet_fixed_beta_breslow")
    expect_equal(
      result$predictions$survival_probability[rows],
      reconstructed_probability,
      tolerance = 1e-12
    )
  }
})

test_that("assessment outcomes cannot change fitted survival probabilities", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  predictors <- c("gene_a", "gene_b")
  make_task <- function(data) {
    sb_surv_task(
      data,
      "patient_id",
      "time",
      "event",
      predictors,
      "month",
      prediction_horizon = 5
    )
  }
  split <- list(
    fold = 1L,
    analysis_rows = seq_len(48L),
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[seq_len(48L)],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline()
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  original <- sb_run_one_fold(
    make_task(d), pipeline, runtime, split, seed = 101L
  )

  changed_outcomes <- d
  changed_outcomes$time[split$assessment_rows] <-
    changed_outcomes$time[split$assessment_rows] + 1000
  changed_outcomes$event[split$assessment_rows] <-
    1L - changed_outcomes$event[split$assessment_rows]
  outcome_perturbed <- sb_run_one_fold(
    make_task(changed_outcomes), pipeline, runtime, split, seed = 101L
  )

  expect_identical(
    original$baseline_hazard,
    outcome_perturbed$baseline_hazard
  )
  expect_identical(original$coefficients, outcome_perturbed$coefficients)
  expect_equal(
    original$predictions$survival_probability,
    outcome_perturbed$predictions$survival_probability,
    tolerance = 0
  )
  analysis_inputs <- list(list(
    fold = 1L,
    time = d$time[split$analysis_rows],
    event = d$event[split$analysis_rows]
  ))
  original_metric <- sb_survival_ipcw_brier_metric(
    original$predictions,
    analysis_inputs,
    prediction_horizon = 5
  )
  perturbed_metric <- sb_survival_ipcw_brier_metric(
    outcome_perturbed$predictions,
    analysis_inputs,
    prediction_horizon = 5
  )
  expect_identical(
    original_metric$censoring_curves,
    perturbed_metric$censoring_curves
  )
  expect_identical(
    original_metric$reference$survival_curves,
    perturbed_metric$reference$survival_curves
  )
  expect_identical(
    original_metric$reference$per_fold$reference_survival_horizon,
    perturbed_metric$reference$per_fold$reference_survival_horizon
  )
  expect_gt(
    abs(
      original_metric$weighted_loss_sum -
        perturbed_metric$weighted_loss_sum
    ),
    0
  )
  expect_gt(
    abs(
      original_metric$reference$weighted_loss_sum -
        perturbed_metric$reference$weighted_loss_sum
    ),
    0
  )

  changed_predictors <- d
  changed_predictors$gene_a[split$assessment_rows] <-
    changed_predictors$gene_a[split$assessment_rows] + 5
  predictor_perturbed <- sb_run_one_fold(
    make_task(changed_predictors), pipeline, runtime, split, seed = 101L
  )

  expect_identical(original$baseline_hazard, predictor_perturbed$baseline_hazard)
  expect_identical(original$coefficients, predictor_perturbed$coefficients)
  expect_gt(
    max(abs(
      original$predictions$survival_probability -
        predictor_perturbed$predictions$survival_probability
    )),
    0
  )
  predictor_metric <- sb_survival_ipcw_brier_metric(
    predictor_perturbed$predictions,
    analysis_inputs,
    prediction_horizon = 5
  )
  expect_identical(original_metric$censoring_curves, predictor_metric$censoring_curves)
  expect_identical(original_metric$reference, predictor_metric$reference)
  expect_gt(
    abs(
      original_metric$weighted_loss_sum -
        predictor_metric$weighted_loss_sum
    ),
    0
  )
})

test_that("survival-probability failure preserves successful risk evidence", {
  d <- make_sb_fixture(n = 30L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 5
  )
  split <- list(
    fold = 1L,
    analysis_rows = 1:20,
    assessment_rows = 21:30,
    analysis_ids = d$patient_id[1:20],
    assessment_ids = d$patient_id[21:30]
  )
  pipeline <- sb_pipeline()
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  runtime$learner$predict_survival <- function(...) {
    sb_abort(
      "Synthetic survival-probability failure.",
      "stormbreaker_test_survival_probability_failure"
    )
  }

  result <- sb_run_one_fold(task, pipeline, runtime, split, seed = 1L)

  expect_identical(result$run$status, "success")
  expect_true(all(is.finite(result$predictions$risk)))
  expect_identical(unique(result$predictions$status), "success")
  expect_true(all(is.na(result$predictions$survival_probability)))
  expect_identical(
    unique(result$predictions$survival_probability_status),
    "failed"
  )
  expect_identical(result$coefficients$status, "success")
  expect_identical(result$baseline_hazard$status, "failed")
  expect_identical(
    result$baseline_hazard$error_class,
    "stormbreaker_test_survival_probability_failure"
  )
  metric <- sb_survival_ipcw_brier_metric(
    result$predictions,
    list(list(
      fold = 1L,
      time = d$time[split$analysis_rows],
      event = d$event[split$analysis_rows]
    )),
    prediction_horizon = 5
  )
  expect_identical(metric$status, "incomplete")
  expect_identical(metric$reference$status, "success")
  expect_identical(metric$reference$n_evaluable_folds, 1L)
  expect_identical(metric$skill_score$status, "incomplete")
  expect_identical(is.na(metric$skill_score$estimate), TRUE)
  uno <- sb_survival_uno_c_metric(
    result$predictions,
    list(list(
      fold = 1L,
      time = d$time[split$analysis_rows],
      event = d$event[split$analysis_rows]
    )),
    prediction_horizon = 5,
    censoring_curves = metric$censoring_curves
  )
  expect_identical(uno$status, "success")
  expect_true(is.finite(uno$estimate))
})

test_that("an unavailable requested Uno C preserves all other CV evidence", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  original_uno <- sb_survival_uno_c_metric
  testthat::local_mocked_bindings(
    sb_survival_uno_c_metric = function(...) {
      metric <- original_uno(...)
      metric$status <- "incomplete"
      metric$estimate <- NA_real_
      metric
    }
  )

  result <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(result$complete_survival_probability_oof, TRUE)
  expect_identical(result$metrics$harrell_c$status, "success")
  expect_identical(result$metrics$ipcw_brier$status, "success")
  expect_identical(result$metrics$time_auc$status, "success")
  expect_identical(result$metrics$uno_c$status, "incomplete")
  expect_true(is.na(result$metrics$uno_c$estimate))
})

test_that("analysis and assessment patients never overlap", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )
  result <- sb_run_cv(task, folds = 5, seed = 17)

  for (split in result$folds) {
    expect_length(
      intersect(split$analysis_ids, split$assessment_ids),
      0L
    )
    expect_setequal(
      c(split$analysis_ids, split$assessment_ids),
      d$patient_id
    )
  }
})

test_that("OOF risk is reproduced by fitting only the recorded analysis rows", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )
  result <- sb_run_cv(task, folds = 5, seed = 88)

  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    audit <- result$coefficients[[fold]]
    analysis <- d[d$patient_id %in% split$analysis_ids, , drop = FALSE]
    assessment <- d[d$patient_id %in% split$assessment_ids, , drop = FALSE]
    manual <- survival::coxph(
      survival::Surv(time, event) ~ gene_a + gene_b,
      data = analysis,
      ties = "efron",
      singular.ok = FALSE,
      x = TRUE,
      y = TRUE
    )
    expected <- unname(stats::predict(
      manual,
      newdata = assessment,
      type = "lp",
      reference = "zero"
    ))
    reconstructed <- unname(as.vector(
      as.matrix(assessment[audit$values$feature]) %*%
        audit$values$coefficient
    ))
    observed <- result$predictions$risk[
      match(assessment$patient_id, result$predictions$patient_id)
    ]
    expect_identical(audit$fold, as.integer(fold))
    expect_identical(audit$learner, "coxph")
    expect_identical(audit$status, "success")
    expect_true(is.na(audit$lambda))
    expect_true(is.na(audit$alpha))
    expect_identical(audit$values$feature, c("gene_a", "gene_b"))
    expect_equal(
      audit$values$coefficient,
      unname(stats::coef(manual)[audit$values$feature]),
      tolerance = 1e-15
    )
    expect_identical(audit$n_nonzero, as.integer(sum(audit$values$nonzero)))
    expect_identical(
      audit$nonzero_features,
      audit$values$feature[audit$values$nonzero]
    )
    expect_equal(observed, expected, tolerance = 1e-12)
    expect_equal(observed, reconstructed, tolerance = 1e-12)
  }
})

test_that("CV declares the common zero-reference linear risk scale", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  result <- sb_run_cv(task, folds = 5, seed = 88)

  expect_identical(result$pipeline$risk_type, "linear_predictor")
  expect_identical(result$pipeline$risk_reference, "zero")
  expect_identical(result$pipeline$risk_formula, "x_beta")
  expect_identical(
    result$pipeline$risk_direction,
    "higher_is_higher_hazard"
  )
})

test_that("the same seed reproduces folds and OOF risks", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )
  first <- sb_run_cv(task, folds = 5, seed = 101)
  second <- sb_run_cv(task, folds = 5, seed = 101)

  expect_identical(first$predictions, second$predictions)
  expect_identical(first$metrics, second$metrics)
  expect_identical(first$runs, second$runs)
  expect_identical(first$folds, second$folds)
})

test_that("fold assignment does not change the caller's random stream", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )

  set.seed(404)
  expected <- stats::runif(3)
  set.seed(404)
  sb_run_cv(task, folds = 5, seed = 101)
  observed <- stats::runif(3)

  expect_identical(observed, expected)
})

test_that("failed Cox fits remain visible and never use a fallback", {
  d <- make_sb_fixture(n = 60)
  d$gene_b <- d$gene_a
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )
  result <- sb_run_cv(task, folds = 5, seed = 9)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(unique(result$runs$learner_requested), "coxph")
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(all(nzchar(result$runs$error_message)), TRUE)
  expect_identical(all(is.na(result$predictions$risk)), TRUE)
  expect_identical(result$metrics$harrell_c$status, "incomplete")
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_fitted", 5L)
  )
  expect_identical(
    vapply(result$coefficients, function(audit) nrow(audit$values), integer(1)),
    rep(0L, 5L)
  )
  expect_identical(is.na(result$metrics$harrell_c$estimate), TRUE)
  expect_identical(
    result$metrics$harrell_c$per_fold$status,
    rep("prediction_failed", 5L)
  )
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(
    printed,
    "Harrell C: unavailable (incomplete)",
    fixed = TRUE
  )
  expect_equal(result$predictions$time, d$time)
  expect_identical(result$predictions$event, d$event)
  expect_identical(unique(result$predictions$status), "failed")
})

test_that("folds require enough events and censored patients", {
  d <- make_sb_fixture(n = 20)
  d$event <- c(rep(1L, 4), rep(0L, 16))
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  error <- capture_sb_error(sb_run_cv(task, folds = 5, seed = 1))

  expect_s3_class(error, "stormbreaker_infeasible_folds")
})

test_that("logical event coding runs through stratified CV", {
  d <- make_sb_fixture(n = 20)
  d$event <- as.logical(d$event)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  result <- sb_run_cv(task, folds = 2, seed = 1)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_setequal(unique(result$predictions$event), c(0L, 1L))
})

test_that("run arguments must use STORMBREAKER contracts", {
  d <- make_sb_fixture(n = 20)
  task_error <- capture_sb_error(sb_run_cv(d, folds = 2, seed = 1))
  expect_s3_class(task_error, "stormbreaker_bad_task")

  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  pipeline_error <- capture_sb_error(
    sb_run_cv(task, pipeline = list(learner = "coxph"), folds = 2, seed = 1)
  )
  expect_s3_class(pipeline_error, "stormbreaker_bad_pipeline")
})

test_that("mutated pipeline metadata cannot masquerade as the core pipeline", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  pipeline <- sb_pipeline()
  pipeline$id <- "forged_pipeline"
  error <- capture_sb_error(
    sb_run_cv(task, pipeline = pipeline, folds = 2, seed = 1)
  )

  expect_s3_class(error, "stormbreaker_unresolved_pipeline")
})

test_that("mutated risk scale metadata cannot masquerade as the core pipeline", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  pipeline <- sb_pipeline()
  pipeline$risk_reference <- "sample"
  error <- capture_sb_error(
    sb_run_cv(task, pipeline = pipeline, folds = 2, seed = 1)
  )

  expect_s3_class(error, "stormbreaker_unresolved_pipeline")
})

test_that("mutated task metadata is rejected before resampling", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  task$n <- task$n - 1L
  error <- capture_sb_error(sb_run_cv(task, folds = 2, seed = 1))

  expect_s3_class(error, "stormbreaker_mutated_task")
})

test_that("Top-K CV stores one complete selection audit per fold", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    selector_n_features = 2L
  )
  result <- sb_run_cv(task, pipeline, folds = 5, seed = 101)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_length(result$selections, 5L)
  expect_length(result$coefficients, 5L)
  expect_identical(result$metrics$harrell_c$status, "success")
  for (selection in result$selections) {
    expect_identical(selection$status, "success")
    expect_length(selection$selected_features, 2L)
    expect_identical(nrow(selection$scores), 3L)
  }
  for (fold in seq_along(result$folds)) {
    selection <- result$selections[[fold]]
    audit <- result$coefficients[[fold]]
    assessment_rows <- result$folds[[fold]]$assessment_rows
    reconstructed <- unname(as.vector(
      as.matrix(d[
        assessment_rows,
        audit$values$feature,
        drop = FALSE
      ]) %*% audit$values$coefficient
    ))
    observed <- result$predictions$risk[
      match(
        result$folds[[fold]]$assessment_ids,
        result$predictions$patient_id
      )
    ]

    expect_identical(audit$status, "success")
    expect_identical(audit$values$feature, selection$selected_features)
    expect_equal(reconstructed, observed, tolerance = 1e-12)
  }
  stability <- result$selection_stability
  expected_count <- vapply(task$predictors, function(feature) {
    as.integer(sum(vapply(
      result$selections,
      function(selection) feature %in% selection$selected_features,
      logical(1)
    )))
  }, integer(1))
  expected_jaccard <- vapply(
    seq_len(nrow(stability$pairwise)),
    function(row) {
      selected_i <- result$selections[[
        stability$pairwise$fold_i[[row]]
      ]]$selected_features
      selected_j <- result$selections[[
        stability$pairwise$fold_j[[row]]
      ]]$selected_features
      length(intersect(selected_i, selected_j)) /
        length(union(selected_i, selected_j))
    },
    numeric(1)
  )
  expect_identical(stability$status, "success")
  expect_identical(stability$n_successful_folds, 5L)
  expect_identical(stability$feature_frequency$feature, task$predictors)
  expect_identical(
    stability$feature_frequency$selected_count,
    unname(expected_count)
  )
  expect_equal(stability$pairwise$jaccard, expected_jaccard, tolerance = 0)
  expect_equal(
    stability$mean_pairwise_jaccard,
    mean(expected_jaccard),
    tolerance = 0
  )
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(printed, "Selection stability", fixed = TRUE)
})

test_that("failed survival selector folds keep stability incomplete", {
  d <- make_sb_fixture(n = 40L, seed = 20260819L)
  d$constant_a <- 1
  d$constant_b <- 2
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("constant_a", "constant_b"),
    "month"
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    selector_n_features = 1L
  )

  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 20260819L)
  stability <- result$selection_stability

  expect_identical(result$status, "failed")
  expect_identical(
    vapply(result$selections, `[[`, "", "status"),
    rep("failed", 4L)
  )
  expect_identical(stability$status, "incomplete")
  expect_identical(stability$n_successful_folds, 0L)
  expect_identical(stability$n_evaluable_pairs, 0L)
  expect_true(is.na(stability$mean_pairwise_jaccard))
  expect_identical(stability$feature_frequency$selected_count, c(0L, 0L))
  expect_true(all(is.na(stability$feature_frequency$selection_frequency)))
  expect_true(all(stability$per_fold$status == "failed"))
  expect_true(all(stability$pairwise$status == "incomplete"))
})

test_that("nested glmnet Cox tunes once inside every outer analysis fold", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260813)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1701L
  )
  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 101L)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(result$metrics$harrell_c$status, "success")
  expect_identical(unique(result$runs$learner_used), "glmnet_cox")
  expect_length(result$tuning, 5L)
  expect_identical(
    vapply(result$tuning, function(x) x$status, ""),
    rep("success", 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) x$selection_rule, ""),
    rep("lambda.min", 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) x$inner_folds, integer(1)),
    rep(3L, 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) is.finite(x$lambda), logical(1)),
    rep(TRUE, 5L)
  )
  for (fold in seq_along(result$tuning)) {
    expect_setequal(
      names(result$tuning[[fold]]$fold_id),
      result$folds[[fold]]$analysis_ids
    )
  }
  expect_identical(
    result$provenance$glmnet_version,
    as.character(utils::packageVersion("glmnet"))
  )
})

test_that("Elastic Net settings are preserved in every outer-fold audit", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260817)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 1701L,
    glmnet_standardize = FALSE,
    glmnet_lambda_rule = "lambda.1se",
    glmnet_tuning_measure = "C"
  )

  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 101L)

  expect_identical(result$status, "success")
  expect_identical(result$pipeline$learner_params$alpha, 0.5)
  expect_identical(result$pipeline$learner_params$standardize, FALSE)
  expect_identical(
    result$pipeline$learner_params$selection_rule,
    "lambda.1se"
  )
  expect_identical(
    vapply(result$tuning, function(x) x$alpha, numeric(1)),
    rep(0.5, 5L)
  )
  expect_identical(
    vapply(result$coefficients, function(x) x$alpha, numeric(1)),
    rep(0.5, 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) x$standardize, logical(1)),
    rep(FALSE, 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) x$selection_rule, character(1)),
    rep("lambda.1se", 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) x$measure, character(1)),
    rep("C", 5L)
  )
})

test_that("Ridge alpha zero reaches tuning and final coefficient audits", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260818)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0,
    inner_folds = 3L,
    inner_seed = 1701L
  )

  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 101L)

  expect_identical(result$status, "success")
  expect_identical(
    vapply(result$tuning, function(x) x$alpha, numeric(1)),
    rep(0, 5L)
  )
  expect_identical(
    vapply(result$coefficients, function(x) x$alpha, numeric(1)),
    rep(0, 5L)
  )
})

test_that("nested glmnet Cox records outer-fold coefficient audits", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260813)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month"
  )
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1701L
  )
  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 101L)

  expect_length(result$coefficients, 5L)
  for (fold in seq_along(result$coefficients)) {
    audit <- result$coefficients[[fold]]
    expect_identical(audit$fold, as.integer(fold))
    expect_identical(audit$learner, "glmnet_cox")
    expect_identical(audit$status, "success")
    expect_equal(audit$lambda, result$tuning[[fold]]$lambda, tolerance = 0)
    expect_identical(audit$values$feature, predictors)
    expect_identical(nrow(audit$values), length(predictors))
    expect_identical(
      audit$n_nonzero,
      as.integer(sum(audit$values$nonzero))
    )
    expect_identical(
      audit$nonzero_features,
      audit$values$feature[audit$values$nonzero]
    )
    assessment <- d[
      result$folds[[fold]]$assessment_rows,
      predictors,
      drop = FALSE
    ]
    reconstructed_risk <- unname(as.vector(
      as.matrix(assessment) %*% audit$values$coefficient
    ))
    observed_risk <- result$predictions$risk[
      match(
        result$folds[[fold]]$assessment_ids,
        result$predictions$patient_id
      )
    ]
    expect_equal(reconstructed_risk, observed_risk, tolerance = 1e-12)
  }
})

test_that("outer assessment patients cannot change fitted coefficients", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260813)
  predictors <- c("gene_a", "gene_b", "gene_c")
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1701L
  )
  make_task <- function(data) {
    sb_surv_task(
      data, "patient_id", "time", "event", predictors, "month"
    )
  }
  original_task <- make_task(d)
  split <- sb_make_stratified_folds(
    original_task,
    folds = 5L,
    seed = 101L
  )$splits[[1]]
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  original <- sb_run_one_fold(
    original_task, pipeline, runtime, split, seed = 101L
  )

  changed <- d
  rows <- split$assessment_rows
  changed[rows, predictors] <- changed[rows, predictors] * 10000
  changed$time[rows] <- changed$time[rows] + 10000
  changed$event[rows] <- 1L - changed$event[rows]
  perturbed <- sb_run_one_fold(
    make_task(changed), pipeline, runtime, split, seed = 101L
  )

  expect_identical(original$coefficients, perturbed$coefficients)
})

test_that("a final glmnet fit failure cannot invent coefficients", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260813)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1701L
  )
  original_adapter <- sb_glmnet_cox_learner_adapter
  local_mocked_bindings(
    sb_glmnet_cox_learner_adapter = function(...) {
      adapter <- original_adapter(...)
      adapter$fit <- function(x, time, event) {
        sb_abort(
          "Synthetic final glmnet fit failure.",
          "stormbreaker_test_glmnet_final_fit_failure"
        )
      }
      adapter
    }
  )

  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 101L)

  expect_identical(
    vapply(result$tuning, function(x) x$status, ""),
    rep("success", 5L)
  )
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_test_glmnet_final_fit_failure"
  )
  expect_identical(
    vapply(result$coefficients, function(x) x$status, ""),
    rep("not_fitted", 5L)
  )
  expect_identical(all(is.na(result$predictions$risk)), TRUE)
})

test_that("each nested lambda ignores its outer assessment patients", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260813)
  predictors <- c("gene_a", "gene_b", "gene_c")
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1701L
  )
  make_task <- function(data) {
    sb_surv_task(
      data, "patient_id", "time", "event", predictors, "month"
    )
  }
  original_task <- make_task(d)
  split <- sb_make_stratified_folds(
    original_task,
    folds = 5L,
    seed = 101L
  )$splits[[1]]
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  original <- sb_run_one_fold(
    original_task,
    pipeline,
    runtime,
    split,
    seed = 101L
  )

  changed <- d
  rows <- split$assessment_rows
  changed[rows, predictors] <- changed[rows, predictors] * 10000
  changed$time[rows] <- changed$time[rows] + 10000
  changed$event[rows] <- 1L - changed$event[rows]
  perturbed <- sb_run_one_fold(
    make_task(changed),
    pipeline,
    runtime,
    split,
    seed = 101L
  )

  expect_equal(
    original$tuning$lambda,
    perturbed$tuning$lambda,
    tolerance = 0
  )
  expect_identical(
    original$tuning$fold_id,
    perturbed$tuning$fold_id
  )
  expect_identical(
    original$tuning$curve,
    perturbed$tuning$curve
  )
})

test_that("nested tuning failures remain visible without fallback", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260813)
  d$constant_a <- 1
  d$constant_b <- 2
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("constant_a", "constant_b"), "month"
  )
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1701L
  )
  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 101L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_glmnet_tuning_failure"
  )
  expect_identical(all(is.na(result$predictions$risk)), TRUE)
  expect_identical(
    vapply(result$tuning, function(x) x$status, ""),
    rep("failed", 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) x$error_class, ""),
    rep("stormbreaker_glmnet_tuning_failure", 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) nzchar(x$error_message), logical(1)),
    rep(TRUE, 5L)
  )
  expect_identical(
    vapply(result$coefficients, function(x) x$status, ""),
    rep("not_fitted", 5L)
  )
  expect_identical(
    vapply(result$coefficients, function(x) nrow(x$values), integer(1)),
    rep(0L, 5L)
  )
})

test_that("failed tail-warning validation retains the original warning", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260816)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1701L
  )
  tail_warning <- paste0(
    "from glmnet C++ code (error code -27); Convergence for 27th ",
    "lambda value not reached after maxit=100000 iterations; ",
    "solutions for larger lambdas returned"
  )
  local_mocked_bindings(
    cv.glmnet = function(...) {
      warning(tail_warning, call. = FALSE)
      list(
        lambda = c(0.2, 0.1),
        cvm = c(2, 1),
        cvsd = c(0.1, 0.1),
        nzero = c(0L, 1L),
        lambda.min = 0.1
      )
    },
    glmnet = function(...) {
      structure(list(lambda = 0.5), class = "glmnet")
    },
    .package = "glmnet"
  )

  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 101L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_glmnet_selected_lambda_not_converged"
  )
  expect_identical(
    vapply(
      result$tuning,
      function(x) identical(x$warnings, tail_warning),
      logical(1)
    ),
    rep(TRUE, 3L)
  )
  expect_identical(
    vapply(result$tuning, function(x) is.null(x$warning_validation),
           logical(1)),
    rep(TRUE, 3L)
  )
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(all(is.na(result$predictions$risk)), TRUE)
})

test_that("survival CV DCA is reconstructed from fixed OOF probabilities", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  thresholds <- c(0.20, 0.40, 0.60)
  plain_task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month",
    prediction_horizon = 7
  )
  threshold_task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month",
    prediction_horizon = 7,
    decision_thresholds = thresholds
  )
  pipeline <- sb_pipeline()
  plain <- sb_run_cv(plain_task, pipeline, folds = 5L, seed = 20260819L)
  observed <- sb_run_cv(
    threshold_task,
    pipeline,
    folds = 5L,
    seed = 20260819L
  )
  metric <- observed$metrics$decision_curve

  expect_identical(observed$status, "success")
  expect_identical(metric$status, "success")
  expect_identical(metric$scope, "pooled_patient_level_out_of_fold")
  expect_identical(metric$curve$threshold_probability, thresholds)
  expect_identical(observed$task$decision_thresholds, thresholds)
  expect_identical(observed$predictions, plain$predictions)
  expect_identical(observed$folds, plain$folds)
  expect_false("decision_curve" %in% names(plain$metrics))
  expect_false(any(c(
    "decision", "action", "treat", "selected"
  ) %in% names(observed$predictions)))

  threshold <- thresholds[[1L]]
  selected <- 1 - observed$predictions$survival_probability >= threshold
  manual_survival <- function(time, event, horizon) {
    event_times <- sort(unique(time[event == 1L & time <= horizon]))
    if (length(event_times) == 0L) {
      return(1)
    }
    prod(vapply(event_times, function(value) {
      1 - sum(time == value & event == 1L) / sum(time >= value)
    }, numeric(1)))
  }
  selected_survival <- manual_survival(
    observed$predictions$time[selected],
    observed$predictions$event[selected],
    7
  )
  overall_survival <- manual_survival(
    observed$predictions$time,
    observed$predictions$event,
    7
  )
  selected_fraction <- sum(selected) / nrow(observed$predictions)
  odds <- threshold / (1 - threshold)
  expected_model <- selected_fraction * (1 - selected_survival) -
    selected_fraction * selected_survival * odds
  expected_all <- (1 - overall_survival) - overall_survival * odds
  row <- metric$curve[1L, ]

  expect_equal(row$selected_survival_probability, selected_survival)
  expect_equal(row$model_net_benefit, expected_model)
  expect_equal(row$treat_all_net_benefit, expected_all)
  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "OOF survival decision curve: 3 declared thresholds")
})

test_that("survival CV reuses an exact patient-keyed outer assignment", {
  d <- make_sb_fixture(n = 60L, seed = 6101L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 6
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "coxph",
    selector_n_features = 2L
  )
  generated <- sb_run_cv(task, pipeline, folds = 5L, seed = 6102L)
  assignment <- setNames(
    generated$predictions$fold,
    generated$predictions$patient_id
  )
  supplied <- assignment[rev(seq_along(assignment))]

  reused <- sb_run_cv(task, pipeline, folds = supplied, seed = 999L)
  task_ids <- as.character(task$data[[task$id_col]])
  expected_assignment <- data.frame(
    patient_id = task_ids,
    fold = as.integer(supplied[match(task_ids, names(supplied))]),
    stringsAsFactors = FALSE
  )

  expect_identical(reused$folds, generated$folds)
  expect_identical(reused$predictions, generated$predictions)
  expect_identical(reused$metrics, generated$metrics)
  expect_identical(reused$coefficients, generated$coefficients)
  expect_identical(reused$baseline_hazards, generated$baseline_hazards)
  expect_identical(reused$selections, generated$selections)
  expect_identical(
    reused$selection_stability,
    generated$selection_stability
  )
  expect_identical(
    reused$runs[setdiff(names(reused$runs), "seed")],
    generated$runs[setdiff(names(generated$runs), "seed")]
  )
  expect_identical(generated$runs$seed, rep(6102L, 5L))
  expect_identical(reused$runs$seed, rep(NA_integer_, 5L))
  expect_identical(generated$provenance$fold_source, "generated_stratified")
  expect_identical(reused$provenance$fold_source, "preassigned_patient_id")
  expect_identical(generated$provenance$seed, 6102L)
  expect_identical(reused$provenance$seed, NA_integer_)
  expect_identical(generated$provenance$folds, 5L)
  expect_identical(reused$provenance$folds, 5L)
  expect_identical(reused$provenance$fold_assignment, expected_assignment)
  expect_identical(
    generated$provenance$fold_assignment,
    expected_assignment
  )
  printed <- paste(capture.output(print(reused)), collapse = "\n")
  expect_match(printed, "preassigned patient-ID folds", fixed = TRUE)
})

test_that("preassigned survival folds keep glmnet inner tuning in analysis", {
  skip_if_not_installed("glmnet")
  d <- make_sb_fixture(n = 60L, seed = 6111L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 6112L
  )
  generated_folds <- sb_make_stratified_folds(task, 3L, 6113L)
  ids <- as.character(task$data[[task$id_col]])
  supplied <- setNames(generated_folds$assignment, ids)

  result <- sb_run_cv(task, pipeline, folds = supplied, seed = 999L)

  expect_identical(result$status, "success")
  expect_identical(result$provenance$fold_source, "preassigned_patient_id")
  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    inner_ids <- names(result$tuning[[fold]]$fold_id)
    expect_setequal(inner_ids, split$analysis_ids)
    expect_length(intersect(inner_ids, split$assessment_ids), 0L)
  }
})
})

# ---- run-binary-cv ----
local({
test_that("binary CV gives every patient one finite OOF probability", {
  d <- make_sb_binary_fixture(n = 60)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )
  result <- sb_run_cv(task, folds = 5, seed = 20260818)

  expect_s3_class(result, "stormbreaker_binary_cv_result")
  expect_s3_class(result, "stormbreaker_cv_result")
  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(result$pipeline$learner, "logistic")
  expect_identical(result$pipeline$task_type, "binary")
  expect_identical(nrow(result$predictions), nrow(d))
  expect_identical(anyDuplicated(result$predictions$patient_id), 0L)
  expect_equal(result$predictions$outcome, d$outcome)
  expect_identical(
    all(result$predictions$probability >= 0 &
          result$predictions$probability <= 1),
    TRUE
  )
  expect_named(
    result$metrics,
    c(
      "roc_auc", "average_precision", "brier_score", "log_loss",
      "calibration", "moderate_calibration"
    )
  )
  expect_identical(result$metrics$roc_auc$status, "success")
  expect_identical(result$metrics$average_precision$status, "success")
  expect_equal(
    result$metrics$average_precision$estimate,
    sum(
      result$metrics$average_precision$per_fold$n_positive *
        result$metrics$average_precision$per_fold$estimate
    ) / result$metrics$average_precision$n_positive,
    tolerance = 1e-15
  )
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(
    printed,
    "Fold-stratified average precision:",
    fixed = TRUE
  )
  expect_identical(result$metrics$brier_score$status, "success")
  expect_equal(
    result$metrics$brier_score$estimate,
    mean(
      (result$predictions$outcome - result$predictions$probability)^2
    ),
    tolerance = 0
  )
  case_loss <- ifelse(
    result$predictions$outcome == 1L,
    -log(result$predictions$probability),
    -log1p(-result$predictions$probability)
  )
  expect_identical(result$metrics$log_loss$status, "success")
  expect_equal(
    result$metrics$log_loss$estimate,
    mean(case_loss),
    tolerance = 0
  )
  analysis_prevalence <- result$runs$n_analysis_positive /
    result$runs$n_analysis
  reference_probability <- analysis_prevalence[result$predictions$fold]
  reference_case_loss <- ifelse(
    result$predictions$outcome == 1L,
    -log(reference_probability),
    -log1p(-reference_probability)
  )
  expect_equal(
    result$metrics$log_loss$reference$estimate,
    mean(reference_case_loss),
    tolerance = 0
  )
  expect_equal(
    result$metrics$log_loss$skill_score$estimate,
    1 - sum(case_loss) / sum(reference_case_loss),
    tolerance = 0
  )
  expect_match(printed, "OOF log loss:", fixed = TRUE)
  expect_match(
    printed,
    "Log-loss skill vs analysis prevalence:",
    fixed = TRUE
  )
  logit_probability <- stats::qlogis(result$predictions$probability)
  direct_calibration_intercept <- stats::glm(
    result$predictions$outcome ~ 1,
    family = stats::binomial(),
    offset = logit_probability
  )
  direct_calibration_slope <- stats::glm(
    result$predictions$outcome ~ logit_probability,
    family = stats::binomial()
  )
  expect_identical(result$metrics$calibration$status, "success")
  expect_equal(
    result$metrics$calibration$calibration_intercept,
    unname(stats::coef(direct_calibration_intercept)[[1L]]),
    tolerance = 1e-12
  )
  expect_equal(
    result$metrics$calibration$calibration_slope,
    unname(stats::coef(direct_calibration_slope)[[2L]]),
    tolerance = 1e-12
  )
  expect_match(printed, "OOF calibration intercept:", fixed = TRUE)
  direct_moderate_calibration <- sb_binary_moderate_calibration_metric(
    result$predictions,
    "pooled_patient_level_out_of_fold"
  )
  expect_identical(result$metrics$moderate_calibration$status, "success")
  expect_identical(
    result$metrics$moderate_calibration$patient_calibration$patient_id,
    result$predictions$patient_id
  )
  expect_equal(
    result$metrics$moderate_calibration[c("ici", "e50", "e90")],
    direct_moderate_calibration[c("ici", "e50", "e90")],
    tolerance = 0
  )
  expect_match(
    printed,
    "OOF binary moderate calibration ICI:",
    fixed = TRUE
  )
  expect_length(result$coefficients, 5L)
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("success", 5L)
  )
})

test_that("binary OOF probabilities reproduce analysis-only glm fits", {
  d <- make_sb_binary_fixture(n = 60)
  predictors <- c("gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  result <- sb_run_cv(task, folds = 5, seed = 101)

  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    analysis_rows <- match(split$analysis_ids, d$patient_id)
    assessment_rows <- match(split$assessment_ids, d$patient_id)
    manual <- stats::glm(
      outcome ~ gene_a + gene_b,
      data = d[analysis_rows, , drop = FALSE],
      family = stats::binomial(link = "logit")
    )
    expected <- unname(stats::predict(
      manual,
      newdata = d[assessment_rows, , drop = FALSE],
      type = "response"
    ))
    observed <- result$predictions$probability[
      match(split$assessment_ids, result$predictions$patient_id)
    ]
    expect_equal(observed, expected, tolerance = 1e-12)

    audit <- result$coefficients[[fold]]$values
    design <- stats::model.matrix(
      ~ gene_a + gene_b,
      data = d[assessment_rows, , drop = FALSE]
    )
    reconstructed <- unname(drop(stats::plogis(
      design %*% audit$coefficient[match(colnames(design), audit$term)]
    )))
    expect_equal(observed, reconstructed, tolerance = 1e-12)
  }
})

test_that("binary CV is reproducible and folds never overlap", {
  d <- make_sb_binary_fixture(n = 60)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )
  first <- sb_run_cv(task, folds = 5, seed = 77)
  second <- sb_run_cv(task, folds = 5, seed = 77)

  expect_identical(first$predictions, second$predictions)
  expect_identical(first$metrics, second$metrics)
  expect_identical(first$runs, second$runs)
  expect_identical(first$folds, second$folds)
  for (split in first$folds) {
    expect_length(intersect(split$analysis_ids, split$assessment_ids), 0L)
    expect_setequal(
      c(split$analysis_ids, split$assessment_ids),
      d$patient_id
    )
  }
})

test_that("binary fit failures remain visible without fallback", {
  d <- make_sb_binary_fixture(n = 60)
  d$gene_b <- d$gene_a
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )
  result <- sb_run_cv(task, folds = 5, seed = 9)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(all(is.na(result$predictions$probability)), TRUE)
  expect_identical(result$metrics$roc_auc$status, "incomplete")
  expect_identical(is.na(result$metrics$roc_auc$estimate), TRUE)
  expect_identical(
    result$metrics$average_precision$status,
    "incomplete"
  )
  expect_identical(
    is.na(result$metrics$average_precision$estimate),
    TRUE
  )
  expect_identical(result$metrics$brier_score$status, "incomplete")
  expect_identical(is.na(result$metrics$brier_score$estimate), TRUE)
  expect_identical(result$metrics$log_loss$status, "incomplete")
  expect_identical(is.na(result$metrics$log_loss$estimate), TRUE)
  expect_identical(result$metrics$log_loss$reference$status, "success")
  expect_identical(result$metrics$log_loss$skill_score$status, "incomplete")
  expect_identical(result$metrics$calibration$status, "incomplete")
  expect_identical(
    result$metrics$calibration$input_status,
    "prediction_failed"
  )
  expect_identical(
    result$metrics$moderate_calibration$status,
    "incomplete"
  )
  expect_identical(
    result$metrics$moderate_calibration$input_status,
    "prediction_failed"
  )
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(
    printed,
    "Fold-stratified average precision: unavailable (incomplete)",
    fixed = TRUE
  )
  expect_match(
    printed,
    "OOF log loss: unavailable (incomplete)",
    fixed = TRUE
  )
  expect_match(
    printed,
    "Log-loss skill: unavailable (incomplete)",
    fixed = TRUE
  )
  expect_match(
    printed,
    "OOF calibration: unavailable (incomplete)",
    fixed = TRUE
  )
  expect_match(
    printed,
    "OOF binary moderate calibration: unavailable",
    fixed = TRUE
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_fitted", 5L)
  )
})

test_that("binary CV fits one optional moderate curve over pooled OOF rows", {
  d <- make_sb_binary_fixture(n = 60L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  original_metric <- sb_binary_moderate_calibration_metric
  calls <- 0L
  observed_scored <- NULL
  testthat::local_mocked_bindings(
    sb_binary_moderate_calibration_metric = function(scored, scope) {
      calls <<- calls + 1L
      observed_scored <<- scored
      original_metric(scored, scope)
    }
  )

  result <- sb_run_cv(task, folds = 5L, seed = 20260820L)

  expect_identical(calls, 1L)
  expect_identical(observed_scored, result$predictions)
  expect_identical(result$metrics$moderate_calibration$status, "success")
  expect_identical(
    result$metrics$moderate_calibration$per_fold$status,
    rep("success", 5L)
  )
  expect_identical(
    result$metrics$moderate_calibration$n_folds,
    5L
  )
})

test_that("a binary moderate-calibration error does not erase OOF evidence", {
  d <- make_sb_binary_fixture(n = 60L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  control_result <- sb_run_cv(task, folds = 5L, seed = 114L)
  testthat::local_mocked_bindings(
    sb_binary_moderate_calibration_metric = function(...) {
      condition <- simpleError("injected binary moderate calibration failure")
      class(condition) <- c(
        "test_binary_moderate_calibration_error",
        class(condition)
      )
      stop(condition)
    }
  )

  observed <- sb_run_cv(task, folds = 5L, seed = 114L)

  expect_identical(observed$status, control_result$status)
  expect_identical(observed$complete_oof, control_result$complete_oof)
  expect_identical(observed$predictions, control_result$predictions)
  expect_identical(observed$runs, control_result$runs)
  expect_identical(
    observed$metrics[names(observed$metrics) != "moderate_calibration"],
    control_result$metrics[
      names(control_result$metrics) != "moderate_calibration"
    ]
  )
  expect_identical(
    observed$metrics$moderate_calibration$status,
    "metric_failed"
  )
  expect_identical(
    observed$metrics$moderate_calibration$error_class,
    "test_binary_moderate_calibration_error"
  )
})

test_that("binary folds require enough patients in both classes", {
  d <- make_sb_binary_fixture(n = 20)
  d$outcome <- c(rep(1L, 4), rep(0L, 16))
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")

  error <- capture_sb_error(sb_run_cv(task, folds = 5, seed = 1))

  expect_s3_class(error, "stormbreaker_infeasible_folds")
})

test_that("task and pipeline outcome types cannot be mixed", {
  binary_data <- make_sb_binary_fixture(n = 20)
  binary_task <- sb_binary_task(
    binary_data, "patient_id", "outcome", "gene_a"
  )
  survival_data <- make_sb_fixture(n = 20)
  survival_task <- sb_surv_task(
    survival_data, "patient_id", "time", "event", "gene_a", "month"
  )

  binary_error <- capture_sb_error(
    sb_run_cv(binary_task, sb_pipeline(), folds = 2, seed = 1)
  )
  survival_error <- capture_sb_error(
    sb_run_cv(
      survival_task,
      sb_pipeline(learner = "logistic"),
      folds = 2,
      seed = 1
    )
  )

  expect_s3_class(binary_error, "stormbreaker_pipeline_task_mismatch")
  expect_s3_class(survival_error, "stormbreaker_pipeline_task_mismatch")
})

test_that("assessment outcomes cannot change fixed-split probabilities", {
  d <- make_sb_binary_fixture(n = 60)
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline(learner = "logistic")
  make_task <- function(data) {
    sb_binary_task(
      data, "patient_id", "outcome", c("gene_a", "gene_b")
    )
  }
  original <- sb_run_one_binary_fold(
    make_task(d), pipeline, sb_resolve_pipeline_runtime(pipeline), split,
    seed = 1L
  )
  changed <- d
  changed$outcome[49:60] <- 1L - changed$outcome[49:60]
  perturbed <- sb_run_one_binary_fold(
    make_task(changed), pipeline, sb_resolve_pipeline_runtime(pipeline),
    split, seed = 1L
  )

  expect_equal(
    original$predictions$probability,
    perturbed$predictions$probability,
    tolerance = 0
  )

  brier_from_fold <- function(fold_result) {
    reference <- data.frame(
      fold = fold_result$run$fold,
      analysis_prevalence =
        fold_result$run$n_analysis_positive / fold_result$run$n_analysis
    )
    sb_binary_brier_metric(fold_result$predictions, reference)
  }
  original_brier <- brier_from_fold(original)
  perturbed_brier <- brier_from_fold(perturbed)
  expect_identical(
    original_brier$reference$per_fold$analysis_prevalence,
    perturbed_brier$reference$per_fold$analysis_prevalence
  )
  reference_probability <-
    original_brier$reference$per_fold$analysis_prevalence
  expect_equal(
    original_brier$reference$estimate,
    mean((original$predictions$outcome - reference_probability)^2),
    tolerance = 0
  )
  expect_equal(
    perturbed_brier$reference$estimate,
    mean((perturbed$predictions$outcome - reference_probability)^2),
    tolerance = 0
  )

  log_loss_from_fold <- function(fold_result) {
    reference <- data.frame(
      fold = fold_result$run$fold,
      analysis_prevalence =
        fold_result$run$n_analysis_positive / fold_result$run$n_analysis
    )
    sb_binary_log_loss_metric(fold_result$predictions, reference)
  }
  original_log_loss <- log_loss_from_fold(original)
  perturbed_log_loss <- log_loss_from_fold(perturbed)
  expect_identical(
    original_log_loss$reference$per_fold$analysis_prevalence,
    perturbed_log_loss$reference$per_fold$analysis_prevalence
  )
  original_reference_loss <- ifelse(
    original$predictions$outcome == 1L,
    -log(reference_probability),
    -log1p(-reference_probability)
  )
  perturbed_reference_loss <- ifelse(
    perturbed$predictions$outcome == 1L,
    -log(reference_probability),
    -log1p(-reference_probability)
  )
  expect_equal(
    original_log_loss$reference$estimate,
    mean(original_reference_loss),
    tolerance = 0
  )
  expect_equal(
    perturbed_log_loss$reference$estimate,
    mean(perturbed_reference_loss),
    tolerance = 0
  )
})

test_that("binary probability references use each outer analysis prevalence", {
  d <- make_sb_binary_fixture(n = 60L, seed = 2101L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )
  result <- sb_run_cv(task, folds = 5L, seed = 2201L)
  reference <- result$metrics$brier_score$reference
  skill <- result$metrics$brier_score$skill_score
  analysis_prevalence <- vapply(result$folds, function(split) {
    mean(d$outcome[match(split$analysis_ids, d$patient_id)])
  }, numeric(1))
  patient_reference <- analysis_prevalence[result$predictions$fold]
  reference_loss <- (
    result$predictions$outcome - patient_reference
  )^2
  model_loss <- (
    result$predictions$outcome - result$predictions$probability
  )^2

  expect_identical(reference$status, "success")
  expect_equal(
    reference$per_fold$analysis_prevalence,
    analysis_prevalence,
    tolerance = 0
  )
  expect_equal(reference$estimate, mean(reference_loss), tolerance = 0)
  expect_identical(skill$status, "success")
  expect_equal(
    skill$estimate,
    1 - sum(model_loss) / sum(reference_loss),
    tolerance = 0
  )

  log_loss <- result$metrics$log_loss
  reference_log_loss <- ifelse(
    result$predictions$outcome == 1L,
    -log(patient_reference),
    -log1p(-patient_reference)
  )
  model_log_loss <- ifelse(
    result$predictions$outcome == 1L,
    -log(result$predictions$probability),
    -log1p(-result$predictions$probability)
  )
  expect_identical(log_loss$reference$status, "success")
  expect_equal(
    log_loss$reference$per_fold$analysis_prevalence,
    analysis_prevalence,
    tolerance = 0
  )
  expect_equal(
    log_loss$reference$estimate,
    mean(reference_log_loss),
    tolerance = 0
  )
  expect_identical(log_loss$skill_score$status, "success")
  expect_equal(
    log_loss$skill_score$estimate,
    1 - sum(model_log_loss) / sum(reference_log_loss),
    tolerance = 1e-14
  )
})

test_that("the binary result schema is explicit", {
  d <- make_sb_binary_fixture(n = 20)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  result <- sb_run_cv(task, folds = 2, seed = 2)

  expect_named(
    result,
    c(
      "status", "complete_oof", "predictions", "metrics", "runs",
      "selection_stability", "coefficients", "folds", "task", "pipeline",
      "provenance", "glmnet_active_set_stability"
    ),
    ignore.order = FALSE
  )
  expect_identical(result$selection_stability$status, "not_applicable")
  expect_identical(result$selection_stability$n_successful_folds, 0L)
  expect_true(is.na(result$selection_stability$mean_pairwise_jaccard))
  expect_identical(
    result$glmnet_active_set_stability$status,
    "not_applicable"
  )
  expect_true(is.na(result$glmnet_active_set_stability$alpha))
  expect_named(
    result$predictions,
    c(
      "patient_id", "fold", "outcome", "positive_class", "probability",
      "status", "pipeline_id", "learner"
    ),
    ignore.order = FALSE
  )
})

test_that("mutated binary task metadata is rejected", {
  d <- make_sb_binary_fixture(n = 20)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  task$n_positive <- task$n_positive + 1L

  error <- capture_sb_error(sb_run_cv(task, folds = 2, seed = 1))

  expect_s3_class(error, "stormbreaker_mutated_task")
})

test_that("nested glmnet logistic supports predictors exceeding patients", {
  skip_if_not_installed("glmnet")

  d <- make_sb_high_dim_binary_fixture(n = 60, p = 80)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 1701L
  )
  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 101L)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(nrow(result$predictions), 60L)
  expect_identical(
    all(is.finite(result$predictions$probability)),
    TRUE
  )
  expect_length(result$tuning, 5L)
  expect_identical(
    vapply(result$tuning, `[[`, "", "status"),
    rep("success", 5L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("success", 5L)
  )
  expect_identical(
    result$provenance$glmnet_version,
    as.character(utils::packageVersion("glmnet"))
  )

  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    audit <- result$tuning[[fold]]
    expect_setequal(names(audit$fold_id), split$analysis_ids)
    expect_length(intersect(names(audit$fold_id), split$assessment_ids), 0L)
    expect_identical(all(is.finite(audit$curve$mean_deviance)), TRUE)
  }
})

test_that("nested binary OOF probabilities match direct glmnet execution", {
  skip_if_not_installed("glmnet")

  d <- make_sb_high_dim_binary_fixture(n = 60, p = 20)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0.4,
    inner_folds = 3L,
    inner_seed = 301L,
    glmnet_standardize = FALSE,
    glmnet_lambda_rule = "lambda.1se",
    glmnet_tuning_measure = "auc"
  )
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 401L)
  split <- result$folds[[1L]]
  analysis_rows <- match(split$analysis_ids, d$patient_id)
  assessment_rows <- match(split$assessment_ids, d$patient_id)
  tuning <- result$tuning[[1L]]
  manual_cv <- glmnet::cv.glmnet(
    x = as.matrix(d[analysis_rows, predictors, drop = FALSE]),
    y = d$outcome[analysis_rows],
    family = "binomial",
    alpha = 0.4,
    foldid = unname(tuning$fold_id[split$analysis_ids]),
    type.measure = "auc",
    standardize = FALSE
  )
  manual_fit <- glmnet::glmnet(
    x = as.matrix(d[analysis_rows, predictors, drop = FALSE]),
    y = d$outcome[analysis_rows],
    family = "binomial",
    alpha = 0.4,
    standardize = FALSE
  )
  expected <- unname(drop(stats::predict(
    manual_fit,
    newx = as.matrix(d[assessment_rows, predictors, drop = FALSE]),
    s = manual_cv$lambda.1se,
    type = "response"
  )))
  observed <- result$predictions$probability[
    match(split$assessment_ids, result$predictions$patient_id)
  ]

  expect_equal(tuning$lambda, unname(manual_cv$lambda.1se), tolerance = 0)
  expect_identical(tuning$selection_rule, "lambda.1se")
  expect_identical(tuning$standardize, FALSE)
  expect_identical(tuning$measure, "auc")
  expect_equal(observed, expected, tolerance = 1e-12)

  coefficient_audit <- result$coefficients[[1L]]$values
  design <- cbind(
    "(Intercept)" = 1,
    as.matrix(d[assessment_rows, predictors, drop = FALSE])
  )
  reconstructed <- unname(drop(stats::plogis(
    design %*% coefficient_audit$coefficient[
      match(colnames(design), coefficient_audit$term)
    ]
  )))
  expect_equal(observed, reconstructed, tolerance = 1e-12)
})

test_that("outer assessment changes cannot alter nested binary fitting", {
  skip_if_not_installed("glmnet")

  d <- make_sb_high_dim_binary_fixture(n = 60, p = 20)
  predictors <- grep("^gene_", names(d), value = TRUE)
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 1,
    inner_folds = 3L,
    inner_seed = 501L
  )
  make_task <- function(data) {
    sb_binary_task(data, "patient_id", "outcome", predictors)
  }
  original <- sb_run_one_binary_fold(
    make_task(d), pipeline, sb_resolve_pipeline_runtime(pipeline), split,
    seed = 1L
  )
  changed <- d
  changed$outcome[49:60] <- 1L - changed$outcome[49:60]
  perturbed <- sb_run_one_binary_fold(
    make_task(changed), pipeline, sb_resolve_pipeline_runtime(pipeline),
    split, seed = 1L
  )

  expect_identical(original$tuning, perturbed$tuning)
  expect_equal(
    original$predictions$probability,
    perturbed$predictions$probability,
    tolerance = 0
  )
})

test_that("nested binary tuning failures remain visible without fallback", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 30)
  d$constant_a <- 1
  d$constant_b <- 2
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("constant_a", "constant_b")
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 1L
  )
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 1L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_glmnet_logistic_tuning_failure"
  )
  expect_identical(
    vapply(result$tuning, `[[`, "", "status"),
    rep("failed", 3L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_fitted", 3L)
  )
  expect_identical(all(is.na(result$predictions$probability)), TRUE)
  expect_identical(result$metrics$roc_auc$status, "incomplete")
  expect_identical(result$metrics$brier_score$status, "incomplete")
})

test_that("binary AUC Top-K selects and audits features inside every outer fold", {
  d <- make_sb_high_dim_binary_fixture(n = 80L, p = 30L, seed = 1201L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 3L
  )
  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 1301L)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_length(result$selections, 4L)
  expect_length(result$coefficients, 4L)
  expect_identical(
    vapply(result$selections, `[[`, "", "status"),
    rep("success", 4L)
  )
  stability <- result$selection_stability
  expected_count <- vapply(predictors, function(feature) {
    as.integer(sum(vapply(
      result$selections,
      function(selection) feature %in% selection$selected_features,
      logical(1)
    )))
  }, integer(1))
  expected_jaccard <- vapply(
    seq_len(nrow(stability$pairwise)),
    function(row) {
      selected_i <- result$selections[[
        stability$pairwise$fold_i[[row]]
      ]]$selected_features
      selected_j <- result$selections[[
        stability$pairwise$fold_j[[row]]
      ]]$selected_features
      length(intersect(selected_i, selected_j)) /
        length(union(selected_i, selected_j))
    },
    numeric(1)
  )
  expect_identical(stability$status, "success")
  expect_identical(stability$n_successful_folds, 4L)
  expect_identical(stability$feature_frequency$feature, predictors)
  expect_identical(
    stability$feature_frequency$selected_count,
    unname(expected_count)
  )
  expect_equal(stability$pairwise$jaccard, expected_jaccard, tolerance = 0)
  expect_equal(
    stability$mean_pairwise_jaccard,
    mean(expected_jaccard),
    tolerance = 0
  )

  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    audit <- result$selections[[fold]]
    analysis_rows <- match(split$analysis_ids, d$patient_id)
    assessment_rows <- match(split$assessment_ids, d$patient_id)
    direct <- sb_fit_univ_binary_auc_top_k(
      d[analysis_rows, predictors, drop = FALSE],
      d$outcome[analysis_rows],
      n_features = 3L
    )

    expect_identical(audit$fold, fold)
    expect_identical(audit$selector, "univ_binary_auc_top_k")
    expect_identical(audit$selected_features, direct$selected_features)
    expect_identical(audit$scores, direct$scores)
    expect_length(audit$selected_features, 3L)
    expect_identical(sum(audit$scores$selected), 3L)
    expect_identical(nrow(audit$scores), length(predictors))

    coefficients <- result$coefficients[[fold]]$values
    expect_identical(
      coefficients$term,
      c("(Intercept)", audit$selected_features)
    )
    design <- cbind(
      "(Intercept)" = 1,
      as.matrix(d[assessment_rows, audit$selected_features, drop = FALSE])
    )
    rebuilt <- unname(drop(stats::plogis(
      design %*% coefficients$coefficient[
        match(colnames(design), coefficients$term)
      ]
    )))
    observed <- result$predictions$probability[
      match(split$assessment_ids, result$predictions$patient_id)
    ]
    expect_equal(observed, rebuilt, tolerance = 1e-12)
  }
})

test_that("outer assessment data cannot alter binary AUC selection or fitting", {
  d <- make_sb_high_dim_binary_fixture(n = 60L, p = 12L, seed = 1401L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 3L
  )
  run_fixed <- function(data) {
    task <- sb_binary_task(data, "patient_id", "outcome", predictors)
    sb_run_one_binary_fold(
      task,
      pipeline,
      sb_resolve_pipeline_runtime(pipeline),
      split,
      seed = 1L
    )
  }
  original <- run_fixed(d)

  changed_outcome <- d
  changed_outcome$outcome[split$assessment_rows] <-
    1L - changed_outcome$outcome[split$assessment_rows]
  outcome_run <- run_fixed(changed_outcome)

  expect_identical(original$selection, outcome_run$selection)
  expect_identical(original$coefficients, outcome_run$coefficients)
  expect_equal(
    original$predictions$probability,
    outcome_run$predictions$probability,
    tolerance = 0
  )

  changed_predictors <- d
  selected <- original$selection$selected_features[[1L]]
  changed_predictors[split$assessment_rows, selected] <-
    changed_predictors[split$assessment_rows, selected] + 100
  predictor_run <- run_fixed(changed_predictors)

  expect_identical(original$selection, predictor_run$selection)
  expect_identical(original$coefficients, predictor_run$coefficients)
  expect_false(isTRUE(all.equal(
    original$predictions$probability,
    predictor_run$predictions$probability,
    tolerance = 0
  )))
})

test_that("learned binary selector results expose selection audit explicitly", {
  d <- make_sb_binary_fixture(n = 40L, seed = 1501L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 1L
  )
  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 1601L)

  expect_named(
    result,
    c(
      "status", "complete_oof", "predictions", "metrics", "runs",
      "selections", "selection_stability", "coefficients", "folds", "task",
      "pipeline", "provenance", "glmnet_active_set_stability"
    ),
    ignore.order = FALSE
  )
})

test_that("binary selector failures remain visible without fallback", {
  d <- make_sb_binary_fixture(n = 40L, seed = 1901L)
  d$constant_a <- 1
  d$constant_b <- 2
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("constant_a", "constant_b")
  )
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 1L
  )
  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 2001L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_insufficient_valid_features"
  )
  expect_identical(
    vapply(result$selections, `[[`, "", "status"),
    rep("failed", 4L)
  )
  expect_identical(
    vapply(result$selections, `[[`, "", "error_class"),
    rep("stormbreaker_insufficient_valid_features", 4L)
  )
  expect_identical(result$selection_stability$status, "incomplete")
  expect_identical(result$selection_stability$n_successful_folds, 0L)
  expect_identical(result$selection_stability$n_evaluable_pairs, 0L)
  expect_true(is.na(result$selection_stability$mean_pairwise_jaccard))
  expect_identical(
    result$selection_stability$per_fold$error_class,
    rep("stormbreaker_insufficient_valid_features", 4L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_fitted", 4L)
  )
  expect_identical(all(is.na(result$predictions$probability)), TRUE)
  expect_identical(result$metrics$roc_auc$status, "incomplete")
  expect_identical(result$metrics$brier_score$status, "incomplete")
})

test_that("Gaussian Naive Bayes supports high-dimensional numeric candidates", {
  d <- make_sb_high_dim_binary_fixture(n = 60L, p = 80L, seed = 2301L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(
    learner = "gaussian_naive_bayes",
    naive_bayes_var_smoothing = 1e-9
  )
  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 2401L)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(all(is.finite(result$predictions$probability)), TRUE)
  expect_length(result$distribution_parameters, 5L)
  expect_identical(
    vapply(result$distribution_parameters, `[[`, "", "status"),
    rep("success", 5L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_applicable", 5L)
  )
  expect_identical(result$metrics$roc_auc$status, "success")
  expect_identical(result$metrics$brier_score$status, "success")
  expect_identical(
    result$metrics$brier_score$reference$status,
    "success"
  )

  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    audit <- result$distribution_parameters[[fold]]
    assessment_rows <- match(split$assessment_ids, d$patient_id)
    means <- do.call(rbind, lapply(c(0L, 1L), function(class) {
      values <- audit$values[audit$values$class == class, , drop = FALSE]
      stats::setNames(values$mean, values$feature)
    }))
    rownames(means) <- c("0", "1")
    variances <- do.call(rbind, lapply(c(0L, 1L), function(class) {
      values <- audit$values[audit$values$class == class, , drop = FALSE]
      stats::setNames(values$variance, values$feature)
    }))
    rownames(variances) <- c("0", "1")
    priors <- stats::setNames(audit$priors$prior, audit$priors$class)
    rebuilt <- manual_gaussian_nb_probability(
      as.matrix(d[assessment_rows, predictors, drop = FALSE]),
      priors,
      means,
      variances
    )
    observed <- result$predictions$probability[
      match(split$assessment_ids, result$predictions$patient_id)
    ]
    expect_equal(observed, rebuilt, tolerance = 1e-14)
    expect_identical(nrow(audit$values), 2L * length(predictors))
  }
})

test_that("assessment data cannot alter Gaussian Naive Bayes fitting", {
  d <- make_sb_high_dim_binary_fixture(n = 60L, p = 12L, seed = 2501L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline(learner = "gaussian_naive_bayes")
  run_fixed <- function(data) {
    task <- sb_binary_task(data, "patient_id", "outcome", predictors)
    sb_run_one_binary_fold(
      task,
      pipeline,
      sb_resolve_pipeline_runtime(pipeline),
      split,
      seed = 1L
    )
  }
  original <- run_fixed(d)

  changed_outcome <- d
  changed_outcome$outcome[split$assessment_rows] <-
    1L - changed_outcome$outcome[split$assessment_rows]
  outcome_run <- run_fixed(changed_outcome)
  expect_identical(
    original$distribution_parameters,
    outcome_run$distribution_parameters
  )
  expect_equal(
    original$predictions$probability,
    outcome_run$predictions$probability,
    tolerance = 0
  )

  changed_predictors <- d
  changed_predictors[split$assessment_rows, predictors[[1L]]] <-
    changed_predictors[split$assessment_rows, predictors[[1L]]] + 100
  predictor_run <- run_fixed(changed_predictors)
  expect_identical(
    original$distribution_parameters,
    predictor_run$distribution_parameters
  )
  expect_false(isTRUE(all.equal(
    original$predictions$probability,
    predictor_run$predictions$probability,
    tolerance = 0
  )))
})

test_that("Gaussian Naive Bayes variance failures have no fallback", {
  d <- make_sb_binary_fixture(n = 40L, seed = 2601L)
  d$constant_a <- 1
  d$constant_b <- 2
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("constant_a", "constant_b")
  )
  pipeline <- sb_pipeline(learner = "gaussian_naive_bayes")
  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 2701L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_invalid_naive_bayes_variance"
  )
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(
    vapply(result$distribution_parameters, `[[`, "", "status"),
    rep("not_fitted", 4L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_applicable", 4L)
  )
  expect_identical(all(is.na(result$predictions$probability)), TRUE)
  expect_identical(result$metrics$brier_score$status, "incomplete")
  expect_identical(
    result$metrics$brier_score$reference$status,
    "success"
  )
  expect_identical(
    result$metrics$brier_score$skill_score$status,
    "incomplete"
  )
})

test_that("fixed LightGBM supports high-dimensional numeric candidates", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  d <- make_sb_high_dim_binary_fixture(n = 48L, p = 64L, seed = 3201L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(
    learner = "lightgbm_binary",
    lightgbm_nrounds = 8L,
    lightgbm_learning_rate = 0.1,
    lightgbm_num_leaves = 7L,
    lightgbm_min_data_in_leaf = 3L,
    lightgbm_seed = 3301L
  )
  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 3401L)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(all(is.finite(result$predictions$probability)), TRUE)
  expect_length(result$tree_models, 4L)
  expect_identical(
    vapply(result$tree_models, `[[`, "", "status"),
    rep("success", 4L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_applicable", 4L)
  )
  expect_identical(result$metrics$roc_auc$status, "success")
  expect_identical(result$metrics$brier_score$status, "success")
  expect_identical(
    result$metrics$brier_score$reference$status,
    "success"
  )

  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    audit <- result$tree_models[[fold]]
    assessment_rows <- match(split$assessment_ids, d$patient_id)
    restored <- lightgbm::lgb.load(model_str = audit$model_string)
    rebuilt <- as.numeric(stats::predict(
      restored,
      as.matrix(d[assessment_rows, audit$features, drop = FALSE]),
      type = "response",
      num_iteration = audit$nrounds,
      params = list(num_threads = 1L)
    ))
    observed <- result$predictions$probability[
      match(split$assessment_ids, result$predictions$patient_id)
    ]
    expect_equal(observed, rebuilt, tolerance = 0)
    expect_identical(audit$features, predictors)
    expect_identical(audit$feature_importance$feature, predictors)
  }
})

test_that("assessment data cannot alter fixed LightGBM fitting", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  d <- make_sb_high_dim_binary_fixture(n = 60L, p = 12L, seed = 3501L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline(
    learner = "lightgbm_binary",
    lightgbm_nrounds = 8L,
    lightgbm_num_leaves = 7L,
    lightgbm_min_data_in_leaf = 3L,
    lightgbm_seed = 3601L
  )
  run_fixed <- function(data) {
    task <- sb_binary_task(data, "patient_id", "outcome", predictors)
    sb_run_one_binary_fold(
      task,
      pipeline,
      sb_resolve_pipeline_runtime(pipeline),
      split,
      seed = 1L
    )
  }
  original <- run_fixed(d)

  changed_outcome <- d
  changed_outcome$outcome[split$assessment_rows] <-
    1L - changed_outcome$outcome[split$assessment_rows]
  outcome_run <- run_fixed(changed_outcome)
  expect_identical(original$tree_model, outcome_run$tree_model)
  expect_equal(
    original$predictions$probability,
    outcome_run$predictions$probability,
    tolerance = 0
  )

  changed_predictors <- d
  changed_predictors[split$assessment_rows, predictors[[1L]]] <-
    changed_predictors[split$assessment_rows, predictors[[1L]]] + 100
  predictor_run <- run_fixed(changed_predictors)
  expect_identical(original$tree_model, predictor_run$tree_model)
  expect_false(isTRUE(all.equal(
    original$predictions$probability,
    predictor_run$predictions$probability,
    tolerance = 0
  )))
})

test_that("missing fixed LightGBM dependency has no fallback", {
  d <- make_sb_binary_fixture(n = 40L, seed = 3701L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )
  pipeline <- sb_pipeline(
    learner = "lightgbm_binary",
    lightgbm_nrounds = 5L,
    lightgbm_min_data_in_leaf = 2L
  )
  local_mocked_bindings(sb_lightgbm_available = function() FALSE)

  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 3801L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_missing_lightgbm"
  )
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(
    vapply(result$tree_models, `[[`, "", "status"),
    rep("not_fitted", 4L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_applicable", 4L)
  )
  expect_identical(all(is.na(result$predictions$probability)), TRUE)
  expect_identical(result$metrics$roc_auc$status, "incomplete")
  expect_identical(result$metrics$brier_score$status, "incomplete")
})

test_that("fixed CatBoost supports high-dimensional numeric candidates", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_high_dim_binary_fixture(n = 48L, p = 64L, seed = 4601L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(
    learner = "catboost_binary",
    catboost_iterations = 8L,
    catboost_learning_rate = 0.1,
    catboost_depth = 3L,
    catboost_l2_leaf_reg = 3,
    catboost_seed = 4701L
  )
  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 4801L)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(all(is.finite(result$predictions$probability)), TRUE)
  expect_length(result$tree_models, 4L)
  expect_identical(
    vapply(result$tree_models, `[[`, "", "status"),
    rep("success", 4L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_applicable", 4L)
  )
  expect_identical(result$metrics$roc_auc$status, "success")
  expect_identical(result$metrics$brier_score$status, "success")
  expect_identical(
    result$metrics$brier_score$reference$status,
    "success"
  )
  expect_identical(
    result$provenance$catboost_version,
    as.character(utils::packageVersion("catboost"))
  )

  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    audit <- result$tree_models[[fold]]
    assessment_rows <- match(split$assessment_ids, d$patient_id)
    restored <- unserialize(audit$model_raw)
    restored_pool <- catboost::catboost.load_pool(as.matrix(
      d[assessment_rows, audit$features, drop = FALSE]
    ))
    rebuilt <- as.numeric(catboost::catboost.predict(
      restored,
      restored_pool,
      prediction_type = "Probability",
      ntree_start = 0L,
      ntree_end = audit$iterations,
      thread_count = 1L
    ))
    observed <- result$predictions$probability[
      match(split$assessment_ids, result$predictions$patient_id)
    ]
    expect_equal(observed, rebuilt, tolerance = 0)
    expect_identical(audit$tree_count, 8L)
    expect_identical(audit$features, predictors)
    expect_identical(audit$feature_importance$feature, predictors)
  }
})

test_that("assessment data cannot alter fixed CatBoost fitting", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_high_dim_binary_fixture(n = 60L, p = 12L, seed = 4901L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline(
    learner = "catboost_binary",
    catboost_iterations = 8L,
    catboost_depth = 3L,
    catboost_seed = 5001L
  )
  run_fixed <- function(data) {
    task <- sb_binary_task(data, "patient_id", "outcome", predictors)
    sb_run_one_binary_fold(
      task,
      pipeline,
      sb_resolve_pipeline_runtime(pipeline),
      split,
      seed = 1L
    )
  }
  predict_saved <- function(audit, data) {
    model <- unserialize(audit$model_raw)
    pool <- catboost::catboost.load_pool(as.matrix(
      data[split$assessment_rows, audit$features, drop = FALSE]
    ))
    as.numeric(catboost::catboost.predict(
      model,
      pool,
      prediction_type = "Probability",
      ntree_start = 0L,
      ntree_end = audit$iterations,
      thread_count = 1L
    ))
  }
  original <- run_fixed(d)

  changed_outcome <- d
  changed_outcome$outcome[split$assessment_rows] <-
    1L - changed_outcome$outcome[split$assessment_rows]
  outcome_run <- run_fixed(changed_outcome)
  expect_equal(
    original$predictions$probability,
    outcome_run$predictions$probability,
    tolerance = 0
  )
  expect_equal(
    predict_saved(original$tree_model, d),
    predict_saved(outcome_run$tree_model, d),
    tolerance = 0
  )
  expect_identical(
    original$tree_model$feature_importance,
    outcome_run$tree_model$feature_importance
  )

  changed_predictors <- d
  changed_predictors[split$assessment_rows, predictors] <-
    -changed_predictors[split$assessment_rows, predictors] + 100
  predictor_run <- run_fixed(changed_predictors)
  expect_equal(
    predict_saved(original$tree_model, d),
    predict_saved(predictor_run$tree_model, d),
    tolerance = 0
  )
  expect_identical(
    original$tree_model$feature_importance,
    predictor_run$tree_model$feature_importance
  )
  expect_false(isTRUE(all.equal(
    original$predictions$probability,
    predictor_run$predictions$probability,
    tolerance = 0
  )))
})

test_that("missing fixed CatBoost dependency has no fallback", {
  d <- make_sb_binary_fixture(n = 40L, seed = 5101L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c"),
    decision_thresholds = c(0.2, 0.5)
  )
  pipeline <- sb_pipeline(
    learner = "catboost_binary",
    catboost_iterations = 5L,
    catboost_depth = 2L
  )
  local_mocked_bindings(sb_catboost_available = function() FALSE)

  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 5201L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_missing_catboost"
  )
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(
    vapply(result$tree_models, `[[`, "", "status"),
    rep("not_fitted", 4L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_applicable", 4L)
  )
  expect_identical(all(is.na(result$predictions$probability)), TRUE)
  expect_identical(result$metrics$roc_auc$status, "incomplete")
  expect_identical(result$metrics$brier_score$status, "incomplete")
  expect_identical(result$metrics$decision_curve$status, "incomplete")
  expect_identical(result$metrics$decision_curve$n_scored, 0L)
})

test_that("binary CV decision curve is reconstructed from OOF patients", {
  d <- make_sb_binary_fixture(n = 60L, seed = 5601L)
  thresholds <- c(0.25, 0.50, 0.75)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c"),
    decision_thresholds = thresholds
  )

  result <- sb_run_cv(
    task,
    sb_pipeline(learner = "logistic"),
    folds = 3L,
    seed = 5701L
  )
  metric <- result$metrics$decision_curve

  expect_identical(result$status, "success")
  expect_identical(metric$status, "success")
  expect_identical(metric$scope, "pooled_patient_level_out_of_fold")
  expect_identical(metric$n_scored, task$n)
  expect_identical(result$task$decision_thresholds, thresholds)
  expect_identical(metric$curve$threshold_probability, thresholds)

  for (index in seq_along(thresholds)) {
    threshold <- thresholds[[index]]
    selected <- result$predictions$probability >= threshold
    true_positive <- as.integer(sum(
      selected & result$predictions$outcome == 1L
    ))
    false_positive <- as.integer(sum(
      selected & result$predictions$outcome == 0L
    ))
    true_negative <- as.integer(sum(
      !selected & result$predictions$outcome == 0L
    ))
    false_negative <- as.integer(sum(
      !selected & result$predictions$outcome == 1L
    ))
    expected_net_benefit <- true_positive / task$n -
      false_positive / task$n * threshold / (1 - threshold)

    expect_identical(metric$curve$true_positive[[index]], true_positive)
    expect_identical(metric$curve$false_positive[[index]], false_positive)
    expect_identical(metric$curve$true_negative[[index]], true_negative)
    expect_identical(metric$curve$false_negative[[index]], false_negative)
    expect_equal(
      metric$curve$model_net_benefit[[index]],
      expected_net_benefit,
      tolerance = 1e-15
    )
  }
})

test_that("binary CV omits decision curve without a threshold declaration", {
  d <- make_sb_binary_fixture(n = 60L, seed = 5801L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )

  result <- sb_run_cv(
    task,
    sb_pipeline(learner = "logistic"),
    folds = 3L,
    seed = 5901L
  )

  expect_false("decision_curve" %in% names(result$metrics))
  expect_false("decision_thresholds" %in% names(result$task))
})

test_that("binary CV reuses an exact patient-keyed outer assignment", {
  d <- make_sb_binary_fixture(n = 60L, seed = 6201L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 2L
  )
  generated <- sb_run_cv(task, pipeline, folds = 5L, seed = 6202L)
  assignment <- setNames(
    generated$predictions$fold,
    generated$predictions$patient_id
  )
  supplied <- assignment[c(
    seq(2L, length(assignment), by = 2L),
    seq(1L, length(assignment), by = 2L)
  )]

  reused <- sb_run_cv(task, pipeline, folds = supplied, seed = 999L)
  task_ids <- as.character(task$data[[task$id_col]])
  expected_assignment <- data.frame(
    patient_id = task_ids,
    fold = as.integer(supplied[match(task_ids, names(supplied))]),
    stringsAsFactors = FALSE
  )

  expect_identical(reused$folds, generated$folds)
  expect_identical(reused$predictions, generated$predictions)
  expect_identical(reused$metrics, generated$metrics)
  expect_identical(reused$coefficients, generated$coefficients)
  expect_identical(reused$selections, generated$selections)
  expect_identical(
    reused$selection_stability,
    generated$selection_stability
  )
  expect_identical(
    reused$runs[setdiff(names(reused$runs), "seed")],
    generated$runs[setdiff(names(generated$runs), "seed")]
  )
  expect_identical(generated$runs$seed, rep(6202L, 5L))
  expect_identical(reused$runs$seed, rep(NA_integer_, 5L))
  expect_identical(generated$provenance$fold_source, "generated_stratified")
  expect_identical(reused$provenance$fold_source, "preassigned_patient_id")
  expect_identical(generated$provenance$seed, 6202L)
  expect_identical(reused$provenance$seed, NA_integer_)
  expect_identical(generated$provenance$folds, 5L)
  expect_identical(reused$provenance$folds, 5L)
  expect_identical(reused$provenance$fold_assignment, expected_assignment)
  expect_identical(
    generated$provenance$fold_assignment,
    expected_assignment
  )
  printed <- paste(capture.output(print(reused)), collapse = "\n")
  expect_match(printed, "preassigned patient-ID folds", fixed = TRUE)
})

test_that("preassigned binary folds are not silently repaired", {
  d <- make_sb_binary_fixture(n = 20L, seed = 6211L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  supplied <- setNames(
    ifelse(d$outcome == 1L, 1L, 2L),
    d$patient_id
  )

  result <- sb_run_cv(
    task,
    sb_pipeline(learner = "logistic"),
    folds = supplied,
    seed = 6212L
  )

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(result$provenance$fold_source, "preassigned_patient_id")
  expect_identical(result$provenance$seed, NA_integer_)
  expect_identical(
    result$provenance$fold_assignment$fold,
    as.integer(supplied)
  )
  expect_identical(
    result$runs$error_class,
    rep("stormbreaker_logistic_single_class", 2L)
  )
  expect_identical(result$runs$status, rep("failed", 2L))
})
})

# ---- cv-comparison ----
local({
# ---- compare-cv ----

test_that("comparison uses identical outer folds and copied OOF metrics", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260816)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    ),
    sb_pipeline(
      learner = "glmnet_cox",
      inner_folds = 3L,
      inner_seed = 1701L
    )
  )

  comparison <- sb_compare_cv(
    task,
    pipelines,
    folds = 5L,
    seed = 101L
  )

  expect_s3_class(comparison, "stormbreaker_cv_comparison")
  expect_identical(comparison$status, "success")
  expect_identical(comparison$n_patients, 60L)
  expect_identical(comparison$folds, 5L)
  expect_identical(comparison$seed, 101L)
  required_columns <- c(
    "pipeline_id", "selector", "learner", "status", "complete_oof",
    "harrell_c_status", "harrell_c", "comparable_pairs",
    "selection_stability_status", "mean_pairwise_jaccard",
    "n_successful_selection_folds",
    "n_folds", "n_successful_folds",
    "glmnet_active_set_stability_status",
    "glmnet_active_set_mean_pairwise_jaccard",
    "n_successful_glmnet_active_set_folds"
  )
  expect_length(
    setdiff(required_columns, names(comparison$comparison)),
    0L
  )
  expected_ids <- vapply(pipelines, function(x) x$id, character(1))
  expect_identical(comparison$comparison$pipeline_id, expected_ids)
  expect_identical(names(comparison$results), expected_ids)
  expect_identical(comparison$comparison$complete_oof, rep(TRUE, 3L))
  expect_identical(
    comparison$comparison$harrell_c_status,
    rep("success", 3L)
  )
  expect_identical(
    comparison$comparison$selection_stability_status,
    c("not_applicable", "success", "not_applicable")
  )
  expect_identical(
    is.finite(comparison$comparison$mean_pairwise_jaccard),
    c(FALSE, TRUE, FALSE)
  )
  expect_identical(
    comparison$comparison$n_successful_selection_folds,
    c(0L, 5L, 0L)
  )
  expect_equal(
    comparison$comparison$mean_pairwise_jaccard,
    unname(vapply(
      comparison$results,
      function(result) result$selection_stability$mean_pairwise_jaccard,
      numeric(1)
    )),
    tolerance = 0
  )
  expected_active_status <- unname(vapply(
    comparison$results,
    function(result) result$glmnet_active_set_stability$status,
    character(1)
  ))
  expect_identical(
    comparison$comparison$glmnet_active_set_stability_status,
    expected_active_status
  )
  expect_identical(
    expected_active_status[1:2],
    rep("not_applicable", 2L)
  )
  expect_match(expected_active_status[[3L]], "^(success|not_evaluable)$")
  expect_equal(
    comparison$comparison$glmnet_active_set_mean_pairwise_jaccard,
    unname(vapply(
      comparison$results,
      function(result) {
        result$glmnet_active_set_stability$mean_pairwise_jaccard
      },
      numeric(1)
    )),
    tolerance = 0
  )
  expect_identical(
    comparison$comparison$n_successful_glmnet_active_set_folds,
    c(0L, 0L, 5L)
  )
  expect_identical(
    vapply(comparison$results[[1L]]$coefficients, `[[`, "", "status"),
    rep("success", 5L)
  )
  expect_identical(
    vapply(comparison$results[[2L]]$coefficients, `[[`, "", "status"),
    rep("success", 5L)
  )

  assignments <- lapply(
    comparison$results,
    function(x) x$predictions[c("patient_id", "fold")]
  )
  for (assignment in assignments[-1L]) {
    expect_identical(assignment, assignments[[1L]])
  }
  expected_c <- vapply(
    comparison$results,
    function(x) x$metrics$harrell_c$estimate,
    numeric(1)
  )
  expected_pairs <- vapply(
    comparison$results,
    function(x) x$metrics$harrell_c$comparable_pairs,
    numeric(1)
  )
  expect_equal(
    comparison$comparison$harrell_c,
    unname(expected_c),
    tolerance = 0
  )
  expect_equal(
    comparison$comparison$comparable_pairs,
    unname(expected_pairs),
    tolerance = 0
  )

  printed <- paste(capture.output(print(comparison)), collapse = "\n")
  expect_match(
    printed,
    "No pipeline is selected by this comparison.",
    fixed = TRUE
  )
})

test_that("survival comparison copies declared-horizon IPCW Brier", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    )
  )

  compared <- sb_compare_cv(
    task,
    pipelines,
    folds = 5L,
    seed = 20260819L
  )

  expect_identical(compared$status, "success")
  expect_identical(
    compared$comparison$complete_survival_probability_oof,
    c(TRUE, TRUE)
  )
  expect_identical(compared$comparison$prediction_horizon, c(7, 7))
  expect_identical(
    compared$comparison$ipcw_brier_status,
    c("success", "success")
  )
  expect_identical(
    compared$comparison$reference_ipcw_brier_status,
    c("success", "success")
  )
  expect_identical(
    compared$comparison$ipcw_brier_skill_status,
    c("success", "success")
  )
  expect_identical(
    compared$comparison$time_auc_status,
    c("success", "success")
  )
  expect_identical(
    compared$comparison$uno_c_status,
    c("success", "success")
  )
  expect_identical(
    compared$comparison$calibration_status,
    c("success", "success")
  )
  expect_identical(
    compared$comparison$mean_calibration_status,
    c("success", "success")
  )
  expect_identical(
    compared$comparison$calibration_slope_status,
    c("success", "success")
  )
  expected <- unname(vapply(
    compared$results,
    function(result) result$metrics$ipcw_brier$estimate,
    numeric(1)
  ))
  expect_equal(compared$comparison$ipcw_brier, expected, tolerance = 0)
  expected_reference <- unname(vapply(
    compared$results,
    function(result) result$metrics$ipcw_brier$reference$estimate,
    numeric(1)
  ))
  expected_skill <- unname(vapply(
    compared$results,
    function(result) result$metrics$ipcw_brier$skill_score$estimate,
    numeric(1)
  ))
  expect_equal(
    compared$comparison$reference_ipcw_brier,
    expected_reference,
    tolerance = 0
  )
  expect_equal(
    compared$comparison$ipcw_brier_skill_score,
    expected_skill,
    tolerance = 0
  )
  expected_time_auc <- unname(vapply(
    compared$results,
    function(result) result$metrics$time_auc$estimate,
    numeric(1)
  ))
  expected_time_auc_denominator <- unname(vapply(
    compared$results,
    function(result) result$metrics$time_auc$pair_denominator,
    numeric(1)
  ))
  expect_equal(
    compared$comparison$time_auc,
    expected_time_auc,
    tolerance = 0
  )
  expect_equal(
    compared$comparison$time_auc_pair_denominator,
    expected_time_auc_denominator,
    tolerance = 0
  )
  expected_uno_c <- unname(vapply(
    compared$results,
    function(result) result$metrics$uno_c$estimate,
    numeric(1)
  ))
  expected_uno_c_denominator <- unname(vapply(
    compared$results,
    function(result) result$metrics$uno_c$weighted_pair_denominator,
    numeric(1)
  ))
  expect_equal(compared$comparison$uno_c, expected_uno_c, tolerance = 0)
  expect_equal(
    compared$comparison$uno_c_weighted_pair_denominator,
    expected_uno_c_denominator,
    tolerance = 0
  )
  expected_calibration_ici <- unname(vapply(
    compared$results,
    function(result) result$metrics$calibration$ici,
    numeric(1)
  ))
  expected_calibration_e50 <- unname(vapply(
    compared$results,
    function(result) result$metrics$calibration$e50,
    numeric(1)
  ))
  expected_calibration_e90 <- unname(vapply(
    compared$results,
    function(result) result$metrics$calibration$e90,
    numeric(1)
  ))
  expect_identical(
    compared$comparison$calibration_ici,
    expected_calibration_ici
  )
  expect_identical(
    compared$comparison$calibration_e50,
    expected_calibration_e50
  )
  expect_identical(
    compared$comparison$calibration_e90,
    expected_calibration_e90
  )
  expected_mean_calibration <- unname(vapply(
    compared$results,
    function(result) result$metrics$weak_calibration$mean_calibration_ratio,
    numeric(1)
  ))
  expected_calibration_slope <- unname(vapply(
    compared$results,
    function(result) result$metrics$weak_calibration$calibration_slope,
    numeric(1)
  ))
  expect_identical(
    compared$comparison$mean_calibration_ratio,
    expected_mean_calibration
  )
  expect_identical(
    compared$comparison$calibration_slope,
    expected_calibration_slope
  )
  first_references <- compared$results[[1L]]$metrics$ipcw_brier$reference
  for (result in compared$results[-1L]) {
    observed_reference <- result$metrics$ipcw_brier$reference
    expect_identical(
      observed_reference$per_fold$reference_survival_horizon,
      first_references$per_fold$reference_survival_horizon
    )
    expect_identical(
      observed_reference$per_fold$weight_sum,
      first_references$per_fold$weight_sum
    )
    expect_identical(
      observed_reference$per_fold$weighted_loss_sum,
      first_references$per_fold$weighted_loss_sum
    )
    expect_identical(
      observed_reference$survival_curves,
      first_references$survival_curves
    )
  }
  printed <- paste(capture.output(print(compared)), collapse = "\n")
  expect_match(printed, "ipcw_brier", fixed = TRUE)
  expect_match(printed, "ipcw_brier_skill_score", fixed = TRUE)
  expect_match(printed, "time_auc", fixed = TRUE)
  expect_match(printed, "calibration_ici", fixed = TRUE)
  expect_match(printed, "calibration_e50", fixed = TRUE)
  expect_match(printed, "calibration_e90", fixed = TRUE)
  expect_match(printed, "mean_calibration_ratio", fixed = TRUE)
  expect_match(printed, "calibration_slope", fixed = TRUE)
})

test_that("survival comparison copies saved calibration without recomputing it", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    )
  )
  saved_results <- lapply(
    pipelines,
    function(pipeline) {
      sb_run_cv(task, pipeline, folds = 5L, seed = 20260819L)
    }
  )
  next_result <- 0L
  testthat::local_mocked_bindings(
    sb_run_cv = function(...) {
      next_result <<- next_result + 1L
      saved_results[[next_result]]
    },
    sb_survival_calibration_metric = function(...) {
      stop("comparison must not refit survival calibration")
    },
    sb_survival_weak_calibration_metric = function(...) {
      stop("comparison must not refit survival weak calibration")
    }
  )

  compared <- sb_compare_cv(
    task,
    pipelines,
    folds = 5L,
    seed = 20260819L
  )

  expect_identical(compared$status, "success")
  expect_identical(next_result, 2L)
  expect_identical(
    compared$comparison$pipeline_id,
    vapply(pipelines, `[[`, "", "id")
  )
  expect_false("winner" %in% names(compared$comparison))
  expect_identical(
    compared$comparison$calibration_ici,
    unname(vapply(saved_results, function(x) x$metrics$calibration$ici, 0))
  )
  expect_identical(
    compared$comparison$mean_calibration_status,
    unname(vapply(
      saved_results,
      function(x) x$metrics$weak_calibration$mean_calibration$status,
      ""
    ))
  )
  expect_identical(
    compared$comparison$mean_calibration_ratio,
    unname(vapply(
      saved_results,
      function(x) x$metrics$weak_calibration$mean_calibration_ratio,
      0
    ))
  )
  expect_identical(
    compared$comparison$calibration_slope_status,
    unname(vapply(
      saved_results,
      function(x) x$metrics$weak_calibration$slope_fit$status,
      ""
    ))
  )
  expect_identical(
    compared$comparison$calibration_slope,
    unname(vapply(
      saved_results,
      function(x) x$metrics$weak_calibration$calibration_slope,
      0
    ))
  )
})

test_that("all survival learners share the fold-local KM reference", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    ),
    sb_pipeline(
      learner = "glmnet_cox",
      inner_folds = 3L,
      inner_seed = 1701L
    )
  )

  compared <- sb_compare_cv(
    task,
    pipelines,
    folds = 5L,
    seed = 20260819L
  )
  reference <- compared$results[[1L]]$metrics$ipcw_brier$reference

  expect_identical(compared$status, "success")
  for (result in compared$results[-1L]) {
    observed <- result$metrics$ipcw_brier$reference
    expect_identical(observed$per_fold, reference$per_fold)
    expect_identical(observed$survival_curves, reference$survival_curves)
    expect_identical(observed$weighted_loss_sum, reference$weighted_loss_sum)
    expect_identical(observed$estimate, reference$estimate)
  }
})

test_that("unsupported survival horizon stays incomplete without metric fallback", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 1e6
  )

  compared <- sb_compare_cv(
    task,
    list(sb_pipeline()),
    folds = 5L,
    seed = 20260819L
  )

  expect_identical(compared$status, "incomplete")
  expect_identical(compared$comparison$status, "failed")
  expect_identical(compared$comparison$complete_oof, TRUE)
  expect_identical(
    compared$comparison$complete_survival_probability_oof,
    FALSE
  )
  expect_identical(compared$comparison$harrell_c_status, "success")
  expect_identical(compared$comparison$ipcw_brier_status, "incomplete")
  expect_identical(compared$comparison$time_auc_status, "incomplete")
  expect_identical(compared$comparison$calibration_status, "incomplete")
  expect_identical(
    compared$comparison$mean_calibration_status,
    "not_evaluated"
  )
  expect_identical(
    compared$comparison$calibration_slope_status,
    "not_fitted"
  )
  expect_identical(
    compared$comparison$reference_ipcw_brier_status,
    "incomplete"
  )
  expect_identical(
    compared$comparison$ipcw_brier_skill_status,
    "incomplete"
  )
  expect_true(is.na(compared$comparison$ipcw_brier))
  expect_true(is.na(compared$comparison$reference_ipcw_brier))
  expect_true(is.na(compared$comparison$ipcw_brier_skill_score))
  expect_true(is.na(compared$comparison$time_auc))
  expect_true(is.na(compared$comparison$calibration_ici))
  expect_true(is.na(compared$comparison$calibration_e50))
  expect_true(is.na(compared$comparison$calibration_e90))
  expect_true(is.na(compared$comparison$mean_calibration_ratio))
  expect_true(is.na(compared$comparison$calibration_slope))
  expect_true(is.finite(compared$comparison$harrell_c))
})

test_that("comparison copies unavailable Uno evidence without ranking", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 7
  )
  original_uno <- sb_survival_uno_c_metric
  testthat::local_mocked_bindings(
    sb_survival_uno_c_metric = function(...) {
      metric <- original_uno(...)
      metric$status <- "incomplete"
      metric$estimate <- NA_real_
      metric
    }
  )

  compared <- sb_compare_cv(
    task,
    list(sb_pipeline()),
    folds = 5L,
    seed = 20260819L
  )

  expect_identical(compared$status, "incomplete")
  expect_identical(compared$comparison$status, "failed")
  expect_identical(compared$comparison$harrell_c_status, "success")
  expect_identical(compared$comparison$ipcw_brier_status, "success")
  expect_identical(compared$comparison$time_auc_status, "success")
  expect_identical(compared$comparison$uno_c_status, "incomplete")
  expect_true(is.na(compared$comparison$uno_c))
  expect_false("winner" %in% names(compared$comparison))
})

test_that("comparison keeps a failed pipeline visible without fallback", {
  d <- make_sb_fixture(n = 60, seed = 20260816)
  d$gene_b <- d$gene_a
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month"
  )
  pipelines <- list(
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 1L
    ),
    sb_pipeline()
  )

  comparison <- sb_compare_cv(task, pipelines, folds = 5L, seed = 101L)

  expect_identical(comparison$status, "incomplete")
  expect_identical(
    comparison$comparison$pipeline_id,
    c("univ_cox_top_k__coxph__k1", "none__coxph")
  )
  expect_identical(comparison$comparison$status, c("success", "failed"))
  expect_identical(comparison$comparison$complete_oof, c(TRUE, FALSE))
  expect_identical(
    comparison$comparison$harrell_c_status,
    c("success", "incomplete")
  )
  expect_identical(is.finite(comparison$comparison$harrell_c), c(TRUE, FALSE))
  failed_result <- comparison$results[["none__coxph"]]
  expect_identical(all(is.na(failed_result$runs$learner_used)), TRUE)
  expect_identical(all(is.na(failed_result$predictions$risk)), TRUE)
})

test_that("comparison admits different alphas on identical outer folds", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260816)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipelines <- list(
    sb_pipeline(
      learner = "glmnet_cox",
      glmnet_alpha = 0.3,
      inner_folds = 3L,
      inner_seed = 101L
    ),
    sb_pipeline(
      learner = "glmnet_cox",
      glmnet_alpha = 0.7,
      inner_folds = 3L,
      inner_seed = 101L
    )
  )

  comparison <- sb_compare_cv(task, pipelines, folds = 5L, seed = 101L)

  expect_identical(comparison$status, "success")
  expect_identical(
    comparison$comparison$pipeline_id,
    c(
      "none__glmnet_cox__alpha0.3__inner3__seed101",
      "none__glmnet_cox__alpha0.7__inner3__seed101"
    )
  )
  expect_length(comparison$results, 2L)
  assignments <- lapply(
    comparison$results,
    function(x) x$predictions[c("patient_id", "fold")]
  )
  expect_identical(assignments[[1L]], assignments[[2L]])
})

test_that("comparison requires a usable collection of unique pipelines", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )

  empty_error <- capture_sb_error(
    sb_compare_cv(task, list(), folds = 2L, seed = 1L)
  )
  expect_s3_class(empty_error, "stormbreaker_bad_pipeline_collection")

  member_error <- capture_sb_error(
    sb_compare_cv(task, list("not a pipeline"), folds = 2L, seed = 1L)
  )
  expect_s3_class(member_error, "stormbreaker_bad_pipeline_collection")

  duplicate_error <- capture_sb_error(
    sb_compare_cv(
      task,
      list(sb_pipeline(), sb_pipeline()),
      folds = 2L,
      seed = 1L
    )
  )
  expect_s3_class(duplicate_error, "stormbreaker_duplicate_pipeline_id")
})

test_that("binary comparison uses identical outer folds and saved AUC", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 60)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes"),
    sb_pipeline(
      learner = "glmnet_logistic",
      glmnet_alpha = 0,
      inner_folds = 3L,
      inner_seed = 101L
    ),
    sb_pipeline(
      learner = "glmnet_logistic",
      glmnet_alpha = 1,
      inner_folds = 3L,
      inner_seed = 101L
    )
  )
  comparison <- sb_compare_cv(
    task, pipelines, folds = 3L, seed = 20260819L
  )

  expect_identical(comparison$status, "success")
  expect_identical(comparison$task_type, "binary")
  expect_identical(
    comparison$comparison$pipeline_id,
    vapply(pipelines, `[[`, "", "id")
  )
  expect_named(
    comparison$comparison,
    c(
      "pipeline_id", "selector", "learner", "status", "complete_oof",
      "roc_auc_status", "roc_auc", "comparable_pairs",
      "average_precision_status", "average_precision", "brier_status",
      "brier_score", "n_scored", "reference_brier_status",
      "reference_brier", "brier_skill_status", "brier_skill_score",
      "log_loss_status", "log_loss", "reference_log_loss_status",
      "reference_log_loss", "log_loss_skill_status", "log_loss_skill_score",
      "calibration_status", "calibration_intercept", "calibration_slope",
      "moderate_calibration_status", "moderate_calibration_ici",
      "moderate_calibration_e50", "moderate_calibration_e90",
      "selection_stability_status", "mean_pairwise_jaccard",
      "n_successful_selection_folds",
      "n_folds", "n_successful_folds",
      "glmnet_active_set_stability_status",
      "glmnet_active_set_mean_pairwise_jaccard",
      "n_successful_glmnet_active_set_folds", "prediction_variation"
    ),
    ignore.order = FALSE
  )
  assignments <- lapply(
    comparison$results,
    function(x) x$predictions[c("patient_id", "fold")]
  )
  expect_identical(assignments[[1L]], assignments[[2L]])
  expect_identical(assignments[[1L]], assignments[[3L]])
  expect_identical(assignments[[1L]], assignments[[4L]])
  expected_auc <- vapply(
    comparison$results,
    function(x) x$metrics$roc_auc$estimate,
    numeric(1)
  )
  expect_equal(
    comparison$comparison$roc_auc,
    unname(expected_auc),
    tolerance = 0
  )
  expected_average_precision <- vapply(
    comparison$results,
    function(x) x$metrics$average_precision$estimate,
    numeric(1)
  )
  expect_equal(
    comparison$comparison$average_precision,
    unname(expected_average_precision),
    tolerance = 0
  )
  expected_brier <- vapply(
    comparison$results,
    function(x) x$metrics$brier_score$estimate,
    numeric(1)
  )
  manual_brier <- vapply(
    comparison$results,
    function(x) {
      mean((x$predictions$outcome - x$predictions$probability)^2)
    },
    numeric(1)
  )
  expect_equal(expected_brier, manual_brier, tolerance = 0)
  expect_equal(
    comparison$comparison$brier_score,
    unname(expected_brier),
    tolerance = 0
  )
  expected_reference <- vapply(
    comparison$results,
    function(x) x$metrics$brier_score$reference$estimate,
    numeric(1)
  )
  expected_skill <- vapply(
    comparison$results,
    function(x) x$metrics$brier_score$skill_score$estimate,
    numeric(1)
  )
  expect_equal(
    comparison$comparison$reference_brier,
    unname(expected_reference),
    tolerance = 0
  )
  expect_length(unique(comparison$comparison$reference_brier), 1L)
  expect_equal(
    comparison$comparison$brier_skill_score,
    unname(expected_skill),
    tolerance = 0
  )
  expected_log_loss <- vapply(
    comparison$results,
    function(x) x$metrics$log_loss$estimate,
    numeric(1)
  )
  expected_reference_log_loss <- vapply(
    comparison$results,
    function(x) x$metrics$log_loss$reference$estimate,
    numeric(1)
  )
  expected_log_loss_skill <- vapply(
    comparison$results,
    function(x) x$metrics$log_loss$skill_score$estimate,
    numeric(1)
  )
  expect_equal(
    comparison$comparison$log_loss,
    unname(expected_log_loss),
    tolerance = 0
  )
  expect_equal(
    comparison$comparison$reference_log_loss,
    unname(expected_reference_log_loss),
    tolerance = 0
  )
  expect_equal(
    comparison$comparison$log_loss_skill_score,
    unname(expected_log_loss_skill),
    tolerance = 0
  )
  expected_calibration_intercept <- vapply(
    comparison$results,
    function(x) x$metrics$calibration$calibration_intercept,
    numeric(1)
  )
  expected_calibration_slope <- vapply(
    comparison$results,
    function(x) x$metrics$calibration$calibration_slope,
    numeric(1)
  )
  expect_equal(
    comparison$comparison$calibration_intercept,
    unname(expected_calibration_intercept),
    tolerance = 0
  )
  expect_equal(
    comparison$comparison$calibration_slope,
    unname(expected_calibration_slope),
    tolerance = 0
  )
  expected_moderate_status <- vapply(
    comparison$results,
    function(x) x$metrics$moderate_calibration$status,
    character(1)
  )
  expect_identical(
    comparison$comparison$moderate_calibration_status,
    unname(expected_moderate_status)
  )
  for (metric in c("ici", "e50", "e90")) {
    expected <- unname(vapply(
      comparison$results,
      function(x) x$metrics$moderate_calibration[[metric]],
      numeric(1)
    ))
    expect_equal(
      comparison$comparison[[paste0("moderate_calibration_", metric)]],
      expected,
      tolerance = 0
    )
  }
  expect_identical(
    comparison$comparison$selection_stability_status,
    rep("not_applicable", length(pipelines))
  )
  expect_true(all(is.na(comparison$comparison$mean_pairwise_jaccard)))
  expect_identical(
    comparison$comparison$n_successful_selection_folds,
    rep(0L, length(pipelines))
  )
  expected_active_status <- unname(vapply(
    comparison$results,
    function(result) result$glmnet_active_set_stability$status,
    character(1)
  ))
  expect_identical(
    comparison$comparison$glmnet_active_set_stability_status,
    expected_active_status
  )
  expect_identical(
    expected_active_status[1:3],
    rep("not_applicable", 3L)
  )
  expect_match(expected_active_status[[4L]], "^(success|not_evaluable)$")
  expect_equal(
    comparison$comparison$glmnet_active_set_mean_pairwise_jaccard,
    unname(vapply(
      comparison$results,
      function(result) {
        result$glmnet_active_set_stability$mean_pairwise_jaccard
      },
      numeric(1)
    )),
    tolerance = 0
  )
  expect_identical(
    comparison$comparison$n_successful_glmnet_active_set_folds,
    c(0L, 0L, 0L, 3L)
  )
  printed <- paste(capture.output(print(comparison)), collapse = "\n")
  expect_match(printed, "average_precision", fixed = TRUE)
  expect_match(printed, "brier_score", fixed = TRUE)
  expect_match(printed, "brier_skill_score", fixed = TRUE)
  expect_match(printed, "log_loss", fixed = TRUE)
  expect_match(printed, "log_loss_skill_score", fixed = TRUE)
  expect_match(printed, "calibration_intercept", fixed = TRUE)
  expect_match(printed, "calibration_slope", fixed = TRUE)
  expect_match(printed, "moderate_calibration_ici", fixed = TRUE)
  expect_match(printed, "moderate_calibration_e50", fixed = TRUE)
  expect_match(printed, "moderate_calibration_e90", fixed = TRUE)
  expect_match(printed, "No pipeline is selected", fixed = TRUE)
})

test_that("binary comparison copies saved moderate calibration without refitting", {
  d <- make_sb_binary_fixture(n = 60L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )
  saved_results <- lapply(
    pipelines,
    function(pipeline) {
      sb_run_cv(task, pipeline, folds = 5L, seed = 20260820L)
    }
  )
  saved_results[[1L]]$metrics$moderate_calibration$ici <- 0.123456789
  next_result <- 0L
  testthat::local_mocked_bindings(
    sb_run_cv = function(...) {
      next_result <<- next_result + 1L
      saved_results[[next_result]]
    },
    sb_binary_moderate_calibration_metric = function(...) {
      stop("comparison must not refit binary moderate calibration")
    }
  )

  compared <- sb_compare_cv(
    task,
    pipelines,
    folds = 5L,
    seed = 20260820L
  )

  expect_identical(next_result, 2L)
  expect_identical(
    compared$comparison$moderate_calibration_ici[[1L]],
    0.123456789
  )
  expect_identical(
    compared$comparison$moderate_calibration_e50,
    unname(vapply(
      saved_results,
      function(x) x$metrics$moderate_calibration$e50,
      numeric(1)
    ))
  )
  expect_false("winner" %in% names(compared$comparison))
})

test_that("binary comparison continues after an earlier failed logistic", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 60)
  d$gene_b <- d$gene_a
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(
      learner = "glmnet_logistic",
      inner_folds = 3L,
      inner_seed = 1L
    )
  )
  comparison <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1L)

  expect_identical(comparison$status, "incomplete")
  expect_identical(
    names(comparison$results),
    vapply(pipelines, `[[`, character(1L), "id")
  )
  expect_identical(comparison$comparison$status, c("failed", "success"))
  expect_identical(
    comparison$comparison$roc_auc_status,
    c("incomplete", "success")
  )
  expect_identical(
    comparison$comparison$average_precision_status,
    c("incomplete", "success")
  )
  expect_identical(
    comparison$comparison$brier_status,
    c("incomplete", "success")
  )
  expect_identical(
    comparison$comparison$reference_brier_status,
    c("success", "success")
  )
  expect_identical(
    comparison$comparison$brier_skill_status,
    c("incomplete", "success")
  )
  expect_identical(
    comparison$comparison$log_loss_status,
    c("incomplete", "success")
  )
  expect_identical(
    comparison$comparison$reference_log_loss_status,
    c("success", "success")
  )
  expect_identical(
    comparison$comparison$log_loss_skill_status,
    c("incomplete", "success")
  )
  expect_identical(
    comparison$comparison$calibration_status,
    c("incomplete", "success")
  )
  expect_identical(
    comparison$comparison$moderate_calibration_status,
    c("incomplete", "success")
  )
  failed <- comparison$results[[pipelines[[1L]]$id]]
  succeeded <- comparison$results[[pipelines[[2L]]$id]]
  expect_identical(all(is.na(failed$runs$learner_used)), TRUE)
  expect_identical(all(is.na(failed$predictions$probability)), TRUE)
  expect_identical(succeeded$status, "success")
  expect_true(all(is.finite(succeeded$predictions$probability)))
})

test_that("binary comparison admits fold-local AUC selection without ranking", {
  d <- make_sb_high_dim_binary_fixture(n = 60L, p = 8L, seed = 1701L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "logistic",
      selector_n_features = 2L
    ),
    sb_pipeline(learner = "gaussian_naive_bayes"),
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "gaussian_naive_bayes",
      selector_n_features = 2L
    )
  )

  comparison <- sb_compare_cv(
    task, pipelines, folds = 3L, seed = 1801L
  )

  expect_identical(comparison$status, "success")
  expect_identical(
    comparison$comparison$pipeline_id,
    vapply(pipelines, `[[`, "", "id")
  )
  fold_assignments <- lapply(
    comparison$results,
    function(result) result$predictions[c("patient_id", "fold")]
  )
  for (index in 2:length(fold_assignments)) {
    expect_identical(fold_assignments[[1L]], fold_assignments[[index]])
  }
  expect_equal(
    comparison$comparison$roc_auc,
    unname(vapply(
      comparison$results,
      function(result) result$metrics$roc_auc$estimate,
      numeric(1)
    )),
    tolerance = 0
  )
  expect_equal(
    comparison$comparison$brier_score,
    unname(vapply(
      comparison$results,
      function(result) result$metrics$brier_score$estimate,
      numeric(1)
    )),
    tolerance = 0
  )
  expect_length(comparison$results[[2L]]$selections, 3L)
  expect_length(comparison$results[[4L]]$selections, 3L)
  expect_identical(
    comparison$comparison$selection_stability_status,
    c("not_applicable", "success", "not_applicable", "success")
  )
  expect_identical(
    is.finite(comparison$comparison$mean_pairwise_jaccard),
    c(FALSE, TRUE, FALSE, TRUE)
  )
  expect_identical(
    comparison$comparison$n_successful_selection_folds,
    c(0L, 3L, 0L, 3L)
  )
  expect_equal(
    comparison$comparison$mean_pairwise_jaccard,
    unname(vapply(
      comparison$results,
      function(result) result$selection_stability$mean_pairwise_jaccard,
      numeric(1)
    )),
    tolerance = 0
  )
  selected_audit_matches <- Map(
    function(audit, selection) {
      length(unique(audit$values$feature)) == 2L &&
        setequal(unique(audit$values$feature), selection$selected_features)
    },
    comparison$results[[4L]]$distribution_parameters,
    comparison$results[[4L]]$selections
  )
  expect_true(all(unlist(selected_audit_matches, use.names = FALSE)))
  printed <- paste(capture.output(print(comparison)), collapse = "\n")
  expect_match(printed, "No pipeline is selected", fixed = TRUE)
})

test_that("fixed LightGBM comparison reuses folds with and without Top-K", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  d <- make_sb_high_dim_binary_fixture(n = 48L, p = 10L, seed = 3901L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  fixed_args <- list(
    lightgbm_nrounds = 6L,
    lightgbm_learning_rate = 0.1,
    lightgbm_num_leaves = 7L,
    lightgbm_min_data_in_leaf = 3L,
    lightgbm_seed = 4001L
  )
  pipelines <- list(
    do.call(
      sb_pipeline,
      c(list(learner = "lightgbm_binary"), fixed_args)
    ),
    do.call(
      sb_pipeline,
      c(
        list(
          selector = "univ_binary_auc_top_k",
          learner = "lightgbm_binary",
          selector_n_features = 3L
        ),
        fixed_args
      )
    )
  )

  comparison <- sb_compare_cv(
    task, pipelines, folds = 3L, seed = 4101L
  )

  expect_identical(comparison$status, "success")
  expect_identical(
    comparison$comparison$pipeline_id,
    vapply(pipelines, `[[`, "", "id")
  )
  fold_assignments <- lapply(
    comparison$results,
    function(result) result$predictions[c("patient_id", "fold")]
  )
  expect_identical(fold_assignments[[1L]], fold_assignments[[2L]])
  expect_length(comparison$results[[1L]]$tree_models, 3L)
  expect_length(comparison$results[[2L]]$tree_models, 3L)
  expect_length(comparison$results[[2L]]$selections, 3L)
  for (fold in seq_len(3L)) {
    expect_setequal(
      comparison$results[[2L]]$tree_models[[fold]]$features,
      comparison$results[[2L]]$selections[[fold]]$selected_features
    )
  }
  printed <- paste(capture.output(print(comparison)), collapse = "\n")
  expect_match(printed, "No pipeline is selected", fixed = TRUE)
})

test_that("fixed CatBoost comparison reuses folds with and without Top-K", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_high_dim_binary_fixture(n = 48L, p = 10L, seed = 5301L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  fixed_args <- list(
    catboost_iterations = 6L,
    catboost_learning_rate = 0.1,
    catboost_depth = 3L,
    catboost_l2_leaf_reg = 3,
    catboost_seed = 5401L
  )
  pipelines <- list(
    do.call(
      sb_pipeline,
      c(list(learner = "catboost_binary"), fixed_args)
    ),
    do.call(
      sb_pipeline,
      c(
        list(
          selector = "univ_binary_auc_top_k",
          learner = "catboost_binary",
          selector_n_features = 3L
        ),
        fixed_args
      )
    )
  )

  comparison <- sb_compare_cv(
    task, pipelines, folds = 3L, seed = 5501L
  )

  expect_identical(comparison$status, "success")
  expect_identical(
    comparison$comparison$pipeline_id,
    vapply(pipelines, `[[`, "", "id")
  )
  fold_assignments <- lapply(
    comparison$results,
    function(result) result$predictions[c("patient_id", "fold")]
  )
  expect_identical(fold_assignments[[1L]], fold_assignments[[2L]])
  expect_length(comparison$results[[1L]]$tree_models, 3L)
  expect_length(comparison$results[[2L]]$tree_models, 3L)
  expect_length(comparison$results[[2L]]$selections, 3L)
  for (fold in seq_len(3L)) {
    expect_setequal(
      comparison$results[[2L]]$tree_models[[fold]]$features,
      comparison$results[[2L]]$selections[[fold]]$selected_features
    )
  }
  printed <- paste(capture.output(print(comparison)), collapse = "\n")
  expect_match(printed, "No pipeline is selected", fixed = TRUE)
})

test_that("binary comparison copies declared decision curves without selection", {
  d <- make_sb_binary_fixture(n = 60L, seed = 5601L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c"),
    decision_thresholds = c(0.2, 0.5, 0.8)
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )

  compared <- sb_compare_cv(
    task,
    pipelines,
    folds = 3L,
    seed = 5701L
  )

  expect_identical(
    names(compared$decision_curves),
    vapply(pipelines, `[[`, character(1), "id")
  )
  for (index in seq_along(pipelines)) {
    expect_identical(
      compared$decision_curves[[index]],
      compared$results[[index]]$metrics$decision_curve
    )
  }
  expect_false(any(grepl(
    "best|winner|selected_threshold",
    names(compared),
    ignore.case = TRUE
  )))
  expect_false(any(grepl(
    "decision",
    names(compared$comparison),
    ignore.case = TRUE
  )))
})

test_that("binary comparison omits decision curves when undeclared", {
  d <- make_sb_binary_fixture(n = 60L, seed = 5801L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )
  compared <- sb_compare_cv(
    task,
    list(
      sb_pipeline(learner = "logistic"),
      sb_pipeline(learner = "gaussian_naive_bayes")
    ),
    folds = 3L,
    seed = 5901L
  )

  expect_false("decision_curves" %in% names(compared))
})

test_that("survival comparison copies declared curves without selection", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7,
    decision_thresholds = c(0.20, 0.40, 0.60)
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    )
  )
  compared <- sb_compare_cv(
    task,
    pipelines,
    folds = 5L,
    seed = 20260819L
  )

  expect_identical(compared$status, "success")
  expect_named(compared$decision_curves, vapply(pipelines, `[[`, "", "id"))
  for (index in seq_along(pipelines)) {
    expect_identical(
      compared$decision_curves[[index]],
      compared$results[[index]]$metrics$decision_curve
    )
  }
  expect_false(any(grepl(
    "decision|winner|best|threshold",
    names(compared$comparison),
    ignore.case = TRUE
  )))
})

test_that("survival comparison reuses one patient-keyed outer assignment", {
  d <- make_sb_fixture(n = 30L, seed = 6301L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    )
  )
  generated <- sb_make_stratified_folds(task, 3L, 6302L)
  supplied <- setNames(generated$assignment, d$patient_id)
  supplied <- supplied[c(30L, seq_len(29L))]

  first <- sb_compare_cv(task, pipelines, folds = supplied, seed = 1L)
  second <- sb_compare_cv(task, pipelines, folds = supplied, seed = 999L)

  expect_identical(first, second)
  expect_identical(first$fold_source, "preassigned_patient_id")
  expect_identical(first$seed, NA_integer_)
  expect_identical(first$folds, 3L)
  expect_identical(
    first$fold_assignment,
    data.frame(
      patient_id = d$patient_id,
      fold = as.integer(supplied[match(d$patient_id, names(supplied))]),
      stringsAsFactors = FALSE
    )
  )
  for (index in seq_along(pipelines)) {
    direct <- sb_run_cv(
      task,
      pipelines[[index]],
      folds = supplied,
      seed = 123L
    )
    expect_identical(first$results[[index]], direct)
    expect_identical(
      first$results[[index]]$provenance$fold_assignment,
      first$fold_assignment
    )
  }
  printed <- paste(capture.output(print(first)), collapse = "\n")
  expect_match(printed, "preassigned patient-ID folds", fixed = TRUE)
})

test_that("binary comparison reuses one patient-keyed outer assignment", {
  d <- make_sb_binary_fixture(n = 30L, seed = 6311L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )
  generated <- sb_make_stratified_folds(task, 3L, 6312L)
  supplied <- setNames(generated$assignment, d$patient_id)
  supplied <- supplied[rev(seq_along(supplied))]

  first <- sb_compare_cv(task, pipelines, folds = supplied, seed = 1L)
  second <- sb_compare_cv(task, pipelines, folds = supplied, seed = 999L)

  expect_identical(first, second)
  expect_identical(first$fold_source, "preassigned_patient_id")
  expect_identical(first$seed, NA_integer_)
  expect_identical(first$folds, 3L)
  expect_identical(
    first$fold_assignment,
    data.frame(
      patient_id = d$patient_id,
      fold = as.integer(supplied[match(d$patient_id, names(supplied))]),
      stringsAsFactors = FALSE
    )
  )
  for (index in seq_along(pipelines)) {
    direct <- sb_run_cv(
      task,
      pipelines[[index]],
      folds = supplied,
      seed = 123L
    )
    expect_identical(first$results[[index]], direct)
    expect_identical(
      first$results[[index]]$provenance$fold_assignment,
      first$fold_assignment
    )
  }
  printed <- paste(capture.output(print(first)), collapse = "\n")
  expect_match(printed, "preassigned patient-ID folds", fixed = TRUE)
})

# ---- compare-repeated-cv ----

test_that("repeated binary comparison preserves every seeded comparison", {
  d <- make_sb_binary_fixture(n = 48L, seed = 7101L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )
  seeds <- c(7102L, 7103L)

  observed <- sb_compare_repeated_cv(
    task,
    pipelines,
    folds = 3L,
    seeds = seeds
  )

  expect_s3_class(observed, "stormbreaker_repeated_cv_comparison")
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$n_patients, 48L)
  expect_identical(observed$n_repeats, 2L)
  expect_identical(names(observed$repetitions), c("repeat_1", "repeat_2"))
  expect_identical(
    observed$comparison$repeat_id,
    rep(1:2, each = length(pipelines))
  )
  expect_identical(
    observed$comparison$seed,
    rep(seeds, each = length(pipelines))
  )
  expect_identical(
    observed$comparison$fold_source,
    rep("generated_stratified", 4L)
  )
  expect_identical(
    observed$fold_assignments$repeat_id,
    rep(1:2, each = nrow(d))
  )

  for (repeat_id in seq_along(seeds)) {
    direct <- sb_compare_cv(
      task,
      pipelines,
      folds = 3L,
      seed = seeds[[repeat_id]]
    )
    expect_identical(observed$repetitions[[repeat_id]], direct)
    expected_assignment <- data.frame(
      repeat_id = repeat_id,
      direct$fold_assignment,
      stringsAsFactors = FALSE
    )
    actual_assignment <- observed$fold_assignments[
      observed$fold_assignments$repeat_id == repeat_id,
      ,
      drop = FALSE
    ]
    rownames(actual_assignment) <- NULL
    expect_identical(actual_assignment, expected_assignment)
  }

  expect_length(
    grep(
      "mean|standard|interval|rank|winner|best",
      names(observed),
      ignore.case = TRUE,
      value = TRUE
    ),
    0L
  )
  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(
    printed,
    "No cross-repeat mean, confidence interval, ranking, or winner",
    fixed = TRUE
  )
})

test_that("repeated survival comparison preserves explicit fold plans", {
  d <- make_sb_fixture(n = 36L, seed = 7201L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    )
  )
  first <- sb_make_stratified_folds(task, 3L, 7202L)
  second <- sb_make_stratified_folds(task, 3L, 7203L)
  fold_plans <- list(
    setNames(first$assignment, d$patient_id),
    rev(setNames(second$assignment, d$patient_id))
  )

  observed <- sb_compare_repeated_cv(
    task,
    pipelines,
    folds = fold_plans,
    seeds = NULL
  )

  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "survival")
  expect_identical(observed$n_repeats, 2L)
  expect_identical(observed$comparison$seed, rep(NA_integer_, 4L))
  expect_identical(
    observed$comparison$fold_source,
    rep("preassigned_patient_id", 4L)
  )

  for (repeat_id in seq_along(fold_plans)) {
    direct <- sb_compare_cv(
      task,
      pipelines,
      folds = fold_plans[[repeat_id]],
      seed = 999L
    )
    expect_identical(observed$repetitions[[repeat_id]], direct)
    expected_fold <- as.integer(fold_plans[[repeat_id]][match(
      d$patient_id,
      names(fold_plans[[repeat_id]])
    )])
    actual <- observed$fold_assignments[
      observed$fold_assignments$repeat_id == repeat_id,
      ,
      drop = FALSE
    ]
    expect_identical(actual$patient_id, d$patient_id)
    expect_identical(actual$fold, expected_fold)
  }
})

test_that("repeated comparison preserves incomplete pipeline evidence", {
  d <- make_sb_fixture(n = 48L, seed = 7251L)
  d$gene_b <- d$gene_a
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month"
  )
  pipelines <- list(
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 1L
    ),
    sb_pipeline()
  )

  observed <- sb_compare_repeated_cv(
    task,
    pipelines,
    folds = 3L,
    seeds = c(7252L, 7253L)
  )

  expect_identical(observed$status, "incomplete")
  expect_identical(
    unname(vapply(
      observed$repetitions,
      `[[`,
      character(1L),
      "status"
    )),
    rep("incomplete", 2L)
  )
  expect_identical(
    observed$comparison$status,
    rep(c("success", "failed"), 2L)
  )
  expect_identical(
    observed$comparison$complete_oof,
    rep(c(TRUE, FALSE), 2L)
  )
})

test_that("repeated comparison keeps generated and explicit modes unambiguous", {
  d <- make_sb_binary_fixture(n = 24L, seed = 7301L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(sb_pipeline(learner = "logistic"))
  generated <- sb_make_stratified_folds(task, 3L, 7302L)
  plan <- setNames(generated$assignment, d$patient_id)

  missing_seeds <- capture_sb_error(sb_compare_repeated_cv(
    task, pipelines, folds = 3L, seeds = NULL
  ))
  conflicting_seeds <- capture_sb_error(sb_compare_repeated_cv(
    task, pipelines, folds = list(plan), seeds = 7303L
  ))
  bare_plan <- capture_sb_error(sb_compare_repeated_cv(
    task, pipelines, folds = plan, seeds = NULL
  ))
  invalid_seed <- capture_sb_error(sb_compare_repeated_cv(
    task, pipelines, folds = 3L, seeds = c(7303L, NA_real_)
  ))
  empty_plans <- capture_sb_error(sb_compare_repeated_cv(
    task, pipelines, folds = list(), seeds = NULL
  ))
  scalar_plan <- capture_sb_error(sb_compare_repeated_cv(
    task, pipelines, folds = list(3L), seeds = NULL
  ))
  mixed_plans <- capture_sb_error(sb_compare_repeated_cv(
    task, pipelines, folds = list(plan, 3L), seeds = NULL
  ))

  expect_s3_class(missing_seeds, "stormbreaker_bad_repeat_seeds")
  expect_s3_class(
    conflicting_seeds,
    "stormbreaker_repeated_fold_seed_conflict"
  )
  expect_s3_class(bare_plan, "stormbreaker_bad_repeated_folds")
  expect_s3_class(invalid_seed, "stormbreaker_bad_repeat_seeds")
  expect_s3_class(empty_plans, "stormbreaker_bad_repeated_folds")
  expect_s3_class(scalar_plan, "stormbreaker_bad_repeated_folds")
  expect_s3_class(mixed_plans, "stormbreaker_bad_repeated_folds")
})

test_that("repeated comparison preserves complete generated stack runs", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 48L, seed = 7351L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 7352L)
  seeds <- c(7353L, 7354L)

  observed <- sb_compare_repeated_cv(
    task,
    stack,
    folds = 3L,
    seeds = seeds
  )

  expect_identical(observed$status, "success")
  expect_identical(observed$strategy_type, "stack")
  expect_identical(observed$n_repeats, 2L)
  expect_identical(
    observed$comparison$pipeline_id,
    rep(c(stack$pipeline_ids, stack$id), 2L)
  )
  for (repeat_id in seq_along(seeds)) {
    direct <- sb_run_stack_cv(
      task,
      stack,
      folds = 3L,
      seed = seeds[[repeat_id]]
    )
    expect_identical(observed$repetitions[[repeat_id]], direct)
    actual_assignment <- observed$fold_assignments[
      observed$fold_assignments$repeat_id == repeat_id,
      c("patient_id", "fold")
    ]
    rownames(actual_assignment) <- NULL
    expect_identical(
      actual_assignment,
      direct$provenance$fold_assignment
    )
  }
})

test_that("repeated stack comparison keeps explicit patient fold plans", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 42L, seed = 7361L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 7362L)
  generated <- sb_make_stratified_folds(task, 3L, 7363L)
  plan <- rev(stats::setNames(generated$assignment, data$patient_id))

  observed <- sb_compare_repeated_cv(
    task,
    stack,
    folds = list(plan),
    seeds = NULL
  )
  expected <- sb_run_stack_cv(task, stack, folds = plan, seed = 1L)

  expect_identical(observed$status, "success")
  expect_identical(observed$repetitions[[1L]], expected)
  expect_identical(observed$comparison$seed, rep(NA_integer_, 3L))
  expect_identical(
    observed$comparison$fold_source,
    rep("preassigned_patient_id", 3L)
  )
  expect_identical(
    observed$fold_assignments$patient_id,
    data$patient_id
  )
  expect_identical(
    observed$fold_assignments$fold,
    as.integer(plan[match(data$patient_id, names(plan))])
  )
})

# ---- comparison-progress ----

test_that("CV progress describes actual folds without changing results", {
  d <- sb_example_binary(45, 10)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  pipelines <- list(sb_pipeline(learner = "logistic"))
  silent <- sb_compare_cv(task, pipelines, folds = 3, seed = 5, progress = FALSE)
  messages <- capture.output(
    verbose <- sb_compare_cv(task, pipelines, folds = 3, seed = 5, progress = TRUE),
    type = "message"
  )
  expect_identical(verbose, silent)
  expect_match(paste(messages, collapse = " "), "CV model 1/1")
  expect_match(paste(messages, collapse = " "), "Fold 3/3")
  expect_match(paste(messages, collapse = " "), "complete")
})

test_that("cohort progress includes final fit, prediction and evaluation", {
  d <- sb_example_binary(45, 10)
  v <- sb_example_binary(25, 11)
  v$patient_id <- paste0("v", v$patient_id)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  val <- sb_binary_task(v, "patient_id", "outcome", c("gene_a", "gene_b"))
  messages <- capture.output(
    result <- sb_compare_cohorts(task, list(Validation = val),
      roles = c(Validation = "validation"), models = "logistic", folds = 3,
      view = FALSE, progress = TRUE), type = "message"
  )
  expect_s3_class(result, "stormbreaker_cohort_comparison")
  text <- paste(messages, collapse = " ")
  expect_match(text, "Final fit 1/1")
  expect_match(text, "Predict Validation")
  expect_match(text, "Evaluate Validation")
})

test_that("survival progress preserves results and invalid flags are rejected", {
  d <- sb_example_survival(45, 3)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month")
  silent <- sb_run_cv(task, folds = 3, progress = FALSE)
  messages <- capture.output(verbose <- sb_run_cv(task, folds = 3, progress = TRUE), type = "message")
  expect_identical(verbose, silent)
  expect_match(paste(messages, collapse = " "), "Fold 3/3 complete")
  expect_error(sb_run_cv(task, progress = NA), "TRUE or FALSE", class = "stormbreaker_bad_progress")
})

# ---- repeated-cv-parallel ----

test_that("two workers preserve repeated binary comparisons and caller RNG", {
  d <- sb_example_binary(60L, 1101L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  pp <- list(sb_pipeline("none", "logistic"), sb_pipeline("none", "gaussian_naive_bayes"))
  withr::local_seed(119L)
  before <- .Random.seed
  serial <- sb_compare_repeated_cv(task, pp, 3L, c(19L, 20L, 21L))
  actual <- sb_compare_repeated_cv(task, pp, 3L, c(19L, 20L, 21L), workers = 2L)
  expect_identical(actual, serial)
  expect_identical(.Random.seed, before)
  plans <- lapply(serial$repetitions, function(x) setNames(x$fold_assignment$fold, x$fold_assignment$patient_id))
  expect_identical(sb_compare_repeated_cv(task, pp, plans, NULL, workers = 2L),
    sb_compare_repeated_cv(task, pp, plans, NULL))
})

test_that("two workers preserve survival horizon and failed routes", {
  d <- sb_example_survival(72L, 11032L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"),
    "month", prediction_horizon = 7)
  pp <- list(sb_pipeline("none", "coxph"), sb_pipeline("univ_cox_top_k", "coxph", selector_n_features = 2L))
  expect_identical(sb_compare_repeated_cv(task, pp, 3L, c(19L, 20L), workers = 2L),
    sb_compare_repeated_cv(task, pp, 3L, c(19L, 20L)))
  bad <- task
  bad$data$gene_b <- bad$data$gene_a
  failed <- sb_compare_repeated_cv(bad, pp[1], 3L, c(19L, 20L))
  expect_true(any(failed$comparison$status != "success"))
  expect_identical(sb_compare_repeated_cv(bad, pp[1], 3L, c(19L, 20L), workers = 2L),
    failed)
})

test_that("parallel repeated CV preserves absent caller RNG", {
  d <- sb_example_binary(48L, 1101L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  withr::local_preserve_seed()
  if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  result <- sb_compare_repeated_cv(task, list(sb_pipeline("none", "logistic")),
    3L, c(19L, 20L), workers = 2L)
  expect_identical(result$n_repeats, 2L)
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
})

test_that("worker count is bounded and worker errors keep their class", {
  d <- sb_example_binary(48L, 12L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  pp <- list(sb_pipeline("none", "logistic"))
  for (workers in list(0L, 3L, NA_real_, c(1L, 2L)))
    expect_error(sb_compare_repeated_cv(task, pp, workers = workers), class = "stormbreaker_bad_workers")
  connections <- nrow(showConnections(all = TRUE))
  expect_error(sb_compare_repeated_cv(task, list("invalid"), seeds = 1:2, workers = 2L),
    class = "stormbreaker_bad_pipeline_collection")
  expect_identical(nrow(showConnections(all = TRUE)), connections)
  expect_error(sb_compare_repeated_cv(task, list(sb_pipeline("none", "adaboost_binary")),
    seeds = 1:2, workers = 2L), class = "stormbreaker_parallel_sequential_rng")
})

test_that("parallel stack repetitions keep their complete results", {
  skip_if_not_installed("glmnet")
  d <- sb_example_binary(72L, 1101L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  stack <- sb_stack(list(sb_pipeline("none", "logistic"),
    sb_pipeline("none", "gaussian_naive_bayes")), inner_folds = 3L, inner_seed = 41L)
  expect_identical(sb_compare_repeated_cv(task, stack, 3L, c(19L, 20L), workers = 2L),
    sb_compare_repeated_cv(task, stack, 3L, c(19L, 20L)))
  d <- sb_example_survival(72L, 11032L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"),
    "month", prediction_horizon = 7)
  stack <- sb_stack(list(sb_pipeline("none", "coxph"),
    sb_pipeline("univ_cox_top_k", "coxph", selector_n_features = 2L)),
    inner_folds = 3L, inner_seed = 41L)
  expect_identical(sb_compare_repeated_cv(task, stack, 3L, c(19L, 20L), workers = 2L),
    sb_compare_repeated_cv(task, stack, 3L, c(19L, 20L)))
})
})

# ---- cohort-comparison ----
local({
# ---- internal-external-cv ----

test_that("IECV preserves group roles and rejects patients shared by cohorts", {
  tasks <- lapply(c("a", "b"), function(cohort) {
    d <- sb_example_binary(40L, 3754L)
    d$patient_id <- paste0(cohort, d$patient_id)
    group_col <- paste0(cohort, "_group")
    d[[group_col]] <- paste0(cohort, rep(1:20, each = 2))
    sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"),
                   group_col = group_col)
  })
  names(tasks) <- c("a", "b")
  validated <- sb_validate_iecv_cohorts(tasks)
  pooled <- sb_iecv_development_task(validated$tasks, tasks[[1]], "binary")
  expect_identical(pooled$group_col, "a_group")
  expect_identical(pooled$data$a_group,
    c(tasks$a$data$a_group, tasks$b$data$b_group))
  tasks$b$data$b_group[1:2] <- tasks$a$data$a_group[1]
  error <- capture_sb_error(sb_validate_iecv_cohorts(tasks))
  expect_s3_class(error, "stormbreaker_internal_external_patient_overlap")
})

make_sb_iecv_binary_task <- function(cohort_id, seed, predictors,
                                     id_col = "patient_id",
                                     outcome_col = "outcome") {
  data <- make_sb_binary_fixture(n = 36L, seed = seed)
  data$patient_id <- paste0(cohort_id, "_", data$patient_id)
  names(data)[names(data) == "patient_id"] <- id_col
  names(data)[names(data) == "outcome"] <- outcome_col
  sb_binary_task(
    data,
    id_col,
    outcome_col,
    predictors,
    decision_thresholds = c(0.2, 0.4)
  )
}

make_sb_iecv_survival_task <- function(cohort_id, seed, predictors,
                                       horizon = 6) {
  data <- make_sb_fixture(n = 36L, seed = seed)
  data$patient_id <- paste0(cohort_id, "_", data$patient_id)
  sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = horizon,
    decision_thresholds = c(0.2, 0.4)
  )
}

test_that("binary IECV holds out each cohort from fitting", {
  cohorts <- list(
    cohort_a = make_sb_iecv_binary_task(
      "a", 9101L, c("gene_a", "gene_b", "gene_c")
    ),
    cohort_b = make_sb_iecv_binary_task(
      "b", 9102L, c("gene_c", "gene_a", "gene_b"),
      id_col = "subject_id", outcome_col = "label"
    ),
    cohort_c = make_sb_iecv_binary_task(
      "c", 9103L, c("gene_b", "gene_c", "gene_a")
    )
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )

  observed <- sb_run_internal_external_cv(cohorts, pipelines)
  pipeline_ids <- vapply(pipelines, `[[`, character(1L), "id")

  expect_s3_class(observed, "stormbreaker_internal_external_cv")
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(
    observed$validation_scheme,
    "leave_one_cohort_out_internal_external_cv"
  )
  expect_identical(observed$pipeline_ids, pipeline_ids)
  expect_identical(names(observed$folds), names(cohorts))
  expect_identical(nrow(observed$results), 6L)
  expect_identical(
    observed$results$heldout_cohort_id,
    rep(names(cohorts), each = 2L)
  )
  expect_identical(observed$results$pipeline_id, rep(pipeline_ids, 3L))
  expect_identical(
    observed$results$n_development_patients,
    rep(72L, 6L)
  )

  all_ids <- lapply(cohorts, function(task) {
    as.character(task$data[[task$id_col]])
  })
  for (heldout_id in names(cohorts)) {
    fold <- observed$folds[[heldout_id]]
    development_ids <- unlist(
      all_ids[fold$development_cohort_ids],
      use.names = FALSE
    )
    heldout_ids <- all_ids[[heldout_id]]

    expect_identical(fold$heldout_cohort_id, heldout_id)
    expect_setequal(fold$development_patient_ids, development_ids)
    expect_length(intersect(fold$development_patient_ids, heldout_ids), 0L)
    expect_identical(fold$heldout_patient_ids, heldout_ids)
    expect_identical(fold$comparison$status, "success")
    expect_identical(names(fold$fits), pipeline_ids)
    expect_identical(names(fold$predictions), pipeline_ids)
    for (pipeline_id in pipeline_ids) {
      expect_identical(
        fold$fits[[pipeline_id]]$development$patient_ids,
        fold$development_patient_ids
      )
      expect_identical(
        fold$predictions[[pipeline_id]]$predictions$patient_id,
        heldout_ids
      )
    }
  }

  expect_false(any(c(
    "pooled_metric", "ranking", "winner", "final_fit"
  ) %in% names(observed)))
  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "different fitted model", fixed = TRUE)
  expect_match(
    printed,
    "not final independent external validation",
    fixed = TRUE
  )
})

test_that("held-out binary outcomes cannot change their fold predictions", {
  cohorts <- list(
    cohort_a = make_sb_iecv_binary_task(
      "a", 9201L, c("gene_a", "gene_b")
    ),
    cohort_b = make_sb_iecv_binary_task(
      "b", 9202L, c("gene_b", "gene_a")
    ),
    cohort_c = make_sb_iecv_binary_task(
      "c", 9203L, c("gene_a", "gene_b")
    )
  )
  pipeline <- sb_pipeline(learner = "logistic")
  original <- sb_run_internal_external_cv(cohorts, list(pipeline))

  changed_data <- cohorts$cohort_c$data
  changed_data$outcome <- 1L - changed_data$outcome
  changed <- cohorts
  changed$cohort_c <- sb_binary_task(
    changed_data,
    "patient_id",
    "outcome",
    cohorts$cohort_c$predictors,
    decision_thresholds = cohorts$cohort_c$decision_thresholds
  )
  altered <- sb_run_internal_external_cv(changed, list(pipeline))

  original_fold <- original$folds$cohort_c
  altered_fold <- altered$folds$cohort_c
  pipeline_id <- pipeline$id
  expect_identical(
    original_fold$development_patient_ids,
    altered_fold$development_patient_ids
  )
  expect_identical(
    original_fold$fits[[pipeline_id]]$learner_state,
    altered_fold$fits[[pipeline_id]]$learner_state
  )
  expect_equal(
    original_fold$predictions[[pipeline_id]]$predictions$probability,
    altered_fold$predictions[[pipeline_id]]$predictions$probability,
    tolerance = 0
  )
})

test_that("binary stack IECV fits only the non-held-out cohorts", {
  skip_if_not_installed("glmnet")

  cohorts <- list(
    cohort_a = make_sb_iecv_binary_task(
      "a", 9251L, c("gene_a", "gene_b")
    ),
    cohort_b = make_sb_iecv_binary_task(
      "b", 9252L, c("gene_b", "gene_a")
    ),
    cohort_c = make_sb_iecv_binary_task(
      "c", 9253L, c("gene_a", "gene_b")
    )
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 9254L)

  observed <- sb_run_internal_external_cv(cohorts, stack)

  expect_identical(observed$status, "success")
  expect_identical(observed$strategy_type, "stack")
  expect_identical(observed$pipeline_ids, stack$id)
  expect_identical(observed$results$pipeline_id, rep(stack$id, 3L))
  for (heldout_id in names(cohorts)) {
    fold <- observed$folds[[heldout_id]]
    heldout_ids <- as.character(
      cohorts[[heldout_id]]$data[[cohorts[[heldout_id]]$id_col]]
    )
    development_ids <- unlist(lapply(
      cohorts[names(cohorts) != heldout_id],
      function(task) as.character(task$data[[task$id_col]])
    ), use.names = FALSE)
    fit <- fold$fits[[stack$id]]
    prediction <- fold$predictions[[stack$id]]

    expect_identical(names(fold$fits), stack$id)
    expect_identical(names(fold$predictions), stack$id)
    expect_setequal(fit$development$patient_ids, development_ids)
    expect_length(intersect(fit$development$patient_ids, heldout_ids), 0L)
    expect_identical(names(fit$member_fits), stack$pipeline_ids)
    expect_identical(
      prediction$predictions$patient_id,
      heldout_ids
    )
    expect_identical(fold$comparison$comparison$pipeline_id, stack$id)
  }

  changed_data <- cohorts$cohort_c$data
  changed_data$outcome <- 1L - changed_data$outcome
  changed <- cohorts
  changed$cohort_c <- sb_binary_task(
    changed_data,
    cohorts$cohort_c$id_col,
    cohorts$cohort_c$outcome_col,
    cohorts$cohort_c$predictors,
    decision_thresholds = cohorts$cohort_c$decision_thresholds
  )
  altered <- sb_run_internal_external_cv(changed, stack)
  original_fold <- observed$folds$cohort_c
  altered_fold <- altered$folds$cohort_c

  expect_identical(
    original_fold$fits[[stack$id]]$meta$state,
    altered_fold$fits[[stack$id]]$meta$state
  )
  expect_identical(
    lapply(
      original_fold$fits[[stack$id]]$member_fits,
      `[[`,
      "learner_state"
    ),
    lapply(
      altered_fold$fits[[stack$id]]$member_fits,
      `[[`,
      "learner_state"
    )
  )
  expect_equal(
    original_fold$predictions[[stack$id]]$predictions$probability,
    altered_fold$predictions[[stack$id]]$predictions$probability,
    tolerance = 0
  )
})

test_that("survival IECV retains one declared time target", {
  cohorts <- list(
    cohort_a = make_sb_iecv_survival_task(
      "a", 9301L, c("gene_a", "gene_b", "gene_c")
    ),
    cohort_b = make_sb_iecv_survival_task(
      "b", 9302L, c("gene_c", "gene_b", "gene_a")
    ),
    cohort_c = make_sb_iecv_survival_task(
      "c", 9303L, c("gene_b", "gene_a", "gene_c")
    )
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "coxph",
    selector_n_features = 2L
  )

  observed <- sb_run_internal_external_cv(cohorts, list(pipeline))

  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "survival")
  expect_identical(observed$cohorts$time_unit, rep("month", 3L))
  expect_identical(observed$cohorts$prediction_horizon, rep(6, 3L))
  expect_identical(
    c("ipcw_brier", "time_auc", "uno_c") %in% names(observed$results),
    rep(TRUE, 3L)
  )
  for (fold in observed$folds) {
    expect_identical(fold$comparison$prediction_target$has_horizon, TRUE)
    expect_identical(
      fold$comparison$prediction_target$prediction_horizon,
      6
    )
    expect_identical(
      fold$comparison$prediction_target$time_unit,
      "month"
    )
    expect_identical(
      fold$predictions[[pipeline$id]]$decision_thresholds,
      c(0.2, 0.4)
    )
  }
})

test_that("survival stack IECV keeps its baseline inside development cohorts", {
  skip_if_not_installed("glmnet")

  cohorts <- list(
    cohort_a = make_sb_iecv_survival_task(
      "sa", 9321L, c("gene_a", "gene_b", "gene_c")
    ),
    cohort_b = make_sb_iecv_survival_task(
      "sb", 9322L, c("gene_c", "gene_a", "gene_b")
    ),
    cohort_c = make_sb_iecv_survival_task(
      "sc", 9323L, c("gene_b", "gene_c", "gene_a")
    )
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "coxph"),
    sb_pipeline(learner = "survreg_weibull")
  ), inner_folds = 3L, inner_seed = 9324L)

  observed <- sb_run_internal_external_cv(cohorts, stack)

  expect_identical(observed$status, "success")
  expect_identical(observed$strategy_type, "stack")
  for (heldout_id in names(cohorts)) {
    fold <- observed$folds[[heldout_id]]
    fit <- fold$fits[[stack$id]]
    prediction <- fold$predictions[[stack$id]]
    heldout_ids <- as.character(
      cohorts[[heldout_id]]$data[[cohorts[[heldout_id]]$id_col]]
    )

    expect_length(intersect(fit$development$patient_ids, heldout_ids), 0L)
    expect_identical(
      fit$meta$state$baseline_hazard$prediction_horizon,
      6
    )
    expect_identical(
      prediction$predictions$prediction_horizon,
      rep(6, length(heldout_ids))
    )
    expect_false(any(!is.finite(prediction$predictions$risk)))
    expect_false(any(!is.finite(
      prediction$predictions$survival_probability
    )))
  }
})

test_that("risk-only survival IECV does not invent a horizon", {
  make_risk_task <- function(cohort_id, seed) {
    data <- make_sb_fixture(n = 30L, seed = seed)
    data$patient_id <- paste0(cohort_id, "_", data$patient_id)
    sb_surv_task(
      data,
      "patient_id",
      "time",
      "event",
      c("gene_a", "gene_b"),
      "month"
    )
  }
  cohorts <- list(
    cohort_a = make_risk_task("a", 9351L),
    cohort_b = make_risk_task("b", 9352L)
  )

  observed <- sb_run_internal_external_cv(
    cohorts,
    list(sb_pipeline())
  )

  expect_identical(observed$status, "success")
  expect_identical(
    "prediction_horizon" %in% names(observed$cohorts),
    FALSE
  )
  expect_identical(
    c("ipcw_brier", "time_auc", "uno_c") %in% names(observed$results),
    rep(FALSE, 3L)
  )
  for (fold in observed$folds) {
    expect_identical(fold$comparison$prediction_target$has_horizon, FALSE)
    expect_identical(
      fold$comparison$prediction_target$prediction_horizon,
      NULL
    )
  }
})

test_that("IECV validates cohort isolation and common targets", {
  binary_a <- make_sb_iecv_binary_task(
    "a", 9401L, c("gene_a", "gene_b")
  )
  binary_b <- make_sb_iecv_binary_task(
    "b", 9402L, c("gene_a", "gene_b")
  )
  different_features <- make_sb_iecv_binary_task(
    "c", 9403L, c("gene_a", "gene_c")
  )
  different_thresholds <- binary_b
  different_thresholds$decision_thresholds <- c(0.3, 0.5)
  different_thresholds <- sb_binary_task(
    different_thresholds$data,
    different_thresholds$id_col,
    different_thresholds$outcome_col,
    different_thresholds$predictors,
    different_thresholds$decision_thresholds
  )
  overlapping <- binary_b
  overlapping$data[[overlapping$id_col]] <- binary_a$data[[binary_a$id_col]]
  overlapping <- sb_binary_task(
    overlapping$data,
    overlapping$id_col,
    overlapping$outcome_col,
    overlapping$predictors,
    overlapping$decision_thresholds
  )
  survival <- make_sb_iecv_survival_task(
    "s", 9404L, c("gene_a", "gene_b")
  )
  pipeline <- sb_pipeline(learner = "logistic")

  expect_s3_class(
    capture_sb_error(sb_run_internal_external_cv(
      list(binary_a), list(pipeline)
    )),
    "stormbreaker_bad_internal_external_cohorts"
  )
  expect_s3_class(
    capture_sb_error(sb_run_internal_external_cv(
      list(a = binary_a, s = survival), list(pipeline)
    )),
    "stormbreaker_internal_external_task_mismatch"
  )
  expect_s3_class(
    capture_sb_error(sb_run_internal_external_cv(
      list(a = binary_a, c = different_features), list(pipeline)
    )),
    "stormbreaker_internal_external_predictor_mismatch"
  )
  expect_s3_class(
    capture_sb_error(sb_run_internal_external_cv(
      list(a = binary_a, b = different_thresholds), list(pipeline)
    )),
    "stormbreaker_internal_external_target_mismatch"
  )
  expect_s3_class(
    capture_sb_error(sb_run_internal_external_cv(
      list(a = binary_a, b = overlapping), list(pipeline)
    )),
    "stormbreaker_internal_external_patient_overlap"
  )
  expect_s3_class(
    capture_sb_error(sb_run_internal_external_cv(
      list(a = binary_a, b = binary_b), list(pipeline, pipeline)
    )),
    "stormbreaker_duplicate_pipeline_id"
  )
})

test_that("IECV continues later pipelines after each earlier fold-fit failure", {
  make_high_dimensional_task <- function(cohort_id, seed) {
    data <- make_sb_fixture(n = 14L, seed = seed)
    data$patient_id <- paste0(cohort_id, "_", data$patient_id)
    withr::with_seed(seed + 1000L, {
      for (index in 4:31) {
        data[[sprintf("gene_%02d", index)]] <- stats::rnorm(nrow(data))
      }
    })
    predictors <- c(
      "gene_a", "gene_b", "gene_c", sprintf("gene_%02d", 4:31)
    )
    sb_surv_task(
      data, "patient_id", "time", "event", predictors, "month"
    )
  }
  cohorts <- list(
    cohort_a = make_high_dimensional_task("a", 9501L),
    cohort_b = make_high_dimensional_task("b", 9502L),
    cohort_c = make_high_dimensional_task("c", 9503L)
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      learner = "coxph",
      selector_n_features = 2L
    )
  )

  observed <- sb_run_internal_external_cv(cohorts, pipelines)
  pipeline_ids <- vapply(pipelines, `[[`, character(1L), "id")

  expect_identical(observed$status, "incomplete")
  expect_identical(
    observed$results$status,
    rep(c("failed", "success"), 3L)
  )
  for (fold in observed$folds) {
    expect_identical(fold$fits[[pipeline_ids[[1L]]]]$status, "failed")
    expect_identical(
      fold$predictions[[pipeline_ids[[1L]]]]$status,
      "failed"
    )
    expect_identical(fold$fits[[pipeline_ids[[2L]]]]$status, "success")
    expect_identical(
      fold$predictions[[pipeline_ids[[2L]]]]$status,
      "success"
    )
  }
})

test_that("IECV keeps a failed stack intact without member fallback", {
  data_a <- make_sb_iecv_binary_task(
    "fa", 9551L, c("gene_a", "gene_b")
  )
  data_b <- make_sb_iecv_binary_task(
    "fb", 9552L, c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 9553L)
  local_mocked_bindings(sb_glmnet_available = function() FALSE)

  observed <- sb_run_internal_external_cv(
    list(cohort_a = data_a, cohort_b = data_b),
    stack
  )

  expect_identical(observed$status, "incomplete")
  expect_identical(observed$results$status, rep("failed", 2L))
  for (fold in observed$folds) {
    expect_identical(fold$fits[[stack$id]]$status, "failed")
    expect_identical(fold$fits[[stack$id]]$error_stage, "meta_fit")
    expect_length(fold$fits[[stack$id]]$member_fits, 0L)
    expect_identical(fold$predictions[[stack$id]]$status, "failed")
    expect_identical(
      fold$predictions[[stack$id]]$error_stage,
      "fit"
    )
  }
})

# ---- multicohort-result-page ----

page_tasks <- function() {
  make <- function(prefix, seed, n) {
    d <- sb_example_binary(n = n, seed = seed)
    d$patient_id <- paste0(prefix, "_", d$patient_id)
    sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  }
  list(development = make("dev", 312, 48),
       validation = make("val", 313, 30),
       test = make("test", 314, 30))
}

test_that("cohort comparison reuses the development fits and saved metrics", {
  tasks <- page_tasks()
  p <- list(sb_pipeline(learner = "logistic"),
            sb_pipeline(learner = "gaussian_naive_bayes"))
  x <- sb_compare_cohorts(tasks$development, tasks[c("validation", "test")],
    roles = c(validation = "validation", test = "test"), pipelines = p,
    folds = 3, seed = 52, view = FALSE)
  direct <- sb_compare_cv(tasks$development, p, folds = 3, seed = 52)
  expect_equal(x$development$comparison, direct$comparison, tolerance = 0)
  expect_identical(x$cohort_roles, c(Development = "development_cv",
    validation = "validation", test = "test"))
  for (i in seq_along(p)) {
    fit <- sb_fit_final(tasks$development, p[[i]])
    expect_equal(x$fits[[i]], fit, tolerance = 0)
    pred <- sb_predict_final(fit, tasks$test$data, "patient_id")
    expect_equal(x$predictions$test[[i]], pred, tolerance = 0)
  }
  page <- sb_result_page_data(x)
  expect_identical(page$default_sort, "Development")
  expect_identical(vapply(page$datasets, `[[`, character(1), "role"),
    c("development_cv", "validation", "test"))
  expect_length(page$models, 2)
  auc <- Filter(function(cell) cell$metric == "roc_auc" && cell$dataset == "test", page$cells)
  expect_equal(vapply(auc, `[[`, numeric(1), "estimate"),
    x$cohorts$test$comparison$roc_auc, tolerance = 0)
})

test_that("held-out outcomes cannot alter the development models or predictions", {
  tasks <- page_tasks()
  args <- list(development = tasks$development, cohorts = list(test = tasks$test),
    roles = c(test = "test"), models = "logistic", folds = 3, view = FALSE)
  a <- do.call(sb_compare_cohorts, args)
  changed <- tasks$test$data
  changed$outcome <- 1L - changed$outcome
  args$cohorts$test <- sb_binary_task(changed, "patient_id", "outcome", c("gene_a", "gene_b"))
  b <- do.call(sb_compare_cohorts, args)
  expect_equal(a$fits, b$fits, tolerance = 0)
  expect_equal(a$predictions, b$predictions, tolerance = 0)
  expect_equal(a$development, b$development, tolerance = 0)
})

test_that("friendly labels describe algorithms without changing configuration identity", {
  forest <- sb_display_learner("ranger_binary")
  expect_match(forest, "ranger", fixed = TRUE)
  expect_false(grepl("fixed|_binary", forest))
  expect_identical(sb_display_selector("none"), "None")
  expect_identical(sb_display_learner("glmnet_logistic", list(alpha = 0)), "Ridge logistic")
  expect_identical(sb_display_learner("glmnet_logistic", list(alpha = 1)), "Lasso logistic")
  expect_false(identical(sb_display_learner("xgboost_binary"),
    sb_display_learner("xgboost_binary_dart")))
  catalog <- sb_model_catalog("binary")
  expect_identical(catalog$model_label[catalog$model_id == "ridge_logistic"], "Ridge logistic")
  expect_identical(catalog$model_label[catalog$model_id == "glmnet_logistic"], "Lasso logistic")
  expect_false(any(grepl("fixed|_binary|_survival", catalog$model_label)))
  expect_false(any(grepl("fixed|_binary|_survival", sb_model_catalog("survival")$model_label)))
})

test_that("result pages contain summary data with escaped labels and no patient rows", {
  tasks <- page_tasks()
  cv <- sb_compare_models(tasks$development, models = "logistic", folds = 3, view = FALSE)
  label <- "</script><img src=x onerror=alert(1)>"
  named <- stats::setNames(list(cv), label)
  path <- tempfile(fileext = ".html")
  expect_identical(sb_view_results(named, path = path, open = FALSE), normalizePath(path, winslash = "/"))
  html <- paste(readLines(path, encoding = "UTF-8"), collapse = "\n")
  expect_false(grepl(label, html, fixed = TRUE))
  expect_false(grepl(tasks$development$data$patient_id[1], html, fixed = TRUE))
  expect_match(html, "Logistic regression", fixed = TRUE)
  expect_error(sb_view_results(named, path = path, open = FALSE),
    class = "stormbreaker_result_page_exists")
  expect_identical(sb_view_results(named, path = path, open = FALSE, overwrite = TRUE),
    normalizePath(path, winslash = "/"))
})

test_that("report opening is optional and never discards a completed comparison", {
  tasks <- page_tasks()
  opened <- character()
  local_mocked_bindings(sb_open_result_page = function(path) opened <<- c(opened, path))
  x <- sb_compare_models(tasks$development, models = "logistic", folds = 3, view = TRUE)
  expect_length(opened, 1)
  expect_identical(x$report_file, opened[[1]])
  expect_true(file.exists(x$report_file))
  opened <- character()
  y <- sb_compare_models(tasks$development, models = "logistic", folds = 3, view = FALSE)
  expect_length(opened, 0)
  expect_null(y$report_file)
  expect_equal(x$comparison, y$comparison, tolerance = 0)
})

test_that("roles and cross-cohort identity checks are explicit", {
  tasks <- page_tasks()
  expect_error(sb_compare_cohorts(tasks$development, list(test = tasks$test),
    roles = c(test = "typo"), models = "logistic", view = FALSE),
    class = "stormbreaker_bad_cohort_roles")
  expect_error(sb_compare_cohorts(tasks$development, list(test = tasks$development),
    roles = c(test = "test"), models = "logistic", view = FALSE),
    class = "stormbreaker_internal_external_patient_overlap")
})

test_that("survival cells preserve horizon-specific estimates and unavailable values", {
  make <- function(prefix, seed) {
    d <- sb_example_survival(45, seed)
    d$patient_id <- paste0(prefix, d$patient_id)
    sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"),
      time_unit = "month", prediction_horizon = 6)
  }
  dev <- make("d", 88)
  ext <- make("e", 89)
  x <- sb_compare_cohorts(dev, list(Test = ext), roles = c(Test = "test"),
    models = c("coxph", "survreg_weibull"), folds = 3, view = FALSE)
  page <- sb_result_page_data(x)
  cells <- Filter(function(cell) cell$metric == "time_auc" && cell$dataset == "Test", page$cells)
  expect_equal(vapply(cells, function(cell) cell$estimate, numeric(1)),
    x$cohorts$Test$comparison$time_auc, tolerance = 0)
  expect_identical(vapply(cells, function(cell) cell$horizon, numeric(1)), c(6, 6))
  expect_identical(vapply(cells, function(cell) cell$time_unit, character(1)), c("month", "month"))
  x$cohorts$Test$evaluations[[1]]$metrics$time_auc$estimate <- Inf
  x$results$Test <- x$cohorts$Test
  missing <- Filter(function(cell) cell$metric == "time_auc" && cell$dataset == "Test",
    sb_result_page_data(x)$cells)
  expect_false(missing[[1]]$available)
  expect_match(sb_page_json(missing[[1]]), '"estimate":null', fixed = TRUE)
})

test_that("JSON encoding preserves strings and blocks script termination", {
  expect_identical(sb_page_json(list(text = '"\\\n</script>&')),
    '{"text":"\\"\\\\\\u000a\\u003c/script\\u003e\\u0026"}')
  expect_identical(sb_page_json(list(value = NA_real_, values = list(1, TRUE, NULL))),
    '{"value":null,"values":[1,true,null]}')
})

test_that("opening failure retains the saved page and completed results", {
  tasks <- page_tasks()
  local_mocked_bindings(sb_open_result_page = function(path) stop("browser unavailable"))
  expect_warning(
    x <- sb_compare_models(tasks$development, models = "logistic", folds = 3, view = TRUE),
    "page was saved"
  )
  expect_true(file.exists(x$report_file))
  expect_identical(x$comparison$status, "success")
})

test_that("a saved failed run remains a visible row", {
  tasks <- page_tasks()
  cv <- sb_compare_models(tasks$development, models = "logistic", folds = 3, view = FALSE)
  cv$results[[1]]$status <- "failed"
  page <- sb_result_page_data(cv)
  expect_length(page$models, 1)
  expect_false(any(vapply(page$cells, function(cell) cell$available, logical(1))))
  expect_identical(page$cells[[1]]$run_status, "failed")
})
})

test_that("binary holdout preserves outcome strata, row IDs and RNG state", {
  d <- sb_example_binary(20L, 72L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"),
    decision_thresholds = c(.25, .5))
  withr::local_seed(91L)
  state <- .Random.seed
  first <- sb_split_binary(task, seed = 6L)
  expect_identical(first, sb_split_binary(task, .7, 6L))
  expect_identical(.Random.seed, state)
  expect_equal(first$counts$n_rows, c(14, 6))
  expect_equal(first$counts$n_positive, c(7, 3))
  expect_equal(first$counts$n_negative, c(7, 3))
  expect_identical(first$assignment$observation_id, d$patient_id)
  expect_length(intersect(first$development$data$patient_id, first$assessment$patient_id), 0L)
  expect_identical(first$development$decision_thresholds, c(.25, .5))
  expect_identical(attr(first$assessment, "stormbreaker_assessment_role"), "internal_holdout")
  expect_error(sb_split_binary(task, train_prop = 1), class = "stormbreaker_bad_split_proportion")
})

test_that("a mixed-outcome patient stays intact during the split", {
  d <- data.frame(obs = letters[1:4], patient = c("mixed", "mixed", "negative", "positive"),
    y = c(0L, 1L, 0L, 1L), x = c(2, 4, 6, 8))
  task <- sb_binary_task(d, "obs", "y", "x", group_col = "patient")
  split <- sb_split_binary(task, .7)
  expect_length(intersect(split$development$data$patient, split$assessment$patient), 0L)
  expect_equal(split$counts$n_rows, c(2, 2))
  expect_equal(split$counts$n_positive, c(1, 1))
  expect_equal(split$counts$n_negative, c(1, 1))
  impossible <- sb_binary_task(d[1:3, ], "obs", "y", "x", group_col = "patient")
  expect_error(sb_split_binary(impossible), class = "stormbreaker_insufficient_split_groups")
})
