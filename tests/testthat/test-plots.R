# ---- plot-calibration ----
local({
# ---- clinical-thresholds ----

clinical_threshold_fixture <- function() {
  structure(list(status = "success", task_type = "binary", pipeline_id = "clinical",
    scored_predictions = data.frame(patient_id = letters[1:4], outcome = c(1L, 0L, 1L, 0L),
      probability = c(.9, .7, .7, .1), status = "success"), metrics = list()),
    class = "stormbreaker_external_evaluation")
}

test_that("Kappa uses observed margins and preserves negative agreement", {
  x <- clinical_threshold_fixture()
  expect_equal(sb_threshold_metrics(x, c(0, .7, 1))$table$kappa, c(0, .5, 0))
  x$scored_predictions$probability <- x$scored_predictions$outcome
  expect_equal(sb_threshold_metrics(x, .5)$table$kappa, 1)
  x$scored_predictions$probability <- 1 - x$scored_predictions$outcome
  r <- sb_threshold_metrics(x, .5)
  expect_equal(r$table$kappa, -1)
  p <- sb_plot_thresholds(r, metrics = "kappa")
  expect_equal(ggplot2::ggplot_build(p)$data[[2]]$y, -1)
  x$scored_predictions$outcome[] <- 1L
  x$scored_predictions$probability[] <- 1
  expect_equal(sb_threshold_metrics(x, .5)$table$kappa, NA_real_)
})

test_that("clinical thresholds preserve ties and all confusion counts", {
  x <- clinical_threshold_fixture()
  before <- x
  result <- sb_threshold_metrics(x, c(.7, 1, 0, .7))
  expect_equal(result$table$threshold, c(0, .7, 1))
  expect_equal(result$table$tp, c(2, 2, 0))
  expect_equal(result$table$fp, c(2, 1, 0))
  expect_equal(result$table$tn, c(0, 1, 2))
  expect_equal(result$table$fn, c(0, 0, 2))
  expect_equal(result$table$ppv, c(.5, 2 / 3, NA_real_))
  expect_equal(result$table$npv, c(NA_real_, 1, .5))
  expect_identical(x, before)
  expect_length(intersect("patient_id", names(result$table)), 0)
})

test_that("single classes retain defined metrics and invalid rows are not dropped", {
  x <- clinical_threshold_fixture()
  x$scored_predictions$outcome[] <- 1L
  r <- sb_threshold_metrics(x, .5)
  expect_equal(r$table$sensitivity, .75)
  expect_equal(r$table$specificity, NA_real_)
  x$scored_predictions$probability[1] <- NA_real_
  r <- sb_threshold_metrics(x, .5)
  expect_equal(nrow(r$table), 0L)
  expect_identical(r$diagnostics$status, "invalid_predictions")
})

test_that("threshold plotting uses all requested points", {
  skip_if_not_installed("ggplot2")
  p <- sb_plot_thresholds(sb_threshold_metrics(clinical_threshold_fixture(), c(.2, .7, .95)))
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 12L)
  expect_length(unique(p$data$metric), 4L)
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 1L)
})

test_that("confusion cells preserve class orientation and counts", {
  skip_if_not_installed("ggplot2")
  x <- sb_threshold_metrics(clinical_threshold_fixture(), .7)
  p <- sb_plot_confusion(x, .7)
  expect_equal(p$data$count, c(1, 1, 0, 2))
  expect_equal(as.character(p$data$observed), c("0", "0", "1", "1"))
  q <- sb_plot_confusion(x, .7, normalize = "observed")
  expect_equal(q$data$value, c(.5, .5, 0, 1))
  expect_s3_class(ggplot2::ggplot_build(q), "ggplot_built")
})

test_that("different models never share confusion-cell text positions", {
  skip_if_not_installed("ggplot2")
  a <- clinical_threshold_fixture()
  b <- a
  b$pipeline_id <- "second_model"
  b$scored_predictions$probability <- 1 - b$scored_predictions$probability
  input <- structure(list(evaluations = list(a, b)), class = "stormbreaker_external_comparison")
  thresholds <- sb_threshold_metrics(input, .5)
  p <- sb_plot_confusion(thresholds, .5)
  built <- ggplot2::ggplot_build(p)
  labels <- built$data[[2L]]
  expect_equal(nrow(built$layout$layout), 2L)
  expect_equal(nrow(labels), 8L)
  expect_equal(anyDuplicated(labels[c("PANEL", "x", "y")]), 0L)
  expect_equal(sum(p$data$count), 8L)
})

# ---- plot-calibration ----

test_that("binary and survival calibration coordinates are copied exactly", {
  skip_if_not_installed("ggplot2")
  for (type in c("binary", "survival")) {
    fixture <- make_sb_curve_view_fixture(type)
    for (boundary in c("cv", "external")) {
      x <- fixture[[boundary]]
      before <- x
      p <- sb_plot_calibration(x, cohort_label = "Demonstration cohort")
      key <- if (type == "binary") "moderate_calibration" else "calibration"
      metric <- x$metrics[[key]]
      x_column <- if (type == "binary") "probability" else "predicted_event_probability"
      y_column <- if (type == "binary") "calibrated_probability" else "calibrated_event_probability"
      expect_s3_class(p, "ggplot")
      expect_null(p$labels$caption)
      expect_null(p$labels$subtitle)
      expect_identical(p$data$x, metric$curve[[x_column]])
      expect_identical(p$data$y, metric$curve[[y_column]])
      expect_identical(p$data$source_row, seq_len(nrow(metric$curve)))
      info <- attr(p, "stormbreaker_curve_status")
      expect_identical(info$scope, metric$scope)
      expect_identical(info$metric_status, metric$status)
      expect_identical(info$run_status, x$status)
      expect_identical(info$cohort, "Demonstration cohort")
      expect_identical(info$validation,
                       if (boundary == "cv") "development_cv" else "external_evaluation")
      if (type == "survival") {
        expect_identical(info$prediction_horizon, 6)
        expect_identical(info$time_unit, "month")
        expect_match(as.character(p$data$panel[[1]]), "6 month", fixed = TRUE)
        expect_match(p$labels$x, "event probability", fixed = TRUE)
      }
      built <- ggplot2::ggplot_build(p)
      expect_equal(built$data[[2]]$x, p$data$x)
      expect_equal(built$data[[2]]$y, p$data$y)
      expect_identical(x, before)
    }
  }
})

test_that("both comparison classes preserve panels and all original curves", {
  skip_if_not_installed("ggplot2")
  for (type in c("binary", "survival")) {
    f <- make_sb_curve_view_fixture(type, comparison = TRUE)
    for (boundary in c("cv", "external")) {
      x <- f[[boundary]]
      results <- if (boundary == "cv") x$results else x$evaluations
      ids <- vapply(results, function(z) {
        if (boundary == "cv") z$pipeline$id else z$pipeline_id
      }, character(1))
      p <- sb_plot_calibration(x)
      info <- attr(p, "stormbreaker_curve_status")
      expect_identical(info$pipeline_id, unname(ids))
      expect_identical(levels(p$data$panel), info$panel)
      expect_identical(nrow(ggplot2::ggplot_build(p)$layout$layout), 2L)
      key <- if (type == "binary") "moderate_calibration" else "calibration"
      column <- if (type == "binary") "probability" else "predicted_event_probability"
      for (i in seq_along(results)) {
        expect_identical(p$data$x[p$data$pipeline_id == ids[[i]]],
                         results[[i]]$metrics[[key]]$curve[[column]])
      }
    }
  }
})

test_that("explicit audit view keeps absent and failed calibration panels", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_survival(60, seed = 1)
  task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month")
  x <- sb_run_cv(task, folds = 3, seed = 1)
  p <- sb_plot_calibration(x, show_unavailable = TRUE)
  expect_identical(attr(p, "stormbreaker_curve_status")$metric_status, "not_requested")
  expect_identical(nrow(ggplot2::ggplot_build(p)$layout$layout), 1L)
  expect_identical(nrow(ggplot2::ggplot_build(p)$data[[2]]), 0L)
  expect_identical(p$data$source_row, NA_integer_)

  x <- make_sb_curve_view_fixture()$cv
  original <- x$metrics$moderate_calibration$curve$calibrated_probability
  x$status <- "failed"
  p <- sb_plot_calibration(x, show_unavailable = TRUE)
  expect_identical(p$data$y, original)
  expect_identical(attr(p, "stormbreaker_curve_status")$run_status, "failed")
  expect_identical(nrow(ggplot2::ggplot_build(p)$data[[2]]), 0L)
  expect_match(p$layers[[3]]$data$label, "No curve", fixed = TRUE)
})

test_that("calibration drawing never refits a model or calibration equation", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture()$external
  before <- x
  local_mocked_bindings(
    sb_fit_final = function(...) stop("No fitting"),
    sb_evaluate_external = function(...) stop("No evaluation"),
    sb_binary_moderate_calibration_metric = function(...) stop("No calibration"),
    sb_survival_calibration_metric = function(...) stop("No calibration")
  )
  expect_s3_class(sb_plot_calibration(x), "ggplot")
  expect_identical(x, before)
})

test_that("curve plotting rejects unsupported containers and invalid cohort labels", {
  skip_if_not_installed("ggplot2")
  expect_error(sb_plot_calibration(list()), class = "stormbreaker_bad_curve_view")
  expect_error(sb_plot_calibration(structure(list(results = list()),
                                            class = "stormbreaker_cv_comparison")),
               class = "stormbreaker_bad_curve_view")
  x <- make_sb_curve_view_fixture()$cv
  for (label in list(NA_character_, character(), c("a", "b"), "", 1)) {
    expect_error(sb_plot_calibration(x, cohort_label = label),
                 class = "stormbreaker_bad_cohort_label")
  }
})

# ---- plot-dca ----

test_that("DCA preserves every threshold and all three saved strategies", {
  skip_if_not_installed("ggplot2")
  for (type in c("binary", "survival")) {
    f <- make_sb_curve_view_fixture(type)
    for (boundary in c("cv", "external")) {
      x <- f[[boundary]]
      before <- x
      p <- sb_plot_dca(x, cohort_label = "Example cohort")
      curve <- x$metrics$decision_curve$curve
      for (strategy in c("model", "treat_all", "treat_none")) {
        rows <- p$data$strategy == strategy
        expect_identical(p$data$x[rows], curve$threshold_probability)
        expect_identical(p$data$y[rows], curve[[paste0(strategy, "_net_benefit")]])
        expect_identical(p$data$source_row[rows], seq_len(nrow(curve)))
      }
      info <- attr(p, "stormbreaker_curve_status")
      expect_identical(info$scope, x$metrics$decision_curve$scope)
      expect_identical(info$cohort, "Example cohort")
      expect_identical(info$metric_status, "success")
      if (type == "survival") {
        expect_identical(info$prediction_horizon, 6)
        expect_identical(info$time_unit, "month")
      }
      built <- ggplot2::ggplot_build(p)
      expect_equal(sort(built$data[[1]]$y), sort(p$data$y))
      expect_equal(sort(built$data[[2]]$x), sort(p$data$x))
      expect_identical(x, before)
    }
  }
})

test_that("DCA retains negative benefits and single-threshold points", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture()$cv
  x$metrics$decision_curve$curve$model_net_benefit[[1]] <- -0.4
  p <- sb_plot_dca(x)
  expect_identical(p$data$y[[1]], -0.4)
  expect_equal(min(ggplot2::ggplot_build(p)$data[[1]]$y), min(p$data$y))
  x$metrics$decision_curve$curve <- x$metrics$decision_curve$curve[1, , drop = FALSE]
  p <- sb_plot_dca(x)
  expect_identical(nrow(p$data), 3L)
  expect_identical(p$scales$get_scales("x")$limits, c(0, 1))
  expect_no_warning(ggplot2::ggplot_build(p))
  expect_identical(nrow(ggplot2::ggplot_build(p)$data[[2]]), 3L)
})

test_that("DCA comparisons keep reference strategies within each pipeline panel", {
  skip_if_not_installed("ggplot2")
  for (type in c("binary", "survival")) {
    f <- make_sb_curve_view_fixture(type, comparison = TRUE)
    for (boundary in c("cv", "external")) {
      x <- f[[boundary]]
      results <- if (boundary == "cv") x$results else x$evaluations
      p <- sb_plot_dca(x)
      info <- attr(p, "stormbreaker_curve_status")
      expect_identical(nrow(info), 2L)
      expect_identical(levels(p$data$panel), info$panel)
      expect_identical(nrow(ggplot2::ggplot_build(p)$layout$layout), 2L)
      for (i in seq_along(results)) {
        rows <- p$data$pipeline_id == info$pipeline_id[[i]] & p$data$strategy == "treat_all"
        expect_identical(p$data$y[rows],
                         results[[i]]$metrics$decision_curve$curve$treat_all_net_benefit)
      }
    }
  }
})

test_that("audit view distinguishes not-requested DCA from failed curves", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_binary(60, seed = 1)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  x <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1)
  p <- sb_plot_dca(x, show_unavailable = TRUE)
  expect_identical(attr(p, "stormbreaker_curve_status")$metric_status, "not_requested")
  expect_identical(p$scales$get_scales("x")$limits, c(0, 1))
  expect_no_warning(ggplot2::ggplot_build(p))
  expect_identical(nrow(ggplot2::ggplot_build(p)$data[[2]]), 0L)
  x$status <- "failed"
  p <- sb_plot_dca(x, show_unavailable = TRUE)
  expect_identical(attr(p, "stormbreaker_curve_status")$metric_status, "not_saved")

  x <- make_sb_curve_view_fixture("survival")$external
  x$metrics$decision_curve$status <- "incomplete"
  x$metrics$decision_curve$curve$model_net_benefit[[2]] <- NA_real_
  p <- sb_plot_dca(x, show_unavailable = TRUE)
  expect_identical(attr(p, "stormbreaker_curve_status")$metric_status, "incomplete")
  expect_no_warning(ggplot2::ggplot_build(p))
  expect_identical(p$data$y[[2]], NA_real_)
  expect_identical(nrow(ggplot2::ggplot_build(p)$data[[1]]), 0L)
})

test_that("unusable DCA rows do not create an interpolated line across a gap", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("survival")$cv
  x$metrics$decision_curve$curve$status[[3]] <- "not_evaluable"
  p <- sb_plot_dca(x)
  expect_identical(sum(p$data$draw), 12L)
  model_line <- p$layers[[1]]$data[p$layers[[1]]$data$strategy == "model", ]
  expect_identical(length(unique(model_line$segment)), 2L)
  expect_identical(p$data$y[3], x$metrics$decision_curve$curve$model_net_benefit[3])
})

test_that("DCA plotting cannot choose thresholds or recompute clinical utility", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture()$external
  before <- x
  local_mocked_bindings(
    sb_fit_final = function(...) stop("No fitting"),
    sb_evaluate_external = function(...) stop("No evaluation"),
    sb_binary_decision_curve_metric = function(...) stop("No recomputation"),
    sb_survival_decision_curve_metric = function(...) stop("No recomputation")
  )
  expect_s3_class(sb_plot_dca(x), "ggplot")
  expect_identical(x, before)
  expect_identical(names(formals(sb_plot_dca)), c("x", "cohort_label", "show_unavailable",
    "standardized", "cost_benefit_axis"))
})
})

# ---- plot-clinical ----
local({
# ---- plot-clinician-panel ----

test_that("combined panels retain ordered evidence and separate explanations", {
  skip_if_not_installed("patchwork")
  result <- make_sb_curve_view_fixture("binary")$cv
  x <- sb_clinician_views(result)
  before <- serialize(x, NULL)
  selected <- c("Calibration", "ROC")
  combined <- sb_plot_clinician_panel(x, selected)
  expect_s3_class(combined, "patchwork")
  for (i in seq_along(selected)) {
    original <- x$plots[[selected[i]]]
    expect_identical(combined[[i]]$data, original$data)
    expect_identical(combined[[i]]$layers, original$layers)
    expect_identical(combined[[i]]$labels$title, selected[i])
    expect_null(combined[[i]]$labels$subtitle)
    expect_null(combined[[i]]$labels$caption)
  }
  expect_identical(attr(combined, "stormbreaker_captions")$view, selected)
  expect_identical(serialize(x, NULL), before)
  expect_s3_class(patchwork::patchworkGrob(combined), "gtable")
  expect_error(sb_plot_clinician_panel(x, "Missing"), class = "stormbreaker_bad_report_input")
  expect_error(sb_plot_clinician_panel(x, c("ROC", "ROC")), class = "stormbreaker_bad_report_input")
})

# ---- probability-distribution ----

probability_distribution_fixture <- function() {
  structure(list(status = "success", task_type = "binary", pipeline_id = "distribution",
    scored_predictions = data.frame(patient_id = letters[1:6], outcome = c(1L, 0L, 1L, 0L, 1L, 0L),
      probability = c(0, .5, 1, 1, .5, 0), status = "success"), metrics = list()),
    class = "stormbreaker_external_evaluation")
}

test_that("probability bins preserve zero, one, ties and absent classes", {
  x <- probability_distribution_fixture()
  before <- x
  r <- sb_probability_distribution(x, c(0, .5, 1))
  expect_equal(r$table$count, c(1, 2, 1, 2))
  expect_equal(r$table$fraction, c(1, 2, 1, 2) / 3)
  expect_equal(sum(r$table$count), 6)
  expect_identical(x, before)
  expect_length(intersect("patient_id", names(r$table)), 0)
  expect_equal(sb_probability_distribution(x, c(0, 1))$table$count, c(3, 3))
  x$scored_predictions$outcome[] <- 1L
  q <- sb_probability_distribution(x, c(0, .5, 1))
  expect_equal(q$table$count, c(0, 0, 2, 4))
  expect_equal(q$table$fraction[1:2], c(NA_real_, NA_real_))
})

test_that("distributions preserve invalid groups and cohort separation", {
  x <- probability_distribution_fixture()
  bad <- x
  bad$scored_predictions$probability[1] <- NA_real_
  r <- sb_probability_distribution(list(A = x, B = bad), c(0, .5, 1))
  expect_equal(sum(r$table$count), 6)
  expect_equal(nrow(r$diagnostics), 2)
  expect_equal(sort(r$diagnostics$status), c("invalid_predictions", "success"))
  both <- sb_probability_distribution(list(A = x, B = x), c(0, .5, 1))
  expect_equal(length(unique(both$table$panel)), 2)
  expect_equal(sum(both$table$count), 12)
  for (b in list(c(0, 0, 1), c(.1, 1), c(0, NA, 1), c(0, .8, .5, 1))) {
    err <- tryCatch(sb_probability_distribution(x, b), error = identity)
    expect_s3_class(err, "stormbreaker_bad_visualization_argument")
  }
})

test_that("probability graphics use aggregate bins including empty bins", {
  skip_if_not_installed("ggplot2")
  r <- sb_probability_distribution(probability_distribution_fixture(), c(0, .2, .5, 1))
  p <- sb_plot_probability_distribution(r)
  expect_equal(nrow(p$data), 6)
  expect_equal(sum(p$data$value), 6)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(sb_plot_probability_distribution(r, "fraction")), "ggplot_built")
})
})

# ---- plot-comparison ----
local({
# ---- plot-performance-heatmap ----

make_sb_heatmap_fixture <- function(id = "model A", estimate = 0.75, status = "success") {
  structure(list(status = status, task_type = "binary", pipeline_id = id,
    patient_match = list(status = "success"), error_stage = NA_character_,
    error_class = NA_character_, error_message = NA_character_,
    metrics = list(roc_auc = list(estimate = estimate, status = "success",
      scope = "external_cohort", aggregation = "patient_pairs"))),
    class = "stormbreaker_external_evaluation")
}

test_that("performance heatmaps preserve named cohorts, model order and missing cells", {
  skip_if_not_installed("ggplot2")
  a <- make_sb_heatmap_fixture()
  b <- make_sb_heatmap_fixture("model B", 0.6)
  comparison <- structure(list(evaluations = list(a, b)), class = "stormbreaker_external_comparison")
  x <- list(Cohort1 = comparison, Cohort2 = a)
  before <- x
  p <- sb_plot_performance_heatmap(x, "roc_auc")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 4L)
  expect_identical(levels(p$data$pipeline_id), c("model B", "model A"))
  expect_equal(sum(p$data$cell_status == "not_supplied"), 1L)
  expect_equal(sort(p$data$estimate[is.finite(p$data$estimate)]), c(0.6, 0.75, 0.75))
  expect_identical(x, before)
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 4L)
})

test_that("failed estimates are retained but never painted as success", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_heatmap_fixture(status = "failed")
  p <- sb_plot_performance_heatmap(list(Cohort = x), "roc_auc")
  expect_equal(p$data$estimate, 0.75)
  expect_identical(p$data$fill_value, NA_real_)
  expect_identical(p$data$cell_status, "failed_run")
  expect_match(p$labels$subtitle, "No usable", fixed = TRUE)
  expect_error(sb_plot_performance_heatmap(list(a = x, a = x), "roc_auc"),
    class = "stormbreaker_bad_visualization_input")
})

test_that("time targets and validation scopes are not conflated", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_survival(45, seed = 14)
  results <- lapply(c(3, 6), function(horizon) {
    task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month",
      prediction_horizon = horizon)
    sb_run_cv(task, folds = 3, seed = 9)
  })
  p <- sb_plot_performance_heatmap(list(ThreeMonths = results[[1]], SixMonths = results[[2]]), "time_auc")
  expect_equal(length(unique(p$data$target_label)), 2L)
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 2L)
  expect_identical(sort(unique(p$data$prediction_horizon)), c(3, 6))
})

# ---- plot-performance ----

test_that("performance points use saved estimates in input order", {
  skip_if_not_installed("ggplot2")
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  x <- sb_compare_cv(task, list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), folds = 3, seed = 1)
  before <- x
  p <- sb_plot_performance(x, "roc_auc")
  expect_s3_class(p, "ggplot")
  expect_identical(p$data$estimate, sb_result_table(x, metric = "roc_auc")$estimate)
  expect_identical(levels(p$data$pipeline_id), names(x$results))
  expect_identical(p$scales$get_scales("x")$limits, rev(names(x$results)))
  expect_length(p$layers, 1L)
  expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
  expect_match(p$labels$subtitle, x$results[[1]]$metrics$roc_auc$aggregation,
               fixed = TRUE)
  expect_match(p$labels$caption, "No confidence intervals", fixed = TRUE)
  expect_identical(x, before)
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$data[[1]]$y, p$data$estimate)
})

test_that("failed runs and unavailable metrics remain visible without zero filling", {
  skip_if_not_installed("ggplot2")
  d <- make_sb_fixture()
  task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month")
  a <- sb_run_cv(task, folds = 3, seed = 1)
  b <- a
  b$pipeline$id <- "failed_run"
  b$status <- "failed"
  c <- a
  c$pipeline$id <- "unavailable_metric"
  c$metrics$harrell_c$status <- "not_evaluable"
  c$metrics$harrell_c$estimate <- NA_real_
  x <- structure(list(results = list(b, a, c)),
                 class = "stormbreaker_cv_comparison")
  p <- sb_plot_performance(x, "harrell_c")
  expect_identical(nrow(p$data), 3L)
  expect_identical(nrow(p$layers[[1]]$data), 1L)
  expect_identical(p$data$estimate,
                   c(a$metrics$harrell_c$estimate, a$metrics$harrell_c$estimate, NA_real_))
  expect_identical(levels(p$data$pipeline_id), c("failed_run", a$pipeline$id,
                                              "unavailable_metric"))
  expect_match(p$labels$caption, "2 pipeline", fixed = TRUE)
  x$results <- list(b, c)
  empty <- sb_plot_performance(x, "harrell_c", show_unavailable = TRUE)
  expect_identical(nrow(empty$layers[[1]]$data), 0L)
  expect_identical(nrow(empty$data), 2L)
  expect_identical(nrow(ggplot2::ggplot_build(empty)$data[[1]]), 0L)
})

test_that("plotting keeps infinite loss in its table but does not draw it", {
  skip_if_not_installed("ggplot2")
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  x <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1)
  x$metrics$log_loss$estimate <- Inf
  p <- sb_plot_performance(x, "log_loss", show_unavailable = TRUE)
  expect_identical(p$data$estimate, Inf)
  expect_identical(nrow(p$layers[[1]]$data), 0L)
  expect_match(p$labels$caption, "1 pipeline", fixed = TRUE)
})

test_that("plotting is read-only and works for each supported binary metric", {
  skip_if_not_installed("ggplot2")
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  x <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1)
  before <- x
  local_mocked_bindings(
    sb_run_cv = function(...) stop("must not run CV"),
    sb_fit_final = function(...) stop("must not fit"),
    sb_binary_auc_metric = function(...) stop("must not recompute")
  )
  for (metric in c("roc_auc", "average_precision", "brier_score", "log_loss")) {
    p <- sb_plot_performance(x, metric)
    expect_identical(p$data$estimate, x$metrics[[metric]]$estimate)
    expect_identical(nrow(ggplot2::ggplot_build(p)$data[[1]]), 1L)
  }
  expect_identical(x, before)
})

# ---- plot-selection-stability ----

make_sb_selection_plot_fixture <- function(failed = FALSE) {
  selections <- list(
    list(fold = 1L, status = "success", selected_features = c("a", "b")),
    list(fold = 2L, status = "success", selected_features = c("b", "c")),
    list(fold = 3L, status = if (failed) "failed" else "success", selected_features = c("a", "b")))
  structure(list(status = if (failed) "failed" else "success", pipeline = list(id = "selected"),
    selections = selections, selection_stability = sb_selection_stability("test", selections, c("a", "b", "c", "d"))),
    class = "stormbreaker_binary_cv_result")
}

test_that("selection frequency copies saved denominators and labels incomplete selections", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_selection_plot_fixture(TRUE)
  before <- x
  p <- sb_plot_selection_frequency(x)
  expect_equal(p$data$selection_frequency, c(0.5, 1, 0.5, 0))
  expect_identical(unique(p$data$n_successful_folds), 2L)
  expect_match(as.character(p$data$panel[1]), "2/3", fixed = TRUE)
  expect_match(as.character(p$data$panel[1]), "incomplete", fixed = TRUE)
  expect_identical(x, before)
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 4L)
})

test_that("Jaccard heatmap uses saved fold pairs without invented diagonal", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_selection_plot_fixture()
  p <- sb_plot_selection_jaccard(x)
  expect_equal(p$data$jaccard, x$selection_stability$pairwise$jaccard)
  expect_identical(p$data$fold_i, x$selection_stability$pairwise$fold_i)
  expect_equal(nrow(p$data), 3L)
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 3L)
})

test_that("UpSet intersections come from exact saved successful fold sets", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_selection_plot_fixture()
  p <- sb_plot_selection_upset(x)
  intersections <- attr(p, "stormbreaker_intersections")
  expect_equal(sum(intersections$count), 3L)
  expect_setequal(intersections$membership, c("1,3", "1,2,3", "2"))
  expect_setequal(unlist(intersections$features, use.names = FALSE), c("a", "b", "c"))
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  p <- sb_plot_selection_upset(make_sb_selection_plot_fixture(TRUE))
  expect_equal(attr(p, "stormbreaker_selection_status")$n_successful_folds, 2L)
  expect_match(p$labels$caption, "successful folds", fixed = TRUE)
})

test_that("non-applicable stability does not invent a selected set", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_selection_plot_fixture()
  x$selection_stability$status <- "not_applicable"
  for (fun in list(sb_plot_selection_frequency, sb_plot_selection_jaccard, sb_plot_selection_upset)) {
    p <- fun(x)
    expect_match(p$labels$subtitle, "No usable", fixed = TRUE)
    expect_identical(attr(p, "stormbreaker_selection_status")$status, "not_applicable")
  }
})
})

# ---- plot-discrimination ----
local({
# ---- plot-discrimination ----

make_sb_discrimination_fixture <- function(external = FALSE) {
  predictions <- data.frame(patient_id = paste0("p", 1:4), outcome = c(1L, 0L, 1L, 0L),
    probability = c(0.9, 0.7, 0.7, 0.1), status = "success", stringsAsFactors = FALSE)
  metrics <- list(
    roc_auc = list(estimate = 0.875, status = "success", scope = "fixture", aggregation = "fixture"),
    average_precision = list(estimate = 5/6, status = "success", scope = "fixture", aggregation = "fixture"))
  if (external) return(structure(list(status = "success", task_type = "binary",
    pipeline_id = "test", scored_predictions = predictions, metrics = metrics),
    class = "stormbreaker_external_evaluation"))
  predictions$fold <- 1L
  predictions$pipeline_id <- "test"
  structure(list(status = "success", task_type = "binary", pipeline = list(id = "test"),
    predictions = predictions, metrics = metrics), class = "stormbreaker_binary_cv_result")
}

test_that("ROC uses tied thresholds and preserves saved estimands", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_discrimination_fixture()
  before <- x
  p <- sb_plot_roc(x)
  expect_s3_class(p, "ggplot")
  expect_equal(p$data$x, c(0, 0, 0.5, 1))
  expect_equal(p$data$y, c(0, 0.5, 1, 1))
  expect_equal(sum(diff(p$data$x) * (head(p$data$y, -1) + tail(p$data$y, -1))/2), 0.875)
  expect_identical(attr(p, "stormbreaker_saved_metrics")[[1]]$estimate, 0.875)
  expect_identical(x, before)
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 1L)
})

test_that("PR uses noninterpolated AP step geometry", {
  skip_if_not_installed("ggplot2")
  p <- sb_plot_pr(make_sb_discrimination_fixture(external = TRUE))
  expect_equal(p$data$x, c(0, 0.5, 1, 1))
  expect_equal(p$data$y, c(1, 1, 2/3, 0.5))
  expect_equal(sum(diff(p$data$x) * tail(p$data$y, -1)), 5/6)
  expect_match(p$labels$caption, "not pooled", fixed = TRUE)
  expect_identical(p$layers[[1]]$geom_params$direction, "vh")
})

test_that("optional per-fold ROC preserves folds and invalid-row diagnostics", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_discrimination_fixture()
  second <- x$predictions
  second$patient_id <- paste0("q", 1:4)
  second$fold <- 2L
  second$probability <- rev(second$probability)
  x$predictions <- rbind(x$predictions, second)
  p <- sb_plot_roc(x, by_fold = TRUE)
  expect_equal(length(unique(p$data$fold)), 2L)
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 2L)
  x$predictions$probability[1] <- NA_real_
  p <- sb_plot_roc(x, by_fold = TRUE)
  status <- attr(p, "stormbreaker_curve_status")
  expect_identical(status$curve_status, c("invalid_predictions", "success"))
  expect_identical(unique(p$data$fold), 2L)
  expect_match(attr(p, "stormbreaker_figure_notes"), "1 unavailable", fixed = TRUE)
})

test_that("failed and single-class curves remain explicit", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_discrimination_fixture()
  x$status <- "failed"
  p <- sb_plot_roc(x)
  expect_identical(attr(p, "stormbreaker_curve_status")$curve_status, "failed_run")
  expect_match(p$labels$subtitle, "No usable", fixed = TRUE)
  x$status <- "success"
  x$predictions$outcome <- 1L
  expect_identical(attr(sb_plot_roc(x), "stormbreaker_curve_status")$curve_status, "single_class")
  expect_equal(tail(sb_plot_pr(x)$data$y, 1), 1)
})

test_that("named cohorts and composite repetitions retain their boundaries", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_discrimination_fixture()
  comparison <- structure(list(results = list(x)), class = "stormbreaker_cv_comparison")
  repeated <- structure(list(repetitions = list(comparison, comparison)),
    class = "stormbreaker_repeated_cv_comparison")
  p <- sb_plot_roc(list(Development = repeated, External = make_sb_discrimination_fixture(TRUE)))
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 3L)
  expect_identical(unique(p$data$repeat_id), c(1L, 2L, NA_integer_))
  expect_equal(length(unique(p$data$series_id)), 3L)
  expect_error(sb_plot_roc(list(x, x)), class = "stormbreaker_bad_visualization_input")
})

# ---- plot-prediction-confusion ----

test_that("direct confusion counts and cell orientation match supplied predictions", {
  y <- c(rep(0, 72), rep(1, 199))
  prob <- c(rep(.1, 60), rep(.9, 12), rep(.1, 8), rep(.9, 191))
  p <- sb_plot_prediction_confusion(y, prob, "apparent", cohort = "Train",
    class_labels = c("Control", "Treat"))
  expect_equal(p$data$count, c(60, 12, 8, 191))
  b <- ggplot2::ggplot_build(p)
  expect_equal(as.numeric(b$data[[1]]$x), c(1, 1, 2, 2))
  expect_equal(as.numeric(b$data[[1]]$y), c(2, 1, 2, 1))
  expect_equal(nrow(b$data[[2]]), 4L)
  expect_identical(anyDuplicated(b$data[[2]][c("x", "y")]), 0L)
  context <- attr(p, "stormbreaker_prediction_context")
  expect_identical(context$kind, "apparent")
  expect_equal(context$n, 271)
  expect_equal(context$counts$tp, 191)
  expect_null(p$labels$subtitle)
  expect_null(p$labels$caption)
  expect_identical(p$labels$title, "Train")
})

test_that("threshold ties and absent classes retain all cells", {
  p <- sb_plot_prediction_confusion(c(0, 1), c(.5, .5), "external")
  expect_equal(p$data$count, c(0, 1, 0, 1))
  q <- sb_plot_prediction_confusion(c(0, 0), c(.1, .9), "oof", normalize = "observed")
  expect_equal(q$data$value, c(.5, .5, NA, NA))
  expect_equal(nrow(ggplot2::ggplot_build(q)$data[[2]]), 4)
  expect_equal(sb_plot_prediction_confusion(c(0, 1), c(0, 1), "external", threshold = 1)$data$count,
    c(1, 0, 0, 1))
})

test_that("direct confusion rejects invalid inputs rather than dropping rows", {
  expect_error(sb_plot_prediction_confusion(c(0, 1), .5, "apparent"), "aligned")
  expect_error(sb_plot_prediction_confusion(c(0, NA), c(.1, .2), "apparent"), "aligned")
  expect_error(sb_plot_prediction_confusion(c(0, 1), c(.1, 2), "apparent"), "aligned")
  expect_error(sb_plot_prediction_confusion(c(0, 1), c(.1, .2), "apparent", threshold = NA_real_), "threshold")
  expect_error(sb_plot_prediction_confusion(c(0, 1), c(.1, .2), "apparent", class_labels = c("A", "A")), "distinct")
})
})

# ---- plot-external-uncertainty ----
local({
make_sb_plot_bootstrap <- function(id = "model", metric = "roc_auc") {
  structure(list(status = "success", pipeline_id = id, task_type = "binary",
    confidence_level = 0.95, prediction_target = list(has_horizon = FALSE),
    intervals = data.frame(metric = metric, point_estimate = 0.7, lower = 0.6, upper = 0.8,
      confidence_level = 0.95, status = "success", interval_available = TRUE,
      replicate_completeness = "complete", n_requested = 100L, n_finite_successful = 100L)),
    class = "stormbreaker_external_bootstrap")
}

test_that("external interval display copies endpoints and bootstrap completeness", {
  skip_if_not_installed("ggplot2")
  a <- make_sb_plot_bootstrap()
  b <- a
  b$status <- "incomplete"
  b$intervals$replicate_completeness <- "partial"
  b$intervals$n_finite_successful <- 80L
  b$intervals$status <- "incomplete"
  before <- b
  p <- sb_plot_external_intervals(list(A = a, B = b), "roc_auc")
  expect_s3_class(p, "ggplot")
  expect_identical(p$data$lower, c(0.6, 0.6))
  expect_identical(p$data$upper, c(0.8, 0.8))
  expect_match(as.character(p$data$row_label[2]), "80/100", fixed = TRUE)
  expect_identical(b, before)
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 2L)
})

test_that("paired differences retain reference direction and failed candidates", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_plot_bootstrap()
  class(x) <- "stormbreaker_external_bootstrap_comparison"
  x$reference_pipeline_id <- "reference"
  x$intervals$candidate_pipeline_id <- "candidate"
  x$intervals$reference_pipeline_id <- "reference"
  x$intervals$metric_direction <- "higher_is_better"
  x$candidate_status <- data.frame(candidate_pipeline_id = c("candidate", "broken"), status = c("success", "failed"))
  p <- sb_plot_external_difference(x, "roc_auc")
  expect_equal(nrow(p$data), 2L)
  expect_identical(p$data$reference_pipeline_id, c("reference", "reference"))
  expect_identical(p$data$status, c("success", "failed"))
  expect_match(p$labels$x, "candidate - reference", fixed = TRUE)
  expect_equal(ggplot2::ggplot_build(p)$data[[3]]$xintercept, 0)
})

test_that("forest plots distinguish saved cohort CI, pooled CI and prediction interval", {
  skip_if_not_installed("ggplot2")
  a <- make_sb_plot_bootstrap()
  b <- a
  b$intervals <- NULL
  x <- structure(list(status = "success", task_type = "binary", metric = "roc_auc",
    confidence_level = 0.95, pipeline_id = "model", metric_target = list(has_horizon = FALSE),
    cohorts = data.frame(cohort_id = c("A", "B"), point_estimate = c(0.7, 0.72)),
    summary = data.frame(estimate = 0.71, lower = 0.63, upper = 0.78,
      prediction_lower = 0.55, prediction_upper = 0.83),
    bootstraps = list(A = a, B = b), prediction_interval_caution = "Only two cohorts."),
    class = "stormbreaker_external_discrimination_meta")
  before <- x
  p <- sb_plot_meta_forest(x)
  expect_equal(p$data$point_estimate, c(0.7, 0.72, 0.71, 0.71))
  expect_equal(p$data$lower, c(0.6, NA, 0.63, 0.55))
  expect_identical(p$data$interval_type, c("Cohort bootstrap CI", "Cohort bootstrap CI", "Pooled CI", "Prediction interval"))
  expect_identical(p$data$status[2], "interval_not_saved")
  expect_match(p$labels$caption, "Only two cohorts", fixed = TRUE)
  expect_identical(x, before)
  expect_equal(nrow(sb_plot_meta_forest(x, show_prediction_interval = FALSE)$data), 3L)
})

test_that("real external bootstrap output can be plotted without re-evaluation", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_binary(60, seed = 17)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  evaluation <- sb_evaluate_external(sb_predict_final(fit, d, "patient_id"),
    sb_binary_endpoint(d, "patient_id", "outcome"))
  b <- sb_bootstrap_external(evaluation, reps = 10, seed = 9)
  local_mocked_bindings(sb_bootstrap_external = function(...) stop("Do not bootstrap"),
    sb_evaluate_external = function(...) stop("Do not evaluate"))
  p <- sb_plot_external_intervals(b, "brier_score")
  expected <- b$intervals[b$intervals$metric == "brier_score", ]
  expect_identical(p$data$lower, expected$lower)
  expect_identical(p$data$point_estimate, expected$point_estimate)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

make_sb_plot_missing_bootstrap <- function(id = "missing", metric = "roc_auc") {
  x <- make_sb_plot_bootstrap(id, metric)
  x$status <- "incomplete"
  x$intervals[c("point_estimate", "lower", "upper")] <- NA_real_
  x$intervals$status <- "incomplete"
  x$intervals$interval_available <- FALSE
  x$intervals$replicate_completeness <- "none"
  x$intervals$n_finite_successful <- 0L
  x
}

test_that("all unavailable intervals retain visible cohort and status rows", {
  skip_if_not_installed("ggplot2")
  withr::local_pdf(NULL)
  inputs <- list(Alpha = make_sb_plot_missing_bootstrap("first"),
    Beta = make_sb_plot_missing_bootstrap("second"))
  before <- inputs
  p <- sb_plot_external_intervals(inputs, "roc_auc")
  built <- ggplot2::ggplot_build(p)
  labels <- as.character(built$layout$panel_params[[1]]$y$get_labels())
  expect_setequal(labels, sb_view_wrap_label(as.character(p$data$row_label)))
  expect_length(labels, 2L)
  expect_equal(sum(grepl("0/100", labels, fixed = TRUE)), 2L)
  expect_equal(sum(grepl("incomplete", labels, fixed = TRUE)), 2L)
  expect_equal(nrow(built$data[[1]]), 0L)
  expect_equal(nrow(built$data[[2]]), 0L)
  expect_identical(p$data$point_estimate, c(NA_real_, NA_real_))
  expect_identical(inputs, before)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})

test_that("unavailable interval rows stay in their own target panels", {
  skip_if_not_installed("ggplot2")
  good <- make_sb_plot_bootstrap("available", "calibration_slope")
  missing <- make_sb_plot_missing_bootstrap("unavailable", "calibration_slope")
  survival <- missing
  survival$task_type <- "survival"
  survival$prediction_target <- list(has_horizon = TRUE,
    prediction_horizon = 24, time_unit = "month")
  p <- sb_plot_external_intervals(list(Binary = good,
    Binary_missing = missing, Survival_missing = survival), "calibration_slope")
  built <- ggplot2::ggplot_build(p)
  panels <- built$layout$layout
  expect_equal(nrow(panels), 2L)
  for (i in seq_len(nrow(panels))) {
    rows <- as.character(p$data$target_label) == as.character(panels$target_label[i])
    labels <- as.character(built$layout$panel_params[[i]]$y$get_labels())
    expect_setequal(labels, sb_view_wrap_label(as.character(p$data$row_label[rows])))
  }
  available <- ggplot2::ggplot_build(sb_plot_external_intervals(good, "calibration_slope"))
  expect_equal(built$layout$panel_params[[1]]$x.range,
    available$layout$panel_params[[1]]$x.range)
  expect_equal(built$data[[1]]$x, 0.6)
  expect_equal(built$data[[1]]$xend, 0.8)
  expect_equal(built$data[[2]]$x, 0.7)
})

test_that("all failed paired candidates keep labels and only the reference line", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_plot_missing_bootstrap()
  class(x) <- "stormbreaker_external_bootstrap_comparison"
  x$reference_pipeline_id <- "reference"
  x$intervals$candidate_pipeline_id <- "first"
  x$intervals$reference_pipeline_id <- "reference"
  x$intervals$status <- "failed"
  x$candidate_status <- data.frame(candidate_pipeline_id = c("first", "second"),
    status = c("failed", "failed"))
  p <- sb_plot_external_difference(x, "roc_auc")
  built <- ggplot2::ggplot_build(p)
  labels <- as.character(built$layout$panel_params[[1]]$y$get_labels())
  expect_setequal(labels, sb_view_wrap_label(as.character(p$data$row_label)))
  expect_length(labels, 2L)
  expect_equal(sum(grepl("Reference: reference", labels, fixed = TRUE)), 2L)
  expect_equal(nrow(built$data[[1]]), 0L)
  expect_equal(nrow(built$data[[2]]), 0L)
  expect_equal(built$data[[3]]$xintercept, 0)
})
})

# ---- plot-pages ----
local({
test_that("empty plot collections remain empty", {
  expect_identical(sb_plot_pages(list()), list())
})

test_that("two-column feature facets paginate without losing support rows", {
  d <- expand.grid(feature = paste0("gene_", 1:12), value = 1:5)
  d$effect <- d$value / 2
  p <- ggplot2::ggplot(d, ggplot2::aes(value, effect)) + ggplot2::geom_line() +
    ggplot2::geom_rug(data = d) + ggplot2::facet_wrap(~feature, ncol = 2)
  pages <- sb_plot_pages(p)
  expect_length(pages, 2)
  expect_equal(sum(vapply(pages, function(z) nrow(z$data), integer(1))), nrow(d))
  expect_equal(sum(vapply(pages, function(z) nrow(z$layers[[2]]$data), integer(1))), nrow(d))
  expect_equal(vapply(pages, function(z) z$facet$params$ncol, numeric(1)), c(plot1 = 2, plot2 = 2), ignore_attr = TRUE)
})

test_that("query pagination retains every observation and ribbon row", {
  d <- data.frame(query = rep(1:15, each = 3), time = rep(1:3, 15), value = seq_len(45))
  ribbon <- transform(d, lower = value - 1, upper = value + 1)
  p <- ggplot2::ggplot(d, ggplot2::aes(time, value)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(data = ribbon, ggplot2::aes(ymin = lower, ymax = upper)) +
    ggplot2::facet_wrap(~query)
  pages <- sb_plot_pages(p)
  expect_length(pages, 3)
  expect_equal(do.call(rbind, lapply(pages, function(z) z$data))$value, d$value)
  expect_equal(sum(vapply(pages, function(z) nrow(z$layers[[2]]$data), integer(1))), nrow(ribbon))
  expect_equal(nrow(p$layers[[2]]$data), 45L)
  expect_equal(vapply(pages, function(z) length(unique(z$data$query)), integer(1)),
    c("plot / 1" = 6L, "plot / 2" = 6L, "plot / 3" = 3L))
  expect_length(lapply(pages, ggplot2::ggplot_build), 3)
})

test_that("labelled posterior plots paginate and retain separate notes", {
  d <- data.frame(query = 1:60, mean = seq_len(60) / 61)
  p <- ggplot2::ggplot(d, ggplot2::aes(factor(query), mean)) + ggplot2::geom_point() +
    ggplot2::labs(caption = "Pointwise 95% posterior credible intervals")
  p <- sb_compact_plot(p)
  pages <- sb_plot_pages(p)
  expect_length(pages, 3)
  expect_null(p$labels$caption)
  expect_identical(attr(pages[[1]], "stormbreaker_notes"), "Pointwise 95% posterior credible intervals")
  expect_equal(sum(vapply(pages, function(z) nrow(z$data), integer(1))), 60L)
})
})

# ---- survival-views ----
local({
# ---- paper-figures ----

test_that("pooled ROC has one curve and direct scores have reproducible intervals", {
  d <- sb_example_binary(80, seed = 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  cv <- sb_run_cv(task, folds = 3)
  p <- sb_plot_roc(cv, ci_reps = 20)
  expect_length(unique(p$data$panel), 1L)
  expect_equal(attr(p, "stormbreaker_roc_statistics")$auc,
    sb_roc_interval(cv$predictions$outcome, cv$predictions$probability)$auc)
  expect_length(unique(sb_plot_roc(cv, by_fold = TRUE)$data$panel), 3L)
  before <- cv
  cv$predictions$probability[1] <- NA_real_
  expect_identical(attr(sb_plot_roc(cv), "stormbreaker_curve_status")$curve_status, "invalid_predictions")
  y <- c(1, 0, 1, 0); score <- c(.9, .7, .7, .1)
  a <- sb_plot_prediction_roc(y, score, kind = "external", ci_reps = 40)
  b <- sb_plot_prediction_roc(y, score, kind = "external", ci_reps = 40)
  expect_equal(attr(a, "stormbreaker_roc_statistics")$auc, .875)
  expect_identical(attr(a, "stormbreaker_roc_statistics"), attr(b, "stormbreaker_roc_statistics"))
  expect_identical(attr(a, "prediction_origin"), "external")
  expect_equal(nrow(ggplot2::ggplot_build(a)$layout$layout), 1)
})

test_that("ranking uses only declared cohorts and leaves incomplete models unranked", {
  skip_if_not_installed("patchwork")
  d <- data.frame(model = rep(c("a", "b", "c"), each = 3),
    cohort = rep(c("Train", "Validation", "Test"), 3),
    auc = c(.9, .8, .5, .8, .9, .4, .99, NA, 1))
  z <- sb_plot_model_ranking(d, "Validation", top_n = 2)
  expect_identical(z$ranking$model, c("b", "a", "c"))
  expect_equal(z$ranking$mean_auc, c(.9, .8, NA))
  expect_identical(z$ranking$status[3], "incomplete")
  expect_s3_class(z$plot, "patchwork")
  expect_equal(nrow(z$data), 6)
})

test_that("bootstrap calibration refits each sample and computes the correction", {
  d <- sb_example_binary(100, seed = 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  original <- sb_fit_final
  calls <- 0L
  local_mocked_bindings(sb_fit_final = function(...) {
    calls <<- calls + 1L
    original(...)
  })
  z <- sb_bootstrap_calibration(task, sb_pipeline(learner = "logistic"), reps = 12, seed = 7)
  expect_identical(calls, 13L)
  expect_identical(z$status$status, rep("success", 12))
  expect_equal(z$curve$corrected, z$curve$apparent - z$curve$optimism)
  expect_equal(z$curve$effective_reps, colSums(is.finite(z$replicates)))
  expect_gt(sum(is.finite(z$curve$corrected)), 20)
  expect_gt(max(abs(z$curve$optimism), na.rm = TRUE), 0)
  again <- sb_bootstrap_calibration(task, sb_pipeline(learner = "logistic"), reps = 12, seed = 7)
  expect_identical(z$curve, again$curve)
  expect_s3_class(sb_plot_bootstrap_calibration(z), "ggplot")
})

test_that("nomogram density and patient marks preserve original probabilities", {
  d <- sb_example_binary(100, seed = 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  nom <- sb_nomogram(fit, d)
  p <- sb_plot_nomogram(nom, reference_data = d, patient = d[1, , drop = FALSE])
  expect_equal(attr(p, "stormbreaker_nomogram_patient")$probability,
    sb_predict_final(fit, d[1, c("patient_id", "gene_a", "gene_b")], "patient_id")$predictions$probability)
  expect_gt(nrow(attr(p, "stormbreaker_nomogram_density")), 100)
  expect_null(p$labels$caption)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that("standardized binary DCA divides by prevalence without changing inputs", {
  d <- sb_example_binary(80, seed = 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"), decision_thresholds = c(.1, .2, .5, .8))
  cv <- sb_run_cv(task, folds = 3)
  before <- serialize(cv, NULL)
  ordinary <- sb_plot_dca(cv)
  standard <- sb_plot_dca(cv, standardized = TRUE, cost_benefit_axis = TRUE)
  expect_equal(standard$data$y, ordinary$data$y / mean(d$outcome))
  expect_identical(serialize(cv, NULL), before)
  expect_identical(standard$labels$y, "Standardized net benefit")
  expect_s3_class(ggplot2::ggplot_build(standard), "ggplot_built")
})

# ---- survival-cohort-summary ----

test_that("cohort bars preserve supplied metric values", {
  d <- data.frame(cohort = c("A", "B"), c_index = c(.68, .75))
  p <- sb_plot_cohort_metrics(d)
  expect_equal(p$data$value, d$c_index)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})

test_that("risk meta agrees with native Cox and inverse-variance estimates", {
  withr::local_seed(102)
  cohorts <- lapply(1:3, function(i) {
    risk <- rnorm(120)
    d <- data.frame(patient_id = paste(i, 1:120), risk = risk,
      time = rexp(120, exp(.5 * risk)), event = rbinom(120, 1, .8),
      group = ifelse(risk > 0, "High", "Low"))
    list(patients = d, cutoff = 0)
  })
  names(cohorts) <- c("A", "B", "C")
  for (type in c("group", "continuous")) {
    m <- surv_meta_risk(cohorts, type)
    expect_identical(m$status, "success")
    e <- m$estimates
    expect_equal(m$pooled$log_hr[1], weighted.mean(e$log_hr, 1/e$se^2))
    expect_equal(sum(e$weight_random), 100)
    expect_s3_class(ggplot2::ggplotGrob(m$plot), "gtable")
  }
  d <- cohorts[[1]]$patients
  native <- survival::coxph(survival::Surv(time, event) ~ risk, d)
  expect_equal(sb_risk_cox_estimate(d, "continuous")$log_hr, unname(coef(native)))
  cohorts[[1]]$patients$event <- 0
  expect_identical(surv_meta_risk(cohorts)$status, "not_evaluable")
  cohorts[[1]] <- cohorts[[2]]
  expect_error(surv_meta_risk(cohorts), "overlap")
})

# ---- survival-metric-views ----

test_that("horizon metric tables preserve native CV and external values", {
  fixture <- make_sb_curve_view_fixture("survival", comparison = TRUE)
  before <- fixture
  for (name in c("cv", "external")) {
    x <- fixture[[name]]
    results <- if (name == "cv") x$results else x$evaluations
    for (metric in c("time_auc", "uno_c", "ipcw_brier")) {
      tab <- sb_result_table(x, metric = metric)
      expect_identical(tab$estimate, unname(vapply(results, function(z) z$metrics[[metric]]$estimate, numeric(1))))
      expect_identical(tab$status, unname(vapply(results, function(z) z$metrics[[metric]]$status, character(1))))
      expect_identical(tab$scope, unname(vapply(results, function(z) z$metrics[[metric]]$scope, character(1))))
      expect_identical(tab$aggregation, unname(vapply(results, function(z) z$metrics[[metric]]$aggregation, character(1))))
      expect_identical(tab$prediction_horizon, rep(6, 2))
      expect_identical(tab$time_unit, rep("month", 2))
      expect_identical(tab$metric, rep(metric, 2))
    }
  }
  expect_identical(fixture, before)
})

test_that("horizon-free requests and unavailable saved results stay distinct", {
  d <- sb_example_survival(60, seed = 1)
  task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month")
  x <- sb_run_cv(task, folds = 3, seed = 1)
  for (metric in c("time_auc", "uno_c", "ipcw_brier")) {
    tab <- sb_result_table(x, metric = metric)
    expect_identical(tab$status, "not_requested")
    expect_identical(tab$estimate, NA_real_)
    expect_identical(tab$prediction_horizon, NA_real_)
    expect_identical(tab$time_unit, "month")
  }
  x$task$prediction_horizon <- 6
  expect_identical(sb_result_table(x, metric = "time_auc")$status, "not_saved")
  x$status <- "failed"
  expect_identical(sb_result_table(x, metric = "time_auc")$run_status, "failed")
})

test_that("horizon views never fit or recompute metrics", {
  x <- make_sb_curve_view_fixture("survival")
  before <- x
  local_mocked_bindings(
    sb_run_cv = function(...) stop("must not run CV"),
    sb_fit_final = function(...) stop("must not fit"),
    sb_evaluate_external = function(...) stop("must not evaluate"),
    sb_survival_time_auc_metric = function(...) stop("must not recompute"),
    sb_survival_uno_c_metric = function(...) stop("must not recompute"),
    sb_survival_ipcw_brier_metric = function(...) stop("must not recompute")
  )
  for (metric in c("time_auc", "uno_c", "ipcw_brier")) {
    expect_identical(sb_result_table(x$cv, metric = metric)$estimate, x$cv$metrics[[metric]]$estimate)
    expect_identical(sb_result_table(x$external, metric = metric)$estimate, x$external$metrics[[metric]]$estimate)
  }
  expect_identical(x, before)
})

test_that("horizon plots retain values above one and distinct targets", {
  skip_if_not_installed("ggplot2")
  a <- make_sb_curve_view_fixture("survival")$cv
  a$metrics$ipcw_brier$estimate <- 1.25
  b <- a
  b$pipeline$id <- "second_target"
  b$task$prediction_horizon <- 12
  b$metrics$ipcw_brier$prediction_horizon <- 12
  b$metrics$ipcw_brier$estimate <- 0.5
  x <- structure(list(results = list(a, b)), class = "stormbreaker_cv_comparison")
  before <- x
  p <- sb_plot_performance(x, "ipcw_brier")
  built <- ggplot2::ggplot_build(p)
  expect_equal(sort(built$data[[1]]$y), c(0.5, 1.25))
  expect_identical(p$labels$y, "IPCW Brier score")
  expect_identical(p$data$prediction_horizon, c(6, 12))
  labels <- p$scales$get_scales("x")$limits
  expect_identical(any(grepl("6 month", labels, fixed = TRUE)), TRUE)
  expect_identical(any(grepl("12 month", labels, fixed = TRUE)), TRUE)
  expect_identical(x, before)
})

test_that("horizon plots omit unavailable methods while retaining diagnostics", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("survival", comparison = TRUE)$cv
  x$results[[1]]$metrics$time_auc$status <- "incomplete"
  x$results[[1]]$metrics$time_auc$estimate <- NA_real_
  p <- sb_plot_performance(x, "time_auc")
  expect_identical(nrow(p$data), 2L)
  expect_identical(length(p$scales$get_scales("x")$limits), 1L)
  audit <- sb_plot_performance(x, "time_auc", show_unavailable = TRUE)
  expect_identical(length(audit$scales$get_scales("x")$limits), 2L)
  x$results[[2]]$metrics$time_auc$estimate <- Inf
  empty <- ggplot2::ggplot_build(sb_plot_performance(x, "time_auc"))
  expect_match(empty$data[[1]]$label, "No usable results", fixed = TRUE)
})

test_that("report includes saved horizon metrics and leaves source failures intact", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("survival")
  sections <- list(CV = x$cv, External = x$external)
  before <- serialize(sections, NULL)
  failure_before <- sb_failure_table(sections)
  output <- tempfile(fileext = ".html")
  sb_render_report(sections, output)
  html <- paste(readLines(output, warn = FALSE), collapse = "\n")
  for (metric in c("time_auc", "uno_c", "ipcw_brier")) {
    expect_identical(grepl(paste0("<summary>", metric, "</summary>"), html, fixed = TRUE), TRUE)
  }
  expect_match(html, "prediction_horizon", fixed = TRUE)
  expect_match(html, "time_unit", fixed = TRUE)
  expect_identical(serialize(sections, NULL), before)
  expect_identical(sb_failure_table(sections), failure_before)
  expect_identical(file.exists(sub(".html$", "_failures.csv", output)), TRUE)
})

test_that("external missing horizon values and original errors are preserved", {
  fixture <- make_sb_curve_view_fixture("survival")
  external <- sb_example_survival(60, seed = 3)
  external$patient_id <- paste0("NEW_", external$patient_id)
  fit <- sb_fit_final(fixture$task, sb_pipeline())
  prediction <- sb_predict_final(fit, external[c("patient_id", "gene_a", "gene_b")], "patient_id")
  endpoint <- sb_surv_endpoint(external, "patient_id", "time", "event", "month")
  evaluated <- sb_evaluate_external(prediction, endpoint)
  expect_identical(is.null(evaluated$endpoint$prediction_horizon), TRUE)
  for (metric in c("time_auc", "uno_c", "ipcw_brier")) {
    expect_identical(sb_result_table(evaluated, metric = metric)$prediction_horizon, 6)
  }
  prediction$predictions$patient_id[1] <- "missing_external_patient"
  failed <- sb_evaluate_external(prediction, endpoint)
  before <- failed
  expect_identical(failed$status, "failed")
  tab <- sb_result_table(failed, metric = "time_auc")
  expect_identical(tab$status, "not_saved")
  expect_identical(tab$estimate, NA_real_)
  expect_identical(tab$error_message, failed$error_message)
  expect_identical(tab$error_class, failed$error_class)
  expect_identical(failed, before)
})

test_that("horizon metrics remain survival-only and distinguish units", {
  skip_if_not_installed("ggplot2")
  binary <- make_sb_curve_view_fixture("binary")
  for (x in binary[c("cv", "external")]) {
    error <- tryCatch(sb_result_table(x, metric = "time_auc"), error = identity)
    expect_s3_class(error, "stormbreaker_result_metric_not_saved")
  }
  a <- make_sb_curve_view_fixture("survival")$external
  b <- a
  b$metrics$time_auc$time_unit <- "day"
  x <- structure(list(evaluations = list(a, b)), class = "stormbreaker_external_comparison")
  tab <- sb_result_table(x, metric = "time_auc")
  expect_identical(tab$time_unit, c("month", "day"))
  p <- sb_plot_performance(x, "time_auc")
  expect_length(p$scales$get_scales("x")$limits, 2)
  expect_identical(nrow(ggplot2::ggplot_build(p)$data[[1]]), 2L)
})

test_that("horizon-free HTML does not invent time-specific sections", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_survival(60, seed = 1)
  task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month")
  cv <- sb_run_cv(task, folds = 3, seed = 1)
  fit <- sb_fit_final(task, sb_pipeline())
  ext <- d
  ext$patient_id <- paste0("EXT_", ext$patient_id)
  pred <- sb_predict_final(fit, ext[c("patient_id", "gene_a")], "patient_id")
  endpoint <- sb_surv_endpoint(ext, "patient_id", "time", "event", "month")
  evaluated <- sb_evaluate_external(pred, endpoint)
  expect_identical(sb_result_table(evaluated, metric = "time_auc")$status, "not_requested")
  file <- tempfile(fileext = ".html")
  sb_render_report(list(CV = cv, External = evaluated), file)
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  for (metric in c("time_auc", "uno_c", "ipcw_brier")) {
    expect_identical(grepl(paste0("<summary>", metric, "</summary>"), html, fixed = TRUE), FALSE)
  }
})

# ---- survival-paper-figures ----

paper_survival_task <- function() {
  withr::local_seed(32)
  d <- as.data.frame(matrix(rnorm(180 * 6), 180, 6))
  names(d) <- paste0("g", 1:6)
  d$id <- paste0("p", seq_len(nrow(d)))
  d$time <- rexp(nrow(d), exp(d$g1 - d$g2 + .7 * d$g3))
  d$event <- rbinom(nrow(d), 1, .8)
  sb_surv_task(d, "id", "time", "event", paste0("g", 1:6), "month")
}

test_that("C-index ranking reuses complete-cohort means without AUC labels", {
  d <- data.frame(model = rep(c("A", "B", "C"), each = 2), cohort = rep(c("V1", "V2"), 3),
    c_index = c(.8, .7, .9, .5, .99, NA))
  z <- sb_plot_model_ranking(d, c("V1", "V2"), metric = "c_index")
  expect_equal(z$ranking$mean_c_index, c(.75, .7, NA))
  expect_identical(z$ranking$status, c("ranked", "ranked", "incomplete"))
  expect_identical(z$heatmap$labels$fill, "C-index")
  expect_identical(z$bar$labels$x, "Mean C-index")
  expect_s3_class(ggplot2::ggplot_build(z$bar), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(z$heatmap), "ggplot_built")
})

test_that("survival Enet matches native coefficients and keeps patient groups", {
  task <- paper_survival_task()
  x <- task$data[task$predictors]
  groups <- rep(1:90, each = 2)
  adapter <- sb_resolve_selector_adapter("enet_cox_nonzero", list(inner_folds = 3L, seed = 1L))
  s <- adapter$fit(x, task$data$time, task$data$event, group_ids = groups)
  native <- glmnet::glmnet(as.matrix(x), survival::Surv(task$data$time, task$data$event),
    family = "cox", cox.ties = "breslow", alpha = .5, standardize = TRUE)
  expect_equal(s$scores$coefficient, as.numeric(stats::coef(native, s = s$lambda)))
  expect_equal(s$alpha, .5)
  expect_identical(all(vapply(split(s$fold_id, groups), function(g) length(unique(g)) == 1L, logical(1))), TRUE)
  expect_identical(surv_selector("Enet"), surv_selector("ElasticNet"))
  for (p in list(surv_pipeline("Cox", surv_selector("Enet")),
    surv_pipeline("ElasticNet", surv_selector("Enet")),
    surv_pipeline("Cox", surv_selector(c("Enet", "StepCox"))))) {
    expect_identical(sb_lint_pipeline(p)$status, "admitted")
    fit <- sb_fit_final(task, p)
    expect_identical(fit$status, "success")
    expect_identical(sb_run_cv(task, p, folds = 3)$status, "success")
  }
})

test_that("risk panels share order and never recalibrate on validation outcomes", {
  task <- paper_survival_task()
  fit <- sb_fit_final(task, surv_pipeline("Cox", surv_selector("Enet")))
  ref <- surv_risk_reference(fit, task)
  z <- surv_risk_plots(ref, task)
  sm <- summary(z$km, censored = TRUE)
  expect_equal(z$plots$km$data$lower[z$plots$km$data$time > 0], sm$lower)
  expect_equal(z$plots$km$data$upper[z$plots$km$data$time > 0], sm$upper)
  expect_equal(z$cox_group$log_hr, unname(stats::coef(survival::coxph(
    survival::Surv(time, event) ~ group, z$patients))))
  scores <- sb_risk_predict(fit, task)
  expect_equal(ref$cutoff, stats::median(scores))
  expect_identical(z$patients$patient_id, task$data$id[order(scores)])
  expect_identical(rownames(z$expression), z$patients$patient_id)
  expected <- scale(as.matrix(task$data[order(scores), ref$genes, drop = FALSE]), center = ref$center, scale = ref$scale)
  expect_equal(unname(z$expression), unname(expected), ignore_attr = TRUE)
  test <- survival::survdiff(survival::Surv(time, event) ~ group, data = z$patients)
  expect_equal(z$logrank_p, stats::pchisq(test$chisq, 1, lower.tail = FALSE))
  expect_equal(z$at_risk$n_risk, mapply(function(t, g) sum(z$patients$time >= t & z$patients$group == g),
    z$at_risk$time, z$at_risk$group))
  altered <- task$data
  altered$time <- rev(altered$time)
  altered$event <- rev(altered$event)
  validation <- sb_surv_task(altered, "id", "time", "event", task$predictors, "month")
  changed <- surv_risk_plots(ref, validation)
  expect_identical(changed$patients[c("patient_id", "risk", "group")], z$patients[c("patient_id", "risk", "group")])
  expect_identical(changed$expression, z$expression)
  for (p in z$plots) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  single <- ref; single$cutoff <- Inf
  expect_identical(surv_risk_plots(single, task)$status, "single_group")
})

test_that("survival ROC handles ties, censoring, ID order and unavailable times", {
  d <- data.frame(id = letters[1:5], time = c(1, 2, 3, 5, 6), event = c(1, 0, 1, 1, 0), x = 1:5)
  task <- sb_surv_task(d, "id", "time", "event", "x", "month")
  risk <- c(3, 0, 1, 1, 2)
  z <- surv_plot_roc(task, risk, c(.5, 4, 7))
  expect_identical(z$statistics$status, c("not_evaluable", "success", "not_evaluable"))
  expect_equal(z$statistics$n_early_censored[2], 1)
  # At t=4 case weights are 1 and 4/3; the second case earns 1/4 pair credit.
  expect_equal(z$statistics$auc[2], (1 + (4/3) * .25) / (1 + 4/3))
  shuffled <- setNames(rev(risk), rev(d$id))
  expect_equal(surv_plot_roc(task, shuffled, 4)$statistics$auc, z$statistics$auc[2])
  tied <- surv_plot_roc(task, rep(1, 5), 4)
  expect_equal(tied$statistics$auc, .5)
  q <- sb_time_auc_stratum(d$time, d$event, plogis(-risk), rep(TRUE, 5),
    sb_fit_censoring_km(d$time, d$event), 4)
  expect_equal(z$statistics$auc[2], q$estimate)
  expect_s3_class(ggplot2::ggplot_build(z$plot), "ggplot_built")
  expect_error(surv_plot_roc(task, risk, c(4, 4)), "unique positive")
})
})
