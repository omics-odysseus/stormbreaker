# ---- model-glmnet-tuning ----
local({
# ---- algorithm-glmnet-tuning ----

test_that("all six glmnet task/mixing variants expose unchanged saved tuning", {
  skip_if_not_installed("glmnet")
  for (kind in c("binary", "survival")) {
    d <- if (kind == "binary") make_sb_binary_fixture(n = 90L, seed = 108L) else
      make_sb_fixture(n = 90L, seed = 8612L)
    task <- if (kind == "binary") sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")) else
      sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
    for (alpha in c(0, .5, 1)) {
      learner <- if (kind == "binary") "glmnet_logistic" else "glmnet_cox"
      fit <- sb_fit_final(task, sb_pipeline(learner = learner,
        inner_folds = 3L, inner_seed = 5L, glmnet_alpha = alpha))
      expect_identical(fit$status, "success")
      expect_null(fit$diagnostic_artifacts)
      before <- fit
      for (type in c("cv_curve", "nonzero_path")) {
        q <- sb_algorithm_diagnostics(fit, type)
        expect_identical(q$table[names(fit$tuning$curve)], fit$tuning$curve)
        expect_equal(q$table$lambda[q$table$selected], fit$tuning$lambda)
        expect_identical(q$alpha, alpha)
        if (requireNamespace("ggplot2", quietly = TRUE)) {
          expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
        }
      }
      expect_identical(fit, before)
    }
  }
})

test_that("tuning diagnostic rejects missing or corrupt native uncertainty", {
  skip_if_not_installed("glmnet")
  d <- make_sb_binary_fixture(n = 90L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "glmnet_logistic", inner_folds = 3L, inner_seed = 5L))
  fit$tuning$curve$standard_error[1] <- -1
  err <- tryCatch(sb_algorithm_diagnostics(fit, "cv_curve"), error = identity)
  expect_s3_class(err, "stormbreaker_invalid_diagnostic")
  fit$tuning$curve <- NULL
  err <- tryCatch(sb_algorithm_diagnostics(fit, "nonzero_path"), error = identity)
  expect_s3_class(err, "stormbreaker_invalid_diagnostic")
})

# ---- glmnet-active-set-stability ----

test_that("glmnet active sets are summarized from outer-fold coefficients", {
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 1,
    inner_folds = 3L,
    inner_seed = 1L
  )
  coefficients <- list(
    list(
      fold = 2L,
      status = "success",
      nonzero_features = c("gene_b", "gene_c")
    ),
    list(
      fold = 1L,
      status = "success",
      nonzero_features = c("gene_a", "gene_b")
    ),
    list(
      fold = 3L,
      status = "success",
      nonzero_features = "gene_b"
    )
  )

  observed <- sb_glmnet_active_set_stability(
    pipeline,
    coefficients,
    c("gene_a", "gene_b", "gene_c")
  )

  expect_identical(observed$learner, "glmnet_logistic")
  expect_identical(observed$alpha, 1)
  expect_identical(observed$status, "success")
  expect_identical(observed$per_fold$fold, 1:3)
  expect_identical(observed$per_fold$n_selected, c(2L, 2L, 1L))
  expect_equal(observed$mean_pairwise_jaccard, 4 / 9, tolerance = 0)
  expect_identical(
    observed$feature_frequency$selected_count,
    c(1L, 3L, 1L)
  )
  expect_equal(
    observed$feature_frequency$selection_frequency,
    c(1 / 3, 1, 1 / 3),
    tolerance = 0
  )
})

test_that("ridge and non-glmnet fits have no sparse active-set stability", {
  coefficients <- list(
    list(
      fold = 1L,
      status = "success",
      nonzero_features = c("gene_a", "gene_b")
    ),
    list(
      fold = 2L,
      status = "success",
      nonzero_features = c("gene_a", "gene_b")
    )
  )
  ridge <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0,
    inner_folds = 3L,
    inner_seed = 1L
  )
  ordinary <- sb_pipeline()

  ridge_audit <- sb_glmnet_active_set_stability(
    ridge,
    coefficients,
    c("gene_a", "gene_b")
  )
  ordinary_audit <- sb_glmnet_active_set_stability(
    ordinary,
    coefficients,
    c("gene_a", "gene_b")
  )

  expect_identical(ridge_audit$status, "not_applicable")
  expect_identical(ridge_audit$alpha, 0)
  expect_identical(ridge_audit$n_successful_folds, 0L)
  expect_true(all(is.na(
    ridge_audit$feature_frequency$selection_frequency
  )))
  expect_identical(ordinary_audit$status, "not_applicable")
  expect_true(is.na(ordinary_audit$alpha))
})

test_that("empty successful glmnet active sets remain visibly unevaluable", {
  pipeline <- sb_pipeline(
    selector = "fixed_features",
    learner = "glmnet_logistic",
    glmnet_alpha = 1,
    inner_folds = 3L,
    inner_seed = 1L,
    selector_features = c("gene_b", "gene_a"),
    feature_set_id = "pair"
  )
  coefficients <- lapply(1:3, function(fold) {
    list(
      fold = fold,
      status = "success",
      nonzero_features = character()
    )
  })

  observed <- sb_glmnet_active_set_stability(
    pipeline,
    coefficients,
    c("gene_a", "gene_b", "gene_c")
  )

  expect_identical(observed$status, "not_evaluable")
  expect_identical(
    observed$feature_frequency$feature,
    c("gene_b", "gene_a")
  )
  expect_identical(observed$per_fold$n_selected, rep(0L, 3L))
  expect_true(is.na(observed$mean_pairwise_jaccard))
})

test_that("binary and survival CV expose glmnet active-set audits", {
  skip_if_not_installed("glmnet")

  binary_data <- make_sb_binary_fixture(n = 60L, seed = 2404L)
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  binary <- sb_run_cv(
    binary_task,
    sb_pipeline(
      learner = "glmnet_logistic",
      glmnet_alpha = 1,
      inner_folds = 3L,
      inner_seed = 2405L
    ),
    folds = 3L,
    seed = 2406L
  )

  survival_data <- make_sb_fixture(n = 60L, seed = 2407L)
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  survival <- sb_run_cv(
    survival_task,
    sb_pipeline(
      learner = "glmnet_cox",
      glmnet_alpha = 1,
      inner_folds = 3L,
      inner_seed = 2408L
    ),
    folds = 3L,
    seed = 2409L
  )

  for (result in list(binary, survival)) {
    audit <- result$glmnet_active_set_stability
    expect_identical(audit$learner, result$pipeline$learner)
    expect_identical(audit$alpha, 1)
    expect_identical(
      audit$per_fold$n_selected,
      vapply(result$coefficients, function(x) x$n_nonzero, integer(1))
    )
    expect_identical(
      audit$n_successful_folds,
      as.integer(sum(vapply(
        result$coefficients,
        function(x) identical(x$status, "success"),
        logical(1)
      )))
    )
  }
})

# ---- glmnet-cox-tuning ----

test_that("inner CV selects lambda reproducibly from its own patients", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 80, seed = 20260813)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  first <- sb_select_glmnet_cox_lambda(
    x, d$time, d$event, inner_folds = 5L, seed = 1701L
  )
  repeated <- sb_select_glmnet_cox_lambda(
    x, d$time, d$event, inner_folds = 5L, seed = 1701L
  )

  expect_identical(first$fold_id, repeated$fold_id)
  expect_equal(first$lambda, repeated$lambda, tolerance = 0)
  expect_equal(
    first$lambda,
    first$curve$lambda[which.min(first$curve$mean_deviance)],
    tolerance = 0
  )
  expect_identical(first$selection_rule, "lambda.min")

  fold_has_both_outcomes <- vapply(seq_len(5L), function(fold) {
    setequal(unique(d$event[first$fold_id == fold]), c(0L, 1L))
  }, logical(1))
  expect_identical(fold_has_both_outcomes, rep(TRUE, 5L))
})

test_that("Cox inner tuning passes the declared standardization setting", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 80, seed = 1)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  tuned <- sb_select_glmnet_cox_lambda(
    x,
    d$time,
    d$event,
    inner_folds = 3L,
    seed = 1701L,
    standardize = FALSE,
    selection_rule = "lambda.1se",
    measure = "C"
  )
  manual <- glmnet::cv.glmnet(
    x = as.matrix(x),
    y = survival::Surv(d$time, d$event),
    family = "cox",
    cox.ties = "breslow",
    alpha = 1,
    foldid = tuned$fold_id,
    type.measure = "C",
    standardize = FALSE
  )

  expect_identical(tuned$standardize, FALSE)
  expect_identical(tuned$selection_rule, "lambda.1se")
  expect_identical(tuned$measure, "C")
  expect_equal(tuned$lambda, unname(manual$lambda.1se), tolerance = 0)
  expect_gt(tuned$lambda, unname(manual$lambda.min))
  expect_equal(tuned$curve$lambda, manual$lambda, tolerance = 0)
  expect_equal(tuned$curve$mean_measure, manual$cvm, tolerance = 0)
})

test_that("outer assessment changes cannot alter inner lambda selection", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 90, seed = 20260813)
  analysis_rows <- 1:70
  assessment_rows <- 71:90
  predictors <- c("gene_a", "gene_b", "gene_c")
  select_from_analysis <- function(data) {
    sb_select_glmnet_cox_lambda(
      data[analysis_rows, predictors, drop = FALSE],
      data$time[analysis_rows],
      data$event[analysis_rows],
      inner_folds = 5L,
      seed = 1701L
    )
  }
  original <- select_from_analysis(d)

  changed <- d
  changed$time[assessment_rows] <- changed$time[assessment_rows] + 10000
  changed$event[assessment_rows] <- 1L - changed$event[assessment_rows]
  changed$gene_a[assessment_rows] <- seq_along(assessment_rows) * 10000
  perturbed <- select_from_analysis(changed)

  expect_equal(original$lambda, perturbed$lambda, tolerance = 0)
  expect_identical(original$fold_id, perturbed$fold_id)
  expect_identical(original$curve, perturbed$curve)
})

test_that("inner CV reports missing glmnet and infeasible folds", {
  d <- make_sb_fixture(n = 20, seed = 20260813)
  x <- d[c("gene_a", "gene_b")]
  local_mocked_bindings(sb_glmnet_available = function() FALSE)

  missing_error <- capture_sb_error(sb_select_glmnet_cox_lambda(
    x, d$time, d$event, inner_folds = 5L, seed = 1701L
  ))
  expect_s3_class(missing_error, "stormbreaker_missing_glmnet")

  local_mocked_bindings(sb_glmnet_available = function() TRUE)
  small <- make_sb_fixture(n = 8, seed = 20260813)
  fold_error <- capture_sb_error(sb_select_glmnet_cox_lambda(
    small[c("gene_a", "gene_b")],
    small$time,
    small$event,
    inner_folds = 5L,
    seed = 1701L
  ))
  expect_s3_class(fold_error, "stormbreaker_infeasible_inner_folds")
})

test_that("inner CV records a glmnet fitting error explicitly", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 20, seed = 20260813)
  x <- data.frame(constant_a = 1, constant_b = 2)[rep(1L, 20L), ]

  error <- capture_sb_error(sb_select_glmnet_cox_lambda(
    x, d$time, d$event, inner_folds = 5L, seed = 1701L
  ))

  expect_s3_class(error, "stormbreaker_glmnet_tuning_failure")
  expect_s3_class(error$parent, "simpleError")
  expect_match(conditionMessage(error), "zero variance")
})

test_that("inner CV preserves its dedicated warning failure class", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 20, seed = 20260813)
  original_cv_glmnet <- glmnet::cv.glmnet
  local_mocked_bindings(
    cv.glmnet = function(...) {
      warning("synthetic cv warning from glmnet", call. = FALSE)
      original_cv_glmnet(...)
    },
    .package = "glmnet"
  )

  error <- capture_sb_error(sb_select_glmnet_cox_lambda(
    d[c("gene_a", "gene_b")],
    d$time,
    d$event,
    inner_folds = 5L,
    seed = 1701L
  ))

  expect_s3_class(error, "stormbreaker_glmnet_tuning_warning")
  expect_s3_class(error$parent, "simpleWarning")
  expect_match(conditionMessage(error), "synthetic cv warning from glmnet")
})

test_that("a validated convergence-tail warning remains visible", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260813)
  tail_warning <- paste0(
    "from glmnet C++ code (error code -27); Convergence for 27th ",
    "lambda value not reached after maxit=100000 iterations; ",
    "solutions for larger lambdas returned"
  )
  original_cv_glmnet <- glmnet::cv.glmnet
  original_glmnet <- glmnet::glmnet
  observed_cv_ties <- character()
  observed_path_ties <- character()
  local_mocked_bindings(
    cv.glmnet = function(..., cox.ties = NULL) {
      observed_cv_ties <<- c(observed_cv_ties, cox.ties)
      warning(tail_warning, call. = FALSE)
      original_cv_glmnet(..., cox.ties = cox.ties)
    },
    glmnet = function(..., cox.ties = NULL) {
      observed_path_ties <<- c(observed_path_ties, cox.ties)
      original_glmnet(..., cox.ties = cox.ties)
    },
    .package = "glmnet"
  )

  tuned <- sb_select_glmnet_cox_lambda(
    d[c("gene_a", "gene_b", "gene_c")],
    d$time,
    d$event,
    inner_folds = 5L,
    seed = 1701L
  )
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = tuned$alpha,
    lambda = tuned$lambda,
    standardize = tuned$standardize
  )
  learner$fit(
    d[c("gene_a", "gene_b", "gene_c")],
    d$time,
    d$event
  )

  expect_identical(tuned$warnings, tail_warning)
  expect_identical(observed_cv_ties, "breslow")
  expect_gte(length(observed_path_ties), 6L)
  expect_identical(unique(observed_path_ties), "breslow")
  expect_identical(nrow(tuned$warning_validation), 5L)
  expect_identical(tuned$warning_validation$fold, seq_len(5L))
  expect_true(all(tuned$lambda >= tuned$warning_validation$min_lambda))
  expect_true(all(tuned$warning_validation$assessment_risk_finite))
})

test_that("Cox refits accept only usable, recorded convergence tails", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 80, seed = 20260813)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  original_glmnet <- glmnet::glmnet
  native <- original_glmnet(
    as.matrix(x), survival::Surv(d$time, d$event),
    family = "cox", cox.ties = "breslow"
  )
  tail_warning <- paste0(
    "from glmnet C++ code (error code -88); Convergence for 88th ",
    "lambda value not reached after maxit=100000 iterations; ",
    "solutions for larger lambdas returned"
  )
  local_mocked_bindings(
    glmnet = function(...) {
      warning(tail_warning, call. = FALSE)
      original_glmnet(...)
    },
    .package = "glmnet"
  )
  learner <- sb_glmnet_cox_learner_adapter(1, native$lambda[[5L]])
  fitted <- learner$fit(x, d$time, d$event)
  expect_identical(attr(fitted, "stormbreaker_fit_warnings"), tail_warning)
  expect_equal(
    learner$predict(fitted, x),
    unname(drop(predict(
      native, as.matrix(x), s = native$lambda[[5L]], type = "link"
    ))),
    tolerance = 0
  )

  outside <- sb_glmnet_cox_learner_adapter(1, min(native$lambda) / 2)
  error <- capture_sb_error(outside$fit(x, d$time, d$event))
  expect_s3_class(error, "stormbreaker_glmnet_lambda_outside_path")

  local_mocked_bindings(
    glmnet = function(...) warning("unrecognized fit warning", call. = FALSE),
    .package = "glmnet"
  )
  error <- capture_sb_error(learner$fit(x, d$time, d$event))
  expect_s3_class(error, "stormbreaker_glmnet_fit_warning")
})

test_that("Cox refit warnings reach CV and final-fit results", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 90, seed = 20260813)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_cox", inner_folds = 3L, inner_seed = 1701L
  )
  original_glmnet <- glmnet::glmnet
  tail_warning <- paste0(
    "from glmnet C++ code (error code -88); Convergence for 88th ",
    "lambda value not reached after maxit=100000 iterations; ",
    "solutions for larger lambdas returned"
  )
  local_mocked_bindings(
    glmnet = function(...) {
      warning(tail_warning, call. = FALSE)
      original_glmnet(...)
    },
    .package = "glmnet"
  )
  cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 1701L)
  final <- sb_fit_final(task, pipeline)
  expect_identical(cv$status, "success")
  expect_identical(final$status, "success")
  expect_identical(
    vapply(
      cv$tuning,
      function(z) any(grepl("Learner refit:", z$warnings)),
      logical(1)
    ),
    rep(TRUE, 3L)
  )
  expect_match(tail(final$tuning$warnings, 1L), "Learner refit:")
})

test_that("an inner-selected lambda refits and predicts held-out patients", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 90, seed = 20260813)
  analysis_rows <- 1:70
  assessment_rows <- 71:90
  predictors <- c("gene_a", "gene_b", "gene_c")
  tuned <- sb_select_glmnet_cox_lambda(
    d[analysis_rows, predictors, drop = FALSE],
    d$time[analysis_rows],
    d$event[analysis_rows],
    inner_folds = 5L,
    seed = 1701L
  )
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = tuned$alpha,
    lambda = tuned$lambda,
    standardize = tuned$standardize
  )
  model <- learner$fit(
    d[analysis_rows, predictors, drop = FALSE],
    d$time[analysis_rows],
    d$event[analysis_rows]
  )
  risk <- learner$predict(
    model,
    d[assessment_rows, predictors, drop = FALSE]
  )

  expect_length(risk, length(assessment_rows))
  expect_identical(all(is.finite(risk)), TRUE)
  expect_gt(stats::sd(risk), 0)
})

test_that("glmnet Cox preserves native predictions at rounded path endpoints", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 80, seed = 20260813)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  native <- glmnet::glmnet(as.matrix(x), survival::Surv(d$time, d$event),
    family = "cox", cox.ties = "breslow")
  endpoints <- range(native$lambda)
  near <- endpoints * c(1 - 2 * .Machine$double.eps,
    1 + 2 * .Machine$double.eps)
  for (lambda in near) {
    learner <- sb_glmnet_cox_learner_adapter(1, lambda)
    model <- learner$fit(x, d$time, d$event)
    expect_identical(model$lambda, native$lambda)
    expect_identical(attr(model, "stormbreaker_lambda"), lambda)
    expect_equal(learner$predict(model, x),
      as.numeric(predict(native, as.matrix(x), s = lambda, type = "link")),
      tolerance = 0)
  }
  for (lambda in endpoints * c(1 - 1e-6, 1 + 1e-6)) {
    learner <- sb_glmnet_cox_learner_adapter(1, lambda)
    error <- capture_sb_error(learner$fit(x, d$time, d$event))
    expect_s3_class(error, "stormbreaker_glmnet_lambda_outside_path")
  }
})

test_that("Cox convergence validation distinguishes endpoint roundoff from a missing solution", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 80, seed = 20260813)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  folds <- rep(1:4, length.out = nrow(d))
  minima <- vapply(1:4, function(f) {
    ix <- folds != f
    model <- glmnet::glmnet(as.matrix(x[ix, ]),
      survival::Surv(d$time[ix], d$event[ix]), family = "cox", cox.ties = "breslow")
    min(model$lambda)
  }, numeric(1))
  lambda <- max(minima) * (1 - 2 * .Machine$double.eps)
  validation <- sb_validate_glmnet_selected_lambda(x, d$time, d$event,
    folds, lambda, alpha = 1, standardize = TRUE)
  expect_equal(nrow(validation), 4L)
  expect_true(all(validation$assessment_risk_finite))
  error <- capture_sb_error(sb_validate_glmnet_selected_lambda(
    x, d$time, d$event, folds, max(minima) * (1 - 1e-6),
    alpha = 1, standardize = TRUE))
  expect_s3_class(error, "stormbreaker_glmnet_selected_lambda_not_converged")
})

# ---- glmnet-fold-local-selection ----

test_that("grouped glmnet keeps row identities and patients at every fit level", {
  skip_if_not_installed("glmnet", minimum_version = "4.1.2")
  predictors <- c("gene_a", "gene_b", "gene_c")
  for (type in c("binary", "survival")) {
    d <- if (type == "binary") sb_example_binary(120L, 3750L) else {
      sb_example_survival(120L, 3751L)
    }
    d$group <- rep(paste0("patient", 1:60), each = 2)
    if (type == "binary") d$outcome <- rep(c(0L, 1L), 60)
    task <- if (type == "binary") {
      sb_binary_task(d, "patient_id", "outcome", predictors, group_col = "group")
    } else {
      sb_surv_task(d, "patient_id", "time", "event", predictors, "month",
                   group_col = "group")
    }
    expect_identical(sb_subset_task(task, 1:100)$group_col, "group")
    for (selected in c(FALSE, TRUE)) {
      selector <- if (!selected) "none" else if (type == "binary") {
        "univ_binary_auc_top_k"
      } else "univ_cox_top_k"
      pipeline <- if (selected) {
        sb_pipeline(selector = selector, selector_n_features = 2L,
          learner = if (type == "binary") "glmnet_logistic" else "glmnet_cox",
          inner_folds = 3L, inner_seed = 3752L)
      } else sb_pipeline(
        learner = if (type == "binary") "glmnet_logistic" else "glmnet_cox",
        inner_folds = 3L, inner_seed = 3752L)
      cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 3753L)
      final <- sb_fit_final(task, pipeline)
      expect_identical(cv$status, "success")
      expect_identical(cv$complete_oof, TRUE)
      expect_identical(final$status, "success")
      expect_identical(final$development$group_col, "group")
      expect_identical(unname(final$development$group_ids), d$group)
      expect_identical(cv$predictions$patient_id, d$patient_id)
      outer <- cv$provenance$fold_assignment
      for (i in seq_along(cv$tuning)) {
        expect_identical(sort(names(cv$tuning[[i]]$fold_id)),
          sort(d$patient_id[outer$fold != i]))
      }
      audits <- c(cv$tuning, list(final$tuning))
      for (audit in audits) {
        rows <- match(names(audit$fold_id), d$patient_id)
        expect_identical(anyNA(rows), FALSE)
        expect_equal(anyDuplicated(names(audit$fold_id)), 0L)
        by_group <- split(unname(audit$fold_id), d$group[rows])
        expect_identical(unname(vapply(by_group, function(x) {
          length(unique(x))
        }, integer(1))), rep(1L, length(by_group)))
        expect_length(audit$inner_selections, if (selected) 3L else 0L)
        if (selected) for (inner in audit$inner_selections) {
          train <- rows[audit$fold_id != inner$fold]
          assess <- rows[audit$fold_id == inner$fold]
          expect_identical(inner$analysis_patient_ids, d$patient_id[train])
          expect_length(intersect(d$group[train], d$group[assess]), 0L)
          expected <- if (type == "binary") {
            sb_fit_univ_binary_auc_top_k(d[train, predictors], d$outcome[train], 2L)
          } else sb_fit_univ_cox_top_k(
            d[train, predictors], d$time[train], d$event[train], 2L)
          expect_identical(inner$selected_features, expected$selected_features)
        }
      }
      if (!selected) {
        args <- list(x = as.matrix(d[predictors]),
          y = if (type == "binary") d$outcome else survival::Surv(d$time, d$event),
          family = if (type == "binary") "binomial" else "cox",
          alpha = pipeline$learner_params$alpha,
          standardize = pipeline$learner_params$standardize,
          foldid = unname(final$tuning$fold_id), type.measure = "deviance")
        if (type == "survival") args$cox.ties <- "breslow"
        native <- do.call(glmnet::cv.glmnet, args)
        expect_identical(final$tuning$lambda, unname(native$lambda.min))
        expect_identical(final$tuning$curve$mean_deviance, native$cvm)
      }
    }
  }
})

test_that("Top-K glmnet combinations are admitted with inner CV", {
  survival_pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "glmnet_cox",
    selector_n_features = 2L,
    inner_folds = 3L,
    inner_seed = 821L
  )
  binary_pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "glmnet_logistic",
    selector_n_features = 2L,
    inner_folds = 3L,
    inner_seed = 822L
  )

  expect_identical(survival_pipeline$task_type, "survival")
  expect_identical(binary_pipeline$task_type, "binary")
  expect_identical(survival_pipeline$learner_params$inner_folds, 3L)
  expect_identical(binary_pipeline$learner_params$inner_folds, 3L)
  expect_identical(
    sb_lint_pipeline(survival_pipeline)$status,
    "admitted"
  )
  expect_identical(
    sb_lint_pipeline(binary_pipeline)$status,
    "admitted"
  )
})

test_that("Cox Top-K is refitted inside every glmnet inner analysis fold", {
  skip_if_not_installed("glmnet", minimum_version = "4.1.2")

  d <- sb_example_survival(n = 90L, seed = 20260822L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  x <- d[predictors]
  patient_ids <- d$patient_id
  fold_id <- sb_make_stratified_fold_id(d$event, 3L, 823L)
  selector <- sb_univ_cox_top_k_selector_adapter(2L)
  fit_outer <- function(data) {
    sb_fit_runtime_selector(
      selector,
      data[predictors],
      list(
        task_type = "survival",
        time = data$time,
        event = data$event
      )
    )
  }
  tune <- function(data) {
    sb_select_glmnet_cox_lambda(
      data[predictors],
      data$time,
      data$event,
      inner_folds = 3L,
      seed = 823L,
      selector = selector,
      fitted_selector = fit_outer(data),
      patient_ids = patient_ids
    )
  }
  local_mocked_bindings(
    sb_make_stratified_fold_id = function(...) fold_id
  )
  original_cv_glmnet <- glmnet::cv.glmnet
  observed_cox_ties <- character()
  local_mocked_bindings(
    cv.glmnet = function(..., cox.ties = NULL) {
      observed_cox_ties <<- c(observed_cox_ties, cox.ties)
      original_cv_glmnet(..., cox.ties = cox.ties)
    },
    .package = "glmnet"
  )

  original <- tune(d)
  expect_identical(original$inner_selection_method, "glmnet_exclude_function")
  expect_identical(original$alignment, "lambda")
  expect_length(original$inner_selections, 3L)
  for (fold in seq_len(3L)) {
    direct <- sb_fit_runtime_selector(
      selector,
      x[fold_id != fold, , drop = FALSE],
      list(
        task_type = "survival",
        time = d$time[fold_id != fold],
        event = d$event[fold_id != fold]
      )
    )
    audit <- original$inner_selections[[fold]]
    expect_identical(audit$fold, fold)
    expect_identical(
      audit$analysis_patient_ids,
      patient_ids[fold_id != fold]
    )
    expect_identical(audit$selected_features, direct$selected_features)
  }

  changed <- d
  assessment_rows <- fold_id == 1L
  changed$time[assessment_rows] <- changed$time[assessment_rows] + 1
  changed$event[assessment_rows] <- 1L - changed$event[assessment_rows]
  changed[assessment_rows, predictors] <-
    changed[assessment_rows, predictors] + 100
  perturbed <- tune(changed)
  expect_identical(
    perturbed$inner_selections[[1L]],
    original$inner_selections[[1L]]
  )
  expect_identical(observed_cox_ties, rep("breslow", 2L))
})

test_that("Top-K glmnet pipelines run and retain compact inner audits", {
  skip_if_not_installed("glmnet", minimum_version = "4.1.2")

  survival_data <- sb_example_survival(n = 90L, seed = 824L)
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  survival_pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "glmnet_cox",
    selector_n_features = 2L,
    inner_folds = 3L,
    inner_seed = 825L
  )
  survival_cv <- sb_run_cv(
    survival_task,
    survival_pipeline,
    folds = 3L,
    seed = 826L
  )
  survival_fit <- sb_fit_final(survival_task, survival_pipeline)

  expect_identical(survival_cv$status, "success")
  expect_true(survival_cv$complete_oof)
  expect_true(all(vapply(survival_cv$tuning, function(audit) {
    identical(audit$status, "success") &&
      identical(audit$inner_selection_method, "glmnet_exclude_function") &&
      length(audit$inner_selections) == 3L
  }, logical(1))))
  expect_identical(survival_fit$status, "success")
  expect_length(survival_fit$tuning$inner_selections, 3L)
  expect_length(survival_fit$required_features, 2L)

  binary_data <- sb_example_binary(n = 90L, seed = 827L)
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  binary_pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "glmnet_logistic",
    selector_n_features = 2L,
    inner_folds = 3L,
    inner_seed = 828L
  )
  binary_cv <- sb_run_cv(
    binary_task,
    binary_pipeline,
    folds = 3L,
    seed = 829L
  )
  binary_fit <- sb_fit_final(binary_task, binary_pipeline)
  prediction <- sb_predict_final(
    binary_fit,
    binary_data[1:4, c("patient_id", "gene_c", "gene_a", "gene_b")],
    "patient_id"
  )

  expect_identical(binary_cv$status, "success")
  expect_true(binary_cv$complete_oof)
  expect_true(all(vapply(binary_cv$tuning, function(audit) {
    identical(audit$status, "success") &&
      identical(audit$inner_selection_method, "glmnet_exclude_function") &&
      length(audit$inner_selections) == 3L
  }, logical(1))))
  expect_identical(binary_fit$status, "success")
  expect_length(binary_fit$tuning$inner_selections, 3L)
  expect_length(binary_fit$required_features, 2L)
  expect_identical(prediction$status, "success")
  expect_true(all(is.finite(prediction$predictions$probability)))
})

test_that("an inner Top-K failure remains an explicit tuning failure", {
  skip_if_not_installed("glmnet", minimum_version = "4.1.2")

  n <- 30L
  outcome <- rep(c(0L, 1L), length.out = n)
  fold_id <- rep(seq_len(3L), length.out = n)
  x <- data.frame(
    signal = seq_len(n),
    fold_one_only = ifelse(fold_id == 1L, seq_len(n), 0)
  )
  patient_ids <- sprintf("P%03d", seq_len(n))
  selector <- sb_univ_binary_auc_top_k_selector_adapter(2L)
  outer <- sb_fit_runtime_selector(
    selector,
    x,
    list(task_type = "binary", outcome = outcome)
  )
  local_mocked_bindings(
    sb_make_stratified_fold_id = function(...) fold_id
  )

  error <- capture_sb_error(sb_select_glmnet_logistic_lambda(
    x,
    outcome,
    inner_folds = 3L,
    seed = 830L,
    selector = selector,
    fitted_selector = outer,
    patient_ids = patient_ids
  ))

  expect_s3_class(error, "stormbreaker_insufficient_valid_features")
  expect_match(conditionMessage(error), "only 1 features")
})

test_that("outer survival selection stays visible when inner tuning fails", {
  d <- sb_example_survival(n = 60L, seed = 831L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "glmnet_cox",
    selector_n_features = 2L,
    inner_folds = 3L,
    inner_seed = 832L
  )
  local_mocked_bindings(
    sb_select_glmnet_cox_lambda = function(...) {
      sb_abort(
        "synthetic inner-selection failure",
        "stormbreaker_synthetic_inner_failure"
      )
    }
  )

  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 833L)

  expect_identical(result$status, "failed")
  expect_true(all(vapply(result$selections, function(audit) {
    identical(audit$status, "success") &&
      length(audit$selected_features) == 2L
  }, logical(1))))
  expect_true(all(vapply(result$tuning, function(audit) {
    identical(audit$status, "failed") &&
      identical(
        audit$error_class,
        "stormbreaker_synthetic_inner_failure"
      )
  }, logical(1))))
})

# ---- glmnet-logistic-adapter ----

test_that("glmnet logistic can use a converged prefix of a truncated path", {
  skip_if_not_installed("glmnet")
  x <- withr::with_seed(5106, as.data.frame(matrix(stats::rnorm(60 * 6), 60, 6)))
  y <- rep(0:1, 30L)
  native <- glmnet::glmnet
  messages <- character()
  reference <- withCallingHandlers(
    native(as.matrix(x), y, family = "binomial", maxit = 20L),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_gt(length(messages), 0L)
  expect_gte(length(reference$lambda), 3L)
  lambda <- reference$lambda[[3L]]
  local_mocked_bindings(glmnet = function(...) {
    args <- list(...)
    args$maxit <- 20L
    do.call(native, args)
  }, .package = "glmnet")
  adapter <- sb_glmnet_logistic_learner_adapter(alpha = 1, lambda = lambda)
  fitted <- adapter$fit(x, y)
  expect_identical(attr(fitted, "stormbreaker_fit_warnings"), messages)
  expect_identical(adapter$predict(fitted, x), unname(drop(stats::predict(
    reference, newx = as.matrix(x), s = lambda, type = "response"))))
  outside <- sb_glmnet_logistic_learner_adapter(
    alpha = 1, lambda = min(reference$lambda) / 10)
  expect_s3_class(tryCatch(outside$fit(x, y), error = identity),
    "stormbreaker_glmnet_lambda_outside_path")
})

test_that("the glmnet logistic adapter matches a direct path fit", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 60)
  analysis_rows <- 1:45
  assessment_rows <- 46:60
  predictors <- c("gene_a", "gene_b", "gene_c")
  x_analysis <- d[analysis_rows, predictors, drop = FALSE]
  x_assessment <- d[assessment_rows, predictors, drop = FALSE]
  manual <- glmnet::glmnet(
    x = as.matrix(x_analysis),
    y = d$outcome[analysis_rows],
    family = "binomial",
    alpha = 0.5,
    standardize = TRUE
  )
  lambda <- manual$lambda[min(20L, length(manual$lambda))]
  learner <- sb_glmnet_logistic_learner_adapter(
    alpha = 0.5,
    lambda = lambda,
    standardize = TRUE
  )
  model <- learner$fit(x_analysis, d$outcome[analysis_rows])
  observed <- learner$predict(model, x_assessment)
  expected <- unname(drop(stats::predict(
    manual,
    newx = as.matrix(x_assessment),
    s = lambda,
    type = "response"
  )))

  expect_identical(learner$id, "glmnet_logistic")
  expect_equal(observed, expected, tolerance = 1e-12)
  expect_identical(all(observed >= 0 & observed <= 1), TRUE)

  coefficients <- learner$extract_coefficients(model)
  expected_coefficients <- as.matrix(stats::coef(manual, s = lambda))
  expect_identical(
    coefficients$term,
    c("(Intercept)", predictors)
  )
  expect_equal(
    coefficients$coefficient,
    unname(expected_coefficients[, 1L]),
    tolerance = 1e-12
  )
  design <- cbind("(Intercept)" = 1, as.matrix(x_assessment))
  reconstructed <- unname(drop(stats::plogis(
    design %*% coefficients$coefficient[
      match(colnames(design), coefficients$term)
    ]
  )))
  expect_equal(observed, reconstructed, tolerance = 1e-12)
})

test_that("the glmnet logistic adapter locks schema and settings", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 60)
  manual <- glmnet::glmnet(
    x = as.matrix(d[c("gene_a", "gene_b")]),
    y = d$outcome,
    family = "binomial"
  )
  lambda <- manual$lambda[min(10L, length(manual$lambda))]
  learner <- sb_glmnet_logistic_learner_adapter(1, lambda)
  model <- learner$fit(d[c("gene_a", "gene_b")], d$outcome)

  schema_error <- capture_sb_error(
    learner$predict(model, d[c("gene_b", "gene_a")])
  )
  expect_s3_class(schema_error, "stormbreaker_learner_schema_mismatch")

  alpha_error <- capture_sb_error(
    sb_glmnet_logistic_learner_adapter(1.2, lambda)
  )
  expect_s3_class(alpha_error, "stormbreaker_bad_glmnet_alpha")

  lambda_error <- capture_sb_error(
    sb_glmnet_logistic_learner_adapter(1, 0)
  )
  expect_s3_class(lambda_error, "stormbreaker_bad_glmnet_lambda")
})

test_that("the glmnet logistic adapter exposes dependency and warnings", {
  d <- make_sb_binary_fixture(n = 40)
  local_mocked_bindings(sb_glmnet_available = function() FALSE)
  learner <- sb_glmnet_logistic_learner_adapter(1, 0.05)
  missing_error <- capture_sb_error(
    learner$fit(d[c("gene_a", "gene_b")], d$outcome)
  )
  expect_s3_class(missing_error, "stormbreaker_missing_glmnet")

  skip_if_not_installed("glmnet")
  local_mocked_bindings(sb_glmnet_available = function() TRUE)
  original_glmnet <- glmnet::glmnet
  local_mocked_bindings(
    glmnet = function(...) {
      warning("synthetic logistic glmnet warning", call. = FALSE)
      original_glmnet(...)
    },
    .package = "glmnet"
  )
  warning_error <- capture_sb_error(
    learner$fit(d[c("gene_a", "gene_b")], d$outcome)
  )
  expect_s3_class(
    warning_error,
    "stormbreaker_glmnet_logistic_fit_warning"
  )
  expect_match(
    conditionMessage(warning_error),
    "synthetic logistic glmnet warning"
  )
})

test_that("glmnet logistic rejects an ambiguous intercept feature", {
  skip_if_not_installed("glmnet")
  d <- make_sb_binary_fixture(n = 40)
  x <- d[c("gene_a", "gene_b")]
  names(x)[1L] <- "(Intercept)"
  learner <- sb_glmnet_logistic_learner_adapter(1, 0.05)

  error <- capture_sb_error(learner$fit(x, d$outcome))

  expect_s3_class(error, "stormbreaker_reserved_intercept_name")
})

test_that("glmnet logistic tolerates only roundoff at path endpoints", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 60)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  native <- glmnet::glmnet(as.matrix(x), d$outcome, family = "binomial")
  endpoints <- range(native$lambda)
  near <- endpoints * c(1 - 2 * .Machine$double.eps,
    1 + 2 * .Machine$double.eps)
  for (lambda in near) {
    learner <- sb_glmnet_logistic_learner_adapter(1, lambda)
    model <- learner$fit(x, d$outcome)
    expect_identical(model$lambda, native$lambda)
    expect_identical(attr(model, "stormbreaker_lambda"), lambda)
    expected <- as.numeric(predict(native, as.matrix(x), s = lambda,
      type = "response"))
    expect_equal(learner$predict(model, x), expected, tolerance = 0)
    state <- learner$extract_state(model)
    expect_equal(learner$predict_state(state, x)$probability, expected,
      tolerance = 1e-12)
  }
  for (lambda in endpoints * c(1 - 1e-6, 1 + 1e-6)) {
    learner <- sb_glmnet_logistic_learner_adapter(1, lambda)
    error <- capture_sb_error(learner$fit(x, d$outcome))
    expect_s3_class(error, "stormbreaker_glmnet_lambda_outside_path")
  }
  for (scale in c(1e-200, 1, 1e200)) {
    path <- c(0.01, 1) * scale
    expect_true(sb_glmnet_lambda_in_path(path[[1L]] *
      (1 - 2 * .Machine$double.eps), path))
    expect_true(sb_glmnet_lambda_in_path(path[[2L]] *
      (1 + 2 * .Machine$double.eps), path))
    expect_false(sb_glmnet_lambda_in_path(path[[1L]] * (1 - 1e-6), path))
    expect_false(sb_glmnet_lambda_in_path(path[[2L]] * (1 + 1e-6), path))
  }
})

# ---- glmnet-logistic-tuning ----

test_that("binary inner CV selects lambda reproducibly", {
  skip_if_not_installed("glmnet")

  d <- make_sb_high_dim_binary_fixture(n = 60, p = 30)
  predictors <- grep("^gene_", names(d), value = TRUE)
  first <- sb_select_glmnet_logistic_lambda(
    d[predictors], d$outcome,
    inner_folds = 5L, seed = 1701L, alpha = 0.5
  )
  repeated <- sb_select_glmnet_logistic_lambda(
    d[predictors], d$outcome,
    inner_folds = 5L, seed = 1701L, alpha = 0.5
  )

  expect_identical(first$fold_id, repeated$fold_id)
  expect_equal(first$lambda, repeated$lambda, tolerance = 0)
  expect_identical(first$curve, repeated$curve)
  expect_equal(
    first$lambda,
    first$curve$lambda[which.min(first$curve$mean_deviance)],
    tolerance = 0
  )
  expect_identical(first$selection_rule, "lambda.min")
  expect_identical(first$measure, "deviance")
  fold_has_both_classes <- vapply(seq_len(5L), function(fold) {
    setequal(unique(d$outcome[first$fold_id == fold]), c(0L, 1L))
  }, logical(1))
  expect_identical(fold_has_both_classes, rep(TRUE, 5L))
})

test_that("binary inner tuning matches direct cv.glmnet", {
  skip_if_not_installed("glmnet")

  d <- make_sb_high_dim_binary_fixture(n = 60, p = 20)
  predictors <- grep("^gene_", names(d), value = TRUE)
  tuned <- sb_select_glmnet_logistic_lambda(
    d[predictors], d$outcome,
    inner_folds = 5L, seed = 101L, alpha = 0.3,
    selection_rule = "lambda.1se", measure = "auc"
  )
  manual <- glmnet::cv.glmnet(
    x = as.matrix(d[predictors]),
    y = d$outcome,
    family = "binomial",
    alpha = 0.3,
    foldid = tuned$fold_id,
    type.measure = "auc",
    standardize = TRUE
  )

  expect_identical(tuned$selection_rule, "lambda.1se")
  expect_identical(tuned$measure, "auc")
  expect_equal(tuned$lambda, unname(manual$lambda.1se), tolerance = 0)
  expect_gt(tuned$lambda, unname(manual$lambda.min))
  expect_equal(tuned$curve$lambda, manual$lambda, tolerance = 0)
  expect_equal(tuned$curve$mean_measure, manual$cvm, tolerance = 0)
  expect_equal(tuned$curve$standard_error, manual$cvsd, tolerance = 0)
  expect_equal(tuned$curve$nonzero, unname(manual$nzero), tolerance = 0)
})

test_that("outer assessment data cannot alter binary inner tuning", {
  skip_if_not_installed("glmnet")

  d <- make_sb_high_dim_binary_fixture(n = 70, p = 20)
  analysis_rows <- 1:55
  assessment_rows <- 56:70
  predictors <- grep("^gene_", names(d), value = TRUE)
  select_from_analysis <- function(data) {
    sb_select_glmnet_logistic_lambda(
      data[analysis_rows, predictors, drop = FALSE],
      data$outcome[analysis_rows],
      inner_folds = 5L,
      seed = 101L,
      alpha = 1
    )
  }
  original <- select_from_analysis(d)
  changed <- d
  changed$outcome[assessment_rows] <- 1L - changed$outcome[assessment_rows]
  changed[assessment_rows, predictors] <-
    changed[assessment_rows, predictors] * 10000
  perturbed <- select_from_analysis(changed)

  expect_identical(original, perturbed)
})

test_that("binary inner tuning records dependency, folds, and warnings", {
  d <- make_sb_binary_fixture(n = 20)
  local_mocked_bindings(sb_glmnet_available = function() FALSE)
  missing_error <- capture_sb_error(sb_select_glmnet_logistic_lambda(
    d[c("gene_a", "gene_b")], d$outcome,
    inner_folds = 5L, seed = 1L
  ))
  expect_s3_class(missing_error, "stormbreaker_missing_glmnet")

  local_mocked_bindings(sb_glmnet_available = function() TRUE)
  small_outcome <- c(rep(1L, 2), rep(0L, 8))
  fold_error <- capture_sb_error(sb_select_glmnet_logistic_lambda(
    data.frame(a = 1:10, b = 11:20), small_outcome,
    inner_folds = 3L, seed = 1L
  ))
  expect_s3_class(fold_error, "stormbreaker_infeasible_inner_folds")

  skip_if_not_installed("glmnet")
  original_cv_glmnet <- glmnet::cv.glmnet
  local_mocked_bindings(
    cv.glmnet = function(...) {
      warning("synthetic binary tuning warning", call. = FALSE)
      original_cv_glmnet(...)
    },
    .package = "glmnet"
  )
  warning_error <- capture_sb_error(sb_select_glmnet_logistic_lambda(
    d[c("gene_a", "gene_b")], d$outcome,
    inner_folds = 5L, seed = 1L
  ))
  expect_s3_class(
    warning_error,
    "stormbreaker_glmnet_logistic_tuning_warning"
  )
  expect_match(conditionMessage(warning_error), "synthetic binary tuning")
})

test_that("binary inner tuning exposes glmnet fitting errors", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 20)
  x <- data.frame(constant_a = rep(1, 20), constant_b = rep(2, 20))
  error <- capture_sb_error(sb_select_glmnet_logistic_lambda(
    x, d$outcome, inner_folds = 5L, seed = 1L
  ))

  expect_s3_class(error, "stormbreaker_glmnet_logistic_tuning_failure")
  expect_s3_class(error$parent, "simpleError")
})
})

# ---- model-sparse ----
local({
# ---- abess-adapters ----

test_that("ABESS components declare one fixed GIC specification per task", {
  expect_identical(
    sb_abess_learner_ids(),
    c("abess_binary_gic", "abess_cox_gic")
  )
  expect_identical(
    sb_abess_selector_ids(),
    c("abess_binary_gic_nonzero", "abess_cox_gic_nonzero")
  )

  binary <- sb_abess_learner_params("abess_binary_gic")
  survival <- sb_abess_learner_params("abess_cox_gic")
  expect_identical(binary$family, "binomial")
  expect_identical(binary$normalize, 2L)
  expect_identical(survival$family, "cox")
  expect_identical(survival$normalize, 3L)
  expect_identical(binary$tune_path, "sequence")
  expect_identical(binary$tune_type, "gic")
  expect_identical(binary$minimum_support, 0L)
  expect_identical(binary$maximum_support, 10L)
  expect_identical(binary$num_threads, 1L)
  expect_identical(binary$seed, 1L)
  expect_identical(sb_abess_dependency()$package, "abess")
})

test_that("ABESS learners match native predictions and portable coefficients", {
  skip_if_not_installed("abess", minimum_version = "0.4.11")

  binary_data <- make_sb_binary_fixture(n = 84L, seed = 8611L)
  binary_x <- binary_data[c("gene_a", "gene_b", "gene_c")]
  binary_new <- binary_x[1:7, , drop = FALSE]
  binary_adapter <- sb_abess_learner_adapter("abess_binary_gic")

  set.seed(17L)
  caller_rng <- .Random.seed
  binary_model <- binary_adapter$fit(binary_x, binary_data$outcome)
  expect_identical(.Random.seed, caller_rng)
  expect_s3_class(binary_model$engine, "abess")
  expect_null(binary_model$engine$call)

  binary_observed <- binary_adapter$predict(binary_model, binary_new)
  binary_native <- stats::predict(
    binary_model$engine,
    as.matrix(binary_new),
    type = "response"
  )
  binary_coefficients <- binary_adapter$extract_coefficients(binary_model)
  binary_manual <- stats::plogis(
    binary_coefficients$coefficient[[1L]] +
      as.vector(as.matrix(binary_new) %*%
        binary_coefficients$coefficient[-1L])
  )
  binary_state <- binary_adapter$extract_state(binary_model, list())
  binary_replay <- binary_adapter$predict_state(binary_state, binary_new)

  expect_equal(binary_observed, unname(as.numeric(binary_native)), tolerance = 1e-12)
  expect_equal(binary_observed, unname(binary_manual), tolerance = 1e-12)
  expect_equal(binary_replay$probability, binary_observed, tolerance = 1e-12)
  expect_null(binary_state$engine)
  expect_null(binary_state$call)

  survival_data <- make_sb_fixture(n = 90L, seed = 8612L)
  survival_x <- survival_data[c("gene_a", "gene_b", "gene_c")]
  survival_new <- survival_x[1:7, , drop = FALSE]
  survival_adapter <- sb_abess_learner_adapter("abess_cox_gic")
  survival_model <- survival_adapter$fit(
    survival_x,
    survival_data$time,
    survival_data$event
  )
  survival_observed <- survival_adapter$predict(survival_model, survival_new)
  survival_native <- stats::predict(
    survival_model$engine,
    as.matrix(survival_new),
    type = "link"
  )
  horizon <- stats::median(survival_data$time)
  survival_probability <- survival_adapter$predict_survival(
    survival_model,
    survival_x,
    survival_data$time,
    survival_data$event,
    survival_new,
    horizon
  )
  survival_state <- survival_adapter$extract_state(
    survival_model,
    list(
      x = survival_x,
      time = survival_data$time,
      event = survival_data$event,
      prediction_horizon = horizon,
      time_unit = "month"
    )
  )
  survival_replay <- survival_adapter$predict_state(
    survival_state,
    survival_new
  )

  expect_equal(survival_observed, unname(as.numeric(survival_native)), tolerance = 1e-12)
  expect_equal(survival_replay$risk, survival_observed, tolerance = 1e-12)
  expect_equal(
    survival_replay$survival_probability,
    survival_probability$probability,
    tolerance = 1e-12
  )
  expect_identical(survival_state$learner, "abess_cox_gic")
  expect_null(survival_state$engine)
  expect_null(survival_state$call)
})

test_that("ABESS Cox reports its native one-input boundary exactly", {
  skip_if_not_installed("abess", minimum_version = "0.4.11")
  survival_data <- make_sb_fixture(n = 90L, seed = 8613L)
  binary_data <- make_sb_binary_fixture(n = 84L, seed = 8614L)

  expect_error(
    sb_abess_learner_adapter("abess_cox_gic")$fit(
      survival_data["gene_a"],
      survival_data$time,
      survival_data$event
    ),
    class = "stormbreaker_abess_cox_too_few_predictors"
  )
  expect_s3_class(
    sb_abess_learner_adapter("abess_binary_gic")$fit(
      binary_data["gene_a"], binary_data$outcome
    ),
    "stormbreaker_abess_model"
  )
})

test_that("ABESS selectors retain only fold-local nonzero feature state", {
  skip_if_not_installed("abess", minimum_version = "0.4.11")

  binary_data <- make_sb_binary_fixture(n = 90L, seed = 8621L)
  binary_x <- binary_data[c("gene_a", "gene_b", "gene_c")]
  binary_adapter <- sb_abess_selector_adapter(
    "abess_binary_gic_nonzero"
  )
  binary_fit <- binary_adapter$fit(binary_x, binary_data$outcome)

  expect_s3_class(binary_fit, "stormbreaker_selector_fit")
  expect_gte(length(binary_fit$selected_features), 2L)
  expect_identical(
    binary_fit$selected_features,
    binary_fit$scores$feature[binary_fit$scores$selected]
  )
  expect_identical(
    names(binary_adapter$apply(binary_fit, binary_x)),
    binary_fit$selected_features
  )
  expect_null(binary_fit$engine)
  expect_null(binary_fit$call)

  survival_data <- make_sb_fixture(n = 90L, seed = 8622L)
  survival_x <- survival_data[c("gene_a", "gene_b", "gene_c")]
  survival_adapter <- sb_abess_selector_adapter(
    "abess_cox_gic_nonzero"
  )
  survival_fit <- survival_adapter$fit(
    survival_x,
    survival_data$time,
    survival_data$event
  )

  expect_gte(length(survival_fit$selected_features), 2L)
  expect_identical(
    names(survival_adapter$apply(survival_fit, survival_x)),
    survival_fit$selected_features
  )
  expect_null(survival_fit$engine)
  expect_null(survival_fit$call)
})

test_that("ABESS routes expand the registry without redundant self chains", {
  registry <- sb_pipeline_registry()
  restricted_internal_cv <- c(
    setdiff(sb_rule_aft_learner_ids(), "bart_aft_survival"),
    sb_slope_binary_learner_ids(),
    "cancerclass_binary"
  )
  base_binary <- registry$learner[
    registry$task_type == "binary" & registry$selector == "none"
  ]
  base_survival <- registry$learner[
    registry$task_type == "survival" & registry$selector == "none"
  ]
  expect_setequal(
    registry$learner[registry$selector == "abess_binary_gic_nonzero"],
    setdiff(base_binary, c("abess_binary_gic", restricted_internal_cv))
  )
  expect_setequal(
    registry$learner[registry$selector == "abess_cox_gic_nonzero"],
    setdiff(
      base_survival,
      c(
        "abess_cox_gic",
        sb_analysis_tuned_survival_learner_ids(),
        restricted_internal_cv
      )
    )
  )
  expect_identical(
    nrow(registry[
      registry$selector == "abess_binary_gic_nonzero" &
        registry$learner %in% c(
          "abess_binary_gic",
          restricted_internal_cv
        ),
      ,
      drop = FALSE
    ]),
    0L
  )
  expect_identical(
    nrow(registry[
      registry$selector == "abess_cox_gic_nonzero" &
        registry$learner %in% c(
          "abess_cox_gic",
          sb_analysis_tuned_survival_learner_ids(),
          restricted_internal_cv
        ),
      ,
      drop = FALSE
    ]),
    0L
  )

  binary <- sb_pipeline(
    selector = "abess_binary_gic_nonzero",
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 1L
  )
  survival <- sb_pipeline(
    selector = "abess_cox_gic_nonzero",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1L
  )
  expect_null(binary$selector_params)
  expect_null(survival$selector_params)
  expect_match(binary$id, "abess_binary_gic_nonzero", fixed = TRUE)
  expect_match(survival$id, "abess_cox_gic_nonzero", fixed = TRUE)
})

test_that("ABESS learner and selector roles run through final fitting", {
  skip_if_not_installed("abess", minimum_version = "0.4.11")

  binary_data <- make_sb_binary_fixture(n = 84L, seed = 8631L)
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  binary_fit <- sb_fit_final(
    binary_task,
    sb_pipeline(learner = "abess_binary_gic")
  )
  binary_selected_fit <- sb_fit_final(
    binary_task,
    sb_pipeline(
      selector = "abess_binary_gic_nonzero",
      learner = "logistic"
    )
  )

  expect_identical(binary_fit$status, "success")
  expect_identical(binary_selected_fit$status, "success")
  expect_gte(length(binary_selected_fit$required_features), 2L)

  survival_data <- make_sb_fixture(n = 90L, seed = 8632L)
  horizon <- stats::median(survival_data$time)
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = horizon
  )
  survival_fit <- sb_fit_final(
    survival_task,
    sb_pipeline(learner = "abess_cox_gic")
  )
  survival_selected_fit <- sb_fit_final(
    survival_task,
    sb_pipeline(
      selector = "abess_cox_gic_nonzero",
      learner = "coxph"
    )
  )

  expect_identical(survival_fit$status, "success")
  expect_identical(survival_selected_fit$status, "success")
  expect_gte(length(survival_selected_fit$required_features), 2L)
})

test_that("ABESS selection is refit inside every glmnet inner analysis fold", {
  skip_if_not_installed("abess", minimum_version = "0.4.11")
  skip_if_not_installed("glmnet", minimum_version = "4.1.2")

  binary_data <- make_sb_binary_fixture(n = 90L, seed = 8641L)
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  binary_result <- sb_run_cv(
    binary_task,
    sb_pipeline(
      selector = "abess_binary_gic_nonzero",
      learner = "glmnet_logistic",
      inner_folds = 3L,
      inner_seed = 8642L
    ),
    folds = 3L,
    seed = 8643L
  )

  survival_data <- make_sb_fixture(n = 90L, seed = 8644L)
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  survival_result <- sb_run_cv(
    survival_task,
    sb_pipeline(
      selector = "abess_cox_gic_nonzero",
      learner = "glmnet_cox",
      inner_folds = 3L,
      inner_seed = 8645L
    ),
    folds = 3L,
    seed = 8646L
  )

  for (result in list(binary_result, survival_result)) {
    expect_identical(result$status, "success")
    expect_true(result$complete_oof)
    expect_true(all(vapply(result$tuning, function(audit) {
      identical(audit$status, "success") &&
        identical(
          audit$inner_selection_method,
          "glmnet_exclude_function"
        ) &&
        length(audit$inner_selections) == 3L &&
        all(vapply(audit$inner_selections, function(selection) {
          length(selection$selected_features) >= 2L &&
            length(selection$analysis_patient_ids) > 0L
        }, logical(1)))
    }, logical(1))))
    expect_identical(
      result$provenance$abess_version,
      as.character(utils::packageVersion("abess"))
    )
  }
})

# ---- algorithm-ncvreg-path ----

test_that("MCP and SCAD expose their actual paths for both task types", {
  skip_if_not_installed("ncvreg")
  for (kind in c("binary", "survival")) {
    d <- if (kind == "binary") make_sb_binary_fixture(n = 140L, seed = 108L) else
      make_sb_fixture(n = 160L, seed = 8612L)
    task <- if (kind == "binary") sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")) else
      sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
    for (penalty in c("mcp", "scad")) {
      learner <- paste0("ncvreg_", penalty, if (kind == "binary") "_binary_bic" else "_cox_bic")
      fit <- sb_fit_final(task, sb_pipeline(learner = learner), capture_diagnostics = TRUE)
      expect_identical(fit$status, "success")
      native <- fit$diagnostic_artifacts$learner$diagnostic_path
      q <- sb_algorithm_diagnostics(fit, "coefficient_path")
      expect_equal(q$table$coefficient, as.vector(as.matrix(native$beta)))
      expect_equal(unique(q$table$lambda), native$lambda)
      b <- sb_algorithm_diagnostics(fit, "bic_path")
      expect_equal(b$table$bic, as.numeric(native$bic))
      expect_equal(b$selected_lambda, b$table$lambda[which.min(b$table$bic)])
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
        expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(b)), "ggplot_built")
      }
    }
  }
})

# ---- ncvreg-adapters ----

make_ncvreg_binary_fixture <- function(
  n = 140L, p = 6L, seed = 42101L, correlated = FALSE
) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    if (correlated && p > 1L) {
      rho <- 0.8
      sigma <- outer(seq_len(p), seq_len(p), function(i, j) {
        rho^abs(i - j)
      })
      x <- x %*% chol(sigma)
    }
    beta <- numeric(p)
    beta[seq_len(min(p, 4L))] <- c(1.0, -0.8, 0.6, -0.4)[
      seq_len(min(p, 4L))
    ]
    outcome <- stats::rbinom(n, 1L, stats::plogis(drop(x %*% beta)))
    names <- if (p == 1L) {
      "gene one"
    } else {
      c("gene one", "2gene", paste0("gene-", seq_len(p - 2L)))
    }
    data.frame(
      patient_id = sprintf("NB%03d", seq_len(n)),
      outcome = outcome,
      x,
      check.names = FALSE
    ) |>
      stats::setNames(c("patient_id", "outcome", names))
  })
}

make_ncvreg_survival_fixture <- function(
  n = 160L, p = 6L, seed = 42201L, tied = FALSE,
  censor_rate = 0.04, correlated = FALSE
) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    if (correlated && p > 1L) {
      rho <- 0.8
      sigma <- outer(seq_len(p), seq_len(p), function(i, j) {
        rho^abs(i - j)
      })
      x <- x %*% chol(sigma)
    }
    beta <- numeric(p)
    beta[seq_len(min(p, 4L))] <- c(0.9, -0.7, 0.5, -0.3)[
      seq_len(min(p, 4L))
    ]
    event_time <- stats::rexp(n, rate = 0.06 * exp(drop(x %*% beta)))
    censor_time <- stats::rexp(n, rate = censor_rate)
    time <- pmin(event_time, censor_time)
    if (tied) {
      time <- pmax(round(time, 1), 0.1)
    }
    names <- if (p == 1L) {
      "gene one"
    } else {
      c("gene one", "2gene", paste0("gene-", seq_len(p - 2L)))
    }
    data.frame(
      patient_id = sprintf("NS%03d", seq_len(n)),
      time = time,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    ) |>
      stats::setNames(c("patient_id", "time", "event", names))
  })
}

ncvreg_native_fit <- function(x, endpoint, params) {
  matrix_x <- as.matrix(x)
  storage.mode(matrix_x) <- "double"
  lambda_min <- if (nrow(matrix_x) > ncol(matrix_x)) {
    params$lambda_min_ratio_n_gt_p
  } else {
    params$lambda_min_ratio_n_le_p
  }
  dfmax <- ncol(matrix_x) + if (identical(params$family, "binomial")) 1L else 0L
  if (identical(params$family, "binomial")) {
    return(ncvreg::ncvreg(
      X = matrix_x,
      y = as.numeric(endpoint$outcome),
      family = "binomial",
      penalty = params$penalty,
      gamma = params$gamma,
      alpha = params$alpha,
      lambda.min = lambda_min,
      nlambda = params$nlambda,
      eps = params$eps,
      max.iter = params$max_iter,
      convex = params$convex_diagnostics,
      dfmax = dfmax,
      penalty.factor = rep(1, ncol(matrix_x)),
      warn = TRUE,
      returnX = FALSE
    ))
  }
  ncvreg::ncvsurv(
    X = matrix_x,
    y = cbind(endpoint$time, endpoint$event),
    penalty = params$penalty,
    gamma = params$gamma,
    alpha = params$alpha,
    lambda.min = lambda_min,
    nlambda = params$nlambda,
    eps = params$eps,
    max.iter = params$max_iter,
    convex = params$convex_diagnostics,
    dfmax = dfmax,
    penalty.factor = rep(1, ncol(matrix_x)),
    warn = TRUE,
    returnX = FALSE
  )
}

test_that("four MCP and SCAD BIC learners are admitted without grids", {
  expected <- c(
    "ncvreg_mcp_binary_bic", "ncvreg_scad_binary_bic",
    "ncvreg_mcp_cox_bic", "ncvreg_scad_cox_bic"
  )
  expect_identical(sb_ncvreg_learner_ids(), expected)

  for (learner in expected) {
    params <- sb_ncvreg_learner_params(learner)
    expect_identical(params$alpha, 1)
    expect_identical(params$nlambda, 100L)
    expect_identical(params$selection, "analysis_fold_bic")
    expect_identical(params$warning_policy, "fail_fit")
    expect_identical(params$return_x, FALSE)
    if (grepl("_mcp_", learner, fixed = TRUE)) {
      expect_identical(params$penalty, "MCP")
      expect_identical(params$gamma, 3)
    } else {
      expect_identical(params$penalty, "SCAD")
      expect_identical(params$gamma, 3.7)
    }
  }

  binary_selectors <- expected_fixed_learner_selectors("binary")
  survival_selectors <- expected_fixed_learner_selectors("survival")
  rows <- sb_pipeline_registry()
  rows <- rows[rows$learner %in% expected, , drop = FALSE]
  expect_identical(
    nrow(rows),
    2L * length(binary_selectors) + 2L * length(survival_selectors)
  )
  expect_true(all(rows$tuning_semantics == "analysis_fold_bic"))
  expect_true(all(grepl("ncvreg", rows$optional_dependency, fixed = TRUE)))
  expect_setequal(
    unique(rows$selector[rows$task_type == "binary"]),
    binary_selectors
  )
  expect_setequal(
    unique(rows$selector[rows$task_type == "survival"]),
    survival_selectors
  )
})

test_that("binary MCP and SCAD reproduce the native BIC-selected model", {
  skip_if_not_installed("ncvreg", minimum_version = "3.16.0")

  data <- make_ncvreg_binary_fixture(correlated = TRUE)
  features <- setdiff(names(data), c("patient_id", "outcome"))
  x <- data[features]
  assessment <- x[c(2L, 31L, 88L), , drop = FALSE]

  for (learner in sb_ncvreg_binary_learner_ids()) {
    adapter <- sb_ncvreg_learner_adapter(learner)
    set.seed(42102L)
    caller_rng <- .Random.seed
    model <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, caller_rng)

    native <- ncvreg_native_fit(
      x,
      list(outcome = data$outcome),
      adapter$params
    )
    bic <- stats::BIC(native)
    selected <- which.min(bic)
    native_coefficient <- drop(stats::coef(
      native,
      lambda = native$lambda[[selected]]
    ))
    native_probability <- drop(stats::predict(
      native,
      X = as.matrix(assessment),
      type = "response",
      lambda = native$lambda[[selected]]
    ))

    expect_equal(model$intercept, native_coefficient[[1L]], tolerance = 0)
    expect_equal(
      unname(model$coefficient),
      unname(native_coefficient[-1L]),
      tolerance = 0
    )
    expect_equal(
      adapter$predict(model, assessment),
      unname(native_probability),
      tolerance = 1e-14
    )
    expect_identical(model$effective$selected_index, as.integer(selected))
    expect_equal(
      model$effective$selected_lambda,
      native$lambda[[selected]],
      tolerance = 0
    )
    expect_false(any(c(
      "engine", "x", "X", "y", "outcome", "time", "event", "fail"
    ) %in% names(model)))

    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, assessment)
    expect_equal(
      replay$probability,
      unname(native_probability),
      tolerance = 1e-14
    )
    expect_false(any(c(
      "engine", "x", "X", "y", "outcome", "time", "event", "fail"
    ) %in% names(state)))
  }
})

test_that("survival MCP and SCAD reproduce native risk and portable survival", {
  skip_if_not_installed("ncvreg", minimum_version = "3.16.0")

  data <- make_ncvreg_survival_fixture(tied = TRUE, correlated = TRUE)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[c(3L, 41L, 109L), , drop = FALSE]
  horizons <- as.numeric(stats::quantile(data$time, c(0.35, 0.65)))

  for (learner in sb_ncvreg_survival_learner_ids()) {
    adapter <- sb_ncvreg_learner_adapter(learner)
    set.seed(42202L)
    caller_rng <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, caller_rng)

    native <- ncvreg_native_fit(
      x,
      list(time = data$time, event = data$event),
      adapter$params
    )
    selected <- which.min(stats::BIC(native))
    native_risk <- drop(stats::predict(
      native,
      X = as.matrix(assessment),
      type = "link",
      lambda = native$lambda[[selected]]
    ))
    risk <- adapter$predict(model, assessment)
    early <- adapter$predict_survival(
      model, x, data$time, data$event, assessment, horizons[[1L]]
    )
    late <- adapter$predict_survival(
      model, x, data$time, data$event, assessment, horizons[[2L]]
    )

    expect_equal(risk, unname(native_risk), tolerance = 1e-14)
    expect_true(all(early$probability >= late$probability))
    expect_true(all(early$probability >= 0 & early$probability <= 1))
    expect_identical(early$ties, "breslow")
    expect_identical(
      early$source,
      paste0(learner, "_breslow_native_link")
    )

    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = horizons[[1L]],
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, assessment)
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      early$probability,
      tolerance = 0
    )
    expect_false(any(c(
      "engine", "x", "X", "y", "outcome", "time", "event", "fail"
    ) %in% names(state)))
  }
})

test_that("one-feature and heavy-censoring fits stay finite", {
  skip_if_not_installed("ncvreg", minimum_version = "3.16.0")

  binary <- make_ncvreg_binary_fixture(p = 1L, seed = 42301L)
  survival <- make_ncvreg_survival_fixture(
    p = 1L,
    seed = 42302L,
    censor_rate = 0.22
  )
  binary_x <- binary[setdiff(names(binary), c("patient_id", "outcome"))]
  survival_x <- survival[setdiff(
    names(survival),
    c("patient_id", "time", "event")
  )]

  for (learner in sb_ncvreg_binary_learner_ids()) {
    adapter <- sb_ncvreg_learner_adapter(learner)
    model <- adapter$fit(binary_x, binary$outcome)
    probability <- adapter$predict(model, binary_x[1:5, , drop = FALSE])
    expect_true(all(is.finite(probability)))
    expect_true(all(probability >= 0 & probability <= 1))
  }
  for (learner in sb_ncvreg_survival_learner_ids()) {
    adapter <- sb_ncvreg_learner_adapter(learner)
    model <- adapter$fit(survival_x, survival$time, survival$event)
    risk <- adapter$predict(model, survival_x[1:5, , drop = FALSE])
    expect_true(all(is.finite(risk)))
  }
})

test_that("native path warnings are preserved as typed fit failures", {
  skip_if_not_installed("ncvreg", minimum_version = "3.16.0")

  cases <- list(
    binary = make_ncvreg_binary_fixture(n = 80L, p = 120L, seed = 41001L),
    survival = make_ncvreg_survival_fixture(
      n = 80L,
      p = 120L,
      seed = 42001L
    )
  )
  for (task in names(cases)) {
    data <- cases[[task]]
    endpoints <- if (identical(task, "binary")) {
      c("patient_id", "outcome")
    } else {
      c("patient_id", "time", "event")
    }
    x <- data[setdiff(names(data), endpoints)]
    learners <- if (identical(task, "binary")) {
      sb_ncvreg_binary_learner_ids()
    } else {
      sb_ncvreg_survival_learner_ids()
    }
    endpoint <- if (identical(task, "binary")) {
      list(outcome = data$outcome)
    } else {
      list(time = data$time, event = data$event)
    }
    for (learner in learners) {
      params <- sb_ncvreg_learner_params(learner)
      native_warnings <- character()
      native <- withCallingHandlers(
        ncvreg_native_fit(x, endpoint, params),
        warning = function(warning) {
          native_warnings <<- c(native_warnings, conditionMessage(warning))
          invokeRestart("muffleWarning")
        }
      )
      expect_s3_class(native, "ncvreg")

      adapter <- sb_ncvreg_learner_adapter(learner)
      fit_call <- if (identical(task, "binary")) {
        quote(adapter$fit(x, data$outcome))
      } else {
        quote(adapter$fit(x, data$time, data$event))
      }
      if (length(native_warnings)) {
        expect_error(
          eval(fit_call),
          class = "stormbreaker_ncvreg_fit_warning"
        )
      } else {
        expect_s3_class(eval(fit_call), "stormbreaker_ncvreg_model")
      }
    }
  }
})

test_that("MCP and SCAD execute through CV and final prediction", {
  skip_if_not_installed("ncvreg", minimum_version = "3.16.0")

  binary <- make_ncvreg_binary_fixture(seed = 42401L)
  binary_features <- setdiff(names(binary), c("patient_id", "outcome"))
  binary_task <- sb_binary_task(
    binary,
    "patient_id",
    "outcome",
    binary_features
  )
  survival <- make_ncvreg_survival_fixture(seed = 42402L)
  survival_features <- setdiff(
    names(survival),
    c("patient_id", "time", "event")
  )
  survival_task <- sb_surv_task(
    survival,
    "patient_id",
    "time",
    "event",
    survival_features,
    "month",
    prediction_horizon = as.numeric(stats::quantile(survival$time, 0.5))
  )

  for (learner in sb_ncvreg_binary_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 4L
    )
    comparison <- sb_run_cv(binary_task, pipeline, folds = 3L, seed = 42403L)
    fit <- sb_fit_final(binary_task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      binary[1:6, c("patient_id", rev(binary_features)), drop = FALSE],
      "patient_id"
    )
    expect_identical(comparison$status, "success")
    expect_true(comparison$complete_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(prediction$predictions$probability >= 0 &
      prediction$predictions$probability <= 1))
  }

  for (learner in sb_ncvreg_survival_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 4L
    )
    comparison <- sb_run_cv(
      survival_task,
      pipeline,
      folds = 3L,
      seed = 42404L
    )
    fit <- sb_fit_final(survival_task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      survival[1:6, c("patient_id", rev(survival_features)), drop = FALSE],
      "patient_id"
    )
    expect_identical(comparison$status, "success")
    expect_true(comparison$complete_oof)
    expect_true(comparison$complete_survival_probability_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})
})

# ---- model-slope-randomglm ----
local({
# ---- algorithm-randomglm ----

test_that("randomGLM diagnostics preserve original bags and probability averaging", {
  skip_if_not_installed("randomGLM")
  d <- sb_example_binary(n = 160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features),
    sb_pipeline(learner = "random_glm_binary"))
  expect_identical(fit$status, "success")
  expect_null(fit$diagnostic_artifacts)
  before <- serialize(fit, NULL, version = 2L)
  state <- fit$learner_state
  counts <- sb_randomglm_diagnostics(fit)
  expected_counts <- vapply(seq_along(features), function(j)
    sum(vapply(state$bag_features, function(index) j %in% index, logical(1))), numeric(1))
  expect_equal(counts$table$bags, expected_counts)
  expect_equal(counts$table$fraction, expected_counts / length(state$intercepts))
  expect_identical(counts$table$feature, features)
  expect_identical(counts$intercepts, state$intercepts)
  coefficients <- sb_randomglm_diagnostics(fit, "bag_coefficients")
  expect_equal(nrow(coefficients$table), length(state$intercepts) * length(features))
  beta <- matrix(coefficients$table$coefficient, nrow = length(state$intercepts))
  for (i in seq_along(state$intercepts)) {
    expect_equal(unname(beta[i, state$bag_features[[i]]]), state$bag_coefficients[[i]])
  }
  score <- sweep(as.matrix(d[features]) %*% t(beta), 2L, state$intercepts, "+")
  expect_equal(rowMeans(plogis(score)), sb_predict_final(fit, d, "patient_id")$predictions$probability,
    tolerance = 1e-12)
  sizes <- sb_randomglm_diagnostics(fit, "bag_size")
  expect_equal(sizes$table$selected_features, unname(lengths(state$bag_features)))
  expect_identical(sizes$table$bag, seq_along(state$bag_features))
  expect_error(sb_nomogram(fit, d), class = "stormbreaker_unsupported_coefficients")
  expect_error(sb_randomglm_diagnostics(fit, "cv_curve"), class = "stormbreaker_unsupported_diagnostic")
  bundle <- sb_native_diagnostics(fit, complete = TRUE)
  expect_identical(bundle$status$status, rep("available", 3L))
  if (requireNamespace("ggplot2", quietly = TRUE)) for (view in bundle$diagnostics) {
    plot <- sb_plot_algorithm_diagnostic(view)
    expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
    expect_null(plot$labels$subtitle)
    expect_null(plot$labels$caption)
  }
  expect_identical(serialize(fit, NULL, version = 2L), before)

  synthetic <- fit
  synthetic$learner_state$intercepts <- c(0, 0, 0)
  synthetic$learner_state$bag_features <- list(integer(), 2L, 1L)
  synthetic$learner_state$bag_coefficients <- list(numeric(), 0, 1)
  synthetic$learner_state$params$n_bags <- 3L
  expect_equal(sb_randomglm_diagnostics(synthetic)$table$bags, c(1, 1, 0))
  expect_equal(sb_randomglm_diagnostics(synthetic, "bag_size")$table$selected_features, c(0, 1, 1))
  tab <- sb_randomglm_diagnostics(synthetic, "bag_coefficients")$table
  expect_true(tab$included[tab$bag == 2L & tab$feature == features[2L]])
  expect_equal(tab$coefficient[tab$bag == 2L & tab$feature == features[2L]], 0)
  expect_false(any(tab$included[tab$bag == 1L]))
})

# ---- slope-randomglm-binary-adapters ----

make_slope_randomglm_fixture <- function(n = 90L, p = 100L, seed = 2901L) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("gene_", seq_len(p))
  linear_predictor <- 1.3 * x[, 1L] - x[, 2L] + 0.6 * x[, 3L]
  outcome <- stats::rbinom(n, 1L, stats::plogis(linear_predictor))
  data <- data.frame(
    patient_id = paste0("patient_", seq_len(n)),
    outcome = outcome,
    x,
    check.names = FALSE
  )
  list(data = data, x = data[colnames(x)], outcome = outcome)
}


test_that("missing SLOPE and randomGLM engines fail explicitly", {
  fixture <- make_slope_randomglm_fixture(n = 30L, p = 5L)

  for (learner in sb_slope_randomglm_binary_learner_ids()) {
    package <- if (identical(learner, "slope_binary")) "SLOPE" else "randomGLM"
    expected_class <- if (identical(learner, "slope_binary")) {
      "stormbreaker_missing_slope"
    } else {
      "stormbreaker_missing_randomglm"
    }
    local_mocked_bindings(
      sb_slope_randomglm_binary_dependency = function(learner) {
        list(package = package, available = FALSE)
      }
    )
    expect_error(
      sb_slope_randomglm_binary_learner_adapter(learner)$fit(
        fixture$x,
        fixture$outcome
      ),
      class = expected_class
    )
  }
})


test_that("compact randomGLM probabilities equal its native ensemble", {
  skip_if_not_installed("randomGLM")
  fixture <- make_slope_randomglm_fixture(n = 70L, p = 100L, seed = 2903L)
  params <- sb_slope_randomglm_binary_params("random_glm_binary")
  params$n_bags <- 10L
  native <- randomGLM::randomGLM(
    x = as.matrix(fixture$x),
    y = factor(fixture$outcome, levels = c(0L, 1L)),
    type = "binary",
    maxInteractionOrder = params$max_interaction_order,
    nBags = params$n_bags,
    replace = params$replace,
    nCandidateCovariates = params$n_candidate_covariates,
    keepModels = TRUE,
    randomSeed = params$seed,
    nThreads = params$threads,
    verbose = 0
  )
  compact <- sb_compact_random_glm_binary(native, names(fixture$x), params)

  expected <- unname(stats::predict(
    native,
    newdata = as.matrix(fixture$x),
    type = "response"
  )[, "1"])
  observed <- sb_predict_slope_randomglm_binary(compact, fixture$x)
  state <- sb_extract_slope_randomglm_binary_state(compact)

  expect_equal(observed, expected, tolerance = 1e-14)
  expect_equal(
    sb_predict_slope_randomglm_binary_state(state, fixture$x)$probability,
    observed,
    tolerance = 0
  )
  expect_false(any(c("models", "training_data", "outcome") %in% names(state)))
})

test_that("randomGLM fitting preserves an absent caller RNG state", {
  skip_if_not_installed("randomGLM")
  fixture <- make_slope_randomglm_fixture(n = 50L, p = 30L, seed = 2906L)
  params <- sb_slope_randomglm_binary_params("random_glm_binary")
  params$n_bags <- 3L
  previous_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit(
    {
      if (is.null(previous_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  sb_fit_random_glm_binary(
    list(x = fixture$x, outcome = fixture$outcome),
    params
  )

  expect_false(exists(".Random.seed", envir = .GlobalEnv))
})

test_that("both learners run through binary CV and final prediction", {
  skip_if_not_installed("SLOPE")
  skip_if_not_installed("randomGLM")
  production_params <- sb_slope_randomglm_binary_params
  local_mocked_bindings(
    sb_slope_randomglm_binary_params = function(learner) {
      params <- production_params(learner)
      if (identical(learner, "slope_binary")) {
        params$q <- c(0.1, 0.2)
        params$inner_folds <- 3L
        params$path_length <- 10L
      } else {
        params$n_bags <- 10L
      }
      params
    }
  )
  fixture <- make_slope_randomglm_fixture(n = 90L, p = 30L, seed = 2904L)
  task <- sb_binary_task(
    fixture$data,
    "patient_id",
    "outcome",
    names(fixture$x)
  )

  for (learner in sb_slope_randomglm_binary_learner_ids()) {
    pipeline <- sb_pipeline(learner = learner)
    cv <- sb_run_cv(task, pipeline = pipeline, folds = 3L, seed = 2905L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      fixture$data[1:6, c("patient_id", names(fixture$x))],
      "patient_id"
    )
    version_field <- if (identical(learner, "slope_binary")) {
      "slope_version"
    } else {
      "randomglm_version"
    }

    expect_identical(cv$status, "success", info = learner)
    expect_identical(cv$complete_oof, TRUE, info = learner)
    expect_true(all(is.finite(cv$predictions$probability)), info = learner)
    expect_true(version_field %in% names(cv$provenance), info = learner)
    expect_identical(fit$status, "success", info = learner)
    expect_identical(prediction$status, "success", info = learner)
    expect_true(all(is.finite(prediction$predictions$probability)), info = learner)
    expect_true(version_field %in% names(fit$provenance), info = learner)
  }
})
})
