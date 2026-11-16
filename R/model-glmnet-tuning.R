# ---- model-glmnet-tuning ----

# ---- algorithm-glmnet-tuning ----

sb_glmnet_tuning_diagnostic <- function(fit, type) {
  if (!fit$pipeline$learner %in% c("glmnet_logistic", "glmnet_cox")) sb_abort(
    "This diagnostic requires a saved glmnet inner-validation curve.",
    "stormbreaker_unsupported_diagnostic")
  tuning <- fit$tuning
  tab <- tuning$curve
  value_col <- if (identical(tuning$measure, "deviance")) "mean_deviance" else "mean_measure"
  required <- c("lambda", value_col, "standard_error", "nonzero")
  valid <- identical(tuning$status, "success") && is.data.frame(tab) && nrow(tab) > 0 &&
    all(required %in% names(tab)) && all(vapply(tab[required], function(z) {
      is.numeric(z) && all(is.finite(z))
    }, logical(1))) && all(tab$lambda > 0) && !anyDuplicated(tab$lambda) &&
    all(tab$standard_error >= 0) && all(tab$nonzero >= 0) &&
    all(tab$nonzero == floor(tab$nonzero)) &&
    is.numeric(tuning$lambda) && length(tuning$lambda) == 1L &&
    is.finite(tuning$lambda) && tuning$lambda %in% tab$lambda
  if (!valid) sb_abort("Saved glmnet tuning curve is missing or invalid.",
    "stormbreaker_invalid_diagnostic")
  tab$mean_value <- tab[[value_col]]
  tab$selected <- tab$lambda == tuning$lambda
  labels <- c(deviance = "Deviance", auc = "ROC AUC", C = "Concordance",
    class = "Classification error", mse = "Mean squared error", mae = "Mean absolute error")
  measure <- unname(labels[tuning$measure])
  if (length(measure) != 1L || is.na(measure)) measure <- tuning$measure
  structure(list(type = type, table = tab, learner = fit$pipeline$learner,
    pipeline_id = fit$pipeline$id, measure = measure,
    selected_lambda = tuning$lambda, selection_rule = tuning$selection_rule,
    alpha = tuning$alpha, source = "saved final-development inner cross-validation",
    uncertainty = "one saved standard error; not a confidence interval"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_glmnet_tuning_diagnostic <- function(x) {
  lambda <- mean_value <- standard_error <- nonzero <- NULL
  p <- ggplot2::ggplot(x$table, ggplot2::aes(lambda))
  if (x$type == "cv_curve") {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_value - standard_error,
      ymax = mean_value + standard_error), width = 0, colour = "#A6A6A6") +
      ggplot2::geom_line(ggplot2::aes(y = mean_value), colour = "#2878A0") +
      ggplot2::geom_point(ggplot2::aes(y = mean_value), size = 1, colour = "#2878A0") +
      ggplot2::labs(title = "CV error", y = x$measure,
        caption = "Bars: saved mean +/- one standard error, not a confidence interval.\nDevelopment inner validation, not untouched-cohort performance; no refitting.")
  } else {
    p <- p + ggplot2::geom_step(ggplot2::aes(y = nonzero), colour = "#2878A0") +
      ggplot2::geom_point(ggplot2::aes(y = nonzero), size = 1, colour = "#2878A0") +
      ggplot2::scale_y_continuous(limits = c(0, max(1, x$table$nonzero)),
        breaks = function(limits) unique(round(pretty(limits, n = 5)))) +
      ggplot2::labs(title = "Model size",
        y = "Nonzero predictors",
        caption = "Counts retained from the original glmnet tuning object.\nThis is model complexity, not validation performance; no refitting.")
  }
  p + ggplot2::geom_vline(xintercept = x$selected_lambda, linetype = 2, colour = "#D77935") +
    ggplot2::scale_x_log10() + ggplot2::labs(x = "Lambda (log scale)",
      subtitle = paste0("Original rule: ", x$selection_rule, "; alpha = ", x$alpha,
        "; selected lambda = ", format(signif(x$selected_lambda, 4)))) + sb_viz_theme()
}

# ---- glmnet-cox-tuning ----

sb_validate_glmnet_lambda_rule <- function(value,
                                           argument = "selection_rule") {
  valid <- is.character(value) && length(value) == 1L && !is.na(value) &&
    value %in% c("lambda.min", "lambda.1se")
  if (!valid) {
    sb_abort(
      paste0(
        "`", argument,
        "` must be either `\"lambda.min\"` or `\"lambda.1se\"`."
      ),
      "stormbreaker_bad_glmnet_lambda_rule"
    )
  }
  unname(as.character(value))
}

sb_validate_glmnet_tuning_measure <- function(value, task_type,
                                              argument = "measure") {
  allowed <- if (identical(task_type, "binary")) {
    c("deviance", "auc")
  } else {
    c("deviance", "C")
  }
  valid <- is.character(value) && length(value) == 1L && !is.na(value) &&
    value %in% allowed
  if (!valid) {
    sb_abort(
      paste0(
        "`", argument, "` must be one of ",
        paste0("`\"", allowed, "\"`", collapse = ", "),
        " for a ", task_type, " glmnet learner."
      ),
      "stormbreaker_bad_glmnet_tuning_measure"
    )
  }
  unname(as.character(value))
}

sb_make_stratified_fold_id <- function(event, folds, seed,
                                       stratum_label = "event",
                                       group_ids = NULL) {
  if (!is.null(group_ids)) {
    group_ids <- sb_validate_group_ids(group_ids, length(event))
    if (anyDuplicated(group_ids)) {
      return(sb_make_group_stratified_fold_id(
        event, group_ids, folds, seed, "stormbreaker_infeasible_inner_folds"
      ))
    }
  }
  strata_counts <- table(factor(event, levels = c(0, 1)))
  if (length(event) < folds || any(strata_counts < folds)) {
    sb_abort(
      paste0(
        "Cannot create ", folds,
        " ", stratum_label,
        "-stratified inner folds: each outcome group must have ",
        "at least ", folds, " patients."
      ),
      "stormbreaker_infeasible_inner_folds"
    )
  }

  fold_id <- integer(length(event))
  withr::with_seed(seed, {
    for (outcome in c(0, 1)) {
      rows <- which(event == outcome)
      rows <- sample(rows, length(rows), replace = FALSE)
      fold_id[rows] <- rep(seq_len(folds), length.out = length(rows))
    }
  })
  fold_id
}

sb_is_glmnet_tail_convergence_warning <- function(warning) {
  message <- if (inherits(warning, "condition")) {
    conditionMessage(warning)
  } else {
    as.character(warning)
  }
  grepl(
    paste0(
      "^from glmnet C\\+\\+ code \\(error code -[0-9]+\\); ",
      "Convergence for [0-9]+(st|nd|rd|th) lambda value not reached ",
      "after maxit=.* iterations; solutions for larger lambdas returned$"
    ),
    message
  )
}

sb_validate_glmnet_selected_lambda <- function(x, time, event, fold_id,
                                               lambda, alpha,
                                               standardize,
                                               inner_selections = list()) {
  x_matrix <- as.matrix(x)
  validation_rows <- lapply(sort(unique(fold_id)), function(fold) {
    analysis_rows <- fold_id != fold
    assessment_rows <- fold_id == fold
    selected_features <- colnames(x_matrix)
    if (length(inner_selections) > 0L) {
      selected_features <- inner_selections[[fold]][["selected_features"]]
    }
    path_warnings <- character()
    model <- tryCatch(
      withCallingHandlers(
        glmnet::glmnet(
          x = x_matrix[
            analysis_rows,
            selected_features,
            drop = FALSE
          ],
          y = survival::Surv(
            time[analysis_rows],
            event[analysis_rows]
          ),
          family = "cox",
          cox.ties = "breslow",
          alpha = alpha,
          standardize = standardize
        ),
        warning = function(warning) {
          if (sb_is_glmnet_tail_convergence_warning(warning)) {
            path_warnings <<- c(
              path_warnings,
              conditionMessage(warning)
            )
            invokeRestart("muffleWarning")
          }
          sb_abort(
            paste0(
              "Selected-lambda validation raised a warning in inner ",
              "fold ", fold, ": ", conditionMessage(warning)
            ),
            "stormbreaker_glmnet_tuning_warning",
            parent = warning
          )
        }
      ),
      error = function(error) {
        if (inherits(error, "stormbreaker_error")) {
          stop(error)
        }
        sb_abort(
          paste0(
            "Selected-lambda validation failed in inner fold ", fold,
            ": ", conditionMessage(error)
          ),
          "stormbreaker_glmnet_tuning_failure",
          parent = error
        )
      }
    )
    min_lambda <- min(model$lambda)
    max_lambda <- max(model$lambda)
    if (lambda < min_lambda &&
        !sb_glmnet_lambda_in_path(lambda, model$lambda)) {
      sb_abort(
        paste0(
          "The selected lambda ", format(lambda),
          " is below the smallest converged lambda ",
          format(min_lambda), " in inner fold ", fold, "."
        ),
        "stormbreaker_glmnet_selected_lambda_not_converged"
      )
    }
    risk <- tryCatch(
      drop(stats::predict(
        model,
        newx = x_matrix[
          assessment_rows,
          selected_features,
          drop = FALSE
        ],
        s = lambda,
        type = "link"
      )),
      error = function(error) {
        sb_abort(
          paste0(
            "Selected-lambda prediction failed in inner fold ", fold,
            ": ", conditionMessage(error)
          ),
          "stormbreaker_invalid_glmnet_tuning",
          parent = error
        )
      }
    )
    risk_finite <- length(risk) == sum(assessment_rows) &&
      !anyNA(risk) && all(is.finite(risk))
    if (!risk_finite) {
      sb_abort(
        paste0(
          "Selected-lambda prediction was incomplete or non-finite in ",
          "inner fold ", fold, "."
        ),
        "stormbreaker_invalid_glmnet_tuning"
      )
    }
    data.frame(
      fold = as.integer(fold),
      n_lambda = as.integer(length(model$lambda)),
      min_lambda = min_lambda,
      max_lambda = max_lambda,
      path_warning_count = as.integer(length(path_warnings)),
      path_warnings = paste(path_warnings, collapse = " | "),
      assessment_risk_finite = risk_finite,
      stringsAsFactors = FALSE
    )
  })
  validation <- do.call(rbind, validation_rows)
  rownames(validation) <- NULL
  validation
}

sb_select_glmnet_cox_lambda <- function(x, time, event,
                                        inner_folds = 5L, seed = 1L,
                                        alpha = 1,
                                        standardize = TRUE,
                                        selector = NULL,
                                        fitted_selector = NULL,
                                        patient_ids = NULL,
                                        selection_rule = "lambda.min",
                                        measure = "deviance",
                                        group_ids = NULL,
                                        retain_analysis_fit = FALSE) {
  if (!sb_glmnet_available()) {
    sb_abort(
      "Package `glmnet` is required for inner Cox tuning.",
      "stormbreaker_missing_glmnet"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  inner_folds <- sb_check_positive_integer(
    inner_folds,
    "inner_folds",
    minimum = 3L
  )
  seed <- sb_check_positive_integer(seed, "seed", minimum = 1L)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha < 0 || alpha > 1) {
    sb_abort(
      "`alpha` must be one finite number between 0 and 1.",
      "stormbreaker_bad_glmnet_alpha"
    )
  }
  if (!is.logical(standardize) || length(standardize) != 1L ||
      is.na(standardize)) {
    sb_abort(
      "`standardize` must be `TRUE` or `FALSE`.",
      "stormbreaker_bad_glmnet_standardize"
    )
  }
  selection_rule <- sb_validate_glmnet_lambda_rule(selection_rule)
  measure <- sb_validate_glmnet_tuning_measure(measure, "survival")
  if (any(time <= 0)) {
    sb_abort(
      paste0(
        "Inner glmnet Cox tuning requires strictly positive ",
        "survival times."
      ),
      "stormbreaker_glmnet_nonpositive_time"
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }

  event <- as.integer(event)
  fold_id <- sb_make_stratified_fold_id(
    event, inner_folds, seed, group_ids = group_ids
  )
  uses_fold_local_selection <- !is.null(selector) ||
    !is.null(fitted_selector) || !is.null(patient_ids)
  if (uses_fold_local_selection &&
      (is.null(selector) || is.null(fitted_selector) ||
       is.null(patient_ids))) {
    sb_abort(
      paste0(
        "Fold-local glmnet selection requires `selector`, ",
        "`fitted_selector`, and `patient_ids` together."
      ),
      "stormbreaker_bad_glmnet_selection_context"
    )
  }
  x_matrix <- as.matrix(x)
  fold_filter <- NULL
  if (uses_fold_local_selection) {
    rownames(x_matrix) <- as.character(patient_ids)
    fold_filter <- sb_make_glmnet_fold_selector_filter(
      selector,
      fitted_selector,
      patient_ids,
      fold_id,
      task_type = "survival",
      group_ids = group_ids
    )
  }
  tail_warnings <- character()
  cv_fit <- tryCatch(
    withCallingHandlers(
      if (uses_fold_local_selection) {
        glmnet::cv.glmnet(
          x = x_matrix,
          y = survival::Surv(time, event),
          family = "cox",
          cox.ties = "breslow",
          alpha = as.numeric(alpha),
          foldid = fold_id,
          type.measure = measure,
          standardize = standardize,
          exclude = fold_filter$exclude,
          alignment = "lambda",
          parallel = FALSE
        )
      } else {
        glmnet::cv.glmnet(
          x = x_matrix,
          y = survival::Surv(time, event),
          family = "cox",
          cox.ties = "breslow",
          alpha = as.numeric(alpha),
          foldid = fold_id,
          type.measure = measure,
          standardize = standardize
        )
      },
      warning = function(warning) {
        if (sb_is_glmnet_tail_convergence_warning(warning)) {
          tail_warnings <<- c(
            tail_warnings,
            conditionMessage(warning)
          )
          invokeRestart("muffleWarning")
        }
        sb_abort(
          paste0(
            "Inner glmnet Cox tuning raised a warning: ",
            conditionMessage(warning)
          ),
          "stormbreaker_glmnet_tuning_warning",
          parent = warning
        )
      }
    ),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "Inner glmnet Cox tuning failed: ",
          conditionMessage(error)
        ),
        "stormbreaker_glmnet_tuning_failure",
        parent = error
      )
    }
  )
  inner_selections <- if (uses_fold_local_selection) {
    fold_filter$inner_selections()
  } else {
    list()
  }

  curve <- data.frame(
    lambda = cv_fit$lambda,
    mean_value = cv_fit$cvm,
    standard_error = cv_fit$cvsd,
    nonzero = cv_fit$nzero,
    stringsAsFactors = FALSE
  )
  names(curve)[[2L]] <- if (identical(measure, "deviance")) {
    "mean_deviance"
  } else {
    "mean_measure"
  }
  if (anyNA(cv_fit$cvm) || any(!is.finite(cv_fit$cvm))) {
    sb_abort(
      "Inner glmnet Cox tuning produced an incomplete error curve.",
      "stormbreaker_invalid_glmnet_tuning"
    )
  }

  selected_lambda <- unname(cv_fit[[selection_rule]])
  if (length(selected_lambda) != 1L || !is.finite(selected_lambda) ||
      selected_lambda <= 0 || !selected_lambda %in% curve$lambda) {
    sb_abort(
      "Inner glmnet Cox tuning produced an invalid selected lambda.",
      "stormbreaker_invalid_glmnet_tuning"
    )
  }
  warning_validation <- data.frame(
    fold = integer(),
    n_lambda = integer(),
    min_lambda = numeric(),
    max_lambda = numeric(),
    path_warning_count = integer(),
    path_warnings = character(),
    assessment_risk_finite = logical(),
    stringsAsFactors = FALSE
  )
  if (length(tail_warnings) > 0L) {
    warning_validation <- tryCatch(
      sb_validate_glmnet_selected_lambda(
        x,
        time,
        event,
        fold_id,
        lambda = selected_lambda,
        alpha = as.numeric(alpha),
        standardize = standardize,
        inner_selections = inner_selections
      ),
      error = function(error) {
        error$glmnet_warnings <- tail_warnings
        stop(error)
      }
    )
  }

  result <- list(
    lambda = selected_lambda,
    selection_rule = selection_rule,
    measure = measure,
    alpha = as.numeric(alpha),
    standardize = standardize,
    inner_folds = inner_folds,
    seed = seed,
    fold_id = fold_id,
    alignment = "lambda",
    inner_selection_method = if (uses_fold_local_selection) {
      "glmnet_exclude_function"
    } else {
      "none"
    },
    inner_selections = inner_selections,
    curve = curve,
    warnings = tail_warnings,
    warning_validation = warning_validation
  )
  if (retain_analysis_fit) result$analysis_fit <- cv_fit$glmnet.fit
  result
}

# ---- glmnet-logistic-tuning ----

sb_select_glmnet_logistic_lambda <- function(x, outcome,
                                             inner_folds = 5L,
                                             seed = 1L,
                                             alpha = 1,
                                             standardize = TRUE,
                                             selector = NULL,
                                             fitted_selector = NULL,
                                             patient_ids = NULL,
                                             selection_rule = "lambda.min",
                                             measure = "deviance",
                                             group_ids = NULL,
                                             capture_path = FALSE) {
  if (!sb_glmnet_available()) {
    sb_abort(
      "Package `glmnet` is required for inner logistic tuning.",
      "stormbreaker_missing_glmnet"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  inner_folds <- sb_check_positive_integer(
    inner_folds,
    "inner_folds",
    minimum = 3L
  )
  seed <- sb_check_positive_integer(seed, "seed", minimum = 1L)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha < 0 || alpha > 1) {
    sb_abort(
      "`alpha` must be one finite number between 0 and 1.",
      "stormbreaker_bad_glmnet_alpha"
    )
  }
  if (!is.logical(standardize) || length(standardize) != 1L ||
      is.na(standardize)) {
    sb_abort(
      "`standardize` must be `TRUE` or `FALSE`.",
      "stormbreaker_bad_glmnet_standardize"
    )
  }
  selection_rule <- sb_validate_glmnet_lambda_rule(selection_rule)
  measure <- sb_validate_glmnet_tuning_measure(measure, "binary")
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "Inner glmnet logistic tuning requires both outcome classes.",
      "stormbreaker_logistic_single_class"
    )
  }

  fold_id <- sb_make_stratified_fold_id(
    outcome,
    inner_folds,
    seed,
    stratum_label = "outcome",
    group_ids = group_ids
  )
  uses_fold_local_selection <- !is.null(selector) ||
    !is.null(fitted_selector) || !is.null(patient_ids)
  if (uses_fold_local_selection &&
      (is.null(selector) || is.null(fitted_selector) ||
       is.null(patient_ids))) {
    sb_abort(
      paste0(
        "Fold-local glmnet selection requires `selector`, ",
        "`fitted_selector`, and `patient_ids` together."
      ),
      "stormbreaker_bad_glmnet_selection_context"
    )
  }
  x_matrix <- as.matrix(x)
  fold_filter <- NULL
  if (uses_fold_local_selection) {
    rownames(x_matrix) <- as.character(patient_ids)
    fold_filter <- sb_make_glmnet_fold_selector_filter(
      selector,
      fitted_selector,
      patient_ids,
      fold_id,
      task_type = "binary", group_ids = group_ids
    )
  }
  cv_fit <- tryCatch(
    withCallingHandlers(
      if (uses_fold_local_selection) {
        glmnet::cv.glmnet(
          x = x_matrix,
          y = outcome,
          family = "binomial",
          alpha = as.numeric(alpha),
          foldid = fold_id,
          type.measure = measure,
          standardize = standardize,
          exclude = fold_filter$exclude,
          alignment = "lambda",
          parallel = FALSE
        )
      } else {
        glmnet::cv.glmnet(
          x = x_matrix,
          y = outcome,
          family = "binomial",
          alpha = as.numeric(alpha),
          foldid = fold_id,
          type.measure = measure,
          standardize = standardize
        )
      },
      warning = function(warning) {
        sb_abort(
          paste0(
            "Inner glmnet logistic tuning raised a warning: ",
            conditionMessage(warning)
          ),
          "stormbreaker_glmnet_logistic_tuning_warning",
          parent = warning
        )
      }
    ),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "Inner glmnet logistic tuning failed: ",
          conditionMessage(error)
        ),
        "stormbreaker_glmnet_logistic_tuning_failure",
        parent = error
      )
    }
  )
  inner_selections <- if (uses_fold_local_selection) {
    fold_filter$inner_selections()
  } else {
    list()
  }

  curve <- data.frame(
    lambda = cv_fit$lambda,
    mean_value = cv_fit$cvm,
    standard_error = cv_fit$cvsd,
    nonzero = cv_fit$nzero,
    stringsAsFactors = FALSE
  )
  names(curve)[[2L]] <- if (identical(measure, "deviance")) {
    "mean_deviance"
  } else {
    "mean_measure"
  }
  valid_curve <- nrow(curve) > 0L &&
    all(is.finite(curve$lambda)) && all(curve$lambda > 0) &&
    all(is.finite(cv_fit$cvm)) &&
    all(is.finite(curve$standard_error)) &&
    all(is.finite(curve$nonzero))
  if (!valid_curve) {
    sb_abort(
      "Inner glmnet logistic tuning produced an incomplete error curve.",
      "stormbreaker_invalid_glmnet_logistic_tuning"
    )
  }
  selected_lambda <- unname(cv_fit[[selection_rule]])
  if (length(selected_lambda) != 1L || !is.finite(selected_lambda) ||
      selected_lambda <= 0 || !selected_lambda %in% curve$lambda) {
    sb_abort(
      "Inner glmnet logistic tuning produced an invalid selected lambda.",
      "stormbreaker_invalid_glmnet_logistic_tuning"
    )
  }

  list(
    lambda = selected_lambda,
    selection_rule = selection_rule,
    measure = measure,
    alpha = as.numeric(alpha),
    standardize = standardize,
    inner_folds = inner_folds,
    seed = seed,
    fold_id = fold_id,
    alignment = "lambda",
    inner_selection_method = if (uses_fold_local_selection) {
      "glmnet_exclude_function"
    } else {
      "none"
    },
    inner_selections = inner_selections,
    curve = curve,
    coefficient_path = if (capture_path) list(lambda = cv_fit$glmnet.fit$lambda,
      beta = cv_fit$glmnet.fit$beta) else NULL,
    warnings = character()
  )
}
