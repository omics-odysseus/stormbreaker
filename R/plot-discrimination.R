# ---- plot-discrimination ----

# ---- paper-roc ----

sb_roc_interval <- function(y, p, reps = 0L) {
  auc <- function(y, p) {
    n1 <- sum(y == 1)
    n0 <- sum(y == 0)
    if (!n1 || !n0) return(NA_real_)
    (sum(rank(p)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  values <- if (reps > 0) replicate(reps, {
    i <- sample.int(length(y), replace = TRUE)
    auc(y[i], p[i])
  }) else numeric()
  ci <- if (sum(is.finite(values)) >= 2) stats::quantile(values, c(.025, .975), na.rm = TRUE) else c(NA_real_, NA_real_)
  data.frame(auc = auc(y, p), lower = unname(ci[1]), upper = unname(ci[2]),
    bootstrap_success = sum(is.finite(values)), bootstrap_requested = reps)
}

#' @rdname sb_plot_performance
#' @export
sb_plot_prediction_roc <- function(outcome, probability, kind, cohort = "Cohort",
                                    ci_reps = 1000L, seed = 1L) {
  kind <- match.arg(kind, c("apparent", "external", "oof"))
  sb_require_curve_plotting()
  if (is.numeric(probability) && is.null(dim(probability))) probability <- matrix(probability, ncol = 1, dimnames = list(NULL, "Model"))
  if (!is.matrix(probability) || !is.numeric(probability) || !ncol(probability) ||
      nrow(probability) != length(outcome) || any(!is.finite(probability)) ||
      !is.numeric(outcome) || anyNA(outcome) || !all(outcome %in% c(0, 1)) ||
      length(unique(outcome)) != 2L) stop("Supply aligned finite scores and both 0/1 classes; no rows are dropped.")
  if (is.null(colnames(probability)) || anyDuplicated(colnames(probability)) || any(!nzchar(colnames(probability)))) stop("Score columns need unique model names.")
  if (length(ci_reps) != 1L || !is.finite(ci_reps) || ci_reps < 0 || ci_reps != as.integer(ci_reps) || ci_reps == 1) stop("ci_reps must be 0 or an integer of at least 2.")
  withr::local_seed(seed)
  tables <- lapply(seq_len(ncol(probability)), function(i) {
    z <- sb_binary_threshold_coordinates(outcome, probability[, i], "roc")
    z$model <- colnames(probability)[i]
    z
  })
  a <- do.call(rbind, lapply(seq_len(ncol(probability)), function(i) sb_roc_interval(outcome, probability[, i], ci_reps)))
  a$model <- colnames(probability)
  labels <- paste0(a$model, ": AUC ", sprintf("%.3f", a$auc),
    if (ci_reps > 0) paste0(" (95% CI ", sprintf("%.3f", a$lower), "\u2013", sprintf("%.3f", a$upper), ")") else "")
  x <- y <- model <- NULL
  p <- ggplot2::ggplot(do.call(rbind, tables), ggplot2::aes(x, y, colour = model)) +
    ggplot2::geom_path(linewidth = .7) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey65") +
    ggplot2::scale_colour_discrete(breaks = a$model, labels = labels) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(title = cohort, x = "1 \u2212 Specificity", y = "Sensitivity", colour = NULL) +
    sb_viz_theme() + ggplot2::theme(legend.position = "bottom")
  attr(p, "stormbreaker_roc_statistics") <- a
  attr(p, "prediction_origin") <- kind
  p
}

# ---- plot-confusion ----

#' @rdname sb_plot_performance
#' @export
sb_plot_confusion <- function(x, threshold, normalize = c("none", "observed"), ncol = 2L) {
  sb_require_curve_plotting()
  normalize <- match.arg(normalize)
  ncol <- sb_viz_positive_integer(ncol, "ncol")
  if (!inherits(x, "stormbreaker_threshold_metrics") || !is.numeric(threshold) ||
      length(threshold) != 1L || !is.finite(threshold) || threshold < 0 || threshold > 1) {
    sb_abort("Supply threshold results and one explicit threshold in [0, 1].",
      "stormbreaker_bad_visualization_argument")
  }
  if (!nrow(x$table)) {
    p <- sb_empty_display_plot("Confusion matrix", "No usable binary predictions.",
      "Unavailable groups remain in diagnostics.")
  } else {
    if (!threshold %in% x$table$threshold) sb_abort(
      "This threshold was not evaluated; call sb_threshold_metrics first.",
      "stormbreaker_bad_visualization_argument")
    tab <- x$table[x$table$threshold == threshold, , drop = FALSE]
    tab$panel_key <- interaction(tab$view_id, tab$pipeline_id,
      addNA(factor(tab$fold)), drop = TRUE)
    data <- do.call(rbind, lapply(seq_len(nrow(tab)), function(i) {
      r <- tab[i, , drop = FALSE]
      d <- r[rep(1L, 4L), , drop = FALSE]
      d$observed <- factor(c(0, 0, 1, 1), levels = c(0, 1))
      d$predicted <- factor(c(0, 1, 0, 1), levels = c(0, 1))
      d$count <- c(r$tn, r$fp, r$fn, r$tp)
      denominator <- c(rep(r$tn + r$fp, 2), rep(r$fn + r$tp, 2))
      d$value <- if (normalize == "none") d$count else
        ifelse(denominator == 0, NA_real_, d$count / denominator)
      d$label <- if (normalize == "none") as.character(d$count) else
        ifelse(is.na(d$value), paste0(d$count, "\nUndefined"),
          paste0(d$count, "\n", sprintf("%.1f%%", 100 * d$value)))
      d
    }))
    predicted <- observed <- value <- label <- panel_key <- NULL
    p <- ggplot2::ggplot(data, ggplot2::aes(predicted, observed, fill = value)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 3.5) +
      ggplot2::scale_fill_gradient(low = "#F2F6F9", high = "#80B1CF", na.value = "grey90") +
      ggplot2::facet_wrap(~panel_key, ncol = ncol,
        labeller = ggplot2::as_labeller(stats::setNames(tab$panel, tab$panel_key))) +
      ggplot2::labs(title = "Which patients are classified correctly?",
        subtitle = paste("Explicit threshold:", threshold), x = "Predicted class", y = "Observed class",
        fill = if (normalize == "none") "Count" else "Within observed class",
        caption = "Class 1 is positive; probability >= threshold. Separate folds/cohorts. No threshold selected.") +
      sb_viz_theme()
    attr(p, "stormbreaker_suggested_size") <- c(width = 11,
      height = max(4, 3 * ceiling(nlevels(tab$panel_key) / ncol)))
  }
  attr(p, "stormbreaker_threshold_diagnostics") <- x$diagnostics
  p
}

# ---- plot-discrimination ----

#' @rdname sb_plot_performance
#' @export
sb_plot_roc <- function(x, cohort_label = NULL, ncol = 2L,
                        by_fold = FALSE, annotate = TRUE, ci_reps = 0L, seed = 1L) {
  sb_plot_binary_discrimination(x, "roc", cohort_label, ncol,
    by_fold, annotate, ci_reps, seed)
}

#' @rdname sb_plot_performance
#' @export
sb_plot_pr <- function(x, cohort_label = NULL, ncol = 2L) {
  sb_plot_binary_discrimination(x, "pr", cohort_label, ncol)
}

sb_binary_threshold_coordinates <- function(outcome, probability, kind) {
  ord <- order(probability, decreasing = TRUE)
  prob <- probability[ord]
  y <- outcome[ord]
  end <- c(which(diff(prob) != 0), length(prob))
  tp <- c(0, cumsum(y == 1)[end])
  fp <- c(0, cumsum(y == 0)[end])
  positive <- sum(y == 1)
  negative <- sum(y == 0)
  data.frame(threshold = c(Inf, prob[end]), true_positive = tp, false_positive = fp,
    x = if (kind == "roc") fp / negative else tp / positive,
    y = if (kind == "roc") tp / positive else c(1, utils::tail(tp, -1)/(utils::tail(tp, -1) + utils::tail(fp, -1))))
}

sb_plot_binary_discrimination <- function(x, kind, cohort_label, ncol,
                                         by_fold = TRUE, annotate = FALSE,
                                         ci_reps = 0L, seed = 1L) {
  sb_require_curve_plotting()
  ncol <- sb_viz_positive_integer(ncol, "ncol")
  leaves <- sb_viz_leaves(x, cohort_label)
  coordinates <- diagnostics <- saved <- annotations <- list()
  if (!is.logical(by_fold) || length(by_fold) != 1L || is.na(by_fold) ||
      !is.logical(annotate) || length(annotate) != 1L || is.na(annotate)) stop("by_fold and annotate must be TRUE or FALSE.")
  if (!is.numeric(ci_reps) || length(ci_reps) != 1L || !is.finite(ci_reps) ||
      ci_reps < 0 || ci_reps != as.integer(ci_reps) || ci_reps == 1) stop("ci_reps must be 0 or an integer of at least 2.")
  withr::local_seed(seed)
  for (leaf in leaves) {
    value <- leaf$value
    external <- inherits(value, "stormbreaker_external_evaluation")
    binary <- if (external) identical(value$task_type, "binary") else
      inherits(value, "stormbreaker_binary_cv_result")
    if (!binary) sb_abort("ROC/PR requires binary saved predictions, not survival scores.",
      "stormbreaker_bad_visualization_task")
    id <- sb_viz_pipeline_id(value)
    pred <- if (external) value$scored_predictions else value$predictions
    metric <- if (kind == "roc") "roc_auc" else "average_precision"
    saved[[leaf$context$view_id]] <- value$metrics[[metric]]
    groups <- if (by_fold && !external && is.data.frame(pred) && nrow(pred) && "fold" %in% names(pred)) {
      unique(pred$fold)
    } else NA_integer_
    for (fold in groups) {
      group <- if (!is.data.frame(pred)) data.frame() else if (external || is.na(fold)) pred else
        pred[!is.na(pred$fold) & pred$fold == fold, , drop = FALSE]
      panel <- if (identical(leaf$context$validation_scheme, "internal_holdout")) {
        sb_viz_panel(leaf)
      } else paste(sb_viz_panel(leaf), if (external) "External patients" else
        if (by_fold) paste("Assessment fold", fold) else "Out-of-fold", sep = "\n")
      status <- if (!identical(value$status, "success")) "failed_run" else if (!nrow(group)) "no_predictions" else "success"
      if (status == "success") {
        columns <- c("outcome", "probability", "status", "patient_id")
        valid <- all(columns %in% names(group)) && is.numeric(group$outcome) &&
          is.numeric(group$probability) && all(is.finite(group$outcome)) &&
          all(group$outcome %in% c(0, 1)) && all(is.finite(group$probability)) &&
          all(group$probability >= 0 & group$probability <= 1) &&
          !anyNA(group$status) && all(group$status == "success") &&
          !anyNA(group$patient_id) && !anyDuplicated(group$patient_id) &&
          (external || ("fold" %in% names(group) && !anyNA(group$fold)))
        if (!valid) status <- "invalid_predictions" else if (!sum(group$outcome == 1) ||
            (kind == "roc" && !sum(group$outcome == 0))) status <- "single_class"
      }
      diagnostic <- data.frame(pipeline_id = id, dataset = leaf$dataset, fold = as.integer(fold),
        panel = panel, curve_status = status, run_status = value$status,
        n = nrow(group), stringsAsFactors = FALSE)
      diagnostic <- sb_view_annotate(diagnostic, leaf$context)
      diagnostics[[length(diagnostics) + 1L]] <- diagnostic
      if (status != "success") next
      tab <- sb_binary_threshold_coordinates(group$outcome, group$probability, kind)
      for (name in names(diagnostic)) tab[[name]] <- rep(diagnostic[[name]], nrow(tab))
      tab$series_id <- paste(leaf$context$view_id, fold, sep = "/")
      coordinates[[length(coordinates) + 1L]] <- tab
      if (annotate && kind == "roc") {
        a <- sb_roc_interval(group$outcome, group$probability, ci_reps)
        a$pipeline_id <- id
        a$panel <- panel
        annotations[[length(annotations) + 1L]] <- a
      }
    }
  }
  status <- sb_view_bind(diagnostics)
  caption <- paste0(if (by_fold) "Separate assessment folds/external cohorts; not pooled. " else
    "One pooled out-of-fold curve per model and resampling context; not apparent training performance. ", "Saved aggregate ",
    if (kind == "roc") "AUC" else "AP", " is unchanged. ",
    sum(status$curve_status != "success"), " unavailable curve group(s). ",
    "AUC labels describe the displayed curve. Patient-bootstrap intervals condition on saved predictions, not development uncertainty.")
  title <- if (kind == "roc") "ROC" else "Precision-recall"
  if (!length(coordinates)) {
    p <- sb_empty_display_plot(title, "No usable saved binary predictions for these curves.", caption)
    p$data <- data.frame()
  } else {
    data <- sb_view_bind(coordinates)
    data$panel <- factor(data$panel, levels = unique(status$panel[status$curve_status == "success"]))
    data$pipeline_id <- factor(data$pipeline_id, levels = unique(status$pipeline_id))
    x <- y <- pipeline_id <- series_id <- panel <- NULL
    p <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = y, colour = pipeline_id, group = series_id))
    p <- p + if (kind == "roc") ggplot2::geom_path(linewidth = 0.65) else
      ggplot2::geom_step(direction = "vh", linewidth = 0.65)
    if (kind == "roc") p <- p + ggplot2::geom_abline(slope = 1, intercept = 0,
      colour = "grey65", linetype = 2, inherit.aes = FALSE)
    p <- p + ggplot2::facet_wrap(~panel, ncol = ncol, labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
      ggplot2::scale_x_continuous(limits = c(0, 1)) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = title, x = if (kind == "roc") "False positive rate" else "Recall",
        y = if (kind == "roc") "Sensitivity" else "Precision", colour = "Pipeline", caption = caption) +
      sb_viz_theme()
    if (length(annotations)) {
      a <- sb_view_bind(annotations)
      a$label <- paste0("AUC ", sprintf("%.3f", a$auc),
        if (ci_reps > 0) paste0(" (95% CI ", sprintf("%.3f", a$lower), "\u2013", sprintf("%.3f", a$upper), ")") else "")
      a$label_y <- .04 + .055 * (stats::ave(seq_len(nrow(a)), a$panel, FUN = seq_along) - 1)
      label <- label_y <- NULL
      p <- p + ggplot2::geom_text(data = a,
        ggplot2::aes(x = .98, y = label_y, label = label, colour = pipeline_id),
        inherit.aes = FALSE, show.legend = FALSE, hjust = 1, size = 3)
      attr(p, "stormbreaker_roc_statistics") <- a
    }
  }
  attr(p, "stormbreaker_curve_status") <- status
  attr(p, "stormbreaker_saved_metrics") <- saved
  if (kind == "roc" && nrow(p$data)) {
    p <- sb_clinician_clean_plot(p, sb_clinician_models(leaves))
    attr(p, "stormbreaker_figure_notes") <- caption
    p <- p + ggplot2::labs(subtitle = NULL, caption = NULL)
  }
  p
}

# ---- plot-prediction-confusion ----

#' @rdname sb_plot_performance
#' @export
sb_plot_prediction_confusion <- function(outcome, probability, kind,
                                          threshold = .5, cohort = NULL,
                                          class_labels = c("0", "1"),
                                          normalize = c("none", "observed")) {
  kind <- match.arg(kind, c("apparent", "external", "oof"))
  normalize <- match.arg(normalize)
  if (!is.numeric(outcome) || !is.null(dim(outcome)) || !length(outcome) ||
      anyNA(outcome) || !all(outcome %in% c(0, 1)) ||
      !is.numeric(probability) || !is.null(dim(probability)) ||
      length(probability) != length(outcome) || any(!is.finite(probability)) ||
      any(probability < 0 | probability > 1)) stop("Supply aligned 0/1 labels and finite probabilities in [0,1].")
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) || threshold < 0 || threshold > 1) stop("threshold must be in [0,1].")
  if (!is.character(class_labels) || length(class_labels) != 2L || anyNA(class_labels) ||
      any(!nzchar(trimws(class_labels))) || anyDuplicated(class_labels)) stop("Supply two distinct class labels for 0 and 1.")
  if (!is.null(cohort) && (!is.character(cohort) || length(cohort) != 1L || is.na(cohort) || !nzchar(trimws(cohort)))) stop("cohort must be a short nonempty label or NULL.")
  positive <- probability >= threshold
  counts <- data.frame(threshold = threshold,
    tp = sum(positive & outcome == 1), fp = sum(positive & outcome == 0),
    tn = sum(!positive & outcome == 0), fn = sum(!positive & outcome == 1))
  table <- cbind(counts, view_id = "predictions", pipeline_id = "Model",
    fold = NA_integer_, panel = "Confusion matrix")
  metrics <- structure(list(table = table, diagnostics = data.frame(), saved_metrics = list()),
    class = "stormbreaker_threshold_metrics")
  p <- sb_plot_confusion(metrics, threshold, normalize = normalize)
  observed <- predicted <- value <- NULL
  p$mapping <- ggplot2::aes(observed, predicted, fill = value)
  labels <- stats::setNames(class_labels, c("0", "1"))
  p <- p + ggplot2::scale_x_discrete(limits = c("0", "1"), labels = labels, drop = FALSE) +
    ggplot2::scale_y_discrete(limits = c("1", "0"), labels = labels, drop = FALSE) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = if (is.null(cohort)) "Confusion matrix" else cohort,
      subtitle = NULL, caption = NULL, x = "Observed", y = "Predicted") +
    ggplot2::theme(strip.text = ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
      legend.position = "none")
  attr(p, "stormbreaker_prediction_context") <- list(kind = kind, cohort = cohort,
    threshold = threshold, n = length(outcome), class_labels = class_labels, counts = counts)
  attr(p, "stormbreaker_suggested_size") <- c(width = 4, height = 4)
  p
}

# ---- survival-views ----

# ---- survival-cohort-summary ----

sb_risk_cox_estimate <- function(d, type = "group") {
  tryCatch({
    value <- if (type == "group") as.numeric(d$group == "High") else d$risk
    if (length(unique(value)) < 2L || !any(d$event == 1)) stop("Insufficient groups, variation or events.")
    model <- withCallingHandlers(survival::coxph(survival::Surv(time, event) ~ value, data = d),
      warning = function(w) stop(conditionMessage(w)))
    te <- unname(stats::coef(model)[1]); se <- sqrt(stats::vcov(model)[1, 1])
    if (!is.finite(te) || !is.finite(se) || se <= 0) stop("Cox estimate is not finite.")
    data.frame(log_hr = te, se = se, hr = exp(te), lower = exp(te - stats::qnorm(.975) * se),
      upper = exp(te + stats::qnorm(.975) * se), p = 2 * stats::pnorm(-abs(te / se)),
      status = "success", reason = "")
  }, error = function(e) data.frame(log_hr = NA_real_, se = NA_real_, hr = NA_real_,
    lower = NA_real_, upper = NA_real_, p = NA_real_, status = "not_evaluable", reason = conditionMessage(e)))
}

#' @rdname sb_plot_performance
#' @export
sb_plot_cohort_metrics <- function(performance, metric = c("c_index", "auc"), title = NULL) {
  metric <- match.arg(metric)
  sb_require_curve_plotting()
  if (!is.data.frame(performance) || !all(c("cohort", metric) %in% names(performance)) ||
      !nrow(performance) || anyNA(performance$cohort) || anyDuplicated(performance$cohort) ||
      !is.numeric(performance[[metric]]) || any(!is.finite(performance[[metric]])) ||
      any(performance[[metric]] < 0 | performance[[metric]] > 1)) stop("Supply unique cohorts and finite metric values in [0,1].")
  d <- data.frame(cohort = factor(performance$cohort, levels = rev(performance$cohort)), value = performance[[metric]])
  cohort <- value <- NULL
  ggplot2::ggplot(d, ggplot2::aes(value, cohort, fill = cohort)) +
    ggplot2::geom_col(width = .72, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", value)), hjust = -.1, size = 3.5) +
    ggplot2::scale_x_continuous(limits = c(0, 1.08), breaks = seq(0, 1, .2)) +
    ggplot2::labs(x = if (metric == "c_index") "C-index" else "AUC", y = NULL, title = title) + sb_viz_theme()
}

#' @rdname sb_meta_external_discrimination
#' @export
surv_meta_risk <- function(cohorts, type = c("group", "continuous")) {
  type <- match.arg(type)
  sb_require_curve_plotting()
  if (!is.list(cohorts) || length(cohorts) < 2L || is.null(names(cohorts)) ||
      anyNA(names(cohorts)) || any(!nzchar(names(cohorts))) || anyDuplicated(names(cohorts))) stop("Supply at least two named independent cohorts.")
  required <- c("patient_id", "risk", "time", "event", "group")
  for (z in cohorts) if (!is.data.frame(z$patients) || !all(required %in% names(z$patients))) stop("Use surv_risk_plots results.")
  ids <- unlist(lapply(cohorts, function(z) as.character(z$patients$patient_id)), use.names = FALSE)
  if (anyNA(ids) || anyDuplicated(ids)) stop("Patient IDs overlap: do not pool combined and constituent cohorts.")
  if (type == "group" && length(unique(vapply(cohorts, function(z) z$cutoff, numeric(1)))) != 1L) stop("Use the same training cutoff in every cohort.")
  estimates <- do.call(rbind, lapply(seq_along(cohorts), function(i) {
    cbind(cohort = names(cohorts)[i], sb_risk_cox_estimate(cohorts[[i]]$patients, type))
  }))
  if (any(estimates$status != "success")) return(list(estimates = estimates,
    status = "not_evaluable", plot = NULL, reason = "Review failed cohort estimates before pooling."))
  if (!requireNamespace("metafor", quietly = TRUE)) stop("Install metafor for Cox meta-analysis.")
  fixed <- metafor::rma.uni(yi = estimates$log_hr, sei = estimates$se, method = "FE")
  random <- metafor::rma.uni(yi = estimates$log_hr, sei = estimates$se, method = "DL")
  estimates$weight_fixed <- as.numeric(stats::weights(fixed))
  estimates$weight_random <- as.numeric(stats::weights(random))
  pooled <- do.call(rbind, lapply(list(Fixed = fixed, Random = random), function(m) {
    data.frame(log_hr = as.numeric(m$b), se = m$se, hr = exp(as.numeric(m$b)),
      lower = exp(m$ci.lb), upper = exp(m$ci.ub), p = m$pval)
  }))
  pooled$cohort <- rownames(pooled)
  d <- rbind(estimates[c("cohort", "hr", "lower", "upper")], pooled[c("cohort", "hr", "lower", "upper")])
  d$cohort <- factor(d$cohort, levels = rev(d$cohort))
  cohort <- hr <- lower <- upper <- NULL
  plot <- ggplot2::ggplot(d, ggplot2::aes(hr, cohort)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper), width = .15, orientation = "y") +
    ggplot2::geom_point(size = 2.5) + ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Hazard ratio (95% CI)", y = NULL) + sb_viz_theme()
  list(estimates = estimates, pooled = pooled, plot = plot, type = type, status = "success",
    heterogeneity = data.frame(Q = random$QE, p = random$QEp, I2 = random$I2, tau2 = random$tau2))
}

# ---- survival-metric-view ----

sb_horizon_view_metrics <- function() c("time_auc", "uno_c", "ipcw_brier")

sb_scalar_view_metrics <- function() {
  c("roc_auc", "average_precision", "brier_score", "log_loss", "harrell_c",
    sb_horizon_view_metrics())
}

sb_horizon_metric_row <- function(result, metric, task_type, external = FALSE,
                                 cohort_label = NULL) {
  id <- if (external) result$pipeline_id else result$pipeline$id
  if (!identical(task_type, "survival")) {
    sb_abort(paste0("Metric `", metric, "` is not saved in pipeline `", id,
      "`. No metric is recomputed by this function."),
      "stormbreaker_result_metric_not_saved")
  }
  value_or <- function(value, default) if (is.null(value)) default else value
  descriptor <- if (external) result$endpoint else result$task
  value <- result$metrics[[metric]]
  horizon <- value_or(value$prediction_horizon,
    value_or(descriptor$prediction_horizon, NA_real_))
  unit <- value_or(value$time_unit, value_or(descriptor$time_unit, NA_character_))
  status <- if (!is.null(value)) value$status else if (
    identical(result$status, "success") && is.na(horizon)) "not_requested" else
      "not_saved"
  output <- data.frame(
    pipeline_id = id, task_type = task_type, metric = metric,
    estimate = value_or(value$estimate, NA_real_), status = status,
    run_status = result$status, scope = value_or(value$scope, NA_character_),
    aggregation = value_or(value$aggregation, NA_character_),
    prediction_horizon = as.numeric(horizon), time_unit = unit,
    metric_error_class = value_or(value$error_class, NA_character_),
    metric_error_message = value_or(value$error_message, NA_character_),
    stringsAsFactors = FALSE
  )
  if (external) {
    output$cohort_label <- sb_result_cohort_label(cohort_label)
    output$patient_match_status <- result$patient_match$status
    output$error_stage <- result$error_stage
    output$error_class <- result$error_class
    output$error_message <- result$error_message
  }
  output
}

sb_horizon_view_labels <- function(tab) {
  horizon <- vapply(tab$prediction_horizon, function(value) {
    if (is.na(value)) "horizon not saved" else format(value, digits = 17, trim = TRUE)
  }, character(1))
  unit <- tab$time_unit
  unit[is.na(unit)] <- "unit not saved"
  paste0(as.character(tab$pipeline_id), " [", horizon, " ", unit, "]")
}

# ---- survival-risk-panels ----

sb_risk_plot_task <- function(task) {
  task <- sb_validate_task_contract(task)
  if (!inherits(task, "stormbreaker_surv_task") || !is.null(task$group_col)) {
    stop("Supply a survival task with one independent row per patient.")
  }
  task
}

sb_risk_plot_scores <- function(risk, task) {
  ids <- as.character(task$data[[task$id_col]])
  if (!is.numeric(risk) || !is.null(dim(risk)) || length(risk) != length(ids) || any(!is.finite(risk))) {
    stop("Supply one finite numeric risk score per patient, higher means worse prognosis.")
  }
  if (!is.null(names(risk))) {
    if (anyDuplicated(names(risk)) || !setequal(names(risk), ids)) stop("Named risk IDs must match task patients.")
    risk <- risk[match(ids, names(risk))]
  }
  unname(risk)
}

sb_risk_predict <- function(fit, task) {
  p <- sb_predict_final(fit, task$data, task$id_col)
  if (p$status != "success" || !"risk" %in% names(p$predictions)) stop("Successful survival risk predictions are required.")
  risk <- stats::setNames(p$predictions$risk, p$predictions$patient_id)
  sb_risk_plot_scores(risk, task)
}

#' @rdname sb_survival_curves
#' @export
surv_risk_reference <- function(fit, development) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  task <- sb_risk_plot_task(development)
  ids <- as.character(task$data[[task$id_col]])
  if (fit$pipeline$task_type != "survival" || !setequal(ids, fit$development$patient_ids)) {
    stop("Supply the original development cohort for this survival fit.")
  }
  risk <- sb_risk_predict(fit, task)
  genes <- fit$selection$selected_features
  x <- task$data[genes]
  scale <- vapply(x, stats::sd, numeric(1))
  scale[scale == 0] <- 1
  structure(list(fit = fit, cutoff = stats::median(risk), genes = genes,
    center = vapply(x, mean, numeric(1)), scale = scale, time_unit = task$time_unit),
    class = "stormbreaker_risk_reference")
}

#' @rdname sb_survival_curves
#' @export
surv_risk_plots <- function(reference, task, genes = NULL) {
  sb_require_curve_plotting()
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install patchwork.")
  if (!inherits(reference, "stormbreaker_risk_reference")) stop("Use surv_risk_reference first.")
  task <- sb_risk_plot_task(task)
  if (!identical(task$time_unit, reference$time_unit)) stop("Time units must match the development reference.")
  if (is.null(genes)) genes <- reference$genes
  if (!is.character(genes) || !length(genes) || anyNA(genes) || anyDuplicated(genes) || !all(genes %in% reference$genes)) stop("Choose unique selected genes.")
  risk <- sb_risk_predict(reference$fit, task)
  d <- data.frame(patient_id = as.character(task$data[[task$id_col]]), risk = risk,
    time = task$data[[task$time_col]], event = task$data[[task$event_col]],
    group = factor(ifelse(risk > reference$cutoff, "High", "Low"), levels = c("Low", "High")))
  ord <- order(d$risk, seq_len(nrow(d)))
  d <- d[ord, ]; d$position <- seq_len(nrow(d))
  d$outcome <- factor(ifelse(d$event == 1, "Event", "Censored"), levels = c("Censored", "Event"))
  z <- sweep(as.matrix(task$data[ord, genes, drop = FALSE]), 2, reference$center[genes], "-")
  z <- sweep(z, 2, reference$scale[genes], "/")
  rownames(z) <- d$patient_id
  heat <- data.frame(position = rep(d$position, times = length(genes)),
    gene = factor(rep(genes, each = nrow(d)), levels = rev(genes)), value = as.vector(z))
  position <- gene <- value <- group <- time <- outcome <- survival <- n_risk <- lower <- upper <- NULL
  colours <- c(Low = "#2878B5", High = "#D45064")
  limits <- c(.5, nrow(d) + .5)
  common_x <- ggplot2::scale_x_continuous(limits = limits, expand = c(0, 0))
  plots <- list(
    expression = ggplot2::ggplot(heat, ggplot2::aes(position, gene, fill = value)) +
      ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#2878B5", mid = "white", high = "#D45064", name = "Z-score") +
      common_x + ggplot2::labs(x = NULL, y = NULL) + sb_viz_theme(),
    risk = ggplot2::ggplot(d, ggplot2::aes(position, risk, colour = group)) +
      ggplot2::geom_point(size = 1.5) + ggplot2::geom_hline(yintercept = reference$cutoff, linetype = 2) +
      common_x + ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
      ggplot2::labs(x = "Patient order", y = "Risk score", colour = "Group") + sb_viz_theme(),
    outcome = ggplot2::ggplot(d, ggplot2::aes(position, time, colour = outcome)) +
      ggplot2::geom_point(size = 1.5) + common_x +
      ggplot2::scale_colour_manual(values = c(Censored = "#2878B5", Event = "#D45064"), drop = FALSE) +
      ggplot2::labs(x = "Patient order", y = paste0("Follow-up (", task$time_unit, ")"), colour = "Outcome") + sb_viz_theme())
  km <- survival::survfit(survival::Surv(time, event) ~ group, data = d)
  sm <- summary(km, censored = TRUE)
  groups <- levels(droplevels(d$group))
  kmdata <- data.frame(time = sm$time, survival = sm$surv, lower = sm$lower, upper = sm$upper,
    group = if (is.null(sm$strata)) rep(groups, length(sm$time)) else sub("^group=", "", as.character(sm$strata)))
  kmdata <- rbind(data.frame(time = 0, survival = 1, lower = 1, upper = 1, group = groups), kmdata)
  kmdata <- kmdata[order(match(kmdata$group, groups), kmdata$time), ]
  ticks <- pretty(c(0, max(d$time)), n = 5)
  ticks <- sort(unique(c(0, ticks[ticks >= 0 & ticks <= max(d$time)])))
  at_risk <- expand.grid(time = ticks, group = groups, stringsAsFactors = FALSE)
  at_risk$n_risk <- mapply(function(t, g) sum(d$time >= t & d$group == g), at_risk$time, at_risk$group)
  band <- do.call(rbind, lapply(split(kmdata, kmdata$group), function(z) {
    n <- nrow(z)
    if (n < 2L) return(z)
    out <- z[rep(seq_len(n), each = 2)[- (2 * n)], ]
    out$time <- c(z$time[1], rep(z$time[-1], each = 2))
    out
  }))
  plots$km <- ggplot2::ggplot(kmdata, ggplot2::aes(time, survival, colour = group)) +
    ggplot2::geom_ribbon(data = band, ggplot2::aes(ymin = lower, ymax = upper, fill = group),
      alpha = .18, colour = NA, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::geom_step() + ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_x_continuous(limits = c(0, max(d$time)), breaks = ticks) +
    ggplot2::labs(x = paste0("Time (", task$time_unit, ")"), y = "Survival probability", colour = "Group") + sb_viz_theme()
  plots$at_risk <- ggplot2::ggplot(at_risk, ggplot2::aes(time, group, label = n_risk, colour = group)) +
    ggplot2::geom_text(show.legend = FALSE) + ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_x_continuous(limits = c(0, max(d$time)), breaks = ticks) +
    ggplot2::labs(x = NULL, y = "At risk") + sb_viz_theme()
  logrank <- sb_risk_logrank(d)
  list(plots = plots, plot = patchwork::wrap_plots(plots, ncol = 1, heights = c(2, 1, 1, 2, .6)),
    patients = d, expression = z, km = km, at_risk = at_risk,
    logrank_p = logrank$p, logrank_status = logrank$status, logrank_error = logrank$error,
    cox_group = sb_risk_cox_estimate(d),
    cutoff = reference$cutoff, status = if (length(groups) == 2) "two_groups" else "single_group")
}

sb_risk_logrank <- function(d) {
  if (length(unique(d$group)) != 2L) {
    return(list(p = NA_real_, status = "single_group", error = NA_character_))
  }
  if (!any(d$event == 1)) {
    return(list(p = NA_real_, status = "no_events", error = NA_character_))
  }
  tryCatch({
    test <- survival::survdiff(survival::Surv(time, event) ~ group, data = d)
    p <- stats::pchisq(test$chisq, df = 1, lower.tail = FALSE)
    if (length(p) != 1L || !is.finite(p) || p < 0 || p > 1) {
      stop("The log-rank test did not return a finite P value in [0, 1].")
    }
    list(p = p, status = "success", error = NA_character_)
  }, error = function(e) {
    message <- conditionMessage(e)
    warning(paste0("Log-rank test failed: ", message), call. = FALSE)
    list(p = NA_real_, status = "failed", error = message)
  })
}

# ---- survival-roc-plot ----

#' @rdname sb_survival_curves
#' @export
surv_plot_roc <- function(task, risk, times) {
  sb_require_curve_plotting()
  task <- sb_risk_plot_task(task)
  risk <- sb_risk_plot_scores(risk, task)
  if (!is.numeric(times) || !length(times) || any(!is.finite(times)) || any(times <= 0) || anyDuplicated(times)) stop("Supply unique positive times in the task time unit.")
  time <- task$data[[task$time_col]]
  event <- task$data[[task$event_col]]
  censoring <- sb_fit_censoring_km(time, event)
  curves <- list()
  statistics <- lapply(sort(times), function(t) {
    cases <- event == 1 & time <= t
    controls <- time > t
    row <- data.frame(time = t, auc = NA_real_, n_cases = sum(cases), n_controls = sum(controls),
      n_early_censored = sum(event == 0 & time <= t), status = "not_evaluable", reason = "")
    if (!any(cases) || !any(controls)) {
      row$reason <- "No cases or no controls at this time."
      return(row)
    }
    gcase <- vapply(time[cases], function(s) sb_censoring_survival_at(censoring, s, left_limit = TRUE), numeric(1))
    gcontrol <- sb_censoring_survival_at(censoring, t)
    if (any(!is.finite(c(gcase, gcontrol))) || any(c(gcase, gcontrol) <= 0)) {
      row$reason <- "Censoring survival is zero or unavailable."
      return(row)
    }
    weights_case <- weights_control <- numeric(length(risk))
    weights_case[cases] <- 1 / gcase
    weights_control[controls] <- 1 / gcontrol
    ord <- order(risk, decreasing = TRUE)
    last <- c(which(diff(risk[ord]) != 0), length(risk))
    sensitivity <- c(0, cumsum(weights_case[ord])[last] / sum(weights_case))
    fpr <- c(0, cumsum(weights_control[ord])[last] / sum(weights_control))
    row$auc <- sum(diff(fpr) * (utils::head(sensitivity, -1) + utils::tail(sensitivity, -1)) / 2)
    row$status <- "success"
    curves[[as.character(t)]] <<- data.frame(time = t, threshold = c(Inf, risk[ord][last]),
      fpr = fpr, sensitivity = sensitivity, label = sprintf("%g %s: AUC %.3f", t, task$time_unit, row$auc))
    row
  })
  statistics <- do.call(rbind, statistics)
  coordinates <- if (length(curves)) do.call(rbind, curves) else
    data.frame(time = numeric(), threshold = numeric(), fpr = numeric(), sensitivity = numeric(), label = character())
  fpr <- sensitivity <- label <- NULL
  p <- ggplot2::ggplot(coordinates, ggplot2::aes(fpr, sensitivity, colour = label)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70") +
    ggplot2::geom_path(linewidth = .7) + ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "1 - Specificity", y = "Sensitivity", colour = NULL) + sb_viz_theme()
  list(plot = p, coordinates = coordinates, statistics = statistics,
    method = "cumulative_dynamic_ipcw_marginal_reverse_km", time_unit = task$time_unit)
}
