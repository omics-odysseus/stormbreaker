#' Survival curves, risk views and nomograms
#' @name sb_survival_curves
#' @description Recover supported survival curves and linear-model nomograms, or draw saved risk summaries.
#' @rawRd \arguments{
#' \item{fit}{
#' \strong{For \code{sb_survival_curves()}:} A final fit augmented by \code{\link[=sb_add_survival_curve]{sb_add_survival_curve()}}.
#' 
#' \strong{For \code{sb_add_survival_curve()}:} Successful portable linear Cox final fit with a saved horizon.
#' 
#' \strong{For \code{sb_nomogram()}:} Successful final model supported by \code{\link[=sb_model_coefficients]{sb_model_coefficients()}}.
#' 
#' \strong{For \code{surv_risk_reference()}:} Successful survival final fit.
#' }
#' \item{new_data}{
#' \strong{For \code{sb_survival_curves()}:} New feature data with unique identifiers; other columns ignored.
#' 
#' \strong{For \code{sb_nomogram_values()}:} Original-scale predictors. No outcomes are read.
#' }
#' \item{id_col}{
#' Unique identifier column name.
#' }
#' \item{times}{
#' \strong{For \code{sb_survival_curves()}:} Optional finite nonnegative time points within training follow-up.
#' Sorted unique points are returned. NULL uses the full saved step grid.
#' 
#' \strong{For \code{surv_plot_roc()}:} Positive time points in the task's declared time unit.
#' }
#' \item{development}{
#' \strong{For \code{sb_add_survival_curve()}:} The original, unchanged development survival task.
#' Test or validation outcomes must never be supplied. IDs, sample count,
#' time units and the saved horizon hazard are checked, but those checks
#' cannot prove that all supplied values are the original training values.
#' 
#' \strong{For \code{surv_risk_reference()}:} Original development survival task; one row per patient.
#' }
#' \item{x}{
#' \strong{For \code{sb_plot_survival_curves()}:} A result from \code{\link[=sb_survival_curves]{sb_survival_curves()}}.
#' 
#' \strong{For \code{sb_nomogram_values()}, \code{sb_plot_nomogram()}:} A result from \code{\link[=sb_nomogram]{sb_nomogram()}}.
#' }
#' \item{rows}{
#' Patient row numbers, or integer(0) for baseline only.
#' }
#' \item{baseline}{
#' Whether to include the zero-score baseline.
#' }
#' \item{reference_data}{
#' \strong{For \code{sb_nomogram()}:} Original-scale training predictors defining display
#' ranges. No outcomes are used and no model is refitted.
#' 
#' \strong{For \code{sb_plot_nomogram()}:} Optional original-scale reference observations for
#' density overlays. They must lie within the nomogram's existing ranges.
#' }
#' \item{ticks}{
#' Approximate ticks per feature axis.
#' }
#' \item{show_zero}{
#' Show axes for exactly zero coefficients. By default these
#' contribute no points and are omitted from the display only; their names
#' remain in the plot metadata and the original nomogram table.
#' }
#' \item{patient}{
#' Optional one-row feature frame for exact patient point and
#' probability markers. Neither option fits or changes the model.
#' }
#' \item{task}{
#' \strong{For \code{sb_cox_forest()}:} A survival task with one independent patient per row.
#' 
#' \strong{For \code{surv_risk_plots()}:} Development or validation survival task, with the same time unit.
#' 
#' \strong{For \code{surv_plot_roc()}:} Survival task with one independent patient per row.
#' }
#' \item{conf_level}{
#' Confidence level, default 0.95.
#' }
#' \item{reference}{
#' Result of surv_risk_reference, locked on development data.
#' }
#' \item{genes}{
#' Optional subset of selected genes for the heatmap. NULL uses all.
#' }
#' \item{risk}{
#' Finite numeric prognostic scores, higher means worse prognosis.
#' Unnamed scores follow task row order; named scores are aligned by patient
#' IDs. Use one saved model's score, not separately selected models per time.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_survival_curves()}}{
#' Predict full survival step curves without reading patient outcomes
#' 
#' \strong{Returns:} A \code{stormbreaker_survival_curves} list with baseline table and
#' time-by-patient survival matrix, IDs, horizon and unchanged final prediction.
#' Full step-grid values are retained separately for faithful plotting even
#' when only a few prediction times were requested.
#' No uncertainty intervals are inferred from point estimates.
#' }
#' 
#' \subsection{\code{sb_add_survival_curve()}}{
#' Add a full baseline curve to a saved linear survival model
#' 
#' \strong{Returns:} A copy of fit with \code{survival_curve} metadata and time/hazard table.
#' Existing coefficients, learner state, and horizon predictions are unchanged.
#' 
#' \strong{Details:} Estimates only the missing baseline at fixed saved coefficients,
#' using the saved Breslow or Efron tie method. Unsupported baseline methods
#' are refused. No variable selection, tuning or coefficient fitting occurs.
#' Store the returned model in a new file to preserve the original artifact.
#' }
#' 
#' \subsection{\code{sb_plot_survival_curves()}}{
#' Plot a baseline and selected individual survival curves
#' 
#' \strong{Returns:} A ggplot step curve with the saved-horizon probability markers.
#' Patient row numbers replace original IDs, but predictions remain sensitive.
#' }
#' 
#' \subsection{\code{sb_nomogram()}}{
#' Construct an exact points representation of an existing linear model
#' 
#' \strong{Returns:} A \code{stormbreaker_nomogram} with the unchanged model coefficients,
#' display ranges, point conversion and optional saved survival baseline.
#' Models without an absolute survival baseline retain a risk-score axis.
#' AFT models use their original distribution and scale at the saved horizon,
#' not a Cox baseline. Without a saved horizon they retain their original
#' negative-log-time-location score axis.
#' }
#' 
#' \subsection{\code{sb_nomogram_values()}}{
#' Recover the original prediction from nomogram points
#' 
#' \strong{Returns:} List with per-feature points, total points, original score, and
#' probability when identified. Outside-display-range values are flagged,
#' never clipped; such points cannot be read from the displayed axes.
#' }
#' 
#' \subsection{\code{sb_plot_nomogram()}}{
#' Draw a nomogram without fitting a surrogate model
#' 
#' \strong{Returns:} A ggplot with separate feature, points, total-points and outcome
#' axes. Total points use their own labelled scale, as in a conventional
#' nomogram. Short feature axes stagger endpoint labels to avoid collisions.
#' Instructions are retained in the stormbreaker_notes attribute, not drawn
#' as a subtitle or figure footer.
#' }
#' 
#' \subsection{\code{sb_cox_forest()}}{
#' Fit single-gene Cox models and draw a hazard-ratio forest plot
#' 
#' \strong{Returns:} List with plot, complete results table and confidence level.
#' Failed fits remain in the table and are not plotted. No significance
#' filtering is performed. P values and BH-adjusted P values are returned.
#' 
#' \strong{Details:} Each gene is fitted separately using Efron ties. Hazard ratios
#' describe a one-unit increase on the supplied expression scale, not an
#' adjusted or causal effect. Wald intervals are pointwise, not simultaneous.
#' Grouped tasks are not supported. This analysis does not select genes for
#' a validated pipeline; screening must be repeated within training folds.
#' }
#' 
#' \subsection{\code{surv_risk_plots()}}{
#' Draw linked expression, risk, outcome and Kaplan-Meier panels
#' 
#' \strong{Returns:} List with plots, a combined plot, patient table in risk order,
#' standardized expression matrix, KM estimates, at-risk table, log-rank P
#' and grouping status. \code{logrank_status} is success, single_group, no_events
#' or failed. \code{logrank_error} retains the original error on failure, otherwise
#' it is NA. Failed tests warn without discarding usable plots.
#' KM curves include pointwise 95\% log confidence bands;
#' cox_group contains the High-versus-Low Cox HR, Wald CI and P value.
#' Patient-level results are sensitive; review sharing.
#' 
#' \strong{Details:} Each cohort is sorted by increasing saved-model risk, breaking ties
#' by input row order. Expression, risk and outcome panels share this order.
#' Grouping and gene scaling never use validation outcomes. KM and log-rank
#' describe the supplied cohort; apparent development separation is not
#' independent validation. Log-rank P is returned, not printed as tiny text.
#' }
#' 
#' \subsection{\code{surv_risk_reference()}}{
#' Lock risk grouping and expression scaling from the development cohort
#' 
#' \strong{Returns:} Reference containing the saved fit, training-median risk threshold,
#' selected genes, training expression means/SDs and time unit. No refitting.
#' 
#' \strong{Details:} Use the original unchanged development data. IDs are checked but
#' cannot prove that supplied feature values are unchanged. Scores strictly
#' above the training median are High; all other scores are Low. Constant
#' training genes use scale 1 for display. No outcome-based cutpoint search.
#' }
#' 
#' \subsection{\code{surv_plot_roc()}}{
#' Plot cumulative/dynamic survival ROC at multiple time points
#' 
#' \strong{Returns:} List with plot, coordinates and per-time AUC/status table.
#' 
#' \strong{Details:} Cases have observed events at or before t; controls survive beyond
#' t. Earlier censored patients are neither cases nor controls. IPCW uses
#' the supplied cohort's marginal reverse Kaplan-Meier censoring estimate:
#' cases receive 1/G(T-), controls 1/G(t). Tied scores enter together and AUC
#' uses trapezoidal integration (half credit for ties). This assumes marginal
#' independent censoring. No confidence intervals or automatic model ranking
#' are produced; unsupported horizons remain in the status table. Do not pool
#' incompatible fold-specific score scales. The score is not a probability.
#' }
#' }
#' @rawRd \examples{
#' # sb_nomogram
#' d <- bi_example(40, seed = 42)
#' genes <- c("gene_a", "gene_b")
#' fit <- fit_final(bi_task(d, "patient_id", "outcome", genes),
#'   bi_pipeline("Logistic"))
#' points <- nomogram(fit, d[genes])
#' points$features
#' 
#' # sb_nomogram_values
#' d <- bi_example(40, seed = 42)
#' genes <- c("gene_a", "gene_b")
#' fit <- fit_final(bi_task(d, "patient_id", "outcome", genes),
#'   bi_pipeline("Logistic"))
#' points <- nomogram(fit, d[genes])
#' nom_values(points, d[1:3, genes])
#' 
#' # sb_plot_nomogram
#' \dontshow{if (requireNamespace("ggplot2", quietly = TRUE)) withAutoprint(\{ # examplesIf}
#' d <- bi_example(40, seed = 42)
#' genes <- c("gene_a", "gene_b")
#' fit <- fit_final(bi_task(d, "patient_id", "outcome", genes),
#'   bi_pipeline("Logistic"))
#' points <- nomogram(fit, d[genes])
#' plot_nomogram(points, patient = d[1, genes, drop = FALSE])
#' \dontshow{\}) # examplesIf}
#' 
#' # surv_risk_plots
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE)) {
#'   d <- sb_example_survival(60, seed = 42)
#'   task <- sb_surv_task(d, "patient_id", "time", "event",
#'     c("gene_a", "gene_b"), "month")
#'   fit <- sb_fit_final(task, surv_pipeline("Cox"))
#'   reference <- surv_risk_reference(fit, task)
#'   panels <- surv_risk_plots(reference, task)
#'   panels[c("logrank_p", "logrank_status", "logrank_error")]
#' }
#' }
NULL

# ---- model-curves ----

# ---- model-survival-curves ----

#' @rdname sb_survival_curves
#' @export
sb_add_survival_curve <- function(fit, development) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  state <- fit$learner_state
  baseline <- state$baseline_hazard
  if (!state$learner %in% sb_linear_survival_learner_ids() ||
    is.null(baseline) || !baseline$ties %in% c("breslow", "efron")) {
    sb_abort(
      "A linear Cox model with a saved Breslow/Efron baseline is required.",
      "stormbreaker_unsupported_survival_curve"
    )
  }
  development <- sb_validate_task_contract(development)
  ids <- as.character(development$data[[development$id_col]])
  if (!setequal(ids, fit$development$patient_ids) ||
    length(ids) != length(fit$development$patient_ids) ||
    !identical(development$time_unit, baseline$time_unit)) {
    sb_abort(
      "Supply the original development task with matching IDs and time units.",
      "stormbreaker_curve_development_mismatch"
    )
  }
  x <- sb_select_final_feature_frame(
    development$data, state$features,
    "Development data", "development fitted features"
  )
  risk <- drop(as.matrix(x) %*% state$coefficient)
  time <- development$data[[development$time_col]]
  event <- development$data[[development$event_col]]
  horizon <- baseline$prediction_horizon
  sb_validate_prediction_horizon(horizon, time)
  event_times <- sort(unique(time[event == 1L]))
  increments <- vapply(event_times, function(t) {
    at_risk <- time >= t
    deaths <- time == t & event == 1L
    n_deaths <- sum(deaths)
    maximum <- max(risk[at_risk])
    denominator <- sum(exp(risk[at_risk] - maximum))
    if (baseline$ties == "breslow") {
      return(exp(log(n_deaths) - maximum - log(denominator)))
    }
    death_risk <- sum(exp(risk[deaths] - maximum))
    sum(exp(-maximum - log(denominator - (seq_len(n_deaths) - 1) / n_deaths * death_risk)))
  }, numeric(1))
  hazard <- cumsum(increments)
  if (any(!is.finite(hazard)) || any(increments <= 0)) {
    sb_abort(
      "Recovered baseline hazard is not finite and nondecreasing.",
      "stormbreaker_invalid_baseline_hazard"
    )
  }
  at_horizon <- sum(increments[event_times <= horizon])
  if (!isTRUE(all.equal(at_horizon, baseline$cumulative_baseline_hazard,
    tolerance = 1e-8, check.attributes = FALSE
  ))) {
    sb_abort(
      "Recovered baseline does not match the saved horizon; check original training data and method.",
      "stormbreaker_curve_baseline_mismatch"
    )
  }
  grid <- sort(unique(c(0, event_times, horizon, max(time))))
  values <- c(0, hazard)[findInterval(grid, event_times) + 1L]
  fit$survival_curve <- list(
    table = data.frame(
      time = grid,
      cumulative_hazard = values, baseline_survival = exp(-values)
    ),
    time_unit = baseline$time_unit, prediction_horizon = horizon,
    max_followup = max(time), source = baseline$source, ties = baseline$ties,
    reference = "x_beta = 0; not an average patient",
    development_n = length(ids), method = "original_development_fixed_coefficients"
  )
  fit
}

#' @rdname sb_survival_curves
#' @export
sb_survival_curves <- function(fit, new_data, id_col, times = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  curve <- fit$survival_curve
  if (is.null(curve)) {
    sb_abort(
      "Add the original-training baseline with sb_add_survival_curve first.",
      "stormbreaker_missing_survival_curve"
    )
  }
  ids <- sb_extract_final_prediction_ids(new_data, id_col)
  x <- sb_select_final_prediction_features(new_data, id_col, fit$required_features)
  if (is.null(times)) times <- curve$table$time
  if (!is.numeric(times) || !length(times) || any(!is.finite(times)) ||
    any(times < 0 | times > curve$max_followup)) {
    sb_abort(
      "Times must be finite, nonnegative and within original training follow-up.",
      "stormbreaker_curve_time_unsupported"
    )
  }
  times <- sort(unique(times))
  hazard <- curve$table$cumulative_hazard[findInterval(times, curve$table$time)]
  risk <- drop(as.matrix(x) %*% fit$learner_state$coefficient)
  probability <- vapply(
    hazard, function(h) sb_survival_probability_from_hazard(risk, h),
    numeric(length(ids))
  )
  probability <- matrix(probability, nrow = length(ids), ncol = length(times))
  full_probability <- vapply(
    curve$table$cumulative_hazard,
    function(h) sb_survival_probability_from_hazard(risk, h), numeric(length(ids))
  )
  full_probability <- matrix(full_probability, nrow = length(ids), ncol = nrow(curve$table))
  structure(
    list(
      baseline = data.frame(
        time = times, cumulative_hazard = hazard,
        survival = exp(-hazard)
      ), survival = t(probability), patient_ids = ids,
      plot_baseline = data.frame(
        time = curve$table$time,
        survival = curve$table$baseline_survival
      ), plot_survival = t(full_probability),
      time_unit = curve$time_unit, prediction_horizon = curve$prediction_horizon,
      max_followup = curve$max_followup, reference = curve$reference,
      prediction = sb_predict_final(fit, new_data, id_col)
    ),
    class = "stormbreaker_survival_curves"
  )
}

#' @rdname sb_survival_curves
#' @export
sb_plot_survival_curves <- function(x, rows = 1L, baseline = TRUE) {
  sb_require_curve_plotting()
  if (!inherits(x, "stormbreaker_survival_curves")) {
    sb_abort(
      "Supply an sb_survival_curves result.", "stormbreaker_bad_survival_curves"
    )
  }
  if (!is.numeric(rows) || anyNA(rows) || any(rows != floor(rows)) ||
    any(rows < 1 | rows > ncol(x$survival)) || anyDuplicated(rows) ||
    !is.logical(baseline) || length(baseline) != 1L || is.na(baseline) ||
    (!length(rows) && !baseline)) {
    sb_abort(
      "Choose valid unique patient row numbers and/or the baseline.",
      "stormbreaker_bad_survival_curves"
    )
  }
  d <- data.frame(time = numeric(), survival = numeric(), series = character())
  if (baseline) {
    d <- data.frame(
      time = x$plot_baseline$time,
      survival = x$plot_baseline$survival, series = "Baseline (x beta = 0)"
    )
  }
  for (i in rows) {
    d <- rbind(d, data.frame(
      time = x$plot_baseline$time,
      survival = x$plot_survival[, i], series = paste("Observation", i)
    ))
  }
  mark <- d[d$time == x$prediction_horizon, , drop = FALSE]
  labels <- paste(mark$series, sprintf("%.1f%%", 100 * mark$survival), collapse = " | ")
  palette <- if (length(rows)) {
    stats::setNames(
      rep(c("#24547A", "#B86236", "#3A806E", "#805A8D"), length.out = length(rows)),
      paste("Observation", rows)
    )
  } else {
    character()
  }
  if (baseline) palette <- c("Baseline (x beta = 0)" = "#838B91", palette)
  time <- survival <- series <- NULL
  ggplot2::ggplot(d, ggplot2::aes(x = time, y = survival, colour = series)) +
    ggplot2::geom_step(direction = "hv", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = x$prediction_horizon, linetype = 3, colour = "#7A8791") +
    ggplot2::geom_point(data = mark, size = 2) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(
      title = "Model-predicted survival", colour = NULL,
      x = paste("Time |", x$time_unit), y = "Survival probability",
      subtitle = sb_viz_wrap(paste("Horizon:", x$prediction_horizon, x$time_unit, "|", labels)),
      caption = sb_viz_wrap(paste(
        "Fixed model coefficients and original-training baseline.",
        "The zero-score baseline is not an average patient. No uncertainty intervals.",
        "No extrapolation beyond observed training follow-up."
      ))
    ) +
    sb_viz_theme()
}

# ---- stable-survival ----

# Preserve a finite scale already defined by the fitted ordinary Cox model.
# No clipping, external-data refitting, or conversion of arbitrary risk scores.
sb_final_cox_log_cumulative_hazard <- function(fit, values) {
  if (is.null(values) || !identical(fit$learner_state$learner, "coxph") ||
      !any(values$survival_probability %in% c(0, 1))) return(NULL)
  h0 <- fit$learner_state$baseline_hazard$cumulative_baseline_hazard
  if (!is.numeric(h0) || length(h0) != 1L || !is.finite(h0) || h0 <= 0) return(NULL)
  stable <- log(h0) + values$risk
  if (!all(is.finite(stable))) return(NULL)
  # Emit only when this additional representation agrees with the saved
  # probability. Do not make previously valid predictions fail here.
  if (any(abs(exp(-exp(stable)) - values$survival_probability) >
      64 * .Machine$double.eps)) return(NULL)
  sb_validate_stable_survival_scale(stable, values$survival_probability)
}

sb_validate_stable_survival_scale <- function(value, probability) {
  if (!is.numeric(value) || length(value) != length(probability) ||
      anyNA(value) || any(!is.finite(value)) || !is.numeric(probability) ||
      anyNA(probability) || any(!is.finite(probability)) ||
      any(probability < 0 | probability > 1)) {
    sb_abort("The saved log cumulative hazard must be finite and patient-aligned.",
      "stormbreaker_bad_external_prediction")
  }
  reconstructed <- exp(-exp(value))
  if (any(abs(reconstructed - probability) > 64 * .Machine$double.eps)) {
    sb_abort("Saved log cumulative hazard and survival probability disagree.",
      "stormbreaker_bad_external_prediction")
  }
  unname(as.numeric(value))
}

sb_external_stable_scale <- function(rows, probability) {
  if (!"log_cumulative_hazard" %in% names(rows)) return(NULL)
  value <- sb_external_column(rows, "log_cumulative_hazard")
  sb_validate_stable_survival_scale(value, probability)
}

# ---- survival-probability ----

sb_validate_prediction_horizon <- function(prediction_horizon, analysis_time) {
  valid <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) &&
    is.finite(prediction_horizon) &&
    prediction_horizon > 0
  if (!valid) {
    sb_abort(
      "`prediction_horizon` must be one finite positive number.",
      "stormbreaker_bad_prediction_horizon"
    )
  }
  if (prediction_horizon > max(analysis_time)) {
    sb_abort(
      paste0(
        "The prediction horizon exceeds this outer analysis fold's ",
        "observed follow-up."
      ),
      "stormbreaker_prediction_horizon_unsupported"
    )
  }
  as.numeric(prediction_horizon)
}

sb_survival_probability_from_hazard <- function(risk,
                                                 cumulative_baseline_hazard) {
  valid_hazard <- is.numeric(cumulative_baseline_hazard) &&
    length(cumulative_baseline_hazard) == 1L &&
    !is.na(cumulative_baseline_hazard) &&
    is.finite(cumulative_baseline_hazard) &&
    cumulative_baseline_hazard >= 0
  if (!valid_hazard) {
    sb_abort(
      "The fitted baseline cumulative hazard is invalid.",
      "stormbreaker_invalid_baseline_hazard"
    )
  }
  if (length(risk) < 1L || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "The fitted linear predictors are incomplete or non-finite.",
      "stormbreaker_invalid_survival_prediction"
    )
  }

  probability <- if (cumulative_baseline_hazard == 0) {
    rep(1, length(risk))
  } else {
    log_patient_hazard <- log(cumulative_baseline_hazard) + risk
    patient_hazard <- exp(log_patient_hazard)
    exp(-patient_hazard)
  }
  if (length(probability) != length(risk) || anyNA(probability) ||
      any(!is.finite(probability)) ||
      any(probability < 0 | probability > 1)) {
    sb_abort(
      "The fitted survival probabilities are incomplete or invalid.",
      "stormbreaker_invalid_survival_prediction"
    )
  }
  unname(probability)
}

sb_breslow_baseline_hazard <- function(time, event, risk,
                                       prediction_horizon) {
  event_times <- sort(unique(time[event == 1L & time <= prediction_horizon]))
  if (length(event_times) == 0L) {
    return(0)
  }

  increments <- vapply(event_times, function(event_time) {
    risk_set <- risk[time >= event_time]
    maximum <- max(risk_set)
    log_denominator <- maximum + log(sum(exp(risk_set - maximum)))
    n_events <- sum(event == 1L & time == event_time)
    exp(log(n_events) - log_denominator)
  }, numeric(1))
  if (anyNA(increments) || any(!is.finite(increments)) ||
      any(increments <= 0)) {
    sb_abort(
      "The fixed-coefficient Breslow baseline hazard is invalid.",
      "stormbreaker_invalid_baseline_hazard"
    )
  }
  cumulative_hazard <- unname(sum(increments))
  if (!is.finite(cumulative_hazard) || cumulative_hazard < 0) {
    sb_abort(
      "The fixed-coefficient Breslow baseline hazard is invalid.",
      "stormbreaker_invalid_baseline_hazard"
    )
  }
  cumulative_hazard
}

sb_survival_prediction_result <- function(prediction_horizon,
                                           cumulative_baseline_hazard,
                                           source, ties, risk) {
  probability <- sb_survival_probability_from_hazard(
    risk,
    cumulative_baseline_hazard
  )
  list(
    prediction_horizon = as.numeric(prediction_horizon),
    cumulative_baseline_hazard = cumulative_baseline_hazard,
    baseline_survival = unname(exp(-cumulative_baseline_hazard)),
    source = source,
    ties = ties,
    probability = probability
  )
}

# ---- nomogram ----

# ---- nomogram-aft ----

sb_nomogram_aft_probability <- function(aft, score) {
  z <- (log(aft$prediction_horizon) + score) / aft$scale
  switch(aft$distribution,
    weibull = exp(-exp(z)), exponential = exp(-exp(z)),
    lognormal = stats::pnorm(z, lower.tail = FALSE),
    loglogistic = stats::plogis(z, lower.tail = FALSE))
}

sb_nomogram_aft_score <- function(aft, probability) {
  z <- switch(aft$distribution,
    weibull = log(-log(probability)), exponential = log(-log(probability)),
    lognormal = stats::qnorm(probability, lower.tail = FALSE),
    loglogistic = stats::qlogis(probability, lower.tail = FALSE))
  aft$scale * z - log(aft$prediction_horizon)
}

# ---- nomogram-overlays ----

sb_nomogram_overlays <- function(p, nomogram, reference_data, patient) {
  axes <- attr(p, "stormbreaker_nomogram_axes")
  features <- nomogram$features$feature
  shown <- features[features %in% axes$name]
  raised <- p$data$label_y > p$data$y
  p$data$label_y[raised] <- p$data$y[raised] - .45
  density_rows <- list()
  if (!is.null(reference_data)) {
    values <- sb_nomogram_values(nomogram, reference_data)
    if (any(values$outside_display_range)) stop("Reference observations extend beyond the nomogram display ranges.")
    add_density <- function(v, y) {
      if (length(v) < 2L || length(unique(v)) < 2L) return(NULL)
      d <- stats::density(v, from = min(v), to = max(v), n = 128)
      data.frame(x = d$x, y = y + .55 * d$y / max(d$y), axis = y)
    }
    for (f in shown) density_rows[[length(density_rows) + 1L]] <- add_density(values$points[, f], axes$y[match(f, axes$name)])
    if (nomogram$total_points > 0) density_rows[[length(density_rows) + 1L]] <- add_density(values$total_points / nomogram$total_points * 100, 2)
  }
  x <- y <- axis <- label <- NULL
  if (length(density_rows)) {
    d <- do.call(rbind, density_rows)
    if (!is.null(d)) p <- p + ggplot2::geom_line(data = d, ggplot2::aes(x, y, group = axis),
      inherit.aes = FALSE, colour = "#4477cc", linewidth = .45)
    attr(p, "stormbreaker_nomogram_density") <- d
  }
  if (!is.null(patient)) {
    if (!is.data.frame(patient) || nrow(patient) != 1L) stop("Supply exactly one patient row.")
    value <- sb_nomogram_values(nomogram, patient)
    if (value$outside_display_range) stop("Patient is outside the displayed ranges; recreate the nomogram with suitable ranges.")
    markers <- data.frame(x = as.numeric(value$points[1, shown]), y = axes$y[match(shown, axes$name)])
    p <- p + ggplot2::geom_segment(data = markers, ggplot2::aes(x = x, xend = x, y = y, yend = max(axes$y)),
      inherit.aes = FALSE, colour = "#d73027", linetype = 3, linewidth = .35)
    total_x <- if (nomogram$total_points > 0) value$total_points / nomogram$total_points * 100 else 0
    markers <- rbind(markers, data.frame(x = total_x, y = c(2, 1)))
    labels <- data.frame(x = total_x, y = c(2, 1) + .20,
      label = c(sprintf("%.1f", value$total_points), sprintf("%.3f", if (is.null(value$probability)) value$score else value$probability)))
    p <- p + ggplot2::geom_point(data = markers, ggplot2::aes(x, y), inherit.aes = FALSE, colour = "#d73027", size = 2) +
      ggplot2::geom_text(data = labels, ggplot2::aes(x, y, label = label), inherit.aes = FALSE, colour = "#d73027", size = 3)
    attr(p, "stormbreaker_nomogram_patient") <- value
  }
  p + ggplot2::labs(title = "Nomogram", subtitle = NULL, caption = NULL)
}

# ---- nomogram ----

#' @rdname sb_survival_curves
#' @export
sb_nomogram <- function(fit, reference_data) {
  parameters <- sb_model_coefficients(fit)
  aft <- NULL
  if (identical(parameters$score_scale, "negative_log_time_location") &&
      !is.null(fit$learner_state$survival_probability)) {
    aft <- sb_validate_additional_survival_probability_state(fit$learner_state, fit$learner_state$learner)
    aft$distribution <- fit$learner_state$params$distribution
    aft$scale <- fit$learner_state$scale
  }
  tab <- parameters$table
  features <- tab$feature
  values <- sb_select_final_feature_frame(reference_data, features,
    "Reference data", "nomogram feature ranges")
  if (!nrow(values)) sb_abort("Reference data must contain observations.",
    "stormbreaker_bad_reference_data")
  tab$minimum <- vapply(values, min, numeric(1))
  tab$maximum <- vapply(values, max, numeric(1))
  tab$minimum_contribution <- pmin(tab$minimum * tab$coefficient, tab$maximum * tab$coefficient)
  tab$span <- abs(tab$coefficient) * (tab$maximum - tab$minimum)
  largest <- max(tab$span)
  unit <- if (largest > 0) largest / 100 else 1
  intercept <- if (is.null(parameters$intercept)) 0 else parameters$intercept
  tab$maximum_points <- tab$span / unit
  structure(list(features = tab, intercept = intercept,
    score_origin = intercept + sum(tab$minimum_contribution), score_per_point = unit,
    total_points = sum(tab$maximum_points), score_scale = parameters$score_scale,
    binary = identical(parameters$score_scale, "class_1_log_odds"),
    baseline = fit$learner_state$baseline_hazard, aft = aft, pipeline_id = parameters$pipeline_id,
    reference_n = nrow(values)), class = "stormbreaker_nomogram")
}

#' @rdname sb_survival_curves
#' @export
sb_nomogram_values <- function(x, new_data) {
  if (!inherits(x, "stormbreaker_nomogram")) sb_abort(
    "Supply an sb_nomogram result.", "stormbreaker_bad_visualization_argument")
  tab <- x$features
  values <- as.matrix(sb_select_final_feature_frame(new_data, tab$feature,
    "New data", "nomogram predictors"))
  contribution <- sweep(values, 2L, tab$coefficient, "*")
  points <- sweep(contribution, 2L, tab$minimum_contribution, "-") / x$score_per_point
  total <- unname(rowSums(points))
  score <- x$score_origin + total * x$score_per_point
  outside <- unname(rowSums(sweep(values, 2L, tab$minimum, "<") |
    sweep(values, 2L, tab$maximum, ">")) > 0)
  probability <- NULL
  if (x$binary) probability <- stats::plogis(score) else if (!is.null(x$aft)) {
    probability <- sb_nomogram_aft_probability(x$aft, score)
  } else if (!is.null(x$baseline)) {
    hazard <- x$baseline$cumulative_baseline_hazard
    if (!is.numeric(hazard) || length(hazard) != 1L || !is.finite(hazard) || hazard < 0) {
      sb_abort("Saved survival baseline is invalid.", "stormbreaker_invalid_baseline_hazard")
    }
    probability <- if (hazard == 0) rep(1, length(score)) else exp(-exp(log(hazard) + score))
  }
  list(points = points, total_points = total, score = score,
    probability = probability, probability_type = if (x$binary) "class_1" else
      if (!is.null(probability)) "survival_at_saved_horizon" else NULL,
    outside_display_range = outside)
}

#' @rdname sb_survival_curves
#' @export
sb_plot_nomogram <- function(x, ticks = 5L, show_zero = FALSE,
                             reference_data = NULL, patient = NULL) {
  sb_require_curve_plotting()
  if (!inherits(x, "stormbreaker_nomogram")) sb_abort(
    "Supply an sb_nomogram result.", "stormbreaker_bad_visualization_argument")
  ticks <- sb_viz_positive_integer(ticks, "ticks")
  if (!is.logical(show_zero) || length(show_zero) != 1L || is.na(show_zero)) sb_abort(
    "show_zero must be TRUE or FALSE.", "stormbreaker_bad_visualization_argument")
  tab <- x$features
  omitted <- if (show_zero) character() else tab$feature[tab$coefficient == 0]
  if (!show_zero) tab <- tab[tab$coefficient != 0, , drop = FALSE]
  k <- nrow(tab)
  axis <- data.frame(name = c("Points", tab$feature, "Total points", "Model score"),
    y = rev(seq_len(k + 3)), start = 0,
    end = c(100, tab$maximum_points, 100, 100))
  if (identical(x$score_scale, "negative_log_time_location")) axis$name[k + 3] <- "Negative log-time location"
  marks <- list(data.frame(x = seq(0, 100, 20), y = k + 3,
    label = as.character(seq(0, 100, 20)), offset = -.25))
  for (j in seq_len(k)) {
    count <- max(2L, min(ticks, floor(tab$maximum_points[j] / 15) + 1L))
    v <- unique(seq(tab$minimum[j], tab$maximum[j], length.out = count))
    pos <- (v * tab$coefficient[j] - tab$minimum_contribution[j]) / x$score_per_point
    if (tab$coefficient[j] == 0) {
      v <- paste(format(signif(tab$minimum[j], 3)), "to", format(signif(tab$maximum[j], 3)), "(0 points)")
      pos <- 0
    } else v <- format(signif(v, 3), trim = TRUE)
    offset <- if (length(pos) == 2L && abs(diff(pos)) < 12) c(.25, -.3) else -.25
    marks[[length(marks) + 1L]] <- data.frame(x = pos, y = k + 3 - j,
      label = v, offset = offset)
  }
  total <- unique(seq(0, x$total_points, length.out = max(2L, ticks)))
  pos <- if (x$total_points > 0) total / x$total_points * 100 else rep(0, length(total))
  marks[[length(marks) + 1L]] <- data.frame(x = pos, y = 2,
    label = format(signif(total, 3), trim = TRUE), offset = -.25)
  score <- x$score_origin + total * x$score_per_point
  label <- format(signif(score, 3), trim = TRUE)
  if (x$binary || !is.null(x$baseline) || !is.null(x$aft)) {
    probability <- c(.001, .01, .05, .1, .2, .4, .5, .6, .8, .9, .95, .99, .999)
    desired_score <- if (x$binary) stats::qlogis(probability) else if (!is.null(x$aft))
      sb_nomogram_aft_score(x$aft, probability) else
      log(-log(probability)) - log(x$baseline$cumulative_baseline_hazard)
    pts <- (desired_score - x$score_origin) / x$score_per_point
    keep <- is.finite(pts) & pts >= 0 & pts <= x$total_points
    if (any(keep) && x$total_points > 0) {
      pos <- pts[keep] / x$total_points * 100
      label <- format(probability[keep], trim = TRUE)
      # Keep readable tick spacing, prioritizing the middle probability.
      priority <- order(abs(probability[keep] - .5))
      chosen <- integer()
      for (i in priority) {
        if (!length(chosen) || all(abs(pos[i] - pos[chosen]) >= 8)) chosen <- c(chosen, i)
      }
      chosen <- chosen[order(pos[chosen])]
      pos <- pos[chosen]
      label <- label[chosen]
    } else {
      pos <- 0
      p <- if (x$binary) stats::plogis(x$score_origin) else if (!is.null(x$aft))
        sb_nomogram_aft_probability(x$aft, x$score_origin) else
        exp(-x$baseline$cumulative_baseline_hazard * exp(x$score_origin))
      label <- format(signif(p, 4), trim = TRUE)
    }
    axis$name[k + 3] <- if (x$binary) "Probability (class 1)" else if (!is.null(x$aft))
      paste("Survival at", x$aft$prediction_horizon, x$aft$time_unit) else
      paste("Survival at", x$baseline$prediction_horizon, x$baseline$time_unit)
  }
  marks[[length(marks) + 1L]] <- data.frame(x = pos, y = 1, label = label, offset = -.25)
  data <- do.call(rbind, marks)
  data$label_y <- data$y + data$offset
  start <- end <- y <- label <- label_y <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(x, y)) +
    ggplot2::geom_segment(data = axis, ggplot2::aes(x = start, xend = end, y = y, yend = y), inherit.aes = FALSE) +
    ggplot2::geom_segment(ggplot2::aes(xend = x, yend = y - .10)) +
    ggplot2::geom_text(ggplot2::aes(y = label_y, label = label), size = 3) +
    ggplot2::scale_y_continuous(breaks = axis$y, labels = axis$name, limits = c(.4, k + 3.5)) +
    ggplot2::scale_x_continuous(limits = c(-3, 103)) +
    ggplot2::labs(title = "Nomogram", x = NULL, y = NULL) +
    sb_viz_theme() + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank())
  attr(p, "stormbreaker_nomogram_axes") <- axis
  attr(p, "stormbreaker_nomogram_omitted_zero_terms") <- omitted
  attr(p, "stormbreaker_notes") <- c(
    "Sum feature points, then read total-points and outcome axes. Short axes have staggered endpoint labels.",
    "Original score and probability mapping; no refitting. Display ranges come from reference data.",
    if (length(omitted)) paste(length(omitted),
      "exactly zero-coefficient terms omitted from axes; retained in model/table.")
  )
  if (!is.null(reference_data) || !is.null(patient)) p <- sb_nomogram_overlays(p, x, reference_data, patient)
  p
}
