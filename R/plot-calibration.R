#' Calibration and clinical plots
#' @name sb_plot_calibration
#' @description Plot calibration, decision curves, thresholds and saved clinical summaries.
#' @rawRd \arguments{
#' \item{x}{
#' \strong{For \code{sb_plot_calibration()}:} A result from \code{\link[=sb_run_cv]{sb_run_cv()}}, \code{\link[=sb_compare_cv]{sb_compare_cv()}},
#' \code{\link[=sb_evaluate_external]{sb_evaluate_external()}}, or \code{\link[=sb_compare_external]{sb_compare_external()}}. Each comparison is
#' displayed in supplied order with one panel per usable pipeline. Whole
#' repeated CV, IECV and stacking CV are supported with separate context
#' panels. IECV remains development validation, not final external validation.
#' 
#' \strong{For \code{sb_plot_dca()}:} A result from \code{\link[=sb_run_cv]{sb_run_cv()}}, \code{\link[=sb_compare_cv]{sb_compare_cv()}},
#' \code{\link[=sb_evaluate_external]{sb_evaluate_external()}}, or \code{\link[=sb_compare_external]{sb_compare_external()}}. Comparisons keep one
#' panel per usable pipeline in supplied order, with references from that same
#' pipeline's saved curve. Whole repeated CV, IECV and stacking CV retain
#' separate context panels and reference curves; they are never pooled.
#' 
#' \strong{For \code{sb_plot_thresholds()}:} A result from \code{\link[=sb_threshold_metrics]{sb_threshold_metrics()}}.
#' 
#' \strong{For \code{sb_plot_clinician_panel()}:} A result from \code{\link[=sb_clinician_views]{sb_clinician_views()}}.
#' 
#' \strong{For \code{sb_plot_probability_distribution()}:} An \code{\link[=sb_probability_distribution]{sb_probability_distribution()}} result.
#' 
#' \strong{For \code{sb_clinician_views()}:} Saved CV or external evaluation results, or a supported collection.
#' Results from \code{\link[=bi_run]{bi_run()}} and \code{\link[=surv_run]{surv_run()}} display development CV;
#' \code{\link[=sb_complete_run]{sb_complete_run()}} results display saved external evaluations by cohort.
#' }
#' \item{cohort_label}{
#' \strong{For \code{sb_plot_calibration()}:} Optional non-empty display label for the cohort. The
#' development/external evaluation boundary is retained in returned metadata.
#' Default figures omit process/status text and explanatory footnotes. A label does
#' not establish independence or external validity.
#' 
#' \strong{For \code{sb_plot_dca()}:} Optional non-empty cohort display label. Development
#' versus external evaluation remains explicit; a name does not prove
#' cohort independence.
#' }
#' \item{show_unavailable}{
#' \strong{For \code{sb_plot_calibration()}:} Logical; default FALSE hides unavailable panels.
#' TRUE retains labeled audit panels for failed or unrequested curves.
#' 
#' \strong{For \code{sb_plot_dca()}:} Logical; default FALSE hides unavailable panels.
#' }
#' \item{standardized}{
#' Divide binary net benefit by the cohort event prevalence.
#' Not implemented for survival outcomes. No model is refitted.
#' }
#' \item{cost_benefit_axis}{
#' Add a secondary threshold-odds cost:benefit scale.
#' TRUE retains the previous labeled audit panels.
#' }
#' \item{metrics}{
#' Metrics to display. Defaults to sensitivity, specificity,
#' positive predictive value and negative predictive value.
#' }
#' \item{ncol}{
#' \strong{For \code{sb_plot_thresholds()}, \code{sb_plot_probability_distribution()}:} Number of panel columns.
#' 
#' \strong{For \code{sb_plot_clinician_panel()}:} Number of columns in the combined figure.
#' }
#' \item{views}{
#' Explicit available plot names in the desired panel order.
#' }
#' \item{scale}{
#' Display count or fraction within each observed class.
#' }
#' \item{thresholds}{
#' Optional explicit binary probability thresholds. No optimum
#' is selected; saved decision curves keep their own original thresholds.
#' }
#' \item{confusion_threshold}{
#' Optional single threshold for the confusion matrix.
#' It must be included in thresholds.
#' }
#' \item{native}{
#' Optional \code{\link[=sb_native_diagnostics]{sb_native_diagnostics()}} result or uniquely named list
#' of extracted algorithm diagnostics.
#' These are already fitted-model summaries; no extraction refit is triggered.
#' }
#' \item{shap}{
#' Optional saved \code{\link[=sb_additive_shap]{sb_additive_shap()}} result. Adds global importance
#' and beeswarm views; never recomputes contributions or creates a surrogate.
#' }
#' \item{feature_values}{
#' Optional feature frame and explicit ID column for
#' SHAP colouring and dependence. IDs must match explained observations exactly.
#' }
#' \item{id_col}{
#' Optional feature frame and explicit ID column for
#' SHAP colouring and dependence. IDs must match explained observations exactly.
#' }
#' \item{dependence_features}{
#' Explicit explained feature names for dependence
#' plots. NULL adds none; no automatic mass per-feature export.
#' }
#' \item{nomogram}{
#' Optional exact \code{\link[=sb_nomogram]{sb_nomogram()}} result. No new model is fitted.
#' }
#' \item{selectors}{
#' Optional uniquely named list of saved
#' \code{\link[=sb_selector_diagnostics]{sb_selector_diagnostics()}} results, shown separately from model evidence.
#' }
#' \item{detail}{
#' Default \code{"brief"} requests ROC, calibration and saved decision
#' curves only. \code{"full"} additionally requests precision-recall and probability
#' distributions. Explicit diagnostics and explanations remain opt-in in both.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_plot_calibration()}}{
#' Display existing descriptive calibration curves without fitting a new
#' calibration equation or changing model predictions.
#' 
#' \strong{Returns:} A \code{ggplot} object, requiring optional \code{ggplot2}. \code{plot$data} contains
#' copied \code{x}/\code{y} coordinates, source row numbers, pipeline/cohort/scope and
#' status information. The \code{stormbreaker_curve_status} attribute is a data
#' frame with one diagnostic row per pipeline. No file is written.
#' \code{stormbreaker_curve_data} retains all source coordinates, including hidden
#' pipelines; \code{plot$data} contains only displayed pipelines by default.
#' 
#' \strong{Details:} Binary coordinates come from the saved moderate-calibration curve.
#' Survival coordinates are the saved predicted and calibrated \emph{event}
#' probabilities at the declared horizon, not survival probabilities. The
#' axis target and panel horizon/unit identify this distinction. The
#' calibrated coordinate is the existing smoothed descriptive estimate of
#' observed probability, not a new model prediction to apply to patients.
#' 
#' Both metric and complete-run status must be successful for a curve to be
#' drawn, and at least one saved line segment must be drawable. Otherwise its
#' panel is hidden by default; source coordinates remain in diagnostic
#' attributes and the original result. If none are usable, a no-results
#' message is shown. Use \code{\link[=sb_export_failures]{sb_export_failures()}} for saved failure details.
#' Missing/non-finite coordinates are not filled, and lines do not bridge
#' excluded coordinates. An absent survival horizon is shown as not requested
#' for an otherwise successful run. Missing information after a failure is
#' not relabeled as an unrequested analysis.
#' 
#' The dashed identity line is an ideal-calibration reference. No confidence
#' band, refitting, smoothing, ranking or threshold selection is performed.
#' Curves remain development or external evaluation as recorded by the source.
#' Comparisons can require a taller output image for readable panel labels.
#' }
#' 
#' \subsection{\code{sb_plot_dca()}}{
#' Display saved model, treat-all and treat-none net benefits at the original
#' study-declared threshold probabilities, without recalculating utility.
#' 
#' \strong{Returns:} A \code{ggplot} object using optional \code{ggplot2}. \code{plot$data} preserves
#' original threshold (\code{x}) and net-benefit (\code{y}) coordinates for each
#' \code{strategy}, with source row numbers and provenance. The attribute
#' \code{stormbreaker_curve_status} contains one diagnostic row per pipeline.
#' \code{stormbreaker_curve_data} retains all original coordinates, including
#' hidden pipelines; \code{plot$data} contains displayed pipelines by default.
#' No file is written.
#' 
#' \strong{Details:} DCA must already have been requested using \code{decision_thresholds}
#' when declaring the development task. This plotting function accepts no
#' new thresholds and cannot compute missing curves. Survival thresholds
#' apply to event probability by the saved horizon; its time unit and horizon
#' are displayed. The saved censoring-aware estimator is not replaced.
#' 
#' Both curve and complete-run status must be successful for points to be
#' drawn. A pipeline also needs at least one usable model net-benefit point.
#' Unrequested and failed curves are hidden by default but retain diagnostic
#' attributes. If no pipelines are usable, an explicit no-results message
#' replaces the plot. Use \code{\link[=sb_export_failures]{sb_export_failures()}} for failure details.
#' Non-finite or unsuccessful row-level coordinates are not
#' drawn; separate line segments prevent bridging their gaps. Original values
#' remain in \code{plot$data}, including negative net benefits. Empty panels use
#' explicit missing-coordinate placeholders rather than invented values.
#' 
#' Points identify the actual saved thresholds; connecting segments are
#' visual guides, not newly evaluated thresholds or a smoothed curve. No
#' confidence interval, optimal threshold, ranking or clinical recommendation
#' is inferred. Reference lines remain within their original pipeline panels.
#' Use a taller output image for comparisons with many pipelines.
#' }
#' 
#' \subsection{\code{sb_plot_thresholds()}}{
#' Plot explicit threshold performance
#' 
#' \strong{Returns:} A ggplot object with group availability attached. Undefined values
#' remain missing. No model fitting or threshold optimization occurs.
#' }
#' 
#' \subsection{\code{sb_plot_clinician_panel()}}{
#' Combine selected clinical views into one figure
#' 
#' \strong{Returns:} A patchwork plot, suitable for ggplot2::ggsave(). The
#' \code{stormbreaker_captions} attribute is a data frame of view names and separate
#' figure legends in panel order. Requires the optional patchwork package.
#' 
#' \strong{Details:} Only layout changes: original plot data, layers, axes, legends and
#' evaluation contexts are retained. Titles use the short selected view names;
#' original titles and explanations stay in the separate caption table.
#' Choose related views and export at a readable physical size. For many
#' cohorts or complex faceted diagnostics, use fewer panels or separate figures
#' rather than shrinking labels. No files are written and no models are fitted.
#' }
#' 
#' \subsection{\code{sb_plot_probability_distribution()}}{
#' Plot saved probability distributions without patient identifiers
#' 
#' \strong{Returns:} A ggplot with adjacent class histograms for each fold or cohort.
#' Missing classes have zero counts and undefined fractions, not fabricated
#' density estimates. All group diagnostics are attached to the plot.
#' }
#' 
#' \subsection{\code{sb_clinician_views()}}{
#' Assemble clinical views from saved evaluation results
#' 
#' \strong{Returns:} A list with named plots, topic groups, separate captions, an availability
#' table and original diagnostic tables. Unavailable views are recorded, not
#' replaced with empty figures or fabricated metrics. Partial-context plots
#' retain their original fold/cohort diagnostics. No model is fitted or selected.
#' 
#' \strong{Details:} ROC, precision-recall and probability distributions require binary
#' predictions; survival inputs retain their saved calibration and decision
#' curves. Unsupported binary views are recorded as unavailable. Missing DCA
#' is not recomputed. SHAP and nomograms must be supplied from their exact-model
#' extractors; no patient card or subgroup module is added. Individual SHAP
#' values remain sensitive even though original identifiers are not plotted.
#' Technical diagnostics are opt-in: call an extractor separately and supply
#' its result through native or selectors. Report assembly never discovers
#' diagnostics from a model, including when all report topics are expanded.
#' Figures use readable display labels; original IDs and status remain in
#' their data and diagnostics. Model settings are kept separately in models.
#' }
#' }
#' @rawRd \seealso{
#' \code{\link[=sb_plot_dca]{sb_plot_dca()}}, \code{\link[=sb_plot_performance]{sb_plot_performance()}}
#' 
#' \code{\link[=sb_plot_calibration]{sb_plot_calibration()}}, \code{\link[=sb_plot_performance]{sb_plot_performance()}}
#' }
#' @rawRd \examples{
#' # sb_plot_calibration
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   d <- sb_example_binary(120, seed = 1)
#'   task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
#'   result <- sb_run_cv(
#'     task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1
#'   )
#'   sb_plot_calibration(result)
#' }
#' 
#' # sb_plot_dca
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   d <- sb_example_binary(120, seed = 1)
#'   task <- sb_binary_task(
#'     d, "patient_id", "outcome", "gene_a",
#'     decision_thresholds = c(0.1, 0.2, 0.3, 0.5)
#'   )
#'   result <- sb_run_cv(
#'     task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1
#'   )
#'   sb_plot_dca(result)
#' }
#' 
#' # sb_clinician_views
#' \dontshow{if (requireNamespace("ggplot2", quietly = TRUE)) withAutoprint(\{ # examplesIf}
#' development <- bi_example(40, seed = 42)
#' run <- bi_run(development, outcome = "outcome", id = "patient_id",
#'   predictors = c("gene_a", "gene_b"), models = "Logistic", folds = 2,
#'   seed = 42, plots = FALSE)
#' views <- clinical_plots(run)
#' names(views$plots)
#' views$plots$ROC
#' \dontshow{\}) # examplesIf}
#' }
NULL

# ---- plot-calibration ----

# ---- clinical-thresholds ----

#' @rdname sb_evaluate_external
#' @export
sb_threshold_metrics <- function(x, thresholds, cohort_label = NULL) {
  if (!is.numeric(thresholds) || !length(thresholds) ||
      any(!is.finite(thresholds)) || any(thresholds < 0 | thresholds > 1)) {
    sb_abort("Supply explicit finite thresholds in [0, 1].",
      "stormbreaker_bad_visualization_argument")
  }
  thresholds <- sort(unique(thresholds))
  leaves <- sb_viz_leaves(x, cohort_label)
  tables <- diagnostics <- saved <- list()
  ratio <- function(a, b) if (b == 0) NA_real_ else a / b
  for (leaf in leaves) {
    value <- leaf$value
    external <- inherits(value, "stormbreaker_external_evaluation")
    binary <- if (external) identical(value$task_type, "binary") else
      inherits(value, "stormbreaker_binary_cv_result")
    if (!binary) sb_abort("Decision thresholds require binary probabilities.",
      "stormbreaker_bad_visualization_task")
    pred <- if (external) value$scored_predictions else value$predictions
    saved[[leaf$context$view_id]] <- value$metrics
    folds <- if (!external && is.data.frame(pred) && nrow(pred) &&
      "fold" %in% names(pred)) unique(pred$fold) else NA_integer_
    for (fold in folds) {
      group <- if (!is.data.frame(pred)) data.frame() else if (external || is.na(fold)) pred else
        pred[!is.na(pred$fold) & pred$fold == fold, , drop = FALSE]
      status <- if (!identical(value$status, "success")) "failed_run" else if (!nrow(group)) "no_predictions" else "success"
      if (status == "success") {
        valid <- all(c("outcome", "probability", "status", "patient_id") %in% names(group)) &&
          is.numeric(group$outcome) && is.numeric(group$probability) &&
          all(is.finite(group$outcome)) && all(group$outcome %in% c(0, 1)) &&
          all(is.finite(group$probability)) && all(group$probability >= 0 & group$probability <= 1) &&
          !anyNA(group$status) && all(group$status == "success") &&
          !anyNA(group$patient_id) && !anyDuplicated(group$patient_id) &&
          (external || ("fold" %in% names(group) && !anyNA(group$fold)))
        if (!valid) status <- "invalid_predictions"
      }
      d <- data.frame(pipeline_id = sb_viz_pipeline_id(value), dataset = leaf$dataset,
        fold = as.integer(fold), panel = paste(sb_viz_panel(leaf),
          if (external) "External patients" else paste("Assessment fold", fold), sep = "\n"),
        status = status, n = nrow(group), stringsAsFactors = FALSE)
      d <- sb_view_annotate(d, leaf$context)
      d$panel <- paste(c(d$pipeline_id,
        if (leaf$dataset != "Unnamed cohort") leaf$dataset,
        if (!is.na(leaf$context$repeat_id)) paste("Repeat", leaf$context$repeat_id),
        if (!is.na(leaf$context$heldout_cohort_id)) leaf$context$heldout_cohort_id,
        if (identical(leaf$context$validation_scheme, "internal_holdout")) "Internal holdout" else
          if (external) "External" else paste("CV fold", fold)), collapse = "\n")
      diagnostics[[length(diagnostics) + 1L]] <- d
      if (status != "success") next
      rows <- lapply(thresholds, function(threshold) {
        positive <- group$probability >= threshold
        tp <- sum(positive & group$outcome == 1)
        fp <- sum(positive & group$outcome == 0)
        tn <- sum(!positive & group$outcome == 0)
        fn <- sum(!positive & group$outcome == 1)
        n <- nrow(group)
        observed <- (tp + tn) / n
        expected <- ((tp + fp) / n) * ((tp + fn) / n) +
          ((tn + fn) / n) * ((tn + fp) / n)
        data.frame(threshold = threshold, tp = tp, fp = fp, tn = tn, fn = fn,
          sensitivity = ratio(tp, tp + fn), specificity = ratio(tn, tn + fp),
          ppv = ratio(tp, tp + fp), npv = ratio(tn, tn + fn),
          accuracy = ratio(tp + tn, nrow(group)), f1 = ratio(2 * tp, 2 * tp + fp + fn),
          kappa = ratio(observed - expected, 1 - expected),
          positive_fraction = ratio(tp + fp, nrow(group)))
      })
      tab <- do.call(rbind, rows)
      for (name in names(d)) tab[[name]] <- rep(d[[name]], nrow(tab))
      tables[[length(tables) + 1L]] <- tab
    }
  }
  structure(list(table = if (length(tables)) sb_view_bind(tables) else data.frame(),
    diagnostics = if (length(diagnostics)) sb_view_bind(diagnostics) else data.frame(),
    saved_metrics = saved, threshold_rule = "probability >= threshold"),
    class = "stormbreaker_threshold_metrics")
}

#' @rdname sb_plot_calibration
#' @export
sb_plot_thresholds <- function(x,
                               metrics = c("sensitivity", "specificity", "ppv", "npv"),
                               ncol = 2L) {
  sb_require_curve_plotting()
  if (!inherits(x, "stormbreaker_threshold_metrics")) sb_abort(
    "Supply an sb_threshold_metrics result.", "stormbreaker_bad_visualization_argument")
  allowed <- c("sensitivity", "specificity", "ppv", "npv", "accuracy", "f1", "kappa", "positive_fraction")
  if (!is.character(metrics) || !length(metrics) || anyNA(metrics) ||
      anyDuplicated(metrics) || any(!metrics %in% allowed)) sb_abort(
    "Supply distinct supported threshold metrics.", "stormbreaker_bad_visualization_argument")
  ncol <- sb_viz_positive_integer(ncol, "ncol")
  if (!nrow(x$table)) {
    p <- sb_empty_display_plot("Threshold performance", "No usable binary predictions.",
      "Unavailable groups remain in diagnostics.")
  } else {
    data <- do.call(rbind, lapply(metrics, function(m) {
      d <- x$table
      d$metric <- m
      d$value <- d[[m]]
      d
    }))
    threshold <- value <- metric <- panel <- NULL
    p <- ggplot2::ggplot(data, ggplot2::aes(threshold, value, colour = metric)) +
      ggplot2::geom_line(na.rm = TRUE) + ggplot2::geom_point(size = 1, na.rm = TRUE) +
      ggplot2::facet_wrap(~panel, ncol = ncol, labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
      ggplot2::scale_y_continuous(limits = c(if ("kappa" %in% metrics) -1 else 0, 1)) +
      ggplot2::labs(title = "How does the decision threshold change performance?",
        x = "Explicit probability threshold", y = "Metric value", colour = "Metric",
        caption = "Class 1 is positive; probability >= threshold. Separate folds/cohorts. No threshold selected. Undefined ratios are missing.") +
      sb_viz_theme()
  }
  attr(p, "stormbreaker_threshold_diagnostics") <- x$diagnostics
  p
}

# ---- plot-calibration ----

#' @rdname sb_plot_calibration
#' @export
sb_plot_calibration <- function(x, cohort_label = NULL, show_unavailable = FALSE) {
  sb_require_curve_plotting()
  original <- sb_curve_view_data(x, "calibration", cohort_label)
  view <- sb_curve_display_view(original, "calibration", show_unavailable)
  if (!nrow(view$status)) {
    p <- sb_empty_display_plot("Saved calibration curve",
      "No successful saved calibration curves are available.",
      "Use sb_export_failures() for details; unrequested curves are not failures.")
    p$data <- view$data
    attr(p, "stormbreaker_curve_status") <- original$status
    attr(p, "stormbreaker_curve_data") <- original$data
    return(p)
  }
  lines <- sb_curve_view_lines(view$data)
  annotations <- sb_curve_view_annotations(view)
  survival <- all(view$status$task_type == "survival")
  target <- if (survival) "event probability" else "probability"
  panel_label <- function(labels) {
    if (nrow(view$status) == 1L) return(rep(if (survival)
      paste(view$status$prediction_horizon, view$status$time_unit) else "", length(labels)))
    sub("\ncurve:.*", "", labels)
  }
  x <- y <- panel <- label <- line_group <- NULL
  p <- ggplot2::ggplot(view$data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey55") +
    ggplot2::geom_path(data = lines, ggplot2::aes(group = line_group),
                       colour = "#24547A") +
    ggplot2::geom_text(data = annotations,
                       ggplot2::aes(x = 0.5, y = 0.5, label = label),
                       inherit.aes = FALSE, colour = "grey35", size = 3.5) +
    ggplot2::facet_wrap(~panel, ncol = 1, drop = FALSE,
      labeller = ggplot2::as_labeller(panel_label)) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = "Calibration",
      x = paste("Predicted", target), y = paste("Observed", target)
    ) +
    sb_viz_theme()
  attr(p, "stormbreaker_curve_status") <- original$status
  attr(p, "stormbreaker_curve_data") <- original$data
  p
}

# ---- plot-dca ----

#' @rdname sb_plot_calibration
#' @export
sb_plot_dca <- function(x, cohort_label = NULL, show_unavailable = FALSE,
                       standardized = FALSE, cost_benefit_axis = FALSE) {
  sb_require_curve_plotting()
  for (flag in list(standardized, cost_benefit_axis)) if (!is.logical(flag) ||
      length(flag) != 1L || is.na(flag)) stop("DCA display options must be TRUE or FALSE.")
  original <- sb_curve_view_data(x, "dca", cohort_label)
  models <- sb_clinician_models(sb_viz_leaves(x, cohort_label))
  if (standardized) {
    leaves <- sb_viz_leaves(x, cohort_label)
    if (!all(vapply(leaves, function(z) identical(z$value$task_type, "binary") ||
      inherits(z$value, "stormbreaker_binary_cv_result"), logical(1)))) stop("Standardized DCA currently supports binary outcomes only.")
    d <- original$data
    keys <- intersect(c("view_id", "pipeline_id", "panel"), names(d))
    key <- interaction(d[keys], drop = TRUE)
    for (k in levels(key)) {
      rows <- which(key == k)
      ref <- rows[d$strategy[rows] == "treat_all" & is.finite(d$y[rows]) & is.finite(d$x[rows])]
      prevalence <- d$y[ref] * (1 - d$x[ref]) + d$x[ref]
      if (!length(prevalence)) {
        if (any(is.finite(d$y[rows]))) stop("Standardization requires the saved treat-all reference.")
        next
      }
      if (any(prevalence <= 0) || diff(range(prevalence)) > 1e-8) stop("A positive common prevalence is required for standardized net benefit.")
      d$y[rows] <- d$y[rows] / mean(prevalence)
    }
    original$data <- d
  }
  view <- sb_curve_display_view(original, "dca", show_unavailable)
  if (!nrow(view$status)) {
    p <- sb_empty_display_plot("Saved decision curve",
      "No successful saved decision curves are available.",
      "Use sb_export_failures() for details; unrequested curves are not failures.")
    p$data <- view$data
    attr(p, "stormbreaker_curve_status") <- original$status
    attr(p, "stormbreaker_curve_data") <- original$data
    return(p)
  }
  lines <- sb_curve_view_lines(view$data)
  points <- view$data[view$data$draw, , drop = FALSE]
  annotations <- sb_curve_view_annotations(view)
  thresholds <- view$data$x[is.finite(view$data$x)]
  threshold_range <- if (length(unique(thresholds)) < 2L) c(0, 1) else range(thresholds)
  annotations$annotation_x <- rep(mean(threshold_range), nrow(annotations))
  x <- y <- strategy <- line_group <- panel <- label <- annotation_x <- NULL
  p <- ggplot2::ggplot(view$data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_path(data = lines,
                       ggplot2::aes(colour = strategy, linetype = strategy,
                                    group = line_group)) +
    ggplot2::geom_point(data = points, ggplot2::aes(colour = strategy), size = 1.5) +
    ggplot2::geom_text(data = annotations,
                       ggplot2::aes(x = annotation_x, y = 0, label = label),
                       inherit.aes = FALSE, colour = "grey35", size = 3.5) +
    ggplot2::facet_wrap(~panel, ncol = 1, drop = FALSE,
      labeller = if ("view_id" %in% names(view$status))
        ggplot2::as_labeller(sb_view_wrap_label) else "label_value") +
    ggplot2::scale_x_continuous(limits = threshold_range,
      sec.axis = if (cost_benefit_axis) ggplot2::dup_axis(name = "Cost:Benefit",
        breaks = c(.01, .1, .2, .5, .8), labels = c("1:99", "1:9", "1:4", "1:1", "4:1")) else ggplot2::waiver()) +
    ggplot2::scale_colour_manual(
      name = NULL, values = c(model = "#24547A", treat_all = "grey50", treat_none = "grey15"),
      limits = c("model", "treat_all", "treat_none"),
      breaks = c("model", "treat_all", "treat_none"),
      labels = c("Model", "Treat all", "Treat none"), drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      name = NULL, values = c(model = 1, treat_all = 2, treat_none = 3),
      limits = c("model", "treat_all", "treat_none"),
      breaks = c("model", "treat_all", "treat_none"),
      labels = c("Model", "Treat all", "Treat none"), drop = FALSE
    ) +
    ggplot2::labs(
      title = "Saved decision curve",
      subtitle = sb_curve_view_scope_label(view),
      x = "Threshold probability", y = "Net benefit",
      caption = paste0(
        "Points: saved thresholds; lines: connecting guides. No threshold selection or confidence intervals.\n",
        sum(!view$data$draw & !is.na(view$data$source_row)),
        " saved coordinate(s) not drawn; inspect attr(plot, 'stormbreaker_curve_status')."
      )
    ) +
    sb_viz_theme() +
    ggplot2::theme(legend.position = "bottom")
  attr(p, "stormbreaker_curve_status") <- original$status
  attr(p, "stormbreaker_curve_data") <- original$data
  if (standardized) p <- p + ggplot2::labs(y = "Standardized net benefit")
  attr(p, "stormbreaker_figure_notes") <- p$labels$caption
  p <- sb_clinician_clean_plot(p, models)
  p <- p + ggplot2::labs(title = "Decision curve", subtitle = NULL, caption = NULL)
  p
}

# ---- plot-clinical ----

# ---- curve-view-data ----

sb_require_curve_plotting <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    sb_abort(
      "Plotting requires ggplot2; install it with install.packages('ggplot2').",
      "stormbreaker_missing_plot_dependency"
    )
  }
}

sb_curve_view_data <- function(x, kind, cohort_label = NULL) {
  x <- sb_viz_result_input(x)
  if (is.list(x) && !is.object(x) && length(x) && !is.null(names(x))) {
    parts <- lapply(sb_viz_leaves(x, cohort_label), function(leaf) {
      value <- leaf$value
      attr(value, "stormbreaker_view_context") <- leaf$context
      sb_curve_view_data(value, kind, leaf$dataset)
    })
    data <- sb_view_bind(lapply(parts, `[[`, "data"))
    status <- sb_view_bind(lapply(parts, `[[`, "status"))
    data$panel <- factor(data$panel, levels = unique(status$panel))
    return(list(data = data, status = status))
  }
  if (sb_is_composite_view(x) || !is.null(attr(x, "stormbreaker_view_context"))) {
    return(sb_composite_curves(x, kind, cohort_label))
  }
  is_cv <- function(value) {
    inherits(value, c("stormbreaker_cv_result", "stormbreaker_binary_cv_result"))
  }
  is_external <- function(value) inherits(value, "stormbreaker_external_evaluation")
  if (inherits(x, "stormbreaker_cv_comparison")) {
    results <- x$results
    external <- FALSE
  } else if (inherits(x, "stormbreaker_external_comparison")) {
    results <- x$evaluations
    external <- TRUE
  } else if (is_cv(x)) {
    results <- list(x)
    external <- FALSE
  } else if (is_external(x)) {
    results <- list(x)
    external <- TRUE
  } else {
    sb_abort(
      "Supply a CV result/comparison or an external evaluation/comparison.",
      "stormbreaker_bad_curve_view"
    )
  }
  predicate <- if (external) is_external else is_cv
  if (!is.list(results) || length(results) == 0L ||
      !all(vapply(results, predicate, logical(1)))) {
    sb_abort("The comparison contains no supported results.",
             "stormbreaker_bad_curve_view")
  }
  if (!is.null(cohort_label) &&
      (!is.character(cohort_label) || length(cohort_label) != 1L ||
       is.na(cohort_label) || !nzchar(trimws(cohort_label)))) {
    sb_abort("`cohort_label` must be NULL or one non-empty string.",
             "stormbreaker_bad_cohort_label")
  }
  boundary <- if (external && all(vapply(results, function(result) {
    identical(result$assessment_role, "internal_holdout")
  }, logical(1)))) "Internal holdout" else if (external) "External evaluation" else "Development CV"
  cohort <- if (is.null(cohort_label)) boundary else cohort_label
  value_or <- function(value, default) if (is.null(value)) default else value

  parts <- lapply(results, function(result) {
    task_type <- if (external) {
      result$task_type
    } else if (inherits(result, "stormbreaker_binary_cv_result")) {
      "binary"
    } else {
      "survival"
    }
    pipeline_id <- if (external) result$pipeline_id else result$pipeline$id
    descriptor <- if (external) result$endpoint else result$task
    key <- if (identical(kind, "dca")) {
      "decision_curve"
    } else if (identical(task_type, "binary")) {
      "moderate_calibration"
    } else {
      "calibration"
    }
    metric <- result$metrics[[key]]
    horizon <- value_or(metric$prediction_horizon,
                        value_or(descriptor$prediction_horizon, NA_real_))
    unit <- if (identical(task_type, "survival")) {
      value_or(metric$time_unit, value_or(descriptor$time_unit, NA_character_))
    } else {
      NA_character_
    }
    status <- if (!is.null(metric)) {
      metric$status
    } else if (identical(result$status, "success") &&
               (identical(kind, "dca") ||
                (identical(task_type, "survival") && is.na(horizon)))) {
      "not_requested"
    } else {
      "not_saved"
    }
    scope <- value_or(metric$scope, NA_character_)
    target <- if (identical(task_type, "survival")) {
      if (is.na(horizon)) "horizon not declared" else paste(horizon, unit)
    } else {
      "class 1"
    }
    panel <- paste0(
      pipeline_id, "\n", cohort,
      if (!identical(cohort, boundary)) paste0(" | ", boundary) else "",
      " | ", target, "\ncurve: ", status, " | run: ", result$status
    )
    curve <- metric$curve
    n_rows <- if (is.data.frame(curve)) nrow(curve) else 0L
    if (n_rows == 0L) {
      coords <- data.frame(
        source_row = NA_integer_, strategy = "model", x = NA_real_, y = NA_real_,
        row_usable = FALSE, stringsAsFactors = FALSE
      )
    } else if (identical(kind, "dca")) {
      columns <- c("threshold_probability", "model_net_benefit",
                   "treat_all_net_benefit", "treat_none_net_benefit")
      if (!all(columns %in% names(curve))) {
        sb_abort("The saved DCA curve lacks its documented coordinate columns.",
                 "stormbreaker_bad_saved_curve")
      }
      row_usable <- if ("status" %in% names(curve)) {
        !is.na(curve$status) & curve$status == "success"
      } else {
        rep(TRUE, n_rows)
      }
      coords <- do.call(rbind, lapply(c("model", "treat_all", "treat_none"),
                                     function(strategy) {
        data.frame(
          source_row = seq_len(n_rows), strategy = strategy,
          x = curve$threshold_probability,
          y = curve[[paste0(strategy, "_net_benefit")]],
          row_usable = row_usable, stringsAsFactors = FALSE
        )
      }))
    } else {
      columns <- if (identical(task_type, "binary")) {
        c("probability", "calibrated_probability")
      } else {
        c("predicted_event_probability", "calibrated_event_probability")
      }
      if (!all(columns %in% names(curve))) {
        sb_abort("The saved calibration curve lacks its documented coordinate columns.",
                 "stormbreaker_bad_saved_curve")
      }
      coords <- data.frame(
        source_row = seq_len(n_rows), strategy = "model",
        x = curve[[columns[[1]]]], y = curve[[columns[[2]]]],
        row_usable = TRUE, stringsAsFactors = FALSE
      )
    }
    coords$draw <- identical(status, "success") &&
      identical(result$status, "success")
    coords$draw <- coords$draw & coords$row_usable &
      is.finite(coords$x) & is.finite(coords$y)
    coords$segment <- cumsum(!coords$draw)
    info <- data.frame(
      pipeline_id = pipeline_id, task_type = task_type,
      validation = if (external) "external_evaluation" else "development_cv",
      cohort = cohort, curve_key = key, metric_status = status,
      run_status = result$status, scope = scope,
      prediction_horizon = horizon, time_unit = unit,
      error_class = value_or(metric$error_class,
                             value_or(result$error_class, NA_character_)),
      error_message = value_or(metric$error_message,
                               value_or(result$error_message, NA_character_)),
      n_curve_rows = as.integer(n_rows), n_drawn_coordinates = as.integer(sum(coords$draw)),
      panel = panel, stringsAsFactors = FALSE
    )
    data <- cbind(info[rep(1L, nrow(coords)), , drop = FALSE], coords)
    list(data = data, status = info)
  })
  data <- do.call(rbind, lapply(parts, `[[`, "data"))
  status <- do.call(rbind, lapply(parts, `[[`, "status"))
  rownames(data) <- rownames(status) <- NULL
  data$panel <- factor(data$panel, levels = status$panel)
  data$line_group <- interaction(data$pipeline_id, data$strategy, data$segment,
                                 drop = TRUE)
  list(data = data, status = status)
}

sb_curve_view_lines <- function(data) {
  rows <- data[data$draw, , drop = FALSE]
  counts <- table(rows$line_group)
  rows[as.integer(counts[as.character(rows$line_group)]) >= 2L, , drop = FALSE]
}

sb_curve_view_annotations <- function(view) {
  info <- view$status
  info <- info[info$n_drawn_coordinates == 0L, , drop = FALSE]
  info$label <- if (nrow(info) > 0L) {
    paste0("No curve\n", info$metric_status, " / run: ", info$run_status)
  } else {
    character()
  }
  info$panel <- factor(info$panel, levels = levels(view$data$panel))
  info
}

sb_curve_view_scope_label <- function(view) {
  scopes <- unique(view$status$scope[!is.na(view$status$scope)])
  if (length(scopes) == 0L) "No saved curve scope" else paste(scopes, collapse = " | ")
}

# ---- plot-clinician-panel ----

#' @rdname sb_plot_calibration
#' @export
sb_plot_clinician_panel <- function(x, views, ncol = 2L) {
  if (!inherits(x, "stormbreaker_clinician_views") || !is.character(views) ||
      !length(views) || anyNA(views) || anyDuplicated(views) ||
      !all(views %in% names(x$plots))) sb_abort(
    "Choose unique available plot names from names(x$plots); inspect x$status for unavailable views.",
    "stormbreaker_bad_report_input")
  if (!is.numeric(ncol) || length(ncol) != 1L || !is.finite(ncol) ||
      ncol < 1 || ncol != floor(ncol)) sb_abort(
    "ncol must be a positive integer.", "stormbreaker_bad_report_input")
  if (!requireNamespace("patchwork", quietly = TRUE)) sb_abort(
    "Install patchwork to combine clinical figures.", "stormbreaker_missing_dependency")
  indices <- match(views, names(x$plots))
  captions <- vapply(seq_along(views), function(i) {
    p <- x$plots[[indices[i]]]
    notes <- unique(c(p$labels$title, p$labels$subtitle, p$labels$caption, x$captions[indices[i]]))
    paste(notes[!is.na(notes) & nzchar(notes) & notes != views[i]], collapse = "\n")
  }, character(1))
  plots <- lapply(seq_along(views), function(i) {
    x$plots[[indices[i]]] + ggplot2::labs(title = views[i], subtitle = NULL, caption = NULL)
  })
  result <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(tag_levels = "A")
  attr(result, "stormbreaker_captions") <- data.frame(view = views, caption = captions)
  result
}

# ---- probability-distribution ----

#' @rdname sb_evaluate_external
#' @export
sb_probability_distribution <- function(x, breaks = seq(0, 1, .1),
                                        cohort_label = NULL) {
  if (!is.numeric(breaks) || length(breaks) < 2L || any(!is.finite(breaks)) ||
      breaks[1] != 0 || utils::tail(breaks, 1) != 1 || any(diff(breaks) <= 0)) sb_abort(
    "Supply strictly increasing bin boundaries from zero to one.",
    "stormbreaker_bad_visualization_argument")
  thresholds <- sb_threshold_metrics(x, breaks, cohort_label)
  tab <- thresholds$table
  rows <- list()
  if (nrow(tab)) {
    for (start in seq.int(1L, nrow(tab), by = length(breaks))) {
      block <- tab[start + seq_along(breaks) - 1L, , drop = FALSE]
      for (observed in c(0L, 1L)) {
        cumulative <- if (observed == 1L) block$tp else block$fp
        # Last bin includes probability == 1 instead of subtracting it.
        count <- if (length(breaks) == 2L) cumulative[1L] else
          utils::head(cumulative, -1L) - c(cumulative[seq.int(2L, length(breaks) - 1L)], 0)
        metadata <- block[rep(1L, length(breaks) - 1L),
          names(thresholds$diagnostics), drop = FALSE]
        metadata$lower <- utils::head(breaks, -1L)
        metadata$upper <- utils::tail(breaks, -1L)
        metadata$observed_class <- observed
        metadata$count <- count
        metadata$class_n <- cumulative[1L]
        metadata$fraction <- if (cumulative[1L] == 0) NA_real_ else count / cumulative[1L]
        rows[[length(rows) + 1L]] <- metadata
      }
    }
  }
  structure(list(table = if (length(rows)) sb_view_bind(rows) else data.frame(),
    diagnostics = thresholds$diagnostics, breaks = breaks,
    bin_rule = "[lower, upper), except final bin [lower, 1]",
    saved_metrics = thresholds$saved_metrics),
    class = "stormbreaker_probability_distribution")
}

#' @rdname sb_plot_calibration
#' @export
sb_plot_probability_distribution <- function(x, scale = c("count", "fraction"),
                                             ncol = 2L) {
  sb_require_curve_plotting()
  if (!inherits(x, "stormbreaker_probability_distribution")) sb_abort(
    "Supply an sb_probability_distribution result.", "stormbreaker_bad_visualization_argument")
  scale <- match.arg(scale)
  ncol <- sb_viz_positive_integer(ncol, "ncol")
  d <- x$table
  if (!nrow(d)) {
    p <- sb_empty_display_plot("Prediction probability distribution",
      "No usable binary predictions.", "Unavailable groups remain in diagnostics.")
  } else {
    d$value <- d[[scale]]
    midpoint <- (d$lower + d$upper) / 2
    gap <- (d$upper - d$lower) * .03
    d$xmin <- ifelse(d$observed_class == 0, d$lower + gap, midpoint + gap)
    d$xmax <- ifelse(d$observed_class == 0, midpoint - gap, d$upper - gap)
    d$class <- factor(d$observed_class, levels = 0:1, labels = c("Observed 0", "Observed 1"))
    xmin <- xmax <- value <- class <- panel <- NULL
    p <- ggplot2::ggplot(d, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = value, fill = class)) +
      ggplot2::geom_rect(na.rm = TRUE) +
      ggplot2::facet_wrap(~panel, ncol = ncol, labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
      ggplot2::scale_fill_manual(values = c("#2878A0", "#D77935"), drop = FALSE) +
      ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
      ggplot2::labs(title = "Prediction probability by observed outcome",
        x = "Predicted probability of class 1",
        y = if (scale == "count") "Patients per bin" else "Fraction within observed class",
        fill = NULL,
        caption = "Adjacent bars share each probability bin. Last bin includes 1.\nSeparate folds/cohorts; no smoothing. Absent-class fractions are undefined.") +
      sb_viz_theme()
  }
  attr(p, "stormbreaker_probability_diagnostics") <- x$diagnostics
  p
}

# ---- plot-external-uncertainty ----

#' @rdname sb_bootstrap_external
#' @export
sb_plot_external_intervals <- function(x, metric, cohort_label = NULL) {
  sb_plot_bootstrap_intervals(x, metric, cohort_label, paired = FALSE)
}

#' @rdname sb_bootstrap_external
#' @export
sb_plot_external_difference <- function(x, metric, cohort_label = NULL) {
  sb_plot_bootstrap_intervals(x, metric, cohort_label, paired = TRUE)
}

sb_viz_uncertainty_inputs <- function(x, classes, cohort_label) {
  if (inherits(x, classes)) {
    label <- sb_result_cohort_label(cohort_label)
    return(stats::setNames(list(x), if (is.na(label)) "Unnamed cohort" else label))
  }
  if (!is.null(cohort_label) || !is.list(x) || !length(x) || is.null(names(x)) ||
      anyNA(names(x)) || any(!nzchar(trimws(names(x)))) || anyDuplicated(names(x)) ||
      !all(vapply(x, inherits, logical(1), what = classes))) {
    sb_abort("Supply a supported uncertainty object or a uniquely named list of them.",
      "stormbreaker_bad_visualization_input")
  }
  x
}

sb_viz_uncertainty_target <- function(x, metric) {
  target <- x$metric_target
  if (is.null(target)) target <- x$prediction_target
  horizon_metric <- !metric %in% c("roc_auc", "average_precision", "brier_score", "log_loss", "harrell_c") &&
    identical(x$task_type, "survival")
  label <- if (!horizon_metric) "horizon-free target" else if (isTRUE(target$has_horizon)) {
    paste(target$prediction_horizon, target$time_unit)
  } else "horizon not saved"
  paste(x$task_type, metric, label, sep = " | ")
}

sb_viz_empty_interval <- function(metric, status = "interval_not_saved") {
  data.frame(metric = metric, point_estimate = NA_real_, lower = NA_real_, upper = NA_real_,
    confidence_level = NA_real_, status = status, interval_available = FALSE,
    replicate_completeness = "not_saved", n_requested = NA_integer_,
    n_finite_successful = NA_integer_, stringsAsFactors = FALSE)
}

sb_viz_interval_label <- function(tab) {
  level <- ifelse(is.finite(tab$confidence_level), paste0(format(100 * tab$confidence_level, trim = TRUE), "% CI"), "CI not saved")
  counts <- ifelse(is.finite(tab$n_requested) & is.finite(tab$n_finite_successful),
    paste0(tab$n_finite_successful, "/", tab$n_requested, " draws; ", tab$replicate_completeness),
    "replicate count not saved")
  paste(level, counts, tab$status, sep = " | ")
}

sb_plot_bootstrap_intervals <- function(x, metric, cohort_label, paired) {
  sb_require_curve_plotting()
  if (!is.character(metric) || length(metric) != 1L || is.na(metric) || !nzchar(metric)) {
    sb_abort("Explicitly select one saved metric.", "stormbreaker_bad_result_metric")
  }
  classes <- if (paired) "stormbreaker_external_bootstrap_comparison" else "stormbreaker_external_bootstrap"
  values <- sb_viz_uncertainty_inputs(x, classes, cohort_label)
  rows <- lapply(seq_along(values), function(i) {
    value <- values[[i]]
    all <- value$intervals
    tab <- if (is.data.frame(all) && "metric" %in% names(all)) all[!is.na(all$metric) & all$metric == metric, , drop = FALSE] else NULL
    if (!is.null(tab) && nrow(tab) > 0L) {
      required <- c("point_estimate", "lower", "upper", "confidence_level", "status", "interval_available",
        "replicate_completeness", "n_requested", "n_finite_successful")
      if (!all(required %in% names(tab))) sb_abort("Saved interval fields are incomplete.", "stormbreaker_bad_visualization_input")
    }
    if (paired) {
      candidates <- value$candidate_status
      ids <- candidates$candidate_pipeline_id
      if (!length(ids) || anyNA(ids) || anyDuplicated(ids)) {
        sb_abort("A paired comparison needs distinct saved candidate IDs.", "stormbreaker_bad_visualization_input")
      }
      pieces <- lapply(seq_along(ids), function(j) {
        row <- if (!is.null(tab) && nrow(tab)) tab[tab$candidate_pipeline_id == ids[[j]], , drop = FALSE] else NULL
        if (is.null(row) || !nrow(row)) {
          row <- sb_viz_empty_interval(metric, if (candidates$status[[j]] == "failed") "failed" else "interval_not_saved")
          row$candidate_pipeline_id <- ids[[j]]
        }
        if (nrow(row) != 1L) sb_abort("Duplicate candidate/metric interval rows.", "stormbreaker_bad_visualization_input")
        row$reference_pipeline_id <- value$reference_pipeline_id
        row$pipeline_id <- ids[[j]]
        row
      })
      tab <- sb_view_bind(pieces)
    } else {
      if (is.null(tab) || !nrow(tab)) tab <- sb_viz_empty_interval(metric,
        if (identical(value$status, "failed")) "failed" else "interval_not_saved")
      if (nrow(tab) != 1L) sb_abort("Duplicate model/metric interval rows.", "stormbreaker_bad_visualization_input")
      tab$pipeline_id <- value$pipeline_id
    }
    tab$dataset <- names(values)[[i]]
    tab$target_label <- sb_viz_uncertainty_target(value, metric)
    tab$interval_type <- "External patient-bootstrap CI"
    tab$row_label <- paste(tab$dataset, tab$pipeline_id,
      if (paired) paste("Reference:", value$reference_pipeline_id) else "Locked model",
      sb_viz_interval_label(tab), sep = "\n")
    tab
  })
  data <- sb_view_bind(rows)
  sb_viz_interval_plot(data,
    title = if (paired) "Paired external performance differences" else "External performance uncertainty",
    xlab = if (paired) paste0(metric, " (candidate - reference)") else metric,
    caption = paste("Conditional on locked predictions; no development-fitting uncertainty or model selection.",
      if (paired) "Same patient draws within each comparison; positive favours candidates only for higher-is-better metrics." else
        "Saved bootstrap intervals and replicate completeness; no new resampling."),
    reference = if (paired) 0 else NULL)
}

sb_viz_interval_plot <- function(data, title, xlab, caption, reference = NULL) {
  data$point_available <- is.finite(data$point_estimate) & !data$status %in% c("failed", "failed_run")
  data$draw_interval <- !is.na(data$interval_available) & data$interval_available &
    is.finite(data$lower) & is.finite(data$upper) & data$lower <= data$upper & data$point_available
  data$row_label <- factor(data$row_label, levels = rev(unique(data$row_label)))
  data$target_label <- factor(data$target_label, levels = unique(data$target_label))
  lower <- upper <- row_label <- point_estimate <- interval_type <- target_label <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(y = row_label)) +
    ggplot2::geom_segment(data = data[data$draw_interval, , drop = FALSE],
      ggplot2::aes(x = lower, xend = upper, yend = row_label, linetype = interval_type),
      colour = "#24547A", linewidth = 0.7) +
    ggplot2::geom_point(data = data[data$point_available, , drop = FALSE],
      ggplot2::aes(x = point_estimate, shape = interval_type), colour = "#24547A", size = 2.4)
  if (!is.null(reference)) p <- p + ggplot2::geom_vline(xintercept = reference, linetype = 2, colour = "grey60")
  # Train row labels from all saved rows without adding estimates or marks.
  p + ggplot2::geom_blank() +
    ggplot2::facet_wrap(~target_label, ncol = 1, scales = "free",
      labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
    ggplot2::scale_y_discrete(labels = sb_view_wrap_label) +
    ggplot2::labs(title = title, subtitle = if (!any(data$point_available)) "No usable point estimates; unavailable rows retained." else NULL,
      x = xlab, y = NULL, linetype = "Interval", shape = "Interval", caption = sb_viz_wrap(caption)) + sb_viz_theme()
}

#' @rdname sb_meta_external_discrimination
#' @export
sb_plot_meta_forest <- function(x, show_prediction_interval = TRUE) {
  sb_require_curve_plotting()
  classes <- c("stormbreaker_external_discrimination_meta", "stormbreaker_external_mean_calibration_meta",
    "stormbreaker_external_calibration_slope_meta", "stormbreaker_external_discrimination_difference_meta",
    "stormbreaker_external_prediction_error_difference_meta")
  if (!inherits(x, classes) || !identical(x$status, "success")) {
    sb_abort("Supply a successful saved external Meta-analysis object.", "stormbreaker_bad_visualization_input")
  }
  if (!is.logical(show_prediction_interval) || length(show_prediction_interval) != 1L || is.na(show_prediction_interval)) {
    sb_abort("`show_prediction_interval` must be TRUE or FALSE.", "stormbreaker_bad_visualization_argument")
  }
  paired <- inherits(x, classes[4:5])
  cohorts <- x$cohorts
  sources <- if (paired) x$comparisons else x$bootstraps
  if (!is.data.frame(cohorts) || !nrow(cohorts) || anyDuplicated(cohorts$cohort_id)) {
    sb_abort("Saved cohort rows are missing or duplicated.", "stormbreaker_bad_visualization_input")
  }
  rows <- lapply(seq_len(nrow(cohorts)), function(i) {
    id <- cohorts$cohort_id[[i]]
    source <- sources[[id]]
    intervals <- source$intervals
    row <- NULL
    if (is.data.frame(intervals) && nrow(intervals)) {
      take <- !is.na(intervals$metric) & intervals$metric == x$metric
      if (paired) take <- take & intervals$candidate_pipeline_id == x$candidate_pipeline_id
      row <- intervals[take, , drop = FALSE]
    }
    if (is.null(row) || !nrow(row)) row <- sb_viz_empty_interval(x$metric)
    if (nrow(row) != 1L) sb_abort("Multiple saved intervals for one cohort and target.", "stormbreaker_bad_visualization_input")
    row$point_estimate <- if (paired) cohorts$point_difference[[i]] else cohorts$point_estimate[[i]]
    row$cohort_id <- id
    row$interval_type <- "Cohort bootstrap CI"
    row$row_label <- paste(id, sb_viz_interval_label(row), sep = "\n")
    row
  })
  summary <- x$summary
  if (!is.data.frame(summary) || nrow(summary) != 1L) sb_abort("One saved pooled summary is required.", "stormbreaker_bad_visualization_input")
  pooled <- sb_viz_empty_interval(x$metric, "success")
  pooled$point_estimate <- summary$estimate
  pooled$lower <- summary$lower
  pooled$upper <- summary$upper
  pooled$confidence_level <- x$confidence_level
  pooled$interval_available <- all(is.finite(c(pooled$lower, pooled$upper)))
  pooled$interval_type <- "Pooled CI"
  pooled$row_label <- paste0("Pooled estimate\n", 100 * x$confidence_level, "% CI")
  rows[[length(rows) + 1L]] <- pooled
  if (show_prediction_interval) {
    pi <- pooled
    pi$lower <- summary$prediction_lower
    pi$upper <- summary$prediction_upper
    pi$interval_available <- all(is.finite(c(pi$lower, pi$upper)))
    pi$interval_type <- "Prediction interval"
    pi$row_label <- paste0("Between-cohort prediction interval\n", 100 * x$confidence_level, "% PI")
    rows[[length(rows) + 1L]] <- pi
  }
  data <- sb_view_bind(rows)
  data$target_label <- sb_viz_uncertainty_target(x, x$metric)
  xlab <- if (paired) paste0(x$metric, " (candidate - reference)") else x$metric
  reference <- if (paired) 0 else if (x$metric %in% c("mean_calibration_ratio", "calibration_slope")) 1 else NULL
  sb_viz_interval_plot(data, "External validation forest plot", xlab,
    paste("Saved cohort bootstrap CIs, pooled CI and between-cohort PI; no refitting.",
      x$prediction_interval_caution), reference)
}
