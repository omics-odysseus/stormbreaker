#' Discrimination and comparison plots
#' @name sb_plot_performance
#' @description Plot saved discrimination, confusion matrices and model-by-cohort comparisons.
#' @rawRd \arguments{
#' \item{x}{
#' \strong{For \code{sb_plot_performance()}:} A CV result/comparison or external evaluation/single-cohort
#' comparison, or whole repeated CV, IECV or stacking CV result. Contexts
#' have distinct labels, without pooling. See \code{\link[=sb_result_table]{sb_result_table()}}.
#' 
#' \strong{For \code{sb_plot_performance_heatmap()}:} A saved result/comparison, supported composite CV object, or a
#' uniquely named list of these; names identify cohorts or result sets.
#' 
#' \strong{For \code{sb_plot_confusion()}:} A result from \code{\link[=sb_threshold_metrics]{sb_threshold_metrics()}}.
#' 
#' \strong{For \code{sb_plot_roc()}:} A binary CV or external evaluation/comparison, a supported
#' repeated/IECV/stacking CV object, or a uniquely named list of these.
#' }
#' \item{metric}{
#' \strong{For \code{sb_plot_performance()}:} One supported scalar metric, as documented in
#' \code{\link[=sb_result_table]{sb_result_table()}}. It must already exist in every supplied result.
#' Failed external evaluations without a saved metric remain in diagnostics.
#' 
#' \strong{For \code{sb_plot_performance_heatmap()}:} A scalar metric supported by \code{\link[=sb_result_table]{sb_result_table()}}.
#' 
#' \strong{For \code{sb_plot_model_ranking()}:} "auc" (default) or "c_index"; supply the matching numeric
#' column. All C-index cells must use the same estimator (e.g. Harrell),
#' endpoint and comparable evaluation scope; do not mix Harrell and Uno.
#' 
#' \strong{For \code{sb_plot_cohort_metrics()}:} Metric column to display.
#' }
#' \item{cohort_label}{
#' \strong{For \code{sb_plot_performance()}:} Optional non-empty display label supplied by the caller,
#' forwarded to \code{\link[=sb_result_table]{sb_result_table()}}. No cohort identity is inferred.
#' 
#' \strong{For \code{sb_plot_performance_heatmap()}:} Optional cohort label for one object, not a named list.
#' 
#' \strong{For \code{sb_plot_roc()}:} Optional cohort display label for a single object.
#' For a named list, labels come from the list names.
#' }
#' \item{show_unavailable}{
#' Logical; default FALSE excludes unavailable pipeline
#' labels. TRUE retains the previous audit view with labels but no points.
#' }
#' \item{digits}{
#' Number of decimal places printed in cells (zero to six).
#' }
#' \item{annotate}{
#' \strong{For \code{sb_plot_performance_heatmap()}:} Logical; print estimates or unavailable status in cells.
#' 
#' \strong{For \code{sb_plot_roc()}:} Add displayed-curve AUC labels to ROC panels.
#' }
#' \item{performance}{
#' \strong{For \code{sb_plot_model_ranking()}:} Data frame with model, cohort and auc columns; one row per
#' model/cohort. Failed or missing AUCs must be NA, never zero.
#' 
#' \strong{For \code{sb_plot_cohort_metrics()}:} Data frame with unique cohort names and a numeric column
#' named c_index or auc. Supply one model and one time horizon at a time.
#' }
#' \item{ranking_cohorts}{
#' Explicit cohort names defining an unweighted mean AUC.
#' Use development or designated validation cohorts, not a held-out final test.
#' }
#' \item{top_n}{
#' Maximum number of complete models to display; Inf shows all.
#' }
#' \item{title}{
#' Optional short title, for example '1-year AUC'.
#' }
#' \item{threshold}{
#' \strong{For \code{sb_plot_confusion()}:} One threshold already evaluated in x.
#' 
#' \strong{For \code{sb_plot_prediction_confusion()}:} Fixed decision threshold, default 0.5. Probability greater
#' than or equal to threshold predicts class 1. No threshold optimization.
#' }
#' \item{normalize}{
#' \strong{For \code{sb_plot_confusion()}:} Show counts or proportions within each observed class.
#' 
#' \strong{For \code{sb_plot_prediction_confusion()}:} "none" for counts, "observed" for proportions within each
#' observed class. A missing observed class has undefined proportions.
#' }
#' \item{ncol}{
#' \strong{For \code{sb_plot_confusion()}:} Panel columns.
#' 
#' \strong{For \code{sb_plot_roc()}:} Positive number of panel columns. Defaults to two.
#' }
#' \item{outcome}{
#' \strong{For \code{sb_plot_prediction_confusion()}:} Numeric vector of observed 0/1 labels in patient row order.
#' 
#' \strong{For \code{sb_plot_prediction_roc()}:} Numeric 0/1 outcomes, one independent patient per row.
#' }
#' \item{probability}{
#' \strong{For \code{sb_plot_prediction_confusion()}:} Numeric vector of aligned class-1 probabilities in \verb{[0,1]}.
#' Rows are matched by position, not vector names. No rows are dropped.
#' 
#' \strong{For \code{sb_plot_prediction_roc()}:} Numeric score vector or matrix with named model columns.
#' Scores need not be probabilities; higher values must indicate class 1.
#' }
#' \item{kind}{
#' \strong{For \code{sb_plot_prediction_confusion()}:} Prediction origin: "apparent" training, "external", or "oof".
#' This records the caller's declaration; it does not verify independence.
#' 
#' \strong{For \code{sb_plot_prediction_roc()}:} Explicit prediction origin: apparent training, external, or OOF.
#' }
#' \item{cohort}{
#' \strong{For \code{sb_plot_prediction_confusion()}:} Optional short cohort title.
#' 
#' \strong{For \code{sb_plot_prediction_roc()}:} Short cohort label.
#' }
#' \item{class_labels}{
#' Two distinct display labels for classes 0 and 1.
#' }
#' \item{ci_reps}{
#' \strong{For \code{sb_plot_prediction_roc()}:} Patient-bootstrap replicates; 0 omits intervals.
#' 
#' \strong{For \code{sb_plot_roc()}:} Optional patient-bootstrap replicates for ROC AUC intervals;
#' 0 omits intervals. Use 1000 for automatic 95 percent interval labels.
#' Requires independent patient rows, not repeated observations per patient.
#' }
#' \item{seed}{
#' \strong{For \code{sb_plot_prediction_roc()}:} Random seed.
#' 
#' \strong{For \code{sb_plot_roc()}:} Random seed for ROC intervals.
#' }
#' \item{by_fold}{
#' For ROC, show separate assessment folds instead of a single
#' pooled out-of-fold curve. Repetitions and external cohorts remain separate.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_plot_performance()}}{
#' Draw one selected scalar performance metric from saved results. This is
#' a descriptive display, not model selection or a new evaluation.
#' 
#' \strong{Returns:} A \code{ggplot} object. Its \code{data} retains every pipeline and the original
#' metric values, with \code{pipeline_id} represented as an input-ordered factor.
#' No file is written. The optional package \code{ggplot2} is required only when
#' this function is called.
#' 
#' \strong{Details:} Pipeline labels run from top to bottom in supplied order, not in
#' performance order. Points require both the selected metric and the complete
#' run to have status \code{"success"}, and the estimate to be finite. Otherwise
#' the pipeline is omitted from the default display, and the caption counts
#' omissions. If none are usable, an explicit no-results message is shown.
#' Inspect \code{\link[=sb_result_table]{sb_result_table()}} for metric and run status separately, or export
#' \code{\link[=sb_export_failures]{sb_export_failures()}} for saved failure details.
#' Estimates in the underlying table are never replaced with zero or clipped.
#' 
#' The subtitle reports the saved metric scope and aggregation. The plot does
#' distinguish external evaluation from development CV and displays a supplied
#' cohort label. An omitted external cohort label is explicitly unnamed. It does
#' not fit models, recompute metrics, infer confidence intervals, or rank
#' pipelines. The probability-scale metrics and concordance are displayed on
#' a zero-to-one axis, except that log loss and IPCW Brier have a zero lower
#' bound and no fixed upper bound. The IPCW weighting can yield estimates
#' above one, which remain visible. No error bars are drawn.
#' For \code{"time_auc"}, \code{"uno_c"}, and \code{"ipcw_brier"}, each pipeline label
#' includes its saved horizon and time unit. Distinct targets remain distinct
#' labels, not pooled or converted; they are not like-for-like comparisons.
#' These plots also retain a \code{display_label} column in their data.
#' Large comparisons may require a taller figure or fewer predeclared
#' candidate pipelines for readable labels.
#' }
#' 
#' \subsection{\code{sb_plot_performance_heatmap()}}{
#' Display one explicitly selected metric without averaging, ranking,
#' reevaluating patients or changing the saved estimates.
#' 
#' \strong{Returns:} A ggplot object. \code{data} retains original estimates and status,
#' a separate \code{fill_value} used only for successful cells, and explicit
#' \code{not_supplied} rows for missing model/cohort combinations. No patient
#' rows are included. Optional ggplot2 is required.
#' 
#' \strong{Details:} Rows retain model input order. Columns retain cohort/context input
#' order, including repetition and stacking role. Different tasks, horizons,
#' time units, metric scopes or aggregations are separated into panels.
#' Cohort names do not prove independent validation. No automatic winner or
#' clustered ordering is inferred. A cell needs successful metric and run
#' status and a finite estimate to be coloured. Unavailable cells are grey,
#' never zero; original values remain in the data. Colour scales use the
#' observed successful estimates and retain their numerical values; lighter
#' colours indicate preferable performance. No confidence interval is
#' inferred. Large grids require a suitably sized output image.
#' }
#' 
#' \subsection{\code{sb_plot_model_ranking()}}{
#' Rank models across explicitly selected cohorts and draw a Top-N panel
#' 
#' \strong{Returns:} List with full ranking, displayed data, heatmap, mean-AUC bar plot and
#' a combined patchwork plot. Models missing any ranking cohort are unranked;
#' other missing cells remain missing. Ties retain input model order.
#' }
#' 
#' \subsection{\code{sb_plot_cohort_metrics()}}{
#' Plot one model's performance across cohorts
#' 
#' \strong{Returns:} A ggplot. Estimates are supplied, not recalculated or ranked.
#' }
#' 
#' \subsection{\code{sb_plot_confusion()}}{
#' Plot a binary confusion matrix at an explicit threshold
#' 
#' \strong{Returns:} A ggplot object. All four cells are retained. A missing observed
#' class has undefined normalized proportions, not artificial zeros.
#' }
#' 
#' \subsection{\code{sb_plot_prediction_confusion()}}{
#' Plot a binary confusion matrix directly from labels and probabilities
#' 
#' \strong{Returns:} A ggplot, with observed classes on x and predicted classes on y.
#' The stormbreaker_prediction_context attribute records kind, threshold,
#' sample count and TP/FP/TN/FN. Plot data contain counts, not patient records.
#' 
#' \strong{Details:} Use a prespecified or development-derived threshold unchanged on
#' validation cohorts. Apparent training predictions are not cross-validation
#' predictions. This function does not fit models or enter automatic reports.
#' }
#' 
#' \subsection{\code{sb_plot_prediction_roc()}}{
#' Plot one binary ROC per supplied score column
#' 
#' \strong{Returns:} A ggplot with AUC/95 percent percentile intervals in its legend and
#' a stormbreaker_roc_statistics attribute. Intervals condition on supplied
#' scores, do not include model development uncertainty, and use ordinary
#' patient resampling. Single-class resamples are counted as unsuccessful.
#' }
#' 
#' \subsection{\code{sb_plot_roc()}}{
#' Compute descriptive threshold coordinates without fitting models or
#' replacing saved performance metrics. ROC defaults to one pooled out-of-fold
#' curve per model/context; precision-recall retains separate assessment folds.
#' 
#' \strong{Returns:} A ggplot object. Its data contains aggregate threshold coordinates,
#' not patient IDs. Attributes \code{stormbreaker_curve_status} and
#' \code{stormbreaker_saved_metrics} retain availability and unchanged saved
#' metrics for each entry. Optional ggplot2 is required.
#' 
#' \strong{Details:} Both functions treat class 1 as positive and higher probabilities
#' as more positive. Equal probabilities enter together. ROC joins the
#' resulting points linearly; PR uses vertical-then-horizontal steps so
#' areas correspond to noninterpolated average precision within a panel.
#' The starting PR precision at recall zero is defined as one.
#' 
#' Repetitions, held-out IECV cohorts and stacking roles remain separate.
#' External evaluation has no fictitious folds. Saved aggregate AUC/AP objects
#' are available as an attribute and are not replaced by displayed-curve AUC.
#' Pooled OOF is not apparent training performance. Bootstrap intervals
#' condition on saved predictions; they do not include development uncertainty.
#' There is no threshold selection or confidence band. A failed run or a
#' group with invalid/incomplete predictions is not
#' drawn; no invalid rows are silently removed. ROC needs both classes;
#' PR needs at least one positive. Unavailable groups remain in diagnostics.
#' Survival time-dependent ROC is not this binary interface.
#' }
#' }
#' @rawRd \seealso{
#' \code{\link[=sb_result_table]{sb_result_table()}}
#' }
#' @rawRd \examples{
#' # sb_plot_performance
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   d <- sb_example_survival(n = 60, seed = 1)
#'   task <- sb_surv_task(
#'     d, "patient_id", "time", "event", "gene_a", "month"
#'   )
#'   result <- sb_run_cv(task, folds = 3, seed = 1)
#'   sb_plot_performance(result, "harrell_c")
#' }
#' 
#' # sb_plot_performance_heatmap
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   d <- sb_example_binary(60, seed = 1)
#'   task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
#'   cv <- sb_run_cv(task, folds = 3, seed = 1)
#'   sb_plot_performance_heatmap(list(Development = cv), "roc_auc")
#' }
#' 
#' # sb_plot_roc
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   d <- sb_example_binary(60, seed = 1)
#'   task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
#'   cv <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1)
#'   sb_plot_roc(cv)
#'   sb_plot_pr(cv)
#' }
#' }
NULL

# ---- plot-comparison ----

# ---- plot-model-ranking ----

#' @rdname sb_plot_performance
#' @export
sb_plot_model_ranking <- function(performance, ranking_cohorts, top_n = 20L,
                                   metric = c("auc", "c_index")) {
  metric <- match.arg(metric)
  if (metric == "c_index") {
    if (!is.data.frame(performance) || !"c_index" %in% names(performance)) stop("Supply a c_index column.")
    d <- performance[c("model", "cohort", "c_index")]
    names(d)[3] <- "auc"
    result <- sb_plot_model_ranking(d, ranking_cohorts, top_n, metric = "auc")
    names(result$ranking)[names(result$ranking) == "mean_auc"] <- "mean_c_index"
    names(result$data)[names(result$data) == "auc"] <- "c_index"
    names(result$heatmap$data)[names(result$heatmap$data) == "auc"] <- "c_index"
    names(result$bar$data)[names(result$bar$data) == "mean_auc"] <- "mean_c_index"
    c_index <- mean_c_index <- NULL
    result$heatmap <- result$heatmap + ggplot2::aes(fill = c_index) + ggplot2::labs(fill = "C-index")
    result$bar <- result$bar + ggplot2::aes(x = mean_c_index) + ggplot2::labs(x = "Mean C-index")
    result$bar$layers[[2]]$mapping$label <- ggplot2::aes(label = sprintf("%.3f", mean_c_index))$label
    result$plot <- patchwork::wrap_plots(result$heatmap, result$bar, widths = c(4, 1.5))
    result$metric <- metric
    return(result)
  }
  sb_require_curve_plotting()
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install patchwork for the combined ranking plot.")
  d <- performance
  if (!is.data.frame(d) || !all(c("model", "cohort", "auc") %in% names(d)) || !nrow(d)) stop("Supply model, cohort and auc columns.")
  d$model <- as.character(d$model)
  d$cohort <- as.character(d$cohort)
  if (anyNA(d[c("model", "cohort")]) || any(!nzchar(d$model)) || any(!nzchar(d$cohort)) ||
      anyDuplicated(d[c("model", "cohort")]) || !is.numeric(d$auc) ||
      any(!is.na(d$auc) & (!is.finite(d$auc) | d$auc < 0 | d$auc > 1))) stop("Use unique model/cohort cells and AUC in [0,1] or NA.")
  if (!is.character(ranking_cohorts) || !length(ranking_cohorts) || anyNA(ranking_cohorts) ||
      anyDuplicated(ranking_cohorts) || !all(ranking_cohorts %in% d$cohort)) stop("Choose explicit unique ranking cohorts from the table.")
  if (length(top_n) != 1L || is.na(top_n) || top_n <= 0 || (is.finite(top_n) && top_n != floor(top_n))) stop("top_n must be a positive integer or Inf.")
  models <- unique(d$model)
  scores <- vapply(models, function(m) {
    z <- d[d$model == m, ]
    a <- z$auc[match(ranking_cohorts, z$cohort)]
    if (anyNA(a)) NA_real_ else mean(a)
  }, numeric(1))
  ranking <- data.frame(model = models, mean_auc = unname(scores), status = ifelse(is.na(scores), "incomplete", "ranked"))
  ranking <- ranking[order(-ranking$mean_auc, seq_along(models), na.last = TRUE), ]
  selected <- utils::head(ranking$model[ranking$status == "ranked"], if (is.infinite(top_n)) length(models) else top_n)
  if (!length(selected)) stop("No model has complete AUCs for every ranking cohort.")
  grid <- expand.grid(model = selected, cohort = unique(d$cohort), stringsAsFactors = FALSE)
  z <- merge(grid, d, by = c("model", "cohort"), all.x = TRUE, sort = FALSE)
  z$model <- factor(z$model, levels = rev(selected))
  z$cohort <- factor(z$cohort, levels = unique(d$cohort))
  z$label <- ifelse(is.na(z$auc), "", sprintf("%.3f", z$auc))
  bars <- ranking[match(selected, ranking$model), ]
  bars$model <- factor(bars$model, levels = rev(selected))
  model <- cohort <- auc <- label <- mean_auc <- NULL
  heatmap <- ggplot2::ggplot(z, ggplot2::aes(cohort, model, fill = auc)) +
    ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = label), size = 3) +
    ggplot2::scale_fill_gradient(low = "#eff3ff", high = "#6a51a3", limits = c(0, 1), na.value = "grey90") +
    ggplot2::labs(x = NULL, y = NULL, fill = "AUC") + sb_viz_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  bar <- ggplot2::ggplot(bars, ggplot2::aes(mean_auc, model)) +
    ggplot2::geom_col(fill = "#a63683", width = .85) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", mean_auc)), hjust = -.1, size = 3) +
    ggplot2::scale_x_continuous(limits = c(0, 1.12), breaks = c(0, .5, 1)) +
    ggplot2::labs(x = "Mean AUC", y = NULL) + sb_viz_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank())
  list(ranking = ranking, data = z, ranking_cohorts = ranking_cohorts,
    heatmap = heatmap, bar = bar, plot = patchwork::wrap_plots(heatmap, bar, widths = c(4, 1.5)))
}

# ---- plot-performance-heatmap ----

#' @rdname sb_plot_performance
#' @export
sb_plot_performance_heatmap <- function(x, metric, cohort_label = NULL,
                                        digits = 3L, annotate = TRUE) {
  sb_require_curve_plotting()
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits != floor(digits) || digits < 0 || digits > 6) {
    sb_abort("`digits` must be an integer from zero to six.", "stormbreaker_bad_visualization_argument")
  }
  if (!is.logical(annotate) || length(annotate) != 1L || is.na(annotate)) {
    sb_abort("`annotate` must be TRUE or FALSE.", "stormbreaker_bad_visualization_argument")
  }
  leaves <- sb_viz_leaves(x, cohort_label)
  rows <- lapply(leaves, function(leaf) {
    tab <- sb_result_table(leaf$value, metric = metric)
    tab <- sb_view_annotate(tab, leaf$context)
    tab$dataset <- leaf$dataset
    tab$column_label <- sb_viz_panel(leaf)
    horizon <- if (metric %in% sb_horizon_view_metrics()) {
      paste(tab$prediction_horizon, tab$time_unit)
    } else "horizon-free target"
    tab$target_label <- paste(tab$task_type, metric, horizon, tab$scope, tab$aggregation, sep = " | ")
    tab$cell_status <- ifelse(tab$run_status != "success" | is.na(tab$run_status), "failed_run",
      ifelse(tab$status != "success" | is.na(tab$status), "unavailable_metric",
        ifelse(is.finite(tab$estimate), "success", "nonfinite_estimate")))
    tab
  })
  data <- sb_view_bind(rows)
  model_order <- unique(data$pipeline_id)
  column_order <- unique(data$column_label)
  target_order <- unique(data$target_label)
  key <- paste(match(data$pipeline_id, model_order), match(data$column_label, column_order),
    match(data$target_label, target_order), sep = "/")
  if (anyDuplicated(key)) sb_abort("Multiple saved entries occupy the same model/cohort/target cell; give them distinct result-set labels.",
    "stormbreaker_ambiguous_heatmap_cell")
  additions <- list()
  for (target in target_order) {
    local <- data[data$target_label == target, , drop = FALSE]
    grid <- expand.grid(pipeline_id = unique(local$pipeline_id), column_label = unique(local$column_label),
      stringsAsFactors = FALSE)
    present <- paste(match(local$pipeline_id, model_order), match(local$column_label, column_order), sep = "/")
    missing <- !paste(match(grid$pipeline_id, model_order), match(grid$column_label, column_order), sep = "/") %in% present
    if (!any(missing)) next
    absent <- grid[missing, , drop = FALSE]
    absent$target_label <- target
    absent$cell_status <- "not_supplied"
    absent$estimate <- NA_real_
    absent$metric <- metric
    additions[[length(additions) + 1L]] <- absent
  }
  if (length(additions)) data <- sb_view_bind(c(list(data), additions))
  data$fill_value <- ifelse(data$cell_status == "success", data$estimate, NA_real_)
  data$annotation <- ifelse(data$cell_status == "success", formatC(data$estimate, format = "f", digits = digits),
    "")
  labels <- c(failed_run = "failed", unavailable_metric = "unavailable", nonfinite_estimate = "nonfinite", not_supplied = "not supplied")
  unavailable <- data$cell_status != "success"
  data$annotation[unavailable] <- unname(labels[data$cell_status[unavailable]])
  data$pipeline_id <- factor(data$pipeline_id, levels = rev(model_order))
  data$column_label <- factor(data$column_label, levels = column_order)
  data$target_label <- factor(data$target_label, levels = target_order)
  lower_better <- metric %in% c("brier_score", "log_loss", "ipcw_brier")
  pipeline_id <- column_label <- fill_value <- annotation <- target_label <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(x = column_label, y = pipeline_id, fill = fill_value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_viridis_c(option = "C", direction = if (lower_better) -1 else 1,
      na.value = "grey90", name = metric) +
    ggplot2::facet_wrap(~target_label, scales = "free", ncol = 1,
      labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
    ggplot2::scale_x_discrete(labels = sb_view_wrap_label) +
    ggplot2::scale_y_discrete(labels = sb_view_wrap_label) +
    ggplot2::labs(title = "Saved model-by-cohort performance",
      subtitle = if (all(unavailable)) "No usable successful estimates; availability is shown." else
        paste(if (lower_better) "Lower" else "Higher", "is better; input order, no model selection."),
      x = NULL, y = "Pipeline",
      caption = "Separate targets and validation contexts. Grey cells are unavailable, not zero. No pooling or inferred intervals.") +
    sb_viz_theme() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 9))
  if (annotate) p <- p + ggplot2::geom_label(ggplot2::aes(label = annotation),
    size = 2.8, fill = "white", alpha = 0.85)
  p
}

# ---- plot-performance ----

#' @rdname sb_plot_performance
#' @export
sb_plot_performance <- function(x, metric, cohort_label = NULL,
                                show_unavailable = FALSE) {
  sb_validate_show_unavailable(show_unavailable)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    sb_abort(
      "Plotting requires ggplot2; install it with install.packages('ggplot2').",
      "stormbreaker_missing_plot_dependency"
    )
  }
  tab <- sb_result_table(x, metric = metric, cohort_label = cohort_label)
  external <- inherits(x, c("stormbreaker_external_evaluation",
                             "stormbreaker_external_comparison"))
  tab$pipeline_id <- factor(tab$pipeline_id, levels = unique(tab$pipeline_id))
  horizon_metric <- metric %in% sb_horizon_view_metrics()
  composite <- "view_id" %in% names(tab)
  separate_labels <- horizon_metric || composite
  if (separate_labels) {
    tab$display_label <- if (horizon_metric) sb_horizon_view_labels(tab) else as.character(tab$pipeline_id)
    if (composite) {
      tab$display_label <- vapply(seq_len(nrow(tab)), function(i) {
        paste(tab$display_label[i], sb_view_label(as.list(tab[i, ])),
              paste("Saved entry:", tab$view_id[i]), sep = "\n")
      }, character(1))
    }
    tab$display_label <- factor(tab$display_label, levels = unique(tab$display_label))
  }
  drawable <- sb_performance_display_rows(tab)
  points <- tab[drawable, , drop = FALSE]
  omitted <- sum(!drawable)
  title <- if (external) "External performance" else "CV performance"
  subtitle <- paste(
    c(if (!is.null(cohort_label)) cohort_label else if (external)
      "External cohort (not named)",
      paste(unique(tab$scope), unique(tab$aggregation), sep = "\n", collapse = "\n")),
    collapse = "\n"
  )
  if (composite) {
    title <- "Performance by context"
    subtitle <- paste(subtitle, paste(unique(vapply(seq_len(nrow(tab)), function(i) {
      sb_view_boundary(as.list(tab[i, ]))
    }, character(1))), collapse = "\n"), "No pooling across repetitions or cohorts.", sep = "\n")
  }
  caption <- paste0(
    omitted, " pipeline(s) omitted from points (run/metric status or nonfinite estimate).\n",
    "Use sb_export_failures() for details. No confidence intervals."
  )
  if (horizon_metric) {
    subtitle <- paste(subtitle,
      "Saved horizon and unit shown per pipeline; different targets are not equivalent.",
      sep = "\n")
  }
  if (!show_unavailable && !any(drawable)) {
    p <- sb_empty_display_plot(title, subtitle, caption)
    p$data <- tab
    return(p)
  }
  metric_label <- switch(metric,
    roc_auc = "ROC AUC",
    average_precision = "Average precision",
    brier_score = "Brier score",
    log_loss = "Log loss",
    harrell_c = "Harrell C",
    time_auc = "Cumulative/dynamic time-AUC",
    uno_c = "Truncated Uno C",
    ipcw_brier = "IPCW Brier score"
  )
  pipeline_id <- estimate <- display_label <- NULL
  mapping <- if (separate_labels) ggplot2::aes(x = display_label, y = estimate) else
    ggplot2::aes(x = pipeline_id, y = estimate)
  labels <- if (separate_labels) tab$display_label else tab$pipeline_id
  visible_labels <- if (separate_labels) points$display_label else points$pipeline_id
  p <- ggplot2::ggplot(tab, mapping) +
    ggplot2::geom_point(data = points, size = 2.5, colour = "#24547A") +
    ggplot2::scale_x_discrete(
      limits = if (show_unavailable) rev(levels(labels)) else
        rev(unique(as.character(visible_labels))), drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      limits = if (metric %in% c("log_loss", "ipcw_brier")) c(0, NA_real_) else c(0, 1)
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "Model",
      y = metric_label,
      title = title, subtitle = subtitle, caption = caption
    ) +
    sb_viz_theme()
  sb_clinician_clean_plot(p, sb_clinician_models(sb_viz_leaves(x, cohort_label)))
}

# ---- plot-selection-stability ----

#' @rdname sb_analyze_methods
#' @export
sb_plot_selection_frequency <- function(x, source = c("selector", "glmnet_active_set"),
                                         top_n = 20L, cohort_label = NULL) {
  source <- match.arg(source)
  top_n <- sb_viz_positive_integer(top_n, "top_n")
  view <- sb_viz_selection_view(x, source, cohort_label)
  rows <- lapply(view$entries, function(entry) {
    if (!entry$usable) return(NULL)
    tab <- entry$summary$feature_frequency
    if (!is.data.frame(tab) || !nrow(tab)) return(NULL)
    ranking <- order(-tab$selection_frequency, seq_len(nrow(tab)), na.last = TRUE)
    chosen <- utils::head(ranking[is.finite(tab$selection_frequency[ranking])], top_n)
    tab <- tab[sort(chosen), , drop = FALSE]
    if (!nrow(tab)) return(NULL)
    tab$panel <- entry$panel
    tab$n_successful_folds <- entry$summary$n_successful_folds
    tab$n_folds <- entry$summary$n_folds
    tab$feature_key <- paste(entry$key, tab$feature, sep = "::")
    tab$display_order <- match(tab$feature, entry$summary$feature_frequency$feature[chosen])
    tab
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(sb_viz_empty_selection(view, "Feature selection frequency"))
  data <- sb_view_bind(rows)
  levels <- unlist(lapply(rows, function(row) row$feature_key[order(row$display_order)]), use.names = FALSE)
  data$feature_key <- factor(data$feature_key, levels = rev(levels))
  labels <- stats::setNames(data$feature, as.character(data$feature_key))
  selection_frequency <- feature_key <- panel <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(x = selection_frequency, y = feature_key)) +
    ggplot2::geom_col(fill = "#24547A", width = 0.7) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_discrete(labels = labels) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y", labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
    ggplot2::labs(title = "Feature selection frequency", x = "Fraction of successful selection folds", y = NULL,
      caption = paste0("Top ", top_n, " per panel by saved frequency; ties retain input order. Not biological importance. ",
        sum(!view$status$usable), " unavailable entry(s).")) + sb_viz_theme()
  attr(p, "stormbreaker_selection_status") <- view$status
  p
}

#' @rdname sb_analyze_methods
#' @export
sb_plot_selection_jaccard <- function(x, source = c("selector", "glmnet_active_set"), cohort_label = NULL) {
  source <- match.arg(source)
  view <- sb_viz_selection_view(x, source, cohort_label)
  rows <- lapply(view$entries, function(entry) {
    if (!entry$usable) return(NULL)
    tab <- entry$summary$pairwise
    if (!is.data.frame(tab) || !nrow(tab)) return(NULL)
    tab$panel <- entry$panel
    tab$fill_value <- ifelse(tab$status == "success" & is.finite(tab$jaccard), tab$jaccard, NA_real_)
    tab$annotation <- ifelse(is.finite(tab$fill_value), formatC(tab$jaccard, digits = 2, format = "f"), "NA")
    tab
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(sb_viz_empty_selection(view, "Fold-pair selection Jaccard"))
  data <- sb_view_bind(rows)
  fold_i <- fold_j <- fill_value <- annotation <- panel <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(x = factor(fold_j), y = factor(fold_i), fill = fill_value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_label(ggplot2::aes(label = annotation), size = 3, fill = "white", alpha = 0.85) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90", name = "Jaccard") +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free", labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
    ggplot2::labs(title = "Fold-pair selection Jaccard", x = "Fold j", y = "Fold i",
      caption = "Saved unique fold pairs only; no diagonal or missing pair is invented. Grey is unavailable, not zero.") +
    sb_viz_theme() + ggplot2::theme(panel.grid = ggplot2::element_blank())
  attr(p, "stormbreaker_selection_status") <- view$status
  p
}

#' @rdname sb_analyze_methods
#' @export
sb_plot_selection_upset <- function(x, source = c("selector", "glmnet_active_set"),
                                     max_intersections = 20L, cohort_label = NULL) {
  source <- match.arg(source)
  max_intersections <- sb_viz_positive_integer(max_intersections, "max_intersections")
  view <- sb_viz_selection_view(x, source, cohort_label)
  pieces <- lapply(view$entries, function(entry) {
    if (!entry$usable) return(NULL)
    records <- entry$records
    if (!is.list(records) || !length(records)) return(NULL)
    successful <- vapply(records, function(record) identical(record$status, "success"), logical(1))
    records <- records[successful]
    if (!length(records)) return(NULL)
    folds <- vapply(records, function(record) as.integer(record$fold), integer(1))
    if (anyNA(folds) || anyDuplicated(folds)) sb_abort("Saved selection folds are missing or duplicated.", "stormbreaker_bad_visualization_input")
    ord <- order(folds)
    records <- records[ord]
    folds <- folds[ord]
    sets <- lapply(records, function(record) unique(as.character(record$selected_features)))
    features <- unique(unlist(sets, use.names = FALSE))
    if (!length(features)) return(NULL)
    membership <- vapply(features, function(feature) paste(folds[vapply(sets, function(set) feature %in% set, logical(1))], collapse = ","), character(1))
    patterns <- unique(membership)
    groups <- lapply(patterns, function(pattern) features[membership == pattern])
    counts <- lengths(groups)
    ord <- order(-counts, seq_along(patterns))
    tab <- data.frame(panel = entry$panel, membership = patterns[ord], count = as.integer(counts[ord]),
      intersection = seq_along(ord), stringsAsFactors = FALSE)
    tab$features <- I(groups[ord])
    tab$shown <- tab$intersection <= max_intersections
    list(intersections = tab, folds = folds)
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(sb_viz_empty_selection(view, "Selection-set intersections"))
  intersections <- sb_view_bind(lapply(pieces, `[[`, "intersections"))
  displayed <- intersections[intersections$shown, , drop = FALSE]
  all_folds <- sort(unique(unlist(lapply(pieces, `[[`, "folds"))))
  n_sets <- length(all_folds)
  baseline <- n_sets + 1.5
  height <- max(3, n_sets * 0.65)
  displayed$bar_bottom <- baseline
  displayed$bar_top <- baseline + height * displayed$count / max(displayed$count)
  dots <- sb_view_bind(lapply(pieces, function(piece) {
    tab <- piece$intersections[piece$intersections$shown, , drop = FALSE]
    sb_view_bind(lapply(seq_len(nrow(tab)), function(i) {
      active <- as.integer(strsplit(tab$membership[[i]], ",", fixed = TRUE)[[1]])
      data.frame(panel = tab$panel[[i]], intersection = tab$intersection[[i]], fold = piece$folds,
        y = match(piece$folds, all_folds), active = piece$folds %in% active)
    }))
  }))
  links <- sb_view_bind(lapply(seq_len(nrow(displayed)), function(i) {
    active <- as.integer(strsplit(displayed$membership[[i]], ",", fixed = TRUE)[[1]])
    data.frame(panel = displayed$panel[[i]], intersection = displayed$intersection[[i]],
      bottom = min(match(active, all_folds)), top = max(match(active, all_folds)))
  }))
  intersection <- y <- active <- panel <- bar_bottom <- bar_top <- count <- bottom <- top <- NULL
  p <- ggplot2::ggplot(dots, ggplot2::aes(x = intersection, y = y)) +
    ggplot2::geom_segment(data = displayed, ggplot2::aes(x = intersection, xend = intersection,
      y = bar_bottom, yend = bar_top), inherit.aes = FALSE, linewidth = 5, colour = "#24547A") +
    ggplot2::geom_text(data = displayed, ggplot2::aes(x = intersection, y = bar_top, label = count),
      inherit.aes = FALSE, vjust = -0.5, size = 3) +
    ggplot2::geom_segment(data = links, ggplot2::aes(x = intersection, xend = intersection, y = bottom, yend = top),
      inherit.aes = FALSE, colour = "#24547A") +
    ggplot2::geom_point(ggplot2::aes(colour = active), size = 2.7) +
    ggplot2::scale_colour_manual(values = c("FALSE" = "grey85", "TRUE" = "#24547A"), guide = "none") +
    ggplot2::scale_y_continuous(breaks = seq_along(all_folds), labels = paste("Fold", all_folds),
      limits = c(0.5, baseline + height + 1)) +
    ggplot2::scale_x_continuous(breaks = sort(unique(displayed$intersection))) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_x", labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
    ggplot2::labs(title = "Selection-set intersections", x = "Exact intersection (size order)", y = NULL,
      caption = paste0("Exact membership among successful folds only; bars label feature counts. ",
        "Never-selected features are outside the union. ", sum(!intersections$shown),
        " intersection(s) omitted from display, retained in attributes.")) +
    sb_viz_theme() + ggplot2::theme(panel.grid = ggplot2::element_blank())
  attr(p, "stormbreaker_selection_status") <- view$status
  attr(p, "stormbreaker_intersections") <- intersections
  p
}

sb_viz_selection_view <- function(x, source, cohort_label) {
  sb_require_curve_plotting()
  leaves <- sb_viz_leaves(x, cohort_label)
  entries <- lapply(leaves, function(leaf) {
    value <- leaf$value
    if (inherits(value, "stormbreaker_external_evaluation")) sb_abort("External evaluations do not contain training-fold selection stability.",
      "stormbreaker_bad_visualization_task")
    summary <- if (source == "selector") value$selection_stability else value$glmnet_active_set_stability
    records <- if (source == "selector") value$selections else lapply(value$coefficients, function(record) {
      list(fold = record$fold, status = record$status, selected_features = record$nonzero_features)
    })
    if (is.null(summary)) summary <- list(status = "not_saved", n_successful_folds = 0L, n_folds = NA_integer_)
    usable <- !summary$status %in% c("not_saved", "not_applicable") && isTRUE(summary$n_successful_folds > 0)
    panel <- paste(sb_viz_panel(leaf), sb_viz_pipeline_id(value), source,
      paste0(summary$n_successful_folds, "/", summary$n_folds, " successful folds | ", summary$status), sep = "\n")
    list(summary = summary, records = records, usable = usable, panel = panel, key = leaf$context$view_id)
  })
  status <- sb_view_bind(lapply(entries, function(entry) data.frame(panel = entry$panel,
    status = entry$summary$status, n_successful_folds = entry$summary$n_successful_folds,
    n_folds = entry$summary$n_folds, usable = entry$usable)))
  list(entries = entries, status = status)
}

sb_viz_empty_selection <- function(view, title) {
  p <- sb_empty_display_plot(title, "No usable saved selection records for this view.",
    "Inspect stormbreaker_selection_status; missing or non-applicable sets are not invented.")
  p$data <- data.frame()
  attr(p, "stormbreaker_selection_status") <- view$status
  p
}
