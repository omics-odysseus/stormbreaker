#' Results, reports and exports
#' @name sb_render_report
#' @description Extract saved results, export tables and figures, and render offline reports.
#' @rawRd \arguments{
#' \item{x}{
#' \strong{For \code{sb_render_report()}:} One result supported by \code{\link[=sb_result_table]{sb_result_table()}}, or a uniquely named
#' list of such results. Whole repeated CV, IECV and stacking CV objects are
#' expanded into separate saved-entry sections with their context and role.
#' IECV is labeled development validation, not final independent validation.
#' Each ordinary list member is a separate report section;
#' names are caller-supplied cohort labels, not verified cohort identities.
#' 
#' \strong{For \code{sb_render_clinician_report()}:} A result from \code{\link[=sb_clinician_views]{sb_clinician_views()}}.
#' 
#' \strong{For \code{sb_export_run()} and \code{sb_run_figure_keys()}:}
#' Result from \code{\link[=sb_complete_run]{sb_complete_run()}}.
#' 
#' \strong{For \code{sb_export_failures()}:} A result or named result list accepted by \code{\link[=sb_failure_table]{sb_failure_table()}}.
#' 
#' \strong{For \code{sb_result_table()}:} A single result from \code{\link[=sb_run_cv]{sb_run_cv()}}, a shared-fold comparison from
#' \code{\link[=sb_compare_cv]{sb_compare_cv()}}, an evaluation from \code{\link[=sb_evaluate_external]{sb_evaluate_external()}}, or a
#' single-cohort comparison from \code{\link[=sb_compare_external]{sb_compare_external()}}. Also accepts whole
#' repeated CV, internal-external CV (IECV), and stacking CV results.
#' Arbitrary cross-cohort collections are not accepted by this table function.
#' 
#' \strong{For \code{sb_failure_table()}:} A CV result/comparison, external evaluation/comparison, whole
#' repeated CV, IECV or stacking CV result, or a uniquely
#' named list of these objects representing distinct report sections.
#' 
#' \strong{For \code{sb_view_results()}:} A saved CV/evaluation result, comparison, external collection,
#' \code{\link[=sb_compare_cohorts]{sb_compare_cohorts()}} result, or uniquely named list of saved results.
#' }
#' \item{file}{
#' \strong{For \code{sb_render_report()}:} A new \code{.html} file in an existing directory. A companion
#' \verb{<stem>_failures.csv} is written in the same directory. Both paths must be
#' new; existing files are refused and there is no overwrite option.
#' 
#' \strong{For \code{sb_render_clinician_report()}:} New HTML filename; existing outputs are refused.
#' 
#' \strong{For \code{sb_render_model_report()}:} New HTML filename in an existing directory. Existing files refused.
#' 
#' \strong{For \code{sb_render_visual_report()}:} New \code{.html} path in an existing directory. Existing files are
#' refused; there is no overwrite option.
#' 
#' \strong{For \code{sb_export_failures()}:} A new CSV path in an existing directory. Existing files are refused.
#' }
#' \item{title}{
#' \strong{For \code{sb_render_report()}:} One non-empty report title.
#' 
#' \strong{For \code{sb_render_clinician_report()}, \code{sb_render_model_report()}:} Report title.
#' 
#' \strong{For \code{sb_render_visual_report()}:} Non-empty report title.
#' 
#' \strong{For \code{sb_view_results()}:} Page title.
#' }
#' \item{complete}{
#' Expand all topics initially instead of only the first.
#' }
#' \item{fit}{
#' A successful supported linear final fit.
#' }
#' \item{new_data}{
#' Feature data for every sample to include; outcome columns
#' are ignored. All rows receive individual plots, without sampling.
#' }
#' \item{id_col}{
#' Unique identifier column in new_data.
#' }
#' \item{background_data}{
#' Original-training feature reference, held fixed for
#' every SHAP view and coefficient-SD calculation. The caller is responsible
#' for supplying the original training reference, not explained test patients.
#' }
#' \item{development}{
#' Optional original survival development task, needed only
#' if the saved model has no full baseline. See \code{\link[=sb_add_survival_curve]{sb_add_survival_curve()}}.
#' }
#' \item{top_n}{
#' Maximum terms in global SHAP and waterfall displays. Remaining
#' waterfall terms are aggregated by the existing plotting implementation.
#' }
#' \item{selection_results}{
#' Optional saved development CV results for gene
#' selection rankings. No selection counts are inferred from the final fit.
#' }
#' \item{compact}{
#' \strong{For \code{sb_render_model_report()}:} Remove repeated plot titles, subtitles and long captions.
#' 
#' \strong{For \code{sb_render_visual_report()}:} Remove plot titles, subtitles and captions from the image and
#' retain their text in a separate report-body section. Axes, legends, data
#' and statistical layers stay.
#' }
#' \item{plots}{
#' Non-empty, uniquely named list of ggplot objects. Names become
#' section headings and alternative image text.
#' }
#' \item{width}{
#' Positive finite plot dimensions in inches, either one
#' value for all plots or one value per plot in list order.
#' }
#' \item{height}{
#' Positive finite plot dimensions in inches, either one
#' value for all plots or one value per plot in list order.
#' }
#' \item{groups}{
#' Optional character vector, one topic per plot in list order.
#' Topics appear in first-occurrence order as expandable sections. NULL keeps
#' the flat gallery. Grouping organizes existing plots; it does not compute
#' missing clinical or algorithm diagnostics.
#' }
#' \item{expanded}{
#' Topics initially expanded. NULL expands only the first topic;
#' character(0) collapses all. Supply all unique groups for complete detail.
#' }
#' \item{captions}{
#' Optional character vector of additional figure legends, one
#' per plot. Use this for diagnostic interpretation stored separately from
#' plot labels. All text is escaped, not interpreted as HTML.
#' }
#' \item{availability}{
#' Optional data frame of view availability and reasons,
#' displayed as a selectable table before the gallery.
#' }
#' \item{directory}{
#' A new output directory. Existing paths are never overwritten.
#' }
#' \item{figure_format}{
#' For \code{sb_export_run()}, \code{"pdf"} (default) or \code{"png"}.
#' }
#' \item{dpi}{
#' Positive integer resolution for PNG figures; defaults to 300. PDF is vector.
#' }
#' \item{font_family}{
#' Optional installed font family, e.g. \code{"Times New Roman"}, for ordinary
#' plot text. NULL preserves the original style. Explicitly mapped fonts and
#' embedded text or custom grobs are not rewritten.
#' }
#' \item{figure_keys}{
#' Exact keys from \code{sb_run_figure_keys(x)} in the desired order. NULL
#' exports all saved figures; \code{character(0)} exports none. Requested
#' figures keep all their pages if paginated. No absent figures are built.
#' }
#' \item{include_data}{
#' Export the workflow's tables and RDS results as well as the figure index
#' and notes (default TRUE). Set FALSE to leave out data files, including
#' patient-level prediction tables and saved explanations.
#' }
#' \item{cohort_label}{
#' \strong{For \code{sb_export_failures()}:} Optional label, forwarded to \code{\link[=sb_failure_table]{sb_failure_table()}}.
#' 
#' \strong{For \code{sb_result_table()}:} Optional non-empty display label supplied by the caller.
#' A cohort identity is never inferred. External tables always include this
#' column (missing when omitted); CV tables include it only when supplied.
#' 
#' \strong{For \code{sb_failure_table()}:} Optional non-empty cohort display label for a single
#' result or comparison. For a named list, leave NULL and use its names.
#' }
#' \item{component}{
#' One of \code{"metrics"}, \code{"predictions"}, or \code{"runs"}.
#' }
#' \item{metric}{
#' For \code{component = "metrics"}, one explicitly selected saved
#' scalar metric: \code{"roc_auc"}, \code{"average_precision"}, \code{"brier_score"},
#' \code{"log_loss"} (binary), or \code{"harrell_c"}, \code{"time_auc"}, \code{"uno_c"},
#' \code{"ipcw_brier"} (survival). The last three are saved fixed-horizon
#' cumulative/dynamic AUC, truncated Uno concordance, and IPCW Brier score.
#' Leave this argument \code{NULL} for prediction and run tables. Nested reference,
#' skill, calibration and decision-curve objects are not flattened here.
#' }
#' \item{roles}{
#' Optional named dataset roles: \code{"development_cv"}, \code{"training"},
#' \code{"validation"}, \code{"test"} or \code{"external"}. CV objects retain their CV role.
#' Unlabeled external evaluations remain external, not assumed validation.
#' }
#' \item{path}{
#' New HTML file path. Defaults to a temporary report.
#' }
#' \item{open}{
#' Open the page in the configured browser; defaults to \code{interactive()}.
#' }
#' \item{overwrite}{
#' Replace an existing report only when explicitly TRUE.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_render_report()}}{
#' Create an offline HTML report without fitting models or recomputing metrics.
#' 
#' \strong{Returns:} Invisibly, the normalized absolute output path. Requires optional
#' \code{ggplot2} and SVG support in R, but no Pandoc or network access.
#' 
#' \strong{Details:} The report contains saved scalar metric and status tables,
#' performance plots, calibration and decision curves. CV and external
#' sections remain separate. Pipeline order is preserved; no winner,
#' confidence interval, missing metric or missing curve is inferred.
#' Survival sections with a saved horizon or fixed-horizon metric also show
#' \code{time_auc}, \code{uno_c}, and \code{ipcw_brier}, retaining their horizon and unit.
#' Horizon-free sections do not request these metrics. In mixed sections,
#' unavailable values remain unavailable and different targets are labeled
#' separately, never pooled or converted.
#' Unsupported or malformed inputs are errors, not silently omitted sections.
#' Default tables and graphs show usable results only. Failure/issue counts
#' and a companion CSV link remain visible; all-failed sections explicitly
#' state that no usable results are available. The CSV is header-only when
#' there are no issues. See \code{\link[=sb_failure_table]{sb_failure_table()}} for its record definitions.
#' 
#' Patient-level rows, raw result objects and free-text error messages are
#' never embedded. Caller-supplied titles, cohort labels and pipeline labels
#' remain visible and must be reviewed before sharing. This is not a guarantee
#' of anonymization: plotted calibration curves describe the supplied cohort.
#' Saved numeric summaries are formatted to eight significant digits.
#' The companion CSV includes saved error text, which may be sensitive;
#' review it separately before sharing. Patient-level rows are not exported.
#' If an I/O error interrupts the two-file write, a completed CSV may remain;
#' it is not deleted or overwritten automatically. Inspect the reported paths.
#' }
#' 
#' \subsection{\code{sb_render_clinician_report()}}{
#' Render grouped clinical views from saved results
#' 
#' \strong{Returns:} Invisibly the absolute report filename. Only unavailable or partial
#' views and errors appear in the status table; normal status rows are hidden.
#' Model settings appear in the report body. Original per-context diagnostics remain
#' in x. No model is refitted and no extra metrics are requested.
#' }
#' 
#' \subsection{\code{sb_render_model_report()}}{
#' Render an offline interactive explanation of a saved linear model
#' 
#' \strong{Returns:} Invisibly, the absolute HTML path. This single offline file embeds
#' R-rendered SVG plots and selectable coefficient/contribution tables.
#' 
#' \strong{Details:} Reuses \code{\link[=sb_additive_shap]{sb_additive_shap()}}, \code{\link[=sb_plot_shap_beeswarm]{sb_plot_shap_beeswarm()}},
#' \code{\link[=sb_plot_shap_importance]{sb_plot_shap_importance()}} and \code{\link[=sb_plot_shap_waterfall]{sb_plot_shap_waterfall()}}; it does not
#' introduce a second SHAP algorithm. One sample selector switches both the
#' waterfall and survival view. SHAP is on log-odds/log-hazard score scales,
#' never probability effects. The baseline is zero-score, not average-patient
#' survival. No external service or resource is contacted. Patient IDs are
#' replaced by row numbers in the report, but plotted values and predictions
#' remain sensitive; this is not anonymization. Rendering all patient views
#' can take time and produce a large file. Use the separate plot functions
#' when an interactive report is not required.
#' Optional svglite produces compact text-based SVG; otherwise the R SVG
#' device is used. Both preserve the same plotted data and require no network.
#' }
#' 
#' \subsection{\code{sb_render_visual_report()}}{
#' Export a gallery of saved plots to offline HTML
#' 
#' \strong{Returns:} Invisibly, the normalized absolute output path. Requires optional
#' ggplot2 and an R build with SVG support; no Pandoc or internet is needed.
#' 
#' \strong{Details:} Embeds rendered SVG images, not R objects or patient-level tables.
#' Plot titles, labels and visible points remain present, including any
#' individual values in SHAP plots. This is not an anonymization service:
#' review all supplied plots before sharing. No metrics are recomputed and
#' no model is fitted. Rendered images are selectable only as images; inspect
#' the original ggplot data locally for numeric values. Device state is
#' restored even if drawing fails. A failed output write may leave a partial
#' file, which is not automatically removed or overwritten.
#' Grouped sections use native HTML disclosure controls and work without
#' JavaScript. Every supplied image remains embedded, including collapsed
#' details; collapsing is not privacy protection or file-size reduction.
#' Expand the sections you want before printing. Extra legends are selectable
#' text in the report body. Supply \code{compact = FALSE} to keep original labels.
#' }
#' 
#' \subsection{\code{sb_export_run()}}{
#' Export a completed workflow into a new directory
#' 
#' \strong{Returns:} Invisible absolute output directory. Tables use CSV; figures
#' use PDF by default or PNG when requested.
#' Cohort/plot indexes map numbered files to their original names. Patient
#' predictions are identifiable data; store this directory appropriately.
#' Original model predictions and evaluated predictions are separate files.
#' Cohort processing records and original metrics accompany evaluations.
#' Export uses existing result and plot objects; it does not refit the model.
#' }
#' 
#' \subsection{\code{sb_run_figure_keys()}}{
#' List the exact keys of saved plots in a completed workflow, without
#' drawing or writing them. The keys refer to whole plots before pagination.
#' }
#' 
#' \subsection{\code{sb_export_failures()}}{
#' Export saved failures to a separate CSV
#' 
#' \strong{Returns:} Invisibly, the normalized absolute CSV path. A header-only CSV is
#' written when there are no issues. Uses base R and does not require ggplot2.
#' 
#' \strong{Details:} The CSV preserves saved error text, including quoted commas and
#' line breaks, plus repetition, held-out cohort and model-role columns for
#' composite objects. It excludes patient-level prediction rows but raw error text
#' still requires privacy review. Results are not retrained or recomputed.
#' }
#' 
#' \subsection{\code{sb_result_table()}}{
#' Read saved metrics, patient predictions, or fold-run records without fitting,
#' evaluating, ranking, or modifying the supplied result.
#' 
#' \strong{Returns:} A data frame in supplied pipeline order. Metric tables contain
#' \code{pipeline_id}, \code{task_type}, \code{metric}, \code{estimate}, \code{status}, \code{run_status},
#' \code{scope}, and \code{aggregation}. \code{status} belongs to the selected metric;
#' \code{run_status} belongs to the complete CV run or external evaluation.
#' External metric tables also retain evaluation errors, metric errors and
#' patient-match status. Prediction and run tables
#' preserve saved columns and values and append \code{task_type}. Rows from each
#' pipeline remain in saved order; output row names are reset.
#' For external objects, predictions are saved scored rows including outcomes;
#' runs contain one evaluation summary per pipeline, not fictitious CV folds.
#' The three fixed-horizon metrics additionally retain \code{prediction_horizon},
#' \code{time_unit}, and saved metric errors. Other metric tables keep their
#' existing columns. Different horizons or units are distinct targets; no
#' conversion or pooling is performed.
#' 
#' \strong{Details:} Composite objects add \code{repeat_id}, \code{heldout_cohort_id}, \code{model_role},
#' \code{validation_scheme}, and \code{view_id} (a structural entry key). These identify
#' saved repetitions, cohorts, and member/stack outputs without changing
#' \code{pipeline_id}. Unused context fields remain missing. Rows stay in stored
#' order; no averaging, pooling, or confidence interval is inferred. IECV is
#' method development, not final independent external validation. Standalone
#' stacking includes outer-CV members and the stack, not inner training folds.
#' When a stack descriptor omits its horizon or unit, an unambiguous value
#' saved in its prediction rows is retained. Conflicting or absent row
#' metadata is not guessed; metric-level metadata still takes precedence.
#' Failed methods and missing or infinite estimates remain in the
#' table. A metric can be evaluable even when the complete run failed at
#' another requested output. No estimate is filled, clipped, recalibrated, or
#' recomputed. Binary ROC AUC, average precision, and survival Harrell C keep
#' their saved fold-aware aggregation; they are not replaced by a pooled
#' calculation. Confidence intervals and a winning model are not inferred.
#' External metrics retain their saved external-cohort scope and aggregation.
#' A failed external evaluation with no saved supported metric has an \code{NA}
#' estimate and metric status \code{"not_saved"}; its failure information remains.
#' An unsupported metric or a missing metric in a nonfailed result is an error.
#' For the three fixed-horizon metrics, an absent value in a successful
#' horizon-free result is instead \code{"not_requested"}; an absent requested
#' value or an absent value in a failed result is \code{"not_saved"}. Its estimate
#' remains \code{NA}. A horizon stored in the metric takes precedence over the
#' task/endpoint descriptor; unavailable metadata remains missing.
#' Failed evaluations with no scored rows contribute no prediction rows;
#' inspect the run table for their failures. Cohort labels are descriptive,
#' not evidence of independence. No cohorts are pooled by this function.
#' 
#' Run records retain their existing error information. An error stage is not
#' invented when the source did not save one. Prediction tables can contain
#' patient identifiers and outcomes; review them before sharing or exporting.
#' }
#' 
#' \subsection{\code{sb_failure_table()}}{
#' Read failure records without executing, modifying or dropping model results.
#' 
#' \strong{Returns:} A data frame with one row per saved pipeline/run/metric issue.
#' Pipeline summaries and detailed fold or metric issues are distinct rows;
#' the number of rows is not the number of failed pipelines. Includes cohort,
#' validation boundary, pipeline, issue type, structural stage, fold, metric,
#' source status, complete-run status, source error stage/class/message, scalar
#' estimate and saved horizon/unit. No patient prediction rows are included.
#' 
#' \strong{Details:} Explicitly unrequested/not-applicable components are not failures.
#' Composite results add repetition, held-out cohort, member/stack role,
#' validation scheme and structural \code{view_id} columns. IECV is development
#' validation. Mixed named lists keep these columns with missing values for
#' ordinary results. Saved failures are neither averaged nor deduplicated.
#' A saved success status with a nonfinite scalar estimate is labeled through
#' issue_type = "metric_nonfinite" without changing its source status.
#' Missing source stages and reasons remain NA. Structural stage "cv_fold"
#' identifies a fold record, not an inferred fitting or prediction error stage.
#' Empty results return the same columns with zero rows. Raw error messages
#' can contain sensitive information: review before sharing or opening in
#' software that may interpret cell contents as formulas.
#' }
#' 
#' \subsection{\code{sb_view_results()}}{
#' Create a standalone HTML page from saved results. Rows show feature selectors
#' and learner names; columns show development, validation and test cohorts.
#' Metrics are read from the supplied results without refitting models.
#' 
#' \strong{Returns:} Invisibly, the normalized path of the saved HTML file.
#' 
#' \strong{Details:} The page contains summary metrics and model metadata, not patient
#' rows or raw error messages. Original model IDs and parameter settings
#' remain available in details. Default sorting uses development CV when
#' present; it never selects or refits a winning model. Optional validation
#' means exclude training/test cohorts and require all selected validation
#' cells to be available for the same metric and horizon.
#' }
#' }
#' @rawRd \seealso{
#' \code{\link[=sb_result_table]{sb_result_table()}}, \code{\link[=sb_plot_performance]{sb_plot_performance()}}, \code{\link[=sb_plot_calibration]{sb_plot_calibration()}},
#' \code{\link[=sb_plot_dca]{sb_plot_dca()}}
#' 
#' \code{\link[=sb_render_report]{sb_render_report()}}, \code{\link[=sb_plot_roc]{sb_plot_roc()}}, \code{\link[=sb_plot_meta_forest]{sb_plot_meta_forest()}},
#' \code{\link[=sb_plot_selection_upset]{sb_plot_selection_upset()}}, \code{\link[=sb_plot_shap_waterfall]{sb_plot_shap_waterfall()}}
#' 
#' \code{\link[=sb_plot_performance]{sb_plot_performance()}}
#' 
#' \code{\link[=sb_export_failures]{sb_export_failures()}}
#' }
#' @rawRd \examples{
#' # sb_export_run
#' development <- bi_example(40, seed = 42)
#' run <- bi_run(development, outcome = "outcome", id = "patient_id",
#'   predictors = c("gene_a", "gene_b"), models = "Logistic", folds = 2,
#'   seed = 42, plots = FALSE)
#' finished <- finish_run(run, "Logistic", plots = FALSE)
#' destination <- tempfile("stormbreaker-export-")
#' export_run(finished, destination)
#' list.files(destination)
#' 
#' # sb_result_table
#' d <- sb_example_binary(n = 60, seed = 1)
#' task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
#' result <- sb_run_cv(
#'   task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1
#' )
#' sb_result_table(result, metric = "roc_auc")
#' head(sb_result_table(result, "predictions"))
#' sb_result_table(result, "runs")
#' 
#' # sb_view_results
#' d <- sb_example_binary(50, seed = 1)
#' task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#' result <- sb_compare_models(task, models = "logistic", folds = 3, view = FALSE)
#' page <- sb_view_results(result, open = FALSE)
#' file.exists(page)
#' }
NULL

# ---- report-clinical ----

# ---- render-clinician-report ----

#' @rdname sb_render_report
#' @export
sb_render_clinician_report <- function(x, file, title = "Clinical model evaluation", complete = FALSE) {
  if (!inherits(x, "stormbreaker_clinician_views") || !length(x$plots)) sb_abort(
    "Supply clinician views with at least one available plot; inspect status for unavailable results.",
    "stormbreaker_bad_report_input")
  if (!is.logical(complete) || length(complete) != 1L || is.na(complete)) sb_abort(
    "complete must be TRUE or FALSE.", "stormbreaker_bad_report_input")
  issues <- x$status[x$status$status != "available", , drop = FALSE]
  captions <- x$captions
  if (!is.null(x$models) && nrow(x$models)) {
    settings <- vapply(seq_len(nrow(x$models)), function(i) {
      m <- x$models[i, ]
      parts <- c(m$model, if (m$selector != "None") m$selector,
        m$parameters, m$selection_parameters)
      paste(parts[nzchar(parts)], collapse = "; ")
    }, character(1))
    captions[1L] <- paste(c(captions[1L], "Model settings:", settings), collapse = "\n")
  }
  sb_render_visual_report(x$plots, file, title = title, groups = x$groups,
    expanded = if (complete) unique(x$groups) else NULL, captions = captions,
    availability = if (nrow(issues)) issues else NULL, width = x$width, height = x$height)
}

# ---- render-visual-report ----

#' @rdname sb_render_report
#' @export
sb_render_visual_report <- function(plots, file,
                                    title = "Saved-result visualization gallery",
                                    width = 11, height = 7, groups = NULL,
                                    expanded = NULL, compact = TRUE, captions = NULL,
                                    availability = NULL) {
  text_ok <- function(x) {
    is.character(x) && length(x) == 1L &&
      !is.na(x) && nzchar(trimws(x))
  }
  if (!is.list(plots) || !length(plots) || is.null(names(plots)) ||
    anyNA(names(plots)) || anyDuplicated(names(plots)) ||
    any(!nzchar(trimws(names(plots)))) ||
    !all(vapply(plots, inherits, logical(1), what = "ggplot"))) {
    sb_abort(
      "Supply a non-empty uniquely named list of ggplot objects.",
      "stormbreaker_bad_report_input"
    )
  }
  if (!text_ok(file) || !grepl("\\.html$", file, ignore.case = TRUE) ||
    !dir.exists(dirname(file)) || !text_ok(title)) {
    sb_abort(
      "Supply a title and a new .html path in an existing directory.",
      "stormbreaker_bad_report_input"
    )
  }
  dimensions_ok <- function(x) {
    is.numeric(x) && length(x) %in% c(1L, length(plots)) &&
      all(is.finite(x) & x > 0)
  }
  if (!dimensions_ok(width) || !dimensions_ok(height)) {
    sb_abort(
      "Plot dimensions must be positive finite numbers, length one or number of plots.",
      "stormbreaker_bad_report_input"
    )
  }
  if (!is.logical(compact) || length(compact) != 1L || is.na(compact)) sb_abort(
    "compact must be TRUE or FALSE.", "stormbreaker_bad_report_input")
  vector_ok <- function(x) is.character(x) && length(x) == length(plots) && !anyNA(x)
  if (!is.null(groups) && (!vector_ok(groups) || any(!nzchar(trimws(groups))))) sb_abort(
    "Supply one non-empty group name per plot.", "stormbreaker_bad_report_input")
  if (!is.null(captions) && !vector_ok(captions)) sb_abort(
    "Supply one caption per plot; use an empty string when none is needed.", "stormbreaker_bad_report_input")
  topics <- unique(groups)
  if (!is.null(availability) && (!is.data.frame(availability) || !ncol(availability) ||
      !all(vapply(availability, is.atomic, logical(1))))) sb_abort(
    "Availability must be a data frame with atomic columns.", "stormbreaker_bad_report_input")
  if (!is.null(expanded) && (!is.character(expanded) || anyNA(expanded) ||
      anyDuplicated(expanded) || !all(expanded %in% topics))) sb_abort(
    "Expanded topics must be unique supplied group names.", "stormbreaker_bad_report_input")
  if (is.null(expanded)) expanded <- utils::head(topics, 1L)
  refuse_existing <- function() {
    if (file.exists(file)) {
      sb_abort(
        "The report path already exists; choose a new filename.",
        "stormbreaker_report_exists"
      )
    }
  }
  refuse_existing()
  sb_require_curve_plotting()
  width <- rep(width, length.out = length(plots))
  height <- rep(height, length.out = length(plots))
  notes <- vector("list", length(plots))
  draw <- function(i) {
    path <- tempfile(fileext = ".svg")
    grDevices::svg(path, width = width[i], height = height[i])
    device <- grDevices::dev.cur()
    on.exit(
      {
        if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
      },
      add = TRUE
    )
    plot <- plots[[i]]
    legend <- captions[i]
    if (compact) {
      original_labels <- plot$labels[c("title", "subtitle", "caption")]
      label_text <- vapply(original_labels, function(z) paste(as.character(z), collapse = " "), character(1))
      legend <- unique(c(label_text[label_text != names(plots)[i]], legend))
      plot <- plot + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
    }
    legend <- legend[!is.na(legend) & nzchar(legend)]
    notes[[i]] <<- legend
    print(plot)
    grDevices::dev.off(device)
    svg <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    label <- sb_report_escape(names(plots)[i])
    paste0(
      "<section id='plot-", i, "'><h2>", label,
      "</h2><figure><img alt='", label, "' src='data:image/svg+xml;charset=utf-8,",
      utils::URLencode(svg, reserved = TRUE), "'></figure></section>"
    )
  }
  sections <- vapply(seq_along(plots), draw, character(1))
  nav <- paste0("<li><a href='#plot-", seq_along(plots), "'>",
    sb_report_escape(names(plots)), "</a></li>",
    collapse = ""
  )
  if (length(topics)) {
    sections <- vapply(seq_along(topics), function(i) {
      indices <- which(groups == topics[i])
      paste0("<details id='topic-", i, "'", if (topics[i] %in% expanded) " open" else "",
        "><summary>", sb_report_escape(topics[i]), " <span>(", length(indices), ")</span></summary>",
        paste(sections[indices], collapse = "\n"), "</details>")
    }, character(1))
    nav <- paste0("<li><a href='#topic-", seq_along(topics), "'>",
      sb_report_escape(topics), "</a></li>", collapse = "")
  }
  html <- paste0(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'><title>",
    sb_report_escape(title), "</title><style>",
    "*{box-sizing:border-box}body{margin:0;color:#172b3a;background:#f5f7f9;font:16px/1.55 system-ui,sans-serif}",
    "main{max-width:1280px;margin:auto;padding:28px}header{border-bottom:4px solid #24547a}h1{line-height:1.2}",
    "section{background:white;padding:20px;margin:24px 0}h1,h2,a{overflow-wrap:anywhere}",
    "details{background:white;margin:24px 0;border:1px solid #ccd6df}summary{cursor:pointer;padding:18px;font-size:1.2rem;font-weight:600;overflow-wrap:anywhere}",
    "summary span{font-weight:400;color:#465966}summary:focus-visible{outline:3px solid #bf6b00}details section{margin:0;border-top:1px solid #e2e8ed}",
    "figure{margin:0}img{display:block;width:100%;height:auto}figcaption,footer{font-size:.9rem;color:#465966}figcaption{white-space:pre-line;overflow-wrap:anywhere}",
    "table{border-collapse:collapse;width:100%;font-size:.9rem}td,th{padding:8px;text-align:left;border-bottom:1px solid #ccd6df;overflow-wrap:anywhere}",
    "a{color:#24547a}a:focus-visible{outline:3px solid #bf6b00}footer{border-top:1px solid #ccd6df;padding:18px 0}",
    "@media(max-width:600px){main{padding:12px}section{padding:8px}h1{font-size:1.6rem}}",
    "@media print{body{background:white}main{padding:0;max-width:none}nav{display:none}section{padding:0}figure{break-inside:avoid}}",
    "</style></head><body><main><header><h1>", sb_report_escape(title),
    "</h1></header>",
    "<nav aria-label='Gallery sections'><ol>", nav, "</ol></nav>",
    if (!is.null(availability)) paste0("<details><summary>View availability</summary>", sb_model_report_table(availability), "</details>") else "",
    paste(sections, collapse = "\n"),
    if (any(lengths(notes))) paste0("<section id='report-notes'><h2>Methods and interpretation</h2>",
      paste(vapply(which(lengths(notes) > 0L), function(i) paste0("<h3>",
        sb_report_escape(names(plots)[i]), "</h3><p>",
        paste(sb_report_escape(unique(notes[[i]])), collapse = "<br>"), "</p>"), character(1)), collapse = "\n"),
      "</section>") else "",
    "</main></body></html>"
  )
  refuse_existing()
  writeLines(enc2utf8(html), file, useBytes = TRUE)
  invisible(normalizePath(file, winslash = "/", mustWork = TRUE))
}

# ---- report-export ----

# ---- render-model-report ----

sb_model_report_image <- function(plot, label, width = 11, height = 7) {
  path <- tempfile(fileext = ".svg")
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(path, width = width, height = height)
  } else {
    grDevices::svg(path, width = width, height = height)
  }
  device <- grDevices::dev.cur()
  on.exit(
    {
      if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
    },
    add = TRUE
  )
  print(plot)
  grDevices::dev.off(device)
  svg <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  paste0(
    "<img alt='", sb_report_escape(label),
    "' src='data:image/svg+xml;charset=utf-8,", utils::URLencode(svg, reserved = TRUE), "'>"
  )
}

sb_model_report_table <- function(data) {
  cells <- lapply(data, function(column) {
    if (is.numeric(column)) column <- format(signif(column, 6), trim = TRUE)
    sb_report_escape(as.character(column))
  })
  rows <- vapply(seq_len(nrow(data)), function(i) {
    paste0(
      "<tr>",
      paste0("<td>", vapply(cells, `[`, character(1), i), "</td>", collapse = ""), "</tr>"
    )
  }, character(1))
  paste0("<table><thead><tr>", paste0("<th scope='col'>", sb_report_escape(names(data)),
    "</th>",
    collapse = ""
  ), "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table>")
}

#' @rdname sb_render_report
#' @export
sb_render_model_report <- function(fit, new_data, id_col, background_data, file,
                                   development = NULL,
                                   title = "Saved model: parameters and explanations",
                                   top_n = 20L, selection_results = NULL, compact = TRUE) {
  text_ok <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
  if (!text_ok(file) || !grepl("\\.html$", file, ignore.case = TRUE) ||
    !dir.exists(dirname(file)) || !text_ok(title)) {
    sb_abort(
      "Supply a title and a new .html path in an existing directory.", "stormbreaker_bad_report_input"
    )
  }
  if (file.exists(file)) {
    sb_abort(
      "The report path already exists; choose a new filename.",
      "stormbreaker_report_exists"
    )
  }
  top_n <- sb_viz_positive_integer(top_n, "top_n")
  if (!is.logical(compact) || length(compact) != 1L || is.na(compact)) {
    sb_abort("compact must be TRUE or FALSE.", "stormbreaker_bad_report_input")
  }
  fit <- sb_validate_final_fit_for_prediction(fit)
  coefficients <- sb_model_coefficients(fit, background_data)
  shap <- sb_additive_shap(fit, new_data, id_col, background_data)
  curves <- NULL
  if (fit$task_type == "survival") {
    if (is.null(fit$survival_curve)) {
      if (is.null(development)) {
        sb_abort(
          "Supply the original development task or an augmented survival model.",
          "stormbreaker_missing_survival_curve"
        )
      }
      fit <- sb_add_survival_curve(fit, development)
    }
    curves <- sb_survival_curves(fit, new_data, id_col)
  }
  image <- function(plot, label, width = 11, height = 7) {
    if (compact) plot <- sb_compact_report_plot(plot)
    sb_model_report_image(plot, label, width, height)
  }
  genes <- if (is.null(selection_results)) "" else sb_model_report_gene_section(selection_results, image)
  shap_height <- max(4.2, 1.7 + 0.30 * min(top_n, ncol(shap$contributions)))
  waterfall_height <- max(4.6, shap_height + 0.6)
  raw <- image(sb_plot_coefficients(coefficients), "Raw signed coefficients")
  scaled <- image(sb_plot_coefficients(coefficients, scale = "reference_sd"), "Coefficients per training SD")
  swarm <- image(sb_plot_shap_beeswarm(shap, top_n, new_data, id_col), "Global SHAP beeswarm", height = shap_height)
  importance <- image(sb_plot_shap_importance(shap, top_n), "Global mean absolute SHAP importance", height = shap_height)
  base <- if (!is.null(curves)) {
    paste0(
      "<details><summary>Baseline survival alone</summary>",
      image(sb_plot_survival_curves(curves, rows = integer()), "Original-training zero-score baseline"),
      "</details>"
    )
  } else {
    ""
  }
  n <- nrow(shap$contributions)
  samples <- vapply(seq_len(n), function(i) {
    summary <- paste0(
      "<p><strong>Observation ", i, "</strong> | Model score: ",
      format(signif(shap$scores$score[i], 6)), " | Background score: ",
      format(signif(shap$reference_value, 6)), "</p>"
    )
    if (is.null(curves)) {
      summary <- paste0(
        summary, "<p>Class-1 probability: ",
        sprintf("%.2f%%", 100 * stats::plogis(shap$scores$score[i])), "</p>"
      )
    } else {
      probability <- curves$survival[curves$baseline$time == curves$prediction_horizon, i]
      summary <- paste0(
        summary, "<p>Survival at ", curves$prediction_horizon, " ",
        sb_report_escape(curves$time_unit), ": <strong>", sprintf("%.2f%%", probability * 100), "</strong></p>"
      )
    }
    contributions <- data.frame(
      Feature = colnames(shap$contributions),
      Contribution = unname(shap$contributions[i, ])
    )
    paste0(
      summary,
      image(sb_plot_shap_waterfall(shap, i, top_n), paste("SHAP waterfall for observation", i), height = waterfall_height),
      "<details><summary>SHAP values</summary>",
      sb_model_report_table(contributions), "</details>",
      if (!is.null(curves)) {
        image(
          sb_plot_survival_curves(curves, rows = i),
          paste("Predicted survival for observation", i)
        )
      } else {
        ""
      }
    )
  }, character(1))
  options <- paste0("<option value='", seq_len(n), "'>Observation ", seq_len(n), "</option>", collapse = "")
  templates <- paste0("<template id='sample-", seq_len(n), "'>", samples, "</template>", collapse = "\n")
  html <- paste0(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'><title>", sb_report_escape(title),
    "</title><style>*{box-sizing:border-box}body{margin:0;background:#f5f7f9;color:#172b3a;font:16px/1.6 system-ui,sans-serif}",
    "main{max-width:1240px;margin:auto;padding:28px}header{border-bottom:4px solid #24547a}h1{line-height:1.2}",
    "h1,h2,p,a,td,th{overflow-wrap:anywhere}section{background:white;margin:26px 0;padding:22px}",
    "img{display:block;width:100%;height:auto}nav{display:flex;gap:20px;flex-wrap:wrap;margin:20px 0}",
    "a{color:#24547a}select,button{font:inherit;padding:7px 12px;max-width:100%;background:white;border:1px solid #8496a2;border-radius:4px}",
    "label{display:inline-block;margin-right:10px;font-weight:600}select:focus-visible,summary:focus-visible,a:focus-visible{outline:3px solid #bf6b00;outline-offset:3px}",
    "details{margin:16px 0}summary{cursor:pointer;color:#24547a}table{border-collapse:collapse;width:100%;font-size:.9rem;table-layout:fixed}",
    "th,td{padding:8px;border-bottom:1px solid #dce3e8;text-align:left}th{background:#edf2f5}",
    "[hidden]{display:none!important}.note,footer{color:#465966;font-size:.92rem}footer{border-top:1px solid #ccd6df;padding:20px 0}",
    "@media(max-width:600px){main{padding:12px}section{padding:10px}h1{font-size:1.6rem}th,td{padding:5px;font-size:.8rem}}",
    "@media print{body{background:white}main{max-width:none;padding:0}.controls,nav{display:none}section{break-inside:avoid;padding:0}details>table{display:table}}",
    "</style></head><body><main><header><h1>", sb_report_escape(title), "</h1>",
    "<p>", n, " samples | ", nrow(coefficients$table), " model features | ",
    sb_report_escape(sb_display_learner(fit$pipeline$learner, fit$pipeline$learner_params)),
    "</p><p>SHAP: ", sb_report_escape(sb_report_score_label(shap$score_scale)),
    " | Training reference: ", shap$background$n, "</p></header>",
    "<nav aria-label='Report sections'>", if (nzchar(genes)) "<a href='#genes'>Gene selection</a>",
    "<a href='#parameters'>Coefficients</a><a href='#global'>Global SHAP</a><a href='#individual'>Individual prediction</a></nav>",
    genes,
    "<section id='parameters'><h2>Model coefficients</h2>",
    if (!is.null(coefficients$intercept)) paste0("<p>Intercept: ", format(signif(coefficients$intercept, 6)), "</p>"),
    "<div class='controls'>",
    "<label for='coefficient-view'>Coefficient units</label><select id='coefficient-view'><option value='raw'>Original input units</option><option value='sd'>Per training standard deviation</option></select></div>",
    "<div id='coef-raw'>", raw, "</div><div id='coef-sd' hidden>", scaled, "</div>",
    "<details><summary>All coefficients</summary>", sb_model_report_table(coefficients$table),
    "</details></section>",
    "<section id='global'><h2>Global SHAP</h2><div class='controls'>",
    "<label for='shap-view'>Global view</label><select id='shap-view'><option value='swarm'>Beeswarm</option><option value='importance'>Mean absolute contribution</option></select></div>",
    "<div id='shap-swarm'>", swarm, "</div><div id='shap-importance' hidden>", importance,
    "</div></section>",
    "<section id='individual'><h2>Individual prediction</h2><div class='controls'><label for='sample-select'>Observation</label>",
    "<select id='sample-select'>", options, "</select></div>",
    "<p id='sample-status' role='status' aria-live='polite'>Showing observation 1 of ", n, "</p><div id='sample-panel'>",
    samples[1], "</div>", base, "</section>", templates,
    "<noscript><p>Enable JavaScript to switch views.</p></noscript></main>",
    "<script>'use strict';const el=id=>document.getElementById(id);",
    "el('coefficient-view').addEventListener('change',e=>{el('coef-raw').hidden=e.target.value!=='raw';el('coef-sd').hidden=e.target.value!=='sd';});",
    "el('shap-view').addEventListener('change',e=>{el('shap-swarm').hidden=e.target.value!=='swarm';el('shap-importance').hidden=e.target.value!=='importance';});",
    if (nzchar(genes)) "el('gene-view').addEventListener('change',e=>{el('gene-selector').hidden=e.target.value!=='selector';el('gene-glmnet_active_set').hidden=e.target.value!=='glmnet_active_set';});",
    "el('sample-select').addEventListener('change',e=>{const t=el('sample-'+e.target.value);el('sample-panel').replaceChildren(t.content.cloneNode(true));",
    "el('sample-status').textContent='Showing observation '+e.target.value+' of ", n, "';});</script></body></html>"
  )
  if (file.exists(file)) {
    sb_abort(
      "The report path already exists; choose a new filename.",
      "stormbreaker_report_exists"
    )
  }
  writeLines(enc2utf8(html), file, useBytes = TRUE)
  invisible(normalizePath(file, winslash = "/", mustWork = TRUE))
}

sb_report_score_label <- function(x) {
  gsub("class_1_log_odds", "log odds (class 1)",
    gsub("zero_reference_log_hazard_ratio", "log hazard ratio", x, fixed = TRUE),
    fixed = TRUE
  )
}

sb_compact_report_plot <- function(plot) {
  xlabel <- plot$labels$x
  plot <- plot + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
  if (is.character(xlabel)) plot <- plot + ggplot2::labs(x = sb_report_score_label(xlabel))
  plot
}

sb_model_report_gene_section <- function(results, image) {
  views <- vapply(c("selector", "glmnet_active_set"), function(source) {
    ranking <- sb_gene_selection_ranking(results, source = source)
    tab <- ranking$ranking
    display <- data.frame(
      Rank = tab$rank, Gene = tab$feature,
      Selected = tab$selected_count, Eligible = tab$eligible_count,
      Frequency = ifelse(is.finite(tab$selection_frequency), sprintf("%.1f%%", 100 * tab$selection_frequency), "NA")
    )
    paste0(
      "<div id='gene-", source, "'", if (source != "selector") " hidden", ">",
      if (ranking$status == "incomplete") "<p><strong>Incomplete selection records</strong></p>",
      if (ranking$status == "not_available") {
        "<p>No saved learned-selection counts.</p>"
      } else {
        image(sb_plot_gene_selection_ranking(ranking), paste("Gene selection ranking", source), height = 7)
      },
      "<details><summary>All genes (", nrow(tab), ")</summary>",
      sb_model_report_table(display), "</details></div>"
    )
  }, character(1))
  paste0(
    "<section id='genes'><h2>Gene selection ranking</h2><div class='controls'>",
    "<label for='gene-view'>Count</label><select id='gene-view'><option value='selector'>Selector selections</option>",
    "<option value='glmnet_active_set'>LASSO / Elastic Net nonzero</option></select></div>",
    "<p>Successful CV folds; shared selector copies count once.</p>", paste(views, collapse = ""), "</section>"
  )
}

# ---- render-report ----

sb_report_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  for (pair in list(
    c("&", "&amp;"), c("<", "&lt;"), c(">", "&gt;"),
    c('"', "&quot;"), c("'", "&#39;")
  )) {
    x <- gsub(pair[1], pair[2], x, fixed = TRUE)
  }
  x
}

sb_report_table <- function(tab, columns, caption) {
  columns <- unique(c(intersect(c(
    "repeat_id", "heldout_cohort_id", "model_role",
    "validation_scheme", "view_id"
  ), names(tab)), columns))
  tab <- tab[intersect(columns, names(tab))]
  header <- paste0("<th scope='col'>", sb_report_escape(names(tab)), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(tab)), function(i) {
    cells <- vapply(tab, function(column) {
      value <- column[i]
      if (is.numeric(value)) value <- format(value, digits = 8, trim = TRUE)
      sb_report_escape(value)
    }, character(1))
    paste0("<tr>", paste0("<td>", cells, "</td>", collapse = ""), "</tr>")
  }, character(1))
  paste0(
    "<div class='table-scroll' tabindex='0' role='region' aria-label='",
    sb_report_escape(caption), "'><table><caption>", sb_report_escape(caption),
    "</caption><thead><tr>", header, "</tr></thead><tbody>",
    paste(rows, collapse = ""), "</tbody></table></div>"
  )
}

sb_report_plot <- function(plot, label, panels) {
  plot <- sb_compact_report_plot(plot)
  path <- tempfile(fileext = ".svg")
  grDevices::svg(path, width = 10, height = max(5, 3 * panels))
  device <- grDevices::dev.cur()
  on.exit(
    {
      if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
    },
    add = TRUE
  )
  print(plot)
  grDevices::dev.off(device)
  svg <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  paste0(
    "<figure><img alt='", sb_report_escape(label),
    "' src='data:image/svg+xml;charset=utf-8,", utils::URLencode(svg, reserved = TRUE),
    "'><figcaption>", sb_report_escape(label), "</figcaption></figure>"
  )
}

sb_report_section <- function(x, label, index) {
  external <- inherits(x, c("stormbreaker_external_evaluation", "stormbreaker_external_comparison"))
  results <- if (inherits(x, "stormbreaker_cv_comparison")) {
    x$results
  } else if (
    inherits(x, "stormbreaker_external_comparison")) {
    x$evaluations
  } else {
    list(x)
  }
  if (!length(results)) sb_abort("Report sections cannot be empty.", "stormbreaker_bad_report_input")
  runs <- sb_result_table(x, "runs", cohort_label = label)
  types <- unique(runs$task_type)
  if (length(types) != 1L || !types %in% c("binary", "survival")) {
    sb_abort("Each section must contain one supported task type.", "stormbreaker_bad_report_input")
  }
  metrics <- if (types == "binary") c("roc_auc", "average_precision", "brier_score", "log_loss") else "harrell_c"
  if (types == "survival") {
    requested <- any(vapply(results, function(result) {
      descriptor <- if (external) result$endpoint else result$task
      !is.null(descriptor$prediction_horizon) ||
        any(sb_horizon_view_metrics() %in% names(result$metrics))
    }, logical(1)))
    if (requested) metrics <- c(metrics, sb_horizon_view_metrics())
  }
  boundary <- if (external) "External evaluation" else "Development CV"
  context <- attr(x, "stormbreaker_view_context")
  if (!is.null(context)) boundary <- sb_view_boundary(context)
  successful <- vapply(results, function(result) identical(result$status, "success"), logical(1))
  successful_ids <- vapply(results[successful], function(result) {
    if (external) result$pipeline_id else result$pipeline$id
  }, character(1))
  runs <- runs[runs$pipeline_id %in% successful_ids & runs$status == "success", , drop = FALSE]
  body <- c(
    paste0(
      "<section id='section-", index, "'><h2>", sb_report_escape(label),
      "</h2><p class='boundary'>", boundary, " &middot; ", types, "</p>"
    ),
    paste0(
      "<p>", sum(successful), " of ", length(results),
      " pipelines completed", if (any(!successful)) paste0("; ", sum(!successful), " failed/incomplete"), ".</p>"
    ),
    if (!any(successful)) "<p><strong>No usable results in this section. This is not a successful analysis.</strong></p>",
    sb_report_table(runs, c(
      "pipeline_id", "task_type", "status", "fold",
      "n", "patient_match_status"
    ), "Successful run status")
  )
  for (metric in metrics) {
    tab <- sb_result_table(x, metric = metric, cohort_label = label)
    tab <- tab[sb_performance_display_rows(tab), , drop = FALSE]
    p <- sb_plot_performance(x, metric, cohort_label = label)
    body <- c(
      body, paste0("<details open><summary>", metric, "</summary>"),
      sb_report_table(tab, c(
        "pipeline_id", "estimate", "status", "run_status",
        "scope", "aggregation", "prediction_horizon",
        "time_unit"
      ), paste("Saved", metric)),
      sb_report_plot(p, paste(boundary, label, metric), max(1, ceiling(nrow(tab) / 8))), "</details>"
    )
  }
  for (kind in c("calibration", "dca")) {
    heading <- if (kind == "calibration") "Calibration" else "Decision curves"
    p <- if (kind == "calibration") {
      sb_plot_calibration(x, cohort_label = label)
    } else {
      sb_plot_dca(x, cohort_label = label)
    }
    status <- attr(p, "stormbreaker_curve_status")
    status_counts <- as.data.frame(table(status$metric_status), stringsAsFactors = FALSE)
    names(status_counts) <- c("saved_metric_status", "pipelines")
    displayed <- unique(as.character(p$data$pipeline_id))
    status <- status[status$pipeline_id %in% displayed, , drop = FALSE]
    body <- c(
      body, paste0("<details open><summary>", heading, "</summary>"),
      sb_report_table(status_counts, names(status_counts), paste(heading, "availability counts")),
      sb_report_table(status, c(
        "pipeline_id", "metric_status", "run_status", "scope",
        "prediction_horizon", "time_unit"
      ), paste(heading, "status")),
      sb_report_plot(p, paste(boundary, label, heading), max(1, nrow(status))), "</details>"
    )
  }
  paste(c(body, "</section>"), collapse = "\n")
}

#' @rdname sb_render_report
#' @export
sb_render_report <- function(x, file, title = "Saved prediction results") {
  valid_text <- function(value) {
    is.character(value) && length(value) == 1L &&
      !is.na(value) && nzchar(trimws(value))
  }
  if (!valid_text(file) || !grepl("\\.html$", file, ignore.case = TRUE) ||
    !dir.exists(dirname(file)) || !valid_text(title)) {
    sb_abort(
      "Supply a non-empty title and a new .html path in an existing directory.",
      "stormbreaker_bad_report_input"
    )
  }
  if (file.exists(file)) {
    sb_abort(
      "The report path already exists; choose a new filename.",
      "stormbreaker_report_exists"
    )
  }
  failure_file <- sb_new_output_path(sub("\\.html$", "_failures.csv", file,
    ignore.case = TRUE
  ), "csv")
  supported <- function(value) {
    inherits(value, c(
      "stormbreaker_binary_cv_result",
      "stormbreaker_cv_result", "stormbreaker_cv_comparison", "stormbreaker_external_evaluation",
      "stormbreaker_external_comparison", "stormbreaker_repeated_cv_comparison",
      "stormbreaker_internal_external_cv", "stormbreaker_stack_cv_result"
    ))
  }
  if (supported(x)) {
    sections <- list(Result = x)
  } else {
    sections <- x
    if (!is.list(sections) || !length(sections) || is.null(names(sections)) ||
      anyDuplicated(names(sections)) || anyNA(names(sections)) ||
      any(!nzchar(trimws(names(sections)))) || !all(vapply(sections, supported, logical(1)))) {
      sb_abort(
        "Supply a supported result or a non-empty uniquely named list of results.",
        "stormbreaker_bad_report_input"
      )
    }
  }
  failures <- sb_failure_table(sections)
  expanded <- list()
  for (i in seq_along(sections)) {
    if (!sb_is_composite_view(sections[[i]])) {
      expanded <- c(expanded, stats::setNames(list(sections[[i]]), names(sections)[i]))
    } else {
      leaves <- sb_view_leaves(sections[[i]])
      for (leaf in leaves) {
        value <- leaf$value
        attr(value, "stormbreaker_view_context") <- leaf$context
        label <- paste(names(sections)[i], sb_view_label(leaf$context),
          leaf$context$view_id,
          sep = " | "
        )
        expanded <- c(expanded, stats::setNames(list(value), label))
      }
    }
  }
  sections <- expanded
  sb_require_curve_plotting()
  body <- vapply(seq_along(sections), function(i) {
    sb_report_section(sections[[i]], names(sections)[i], i)
  }, character(1))
  nav <- paste0("<li><a href='#section-", seq_along(sections), "'>",
    sb_report_escape(names(sections)), "</a></li>",
    collapse = ""
  )
  html <- paste0(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'><title>",
    sb_report_escape(title), "</title><style>",
    "*{box-sizing:border-box}body{margin:0;background:#f5f7f9;color:#172b3a;font:16px/1.55 system-ui,sans-serif}",
    "main{max-width:1180px;margin:auto;padding:28px}header{border-bottom:4px solid #24547a;padding-bottom:20px}",
    "h1{font-size:2rem;line-height:1.2}h2{margin-bottom:.3em}section{margin:32px 0;padding:22px;background:white}",
    ".boundary{color:#24547a;font-weight:650}details{margin:18px 0;border-top:1px solid #ccd6df;padding-top:12px}",
    "summary{cursor:pointer;font-weight:650;font-size:1.2rem}a{color:#24547a}a:focus-visible,summary:focus-visible,.table-scroll:focus-visible{outline:3px solid #bf6b00}",
    ".table-scroll{overflow-x:auto;max-width:100%}table{border-collapse:collapse;width:100%;font-size:.88rem}",
    "caption{text-align:left;font-weight:bold;padding:10px 0}td,th{padding:8px;text-align:left;border-bottom:1px solid #dbe2e8;overflow-wrap:anywhere}",
    "th{background:#edf2f5}figure{margin:18px 0}img{display:block;width:100%;height:auto}figcaption{font-size:.85rem;color:#465966}",
    "footer{font-size:.85rem;border-top:1px solid #ccd6df;padding:18px 0}@media(max-width:600px){main{padding:12px}section{padding:12px}h1{font-size:1.6rem}td,th{padding:5px}}",
    "@media print{body{background:white}main{max-width:none;padding:0}nav{display:none}details::details-content{content-visibility:visible!important}section{padding:0}figure{break-inside:avoid}.table-scroll{overflow:visible}}",
    "</style></head><body><main><header><h1>", sb_report_escape(title),
    "</h1>",
    "<p><a href='", sb_report_escape(utils::URLencode(basename(failure_file), reserved = TRUE)),
    "' download>Download failures CSV</a>: ", nrow(failures),
    " issue records.</p></header>",
    "<nav aria-label='Report sections'><ol>", nav, "</ol></nav>", paste(body, collapse = "\n"),
    "</main></body></html>"
  )
  if (file.exists(file)) {
    sb_abort(
      "The report path now exists; choose a new filename.",
      "stormbreaker_report_exists"
    )
  }
  failure_file <- sb_new_output_path(failure_file, "csv")
  utils::write.csv(failures, failure_file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  writeLines(enc2utf8(html), file, useBytes = TRUE)
  invisible(normalizePath(file, winslash = "/", mustWork = TRUE))
}
