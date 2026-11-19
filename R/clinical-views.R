#' Native diagnostic workflow
#' @name sb_native_diagnostics
#' @description Collect available native diagnostics and draw or paginate their saved evidence.
#' @rawRd \arguments{
#' \item{fit}{
#' \strong{For \code{sb_native_diagnostics()}:} A successful final fit with the materials needed by the requested
#' diagnostics. Extraction uses the saved model without refitting.
#' 
#' \strong{For \code{sb_plots()}:} Successful final fit (capture_diagnostics=TRUE), or a saved
#' multiclass fit from sb_analyze_methods(). Multiclass SHAP is not supported.
#' }
#' \item{complete}{
#' FALSE stops after the first available registered diagnostic;
#' TRUE attempts every registered native view and records its status.
#' }
#' \item{newdata}{
#' Optional explicit query rows for query-dependent diagnostics,
#' including nearest-neighbor, kernel-survival, BART, Gaussian-process,
#' local-LDA, FDA/PDA/MDA and rotation-forest views. Supply these rows when
#' the requested view needs queries; NULL leaves queries unspecified.
#' }
#' \item{x}{
#' A \code{stormbreaker_algorithm_diagnostic} result, from
#' \code{\link[=sb_algorithm_diagnostics]{sb_algorithm_diagnostics()}} or a family-specific diagnostic extractor such
#' as \code{\link[=sb_naive_bayes_diagnostics]{sb_naive_bayes_diagnostics()}} or \code{\link[=sb_survival_svm_diagnostics]{sb_survival_svm_diagnostics()}}.
#' }
#' \item{new_data}{
#' Observations for local diagnostics or SHAP. Required for
#' multiclass prediction neighborhoods; training rows are not substituted.
#' }
#' \item{id_col}{
#' Observation ID column for SHAP.
#' }
#' \item{background_data}{
#' Optional development reference data. Supplying this
#' and new_data explicitly requests SHAP (approximate for supported trees,
#' LogitBoost and RDA).
#' }
#' \item{plots}{
#' A ggplot or named list of ggplots, such as sb_plots(fit)$plots.
#' }
#' \item{panels_per_page}{
#' Maximum number of query or feature facets on each page.
#' }
#' \item{observations_per_page}{
#' Maximum number of labelled query observations
#' on a plot without query facets.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_native_diagnostics()}}{
#' Extract available diagnostics from a saved model
#' 
#' \strong{Returns:} A \code{stormbreaker_native_diagnostics} list with extracted \code{diagnostics}
#' and a \code{status} table with reasons and error classes. Statuses distinguish
#' \code{available}, \code{not_implemented}, \code{missing_capture}, \code{unavailable} and \code{error}.
#' Diagnostic tables may contain sample-level data.
#' }
#' 
#' \subsection{\code{sb_plot_algorithm_diagnostic()}}{
#' Plot a saved native algorithm diagnostic
#' 
#' \strong{Returns:} A ggplot displaying the requested saved native diagnostic, without
#' fitting or changing the model. Consult the extractor's documentation and
#' result metadata for the meaning of the values and any native uncertainty;
#' a development diagnostic is not external validation performance.
#' }
#' 
#' \subsection{\code{sb_plots()}}{
#' Draw available algorithm and selector evidence through one entrypoint
#' 
#' \strong{Returns:} List with plots, diagnostics/status and optional SHAP evidence.
#' Nothing is refitted. Missing evidence is recorded, never invented.
#' Explanatory metadata stays in returned results, not tiny plot footnotes.
#' }
#' 
#' \subsection{\code{sb_plot_pages()}}{
#' Paginate algorithm plots without dropping observations
#' 
#' \strong{Returns:} Named list of ggplots. All query rows and corresponding layer rows
#' are retained across pages. Explanatory notes remain in the
#' stormbreaker_notes attribute, separate from the visible panel.
#' }
#' }
#' @rawRd \seealso{
#' \code{\link[=sb_clinician_views]{sb_clinician_views()}}
#' }
NULL

# ---- clinician-views ----

#' @rdname sb_plot_calibration
#' @export
sb_clinician_views <- function(x, thresholds = NULL, confusion_threshold = NULL, native = list(),
                               shap = NULL, feature_values = NULL, id_col = NULL,
                               dependence_features = NULL, nomogram = NULL, selectors = list(),
                               detail = c("brief", "full")) {
  detail <- match.arg(detail)
  sb_require_curve_plotting()
  if (!is.list(selectors) || (length(selectors) && (is.null(names(selectors)) ||
      anyNA(names(selectors)) || anyDuplicated(names(selectors)) || any(!nzchar(trimws(names(selectors)))) ||
      !all(vapply(selectors, inherits, logical(1), what = "stormbreaker_selector_diagnostic"))))) sb_abort(
    "Supply a uniquely named list of extracted selector diagnostics.", "stormbreaker_bad_report_input")
  native_status <- NULL
  if (inherits(native, "stormbreaker_native_diagnostics")) {
    native_status <- native$status
    native <- native$diagnostics
  }
  if (!is.list(native) || (length(native) && (is.null(names(native)) || anyNA(names(native)) ||
      anyDuplicated(names(native)) || any(!nzchar(trimws(names(native))))))) sb_abort(
    "Native diagnostics must be a uniquely named list.", "stormbreaker_bad_report_input")
  if (length(native) && !all(vapply(native, inherits, logical(1), what = "stormbreaker_algorithm_diagnostic"))) sb_abort(
    "Supply extracted algorithm diagnostics, not fits or arbitrary objects.", "stormbreaker_bad_report_input")
  if (!is.null(thresholds) && (!is.numeric(thresholds) || !length(thresholds) ||
      any(!is.finite(thresholds)) || any(thresholds < 0 | thresholds > 1))) sb_abort(
    "Supply explicit thresholds in [0, 1].", "stormbreaker_bad_report_input")
  if (!is.null(confusion_threshold) && (!is.numeric(confusion_threshold) ||
      length(confusion_threshold) != 1L || !is.finite(confusion_threshold) ||
      !confusion_threshold %in% thresholds)) sb_abort(
    "The confusion threshold must be one of the supplied thresholds.", "stormbreaker_bad_report_input")
  # Validate the result container before treating per-view unavailability as data.
  models <- sb_clinician_models(sb_viz_leaves(x, NULL))
  if (!is.null(nomogram) && !inherits(nomogram, "stormbreaker_nomogram")) sb_abort(
    "Supply an exact sb_nomogram result.", "stormbreaker_bad_report_input")
  if (!is.null(shap)) sb_viz_shap_matrix(shap)
  if (!is.null(dependence_features) && (is.null(shap) || !is.character(dependence_features) ||
      !length(dependence_features) || anyNA(dependence_features) || anyDuplicated(dependence_features) ||
      !all(dependence_features %in% colnames(shap$contributions)) || is.null(feature_values) || is.null(id_col))) sb_abort(
    "Dependence plots need saved SHAP, unique explained feature names, feature values and an ID column.",
    "stormbreaker_bad_report_input")
  plots <- groups <- captions <- list()
  widths <- heights <- numeric()
  diagnostics <- list()
  status <- list()
  add <- function(label, group, draw, caption = "", height = 7) {
    value <- tryCatch(draw(), error = identity)
    if (inherits(value, "error")) {
      unavailable <- inherits(value, c("stormbreaker_bad_visualization_task", "stormbreaker_unavailable_diagnostic",
        "stormbreaker_missing_diagnostic_artifacts", "stormbreaker_unsupported_diagnostic"))
      status[[length(status) + 1L]] <<- data.frame(view = label, group = group,
        status = if (unavailable) "unavailable" else "error", reason = conditionMessage(value))
      diagnostics[[label]] <<- list(error_class = class(value), message = conditionMessage(value))
      return(invisible(NULL))
    }
    attrs <- attributes(value)
    notes <- unique(c(caption, attr(value, "stormbreaker_notes"), attr(value, "stormbreaker_figure_notes"),
      value$labels$subtitle, value$labels$caption))
    caption <- paste(notes[!is.na(notes) & nzchar(notes)], collapse = "\n")
    caption <- gsub("0 unavailable curve group(s). ", "", caption, fixed = TRUE)
    lines <- strsplit(caption, "\n", fixed = TRUE)[[1L]]
    caption <- paste(lines[!grepl("^0 saved coordinate\\(s\\) not drawn;", lines)], collapse = "\n")
    value <- value + ggplot2::labs(caption = NULL, subtitle = NULL)
    value <- sb_clinician_clean_plot(value, models)
    details <- attrs[grep("^stormbreaker_.*(status|diagnostics)$", names(attrs))]
    diagnostics[[label]] <<- details
    usable <- is.data.frame(value$data) && nrow(value$data) > 0L
    if (usable && "draw" %in% names(value$data)) usable <- any(value$data$draw %in% TRUE)
    issues <- unique(unlist(lapply(details, function(d) {
      if (!is.data.frame(d)) return(character())
      fields <- intersect(c("curve_status", "metric_status", "run_status", "status"), names(d))
      values <- as.character(unlist(d[fields], use.names = FALSE))
      values[!is.na(values) & values != "success"]
    }), use.names = FALSE))
    reason <- if (length(issues)) paste(issues, collapse = "; ") else
      if (usable) "" else "No usable coordinates; see original diagnostic tables."
    status[[length(status) + 1L]] <<- data.frame(view = label, group = group,
      status = if (!usable) "unavailable" else if (length(issues)) "partial" else "available", reason = reason)
    if (usable) {
      plots[[label]] <<- value
      groups[[label]] <<- group
      captions[[label]] <<- caption
      dimensions <- attr(value, "stormbreaker_suggested_size")
      widths[[label]] <<- if (is.null(dimensions)) 11 else dimensions[["width"]]
      heights[[label]] <<- if (is.null(dimensions)) height else dimensions[["height"]]
    }
  }
  add("ROC", "Performance", function() sb_plot_roc(x))
  if (detail == "full") add("Precision-recall", "Performance", function() sb_plot_pr(x))
  add("Calibration", "Performance", function() sb_plot_calibration(x))
  add("Decision curve", "Clinical decisions", function() sb_plot_dca(x))
  if (detail == "full") add("Probability distribution", "Clinical decisions", function() sb_plot_probability_distribution(sb_probability_distribution(x)))
  if (!is.null(thresholds)) {
    add("Threshold performance", "Clinical decisions", function() sb_plot_thresholds(sb_threshold_metrics(x, thresholds)))
    if (!is.null(confusion_threshold)) add("Confusion matrix", "Clinical decisions",
      function() sb_plot_confusion(sb_threshold_metrics(x, thresholds), confusion_threshold))
  }
  for (label in names(native)) {
    q <- native[[label]]
    add(paste0("Model: ", label), "Model diagnostics", function() sb_plot_algorithm_diagnostic(q),
      paste(c(q$source, q$interpretation), collapse = "\n"))
  }
  for (label in names(selectors)) {
    q <- selectors[[label]]
    add(paste0("Selector: ", label), "Feature selection", function() sb_plot_selector_diagnostic(q),
      q$source, height = if (q$type %in% c("vsurf_importance", "boruta_summary", "boruta_hits"))
        max(5, 2 + .3 * nrow(q$table)) else 7)
  }
  if (!is.null(native_status)) {
    diagnostics[["Native extraction"]] <- native_status
    missing <- native_status[native_status$status != "available", , drop = FALSE]
    if (nrow(missing)) status[[length(status) + 1L]] <- data.frame(
      view = paste0("Model: ", missing$view), group = "Model diagnostics",
      status = missing$status, reason = missing$reason)
  }
  if (!is.null(shap)) {
    height <- max(4.2, 1.8 + .32 * min(20L, ncol(shap$contributions)))
    add("SHAP importance", "Model explanations", function() sb_plot_shap_importance(shap), height = height)
    add("SHAP beeswarm", "Model explanations",
      function() sb_plot_shap_beeswarm(shap, feature_values = feature_values, id_col = id_col), height = height)
    for (feature in dependence_features) add(paste0("SHAP dependence: ", feature), "Model explanations",
      function() sb_plot_shap_dependence(shap, feature, feature_values, id_col), height = 5)
  }
  if (!is.null(nomogram)) add("Original-model nomogram", "Model explanations",
    function() sb_plot_nomogram(nomogram), height = max(5, 2 + .45 * nrow(nomogram$features)))
  structure(list(plots = plots, groups = unlist(groups, use.names = FALSE),
    captions = unlist(captions, use.names = FALSE), status = do.call(rbind, status),
    diagnostics = diagnostics, models = models, width = unname(widths), height = unname(heights)), class = "stormbreaker_clinician_views")
}

# ---- doctor-presentation ----

# ---- clinician-labels ----

sb_clinician_models <- function(leaves) {
  rows <- lapply(leaves, function(leaf) {
    m <- sb_page_model(leaf$value, leaf$value$task_type)
    data.frame(id = m$id, model = m$learner, selector = m$selector,
      parameters = m$parameters, selection_parameters = m$selection_parameters)
  })
  models <- unique(do.call(rbind, rows))
  models <- models[!duplicated(models$id), , drop = FALSE]
  repeated <- duplicated(models$model) | duplicated(models$model, fromLast = TRUE)
  models$model[repeated] <- paste0(models$model[repeated], " (", seq_len(sum(repeated)), ")")
  models
}

sb_clinician_labeler <- function(models) {
  function(labels) vapply(as.character(labels), function(label) {
    for (i in order(nchar(models$id), decreasing = TRUE)) {
      label <- gsub(models$id[i], models$model[i], label, fixed = TRUE)
    }
    label <- gsub("Logistic regression", "Logistic", label, fixed = TRUE)
    label <- gsub("Gaussian naive Bayes", "Gaussian NB", label, fixed = TRUE)
    lines <- trimws(strsplit(label, "\n", fixed = TRUE)[[1L]])
    lines <- lines[!grepl("^(Unnamed cohort|Role: pipeline|Saved entry:|curve: success run: success)$", lines) &
      !grepl("^Saved entry:", lines)]
    lines <- gsub("Role: pipeline", "", lines, fixed = TRUE)
    lines <- gsub("Assessment fold", "Fold", lines, fixed = TRUE)
    lines <- gsub("Development CV", "CV", lines, fixed = TRUE)
    lines <- gsub("Cross-validation", "CV", lines, fixed = TRUE)
    lines <- gsub(" | class 1", "", lines, fixed = TRUE)
    lines <- gsub("External patients", "", lines, fixed = TRUE)
    lines <- gsub("curve: success", "", lines, fixed = TRUE)
    lines <- gsub("run: success", "", lines, fixed = TRUE)
    lines <- gsub("(^[ |]+|[ |]+$)", "", lines)
    if (any(lines == "CV")) lines <- lines[lines != "Out-of-fold"]
    paste(unique(lines[nzchar(lines)]), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

sb_clinician_clean_plot <- function(plot, models) {
  clean <- sb_clinician_labeler(models)
  if (length(intersect(all.vars(plot$mapping$x), c("pipeline_id", "display_label")))) {
    plot$scales <- plot$scales$clone()
    axis <- plot$scales$get_scales("x")
    if (!is.null(axis)) axis$labels <- clean
  }
  original <- plot$facet$params$labeller
  if (!is.null(original)) {
    facet <- plot$facet
    plot$facet <- ggplot2::ggproto(NULL, facet)
    plot$facet$params$labeller <- function(labels) {
      if ("panel_key" %in% names(labels)) return(lapply(original(labels), clean))
      lapply(labels, clean)
    }
  }
  if (identical(all.vars(plot$mapping$colour), "pipeline_id")) {
    if (is.null(plot$scales$get_scales("colour"))) {
      plot <- plot + ggplot2::scale_colour_discrete(name = NULL, labels = clean)
    } else {
      plot$scales <- plot$scales$clone()
      scale <- plot$scales$get_scales("colour")
      scale$labels <- clean
      scale$name <- NULL
    }
  }
  plot
}

# ---- doctor-export ----

#' @rdname sb_render_report
#' @export
sb_run_figure_keys <- function(x) {
  as.character(names(sb_export_plot_list(x$plots)))
}

#' @rdname sb_render_report
#' @export
sb_export_run <- function(x, directory, figure_format = c("pdf", "png"),
                          dpi = 300L, font_family = NULL, figure_keys = NULL,
                          include_data = TRUE) {
  figure_format <- match.arg(figure_format)
  if (figure_format == "png") dpi <- sb_check_positive_integer(dpi, "dpi")
  if (!is.null(font_family) && (!is.character(font_family) ||
      length(font_family) != 1L || is.na(font_family) || !nzchar(trimws(font_family)))) {
    stop("font_family must be NULL or one non-empty font family.", call. = FALSE)
  }
  if (!is.logical(include_data) || length(include_data) != 1L || is.na(include_data)) {
    stop("include_data must be TRUE or FALSE.", call. = FALSE)
  }
  figures <- sb_export_plot_list(x$plots)
  if (!is.null(figure_keys)) {
    if (!is.character(figure_keys) || anyNA(figure_keys) || anyDuplicated(figure_keys)) {
      stop("figure_keys must be NULL or a character vector of unique figure keys.", call. = FALSE)
    }
    unknown <- setdiff(figure_keys, names(figures))
    if (length(unknown)) {
      stop("Unknown figure key(s): ", paste(unknown, collapse = ", "),
        ". See sb_run_figure_keys(x).", call. = FALSE)
    }
    figures <- figures[figure_keys]
  }
  figures <- sb_plot_pages(figures)
  if (file.exists(directory)) stop("Choose a new output directory; existing paths are not overwritten.")
  if (!dir.create(directory, recursive = TRUE)) stop("Could not create output directory.")
  write <- function(d, name) utils::write.csv(d, file.path(directory, name), row.names = FALSE)
  if (include_data) sb_export_run_tables(x, directory, write)
  plot_index <- data.frame(plot = as.character(names(figures)),
    file = sprintf(paste0("figure-%03d.", figure_format), seq_along(figures)))
  for (i in seq_len(nrow(plot_index))) {
    plot <- figures[[i]]
    size <- attr(plot, "stormbreaker_suggested_size")
    if (is.null(size)) {
      size <- c(width = 11, height = if (inherits(plot, "patchwork")) 12 else 8.5)
    }
    if (!is.null(font_family)) plot <- sb_export_plot_font(plot, font_family)
    path <- file.path(directory, plot_index$file[i])
    if (figure_format == "pdf") {
      ggplot2::ggsave(path, plot, device = grDevices::cairo_pdf,
        width = size[[1]], height = size[[2]])
    } else {
      ggplot2::ggsave(path, plot, device = "png", dpi = dpi, bg = "white",
        width = size[[1]], height = size[[2]])
    }
  }
  write(plot_index, "figures.csv")
  notes <- vapply(figures, function(p) {
    paste(unique(c(
      attr(p, "stormbreaker_notes"),
      p$labels$caption
    )), collapse = "\n")
  }, character(1))
  write(data.frame(plot = as.character(names(figures)), note = notes), "figure-notes.csv")
  invisible(normalizePath(directory, winslash = "/"))
}

sb_export_run_tables <- function(x, directory, write) {
  write(x$comparison, "model-comparison.csv")
  write(x$genes, "genes-by-stage.csv")
  index <- data.frame(
    cohort = as.character(names(x$evaluations)), file = sprintf("cohort-%03d", seq_along(x$evaluations)),
    status = vapply(x$evaluations, function(z) z$status, character(1))
  )
  is_holdout <- vapply(x$evaluations, function(z) {
    identical(z$assessment_role, "internal_holdout")
  }, logical(1))
  if (any(is_holdout)) index$assessment_role <- ifelse(is_holdout, "internal_holdout", "external")
  write(index, "cohorts.csv")
  if (length(x$risk_groups)) {
    write(do.call(rbind, lapply(names(x$risk_groups), function(name) {
      data.frame(cohort = name, x$risk_groups[[name]], row.names = NULL)
    })), "risk-groups.csv")
  }
  for (i in seq_len(nrow(index))) {
    name <- index$cohort[i]
    write(x$predictions[[name]], paste0(index$file[i], "-predictions.csv"))
    evaluation <- x$evaluations[[name]]
    if (is.data.frame(evaluation$scored_predictions)) {
      write(evaluation$scored_predictions, paste0(index$file[i], "-evaluated-predictions.csv"))
    }
    if (is.data.frame(evaluation$direction_correction)) {
      write(evaluation$direction_correction, paste0(index$file[i], "-direction-correction.csv"))
      original <- evaluation$original_metrics
      original_tab <- do.call(rbind, lapply(names(original), function(m) {
        z <- original[[m]]
        if (!is.list(z) || is.null(z$estimate) || length(z$estimate) != 1L) return(NULL)
        data.frame(metric = m, estimate = z$estimate,
          status = if (is.null(z$status)) NA_character_ else z$status)
      }))
      if (!is.null(original_tab)) write(original_tab, paste0(index$file[i], "-original-metrics.csv"))
    }
    metrics <- x$metrics[[name]]
    tab <- do.call(rbind, lapply(names(metrics), function(m) {
      z <- metrics[[m]]
      if (!is.list(z) || is.null(z$estimate) || length(z$estimate) != 1L) {
        return(NULL)
      }
      data.frame(
        metric = m, estimate = z$estimate,
        status = if (is.null(z$status)) NA_character_ else z$status
      )
    }))
    if (!is.null(tab)) write(tab, paste0(index$file[i], "-metrics.csv"))
  }
  write(sb_export_issue_rows(x$issues), "issues.csv")
  saveRDS(x$issues, file.path(directory, "issues.rds"))
  if (!is.null(x$explanation)) saveRDS(x$explanation, file.path(directory, "explanation.rds"))
}

sb_export_plot_font <- function(plot, family) {
  style <- ggplot2::theme(text = ggplot2::element_text(family = family))
  if (inherits(plot, "patchwork")) return(plot & style)
  plot$layers <- lapply(plot$layers, function(layer) {
    if ("family" %in% names(layer$geom$default_aes) &&
        !"family" %in% names(layer$mapping) &&
        !(isTRUE(layer$inherit.aes) && "family" %in% names(plot$mapping))) {
      params <- layer$aes_params
      params$family <- family
      return(ggplot2::ggproto(NULL, layer, aes_params = params))
    }
    layer
  })
  plot + style
}

sb_export_plot_list <- function(plots, prefix = "") {
  if (inherits(plots, "ggplot")) {
    return(stats::setNames(list(plots), prefix))
  }
  collected <- list()
  if (is.list(plots)) {
    for (name in names(plots)) {
      collected <- c(collected, sb_export_plot_list(plots[[name]], paste(prefix, name, sep = "/")))
    }
  }
  collected
}

sb_export_issue_rows <- function(issues, scope = "") {
  empty <- data.frame(scope = character(), message = character())
  if (!length(issues)) {
    return(empty)
  }
  if (is.character(issues)) {
    return(data.frame(scope = scope, message = paste(issues, collapse = "; ")))
  }
  if (is.list(issues) && !is.null(issues$status)) {
    return(data.frame(scope = scope, message = paste(c(issues$status, issues$error_message), collapse = "; ")))
  }
  if (!is.list(issues) || is.null(names(issues))) {
    return(empty)
  }
  rows <- lapply(names(issues), function(name) {
    sb_export_issue_rows(issues[[name]], paste(scope, name, sep = "/"))
  })
  do.call(rbind, rows)
}

# ---- doctor-plots ----

sb_compact_plot <- function(p) {
  notes <- unique(c(attr(p, "stormbreaker_notes"), attr(p, "stormbreaker_figure_notes"),
    p$labels$subtitle, p$labels$caption))
  title <- p$labels$title
  if (is.character(title)) {
    title <- sub("^(Original|Native|Saved) ", "", title)
    title <- paste0(toupper(substr(title, 1, 1)), substring(title, 2))
  }
  p <- p + ggplot2::labs(title = title, subtitle = NULL, caption = NULL)
  attr(p, "stormbreaker_notes") <- notes
  p
}

sb_small_model_plots <- function(fit) {
  id <- fit$pipeline$learner
  if (!id %in% c("nearest_centroid_binary", "bernoulli_naive_bayes_binary")) return(list())
  model <- fit$diagnostic_artifacts$learner
  if (is.null(model)) stop("Refit with capture_diagnostics=TRUE to retain algorithm evidence.")
  m <- if (id == "nearest_centroid_binary") model$engine$centroids else model$engine$feature_probability
  d <- data.frame(feature = rep(colnames(m), each = 2), class = rep(c("0", "1"), ncol(m)), value = as.vector(m))
  feature <- class <- value <- NULL
  list(class_profiles = ggplot2::ggplot(d, ggplot2::aes(value, feature, colour = class)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = .4), size = 2.4) +
    ggplot2::labs(x = if (id == "nearest_centroid_binary") "Centroid" else "P(feature = 1 | class)", y = NULL, colour = "Class") + sb_viz_theme())
}

#' @rdname sb_native_diagnostics
#' @export
sb_plots <- function(fit, new_data = NULL, id_col = NULL, background_data = NULL) {
  if (inherits(fit, "stormbreaker_multiclass_fit")) {
    if (!is.null(background_data)) stop("Multiclass SHAP is not implemented through sb_plots.", call. = FALSE)
    drawn <- sb_multiclass_plots(fit, sb_multiclass_evidence(fit), new_data)
    plots <- lapply(drawn$plots, sb_compact_plot)
    return(list(plots = plots, pages = sb_plot_pages(plots), diagnostics = drawn$tables,
      notes = lapply(plots, attr, which = "stormbreaker_notes")))
  }
  fit <- sb_validate_final_fit_for_prediction(fit)
  sb_require_curve_plotting()
  plots <- list()
  status <- list()
  add <- function(name, expression) {
    p <- tryCatch(force(expression), error = identity)
    if (inherits(p, "error")) status[[name]] <<- conditionMessage(p) else plots[[name]] <<- sb_compact_plot(p)
  }
  if (fit$pipeline$learner %in% c("nearest_centroid_binary", "bernoulli_naive_bayes_binary")) {
    small <- tryCatch(sb_small_model_plots(fit), error = identity)
    if (inherits(small, "error")) status$native <- conditionMessage(small) else plots <- small
    native <- NULL
  } else {
    diagnostic_data <- if (is.null(new_data)) fit$diagnostic_artifacts$predictors else new_data
    native <- sb_native_diagnostics(fit, complete = TRUE, newdata = diagnostic_data)
    missing <- native$status[native$status$status != "available", , drop = FALSE]
    if (nrow(missing)) status$native <- missing
    for (name in names(native$diagnostics)) add(name, sb_plot_algorithm_diagnostic(native$diagnostics[[name]]))
  }
  state <- fit$selection$state
  stages <- if (!is.null(state$stages)) state$stages else list(state)
  for (i in seq_along(stages)) {
    s <- stages[[i]]
    if (is.null(s) || s$id == "none") next
    if (s$id %in% sb_abess_selector_ids()) {
      criterion <- s$criterion_path
      criterion$selected <- criterion$support_size == s$selected_support_size
      diagnostic <- structure(list(type = "support_criterion", table = criterion,
        selected_support_size = s$selected_support_size), class = "stormbreaker_algorithm_diagnostic")
      add(paste0("selector_", i, "_criterion"), sb_plot_algorithm_diagnostic(diagnostic))
      if (is.null(s$coefficient_path)) status[[paste0("selector_", i, "_coefficient_path")]] <- "ABESS path was not retained; refit with capture_diagnostics=TRUE."
      else {
        b <- as.matrix(s$coefficient_path$beta)
        diagnostic$type <- "abess_coefficient_path"
        diagnostic$table <- data.frame(support_size = rep(s$coefficient_path$support_size, each = nrow(b)),
          feature = rep(s$source_features, ncol(b)), coefficient = as.vector(b))
        add(paste0("selector_", i, "_coefficient_path"), sb_plot_algorithm_diagnostic(diagnostic))
      }
    } else if (s$id %in% sb_survival_weight_selector_ids()) {
      d <- s$scores
      d$selected <- ifelse(d$selected, "Yes", "No")
      feature <- score <- selected <- NULL
      add(paste0("selector_", i), ggplot2::ggplot(d,
        ggplot2::aes(score, stats::reorder(feature, score), colour = selected)) +
        ggplot2::geom_point(size = 2.3) +
        ggplot2::scale_colour_manual(values = c(Yes = "#2878A0", No = "grey65")) +
        ggplot2::labs(x = "Screening weight", y = NULL, colour = "Selected") + sb_viz_theme())
      if (s$id == "ridge_cox_weights") {
        tab <- s$evidence$tuning_curve
        tab$mean_value <- tab$mean_deviance
        diagnostic <- list(type = "cv_curve", table = tab, measure = "CV deviance",
          selected_lambda = s$evidence$lambda, selection_rule = "lambda.min", alpha = 0)
        add(paste0("selector_", i, "_cv"), sb_plot_glmnet_tuning_diagnostic(diagnostic))
        add(paste0("selector_", i, "_coefficient_path"),
          sb_selector_path_plot(s$evidence$coefficient_path, s$evidence$lambda))
      }
      if (s$id == "plsrcox_weights" && !is.null(s$evidence$loadings)) {
        m <- s$evidence$loadings
        tab <- data.frame(feature = rep(s$source_features, ncol(m)),
          component = rep(seq_len(ncol(m)), each = nrow(m)), loading = as.vector(m))
        add(paste0("selector_", i, "_loadings"), sb_plot_pls_diagnostic(list(
          type = "pls_loadings", table = tab, features = s$source_features,
          loading_definition = "PLS loadings", interpretation = NULL, learner = "PLS-Cox")))
      }
    } else if (s$id %in% sb_survival_boosting_selector_ids()) {
      d <- s$scores
      d$selected <- ifelse(d$selected, "Yes", "No")
      feature <- score <- selected <- step <- loss <- NULL
      add(paste0("selector_", i), ggplot2::ggplot(d,
        ggplot2::aes(score, stats::reorder(feature, score), colour = selected)) +
        ggplot2::geom_point(size = 2.3, show.legend = length(unique(d$selected)) > 1L) +
        ggplot2::scale_colour_manual(values = c(Yes = "#2878A0", No = "grey65")) +
        ggplot2::labs(y = NULL, colour = "Selected",
          x = if (s$id == "xgboost_cox_importance") "Gain" else if (s$id == "gbm_cox_importance") "Importance" else
            if (s$id == "superpc_survival_scores") "SuperPC score" else
              if (s$id == "gamboost_cox_nonzero") "Spline coefficient norm" else "Coefficient") + sb_viz_theme())
      if (s$id == "superpc_survival_scores") {
        name <- paste0("selector_", i)
        plots[[name]] <- plots[[name]] + ggplot2::geom_vline(xintercept = c(-1, 1) * s$params$threshold, linetype = 2)
      } else if (nrow(s$evidence$training_loss)) add(paste0("selector_", i, "_training_loss"), ggplot2::ggplot(s$evidence$training_loss,
        ggplot2::aes(step, loss)) + ggplot2::geom_line() +
        ggplot2::labs(x = "Boosting step", y = "Training loss") + sb_viz_theme())
    } else if (s$id == "svm_rfe_binary") {
      d <- s$evidence$curve
      n_features <- value <- NULL
      for (measure in c("accuracy", "error")) {
        d$value <- d[[measure]]
        add(paste0("selector_", i, "_cv_", measure),
          ggplot2::ggplot(d, ggplot2::aes(n_features, value)) +
            ggplot2::geom_line(colour = "#2878A0") + ggplot2::geom_point() +
            ggplot2::geom_vline(xintercept = d$n_features[d$selected], linetype = 2) +
            ggplot2::labs(x = "Number of features", y = paste("CV", measure)) + sb_viz_theme())
      }
    } else if (s$id %in% c(sb_binary_model_selector_ids(), sb_pamr_selector_ids(), "lasso_nonzero", "enet_cox_nonzero", "coxboost_nonzero", "stepcox")) {
      d <- s$scores
      if (!"score" %in% names(d)) d$score <- d$coefficient
      d <- d[is.finite(d$score), ]
      d$selected <- ifelse(d$selected, "Yes", "No")
      feature <- score <- selected <- NULL
      add(paste0("selector_", i), ggplot2::ggplot(d, ggplot2::aes(score, stats::reorder(feature, score), colour = selected)) +
        ggplot2::geom_point(size = 2.3) + ggplot2::labs(x = if (s$id %in% sb_pamr_selector_ids()) "Shrunken-centroid score" else if (s$id == "plsrglm_binary_weights") "|Coefficient| \u00d7 SD" else if (s$id == "xgboost_binary_importance") "Gain" else if (s$id %in% c("rf_binary_importance", "gbm_binary_importance")) "Importance" else if (s$id == "gamboost_binary_nonzero") "Spline coefficient norm" else "Coefficient", y = NULL, colour = "Selected") + sb_viz_theme())
      if (s$id == "gbm_binary_importance") {
        trees <- loss <- NULL
        add(paste0("selector_", i, "_training_loss"), ggplot2::ggplot(s$evidence$training_loss, ggplot2::aes(trees, loss)) +
          ggplot2::geom_line() + ggplot2::labs(x = "Trees", y = "Training loss") + sb_viz_theme())
      }
      if (s$id %in% c("stepglm_binary_aic", "glmboost_binary_nonzero", "gamboost_binary_nonzero")) {
        history <- if (s$id == "stepglm_binary_aic") s$evidence$selection_path else s$evidence$training_loss
        if (is.null(history)) status[[paste0("selector_", i, "_history")]] <- "History was not retained; refit with capture_diagnostics=TRUE."
        else {
          d <- data.frame(step = seq_len(nrow(history)) - 1L,
            value = if (s$id == "stepglm_binary_aic") history$AIC else history$loss)
          step <- value <- NULL
          add(paste0("selector_", i, "_history"), ggplot2::ggplot(d, ggplot2::aes(step, value)) +
            ggplot2::geom_line() + ggplot2::geom_point() +
            ggplot2::labs(x = "Step", y = if (s$id == "stepglm_binary_aic") "AIC" else "Training loss") + sb_viz_theme())
        }
      }
      curve <- if (!is.null(s$tuning_curve)) s$tuning_curve else s$evidence$tuning_curve
      if (s$id == "coxboost_nonzero") {
        diagnostic <- structure(list(type = "coxboost_cv_curve", table = curve,
          stopping_iteration = s$optimal_step, learner = s$id, measure = "CoxBoost CV"),
          class = "stormbreaker_algorithm_diagnostic")
        add(paste0("selector_", i, "_cv"), sb_plot_coxboost_diagnostic(diagnostic))
        if (is.null(s$coefficient_path)) status[[paste0("selector_", i, "_coefficient_path")]] <- "CoxBoost path was not retained; refit with capture_diagnostics=TRUE."
        else {
          b <- s$coefficient_path
          diagnostic$type <- "coxboost_coefficient_path"
          diagnostic$measure <- "CoxBoost coefficient path"
          diagnostic$table <- data.frame(step = rep(seq_len(nrow(b)) - 1L, ncol(b)),
            feature = rep(s$source_features, each = nrow(b)), coefficient = as.vector(b))
          add(paste0("selector_", i, "_coefficient_path"), sb_plot_coxboost_diagnostic(diagnostic))
        }
      }
      selected_lambda <- if (!is.null(s$lambda)) s$lambda else s$evidence$lambda
      if (s$id %in% c("lasso_binary_nonzero", "enet_binary_nonzero", "lasso_nonzero", "enet_cox_nonzero")) {
        path <- if (!is.null(s$coefficient_path)) s$coefficient_path else s$evidence$coefficient_path
        add(paste0("selector_", i, "_coefficient_path"), sb_selector_path_plot(path, selected_lambda))
      }
      if (s$id == "rf_binary_importance") add(paste0("selector_", i, "_oob_error"),
        sb_selector_oob_plot(s$evidence$oob_error))
      if (is.data.frame(curve) && "lambda" %in% names(curve)) {
        y <- intersect(c("mean_deviance", "mean_auc", "mean_value"), names(curve))[1L]
        if (!is.na(y)) {
          cvdata <- data.frame(lambda = curve$lambda, mean = curve[[y]], se = curve$standard_error)
          lambda <- mean <- se <- NULL
          ticks <- unique(round(seq(1, nrow(curve), length.out = min(7L, nrow(curve)))))
          axis <- if ("nonzero" %in% names(curve)) ggplot2::scale_x_log10(
            sec.axis = ggplot2::dup_axis(breaks = curve$lambda[ticks],
              labels = curve$nonzero[ticks], name = "Nonzero coefficients")) else ggplot2::scale_x_log10()
          add(paste0("selector_", i, "_cv"), ggplot2::ggplot(cvdata, ggplot2::aes(lambda, mean)) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0, colour = "grey70") +
            ggplot2::geom_point(colour = "#2878B5", size = 1.5) + axis +
            ggplot2::geom_vline(xintercept = selected_lambda, linetype = 2, colour = "grey50") +
            ggplot2::labs(x = "Lambda", y = if (y == "mean_auc") "CV AUC" else "CV loss") + sb_viz_theme())
        }
      }
      if (s$id == "stepcox" && is.data.frame(s$selection_path)) {
        path <- data.frame(step = seq_len(nrow(s$selection_path)) - 1L, AIC = s$selection_path$AIC)
        step <- AIC <- NULL
        add(paste0("selector_", i, "_aic"), ggplot2::ggplot(path, ggplot2::aes(step, AIC)) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(x = "Step", y = "AIC") + sb_viz_theme())
      }
    } else {
      views <- sb_selector_score_views(s$id)
      for (v in names(views)) add(paste0("selector_", i, "_", v), sb_plot_selector_scores(s, v))
      types <- if (s$id %in% sb_boruta_selector_ids()) c("boruta_summary", "boruta_history", "boruta_boxplot") else
        if (s$id == "varselrf_binary_oob_1se") "selection_history" else
          if (s$id == "vsurf_binary_prediction") c("vsurf_importance", "vsurf_stages") else character()
      for (t in types) add(paste0("selector_", i, "_", t), sb_plot_selector_diagnostic(sb_selector_diagnostics(s, t)))
      if (!length(types) && !length(views)) status[[paste0("selector_", i)]] <- "No specialized selector plot route; selected genes remain in fit$selection."
    }
  }
  shap <- NULL
  if (!is.null(background_data)) {
    shap <- tryCatch(if (fit$pipeline$learner %in% sb_tree_shap_ids()) {
      sb_tree_shap(fit, new_data, id_col, background_data)
    } else if (fit$pipeline$learner %in% c("logitboost_binary", "rda_binary")) {
      sb_model_shap(fit, new_data, id_col, background_data)
    } else sb_additive_shap(fit, new_data, id_col, background_data), error = identity)
    if (inherits(shap, "error")) { status$shap <- conditionMessage(shap); shap <- NULL } else {
      add("shap_importance", sb_plot_shap_importance(shap))
      add("shap_beeswarm", sb_plot_shap_beeswarm(shap, feature_values = new_data, id_col = id_col))
      add("shap_waterfall", sb_plot_shap_waterfall(shap))
      plots$shap_importance <- plots$shap_importance + ggplot2::labs(title = "SHAP importance", x = "Mean |SHAP|")
      plots$shap_beeswarm <- plots$shap_beeswarm + ggplot2::labs(title = "SHAP summary", x = "SHAP value", colour = "Feature value")
      plots$shap_waterfall <- plots$shap_waterfall + ggplot2::labs(title = "Individual prediction",
        x = if (shap$score_scale == "class_1_probability") "Predicted probability" else "Model score")
      top <- utils::head(names(sort(colMeans(abs(shap$contributions)), decreasing = TRUE)), 6L)
      for (gene in top) add(paste0("shap_dependence_", gene), sb_plot_shap_dependence(shap, gene, new_data, id_col))
    }
  }
  notes <- lapply(plots, function(p) attr(p, "stormbreaker_notes"))
  if (!is.null(native)) for (name in names(native$diagnostics)) {
    notes[[name]] <- unique(c(notes[[name]], native$diagnostics[[name]]$interpretation))
    if (name %in% names(plots)) attr(plots[[name]], "stormbreaker_notes") <- notes[[name]]
  }
  list(plots = plots, pages = sb_plot_pages(plots), notes = notes,
    unavailable = status, native = native, shap = shap)
}

# ---- successful-display ----

sb_validate_show_unavailable <- function(value) {
  if(!is.logical(value) || length(value)!=1L || is.na(value)) {
    sb_abort("`show_unavailable` must be TRUE or FALSE.","stormbreaker_bad_display_option")
  }
  invisible(value)
}

sb_empty_display_plot <- function(title,subtitle=NULL,caption=NULL) {
  ggplot2::ggplot()+
    ggplot2::annotate("text",x=0,y=0,label="No usable results for this display",size=4)+
    ggplot2::labs(title=title,subtitle=subtitle,caption=caption)+
    ggplot2::theme_void(base_size=11)+
    ggplot2::theme(plot.caption=ggplot2::element_text(hjust=0))
}

sb_performance_display_rows <- function(tab) {
  !is.na(tab$status) & !is.na(tab$run_status) & tab$status=="success" &
    tab$run_status=="success" & is.finite(tab$estimate)
}

sb_curve_display_view <- function(view,kind,show_unavailable=FALSE) {
  sb_validate_show_unavailable(show_unavailable)
  if(show_unavailable) return(view)
  usable <- if(kind=="calibration") {
    unique(sb_curve_view_lines(view$data)$pipeline_id)
  } else {
    unique(view$data$pipeline_id[view$data$draw & view$data$strategy=="model"])
  }
  view$status <- view$status[view$status$pipeline_id %in% usable,,drop=FALSE]
  view$data <- view$data[view$data$pipeline_id %in% usable,,drop=FALSE]
  view$data$panel <- factor(view$data$panel,levels=view$status$panel)
  view
}

# ---- diagnostics-core ----

# ---- algorithm-diagnostics ----

#' @rdname sb_model_coefficients
#' @export
sb_algorithm_diagnostics <- function(fit, type = "coefficient_path") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("coefficient_path", "support_criterion", "cv_curve", "nonzero_path", "bic_path")) sb_abort(
    "This diagnostic type has not been implemented.", "stormbreaker_unsupported_diagnostic")
  if (type %in% c("cv_curve", "nonzero_path")) return(sb_glmnet_tuning_diagnostic(fit, type))
  if (is.null(fit$diagnostic_artifacts)) sb_abort(
    "Native artifacts were not captured. Explicitly refit on the original development task with capture_diagnostics = TRUE if this diagnostic is needed.",
    "stormbreaker_missing_diagnostic_artifacts")
  model <- fit$diagnostic_artifacts$learner
  if (inherits(model, "stormbreaker_ncvreg_model") &&
      type %in% c("coefficient_path", "bic_path")) return(sb_ncvreg_path_diagnostic(fit, type))
  if (identical(type, "bic_path")) sb_abort(
    "BIC-path diagnostics currently require an MCP/SCAD learner.",
    "stormbreaker_unsupported_diagnostic")
  engine <- if (is.list(model) && !is.null(model$engine)) model$engine else model
  if (type == "coefficient_path" && inherits(model, "stormbreaker_abess_model")) {
    b <- as.matrix(engine$beta)
    d <- data.frame(support_size = rep(engine$support.size, each = nrow(b)),
      feature = rep(model$features, ncol(b)), coefficient = as.vector(b))
    return(structure(list(type = "abess_coefficient_path", table = d,
      selected_support_size = model$selected_support_size, learner = fit$pipeline$learner,
      pipeline_id = fit$pipeline$id, source = "native ABESS support-size path"),
      class = "stormbreaker_algorithm_diagnostic"))
  }
  if (identical(type, "support_criterion")) {
    if (!inherits(model, "stormbreaker_abess_model")) sb_abort(
      "Support criterion diagnostics require the original ABESS model.",
      "stormbreaker_unsupported_diagnostic")
    path <- model$criterion_path
    selected <- model$selected_support_size
    if (!is.data.frame(path) || !all(c("support_size", "gic", "deviance") %in% names(path)) ||
        !nrow(path) || !all(vapply(path[c("support_size", "gic", "deviance")],
          function(z) is.numeric(z) && all(is.finite(z)), logical(1))) ||
        anyDuplicated(path$support_size) || any(path$support_size < 0) ||
        any(path$support_size != floor(path$support_size)) ||
        length(selected) != 1L || is.na(selected) || !selected %in% path$support_size) sb_abort(
      "Native ABESS selection path is incomplete or invalid.",
      "stormbreaker_invalid_diagnostic")
    path$selected <- path$support_size == selected
    return(structure(list(type = type, table = path,
      learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
      selected_support_size = selected,
      source = "native final-fit ABESS GIC path", engine_class = class(engine)),
      class = "stormbreaker_algorithm_diagnostic"))
  }
  if (!inherits(engine, c("glmnet", "ncvreg", "ncvsurv"))) sb_abort(
    "This learner does not have a supported native regularization path.",
    "stormbreaker_unsupported_diagnostic")
  beta <- as.matrix(engine$beta)
  lambda <- engine$lambda
  if (!is.numeric(lambda) || ncol(beta) != length(lambda) ||
      !length(lambda) || any(!is.finite(lambda)) || any(lambda <= 0) ||
      any(!is.finite(beta))) sb_abort(
    "Native regularization path is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
  labels <- rownames(beta)
  if (is.null(labels)) sb_abort("Native path feature labels are missing.",
    "stormbreaker_invalid_diagnostic")
  table <- data.frame(lambda = rep(lambda, each = nrow(beta)),
    feature = rep(labels, length(lambda)), coefficient = as.vector(beta))
  structure(list(type = type, table = table, learner = fit$pipeline$learner,
    pipeline_id = fit$pipeline$id, source = "native final-fit regularization path",
    engine_class = class(engine), tuning = fit$tuning),
    class = "stormbreaker_algorithm_diagnostic")
}

#' @rdname sb_native_diagnostics
#' @export
sb_plot_algorithm_diagnostic <- function(x) {
  sb_require_curve_plotting()
  if (inherits(x, "stormbreaker_algorithm_diagnostic") &&
      x$type %in% c("akritas_neighborhood", "akritas_survival_curves")) return(sb_plot_akritas_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type == "abess_coefficient_path") {
    support_size <- coefficient <- feature <- NULL
    return(ggplot2::ggplot(x$table, ggplot2::aes(support_size, coefficient, colour = feature)) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(xintercept = x$selected_support_size, linetype = 2, colour = "grey50") +
      ggplot2::labs(x = "Number of selected predictors", y = "Coefficient", title = "ABESS coefficient path") + sb_viz_theme())
  }
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% c("cancerclass_centroids", "cancerclass_score_mapping")) return(sb_plot_cancerclass_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type == "loclda_neighborhood") return(sb_plot_loclda_neighborhood(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% paste0("high_dimensional_", c("loadings", "densities", "spectrum", "centers"))) return(sb_plot_high_dimensional_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% paste0("mclust_da_", c("centers", "correlations", "mixture_weights"))) return(sb_plot_mclust_da_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% paste0("shrinkage_discriminant_", c("coefficients", "regularization", "loadings", "class_means"))) return(sb_plot_shrinkage_discriminant_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% c("fda_centers", "fda_query_variates", "fda_mixture_weights")) return(sb_plot_fda_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% c("gaussian_process_basis_weights", "gaussian_process_query_contributions")) return(sb_plot_gaussian_process_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% c("ferns_split_usage", "ferns_thresholds")) return(sb_plot_ferns_diagnostic(x))
  if(inherits(x,"stormbreaker_algorithm_diagnostic") && x$type %in% c("rotation_forest_loadings","rotation_forest_tree_size","rotation_forest_tree_predictions")) return(sb_plot_rotation_forest_diagnostic(x))
  if(inherits(x,"stormbreaker_algorithm_diagnostic") && x$type %in% c("c50_usage","c50_splits","c50_size")) return(sb_plot_c50_diagnostic(x))
  if (inherits(x,"stormbreaker_algorithm_diagnostic") && x$type == "classification_tree_leaf_probabilities") return(sb_plot_classification_tree_probabilities(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type == "bart_posterior") return(sb_plot_bart_posterior(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% c("bagged_knn_bootstrap", "bagged_knn_bag_curves")) return(sb_plot_bagged_knn_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type == "knn_neighborhood") return(sb_plot_knn_neighborhood(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("ahaz_coefficients", "ahaz_bic_path", "ahaz_coefficient_path")) return(sb_plot_ahaz_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("bujar_coefficients", "bujar_residual_survival")) return(sb_plot_bujar_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("penaft_coefficients", "penaft_cv_curve", "penaft_coefficient_path")) return(sb_plot_penaft_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type == "pamr_centroids") return(sb_plot_pamr_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("wsrf_importance", "wsrf_tree_nodes", "wsrf_tree_oob")) return(sb_plot_wsrf_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("rlt_split_usage", "rlt_tree_size")) return(sb_plot_rlt_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("aorsf_tree_size", "aorsf_oob_score")) return(sb_plot_aorsf_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("cforest_split_usage", "cforest_tree_size")) return(sb_plot_cforest_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("bagging_split_usage", "bagging_tree_size", "bagging_importance")) return(sb_plot_bagging_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("timereg_aalen_effects", "timereg_aalen_baseline", "timereg_odds_effects", "timereg_odds_baseline")) return(sb_plot_timereg_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("mtlr_feature_weights", "mtlr_biases")) return(sb_plot_mtlr_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && identical(x$type, "ktsp_pair_scores")) return(sb_plot_ktsp_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("flexsurv_location_coefficients", "flexsurv_spline_baseline")) return(sb_plot_flexsurv_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("aft_time_ratios", "aft_log_time_coefficients")) return(sb_plot_aft_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && identical(x$type, "ctree_structure")) return(sb_plot_ctree_structure(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("earth_basis_coefficients", "earth_pruning_gcv")) return(sb_plot_earth_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("neural_topology", "neural_weights")) return(sb_plot_neural_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("randomglm_inclusion_frequency", "randomglm_bag_coefficients", "randomglm_bag_size")) return(sb_plot_randomglm_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type == "slope_cv_curve") return(sb_plot_slope_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("lvq_prototypes", "lvq_prototype_distances")) return(sb_plot_lvq_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("discriminant_class_means", "discriminant_class_correlation")) return(sb_plot_discriminant_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type == "linear_coefficients") {
    p <- sb_plot_coefficients(x$parameters, include_zero = TRUE) +
      ggplot2::labs(subtitle = NULL, caption = NULL, x = if (x$parameters$score_scale == "class_1_log_odds")
        "Coefficient (log odds/unit)" else "Coefficient (log HR/unit)")
    attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = max(4, 2 + .3 * nrow(x$table)))
    return(p)
  }
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("naive_bayes_conditional_density", "naive_bayes_log_likelihood_ratio")) return(sb_plot_naive_bayes_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("survival_svm_score_distribution", "survival_svm_basis_weights", "survival_svm_coefficients")) return(sb_plot_survival_svm_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") &&
      x$type %in% c("svm_support_counts", "svm_decision_distribution", "svm_coefficients")) return(sb_plot_svm_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") &&
      x$type %in% c("pls_loadings", "pls_score_distribution", "pls_cv_curve")) return(sb_plot_pls_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") &&
      x$type %in% c("superpc_feature_scores", "superpc_loadings", "superpc_cv_curve")) {
    return(sb_plot_superpc_diagnostic(x))
  }
  if (inherits(x, "stormbreaker_algorithm_diagnostic") &&
      x$type %in% c("hinge_training_loss", "hinge_component_usage", "hinge_coefficient_path")) {
    return(sb_plot_hinge_boosting_diagnostic(x))
  }
  if (inherits(x, "stormbreaker_algorithm_diagnostic") &&
      x$type %in% c("classical_boosting_split_usage", "classical_boosting_tree_weights",
        "classical_boosting_tree_size", "classical_boosting_stump_sequence")) {
    return(sb_plot_classical_boosting_diagnostic(x))
  }
  if (inherits(x, "stormbreaker_algorithm_diagnostic") &&
      x$type %in% c("coxboost_coefficient_path", "coxboost_nonzero_path", "coxboost_cv_curve")) {
    return(sb_plot_coxboost_diagnostic(x))
  }
  if (!inherits(x, "stormbreaker_algorithm_diagnostic") ||
      !x$type %in% c("coefficient_path", "support_criterion", "cv_curve", "nonzero_path", "bic_path", "gam_effects", "tree_structure", "rule_coefficients", "rule_cv_curve", "forest_importance", "forest_tree_size", "forest_oob_error", "boosting_training_loss", "boosting_component_usage", "boosted_tree_importance", "boosted_tree_training_loss", "boosted_tree_cv_curve", "light_cat_importance", "light_cat_tree_size", "light_cat_interaction")) sb_abort(
    "Supply a supported algorithm diagnostic.", "stormbreaker_bad_visualization_argument")
  lambda <- coefficient <- feature <- NULL
  if (x$type %in% c("light_cat_importance", "light_cat_tree_size", "light_cat_interaction")) return(sb_plot_light_cat_diagnostic(x))
  if (x$type %in% c("boosted_tree_importance", "boosted_tree_training_loss", "boosted_tree_cv_curve")) return(sb_plot_gbm_xgboost_diagnostic(x))
  if (x$type %in% c("boosting_training_loss", "boosting_component_usage")) return(sb_plot_boosting_diagnostic(x))
  if (x$type %in% c("forest_importance", "forest_tree_size", "forest_oob_error")) return(sb_plot_forest_diagnostic(x))
  if (x$type %in% c("rule_coefficients", "rule_cv_curve")) return(sb_plot_pre_diagnostic(x))
  if (identical(x$type, "tree_structure")) return(sb_plot_rpart_tree(x))
  if (identical(x$type, "gam_effects")) return(sb_plot_gam_effects(x))
  if (identical(x$type, "bic_path")) {
    bic <- NULL
    return(ggplot2::ggplot(x$table, ggplot2::aes(lambda, bic)) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::geom_vline(xintercept = x$selected_lambda, linetype = 2, colour = "#D77935") +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = "MCP/SCAD BIC path", x = "Lambda (log scale)",
        y = "BIC (lower is better)",
        caption = "Dashed line: original selected lambda. Development-fit criterion, not cross-validation.\nNo refitting or replacement by an unpenalized model.") + sb_viz_theme())
  }
  if (x$type %in% c("cv_curve", "nonzero_path")) return(sb_plot_glmnet_tuning_diagnostic(x))
  if (identical(x$type, "support_criterion")) {
    support_size <- gic <- NULL
    return(ggplot2::ggplot(x$table, ggplot2::aes(support_size, gic)) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::geom_point(colour = "#2878A0") +
      ggplot2::geom_point(data = x$table[x$table$selected, , drop = FALSE],
        shape = 21, fill = "#D77935", size = 3) +
      ggplot2::scale_x_continuous(breaks = x$table$support_size) +
      ggplot2::labs(title = "ABESS support-size selection",
        subtitle = paste("Original selected support size:", x$selected_support_size),
        x = "Number of selected predictors", y = "GIC (lower is better)",
        caption = "Selection criterion from the original fit on the development data.") +
      sb_viz_theme())
  }
  ggplot2::ggplot(x$table, ggplot2::aes(lambda, coefficient, colour = feature)) +
    ggplot2::geom_line() + ggplot2::scale_x_log10() +
    ggplot2::labs(title = "Regularization path", x = "Lambda (log scale)",
      y = "Coefficient", colour = "Feature",
      caption = "Original final-fit path; no refitting. Not a coefficient confidence interval or a cross-validation performance curve.") +
    sb_viz_theme()
}

# ---- native-diagnostics ----

sb_native_diagnostic_plan <- function(id) {
  route <- function(extractor, types) data.frame(extractor = extractor, type = types)
  if (id == "akritas_survival") return(route("sb_akritas_diagnostics", c("neighborhood", "survival_curves")))
  if (id == "cancerclass_binary") return(route("sb_cancerclass_diagnostics", c("centroids", "score_mapping")))
  if (id == "loclda_binary") return(route("sb_loclda_neighborhood", ""))
  if (id %in% c("hda_binary", "hdda_binary_cattell")) return(route("sb_high_dimensional_diagnostics", if (id == "hda_binary") c("loadings", "densities") else c("loadings", "spectrum", "centers")))
  if (id %in% sb_mclust_da_learner_ids()) return(route("sb_mclust_da_diagnostics", c("centers", "correlations", "mixture_weights")))
  if (id %in% c("sda_binary", "slda_binary")) return(route("sb_shrinkage_discriminant_diagnostics", if (id == "sda_binary") c("coefficients", "regularization") else c("loadings", "class_means")))
  if (id %in% c("fda_binary", "pda_binary", "mda_binary")) return(route("sb_fda_diagnostics", c("centers", if (id == "mda_binary") "mixture_weights", "query_variates")))
  if (id == "gaussian_process_binary") return(route("sb_gaussian_process_diagnostics", c("basis_weights", "query_contributions")))
  if (id == "random_ferns_binary") return(route("sb_ferns_diagnostics", c("split_usage", "thresholds")))
  if(id=="rotation_forest_binary") return(route("sb_rotation_forest_diagnostics",c("loadings","tree_size","tree_predictions")))
  if (id %in% c("c50_binary","c50_rules_binary")) return(route("sb_c50_diagnostics",c("usage","splits","size")))
  if (id == "tree_binary") return(route("sb_classification_tree_diagnostics", c("structure","leaf_probabilities")))
  if (id == "evtree_binary") return(route("sb_evtree_structure", ""))
  if (id %in% c("bart_binary", "bart_survival", "bart_aft_survival")) return(route("sb_bart_posterior", ""))
  if (id == "bagged_knn_survival") return(route("sb_bagged_knn_diagnostics", c("bootstrap", "bag_curves")))
  if (id %in% c("knn_binary", "kernel_weighted_knn_binary")) return(route("sb_knn_neighborhood", ""))
  if (id == "bagging_binary") return(route("sb_bagging_diagnostics", c("split_usage", "tree_size", "importance")))
  if (id %in% c("cforest_binary", "cforest_survival")) return(route("sb_cforest_diagnostics", c("split_usage", "tree_size")))
  if (id %in% c("aorsf_binary", "aorsf_survival")) return(route("sb_aorsf_diagnostics", c("tree_size", "oob_score")))
  if (id %in% sb_rlt_learner_ids()) return(route("sb_rlt_diagnostics", c("split_usage", "tree_size")))
  if (id == "wsrf_binary") return(route("sb_wsrf_diagnostics", c("importance", "tree_nodes", "tree_oob")))
  if (id %in% c("pamr_binary", "pamr_survival")) return(route("sb_pamr_diagnostics", ""))
  if (id == "penaft_gehan_elastic_net") return(route("sb_penaft_diagnostics", c("coefficients", "cv_curve", "coefficient_path")))
  if (id == sb_bujar_linear_learner_ids()) return(route("sb_bujar_diagnostics", c("coefficients", "residual_survival")))
  if (id %in% sb_ahaz_survival_learner_ids()) return(route("sb_ahaz_diagnostics",
    if (id == "ahaz_lin_ying") "coefficients" else c("coefficients", "bic_path", "coefficient_path")))
  if (id %in% sb_timereg_survival_learner_ids()) return(route("sb_timereg_diagnostics", c("effects", "baseline")))
  if (id == "mtlr_survival") return(route("sb_mtlr_diagnostics", c("feature_weights", "biases")))
  if (id == "ktsp_binary") return(route("sb_ktsp_diagnostics", ""))
  if (id %in% sb_flexsurv_parametric_learner_ids()) return(route("sb_flexsurv_diagnostics", "location_coefficients"))
  if (id %in% sb_flexsurvspline_learner_ids()) return(route("sb_flexsurv_diagnostics", c("location_coefficients", "spline_baseline")))
  if (id %in% sb_survreg_learner_ids()) return(route("sb_aft_diagnostics", c("time_ratios", "log_time_coefficients")))
  if (id %in% c("ctree_binary", "ctree_survival")) return(route("sb_ctree_structure", ""))
  if (id == "earth_binary") return(route("sb_earth_diagnostics", c("basis_coefficients", "pruning_gcv")))
  if (id %in% c("nnet_binary", "neuralnet_multilayer_binary")) return(route("sb_neural_diagnostics", c("topology", "weights")))
  if (id == "random_glm_binary") return(route("sb_randomglm_diagnostics", c("inclusion_frequency", "bag_coefficients", "bag_size")))
  if (id == "slope_binary") return(rbind(route("sb_slope_diagnostics", ""),
    route("sb_linear_coefficients_diagnostic", "")))
  if (id %in% sb_lvq_binary_learner_ids()) return(route("sb_lvq_diagnostics", c("prototypes", "prototype_distances")))
  if (id %in% c("lda_binary", "qda_binary", "rda_binary")) return(rbind(
    if (id == "lda_binary") route("sb_linear_coefficients_diagnostic", ""),
    route("sb_discriminant_diagnostics", c("class_means", if (id != "lda_binary") "class_correlation"))))
  if (id %in% c("gaussian_naive_bayes", "kernel_naive_bayes_binary")) return(
    route("sb_naive_bayes_diagnostics", c("conditional_density", if (id == "gaussian_naive_bayes") "log_likelihood_ratio")))
  if (id %in% c("glmnet_logistic", "glmnet_cox")) return(
    route("sb_algorithm_diagnostics", c("cv_curve", "nonzero_path", "coefficient_path")))
  if (id %in% sb_ncvreg_learner_ids()) return(route("sb_algorithm_diagnostics", c("coefficient_path", "bic_path")))
  if (id %in% sb_abess_learner_ids()) return(route("sb_algorithm_diagnostics", c("support_criterion", "coefficient_path")))
  if (id %in% sb_mgcv_gam_learner_ids()) return(route("sb_gam_effects", ""))
  if (id %in% c("decision_tree_binary", "rpart_survival")) return(route("sb_tree_structure", ""))
  if (id %in% sb_pre_rule_ensemble_learner_ids()) return(route("sb_rule_ensemble", c("rule_coefficients", "rule_cv_curve")))
  if (id %in% c(sb_ranger_learner_ids(), "random_forest_binary", "rrf_binary", "rsf_survival")) return(
    route("sb_forest_diagnostics", c("importance", "tree_size", "oob_error")))
  if (id %in% c("plsrglm_binary", "plsda_binary", "splsda_binary", "plsrcox", "plsrcox_cv")) return(
    rbind(route("sb_pls_diagnostics", c("loadings", "score_distribution", if (id == "plsrcox_cv") "cv_curve")),
      if (id == "plsrglm_binary") route("sb_linear_coefficients_diagnostic", "")))
  if (id %in% c("superpc", "superpc_cv")) return(route("sb_superpc_diagnostics",
    c("feature_scores", "loadings", if (id == "superpc_cv") "cv_curve")))
  if (id %in% sb_survivalsvm_learner_ids()) return(route("sb_survival_svm_diagnostics",
    c("score_distribution", "basis_weights", "coefficients")))
  if (id %in% sb_linear_margin_diagnostic_ids()) return(route("sb_svm_diagnostics", c("coefficients", "decision_distribution")))
  if (id == "lssvm_radial_binary") return(route("sb_svm_diagnostics", "decision_distribution"))
  if (id %in% c("svm_linear_binary", "svm_radial_binary", "svm_radial_weighted_binary",
      "svm_polynomial_binary", "svm_sigmoid_binary")) return(
    route("sb_svm_diagnostics", c("support_counts", "decision_distribution", if (id == "svm_linear_binary") "coefficients")))
  if (id == "bst_hinge_binary") return(route("sb_boosting_diagnostics", c("training_loss", "component_usage", "coefficient_path")))
  if (id %in% sb_classical_binary_boosting_learner_ids()) return(route("sb_boosting_diagnostics",
    c("split_usage", if (id == "adaboost_binary") c("tree_weights", "tree_size") else "stump_sequence")))
  if (id %in% c("coxboost", "coxboost_cv")) return(route("sb_boosting_diagnostics",
    c("coefficient_path", "nonzero_path", if (id == "coxboost_cv") "cv_curve")))
  if (id %in% sb_lightgbm_binary_learner_ids()) return(route("sb_boosting_diagnostics", c("importance", "tree_size")))
  if (id %in% sb_catboost_diagnostic_ids()) return(route("sb_boosting_diagnostics", c("importance", "interaction")))
  if (id %in% sb_xgboost_learner_ids()) return(route("sb_boosting_diagnostics", "importance"))
  if (id %in% sb_gbm_diagnostic_ids()) return(route("sb_boosting_diagnostics",
    c("importance", "training_loss", if (id == "gbm_cox_cv") "cv_curve")))
  if (id %in% sb_mboost_diagnostic_ids()) return(route("sb_boosting_diagnostics", c("training_loss", "component_usage")))
  if (id %in% c("logistic", sb_extended_logistic_learner_ids(), sb_linear_survival_learner_ids())) return(
    route("sb_linear_coefficients_diagnostic", ""))
  data.frame(extractor = character(), type = character())
}

#' @rdname sb_native_diagnostics
#' @export
sb_native_diagnostics <- function(fit, complete = FALSE, newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!is.logical(complete) || length(complete) != 1L || is.na(complete)) sb_abort(
    "complete must be TRUE or FALSE.", "stormbreaker_bad_report_input")
  plan <- sb_native_diagnostic_plan(fit$pipeline$learner)
  diagnostics <- list()
  rows <- list()
  if (!nrow(plan)) rows[[1L]] <- data.frame(view = fit$pipeline$learner, extractor = "", type = "",
    status = "not_implemented", reason = "No native diagnostic route is implemented for this learner.", error_class = "")
  for (i in seq_len(nrow(plan))) {
    extractor <- get(plan$extractor[i], envir = asNamespace("stormbreaker"))
    args <- list(fit = fit)
    if ("newdata" %in% names(formals(extractor))) args$newdata <- newdata
    if (nzchar(plan$type[i])) args$type <- plan$type[i]
    value <- tryCatch(do.call(extractor, args), error = identity)
    label <- paste(fit$pipeline$learner, if (nzchar(plan$type[i])) plan$type[i] else sub("sb_", "", plan$extractor[i]), sep = ": ")
    status <- "available"
    reason <- error_class <- ""
    if (inherits(value, "error")) {
      status <- if (inherits(value, "stormbreaker_missing_diagnostic_artifacts")) "missing_capture" else
        if (inherits(value, c("stormbreaker_unavailable_diagnostic", "stormbreaker_unsupported_diagnostic"))) "unavailable" else "error"
      reason <- conditionMessage(value)
      error_class <- class(value)[1L]
    } else diagnostics[[label]] <- value
    rows[[i]] <- data.frame(view = label, extractor = plan$extractor[i], type = plan$type[i],
      status = status, reason = reason, error_class = error_class)
    if (!complete && status == "available") break
  }
  structure(list(diagnostics = diagnostics, status = do.call(rbind, rows), learner = fit$pipeline$learner,
    pipeline_id = fit$pipeline$id), class = "stormbreaker_native_diagnostics")
}
