#' Feature-selector diagnostics
#' @name sb_selector_diagnostics
#' @description Inspect and plot feature-selection evidence saved during fitting.
#' @rawRd \arguments{
#' \item{x}{
#' \strong{For \code{sb_selector_diagnostics()}:} A final fit or saved selector state.
#' 
#' \strong{For \code{sb_plot_selector_diagnostic()}:} An \code{\link[=sb_selector_diagnostics]{sb_selector_diagnostics()}} result.
#' 
#' \strong{For \code{sb_plot_selector_scores()}:} Saved selector state or a final fit with one selector. For a chain,
#' pass an individual stage, or use sb_plots to draw every stage.
#' }
#' \item{type}{
#' One of boruta_summary, boruta_hits, boruta_history, boruta_boxplot, or
#' selection_history (varSelRF's saved OOB sequence), vsurf_importance or
#' vsurf_stages. VSURF ranges are native between-forest SD, not confidence
#' intervals. Stage counts do not reconstruct unsaved stage membership.
#' }
#' \item{ncol}{
#' Columns for Boruta history facets; all predictors and shadow
#' references are shown, without top-N selection.
#' }
#' \item{metric}{
#' NULL uses the selector's default diagnostic. RSF also accepts
#' importance in addition to minimal_depth.
#' }
#' \item{top_n}{
#' Maximum number of genes shown, ordered by native rank when
#' available, otherwise by score (absolute LPC score, smallest RSF depth).
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_selector_diagnostics()}}{
#' Extract saved feature-selector evidence without fitting a classifier
#' 
#' \strong{Returns:} A selector diagnostic with the original evidence table and audit.
#' Boruta history requires final fitting with capture_diagnostics = TRUE.
#' Boruta's negative-infinity sentinel after rejection is preserved in
#' raw_importance and represented as missing in the plotted importance.
#' Missing importance iterations remain missing. Boruta importance and hits
#' are selection evidence, not classifier accuracy or causal importance.
#' }
#' 
#' \subsection{\code{sb_plot_selector_diagnostic()}}{
#' Plot original feature-selection evidence
#' 
#' \strong{Returns:} A ggplot. Boruta ranges are iteration minima/maxima, not confidence
#' intervals. varSelRF error bars are the native OOB standard deviations.
#' }
#' 
#' \subsection{\code{sb_plot_selector_scores()}}{
#' Plot saved screening scores or stability-selection probabilities
#' 
#' \strong{Returns:} A ggplot with the displayed native values in its data. The complete
#' score table is retained in attribute selector_scores. Unevaluated scores
#' remain missing and are never replaced by zero.
#' 
#' \strong{Details:} SIS supplies ranks, not invented importance values. Information
#' criteria are sequential selection scores, not classifier performance or
#' causal effects. Stability probabilities are the saved subsampling results;
#' the line is the original selection cutoff, not a newly optimized threshold.
#' Univariable binary screening displays the saved direction-adjusted AUC.
#' Univariable Cox screening displays signed Wald Z, ordered by the original
#' p-value ranking. Neither view refits the selector.
#' }
#' }
NULL

# ---- selector-views ----

# ---- selector-captured-plots ----

sb_captured_selector_adapter <- function(id, params) {
  if (id %in% sb_survival_weight_selector_ids()) return(sb_survival_weight_selector_adapter(id, params, TRUE))
  if (id == "ordered_chain") return(sb_chain_selector_adapter(params, TRUE))
  if (id %in% sb_boruta_selector_ids()) return(sb_forest_selector_adapter(id, capture_history = TRUE, options = params))
  if (id %in% sb_binary_model_selector_ids()) return(sb_binary_model_selector_adapter(id, params, TRUE))
  if (id %in% sb_abess_selector_ids()) return(sb_abess_selector_adapter(id, TRUE))
  if (id == "coxboost_nonzero") return(sb_coxboost_nonzero_selector_adapter(params$inner_folds, params$seed, TRUE))
  if (id == "coxboost_then_stepcox") {
    adapter <- sb_sequential_survival_selector_adapter(
      sb_coxboost_nonzero_selector_adapter(params$inner_folds, params$seed, TRUE),
      sb_stepcox_selector_adapter(params$direction), id)
    adapter$params <- params
    return(adapter)
  }
  if (id %in% c("lasso_nonzero", "enet_cox_nonzero")) return(
    sb_lasso_nonzero_selector_adapter(params$inner_folds, params$seed, id, TRUE))
  if (id == "lasso_then_stepcox") {
    adapter <- sb_sequential_survival_selector_adapter(
      sb_lasso_nonzero_selector_adapter(params$inner_folds, params$seed, capture_history = TRUE),
      sb_stepcox_selector_adapter(params$direction), id)
    adapter$params <- params
    return(adapter)
  }
  sb_resolve_selector_adapter(id, params)
}

sb_selector_path_plot <- function(path, selected_lambda) {
  if (is.null(path) || is.null(path$beta) || !length(path$lambda)) {
    stop("Coefficient path was not retained. Refit once with capture_diagnostics=TRUE.")
  }
  b <- as.matrix(path$beta)
  d <- data.frame(lambda = rep(path$lambda, each = nrow(b)),
    feature = rep(rownames(b), ncol(b)), coefficient = as.vector(b))
  k <- unique(round(seq(1, ncol(b), length.out = min(7L, ncol(b)))))
  lambda <- feature <- coefficient <- NULL
  ggplot2::ggplot(d, ggplot2::aes(lambda, coefficient, colour = feature)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = selected_lambda, linetype = 2, colour = "grey50") +
    ggplot2::scale_x_log10(sec.axis = ggplot2::dup_axis(breaks = path$lambda[k],
      labels = colSums(b != 0)[k], name = "Nonzero coefficients")) +
    ggplot2::labs(x = "Lambda (log scale)", y = "Coefficient", colour = "Gene") + sb_viz_theme()
}

sb_selector_oob_plot <- function(d) {
  if (is.null(d)) stop("OOB history was not retained. Refit once with capture_diagnostics=TRUE.")
  classes <- setdiff(names(d), "trees")
  long <- data.frame(trees = rep(d$trees, length(classes)),
    series = rep(classes, each = nrow(d)), error = unlist(d[classes], use.names = FALSE))
  trees <- series <- error <- NULL
  colours <- stats::setNames(grDevices::hcl.colors(length(classes), "Dark 3"), classes)
  known <- c("0" = "#2878A0", "1" = "#B86236", "oob" = "#686868", "OOB" = "#686868")
  common <- intersect(names(colours), names(known))
  colours[common] <- known[common]
  ggplot2::ggplot(long, ggplot2::aes(trees, error, colour = series)) +
    ggplot2::geom_line(linewidth = .5) + ggplot2::scale_colour_manual(values = colours) +
    ggplot2::labs(x = "Trees", y = "OOB error", colour = NULL) + sb_viz_theme()
}

# ---- selector-diagnostics ----

#' @rdname sb_selector_diagnostics
#' @export
sb_selector_diagnostics <- function(x, type = "boruta_summary") {
  if (identical(type, "boruta_boxplot")) {
    result <- sb_selector_diagnostics(x, "boruta_history")
    result$type <- "boruta_boxplot"
    return(result)
  }
  state <- if (inherits(x, "stormbreaker_final_fit")) x$selection$state else x
  if (!inherits(state, "stormbreaker_selector_fit")) sb_abort(
    "Supply a final fit with saved selector state, or the selector state itself.",
    "stormbreaker_bad_visualization_argument")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("boruta_summary", "boruta_hits", "boruta_history", "selection_history", "vsurf_importance", "vsurf_stages")) sb_abort(
    "Unsupported selector diagnostic type.", "stormbreaker_unsupported_diagnostic")
  if (type %in% c("boruta_summary", "boruta_hits", "boruta_history")) {
    if (!state$id %in% c(sb_boruta_selector_ids(), "boruta_multiclass_roughfix")) sb_abort(
      "Boruta diagnostics require an actual Boruta selector.", "stormbreaker_unsupported_diagnostic")
    tab <- state$scores
    required <- c("feature", "median_importance", "minimum_importance", "maximum_importance",
      "normalized_hits", "decision", "selected")
    if (!is.data.frame(tab) || !all(required %in% names(tab)) ||
        !identical(tab$feature, state$source_features)) sb_abort(
      "Saved Boruta feature statistics are incomplete.", "stormbreaker_invalid_diagnostic")
    if (type == "boruta_history") {
      history <- state$diagnostic_history$importance
      if (is.null(history)) sb_abort(
        "Boruta history was not captured; enable capture_diagnostics in the original final fit.",
        "stormbreaker_missing_diagnostic_artifacts")
      if (!is.matrix(history) || !is.numeric(history) || !nrow(history) ||
          is.null(colnames(history)) || anyDuplicated(colnames(history)) ||
          !all(state$source_features %in% colnames(history)) || any(history == Inf, na.rm = TRUE)) sb_abort(
        "Saved Boruta importance history is invalid.", "stormbreaker_invalid_diagnostic")
      feature <- rep(colnames(history), each = nrow(history))
      raw <- as.vector(history)
      stopped <- !is.na(raw) & raw == -Inf
      final_decision <- state$scores$decision[match(feature, state$scores$feature)]
      if (any(stopped & (is.na(final_decision) | final_decision != "Rejected"))) sb_abort(
        "Boruta rejection sentinels do not match final rejected features.", "stormbreaker_invalid_diagnostic")
      tab <- data.frame(iteration = rep(seq_len(nrow(history)), ncol(history)),
        feature = feature, raw_importance = raw,
        importance = ifelse(stopped, NA_real_, raw), not_evaluated_after_rejection = stopped,
        kind = ifelse(feature %in% state$source_features, "Predictor", "Shadow reference"),
        final_decision = final_decision)
    }
  } else if (type %in% c("vsurf_importance", "vsurf_stages")) {
    if (!state$id %in% c("vsurf_binary_prediction", "vsurf_multiclass_prediction")) sb_abort(
      "VSURF diagnostics require the original VSURF selector.", "stormbreaker_unsupported_diagnostic")
    tab <- if (type == "vsurf_importance") state$scores else state$selection_history
    required <- if (type == "vsurf_importance") c("feature", "importance", "importance_sd", "rank", "selected") else c("stage", "n_features")
    if (!is.data.frame(tab) || !nrow(tab) || !all(required %in% names(tab))) sb_abort(
      "Saved VSURF evidence is incomplete.", "stormbreaker_invalid_diagnostic")
    if (type == "vsurf_importance" && (!identical(tab$feature, state$source_features) ||
        any(!is.finite(tab$importance)) || any(!is.finite(tab$importance_sd)) || any(tab$importance_sd < 0))) sb_abort(
      "Saved VSURF importance is invalid.", "stormbreaker_invalid_diagnostic")
    if (type == "vsurf_stages" && (!identical(tab$stage, c("threshold", "interpretation", "prediction")) ||
        any(!is.finite(tab$n_features)) || any(tab$n_features < 0))) sb_abort(
      "Saved VSURF stages are invalid.", "stormbreaker_invalid_diagnostic")
  } else {
    if (!state$id %in% c(sb_varselrf_selector_ids(), "varselrf_multiclass_oob_1se")) sb_abort(
      "This selection-history adapter requires varSelRF.", "stormbreaker_unsupported_diagnostic")
    tab <- state$selection_history
    required <- c("number_of_variables", "oob_error", "oob_standard_error")
    if (!is.data.frame(tab) || !nrow(tab) || !all(required %in% names(tab)) ||
        !all(vapply(tab[required], function(z) is.numeric(z) && all(is.finite(z)), logical(1)))) sb_abort(
      "Saved varSelRF selection history is incomplete.", "stormbreaker_invalid_diagnostic")
  }
  structure(list(type = type, selector = state$id, table = tab, audit = state$audit,
    selected_features = state$selected_features,
    decision_before_rough_fix = state$diagnostic_history$decision_before_rough_fix,
    decision_after_rough_fix = state$diagnostic_history$decision_after_rough_fix,
    source = "original selector fit; no fitting or predictive performance calculation"),
    class = "stormbreaker_selector_diagnostic")
}

#' @rdname sb_selector_diagnostics
#' @export
sb_plot_selector_diagnostic <- function(x, ncol = 3L) {
  sb_require_curve_plotting()
  if (!inherits(x, "stormbreaker_selector_diagnostic")) sb_abort(
    "Supply a selector diagnostic.", "stormbreaker_bad_visualization_argument")
  ncol <- sb_viz_positive_integer(ncol, "ncol")
  d <- x$table
  feature <- median_importance <- minimum_importance <- maximum_importance <- NULL
  normalized_hits <- decision <- iteration <- importance <- kind <- NULL
  number_of_variables <- oob_error <- oob_standard_error <- importance_sd <- selected <- stage <- n_features <- NULL
  if (x$type %in% c("boruta_summary", "boruta_hits")) {
    d$feature <- factor(d$feature, levels = d$feature[order(d$median_importance)])
    p <- ggplot2::ggplot(d, ggplot2::aes(y = feature, colour = decision))
    if (x$type == "boruta_summary") {
      p <- p + ggplot2::geom_segment(ggplot2::aes(x = minimum_importance,
        xend = maximum_importance, yend = feature)) +
        ggplot2::geom_point(ggplot2::aes(x = median_importance), size = 2) +
        ggplot2::labs(x = "Native importance", title = "Boruta importance distribution summary",
          caption = "Point: median; line: observed iteration minimum to maximum, not a confidence interval.\nAll original features; decisions after any declared tentative rough fix.")
    } else {
      p <- p + ggplot2::geom_point(ggplot2::aes(x = normalized_hits), size = 2) +
        ggplot2::scale_x_continuous(limits = c(0, 1)) +
        ggplot2::labs(x = "Native normalized hit rate", title = "Boruta feature hits",
          caption = "Selection evidence from original Boruta attStats; not prediction accuracy.\nAll original features; decisions after any declared tentative rough fix.")
    }
    p <- p + ggplot2::scale_colour_manual(values = c(Confirmed = "#2878A0",
      Tentative = "#D77935", Rejected = "#909090"), drop = FALSE) +
      ggplot2::labs(y = NULL, colour = "Final decision")
  } else if (x$type == "boruta_boxplot") {
    d <- d[is.finite(d$importance), , drop = FALSE]
    if (!nrow(d)) stop("No observed Boruta importance history is available.")
    d$decision <- ifelse(d$kind == "Shadow reference", "Shadow", d$final_decision)
    medians <- tapply(d$importance, d$feature, stats::median)
    d$feature <- factor(d$feature, levels = names(sort(medians)))
    p <- ggplot2::ggplot(d, ggplot2::aes(importance, feature, fill = decision)) +
      ggplot2::geom_boxplot(width = .65, outlier.size = .7) +
      ggplot2::scale_fill_manual(values = c(Confirmed = "#2CA25F", Tentative = "#F0C75E",
        Rejected = "#D95F59", Shadow = "#4C78A8")) +
      ggplot2::labs(x = "Importance", y = NULL, fill = NULL, title = "Boruta")
  } else if (x$type == "boruta_history") {
    p <- ggplot2::ggplot(d, ggplot2::aes(iteration, importance, colour = kind)) +
      ggplot2::geom_line(na.rm = TRUE) +
      ggplot2::facet_wrap(~feature, ncol = ncol) +
      ggplot2::scale_colour_manual(values = c(Predictor = "#2878A0", "Shadow reference" = "#909090")) +
      ggplot2::labs(x = "Original iteration", y = "Native importance",
        title = "Boruta importance history", colour = NULL,
        caption = "Original predictor and shadow histories; iterations after rejection are not plotted as scores.\nMissing entries remain missing. Selection process evidence, not classifier performance.")
  } else if (x$type == "selection_history") {
    p <- ggplot2::ggplot(d, ggplot2::aes(number_of_variables, oob_error)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = oob_error - oob_standard_error,
        ymax = oob_error + oob_standard_error), width = 0, colour = "#909090") +
      ggplot2::geom_line(colour = "#2878A0") + ggplot2::geom_point(colour = "#2878A0") +
      ggplot2::geom_vline(xintercept = length(x$selected_features), linetype = 2) +
      ggplot2::labs(x = "Number of variables", y = "OOB error",
        title = "varSelRF elimination history",
        caption = "Bars: native OOB standard deviation, not confidence intervals. Dashed: selected subset size.\nSelection-stage OOB error, not untouched-cohort performance.")
  } else if (x$type == "vsurf_importance") {
    d$feature <- factor(d$feature, levels = rev(d$feature[order(d$rank)]))
    p <- ggplot2::ggplot(d, ggplot2::aes(importance, feature, colour = selected)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey75", linewidth = .3) +
      ggplot2::geom_segment(ggplot2::aes(x = importance - importance_sd,
        xend = importance + importance_sd, yend = feature), linewidth = .4) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_colour_manual(values = c("FALSE" = "#909090", "TRUE" = "#2878A0"),
        breaks = c(FALSE, TRUE), labels = c("Not retained", "Retained")) +
      ggplot2::labs(title = "VSURF importance", x = "Native importance", y = NULL, colour = NULL,
        caption = "Points: original forest-averaged importance; bars: plus/minus one native between-forest SD, not confidence intervals.\nAll input features shown; colour indicates the final prediction-stage subset, not stage-specific membership or clinical importance.")
  } else if (x$type == "vsurf_stages") {
    d$stage <- factor(d$stage, levels = c("threshold", "interpretation", "prediction"),
      labels = c("Threshold", "Interpretation", "Prediction"))
    p <- ggplot2::ggplot(d, ggplot2::aes(stage, n_features)) +
      ggplot2::geom_col(fill = "#2878A0", width = .6) +
      ggplot2::scale_y_continuous(breaks = function(limits) unique(round(pretty(limits)))) +
      ggplot2::labs(title = "VSURF selection stages", x = NULL, y = "Retained features",
        caption = "Original threshold, interpretation and prediction-stage subset sizes.\nCounts alone do not identify each intermediate subset; no missing membership or OOB sequence is reconstructed.")
  } else sb_abort("Unsupported selector diagnostic.", "stormbreaker_unsupported_diagnostic")
  p + sb_viz_theme()
}

# ---- selector-score-plots ----

sb_selector_score_views <- function(id) {
  if (id == "univ_binary_auc_top_k") return(c(discrimination = "AUC (direction-adjusted)"))
  if (id == "univ_cox_top_k") return(c(wald_z = "Wald Z"))
  if (id == "rsf_minimal_depth") return(c(minimal_depth = "Minimal depth", importance = "Importance"))
  if (id == "sis_binary_nlogn") return(c(rank = "Screening rank"))
  if (id == "mrmr_binary_nlogn") return(c(mrmr_score = "mRMR score"))
  if (id %in% c("jmi_binary_nlogn", "cmim_binary_nlogn", "disr_binary_nlogn", "jmim_binary_nlogn")) return(c(information_score = "Selection score"))
  if (id %in% c("lpc_binary_nlogn", "lpc_cox_nlogn")) return(c(lpc_score = "LPC score"))
  if (id %in% c("stabsel_binary_lasso", "stabsel_cox_lasso")) return(c(selection_probability = "Selection probability"))
  character()
}

#' @rdname sb_selector_diagnostics
#' @export
sb_plot_selector_scores <- function(x, metric = NULL, top_n = 20L) {
  sb_require_curve_plotting()
  s <- if (inherits(x, "stormbreaker_final_fit")) x$selection$state else x
  if (!inherits(s, "stormbreaker_selector_fit")) stop("Supply a saved selector state.")
  views <- sb_selector_score_views(s$id)
  if (is.null(metric)) metric <- names(views)[1L]
  if (length(metric) != 1L || is.na(metric) || !metric %in% names(views)) stop("No supported score view for this selector.")
  top_n <- sb_check_positive_integer(top_n, "top_n", minimum = 1L)
  all <- s$scores
  if (!all(c("feature", "selected", metric) %in% names(all))) stop("The original score table is incomplete.")
  d <- all[is.finite(all[[metric]]), , drop = FALSE]
  if (!nrow(d)) stop("No finite native scores are available.")
  priority <- if ("rank" %in% names(d) && metric != "importance") d$rank else
    if (metric == "minimal_depth") d[[metric]] else
      if (metric == "lpc_score") -abs(d[[metric]]) else -d[[metric]]
  d <- utils::head(d[order(priority, seq_len(nrow(d)), na.last = TRUE), , drop = FALSE], top_n)
  d$value <- d[[metric]]
  d$feature <- factor(d$feature, levels = rev(d$feature))
  d$selection <- ifelse(d$selected, "Selected", "Not selected")
  value <- feature <- selection <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(value, feature, colour = selection)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(values = c(Selected = "#2878B5", `Not selected` = "grey65")) +
    ggplot2::labs(x = unname(views[metric]), y = NULL, colour = NULL) + sb_viz_theme()
  if (metric == "selection_probability") p <- p +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::geom_vline(xintercept = s$params$cutoff, linetype = 2, colour = "grey50")
  attr(p, "selector_scores") <- all
  p
}

# ---- selection-summary ----

# ---- gene-overlap ----

#' @rdname sb_analyze_methods
#' @export
sb_gene_overlap <- function(sets) {
  if (!is.list(sets) || length(sets) < 2L || is.null(names(sets)) ||
      anyNA(names(sets)) || anyDuplicated(names(sets)) || any(!nzchar(trimws(names(sets))))) {
    stop("Supply at least two uniquely named gene sets.")
  }
  if (!all(vapply(sets, function(s) is.character(s) && !anyNA(s) &&
      all(nzchar(trimws(s))), logical(1)))) stop("Gene sets must contain nonmissing, nonempty character IDs.")
  sets <- lapply(sets, unique)
  genes <- unique(unlist(sets, use.names = FALSE))
  membership <- matrix(FALSE, length(genes), length(sets), dimnames = list(genes, names(sets)))
  for (j in seq_along(sets)) membership[, j] <- genes %in% sets[[j]]
  keys <- if (length(genes)) apply(membership, 1L, paste0, collapse = ":") else character()
  members <- split(genes, keys)
  combinations <- lapply(members, function(g) names(sets)[membership[g[1L], ]])
  intersections <- data.frame(intersection = vapply(combinations, paste, character(1), collapse = " & "),
    n = lengths(members), row.names = NULL)
  intersections$sets <- unname(combinations)
  intersections$genes <- unname(members)
  intersections <- intersections[order(-intersections$n, intersections$intersection), , drop = FALSE]
  list(sets = sets, membership = membership, intersections = intersections,
    common = Reduce(intersect, sets))
}

#' @rdname sb_analyze_methods
#' @export
sb_plot_gene_overlap <- function(x, type = c("auto", "venn", "bar")) {
  type <- match.arg(type)
  x <- sb_gene_overlap(x$sets)
  if (type == "auto") type <- if (length(x$sets) <= 4L) "venn" else "bar"
  if (type == "venn") {
    if (length(x$sets) > 4L) stop("Use bar for more than four sets.")
    if (!requireNamespace("ggvenn", quietly = TRUE)) stop("Install ggvenn for Venn diagrams, or use type='bar'.")
    p <- ggvenn::ggvenn(x$sets, show_percentage = FALSE,
      fill_color = c("#6497B1", "#E69F80", "#8FBF9F", "#B4A7D6"), stroke_size = .5)
  } else {
    sb_require_curve_plotting()
    d <- x$intersections
    d$intersection <- factor(d$intersection, levels = rev(d$intersection))
    intersection <- n <- NULL
    p <- ggplot2::ggplot(d, ggplot2::aes(n, intersection)) +
      ggplot2::geom_col(fill = "#2878A0", width = .65) +
      ggplot2::labs(x = "Genes", y = NULL) + sb_viz_theme()
  }
  attr(p, "gene_overlap") <- x
  p
}

# ---- gene-selection-ranking ----

#' @rdname sb_analyze_methods
#' @export
sb_gene_selection_ranking <- function(x, source = c("selector", "glmnet_active_set"),
                                      candidate_features = NULL, cohort_label = NULL) {
  source <- match.arg(source)
  if (inherits(x, "stormbreaker_cohort_comparison")) x <- x$development
  supported <- c(
    "stormbreaker_binary_cv_result", "stormbreaker_cv_result",
    "stormbreaker_cv_comparison", "stormbreaker_repeated_cv_comparison"
  )
  if (!inherits(x, supported)) {
    sb_abort("Supply one development CV result or comparison.", "stormbreaker_bad_gene_ranking")
  }
  leaves <- sb_viz_leaves(x, cohort_label)
  tasks <- unique(vapply(leaves, function(leaf) {
    if (inherits(leaf$value, "stormbreaker_binary_cv_result")) "binary" else "survival"
  }, character(1)))
  if (length(tasks) != 1L) sb_abort("Rank binary and survival tasks separately.", "stormbreaker_bad_gene_ranking")
  valid_genes <- function(z) {
    is.character(z) && length(z) > 0L &&
      !anyNA(z) && all(nzchar(trimws(z))) && !anyDuplicated(z)
  }
  saved_candidates <- unique(unlist(lapply(leaves, function(leaf) leaf$value$task$predictors), use.names = FALSE))
  if (is.null(candidate_features)) candidate_features <- saved_candidates
  if (!valid_genes(candidate_features)) sb_abort("Supply distinct non-empty candidate gene names.", "stormbreaker_bad_gene_ranking")
  groups <- list()
  pipeline_rows <- list()
  for (i in seq_along(leaves)) {
    leaf <- leaves[[i]]
    result <- leaf$value
    pipeline <- result$pipeline
    summary <- if (source == "selector") result$selection_stability else result$glmnet_active_set_stability
    applicable <- if (source == "selector") {
      !is.null(pipeline$selector) && !pipeline$selector %in% c("none", "fixed_features")
    } else {
      pipeline$learner %in% c("glmnet_cox", "glmnet_logistic") &&
        isTRUE(pipeline$learner_params$alpha > 0)
    }
    status <- if (!applicable) "not_applicable" else if (is.null(summary)) "not_saved" else summary$status
    tab <- if (applicable) summary$feature_frequency else NULL
    n <- if (applicable) summary$n_successful_folds else 0L
    valid <- is.data.frame(tab) && all(c("feature", "selected_count") %in% names(tab)) &&
      valid_genes(tab$feature) && is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 0 && n == floor(n) &&
      is.numeric(tab$selected_count) && all(is.finite(tab$selected_count)) &&
      all(tab$selected_count >= 0 & tab$selected_count <= n & tab$selected_count == floor(tab$selected_count))
    if (applicable && !is.null(tab) && !valid) {
      sb_abort("Saved selection counts are malformed.", "stormbreaker_bad_gene_ranking")
    }
    usable <- applicable && valid && n > 0
    group_id <- NA_integer_
    duplicate <- FALSE
    if (usable) {
      tab <- tab[order(tab$feature), c("feature", "selected_count"), drop = FALSE]
      rownames(tab) <- NULL
      assignment <- result$provenance$fold_assignment
      if (is.data.frame(assignment) && all(c("patient_id", "fold") %in% names(assignment))) {
        assignment <- assignment[order(assignment$patient_id), c("patient_id", "fold"), drop = FALSE]
        rownames(assignment) <- NULL
      } else {
        assignment <- result$folds
      }
      identity_known <- !is.null(assignment) && length(assignment) > 0L
      key <- list(
        source = source, selector = pipeline$selector,
        selector_params = pipeline$selector_params,
        learner = if (source == "glmnet_active_set") pipeline$learner else NULL,
        learner_params = if (source == "glmnet_active_set") pipeline$learner_params else NULL,
        candidates = sort(result$task$predictors), assignment = assignment,
        repeat_id = leaf$context$repeat_id,
        separate_entry = if (!identity_known) i else NULL
      )
      matches <- which(vapply(groups, function(group) identical(group$key, key), logical(1)))
      counts <- list(
        table = tab, denominator = as.integer(n),
        per_fold = summary$per_fold[intersect(c("fold", "status", "n_selected"), names(summary$per_fold))]
      )
      if (length(matches)) {
        group_id <- matches[1]
        duplicate <- TRUE
        groups[[group_id]]$pipeline_ids <- c(groups[[group_id]]$pipeline_ids, pipeline$id)
        if (!identical(groups[[group_id]]$counts, counts)) groups[[group_id]]$status <- "conflicting_counts"
      } else {
        group_id <- length(groups) + 1L
        groups[[group_id]] <- list(
          key = key, counts = counts, status = "success",
          selector = pipeline$selector, learner = pipeline$learner,
          pipeline_ids = pipeline$id, split_identity_saved = identity_known
        )
      }
    }
    j <- if (valid) match(candidate_features, tab$feature) else rep(NA_integer_, length(candidate_features))
    pipeline_rows[[i]] <- data.frame(
      pipeline_id = pipeline$id, source = source,
      feature = candidate_features, selected_count = if (usable) tab$selected_count[j] else NA_integer_,
      eligible_count = if (usable) ifelse(is.na(j), NA_integer_, as.integer(n)) else NA_integer_,
      status = if (is.null(status)) "not_saved" else status,
      group_id = group_id, duplicate_copy = duplicate, stringsAsFactors = FALSE
    )
  }
  successful_groups <- Filter(function(group) identical(group$status, "success"), groups)
  selected <- eligible <- methods <- integer(length(candidate_features))
  for (group in successful_groups) {
    j <- match(candidate_features, group$counts$table$feature)
    present <- !is.na(j)
    count <- group$counts$table$selected_count[j[present]]
    selected[present] <- selected[present] + count
    eligible[present] <- eligible[present] + group$counts$denominator
    methods[present] <- methods[present] + as.integer(count > 0)
  }
  selected[eligible == 0] <- NA_integer_
  methods[eligible == 0] <- NA_integer_
  ranking <- data.frame(
    rank = as.integer(rank(-selected, ties.method = "min", na.last = "keep")),
    feature = candidate_features, selected_count = selected, eligible_count = eligible,
    selection_frequency = ifelse(eligible > 0, selected / eligible, NA_real_),
    selected_methods = methods, stringsAsFactors = FALSE
  )
  ranking <- ranking[order(-ranking$selected_count, seq_len(nrow(ranking)), na.last = TRUE), , drop = FALSE]
  rownames(ranking) <- NULL
  group_table <- if (length(groups)) {
    sb_view_bind(lapply(seq_along(groups), function(i) {
      group <- groups[[i]]
      data.frame(
        group_id = i, source = source, selector = group$selector,
        learner = if (source == "glmnet_active_set") group$learner else NA_character_,
        status = group$status, successful_folds = group$counts$denominator,
        pipeline_copies = length(group$pipeline_ids),
        split_identity_saved = group$split_identity_saved,
        pipeline_ids = paste(unique(group$pipeline_ids), collapse = ";"), stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(
      group_id = integer(), source = character(), selector = character(),
      learner = character(), status = character(), successful_folds = integer(),
      pipeline_copies = integer(), split_identity_saved = logical(), pipeline_ids = character()
    )
  }
  by_pipeline <- sb_view_bind(pipeline_rows)
  status <- if (!length(successful_groups)) {
    "not_available"
  } else if (
    any(group_table$status != "success") || any(by_pipeline$status %in% c("incomplete", "not_saved", "failed"))) {
    "incomplete"
  } else {
    "success"
  }
  structure(list(
    status = status, source = source, task_type = tasks,
    cohort = leaves[[1]]$dataset, counting_unit = "unique_configuration_successful_outer_folds",
    ranking = ranking, by_pipeline = by_pipeline, groups = group_table
  ), class = "stormbreaker_gene_selection_ranking")
}

#' @rdname sb_analyze_methods
#' @export
sb_plot_gene_selection_ranking <- function(x, top_n = 20L) {
  if (!inherits(x, "stormbreaker_gene_selection_ranking")) {
    sb_abort("Supply a gene-selection ranking.", "stormbreaker_bad_gene_ranking")
  }
  sb_require_curve_plotting()
  top_n <- sb_viz_positive_integer(top_n, "top_n")
  d <- utils::head(x$ranking[is.finite(x$ranking$selected_count), , drop = FALSE], top_n)
  if (!nrow(d)) {
    return(sb_empty_display_plot("Gene selection counts", "No saved learned-selection counts.", ""))
  }
  d$feature_axis <- factor(d$feature, levels = rev(d$feature))
  d$count_label <- paste(d$selected_count, "/", d$eligible_count)
  feature_axis <- selected_count <- count_label <- NULL
  ggplot2::ggplot(d, ggplot2::aes(x = selected_count, y = feature_axis)) +
    ggplot2::geom_col(fill = "#2878B5", width = 0.62) +
    ggplot2::geom_text(ggplot2::aes(label = count_label), hjust = -0.15, family = "Arial", size = 3.5) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.23)),
      breaks = function(limits) {
        b <- pretty(limits)
        b[b >= 0 & b == floor(b)]
      }
    ) +
    ggplot2::labs(
      title = if (x$source == "selector") "Gene selection counts" else "LASSO / Elastic Net nonzero counts",
      x = "Selected / eligible folds", y = NULL,
      subtitle = if (x$status == "incomplete") "Incomplete records" else NULL
    ) +
    sb_shap_theme()
}

# ---- selection-stability ----

sb_empty_selection_pairwise <- function() {
  data.frame(
    fold_i = integer(),
    fold_j = integer(),
    status = character(),
    n_i = integer(),
    n_j = integer(),
    intersection = integer(),
    union = integer(),
    jaccard = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_selection_stability <- function(selector, selections, candidate_features) {
  fold_ids <- vapply(
    selections,
    function(selection) as.integer(selection$fold[[1L]]),
    integer(1)
  )
  ordered <- order(fold_ids)
  fold_ids <- fold_ids[ordered]
  selections <- selections[ordered]
  n_folds <- as.integer(length(selections))
  not_applicable <- selector %in% c("none", "fixed_features")

  saved_success <- vapply(
    selections,
    function(selection) identical(selection$status, "success"),
    logical(1)
  )
  successful <- if (not_applicable) {
    rep(FALSE, n_folds)
  } else {
    saved_success
  }
  selected_sets <- lapply(
    selections,
    function(selection) unique(as.character(selection$selected_features))
  )

  error_field <- function(selection, field) {
    value <- selection[[field]]
    if (is.null(value) || length(value) == 0L || is.na(value[[1L]])) {
      return(NA_character_)
    }
    as.character(value[[1L]])
  }
  per_fold <- data.frame(
    fold = as.integer(fold_ids),
    status = if (not_applicable) {
      rep("not_applicable", n_folds)
    } else {
      vapply(selections, function(selection) selection$status, character(1))
    },
    n_selected = as.integer(vapply(seq_along(selections), function(index) {
      if (successful[[index]]) length(selected_sets[[index]]) else NA_integer_
    }, integer(1))),
    error_class = vapply(
      selections,
      error_field,
      character(1),
      field = "error_class"
    ),
    error_message = vapply(
      selections,
      error_field,
      character(1),
      field = "error_message"
    ),
    stringsAsFactors = FALSE
  )
  if (not_applicable) {
    per_fold$error_class <- NA_character_
    per_fold$error_message <- NA_character_
  }

  n_successful <- as.integer(sum(successful))
  selected_count <- vapply(candidate_features, function(feature) {
    as.integer(sum(vapply(
      selected_sets[successful],
      function(selected) feature %in% selected,
      logical(1)
    )))
  }, integer(1))
  feature_frequency <- data.frame(
    feature = as.character(candidate_features),
    selected_count = as.integer(selected_count),
    selection_frequency = if (n_successful > 0L) {
      unname(selected_count / n_successful)
    } else {
      rep(NA_real_, length(candidate_features))
    },
    stringsAsFactors = FALSE
  )

  pair_indices <- if (n_folds >= 2L) {
    t(utils::combn(seq_len(n_folds), 2L))
  } else {
    matrix(integer(), ncol = 2L)
  }
  pairwise <- if (nrow(pair_indices) == 0L) {
    sb_empty_selection_pairwise()
  } else {
    rows <- lapply(seq_len(nrow(pair_indices)), function(row) {
      index_i <- pair_indices[row, 1L]
      index_j <- pair_indices[row, 2L]
      status <- if (not_applicable) {
        "not_applicable"
      } else if (successful[[index_i]] && successful[[index_j]]) {
        "success"
      } else {
        "incomplete"
      }
      n_i <- per_fold$n_selected[[index_i]]
      n_j <- per_fold$n_selected[[index_j]]
      intersection_size <- NA_integer_
      union_size <- NA_integer_
      jaccard <- NA_real_
      if (identical(status, "success")) {
        intersection_size <- as.integer(length(intersect(
          selected_sets[[index_i]],
          selected_sets[[index_j]]
        )))
        union_size <- as.integer(length(union(
          selected_sets[[index_i]],
          selected_sets[[index_j]]
        )))
        if (union_size > 0L) {
          jaccard <- unname(intersection_size / union_size)
        } else {
          status <- "not_evaluable"
        }
      }
      data.frame(
        fold_i = as.integer(fold_ids[[index_i]]),
        fold_j = as.integer(fold_ids[[index_j]]),
        status = status,
        n_i = as.integer(n_i),
        n_j = as.integer(n_j),
        intersection = intersection_size,
        union = union_size,
        jaccard = jaccard,
        stringsAsFactors = FALSE
      )
    })
    combined <- do.call(rbind, rows)
    rownames(combined) <- NULL
    combined
  }

  n_pairs <- as.integer(nrow(pairwise))
  n_evaluable_pairs <- as.integer(sum(pairwise$status == "success"))
  complete <- !not_applicable && n_folds >= 2L && all(successful) &&
    n_evaluable_pairs == n_pairs
  status <- if (not_applicable) {
    "not_applicable"
  } else if (complete) {
    "success"
  } else {
    "incomplete"
  }
  error_message <- if (not_applicable || complete) {
    NA_character_
  } else if (n_folds < 2L) {
    "At least two outer-fold selections are required for stability."
  } else if (!all(successful)) {
    paste0(
      "Selection stability is incomplete because not every outer-fold ",
      "selection succeeded."
    )
  } else {
    "Selection stability is incomplete because a fold pair was not evaluable."
  }

  list(
    selector = selector,
    status = status,
    scope = "outer_analysis_fold_selected_feature_sets",
    aggregation = "mean_of_unique_fold_pair_jaccard",
    frequency_denominator = "successful_selection_folds",
    n_candidate_features = as.integer(length(candidate_features)),
    n_folds = n_folds,
    n_successful_folds = n_successful,
    n_pairs = n_pairs,
    n_evaluable_pairs = n_evaluable_pairs,
    mean_pairwise_jaccard = if (complete) {
      unname(mean(pairwise$jaccard))
    } else {
      NA_real_
    },
    feature_frequency = feature_frequency,
    per_fold = per_fold,
    pairwise = pairwise,
    error_message = error_message
  )
}
