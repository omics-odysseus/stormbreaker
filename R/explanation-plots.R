#' Linear-model coefficients and paths
#' @name sb_model_coefficients
#' @description Inspect saved coefficients, regularization paths and related native tuning evidence.
#' @rawRd \arguments{
#' \item{fit}{
#' \strong{For \code{sb_model_coefficients()}:} Successful portable linear final fit.
#' 
#' \strong{For \code{sb_algorithm_diagnostics()}:} Final model. Native paths require capture_diagnostics = TRUE;
#' saved glmnet tuning curves do not require extra capture.
#' 
#' \strong{For \code{sb_slope_diagnostics()}:} A successful SLOPE final fit captured with capture_diagnostics.
#' 
#' \strong{For \code{sb_randomglm_diagnostics()}:} A successful randomGLM final fit.
#' }
#' \item{reference_data}{
#' Optional original-training feature data. Only fitted
#' numeric features are read. At least two rows are needed for sample SD.
#' }
#' \item{x}{
#' A result from \code{\link[=sb_model_coefficients]{sb_model_coefficients()}}.
#' }
#' \item{scale}{
#' Raw input-unit coefficients or coefficients multiplied by the
#' supplied reference sample SD. Predictive models are never modified.
#' }
#' \item{top_n}{
#' Maximum number of terms shown, or NULL for all.
#' }
#' \item{include_zero}{
#' Whether to display zero coefficients (default FALSE).
#' }
#' \item{type}{
#' \strong{For \code{sb_algorithm_diagnostics()}:} Diagnostic type: coefficient_path for native glmnet, ABESS or ncvreg
#' engines, support_criterion for the ABESS GIC selection path, or cv_curve
#' and nonzero_path for the saved glmnet inner-validation sequence, or
#' bic_path for MCP/SCAD's original BIC selection sequence.
#' 
#' \strong{For \code{sb_randomglm_diagnostics()}:} inclusion_frequency, bag_coefficients or bag_size.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_model_coefficients()}}{
#' Inspect the parameters of a saved linear model
#' 
#' \strong{Returns:} A \code{stormbreaker_coefficients} list with a feature table, intercept
#' (NULL for Cox), score scale and reference sample count. Zero coefficients
#' remain in the table. Scaled coefficients are beta times sample SD, not a
#' replacement model and not a SHAP importance or causal effect.
#' Stepwise GLM, GLMBoost and plsRglm retain their fitted slopes and convert
#' the original centered intercept to input units. Bias-reduced and Cauchy-prior
#' logistic fits retain their original estimates, not a new ordinary GLM.
#' LDA retains the exact posterior log-odds equation from its original class
#' priors, means and discriminant scaling; no logistic model is fitted.
#' SLOPE uses its original selected penalized logistic coefficients in input units.
#' LIBLINEAR L1/L2 logistic models undo their original standardization and
#' orient native weights toward class 1, including the original bias term.
#' Their coefficients are penalized log-odds, without ordinary GLM intervals.
#' Robust Mqle logistic models retain their original robust estimates;
#' neither ordinary logistic refitting nor robust weight reconstruction occurs.
#' Survreg AFT models use the original ranking score, negative log-time
#' location: both slopes and intercept are the negatives of native AFT values.
#' }
#' 
#' \subsection{\code{sb_plot_coefficients()}}{
#' Plot signed coefficients, optionally per reference standard deviation
#' 
#' \strong{Returns:} A ggplot object. The table retains all coefficients regardless of
#' display choices. Direction describes conditional association, not causation.
#' }
#' 
#' \subsection{\code{sb_algorithm_diagnostics()}}{
#' Extract native regularization paths from a captured final model
#' 
#' \strong{Returns:} A \code{stormbreaker_algorithm_diagnostic} list with the native path
#' table and origin metadata. No fitting or parameter search occurs.
#' }
#' 
#' \subsection{\code{sb_slope_diagnostics()}}{
#' Inspect the original SLOPE cross-validation search
#' 
#' \strong{Returns:} An algorithm diagnostic retaining the full native q/gamma/alpha
#' summary, means, standard errors, bounds and selected optimum. Paths may
#' stop early; no uncomputed alpha values are added. Bounds are the native
#' mean plus/minus qt(0.975, score-row-count - 1) times native SE, not external
#' validation uncertainty. Training observations are not retained in this
#' added diagnostic summary. Historical compact models have only their saved
#' optimum and cannot supply a complete CV path. Original coefficients remain
#' available through \code{\link[=sb_model_coefficients]{sb_model_coefficients()}}, including without capture.
#' }
#' 
#' \subsection{\code{sb_randomglm_diagnostics()}}{
#' Inspect original randomGLM ensemble members
#' 
#' \strong{Returns:} An algorithm diagnostic from compact saved ensemble members. Every
#' fitted feature and bag is retained, including intercept-only bags. Inclusion
#' means membership in a bag's selected feature list, not nonzero coefficient.
#' Coefficients are individual bag log-odds slopes in original input units;
#' omitted features have zero contribution and a separate included flag.
#' Bag frequency is not cross-validation stability, population prevalence,
#' uncertainty or importance. The ensemble averages probabilities, not slopes:
#' these coefficients must not be averaged into a surrogate nomogram.
#' No training data, refitting, bag resampling or new model is required.
#' }
#' }
NULL

# ---- plot-shap ----

# ---- plot-shap-dependence ----

#' @rdname sb_model_shap
#' @export
sb_plot_shap_dependence <- function(x, feature, feature_values, id_col, colour_by = NULL) {
  m <- sb_viz_shap_matrix(x)
  valid_name <- function(z) is.character(z) && length(z) == 1L && !is.na(z) && z %in% colnames(m)
  if (!valid_name(feature) || (!is.null(colour_by) && !valid_name(colour_by))) {
    sb_abort("Choose named explained features.", "stormbreaker_bad_visualization_argument")
  }
  needed <- unique(c(feature, colour_by))
  if (!is.data.frame(feature_values) || !is.character(id_col) || length(id_col) != 1L ||
      is.na(id_col) || !all(c(id_col, needed) %in% names(feature_values))) {
    sb_abort("Supply feature values with an explicit ID column and requested features.",
      "stormbreaker_bad_visualization_argument")
  }
  ids <- as.character(feature_values[[id_col]])
  if (anyNA(ids) || anyDuplicated(ids) || length(ids) != nrow(m) || !setequal(ids, rownames(m))) {
    sb_abort("Feature-value IDs must match explained observations exactly.",
      "stormbreaker_bad_visualization_argument")
  }
  values <- feature_values[match(rownames(m), ids), needed, drop = FALSE]
  if (!all(vapply(values, function(v) is.numeric(v) && all(is.finite(v)), logical(1)))) {
    sb_abort("Feature values must be finite numeric columns.", "stormbreaker_bad_visualization_argument")
  }
  d <- data.frame(observation = seq_len(nrow(m)), feature_value = values[[feature]],
    contribution = unname(m[, feature]))
  feature_value <- contribution <- colour_value <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(feature_value, contribution)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey65", linetype = 2)
  if (is.null(colour_by)) {
    p <- p + ggplot2::geom_point(colour = "#2878B5", alpha = .7)
  } else {
    d$colour_value <- values[[colour_by]]
    p$data <- d
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = colour_value), alpha = .7) +
      ggplot2::scale_colour_gradient(low = "#2878B5", high = "#D64E66", name = colour_by)
  }
  p + ggplot2::geom_rug(sides = "b", alpha = .2) +
    ggplot2::labs(title = paste("SHAP:", feature),
      x = feature, y = paste0("SHAP (", sb_viz_score_label(x$score_scale), ")"),
      subtitle = paste(nrow(m), "saved observations | original feature units"),
      caption = "Predictive association, not causality. No clinical cutoff inferred. Colouring is not an interaction estimate. Individual values may be sensitive.") +
    sb_shap_theme()
}

# ---- plot-shap ----

#' @rdname sb_model_shap
#' @export
sb_plot_shap_importance <- function(x, top_n = 20L, groupby = NULL) {
  m <- sb_viz_shap_matrix(x)
  top_n <- sb_viz_positive_integer(top_n, "top_n")
  importance <- colMeans(abs(m))
  chosen <- utils::head(order(-importance, seq_along(importance)), top_n)
  data <- data.frame(
    feature = colnames(m)[chosen], mean_absolute_shap = unname(importance[chosen]),
    n_observations = nrow(m), stringsAsFactors = FALSE
  )
  features <- data$feature
  if (!is.null(groupby)) {
    if (!is.atomic(groupby) || !is.null(dim(groupby)) || length(groupby) != nrow(m) || anyNA(groupby)) {
      sb_abort("groupby must supply one nonmissing group per explained observation.", "stormbreaker_bad_visualization_argument")
    }
    ids <- names(groupby)
    if (!is.null(ids)) {
      if (anyNA(ids) || anyDuplicated(ids) || !setequal(ids, rownames(m))) {
        sb_abort("Named groupby IDs must match explained observations exactly.", "stormbreaker_bad_visualization_argument")
      }
      groupby <- groupby[match(rownames(m), ids)]
    }
    groups <- as.character(groupby)
    if (any(!nzchar(trimws(groups)))) sb_abort("groupby labels cannot be empty.", "stormbreaker_bad_visualization_argument")
    levels <- if (is.factor(groupby)) levels(droplevels(groupby)) else unique(groups)
    data <- do.call(rbind, lapply(levels, function(g) {
      rows <- groups == g
      data.frame(feature = features,
        mean_absolute_shap = unname(colMeans(abs(m[rows, chosen, drop = FALSE]))),
        n_observations = sum(rows), group = g, stringsAsFactors = FALSE)
    }))
    data$group <- factor(data$group, levels = levels)
  }
  data$feature_axis <- factor(data$feature, levels = rev(features))
  data$value_label <- sb_shap_number(data$mean_absolute_shap)
  feature_axis <- mean_absolute_shap <- value_label <- group <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(x = mean_absolute_shap, y = feature_axis))
  if (is.null(groupby)) {
    p <- p + ggplot2::geom_col(fill = "#2878B5", width = 0.62)
  } else {
    p <- p + ggplot2::aes(fill = group, group = group) +
      ggplot2::geom_col(width = 0.72, position = ggplot2::position_dodge(width = 0.8, orientation = "y")) +
      ggplot2::labs(fill = "Group")
  }
  p +
    ggplot2::geom_text(ggplot2::aes(label = value_label),
      position = if (is.null(groupby)) "identity" else ggplot2::position_dodge(width = 0.8, orientation = "y"),
      hjust = -0.2,
      family = "Arial", size = 3.5, colour = "#253746"
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(
      title = "SHAP importance", x = paste0("Mean |SHAP| (", sb_viz_score_label(x$score_scale), ")"), y = NULL
    ) +
    sb_shap_theme()
}

#' @rdname sb_model_shap
#' @export
sb_plot_shap_heatmap <- function(x, top_n = 20L) {
  m <- sb_viz_shap_matrix(x)
  top_n <- sb_viz_positive_integer(top_n, "top_n")
  importance <- colMeans(abs(m))
  chosen <- utils::head(order(-importance, seq_along(importance)), top_n)
  features <- colnames(m)[chosen]
  data <- data.frame(
    observation = rep(seq_len(nrow(m)), length(features)),
    feature = factor(rep(features, each = nrow(m)), levels = rev(features)),
    contribution = as.vector(m[, chosen, drop = FALSE])
  )
  observation <- feature <- contribution <- NULL
  ggplot2::ggplot(data, ggplot2::aes(observation, feature, fill = contribution)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient2(low = "#2878B5", mid = "#FFFFFF", high = "#BE3354",
      name = paste0("SHAP (", sb_viz_score_label(x$score_scale), ")")) +
    ggplot2::labs(title = "SHAP heatmap", x = "Observation", y = NULL) +
    sb_shap_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' @rdname sb_model_shap
#' @export
sb_plot_shap_beeswarm <- function(x, top_n = 20L, feature_values = NULL, id_col = NULL) {
  m <- sb_viz_shap_matrix(x)
  top_n <- sb_viz_positive_integer(top_n, "top_n")
  importance <- colMeans(abs(m))
  chosen <- utils::head(order(-importance, seq_along(importance)), top_n)
  values <- NULL
  if (!is.null(feature_values)) {
    valid <- is.data.frame(feature_values) && is.character(id_col) && length(id_col) == 1L &&
      !is.na(id_col) && id_col %in% names(feature_values) && all(colnames(m) %in% names(feature_values))
    if (!valid) {
      sb_abort(
        "Provide feature_values with an explicit id_col and all explained feature columns.",
        "stormbreaker_bad_visualization_argument"
      )
    }
    ids <- as.character(feature_values[[id_col]])
    if (anyNA(ids) || anyDuplicated(ids) || !setequal(ids, rownames(m)) || length(ids) != nrow(m)) {
      sb_abort("Feature-value IDs must match explained observations exactly.", "stormbreaker_bad_visualization_argument")
    }
    values <- feature_values[match(rownames(m), ids), colnames(m), drop = FALSE]
    if (!all(vapply(values, function(value) is.numeric(value) && all(is.finite(value)), logical(1)))) {
      sb_abort("Explained feature values must be finite numeric columns.", "stormbreaker_bad_visualization_argument")
    }
  } else if (!is.null(id_col)) sb_abort("id_col is only used with feature_values.", "stormbreaker_bad_visualization_argument")
  rows <- lapply(seq_along(chosen), function(i) {
    j <- chosen[[i]]
    contribution <- m[, j]
    span <- diff(range(contribution))
    bins <- if (span == 0) {
      rep(1L, length(contribution))
    } else {
      floor(30 * (contribution - min(contribution)) / span)
    }
    max_bin <- max(tabulate(match(bins, unique(bins))))
    offset <- numeric(length(contribution))
    for (bin in unique(bins)) {
      idx <- which(bins == bin)
      idx <- idx[order(contribution[idx], idx)]
      sequence <- seq_along(idx) - 1L
      offset[idx] <- ceiling(sequence / 2) * ifelse(sequence %% 2 == 1, 1, -1) * 0.7 / max(2, max_bin)
    }
    tab <- data.frame(
      observation = seq_len(nrow(m)), feature = colnames(m)[j],
      contribution = unname(contribution), y = length(chosen) - i + 1 + offset,
      feature_value = NA_real_, colour_value = NA_real_, stringsAsFactors = FALSE
    )
    if (!is.null(values)) {
      v <- values[[j]]
      width <- diff(range(v))
      tab$feature_value <- v
      tab$colour_value <- if (width == 0) rep(0.5, length(v)) else (v - min(v)) / width
    }
    tab
  })
  data <- sb_view_bind(rows)
  contribution <- y <- colour_value <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(x = contribution, y = y)) +
    ggplot2::geom_hline(yintercept = seq_along(chosen), colour = "#E3E8EC", linetype = 3, linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, colour = "#9AA7B0", linewidth = 0.45)
  if (is.null(values)) {
    p <- p + ggplot2::geom_point(colour = "#2878B5", size = 1.9, alpha = 0.8)
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = colour_value), size = 1.9, alpha = 0.8) +
      ggplot2::scale_colour_gradientn(
        colours = c("#2878B5", "#AD92B5", "#D45064"), limits = c(0, 1),
        breaks = c(0, 1), labels = c("Low", "High"), name = "Feature value\n(within gene)",
        guide = ggplot2::guide_colourbar(
          barheight = grid::unit(38, "mm"), barwidth = grid::unit(3, "mm"),
          theme = ggplot2::theme(legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 8)))
        )
      )
  }
  p + ggplot2::scale_y_continuous(breaks = seq_along(chosen), labels = rev(colnames(m)[chosen])) +
    ggplot2::labs(
      title = "SHAP summary", x = paste0("SHAP (", sb_viz_score_label(x$score_scale), ")"), y = NULL,
      subtitle = paste(nrow(m), "observations | one dot per sample per gene | left: lower score; right: higher score"),
      caption = sb_viz_wrap(sb_viz_shap_caption(x), width = 105)
    ) + sb_shap_theme()
}

#' @rdname sb_model_shap
#' @export
sb_plot_shap_waterfall <- function(x, row = 1L, top_n = 10L) {
  m <- sb_viz_shap_matrix(x)
  row <- sb_viz_positive_integer(row, "row")
  top_n <- sb_viz_positive_integer(top_n, "top_n")
  if (row > nrow(m)) sb_abort("row exceeds the number of explained observations.", "stormbreaker_bad_visualization_argument")
  terms <- stats::setNames(as.numeric(m[row, , drop = FALSE]), colnames(m))
  chosen <- utils::head(order(-abs(terms), seq_along(terms)), top_n)
  other <- setdiff(seq_along(terms), chosen)
  contributions <- unname(terms[chosen])
  feature <- names(terms)[chosen]
  kind <- rep("feature", length(chosen))
  if (length(other)) {
    contributions <- c(contributions, sum(terms[other]))
    feature <- c(feature, paste0("Other features (", length(other), ")"))
    kind <- c(kind, "remainder")
  }
  end <- x$reference_value + cumsum(contributions)
  start <- c(x$reference_value, utils::head(end, -1))
  score <- x$reference_value + sum(terms)
  data <- data.frame(
    feature = c("Background reference", feature, "Final model score"),
    kind = c("reference", kind, "prediction"),
    contribution = c(x$reference_value, contributions, score),
    start = c(0, start, 0), end = c(x$reference_value, end, score), stringsAsFactors = FALSE
  )
  data$y <- rev(seq_len(nrow(data)))
  data$xmin <- pmin(data$start, data$end)
  data$xmax <- pmax(data$start, data$end)
  data$direction <- ifelse(data$kind == "reference", "Reference", ifelse(data$kind == "prediction", "Final score",
    ifelse(data$contribution >= 0, "Increase", "Decrease")
  ))
  data$value_label <- sb_shap_number(data$contribution,
    signed = !data$kind %in% c("reference", "prediction")
  )
  data$label <- data$value_label
  extent <- range(c(x$reference_value, start, end, score))
  span <- diff(extent)
  if (span == 0) span <- max(1, abs(score) * 0.1)
  data$label_x <- extent[2] + span * 0.12
  terms_data <- data[!data$kind %in% c("reference", "prediction"), , drop = FALSE]
  endpoints <- data[data$kind %in% c("reference", "prediction"), , drop = FALSE]
  polygons <- lapply(seq_len(nrow(terms_data)), function(i) {
    d <- terms_data[i, ]
    tip <- min(abs(d$end - d$start) * 0.4, span * 0.016)
    neck <- d$end - sign(d$contribution) * tip
    data.frame(
      px = c(d$start, neck, d$end, neck, d$start),
      py = d$y + c(-0.27, -0.27, 0, 0.27, 0.27),
      group = i, direction = d$direction
    )
  })
  polygons <- sb_view_bind(polygons)
  links <- data[-nrow(data), , drop = FALSE]
  links$next_y <- utils::tail(data$y, -1)
  xmin <- xmax <- y <- direction <- end <- next_y <- label <- NULL
  px <- py <- group <- label_x <- value_label <- NULL
  ggplot2::ggplot(data) +
    ggplot2::geom_segment(ggplot2::aes(y = y, yend = y),
      x = extent[1] - span * 0.06, xend = extent[2] + span * 0.07,
      colour = "#E3E8EC", linetype = 3, linewidth = 0.3
    ) +
    ggplot2::geom_polygon(
      data = polygons,
      ggplot2::aes(x = px, y = py, group = group, fill = direction)
    ) +
    ggplot2::geom_point(
      data = endpoints, ggplot2::aes(x = end, y = y),
      shape = 18, size = 3, colour = "#536777"
    ) +
    ggplot2::geom_segment(
      data = terms_data[terms_data$contribution == 0, , drop = FALSE],
      ggplot2::aes(x = end, xend = end, y = y - 0.2, yend = y + 0.2), colour = "#9AA7B0"
    ) +
    ggplot2::geom_segment(
      data = links, ggplot2::aes(x = end, xend = end, y = y - 0.3, yend = next_y + 0.3),
      colour = "#9AA7B0", linetype = 3, linewidth = 0.4
    ) +
    ggplot2::geom_text(ggplot2::aes(x = label_x, y = y, label = value_label),
      hjust = 0, size = 3.5, family = "Arial", colour = "#253746"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(extent[1] - span * 0.06, extent[2] + span * 0.32),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_continuous(breaks = data$y, labels = data$feature) +
    ggplot2::scale_fill_manual(
      values = c(Increase = "#D45064", Decrease = "#2878B5"),
      breaks = c("Increase", "Decrease"), labels = c("Increases score", "Decreases score"), drop = FALSE
    ) +
    ggplot2::labs(
      title = "Individual prediction",
      subtitle = paste("Observation", row, "| contributions ranked by magnitude | numbers: signed contribution or endpoint score"),
      x = paste0("Score (", sb_viz_score_label(x$score_scale), ")"), y = NULL, fill = NULL,
      caption = sb_viz_wrap(paste(sb_viz_shap_caption(x), "Omitted terms are summed, not discarded."), width = 105)
    ) +
    sb_shap_theme() +
    ggplot2::theme(legend.position = "top", legend.justification = "left")
}

sb_shap_number <- function(x, signed = FALSE) {
  label <- vapply(x, function(value) format(signif(value, 3), trim = TRUE), character(1))
  ifelse(rep_len(signed, length(x)) & x >= 0, paste0("+", label), label)
}

sb_shap_theme <- function() {
  ggplot2::theme_classic(base_size = 11, base_family = "Arial") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold", colour = "#172B3A", margin = ggplot2::margin(b = 8)),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "#253746"),
      axis.text.y = ggplot2::element_text(size = 11),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12)),
      axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(colour = "#9AA7B0", linewidth = 0.4),
      axis.ticks.x = ggplot2::element_line(colour = "#9AA7B0", linewidth = 0.4),
      legend.title = ggplot2::element_text(size = 10), legend.text = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(16, 22, 12, 12), plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

sb_viz_shap_matrix <- function(x) {
  sb_require_curve_plotting()
  m <- x$contributions
  valid <- inherits(x, c("stormbreaker_additive_shap", "stormbreaker_tree_shap", "stormbreaker_model_shap")) && identical(x$status, "success") &&
    is.matrix(m) && is.numeric(m) && nrow(m) > 0 && ncol(m) > 0 && all(is.finite(m)) &&
    !is.null(colnames(m)) && !is.null(rownames(m)) && !anyNA(colnames(m)) && !anyNA(rownames(m)) &&
    !anyDuplicated(colnames(m)) && !anyDuplicated(rownames(m)) &&
    is.numeric(x$reference_value) && length(x$reference_value) == 1L && is.finite(x$reference_value) &&
    is.character(x$score_scale) && length(x$score_scale) == 1L && !is.na(x$score_scale)
  if (!valid) sb_abort("Supply a valid successful SHAP result with finite named contributions.", "stormbreaker_bad_shap_visualization")
  if (!is.data.frame(x$scores) || !identical(as.character(x$scores$patient_id), rownames(m)) ||
    !isTRUE(all.equal(unname(x$scores$score), unname(x$reference_value + rowSums(m)), tolerance = 1e-10))) {
    sb_abort("Saved SHAP scores do not reconstruct from the reference and contributions.", "stormbreaker_bad_shap_visualization")
  }
  m
}

sb_viz_shap_caption <- function(x) {
  if (inherits(x, c("stormbreaker_tree_shap", "stormbreaker_model_shap"))) return(paste0(
    "Approximate permutation interventional SHAP; ", x$permutations,
    " orders, background n = ", x$background$n, "; scale = ", x$score_scale, ". Not exact TreeSHAP."))
  paste0(
    "Exact empirical interventional contributions on the saved score scale, not probability effects. Background n = ",
    x$background$n, "; reference = ", format(x$reference_value, digits = 4), ". Descriptive, not causal."
  )
}

# ---- model-coefficients ----

# ---- cox-forest ----

#' @rdname sb_survival_curves
#' @export
sb_cox_forest <- function(task, conf_level = 0.95) {
  task <- sb_validate_task_contract(task)
  if (!is.null(task$group_col)) stop("Use one independent patient per row for this Cox forest.")
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) stop("conf_level must be in (0,1).")
  sb_require_curve_plotting()
  x <- task$data[task$predictors]
  response <- survival::Surv(task$data[[task$time_col]], task$data[[task$event_col]])
  tab <- do.call(rbind, lapply(seq_along(x), function(i) {
    sb_score_univariable_cox_feature(names(x)[i], x, response,
      rep(1, nrow(x)), survival::coxph.control(), i)
  }))
  z <- stats::qnorm((1 + conf_level) / 2)
  tab$hr <- exp(tab$coefficient)
  tab$lower <- exp(tab$coefficient - z * tab$standard_error)
  tab$upper <- exp(tab$coefficient + z * tab$standard_error)
  tab$p_adjusted <- stats::p.adjust(tab$p_value, method = "BH")
  tab$plot_available <- tab$status == "success" & is.finite(tab$hr) &
    is.finite(tab$lower) & is.finite(tab$upper) & tab$lower > 0
  shown <- tab[tab$plot_available, ]
  shown$feature <- factor(shown$feature, levels = rev(tab$feature))
  shown$direction <- ifelse(shown$hr >= 1, "HR > 1", "HR < 1")
  shown$direction[shown$hr == 1] <- "HR = 1"
  feature <- hr <- lower <- upper <- direction <- NULL
  p <- ggplot2::ggplot(shown, ggplot2::aes(hr, feature, colour = direction)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper, yend = feature), linewidth = .6) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c("HR > 1" = "#C87558", "HR < 1" = "#287E8E", "HR = 1" = "grey50")) +
    ggplot2::labs(x = "Hazard ratio", y = NULL) + sb_viz_theme() +
    ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank())
  list(plot = p, table = tab, conf_level = conf_level)
}

# ---- model-coefficients ----

#' @rdname sb_model_coefficients
#' @export
sb_model_coefficients <- function(fit, reference_data = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  state <- fit$learner_state
  supported <- c(
    "logistic", "glmnet_logistic", "lda_binary", "slope_binary", sb_abess_binary_learner_ids(),
    sb_ncvreg_binary_learner_ids(), sb_linear_survival_learner_ids(), sb_extended_logistic_learner_ids(),
    sb_survreg_learner_ids()
  )
  if (!state$learner %in% supported) {
    sb_abort(
      "This learner has no supported linear coefficient view.",
      "stormbreaker_unsupported_coefficients"
    )
  }
  reference_sd <- rep(NA_real_, length(state$features))
  reference_n <- 0L
  if (!is.null(reference_data)) {
    reference <- sb_select_final_feature_frame(
      reference_data, state$features,
      "Reference data", "reference fitted features"
    )
    if (nrow(reference) < 2L) {
      sb_abort(
        "Reference data need at least two rows.",
        "stormbreaker_bad_reference_data"
      )
    }
    reference_sd <- unname(vapply(reference, stats::sd, numeric(1)))
    reference_n <- nrow(reference)
  }
  survival <- state$learner %in% sb_linear_survival_learner_ids()
  aft <- state$learner %in% sb_survreg_learner_ids()
  if (aft) sb_restore_additional_survival_model(state, state$learner)
  parameters <- if (state$learner == "lda_binary") sb_discriminant_parameters(state) else
    if (state$learner %in% sb_extended_logistic_learner_ids()) sb_extended_logistic_parameters(state) else
    if (aft) list(coefficient = -state$coefficient, intercept = -state$intercept) else
    list(coefficient = state$coefficient, intercept = if (survival) NULL else state$intercept)
  structure(list(
    table = data.frame(
      feature = state$features, coefficient = parameters$coefficient,
      reference_sd = reference_sd,
      scaled_coefficient = parameters$coefficient * reference_sd
    ),
    intercept = parameters$intercept,
    score_scale = if (aft) "negative_log_time_location" else
      if (survival) "zero_reference_log_hazard_ratio" else "class_1_log_odds",
    reference_n = reference_n, learner = state$learner,
    pipeline_id = fit$pipeline$id
  ), class = "stormbreaker_coefficients")
}

#' @rdname sb_model_coefficients
#' @export
sb_plot_coefficients <- function(x, scale = c("raw", "reference_sd"),
                                 top_n = NULL, include_zero = FALSE) {
  sb_require_curve_plotting()
  if (!inherits(x, "stormbreaker_coefficients")) {
    sb_abort(
      "Supply an sb_model_coefficients result.", "stormbreaker_bad_coefficients"
    )
  }
  scale <- match.arg(scale)
  if (!is.logical(include_zero) || length(include_zero) != 1L || is.na(include_zero)) {
    sb_abort("include_zero must be TRUE or FALSE.", "stormbreaker_bad_coefficients")
  }
  if (scale == "reference_sd" && x$reference_n < 2L) {
    sb_abort(
      "Supply original-training reference data to use SD scaling.",
      "stormbreaker_bad_reference_data"
    )
  }
  d <- x$table
  d$value <- if (scale == "raw") d$coefficient else d$scaled_coefficient
  if (!include_zero) d <- d[d$coefficient != 0, , drop = FALSE]
  d <- d[order(-abs(d$value), seq_len(nrow(d))), , drop = FALSE]
  if (!is.null(top_n)) d <- utils::head(d, sb_viz_positive_integer(top_n, "top_n"))
  d$feature_axis <- factor(d$feature, levels = rev(d$feature))
  d$direction <- ifelse(d$value >= 0, "Positive", "Negative")
  feature_axis <- value <- direction <- NULL
  intercept <- if (is.null(x$intercept)) {
    "Survival score: no intercept; baseline defines absolute survival."
  } else {
    paste0("Intercept = ", format(signif(x$intercept, 5)), " (", x$score_scale, ").")
  }
  score_label <- sb_viz_score_label(x$score_scale)
  ggplot2::ggplot(d, ggplot2::aes(x = value, y = feature_axis, fill = direction)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#9AA6AF", linewidth = 0.4) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::scale_fill_manual(values = c(Negative = "#24547A", Positive = "#B86236")) +
    ggplot2::labs(
      title = "Coefficients", y = NULL, fill = NULL,
      subtitle = paste(nrow(d), "displayed /", nrow(x$table), "fitted terms |", x$learner),
      x = paste0("Coefficient (", score_label,
        if (scale == "raw") "/unit)" else "/SD)"),
      caption = sb_viz_wrap(paste(
        intercept,
        "Conditional associations, not causal effects or SHAP importance.",
        if (scale == "reference_sd") {
          paste("Reference:", x$reference_n, "caller-supplied training rows.")
        } else {
          "Input units are those of the supplied model features."
        }
      ))
    ) +
    sb_viz_theme()
}
