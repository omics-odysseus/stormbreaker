#' Feature analysis and selection summaries
#' @name sb_analyze_methods
#' @description Analyze feature-selection methods and summarize saved selections, rankings and overlaps.
#' @rawRd \arguments{
#' \item{task}{
#' A \code{\link[=sb_binary_task]{sb_binary_task()}}, \code{\link[=sb_multiclass_task]{sb_multiclass_task()}} or \code{\link[=sb_surv_task]{sb_surv_task()}} containing the data,
#' outcome and candidate features. Each method receives all candidate features.
#' }
#' \item{methods}{
#' Explicit names from \code{\link[=sb_method_catalog]{sb_method_catalog()}}. Names are
#' case-insensitive; spaces and hyphens may replace underscores.
#' }
#' \item{inner_folds}{
#' Inner-CV settings for LASSO and Elastic Net.
#' Other engines retain their existing documented parameters and seeds.
#' }
#' \item{inner_seed}{
#' Inner-CV settings for LASSO and Elastic Net.
#' Other engines retain their existing documented parameters and seeds.
#' }
#' \item{elastic_net_alpha}{
#' Elastic Net mixing fraction, strictly between 0 and 1.
#' }
#' \item{lambda_rule}{
#' Glmnet lambda rule: lambda.1se or lambda.min.
#' }
#' \item{task_type}{
#' Binary, multiclass classification or right-censored survival.
#' }
#' \item{x}{
#' \strong{For \code{sb_method_features()}:} An \code{\link[=sb_analyze_methods]{sb_analyze_methods()}} result.
#' 
#' \strong{For \code{sb_gene_selection_ranking()}:} A saved CV result, CV comparison, repeated CV comparison, or a
#' cohort comparison (only its development CV is used).
#' 
#' \strong{For \code{sb_plot_gene_overlap()}:} Result of sb_gene_overlap.
#' 
#' \strong{For \code{sb_plot_gene_selection_ranking()}:} A \code{\link[=sb_gene_selection_ranking]{sb_gene_selection_ranking()}} object.
#' 
#' \strong{For \code{sb_plot_selection_frequency()}:} A CV result/comparison or supported composite CV result; a uniquely
#' named list may identify result sets. External evaluations have no training
#' fold selection records and are not accepted.
#' }
#' \item{method}{
#' One method present in the result.
#' }
#' \item{class}{
#' Optional multiclass label. NULL returns overall features;
#' non-NULL returns the saved class-specific coefficients or importance.
#' }
#' \item{sets}{
#' Named list of character gene vectors. Duplicate genes within a
#' set are counted once; empty sets are retained. Identifiers must already
#' use the same naming system. No gene annotation or model fitting occurs.
#' }
#' \item{source}{
#' \strong{For \code{sb_gene_selection_ranking()}:} Learned selector counts or sparse glmnet nonzero counts.
#' 
#' \strong{For \code{sb_plot_selection_frequency()}:} Which saved sets to show: \code{"selector"} or
#' \code{"glmnet_active_set"}. These are different notions of selection.
#' }
#' \item{candidate_features}{
#' Optional complete input gene list. Evaluated genes
#' with no selections receive zero; genes absent from the saved candidate
#' sets receive NA, not zero. Default: saved task predictors in input order.
#' }
#' \item{cohort_label}{
#' \strong{For \code{sb_gene_selection_ranking()}:} Optional name for this development cohort.
#' 
#' \strong{For \code{sb_plot_selection_frequency()}:} Optional label for a single object.
#' }
#' \item{type}{
#' auto uses Venn for two to four sets, otherwise bars.
#' }
#' \item{top_n}{
#' \strong{For \code{sb_plot_gene_selection_ranking()}:} Maximum genes to display; the ranking table retains all genes.
#' 
#' \strong{For \code{sb_plot_selection_frequency()}:} Maximum features shown in each frequency panel. Selection for
#' display uses decreasing saved frequency, with input-order ties.
#' }
#' \item{max_intersections}{
#' Maximum exact intersections displayed per UpSet
#' panel, ordered by decreasing size with stable ties. All intersections
#' remain in the plot attribute. Bars share a scale across panels and their
#' actual feature counts are labelled.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_analyze_methods()}}{
#' Independently analyze named feature-selection methods
#' 
#' \strong{Returns:} A \code{stormbreaker_method_analysis} with named results, a method status
#' table and candidate features. Each successful result retains its fit,
#' feature table, native evidence, parameters, diagnostics, ggplots and separate
#' explanatory notes. Fit and plot statuses are reported separately. The
#' method status table records failures while the other methods run. Results
#' are returned in memory.
#' 
#' Each method uses all supplied samples. CV curves describe tuning and OOB
#' curves describe the fitted forest. To evaluate feature selection on
#' held-out data, use pipeline comparison so selection runs inside each
#' training fold. Saved fits can contain patient-derived data. The survival
#' random forest adapter reports ranger impurity rankings rather than OOB curves.
#' }
#' 
#' \subsection{\code{sb_method_catalog()}}{
#' List standalone feature-analysis methods
#' 
#' \strong{Returns:} A data frame of method names, engines and feature-result semantics.
#' Unsupported methods remain visible with supported = FALSE. Availability
#' of optional dependencies is checked when the requested method executes.
#' }
#' 
#' \subsection{\code{sb_method_features()}}{
#' Retrieve saved gene selection or ranking without refitting
#' 
#' \strong{Returns:} \code{sb_method_features()} returns the complete feature table, including
#' zero coefficients and unselected genes. For random forest, selected is NA
#' and rank is decreasing native importance (ties share a rank).
#' \code{sb_selected_features()} returns a character vector, possibly empty, for
#' methods with a defined selection rule. It errors for ranking-only methods
#' rather than silently treating all or top-ranked genes as selected.
#' Multiclass LDA, QDA and Gaussian naive Bayes return model features and
#' optional class means/SD, not selected genes. A decision tree returns native
#' importance and a separate used_for_split flag for actual primary splits.
#' }
#' 
#' \subsection{\code{sb_gene_overlap()}}{
#' Summarize exact intersections of named gene sets
#' 
#' \strong{Returns:} List with sets, membership, exclusive intersections (including gene
#' members), and common genes present in every supplied set.
#' }
#' 
#' \subsection{\code{sb_gene_selection_ranking()}}{
#' Rank candidate genes by saved cross-validation selection counts
#' 
#' \strong{Returns:} An object with \code{ranking}, \code{by_pipeline}, \code{groups} and \code{status}.
#' Ranking columns include selected_count, eligible_count, selection_frequency
#' and selected_methods. All supplied genes remain present. Tied counts share
#' a rank and retain input order. No model is fitted or selected.
#' 
#' \strong{Details:} Each successful outer selection fold is one opportunity. Copies of
#' the same selector configuration on the same saved split assignment are
#' counted once, even if several downstream learners use them. Conflicting
#' saved counts for the same configuration are marked unavailable. Missing
#' split identity prevents deduplication and is recorded in groups.
#' A failed learner does not erase a separately successful selector record.
#' None/fixed feature panels are not learned selections. Glmnet active sets
#' require alpha > 0 and remain separate from selector counts. Repeated CV
#' repetitions remain distinct counting units. External outcomes are not read.
#' One call covers one development cohort and task; do not pool cohorts.
#' }
#' 
#' \subsection{\code{sb_plot_gene_overlap()}}{
#' Plot gene-set overlap as a Venn diagram or intersection bars
#' 
#' \strong{Returns:} A plot with exact membership in its gene_overlap attribute. Venn
#' regions are schematic, not proportional to set sizes. Bars show exclusive
#' intersections: genes in the named sets and absent from every other set.
#' }
#' 
#' \subsection{\code{sb_plot_gene_selection_ranking()}}{
#' Plot ranked gene selection counts
#' 
#' \strong{Returns:} A ggplot with counts and successful-fold denominators.
#' }
#' 
#' \subsection{\code{sb_plot_selection_frequency()}}{
#' Visualize saved feature-selection stability
#' 
#' \strong{Returns:} A ggplot object with saved frequency or Jaccard data and a
#' \code{stormbreaker_selection_status} attribute. The UpSet variant also attaches
#' all exact intersections as \code{stormbreaker_intersections}, including features
#' omitted from the displayed top intersections.
#' 
#' \strong{Details:} Frequencies use the saved successful-fold denominator, not all
#' requested folds. Incomplete selections are explicitly labelled and do not
#' become complete merely because the plot draws. Jaccard shows only saved
#' unique fold pairs; no diagonal, missing pair or zero is invented.
#' A failed model fit does not erase separately successful selector records.
#' Repetitions, cohorts and models remain separate panels. These are
#' descriptive stability displays, not biological importance or causal
#' evidence. An unscreened or fixed feature panel is not reclassified as a
#' learned feature-selection result.
#' }
#' }
#' @rawRd \examples{
#' # sb_analyze_methods
#' \dontrun{
#' task <- sb_binary_task(data, "patient_id", "outcome", predictors = genes)
#' analysis <- sb_analyze_methods(task, c("LASSO", "Boruta", "random_forest"))
#' sb_selected_features(analysis, "LASSO")
#' analysis$results$lasso$plots$coefficient_path
#' sb_method_features(analysis, "random_forest")
#' }
#' }
NULL

# ---- method-analysis ----

#' @rdname sb_analyze_methods
#' @export
sb_method_catalog <- function(task_type = c("binary", "survival", "multiclass")) {
  task_type <- match.arg(task_type)
  if (task_type == "multiclass") {
    catalog <- sb_method_catalog("binary")
    catalog$supported <- catalog$method %in% c("lasso", "elastic_net", "random_forest",
      "boruta", "varselrf", "vsurf", "abess")
    catalog$engine <- NA_character_
    catalog$engine[catalog$method %in% c("lasso", "elastic_net")] <- "glmnet_multinomial"
    catalog$engine[catalog$method == "random_forest"] <- "randomForest_multiclass"
    catalog$engine[catalog$method == "boruta"] <- "boruta_multiclass_roughfix"
    catalog$engine[catalog$method == "varselrf"] <- "varselrf_multiclass_oob_1se"
    catalog$engine[catalog$method == "vsurf"] <- "vsurf_multiclass_prediction"
    catalog$engine[catalog$method == "abess"] <- "abess_multinomial_gic"
    catalog$selection_rule[!catalog$supported] <- "Not implemented for multiclass"
    catalog$selection_rule[catalog$method %in% c("lasso", "elastic_net")] <-
      "Union of class-specific nonzero coefficients at selected lambda"
    catalog$selection_rule[catalog$method == "abess"] <-
      "Union of class-specific nonzero coefficients at GIC-selected support size"
    ridge <- data.frame(method = "ridge", label = "Ridge", supported = TRUE,
      engine = "glmnet_multinomial", feature_result = "model_features",
      selection_rule = "Ridge shrinks coefficients; it does not select genes")
    return(rbind(catalog, sb_multiclass_classical_catalog(), ridge,
      sb_multiclass_native_catalog(), sb_multiclass_information_catalog(), sb_multiclass_omics_catalog()))
  }
  binary <- task_type == "binary"
  data.frame(
    method = c("boruta", "varselrf", "vsurf", "lasso", "elastic_net",
      "mcp", "scad", "abess", "random_forest"),
    label = c("Boruta", "varSelRF", "VSURF", "LASSO", "Elastic Net",
      "MCP", "SCAD", "ABESS", "Random forest"),
    supported = c(TRUE, binary, binary, rep(TRUE, 6)),
    engine = c(paste0("boruta_", task_type, "_roughfix"),
      "varselrf_binary_oob_1se", "vsurf_binary_prediction",
      rep(if (binary) "glmnet_logistic" else "glmnet_cox", 2),
      paste0("ncvreg_", c("mcp", "scad"), if (binary) "_binary_bic" else "_cox_bic"),
      if (binary) "abess_binary_gic" else "abess_cox_gic",
      if (binary) "random_forest_binary" else "ranger_survival"),
    feature_result = c(rep("selection", 8), "ranking"),
    selection_rule = c("Boruta confirmed after TentativeRoughFix",
      "Smallest subset within one SD of minimum OOB error",
      "VSURF prediction stage", rep("Nonzero coefficients at selected lambda", 2),
      rep("Nonzero coefficients at minimum BIC", 2),
      "Nonzero coefficients in GIC-selected subset", "Not defined")
  )
}

sb_method_names <- function(methods) {
  if (!is.character(methods) || !length(methods) || anyNA(methods)) {
    sb_abort("Supply one or more method names.", "stormbreaker_bad_method")
  }
  gsub("[ -]+", "_", tolower(trimws(methods)))
}

sb_method_plot <- function(plot, type, label) {
  titles <- c(boruta_summary = "Importance summary", boruta_hits = "Feature hits",
    boruta_history = "Importance history", selection_history = "Feature elimination",
    vsurf_importance = "Importance", vsurf_stages = "Selection stages",
    coefficient_path = "Coefficient path", cv_curve = "Cross-validation",
    nonzero_path = "Nonzero coefficients", bic_path = "BIC",
    support_criterion = "Subset selection", importance = "Importance",
    oob_error = "OOB error", tree_size = "Tree size")
  shorten <- function(label) {
    if (is.null(label)) return(NULL)
    label <- sub("^(Native|Original) ", "", label)
    sub(" \\(lower is better\\)$", "", label)
  }
  plot + ggplot2::labs(title = paste(label, titles[[type]], sep = ": "),
    x = shorten(plot$labels$x), y = shorten(plot$labels$y),
    subtitle = NULL, caption = NULL)
}

sb_analyze_one_method <- function(task, spec, inner_folds, inner_seed,
                                  elastic_net_alpha, lambda_rule) {
  if (inherits(task, "stormbreaker_multiclass_task")) return(
    sb_analyze_multiclass_method(task, spec, inner_folds, inner_seed, elastic_net_alpha, lambda_rule))
  method <- spec$method
  if (!spec$supported) sb_abort(
    paste(spec$label, "is not implemented for this task type."),
    "stormbreaker_unsupported_method_task")
  selector <- method %in% c("boruta", "varselrf", "vsurf")
  if (selector) {
    x <- task$data[task$predictors]
    binary <- inherits(task, "stormbreaker_binary_task")
    endpoint <- if (binary) list(outcome = as.integer(task$data[[task$outcome_col]])) else
      list(time = task$data[[task$time_col]], event = as.integer(task$data[[task$event_col]]))
    fit <- switch(method,
      boruta = sb_fit_boruta_selector(x, endpoint, spec$engine, capture_history = TRUE),
      varselrf = sb_fit_varselrf_selector(x, endpoint$outcome, spec$engine),
      vsurf = sb_fit_vsurf_binary_selector(x, endpoint$outcome))
    features <- data.frame(feature = task$predictors,
      selected = task$predictors %in% fit$selected_features)
    evidence <- fit$scores
    parameters <- fit$params
    types <- switch(method,
      boruta = c("boruta_summary", "boruta_hits", "boruta_history"),
      varselrf = "selection_history", vsurf = c("vsurf_importance", "vsurf_stages"))
    extract <- function(type) sb_selector_diagnostics(fit, type)
    draw <- sb_plot_selector_diagnostic
  } else {
    args <- list(selector = "none", learner = spec$engine)
    if (method %in% c("lasso", "elastic_net")) args <- c(args, list(
      inner_folds = inner_folds, inner_seed = inner_seed,
      glmnet_alpha = if (method == "lasso") 1 else elastic_net_alpha,
      glmnet_lambda_rule = lambda_rule))
    pipeline <- do.call(sb_pipeline, args)
    fit <- sb_fit_final(task, pipeline, capture_diagnostics = TRUE)
    if (fit$status != "success") sb_abort(fit$error_message, fit$error_class)
    parameters <- list(pipeline = pipeline, tuning = fit$tuning)
    if (method == "random_forest") {
      evidence <- sb_forest_diagnostics(fit, "importance")$table
      features <- evidence
      features$rank <- rank(-features$importance, ties.method = "min")
      features$selected <- NA
      types <- if (inherits(task, "stormbreaker_binary_task"))
        c("importance", "oob_error", "tree_size") else c("importance", "tree_size")
      extract <- function(type) sb_forest_diagnostics(fit, type)
    } else {
      evidence <- sb_model_coefficients(fit)$table
      features <- evidence[c("feature", "coefficient")]
      features$selected <- features$coefficient != 0
      types <- switch(method, lasso = , elastic_net =
        c("coefficient_path", "cv_curve", "nonzero_path"),
        mcp = , scad = c("coefficient_path", "bic_path"),
        abess = "support_criterion")
      extract <- function(type) sb_algorithm_diagnostics(fit, type)
    }
    draw <- sb_plot_algorithm_diagnostic
  }
  diagnostics <- plots <- notes <- list()
  status <- lapply(types, function(type) {
    value <- tryCatch({
      diagnostic <- extract(type)
      diagnostics[[type]] <<- diagnostic
      plot <- draw(diagnostic)
      notes[[type]] <<- list(labels = plot$labels, source = diagnostic$source)
      plots[[type]] <<- sb_method_plot(plot, type, spec$label)
      NULL
    }, error = identity)
    data.frame(view = type, status = if (is.null(value)) "available" else "error",
      reason = if (is.null(value)) "" else conditionMessage(value))
  })
  list(status = "success", method = method, label = spec$label,
    feature_result = spec$feature_result, selection_rule = spec$selection_rule,
    features = features, evidence = evidence, parameters = parameters,
    fit = fit, diagnostics = diagnostics, plots = plots, notes = notes,
    plot_status = do.call(rbind, status))
}

#' @rdname sb_analyze_methods
#' @export
sb_analyze_methods <- function(task, methods, inner_folds = 5L,
                               inner_seed = 1701L, elastic_net_alpha = 0.5,
                               lambda_rule = "lambda.1se") {
  binary <- inherits(task, "stormbreaker_binary_task")
  multiclass <- inherits(task, "stormbreaker_multiclass_task")
  task <- if (binary) sb_validate_binary_task_contract(task) else
    if (multiclass) sb_validate_multiclass_task(task) else sb_validate_task_contract(task)
  task_type <- if (binary) "binary" else if (multiclass) "multiclass" else "survival"
  catalog <- sb_method_catalog(task_type)
  methods <- sb_method_names(methods)
  if (anyDuplicated(methods) || any(!methods %in% catalog$method)) sb_abort(
    "Use distinct method names from sb_method_catalog().", "stormbreaker_bad_method")
  if ("elastic_net" %in% methods && (!is.numeric(elastic_net_alpha) ||
      length(elastic_net_alpha) != 1L || is.na(elastic_net_alpha) ||
      elastic_net_alpha <= 0 || elastic_net_alpha >= 1)) sb_abort(
    "elastic_net_alpha must be strictly between 0 and 1.", "stormbreaker_bad_method")
  results <- lapply(methods, function(method) {
    spec <- catalog[catalog$method == method, ]
    tryCatch(sb_analyze_one_method(task, spec, inner_folds, inner_seed,
      elastic_net_alpha, lambda_rule), error = function(e) list(
        status = "failed", method = method, label = spec$label,
        feature_result = spec$feature_result, selection_rule = spec$selection_rule,
        error_class = class(e)[1L], error_message = conditionMessage(e)))
  })
  names(results) <- methods
  status <- do.call(rbind, lapply(results, function(result) data.frame(
    method = result$method, status = result$status,
    feature_result = result$feature_result,
    n_plots = length(result$plots),
    n_plot_errors = sum(result$plot_status$status == "error"),
    reason = if (result$status == "failed") result$error_message else "")))
  rownames(status) <- NULL
  structure(list(results = results, status = status,
    candidate_features = task$predictors, task_type = task_type,
    interpretation = "Independent method analysis using all supplied samples."),
    class = "stormbreaker_method_analysis")
}

sb_method_result <- function(x, method) {
  method <- sb_method_names(method)
  if (!inherits(x, "stormbreaker_method_analysis") || length(method) != 1L ||
      !method %in% names(x$results)) sb_abort(
    "Choose one method present in this analysis.", "stormbreaker_bad_method")
  result <- x$results[[method]]
  if (result$status != "success") sb_abort(result$error_message, result$error_class)
  result
}

#' @rdname sb_analyze_methods
#' @export
sb_method_features <- function(x, method, class = NULL) {
  result <- sb_method_result(x, method)
  if (is.null(class)) return(result$features)
  if (!is.character(class) || length(class) != 1L || is.na(class) ||
      is.null(result$class_features) || !class %in% result$class_features$class) sb_abort(
    "Choose a class present in this multiclass method result.", "stormbreaker_bad_multiclass_outcome")
  result$class_features[result$class_features$class == class, , drop = FALSE]
}

#' @rdname sb_analyze_methods
#' @export
sb_selected_features <- function(x, method, class = NULL) {
  result <- sb_method_result(x, method)
  if (result$feature_result != "selection") sb_abort(
    "This method has no gene-selection rule. Use sb_method_features() for saved evidence.",
    "stormbreaker_selection_not_defined")
  features <- sb_method_features(x, method, class)
  features$feature[features$selected]
}
