# Source this file to define the examples; no model is trained on loading.
# List names with advanced_examples(); run one with run_advanced_example(name).
# Each call returns its local objects. Choose a fresh working/output directory
# for examples that explicitly export files. Synthetic data regeneration is opt-in.

# ---- regenerate-demo-data ----

# Synthetic inputs only. Source this file, then call regenerate_demo_data("new-directory").
regenerate_demo_data <- function(directory) {
  if (file.exists(directory)) stop("Choose a new directory; existing files are not overwritten.")
  specs <- data.frame(
    file = c("demo-binary-development.csv", "demo-binary-external-a.csv",
      "demo-binary-external-b.csv", "demo-survival-development.csv",
      "demo-survival-external.csv"),
    task = c("binary", "binary", "binary", "survival", "survival"),
    n = c(120, 80, 70, 160, 100), seed = c(42, 17, 23, 42, 17),
    prefix = c("DEV_", "EXTA_", "EXTB_", "SDEV_", "SEXT_")
  )
  if (!dir.create(directory, recursive = TRUE)) stop("Could not create the requested directory.")
  for (i in seq_len(nrow(specs))) {
    d <- if (specs$task[i] == "binary") {
      stormbreaker::bi_example(specs$n[i], specs$seed[i])
    } else {
      stormbreaker::surv_example(specs$n[i], specs$seed[i])
    }
    d$patient_id <- paste0(specs$prefix[i], d$patient_id)
    utils::write.csv(d, file.path(directory, specs$file[i]),
      row.names = FALSE, fileEncoding = "UTF-8")
  }
  invisible(specs)
}

# ---- user-workflow ----

# Source this file to define helpers; sourcing does not fit or select a model.
# All generated templates/demo results are SYNTHETIC, not clinical evidence.
# Requires an installed stormbreaker package. See guides/user-guide.md.

workflow_new_directory <- function(path) {
  # Output files can contain patient data; never replace an existing run.
  if (file.exists(path)) stop("Choose a new output directory: ", path)
  if (!dir.create(path, recursive = TRUE)) stop("Cannot create: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

workflow_read_csv <- function(path) {
  # Preserve IDs such as 00123 and original feature column names.
  header <- names(utils::read.csv(path, nrows = 0, check.names = FALSE))
  classes <- rep(NA_character_, length(header))
  classes[header == "patient_id"] <- "character"
  utils::read.csv(path,
    colClasses = classes, check.names = FALSE,
    stringsAsFactors = FALSE, fileEncoding = "UTF-8"
  )
}

workflow_make_templates <- function(directory) {
  directory <- workflow_new_directory(directory)
  for (kind in c("binary", "survival")) {
    generator <- if (kind == "binary") {
      stormbreaker::bi_example
    } else {
      stormbreaker::surv_example
    }
    development <- generator(n = 100, seed = 20260819)
    external <- generator(n = 40, seed = 20260821)
    development$patient_id <- sprintf("%s_dev_%03d", kind, seq_len(100))
    external$patient_id <- sprintf("%s_ext_%03d", kind, seq_len(40))
    outcomes <- if (kind == "binary") "outcome" else c("time", "event")
    tables <- list(
      development = development,
      external_predictors = external[c("patient_id", "gene_a", "gene_b", "gene_c")],
      external_outcomes = external[c("patient_id", outcomes)]
    )
    for (role in names(tables)) {
      utils::write.csv(tables[[role]],
        file.path(directory, paste0(kind, "_", role, ".csv")),
        row.names = FALSE
      )
    }
  }
  invisible(directory)
}

workflow_develop <- function(development_csv, output_dir,
                             kind = c("binary", "survival"),
                             predictors = c("gene_a", "gene_b", "gene_c"),
                             time_unit = "month", prediction_horizon = 6,
                             folds = 5L, seed = 5201314L,
                             near_tie = 0.01) {
  kind <- match.arg(kind)
  if (!is.numeric(near_tie) || length(near_tie) != 1L ||
    !is.finite(near_tie) || near_tie < 0) {
    stop("near_tie must be nonnegative.")
  }
  directory <- workflow_new_directory(output_dir)
  data <- workflow_read_csv(development_csv)
  task <- if (kind == "binary") {
    stormbreaker::bi_task(
      data,
      "patient_id", "outcome", predictors
    )
  } else {
    stormbreaker::surv_task(data,
      "patient_id", "time", "event", predictors, time_unit,
      prediction_horizon = prediction_horizon
    )
  }
  # Declare BOTH candidates and the decision rule before running CV.
  # Here 'simpler' means the predeclared modelling procedure, not fitted size.
  pipelines <- if (kind == "binary") {
    list(
      logistic = stormbreaker::sb_pipeline(learner = "logistic"),
      naive_bayes = stormbreaker::sb_pipeline(learner = "gaussian_naive_bayes")
    )
  } else {
    list(
      all_features_cox = stormbreaker::sb_pipeline(learner = "coxph"),
      top2_cox = stormbreaker::sb_pipeline(
        learner = "coxph",
        selector = "univ_cox_top_k", selector_n_features = 2L
      )
    )
  }
  metric <- if (kind == "binary") "roc_auc" else "harrell_c"
  rule <- list(
    primary_metric = metric, direction = "higher",
    near_tie = near_tie, simplicity_order = names(pipelines),
    stability_review = paste(
      "Review fold failures, selection stability,",
      "calibration and secondary losses before confirming.",
      "No universal Jaccard cutoff; not_applicable is not failure."
    )
  )
  saveRDS(rule, file.path(directory, "selection_rule_before_cv.rds"))
  compared <- stormbreaker::compare_cv(task, pipelines, folds, seed)
  tab <- compared$comparison
  eligible <- tab$status == "success" & tab$complete_oof &
    tab[[paste0(metric, "_status")]] == "success" & is.finite(tab[[metric]])
  eligible[is.na(eligible)] <- FALSE
  candidate_ids <- vapply(pipelines, function(x) x$id, character(1))
  tab$candidate_name <- names(candidate_ids)[match(tab$pipeline_id, candidate_ids)]
  tab$eligible_for_primary_rule <- eligible
  tab$near_best <- FALSE
  suggested <- NA_character_
  if (any(eligible)) {
    tab$near_best <- eligible & tab[[metric]] >= max(tab[[metric]][eligible]) - near_tie
    tab$near_best[is.na(tab$near_best)] <- FALSE
    suggested <- names(pipelines)[which(names(pipelines) %in%
      tab$candidate_name[tab$near_best])[1L]]
  }
  utils::write.csv(tab, file.path(directory, "candidate_comparison.csv"), row.names = FALSE)
  utils::write.csv(compared$fold_assignment, file.path(directory, "shared_folds.csv"),
    row.names = FALSE
  )
  for (source in c("selector", "glmnet_active_set")) {
    ranking <- stormbreaker::gene_rank(compared, source,
      candidate_features = predictors
    )
    stem <- if (source == "selector") "gene_selection" else "gene_nonzero"
    saveRDS(ranking, file.path(directory, paste0(stem, "_ranking.rds")))
    utils::write.csv(ranking$ranking, file.path(directory, paste0(stem, "_ranking.csv")), row.names = FALSE)
    utils::write.csv(ranking$by_pipeline, file.path(directory, paste0(stem, "_by_pipeline.csv")), row.names = FALSE)
    utils::write.csv(ranking$groups, file.path(directory, paste0(stem, "_counting_groups.csv")), row.names = FALSE)
    if (ranking$status != "not_available" && requireNamespace("ggplot2", quietly = TRUE)) {
      ggplot2::ggsave(file.path(directory, paste0(stem, "_ranking.png")),
        stormbreaker::plot_gene_rank(ranking),
        width = 9, height = 7, dpi = 300
      )
    }
  }
  state <- list(
    kind = kind, task = task, pipelines = pipelines, cv = compared,
    comparison = tab, rule = rule, suggested_candidate = suggested,
    predictors = predictors, time_unit = time_unit,
    prediction_horizon = prediction_horizon
  )
  saveRDS(state, file.path(directory, "development.rds"))
  message(
    "Development only. Review candidate_comparison.csv. Suggested for review: ",
    suggested, ". Nothing has been selected or fitted as a final model."
  )
  invisible(state)
}

workflow_finalize <- function(development_rds, selected_candidate, selection_reason,
                              external_predictors_csv, external_outcomes_csv,
                              output_dir, render_report = TRUE) {
  if (missing(selected_candidate) || missing(selection_reason) ||
    !is.character(selection_reason) || length(selection_reason) != 1L ||
    is.na(selection_reason) || !nzchar(trimws(selection_reason))) {
    stop("Explicit selected_candidate and selection_reason are required.")
  }
  state <- readRDS(development_rds)
  if (!is.character(selected_candidate) || length(selected_candidate) != 1L ||
    is.na(selected_candidate) || !selected_candidate %in% names(state$pipelines)) {
    stop("selected_candidate must name a declared candidate.")
  }
  row <- state$comparison[state$comparison$candidate_name == selected_candidate, ]
  if (!isTRUE(row$eligible_for_primary_rule)) stop("Selected candidate did not complete primary evaluation.")
  directory <- workflow_new_directory(output_dir)
  decision <- list(
    selected_candidate = selected_candidate, reason = selection_reason,
    suggested_candidate = state$suggested_candidate, rule = state$rule,
    chosen_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package_version = as.character(utils::packageVersion("stormbreaker"))
  )
  saveRDS(decision, file.path(directory, "selection_decision.rds"))
  fit <- stormbreaker::fit_final(state$task, state$pipelines[[selected_candidate]])
  if (!identical(fit$status, "success")) stop("Final fit failed; no replacement model.")
  saveRDS(fit, file.path(directory, "final_model.rds"))
  # Reload the artifact to demonstrate that prediction does not require refitting.
  fit <- readRDS(file.path(directory, "final_model.rds"))
  new_data <- workflow_read_csv(external_predictors_csv)
  # Deliberately supply only the identifier and the declared predictors.
  new_data <- new_data[c("patient_id", state$predictors)]
  prediction <- stormbreaker::predict_final(fit, new_data, id_col = "patient_id")
  saveRDS(prediction, file.path(directory, "external_predictions.rds"))
  utils::write.csv(prediction$predictions,
    file.path(directory, "external_predictions.csv"),
    row.names = FALSE
  )
  # External outcomes enter for the first time AFTER the model and predictions exist.
  outcomes <- workflow_read_csv(external_outcomes_csv)
  endpoint <- if (state$kind == "binary") {
    stormbreaker::bi_endpoint(
      outcomes, "patient_id", "outcome"
    )
  } else {
    stormbreaker::surv_endpoint(
      outcomes, "patient_id", "time", "event", state$time_unit,
      prediction_horizon = state$prediction_horizon
    )
  }
  evaluation <- stormbreaker::eval_external(prediction, endpoint)
  saveRDS(evaluation, file.path(directory, "external_evaluation.rds"))
  capture.output(print(evaluation), file = file.path(directory, "evaluation_summary.txt"))
  if (isTRUE(render_report)) {
    stormbreaker::sb_render_report(
      list(development_cv = state$cv, external_validation = evaluation),
      file = file.path(directory, "report.html")
    )
  }
  invisible(list(
    decision = decision, fit = fit, prediction = prediction,
    evaluation = evaluation, directory = directory
  ))
}

workflow_demo <- function(output_dir) {
  directory <- workflow_new_directory(output_dir)
  input <- file.path(directory, "synthetic_inputs")
  workflow_make_templates(input)
  results <- list()
  for (kind in c("binary", "survival")) {
    development_dir <- file.path(directory, paste0(kind, "_development"))
    state <- workflow_develop(file.path(input, paste0(kind, "_development.csv")),
      development_dir,
      kind = kind
    )
    # Synthetic demonstration ONLY: enact the predeclared suggestion.
    # Real projects stop after workflow_develop for human review and confirmation.
    results[[kind]] <- workflow_finalize(file.path(development_dir, "development.rds"),
      selected_candidate = state$suggested_candidate,
      selection_reason = "SYNTHETIC DEMO ONLY: enact predeclared near-tie/simplicity rule.",
      external_predictors_csv = file.path(input, paste0(kind, "_external_predictors.csv")),
      external_outcomes_csv = file.path(input, paste0(kind, "_external_outcomes.csv")),
      output_dir = file.path(directory, paste0(kind, "_final"))
    )
  }
  invisible(results)
}

# ---- visualization-core ----

# STORMBREAKER 0.3.0 visualization tutorial. Synthetic data, not clinical evidence.
# Requires stormbreaker and optional ggplot2 and metafor (for Meta-analysis).
# Source this file, then call stormbreaker_visualization_example().
# Pass a NEW output directory to export all figures and an offline gallery.
stormbreaker_visualization_example <- function(output_dir = NULL, reps = 40L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 to draw the gallery.")
  if (!requireNamespace("metafor", quietly = TRUE)) stop("Install metafor to run the five Meta-analysis examples.")
  library(stormbreaker)
  features <- c("gene_a", "gene_b", "gene_c")
  development <- bi_example(120, seed = 301)
  task <- bi_task(development, "patient_id", "outcome", features)
  pipelines <- list(Logistic = sb_pipeline(learner = "logistic"),
    Top2 = sb_pipeline(selector = "univ_binary_auc_top_k", learner = "logistic",
      selector_n_features = 2L))
  cv <- compare_cv(task, pipelines, folds = 3L, seed = 302)
  # No winning model is chosen from the comparison; both are predeclared.
  fits <- lapply(pipelines, function(p) fit_final(task, p))
  cohorts <- lapply(1:3, function(i) {
    d <- bi_example(100, seed = 310 + i)
    d$patient_id <- paste0("external_", i, "_", seq_len(nrow(d)))
    d
  })
  names(cohorts) <- c("Cohort A", "Cohort B", "Cohort C")
  evaluations <- lapply(cohorts, function(d) {
    predictions <- lapply(fits, function(fit)
      predict_final(fit, d[c("patient_id", features)], "patient_id"))
    compare_external(predictions, bi_endpoint(d, "patient_id", "outcome"))
  })
  bootstraps <- lapply(evaluations, function(e)
    sb_bootstrap_external(e$evaluations[[1]], reps = reps, seed = 320))
  paired <- lapply(evaluations, function(e)
    sb_bootstrap_external_comparison(e, pipelines[[1]]$id, reps = reps, seed = 321))
  meta <- list(
    Discrimination = sb_meta_external_discrimination(bootstraps),
    Calibration_slope = sb_meta_external_calibration_slope(bootstraps),
    Mean_calibration = sb_meta_external_mean_calibration(bootstraps),
    Discrimination_difference = sb_meta_external_discrimination_difference(paired, pipelines[[2]]$id),
    Prediction_error_difference = sb_meta_external_prediction_error_difference(paired, pipelines[[2]]$id)
  )
  explanation <- additive_shap(fits[[1]], cohorts[[1]][1:30, ],
    "patient_id", development[features])
  plots <- list(
    ROC_by_assessment_fold = plot_roc(cv),
    PR_by_assessment_fold = plot_pr(cv),
    External_AUC_heatmap = plot_heatmap(evaluations, "roc_auc"),
    External_AUC_intervals = plot_ext_ci(bootstraps, "roc_auc"),
    Paired_AUC_differences = plot_ext_diff(paired, "roc_auc"),
    Selection_frequency = plot_freq(cv$results[[2]]),
    Selection_Jaccard = plot_jaccard(cv$results[[2]]),
    Selection_intersections = plot_upset(cv$results[[2]]),
    SHAP_importance = plot_shap_imp(explanation),
    SHAP_beeswarm = plot_shap(explanation,
      feature_values = cohorts[[1]][1:30, ], id_col = "patient_id"),
    SHAP_waterfall = plot_shap_one(explanation, top_n = 2L)
  )
  for (name in names(meta)) plots[[paste0("Meta_", name)]] <- plot_forest(meta[[name]])
  # Fixed-horizon survival evidence stays separate from binary evidence.
  survival <- surv_example(120, seed = 330)
  survival_task <- surv_task(survival, "patient_id", "time", "event", features,
    time_unit = "month", prediction_horizon = 6)
  survival_fit <- fit_final(survival_task, sb_pipeline(learner = "coxph"))
  survival_external <- surv_example(100, seed = 331)
  survival_external$patient_id <- paste0("surv_external_", seq_len(nrow(survival_external)))
  survival_evaluation <- eval_external(
    predict_final(survival_fit, survival_external[c("patient_id", features)], "patient_id"),
    surv_endpoint(survival_external, "patient_id", "time", "event",
      time_unit = "month", prediction_horizon = 6))
  survival_bootstrap <- sb_bootstrap_external(survival_evaluation, reps = reps, seed = 332)
  plots$Survival_horizon_heatmap <- plot_heatmap(survival_evaluation, "ipcw_brier")
  # Bootstrap interval keys retain their saved name, ipcw_brier_score.
  plots$Survival_horizon_intervals <- plot_ext_ci(survival_bootstrap, "ipcw_brier_score")
  plots$No_selection_available <- plot_freq(cv$results[[1]])
  # The small number of bootstrap replicates is for a fast software tutorial,
  # not a recommended precision for research. Only three synthetic cohorts:
  # pooled intervals and heterogeneity estimates are highly uncertain.
  if (!is.null(output_dir)) {
    if (file.exists(output_dir)) stop("Choose a new output directory; existing output is never overwritten.")
    if (!dir.create(output_dir, recursive = TRUE)) stop("Cannot create output directory.")
    heights <- rep(7, length(plots))
    heights[1:2] <- 9
    for (i in seq_along(plots)) {
      ggplot2::ggsave(file.path(output_dir, paste0(names(plots)[i], ".png")),
        plots[[i]], width = 12, height = heights[i], dpi = 140, bg = "white")
    }
    sb_render_visual_report(plots, file.path(output_dir, "gallery.html"),
      title = "STORMBREAKER 0.3.0 | Synthetic visualization gallery",
      width = 12, height = heights)
  }
  invisible(list(plots = plots, cv = cv, evaluations = evaluations,
    bootstraps = bootstraps, paired = paired, meta = meta, explanation = explanation,
    survival_evaluation = survival_evaluation, survival_bootstrap = survival_bootstrap))
}

# ---- clinical-binary ----

stormbreaker_clinical_example <- function(output_dir = NULL, figure_format = "png",
                                           font_family = NULL) {
  required <- c("randomForest", "gbm", "xgboost", "lightgbm", "ranger", "ggplot2")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("This six-family example needs: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!is.null(output_dir) && file.exists(output_dir)) {
    stop("Choose a new output directory.", call. = FALSE)
  }
  d <- bi_example(240L, seed = 42L)
  i <- seq_len(nrow(d))
  d$gene_a_proxy <- d$gene_a + sin(i) / 30
  d$gene_b_proxy <- -d$gene_b + cos(i) / 30
  d$gene_d <- sin(i / 3)
  d$gene_e <- cos(i / 5)
  features <- setdiff(names(d), c("patient_id", "outcome"))
  thresholds <- c(.1, .2, .3, .4, .5, .6, .7, .8, .9)
  threshold <- .5
  task <- bi_task(d, "patient_id", "outcome", features, decision_thresholds = thresholds)
  split <- sb_split_binary(task, train_prop = .7, seed = 82L)

  # A small, declared search space. The first four are fixed configurations.
  controls <- list(Logistic = list(learner = "logistic"),
    KNN = list(learner = "knn_binary"),
    GBDT = list(learner = "gbm_binary"),
    XGBoost = list(learner = "xgboost_binary"))
  for (rounds in c(60L, 120L)) for (leaves in c(3L, 7L)) {
    controls[[paste0("LightGBM-", rounds, "-leaves-", leaves)]] <- list(
      learner = "lightgbm_binary", lightgbm_nrounds = rounds, lightgbm_num_leaves = leaves)
  }
  for (trees in c(100L, 200L)) {
    controls[[paste0("ExtraTrees-", trees)]] <- list(
      learner = "ranger_extratrees_binary", ranger_num_trees = trees)
  }
  candidates <- lapply(controls, function(control) {
    do.call(sb_pipeline, c(list(selector = "rf_binary_importance",
      selector_options = list(top_n = 5L, spearman_cutoff = .7, trees = 100L, seed = 4L)),
      control))
  })
  run <- bi_run(split$development, pipelines = candidates, folds = 5L, seed = 11L, plots = FALSE)
  ids <- vapply(candidates, function(x) x$id, character(1))
  cv <- run$comparison$comparison[match(ids, run$comparison$comparison$pipeline_id), , drop = FALSE]
  eligible <- which(cv$status == "success" & cv$roc_auc_status == "success" & is.finite(cv$roc_auc))
  if (!length(eligible)) stop("No candidate has a usable development-CV AUC.", call. = FALSE)
  # which.max uses declared candidate order for ties. Do not pick on assessment data.
  locked <- names(candidates)[eligible[which.max(cv$roc_auc[eligible])]]
  finished <- finish_run(run, locked, external = list(Holdout = split$assessment),
    plots = TRUE, reverse = FALSE)
  evaluation <- finished$evaluations$Holdout
  if (!identical(finished$fit$status, "success") || !identical(evaluation$status, "success")) {
    stop("The locked candidate did not complete holdout evaluation; inspect its saved status.",
      call. = FALSE)
  }
  views <- clinical_plots(evaluation, thresholds = thresholds,
    confusion_threshold = threshold, detail = "full")
  # The interval is for AUC, conditional on the locked model, not a curve band.
  views$plots$ROC <- plot_roc(evaluation, ci_reps = 500L, seed = 19L)
  auc_interval <- attr(views$plots$ROC, "stormbreaker_roc_statistics")
  at_threshold <- sb_threshold_metrics(evaluation, threshold)$table
  scores <- c(roc_auc = evaluation$metrics$roc_auc$estimate,
    average_precision = evaluation$metrics$average_precision$estimate,
    accuracy = at_threshold$accuracy, precision = at_threshold$ppv,
    recall = at_threshold$sensitivity, f1 = at_threshold$f1,
    brier_score = evaluation$metrics$brier_score$estimate)
  metrics <- data.frame(metric = names(scores), estimate = unname(scores))

  selection <- finished$fit$selection$state
  screen <- selection$evidence$rf_screen
  ranking <- selection$scores[match(screen$top_features, selection$scores$feature), ]
  ranking$feature <- factor(ranking$feature, levels = rev(screen$top_features))
  rho <- as.data.frame(as.table(screen$correlation), responseName = "rho")
  names(rho)[1:2] <- c("feature_x", "feature_y")
  distribution <- rbind(
    data.frame(split = split$counts$role, outcome = "0", n = split$counts$n_negative),
    data.frame(split = split$counts$role, outcome = "1", n = split$counts$n_positive))
  exploratory <- list(
    outcome = ggplot2::ggplot(distribution, ggplot2::aes(split, n, fill = outcome)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::scale_fill_manual(values = c("0" = "#497DA1", "1" = "#B34656")) +
      ggplot2::labs(title = "Outcome distribution", x = NULL, y = "Patients", fill = "Outcome") +
      ggplot2::theme_minimal(),
    importance = ggplot2::ggplot(ranking, ggplot2::aes(score, feature, fill = selected)) +
      ggplot2::geom_col(width = .65) +
      ggplot2::scale_fill_manual(values = c("TRUE" = "#497DA1", "FALSE" = "#A7AFB6")) +
      ggplot2::labs(title = "RF importance", x = "Permutation importance", y = NULL,
        fill = "Retained") +
      ggplot2::theme_minimal(),
    correlation = ggplot2::ggplot(rho, ggplot2::aes(feature_x, feature_y, fill = rho)) +
      ggplot2::geom_tile() + ggplot2::coord_equal() +
      ggplot2::scale_fill_gradient2(low = "#2878B5", mid = "#FFFFFF", high = "#BE3354",
        limits = c(-1, 1), name = "Spearman") +
      ggplot2::labs(title = "Development correlations", x = NULL, y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  )
  predictor_columns <- c("patient_id", features)
  # Demonstration subset only. Explanation never uses assessment outcomes.
  query <- utils::head(split$assessment[predictor_columns], 16L)
  shap <- model_shap(finished$fit, query, "patient_id",
    background_data = split$development$data[predictor_columns],
    permutations = 16L, background_n = 12L, seed = 38L)
  finished$explanation <- shap
  finished$plots$external <- NULL
  finished$plots$assessment <- views$plots
  finished$plots$training_data <- exploratory
  finished$plots$explanation <- list(importance = plot_shap_imp(shap),
    beeswarm = plot_shap(shap, feature_values = query, id_col = "patient_id"),
    waterfall = plot_shap_one(shap),
    dependence = plot_shap_dep(shap, colnames(shap$contributions)[1L], query, "patient_id"),
    heatmap = sb_plot_shap_heatmap(shap))
  cv$candidate <- names(candidates)
  cv$learner <- vapply(candidates, function(x) x$learner, character(1))
  if (!is.null(output_dir)) {
    export_run(finished, output_dir, figure_format = figure_format,
      dpi = 300L, font_family = font_family)
    utils::write.csv(split$counts, file.path(output_dir, "split-counts.csv"), row.names = FALSE)
    utils::write.csv(cv, file.path(output_dir, "development-candidates.csv"), row.names = FALSE)
    utils::write.csv(metrics, file.path(output_dir, "holdout-summary.csv"), row.names = FALSE)
    utils::write.csv(auc_interval, file.path(output_dir, "holdout-auc-interval.csv"), row.names = FALSE)
    utils::write.csv(at_threshold, file.path(output_dir, "holdout-threshold.csv"), row.names = FALSE)
    saveRDS(finished$fit, file.path(output_dir, "model.rds"))
  }
  invisible(list(split = split, candidates = cv, selected = locked, finished = finished,
    metrics = metrics, threshold = threshold, auc_interval = auc_interval,
    views = views$status))
}

.stormbreaker_advanced_examples <- list(
  "clinical-binary" = function() stormbreaker_clinical_example(),
  "doctor-compare" = function() {
data <- bi_example(n = 100, seed = 42)
genes <- c("gene_a", "gene_b", "gene_c")
task <- bi_task(data, "patient_id", "outcome", predictors = genes)


run <- bi_run(task, models = c("Logistic", "Lasso", "RF"), folds = 3, plots = FALSE)
comparison <- run$comparison
comparison$comparison[c("status", "roc_auc", "brier_score")]

views <- clinical_plots(comparison)
print(views$plots$ROC)


route <- bi_pipeline(c("Lasso", "RF", "GBM"))
# bi_run(task, pipelines = list(Lasso_RF_GBM = route), folds = 3)
# bi_run(task, models = c("Lasso", "RF", "GBM"), stages = c(1, 2, 3))
# finished <- finish_run(run, "RF", plots = FALSE)
# finished$genes
    invisible(as.list(environment(), all.names = FALSE))
  },
  "doctor-select" = function() {
# A predeclared Lasso-selection + Logistic route; no automatic winner selection.
data <- bi_example(n = 100, seed = 42)
genes <- c("gene_a", "gene_b", "gene_c")
r <- bi_run(data, outcome = "outcome", id = "patient_id", predictors = genes,
            pipelines = list("Lasso + Logistic" = bi_pipeline(c("Lasso", "Logistic"))),
            folds = 3, plots = FALSE)
finished <- finish_run(r, "Lasso + Logistic", plots = FALSE)
finished$genes
    invisible(as.list(environment(), all.names = FALSE))
  },
  "doctor-validate" = function() {
development <- bi_example(n = 100, seed = 42)
genes <- c("gene_a", "gene_b", "gene_c")
task <- bi_task(development, "patient_id", "outcome", predictors = genes)
fit <- fit_final(task, bi_pipeline("Logistic"))


external <- bi_example(n = 80, seed = 43)
external$patient_id <- paste0("external_", external$patient_id)
prediction <- predict_final(fit, external[c("patient_id", genes)], "patient_id")
endpoint <- bi_endpoint(external, "patient_id", "outcome")
evaluation <- eval_external(prediction, endpoint)
evaluation$metrics
    invisible(as.list(environment(), all.names = FALSE))
  },
  "method-analysis" = function() {
# Runnable synthetic example; no files are exported automatically.
set.seed(881)
n <- 120L
d <- data.frame(patient_id = paste0("s", seq_len(n)))
genes <- paste0("gene_", letters[1:8])
for (gene in genes) d[[gene]] <- rnorm(n)
d$outcome <- rbinom(n, 1, plogis(2 * d$gene_a - 1.3 * d$gene_b))
task <- bi_task(d, "patient_id", "outcome", genes)

# Each requested method is fitted independently, once.
analysis <- analyze_methods(task, c("LASSO", "Boruta", "random_forest"))
analysis$status
selected_genes(analysis, "LASSO")
sb_method_features(analysis, "random_forest")
analysis$results$lasso$plots$coefficient_path

# Optional: combine saved plots. No method is fitted again.
if (requireNamespace("patchwork", quietly = TRUE)) {
  panel <- patchwork::wrap_plots(list(
    analysis$results$random_forest$plots$oob_error,
    analysis$results$random_forest$plots$importance,
    analysis$results$boruta$plots$boruta_summary,
    analysis$results$boruta$plots$boruta_history,
    analysis$results$lasso$plots$coefficient_path,
    analysis$results$lasso$plots$cv_curve
  ), ncol = 2) + patchwork::plot_annotation(tag_levels = "A")
  print(panel)
}
    invisible(as.list(environment(), all.names = FALSE))
  },
  "minimum-core" = function() {
# STORMBREAKER minimum end-to-end example
# Synthetic data only: no result below is clinical evidence.
# ---- Right-censored survival workflow --------------------------------------

survival_development <- surv_example(n = 100, seed = 20260812)
survival_task <- surv_task(
  data = survival_development,
  id_col = "patient_id",
  time_col = "time",
  event_col = "event",
  predictors = c("gene_a", "gene_b", "gene_c"),
  time_unit = "month",
  prediction_horizon = 6
)

# Both specifications are declared before resampling and use identical folds.
survival_pipelines <- list(
  all_features_cox = sb_pipeline(
    selector = "none",
    learner = "coxph"
  ),
  top2_cox = sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "coxph",
    selector_n_features = 2L
  )
)
survival_cv <- compare_cv(
  survival_task,
  survival_pipelines,
  folds = 5L,
  seed = 5201314L
)
print(survival_cv)
survival_cv$comparison

# STORMBREAKER does not choose a winner. The study team supplies the final
# specification; this synthetic example uses the predeclared all-feature Cox.
survival_final_pipeline <- survival_pipelines[["all_features_cox"]]
survival_fit <- fit_final(survival_task, survival_final_pipeline)

survival_external <- surv_example(n = 40, seed = 20260820)
survival_external$patient_id <- paste0(
  "survival_external_",
  seq_len(nrow(survival_external))
)
survival_new_data <- survival_external[
  c("patient_id", "gene_a", "gene_b", "gene_c")
]

# Prediction is outcome-blind: time and event are not supplied here.
survival_prediction <- predict_final(
  survival_fit,
  survival_new_data,
  id_col = "patient_id"
)

# Outcomes enter only after predictions have been saved.
survival_endpoint <- surv_endpoint(
  survival_external,
  id_col = "patient_id",
  time_col = "time",
  event_col = "event",
  time_unit = "month",
  prediction_horizon = 6
)
survival_evaluation <- eval_external(
  survival_prediction,
  survival_endpoint
)
print(survival_evaluation)
survival_evaluation$metrics$harrell_c
survival_evaluation$metrics$ipcw_brier

# ---- Binary workflow -------------------------------------------------------

binary_development <- bi_example(n = 100, seed = 20260819)
binary_task <- bi_task(
  data = binary_development,
  id_col = "patient_id",
  outcome_col = "outcome",
  predictors = c("gene_a", "gene_b", "gene_c")
)

binary_pipelines <- list(
  logistic = sb_pipeline(learner = "logistic"),
  gaussian_naive_bayes = sb_pipeline(learner = "gaussian_naive_bayes")
)
binary_cv <- compare_cv(
  binary_task,
  binary_pipelines,
  folds = 5L,
  seed = 5201314L
)
print(binary_cv)
binary_cv$comparison

# Again, the package does not select from the comparison automatically.
binary_final_pipeline <- binary_pipelines[["logistic"]]
binary_fit <- fit_final(binary_task, binary_final_pipeline)

binary_external <- bi_example(n = 40, seed = 20260821)
binary_external$patient_id <- paste0(
  "binary_external_",
  seq_len(nrow(binary_external))
)
binary_new_data <- binary_external[
  c("patient_id", "gene_a", "gene_b", "gene_c")
]
binary_prediction <- predict_final(
  binary_fit,
  binary_new_data,
  id_col = "patient_id"
)
binary_endpoint <- bi_endpoint(
  binary_external,
  id_col = "patient_id",
  outcome_col = "outcome"
)
binary_evaluation <- eval_external(
  binary_prediction,
  binary_endpoint
)
print(binary_evaluation)
binary_evaluation$metrics$roc_auc
binary_evaluation$metrics$brier_score
binary_evaluation$metrics$calibration

# Minimal integration assertions. Detailed statistical reconstruction tests
# belong to tests/testthat rather than this user-facing example.
stopifnot(
  identical(survival_cv$status, "success"),
  identical(survival_fit$status, "success"),
  identical(survival_prediction$status, "success"),
  identical(survival_evaluation$status, "success"),
  !any(c("time", "event") %in% names(survival_new_data)),
  identical(binary_cv$status, "success"),
  identical(binary_fit$status, "success"),
  identical(binary_prediction$status, "success"),
  identical(binary_evaluation$status, "success"),
  !"outcome" %in% names(binary_new_data)
)
    invisible(as.list(environment(), all.names = FALSE))
  },
  "multi-cohort" = function() {
make_binary_task <- function(prefix, n, seed) {
  d <- bi_example(n = n, seed = seed)
  d$patient_id <- paste0(prefix, d$patient_id)
  bi_task(d, "patient_id", "outcome",
    predictors = c("gene_a", "gene_b", "gene_c"))
}

multi_cohort_result <- sb_compare_cohorts(
  development = make_binary_task("dev_", 100, 42),
  cohorts = list(
    Validation1 = make_binary_task("val1_", 60, 43),
    Validation2 = make_binary_task("val2_", 60, 44),
    Test = make_binary_task("test_", 80, 45)
  ),
  roles = c(Validation1 = "validation", Validation2 = "validation", Test = "test"),
  models = c("logistic", "gaussian_naive_bayes"),
  development_name = "Training", folds = 5, seed = 42
)

print(multi_cohort_result)
multi_cohort_result$cohorts$Test$comparison

# For a persistent report, choose a new filename:
# sb_view_results(multi_cohort_result, path = "cohort-results.html", open = FALSE)
    invisible(as.list(environment(), all.names = FALSE))
  },
  "multiclass-classical" = function() {
# Iris is a public demonstration dataset, not clinical validation evidence.
d <- datasets::iris
d$id <- seq_len(nrow(d))
task <- multi_task(d, "id", "Species", names(d)[1:4],
  classes = c("virginica", "setosa", "versicolor"))
methods <- c("lda", "qda", "gaussian_naive_bayes", "decision_tree")
analysis <- analyze_methods(task, methods)
analysis$status
sb_method_features(analysis, "gaussian_naive_bayes", class = "setosa")
subset(sb_method_features(analysis, "decision_tree"), used_for_split)
analysis$results$gaussian_naive_bayes$plots$class_statistics
sb_predict_multiclass(analysis$results$qda$fit, d[1:3, ])
comparison <- sb_compare_multiclass(task, methods, folds = 3)
comparison$comparison
# No automatic files, reports or confusion plots.
    invisible(as.list(environment(), all.names = FALSE))
  },
  "multiclass-quickstart" = function() {
# Three-class example. Replace iris with a numeric gene/clinical feature table.
data <- iris
data$patient_id <- seq_len(nrow(data))
genes <- names(data)[1:4]
task <- multi_task(data, "patient_id", "Species", predictors = genes)

# Independent full-data feature exploration, not held-out performance.
analysis <- analyze_methods(task, c("Boruta", "varSelRF", "VSURF", "ABESS"))
analysis$status
selected_genes(analysis, "Boruta")
selected_genes(analysis, "ABESS")
analysis$results$abess$plots$support_criterion

# Each learner receives the original candidate features inside its training fold.
# Do not substitute genes selected on the full data above.
comparison <- sb_compare_multiclass(task,
  c("svm_radial", "ranger", "xgboost", "lightgbm", "catboost"), folds = 3)
comparison$comparison

# Fit one chosen model and predict from its saved native fit.
model <- analyze_methods(task, "ranger")
probability <- sb_predict_multiclass(model$results$ranger$fit, data[1:3, ])
probability
# This last prediction uses training rows only to demonstrate the API.
# For clinical use, pass genuinely new patients with the same predictor columns.
    invisible(as.list(environment(), all.names = FALSE))
  },
  "multiclass" = function() {
# Public iris data demonstrate the interface; these are not clinical results.
d <- datasets::iris
d$sample_id <- seq_len(nrow(d))
task <- multi_task(d, "sample_id", "Species", names(d)[1:4],
  classes = c("virginica", "setosa", "versicolor"))

analysis <- analyze_methods(task, c("lasso", "elastic_net", "random_forest"),
  inner_folds = 3)
analysis$status
selected_genes(analysis, "lasso", class = "setosa")
analysis$results$lasso$plots$coefficient_path
sb_predict_multiclass(analysis$results$lasso$fit, d[1:3, ])

comparison <- sb_compare_multiclass(task, c("lasso", "elastic_net", "random_forest"),
  folds = 3, inner_folds = 3)
comparison$comparison
comparison$results$lasso$metrics$per_class
# Optional only: plot_multi_cm(comparison$results$lasso$metrics)
    invisible(as.list(environment(), all.names = FALSE))
  },
  "paper-figures" = function() {
d <- bi_example(n = 160, seed = 42)
genes <- c("gene_a", "gene_b", "gene_c")
task <- bi_task(d, "patient_id", "outcome", genes,
  decision_thresholds = seq(.05, .8, by = .05))
pipeline <- sb_pipeline(learner = "logistic")
fit <- fit_final(task, pipeline)
pred <- predict_final(fit, d[c("patient_id", genes)], "patient_id")
training_roc <- plot_pred_roc(d$outcome, pred$predictions$probability,
  kind = "apparent", cohort = "Training", ci_reps = 1000)

external <- bi_example(n = 120, seed = 43)
external$patient_id <- paste0("external_", external$patient_id)
ep <- predict_final(fit, external[c("patient_id", genes)], "patient_id")
evaluation <- eval_external(ep, bi_endpoint(external, "patient_id", "outcome"))
external_roc <- plot_roc(evaluation, cohort_label = "Validation", ci_reps = 1000)

cv <- run_cv(task, pipeline, folds = 3, seed = 42)
oof_roc <- plot_roc(cv)
nomogram <- plot_nomogram(nomogram(fit, d), reference_data = d,
  patient = d[1, , drop = FALSE])
calibration <- sb_bootstrap_calibration(task, pipeline, reps = 1000, seed = 42)
calibration_plot <- plot_boot_calib(calibration)
dca <- plot_dca(cv, standardized = TRUE, cost_benefit_axis = TRUE)


performance <- data.frame(model = rep(c("Model A", "Model B", "Model C"), each = 3),
  cohort = rep(c("Development", "Validation1", "Validation2"), 3),
  auc = c(.85, .80, .82, .90, .78, .76, .82, .79, .80))
ranking <- plot_rank(performance, c("Validation1", "Validation2"))


training_roc
ranking$plot
nomogram
calibration_plot
dca
    invisible(as.list(environment(), all.names = FALSE))
  },
  "quick-start" = function() {
genes <- c("gene_a", "gene_b", "gene_c")

binary_data <- bi_example(n = 100, seed = 42)
binary_task <- bi_task(
  binary_data, "patient_id", "outcome", predictors = genes
)
binary_result <- sb_compare_models(
  binary_task, models = c("logistic", "gaussian_naive_bayes"),
  folds = 5, seed = 42
)
binary_result$comparison[c("model_id", "status", "roc_auc", "brier_score")]

survival_data <- surv_example(n = 100, seed = 42)
survival_task <- surv_task(
  survival_data, "patient_id", "time", "event",
  predictors = genes, time_unit = "month"
)
survival_result <- sb_compare_models(
  survival_task, models = c("coxph", "survreg_weibull"),
  folds = 5, seed = 42
)
survival_result$comparison[c("model_id", "status", "harrell_c")]
    invisible(as.list(environment(), all.names = FALSE))
  },
  "short-quickstart" = function() {
development <- bi_example(120, seed = 42)
validation <- bi_example(80, seed = 17)
genes <- c("gene_a", "gene_b")

run <- bi_run(development, outcome = "outcome", predictors = genes,
              id = "patient_id", models = "Logistic", plots = FALSE)
finished <- finish_run(run, "Logistic", external = list(validation = validation))

views <- clinical_plots(finished)
roc <- views$plots$ROC
metrics <- result_table(finished$evaluations$validation, metric = "roc_auc")

if (interactive()) print(roc)
    invisible(as.list(environment(), all.names = FALSE))
  },
  "specialty-short-names" = function() {
## ---- specialty-lasso ----
development <- read.csv(system.file("extdata", "demo-binary-development.csv",
  package = "stormbreaker"), check.names = FALSE)
genes <- c("gene_a", "gene_b", "gene_c")
task <- bi_task(development, "patient_id", "outcome", genes)
fit <- fit_final(task, bi_pipeline("Lasso"), capture_diagnostics = TRUE)
cv <- path_diag(fit, type = "cv_curve")
path <- path_diag(fit, type = "coefficient_path")
native <- native_diag(fit, complete = TRUE)
figures <- list(CV = plot_diag(cv), Path = plot_diag(path))
example_result <- list(input = development, task = task, fit = fit,
  cv = cv, path = path, native = native, figures = figures)

if (interactive()) print(figures$CV)
    invisible(as.list(environment(), all.names = FALSE))
  },
  "user-models" = function() {
# User-controlled models and genes. These functions do not infer the endpoint.
# Replace paths, endpoint columns and genes with YOUR predeclared study inputs.
# Lines below are illustrative, not executed automatically by the package.
# Read your own gene list: one exact predictor column name per line.
genes <- readLines("D:/my_study/my_genes.txt", encoding = "UTF-8")
development <- read.csv("D:/my_study/development.csv", check.names = FALSE,
  colClasses = c(patient_id = "character"))

# Choose ONE task definition yourself:
task <- bi_task(development, "patient_id", "outcome", predictors = genes)
# For survival use this INSTEAD:
# task <- surv_task(development, "patient_id", "time", "event",
#   predictors = genes, time_unit = "month")

sb_model_catalog("binary") # use "survival" for a survival task

# Specific requirements: run just these models, on the same user gene set.
compared <- sb_compare_models(task,
  models = c("logistic", "ridge_logistic", "ranger_binary"),
  folds = 5L, seed = 2026L)

# No model preference: use this INSTEAD to run the complete task-specific menu.
# compared <- sb_compare_models(task, models = NULL, folds = 5L, seed = 2026L)

compared$comparison
compared$elapsed_seconds
# Use new file paths for each run; do not replace past patient results.
# write.csv(compared$comparison, "D:/my_study/run_01/comparison.csv", row.names = FALSE)
# saveRDS(compared, "D:/my_study/run_01/comparison.rds")

# After development-only review, explicitly choose a compared specification.
# selected <- "logistic"  # example, NOT an automatically selected winner
# final_fit <- fit_final(task, compared$model_pipelines[[selected]])
# saveRDS(final_fit, "D:/my_study/run_01/model.rds")
# Predict external predictors first, then supply outcomes for external evaluation;
# see guides/user-guide.md for the outcome-blind sequence.
    invisible(as.list(environment(), all.names = FALSE))
  }
)

advanced_examples <- function() names(.stormbreaker_advanced_examples)

run_advanced_example <- function(name) {
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !name %in% advanced_examples()) {
    stop("Choose one name returned by advanced_examples().")
  }
  library(stormbreaker)
  .stormbreaker_advanced_examples[[name]]()
}
