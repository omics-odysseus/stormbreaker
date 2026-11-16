# ---- method-analysis ----
local({
method_analysis_task <- function() {
  withr::with_seed(881L, {
    n <- 90L
    d <- data.frame(id = seq_len(n), a = rnorm(n), b = rnorm(n), c = rnorm(n))
    d$y <- rbinom(n, 1, plogis(2 * d$a - d$b))
    sb_binary_task(d, "id", "y", c("a", "b", "c"))
  })
}

test_that("catalog distinguishes task support and selection from ranking", {
  binary <- sb_method_catalog("binary")
  survival <- sb_method_catalog("survival")
  expect_identical(nrow(binary), 9L)
  expect_identical(binary$supported, rep(TRUE, 9))
  expect_identical(survival$method[!survival$supported], c("varselrf", "vsurf"))
  expect_identical(binary$method[binary$feature_result == "ranking"], "random_forest")
  expect_identical(sb_method_names(c("LASSO", "Elastic Net", "random-forest")),
    c("lasso", "elastic_net", "random_forest"))
})

test_that("actual LASSO and forest runs retain native results and reusable plots", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("ggplot2")
  task <- method_analysis_task()
  x <- sb_analyze_methods(task, c("LASSO", "random_forest"), inner_folds = 3L)
  expect_identical(x$status$status, rep("success", 2))
  lasso <- x$results$lasso
  expect_identical(lasso$fit$pipeline$selector, "none")
  expect_identical(lasso$fit$development$candidate_features, task$predictors)
  beta <- as.matrix(stats::coef(lasso$fit$diagnostic_artifacts$learner,
    s = lasso$fit$tuning$lambda))
  expected <- rownames(beta)[beta[, 1] != 0 & rownames(beta) != "(Intercept)"]
  expect_setequal(sb_selected_features(x, "lasso"), expected)
  rf <- sb_method_features(x, "random_forest")
  native <- unserialize(x$results$random_forest$fit$learner_state$engine_raw)
  expect_equal(rf$importance, unname(native$importance[rf$feature, "MeanDecreaseGini"]))
  expect_identical(rf$selected, rep(NA, nrow(rf)))
  expect_equal(rf$rank, rank(-rf$importance, ties.method = "min"))
  for (result in x$results) {
    expect_identical(result$plot_status$status, rep("available", length(result$plots)))
    for (plot in result$plots) {
      expect_null(plot$labels$subtitle)
      expect_null(plot$labels$caption)
      expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
    }
  }
  local_mocked_bindings(sb_fit_final = function(...) stop("unexpected refit"),
    sb_analyze_one_method = function(...) stop("unexpected analysis"))
  restored <- unserialize(serialize(x, NULL))
  expect_identical(sb_selected_features(restored, "LASSO"), sb_selected_features(x, "lasso"))
  expect_identical(sb_method_features(restored, "random_forest"), rf)
  expect_s3_class(tryCatch(sb_selected_features(restored, "random_forest"), error = identity),
    "stormbreaker_selection_not_defined")
})

test_that("method failures are retained and later methods still run", {
  called <- character()
  local_mocked_bindings(sb_analyze_one_method = function(task, spec, ...) {
    called <<- c(called, spec$method)
    if (spec$method == "boruta") stop("native failure")
    list(status = "success", method = spec$method, feature_result = spec$feature_result)
  })
  x <- sb_analyze_methods(method_analysis_task(), c("boruta", "lasso"))
  expect_identical(called, c("boruta", "lasso"))
  expect_identical(x$status$status, c("failed", "success"))
  expect_identical(x$status$reason, c("native failure", ""))
  expect_identical(x$results$boruta$error_message, "native failure")
  expect_s3_class(tryCatch(sb_selected_features(x, "boruta"), error = identity), "error")
})

test_that("ambiguous requests fail before fitting and empty selections remain empty", {
  task <- method_analysis_task()
  for (methods in list(c("LASSO", "lasso"), "svm_rfe", character(), NA_character_)) {
    expect_s3_class(tryCatch(sb_analyze_methods(task, methods), error = identity),
      "stormbreaker_bad_method")
  }
  for (alpha in c(0, 1, NA_real_)) {
    expect_s3_class(tryCatch(sb_analyze_methods(task, "elastic_net", elastic_net_alpha = alpha),
      error = identity), "stormbreaker_bad_method")
  }
  x <- structure(list(results = list(lasso = list(status = "success",
    feature_result = "selection", features = data.frame(feature = "a", selected = FALSE)))),
    class = "stormbreaker_method_analysis")
  expect_identical(sb_selected_features(x, "lasso"), character())
})
})

# ---- selection-summary ----
local({
# ---- gene-overlap ----

test_that("gene intersections preserve exact exclusive members", {
  x <- sb_gene_overlap(list(A = c("a", "b", "b"), B = c("b", "c"), C = character()))
  expect_identical(x$common, character())
  expect_equal(sum(x$intersections$n), 3)
  expect_identical(x$intersections$genes[[match("A & B", x$intersections$intersection)]], "b")
  expect_identical(colnames(x$membership), c("A", "B", "C"))
  expect_s3_class(ggplot2::ggplotGrob(sb_plot_gene_overlap(x, "bar")), "gtable")
  empty <- sb_gene_overlap(list(A = character(), B = character()))
  expect_equal(nrow(empty$intersections), 0)
  expect_equal(dim(empty$membership), c(0L, 2L))
})

test_that("Venn diagrams render two to four gene sets", {
  skip_if_not_installed("ggvenn")
  sets <- list(A = c("a", "b"), B = c("b", "c"), C = c("b", "d"), D = "b")
  for (k in 2:4) {
    x <- sb_gene_overlap(sets[seq_len(k)])
    p <- sb_plot_gene_overlap(x)
    expect_identical(attr(p, "gene_overlap")$common, "b")
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  }
})

# ---- gene-selection-ranking ----

make_gene_rank_result <- function(id = 'p1', selector = 'screen', counts = c(a=2L,b=1L,c=0L), n = 2L, learner = 'logistic') {
  summary <- list(status = if(n == 2L) 'success' else 'incomplete', n_successful_folds = n,
    feature_frequency = data.frame(feature=names(counts), selected_count=unname(counts), selection_frequency=unname(counts)/n),
    per_fold=data.frame(fold=1:2,status=c('success',if(n==2L)'success' else 'failed'),n_selected=c(2L,1L)))
  structure(list(status='success', task=list(predictors=names(counts)),
    pipeline=list(id=id, selector=selector, selector_params=list(k=2L), learner=learner, learner_params=list(alpha=1)),
    selection_stability=summary, glmnet_active_set_stability=summary,
    provenance=list(fold_assignment=data.frame(patient_id=c('secret-patient-1','secret-patient-2'),fold=1:2))),
    class='stormbreaker_binary_cv_result')
}
make_gene_rank_comparison <- function(...) structure(list(results=list(...)),class='stormbreaker_cv_comparison')

test_that('shared selector copies count once and retain zero and unobserved genes', {
  a <- make_gene_rank_result()
  b <- a; b$pipeline$id <- 'p2'; b$pipeline$learner <- 'another_learner'
  c <- make_gene_rank_result('p3','other',c(a=1L,b=2L,c=0L))
  r <- sb_gene_selection_ranking(make_gene_rank_comparison(a,b,c), candidate_features=c('a','b','c','missing'))
  expect_identical(r$ranking$feature,c('a','b','c','missing'))
  expect_equal(r$ranking$selected_count,c(3,3,0,NA))
  expect_equal(r$ranking$eligible_count,c(4,4,4,0))
  expect_equal(r$ranking$rank,c(1,1,3,NA))
  expect_equal(r$ranking$selection_frequency,c(.75,.75,0,NA))
  expect_equal(r$groups$pipeline_copies,c(2L,1L))
  expect_equal(r$status,'success')
  expect_identical(any(grepl('secret-patient',unlist(r))),FALSE)
  expect_s3_class(sb_plot_gene_selection_ranking(r), 'ggplot')
})

test_that('failed folds do not enter the eligible denominator', {
  a <- make_gene_rank_result(counts=c(a=1L,b=0L,c=0L),n=1L)
  a$status <- 'failed'
  r <- sb_gene_selection_ranking(a)
  expect_equal(r$ranking$selected_count,c(1,0,0))
  expect_equal(r$ranking$eligible_count,c(1,1,1))
  expect_equal(r$ranking$selection_frequency,c(1,0,0))
  expect_identical(r$status,'incomplete')
})

test_that('none fixed panels and all-failed summaries do not invent zero evidence', {
  for(selector in c('none','fixed_features')) {
    r <- sb_gene_selection_ranking(make_gene_rank_result(selector=selector))
    expect_identical(r$status,'not_available')
    expect_identical(all(is.na(r$ranking$selected_count)),TRUE)
  }
  a <- make_gene_rank_result(counts=c(a=0L,b=0L,c=0L),n=0L)
  r <- sb_gene_selection_ranking(a)
  expect_identical(all(is.na(r$ranking$selected_count)),TRUE)
})

test_that('conflicting saved counts are visible and not arbitrarily selected', {
  a <- make_gene_rank_result()
  b <- make_gene_rank_result('p2',counts=c(a=1L,b=1L,c=0L))
  c <- make_gene_rank_result('p3','other')
  r <- sb_gene_selection_ranking(make_gene_rank_comparison(a,b,c))
  expect_identical(r$groups$status,c('conflicting_counts','success'))
  expect_equal(r$ranking$selected_count,c(2,1,0))
  expect_equal(r$ranking$eligible_count,c(2,2,2))
  expect_identical(r$status,'incomplete')
})

test_that('different split assignments and CV repetitions remain separate', {
  a <- make_gene_rank_result()
  b <- a; b$pipeline$id <- 'p2'; b$provenance$fold_assignment$fold <- 2:1
  r <- sb_gene_selection_ranking(make_gene_rank_comparison(a,b))
  expect_equal(r$ranking$selected_count,c(4,2,0))
  repeated <- structure(list(repetitions=list(make_gene_rank_comparison(a),make_gene_rank_comparison(a))),class='stormbreaker_repeated_cv_comparison')
  q <- sb_gene_selection_ranking(repeated)
  expect_equal(q$ranking$selected_count,c(4,2,0))
  expect_equal(nrow(q$groups),2L)
})

test_that('glmnet nonzero counts exclude ridge and stay separate from selectors', {
  a <- make_gene_rank_result(selector='none',learner='glmnet_logistic')
  b <- a; b$pipeline$id <- 'ridge'; b$pipeline$learner_params$alpha <- 0
  r <- sb_gene_selection_ranking(make_gene_rank_comparison(a,b),source='glmnet_active_set')
  expect_equal(r$ranking$selected_count,c(2,1,0))
  expect_equal(r$ranking$eligible_count,c(2,2,2))
  expect_identical(sb_gene_selection_ranking(a)$status,'not_available')
})

test_that('malformed counts and external evaluation inputs are rejected', {
  a <- make_gene_rank_result(counts=c(a=3L,b=1L,c=0L))
  expect_error(sb_gene_selection_ranking(a),class='stormbreaker_bad_gene_ranking')
  expect_error(sb_gene_selection_ranking(structure(list(),class='stormbreaker_external_evaluation')),class='stormbreaker_bad_gene_ranking')
})

# ---- selection-stability ----

make_selection_audit <- function(
    fold,
    selected_features,
    status = "success",
    error_class = NA_character_,
    error_message = NA_character_) {
  list(
    fold = as.integer(fold),
    selector = "test_selector",
    status = status,
    selected_features = selected_features,
    scores = NULL,
    error_class = error_class,
    error_message = error_message
  )
}

test_that("selection frequency and Jaccard match a hand calculation", {
  selections <- list(
    make_selection_audit(1L, c("a", "b")),
    make_selection_audit(2L, c("a", "c")),
    make_selection_audit(3L, c("a", "b"))
  )

  audit <- sb_selection_stability(
    "test_selector",
    selections,
    c("d", "b", "a", "c")
  )

  expect_named(
    audit,
    c(
      "selector", "status", "scope", "aggregation",
      "frequency_denominator", "n_candidate_features", "n_folds",
      "n_successful_folds", "n_pairs", "n_evaluable_pairs",
      "mean_pairwise_jaccard", "feature_frequency", "per_fold",
      "pairwise", "error_message"
    ),
    ignore.order = FALSE
  )
  expect_identical(audit$status, "success")
  expect_identical(audit$scope, "outer_analysis_fold_selected_feature_sets")
  expect_identical(
    audit$aggregation,
    "mean_of_unique_fold_pair_jaccard"
  )
  expect_identical(audit$frequency_denominator, "successful_selection_folds")
  expect_identical(audit$n_candidate_features, 4L)
  expect_identical(audit$n_folds, 3L)
  expect_identical(audit$n_successful_folds, 3L)
  expect_identical(audit$n_pairs, 3L)
  expect_identical(audit$n_evaluable_pairs, 3L)
  expect_identical(audit$feature_frequency$feature, c("d", "b", "a", "c"))
  expect_identical(audit$feature_frequency$selected_count, c(0L, 2L, 3L, 1L))
  expect_equal(
    audit$feature_frequency$selection_frequency,
    c(0, 2 / 3, 1, 1 / 3),
    tolerance = 0
  )
  expect_identical(audit$pairwise$fold_i, c(1L, 1L, 2L))
  expect_identical(audit$pairwise$fold_j, c(2L, 3L, 3L))
  expect_identical(audit$pairwise$intersection, c(1L, 2L, 1L))
  expect_identical(audit$pairwise$union, c(3L, 2L, 3L))
  expect_equal(audit$pairwise$jaccard, c(1 / 3, 1, 1 / 3), tolerance = 0)
  expect_equal(audit$mean_pairwise_jaccard, 5 / 9, tolerance = 0)
  expect_true(all(audit$pairwise$status == "success"))
  expect_true(is.na(audit$error_message))
})

test_that("fold and selected-feature order cannot alter stability", {
  original <- list(
    make_selection_audit(1L, c("a", "b")),
    make_selection_audit(2L, c("a", "c")),
    make_selection_audit(3L, c("a", "b"))
  )
  reordered <- list(
    make_selection_audit(3L, c("b", "a")),
    make_selection_audit(1L, c("b", "a")),
    make_selection_audit(2L, c("c", "a"))
  )

  expected <- sb_selection_stability(
    "test_selector",
    original,
    c("d", "b", "a", "c")
  )
  observed <- sb_selection_stability(
    "test_selector",
    reordered,
    c("d", "b", "a", "c")
  )

  expect_identical(observed, expected)
})

test_that("no learned selector is explicitly not applicable", {
  selections <- list(
    make_selection_audit(1L, c("a", "b")),
    make_selection_audit(2L, c("a", "b"))
  )

  audit <- sb_selection_stability("none", selections, c("a", "b"))

  expect_identical(audit$status, "not_applicable")
  expect_identical(audit$n_folds, 2L)
  expect_identical(audit$n_successful_folds, 0L)
  expect_identical(audit$n_pairs, 1L)
  expect_identical(audit$n_evaluable_pairs, 0L)
  expect_true(is.na(audit$mean_pairwise_jaccard))
  expect_identical(audit$feature_frequency$selected_count, c(0L, 0L))
  expect_true(all(is.na(audit$feature_frequency$selection_frequency)))
  expect_true(all(audit$per_fold$status == "not_applicable"))
  expect_true(all(audit$pairwise$status == "not_applicable"))
})

test_that("failed selector folds retain partial evidence without a headline", {
  selections <- list(
    make_selection_audit(3L, c("a", "c")),
    make_selection_audit(
      2L,
      character(),
      status = "failed",
      error_class = "stormbreaker_selection_failure",
      error_message = "selector failed"
    ),
    make_selection_audit(1L, c("a", "b"))
  )

  audit <- sb_selection_stability(
    "test_selector",
    selections,
    c("a", "b", "c", "d")
  )

  expect_identical(audit$status, "incomplete")
  expect_identical(audit$n_successful_folds, 2L)
  expect_identical(audit$n_pairs, 3L)
  expect_identical(audit$n_evaluable_pairs, 1L)
  expect_true(is.na(audit$mean_pairwise_jaccard))
  expect_identical(audit$feature_frequency$selected_count, c(2L, 1L, 1L, 0L))
  expect_equal(
    audit$feature_frequency$selection_frequency,
    c(1, 1 / 2, 1 / 2, 0),
    tolerance = 0
  )
  expect_identical(audit$per_fold$fold, c(1L, 2L, 3L))
  expect_identical(audit$per_fold$status, c("success", "failed", "success"))
  expect_identical(
    audit$per_fold$error_class,
    c(NA_character_, "stormbreaker_selection_failure", NA_character_)
  )
  expect_identical(audit$per_fold$error_message[[2L]], "selector failed")
  expect_true(is.na(audit$per_fold$n_selected[[2L]]))
  expect_identical(audit$pairwise$status, c("incomplete", "success", "incomplete"))
  expect_equal(audit$pairwise$jaccard, c(NA, 1 / 3, NA), tolerance = 0)
  expect_match(audit$error_message, "not every outer-fold selection succeeded")
})
})

# ---- selector-binary-models ----
local({
binary_model_selector_task <- function(grouped = FALSE) {
  withr::local_seed(5)
  d <- as.data.frame(matrix(rnorm(180 * 8), 180, 8))
  names(d) <- paste0("gene", 1:8)
  d$id <- paste0("p", 1:180)
  d$y <- rbinom(180, 1, plogis(1.5 * d$gene1 - d$gene2))
  if (grouped) {
    d$patient <- d$id
    d <- d[rep(seq_len(nrow(d)), each = 2), ]
    d$id <- paste0("row", seq_len(nrow(d)))
  }
  sb_binary_task(d, "id", "y", paste0("gene", 1:8), group_col = if (grouped) "patient" else NULL)
}

test_that("binary GamBoost selection matches native active spline terms", {
  skip_if_not_installed("mboost")
  task <- binary_model_selector_task()
  x <- task$data[task$predictors]
  names(x) <- paste0("gene ", seq_len(ncol(x)))
  p <- sb_pipeline("gamboost_binary_nonzero", "logistic")
  selected <- sb_fit_binary_model_selector(x, task$data$y, p$selector, p$selector_params)
  work <- x; names(work) <- paste0("feature_", seq_len(ncol(x)))
  native <- sb_fit_classical_binary(work, task$data$y, "gamboost_binary",
    sb_classical_binary_params("gamboost_binary"))$engine
  beta <- stats::coef(native)
  active <- names(beta)[vapply(beta, function(b) any(b != 0), logical(1))]
  index <- match(sub("^bbs\\(([^,)]+).*", "\\1", active), names(work))
  expect_setequal(selected$selected_features, names(x)[index])
  expect_identical(sb_doctor_selector_menu("binary")[["GamBoost"]], p$selector)
  expect_null(selected$engine)
  expect_identical(selected, unserialize(serialize(selected, NULL)))
})

test_that("binary model selectors round-trip and retain genuine selected genes", {
  task <- binary_model_selector_task()
  for (id in sb_binary_model_selector_ids()) {
    pipeline <- sb_pipeline(id, "gbm_binary")
    expect_identical(sb_lint_pipeline(pipeline)$status, "admitted")
    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    if (id == "svm_rfe_binary") {
      evidence <- fit$selection$state$evidence
      expect_identical(fit$selection$selected_features,
        evidence$subsets[[which.max(evidence$curve$accuracy)]])
    } else expect_identical(fit$selection$selected_features,
      fit$selection$scores$feature[fit$selection$scores$selected])
    expect_identical(all(fit$selection$selected_features %in% task$predictors), TRUE)
    copied <- unserialize(serialize(fit, NULL))
    pred <- sb_predict_final(copied, task$data[c("id", task$predictors)], "id")
    expect_identical(pred$status, "success")
  }
  a <- sb_pipeline("stepglm_binary_aic", "logistic", selector_options = list(direction = "forward"))
  b <- sb_pipeline("stepglm_binary_aic", "logistic", selector_options = list(direction = "backward"))
  expect_identical(a$id == b$id, FALSE)
  expect_identical(sb_lint_pipeline(a)$status, "admitted")
  bad <- tryCatch(sb_pipeline("rf_binary_importance", "logistic", selector_options = list(thresold = 1)), error = identity)
  expect_s3_class(bad, "stormbreaker_bad_selector_options")
  empty <- sb_fit_final(task, sb_pipeline("rf_binary_importance", "logistic", selector_options = list(threshold = 1e6)))
  expect_identical(empty$status, "failed")
  expect_identical(empty$error_stage, "selection")
})

test_that("selectors execute within outer and learner-inner training folds", {
  task <- binary_model_selector_task()
  original <- sb_fit_binary_model_selector
  seen <- list()
  local_mocked_bindings(sb_fit_binary_model_selector = function(x, outcome, id, params, group_ids = NULL) {
    ids <- rownames(x)
    if (all(ids %in% rownames(task$data))) ids <- task$data$id[match(ids, rownames(task$data))]
    seen[[length(seen) + 1L]] <<- ids
    original(x, outcome, id, params, group_ids)
  })
  for (id in sb_binary_model_selector_ids()) {
    seen <- list()
    cv <- sb_run_cv(task, sb_pipeline(id, "glmnet_logistic", inner_folds = 3, inner_seed = 7,
      selector_options = if (id %in% c("lasso_binary_nonzero", "enet_binary_nonzero")) list(lambda_rule = "lambda.min") else NULL), folds = 3)
    if (id == "svm_rfe_binary") {
      # This fixture selects one feature in a nested fold; native glmnet
      # cannot fit that design. Failure must remain visible, without padding.
      expect_identical(cv$status, "failed")
      failed <- cv$runs[cv$runs$status == "failed", ]
      expect_identical(failed$error_class, "stormbreaker_glmnet_logistic_tuning_failure")
      expect_match(failed$error_message, "cannot retain 1 or less variables", fixed = TRUE)
    } else expect_identical(cv$status, "success")
    expect_gt(length(seen), 3)
    expect_identical(all(lengths(seen) < task$n), TRUE)
    expect_identical(any(lengths(seen) < floor(task$n * 2/3)), TRUE)
    for (ids in seen) expect_identical(any(vapply(cv$folds,
      function(fold) all(ids %in% fold$analysis_ids), logical(1))), TRUE)
  }
})

test_that("LASSO inner folds keep repeated patients together", {
  task <- binary_model_selector_task(TRUE)
  original <- sb_fit_binary_model_selector
  grouped_calls <- 0L
  local_mocked_bindings(sb_fit_binary_model_selector = function(x, outcome, id, params, group_ids = NULL) {
    result <- original(x, outcome, id, params, group_ids)
    if (id == "lasso_binary_nonzero") {
      expect_length(group_ids, nrow(x))
      expect_identical(all(vapply(split(result$evidence$fold_id, group_ids),
        function(z) length(unique(z)) == 1L, logical(1))), TRUE)
      grouped_calls <<- grouped_calls + 1L
    }
    result
  })
  p <- sb_pipeline("lasso_binary_nonzero", "glmnet_logistic", inner_folds = 3, inner_seed = 7,
    selector_options = list(lambda_rule = "lambda.min"))
  fit <- sb_fit_final(task, p)
  expect_identical(fit$status, "success")
  cv <- sb_run_cv(task, p, folds = 3)
  expect_identical(cv$status, "success")
  expect_gt(grouped_calls, 5)
})

test_that("new selectors refit in SVM calibration and preserve patient groups", {
  task <- binary_model_selector_task(TRUE)
  original <- sb_fit_binary_model_selector
  calls <- 0L
  local_mocked_bindings(sb_fit_binary_model_selector = function(x, outcome, id, params, group_ids = NULL) {
    result <- original(x, outcome, id, params, group_ids)
    calls <<- calls + 1L
    if (id == "lasso_binary_nonzero") {
      expect_length(group_ids, nrow(x))
      expect_identical(all(vapply(split(result$evidence$fold_id, group_ids),
        function(z) length(unique(z)) == 1L, logical(1))), TRUE)
    }
    result
  })
  for (id in sb_binary_model_selector_ids()) {
    calls <- 0L
    fit <- sb_fit_final(task, sb_pipeline(id, "svm_linear_binary",
      selector_options = if (id == "lasso_binary_nonzero") list(lambda_rule = "lambda.min") else NULL))
    expect_identical(fit$status, "success")
    expect_gt(calls, 1L)
    expect_length(fit$tuning$inner_selections, 5L)
  }
})

test_that("RF scores match native importance with original gene names", {
  task <- binary_model_selector_task()
  x <- task$data[task$predictors]
  names(x)[1:2] <- c("HLA-DRA", "gene with space")
  params <- sb_binary_model_selector_options("rf_binary_importance")
  selected <- sb_fit_binary_model_selector(x, task$data$y, "rf_binary_importance", params)
  work <- x
  names(work) <- paste0("feature_", seq_len(ncol(x)))
  withr::local_seed(params$seed)
  native <- randomForest::randomForest(x = work, y = factor(task$data$y, levels = c(0, 1)),
    ntree = params$trees, importance = TRUE, keep.forest = FALSE)
  score <- as.numeric(randomForest::importance(native, type = 1, scale = FALSE)[, 1])
  expect_equal(selected$scores$score, score)
  expect_identical(selected$selected_features, names(x)[score > params$threshold])
  for (id in c("stepglm_binary_aic", "glmboost_binary_nonzero")) {
    result <- sb_fit_binary_model_selector(x, task$data$y, id, sb_binary_model_selector_options(id))
    expect_identical(all(result$selected_features %in% names(x)), TRUE)
    expect_identical(all(result$evidence$selected_features %in% names(x)), TRUE)
  }
})
})

# ---- selector-forest ----
local({
weight_screen_data <- function() {
  withr::local_seed(92)
  d <- as.data.frame(matrix(rnorm(160 * 8), 160, 8))
  names(d) <- paste0("gene", 1:8)
  d$id <- paste0("p", 1:160)
  d$time <- rexp(160, exp(d$gene1 - d$gene2) / 10) + .01
  d$event <- rbinom(160, 1, .8)
  d
}

test_that("default survival weight screens select original genes and render", {
  d <- weight_screen_data()
  t <- sb_surv_task(d, "id", "time", "event", paste0("gene", 1:8), "month")
  for (method in c("Ridge", "PLSR", "SVM")) {
    p <- surv_pipeline("Cox", surv_selector(method))
    f <- sb_fit_final(t, p, capture_diagnostics = TRUE)
    expect_identical(f$status, "success")
    s <- f$selection$state
    expect_length(s$selected_features, 4)
    expect_equal(s$scores$score, abs(s$scores$coefficient * s$scores$training_sd))
    expect_setequal(s$selected_features, s$scores$feature[order(-s$scores$score)[1:4]])
    expect_identical(sb_predict_final(f, d[1:5, ], "id")$status, "success")
    z <- sb_plots(f)
    expect_length(z$unavailable, 0)
    expect_s3_class(z$plots$selector_1, "ggplot")
    for (plot in z$plots) expect_s3_class(ggplot2::ggplotGrob(plot +
      ggplot2::theme(text = ggplot2::element_text(family = "sans"))), "gtable")
    cv <- sb_run_cv(t, p, folds = 3, seed = 12)
    expect_identical(cv$status, "success")
  }
  f <- sb_fit_final(t, surv_pipeline("Cox", surv_selector(c("Ridge", "PLSR"))))
  expect_identical(f$status, "success")
  expect_length(f$selection$selected_features, 2)
})

test_that("Ridge screening keeps repeated patients in the same tuning fold", {
  d <- weight_screen_data()
  groups <- rep(seq_len(80), each = 2)
  p <- sb_survival_weight_selector_params("ridge_cox_weights")
  s <- sb_survival_weight_selector_adapter("ridge_cox_weights", p)$fit(
    d[paste0("gene", 1:8)], d$time, d$event, group_ids = groups)
  expect_equal(vapply(split(s$evidence$fold_id, groups), function(x) length(unique(x)), integer(1)),
    stats::setNames(rep(1L, 80), as.character(1:80)))
})
})

# ---- selector-native ----
local({
# ---- boruta-worker-reuse ----

test_that("Boruta worker handoff preserves matching and analysis boundaries", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    list(selected_features = names(x), value = calls)
  })
  selector <- sb_resolve_selector_adapter("boruta_binary_roughfix", list())
  x <- data.frame(a = 1:8, b = 8:1)
  endpoint <- list(task_type = "binary", outcome = rep(0:1, 4))
  saved <- sb_with_selector_reuse({
    fit <- sb_fit_runtime_selector(selector, x, endpoint)
    list(fit = fit, checkpoint = sb_boruta_worker_checkpoint(selector$id))
  })
  saved <- unserialize(serialize(saved, NULL))
  result <- sb_with_selector_reuse({
    sb_restore_boruta_worker_checkpoint(saved$checkpoint)
    expect_identical(sb_fit_runtime_selector(selector, x, endpoint), saved$fit)
    sb_fit_runtime_selector(selector, x[-1, ],
      list(task_type = "binary", outcome = endpoint$outcome[-1]))
  }, diagnostics = TRUE)
  expect_identical(calls, 2L)
  expect_identical(result$reuse$hits, 1L)
  expect_identical(result$reuse$misses, 1L)
  expect_null(.sb_selector_reuse_context$active)
})

test_that("failed consumers do not alter saved Boruta entries", {
  selector <- sb_resolve_selector_adapter("boruta_binary_roughfix", list())
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(...) list(selected_features = "a"))
  x <- data.frame(a = 1:6)
  endpoint <- list(task_type = "binary", outcome = rep(0:1, 3))
  saved <- sb_with_selector_reuse({
    sb_fit_runtime_selector(selector, x, endpoint)
    sb_boruta_worker_checkpoint(selector$id)
  })
  expect_error(sb_with_selector_reuse({
    sb_restore_boruta_worker_checkpoint(saved)
    stop("downstream failure")
  }), "downstream failure")
  result <- sb_with_selector_reuse({
    sb_restore_boruta_worker_checkpoint(saved)
    sb_fit_runtime_selector(selector, x, endpoint)
  }, diagnostics = TRUE)
  expect_identical(result$reuse$hits, 1L)
  expect_null(.sb_selector_reuse_context$active)
})

# ---- enet-gbm-selectors ----

enet_gbm_task <- function() {
  withr::local_seed(42)
  d <- as.data.frame(matrix(rnorm(240 * 8), 240, 8))
  names(d) <- c("HLA-DRA", "gene with space", paste0("g", 3:8))
  d$id <- paste0("p", seq_len(nrow(d)))
  d$y <- rbinom(nrow(d), 1, plogis(2 * d[[1]] - 1.5 * d[[2]] + d[[3]]))
  sb_binary_task(d, "id", "y", names(d)[1:8])
}

test_that("Enet selection uses its alpha and preserves patient-level folds", {
  task <- enet_gbm_task()
  x <- task$data[task$predictors]
  groups <- rep(seq_len(120), each = 2)
  params <- sb_binary_model_selector_options("enet_binary_nonzero")
  state <- sb_fit_binary_model_selector(x, task$data$y, "enet_binary_nonzero", params, groups)
  work <- x
  names(work) <- paste0("feature_", seq_len(ncol(x)))
  native <- glmnet::glmnet(as.matrix(work), task$data$y, family = "binomial",
    alpha = .5, lambda = state$evidence$lambda, standardize = TRUE)
  expect_equal(state$scores$score, as.numeric(stats::coef(native))[-1])
  expect_identical(state$selected_features, names(x)[state$scores$score != 0])
  expect_identical(all(vapply(split(state$evidence$fold_id, groups),
    function(g) length(unique(g)) == 1L, logical(1))), TRUE)
  expect_identical(bi_selector("ElasticNet"), bi_selector("Enet"))
  expect_error(sb_binary_model_selector_options("enet_binary_nonzero", list(alpha = 0)),
    class = "stormbreaker_bad_selector_options")
})

test_that("GBM selection matches native relative influence without invented cutoffs", {
  task <- enet_gbm_task()
  x <- task$data[task$predictors]
  state <- sb_fit_binary_model_selector(x, task$data$y, "gbm_binary_importance",
    sb_binary_model_selector_options("gbm_binary_importance"))
  work <- x
  names(work) <- paste0("feature_", seq_len(ncol(x)))
  withr::local_seed(1)
  native <- gbm::gbm.fit(as.matrix(work), task$data$y, distribution = "bernoulli",
    n.trees = 100, interaction.depth = 1, n.minobsinnode = 10, shrinkage = .1,
    bag.fraction = 1, nTrain = nrow(x), keep.data = FALSE, verbose = FALSE,
    var.names = names(work))
  imp <- summary(native, n.trees = 100, plotit = FALSE)
  expect_equal(state$scores$score, imp$rel.inf[match(names(work), imp$var)])
  expect_identical(state$selected_features, names(x)[state$scores$score > 0])
  expect_equal(state$evidence$training_loss$loss, native$train.error)
})

test_that("new selectors run in ordered chains with saved genes and native plots", {
  task <- enet_gbm_task()
  for (selectors in list("Enet", "GBM", c("Enet", "GBM"), c("GBM", "Lasso"))) {
    p <- bi_pipeline("RF", bi_selector(selectors))
    fit <- sb_fit_final(task, p, capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    expect_identical(sb_lint_pipeline(p)$status, "admitted")
    z <- sb_plots(fit)
    expect_contains(names(z$plots), "selector_1")
    if (selectors[1] == "GBM") {
      expect_contains(names(z$plots), "selector_1_training_loss")
      expect_identical(z$plots$selector_1$labels$x, "Importance")
    }
    if (length(selectors) == 2) {
      states <- fit$selection$state$stages
      expect_identical(states[[2]]$source_features, states[[1]]$selected_features)
    }
    cv <- sb_run_cv(task, p, folds = 3)
    expect_identical(cv$status, "success")
    expect_equal(sb_predict_final(fit, task$data, "id")$predictions,
      sb_predict_final(unserialize(serialize(fit, NULL)), task$data, "id")$predictions)
  }
})

# ---- forest-selectors ----

make_forest_selector_binary_fixture <- function(n = 100L, seed = 8711L) {
  withr::with_seed(seed, {
    outcome <- sample(rep(c(0L, 1L), length.out = n))
    data.frame(
      patient_id = sprintf("FSB%03d", seq_len(n)),
      outcome = outcome,
      gene_a = stats::rnorm(n, 1.2 * outcome),
      gene_b = stats::rnorm(n, -0.9 * outcome),
      gene_c = stats::rnorm(n, 0.7 * outcome),
      gene_d = stats::rnorm(n),
      gene_e = stats::rnorm(n),
      gene_f = stats::rnorm(n),
      check.names = FALSE
    )
  })
}

make_forest_selector_survival_fixture <- function(n = 100L, seed = 8712L) {
  withr::with_seed(seed, {
    gene_a <- stats::rnorm(n)
    gene_b <- stats::rnorm(n)
    gene_c <- stats::rnorm(n)
    gene_d <- stats::rnorm(n)
    gene_e <- stats::rnorm(n)
    gene_f <- stats::rnorm(n)
    event_time <- stats::rexp(
      n,
      rate = 0.06 * exp(1.1 * gene_a - 0.9 * gene_b + 0.5 * gene_c)
    ) + 0.1
    cutoff <- stats::median(event_time)
    data.frame(
      patient_id = sprintf("FSS%03d", seq_len(n)),
      time = round(pmin(event_time, cutoff), 4L),
      event = as.integer(event_time <= cutoff),
      gene_a = gene_a,
      gene_b = gene_b,
      gene_c = gene_c,
      gene_d = gene_d,
      gene_e = gene_e,
      gene_f = gene_f,
      check.names = FALSE
    )
  })
}

test_that("forest selectors expose fixed native specifications", {
  expect_identical(
    sb_forest_selector_ids(),
    c(
      "boruta_binary_roughfix",
      "boruta_survival_roughfix",
      "varselrf_binary_oob_1se"
    )
  )

  boruta <- sb_forest_selector_params("boruta_binary_roughfix")
  expect_identical(boruta$method, "Boruta")
  expect_identical(boruta$p_value, 0.01)
  expect_true(boruta$multiple_comparison_adjustment)
  expect_identical(boruta$maximum_runs, 100L)
  expect_identical(boruta$num_trees, 500L)
  expect_identical(boruta$num_threads, 1L)
  expect_identical(boruta$tentative_rule, "TentativeRoughFix")
  expect_identical(boruta$seed, 1L)

  varselrf <- sb_forest_selector_params("varselrf_binary_oob_1se")
  expect_identical(varselrf$method, "varSelRF")
  expect_identical(varselrf$c_sd, 1)
  expect_identical(varselrf$initial_trees, 5000L)
  expect_identical(varselrf$iteration_trees, 2000L)
  expect_identical(varselrf$drop_fraction, 0.2)
  expect_true(varselrf$whole_range)
  expect_false(varselrf$recompute_importance)
  expect_identical(varselrf$seed, 1L)
  expect_identical(
    names(sb_forest_selector_provenance("boruta_binary_roughfix")),
    c("boruta_version", "ranger_version")
  )
  expect_identical(
    names(sb_forest_selector_provenance("varselrf_binary_oob_1se")),
    c("varselrf_version", "randomforest_version")
  )
})

test_that("Boruta calls pin the ranger normalized-permutation backend", {
  skip_if_not_installed("Boruta", minimum_version = "9.0.0")

  captured_importance <- NULL
  observed <- testthat::with_mocked_bindings(
    sb_call_boruta(
      data.frame(gene_a = 1:4, gene_b = 4:1),
      factor(c(0L, 0L, 1L, 1L)),
      sb_forest_selector_params("boruta_binary_roughfix")
    ),
    Boruta = function(..., getImp) {
      captured_importance <<- getImp
      "captured"
    },
    .package = "Boruta"
  )

  expect_identical(observed, "captured")
  expect_identical(captured_importance, Boruta::getImpRfZ)
})

test_that("Boruta selectors reproduce native rough-fixed decisions", {
  skip_if_not_installed("Boruta", minimum_version = "9.0.0")
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  binary <- make_forest_selector_binary_fixture()
  x <- binary[paste0("gene_", letters[1:6])]
  adapter <- sb_forest_selector_adapter("boruta_binary_roughfix")
  params <- sb_forest_selector_params("boruta_binary_roughfix")

  set.seed(8713L)
  caller_rng <- .Random.seed
  fitted <- adapter$fit(x, binary$outcome)
  expect_identical(.Random.seed, caller_rng)

  native <- withr::with_seed(params$seed, Boruta::Boruta(
    x,
    factor(binary$outcome, levels = c(0L, 1L)),
    getImp = Boruta::getImpRfZ,
    pValue = params$p_value,
    mcAdj = params$multiple_comparison_adjustment,
    maxRuns = params$maximum_runs,
    doTrace = 0L,
    holdHistory = TRUE,
    num.trees = params$num_trees,
    num.threads = params$num_threads
  ))
  if (any(native$finalDecision == "Tentative")) {
    native <- Boruta::TentativeRoughFix(native)
  }

  expect_identical(
    fitted$selected_features,
    Boruta::getSelectedAttributes(native, withTentative = FALSE)
  )
  expect_identical(
    fitted$selected_features,
    fitted$scores$feature[fitted$scores$selected]
  )
  expect_null(fitted$engine)
  expect_null(fitted$call)

  survival <- make_forest_selector_survival_fixture()
  survival_x <- survival[paste0("gene_", letters[1:6])]
  survival_fit <- sb_forest_selector_adapter(
    "boruta_survival_roughfix"
  )$fit(survival_x, survival$time, survival$event)
  expect_gte(length(survival_fit$selected_features), 1L)
  expect_identical(
    names(sb_apply_forest_selector(survival_fit, survival_x)),
    survival_fit$selected_features
  )
  expect_null(survival_fit$engine)
  expect_null(survival_fit$call)
})

test_that("varSelRF selector reproduces its native one-SE OOB subset", {
  skip_if_not_installed("varSelRF", minimum_version = "0.7-9")
  skip_if_not_installed("randomForest", minimum_version = "4.7")

  binary <- make_forest_selector_binary_fixture(seed = 8721L)
  x <- binary[paste0("gene_", letters[1:6])]
  params <- sb_forest_selector_params("varselrf_binary_oob_1se")
  adapter <- sb_forest_selector_adapter("varselrf_binary_oob_1se")

  set.seed(8722L)
  caller_rng <- .Random.seed
  fitted <- adapter$fit(x, binary$outcome)
  expect_identical(.Random.seed, caller_rng)

  native <- withr::with_seed(params$seed, varSelRF::varSelRF(
    xdata = x,
    Class = factor(binary$outcome, levels = c(0L, 1L)),
    c.sd = params$c_sd,
    mtryFactor = params$mtry_factor,
    ntree = params$initial_trees,
    ntreeIterat = params$iteration_trees,
    vars.drop.num = NULL,
    vars.drop.frac = params$drop_fraction,
    whole.range = params$whole_range,
    recompute.var.imp = params$recompute_importance,
    verbose = FALSE,
    returnFirstForest = FALSE,
    keep.forest = FALSE
  ))

  expect_identical(fitted$selected_features, native$selected.vars)
  expect_gte(length(fitted$selected_features), 2L)
  expect_identical(nrow(fitted$selection_history), nrow(native$selec.history))
  expect_false("Vars.in.Forest" %in% names(fitted$selection_history))
  expect_null(fitted$engine)
  expect_null(fitted$call)
})

test_that("forest selectors add only legal generated routes", {
  registry <- sb_pipeline_registry()
  restricted_internal_cv <- c(
    setdiff(sb_rule_aft_learner_ids(), "bart_aft_survival"),
    sb_slope_binary_learner_ids(),
    "cancerclass_binary"
  )
  binary_learners <- registry$learner[
    registry$task_type == "binary" & registry$selector == "none" &
      !registry$learner %in% restricted_internal_cv
  ]
  survival_learners <- registry$learner[
    registry$task_type == "survival" & registry$selector == "none" &
      !registry$learner %in% c(
        sb_analysis_tuned_survival_learner_ids(),
        restricted_internal_cv
      )
  ]

  expect_setequal(
    registry$learner[registry$selector == "boruta_binary_roughfix"],
    binary_learners
  )
  expect_setequal(
    registry$learner[registry$selector == "varselrf_binary_oob_1se"],
    binary_learners
  )
  expect_setequal(
    registry$learner[registry$selector == "boruta_survival_roughfix"],
    survival_learners
  )
  expect_identical(
    nrow(registry[
      registry$selector == "boruta_survival_roughfix" &
        registry$learner %in% sb_analysis_tuned_survival_learner_ids(),
      ,
      drop = FALSE
    ]),
    0L
  )
})

test_that("forest selectors are repeated inside glmnet inner folds", {
  skip_if_not_installed("Boruta", minimum_version = "9.0.0")
  skip_if_not_installed("ranger", minimum_version = "0.18.0")
  skip_if_not_installed("varSelRF", minimum_version = "0.7-9")
  skip_if_not_installed("randomForest", minimum_version = "4.7")
  skip_if_not_installed("glmnet", minimum_version = "4.1.2")

  binary <- make_forest_selector_binary_fixture(n = 110L, seed = 8731L)
  binary_task <- sb_binary_task(
    binary,
    "patient_id",
    "outcome",
    paste0("gene_", letters[1:6])
  )
  survival <- make_forest_selector_survival_fixture(n = 110L, seed = 8732L)
  survival_task <- sb_surv_task(
    survival,
    "patient_id",
    "time",
    "event",
    paste0("gene_", letters[1:6]),
    "month"
  )

  fits <- list(
    sb_fit_final(binary_task, sb_pipeline(
      selector = "boruta_binary_roughfix",
      learner = "glmnet_logistic",
      inner_folds = 3L,
      inner_seed = 8733L
    )),
    sb_fit_final(binary_task, sb_pipeline(
      selector = "varselrf_binary_oob_1se",
      learner = "glmnet_logistic",
      inner_folds = 3L,
      inner_seed = 8734L
    )),
    sb_fit_final(survival_task, sb_pipeline(
      selector = "boruta_survival_roughfix",
      learner = "glmnet_cox",
      inner_folds = 3L,
      inner_seed = 8735L
    ))
  )

  for (fit in fits) {
    expect_identical(fit$status, "success")
    expect_identical(
      fit$tuning$inner_selection_method,
      "glmnet_exclude_function"
    )
    expect_length(fit$tuning$inner_selections, 3L)
    expect_true(all(vapply(fit$tuning$inner_selections, function(selection) {
      length(selection$analysis_patient_ids) > 0L &&
        length(selection$selected_features) > 0L
    }, logical(1))))
  }
  expect_identical(
    fits[[1L]]$provenance$boruta_version,
    as.character(utils::packageVersion("Boruta"))
  )
  expect_identical(
    fits[[1L]]$provenance$ranger_version,
    as.character(utils::packageVersion("ranger"))
  )
  expect_identical(
    fits[[2L]]$provenance$varselrf_version,
    as.character(utils::packageVersion("varSelRF"))
  )
  expect_identical(
    fits[[2L]]$provenance$randomforest_version,
    as.character(utils::packageVersion("randomForest"))
  )
  expect_identical(
    fits[[3L]]$provenance$boruta_version,
    as.character(utils::packageVersion("Boruta"))
  )
})

# ---- selector-pamr ----

test_that("PAMR selectors use native nonzero centroids for binary and survival tasks", {
  skip_if_not_installed("pamr")
  d <- sb_example_survival(180, 82)
  d$outcome <- d$event
  x <- d[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene A", "gene-B", "gene C")
  for (surv in c(FALSE, TRUE)) {
    id <- if (surv) "pamr_survival_nonzero" else "pamr_binary_nonzero"
    params <- sb_pamr_selector_params(id)
    a <- sb_pamr_selector_adapter(id, params)
    set.seed(19); before <- .Random.seed
    selected <- if (surv) a$fit(x, d$time, d$event) else a$fit(x, d$outcome)
    expect_identical(.Random.seed, before)
    data <- list(x = t(as.matrix(x)), geneid = names(x))
    if (surv) {
      data$survival.time <- d$time; data$censoring.status <- d$event
    } else data$y <- factor(d$outcome, levels = c(0, 1))
    invisible(capture.output(native <- pamr::pamr.train(data, threshold = 1,
      scale.sd = TRUE, offset.percent = 50, remove.zeros = FALSE,
      prior = if (surv) NULL else as.numeric(table(data$y) / nrow(x)),
      ngroup.survival = 2, sign.contrast = "both")))
    expected <- pamr::pamr.predict(native, data$x, threshold = 1, type = "nonzero")
    expect_identical(selected$selected_features, names(x)[expected])
    expect_identical(a$apply(unserialize(serialize(selected, NULL)), x), x[expected])
    expect_null(selected$engine)
    expect_null(selected$time)
    expect_null(selected$outcome)
    p <- sb_pipeline(id, if (surv) "coxph" else "logistic")
    expect_identical(sb_validate_pipeline_contract(p, p$task_type), p)
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (surv) sb_surv_task(d, "patient_id", "time", "event", features, "month") else
      sb_binary_task(d, "patient_id", "outcome", features)
    cv <- sb_run_cv(task, p, folds = 3, seed = 82)
    expect_identical(cv$status, "success")
    expect_identical(cv$provenance$pamr_version, as.character(packageVersion("pamr")))
    fit <- sb_fit_final(task, p)
    expect_identical(fit$status, "success")
    expect_identical(fit$provenance$pamr_version, as.character(packageVersion("pamr")))
    pred <- sb_predict_final(unserialize(serialize(fit, NULL)), d[c("patient_id", features)], "patient_id")
    expect_identical(pred$status, "success")
  }
})

test_that("PAMR selection is refitted within glmnet inner analysis folds", {
  skip_if_not_installed("pamr")
  skip_if_not_installed("glmnet")
  d <- withr::with_seed(82, {
    x <- matrix(rnorm(540), 180, 3)
    eta <- rowSums(x)
    t <- rexp(180, exp(eta) / 10)
    censor <- rexp(180, 0.02)
    data.frame(patient_id = paste0("p", seq_len(180)),
      time = pmin(t, censor), event = as.integer(t <= censor),
      outcome = rbinom(180, 1, plogis(eta)),
      gene_a = x[, 1], gene_b = x[, 2], gene_c = x[, 3])
  })
  features <- c("gene_a", "gene_b", "gene_c")
  original <- sb_pamr_selector_state
  seen <- list()
  local_mocked_bindings(sb_pamr_selector_state = function(engine, x, id, params) {
    seen[[length(seen) + 1L]] <<- d$patient_id[match(x[[1]], d$gene_a)]
    original(engine, x, id, params)
  })
  for (surv in c(FALSE, TRUE)) {
    task <- if (surv) sb_surv_task(d, "patient_id", "time", "event", features, "month") else
      sb_binary_task(d, "patient_id", "outcome", features)
    p <- sb_pipeline(if (surv) "pamr_survival_nonzero" else "pamr_binary_nonzero",
      if (surv) "glmnet_cox" else "glmnet_logistic", inner_folds = 3, inner_seed = 7)
    seen <- list()
    cv <- sb_run_cv(task, p, folds = 3, seed = 82)
    expect_gt(length(seen), 3L)
    expect_identical(any(lengths(seen) < 120L), TRUE)
    for (ids in seen) expect_identical(any(vapply(cv$folds,
      function(fold) all(ids %in% fold$analysis_ids), logical(1))), TRUE)
    expect_identical(cv$status, "success")
  }
})

test_that("PAMR does not pad a one-feature inner selection for glmnet", {
  skip_if_not_installed("pamr")
  skip_if_not_installed("glmnet")
  d <- sb_example_survival(180, 82); d$outcome <- d$event
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  cv <- sb_run_cv(task, sb_pipeline("pamr_binary_nonzero", "glmnet_logistic",
    inner_folds = 3, inner_seed = 7), folds = 3, seed = 82)
  expect_identical(cv$status, "failed")
  expect_identical(unique(cv$runs$error_class), "stormbreaker_glmnet_logistic_tuning_failure")
  expect_match(unique(cv$runs$error_message), "cannot retain 1 or less variables", fixed = TRUE)
})

# ---- selector-plsrglm ----

test_that("PLS logistic selection uses native coefficients in the analysis scale", {
  skip_if_not_installed("plsRglm")
  d <- sb_example_survival(180, 82); d$outcome <- d$event
  x <- d[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene A", "gene-B", "gene C")
  p <- sb_pipeline("plsrglm_binary_weights", "logistic")
  a <- sb_resolve_selector_adapter(p$selector, p$selector_params)
  set.seed(19); before <- .Random.seed
  selected <- a$fit(x, d$outcome)
  expect_identical(.Random.seed, before)
  native <- plsRglm::PLS_glm(dataY = d$outcome, dataX = x, nt = 2,
    modele = "pls-glm-logistic", typeVC = "none", scaleX = TRUE,
    pvals.expli = FALSE, MClassed = FALSE, sparse = FALSE, naive = FALSE, verbose = FALSE)
  score <- abs(drop(native$Coeffs)[-1] * vapply(x, sd, numeric(1)))
  expect_equal(selected$scores$score, unname(score), tolerance = 1e-12)
  expect_setequal(selected$selected_features, names(x)[order(-score)[1:2]])
  expect_null(selected$engine)
  expect_identical(a$apply(unserialize(serialize(selected, NULL)), x), x[selected$selected_features])
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  cv <- sb_run_cv(task, p, folds = 3, seed = 82)
  expect_identical(cv$status, "success")
  expect_identical(cv$provenance$plsrglm_version, as.character(packageVersion("plsRglm")))
  fit <- sb_fit_final(task, p)
  expect_identical(fit$status, "success")
  expect_identical(fit$provenance$plsrglm_version, as.character(packageVersion("plsRglm")))
  pred <- sb_predict_final(unserialize(serialize(fit, NULL)), d[c("patient_id", features)], "patient_id")
  expect_identical(pred$status, "success")
})

# ---- selector-svm-rfe ----

test_that("SVM-RFE short entry fits, predicts and draws its native curves", {
  skip_if_not_installed("e1071")
  set.seed(20)
  d <- data.frame(id = paste0("p", 1:120), y = rep(0:1, 60),
    a = rnorm(120), b = rnorm(120), c = rnorm(120))
  d$a <- d$a + 2 * d$y
  task <- sb_binary_task(d, "id", "y", c("a", "b", "c"))
  pipeline <- bi_pipeline("Logistic", bi_selector("SVM_RFE"))
  fit <- sb_fit_final(task, pipeline, capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  plots <- sb_plots(fit)
  expect_equal(length(plots$unavailable), 0)
  expect_equal(all(c("selector_1_cv_accuracy", "selector_1_cv_error") %in% names(plots$plots)), TRUE)
  expect_s3_class(ggplot2::ggplot_build(plots$plots$selector_1_cv_accuracy), "ggplot_built")
  expect_identical(sb_run_cv(task, pipeline, folds = 3, seed = 1)$status, "success")
  chain <- bi_pipeline("Logistic", bi_selector(c("GBM", "SVM_RFE")))
  chained <- sb_fit_final(task, chain, capture_diagnostics = TRUE)
  expect_identical(chained$status, "success")
  stages <- chained$selection$state$stages
  evidence <- stages[[2]]$evidence
  rows <- evidence$fold_id != 1L
  upstream <- sb_fit_binary_model_selector(d[rows, c("a", "b", "c")], d$y[rows],
    "gbm_binary_importance", sb_binary_model_selector_options("gbm_binary_importance"))
  expect_setequal(evidence$fold_subsets[[1]][[length(evidence$fold_subsets[[1]])]],
    upstream$selected_features)
  expect_identical(sb_run_cv(task, chain, folds = 3, seed = 1)$status, "success")
})

test_that("SVM-RFE refits nested subsets and respects patient folds", {
  skip_if_not_installed("e1071")
  set.seed(19)
  x <- as.data.frame(matrix(rnorm(120 * 4), 120, 4))
  y <- rep(rep(c(0L, 1L), each = 30), each = 2)
  x$V1 <- x$V1 + 3 * y
  groups <- rep(seq_len(60), each = 2)
  fit <- sb_fit_svm_rfe(x, y, groups)
  e <- fit$evidence
  expect_equal(vapply(split(e$fold_id, groups), function(z) length(unique(z)), integer(1)),
    setNames(rep(1L, 60), as.character(seq_len(60))))
  expect_equal(e$curve$accuracy, colSums(e$correct) / nrow(x))
  expect_equal(e$curve$error, 1 - e$curve$accuracy)
  expect_length(fit$selected_features, which.max(e$curve$accuracy))
  for (k in 2:4) expect_equal(setdiff(e$subsets[[k - 1]], e$subsets[[k]]), character())
  train <- e$fold_id != 1
  replay <- sb_svm_rfe_path(x[train, ], y[train], x[!train, ])
  expect_identical(replay$subsets, e$fold_subsets[[1]])
})

# ---- selector-univ-binary-auc-top-k ----

manual_binary_auc <- function(outcome, score) {
  positive <- score[outcome == 1L]
  negative <- score[outcome == 0L]
  credits <- unlist(lapply(positive, function(pos) {
    (pos > negative) + 0.5 * (pos == negative)
  }))
  mean(credits)
}

test_that("binary AUC Top-K matches independently enumerated pairs", {
  outcome <- c(0L, 0L, 1L, 1L)
  x <- data.frame(
    higher = c(0, 1, 2, 3),
    lower = c(3, 2, 1, 0),
    tied = c(0, 1, 0, 1),
    constant = rep(4, 4),
    check.names = FALSE
  )
  fitted <- sb_fit_univ_binary_auc_top_k(x, outcome, n_features = 2L)
  scores <- fitted$scores

  expect_s3_class(fitted, "stormbreaker_selector_fit")
  expect_identical(fitted$id, "univ_binary_auc_top_k")
  expect_identical(fitted$source_features, names(x))
  expect_identical(fitted$selected_features, c("higher", "lower"))
  expect_identical(scores$feature, names(x))
  expect_equal(
    scores$auc[1:3],
    unname(vapply(
      x[1:3], manual_binary_auc, numeric(1), outcome = outcome
    )),
    tolerance = 0
  )
  expect_equal(scores$auc[1:3], c(1, 0, 0.5), tolerance = 0)
  expect_equal(scores$discrimination[1:3], c(1, 1, 0.5), tolerance = 0)
  expect_identical(
    scores$direction[1:3],
    c("higher_in_positive", "lower_in_positive", "tied")
  )
  expect_identical(scores$comparable_pairs[1:3], rep(4, 3L))
  expect_identical(scores$rank[1:3], 1:3)
  expect_identical(scores$selected, c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(scores$status[1:3], rep("success", 3L))
  expect_identical(scores$status[4], "failed")
  expect_identical(
    scores$error_class[4],
    "stormbreaker_selector_constant_feature"
  )
})

test_that("binary AUC Top-K applies only its fitted columns in fitted order", {
  x <- data.frame(
    reverse = c(4, 3, 2, 1, 0, -1),
    weak = c(0, 1, 1, 0, 1, 0),
    direct = c(-1, 0, 1, 2, 3, 4),
    check.names = FALSE
  )
  outcome <- c(0L, 0L, 0L, 1L, 1L, 1L)
  adapter <- sb_univ_binary_auc_top_k_selector_adapter(2L)
  fitted <- adapter$fit(x, outcome)
  applied <- adapter$apply(fitted, x)

  expect_identical(adapter$id, "univ_binary_auc_top_k")
  expect_identical(adapter$params, list(n_features = 2L))
  expect_identical(adapter$needs_outcome, TRUE)
  expect_identical(names(applied), fitted$selected_features)
  expect_identical(applied, x[fitted$selected_features])

  schema_error <- capture_sb_error(
    adapter$apply(fitted, x[c("weak", "reverse", "direct")])
  )
  expect_s3_class(schema_error, "stormbreaker_selector_schema_mismatch")
})

test_that("complementary AUC directions break ties by original order", {
  outcome <- c(1L, 0L, 0L, 0L)
  first <- c(3, 1, 2, 4)
  x <- data.frame(
    first = first,
    second = -first,
    check.names = FALSE
  )

  fitted <- sb_fit_univ_binary_auc_top_k(x, outcome, n_features = 1L)

  expect_equal(fitted$scores$auc, c(2 / 3, 1 / 3), tolerance = 1e-15)
  expect_identical(
    fitted$scores$discrimination[[1L]],
    fitted$scores$discrimination[[2L]]
  )
  expect_identical(fitted$scores$rank, 1:2)
  expect_identical(fitted$selected_features, "first")
})

test_that("binary AUC Top-K fails when fewer than K scores are valid", {
  x <- data.frame(
    informative = c(0, 1, 2, 3, 4, 5),
    constant_a = rep(1, 6),
    constant_b = rep(2, 6)
  )
  outcome <- c(0L, 0L, 0L, 1L, 1L, 1L)

  error <- capture_sb_error(
    sb_fit_univ_binary_auc_top_k(x, outcome, n_features = 2L)
  )

  expect_s3_class(error, "stormbreaker_insufficient_valid_features")
  expect_match(conditionMessage(error), "only 1 features")
})

test_that("binary AUC Top-K resolves through the admitted runtime", {
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 2L
  )
  runtime <- sb_resolve_pipeline_runtime(pipeline)

  expect_identical(runtime$id, pipeline$id)
  expect_identical(runtime$selector$id, pipeline$selector)
  expect_identical(runtime$selector$params, pipeline$selector_params)
  expect_type(runtime$selector$fit, "closure")
  expect_type(runtime$selector$apply, "closure")
  expect_identical(runtime$learner$id, "logistic")
  expect_identical(sb_validate_pipeline_runtime(runtime, pipeline), runtime)
})

# ---- selector-vsurf ----

test_that("VSURF saved evidence plots integrate into a separate clinician group", {
  skip_if_not_installed("VSURF")
  fixture <- withr::with_seed(9211L, {
    y <- sample(rep(c(0L, 1L), length.out = 72L))
    x <- matrix(rnorm(72L * 30L), nrow = 72L)
    x[, 1] <- x[, 1] + 1.4 * y
    x[, 2] <- x[, 2] - y
    x[, 3] <- x[, 3] + .7 * y
    x <- as.data.frame(x)
    names(x) <- paste0("gene_", seq_len(30L))
    list(x = x, y = y)
  })
  state <- sb_fit_vsurf_binary_selector(fixture$x, fixture$y)
  before <- serialize(state, NULL, version = 2L)
  importance <- sb_selector_diagnostics(state, "vsurf_importance")
  stages <- sb_selector_diagnostics(state, "vsurf_stages")
  expect_identical(importance$table, state$scores)
  expect_identical(stages$table, state$selection_history)
  expect_identical(importance$selected_features, state$selected_features)
  expect_equal(nrow(importance$table), 30L)
  ip <- sb_plot_selector_diagnostic(importance)
  sp <- sb_plot_selector_diagnostic(stages)
  expect_equal(ggplot2::ggplot_build(ip)$data[[2]]$x,
    state$scores$importance - state$scores$importance_sd)
  expect_equal(ggplot2::ggplot_build(ip)$data[[2]]$xend,
    state$scores$importance + state$scores$importance_sd)
  expect_identical(levels(ip$data$feature), rev(state$scores$feature[order(state$scores$rank)]))
  expect_equal(ggplot2::ggplot_build(sp)$data[[1]]$y, state$selection_history$n_features)
  result <- make_sb_curve_view_fixture("binary")$cv
  views <- sb_clinician_views(result, selectors = list(Importance = importance, Stages = stages))
  expect_identical(views$groups[grepl("^Selector:", names(views$plots))], rep("Feature selection", 2))
  expect_identical(views$status$status[grepl("^Selector:", views$status$view)], rep("available", 2))
  expect_null(views$plots[["Selector: Importance"]]$labels$caption)
  expect_match(views$captions[match("Selector: Importance", names(views$plots))], "between-forest SD", fixed = TRUE)
  file <- tempfile(fileext = ".html")
  sb_render_clinician_report(views, file)
  html <- paste(readLines(file), collapse = "\n")
  expect_match(html, "Feature selection", fixed = TRUE)
  expect_match(html, "between-forest SD", fixed = TRUE)
  expect_identical(serialize(state, NULL, version = 2L), before)
})

# ---- selector-xgboost ----

test_that("XGBoost selectors retain exactly native positive-gain features in both tasks", {
  skip_if_not_installed("xgboost")
  d <- sb_example_survival(150, 82)
  features <- c("gene_a", "gene_b", "gene_c")
  x <- d[features]; names(x) <- c("gene A", "gene-B", "gene C")
  d$outcome <- d$event
  original <- sb_fit_xgboost
  seen <- list()
  local_mocked_bindings(sb_fit_xgboost = function(x, endpoint, learner, params) {
    seen[[length(seen) + 1L]] <<- d$patient_id[match(x[[1]], d$gene_a)]
    original(x, endpoint, learner, params)
  })
  for (surv in c(FALSE, TRUE)) {
    id <- if (surv) "xgboost_cox_importance" else "xgboost_binary_importance"
    learner <- if (surv) "xgboost_cox" else "xgboost_binary"
    params <- if (surv) sb_survival_boosting_selector_params(id) else sb_binary_model_selector_options(id)
    a <- sb_resolve_selector_adapter(id, params)
    set.seed(11); before <- .Random.seed
    s <- if (surv) a$fit(x, d$time, d$event) else a$fit(x, d$outcome)
    expect_identical(.Random.seed, before)
    work <- x; names(work) <- paste0("feature_", 1:3)
    native <- sb_fit_xgboost(work, if (surv) list(time = d$time, event = d$event) else d$outcome,
      learner, sb_xgboost_params(learner))
    imp <- xgboost::xgb.importance(feature_names = names(work), model = native$booster)
    expect_setequal(s$selected_features, names(x)[match(imp$Feature[imp$Gain > 0], names(work))])
    expect_null(s$engine)
    expect_identical(a$apply(unserialize(serialize(s, NULL)), x), x[s$selected_features])
    task <- if (surv) sb_surv_task(d, "patient_id", "time", "event", features, "month") else
      sb_binary_task(d, "patient_id", "outcome", features)
    p <- sb_pipeline(id, if (surv) "glmnet_cox" else "glmnet_logistic", inner_folds = 3, inner_seed = 7)
    seen <- list()
    cv <- sb_run_cv(task, p, folds = 3, seed = 82)
    expect_identical(cv$status, "success")
    expect_gt(length(seen), 3L)
    expect_identical(any(lengths(seen) < 100L), TRUE)
    for (ids in seen) expect_identical(any(vapply(cv$folds,
      function(fold) all(ids %in% fold$analysis_ids), logical(1))), TRUE)
    expect_identical(cv$provenance$xgboost_version, as.character(packageVersion("xgboost")))
    fit <- sb_fit_final(task, p)
    expect_identical(fit$status, "success")
    expect_identical(fit$provenance$xgboost_version, as.character(packageVersion("xgboost")))
    pred <- sb_predict_final(unserialize(serialize(fit, NULL)), d[c("patient_id", features)], "patient_id")
    expect_identical(pred$status, "success")
  }
})
})

# ---- selector-specification ----
local({
test_that("fixed feature pipelines retain an explicit feature block", {
  binary <- sb_pipeline(
    selector = "fixed_features",
    learner = "logistic",
    selector_features = c("age", "gene_a"),
    feature_set_id = "clinical_omics"
  )
  survival <- sb_pipeline(
    selector = "fixed_features",
    learner = "coxph",
    selector_features = c("gene_b", "gene_a"),
    feature_set_id = "two_gene"
  )

  expect_identical(binary$id, "fixed_features__logistic__setclinical_omics")
  expect_identical(binary$task_type, "binary")
  expect_identical(binary$selector_params, list(
    features = c("age", "gene_a"),
    feature_set_id = "clinical_omics",
    n_features = 2L
  ))
  expect_identical(survival$task_type, "survival")
  expect_identical(
    sb_lint_pipeline(binary, task_type = "binary")$status,
    "admitted"
  )

  registry <- sb_pipeline_registry()
  fixed <- registry[registry$selector == "fixed_features", , drop = FALSE]
  expect_identical(nrow(fixed), length(unique(registry$learner)))
  for (task_type in c("binary", "survival")) {
    expect_identical(
      sum(fixed$task_type == task_type),
      length(unique(registry$learner[registry$task_type == task_type]))
    )
  }
})

test_that("fixed feature arguments preserve the existing positional API", {
  positional <- sb_pipeline("none", "glmnet_cox", NULL, 3L, 101L)
  named <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    selector_n_features = NULL,
    inner_folds = 3L,
    inner_seed = 101L
  )

  expect_identical(positional, named)
  expect_identical(
    positional$id,
    "none__glmnet_cox__alpha1__inner3__seed101"
  )

  ranger_positional <- do.call(
    sb_pipeline,
    c(
      list("none", "ranger_survival"),
      rep(list(NULL), 15L),
      list(100L, 2L, 3L, 7L)
    )
  )
  ranger_named <- sb_pipeline(
    learner = "ranger_survival",
    ranger_num_trees = 100L,
    ranger_mtry = 2L,
    ranger_min_node_size = 3L,
    ranger_seed = 7L
  )

  expect_identical(ranger_positional, ranger_named)
  expect_identical(
    ranger_positional$id,
    "none__ranger_survival__trees100__mtry2__minnode3__seed7"
  )

  gbm_positional <- do.call(
    sb_pipeline,
    c(
      list("none", "gbm_cox"),
      rep(list(NULL), 25L),
      list(100L, 2L, 10L, 0.1, 1, 1L)
    )
  )
  gbm_named <- sb_pipeline(
    learner = "gbm_cox",
    gbm_n_trees = 100L,
    gbm_interaction_depth = 2L,
    gbm_min_obs_node = 10L,
    gbm_shrinkage = 0.1,
    gbm_bag_fraction = 1,
    gbm_seed = 1L
  )

  expect_identical(gbm_positional, gbm_named)
  expect_identical(
    gbm_positional$id,
    paste0(
      "none__gbm_cox__trees100__depth2__minnode10__",
      "shrink0.1__bag1__seed1"
    )
  )
})

test_that("fixed feature specifications reject incomplete or ambiguous input", {
  missing_features <- capture_sb_error(sb_pipeline(
    selector = "fixed_features",
    learner = "logistic",
    feature_set_id = "set_a"
  ))
  missing_id <- capture_sb_error(sb_pipeline(
    selector = "fixed_features",
    learner = "logistic",
    selector_features = c("gene_a", "gene_b")
  ))
  duplicate_features <- capture_sb_error(sb_pipeline(
    selector = "fixed_features",
    learner = "logistic",
    selector_features = c("gene_a", "gene_a"),
    feature_set_id = "set_a"
  ))
  blank_id <- capture_sb_error(sb_pipeline(
    selector = "fixed_features",
    learner = "logistic",
    selector_features = c("gene_a", "gene_b"),
    feature_set_id = " "
  ))
  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    selector_features = c("gene_a", "gene_b"),
    feature_set_id = "set_a"
  ))
  one_feature_glmnet <- capture_sb_error(sb_pipeline(
    selector = "fixed_features",
    learner = "glmnet_logistic",
    selector_features = "gene_a",
    feature_set_id = "one_gene",
    inner_folds = 3L,
    inner_seed = 1L
  ))

  expect_s3_class(
    missing_features,
    "stormbreaker_missing_selector_parameter"
  )
  expect_s3_class(missing_id, "stormbreaker_missing_selector_parameter")
  expect_s3_class(
    duplicate_features,
    "stormbreaker_bad_fixed_features"
  )
  expect_s3_class(blank_id, "stormbreaker_bad_feature_set_id")
  expect_s3_class(
    irrelevant,
    "stormbreaker_irrelevant_selector_parameter"
  )
  expect_s3_class(
    one_feature_glmnet,
    "stormbreaker_glmnet_too_few_selected_features"
  )
})

test_that("fixed feature selection is outcome-blind and order-preserving", {
  x <- data.frame(
    gene_a = c(1, 2, 3),
    gene_b = c(4, 5, 6),
    gene_c = c(7, 8, 9)
  )
  adapter <- sb_resolve_selector_adapter(
    "fixed_features",
    list(
      features = c("gene_c", "gene_a"),
      feature_set_id = "ordered",
      n_features = 2L
    )
  )
  first <- sb_fit_runtime_selector(
    adapter,
    x,
    list(task_type = "binary", outcome = c(0L, 1L, 0L))
  )
  changed <- sb_fit_runtime_selector(
    adapter,
    x,
    list(task_type = "binary", outcome = c(1L, 0L, 1L))
  )
  selected <- adapter$apply(first, x)

  expect_identical(adapter$needs_outcome, FALSE)
  expect_identical(changed, first)
  expect_identical(first$feature_set_id, "ordered")
  expect_identical(first$selected_features, c("gene_c", "gene_a"))
  expect_identical(names(selected), c("gene_c", "gene_a"))
  expect_identical(selected$gene_c, x$gene_c)
  expect_s3_class(
    capture_sb_error(adapter$fit(x[c("gene_a", "gene_b")])),
    "stormbreaker_missing_fixed_feature"
  )
  expect_s3_class(
    capture_sb_error(adapter$apply(first, x[c("gene_b", "gene_a", "gene_c")])),
    "stormbreaker_selector_schema_mismatch"
  )
})

test_that("fixed feature blocks compare on identical outer folds", {
  data <- make_sb_binary_fixture(n = 60L, seed = 9961L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(
      selector = "fixed_features",
      learner = "logistic",
      selector_features = c("gene_a", "gene_b"),
      feature_set_id = "omics_pair"
    ),
    sb_pipeline(
      selector = "fixed_features",
      learner = "logistic",
      selector_features = "gene_c",
      feature_set_id = "single_numeric"
    )
  )

  observed <- sb_compare_cv(task, pipelines, folds = 3L, seed = 9962L)
  assignments <- lapply(observed$results, function(result) {
    result$predictions[c("patient_id", "fold")]
  })

  expect_identical(observed$status, "success")
  expect_identical(assignments[[1L]], assignments[[2L]])
  expect_identical(
    observed$comparison$selection_stability_status,
    rep("not_applicable", 2L)
  )
  expect_identical(
    observed$comparison$n_successful_selection_folds,
    rep(0L, 2L)
  )
  expect_identical(
    lapply(observed$results[[1L]]$selections, function(selection) {
      selection$selected_features
    }),
    rep(list(c("gene_a", "gene_b")), 3L)
  )
  expect_identical(
    lapply(observed$results[[2L]]$selections, function(selection) {
      selection$selected_features
    }),
    rep(list("gene_c"), 3L)
  )
})

test_that("fixed subsets feed final binary and survival prediction", {
  binary_data <- make_sb_binary_fixture(n = 60L, seed = 9971L)
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  binary_pipeline <- sb_pipeline(
    selector = "fixed_features",
    learner = "logistic",
    selector_features = c("gene_c", "gene_a"),
    feature_set_id = "binary_pair"
  )
  binary_fit <- sb_fit_final(binary_task, binary_pipeline)
  binary_new <- binary_data[1:6, c("patient_id", "gene_c", "gene_a")]
  binary_prediction <- sb_predict_final(
    binary_fit,
    binary_new,
    "patient_id"
  )
  binary_direct_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    c("gene_c", "gene_a")
  )
  binary_direct_fit <- sb_fit_final(
    binary_direct_task,
    sb_pipeline(learner = "logistic")
  )
  binary_direct_prediction <- sb_predict_final(
    binary_direct_fit,
    binary_new,
    "patient_id"
  )

  survival_data <- make_sb_fixture(n = 60L, seed = 9972L)
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  survival_pipeline <- sb_pipeline(
    selector = "fixed_features",
    learner = "coxph",
    selector_features = c("gene_b", "gene_a"),
    feature_set_id = "survival_pair"
  )
  survival_fit <- sb_fit_final(survival_task, survival_pipeline)
  survival_new <- survival_data[1:6, c("patient_id", "gene_b", "gene_a")]
  survival_prediction <- sb_predict_final(
    survival_fit,
    survival_new,
    "patient_id"
  )
  survival_direct_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    c("gene_b", "gene_a"),
    "month"
  )
  survival_direct_fit <- sb_fit_final(
    survival_direct_task,
    sb_pipeline(learner = "coxph")
  )
  survival_direct_prediction <- sb_predict_final(
    survival_direct_fit,
    survival_new,
    "patient_id"
  )
  survival_cv <- sb_run_cv(
    survival_task,
    survival_pipeline,
    folds = 3L,
    seed = 9973L
  )

  expect_identical(binary_fit$status, "success")
  expect_identical(binary_fit$required_features, c("gene_c", "gene_a"))
  expect_identical(binary_prediction$status, "success")
  expect_equal(
    binary_prediction$predictions$probability,
    binary_direct_prediction$predictions$probability,
    tolerance = 0
  )
  expect_identical(
    binary_prediction$predictions$patient_id,
    binary_new$patient_id
  )
  expect_identical(survival_fit$status, "success")
  expect_identical(
    survival_fit$required_features,
    c("gene_b", "gene_a")
  )
  expect_identical(survival_prediction$status, "success")
  expect_equal(
    survival_prediction$predictions$risk,
    survival_direct_prediction$predictions$risk,
    tolerance = 0
  )
  expect_identical(
    survival_prediction$predictions$patient_id,
    survival_new$patient_id
  )
  expect_identical(survival_cv$status, "success")
  expect_identical(
    survival_cv$selection_stability$status,
    "not_applicable"
  )
  expect_identical(
    lapply(survival_cv$selections, function(selection) {
      selection$selected_features
    }),
    rep(list(c("gene_b", "gene_a")), 3L)
  )
})

test_that("fixed glmnet subsets do not masquerade as learned inner selection", {
  skip_if_not_installed("glmnet")

  binary_data <- make_sb_binary_fixture(n = 60L, seed = 9981L)
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  binary_pipeline <- sb_pipeline(
    selector = "fixed_features",
    learner = "glmnet_logistic",
    selector_features = c("gene_b", "gene_a"),
    feature_set_id = "binary_glmnet_pair",
    glmnet_alpha = 0,
    inner_folds = 3L,
    inner_seed = 9982L
  )
  binary_fit <- sb_fit_final(binary_task, binary_pipeline)

  survival_data <- make_sb_fixture(n = 60L, seed = 9983L)
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  survival_pipeline <- sb_pipeline(
    selector = "fixed_features",
    learner = "glmnet_cox",
    selector_features = c("gene_c", "gene_a"),
    feature_set_id = "survival_glmnet_pair",
    glmnet_alpha = 0,
    inner_folds = 3L,
    inner_seed = 9984L
  )
  survival_fit <- sb_fit_final(survival_task, survival_pipeline)

  expect_identical(binary_fit$status, "success")
  expect_identical(
    binary_fit$tuning$inner_selection_method,
    "predeclared_feature_subset"
  )
  expect_length(binary_fit$tuning$inner_selections, 0L)
  expect_identical(survival_fit$status, "success")
  expect_identical(
    survival_fit$tuning$inner_selection_method,
    "predeclared_feature_subset"
  )
  expect_length(survival_fit$tuning$inner_selections, 0L)
})
})

# ---- selector-survival ----
local({
# ---- selector-mime-survival ----

make_mime_selector_fixture <- function(n = 90L, seed = 8601L) {
  set.seed(seed)
  signal_a <- stats::rnorm(n)
  signal_b <- stats::rnorm(n)
  signal_c <- stats::rnorm(n)
  noise_a <- stats::rnorm(n)
  noise_b <- stats::rnorm(n)
  noise_c <- stats::rnorm(n)
  linear_predictor <- 1.2 * signal_a - 0.9 * signal_b + 0.6 * signal_c
  event_time <- stats::rexp(n, rate = exp(linear_predictor) / 15)
  censor_time <- stats::rexp(n, rate = 1 / 20)
  data.frame(
    patient_id = sprintf("mime_selector_%03d", seq_len(n)),
    time = pmax(pmin(event_time, censor_time), 0.05),
    event = as.integer(event_time <= censor_time),
    signal_a = signal_a,
    signal_b = signal_b,
    signal_c = signal_c,
    noise_a = noise_a,
    noise_b = noise_b,
    noise_c = noise_c,
    stringsAsFactors = FALSE
  )
}

test_that("MIME selection failures report how many predictors were returned", {
  source_features <- c("gene_a", "gene_b", "gene_c")
  for (selected in list(character(), NA_character_, c("gene_a", "gene_a"))) {
    error <- capture_sb_error(
      sb_validate_mime_selector_features(selected, source_features, "Lasso")
    )
    expect_s3_class(error, "stormbreaker_mime_selector_too_few_features")
    expect_match(
      conditionMessage(error),
      paste0("returned ", length(selected), " of 3 candidate predictors"),
      fixed = TRUE
    )
  }
})

test_that("Lasso can retain one predictor for a single-input survival learner", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("randomForestSRC")
  d <- sb_example_survival(n = 200L, seed = 20260905L)
  noise <- withr::with_seed(20260906L, matrix(rnorm(200L * 17L), nrow = 200L))
  for (j in seq_len(ncol(noise))) d[[paste0("noise_", j)]] <- noise[, j]
  features <- setdiff(names(d), c("patient_id", "time", "event"))
  task <- sb_surv_task(d[1:160, ], "patient_id", "time", "event", features,
                       "month", prediction_horizon = 6)
  pipeline <- sb_pipeline(
    selector = "lasso_nonzero", learner = "rsf_survival",
    selector_inner_folds = 3L, selector_seed = 20260905L
  )
  fit <- sb_fit_final(task, pipeline)
  expect_identical(fit$status, "success")
  expect_identical(fit$selection$selected_features, "gene_a")
  prediction <- sb_predict_final(fit, d[161:200, c("patient_id", features)], "patient_id")
  expect_identical(prediction$status, "success")
  expect_equal(nrow(prediction$predictions), 40L)
  expect_equal(sum(is.finite(prediction$predictions$risk)), 40L)
  expect_identical(prediction$predictions$patient_id, d$patient_id[161:200])
})

test_that("RSF minimal-depth selection uses native top variables", {
  skip_if_not_installed("randomForestSRC")
  data <- make_mime_selector_fixture()
  x <- data[4:9]
  adapter <- sb_rsf_minimal_depth_selector_adapter(
    nodesize = 5L,
    seed = 8602L
  )

  fitted <- adapter$fit(x, data$time, data$event)
  selected <- adapter$apply(fitted, x[rev(seq_len(nrow(x))), ])

  expect_s3_class(fitted, "stormbreaker_selector_fit")
  expect_identical(fitted$id, "rsf_minimal_depth")
  expect_identical(fitted$source_features, names(x))
  expect_gt(length(fitted$selected_features), 1L)
  expect_true(all(fitted$selected_features %in% names(x)))
  expect_identical(names(selected), fitted$selected_features)
  expect_identical(
    selected,
    x[rev(seq_len(nrow(x))), fitted$selected_features, drop = FALSE]
  )
})

test_that("CoxBoost nonzero selection records tuned native coefficients", {
  skip_if_not_installed("CoxBoost")
  data <- make_mime_selector_fixture(seed = 8611L)
  x <- data[4:9]
  adapter <- sb_coxboost_nonzero_selector_adapter(
    inner_folds = 3L,
    seed = 8612L
  )

  fitted <- adapter$fit(x, data$time, data$event)

  expect_s3_class(fitted, "stormbreaker_selector_fit")
  expect_identical(fitted$id, "coxboost_nonzero")
  expect_gt(length(fitted$selected_features), 1L)
  expect_true(all(fitted$selected_features %in% names(x)))
  expect_true(is.finite(fitted$penalty))
  expect_true(is.finite(fitted$optimal_step))
  expect_identical(
    fitted$scores$selected,
    fitted$scores$coefficient != 0
  )
})

test_that("Lasso nonzero selection records fold-local lambda", {
  skip_if_not_installed("glmnet")
  data <- make_mime_selector_fixture(seed = 8621L)
  x <- data[4:9]
  adapter <- sb_lasso_nonzero_selector_adapter(
    inner_folds = 3L,
    seed = 8622L
  )

  fitted <- adapter$fit(x, data$time, data$event)

  expect_s3_class(fitted, "stormbreaker_selector_fit")
  expect_identical(fitted$id, "lasso_nonzero")
  expect_gt(length(fitted$selected_features), 1L)
  expect_true(all(fitted$selected_features %in% names(x)))
  expect_true(is.finite(fitted$lambda))
  expect_identical(fitted$alpha, 1)
  expect_identical(fitted$selection_rule, "lambda.min")
  expect_identical(
    fitted$scores$selected,
    fitted$scores$coefficient != 0
  )
})

test_that("sequential selection fits StepCox only on stage-one features", {
  skip_if_not_installed("randomForestSRC")
  data <- make_mime_selector_fixture(n = 120L, seed = 8631L)
  x <- data[4:9]
  adapter <- sb_rsf_then_stepcox_selector_adapter(
    nodesize = 5L,
    seed = 8632L,
    direction = "forward"
  )

  fitted <- adapter$fit(x, data$time, data$event)
  selected <- adapter$apply(fitted, x)

  expect_s3_class(fitted, "stormbreaker_selector_fit")
  expect_identical(fitted$id, "rsf_then_stepcox")
  expect_identical(
    fitted$stages$second$source_features,
    fitted$stages$first$selected_features
  )
  expect_true(all(
    fitted$selected_features %in% fitted$stages$first$selected_features
  ))
  expect_identical(names(selected), fitted$selected_features)
  expect_identical(fitted$stages$second$direction, "forward")
})

test_that("all MIME first-stage selectors compose with StepCox", {
  adapters <- list(
    rsf_then_stepcox = sb_rsf_then_stepcox_selector_adapter(
      nodesize = 5L,
      seed = 8641L,
      direction = "both"
    ),
    coxboost_then_stepcox = sb_coxboost_then_stepcox_selector_adapter(
      inner_folds = 3L,
      seed = 8642L,
      direction = "backward"
    ),
    lasso_then_stepcox = sb_lasso_then_stepcox_selector_adapter(
      inner_folds = 3L,
      seed = 8643L,
      direction = "forward"
    )
  )

  expect_identical(
    unname(vapply(adapters, `[[`, character(1), "id")),
    names(adapters)
  )
  expect_true(all(vapply(adapters, `[[`, logical(1), "needs_outcome")))
  expect_true(all(vapply(adapters, function(adapter) {
    is.function(adapter$fit) && is.function(adapter$apply)
  }, logical(1))))
})

test_that("MIME survival selector controls reconstruct registered pipelines", {
  pipelines <- list(
    sb_pipeline(
      selector = "rsf_minimal_depth",
      learner = "coxboost_cv",
      selector_rsf_nodesize = 5L,
      selector_seed = 8651L,
      learner_tuning_folds = 3L,
      learner_seed = 8652L
    ),
    sb_pipeline(
      selector = "coxboost_nonzero",
      learner = "glmnet_cox",
      selector_inner_folds = 3L,
      selector_seed = 8652L,
      inner_folds = 3L,
      inner_seed = 8653L,
      glmnet_alpha = 0.5
    ),
    sb_pipeline(
      selector = "lasso_then_stepcox",
      learner = "coxph",
      selector_inner_folds = 3L,
      selector_seed = 8654L,
      stepcox_direction = "both"
    )
  )

  expect_length(unique(vapply(pipelines, `[[`, character(1), "id")), 3L)
  for (pipeline in pipelines) {
    runtime <- sb_resolve_pipeline_runtime(pipeline)
    expect_identical(runtime$selector$params, pipeline$selector_params)
    expect_identical(
      sb_validate_pipeline_contract(pipeline, "survival"),
      pipeline
    )
  }

  missing_seed <- capture_sb_error(sb_pipeline(
    selector = "lasso_nonzero",
    learner = "coxboost_cv",
    selector_inner_folds = 3L
  ))
  irrelevant_control <- capture_sb_error(sb_pipeline(
    learner = "coxph",
    selector_seed = 8655L
  ))
  admitted_shape <- sb_pipeline(
    selector = "rsf_minimal_depth",
    learner = "bart_survival",
    selector_rsf_nodesize = 5L,
    selector_seed = 8656L
  )

  expect_s3_class(missing_seed, "stormbreaker_missing_selector_parameter")
  expect_s3_class(
    irrelevant_control,
    "stormbreaker_irrelevant_selector_parameter"
  )
  expect_s3_class(admitted_shape, "stormbreaker_pipeline")
  expect_identical(admitted_shape$selector, "rsf_minimal_depth")
  expect_identical(admitted_shape$learner, "bart_survival")
})

# ---- selector-stepcox ----

make_stepcox_fixture <- function(n = 120L, seed = 2401L) {
  set.seed(seed)
  signal_a <- stats::rnorm(n)
  signal_b <- stats::rnorm(n)
  noise <- stats::rnorm(n)
  event_time <- stats::rexp(
    n,
    rate = exp(1.5 * signal_a - 1.2 * signal_b)
  )
  censor_time <- stats::rexp(n, rate = 0.12)
  data.frame(
    patient_id = sprintf("stepcox_%03d", seq_len(n)),
    time = pmin(event_time, censor_time),
    event = as.integer(event_time <= censor_time),
    signal_a = signal_a,
    signal_b = signal_b,
    noise = noise,
    stringsAsFactors = FALSE
  )
}

test_that("StepCox directions select from analysis outcomes and preserve names", {
  data <- make_stepcox_fixture()
  x <- data[c("signal_a", "signal_b", "noise")]
  names(x) <- c("signal A", "signal-B", "noise")

  for (direction in c("forward", "backward", "both")) {
    adapter <- sb_stepcox_selector_adapter(direction)
    fitted <- adapter$fit(x, data$time, data$event)
    selected <- adapter$apply(fitted, x)

    expect_s3_class(fitted, "stormbreaker_selector_fit")
    expect_identical(fitted$id, "stepcox")
    expect_identical(fitted$direction, direction)
    expect_identical(fitted$source_features, names(x))
    expect_setequal(
      intersect(fitted$selected_features, c("signal A", "signal-B")),
      c("signal A", "signal-B")
    )
    expect_identical(names(selected), fitted$selected_features)
    expect_identical(selected, x[fitted$selected_features])
    expect_identical(adapter$params, list(direction = direction))
    expect_identical(adapter$needs_outcome, TRUE)
    expect_identical(fitted$scores$feature, names(x))
    expect_identical(
      fitted$scores$selected,
      fitted$scores$feature %in% fitted$selected_features
    )
    expect_true(all(is.finite(
      fitted$scores$coefficient[fitted$scores$selected]
    )))
    expect_true(is.finite(fitted$final_aic))
    expect_s3_class(fitted$selection_path, "data.frame")
    expect_false(any(grepl(
      ".stormbreaker_stepcox_",
      fitted$selection_path$Step,
      fixed = TRUE
    )))
  }

  fitted <- sb_stepcox_selector_adapter("forward")$fit(
    x,
    data$time,
    data$event
  )
  schema_error <- capture_sb_error(sb_apply_stepcox(
    fitted,
    x[c("noise", "signal A", "signal-B")]
  ))
  expect_s3_class(schema_error, "stormbreaker_selector_schema_mismatch")
})

test_that("StepCox pipeline parameters are explicit and task-specific", {
  pipelines <- lapply(c("forward", "backward", "both"), function(direction) {
    sb_pipeline(
      selector = "stepcox",
      learner = "coxph",
      stepcox_direction = direction
    )
  })

  expect_identical(
    vapply(pipelines, function(pipeline) pipeline$id, character(1)),
    paste0("stepcox__coxph__direction", c("forward", "backward", "both"))
  )
  expect_identical(
    lapply(pipelines, function(pipeline) pipeline$selector_params),
    lapply(
      c("forward", "backward", "both"),
      function(direction) list(direction = direction)
    )
  )

  missing <- capture_sb_error(sb_pipeline(
    selector = "stepcox",
    learner = "coxph"
  ))
  invalid <- capture_sb_error(sb_pipeline(
    selector = "stepcox",
    learner = "coxph",
    stepcox_direction = "sideways"
  ))
  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "coxph",
    stepcox_direction = "forward"
  ))
  incompatible <- capture_sb_error(sb_pipeline(
    selector = "stepcox",
    learner = "logistic",
    stepcox_direction = "forward"
  ))

  expect_s3_class(missing, "stormbreaker_missing_selector_parameter")
  expect_s3_class(invalid, "stormbreaker_bad_stepcox_direction")
  expect_s3_class(irrelevant, "stormbreaker_irrelevant_selector_parameter")
  expect_s3_class(incompatible, "stormbreaker_incompatible_binary_selector")
})

test_that("StepCox canonicalizes named and classed direction values", {
  data <- make_stepcox_fixture(n = 120L, seed = 2405L)
  x <- data[c("signal_a", "signal_b", "noise")]
  plain_direction <- "forward"
  named_direction <- c(requested = "forward")
  classed_direction <- structure(
    "forward",
    class = c("stormbreaker_test_direction", "character")
  )

  plain_pipeline <- sb_pipeline(
    selector = "stepcox",
    learner = "coxph",
    stepcox_direction = plain_direction
  )
  named_pipeline <- sb_pipeline(
    selector = "stepcox",
    learner = "coxph",
    stepcox_direction = named_direction
  )
  classed_pipeline <- sb_pipeline(
    selector = "stepcox",
    learner = "coxph",
    stepcox_direction = classed_direction
  )
  plain_fit <- sb_fit_stepcox(
    x,
    data$time,
    data$event,
    plain_direction
  )
  named_fit <- sb_fit_stepcox(
    x,
    data$time,
    data$event,
    named_direction
  )
  classed_fit <- sb_fit_stepcox(
    x,
    data$time,
    data$event,
    classed_direction
  )

  expect_identical(named_pipeline, plain_pipeline)
  expect_identical(classed_pipeline, plain_pipeline)
  expect_identical(named_fit, plain_fit)
  expect_identical(classed_fit, plain_fit)
  expect_identical(plain_fit$start_model, "null")
})

test_that("StepCox runs fold-locally and supports final Cox prediction", {
  data <- make_stepcox_fixture(n = 150L, seed = 2411L)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    c("signal_a", "signal_b", "noise"),
    "month"
  )
  pipeline <- sb_pipeline(
    selector = "stepcox",
    learner = "coxph",
    stepcox_direction = "forward"
  )

  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 2412L)
  fit <- sb_fit_final(task, pipeline)
  new_data <- data[1:12, c(
    "patient_id",
    fit$required_features
  )]
  prediction <- sb_predict_final(fit, new_data, "patient_id")

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(result$selection_stability$selector, "stepcox")
  expect_identical(result$selection_stability$status, "success")
  expect_true(all(vapply(result$selections, function(selection) {
    selection$status == "success" &&
      all(c("signal_a", "signal_b") %in% selection$selected_features)
  }, logical(1))))
  expect_identical(fit$status, "success")
  expect_identical(fit$selection$state$direction, "forward")
  expect_setequal(
    intersect(fit$required_features, c("signal_a", "signal_b")),
    c("signal_a", "signal_b")
  )
  expect_identical(prediction$status, "success")
  expect_identical(
    prediction$predictions$patient_id,
    new_data$patient_id
  )
  expect_true(all(is.finite(prediction$predictions$risk)))
})

test_that("outer assessment outcomes cannot change StepCox selection", {
  data <- make_stepcox_fixture(n = 120L, seed = 2415L)
  split <- list(
    fold = 1L,
    analysis_rows = 1:90,
    assessment_rows = 91:120,
    analysis_ids = data$patient_id[1:90],
    assessment_ids = data$patient_id[91:120]
  )
  pipeline <- sb_pipeline(
    selector = "stepcox",
    learner = "coxph",
    stepcox_direction = "forward"
  )
  make_task <- function(frame) {
    sb_surv_task(
      frame,
      "patient_id",
      "time",
      "event",
      c("signal_a", "signal_b", "noise"),
      "month"
    )
  }

  original <- sb_run_one_fold(
    make_task(data),
    pipeline,
    sb_resolve_pipeline_runtime(pipeline),
    split,
    seed = 2416L
  )
  changed <- data
  changed$time[91:120] <- changed$time[91:120] + 1000
  changed$event[91:120] <- 1L - changed$event[91:120]
  perturbed <- sb_run_one_fold(
    make_task(changed),
    pipeline,
    sb_resolve_pipeline_runtime(pipeline),
    split,
    seed = 2416L
  )

  expect_identical(
    original$selection$selected_features,
    perturbed$selection$selected_features
  )
  expect_identical(original$selection$scores, perturbed$selection$scores)
  expect_equal(
    original$predictions$risk,
    perturbed$predictions$risk,
    tolerance = 0
  )
})

test_that("StepCox is repeated inside glmnet inner analysis folds", {
  skip_if_not_installed("glmnet")

  data <- make_stepcox_fixture(n = 150L, seed = 2421L)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    c("signal_a", "signal_b", "noise"),
    "month"
  )
  pipeline <- sb_pipeline(
    selector = "stepcox",
    learner = "glmnet_cox",
    stepcox_direction = "backward",
    glmnet_alpha = 0,
    inner_folds = 3L,
    inner_seed = 2422L
  )

  fit <- sb_fit_final(task, pipeline)

  expect_identical(fit$status, "success")
  expect_identical(
    fit$tuning$inner_selection_method,
    "glmnet_exclude_function"
  )
  expect_length(fit$tuning$inner_selections, 3L)
  expect_true(all(vapply(fit$tuning$inner_selections, function(selection) {
    length(selection$analysis_patient_ids) < nrow(data) &&
      all(c("signal_a", "signal_b") %in% selection$selected_features)
  }, logical(1))))
})

# ---- selector-survival-boosting ----

test_that("survival boosting selectors fit, predict and expose their evidence", {
  d <- sb_example_survival(150, 82)
  t <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
  for (method in c("GBM", "glmBoost", "SuperPC", "GamBoost")) {
    p <- surv_pipeline("Cox", surv_selector(method))
    expect_equal(sb_validate_pipeline_contract(p, "survival"), p)
    f <- sb_fit_final(t, p, capture_diagnostics = TRUE)
    expect_identical(f$status, "success")
    s <- f$selection$state
    expect_identical(s$selected_features, s$scores$feature[s$scores$selected])
    expect_equal(nrow(s$scores), 3)
    if (method != "SuperPC") expect_gt(nrow(s$evidence$training_loss), 0)
    expect_identical(sb_predict_final(f, d, "patient_id")$status, "success")
    plots <- sb_plots(f)
    expect_length(plots$unavailable, 0)
    expect_s3_class(plots$plots$selector_1, "ggplot")
    if (method != "SuperPC") expect_s3_class(plots$plots$selector_1_training_loss, "ggplot")
    cv <- sb_run_cv(t, p, folds = 3, seed = 82)
    expect_identical(cv$status, "success")
  }
  chain <- sb_fit_final(t, surv_pipeline("Cox", surv_selector(c("GBM", "glmBoost"))))
  expect_identical(chain$status, "success")
  expect_length(chain$selection$state$stages, 2)
})

test_that("GamBoost selects native active spline terms, not linear glmBoost terms", {
  skip_if_not_installed("mboost")
  d <- sb_example_survival(150, 82)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene A", "gene-B", "gene C")
  p <- sb_pipeline("gamboost_cox_nonzero", "coxph")
  adapter <- sb_survival_boosting_selector_adapter(p$selector, p$selector_params)
  fitted <- adapter$fit(x, d$time, d$event)
  work <- x
  names(work) <- paste0("feature_", seq_len(ncol(x)))
  native <- sb_fit_mboost_cox(work, d$time, d$event,
    "gamboost_cox", p$selector_params)$engine
  beta <- stats::coef(native)
  active <- names(beta)[vapply(beta, function(b) any(b != 0), logical(1))]
  native_features <- sub("^bbs\\(([^,)]+).*", "\\1", active)
  expect_setequal(fitted$selected_features, names(x)[match(native_features, names(work))])
  expect_identical(adapter$apply(fitted, x), x[fitted$selected_features])
  expect_identical(unserialize(serialize(fitted, NULL)), fitted)
  expect_null(fitted$engine)
  expect_identical(sb_selector_is_learned(p$selector), TRUE)
})

test_that("GamBoost is refitted inside glmnet outer and inner analysis folds", {
  skip_if_not_installed("mboost")
  d <- sb_example_survival(150, 82)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month")
  original <- sb_fit_mboost_cox
  seen <- list()
  local_mocked_bindings(sb_fit_mboost_cox = function(x, time, event, learner, params) {
    seen[[length(seen) + 1L]] <<- d$patient_id[match(x[[1L]], d$gene_a)]
    original(x, time, event, learner, params)
  })
  p <- sb_pipeline("gamboost_cox_nonzero", "glmnet_cox", inner_folds = 3, inner_seed = 7)
  cv <- sb_run_cv(task, p, folds = 3, seed = 82)
  expect_identical(cv$status, "success")
  expect_gt(length(seen), 3L)
  expect_identical(all(lengths(seen) < nrow(d)), TRUE)
  expect_identical(any(lengths(seen) < 100L), TRUE)
  for (ids in seen) expect_identical(any(vapply(cv$folds,
    function(fold) all(ids %in% fold$analysis_ids), logical(1))), TRUE)
})
})

# ---- selector-views ----
local({
# ---- selector-captured-plots ----

test_that("captured selectors retain classic paths without changing selections", {
  withr::local_seed(81)
  x <- as.data.frame(matrix(rnorm(300 * 6), 300))
  names(x) <- paste0("gene", 1:6)
  y <- rbinom(300, 1, plogis(2 * x$gene1 - 2 * x$gene2))
  for (id in c("lasso_binary_nonzero", "enet_binary_nonzero", "rf_binary_importance")) {
    p <- sb_binary_model_selector_options(id)
    compact <- sb_fit_binary_model_selector(x, y, id, p)
    captured <- sb_captured_selector_adapter(id, p)$fit(x, y)
    expect_identical(compact$selected_features, captured$selected_features)
    expect_equal(compact$scores, captured$scores)
    if (id == "rf_binary_importance") {
      expect_equal(nrow(captured$evidence$oob_error), p$trees)
      expect_null(compact$evidence$oob_error)
      plot <- sb_selector_oob_plot(captured$evidence$oob_error)
    } else {
      path <- captured$evidence$coefficient_path
      expect_identical(rownames(path$beta), names(x))
      expect_gt(length(path$lambda), 10)
      expect_null(compact$evidence$coefficient_path)
      plot <- sb_selector_path_plot(path, captured$evidence$lambda)
    }
    expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
  }
})

test_that("Cox Lasso and Enet capture the fitted native path", {
  d <- sb_example_survival(180, 41)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  for (id in c("lasso_nonzero", "enet_cox_nonzero")) {
    s <- sb_captured_selector_adapter(id, list(inner_folds = 5L, seed = 1L))$fit(x, d$time, d$event)
    expect_gt(length(s$coefficient_path$lambda), 10)
    expect_identical(rownames(s$coefficient_path$beta), names(x))
    expect_s3_class(ggplot2::ggplotGrob(sb_selector_path_plot(s$coefficient_path, s$lambda)), "gtable")
  }
})

test_that("additional classic selector plots reach the doctor entrypoint", {
  withr::local_seed(81)
  d <- as.data.frame(matrix(rnorm(300 * 6), 300)); names(d) <- paste0("gene", 1:6)
  d$id <- paste0("p", 1:300); d$y <- rbinom(300, 1, plogis(2*d$gene1 - 2*d$gene2))
  task <- sb_binary_task(d, "id", "y", paste0("gene", 1:6))
  for (selector in c("ABESS", "StepGLM", "glmBoost")) {
    f <- sb_fit_final(task, bi_pipeline("RF", bi_selector(selector)), capture_diagnostics = TRUE)
    expect_identical(f$status, "success")
    z <- sb_plots(f)
    names <- grep("^selector_", names(z$plots), value = TRUE)
    expect_gte(length(names), 2)
    for (nm in names) expect_s3_class(ggplot2::ggplotGrob(z$plots[[nm]]), "gtable")
  }
  d <- sb_example_survival(180, 41)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
  f <- sb_fit_final(task, sb_pipeline("coxboost_then_stepcox", "coxph", stepcox_direction = "both",
    selector_seed = 1L, selector_inner_folds = 5L), capture_diagnostics = TRUE)
  expect_identical(f$status, "success")
  z <- sb_plots(f)
  for (nm in c("selector_1_cv", "selector_1_coefficient_path")) {
    expect_s3_class(z$plots[[nm]], "ggplot")
    expect_s3_class(ggplot2::ggplotGrob(z$plots[[nm]]), "gtable")
  }
})

# ---- selector-diagnostics ----

selector_diagnostic_fixture <- function(kind) {
  withr::with_seed(881L, {
    n <- 120L
    d <- data.frame(patient_id = paste0("s", seq_len(n)),
      gene_a = rnorm(n), gene_b = rnorm(n), gene_c = rnorm(n),
      gene_d = rnorm(n), gene_e = rnorm(n), gene_f = rnorm(n))
    if (kind == "binary") {
      d$outcome <- rbinom(n, 1, plogis(2 * d$gene_a - 1.3 * d$gene_b))
      sb_binary_task(d, "patient_id", "outcome", paste0("gene_", letters[1:6]))
    } else {
      times <- rexp(n, rate = .08 * exp(1.8 * d$gene_a - 1.3 * d$gene_b)) + .1
      d$time <- pmin(times, median(times))
      d$event <- as.integer(times <= median(times))
      sb_surv_task(d, "patient_id", "time", "event", paste0("gene_", letters[1:6]), "month")
    }
  })
}

test_that("Boruta capture preserves selection and real predictor/shadow histories", {
  skip_if_not_installed("Boruta")
  skip_if_not_installed("ranger")
  for (kind in c("binary", "survival")) {
    task <- selector_diagnostic_fixture(kind)
    pipeline <- sb_pipeline(selector = paste0("boruta_", kind, "_roughfix"),
      learner = if (kind == "binary") "logistic" else "coxph")
    plain <- sb_fit_final(task, pipeline)
    captured <- sb_fit_final(task, pipeline, capture_diagnostics = TRUE)
    expect_identical(plain$status, "success")
    expect_identical(captured$status, "success")
    expect_identical(plain$learner_state, captured$learner_state)
    expect_identical(plain$selection$scores, captured$selection$scores)
    expect_null(plain$selection$state$diagnostic_history)
    q <- sb_selector_diagnostics(captured, "boruta_summary")
    expect_identical(q$table, captured$selection$scores)
    h <- sb_selector_diagnostics(captured, "boruta_history")
    native <- captured$selection$state$diagnostic_history$importance
    expect_identical(h$table$raw_importance, as.vector(native))
    expect_equal(sum(h$table$not_evaluated_after_rejection), sum(native == -Inf, na.rm = TRUE))
    expect_equal(nrow(h$table), length(native))
    expect_identical(sort(unique(h$table$kind)), c("Predictor", "Shadow reference"))
    err <- tryCatch(sb_selector_diagnostics(plain, "boruta_history"), error = identity)
    expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      for (type in c("boruta_summary", "boruta_hits", "boruta_history", "boruta_boxplot")) {
        expect_s3_class(ggplot2::ggplot_build(sb_plot_selector_diagnostic(
          sb_selector_diagnostics(captured, type))), "ggplot_built")
      }
    }
  }
})

test_that("varSelRF exposes saved OOB sequence without relabelling it CV", {
  skip_if_not_installed("varSelRF")
  fit <- sb_fit_final(selector_diagnostic_fixture("binary"),
    sb_pipeline(selector = "varselrf_binary_oob_1se", learner = "logistic"))
  expect_identical(fit$status, "success")
  q <- sb_selector_diagnostics(fit, "selection_history")
  expect_identical(q$table, fit$selection$state$selection_history)
  expect_identical(q$selected_features, fit$required_features)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(ggplot2::ggplot_build(sb_plot_selector_diagnostic(q)), "ggplot_built")
  }
})

# ---- selector-score-plots ----

test_that("saved score plots preserve native values and missing evidence", {
  ids <- c("rsf_minimal_depth", "sis_binary_nlogn", "mrmr_binary_nlogn",
    "jmi_binary_nlogn", "cmim_binary_nlogn", "disr_binary_nlogn", "jmim_binary_nlogn",
    "lpc_binary_nlogn", "lpc_cox_nlogn", "stabsel_binary_lasso", "stabsel_cox_lasso")
  for (id in ids) for (metric in names(sb_selector_score_views(id))) {
    scores <- data.frame(feature = c("a", "b", "c"), selected = c(TRUE, FALSE, FALSE))
    scores[[metric]] <- c(.9, .2, NA_real_)
    if (metric == "rank") scores[[metric]] <- c(1L, 2L, NA_integer_)
    s <- structure(list(id = id, scores = scores, params = list(cutoff = .75)), class = "stormbreaker_selector_fit")
    p <- sb_plot_selector_scores(s, metric)
    expect_identical(attr(p, "selector_scores"), scores)
    expect_setequal(p$data$value, scores[[metric]][1:2])
    expect_equal(nrow(p$data), 2)
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
    if (metric == "selection_probability") expect_equal(p$layers[[2]]$data$xintercept, .75)
  }
})

test_that("doctor plots expose screening views without refitting", {
  d <- sb_example_binary(100, 11)
  t <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  f <- sb_fit_final(t, bi_pipeline("RF", bi_selector("mRMR")), capture_diagnostics = TRUE)
  expect_identical(f$status, "success")
  before <- f$selection$state
  p <- sb_plots(f)
  expect_s3_class(p$plots$selector_1_mrmr_score, "ggplot")
  expect_identical(f$selection$state, before)
  expect_null(p$unavailable$selector_1)
})

test_that("univariable selector plots expose the saved ranking statistics", {
  fits <- list()
  for (type in c("binary", "survival")) {
    d <- if (type == "binary") sb_example_binary(100, 21) else sb_example_survival(100, 21)
    genes <- c("gene_a", "gene_b", "gene_c")
    task <- if (type == "binary") sb_binary_task(d, "patient_id", "outcome", genes) else
      sb_surv_task(d, "patient_id", "time", "event", genes, "month")
    p <- sb_pipeline(selector = if (type == "binary") "univ_binary_auc_top_k" else "univ_cox_top_k",
      learner = if (type == "binary") "logistic" else "coxph", selector_n_features = 2L)
    fits[[type]] <- sb_fit_final(task, p, capture_diagnostics = TRUE)
  }
  local_mocked_bindings(sb_fit_final = function(...) stop("Plotting must not refit"))
  for (type in names(fits)) {
    fit <- fits[[type]]
    expect_identical(fit$status, "success")
    metric <- if (type == "binary") "discrimination" else "wald_z"
    key <- paste0("selector_1_", metric)
    out <- sb_plots(fit)
    expect_s3_class(out$plots[[key]], "ggplot")
    expect_identical(attr(out$plots[[key]], "selector_scores"), fit$selection$state$scores)
    expect_setequal(out$plots[[key]]$data$value, fit$selection$state$scores[[metric]])
    expect_null(out$unavailable$selector_1)
    expect_s3_class(ggplot2::ggplotGrob(out$plots[[key]]), "gtable")
  }
})
})

test_that("RF Top-K removes lower-priority positive and negative correlations without refilling", {
  d <- data.frame(a = 1:8, b = -(1:8), c = (1:8)^2, d = c(7, 2, 6, 1, 4, 8, 3, 5),
    constant = rep(1, 8), outside = c(4, 7, 2, 1, 8, 6, 3, 5))
  result <- sb_rf_ranked_screen(d, c(5, 4, 4, 3, 2, 1), top_n = 5L, spearman_cutoff = .7)
  expect_identical(result$top_features, c("a", "b", "c", "d", "constant"))
  expect_identical(result$constant_features, "constant")
  expect_identical(result$retained_features, c("a", "d"))
  expect_equal(unname(result$correlation["a", c("b", "c")]), c(-1, 1))
  default <- sb_pipeline("rf_binary_importance", "logistic")
  expect_named(default$selector_params, c("trees", "seed", "importance", "threshold"))
  expect_error(sb_pipeline("rf_binary_importance", "logistic",
    selector_options = list(spearman_cutoff = .7)), class = "stormbreaker_bad_selector_options")
  expect_error(sb_pipeline("rf_binary_importance", "logistic",
    selector_options = list(top_n = 5L, threshold = 1)), class = "stormbreaker_bad_selector_options")
})

test_that("rank and correlation pruning use only analysis rows", {
  skip_if_not_installed("randomForest")
  d <- sb_example_binary(80L, 27L)
  d$proxy <- -d$gene_a
  features <- c("gene_a", "gene_b", "gene_c", "proxy")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  pipeline <- sb_pipeline("rf_binary_importance", "logistic",
    selector_options = list(trees = 40L, seed = 2L, top_n = 4L, spearman_cutoff = .7))
  expect_identical(sb_validate_pipeline_contract(pipeline, "binary"), pipeline)
  folds <- sb_make_stratified_folds(task, 3L, 11L)
  assignment <- stats::setNames(folds$assignment, d$patient_id)
  a <- sb_run_cv(task, pipeline, folds = assignment)
  changed <- d
  idx <- folds$splits[[1L]]$assessment_rows
  changed$gene_a[idx] <- changed$gene_a[idx] * -50
  changed$proxy[idx] <- changed$gene_b[idx]
  changed$outcome[idx] <- 1L - changed$outcome[idx]
  b <- sb_run_cv(sb_binary_task(changed, "patient_id", "outcome", features),
    pipeline, folds = assignment)
  expect_identical(a$selections[[1L]]$status, "success")
  expect_identical(a$selections[[1L]], b$selections[[1L]])
  fit <- sb_fit_final(task, pipeline)
  expect_identical(fit$status, "success")
  expect_length(intersect(fit$required_features, c("gene_a", "proxy")), 1L)
  expect_equal(fit$selection$state$evidence$rf_screen$correlation["gene_a", "proxy"], -1)
})
