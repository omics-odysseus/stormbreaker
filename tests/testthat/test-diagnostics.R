# ---- diagnostics-boosting ----
local({
# ---- algorithm-classical-boosting ----

test_that("classical boosting views retain original ensembles and exact native counts", {
  skip_if_not_installed("fastAdaboost")
  skip_if_not_installed("caTools")
  d <- make_sb_binary_fixture(n = 100L, seed = 108L)
  d$unused <- 0
  features <- c("gene_c", "unused", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  for (id in sb_classical_binary_boosting_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    original <- serialize(fit, NULL)
    engine <- unserialize(fit$learner_state$model_raw)
    usage <- sb_boosting_diagnostics(fit, "split_usage")
    expect_identical(usage$native_table, fit$learner_state$feature_importance)
    selected <- if (id == "adaboost_binary") unlist(lapply(engine$trees,
      function(tree) tree$frame$var[tree$frame$var != "<leaf>"])) else features[engine$Stump[, 1]]
    expect_equal(usage$table$count, as.numeric(tabulate(match(selected, features), nbins = length(features))))
    expect_equal(usage$table$count[usage$table$feature == "unused"], 0)
    if (id == "adaboost_binary") {
      weights <- sb_boosting_diagnostics(fit, "tree_weights")
      expect_identical(weights$native_table, engine$weights)
      expect_equal(weights$table$weight, as.numeric(engine$weights))
      sizes <- sb_boosting_diagnostics(fit, "tree_size")
      expect_equal(sizes$table$terminal_nodes, unname(vapply(engine$trees,
        function(tree) as.integer((nrow(tree$frame) + 1L) / 2L), integer(1))))
      expect_equal(nrow(weights$table), length(engine$trees))
      plots <- list(usage, weights, sizes)
    } else {
      stumps <- sb_boosting_diagnostics(fit, "stump_sequence")
      expect_identical(stumps$native_table, engine$Stump)
      expect_equal(stumps$table$feature, features[engine$Stump[, 1]])
      expect_equal(stumps$table$threshold, as.numeric(engine$Stump[, 2]))
      # Reconstruct class-1 probabilities independently from the returned rules.
      values <- d[features]
      margin <- numeric(nrow(values))
      for (i in seq_len(nrow(stumps$table))) {
        rule <- stumps$table[i, ]
        margin <- margin + ifelse(values[[rule$feature]] <= rule$threshold,
          rule$contribution_at_or_below, rule$contribution_above)
      }
      expected <- stats::predict(engine, as.matrix(values), type = "raw")[, "1"]
      expect_equal(stats::plogis(margin), unname(expected))
      plots <- list(usage, stumps)
    }
    if (requireNamespace("ggplot2", quietly = TRUE)) for (q in plots) expect_s3_class(
      ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    err <- tryCatch(sb_boosting_diagnostics(fit, "training_loss"), error = identity)
    expect_s3_class(err, "stormbreaker_unavailable_diagnostic")
    expect_identical(serialize(fit, NULL), original)
  }
})

test_that("AdaBoost early convergence does not produce invented extra trees", {
  skip_if_not_installed("fastAdaboost")
  y <- rep(c(0L, 1L), each = 10L)
  d <- data.frame(id = seq_along(y), outcome = y, separating = ifelse(y == 1L, 1, -1), noise = seq_along(y))
  fit <- sb_fit_final(sb_binary_task(d, "id", "outcome", c("separating", "noise")),
    sb_pipeline(learner = "adaboost_binary", adaboost_iterations = 5L))
  expect_identical(fit$status, "success")
  q <- sb_boosting_diagnostics(fit, "tree_weights")
  expect_equal(q$fitted_iterations, 1L)
  expect_equal(q$requested_iterations, 5L)
  expect_equal(nrow(q$table), 1L)
  expect_equal(q$table$weight, as.numeric(unserialize(fit$learner_state$model_raw)$weights))
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(
    ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
})

# ---- algorithm-hinge-boosting ----

test_that("captured hinge paths reproduce native updates without changing portable state", {
  skip_if_not_installed("bst")
  d <- make_sb_binary_fixture(n = 100L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  pipeline <- sb_pipeline(learner = "bst_hinge_binary")
  compact <- sb_fit_final(task, pipeline)
  fit <- sb_fit_final(task, pipeline, capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  expect_identical(compact$status, "success")
  expect_identical(fit$learner_state, compact$learner_state)
  expect_null(fit$learner_state$engine$diagnostic_path)
  err <- tryCatch(sb_boosting_diagnostics(compact, "training_loss"), error = identity)
  expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
  original <- serialize(fit, NULL)
  params <- pipeline$learner_params
  native <- bst::bst(as.matrix(d[features]), ifelse(d$outcome == 1, 1, -1),
    cost = params$cost, family = params$family, learner = params$base_learner,
    ctrl = bst::bst_control(mstop = params$iterations, nu = params$learning_rate,
      center = params$center, intercept = params$fit_intercept, trace = FALSE))
  loss <- sb_boosting_diagnostics(fit, "training_loss")
  usage <- sb_boosting_diagnostics(fit, "component_usage")
  path <- sb_boosting_diagnostics(fit, "coefficient_path")
  expect_equal(loss$table$loss, native$risk)
  expect_equal(loss$table$iteration, seq_len(params$iterations))
  expect_equal(usage$table$updates, tabulate(native$ensemble, nbins = length(features)))
  expect_equal(path$selections$native_feature_index, native$ensemble)
  expect_equal(path$selections$native_update, native$coef)
  expect_equal(path$table$feature[path$table$iteration == 0], c("(Intercept)", features))
  reconstructed_loss <- vapply(seq_len(params$iterations), function(i) {
    beta <- path$table$coefficient[path$table$iteration == i]
    margin_i <- beta[1] + as.numeric(as.matrix(d[features]) %*% beta[-1])
    mean(bst::loss(ifelse(d$outcome == 1, 1, -1), margin_i, cost = params$cost, family = "hinge"))
  }, numeric(1))
  expect_equal(reconstructed_loss, native$risk, tolerance = 1e-12)
  endpoint <- path$table[path$table$iteration == params$iterations, ]$coefficient
  margin <- endpoint[1] + as.numeric(as.matrix(d[features]) %*% endpoint[-1])
  expect_equal(margin, as.numeric(native$yhat), tolerance = 1e-12)
  expect_equal(stats::plogis(2 * margin), sb_predict_omics_specialist_binary_state(
    fit$learner_state, d[features])$probability, tolerance = 1e-12)
  if (requireNamespace("ggplot2", quietly = TRUE)) for (q in list(loss, usage, path)) expect_s3_class(
    ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
  expect_identical(serialize(fit, NULL), original)
})
})

# ---- diagnostics-core ----
local({
# ---- algorithm-diagnostics ----

test_that("optional diagnostic capture preserves original model predictions", {
  skip_if_not_installed("glmnet")
  d <- make_sb_binary_fixture(n = 90L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  pipeline <- sb_pipeline(learner = "glmnet_logistic", inner_folds = 3L, inner_seed = 5L)
  a <- sb_fit_final(task, pipeline)
  b <- sb_fit_final(task, pipeline, capture_diagnostics = TRUE)
  expect_identical(b$status, "success")
  expect_null(a$diagnostic_artifacts)
  expect_identical(a$learner_state, b$learner_state)
  expect_identical(sb_predict_final(a, d, "patient_id"), sb_predict_final(b, d, "patient_id"))
  q <- sb_algorithm_diagnostics(b)
  expect_equal(q$table$coefficient, as.vector(as.matrix(b$diagnostic_artifacts$learner$beta)))
  expect_equal(length(unique(q$table$feature)), 2)
  skip_if_not_installed("ggplot2")
  expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
})

test_that("ABESS diagnostics preserve the native GIC path and chosen support", {
  skip_if_not_installed("abess", minimum_version = "0.4.11")
  for (kind in c("binary", "survival")) {
    if (kind == "binary") {
      d <- make_sb_binary_fixture(n = 84L, seed = 8611L)
      task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
      learner <- "abess_binary_gic"
    } else {
      d <- make_sb_fixture(n = 90L, seed = 8612L)
      task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
      learner <- "abess_cox_gic"
    }
    fit <- sb_fit_final(task, sb_pipeline(learner = learner), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- fit$learner_state
    q <- sb_algorithm_diagnostics(fit, "support_criterion")
    native <- fit$diagnostic_artifacts$learner
    coefficient_path <- sb_algorithm_diagnostics(fit, "coefficient_path")
    expect_equal(coefficient_path$table$coefficient, as.vector(as.matrix(native$engine$beta)))
    expect_s3_class(ggplot2::ggplotGrob(sb_plot_algorithm_diagnostic(coefficient_path)), "gtable")
    expect_identical(q$table[names(native$criterion_path)], native$criterion_path)
    expect_identical(q$table$support_size[q$table$selected], native$selected_support_size)
    expect_identical(fit$learner_state, before)
    expect_equal(sum(q$table$selected), 1)
    damaged <- fit
    damaged$diagnostic_artifacts$learner$criterion_path$gic[1] <- NA_real_
    err <- tryCatch(sb_algorithm_diagnostics(damaged, "support_criterion"), error = identity)
    expect_s3_class(err, "stormbreaker_invalid_diagnostic")
    missing <- fit
    missing$diagnostic_artifacts <- NULL
    err <- tryCatch(sb_algorithm_diagnostics(missing, "support_criterion"), error = identity)
    expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    }
  }
})

# ---- native-diagnostics ----

test_that("automatic extraction executes native routes and preserves unavailable reasons", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_binary(n = 100L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  for (id in c("gaussian_naive_bayes", "logistic", "decision_tree_binary", "knn_binary")) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    auto <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(serialize(fit, NULL, version = 2L), before)
    if (id == "gaussian_naive_bayes") {
      expect_identical(auto$status$status, c("available", "available"))
      expect_identical(auto$diagnostics[[1]], sb_naive_bayes_diagnostics(fit))
      expect_length(sb_native_diagnostics(fit)$diagnostics, 1L)
    } else if (id == "logistic") {
      expect_identical(auto$status$status, "available")
      expect_identical(auto$diagnostics[[1]]$table, sb_model_coefficients(fit)$table)
    } else {
      expected <- if (id == "knn_binary") "unavailable" else "missing_capture"
      expect_identical(auto$status$status, expected)
      expect_length(auto$diagnostics, 0L)
      result <- sb_run_cv(task, sb_pipeline(learner = id), folds = 3L, seed = 41L)
      views <- sb_clinician_views(result, native = auto)
      expect_true(expected %in% views$status$status)
      expect_identical(views$diagnostics[["Native extraction"]], auto$status)
    }
  }
  fit <- sb_fit_final(task, sb_pipeline(learner = "decision_tree_binary"), capture_diagnostics = TRUE)
  expect_identical(sb_native_diagnostics(fit)$diagnostics[[1]], sb_tree_structure(fit))
})

test_that("all catalog routes name real compatible extractors without claiming availability", {
  plans <- list()
  for (type in c("binary", "survival")) {
    catalog <- sb_model_catalog(type)
    expect_gt(nrow(catalog), 0L)
    plans[[type]] <- do.call(rbind, lapply(unique(catalog$learner), sb_native_diagnostic_plan))
  }
  plan <- do.call(rbind, plans)
  for (name in unique(plan$extractor)) {
    extractor <- get(name, envir = asNamespace("stormbreaker"))
    expect_true(is.function(extractor), info = name)
    expect_true("fit" %in% names(formals(extractor)), info = name)
    if (any(nzchar(plan$type[plan$extractor == name]))) {
      expect_true("type" %in% names(formals(extractor)), info = name)
    }
  }
})

test_that("automatic glmnet detail distinguishes compact tuning from missing native paths", {
  skip_if_not_installed("glmnet")
  d <- sb_example_binary(n = 100L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "glmnet_logistic", inner_folds = 3L, inner_seed = 5L))
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  auto <- sb_native_diagnostics(fit, complete = TRUE)
  expect_identical(auto$status$type, c("cv_curve", "nonzero_path", "coefficient_path"))
  expect_identical(auto$status$status, c("available", "available", "missing_capture"))
  expect_identical(auto$diagnostics[[1]], sb_algorithm_diagnostics(fit, "cv_curve"))
  expect_identical(serialize(fit, NULL, version = 2L), before)
})
})

# ---- diagnostics-forests ----
local({
# ---- algorithm-aorsf ----

test_that("aorsf views preserve native leaf counts and final-only OOB scores", {
  skip_if_not_installed("aorsf")
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in c("aorsf_binary", "aorsf_survival")) {
    binary <- id == "aorsf_binary"
    d <- if (binary) sb_example_binary(100L, seed = 108L) else sb_example_survival(100L, seed = 8612L)
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    prediction <- sb_predict_final(fit, d, "patient_id")
    e <- unserialize(fit$learner_state$engine_raw)
    sizes <- sb_aorsf_diagnostics(fit)
    oob <- sb_aorsf_diagnostics(fit, "oob_score")
    expect_equal(mean(sizes$table$terminal_nodes), e$get_mean_leaves_per_tree(), tolerance = 1e-12)
    expect_equal(sizes$table$terminal_nodes,
      vapply(e$forest$leaf_summary, function(z) length(which(z != 0)), integer(1)))
    expect_equal(nrow(sizes$table), 500L)
    expect_equal(oob$table$score, as.numeric(e$eval_oobag$stat_values), tolerance = 0)
    expect_identical(oob$measure, if (binary) "AUC-ROC" else "Harrell's C-index")
    expect_equal(oob$table$trees, 500)
    expect_equal(e$oobag_eval_every, 500)
    expect_equal(nrow(oob$table), 1L)
    expect_null(fit$diagnostic_artifacts)
    expect_null(e$importance)
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 2L))
    for (v in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(v)
      b <- ggplot2::ggplot_build(p)
      expect_s3_class(b, "ggplot_built")
      expect_null(p$labels$caption)
      expect_null(p$labels$subtitle)
      if (v$type == "aorsf_tree_size") expect_equal(sum(p$data$trees), 500L) else {
        expect_length(b$data, 1L)
        expect_equal(b$data[[1]]$y, oob$table$score, tolerance = 0)
        expect_equal(nrow(b$data[[1]]), 1L)
      }
    }
    expect_identical(sb_predict_final(fit, d, "patient_id"), prediction)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})

# ---- algorithm-bagging ----

test_that("bagging views preserve every original tree and native importance entry", {
  skip_if_not_installed("ipred")
  d <- sb_example_binary(100L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = "bagging_binary"))
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  prediction <- sb_predict_final(fit, d, "patient_id")
  e <- unserialize(fit$learner_state$engine_raw)
  sizes <- sb_bagging_diagnostics(fit, "tree_size")
  usage <- sb_bagging_diagnostics(fit)
  imp <- sb_bagging_diagnostics(fit, "importance")
  for (i in seq_along(e$mtrees)) {
    tree <- e$mtrees[[i]]$btree
    ids <- as.integer(rownames(tree$frame))
    expect_equal(sizes$table$terminal_nodes[i], sum(!((2L * ids) %in% ids)))
    expect_equal(sizes$table$nodes[i], length(ids))
    rows <- imp$table[imp$table$tree == i, ]
    expect_identical(rows$feature, features)
    native <- tree$variable.importance
    for (j in seq_along(features)) expect_equal(rows$importance[j],
      if (features[j] %in% names(native)) unname(native[features[j]]) else 0, tolerance = 0)
  }
  expect_equal(nrow(imp$table), 300L)
  expect_equal(sum(usage$table$splits), sum(sizes$table$nodes - sizes$table$terminal_nodes))
  vars <- unlist(lapply(e$mtrees, function(x) x$btree$frame$var))
  expect_equal(usage$table$splits, as.numeric(table(factor(vars, levels = features))))
  expect_null(fit$diagnostic_artifacts)
  bundle <- sb_native_diagnostics(fit, complete = TRUE)
  expect_identical(bundle$status$status, rep("available", 3L))
  for (v in bundle$diagnostics) {
    p <- sb_plot_algorithm_diagnostic(v)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    expect_null(p$labels$caption)
    expect_null(p$labels$subtitle)
    if (v$type == "bagging_tree_size") expect_equal(sum(p$data$trees), 100L)
    if (v$type == "bagging_importance") expect_equal(p$data$value,
      as.numeric(tapply(imp$table$importance, imp$table$feature, mean)), tolerance = 0)
  }
  expect_identical(sb_predict_final(fit, d, "patient_id"), prediction)
  expect_identical(serialize(fit, NULL, version = 2L), before)
})

# ---- algorithm-cforest ----

test_that("cforest diagnostics map native split columns and retain every tree", {
  skip_if_not_installed("partykit")
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in c("cforest_binary", "cforest_survival")) {
    binary <- id == "cforest_binary"
    d <- if (binary) sb_example_binary(100L, seed = 108L) else sb_example_survival(100L, seed = 8612L)
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    prediction <- sb_predict_final(fit, d, "patient_id")
    e <- unserialize(fit$learner_state$engine_raw)
    sizes <- sb_cforest_diagnostics(fit, "tree_size")
    usage <- sb_cforest_diagnostics(fit)
    expect_equal(sizes$table$terminal_nodes, vapply(e$nodes, function(z) length(partykit::nodeids(z, terminal = TRUE)), integer(1)))
    expect_equal(sizes$table$nodes, vapply(e$nodes, function(z) length(partykit::nodeids(z)), integer(1)))
    native_features <- unlist(lapply(e$nodes, function(z) {
      ids <- setdiff(partykit::nodeids(z), partykit::nodeids(z, terminal = TRUE))
      unlist(partykit::nodeapply(z, ids = ids, FUN = function(n) {
        partykit::character_split(partykit::split_node(n), data = e$data)$name
      }), use.names = FALSE)
    }), use.names = FALSE)
    expect_equal(usage$table$splits, as.numeric(table(factor(native_features, levels = features))))
    expect_equal(usage$table$feature, features)
    expect_equal(sum(usage$table$splits), sum(sizes$table$nodes - sizes$table$terminal_nodes))
    expect_equal(nrow(sizes$table), 500L)
    expect_null(fit$diagnostic_artifacts)
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 2L))
    for (v in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(v)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_null(p$labels$caption)
      expect_null(p$labels$subtitle)
      if (v$type == "cforest_tree_size") expect_equal(sum(p$data$trees), 500L)
    }
    expect_identical(sb_predict_final(fit, d, "patient_id"), prediction)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})

# ---- algorithm-ferns ----

test_that("fern split mapping reproduces native scores without changing fits", {
  skip_if_not_installed("rFerns")
  for (features in list("gene_c", c("gene_c", "gene_a", "gene_b"))) {
    d <- sb_example_binary(100L, seed = 108L)
    fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features),
      sb_pipeline(learner = "random_ferns_binary"))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    e <- unserialize(fit$learner_state$engine_raw)
    usage <- sb_ferns_diagnostics(fit)
    thresholds <- sb_ferns_diagnostics(fit, "thresholds")
    expect_equal(usage$table$splits, as.integer(table(factor(thresholds$table$feature, levels = features))))
    expect_identical(thresholds$table$threshold, e$model$threReal)
    query <- d[c(1L, 8L, 13L), features, drop = FALSE]
    query[[1L]] <- query[[1L]] + .137
    depth <- unname(e$parameters["depth"])
    nferns <- unname(e$parameters["ferns"])
    scores <- array(e$model$scores, c(length(e$classLabels), 2^depth, nferns))
    reconstructed <- matrix(0, nrow(query), length(e$classLabels))
    for (f in seq_len(nferns)) {
      z <- thresholds$table[thresholds$table$fern == f, ]
      leaves <- vapply(seq_len(nrow(query)), function(i) {
        sum((as.numeric(query[i, z$feature, drop = FALSE]) < z$threshold) * 2^(z$split - 1L)) + 1L
      }, numeric(1))
      reconstructed <- reconstructed + t(scores[, leaves, f, drop = FALSE][, , 1L])
    }
    native <- get("predict.rFerns", asNamespace("rFerns"))(e, query, scores = TRUE)
    expect_equal(reconstructed / nferns, unname(as.matrix(native)), tolerance = 1e-12)
    for (z in list(usage, thresholds)) {
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    }
    expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status, rep("available", 2L))
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})

# ---- algorithm-rotation-forest ----

test_that("rotation forest diagnostics preserve original bases and native predictions", {
  skip_if_not_installed("rotationForest")
  d <- sb_example_binary(160L,seed=108L)
  features <- c("gene_c","gene_a","gene_b")
  fit <- sb_fit_final(sb_binary_task(d,"patient_id","outcome",features),sb_pipeline(learner="rotation_forest_binary"))
  expect_identical(fit$status,"success")
  before <- serialize(fit,NULL,version=2L)
  e <- unserialize(fit$learner_state$engine_raw)
  loadings <- sb_rotation_forest_diagnostics(fit)
  sizes <- sb_rotation_forest_diagnostics(fit,"tree_size")
  query <- d[c(1L,7L), ]; query$gene_c <- query$gene_c+.17
  predictions <- sb_rotation_forest_diagnostics(fit,"tree_predictions",query)
  for(i in seq_along(e$models)) {
    tab <- loadings$table[loadings$table$tree==i,]
    m <- matrix(tab$loading,nrow=length(features))
    expect_equal(m,unname(as.matrix(e$loadings[[i]])))
    expect_identical(tab$feature,rep(features,length(features)))
    ids <- as.integer(rownames(e$models[[i]]$frame))
    expect_equal(sizes$table$terminal_nodes[i],sum(!((2L*ids) %in% ids)))
    rotated <- as.data.frame(as.matrix(query[features]) %*% as.matrix(e$loadings[[i]]))
    native <- predict(e$models[[i]],rotated,type="prob")[,"1"]
    expect_equal(predictions$table$probability[predictions$table$tree==i],unname(native))
  }
  expect_equal(as.numeric(tapply(predictions$table$probability,predictions$table$query,mean)),
    sb_predict_final(fit,query,"patient_id")$predictions$probability,tolerance=1e-12)
  for(z in list(loadings,sizes,predictions)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)),"ggplot_built")
  expect_identical(sb_native_diagnostics(fit,complete=TRUE)$status$status,c("available","available","unavailable"))
  expect_identical(sb_native_diagnostics(fit,complete=TRUE,newdata=query)$status$status,rep("available",3))
  expect_identical(serialize(fit,NULL,version=2L),before)
})
})

# ---- diagnostics-kernels ----
local({
# ---- algorithm-gaussian-process ----

test_that("Gaussian-process kernel terms preserve native probabilities", {
  skip_if_not_installed("kernlab")
  d <- sb_example_binary(100L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features),
    sb_pipeline(learner = "gaussian_process_binary"))
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  e <- unserialize(fit$learner_state$engine_raw)
  weights <- sb_gaussian_process_diagnostics(fit)
  expect_equal(weights$table$weight, -as.numeric(e@alpha[[1L]]))
  expect_identical(weights$table$reference_row, e@alphaindex[[1L]])
  for (rows in list(8L, c(1L, 8L, 19L))) {
    q <- d[rows, ]; q$gene_c <- q$gene_c * 1.3 + .71
    z <- sb_gaussian_process_diagnostics(fit, "query_contributions", q)
    native <- kernlab::predict(e, as.matrix(q[features]), type = "probabilities")
    expect_equal(z$summary$probability, as.numeric(native[, "1"]), tolerance = 1e-12)
    expect_equal(plogis(z$summary$latent_score), z$summary$probability, tolerance = 1e-12)
    for (i in seq_along(rows)) {
      scaled <- (as.numeric(q[i, features]) - e@scaling$x.scale[["scaled:center"]]) /
        e@scaling$x.scale[["scaled:scale"]]
      terms <- vapply(seq_along(e@alphaindex[[1L]]), function(j) {
        as.numeric(e@kernelf(scaled, e@xmatrix[e@alphaindex[[1L]][j], ])) * -e@alpha[[1L]][j]
      }, numeric(1))
      expect_equal(z$table$contribution[z$table$query == i], terms)
    }
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    expect_identical(sb_native_diagnostics(fit, complete = TRUE, newdata = q)$status$status, rep("available", 2L))
  }
  expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(weights)), "ggplot_built")
  expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status, c("available", "unavailable"))
  expect_identical(serialize(fit, NULL, version = 2L), before)
})

# ---- algorithm-lssvm ----

test_that("LS-SVM preserves native kernel scores rather than inventing support vectors", {
  skip_if_not_installed("kernlab")
  d <- sb_example_binary(n = 100L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  fit <- sb_fit_final(task, sb_pipeline(learner = "lssvm_radial_binary"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  q <- sb_svm_diagnostics(fit, "decision_distribution")
  e <- fit$diagnostic_artifacts$learner$engine
  z <- scale(data.matrix(d[features]), center = e$center, scale = e$scale)
  K <- kernlab::kernelMatrix(kernlab::rbfdot(sigma = e$sigma), z, e$scaled_training)
  native <- sweep(as.matrix(K) %*% e$alpha, 2L, e$intercept, "+")
  expect_equal(q$native_table, unname(native), tolerance = 1e-10)
  expect_equal(q$table$decision, as.numeric(native[, 2L] - native[, 1L]), tolerance = 1e-10)
  expect_equal(plogis(q$table$decision), sb_predict_lssvm_radial_binary(e, data.matrix(d[features])))
  expect_true(q$sensitive_training_data)
  expect_null(rownames(q$native_table))
  expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
  for (type in c("support_counts", "coefficients")) expect_error(sb_svm_diagnostics(fit, type),
    class = "stormbreaker_unavailable_diagnostic")
  expect_true(identical(serialize(fit, NULL, version = 2L), before))
  fit$diagnostic_artifacts <- NULL
  expect_error(sb_svm_diagnostics(fit, "decision_distribution"), class = "stormbreaker_missing_diagnostic_artifacts")
})

# ---- algorithm-survival-svm ----

test_that("all twelve survival-SVM variants retain native scores and basis weights", {
  skip_if_not_installed("survivalsvm")
  task <- make_survival_svm_diagnostic_task()
  for (id in sb_survivalsvm_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success", info = id)
    if (fit$status != "success") next
    before <- serialize(fit, NULL, version = 2L)
    model <- fit$diagnostic_artifacts$learner
    a <- model$audit
    x <- fit$diagnostic_artifacts$predictors
    z <- as.data.frame(scale(as.matrix(x), center = a$center, scale = a$scale))
    score <- sb_survival_svm_diagnostics(fit)
    native <- as.numeric(stats::predict(model$engine, newdata = z)$predicted)
    expect_equal(score$table$value, native, tolerance = 1e-10)
    q <- sb_survival_svm_diagnostics(fit, "basis_weights")
    m <- model$engine$model.fit
    w <- if (model$params$type == "regression") as.numeric(m$Beta) else {
      v <- as.numeric(t(m$DifMat$Mat) %*% m$Alpha)
      if (model$params$type == "hybrid") v <- v + m$Beta - m$Delta * m$Betarstar
      as.numeric(v)
    }
    expect_equal(q$table$value, w)
    expect_true(q$sensitive_training_data)
    constant <- q
    constant$table$value[] <- 0
    constant_plot <- ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(constant))
    expect_equal(constant_plot$data[[2L]]$y, 1)
    expect_true(score$sensitive_training_data)
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(score)), "ggplot_built")
    if (id %in% sb_survivalsvm_linear_learner_ids()) {
      beta <- sb_survival_svm_diagnostics(fit, "coefficients")
      expect_equal(as.numeric(as.matrix(x) %*% beta$table$value + beta$intercept), native, tolerance = 1e-10)
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(beta)), "ggplot_built")
    } else expect_error(sb_survival_svm_diagnostics(fit, "coefficients"), class = "stormbreaker_unavailable_diagnostic")
    expect_true(identical(serialize(fit, NULL, version = 2L), before))
    fit$diagnostic_artifacts <- NULL
    expect_error(sb_survival_svm_diagnostics(fit), class = "stormbreaker_missing_diagnostic_artifacts")
  }
})

# ---- algorithm-svm ----

test_that("five LIBSVM learners expose native counts and oriented decisions", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(n = 100L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in c("svm_linear_binary", "svm_radial_binary",
      "svm_radial_weighted_binary", "svm_polynomial_binary", "svm_sigmoid_binary")) {
    task <- sb_binary_task(d, "patient_id", "outcome", features)
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    original <- serialize(fit, NULL, version = 2L)
    e <- fit$diagnostic_artifacts$learner$engine
    counts <- sb_svm_diagnostics(fit)
    native_counts <- table(factor(d$outcome[e$index], levels = c(0, 1)))
    expect_equal(counts$table$support_vectors[match(c("0", "1"), counts$table$class)], as.integer(native_counts))
    q <- sb_svm_diagnostics(fit, "decision_distribution")
    native <- stats::predict(e, data.matrix(d[features]), decision.values = TRUE)
    expect_identical(q$native_table, unname(attr(native, "decision.values")))
    expect_equal(as.integer(q$table$decision > 0), as.integer(as.character(native)))
    expect_identical(q$table$class, as.character(d$outcome))
    expect_true(q$sensitive_training_data)
    expect_null(rownames(q$native_table))
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(counts)), "ggplot_built")
    compact <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(compact$status, "success")
    expect_identical(compact$learner_state, fit$learner_state)
    expect_identical(sb_svm_diagnostics(compact), counts)
    expect_error(sb_svm_diagnostics(compact, "decision_distribution"), class = "stormbreaker_missing_diagnostic_artifacts")
    bad <- fit
    bad$diagnostic_artifacts$predictors <- bad$diagnostic_artifacts$predictors[rev(features)]
    expect_error(sb_svm_diagnostics(bad, "decision_distribution"), class = "stormbreaker_invalid_diagnostic")
    if (id != "svm_linear_binary") expect_error(sb_svm_diagnostics(fit, "coefficients"),
      class = "stormbreaker_unavailable_diagnostic")
    expect_true(identical(serialize(fit, NULL, version = 2L), original))
  }
})

test_that("linear SVM weights reproduce native decisions for either native label order", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(n = 100L, seed = 108L)
  d$gene_a <- 30 + 15 * d$gene_a
  features <- c("gene_b", "gene_a", "gene_c")
  orientations <- numeric()
  for (descending in c(FALSE, TRUE)) {
    data <- d[order(d$outcome, decreasing = descending), ]
    task <- sb_binary_task(data, "patient_id", "outcome", features)
    fit <- sb_fit_final(task, sb_pipeline(learner = "svm_linear_binary"), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    q <- sb_svm_diagnostics(fit, "coefficients")
    scores <- sb_svm_diagnostics(fit, "decision_distribution")
    expect_equal(as.numeric(data.matrix(data[features]) %*% q$table$coefficient + q$intercept),
      scores$table$decision, tolerance = 1e-10)
    expect_identical(q$table$feature, features)
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    orientations <- c(orientations, q$orientation)
  }
  expect_setequal(orientations, c(-1, 1))
})
})

# ---- diagnostics-linear ----
local({
# ---- algorithm-discriminant ----

test_that("discriminant diagnostics preserve native prediction equations and compact state", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("klaR")
  d <- sb_example_binary(n = 160L, seed = 108L)
  d$gene_a <- 10 + d$gene_a
  d$gene_b <- 3 * d$gene_b - 5
  for (features in list(c("gene_c", "gene_a", "gene_b"), "gene_a")) {
    task <- sb_binary_task(d, "patient_id", "outcome", features)
    for (id in c("lda_binary", "qda_binary", "rda_binary")) {
      fit <- sb_fit_final(task, sb_pipeline(learner = id))
      expect_identical(fit$status, "success", info = id)
      expect_null(fit$diagnostic_artifacts)
      before <- serialize(fit, NULL, version = 2L)
      diagnostic <- sb_discriminant_diagnostics(fit)
      p <- diagnostic$parameters
      x <- as.matrix(d[features])
      score <- p$intercept + drop(x %*% p$coefficient) + rowSums((x %*% p$quadratic) * x)
      native <- sb_predict_final(fit, d, "patient_id")$predictions$probability
      expect_equal(plogis(score), native, tolerance = 1e-12, info = id)
      e <- unserialize(fit$learner_state$engine_raw)
      expect_equal(p$means, if (id == "rda_binary") t(e$means) else e$means)
      expect_equal(diagnostic$table$value, as.vector(p$means))
      expect_equal(p$priors, e$prior)
      expect_identical(names(p$coefficient), features)
      expect_identical(dimnames(p$quadratic), list(features, features))
      if (id == "lda_binary") {
        expect_true(all(p$quadratic == 0))
        expect_error(sb_discriminant_diagnostics(fit, "class_correlation"),
          class = "stormbreaker_unavailable_diagnostic")
        coefficients <- sb_model_coefficients(fit, d)
        expect_equal(unname(coefficients$table$coefficient), unname(p$coefficient))
        nomogram <- sb_nomogram(fit, d)
        expect_equal(sb_nomogram_values(nomogram, d)$score, score, tolerance = 1e-12)
        shap <- sb_additive_shap(fit, d, "patient_id", d)
        expect_equal(unname(rowSums(shap$contributions)) + shap$reference_value, score, tolerance = 1e-12)
      } else {
        correlation <- sb_discriminant_diagnostics(fit, "class_correlation")
        expect_equal(nrow(correlation$table), 2 * length(features)^2)
        for (class in c("0", "1")) expect_equal(correlation$table$value[correlation$table$class == class],
          as.vector(stats::cov2cor(p$covariance[[class]])))
        expect_error(sb_nomogram(fit, d), class = "stormbreaker_unsupported_coefficients")
        expect_error(sb_additive_shap(fit, d, "patient_id", d), class = "stormbreaker_unsupported_additive_shap")
      }
      expect_error(sb_discriminant_diagnostics(fit, "unknown"), class = "stormbreaker_unsupported_diagnostic")
      bundle <- sb_native_diagnostics(fit, complete = TRUE)
      expect_identical(bundle$status$status, rep("available", 2L))
      if (requireNamespace("ggplot2", quietly = TRUE)) for (view in bundle$diagnostics) {
        plot <- sb_plot_algorithm_diagnostic(view)
        expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
        expect_null(plot$labels$subtitle)
        expect_null(plot$labels$caption)
      }
      expect_identical(serialize(fit, NULL, version = 2L), before)
    }
  }
})

# ---- algorithm-fda ----

test_that("FDA/PDA/MDA native geometry reconstructs posterior probabilities", {
  skip_if_not_installed("mda")
  d <- sb_example_binary(160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in c("fda_binary", "pda_binary", "mda_binary")) {
    fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    e <- unserialize(fit$learner_state$engine_raw)
    centers <- sb_fda_diagnostics(fit)
    means <- matrix(centers$table$coordinate, nrow = nrow(e$means))
    expect_equal(means, unname(as.matrix(e$means)), ignore_attr = TRUE)
    mapping <- centers$table[centers$table$axis == 1L, ]
    weights <- if (id == "mda_binary") sb_fda_diagnostics(fit, "mixture_weights") else NULL
    if (!is.null(weights)) {
      expect_equal(weights$table$weight, as.numeric(unlist(e$sub.prior)))
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(weights)), "ggplot_built")
    }
    for (rows in list(8L, c(1L, 8L, 13L))) {
      q <- d[rows, ]; q$gene_c <- q$gene_c + .71
      z <- sb_fda_diagnostics(fit, "query_variates", q)
      coordinates <- matrix(z$table$coordinate, nrow = nrow(q))
      batch <- rbind(q[features], d[22L, features, drop = FALSE])
      native_coordinates <- predict(e, batch, type = "variates", dimension = e$dimension)[seq_len(nrow(q)), , drop = FALSE]
      expect_equal(coordinates, unname(as.matrix(native_coordinates)), ignore_attr = TRUE)
      mass <- matrix(0, nrow(q), 2L)
      for (i in seq_len(nrow(means))) {
        cls <- match(mapping$class[i], c("0", "1"))
        weight <- if (is.null(weights)) 1 else weights$table$weight[i]
        distance <- rowSums(sweep(coordinates, 2L, means[i, ], "-")^2)
        mass[, cls] <- mass[, cls] + exp(-distance / 2) * weight * as.numeric(e$prior[cls])
      }
      posterior <- mass[, 2L] / rowSums(mass)
      expect_equal(posterior, sb_predict_final(fit, q, "patient_id")$predictions$probability, tolerance = 1e-12)
      native_posterior <- predict(e, batch, type = "posterior")[seq_len(nrow(q)), "1"]
      expect_equal(posterior, unname(native_posterior), tolerance = 1e-12)
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    }
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(centers)), "ggplot_built")
    expect_identical(tail(sb_native_diagnostics(fit, complete = TRUE)$status$status, 1L), "unavailable")
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})

# ---- algorithm-linear-margin ----

test_that("LIBLINEAR SVC and DWD weights preserve the original decision function", {
  skip_if_not_installed("LiblineaR")
  skip_if_not_installed("DWDLargeR")
  d <- sb_example_binary(n = 100L, seed = 108L)
  d$gene_a <- 20 + 5 * d$gene_a
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in sb_linear_margin_diagnostic_ids()) {
    task <- sb_binary_task(d, "patient_id", "outcome", features)
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    q <- sb_svm_diagnostics(fit, "coefficients")
    scores <- sb_svm_diagnostics(fit, "decision_distribution")
    expect_identical(q$table$feature, features)
    expect_equal(as.numeric(data.matrix(d[features]) %*% q$table$coefficient + q$intercept), scores$table$decision)
    e <- fit$diagnostic_artifacts$learner$engine
    if (id != "dwd_binary") {
      x <- sweep(sweep(data.matrix(d[features]), 2L, e$center, "-"), 2L, e$scale, "/")
      native <- stats::predict(e$model, x, decisionValues = TRUE)$decisionValues
      expect_equal(scores$table$decision, as.numeric(native[, "1"] - native[, "0"]), tolerance = 1e-10)
      p <- sb_predict_liblinear_binary(e, data.matrix(d[features]))
    } else {
      expect_identical(q$table$coefficient, as.numeric(e$coefficient))
      p <- plogis(2 * (e$intercept + data.matrix(d[features]) %*% e$coefficient))
    }
    factor <- if (id == "dwd_binary") 2 else 1
    expect_equal(plogis(factor * scores$table$decision), as.numeric(p), tolerance = 1e-10)
    for (plot in list(q, scores)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(plot)), "ggplot_built")
    expect_true(scores$sensitive_training_data)
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_identical(sb_svm_diagnostics(compact, "coefficients"), q)
    expect_error(sb_svm_diagnostics(compact, "decision_distribution"), class = "stormbreaker_missing_diagnostic_artifacts")
    expect_error(sb_svm_diagnostics(fit, "support_counts"), class = "stormbreaker_unavailable_diagnostic")
    expect_true(identical(serialize(fit, NULL, version = 2L), before))
  }
})

# ---- algorithm-shrinkage-discriminant ----

test_that("SDA and SLDA diagnostics preserve native prediction geometry", {
  skip_if_not_installed("sda")
  skip_if_not_installed("ipred")
  d <- sb_example_binary(160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in c("sda_binary", "slda_binary")) {
    fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    e <- unserialize(fit$learner_state$engine_raw)
    first <- sb_shrinkage_discriminant_diagnostics(fit)
    second <- sb_shrinkage_discriminant_diagnostics(fit, if (id == "sda_binary") "regularization" else "class_means")
    if (id == "sda_binary") {
      expect_identical(first$table$feature, features)
      expect_equal(first$table$coefficient, as.numeric(e$beta[2L, ] - e$beta[1L, ]))
      expect_identical(second$table$value, as.numeric(e$regularization))
    } else {
      loadings <- matrix(first$table$loading, nrow = length(features))
      expect_equal(loadings, unname(e$scores))
      expect_identical(first$table$feature, rep(features, ncol(e$scores)))
      expect_equal(matrix(second$table$mean, nrow = 2L), unname(e$mylda$means))
    }
    for (rows in list(8L, c(1L, 8L, 13L))) {
      q <- d[rows, ]; q$gene_c <- q$gene_c + .71
      if (id == "sda_binary") {
        p <- plogis(as.numeric(as.matrix(q[features]) %*% first$table$coefficient) + first$intercept)
        reconstructed <- zapsmall(cbind(1 - p, p))[, 2L]
        native <- get("predict.sda", asNamespace("sda"))(e, as.matrix(q[features]), verbose = FALSE)$posterior[, "1"]
      } else {
        projected <- as.data.frame(as.matrix(q[features]) %*% loadings)
        names(projected) <- colnames(e$mylda$means)
        reconstructed <- predict(e$mylda, projected)$posterior[, "1"]
        native <- sb_predict_slda_binary(e, as.matrix(q[features]))
      }
      expect_equal(as.numeric(reconstructed), as.numeric(native), tolerance = 1e-12)
      expect_equal(as.numeric(reconstructed), sb_predict_final(fit, q, "patient_id")$predictions$probability, tolerance = 1e-12)
    }
    for (z in list(first, second)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status, rep("available", 2L))
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
})

# ---- diagnostics-neighbors ----
local({
# ---- algorithm-knn ----

test_that("explicit KNN neighborhoods reproduce original predictions", {
  skip_if_not_installed("kknn")
  d <- sb_example_binary(n = 100L, seed = 108L)
  d$gene_c <- 2 * d$gene_c + 4
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  for (id in c("knn_binary", "kernel_weighted_knn_binary")) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    e <- unserialize(fit$learner_state$engine_raw)
    for (rows in list(1L, c(1L, 7L, 11L))) {
      query <- d[rows, ]
      if (length(rows) > 1L) query$gene_c <- query$gene_c + 0.17
      z <- sb_knn_neighborhood(fit, query)
      expect_equal(nrow(z$table), length(rows) * e$k)
      expect_equal(z$predictions$probability,
        sb_predict_final(fit, query, "patient_id")$predictions$probability, tolerance = 1e-12)
      expect_equal(as.numeric(tapply(z$table$voting_weight * (z$table$class == "1"),
        z$table$query, sum)), z$predictions$probability, tolerance = 1e-12)
      if (id == "knn_binary") {
        for (i in seq_along(rows)) {
          standardized <- (as.numeric(query[i, features]) - e$center) / e$scale
          distances <- apply(e$training, 1L, function(row) sum((row - standardized)^2))
          chosen <- order(distances, seq_along(distances))[seq_len(e$k)]
          tab <- z$table[z$table$query == i, ]
          expect_identical(tab$training_row, chosen)
          expect_equal(tab$distance, sqrt(unname(distances[chosen])), tolerance = 1e-12)
          expect_equal(tab$voting_weight, rep(1 / e$k, e$k))
        }
      } else {
        training <- e$training
        training$y <- e$outcome
        native <- kknn::kknn(y ~ ., training, query[features], k = e$k,
          distance = e$distance, kernel = e$kernel, scale = e$scale)
        expect_equal(z$table$training_row, as.integer(t(native$C)))
        expect_equal(z$table$distance, as.numeric(t(native$D)))
        expect_equal(z$table$native_weight, as.numeric(t(native$W)))
        expect_identical(z$table$class, as.character(t(native$CL)))
      }
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
      expect_identical(sb_native_diagnostics(fit, newdata = query)$status$status, "available")
    }
    expect_identical(sb_native_diagnostics(fit)$status$status, "unavailable")
    expect_match(sb_native_diagnostics(fit)$status$reason, "explicit query")
    if (id == "knn_binary") {
      tied <- fit
      e$training[,] <- 0
      tied$learner_state$engine_raw <- serialize(e, NULL, version = 2L)
      tie_query <- as.data.frame(as.list(e$center))
      explained <- sb_knn_neighborhood(tied, tie_query)
      expect_identical(explained$table$training_row, seq_len(e$k))
      expect_equal(explained$table$distance, rep(0, e$k))
      expect_equal(explained$predictions$probability, mean(e$outcome[seq_len(e$k)]))
    }
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})

# ---- algorithm-loclda ----

test_that("local-LDA neighborhoods reproduce native local posterior calculations", {
  skip_if_not_installed("klaR")
  d <- sb_example_binary(160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = "loclda_binary"))
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  e <- unserialize(fit$learner_state$engine_raw)
  for (rows in list(8L, c(1L, 8L, 13L))) {
    q <- d[rows, ]; q$gene_c <- q$gene_c + .71
    z <- sb_loclda_neighborhood(fit, q)
    for (i in seq_len(nrow(q))) {
      t <- z$table[z$table$query == i, ]
      distances <- as.numeric(dist(rbind(as.matrix(q[i, features]), e$learn)))[seq_len(nrow(e$learn))]
      ids <- which(distances <= sort(distances)[e$k])
      expect_identical(t$reference_row, ids)
      expect_equal(t$distance, distances[ids])
      expect_equal(t$weight, vapply(distances[ids] / sort(distances)[e$k], e$weight.func, numeric(1)))
      t <- t[t$used, ]
      pooled <- matrix(0, length(features), length(features))
      means <- matrix(0, 2L, length(features))
      prior <- numeric(2L)
      for (j in 1:2) {
        tt <- t[t$class == as.character(j - 1L), ]
        X <- e$learn[tt$reference_row, , drop = FALSE]
        means[j, ] <- colSums(X * tt$weight) / sum(tt$weight)
        pooled <- pooled + cov.wt(X, wt = tt$weight)$cov * nrow(X)
        prior[j] <- sum(tt$weight) / sum(t$weight)
      }
      pooled <- pooled / (nrow(t) - 2L)
      delta <- sweep(means, 2L, as.numeric(q[i, features]))
      score <- -.5 * rowSums((delta %*% solve(pooled)) * delta) + log(prior)
      expect_equal(plogis(score[2L] - score[1L]), z$predictions$probability[i], tolerance = 1e-12)
    }
    native <- get("predict.loclda", asNamespace("klaR"))(e, as.matrix(q[features]))$posterior[, 2L]
    expect_equal(z$predictions$probability, as.numeric(native))
    expect_equal(z$predictions$probability, sb_predict_final(fit, q, "patient_id")$predictions$probability)
    built <- ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z))
    expect_s3_class(built, "ggplot_built")
    expect_identical(as.character(built$plot$scales$get_scales("shape")$get_labels()), "Used")
    expect_identical(sb_native_diagnostics(fit, newdata = q)$status$status, "available")
  }
  expect_identical(sb_native_diagnostics(fit)$status$status, "unavailable")
  expect_identical(serialize(fit, NULL, version = 2L), before)
})

# ---- algorithm-lvq ----

test_that("all LVQ diagnostic views retain original codebooks and prediction geometry", {
  skip_if_not_installed("class")
  d <- sb_example_binary(n = 160L, seed = 108L)
  d$gene_a <- d$gene_a + 10
  d$gene_b <- 4 * d$gene_b - 5
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  for (id in sb_lvq_binary_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success", info = id)
    expect_null(fit$diagnostic_artifacts)
    before <- serialize(fit, NULL, version = 2L)
    prediction_before <- sb_predict_final(fit, d, "patient_id")
    e <- unserialize(fit$learner_state$engine_raw)
    diagnostic <- sb_lvq_diagnostics(fit)
    expect_identical(diagnostic$table$value, as.vector(e$prototypes))
    expect_equal(diagnostic$table$original_value,
      as.vector(sweep(sweep(e$prototypes, 2L, e$scale, "*"), 2L, e$center, "+")), tolerance = 0)
    expect_identical(diagnostic$table$feature, rep(features, each = nrow(e$prototypes)))
    expect_identical(diagnostic$table$class, rep(e$prototype_class, length(features)))
    expect_equal(diagnostic$prototype_counts$count,
      as.integer(table(factor(e$prototype_class, levels = c("0", "1")))))
    expect_identical(diagnostic$iterations, e$iterations)
    expect_identical(names(diagnostic$center), features)
    expect_equal(unname(diagnostic$scale), e$scale)
    distances <- sb_lvq_diagnostics(fit, "prototype_distances")
    expected <- vapply(seq_len(nrow(e$prototypes)), function(i)
      rowSums(sweep(e$prototypes, 2L, e$prototypes[i, ], "-")^2), numeric(nrow(e$prototypes)))
    expect_equal(distances$table$value, as.vector(expected), tolerance = 1e-12)
    expect_equal(sum(diag(matrix(distances$table$value, nrow(e$prototypes)))), 0)
    expect_error(sb_lvq_diagnostics(fit, "training_loss"), class = "stormbreaker_unsupported_diagnostic")
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 2L))
    if (requireNamespace("ggplot2", quietly = TRUE)) for (view in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(view)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_null(p$labels$subtitle)
      expect_null(p$labels$caption)
    }
    expect_identical(sb_predict_final(fit, d, "patient_id"), prediction_before)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
})

# ---- diagnostics-omics ----
local({
# ---- algorithm-cancerclass ----

test_that("Cancerclass views preserve centroid order and the saved probability map", {
  skip_if_not_installed("cancerclass")
  skip_if_not_installed("Biobase")
  d <- sb_example_binary(160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = "cancerclass_binary"))
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  centroids <- sb_cancerclass_diagnostics(fit)
  mapping <- sb_cancerclass_diagnostics(fit, "score_mapping", c(-2, 3))
  m <- matrix(centroids$table$mean, nrow = length(features))
  expect_identical(centroids$table$feature, rep(features, 2L))
  expect_equal(m, unname(fit$learner_state$predictor[features, c("0", "1")]))
  expect_identical(mapping$calibration_coefficients, fit$learner_state$calibration_coefficients)
  expect_equal(range(mapping$table$score), c(-2, 3))
  b <- mapping$calibration_coefficients
  expect_equal(mapping$table$probability, plogis(b[1L] + b[2L] * mapping$table$score))
  q <- d[c(1L, 8L, 13L), features]; q$gene_c <- q$gene_c + .71
  z <- (cancerclass::get.d(t(as.matrix(q)), m[, 1L], method = "cor") -
    cancerclass::get.d(t(as.matrix(q)), m[, 2L], method = "cor")) /
    cancerclass::get.d(m[, 1L], m[, 2L], method = "cor")
  expect_equal(as.numeric(plogis(b[1L] + b[2L] * z)),
    sb_predict_final(fit, data.frame(patient_id = 1:3, q), "patient_id")$predictions$probability,
    tolerance = 1e-12)
  for (view in list(centroids, mapping)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(view)), "ggplot_built")
  expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status, rep("available", 2L))
  expect_identical(serialize(fit, NULL, version = 2L), before)
})

# ---- algorithm-high-dimensional ----

test_that("HDA and HDDA original geometry reproduces native probabilities", {
  skip_if_not_installed("hda")
  skip_if_not_installed("HDclassif")
  d <- sb_example_binary(160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in c("hda_binary", "hdda_binary_cattell")) {
    fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    e <- fit$learner_state$engine
    loadings <- sb_high_dimensional_diagnostics(fit)
    q <- d[c(1L, 8L, 13L), ]; q$gene_c <- q$gene_c + .71
    mass <- matrix(0, nrow(q), 2L)
    if (id == "hda_binary") {
      densities <- sb_high_dimensional_diagnostics(fit, "densities")
      m <- matrix(loadings$table$loading, nrow = length(features))
      expect_equal(m, unname(e$hda.loadings[, seq_len(e$reduced.dimension), drop = FALSE]))
      expect_identical(loadings$table$feature, rep(features, e$reduced.dimension))
      projected <- as.matrix(q[features]) %*% m
      native_projected <- get("predict.hda", asNamespace("hda"))(e, as.matrix(q[features]), task = "dr")
      expect_equal(unname(projected), unname(native_projected))
      for (i in 1:2) {
        par <- densities$parameters[densities$parameters$class == as.character(i - 1L), ]
        expect_equal(par$mean, unname(e$naivebayes$tables[[1L]][i, 1L]))
        mass[, i] <- dnorm(projected[, 1L], par$mean, par$sd) * as.numeric(e$naivebayes$apriori[i])
      }
      views <- list(loadings, densities)
    } else {
      spectrum <- sb_high_dimensional_diagnostics(fit, "spectrum")
      centers <- sb_high_dimensional_diagnostics(fit, "centers")
      for (i in 1:2) {
        cls <- e$kname[i]
        rows <- loadings$table$class == cls
        Q <- matrix(loadings$table$loading[rows], nrow = length(features))
        expect_equal(Q, e$Q[[i]])
        spec <- spectrum$table[spectrum$table$class == cls, ]
        expect_equal(spec$eigenvalue, as.numeric(e$ev[i, ]))
        expect_equal(sum(spec$retained), as.integer(e$d[i]))
        k <- ncol(Q)
        a <- spec$model_variance[seq_len(k)]
        b <- spec$model_variance[k + 1L]
        cov <- Q %*% diag(a - b, k) %*% t(Q) + diag(b, length(features))
        cov <- cov * outer(e$scaling$sd, e$scaling$sd)
        mu <- centers$table$mean[centers$table$class == cls]
        expect_equal(mu, as.numeric(e$mu[i, ]) * as.numeric(e$scaling$sd) + as.numeric(e$scaling$mu))
        delta <- sweep(as.matrix(q[features]), 2L, mu)
        mass[, i] <- exp(-rowSums((delta %*% solve(cov)) * delta) / 2) /
          sqrt((2 * pi)^length(features) * det(cov)) * as.numeric(e$prop[i])
      }
      views <- list(loadings, spectrum, centers)
    }
    expect_equal(as.numeric(mass[, 2L] / rowSums(mass)), sb_predict_final(fit, q, "patient_id")$predictions$probability, tolerance = 1e-12)
    for (z in views) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status, rep("available", length(views)))
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})

# ---- algorithm-ktsp ----

test_that("k-TSP diagnostics preserve ordered pair rules and tie votes", {
  skip_if_not_installed("switchBox")
  d <- withr::with_seed(1701L, {
    x <- matrix(rnorm(72 * 8), nrow = 72)
    colnames(x) <- paste0("gene_", seq_len(8))
    y <- as.integer(x[, 1] - .7 * x[, 2] + .3 * x[, 3] + rnorm(72, sd = .6) > 0)
    data.frame(patient_id = paste0("NM", seq_len(72)), outcome = y, x)
  })
  features <- paste0("gene_", c(3, 1, 2, 4:8))
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features),
    sb_pipeline(learner = "ktsp_binary"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  e <- fit$learner_state$engine
  view <- sb_ktsp_diagnostics(fit)
  expect_identical(view$table$first, e$pairs[, 1])
  expect_identical(view$table$second, e$pairs[, 2])
  expect_equal(view$table$native_score, e$pair_scores, tolerance = 0)
  expect_equal(view$table$pair, seq_len(e$selected_pairs))
  expect_identical(view$unused_features, setdiff(features, as.vector(e$pairs)))
  expect_identical(view$vote_aggregation, "unweighted_mean")
  expect_equal(unname(view$vote_rule), c(1, .5, 0), tolerance = 0)
  new <- d[1:3, , drop = FALSE]
  for (feature in features) new[[feature]] <- 0
  for (first in e$pairs[, 1]) new[[first]] <- c(1, 0, -1)
  expect_equal(sb_predict_final(fit, new, "patient_id")$predictions$probability, c(1, .5, 0), tolerance = 0)
  compact <- fit
  compact$diagnostic_artifacts <- NULL
  expect_identical(sb_ktsp_diagnostics(compact), view)
  expect_identical(sb_native_diagnostics(compact)$status$status, "available")
  p <- sb_plot_algorithm_diagnostic(view)
  expect_identical(levels(p$data$rule), rev(paste0("Pair ", view$table$pair, ": ", view$table$class_1_rule)))
  expect_equal(ggplot2::ggplot_build(p)$data[[1]]$x, view$table$native_score, tolerance = 0)
  expect_null(p$labels$subtitle)
  expect_null(p$labels$caption)
  expect_identical(serialize(fit, NULL, version = 2L), before)
})

test_that("a multi-pair k-TSP fit keeps every rule and averages rather than weights votes", {
  skip_if_not_installed("switchBox")
  d <- withr::with_seed(1001L, {
    y <- sample(rep(c(0L, 1L), 120L))
    x <- matrix(rnorm(240 * 12), nrow = 240)
    for (j in seq(1, 12, 2)) x[, j] <- x[, j] + .9 * y
    colnames(x) <- paste0("gene_", seq_len(12))
    data.frame(patient_id = paste0("K", seq_len(240)), outcome = y, x)
  })
  features <- paste0("gene_", seq_len(12))
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = "ktsp_binary"))
  view <- sb_ktsp_diagnostics(fit)
  expect_gt(nrow(view$table), 1L)
  expect_equal(view$table$native_score, fit$learner_state$engine$pair_scores, tolerance = 0)
  expect_equal(nrow(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(view))$data[[1]]), nrow(view$table))
  new <- d[1:3, , drop = FALSE]
  for (feature in features) new[[feature]] <- 0
  for (first in view$table$first) new[[first]] <- c(-1, 0, 1)
  new[[view$table$first[1]]][1] <- 1
  expect_equal(sb_predict_final(fit, new, "patient_id")$predictions$probability,
    c(1 / nrow(view$table), .5, 1), tolerance = 1e-12)
})

# ---- algorithm-superpc ----

test_that("SuperPC screening, loadings and CV are read from the original objects", {
  skip_if_not_installed("superpc")
  for (id in c("superpc", "superpc_cv")) {
    pipeline <- if (id == "superpc_cv") sb_pipeline(learner = id,
      learner_tuning_folds = 3L, learner_seed = 8903L) else sb_pipeline(learner = id)
    fit <- sb_fit_final(make_component_survival_task(), pipeline, capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    original <- serialize(fit, NULL)
    model <- fit$diagnostic_artifacts$learner
    engine <- model$engine
    screen <- sb_superpc_diagnostics(fit)
    loadings <- sb_superpc_diagnostics(fit, "loadings")
    expect_identical(screen$native_table, engine$trained$feature.scores)
    expect_equal(screen$table$selected, unname(engine$prediction$which.features))
    expect_identical(loadings$native_table, engine$prediction$u)
    expect_equal(loadings$table$loading, as.vector(engine$prediction$u))
    expect_equal(loadings$table$feature, model$features[engine$prediction$which.features])
    expect_equal(screen$table$selected, abs(screen$table$score) >= screen$threshold)
    if (requireNamespace("ggplot2", quietly = TRUE)) for (q in list(screen, loadings)) expect_s3_class(
      ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    cv <- tryCatch(sb_superpc_diagnostics(fit, "cv_curve"), error = identity)
    if (id == "superpc_cv") {
      native <- engine$cross_validation
      expect_equal(cv$table$threshold, rep(native$thresholds, each = nrow(native$scor)))
      expect_equal(cv$table$score, as.vector(native$scor))
      expect_equal(cv$table$lower, as.vector(native$scor.lower))
      expect_equal(cv$table$upper, as.vector(native$scor.upper))
      expect_equal(cv$table$component[cv$table$selected], 1L)
      expect_equal(cv$table$threshold[cv$table$selected], native$thresholds[which.max(native$scor[1, ])])
      expect_equal(unname(cv$folds), c(model$audit$requested_folds, native$n.fold))
      if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(
        ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(cv)), "ggplot_built")
    } else expect_s3_class(cv, "stormbreaker_unavailable_diagnostic")
    expect_identical(serialize(fit, NULL), original)
    fit$diagnostic_artifacts <- NULL
    err <- tryCatch(sb_superpc_diagnostics(fit), error = identity)
    expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
  }
})

test_that("SuperPC preserves its one-selected-feature singular vector", {
  skip_if_not_installed("superpc")
  fit <- sb_fit_final(make_component_survival_task(),
    sb_pipeline(learner = "superpc", superpc_threshold = 2.5), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  q <- sb_superpc_diagnostics(fit, "loadings")
  expect_equal(nrow(q$table), 1L)
  expect_identical(q$native_table, fit$diagnostic_artifacts$learner$engine$prediction$u)
  expect_equal(q$table$feature, "gene_01")
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(
    ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
})
})

# ---- diagnostics-projections ----
local({
# ---- algorithm-pls-cv ----

test_that("PLS-Cox CV plots the original higher-is-better criterion without capture", {
  skip_if_not_installed("plsRcox")
  task <- make_component_survival_task()
  fit <- sb_fit_final(task, sb_pipeline(learner = "plsrcox_cv", learner_tuning_folds = 3L,
    learner_seed = 8903L))
  expect_identical(fit$status, "success")
  original <- serialize(fit, NULL, version = 2L)
  q <- sb_pls_diagnostics(fit, "cv_curve")
  a <- fit$learner_state$audit
  expect_null(fit$diagnostic_artifacts)
  expect_identical(q$table$iauc, a$cv_error)
  expect_identical(q$table$native_spread, a$cv_standard_error)
  expect_identical(q$table$component, seq.int(0L, a$maximum_components))
  expect_equal(q$table$component[q$table$selected], a$selected_components)
  expect_equal(q$selected_components, which.max(q$table$iauc) - 1L)
  expect_identical(q$folds, 3L)
  p <- sb_plot_algorithm_diagnostic(q)
  expect_match(p$labels$y, "higher is better", fixed = TRUE)
  expect_match(p$labels$caption, "not conventional SE", fixed = TRUE)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  expect_true(identical(serialize(fit, NULL, version = 2L), original))
  # Unavailable components must not be imputed or cause the remaining curve to
  # be reselected; synthetic corruption and missingness exercise this boundary.
  missing <- fit
  missing$learner_state$audit$cv_error[2L] <- NA_real_
  missing$learner_state$audit$cv_standard_error[2L] <- NA_real_
  q <- sb_pls_diagnostics(missing, "cv_curve")
  expect_true(is.na(q$table$iauc[2L]))
  expect_identical(q$selected_components, a$selected_components)
  expect_warning(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), NA)
  bad <- fit
  bad$learner_state$audit$cv_standard_error[1L] <- -1
  expect_error(sb_pls_diagnostics(bad, "cv_curve"), class = "stormbreaker_invalid_diagnostic")
  bad <- fit
  bad$learner_state$audit$cv_error <- 0.5
  expect_error(sb_pls_diagnostics(bad, "cv_curve"), class = "stormbreaker_invalid_diagnostic")
  single <- fit
  single$learner_state$audit$component_selection <- "single_predictor_has_one_feasible_component"
  expect_error(sb_pls_diagnostics(single, "cv_curve"), class = "stormbreaker_unavailable_diagnostic")
  fixed <- sb_fit_final(task, sb_pipeline(learner = "plsrcox"))
  expect_identical(fixed$status, "success")
  expect_error(sb_pls_diagnostics(fixed, "cv_curve"), class = "stormbreaker_unavailable_diagnostic")
})

# ---- algorithm-pls ----

test_that("five PLS learners retain exact native components without enlarging portable state", {
  skip_if_not_installed("mixOmics")
  skip_if_not_installed("plsRglm")
  skip_if_not_installed("plsRcox")
  d <- make_sb_binary_fixture(n = 100L, seed = 108L)
  binary <- sb_binary_task(d, "patient_id", "outcome", c("gene_c", "gene_a", "gene_b"))
  for (id in c("plsrglm_binary", "plsda_binary", "splsda_binary", "plsrcox", "plsrcox_cv")) {
    task <- if (grepl("binary", id)) binary else make_component_survival_task()
    pipeline <- if (id == "plsrcox_cv") sb_pipeline(learner = id, learner_tuning_folds = 3L,
      learner_seed = 8903L) else sb_pipeline(learner = id)
    fit <- sb_fit_final(task, pipeline, capture_diagnostics = TRUE)
    compact <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    expect_identical(compact$status, "success")
    expect_identical(fit$learner_state, compact$learner_state)
    expect_null(fit$learner_state$diagnostic_components)
    # Version 2 compares the full saved state without ALTREP representation
    # changes when native deferred score row names are materialized by unname().
    original <- serialize(fit, NULL, version = 2L)
    loadings <- sb_pls_diagnostics(fit)
    scores <- sb_pls_diagnostics(fit, "score_distribution")
    model <- fit$diagnostic_artifacts$learner
    native <- if (startsWith(id, "plsrcox")) list(loadings = model$engine$pp,
      scores = model$engine$tt) else model$diagnostic_components
    expect_identical(loadings$native_table, native$loadings)
    expect_identical(scores$native_table, unname(native$scores))
    expect_equal(loadings$table$loading, as.vector(native$loadings))
    expect_equal(scores$table$score, as.vector(native$scores))
    expect_equal(loadings$table$feature[loadings$table$component == 1L], model$features)
    expect_identical(scores$sensitive_training_data, TRUE)
    expect_null(rownames(scores$native_table))
    if (requireNamespace("ggplot2", quietly = TRUE)) for (q in list(loadings, scores)) expect_s3_class(
      ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    expect_true(identical(serialize(fit, NULL, version = 2L), original), info = id)
    err <- tryCatch(sb_pls_diagnostics(compact), error = identity)
    expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
  }
})

test_that("PLS diagnostics preserve sparse zeros and actual one-component fits", {
  skip_if_not_installed("mixOmics")
  set.seed(108L)
  x <- as.data.frame(matrix(rnorm(100L * 12L), 100L))
  names(x) <- paste0("gene_", seq_len(12L))
  d <- data.frame(id = seq_len(100L), outcome = rbinom(100L, 1L,
    plogis(x[[1]] - x[[2]])), x)
  task <- sb_binary_task(d, "id", "outcome", names(x))
  fit <- sb_fit_final(task, sb_pipeline(learner = "splsda_binary"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  q <- sb_pls_diagnostics(fit)
  expect_equal(nrow(q$table), 24L)
  expect_true(any(q$table$loading == 0))
  expect_identical(q$native_table, fit$diagnostic_artifacts$learner$diagnostic_components$loadings)
  task <- sb_binary_task(d, "id", "outcome", "gene_1")
  fit <- sb_fit_final(task, sb_pipeline(learner = "plsda_binary"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  q <- sb_pls_diagnostics(fit)
  expect_identical(unique(q$table$component), 1L)
  expect_equal(nrow(q$table), 1L)
  fit$diagnostic_artifacts$learner$diagnostic_components <- NULL
  expect_error(sb_pls_diagnostics(fit), class = "stormbreaker_missing_diagnostic_artifacts")
})
})

# ---- diagnostics-smooth ----
local({
# ---- algorithm-bart ----

test_that("BART posterior views use original draws and original score scales", {
  skip_if_not_installed("BART")
  for (id in c("bart_binary", "bart_survival", "bart_aft_survival")) {
    d <- if (id == "bart_binary") sb_example_binary(100L, seed = 108L) else sb_example_survival(100L, seed = 8612L)
    features <- c("gene_c", "gene_a", "gene_b")
    task <- if (id == "bart_binary") sb_binary_task(d, "patient_id", "outcome", features) else
      sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    query <- d[c(1L, 7L), ]
    query$gene_c <- query$gene_c + .17
    z <- sb_bart_posterior(fit, query)
    e <- if (is.null(fit$learner_state$engine_raw)) fit$learner_state$engine else unserialize(fit$learner_state$engine_raw)
    matrix_x <- as.matrix(query[features])
    if (id == "bart_aft_survival") {
      invisible(capture.output(native_draws <- BART::pwbart(matrix_x[, e$kept_columns, drop = FALSE],
        e$treedraws, mu = e$offset, mc.cores = 1L, dodraws = TRUE)))
    } else {
      if (id == "bart_survival") matrix_x <- cbind(t = rep(e$times, 2),
        matrix_x[rep(1:2, each = length(e$times)), , drop = FALSE])
      if (length(e$rm.const)) matrix_x <- matrix_x[, e$rm.const, drop = FALSE]
      fn <- get(if(id == "bart_binary") "predict.pbart" else "predict.survbart", asNamespace("BART"))
      invisible(capture.output(original <- fn(e, matrix_x, mc.cores = 1L, openmp = FALSE)))
      native_draws <- if(id == "bart_binary") original$prob.test else original$surv.test
    }
    expect_equal(z$draws, native_draws, tolerance = 1e-12)
    expect_equal(z$table$mean, unname(colMeans(z$draws)))
    expect_equal(z$table$lower, unname(apply(z$draws, 2L, quantile, .025)))
    expect_equal(z$table$upper, unname(apply(z$draws, 2L, quantile, .975)))
    if (id == "bart_binary") {
      expect_equal(z$table$mean, sb_predict_final(fit, query, "patient_id")$predictions$probability, tolerance = 1e-12)
    } else if (id == "bart_aft_survival") {
      expect_equal(-z$table$mean, sb_predict_final(fit, query, "patient_id")$predictions$risk, tolerance = 1e-12)
    } else {
      model <- sb_restore_additional_survival_model(fit$learner_state, id)
      native <- sb_bart_survival_curves(model, query[features])
      expect_equal(z$table$mean, as.numeric(t(native$survival)), tolerance = 1e-12)
      expect_equal(z$table$time, rep(native$time, 2))
    }
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    expect_identical(sb_native_diagnostics(fit)$status$status, "unavailable")
    expect_identical(sb_native_diagnostics(fit, newdata = query)$status$status, "available")
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})

# ---- algorithm-earth ----

test_that("earth diagnostics retain logistic rather than least-squares coefficients", {
  skip_if_not_installed("earth")
  d <- sb_example_binary(n = 160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  fit <- sb_fit_final(task, sb_pipeline(learner = "earth_binary"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  native <- fit$diagnostic_artifacts$learner$engine
  diagnostic <- sb_earth_diagnostics(fit)
  expect_equal(diagnostic$table$coefficient, as.numeric(native$glm.coefficients), tolerance = 0)
  expect_false(isTRUE(all.equal(diagnostic$table$coefficient, as.numeric(native$coefficients))))
  expect_identical(diagnostic$table$term, rownames(native$glm.coefficients))
  expect_equal(diagnostic$table$basis_index, native$selected.terms)
  expect_equal(diagnostic$cuts, native$cuts, tolerance = 0)
  expect_equal(diagnostic$directions, native$dirs, tolerance = 0)
  basis <- get("model.matrix.earth", asNamespace("earth"))(native, d[features])
  expect_equal(plogis(drop(basis %*% diagnostic$table$coefficient)),
    sb_predict_final(fit, d, "patient_id")$predictions$probability, tolerance = 1e-12)
  shap <- sb_additive_shap(fit, d, "patient_id", d)
  expect_equal(unname(rowSums(shap$contributions)) + shap$reference_value,
    drop(basis %*% diagnostic$table$coefficient), tolerance = 1e-12)
  expect_identical(colnames(shap$contributions), features)
  expect_identical(shap$score_scale, "class_1_log_odds")
  pruning <- sb_earth_diagnostics(fit, "pruning_gcv")
  expect_equal(pruning$table$gcv, as.numeric(native$gcv.per.subset), tolerance = 0)
  expect_equal(pruning$table$terms[pruning$table$selected], length(native$selected.terms))
  compact <- fit
  compact$diagnostic_artifacts <- NULL
  expect_equal(sb_earth_diagnostics(compact)$table, diagnostic$table)
  expect_error(sb_earth_diagnostics(compact, "pruning_gcv"), class = "stormbreaker_missing_diagnostic_artifacts")
  expect_identical(sb_native_diagnostics(compact, complete = TRUE)$status$status, c("available", "missing_capture"))
  expect_error(sb_nomogram(fit, d), class = "stormbreaker_unsupported_coefficients")
  bundle <- sb_native_diagnostics(fit, complete = TRUE)
  expect_identical(bundle$status$status, rep("available", 2L))
  if (requireNamespace("ggplot2", quietly = TRUE)) for (view in bundle$diagnostics) {
    p <- sb_plot_algorithm_diagnostic(view)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    expect_null(p$labels$subtitle)
    expect_null(p$labels$caption)
  }
  expect_identical(serialize(fit, NULL, version = 2L), before)
})

# ---- algorithm-neural ----

test_that("neural diagnostics preserve every native connection and forward prediction", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("neuralnet")
  d <- sb_example_binary(n = 160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  for (id in c("nnet_binary", "neuralnet_multilayer_binary")) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success", info = id)
    before <- serialize(fit, NULL, version = 2L)
    diagnostic <- sb_neural_diagnostics(fit)
    parameters <- diagnostic$parameters
    e <- unserialize(fit$learner_state$engine_raw)
    if (id == "nnet_binary") {
      expect_equal(unlist(parameters$weights, use.names = FALSE), e$model$wts, tolerance = 0)
    } else expect_equal(parameters$weights, e$model$weights[[1L]], tolerance = 0)
    expect_equal(nrow(diagnostic$table), sum(vapply(parameters$weights, length, integer(1))))
    expect_equal(sum(diagnostic$nodes$bias), length(parameters$weights))
    expect_identical(diagnostic$nodes$label[diagnostic$nodes$layer == 0 & !diagnostic$nodes$bias], features)
    value <- sweep(sweep(as.matrix(d[features]), 2L, parameters$center, "-"), 2L, parameters$scale, "/")
    for (weights in parameters$weights) {
      score <- cbind(1, value) %*% weights
      value <- if (id == "nnet_binary") ifelse(score < -15, 0,
        ifelse(score > 15, 1, 1 / (1 + exp(-score)))) else e$model$act.fct(score)
    }
    expect_equal(as.numeric(value), sb_predict_final(fit, d, "patient_id")$predictions$probability, tolerance = 1e-12)
    expect_error(sb_neural_diagnostics(fit, "learning_curve"), class = "stormbreaker_unsupported_diagnostic")
    expect_error(sb_nomogram(fit, d), class = "stormbreaker_unsupported_coefficients")
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 2L))
    if (requireNamespace("ggplot2", quietly = TRUE)) for (view in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(view)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_null(p$labels$subtitle)
      expect_null(p$labels$caption)
    }
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
})

# ---- diagnostics-survival ----
local({
# ---- algorithm-coxboost ----

make_coxboost_diagnostic_task <- function() {
  set.seed(8901L)
  x <- as.data.frame(matrix(stats::rnorm(90L * 12L), nrow = 90L))
  names(x) <- sprintf("gene_%02d", seq_len(12L))
  lp <- x[[1]] - x[[2]] + .6 * x[[3]]
  et <- stats::rexp(90L, rate = exp(lp) / 15)
  ct <- stats::rexp(90L, rate = 1 / 25)
  d <- data.frame(id = seq_len(90L), time = pmax(pmin(et, ct), .05),
    event = as.integer(et <= ct), x)
  sb_surv_task(d, "id", "time", "event", rev(names(x)), "month")
}

test_that("both CoxBoost variants expose their exact captured coefficient paths", {
  skip_if_not_installed("CoxBoost")
  for (id in c("coxboost", "coxboost_cv")) {
    pipeline <- if (id == "coxboost_cv") sb_pipeline(learner = id,
      learner_tuning_folds = 3L, learner_seed = 8902L) else sb_pipeline(learner = id)
    fit <- sb_fit_final(make_coxboost_diagnostic_task(), pipeline, capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    original <- serialize(fit, NULL)
    model <- fit$diagnostic_artifacts$learner
    engine <- model$engine
    q <- sb_boosting_diagnostics(fit, "coefficient_path")
    native <- as.matrix(stats::coef(engine, at.step = seq.int(0L, engine$stepno), scaled = FALSE))
    expect_equal(unname(q$native_table), unname(native))
    expected <- sweep(native, 2L, engine$sdx, "/")
    expect_equal(q$table$coefficient, as.vector(expected))
    endpoint <- q$table[q$table$step == engine$stepno, ]
    expect_identical(endpoint$feature, fit$learner_state$features)
    expect_equal(endpoint$coefficient, fit$learner_state$coefficient)
    centered <- sweep(as.matrix(fit$diagnostic_artifacts$predictors), 2L, engine$meanx, "-")
    expect_equal(as.numeric(centered %*% endpoint$coefficient),
      as.numeric(engine$linear.predictor[nrow(engine$linear.predictor), ]), tolerance = 1e-10)
    count <- sb_boosting_diagnostics(fit, "nonzero_path")
    expect_equal(count$table$nonzero, rowSums(native != 0))
    expect_equal(count$table$step, seq.int(0L, engine$stepno))
    expect_identical(serialize(fit, NULL), original)
    if (requireNamespace("ggplot2", quietly = TRUE)) for (z in list(q, count)) {
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    }
    fit$diagnostic_artifacts <- NULL
    err <- tryCatch(sb_boosting_diagnostics(fit, "coefficient_path"), error = identity)
    expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
    cv <- tryCatch(sb_boosting_diagnostics(fit, "cv_curve"), error = identity)
    if (id == "coxboost_cv") {
      expect_identical(cv$table, model$audit$tuning_curve)
      expect_equal(cv$stopping_iteration, model$audit$optimal_step)
      expect_equal(cv$stopping_iteration, cv$table$step[which.max(cv$table$mean_log_likelihood)])
      if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(
        ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(cv)), "ggplot_built")
    } else expect_s3_class(cv, "stormbreaker_unavailable_diagnostic")
  }
})

test_that("single-feature CoxBoost paths retain a matrix and all zero-step rows", {
  skip_if_not_installed("CoxBoost")
  d <- make_sb_fixture(n = 100L, seed = 8612L)
  fit <- sb_fit_final(sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month"),
    sb_pipeline(learner = "coxboost", coxboost_stepno = 1L), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  q <- sb_boosting_diagnostics(fit, "coefficient_path")
  expect_equal(dim(q$native_table), c(2L, 1L))
  expect_equal(q$table$step, c(0L, 1L))
  expect_equal(q$table$coefficient[1], 0)
  expect_equal(q$table$coefficient[2], fit$learner_state$coefficient)
  broken <- fit
  broken$learner_state$coefficient <- broken$learner_state$coefficient + 1
  err <- tryCatch(sb_boosting_diagnostics(broken, "coefficient_path"), error = identity)
  expect_s3_class(err, "stormbreaker_invalid_diagnostic")
})

# ---- algorithm-penaft ----

test_that("penAFT diagnostics preserve native CV and input-unit path predictions", {
  skip_if_not_installed("penAFT")
  d <- sb_example_survival(100L, seed = 8612L)
  d$gene_c <- d$gene_c * 2 + 4
  features <- c("gene_c", "gene_a", "gene_b")
  x <- as.matrix(d[features])
  fold <- sb_make_stratified_fold_id(d$event, 5L, 1L, stratum_label = "event")
  invisible(capture.output(native <- withr::with_seed(1L, penAFT::penAFT.cv(
    X = x, logY = log(d$time), delta = d$event, nlambda = 50L,
    lambda.ratio.min = .1, penalty = "EN", alpha = .5, tol.abs = 1e-8,
    tol.rel = 2.5e-4, standardize = TRUE,
    cv.index = lapply(1:5, function(i) which(fold == i)), admm.max.iter = 10000L, quiet = TRUE))))
  fit <- sb_fit_final(sb_surv_task(d, "patient_id", "time", "event", features, "month"),
    sb_pipeline(learner = "penaft_gehan_elastic_net"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  cv <- sb_penaft_diagnostics(fit, "cv_curve")
  path <- sb_penaft_diagnostics(fit, "coefficient_path")
  co <- sb_penaft_diagnostics(fit)
  f <- native$full.fit
  expect_equal(cv$table$lambda, f$lambda, tolerance = 0)
  expect_equal(cv$table$score, native$cv.err.linPred, tolerance = 0)
  expect_equal(which(cv$table$selected), which.min(native$cv.err.linPred))
  expect_equal(path$native_path$full_fit$beta, f$beta, tolerance = 0)
  raw <- sweep(as.matrix(f$beta), 1, f$X.sd, "/")
  expect_equal(path$table$coefficient, as.numeric(t(raw)), tolerance = 0)
  expect_equal(path$table$feature, rep(features, each = length(f$lambda)))
  expect_equal(path$intercept_path, -as.numeric(crossprod(f$X.mean, raw)), tolerance = 0)
  native_lp <- scale(x, center = f$X.mean, scale = f$X.sd) %*% as.matrix(f$beta)
  replay <- sweep(x %*% raw, 2, path$intercept_path, "+")
  expect_equal(unname(replay), unname(native_lp), tolerance = 1e-12)
  j <- which.min(native$cv.err.linPred)
  expect_equal(co$table$coefficient, as.numeric(raw[, j]), tolerance = 0)
  expect_equal(co$intercept, path$intercept_path[j], tolerance = 0)
  expect_null(fit$learner_state$diagnostic_path)
  compact <- fit
  compact$diagnostic_artifacts <- NULL
  expect_identical(sb_penaft_diagnostics(compact), co)
  expect_identical(sb_native_diagnostics(compact, complete = TRUE)$status$status,
    c("available", "missing_capture", "missing_capture"))
  bundle <- sb_native_diagnostics(fit, complete = TRUE)
  expect_identical(bundle$status$status, rep("available", 3L))
  for (v in bundle$diagnostics) {
    p <- sb_plot_algorithm_diagnostic(v)
    b <- ggplot2::ggplot_build(p)
    expect_s3_class(b, "ggplot_built")
    expect_null(p$labels$caption)
    expect_null(p$labels$subtitle)
    if (v$type == "penaft_coefficients") expect_equal(b$data[[2]]$x, v$table$coefficient, tolerance = 0)
    if (v$type == "penaft_cv_curve") expect_equal(b$data[[2]]$x, log10(f$lambda[j]), tolerance = 0)
    if (v$type == "penaft_coefficient_path") expect_equal(sort(b$data[[1]]$y), sort(v$table$coefficient), tolerance = 0)
  }
  expect_identical(serialize(fit, NULL, version = 2L), before)
})
})

# ---- diagnostics-trees ----
local({
# ---- algorithm-c50 ----

test_that("C5.0 native feature summaries and trial sizes retain original values", {
  skip_if_not_installed("C50")
  d <- sb_example_binary(160L,seed=108L); d$gene_a <- d$gene_a+1.5*d$outcome
  task <- sb_binary_task(d,"patient_id","outcome",c("gene_c","gene_a","gene_b"))
  for(id in c("c50_binary","c50_rules_binary")) {
    fit <- sb_fit_final(task,sb_pipeline(learner=id),capture_diagnostics=TRUE)
    expect_identical(fit$status,"success")
    before <- serialize(fit,NULL,version=2L)
    prediction <- sb_predict_final(fit,d,"patient_id")
    e <- fit$diagnostic_artifacts$learner$engine
    for(type in c("usage","splits","size")) {
      z <- sb_c50_diagnostics(fit,type)
      if(type=="size") expect_equal(z$table$size,e$size) else
        expect_equal(z$table$value,C50::C5imp(e,metric=type)[fit$learner_state$features,"Overall"])
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)),"ggplot_built")
    }
    expect_identical(sb_native_diagnostics(fit,complete=TRUE)$status$status,rep("available",3))
    compact <- fit; compact$diagnostic_artifacts <- NULL
    expect_identical(sb_native_diagnostics(compact,complete=TRUE)$status$status,rep("missing_capture",3))
    expect_identical(sb_predict_final(compact,d,"patient_id"),prediction)
    expect_identical(serialize(fit,NULL,version=2L),before)
  }
})

test_that("unsplit C5.0 trials do not invent feature importance", {
  skip_if_not_installed("C50")
  d <- sb_example_binary(160L,seed=108L)
  fit <- sb_fit_final(sb_binary_task(d,"patient_id","outcome",c("gene_c","gene_a","gene_b")),
    sb_pipeline(learner="c50_binary"),capture_diagnostics=TRUE)
  expect_identical(fit$status,"success")
  expect_identical(sb_native_diagnostics(fit,complete=TRUE)$status$status,c("unavailable","unavailable","available"))
  expect_equal(sb_c50_diagnostics(fit,"size")$table$size,rep(1,10))
})

# ---- algorithm-classification-tree ----

test_that("tree-package diagnostics retain native nodes, rules and probabilities", {
  skip_if_not_installed("tree")
  d <- sb_example_binary(160L,seed=108L)
  d$gene_a <- d$gene_a + 1.5*d$outcome
  features <- c("gene_c","gene_a","gene_b")
  fit <- sb_fit_final(sb_binary_task(d,"patient_id","outcome",features),sb_pipeline(learner="tree_binary"))
  expect_identical(fit$status,"success")
  before <- serialize(fit,NULL,version=2L)
  e <- unserialize(fit$learner_state$engine_raw)
  f <- e$frame
  z <- sb_classification_tree_diagnostics(fit)
  expect_identical(z$native_frame,f)
  expect_equal(z$table$n,f$n)
  expect_equal(z$table$probability,unname(f$yprob[,"1"]))
  expect_equal(z$table$node,as.integer(rownames(f)))
  for(i in which(!z$table$terminal)) {
    kids <- match(z$table$node[i]*2L+c(0L,1L),z$table$node)
    expect_equal(sum(z$table$n[kids]),z$table$n[i])
    expect_identical(z$table$incoming_rule[kids],paste(z$table$split_feature[i],f$splits[i,]))
  }
  native_x <- d[features]; names(native_x) <- sb_tree_binary_internal_names(3L)
  leaves <- predict(e,native_x,type="where")
  expect_equal(z$table$probability[leaves],sb_predict_final(fit,d,"patient_id")$predictions$probability)
  p <- sb_classification_tree_diagnostics(fit,"leaf_probabilities")
  expect_equal(p$table$probability,unname(f$yprob[f$var=="<leaf>","1"]))
  expect_equal(sum(p$table$n),nrow(d))
  expect_identical(sb_native_diagnostics(fit,complete=TRUE)$status$status,c("available","available"))
  for(view in list(z,p)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(view)),"ggplot_built")
  expect_identical(serialize(fit,NULL,version=2L),before)
})

# ---- algorithm-ctree ----

test_that("conditional trees retain native branches, membership and compact state", {
  skip_if_not_installed("partykit")
  for (id in c("ctree_binary", "ctree_survival")) {
    binary <- id == "ctree_binary"
    d <- if (binary) sb_example_binary(160L, seed = 108L) else sb_example_survival(160L, seed = 8612L)
    if (binary) d$gene_a <- d$gene_a + 1.5 * d$outcome
    features <- c("gene_c", "gene_a", "gene_b")
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else
      sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    e <- fit$diagnostic_artifacts$learner$engine
    view <- sb_ctree_structure(fit)
    tab <- view$table
    expect_equal(tab$node, partykit::nodeids(e))
    expect_equal(tab$node[tab$terminal], partykit::nodeids(e, terminal = TRUE))
    expect_gt(nrow(tab), 1L)
    membership <- partykit::fitted_node(partykit::node_party(e), partykit::data_party(e))
    expect_equal(tab$n[tab$terminal], as.integer(table(factor(membership, levels = tab$node[tab$terminal]))))
    expect_equal(tab$n[1], nrow(d))
    expect_equal(sum(tab$n[tab$terminal]), nrow(d))
    inspect <- function(node) {
      key <- as.character(partykit::id_node(node))
      split <- partykit::split_node(node)
      expect_identical(view$native_splits[[key]], split)
      expect_identical(view$native_node_information[[key]], partykit::info_node(node))
      kids <- partykit::kids_node(node)
      if (length(kids)) {
        labels <- partykit::character_split(split, partykit::data_party(e), digits = 7)
        children <- match(vapply(kids, partykit::id_node, 1L), tab$node)
        expect_equal(tab$parent[children], rep(as.integer(key), length(kids)))
        expect_identical(tab$incoming_rule[children], paste(labels$name, labels$levels))
        expect_equal(sum(tab$n[children]), tab$n[match(as.integer(key), tab$node)])
        invisible(lapply(kids, inspect))
      }
    }
    inspect(partykit::node_party(e))
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_identical(sb_ctree_structure(compact), view)
    expect_identical(sb_native_diagnostics(compact)$status$status, "available")
    p <- sb_plot_algorithm_diagnostic(view)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    expect_null(p$labels$subtitle)
    expect_null(p$labels$caption)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})

test_that("an unsplit conditional tree remains a valid one-node view", {
  skip_if_not_installed("partykit")
  d <- sb_example_binary(160L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_c", "gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "ctree_binary"))
  view <- sb_ctree_structure(fit)
  expect_equal(nrow(view$table), 1L)
  expect_identical(view$table$terminal, TRUE)
  expect_identical(view$table$incoming_rule, "Root")
  expect_equal(view$table$n, 160L)
  expect_null(view$native_splits[[1]])
  expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(view)), "ggplot_built")
})

# ---- algorithm-rpart-tree ----

test_that("single-tree diagnostics preserve all nodes, exact splits and native predictions", {
  skip_if_not_installed("rpart")
  for (kind in c("binary", "survival")) {
    d <- if (kind == "binary") make_sb_binary_fixture(n = 150L, seed = 108L) else
      make_sb_fixture(n = 150L, seed = 8612L)
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (kind == "binary") sb_binary_task(d, "patient_id", "outcome", features) else
      sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = if (kind == "binary") "decision_tree_binary" else "rpart_survival"),
      capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    original <- fit
    q <- sb_tree_structure(fit)
    engine <- fit$diagnostic_artifacts$learner$engine
    expect_equal(q$table$node, as.numeric(rownames(engine$frame)))
    expect_identical(q$native_splits, engine$splits)
    expect_equal(q$table$n, engine$frame$n)
    expect_equal(q$table$parent[-1], floor(q$table$node[-1] / 2))
    expected <- if (kind == "binary") engine$frame$yval2[, 5] else engine$frame$yval
    expect_equal(q$table$prediction, unname(expected))
    expect_equal(sum(q$table$n[q$table$terminal]), nrow(d))
    expect_identical(fit, original)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- sb_plot_algorithm_diagnostic(q)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_equal(nrow(p$data), nrow(engine$frame))
    }
  }
})

test_that("single-node trees and missing captures are explicit", {
  skip_if_not_installed("rpart")
  d <- make_sb_binary_fixture(n = 60L, seed = 108L)
  d$gene_a[] <- 0
  d$gene_b[] <- 0
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "decision_tree_binary"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  q <- sb_tree_structure(fit)
  expect_equal(nrow(q$table), 1)
  expect_identical(q$table$incoming_rule, "root")
  expect_equal(q$table$prediction, mean(d$outcome))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
  }
  fit$diagnostic_artifacts <- NULL
  err <- tryCatch(sb_tree_structure(fit), error = identity)
  expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
})
})
