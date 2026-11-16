# ---- multiclass-boosting ----
local({
test_that("boosting objectives produce portable native multiclass probabilities", {
  for (p in c("xgboost", "lightgbm", "catboost")) skip_if_not_installed(p)
  methods <- c("xgboost", "xgboost_dart", "xgboost_linear", "lightgbm",
    "lightgbm_dart", "lightgbm_goss", "lightgbm_rf", "catboost")
  for (k in c(3L, 4L)) {
    withr::local_seed(971)
    code <- rep(seq_len(k), each = 40)
    d <- data.frame(id = seq_along(code), type = paste0("Type ", code))
    for (j in seq_len(5)) d[[paste0("gene ", j)]] <- 3 * (code == j) + rnorm(length(code))
    task <- sb_multiclass_task(d, "id", "type", names(d)[-(1:2)],
      classes = rev(paste0("Type ", seq_len(k))))
    a <- sb_analyze_methods(task, methods)
    expect_identical(a$status$status, rep("success", length(methods)))
    for (method in methods) {
      fit <- a$results[[method]]$fit
      p <- sb_predict_multiclass(fit, d)
      expect_identical(colnames(p), task$classes)
      expect_equal(unname(rowSums(p)), rep(1, nrow(d)), tolerance = 1e-6)
      expect_gt(mean(task$classes[max.col(p)] == d$type), .8)
      rows <- c(nrow(d), 1, 41)
      expect_equal(unname(sb_predict_multiclass(fit, d[rows, ])), unname(p[rows, , drop = FALSE]))
      expect_equal(unname(sb_predict_multiclass(fit, d[1, , drop = FALSE])), unname(p[1, , drop = FALSE]))
      expect_equal(unname(sb_predict_multiclass(unserialize(serialize(fit, NULL)), d[rows, ])),
        unname(p[rows, , drop = FALSE]))
      x <- as.matrix(d[rows, task$predictors])
      colnames(x) <- fit$feature_codes
      native <- if (fit$package == "xgboost") {
        predict(xgboost::xgb.load.raw(fit$engine), x, reshape = TRUE)
      } else if (fit$package == "lightgbm") {
        predict(lightgbm::lgb.load(model_str = fit$engine), x)
      } else catboost::catboost.predict(unserialize(fit$engine), catboost::catboost.load_pool(x),
        prediction_type = "Probability", thread_count = 1L)
      expect_equal(unname(p[rows, , drop = FALSE]), unname(native), ignore_attr = TRUE)
      for (plot in a$results[[method]]$plots) expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
    }
  }
  comparison <- sb_compare_multiclass(task, methods, folds = 3)
  expect_identical(comparison$comparison$status, rep("success", length(methods)))
})
})

# ---- multiclass-classical ----
local({
# ---- multiclass-classical ----

classical_multiclass_task <- function(k = 3L) {
  withr::with_seed(816, {
    n <- 40L * k
    code <- rep(seq_len(k), each = 40L)
    d <- data.frame(id = seq_len(n), a = rnorm(n) + code,
      b = rnorm(n) - code / 2, c = rnorm(n), label = paste0("Type ", code))
    sb_multiclass_task(d, "id", "label", c("a", "b", "c"),
      classes = rev(paste0("Type ", seq_len(k))))
  })
}

test_that("four classical adapters retain native multiclass probabilities and evidence", {
  for (package in c("MASS", "e1071", "rpart", "ggplot2")) skip_if_not_installed(package)
  methods <- c("lda", "qda", "gaussian_naive_bayes", "decision_tree")
  for (k in c(3L, 4L)) {
    task <- classical_multiclass_task(k)
    analysis <- sb_analyze_methods(task, methods)
    expect_identical(analysis$status$status, rep("success", 4))
    for (method in methods) {
      result <- analysis$results[[method]]
      fit <- result$fit
      rows <- c(task$n, 1, 41)
      newdata <- task$data[rows, task$predictors]
      p <- sb_predict_multiclass(fit, newdata)
      native <- if (method %in% c("lda", "qda")) predict(fit$engine, newdata)$posterior else
        if (method == "gaussian_naive_bayes") predict(fit$engine, newdata, type = "raw") else {
          names(newdata) <- fit$feature_codes
          predict(fit$engine, newdata, type = "prob")
        }
      expect_equal(unname(p), unname(native[, fit$codes, drop = FALSE]))
      expect_identical(colnames(p), task$classes)
      expect_equal(unname(rowSums(p)), rep(1, length(rows)))
      expect_equal(unname(sb_predict_multiclass(fit, task$data[rows[1], , drop = FALSE])),
        unname(p[1, , drop = FALSE]))
      expect_s3_class(tryCatch(sb_selected_features(analysis, method), error = identity),
        "stormbreaker_selection_not_defined")
      if (method != "decision_tree") {
        for (label in task$classes) {
          tab <- sb_method_features(analysis, method, class = label)
          original <- task$data[task$data$label == label, task$predictors]
          expect_equal(unname(tab$mean), unname(colMeans(original)))
          if (method == "gaussian_naive_bayes") expect_equal(unname(tab$sd),
            unname(vapply(original, stats::sd, numeric(1))))
        }
      } else {
        tab <- sb_method_features(analysis, method)
        expect_identical(tab$used_for_split, fit$feature_codes %in% fit$engine$frame$var)
        expect_identical(fit$parameters$xval, 0L)
      }
      for (plot in result$plots) {
        expect_null(plot$labels$subtitle)
        expect_null(plot$labels$caption)
        expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
      }
      if (method == "decision_tree") {
        expect_match(result$notes$tree_structure, "Every original node", all = FALSE)
      }
    }
  }
  local_mocked_bindings(sb_fit_multiclass_classical = function(...) stop("Unexpected refit"))
  restored <- unserialize(serialize(analysis, NULL))
  expect_identical(sb_method_features(restored, "lda", "Type 2"),
    sb_method_features(analysis, "lda", "Type 2"))
  expect_equal(sb_predict_multiclass(restored$results$qda$fit, task$data),
    sb_predict_multiclass(analysis$results$qda$fit, task$data))
  expect_equal(sb_predict_multiclass(restored$results$decision_tree$fit, task$data),
    sb_predict_multiclass(analysis$results$decision_tree$fit, task$data))
  expect_identical(attr(restored$results$decision_tree$fit$engine$terms, ".Environment"), baseenv())
})

test_that("shared CV supports new methods and retains singular QDA failures", {
  for (package in c("MASS", "e1071", "rpart")) skip_if_not_installed(package)
  task <- classical_multiclass_task()
  methods <- c("lda", "qda", "gaussian_naive_bayes", "decision_tree")
  result <- sb_compare_multiclass(task, methods, folds = 3)
  expect_identical(result$comparison$status, rep("success", 4))
  for (method in methods) expect_equal(unname(rowSums(result$results[[method]]$probability)),
    rep(1, task$n))
  task$data$c <- task$data$a
  singular <- sb_compare_multiclass(task, c("qda", "gaussian_naive_bayes"), folds = 3)
  expect_identical(singular$comparison$status, c("failed", "success"))
  expect_null(singular$results$qda$metrics)
  expect_identical(is.na(singular$comparison$accuracy[1]), TRUE)
})

# ---- multiclass-distribution-plots ----

test_that("native distribution plots preserve densities and fitted means", {
  skip_if_not_installed("klaR")
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  for (method in c("kernel_naive_bayes", "rda")) {
    f <- sb_fit_multiclass(task, method, 3, 1)
    p <- sb_multiclass_native_plots(f)
    if (method == "kernel_naive_bayes") {
      expected <- unlist(lapply(f$engine$tables, function(t) unlist(lapply(t, `[[`, "y"))), use.names = FALSE)
      expect_equal(p$tables$class_densities$density, expected)
      expect_equal(nrow(p$tables$class_densities), 4L * 3L * 512L)
    } else {
      expect_equal(p$tables$class_means$value, as.vector(f$engine$means))
      co <- p$tables$class_correlations
      expect_equal(nrow(co), 4L * 4L * 3L)
      expect_equal(co$correlation[co$feature == co$other], rep(1, 12))
      expect_setequal(co$class, f$classes)
    }
    for (plot in p$plots) expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  }
})

# ---- multiclass-gp-plots ----

test_that("Gaussian-process figures preserve all native pairwise basis weights", {
  skip_if_not_installed("kernlab")
  d <- iris
  d$id <- seq_len(nrow(d))
  f <- sb_fit_multiclass(sb_multiclass_task(d, "id", "Species", names(iris)[1:4]), "gaussian_process", 3, 1)
  out <- sb_multiclass_native_plots(f)
  tab <- out$tables$kernel_weights
  expect_equal(tab$weight, as.numeric(unlist(f$engine@alpha)))
  expect_equal(tab$reference_row, unlist(f$engine@alphaindex, use.names = FALSE))
  pairs <- apply(utils::combn(f$classes, 2), 2, paste, collapse = " vs ")
  expect_equal(unique(tab$pair), pairs)
  expect_equal(nrow(tab), sum(lengths(f$engine@alpha)))
  expect_s3_class(ggplot2::ggplot_build(out$plots$kernel_weights), "ggplot_built")
})

# ---- multiclass-neighbor-plots ----

test_that("multiclass neighbor weights recover native predictions and paginate all queries", {
  skip_if_not_installed("kknn")
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  for (method in c("knn", "kernel_weighted_knn")) {
    f <- sb_fit_multiclass(task, method, 3, 1)
    expect_error(sb_plots(f), "Supply new_data")
    q <- d[1:9, ]
    out <- sb_plots(f, new_data = q)
    tab <- out$diagnostics$neighborhoods
    predicted <- sb_predict_multiclass(f, q)
    votes <- matrix(0, nrow(q), length(f$classes))
    for (i in seq_len(nrow(tab))) {
      cl <- match(tab$class[i], f$classes)
      votes[tab$query[i], cl] <- votes[tab$query[i], cl] + tab$weight[i]
    }
    expect_equal(unname(predicted), votes, tolerance = 1e-10)
    expect_equal(rowSums(votes), rep(1, nrow(q)))
    expect_equal(sum(vapply(out$pages, function(p) nrow(p$data), 1L)), nrow(tab))
    expect_equal(length(out$pages), 2L)
    for (p in out$pages) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  }
})

# ---- multiclass-neural-plots ----

test_that("multiclass network heatmaps preserve every native connection", {
  skip_if_not_installed("nnet")
  x <- scale(as.matrix(iris[1:4]))
  target <- stats::model.matrix(~ Species - 1, iris)
  m <- withr::with_seed(41, nnet::nnet(x, target, size = 3, softmax = TRUE,
    maxit = 500, decay = 0.1, trace = FALSE))
  fit <- list(method = "nnet", engine = m, features = names(iris)[1:4],
    classes = levels(iris$Species))
  result <- sb_multiclass_native_plots(fit)
  d <- result$tables$connection_weights
  expect_equal(d$weight, m$wts)
  expect_equal(nrow(d), length(m$wts))
  expect_setequal(as.character(d$destination), c(paste("Hidden", 1:3), fit$classes))
  expect_equal(sum(d$source == "Bias"), 6L)
  expect_setequal(d$layer, c("Input to hidden", "Hidden to output"))
  expect_s3_class(ggplot2::ggplot_build(result$plots$connection_weights), "ggplot_built")
})
})

# ---- multiclass-discriminants ----
local({
# ---- multiclass-discriminants ----

test_that("remaining classical engines preserve multiclass labels and saved predictions", {
  for (p in c("mda", "ipred", "sda", "klaR", "pamr", "kernlab", "kknn", "LiblineaR")) skip_if_not_installed(p)
  methods <- c("fda", "mda", "pda", "slda", "sda", "loclda", "pamr", "gaussian_process",
    "kernel_weighted_knn", "liblinear_l2_logistic", "liblinear_l1_logistic")
  for (k in c(3L, 4L)) {
    withr::local_seed(937)
    label <- rep(seq_len(k), each = 40)
    d <- data.frame(id = seq_along(label), outcome = paste0("Type ", label))
    for (j in seq_len(5)) d[[paste0("gene ", j)]] <- 2 * (label == j) + rnorm(length(label))
    task <- sb_multiclass_task(d, "id", "outcome", names(d)[-(1:2)],
      classes = rev(paste0("Type ", seq_len(k))))
    a <- sb_analyze_methods(task, methods)
    expect_identical(a$status$status, rep("success", length(methods)))
    for (method in methods) {
      fit <- a$results[[method]]$fit
      p <- sb_predict_multiclass(fit, d)
      expect_identical(colnames(p), task$classes)
      expect_equal(unname(rowSums(p)), rep(1, nrow(d)), tolerance = 1e-6)
      expect_gt(mean(task$classes[max.col(p)] == d$outcome), .7)
      rows <- c(nrow(d), 1, 41)
      expect_equal(unname(sb_predict_multiclass(fit, d[rows, ])), unname(p[rows, , drop = FALSE]))
      expect_equal(unname(sb_predict_multiclass(fit, d[rows[1], , drop = FALSE])),
        unname(p[rows[1], , drop = FALSE]))
      expect_equal(sb_predict_multiclass(unserialize(serialize(fit, NULL)), d[rows, ]),
        sb_predict_multiclass(fit, d[rows, ]))
    }
  }
  comparison <- sb_compare_multiclass(task, methods, folds = 3)
  expect_identical(comparison$comparison$status, rep("success", length(methods)))
})

# ---- multiclass-loclda-plots ----

test_that("local LDA retains radius ties and excluded class neighbors", {
  skip_if_not_installed("klaR")
  d <- iris
  d$id <- seq_len(nrow(d))
  f <- sb_fit_multiclass(sb_multiclass_task(d, "id", "Species", names(iris)[1:4]), "loclda", 3, 1)
  expect_error(sb_plots(f), "Supply new_data")
  out <- sb_plots(f, new_data = d[1:7, ])
  expect_equal(length(out$pages), 2L)
  tab <- out$diagnostics$neighborhoods
  for (q in 1:7) {
    distances <- sqrt(rowSums(sweep(f$engine$learn, 2, as.numeric(d[q, 1:4]), "-")^2))
    radius <- sort(distances)[f$engine$k]
    selected <- which(distances <= radius)
    actual <- tab[tab$query == q, ]
    expect_equal(actual$reference_row, selected)
    expect_equal(actual$weight, vapply(distances[selected] / radius, f$engine$weight.func, 1), ignore_attr = TRUE)
  }
  e <- list(learn = matrix(c(0, 1, -1, 2, 5), ncol = 1), k = 3,
    grouping = factor(c("a", "a", "b", "c", "c")), lev = c("a", "b", "c"), weight.func = function(x) 1 - x)
  ties <- sb_loclda_neighbor_table(e, matrix(0, 1, 1))
  expect_equal(ties$reference_row, 1:3)
  expect_equal(ties$used, c(TRUE, TRUE, FALSE))
  expect_equal(ties$effective_weight, c(1, 0, 0))
  for (p in out$pages) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

# ---- multiclass-rda-dispatch ----

test_that("multiclass RDA uses klaR even when another package registers predict.rda", {
  skip_if_not_installed("klaR")
  d <- iris
  d$id <- seq_len(nrow(d))
  f <- sb_fit_multiclass(sb_multiclass_task(d, "id", "Species", names(iris)[1:4]), "rda", 3, 1)
  expected <- get("predict.rda", asNamespace("klaR"))(f$engine, d[1:5, 1:4])$posterior
  previous <- utils::getS3method("predict", "rda")
  withr::defer(registerS3method("predict", "rda", previous, envir = asNamespace("stats")))
  registerS3method("predict", "rda", function(...) stop("Unrelated RDA method"), envir = asNamespace("stats"))
  expect_error(stats::predict(f$engine, d[1:5, 1:4]), "Unrelated RDA")
  expect_equal(unname(sb_predict_multiclass(f, d[1:5, ])), unname(expected[, f$codes]))
})
})

# ---- multiclass-linear ----
local({
# ---- multiclass-linear-weight-plots ----

test_that("linear multiclass weight plots retain native parameters", {
  skip_if_not_installed("LiblineaR")
  skip_if_not_installed("xgboost")
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  for (method in c("liblinear_l2_logistic", "liblinear_l1_logistic", "xgboost_linear")) {
    f <- sb_fit_multiclass(task, method, 3, 1)
    p <- sb_multiclass_native_plots(f)
    tab <- p$tables$linear_weights
    expect_equal(nrow(tab), (length(f$features) + 1L) * length(f$classes))
    expect_setequal(tab$class, f$classes)
    if (method != "xgboost_linear") expect_equal(tab$weight, as.vector(f$engine$W)) else {
      b <- sb_load_xgboost_state_model(list(model_raw = f$engine, params = f$parameters),
        "xgboost_binary_linear")$booster
      x <- as.matrix(d[1:4])
      colnames(x) <- f$feature_codes
      native <- predict(b, x, outputmargin = TRUE, reshape = TRUE)
      beta <- matrix(tab$weight, nrow = 3)
      reconstructed <- cbind(x, 1) %*% t(beta)
      expect_equal(sweep(reconstructed, 1, reconstructed[, 1]),
        sweep(native, 1, native[, 1]), tolerance = 1e-5, ignore_attr = TRUE)
    }
    expect_s3_class(ggplot2::ggplot_build(p$plots$linear_weights), "ggplot_built")
  }
})

# ---- multiclass-logistic-plots ----

test_that("multinomial plots retain coefficients and native AIC history", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("brglm2")
  d <- iris
  d$Species <- withr::with_seed(91, sample(d$Species))
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  for (method in c("logistic", "stepglm", "brglm2")) {
    fit <- sb_fit_multiclass(task, method, 3, 1)
    plots <- sb_multiclass_native_plots(fit)
    expect_equal(plots$tables$coefficients$coefficient, as.vector(stats::coef(fit$engine)))
    expect_length(grep("feature[0-9]", plots$tables$coefficients$feature), 0L)
    expect_setequal(unique(plots$tables$coefficients$contrast),
      paste(fit$classes[-1], "vs", fit$classes[1]))
    expect_s3_class(ggplot2::ggplot_build(plots$plots$coefficients), "ggplot_built")
    if (method == "stepglm") expect_equal(plots$tables$selection_path$AIC, fit$engine$anova$AIC)
  }
})

# ---- multiclass-multinomial ----

test_that("multinomial model families preserve class order and native saved predictions", {
  for (p in c("nnet", "brglm2", "mgcv", "pre")) skip_if_not_installed(p)
  methods <- c("logistic", "stepglm", "brglm2", "gam", "pre")
  for (k in c(3L, 4L)) {
    withr::local_seed(543)
    labels <- rep(seq_len(k), each = 35)
    d <- data.frame(id = seq_along(labels), label = paste0("Type ", labels))
    for (j in seq_len(3)) d[[paste0("gene ", j)]] <- 2 * (labels == j) + rnorm(length(labels))
    task <- sb_multiclass_task(d, "id", "label", names(d)[-(1:2)], classes = rev(paste0("Type ", seq_len(k))))
    a <- sb_analyze_methods(task, methods, inner_folds = 3)
    expect_identical(a$status$status, rep("success", length(methods)),
      info = paste(a$status$method, a$status$reason, collapse = "; "))
    for (method in methods) {
      fit <- a$results[[method]]$fit
      p <- sb_predict_multiclass(fit, d)
      expect_identical(colnames(p), task$classes)
      expect_equal(unname(rowSums(p)), rep(1, nrow(d)), tolerance = 1e-6)
      expect_gt(mean(task$classes[max.col(p)] == d$label), .65)
      rows <- c(nrow(d), 1, 36)
      expect_equal(unname(sb_predict_multiclass(fit, d[rows, ])), unname(p[rows, , drop = FALSE]))
      expect_equal(unname(sb_predict_multiclass(fit, d[1, , drop = FALSE])), unname(p[1, , drop = FALSE]))
      expect_equal(unname(sb_predict_multiclass(unserialize(serialize(fit, NULL)), d[rows, ])),
        unname(p[rows, , drop = FALSE]))
    }
  }
  cv <- sb_compare_multiclass(task, methods, folds = 3, inner_folds = 3)
  expect_identical(cv$comparison$status, rep("success", length(methods)),
    info = paste(cv$comparison$method, cv$comparison$reason, collapse = "; "))
})

test_that("unverified SLOPE is not advertised or fitted", {
  catalog <- sb_method_catalog("multiclass")
  expect_identical(catalog$supported[catalog$method == "slope"], FALSE)
  d <- datasets::iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(d)[1:4])
  a <- sb_analyze_methods(task, "slope")
  expect_identical(a$status$status, "failed")
  expect_identical(a$results$slope$error_class, "stormbreaker_unsupported_method_task")
})
})

# ---- multiclass-omics ----
local({
# ---- multiclass-omics ----

test_that("multiclass LPC retains native scores and selected genes", {
  skip_if_not_installed("lpc")
  for (k in c(3L, 4L)) {
    withr::local_seed(190)
    labels <- rep(seq_len(k), each = 20)
    d <- as.data.frame(matrix(rnorm(length(labels) * 60), ncol = 60))
    d[1:4] <- d[1:4] + labels
    genes <- names(d)
    d$id <- seq_along(labels)
    d$outcome <- paste0("Type ", labels)
    task <- sb_multiclass_task(d, "id", "outcome", genes,
      classes = rev(paste0("Type ", seq_len(k))))
    a <- sb_analyze_methods(task, "lpc")
    expect_identical(a$status$status, "success", info = a$status$reason)
    native <- NULL
    invisible(capture.output(native <- withr::with_seed(1,
      lpc::LPC(t(as.matrix(d[genes])), as.integer(factor(d$outcome, levels = task$classes)),
        type = "multiclass"))))
    tab <- sb_method_features(a, "lpc")
    expect_equal(tab$lpc_score, as.numeric(native$lpcscores))
    expect_equal(tab$base_score, as.numeric(native$tscores))
    ranked <- order(-abs(native$lpcscores), seq_along(genes))
    expected <- genes[head(ranked, sb_high_dimensional_nlogn_size(nrow(d), length(genes)))]
    expect_setequal(sb_selected_features(a, "lpc"), expected)
    expect_identical(sb_selected_features(unserialize(serialize(a, NULL)), "lpc"),
      sb_selected_features(a, "lpc"))
    expect_identical(a$results$lpc$plot_status$status, "available")
  }
})

# ---- multiclass-pamr-plots ----

test_that("PAMR displays native shrunken contrasts without dropping genes", {
  skip_if_not_installed("pamr")
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  fit <- sb_fit_multiclass(task, "pamr", 3, 1)
  result <- sb_multiclass_native_plots(fit)
  native <- pamr::pamr.predict(fit$engine, newx = NULL,
    threshold = fit$parameters$threshold, type = "centroid")
  expected <- (native - fit$engine$centroid.overall) / fit$engine$sd
  tab <- result$tables$shrunken_centroids
  expect_equal(tab$contrast, as.vector(expected))
  expect_equal(nrow(tab), 12L)
  expect_setequal(as.character(tab$feature), names(iris)[1:4])
  expect_setequal(tab$class, levels(iris$Species))
  expect_s3_class(ggplot2::ggplot_build(result$plots$shrunken_centroids), "ggplot_built")
  fit$parameters$threshold <- 1e6
  zero <- sb_multiclass_native_plots(fit)$tables$shrunken_centroids
  expect_equal(nrow(zero), 12L)
  expect_true(all(zero$contrast == 0))
})

# ---- multiclass-projection-plots ----

test_that("projection plots retain native coordinates and component identities", {
  for (package in c("mda", "ipred", "sda")) skip_if_not_installed(package)
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  for (method in c("fda", "mda", "pda", "slda", "sda")) {
    f <- sb_fit_multiclass(task, method, 3, 1)
    p <- sb_multiclass_native_plots(f)
    if (method %in% c("fda", "mda", "pda")) {
      expect_equal(p$tables$projected_centroids$value, as.vector(f$engine$means))
      expect_length(unique(p$tables$projected_centroids$label), nrow(f$engine$means))
    } else if (method == "slda") {
      expect_equal(p$tables$projection_loadings$value, as.vector(f$engine$scores))
      expect_equal(p$tables$projected_centroids$value, as.vector(f$engine$mylda$means))
    } else {
      expect_equal(p$tables$discriminant_coefficients$value, as.vector(t(f$engine$beta)))
      expect_equal(p$tables$shrinkage$intensity, unname(f$engine$regularization))
    }
    if (method == "mda") expect_equal(p$tables$subclass_weights$probability,
      unname(unlist(f$engine$sub.prior)))
    for (plot in p$plots) expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  }
})
})

# ---- multiclass-specialists ----
local({
# ---- multiclass-gam-plots ----

test_that("GAM smooths match native prediction matrices and retain training rugs", {
  skip_if_not_installed("mgcv")
  d <- iris
  d$id <- seq_len(nrow(d))
  set.seed(27)
  d$Species <- sample(d$Species)
  f <- sb_fit_multiclass(sb_multiclass_task(d, "id", "Species", names(iris)[1:2]), "gam", 3, 1)
  out <- sb_multiclass_native_plots(f)
  expect_equal(nrow(out$tables$training_rug), nrow(d) * 4L)
  expect_setequal(out$tables$smooths$class, levels(d$Species)[2:3])
  s <- f$engine$smooth[[1]]
  p <- s$first.para:s$last.para
  tab <- out$tables$smooths[1:100, ]
  nd <- f$engine$model[rep(1, 100), , drop = FALSE]
  nd[[s$term]] <- tab$value
  X <- stats::predict(f$engine, nd, type = "lpmatrix")[, p, drop = FALSE]
  expect_equal(tab$estimate, as.numeric(X %*% stats::coef(f$engine)[p]), ignore_attr = TRUE)
  expect_equal(tab$se, sqrt(diag(X %*% f$engine$Vp[p, p] %*% t(X))), ignore_attr = TRUE)
  pages <- sb_plot_pages(out$plots)
  for (plot in pages) expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
})

# ---- multiclass-pre-plots ----

test_that("PRE preserves native multiclass coefficients and coefficient-only CV", {
  skip_if_not_installed("pre")
  d <- iris
  d$id <- seq_len(nrow(d))
  f <- sb_fit_multiclass(sb_multiclass_task(d, "id", "Species", names(iris)[1:4]), "pre", 3, 1)
  out <- sb_multiclass_native_plots(f)
  native <- stats::coef(f$engine, penalty.par.val = "lambda.1se")
  expect_equal(nrow(out$tables$coefficients), nrow(native) * 3L)
  for (k in 1:3) expect_equal(out$tables$coefficients$coefficient[out$tables$coefficients$class == f$classes[k]],
    native[[paste0("coefficient.", f$codes[k])]])
  expect_equal(out$tables$coefficient_cv$mean_value, f$engine$glmnet.fit$cvm)
  expect_equal(out$tables$coefficient_cv$standard_error, f$engine$glmnet.fit$cvsd)
  expect_false(any(grepl("feature[0-9]", out$tables$coefficients$description)))
  for (p in sb_plot_pages(out$plots)) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

# ---- multiclass-specialist-plots ----

test_that("specialist figures preserve native means, spectra and loadings", {
  for (package in c("mclust", "HDclassif", "hda")) skip_if_not_installed(package)
  d <- withr::with_seed(701, data.frame(id = 1:90,
    x = rnorm(90), y = rnorm(90), z = rnorm(90), w = rnorm(90),
    label = rep(c("A", "B", "C"), 30)))
  task <- sb_multiclass_task(d, "id", "label", c("x", "y", "z", "w"))
  for (method in c("edda", "mclustda", "hdda", "hda")) {
    f <- sb_fit_multiclass(task, method, 3, 1)
    p <- sb_multiclass_native_plots(f)
    if (method %in% c("edda", "mclustda")) {
      expect_equal(p$tables$component_means$value,
        unlist(lapply(f$engine$models, function(x) as.vector(x$parameters$mean)), use.names = FALSE))
      expect_setequal(p$tables$component_means$class, f$classes)
    } else if (method == "hdda") {
      expect_equal(p$tables$eigenvalue_spectrum$value, as.vector(f$engine$ev))
      expect_equal(sum(p$tables$eigenvalue_spectrum$retained), sum(f$engine$d))
    } else {
      expect_equal(p$tables$projection_loadings$value, as.vector(f$engine$hda.loadings))
      expect_equal(sum(p$tables$projection_loadings$retained), length(f$features) * f$engine$reduced.dimension)
    }
    for (plot in p$plots) expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  }
})

# ---- multiclass-specialists ----

test_that("specialist engines return aligned portable multiclass predictions", {
  for (p in c("mclust", "HDclassif", "hda", "RLT", "aorsf", "wsrf")) skip_if_not_installed(p)
  methods <- c("edda", "mclustda", "hdda", "hda", "rlt", "aorsf", "wsrf")
  for (k in c(3L, 4L)) {
    withr::local_seed(934)
    label <- rep(seq_len(k), each = 45)
    d <- data.frame(id = seq_along(label), outcome = paste0("Type ", label))
    for (j in seq_len(6)) d[[paste0("gene ", j)]] <- 2 * (label == j) + rnorm(length(label))
    task <- sb_multiclass_task(d, "id", "outcome", names(d)[-(1:2)],
      classes = rev(paste0("Type ", seq_len(k))))
    a <- sb_analyze_methods(task, methods)
    expect_identical(a$status$status, rep("success", length(methods)))
    for (method in methods) {
      fit <- a$results[[method]]$fit
      p <- sb_predict_multiclass(fit, d)
      expect_identical(colnames(p), task$classes)
      expect_equal(unname(rowSums(p)), rep(1, nrow(d)), tolerance = 1e-6)
      if (method == "hda") {
        x <- as.matrix(d[task$predictors])
        colnames(x) <- fit$feature_codes
        native <- predict(fit$engine, newdata = x, task = "c")
        expect_equal(unname(p), unname(native$posteriors[, fit$codes, drop = FALSE]))
        expect_identical(fit$codes[max.col(p)], as.character(native$prediction))
      } else expect_gt(mean(task$classes[max.col(p)] == d$outcome), .65)
      rows <- c(nrow(d), 1, 46)
      expect_equal(unname(sb_predict_multiclass(fit, d[rows, ])), unname(p[rows, , drop = FALSE]))
      expect_equal(unname(sb_predict_multiclass(fit, d[1, , drop = FALSE])), unname(p[1, , drop = FALSE]))
      expect_equal(unname(sb_predict_multiclass(unserialize(serialize(fit, NULL)), d[rows, ])),
        unname(p[rows, , drop = FALSE]))
    }
  }
  cv <- sb_compare_multiclass(task, methods, folds = 3)
  expect_identical(cv$comparison$status, rep("success", length(methods)))
})

# ---- multiclass-svm-plots ----

test_that("multiclass SVM figures retain native pairwise decision values", {
  skip_if_not_installed("e1071")
  y <- factor(paste0("class", as.integer(iris$Species)))
  engine <- e1071::svm(iris[1:4], y, kernel = "radial", probability = TRUE)
  fit <- list(method = "svm_radial", engine = engine, codes = levels(y), classes = levels(iris$Species))
  drawn <- sb_multiclass_native_plots(fit)
  expect_named(drawn$plots, c("support_vectors", "decision_distribution"))
  expect_equal(drawn$tables$decision_distribution$decision, as.vector(engine$decision.values))
  expect_equal(nrow(drawn$tables$decision_distribution), 450L)
  expect_length(unique(drawn$tables$decision_distribution$pair), 3)
  expect_s3_class(ggplot2::ggplot_build(drawn$plots$decision_distribution), "ggplot_built")
})
})

# ---- multiclass-trees ----
local({
# ---- multiclass-forest-structure-plots ----

test_that("forest structure figures count every original tree", {
  for (package in c("partykit", "ipred", "RLT", "aorsf", "wsrf")) skip_if_not_installed(package)
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  for (method in c("cforest", "bagging", "rlt", "aorsf", "wsrf")) {
    f <- sb_fit_multiclass(task, method, 3, 1)
    p <- sb_multiclass_native_plots(f)
    trees <- switch(method, cforest = length(f$engine$nodes), bagging = length(f$engine$mtrees),
      rlt = length(f$engine$FittedTrees), aorsf = f$engine$n_tree, wsrf = length(f$engine$trees))
    expect_equal(nrow(p$tables$tree_size), trees)
    expect_equal(sum(p$plots$tree_size$data$trees), trees)
    axes <- ggplot2::ggplot_build(p$plots$tree_size)$layout$panel_params[[1]]
    for (axis in c("x", "y")) {
      ticks <- stats::na.omit(axes[[axis]]$get_breaks())
      expect_equal(as.numeric(ticks), floor(as.numeric(ticks)))
    }
    if (method %in% c("cforest", "bagging", "rlt")) {
      expect_equal(sum(p$tables$split_usage$splits), sum(p$tables$tree_size$nodes - 1L))
      expect_setequal(as.character(p$tables$split_usage$feature), f$features)
    }
    if (method == "wsrf") expect_equal(p$tables$individual_tree_oob$error,
      as.numeric(f$engine$treeOOBErrorRates))
    for (plot in p$plots) expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  }
})

# ---- multiclass-party-plots ----

test_that("multiclass conditional trees retain original nodes and memberships", {
  skip_if_not_installed("partykit")
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  f <- sb_fit_multiclass(task, "ctree", 3, 1)
  p <- sb_multiclass_native_plots(f)
  tab <- p$tables$tree_structure
  expect_setequal(tab$node, partykit::nodeids(f$engine))
  expect_equal(sum(tab$n[tab$terminal]), nrow(d))
  expect_equal(tab$n[is.na(tab$parent)], nrow(d))
  expect_length(grep("feature[0-9]", tab$incoming_rule), 0L)
  expect_s3_class(ggplot2::ggplot_build(p$plots$tree_structure), "ggplot_built")
})

# ---- multiclass-tree-native-plots ----

test_that("C50 and regularized forest plots use saved native quantities", {
  skip_if_not_installed("C50")
  skip_if_not_installed("RRF")
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  for (method in c("c50", "c50_rules", "rrf")) {
    fit <- sb_fit_multiclass(task, method, 3, 1)
    p <- sb_multiclass_native_plots(fit)
    native <- if (method == "rrf") fit$engine$importance else C50::C5imp(fit$engine, metric = "usage", pct = TRUE)
    expect_equal(p$tables$importance$importance, as.numeric(native[, 1]))
    expect_setequal(as.character(p$tables$importance$feature), fit$features[match(rownames(native), fit$feature_codes)])
    expect_s3_class(ggplot2::ggplot_build(p$plots$importance), "ggplot_built")
    if (method == "rrf") {
      expect_equal(p$tables$oob_error$error, as.vector(fit$engine$err.rate))
      expect_setequal(p$tables$oob_error$series, c("OOB", fit$classes))
      expect_s3_class(ggplot2::ggplot_build(p$plots$oob_error), "ggplot_built")
    }
  }
})

# ---- multiclass-tree-package-plots ----

test_that("tree package plotting retains the complete fitted tree", {
  skip_if_not_installed("tree")
  d <- iris
  d$id <- seq_len(nrow(d))
  fit <- sb_fit_multiclass(sb_multiclass_task(d, "id", "Species", names(iris)[1:4]), "tree", 3, 1)
  out <- sb_multiclass_native_plots(fit)
  tab <- out$tables$tree_structure
  frame <- fit$engine$frame
  expect_equal(tab$node, as.numeric(rownames(frame)))
  expect_equal(tab$n, frame$n)
  expect_equal(tab$prediction, fit$classes[match(as.character(frame$yval), fit$codes)])
  expect_equal(sum(tab$n[tab$terminal]), nrow(d))
  expect_true(all(tab$parent[-1] %in% tab$node))
  expect_false(any(grepl("feature[0-9]", tab$incoming_rule)))
  built <- ggplot2::ggplot_build(out$plots$tree_structure)
  expect_equal(nrow(built$data[[1]]), nrow(frame) - 1L)
  expect_equal(nrow(built$data[[2]]), nrow(frame))
})

# ---- multiclass-tree-structure ----

test_that("multiclass tree displays every original node and class", {
  skip_if_not_installed("rpart")
  d <- iris
  names(d)[1:4] <- paste0("feature", 1:4)
  d$Species <- factor(paste0("class", as.integer(d$Species)))
  engine <- rpart::rpart(Species ~ ., d, method = "class", control = rpart::rpart.control(xval = 0))
  fit <- list(engine = engine, feature_codes = paste0("feature", 1:4),
    features = names(iris)[1:4], codes = paste0("class", 1:3), classes = levels(iris$Species))
  tree <- sb_multiclass_tree_structure(fit)
  expect_equal(nrow(tree$table), nrow(engine$frame))
  expect_equal(tree$table$prediction, fit$classes[engine$frame$yval])
  expect_length(grep("feature[1-4]", tree$table$incoming_rule), 0)
  plot <- sb_plot_rpart_tree(tree)
  expect_equal(nrow(ggplot2::ggplot_build(plot)$data[[2]]), nrow(engine$frame))
})

# ---- multiclass-trees ----

test_that("native tree variants preserve multiclass probabilities", {
  for (p in c("ranger", "C50", "partykit", "ipred", "RRF", "tree", "evtree")) skip_if_not_installed(p)
  methods <- c("ranger", "ranger_extratrees", "c50", "c50_rules", "ctree",
    "cforest", "bagging", "rrf", "tree", "evtree")
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(d)[1:4],
    classes = c("virginica", "setosa", "versicolor"))
  a <- sb_analyze_methods(task, methods)
  expect_identical(a$status$status, rep("success", length(methods)))
  for (method in methods) {
    fit <- a$results[[method]]$fit
    rows <- c(120, 2, 61)
    p <- sb_predict_multiclass(fit, d[rows, ])
    expect_identical(colnames(p), task$classes)
    expect_equal(unname(rowSums(p)), rep(1, length(rows)), tolerance = 1e-6)
    expect_equal(unname(sb_predict_multiclass(fit, d[rows[1], , drop = FALSE])),
      unname(p[1, , drop = FALSE]))
    expect_equal(sb_predict_multiclass(unserialize(serialize(fit, NULL)), d[rows, ]), p)
    allp <- sb_predict_multiclass(fit, d)
    expect_gt(mean(task$classes[max.col(allp)] == d$Species), .8)
    expect_s3_class(tryCatch(sb_selected_features(a, method), error = identity),
      "stormbreaker_selection_not_defined")
  }
  cv <- sb_compare_multiclass(task, methods, folds = 3)
  expect_identical(cv$comparison$status, rep("success", length(methods)))
})
})

# ---- multiclass-workflow ----
local({
# ---- multiclass-information ----

test_that("multiclass information selectors reproduce native selected genes", {
  skip_if_not_installed("praznik")
  methods <- c("mrmr", "jmi", "cmim", "disr", "jmim")
  for (k in c(3L, 4L)) {
    withr::local_seed(122)
    y <- rep(seq_len(k), each = 15)
    d <- as.data.frame(matrix(rnorm(length(y) * 24), ncol = 24))
    d[[1]] <- y + rnorm(length(y), sd = .1)
    genes <- names(d)
    d$id <- seq_along(y)
    d$outcome <- paste0("Type ", y)
    task <- sb_multiclass_task(d, "id", "outcome", genes,
      classes = rev(paste0("Type ", seq_len(k))))
    a <- sb_analyze_methods(task, methods)
    expect_identical(a$status$status, rep("success", 5))
    for (method in methods) {
      f <- getExportedValue("praznik", toupper(method))
      args <- list(X = d[genes], Y = factor(d$outcome, levels = task$classes),
        k = sb_high_dimensional_nlogn_size(nrow(d), length(genes)), threads = 1L)
      if (method == "mrmr") args$positive <- FALSE
      native <- do.call(f, args)
      tab <- sb_method_features(a, method)
      expect_identical(tab$feature[order(tab$rank, na.last = NA)], genes[native$selection])
      expect_setequal(sb_selected_features(a, method), genes[native$selection])
      expect_identical(sb_selected_features(unserialize(serialize(a, NULL)), method),
        sb_selected_features(a, method))
      expect_identical(length(a$results[[method]]$plots), 1L)
    }
  }
})

# ---- multiclass-native ----

test_that("native probability engines support three and four classes", {
  for (p in c("e1071", "klaR", "nnet", "glmnet")) skip_if_not_installed(p)
  methods <- c("svm_linear", "svm_radial", "svm_radial_weighted", "svm_polynomial",
    "svm_sigmoid", "kernel_naive_bayes", "rda", "nnet", "knn", "ridge")
  for (k in c(3L, 4L)) {
    withr::local_seed(991)
    label <- rep(seq_len(k), each = 35L)
    d <- data.frame(id = seq_along(label), a = label + rnorm(length(label)),
      b = sin(label) + rnorm(length(label)), c = rnorm(length(label)),
      outcome = paste0("Type ", label))
    task <- sb_multiclass_task(d, "id", "outcome", c("a", "b", "c"),
      classes = rev(paste0("Type ", seq_len(k))))
    a <- sb_analyze_methods(task, methods, inner_folds = 3)
    expect_identical(a$status$status, rep("success", length(methods)))
    for (method in methods) {
      fit <- a$results[[method]]$fit
      rows <- c(nrow(d), 1, 36)
      p <- sb_predict_multiclass(fit, d[rows, ])
      expect_identical(colnames(p), task$classes)
      expect_equal(unname(rowSums(p)), rep(1, 3), tolerance = 1e-6)
      expect_equal(unname(sb_predict_multiclass(fit, d[rows[1], , drop = FALSE])),
        unname(p[1, , drop = FALSE]))
      expect_equal(sb_predict_multiclass(unserialize(serialize(fit, NULL)), d[rows, ]), p)
      expect_s3_class(tryCatch(sb_selected_features(a, method), error = identity),
        "stormbreaker_selection_not_defined")
      if (startsWith(method, "svm_")) {
        native <- attr(predict(fit$engine, d[rows, task$predictors], probability = TRUE), "probabilities")
        expect_equal(unname(p), unname(native[, fit$codes, drop = FALSE]))
        expect_equal(sum(a$results[[method]]$diagnostics$support_vectors$count), fit$engine$tot.nSV)
        counts <- table(factor(d$outcome[fit$engine$index], levels = task$classes))
        tab <- a$results[[method]]$diagnostics$support_vectors
        expect_equal(unname(tab$count[match(task$classes, as.character(tab$class))]), unname(as.integer(counts)))
      } else if (method == "rda") {
        native <- get("predict.rda", asNamespace("klaR"))(fit$engine, d[rows, task$predictors])$posterior
        expect_equal(unname(p), unname(native[, fit$codes, drop = FALSE]))
      } else if (method == "kernel_naive_bayes") {
        native <- predict(fit$engine, d[rows, task$predictors])$posterior
        expect_equal(unname(p), unname(native[, fit$codes, drop = FALSE]))
      } else if (method == "knn") {
        z <- sweep(sweep(as.matrix(d[rows, task$predictors]), 2, fit$center, "-"), 2, fit$scale, "/")
        for (i in seq_along(rows)) {
          distance <- rowSums((fit$engine$training - matrix(z[i, ], nrow(d), 3, byrow = TRUE))^2)
          neighbors <- order(distance, seq_along(distance))[1:5]
          expected <- as.numeric(table(factor(d$outcome[neighbors], levels = task$classes))) / 5
          expect_equal(unname(p[i, ]), expected)
        }
      } else if (method == "nnet") {
        z <- sweep(sweep(as.matrix(d[rows, task$predictors]), 2, fit$center, "-"), 2, fit$scale, "/")
        expect_equal(unname(p), unname(predict(fit$engine, z, type = "raw")))
        expect_equal(fit$center, vapply(d[task$predictors], mean, numeric(1)))
      }
    }
  }
  cv <- sb_compare_multiclass(task, c("svm_radial", "nnet", "knn", "ridge"),
    folds = 3, inner_folds = 3)
  expect_identical(cv$comparison$status, rep("success", 4))
})

# ---- multiclass-selection ----

multiclass_selection_task <- function(k = 3L) {
  withr::with_seed(1717, {
    label <- rep(seq_len(k), each = 30L)
    d <- data.frame(id = seq_along(label), group = paste0("Type ", label))
    for (j in seq_len(k)) d[[paste0("gene ", j)]] <- 3 * (label == j) + rnorm(length(label))
    d$noise <- rnorm(length(label))
    sb_multiclass_task(d, "id", "group", names(d)[-(1:2)],
      classes = rev(paste0("Type ", seq_len(k))))
  })
}

test_that("multiclass forest selectors return native selections and saved plots", {
  for (p in c("Boruta", "varSelRF", "VSURF")) skip_if_not_installed(p)
  task <- multiclass_selection_task()
  methods <- c("boruta", "varselrf", "vsurf")
  analysis <- sb_analyze_methods(task, methods)
  expect_identical(analysis$status$status, rep("success", 3))
  for (method in methods) {
    result <- analysis$results[[method]]
    expect_setequal(sb_selected_features(analysis, method), result$fit$selected_features)
    expect_identical(result$fit$classes, task$classes)
    expect_identical(result$features$selected,
      task$predictors %in% result$fit$selected_features)
    expect_identical(result$plot_status$status, rep("available", nrow(result$plot_status)))
    for (p in result$plots) {
      expect_null(p$labels$caption)
      expect_null(p$labels$subtitle)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    }
  }
  restored <- unserialize(serialize(analysis, NULL))
  expect_identical(sb_selected_features(restored, "boruta"), sb_selected_features(analysis, "boruta"))
  expect_s3_class(tryCatch(sb_compare_multiclass(task, "boruta", folds = 3), error = identity),
    "stormbreaker_bad_method")
})

test_that("ABESS uses native multinomial class order and GIC selection", {
  skip_if_not_installed("abess")
  for (k in c(3L, 4L)) {
    task <- multiclass_selection_task(k)
    a <- sb_analyze_methods(task, "abess")
    expect_identical(a$status$status, "success")
    fit <- a$results$abess$fit
    rows <- c(task$n, 1, 31)
    p <- sb_predict_multiclass(fit, task$data[rows, ])
    native <- predict(fit$engine, as.matrix(task$data[rows, task$predictors]),
      support.size = fit$engine$best.size, type = "response")[[1]]
    expect_equal(unname(p), unname(native))
    expect_identical(colnames(p), task$classes)
    expect_equal(unname(rowSums(p)), rep(1, length(rows)))
    expect_equal(unname(sb_predict_multiclass(fit, task$data[rows[1], , drop = FALSE])),
      unname(p[1, , drop = FALSE]))
    allp <- sb_predict_multiclass(fit, task$data)
    expect_gt(mean(task$classes[max.col(allp)] == task$data$group), .85)
    beta <- as.matrix(coef(fit$engine, support.size = fit$engine$best.size)[[1]])[-1, , drop = FALSE]
    expect_setequal(sb_selected_features(a, "abess"), task$predictors[rowSums(beta != 0) > 0])
    reference <- sb_method_features(a, "abess", class = tail(task$classes, 1))
    expect_equal(reference$coefficient, rep(0, length(task$predictors)))
    expect_identical(a$results$abess$plot_status$status, c("available", "available"))
    for (plot in a$results$abess$plots) expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
    restored <- unserialize(serialize(fit, NULL))
    expect_equal(sb_predict_multiclass(restored, task$data), allp)
  }
  cv <- sb_compare_multiclass(task, "abess", folds = 3)
  expect_identical(cv$comparison$status, "success")
  expect_equal(unname(rowSums(cv$results$abess$probability)), rep(1, task$n))
})

# ---- multiclass ----

multiclass_test_task <- function(grouped = FALSE) {
  d <- datasets::iris
  d$id <- seq_len(nrow(d))
  if (grouped) {
    d$patient <- d$id
    d <- d[rep(seq_len(nrow(d)), each = 2), ]
    d$id <- seq_len(nrow(d))
  }
  sb_multiclass_task(d, "id", "Species", names(d)[1:4],
    classes = c("virginica", "setosa", "versicolor"),
    group_col = if (grouped) "patient" else NULL)
}

test_that("multiclass task preserves labels and rejects ambiguous roles", {
  task <- multiclass_test_task()
  expect_identical(task$classes, c("virginica", "setosa", "versicolor"))
  expect_identical(task$data$Species, datasets::iris$Species)
  bad <- task$data
  bad$Species[1] <- NA
  expect_s3_class(tryCatch(sb_multiclass_task(bad, "id", "Species", task$predictors),
    error = identity), "stormbreaker_bad_multiclass_outcome")
  expect_s3_class(tryCatch(sb_multiclass_task(task$data, "id", "Species", "Species"),
    error = identity), "stormbreaker_bad_predictors")
  bad <- task$data[task$data$Species != "setosa", ]
  expect_s3_class(tryCatch(sb_multiclass_task(bad, "id", "Species", task$predictors),
    error = identity), "stormbreaker_bad_multiclass_outcome")
  expect_identical(sb_method_catalog("multiclass")$method[
    sb_method_catalog("multiclass")$supported], c("boruta", "varselrf", "vsurf",
      "lasso", "elastic_net", "abess", "random_forest",
      "lda", "qda", "gaussian_naive_bayes", "decision_tree", "ridge",
      "svm_linear", "svm_radial", "svm_radial_weighted", "svm_polynomial", "svm_sigmoid",
      "kernel_naive_bayes", "rda", "nnet", "knn", "ranger", "ranger_extratrees",
      "c50", "c50_rules", "ctree", "cforest", "bagging", "rrf", "tree", "evtree",
      "fda", "mda", "pda", "slda", "sda", "loclda", "pamr", "gaussian_process",
      "kernel_weighted_knn", "liblinear_l2_logistic", "liblinear_l1_logistic",
      "xgboost", "xgboost_dart", "xgboost_linear", "lightgbm", "lightgbm_dart",
      "lightgbm_goss", "lightgbm_rf", "catboost",
      "edda", "mclustda", "hdda", "hda", "rlt", "aorsf", "wsrf",
      "logistic", "stepglm", "brglm2", "gam", "pre",
      "mrmr", "jmi", "cmim", "disr", "jmim", "lpc"))
})

test_that("native three-class fits preserve probabilities and selection per class", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("ggplot2")
  task <- multiclass_test_task()
  analysis <- sb_analyze_methods(task, c("lasso", "elastic_net", "random_forest", "mcp"), inner_folds = 3)
  expect_identical(analysis$status$status, c(rep("success", 3), "failed"))
  expect_identical(analysis$results$mcp$error_class, "stormbreaker_unsupported_method_task")
  expect_identical(analysis$task_type, "multiclass")
  rows <- c(121, 2, 61)
  for (method in c("lasso", "elastic_net", "random_forest")) {
    result <- analysis$results[[method]]
    fit <- result$fit
    p <- sb_predict_multiclass(fit, task$data[rows, ])
    expect_identical(colnames(p), task$classes)
    expect_equal(unname(rowSums(p)), rep(1, 3), tolerance = 1e-8)
    one <- sb_predict_multiclass(fit, task$data[rows[1], , drop = FALSE])
    expect_equal(unname(one), unname(p[1, , drop = FALSE]))
    if (method == "random_forest") {
      native <- predict(fit$engine, task$data[rows, task$predictors], type = "prob")
      expect_equal(unname(p), unname(as.matrix(native[, fit$codes])))
      expect_identical(sb_method_features(analysis, method)$selected, rep(NA, 4))
      expect_s3_class(tryCatch(sb_selected_features(analysis, method), error = identity),
        "stormbreaker_selection_not_defined")
    } else {
      native <- predict(fit$engine, as.matrix(task$data[rows, task$predictors]),
        s = fit$lambda, type = "response")
      expect_equal(unname(p), unname(native[, fit$codes, 1]))
      beta <- coef(fit$engine, s = fit$lambda)
      union <- character()
      for (k in seq_along(task$classes)) {
        b <- as.matrix(beta[[fit$codes[k]]])
        selected <- rownames(b)[b[, 1] != 0 & rownames(b) != "(Intercept)"]
        expect_setequal(sb_selected_features(analysis, method, class = task$classes[k]), selected)
        union <- union(union, selected)
      }
      expect_setequal(sb_selected_features(analysis, method), union)
    }
    expect_identical(result$plot_status$status, rep("available", length(result$plots)))
    class_plot <- if (method == "random_forest") result$plots$class_importance else
      result$plots$coefficient_path
    layout <- ggplot2::ggplot_build(class_plot)$layout$layout
    expect_identical(as.character(layout$class), task$classes)
    expect_identical(any(grepl("confusion", names(result$plots))), FALSE)
    for (plot in result$plots) {
      expect_null(plot$labels$subtitle)
      expect_null(plot$labels$caption)
      expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
    }
  }
  local_mocked_bindings(sb_fit_multiclass = function(...) stop("Unexpected refit"))
  saved <- unserialize(serialize(analysis, NULL))
  expect_identical(sb_selected_features(saved, "lasso", "setosa"),
    sb_selected_features(analysis, "lasso", "setosa"))
  expect_equal(sb_predict_multiclass(saved$results$lasso$fit, task$data[rows, ]),
    sb_predict_multiclass(analysis$results$lasso$fit, task$data[rows, ]))
})

test_that("metrics match hand calculations and confusion counts do not overlap", {
  truth <- c("a", "b", "c", "a", "b", "c")
  p <- matrix(c(.8,.1,.1, .1,.8,.1, .1,.1,.8,
    .1,.8,.1, .1,.8,.1, .1,.8,.1), byrow = TRUE, ncol = 3,
    dimnames = list(NULL, c("a", "b", "c")))
  m <- sb_multiclass_metrics(truth, p)
  expect_equal(m$overall$accuracy, 4 / 6)
  expect_equal(m$overall$balanced_accuracy, 2 / 3)
  expect_equal(m$overall$macro_f1, 2 / 3)
  expect_equal(m$overall$log_loss, -(4 * log(.8) + 2 * log(.1)) / 6)
  expect_equal(m$overall$brier, (4 * .06 + 2 * 1.46) / 6)
  absent <- sb_multiclass_metrics(truth[1:2], p[1:2, ])
  expect_identical(is.na(absent$overall$macro_auc), TRUE)
  expect_identical(is.na(absent$overall$balanced_accuracy), TRUE)
  invalid <- p
  invalid[1, 1] <- .9
  expect_s3_class(tryCatch(sb_multiclass_metrics(truth, invalid), error = identity),
    "stormbreaker_bad_multiclass_probabilities")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    text <- ggplot2::ggplot_build(sb_plot_multiclass_confusion(m))$data[[2]]
    expect_identical(nrow(text), 9L)
    expect_identical(anyDuplicated(text[c("x", "y")]), 0L)
    expect_equal(sum(as.numeric(text$label)), length(truth))
  }
})

test_that("four classes with engine-like labels keep their original mapping", {
  skip_if_not_installed("randomForest")
  skip_if_not_installed("glmnet")
  d <- withr::with_seed(15, {
    x <- data.frame(id = seq_len(120), x1 = rnorm(120), x2 = rnorm(120), x3 = rnorm(120))
    x$label <- rep(c("Overall", "MeanDecreaseGini", "A type", "class2"), each = 30)
    x
  })
  classes <- c("class2", "Overall", "A type", "MeanDecreaseGini")
  task <- sb_multiclass_task(d, "id", "label", c("x1", "x2", "x3"), classes)
  result <- sb_analyze_methods(task, c("random_forest", "lasso"), inner_folds = 3)
  expect_identical(result$status$status, c("success", "success"))
  for (fit in lapply(result$results, `[[`, "fit")) {
    p <- sb_predict_multiclass(fit, d[c(91, 1, 61, 31), ])
    expect_identical(colnames(p), classes)
    expect_equal(unname(rowSums(p)), rep(1, 4))
  }
})

test_that("nested folds keep repeated patients together and exclude held-out rows", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("randomForest")
  task <- multiclass_test_task(grouped = TRUE)
  v <- sb_compare_multiclass(task, c("lasso", "random_forest"), folds = 3, inner_folds = 3)
  expect_identical(v$comparison$status, c("success", "success"))
  groups <- task$data$patient
  expect_identical(all(vapply(split(v$assignment, groups), function(x) length(unique(x)) == 1, logical(1))), TRUE)
  for (r in v$results) {
    expect_identical(rownames(r$probability), as.character(task$data$id))
    expect_equal(unname(rowSums(r$probability)), rep(1, task$n))
    for (fold in r$folds) {
      expect_length(intersect(groups[fold$training_rows], groups[fold$assessment_rows]), 0)
      if (!is.null(fold$inner_fold_id)) {
        expect_length(fold$inner_fold_id, length(fold$training_rows))
        sizes <- vapply(split(fold$inner_fold_id, groups[fold$training_rows]),
          function(x) length(unique(x)), integer(1))
        expect_identical(all(sizes == 1L), TRUE)
      }
    }
  }
  failed <- sb_compare_multiclass(multiclass_test_task(), c("lasso", "random_forest"),
    folds = 3, inner_folds = 100)
  expect_identical(failed$comparison$status, c("failed", "success"))
  expect_null(failed$results$lasso$metrics)
  expect_identical(is.na(failed$comparison$accuracy[1]), TRUE)
})
})
