# ---- model-classical-binary ----
local({
# ---- classic-tree-rule-adapters ----

test_that("classic tree and C5.0 rules declare distinct native specifications", {
  expected_dependencies <- c(
    tree_binary = "tree",
    c50_rules_binary = "C50"
  )

  expect_setequal(
    intersect(
      names(expected_dependencies),
      sb_classical_binary_learner_ids()
    ),
    names(expected_dependencies)
  )
  for (learner in names(expected_dependencies)) {
    dependency <- sb_classical_binary_dependency(learner)
    expect_identical(
      dependency$package,
      unname(expected_dependencies[[learner]])
    )
    expect_identical(
      dependency$available,
      requireNamespace(dependency$package, quietly = TRUE)
    )
  }

  expect_identical(
    sb_classical_binary_params("tree_binary"),
    list(
      split = "deviance",
      mincut = 5L,
      minsize = 10L,
      mindev = 0.01,
      keep_model = FALSE,
      keep_x = FALSE,
      keep_y = FALSE,
      keep_weights = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("c50_rules_binary"),
    list(
      trials = 10L,
      rules = TRUE,
      winnow = FALSE,
      early_stopping = FALSE,
      seed = 1L
    )
  )

  selectors <- expected_fixed_learner_selectors("binary")
  registry <- sb_pipeline_registry()
  for (learner in names(expected_dependencies)) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(nrow(rows), length(selectors))
    expect_setequal(rows$selector, selectors)
    expect_identical(unique(rows$task_type), "binary")
    expect_identical(unique(rows$tuning_semantics), "none")
  }
})

test_that("classic tree and C5.0 rules reproduce native probabilities", {
  skip_if_not_installed("tree", minimum_version = "1.0-40")
  skip_if_not_installed("C50", minimum_version = "0.2.0")

  data <- make_sb_binary_fixture(n = 96L, seed = 8631L)
  names(data)[names(data) == "gene_a"] <- "gene A"
  names(data)[names(data) == "gene_b"] <- "gene-B"
  predictors <- setdiff(names(data), c("patient_id", "outcome"))
  x <- data[predictors]
  grouping <- factor(data$outcome, levels = c(0L, 1L))
  assessment <- x[1:9, , drop = FALSE]

  tree_adapter <- sb_classical_binary_learner_adapter(
    "tree_binary"
  )
  tree_model <- tree_adapter$fit(x, data$outcome)
  formula_data <- sb_tree_binary_formula_data(x, data$outcome)
  native_tree <- tree::tree(
    formula_data$formula,
    data = formula_data$data,
    control = tree::tree.control(
      nobs = nrow(x),
      mincut = 5L,
      minsize = 10L,
      mindev = 0.01
    ),
    split = "deviance",
    model = FALSE,
    x = FALSE,
    y = FALSE,
    wts = FALSE
  )
  tree_assessment <- assessment
  names(tree_assessment) <- sb_tree_binary_internal_names(ncol(assessment))
  expected_tree <- stats::predict(
    native_tree,
    newdata = tree_assessment,
    type = "vector"
  )[, "1"]

  expect_s3_class(tree_model$engine, "tree")
  expect_equal(
    tree_adapter$predict(tree_model, assessment),
    unname(expected_tree),
    tolerance = 0
  )

  rules_adapter <- sb_classical_binary_learner_adapter(
    "c50_rules_binary"
  )
  rules_model <- rules_adapter$fit(x, data$outcome)
  native_rules <- C50::C5.0(
    x = x,
    y = grouping,
    trials = 10L,
    rules = TRUE,
    control = C50::C5.0Control(
      winnow = FALSE,
      seed = 1L,
      earlyStopping = FALSE
    )
  )
  expected_rules <- get(
    "predict.C5.0",
    envir = asNamespace("C50"),
    inherits = FALSE
  )(native_rules, newdata = assessment, type = "prob")[, "1"]

  expect_s3_class(rules_model$engine, "C5.0")
  expect_type(rules_model$engine$rules, "character")
  expect_true(nzchar(rules_model$engine$rules))
  expect_equal(
    rules_adapter$predict(rules_model, assessment),
    unname(expected_rules),
    tolerance = 0
  )
})

test_that("classic tree and C5.0 rule states replay without fitted rows", {
  skip_if_not_installed("tree", minimum_version = "1.0-40")
  skip_if_not_installed("C50", minimum_version = "0.2.0")

  data <- make_sb_binary_fixture(n = 84L, seed = 8632L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  x <- data[predictors]
  assessment <- x[1:7, , drop = FALSE]

  for (learner in c("tree_binary", "c50_rules_binary")) {
    adapter <- sb_classical_binary_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    expected <- adapter$predict(model, assessment)
    state <- adapter$extract_state(model)
    engine <- unserialize(state$engine_raw)
    replay <- adapter$predict_state(
      unserialize(serialize(state, NULL)),
      assessment
    )$probability

    expect_equal(replay, expected, tolerance = 0)
    expect_null(state$training_data_retention)
    if (identical(learner, "tree_binary")) {
      expect_identical(names(engine), c("frame", "terms"))
      expect_identical(environment(engine$terms), baseenv())
    } else {
      expect_identical(
        intersect(c("call", "output", "boostResults", "size"), names(engine)),
        character()
      )
    }
  }
})

test_that("classic tree and C5.0 rules execute through CV and final prediction", {
  skip_if_not_installed("tree", minimum_version = "1.0-40")
  skip_if_not_installed("C50", minimum_version = "0.2.0")

  data <- make_sb_binary_fixture(n = 72L, seed = 8633L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    predictors = c("gene_a", "gene_b", "gene_c")
  )

  for (learner in c("tree_binary", "c50_rules_binary")) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 2L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 8634L)
    final <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      final,
      data[1:8, c("patient_id", "gene_a", "gene_b", "gene_c")],
      id_col = "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(final$status, "success")
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$probability, 8L)
  }
})

# ---- classical-binary-boosting-fixed-adapters ----

test_that("classical binary boosting settings are fixed and explicit", {
  adaboost <- sb_pipeline(learner = "adaboost_binary")
  logitboost <- sb_pipeline(learner = "logitboost_binary")

  expect_identical(
    adaboost$id,
    "none__adaboost_binary__iterations50"
  )
  expect_identical(
    logitboost$id,
    "none__logitboost_binary__iterations51"
  )
  expect_identical(adaboost$learner_params$method, "AdaBoost.M1")
  expect_identical(adaboost$learner_params$rpart_cp, 0)
  expect_identical(logitboost$learner_params$method, "LogitBoost")
  expect_identical(logitboost$learner_params$base_learner, "decision_stump")
  expect_identical(
    adaboost$learner_params$iteration_semantics,
    "maximum_until_native_convergence"
  )
  expect_identical(logitboost$learner_params$iteration_semantics, "exact")
  expect_identical(adaboost$learner_params$tuning, "none")
  expect_identical(logitboost$learner_params$tuning, "none")

  expect_error(
    sb_pipeline(
      learner = "adaboost_binary",
      adaboost_iterations = 0L
    ),
    class = "stormbreaker_bad_integer"
  )
  expect_error(
    sb_pipeline(learner = "logistic", logitboost_iterations = 11L),
    class = "stormbreaker_irrelevant_classical_binary_boosting_parameter"
  )
  expect_error(
    sb_pipeline(
      learner = "adaboost_binary",
      logitboost_iterations = 11L
    ),
    class = "stormbreaker_irrelevant_classical_binary_boosting_parameter"
  )
})

test_that("classical binary boosting preserves native class-one probability", {
  skip_if_not_installed("fastAdaboost", minimum_version = "1.0.0")
  skip_if_not_installed("caTools", minimum_version = "1.18.3")

  d <- make_sb_binary_fixture(n = 54L, seed = 7911L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(4L, 19L, 37L), , drop = FALSE]
  adapters <- list(
    adaboost_binary =
      sb_classical_binary_boosting_learner_adapter(
        "adaboost_binary",
        9L
      ),
    logitboost_binary =
      sb_classical_binary_boosting_learner_adapter(
        "logitboost_binary",
        11L
      )
  )

  for (learner in names(adapters)) {
    adapter <- adapters[[learner]]
    model <- adapter$fit(x, d$outcome)
    observed <- adapter$predict(model, new_x)
    expected <- if (identical(learner, "adaboost_binary")) {
      native <- stats::predict(model$engine, newdata = new_x)$prob
      unname(native[, match("1", unname(model$engine$classnames))])
    } else {
      native <- stats::predict(
        model$engine,
        as.matrix(new_x),
        type = "raw",
        nIter = model$iterations
      )
      unname(native[, match("1", as.character(model$engine$lablist))])
    }
    audit <- adapter$extract_tree_model(model)
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    restored <- adapter$predict_state(state, new_x)$probability

    expect_equal(observed, expected, tolerance = 0)
    expect_equal(restored, expected, tolerance = 0)
    expect_identical(state$learner, learner)
    expect_false("engine" %in% names(state))
    expect_type(state$model_raw, "raw")
    expect_gt(length(state$model_raw), 0L)
    expect_identical(audit$iterations, adapter$params$iterations)
    expect_identical(audit$fitted_iterations, model$fitted_iterations)
    expect_lte(audit$fitted_iterations, audit$iterations)
    expect_identical(audit$importance_type, "split_count")
    expect_identical(audit$feature_importance$feature, names(x))
    expect_true(all(audit$feature_importance$importance >= 0))
  }
})

test_that("AdaBoost records native convergence before its iteration maximum", {
  skip_if_not_installed("fastAdaboost", minimum_version = "1.0.0")

  outcome <- rep(c(0L, 1L), each = 10L)
  x <- data.frame(
    separating_feature = ifelse(outcome == 1L, 1, -1),
    noise = seq_along(outcome)
  )
  adapter <- sb_classical_binary_boosting_learner_adapter(
    "adaboost_binary",
    5L
  )
  model <- adapter$fit(x, outcome)
  audit <- adapter$extract_tree_model(model)

  expect_identical(model$iterations, 5L)
  expect_identical(model$fitted_iterations, 1L)
  expect_identical(audit$fitted_iterations, 1L)
  expect_true(all(adapter$predict(model, x) %in% c(0, 1)))
})

test_that("classical binary boosting runs CV and portable final prediction", {
  skip_if_not_installed("fastAdaboost", minimum_version = "1.0.0")
  skip_if_not_installed("caTools", minimum_version = "1.18.3")

  d <- make_sb_binary_fixture(n = 72L, seed = 7921L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipelines <- list(
    sb_pipeline(
      learner = "adaboost_binary",
      adaboost_iterations = 9L
    ),
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "logitboost_binary",
      selector_n_features = 2L,
      logitboost_iterations = 11L
    )
  )
  version_fields <- c("fastadaboost_version", "catools_version")

  for (index in seq_along(pipelines)) {
    pipeline <- pipelines[[index]]
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 7922L)
    fit <- sb_fit_final(task, pipeline)
    new_data <- d[
      c(5L, 29L, 61L),
      c("patient_id", rev(predictors)),
      drop = FALSE
    ]
    prediction <- sb_predict_final(fit, new_data, "patient_id")
    direct <- sb_predict_classical_binary_boosting_learner_state(
      fit$learner_state,
      new_data[fit$required_features]
    )

    expect_identical(result$status, "success")
    expect_identical(result$complete_oof, TRUE)
    expect_length(result$tree_models, 3L)
    expect_identical(
      vapply(result$tree_models, `[[`, "", "status"),
      rep("success", 3L)
    )
    expect_identical(
      vapply(result$coefficients, `[[`, "", "status"),
      rep("not_applicable", 3L)
    )
    expect_true(version_fields[[index]] %in% names(result$provenance))
    expect_identical(fit$status, "success")
    expect_true(version_fields[[index]] %in% names(fit$provenance))
    expect_identical(prediction$status, "success")
    expect_equal(
      prediction$predictions$probability,
      direct$probability,
      tolerance = 0
    )
  }
})

test_that("classical binary boosting has no dependency fallback", {
  d <- make_sb_binary_fixture(n = 30L, seed = 7931L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adaboost <- sb_classical_binary_boosting_learner_adapter(
    "adaboost_binary",
    3L
  )
  logitboost <- sb_classical_binary_boosting_learner_adapter(
    "logitboost_binary",
    3L
  )
  local_mocked_bindings(
    sb_fastadaboost_available = function() FALSE,
    sb_catools_available = function() FALSE
  )

  expect_error(
    adaboost$fit(x, d$outcome),
    class = "stormbreaker_missing_fastadaboost"
  )
  expect_error(
    logitboost$fit(x, d$outcome),
    class = "stormbreaker_missing_catools"
  )
})

test_that("classical binary boosting rejects invalid fits despite warnings", {
  skip_if_not_installed("fastAdaboost", minimum_version = "1.0.0")
  skip_if_not_installed("caTools", minimum_version = "1.18.3")

  d <- make_sb_binary_fixture(n = 30L, seed = 7941L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adaboost <- sb_classical_binary_boosting_learner_adapter(
    "adaboost_binary",
    3L
  )
  logitboost <- sb_classical_binary_boosting_learner_adapter(
    "logitboost_binary",
    3L
  )
  local_mocked_bindings(
    sb_fastadaboost_fit = function(...) {
      warning("injected AdaBoost warning")
    },
    sb_catools_logitboost_fit = function(...) {
      warning("injected LogitBoost warning")
    }
  )

  expect_error(
    suppressWarnings(adaboost$fit(x, d$outcome)),
    class = "stormbreaker_invalid_classical_binary_boosting_fit"
  )
  expect_error(
    suppressWarnings(logitboost$fit(x, d$outcome)),
    class = "stormbreaker_invalid_classical_binary_boosting_fit"
  )
})

# ---- classical-binary-fixed-adapters ----

test_that("additional fixed binary learners keep distinct native methods", {
  learners <- c(
    "kernel_naive_bayes_binary",
    "rda_binary",
    "nnet_binary",
    "svm_polynomial_binary",
    "svm_sigmoid_binary",
    "bart_binary",
    "c50_binary",
    "kernel_weighted_knn_binary",
    "gaussian_process_binary",
    "pamr_binary",
    "neuralnet_multilayer_binary"
  )
  expect_setequal(
    intersect(sb_classical_binary_learner_ids(), learners),
    learners
  )
  expect_identical(
    sb_classical_binary_params(
      "kernel_naive_bayes_binary"
    ),
    list(
      prior = "empirical",
      usekernel = TRUE,
      fL = 0,
      adjust = 1
    )
  )
  expect_identical(
    sb_classical_binary_params("rda_binary"),
    list(
      prior = "empirical",
      gamma = 0.5,
      lambda = 0.5,
      estimate_error = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("nnet_binary"),
    list(
      size = 5L,
      decay = 0.01,
      maxit = 2000L,
      rang = 0.1,
      entropy = TRUE,
      standardize = TRUE,
      max_weights = 1000000L,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("svm_polynomial_binary"),
    list(
      kernel = "polynomial",
      cost = 1,
      gamma = "one_over_p",
      degree = 3L,
      coef0 = 0,
      standardize = TRUE,
      probability = TRUE,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("svm_sigmoid_binary"),
    list(
      kernel = "sigmoid",
      cost = 1,
      gamma = "one_over_p",
      coef0 = 0,
      standardize = TRUE,
      probability = TRUE,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("bart_binary"),
    list(
      num_trees = 50L,
      num_cutpoints = 100L,
      posterior_draws = 250L,
      burn_in = 100L,
      keep_every = 1L,
      prior_k = 2,
      tree_power = 2,
      tree_base = 0.95,
      sparse = FALSE,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("c50_binary"),
    list(
      trials = 10L,
      rules = FALSE,
      winnow = FALSE,
      early_stopping = FALSE,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("kernel_weighted_knn_binary"),
    list(k = 7L, distance = 2, kernel = "optimal", scale = TRUE)
  )
  expect_identical(
    sb_classical_binary_params("gaussian_process_binary"),
    list(
      kernel = "rbfdot",
      kernel_parameters = "automatic",
      scaled = TRUE,
      tolerance = 0.0005,
      cross_validation = 0L,
      keep_fitted = FALSE,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("pamr_binary"),
    list(
      threshold = 1,
      scale_sd = TRUE,
      offset_percent = 50,
      remove_zero_thresholds = FALSE,
      prior = "empirical"
    )
  )
  expect_identical(
    sb_classical_binary_params("neuralnet_multilayer_binary"),
    list(
      hidden = c(3L, 2L),
      threshold = 0.01,
      stepmax = 500000L,
      algorithm = "rprop+",
      error_function = "sse",
      activation = "logistic",
      linear_output = FALSE,
      standardize = TRUE,
      seed = 1L
    )
  )
})

test_that("additional fixed binary probabilities match native engines", {
  data <- sb_example_binary(n = 72L, seed = 6810L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:8, , drop = FALSE]
  learners <- c(
    "kernel_naive_bayes_binary",
    "rda_binary",
    "nnet_binary",
    "svm_polynomial_binary",
    "svm_sigmoid_binary"
  )

  for (learner in learners) {
    dependency <- sb_classical_binary_dependency(learner)
    if (!dependency$available) {
      next
    }
    params <- sb_classical_binary_params(learner)
    if (identical(learner, "cforest_binary")) {
      params$num_trees <- 20L
    }
    if (identical(learner, "bagging_binary")) {
      params$num_trees <- 20L
    }
    if (learner %in% c("gamboost_binary", "blackboost_binary")) {
      params$mstop <- 20L
    }
    if (identical(learner, "aorsf_binary")) {
      params$num_trees <- 20L
    }
    adapter <- sb_classical_binary_learner_adapter(learner, params)
    model <- adapter$fit(x, data$outcome)
    observed <- adapter$predict(model, assessment)
    native_engine <- if (identical(learner, "nnet_binary")) {
      model$engine$model
    } else {
      model$engine
    }
    expected_class <- switch(learner,
      kernel_naive_bayes_binary = "NaiveBayes",
      rda_binary = "rda",
      nnet_binary = "nnet",
      svm_polynomial_binary = "svm",
      svm_sigmoid_binary = "svm"
    )
    expect_s3_class(native_engine, expected_class)
    native <- switch(learner,
      kernel_naive_bayes_binary = stats::predict(
        model$engine,
        newdata = assessment,
        type = "raw"
      )$posterior[, "1"],
      rda_binary = get("predict.rda", envir = asNamespace("klaR"), inherits = FALSE)(
        model$engine,
        newdata = assessment
      )$posterior[, "1"],
      nnet_binary = {
        matrix_x <- as.matrix(assessment)
        scaled <- sweep(
          sweep(matrix_x, 2L, model$engine$center, "-"),
          2L,
          model$engine$scale,
          "/"
        )
        as.numeric(stats::predict(
          model$engine$model,
          newdata = scaled,
          type = "raw"
        ))
      },
      svm_polynomial_binary = {
        prediction <- stats::predict(
          model$engine,
          newdata = as.matrix(assessment),
          probability = TRUE
        )
        attr(prediction, "probabilities", exact = TRUE)[, "1"]
      },
      svm_sigmoid_binary = {
        prediction <- stats::predict(
          model$engine,
          newdata = as.matrix(assessment),
          probability = TRUE
        )
        attr(prediction, "probabilities", exact = TRUE)[, "1"]
      }
    )
    expect_equal(observed, unname(native), tolerance = 1e-14)
    if (identical(learner, "kernel_naive_bayes_binary")) {
      expect_identical(model$engine$usekernel, TRUE)
    }
    if (identical(learner, "rda_binary")) {
      expect_equal(
        unname(model$engine$regularization),
        c(0.5, 0.5),
        tolerance = 0
      )
    }
    if (identical(learner, "nnet_binary")) {
      expect_identical(model$engine$model$convergence, 0L)
    }
    if (identical(learner, "svm_polynomial_binary")) {
      expect_identical(model$engine$kernel, 1)
    }
    if (identical(learner, "svm_sigmoid_binary")) {
      expect_identical(model$engine$kernel, 3)
    }
  }
})

test_that("sixth-batch fixed binary probabilities match native engines", {
  data <- sb_example_binary(n = 72L, seed = 6815L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:8, , drop = FALSE]
  learners <- c(
    "bart_binary",
    "c50_binary",
    "kernel_weighted_knn_binary",
    "gaussian_process_binary",
    "pamr_binary",
    "neuralnet_multilayer_binary"
  )

  for (learner in learners) {
    dependency <- sb_classical_binary_dependency(learner)
    if (!dependency$available) {
      next
    }
    params <- sb_classical_binary_params(learner)
    if (identical(learner, "bart_binary")) {
      params$posterior_draws <- 40L
      params$burn_in <- 20L
    }
    adapter <- sb_classical_binary_learner_adapter(learner, params)
    model <- adapter$fit(x, data$outcome)
    observed <- adapter$predict(model, assessment)
    grouping <- factor(data$outcome, levels = c(0L, 1L))

    native <- switch(learner,
      bart_binary = {
        prediction <- NULL
        method <- get(
          "predict.pbart",
          envir = asNamespace("BART"),
          inherits = FALSE
        )
        invisible(utils::capture.output(
          prediction <- method(
            model$engine,
            sb_bart_prediction_matrix(model$engine, as.matrix(assessment)),
            mc.cores = 1L,
            openmp = FALSE
          ),
          type = "output"
        ))
        expect_s3_class(model$engine, "pbart")
        prediction$prob.test.mean
      },
      c50_binary = {
        expect_s3_class(model$engine, "C5.0")
        expect_identical(
          as.integer(unname(model$engine$trials["Actual"])),
          params$trials
        )
        get(
          "predict.C5.0",
          envir = asNamespace("C50"),
          inherits = FALSE
        )(model$engine, assessment, type = "prob")[, "1"]
      },
      kernel_weighted_knn_binary = {
        training <- x
        training$.stormbreaker_outcome <- grouping
        fitted <- kknn::kknn(
          .stormbreaker_outcome ~ .,
          train = training,
          test = assessment,
          k = params$k,
          distance = params$distance,
          kernel = params$kernel,
          scale = params$scale
        )
        expect_s3_class(fitted, "kknn")
        fitted$prob[, "1"]
      },
      gaussian_process_binary = {
        expect_s4_class(model$engine, "gausspr")
        kernlab::predict(
          model$engine,
          as.matrix(assessment),
          type = "probabilities"
        )[, "1"]
      },
      pamr_binary = {
        expect_s3_class(model$engine, "pamrtrained")
        expect_equal(model$engine$prior, as.numeric(table(grouping) / 72))
        pamr::pamr.predict(
          model$engine,
          newx = t(as.matrix(assessment)),
          threshold = params$threshold,
          type = "posterior"
        )[, "1"]
      },
      neuralnet_multilayer_binary = {
        expect_s3_class(model$engine$model, "nn")
        expect_length(model$engine$model$weights[[1L]], 3L)
        matrix_x <- as.matrix(assessment)
        standardized <- sweep(
          sweep(matrix_x, 2L, model$engine$center, "-"),
          2L,
          model$engine$scale,
          "/"
        )
        standardized <- as.data.frame(standardized)
        names(standardized) <- names(assessment)
        as.numeric(get(
          "predict.nn",
          envir = asNamespace("neuralnet"),
          inherits = FALSE
        )(model$engine$model, standardized))
      }
    )
    expect_equal(observed, unname(native), tolerance = 1e-14)
  }
})

test_that("binary BART replays after native constant-column removal", {
  skip_if_not_installed("BART", minimum_version = "2.9.0")

  data <- sb_example_binary(n = 60L, seed = 6816L)
  data$constant <- 1
  x <- data[c("gene_a", "constant", "gene_b")]
  params <- sb_classical_binary_params("bart_binary")
  params$posterior_draws <- 30L
  params$burn_in <- 15L
  adapter <- sb_classical_binary_learner_adapter(
    "bart_binary",
    params
  )
  model <- adapter$fit(x, data$outcome)
  assessment <- x[1:6, , drop = FALSE]

  observed <- adapter$predict(model, assessment)
  replay <- adapter$predict_state(
    adapter$extract_state(model),
    assessment
  )$probability

  expect_true(length(model$engine$rm.const) > 0L)
  expect_equal(replay, observed, tolerance = 0)
})

test_that("fixed classical binary learners fit and preserve predictions", {
  data <- sb_example_binary(n = 72L, seed = 6801L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  analysis <- data[predictors]
  assessment <- data[1:8, predictors, drop = FALSE]

  for (learner in sb_classical_binary_learner_ids()) {
    dependency <- sb_classical_binary_dependency(learner)
    if (!dependency$available) {
      next
    }
    adapter <- sb_classical_binary_learner_adapter(learner)
    model <- adapter$fit(analysis, data$outcome)
    probability <- adapter$predict(model, assessment)
    state_probability <- adapter$predict_state(
      adapter$extract_state(model),
      assessment
    )$probability

    expect_identical(adapter$id, learner)
    expect_length(probability, nrow(assessment))
    expect_identical(all(is.finite(probability)), TRUE)
    expect_identical(all(probability >= 0 & probability <= 1), TRUE)
    expect_equal(state_probability, probability, tolerance = 0)
  }
})

test_that("fixed classical binary learners execute through CV and final fit", {
  data <- sb_example_binary(n = 72L, seed = 6802L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    predictors = predictors
  )
  original_params <- sb_classical_binary_params
  local_mocked_bindings(
    sb_classical_binary_params = function(learner) {
      params <- original_params(learner)
      if (identical(learner, "cforest_binary")) {
        params$num_trees <- 20L
      }
      if (identical(learner, "bagging_binary")) {
        params$num_trees <- 20L
      }
      if (learner %in% c(
        "gamboost_binary",
        "blackboost_binary"
      )) {
        params$mstop <- 20L
      }
      if (identical(learner, "aorsf_binary")) {
        params$num_trees <- 20L
      }
      params
    }
  )

  for (learner in sb_classical_binary_learner_ids()) {
    dependency <- sb_classical_binary_dependency(learner)
    if (!dependency$available) {
      next
    }
    pipeline <- sb_pipeline(learner = learner)
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 6803L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:6, c("patient_id", predictors), drop = FALSE],
      id_col = "patient_id"
    )

    expect_identical(result$status, "success")
    expect_identical(result$complete_oof, TRUE)
    expect_identical(result$runs$learner_used, rep(learner, 3L))
    if (!is.na(dependency$package)) {
      field <- paste0(tolower(dependency$package), "_version")
      expect_identical(
        result$provenance[[field]],
        as.character(utils::packageVersion(dependency$package))
      )
      expect_identical(
        fit$provenance[[field]],
        as.character(utils::packageVersion(dependency$package))
      )
    }
    expect_identical(fit$status, "success")
    expect_identical(fit$learner_state$learner, learner)
    expect_identical(prediction$status, "success")
    expect_identical(
      prediction$predictions$learner,
      rep(learner, 6L)
    )
  }
})

test_that("fixed classical binary learners enforce prediction schema", {
  data <- sb_example_binary(n = 60L, seed = 6804L)
  adapter <- sb_classical_binary_learner_adapter("lda_binary")
  model <- adapter$fit(data[c("gene_a", "gene_b")], data$outcome)

  error <- capture_sb_error(
    adapter$predict(model, data[c("gene_b", "gene_a")])
  )

  expect_s3_class(error, "stormbreaker_learner_schema_mismatch")
})

test_that("weighted radial SVM derives fold-local inverse-frequency weights", {
  skip_if_not_installed("e1071", minimum_version = "1.7.0")

  data <- sb_example_binary(n = 72L, seed = 6805L)
  data$outcome <- c(rep(0L, 48L), rep(1L, 24L))
  x <- data[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_classical_binary_learner_adapter(
    "svm_radial_weighted_binary"
  )
  original_svm <- e1071::svm
  observed_weights <- NULL

  testthat::with_mocked_bindings(
    svm = function(..., class.weights = NULL) {
      observed_weights <<- class.weights
      original_svm(..., class.weights = class.weights)
    },
    adapter$fit(x, data$outcome),
    .package = "e1071"
  )

  counts <- table(factor(data$outcome, levels = c(0L, 1L)))
  expected <- setNames(
    length(data$outcome) / (2 * as.numeric(counts)),
    names(counts)
  )
  expect_equal(observed_weights, expected, tolerance = 0)
})

test_that("LDA collinearity is a typed learner failure", {
  data <- sb_example_binary(n = 60L, seed = 6805L)
  data$duplicate_gene_a <- data$gene_a
  predictors <- c("gene_a", "gene_b", "duplicate_gene_a")
  adapter <- sb_classical_binary_learner_adapter("lda_binary")

  error <- capture_sb_error(
    adapter$fit(data[predictors], data$outcome)
  )
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    predictors = predictors
  )
  pipeline <- sb_pipeline(learner = "lda_binary")
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 6806L)
  fit <- sb_fit_final(task, pipeline)

  expect_s3_class(error, "stormbreaker_lda_collinear_predictors")
  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_lda_collinear_predictors"
  )
  expect_identical(fit$status, "failed")
  expect_identical(fit$error_stage, "learner_fit")
  expect_identical(
    fit$error_class,
    "stormbreaker_lda_collinear_predictors"
  )
})

test_that("fixed neural-network nonconvergence is a typed learner failure", {
  skip_if_not_installed("nnet", minimum_version = "7.3")

  data <- sb_example_binary(n = 60L, seed = 6811L)
  params <- sb_classical_binary_params("nnet_binary")
  params$maxit <- 1L
  error <- capture_sb_error(sb_fit_classical_binary(
    data[c("gene_a", "gene_b", "gene_c")],
    data$outcome,
    "nnet_binary",
    params
  ))

  expect_s3_class(error, "stormbreaker_nnet_nonconvergence")
})

test_that("multilayer neural-network nonconvergence is typed", {
  skip_if_not_installed("neuralnet", minimum_version = "1.44.0")

  data <- sb_example_binary(n = 60L, seed = 6817L)
  params <- sb_classical_binary_params(
    "neuralnet_multilayer_binary"
  )
  params$stepmax <- 1L
  error <- capture_sb_error(sb_fit_classical_binary(
    data[c("gene_a", "gene_b", "gene_c")],
    data$outcome,
    "neuralnet_multilayer_binary",
    params
  ))

  expect_s3_class(error, "stormbreaker_neuralnet_nonconvergence")
})

test_that("RDA prediction is isolated from conflicting S3 registrations", {
  skip_if_not_installed("klaR", minimum_version = "1.7.0")

  data <- sb_example_binary(n = 60L, seed = 6812L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_classical_binary_learner_adapter("rda_binary")
  model <- adapter$fit(x, data$outcome)
  expected <- get(
    "predict.rda",
    envir = asNamespace("klaR"),
    inherits = FALSE
  )(model$engine, newdata = x[1:6, , drop = FALSE])$posterior[, "1"]
  original <- getS3method("predict", "rda")
  withr::defer(registerS3method(
    "predict",
    "rda",
    original,
    envir = globalenv()
  ))
  registerS3method(
    "predict",
    "rda",
    function(object, ...) stop("conflicting predict.rda", call. = FALSE),
    envir = globalenv()
  )
  generic_error <- tryCatch(
    stats::predict(model$engine, newdata = x[1:6, , drop = FALSE]),
    error = identity
  )
  observed <- adapter$predict(model, x[1:6, , drop = FALSE])

  expect_s3_class(generic_error, "error")
  expect_equal(observed, unname(expected), tolerance = 1e-14)
})

test_that("saved fixed learner engines omit non-prediction training payloads", {
  data <- sb_example_binary(n = 60L, seed = 6807L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  analysis <- data[predictors]
  engine_for <- function(learner) {
    adapter <- sb_classical_binary_learner_adapter(learner)
    model <- adapter$fit(analysis, data$outcome)
    unserialize(adapter$extract_state(model)$engine_raw)
  }

  knn <- engine_for("knn_binary")
  expect_identical(nrow(knn$training), nrow(analysis))
  expect_identical(knn$outcome, data$outcome)

  if (sb_classical_binary_dependency(
    "decision_tree_binary"
  )$available) {
    tree <- engine_for("decision_tree_binary")
    expect_identical(
      intersect(c("model", "x", "y", "where"), names(tree)),
      character()
    )
    expect_identical(environment(tree$terms), baseenv())
  }

  if (sb_classical_binary_dependency(
    "random_forest_binary"
  )$available) {
    forest <- engine_for("random_forest_binary")
    expect_identical(
      intersect(
        c(
          "predicted", "err.rate", "confusion", "votes", "oob.times",
          "y", "test", "inbag", "localImportance", "proximity"
        ),
        names(forest)
      ),
      character()
    )
  }

  if (sb_classical_binary_dependency(
    "svm_linear_binary"
  )$available) {
    svm <- engine_for("svm_linear_binary")
    expect_identical(
      intersect(
        c("call", "index", "na.action", "fitted", "decision.values"),
        names(svm)
      ),
      character()
    )
  }

  if (sb_classical_binary_dependency(
    "kernel_naive_bayes_binary"
  )$available) {
    kernel_nb <- engine_for("kernel_naive_bayes_binary")
    expect_identical(nrow(kernel_nb$x), 0L)
    expect_identical("call" %in% names(kernel_nb), FALSE)
  }

  if (sb_classical_binary_dependency("rda_binary")$available) {
    rda <- engine_for("rda_binary")
    expect_identical(
      intersect(c("call", "error.rate", "converged", "iter"), names(rda)),
      character()
    )
  }

  if (sb_classical_binary_dependency("nnet_binary")$available) {
    neural_net <- engine_for("nnet_binary")
    expect_identical(
      intersect(
        c("call", "fitted.values", "residuals"),
        names(neural_net$model)
      ),
      character()
    )
  }

  if (sb_classical_binary_dependency("bart_binary")$available) {
    bart <- engine_for("bart_binary")
    expect_true(all(names(bart) %in% c(
      "treedraws", "rm.const", "binaryOffset", "offset"
    )))
    expect_false(any(c("yhat.train", "prob.train", "varcount") %in%
      names(bart)))
  }

  if (sb_classical_binary_dependency("c50_binary")$available) {
    c50 <- engine_for("c50_binary")
    expect_identical(
      intersect(c("call", "output", "boostResults", "size"), names(c50)),
      character()
    )
  }

  if (sb_classical_binary_dependency(
    "kernel_weighted_knn_binary"
  )$available) {
    weighted_knn <- engine_for("kernel_weighted_knn_binary")
    expect_identical(nrow(weighted_knn$training), nrow(analysis))
    expect_identical(as.integer(weighted_knn$outcome), data$outcome + 1L)
  }

  if (sb_classical_binary_dependency(
    "gaussian_process_binary"
  )$available) {
    gaussian_process <- engine_for("gaussian_process_binary")
    expect_s4_class(gaussian_process, "gausspr")
    expect_identical(methods::slot(gaussian_process, "fitted"), NA)
  }

  if (sb_classical_binary_dependency("pamr_binary")$available) {
    pamr <- engine_for("pamr_binary")
    expect_identical(
      intersect(c("y", "yhat", "prob", "errors", "call"), names(pamr)),
      character()
    )
  }

  if (sb_classical_binary_dependency(
    "neuralnet_multilayer_binary"
  )$available) {
    neuralnet <- engine_for("neuralnet_multilayer_binary")
    expect_setequal(
      names(neuralnet$model),
      c("model.list", "weights", "act.fct", "linear.output")
    )
  }

  if (sb_classical_binary_dependency("gbm_binary")$available) {
    gbm <- engine_for("gbm_binary")
    expect_identical(
      intersect(
        c("fit", "train.error", "valid.error", "oobag.improve"),
        names(gbm)
      ),
      character()
    )
  }
})

test_that("CForest safely predicts more rows than its training data", {
  skip_if_not_installed("partykit", minimum_version = "1.2.0")

  data <- sb_example_binary(n = 72L, seed = 6821L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  params <- sb_classical_binary_params("cforest_binary")
  params$num_trees <- 20L
  adapter <- sb_classical_binary_learner_adapter(
    "cforest_binary", params
  )
  model <- adapter$fit(x, data$outcome)
  training_data <- model$engine$data
  rows <- rep(rev(seq_len(nrow(x))), length.out = 2L * nrow(x) + 5L)
  assessment <- x[rows, , drop = FALSE]
  chunks <- split(seq_len(nrow(assessment)), ceiling(
    seq_len(nrow(assessment)) / nrow(x)
  ))
  native_predict <- stats::predict
  expected <- unlist(lapply(chunks, function(index) {
    native_predict(
      model$engine, newdata = assessment[index, , drop = FALSE], type = "prob"
    )[, "1"]
  }), use.names = FALSE)
  native_batch_sizes <- integer()
  local_mocked_bindings(
    predict = function(object, newdata, ...) {
      if (inherits(object, "cforest")) {
        native_batch_sizes <<- c(native_batch_sizes, nrow(newdata))
        if (nrow(newdata) > nrow(object$data)) {
          stop("The synthetic test intercepted an unsafe native CForest batch.")
        }
      }
      native_predict(object, newdata = newdata, ...)
    },
    .package = "stats"
  )

  observed <- adapter$predict(model, assessment)
  state <- sb_extract_classical_binary_learner_state(model)
  replayed <- sb_predict_classical_binary_learner_state(
    state, assessment
  )$probability
  expect_equal(observed, expected, tolerance = 0)
  expect_equal(replayed, expected, tolerance = 0)
  expect_identical(native_batch_sizes, rep(unname(lengths(chunks)), 2L))
  expect_identical(model$params, params)
  expect_identical(model$engine$data, training_data)
  expect_equal(
    adapter$predict(model, assessment[1:8, , drop = FALSE]),
    expected[1:8],
    tolerance = 0
  )
  expect_identical(tail(native_batch_sizes, 1L), 8L)
})

test_that("new fixed learners run with learned and fixed feature selectors", {
  data <- sb_example_binary(n = 72L, seed = 6808L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    predictors = c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "lda_binary",
      selector_n_features = 2L
    ),
    sb_pipeline(
      selector = "fixed_features",
      learner = "decision_tree_binary",
      selector_features = c("gene_b", "gene_a"),
      feature_set_id = "two_genes"
    ),
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "kernel_naive_bayes_binary",
      selector_n_features = 2L
    ),
    sb_pipeline(
      selector = "fixed_features",
      learner = "nnet_binary",
      selector_features = c("gene_b", "gene_a"),
      feature_set_id = "neural_two_genes"
    ),
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "bart_binary",
      selector_n_features = 2L
    ),
    sb_pipeline(
      selector = "fixed_features",
      learner = "pamr_binary",
      selector_features = c("gene_b", "gene_a"),
      feature_set_id = "pamr_two_genes"
    )
  )

  for (pipeline in pipelines) {
    dependency <- sb_classical_binary_dependency(pipeline$learner)
    if (!dependency$available) {
      next
    }
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 6809L)
    fit <- sb_fit_final(task, pipeline)

    expect_identical(result$status, "success")
    expect_identical(fit$status, "success")
    expect_length(fit$required_features, 2L)
  }
})

# ---- classical-expansion-batch9-adapters ----

test_that("batch-nine binary learners declare distinct native methods", {
  expected <- c(
    earth_binary = "earth",
    fda_binary = "mda",
    mda_binary = "mda",
    liblinear_l2_logistic_binary = "LiblineaR",
    evtree_binary = "evtree",
    sda_binary = "sda",
    rotation_forest_binary = "rotationForest",
    rrf_binary = "RRF"
  )

  expect_true(all(names(expected) %in% sb_classical_binary_learner_ids()))
  for (learner in names(expected)) {
    dependency <- sb_classical_binary_dependency(learner)
    expect_identical(dependency$package, unname(expected[[learner]]))
    expect_identical(
      dependency$available,
      requireNamespace(dependency$package, quietly = TRUE)
    )
  }

  expect_identical(
    sb_classical_binary_params("earth_binary"),
    list(
      degree = 1L,
      nk_rule = "max_21_2p_plus1",
      nprune = 15L,
      pmethod = "backward",
      glm_family = "binomial_logit",
      nfold = 0L,
      keep_xy = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("fda_binary"),
    list(
      method = "polyreg",
      dimension = 1L,
      degree = 1L,
      monomial = FALSE,
      keep_fitted = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("mda_binary"),
    list(
      subclasses_per_class = 2L,
      dimension = 3L,
      iterations = 10L,
      method = "polyreg",
      degree = 1L,
      keep_fitted = FALSE,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params(
      "liblinear_l2_logistic_binary"
    ),
    list(
      solver_type = 0L,
      cost = 1,
      tolerance = 0.01,
      bias = 1,
      standardize = TRUE,
      probability = TRUE
    )
  )
  expect_identical(
    sb_classical_binary_params("evtree_binary"),
    list(
      minbucket = 7L,
      minsplit = 20L,
      maxdepth = 5L,
      iterations = 1000L,
      population_size = 50L,
      complexity_weight = 1,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("sda_binary"),
    list(
      covariance_shrinkage = "analytic",
      variance_shrinkage = "analytic",
      class_frequency_shrinkage = "analytic",
      diagonal = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("rotation_forest_binary"),
    list(
      num_trees = 10L,
      subset_rule = "round_p_over_3_bounded_by_1_and_floor_p_over_2",
      cp = 0.01,
      minsplit = 10L,
      maxdepth = 10L,
      xval = 0L,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("rrf_binary"),
    list(
      num_trees = 500L,
      mtry_rule = "floor_sqrt_p",
      min_node_size = 1L,
      replace = TRUE,
      regularization_mode = 1L,
      regularization_coefficient = 0.5,
      seed = 1L
    )
  )
})

test_that("batch-nine binary predictions equal native engine predictions", {
  data <- sb_example_binary(n = 84L, seed = 8511L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:8, , drop = FALSE]
  learners <- c(
    "earth_binary",
    "fda_binary",
    "mda_binary",
    "liblinear_l2_logistic_binary",
    "evtree_binary",
    "sda_binary",
    "rotation_forest_binary",
    "rrf_binary"
  )

  learners <- learners[vapply(learners, function(learner) {
    sb_classical_binary_dependency(learner)$available
  }, logical(1))]
  if (!length(learners)) skip("No optional engines for this batch are installed")

  for (learner in learners) {
    params <- sb_classical_binary_params(learner)
    if (identical(learner, "evtree_binary")) {
      params$iterations <- 200L
      params$population_size <- 20L
    }
    if (identical(learner, "rotation_forest_binary")) {
      params$num_trees <- 3L
    }
    if (identical(learner, "rrf_binary")) {
      params$num_trees <- 40L
    }
    adapter <- sb_classical_binary_learner_adapter(learner, params)
    model <- adapter$fit(x, data$outcome)
    observed <- adapter$predict(model, assessment)

    native <- switch(learner,
      earth_binary = {
        expect_s3_class(model$engine, "earth")
        stats::predict(model$engine, assessment, type = "response")
      },
      fda_binary = {
        expect_s3_class(model$engine, "fda")
        stats::predict(model$engine, assessment, type = "posterior")[, "1"]
      },
      mda_binary = {
        expect_s3_class(model$engine, "mda")
        stats::predict(model$engine, assessment, type = "posterior")[, "1"]
      },
      liblinear_l2_logistic_binary = {
        expect_s3_class(model$engine$model, "LiblineaR")
        matrix_x <- as.matrix(assessment)
        standardized <- sweep(
          sweep(matrix_x, 2L, model$engine$center, "-"),
          2L,
          model$engine$scale,
          "/"
        )
        stats::predict(
          model$engine$model,
          standardized,
          proba = TRUE,
          decisionValues = FALSE
        )$probabilities[, "1"]
      },
      evtree_binary = {
        expect_s3_class(model$engine, "constparty")
        stats::predict(model$engine, assessment, type = "prob")[, "1"]
      },
      sda_binary = {
        expect_s3_class(model$engine, "sda")
        sda::predict.sda(
          model$engine,
          as.matrix(assessment),
          verbose = FALSE
        )$posterior[, "1"]
      },
      rotation_forest_binary = {
        expect_s3_class(model$engine, "rotationForest")
        stats::predict(model$engine, assessment)
      },
      rrf_binary = {
        expect_s3_class(model$engine, "RRF")
        stats::predict(
          model$engine,
          as.matrix(assessment),
          type = "prob"
        )[, "1"]
      }
    )
    state <- adapter$extract_state(model)
    replay <- adapter$predict_state(
      unserialize(serialize(state, NULL)),
      assessment
    )$probability

    expect_equal(observed, unname(as.numeric(native)), tolerance = 1e-13)
    expect_equal(replay, observed, tolerance = 0)
  }
})

test_that("batch-nine prediction states omit avoidable fitted-row payloads", {
  for (package in c("earth", "mda", "rotationForest", "RRF", "evtree")) {
    skip_if_not_installed(package)
  }
  data <- sb_example_binary(n = 72L, seed = 8512L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  state_engine <- function(learner, params = NULL) {
    if (is.null(params)) {
      params <- sb_classical_binary_params(learner)
    }
    adapter <- sb_classical_binary_learner_adapter(learner, params)
    model <- adapter$fit(x, data$outcome)
    list(
      state = adapter$extract_state(model),
      engine = unserialize(adapter$extract_state(model)$engine_raw)
    )
  }

  earth <- state_engine("earth_binary")
  expect_identical(
    intersect(
      c("bx", "fitted.values", "residuals", "glm.list", "weights", "call"),
      names(earth$engine)
    ),
    character()
  )

  fda <- state_engine("fda_binary")
  expect_identical(
    intersect(c("x", "residuals", "call"), names(fda$engine$fit)),
    character()
  )

  mda <- state_engine("mda_binary")
  expect_false("weights" %in% names(mda$engine))

  params <- sb_classical_binary_params("rotation_forest_binary")
  params$num_trees <- 3L
  rotation <- state_engine("rotation_forest_binary", params)
  expect_true(all(vapply(rotation$engine$models, function(tree) {
    identical(intersect(c("model", "x", "y", "where", "call"), names(tree)),
      character()) && identical(environment(tree$terms), baseenv())
  }, logical(1))))

  params <- sb_classical_binary_params("rrf_binary")
  params$num_trees <- 30L
  rrf <- state_engine("rrf_binary", params)
  expect_identical(
    intersect(
      c("predicted", "votes", "oob.times", "y", "inbag", "call"),
      names(rrf$engine)
    ),
    character()
  )

  evtree <- state_engine("evtree_binary")
  expect_identical(
    evtree$state$training_data_retention,
    "required_by_evtree_terminal_class_probability_state"
  )
})

test_that("MDA fixes local initialization without changing caller RNG", {
  skip_if_not_installed("mda")
  data <- sb_example_binary(n = 84L, seed = 8513L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:8, , drop = FALSE]
  adapter <- sb_classical_binary_learner_adapter("mda_binary")

  set.seed(19L)
  caller_state <- .Random.seed
  first <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, caller_state)

  set.seed(991L)
  second <- adapter$fit(x, data$outcome)
  expect_equal(
    adapter$predict(first, assessment),
    adapter$predict(second, assessment),
    tolerance = 0
  )
  expect_identical(first$params$seed, 1L)
})

# ---- classical-native-binary-adapters ----

test_that("batch-eight binary learners declare distinct native methods", {
  expected <- c(
    ctree_binary = "partykit",
    cforest_binary = "partykit",
    bagging_binary = "ipred",
    gamboost_binary = "mboost",
    blackboost_binary = "mboost",
    loclda_binary = "klaR",
    aorsf_binary = "aorsf"
  )

  expect_true(all(names(expected) %in% sb_classical_binary_learner_ids()))
  for (learner in names(expected)) {
    dependency <- sb_classical_binary_dependency(learner)
    expect_identical(dependency$package, unname(expected[[learner]]))
    expect_identical(
      dependency$available,
      requireNamespace(dependency$package, quietly = TRUE)
    )
  }

  expect_identical(
    sb_classical_binary_params("ctree_binary"),
    list(
      teststat = "quad",
      testtype = "Bonferroni",
      mincriterion = 0.95,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = Inf,
      saveinfo = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("cforest_binary"),
    list(
      num_trees = 500L,
      mtry_rule = "ceiling_sqrt_p",
      replace = TRUE,
      sample_fraction = 1,
      teststat = "quad",
      testtype = "Univ",
      mincriterion = 0,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = Inf,
      saveinfo = FALSE,
      seed = 1L,
      cores = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("bagging_binary"),
    list(
      num_trees = 100L,
      cp = 0,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = 30L,
      xval = 0L,
      seed = 1L,
      coob = FALSE,
      keep_x = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("gamboost_binary"),
    list(
      mstop = 100L,
      nu = 0.1,
      risk = "inbag",
      family = "Binomial_logit",
      baselearner = "bbs",
      dfbase = 4,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("blackboost_binary"),
    list(
      mstop = 100L,
      nu = 0.1,
      risk = "inbag",
      family = "Binomial_logit",
      teststat = "quad",
      testtype = "Teststatistic",
      mincriterion = 0,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = 2L,
      saveinfo = FALSE,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("loclda_binary"),
    list(
      k_max = 30L,
      k_rule = "dimension_aware",
      weighted_apriori = TRUE,
      prior = "local_empirical"
    )
  )
  expect_identical(
    sb_classical_binary_params("aorsf_binary"),
    list(
      num_trees = 500L,
      num_split = 5L,
      num_retry = 3L,
      mtry_rule = "ceiling_sqrt_p",
      sample_with_replacement = TRUE,
      sample_fraction = 0.632,
      control_method = "glm",
      control_scale_x = TRUE,
      control_max_iter = 1L,
      seed = 1L,
      threads = 1L,
      importance = "none",
      attach_data = FALSE
    )
  )
})

test_that("batch-eight binary predictions equal native engine predictions", {
  data <- sb_example_binary(n = 72L, seed = 8421L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:8, , drop = FALSE]
  learners <- c(
    "ctree_binary",
    "cforest_binary",
    "bagging_binary",
    "gamboost_binary",
    "blackboost_binary",
    "loclda_binary",
    "aorsf_binary"
  )

  learners <- learners[vapply(learners, function(learner) {
    sb_classical_binary_dependency(learner)$available
  }, logical(1))]
  if (!length(learners)) skip("No optional engines for this batch are installed")

  for (learner in learners) {
    params <- sb_classical_binary_params(learner)
    if (identical(learner, "cforest_binary")) {
      params$num_trees <- 30L
    }
    if (identical(learner, "bagging_binary")) {
      params$num_trees <- 20L
    }
    if (learner %in% c("gamboost_binary", "blackboost_binary")) {
      params$mstop <- 20L
    }
    if (identical(learner, "aorsf_binary")) {
      params$num_trees <- 30L
    }
    adapter <- sb_classical_binary_learner_adapter(learner, params)
    model <- adapter$fit(x, data$outcome)
    observed <- adapter$predict(model, assessment)

    native <- switch(learner,
      ctree_binary = {
        expect_s3_class(model$engine, "constparty")
        stats::predict(model$engine, assessment, type = "prob")[, "1"]
      },
      cforest_binary = {
        expect_s3_class(model$engine, "cforest")
        stats::predict(model$engine, assessment, type = "prob")[, "1"]
      },
      bagging_binary = {
        expect_s3_class(model$engine, "classbagg")
        stats::predict(
          model$engine,
          assessment,
          type = "prob",
          aggregation = "average"
        )[, "1"]
      },
      gamboost_binary = {
        expect_s3_class(model$engine, "gamboost")
        stats::predict(model$engine, assessment, type = "response")
      },
      blackboost_binary = {
        expect_s3_class(model$engine, "blackboost")
        stats::predict(model$engine, assessment, type = "response")
      },
      loclda_binary = {
        expect_s3_class(model$engine, "loclda")
        get(
          "predict.loclda",
          envir = asNamespace("klaR"),
          inherits = FALSE
        )(model$engine, assessment)$posterior[, "1"]
      },
      aorsf_binary = {
        expect_s3_class(model$engine, "ObliqueForestClassification")
        stats::predict(
          model$engine,
          new_data = assessment,
          pred_type = "prob",
          n_thread = 1L
        )[, "1"]
      }
    )
    state <- adapter$extract_state(model)
    replay <- adapter$predict_state(
      unserialize(serialize(state, NULL)),
      assessment
    )$probability

    expect_equal(observed, unname(as.numeric(native)), tolerance = 1e-14)
    expect_equal(replay, observed, tolerance = 0)
  }
})

test_that("lazy native binary states disclose required analysis data", {
  for (package in c("partykit", "mboost", "klaR", "aorsf")) {
    skip_if_not_installed(package)
  }
  data <- sb_example_binary(n = 60L, seed = 8422L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  expected <- c(
    ctree_binary =
      "required_by_partykit_class_probability_state",
    cforest_binary =
      "required_by_partykit_class_probability_state",
    gamboost_binary =
      "required_by_mboost_prediction_state",
    blackboost_binary =
      "required_by_mboost_prediction_state",
    loclda_binary =
      "required_by_klar_localized_lda_prediction_state"
  )

  for (learner in names(expected)) {
    params <- sb_classical_binary_params(learner)
    if (identical(learner, "cforest_binary")) {
      params$num_trees <- 20L
    }
    if (learner %in% c("gamboost_binary", "blackboost_binary")) {
      params$mstop <- 15L
    }
    adapter <- sb_classical_binary_learner_adapter(learner, params)
    state <- adapter$extract_state(adapter$fit(x, data$outcome))
    expect_identical(state$training_data_retention, unname(expected[[learner]]))
  }

  params <- sb_classical_binary_params("aorsf_binary")
  params$num_trees <- 20L
  adapter <- sb_classical_binary_learner_adapter(
    "aorsf_binary",
    params
  )
  model <- adapter$fit(x, data$outcome)
  state <- adapter$extract_state(model)
  expect_null(model$engine$data)
  expect_false("training_data_retention" %in% names(state))
})

test_that("bagging state drops fitted-row outcomes and replays exactly", {
  skip_if_not_installed("ipred")
  data <- sb_example_binary(n = 60L, seed = 8423L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:8, , drop = FALSE]
  params <- sb_classical_binary_params("bagging_binary")
  params$num_trees <- 20L
  adapter <- sb_classical_binary_learner_adapter(
    "bagging_binary",
    params
  )
  model <- adapter$fit(x, data$outcome)
  expected <- adapter$predict(model, assessment)
  state <- adapter$extract_state(model)
  engine <- unserialize(state$engine_raw)
  replay <- adapter$predict_state(state, assessment)$probability
  has_full_row_payload <- function(value, n) {
    if (is.environment(value) || is.function(value)) {
      return(FALSE)
    }
    if ((is.data.frame(value) || is.matrix(value)) && nrow(value) == n) {
      return(TRUE)
    }
    if (is.atomic(value) && is.null(dim(value)) && length(value) == n) {
      return(TRUE)
    }
    if (is.list(value) || is.pairlist(value) || is.call(value)) {
      return(any(vapply(as.list(value), has_full_row_payload, logical(1), n)))
    }
    FALSE
  }

  expect_identical(levels(engine$y), c("0", "1"))
  expect_length(engine$y, 0L)
  expect_null(engine$X)
  expect_null(engine$OOB)
  expect_null(engine$call)
  expect_true(all(vapply(engine$mtrees, function(tree) {
    is.null(tree$bindx) && is.null(tree$btree$y) &&
      is.null(tree$btree$where) && is.null(tree$btree$call) &&
      identical(environment(tree$btree$terms), baseenv())
  }, logical(1))))
  expect_false(has_full_row_payload(engine, nrow(x)))
  expect_false("training_data_retention" %in% names(state))
  expect_equal(replay, expected, tolerance = 0)
})

# ---- lvq-binary-fixed ----

test_that("four LVQ update algorithms are admitted as binary learners", {
  learners <- c(
    "lvq1_binary",
    "olvq1_binary",
    "lvq21_binary",
    "lvq3_binary"
  )
  expect_setequal(sb_lvq_binary_learner_ids(), learners)

  expected_methods <- c("lvq1", "olvq1", "lvq2.1", "lvq3")
  expected_iterations <- c(100L, 40L, 100L, 100L)
  expected_alpha <- c(0.03, 0.3, 0.03, 0.03)
  expected_selectors <- expected_fixed_learner_selectors("binary")
  registry <- sb_pipeline_registry()

  for (index in seq_along(learners)) {
    learner <- learners[[index]]
    params <- sb_classical_binary_params(learner)
    expect_identical(params$method, expected_methods[[index]])
    expect_identical(
      params$iterations_per_prototype,
      expected_iterations[[index]]
    )
    expect_identical(params$alpha, expected_alpha[[index]])
    expect_identical(params$codebook_size_rule, "max_two_native_default")
    expect_identical(params$prototype_allocation, "equal_prior_target")
    expect_identical(params$initialization_k, 5L)
    expect_identical(params$standardize, TRUE)
    expect_identical(params$seed, 1L)
    expect_identical(
      params$score_mapping,
      "logistic_nearest_prototype_squared_distance_margin"
    )
    if (learner %in% c("lvq21_binary", "lvq3_binary")) {
      expect_identical(params$window, 0.3)
    } else {
      expect_null(params$window)
    }
    if (identical(learner, "lvq3_binary")) {
      expect_identical(params$epsilon, 0.1)
    } else {
      expect_null(params$epsilon)
    }

    dependency <- sb_classical_binary_dependency(learner)
    expect_identical(dependency$package, "class")
    expect_identical(dependency$available, TRUE)

    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(unique(rows$task_type), "binary")
    expect_setequal(rows$selector, expected_selectors)
    expect_identical(unique(rows$tuning_semantics), "none")
    expect_identical(
      all(grepl("class", rows$optional_dependency, fixed = TRUE)),
      TRUE
    )
  }
})

test_that("portable LVQ states preserve the native nearest-prototype boundary", {
  skip_if_not_installed("class", minimum_version = "7.3")

  data <- make_sb_binary_fixture(n = 84L, seed = 11101L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene one", "gene-two", "gene_c")
  assessment <- x[c(3L, 19L, 47L, 72L), , drop = FALSE]
  learners <- sb_lvq_binary_learner_ids()

  withr::local_seed(11102L)
  rng_before <- .Random.seed
  for (learner in learners) {
    adapter <- sb_classical_binary_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, rng_before)
    repeated <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, rng_before)
    expect_equal(repeated$engine$prototypes, model$engine$prototypes)

    params <- sb_classical_binary_params(learner)
    standardized_training <- sweep(
      sweep(as.matrix(x), 2L, model$engine$center, "-"),
      2L,
      model$engine$scale,
      "/"
    )
    native <- withr::with_seed(params$seed, {
      size <- max(
        2L,
        as.integer(min(
          round(0.4 * 2 * (2 - 1 + ncol(x) / 2), 0),
          nrow(x)
        ))
      )
      initial <- class::lvqinit(
        standardized_training,
        factor(data$outcome, levels = c(0L, 1L)),
        size = size,
        prior = c(0.5, 0.5),
        k = params$initialization_k
      )
      iterations <- params$iterations_per_prototype * nrow(initial$x)
      switch(params$method,
        lvq1 = class::lvq1(
          standardized_training,
          factor(data$outcome, levels = c(0L, 1L)),
          initial,
          niter = iterations,
          alpha = params$alpha
        ),
        olvq1 = class::olvq1(
          standardized_training,
          factor(data$outcome, levels = c(0L, 1L)),
          initial,
          niter = iterations,
          alpha = params$alpha
        ),
        `lvq2.1` = class::lvq2(
          standardized_training,
          factor(data$outcome, levels = c(0L, 1L)),
          initial,
          niter = iterations,
          alpha = params$alpha,
          win = params$window
        ),
        lvq3 = class::lvq3(
          standardized_training,
          factor(data$outcome, levels = c(0L, 1L)),
          initial,
          niter = iterations,
          alpha = params$alpha,
          win = params$window,
          epsilon = params$epsilon
        )
      )
    })
    expect_equal(
      model$engine$prototypes,
      unname(as.matrix(native$x)),
      tolerance = 0
    )
    expect_identical(
      model$engine$prototype_class,
      unname(as.character(native$cl))
    )

    observed <- adapter$predict(model, assessment)
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, assessment)$probability
    engine <- unserialize(state$engine_raw)
    standardized <- sweep(
      sweep(as.matrix(assessment), 2L, engine$center, "-"),
      2L,
      engine$scale,
      "/"
    )
    codebook <- list(
      x = engine$prototypes,
      cl = factor(engine$prototype_class, levels = c("0", "1"))
    )
    native_class <- class::lvqtest(codebook, standardized)
    squared_distance <- outer(
      rowSums(standardized^2),
      rowSums(engine$prototypes^2),
      "+"
    ) - 2 * tcrossprod(standardized, engine$prototypes)
    squared_distance <- pmax(squared_distance, 0)
    distance_zero <- apply(
      squared_distance[, engine$prototype_class == "0", drop = FALSE],
      1L,
      min
    )
    distance_one <- apply(
      squared_distance[, engine$prototype_class == "1", drop = FALSE],
      1L,
      min
    )
    expected <- stats::plogis(distance_zero - distance_one)

    expect_equal(observed, unname(expected), tolerance = 1e-12)
    expect_equal(replay, observed, tolerance = 0)
    expect_identical(
      as.integer(observed > 0.5),
      as.integer(as.character(native_class))
    )
    expect_setequal(
      names(engine),
      c(
        "center", "scale", "prototypes", "prototype_class",
        "codebook_size", "iterations"
      )
    )
    expect_null(rownames(engine$prototypes))
    expect_identical(
      any(c("outcome", "training", "native_fit") %in% names(engine)),
      FALSE
    )
    expect_identical("training_data_retention" %in% names(state), FALSE)
  }
})

test_that("the LVQ initialization rule supports a one-feature analysis fit", {
  skip_if_not_installed("class", minimum_version = "7.3")

  outcome <- rep(c(0L, 1L), each = 30L)
  x <- data.frame(signal = c(
    seq(-3, -1, length.out = 30L),
    seq(1, 3, length.out = 30L)
  ))
  for (learner in sb_lvq_binary_learner_ids()) {
    adapter <- sb_classical_binary_learner_adapter(learner)
    model <- adapter$fit(x, outcome)
    probability <- adapter$predict(model, x[c(5L, 55L), , drop = FALSE])

    expect_identical(model$engine$codebook_size, 2L)
    expect_setequal(model$engine$prototype_class, c("0", "1"))
    expect_length(probability, 2L)
    expect_identical(all(is.finite(probability)), TRUE)
  }
})

test_that("LVQ learners execute through fold-local CV and final fit", {
  skip_if_not_installed("class", minimum_version = "7.3")

  data <- make_sb_binary_fixture(n = 72L, seed = 11111L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )

  for (learner in sb_lvq_binary_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 2L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 11112L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:6, c("patient_id", "gene_c", "gene_a", "gene_b")],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(sum(comparison$runs$status == "success"), 3L)
    expect_identical("class_version" %in% names(comparison$provenance), TRUE)
    expect_identical(fit$status, "success")
    expect_identical(fit$learner_state$learner, learner)
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$probability, 6L)
  }
})
})

# ---- model-logistic ----
local({
# ---- extended-logistic-parameters ----

test_that("five additional logistic predictors retain exact original-unit explanations", {
  d <- sb_example_binary(n = 160L, seed = 108L)
  d$gene_a <- d$gene_a + 10
  d$gene_b <- 3 * d$gene_b - 5
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  dependencies <- c(stepglm_binary = "stats", glmboost_binary = "mboost",
    plsrglm_binary = "plsRglm", brglm2_logistic = "brglm2", bayesglm_logistic_cauchy = "arm")
  for (id in names(dependencies)) {
    skip_if_not_installed(dependencies[[id]])
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success", info = id)
    before <- serialize(fit, NULL, version = 2L)
    p <- sb_model_coefficients(fit, d)
    score <- p$intercept + as.numeric(as.matrix(d[features]) %*% p$table$coefficient)
    prediction <- sb_predict_final(fit, d, "patient_id")
    expect_equal(plogis(score), prediction$predictions$probability, tolerance = 1e-12, info = id)
    nomogram <- sb_nomogram(fit, d)
    expect_equal(sb_nomogram_values(nomogram, d)$score, score, tolerance = 1e-12, info = id)
    shap <- sb_additive_shap(fit, d, "patient_id", d)
    expect_identical(shap$score_scale, "class_1_log_odds")
    expect_equal(unname(rowSums(shap$contributions)) + shap$reference_value, score, tolerance = 1e-12, info = id)
    expect_identical(serialize(fit, NULL, version = 2L), before)
    expect_equal(unname(p$table$scaled_coefficient), unname(p$table$coefficient * vapply(d[features], sd, numeric(1))))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      expect_s3_class(ggplot2::ggplot_build(sb_plot_nomogram(nomogram)), "ggplot_built")
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(sb_linear_coefficients_diagnostic(fit))), "ggplot_built")
    }
  }
})

# ---- liblinear-logistic-explanations ----

test_that("LIBLINEAR logistic explanations match native probabilities in input units", {
  skip_if_not_installed("LiblineaR")
  d <- sb_example_binary(n = 100L, seed = 108L)
  d$gene_c <- d$gene_c * 2 + 4
  d$gene_a <- d$gene_a * 0.3 - 7
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  for (id in c("liblinear_l2_logistic_binary", "liblinear_l1_logistic_binary")) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success", info = id)
    before <- serialize(fit, NULL, version = 2L)
    e <- unserialize(fit$learner_state$engine_raw)
    x <- as.matrix(d[features])
    native_x <- sweep(sweep(x, 2L, e$center, "-"), 2L, e$scale, "/")
    native <- predict(e$model, native_x, proba = TRUE)$probabilities[, "1"]
    parameters <- sb_model_coefficients(fit, d)
    score <- as.numeric(x %*% parameters$table$coefficient + parameters$intercept)
    expect_identical(parameters$table$feature, features)
    expect_equal(plogis(score), unname(native), tolerance = 1e-12)
    expect_equal(sb_predict_final(fit, d, "patient_id")$predictions$probability,
      unname(native), tolerance = 1e-12)
    nomogram <- sb_nomogram(fit, d)
    expect_equal(sb_nomogram_values(nomogram, d)$score, score, tolerance = 1e-12)
    shap <- sb_additive_shap(fit, d, "patient_id", d)
    expect_identical(shap$score_scale, "class_1_log_odds")
    expect_equal(unname(rowSums(shap$contributions)) + shap$reference_value,
      score, tolerance = 1e-12)
    expected_shap <- sweep(sweep(x, 2L, colMeans(x), "-"),
      2L, parameters$table$coefficient, "*")
    expect_equal(unname(shap$contributions), unname(expected_shap), tolerance = 1e-12)
    expect_s3_class(ggplot2::ggplot_build(sb_plot_nomogram(nomogram)), "ggplot_built")
    diagnostic <- sb_linear_coefficients_diagnostic(fit)
    expect_equal(diagnostic$table$coefficient, parameters$table$coefficient)
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(diagnostic)), "ggplot_built")
    expect_identical(serialize(fit, NULL, version = 2L), before)
    # Equivalent native encoding with the opposite first class tests sign handling.
    reversed <- fit
    e$model$ClassNames <- rev(e$model$ClassNames)
    e$model$W <- -e$model$W
    reversed$learner_state$engine_raw <- serialize(e, NULL, version = 2L)
    rp <- sb_model_coefficients(reversed)
    expect_equal(rp$table$coefficient, parameters$table$coefficient, tolerance = 1e-12)
    expect_equal(rp$intercept, parameters$intercept, tolerance = 1e-12)
    expect_equal(unname(predict(e$model, native_x, proba = TRUE)$probabilities[, "1"]),
      unname(native), tolerance = 1e-12)
  }
})

# ---- logistic-adapter ----

test_that("the logistic adapter matches an independent glm fit", {
  d <- make_sb_binary_fixture(n = 60)
  analysis_rows <- 1:48
  assessment_rows <- 49:60
  predictors <- c("gene_a", "gene_b", "gene_c")
  learner <- sb_logistic_learner_adapter()

  model <- learner$fit(
    d[analysis_rows, predictors, drop = FALSE],
    d$outcome[analysis_rows]
  )
  observed <- learner$predict(
    model,
    d[assessment_rows, predictors, drop = FALSE]
  )
  manual <- stats::glm(
    outcome ~ gene_a + gene_b + gene_c,
    data = d[analysis_rows, , drop = FALSE],
    family = stats::binomial(link = "logit")
  )
  expected <- unname(stats::predict(
    manual,
    newdata = d[assessment_rows, , drop = FALSE],
    type = "response"
  ))

  expect_identical(learner$id, "logistic")
  expect_equal(observed, expected, tolerance = 1e-12)
  expect_identical(all(observed >= 0 & observed <= 1), TRUE)

  coefficients <- learner$extract_coefficients(model)
  expect_identical(
    coefficients$term,
    c("(Intercept)", predictors)
  )
  expect_equal(
    coefficients$coefficient,
    unname(stats::coef(manual)),
    tolerance = 1e-12
  )
})

test_that("the logistic adapter locks its numeric predictor schema", {
  d <- make_sb_binary_fixture(n = 40)
  learner <- sb_logistic_learner_adapter()
  model <- learner$fit(d[c("gene_a", "gene_b")], d$outcome)

  error <- capture_sb_error(
    learner$predict(model, d[c("gene_b", "gene_a")])
  )

  expect_s3_class(error, "stormbreaker_learner_schema_mismatch")
})

test_that("rank-deficient logistic fits fail without dropping terms", {
  d <- make_sb_binary_fixture(n = 40)
  x <- data.frame(first = d$gene_a, duplicate = d$gene_a)
  learner <- sb_logistic_learner_adapter()

  error <- capture_sb_error(learner$fit(x, d$outcome))

  expect_s3_class(error, "stormbreaker_invalid_logistic_fit")
})

test_that("the logistic adapter rejects an ambiguous intercept feature", {
  d <- make_sb_binary_fixture(n = 40)
  x <- d[c("gene_a", "gene_b")]
  names(x)[1L] <- "(Intercept)"
  learner <- sb_logistic_learner_adapter()

  error <- capture_sb_error(learner$fit(x, d$outcome))

  expect_s3_class(error, "stormbreaker_reserved_intercept_name")
})

test_that("saturated probabilities alone do not invalidate an overlapping fit", {
  d <- data.frame(
    gene = c(rep(-1, 20), rep(1, 20), -100, 100),
    outcome = c(rep(0, 18), rep(1, 2), rep(0, 2), rep(1, 18), 0, 1)
  )
  native <- suppressWarnings(stats::glm(
    outcome ~ gene, d, family = stats::binomial()
  ))
  model <- sb_fit_logistic(d["gene"], d$outcome)
  expect_equal(stats::coef(model), stats::coef(native), tolerance = 1e-12)
  expect_equal(
    sb_predict_logistic(model, d["gene"]),
    unname(stats::predict(native, d, type = "response")),
    tolerance = 1e-12
  )
  expect_match(
    attr(model, "stormbreaker_fit_warnings"),
    "fitted probabilities numerically 0 or 1 occurred", fixed = TRUE
  )
})

test_that("a probability warning does not admit a separated logistic fit", {
  d <- data.frame(gene = rep(c(-1, 1), each = 20), outcome = rep(0:1, each = 20))
  native <- suppressWarnings(stats::glm(
    outcome ~ gene, d, family = stats::binomial(), x = TRUE, y = TRUE
  ))
  native$converged <- TRUE
  local_mocked_bindings(
    glm = function(...) {
      warning("glm.fit: fitted probabilities numerically 0 or 1 occurred")
      native
    },
    .package = "stats"
  )
  error <- capture_sb_error(sb_fit_logistic(d["gene"], d$outcome))
  expect_s3_class(error, "stormbreaker_logistic_separation")
})

test_that("logistic still rejects unrelated warnings and inconclusive separation", {
  d <- make_sb_binary_fixture(n = 40)
  original_glm <- stats::glm
  message <- "unexpected fitting warning"
  local_mocked_bindings(
    glm = function(...) {
      warning(message)
      original_glm(...)
    },
    .package = "stats"
  )
  error <- capture_sb_error(sb_fit_logistic(d["gene_a"], d$outcome))
  expect_s3_class(error, "stormbreaker_logistic_warning")

  message <- "glm.fit: fitted probabilities numerically 0 or 1 occurred"
  local_mocked_bindings(
    sb_binary_moderate_calibration_separation = function(...) {
      list(status = "failed", detected = NA, error_message = "solver failed")
    }
  )
  error <- capture_sb_error(sb_fit_logistic(d["gene_a"], d$outcome))
  expect_s3_class(error, "stormbreaker_logistic_separation_check_failed")
})

test_that("logistic fit warnings remain in CV and final-fit results", {
  d <- make_sb_binary_fixture(n = 90)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  original_fit <- sb_fit_logistic
  message <- "glm.fit: fitted probabilities numerically 0 or 1 occurred"
  local_mocked_bindings(
    sb_fit_logistic = function(x, outcome) {
      model <- original_fit(x, outcome)
      attr(model, "stormbreaker_fit_warnings") <- message
      model
    }
  )
  cv <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3L, seed = 1L)
  final <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  expect_identical(cv$status, "success")
  expect_identical(final$status, "success")
  expect_identical(
    vapply(
      cv$tuning, function(x) any(grepl("Learner refit:", x$warnings)), logical(1)
    ),
    rep(TRUE, 3L)
  )
  expect_match(tail(final$tuning$warnings, 1L), "Learner refit:")
})

# ---- robust-logistic-explanations ----

test_that("robust logistic explanations preserve the native Mqle model", {
  skip_if_not_installed("robustbase")
  d <- sb_example_binary(n = 160L, seed = 108L)
  d$gene_c <- d$gene_c * 2 + 4
  d$gene_a <- d$gene_a * 0.3 - 7
  features <- c("gene_c", "gene_a", "gene_b")
  for (selected in list(features, "gene_c")) {
    task <- sb_binary_task(d, "patient_id", "outcome", selected)
    fit <- sb_fit_final(task, sb_pipeline(learner = "robust_logistic_mqle"))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    native_data <- d[selected]
    names(native_data) <- paste0(".stormbreaker_robust_x", seq_along(selected))
    native_data$.stormbreaker_robust_outcome <- factor(d$outcome, levels = c(0, 1))
    native <- robustbase::glmrob(
      reformulate(names(native_data)[seq_along(selected)], ".stormbreaker_robust_outcome"),
      family = binomial("logit"), data = native_data, method = "Mqle", weights.on.x = "none",
      control = robustbase::glmrobMqle.control(acc = 1e-4, test.acc = "coef", maxit = 50L, tcc = 3.5),
      model = FALSE, x = FALSE)
    expect_identical(native$converged, TRUE)
    parameters <- sb_model_coefficients(fit, d)
    expect_identical(parameters$table$feature, selected)
    expect_equal(unname(parameters$table$coefficient), unname(coef(native)[-1L]), tolerance = 1e-12)
    expect_equal(parameters$intercept, unname(coef(native)[1L]), tolerance = 1e-12)
    score <- as.numeric(as.matrix(d[selected]) %*% parameters$table$coefficient + parameters$intercept)
    expect_equal(plogis(score), unname(predict(native, native_data, type = "response")), tolerance = 1e-12)
    expect_equal(sb_predict_final(fit, d, "patient_id")$predictions$probability,
      plogis(score), tolerance = 1e-12)
    nomogram <- sb_nomogram(fit, d)
    expect_equal(sb_nomogram_values(nomogram, d)$score, score, tolerance = 1e-12)
    shap <- sb_additive_shap(fit, d, "patient_id", d)
    expected <- sweep(sweep(as.matrix(d[selected]), 2L, colMeans(d[selected]), "-"),
      2L, unname(coef(native)[-1L]), "*")
    expect_equal(unname(shap$contributions), unname(expected), tolerance = 1e-12)
    expect_equal(unname(rowSums(shap$contributions)) + shap$reference_value, score, tolerance = 1e-12)
    expect_identical(shap$score_scale, "class_1_log_odds")
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, "available")
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(bundle$diagnostics[[1L]])), "ggplot_built")
    expect_s3_class(ggplot2::ggplot_build(sb_plot_nomogram(nomogram)), "ggplot_built")
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
})

# ---- model-naive-bayes ----
local({
# ---- algorithm-naive-bayes ----

test_that("Naive Bayes curves preserve fitted densities and posterior decomposition", {
  skip_if_not_installed("klaR")
  d <- sb_example_binary(n = 100L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  for (id in c("gaussian_naive_bayes", "kernel_naive_bayes_binary")) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    q <- sb_naive_bayes_diagnostics(fit)
    expect_null(fit$diagnostic_artifacts)
    expect_identical(q$features, features)
    p <- sb_plot_algorithm_diagnostic(q)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    expect_null(p$labels$caption)
    expect_null(p$labels$subtitle)
    expect_identical(p$facet$params$free, list(x = TRUE, y = TRUE))
    if (id == "gaussian_naive_bayes") {
      m <- fit$learner_state$payload
      for (feature in features) for (class in c("0", "1")) {
        tab <- q$table[q$table$feature == feature & q$table$class == class, ]
        expect_equal(tab$value, dnorm(tab$x, m$means[class, feature], sqrt(m$variances[class, feature])))
      }
      ratio <- sb_naive_bayes_diagnostics(fit, "log_likelihood_ratio")
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(ratio)), "ggplot_built")
      expect_identical(sb_plot_algorithm_diagnostic(ratio)$facet$params$free, list(x = TRUE, y = FALSE))
      # Pick one stored grid point per feature to reconstruct the exact model
      # posterior, including class priors and the original variance smoothing.
      row <- setNames(as.data.frame(as.list(vapply(features, function(f) ratio$table$x[ratio$table$feature == f][40L], numeric(1)))), features)
      terms <- vapply(features, function(f) ratio$table$value[ratio$table$feature == f][40L], numeric(1))
      original <- sb_predict_naive_bayes_learner_state(fit$learner_state, row)
      expect_equal(plogis(ratio$prior_log_odds + sum(terms)), as.numeric(original), tolerance = 1e-12)
    } else {
      e <- unserialize(fit$learner_state$engine_raw)
      for (feature in features) for (class in c("0", "1")) {
        tab <- q$table[q$table$feature == feature & q$table$class == class, ]
        expect_identical(tab$x, e$tables[[feature]][[class]]$x)
        expect_identical(tab$value, e$tables[[feature]][[class]]$y)
      }
      expect_error(sb_naive_bayes_diagnostics(fit, "log_likelihood_ratio"), class = "stormbreaker_unavailable_diagnostic")
    }
    single <- sb_naive_bayes_diagnostics(fit, features = features[2L])
    expect_identical(unique(single$table$feature), features[2L])
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(single)), "ggplot_built")
    expect_error(sb_naive_bayes_diagnostics(fit, features = "missing"), class = "stormbreaker_invalid_diagnostic")
    expect_true(identical(serialize(fit, NULL, version = 2L), before))
  }
})

# ---- centroid-bernoulli ----

test_that("centroid scores match Euclidean geometry including ties", {
  x <- data.frame(a = c(0, 2, 4, 6), constant = 1)
  y <- c(0L, 0L, 1L, 1L)
  id <- "nearest_centroid_binary"
  m <- sb_fit_classical_binary(x, y, id, sb_classical_binary_params(id))
  expect_equal(unname(m$engine$centroids), matrix(c(1, 5, 1, 1), 2))
  q <- data.frame(a = c(1, 3, 5), constant = 1)
  p <- sb_predict_classical_binary(m, q)
  expect_equal(p, plogis(c(-8, 0, 8)))
  saved <- sb_extract_classical_binary_learner_state(m)
  expect_equal(sb_predict_classical_binary_learner_state(saved, q)$probability, p)
})

test_that("Bernoulli NB includes absent features, priors and Laplace smoothing", {
  x <- data.frame(a = c(-2, 0, 1, 2), b = c(0, 1, 0, 1), constant = 0)
  id <- "bernoulli_naive_bayes_binary"
  m <- sb_fit_classical_binary(x, c(0, 0, 0, 1), id,
    sb_classical_binary_params(id))
  expect_equal(m$engine$prior, c(.75, .25))
  expect_equal(unname(m$engine$feature_probability), rbind(c(2/5, 2/5, 1/5), c(2/3, 2/3, 1/3)))
  q <- data.frame(a = 0, b = 1, constant = 0)
  likelihood <- c(.75 * 3/5 * 2/5 * 4/5, .25 * 1/3 * 2/3 * 2/3)
  expect_equal(sb_predict_classical_binary(m, q), likelihood[2]/sum(likelihood))
})

test_that("new learners support CV, selection, saved prediction and external AUC", {
  d <- sb_example_binary(100, 17)
  genes <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", genes)
  external <- sb_example_binary(50, 18)
  for (id in c("nearest_centroid_binary", "bernoulli_naive_bayes_binary")) {
    for (selector in c("none", "univ_binary_auc_top_k", "lasso_binary_nonzero")) {
      p <- sb_pipeline(selector, id,
        selector_n_features = if (selector == "univ_binary_auc_top_k") 2L else NULL,
        selector_options = if (selector == "lasso_binary_nonzero") list(lambda_rule = "lambda.min") else NULL)
      expect_identical(sb_lint_pipeline(p)$status, "admitted")
      fit <- sb_fit_final(task, p)
      expect_identical(fit$status, "success")
      pred <- sb_predict_final(unserialize(serialize(fit, NULL)), external[c("patient_id", genes)], "patient_id")
      ev <- sb_evaluate_external(pred, sb_binary_endpoint(external, "patient_id", "outcome"))
      expect_identical(ev$status, "success")
      if (selector != "lasso_binary_nonzero") {
        cv <- sb_run_cv(task, p, folds = 3, seed = 17)
        expect_identical(cv$status, "success")
        expect_s3_class(sb_plot_roc(cv), "ggplot")
      }
    }
  }
})

test_that("CV learns class statistics only from analysis rows", {
  d <- sb_example_binary(90, 22)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  original <- sb_fit_classical_binary
  seen <- list()
  local_mocked_bindings(sb_fit_classical_binary = function(x, outcome, learner, params) {
    seen[[length(seen) + 1L]] <<- rownames(x)
    original(x, outcome, learner, params)
  })
  for (id in c("nearest_centroid_binary", "bernoulli_naive_bayes_binary")) {
    seen <- list()
    cv <- sb_run_cv(task, sb_pipeline(learner = id), folds = 3, seed = 22)
    expect_identical(cv$status, "success")
    expect_length(seen, 3L)
    for (i in seq_along(seen)) {
      expect_setequal(seen[[i]], rownames(d)[cv$folds[[i]]$analysis_rows])
    }
  }
})

# ---- gaussian-naive-bayes-adapter ----

test_that("Gaussian Naive Bayes matches independent population calculations", {
  x <- data.frame(
    gene_a = c(0, 1, 2, 2, 4, 6),
    gene_b = c(1, 3, 5, 0, 1, 2),
    check.names = FALSE
  )
  outcome <- c(0L, 0L, 0L, 1L, 1L, 1L)
  new_x <- data.frame(
    gene_a = c(1.5, 5),
    gene_b = c(2, 1),
    check.names = FALSE
  )
  var_smoothing <- 0.1
  epsilon <- var_smoothing * max(vapply(x, population_variance, numeric(1)))
  means <- rbind(
    `0` = vapply(x[outcome == 0L, , drop = FALSE], mean, numeric(1)),
    `1` = vapply(x[outcome == 1L, , drop = FALSE], mean, numeric(1))
  )
  variances <- rbind(
    `0` = vapply(
      x[outcome == 0L, , drop = FALSE], population_variance, numeric(1)
    ) + epsilon,
    `1` = vapply(
      x[outcome == 1L, , drop = FALSE], population_variance, numeric(1)
    ) + epsilon
  )
  priors <- c(`0` = 0.5, `1` = 0.5)
  adapter <- sb_gaussian_naive_bayes_learner_adapter(var_smoothing)
  model <- adapter$fit(x, outcome)
  observed <- adapter$predict(model, new_x)
  one_at_a_time <- vapply(seq_len(nrow(new_x)), function(row) {
    adapter$predict(model, new_x[row, , drop = FALSE])
  }, numeric(1))
  reversed <- adapter$predict(model, new_x[2:1, , drop = FALSE])[2:1]
  expected <- manual_gaussian_nb_probability(
    as.matrix(new_x), priors, means, variances
  )
  audit <- adapter$extract_parameters(model)

  expect_s3_class(model, "stormbreaker_gaussian_naive_bayes")
  expect_equal(model$epsilon, epsilon, tolerance = 0)
  expect_equal(model$means, means, tolerance = 0)
  expect_equal(model$variances, variances, tolerance = 0)
  expect_equal(observed, expected, tolerance = 1e-15)
  expect_equal(observed, one_at_a_time, tolerance = 0)
  expect_equal(observed, reversed, tolerance = 0)
  expect_identical(all(observed >= 0 & observed <= 1), TRUE)
  expect_identical(
    audit$priors,
    data.frame(
      class = c(0L, 1L),
      n = c(3L, 3L),
      prior = c(0.5, 0.5),
      stringsAsFactors = FALSE
    )
  )
  expect_identical(audit$var_smoothing, var_smoothing)
  expect_equal(audit$epsilon, epsilon, tolerance = 0)
  expect_identical(
    audit$values$feature,
    rep(names(x), times = 2L)
  )
  expect_identical(audit$values$class, rep(c(0L, 1L), each = ncol(x)))
  expect_equal(
    audit$values$mean,
    unname(c(means["0", ], means["1", ])),
    tolerance = 0
  )
  expect_equal(
    audit$values$variance,
    unname(c(variances["0", ], variances["1", ])),
    tolerance = 0
  )
  expect_equal(
    audit$values$standard_deviation,
    sqrt(audit$values$variance),
    tolerance = 0
  )
})

test_that("Gaussian Naive Bayes locks predictor schema and smoothing", {
  x <- data.frame(
    gene_a = c(0, 1, 2, 3),
    gene_b = c(1, 3, 2, 4)
  )
  outcome <- c(0L, 0L, 1L, 1L)
  adapter <- sb_gaussian_naive_bayes_learner_adapter(1e-9)
  model <- adapter$fit(x, outcome)

  schema_error <- capture_sb_error(
    adapter$predict(model, x[c("gene_b", "gene_a")])
  )
  expect_s3_class(schema_error, "stormbreaker_learner_schema_mismatch")

  invalid_values <- list(-1, NA_real_, Inf, "1e-9", c(0, 1e-9))
  for (value in invalid_values) {
    error <- capture_sb_error(
      sb_gaussian_naive_bayes_learner_adapter(value)
    )
    expect_s3_class(
      error,
      "stormbreaker_bad_naive_bayes_var_smoothing"
    )
  }
})

test_that("Gaussian Naive Bayes predicts a one-patient assessment fold", {
  x <- data.frame(
    gene_a = c(0, 1, 2, 3, 4, 5),
    gene_b = c(1, 3, 5, 0, 2, 4)
  )
  outcome <- rep(c(0L, 1L), each = 3L)
  new_x <- x[4L, , drop = FALSE]
  adapter <- sb_gaussian_naive_bayes_learner_adapter(1e-9)
  model <- adapter$fit(x, outcome)

  probability <- adapter$predict(model, new_x)

  expect_length(probability, 1L)
  expect_true(is.finite(probability))
  expect_true(probability >= 0 && probability <= 1)
})

test_that("Gaussian Naive Bayes exposes invalid fitted variances", {
  all_constant <- data.frame(
    gene_a = rep(1, 6),
    gene_b = rep(2, 6)
  )
  outcome <- rep(c(0L, 1L), each = 3L)
  adapter <- sb_gaussian_naive_bayes_learner_adapter(1e-9)

  constant_error <- capture_sb_error(adapter$fit(all_constant, outcome))
  expect_s3_class(
    constant_error,
    "stormbreaker_invalid_naive_bayes_variance"
  )

  within_class_constant <- data.frame(
    gene_a = c(rep(0, 3), rep(1, 3)),
    gene_b = c(rep(2, 3), rep(3, 3))
  )
  no_smoothing <- sb_gaussian_naive_bayes_learner_adapter(0)
  zero_error <- capture_sb_error(
    no_smoothing$fit(within_class_constant, outcome)
  )
  expect_s3_class(
    zero_error,
    "stormbreaker_invalid_naive_bayes_variance"
  )
})
})

# ---- model-mclust ----
local({
# ---- algorithm-mclust-da ----

test_that("model-based discriminant views preserve full component distributions", {
  skip_if_not_installed("mclust")
  d <- sb_example_binary(180L, seed = 108L)
  d$gene_c <- d$gene_c + rep(c(-5, 5), 90L)
  d$gene_b <- .4 * d$gene_c + d$gene_b
  for (features in list("gene_c", c("gene_c", "gene_a", "gene_b"))) {
    for (id in sb_mclust_da_learner_ids()) {
      fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = id))
      expect_identical(fit$status, "success")
      before <- serialize(fit, NULL, version = 2L)
      e <- fit$learner_state$engine
      centers <- sb_mclust_da_diagnostics(fit)
      correlations <- sb_mclust_da_diagnostics(fit, "correlations")
      weights <- sb_mclust_da_diagnostics(fit, "mixture_weights")
      if (id == "mclustda_binary_bic") expect_gt(max(weights$table$component), 1L)
      q <- d[c(1L, 8L, 13L), ]; q$gene_c <- q$gene_c + .71
      mass <- matrix(0, nrow(q), 2L)
      for (cls in c("0", "1")) {
        par <- e$models[[cls]]$parameters
        w <- weights$table[weights$table$class == cls, ]
        expect_equal(w$weight, as.numeric(par$pro))
        for (j in seq_along(par$pro)) {
          rows <- centers$table$class == cls & centers$table$component == j
          mu <- centers$table$mean[rows]
          expect_identical(centers$table$feature[rows], features)
          expect_equal(mu, as.numeric(matrix(par$mean, nrow = length(features))[, j]))
          covariance <- if (length(features) == 1L) matrix(rep(par$variance$sigmasq, length.out = length(par$pro))[j]) else par$variance$sigma[, , j]
          cor_rows <- correlations$table$class == cls & correlations$table$component == j
          cor <- matrix(correlations$table$correlation[cor_rows], nrow = length(features))
          expect_equal(cor, unname(stats::cov2cor(covariance)))
          sd <- sqrt(diag(covariance))
          recovered_cov <- cor * outer(sd, sd)
          delta <- sweep(as.matrix(q[features]), 2L, mu)
          density <- exp(-rowSums((delta %*% solve(recovered_cov)) * delta) / 2) /
            sqrt((2 * pi)^length(features) * det(recovered_cov))
          k <- match(cls, c("0", "1"))
          mass[, k] <- mass[, k] + density * w$weight[j] * w$class_prior[j]
        }
      }
      expect_equal(mass[, 2L] / rowSums(mass), sb_predict_final(fit, q, "patient_id")$predictions$probability, ignore_attr = TRUE, tolerance = 1e-12)
      for (z in list(centers, correlations, weights)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
      expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status, rep("available", 3L))
      expect_identical(serialize(fit, NULL, version = 2L), before)
    }
  }
})

# ---- mclust-da-adapters ----

make_mclust_da_fixture <- function(n = 160L, seed = 9201L) {
  withr::with_seed(seed, {
    gene_a <- stats::rnorm(n)
    gene_b <- stats::rnorm(n)
    gene_c <- stats::rnorm(n)
    gene_d <- stats::rnorm(n)
    probability <- stats::plogis(
      1.1 * gene_a - 0.7 * gene_b + 0.35 * gene_c
    )
    data.frame(
      patient_id = sprintf("MCLUST%03d", seq_len(n)),
      outcome = stats::rbinom(n, 1L, probability),
      check.names = FALSE,
      "gene A" = gene_a,
      "gene-B" = gene_b,
      "gene.C" = gene_c,
      "gene D" = gene_d
    )
  })
}

native_mclust_da_fit <- function(x, outcome, params) {
  environment <- new.env(parent = asNamespace("mclust"))
  environment$data <- data.matrix(x)
  environment$class <- factor(outcome, levels = c(0L, 1L))
  environment$G <- params$components
  environment$modelNames <- if (ncol(x) == 1L) {
    params$univariate_covariance_models
  } else {
    params$multivariate_covariance_models
  }
  environment$modelType <- params$model_type
  environment$control <- do.call(mclust::emControl, params$em_control)
  eval(
    quote(MclustDA(
      data = data,
      class = class,
      G = G,
      modelNames = modelNames,
      modelType = modelType,
      prior = NULL,
      control = control,
      warn = FALSE,
      verbose = FALSE
    )),
    envir = environment
  )
}

test_that("mclust DA learners declare one analysis-fold BIC specification", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")

  expect_identical(
    sb_mclust_da_learner_ids(),
    c("edda_binary_bic", "mclustda_binary_bic")
  )
  edda <- sb_mclust_da_params("edda_binary_bic")
  mixture <- sb_mclust_da_params("mclustda_binary_bic")

  expect_identical(edda$model_type, "EDDA")
  expect_identical(edda$components, 1L)
  expect_identical(mixture$model_type, "MclustDA")
  expect_identical(mixture$components, 1:3)
  expect_identical(edda$univariate_covariance_models, c("E", "V"))
  expect_identical(
    edda$multivariate_covariance_models,
    mclust::mclust.options("emModelNames")
  )
  expect_identical(edda$class_prior, "analysis_fold_empirical")
  expect_identical(edda$parameter_prior, "none")
  expect_identical(edda$scale, "original")
  expect_identical(edda$criterion, "BIC")
  expect_identical(edda$native_candidate_warnings, FALSE)
  expect_identical(edda$all_model_failure, "typed")
  expect_identical(
    edda$initialization,
    list(method = "known_class_labels")
  )
  expect_identical(
    mixture$initialization,
    list(
      method = "native_hierarchical",
      subset_max = 2000L,
      hc_model_name = "VVV",
      hc_use = "SVD",
      seed = 1L
    )
  )

  expect_identical(
    sb_pipeline(learner = "edda_binary_bic")$learner_params,
    edda
  )
  expect_identical(
    sb_pipeline(learner = "mclustda_binary_bic")$learner_params,
    mixture
  )
})

test_that("mclust DA probabilities match the native posterior", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  data <- make_mclust_da_fixture()
  features <- c("gene A", "gene-B", "gene.C", "gene D")
  x <- data[features]
  new_x <- x[1:13, , drop = FALSE]

  for (learner in sb_mclust_da_learner_ids()) {
    params <- sb_mclust_da_params(learner)
    adapter <- sb_mclust_da_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    native <- native_mclust_da_fit(x, data$outcome, params)
    expected <- unname(stats::predict(
      native,
      newdata = data.matrix(new_x)
    )$z[, "1"])

    expect_equal(adapter$predict(model, new_x), expected, tolerance = 1e-12)
    expect_identical(model$audit$class_prior, native$prop)
    expect_identical(model$audit$selected$class, c("0", "1"))
  }
})

test_that("mclust DA uses the native univariate model grid", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  data <- make_mclust_da_fixture(seed = 9202L)
  x <- data["gene A"]

  for (learner in sb_mclust_da_learner_ids()) {
    adapter <- sb_mclust_da_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    probability <- adapter$predict(model, x[1:9, , drop = FALSE])

    expect_length(probability, 9L)
    expect_equal(range(probability), pmin(pmax(range(probability), 0), 1))
    expect_identical(model$params$active_covariance_models, c("E", "V"))
  }
})

test_that("mclust DA rejects a native silent reduction of the BIC grid", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  data <- make_mclust_da_fixture(seed = 9203L)
  x <- data[c("gene A", "gene-B", "gene.C", "gene D")]
  params <- sb_mclust_da_params("mclustda_binary_bic")
  engine <- native_mclust_da_fit(x, data$outcome, params)
  engine$models[[1L]]$BIC <- engine$models[[1L]]$BIC[-3L, , drop = FALSE]

  error <- tryCatch(
    {
      sb_validate_mclust_da_search_grid(engine, params)
      NULL
    },
    error = identity
  )
  expect_s3_class(error, "stormbreaker_mclust_da_search_reduction")
})

test_that("EDDA rejects its real dimension-based covariance-grid reduction", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  x <- withr::with_seed(92031L, {
    as.data.frame(matrix(
      stats::rnorm(20L * 10L),
      nrow = 20L,
      ncol = 10L
    ))
  })
  names(x) <- sprintf("gene_%02d", seq_len(ncol(x)))
  outcome <- rep(0:1, each = 10L)

  error <- tryCatch(
    {
      sb_mclust_da_learner_adapter("edda_binary_bic")$fit(x, outcome)
      NULL
    },
    error = identity
  )

  expect_s3_class(error, "stormbreaker_mclust_da_search_reduction")
})

test_that("EDDA keeps the full grid when one class exceeds the dimension", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  outcome <- c(rep(0L, 10L), rep(1L, 20L))
  x <- withr::with_seed(92032L, {
    matrix_x <- matrix(
      stats::rnorm(30L * 15L),
      nrow = 30L,
      ncol = 15L
    )
    matrix_x[, 1L] <- matrix_x[, 1L] + 0.6 * outcome
    as.data.frame(matrix_x)
  })
  names(x) <- sprintf("gene_%02d", seq_len(ncol(x)))

  model <- sb_mclust_da_learner_adapter("edda_binary_bic")$fit(x, outcome)

  expect_s3_class(model, "stormbreaker_mclust_da")
  expect_identical(
    model$audit$candidate_summary$requested,
    c(14L, 14L)
  )
})

test_that("mclust DA portable states omit analysis rows and replay", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  data <- make_mclust_da_fixture(seed = 9204L)
  features <- c("gene A", "gene-B", "gene.C", "gene D")
  x <- data[features]
  new_x <- x[1:11, , drop = FALSE]

  for (learner in sb_mclust_da_learner_ids()) {
    adapter <- sb_mclust_da_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    state <- unserialize(serialize(adapter$extract_state(model, list()), NULL))

    expect_equal(
      adapter$predict_state(state, new_x)$probability,
      adapter$predict(model, new_x),
      tolerance = 1e-12
    )
    expect_setequal(names(state$engine), c("type", "prop", "models"))
    for (class_model in state$engine$models) {
      expect_setequal(names(class_model), c("modelName", "parameters"))
    }
    expect_null(state$data)
    expect_null(state$class)
    expect_null(state$outcome)
    expect_null(state$engine$data)
    expect_null(state$engine$class)
    expect_null(state$engine$call)
  }
})

test_that("mclust DA does not alter the caller RNG state", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  data <- make_mclust_da_fixture(seed = 9205L)
  x <- data[c("gene A", "gene-B", "gene.C", "gene D")]

  set.seed(9206L)
  before <- .Random.seed
  for (learner in sb_mclust_da_learner_ids()) {
    sb_mclust_da_learner_adapter(learner)$fit(x, data$outcome)
  }
  expect_identical(.Random.seed, before)
})

test_that("MclustDA contains and audits sampled hierarchical initialization", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  if (!"package:mclust" %in% search()) {
    suppressWarnings(withr::local_package("mclust"))
  }
  original_options <- list(
    subset = mclust::mclust.options("subset"),
    hcModelName = mclust::mclust.options("hcModelName"),
    hcUse = mclust::mclust.options("hcUse")
  )
  withr::defer(do.call(mclust::mclust.options, original_options))

  n <- 240L
  outcome <- rep(0:1, each = n / 2L)
  x <- withr::with_seed(92061L, {
    data.frame(
      gene_a = stats::rnorm(n) + 0.8 * outcome,
      gene_b = stats::rnorm(n) - 0.5 * outcome,
      gene_c = stats::rnorm(n),
      gene_d = stats::rnorm(n)
    )
  })
  params <- sb_mclust_da_params("mclustda_binary_bic")
  params$initialization$subset_max <- 50L
  adapter <- sb_mclust_da_learner_adapter(
    "mclustda_binary_bic",
    params = params
  )
  new_x <- x[1:15, , drop = FALSE]

  mclust::mclust.options(
    subset = 40,
    hcModelName = "EII",
    hcUse = "VARS"
  )
  set.seed(92062L)
  rng_before_first <- .Random.seed
  first <- adapter$fit(x, outcome)
  first_probability <- adapter$predict(first, new_x)

  expect_identical(.Random.seed, rng_before_first)
  expect_identical(mclust::mclust.options("subset"), 40)
  expect_identical(mclust::mclust.options("hcModelName"), "EII")
  expect_identical(mclust::mclust.options("hcUse"), "VARS")
  expect_identical(
    first$audit$candidate_summary$initialization_subset_size,
    c(50L, 50L)
  )

  mclust::mclust.options(
    subset = 30,
    hcModelName = "VVV",
    hcUse = "SVD"
  )
  set.seed(92063L)
  rng_before_second <- .Random.seed
  second <- adapter$fit(x, outcome)
  second_probability <- adapter$predict(second, new_x)

  expect_identical(.Random.seed, rng_before_second)
  expect_identical(mclust::mclust.options("subset"), 30)
  expect_identical(mclust::mclust.options("hcModelName"), "VVV")
  expect_identical(mclust::mclust.options("hcUse"), "SVD")
  expect_identical(
    second$audit$candidate_summary$initialization_subset_size,
    c(50L, 50L)
  )
  expect_identical(first$audit$selected, second$audit$selected)
  expect_equal(first_probability, second_probability, tolerance = 1e-12)
})

test_that("mclust DA learners run through CV and outcome-blind prediction", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  data <- make_mclust_da_fixture(n = 180L, seed = 9207L)
  features <- c("gene A", "gene-B", "gene.C", "gene D")
  task <- sb_binary_task(data, "patient_id", "outcome", features)
  new_data <- data[1:9, c("patient_id", features), drop = FALSE]

  for (learner in sb_mclust_da_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 3L
    )
    cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 9208L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(fit, new_data, "patient_id")

    expect_identical(cv$status, "success")
    expect_identical(cv$complete_oof, TRUE)
    expect_identical(length(cv$selections), 3L)
    expect_identical(
      vapply(cv$selections, function(selection) {
        length(selection$selected_features)
      }, integer(1)),
      rep.int(3L, 3L)
    )
    expect_identical(
      unique(vapply(cv$coefficients, `[[`, "", "status")),
      "not_applicable"
    )
    expect_identical(
      cv$provenance$mclust_version,
      as.character(utils::packageVersion("mclust"))
    )
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_identical(
      fit$provenance$mclust_version,
      as.character(utils::packageVersion("mclust"))
    )
    expect_equal(nrow(prediction$predictions), nrow(new_data))
    expect_equal(
      range(prediction$predictions$probability),
      pmin(pmax(range(prediction$predictions$probability), 0), 1)
    )
  }
})

test_that("mclust DA learners retain the supported selector routes", {
  registry <- sb_pipeline_registry()
  rows <- registry[
    registry$learner %in% c("edda_binary_bic", "mclustda_binary_bic"), ,
    drop = FALSE
  ]

  expected_selectors <- expected_fixed_learner_selectors("binary")
  expect_identical(nrow(rows), 2L * length(expected_selectors))
  for (learner in c("edda_binary_bic", "mclustda_binary_bic")) {
    expect_setequal(rows$selector[rows$learner == learner], expected_selectors)
  }
  expect_identical(unique(rows$tuning_semantics), "analysis_fold_bic")
  expect_identical(
    all(grepl("mclust", rows$optional_dependency, fixed = TRUE)),
    TRUE
  )
})
})

# ---- model-linear-margin ----
local({
# ---- liblinear-formulations ----

make_liblinear_formulation_fixture <- function(n = 72L,
                                               p = 5L,
                                               seed = 9311L) {
  withr::with_seed(seed, {
    outcome <- sample(rep(c(0L, 1L), length.out = n))
    x <- matrix(stats::rnorm(n * p), nrow = n)
    x[, 1L] <- x[, 1L] + 1.5 * outcome
    if (p >= 2L) {
      x[, 2L] <- x[, 2L] - 1.1 * outcome
    }
    names <- c(
      "gene one",
      if (p >= 2L) "2gene" else NULL,
      if (p >= 3L) paste0("gene-", seq_len(p - 2L)) else NULL
    )
    colnames(x) <- names
    data.frame(
      patient_id = sprintf("LL%03d", seq_len(n)),
      outcome = outcome,
      x,
      check.names = FALSE
    )
  })
}

test_that("LIBLINEAR formulations declare four distinct new objectives", {
  expected <- c(
    "liblinear_l2_logistic_binary",
    "liblinear_l1_logistic_binary",
    "liblinear_l2_l2loss_svc_binary",
    "liblinear_l2_l1loss_svc_binary",
    "liblinear_l1_l2loss_svc_binary"
  )
  expect_identical(sb_liblinear_binary_learner_ids(), expected)

  expected_types <- c(0L, 6L, 1L, 3L, 5L)
  expected_tolerances <- c(0.01, 0.01, 0.1, 0.1, 0.01)
  expected_probabilities <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
  for (i in seq_along(expected)) {
    params <- sb_classical_binary_params(expected[[i]])
    expect_identical(params$solver_type, expected_types[[i]])
    expect_identical(params$tolerance, expected_tolerances[[i]])
    expect_identical(params$probability, expected_probabilities[[i]])
    expect_identical(params$cost, 1)
    expect_identical(params$bias, 1)
    expect_identical(params$standardize, TRUE)
  }
  for (learner in expected[-1L]) {
    expect_identical(
      sb_classical_binary_params(learner)$seed,
      1L
    )
  }

  registry <- sb_pipeline_registry()
  expected_selectors <- expected_fixed_learner_selectors("binary")
  rows <- registry[
    registry$learner %in% expected[-1L], ,
    drop = FALSE
  ]
  expect_identical(
    nrow(rows),
    length(expected[-1L]) * length(expected_selectors)
  )
  expect_identical(
    sort(unique(rows$selector)),
    sort(expected_selectors)
  )
  expect_true(all(grepl("LiblineaR", rows$optional_dependency, fixed = TRUE)))
})

test_that("LIBLINEAR formulations reproduce native probability or margin", {
  skip_if_not_installed("LiblineaR", minimum_version = "2.10.0")

  data <- make_liblinear_formulation_fixture()
  features <- setdiff(names(data), c("patient_id", "outcome"))
  x <- data[features]
  assessment <- x[c(2L, 19L, 51L), , drop = FALSE]

  for (learner in sb_liblinear_binary_learner_ids()[-1L]) {
    adapter <- sb_classical_binary_learner_adapter(learner)

    set.seed(9312L)
    caller_rng <- .Random.seed
    model <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, caller_rng)

    set.seed(9313L)
    second_caller_rng <- .Random.seed
    repeated <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, second_caller_rng)
    expect_identical(
      serialize(model$engine, NULL),
      serialize(repeated$engine, NULL)
    )

    standardized <- sweep(
      sweep(as.matrix(assessment), 2L, model$engine$center, "-"),
      2L,
      model$engine$scale,
      "/"
    )
    native <- if (isTRUE(adapter$params$probability)) {
      stats::predict(
        model$engine$model,
        standardized,
        proba = TRUE,
        decisionValues = FALSE
      )$probabilities[, "1"]
    } else {
      prediction <- stats::predict(
        model$engine$model,
        standardized,
        proba = FALSE,
        decisionValues = TRUE
      )
      decision <- prediction$decisionValues
      mapped <- stats::plogis(decision[, "1"] - decision[, "0"])
      expect_identical(
        as.integer(mapped >= 0.5),
        as.integer(prediction$predictions)
      )
      mapped
    }
    observed <- adapter$predict(model, assessment)
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, assessment)$probability

    expect_identical(model$engine$model$Type, adapter$params$solver_type)
    expect_equal(observed, unname(native), tolerance = 0)
    expect_equal(replay, observed, tolerance = 0)
    expect_type(state$engine_raw, "raw")
    expect_null(state$training_data_retention)
  }
})

test_that("LIBLINEAR formulations support p greater than n and one feature", {
  skip_if_not_installed("LiblineaR", minimum_version = "2.10.0")

  cases <- list(
    high_dimensional = make_liblinear_formulation_fixture(
      n = 40L,
      p = 48L,
      seed = 9314L
    ),
    one_feature = make_liblinear_formulation_fixture(
      n = 48L,
      p = 1L,
      seed = 9315L
    )
  )
  for (data in cases) {
    features <- setdiff(names(data), c("patient_id", "outcome"))
    x <- data[features]
    for (learner in sb_liblinear_binary_learner_ids()[-1L]) {
      adapter <- sb_classical_binary_learner_adapter(learner)
      model <- adapter$fit(x, data$outcome)
      probability <- adapter$predict(model, x[1:5, , drop = FALSE])

      expect_length(probability, 5L)
      expect_true(all(is.finite(probability)))
      expect_true(all(probability >= 0 & probability <= 1))
    }
  }
})

test_that("LIBLINEAR formulations execute through CV and final prediction", {
  skip_if_not_installed("LiblineaR", minimum_version = "2.10.0")

  data <- make_liblinear_formulation_fixture(n = 72L, seed = 9316L)
  features <- setdiff(names(data), c("patient_id", "outcome"))
  task <- sb_binary_task(data, "patient_id", "outcome", features)

  for (learner in sb_liblinear_binary_learner_ids()[-1L]) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 3L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 9317L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      data[1:6, c("patient_id", rev(features)), drop = FALSE],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_identical(fit$learner_state$learner, learner)
  }
})

# ---- random-ferns-robust-logistic-lssvm ----

test_that("random ferns, robust logistic, and LS-SVM are distinct learners", {
  learners <- c(
    "random_ferns_binary",
    "robust_logistic_mqle",
    "lssvm_radial_binary"
  )

  expect_true(all(learners %in% sb_classical_binary_learner_ids()))
  expect_identical(
    sb_classical_binary_params("random_ferns_binary"),
    list(
      depth = 5L,
      num_ferns = 500L,
      importance = "none",
      save_forest = TRUE,
      threads = 1L,
      seed = 1L,
      score_mapping = "binary_softmax"
    )
  )
  expect_identical(
    sb_classical_binary_params("robust_logistic_mqle"),
    list(
      method = "Mqle",
      weights_on_x = "none",
      tuning_constant = 3.5,
      convergence_tolerance = 1e-4,
      convergence_test = "coef",
      max_iterations = 50L,
      family = "binomial_logit"
    )
  )
  expect_identical(
    sb_classical_binary_params("lssvm_radial_binary"),
    list(
      kernel = "rbfdot",
      kernel_parameters = "automatic",
      sigma_seed = 1L,
      tau = 0.01,
      scaled = TRUE,
      reduced = FALSE,
      cross_validation = 0L,
      keep_fitted = FALSE,
      score_mapping = "binary_softmax"
    )
  )

  dependencies <- vapply(
    learners,
    function(learner) {
      sb_classical_binary_dependency(learner)$package
    },
    character(1)
  )
  expect_identical(
    unname(dependencies),
    c("rFerns", "robustbase", "kernlab")
  )

  registry <- sb_pipeline_registry()
  expected_selectors <- expected_fixed_learner_selectors("binary")
  for (learner in learners) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(unique(rows$task_type), "binary")
    expect_setequal(rows$selector, expected_selectors)
  }
  expect_identical(
    unique(registry$tuning_semantics[
      registry$learner == "lssvm_radial_binary"
    ]),
    "analysis_fold_predictor_only_sigest"
  )
})

test_that("new binary learners replay their declared native predictions", {
  skip_if_not_installed("rFerns", minimum_version = "2.0.0")
  skip_if_not_installed("robustbase", minimum_version = "0.99.0")
  skip_if_not_installed("kernlab", minimum_version = "0.9.0")

  data <- make_sb_binary_fixture(n = 84L, seed = 8811L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene one", "gene-two", "gene_c")
  assessment <- x[c(3L, 19L, 47L, 72L), , drop = FALSE]

  learners <- c(
    "random_ferns_binary",
    "robust_logistic_mqle",
    "lssvm_radial_binary"
  )
  for (learner in learners) {
    params <- sb_classical_binary_params(learner)
    if (identical(learner, "random_ferns_binary")) {
      params$num_ferns <- 120L
    }
    adapter <- sb_classical_binary_learner_adapter(learner, params)
    model <- adapter$fit(x, data$outcome)
    observed <- adapter$predict(model, assessment)
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, assessment)$probability

    native <- switch(learner,
      random_ferns_binary = {
        fit <- withr::with_seed(
          params$seed,
          rFerns::rFerns(
            x = x,
            y = factor(data$outcome, levels = c(0L, 1L)),
            depth = params$depth,
            ferns = params$num_ferns,
            importance = params$importance,
            saveForest = params$save_forest,
            threads = params$threads
          )
        )
        score <- as.matrix(stats::predict(
          fit,
          x = assessment,
          scores = TRUE
        ))
        stats::plogis(score[, "1"] - score[, "0"])
      },
      robust_logistic_mqle = {
        internal_names <- paste0(
          ".stormbreaker_robust_x",
          seq_len(ncol(x))
        )
        analysis <- x
        names(analysis) <- internal_names
        analysis$.stormbreaker_robust_outcome <- factor(
          data$outcome,
          levels = c(0L, 1L)
        )
        fit <- robustbase::glmrob(
          stats::reformulate(
            internal_names,
            response = ".stormbreaker_robust_outcome",
            env = baseenv()
          ),
          data = analysis,
          family = stats::binomial(link = "logit"),
          method = params$method,
          weights.on.x = params$weights_on_x,
          control = robustbase::glmrobMqle.control(
            acc = params$convergence_tolerance,
            test.acc = params$convergence_test,
            maxit = params$max_iterations,
            tcc = params$tuning_constant
          ),
          model = FALSE,
          x = FALSE
        )
        native_new <- assessment
        names(native_new) <- internal_names
        unname(stats::predict(fit, native_new, type = "response"))
      },
      lssvm_radial_binary = {
        fit <- withr::with_seed(params$sigma_seed, {
          invisible(utils::capture.output(
            native_fit <- kernlab::lssvm(
              x = as.matrix(x),
              y = factor(data$outcome, levels = c(0L, 1L)),
              scaled = params$scaled,
              kernel = params$kernel,
              kpar = params$kernel_parameters,
              type = "classification",
              tau = params$tau,
              reduced = params$reduced,
              cross = params$cross_validation,
              fit = params$keep_fitted
            )
          ))
          native_fit
        })
        score <- kernlab::predict(
          fit,
          as.matrix(assessment),
          type = "decision"
        )
        stats::plogis(score[, 2L] - score[, 1L])
      }
    )

    expect_equal(observed, unname(as.numeric(native)), tolerance = 1e-12)
    expect_equal(replay, observed, tolerance = 0)
    expect_true(all(observed > 0 & observed < 1))
  }

  fern_params <- sb_classical_binary_params(
    "random_ferns_binary"
  )
  fern_params$num_ferns <- 120L
  fern_adapter <- sb_classical_binary_learner_adapter(
    "random_ferns_binary",
    fern_params
  )
  fern_state <- fern_adapter$extract_state(
    fern_adapter$fit(x, data$outcome)
  )
  fern_engine <- unserialize(fern_state$engine_raw)
  expect_setequal(
    names(fern_engine),
    c("model", "classLabels", "isStruct", "parameters", "type")
  )
  expect_false("training_data_retention" %in% names(fern_state))

  robust_adapter <- sb_classical_binary_learner_adapter(
    "robust_logistic_mqle"
  )
  robust_state <- robust_adapter$extract_state(
    robust_adapter$fit(x, data$outcome)
  )
  robust_engine <- unserialize(robust_state$engine_raw)
  expect_setequal(
    names(robust_engine),
    c("intercept", "coefficient", "iterations")
  )
  expect_false("training_data_retention" %in% names(robust_state))

  lssvm_adapter <- sb_classical_binary_learner_adapter(
    "lssvm_radial_binary"
  )
  lssvm_state <- lssvm_adapter$extract_state(
    lssvm_adapter$fit(x, data$outcome)
  )
  lssvm_engine <- unserialize(lssvm_state$engine_raw)
  expect_identical(
    lssvm_state$training_data_retention,
    "required_scaled_analysis_predictors_for_radial_kernel_prediction"
  )
  expect_identical(nrow(lssvm_engine$scaled_training), nrow(x))
  expect_false("outcome" %in% names(lssvm_engine))
})

test_that("fixed Random Ferns isolates its deterministic seed", {
  skip_if_not_installed("rFerns", minimum_version = "2.0.0")

  data <- make_sb_binary_fixture(n = 72L, seed = 8812L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  params <- sb_classical_binary_params("random_ferns_binary")
  params$num_ferns <- 100L
  adapter <- sb_classical_binary_learner_adapter(
    "random_ferns_binary",
    params
  )
  withr::local_seed(8813L)
  rng_before <- .Random.seed
  first <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, rng_before)
  second <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, rng_before)
  expect_identical(first$engine$model, second$engine$model)
  expect_equal(
    adapter$predict(first, x[1:8, , drop = FALSE]),
    adapter$predict(second, x[1:8, , drop = FALSE]),
    tolerance = 0
  )
})

test_that("fixed radial LS-SVM isolates its predictor-only sigest seed", {
  skip_if_not_installed("kernlab", minimum_version = "0.9.0")

  data <- make_sb_binary_fixture(n = 72L, seed = 8814L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_classical_binary_learner_adapter(
    "lssvm_radial_binary"
  )
  withr::local_seed(8815L)
  rng_before <- .Random.seed
  first <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, rng_before)
  second <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, rng_before)
  expect_equal(first$engine$sigma, second$engine$sigma, tolerance = 0)
  expect_gt(first$engine$sigma, 0)
  expect_equal(
    adapter$predict(first, x[1:8, , drop = FALSE]),
    adapter$predict(second, x[1:8, , drop = FALSE]),
    tolerance = 0
  )
})

test_that("new binary learners execute through fold-local CV and final fit", {
  skip_if_not_installed("rFerns", minimum_version = "2.0.0")
  skip_if_not_installed("robustbase", minimum_version = "0.99.0")
  skip_if_not_installed("kernlab", minimum_version = "0.9.0")

  data <- make_sb_binary_fixture(n = 72L, seed = 8821L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  dependencies <- c(
    random_ferns_binary = "rferns_version",
    robust_logistic_mqle = "robustbase_version",
    lssvm_radial_binary = "kernlab_version"
  )
  for (learner in names(dependencies)) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 2L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 8822L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:6, c("patient_id", "gene_c", "gene_a", "gene_b")],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_true(all(comparison$runs$status == "success"))
    expect_true(unname(dependencies[[learner]]) %in% names(
      comparison$provenance
    ))
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$probability, 6L)
  }
})

# ---- specialized-binary-linear-adapters ----

test_that("specialized binary linear specifications are explicit", {
  expect_identical(
    sb_specialized_binary_linear_learner_ids(),
    c(
      "stepglm_binary",
      "glmboost_binary",
      "plsrglm_binary"
    )
  )
  expect_identical(
    sb_specialized_binary_linear_params(
      "stepglm_binary",
      stepglm_direction = "forward"
    ),
    list(
      direction = "forward",
      criterion = "AIC",
      k = 2,
      max_steps = 1000L
    )
  )
  expect_identical(
    sb_specialized_binary_linear_params("glmboost_binary"),
    list(
      iterations = 100L,
      shrinkage = 0.1,
      center = FALSE,
      family = "binomial_logit"
    )
  )
  expect_identical(
    sb_specialized_binary_linear_params("plsrglm_binary"),
    list(
      components = 2L,
      scale_x = TRUE,
      model = "pls-glm-logistic"
    )
  )
})

test_that("Stepglm direction is part of the pipeline specification", {
  default <- sb_pipeline(learner = "stepglm_binary")
  backward <- sb_pipeline(
    learner = "stepglm_binary",
    stepglm_direction = "backward"
  )
  bad <- capture_sb_error(sb_pipeline(
    learner = "stepglm_binary",
    stepglm_direction = "sideways"
  ))
  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    stepglm_direction = "forward"
  ))

  expect_identical(default$learner_params$direction, "forward")
  expect_identical(
    default$id,
    "none__stepglm_binary__directionforward"
  )
  expect_identical(
    backward$id,
    "none__stepglm_binary__directionbackward"
  )
  expect_s3_class(bad, "stormbreaker_bad_stepglm_direction")
  expect_s3_class(
    irrelevant,
    "stormbreaker_irrelevant_stepglm_parameter"
  )
})

test_that("Stepglm normalizes equivalent character direction objects", {
  data <- sb_example_binary(n = 72L, seed = 8001L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  directions <- list(
    "forward",
    c(requested = "forward"),
    structure("forward", class = "declared_direction")
  )
  fits <- lapply(directions, function(direction) {
    params <- sb_specialized_binary_linear_params(
      "stepglm_binary",
      stepglm_direction = direction
    )
    adapter <- sb_specialized_binary_linear_learner_adapter(
      "stepglm_binary",
      params
    )
    model <- adapter$fit(x, data$outcome)
    list(
      params = params,
      model = model,
      probability = adapter$predict(model, x)
    )
  })

  expect_identical(
    lapply(fits, function(fit) fit$params$direction),
    rep(list("forward"), 3L)
  )
  expect_identical(
    lapply(fits, function(fit) fit$model$effective$selected_features),
    rep(list(fits[[1L]]$model$effective$selected_features), 3L)
  )
  expect_identical(
    lapply(fits, function(fit) fit$probability),
    rep(list(fits[[1L]]$probability), 3L)
  )
})

test_that("specialized binary linear adapters preserve fitted predictions", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("plsRglm", minimum_version = "1.6.0")

  data <- sb_example_binary(n = 72L, seed = 8001L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(4L, 29L, 61L), , drop = FALSE]

  for (direction in c("forward", "backward", "both")) {
    params <- sb_specialized_binary_linear_params(
      "stepglm_binary",
      stepglm_direction = direction
    )
    adapter <- sb_specialized_binary_linear_learner_adapter(
      "stepglm_binary",
      params
    )
    model <- adapter$fit(x, data$outcome)
    probability <- adapter$predict(model, new_x)
    replay <- adapter$predict_state(adapter$extract_state(model), new_x)

    expect_identical(model$params$direction, direction)
    expect_equal(replay$probability, probability, tolerance = 0)
  }

  glmboost_params <- sb_specialized_binary_linear_params(
    "glmboost_binary"
  )
  glmboost_adapter <- sb_specialized_binary_linear_learner_adapter(
    "glmboost_binary",
    glmboost_params
  )
  glmboost_model <- glmboost_adapter$fit(x, data$outcome)
  glmboost_native <- mboost::glmboost(
    x = as.matrix(x),
    y = factor(data$outcome, levels = c(0L, 1L)),
    family = mboost::Binomial(type = "glm", link = "logit"),
    control = mboost::boost_control(
      mstop = glmboost_params$iterations,
      nu = glmboost_params$shrinkage,
      risk = "inbag",
      stopintern = FALSE,
      trace = FALSE
    ),
    center = glmboost_params$center
  )
  expect_equal(
    glmboost_adapter$predict(glmboost_model, new_x),
    as.numeric(stats::predict(
      glmboost_native,
      newdata = as.matrix(new_x),
      type = "response"
    )),
    tolerance = 1e-14
  )

  plsrglm_params <- sb_specialized_binary_linear_params(
    "plsrglm_binary"
  )
  plsrglm_adapter <- sb_specialized_binary_linear_learner_adapter(
    "plsrglm_binary",
    plsrglm_params
  )
  plsrglm_model <- plsrglm_adapter$fit(x, data$outcome)
  plsrglm_native <- plsRglm::PLS_glm(
    dataY = data$outcome,
    dataX = x,
    nt = plsrglm_params$components,
    modele = plsrglm_params$model,
    typeVC = "none",
    scaleX = plsrglm_params$scale_x,
    pvals.expli = FALSE,
    MClassed = FALSE,
    sparse = FALSE,
    naive = FALSE,
    verbose = FALSE
  )
  native_coefficients <- drop(plsrglm_native$Coeffs)
  native_probability <- stats::plogis(
    unname(native_coefficients[["Intercept"]]) +
      as.vector(
        as.matrix(new_x) %*% native_coefficients[names(new_x)]
      )
  )
  expect_equal(
    plsrglm_adapter$predict(plsrglm_model, new_x),
    native_probability,
    tolerance = 1e-14
  )

  for (model_adapter in list(
    list(model = glmboost_model, adapter = glmboost_adapter),
    list(model = plsrglm_model, adapter = plsrglm_adapter)
  )) {
    state <- model_adapter$adapter$extract_state(model_adapter$model)
    replay <- model_adapter$adapter$predict_state(
      unserialize(serialize(state, NULL)),
      new_x
    )
    direct <- model_adapter$adapter$predict(model_adapter$model, new_x)

    expect_identical("engine_raw" %in% names(state), FALSE)
    expect_equal(replay$probability, direct, tolerance = 0)
  }
})

test_that("specialized binary linear learners run CV and final prediction", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("plsRglm", minimum_version = "1.6.0")

  data <- sb_example_binary(n = 72L, seed = 8002L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    predictors = c("gene_a", "gene_b", "gene_c")
  )
  learners <- sb_specialized_binary_linear_learner_ids()

  for (learner in learners) {
    pipeline <- if (identical(learner, "stepglm_binary")) {
      sb_pipeline(learner = learner, stepglm_direction = "both")
    } else {
      sb_pipeline(learner = learner)
    }
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 8003L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:6, c("patient_id", "gene_c", "gene_a", "gene_b")],
      "patient_id"
    )

    expect_identical(result$status, "success")
    expect_identical(result$complete_oof, TRUE)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
  }
})

test_that("specialized binary linear learners use legal feature selectors", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")

  data <- sb_example_binary(n = 72L, seed = 8004L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    predictors = c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "stepglm_binary",
      selector_n_features = 2L,
      stepglm_direction = "forward"
    ),
    sb_pipeline(
      selector = "fixed_features",
      learner = "glmboost_binary",
      selector_features = c("gene_b", "gene_a"),
      feature_set_id = "two_genes"
    )
  )

  for (pipeline in pipelines) {
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 8005L)
    fit <- sb_fit_final(task, pipeline)

    expect_identical(result$status, "success")
    expect_identical(fit$status, "success")
    expect_length(fit$required_features, 2L)
  }
})
})

# ---- model-svm ----
local({
# ---- survivalsvm-additive-kernel ----

test_that("additive kernel query batching preserves feature and sample ordering", {
  basis <- matrix(c(-2, 0, 3, 1, 4, 7), nrow = 3L)
  colnames(basis) <- c("gene A", "gene-B")
  x <- matrix(c(-4, 2, 6, 0, 5, 9), nrow = 3L)
  colnames(x) <- colnames(basis)
  for (columns in list(1L, 1:2)) {
    b <- basis[, columns, drop = FALSE]
    queries <- x[, columns, drop = FALSE]
    result <- sb_survivalsvm_additive_kernel_matrix(b, queries)
    individual <- do.call(cbind, lapply(seq_len(nrow(queries)), function(i) {
      sb_survivalsvm_additive_kernel_matrix(b, queries[i, , drop = FALSE])
    }))
    expect_identical(result, individual)
    expect_identical(
      sb_survivalsvm_additive_kernel_matrix(b, queries[3:1, , drop = FALSE]),
      result[, 3:1, drop = FALSE]
    )
  }
})

make_survivalsvm_additive_fixture <- function(n = 72L, seed = 10032L) {
  withr::with_seed(seed, {
    gene_a <- stats::rnorm(n)
    gene_b <- stats::rnorm(n)
    gene_c <- stats::rnorm(n)
    gene_d <- stats::rnorm(n)
    linear_predictor <- 0.8 * gene_a - 0.6 * gene_b +
      0.35 * gene_c + 0.2 * gene_d
    event_time <- stats::rexp(n, rate = 0.08 * exp(linear_predictor))
    censor_time <- stats::rexp(n, rate = 0.03)
    data.frame(
      patient_id = sprintf("ASVM%03d", seq_len(n)),
      time = pmax(pmin(event_time, censor_time), 0.01),
      event = as.integer(event_time <= censor_time),
      `gene A` = gene_a,
      `gene-B` = gene_b,
      gene_c = gene_c,
      gene_d = gene_d,
      check.names = FALSE
    )
  })
}

test_that("additive survival-SVM formulations are registered native learners", {
  learners <- c(
    "survivalsvm_additive",
    "survivalsvm_vanbelle1_additive",
    "survivalsvm_vanbelle2_additive",
    "survivalsvm_hybrid_additive"
  )
  expect_setequal(sb_survivalsvm_additive_learner_ids(), learners)
  expect_setequal(
    sb_survivalsvm_learner_ids(),
    c(
      sb_survivalsvm_linear_learner_ids(),
      sb_survivalsvm_kernel_learner_ids()
    )
  )
  expect_length(intersect(learners, sb_linear_survival_learner_ids()), 0L)

  expected_types <- c("regression", "vanbelle1", "vanbelle2", "hybrid")
  expected_selectors <- expected_fixed_learner_selectors("survival")
  registry <- sb_pipeline_registry()
  for (index in seq_along(learners)) {
    learner <- learners[[index]]
    params <- sb_survivalsvm_params(0.5, learner)
    expect_identical(params$type, expected_types[[index]])
    expect_identical(params$kernel, "add_kernel")
    expect_identical(params$gamma_mu, 0.5)
    expect_identical(params$standardize, TRUE)
    expect_identical(params$opt_meth, "quadprog")
    if (identical(params$type, "regression")) {
      expect_false("diff_meth" %in% names(params))
    } else {
      expect_identical(params$diff_meth, "makediff3")
    }

    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_setequal(rows$selector, expected_selectors)
    expect_true(all(grepl(
      "survivalsvm", rows$optional_dependency,
      fixed = TRUE
    )))
    expect_identical(unique(rows$tuning_semantics), "none")

    pipeline <- sb_pipeline(
      learner = learner,
      survivalsvm_gamma_mu = 0.5
    )
    expect_identical(pipeline$learner_params$kernel, "add_kernel")
    expect_identical(
      pipeline$risk_type,
      "calibrated_survival_svm_kernel_score"
    )
    expect_identical(
      pipeline$risk_reference,
      "survivalsvm_additive_kernel_native"
    )
    expect_identical(
      pipeline$risk_formula,
      "negative_native_kernel_score_times_analysis_cox_coefficient"
    )
  }
})

test_that("additive survival-SVM states exactly replay native kernel scores", {
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  data <- make_survivalsvm_additive_fixture()
  predictors <- c("gene A", "gene-B", "gene_c", "gene_d")
  x <- data[predictors]
  rownames(x) <- data$patient_id
  assessment <- x[c(3L, 18L, 41L, 67L), , drop = FALSE]
  learners <- sb_survivalsvm_additive_learner_ids()
  risks <- list()
  withr::local_seed(10033L)

  for (learner in learners) {
    adapter <- sb_survivalsvm_learner_adapter(
      gamma_mu = 1,
      learner = learner
    )
    rng_before <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, rng_before)

    scaled_assessment <- as.data.frame(
      sweep(
        sweep(as.matrix(assessment), 2L, model$audit$center, "-"),
        2L,
        model$audit$scale,
        "/"
      ),
      optional = TRUE,
      stringsAsFactors = FALSE
    )
    names(scaled_assessment) <- predictors
    native <- -as.numeric(stats::predict(
      model$engine,
      newdata = scaled_assessment
    )$predicted) * model$audit$calibration_coefficient
    observed <- adapter$predict(model, assessment)
    expect_equal(observed, native, tolerance = 1e-10)

    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = 7,
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, assessment)
    expected_survival <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      assessment,
      prediction_horizon = 7
    )$probability

    expect_s3_class(state, "stormbreaker_survivalsvm_additive_learner_state")
    expect_false("engine" %in% names(state))
    expect_null(rownames(state$basis))
    expect_identical(
      state$training_data_retention,
      "required_scaled_kernel_basis_without_outcomes"
    )
    expect_identical(dim(state$basis), c(length(state$row_weight), 4L))
    expect_equal(replay$risk, observed, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      expected_survival,
      tolerance = 0
    )
    risks[[learner]] <- observed
  }

  for (pair in utils::combn(names(risks), 2L, simplify = FALSE)) {
    expect_gt(max(abs(risks[[pair[[1L]]]] - risks[[pair[[2L]]]])), 1e-6)
  }
})

test_that("additive survival-SVM formulations execute fold-local CV and final fit", {
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  data <- make_survivalsvm_additive_fixture(seed = 10034L)
  predictors <- c("gene A", "gene-B", "gene_c", "gene_d")
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )

  for (learner in sb_survivalsvm_additive_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 3L,
      survivalsvm_gamma_mu = 1
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 10035L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:6, c("patient_id", rev(predictors)), drop = FALSE],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(comparison$complete_survival_probability_oof, TRUE)
    expect_identical(sum(comparison$runs$status == "success"), 3L)
    expect_identical(fit$status, "success")
    expect_identical(fit$learner_state$learner, learner)
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$risk, 6L)
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})

# ---- survivalsvm-centered-baseline ----

test_that("kernel probabilities are invariant to large risk offsets", {
  x <- data.frame(a = c(-1, 0, 1, 2), b = c(1, -1, 0, 2))
  time <- c(1, 2, 3, 4)
  event <- c(1, 1, 0, 1)
  for (learner in c("survivalsvm_rbf", "survivalsvm_additive")) {
    params <- sb_survivalsvm_params(1, learner)
    model <- sb_new_survivalsvm_kernel_model(
      learner, names(x), as.matrix(x), c(.1, -.1, .2, -.2), 0,
      1, if (params$kernel == "rbf_kernel") 1 else NULL,
      params, list(center = c(0, 0), scale = c(1, 1)), NULL, "test"
    )
    raw <- sb_predict_survivalsvm_kernel(model, x)
    expected <- exp(-vapply(raw, function(r) {
      sum(vapply(c(1, 2), function(t) {
        sum(event == 1 & time == t) / sum(exp(raw[time >= t] - r))
      }, numeric(1)))
    }, numeric(1)))
    context <- list(x = x, time = time, event = event,
                    prediction_horizon = 2, time_unit = "day")
    for (offset in c(0, -2000, 2000)) {
      shifted <- model
      shifted$intercept <- offset
      direct <- sb_predict_survivalsvm_kernel_survival(shifted, x, time, event, x, 2)
      expect_equal(direct$probability, expected, tolerance = 1e-10)
      state <- sb_extract_survivalsvm_kernel_state(shifted, context)
      saved <- sb_predict_survivalsvm_kernel_state(unserialize(serialize(state, NULL)), x)
      expect_equal(saved$risk, raw - offset)
      expect_equal(saved$survival_probability, expected, tolerance = 1e-10)
      expect_equal(state$baseline_hazard$risk_reference, mean(raw - offset))
      for (bad in list(NA_real_, Inf, "zero", c(0, 1))) {
        broken <- state
        broken$baseline_hazard$risk_reference <- bad
        expect_error(sb_predict_survivalsvm_kernel_state(broken, x),
                     class = "stormbreaker_bad_learner_state")
      }
      broken <- state
      broken$baseline_hazard$prediction_horizon <- -1
      expect_error(sb_predict_survivalsvm_kernel_state(broken, x))
      broken <- state
      broken$baseline_hazard$cumulative_baseline_hazard <- -1
      expect_error(sb_predict_survivalsvm_kernel_state(broken, x))
    }
    legacy <- sb_extract_survivalsvm_kernel_state(model, context)
    legacy$baseline_hazard$risk_reference <- NULL
    legacy$baseline_hazard$cumulative_baseline_hazard <-
      sb_breslow_baseline_hazard(time, event, raw, 2)
    expect_equal(sb_predict_survivalsvm_kernel_state(legacy, x)$survival_probability,
                 expected, tolerance = 1e-10)
  }
})

# ---- survivalsvm-formulations ----

make_survivalsvm_formulation_fixture <- function(n = 72L, seed = 9303L) {
  withr::with_seed(seed, {
    gene_a <- stats::rnorm(n)
    gene_b <- stats::rnorm(n)
    gene_c <- stats::rnorm(n)
    linear_predictor <- 0.7 * gene_a - 0.5 * gene_b + 0.25 * gene_c
    event_time <- stats::rexp(n, rate = 0.08 * exp(linear_predictor))
    censor_time <- stats::rexp(n, rate = 0.035)
    data.frame(
      patient_id = sprintf("SVM%03d", seq_len(n)),
      time = pmax(pmin(event_time, censor_time), 0.01),
      event = as.integer(event_time <= censor_time),
      `gene A` = gene_a,
      `gene-B` = gene_b,
      gene_c = gene_c,
      check.names = FALSE
    )
  })
}

test_that("survival-SVM formulations are distinct fixed native models", {
  learners <- c(
    "survivalsvm_vanbelle1",
    "survivalsvm_vanbelle2",
    "survivalsvm_hybrid"
  )
  expect_setequal(
    sb_survivalsvm_linear_learner_ids(),
    c("survivalsvm", learners)
  )

  expected_types <- c("vanbelle1", "vanbelle2", "hybrid")
  for (index in seq_along(learners)) {
    params <- sb_survivalsvm_params(
      gamma_mu = 0.5,
      learner = learners[[index]]
    )
    expect_identical(params$gamma_mu, 0.5)
    expect_identical(params$type, expected_types[[index]])
    expect_identical(params$diff_meth, "makediff3")
    expect_identical(params$kernel, "lin_kernel")
    expect_identical(params$opt_meth, "quadprog")
    expect_identical(params$standardize, TRUE)
    expect_identical(params$calibration_ties, "breslow")
  }

  regression <- sb_survivalsvm_params(gamma_mu = 0.5)
  expect_identical(regression$type, "regression")
  expect_false("diff_meth" %in% names(regression))

  registry <- sb_pipeline_registry()
  expected_selectors <- expected_fixed_learner_selectors("survival")
  for (index in seq_along(learners)) {
    learner <- learners[[index]]
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(unique(rows$task_type), "survival")
    expect_setequal(rows$selector, expected_selectors)
    expect_true(all(grepl(
      "survivalsvm", rows$optional_dependency,
      fixed = TRUE
    )))
    expect_identical(unique(rows$tuning_semantics), "none")

    pipeline <- sb_pipeline(
      learner = learner,
      survivalsvm_gamma_mu = 0.5
    )
    expect_identical(pipeline$learner, learner)
    expect_identical(pipeline$learner_params$type, expected_types[[index]])
    expect_identical(pipeline$learner_params$diff_meth, "makediff3")
    expect_identical(
      pipeline$id,
      paste0("none__", learner, "__gamma0.5")
    )
  }
})

test_that("new survival-SVM formulations preserve native linear scores", {
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  data <- make_survivalsvm_formulation_fixture()
  predictors <- c("gene A", "gene-B", "gene_c")
  x <- data[predictors]
  assessment <- x[c(3L, 18L, 41L, 67L), , drop = FALSE]
  learners <- c(
    "survivalsvm_vanbelle1",
    "survivalsvm_vanbelle2",
    "survivalsvm_hybrid"
  )
  native_types <- c("vanbelle1", "vanbelle2", "hybrid")
  risks <- list()
  withr::local_seed(9304L)

  for (index in seq_along(learners)) {
    learner <- learners[[index]]
    adapter <- sb_survivalsvm_learner_adapter(
      gamma_mu = 1,
      learner = learner
    )
    rng_before <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, rng_before)
    repeated <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, rng_before)
    expect_identical(model$coefficient, repeated$coefficient)

    scaled_assessment <- as.data.frame(
      sweep(
        sweep(as.matrix(assessment), 2L, model$audit$center, "-"),
        2L,
        model$audit$scale,
        "/"
      ),
      optional = TRUE,
      stringsAsFactors = FALSE
    )
    names(scaled_assessment) <- predictors
    native <- -as.numeric(stats::predict(
      model$engine,
      newdata = scaled_assessment
    )$predicted) * model$audit$calibration_coefficient
    observed <- adapter$predict(model, assessment)
    centered <- function(value) value - value[[1L]]

    expect_equal(centered(observed), centered(native), tolerance = 1e-10)
    expect_identical(model$learner, learner)
    expect_identical(model$params$type, native_types[[index]])
    expect_identical(model$audit$diff_meth, "makediff3")

    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = 7,
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, assessment)
    expected_survival <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      assessment,
      prediction_horizon = 7
    )$probability

    expect_identical(state$learner, learner)
    expect_false("engine" %in% names(state))
    expect_equal(replay$risk, observed, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      expected_survival,
      tolerance = 0
    )
    risks[[learner]] <- observed
  }

  for (pair in utils::combn(names(risks), 2L, simplify = FALSE)) {
    expect_gt(max(abs(risks[[pair[[1L]]]] - risks[[pair[[2L]]]])), 1e-6)
  }
})

test_that("new survival-SVM formulations execute fold-local CV and final fit", {
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  data <- make_survivalsvm_formulation_fixture(seed = 9305L)
  predictors <- c("gene A", "gene-B", "gene_c")
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )
  learners <- c(
    "survivalsvm_vanbelle1",
    "survivalsvm_vanbelle2",
    "survivalsvm_hybrid"
  )

  for (learner in learners) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 2L,
      survivalsvm_gamma_mu = 1
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 9306L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:6, c("patient_id", rev(predictors)), drop = FALSE],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(comparison$complete_survival_probability_oof, TRUE)
    expect_identical(sum(comparison$runs$status == "success"), 3L)
    expect_identical(
      comparison$provenance$survivalsvm_version,
      as.character(utils::packageVersion("survivalsvm"))
    )
    expect_identical(fit$status, "success")
    expect_identical(fit$learner_state$learner, learner)
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$risk, 6L)
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})

# ---- survivalsvm-rbf-kernel ----

make_survivalsvm_rbf_fixture <- function(n = 72L, seed = 11032L) {
  withr::with_seed(seed, {
    gene_a <- stats::rnorm(n)
    gene_b <- stats::rnorm(n)
    gene_c <- stats::rnorm(n)
    gene_d <- stats::rnorm(n)
    nonlinear_predictor <- 0.7 * gene_a - 0.5 * gene_b +
      0.35 * gene_c^2 - 0.25 * gene_a * gene_d
    event_time <- stats::rexp(
      n,
      rate = 0.07 * exp(nonlinear_predictor)
    )
    censor_time <- stats::rexp(n, rate = 0.03)
    data.frame(
      patient_id = sprintf("RBF%03d", seq_len(n)),
      time = pmax(pmin(event_time, censor_time), 0.01),
      event = as.integer(event_time <= censor_time),
      `gene A` = gene_a,
      `gene-B` = gene_b,
      gene_c = gene_c,
      gene_d = gene_d,
      check.names = FALSE
    )
  })
}

test_that("kernel survival-SVM one-input boundaries are exact", {
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")
  data <- sb_example_survival(n = 200L, seed = 20260905L)[1:160, ]
  x <- data["gene_a"]
  unsupported <- c(
    sb_survivalsvm_additive_learner_ids(),
    setdiff(
      sb_survivalsvm_rbf_learner_ids(),
      "survivalsvm_rbf"
    )
  )

  for (learner in unsupported) {
    expect_error(
      sb_survivalsvm_learner_adapter(1, learner)$fit(
        x, data$time, data$event
      ),
      class = "stormbreaker_survivalsvm_too_few_features",
      info = learner
    )
  }

  expect_s3_class(
    sb_survivalsvm_learner_adapter(
      1, "survivalsvm_rbf"
    )$fit(x, data$time, data$event),
    "stormbreaker_survivalsvm_rbf"
  )
  expect_s3_class(
    sb_survivalsvm_learner_adapter(
      1, "survivalsvm_vanbelle1"
    )$fit(x, data$time, data$event),
    "stormbreaker_survivalsvm_vanbelle1"
  )
})

test_that("RBF survival-SVM formulations are registered native learners", {
  learners <- c(
    "survivalsvm_rbf",
    "survivalsvm_vanbelle1_rbf",
    "survivalsvm_vanbelle2_rbf",
    "survivalsvm_hybrid_rbf"
  )
  expect_setequal(sb_survivalsvm_rbf_learner_ids(), learners)
  expect_setequal(
    sb_survivalsvm_kernel_learner_ids(),
    c(sb_survivalsvm_additive_learner_ids(), learners)
  )
  expect_setequal(
    sb_survivalsvm_learner_ids(),
    c(
      sb_survivalsvm_linear_learner_ids(),
      sb_survivalsvm_kernel_learner_ids()
    )
  )
  expect_length(intersect(learners, sb_linear_survival_learner_ids()), 0L)

  expected_types <- c("regression", "vanbelle1", "vanbelle2", "hybrid")
  expected_selectors <- expected_fixed_learner_selectors("survival")
  registry <- sb_pipeline_registry()
  for (index in seq_along(learners)) {
    learner <- learners[[index]]
    params <- sb_survivalsvm_params(0.5, learner)
    expect_identical(params$type, expected_types[[index]])
    expect_identical(params$kernel, "rbf_kernel")
    expect_null(params$kernel_pars)
    expect_identical(
      params$kernel_parameter_rule,
      "analysis_median_squared_distance"
    )
    expect_identical(params$gamma_mu, 0.5)
    expect_identical(params$standardize, TRUE)
    expect_identical(params$opt_meth, "quadprog")
    if (identical(params$type, "regression")) {
      expect_false("diff_meth" %in% names(params))
    } else {
      expect_identical(params$diff_meth, "makediff3")
    }

    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_setequal(rows$selector, expected_selectors)
    expect_true(all(grepl(
      "survivalsvm", rows$optional_dependency,
      fixed = TRUE
    )))
    expect_identical(unique(rows$tuning_semantics), "none")
    expect_true(all(grepl(
      "analysis-median-squared-distance", rows$parameter_semantics,
      fixed = TRUE
    )))

    pipeline <- sb_pipeline(
      learner = learner,
      survivalsvm_gamma_mu = 0.5
    )
    expect_identical(pipeline$learner_params$kernel, "rbf_kernel")
    expect_identical(
      pipeline$risk_type,
      "calibrated_survival_svm_kernel_score"
    )
    expect_identical(
      pipeline$risk_reference,
      "survivalsvm_rbf_kernel_native"
    )
    expect_identical(
      pipeline$risk_formula,
      "negative_native_kernel_score_times_analysis_cox_coefficient"
    )
  }
})

test_that("RBF survival-SVM kernel matrix follows the native formula", {
  basis <- matrix(c(0, 0, 1, -1), nrow = 2L, byrow = TRUE)
  assessment <- matrix(c(0, 1, 2, -1), nrow = 2L, byrow = TRUE)
  observed <- sb_survivalsvm_rbf_kernel_matrix(
    basis,
    assessment,
    kernel_pars = c(0.5, 0.5)
  )
  expected <- matrix(
    c(exp(-1), exp(-5), exp(-5), exp(-1)),
    nrow = 2L,
    byrow = TRUE
  )

  expect_equal(observed, expected, tolerance = 1e-15)
})

test_that("RBF survival-SVM states exactly replay native kernel scores", {
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  data <- make_survivalsvm_rbf_fixture()
  predictors <- c("gene A", "gene-B", "gene_c", "gene_d")
  x <- data[predictors]
  rownames(x) <- data$patient_id
  assessment <- x[c(4L, 19L, 42L, 68L), , drop = FALSE]
  learners <- sb_survivalsvm_rbf_learner_ids()
  risks <- list()
  withr::local_seed(11033L)

  for (learner in learners) {
    adapter <- sb_survivalsvm_learner_adapter(
      gamma_mu = 1,
      learner = learner
    )
    rng_before <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, rng_before)

    scaled_assessment <- as.data.frame(
      sweep(
        sweep(as.matrix(assessment), 2L, model$audit$center, "-"),
        2L,
        model$audit$scale,
        "/"
      ),
      optional = TRUE,
      stringsAsFactors = FALSE
    )
    names(scaled_assessment) <- predictors
    native <- -as.numeric(stats::predict(
      model$engine,
      newdata = scaled_assessment
    )$predicted) * model$audit$calibration_coefficient
    observed <- adapter$predict(model, assessment)
    expect_equal(observed, native, tolerance = 1e-10)
    expect_equal(
      model$audit$kernel_pars,
      stats::median(as.numeric(stats::dist(scale(x)))^2),
      tolerance = 0
    )

    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = 7,
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, assessment)
    expected_survival <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      assessment,
      prediction_horizon = 7
    )$probability

    expect_s3_class(state, "stormbreaker_survivalsvm_rbf_learner_state")
    expect_false("engine" %in% names(state))
    expect_null(rownames(state$basis))
    expect_identical(
      state$training_data_retention,
      "required_scaled_kernel_basis_without_outcomes"
    )
    expect_identical(dim(state$basis), c(length(state$row_weight), 4L))
    expect_equal(state$kernel_pars, model$audit$kernel_pars, tolerance = 0)
    expect_equal(replay$risk, observed, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      expected_survival,
      tolerance = 0
    )
    risks[[learner]] <- observed
  }

  for (pair in utils::combn(names(risks), 2L, simplify = FALSE)) {
    expect_gt(max(abs(risks[[pair[[1L]]]] - risks[[pair[[2L]]]])), 1e-6)
  }
})

test_that("RBF survival-SVM formulations execute fold-local CV and final fit", {
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  data <- make_survivalsvm_rbf_fixture(seed = 11034L)
  predictors <- c("gene A", "gene-B", "gene_c", "gene_d")
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )

  for (learner in sb_survivalsvm_rbf_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 3L,
      survivalsvm_gamma_mu = 1
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 11035L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:6, c("patient_id", rev(predictors)), drop = FALSE],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(comparison$complete_survival_probability_oof, TRUE)
    expect_identical(sum(comparison$runs$status == "success"), 3L)
    expect_identical(fit$status, "success")
    expect_identical(fit$learner_state$learner, learner)
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$risk, 6L)
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})

# ---- svm-probability-calibration ----

test_that("SVM calibration keeps paired patients and fold-local Top-K", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(100L, 3800L)
  d$outcome <- rep(c(0L, 1L), 50L)
  d$group <- rep(paste0("pair", seq_len(50L)), each = 2L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors,
                         group_col = "group")
  pipeline <- sb_pipeline(selector = "univ_binary_auc_top_k",
    selector_n_features = 2L, learner = "svm_linear_binary")
  fit <- sb_fit_final(task, pipeline)
  expect_identical(fit$status, "success")
  expect_identical(fit$tuning$status, "success")
  expect_identical(fit$tuning$inner_folds, 5L)
  expect_length(fit$tuning$inner_selections, 5L)
  for (audit in fit$tuning$inner_selections) {
    train <- match(audit$analysis_patient_ids, d$patient_id)
    assessment <- match(audit$assessment_patient_ids, d$patient_id)
    expect_length(intersect(d$group[train], d$group[assessment]), 0L)
    direct <- sb_fit_univ_binary_auc_top_k(d[train, predictors],
                                           d$outcome[train], 2L)
    expect_identical(audit$selected_features, direct$selected_features)
  }
  cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 3801L)
  expect_identical(cv$status, "success")
  expect_length(cv$tuning, 3L)
  for (audit in c(cv$tuning, list(fit$tuning))) {
    rows <- match(names(audit$fold_id), d$patient_id)
    expect_identical(anyNA(rows), FALSE)
    group_folds <- split(unname(audit$fold_id), d$group[rows])
    expect_identical(unname(vapply(group_folds, function(x) {
      length(unique(x))
    }, integer(1))), rep(1L, length(group_folds)))
    expect_identical(audit$calibration$method, "libsvm_platt_oof")
  }
})

test_that("unscaled unscreened calibration reproduces native LIBSVM", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(100L, 3802L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  for (first in c(0L, 1L)) for (learner in c(
    "svm_linear_binary", "svm_radial_binary",
    "svm_polynomial_binary", "svm_sigmoid_binary"
  )) {
    rows <- c(which(d$outcome == first), which(d$outcome != first))
    params <- sb_classical_binary_params(learner)
    params$standardize <- FALSE
    native <- sb_fit_classical_binary(x[rows, ], d$outcome[rows],
                                             learner, params)
    tuned <- sb_calibrate_svm_probability(x[rows, ], d$outcome[rows],
      learner, params, d$patient_id[rows])
    expect_equal(tuned$calibration$A, unname(native$engine$probA),
                 tolerance = 1e-8)
    expect_equal(tuned$calibration$B, unname(native$engine$probB),
                 tolerance = 1e-8)
    expect_identical(tuned$calibration$positive_class, as.character(first))
  }
})

test_that("all five SVM variants use native fold-local margins and portable probabilities", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(100L, 3803L)
  d$outcome <- rep(c(0L, 1L, 0L, 0L), 25L)
  d$group <- rep(paste0("pair", seq_len(50L)), each = 2L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors,
                         group_col = "group")
  native_fit <- function(rows, features, params) {
    args <- list(x = as.matrix(d[rows, features, drop = FALSE]),
      y = factor(d$outcome[rows], levels = c(0L, 1L)), type = "C-classification",
      kernel = params$kernel, cost = params$cost, scale = TRUE,
      probability = FALSE, fitted = FALSE)
    if (!is.null(params$gamma)) args$gamma <- 1 / length(features)
    if (!is.null(params$degree)) args$degree <- params$degree
    if (!is.null(params$coef0)) args$coef0 <- params$coef0
    if (!is.null(params$class_weights)) {
      counts <- table(args$y)
      args$class.weights <- setNames(length(rows) / (2 * as.numeric(counts)),
                                      names(counts))
    }
    do.call(e1071::svm, args)
  }
  set.seed(3804L)
  before <- .Random.seed
  for (learner in sb_svm_probability_learner_ids()) {
    pipeline <- sb_pipeline(selector = "univ_binary_auc_top_k",
      selector_n_features = 2L, learner = learner)
    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    expect_identical(.Random.seed, before)
    expected_decision <- numeric(nrow(d))
    for (inner in fit$tuning$inner_selections) {
      train <- match(inner$analysis_patient_ids, d$patient_id)
      assess <- match(inner$assessment_patient_ids, d$patient_id)
      features <- sb_fit_univ_binary_auc_top_k(d[train, predictors],
                                                d$outcome[train], 2L)$selected_features
      native <- native_fit(train, features, pipeline$learner_params)
      decision <- attr(predict(native, as.matrix(d[assess, features]),
                                decision.values = TRUE), "decision.values")
      sign <- if (native$levels[native$labels[1L]] ==
                  fit$tuning$calibration$positive_class) 1 else -1
      expected_decision[assess] <- sign * as.numeric(decision)
      expect_identical(inner$selected_features, features)
      expect_equal(inner$scale_center, native$x.scale[["scaled:center"]],
                   tolerance = 1e-14)
      expect_equal(inner$scale_scale, native$x.scale[["scaled:scale"]],
                   tolerance = 1e-14)
      if (!is.null(pipeline$learner_params$class_weights)) {
        counts <- table(factor(d$outcome[train], levels = c(0L, 1L)))
        expect_equal(inner$class_weights,
          setNames(length(train) / (2 * as.numeric(counts)), names(counts)),
          tolerance = 0)
      }
    }
    expect_equal(fit$tuning$calibration$out_of_fold$decision, expected_decision,
                 tolerance = 1e-12)
    positive <- d$outcome == as.integer(fit$tuning$calibration$positive_class)
    target <- ifelse(positive, (sum(positive) + 1) / (sum(positive) + 2),
                      1 / (sum(!positive) + 2))
    # Some native SVMs have almost constant margins. Condition the independent
    # optimizer's coordinates; this does not change production calibration.
    centered <- (expected_decision - mean(expected_decision)) / sd(expected_decision)
    objective <- function(beta) {
      z <- beta[1L] * centered + beta[2L]
      sum(pmax(z, 0) - (1 - target) * z + log1p(exp(-abs(z))))
    }
    gradient <- function(beta) {
      z <- beta[1L] * centered + beta[2L]
      residual <- target - plogis(-z)
      c(sum(centered * residual), sum(residual))
    }
    optimum <- optim(c(0, 0), objective, gradient, method = "BFGS",
                      control = list(reltol = 1e-13))
    expect_identical(optimum$convergence, 0L)
    transformed <- c(fit$tuning$calibration$A * sd(expected_decision),
      fit$tuning$calibration$B + fit$tuning$calibration$A * mean(expected_decision))
    expect_equal(transformed, optimum$par, tolerance = 1e-5)
    expect_equal(objective(transformed), optimum$value, tolerance = 1e-9)
    native <- native_fit(seq_len(nrow(d)), fit$required_features,
                           pipeline$learner_params)
    native$probA <- fit$tuning$calibration$A
    native$probB <- fit$tuning$calibration$B
    native$compprob <- TRUE
    expected <- unname(attr(predict(native, as.matrix(d[fit$required_features]),
                                     probability = TRUE), "probabilities")[, "1"])
    portable <- unserialize(serialize(fit, NULL))
    actual <- sb_predict_final(portable, d[c("patient_id", predictors)],
                                "patient_id")
    expect_identical(actual$status, "success")
    expect_equal(actual$predictions$probability, expected, tolerance = 1e-14)
    expect_identical(portable$learner_state$params, pipeline$learner_params)
  }
})

test_that("calibration training is unchanged by its held-out fold", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(100L, 3805L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  learner <- "svm_radial_binary"
  params <- sb_classical_binary_params(learner)
  selector <- sb_univ_binary_auc_top_k_selector_adapter(2L)
  declared <- sb_make_svm_calibration_folds(d$outcome, params$seed)
  fold <- if (declared$fold_id[1L] == 1L) 2L else 1L
  local_mocked_bindings(sb_make_svm_calibration_folds = function(...) declared)
  models <- list()
  original <- sb_fit_classical_binary
  local_mocked_bindings(sb_fit_classical_binary = function(...) {
    model <- original(...)
    models[[length(models) + 1L]] <<- model
    model
  })
  baseline <- sb_calibrate_svm_probability(d[predictors], d$outcome, learner,
    params, d$patient_id, selector = selector)
  changed <- d
  assess <- declared$fold_id == fold
  changed[assess, predictors] <- changed[assess, predictors] + 500
  changed$outcome[assess] <- 1L - changed$outcome[assess]
  perturbed <- sb_calibrate_svm_probability(changed[predictors], changed$outcome,
    learner, params, d$patient_id, selector = selector)
  expect_identical(perturbed$inner_selections[[fold]],
                    baseline$inner_selections[[fold]])
  first <- models[[fold]]
  second <- models[[fold + 5L]]
  probe <- as.matrix(d[first$features])
  expect_equal(attr(predict(first$engine, probe, decision.values = TRUE),
                     "decision.values"),
    attr(predict(second$engine, probe, decision.values = TRUE), "decision.values"),
    tolerance = 0)
})

test_that("SVM calibration preserves fixed blocks and singleton native folds", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(60L, 3806L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  for (selector in c("none", "fixed_features")) {
    pipeline <- if (selector == "none") {
      sb_pipeline(learner = "svm_radial_binary")
    } else sb_pipeline(selector = "fixed_features",
      selector_features = c("gene_c", "gene_a"), feature_set_id = "predeclared",
      learner = "svm_radial_binary")
    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    for (inner in fit$tuning$inner_selections) {
      expect_identical(inner$selected_features, fit$required_features)
    }
    expect_identical(fit$tuning$inner_selection_method,
      if (selector == "none") "none" else "predeclared_feature_subset")
  }
  plain <- sb_make_svm_calibration_folds(d$outcome, 1L)
  singleton <- sb_make_svm_calibration_folds(d$outcome, 1L, d$patient_id)
  expect_identical(singleton, plain)
})

test_that("single-class calibration subproblems retain native constant margins", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(10L, 3807L)
  d$outcome <- c(1L, rep(0L, 9L))
  learner <- "svm_linear_binary"
  params <- sb_classical_binary_params(learner)
  params$standardize <- FALSE
  x <- d[c("gene_a", "gene_b")]
  native <- sb_fit_classical_binary(x, d$outcome, learner, params)
  tuned <- sb_calibrate_svm_probability(x, d$outcome, learner, params,
                                        d$patient_id)
  expect_equal(tuned$calibration$A, unname(native$engine$probA), tolerance = 1e-8)
  expect_equal(tuned$calibration$B, unname(native$engine$probB), tolerance = 1e-8)
  expect_equal(sum(vapply(tuned$inner_selections, function(x) {
    x$status == "single_class_constant"
  }, logical(1))), 1L)
})

test_that("infeasible grouped calibration is not silently replaced by row folds", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(20L, 3808L)
  d$outcome <- rep(c(0L, 1L), 10L)
  d$group <- rep(paste0("group", seq_len(4L)), each = 5L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"),
                         group_col = "group")
  fit <- sb_fit_final(task, sb_pipeline(learner = "svm_linear_binary"))
  expect_identical(fit$status, "failed")
  expect_identical(fit$error_stage, "tuning")
  expect_identical(fit$error_class, "stormbreaker_infeasible_inner_folds")
  expect_identical(fit$tuning$status, "failed")
  expect_identical(fit$tuning$inner_folds, 5L)
  expect_null(fit$learner_state)
})

test_that("SVM calibration failures remain visible in CV and final fitting", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(60L, 3809L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  pipeline <- sb_pipeline(learner = "svm_linear_binary")
  local_mocked_bindings(sb_calibrate_svm_probability = function(...) {
    sb_abort("Synthetic calibration failure", "stormbreaker_synthetic_calibration")
  })
  final <- sb_fit_final(task, pipeline)
  cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 3810L)
  expect_identical(final$status, "failed")
  expect_identical(final$error_stage, "tuning")
  expect_identical(cv$status, "failed")
  expect_length(cv$tuning, 3L)
  for (audit in c(cv$tuning, list(final$tuning))) {
    expect_identical(audit$status, "failed")
    expect_identical(audit$error_class, "stormbreaker_synthetic_calibration")
    expect_identical(audit$inner_folds, 5L)
  }
})

test_that("rare classes do not block otherwise feasible grouped calibration", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(40L, 3811L)
  d$outcome <- c(1L, 1L, rep(0L, 38L))
  d$group <- rep(paste0("group", seq_len(20L)), each = 2L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"),
                         group_col = "group")
  fit <- sb_fit_final(task, sb_pipeline(learner = "svm_linear_binary"))
  expect_identical(fit$status, "success")
  expect_identical(sort(unique(unname(fit$tuning$fold_id))), seq_len(5L))
  for (inner in fit$tuning$inner_selections) {
    train <- match(inner$analysis_patient_ids, d$patient_id)
    assess <- match(inner$assessment_patient_ids, d$patient_id)
    expect_length(intersect(d$group[train], d$group[assess]), 0L)
    if (length(unique(d$outcome[train])) == 1L) {
      expect_identical(inner$status, "single_class_constant")
      oof <- fit$tuning$calibration$out_of_fold
      expect_identical(oof$decision[match(inner$assessment_patient_ids, oof$patient_id)],
                         rep(-1, length(assess)))
    }
  }
  # Exactly K groups must also yield K nonempty folds when only one has class 1.
  for (seed in 1:10) {
    y <- c(1L, 1L, rep(0L, 8L))
    folds <- sb_make_svm_calibration_folds(y, seed,
      rep(paste0("g", seq_len(5L)), each = 2L))
    expect_identical(sort(unique(folds$fold_id)), seq_len(5L))
    expect_identical(unname(tabulate(folds$fold_id)), rep(2L, 5L))
  }
  error <- capture_sb_error(sb_make_group_stratified_fold_id(
    d$outcome, d$group, 5L, 1L, "stormbreaker_infeasible_inner_folds"))
  expect_s3_class(error, "stormbreaker_infeasible_inner_folds")
})
})

# ---- loclda-posterior-extraction ----
local({
test_that("locLDA unnamed posterior follows native grouping order", {
  skip_if_not_installed("klaR")
  posterior <- matrix(c(0.8, 0.3, 0.2, 0.7), nrow = 2L)
  local_mocked_bindings(
    predict.loclda = function(object, newdata) list(posterior = posterior),
    .package = "klaR"
  )
  engine <- list(grouping = factor(c(0, 1, 0), levels = c(0, 1)))
  expect_equal(sb_predict_loclda_binary(engine, matrix(1, 2, 1)),
               c(0.2, 0.7))
  engine$grouping <- factor(c(0, 1, 0), levels = c(1, 0))
  expect_equal(sb_predict_loclda_binary(engine, matrix(1, 2, 1)),
               c(0.8, 0.3))
})

test_that("locLDA named columns and invalid native values are preserved", {
  skip_if_not_installed("klaR")
  posterior <- matrix(c(0.2, NaN, 0.8, NaN), nrow = 2L,
                      dimnames = list(NULL, c("1", "0")))
  local_mocked_bindings(
    predict.loclda = function(object, newdata) list(posterior = posterior),
    .package = "klaR"
  )
  engine <- list(grouping = factor(c(0, 1)))
  expect_equal(sb_predict_loclda_binary(engine, matrix(1, 2, 1)),
               c(0.2, NaN))
  colnames(posterior) <- NULL
  expect_equal(sb_predict_loclda_binary(engine, matrix(1, 2, 1)),
               c(0.8, NaN))
})

test_that("locLDA native finite unnamed probabilities remain usable", {
  skip_if_not_installed("klaR")
  x <- withr::with_seed(1907L, matrix(rnorm(160L), ncol = 2L))
  grouping <- factor(rep(c(0, 1), 40L), levels = c(0, 1))
  engine <- klaR::loclda(x, grouping, k = 30L, weighted.apriori = TRUE)
  assessment <- matrix(c(1000, 1001, 1000, 1001), ncol = 2L)
  native <- get("predict.loclda", asNamespace("klaR"))(
    engine, assessment
  )$posterior
  expect_null(colnames(native))
  expect_equal(unname(native[, 2L]), c(0.5, 0.5))
  expect_equal(unname(sb_predict_loclda_binary(engine, assessment)),
               unname(native[, 2L]))
})

test_that("locLDA does not guess unidentified posterior classes", {
  skip_if_not_installed("klaR")
  local_mocked_bindings(
    predict.loclda = function(object, newdata) {
      list(posterior = matrix(0.5, nrow = 2L, ncol = 2L))
    },
    .package = "klaR"
  )
  expect_error(
    sb_predict_loclda_binary(list(grouping = factor("0")),
                                   matrix(1, 2, 1)),
    class = "stormbreaker_invalid_classical_binary_prediction"
  )
})
})

# ---- nnet-iteration-budget ----
local({
test_that("nnet gains optimization time without changing saved budgets", {
  p <- sb_pipeline(learner = "nnet_binary")
  expect_identical(p$learner_params$maxit, 2000L)
  expect_identical(sb_validate_pipeline_contract(p, "binary"), p)
  old <- p
  old$learner_params$maxit <- 200L
  expect_identical(sb_validate_pipeline_contract(old, "binary"), old)
  old$learner_params$decay <- 10
  error <- tryCatch(sb_validate_pipeline_contract(old, "binary"), error = identity)
  expect_s3_class(error, "stormbreaker_unresolved_pipeline")
})
})

# ---- neural-weight-scale ----
local({
test_that("neural color mapping preserves extreme and small signed weights", {
  d <- data.frame(unit = 1:5, weight = c(-3000, -1, 0, 1, 2000))
  p <- ggplot2::ggplot(d, ggplot2::aes(unit, 1, fill = weight)) +
    ggplot2::geom_tile() + sb_neural_weight_scale()
  b <- ggplot2::ggplot_build(p)
  expect_identical(p$data, d)
  expect_equal(nrow(b$data[[1]]), nrow(d))
  expect_false(anyNA(b$data[[1]]$fill))
  expect_equal(length(unique(b$data[[1]]$fill)), nrow(d))
  scale <- b$plot$scales$get_scales("fill")
  expect_equal(scale$trans$inverse(scale$trans$transform(d$weight)), d$weight)
})
})
