# ---- algorithm-rrf ----
local({
test_that("RRF forest views retain native importance, all trees and OOB histories", {
  skip_if_not_installed("RRF")
  d <- sb_example_binary(160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features),
    sb_pipeline(learner = "rrf_binary"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  prediction <- sb_predict_final(fit, d, "patient_id")
  e <- fit$diagnostic_artifacts$learner$engine
  importance <- sb_forest_diagnostics(fit)
  sizes <- sb_forest_diagnostics(fit, "tree_size")
  oob <- sb_forest_diagnostics(fit, "oob_error")
  expect_equal(importance$table$feature, rownames(e$importance))
  expect_equal(importance$table$importance, unname(e$importance[, "MeanDecreaseGini"]), tolerance = 0)
  expect_equal(nrow(importance$table), length(features))
  expect_gt(sum(importance$table$importance == 0), 0)
  native_sizes <- vapply(seq_len(e$ntree), function(k) {
    sum(RRF::getTree(e, k = k)[, "status"] == -1)
  }, integer(1))
  expect_equal(sizes$table$terminal_nodes, native_sizes, tolerance = 0)
  expect_equal(oob$table$error, as.numeric(e$err.rate), tolerance = 0)
  expect_equal(oob$table$available, is.finite(as.numeric(e$err.rate)))
  expect_equal(oob$table$trees, rep(seq_len(e$ntree), 3))
  compact <- fit
  compact$diagnostic_artifacts <- NULL
  expect_identical(sb_forest_diagnostics(compact), importance)
  expect_identical(sb_forest_diagnostics(compact, "tree_size"), sizes)
  expect_identical(sb_native_diagnostics(compact, complete = TRUE)$status$status,
    c("available", "available", "missing_capture"))
  bundle <- sb_native_diagnostics(fit, complete = TRUE)
  expect_identical(bundle$status$status, rep("available", 3L))
  for (v in bundle$diagnostics) expect_s3_class(
    ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(v)), "ggplot_built")
  expect_identical(sb_predict_final(fit, d, "patient_id"), prediction)
  expect_identical(serialize(fit, NULL, version = 2L), before)
})
})

# ---- algorithm-evtree ----
local({
test_that("evolutionary tree labels use the original predictor schema", {
  skip_if_not_installed("evtree")
  d <- sb_example_binary(160L, seed = 108L)
  d$gene_a <- d$gene_a + 1.5 * d$outcome
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d,"patient_id","outcome",features)
  fit <- sb_fit_final(task,sb_pipeline(learner="evtree_binary"))
  expect_identical(fit$status,"success")
  before <- serialize(fit,NULL,version=2L)
  prediction <- sb_predict_final(fit,d,"patient_id")
  e <- unserialize(fit$learner_state$engine_raw)
  z <- sb_evtree_structure(fit)
  expect_equal(nrow(e$data),0L)
  expect_equal(z$table$n[1L],nrow(d))
  expect_equal(sum(z$table$n[z$table$terminal]),nrow(d))
  check <- function(node) {
    key <- as.character(partykit::id_node(node))
    expect_identical(z$native_splits[[key]],partykit::split_node(node))
    kids <- partykit::kids_node(node)
    if(length(kids)) {
      split <- partykit::split_node(node)
      expected_feature <- names(e$data)[partykit::varid_split(split)]
      row <- match(as.integer(key),z$table$node)
      expect_identical(z$table$split_feature[row],expected_feature)
      expect_identical(expected_feature %in% features,TRUE)
      branches <- partykit::character_split(split,e$data,digits=7)
      children <- match(vapply(kids,partykit::id_node,1L),z$table$node)
      expect_identical(z$table$incoming_rule[children],paste(expected_feature,branches$levels))
      expect_equal(sum(z$table$n[children]),z$table$n[row])
      invisible(lapply(kids,check))
    }
  }
  check(partykit::node_party(e))
  expect_identical(sb_native_diagnostics(fit)$status$status,"available")
  p <- sb_plot_algorithm_diagnostic(z)
  expect_identical(p$labels$title,"Evolutionary classification tree")
  expect_s3_class(ggplot2::ggplot_build(p),"ggplot_built")
  expect_identical(sb_predict_final(fit,d,"patient_id"),prediction)
  expect_identical(serialize(fit,NULL,version=2L),before)
})
})

# ---- model-forests ----
local({
# ---- aorsf-survival-fixed-adapter ----

test_that("aorsf survival declares one fixed native forest", {
  expect_true(
    "aorsf_survival" %in% sb_additional_survival_learner_ids()
  )
  expect_identical(
    sb_additional_survival_dependency("aorsf_survival"),
    list(package = "aorsf", packages = "aorsf",
         available = requireNamespace("aorsf", quietly = TRUE))
  )
  expect_identical(
    sb_additional_survival_params("aorsf_survival"),
    list(
      num_trees = 500L,
      num_split = 5L,
      num_retry = 3L,
      mtry_rule = "ceiling_sqrt_p",
      sample_with_replacement = TRUE,
      sample_fraction = 0.632,
      control_method = "glm",
      control_scale_x = TRUE,
      control_ties = "efron",
      control_max_iter = 1L,
      seed = 1L,
      threads = 1L,
      importance = "none",
      attach_data = FALSE
    )
  )
})

test_that("aorsf survival risk and horizon survival equal native predictions", {
  skip_if_not_installed("aorsf")
  d <- make_sb_fixture(n = 60L, seed = 8431L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:7, , drop = FALSE]
  params <- sb_additional_survival_params("aorsf_survival")
  params$num_trees <- 40L
  adapter <- sb_additional_survival_learner_adapter(
    "aorsf_survival",
    params
  )
  model <- adapter$fit(x, d$time, d$event)
  risk <- adapter$predict(model, assessment)
  survival <- adapter$predict_survival(
    model,
    x,
    d$time,
    d$event,
    assessment,
    prediction_horizon = 7
  )
  native_risk <- stats::predict(
    model$engine,
    new_data = assessment,
    pred_type = "mort",
    n_thread = 1L
  )
  native_survival <- stats::predict(
    model$engine,
    new_data = assessment,
    pred_type = "surv",
    pred_horizon = 7,
    n_thread = 1L
  )
  state <- adapter$extract_state(model, list(
    x = x,
    time = d$time,
    event = d$event,
    prediction_horizon = 7,
    time_unit = "month"
  ))
  replay <- adapter$predict_state(
    unserialize(serialize(state, NULL)),
    assessment
  )

  expect_s3_class(model$engine, "ObliqueForestSurvival")
  expect_null(model$engine$data)
  expect_equal(risk, unname(as.numeric(native_risk)), tolerance = 1e-14)
  expect_equal(
    survival$probability,
    unname(as.numeric(native_survival)),
    tolerance = 1e-14
  )
  expect_equal(replay$risk, risk, tolerance = 0)
  expect_equal(
    replay$survival_probability,
    survival$probability,
    tolerance = 0
  )
  expect_identical(state$risk_reference, "aorsf_native")
  expect_identical(state$risk_formula, "native_mortality_over_event_times")
  expect_identical(
    state$survival_probability$source,
    "aorsf_oblique_survival_forest_direct"
  )
  expect_false("baseline_hazard" %in% names(Filter(Negate(is.null), state)))
})

test_that("aorsf keeps raw feature names while using unique native aliases", {
  skip_if_not_installed("aorsf")
  d <- make_sb_fixture(n = 60L, seed = 8431L)
  x <- d[c("gene_c", "gene_a", "gene_b")]
  names(x) <- c("NKX3-1", "NKX3.1", "if")
  native_x <- x
  names(native_x) <- make.names(names(x), unique = TRUE)
  params <- sb_additional_survival_params("aorsf_survival")
  params$num_trees <- 40L
  adapter <- sb_additional_survival_learner_adapter("aorsf_survival", params)
  model <- adapter$fit(x, d$time, d$event)
  renamed <- adapter$fit(native_x, d$time, d$event)
  assessment <- x[1:7, , drop = FALSE]
  native_assessment <- native_x[1:7, , drop = FALSE]
  state <- adapter$extract_state(model, list(
    x = x, time = d$time, event = d$event,
    prediction_horizon = 7, time_unit = "month"
  ))
  replay <- adapter$predict_state(unserialize(serialize(state, NULL)), assessment)
  native_risk <- unname(as.numeric(stats::predict(
    renamed$engine, new_data = native_assessment, pred_type = "mort", n_thread = 1L
  )))
  native_survival <- unname(as.numeric(stats::predict(
    renamed$engine, new_data = native_assessment,
    pred_type = "surv", pred_horizon = 7, n_thread = 1L
  )))

  expect_identical(model$features, names(x))
  expect_identical(state$features, names(x))
  expect_identical(model$engine$get_names_x(), names(native_x))
  expect_identical(model$engine$tree_seeds, renamed$engine$tree_seeds)
  expect_identical(adapter$predict(model, assessment), native_risk)
  expect_identical(replay$risk, native_risk)
  expect_identical(replay$survival_probability, native_survival)
  expect_identical(adapter$predict_survival(model, x, d$time, d$event,
    assessment, prediction_horizon = 7)$probability, native_survival)
})

test_that("aorsf survival executes through CV and final prediction", {
  skip_if_not_installed("aorsf")
  d <- make_sb_fixture(n = 60L, seed = 8432L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  original_params <- sb_additional_survival_params
  local_mocked_bindings(
    sb_additional_survival_params = function(learner) {
      params <- original_params(learner)
      if (identical(learner, "aorsf_survival")) {
        params$num_trees <- 25L
      }
      params
    }
  )
  pipeline <- sb_pipeline(learner = "aorsf_survival")
  cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 8433L)
  fit <- sb_fit_final(task, pipeline)
  prediction <- sb_predict_final(
    fit,
    d[1:6, c("patient_id", "gene_c", "gene_a", "gene_b")],
    "patient_id"
  )

  expect_identical(cv$status, "success")
  expect_identical(cv$complete_oof, TRUE)
  expect_identical(cv$complete_survival_probability_oof, TRUE)
  expect_identical(
    cv$provenance$aorsf_version,
    as.character(utils::packageVersion("aorsf"))
  )
  expect_identical(fit$status, "success")
  expect_identical(
    fit$provenance$aorsf_version,
    as.character(utils::packageVersion("aorsf"))
  )
  expect_identical(fit$learner_state$learner, "aorsf_survival")
  expect_identical(prediction$status, "success")
  expect_true(all(is.finite(prediction$predictions$risk)))
  expect_true(all(prediction$predictions$survival_probability >= 0 &
    prediction$predictions$survival_probability <= 1))
})

# ---- bagged-knn-survival-adapter ----

test_that("bagged distance screening preserves exact selected neighbors", {
  withr::local_seed(491)
  for (p in c(1L, 5L, 30L)) {
    for (scenario in c("ordinary", "duplicates", "offset", "ill_conditioned",
                       "near_ties")) {
      x <- matrix(rnorm(100L * p), 100L, p)
      if (scenario == "duplicates") x[11:20, ] <- x[1:10, ]
      if (scenario == "offset") x <- x + 1e9
      if (scenario == "ill_conditioned" && p > 1L) {
        x[, p] <- x[, 1L] + 1e-6 * x[, p]
      }
      if (scenario == "near_ties") {
        x[2L, ] <- x[1L, ] + 1e-14
        x[3L, ] <- x[1L, ] - 1e-14
      }
      rownames(x) <- paste0("row", seq_len(nrow(x)))
      covariance <- cov(x)
      inverse <- solve(covariance, tol = 1e-25)
      distance <- sb_bagged_knn_distance_function(10L)
      queries <- rbind(x[1:4, , drop = FALSE], colMeans(x))
      for (i in seq_len(nrow(queries))) {
        expected <- sort(stats::mahalanobis(x, queries[i, ], inverse,
          inverted = TRUE), index.return = TRUE)
        actual <- sort(distance(x, queries[i, ], covariance, tol = 1e-25),
          index.return = TRUE)
        expect_identical(actual$ix[1:10], expected$ix[1:10])
        expect_identical(actual$x[1:10], expected$x[1:10])
      }
    }
  }
})

test_that("bagged full-neighborhood distances and singular errors stay native", {
  x <- matrix(seq_len(20), ncol = 1L)
  distance <- sb_bagged_knn_distance_function(nrow(x))
  expect_identical(distance(x, 2.5, cov(x), tol = 1e-25),
    stats::mahalanobis(x, 2.5, cov(x), tol = 1e-25))
  singular <- cbind(x, x)
  condition <- function(expr) tryCatch(expr, error = identity)
  native <- condition(stats::mahalanobis(singular, c(1, 1), cov(singular)))
  got <- condition(sb_bagged_knn_distance_function(2L)(
    singular, c(1, 1), cov(singular)))
  expect_identical(class(got), class(native))
  expect_identical(conditionMessage(got), conditionMessage(native))
})

test_that("bagged single-feature single-query native curves remain unchanged", {
  skip_if_not_installed("bnnSurvival", minimum_version = "0.1.5")
  withr::local_seed(731)
  x <- data.frame("gene one" = rnorm(60L), check.names = FALSE)
  time <- rexp(60L) + 0.01
  event <- rep(c(1L, 0L, 1L), 20L)
  adapter <- sb_bagged_knn_survival_learner_adapter("bagged_knn_survival")
  model <- adapter$fit(x, time, event)
  model <- unserialize(serialize(model, NULL))
  query <- x[4L, , drop = FALSE]
  prepared <- sb_bagged_knn_prediction_data(model, query)
  native <- get("predict", asNamespace("bnnSurvival"))
  before <- .Random.seed
  reference <- sb_with_bagged_knn_single_core(
    native(prepared$model$engine, prepared$internal))
  actual <- sb_predict_bagged_knn_survival_curves(model, query)
  expect_identical(actual$survival,
    as.matrix(bnnSurvival::predictions(reference)))
  expect_identical(.Random.seed, before)
})

# ---- cox-control-separation ----

make_cox_control_fixture <- function() {
  withr::with_seed(12413L, {
    n <- 160L
    x <- as.data.frame(matrix(stats::rnorm(n * 4L), nrow = n))
    names(x) <- paste0("gene", seq_len(ncol(x)))
    event_time <- stats::rexp(n, exp(1.1 * x$gene1 - 0.7 * x$gene2))
    censor_time <- stats::rexp(n, 0.2)
    list(x = x, time = pmin(event_time, censor_time),
      event = as.integer(event_time <= censor_time))
  })
}

test_that("StepCox tightens likelihood precision without changing warning tolerance", {
  d <- make_cox_control_fixture()
  fit <- sb_fit_stepcox(d$x, d$time, d$event, "both")
  native <- survival::coxph.control()
  expect_equal(fit$cox_control$eps, 1e-11)
  expect_identical(fit$cox_control$iter.max, 50L)
  expect_identical(fit$cox_control$toler.inf, native$toler.inf)
  expect_identical(fit$cox_control$toler.chol, native$toler.chol)
  expect_gt(length(fit$selected_features), 0L)
})

test_that("PLS-Cox preserves the native infinite-coefficient warning threshold", {
  skip_if_not_installed("plsRcox")
  d <- make_cox_control_fixture()
  fit <- sb_plsrcox_cv_learner_adapter(inner_folds = 3L, seed = 1L,
    maximum_components = 3L)$fit(d$x, d$time, d$event)
  native <- survival::coxph.control()
  expect_equal(fit$audit$cox_control$eps, 1e-11)
  expect_identical(fit$audit$cox_control$iter.max, 50L)
  expect_identical(fit$audit$cox_control$toler.inf, native$toler.inf)
  expect_identical(fit$audit$cox_control$toler.chol, native$toler.chol)
})

test_that("StepCox still rejects complete survival separation", {
  x <- data.frame(separated = rep(c(1, 0), each = 20L))
  error <- tryCatch(sb_fit_stepcox(x, seq_len(40L), rep(1L, 40L), "both"),
    error = identity)
  expect_s3_class(error, "stormbreaker_stepcox_warning")
  expect_match(conditionMessage(error), "coefficient may be infinite", fixed = TRUE)
})

# ---- coxphf-bagged-knn-wsrf ----

make_native_expansion_binary_fixture <- function(
  n = 90L, p = 5L, seed = 44101L
) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    eta <- 0.8 * x[, 1L] - 0.6 * x[, min(2L, p)]
    outcome <- stats::rbinom(n, 1L, stats::plogis(eta))
    data.frame(
      patient_id = sprintf("NB%03d", seq_len(n)),
      outcome = outcome,
      x,
      check.names = FALSE
    ) |>
      stats::setNames(c(
        "patient_id", "outcome", "gene one",
        if (p >= 2L) "2gene",
        if (p > 2L) paste0("gene-", seq_len(p - 2L))
      ))
  })
}

make_native_expansion_survival_fixture <- function(
  n = 100L, p = 5L, seed = 44201L
) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    eta <- 0.7 * x[, 1L] - 0.5 * x[, min(2L, p)]
    event_time <- stats::rexp(n, rate = exp(eta))
    censor_time <- stats::rexp(n, rate = 0.45)
    data.frame(
      patient_id = sprintf("NS%03d", seq_len(n)),
      time = pmin(event_time, censor_time) + 0.02,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    ) |>
      stats::setNames(c(
        "patient_id", "time", "event", "gene one",
        if (p >= 2L) "2gene",
        if (p > 2L) paste0("gene-", seq_len(p - 2L))
      ))
  })
}

test_that("three native learners expose their declared selector routes", {
  ids <- c(
    sb_wsrf_binary_learner_ids(),
    sb_coxphf_learner_ids(),
    sb_bagged_knn_survival_learner_ids()
  )
  binary_selectors <- expected_fixed_learner_selectors("binary")
  survival_selectors <- expected_fixed_learner_selectors("survival")
  registry <- sb_pipeline_registry()
  rows <- registry[registry$learner %in% ids, , drop = FALSE]
  expect_identical(
    nrow(rows),
    length(binary_selectors) + 2L * length(survival_selectors)
  )
  expect_setequal(
    rows$selector[rows$learner == sb_wsrf_binary_learner_ids()],
    binary_selectors
  )
  for (learner in c(
    sb_coxphf_learner_ids(),
    sb_bagged_knn_survival_learner_ids()
  )) {
    expect_setequal(
      rows$selector[rows$learner == learner],
      survival_selectors
    )
  }
  expect_route_dependencies(rows, c(wsrf_binary = "wsrf",
    coxphf_firth = "coxphf", bagged_knn_survival = "bnnSurvival"))

  for (index in seq_len(nrow(rows))) {
    spec <- list(
      selector = rows$selector[[index]],
      learner = rows$learner[[index]]
    )
    if (rows$selector[[index]] %in% c(
      "univ_binary_auc_top_k", "univ_cox_top_k"
    )) {
      spec$selector_n_features <- 3L
    }
    if (identical(rows$selector[[index]], "stepcox")) {
      spec$stepcox_direction <- "both"
    }
    if (identical(rows$selector[[index]], "fixed_features")) {
      spec$selector_features <- c("gene one", "2gene", "gene-1")
      spec$feature_set_id <- "native_expansion_test"
    }
    if (spec$selector %in% c("lasso_nonzero", "enet_cox_nonzero", "coxboost_nonzero")) {
      spec$selector_inner_folds <- 3L
      spec$selector_seed <- 1L
    }
    if (spec$selector == "rsf_minimal_depth") {
      spec$selector_rsf_nodesize <- 5L
      spec$selector_seed <- 1L
    }
    expect_identical(sb_lint_pipeline(spec)$status, "admitted")
  }
})

test_that("Firth Cox matches native coefficients and portable x-beta replay", {
  skip_if_not_installed("coxphf", minimum_version = "1.13.4")
  data <- make_native_expansion_survival_fixture()
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[c(2L, 19L, 73L), , drop = FALSE]
  adapter <- sb_coxphf_learner_adapter("coxphf_firth")

  model <- adapter$fit(x, data$time, data$event)
  prepared <- sb_coxphf_formula_data(list(
    x = x,
    time = data$time,
    event = data$event
  ))
  native <- coxphf::coxphf(
    prepared$formula,
    prepared$data,
    pl = FALSE,
    maxit = 100L,
    maxhs = 5L,
    epsilon = 1e-6,
    gconv = 1e-4,
    maxstep = 0.5,
    firth = TRUE,
    penalty = 0.5
  )
  expect_identical(model$coefficient, unname(native$coefficients))
  risk <- adapter$predict(model, assessment)
  expect_identical(
    risk,
    unname(as.vector(as.matrix(assessment) %*% native$coefficients))
  )

  horizon <- as.numeric(stats::median(data$time))
  probability <- adapter$predict_survival(
    model, x, data$time, data$event, assessment, horizon
  )
  state <- adapter$extract_state(model, list(
    x = x,
    time = data$time,
    event = data$event,
    prediction_horizon = horizon,
    time_unit = "month"
  ))
  replay <- adapter$predict_state(
    unserialize(serialize(state, NULL)),
    assessment
  )
  expect_identical(replay$risk, risk)
  expect_identical(replay$survival_probability, probability$probability)
  expect_identical(probability$ties, "breslow")
  expect_false(any(c("engine", "time", "event", "outcome") %in% names(state)))
  expect_error(
    adapter$predict(model, assessment[rev(names(assessment))]),
    class = "stormbreaker_learner_schema_mismatch"
  )
})

test_that("Firth Cox handles monotone likelihood and exposes warnings", {
  skip_if_not_installed("coxphf", minimum_version = "1.13.4")
  n <- 30L
  x <- data.frame(marker = rev(seq_len(n)))
  adapter <- sb_coxphf_learner_adapter("coxphf_firth")
  model <- adapter$fit(x, seq_len(n) + 0.01, rep(1L, n))
  expect_true(all(is.finite(model$coefficient)))
  expect_true(model$iterations < adapter$params$max_iterations)
  expect_error(
    sb_coxphf_stop_on_warning(
      {
        warning("native convergence warning")
        TRUE
      },
      "fitting"
    ),
    class = "stormbreaker_coxphf_warning"
  )
})

test_that("bagged kNN survival matches native curves and retains required data", {
  skip_if_not_installed("bnnSurvival", minimum_version = "0.1.5")
  data <- make_native_expansion_survival_fixture(seed = 44301L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  rownames(x) <- paste0("PRIVATE_PATIENT_", seq_len(nrow(x)))
  assessment <- x[c(4L, 28L, 77L), , drop = FALSE]
  adapter <- sb_bagged_knn_survival_learner_adapter(
    "bagged_knn_survival"
  )
  old_cores <- getOption("mc.cores")
  set.seed(44302L)
  caller_rng <- .Random.seed
  model <- adapter$fit(x, data$time, data$event)
  expect_identical(.Random.seed, caller_rng)
  expect_identical(getOption("mc.cores"), old_cores)

  prepared <- sb_bagged_knn_formula_data(list(
    x = x,
    time = data$time,
    event = data$event
  ))
  native <- sb_with_bagged_knn_single_core(withr::with_seed(
    1L,
    bnnSurvival::bnnSurvival(
      prepared$formula,
      prepared$analysis,
      k = max(1L, floor(nrow(x) / 10L)),
      num_base_learners = 50L,
      num_features_per_base_learner = length(prepared$internal_features),
      metric = "mahalanobis",
      weighting_function = sb_bagged_knn_constant_weight,
      replace = TRUE,
      sample_fraction = 1
    )
  ))
  predict_generic <- get(
    "predict",
    envir = asNamespace("bnnSurvival"),
    inherits = FALSE
  )
  native_result <- sb_with_bagged_knn_single_core(
    predict_generic(native, sb_bagged_knn_internal_data(assessment)$data)
  )
  curves <- sb_predict_bagged_knn_survival_curves(model, assessment)
  expect_identical(
    curves$survival,
    as.matrix(bnnSurvival::predictions(native_result))
  )
  expect_true(all(apply(curves$survival, 1L, function(values) {
    all(diff(values) <= 1e-10)
  })))

  horizon <- as.numeric(stats::median(data$time))
  probability <- adapter$predict_survival(
    model, x, data$time, data$event, assessment, horizon
  )
  state <- adapter$extract_state(model, list(
    x = x,
    time = data$time,
    event = data$event,
    prediction_horizon = horizon,
    time_unit = "month"
  ))
  replay <- adapter$predict_state(
    unserialize(serialize(state, NULL)),
    assessment
  )
  expect_identical(replay$risk, curves$risk)
  expect_identical(replay$survival_probability, probability$probability)
  expect_identical(
    state$training_data_retention,
    "required_by_bagged_knn_prediction_state"
  )
  restored_engine <- unserialize(state$engine_raw)
  training_data <- methods::slot(restored_engine, "train_data")
  expect_gt(nrow(training_data), 0L)
  expect_false(any(startsWith(
    rownames(training_data),
    "PRIVATE_PATIENT_"
  )))
  expect_length(
    grepRaw("PRIVATE_PATIENT_", state$engine_raw, fixed = TRUE, all = TRUE),
    0L
  )
  expect_identical(
    environment(methods::slot(restored_engine, "weighting_function")),
    environment(sb_bagged_knn_constant_weight)
  )

  one_feature <- x[1L]
  one_model <- adapter$fit(one_feature, data$time, data$event)
  expect_true(one_model$uses_zero_variance_bridge)
  expect_true(all(is.finite(adapter$predict(
    one_model,
    one_feature[1:2, , drop = FALSE]
  ))))
})

test_that("weighted-subspace forest matches native votes and compact replay", {
  skip_if_not_installed("wsrf", minimum_version = "1.7.32")
  data <- make_native_expansion_binary_fixture(seed = 44401L)
  features <- setdiff(names(data), c("patient_id", "outcome"))
  x <- data[features]
  assessment <- x[c(3L, 31L, 82L), , drop = FALSE]
  adapter <- sb_wsrf_binary_learner_adapter("wsrf_binary")

  set.seed(44402L)
  caller_rng <- .Random.seed
  model <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, caller_rng)
  prepared <- sb_wsrf_internal_data(x)
  native <- withr::with_seed(
    1L,
    wsrf::wsrf(
      prepared$data,
      factor(
        data$outcome,
        levels = c(0L, 1L),
        labels = c("control", "case")
      ),
      mtry = floor(log2(ncol(x)) + 1),
      ntree = 500L,
      weights = TRUE,
      parallel = FALSE,
      importance = FALSE,
      nodesize = 2L
    )
  )
  native_probability <- stats::predict(
    native,
    sb_wsrf_internal_data(assessment)$data,
    type = "prob"
  )$prob[, "case"]
  probability <- adapter$predict(model, assessment)
  expect_identical(probability, unname(native_probability))
  expect_null(model$engine$targetData)
  expect_true(all(lengths(model$engine$OOBSets) == 0L))
  expect_true(all(lengths(model$engine$OOBPredictLabels) == 0L))

  state <- adapter$extract_state(model)
  replay <- adapter$predict_state(
    unserialize(serialize(state, NULL)),
    assessment
  )
  expect_identical(replay$probability, probability)
  expect_false(state$training_outcome_retained)
  expect_false(state$oob_patient_indices_retained)
  expect_error(
    adapter$predict(model, assessment[rev(names(assessment))]),
    class = "stormbreaker_learner_schema_mismatch"
  )

  one_feature <- x[1L]
  one_model <- adapter$fit(one_feature, data$outcome)
  expect_true(one_model$uses_zero_variance_bridge)
  expect_true(all(adapter$predict(
    one_model,
    one_feature[1:2, , drop = FALSE]
  ) >= 0))
})

test_that("new native dependency failures remain explicit", {
  local_mocked_bindings(sb_coxphf_available = function() FALSE)
  expect_error(
    sb_require_coxphf_dependency(),
    class = "stormbreaker_missing_coxphf"
  )
  local_mocked_bindings(sb_bagged_knn_survival_available = function() FALSE)
  expect_error(
    sb_require_bagged_knn_survival_dependency(),
    class = "stormbreaker_missing_bnn_survival"
  )
  local_mocked_bindings(sb_wsrf_available = function() FALSE)
  expect_error(sb_require_wsrf_dependency(), class = "stormbreaker_missing_wsrf")
})

test_that("three native learners execute CV and outcome-blind final prediction", {
  skip_if_not_installed("coxphf", minimum_version = "1.13.4")
  skip_if_not_installed("bnnSurvival", minimum_version = "0.1.5")
  skip_if_not_installed("wsrf", minimum_version = "1.7.32")

  binary <- make_native_expansion_binary_fixture(n = 84L, seed = 44501L)
  binary_features <- setdiff(names(binary), c("patient_id", "outcome"))
  binary_task <- sb_binary_task(
    binary,
    "patient_id",
    "outcome",
    binary_features
  )
  binary_pipeline <- sb_pipeline(learner = "wsrf_binary")
  binary_cv <- sb_run_cv(binary_task, binary_pipeline, folds = 3L, seed = 1L)
  binary_fit <- sb_fit_final(binary_task, binary_pipeline)
  binary_prediction <- sb_predict_final(
    unserialize(serialize(binary_fit, NULL)),
    binary[1:5, c("patient_id", rev(binary_features)), drop = FALSE],
    "patient_id"
  )
  expect_identical(binary_cv$status, "success")
  expect_true(binary_cv$complete_oof)
  expect_identical(
    binary_cv$provenance$wsrf_version,
    as.character(utils::packageVersion("wsrf"))
  )
  expect_identical(binary_fit$status, "success")
  expect_identical(binary_prediction$status, "success")

  survival <- make_native_expansion_survival_fixture(
    n = 90L,
    seed = 44502L
  )
  survival_features <- setdiff(
    names(survival),
    c("patient_id", "time", "event")
  )
  horizon <- as.numeric(stats::quantile(survival$time, 0.45))
  survival_task <- sb_surv_task(
    survival,
    "patient_id",
    "time",
    "event",
    survival_features,
    "month",
    prediction_horizon = horizon
  )
  for (learner in c(
    "coxphf_firth",
    "bagged_knn_survival"
  )) {
    pipeline <- sb_pipeline(learner = learner)
    cv <- sb_run_cv(survival_task, pipeline, folds = 3L, seed = 2L)
    fit <- sb_fit_final(survival_task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      survival[1:5, c("patient_id", rev(survival_features)), drop = FALSE],
      "patient_id"
    )
    expect_identical(cv$status, "success")
    expect_true(cv$complete_oof)
    expect_true(cv$complete_survival_probability_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(
      prediction$predictions$survival_probability >= 0 &
        prediction$predictions$survival_probability <= 1
    ))
  }
})

test_that("bagged kNN ordinary warnings preserve native fit and predictions", {
  skip_if_not_installed("bnnSurvival")
  d <- make_native_expansion_survival_fixture(n = 40L, p = 3L)
  x <- d[setdiff(names(d), c("patient_id", "time", "event"))]
  id <- "bagged_knn_survival"
  params <- sb_bagged_knn_survival_params(id)
  params$num_base_learners <- 5L
  adapter <- sb_bagged_knn_survival_learner_adapter(id, params)
  reference <- adapter$fit(x, d$time, d$event)
  expected <- adapter$predict(reference, x)
  native_fit <- bnnSurvival::bnnSurvival
  native_predict <- sb_predict_bagged_knn_engine
  local_mocked_bindings(bnnSurvival = function(...) {
    warning("informational fit warning")
    native_fit(...)
  }, .package = "bnnSurvival")
  local_mocked_bindings(sb_predict_bagged_knn_engine = function(...) { warning("informational prediction warning"); native_predict(...) })
  seen <- character()
  observed <- withCallingHandlers({
    fitted <- adapter$fit(x, d$time, d$event)
    adapter$predict(fitted, x)
  }, warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(observed, expected)
  expect_identical(all(c("informational fit warning", "informational prediction warning") %in% seen), TRUE)
})

test_that("WSRF ordinary warnings preserve native fit and predictions", {
  skip_if_not_installed("wsrf")
  d <- make_native_expansion_binary_fixture(n = 40L, p = 3L)
  x <- d[setdiff(names(d), c("patient_id", "outcome"))]
  id <- "wsrf_binary"
  params <- sb_wsrf_binary_params(id)
  params$num_trees <- 5L
  adapter <- sb_wsrf_binary_learner_adapter(id, params)
  reference <- adapter$fit(x, d$outcome)
  expected <- adapter$predict(reference, x)
  native_fit <- wsrf::wsrf
  native_predict <- stats::predict
  local_mocked_bindings(wsrf = function(...) {
    warning("informational fit warning")
    native_fit(...)
  }, .package = "wsrf")
  local_mocked_bindings(predict = function(...) { warning("informational prediction warning"); native_predict(...) }, .package = "stats")
  seen <- character()
  observed <- withCallingHandlers({
    fitted <- adapter$fit(x, d$outcome)
    adapter$predict(fitted, x)
  }, warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(observed, expected)
  expect_identical(all(c("informational fit warning", "informational prediction warning") %in% seen), TRUE)
})

# ---- rlt-adapters ----

make_rlt_binary_fixture <- function(n = 100L, p = 6L, seed = 43101L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    eta <- x[, 1L] - 0.7 * x[, min(2L, p)]
    outcome <- stats::rbinom(n, 1L, stats::plogis(eta))
    feature_names <- c(
      "gene one",
      if (p >= 2L) "2gene",
      if (p > 2L) paste0("gene-", seq_len(p - 2L))
    )
    data.frame(
      patient_id = sprintf("RB%03d", seq_len(n)),
      outcome = outcome,
      x,
      check.names = FALSE
    ) |>
      stats::setNames(c("patient_id", "outcome", feature_names))
  })
}

make_rlt_survival_fixture <- function(n = 120L, p = 6L, seed = 43201L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    eta <- 0.7 * x[, 1L] - 0.5 * x[, min(2L, p)]
    event_time <- stats::rexp(n, rate = exp(eta))
    censor_time <- stats::rexp(n, rate = 0.45)
    time <- pmin(event_time, censor_time) + 0.01
    feature_names <- c(
      "gene one",
      if (p >= 2L) "2gene",
      if (p > 2L) paste0("gene-", seq_len(p - 2L))
    )
    data.frame(
      patient_id = sprintf("RS%03d", seq_len(n)),
      time = time,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    ) |>
      stats::setNames(c("patient_id", "time", "event", feature_names))
  })
}

rlt_native_fit <- function(x, endpoint, params) {
  matrix_x <- as.matrix(x)
  storage.mode(matrix_x) <- "double"
  arguments <- list(
    x = matrix_x,
    model = params$model,
    print.summary = 0L,
    use.cores = params$cores,
    ntrees = params$num_trees,
    nmin = params$native_nmin,
    reinforcement = params$reinforcement,
    combsplit = params$combination_split_size,
    embed.n.th = params$embedded_stop_size,
    importance = params$importance,
    track.obs = params$track_observations
  )
  if (identical(params$model, "classification")) {
    arguments$y <- factor(
      endpoint$outcome,
      levels = c(0L, 1L),
      labels = c("control", "case")
    )
  } else {
    arguments$y <- endpoint$time
    arguments$censor <- endpoint$event
  }
  withr::with_seed(params$seed, do.call(RLT::RLT, arguments))
}

test_that("two RLT reinforcement learners expose their legal routes", {
  expect_identical(
    sb_rlt_learner_ids(),
    c(
      "rlt_reinforcement_binary",
      "rlt_reinforcement_survival"
    )
  )
  for (learner in sb_rlt_learner_ids()) {
    params <- sb_rlt_learner_params(learner)
    expect_identical(params$num_trees, 100L)
    expect_identical(params$reinforcement, TRUE)
    expect_identical(params$combination_split_size, 1L)
    expect_identical(params$native_nmin, 5L)
    expect_identical(params$embedded_stop_size, 1L)
    expect_identical(params$cores, 1L)
    expect_identical(params$seed, 1L)
  }

  binary_selectors <- expected_fixed_learner_selectors("binary")
  survival_selectors <- expected_fixed_learner_selectors("survival")
  registry <- sb_pipeline_registry()
  rows <- registry[registry$learner %in% sb_rlt_learner_ids(), , drop = FALSE]
  expect_identical(
    nrow(rows),
    length(binary_selectors) + length(survival_selectors)
  )
  expect_true(all(rows$optional_dependency == "RLT" |
    grepl("RLT", rows$optional_dependency, fixed = TRUE)))
  expect_true(all(rows$tuning_semantics == "none"))
  expect_setequal(
    rows$selector[rows$task_type == "binary"],
    binary_selectors
  )
  expect_setequal(
    rows$selector[rows$task_type == "survival"],
    survival_selectors
  )

  for (index in seq_len(nrow(rows))) {
    spec <- list(
      selector = rows$selector[[index]],
      learner = rows$learner[[index]]
    )
    if (rows$selector[[index]] %in% c(
      "univ_binary_auc_top_k", "univ_cox_top_k"
    )) {
      spec$selector_n_features <- 3L
    }
    if (identical(rows$selector[[index]], "stepcox")) {
      spec$stepcox_direction <- "both"
    }
    if (identical(rows$selector[[index]], "fixed_features")) {
      spec$selector_features <- c("gene one", "2gene", "gene-1")
      spec$feature_set_id <- "rlt_test"
    }
    if (spec$selector %in% c("lasso_nonzero", "enet_cox_nonzero", "coxboost_nonzero")) {
      spec$selector_inner_folds <- 3L
      spec$selector_seed <- 1L
    }
    if (spec$selector == "rsf_minimal_depth") {
      spec$selector_rsf_nodesize <- 5L
      spec$selector_seed <- 1L
    }
    expect_identical(sb_lint_pipeline(spec)$status, "admitted")
  }
})

test_that("binary RLT matches native probabilities and portable replay", {
  skip_if_not_installed("RLT", minimum_version = "3.2.6")

  data <- make_rlt_binary_fixture()
  features <- setdiff(names(data), c("patient_id", "outcome"))
  x <- data[features]
  assessment <- x[c(2L, 19L, 61L), , drop = FALSE]
  adapter <- sb_rlt_learner_adapter("rlt_reinforcement_binary")

  set.seed(43102L)
  caller_rng <- .Random.seed
  model <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, caller_rng)
  native <- rlt_native_fit(
    x,
    list(outcome = data$outcome),
    adapter$params
  )
  native_probability <- as.matrix(
    stats::predict(native, assessment)$ProbPrediction
  )[, 2L]
  probability <- adapter$predict(model, assessment)
  expect_identical(probability, unname(native_probability))
  expect_true(all(probability >= 0 & probability <= 1))
  expect_identical(model$engine$ylevels, c("control", "case"))
  expect_false(any(c(
    "AllError", "VarImp", "ObsTrack", "n", "subjectnames",
    "subject.weight", "x", "y", "outcome"
  ) %in% names(model$engine)))

  state <- unserialize(serialize(adapter$extract_state(model), NULL))
  compact <- unserialize(state$engine_raw)
  replay <- adapter$predict_state(state, assessment)
  expect_identical(replay$probability, probability)
  expect_false(any(c(
    "AllError", "VarImp", "ObsTrack", "n", "subjectnames",
    "subject.weight", "x", "y", "outcome"
  ) %in% names(compact)))
  expect_error(
    adapter$predict(model, assessment[rev(names(assessment))]),
    class = "stormbreaker_learner_schema_mismatch"
  )
})

test_that("survival RLT reconstructs native hazards and direct survival", {
  skip_if_not_installed("RLT", minimum_version = "3.2.6")

  data <- make_rlt_survival_fixture()
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[c(3L, 27L, 92L), , drop = FALSE]
  adapter <- sb_rlt_learner_adapter("rlt_reinforcement_survival")

  set.seed(43202L)
  caller_rng <- .Random.seed
  model <- adapter$fit(x, data$time, data$event)
  expect_identical(.Random.seed, caller_rng)
  native <- rlt_native_fit(
    x,
    list(time = data$time, event = data$event),
    adapter$params
  )
  native_prediction <- stats::predict(native, assessment)
  hazard <- as.matrix(native_prediction$SurvPred)
  survival <- matrix(NA_real_, nrow(hazard), ncol(hazard))
  for (row in seq_len(nrow(hazard))) {
    survival[row, ] <- cumprod(1 - hazard[row, ])
  }
  interval <- diff(c(0, native_prediction$timepoints))
  previous <- cbind(1, survival[, -ncol(survival), drop = FALSE])
  native_risk <- -as.vector(previous %*% interval)
  curves <- sb_predict_rlt_survival_curves(model, assessment)
  expect_identical(curves$hazard, hazard)
  expect_identical(curves$survival, survival)
  expect_identical(curves$risk, unname(native_risk))
  expect_true(all(apply(
    curves$survival,
    1L,
    function(values) all(diff(values) <= 0)
  )))

  early_horizon <- min(curves$time) / 2
  early <- adapter$predict_survival(
    model, x, data$time, data$event, assessment, early_horizon
  )
  later_horizon <- stats::median(curves$time)
  later <- adapter$predict_survival(
    model, x, data$time, data$event, assessment, later_horizon
  )
  expect_identical(early$probability, rep(1, nrow(assessment)))
  expect_true(all(later$probability >= 0 & later$probability <= 1))
  expect_identical(later$source, sb_rlt_direct_survival_source())
  expect_identical(later$ties, "not_applicable")

  state <- unserialize(serialize(adapter$extract_state(model, list(
    x = x,
    time = data$time,
    event = data$event,
    prediction_horizon = later_horizon,
    time_unit = "month"
  )), NULL))
  replay <- adapter$predict_state(state, assessment)
  expect_identical(replay$risk, curves$risk)
  expect_identical(replay$survival_probability, later$probability)
  expect_identical(
    state$risk_formula,
    "negative_step_RMST_over_native_event_grid"
  )
  compact <- unserialize(state$engine_raw)
  expect_false(any(c(
    "AllError", "VarImp", "ObsTrack", "n", "subjectnames",
    "subject.weight", "x", "y", "outcome", "time", "event"
  ) %in% names(compact)))
})

test_that("RLT native-domain failures remain typed and do not enter C code", {
  skip_if_not_installed("RLT", minimum_version = "3.2.6")

  data <- make_rlt_survival_fixture(n = 30L, p = 2L, seed = 43301L)
  x <- data[setdiff(names(data), c("patient_id", "time", "event"))]
  adapter <- sb_rlt_learner_adapter("rlt_reinforcement_survival")
  event <- integer(nrow(data))
  event[[1L]] <- 1L
  expect_error(
    adapter$fit(x, data$time, event),
    class = "stormbreaker_rlt_insufficient_event_times"
  )
  event[[2L]] <- 1L
  tied_time <- data$time
  tied_time[1:2] <- tied_time[[1L]]
  expect_error(
    adapter$fit(x, tied_time, event),
    class = "stormbreaker_rlt_insufficient_event_times"
  )
  zero_time <- data$time
  zero_time[[1L]] <- 0
  expect_error(
    adapter$fit(x, zero_time, data$event),
    class = "stormbreaker_bad_rlt_time"
  )
  expect_error(
    sb_rlt_stop_on_warning(
      {
        warning("Subject weights length must be n, reset to equal weights")
        TRUE
      },
      "rlt_reinforcement_binary",
      "fitting"
    ),
    class = "stormbreaker_rlt_warning"
  )
})

test_that("RLT dependency failure is explicit", {
  local_mocked_bindings(sb_rlt_available = function() FALSE)
  expect_error(
    sb_require_rlt_dependency(),
    class = "stormbreaker_missing_rlt"
  )
})

test_that("RLT executes through CV, final fit, and outcome-blind prediction", {
  skip_if_not_installed("RLT", minimum_version = "3.2.6")

  binary <- make_rlt_binary_fixture(n = 90L, seed = 43401L)
  binary_features <- setdiff(names(binary), c("patient_id", "outcome"))
  binary_task <- sb_binary_task(
    binary,
    "patient_id",
    "outcome",
    binary_features
  )
  binary_pipeline <- sb_pipeline(
    learner = "rlt_reinforcement_binary"
  )
  binary_cv <- sb_run_cv(
    binary_task,
    binary_pipeline,
    folds = 3L,
    seed = 43402L
  )
  binary_fit <- sb_fit_final(binary_task, binary_pipeline)
  binary_prediction <- sb_predict_final(
    unserialize(serialize(binary_fit, NULL)),
    binary[1:6, c("patient_id", rev(binary_features)), drop = FALSE],
    "patient_id"
  )
  rlt_version <- as.character(utils::packageVersion("RLT"))
  expect_identical(binary_cv$status, "success")
  expect_true(binary_cv$complete_oof)
  expect_identical(binary_cv$provenance$rlt_version, rlt_version)
  expect_identical(binary_fit$status, "success")
  expect_identical(binary_fit$provenance$rlt_version, rlt_version)
  expect_identical(binary_prediction$status, "success")

  survival <- make_rlt_survival_fixture(n = 105L, seed = 43403L)
  survival_features <- setdiff(
    names(survival),
    c("patient_id", "time", "event")
  )
  horizon <- as.numeric(stats::quantile(survival$time, 0.45))
  survival_task <- sb_surv_task(
    survival,
    "patient_id",
    "time",
    "event",
    survival_features,
    "month",
    prediction_horizon = horizon
  )
  survival_pipeline <- sb_pipeline(
    learner = "rlt_reinforcement_survival"
  )
  survival_cv <- sb_run_cv(
    survival_task,
    survival_pipeline,
    folds = 3L,
    seed = 43404L
  )
  survival_fit <- sb_fit_final(survival_task, survival_pipeline)
  survival_prediction <- sb_predict_final(
    unserialize(serialize(survival_fit, NULL)),
    survival[1:6, c("patient_id", rev(survival_features)), drop = FALSE],
    "patient_id"
  )
  expect_identical(survival_cv$status, "success")
  expect_true(survival_cv$complete_oof)
  expect_true(survival_cv$complete_survival_probability_oof)
  expect_identical(survival_cv$provenance$rlt_version, rlt_version)
  expect_identical(survival_fit$status, "success")
  expect_identical(survival_fit$provenance$rlt_version, rlt_version)
  expect_identical(survival_prediction$status, "success")
  expect_true(all(is.finite(survival_prediction$predictions$risk)))
  expect_true(all(
    survival_prediction$predictions$survival_probability >= 0 &
      survival_prediction$predictions$survival_probability <= 1
  ))
})

test_that("RLT ordinary warnings preserve native fit and predictions", {
  skip_if_not_installed("RLT")
  d <- make_rlt_binary_fixture(n = 40L, p = 3L)
  x <- d[setdiff(names(d), c("patient_id", "outcome"))]
  id <- "rlt_reinforcement_binary"
  params <- sb_rlt_learner_params(id)
  params$num_trees <- 5L
  adapter <- sb_rlt_learner_adapter(id, params)
  reference <- adapter$fit(x, d$outcome)
  expected <- adapter$predict(reference, x)
  native_fit <- RLT::RLT
  native_predict <- get("predict.RLT", asNamespace("RLT"))
  local_mocked_bindings(RLT = function(...) {
    warning("informational fit warning")
    native_fit(...)
  }, .package = "RLT")
  local_mocked_bindings(predict.RLT = function(...) {
    warning("informational prediction warning")
    native_predict(...)
  }, .package = "RLT")
  seen <- character()
  observed <- withCallingHandlers({
    fitted <- adapter$fit(x, d$outcome)
    adapter$predict(fitted, x)
  }, warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(observed, expected)
  expect_identical(all(c("informational fit warning", "informational prediction warning") %in% seen), TRUE)
})

test_that("RLT retains weight reset and feature mismatch protection", {
  messages <- c(
    "Subject weights length must be n, reset to equal weights",
    "Subject weights cannot be 0 or negative, reset to equal weights",
    "Variable weights length must by p, reset to equal weights",
    "Variable weights cannot be negative, reset to equal weights",
    "test data variables names does not match training data...")
  for (message in messages) {
    error <- tryCatch(sb_rlt_stop_on_warning({ warning(message); TRUE },
      "rlt_reinforcement_binary", "fitting"), error = identity)
    expect_s3_class(error, "stormbreaker_rlt_warning")
  }
})
})

# ---- model-ranger ----
local({
# ---- algorithm-forest ----

test_that("all seven ranger variants expose exact saved importance and tree counts", {
  skip_if_not_installed("ranger")
  for (id in sb_ranger_learner_ids()) {
    binary <- id %in% sb_ranger_binary_learner_ids()
    d <- if (binary) make_sb_binary_fixture(n = 100L, seed = 108L) else make_sb_fixture(n = 100L, seed = 8612L)
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    expect_null(fit$diagnostic_artifacts)
    original <- fit
    q <- sb_forest_diagnostics(fit)
    expect_identical(q$table, fit$learner_state$feature_importance)
    sizes <- sb_forest_diagnostics(fit, "tree_size")
    forest <- unserialize(fit$learner_state$model_raw)
    native_object <- structure(list(forest = forest), class = "ranger")
    native <- vapply(seq_len(forest$num.trees), function(i) {
      nrow(subset(ranger::treeInfo(native_object, tree = i), terminal))
    }, integer(1))
    expect_equal(sizes$table$terminal_nodes, native)
    expect_equal(nrow(sizes$table), fit$learner_state$num_trees)
    expect_identical(fit, original)
    err <- tryCatch(sb_forest_diagnostics(fit, "oob_error"), error = identity)
    expect_s3_class(err, "stormbreaker_unavailable_diagnostic")
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      for (z in list(q, sizes)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
      expect_equal(sum(sb_plot_algorithm_diagnostic(sizes)$data$trees), forest$num.trees)
    }
  }
})

test_that("randomForest and randomForestSRC diagnostics preserve native OOB missingness", {
  skip_if_not_installed("randomForest")
  skip_if_not_installed("randomForestSRC")
  for (id in c("random_forest_binary", "rsf_survival")) {
    binary <- id == "random_forest_binary"
    d <- if (binary) make_sb_binary_fixture(n = 100L, seed = 108L) else make_sb_fixture(n = 100L, seed = 8612L)
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    original <- fit
    engine <- fit$diagnostic_artifacts$learner$engine
    sizes <- sb_forest_diagnostics(fit, "tree_size")
    native <- if (binary) (engine$forest$ndbigtree + 1L) / 2L else engine$leaf.count
    expect_equal(sizes$table$terminal_nodes, native)
    q <- sb_forest_diagnostics(fit, "oob_error")
    expect_equal(q$table$error, as.vector(engine$err.rate))
    expect_equal(q$table$available, is.finite(as.vector(engine$err.rate)))
    expect_identical(fit, original)
    if (binary) {
      importance <- sb_forest_diagnostics(fit)
      expect_equal(importance$table$importance, unname(engine$importance[, "MeanDecreaseGini"]))
    } else {
      err <- tryCatch(sb_forest_diagnostics(fit), error = identity)
      expect_s3_class(err, "stormbreaker_unavailable_diagnostic")
    }
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      for (z in list(q, sizes)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
      expect_equal(nrow(sb_plot_algorithm_diagnostic(q)$data), sum(q$table$available))
    }
    fit$diagnostic_artifacts <- NULL
    expect_identical(sb_forest_diagnostics(fit, "tree_size"), sizes)
    err <- tryCatch(sb_forest_diagnostics(fit, "oob_error"), error = identity)
    expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
  }
})

# ---- ranger-binary-fixed-adapter ----

test_that("fixed ranger binary state reconstructs class-one probabilities", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  d <- make_sb_binary_fixture(n = 54L, seed = 7701L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(4L, 19L, 37L), , drop = FALSE]
  adapter <- sb_ranger_binary_learner_adapter(
    num_trees = 18L,
    mtry = NULL,
    min_node_size = 10L,
    seed = 7702L
  )
  model <- adapter$fit(x, d$outcome)
  observed <- adapter$predict(model, new_x)
  audit <- adapter$extract_tree_model(model)
  restored <- unserialize(audit$model_raw)
  oracle <- stats::predict(
    restored,
    data = new_x,
    num.trees = 18L,
    type = "response",
    seed = 0L,
    num.threads = 1L,
    verbose = FALSE
  )$predictions
  expected <- unname(oracle[, "1"])
  state <- adapter$extract_state(model)

  expect_s3_class(model, "stormbreaker_ranger_binary")
  expect_identical(model$effective_mtry, 1L)
  expect_equal(observed, expected, tolerance = 0)
  expect_identical(adapter$params$probability, TRUE)
  expect_identical(adapter$params$splitrule, "gini")
  expect_identical(adapter$params$min_node_size, 10L)
  expect_identical(audit$num_trees, 18L)
  expect_identical(audit$effective_mtry, 1L)
  expect_identical(audit$features, names(x))
  expect_identical(audit$params, adapter$params)
  expect_identical(audit$importance_type, "impurity")
  expect_identical(audit$model_serialization, "R_xdr_v3_ranger.forest")
  expect_type(audit$model_raw, "raw")
  expect_gt(length(audit$model_raw), 0L)
  expect_identical(audit$feature_importance$feature, names(x))
  expect_true(all(is.finite(audit$feature_importance$importance)))
  expect_identical(state$learner, "ranger_binary")
  expect_s3_class(state, "stormbreaker_learner_state")
  expect_false("forest" %in% names(state))
  expect_equal(
    adapter$predict_state(state, new_x)$probability,
    expected,
    tolerance = 0
  )
})

test_that("fixed ranger binary is deterministic and fails explicitly", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  d <- make_sb_binary_fixture(n = 50L, seed = 7711L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_ranger_binary_learner_adapter(
    num_trees = 15L,
    mtry = 2L,
    min_node_size = 10L,
    seed = 7712L
  )
  first <- adapter$fit(x[1:40, , drop = FALSE], d$outcome[1:40])
  second <- adapter$fit(x[1:40, , drop = FALSE], d$outcome[1:40])
  new_x <- x[41:50, , drop = FALSE]
  batch <- adapter$predict(first, new_x)
  one_at_a_time <- vapply(seq_len(nrow(new_x)), function(row) {
    adapter$predict(first, new_x[row, , drop = FALSE])
  }, numeric(1))

  expect_equal(batch, adapter$predict(second, new_x), tolerance = 0)
  expect_equal(batch, one_at_a_time, tolerance = 0)
  expect_error(
    adapter$predict(first, new_x[c("gene_b", "gene_a", "gene_c")]),
    class = "stormbreaker_learner_schema_mismatch"
  )

  local_mocked_bindings(sb_ranger_available = function() FALSE)
  expect_error(
    adapter$fit(x, d$outcome),
    class = "stormbreaker_missing_ranger"
  )
  local_mocked_bindings(sb_ranger_available = function() TRUE)
  expect_error(
    adapter$fit(x, rep(0L, nrow(x))),
    class = "stormbreaker_ranger_binary_single_class"
  )
  too_wide <- sb_ranger_binary_learner_adapter(
    num_trees = 5L,
    mtry = 4L,
    seed = 7713L
  )
  expect_error(
    too_wide$fit(x, d$outcome),
    class = "stormbreaker_ranger_mtry_exceeds_features"
  )
})

test_that("ranger binary runs all-feature and fold-local Top-K workflows", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  d <- make_sb_binary_fixture(n = 72L, seed = 7721L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  all_features <- sb_pipeline(
    learner = "ranger_binary",
    ranger_num_trees = 12L,
    ranger_seed = 7722L
  )
  screened <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "ranger_binary",
    selector_n_features = 2L,
    ranger_num_trees = 12L,
    ranger_mtry = 1L,
    ranger_seed = 7723L
  )
  result <- sb_run_cv(task, all_features, folds = 3L, seed = 7724L)
  screened_result <- sb_run_cv(task, screened, folds = 3L, seed = 7724L)

  expect_identical(all_features$task_type, "binary")
  expect_identical(
    all_features$id,
    "none__ranger_binary__trees12__mtrysqrt__minnode10__seed7722"
  )
  expect_identical(
    screened$id,
    paste0(
      "univ_binary_auc_top_k__ranger_binary__k2",
      "__trees12__mtry1__minnode10__seed7723"
    )
  )
  expect_identical(all_features$learner_params$min_node_size, 10L)
  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_true(all(is.finite(result$predictions$probability)))
  expect_length(result$tree_models, 3L)
  expect_identical(
    vapply(result$tree_models, `[[`, "", "status"),
    rep("success", 3L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_applicable", 3L)
  )
  expect_identical(
    result$provenance$ranger_version,
    as.character(utils::packageVersion("ranger"))
  )
  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    audit <- result$tree_models[[fold]]
    rows <- match(split$assessment_ids, d$patient_id)
    rebuilt <- stats::predict(
      unserialize(audit$model_raw),
      data = d[rows, audit$features, drop = FALSE],
      num.trees = audit$num_trees,
      type = "response",
      seed = 0L,
      num.threads = 1L,
      verbose = FALSE
    )$predictions[, "1"]
    observed <- result$predictions$probability[
      match(split$assessment_ids, result$predictions$patient_id)
    ]
    expect_equal(observed, unname(rebuilt), tolerance = 0)
  }
  expect_identical(screened_result$status, "success")
  expect_true(all(vapply(
    screened_result$selections,
    function(selection) length(selection$selected_features) == 2L,
    logical(1)
  )))

  fit <- sb_fit_final(task, screened)
  rows <- c(5L, 21L, 49L)
  new_data <- d[rows, c("patient_id", rev(predictors)), drop = FALSE]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  direct <- sb_predict_ranger_binary_learner_state(
    fit$learner_state,
    new_data[fit$required_features]
  )

  expect_identical(fit$status, "success")
  expect_identical(fit$learner_state$learner, "ranger_binary")
  expect_type(fit$learner_state$model_raw, "raw")
  expect_identical(prediction$status, "success")
  expect_equal(prediction$predictions$probability, direct$probability, tolerance = 0)
  expect_identical(
    fit$provenance$ranger_version,
    as.character(utils::packageVersion("ranger"))
  )
})

test_that("ranger binary informational warnings retain native results", {
  skip_if_not_installed("ranger")
  x <- data.frame(a = sin(1:30), b = cos(1:30))
  adapter <- sb_ranger_binary_learner_adapter(num_trees = 5L)
  reference <- adapter$fit(x, rep(0:1, 15L))
  native <- ranger::ranger
  local_mocked_bindings(ranger = function(...) {
    warning("informational native warning")
    native(...)
  }, .package = "ranger")
  seen <- character()
  observed <- withCallingHandlers(adapter$fit(x, rep(0:1, 15L)), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(seen, "informational native warning")
  expect_identical(observed$forest, reference$forest)
  expect_identical(adapter$predict(observed, x), adapter$predict(reference, x))
})

test_that("ranger binary still rejects ignored model arguments", {
  skip_if_not_installed("ranger")
  x <- data.frame(a = sin(1:30), b = cos(1:30))
  adapter <- sb_ranger_binary_learner_adapter(num_trees = 5L)
  message <- "Unused arguments: misspelled_parameter"
  local_mocked_bindings(ranger = function(...) {
    warning(message)
    stop("ignored argument warning should stop fitting first")
  }, .package = "ranger")
  expect_s3_class(tryCatch(adapter$fit(x, rep(0:1, 15L)), error = identity), "stormbreaker_ranger_binary_warning")
  message <- "Argument 'num.random.splits' ignored if splitrule is not 'extratrees'."
  expect_s3_class(tryCatch(adapter$fit(x, rep(0:1, 15L)), error = identity), "stormbreaker_ranger_binary_warning")
})

test_that("ranger binary still rejects incomplete forests", {
  skip_if_not_installed("ranger")
  x <- data.frame(a = sin(1:30), b = cos(1:30))
  adapter <- sb_ranger_binary_learner_adapter(num_trees = 5L)
  native <- ranger::ranger
  local_mocked_bindings(ranger = function(...) {
    fitted <- native(...)
    fitted$forest <- NULL
    fitted
  }, .package = "ranger")
  expect_s3_class(tryCatch(adapter$fit(x, rep(0:1, 15L)), error = identity), "stormbreaker_incomplete_ranger_fit")
})

# ---- ranger-splitrule-variants ----

make_ranger_splitrule_binary_fixture <- function(n = 84L, seed = 8841L) {
  withr::with_seed(seed, {
    outcome <- sample(rep(c(0L, 1L), length.out = n))
    data.frame(
      patient_id = sprintf("RSB%03d", seq_len(n)),
      outcome = outcome,
      "gene one" = stats::rnorm(n, 1.1 * outcome),
      "2gene" = stats::rnorm(n, -0.8 * outcome),
      "gene-three" = stats::rnorm(n),
      check.names = FALSE
    )
  })
}

make_ranger_splitrule_survival_fixture <- function(n = 84L, seed = 8842L) {
  withr::with_seed(seed, {
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    x3 <- stats::rnorm(n)
    event_time <- stats::rexp(
      n,
      rate = 0.08 * exp(x1 - 0.7 * x2)
    ) + 0.1
    censor_time <- stats::rexp(n, rate = 0.05) + 0.1
    data.frame(
      patient_id = sprintf("RSS%03d", seq_len(n)),
      time = pmin(event_time, censor_time),
      event = as.integer(event_time <= censor_time),
      "gene one" = x1,
      "2gene" = x2,
      "gene-three" = x3,
      check.names = FALSE
    )
  })
}

test_that("ranger split-rule learners declare five distinct specifications", {
  binary <- c(
    "ranger_binary",
    "ranger_extratrees_binary",
    "ranger_hellinger_binary"
  )
  survival <- c(
    "ranger_survival",
    "ranger_extratrees_survival",
    "ranger_cindex_survival",
    "ranger_maxstat_survival"
  )

  expect_identical(sb_ranger_binary_learner_ids(), binary)
  expect_identical(sb_ranger_survival_learner_ids(), survival)
  expect_identical(sb_ranger_learner_ids(), c(binary, survival))

  binary_splitrules <- c("gini", "extratrees", "hellinger")
  for (i in seq_along(binary)) {
    params <- sb_ranger_binary_params(learner = binary[[i]])
    expect_identical(params$splitrule, binary_splitrules[[i]])
    expect_identical(params$replace, !identical(binary_splitrules[[i]], "extratrees"))
    expect_identical(params$probability, TRUE)
  }
  expect_identical(
    sb_ranger_binary_params(
      learner = "ranger_extratrees_binary"
    )$num_random_splits,
    1L
  )

  survival_splitrules <- c("logrank", "extratrees", "C", "maxstat")
  for (i in seq_along(survival)) {
    params <- sb_ranger_survival_params(learner = survival[[i]])
    expect_identical(params$splitrule, survival_splitrules[[i]])
    expect_identical(
      params$replace,
      !identical(survival_splitrules[[i]], "extratrees")
    )
  }
  extra <- sb_ranger_survival_params(
    learner = "ranger_extratrees_survival"
  )
  maxstat <- sb_ranger_survival_params(
    learner = "ranger_maxstat_survival"
  )
  expect_identical(extra$num_random_splits, 1L)
  expect_identical(maxstat$alpha, 0.5)
  expect_identical(maxstat$minprop, 0.1)

  binary_selectors <- expected_fixed_learner_selectors("binary")
  survival_selectors <- expected_fixed_learner_selectors("survival")
  registry <- sb_pipeline_registry()
  rows <- registry[
    registry$learner %in% c(binary[-1L], survival[-1L]), ,
    drop = FALSE
  ]
  expect_identical(
    nrow(rows),
    length(binary[-1L]) * length(binary_selectors) +
      length(survival[-1L]) * length(survival_selectors)
  )
  expect_identical(
    sort(unique(rows$selector[rows$task_type == "binary"])),
    sort(binary_selectors)
  )
  expect_identical(
    sort(unique(rows$selector[rows$task_type == "survival"])),
    sort(survival_selectors)
  )
  expect_route_dependencies(rows, stats::setNames(
    rep("ranger", length(c(binary[-1L], survival[-1L]))), c(binary[-1L], survival[-1L])))
})

test_that("new ranger split rules support one non-syntactic feature", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  binary_data <- make_ranger_splitrule_binary_fixture(n = 48L, seed = 8849L)
  binary_x <- binary_data["gene one"]
  for (learner in sb_ranger_binary_learner_ids()[-1L]) {
    adapter <- sb_ranger_binary_learner_adapter(
      num_trees = 12L,
      learner = learner
    )
    model <- adapter$fit(binary_x, binary_data$outcome)
    prediction <- adapter$predict(model, binary_x[1:4, , drop = FALSE])

    expect_length(prediction, 4L)
    expect_true(all(is.finite(prediction)))
  }

  survival_data <- make_ranger_splitrule_survival_fixture(
    n = 48L,
    seed = 8850L
  )
  survival_x <- survival_data["gene one"]
  for (learner in sb_ranger_survival_learner_ids()[-1L]) {
    adapter <- sb_ranger_survival_learner_adapter(
      num_trees = 12L,
      learner = learner
    )
    model <- adapter$fit(
      survival_x,
      survival_data$time,
      survival_data$event
    )
    prediction <- adapter$predict(
      model,
      survival_x[1:4, , drop = FALSE]
    )

    expect_length(prediction, 4L)
    expect_true(all(is.finite(prediction)))
  }
})

test_that("binary ranger split rules replay native probabilities", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  data <- make_ranger_splitrule_binary_fixture()
  x <- data[c("gene one", "2gene", "gene-three")]
  new_x <- x[c(3L, 27L, 68L), , drop = FALSE]
  withr::local_seed(8843L)

  for (learner in sb_ranger_binary_learner_ids()[-1L]) {
    adapter <- sb_ranger_binary_learner_adapter(
      num_trees = 24L,
      mtry = 1L,
      min_node_size = 10L,
      seed = 1L,
      learner = learner
    )
    rng_before <- .Random.seed
    model <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, rng_before)
    repeated <- adapter$fit(x, data$outcome)
    expect_identical(
      serialize(model$forest, NULL),
      serialize(repeated$forest, NULL)
    )

    native <- stats::predict(
      model$forest,
      data = new_x,
      num.trees = 24L,
      type = "response",
      seed = 0L,
      num.threads = 1L,
      verbose = FALSE
    )$predictions[, "1"]
    observed <- adapter$predict(model, new_x)
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, new_x)

    expect_identical(model$learner, learner)
    expect_equal(observed, unname(native), tolerance = 0)
    expect_equal(replay$probability, observed, tolerance = 0)
    expect_identical(state$learner, learner)
    expect_identical(state$params$splitrule, adapter$params$splitrule)
    expect_type(state$model_raw, "raw")
    expect_false("forest" %in% names(state))
  }
})

test_that("survival ranger split rules replay native mortality and survival", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  data <- make_ranger_splitrule_survival_fixture()
  x <- data[c("gene one", "2gene", "gene-three")]
  new_x <- x[c(3L, 27L, 68L), , drop = FALSE]
  horizon <- as.numeric(stats::quantile(data$time, 0.4))
  withr::local_seed(8844L)

  for (learner in sb_ranger_survival_learner_ids()[-1L]) {
    adapter <- sb_ranger_survival_learner_adapter(
      num_trees = 24L,
      mtry = 1L,
      min_node_size = 3L,
      seed = 1L,
      learner = learner
    )
    rng_before <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, rng_before)
    repeated <- adapter$fit(x, data$time, data$event)
    expect_identical(
      serialize(model$forest, NULL),
      serialize(repeated$forest, NULL)
    )

    native <- stats::predict(
      model$forest,
      data = new_x,
      num.trees = 24L,
      type = "response",
      seed = 0L,
      num.threads = 1L,
      verbose = FALSE
    )
    risk <- adapter$predict(model, new_x)
    probability <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      new_x,
      horizon
    )
    eligible <- which(native$unique.death.times <= horizon)
    native_probability <- if (length(eligible) == 0L) {
      rep(1, nrow(new_x))
    } else {
      unname(native$survival[, max(eligible)])
    }
    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = horizon,
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, new_x)

    expect_identical(model$learner, learner)
    expect_equal(risk, unname(rowSums(native$chf)), tolerance = 0)
    expect_equal(probability$probability, native_probability, tolerance = 0)
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      probability$probability,
      tolerance = 0
    )
    expect_identical(state$learner, learner)
    expect_identical(state$params$splitrule, adapter$params$splitrule)
    expect_null(state$baseline_hazard)
    expect_type(state$model_raw, "raw")
  }
})

test_that("new ranger split rules execute through CV and final prediction", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  binary_data <- make_ranger_splitrule_binary_fixture(n = 72L, seed = 8845L)
  binary_features <- c("gene one", "2gene", "gene-three")
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    binary_features
  )
  for (learner in sb_ranger_binary_learner_ids()[-1L]) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 2L,
      ranger_num_trees = 18L
    )
    comparison <- sb_run_cv(binary_task, pipeline, folds = 3L, seed = 8846L)
    fit <- sb_fit_final(binary_task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      binary_data[1:5, c("patient_id", rev(binary_features))],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_identical(fit$learner_state$learner, learner)
  }

  survival_data <- make_ranger_splitrule_survival_fixture(
    n = 72L,
    seed = 8847L
  )
  survival_features <- c("gene one", "2gene", "gene-three")
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    survival_features,
    "month",
    prediction_horizon = as.numeric(stats::quantile(survival_data$time, 0.4))
  )
  for (learner in sb_ranger_survival_learner_ids()[-1L]) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 2L,
      ranger_num_trees = 18L
    )
    comparison <- sb_run_cv(
      survival_task,
      pipeline,
      folds = 3L,
      seed = 8848L
    )
    fit <- sb_fit_final(survival_task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      survival_data[1:5, c("patient_id", rev(survival_features))],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(
      comparison$complete_survival_probability_oof,
      TRUE
    )
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_identical(fit$learner_state$learner, learner)
  }
})

# ---- ranger-survival-fixed-adapter ----

test_that("C-index execution threads preserve forest and portable predictions", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")
  d <- make_sb_fixture(n = 54L, seed = 7601L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_ranger_survival_learner_adapter(
    num_trees = 18L, learner = "ranger_cindex_survival")
  one <- withr::with_options(list(stormbreaker.ranger_cindex_threads = 1L),
    adapter$fit(x, d$time, d$event))
  four <- withr::with_options(list(stormbreaker.ranger_cindex_threads = 4L),
    adapter$fit(x, d$time, d$event))
  expect_identical(one$params$num_threads, 1L)
  expect_identical(four$params$num_threads, 4L)
  expect_identical(one$forest, four$forest)
  expect_identical(adapter$predict(one, x), adapter$predict(four, x))
  context <- list(x = x, time = d$time, event = d$event,
    prediction_horizon = 7, time_unit = "month")
  s1 <- adapter$extract_state(one, context)
  s4 <- unserialize(serialize(adapter$extract_state(four, context), NULL))
  expect_identical(adapter$predict_state(s1, x), adapter$predict_state(s4, x))
  expect_equal(one$feature_importance, four$feature_importance, tolerance = 1e-12)
})

test_that("fixed ranger survival state reconstructs risk and direct survival", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  d <- make_sb_fixture(n = 54L, seed = 7601L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(4L, 19L, 37L), , drop = FALSE]
  adapter <- sb_ranger_survival_learner_adapter(
    num_trees = 18L,
    mtry = NULL,
    min_node_size = 3L,
    seed = 7602L
  )
  model <- adapter$fit(x, d$time, d$event)
  observed_risk <- adapter$predict(model, new_x)
  observed_survival <- adapter$predict_survival(
    model,
    x,
    d$time,
    d$event,
    new_x,
    prediction_horizon = 7
  )
  audit <- adapter$extract_tree_model(model)

  restored <- unserialize(audit$model_raw)
  oracle <- stats::predict(
    restored,
    data = new_x,
    num.trees = 18L,
    type = "response",
    seed = 0L,
    num.threads = 1L,
    verbose = FALSE
  )
  expected_risk <- unname(rowSums(oracle$chf))
  eligible <- which(oracle$unique.death.times <= 7)
  expected_survival <- if (length(eligible) == 0L) {
    rep(1, nrow(new_x))
  } else {
    unname(oracle$survival[, max(eligible)])
  }
  state <- adapter$extract_state(model, list(
    x = x,
    time = d$time,
    event = d$event,
    prediction_horizon = 7,
    time_unit = "month"
  ))
  state_prediction <- adapter$predict_state(state, new_x)

  expect_s3_class(model, "stormbreaker_ranger_survival")
  expect_identical(model$effective_mtry, 1L)
  expect_equal(observed_risk, expected_risk, tolerance = 0)
  expect_equal(observed_survival$probability, expected_survival, tolerance = 0)
  expect_identical(observed_survival$prediction_horizon, 7)
  expect_identical(
    observed_survival$source,
    "ranger_survival_forest_direct"
  )
  expect_identical(observed_survival$ties, "not_applicable")
  expect_false("cumulative_baseline_hazard" %in% names(observed_survival))
  expect_identical(audit$num_trees, 18L)
  expect_identical(audit$effective_mtry, 1L)
  expect_identical(audit$features, names(x))
  expect_identical(audit$params, adapter$params)
  expect_identical(audit$importance_type, "impurity")
  expect_identical(audit$model_serialization, "R_xdr_v3_ranger.forest")
  expect_type(audit$model_raw, "raw")
  expect_gt(length(audit$model_raw), 0L)
  expect_identical(audit$feature_importance$feature, names(x))
  expect_true(all(is.finite(audit$feature_importance$importance)))
  expect_identical(state$learner, "ranger_survival")
  expect_false("forest" %in% names(state))
  expect_identical(
    state$survival_probability$source,
    "ranger_survival_forest_direct"
  )
  expect_equal(state_prediction$risk, expected_risk, tolerance = 0)
  expect_equal(
    state_prediction$survival_probability,
    expected_survival,
    tolerance = 0
  )
  expect_identical(state_prediction$prediction_horizon, 7)
  expect_identical(state_prediction$time_unit, "month")

  risk_only_state <- adapter$extract_state(model, list())
  expect_null(risk_only_state$survival_probability)
  expect_identical(
    names(adapter$predict_state(risk_only_state, new_x)),
    "risk"
  )
})

test_that("fixed ranger survival is deterministic and batching independent", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  d <- make_sb_fixture(n = 50L, seed = 7611L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_ranger_survival_learner_adapter(
    num_trees = 15L,
    mtry = 2L,
    min_node_size = 3L,
    seed = 7612L
  )
  first <- adapter$fit(x[1:40, , drop = FALSE], d$time[1:40], d$event[1:40])
  second <- adapter$fit(x[1:40, , drop = FALSE], d$time[1:40], d$event[1:40])
  new_x <- x[41:50, , drop = FALSE]

  batch <- adapter$predict(first, new_x)
  repeated <- adapter$predict(second, new_x)
  one_at_a_time <- vapply(seq_len(nrow(new_x)), function(row) {
    adapter$predict(first, new_x[row, , drop = FALSE])
  }, numeric(1))
  first_survival <- adapter$predict_survival(
    first,
    x[1:40, , drop = FALSE],
    d$time[1:40],
    d$event[1:40],
    new_x,
    7
  )$probability
  second_survival <- adapter$predict_survival(
    second,
    x[1:40, , drop = FALSE],
    d$time[1:40],
    d$event[1:40],
    new_x,
    7
  )$probability

  expect_equal(batch, repeated, tolerance = 0)
  expect_equal(batch, one_at_a_time, tolerance = 0)
  expect_equal(first_survival, second_survival, tolerance = 0)
  expect_error(
    adapter$predict(first, new_x[c("gene_b", "gene_a", "gene_c")]),
    class = "stormbreaker_learner_schema_mismatch"
  )
})

test_that("fixed ranger survival failures remain explicit", {
  d <- make_sb_fixture(n = 24L, seed = 7621L)
  x <- d[c("gene_a", "gene_b")]
  adapter <- sb_ranger_survival_learner_adapter(
    num_trees = 5L,
    min_node_size = 3L,
    seed = 7622L
  )

  local_mocked_bindings(sb_ranger_available = function() FALSE)
  expect_error(
    adapter$fit(x, d$time, d$event),
    class = "stormbreaker_missing_ranger"
  )

  local_mocked_bindings(sb_ranger_available = function() TRUE)
  expect_error(
    adapter$fit(x, d$time, rep(0L, nrow(d))),
    class = "stormbreaker_ranger_survival_no_event"
  )
  too_wide <- sb_ranger_survival_learner_adapter(
    num_trees = 5L,
    mtry = 3L,
    min_node_size = 3L,
    seed = 7623L
  )
  expect_error(
    too_wide$fit(x, d$time, d$event),
    class = "stormbreaker_ranger_mtry_exceeds_features"
  )
})

test_that("ranger survival runs audited CV and final direct prediction", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  d <- make_sb_fixture(n = 72L, seed = 7631L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  horizon_task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )
  risk_task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month"
  )
  pipeline <- sb_pipeline(
    learner = "ranger_survival",
    ranger_num_trees = 10L,
    ranger_min_node_size = 3L,
    ranger_seed = 7632L
  )
  screened <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "ranger_survival",
    selector_n_features = 2L,
    ranger_num_trees = 10L,
    ranger_mtry = 1L,
    ranger_min_node_size = 3L,
    ranger_seed = 7633L
  )
  result <- sb_run_cv(horizon_task, pipeline, folds = 3L, seed = 7634L)
  screened_result <- sb_run_cv(
    horizon_task,
    screened,
    folds = 3L,
    seed = 7634L
  )
  risk_result <- sb_run_cv(risk_task, pipeline, folds = 3L, seed = 7634L)

  expect_identical(pipeline$task_type, "survival")
  expect_identical(pipeline$risk_type, "survival_forest_mortality")
  expect_identical(pipeline$risk_reference, "ranger_native")
  expect_identical(pipeline$risk_formula, "row_sums_cumulative_hazard")
  expect_identical(pipeline$learner_params$mtry_rule, "floor_sqrt_p")
  expect_identical(result$status, "success")
  expect_true(result$complete_oof)
  expect_true(result$complete_survival_probability_oof)
  expect_true(all(is.finite(result$predictions$risk)))
  expect_true(all(is.finite(result$predictions$survival_probability)))
  expect_identical(
    unique(result$predictions$learner),
    "ranger_survival"
  )
  expect_identical(result$metrics$harrell_c$status, "success")
  expect_identical(result$metrics$ipcw_brier$status, "success")
  expect_length(result$tree_models, 3L)
  expect_true(all(vapply(
    result$tree_models,
    function(audit) identical(audit$status, "success") &&
      is.raw(audit$model_raw) && length(audit$model_raw) > 0L,
    logical(1)
  )))
  expect_true(all(vapply(
    result$baseline_hazards,
    function(audit) identical(audit$status, "not_applicable") &&
      identical(audit$source, "ranger_survival_forest_direct") &&
      is.na(audit$cumulative_baseline_hazard),
    logical(1)
  )))
  expect_identical(
    result$provenance$ranger_version,
    as.character(utils::packageVersion("ranger"))
  )
  expect_identical(screened_result$status, "success")
  expect_length(screened_result$selections, 3L)
  expect_true(all(vapply(
    screened_result$selections,
    function(selection) length(selection$selected_features) == 2L,
    logical(1)
  )))
  expect_identical(risk_result$status, "success")
  expect_true(risk_result$complete_oof)
  expect_false("survival_probability" %in% names(risk_result$predictions))
  expect_false("baseline_hazards" %in% names(risk_result))

  fit <- sb_fit_final(horizon_task, pipeline)
  rows <- c(5L, 21L, 49L)
  new_data <- d[rows, c("patient_id", rev(predictors)), drop = FALSE]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  selected_x <- new_data[fit$required_features]
  direct <- sb_predict_ranger_survival_learner_state(
    fit$learner_state,
    selected_x
  )

  expect_identical(fit$status, "success")
  expect_s3_class(fit$learner_state, "stormbreaker_learner_state")
  expect_identical(fit$learner_state$learner, "ranger_survival")
  expect_type(fit$learner_state$model_raw, "raw")
  expect_null(fit$learner_state$baseline_hazard)
  expect_identical(prediction$status, "success")
  expect_equal(prediction$predictions$risk, direct$risk, tolerance = 0)
  expect_equal(
    prediction$predictions$survival_probability,
    direct$survival_probability,
    tolerance = 0
  )
  expect_identical(
    fit$provenance$ranger_version,
    as.character(utils::packageVersion("ranger"))
  )
})

test_that("ranger survival informational warnings retain native results", {
  skip_if_not_installed("ranger")
  x <- data.frame(a = sin(1:30), b = cos(1:30))
  adapter <- sb_ranger_survival_learner_adapter(num_trees = 5L)
  reference <- adapter$fit(x, seq_len(30), rep(0:1, 15L))
  native <- ranger::ranger
  local_mocked_bindings(ranger = function(...) {
    warning("informational native warning")
    native(...)
  }, .package = "ranger")
  seen <- character()
  observed <- withCallingHandlers(adapter$fit(x, seq_len(30), rep(0:1, 15L)), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(seen, "informational native warning")
  expect_identical(observed$forest, reference$forest)
  expect_identical(adapter$predict(observed, x), adapter$predict(reference, x))
})

test_that("ranger survival still rejects ignored model arguments", {
  skip_if_not_installed("ranger")
  x <- data.frame(a = sin(1:30), b = cos(1:30))
  adapter <- sb_ranger_survival_learner_adapter(num_trees = 5L)
  message <- "Unused arguments: misspelled_parameter"
  local_mocked_bindings(ranger = function(...) {
    warning(message)
    stop("ignored argument warning should stop fitting first")
  }, .package = "ranger")
  expect_s3_class(tryCatch(adapter$fit(x, seq_len(30), rep(0:1, 15L)), error = identity), "stormbreaker_ranger_survival_warning")
  message <- "Argument 'num.random.splits' ignored if splitrule is not 'extratrees'."
  expect_s3_class(tryCatch(adapter$fit(x, seq_len(30), rep(0:1, 15L)), error = identity), "stormbreaker_ranger_survival_warning")
})

test_that("ranger survival still rejects incomplete forests", {
  skip_if_not_installed("ranger")
  x <- data.frame(a = sin(1:30), b = cos(1:30))
  adapter <- sb_ranger_survival_learner_adapter(num_trees = 5L)
  native <- ranger::ranger
  local_mocked_bindings(ranger = function(...) {
    fitted <- native(...)
    fitted$forest <- NULL
    fitted
  }, .package = "ranger")
  expect_s3_class(tryCatch(adapter$fit(x, seq_len(30), rep(0:1, 15L)), error = identity), "stormbreaker_incomplete_ranger_fit")
})
})

# ---- model-rlt ----
local({
test_that("RLT views retain all native trees and axis-aligned split indices", {
  skip_if_not_installed("RLT")
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in sb_rlt_learner_ids()) {
    binary <- id %in% sb_rlt_binary_learner_ids()
    d <- if (binary) sb_example_binary(100L, seed = 108L) else sb_example_survival(100L, seed = 8612L)
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    prediction <- sb_predict_final(fit, d, "patient_id")
    e <- unserialize(fit$learner_state$engine_raw)
    sizes <- sb_rlt_diagnostics(fit, "tree_size")
    usage <- sb_rlt_diagnostics(fit)
    leaves <- vapply(e$FittedTrees, function(z) sum(!is.finite(z["NextLeft", ]) & !is.finite(z["NextRight", ])), integer(1))
    expect_equal(sizes$table$terminal_nodes, leaves)
    expect_equal(sizes$table$nodes, vapply(e$FittedTrees, ncol, integer(1)))
    expect_equal(sizes$table$nodes, 2L * sizes$table$terminal_nodes - 1L)
    vars <- unlist(lapply(e$FittedTrees, function(z) z["SplitVar1", is.finite(z["NextLeft", ])]), use.names = FALSE)
    expect_equal(usage$table$splits, as.numeric(table(factor(vars, levels = seq_along(features)))))
    expect_equal(usage$table$feature, features)
    expect_equal(sum(usage$table$splits), sum(sizes$table$nodes - sizes$table$terminal_nodes))
    expect_equal(nrow(sizes$table), 100L)
    expect_null(fit$diagnostic_artifacts)
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 2L))
    for (v in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(v)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_null(p$labels$caption)
      expect_null(p$labels$subtitle)
      if (v$type == "rlt_tree_size") expect_equal(sum(p$data$trees), 100L) else
        expect_equal(p$data$splits, usage$table$splits)
    }
    expect_identical(sb_predict_final(fit, d, "patient_id"), prediction)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
})

# ---- model-wsrf ----
local({
test_that("WSRF diagnostics never confuse compact placeholders with native metrics", {
  skip_if_not_installed("wsrf")
  d <- sb_example_binary(160L, seed = 108L)
  for (features in list(c("gene_c", "gene_a", "gene_b"), "gene_a")) {
    task <- sb_binary_task(d, "patient_id", "outcome", features)
    fit <- sb_fit_final(task, sb_pipeline(learner = "wsrf_binary"), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    prediction <- sb_predict_final(fit, d, "patient_id")
    x <- d[features]
    names(x) <- paste0(".stormbreaker_wsrf_x", seq_along(features))
    if (length(features) == 1L) x$.stormbreaker_wsrf_zero <- 0
    native <- withr::with_seed(1L, wsrf::wsrf(x,
      factor(d$outcome, levels = c(0, 1), labels = c("control", "case")),
      mtry = as.integer(floor(log2(length(features)) + 1)), ntree = 500L,
      weights = TRUE, parallel = FALSE, importance = FALSE, nodesize = 2L))
    importance <- sb_wsrf_diagnostics(fit)
    nodes <- sb_wsrf_diagnostics(fit, "tree_nodes")
    oob <- sb_wsrf_diagnostics(fit, "tree_oob")
    expect_equal(importance$table$importance, unname(native$importance[, "MeanDecreaseIGR"]), tolerance = 0)
    expect_equal(importance$table$feature[!importance$table$artificial_bridge], features)
    expect_equal(sum(importance$table$artificial_bridge), as.integer(length(features) == 1L))
    expect_equal(nodes$table$nodes, lengths(native$trees))
    expect_equal(oob$table$error, native$treeOOBErrorRates, tolerance = 0)
    expect_equal(oob$forest_oob, native$RFOOBErrorRate, tolerance = 0)
    expect_gt(sum(oob$table$error > 0), 0)
    expect_null(fit$learner_state$diagnostic_metrics)
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_equal(unserialize(compact$learner_state$engine_raw)$treeOOBErrorRates, rep(0, 500))
    expect_identical(sb_wsrf_diagnostics(compact, "tree_nodes"), nodes)
    expect_identical(sb_native_diagnostics(compact, complete = TRUE)$status$status,
      c("missing_capture", "available", "missing_capture"))
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 3))
    for (v in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(v)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_null(p$labels$caption)
      expect_null(p$labels$subtitle)
      if (v$type == "wsrf_tree_nodes") expect_equal(sum(p$data$trees), 500L)
      if (v$type == "wsrf_tree_oob") expect_equal(p$data$error, native$treeOOBErrorRates, tolerance = 0)
    }
    expect_identical(sb_predict_final(fit, d, "patient_id"), prediction)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
})
