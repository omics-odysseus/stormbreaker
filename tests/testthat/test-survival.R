# ---- additional-survival-fixed-adapters ----
local({
test_that("additional survival specifications are explicit", {
  survreg <- sb_additional_survival_params("survreg_weibull")
  rpart <- sb_additional_survival_params("rpart_survival")
  rsf <- sb_additional_survival_params("rsf_survival")
  ctree <- sb_additional_survival_params("ctree_survival")
  cforest <- sb_additional_survival_params("cforest_survival")
  bart <- sb_additional_survival_params("bart_survival")

  expect_identical(survreg, list(
    distribution = "weibull",
    scale = "estimated"
  ))
  expect_identical(rpart$method, "exp")
  expect_identical(rpart[c("cp", "minsplit", "maxdepth", "xval")], list(
    cp = 0.01,
    minsplit = 20L,
    maxdepth = 30L,
    xval = 0L
  ))
  expect_identical(rsf$num_trees, 1000L)
  expect_identical(rsf$mtry_rule, "randomForestSRC_default")
  expect_identical(rsf$min_node_size, 5L)
  expect_identical(rsf$splitrule, "logrank")
  expect_identical(rsf$seed, 1L)
  expect_identical(ctree$mincriterion, 0.95)
  expect_identical(ctree$minsplit, 20L)
  expect_identical(ctree$minbucket, 7L)
  expect_identical(cforest$num_trees, 500L)
  expect_identical(cforest$mtry_rule, "ceiling_sqrt_p")
  expect_identical(cforest$replace, TRUE)
  expect_identical(cforest$sample_fraction, 1)
  expect_identical(cforest$seed, 1L)
  expect_identical(bart, list(
    type = "pbart",
    num_trees = 50L,
    num_cutpoints = 100L,
    time_grid_size = 50L,
    posterior_draws = 250L,
    burn_in = 100L,
    keep_every = 1L,
    prior_k = 2,
    tree_power = 2,
    tree_base = 0.95,
    sparse = FALSE,
    seed = 1L,
    cores = 1L
  ))
})

test_that("Weibull SurvReg nonconvergence is a typed learner failure", {
  value <- seq(-1, 1, length.out = 18L)
  data <- data.frame(
    patient_id = sprintf("patient_%02d", seq_along(value)),
    time = exp(value),
    event = as.integer(!seq_along(value) %in% c(6L, 12L, 18L)),
    gene_a = value
  )
  adapter <- sb_additional_survival_learner_adapter(
    "survreg_weibull"
  )

  error <- capture_sb_error(
    adapter$fit(data["gene_a"], data$time, data$event)
  )
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    "gene_a",
    "month",
    prediction_horizon = 1
  )
  pipeline <- sb_pipeline(learner = "survreg_weibull")
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 7920L)
  fit <- sb_fit_final(task, pipeline)

  expect_s3_class(error, "stormbreaker_survreg_nonconvergence")
  expect_match(
    conditionMessage(error),
    "Ran out of iterations and did not converge",
    fixed = TRUE
  )
  expect_s3_class(error$parent, "warning")
  expect_identical(result$status, "failed")
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_survreg_nonconvergence"
  )
  expect_identical(fit$status, "failed")
  expect_identical(fit$error_stage, "learner_fit")
  expect_identical(fit$error_class, "stormbreaker_survreg_nonconvergence")
})

test_that("survival BART uses native posterior curves and declared risk", {
  skip_if_not_installed("BART", minimum_version = "2.9.0")

  d <- make_sb_fixture(n = 54L, seed = 7922L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(4L, 21L, 47L), , drop = FALSE]
  params <- sb_additional_survival_params("bart_survival")
  params$time_grid_size <- 20L
  params$posterior_draws <- 40L
  params$burn_in <- 20L
  adapter <- sb_additional_survival_learner_adapter(
    "bart_survival",
    params
  )
  model <- adapter$fit(x, d$time, d$event)
  curves <- sb_bart_survival_curves(model, new_x)

  grid <- model$engine$times
  augmented <- cbind(
    t = rep(grid, times = nrow(new_x)),
    as.matrix(new_x)[
      rep(seq_len(nrow(new_x)), each = length(grid)),
      ,
      drop = FALSE
    ]
  )
  if (length(model$engine$rm.const) > 0L) {
    augmented <- augmented[, model$engine$rm.const, drop = FALSE]
  }
  native <- NULL
  method <- get(
    "predict.survbart",
    envir = asNamespace("BART"),
    inherits = FALSE
  )
  invisible(utils::capture.output(
    native <- method(
      model$engine,
      augmented,
      mc.cores = 1L,
      openmp = FALSE
    ),
    type = "output"
  ))
  native_survival <- matrix(
    native$surv.test.mean,
    nrow = nrow(new_x),
    ncol = length(grid),
    byrow = TRUE
  )
  interval <- diff(c(0, grid))
  native_risk <- -as.vector(cbind(
    1,
    native_survival[, -ncol(native_survival), drop = FALSE]
  ) %*% interval)
  horizon <- stats::median(grid)
  horizon_column <- max(which(grid <= horizon))
  probability <- adapter$predict_survival(
    model,
    x,
    d$time,
    d$event,
    new_x,
    horizon
  )

  expect_s3_class(model$engine, "survbart")
  expect_equal(curves$survival, native_survival, tolerance = 1e-14)
  expect_equal(adapter$predict(model, new_x), native_risk, tolerance = 1e-14)
  expect_equal(
    probability$probability,
    native_survival[, horizon_column],
    tolerance = 1e-14
  )
  expect_identical(
    probability$source,
    "BART_discrete_time_posterior_survival_direct"
  )
})

test_that("additional survival adapters return risk, horizon survival, and state", {
  skip_if_not_installed("rpart", minimum_version = "4.1")
  skip_if_not_installed("randomForestSRC", minimum_version = "3.3.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")
  skip_if_not_installed("coin", minimum_version = "1.1")
  skip_if_not_installed("BART", minimum_version = "2.9.0")

  d <- make_sb_fixture(n = 60L, seed = 7921L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(4L, 23L, 49L), , drop = FALSE]
  learners <- sb_additional_survival_learner_ids()

  for (learner in learners) {
    params <- sb_additional_survival_params(learner)
    if (identical(learner, "rsf_survival")) {
      params$num_trees <- 20L
    }
    if (identical(learner, "cforest_survival")) {
      params$num_trees <- 20L
    }
    if (identical(learner, "bart_survival")) {
      params$time_grid_size <- 20L
      params$posterior_draws <- 40L
      params$burn_in <- 20L
    }
    if (identical(learner, "aorsf_survival")) {
      params$num_trees <- 20L
    }
    adapter <- sb_additional_survival_learner_adapter(learner, params)
    model <- adapter$fit(x, d$time, d$event)
    risk <- adapter$predict(model, new_x)
    probability <- adapter$predict_survival(
      model,
      x,
      d$time,
      d$event,
      new_x,
      prediction_horizon = 7
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
      new_x
    )

    expect_identical(state$learner, learner)
    expect_length(risk, nrow(new_x))
    expect_true(all(is.finite(risk)))
    expect_true(all(probability$probability >= 0 &
      probability$probability <= 1))
    expect_equal(replay$risk, risk, tolerance = 1e-14)
    expect_equal(
      replay$survival_probability,
      probability$probability,
      tolerance = 1e-14
    )
  }
})

test_that("additional survival states retain only their required engine form", {
  skip_if_not_installed("rpart", minimum_version = "4.1")
  skip_if_not_installed("partykit", minimum_version = "1.2")
  skip_if_not_installed("coin", minimum_version = "1.1")

  d <- make_sb_fixture(n = 48L, seed = 7931L)
  x <- d[c("gene_a", "gene_b", "gene_c")]

  survreg_adapter <- sb_additional_survival_learner_adapter(
    "survreg_weibull"
  )
  survreg_model <- survreg_adapter$fit(x, d$time, d$event)
  survreg_state <- survreg_adapter$extract_state(survreg_model)
  expect_false("engine_raw" %in% names(survreg_state))
  expect_true(all(is.finite(c(
    survreg_state$intercept,
    survreg_state$coefficient,
    survreg_state$scale
  ))))

  rpart_adapter <- sb_additional_survival_learner_adapter(
    "rpart_survival"
  )
  rpart_model <- rpart_adapter$fit(x, d$time, d$event)
  rpart_state <- rpart_adapter$extract_state(rpart_model)
  rpart_engine <- unserialize(rpart_state$engine_raw)
  expect_identical(
    intersect(c("model", "x", "y", "where"), names(rpart_engine)),
    character()
  )
  expect_identical(environment(rpart_engine$terms), baseenv())

  rsf_params <- sb_additional_survival_params("rsf_survival")
  rsf_params$num_trees <- 10L
  rsf_adapter <- sb_additional_survival_learner_adapter(
    "rsf_survival",
    rsf_params
  )
  rsf_state <- rsf_adapter$extract_state(rsf_adapter$fit(
    x,
    d$time,
    d$event
  ))
  expect_identical(
    rsf_state$training_data_retention,
    "event_time_metadata_required_by_randomForestSRC_forest"
  )
  rsf_engine <- unserialize(rsf_state$engine_raw)
  expect_identical(intersect(c("xvar", "yvar"), names(rsf_engine)), character())

  for (learner in c("ctree_survival", "cforest_survival")) {
    params <- sb_additional_survival_params(learner)
    if (identical(learner, "cforest_survival")) {
      params$num_trees <- 10L
    }
    adapter <- sb_additional_survival_learner_adapter(learner, params)
    state <- adapter$extract_state(adapter$fit(x, d$time, d$event))
    expect_identical(
      state$training_data_retention,
      "required_by_partykit_conditional_distribution_state"
    )
  }

  if (requireNamespace("BART", quietly = TRUE)) {
    params <- sb_additional_survival_params("bart_survival")
    params$time_grid_size <- 20L
    params$posterior_draws <- 30L
    params$burn_in <- 15L
    adapter <- sb_additional_survival_learner_adapter(
      "bart_survival",
      params
    )
    state <- adapter$extract_state(adapter$fit(x, d$time, d$event))
    engine <- unserialize(state$engine_raw)
    expect_true(all(names(engine) %in% c(
      "treedraws", "rm.const", "binaryOffset", "offset", "type", "times", "K"
    )))
    expect_false(any(c("yhat.train", "tx.train", "prob.train") %in%
      names(engine)))
  }
})

test_that("RSF fit and prediction restore the caller core option", {
  skip_if_not_installed("randomForestSRC", minimum_version = "3.3.0")

  d <- make_sb_fixture(n = 42L, seed = 7932L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  params <- sb_additional_survival_params("rsf_survival")
  params$num_trees <- 10L
  adapter <- sb_additional_survival_learner_adapter(
    "rsf_survival",
    params
  )
  withr::local_options(rf.cores = 7L, mc.cores = 6L)
  original_rfsrc <- randomForestSRC::rfsrc
  original_predict <- stats::predict
  fit_options <- integer()
  prediction_options <- integer()

  model <- testthat::with_mocked_bindings(
    rfsrc = function(...) {
      fit_options <<- c(
        rf.cores = getOption("rf.cores"),
        mc.cores = getOption("mc.cores")
      )
      original_rfsrc(...)
    },
    adapter$fit(x, d$time, d$event),
    .package = "randomForestSRC"
  )
  expect_identical(fit_options, c(rf.cores = 1L, mc.cores = 1L))
  expect_identical(
    c(rf.cores = getOption("rf.cores"), mc.cores = getOption("mc.cores")),
    c(rf.cores = 7L, mc.cores = 6L)
  )

  testthat::with_mocked_bindings(
    predict = function(...) {
      prediction_options <<- c(
        rf.cores = getOption("rf.cores"),
        mc.cores = getOption("mc.cores")
      )
      original_predict(...)
    },
    adapter$predict(model, x[1:3, , drop = FALSE]),
    .package = "stats"
  )
  expect_identical(
    prediction_options,
    c(rf.cores = 1L, mc.cores = 1L)
  )
  expect_identical(
    c(rf.cores = getOption("rf.cores"), mc.cores = getOption("mc.cores")),
    c(rf.cores = 7L, mc.cores = 6L)
  )
})

test_that("additional survival routes run CV and final prediction", {
  skip_if_not_installed("rpart", minimum_version = "4.1")
  skip_if_not_installed("randomForestSRC", minimum_version = "3.3.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")
  skip_if_not_installed("coin", minimum_version = "1.1")
  skip_if_not_installed("BART", minimum_version = "2.9.0")

  d <- make_sb_fixture(n = 60L, seed = 7941L)
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
      if (identical(learner, "rsf_survival")) {
        params$num_trees <- 20L
      }
      if (identical(learner, "cforest_survival")) {
        params$num_trees <- 20L
      }
      if (identical(learner, "bart_survival")) {
        params$time_grid_size <- 20L
        params$posterior_draws <- 40L
        params$burn_in <- 20L
      }
      if (identical(learner, "aorsf_survival")) {
        params$num_trees <- 20L
      }
      params
    }
  )

  for (learner in sb_additional_survival_learner_ids()) {
    pipeline <- sb_pipeline(learner = learner)
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 7942L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      d[c(6L, 31L, 52L), c("patient_id", "gene_c", "gene_a", "gene_b")],
      "patient_id"
    )

    expect_identical(result$status, "success")
    expect_true(result$complete_oof)
    expect_true(result$complete_survival_probability_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})
})

# ---- ahaz-bujar-adapters ----
local({
make_ahaz_bujar_fixture <- function(
  n = 90L, p = 6L, seed = 601L, tied = FALSE,
  heavy_censoring = FALSE, signal = c(0.8, -0.6, 0.35)
) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    beta <- rep(0, p)
    beta[seq_len(min(length(signal), p))] <- signal[
      seq_len(min(length(signal), p))
    ]
    risk <- drop(x %*% beta)
    event_time <- stats::rexp(n, exp(risk) / 10)
    censor_time <- stats::rexp(
      n,
      if (heavy_censoring) 0.16 else 0.04
    )
    time <- pmin(event_time, censor_time)
    if (tied) {
      time <- pmax(0.25, round(time * 2) / 2)
    } else {
      time <- time + seq_len(n) * 1e-10
    }
    names <- c(
      "gene one", "2gene", "if", "gene-4",
      paste0("gene_", seq_len(max(0L, p - 4L)))
    )[seq_len(p)]
    data.frame(
      patient_id = paste0("PRIVATE_PATIENT_", seq_len(n)),
      time = time,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    ) |>
      stats::setNames(c("patient_id", "time", "event", names))
  })
}

ahaz_bujar_features <- function(data) {
  setdiff(names(data), c("patient_id", "time", "event"))
}

test_that("four additive-hazard and Buckley-James learners add 32 routes", {
  ids <- c(sb_ahaz_survival_learner_ids(), sb_bujar_linear_learner_ids())
  expect_identical(ids, c(
    "ahaz_lin_ying",
    "ahaz_lasso_bic",
    "ahaz_sscad_bic",
    "bujar_linear_boost"
  ))
  rows <- sb_pipeline_registry()
  rows <- rows[rows$learner %in% ids, , drop = FALSE]
  expect_identical(nrow(rows), 4L * length(expected_fixed_learner_selectors("survival")))
  for (learner in ids) {
    expect_setequal(
      rows$selector[rows$learner == learner],
      expected_fixed_learner_selectors("survival")
    )
  }
  expect_true(all(
    rows$tuning_semantics[
      rows$learner %in% c("ahaz_lasso_bic", "ahaz_sscad_bic")
    ] == "analysis_fold_bic"
  ))
  expect_true(all(
    rows$tuning_semantics[
      rows$learner %in% c(
        "ahaz_lin_ying",
        "bujar_linear_boost"
      )
    ] == "none"
  ))
  expect_route_dependencies(rows, c(ahaz_lin_ying = "ahaz",
    ahaz_lasso_bic = "ahaz", ahaz_sscad_bic = "ahaz", bujar_linear_boost = "bujar"))
})

test_that("ahaz learners reproduce native coefficients and BIC choices", {
  skip_if_not_installed("ahaz", minimum_version = "1.15.1")
  data <- make_ahaz_bujar_fixture()
  features <- ahaz_bujar_features(data)
  x <- data[features]
  assessment <- x[c(2L, 31L, 76L), , drop = FALSE]
  engine_x <- sb_ahaz_engine_x(x)
  surv <- survival::Surv(data$time, data$event)

  for (learner in sb_ahaz_survival_learner_ids()) {
    adapter <- sb_ahaz_survival_learner_adapter(learner)
    model <- adapter$fit(x, data$time, data$event)
    native <- sb_fit_ahaz_engine(surv, engine_x, learner, adapter$params)
    native_coefficient <- unname(as.numeric(stats::coef(native)))
    expect_equal(model$coefficient, native_coefficient, tolerance = 0)
    expect_equal(
      adapter$predict(model, assessment),
      unname(as.vector(as.matrix(assessment) %*% native_coefficient)),
      tolerance = 0
    )
    if (identical(learner, "ahaz_lin_ying")) {
      expect_true(is.na(model$audit$selected_lambda))
    } else {
      expect_equal(
        model$audit$selected_lambda,
        as.numeric(native$lambda.min),
        tolerance = 0
      )
      expect_identical(
        model$audit$path_length,
        as.integer(length(native$lambda))
      )
    }
  }
})

test_that("ahaz horizon output is analysis-fold Cox calibrated and portable", {
  skip_if_not_installed("ahaz", minimum_version = "1.15.1")
  data <- make_ahaz_bujar_fixture(seed = 602L)
  features <- ahaz_bujar_features(data)
  x <- data[features]
  rownames(x) <- data$patient_id
  assessment <- x[c(3L, 40L, 81L), , drop = FALSE]
  horizon <- as.numeric(stats::median(data$time[data$event == 1L]))

  for (learner in sb_ahaz_survival_learner_ids()) {
    adapter <- sb_ahaz_survival_learner_adapter(learner)
    model <- adapter$fit(x, data$time, data$event)
    probability <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      assessment,
      horizon
    )
    expect_identical(
      probability$source,
      paste0(learner, "_analysis_cox_breslow_calibration")
    )
    expect_identical(probability$ties, "breslow")
    expect_gte(probability$calibration_coefficient, 0)
    expect_true(all(probability$probability >= 0 &
      probability$probability <= 1))

    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = horizon,
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, assessment)
    expect_equal(
      replay$risk,
      adapter$predict(model, assessment),
      tolerance = 0
    )
    expect_equal(
      replay$survival_probability,
      probability$probability,
      tolerance = 0
    )
    expect_false(any(c(
      "engine", "x", "X", "y", "outcome", "time", "event"
    ) %in% names(state)))
    expect_length(
      grepRaw(
        "PRIVATE_PATIENT_",
        serialize(state, NULL),
        fixed = TRUE,
        all = TRUE
      ),
      0L
    )
  }
})

test_that("ahaz resolves event ties while retaining other native limitations", {
  skip_if_not_installed("ahaz", minimum_version = "1.15.1")
  tied <- make_ahaz_bujar_fixture(tied = TRUE, seed = 603L)
  tied_x <- tied[ahaz_bujar_features(tied)]
  original <- tied
  for (learner in sb_ahaz_survival_learner_ids()) {
    adapter <- sb_ahaz_survival_learner_adapter(learner)
    model <- adapter$fit(tied_x, tied$time, tied$event)
    fit_time <- sb_ahaz_fit_times(tied$time, tied$event)
    native <- sb_fit_ahaz_engine(
      survival::Surv(fit_time$time, tied$event),
      sb_ahaz_engine_x(tied_x), learner, adapter$params
    )
    expect_equal(model$coefficient, unname(as.numeric(stats::coef(native))))
    expect_equal(all(is.finite(adapter$predict(model, tied_x))), TRUE)
    expect_identical(tied, original)
    expect_identical(model$audit$time_ties, fit_time$audit)
    risk <- adapter$predict(model, tied_x)
    horizon <- stats::median(tied$time)
    expected <- sb_fit_ahaz_calibration(tied$time, tied$event, risk, horizon)
    probability <- sb_predict_ahaz_survival_probability(
      model, tied_x, tied$time, tied$event, tied_x, horizon
    )
    expect_equal(probability$calibration_coefficient, expected$coefficient)
    expect_equal(probability$cumulative_baseline_hazard,
      expected$cumulative_baseline_hazard)
  }

  high_dimensional <- make_ahaz_bujar_fixture(
    n = 45L,
    p = 80L,
    seed = 704L
  )
  high_x <- high_dimensional[ahaz_bujar_features(high_dimensional)]
  expect_error(
    sb_ahaz_survival_learner_adapter("ahaz_lin_ying")$fit(
      high_x,
      high_dimensional$time,
      high_dimensional$event
    ),
    class = "stormbreaker_invalid_ahaz_fit"
  )

  expect_error(
    sb_fit_ahaz_calibration(
      seq_len(20),
      rep(c(1L, 0L), 10L),
      rep(0, 20L),
      10
    ),
    class = "stormbreaker_ahaz_calibration_failed"
  )
})

test_that("Bujar linear boosting is deterministic, portable, and side-effect free", {
  skip_if_not_installed("bujar", minimum_version = "0.2-11")
  data <- make_ahaz_bujar_fixture(seed = 604L)
  features <- ahaz_bujar_features(data)
  x <- data[features]
  rownames(x) <- data$patient_id
  assessment <- x[c(4L, 38L, 79L), , drop = FALSE]
  adapter <- sb_bujar_linear_learner_adapter()

  set.seed(9604L)
  caller_rng <- .Random.seed
  withr::local_dir(withr::local_tempdir())
  model <- adapter$fit(x, data$time, data$event)
  expect_identical(.Random.seed, caller_rng)
  expect_false(file.exists("file"))
  link <- unname(
    model$intercept + as.vector(as.matrix(assessment) %*% model$coefficient)
  )
  expect_equal(adapter$predict(model, assessment), -link, tolerance = 0)

  horizon <- as.numeric(stats::median(data$time[data$event == 1L]))
  probability <- adapter$predict_survival(
    model,
    x,
    data$time,
    data$event,
    assessment,
    horizon
  )
  expect_identical(
    probability$source,
    "bujar_linear_boost_analysis_residual_kaplan_meier"
  )
  expect_identical(probability$ties, "not_applicable")
  expect_true(all(probability$probability >= 0 &
    probability$probability <= 1))

  state <- unserialize(serialize(adapter$extract_state(model, list(
    x = x,
    time = data$time,
    event = data$event,
    prediction_horizon = horizon,
    time_unit = "month"
  )), NULL))
  replay <- adapter$predict_state(state, assessment)
  expect_equal(replay$risk, -link, tolerance = 0)
  expect_equal(
    replay$survival_probability,
    probability$probability,
    tolerance = 0
  )
  expect_false(any(c(
    "engine", "x", "X", "y", "outcome", "time", "event"
  ) %in% names(state)))
  expect_length(
    grepRaw(
      "PRIVATE_PATIENT_",
      serialize(state, NULL),
      fixed = TRUE,
      all = TRUE
    ),
    0L
  )
})

test_that("Bujar nonconvergence is a typed failure without a fallback", {
  skip_if_not_installed("bujar", minimum_version = "0.2-11")
  data <- make_ahaz_bujar_fixture(seed = 603L)
  x <- data[ahaz_bujar_features(data)]
  expect_error(
    sb_bujar_linear_learner_adapter()$fit(x, data$time, data$event),
    class = "stormbreaker_bujar_nonconvergence"
  )
})

test_that("new optional dependency failures are explicit", {
  local_mocked_bindings(sb_ahaz_available = function() FALSE)
  expect_error(sb_require_ahaz_dependency(), class = "stormbreaker_missing_ahaz")
  local_mocked_bindings(sb_bujar_available = function() FALSE)
  expect_error(
    sb_require_bujar_dependency(),
    class = "stormbreaker_missing_bujar"
  )
})

test_that("ahaz and Bujar run through final outcome-blind prediction", {
  skip_if_not_installed("ahaz", minimum_version = "1.15.1")
  skip_if_not_installed("bujar", minimum_version = "0.2-11")
  data <- make_ahaz_bujar_fixture(n = 100L, p = 4L, seed = 609L)
  features <- ahaz_bujar_features(data)
  horizon <- as.numeric(stats::median(data$time[data$event == 1L]))
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    features,
    "month",
    prediction_horizon = horizon
  )
  new_data <- data[1:5, c(
    "patient_id", "time", "event", rev(features)
  )]
  altered <- new_data
  altered$time <- rev(altered$time) + 1000
  altered$event <- 1L - altered$event

  for (learner in c(
    sb_ahaz_survival_learner_ids(),
    sb_bujar_linear_learner_ids()
  )) {
    pipeline <- sb_pipeline(learner = learner)
    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      new_data,
      "patient_id"
    )
    altered_prediction <- sb_predict_final(fit, altered, "patient_id")
    expect_identical(prediction$status, "success")
    expect_equal(
      prediction$predictions$risk,
      altered_prediction$predictions$risk,
      tolerance = 0
    )
    expect_equal(
      prediction$predictions$survival_probability,
      altered_prediction$predictions$survival_probability,
      tolerance = 0
    )
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(
      prediction$predictions$survival_probability >= 0 &
        prediction$predictions$survival_probability <= 1
    ))
  }
})

test_that("ahaz BIC and Bujar execute inside outer analysis folds", {
  skip_if_not_installed("ahaz", minimum_version = "1.15.1")
  skip_if_not_installed("bujar", minimum_version = "0.2-11")
  data <- make_ahaz_bujar_fixture(n = 100L, p = 4L, seed = 609L)
  features <- ahaz_bujar_features(data)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    features,
    "month"
  )

  ahaz <- sb_run_cv(
    task,
    sb_pipeline(learner = "ahaz_lasso_bic"),
    folds = 3L,
    seed = 17L
  )
  expect_identical(ahaz$status, "success")
  expect_true(ahaz$complete_oof)
  expect_identical(
    ahaz$provenance$ahaz_version,
    as.character(utils::packageVersion("ahaz"))
  )
  expect_true(all(vapply(
    ahaz$coefficients,
    function(fold) is.finite(fold$lambda) && fold$lambda > 0,
    logical(1)
  )))

  bujar <- sb_run_cv(
    task,
    sb_pipeline(learner = "bujar_linear_boost"),
    folds = 3L,
    seed = 17L
  )
  expect_identical(bujar$status, "success")
  expect_true(bujar$complete_oof)
  expect_identical(
    bujar$provenance$bujar_version,
    as.character(utils::packageVersion("bujar"))
  )
})

test_that("ahaz tied-event CV retains original assessment endpoints", {
  skip_if_not_installed("ahaz", minimum_version = "1.15.1")
  data <- make_ahaz_bujar_fixture(n = 120L, p = 4L, seed = 609L, tied = TRUE)
  original <- data
  task <- sb_surv_task(data, "patient_id", "time", "event",
    ahaz_bujar_features(data), "month")
  for (learner in sb_ahaz_survival_learner_ids()) {
    result <- sb_run_cv(task, sb_pipeline(learner = learner), folds = 3L, seed = 17L)
    expect_identical(result$status, "success")
    expect_identical(result$complete_oof, TRUE)
    expect_identical(data, original)
    rows <- match(result$predictions$patient_id, data$patient_id)
    expect_equal(result$predictions$time, data$time[rows])
    expect_equal(result$predictions$event, data$event[rows])
  }
})

test_that("ahaz CV audit reconstructs calibrated horizon probabilities", {
  skip_if_not_installed("ahaz", minimum_version = "1.15.1")
  data <- make_ahaz_bujar_fixture(
    n = 160L,
    p = 6L,
    seed = 645L,
    signal = c(1.4, -1.1, 0.8)
  )
  features <- ahaz_bujar_features(data)
  horizon <- as.numeric(stats::median(data$time[data$event == 1L]))
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    features,
    "month",
    prediction_horizon = horizon
  )
  result <- sb_run_cv(
    task,
    sb_pipeline(learner = "ahaz_lasso_bic"),
    folds = 3L,
    seed = 17L
  )

  expect_identical(result$status, "success")
  expect_identical(result$complete_survival_probability_oof, TRUE)
  for (fold in seq_along(result$folds)) {
    rows <- result$predictions$fold == fold
    audit <- result$baseline_hazards[[fold]]
    coefficient <- audit$calibration_coefficient
    expect_length(coefficient, 1L)
    if (length(coefficient) != 1L) {
      next
    }
    reconstructed <- exp(
      -audit$cumulative_baseline_hazard *
        exp(coefficient * result$predictions$risk[rows])
    )
    expect_equal(
      result$predictions$survival_probability[rows],
      reconstructed,
      tolerance = 1e-14
    )
  }
})

for (warning_stage in c("fit", "prediction")) {
  test_that(paste("Bujar ordinary", warning_stage, "warnings preserve predictions"), {
    skip_if_not_installed("bujar")
    d <- make_ahaz_bujar_fixture(n = 40L, p = 3L)
    x <- d[ahaz_bujar_features(d)]
    adapter <- sb_bujar_linear_learner_adapter()
    reference <- adapter$fit(x, d$time, d$event)
    expected <- adapter$predict(reference, x)
    message <- paste("informational", warning_stage, "warning")
    if (warning_stage == "fit") {
      native_fit <- bujar::bujar
      local_mocked_bindings(bujar = function(...) {
        warning(message)
        native_fit(...)
      }, .package = "bujar")
    } else {
      native_predict <- stats::predict
      local_mocked_bindings(predict = function(...) {
        warning(message)
        native_predict(...)
      }, .package = "stats")
    }
    seen <- character()
    observed <- withCallingHandlers({
      fitted <- adapter$fit(x, d$time, d$event)
      adapter$predict(fitted, x)
    }, warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    expect_identical(observed, expected)
    expect_identical(seen, message)
  })
}

test_that("Bujar warnings do not bypass convergence or coefficient checks", {
  skip_if_not_installed("bujar")
  d <- make_ahaz_bujar_fixture(n = 40L, p = 3L)
  x <- d[ahaz_bujar_features(d)]
  adapter <- sb_bujar_linear_learner_adapter()
  native_fit <- bujar::bujar
  defect <- "contype"
  local_mocked_bindings(bujar = function(...) {
    engine <- native_fit(...)
    warning("informational fit warning")
    if (defect == "contype") engine$contype <- 1
    if (defect == "nonconv") engine$nonconv <- TRUE
    if (defect == "coefficient") engine$coef.bj[1] <- NA_real_
    engine
  }, .package = "bujar")
  for (defect in c("contype", "nonconv", "coefficient")) {
    error <- tryCatch(
      suppressWarnings(adapter$fit(x, d$time, d$event)),
      error = identity
    )
    expected <- if (defect == "coefficient") {
      "stormbreaker_invalid_bujar_fit"
    } else {
      "stormbreaker_bujar_nonconvergence"
    }
    expect_s3_class(error, expected)
  }
})
})

# ---- ahaz-ties ----
local({
test_that("ahaz fit-only offsets preserve endpoints, order and caller RNG", {
  time <- c(1, 2, 2, 2, 3, 3, 4)
  event <- c(1L, 1L, 1L, 0L, 0L, 0L, 1L)
  original <- time
  withr::local_seed(62)
  before <- .Random.seed
  got <- sb_ahaz_fit_times(time, event)
  expect_identical(time, original)
  expect_identical(.Random.seed, before)
  expect_identical(got, sb_ahaz_fit_times(time, event))
  expect_identical(got$time[event == 0L], time[event == 0L])
  expect_equal(anyDuplicated(got$time[event == 1L]), 0L)
  expect_equal(all(got$time <= time), TRUE)
  expect_equal(all(got$time[time == 2] > 1), TRUE)
  expect_equal(got$audit$tied_event_groups, 1L)
  expect_equal(got$audit$tied_events, 2L)
  expect_equal(got$audit$changed_events, 1L)
  expect_equal(got$audit$max_abs_shift > 0, TRUE)
})

test_that("ahaz leaves supported times unchanged and breaks all-event ties", {
  time <- c(1, 1, 2, 2)
  event <- c(1L, 0L, 1L, 0L)
  expect_identical(sb_ahaz_fit_times(time, event)$time, time)
  same <- sb_ahaz_fit_times(rep(10, 20), rep(1L, 20))
  expect_equal(anyDuplicated(same$time), 0L)
  expect_equal(all(same$time > 0 & same$time <= 10), TRUE)
  expect_equal(max(10 - same$time) <= 1e-7, TRUE)
  alternative <- sb_ahaz_fit_times(rep(10, 20), rep(1L, 20), seed = 10102L)
  expect_equal(identical(same$time, alternative$time), FALSE)
  for (scale in c(1e-10, 1, 1e10)) {
    x <- c(1, 2, 2, 2.0001) * scale
    z <- sb_ahaz_fit_times(x, rep(1L, 4))$time
    expect_equal(anyDuplicated(z), 0L)
    expect_equal(z[1] < min(z[2:3]) && max(z[2:3]) < z[4], TRUE)
  }
})

test_that("unrepresentable ahaz offsets fail rather than reorder distinct times", {
  time <- c(1, 1 + .Machine$double.eps, 1 + .Machine$double.eps)
  error <- tryCatch(sb_ahaz_fit_times(time, rep(1L, 3)), error = identity)
  expect_s3_class(error, "stormbreaker_ahaz_time_precision")
})
})

# ---- classical-survival-fixed-adapters ----
local({
test_that("classical survival pipeline settings are fixed and explicit", {
  coxboost <- sb_pipeline(learner = "coxboost")
  plsrcox <- sb_pipeline(learner = "plsrcox")
  superpc <- sb_pipeline(learner = "superpc")
  survivalsvm <- sb_pipeline(learner = "survivalsvm")

  expect_identical(
    coxboost$id,
    "none__coxboost__steps100__penaltyevent9"
  )
  expect_identical(plsrcox$id, "none__plsrcox__components2")
  expect_identical(
    superpc$id,
    "none__superpc__threshold1__components1"
  )
  expect_identical(survivalsvm$id, "none__survivalsvm__gamma1")
  expect_identical(coxboost$risk_formula, "x_beta")
  expect_identical(plsrcox$risk_reference, "zero")
  expect_identical(superpc$risk_direction, "higher_is_higher_hazard")
  expect_identical(survivalsvm$learner_params$kernel, "lin_kernel")

  expect_error(
    sb_pipeline(learner = "coxboost", coxboost_penalty = 0),
    class = "stormbreaker_bad_coxboost_penalty"
  )
  expect_error(
    sb_pipeline(learner = "superpc", superpc_threshold = -1),
    class = "stormbreaker_bad_superpc_threshold"
  )
  expect_error(
    sb_pipeline(learner = "superpc", superpc_components = 4L),
    class = "stormbreaker_bad_superpc_components"
  )
  expect_error(
    sb_pipeline(learner = "survivalsvm", survivalsvm_gamma_mu = 0),
    class = "stormbreaker_bad_survivalsvm_gamma_mu"
  )
  expect_error(
    sb_pipeline(learner = "coxph", plsrcox_components = 1L),
    class = "stormbreaker_irrelevant_classical_survival_parameter"
  )
})

test_that("classical survival adapters preserve their native linear scores", {
  skip_if_not_installed("CoxBoost", minimum_version = "1.5.1")
  skip_if_not_installed("plsRcox", minimum_version = "1.8.0")
  skip_if_not_installed("superpc", minimum_version = "1.12")
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  d <- make_sb_fixture(n = 72L, seed = 7811L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(3L, 18L, 41L), , drop = FALSE]

  adapters <- list(
    coxboost = sb_coxboost_learner_adapter(stepno = 10L),
    plsrcox = sb_plsrcox_learner_adapter(components = 2L),
    superpc = sb_superpc_learner_adapter(
      threshold = 0,
      components = 1L
    ),
    survivalsvm = sb_survivalsvm_learner_adapter(gamma_mu = 1)
  )
  models <- lapply(adapters, function(adapter) {
    adapter$fit(x, d$time, d$event)
  })

  coxboost_native <- as.numeric(stats::predict(
    models$coxboost$engine,
    newdata = as.matrix(new_x),
    at.step = 10L,
    type = "lp"
  ))
  plsrcox_native <- as.numeric(stats::predict(
    models$plsrcox$engine,
    newdata = new_x,
    comps = 2L,
    type = "lp",
    verbose = FALSE
  ))
  superpc_native <- as.numeric(
    models$superpc$engine$prediction$v.pred.1df
  )
  svm_center <- models$survivalsvm$audit$center
  svm_scale <- models$survivalsvm$audit$scale
  scaled_new_x <- as.data.frame(
    sweep(sweep(as.matrix(new_x), 2L, svm_center, "-"),
          2L, svm_scale, "/"),
    optional = TRUE
  )
  names(scaled_new_x) <- names(new_x)
  survivalsvm_native <- -as.numeric(stats::predict(
    models$survivalsvm$engine,
    newdata = scaled_new_x
  )$predicted) * models$survivalsvm$audit$calibration_coefficient

  centered <- function(value) value - value[[1L]]
  expect_equal(
    centered(adapters$coxboost$predict(
      models$coxboost,
      new_x
    )),
    centered(coxboost_native),
    tolerance = 1e-10
  )
  expect_equal(
    centered(adapters$plsrcox$predict(
      models$plsrcox,
      new_x
    )),
    centered(plsrcox_native),
    tolerance = 1e-10
  )
  expect_equal(
    centered(adapters$superpc$predict(models$superpc, x)),
    centered(superpc_native),
    tolerance = 1e-10
  )
  expect_equal(
    centered(adapters$survivalsvm$predict(
      models$survivalsvm,
      new_x
    )),
    centered(survivalsvm_native),
    tolerance = 1e-8
  )

  for (learner in names(adapters)) {
    adapter <- adapters[[learner]]
    model <- models[[learner]]
    risk <- adapter$predict(model, new_x)
    probability <- adapter$predict_survival(
      model,
      x,
      d$time,
      d$event,
      new_x,
      prediction_horizon = 7
    )
    state <- adapter$extract_state(model, list(
      x = x,
      time = d$time,
      event = d$event,
      prediction_horizon = 7,
      time_unit = "month"
    ))
    restored <- unserialize(serialize(state, NULL))
    state_prediction <- adapter$predict_state(restored, new_x)

    expect_identical(state$learner, learner)
    expect_false("engine" %in% names(state))
    expect_identical(state$audit, model$audit)
    expect_equal(state_prediction$risk, risk, tolerance = 0)
    expect_equal(
      state_prediction$survival_probability,
      probability$probability,
      tolerance = 0
    )
    expect_true(all(probability$probability >= 0 &
      probability$probability <= 1))
    expect_identical(probability$ties, "breslow")
  }
})

test_that("classical survival learners run CV and portable final prediction", {
  skip_if_not_installed("CoxBoost", minimum_version = "1.5.1")
  skip_if_not_installed("plsRcox", minimum_version = "1.8.0")
  skip_if_not_installed("superpc", minimum_version = "1.12")
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  d <- make_sb_fixture(n = 72L, seed = 7821L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )
  pipelines <- list(
    sb_pipeline(learner = "coxboost", coxboost_stepno = 10L),
    sb_pipeline(learner = "plsrcox", plsrcox_components = 1L),
    sb_pipeline(
      learner = "superpc",
      superpc_threshold = 0,
      superpc_components = 1L
    ),
    sb_pipeline(learner = "survivalsvm", survivalsvm_gamma_mu = 1)
  )
  version_fields <- c(
    "coxboost_version", "plsrcox_version", "superpc_version",
    "survivalsvm_version"
  )

  for (index in seq_along(pipelines)) {
    pipeline <- pipelines[[index]]
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 7822L)
    fit <- sb_fit_final(task, pipeline)
    new_data <- d[c(5L, 29L, 61L), c("patient_id", predictors)]
    prediction <- sb_predict_final(fit, new_data, "patient_id")
    direct <- sb_predict_cox_learner_state(
      fit$learner_state,
      new_data[fit$required_features]
    )
    explanation <- sb_additive_shap(
      fit,
      new_data,
      "patient_id",
      d[1:24, fit$required_features, drop = FALSE]
    )

    expect_identical(result$status, "success")
    expect_true(result$complete_oof)
    expect_true(result$complete_survival_probability_oof)
    expect_length(result$coefficients, 3L)
    expect_true(all(vapply(
      result$coefficients,
      function(value) identical(value$status, "success"),
      logical(1)
    )))
    expect_true(version_fields[[index]] %in% names(result$provenance))
    expect_identical(fit$status, "success")
    expect_true(version_fields[[index]] %in% names(fit$provenance))
    expect_identical(prediction$status, "success")
    expect_equal(prediction$predictions$risk, direct$risk, tolerance = 0)
    expect_equal(
      prediction$predictions$survival_probability,
      direct$survival_probability,
      tolerance = 0
    )
    expect_identical(explanation$learner, pipeline$learner)
    expect_equal(
      explanation$scores$score,
      prediction$predictions$risk,
      tolerance = 1e-14
    )
  }
})

test_that("classical survival fit warnings are explicit failures", {
  skip_if_not_installed("CoxBoost", minimum_version = "1.5.1")
  skip_if_not_installed("plsRcox", minimum_version = "1.8.0")
  skip_if_not_installed("superpc", minimum_version = "1.12")
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  d <- make_sb_fixture(n = 36L, seed = 7829L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  injected_warning <- function(...) {
    warning("INJECTED_CLASSICAL_SURVIVAL_WARNING")
  }

  coxboost_error <- capture_sb_error(testthat::with_mocked_bindings(
    CoxBoost = injected_warning,
    sb_coxboost_learner_adapter(5L)$fit(x, d$time, d$event),
    .package = "CoxBoost"
  ))
  expect_s3_class(coxboost_error, "stormbreaker_coxboost_fit_warning")
  expect_match(
    conditionMessage(coxboost_error),
    "INJECTED_CLASSICAL_SURVIVAL_WARNING",
    fixed = TRUE
  )
  expect_s3_class(coxboost_error$parent, "warning")

  plsrcox_error <- capture_sb_error(testthat::with_mocked_bindings(
    plsRcox = injected_warning,
    sb_plsrcox_learner_adapter(1L)$fit(x, d$time, d$event),
    .package = "plsRcox"
  ))
  expect_s3_class(plsrcox_error, "stormbreaker_plsrcox_fit_warning")
  expect_match(
    conditionMessage(plsrcox_error),
    "INJECTED_CLASSICAL_SURVIVAL_WARNING",
    fixed = TRUE
  )
  expect_s3_class(plsrcox_error$parent, "warning")

  superpc_train_error <- capture_sb_error(testthat::with_mocked_bindings(
    `superpc.train` = injected_warning,
    sb_superpc_learner_adapter(0, 1L)$fit(x, d$time, d$event),
    .package = "superpc"
  ))
  expect_s3_class(
    superpc_train_error,
    "stormbreaker_superpc_train_warning"
  )
  expect_s3_class(superpc_train_error$parent, "warning")

  superpc_predict_error <- capture_sb_error(testthat::with_mocked_bindings(
    `superpc.predict` = injected_warning,
    sb_superpc_learner_adapter(0, 1L)$fit(x, d$time, d$event),
    .package = "superpc"
  ))
  expect_s3_class(
    superpc_predict_error,
    "stormbreaker_superpc_predict_warning"
  )
  expect_s3_class(superpc_predict_error$parent, "warning")

  survivalsvm_error <- capture_sb_error(testthat::with_mocked_bindings(
    survivalsvm = injected_warning,
    sb_survivalsvm_learner_adapter(1)$fit(x, d$time, d$event),
    .package = "survivalsvm"
  ))
  expect_s3_class(
    survivalsvm_error,
    "stormbreaker_survivalsvm_fit_warning"
  )
  expect_match(
    conditionMessage(survivalsvm_error),
    "INJECTED_CLASSICAL_SURVIVAL_WARNING",
    fixed = TRUE
  )
  expect_s3_class(survivalsvm_error$parent, "warning")

  calibration_error <- capture_sb_error(sb_fit_survivalsvm_calibration(
    time = 1:12,
    event = rep(1L, 12L),
    raw_risk = -(1:12),
    ties = "breslow"
  ))
  expect_s3_class(
    calibration_error,
    "stormbreaker_survivalsvm_calibration_warning"
  )
  expect_match(
    conditionMessage(calibration_error),
    "did not converge",
    fixed = TRUE
  )
  expect_s3_class(calibration_error$parent, "warning")
})

test_that("single-feature SuperPC preserves the native one-component score", {
  skip_if_not_installed("superpc")
  d <- sb_example_survival(n = 200L, seed = 20260905L)
  noise <- withr::with_seed(20260906L, matrix(rnorm(200L * 17L), nrow = 200L))
  for (j in seq_len(ncol(noise))) d[[paste0("noise_", j)]] <- noise[, j]
  d <- d[1:160, ]
  x <- d[setdiff(names(d), c("patient_id", "time", "event"))]
  adapter <- sb_superpc_learner_adapter(1, 1L)
  model <- adapter$fit(x, d$time, d$event)

  native_data <- list(x = t(as.matrix(x)), y = d$time, censoring.status = d$event)
  trained <- superpc::superpc.train(native_data, type = "survival", s0.perc = 0.5)
  selected <- abs(trained$feature.scores) >= 1
  expect_equal(sum(selected), 1L)
  # Oracle only: duplicate the same signal so native PCA retains matrix shape.
  # Its normalized rank-one score is unchanged; production adds no predictors.
  native_data$x <- rbind(native_data$x, native_data$x[selected, , drop = FALSE])
  trained$feature.scores <- c(trained$feature.scores, trained$feature.scores[selected])
  native <- superpc::superpc.predict(
    trained, native_data, native_data, threshold = 1,
    n.components = 1L, prediction.type = "continuous"
  )
  risk <- adapter$predict(model, x)
  expected <- as.numeric(native$v.pred.1df)
  expect_equal(risk - risk[1L], expected - expected[1L], tolerance = 1e-10)
  expect_identical(model$audit$selected_features, "gene_a")
  expect_identical(model$params$threshold, 1)
  expect_identical(model$params$components, 1L)
})

test_that("fixed SuperPC reports its native one-input training boundary", {
  skip_if_not_installed("superpc", minimum_version = "1.12")
  d <- make_sb_fixture(n = 72L, seed = 7830L)

  expect_error(
    sb_superpc_learner_adapter(1, 1L)$fit(
      d["gene_a"], d$time, d$event
    ),
    class = "stormbreaker_superpc_too_few_features"
  )
})

test_that("SuperPC native prediction errors explain the selected feature count", {
  skip_if_not_installed("superpc")
  d <- make_sb_fixture(n = 36L, seed = 7829L)
  local_mocked_bindings(
    `superpc.train` = function(...) {
      list(feature.scores = c(1.5, 0.5, 0.2))
    },
    `superpc.predict` = function(...) stop("INJECTED_ARRAY_ERROR"),
    .package = "superpc"
  )

  error <- capture_sb_error(
    sb_superpc_learner_adapter(1, 2L)$fit(
      d[c("gene_a", "gene_b", "gene_c")], d$time, d$event
    )
  )
  expect_s3_class(error, "stormbreaker_superpc_predict_error")
  expect_match(conditionMessage(error), "selecting 1 of 3 predictors", fixed = TRUE)
  expect_match(conditionMessage(error), "threshold=1, components=2", fixed = TRUE)
  expect_s3_class(error$parent, "simpleError")
  expect_identical(conditionMessage(error$parent), "INJECTED_ARRAY_ERROR")
})

test_that("missing classical survival dependencies fail without fallback", {
  d <- make_sb_fixture(n = 30L, seed = 7831L)
  x <- d[c("gene_a", "gene_b", "gene_c")]

  local_mocked_bindings(sb_coxboost_available = function() FALSE)
  expect_error(
    sb_coxboost_learner_adapter(5L)$fit(x, d$time, d$event),
    class = "stormbreaker_missing_coxboost"
  )

  local_mocked_bindings(sb_plsrcox_available = function() FALSE)
  expect_error(
    sb_plsrcox_learner_adapter(1L)$fit(x, d$time, d$event),
    class = "stormbreaker_missing_plsrcox"
  )

  local_mocked_bindings(sb_superpc_available = function() FALSE)
  expect_error(
    sb_superpc_learner_adapter(0, 1L)$fit(x, d$time, d$event),
    class = "stormbreaker_missing_superpc"
  )

  local_mocked_bindings(sb_survivalsvm_available = function() FALSE)
  expect_error(
    sb_survivalsvm_learner_adapter(1)$fit(x, d$time, d$event),
    class = "stormbreaker_missing_survivalsvm"
  )
})
})

# ---- mime-survival-cv-learners ----
local({
make_mime_cv_learner_fixture <- function(n = 90L, p = 12L,
                                         seed = 8901L) {
  set.seed(seed)
  x <- as.data.frame(matrix(stats::rnorm(n * p), nrow = n, ncol = p))
  names(x) <- sprintf("gene_%02d", seq_len(p))
  linear_predictor <- x[[1L]] - x[[2L]] + 0.6 * x[[3L]]
  event_time <- stats::rexp(n, rate = exp(linear_predictor) / 15)
  censor_time <- stats::rexp(n, rate = 1 / 25)
  list(
    x = x,
    time = pmax(pmin(event_time, censor_time), 0.05),
    event = as.integer(event_time <= censor_time)
  )
}

test_that("reference CV survival learners fit through their native engines", {
  skip_if_not_installed("CoxBoost")
  skip_if_not_installed("plsRcox")
  skip_if_not_installed("superpc")
  skip_if_not_installed("gbm")

  data <- make_mime_cv_learner_fixture()
  adapters <- list(
    coxboost_cv = sb_coxboost_cv_learner_adapter(
      inner_folds = 3L,
      seed = 8902L,
      maximum_steps = 80L
    ),
    plsrcox_cv = sb_plsrcox_cv_learner_adapter(
      inner_folds = 3L,
      seed = 8903L,
      maximum_components = 10L
    ),
    superpc_cv = sb_superpc_cv_learner_adapter(
      inner_folds = 3L,
      seed = 8904L
    ),
    gbm_cox_cv = sb_gbm_cox_cv_learner_adapter(
      inner_folds = 3L,
      seed = 8905L,
      maximum_trees = 60L
    )
  )

  models <- lapply(adapters, function(adapter) {
    model <- adapter$fit(data$x, data$time, data$event)
    observed <- adapter$predict(model, data$x[1:5, , drop = FALSE])
    state <- adapter$extract_state(model, list())
    restored <- adapter$predict_state(state, data$x[1:5, , drop = FALSE])
    expect_identical(adapter$id, state$learner)
    expect_equal(restored$risk, observed, tolerance = 0)
    expect_true(all(is.finite(observed)))
    model
  })

  expect_true(is.finite(models$coxboost_cv$audit$effective_penalty))
  expect_gte(models$coxboost_cv$audit$optimal_step, 1L)
  expect_lte(models$coxboost_cv$audit$optimal_step, 80L)
  expect_setequal(models$coxboost_cv$audit$fold_id, 1:3)
  expect_gte(models$plsrcox_cv$audit$selected_components, 1L)
  expect_lte(models$plsrcox_cv$audit$selected_components, 10L)
  expect_length(models$plsrcox_cv$audit$folds, 3L)
  expect_true(is.finite(models$superpc_cv$audit$selected_threshold))
  expect_identical(models$superpc_cv$audit$requested_folds, 3L)
  expect_gte(models$gbm_cox_cv$n_trees, 1L)
  expect_lte(models$gbm_cox_cv$n_trees, 60L)
  expect_length(models$gbm_cox_cv$tuning$cv_error, 60L)
})

test_that("reference CV learner controls reconstruct pipeline contracts", {
  for (learner in sb_analysis_tuned_survival_learner_ids()) {
    pipeline <- sb_pipeline(
      learner = learner,
      learner_tuning_folds = 3L,
      learner_seed = 8911L,
      gbm_cv_max_trees = if (identical(learner, "gbm_cox_cv")) {
        60L
      } else {
        NULL
      }
    )
    expect_identical(sb_lint_pipeline(pipeline)$status, "admitted")
    expect_identical(
      sb_validate_pipeline_contract(pipeline, "survival"),
      pipeline
    )
  }
})

test_that("PLS-Cox CV uses the sole feasible component for one predictor", {
  skip_if_not_installed("plsRcox")
  data <- make_mime_cv_learner_fixture(p = 3L)
  data$x <- data$x[1L]
  adapter <- sb_plsrcox_cv_learner_adapter(
    inner_folds = 3L, seed = 8912L, maximum_components = 10L
  )
  model <- adapter$fit(data$x, data$time, data$event)
  observed <- adapter$predict(model, data$x[1:5, , drop = FALSE])
  native <- plsRcox::plsRcox(
    data$x,
    time = data$time,
    event = data$event,
    nt = 1L,
    scaleX = TRUE,
    scaleY = FALSE,
    sparseStop = FALSE,
    allres = FALSE,
    verbose = FALSE
  )
  expected <- as.numeric(stats::predict(
    native, newdata = data$x[1:5, , drop = FALSE], type = "lp"
  ))
  expect_identical(model$audit$selected_components, 1L)
  expect_identical(
    model$audit$component_selection,
    "single_predictor_has_one_feasible_component"
  )
  expect_length(model$audit$folds, 0L)
  expect_true(all(is.finite(observed)))
  expect_equal(
    observed - mean(observed), expected - mean(expected), tolerance = 1e-12
  )
})

test_that("SuperPC CV reports its two-predictor tuning boundary", {
  skip_if_not_installed("superpc")
  data <- make_mime_cv_learner_fixture(p = 3L)
  data$x <- data$x[1L]
  error <- capture_sb_error(
    sb_superpc_cv_learner_adapter(inner_folds = 3L, seed = 8913L)$fit(
      data$x, data$time, data$event
    )
  )
  expect_s3_class(error, "stormbreaker_superpc_cv_too_few_features")
  expect_match(
    conditionMessage(error),
    "requires at least 2 predictors but received 1",
    fixed = TRUE
  )
})

test_that("PLS-Cox CV retains warnings beyond its selected component count", {
  skip_if_not_installed("plsRcox")
  data <- make_mime_cv_learner_fixture(n = 30L, p = 3L)
  message <- "Loglik converged before variable  3 ; coefficient may be infinite. "
  cv_calls <- 0L
  final_components <- NULL
  native_cv <- function(...) {
    cv_calls <<- cv_calls + 1L
    warning(message)
    warning(message)
    list(lambda.min5 = 2L, folds = split(seq_len(30L), rep(1:3, each = 10L)),
         cv.error5 = c(0.5, 0.6, 0.7, 0.65), cv.se5 = rep(0.01, 4L),
         completed.cv5 = matrix(TRUE, nrow = 4L, ncol = 3L))
  }
  native_fit <- function(Xplan, nt, ...) {
    final_components <<- nt
    list(Coeffs = matrix(0.1, nrow = ncol(Xplan), ncol = 1L,
                         dimnames = list(names(Xplan), NULL)), computed_nt = nt)
  }
  model <- testthat::with_mocked_bindings(
    sb_plsrcox_cv_learner_adapter(inner_folds = 3L, maximum_components = 3L)$fit(
      data$x, data$time, data$event
    ),
    cv.plsRcox = native_cv, plsRcox = native_fit, .package = "plsRcox"
  )
  expect_identical(cv_calls, 1L)
  expect_identical(final_components, 2L)
  expect_identical(model$audit$selected_components, 2L)
  expect_identical(model$audit$component_search_warnings, rep(message, 2L))
  expect_identical(attr(model, "stormbreaker_fit_warnings"), rep(message, 2L))
})

test_that("PLS-Cox CV still rejects selected and unclassified warnings", {
  skip_if_not_installed("plsRcox")
  data <- make_mime_cv_learner_fixture(n = 30L, p = 3L)
  messages <- c(
    "Loglik converged before variable  2 ; coefficient may be infinite. ",
    "Loglik converged before variable  2,3 ; coefficient may be infinite. ",
    "Loglik converged before variable  0 ; coefficient may be infinite. ",
    "Loglik converged before variable  4 ; coefficient may be infinite. ",
    "Ran out of iterations and did not converge",
    "INJECTED_UNKNOWN_PLSRCOX_WARNING"
  )
  for (message in messages) {
    error <- capture_sb_error(testthat::with_mocked_bindings(
      sb_plsrcox_cv_learner_adapter(inner_folds = 3L, maximum_components = 3L)$fit(
        data$x, data$time, data$event
      ),
      cv.plsRcox = function(...) {
        warning(message)
        list(lambda.min5 = 2L)
      },
      .package = "plsRcox"
    ))
    expect_s3_class(error, "stormbreaker_plsrcox_cv_warning")
    expect_match(conditionMessage(error), message, fixed = TRUE)
    expect_s3_class(error$parent, "warning")
  }
})

test_that("PLS-Cox CV does not accept incomplete selected folds or refit warnings", {
  skip_if_not_installed("plsRcox")
  data <- make_mime_cv_learner_fixture(n = 30L, p = 3L)
  message <- "Loglik converged before variable  3 ; coefficient may be infinite. "
  for (complete in c(FALSE, TRUE)) {
    error <- capture_sb_error(testthat::with_mocked_bindings(
      sb_plsrcox_cv_learner_adapter(inner_folds = 3L, maximum_components = 3L)$fit(
        data$x, data$time, data$event
      ),
      cv.plsRcox = function(...) {
        warning(message)
        completed <- matrix(TRUE, nrow = 4L, ncol = 3L)
        completed[3L, 1L] <- complete
        list(lambda.min5 = 2L, completed.cv5 = completed)
      },
      plsRcox = function(...) warning(message),
      .package = "plsRcox"
    ))
    expect_s3_class(error, if (complete) {
      "stormbreaker_plsrcox_cv_warning"
    } else {
      "stormbreaker_invalid_plsrcox_cv"
    })
    expect_match(conditionMessage(error), if (complete) {
      "final fit raised a warning"
    } else {
      "did not complete every fold"
    }, fixed = TRUE)
  }
})

test_that("PLS-Cox CV fits the original 20-feature case without changing selection", {
  skip_if_not_installed("plsRcox", minimum_version = "1.8.0")
  data <- sb_example_survival(n = 200L, seed = 20260905L)
  noise <- withr::with_seed(20260906L, matrix(rnorm(200L * 17L), nrow = 200L))
  for (j in seq_len(ncol(noise))) {
    data[[sprintf("gene_%02d", j + 3L)]] <- noise[, j]
  }
  features <- grep("^gene_", names(data), value = TRUE)
  training <- data[1:160, ]
  task <- sb_surv_task(training, "patient_id", "time", "event", features,
                       "month", prediction_horizon = 6)
  fit <- sb_fit_final(task, sb_pipeline(learner = "plsrcox_cv"))
  expect_identical(fit$status, "success")
  expect_identical(fit$learner_state$audit$selected_components, 4L)
  native_warnings <- character()
  native_cv <- withCallingHandlers(withr::with_seed(1L, plsRcox::cv.plsRcox(
    list(x = training[features], time = training$time, status = training$event),
    nfold = 5L, nt = 10L, plot.it = FALSE, folddetails = TRUE, verbose = FALSE,
    control = survival::coxph.control(eps = 1e-11, iter.max = 50L,
      toler.inf = survival::coxph.control()$toler.inf)
  )), warning = function(w) {
    native_warnings <<- c(native_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(as.integer(native_cv$lambda.min5), 4L)
  expect_identical(fit$learner_state$audit$component_search_warnings, native_warnings)
  expect_identical(fit$tuning$warnings, paste0("Learner refit: ", native_warnings))
  expect_match(fit$tuning$warnings, "coefficient may be infinite", fixed = TRUE)
  assessment <- data[161:200, c("patient_id", features)]
  prediction <- sb_predict_final(fit, assessment, "patient_id")
  expect_identical(prediction$status, "success")
  expect_identical(prediction$predictions$patient_id, assessment$patient_id)
  expect_identical(all(is.finite(prediction$predictions$risk)), TRUE)
  expect_identical(all(prediction$predictions$survival_probability >= 0 &
                         prediction$predictions$survival_probability <= 1), TRUE)
  native <- plsRcox::plsRcox(training[features], time = training$time,
                            event = training$event, nt = 4L, scaleX = TRUE,
                            scaleY = FALSE, sparseStop = FALSE,
                            allres = FALSE, verbose = FALSE)
  native_risk <- as.numeric(stats::predict(native,
    newdata = assessment[features], type = "lp"
  ))
  observed <- prediction$predictions$risk
  expect_equal(observed - mean(observed), native_risk - mean(native_risk),
               tolerance = 1e-10)
})
})

# ---- mime117-execution ----
local({
make_mime117_execution_fixture <- function(n = 150L, p = 15L,
                                           seed = 8801L) {
  set.seed(seed)
  x <- replicate(p, stats::rnorm(n))
  colnames(x) <- sprintf("gene_%02d", seq_len(p))
  linear_predictor <- 1.4 * x[, 1L] - 1.1 * x[, 2L] +
    0.8 * x[, 3L] - 0.5 * x[, 4L]
  event_time <- stats::rexp(n, rate = exp(linear_predictor) / 18)
  censor_time <- stats::rexp(n, rate = 1 / 26)
  data.frame(
    patient_id = sprintf("mime117_%03d", seq_len(n)),
    time = pmax(pmin(event_time, censor_time), 0.05),
    event = as.integer(event_time <= censor_time),
    x,
    check.names = FALSE
  )
}

test_that("representative MIME117 pipelines really fit and predict", {
  skip_if_not_installed("randomForestSRC")
  skip_if_not_installed("CoxBoost")
  skip_if_not_installed("glmnet")

  data <- make_mime117_execution_fixture()
  predictors <- grep("^gene_", names(data), value = TRUE)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    predictors,
    "month"
  )
  suite <- sb_mime117_pipelines(
    inner_folds = 3L,
    inner_seed = 8802L,
    selector_inner_folds = 3L,
    selector_seed = 8803L,
    rsf_nodesize = 5L,
    learner_tuning_folds = 3L,
    plsrcox_tuning_folds = 3L,
    learner_seed = 8804L,
    gbm_cv_max_trees = 50L
  )
  labels <- c(
    "RSF + CoxBoost",
    "CoxBoost + Enet[\u03b1=0.5]",
    "Lasso + StepCox[forward]",
    "StepCox[backward] + RSF"
  )

  fits <- lapply(labels, function(label) {
    index <- match(label, suite$official_model_label)
    expect_false(is.na(index), info = label)
    fit <- sb_fit_final(task, suite$pipeline[[index]])
    expect_identical(fit$status, "success", info = label)
    prediction <- sb_predict_final(
      fit,
      data[1:5, , drop = FALSE],
      "patient_id"
    )
    expect_identical(prediction$status, "success", info = label)
    expect_true(
      all(is.finite(prediction$predictions$risk)),
      info = label
    )
    fit
  })

  nested <- fits[[match("CoxBoost + Enet[\u03b1=0.5]", labels)]]
  expect_identical(
    nested$tuning$inner_selection_method,
    "glmnet_exclude_function"
  )
  expect_length(nested$tuning$inner_selections, 3L)
  expect_true(all(vapply(
    nested$tuning$inner_selections,
    function(selection) {
      length(selection$analysis_patient_ids) < nrow(data) &&
        length(selection$selected_features) > 1L
    },
    logical(1)
  )))
})
})

# ---- model-ahaz-bujar ----
local({
# ---- algorithm-ahaz ----

test_that("ahaz diagnostics preserve native input-unit coefficients and complete BIC paths", {
  skip_if_not_installed("ahaz")
  d <- sb_example_survival(160L, seed = 8612L)
  d$gene_c <- d$gene_c * 2 + 4
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  for (id in sb_ahaz_survival_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    x <- as.matrix(d[features])
    colnames(x) <- paste0("x", seq_along(features))
    native <- if (id == "ahaz_lin_ying") {
      ahaz::ahaz(survival::Surv(d$time, d$event), x, robust = FALSE)
    } else ahaz::tune.ahazpen(survival::Surv(d$time, d$event), x,
      standardize = TRUE, penalty = if (id == "ahaz_lasso_bic")
        ahaz::lasso.control(alpha = 1) else ahaz::sscad.control(a = 3.7, nsteps = 1L),
      tune = ahaz::bic.control(), dfmax = 3L, nlambda = 100L)
    co <- sb_ahaz_diagnostics(fit)
    expect_identical(co$table$feature, features)
    expect_equal(co$table$coefficient, as.numeric(stats::coef(native)), tolerance = 1e-12)
    expect_equal(as.numeric(x %*% co$table$coefficient),
      as.numeric(x %*% stats::coef(native)), tolerance = 1e-12)
    expect_null(fit$learner_state$diagnostic_path)
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_identical(sb_ahaz_diagnostics(compact), co)
    if (id != "ahaz_lin_ying") {
      bic <- sb_ahaz_diagnostics(fit, "bic_path")
      path <- sb_ahaz_diagnostics(fit, "coefficient_path")
      expect_equal(bic$table$lambda, native$lambda, tolerance = 0)
      expect_equal(bic$table$criterion, native$tunem, tolerance = 0)
      expect_equal(bic$table$df, native$df, tolerance = 0)
      expect_equal(bic$table$lambda[bic$table$selected], native$lambda.min, tolerance = 0)
      expect_equal(path$table$coefficient, as.numeric(t(as.matrix(native$fit$beta))), tolerance = 0)
      expect_equal(path$table$feature, rep(features, each = length(native$lambda)))
      expect_equal(path$table$lambda, rep(native$lambda, length(features)), tolerance = 0)
      expect_equal(nrow(path$table), 3L * length(native$lambda))
      expect_identical(sb_native_diagnostics(compact, complete = TRUE)$status$status,
        c("available", "missing_capture", "missing_capture"))
    }
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", if (id == "ahaz_lin_ying") 1L else 3L))
    for (view in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(view)
      built <- ggplot2::ggplot_build(p)
      expect_s3_class(built, "ggplot_built")
      expect_null(p$labels$caption)
      expect_null(p$labels$subtitle)
      if (view$type == "ahaz_coefficients") expect_equal(built$data[[2]]$x, view$table$coefficient, tolerance = 0)
      if (view$type == "ahaz_bic_path") expect_equal(built$data[[2]]$x,
        log10(view$table$lambda[view$table$selected]), tolerance = 0)
      if (view$type == "ahaz_coefficient_path") expect_equal(sort(built$data[[1]]$y), sort(view$table$coefficient), tolerance = 0)
    }
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})

# ---- algorithm-bujar ----

test_that("Bujar views preserve original coefficients and residual event steps", {
  skip_if_not_installed("bujar")
  d <- sb_example_survival(160L, seed = 8612L)
  features <- c("gene_c", "gene_a", "gene_b")
  x <- as.matrix(d[features])
  colnames(x) <- paste0("x", seq_along(features))
  native <- withr::with_seed(1L, bujar::bujar(y = log(d$time), cens = d$event, x = x,
    degree = 1, learner = "linear.regression", center = TRUE, mimpu = NULL,
    iter.bj = 20L, max.cycle = 10L, nu = .1, mstop = 50L, twin = FALSE,
    mstop2 = 100L, tuning = FALSE, cv = FALSE, vimpint = FALSE, rel.inf = FALSE,
    tol = 1e-4, n.cores = 1L, rng = 1L, trace = FALSE))
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
    prediction_horizon = stats::median(d$time))
  fit <- sb_fit_final(task, sb_pipeline(learner = "bujar_linear_boost"))
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  co <- sb_bujar_diagnostics(fit)
  expect_identical(co$features, features)
  expect_equal(co$table$coefficient, unname(native$coef.bj[-1]), tolerance = 0)
  expect_equal(co$intercept, unname(native$coef.bj[1]), tolerance = 0)
  lp <- as.numeric(stats::predict(native, newx = x))
  expect_equal(as.numeric(co$intercept + x %*% co$table$coefficient), lp, tolerance = 0)
  residual <- log(d$time) - lp
  km <- survival::survfit(survival::Surv(residual, d$event) ~ 1)
  keep <- km$n.event > 0
  curve <- sb_bujar_diagnostics(fit, "residual_survival")
  expect_equal(curve$table$residual, unname(km$time[keep]), tolerance = 0)
  expect_equal(curve$table$survival, unname(km$surv[keep]), tolerance = 0)
  expect_gt(sum(curve$table$residual < 0), 0)
  bundle <- sb_native_diagnostics(fit, complete = TRUE)
  expect_identical(bundle$status$status, c("available", "available"))
  for (v in bundle$diagnostics) {
    p <- sb_plot_algorithm_diagnostic(v)
    b <- ggplot2::ggplot_build(p)
    expect_s3_class(b, "ggplot_built")
    expect_null(p$labels$caption)
    expect_null(p$labels$subtitle)
    if (v$type == "bujar_residual_survival") {
      expect_equal(b$data[[1]]$x, v$table$residual, tolerance = 0)
      expect_equal(b$data[[1]]$y, v$table$survival, tolerance = 0)
    } else expect_equal(b$data[[2]]$x, v$table$coefficient, tolerance = 0)
  }
  expect_identical(serialize(fit, NULL, version = 2L), before)
  no_horizon <- sb_fit_final(sb_surv_task(d, "patient_id", "time", "event", features, "month"),
    sb_pipeline(learner = "bujar_linear_boost"))
  expect_identical(no_horizon$status, "success")
  expect_identical(sb_native_diagnostics(no_horizon, complete = TRUE)$status$status,
    c("available", "unavailable"))
})
})

# ---- model-akritas ----
local({
# ---- akritas-survival-adapter ----

test_that("Akritas uses the native empirical CDF and native tied-time survival engine", {
  skip_if_not_installed("distr6")
  skip_if_not_installed("survivalmodels")
  withr::local_seed(296)
  x <- as.data.frame(matrix(rnorm(80 * 3), 80, 3))
  names(x) <- c("gene A", "gene-B", "gene C")
  rownames(x) <- paste0("PRIVATE_", 1:80)
  time <- ceiling(rexp(80, exp(x[[1]]) / 20))
  event <- rbinom(80, 1, .7)
  for (p in c(1, 3)) {
    a <- x[seq_len(p)]
    before <- .Random.seed
    fit <- sb_fit_akritas_survival(a, time, event, "akritas_survival", sb_akritas_params())
    expect_identical(.Random.seed, before)
    query <- a[c(2, 3, 7), , drop = FALSE]
    dist <- if (p == 1) distr6::Empirical$new(as.matrix(a)) else distr6::EmpiricalMV$new(as.matrix(a))
    ord <- order(time, decreasing = TRUE)
    expected <- getFromNamespace("C_Akritas", "survivalmodels")(
      cbind(time, event)[ord, ], fit$grid, fit$grid,
      as.numeric(dist$cdf(data = as.matrix(a)))[ord],
      as.numeric(dist$cdf(data = as.matrix(query))), .5)
    expect_identical(sb_akritas_curves(fit, query), expected)
    horizon <- as.numeric(median(time))
    state <- sb_extract_akritas_state(fit, list(prediction_horizon = horizon, time = time, time_unit = "month"))
    expect_null(rownames(state$basis))
    expect_null(rownames(state$truth))
    expect_null(state$engine)
    replay <- sb_predict_akritas_state(unserialize(serialize(state, NULL)), query)
    expect_identical(replay$risk, sb_predict_akritas_risk(fit, query))
    expect_identical(replay$survival_probability,
      as.numeric(sb_akritas_curves(fit, query, horizon)))
    expected_rmst <- as.numeric(cbind(1, expected)[,seq_along(fit$grid),drop=FALSE] %*% diff(c(0,fit$grid)))
    expect_identical(replay$risk, -expected_rmst)
  }
})

test_that("Akritas computes each distinct query CDF once and restores row order", {
  skip_if_not_installed("distr6")
  skip_if_not_installed("survivalmodels")
  withr::local_seed(961)
  x <- data.frame("gene A" = seq_len(60), "gene-B" = seq_len(60)^2,
    check.names = FALSE)
  time <- ceiling(rexp(60, 0.1))
  event <- rep(c(1L, 0L, 1L), 20L)
  native <- getFromNamespace("C_Akritas", "survivalmodels")
  original_lookup <- utils::getFromNamespace
  observed <- list()
  local_mocked_bindings(getFromNamespace = function(x, ns, ...) {
    f <- original_lookup(x, ns, ...)
    if (identical(x, "C_Akritas") && identical(ns, "survivalmodels")) {
      return(function(...) {
        arguments <- list(...)
        observed[[length(observed) + 1L]] <<- arguments[[5L]]
        f(...)
      })
    }
    f
  }, .package = "utils")
  for (p in c(1L, 2L)) {
    fit <- sb_fit_akritas_survival(x[seq_len(p)], time, event,
      "akritas_survival", sb_akritas_params())
    for (indices in list(c(7, 2, 7, 55, 2, 14), c(55, 7, 14), 7)) {
      query <- x[indices, seq_len(p), drop = FALSE]
      cdf <- sb_akritas_empirical_cdf(fit$basis, as.matrix(query))
      for (times in list(c(0, 2.5, max(time) + 1), 2.5)) {
        expected <- native(fit$truth, as.numeric(times), fit$grid,
          fit$cdf_train, cdf, fit$params$lambda)
        before <- .Random.seed
        actual <- sb_akritas_curves(fit, query, times)
        expect_identical(actual, expected)
        expect_identical(tail(observed, 1L)[[1L]], unique(cdf))
        expect_identical(.Random.seed, before)
      }
    }
  }
})

test_that("Akritas runs public CV and final state-only survival predictions", {
  skip_if_not_installed("distr6"); skip_if_not_installed("survivalmodels")
  d <- sb_example_survival(90, 82)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
    prediction_horizon = as.numeric(median(d$time)))
  for (p in list(sb_pipeline("none", "akritas_survival"),
    sb_pipeline("univ_cox_top_k", "akritas_survival", selector_n_features = 2),
    sb_pipeline("fixed_features", "akritas_survival", selector_features = features[1:2], feature_set_id = "declared_two"))) {
    cv <- sb_run_cv(task, p, folds = 3, seed = 82)
    expect_identical(cv$status, "success")
    expect_identical(cv$complete_oof, TRUE)
    fit <- sb_fit_final(task, p)
    expect_identical(fit$status, "success")
    expect_identical(fit$provenance$distr6_version, as.character(packageVersion("distr6")))
    pred <- sb_predict_final(unserialize(serialize(fit, NULL)), d[c("patient_id", features)], "patient_id")
    expect_identical(pred$status, "success")
  }
})

# ---- algorithm-akritas ----

test_that("Akritas diagnostic neighborhoods and curves match the original engine", {
  skip_if_not_installed("distr6")
  skip_if_not_installed("survivalmodels")
  skip_if_not_installed("ggplot2")
  d <- data.frame(id = paste0("private-", 1:64), time = as.numeric(1:64),
    event = rep(c(1L, 0L), 32), gene_a = as.numeric(1:64), gene_b = sin(1:64))
  for (features in list("gene_a", c("gene_a", "gene_b"))) {
    task <- sb_surv_task(d, "id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = "akritas_survival"))
    expect_identical(fit$status, "success")
    file <- tempfile(fileext = ".rds")
    saveRDS(fit, file)
    restored <- readRDS(file)
    query <- d[c(1, 32, 64, 1), ]
    before <- serialize(restored, NULL, version = 2L)
    q <- sb_akritas_diagnostics(restored, "neighborhood", query)
    x <- as.matrix(d[features])
    empirical <- if (length(features) == 1L) distr6::Empirical$new(x) else distr6::EmpiricalMV$new(x)
    train_cdf <- as.numeric(empirical$cdf(data = x))
    query_cdf <- as.numeric(empirical$cdf(data = as.matrix(query[features])))
    neighbors <- vapply(query_cdf, function(z) sum(abs(z - train_cdf) <= 0.5), integer(1))
    events <- vapply(query_cdf, function(z) sum(d$event[abs(z - train_cdf) <= 0.5]), integer(1))
    expect_equal(q$table$empirical_cdf, query_cdf)
    expect_equal(q$table$neighbors, neighbors)
    expect_equal(q$table$events, events)
    expect_equal(q$table$censored, neighbors - events)
    if (length(features) == 1L) expect_equal(q$table$neighbors, c(33L, 64L, 33L, 33L))
    curves <- sb_akritas_diagnostics(restored, "survival_curves", query)
    times <- c(0, restored$learner_state$grid)
    expected <- getFromNamespace("C_Akritas", "survivalmodels")(
      restored$learner_state$truth, times, restored$learner_state$grid,
      restored$learner_state$cdf_train, query_cdf, 0.5)
    expect_equal(curves$table$survival, as.numeric(t(expected)))
    expect_equal(curves$table$time, rep(times, nrow(query)))
    expect_equal(curves$table$query, rep(seq_len(nrow(query)), each = length(times)))
    expect_identical(sb_predict_final(restored, query[c(1:3), ], "id"),
      sb_predict_final(fit, query[c(1:3), ], "id"))
    expect_identical(serialize(restored, NULL, version = 2L), before)
    for (item in list(q, curves)) {
      expect_s3_class(sb_plot_algorithm_diagnostic(item), "ggplot")
      expect_s3_class(ggplot2::ggplotGrob(sb_plot_algorithm_diagnostic(item)), "gtable")
      expect_identical(any(grepl("private-", unlist(item$table))), FALSE)
    }
    auto <- sb_native_diagnostics(restored, complete = TRUE, newdata = query)
    expect_identical(auto$status$status, c("available", "available"))
    expect_identical(sb_native_diagnostics(restored)$status$status, c("unavailable", "unavailable"))
    drawn <- plots(restored, new_data = query)
    expect_length(drawn$plots, 2L)
    expect_length(drawn$unavailable, 0L)
  }
})

test_that("the common survival workflow collects and exports diagnostics automatically", {
  skip_if_not_installed("distr6")
  skip_if_not_installed("survivalmodels")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  development <- surv_example(60, seed = 28)
  validation <- surv_example(36, seed = 29)
  for (methods in list("akritas", c("Cox", "akritas"))) {
    run <- surv_run(development, time = "time", event = "event",
      id = "patient_id", models = methods, folds = 2)
    finished <- finish_run(run, "akritas", external = list(Validation = validation))
    expect_identical(finished$fit$status, "success")
    expect_identical(finished$evaluations$Validation$status, "success")
    expect_length(finished$plots$algorithm, 2L)
    expect_length(finished$issues$algorithm, 0L)
    if (!is.null(run$details$akritas$fit)) {
      expect_identical(finished$fit, run$details$akritas$fit)
    }
    directory <- tempfile("common-workflow-")
    export_run(finished, directory)
    index <- read.csv(file.path(directory, "figures.csv"))
    algorithm_files <- index$file[grepl("/algorithm/", index$plot, fixed = TRUE)]
    expected_pages <- sb_plot_pages(finished$plots$algorithm)
    expect_length(algorithm_files, length(expected_pages))
    expect_identical(index$plot[grepl("/algorithm/", index$plot, fixed = TRUE)],
      paste0("/algorithm/", names(expected_pages)))
    expect_gt(min(file.info(file.path(directory, algorithm_files))$size), 0)
    saveRDS(finished$fit, file.path(directory, "model.rds"))
    restored <- readRDS(file.path(directory, "model.rds"))
    drawn <- plots(restored)
    expect_identical(names(drawn$plots), names(finished$plots$algorithm))
    expect_equal(lapply(drawn$plots, function(p) p$data),
      lapply(finished$plots$algorithm, function(p) p$data))
    expect_length(drawn$unavailable, 0L)
    expect_identical(predict_final(restored, validation, "patient_id"),
      predict_final(finished$fit, validation, "patient_id"))
    expect_equal(vapply(drawn$native$diagnostics, function(x) length(unique(x$table$query)),
      integer(1)), rep(60L, 2L), ignore_attr = TRUE)
  }
})
})

# ---- model-bagged-knn ----
local({
test_that("bagged KNN diagnostics preserve all original bags and native curves", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(100L, seed = 8612L)
  for (features in list(c("gene_c", "gene_a", "gene_b"), "gene_c")) {
    task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = "bagged_knn_survival"))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    e <- unserialize(fit$learner_state$engine_raw)
    bags <- methods::slot(e, "base_learners")
    z <- sb_bagged_knn_diagnostics(fit)
    expect_equal(nrow(z$table), 50L * nrow(d))
    for (i in seq_along(bags)) {
      index <- methods::slot(bags[[i]], "bootstrap_sample")
      expect_identical(z$table$multiplicity[z$table$bag == i],
        as.integer(table(factor(index, levels = seq_len(nrow(d))))))
    }
    query <- d[c(1L, 7L), features, drop = FALSE]
    query[[1L]] <- query[[1L]] + 0.17
    curves <- sb_bagged_knn_diagnostics(fit, "bag_curves", query)
    model <- sb_restore_bagged_knn_survival_model(fit$learner_state)
    original <- sb_predict_bagged_knn_survival_curves(model, query)
    for (i in seq_len(nrow(query))) {
      summary <- curves$summary[curves$summary$query == i, ]
      expect_equal(summary$time, original$time)
      expect_equal(summary$survival, unname(original$survival[i, ]), tolerance = 1e-12)
    }
    native_predict <- get("predict", asNamespace("bnnSurvival"))
    native <- native_predict(bags[[1L]], methods::slot(e, "train_data"),
      data.matrix(sb_bagged_knn_internal_data(query)$data), methods::slot(e, "timepoints"),
      methods::slot(e, "metric"), methods::slot(e, "weighting_function"), methods::slot(e, "k"))
    expect_equal(curves$table$survival[curves$table$bag == 1L], as.numeric(t(native)))
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(curves)), "ggplot_built")
    expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status,
      c("available", "unavailable"))
    expect_identical(sb_native_diagnostics(fit, complete = TRUE, newdata = query)$status$status,
      c("available", "available"))
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
})

# ---- model-flexsurv ----
local({
# ---- algorithm-flexsurv ----

test_that("flexsurv diagnostics preserve original parameter scales and spline baseline", {
  skip_if_not_installed("flexsurv")
  d <- sb_example_survival(160L, seed = 8612L)
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  for (id in c(sb_flexsurv_parametric_learner_ids(), sb_flexsurvspline_learner_ids())) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    engine <- fit$diagnostic_artifacts$learner$engine
    engine_features <- sprintf(".stormbreaker_x%04d", seq_along(features))
    numerical <- fit$learner_state$numerical_fitting
    original_beta <- unname(stats::coef(engine)[engine_features] / numerical$scale)
    location_shift <- sum(original_beta * numerical$center)
    view <- sb_flexsurv_diagnostics(fit)
    expect_identical(view$table$feature, features)
    expect_equal(view$table$coefficient, original_beta, tolerance = 0)
    spline <- id %in% sb_flexsurvspline_learner_ids()
    if (spline) {
      original_gamma <- stats::coef(engine)[names(view$gamma)]
      original_gamma[[1L]] <- original_gamma[[1L]] - location_shift
      expect_equal(view$gamma, original_gamma, tolerance = 0)
      expect_equal(view$knots, unname(engine$knots), tolerance = 0)
      baseline <- sb_flexsurv_diagnostics(fit, "spline_baseline")
      zero <- as.data.frame(as.list(stats::setNames(-numerical$center / numerical$scale, engine_features)))
      native <- get("summary.flexsurvreg", asNamespace("flexsurv"))(engine,
        newdata = zero, type = "survival", t = baseline$table$time, ci = FALSE)[[1]]
      expect_equal(baseline$table$survival, native$est, tolerance = 1e-12)
      expect_equal(range(log(baseline$table$time)), range(engine$knots), tolerance = 1e-12)
      expect_equal(sum(exp(engine$knots) %in% baseline$table$time), length(engine$knots))
      expect_identical(baseline$params, fit$learner_state$params)
      expect_identical(baseline$time_unit, "month")
      expect_identical(sb_plot_algorithm_diagnostic(baseline)$labels$x, "Time (month)")
    } else {
      location <- fit$learner_state$params$location_parameter
      expect_equal(view$intercept, unname(stats::coef(engine)[location]) - location_shift, tolerance = 0)
      expect_equal(view$ancillary, stats::setNames(engine$res[names(view$ancillary), "est"], names(view$ancillary)), tolerance = 0)
    }
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_identical(sb_flexsurv_diagnostics(compact), view)
    bundle <- sb_native_diagnostics(compact, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", if (spline) 2L else 1L))
    for (diagnostic in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(diagnostic)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_null(p$labels$subtitle)
      expect_null(p$labels$caption)
    }
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})

# ---- flexsurv-parametric-fixed ----

make_flexsurv_parametric_fixture <- function(
  n = 84L, p = 6L, seed = 9901L, tied = FALSE
) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n)
    true_risk <- 0.8 * x[, 1L]
    if (p >= 2L) {
      true_risk <- true_risk - 0.45 * x[, 2L]
    }
    event_time <- flexsurv::rgengamma(
      n,
      mu = 2.3 - true_risk,
      sigma = 0.65,
      Q = 0.55
    )
    censor_time <- stats::rexp(n, rate = 1 / 20)
    time <- pmin(event_time, censor_time)
    if (tied) {
      time <- pmax(round(time, 1), 0.1)
    }
    feature_names <- c(
      "gene one",
      if (p >= 2L) "2gene" else NULL,
      if (p >= 3L) paste0("gene-", seq_len(p - 2L)) else NULL
    )
    data.frame(
      patient_id = sprintf("FSP%03d", seq_len(n)),
      time = time,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    )
  })
}

flexsurv_parametric_expected_params <- function() {
  list(
    flexsurv_gamma = list(
      distribution = "gamma",
      parameterization = "shape_rate",
      location_parameter = "rate",
      risk_transform = "log_location"
    ),
    flexsurv_gompertz = list(
      distribution = "gompertz",
      parameterization = "shape_rate",
      location_parameter = "rate",
      risk_transform = "log_location"
    ),
    flexsurv_gengamma = list(
      distribution = "gengamma",
      parameterization = "stable_mu_sigma_Q",
      location_parameter = "mu",
      risk_transform = "negative_location"
    )
  )
}

test_that("flexsurv admits three nonduplicate parametric survival families", {
  expected <- flexsurv_parametric_expected_params()

  expect_identical(sb_flexsurv_parametric_learner_ids(), names(expected))
  for (learner in names(expected)) {
    expect_identical(
      sb_flexsurv_parametric_params(learner),
      expected[[learner]]
    )
  }
  expect_false(any(c(
    "flexsurv_genf_fixed", "flexsurv_gengamma_orig_fixed",
    "flexsurv_genf_orig_fixed", "flexsurv_weibull_fixed",
    "flexsurv_exponential_fixed", "flexsurv_lognormal_fixed",
    "flexsurv_loglogistic_fixed"
  ) %in% sb_flexsurv_parametric_learner_ids()))

  rows <- sb_pipeline_registry()
  rows <- rows[
    rows$learner %in% sb_flexsurv_parametric_learner_ids(), ,
    drop = FALSE
  ]
  expect_identical(nrow(rows), 3L * length(expected_fixed_learner_selectors("survival")))
  expect_setequal(
    unique(rows$selector),
    expected_fixed_learner_selectors("survival")
  )
  expect_true(all(grepl("flexsurv", rows$optional_dependency, fixed = TRUE)))
})

test_that("flexsurv risk and horizon survival equal native predictions", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  data <- make_flexsurv_parametric_fixture()
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[c(3L, 29L, 71L), , drop = FALSE]
  horizon <- as.numeric(stats::quantile(data$time, 0.45))

  for (learner in sb_flexsurv_parametric_learner_ids()) {
    adapter <- sb_flexsurv_parametric_learner_adapter(learner)
    set.seed(9902L)
    caller_rng <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, caller_rng)

    engine_assessment <- sb_flexsurv_numerical_predictors(assessment, model$numerical_fitting)
    names(engine_assessment) <- model$engine_features
    native_location <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_assessment,
      type = "link",
      conf.int = FALSE
    )$.pred_link)
    native_risk <- if (model$params$risk_transform == "log_location") {
      log(native_location)
    } else {
      -native_location
    }
    native_survival <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_assessment,
      type = "survival",
      times = horizon,
      conf.int = FALSE
    )$.pred_survival)

    risk <- adapter$predict(model, assessment)
    survival <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      assessment,
      horizon
    )
    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = horizon,
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, assessment)

    expect_equal(risk, native_risk, tolerance = 1e-12)
    expect_equal(survival$probability, native_survival, tolerance = 1e-12)
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      survival$probability,
      tolerance = 0
    )
    expect_identical(survival$ties, "not_applicable")
    expect_false("engine_raw" %in% names(state))
    expect_false(any(c("time", "event", "outcome") %in% names(state)))
    expect_identical(
      state$survival_probability$source,
      paste0(
        "flexsurv_", model$params$distribution,
        "_distribution_direct"
      )
    )
  }
})

test_that("flexsurv parametric families support one feature and tied times", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  cases <- list(
    one_feature = make_flexsurv_parametric_fixture(
      n = 64L, p = 1L, seed = 9903L
    ),
    tied_times = make_flexsurv_parametric_fixture(
      n = 84L, p = 5L, seed = 9904L, tied = TRUE
    )
  )
  for (data in cases) {
    features <- setdiff(names(data), c("patient_id", "time", "event"))
    x <- data[features]
    horizons <- as.numeric(stats::quantile(data$time, c(0.3, 0.7)))
    for (learner in sb_flexsurv_parametric_learner_ids()) {
      adapter <- sb_flexsurv_parametric_learner_adapter(learner)
      model <- adapter$fit(x, data$time, data$event)
      early <- adapter$predict_survival(
        model, x, data$time, data$event,
        x[1:5, , drop = FALSE], horizons[[1L]]
      )$probability
      late <- adapter$predict_survival(
        model, x, data$time, data$event,
        x[1:5, , drop = FALSE], horizons[[2L]]
      )$probability

      expect_true(all(is.finite(early)))
      expect_true(all(is.finite(late)))
      expect_true(all(early >= 0 & early <= 1))
      expect_true(all(late >= 0 & late <= 1))
      expect_true(all(early >= late))
    }
  }
})

test_that("flexsurv reports a nonfinite Hessian as a failed fit", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  data <- withr::with_seed(2105L, {
    n <- 44L
    p <- 52L
    x <- matrix(stats::rnorm(n * p), nrow = n)
    true_risk <- 0.8 * x[, 1L] - 0.45 * x[, 2L]
    event_time <- flexsurv::rgompertz(
      n,
      shape = 0.035,
      rate = exp(true_risk) / 14
    )
    censor_time <- stats::rexp(n, rate = 1 / 20)
    x <- as.data.frame(x)
    names(x) <- paste0("gene ", seq_len(p))
    list(
      x = x,
      time = pmin(event_time, censor_time),
      event = as.integer(event_time <= censor_time)
    )
  })

  expect_error(
    sb_flexsurv_parametric_learner_adapter(
      "flexsurv_gompertz"
    )$fit(data$x, data$time, data$event),
    class = "stormbreaker_flexsurv_nonconvergence"
  )
})

test_that("flexsurv parametric methods execute through CV and final prediction", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  data <- make_flexsurv_parametric_fixture(n = 78L, seed = 9906L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    features,
    "month",
    prediction_horizon = as.numeric(stats::quantile(data$time, 0.45))
  )
  new_data <- data[1:6, c("patient_id", rev(features)), drop = FALSE]

  for (learner in sb_flexsurv_parametric_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 3L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 9907L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      new_data,
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_true(comparison$complete_oof)
    expect_true(comparison$complete_survival_probability_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})

# ---- flexsurvspline-fixed ----

make_flexsurvspline_fixture <- function(
  n = 120L, p = 4L, seed = 4111L, tied = FALSE
) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    true_risk <- 0.75 * x[, 1L]
    if (p >= 2L) {
      true_risk <- true_risk - 0.45 * x[, 2L]
    }
    if (p >= 3L) {
      true_risk <- true_risk + 0.20 * x[, 3L]
    }
    event_time <- (
      -log(stats::runif(n)) / exp(true_risk)
    )^(1 / 1.35) * 8
    censor_time <- stats::rexp(n, rate = 1 / 16)
    time <- pmin(event_time, censor_time)
    if (tied) {
      time <- pmax(round(time * 2) / 2, 0.25)
    }
    feature_names <- c(
      "gene one",
      if (p >= 2L) "2gene" else NULL,
      if (p >= 3L) paste0("gene-", seq_len(p - 2L)) else NULL
    )
    data.frame(
      patient_id = sprintf("FSS%03d", seq_len(n)),
      time = time,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    )
  })
}

flexsurvspline_expected_params <- function() {
  list(
    flexsurvspline_hazard_k1 = list(
      scale = "hazard",
      k = 1L,
      timescale = "log",
      spline = "rp",
      knot_rule = "analysis_event_log_time_quantiles"
    ),
    flexsurvspline_odds_k1 = list(
      scale = "odds",
      k = 1L,
      timescale = "log",
      spline = "rp",
      knot_rule = "analysis_event_log_time_quantiles"
    ),
    flexsurvspline_normal_k1 = list(
      scale = "normal",
      k = 1L,
      timescale = "log",
      spline = "rp",
      knot_rule = "analysis_event_log_time_quantiles"
    )
  )
}

test_that("three nonduplicate one-knot Royston-Parmar scales are admitted", {
  expected <- flexsurvspline_expected_params()

  expect_identical(sb_flexsurvspline_learner_ids(), names(expected))
  for (learner in names(expected)) {
    expect_identical(
      sb_flexsurvspline_params(learner),
      expected[[learner]]
    )
  }
  expect_false(any(c(
    "flexsurvspline_hazard_k0_fixed",
    "flexsurvspline_odds_k0_fixed",
    "flexsurvspline_normal_k0_fixed",
    "flexsurvspline_hazard_k2_fixed",
    "flexsurvspline_hazard_splines2ns_k1_fixed"
  ) %in% sb_flexsurvspline_learner_ids()))

  rows <- sb_pipeline_registry()
  rows <- rows[
    rows$learner %in% sb_flexsurvspline_learner_ids(), ,
    drop = FALSE
  ]
  expect_identical(nrow(rows), 3L * length(expected_fixed_learner_selectors("survival")))
  expect_setequal(
    unique(rows$selector),
    expected_fixed_learner_selectors("survival")
  )
  expect_true(all(grepl("flexsurv", rows$optional_dependency, fixed = TRUE)))
})

test_that("spline risk, knots, and direct survival equal native predictions", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  data <- make_flexsurvspline_fixture()
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[c(2L, 37L, 91L), , drop = FALSE]
  horizons <- as.numeric(stats::quantile(data$time, c(0.35, 0.65)))
  expected_knots <- c(
    min(log(data$time[data$event == 1L])),
    as.numeric(stats::quantile(log(data$time[data$event == 1L]), 0.5)),
    max(log(data$time[data$event == 1L]))
  )

  for (learner in sb_flexsurvspline_learner_ids()) {
    adapter <- sb_flexsurvspline_learner_adapter(learner)
    set.seed(4112L)
    caller_rng <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, caller_rng)

    engine_assessment <- sb_flexsurv_numerical_predictors(assessment, model$numerical_fitting)
    names(engine_assessment) <- model$engine_features
    native_risk <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_assessment,
      type = "link",
      conf.int = FALSE
    )$.pred_link)
    native_early <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_assessment,
      type = "survival",
      times = horizons[[1L]],
      conf.int = FALSE
    )$.pred_survival)
    native_late <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_assessment,
      type = "survival",
      times = horizons[[2L]],
      conf.int = FALSE
    )$.pred_survival)

    risk <- adapter$predict(model, assessment)
    early <- adapter$predict_survival(
      model, x, data$time, data$event, assessment, horizons[[1L]]
    )
    late <- adapter$predict_survival(
      model, x, data$time, data$event, assessment, horizons[[2L]]
    )
    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = horizons[[1L]],
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, assessment)

    expect_equal(risk, native_risk, tolerance = 1e-12)
    expect_equal(early$probability, native_early, tolerance = 1e-12)
    expect_equal(late$probability, native_late, tolerance = 1e-12)
    expect_true(all(early$probability >= late$probability))
    expect_equal(model$prediction$knots, expected_knots, tolerance = 1e-12)
    expect_equal(state$knots, expected_knots, tolerance = 1e-12)
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      early$probability,
      tolerance = 0
    )
    expect_false(any(c(
      "engine", "engine_raw", "time", "event", "outcome", "data"
    ) %in% names(state)))
    expect_identical(early$ties, "not_applicable")
    expect_identical(
      state$survival_probability$source,
      paste0(
        "flexsurvspline_", model$params$scale,
        "_k1_rp_distribution_direct"
      )
    )
  }
})

test_that("spline knots are learned only from the supplied analysis fold", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  data <- make_flexsurvspline_fixture(n = 121L, seed = 4120L)
  data$time[[121L]] <- max(data$time) * 50
  data$event[[121L]] <- 1L
  analysis <- data[seq_len(120L), , drop = FALSE]
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  expected_analysis_knots <- c(
    min(log(analysis$time[analysis$event == 1L])),
    as.numeric(stats::quantile(
      log(analysis$time[analysis$event == 1L]), 0.5
    )),
    max(log(analysis$time[analysis$event == 1L]))
  )
  full_boundary <- max(log(data$time[data$event == 1L]))

  for (learner in sb_flexsurvspline_learner_ids()) {
    model <- sb_flexsurvspline_learner_adapter(learner)$fit(
      analysis[features], analysis$time, analysis$event
    )
    expect_equal(
      model$prediction$knots,
      expected_analysis_knots,
      tolerance = 1e-12
    )
    expect_false(isTRUE(all.equal(
      model$prediction$knots[[3L]], full_boundary,
      tolerance = 1e-12
    )))
  }
})

test_that("one-feature and tied-time spline fits remain portable", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  cases <- list(
    one_feature = make_flexsurvspline_fixture(
      n = 96L, p = 1L, seed = 4201L
    ),
    tied_times = make_flexsurvspline_fixture(
      n = 120L, p = 5L, seed = 4202L, tied = TRUE
    )
  )
  for (data in cases) {
    features <- setdiff(names(data), c("patient_id", "time", "event"))
    x <- data[features]
    horizons <- as.numeric(stats::quantile(data$time, c(0.3, 0.7)))
    for (learner in sb_flexsurvspline_learner_ids()) {
      adapter <- sb_flexsurvspline_learner_adapter(learner)
      model <- adapter$fit(x, data$time, data$event)
      early <- adapter$predict_survival(
        model, x, data$time, data$event,
        x[1:5, , drop = FALSE], horizons[[1L]]
      )$probability
      late <- adapter$predict_survival(
        model, x, data$time, data$event,
        x[1:5, , drop = FALSE], horizons[[2L]]
      )$probability
      expect_true(all(is.finite(adapter$predict(
        model, x[1:5, , drop = FALSE]
      ))))
      expect_true(all(early >= late))
    }
  }
})

test_that("one-knot spline methods execute through CV and final prediction", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  data <- make_flexsurvspline_fixture(n = 120L, seed = 4111L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    features,
    "month",
    prediction_horizon = as.numeric(stats::quantile(data$time, 0.45))
  )
  new_data <- data[1:6, c("patient_id", rev(features)), drop = FALSE]

  for (learner in sb_flexsurvspline_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 3L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 4113L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      new_data,
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_true(comparison$complete_oof)
    expect_true(comparison$complete_survival_probability_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})

# ---- survreg-distribution-fixed-adapters ----

test_that("fixed SurvReg learners declare four distribution assumptions", {
  expected <- list(
    survreg_weibull = list(
      distribution = "weibull",
      scale = "estimated"
    ),
    survreg_exponential = list(
      distribution = "exponential",
      scale = "fixed_1"
    ),
    survreg_lognormal = list(
      distribution = "lognormal",
      scale = "estimated"
    ),
    survreg_loglogistic = list(
      distribution = "loglogistic",
      scale = "estimated"
    )
  )

  expect_identical(sb_survreg_learner_ids(), names(expected))
  expect_true(all(names(expected) %in%
    sb_additional_survival_learner_ids()))
  for (learner in names(expected)) {
    expect_identical(
      sb_additional_survival_params(learner),
      expected[[learner]]
    )
    expect_identical(
      sb_additional_survival_dependency(learner),
      list(package = "survival", packages = "survival", available = TRUE)
    )
  }
})

test_that("fixed SurvReg risk and survival equal native distribution results", {
  data <- make_sb_fixture(n = 72L, seed = 8521L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:7, , drop = FALSE]
  horizon <- 7

  for (learner in sb_survreg_learner_ids()) {
    params <- sb_additional_survival_params(learner)
    adapter <- sb_additional_survival_learner_adapter(learner, params)
    model <- adapter$fit(x, data$time, data$event)
    observed_risk <- adapter$predict(model, assessment)
    observed_survival <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      assessment,
      prediction_horizon = horizon
    )
    native_location <- stats::predict(
      model$engine,
      newdata = assessment,
      type = "lp"
    )
    native_survival <- 1 - survival::psurvreg(
      horizon,
      mean = native_location,
      scale = model$engine$scale,
      distribution = params$distribution
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

    expect_equal(observed_risk, -unname(native_location), tolerance = 1e-14)
    expect_equal(
      observed_survival$probability,
      unname(native_survival),
      tolerance = 1e-14
    )
    expect_equal(replay$risk, observed_risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      observed_survival$probability,
      tolerance = 0
    )
    expect_false("engine_raw" %in% names(state))
    expect_identical(
      state$risk_reference,
      paste0(params$distribution, "_aft_native")
    )
    expect_identical(
      state$survival_probability$source,
      paste0("survreg_", params$distribution, "_distribution_direct")
    )
  }
})

test_that("fixed SurvReg preserves quoted predictor names and native predictions", {
  data <- make_sb_fixture(n = 72L, seed = 8521L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("NKX3-1", "2 gene", "if")
  assessment <- x[1:7, , drop = FALSE]

  for (learner in sb_survreg_learner_ids()) {
    params <- sb_additional_survival_params(learner)
    formula_data <- sb_additional_survival_formula_data(x, data$time, data$event)
    native <- survival::survreg(
      formula_data$formula, data = formula_data$data,
      dist = params$distribution, model = FALSE, x = FALSE, y = FALSE
    )
    adapter <- sb_additional_survival_learner_adapter(learner, params)
    model <- adapter$fit(x, data$time, data$event)
    location <- unname(stats::predict(native, newdata = assessment, type = "lp"))
    probability <- unname(1 - survival::psurvreg(
      7, mean = location, scale = native$scale,
      distribution = params$distribution
    ))
    state <- adapter$extract_state(model, list(
      x = x, time = data$time, event = data$event,
      prediction_horizon = 7, time_unit = "month"
    ))
    replay <- adapter$predict_state(unserialize(serialize(state, NULL)), assessment)

    expect_identical(model$features, names(x))
    expect_identical(state$features, names(x))
    expect_identical(stats::coef(model$engine), stats::coef(native))
    expect_identical(model$engine$scale, native$scale)
    expect_equal(adapter$predict(model, assessment), -location, tolerance = 1e-14)
    expect_equal(replay$risk, -location, tolerance = 1e-14)
    expect_equal(replay$survival_probability, probability, tolerance = 1e-14)
  }
})

test_that("fixed SurvReg still rejects misaligned and nonfinite coefficients", {
  data <- make_sb_fixture(n = 72L, seed = 8521L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  formula_data <- sb_additional_survival_formula_data(x, data$time, data$event)
  native <- survival::survreg(
    formula_data$formula, data = formula_data$data, dist = "weibull"
  )
  adapter <- sb_additional_survival_learner_adapter("survreg_weibull")
  bad <- native
  names(bad$coefficients)[2:3] <- rev(names(bad$coefficients)[2:3])
  local_mocked_bindings(sb_fit_survreg = function(...) bad)
  expect_error(adapter$fit(x, data$time, data$event),
    class = "stormbreaker_incomplete_survreg_fit")
  bad <- native
  bad$coefficients[[2L]] <- NA_real_
  expect_error(adapter$fit(x, data$time, data$event),
    class = "stormbreaker_incomplete_survreg_fit")
})
})

# ---- model-mtlr ----
local({
test_that("MTLR diagnostics preserve native time weights and exact input-unit potentials", {
  skip_if_not_installed("MTLR")
  d <- sb_example_survival(160L, seed = 8612L)
  features <- c("gene_c", "gene_a", "gene_b")
  d$gene_c <- d$gene_c * 2 + 4
  fit <- sb_fit_final(sb_surv_task(d, "patient_id", "time", "event", features, "month"),
    sb_pipeline(learner = "mtlr_survival"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  e <- fit$learner_state$engine
  view <- sb_mtlr_diagnostics(fit)
  expect_equal(view$native_weights, e$weight_matrix, tolerance = 0)
  expect_equal(view$table$weight, as.numeric(e$weight_matrix[, -1]), tolerance = 0)
  expect_equal(view$time_points, e$time_points, tolerance = 0)
  expect_equal(view$table$time, rep(e$time_points, length(features)), tolerance = 0)
  expect_identical(view$normalization$feature, features)
  expect_equal(view$normalization$center, unname(e$scale$center), tolerance = 0)
  expect_equal(view$normalization$sd, unname(e$scale$sd), tolerance = 0)
  expect_equal(view$requested_intervals, ceiling(sqrt(nrow(d))))
  expect_equal(view$effective_time_points, nrow(e$weight_matrix))
  expect_lt(view$effective_time_points, view$requested_intervals)
  x <- as.matrix(d[features])
  normalized <- scale(x, center = e$scale$center, scale = e$scale$sd)
  expect_equal(unname(cbind(1, x) %*% t(view$original_unit_weights)),
    unname(cbind(1, normalized) %*% t(e$weight_matrix)), tolerance = 1e-12)
  newdata <- d[features]
  names(newdata) <- fit$learner_state$internal_features
  native_predict <- get("predict.mtlr", asNamespace("MTLR"))
  raw_engine <- e
  raw_engine$weight_matrix <- view$original_unit_weights
  raw_engine$scale <- NULL
  expect_equal(native_predict(raw_engine, newdata, type = "survivalcurve", add_zero = FALSE),
    native_predict(e, newdata, type = "survivalcurve", add_zero = FALSE), tolerance = 1e-12)
  biases <- sb_mtlr_diagnostics(fit, "biases")
  expect_equal(biases$table$bias, as.numeric(e$weight_matrix[, 1]), tolerance = 0)
  compact <- fit
  compact$diagnostic_artifacts <- NULL
  expect_identical(sb_mtlr_diagnostics(compact), view)
  # Synthetic zero-SD state exercises the native ignored-column prediction branch.
  zero <- compact
  zero$learner_state$engine$scale$sd[1] <- 0
  zero_view <- sb_mtlr_diagnostics(zero)
  expect_equal(as.numeric(zero_view$original_unit_weights[, 2]), rep(0, length(e$time_points)))
  zero_raw <- zero$learner_state$engine
  zero_raw$weight_matrix <- zero_view$original_unit_weights
  zero_raw$scale <- NULL
  expect_equal(native_predict(zero_raw, newdata, type = "survivalcurve", add_zero = FALSE),
    native_predict(zero$learner_state$engine, newdata, type = "survivalcurve", add_zero = FALSE), tolerance = 1e-12)
  bundle <- sb_native_diagnostics(compact, complete = TRUE)
  expect_identical(bundle$status$status, rep("available", 2L))
  for (q in bundle$diagnostics) {
    p <- sb_plot_algorithm_diagnostic(q)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    expect_null(p$labels$subtitle)
    expect_null(p$labels$caption)
  }
  expect_identical(serialize(fit, NULL, version = 2L), before)
})
})

# ---- model-rule-aft ----
local({
# ---- algorithm-aft ----

test_that("AFT diagnostic time ratios reproduce native conditional quantile ratios", {
  skip_if_not_installed("survival")
  d <- sb_example_survival(160L, seed = 8612L)
  features <- c("gene_c", "gene_a", "gene_b")
  d$gene_c <- d$gene_c * 2 + 4
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  for (id in sb_survreg_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    engine <- fit$diagnostic_artifacts$learner$engine
    view <- sb_aft_diagnostics(fit)
    expect_identical(view$table$feature, features)
    expect_equal(view$table$coefficient, unname(stats::coef(engine)[-1]), tolerance = 0)
    expect_equal(view$intercept, unname(stats::coef(engine)[1]), tolerance = 0)
    expect_equal(view$distribution_scale, unname(engine$scale), tolerance = 0)
    expect_identical(view$distribution, engine$dist)
    rows <- d[1:5, features, drop = FALSE]
    native_predict <- get("predict.survreg", asNamespace("survival"))
    base <- native_predict(engine, newdata = rows, type = "quantile", p = c(.2, .5, .8))
    for (j in seq_along(features)) {
      shifted <- rows
      shifted[[j]] <- shifted[[j]] + 1
      quantiles <- native_predict(engine, newdata = shifted, type = "quantile", p = c(.2, .5, .8))
      expect_equal(unname(quantiles / base), matrix(view$table$time_ratio[j], 5, 3), tolerance = 1e-12)
    }
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_identical(sb_aft_diagnostics(compact), view)
    bundle <- sb_native_diagnostics(compact, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 2L))
    for (diagnostic in bundle$diagnostics) {
      expect_identical(diagnostic$table, view$table)
      p <- sb_plot_algorithm_diagnostic(diagnostic)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_null(p$labels$subtitle)
      expect_null(p$labels$caption)
    }
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})

# ---- algorithm-pre-rules ----

test_that("PRE diagnostics preserve native terms, predictions and coefficient CV", {
  skip_if_not_installed("pre")
  skip_if_not_installed("mboost")
  for (kind in c("binary", "survival")) {
    d <- make_sb_binary_fixture(n = 150L, seed = 108L)
    d$time <- exp(2 + .8 * d$gene_a - .4 * d$gene_b)
    d$event <- rep(c(1L, 1L, 0L), length.out = nrow(d))
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (kind == "binary") sb_binary_task(d, "patient_id", "outcome", features) else
      sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = paste0("pre_", kind, "_rule_ensemble")))
    expect_identical(fit$status, "success")
    expect_null(fit$diagnostic_artifacts)
    original <- fit
    q <- sb_rule_ensemble(fit)
    engine <- unserialize(fit$learner_state$engine_raw)
    native <- stats::coef(engine, penalty.par.val = fit$learner_state$params$penalty_rule)
    rownames(native) <- NULL
    expect_identical(q$table[names(native)], native)
    expect_equal(q$table$selected, native$coefficient != 0)
    expect_equal(q$table$term_type[q$table$rule %in% engine$rules$rule],
      rep("rule", sum(q$table$rule %in% engine$rules$rule)))
    # These expressions come only from this synthetic native fit, not user text.
    reconstructed <- Reduce(`+`, lapply(seq_len(nrow(q$table)), function(i) {
      value <- if (q$table$term_type[i] == "linear") {
        pmin(pmax(d[[q$table$rule[i]]], q$table$winsor_lower[i]), q$table$winsor_upper[i])
      } else eval(parse(text = q$table$description[i]), envir = d)
      q$table$coefficient[i] * value
    }))
    native_link <- as.numeric(stats::predict(engine, newdata = d[features], type = "link"))
    expect_equal(as.numeric(reconstructed), native_link, tolerance = 1e-10)
    cv <- sb_rule_ensemble(fit, "rule_cv_curve")
    expect_equal(cv$table$lambda, engine$glmnet.fit$lambda)
    expect_equal(cv$table$mean_value, unname(engine$glmnet.fit$cvm))
    expect_equal(cv$table$standard_error, unname(engine$glmnet.fit$cvsd))
    expect_equal(cv$table$nonzero, unname(engine$glmnet.fit$nzero))
    expect_equal(cv$table$lambda[cv$table$selected], engine$glmnet.fit$lambda.1se)
    expect_identical(fit, original)
    expect_equal(length(grepRaw("B001", serialize(q, NULL), fixed = TRUE)), 0)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      for (z in list(q, cv)) {
        p <- sb_plot_algorithm_diagnostic(z)
        expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      }
      p <- sb_plot_algorithm_diagnostic(q)
      expect_equal(nrow(p$data), sum(q$table$selected & q$table$term_type != "intercept"))
      empty <- q
      empty$table$coefficient[] <- 0
      empty$table$selected[] <- FALSE
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(empty)), "ggplot_built")
    }
    bad <- fit
    engine$glmnet.fit$cvsd[1] <- -1
    bad$learner_state$engine_raw <- serialize(engine, NULL)
    err <- tryCatch(sb_rule_ensemble(bad, "rule_cv_curve"), error = identity)
    expect_s3_class(err, "stormbreaker_invalid_diagnostic")
    err <- tryCatch(sb_rule_ensemble(fit, "surrogate_tree"), error = identity)
    expect_s3_class(err, "stormbreaker_unsupported_diagnostic")
  }
})

# ---- rule-aft-adapters ----

make_rule_aft_fixture <- function(n = 72L, p = 12L, seed = 9101L) {
  withr::with_seed(seed, {
    expression <- matrix(stats::rnorm(n * p), nrow = n)
    colnames(expression) <- paste0("gene_", seq_len(p))
    nonlinear <- 1.3 * expression[, 1L] * expression[, 2L] +
      0.8 * expression[, 3L] - 0.5 * expression[, 4L]^2
    outcome <- as.integer(nonlinear > stats::median(nonlinear))
    log_time <- 2 + 0.6 * expression[, 1L] -
      0.5 * expression[, 2L] +
      0.7 * (expression[, 3L] > 0) + stats::rnorm(n, sd = 0.35)
    data.frame(
      patient_id = paste0("RULE_AFT_PRIVATE_", seq_len(n)),
      time = exp(log_time),
      event = rep(c(1L, 1L, 0L, 1L, 0L), length.out = n),
      outcome = outcome,
      as.data.frame(expression, check.names = FALSE),
      check.names = FALSE
    )
  })
}

rule_aft_features <- function(data) {
  setdiff(names(data), c("patient_id", "time", "event", "outcome"))
}

rule_aft_token_params <- function(learner) {
  params <- sb_rule_aft_params(learner)
  if (identical(learner, "bart_aft_survival")) {
    params$num_trees <- 15L
    params$posterior_draws <- 40L
    params$burn_in <- 20L
    params$time_grid_size <- 20L
  } else if (learner %in% sb_pre_rule_ensemble_learner_ids()) {
    params$num_trees <- 80L
    params$internal_folds <- 3L
  } else {
    params$num_lambda <- 8L
    params$internal_folds <- 3L
  }
  params
}

test_that("four rule/AFT learners expose only meaningful routes", {
  ids <- sb_rule_aft_learner_ids()
  expect_identical(ids, c(
    "bart_aft_survival",
    "pre_binary_rule_ensemble",
    "pre_survival_rule_ensemble",
    "penaft_gehan_elastic_net"
  ))
  registry <- sb_pipeline_registry()
  rows <- registry[registry$learner %in% ids, , drop = FALSE]
  expect_setequal(
    rows$selector[rows$learner == "bart_aft_survival"],
    expected_fixed_learner_selectors("survival")
  )
  for (learner in setdiff(ids, "bart_aft_survival")) {
    expect_setequal(
      rows$selector[rows$learner == learner],
      c("none", "fixed_features")
    )
  }
  pre_rows <- rows$learner %in% sb_pre_rule_ensemble_learner_ids()
  expect_identical(
    unique(rows$tuning_semantics[pre_rows]),
    "analysis_fold_native_coefficient_cv_after_rule_induction"
  )
  expect_true(all(grepl(
    "complete analysis fold",
    rows$parameter_semantics[pre_rows],
    fixed = TRUE
  )))
  expect_identical(
    unique(rows$optional_dependency[
      rows$learner == "pre_survival_rule_ensemble"
    ]),
    "pre + mboost"
  )
  expect_error(
    sb_pipeline(
      selector = "univ_cox_top_k",
      learner = "pre_survival_rule_ensemble",
      selector_n_features = 5L
    ),
    class = "stormbreaker_unadmitted_pipeline"
  )
})

test_that("native rule/AFT fits are high-dimensional, portable, and local-RNG", {
  skip_if_not_installed("BART", minimum_version = "2.9.0")
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("pre", minimum_version = "1.1.0")
  skip_if_not_installed("penAFT", minimum_version = "0.3.2")
  data <- make_rule_aft_fixture(n = 54L, p = 70L, seed = 9102L)
  features <- rule_aft_features(data)
  x <- data[features]
  rownames(x) <- data$patient_id
  assessment <- x[c(2L, 19L, 47L), , drop = FALSE]
  horizon <- as.numeric(stats::median(data$time))

  for (learner in sb_rule_aft_learner_ids()) {
    adapter <- sb_rule_aft_learner_adapter(
      learner,
      rule_aft_token_params(learner)
    )
    set.seed(9200L)
    caller_rng <- .Random.seed
    model <- if (identical(learner, "pre_binary_rule_ensemble")) {
      adapter$fit(x, data$outcome)
    } else {
      adapter$fit(x, data$time, data$event)
    }
    expect_identical(.Random.seed, caller_rng)
    prediction <- adapter$predict(model, assessment)
    expect_length(prediction, nrow(assessment))
    expect_equal(sum(is.finite(prediction)), nrow(assessment))
    if (identical(learner, "pre_binary_rule_ensemble")) {
      native <- stats::predict(
        model$engine,
        newdata = assessment,
        type = "response",
        penalty.par.val = model$params$penalty_rule
      )
      expect_equal(prediction, unname(as.numeric(native)), tolerance = 0)
    }
    if (identical(learner, "penaft_gehan_elastic_net")) {
      manual <- -unname(
        model$intercept +
          as.vector(as.matrix(assessment) %*% model$coefficient)
      )
      expect_equal(prediction, manual, tolerance = 0)
    }

    context <- if (learner %in% sb_rule_aft_survival_learner_ids()) {
      list(
        x = x,
        time = data$time,
        event = data$event,
        prediction_horizon = horizon,
        time_unit = "month"
      )
    } else {
      list()
    }
    state <- unserialize(serialize(
      adapter$extract_state(model, context),
      NULL
    ))
    replay <- adapter$predict_state(state, assessment)
    replay_value <- if (identical(learner, "pre_binary_rule_ensemble")) {
      replay$probability
    } else {
      replay$risk
    }
    expect_equal(replay_value, prediction, tolerance = 1e-10)
    if (learner %in% sb_rule_aft_survival_learner_ids()) {
      expect_equal(
        sum(replay$survival_probability >= 0 &
          replay$survival_probability <= 1),
        nrow(assessment)
      )
    }
    expect_length(
      grepRaw(
        "RULE_AFT_PRIVATE_",
        serialize(state, NULL),
        fixed = TRUE,
        all = TRUE
      ),
      0L
    )
  }
})

test_that("AFT-BART and penAFT horizon formulas match their fitted states", {
  skip_if_not_installed("BART", minimum_version = "2.9.0")
  skip_if_not_installed("penAFT", minimum_version = "0.3.2")
  data <- make_rule_aft_fixture(n = 60L, p = 18L, seed = 9103L)
  features <- rule_aft_features(data)
  x <- data[features]
  assessment <- x[c(3L, 28L, 55L), , drop = FALSE]
  horizon <- as.numeric(stats::median(data$time))

  bart <- sb_rule_aft_learner_adapter(
    "bart_aft_survival",
    rule_aft_token_params("bart_aft_survival")
  )
  bart_model <- bart$fit(x, data$time, data$event)
  draws <- sb_bart_aft_log_time_draws(bart_model, assessment)
  manual_bart <- colMeans(stats::pnorm(
    (log(horizon) - draws) / bart_model$engine$sigma,
    lower.tail = FALSE
  ))
  bart_probability <- bart$predict_survival(
    bart_model,
    x,
    data$time,
    data$event,
    assessment,
    horizon
  )
  expect_equal(
    bart_probability$probability,
    unname(manual_bart),
    tolerance = 0
  )

  penaft <- sb_rule_aft_learner_adapter(
    "penaft_gehan_elastic_net",
    rule_aft_token_params("penaft_gehan_elastic_net")
  )
  penaft_model <- penaft$fit(x, data$time, data$event)
  link <- sb_penaft_log_time_link(penaft_model, x)
  curve <- sb_fit_mboost_gehan_residual_survival(
    data$time,
    data$event,
    link
  )
  manual_penaft <- sb_predict_mboost_gehan_residual_survival(
    curve,
    sb_penaft_log_time_link(penaft_model, assessment),
    horizon
  )
  penaft_probability <- penaft$predict_survival(
    penaft_model,
    x,
    data$time,
    data$event,
    assessment,
    horizon
  )
  expect_equal(
    penaft_probability$probability,
    manual_penaft,
    tolerance = 0
  )
})

test_that("rule/AFT learners run through CV and outcome-blind final prediction", {
  skip_on_cran()
  skip_if_not_installed("BART", minimum_version = "2.9.0")
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("pre", minimum_version = "1.1.0")
  skip_if_not_installed("penAFT", minimum_version = "0.3.2")
  data <- make_rule_aft_fixture(n = 75L, p = 6L, seed = 9104L)
  features <- rule_aft_features(data)
  horizon <- as.numeric(stats::median(data$time))

  for (learner in sb_rule_aft_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "fixed_features",
      learner = learner,
      selector_features = features,
      feature_set_id = "rule_aft_test"
    )
    task <- if (identical(learner, "pre_binary_rule_ensemble")) {
      sb_binary_task(
        data,
        "patient_id",
        "outcome",
        features
      )
    } else {
      sb_surv_task(
        data,
        "patient_id",
        "time",
        "event",
        features,
        "month",
        prediction_horizon = horizon
      )
    }
    cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 9301L)
    expect_identical(unique(cv$runs$status), "success")
    expect_identical(unique(cv$predictions$status), "success")
    if (identical(learner, "pre_survival_rule_ensemble")) {
      expect_true(all(c("pre_version", "mboost_version") %in%
        names(cv$provenance)))
      expect_false(anyNA(unlist(cv$provenance[
        c("pre_version", "mboost_version")
      ])))
    }

    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    if (identical(learner, "pre_survival_rule_ensemble")) {
      expect_true(all(c("pre_version", "mboost_version") %in%
        names(fit$provenance)))
      expect_false(anyNA(unlist(fit$provenance[
        c("pre_version", "mboost_version")
      ])))
    }
    new_data <- data[1:5, , drop = FALSE]
    altered <- new_data
    altered$outcome <- 1L - altered$outcome
    altered$time <- altered$time + 1000
    altered$event <- 1L - altered$event
    prediction <- sb_predict_final(fit, new_data, "patient_id")
    changed <- sb_predict_final(fit, altered, "patient_id")
    expect_identical(prediction$status, "success")
    value <- if (identical(learner, "pre_binary_rule_ensemble")) {
      "probability"
    } else {
      "risk"
    }
    expect_equal(
      prediction$predictions[[value]],
      changed$predictions[[value]],
      tolerance = 0
    )
    if (learner %in% sb_rule_aft_survival_learner_ids()) {
      expect_equal(
        prediction$predictions$survival_probability,
        changed$predictions$survival_probability,
        tolerance = 0
      )
    }
  }
})

test_that("rule/AFT dependencies include the survival PRE engine", {
  dependency <- sb_rule_aft_dependency("pre_survival_rule_ensemble")
  expect_identical(dependency$packages, c("pre", "mboost"))

  local_mocked_bindings(
    sb_rule_aft_dependency = function(learner) {
      list(
        package = "pre",
        packages = c("pre", "mboost"),
        installed = c(pre = TRUE, mboost = FALSE),
        available = FALSE
      )
    }
  )
  expect_error(
    sb_require_rule_aft_dependency("pre_survival_rule_ensemble"),
    "mboost",
    class = "stormbreaker_missing_rule_aft_dependency"
  )
})
})

# ---- model-timereg ----
local({
test_that("timereg views retain native cumulative processes and coefficient meanings", {
  skip_if_not_installed("timereg")
  d <- sb_example_survival(160L, seed = 8612L)
  names(d)[names(d) == "gene_c"] <- "gene c"
  names(d)[names(d) == "gene_a"] <- "gene a"
  features <- c("gene c", "gene a", "gene_b")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  for (id in sb_timereg_survival_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    aalen <- id == "aalen_additive_survival"
    native <- withr::with_seed(1L, if (aalen) timereg::aalen(
      survival::Surv(time, event) ~ `gene c` + `gene a` + gene_b, data = d,
      robust = 0, residuals = 0, n.sim = 0, weighted.test = 0, covariance = 0, resample.iid = 0, silent = 1) else
      timereg::prop.odds(timereg::Event(time, event) ~ `gene c` + `gene a` + gene_b, data = d,
        Nit = 20, detail = 0, n.sim = 0, weighted.test = 0, profile = 1, sym = 0, baselinevar = 1))
    effects <- sb_timereg_diagnostics(fit)
    baseline <- sb_timereg_diagnostics(fit, "baseline")
    expect_equal(baseline$table$time, unname(native$cum[, 1]), tolerance = 0)
    expect_equal(baseline$table$value, unname(native$cum[, 2]), tolerance = 0)
    expect_identical(effects$features, features)
    if (aalen) {
      expect_equal(effects$table$value, as.numeric(native$cum[, -(1:2)]), tolerance = 0)
      expect_equal(effects$table$time, rep(unname(native$cum[, 1]), length(features)), tolerance = 0)
      expect_equal(nrow(effects$table), nrow(native$cum) * length(features))
      expect_gt(sum(effects$table$value < 0), 0)
    } else expect_equal(effects$table$value, as.numeric(native$gamma[, 1]), tolerance = 0)
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_identical(sb_timereg_diagnostics(compact), effects)
    bundle <- sb_native_diagnostics(compact, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 2L))
    for (view in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(view)
      built <- ggplot2::ggplot_build(p)
      expect_s3_class(built, "ggplot_built")
      expect_equal(if (view$type == "timereg_odds_effects") built$data[[2]]$x else built$data[[2]]$y,
        view$table$value, tolerance = 0)
      expect_null(p$labels$subtitle)
      expect_null(p$labels$caption)
    }
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
})

# ---- survival-numeric-regressions ----
local({
make_survival_numeric_fixture <- function() {
  withr::with_seed(20260930L, {
    n <- 600L
    x <- as.data.frame(matrix(stats::rnorm(n * 24L), nrow = n))
    names(x) <- sprintf("gene%02d", seq_len(ncol(x)))
    beta <- c(.65, -.60, .55, -.50, .45, -.40, .30, -.25, rep(0, 16L))
    eta <- drop(as.matrix(x) %*% beta)
    event_time <- 365 * (-log(stats::runif(n)) / exp(eta))^(1 / 1.3)
    censor_time <- stats::rexp(n, rate = 1 / 900)
    list(x = x, time = pmin(event_time, censor_time),
         event = as.integer(event_time <= censor_time))
  })
}

test_that("gamboost fits more than 500 rows without sparse QR warnings", {
  skip_if_not_installed("mboost")
  d <- make_survival_numeric_fixture()
  old <- getOption("mboost_useMatrix")
  ids <- c("gamboost_cox", "gamboost_weibull_aft",
    "gamboost_loglogistic_aft", "gamboost_lognormal_aft",
    "gamboost_gehan_aft", "gamboost_unoc")
  for (id in ids) {
    adapter <- sb_resolve_pipeline_runtime(sb_pipeline("none", id))$learner
    model <- adapter$fit(d$x, d$time, d$event)
    risk <- adapter$predict(model, d$x[seq_len(12L), , drop = FALSE])
    expect_identical(all(is.finite(risk)), TRUE)
    expect_gt(length(unique(risk)), 1L)
    expect_equal(mboost::mstop(model$engine), 100)
    expect_identical(getOption("mboost_useMatrix"), old)
  }
})

test_that("Gompertz initialization is stable across time units", {
  skip_if_not_installed("flexsurv")
  d <- make_survival_numeric_fixture()
  adapter <- sb_flexsurv_parametric_learner_adapter("flexsurv_gompertz")
  days <- adapter$fit(d$x, d$time, d$event)
  years <- adapter$fit(d$x, d$time / 365, d$event)
  assessment <- d$x[seq_len(12L), , drop = FALSE]
  p_days <- sb_flexsurv_parametric_direct_probability(days, assessment, 365)
  p_years <- sb_flexsurv_parametric_direct_probability(years, assessment, 1)
  expect_equal(p_days, p_years, tolerance = 1e-5)
  expect_identical(days$engine$opt$convergence, 0L)
  expect_identical(all(is.finite(days$engine$cov)), TRUE)
  native_x <- sb_flexsurv_numerical_predictors(assessment, days$numerical_fitting)
  names(native_x) <- days$engine_features
  native <- as.numeric(stats::predict(days$engine, newdata = native_x,
    type = "survival", times = 365, conf.int = FALSE)$.pred_survival)
  expect_equal(p_days, native, tolerance = 1e-12)
  expect_equal(adapter$predict(years, assessment) - adapter$predict(days, assessment),
    rep(log(365), nrow(assessment)), tolerance = 1e-5)
})

test_that("LPC explains the native fifty-feature automatic-search limitation", {
  skip_if_not_installed("lpc")
  d <- make_survival_numeric_fixture()
  endpoint <- list(task_type = "survival", time = d$time, event = d$event)
  error <- tryCatch(sb_fit_lpc_selector(d$x, endpoint, "lpc_cox_nlogn"), error = identity)
  expect_s3_class(error, "stormbreaker_invalid_lpc_selector")
  expect_s3_class(error, "stormbreaker_lpc_dimension_unsupported")
  expect_match(conditionMessage(error), "at least 50 predictors", fixed = TRUE)
})

test_that("PLS-Cox resolves late-component numerical warnings without lowering the cap", {
  skip_if_not_installed("plsRcox")
  d <- make_survival_numeric_fixture()
  panels <- list(c(1:9, 11, 13, 16, 20, 21, 24), c(1:9, 16))
  for (panel in panels) {
    x <- d$x[panel]
    adapter <- sb_plsrcox_cv_learner_adapter(inner_folds = 5L, seed = 1L,
      maximum_components = 10L)
    model <- adapter$fit(x, d$time, d$event)
    expect_length(model$audit$component_search_warnings, 0L)
    expect_identical(model$audit$maximum_components, 10L)
    expect_equal(model$audit$cox_control$eps, 1e-11)
    expect_identical(model$audit$cox_control$toler.chol,
      survival::coxph.control()$toler.chol)
    expected <- as.numeric(stats::predict(model$engine,
      newdata = x[1:12, , drop = FALSE], type = "lp"))
    actual <- adapter$predict(model, x[1:12, , drop = FALSE])
    expect_equal(actual - mean(actual), expected - mean(expected), tolerance = 1e-8)
    expect_gt(length(unique(actual)), 1L)
  }
})

test_that("PLS-Cox still rejects genuine selected-component search warnings", {
  skip_if_not_installed("plsRcox")
  d <- make_survival_numeric_fixture()
  testthat::local_mocked_bindings(cv.plsRcox = function(...) {
    warning("Loglik converged before variable  1 ; coefficient may be infinite. ")
    list(lambda.min5 = 1L)
  }, .package = "plsRcox")
  error <- tryCatch(sb_plsrcox_cv_learner_adapter()$fit(d$x, d$time, d$event),
    error = identity)
  expect_s3_class(error, "stormbreaker_plsrcox_cv_warning")
  expect_s3_class(error$parent, "warning")
})
})

# ---- timereg-survival-fixed ----
local({
make_timereg_survival_fixture <- function(n = 90L, seed = 20260831L) {
  withr::with_seed(seed, {
    gene_a <- stats::runif(n)
    gene_b <- stats::runif(n)
    gene_c <- stats::runif(n)
    event_time <- stats::rexp(
      n,
      rate = 0.08 + 0.08 * gene_a + 0.05 * gene_b + 0.03 * gene_c
    )
    censor_time <- stats::rexp(n, rate = 0.025)
    data.frame(
      patient_id = sprintf("TR%03d", seq_len(n)),
      time = pmax(pmin(event_time, censor_time), 0.01),
      event = as.integer(event_time <= censor_time),
      gene_a,
      gene_b,
      gene_c,
      check.names = FALSE
    )
  })
}

test_that("timereg contributes two distinct fixed survival learners", {
  expect_identical(
    sb_timereg_survival_learner_ids(),
    c("aalen_additive_survival", "prop_odds_survival")
  )

  registry <- sb_pipeline_registry()
  rows <- registry[
    registry$learner %in% sb_timereg_survival_learner_ids(), ,
    drop = FALSE
  ]
  expect_identical(nrow(rows), 2L * length(expected_fixed_learner_selectors("survival")))
  expect_setequal(
    unique(rows$selector),
    expected_fixed_learner_selectors("survival")
  )
  expect_true(all(grepl("timereg", rows$optional_dependency, fixed = TRUE)))

  aalen <- sb_pipeline(learner = "aalen_additive_survival")
  prop_odds <- sb_pipeline(learner = "prop_odds_survival")
  expect_identical(
    aalen$risk_formula,
    "cumulative_hazard_at_last_analysis_event"
  )
  expect_identical(aalen$learner_params$risk_time_rule, "last_analysis_event")
  expect_identical(prop_odds$risk_formula, "x_beta")
  expect_identical(prop_odds$risk_reference, "zero")
})

test_that("timereg portable predictions reproduce native model formulas", {
  skip_if_not_installed("timereg", minimum_version = "2.0.7")
  d <- make_timereg_survival_fixture(n = 120L, seed = 9431L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:9, , drop = FALSE]

  set.seed(9432L)
  caller_seed <- .Random.seed
  aalen_adapter <- sb_timereg_survival_learner_adapter(
    "aalen_additive_survival"
  )
  aalen <- aalen_adapter$fit(x, d$time, d$event)
  expect_identical(.Random.seed, caller_seed)

  aalen_risk <- aalen_adapter$predict(aalen, assessment)
  expected_aalen_risk <- unname(
    aalen$cumulative_intercept[[length(aalen$time)]] +
      as.vector(
        as.matrix(assessment) %*%
          aalen$cumulative_coefficient[nrow(aalen$cumulative_coefficient), ]
      )
  )
  expect_equal(aalen_risk, expected_aalen_risk, tolerance = 0)
  native_aalen <- withr::with_seed(
    1L,
    timereg::aalen(
      survival::Surv(time, event) ~ gene_a + gene_b + gene_c,
      data = d,
      robust = 0L,
      residuals = 0L,
      n.sim = 0L,
      weighted.test = 0L,
      covariance = 0L,
      resample.iid = 0L,
      silent = 1L
    )
  )
  native_aalen_survival <- stats::predict(
    native_aalen,
    X = cbind(1, as.matrix(assessment)),
    times = aalen$risk_reference_time,
    se = FALSE,
    uniform = FALSE
  )$S0
  expect_equal(
    aalen_risk,
    -log(unname(as.numeric(native_aalen_survival))),
    tolerance = 1e-12
  )
  native_aalen_horizon <- stats::predict(
    native_aalen,
    X = cbind(1, as.matrix(assessment)),
    times = 5,
    se = FALSE,
    uniform = FALSE
  )$S0
  portable_aalen_horizon <- aalen_adapter$predict_survival(
    aalen,
    x,
    d$time,
    d$event,
    assessment,
    5
  )$probability
  expect_equal(
    portable_aalen_horizon,
    unname(as.numeric(native_aalen_horizon)),
    tolerance = 1e-12
  )

  prop_adapter <- sb_timereg_survival_learner_adapter(
    "prop_odds_survival"
  )
  prop_odds <- prop_adapter$fit(x, d$time, d$event)
  prop_risk <- prop_adapter$predict(prop_odds, assessment)
  expect_equal(
    prop_risk,
    unname(as.vector(as.matrix(assessment) %*% prop_odds$coefficient)),
    tolerance = 0
  )

  horizon <- 5
  native_prop <- withr::with_seed(
    1L,
    timereg::prop.odds(
      timereg::Event(time, event) ~ gene_a + gene_b + gene_c,
      data = d,
      Nit = 20L,
      detail = 0L,
      n.sim = 0L,
      weighted.test = 0L,
      profile = 1L,
      sym = 0L,
      baselinevar = 1L
    )
  )
  native_survival <- stats::predict(
    native_prop,
    X = matrix(1, nrow(assessment), 1L),
    Z = as.matrix(assessment),
    times = horizon,
    se = FALSE,
    uniform = FALSE
  )$S0
  portable_survival <- prop_adapter$predict_survival(
    prop_odds,
    x,
    d$time,
    d$event,
    assessment,
    horizon
  )$probability
  expect_equal(
    portable_survival,
    unname(as.numeric(native_survival)),
    tolerance = 1e-12
  )
})

test_that("Aalen native survival is never silently clipped", {
  model <- structure(
    list(
      learner = "aalen_additive_survival",
      features = "gene_a",
      params = sb_timereg_survival_params(
        "aalen_additive_survival"
      ),
      time = c(0, 1),
      cumulative_intercept = c(0, -0.2),
      cumulative_coefficient = matrix(c(0, 0), ncol = 1L),
      risk_reference_time = 1
    ),
    class = "stormbreaker_timereg_survival_fixed"
  )

  expect_error(
    sb_timereg_direct_survival_probability(
      model,
      data.frame(gene_a = 0),
      horizon = 1
    ),
    class = "stormbreaker_invalid_timereg_survival_prediction"
  )
})

test_that("Aalen horizon state validation is invariant to patient row order", {
  model <- structure(
    list(
      learner = "aalen_additive_survival",
      features = "gene_a",
      params = sb_timereg_survival_params(
        "aalen_additive_survival"
      ),
      time = c(0, 1),
      cumulative_intercept = c(0, 0.2),
      cumulative_coefficient = matrix(c(0, 1), ncol = 1L),
      risk_reference_time = 1,
      baseline_cumulative_odds = NULL,
      coefficient = NULL
    ),
    class = "stormbreaker_timereg_survival_fixed"
  )
  context <- list(
    x = data.frame(gene_a = c(0, -1)),
    time = c(1, 2),
    event = c(1L, 0L),
    prediction_horizon = 1,
    time_unit = "month"
  )

  expect_error(
    sb_extract_timereg_survival_state(model, context),
    class = "stormbreaker_invalid_timereg_survival_prediction"
  )
  reversed <- c(2L, 1L)
  context$x <- context$x[reversed, , drop = FALSE]
  context$time <- context$time[reversed]
  context$event <- context$event[reversed]
  expect_error(
    sb_extract_timereg_survival_state(model, context),
    class = "stormbreaker_invalid_timereg_survival_prediction"
  )
})

test_that("timereg learners execute CV and portable final prediction", {
  skip_if_not_installed("timereg", minimum_version = "2.0.7")
  d <- make_timereg_survival_fixture(n = 90L, seed = 9433L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  risk_task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month"
  )

  for (learner in sb_timereg_survival_learner_ids()) {
    pipeline <- sb_pipeline(learner = learner)
    cv <- sb_run_cv(risk_task, pipeline, folds = 3L, seed = 9434L)
    fit <- sb_fit_final(risk_task, pipeline)
    state <- unserialize(serialize(fit$learner_state, NULL))
    fit$learner_state <- state
    new_data <- d[1:6, c("patient_id", "gene_c", "gene_a", "gene_b")]
    with_endpoints <- new_data
    with_endpoints$time <- -999
    with_endpoints$event <- 1L
    prediction <- sb_predict_final(fit, new_data, "patient_id")
    endpoint_prediction <- sb_predict_final(
      fit,
      with_endpoints,
      "patient_id"
    )

    expect_identical(cv$status, "success", info = learner)
    expect_identical(cv$complete_oof, TRUE, info = learner)
    expect_identical(fit$status, "success", info = learner)
    expect_identical(prediction$status, "success", info = learner)
    expect_equal(
      prediction$predictions$risk,
      endpoint_prediction$predictions$risk,
      tolerance = 0,
      info = learner
    )
    expect_identical(
      fit$provenance$timereg_version,
      as.character(utils::packageVersion("timereg")),
      info = learner
    )
    expect_false("engine_raw" %in% names(state), info = learner)
    expect_false(any(c("patient_ids", "time_outcome", "event") %in%
      names(state)), info = learner)
  }

  horizon_task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 5
  )
  fit <- sb_fit_final(
    horizon_task,
    sb_pipeline(learner = "prop_odds_survival")
  )
  prediction <- sb_predict_final(
    fit,
    d[1:6, c("patient_id", predictors)],
    "patient_id"
  )
  expect_identical(fit$status, "success")
  expect_identical(prediction$status, "success")
  expect_true(all(prediction$predictions$survival_probability >= 0 &
    prediction$predictions$survival_probability <= 1))
})
})
