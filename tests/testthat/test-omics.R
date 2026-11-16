# ---- model-omics ----
local({
# ---- cancerclass-binary-adapter ----

test_that("cancerclass reproduces native z scores and fold-local calibration", {
  skip_if_not_installed("cancerclass", minimum_version = "1.52.0")
  skip_if_not_installed("Biobase", minimum_version = "2.68.0")

  d <- make_sb_binary_fixture(n = 60L, seed = 8901L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  outcome <- d$outcome
  sample_names <- paste0("sample_", seq_len(nrow(x)))
  make_eset <- function(rows) {
    expression <- t(data.matrix(x[rows, , drop = FALSE]))
    colnames(expression) <- sample_names[rows]
    phenotype <- data.frame(
      class = factor(outcome[rows], levels = c(0L, 1L)),
      row.names = sample_names[rows]
    )
    Biobase::ExpressionSet(
      assayData = expression,
      phenoData = Biobase::AnnotatedDataFrame(phenotype)
    )
  }

  adapter <- sb_cancerclass_binary_learner_adapter()
  search_before_fit <- search()
  model <- adapter$fit(x, outcome)
  expect_identical(search(), search_before_fit)
  assessment_rows <- c(3L, 14L, 27L, 45L, 58L)
  withr::local_package("generics", quietly = TRUE)
  withr::local_package("BiocGenerics", quietly = TRUE)
  withr::local_package("Biobase", quietly = TRUE)
  native <- suppressMessages(cancerclass::fit(
    make_eset(seq_len(nrow(x))),
    class = "class",
    method = "welch.test",
    hparam = 0.75
  ))
  native_prediction <- suppressMessages(cancerclass::predict(
    native,
    make_eset(assessment_rows),
    positive = "1",
    class = "class",
    ngenes = ncol(x),
    dist = "cor"
  ))
  expected_score <- as.numeric(
    methods::slot(native_prediction, "prediction")[, "z"]
  )
  raw_score <- sb_cancerclass_raw_score(
    model,
    x[assessment_rows, , drop = FALSE]
  )
  probability <- adapter$predict(
    model,
    x[assessment_rows, , drop = FALSE]
  )

  expect_equal(model$predictor, methods::slot(native, "predictor"),
    tolerance = 0
  )
  expect_equal(raw_score, expected_score, tolerance = 1e-14)
  expect_equal(
    probability,
    stats::plogis(
      model$calibration_coefficients[[1L]] +
        model$calibration_coefficients[[2L]] * raw_score
    ),
    tolerance = 0
  )

  expected_fold_id <- sb_make_stratified_fold_id(
    outcome,
    folds = 5L,
    seed = 1L,
    stratum_label = "class"
  )
  expected_oof <- numeric(nrow(x))
  for (fold in seq_len(5L)) {
    analysis_rows <- which(expected_fold_id != fold)
    held_out_rows <- which(expected_fold_id == fold)
    inner_fit <- suppressMessages(cancerclass::fit(
      make_eset(analysis_rows),
      class = "class",
      method = "welch.test",
      hparam = 0.75
    ))
    predictor <- methods::slot(inner_fit, "predictor")
    held_out <- t(data.matrix(x[held_out_rows, , drop = FALSE]))
    held_out <- held_out[rownames(predictor), , drop = FALSE]
    distance_zero <- cancerclass::get.d(
      held_out,
      predictor[, "0"],
      method = "cor"
    )
    distance_one <- cancerclass::get.d(
      held_out,
      predictor[, "1"],
      method = "cor"
    )
    class_distance <- cancerclass::get.d(
      predictor[, "0"],
      predictor[, "1"],
      method = "cor"
    )
    expected_oof[held_out_rows] <-
      (distance_zero - distance_one) / class_distance
  }

  expect_identical(model$calibration_fold_id, expected_fold_id)
  expect_equal(model$calibration_oof_score, expected_oof, tolerance = 0)
  state <- unserialize(serialize(adapter$extract_state(model), NULL))
  replay <- adapter$predict_state(
    state,
    x[assessment_rows, , drop = FALSE]
  )
  expect_equal(replay$probability, probability, tolerance = 0)
  expect_identical(state$training_data_retention, "none")
  expect_identical(
    intersect(
      c("outcome", "calibration_fold_id", "calibration_oof_score"),
      names(state)
    ),
    character()
  )
})

test_that("cancerclass exposes only leakage-safe routes and runs binary CV", {
  skip_if_not_installed("cancerclass", minimum_version = "1.52.0")
  skip_if_not_installed("Biobase", minimum_version = "2.68.0")

  registry <- sb_pipeline_registry()
  rows <- registry[
    registry$learner == "cancerclass_binary",
    ,
    drop = FALSE
  ]
  expect_setequal(rows$selector, c("none", "fixed_features"))
  expect_identical(unique(rows$task_type), "binary")
  expect_identical(unique(rows$optional_dependency), "cancerclass + Biobase")
  unsafe <- capture_sb_error(sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "cancerclass_binary",
    selector_n_features = 2L
  ))
  expect_s3_class(unsafe, "stormbreaker_unadmitted_pipeline")

  d <- make_sb_binary_fixture(n = 90L, seed = 8911L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(learner = "cancerclass_binary")
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 8912L)

  expect_identical(pipeline$prediction_type, "probability")
  expect_identical(pipeline$learner_params$feature_method, "welch.test")
  expect_identical(pipeline$learner_params$calibration, "oof_logistic")
  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_true(all(is.finite(result$predictions$probability)))
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("success", 3L)
  )
  expect_identical(
    result$provenance$cancerclass_version,
    as.character(utils::packageVersion("cancerclass"))
  )
  expect_identical(
    result$provenance$biobase_version,
    as.character(utils::packageVersion("Biobase"))
  )

  fit <- sb_fit_final(task, pipeline)
  new_data <- d[c(7L, 28L, 74L),
    c("patient_id", "outcome", rev(predictors)),
    drop = FALSE
  ]
  altered <- new_data
  altered$outcome <- 1L - altered$outcome
  original_prediction <- sb_predict_final(fit, new_data, "patient_id")
  altered_prediction <- sb_predict_final(fit, altered, "patient_id")

  expect_identical(fit$status, "success")
  expect_identical(fit$learner_state$learner, "cancerclass_binary")
  expect_identical(original_prediction$status, "success")
  expect_equal(
    original_prediction$predictions$probability,
    altered_prediction$predictions$probability,
    tolerance = 0
  )
})

test_that("cancerclass fails explicitly without dependencies or feasible folds", {
  x <- data.frame(
    gene_a = seq_len(12L),
    gene_b = rev(seq_len(12L)),
    gene_c = rep(c(-1, 1), length.out = 12L)
  )
  adapter <- sb_cancerclass_binary_learner_adapter()

  local_mocked_bindings(sb_cancerclass_available = function() FALSE)
  missing <- capture_sb_error(adapter$fit(x, rep(c(0L, 1L), each = 6L)))
  expect_s3_class(missing, "stormbreaker_missing_cancerclass")

  local_mocked_bindings(sb_cancerclass_available = function() TRUE)
  single_class <- capture_sb_error(adapter$fit(x, rep(0L, nrow(x))))
  infeasible <- capture_sb_error(adapter$fit(
    x[seq_len(8L), , drop = FALSE],
    rep(c(0L, 1L), each = 4L)
  ))
  expect_s3_class(single_class, "stormbreaker_cancerclass_single_class")
  expect_s3_class(infeasible, "stormbreaker_infeasible_inner_folds")
})

# ---- high-dimensional-discriminant-adapters ----

test_that("PDA and stabilised LDA declare distinct native specifications", {
  expected_dependencies <- c(
    pda_binary = "mda",
    slda_binary = "ipred"
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
    sb_classical_binary_params("pda_binary"),
    list(
      method = "gen.ridge",
      dimension = 1L,
      lambda = 1,
      prior = "empirical",
      keep_fitted = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("slda_binary"),
    list(
      score_rule = "native_PC_q",
      q = NULL,
      prior = "empirical"
    )
  )

  registry <- sb_pipeline_registry()
  selectors <- expected_fixed_learner_selectors("binary")
  for (learner in names(expected_dependencies)) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(nrow(rows), length(selectors))
    expect_setequal(rows$selector, selectors)
    expect_identical(unique(rows$task_type), "binary")
    expect_identical(unique(rows$tuning_semantics), "none")
  }
})

test_that("PDA and stabilised LDA reproduce native high-dimensional probabilities", {
  skip_if_not_installed("mda", minimum_version = "0.5.0")
  skip_if_not_installed("ipred", minimum_version = "0.9.0")

  data <- make_sb_high_dim_binary_fixture(n = 64L, p = 90L, seed = 8611L)
  names(data)[names(data) == "gene_001"] <- "gene A"
  names(data)[names(data) == "gene_002"] <- "gene-B"
  predictors <- setdiff(names(data), c("patient_id", "outcome"))
  x <- data[predictors]
  assessment <- x[1:9, , drop = FALSE]
  grouping <- factor(data$outcome, levels = c(0L, 1L))
  priors <- as.numeric(table(grouping) / length(grouping))

  pda_adapter <- sb_classical_binary_learner_adapter(
    "pda_binary"
  )
  pda <- pda_adapter$fit(x, data$outcome)
  formula_data <- sb_classical_binary_formula_data(x, data$outcome)
  native_pda <- mda::fda(
    formula_data$formula,
    data = formula_data$data,
    dimension = 1L,
    method = mda::gen.ridge,
    keep.fitted = FALSE,
    lambda = 1
  )
  expected_pda <- stats::predict(
    native_pda,
    newdata = assessment,
    type = "posterior"
  )[, "1"]

  expect_s3_class(pda$engine, "fda")
  expect_s3_class(pda$engine$fit, "gen.ridge")
  expect_equal(as.numeric(pda$engine$prior), priors, tolerance = 0)
  expect_identical(pda$engine$fit$lambda, 1)
  expect_gt(ncol(x), nrow(x))
  expect_equal(
    pda_adapter$predict(pda, assessment),
    unname(expected_pda),
    tolerance = 1e-13
  )

  slda_adapter <- sb_classical_binary_learner_adapter(
    "slda_binary"
  )
  slda <- slda_adapter$fit(x, data$outcome)
  native_slda <- ipred::slda(
    grouping,
    X = as.matrix(x),
    q = NULL,
    prior = priors
  )
  projected <- as.matrix(assessment) %*% native_slda$scores
  colnames(projected) <- paste0("V", seq_len(ncol(projected)))
  expected_slda <- stats::predict(
    native_slda$mylda,
    newdata = as.data.frame(projected)
  )$posterior[, "1"]

  expect_s3_class(slda$engine, "slda")
  expect_s3_class(slda$engine$mylda, "lda")
  expect_identical(
    ncol(slda$engine$scores),
    ncol(native_slda$scores)
  )
  expect_equal(unname(slda$engine$mylda$prior), priors, tolerance = 0)
  expect_equal(
    slda_adapter$predict(slda, assessment),
    unname(expected_slda),
    tolerance = 1e-13
  )
})

test_that("PDA and stabilised LDA states replay without fitted patient rows", {
  skip_if_not_installed("mda", minimum_version = "0.5.0")
  skip_if_not_installed("ipred", minimum_version = "0.9.0")

  data <- make_sb_high_dim_binary_fixture(n = 64L, p = 90L, seed = 8612L)
  names(data)[names(data) == "gene_001"] <- "gene A"
  names(data)[names(data) == "gene_002"] <- "gene-B"
  predictors <- setdiff(names(data), c("patient_id", "outcome"))
  x <- data[predictors]
  assessment <- x[1:7, , drop = FALSE]

  for (learner in c("pda_binary", "slda_binary")) {
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
    if (identical(learner, "pda_binary")) {
      expect_identical(
        intersect(c("call", "confusion"), names(engine)),
        character()
      )
      expect_identical(
        intersect(
          c("fitted.values", "x", "residuals", "call"),
          names(engine$fit)
        ),
        character()
      )
      expect_identical(environment(engine$terms), baseenv())
    } else {
      expect_identical(names(engine), c("scores", "mylda"))
      expect_identical("call" %in% names(engine$mylda), FALSE)
    }
  }
})

test_that("PDA and stabilised LDA execute through CV and final prediction", {
  skip_if_not_installed("mda", minimum_version = "0.5.0")
  skip_if_not_installed("ipred", minimum_version = "0.9.0")

  data <- make_sb_high_dim_binary_fixture(n = 72L, p = 12L, seed = 8613L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    predictors = sprintf("gene_%03d", seq_len(12L))
  )

  for (learner in c("pda_binary", "slda_binary")) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 4L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 8614L)
    final <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      final,
      data[1:8, c("patient_id", sprintf("gene_%03d", seq_len(12L)))],
      id_col = "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(final$status, "success")
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$probability, 8L)
  }
})

test_that("PDA no-discrimination fits fail before prediction", {
  skip_if_not_installed("mda", minimum_version = "0.5.0")

  data <- data.frame(
    patient_id = sprintf("PD%02d", seq_len(30L)),
    outcome = rep(c(0L, 1L), each = 15L),
    constant_gene = 0,
    check.names = FALSE
  )
  adapter <- sb_classical_binary_learner_adapter("pda_binary")
  error <- capture_sb_error(suppressWarnings(adapter$fit(
    data["constant_gene"],
    data$outcome
  )))
  task <- sb_binary_task(
    data,
    id_col = "patient_id",
    outcome_col = "outcome",
    predictors = "constant_gene"
  )
  pipeline <- sb_pipeline(learner = "pda_binary")
  comparison <- suppressWarnings(sb_run_cv(
    task,
    pipeline,
    folds = 3L,
    seed = 8615L
  ))
  final <- suppressWarnings(sb_fit_final(task, pipeline))

  expect_s3_class(error, "stormbreaker_pda_no_discrimination")
  expect_s3_class(error$parent, "warning")
  expect_identical(comparison$status, "failed")
  expect_identical(comparison$complete_oof, FALSE)
  expect_identical(
    unique(comparison$runs$error_class),
    "stormbreaker_pda_no_discrimination"
  )
  expect_identical(final$status, "failed")
  expect_identical(final$error_stage, "learner_fit")
  expect_identical(
    final$error_class,
    "stormbreaker_pda_no_discrimination"
  )
})

# ---- high-dimensional-omics-adapters ----

make_high_dimensional_omics_fixture <- function(n = 72L, p = 48L,
                                                seed = 9201L) {
  withr::with_seed(seed, {
    outcome <- sample(rep(c(0L, 1L), length.out = n))
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    x[, 1L] <- x[, 1L] + 1.4 * outcome
    x[, 2L] <- x[, 2L] - outcome
    x[, 3L] <- x[, 3L] + 0.7 * outcome
    x <- as.data.frame(x, optional = TRUE)
    names(x) <- paste0("gene_", seq_len(p))
    list(x = x, outcome = outcome)
  })
}

make_high_dimensional_survival_fixture <- function(n = 84L, p = 48L,
                                                   seed = 9203L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    event_time <- stats::rexp(
      n,
      rate = exp(1.1 * x[, 1L] - 0.9 * x[, 2L] + 0.6 * x[, 3L])
    )
    censor_time <- stats::rexp(n, rate = 0.3)
    x <- as.data.frame(x, optional = TRUE)
    names(x) <- paste0("gene_", seq_len(p))
    list(
      x = x,
      time = pmin(event_time, censor_time) + 0.01,
      event = as.integer(event_time <= censor_time)
    )
  })
}

test_that("SIS screening uses its fold-local native marginal ranking", {
  skip_if_not_installed("SIS", minimum_version = "1.5")

  fixture <- make_high_dimensional_omics_fixture(p = 80L)
  params <- sb_sis_binary_selector_params()
  expected_size <- sb_sis_binary_screen_size(
    nrow(fixture$x),
    ncol(fixture$x)
  )
  native <- withr::with_seed(params$seed, SIS::SIS(
    x = data.matrix(fixture$x),
    y = fixture$outcome,
    family = params$family,
    penalty = params$penalty,
    tune = params$tuning,
    type.measure = "deviance",
    nsis = expected_size,
    iter = FALSE,
    varISIS = "vanilla",
    perm = FALSE,
    greedy = FALSE,
    seed = params$seed,
    standardize = params$standardize,
    boot_ci = FALSE,
    parallel = FALSE
  ))

  set.seed(9202L)
  caller_rng <- .Random.seed
  fitted <- sb_fit_sis_binary_selector(fixture$x, fixture$outcome)

  expect_identical(.Random.seed, caller_rng)
  expect_identical(
    fitted$selected_features,
    names(fixture$x)[native$sis.ix0]
  )
  expect_identical(length(fitted$selected_features), expected_size)
  expect_identical(
    names(sb_apply_high_dimensional_selector(fitted, fixture$x)),
    fitted$selected_features
  )
})

test_that("VSURF exposes its native prediction-stage subset", {
  skip_if_not_installed("VSURF", minimum_version = "1.2.1")

  fixture <- make_high_dimensional_omics_fixture(p = 30L, seed = 9211L)
  params <- sb_vsurf_binary_selector_params()
  mtry <- max(floor(ncol(fixture$x) / 3), 1L)
  native <- withr::with_seed(params$seed, VSURF::VSURF(
    x = data.matrix(fixture$x),
    y = factor(fixture$outcome, levels = c(0L, 1L)),
    mtry = mtry,
    ntree.thres = params$threshold_trees,
    nfor.thres = params$threshold_forests,
    nmin = params$nmin,
    ntree.interp = params$interpretation_trees,
    nfor.interp = params$interpretation_forests,
    nsd = params$nsd,
    ntree.pred = params$prediction_trees,
    nfor.pred = params$prediction_forests,
    nmj = params$nmj,
    RFimplem = params$random_forest_engine,
    parallel = FALSE,
    verbose = FALSE
  ))

  set.seed(9212L)
  caller_rng <- .Random.seed
  fitted <- sb_fit_vsurf_binary_selector(fixture$x, fixture$outcome)
  diagnostic <- sb_selector_diagnostics(fitted, "vsurf_importance")
  expect_equal(diagnostic$table$importance[native$imp.mean.dec.ind], as.numeric(native$imp.mean.dec), tolerance = 0)
  expect_equal(diagnostic$table$importance_sd[native$imp.mean.dec.ind], as.numeric(native$imp.sd.dec), tolerance = 0)

  expect_identical(.Random.seed, caller_rng)
  expect_identical(
    fitted$selected_features,
    names(fixture$x)[native$varselect.pred]
  )
  expect_identical(
    fitted$selection_history$n_features,
    as.integer(c(
      length(native$varselect.thres),
      length(native$varselect.interp),
      length(native$varselect.pred)
    ))
  )
})

test_that("mRMR preserves its native ordered high-dimensional subset", {
  skip_if_not_installed("praznik", minimum_version = "13.0.0")

  fixture <- make_high_dimensional_omics_fixture(p = 90L, seed = 9213L)
  params <- sb_mrmr_binary_selector_params()
  expected_size <- sb_high_dimensional_nlogn_size(
    nrow(fixture$x),
    ncol(fixture$x)
  )
  native <- praznik::MRMR(
    X = fixture$x,
    Y = factor(fixture$outcome, levels = c(0L, 1L)),
    k = expected_size,
    positive = params$positive_scores_only,
    threads = params$threads
  )

  set.seed(9214L)
  caller_rng <- .Random.seed
  fitted <- sb_fit_mrmr_binary_selector(fixture$x, fixture$outcome)

  expect_identical(.Random.seed, caller_rng)
  expect_identical(
    fitted$selected_features,
    names(fixture$x)[as.integer(native$selection)]
  )
  expect_equal(
    fitted$scores$mrmr_score[
      order(fitted$scores$rank, na.last = NA)
    ],
    as.numeric(native$score),
    tolerance = 0
  )
})

test_that("information filters preserve their native ordered subsets", {
  skip_if_not_installed("praznik", minimum_version = "13.0.0")

  fixture <- make_high_dimensional_omics_fixture(p = 90L, seed = 9215L)
  expected_size <- sb_high_dimensional_nlogn_size(
    nrow(fixture$x),
    ncol(fixture$x)
  )
  native_methods <- c(
    jmi_binary_nlogn = "JMI",
    cmim_binary_nlogn = "CMIM",
    disr_binary_nlogn = "DISR",
    jmim_binary_nlogn = "JMIM"
  )

  for (selector in names(native_methods)) {
    params <- sb_information_binary_selector_params(selector)
    native <- do.call(
      getExportedValue("praznik", native_methods[[selector]]),
      list(
        X = fixture$x,
        Y = factor(fixture$outcome, levels = c(0L, 1L)),
        k = expected_size,
        threads = params$threads
      )
    )

    set.seed(9216L)
    caller_rng <- .Random.seed
    fitted <- sb_fit_information_binary_selector(
      fixture$x,
      fixture$outcome,
      selector
    )

    expect_identical(.Random.seed, caller_rng)
    expect_identical(
      fitted$selected_features,
      names(fixture$x)[as.integer(native$selection)]
    )
    expect_equal(
      fitted$scores$information_score[
        order(fitted$scores$rank, na.last = NA)
      ],
      as.numeric(native$score),
      tolerance = 0
    )
    expect_identical(fitted$params$native_method, native_methods[[selector]])
  }
})

test_that("LPC selectors preserve native gene scores and absolute ranking", {
  skip_if_not_installed("lpc", minimum_version = "1.0.2.1")

  binary <- make_high_dimensional_omics_fixture(
    n = 80L,
    p = 60L,
    seed = 9233L
  )
  survival <- make_high_dimensional_survival_fixture(
    n = 84L,
    p = 60L,
    seed = 9234L
  )
  specifications <- list(
    list(
      selector = "lpc_binary_nlogn",
      x = binary$x,
      response = binary$outcome + 1L,
      censoring_status = NULL,
      fit = function() {
        sb_fit_lpc_binary_selector(binary$x, binary$outcome)
      }
    ),
    list(
      selector = "lpc_cox_nlogn",
      x = survival$x,
      response = survival$time,
      censoring_status = survival$event,
      fit = function() {
        sb_fit_lpc_cox_selector(
          survival$x,
          survival$time,
          survival$event
        )
      }
    )
  )

  for (specification in specifications) {
    params <- sb_lpc_selector_params(specification$selector)
    native <- NULL
    invisible(utils::capture.output(
      native <- withr::with_seed(
        params$seed,
        lpc::LPC(
          x = t(data.matrix(specification$x)),
          y = specification$response,
          type = params$native_type,
          soft.thresh = NULL,
          censoring.status = specification$censoring_status
        )
      )
    ))
    native_score <- as.numeric(native$lpcscores)
    expected_size <- sb_high_dimensional_nlogn_size(
      nrow(specification$x),
      ncol(specification$x)
    )
    expected_order <- order(
      -abs(native_score),
      seq_along(native_score),
      method = "radix"
    )

    set.seed(9235L)
    caller_rng <- .Random.seed
    fitted <- specification$fit()

    expect_identical(.Random.seed, caller_rng)
    expect_identical(
      fitted$selected_features,
      names(specification$x)[expected_order[seq_len(expected_size)]]
    )
    expect_equal(fitted$scores$lpc_score, native_score, tolerance = 0)
    expect_equal(
      fitted$scores$base_score,
      as.numeric(native$tscores),
      tolerance = 0
    )
    expect_identical(
      fitted$params$realized_soft_threshold,
      as.numeric(native$soft.thresh)
    )
  }
})

test_that("LPC fails transparently when its native threshold is unavailable", {
  skip_if_not_installed("lpc", minimum_version = "1.0.2.1")

  fixture <- make_high_dimensional_omics_fixture(
    n = 60L,
    p = 20L,
    seed = 9236L
  )

  expect_error(
    sb_fit_lpc_binary_selector(fixture$x, fixture$outcome),
    class = "stormbreaker_invalid_lpc_selector"
  )
})

test_that("stability selectors reproduce native selected variables", {
  skip_if_not_installed("stabs", minimum_version = "0.7-1")
  skip_if_not_installed("glmnet", minimum_version = "4.1-2")

  binary <- make_high_dimensional_omics_fixture(p = 60L, seed = 9216L)
  survival <- make_high_dimensional_survival_fixture(
    p = 60L,
    seed = 9217L
  )
  specifications <- list(
    list(
      selector = "stabsel_binary_lasso",
      x = binary$x,
      response = binary$outcome,
      fit = function() {
        sb_fit_stabsel_binary_selector(binary$x, binary$outcome)
      }
    ),
    list(
      selector = "stabsel_cox_lasso",
      x = survival$x,
      response = survival::Surv(survival$time, survival$event),
      fit = function() {
        sb_fit_stabsel_cox_selector(
          survival$x,
          survival$time,
          survival$event
        )
      }
    )
  )

  for (specification in specifications) {
    params <- sb_stabsel_selector_params(specification$selector)
    subsample_q <- sb_stabsel_subsample_q(
      nrow(specification$x),
      ncol(specification$x)
    )
    fitfun_args <- list(
      type = params$base_selection_type,
      family = params$family,
      standardize = params$standardize
    )
    if (identical(params$family, "cox")) {
      fitfun_args$cox.ties <- params$cox_ties
    }
    native <- withr::with_seed(params$seed, stabs::stabsel(
      x = data.matrix(specification$x),
      y = specification$response,
      fitfun = stabs::glmnet.lasso,
      args.fitfun = fitfun_args,
      cutoff = params$cutoff,
      q = subsample_q,
      B = params$complementary_pairs,
      sampling.type = params$sampling_type,
      assumption = params$assumption,
      papply = parallel::mclapply,
      mc.cores = 1L,
      verbose = FALSE
    ))

    set.seed(9218L)
    caller_rng <- .Random.seed
    fitted <- specification$fit()

    expect_identical(.Random.seed, caller_rng)
    expect_identical(
      fitted$selected_features,
      names(specification$x)[unname(as.integer(stabs::selected(native)))]
    )
    expect_equal(
      fitted$scores$selection_probability,
      unname(native$max[names(specification$x)]),
      tolerance = 0
    )
  }
})

test_that("Cox stability selection locks glmnet ties to Breslow", {
  skip_if_not_installed("stabs", minimum_version = "0.7-1")
  skip_if_not_installed("glmnet", minimum_version = "4.1-2")

  fixture <- make_high_dimensional_survival_fixture(
    n = 30L,
    p = 6L,
    seed = 9219L
  )
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    stabsel = function(x, y, ..., args.fitfun) {
      captured$args.fitfun <- args.fitfun
      structure(
        list(
          selected = 1L,
          max = stats::setNames(
            c(1, rep(0, ncol(x) - 1L)),
            colnames(x)
          ),
          PFER = 0
        ),
        class = "stabsel"
      )
    },
    .package = "stabs"
  )

  sb_fit_stabsel_cox_selector(
    fixture$x,
    fixture$time,
    fixture$event
  )

  expect_identical(
    captured$args.fitfun,
    list(
      type = "conservative",
      family = "cox",
      standardize = TRUE,
      cox.ties = "breslow"
    )
  )
})

test_that("HDDA predicts native class-one posteriors from p greater than n", {
  skip_if_not_installed("HDclassif", minimum_version = "2.2.2")

  fixture <- make_high_dimensional_omics_fixture(n = 60L, p = 90L)
  adapter <- sb_hdda_binary_learner_adapter("hdda_binary_cattell")
  model <- adapter$fit(fixture$x, fixture$outcome)
  assessment <- fixture$x[1:7, , drop = FALSE]
  native <- stats::predict(model$engine, data = data.matrix(assessment))
  observed <- adapter$predict(model, assessment)

  expect_equal(observed, unname(native$posterior[, "1"]), tolerance = 1e-12)
  state <- unserialize(serialize(adapter$extract_state(model), NULL))
  expect_equal(
    adapter$predict_state(state, assessment)$probability,
    observed,
    tolerance = 1e-12
  )
})

test_that("VSURF and HDDA complete final fit and outcome-blind prediction", {
  skip_if_not_installed("VSURF", minimum_version = "1.2.1")
  skip_if_not_installed("HDclassif", minimum_version = "2.2.2")

  fixture <- make_high_dimensional_omics_fixture(p = 24L, seed = 9215L)
  data <- data.frame(
    patient_id = sprintf("HF%03d", seq_len(nrow(fixture$x))),
    outcome = fixture$outcome,
    fixture$x,
    check.names = FALSE
  )
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    names(fixture$x)
  )
  new_data <- data[1:6, c("patient_id", names(fixture$x))]
  new_data$patient_id <- paste0("NEW_", new_data$patient_id)
  fits <- list(
    sb_fit_final(task, sb_pipeline(
      selector = "vsurf_binary_prediction",
      learner = "logistic"
    )),
    sb_fit_final(task, sb_pipeline(learner = "hdda_binary_cattell"))
  )

  for (fit in fits) {
    expect_identical(fit$status, "success")
    prediction <- sb_predict_final(fit, new_data, "patient_id")
    expect_identical(prediction$status, "success")
    expect_identical(nrow(prediction$predictions), nrow(new_data))
    expect_gte(min(prediction$predictions$probability), 0)
    expect_lte(max(prediction$predictions$probability), 1)
  }
  expect_identical(
    fits[[1L]]$provenance$vsurf_version,
    as.character(utils::packageVersion("VSURF"))
  )
  expect_identical(
    fits[[2L]]$provenance$hdclassif_version,
    as.character(utils::packageVersion("HDclassif"))
  )
})

test_that("information and stability selectors complete final prediction", {
  skip_if_not_installed("praznik", minimum_version = "13.0.0")
  skip_if_not_installed("stabs", minimum_version = "0.7-1")
  skip_if_not_installed("glmnet", minimum_version = "4.1-2")

  binary <- make_high_dimensional_omics_fixture(
    n = 84L,
    p = 40L,
    seed = 9223L
  )
  binary_data <- data.frame(
    patient_id = sprintf("MB%03d", seq_len(nrow(binary$x))),
    outcome = binary$outcome,
    binary$x,
    check.names = FALSE
  )
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    names(binary$x)
  )
  binary_new <- binary_data[1:6, c("patient_id", names(binary$x))]
  binary_new$patient_id <- paste0("NEW_", binary_new$patient_id)

  survival <- make_high_dimensional_survival_fixture(
    n = 90L,
    p = 40L,
    seed = 9224L
  )
  survival_data <- data.frame(
    patient_id = sprintf("MS%03d", seq_len(nrow(survival$x))),
    time = survival$time,
    event = survival$event,
    survival$x,
    check.names = FALSE
  )
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    names(survival$x),
    "month",
    prediction_horizon = stats::median(survival$time)
  )
  survival_new <- survival_data[1:6, c("patient_id", names(survival$x))]
  survival_new$patient_id <- paste0("NEW_", survival_new$patient_id)

  binary_selectors <- c(
    "mrmr_binary_nlogn",
    sb_information_binary_selector_ids(),
    "stabsel_binary_lasso"
  )
  binary_fits <- lapply(binary_selectors, function(selector) {
    sb_fit_final(binary_task, sb_pipeline(
      selector = selector,
      learner = "gaussian_naive_bayes"
    ))
  })
  fits <- c(
    binary_fits,
    list(
      sb_fit_final(survival_task, sb_pipeline(
        selector = "stabsel_cox_lasso",
        learner = "coxph"
      ))
    )
  )
  new_data <- c(
    rep(list(binary_new), length(binary_fits)),
    list(survival_new)
  )

  for (index in seq_along(fits)) {
    expect_identical(fits[[index]]$status, "success")
    prediction <- sb_predict_final(
      fits[[index]],
      new_data[[index]],
      "patient_id"
    )
    expect_identical(prediction$status, "success")
    expect_identical(nrow(prediction$predictions), 6L)
  }
  praznik_indexes <- which(binary_selectors != "stabsel_binary_lasso")
  expect_identical(
    vapply(
      fits[praznik_indexes],
      function(fit) fit$provenance$praznik_version,
      character(1)
    ),
    rep(as.character(utils::packageVersion("praznik")), 5L)
  )
  expect_identical(
    fits[[length(binary_fits)]]$provenance$stabs_version,
    as.character(utils::packageVersion("stabs"))
  )
  expect_identical(
    fits[[length(fits)]]$provenance$glmnet_version,
    as.character(utils::packageVersion("glmnet"))
  )

  cv <- sb_run_cv(
    survival_task,
    pipeline = sb_pipeline(
      selector = "stabsel_cox_lasso",
      learner = "coxph"
    ),
    folds = 3L,
    seed = 9225L
  )
  expect_identical(cv$status, "success")
  expect_identical(
    cv$provenance$stabs_version,
    as.character(utils::packageVersion("stabs"))
  )
})

test_that("LPC selectors complete serializable outcome-blind prediction", {
  skip_if_not_installed("lpc", minimum_version = "1.0.2.1")
  skip_if_not_installed("e1071", minimum_version = "1.7-13")

  binary <- make_high_dimensional_omics_fixture(
    n = 84L,
    p = 60L,
    seed = 9237L
  )
  binary_data <- data.frame(
    patient_id = sprintf("LB%03d", seq_len(nrow(binary$x))),
    outcome = binary$outcome,
    binary$x,
    check.names = FALSE
  )
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    names(binary$x)
  )

  survival <- make_high_dimensional_survival_fixture(
    n = 90L,
    p = 60L,
    seed = 9238L
  )
  survival_data <- data.frame(
    patient_id = sprintf("LS%03d", seq_len(nrow(survival$x))),
    time = survival$time,
    event = survival$event,
    survival$x,
    check.names = FALSE
  )
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    names(survival$x),
    "month",
    prediction_horizon = stats::median(survival$time)
  )

  fits <- list(
    sb_fit_final(binary_task, sb_pipeline(
      selector = "lpc_binary_nlogn",
      learner = "gaussian_naive_bayes"
    )),
    sb_fit_final(survival_task, sb_pipeline(
      selector = "lpc_cox_nlogn",
      learner = "coxph"
    ))
  )
  new_data <- list(
    binary_data[1:6, c("patient_id", names(binary$x))],
    survival_data[1:6, c("patient_id", names(survival$x))]
  )

  for (index in seq_along(fits)) {
    fit <- unserialize(serialize(fits[[index]], NULL))
    prediction <- sb_predict_final(fit, new_data[[index]], "patient_id")

    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_identical(nrow(prediction$predictions), 6L)
    prediction_values <- if (index == 1L) {
      prediction$predictions$probability
    } else {
      prediction$predictions$risk
    }
    expect_true(all(is.finite(prediction_values)))
    if (index == 2L) {
      expect_true(all(is.finite(
        prediction$predictions$survival_probability
      )))
    }
    expect_identical(
      fit$provenance$lpc_version,
      as.character(utils::packageVersion("lpc"))
    )
  }
  expect_true(all(
    fits[[1L]]$selection$selected_features %in% names(binary$x)
  ))
  expect_true(all(
    fits[[2L]]$selection$selected_features %in% names(survival$x)
  ))
})

test_that("high-dimensional methods resolve through admitted binary routes", {
  registry <- sb_pipeline_registry()
  unrestricted <- registry$learner[
    registry$task_type == "binary" & registry$selector == "none" &
      !registry$learner %in% c(
        "pre_binary_rule_ensemble",
        sb_slope_binary_learner_ids(),
        "cancerclass_binary"
      )
  ]

  for (selector in setdiff(
    sb_high_dimensional_binary_selector_ids(),
    "stabsel_binary_lasso"
  )) {
    expect_setequal(
      registry$learner[registry$selector == selector],
      unrestricted
    )
  }
  expect_setequal(
    registry$learner[registry$selector == "stabsel_binary_lasso"],
    setdiff(unrestricted, "glmnet_logistic")
  )
  expect_setequal(
    registry$learner[registry$selector == "lpc_cox_nlogn"],
    registry$learner[registry$selector == "boruta_survival_roughfix"]
  )
  expect_setequal(
    registry$learner[registry$selector == "stabsel_cox_lasso"],
    setdiff(
      registry$learner[registry$selector == "boruta_survival_roughfix"],
      "glmnet_cox"
    )
  )
  expect_setequal(
    registry$selector[registry$learner == "hdda_binary_cattell"],
    expected_fixed_learner_selectors("binary")
  )

  binary_specifications <- lapply(
    sb_high_dimensional_binary_selector_ids(),
    function(selector) {
      sb_pipeline(selector = selector, learner = "logistic")
    }
  )
  specifications <- c(
    binary_specifications,
    list(
      sb_pipeline(selector = "lpc_cox_nlogn", learner = "coxph"),
      sb_pipeline(selector = "stabsel_cox_lasso", learner = "coxph"),
      sb_pipeline(learner = "hdda_binary_cattell")
    )
  )
  for (pipeline in specifications) {
    runtime <- sb_resolve_pipeline_runtime(pipeline)
    expect_identical(sb_validate_pipeline_runtime(runtime, pipeline), runtime)
  }
})

test_that("SIS is refitted inside glmnet inner analysis folds", {
  skip_if_not_installed("SIS", minimum_version = "1.5")
  skip_if_not_installed("glmnet", minimum_version = "4.1.2")

  fixture <- make_high_dimensional_omics_fixture(n = 84L, p = 36L)
  data <- data.frame(
    patient_id = sprintf("HD%03d", seq_len(nrow(fixture$x))),
    outcome = fixture$outcome,
    fixture$x,
    check.names = FALSE
  )
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    names(fixture$x)
  )
  fit <- sb_fit_final(task, sb_pipeline(
    selector = "sis_binary_nlogn",
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 9221L
  ))

  expect_identical(fit$status, "success")
  expect_identical(
    fit$tuning$inner_selection_method,
    "glmnet_exclude_function"
  )
  expect_length(fit$tuning$inner_selections, 3L)
  expect_identical(
    fit$provenance$sis_version,
    as.character(utils::packageVersion("SIS"))
  )
})

test_that("praznik selectors are refitted inside glmnet inner folds", {
  skip_if_not_installed("praznik", minimum_version = "13.0.0")
  skip_if_not_installed("glmnet", minimum_version = "4.1-2")

  binary <- make_high_dimensional_omics_fixture(
    n = 90L,
    p = 36L,
    seed = 9231L
  )
  binary_data <- data.frame(
    patient_id = sprintf("IB%03d", seq_len(nrow(binary$x))),
    outcome = binary$outcome,
    binary$x,
    check.names = FALSE
  )
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    names(binary$x)
  )
  selectors <- c("mrmr_binary_nlogn", sb_information_binary_selector_ids())
  for (selector in selectors) {
    fit <- sb_fit_final(binary_task, sb_pipeline(
      selector = selector,
      learner = "glmnet_logistic",
      inner_folds = 3L,
      inner_seed = 9232L
    ))

    expect_identical(fit$status, "success")
    expect_identical(
      fit$tuning$inner_selection_method,
      "glmnet_exclude_function"
    )
    expect_length(fit$tuning$inner_selections, 3L)
    expect_gt(
      min(vapply(
        fit$tuning$inner_selections,
        function(selection) length(selection$selected_features),
        integer(1)
      )),
      0L
    )
  }
})

test_that("LPC selectors are refitted inside glmnet inner folds", {
  skip_if_not_installed("lpc", minimum_version = "1.0.2.1")
  skip_if_not_installed("glmnet", minimum_version = "4.1-2")

  binary <- make_high_dimensional_omics_fixture(
    n = 90L,
    p = 60L,
    seed = 9239L
  )
  binary_data <- data.frame(
    patient_id = sprintf("LI%03d", seq_len(nrow(binary$x))),
    outcome = binary$outcome,
    binary$x,
    check.names = FALSE
  )
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    names(binary$x)
  )

  survival <- make_high_dimensional_survival_fixture(
    n = 90L,
    p = 60L,
    seed = 9240L
  )
  survival_data <- data.frame(
    patient_id = sprintf("LC%03d", seq_len(nrow(survival$x))),
    time = survival$time,
    event = survival$event,
    survival$x,
    check.names = FALSE
  )
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    names(survival$x),
    "month",
    prediction_horizon = stats::median(survival$time)
  )

  fits <- list(
    sb_fit_final(binary_task, sb_pipeline(
      selector = "lpc_binary_nlogn",
      learner = "glmnet_logistic",
      inner_folds = 3L,
      inner_seed = 9241L
    )),
    sb_fit_final(survival_task, sb_pipeline(
      selector = "lpc_cox_nlogn",
      learner = "glmnet_cox",
      inner_folds = 3L,
      inner_seed = 9242L
    ))
  )

  for (fit in fits) {
    expect_identical(fit$status, "success")
    expect_identical(
      fit$tuning$inner_selection_method,
      "glmnet_exclude_function"
    )
    expect_length(fit$tuning$inner_selections, 3L)
    expect_gt(
      min(vapply(
        fit$tuning$inner_selections,
        function(selection) length(selection$selected_features),
        integer(1)
      )),
      0L
    )
    expect_identical(
      fit$provenance$lpc_version,
      as.character(utils::packageVersion("lpc"))
    )
  }
})

# ---- mixomics-da-adapters ----

test_that("mixOmics PLS-DA learners declare two distinct fixed methods", {
  expected_dependencies <- c(
    plsda_binary = "mixOmics",
    splsda_binary = "mixOmics"
  )

  expect_setequal(
    sb_mixomics_da_learner_ids(),
    names(expected_dependencies)
  )
  expect_identical(
    sb_mixomics_da_params("plsda_binary"),
    list(
      method = "PLS-DA",
      requested_components = 2L,
      component_rule = "min_requested_features_analysis_rows_minus_one",
      scale = TRUE,
      tolerance = 1e-6,
      maximum_iterations = 100L,
      near_zero_variance_filter = FALSE,
      logratio = "none",
      keep_x_rule = "all_features",
      probability_mapping = "binary_softmax_of_native_dummy_scores"
    )
  )
  expect_identical(
    sb_mixomics_da_params("splsda_binary"),
    list(
      method = "sPLS-DA",
      requested_components = 2L,
      component_rule = "min_requested_features_analysis_rows_minus_one",
      scale = TRUE,
      tolerance = 1e-6,
      maximum_iterations = 100L,
      near_zero_variance_filter = FALSE,
      logratio = "none",
      keep_x_rule = "min_10_and_p_per_component",
      probability_mapping = "binary_softmax_of_native_dummy_scores"
    )
  )

  registry <- sb_pipeline_registry()
  for (learner in names(expected_dependencies)) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_setequal(
      rows$selector,
      expected_fixed_learner_selectors("binary")
    )
    expect_true(all(grepl("mixOmics", rows$optional_dependency, fixed = TRUE)))
    expect_length(unique(rows$task_type), 1L)
    expect_identical(unique(rows$task_type), "binary")

    pipelines <- list(
      sb_pipeline(learner = learner),
      sb_pipeline(
        selector = "univ_binary_auc_top_k",
        learner = learner,
        selector_n_features = 3L
      ),
      sb_pipeline(
        selector = "abess_binary_gic_nonzero",
        learner = learner
      ),
      sb_pipeline(
        selector = "boruta_binary_roughfix",
        learner = learner
      ),
      sb_pipeline(
        selector = "varselrf_binary_oob_1se",
        learner = learner
      ),
      sb_pipeline(
        selector = "fixed_features",
        learner = learner,
        selector_features = c("gene_a", "gene_b", "gene_c"),
        feature_set_id = "mixomics_fixed"
      )
    )
    expect_length(pipelines, 6L)
    expect_true(all(vapply(
      pipelines,
      inherits,
      logical(1),
      what = "stormbreaker_pipeline"
    )))
  }
})

test_that("mixOmics adapters reproduce native dummy scores and softmax", {
  skip_if_not_installed("mixOmics", minimum_version = "6.32.0")

  data <- make_sb_high_dim_binary_fixture(n = 90L, p = 12L, seed = 8641L)
  x <- data[, -(1:2), drop = FALSE]
  names(x) <- c("gene A", "gene-B", paste0("gene_", 3:12))
  rownames(x) <- data$patient_id
  assessment <- x[1:9, , drop = FALSE]
  grouping <- factor(data$outcome, levels = c(0L, 1L))

  for (learner in sb_mixomics_da_learner_ids()) {
    adapter <- sb_mixomics_da_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    native <- if (identical(learner, "plsda_binary")) {
      mixOmics::plsda(
        data.matrix(x), grouping,
        ncomp = 2L, scale = TRUE, tol = 1e-6, max.iter = 100L,
        near.zero.var = FALSE, logratio = "none", all.outputs = TRUE
      )
    } else {
      mixOmics::splsda(
        data.matrix(x), grouping,
        ncomp = 2L, keepX = c(10L, 10L), scale = TRUE,
        tol = 1e-6, max.iter = 100L, near.zero.var = FALSE,
        logratio = "none", all.outputs = TRUE
      )
    }
    native_scores <- stats::predict(
      native,
      data.matrix(assessment)
    )$predict[, , "dim2"]
    expected <- stats::plogis(
      native_scores[, "1"] - native_scores[, "0"]
    )

    expect_equal(
      adapter$predict(model, assessment),
      unname(expected),
      tolerance = 1e-12
    )
    expect_identical(
      adapter$predict(model, assessment) > 0.5,
      unname(native_scores[, "1"] > native_scores[, "0"])
    )
    expect_identical(model$audit$effective_components, 2L)
    expect_identical(
      model$audit$probability_mapping,
      "binary_softmax_of_native_dummy_scores"
    )
    if (identical(learner, "splsda_binary")) {
      expect_identical(
        vapply(model$audit$component_features, length, integer(1)),
        c(10L, 10L)
      )
    } else {
      expect_identical(model$audit$component_features, list())
    }
  }
})

test_that("mixOmics component count follows available dimensions", {
  skip_if_not_installed("mixOmics", minimum_version = "6.32.0")

  data <- make_sb_binary_fixture(n = 70L, seed = 8642L)
  x <- data.frame("gene A" = data$gene_a, check.names = FALSE)
  rownames(x) <- data$patient_id

  for (learner in sb_mixomics_da_learner_ids()) {
    adapter <- sb_mixomics_da_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    probability <- adapter$predict(model, x[1:5, , drop = FALSE])

    expect_identical(model$audit$effective_components, 1L)
    expect_identical(dim(model$engine$coefficients), c(1L, 2L))
    expect_length(probability, 5L)
    expect_true(all(is.finite(probability)))
  }
})

test_that("mixOmics states replay without analysis rows", {
  skip_if_not_installed("mixOmics", minimum_version = "6.32.0")

  data <- make_sb_high_dim_binary_fixture(n = 84L, p = 4L, seed = 8643L)
  x <- data[, -(1:2), drop = FALSE]
  names(x) <- c("gene A", "gene-B", "gene_c", "gene_d")
  rownames(x) <- data$patient_id
  assessment <- x[1:7, , drop = FALSE]

  for (learner in sb_mixomics_da_learner_ids()) {
    adapter <- sb_mixomics_da_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    expected <- adapter$predict(model, assessment)
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, assessment)$probability

    expect_equal(replay, expected, tolerance = 0)
    expect_identical(
      names(state$engine),
      c("coefficients", "x_center", "x_scale", "y_center", "y_scale")
    )
    expect_identical(
      dim(state$engine$coefficients),
      c(ncol(x), 2L)
    )
    expect_length(state$engine$x_center, ncol(x))
    expect_length(state$engine$x_scale, ncol(x))
    expect_length(state$engine$y_center, 2L)
    expect_length(state$engine$y_scale, 2L)
  }
})

test_that("mixOmics learners execute through CV and final prediction", {
  skip_if_not_installed("mixOmics", minimum_version = "6.32.0")

  data <- make_sb_high_dim_binary_fixture(n = 72L, p = 4L, seed = 8644L)
  task <- sb_binary_task(
    data = data,
    id_col = "patient_id",
    outcome_col = "outcome",
    predictors = c("gene_001", "gene_002", "gene_003", "gene_004")
  )

  for (learner in sb_mixomics_da_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 3L
    )
    cv <- sb_run_binary_cv(task, pipeline, folds = 3L, seed = 8645L)
    final <- sb_fit_final(task, pipeline)
    predicted <- sb_predict_final(
      final,
      data[1:6, c(
        "patient_id", "gene_001", "gene_002", "gene_003", "gene_004"
      )],
      id_col = "patient_id"
    )

    expect_identical(cv$status, "success")
    expect_identical(cv$complete_oof, TRUE)
    expect_identical(final$status, "success")
    expect_identical(predicted$status, "success")
    expect_length(predicted$predictions$probability, 6L)
    expect_true(all(is.finite(predicted$predictions$probability)))
    expect_true(all(
      predicted$predictions$probability >= 0 &
        predicted$predictions$probability <= 1
    ))
  }
})

# ---- omics-specialist-and-mtlr-adapters ----

make_new_method_binary_fixture <- function(n = 72L, p = 8L, seed = 1701L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n)
    colnames(x) <- paste0("gene_", seq_len(p))
    outcome <- as.integer(
      x[, 1L] - 0.7 * x[, 2L] + 0.3 * x[, 3L] +
        stats::rnorm(n, sd = 0.6) > 0
    )
    data.frame(
      patient_id = sprintf("NM%03d", seq_len(n)),
      outcome = outcome,
      x,
      check.names = FALSE
    )
  })
}

make_new_method_survival_fixture <- function(n = 72L, p = 5L,
                                             seed = 1703L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n)
    colnames(x) <- paste0("gene_", seq_len(p))
    event_time <- stats::rexp(n, rate = exp(0.45 * x[, 1L] - 0.3 * x[, 2L]))
    censor_time <- stats::rexp(n, rate = 0.35)
    data.frame(
      patient_id = sprintf("MS%03d", seq_len(n)),
      time = pmin(event_time, censor_time) + 0.05,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    )
  })
}

test_that("K-TSP keeps the declared class-one orientation", {
  skip_if_not_installed("switchBox", minimum_version = "1.44.0")

  data <- make_new_method_binary_fixture()
  x <- data[setdiff(names(data), c("patient_id", "outcome"))]
  adapter <- sb_omics_specialist_binary_learner_adapter("ktsp_binary")
  model <- adapter$fit(x, data$outcome)
  native <- switchBox::SWAP.Train.KTSP(
    inputMat = t(data.matrix(x)),
    phenoGroup = factor(data$outcome, levels = c(0L, 1L)),
    krange = seq_len(4L),
    FilterFunc = switchBox::SWAP.Filter.Wilcoxon,
    handleTies = FALSE,
    disjoint = TRUE,
    featureNo = 8L,
    verbose = FALSE
  )
  assessment <- x[1:7, , drop = FALSE]
  native_statistic <- switchBox::SWAP.KTSP.Statistics(
    t(data.matrix(assessment)),
    native
  )$statistics
  expected <- unname(native_statistic / nrow(native$TSPs) + 0.5)

  expect_identical(model$engine$labels, c("0", "1"))
  expect_equal(model$engine$pairs, unname(native$TSPs), tolerance = 0)
  expect_equal(model$engine$pair_scores, unname(as.numeric(native$score)), tolerance = 0)
  expect_equal(adapter$predict(model, assessment), expected, tolerance = 0)
  expect_identical(
    adapter$predict(model, assessment) > 0.5,
    unname(switchBox::SWAP.KTSP.Classify(
      t(data.matrix(assessment)),
      native
    ) == "1")
  )
})

test_that("MTLR accepts its native interval-by-coefficient weight matrix", {
  skip_if_not_installed("MTLR", minimum_version = "0.2.2")

  data <- make_new_method_survival_fixture()
  x <- data[setdiff(names(data), c("patient_id", "time", "event"))]
  adapter <- sb_mtlr_survival_learner_adapter()
  model <- adapter$fit(x, data$time, data$event)

  expect_identical(ncol(model$engine$weight_matrix), ncol(x) + 1L)
  expect_identical(
    nrow(model$engine$weight_matrix),
    length(model$engine$time_points)
  )
})

test_that("new methods add their generated legal routes", {
  registry <- sb_pipeline_registry()
  binary_selector_count <- length(expected_fixed_learner_selectors("binary"))
  survival_selector_count <- length(expected_fixed_learner_selectors("survival"))
  expected_task <- c(
    stats::setNames(
      rep("binary", length(sb_omics_specialist_binary_learner_ids())),
      sb_omics_specialist_binary_learner_ids()
    ),
    mtlr_survival = "survival"
  )
  expected_routes <- c(
    stats::setNames(
      rep(
        binary_selector_count,
        length(sb_omics_specialist_binary_learner_ids())
      ),
      sb_omics_specialist_binary_learner_ids()
    ),
    mtlr_survival = survival_selector_count
  )

  for (learner in names(expected_task)) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(nrow(rows), unname(expected_routes[[learner]]))
    expect_identical(unique(rows$task_type), unname(expected_task[[learner]]))
    expect_identical(anyDuplicated(rows$selector), 0L)
    expect_true(all(rows$admission_status == "admitted"))
    expect_s3_class(sb_pipeline(learner = learner), "stormbreaker_pipeline")
  }
})

test_that("bias-reduced logistic remains finite under complete separation", {
  skip_if_not_installed("brglm2", minimum_version = "1.1.0")

  n <- 40L
  separating_feature <- seq(-1, 1, length.out = n)
  data <- data.frame(
    patient_id = sprintf("NM%03d", seq_len(n)),
    outcome = as.integer(separating_feature > 0),
    gene_1 = separating_feature,
    gene_2 = sin(seq_len(n)),
    gene_3 = cos(seq_len(n) / 3),
    gene_4 = rep(c(-1, 1), length.out = n)
  )
  x <- data[setdiff(names(data), c("patient_id", "outcome"))]
  adapter <- sb_omics_specialist_binary_learner_adapter(
    "brglm2_logistic"
  )
  model <- adapter$fit(x, data$outcome)
  probability <- adapter$predict(model, x[1:8, , drop = FALSE])
  state <- unserialize(serialize(adapter$extract_state(model), NULL))

  expect_true(all(is.finite(model$engine$coefficient)))
  expect_true(all(probability > 0 & probability < 1))
  expect_equal(
    adapter$predict_state(state, x[1:8, , drop = FALSE])$probability,
    probability,
    tolerance = 0
  )
  expect_identical(state$training_data_retention, "none")
})

test_that("Bayesian Cauchy logistic replays its native posterior mode", {
  skip_if_not_installed("arm", minimum_version = "1.15.3")

  data <- make_new_method_binary_fixture(n = 48L, p = 70L, seed = 1712L)
  x <- data[setdiff(names(data), c("patient_id", "outcome"))]
  adapter <- sb_omics_specialist_binary_learner_adapter(
    "bayesglm_logistic_cauchy"
  )
  params <- adapter$params
  set.seed(2719)
  rng_before <- .Random.seed
  model <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, rng_before)

  design <- cbind("(Intercept)" = 1, data.matrix(x))
  native <- suppressWarnings(arm::bayesglm.fit(
    x = design,
    y = data$outcome,
    family = stats::binomial(link = params$link),
    control = list(maxit = params$maximum_iterations),
    intercept = TRUE,
    prior.mean = params$prior_mean,
    prior.scale = params$prior_scale,
    prior.df = params$prior_df,
    prior.mean.for.intercept = params$intercept_prior_mean,
    prior.scale.for.intercept = params$intercept_prior_scale,
    prior.df.for.intercept = params$intercept_prior_df,
    scaled = params$scaled,
    Warning = FALSE
  ))
  assessment <- x[1:7, , drop = FALSE]
  expected <- stats::plogis(
    unname(native$coefficients[[1L]]) +
      as.vector(data.matrix(assessment) %*%
        unname(native$coefficients[-1L]))
  )
  state <- unserialize(serialize(adapter$extract_state(model), NULL))

  expect_equal(
    unname(model$engine$coefficient),
    unname(native$coefficients[-1L]),
    tolerance = 0
  )
  expect_equal(adapter$predict(model, assessment), expected, tolerance = 0)
  expect_equal(
    adapter$predict_state(state, assessment)$probability,
    expected,
    tolerance = 0
  )
  expect_identical(state$training_data_retention, "none")
})

test_that("hinge boosting stores a portable native margin", {
  skip_if_not_installed("bst", minimum_version = "0.3-24")

  data <- make_new_method_binary_fixture(seed = 1713L)
  x <- data[setdiff(names(data), c("patient_id", "outcome"))]
  adapter <- sb_omics_specialist_binary_learner_adapter(
    "bst_hinge_binary"
  )
  params <- adapter$params
  model <- adapter$fit(x, data$outcome)
  labels <- ifelse(data$outcome == 1L, 1, -1)
  native <- bst::bst(
    data.matrix(x),
    labels,
    cost = params$cost,
    family = params$family,
    learner = params$base_learner,
    ctrl = bst::bst_control(
      mstop = params$iterations,
      nu = params$learning_rate,
      center = params$center,
      intercept = params$fit_intercept,
      trace = FALSE
    )
  )
  assessment <- x[1:7, , drop = FALSE]
  native_margin <- as.numeric(stats::predict(
    native,
    newdata = data.matrix(assessment),
    type = "response"
  ))
  portable_margin <- model$engine$intercept +
    as.vector(data.matrix(assessment) %*% model$engine$coefficient)
  state <- unserialize(serialize(adapter$extract_state(model), NULL))

  expect_equal(portable_margin, native_margin, tolerance = 1e-12)
  expect_equal(
    adapter$predict(model, assessment),
    stats::plogis(2 * native_margin),
    tolerance = 1e-12
  )
  expect_equal(
    adapter$predict_state(state, assessment)$probability,
    adapter$predict(model, assessment),
    tolerance = 0
  )
  expect_identical(state$training_data_retention, "none")
})

test_that("HDA preserves its native one-dimensional posterior", {
  skip_if_not_installed("hda", minimum_version = "0.2-14")

  data <- make_new_method_binary_fixture(n = 96L, p = 8L, seed = 1714L)
  x <- data[setdiff(names(data), c("patient_id", "outcome"))]
  adapter <- sb_omics_specialist_binary_learner_adapter("hda_binary")
  params <- adapter$params
  model <- adapter$fit(x, data$outcome)
  native <- NULL
  invisible(utils::capture.output(
    native <- hda::hda(
      data.matrix(x),
      grouping = factor(data$outcome, levels = c(0L, 1L)),
      newdim = params$reduced_dimension,
      crule = TRUE,
      reg.lamb = params$regularization_lambda,
      reg.gamm = params$regularization_gamma,
      initial.loadings = NULL,
      noutit = params$outer_iterations,
      ninit = params$inner_iterations,
      verbose = FALSE
    ),
    type = "output"
  ))
  assessment <- x[1:7, , drop = FALSE]
  native_prediction <- suppressWarnings(stats::predict(
    native,
    data.matrix(assessment),
    task = "c"
  ))
  expected <- unname(native_prediction$posteriors[, 2L])
  state <- unserialize(serialize(adapter$extract_state(model), NULL))

  expect_equal(adapter$predict(model, assessment), expected, tolerance = 0)
  expect_equal(
    adapter$predict_state(state, assessment)$probability,
    expected,
    tolerance = 0
  )
  expect_setequal(
    names(model$engine),
    c("hda.loadings", "reduced.dimension", "naivebayes")
  )
  expect_identical(model$engine$naivebayes$levels, c("0", "1"))
  expect_identical(state$training_data_retention, "none")
})

test_that("DWD supports p greater than n and replays its native margin", {
  skip_if_not_installed("DWDLargeR", minimum_version = "0.2.0")

  data <- make_new_method_binary_fixture(n = 28L, p = 90L, seed = 1706L)
  x <- data[setdiff(names(data), c("patient_id", "outcome"))]
  x[] <- sweep(
    data.matrix(x),
    2L,
    10^seq(-3, 3, length.out = ncol(x)),
    "*"
  )
  adapter <- sb_omics_specialist_binary_learner_adapter("dwd_binary")
  set.seed(2718)
  rng_before <- .Random.seed
  model <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, rng_before)
  assessment <- x[1:6, , drop = FALSE]
  native_margin <- model$engine$intercept +
    as.vector(data.matrix(assessment) %*% model$engine$coefficient)
  probability <- adapter$predict(model, assessment)

  expect_equal(probability, stats::plogis(2 * native_margin), tolerance = 0)
  expect_true(is.finite(model$engine$penalty))
  training_probability <- adapter$predict(model, x)
  expect_equal(
    mean((training_probability > 0.5) != data$outcome) * 100,
    model$engine$training_error,
    tolerance = 1e-12
  )
  expect_equal(
    adapter$predict_state(
      unserialize(serialize(adapter$extract_state(model), NULL)),
      assessment
    )$probability,
    probability,
    tolerance = 0
  )
})

test_that("DWD maps removed all-zero features back to the original schema", {
  skip_if_not_installed("DWDLargeR", minimum_version = "0.2.0")

  data <- make_new_method_binary_fixture(n = 28L, p = 45L, seed = 1710L)
  x <- data[setdiff(names(data), c("patient_id", "outcome"))]
  x$gene_10 <- 0
  x$gene_20 <- 0
  adapter <- sb_omics_specialist_binary_learner_adapter("dwd_binary")
  model <- adapter$fit(x, data$outcome)

  expect_identical(
    unname(model$engine$coefficient[c("gene_10", "gene_20")]),
    c(0, 0)
  )
  expect_length(adapter$predict(model, x[1:5, , drop = FALSE]), 5L)
})

test_that("MTLR risk, horizon probability, and lean state match native output", {
  skip_if_not_installed("MTLR", minimum_version = "0.2.2")

  data <- make_new_method_survival_fixture(seed = 1707L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[1:7, , drop = FALSE]
  horizon <- unname(stats::quantile(data$time, 0.35))
  adapter <- sb_mtlr_survival_learner_adapter()
  model <- adapter$fit(x, data$time, data$event)
  internal <- assessment
  names(internal) <- model$internal_features
  native_risk <- -as.numeric(stats::predict(
    model$engine,
    newdata = internal,
    type = "mean_time"
  ))
  native_curve <- stats::predict(
    model$engine,
    newdata = internal,
    type = "prob_times",
    times = horizon
  )
  native_probability <- unname(as.numeric(
    native_curve[nrow(native_curve), -1L]
  ))
  probability <- adapter$predict_survival(
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

  expect_equal(adapter$predict(model, assessment), native_risk, tolerance = 0)
  expect_equal(probability$probability, native_probability, tolerance = 0)
  expect_equal(replay$risk, native_risk, tolerance = 0)
  expect_equal(
    replay$survival_probability,
    native_probability,
    tolerance = 0
  )
  expect_false(any(c("x", "y", "response", "Call") %in%
    names(state$engine)))
  expect_identical(state$training_data_retention, "none")
})

test_that("all four methods execute through CV and final state prediction", {
  binary_data <- make_new_method_binary_fixture(seed = 1708L)
  binary_features <- setdiff(
    names(binary_data),
    c("patient_id", "outcome")
  )
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    binary_features
  )
  dependency <- c(
    brglm2_logistic = "brglm2",
    bayesglm_logistic_cauchy = "arm",
    bst_hinge_binary = "bst",
    hda_binary = "hda",
    dwd_binary = "DWDLargeR",
    ktsp_binary = "switchBox"
  )
  for (learner in names(dependency)) {
    skip_if_not_installed(dependency[[learner]])
    pipeline <- sb_pipeline(learner = learner)
    comparison <- sb_run_binary_cv(
      binary_task,
      pipeline,
      folds = 3L,
      seed = 1709L
    )
    fit <- sb_fit_final(binary_task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      binary_data[1:6, c("patient_id", rev(binary_features)), drop = FALSE],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_true(comparison$complete_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
  }

  skip_if_not_installed("MTLR", minimum_version = "0.2.2")
  survival_data <- make_new_method_survival_fixture(seed = 1710L)
  survival_features <- setdiff(
    names(survival_data),
    c("patient_id", "time", "event")
  )
  horizon <- unname(stats::quantile(survival_data$time, 0.3))
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    survival_features,
    "month",
    prediction_horizon = horizon
  )
  pipeline <- sb_pipeline(learner = "mtlr_survival")
  comparison <- sb_run_cv(survival_task, pipeline, folds = 3L, seed = 1711L)
  fit <- sb_fit_final(survival_task, pipeline)
  prediction <- sb_predict_final(
    unserialize(serialize(fit, NULL)),
    survival_data[
      1:6,
      c("patient_id", rev(survival_features)),
      drop = FALSE
    ],
    "patient_id"
  )

  expect_identical(comparison$status, "success")
  expect_true(comparison$complete_oof)
  expect_true(comparison$complete_survival_probability_oof)
  expect_identical(fit$status, "success")
  expect_identical(prediction$status, "success")
})
})

# ---- model-pamr ----
local({
# ---- algorithm-pamr ----

test_that("PAMR views retain complete original and native shrunken centroids", {
  skip_if_not_installed("pamr")
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in c("pamr_binary", "pamr_survival")) {
    d <- if (id == "pamr_binary") sb_example_binary(160L, seed = 108L) else sb_example_survival(160L, seed = 8612L)
    task <- if (id == "pamr_binary") sb_binary_task(d, "patient_id", "outcome", features) else
      sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    prediction <- sb_predict_final(fit, d, "patient_id")
    e <- if (id == "pamr_binary") unserialize(fit$learner_state$engine_raw) else fit$learner_state$engine
    v <- sb_pamr_diagnostics(fit)
    original <- v$table[v$table$stage == "Before", ]
    after <- v$table[v$table$stage == "After", ]
    expect_equal(original$original_value, as.numeric(e$centroids), tolerance = 0)
    native <- pamr::pamr.predict(e, t(as.matrix(d[features])), threshold = 1, type = "centroid")
    expect_equal(after$original_value, as.numeric(native), tolerance = 0)
    delta <- sweep((e$centroids - e$centroid.overall) / e$sd, 2, e$threshold.scale * e$se.scale, "/")
    soft <- sign(delta) * pmax(abs(delta) - 1, 0)
    shrink <- sweep(soft, 2, e$threshold.scale * e$se.scale, "*")
    expect_equal(after$value, as.numeric(shrink), tolerance = 1e-12)
    expect_equal(after$original_value, as.numeric(e$centroid.overall + shrink * e$sd), tolerance = 1e-12)
    expect_equal(original$feature, rep(features, 2L))
    expect_equal(after$feature, rep(features, 2L))
    expect_equal(nrow(v$table), 12L)
    expect_gt(sum(after$value == 0), 0)
    expect_null(fit$diagnostic_artifacts)
    expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status, "available")
    p <- sb_plot_algorithm_diagnostic(v)
    expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), nrow(v$table))
    expect_equal(p$data$value, v$table$value, tolerance = 0)
    expect_null(p$labels$caption)
    expect_null(p$labels$subtitle)
    expect_identical(sb_predict_final(fit, d, "patient_id"), prediction)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})

# ---- pamr-survival-adapter ----

make_pamr_survival_fixture <- function(n = 90L, p = 12L, seed = 1721L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n)
    colnames(x) <- paste0("gene_", seq_len(p))
    linear_predictor <- 0.8 * x[, 1L] - 0.55 * x[, 2L]
    event_time <- stats::rexp(n, rate = exp(linear_predictor) / 10)
    censor_time <- stats::rexp(n, rate = 1 / 18)
    data.frame(
      patient_id = sprintf("PS%03d", seq_len(n)),
      time = pmin(event_time, censor_time) + 0.05,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    )
  })
}

test_that("PAMR survival canonicalizes equivalent learner identifiers", {
  plain <- sb_pipeline("none", "pamr_survival")
  named <- c(requested = "pamr_survival")
  classed <- structure(
    "pamr_survival",
    class = c("stormbreaker_test_learner", "character")
  )

  expect_identical(sb_pipeline("none", named), plain)
  expect_identical(sb_pipeline("none", classed), plain)
})

test_that("PAMR survival reports its native one-input boundary", {
  skip_if_not_installed("pamr", minimum_version = "1.56")
  data <- make_pamr_survival_fixture(n = 72L, p = 2L, seed = 1720L)

  expect_error(
    sb_pamr_survival_learner_adapter()$fit(
      data["gene_1"], data$time, data$event
    ),
    class = "stormbreaker_pamr_survival_too_few_features"
  )
})

test_that("PAMR survival risk is the native early-survival posterior", {
  skip_if_not_installed("pamr", minimum_version = "1.56")

  data <- make_pamr_survival_fixture(n = 60L, p = 90L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  adapter <- sb_pamr_survival_learner_adapter()
  params <- adapter$params
  set.seed(2721)
  rng_before <- .Random.seed
  model <- adapter$fit(x, data$time, data$event)
  expect_identical(.Random.seed, rng_before)

  native_data <- list(
    x = t(data.matrix(x)),
    survival.time = data$time,
    censoring.status = data$event,
    geneid = features,
    genenames = features
  )
  native <- NULL
  invisible(utils::capture.output(
    native <- pamr::pamr.train(
      native_data,
      threshold = params$threshold,
      scale.sd = params$scale_sd,
      offset.percent = params$offset_percent,
      prior = NULL,
      remove.zeros = params$remove_zero_thresholds,
      sign.contrast = "both",
      ngroup.survival = params$survival_groups
    ),
    type = "output"
  ))
  assessment <- x[1:8, , drop = FALSE]
  native_posterior <- pamr::pamr.predict(
    native,
    t(data.matrix(assessment)),
    threshold = params$threshold,
    type = "posterior"
  )
  risk <- adapter$predict(model, assessment)
  state <- unserialize(serialize(adapter$extract_state(model), NULL))

  expect_equal(risk, unname(native_posterior[, "1"]), tolerance = 0)
  expect_equal(
    adapter$predict_state(state, assessment)$risk,
    risk,
    tolerance = 0
  )
  expect_equal(
    model$audit$survival_cutoff,
    unname(native$cutoffs.survival),
    tolerance = 0
  )
  expect_setequal(
    names(model$engine),
    c(
      "centroids", "centroid.overall", "sd", "threshold.scale",
      "se.scale", "sign.contrast", "prior", "norm.cent"
    )
  )
  expect_identical(state$training_data_retention, "none")
})

test_that("PAMR survival horizon output uses analysis-fold Cox calibration", {
  skip_if_not_installed("pamr", minimum_version = "1.56")

  data <- make_pamr_survival_fixture(seed = 1722L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[1:7, , drop = FALSE]
  horizon <- unname(stats::quantile(data$time, 0.4))
  adapter <- sb_pamr_survival_learner_adapter()
  model <- adapter$fit(x, data$time, data$event)
  analysis_risk <- adapter$predict(model, x)
  calibration <- survival::coxph(
    survival::Surv(data$time, data$event) ~ analysis_risk,
    ties = "breslow",
    singular.ok = FALSE
  )
  coefficient <- unname(as.numeric(stats::coef(calibration)))
  cumulative_hazard <- sb_breslow_baseline_hazard(
    data$time,
    data$event,
    coefficient * analysis_risk,
    horizon
  )
  observed <- adapter$predict_survival(
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

  expect_equal(observed$calibration_coefficient, coefficient, tolerance = 0)
  expect_equal(
    observed$probability,
    sb_survival_probability_from_hazard(
      coefficient * adapter$predict(model, assessment),
      cumulative_hazard
    ),
    tolerance = 0
  )
  expect_equal(replay$risk, adapter$predict(model, assessment), tolerance = 0)
  expect_equal(
    replay$survival_probability,
    observed$probability,
    tolerance = 0
  )
  expect_identical(observed$source, sb_pamr_survival_source())
})

test_that("PAMR survival executes through CV and final prediction", {
  skip_if_not_installed("pamr", minimum_version = "1.56")

  data <- make_pamr_survival_fixture(seed = 1723L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  horizon <- unname(stats::quantile(data$time, 0.35))
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    features,
    "month",
    prediction_horizon = horizon
  )
  pipeline <- sb_pipeline(learner = "pamr_survival")
  comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 1724L)
  fit <- sb_fit_final(task, pipeline)
  prediction <- sb_predict_final(
    unserialize(serialize(fit, NULL)),
    data[1:6, c("patient_id", rev(features)), drop = FALSE],
    "patient_id"
  )

  expect_identical(comparison$status, "success")
  expect_identical(comparison$complete_oof, TRUE)
  expect_identical(comparison$complete_survival_probability_oof, TRUE)
  expect_identical(fit$status, "success")
  expect_identical(prediction$status, "success")
  expect_length(prediction$predictions$risk, 6L)
  expect_length(prediction$predictions$survival_probability, 6L)
  expect_identical(
    fit$provenance$pamr_version,
    as.character(utils::packageVersion("pamr"))
  )
})
})
