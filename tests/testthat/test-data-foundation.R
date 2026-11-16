# ---- adapters ----
local({
test_that("the admitted pipeline resolves explicit component adapters", {
  runtime <- sb_resolve_pipeline_runtime(sb_pipeline())

  expect_s3_class(runtime, "stormbreaker_pipeline_runtime")
  expect_identical(runtime$id, "none__coxph")
  expect_identical(runtime$selector$id, "none")
  expect_type(runtime$selector$fit, "closure")
  expect_type(runtime$selector$apply, "closure")
  expect_identical(runtime$learner$id, "coxph")
  expect_type(runtime$learner$fit, "closure")
  expect_type(runtime$learner$predict, "closure")
  expect_type(runtime$learner$predict_survival, "closure")
  expect_type(runtime$learner$extract_coefficients, "closure")
})

test_that("uniform learner adapter families are distinct and admitted", {
  families <- sb_uniform_learner_adapter_families()
  learner_ids <- unlist(
    lapply(families, function(family) family[[1L]]),
    use.names = FALSE
  )
  admitted <- sb_pipeline_registry(level = "learner")$learner

  expect_identical(anyDuplicated(learner_ids), 0L)
  expect_identical(setdiff(learner_ids, admitted), character())
})

test_that("glmnet logistic resolves an explicit tuning runtime", {
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 101L
  )
  runtime <- sb_resolve_pipeline_runtime(pipeline)

  expect_s3_class(runtime, "stormbreaker_pipeline_runtime")
  expect_identical(runtime$id, pipeline$id)
  expect_identical(runtime$selector$id, "none")
  expect_identical(runtime$learner$id, "glmnet_logistic")
  expect_identical(runtime$learner$params, pipeline$learner_params)
  expect_identical(runtime$learner$requires_tuning, TRUE)
  expect_identical(
    sb_validate_pipeline_runtime(runtime, pipeline),
    runtime
  )
})

test_that("fixed CatBoost resolves an executable tree runtime", {
  pipeline <- sb_pipeline(
    learner = "catboost_binary",
    catboost_iterations = 7L,
    catboost_depth = 3L,
    catboost_seed = 101L
  )
  runtime <- sb_resolve_pipeline_runtime(pipeline)

  expect_s3_class(runtime, "stormbreaker_pipeline_runtime")
  expect_identical(runtime$id, pipeline$id)
  expect_identical(runtime$learner$id, "catboost_binary")
  expect_identical(runtime$learner$params, pipeline$learner_params)
  expect_type(runtime$learner$fit, "closure")
  expect_type(runtime$learner$predict, "closure")
  expect_type(runtime$learner$extract_tree_model, "closure")
  expect_identical(
    sb_validate_pipeline_runtime(runtime, pipeline),
    runtime
  )
})

test_that("the runtime resolver rejects unadmitted component identifiers", {
  pipeline <- sb_pipeline()
  pipeline$selector <- "unknown_selector"
  selector_error <- capture_sb_error(
    sb_resolve_pipeline_runtime(pipeline)
  )
  expect_s3_class(selector_error, "stormbreaker_unresolved_selector_adapter")

  pipeline <- sb_pipeline()
  pipeline$learner <- "unknown_learner"
  learner_error <- capture_sb_error(
    sb_resolve_pipeline_runtime(pipeline)
  )
  expect_s3_class(learner_error, "stormbreaker_unresolved_learner_adapter")
})

test_that("an inconsistent runtime is rejected as a visible fold failure", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  split <- list(
    fold = 1L,
    analysis_rows = 1:10,
    assessment_rows = 11:20,
    analysis_ids = d$patient_id[1:10],
    assessment_ids = d$patient_id[11:20]
  )
  runtime <- sb_resolve_pipeline_runtime(sb_pipeline())
  runtime$learner$id <- "different_learner"

  result <- sb_run_one_fold(
    task, sb_pipeline(), runtime, split, seed = 3L
  )

  expect_identical(result$run$status, "failed")
  expect_identical(
    result$run$error_class,
    "stormbreaker_bad_pipeline_runtime"
  )
  expect_identical(all(is.na(result$predictions$risk)), TRUE)
})

test_that("the none selector locks names without retaining patient data", {
  x <- data.frame(
    gene_b = c(4, 5, 6),
    gene_a = c(1, 2, 3),
    check.names = FALSE
  )
  selector <- sb_none_selector_adapter()
  fitted <- selector$fit(x)

  expect_s3_class(fitted, "stormbreaker_selector_fit")
  expect_named(fitted, c("id", "selected_features"))
  expect_identical(fitted$id, "none")
  expect_identical(fitted$selected_features, names(x))
  expect_identical(any(vapply(fitted, is.data.frame, logical(1))), FALSE)
  expect_identical(any(vapply(fitted, is.matrix, logical(1))), FALSE)
  expect_identical(selector$apply(fitted, x), x)
})

test_that("the none selector enforces its locked predictor schema", {
  x <- data.frame(gene_a = 1:3, gene_b = 4:6)
  selector <- sb_none_selector_adapter()
  fitted <- selector$fit(x)

  reordered_error <- capture_sb_error(
    selector$apply(fitted, x[c("gene_b", "gene_a")])
  )
  expect_s3_class(
    reordered_error,
    "stormbreaker_selector_schema_mismatch"
  )

  missing_error <- capture_sb_error(
    selector$apply(fitted, x["gene_a"])
  )
  expect_s3_class(missing_error, "stormbreaker_selector_schema_mismatch")

  extra <- x
  extra$gene_c <- 7:9
  extra_error <- capture_sb_error(selector$apply(fitted, extra))
  expect_s3_class(extra_error, "stormbreaker_selector_schema_mismatch")
})

test_that("the Cox adapter audits coefficients and reconstructs risk", {
  d <- make_sb_fixture(n = 60)
  analysis_rows <- seq_len(40)
  assessment_rows <- 41:60
  x_analysis <- d[
    analysis_rows,
    c("gene_b", "gene_a"),
    drop = FALSE
  ]
  x_assessment <- d[
    assessment_rows,
    c("gene_b", "gene_a"),
    drop = FALSE
  ]
  learner <- sb_coxph_learner_adapter()
  model <- learner$fit(
    x_analysis,
    d$time[analysis_rows],
    d$event[analysis_rows]
  )
  observed <- learner$predict(model, x_assessment)
  coefficients <- learner$extract_coefficients(model)

  manual <- survival::coxph(
    survival::Surv(time, event) ~ gene_b + gene_a,
    data = d[analysis_rows, , drop = FALSE],
    ties = "efron",
    singular.ok = FALSE,
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
  expected <- unname(stats::predict(
    manual,
    newdata = d[assessment_rows, , drop = FALSE],
    type = "lp",
    reference = "zero"
  ))
  reconstructed <- unname(as.vector(
    as.matrix(x_assessment) %*% coefficients$coefficient
  ))

  expect_named(
    coefficients,
    c("feature", "coefficient", "nonzero", "direction"),
    ignore.order = FALSE
  )
  expect_identical(coefficients$feature, c("gene_b", "gene_a"))
  expect_equal(
    coefficients$coefficient,
    unname(stats::coef(model)[coefficients$feature]),
    tolerance = 0
  )
  expect_identical(coefficients$nonzero, coefficients$coefficient != 0)
  expect_identical(
    coefficients$direction,
    ifelse(
      coefficients$coefficient > 0,
      "higher_hazard",
      ifelse(coefficients$coefficient < 0, "lower_hazard", "zero")
    )
  )
  expect_equal(observed, expected, tolerance = 1e-12)
  expect_equal(observed, reconstructed, tolerance = 1e-12)

  schema_error <- capture_sb_error(
    learner$predict(model, x_assessment[c("gene_a", "gene_b")])
  )
  expect_s3_class(schema_error, "stormbreaker_learner_schema_mismatch")
})

test_that("ordinary Cox survival probabilities retain Efron baseline semantics", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  analysis_rows <- seq_len(40L)
  assessment_rows <- 41:60
  predictors <- c("gene_b", "gene_a")
  x_analysis <- d[analysis_rows, predictors, drop = FALSE]
  x_assessment <- d[assessment_rows, predictors, drop = FALSE]
  horizon <- 12
  learner <- sb_coxph_learner_adapter()
  model <- learner$fit(
    x_analysis,
    d$time[analysis_rows],
    d$event[analysis_rows]
  )

  observed <- learner$predict_survival(
    model,
    x_analysis,
    d$time[analysis_rows],
    d$event[analysis_rows],
    x_assessment,
    horizon
  )
  basehaz <- survival::basehaz(model, centered = FALSE)
  expected_hazard <- tail(basehaz$hazard[basehaz$time <= horizon], 1L)
  expected_survival <- drop(summary(
    survival::survfit(model, newdata = x_assessment, se.fit = FALSE),
    times = horizon,
    extend = TRUE
  )$surv)
  risk <- learner$predict(model, x_assessment)
  reconstructed <- exp(-expected_hazard * exp(risk))

  expect_named(
    observed,
    c(
      "prediction_horizon", "cumulative_baseline_hazard",
      "baseline_survival", "source", "ties", "probability"
    ),
    ignore.order = FALSE
  )
  expect_identical(observed$prediction_horizon, horizon)
  expect_identical(observed$source, "coxph_basehaz_zero_reference")
  expect_identical(observed$ties, "efron")
  expect_equal(
    observed$cumulative_baseline_hazard,
    expected_hazard,
    tolerance = 1e-14
  )
  expect_equal(
    observed$baseline_survival,
    exp(-expected_hazard),
    tolerance = 1e-14
  )
  expect_equal(
    observed$probability,
    unname(expected_survival),
    tolerance = 1e-12
  )
  expect_equal(observed$probability, reconstructed, tolerance = 1e-12)
})

test_that("Cox coefficient audits preserve non-syntactic feature names", {
  d <- make_sb_fixture(n = 40L, seed = 20260819L)
  x_analysis <- d[1:30, c("gene_a", "gene_b"), drop = FALSE]
  x_assessment <- d[31:40, c("gene_a", "gene_b"), drop = FALSE]
  names(x_analysis) <- c("gene a", "gene-b")
  names(x_assessment) <- names(x_analysis)
  learner <- sb_coxph_learner_adapter()

  model <- learner$fit(
    x_analysis,
    d$time[1:30],
    d$event[1:30]
  )
  coefficients <- learner$extract_coefficients(model)
  observed <- learner$predict(model, x_assessment)
  survival_prediction <- learner$predict_survival(
    model,
    x_analysis,
    d$time[1:30],
    d$event[1:30],
    x_assessment,
    prediction_horizon = 12
  )
  reconstructed <- unname(as.vector(
    as.matrix(x_assessment[coefficients$feature]) %*%
      coefficients$coefficient
  ))

  expect_identical(coefficients$feature, c("gene a", "gene-b"))
  expect_true(all(is.finite(coefficients$coefficient)))
  expect_equal(reconstructed, observed, tolerance = 1e-12)
  expect_length(survival_prediction$probability, nrow(x_assessment))
  expect_true(all(is.finite(survival_prediction$probability)))
})

test_that("the unpenalized Cox adapter rejects p greater than or equal to n", {
  x <- data.frame(
    gene_a = c(0.1, 0.4, 0.2, 0.8),
    gene_b = c(1.2, 0.9, 1.5, 1.1),
    gene_c = c(2.1, 2.4, 2.3, 2.8),
    gene_d = c(3.2, 3.4, 3.1, 3.9)
  )

  error <- capture_sb_error(
    sb_fit_coxph(
      x,
      time = c(1, 2, 3, 4),
      event = c(1, 0, 1, 0)
    )
  )

  expect_s3_class(error, "stormbreaker_high_dimensional_cox")
  expect_match(conditionMessage(error), "4 predictors with 4 patients")
  expect_match(conditionMessage(error), "fold-local selector or glmnet Cox")
})

test_that("the fixed glmnet Cox adapter matches a direct path fit", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60)
  analysis_rows <- seq_len(40)
  assessment_rows <- 41:60
  x_analysis <- d[
    analysis_rows,
    c("gene_a", "gene_b", "gene_c"),
    drop = FALSE
  ]
  x_assessment <- d[
    assessment_rows,
    c("gene_a", "gene_b", "gene_c"),
    drop = FALSE
  ]
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = 1,
    lambda = 0.05,
    standardize = TRUE
  )
  model <- learner$fit(
    x_analysis,
    d$time[analysis_rows],
    d$event[analysis_rows]
  )
  observed <- learner$predict(model, x_assessment)

  manual <- glmnet::glmnet(
    x = as.matrix(x_analysis),
    y = survival::Surv(
      d$time[analysis_rows],
      d$event[analysis_rows]
    ),
    family = "cox",
    alpha = 1,
    standardize = TRUE
  )
  expected <- unname(drop(stats::predict(
    manual,
    newx = as.matrix(x_assessment),
    s = 0.05,
    type = "link"
  )))

  expect_identical(learner$id, "glmnet_cox")
  expect_equal(observed, expected, tolerance = 1e-12)
  expect_length(observed, nrow(x_assessment))
  expect_identical(all(is.finite(observed)), TRUE)

  schema_error <- capture_sb_error(
    learner$predict(model, x_assessment[c("gene_b", "gene_a", "gene_c")])
  )
  expect_s3_class(
    schema_error,
    "stormbreaker_learner_schema_mismatch"
  )
})

test_that("glmnet Cox survival probabilities use fixed-beta Breslow ties", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  d$time <- rep(seq_len(10L), each = 6L)
  d$event <- rep(c(1L, 0L, 1L, 0L, 1L, 0L), 10L)
  analysis_rows <- seq_len(45L)
  assessment_rows <- 46:60
  predictors <- c("gene_a", "gene_b", "gene_c")
  x_analysis <- d[analysis_rows, predictors, drop = FALSE]
  x_assessment <- d[assessment_rows, predictors, drop = FALSE]
  time <- d$time[analysis_rows]
  event <- d$event[analysis_rows]
  horizon <- 6
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = 1,
    lambda = 0.05,
    standardize = TRUE
  )
  model <- learner$fit(x_analysis, time, event)

  observed <- learner$predict_survival(
    model,
    x_analysis,
    time,
    event,
    x_assessment,
    horizon
  )
  analysis_risk <- learner$predict(model, x_analysis)
  oracle_data <- data.frame(time = time, event = event, eta = analysis_risk)
  oracle <- survival::coxph(
    survival::Surv(time, event) ~ eta,
    data = oracle_data,
    ties = "breslow",
    init = 1,
    iter = 0,
    singular.ok = FALSE
  )
  oracle_basehaz <- survival::basehaz(oracle, centered = FALSE)
  expected_hazard <- tail(
    oracle_basehaz$hazard[oracle_basehaz$time <= horizon],
    1L
  )
  assessment_risk <- learner$predict(model, x_assessment)
  expected_probability <- exp(-expected_hazard * exp(assessment_risk))
  efron_oracle <- survival::coxph(
    survival::Surv(time, event) ~ eta,
    data = oracle_data,
    ties = "efron",
    init = 1,
    iter = 0,
    singular.ok = FALSE
  )
  efron_basehaz <- survival::basehaz(efron_oracle, centered = FALSE)
  efron_hazard <- tail(
    efron_basehaz$hazard[efron_basehaz$time <= horizon],
    1L
  )

  expect_identical(observed$source, "glmnet_fixed_beta_breslow")
  expect_identical(observed$ties, "breslow")
  expect_equal(
    observed$cumulative_baseline_hazard,
    expected_hazard,
    tolerance = 1e-12
  )
  expect_equal(
    observed$probability,
    expected_probability,
    tolerance = 1e-12
  )
  expect_gt(abs(expected_hazard - efron_hazard), 1e-8)
  expect_equal(
    learner$extract_coefficients(model)$coefficient,
    unname(as.matrix(stats::coef(model, s = 0.05))[, 1L]),
    tolerance = 0
  )
})

test_that("glmnet Cox reports every fitted coefficient at its lambda", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60)
  predictors <- c("gene_a", "gene_b", "gene_c")
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = 1,
    lambda = 0.05,
    standardize = TRUE
  )
  model <- learner$fit(
    d[predictors],
    d$time,
    d$event
  )
  observed <- learner$extract_coefficients(model)
  expected <- drop(as.matrix(stats::coef(model, s = 0.05)))

  expect_named(
    observed,
    c("feature", "coefficient", "nonzero", "direction"),
    ignore.order = FALSE
  )
  expect_identical(observed$feature, predictors)
  expect_equal(observed$coefficient, unname(expected), tolerance = 0)
  expect_identical(observed$nonzero, unname(expected != 0))
  expect_identical(
    observed$direction,
    unname(ifelse(expected > 0, "higher_hazard", ifelse(
      expected < 0,
      "lower_hazard",
      "zero"
    )))
  )
  expect_identical(any(observed$nonzero), TRUE)
  expect_identical(any(!observed$nonzero), TRUE)
})

test_that("the fixed glmnet Cox adapter rejects invalid settings", {
  alpha_error <- capture_sb_error(
    sb_glmnet_cox_learner_adapter(alpha = 1.2, lambda = 0.05)
  )
  expect_s3_class(alpha_error, "stormbreaker_bad_glmnet_alpha")

  lambda_error <- capture_sb_error(
    sb_glmnet_cox_learner_adapter(alpha = 1, lambda = 0)
  )
  expect_s3_class(lambda_error, "stormbreaker_bad_glmnet_lambda")

  standardize_error <- capture_sb_error(
    sb_glmnet_cox_learner_adapter(
      alpha = 1,
      lambda = 0.05,
      standardize = NA
    )
  )
  expect_s3_class(
    standardize_error,
    "stormbreaker_bad_glmnet_standardize"
  )
})

test_that("the fixed glmnet Cox adapter reports a missing dependency", {
  local_mocked_bindings(
    sb_glmnet_available = function() FALSE
  )
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = 1,
    lambda = 0.05,
    standardize = TRUE
  )
  d <- make_sb_fixture(n = 20)

  error <- capture_sb_error(learner$fit(
    d[c("gene_a", "gene_b")],
    d$time,
    d$event
  ))

  expect_s3_class(error, "stormbreaker_missing_glmnet")
})

test_that("the fixed glmnet Cox adapter rejects non-positive time", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60)
  d$time[1] <- 0
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = 1,
    lambda = 0.05,
    standardize = TRUE
  )

  error <- capture_sb_error(learner$fit(
    d[c("gene_a", "gene_b", "gene_c")],
    d$time,
    d$event
  ))

  expect_s3_class(error, "stormbreaker_glmnet_nonpositive_time")
})

test_that("glmnet Cox refuses lambda values outside its fitted path", {
  skip_if_not_installed("glmnet")

  x <- data.frame(
    g1 = c(
      -1.941518376040617833, 1.181180891438521785,
      1.859910861543160987, 1.074012260281989928,
      -0.027346966413741798
    ),
    g2 = c(
      -0.033330342275919272, -1.516067623151099220,
      0.790385343888228009, -0.210734180691775658,
      -0.656742925361075369
    )
  )
  time <- c(
    2.34298038283798737, 0.62733672885224223,
    1.22186382110686331, 0.60436541587114334,
    1.96968345521929300
  )
  event <- c(1, 1, 1, 0, 0)
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = 1,
    lambda = 4.0337729269749e-12,
    standardize = TRUE
  )

  error <- capture_sb_error(learner$fit(x, time, event))

  expect_s3_class(error, "stormbreaker_glmnet_lambda_outside_path")
})

test_that("assessment endpoints cannot change fixed-split adapter risks", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )
  split <- list(
    fold = 1L,
    analysis_rows = seq_len(48),
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[seq_len(48)],
    assessment_ids = d$patient_id[49:60]
  )
  runtime <- sb_resolve_pipeline_runtime(sb_pipeline())
  original <- sb_run_one_fold(
    task, sb_pipeline(), runtime, split, seed = 17L
  )

  changed <- d
  changed$time[split$assessment_rows] <-
    changed$time[split$assessment_rows] + 1000
  changed$event[split$assessment_rows] <-
    1L - changed$event[split$assessment_rows]
  changed_task <- sb_surv_task(
    changed, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )
  perturbed <- sb_run_one_fold(
    changed_task, sb_pipeline(), runtime, split, seed = 17L
  )

  expect_equal(
    original$predictions$risk,
    perturbed$predictions$risk,
    tolerance = 0
  )
  expect_identical(original$coefficients$status, "success")
  expect_identical(original$coefficients, perturbed$coefficients)
})

test_that("each adapter operation receives only its permitted fold data", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )
  split <- list(
    fold = 1L,
    analysis_rows = 1:12,
    assessment_rows = 13:20,
    analysis_ids = d$patient_id[1:12],
    assessment_ids = d$patient_id[13:20]
  )
  received <- new.env(parent = emptyenv())
  runtime <- structure(list(
    id = "none__coxph",
    selector = list(
      id = "none",
      params = list(),
      needs_outcome = FALSE,
      fit = function(x) {
        received$selector_fit <- x
        list(features = names(x))
      },
      apply = function(fitted, x) {
        if (is.null(received$selector_apply)) {
          received$selector_apply <- list()
        }
        received$selector_apply[[length(received$selector_apply) + 1L]] <- x
        x[fitted$features]
      }
    ),
    learner = list(
      id = "coxph",
      fit = function(x, time, event) {
        received$learner_fit <- list(x = x, time = time, event = event)
        list(features = names(x))
      },
      predict = function(model, x) {
        received$learner_predict <- x
        rep.int(0, nrow(x))
      },
      predict_survival = function(model, x_analysis, time, event,
                                  x_assessment, prediction_horizon) {
        stop("not requested by this fixture")
      },
      extract_coefficients = function(model) {
        data.frame(
          feature = model$features,
          coefficient = rep(0, length(model$features)),
          nonzero = rep(FALSE, length(model$features)),
          direction = rep("zero", length(model$features)),
          stringsAsFactors = FALSE
        )
      },
      extract_state = function(model, context = list()) {
        list(features = model$features)
      },
      predict_state = function(state, x) rep.int(0, nrow(x))
    )
  ), class = "stormbreaker_pipeline_runtime")

  result <- sb_run_one_fold(
    task, sb_pipeline(), runtime, split, seed = 23L
  )

  expected_analysis <- d[1:12, c("gene_a", "gene_b"), drop = FALSE]
  expected_assessment <- d[13:20, c("gene_a", "gene_b"), drop = FALSE]
  expect_identical(received$selector_fit, expected_analysis)
  expect_identical(received$selector_apply[[1]], expected_analysis)
  expect_identical(received$selector_apply[[2]], expected_assessment)
  expect_identical(received$learner_fit$x, expected_analysis)
  expect_identical(received$learner_fit$time, d$time[1:12])
  expect_identical(received$learner_fit$event, d$event[1:12])
  expect_identical(received$learner_predict, expected_assessment)
  expect_identical(result$run$status, "success")
})

test_that("an adapter failure is recorded without invoking a fallback", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  split <- list(
    fold = 1L,
    analysis_rows = seq_len(10),
    assessment_rows = 11:20,
    analysis_ids = d$patient_id[seq_len(10)],
    assessment_ids = d$patient_id[11:20]
  )
  learner_fit_calls <- 0L
  runtime <- structure(list(
    id = "none__coxph",
    selector = list(
      id = "none",
      params = list(),
      needs_outcome = FALSE,
      fit = function(x) {
        sb_abort(
          "Synthetic selector failure.",
          "stormbreaker_test_selector_failure"
        )
      },
      apply = function(fitted, x) x
    ),
    learner = list(
      id = "coxph",
      fit = function(x, time, event) {
        learner_fit_calls <<- learner_fit_calls + 1L
        NULL
      },
      predict = function(model, x) rep.int(0, nrow(x)),
      predict_survival = function(model, x_analysis, time, event,
                                  x_assessment, prediction_horizon) {
        stop("not requested by this fixture")
      },
      extract_coefficients = function(model) NULL,
      extract_state = function(model, context = list()) list(),
      predict_state = function(state, x) rep.int(0, nrow(x))
    )
  ), class = "stormbreaker_pipeline_runtime")

  result <- sb_run_one_fold(
    task, sb_pipeline(), runtime, split, seed = 19L
  )

  expect_identical(result$run$status, "failed")
  expect_identical(
    result$run$error_class,
    "stormbreaker_test_selector_failure"
  )
  expect_identical(result$run$learner_used, NA_character_)
  expect_identical(learner_fit_calls, 0L)
  expect_identical(all(is.na(result$predictions$risk)), TRUE)
})

test_that("a prediction failure retains an already fitted Cox audit", {
  d <- make_sb_fixture(n = 30L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month"
  )
  split <- list(
    fold = 1L,
    analysis_rows = 1:20,
    assessment_rows = 21:30,
    analysis_ids = d$patient_id[1:20],
    assessment_ids = d$patient_id[21:30]
  )
  pipeline <- sb_pipeline()
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  runtime$learner$predict <- function(model, x) {
    sb_abort(
      "Synthetic Cox prediction failure.",
      "stormbreaker_test_cox_prediction_failure"
    )
  }

  result <- sb_run_one_fold(task, pipeline, runtime, split, seed = 1L)

  expect_identical(result$run$status, "failed")
  expect_identical(
    result$run$error_class,
    "stormbreaker_test_cox_prediction_failure"
  )
  expect_true(all(is.na(result$predictions$risk)))
  expect_identical(result$coefficients$status, "success")
  expect_identical(result$coefficients$values$feature, task$predictors)
  expect_true(all(is.finite(result$coefficients$values$coefficient)))
})

test_that("univariable Cox Top-K matches independent manual scores", {
  d <- make_sb_fixture(n = 60)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  selector <- sb_univ_cox_top_k_selector_adapter(n_features = 2L)
  fitted <- selector$fit(x, d$time, d$event)

  manual_p <- vapply(names(x), function(feature) {
    fit <- survival::coxph(
      stats::reformulate(
        feature,
        response = "survival::Surv(time, event)"
      ),
      data = d,
      ties = "efron",
      singular.ok = FALSE
    )
    unname(summary(fit)$coefficients[1, "Pr(>|z|)"])
  }, numeric(1))
  expected <- names(x)[order(manual_p, seq_along(manual_p))][1:2]

  expect_s3_class(fitted, "stormbreaker_selector_fit")
  expect_identical(fitted$id, "univ_cox_top_k")
  expect_identical(fitted$selected_features, expected)
  expect_named(
    fitted$scores,
    c(
      "feature", "original_order", "coefficient", "standard_error",
      "wald_z", "p_value", "status", "error_class", "error_message",
      "rank", "selected"
    ),
    ignore.order = FALSE
  )
  expect_equal(
    fitted$scores$p_value,
    unname(manual_p[match(fitted$scores$feature, names(manual_p))]),
    tolerance = 1e-12
  )
  expect_identical(sum(fitted$scores$selected), 2L)
  expect_identical(selector$apply(fitted, x), x[expected])
})

test_that("Top-K breaks exact score ties by original predictor order", {
  d <- make_sb_fixture(n = 60)
  x <- data.frame(
    first = d$gene_a,
    second = d$gene_a,
    third = d$gene_b
  )
  selector <- sb_univ_cox_top_k_selector_adapter(n_features = 1L)
  fitted <- selector$fit(x, d$time, d$event)

  expect_identical(fitted$selected_features, "first")
})

test_that("Top-K records invalid features and fails when K is infeasible", {
  d <- make_sb_fixture(n = 60)
  x <- data.frame(valid = d$gene_a, constant = 1)
  selector <- sb_univ_cox_top_k_selector_adapter(n_features = 1L)
  fitted <- selector$fit(x, d$time, d$event)

  expect_identical(fitted$selected_features, "valid")
  expect_identical(
    fitted$scores$status[fitted$scores$feature == "constant"],
    "failed"
  )
  expect_identical(
    is.na(fitted$scores$p_value[fitted$scores$feature == "constant"]),
    TRUE
  )

  infeasible <- sb_univ_cox_top_k_selector_adapter(n_features = 2L)
  error <- capture_sb_error(infeasible$fit(x, d$time, d$event))
  expect_s3_class(error, "stormbreaker_insufficient_valid_features")
})

test_that("Top-K does not let a single patient define a feature score", {
  d <- make_sb_fixture(n = 60)
  x <- data.frame(
    stable = d$gene_a,
    singleton = c(1, rep(0, 59))
  )
  selector <- sb_univ_cox_top_k_selector_adapter(n_features = 1L)

  fitted <- selector$fit(x, d$time, d$event)
  singleton_score <- fitted$scores[
    fitted$scores$feature == "singleton",
    ,
    drop = FALSE
  ]

  expect_identical(fitted$selected_features, "stable")
  expect_identical(singleton_score$status, "failed")
  expect_identical(
    singleton_score$error_class,
    "stormbreaker_selector_insufficient_feature_support"
  )
  expect_match(singleton_score$error_message, "single patient")
})

test_that("Top-K rejects K larger than the candidate feature count", {
  d <- make_sb_fixture(n = 20)
  x <- d[c("gene_a", "gene_b")]
  selector <- sb_univ_cox_top_k_selector_adapter(n_features = 3L)

  error <- capture_sb_error(selector$fit(x, d$time, d$event))

  expect_s3_class(error, "stormbreaker_insufficient_candidate_features")
})

test_that("a Top-K selector failure never reaches learner fitting", {
  d <- make_sb_fixture(n = 20)
  d$constant_a <- 1
  d$constant_b <- 1
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("constant_a", "constant_b"), "month"
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    selector_n_features = 1L
  )
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  learner_fit_calls <- 0L
  runtime$learner$fit <- function(x, time, event) {
    learner_fit_calls <<- learner_fit_calls + 1L
    NULL
  }
  split <- list(
    fold = 1L,
    analysis_rows = 1:12,
    assessment_rows = 13:20,
    analysis_ids = d$patient_id[1:12],
    assessment_ids = d$patient_id[13:20]
  )

  result <- sb_run_one_fold(task, pipeline, runtime, split, seed = 7L)

  expect_identical(result$run$status, "failed")
  expect_identical(
    result$run$error_class,
    "stormbreaker_insufficient_valid_features"
  )
  expect_identical(learner_fit_calls, 0L)
  expect_identical(all(is.na(result$predictions$risk)), TRUE)
})

test_that("Top-K apply enforces the fitted source schema", {
  d <- make_sb_fixture(n = 60)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  selector <- sb_univ_cox_top_k_selector_adapter(n_features = 2L)
  fitted <- selector$fit(x, d$time, d$event)

  reordered_error <- capture_sb_error(
    selector$apply(fitted, x[c("gene_b", "gene_a", "gene_c")])
  )
  expect_s3_class(
    reordered_error,
    "stormbreaker_selector_schema_mismatch"
  )
})

test_that("Top-K receives only analysis endpoints in fold execution", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    selector_n_features = 2L
  )
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  captured <- new.env(parent = emptyenv())
  original_fit <- runtime$selector$fit
  runtime$selector$fit <- function(x, time, event) {
    captured$x <- x
    captured$time <- time
    captured$event <- event
    original_fit(x, time, event)
  }

  result <- sb_run_one_fold(task, pipeline, runtime, split, seed = 17L)

  expect_identical(captured$x, d[1:48, task$predictors, drop = FALSE])
  expect_identical(captured$time, d$time[1:48])
  expect_identical(captured$event, d$event[1:48])
  expect_identical(result$run$status, "success")
  expect_length(result$selection$selected_features, 2L)
})

test_that("assessment endpoints cannot change fixed-split Top-K selection", {
  d <- make_sb_fixture(n = 60)
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    selector_n_features = 2L
  )
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  make_task <- function(data) {
    sb_surv_task(
      data, "patient_id", "time", "event",
      c("gene_a", "gene_b", "gene_c"), "month"
    )
  }
  original <- sb_run_one_fold(
    make_task(d), pipeline, sb_resolve_pipeline_runtime(pipeline), split,
    seed = 17L
  )
  changed <- d
  changed$time[49:60] <- changed$time[49:60] + 1000
  changed$event[49:60] <- 1L - changed$event[49:60]
  perturbed <- sb_run_one_fold(
    make_task(changed), pipeline, sb_resolve_pipeline_runtime(pipeline),
    split, seed = 17L
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
})

# ---- constant-sum ----
local({
# ---- constant-sum-final-fit ----

test_that("public final fits replay constrained predictions and original-unit views", {
  skip_if_not_installed("brglm2")
  withr::local_seed(410)
  x <- as.data.frame(matrix(stats::rnorm(320), 80, 4))
  names(x) <- paste0("gene", 1:4)
  x[] <- sweep(as.matrix(x), 1L, rowMeans(x), "-")
  d <- x
  d$id <- paste0("p", seq_len(nrow(d)))
  d$outcome <- rep(0:1, 40)
  task <- sb_binary_task(d, "id", "outcome", names(x))
  for (id in c("brglm2_logistic", "lda_binary")) {
    fit <- sb_fit_final(task, sb_model_pipelines("binary", id)[[1L]])
    expect_identical(fit$status, "success")
    prediction <- sb_predict_final(fit, d[c("id", names(x))], "id")
    expect_identical(prediction$status, "success")
    replay <- sb_predict_final(unserialize(serialize(fit, NULL)), d[c("id", names(x))], "id")
    expect_identical(replay$predictions, prediction$predictions)
    coefficients <- sb_model_coefficients(fit)
    expect_identical(coefficients$table$feature, names(x))
    expect_equal(stats::plogis(coefficients$intercept +
      as.vector(as.matrix(x) %*% coefficients$table$coefficient)),
      prediction$predictions$probability, tolerance = 1e-10)
    if (id == "lda_binary") {
      diagnostic <- sb_discriminant_diagnostics(fit)
      expect_equal(nrow(diagnostic$table), 2L * ncol(x))
      expect_equal(unname(diagnostic$parameters$means["0", ]),
        unname(colMeans(x[d$outcome == 0L, ])), tolerance = 1e-12)
    }
  }
})

# ---- constant-sum-representation ----

constant_sum_fixture <- function() {
  withr::with_seed(41, {
    x <- as.data.frame(matrix(stats::rnorm(320), 80, 4))
    names(x) <- paste0("g", 1:4)
    x[] <- sweep(as.matrix(x), 1L, rowMeans(x), "-")
    list(x = x, outcome = rep(0:1, 40))
  })
}

test_that("formal-run constant-sum coordinates preserve native predictions", {
  for (package in c("brglm2", "MASS", "mda", "hda", "klaR")) skip_if_not_installed(package)
  d <- constant_sum_fixture()
  for (id in c("brglm2_logistic", "lda_binary", "fda_binary", "mda_binary", "hda_binary", "loclda_binary")) {
    pipeline <- sb_model_pipelines("binary", id)[[1L]]
    adapter <- sb_resolve_learner_adapter(id, pipeline$learner_params)
    model <- adapter$fit(d$x, d$outcome)
    state <- adapter$extract_state(model)
    expect_s3_class(state, "stormbreaker_learner_state")
    expect_identical(state$features, names(d$x), info = id)
    expect_identical(state$projection$features, names(d$x), info = id)
    prediction <- adapter$predict(model, d$x)
    replay <- adapter$predict_state(state, d$x)$probability
    expect_equal(replay, prediction, tolerance = 1e-12, info = id)
    family <- Filter(function(f) id %in% f[[1L]], sb_uniform_learner_adapter_families())[[1L]]
    native <- family[[2L]](id, state$inner_params)
    native_x <- as.data.frame(sweep(as.matrix(d$x), 2L, state$projection$center, "-") %*% state$projection$basis)
    names(native_x) <- paste0("contrast_", seq_len(ncol(native_x)))
    reference <- native$fit(native_x, d$outcome)
    expect_equal(prediction, native$predict(reference, native_x), tolerance = 1e-12, info = id)
    old <- if (id == "brglm2_logistic") {
      list(variant = "brglm2_constant_sum_v1", projection = state$projection, inner = state$inner)
    } else {
      list(round2_representation = "constant_sum_orthonormal_v1", projection = state$projection,
        inner_params = state$inner_params, inner = state$inner)
    }
    expect_equal(adapter$predict_state(old, d$x)$probability, prediction, tolerance = 1e-12, info = id)
    shifted <- d$x
    shifted$g1 <- shifted$g1 + 0.1
    error <- tryCatch(adapter$predict_state(state, shifted), error = identity)
    expect_s3_class(error, "stormbreaker_constant_sum_mismatch")
  }
})

test_that("ordinary inputs retain native learner states", {
  skip_if_not_installed("brglm2")
  d <- constant_sum_fixture()
  d$x$g1 <- seq_len(nrow(d$x)) / 100
  p <- sb_model_pipelines("binary", "brglm2_logistic")[[1L]]
  adapter <- sb_resolve_learner_adapter(p$learner, p$learner_params)
  model <- adapter$fit(d$x, d$outcome)
  expect_s3_class(adapter$extract_state(model), "stormbreaker_omics_specialist_binary_state")
  expect_null(model$projection)
})

test_that("constant-sum linear parameters stay in original gene units", {
  skip_if_not_installed("brglm2")
  d <- constant_sum_fixture()
  for (id in c("brglm2_logistic", "lda_binary")) {
    pipeline <- sb_model_pipelines("binary", id)[[1L]]
    adapter <- sb_resolve_learner_adapter(id, pipeline$learner_params)
    model <- adapter$fit(d$x, d$outcome)
    state <- adapter$extract_state(model)
    parameters <- if (id == "lda_binary") sb_discriminant_parameters(state) else sb_extended_logistic_parameters(state)
    expect_identical(names(parameters$coefficient), names(d$x))
    expect_equal(stats::plogis(parameters$intercept + as.vector(as.matrix(d$x) %*% parameters$coefficient)),
      adapter$predict(model, d$x), tolerance = 1e-10)
    wrong_view <- tryCatch({
      if (id == "lda_binary") sb_extended_logistic_parameters(state) else sb_discriminant_parameters(state)
    }, error = identity)
    expect_s3_class(wrong_view, if (id == "lda_binary") "stormbreaker_unsupported_coefficients" else
      "stormbreaker_unsupported_diagnostic")
  }
})

# ---- state-parameter-equivalence ----

test_that("flexsurv states accept equivalent parameter representations", {
  skip_if_not_installed("flexsurv")
  withr::local_seed(741L)
  x <- data.frame(x1 = rnorm(60L), x2 = rnorm(60L))
  time <- rgamma(60L, shape = 2, rate = exp(0.2 * x$x1) / 10)
  event <- rep(c(1L, 1L, 0L), 20L)
  for (kind in c("parametric", "spline")) {
    id <- if (kind == "parametric") "flexsurv_gamma" else "flexsurvspline_hazard_k1"
    stem <- if (kind == "parametric") "flexsurv_parametric" else "flexsurvspline"
    params <- get(paste0("sb_", stem, "_fixed_params"))(id)
    model <- get(paste0("sb_fit_", stem, "_fixed"))(x, time, event, id, params)
    state <- get(paste0("sb_extract_", stem, "_state"))(model)
    predict <- get(paste0("sb_predict_", stem, "_state"))
    expected <- predict(state, x)
    changed <- state
    changed$params <- rev(state$params)
    expect_identical(predict(changed, x), expected)
    if (kind == "spline") {
      changed$params$k <- as.numeric(changed$params$k)
      expect_identical(predict(changed, x), expected)
      changed$params$k <- 2
    } else {
      changed$params$distribution <- "gompertz"
    }
    error <- tryCatch(predict(changed, x), error = identity)
    expect_s3_class(error, "stormbreaker_bad_learner_state")
  }
})

test_that("mboost states compare iteration values rather than storage types", {
  skip_if_not_installed("mboost")
  withr::local_seed(742L)
  x <- data.frame(x1 = rnorm(40L), x2 = rnorm(40L))
  time <- rgamma(40L, shape = 2, rate = exp(0.2 * x$x1) / 10)
  event <- rep(c(1L, 1L, 1L, 0L), 10L)
  for (kind in c("aft", "ranking")) {
    id <- if (kind == "aft") "glmboost_lognormal_aft" else "glmboost_gehan_aft"
    stem <- paste0("mboost_", kind)
    params <- get(paste0("sb_", stem, "_fixed_params"))(id)
    params$mstop <- 5L
    model <- get(paste0("sb_fit_", stem, "_fixed"))(x, time, event, id, params)
    state <- get(paste0("sb_extract_", stem, "_state"))(model)
    predict <- get(paste0("sb_predict_", stem, "_state"))
    expected <- predict(state, x)
    changed <- state
    changed$params$mstop <- 5
    expect_identical(predict(changed, x), expected)
    for (wrong in list(6, 5.01, TRUE, NA_real_)) {
      changed$params$mstop <- wrong
      error <- tryCatch(predict(changed, x), error = identity)
      expect_s3_class(error, "stormbreaker_bad_learner_state")
    }
  }
})
})

# ---- core-conditions ----
local({
test_that("retired learners have no registry, constructor or adapter route", {
  retired <- sb_retired_learner_ids()
  registry <- sb_pipeline_registry()
  expect_identical(intersect(registry$learner, retired), character())
  expect_length(unique(registry$learner), 176L)
  for (id in retired) {
    expect_error(sb_pipeline(learner = id), class = "stormbreaker_retired_learner")
    expect_error(sb_pipeline("ordered_chain", id), class = "stormbreaker_retired_learner")
    expect_error(sb_resolve_learner_adapter(id), class = "stormbreaker_retired_learner")
    p <- sb_pipeline(learner = "logistic")
    p$learner <- id
    expect_error(sb_validate_pipeline_contract(p, "binary"),
      class = "stormbreaker_retired_learner")
  }
  expect_identical(sb_mgcv_gam_learner_ids(), character())
  expect_identical(sb_slope_binary_learner_ids(), character())
  expect_identical(sb_slope_randomglm_binary_learner_ids(), "random_glm_binary")
})

test_that("removed native fitting implementations cannot be invoked", {
  expect_error(sb_fit_mgcv_gam(NULL, NULL, "gam_cox_reml", NULL),
    class = "stormbreaker_retired_learner")
  expect_error(sb_fit_mgcv_gam(NULL, NULL, "gam_binary_reml", NULL),
    class = "stormbreaker_retired_learner")
  expect_error(sb_fit_slope_binary(NULL, NULL), class = "stormbreaker_retired_learner")
})

test_that("randomGLM and GAM boosting remain constructible", {
  for (id in c("random_glm_binary", "gamboost_binary", "gamboost_cox",
    "gamboost_gehan_aft", "gamboost_loglogistic_aft", "gamboost_lognormal_aft",
    "gamboost_weibull_aft", "gamboost_unoc")) {
    expect_identical(sb_pipeline(learner = id)$learner, id)
  }
  expect_identical(sb_slope_randomglm_binary_params("random_glm_binary")$n_bags, 100L)
})
})

# ---- feature-specification ----
local({
# ---- feature-mapping ----

test_that("explicit feature mapping is outcome-blind and ordered", {
  new_data <- data.frame(
    patient_id = c("p1", "p2", "p3"),
    outcome = c(0L, 1L, 0L),
    ENSG_B = c(2.0, 2.5, 3.0),
    unrelated = c(10, 20, 30),
    ENSG_A = c(-1.0, 0.0, 1.0),
    check.names = FALSE
  )
  feature_map <- data.frame(
    source_feature = c("ENSG_B", "ENSG_A"),
    target_feature = c("gene_b", "gene_a"),
    annotation_release = c("declared_release", "declared_release"),
    stringsAsFactors = FALSE
  )

  observed <- sb_apply_feature_map(
    new_data,
    "patient_id",
    feature_map,
    target_features = c("gene_a", "gene_b")
  )
  changed <- new_data
  changed$outcome <- 1L - changed$outcome
  changed$unrelated <- -changed$unrelated
  repeated <- sb_apply_feature_map(
    changed,
    "patient_id",
    feature_map,
    target_features = c("gene_a", "gene_b")
  )

  expect_s3_class(observed, "stormbreaker_mapped_features")
  expect_identical(observed$status, "success")
  expect_identical(observed$id_col, "patient_id")
  expect_identical(observed$n, 3L)
  expect_identical(
    observed$mapping$target_feature,
    c("gene_a", "gene_b")
  )
  expect_identical(
    observed$mapping$source_feature,
    c("ENSG_A", "ENSG_B")
  )
  expect_identical(
    observed$mapping$annotation_release,
    rep("declared_release", 2L)
  )
  expect_identical(
    names(observed$data),
    c("patient_id", "gene_a", "gene_b")
  )
  expect_identical(observed$data$patient_id, new_data$patient_id)
  expect_identical(observed$data$gene_a, new_data$ENSG_A)
  expect_identical(observed$data$gene_b, new_data$ENSG_B)
  expect_identical(repeated, observed)

  signature <- sb_linear_signature(
    "mapped_signature",
    "binary",
    coefficients = c(gene_a = 0.5, gene_b = -0.25),
    intercept = -0.1
  )
  prediction <- sb_predict_signature(
    signature,
    observed$data,
    "patient_id"
  )
  expected_score <- -0.1 + 0.5 * new_data$ENSG_A -
    0.25 * new_data$ENSG_B
  expect_identical(prediction$status, "success")
  expect_equal(
    prediction$predictions$probability,
    stats::plogis(expected_score),
    tolerance = 0
  )
})

test_that("mapped features feed an unchanged final prediction", {
  development <- make_sb_binary_fixture(n = 60L, seed = 9951L)
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  direct_data <- make_sb_binary_fixture(n = 6L, seed = 9952L)
  direct_data$patient_id <- paste0("new_", direct_data$patient_id)
  source_data <- direct_data
  source_names <- paste0("source_", seq_along(fit$required_features))
  names(source_data)[match(fit$required_features, names(source_data))] <-
    source_names
  feature_map <- data.frame(
    source_feature = rev(source_names),
    target_feature = rev(fit$required_features),
    stringsAsFactors = FALSE
  )

  mapped <- sb_apply_feature_map(
    source_data,
    "patient_id",
    feature_map,
    fit$required_features
  )
  direct <- sb_predict_final(fit, direct_data, "patient_id")
  through_map <- sb_predict_final(fit, mapped$data, "patient_id")

  expect_identical(through_map, direct)
})

test_that("feature mapping refuses ambiguity and missing sources", {
  data <- data.frame(
    patient_id = c("p1", "p2"),
    source_a = c(1, 2),
    source_b = c(3, 4),
    check.names = FALSE
  )
  one_to_many <- data.frame(
    source_feature = c("source_a", "source_a"),
    target_feature = c("gene_a", "gene_b")
  )
  many_to_one <- data.frame(
    source_feature = c("source_a", "source_b"),
    target_feature = c("gene_a", "gene_a")
  )
  incomplete <- data.frame(
    source_feature = "source_a",
    target_feature = "gene_a"
  )
  missing_source <- data.frame(
    source_feature = c("source_a", "source_c"),
    target_feature = c("gene_a", "gene_b")
  )

  expect_s3_class(
    capture_sb_error(sb_apply_feature_map(
      data, "patient_id", one_to_many, c("gene_a", "gene_b")
    )),
    "stormbreaker_non_bijective_feature_map"
  )
  expect_s3_class(
    capture_sb_error(sb_apply_feature_map(
      data, "patient_id", many_to_one, "gene_a"
    )),
    "stormbreaker_non_bijective_feature_map"
  )
  expect_s3_class(
    capture_sb_error(sb_apply_feature_map(
      data, "patient_id", incomplete, c("gene_a", "gene_b")
    )),
    "stormbreaker_feature_map_target_mismatch"
  )
  expect_s3_class(
    capture_sb_error(sb_apply_feature_map(
      data, "patient_id", missing_source, c("gene_a", "gene_b")
    )),
    "stormbreaker_missing_feature_map_source"
  )
})

test_that("feature mapping never aggregates or repairs predictor values", {
  nonnumeric <- data.frame(
    patient_id = c("p1", "p2"),
    source_a = c("low", "high")
  )
  missing <- data.frame(
    patient_id = c("p1", "p2"),
    source_a = c(1, NA_real_)
  )
  mapping <- data.frame(
    source_feature = "source_a",
    target_feature = "gene_a"
  )

  expect_s3_class(
    capture_sb_error(sb_apply_feature_map(
      nonnumeric, "patient_id", mapping, "gene_a"
    )),
    "stormbreaker_bad_adapter_predictor_type"
  )
  expect_s3_class(
    capture_sb_error(sb_apply_feature_map(
      missing, "patient_id", mapping, "gene_a"
    )),
    "stormbreaker_bad_adapter_predictor_values"
  )
})

# ---- linear-signature ----

test_that("a locked binary signature computes its published logistic formula", {
  signature <- sb_linear_signature(
    signature_id = "published_binary_demo",
    task_type = "binary",
    coefficients = c(gene_b = -0.4, gene_a = 0.7),
    intercept = -0.2,
    source = "Example publication formula",
    input_description = "Preprocessed values on the publication scale"
  )
  new_data <- data.frame(
    patient_id = paste0("p", seq_len(8L)),
    gene_a = seq(-1, 1, length.out = 8L),
    unrelated_outcome = rep(c(0L, 1L), 4L),
    gene_b = seq(0.5, -0.5, length.out = 8L),
    stringsAsFactors = FALSE
  )

  observed <- sb_predict_signature(signature, new_data, "patient_id")
  changed_outcome <- new_data
  changed_outcome$unrelated_outcome <- 1L - changed_outcome$unrelated_outcome
  repeated <- sb_predict_signature(
    signature,
    changed_outcome[c("gene_b", "patient_id", "unrelated_outcome", "gene_a")],
    "patient_id"
  )
  expected_score <- -0.2 - 0.4 * new_data$gene_b + 0.7 * new_data$gene_a

  expect_s3_class(signature, "stormbreaker_linear_signature")
  expect_identical(signature$features, c("gene_b", "gene_a"))
  expect_identical(signature$source, "Example publication formula")
  expect_s3_class(observed, "stormbreaker_final_prediction")
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$pipeline_id, "signature__published_binary_demo")
  expect_identical(observed$learner, "locked_linear_signature")
  expect_equal(
    observed$predictions$probability,
    stats::plogis(expected_score),
    tolerance = 0
  )
  expect_identical(
    repeated$predictions$probability,
    observed$predictions$probability
  )
  expect_identical(
    observed$signature$coefficients,
    c(gene_b = -0.4, gene_a = 0.7)
  )

  endpoint <- sb_binary_endpoint(
    data.frame(
      patient_id = new_data$patient_id,
      outcome = rep(c(0L, 1L), 4L)
    ),
    "patient_id",
    "outcome"
  )
  evaluation <- sb_evaluate_external(observed, endpoint)
  expect_identical(evaluation$status, "success")
  expect_identical(
    evaluation$scored_predictions$probability,
    observed$predictions$probability
  )
})

test_that("binary probabilities require a stated intercept", {
  omitted <- capture_sb_error(sb_linear_signature(
    "omitted_binary_intercept",
    "binary",
    c(gene_a = 1)
  ))
  explicit_zero <- sb_linear_signature(
    "explicit_zero_binary_intercept",
    "binary",
    c(gene_a = 1),
    intercept = 0
  )
  prediction <- sb_predict_signature(
    explicit_zero,
    data.frame(patient_id = "p1", gene_a = 0),
    "patient_id"
  )
  survival_default <- sb_linear_signature(
    "survival_zero_reference",
    "survival",
    c(gene_a = 1)
  )

  expect_s3_class(omitted, "stormbreaker_missing_signature_intercept")
  expect_identical(prediction$predictions$probability, 0.5)
  expect_identical(survival_default$intercept, 0)
})

test_that("a survival signature can remain risk-only or use a stated baseline", {
  coefficients <- c(gene_a = 0.5, gene_b = -0.25)
  risk_signature <- sb_linear_signature(
    "published_survival_risk",
    "survival",
    coefficients,
    intercept = 0.1
  )
  horizon_signature <- sb_linear_signature(
    "published_survival_horizon",
    "survival",
    coefficients,
    intercept = 0.1,
    prediction_horizon = 6,
    time_unit = "month",
    baseline_survival = 0.82
  )
  new_data <- data.frame(
    patient_id = paste0("s", seq_len(16L)),
    gene_b = seq(0.8, -0.8, length.out = 16L),
    event = rep(c(1L, 0L), 8L),
    gene_a = seq(-1.2, 1.2, length.out = 16L),
    stringsAsFactors = FALSE
  )
  expected_risk <- 0.1 + 0.5 * new_data$gene_a - 0.25 * new_data$gene_b

  risk_prediction <- sb_predict_signature(
    risk_signature,
    new_data,
    "patient_id"
  )
  horizon_prediction <- sb_predict_signature(
    horizon_signature,
    new_data,
    "patient_id"
  )

  expect_identical(risk_prediction$status, "success")
  expect_equal(risk_prediction$predictions$risk, expected_risk, tolerance = 0)
  expect_identical(
    "survival_probability" %in% names(risk_prediction$predictions),
    FALSE
  )
  expect_equal(
    horizon_prediction$predictions$risk,
    expected_risk,
    tolerance = 0
  )
  expect_equal(
    horizon_prediction$predictions$survival_probability,
    0.82^exp(expected_risk),
    tolerance = 0
  )
  expect_identical(
    unique(horizon_prediction$predictions$prediction_horizon),
    6
  )
  expect_identical(
    unique(horizon_prediction$predictions$time_unit),
    "month"
  )

  endpoint_data <- data.frame(
    patient_id = new_data$patient_id,
    time = as.numeric(seq_len(16L)),
    event = rep(c(1L, 0L), 8L)
  )
  evaluation <- sb_evaluate_external(
    horizon_prediction,
    sb_surv_endpoint(
      endpoint_data,
      "patient_id",
      "time",
      "event",
      "month",
      prediction_horizon = 6
    )
  )
  expect_identical(evaluation$status, "success")
  expect_identical(evaluation$metrics$ipcw_brier$status, "success")
  expect_identical(evaluation$metrics$time_auc$status, "success")
  expect_identical(evaluation$metrics$uno_c$status, "success")

  zero_baseline <- sb_linear_signature(
    "zero_baseline",
    "survival",
    c(gene_a = 1),
    prediction_horizon = 6,
    time_unit = "month",
    baseline_survival = 0
  )
  zero_prediction <- sb_predict_signature(
    zero_baseline,
    data.frame(patient_id = "extreme", gene_a = -1000),
    "patient_id"
  )
  expect_identical(
    zero_prediction$predictions$survival_probability,
    0
  )
})

test_that("signature prediction failures stay explicit without substitution", {
  signature <- sb_linear_signature(
    "missing_feature_demo",
    "binary",
    c(gene_a = 0.5, gene_b = -0.2),
    intercept = 0
  )
  missing <- data.frame(
    patient_id = c("p1", "p2"),
    gene_a = c(0.1, 0.2)
  )
  nonnumeric <- data.frame(
    patient_id = c("p1", "p2"),
    gene_a = c(0.1, 0.2),
    gene_b = c("low", "high")
  )

  missing_result <- sb_predict_signature(signature, missing, "patient_id")
  nonnumeric_result <- sb_predict_signature(
    signature,
    nonnumeric,
    "patient_id"
  )
  partial_horizon <- capture_sb_error(sb_linear_signature(
    "partial_horizon",
    "survival",
    c(gene_a = 0.5),
    prediction_horizon = 12,
    time_unit = "month"
  ))
  unnamed <- capture_sb_error(sb_linear_signature(
    "unnamed",
    "binary",
    coefficients = c(0.5, -0.2)
  ))

  expect_identical(missing_result$status, "failed")
  expect_identical(missing_result$error_stage, "features")
  expect_identical(
    missing_result$error_class,
    "stormbreaker_missing_signature_features"
  )
  expect_identical(
    missing_result$predictions$status,
    rep("failed", 2L)
  )
  expect_identical(missing_result$learner, "locked_linear_signature")
  expect_identical(nonnumeric_result$status, "failed")
  expect_identical(nonnumeric_result$error_stage, "features")
  expect_s3_class(partial_horizon, "stormbreaker_bad_signature_horizon")
  expect_s3_class(unnamed, "stormbreaker_bad_signature_coefficients")
})

test_that("locked signatures reuse external comparison and paired bootstrap", {
  ids <- sprintf("e%02d", seq_len(40L))
  data <- data.frame(
    patient_id = ids,
    outcome = rep(c(0L, 1L), 20L),
    gene_a = seq(-1, 1, length.out = 40L),
    gene_b = cos(seq_len(40L) / 5),
    stringsAsFactors = FALSE
  )
  signatures <- list(
    sb_linear_signature(
      "published_reference",
      "binary",
      c(gene_a = 0.6, gene_b = -0.2),
      intercept = -0.1
    ),
    sb_linear_signature(
      "published_candidate",
      "binary",
      c(gene_a = 0.4, gene_b = 0.3),
      intercept = 0.2
    )
  )
  predictions <- lapply(signatures, function(signature) {
    sb_predict_signature(
      signature,
      data[c("patient_id", "gene_a", "gene_b")],
      "patient_id"
    )
  })
  comparison <- sb_compare_external(
    predictions,
    sb_binary_endpoint(data, "patient_id", "outcome")
  )
  paired <- sb_bootstrap_external_comparison(
    comparison,
    reference_pipeline_id = predictions[[1L]]$pipeline_id,
    reps = 12L,
    confidence_level = 0.90,
    seed = 9901L
  )

  expect_identical(comparison$status, "success")
  expect_identical(
    comparison$comparison$learner,
    rep("locked_linear_signature", 2L)
  )
  expect_identical(paired$status, "success")
  expect_identical(
    unique(paired$intervals$candidate_pipeline_id),
    predictions[[2L]]$pipeline_id
  )
  expect_identical(
    paired$pairing$draw_sharing,
    "same_patient_indices_for_all_models"
  )
})

test_that("signature predictions print their formula-only origin", {
  signature <- sb_linear_signature(
    "printed_signature",
    "binary",
    c(gene_a = 0.5),
    intercept = 0
  )
  prediction <- sb_predict_signature(
    signature,
    data.frame(patient_id = "patient_1", gene_a = 0.2),
    "patient_id"
  )

  printed <- paste(capture.output(print(prediction)), collapse = "\n")

  expect_match(printed, "locked formula-only prediction", fixed = TRUE)
  expect_match(printed, "printed_signature", fixed = TRUE)
  expect_identical(
    grepl("final state-only prediction", printed, fixed = TRUE),
    FALSE
  )
})

# ---- signature-inventory ----

test_that("a signature inventory retains every literature record in order", {
  records <- list(
    paper_a = list(
      signature_id = "complete_binary",
      task_type = "binary",
      coefficients = c(gene_a = 0.7),
      intercept = -0.2,
      source = "doi:paper-a"
    ),
    paper_b = list(
      signature_id = "binary_without_intercept",
      task_type = "binary",
      coefficients = c(gene_a = 0.5),
      source = "doi:paper-b"
    ),
    paper_c = list(
      task_type = "survival",
      coefficients = c(gene_a = 0.4),
      source = "doi:paper-c"
    ),
    paper_d = sb_linear_signature(
      "complete_survival",
      "survival",
      c(gene_a = 0.3)
    ),
    paper_e = "unparseable record"
  )

  observed <- sb_signature_inventory(records)

  expect_identical(
    names(observed),
    c(
      "record_id", "status", "signature_id", "task_type", "n_features",
      "output", "source", "issue_code", "message"
    )
  )
  expect_identical(observed$record_id, names(records))
  expect_identical(
    observed$status,
    c("admitted", "rejected", "rejected", "admitted", "rejected")
  )
  expect_identical(observed$source[1:3], c(
    "doi:paper-a", "doi:paper-b", "doi:paper-c"
  ))
  expect_identical(
    observed$issue_code[2L],
    "stormbreaker_missing_signature_intercept"
  )
  expect_identical(
    observed$issue_code[3L],
    "stormbreaker_missing_signature_argument"
  )
  expect_true(is.na(observed$signature_id[3L]))
  expect_identical(observed$output[4L], "risk")
  expect_identical(
    observed$issue_code[5L],
    "stormbreaker_bad_signature_spec"
  )
  expect_identical(nrow(observed), length(records))
})

test_that("a signature inventory requires unique non-empty record IDs", {
  valid_spec <- list(
    signature_id = "complete",
    task_type = "survival",
    coefficients = c(gene_a = 0.4)
  )
  empty <- capture_sb_error(sb_signature_inventory(list()))
  unnamed <- capture_sb_error(sb_signature_inventory(list(valid_spec)))
  duplicated <- capture_sb_error(sb_signature_inventory(structure(
    list(valid_spec, valid_spec),
    names = c("paper", "paper")
  )))
  blank <- capture_sb_error(sb_signature_inventory(structure(
    list(valid_spec),
    names = " "
  )))
  bad_container <- capture_sb_error(sb_signature_inventory("not a list"))
  bad_frame <- capture_sb_error(sb_signature_inventory(data.frame(record = 1)))
  single_spec <- capture_sb_error(sb_signature_inventory(valid_spec))
  single_object <- capture_sb_error(sb_signature_inventory(
    sb_linear_signature(
      "single_object",
      "survival",
      c(gene_a = 0.4)
    )
  ))
  legal_core_named_record <- sb_signature_inventory(list(
    task_type = valid_spec
  ))

  expect_s3_class(empty, "stormbreaker_bad_signature_inventory")
  expect_s3_class(unnamed, "stormbreaker_bad_signature_inventory")
  expect_s3_class(duplicated, "stormbreaker_bad_signature_inventory")
  expect_s3_class(blank, "stormbreaker_bad_signature_inventory")
  expect_s3_class(bad_container, "stormbreaker_bad_signature_inventory")
  expect_s3_class(bad_frame, "stormbreaker_bad_signature_inventory")
  expect_s3_class(single_spec, "stormbreaker_bad_signature_inventory")
  expect_s3_class(single_object, "stormbreaker_bad_signature_inventory")
  expect_identical(legal_core_named_record$record_id, "task_type")
  expect_identical(legal_core_named_record$status, "admitted")
})

test_that("a signature inventory never disguises an internal linter error", {
  local_mocked_bindings(
    sb_lint_signature = function(spec) {
      sb_abort(
        "injected inventory failure",
        "stormbreaker_internal_inventory_failure"
      )
    }
  )

  error <- tryCatch(
    sb_signature_inventory(list(record = list(
      signature_id = "complete",
      task_type = "survival",
      coefficients = c(gene_a = 0.4)
    ))),
    error = function(error) error
  )

  expect_s3_class(error, "stormbreaker_internal_inventory_failure")
  expect_identical(conditionMessage(error), "injected inventory failure")
})

# ---- signature-linter ----

test_that("the signature linter admits complete raw and constructed formulas", {
  raw_report <- sb_lint_signature(list(
    signature_id = "published_binary",
    task_type = "binary",
    coefficients = c(gene_a = 0.7, gene_b = -0.4),
    intercept = -0.2,
    source = "DOI supplied by the reviewer"
  ))
  survival <- sb_linear_signature(
    "published_survival",
    "survival",
    c(gene_a = 0.5),
    prediction_horizon = 12,
    time_unit = "month",
    baseline_survival = 0.8
  )
  object_report <- sb_lint_signature(survival)
  risk_report <- sb_lint_signature(list(
    signature_id = "published_risk_score",
    task_type = "survival",
    coefficients = c(gene_a = 0.5)
  ))

  expect_identical(
    names(raw_report),
    c(
      "status", "signature_id", "task_type", "n_features", "output",
      "source", "issue_code", "message"
    )
  )
  expect_identical(raw_report$status, "admitted")
  expect_identical(raw_report$signature_id, "published_binary")
  expect_identical(raw_report$task_type, "binary")
  expect_identical(raw_report$n_features, 2L)
  expect_identical(raw_report$output, "probability")
  expect_identical(raw_report$source, "DOI supplied by the reviewer")
  expect_true(is.na(raw_report$issue_code))
  expect_true(is.na(raw_report$message))

  expect_identical(object_report$status, "admitted")
  expect_identical(object_report$signature_id, "published_survival")
  expect_identical(object_report$n_features, 1L)
  expect_identical(object_report$output, "risk+survival_probability")
  expect_identical(risk_report$status, "admitted")
  expect_identical(risk_report$output, "risk")
})

test_that("the signature linter retains incomplete literature formulas", {
  missing_coefficients <- sb_lint_signature(list(
    signature_id = "missing_coefficients",
    task_type = "survival"
  ))
  missing_binary_intercept <- sb_lint_signature(list(
    signature_id = "missing_binary_intercept",
    task_type = "binary",
    coefficients = c(gene_a = 0.7)
  ))
  partial_horizon <- sb_lint_signature(list(
    signature_id = "partial_horizon",
    task_type = "survival",
    coefficients = c(gene_a = 0.5),
    prediction_horizon = 12,
    time_unit = "month"
  ))
  unknown_argument <- sb_lint_signature(list(
    signature_id = "unknown_argument",
    task_type = "survival",
    coefficients = c(gene_a = 0.5),
    published_cutoff = 0.3
  ))
  bad_time_unit <- sb_lint_signature(list(
    signature_id = "bad_time_unit",
    task_type = "survival",
    coefficients = c(gene_a = 0.5),
    prediction_horizon = 12,
    time_unit = "week",
    baseline_survival = 0.8
  ))
  unnamed <- sb_lint_signature(list(
    "unnamed_specification",
    "survival",
    c(gene_a = 0.5)
  ))
  bad_container <- sb_lint_signature("not a signature specification")

  expect_identical(missing_coefficients$status, "rejected")
  expect_identical(
    missing_coefficients$issue_code,
    "stormbreaker_missing_signature_argument"
  )
  expect_identical(missing_coefficients$signature_id, "missing_coefficients")
  expect_true(is.na(missing_coefficients$n_features))

  expect_identical(missing_binary_intercept$status, "rejected")
  expect_identical(
    missing_binary_intercept$issue_code,
    "stormbreaker_missing_signature_intercept"
  )
  expect_identical(missing_binary_intercept$n_features, 1L)
  expect_identical(missing_binary_intercept$task_type, "binary")

  expect_identical(partial_horizon$status, "rejected")
  expect_identical(
    partial_horizon$issue_code,
    "stormbreaker_bad_signature_horizon"
  )
  expect_identical(unknown_argument$status, "rejected")
  expect_identical(
    unknown_argument$issue_code,
    "stormbreaker_unknown_signature_argument"
  )
  expect_identical(bad_time_unit$status, "rejected")
  expect_identical(bad_time_unit$issue_code, "stormbreaker_bad_time_unit")
  expect_identical(unnamed$status, "rejected")
  expect_identical(unnamed$issue_code, "stormbreaker_bad_signature_spec")
  expect_identical(bad_container$status, "rejected")
  expect_identical(
    bad_container$issue_code,
    "stormbreaker_bad_signature_spec"
  )
})

test_that("the signature linter does not disguise unexpected internal errors", {
  local_mocked_bindings(
    sb_signature_from_lint_spec = function(spec) {
      sb_abort(
        "injected internal signature failure",
        "stormbreaker_internal_signature_failure"
      )
    }
  )

  error <- tryCatch(
    sb_lint_signature(list(
      signature_id = "otherwise_complete",
      task_type = "binary",
      coefficients = c(gene_a = 0.7),
      intercept = -0.2
    )),
    error = function(error) error
  )

  expect_s3_class(error, "stormbreaker_internal_signature_failure")
  expect_identical(
    conditionMessage(error),
    "injected internal signature failure"
  )
})
})

# ---- learner-state ----
local({
sb_state_contains_live_object <- function(value) {
  forbidden_class <- inherits(
    value,
    c(
      "coxph", "glmnetfit", "lognet", "coxnet",
      "lgb.Booster", "catboost.Model"
    )
  )
  if (is.function(value) || is.environment(value) ||
      typeof(value) == "externalptr" || forbidden_class) {
    return(TRUE)
  }
  if (!is.list(value)) {
    return(FALSE)
  }
  any(vapply(unclass(value), sb_state_contains_live_object, logical(1)))
}

test_that("ordinary Cox state reconstructs risk and horizon probability", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  x <- d[c("gene_b", "gene_a")]
  adapter <- sb_coxph_learner_adapter()
  model <- adapter$fit(x, d$time, d$event)
  direct_risk <- adapter$predict(model, x)
  direct_probability <- adapter$predict_survival(
    model,
    x,
    d$time,
    d$event,
    x,
    prediction_horizon = 12
  )$probability
  state <- adapter$extract_state(model, list(
    x = x,
    time = d$time,
    event = as.integer(d$event),
    prediction_horizon = 12,
    time_unit = "month"
  ))
  rm(model)
  observed <- adapter$predict_state(state, x)
  reconstructed_risk <- unname(as.vector(
    as.matrix(x[state$features]) %*% state$coefficient
  ))
  reconstructed_probability <- exp(
    -state$baseline_hazard$cumulative_baseline_hazard *
      exp(reconstructed_risk)
  )

  expect_s3_class(state, "stormbreaker_learner_state")
  expect_identical(state$learner, "coxph")
  expect_identical(inherits(state, "coxph"), FALSE)
  expect_false(sb_state_contains_live_object(state))
  expect_equal(observed$risk, direct_risk, tolerance = 1e-12)
  expect_equal(observed$risk, reconstructed_risk, tolerance = 1e-12)
  expect_equal(
    observed$survival_probability,
    direct_probability,
    tolerance = 1e-12
  )
  expect_equal(
    observed$survival_probability,
    reconstructed_probability,
    tolerance = 1e-12
  )
  expect_identical(
    adapter$predict_state(unserialize(serialize(state, NULL)), x),
    observed
  )
})

test_that("glmnet Cox state predicts without retaining a glmnet fit", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_glmnet_cox_learner_adapter(
    alpha = 0.5,
    lambda = 0.05,
    standardize = TRUE
  )
  model <- adapter$fit(x, d$time, d$event)
  direct_risk <- adapter$predict(model, x)
  direct_probability <- adapter$predict_survival(
    model,
    x,
    d$time,
    d$event,
    x,
    prediction_horizon = 12
  )$probability
  state <- adapter$extract_state(model, list(
    x = x,
    time = d$time,
    event = as.integer(d$event),
    prediction_horizon = 12,
    time_unit = "month"
  ))
  rm(model)
  observed <- adapter$predict_state(state, x)
  reconstructed_risk <- unname(as.vector(
    as.matrix(x[state$features]) %*% state$coefficient
  ))
  reconstructed_probability <- exp(
    -state$baseline_hazard$cumulative_baseline_hazard *
      exp(reconstructed_risk)
  )

  expect_identical(state$learner, "glmnet_cox")
  expect_identical(inherits(state, "glmnetfit"), FALSE)
  expect_false(sb_state_contains_live_object(state))
  expect_equal(observed$risk, direct_risk, tolerance = 1e-12)
  expect_equal(observed$risk, reconstructed_risk, tolerance = 1e-12)
  expect_equal(
    observed$survival_probability,
    direct_probability,
    tolerance = 1e-12
  )
  expect_equal(
    observed$survival_probability,
    reconstructed_probability,
    tolerance = 1e-12
  )
  expect_identical(
    adapter$predict_state(unserialize(serialize(state, NULL)), x),
    observed
  )
})

test_that("logistic states reconstruct class-one probabilities", {
  d <- sb_example_binary(n = 80L, seed = 20260819L)
  x <- d[c("gene_b", "gene_a")]
  adapter <- sb_logistic_learner_adapter()
  model <- adapter$fit(x, d$outcome)
  direct <- adapter$predict(model, x)
  state <- adapter$extract_state(model, list())
  rm(model)
  observed <- adapter$predict_state(state, x)
  linear_predictor <- state$intercept + unname(as.vector(
    as.matrix(x[state$features]) %*% state$coefficient
  ))

  expect_identical(state$learner, "logistic")
  expect_identical(inherits(state, "glm"), FALSE)
  expect_false(sb_state_contains_live_object(state))
  expect_equal(observed$probability, direct, tolerance = 1e-12)
  expect_equal(
    observed$probability,
    stats::plogis(linear_predictor),
    tolerance = 1e-12
  )
})

test_that("glmnet logistic state predicts from saved coefficients only", {
  skip_if_not_installed("glmnet")

  d <- sb_example_binary(n = 80L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_glmnet_logistic_learner_adapter(
    alpha = 0.5,
    lambda = 0.05,
    standardize = TRUE
  )
  model <- adapter$fit(x, d$outcome)
  direct <- adapter$predict(model, x)
  state <- adapter$extract_state(model, list())
  rm(model)
  observed <- adapter$predict_state(state, x)

  expect_identical(state$learner, "glmnet_logistic")
  expect_identical(inherits(state, "glmnetfit"), FALSE)
  expect_false(sb_state_contains_live_object(state))
  expect_equal(observed$probability, direct, tolerance = 1e-12)
  expect_identical(
    adapter$predict_state(unserialize(serialize(state, NULL)), x),
    observed
  )
})

test_that("Gaussian Naive Bayes state retains only portable distributions", {
  d <- sb_example_binary(n = 80L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_gaussian_naive_bayes_learner_adapter(1e-9)
  model <- adapter$fit(x, d$outcome)
  direct <- adapter$predict(model, x)
  state <- adapter$extract_state(model, list())
  rm(model)
  observed <- adapter$predict_state(state, x)

  expect_identical(state$learner, "gaussian_naive_bayes")
  expect_false(sb_state_contains_live_object(state))
  expect_equal(observed$probability, direct, tolerance = 1e-14)
  expect_identical(
    adapter$predict_state(unserialize(serialize(state, NULL)), x),
    observed
  )
})

test_that("LightGBM state reloads its saved model string", {
  skip_if_not_installed("lightgbm")

  d <- sb_example_binary(n = 80L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_lightgbm_binary_learner_adapter(
    nrounds = 5L,
    learning_rate = 0.1,
    num_leaves = 7L,
    min_data_in_leaf = 3L,
    seed = 101L
  )
  model <- adapter$fit(x, d$outcome)
  direct <- adapter$predict(model, x)
  state <- adapter$extract_state(model, list())
  rm(model)
  observed <- adapter$predict_state(state, x)

  expect_identical(state$learner, "lightgbm_binary")
  expect_type(state$model_string, "character")
  expect_identical("booster" %in% names(state), FALSE)
  expect_false(sb_state_contains_live_object(state))
  expect_equal(observed$probability, direct, tolerance = 1e-12)
  expect_equal(
    adapter$predict_state(unserialize(serialize(state, NULL)), x)$probability,
    observed$probability,
    tolerance = 1e-12
  )
})

test_that("CatBoost state reloads its saved raw model", {
  skip_if_not_installed("catboost")

  d <- sb_example_binary(n = 80L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_catboost_binary_learner_adapter(
    iterations = 5L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 3,
    seed = 101L
  )
  model <- adapter$fit(x, d$outcome)
  direct <- adapter$predict(model, x)
  state <- adapter$extract_state(model, list())
  rm(model)
  observed <- adapter$predict_state(state, x)

  expect_identical(state$learner, "catboost_binary")
  expect_type(state$model_raw, "raw")
  expect_identical("booster" %in% names(state), FALSE)
  expect_false(sb_state_contains_live_object(state))
  expect_equal(observed$probability, direct, tolerance = 1e-12)
  expect_equal(
    adapter$predict_state(unserialize(serialize(state, NULL)), x)$probability,
    observed$probability,
    tolerance = 1e-12
  )
})

test_that("invalid portable states fail visibly without fallback", {
  d <- make_sb_fixture(n = 40L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b")]
  adapter <- sb_coxph_learner_adapter()
  model <- adapter$fit(x, d$time, d$event)
  state <- adapter$extract_state(model, list())
  state$coefficient[1L] <- NA_real_
  coefficient_error <- capture_sb_error(adapter$predict_state(state, x))
  expect_s3_class(coefficient_error, "stormbreaker_bad_learner_state")

  if (sb_lightgbm_available()) {
    binary <- sb_example_binary(n = 60L, seed = 20260819L)
    binary_x <- binary[c("gene_a", "gene_b")]
    lightgbm <- sb_lightgbm_binary_learner_adapter(
      nrounds = 3L,
      num_leaves = 5L,
      min_data_in_leaf = 2L
    )
    lightgbm_model <- lightgbm$fit(binary_x, binary$outcome)
    lightgbm_state <- lightgbm$extract_state(lightgbm_model, list())
    lightgbm_state$model_string <- "not a LightGBM model"
    lightgbm_error <- capture_sb_error(
      lightgbm$predict_state(lightgbm_state, binary_x)
    )
    expect_s3_class(
      lightgbm_error,
      "stormbreaker_invalid_lightgbm_state"
    )
  }

  if (sb_catboost_available()) {
    binary <- sb_example_binary(n = 60L, seed = 20260819L)
    binary_x <- binary[c("gene_a", "gene_b")]
    catboost <- sb_catboost_binary_learner_adapter(
      iterations = 3L,
      depth = 2L
    )
    catboost_model <- catboost$fit(binary_x, binary$outcome)
    catboost_state <- catboost$extract_state(catboost_model, list())
    catboost_state$model_raw <- charToRaw("not a serialized model")
    catboost_error <- capture_sb_error(
      catboost$predict_state(catboost_state, binary_x)
    )
    expect_s3_class(catboost_error, "stormbreaker_invalid_catboost_state")
  }
})

test_that("portable tree states still require their declared dependencies", {
  skip_if_not_installed("lightgbm")

  d <- sb_example_binary(n = 60L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b")]
  adapter <- sb_lightgbm_binary_learner_adapter(
    nrounds = 3L,
    num_leaves = 5L,
    min_data_in_leaf = 2L
  )
  model <- adapter$fit(x, d$outcome)
  state <- adapter$extract_state(model, list())
  testthat::local_mocked_bindings(
    sb_lightgbm_available = function() FALSE
  )

  error <- capture_sb_error(adapter$predict_state(state, x))
  expect_s3_class(error, "stormbreaker_missing_lightgbm")
})

test_that("portable CatBoost state still requires CatBoost", {
  skip_if_not_installed("catboost")

  d <- sb_example_binary(n = 60L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b")]
  adapter <- sb_catboost_binary_learner_adapter(
    iterations = 3L,
    depth = 2L
  )
  model <- adapter$fit(x, d$outcome)
  state <- adapter$extract_state(model, list())
  testthat::local_mocked_bindings(
    sb_catboost_available = function() FALSE
  )

  error <- capture_sb_error(adapter$predict_state(state, x))
  expect_s3_class(error, "stormbreaker_missing_catboost")
})
})

# ---- legacy-name-compatibility ----
local({
test_that("the full learner registry uses clean canonical identifiers", {
  registry <- sb_pipeline_registry()
  expect_length(unique(registry$learner), 176L)
  expect_length(grep("_fixed$", registry$learner), 0L)
  for (old in names(sb_legacy_learner_names)) {
    clean <- sb_legacy_learner_names[[old]]
    expect_identical(sb_pipeline(learner = old), sb_pipeline(learner = clean))
  }
})

test_that("name compatibility preserves feature names and parameter values", {
  pipeline <- sb_pipeline(learner = "ranger_binary", ranger_num_trees = 23L)
  legacy <- pipeline
  legacy$learner <- "ranger_binary_fixed"
  legacy$id <- sub("ranger_binary", "ranger_binary_fixed", legacy$id, fixed = TRUE)
  expect_identical(sb_upgrade_legacy_names(legacy), pipeline)
  expect_identical(sb_ranger_binary_fixed_params(), sb_ranger_binary_params())
  genes <- data.frame(ranger_binary_fixed = c(1, 2), check.names = FALSE)
  expect_identical(sb_upgrade_legacy_names(genes), genes)
  expect_identical(sb_canonical_learner_id("invented_fixed"), "invented_fixed")
  expect_identical(sb_canonical_learner_id("fixed_features"), "fixed_features")
})
})

# ---- runtime-dependencies ----
local({
test_that("dependency descriptions separate R packages and Python modules", {
  expect_identical(
    sb_dependency_components("mboost + partykit"),
    list(
      r_packages = c("mboost", "partykit"),
      python_modules = character()
    )
  )
  expect_identical(
    sb_dependency_components(
      "reticulate + Python numpy + pandas + torch + torchtuples + pycox"
    ),
    list(
      r_packages = "reticulate",
      python_modules = c(
        "numpy", "pandas", "torch", "torchtuples", "pycox"
      )
    )
  )
})

test_that("dependency status distinguishes missing R and unchecked Python", {
  local_mocked_bindings(
    sb_r_package_installed = function(package) {
      package %in% c("mboost", "reticulate")
    }
  )

  missing <- sb_dependency_status_row("mboost + partykit")
  python <- sb_dependency_status_row(
    "reticulate + Python numpy + pandas + torch + torchtuples + pycox"
  )

  expect_identical(missing$missing_r_packages, "partykit")
  expect_identical(missing$dependency_status, "missing_r_packages")
  expect_identical(python$missing_r_packages, NA_character_)
  expect_identical(python$dependency_status, "python_not_checked")
})

test_that("installed checks require one logical value", {
  error <- tryCatch(
    sb_pipeline_registry(check_installed = 1),
    error = function(error) error
  )

  expect_s3_class(error, "stormbreaker_bad_registry_check_installed")
})
})

# ---- task-data ----
local({
# ---- task-binary ----

test_that("group roles are explicit without weakening unique observation IDs", {
  d <- sb_example_binary(20L, 3755L)
  d$group <- rep(1:10, each = 2)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a", group_col = "group")
  expect_identical(task$data, d)
  expect_identical(task$n, 20L)
  expect_identical(task$group_col, "group")
  for (column in c("patient_id", "outcome", "gene_a")) {
    error <- capture_sb_error(sb_binary_task(
      d, "patient_id", "outcome", "gene_a", group_col = column
    ))
    expect_s3_class(error, "stormbreaker_overlapping_group_roles")
  }
  missing <- capture_sb_error(sb_binary_task(
    d, "patient_id", "outcome", "gene_a", group_col = "absent"
  ))
  expect_s3_class(missing, "stormbreaker_missing_group_column")
  for (bad in list(rep(NA_character_, 20), rep(" ", 20),
                  rep(Inf, 20), as.list(1:20))) {
    changed <- d
    changed$group <- bad
    error <- capture_sb_error(sb_binary_task(
      changed, "patient_id", "outcome", "gene_a", group_col = "group"
    ))
    expect_s3_class(error, "stormbreaker_bad_group_ids")
  }
  d$patient_id[2] <- d$patient_id[1]
  error <- capture_sb_error(sb_binary_task(
    d, "patient_id", "outcome", "gene_a", group_col = "group"
  ))
  expect_s3_class(error, "stormbreaker_duplicate_id")
})

test_that("a valid binary task records explicit class semantics", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )

  expect_s3_class(task, "stormbreaker_binary_task")
  expect_identical(task$n, nrow(d))
  expect_identical(task$n_positive, sum(d$outcome == 1L))
  expect_identical(task$n_negative, sum(d$outcome == 0L))
  expect_identical(task$positive_class, 1L)
  expect_identical(task$negative_class, 0L)
  expect_identical(task$predictors, c("gene_a", "gene_b"))
  expect_identical(task$independent_unit, "patient")
})

test_that("logical binary outcomes preserve explicit zero-one meaning", {
  d <- make_sb_binary_fixture()
  d$outcome <- as.logical(d$outcome)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")

  expect_identical(task$n_positive, sum(d$outcome))
  expect_identical(task$positive_class, 1L)
})

test_that("binary outcomes require finite zero and one with both classes", {
  d <- make_sb_binary_fixture()
  d$outcome[1] <- 2L
  value_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "gene_a")
  )
  expect_s3_class(value_error, "stormbreaker_bad_binary_outcome")

  d <- make_sb_binary_fixture()
  d$outcome <- 1L
  negative_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "gene_a")
  )
  expect_s3_class(negative_error, "stormbreaker_no_negative_class")

  d$outcome <- 0L
  positive_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "gene_a")
  )
  expect_s3_class(positive_error, "stormbreaker_no_positive_class")
})

test_that("binary task roles and numeric predictor values are explicit", {
  d <- make_sb_binary_fixture()
  role_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "outcome")
  )
  expect_s3_class(role_error, "stormbreaker_reserved_predictor")

  d$group <- factor(rep(c("A", "B"), length.out = nrow(d)))
  type_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "group")
  )
  expect_s3_class(type_error, "stormbreaker_bad_predictor_type")

  d$gene_a[1] <- NA_real_
  value_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "gene_a")
  )
  expect_s3_class(value_error, "stormbreaker_bad_predictor_values")
})

test_that("binary tasks reserve the model intercept name", {
  d <- make_sb_binary_fixture()
  d[["(Intercept)"]] <- d$gene_a

  error <- capture_sb_error(sb_binary_task(
    d, "patient_id", "outcome", c("(Intercept)", "gene_b")
  ))

  expect_s3_class(error, "stormbreaker_reserved_intercept_name")
  expect_match(conditionMessage(error), "reserved for the model intercept")
})

test_that("binary patient identifiers remain unique and complete", {
  d <- make_sb_binary_fixture()
  d$patient_id[2] <- d$patient_id[1]
  duplicate_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "gene_a")
  )
  expect_s3_class(duplicate_error, "stormbreaker_duplicate_id")

  d <- make_sb_binary_fixture()
  d$patient_id[1] <- ""
  blank_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "gene_a")
  )
  expect_s3_class(blank_error, "stormbreaker_bad_id")
})

test_that("binary tasks preserve study-declared decision thresholds", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b"),
    decision_thresholds = c(low = 0.10, high = 0.35)
  )

  expect_identical(task$decision_thresholds, c(0.10, 0.35))
  expect_identical(
    stormbreaker:::sb_validate_binary_task_contract(task)$decision_thresholds,
    c(0.10, 0.35)
  )
})

test_that("binary tasks omit decision thresholds by default", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )

  expect_false("decision_thresholds" %in% names(task))
})

test_that("binary decision thresholds must be unique increasing probabilities", {
  d <- make_sb_binary_fixture()
  bad <- list(
    numeric(),
    c(0, 0.2),
    c(0.2, 1),
    c(NA_real_, 0.2),
    c(Inf, 0.2),
    c(0.2, 0.2),
    c(0.4, 0.2),
    "0.2"
  )

  for (value in bad) {
    error <- capture_sb_error(sb_binary_task(
      d,
      "patient_id",
      "outcome",
      c("gene_a", "gene_b"),
      decision_thresholds = value
    ))
    expect_s3_class(error, "stormbreaker_bad_decision_thresholds")
  }
})

test_that("binary task reconstruction rejects invalid threshold mutation", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b"),
    decision_thresholds = c(0.10, 0.30)
  )
  task$decision_thresholds <- c(0.10, 1)

  error <- capture_sb_error(
    stormbreaker:::sb_validate_binary_task_contract(task)
  )
  expect_s3_class(error, "stormbreaker_bad_decision_thresholds")
})

test_that("default binary tasks still enter CV and final fitting", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )
  pipeline <- sb_pipeline(learner = "logistic")

  cv <- sb_run_cv(task, pipeline, folds = 3, seed = 2026)
  fit <- sb_fit_final(task, pipeline)

  expect_identical(cv$status, "success")
  expect_identical(fit$status, "success")
  expect_false("decision_thresholds" %in% names(cv$task))
  expect_false("decision_thresholds" %in% names(fit$prediction))
})

# ---- task-survival ----

test_that("a valid survival task records endpoint semantics", {
  d <- make_sb_fixture()
  task <- sb_surv_task(
    d,
    id_col = "patient_id",
    time_col = "time",
    event_col = "event",
    predictors = c("gene_a", "gene_b"),
    time_unit = "month"
  )

  expect_s3_class(task, "stormbreaker_surv_task")
  expect_identical(task$n, nrow(d))
  expect_identical(task$n_events, sum(d$event))
  expect_identical(task$predictors, c("gene_a", "gene_b"))
  expect_identical(task$time_unit, "month")
  expect_null(task$prediction_horizon)
  expect_named(
    task,
    c(
      "data", "id_col", "time_col", "event_col", "predictors",
      "time_unit", "n", "n_events", "n_censored", "independent_unit"
    ),
    ignore.order = FALSE
  )
})

test_that("a survival task records one predeclared prediction horizon", {
  d <- make_sb_fixture()
  task <- sb_surv_task(
    d,
    id_col = "patient_id",
    time_col = "time",
    event_col = "event",
    predictors = c("gene_a", "gene_b"),
    time_unit = "month",
    prediction_horizon = 24L
  )

  expect_identical(task$prediction_horizon, 24)
  expect_identical(
    names(task),
    c(
      "data", "id_col", "time_col", "event_col", "predictors",
      "time_unit", "prediction_horizon", "n", "n_events", "n_censored",
      "independent_unit"
    )
  )
  printed <- paste(capture.output(print(task)), collapse = "\n")
  expect_match(printed, "Prediction horizon: 24 \\[ month \\]")
})

test_that("a survival task preserves declared event-risk thresholds", {
  d <- make_sb_fixture()
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 24,
    decision_thresholds = c(low = 0.10, high = 0.35)
  )

  expect_identical(task$decision_thresholds, c(0.10, 0.35))
  expect_identical(
    sb_validate_task_contract(task)$decision_thresholds,
    c(0.10, 0.35)
  )
  printed <- paste(capture.output(print(task)), collapse = "\n")
  expect_match(printed, "Decision thresholds: 2")
})

test_that("survival decision thresholds require a declared horizon", {
  d <- make_sb_fixture()
  error <- capture_sb_error(sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    "gene_a",
    "month",
    decision_thresholds = c(0.10, 0.35)
  ))

  expect_s3_class(
    error,
    "stormbreaker_decision_thresholds_require_prediction_horizon"
  )
})

test_that("survival decision thresholds use the shared strict validator", {
  d <- make_sb_fixture()
  invalid <- list(numeric(), c(0, 0.2), c(0.2, 1), c(0.3, 0.2), c(0.2, 0.2))

  for (value in invalid) {
    error <- capture_sb_error(sb_surv_task(
      d,
      "patient_id",
      "time",
      "event",
      "gene_a",
      "month",
      prediction_horizon = 24,
      decision_thresholds = value
    ))
    expect_s3_class(error, "stormbreaker_bad_decision_thresholds")
  }
})

test_that("survival task reconstruction rejects threshold mutation", {
  d <- make_sb_fixture()
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    "gene_a",
    "month",
    prediction_horizon = 24,
    decision_thresholds = c(0.10, 0.35)
  )
  task$decision_thresholds <- c(0.10, 1)

  error <- capture_sb_error(sb_validate_task_contract(task))
  expect_s3_class(error, "stormbreaker_bad_decision_thresholds")

  no_horizon <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  no_horizon$decision_thresholds <- 0.20
  horizon_error <- capture_sb_error(sb_validate_task_contract(no_horizon))
  expect_s3_class(
    horizon_error,
    "stormbreaker_decision_thresholds_require_prediction_horizon"
  )
})

test_that("a similarly named field cannot activate a prediction horizon", {
  d <- make_sb_fixture()
  task <- sb_surv_task(
    d,
    id_col = "patient_id",
    time_col = "time",
    event_col = "event",
    predictors = c("gene_a", "gene_b"),
    time_unit = "month"
  )
  task$prediction_horizons <- 5

  validated <- sb_validate_task_contract(task)

  expect_identical("prediction_horizon" %in% names(validated), FALSE)
})

test_that("prediction horizon must be one finite positive number", {
  d <- make_sb_fixture()
  invalid <- list(
    0,
    -1,
    NA_real_,
    Inf,
    "24",
    c(12, 24)
  )

  for (value in invalid) {
    error <- capture_sb_error(
      sb_surv_task(
        d, "patient_id", "time", "event", "gene_a", "month",
        prediction_horizon = value
      )
    )
    expect_s3_class(error, "stormbreaker_bad_prediction_horizon")
  }
})

test_that("duplicate patients and invalid endpoints are rejected", {
  d <- make_sb_fixture()
  d$patient_id[2] <- d$patient_id[1]
  duplicate_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "month"
    )
  )
  expect_s3_class(duplicate_error, "stormbreaker_duplicate_id")

  d <- make_sb_fixture()
  d$event[1] <- 2
  event_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "month"
    )
  )
  expect_s3_class(event_error, "stormbreaker_bad_event")

  d <- make_sb_fixture()
  d$time[1] <- -1
  time_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "month"
    )
  )
  expect_s3_class(time_error, "stormbreaker_bad_time")
})

test_that("predictor roles, values and time units are explicit", {
  d <- make_sb_fixture()
  role_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "event", "month"
    )
  )
  expect_s3_class(role_error, "stormbreaker_reserved_predictor")

  d$gene_a[1] <- NA_real_
  missing_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "month"
    )
  )
  expect_s3_class(missing_error, "stormbreaker_bad_predictor_values")

  d <- make_sb_fixture()
  unit_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "months"
    )
  )
  expect_s3_class(unit_error, "stormbreaker_bad_time_unit")
})

test_that("all-event, all-censored, and nonnumeric predictors are rejected", {
  d <- make_sb_fixture()
  d$event <- 1L
  censoring_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "month"
    )
  )
  expect_s3_class(censoring_error, "stormbreaker_no_censoring")

  d <- make_sb_fixture()
  d$event <- 0L
  event_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "month"
    )
  )
  expect_s3_class(event_error, "stormbreaker_no_events")

  d <- make_sb_fixture()
  d$gene_text <- as.character(d$gene_a)
  predictor_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_text", "month"
    )
  )
  expect_s3_class(predictor_error, "stormbreaker_bad_predictor_type")
})
test_that("survival task validation and subsetting retain patient grouping", {
  d <- sb_example_survival(40L, 3756L)
  d$group <- rep(1:20, each = 2)
  task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month",
    prediction_horizon = 6, decision_thresholds = c(0.1, 0.5), group_col = "group")
  expect_identical(sb_validate_task_contract(task), task)
  subset <- sb_subset_task(task, 1:30)
  expect_identical(subset$group_col, "group")
  expect_identical(subset$data$group, d$group[1:30])
  expect_identical(subset$prediction_horizon, 6)
  for (column in c("patient_id", "time", "event", "gene_a")) {
    error <- capture_sb_error(sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "month", group_col = column
    ))
    expect_s3_class(error, "stormbreaker_overlapping_group_roles")
  }
})
})
