# ---- binary-model-workflow ----
local({
test_that("one binary menu contains every existing prediction and feature method", {
  menu <- sb_models("binary")
  expect_equal(sum(menu$can_predict), nrow(sb_model_catalog("binary")))
  expect_equal(sum(menu$can_extract_features), length(sb_doctor_selector_menu("binary")))
  expect_identical(any(menu$model_id %in% sb_retired_learner_ids()), FALSE)
  expect_setequal(menu$model_id[menu$can_predict], sb_model_catalog("binary")$model_id)
  expect_identical(anyDuplicated(tolower(menu$name)), 0L)
  expect_identical(anyDuplicated(menu$model_id), 0L)
  expect_identical(menu$can_predict[menu$name == "Boruta"], FALSE)
  expect_identical(menu$can_extract_features[menu$name == "RF"], TRUE)
  expect_identical(menu$can_extract_features[menu$name == "Ridge"], FALSE)
  expect_identical(menu$r_packages[menu$name == "Boruta"], "Boruta, ranger")
})

test_that("ordered model names preserve native stage definitions and their order", {
  forward <- bi_pipeline(c("Lasso", "RF", "GBM"))
  old <- bi_pipeline("GBM", bi_selector(c("Lasso", "RF")))
  backward <- bi_pipeline(c("RF", "Lasso", "GBM"))
  expect_identical(forward, old)
  expect_identical(forward$id == backward$id, FALSE)
  expect_identical(bi_pipeline(c("PLSR", "Logistic"))$selector, "plsrglm_binary_weights")
  expect_identical(bi_pipeline(c("Enet", "Logistic"))$selector, "enet_binary_nonzero")
})

test_that("methods cannot silently acquire unsupported stage capabilities", {
  terminal <- tryCatch(bi_pipeline(c("Logistic", "Boruta")), error = identity)
  intermediate <- tryCatch(bi_pipeline(c("NaiveBayes", "Logistic")), error = identity)
  unknown <- tryCatch(bi_pipeline("not_a_method"), error = identity)
  expect_s3_class(terminal, "stormbreaker_method_cannot_predict")
  expect_s3_class(intermediate, "stormbreaker_method_cannot_extract")
  expect_s3_class(unknown, "stormbreaker_bad_model_choice")
})

test_that("three-stage enumeration covers the chosen method pool rather than presets", {
  plan <- sb_doctor_pipeline_set("binary", c("Lasso", "RF", "GBM"), c(1, 2, 3))
  expect_length(plan$pipelines, 39L)
  expect_length(plan$unavailable, 0L)
  expect_contains(names(plan$pipelines), "RF + Lasso + GBM")
  expect_contains(names(plan$pipelines), "Lasso + RF + GBM")
  expect_contains(names(plan$pipelines), "GBM + RF + Lasso")
  expect_identical(plan$pipelines[["RF + Lasso + GBM"]], bi_pipeline(c("RF", "Lasso", "GBM")))
})

test_that("feature-only methods are available in the same pool but not as predictors", {
  plan <- sb_doctor_pipeline_set("binary", c("Boruta", "Logistic"), c(1, 2))
  expect_named(plan$pipelines, c("Logistic", "Boruta + Logistic"))
  expect_length(plan$unavailable, 0L)
  expect_identical(plan$pipelines[[2L]]$selector, "boruta_binary_roughfix")
  expect_identical(plan$pipelines[[2L]]$learner, "logistic")
  defaults <- eval(formals(bi_run)$stages, list(models = NULL))
  expect_equal(defaults, c(1L, 2L))
})

test_that("unified ordered routes execute CV and report their saved genes", {
  data <- sb_example_binary(120, seed = 37)
  task <- sb_binary_task(data, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  route <- bi_pipeline(c("Stepglm", "Logistic"))
  run <- bi_run(task, pipelines = list(route = route), folds = 3, seed = 7, plots = FALSE)
  completed <- sb_complete_run(run, "route", plots = FALSE)
  expect_identical(run$comparison$comparison$status, "success")
  expect_identical(completed$fit$status, "success")
  expect_setequal(completed$genes$gene, completed$fit$selection$state$selected_features)
  expect_identical(task$data, data)
})

test_that("the home page has one modeling entry and no separate selection task", {
  file <- sb_start(view = FALSE)
  html <- paste(readLines(file, encoding = "UTF-8"), collapse = "\n")
  expect_identical(grepl('data-action="select"', html, fixed = TRUE), FALSE)
  expect_identical(grepl('id="select"', html, fixed = TRUE), FALSE)
  expect_identical(grepl("@@", html, fixed = TRUE), FALSE)
  expect_match(html, "bi_run", fixed = TRUE)
  expect_match(html, "bi_pipeline", fixed = TRUE)
})

test_that("invalid pools retain explicit failures without inventing a model", {
  alone <- sb_doctor_pipeline_set("binary", "Boruta", 1)
  expect_length(alone$pipelines, 0L)
  expect_match(alone$unavailable$Boruta, "final prediction model", fixed = TRUE)
  repeated <- tryCatch(sb_doctor_pipeline_set("binary", c("RF", "random_forest_binary_fixed"), 1), error = identity)
  expect_s3_class(repeated, "stormbreaker_bad_model_choice")
  stages <- tryCatch(sb_doctor_pipeline_set("binary", "RF", 4), error = identity)
  expect_s3_class(stages, "error")
})

test_that("feature aliases and explicit native IDs resolve to the same route", {
  expect_identical(bi_pipeline(c("ElasticNet", "Logistic")), bi_pipeline(c("Enet", "Logistic")))
  expect_identical(bi_pipeline(c("boruta_binary_roughfix", "logistic")), bi_pipeline(c("Boruta", "Logistic")))
  expect_identical(bi_pipeline("random_forest_binary_fixed"), bi_pipeline("RF"))
})

test_that("previous native admission failures remain visible in the unified pool", {
  expected <- tryCatch(bi_pipeline(c("Boruta", "cancerclass")), error = identity)
  plan <- sb_doctor_pipeline_set("binary", c("Boruta", "cancerclass"), 2)
  expect_s3_class(expected, "stormbreaker_unadmitted_pipeline")
  expect_identical(plan$unavailable[["Boruta + cancerclass"]], conditionMessage(expected))
  expect_length(plan$pipelines, 0L)
})
})

# ---- doctor-workflow ----
local({
# ---- binary-framework-alignment ----

test_that("binary comparison reuses fold-local selectors with identical results", {
  d <- sb_example_binary(120, 37)
  task <- sb_binary_task(d, "patient_id", "outcome",
    c("gene_a", "gene_b", "gene_c"))
  pipelines <- list(
    logistic = sb_pipeline("univ_binary_auc_top_k", "logistic", selector_n_features = 2L),
    nb = sb_pipeline("univ_binary_auc_top_k", "gaussian_naive_bayes", selector_n_features = 2L)
  )
  original <- sb_fit_runtime_selector_uncached
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    if (identical(selector$id, "univ_binary_auc_top_k")) calls <<- calls + 1L
    original(selector, x, endpoint)
  })
  plain <- sb_compare_cv(task, unname(pipelines), folds = 3L, seed = 19L,
    reuse_selectors = FALSE, progress = FALSE)
  plain_calls <- calls
  calls <- 0L
  shared <- sb_compare_cv(task, unname(pipelines), folds = 3L, seed = 19L,
    reuse_selectors = TRUE, progress = FALSE)
  expect_identical(shared, plain)
  expect_identical(shared$comparison$status, rep("success", 2L))
  expect_lt(calls, plain_calls)
  expect_identical(calls, 3L)
})

test_that("binary completion exports usable results and preserves cohort failures", {
  d <- sb_example_binary(120, 38)
  run <- bi_run(d, "outcome", c("gene_a", "gene_b", "gene_c"),
    "patient_id", models = c("Logistic", "NaiveBayes"), folds = 3L,
    plots = FALSE)
  external <- sb_example_binary(60, 39)
  external$patient_id <- paste0("external_", external$patient_id)
  result <- sb_complete_run(run, "Logistic",
    external = list(valid = external, incomplete = external["patient_id"]),
    plots = FALSE)
  expect_identical(result$evaluations$valid$status, "success")
  expect_identical(nrow(result$predictions$valid), 60L)
  expect_identical(result$fit$status, "success")
  expect_type(result$issues$external$incomplete, "character")
  target <- tempfile("binary-alignment-export-")
  sb_export_run(result, target)
  scores <- read.csv(file.path(target, "cohort-001-predictions.csv"))
  expect_identical(nrow(scores), 60L)
  expect_setequal(read.csv(file.path(target, "genes-by-stage.csv"))$gene,
    c("gene_a", "gene_b", "gene_c"))
})

# ---- doctor-complete ----

test_that("completion reuses the chosen fit and isolates external failures", {
  set.seed(31)
  d <- data.frame(id = paste0("d", 1:100), a = rnorm(100), b = rnorm(100))
  d$y <- rbinom(100, 1, plogis(d$a))
  run <- bi_run(d, "y", c("a", "b"), "id", models = "Logistic", folds = 3)
  e <- d
  e$id <- paste0("e", 1:100)
  saved <- run$details$Logistic$fit
  local_mocked_bindings(sb_fit_final = function(...) stop("Unexpected refit"),
    sb_compare_cv = function(...) stop("Unexpected comparison"))
  bundle <- sb_complete_run(run, "Logistic", list(valid = e, broken = e["id"]), explain = "valid")
  expect_identical(bundle$fit, saved)
  expect_identical(bundle$comparison, run$comparison$comparison)
  expect_identical(bundle$evaluations$valid$status, "success")
  expect_equal(nrow(bundle$predictions$valid), 100L)
  expect_setequal(bundle$genes$gene, c("a", "b"))
  expect_type(bundle$issues$external$broken, "character")
  expect_s3_class(bundle$plots$external$valid, "ggplot")
  expect_equal(nrow(bundle$explanation$contributions), 100L)
  path <- tempfile("workflow-export-")
  sb_export_run(bundle, path)
  expect_equal(nrow(read.csv(file.path(path, "cohort-001-predictions.csv"))), 100L)
  expect_equal(read.csv(file.path(path, "genes-by-stage.csv"))$gene, bundle$genes$gene)
  expect_gt(file.info(file.path(path, "figure-001.pdf"))$size, 0)
  expect_equal(any(grepl("broken", read.csv(file.path(path, "issues.csv"))$scope)), TRUE)
  broken_explanation <- sb_complete_run(run, "Logistic", list(valid = e, broken = e["id"]), explain = "broken")
  expect_identical(broken_explanation$evaluations$valid$status, "success")
  expect_type(broken_explanation$issues$algorithm, "character")
})

test_that("survival completion fits only the selected pipeline when plots were skipped", {
  set.seed(32)
  d <- data.frame(id = paste0("s", 1:100), a = rnorm(100), b = rnorm(100),
    time = rexp(100), event = rep(0:1, 50))
  run <- surv_run(d, "time", "event", c("a", "b"), "id", models = "Cox", folds = 3, plots = FALSE)
  e <- d
  e$id <- paste0("v", 1:100)
  local_mocked_bindings(sb_compare_cv = function(...) stop("Unexpected comparison"))
  bundle <- sb_complete_run(run, "Cox", list(validation = e), plots = TRUE)
  expect_identical(bundle$fit$status, "success")
  expect_identical(bundle$evaluations$validation$status, "success")
  expect_equal(nrow(bundle$predictions$validation), 100L)
  expect_s3_class(bundle$plots$external$validation, "patchwork")
})

# ---- doctor-custom-pipelines ----

test_that("explicit mixed-stage lists use one shared comparison", {
  for (type in c("binary", "survival")) {
    binary <- type == "binary"
    d <- if (binary) sb_example_binary(120, 37) else sb_example_survival(120, 37)
    build <- if (binary) bi_pipeline else surv_pipeline
    select <- if (binary) bi_selector else surv_selector
    method <- if (binary) "Stepglm" else "StepCox"
    model <- if (binary) "Logistic" else "Cox"
    pipelines <- list(single = build(model), pair = build(model, select(method)),
      chain = build(model, select(c(method, method))))
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")) else
      sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
    run <- if (binary) bi_run else surv_run
    result <- run(task, pipelines = pipelines, folds = 3, seed = 19, plots = FALSE)
    expected <- sb_compare_cv(task, unname(pipelines), folds = 3, seed = 19, progress = FALSE)
    expect_identical(result$pipelines, pipelines)
    expect_equal(result$comparison, expected)
    expect_length(result$comparison$results, 3L)
    expect_equal(result$comparison$comparison$status, rep("success", 3))
    completed <- sb_complete_run(result, "pair", plots = FALSE)
    expect_identical(completed$fit$status, "success")
  }
})

test_that("custom lists reject ambiguous names and incompatible task types", {
  d <- sb_example_binary(80, 11)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  p <- bi_pipeline("Logistic")
  for (pipelines in list(list(p), setNames(list(p, p), c("same", "same")),
    list(first = p, duplicate = p), list(wrong = surv_pipeline("Cox")))) {
    problem <- tryCatch(bi_run(task, pipelines = pipelines, plots = FALSE), error = identity)
    expect_s3_class(problem, "error")
  }
  problem <- tryCatch(bi_run(task, models = "Logistic", pipelines = list(one = p)), error = identity)
  expect_match(conditionMessage(problem), "not both")
})

test_that("single external risk groups are reported and exported without plotting", {
  d <- sb_example_survival(120, 37)
  run <- surv_run(d, "time", "event", c("gene_a", "gene_b", "gene_c"),
    "patient_id", models = "Cox", folds = 3, plots = FALSE)
  run$details$Cox <- list(fit = sb_fit_final(run$development_task, run$pipelines$Cox))
  e <- d
  e$patient_id <- paste0("external", seq_len(nrow(e)))
  one <- e
  for (gene in c("gene_a", "gene_b", "gene_c")) one[[gene]] <- e[[gene]][1]
  local_mocked_bindings(sb_fit_final = function(...) stop("Unexpected refit"))
  for (plots in c(FALSE, TRUE)) {
    bundle <- sb_complete_run(run, "Cox", list(one = one, two = e), plots = plots)
    expect_identical(bundle$evaluations$one$status, "success")
    expect_identical(bundle$risk_groups$one$status, "single_group")
    expect_equal(bundle$risk_groups$one$n_low + bundle$risk_groups$one$n_high, nrow(e))
    expect_identical(bundle$risk_groups$two$status, "two_groups")
    expect_null(bundle$issues$risk$two)
    expect_match(bundle$issues$risk$one, "single_group")
    expect_equal(bundle$risk_groups$one$cutoff, bundle$risk_groups$two$cutoff)
    path <- tempfile("risk-export-")
    sb_export_run(bundle, path)
    expect_match(paste(read.csv(file.path(path, "issues.csv"))$message, collapse = " "), "single_group")
    expect_equal(read.csv(file.path(path, "risk-groups.csv"))$status, c("single_group", "two_groups"))
  }
})

# ---- doctor-issues-native ----

test_that("native unavailable tables produce one issue per view without stopping runs", {
  local_mocked_bindings(sb_failure_table = function(x) sb_failure_empty())
  native <- data.frame(view = c("path", "importance"),
    status = c("unavailable", "unavailable"),
    reason = c("Original path was not captured", "No native importance"))
  x <- list(comparison = structure(list(), class = "stormbreaker_cv_comparison"),
    pipelines = list(), details = list(model = list(fit = list(status = "success"),
      figures = list(unavailable = list(native = native, extra = "No specialized plot"),
        native = list(status = native)))))
  tab <- sb_run_issues(x)
  expect_equal(nrow(tab), 3L)
  expect_setequal(tab$stage, c("plot path", "plot importance", "plot extra"))
  expect_setequal(tab$reason, c(native$reason, "No specialized plot"))
})

# ---- doctor-issues ----

test_that("issue explanations keep original reasons and distinguish stages", {
  local_mocked_bindings(sb_failure_table = function(x) sb_failure_empty())
  reasons <- c("Install missing_package", "cannot retain 1 or less variables",
    "Original path was not captured", "Predictor columns must match",
    "This nonlinear kernel has no single exact coefficient vector", "Unrecognized engine exception")
  r <- list(comparison = structure(list(), class = "stormbreaker_cv_comparison"),
    pipelines = list(), unavailable = list(), plot_errors = list(), details = list())
  for(i in seq_along(reasons))r$details[[paste0("model",i)]]<-list(
    fit=list(status="failed",error_message=reasons[i]),figures=NULL)
  tab <- sb_run_issues(r)
  expect_identical(tab$category, c("dependency","features","evidence","data","unsupported","other"))
  expect_identical(tab$reason, reasons)
  expect_identical(tab$stage, rep("final_fit",6))
  expect_identical(all(nzchar(tab$explanation)), TRUE)
  r$details <- list(ok=list(fit=list(status="success"),figures=list(unavailable=list(path=reasons[3]))))
  tab <- sb_run_issues(r)
  expect_identical(tab$stage,"plot path")
  expect_identical(r$details$ok$fit$status,"success")
})

test_that("survival doctor runs do not request ordinary binary ROC", {
  d <- sb_example_survival(90,19)
  r <- surv_run(d,time="time",event="event",id="patient_id",models="Cox",folds=3)
  expect_identical(r$comparison$status,"success")
  expect_null(r$plot_errors$roc)
  expect_null(r$plots$roc)
  expect_identical(r$issues,sb_run_issues(r))
  expect_s3_class(r$plots$performance,"ggplot")
})

# ---- doctor-workflow ----

doctor_chain_task <- function(grouped = FALSE) {
  withr::local_seed(81)
  d <- as.data.frame(matrix(rnorm(300 * 6), 300, 6))
  names(d) <- paste0("gene", 1:6)
  d$id <- paste0("p", 1:300)
  d$y <- rbinom(300, 1, plogis(2 * d$gene1 - 2 * d$gene2))
  if (grouped) {
    d$patient <- d$id
    d <- d[rep(seq_len(nrow(d)), each = 2), ]
    d$id <- paste0("r", seq_len(nrow(d)))
  }
  sb_binary_task(d, "id", "y", paste0("gene", 1:6), group_col = if (grouped) "patient" else NULL)
}

test_that("short menus cover all existing binary and survival learners", {
  for (task in c("binary", "survival")) {
    menu <- sb_models(task)
    predictors <- if (task == "binary") menu$model_id[menu$can_predict] else menu$model_id
    expect_setequal(predictors, sb_model_catalog(task)$model_id)
    expect_identical(anyDuplicated(tolower(menu$name)), 0L)
    plan <- sb_doctor_pipeline_set(task, if (task == "binary") c("RF", "GBM") else c("Cox", "CoxBoost"), c(1, 3))
    expect_length(plan$unavailable, 0L)
    expect_gte(length(plan$pipelines), 3L)
    for (p in plan$pipelines) expect_identical(sb_lint_pipeline(p)$status, "admitted")
  }
})

test_that("ordered binary chains keep stage evidence and frozen predictions", {
  task <- doctor_chain_task()
  p <- bi_pipeline("GBM", bi_selector(c("Lasso", "RF")))
  f <- sb_fit_final(task, p, capture_diagnostics = TRUE)
  expect_identical(f$status, "success")
  stages <- f$selection$state$stages
  expect_length(stages, 2L)
  expect_identical(stages[[2]]$source_features, stages[[1]]$selected_features)
  copied <- unserialize(serialize(f, NULL))
  new <- task$data[c("id", task$predictors)]
  expect_equal(sb_predict_final(copied, new, "id")$predictions,
    sb_predict_final(f, new, "id")$predictions)
  z <- sb_plots(f)
  expect_identical(all(c("selector_1", "selector_2") %in% names(z$plots)), TRUE)
  cv <- sb_run_cv(task, p, folds = 3, seed = 7)
  expect_identical(cv$status, "success")
  reverse <- bi_pipeline("GBM", bi_selector(c("RF", "Lasso")))
  expect_identical(p$id == reverse$id, FALSE)
})

test_that("chains refit within glmnet and SVM inner folds with grouped patients", {
  task <- doctor_chain_task(TRUE)
  original <- sb_fit_binary_model_selector
  sizes <- integer()
  local_mocked_bindings(sb_fit_binary_model_selector = function(x, outcome, id, params, group_ids = NULL) {
    result <- original(x, outcome, id, params, group_ids)
    sizes <<- c(sizes, nrow(x))
    if (id == "lasso_binary_nonzero") {
      expect_length(group_ids, nrow(x))
      expect_identical(all(vapply(split(result$evidence$fold_id, group_ids),
        function(v) length(unique(v)) == 1L, logical(1))), TRUE)
    }
    result
  })
  for (model in c("Lasso", "SVM")) {
    sizes <- integer()
    f <- sb_fit_final(task, bi_pipeline(model, bi_selector(c("Lasso", "RF"))))
    expect_identical(f$status, "success")
    expect_gt(length(sizes), 2L)
    expect_identical(any(sizes < task$n), TRUE)
  }
})

test_that("survival selector lists use saved stages in final and CV paths", {
  d <- sb_example_survival(120, 37)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
  p <- surv_pipeline("Cox", surv_selector(c("StepCox", "StepCox")))
  f <- sb_fit_final(task, p)
  expect_identical(f$status, "success")
  expect_length(f$selection$state$stages, 2)
  expect_identical(sb_run_cv(task, p, folds = 3)$status, "success")
})

test_that("doctor runs produce small-model figures without a browser", {
  d <- sb_example_binary(80, 11)
  r <- bi_run(d, "outcome", id = "patient_id", models = c("NearestCentroid", "BernoulliNB"), stages = 1, folds = 3, plots = TRUE)
  expect_length(r$details, 2)
  for (detail in r$details) {
    expect_identical(detail$fit$status, "success")
    expect_length(detail$figures$plots, 1)
    expect_s3_class(detail$figures$plots[[1]], "ggplot")
  }
  expect_length(sb_doctor_pipeline_set("binary", "RF", 2)$pipelines, 1L)
})

test_that("multiclass short entrypoint returns CV and algorithm figures", {
  d <- iris
  r <- multi_run(d, "Species", models = "RF", folds = 3)
  expect_identical(r$details$status$status, "success")
  expect_identical(r$comparison$comparison$status, "success")
  expect_gt(r$details$status$n_plots, 0)
})

# ---- lean-workflow ----

test_that("batch defaults defer final fits but preserve the comparison", {
  d <- sb_example_binary(70, 19)
  args <- list(data = d, outcome = "outcome", predictors = c("gene_a", "gene_b"),
    id = "patient_id", models = c("Logistic", "NaiveBayes"), folds = 3, seed = 11)
  eager <- do.call(bi_run, c(args, list(plots = TRUE)))
  fit_original <- sb_fit_final
  calls <- 0L
  local_mocked_bindings(sb_fit_final = function(...) {
    calls <<- calls + 1L
    fit_original(...)
  })
  deferred <- do.call(bi_run, args)
  expect_identical(calls, 0L)
  expect_length(deferred$details, 0L)
  expect_gt(length(deferred$plots), 0L)
  expect_equal(deferred$comparison$comparison, eager$comparison$comparison)
  local_mocked_bindings(sb_compare_cv = function(...) stop("Comparison must not repeat"))
  chosen <- sb_complete_run(deferred, "Logistic")
  expect_identical(calls, 1L)
  expect_equal(chosen$genes$gene, args$predictors)
  expect_equal(sb_predict_final(chosen$fit, d, "patient_id")$predictions,
    sb_predict_final(eager$details$Logistic$fit, d, "patient_id")$predictions)
  expect_setequal(names(chosen$plots$algorithm), names(eager$details$Logistic$figures$plots))
})

test_that("calculation entrypoints do not request browsers by default", {
  expect_identical(formals(sb_compare_models)$view, FALSE)
  expect_identical(formals(sb_compare_cohorts)$view, FALSE)
})

test_that("survival batch defaults preserve results and defer selection to completion", {
  d <- sb_example_survival(80, 21)
  args <- list(data = d, time = "time", event = "event", predictors = c("gene_a", "gene_b"),
    id = "patient_id", models = c("Cox", "survreg_weibull"), folds = 3, seed = 11)
  eager <- do.call(surv_run, c(args, list(plots = TRUE)))
  original <- sb_fit_final
  calls <- 0L
  local_mocked_bindings(sb_fit_final = function(...) {
    calls <<- calls + 1L
    original(...)
  })
  run <- do.call(surv_run, args)
  expect_identical(calls, 0L)
  expect_length(run$details, 0L)
  expect_equal(run$comparison$comparison, eager$comparison$comparison)
  expect_gt(length(run$plots), 0L)
  local_mocked_bindings(sb_compare_cv = function(...) stop("Unexpected comparison"))
  done <- sb_complete_run(run, "Cox")
  expect_identical(calls, 1L)
  expect_equal(sb_predict_final(done$fit, d, "patient_id")$predictions,
    sb_predict_final(eager$details$Cox$fit, d, "patient_id")$predictions)
  expect_setequal(done$genes$gene, args$predictors)
  expect_setequal(names(done$plots$algorithm), names(eager$details$Cox$figures$plots))
})

test_that("single pipeline still draws details and explicit FALSE skips them", {
  d <- sb_example_binary(60, 21)
  args <- list(data = d, outcome = "outcome", predictors = c("gene_a", "gene_b"),
    id = "patient_id", models = "Logistic", folds = 3)
  one <- do.call(bi_run, args)
  expect_identical(one$details$Logistic$fit$status, "success")
  expect_gt(length(one$details$Logistic$figures$plots), 0L)
  none <- do.call(bi_run, c(args, list(plots = FALSE)))
  expect_length(none$details, 0L)
  expect_length(none$plots, 0L)
  expect_equal(one$comparison$comparison, none$comparison$comparison)
})

test_that("multiclass batch comparisons avoid extra method fits", {
  local_mocked_bindings(sb_analyze_methods = function(...) stop("Unexpected final fits"))
  out <- multi_run(iris, "Species", models = c("lda", "qda"), folds = 3)
  expect_null(out$details)
  expect_identical(out$comparison$comparison$status, c("success", "success"))
})


test_that("default comparisons never open pages but explicit requests still do", {
  d <- sb_example_binary(60, 13)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  opened <- character()
  local_mocked_bindings(sb_open_result_page = function(path) opened <<- c(opened, path))
  x <- sb_compare_models(task, models = "logistic", folds = 3)
  d$patient_id <- paste0("external_", d$patient_id)
  external <- list(validation = sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b")))
  y <- sb_compare_cohorts(task, external, models = "logistic", folds = 3)
  expect_length(opened, 0L)
  expect_null(x$report_file)
  expect_null(y$report_file)
  visible <- sb_compare_cohorts(task, external, models = "logistic", folds = 3, view = TRUE)
  expect_length(opened, 1L)
  expect_equal(visible$development$comparison, y$development$comparison)
})

test_that("deferred screening retains selected genes and selector figures", {
  d <- sb_example_binary(200, 21)
  withr::local_seed(81)
  d$outcome <- rbinom(nrow(d), 1, plogis(2 * d$gene_a - 2 * d$gene_b))
  args <- list(data = d, outcome = "outcome", id = "patient_id",
    predictors = c("gene_a", "gene_b", "gene_c"),
    models = c("Logistic", "NaiveBayes"), selector = bi_selector("Lasso"), folds = 3)
  eager <- do.call(bi_run, c(args, list(plots = TRUE)))
  run <- do.call(bi_run, args)
  name <- names(run$pipelines)[1L]
  done <- sb_complete_run(run, name)
  reference <- sb_complete_run(eager, name)
  expect_identical(done$fit$status, "success")
  expect_equal(done$genes, reference$genes)
  expect_equal(sb_predict_final(done$fit, d, "patient_id")$predictions,
    sb_predict_final(reference$fit, d, "patient_id")$predictions)
  expect_setequal(names(done$plots$algorithm), names(reference$plots$algorithm))
})

# ---- sharing-entrypoints ----

test_that("IECV doctor and stack final loops activate bounded model sharing", {
  skip_if_not_installed("survivalsvm")
  skip_if_not_installed("glmnet")
  d <- sb_example_survival(90L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  pp <- list(first = sb_pipeline("none", "survivalsvm"),
    second = sb_pipeline("fixed_features", "survivalsvm",
      selector_features = features, feature_set_id = "all"))
  native <- sb_fit_survivalsvm
  calls <- 0L
  final_counts <- integer()
  final <- sb_fit_final
  scope <- sb_with_selector_reuse
  reuse_enabled <- TRUE
  local_mocked_bindings(sb_fit_survivalsvm = function(...) {
    calls <<- calls + 1L; native(...)
  }, sb_fit_final = function(...) {
    before <- calls
    value <- final(...)
    final_counts <<- c(final_counts, calls - before)
    value
  }, sb_plots = function(...) list(),
  sb_with_selector_reuse = function(expr, enabled = TRUE, diagnostics = FALSE) {
    scope(expr, enabled = enabled && reuse_enabled, diagnostics = diagnostics)
  })
  for (entry in c("iecv", "doctor", "stack")) {
    run <- function() switch(entry,
      iecv = sb_fit_iecv_strategy(task, sb_validate_iecv_strategy(pp, "survival")),
      doctor = surv_run(task, pipelines = pp, folds = 3L, seed = 21L, plots = TRUE),
      stack = sb_fit_stack_final(task, sb_stack(pp, inner_folds = 3L, inner_seed = 21L)))
    reuse_enabled <- FALSE
    final_counts <- integer()
    plain <- run()
    expect_identical(final_counts, c(1L, 1L), info = entry)
    reuse_enabled <- TRUE
    final_counts <- integer()
    value <- run()
    expect_equal(value, plain, info = entry)
    fits <- switch(entry, iecv = value, doctor = lapply(value$details, `[[`, "fit"), stack = value$member_fits)
    expect_identical(vapply(fits, `[[`, character(1), "status"),
      stats::setNames(rep("success", 2L), names(fits)), info = entry)
    expect_identical(final_counts, c(1L, 0L), info = entry)
    expect_null(.sb_selector_reuse_context$active)
  }
})

# ---- start ----

test_that("start page contains installed recipes and never runs them", {
  local_mocked_bindings(sb_fit_final = function(...) stop("Must not fit"),
    sb_compare_models = function(...) stop("Must not compare"),
    sb_analyze_methods = function(...) stop("Must not analyze"))
  file <- sb_start(view = FALSE)
  html <- paste(readLines(file, encoding = "UTF-8"), collapse = "\n")
  expect_identical(file.exists(file), TRUE)
  expect_identical(grepl("@@", html, fixed = TRUE), FALSE)
  for (action in c("compare", "validate")) {
    expect_match(html, paste0('data-action="', action, '"'), fixed = TRUE)
    expect_match(html, paste0('id="', action, '-code"'), fixed = TRUE)
  }
  expect_match(html, "eval_external(", fixed = TRUE)
  expect_identical(grepl('src="https://', html, fixed = TRUE), FALSE)
})

test_that("both doctor examples execute through real interfaces", {
  for (p in c("glmnet", "randomForest", "ggplot2")) skip_if_not_installed(p)
  withr::local_pdf(tempfile(fileext = ".pdf"))
  results <- lapply(c("compare", "validate"), function(action) {
    e <- new.env(parent = globalenv())
    sys.source(system.file("examples", "advanced-examples.R",
      package = "stormbreaker", mustWork = TRUE), envir = e)
    e$run_advanced_example(paste0("doctor-", action))
  })
  comparison <- results[[1L]]$comparison
  expect_length(results[[1L]]$run$pipelines, 3L)
  expect_identical(comparison$comparison$status, rep("success", 3))
  expect_identical(results[[1L]]$views$status$view, c("ROC", "Calibration", "Decision curve"))
  expect_identical(results[[2L]]$fit$status, "success")
  expect_identical(results[[2L]]$prediction$status, "success")
  expect_identical(results[[2L]]$evaluation$status, "success")
})

test_that("selection example finalizes its explicitly declared successful route", {
  skip_if_not_installed("glmnet")
  e <- new.env(parent = globalenv())
  sys.source(system.file("examples", "advanced-examples.R",
    package = "stormbreaker", mustWork = TRUE), envir = e)
  result <- e$run_advanced_example("doctor-select")
  expect_identical(names(result$r$pipelines), "Lasso + Logistic")
  expect_identical(result$r$comparison$comparison$status, "success")
  expect_identical(result$finished$fit$status, "success")
})

# ---- user-workflow ----

workflow_test_environment <- function() {
  env <- new.env(parent = globalenv())
  source(system.file("examples", "advanced-examples.R", package = "stormbreaker"),
    local = env)
  env
}

test_that("CSV helper preserves IDs and requires explicit model confirmation", {
  helper <- workflow_test_environment()
  directory <- withr::local_tempdir()
  path <- file.path(directory, "id.csv")
  writeLines(c("patient_id,gene_a", "00123,2.5"), path)
  expect_identical(helper$workflow_read_csv(path)$patient_id, "00123")
  expect_error(helper$workflow_finalize("does-not-exist.rds"), "Explicit")
  expect_error(helper$workflow_new_directory(directory), "new output directory")
})

test_that("both workflow templates complete outcome-blind prediction and evaluation", {
  helper <- workflow_test_environment()
  directory <- withr::local_tempdir()
  input <- file.path(directory, "input")
  helper$workflow_make_templates(input)
  expect_length(list.files(input), 6L)
  for (kind in c("binary", "survival")) {
    new_data <- helper$workflow_read_csv(file.path(input,
      paste0(kind, "_external_predictors.csv")))
    expect_false(any(c("outcome", "time", "event") %in% names(new_data)))
    development_dir <- file.path(directory, paste0(kind, "_development"))
    state <- suppressMessages(helper$workflow_develop(file.path(input,
      paste0(kind, "_development.csv")), development_dir, kind))
    expect_true(state$suggested_candidate %in% names(state$pipelines))
    expect_equal(nrow(state$comparison), 2L)
    for (stem in c("gene_selection", "gene_nonzero")) {
      ranking_file <- file.path(development_dir, paste0(stem, "_ranking.csv"))
      expect_identical(file.exists(ranking_file), TRUE)
      ranking <- read.csv(ranking_file)
      expect_setequal(ranking$feature, state$predictors)
      expect_identical(file.exists(file.path(development_dir,
        paste0(stem, "_counting_groups.csv"))), TRUE)
    }
    expect_true(all(state$comparison$eligible_for_primary_rule))
    expect_true(all(vapply(state$cv$results, function(x)
      identical(x$provenance$fold_assignment, state$cv$fold_assignment), logical(1))))
    final_dir <- file.path(directory, paste0(kind, "_final"))
    result <- helper$workflow_finalize(file.path(development_dir, "development.rds"),
      state$suggested_candidate, "Synthetic test of the declared rule.",
      file.path(input, paste0(kind, "_external_predictors.csv")),
      file.path(input, paste0(kind, "_external_outcomes.csv")), final_dir,
      render_report = FALSE)
    expect_identical(result$fit$status, "success")
    expect_identical(result$prediction$status, "success")
    expect_identical(result$evaluation$status, "success")
    expect_equal(nrow(result$prediction$predictions), 40L)
    reloaded <- readRDS(file.path(final_dir, "final_model.rds"))
    expect_identical(sb_predict_final(reloaded, new_data, "patient_id"), result$prediction)
    expect_true(all(file.exists(file.path(final_dir,
      c("selection_decision.rds", "final_model.rds", "external_predictions.csv",
        "external_predictions.rds", "external_evaluation.rds")))))
    # A missing outcome file fails only after outcome-blind predictions are saved.
    blind_dir <- file.path(directory, paste0(kind, "_blind"))
    expect_error(suppressWarnings(helper$workflow_finalize(
      file.path(development_dir, "development.rds"), state$suggested_candidate,
      "Synthetic outcome isolation test.",
      file.path(input, paste0(kind, "_external_predictors.csv")),
      file.path(directory, "absent_outcomes.csv"), blind_dir,
      render_report = FALSE)))
    expect_identical(readRDS(file.path(blind_dir, "external_predictions.rds")),
      result$prediction)
  }
})
})

# ---- doctor-presentation ----
local({
# ---- clinician-explanations ----

test_that("clinical reports reuse exact saved explanations and declared dependence features", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_binary(n = 80L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  pipeline <- sb_pipeline(learner = "logistic")
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 41L)
  fit <- sb_fit_final(task, pipeline)
  shap <- sb_additive_shap(fit, d, "patient_id", d)
  nomogram <- sb_nomogram(fit, d)
  before <- serialize(list(fit, shap, nomogram), NULL, version = 2L)
  values <- d[rev(seq_len(nrow(d))), ]
  views <- sb_clinician_views(result, shap = shap, feature_values = values,
    id_col = "patient_id", dependence_features = "gene_b", nomogram = nomogram)
  expect_identical(serialize(list(fit, shap, nomogram), NULL, version = 2L), before)
  expect_identical(views$plots[["SHAP importance"]]$data, sb_plot_shap_importance(shap)$data)
  expect_identical(views$plots[["SHAP beeswarm"]]$data,
    sb_plot_shap_beeswarm(shap, feature_values = values, id_col = "patient_id")$data)
  expect_identical(views$plots[["SHAP dependence: gene_b"]]$data,
    sb_plot_shap_dependence(shap, "gene_b", values, "patient_id")$data)
  expect_identical(views$plots[["Original-model nomogram"]]$data, sb_plot_nomogram(nomogram)$data)
  expect_null(views$plots[["SHAP dependence: gene_a"]])
  expect_length(views$height, length(views$plots))
  expect_length(views$width, length(views$plots))
  expect_true(all(views$status$status[views$status$group == "Model explanations"] == "available"))
  file <- tempfile(fileext = ".html")
  sb_render_clinician_report(views, file, complete = TRUE)
  html <- paste(readLines(file), collapse = "\n")
  expect_true(grepl("Original-model nomogram", html, fixed = TRUE))
  expect_true(grepl("SHAP dependence: gene_b", html, fixed = TRUE))
  expect_false(grepl("patient card", html, fixed = TRUE))
  broken <- values
  broken$patient_id[1L] <- "wrong-id"
  failed <- sb_clinician_views(result, shap = shap, feature_values = broken,
    id_col = "patient_id", dependence_features = "gene_b")
  expect_identical(failed$status$status[failed$status$view == "SHAP dependence: gene_b"], "error")
  expect_null(failed$plots[["SHAP dependence: gene_b"]])
  expect_error(sb_clinician_views(result, shap = shap, dependence_features = "gene_b"),
    class = "stormbreaker_bad_report_input")
  expect_error(sb_clinician_views(result, nomogram = fit), class = "stormbreaker_bad_report_input")
})

# ---- clinician-labels ----

test_that("clinical labels hide internal status without changing model grouping", {
  skip_if_not_installed("ggplot2")
  fixture <- make_sb_curve_view_fixture("binary", comparison = TRUE)
  views <- sb_clinician_views(fixture$cv)
  clean <- sb_clinician_labeler(views$models)
  expect_identical(clean("Unnamed cohort\nRole: pipeline\nAssessment fold 2\ncurve: success run: success\nSaved entry: root/1"), "Fold 2")
  expect_identical(views$models$model, c("Logistic regression", "Gaussian naive Bayes"))
  expect_identical(clean(views$models$id), c("Logistic", "Gaussian NB"))
  expect_identical(grepl("smooth|none__", paste(clean(views$models$id), collapse = " ")), FALSE)
  expect_null(views$plots[["Confusion matrix"]])
  for (name in c("ROC", "Calibration", "Decision curve")) {
    p <- views$plots[[name]]
    rendered <- unlist(p$facet$params$labeller(data.frame(panel = unique(p$data$panel))))
    expect_identical(any(grepl("Unnamed cohort|Role:|curve: success|run: success|Saved entry:|none__", rendered)), FALSE)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  }
  expect_identical(views$plots$ROC$data, sb_plot_roc(fixture$cv)$data)
})

test_that("same-family configurations get distinct readable labels", {
  leaf <- function(id) list(value = list(task_type = "binary",
    pipeline = list(id = id, learner = "logistic", selector = "none")))
  models <- sb_clinician_models(list(leaf("one"), leaf("two")))
  expect_length(unique(models$model), 2L)
  expect_identical(models$model, c("Logistic regression (1)", "Logistic regression (2)"))
  expect_identical(sb_clinician_labeler(models)(models$id), c("Logistic (1)", "Logistic (2)"))
})

# ---- clinician-simplification ----

test_that("ordinary model comparison defaults to three declared models", {
  local_mocked_bindings(
    sb_compare_cv = function(task, pipelines, ...) list(comparison = data.frame(
      pipeline_id = vapply(pipelines, function(p) p$id, character(1)))),
    sb_attach_result_page = function(x, ...) x)
  for (type in c("binary", "survival")) {
    d <- if (type == "binary") sb_example_binary(60) else sb_example_survival(60)
    task <- if (type == "binary") sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b")) else
      sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month")
    expected <- if (type == "binary") c("logistic", "glmnet_logistic", "random_forest_binary") else
      c("coxph", "glmnet_cox", "ranger_survival")
    result <- sb_compare_models(task, view = FALSE)
    expect_identical(result$requested_models, expected)
    expect_identical(sb_compare_models(task, models = NULL, view = FALSE)$requested_models,
      sb_model_catalog(type)$model_id)
  }
  expect_identical(formals(sb_compare_multiclass)$methods,
    quote(c("lasso", "random_forest", "svm_radial")))
})

test_that("brief clinical views do not compute optional plots", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("binary")$cv
  before <- serialize(x, NULL)
  full <- sb_clinician_views(x, detail = "full")
  local_mocked_bindings(sb_plot_pr = function(...) stop("Not requested"),
    sb_probability_distribution = function(...) stop("Not requested"))
  brief <- sb_clinician_views(x)
  expect_identical(brief$status$view, c("ROC", "Calibration", "Decision curve"))
  expect_identical(brief$plots$ROC$data, full$plots$ROC$data)
  expect_identical(serialize(x, NULL), before)
  expect_identical(all(c("Precision-recall", "Probability distribution") %in% full$status$view), TRUE)
})

# ---- clinician-views ----

test_that("clinical assembly preserves saved DCA and native data without fitting", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_binary(n = 100L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"),
    decision_thresholds = c(.1, .3, .6))
  pipeline <- sb_pipeline(learner = "gaussian_naive_bayes")
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 41L)
  expect_identical(result$status, "success")
  fit <- sb_fit_final(task, pipeline)
  native <- sb_naive_bayes_diagnostics(fit)
  before <- serialize(result, NULL, version = 2L)
  views <- sb_clinician_views(result, c(.2, .5, .8), .5, native = list(Density = native), detail = "full")
  expect_identical(views$plots[["Decision curve"]]$data, sb_plot_dca(result)$data)
  expect_identical(views$plots[["Model: Density"]]$data, sb_plot_algorithm_diagnostic(native)$data)
  for (plot in views$plots) {
    expect_null(plot$labels$subtitle)
    expect_null(plot$labels$caption)
  }
  expect_equal(sort(unique(views$plots[["Threshold performance"]]$data$threshold)), c(.2, .5, .8))
  expect_identical(serialize(result, NULL, version = 2L), before)
  expect_equal(nrow(views$status), 8L)
  expect_false(any(views$status$status == "error"))
  expect_true(all(c("Clinical decisions", "Model diagnostics") %in% views$groups))
  expect_true(anyNA(views$plots[["Threshold performance"]]$data$value))
  file <- tempfile(fileext = ".html")
  sb_render_clinician_report(views, file)
  html <- paste(readLines(file), collapse = "\n")
  expect_identical(grepl("View availability", html, fixed = TRUE), FALSE)
  expect_match(html, "Model diagnostics", fixed = TRUE)
  expect_true(grepl("original Gaussian means", html, fixed = TRUE))
  expect_error(sb_clinician_views(result, .2, .5), class = "stormbreaker_bad_report_input")
  fit_before <- serialize(fit, NULL, version = 2L)
  local_mocked_bindings(sb_fit_final = function(...) stop("Report must not fit models"))
  requested <- sb_clinician_views(result, native = sb_native_diagnostics(fit))
  model_plots <- requested$plots[requested$groups == "Model diagnostics"]
  expect_length(model_plots, 1L)
  expect_identical(model_plots[[1L]]$data, sb_plot_algorithm_diagnostic(native)$data)
  detailed <- sb_clinician_views(result, native = sb_native_diagnostics(fit, complete = TRUE))
  expect_equal(sum(detailed$groups == "Model diagnostics"), 2L)
  expect_identical(serialize(fit, NULL, version = 2L), fit_before)
  expect_identical(serialize(result, NULL, version = 2L), before)
  result$predictions$probability[1L] <- NA_real_
  partial <- sb_clinician_views(result)
  expect_identical(partial$status$status[partial$status$view == "ROC"], "unavailable")
  expect_match(partial$status$reason[partial$status$view == "ROC"], "invalid_predictions", fixed = TRUE)
})

test_that("clinical figures move explanations to captions without changing their evidence", {
  skip_if_not_installed("ggplot2")
  result <- make_sb_curve_view_fixture("survival")$cv
  original <- sb_plot_calibration(result) + ggplot2::labs(
    subtitle = "Saved cohort <A>", caption = "Intervals use saved estimates.")
  local_mocked_bindings(sb_plot_calibration = function(...) original)
  views <- sb_clinician_views(result)
  plot <- views$plots[["Calibration"]]
  expect_identical(plot$data, original$data)
  expect_identical(plot$layers, original$layers)
  expect_identical(plot$labels$title, original$labels$title)
  expect_identical(views$captions[match("Calibration", names(views$plots))],
    "Saved cohort <A>\nIntervals use saved estimates.")
  expect_null(plot$labels$subtitle)
  expect_null(plot$labels$caption)
  file <- tempfile(fileext = ".html")
  sb_render_clinician_report(views, file)
  html <- paste(readLines(file), collapse = "\n")
  expect_match(html, "Saved cohort &lt;A&gt;", fixed = TRUE)
  expect_match(html, "Intervals use saved estimates.", fixed = TRUE)
})

test_that("explicit query diagnostics can be added without automatic extraction", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_binary(n = 60L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "knn_binary"))
  result <- make_sb_curve_view_fixture("binary")$cv
  missing <- sb_clinician_views(result, native = sb_native_diagnostics(fit))
  expect_identical(missing$status$status[missing$status$group == "Model diagnostics"], "unavailable")
  query <- d[c(2L, 9L), ]
  supplied <- sb_clinician_views(result, native = sb_native_diagnostics(fit, newdata = query))
  plot <- supplied$plots[supplied$groups == "Model diagnostics"][[1L]]
  expect_identical(plot$data, sb_plot_algorithm_diagnostic(sb_knn_neighborhood(fit, query))$data)
})

test_that("default and expanded reports never extract or plot technical diagnostics", {
  skip_if_not_installed("ggplot2")
  result <- make_sb_curve_view_fixture("binary")$cv
  unexpected <- function(...) stop("Technical diagnostics were not requested")
  local_mocked_bindings(sb_native_diagnostics = unexpected,
    sb_selector_diagnostics = unexpected, sb_plot_algorithm_diagnostic = unexpected,
    sb_plot_selector_diagnostic = unexpected)
  views <- sb_clinician_views(result)
  expect_equal(sum(views$groups %in% c("Model diagnostics", "Feature selection")), 0L)
  for (complete in c(FALSE, TRUE)) {
    file <- tempfile(fileext = ".html")
    sb_render_clinician_report(views, file, complete = complete)
    html <- paste(readLines(file), collapse = "\n")
    expect_identical(grepl("Model diagnostics|Feature selection", html), FALSE)
  }
})

test_that("survival assembly does not reinterpret censoring as binary outcomes", {
  skip_if_not_installed("ggplot2")
  result <- make_sb_curve_view_fixture("survival")$cv
  views <- sb_clinician_views(result)
  expect_identical(views$status$status[views$status$view == "ROC"], "unavailable")
  expect_null(views$plots[["Probability distribution"]])
  expect_identical(views$plots[["Decision curve"]]$data, sb_plot_dca(result)$data)
  expect_identical(views$plots[["Calibration"]]$data, sb_plot_calibration(result)$data)
})

test_that("missing clinical curves remain explicit and are never synthesized", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_binary(n = 80L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  result <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3L, seed = 41L)
  views <- sb_clinician_views(result)
  expect_identical(views$status$status[views$status$view == "Decision curve"], "unavailable")
  expect_null(views$plots[["Decision curve"]])
  expect_true(length(views$diagnostics[["Decision curve"]]) > 0L)
})

# ---- composite-result-displays ----

test_that("repeated displays and exports preserve separate identities and failures", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture()$cv
  failed <- x
  failed$status <- "failed"
  failed$runs$status[1] <- "failed"
  failed$runs$error_message[1] <- "original, \"reason\"\nnext line"
  repeated <- structure(list(repetitions = list(x, failed)),
                        class = "stormbreaker_repeated_cv_comparison")
  before <- serialize(repeated, NULL)
  p <- sb_plot_performance(repeated, "roc_auc", show_unavailable = TRUE)
  expect_length(levels(p$data$display_label), 2)
  expect_equal(nrow(p$layers[[1]]$data), 1)
  for (fn in list(sb_plot_calibration, sb_plot_dca)) {
    p <- fn(repeated, show_unavailable = TRUE)
    status <- attr(p, "stormbreaker_curve_status")
    expect_equal(status$repeat_id, 1:2)
    expect_length(unique(status$panel), 2)
    expect_silent(ggplot2::ggplot_build(p))
  }
  failures <- sb_failure_table(repeated)
  expect_true(all(failures$repeat_id == 2L))
  expect_true(any(failures$error_message == failed$runs$error_message[1], na.rm = TRUE))
  mixed <- sb_failure_table(list(single = x, repeated = repeated))
  expect_equal(nrow(mixed), nrow(failures))
  out <- withr::local_tempdir()
  html_path <- file.path(out, "report.html")
  sb_render_report(list(repeated = repeated, ordinary = x), html_path)
  html <- paste(readLines(html_path), collapse = "\n")
  expect_true(grepl("Development repeated CV", html, fixed = TRUE))
  csv <- read.csv(file.path(out, "report_failures.csv"))
  expect_equal(csv$repeat_id, failures$repeat_id)
  expect_equal(csv$error_message, failures$error_message)
  expect_error(sb_render_report(repeated, html_path), class = "stormbreaker_report_exists")
  expect_identical(serialize(repeated, NULL), before)
})

test_that("IECV reports and curves label development validation without pooling", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("survival")$external
  iecv <- structure(list(folds = list(
    list(heldout_cohort_id = "Site <A>", comparison = x),
    list(heldout_cohort_id = "Site B", comparison = x))),
    class = "stormbreaker_internal_external_cv")
  for (fn in list(sb_plot_calibration, sb_plot_dca)) {
    p <- fn(iecv, show_unavailable = TRUE)
    status <- attr(p, "stormbreaker_curve_status")
    expect_equal(status$heldout_cohort_id, c("Site <A>", "Site B"))
    expect_true(all(status$validation == "internal_external_cv"))
    expect_true(all(grepl("not final independent", status$panel, fixed = TRUE)))
    expect_silent(ggplot2::ggplot_build(p))
  }
  path <- file.path(withr::local_tempdir(), "iecv.html")
  sb_render_report(iecv, path)
  html <- paste(readLines(path), collapse = "\n")
  expect_true(grepl("Development IECV (not final independent external validation)", html, fixed = TRUE))
  expect_true(grepl("Site &lt;A&gt;", html, fixed = TRUE))
  expect_true(grepl("heldout_cohort_id", html, fixed = TRUE))
  expect_equal(nrow(sb_failure_table(iecv)), 0)
})

test_that("composite strip labels wrap long identifiers without losing characters", {
  labels <- c(paste(rep("long_model_identifier", 8), collapse = "__"),
              "Held-out cohort A\nDevelopment IECV")
  wrapped <- sb_view_wrap_label(labels)
  expect_true(all(nchar(unlist(strsplit(wrapped, "\n", fixed = TRUE))) <= 50L))
  expect_equal(gsub("\n", "", wrapped, fixed = TRUE), gsub("\n", "", labels, fixed = TRUE))
})

# ---- concise-plot-style ----

test_that("shared figure themes keep explanatory prose outside the drawing", {
  skip_if_not_installed("ggplot2")

  expect_s3_class(sb_viz_theme()$plot.caption, "element_blank")
  expect_s3_class(sb_viz_theme()$plot.subtitle, "element_blank")
  expect_s3_class(sb_shap_theme()$plot.caption, "element_blank")
  expect_s3_class(sb_shap_theme()$plot.subtitle, "element_blank")
  expect_equal(sb_viz_theme()$text$size, 11)
})

test_that("compacting keeps explanations and numeric plot content", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(x = c(1, 2), y = c(.2, .8))
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y)) + ggplot2::geom_point() +
    ggplot2::labs(title = "Saved model coefficients", x = "Time (month)",
      subtitle = "Reference population", caption = "Conditional association")
  attr(p, "stormbreaker_notes") <- "Original interpretation"
  before <- ggplot2::ggplot_build(p)$data

  q <- sb_compact_plot(p)

  expect_identical(q$labels$title, "Model coefficients")
  expect_null(q$labels$subtitle)
  expect_null(q$labels$caption)
  expect_setequal(attr(q, "stormbreaker_notes"),
    c("Original interpretation", "Reference population", "Conditional association"))
  expect_identical(q$data, d)
  expect_identical(q$labels$x, "Time (month)")
  expect_identical(ggplot2::ggplot_build(q)$data, before)
})

test_that("short clinical labels retain cohort and fold distinctions", {
  models <- data.frame(id = "model_id", model = "Logistic regression")
  clean <- sb_clinician_labeler(models)

  expect_identical(clean("model_id"), "Logistic")
  expect_identical(clean("Development CV\nOut-of-fold"), "CV")
  expect_identical(clean("Validation A\nAssessment fold 2"), "Validation A\nFold 2")
  expect_identical(clean("Validation B\nAssessment fold 2"), "Validation B\nFold 2")
})

test_that("ROC AUC labels use colour instead of repeating model names", {
  skip_if_not_installed("ggplot2")
  run <- bi_run(sb_example_binary(80, seed = 42), "outcome",
    c("gene_a", "gene_b"), "patient_id", models = c("Logistic", "NaiveBayes"),
    folds = 3, plots = FALSE)
  original <- lapply(run$comparison$results, `[[`, "predictions")
  local_mocked_bindings(sb_fit_final = function(...) stop("Unexpected refit"))

  p <- sb_plot_roc(run)

  statistics <- attr(p, "stormbreaker_roc_statistics")
  expect_identical(statistics$label, paste0("AUC ", sprintf("%.3f", statistics$auc)))
  expect_identical(all.vars(p$layers[[3]]$mapping$colour), "pipeline_id")
  expect_identical(lapply(run$comparison$results, `[[`, "predictions"), original)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  views <- sb_clinician_views(run)
  expect_match(views$captions[match("ROC", names(views$plots))],
    "AUC labels describe the displayed curve", fixed = TRUE)
})

test_that("performance axes name the metric without instruction text", {
  skip_if_not_installed("ggplot2")
  run <- bi_run(sb_example_binary(80, seed = 42), "outcome",
    c("gene_a", "gene_b"), "patient_id", models = "Logistic", folds = 3,
    plots = FALSE)

  p <- sb_plot_performance(run$comparison, "roc_auc")

  expect_identical(p$labels$y, "ROC AUC")
  expect_identical(p$scales$get_scales("x")$get_labels(), "Logistic")
  expect_s3_class(p$theme$plot.caption, "element_blank")
  expect_s3_class(p$theme$plot.subtitle, "element_blank")
  expect_match(p$labels$caption, "No confidence intervals", fixed = TRUE)
})

# ---- doctor-clinical-views ----

clinical_binary_run <- function(thresholds = NULL) {
  bi_run(sb_example_binary(120, seed = 42), "outcome",
    c("gene_a", "gene_b"), "patient_id", models = "Logistic",
    folds = 3, plots = FALSE, decision_thresholds = thresholds
  )
}

test_that("short binary entrypoint retains declared decision thresholds", {
  thresholds <- c(.1, .3, .5)
  run <- clinical_binary_run(thresholds)

  expect_identical(run$development_task$decision_thresholds, thresholds)
  expect_identical(run$comparison$results[[1]]$metrics$decision_curve$status, "success")
})

test_that("clinical views accept a doctor run without fitting again", {
  run <- clinical_binary_run(c(.1, .3, .5))
  saved <- run$comparison$results[[1]]$metrics$decision_curve$curve
  original <- run$comparison$results[[1]]$predictions
  local_mocked_bindings(sb_fit_final = function(...) stop("Unexpected refit"))

  views <- sb_clinician_views(run)

  expect_setequal(names(views$plots), c("ROC", "Calibration", "Decision curve"))
  expect_equal(views$plots[["Decision curve"]]$data$x, rep(saved$threshold_probability, 3))
  expect_equal(sort(views$plots[["Decision curve"]]$data$y),
    sort(c(saved$model_net_benefit, saved$treat_all_net_benefit, saved$treat_none_net_benefit)))
  expect_identical(run$comparison$results[[1]]$predictions, original)
})

test_that("named external cohorts retain saved calibration and decision curves", {
  run <- clinical_binary_run(c(.1, .3, .5))
  external <- sb_example_binary(100, seed = 17)
  completed <- sb_complete_run(run, "Logistic", list(first = external),
    plots = FALSE, reverse = FALSE
  )
  evaluation <- completed$evaluations$first
  cohorts <- list(first = evaluation, second = evaluation)

  views <- sb_clinician_views(cohorts)

  expect_identical(views$status$status, rep("available", 3))
  expect_setequal(names(views$plots), c("ROC", "Calibration", "Decision curve"))
  curve <- sb_plot_calibration(evaluation)$data
  expect_equal(views$plots$Calibration$data$x, rep(curve$x, 2))
  expect_equal(views$plots$Calibration$data$y, rep(curve$y, 2))
  expect_match(paste(views$plots$Calibration$data$panel, collapse = " "), "first")
  expect_match(paste(views$plots$Calibration$data$panel, collapse = " "), "second")
})

test_that("completed doctor results display external rather than development evidence", {
  run <- clinical_binary_run(c(.1, .3, .5))
  completed <- sb_complete_run(run, "Logistic",
    list(validation = sb_example_binary(100, seed = 17)),
    plots = FALSE, reverse = FALSE
  )
  original <- completed$evaluations$validation$scored_predictions
  local_mocked_bindings(sb_fit_final = function(...) stop("Unexpected refit"))

  views <- sb_clinician_views(completed)

  expect_setequal(names(views$plots), c("ROC", "Calibration", "Decision curve"))
  expect_match(paste(views$plots$Calibration$data$panel, collapse = " "), "External evaluation")
  expect_identical(completed$evaluations$validation$scored_predictions, original)
})

test_that("short survival entrypoint retains the clinical horizon and thresholds", {
  thresholds <- c(.1, .3, .5)
  run <- surv_run(sb_example_survival(120, seed = 42), "time", "event",
    c("gene_a", "gene_b"), "patient_id", models = "Cox",
    folds = 3, plots = FALSE, prediction_horizon = 6,
    decision_thresholds = thresholds
  )

  expect_identical(run$development_task$prediction_horizon, 6)
  expect_identical(run$development_task$decision_thresholds, thresholds)
  expect_identical(run$comparison$results[[1]]$metrics$time_auc$status, "success")
  expect_identical(run$comparison$results[[1]]$metrics$decision_curve$status, "success")
})

test_that("missing external evaluations are not replaced with training evidence", {
  run <- clinical_binary_run()
  completed <- sb_complete_run(run, "Logistic", plots = FALSE)

  result <- tryCatch(sb_clinician_views(completed), error = identity)

  expect_s3_class(result, "stormbreaker_bad_visualization_input")
  expect_match(conditionMessage(result), "No saved external evaluations", fixed = TRUE)
})

test_that("doctor comparison plots show model names instead of parameter IDs", {
  run <- bi_run(sb_example_binary(120, seed = 42), "outcome",
    c("gene_a", "gene_b"), "patient_id", models = c("Logistic", "NaiveBayes"),
    folds = 3
  )
  plot <- run$plots$performance

  expect_identical(plot$labels$title, "Model comparison")
  expect_null(plot$labels$subtitle)
  expect_null(plot$labels$caption)
  expect_setequal(plot$scales$get_scales("x")$get_labels(),
    c("Logistic", "Gaussian NB"))
  expect_setequal(as.character(plot$data$pipeline_id),
    vapply(run$pipelines, function(pipeline) pipeline$id, character(1)))
})

# ---- doctor-readability ----

readability_binary_run <- function(plots = FALSE) {
  data <- sb_example_binary(80, seed = 11)
  bi_run(data, "outcome", c("gene_a", "gene_b"), "patient_id",
    models = "Logistic", folds = 3, plots = plots
  )
}

test_that("an optional plot failure does not discard a successful model run", {
  local_mocked_bindings(sb_plots = function(...) stop("Plot device unavailable"))

  run <- readability_binary_run(plots = TRUE)

  expect_identical(run$comparison$status, "success")
  expect_identical(run$details$Logistic$fit$status, "success")
  expect_contains(run$issues$reason, "Plot device unavailable")
  expect_contains(run$issues$stage, "plot algorithm")
})

test_that("completion reuses saved algorithm figures without calling the plotter", {
  run <- readability_binary_run(plots = TRUE)
  local_mocked_bindings(sb_plots = function(...) stop("Unexpected second plot call"))

  completed <- sb_complete_run(run, "Logistic")

  expect_identical(completed$fit, run$details$Logistic$fit)
  expect_identical(completed$plots$algorithm, run$details$Logistic$figures$plots)
  expect_identical(completed$issues$algorithm, run$details$Logistic$figures$unavailable)
})

test_that("printing a comparison shows a short summary instead of patient data", {
  run <- readability_binary_run()
  run$development_task$data$patient_id[1] <- "PRIVATE_PATIENT_731"

  printed <- capture.output(returned <- withVisible(print(run)))

  expect_s3_class(run, "stormbreaker_run")
  expect_identical(returned$value, run)
  expect_identical(returned$visible, FALSE)
  expect_lte(length(printed), 12L)
  expect_match(paste(printed, collapse = "\n"), "sb_complete_run", fixed = TRUE)
  expect_identical(any(grepl("PRIVATE_PATIENT_731", printed, fixed = TRUE)), FALSE)
})

test_that("printing a completed run leaves its detailed results untouched", {
  run <- readability_binary_run()
  external <- run$development_task$data
  external$patient_id <- paste0("PRIVATE_PATIENT_952_", seq_len(nrow(external)))
  completed <- sb_complete_run(run, "Logistic", list(validation = external), plots = FALSE)

  printed <- capture.output(returned <- withVisible(print(completed)))

  expect_s3_class(completed, "stormbreaker_completed_run")
  expect_identical(returned$value, completed)
  expect_identical(returned$visible, FALSE)
  expect_lte(length(printed), 12L)
  expect_match(paste(printed, collapse = "\n"), "Logistic", fixed = TRUE)
  expect_match(paste(printed, collapse = "\n"), "sb_export_run", fixed = TRUE)
  expect_identical(any(grepl("PRIVATE_PATIENT_952", printed, fixed = TRUE)), FALSE)
})

test_that("exported indexes retain their column names when no cohorts or plots exist", {
  completed <- sb_complete_run(readability_binary_run(), "Logistic", plots = FALSE)
  directory <- file.path(withr::local_tempdir(), "results")

  sb_export_run(completed, directory)

  cohorts <- read.csv(file.path(directory, "cohorts.csv"))
  figures <- read.csv(file.path(directory, "figures.csv"))
  notes <- read.csv(file.path(directory, "figure-notes.csv"))
  expect_named(cohorts, c("cohort", "file", "status"))
  expect_named(figures, c("plot", "file"))
  expect_named(notes, c("plot", "note"))
  expect_equal(nrow(cohorts), 0L)
  expect_equal(nrow(figures), 0L)
  expect_equal(nrow(notes), 0L)
})

test_that("skipping completion plots neither redraws nor discards a saved fit", {
  run <- readability_binary_run(plots = TRUE)
  local_mocked_bindings(
    sb_plots = function(...) stop("Unexpected plotting"),
    sb_fit_final = function(...) stop("Unexpected refit")
  )

  completed <- sb_complete_run(run, "Logistic", plots = FALSE)

  expect_length(completed$plots, 0L)
  expect_identical(completed$fit, run$details$Logistic$fit)
  expect_setequal(completed$genes$gene, c("gene_a", "gene_b"))
})

test_that("a failed final fit remains visible in the completed summary", {
  run <- readability_binary_run()
  local_mocked_bindings(sb_fit_final = function(...) {
    list(status = "fit_failed", error_message = "No usable features")
  })

  completed <- sb_complete_run(run, "Logistic", plots = FALSE)

  expect_s3_class(completed, "stormbreaker_completed_run")
  expect_identical(completed$issues$fit$error_message, "No usable features")
  expect_match(paste(capture.output(print(completed)), collapse = "\n"),
    "fit_failed",
    fixed = TRUE
  )
  expect_length(completed$predictions, 0L)
})

test_that("export preserves the original reasons for nested failures", {
  completed <- sb_complete_run(readability_binary_run(), "Logistic", plots = FALSE)
  completed$issues <- list(
    external = list(broken = "Missing outcome column"),
    algorithm = list(path = "Saved path unavailable"),
    fit = list(status = "fit_failed", error_message = "Original engine error")
  )
  directory <- file.path(withr::local_tempdir(), "results")

  sb_export_run(completed, directory)

  issues <- read.csv(file.path(directory, "issues.csv"))
  expect_identical(issues$scope, c("/external/broken", "/algorithm/path", "/fit"))
  expect_identical(issues$message, c(
    "Missing outcome column", "Saved path unavailable", "fit_failed; Original engine error"
  ))
  expect_identical(readRDS(file.path(directory, "issues.rds")), completed$issues)
})

test_that("export never replaces an existing directory or its contents", {
  directory <- withr::local_tempdir()
  marker <- file.path(directory, "keep.txt")
  writeLines("Original contents", marker)

  expect_error(sb_export_run(list(), directory), "existing paths are not overwritten", fixed = TRUE)

  expect_identical(readLines(marker), "Original contents")
  expect_identical(list.files(directory), "keep.txt")
})

test_that("survival completion preserves risk groups while reusing algorithm plots", {
  data <- sb_example_survival(80, seed = 21)
  run <- surv_run(data, "time", "event", c("gene_a", "gene_b"), "patient_id",
    models = "Cox", folds = 3, plots = TRUE
  )
  external <- data
  external$patient_id <- paste0("external-", external$patient_id)
  local_mocked_bindings(sb_plots = function(...) stop("Unexpected second plot call"))

  completed <- sb_complete_run(run, "Cox", list(validation = external))

  expect_s3_class(run, "stormbreaker_run")
  expect_s3_class(completed, "stormbreaker_completed_run")
  expect_identical(completed$plots$algorithm, run$details$Cox$figures$plots)
  expect_identical(completed$evaluations$validation$status, "success")
  expect_equal(completed$risk_groups$validation$n_low + completed$risk_groups$validation$n_high, 80L)
  expect_match(paste(capture.output(print(completed)), collapse = "\n"), "survival", fixed = TRUE)
})

# ---- figure-polish ----

test_that("shared styling and OOB colours preserve recorded values", {
  skip_if_not_installed("ggplot2")
  t <- sb_viz_theme()
  expect_identical(t$text$family, "sans")
  expect_identical(t$axis.text$colour, "#333333")
  d <- data.frame(trees = 1:3, `0` = c(.4, .3, .2), `1` = c(.5, .4, .3),
    oob = c(.45, .35, .25), check.names = FALSE)
  p <- sb_selector_oob_plot(d)
  expect_equal(p$data$error, unlist(d[-1], use.names = FALSE))
  expect_equal(ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")$map(c("0", "1", "oob")),
    c("#2878A0", "#B86236", "#686868"))
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})

test_that("mixture correlation labels fit without changing correlations", {
  skip_if_not_installed("ggplot2")
  d <- expand.grid(feature = paste0("gene_", 1:12), feature2 = paste0("gene_", 1:12),
    class = c("0", "1"), component = 1)
  d$correlation <- as.numeric(d$feature == d$feature2)
  p <- sb_plot_mclust_da_diagnostic(list(type = "mclust_da_correlations", table = d))
  expect_identical(p$data, d)
  expect_equal(p$theme$axis.text.x$angle, 45)
  expect_equal(p$coordinates$ratio, 1)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})

# ---- successful-display ----

test_that("default performance display omits failed labels while preserving source", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("binary",TRUE)$external
  x$evaluations[[1]]$status <- "failed"
  before <- x
  p <- sb_plot_performance(x,"roc_auc")
  expect_identical(p$scales$get_scales("x")$limits,x$evaluations[[2]]$pipeline_id)
  expect_identical(x,before)
  audit <- sb_plot_performance(x,"roc_auc",show_unavailable=TRUE)
  expect_identical(audit$scales$get_scales("x")$limits,rev(names(x$evaluations)))
})

test_that("default calibration and DCA omit failed facets but retain diagnostics", {
  skip_if_not_installed("ggplot2")
  for(type in c("binary","survival")) {
    x <- make_sb_curve_view_fixture(type,TRUE)$external
    x$evaluations[[1]]$status <- "failed"
    before <- x
    for(fun in list(sb_plot_calibration,sb_plot_dca)) {
      p <- fun(x)
      expect_identical(nrow(ggplot2::ggplot_build(p)$layout$layout),1L)
      expect_identical(unique(as.character(p$data$pipeline_id)),x$evaluations[[2]]$pipeline_id)
      expect_identical(nrow(attr(p,"stormbreaker_curve_status")),2L)
      audit <- fun(x,show_unavailable=TRUE)
      expect_identical(nrow(ggplot2::ggplot_build(audit)$layout$layout),2L)
    }
    expect_identical(x,before)
  }
})

test_that("all unavailable displays clearly report absence without pipeline axes", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("binary")$external
  x$status <- "failed"
  for(p in list(sb_plot_performance(x,"roc_auc"),sb_plot_calibration(x),sb_plot_dca(x))) {
    expect_match(ggplot2::ggplot_build(p)$data[[1]]$label,"No usable results",fixed=TRUE)
    expect_identical(nrow(ggplot2::ggplot_build(p)$layout$layout),1L)
  }
})

# ---- visualization-integration ----

test_that("all five real Meta outputs and the synthetic tutorial retain saved endpoints", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("metafor")
  path <- system.file("examples", "advanced-examples.R", package = "stormbreaker")
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
  x <- env$stormbreaker_visualization_example(reps = 40L)
  expect_length(x$plots, 19L)
  for (p in x$plots) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  for (m in x$meta) {
    before <- serialize(m, NULL)
    p <- sb_plot_meta_forest(m)
    expect_identical(tail(p$data$point_estimate, 2), rep(m$summary$estimate, 2))
    expect_identical(tail(p$data$lower, 2), c(m$summary$lower, m$summary$prediction_lower))
    expect_identical(tail(p$data$upper, 2), c(m$summary$upper, m$summary$prediction_upper))
    expect_identical(serialize(m, NULL), before)
  }
  expected <- subset(x$survival_bootstrap$intervals, metric == "ipcw_brier_score")
  p <- sb_plot_external_intervals(x$survival_bootstrap, "ipcw_brier_score")
  expect_identical(p$data$lower, expected$lower)
  expect_identical(p$data$upper, expected$upper)
  expect_match(as.character(p$data$target_label), "6 month", fixed = TRUE)
})

test_that("real repeated and IECV results remain separate plotting panels", {
  skip_if_not_installed("ggplot2")
  tasks <- lapply(1:2, function(i) {
    d <- sb_example_binary(60, seed = 500 + i)
    d$patient_id <- paste0(i, "_", d$patient_id)
    sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  })
  names(tasks) <- c("A", "B")
  pipelines <- list(sb_pipeline(learner = "logistic"))
  repeated <- sb_compare_repeated_cv(tasks[[1]], pipelines, folds = 3, seeds = c(501, 502))
  iecv <- sb_run_internal_external_cv(tasks, pipelines)
  for (x in list(repeated, iecv)) {
    before <- serialize(x, NULL)
    roc <- sb_plot_roc(x, by_fold = TRUE)
    pr <- sb_plot_pr(x)
    expect_equal(length(unique(roc$data$panel)), if (identical(x, repeated)) 6L else 2L)
    expect_equal(length(unique(pr$data$panel)), length(unique(roc$data$panel)))
    expect_s3_class(ggplot2::ggplot_build(sb_plot_performance_heatmap(x, "roc_auc")), "ggplot_built")
    expect_identical(serialize(x, NULL), before)
  }
})

test_that("real glmnet active sets and stack roles are usable without pooling", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("glmnet")
  d <- sb_example_binary(90, seed = 510)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  cv <- sb_run_cv(task, sb_pipeline(learner = "glmnet_logistic", glmnet_alpha = 1,
    inner_folds = 3, inner_seed = 511), folds = 3, seed = 511)
  p <- sb_plot_selection_frequency(cv, source = "glmnet_active_set")
  expect_identical(p$data$selection_frequency,
    cv$glmnet_active_set_stability$feature_frequency$selection_frequency)
  for (fun in list(sb_plot_selection_jaccard, sb_plot_selection_upset)) {
    expect_s3_class(ggplot2::ggplot_build(fun(cv, source = "glmnet_active_set")), "ggplot_built")
  }
  stack <- sb_stack(list(sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")), inner_folds = 3)
  result <- sb_run_stack_cv(task, stack, folds = 3, seed = 512)
  p <- sb_plot_roc(result, by_fold = TRUE)
  expect_setequal(unique(p$data$model_role), c("member", "stack"))
  expect_equal(length(unique(p$data$panel)), 6L)
  expect_s3_class(ggplot2::ggplot_build(sb_plot_performance_heatmap(result, "roc_auc")), "ggplot_built")
})

test_that("real Cox and naive Bayes SHAP contributions reconstruct plotted scores", {
  skip_if_not_installed("ggplot2")
  for (learner in c("coxph", "gaussian_naive_bayes")) {
    d <- if (learner == "coxph") sb_example_survival(80, seed = 520) else sb_example_binary(80, seed = 521)
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (learner == "coxph") sb_surv_task(d, "patient_id", "time", "event", features, time_unit = "month") else
      sb_binary_task(d, "patient_id", "outcome", features)
    fit <- sb_fit_final(task, sb_pipeline(learner = learner))
    explanation <- sb_additive_shap(fit, d[1:10, ], "patient_id", d[features])
    before <- serialize(explanation, NULL)
    for (fun in list(sb_plot_shap_importance, sb_plot_shap_beeswarm, sb_plot_shap_waterfall)) {
      expect_s3_class(ggplot2::ggplot_build(fun(explanation, top_n = 2)), "ggplot_built")
    }
    expect_identical(serialize(explanation, NULL), before)
  }
})
})

# ---- composite-result-views ----
local({
test_that("repeated tables keep repetitions without pooling or changing inputs", {
  d <- sb_example_binary(60, seed = 1)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  x <- sb_compare_repeated_cv(task, list(sb_pipeline(learner = "logistic")),
                              folds = 3, seeds = c(1, 1))
  before <- serialize(x, NULL)
  tab <- sb_result_table(x, metric = "roc_auc")
  expect_equal(tab$repeat_id, 1:2)
  expect_equal(tab$estimate, x$comparison$roc_auc)
  expect_equal(nrow(sb_result_table(x, "predictions")), 120)
  expect_equal(nrow(sb_result_table(x, "runs")), 6)
  expect_identical(serialize(x, NULL), before)
})

test_that("native stack and repeated stack views retain member and stack values", {
  skip_if_not_installed("glmnet")
  for (type in c("binary", "survival")) {
    fixture <- make_sb_curve_view_fixture(type)
    pipelines <- if (type == "binary") list(sb_pipeline(learner = "logistic"),
      sb_pipeline(learner = "gaussian_naive_bayes")) else list(sb_pipeline(),
      sb_pipeline(selector = "univ_cox_top_k", selector_n_features = 1))
    stack <- sb_stack(pipelines, inner_folds = 3, inner_seed = 31)
    x <- sb_run_stack_cv(fixture$task, stack, folds = 3, seed = 37)
    before <- serialize(x, NULL)
    metrics <- if (type == "binary") c("roc_auc", "average_precision", "brier_score", "log_loss") else
      c("harrell_c", "time_auc", "uno_c", "ipcw_brier")
    for (metric in metrics) {
      tab <- sb_result_table(x, metric = metric)
      expect_equal(tab$model_role, c("member", "member", "stack"))
      expect_equal(tail(tab$estimate, 1), x$metrics[[metric]]$estimate)
      expect_equal(head(tab$estimate, 2), sb_result_table(x$members, metric = metric)$estimate)
      if (type == "survival" && metric != "harrell_c") {
        expect_equal(tab$time_unit, rep("month", 3))
        expect_equal(tab$prediction_horizon, rep(6, 3))
      }
    }
    predictions <- sb_result_table(x, "predictions")
    expect_equal(predictions[predictions$model_role == "stack", names(x$predictions)],
                 x$predictions, ignore_attr = TRUE)
    expect_equal(nrow(sb_result_table(x, "runs")), 9)
    repeated <- structure(list(repetitions = list(x, x)),
                          class = "stormbreaker_repeated_cv_comparison")
    tab <- sb_result_table(repeated, metric = metrics[1])
    expect_equal(tab$repeat_id, rep(1:2, each = 3))
    expect_equal(length(unique(tab$view_id)), 6)
    expect_identical(serialize(x, NULL), before)
  }
})

test_that("stack metadata uses saved prediction units only when unambiguous", {
  x <- make_sb_curve_view_fixture("survival")$cv
  stack <- structure(list(status = x$status, task_type = "survival",
    predictions = x$predictions, metrics = x$metrics, runs = x$runs,
    task = list(), stack = list(id = "saved-stack"),
    members = structure(list(results = list(x)), class = "stormbreaker_cv_comparison")),
    class = "stormbreaker_stack_cv_result")
  stack$predictions$time_unit[1] <- "day"
  tab <- sb_result_table(stack, metric = "time_auc")
  expect_true(is.na(tail(tab$time_unit, 1)))
  stack$predictions$time_unit <- NULL
  tab <- sb_result_table(stack, metric = "time_auc")
  expect_true(is.na(tail(tab$time_unit, 1)))
})

test_that("composite readers cannot fit or recalculate and preserve missing values", {
  x <- make_sb_curve_view_fixture("survival")$cv
  y <- x
  y$metrics$ipcw_brier$estimate <- 1.25
  y$metrics$ipcw_brier$prediction_horizon <- 12
  y$metrics$ipcw_brier$time_unit <- "day"
  repeated <- structure(list(repetitions = list(x, y)),
                        class = "stormbreaker_repeated_cv_comparison")
  before <- serialize(repeated, NULL)
  fail <- function(...) stop("Must not execute fitting or metrics")
  local_mocked_bindings(sb_run_cv = fail, sb_fit_final = fail,
                       sb_compare_external = fail, sb_evaluate_external = fail)
  tab <- sb_result_table(repeated, metric = "ipcw_brier")
  expect_equal(tab$prediction_horizon, c(6, 12))
  expect_equal(tab$estimate[2], 1.25)
  expect_equal(tab$time_unit[2], "day")
  expect_identical(serialize(repeated, NULL), before)
})

test_that("native IECV stacks preserve their actual saved evaluation roles", {
  skip_if_not_installed("glmnet")
  cohorts <- lapply(1:2, function(i) {
    d <- sb_example_binary(90, seed = i)
    d$patient_id <- paste0(i, "_", d$patient_id)
    sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  })
  names(cohorts) <- c("A", "B")
  stack <- sb_stack(list(sb_pipeline(learner = "logistic"),
                        sb_pipeline(learner = "gaussian_naive_bayes")), inner_folds = 3)
  x <- sb_run_internal_external_cv(cohorts, stack)
  before <- serialize(x, NULL)
  tab <- sb_result_table(x, metric = "roc_auc")
  expect_equal(tab$model_role, rep("stack", 2))
  expect_equal(tab$heldout_cohort_id, c("A", "B"))
  expect_equal(tab$estimate, x$results$roc_auc)
  expect_identical(serialize(x, NULL), before)
})

test_that("empty prediction leaves and raw evaluation failures keep their identity", {
  x <- make_sb_curve_view_fixture()$external
  failed <- x
  failed$status <- "failed"
  failed$scored_predictions <- x$scored_predictions[0, , drop = FALSE]
  failed$metrics <- list()
  failed$error_message <- "saved ID mismatch"
  iecv <- structure(list(folds = list(
    list(heldout_cohort_id = "A", comparison = x),
    list(heldout_cohort_id = "B", comparison = failed))),
    class = "stormbreaker_internal_external_cv")
  expect_equal(nrow(sb_result_table(iecv, "predictions")), nrow(x$scored_predictions))
  tab <- sb_result_table(iecv, metric = "roc_auc")
  expect_true(is.na(tab$estimate[2]))
  expect_equal(tab$status[2], "not_saved")
  failures <- sb_failure_table(iecv)
  expect_equal(failures$heldout_cohort_id, "B")
  expect_equal(failures$error_message, "saved ID mismatch")
})

test_that("IECV tables identify held-out cohorts and development validation", {
  cohorts <- lapply(1:2, function(i) {
    d <- sb_example_binary(60, seed = i)
    d$patient_id <- paste0(i, "_", d$patient_id)
    sb_binary_task(d, "patient_id", "outcome", "gene_a")
  })
  names(cohorts) <- c("site <A>", "site B")
  x <- sb_run_internal_external_cv(cohorts, list(sb_pipeline(learner = "logistic")))
  tab <- sb_result_table(x, metric = "roc_auc")
  expect_equal(tab$heldout_cohort_id, names(cohorts))
  expect_equal(tab$estimate, x$results$roc_auc)
  expect_true(all(tab$validation_scheme == "internal_external_cv"))
  expect_equal(nrow(sb_result_table(x, "predictions")), 120)
  expect_equal(nrow(sb_result_table(x, "runs")), 2)
})
})

# ---- results ----
local({
# ---- failure-table ----

test_that("failure tables distinguish run and metric failures without mutation", {
  x <- make_sb_curve_view_fixture("binary")$external
  x$status <- "failed"
  x$error_stage <- "patient_match"
  x$error_class <- "saved_match_error"
  x$error_message <- "Saved original reason"
  x$metrics <- list()
  before <- x
  tab <- sb_failure_table(x,cohort_label="External A")
  expect_identical(tab$pipeline_id,x$pipeline_id)
  expect_identical(tab$issue_type,"pipeline")
  expect_identical(tab$stage,"external_evaluation")
  expect_identical(tab$source_error_stage,"patient_match")
  expect_identical(tab$error_message,"Saved original reason")
  expect_identical(tab$cohort_label,"External A")
  expect_identical(x,before)

  y <- make_sb_curve_view_fixture("survival")$external
  y$status <- "success"
  y$metrics <- list(harrell_c=y$metrics$harrell_c)
  y$metrics$harrell_c$status <- "not_evaluable"
  y$metrics$harrell_c$error_message <- "No comparable pairs saved"
  tab <- sb_failure_table(y)
  expect_identical(tab$issue_type,"metric")
  expect_identical(tab$metric,"harrell_c")
  expect_identical(tab$run_status,"success")
  expect_identical(tab$status,"not_evaluable")
  expect_identical(tab$prediction_horizon,6)
  expect_identical(tab$time_unit,"month")
})

test_that("fold records retain saved details and successful metrics are not failures", {
  x <- make_sb_curve_view_fixture("binary")$cv
  x$metrics <- list(roc_auc=x$metrics$roc_auc)
  expect_identical(nrow(sb_failure_table(x)),0L)
  x$status <- "failed"
  x$runs$status[1] <- "failed"
  x$runs$error_class[1] <- "saved_fit_error"
  x$runs$error_message[1] <- "Original fold reason"
  tab <- sb_failure_table(x)
  folds <- tab[tab$issue_type=="run",]
  expect_identical(folds$fold,as.integer(x$runs$fold[1]))
  expect_identical(folds$error_message,"Original fold reason")
  expect_identical(folds$stage,"cv_fold")
  expect_identical(folds$source_error_stage,NA_character_)
  expect_identical("patient_id" %in% names(tab),FALSE)
})

test_that("named cohorts remain distinct and unrequested metrics are not errors", {
  x <- make_sb_curve_view_fixture("binary")$external
  x$metrics <- list(log_loss=x$metrics$log_loss,unused=list(status="not_requested"))
  x$metrics$log_loss$estimate <- Inf
  tab <- sb_failure_table(list(A=x,B=x))
  expect_identical(tab$cohort_label,c("A","B"))
  expect_identical(tab$issue_type,rep("metric_nonfinite",2L))
  expect_identical(tab$status,rep("success",2L))
  expect_identical(tab$estimate,rep(Inf,2L))
})

test_that("failure CSV preserves errors, empty headers, and existing output", {
  x <- make_sb_curve_view_fixture("binary")$external
  x$metrics <- list(roc_auc=x$metrics$roc_auc)
  file <- tempfile(fileext=".csv")
  path <- sb_export_failures(x,file)
  expect_identical(path,normalizePath(file,winslash="/"))
  expect_identical(nrow(read.csv(file)),0L)
  saved <- readLines(file)
  expect_error(sb_export_failures(x,file),class="stormbreaker_output_exists")
  expect_identical(readLines(file),saved)
  x$status <- "failed"; x$error_message <- 'Reason, with "quotes" and newline\nsecond line'
  file <- tempfile(fileext=".csv")
  sb_export_failures(x,file)
  tab <- read.csv(file,stringsAsFactors=FALSE)
  expect_identical(tab$error_message[tab$issue_type=="pipeline"],x$error_message)
  expect_error(sb_failure_table(list(x)),class="stormbreaker_bad_failure_view")
})

# ---- result-table ----

test_that("scalar tables copy all five supported metrics without changing input", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  x <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1)
  before <- x
  for (metric in c("roc_auc", "average_precision", "brier_score", "log_loss")) {
    tab <- sb_result_table(x, metric = metric)
    saved <- x$metrics[[metric]]
    expect_identical(nrow(tab), 1L)
    expect_identical(tab$pipeline_id, x$pipeline$id)
    expect_identical(tab$task_type, "binary")
    expect_identical(tab$metric, metric)
    expect_identical(tab$estimate, saved$estimate)
    expect_identical(tab$status, saved$status)
    expect_identical(tab$run_status, x$status)
    expect_identical(tab$scope, saved$scope)
    expect_identical(tab$aggregation, saved$aggregation)
  }
  expect_identical(x, before)

  s <- make_sb_fixture()
  st <- sb_surv_task(s, "patient_id", "time", "event", "gene_a", "month")
  sx <- sb_run_cv(st, folds = 3, seed = 1)
  tab <- sb_result_table(sx, metric = "harrell_c")
  expect_identical(tab$task_type, "survival")
  expect_identical(tab$estimate, sx$metrics$harrell_c$estimate)
  expect_identical(tab$scope, sx$metrics$harrell_c$scope)
  expect_identical(tab$aggregation, sx$metrics$harrell_c$aggregation)
})

test_that("prediction and run tables preserve saved columns and row order", {
  for (binary in c(TRUE, FALSE)) {
    if (binary) {
      d <- make_sb_binary_fixture()
      task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
      pipeline <- sb_pipeline(learner = "logistic")
    } else {
      d <- make_sb_fixture()
      task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month")
      pipeline <- sb_pipeline()
    }
    x <- sb_run_cv(task, pipeline, folds = 3, seed = 1)
    before <- x
    for (component in c("predictions", "runs")) {
      tab <- sb_result_table(x, component)
      expect_identical(tab[names(x[[component]])], x[[component]])
      expect_identical(tab$task_type, rep(if (binary) "binary" else "survival",
                                         nrow(x[[component]])))
    }
    expect_identical(x, before)
    expect_identical("error_stage" %in% names(sb_result_table(x, "runs")), FALSE)
  }
})

test_that("comparison tables retain every supplied pipeline and fold", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  x <- sb_compare_cv(task, list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), folds = 3, seed = 1)
  before <- x
  tab <- sb_result_table(x, metric = "roc_auc")
  expect_identical(tab$pipeline_id, names(x$results))
  expect_identical(tab$estimate, unname(vapply(
    x$results, function(z) z$metrics$roc_auc$estimate, numeric(1)
  )))
  for (component in c("predictions", "runs")) {
    tab <- sb_result_table(x, component)
    expected <- do.call(rbind, lapply(unname(x$results), `[[`, component))
    rownames(expected) <- NULL
    expect_identical(tab[names(expected)], expected)
  }
  expect_identical(x, before)
})

test_that("incomplete and nonfinite estimates are retained rather than repaired", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  a <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1)
  b <- a
  b$pipeline$id <- "synthetic_failure"
  b$status <- "failed"
  b$metrics$log_loss$status <- "incomplete"
  b$metrics$log_loss$estimate <- NA_real_
  c <- a
  c$pipeline$id <- "synthetic_infinite_loss"
  c$metrics$log_loss$estimate <- Inf
  x <- structure(list(results = list(b, a, c)),
                 class = "stormbreaker_cv_comparison")
  tab <- sb_result_table(x, metric = "log_loss")
  expect_identical(tab$pipeline_id,
                   c("synthetic_failure", a$pipeline$id, "synthetic_infinite_loss"))
  expect_identical(tab$estimate, c(NA_real_, a$metrics$log_loss$estimate, Inf))
  expect_identical(tab$run_status, c("failed", "success", "success"))
  expect_identical(tab$status, c("incomplete", a$metrics$log_loss$status,
                               a$metrics$log_loss$status))
})

test_that("extraction never calls fitting or metric computation", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  x <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1)
  local_mocked_bindings(
    sb_run_cv = function(...) stop("must not run CV"),
    sb_fit_final = function(...) stop("must not fit"),
    sb_binary_auc_metric = function(...) stop("must not recompute")
  )
  expect_identical(sb_result_table(x, metric = "roc_auc")$estimate,
                   x$metrics$roc_auc$estimate)
  expect_identical(nrow(sb_result_table(x, "predictions")), nrow(d))
})

test_that("unsupported requests have explicit errors", {
  d <- make_sb_fixture()
  task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month")
  x <- sb_run_cv(task, folds = 3, seed = 1)
  for (metric in list(NULL, NA_character_, character(), c("harrell_c", "roc_auc"),
                      "calibration", "unsupported_score", 1)) {
    expect_error(sb_result_table(x, metric = metric),
                 class = "stormbreaker_bad_result_metric")
  }
  expect_error(sb_result_table(x, metric = "roc_auc"),
               class = "stormbreaker_result_metric_not_saved")
  expect_error(sb_result_table(x, "predictions", metric = "harrell_c"),
               class = "stormbreaker_unused_result_metric")
  expect_error(sb_result_table(list(), metric = "harrell_c"),
               class = "stormbreaker_bad_result_view")
  expect_error(sb_result_table(structure(list(results = list()),
                                        class = "stormbreaker_cv_comparison"),
                               metric = "harrell_c"),
               class = "stormbreaker_bad_result_view")
})
})

# ---- select-model ----
local({
selection_example <- function() {
  make_task <- function(prefix, seed) {
    d <- sb_example_binary(60, seed)
    d$patient_id <- paste0(prefix, d$patient_id)
    sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  }
  sb_compare_cohorts(make_task("train", 1),
    list(Validation = make_task("val", 2), Test = make_task("test", 3)),
    roles = c(Validation = "validation", Test = "test"),
    models = c("logistic", "gaussian_naive_bayes"), folds = 3, view = FALSE
  )
}

test_that("validation selection returns a directly predictable existing fit", {
  compared <- selection_example()
  table <- compared$cohorts$Validation$comparison
  selected <- sb_select_model(compared)
  id <- table$pipeline_id[which.max(table$roc_auc)]
  index <- match(id, vapply(compared$fits, function(x) x$pipeline$id, character(1)))
  expect_s3_class(selected, "stormbreaker_final_fit")
  expect_identical(selected$pipeline$id, id)
  expect_identical(selected$learner_state, compared$fits[[index]]$learner_state)
  d <- sb_example_binary(20, 4)
  expect_identical(
    sb_predict_final(selected, d, "patient_id")$predictions,
    sb_predict_final(unserialize(serialize(selected, NULL)), d, "patient_id")$predictions
  )
  expect_identical(selected$model_selection$cohorts, "Validation")
})

test_that("validation scores are aligned by model ID and equally weighted", {
  compared <- selection_example()
  compared$cohorts$Second <- compared$cohorts$Validation
  compared$cohort_roles <- c(compared$cohort_roles, Second = "validation")
  compared$cohorts$Validation$comparison$roc_auc <- c(0.9, 0.6)
  compared$cohorts$Second$comparison$roc_auc <- c(0.2, 0.8)
  compared$cohorts$Second$comparison <- compared$cohorts$Second$comparison[2:1, ]
  selected <- sb_select_model(compared)
  expect_identical(selected$pipeline$id, compared$model_pipelines[[2]]$id)
  expect_equal(unname(selected$model_selection$ranking$validation_mean), c(0.7, 0.55))
  expect_identical(selected$model_selection$cohorts, c("Validation", "Second"))
})

test_that("test scores never change selection and ties keep declaration order", {
  compared <- selection_example()
  before <- sb_select_model(compared)$pipeline$id
  compared$cohorts$Test$comparison$roc_auc <- c(999, -999)
  expect_identical(sb_select_model(compared)$pipeline$id, before)
  compared$cohorts$Validation$comparison$roc_auc <- c(0.7, 0.7)
  expect_identical(sb_select_model(compared)$pipeline$id, compared$model_pipelines[[1]]$id)
  compared$fits[[1]]$status <- "failed"
  expect_identical(sb_select_model(compared)$pipeline$id, compared$model_pipelines[[2]]$id)
})

test_that("selection requires validation and complete eligible candidates", {
  compared <- selection_example()
  compared$cohort_roles["Validation"] <- "external"
  expect_error(sb_select_model(compared), "validation", class = "stormbreaker_no_validation_cohort")
  compared$cohort_roles["Validation"] <- "validation"
  compared$cohorts$Validation$comparison$roc_auc <- NA_real_
  expect_error(sb_select_model(compared), "eligible", class = "stormbreaker_no_eligible_model")
})

test_that("survival selection uses Harrell C and preserves the fitted state", {
  make_task <- function(prefix, seed) {
    d <- sb_example_survival(55, seed)
    d$patient_id <- paste0(prefix, d$patient_id)
    sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month")
  }
  compared <- sb_compare_cohorts(make_task("train", 10),
    list(Validation = make_task("val", 11)), roles = c(Validation = "validation"),
    pipelines = list(sb_pipeline(), sb_pipeline(selector = "univ_cox_top_k", selector_n_features = 1)),
    folds = 3, view = FALSE, progress = FALSE)
  selected <- sb_select_model(compared)
  tab <- compared$cohorts$Validation$comparison
  expected <- tab$pipeline_id[which.max(tab$harrell_c)]
  expect_identical(selected$pipeline$id, expected)
  expect_identical(selected$model_selection$metric, "harrell_c")
})
})

# ---- failure-continuation ----
local({
test_that("valid comparisons continue after a model-dependent failure below 100 inputs", {
  withr::local_seed(18)
  for(type in c("binary","survival")) {
    d <- if(type=="binary") sb_example_binary(45,seed=11) else sb_example_survival(45,seed=11)
    predictors <- paste0("input_",seq_len(35))
    for(name in predictors) d[[name]] <- rnorm(nrow(d))
    task <- if(type=="binary") sb_binary_task(d,"patient_id","outcome",predictors) else
      sb_surv_task(d,"patient_id","time","event",predictors,"month")
    pipelines <- if(type=="binary") list(sb_pipeline(learner="logistic"),sb_pipeline(learner="gaussian_naive_bayes")) else
      list(sb_pipeline(),sb_pipeline(selector="univ_cox_top_k",selector_n_features=1))
    x <- sb_compare_cv(task,pipelines,folds=3,seed=4)
    expect_identical(length(predictors),35L)
    expect_identical(unname(vapply(x$results,`[[`,character(1),"status")),c("failed","success"))
    expect_identical(names(x$results),vapply(pipelines,function(p) p$id,character(1)))
    issues <- sb_failure_table(x)
    failed_runs <- issues[issues$issue_type=="run",]
    expect_identical(failed_runs$fold,1:3)
    expect_identical(failed_runs$pipeline_id,rep(pipelines[[1]]$id,3L))
    expect_identical(all(!is.na(failed_runs$error_message)),TRUE)
  }
})
})

test_that("internal holdout labels survive the workflow and PNG export is opt-in", {
  skip_if_not_installed("ggplot2")
  d <- bi_example(60L, 61L)
  s <- sb_split_binary(bi_task(d, "patient_id", "outcome", c("gene_a", "gene_b")))
  run <- bi_run(s$development, models = "Logistic", folds = 3L, plots = FALSE)
  f <- finish_run(run, "Logistic", external = list(Test = s$assessment), plots = FALSE)
  expect_identical(f$evaluations$Test$assessment_role, "internal_holdout")
  expect_identical(unique(plot_roc(f$evaluations$Test)$data$validation_scheme), "internal_holdout")
  expect_identical(unique(sb_threshold_metrics(f$evaluations$Test, .5)$table$validation_scheme),
    "internal_holdout")
  graph <- ggplot2::ggplot(data.frame(x = 1, y = 1, label = "Text"), ggplot2::aes(x, y)) +
    ggplot2::geom_text(ggplot2::aes(label = label), family = "mono")
  typed <- sb_export_plot_font(graph, "serif")
  expect_identical(ggplot2::ggplot_build(typed)$data[[1L]]$family, "serif")
  expect_identical(ggplot2::ggplot_build(graph)$data[[1L]]$family, "mono")
  attr(graph, "stormbreaker_suggested_size") <- c(width = 2, height = 2)
  f$plots <- list(graph = graph)
  destination <- file.path(withr::local_tempdir(), "png")
  export_run(f, destination, figure_format = "png", dpi = 150L, font_family = "serif")
  index <- utils::read.csv(file.path(destination, "figures.csv"))
  expect_identical(index$file, "figure-001.png")
  png_header <- readBin(file.path(destination, index$file), "raw", n = 24L)
  expect_identical(readBin(png_header[17:24], "integer", n = 2L, endian = "big"), c(300L, 300L))
  expect_identical(utils::read.csv(file.path(destination, "cohorts.csv"))$assessment_role,
    "internal_holdout")
  default <- file.path(dirname(destination), "default-pdf")
  export_run(f, default)
  expect_identical(utils::read.csv(file.path(default, "figures.csv"))$file, "figure-001.pdf")
  tables <- file.path(dirname(destination), "tables")
  export_run(f, tables, figure_keys = character(0))
  expect_equal(nrow(utils::read.csv(file.path(tables, "figures.csv"))), 0L)
  expect_true(file.exists(file.path(tables, "cohort-001-metrics.csv")))
  expect_length(list.files(tables, pattern = "\\.(pdf|png)$"), 0L)
})

test_that("the binary plan separates screened predictors, selector-only methods, and restricted routes", {
  skip_if_not_installed("randomForest")
  skip_if_not_installed("glmnet")
  plan <- sb_binary_workflow_plan(models = c("Logistic", "Lasso", "Ridge", "ElasticNet",
    "Boruta", "pre_binary_rule_ensemble", "cancerclass_binary"))
  expect_identical(names(plan$pipelines), c("Logistic", "Lasso", "Ridge", "ElasticNet"))
  expect_identical(plan$candidates$status, c(rep("ready", 4), "selector_only",
    "selector_not_admitted", "selector_not_admitted"))
  expect_identical(plan$candidates$stage_selector[5L], "boruta_binary_roughfix")
  expect_identical(unname(vapply(plan$pipelines[-1L], function(p) p$learner_params$alpha, numeric(1))),
    c(1, 0, .5))
  expect_equal(plan$pipelines$Logistic$selector_params$top_n, 5L)
  expect_equal(plan$pipelines$Logistic$selector_params$spearman_cutoff, .7)
})

test_that("unfiltered binary routes require explicit opt-in and never replace missing engines", {
  testthat::local_mocked_bindings(
    sb_r_package_installed = function(package) package != "randomForest",
    .package = "stormbreaker"
  )
  plan <- sb_binary_workflow_plan(models = c("Logistic", "cancerclass_binary"),
    include_unfiltered = TRUE)
  expect_identical(plan$candidates$status, c("missing_dependency", "ready"))
  expect_identical(plan$candidates$effective_selector, c("rf_binary_importance", "none"))
  expect_identical(plan$candidates$missing_r_packages[1L], "randomForest")
  expect_match(plan$candidates$reason[2L], "distinct route")
  expect_identical(unname(vapply(plan$pipelines, `[[`, character(1), "selector")), "none")
})

test_that("a batch can use an alternate registered selector and only one route per model", {
  plan <- sb_binary_workflow_plan(models = c("Logistic", "NaiveBayes"),
    selector = "fixed_features", selector_features = c("gene_a", "gene_b"),
    feature_set_id = "declared-clinical")
  expect_identical(plan$candidates$status, c("ready", "ready"))
  expect_identical(names(plan$pipelines), c("Logistic", "NaiveBayes"))
  expect_identical(plan$pipelines$NaiveBayes$selector, "fixed_features")
  expect_error(sb_binary_workflow_plan(models = c("Logistic", "logistic")), "only once")
})

test_that("the binary batch keeps shared OOF data and retains deferred diagnostics for one final fit", {
  skip_if_not_installed("randomForest")
  withr::local_dir(withr::local_tempdir())
  d <- bi_example(100L, 61L)
  task <- bi_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  split <- sb_split_binary(task)
  plan <- sb_binary_workflow_plan(models = c("Logistic", "NaiveBayes"))
  run <- bi_run(split$development, pipelines = plan$pipelines, folds = 3L, plots = FALSE)
  expect_length(run$details, 0L)
  expect_length(run$plots, 0L)
  results <- run$comparison$results
  expect_identical(unname(vapply(results, `[[`, character(1), "status")), rep("success", 2L))
  expect_identical(results[[1L]]$provenance$fold_assignment, results[[2L]]$provenance$fold_assignment)
  expect_setequal(results[[1L]]$predictions$patient_id, split$development$data$patient_id)
  finished <- finish_run(run, "Logistic", external = list(Holdout = split$assessment),
    plots = FALSE, capture_diagnostics = TRUE)
  expect_identical(finished$fit$status, "success")
  expect_true(length(finished$fit$diagnostic_artifacts$selector) > 0L)
  expect_identical(finished$evaluations$Holdout$assessment_role, "internal_holdout")
  expect_setequal(finished$predictions$Holdout$patient_id, split$assessment$patient_id)
  expect_length(finished$plots, 0L)
  expect_length(list.files(), 0L)
})

test_that("figure export selects whole plots and can exclude all model data", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(query = seq_len(30L), y = seq_len(30L))
  graph <- ggplot2::ggplot(d, ggplot2::aes(query, y)) + ggplot2::geom_point()
  attr(graph, "stormbreaker_suggested_size") <- c(width = 2, height = 2)
  completed <- list(plots = list(Assessment = list(ROC = graph, Calibration = graph)))
  folder <- withr::local_tempdir()
  wanted <- "/Assessment/ROC"
  destination <- file.path(folder, "figures")
  expect_identical(sb_run_figure_keys(completed), c(wanted, "/Assessment/Calibration"))
  export_run(completed, destination, figure_keys = wanted, include_data = FALSE, figure_format = "png")
  expect_setequal(list.files(destination), c("figures.csv", "figure-notes.csv", "figure-001.png", "figure-002.png"))
  index <- utils::read.csv(file.path(destination, "figures.csv"))
  expect_true(all(startsWith(index$plot, wanted)))
  unknown <- file.path(folder, "unknown")
  expect_error(export_run(completed, unknown, figure_keys = "/Assessment/Typo"), "Unknown figure key")
  expect_false(file.exists(unknown))
})
