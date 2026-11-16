# ---- report-clinical ----
local({
test_that("visual gallery is offline, escaped and refuses existing outputs", {
  skip_if_not_installed("ggplot2")
  skip_if_not(capabilities("cairo"))
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 3:1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  file <- tempfile(fileext = ".html")
  device <- grDevices::dev.cur()
  result <- sb_render_visual_report(list("A < B" = p, Other = p), file,
    title = "<script>not executable</script>", width = c(4, 5), height = 4)
  html <- paste(readLines(file), collapse = "\n")
  expect_identical(result, normalizePath(file, winslash = "/"))
  expect_identical(grDevices::dev.cur(), device)
  expect_match(html, "&lt;script&gt;not executable&lt;/script&gt;", fixed = TRUE)
  expect_match(html, "A &lt; B", fixed = TRUE)
  expect_equal(length(gregexpr("data:image/svg+xml", html, fixed = TRUE)[[1]]), 2L)
  expect_identical(grepl("<script|https?://|<table", html), FALSE)
  expect_error(sb_render_visual_report(list(A = p), file), class = "stormbreaker_report_exists")
  expect_identical(paste(readLines(file), collapse = "\n"), html)
  expect_error(sb_render_visual_report(list(A = p), tempfile(fileext = ".html"), width = NA_real_),
    class = "stormbreaker_bad_report_input")
  expect_error(sb_render_visual_report(list(A = 1), tempfile(fileext = ".html")),
    class = "stormbreaker_bad_report_input")
})

test_that("gallery restores graphics device after a drawing error", {
  skip_if_not_installed("ggplot2")
  skip_if_not(capabilities("cairo"))
  p <- ggplot2::ggplot(data.frame(x = 1), ggplot2::aes(x, missing_y)) +
    ggplot2::geom_point()
  device <- grDevices::dev.cur()
  file <- tempfile(fileext = ".html")
  expect_error(sb_render_visual_report(list(Broken = p), file))
  expect_identical(grDevices::dev.cur(), device)
  expect_identical(file.exists(file), FALSE)
})

test_that("grouped reports retain every figure and move explanations outside images", {
  skip_if_not_installed("ggplot2")
  skip_if_not(capabilities("cairo"))
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 3:1), ggplot2::aes(x, y)) +
    ggplot2::geom_point() + ggplot2::labs(title = "Native result", subtitle = "Saved fit",
      caption = "Original scale < not probability")
  original <- list(data = p$data, labels = p$labels, built = ggplot2::ggplot_build(p)$data)
  file <- tempfile(fileext = ".html")
  sb_render_visual_report(list(A = p, B = p, C = p), file,
    groups = c("Clinical", "Algorithm < details", "Clinical"),
    captions = c("Additional legend", "", ""), width = 4, height = 3)
  html <- paste(readLines(file), collapse = "\n")
  expect_match(html, "<details id='topic-1' open>", fixed = TRUE)
  expect_match(html, "<details id='topic-2'>", fixed = TRUE)
  expect_match(html, "Algorithm &lt; details", fixed = TRUE)
  expect_match(html, "Original scale &lt; not probability", fixed = TRUE)
  expect_match(html, "Additional legend", fixed = TRUE)
  expect_equal(length(gregexpr("data:image/svg+xml", html, fixed = TRUE)[[1]]), 3L)
  expect_identical(p$data, original$data)
  expect_identical(p$labels, original$labels)
  expect_identical(ggplot2::ggplot_build(p)$data, original$built)
  file_all <- tempfile(fileext = ".html")
  sb_render_visual_report(list(A = p, B = p), file_all, groups = c("A", "B"), expanded = c("A", "B"))
  all_html <- paste(readLines(file_all), collapse = "\n")
  expect_match(all_html, "<details id='topic-2' open>", fixed = TRUE)
  expect_error(sb_render_visual_report(list(A = p), tempfile(fileext = ".html"), groups = "A", expanded = "missing"),
    class = "stormbreaker_bad_report_input")
  expect_error(sb_render_visual_report(list(A = p), tempfile(fileext = ".html"), captions = c("a", "b")),
    class = "stormbreaker_bad_report_input")
})
})

# ---- report-export ----
local({
# ---- render-model-report ----

test_that("offline report reuses SHAP with accessible model and sample controls", {
  skip_if_not_installed("ggplot2")
  skip_if_not(capabilities("cairo"))
  d <- make_sb_binary_fixture(n = 50L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  new <- d[1:2, c("patient_id", "gene_a", "gene_b")]
  new$patient_id <- c("private_identifier_one", "private_identifier_two")
  file <- tempfile(fileext = ".html")
  path <- sb_render_model_report(fit, new, "patient_id", d, file,
    title = "Model <script>untrusted</script>"
  )
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(html, "Model &lt;script&gt;untrusted&lt;/script&gt;", fixed = TRUE)
  expect_match(html, "id='coefficient-view'", fixed = TRUE)
  expect_match(html, "id='shap-view'", fixed = TRUE)
  expect_match(html, "id='sample-2'", fixed = TRUE)
  expect_match(html, "replaceChildren", fixed = TRUE)
  expect_match(html, "Training reference: 50", fixed = TRUE)
  expect_match(html, "<table>", fixed = TRUE)
  expect_identical(grepl("private_identifier", html, fixed = TRUE), FALSE)
  expect_identical(grepl("src='https", html, fixed = TRUE), FALSE)
  expect_error(sb_render_model_report(fit, new, "patient_id", d, file),
    class = "stormbreaker_report_exists"
  )
})

test_that("survival report binds the same sample to waterfall and survival", {
  skip_if_not_installed("ggplot2")
  skip_if_not(capabilities("cairo"))
  d <- make_sb_fixture(n = 60L, seed = 107L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 6
  )
  fit <- sb_fit_final(task, sb_pipeline())
  file <- tempfile(fileext = ".html")
  path <- sb_render_model_report(fit, d[1:2, ], "patient_id", d, file, development = task)
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(html, "Predicted survival for observation 2", fixed = TRUE)
  expect_match(html, "SHAP waterfall for observation 2", fixed = TRUE)
  expect_match(html, "Survival at 6 month", fixed = TRUE)
  expect_match(html, "Baseline survival alone", fixed = TRUE)
  expect_error(sb_render_model_report(fit, d, "patient_id", d, tempfile(fileext = ".html")),
    class = "stormbreaker_missing_survival_curve"
  )
})

# ---- render-report ----

test_that("reports preserve boundaries and omit private source data", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("binary")
  x$external$scored_predictions$patient_id[1] <- "PRIVATE_PATIENT_SENTINEL"
  x$external$error_message <- "PRIVATE_ERROR_SENTINEL"
  before <- x
  path <- tempfile(fileext = ".html")
  result <- sb_render_report(list(Development = x$cv, Validation = x$external), path,
                             title = "My <report> & evidence")
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_identical(result, normalizePath(path, winslash = "/"))
  for (text in c("Development CV", "External evaluation", "Calibration", "Decision curves",
                 "data:image/svg+xml", "My &lt;report&gt; &amp; evidence", "roc_auc")) {
    expect_match(html, text, fixed = TRUE)
  }
  for (text in c("PRIVATE_PATIENT_SENTINEL", "PRIVATE_ERROR_SENTINEL", "<script", "<iframe")) {
    expect_identical(grepl(text, html, fixed = TRUE), FALSE)
  }
  expect_identical(x, before)
})

test_that("reports never overwrite an existing file or accept ambiguous sections", {
  x <- make_sb_curve_view_fixture("binary")$cv
  path <- tempfile(fileext = ".html")
  writeLines("PRESERVE_THIS_FILE", path)
  expect_error(sb_render_report(x, path), class = "stormbreaker_report_exists")
  expect_identical(readLines(path), "PRESERVE_THIS_FILE")
  for (bad in list(list(x), list(), structure(list(x,x), names=c("Same","Same")), 1)) {
    target <- tempfile(fileext = ".html")
    expect_error(sb_render_report(bad, target), class = "stormbreaker_bad_report_input")
    expect_identical(file.exists(target), FALSE)
  }
})

test_that("survival report retains saved horizon and failed curves", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("survival", TRUE)$external
  x$evaluations[[1]]$status <- "failed"
  x$evaluations[[1]]$metrics <- list()
  x$evaluations[[1]]$error_message <- "SECRET_FAILURE_DETAIL"
  path <- tempfile(fileext = ".html")
  sb_render_report(x, path)
  html <- paste(readLines(path, warn=FALSE),collapse="\n")
  for (text in c("harrell_c", "prediction_horizon", "month", "not_saved", "failed")) {
    expect_match(html,text,fixed=TRUE)
  }
  expect_identical(grepl("SECRET_FAILURE_DETAIL",html,fixed=TRUE),FALSE)
})

test_that("report rendering cannot fit or evaluate models", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("binary")$external
  local_mocked_bindings(
    sb_fit_final=function(...) stop("must not fit"),
    sb_run_cv=function(...) stop("must not run CV"),
    sb_evaluate_external=function(...) stop("must not evaluate"),
    sb_external_binary_auc_metric=function(...) stop("must not compute")
  )
  path <- tempfile(fileext=".html")
  sb_render_report(x,path)
  expect_identical(file.exists(path),TRUE)
})

test_that("survival CV without requested horizon has explicit missing curves", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_survival(60, seed=1)
  task <- sb_surv_task(d,"patient_id","time","event","gene_a","month")
  x <- sb_run_cv(task,folds=3,seed=1)
  path <- tempfile(fileext=".html")
  sb_render_report(x,path)
  html <- paste(readLines(path,warn=FALSE),collapse="\n")
  expect_match(html,"not_requested",fixed=TRUE)
  expect_match(html,"Development CV",fixed=TRUE)
  expect_match(html,"harrell_c",fixed=TRUE)
})

test_that("report tables faithfully format saved data without raw values hidden", {
  tab <- data.frame(pipeline_id="<unsafe&label>", estimate=c(0.123456789),
                    status="success", raw_patient="HIDDEN_PATIENT")
  html <- sb_report_table(tab,c("pipeline_id","estimate","status"),"Saved <metric>")
  expect_match(html,"0.12345679",fixed=TRUE)
  expect_match(html,"&lt;unsafe&amp;label&gt;",fixed=TRUE)
  expect_match(html,"Saved &lt;metric&gt;",fixed=TRUE)
  expect_identical(grepl("HIDDEN_PATIENT",html,fixed=TRUE),FALSE)
})

test_that("report writes separate failures and excludes failed pipeline names from HTML", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("binary",TRUE)$external
  x$evaluations[[1]]$status <- "failed"
  x$evaluations[[1]]$error_message <- 'PRIVATE reason, with "quotes"'
  before <- x
  path <- tempfile(fileext=".html")
  sb_render_report(x,path)
  csv <- sub("\\.html$","_failures.csv",path)
  html <- paste(readLines(path,warn=FALSE),collapse="\n")
  issues <- read.csv(csv,stringsAsFactors=FALSE)
  expect_identical(grepl(x$evaluations[[1]]$pipeline_id,html,fixed=TRUE),FALSE)
  expect_match(html,"Download failures CSV",fixed=TRUE)
  expect_identical(grepl("PRIVATE reason",html,fixed=TRUE),FALSE)
  expect_identical(issues$error_message[issues$issue_type=="pipeline"],x$evaluations[[1]]$error_message)
  expect_identical(nrow(issues),nrow(sb_failure_table(list(Result=x))))
  expect_identical(x,before)
})

test_that("existing companion CSV prevents either report output from being overwritten", {
  x <- make_sb_curve_view_fixture("binary")$external
  path <- tempfile(fileext=".html")
  csv <- sub("\\.html$","_failures.csv",path)
  writeLines("PRESERVE_CSV",csv)
  expect_error(sb_render_report(x,path),class="stormbreaker_output_exists")
  expect_identical(file.exists(path),FALSE)
  expect_identical(readLines(csv),"PRESERVE_CSV")
})

test_that("all-failed report explicitly states no usable results and exports failures", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("binary")$external
  x$status <- "failed"; x$metrics <- list()
  path <- tempfile(fileext=".html")
  sb_render_report(x,path)
  html <- paste(readLines(path,warn=FALSE),collapse="\n")
  expect_match(html,"No usable results in this section",fixed=TRUE)
  expect_match(html,"This is not a successful analysis",fixed=TRUE)
  expect_identical(grepl(x$pipeline_id,html,fixed=TRUE),FALSE)
  csv <- read.csv(sub("\\.html$","_failures.csv",path),stringsAsFactors=FALSE)
  expect_identical(csv$pipeline_id,x$pipeline_id)
  expect_identical(csv$status,"failed")
})

# ---- report-compact ----

test_that('compact model reports keep numeric semantics but remove repeated prose', {
  skip_if_not_installed('ggplot2')
  d <- make_sb_binary_fixture(n=50L,seed=108L)
  task <- sb_binary_task(d,'patient_id','outcome',c('gene_a','gene_b'))
  fit <- sb_fit_final(task,sb_pipeline(learner='logistic'))
  before <- fit
  file <- tempfile(fileext='.html')
  sb_render_model_report(fit,d[1:2,],'patient_id',d,file)
  html <- paste(readLines(file,warn=FALSE),collapse='\n')
  decoded <- utils::URLdecode(html)
  expect_identical(grepl('Exact empirical interventional contributions',decoded,fixed=TRUE),FALSE)
  expect_identical(grepl('Verify critical decisions',html,fixed=TRUE),FALSE)
  expect_identical(grepl('Coefficient magnitude is not',html,fixed=TRUE),FALSE)
  expect_match(html,'Intercept:',fixed=TRUE)
  expect_match(html,'log odds (class 1)',fixed=TRUE)
  expect_match(html,'Training reference: 50',fixed=TRUE)
  expect_match(html,'Class-1 probability:',fixed=TRUE)
  expect_identical(fit,before)
  detailed <- tempfile(fileext='.html')
  sb_render_model_report(fit,d[1:2,],'patient_id',d,detailed,compact=FALSE)
  detailed_html <- utils::URLdecode(paste(readLines(detailed,warn=FALSE),collapse='\n'))
  expect_identical(grepl('>Individual prediction</text>',decoded,fixed=TRUE),FALSE)
  expect_identical(grepl('>Individual prediction</text>',detailed_html,fixed=TRUE),TRUE)
})

test_that('gene counts appear only when saved development selections are supplied', {
  skip_if_not_installed('ggplot2')
  d <- make_sb_binary_fixture(n=50L,seed=108L)
  task <- sb_binary_task(d,'patient_id','outcome',c('gene_a','gene_b'))
  fit <- sb_fit_final(task,sb_pipeline(learner='logistic'))
  cv <- sb_run_cv(task,sb_pipeline(selector='univ_binary_auc_top_k',selector_n_features=1,learner='logistic'),folds=2,seed=5)
  before <- cv
  file <- tempfile(fileext='.html')
  sb_render_model_report(fit,d[1:2,],'patient_id',d,file,selection_results=cv)
  html <- paste(readLines(file,warn=FALSE),collapse='\n')
  expect_match(html,"id='gene-view'",fixed=TRUE)
  expect_match(html,"id='gene-selector'",fixed=TRUE)
  expect_match(html,"id='gene-glmnet_active_set'",fixed=TRUE)
  expect_match(html,'All genes (2)',fixed=TRUE)
  expect_match(html,'shared selector copies count once',fixed=TRUE)
  expect_identical(cv,before)
})
})

# ---- stacking ----
local({
# ---- stack-prediction-variation ----

test_that("stack comparison reports actual prediction variation for both tasks", {
  skip_if_not_installed("glmnet")
  for (type in c("binary", "survival")) {
    if (type == "binary") {
      d <- make_sb_binary_fixture(n = 60L, seed = 101L)
      task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
      learners <- c("logistic", "gaussian_naive_bayes")
    } else {
      d <- make_sb_fixture(n = 72L, seed = 103L)
      task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
      learners <- c("coxph", "survreg_weibull")
    }
    stack <- sb_stack(lapply(learners, function(id) sb_pipeline(learner = id)),
      inner_folds = 3L, inner_seed = 31L)
    result <- sb_run_stack_cv(task, stack, folds = 3L, seed = 37L)
    expect_identical(result$status, "success")
    expect_identical(tail(result$comparison$prediction_variation, 1L),
      sb_prediction_variation(result$predictions, type)$status)
    value <- if (type == "binary") "probability" else "risk"
    result$predictions[[value]][] <- 0.5
    expect_identical(sb_stack_comparison_row(result)$prediction_variation, "constant")
    result$predictions$status[1L] <- "failed"
    expect_identical(sb_stack_comparison_row(result)$prediction_variation, "incomplete")
  }
})

# ---- stacking-cv ----

test_that("binary stacking learns only from inner OOF member predictions", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 60L, seed = 101L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c"),
    decision_thresholds = c(0.2, 0.5)
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 31L)

  result <- sb_run_stack_cv(task, stack, folds = 3L, seed = 37L)

  expect_s3_class(result, "stormbreaker_stack_cv_result")
  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(result$predictions$patient_id, data$patient_id)
  expect_true(all(is.finite(result$predictions$probability)))
  expect_identical(
    result$provenance$fold_assignment,
    result$members$fold_assignment
  )
  expect_length(result$inner_folds, 3L)
  for (audit in result$inner_folds) {
    analysis_ids <- result$folds[[audit$fold]]$analysis_ids
    expect_identical(rownames(audit$member_matrix), analysis_ids)
    expect_identical(names(audit$member_matrix), stack$pipeline_ids)
    for (pipeline_id in stack$pipeline_ids) {
      predictions <- audit$comparison$results[[pipeline_id]]$predictions
      expect_setequal(predictions$patient_id, analysis_ids)
      expect_identical(anyDuplicated(predictions$patient_id), 0L)
      expect_identical(all(predictions$status == "success"), TRUE)
    }
  }
  expect_identical(
    tail(result$comparison$pipeline_id, 1L),
    stack$id
  )
  expect_named(
    result$metrics,
    c(
      "roc_auc", "average_precision", "brier_score", "log_loss",
      "calibration", "moderate_calibration", "decision_curve"
    )
  )
})

test_that("outer assessment outcomes cannot alter their stacked predictions", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 48L, seed = 102L)
  assignment <- setNames(rep(1:3, length.out = nrow(data)), data$patient_id)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 41L)

  original <- sb_run_stack_cv(task, stack, folds = assignment, seed = 43L)
  changed <- data
  changed_rows <- assignment[data$patient_id] == 1L
  changed$outcome[changed_rows] <- 1L - changed$outcome[changed_rows]
  changed_task <- sb_binary_task(
    changed,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  perturbed <- sb_run_stack_cv(
    changed_task,
    stack,
    folds = assignment,
    seed = 43L
  )

  expect_equal(
    original$predictions$probability[changed_rows],
    perturbed$predictions$probability[changed_rows],
    tolerance = 0
  )
  expect_identical(
    original$inner_folds[[1L]]$member_matrix,
    perturbed$inner_folds[[1L]]$member_matrix
  )
  expect_identical(
    original$inner_folds[[1L]]$meta$state,
    perturbed$inner_folds[[1L]]$meta$state
  )
})

test_that("survival stacking returns risk and locked-horizon evaluation", {
  skip_if_not_installed("glmnet")

  data <- make_sb_fixture(n = 72L, seed = 103L)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 2,
    decision_thresholds = c(0.2, 0.4)
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "coxph"),
    sb_pipeline(learner = "survreg_weibull")
  ), inner_folds = 3L, inner_seed = 47L)

  result <- sb_run_stack_cv(task, stack, folds = 3L, seed = 53L)

  expect_s3_class(result, "stormbreaker_stack_cv_result")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(result$complete_survival_probability_oof, TRUE)
  expect_true(all(is.finite(result$predictions$risk)))
  expect_true(all(
    result$predictions$survival_probability >= 0 &
      result$predictions$survival_probability <= 1
  ))
  expect_named(
    result$metrics,
    c(
      "harrell_c", "ipcw_brier", "time_auc", "uno_c",
      "calibration", "weak_calibration", "decision_curve"
    )
  )
  expect_identical(
    vapply(
      result$inner_folds,
      function(audit) audit$meta$state$baseline_hazard$prediction_horizon,
      numeric(1)
    ),
    rep(2, 3L)
  )
})

test_that("stack folds fail visibly when the ridge dependency is unavailable", {
  data <- make_sb_binary_fixture(n = 48L, seed = 104L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L)
  local_mocked_bindings(sb_glmnet_available = function() FALSE)

  result <- sb_run_stack_cv(task, stack, folds = 3L, seed = 59L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$error_class), "stormbreaker_missing_glmnet")
  expect_true(all(is.na(result$predictions$probability)))
  expect_identical(names(result$members$results), stack$pipeline_ids)
  expect_identical(nrow(result$comparison), length(stack$pipelines) + 1L)
})

# ---- stacking-final ----

test_that("binary stack final fit predicts without external outcomes", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 60L, seed = 111L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c"),
    decision_thresholds = c(0.2, 0.5)
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 61L)

  fit <- sb_fit_stack_final(task, stack)
  new_data <- data[1:8, c("patient_id", "gene_a", "gene_b", "gene_c")]
  prediction <- sb_predict_stack_final(fit, new_data, "patient_id")
  expected <- stats::plogis(
    fit$meta$state$intercept +
      as.vector(
        as.matrix(prediction$member_predictions) %*%
          fit$meta$state$coefficient
      )
  )

  expect_s3_class(fit, "stormbreaker_stack_final_fit")
  expect_identical(fit$status, "success")
  expect_identical(names(fit$member_fits), stack$pipeline_ids)
  expect_identical(fit$meta$state$features, stack$pipeline_ids)
  expect_identical(fit$required_features, task$predictors)
  expect_s3_class(prediction, "stormbreaker_final_prediction")
  expect_identical(prediction$status, "success")
  expect_identical(prediction$predictions$patient_id, new_data$patient_id)
  expect_equal(prediction$predictions$probability, expected, tolerance = 1e-14)
  expect_identical(prediction$decision_thresholds, c(0.2, 0.5))

  endpoint <- sb_binary_endpoint(
    data[1:8, c("patient_id", "outcome")],
    "patient_id",
    "outcome"
  )
  evaluation <- sb_evaluate_external(prediction, endpoint)
  expect_s3_class(evaluation, "stormbreaker_external_evaluation")
  expect_identical(evaluation$status, "success")
})

test_that("named member lists fit and predict like ordinary stack inputs", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 48L, seed = 117L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    first = sb_pipeline(learner = "logistic"),
    second = sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 73L)

  fit <- sb_fit_stack_final(task, stack)
  prediction <- sb_predict_stack_final(
    fit,
    data[1:6, c("patient_id", "gene_a", "gene_b")],
    "patient_id"
  )

  expect_null(names(stack$pipeline_ids))
  expect_identical(fit$status, "success")
  expect_identical(prediction$status, "success")
})

test_that("serialized stack fits reproduce binary predictions", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 54L, seed = 112L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 67L)
  fit <- sb_fit_stack_final(task, stack)
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(fit, path)
  restored <- readRDS(path)
  new_data <- data[1:6, c("patient_id", "gene_a", "gene_b")]

  original <- sb_predict_stack_final(fit, new_data, "patient_id")
  replayed <- sb_predict_stack_final(restored, new_data, "patient_id")

  expect_identical(replayed, original)
})

test_that("survival stack final prediction reuses its saved OOF baseline", {
  skip_if_not_installed("glmnet")

  data <- make_sb_fixture(n = 72L, seed = 113L)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 2,
    decision_thresholds = c(0.2, 0.4)
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "coxph"),
    sb_pipeline(learner = "survreg_weibull")
  ), inner_folds = 3L, inner_seed = 71L)

  fit <- sb_fit_stack_final(task, stack)
  new_data <- data[, c(
    "patient_id", "gene_a", "gene_b", "gene_c"
  )]
  prediction <- sb_predict_stack_final(fit, new_data, "patient_id")
  expected_risk <- as.vector(
    as.matrix(prediction$member_predictions) %*%
      fit$meta$state$coefficient
  )
  expected_probability <- sb_survival_probability_from_hazard(
    expected_risk,
    fit$meta$state$baseline_hazard$cumulative_baseline_hazard
  )

  expect_identical(fit$status, "success")
  expect_identical(prediction$status, "success")
  expect_equal(prediction$predictions$risk, expected_risk, tolerance = 1e-14)
  expect_equal(
    prediction$predictions$survival_probability,
    expected_probability,
    tolerance = 1e-14
  )
  expect_identical(
    prediction$predictions$prediction_horizon,
    rep(2, nrow(new_data))
  )

  endpoint <- sb_surv_endpoint(
    data[, c("patient_id", "time", "event")],
    "patient_id",
    "time",
    "event",
    "month",
    prediction_horizon = 2
  )
  evaluation <- sb_evaluate_external(prediction, endpoint)
  expect_s3_class(evaluation, "stormbreaker_external_evaluation")
  expect_identical(evaluation$status, "success")
})

test_that("final stack fitting reports a missing ridge dependency", {
  data <- make_sb_binary_fixture(n = 48L, seed = 114L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L)
  local_mocked_bindings(sb_glmnet_available = function() FALSE)

  fit <- sb_fit_stack_final(task, stack)

  expect_identical(fit$status, "failed")
  expect_identical(fit$error_stage, "meta_fit")
  expect_identical(fit$error_class, "stormbreaker_missing_glmnet")
  expect_length(fit$member_fits, 0L)
})

test_that("a failed member prediction stops the stack without fallback", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 48L, seed = 115L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L)
  fit <- sb_fit_stack_final(task, stack)
  calls <- 0L
  local_mocked_bindings(
    sb_predict_final = function(...) {
      calls <<- calls + 1L
      list(
        status = "failed",
        error_class = "stormbreaker_synthetic_member_failure",
        error_message = "synthetic member failure"
      )
    }
  )
  new_data <- data[1:4, c("patient_id", "gene_a", "gene_b")]

  prediction <- sb_predict_stack_final(fit, new_data, "patient_id")

  expect_identical(prediction$status, "failed")
  expect_identical(prediction$error_stage, "member_prediction")
  expect_identical(
    prediction$error_class,
    "stormbreaker_stack_member_prediction_failure"
  )
  expect_identical(calls, 1L)
})

# ---- stacking-meta ----

test_that("member matrices realign binary OOF predictions by patient ID", {
  data <- make_sb_binary_fixture(n = 48L, seed = 91L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 17L)
  comparison <- sb_compare_cv(
    task,
    stack$pipelines,
    folds = 3L,
    seed = 17L
  )
  second_id <- stack$pipeline_ids[[2L]]
  comparison$results[[second_id]]$predictions <-
    comparison$results[[second_id]]$predictions[48:1, , drop = FALSE]

  matrix <- sb_stack_member_matrix(comparison, task, stack)

  expect_s3_class(matrix, "data.frame")
  expect_identical(names(matrix), stack$pipeline_ids)
  expect_identical(rownames(matrix), data$patient_id)
  for (pipeline_id in stack$pipeline_ids) {
    predictions <- comparison$results[[pipeline_id]]$predictions
    expected <- predictions$probability[
      match(data$patient_id, predictions$patient_id)
    ]
    expect_equal(matrix[[pipeline_id]], expected, tolerance = 0)
  }
})

test_that("member matrices reject incomplete OOF predictions", {
  data <- make_sb_binary_fixture(n = 48L, seed = 92L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L)
  comparison <- sb_compare_cv(task, stack$pipelines, folds = 3L, seed = 9L)
  pipeline_id <- stack$pipeline_ids[[1L]]

  duplicate <- comparison
  duplicate$results[[pipeline_id]]$predictions$patient_id[[1L]] <-
    duplicate$results[[pipeline_id]]$predictions$patient_id[[2L]]
  duplicate_error <- capture_sb_error(
    sb_stack_member_matrix(duplicate, task, stack)
  )
  expect_s3_class(
    duplicate_error,
    "stormbreaker_incomplete_stack_member_oof"
  )

  failed <- comparison
  failed$results[[pipeline_id]]$predictions$status[[1L]] <- "failed"
  failed_error <- capture_sb_error(
    sb_stack_member_matrix(failed, task, stack)
  )
  expect_s3_class(
    failed_error,
    "stormbreaker_incomplete_stack_member_oof"
  )

  nonfinite <- comparison
  nonfinite$results[[pipeline_id]]$predictions$probability[[1L]] <- Inf
  nonfinite_error <- capture_sb_error(
    sb_stack_member_matrix(nonfinite, task, stack)
  )
  expect_s3_class(
    nonfinite_error,
    "stormbreaker_incomplete_stack_member_oof"
  )
})

test_that("ridge logistic stack state matches direct coefficient prediction", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 60L, seed = 93L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 23L)
  comparison <- sb_compare_cv(
    task,
    stack$pipelines,
    folds = 3L,
    seed = stack$inner_seed
  )
  matrix <- sb_stack_member_matrix(comparison, task, stack)

  fitted <- sb_fit_stack_meta(task, stack, matrix)
  observed <- sb_predict_stack_meta(fitted$state, matrix)
  expected <- stats::plogis(
    fitted$state$intercept +
      as.vector(as.matrix(matrix) %*% fitted$state$coefficient)
  )

  expect_s3_class(fitted, "stormbreaker_stack_meta_fit")
  expect_identical(fitted$meta_learner, "ridge_logistic")
  expect_identical(fitted$tuning$alpha, 0)
  expect_identical(fitted$tuning$selection_rule, "lambda.min")
  expect_named(observed, "probability")
  expect_equal(observed$probability, expected, tolerance = 1e-14)
  expect_identical(fitted$state$features, stack$pipeline_ids)
})

test_that("ridge Cox stack state saves its OOF Breslow horizon", {
  skip_if_not_installed("glmnet")

  data <- make_sb_fixture(n = 60L, seed = 94L)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 2
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "coxph"),
    sb_pipeline(learner = "rpart_survival")
  ), inner_folds = 3L, inner_seed = 29L)
  comparison <- sb_compare_cv(
    task,
    stack$pipelines,
    folds = 3L,
    seed = stack$inner_seed
  )
  matrix <- sb_stack_member_matrix(comparison, task, stack)

  fitted <- sb_fit_stack_meta(task, stack, matrix)
  observed <- sb_predict_stack_meta(fitted$state, matrix)
  expected_risk <- as.vector(
    as.matrix(matrix) %*% fitted$state$coefficient
  )
  expected_probability <- sb_survival_probability_from_hazard(
    expected_risk,
    fitted$state$baseline_hazard$cumulative_baseline_hazard
  )

  expect_s3_class(fitted, "stormbreaker_stack_meta_fit")
  expect_identical(fitted$meta_learner, "ridge_cox")
  expect_identical(fitted$tuning$alpha, 0)
  expect_equal(observed$risk, expected_risk, tolerance = 1e-14)
  expect_equal(
    observed$survival_probability,
    expected_probability,
    tolerance = 1e-14
  )
  expect_identical(observed$prediction_horizon, 2)
  expect_identical(observed$time_unit, "month")
  expect_identical(
    fitted$state$baseline_hazard$source,
    "glmnet_fixed_beta_breslow"
  )
})

# ---- stacking-spec ----

test_that("a stack preserves distinct predeclared binary members", {
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )

  stack <- sb_stack(pipelines, inner_folds = 3L, inner_seed = 41L)

  expect_s3_class(stack, "stormbreaker_stack")
  expect_identical(stack$task_type, "binary")
  expect_identical(stack$pipelines, pipelines)
  expect_identical(stack$pipeline_ids, vapply(
    pipelines,
    function(pipeline) pipeline$id,
    character(1)
  ))
  expect_identical(stack$meta_learner, "ridge_logistic")
  expect_identical(stack$inner_folds, 3L)
  expect_identical(stack$inner_seed, 41L)
  expect_identical(
    stack$id,
    paste0("ridge_stack__", paste(stack$pipeline_ids, collapse = "__"))
  )
  expect_identical(sb_validate_stack_contract(stack, "binary"), stack)
  expect_match(capture_output(print(stack)), "Members: 2")
})

test_that("a survival stack declares ridge Cox meta semantics", {
  stack <- sb_stack(list(
    sb_pipeline(learner = "coxph"),
    sb_pipeline(learner = "rpart_survival")
  ), inner_folds = 3L)

  expect_identical(stack$task_type, "survival")
  expect_identical(stack$meta_learner, "ridge_cox")
  expect_identical(stack$meta_alpha, 0)
  expect_identical(stack$meta_standardize, TRUE)
  expect_identical(stack$meta_selection_rule, "lambda.min")
  expect_identical(stack$meta_measure, "deviance")
})

test_that("stack construction rejects ambiguous member collections", {
  one_member <- capture_sb_error(sb_stack(list(
    sb_pipeline(learner = "logistic")
  )))
  expect_s3_class(one_member, "stormbreaker_bad_stack_pipelines")

  binary <- sb_pipeline(learner = "logistic")
  duplicate <- capture_sb_error(sb_stack(list(binary, binary)))
  expect_s3_class(
    duplicate,
    "stormbreaker_duplicate_stack_pipeline_id"
  )

  mixed <- capture_sb_error(sb_stack(list(
    binary,
    sb_pipeline(learner = "coxph")
  )))
  expect_s3_class(mixed, "stormbreaker_stack_task_mismatch")

  not_pipelines <- capture_sb_error(sb_stack(list(list(), list())))
  expect_s3_class(not_pipelines, "stormbreaker_bad_stack_pipelines")

  bad_folds <- capture_sb_error(sb_stack(
    list(
      sb_pipeline(learner = "logistic"),
      sb_pipeline(learner = "gaussian_naive_bayes")
    ),
    inner_folds = 2L
  ))
  expect_s3_class(bad_folds, "stormbreaker_bad_integer")
})

test_that("task slicing preserves binary roles and exact patient order", {
  data <- make_sb_binary_fixture(n = 40L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b"),
    decision_thresholds = c(0.2, 0.5)
  )
  rows <- sort(c(
    head(which(data$outcome == 0L), 8L),
    head(which(data$outcome == 1L), 8L)
  ))

  sliced <- sb_subset_task(task, rows)

  expect_s3_class(sliced, "stormbreaker_binary_task")
  expect_identical(sliced$data, data[rows, , drop = FALSE])
  expect_identical(sliced$predictors, task$predictors)
  expect_identical(sliced$decision_thresholds, task$decision_thresholds)
})

test_that("task slicing preserves survival horizon metadata", {
  data <- make_sb_fixture(n = 40L)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 2,
    decision_thresholds = c(0.2, 0.4)
  )
  rows <- sort(c(
    head(which(data$event == 0L), 8L),
    head(which(data$event == 1L), 8L)
  ))

  sliced <- sb_subset_task(task, rows)

  expect_s3_class(sliced, "stormbreaker_surv_task")
  expect_identical(sliced$data, data[rows, , drop = FALSE])
  expect_identical(sliced$time_unit, "month")
  expect_identical(sliced$prediction_horizon, 2)
  expect_identical(sliced$decision_thresholds, c(0.2, 0.4))
})
})
