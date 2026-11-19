test_that("route quick start compares both routes and reuses the selected fit", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ggplot2")
  env <- new.env(parent = globalenv())
  sys.source(system.file("examples", "scenarios", "07-route-workflow.R",
    package = "stormbreaker", mustWork = TRUE), envir = env)
  x <- env$example_result

  expect_identical(x$cv_auc$route, c("Direct", "Lasso"))
  expect_identical(x$cv_auc$status, rep("success", 2))
  runs <- result_table(x$run$comparison, "runs")
  expect_identical(runs$status, rep("success", 6))
  expect_identical(x$run$comparison$results[[1]]$provenance$fold_assignment,
    x$run$comparison$results[[2]]$provenance$fold_assignment)
  selected <- x$cv_auc$route[which.max(x$cv_auc$estimate)]
  expect_identical(x$finished$fit$pipeline$id, x$routes[[selected]]$id)
  expect_identical(x$reloaded$pipeline$id, x$finished$fit$pipeline$id)
  expect_identical(x$test_auc$status, "success")
  expect_false("outcome" %in% names(x$test_predictors))
  expect_equal(nrow(x$prediction$predictions), nrow(x$test_predictors))
  expect_identical(x$prediction$predictions, x$finished$predictions$Test)
  expect_s3_class(x$figures$ROC, "ggplot")
})
