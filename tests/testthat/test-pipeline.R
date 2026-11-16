# ---- candidate-menus ----
local({
# ---- binary-menu-pipelines ----

test_that("binary menu covers the registry without fitting or changing defaults", {
  genes <- paste0("gene_", seq_len(55L))
  menu <- sb_binary_menu_pipelines(genes, "fixed55")
  registry <- sb_pipeline_registry("binary")
  expect_setequal(paste(menu$selector, menu$learner),
    paste(registry$selector, registry$learner))
  expect_identical(nrow(menu), nrow(registry))
  expect_identical(anyDuplicated(menu$config_id), 0L)
  expect_identical(anyDuplicated(vapply(menu$pipeline, `[[`, "", "id")), 0L)
  expect_identical(attr(menu, "registry_pairs"), nrow(registry))
  for (p in menu$pipeline) {
    expect_identical(p$task_type, "binary")
    if (p$selector == "fixed_features") {
      expect_identical(p$selector_params$features, genes)
      expect_identical(p$selector_params$feature_set_id, "fixed55")
    }
    if (p$selector == "univ_binary_auc_top_k") {
      expect_identical(p$selector_params$n_features, 10L)
    }
    if (p$learner == "glmnet_logistic") {
      expect_identical(p$learner_params$inner_folds, 5L)
      expect_identical(p$learner_params$inner_seed, 1L)
    }
  }
  reference <- attr(menu, "mime7_family_mapping")
  expect_setequal(reference$official_method,
    c("nb", "svmRadialWeights", "rf", "kknn", "adaboost", "LogitBoost", "cancerclass"))
  expect_identical(anyNA(match(reference$config_id, menu$config_id)), FALSE)
  expect_identical(unique(reference$relationship), "algorithm_family_only")
  expect_identical(sum(nzchar(menu$mime_family)), 7L)
})

test_that("binary menu controls are explicit and invalid inputs are rejected", {
  invalid <- tryCatch(sb_binary_menu_pipelines(c("a", "a"), "duplicate"), error = identity)
  expect_s3_class(invalid, "stormbreaker_bad_fixed_features")
  invalid <- tryCatch(sb_binary_menu_pipelines(c("a", "b"), "valid", inner_folds = 2L), error = identity)
  expect_s3_class(invalid, "error")
  row <- sb_pipeline_registry("binary")
  row <- row[row$selector == "none" & row$learner == "glmnet_logistic", ]
  args <- sb_binary_menu_arguments(row, c("a", "b"), "fixed", 2L, 3L, 42L)
  p <- do.call(sb_pipeline, args)
  expect_identical(p$learner_params$inner_folds, 3L)
  expect_identical(p$learner_params$inner_seed, 42L)
  direct <- sb_pipeline(learner = "glmnet_logistic", inner_folds = 3L, inner_seed = 42L)
  expect_identical(p, direct)
})

# ---- candidate-registry ----

test_that("candidate registry preserves the retained historical catalog", {
  registry <- sb_candidate_registry()

  expect_s3_class(registry, "data.frame")
  expect_identical(nrow(registry), 1175L)
  task_counts <- table(registry$runtime_task_type)
  expect_identical(names(task_counts), c("binary", "survival"))
  expect_identical(as.integer(task_counts), c(584L, 591L))
  expect_identical(anyDuplicated(registry$pipeline_id), 0L)
  expect_identical(registry$pipeline_id[[1L]], "surv_0001")
  expect_identical(registry$pipeline_id[[1175L]], "class_0584")
  expect_identical(unique(registry$adapter_status), "not_verified")
  expect_identical(any(registry$execution_status == "executable"), FALSE)
  expect_identical(
    all(grepl("^(not_admitted|blocked)", registry$admission_status)),
    TRUE
  )
  expect_identical(
    any(registry$runtime_admission_status == "admitted"),
    FALSE
  )
  expect_identical(sum(registry$runtime_route_available), 138L)
})

test_that("candidate registry maps only explicit current runtime roles", {
  registry <- sb_candidate_registry()
  row_by_id <- function(id) {
    registry[registry$pipeline_id == id, , drop = FALSE]
  }

  logistic <- row_by_id("class_0001")
  adaboost <- row_by_id("class_0029")
  stepcox <- row_by_id("surv_0013")
  xgboost <- row_by_id("class_0024")
  decision_tree <- registry[
    registry$task_type == "binary_classification" &
      registry$pipeline_label == "DecisionTree" &
      registry$model_type == "single", ,
    drop = FALSE
  ]
  stepcox_ridge <- registry[
    registry$task_type == "survival_right_censored" &
      registry$pipeline_label == "StepCox[forward] + Ridge", ,
    drop = FALSE
  ]
  blocked <- registry[grepl("^blocked", registry$admission_status), ,
    drop = FALSE
  ][1L, , drop = FALSE]
  unresolved <- registry[
    registry$role_semantics_status == "stage_roles_unresolved", ,
    drop = FALSE
  ][1L, , drop = FALSE]
  nested_stepcox <- registry[
    registry$selector_1 == "StepCox[forward]" &
      registry$learner == "StepCox[both]", ,
    drop = FALSE
  ]

  expect_identical(logistic$runtime_selector, "none")
  expect_identical(logistic$runtime_learner, "logistic")
  expect_identical(logistic$runtime_route_status, "route_available")
  expect_identical(
    logistic$runtime_configuration_status,
    "no_additional_runtime_parameters"
  )
  expect_identical(adaboost$runtime_learner, "adaboost_binary")
  expect_identical(
    adaboost$runtime_configuration_status,
    "additional_parameters_required"
  )
  expect_identical(stepcox$runtime_selector, "stepcox")
  expect_identical(stepcox$runtime_learner, "coxph")
  expect_identical(
    stepcox$runtime_parameters,
    "stepcox_direction=forward"
  )
  expect_identical(nrow(stepcox_ridge), 1L)
  expect_identical(stepcox_ridge$runtime_route_status, "route_available")
  expect_identical(stepcox_ridge$runtime_learner, "glmnet_cox")
  expect_identical(
    stepcox_ridge$runtime_parameters,
    "stepcox_direction=forward; glmnet_alpha=0"
  )
  expect_identical(xgboost$runtime_route_status, "route_available")
  expect_identical(xgboost$runtime_learner, "xgboost_binary")
  expect_true(is.na(xgboost$missing_runtime_components))
  expect_identical(nrow(decision_tree), 1L)
  expect_identical(
    decision_tree$runtime_learner,
    "decision_tree_binary"
  )
  expect_identical(decision_tree$runtime_route_status, "route_available")
  expect_identical(blocked$runtime_route_status, "methodologically_blocked")
  expect_identical(blocked$runtime_route_available, FALSE)
  expect_identical(unresolved$runtime_route_status, "stage_roles_unresolved")
  expect_identical(unresolved$runtime_route_available, FALSE)
  expect_identical(nrow(nested_stepcox), 1L)
  expect_identical(
    nested_stepcox$runtime_route_status,
    "unsupported_pipeline_shape"
  )
})

test_that("new specialized binary labels map to real runtime learners", {
  registry <- sb_candidate_registry()
  labels <- c(
    "Stepglm[forward]", "Stepglm[backward]", "Stepglm[both]",
    "svmRadialWeights", "glmBoost", "plsRglm"
  )
  direct <- registry[
    registry$runtime_task_type == "binary" &
      registry$model_type == "single" &
      registry$pipeline_label %in% labels, ,
    drop = FALSE
  ]
  expected_learners <- c(
    `Stepglm[forward]` = "stepglm_binary",
    `Stepglm[backward]` = "stepglm_binary",
    `Stepglm[both]` = "stepglm_binary",
    svmRadialWeights = "svm_radial_weighted_binary",
    glmBoost = "glmboost_binary",
    plsRglm = "plsrglm_binary"
  )

  expect_identical(nrow(direct), 6L)
  expect_identical(sort(direct$pipeline_label), sort(labels))
  expect_identical(
    direct$runtime_learner,
    unname(expected_learners[direct$pipeline_label])
  )
  expect_identical(unique(direct$runtime_selector), "none")
  expect_identical(unique(direct$runtime_route_status), "route_available")
  expect_identical(unique(direct$runtime_route_available), TRUE)
  expect_identical(
    direct$runtime_parameters[grepl("^Stepglm", direct$pipeline_label)],
    paste0(
      "stepglm_direction=",
      sub("^Stepglm\\[(.*)\\]$", "\\1", direct$pipeline_label[
        grepl("^Stepglm", direct$pipeline_label)
      ])
    )
  )
  expect_identical(
    unique(direct$runtime_admission_status),
    "not_admitted_source_to_spec_unverified"
  )
  expect_identical(any(direct$execution_status == "executable"), FALSE)
})

test_that("three mboost survival labels map to distinct real learners", {
  registry <- sb_candidate_registry()
  labels <- c("GamBoost", "GlmBoost", "BlackBoost")
  direct <- registry[
    registry$runtime_task_type == "survival" &
      registry$model_type == "single" &
      registry$pipeline_label %in% labels, ,
    drop = FALSE
  ]
  expected <- c(
    GamBoost = "gamboost_cox",
    GlmBoost = "glmboost_cox",
    BlackBoost = "blackboost_cox"
  )

  expect_identical(nrow(direct), 3L)
  expect_identical(
    direct$runtime_learner,
    unname(expected[direct$pipeline_label])
  )
  expect_identical(unique(direct$runtime_route_status), "route_available")
  expect_identical(unique(direct$runtime_route_available), TRUE)
  expect_identical(
    unique(direct$runtime_source_spec_status),
    "audited_source_identity_unverified"
  )
  expect_identical(
    unique(direct$runtime_admission_status),
    "not_admitted_source_identity_unverified"
  )
  expect_true(grepl(
    "distinct native gamboost",
    direct$runtime_source_spec_difference[direct$pipeline_label == "GamBoost"],
    fixed = TRUE
  ))
})

test_that("available candidate routes refer only to admitted runtime templates", {
  candidates <- sb_candidate_registry()
  routes <- sb_pipeline_registry()
  available <- candidates[candidates$runtime_route_available, , drop = FALSE]
  candidate_keys <- paste(
    available$runtime_task_type,
    available$runtime_selector,
    available$runtime_learner,
    sep = "::"
  )
  route_keys <- paste(
    routes$task_type,
    routes$selector,
    routes$learner,
    sep = "::"
  )

  expect_gt(nrow(available), 0L)
  expect_identical(all(candidate_keys %in% route_keys), TRUE)
  expect_identical(
    sort(unique(available$runtime_admission_status)),
    sort(c(
      "not_admitted_source_identity_ambiguous",
      "not_admitted_source_identity_unverified",
      "not_admitted_source_to_spec_unverified"
    ))
  )
  expect_identical(
    any(available$execution_status == "executable"),
    FALSE
  )
})

test_that("StepCox source identity ambiguity is explicit without admission", {
  registry <- sb_candidate_registry()
  available <- registry$runtime_route_available
  stepcox <- available & registry$runtime_selector == "stepcox"
  mime_binary <- available &
    registry$runtime_task_type == "binary" &
    registry$model_type == "single" &
    registry$runtime_selector == "none" &
    registry$pipeline_label %in% c("NaiveBayes", "AdaBoost", "LogitBoost")
  survival_glmnet <- available &
    registry$runtime_task_type == "survival" &
    registry$model_type == "single" &
    registry$runtime_selector == "none" &
    registry$runtime_learner == "glmnet_cox"
  binary_linear <- available &
    registry$runtime_task_type == "binary" &
    registry$model_type == "single" &
    registry$runtime_selector == "none" &
    registry$runtime_learner %in% c("logistic", "glmnet_logistic")
  direct_survival_fixed <- available &
    registry$runtime_task_type == "survival" &
    registry$model_type == "single" &
    registry$runtime_selector == "none" &
    registry$runtime_learner %in% c(
      "survivalsvm",
      "coxboost",
      "gbm_cox",
      "superpc",
      "plsrcox",
      "ranger_survival",
      sb_mboost_cox_learner_ids()
    )
  new_binary <- available &
    registry$runtime_task_type == "binary" &
    registry$model_type == "single" &
    registry$runtime_selector == "none" &
    registry$runtime_learner %in% sb_classical_binary_learner_ids()
  new_specialized_binary <- available &
    registry$runtime_task_type == "binary" &
    registry$model_type == "single" &
    registry$runtime_selector == "none" &
    registry$runtime_learner %in%
      sb_specialized_binary_linear_learner_ids()
  new_expansion <- available &
    registry$model_type == "single" &
    registry$runtime_selector == "none" &
    registry$runtime_learner %in% c(
      "xgboost_binary", "xgboost_cox",
      "survreg_weibull", "rpart_survival",
      "rsf_survival", "ctree_survival",
      "cforest_survival"
    )
  former_python <- registry$learner %in% c(
      "CoxTime", "DeepHit", "DeepSurv", "Logistic-Hazard", "PC-Hazard"
    )
  other_available <- available & !stepcox & !mime_binary &
    !survival_glmnet & !binary_linear & !direct_survival_fixed & !new_binary &
    !new_specialized_binary & !new_expansion
  unavailable <- !available

  expect_gt(sum(stepcox), 0L)
  expect_identical(sum(other_available), 0L)
  expect_identical(sum(former_python), 60L)
  expect_true(all(!registry$runtime_route_available[former_python]))
  expect_true(all(
    registry$runtime_route_status[former_python] %in%
      c("component_unavailable", "methodologically_blocked")
  ))
  former_unblocked <- former_python &
    registry$runtime_route_status == "component_unavailable"
  expect_true(all(
    registry$runtime_admission_status[former_unblocked] ==
      "not_admitted_runtime_route_unavailable"
  ))
  expect_true(all(grepl(
    "learner:",
    registry$missing_runtime_components[former_unblocked]
  )))
  expect_identical(sum(new_binary), 8L)
  expect_identical(
    sort(registry$pipeline_label[new_binary]),
    sort(c(
      "DecisionTree", "GBM", "KNN", "LDA", "QDA", "RF", "SVM",
      "svmRadialWeights"
    ))
  )
  expect_identical(sum(new_specialized_binary), 5L)
  expect_identical(
    sort(registry$pipeline_label[new_specialized_binary]),
    sort(c(
      "Stepglm[forward]", "Stepglm[backward]", "Stepglm[both]",
      "glmBoost", "plsRglm"
    ))
  )
  expect_identical(sum(new_expansion), 6L)
  expect_identical(
    unique(registry$runtime_source_spec_status[new_expansion]),
    "not_yet_audited"
  )
  expect_identical(
    unique(registry$runtime_source_spec_status[stepcox]),
    "audited_source_identity_ambiguous"
  )
  expect_identical(
    unique(registry$runtime_source_spec_status[new_binary]),
    "not_yet_audited"
  )
  expect_identical(
    unique(registry$runtime_source_spec_status[new_specialized_binary]),
    "not_yet_audited"
  )
  expect_identical(
    unique(registry$runtime_source_spec_status[unavailable]),
    "not_applicable"
  )
  expect_identical(
    unique(registry$runtime_admission_status[stepcox]),
    "not_admitted_source_identity_ambiguous"
  )
  expect_identical(
    all(nzchar(registry$runtime_source_spec_reference[stepcox])),
    TRUE
  )
  expect_identical(
    all(grepl(
      "Mime1.*Fig.5A",
      registry$runtime_source_spec_reference[stepcox]
    )),
    TRUE
  )
  expect_identical(
    all(grepl(
      "forward.*null",
      registry$runtime_source_spec_difference[stepcox]
    )),
    TRUE
  )
  expect_identical(
    all(registry$adapter_status[stepcox] == "not_verified"),
    TRUE
  )
  expect_identical(
    any(registry$execution_status[stepcox] == "executable"),
    FALSE
  )
})

test_that("direct survival glmnet source identity ambiguity is explicit", {
  registry <- sb_candidate_registry()
  labels <- c(
    "Ridge",
    sprintf("Enet[alpha=%.1f]", seq(0.1, 0.9, by = 0.1)),
    "Lasso"
  )
  direct <- registry[
    registry$runtime_route_available &
      registry$runtime_task_type == "survival" &
      registry$model_type == "single" &
      registry$runtime_selector == "none" &
      registry$runtime_learner == "glmnet_cox", ,
    drop = FALSE
  ]

  expect_identical(nrow(direct), 11L)
  expect_identical(sort(direct$pipeline_label), sort(labels))
  expect_identical(anyDuplicated(direct$pipeline_label), 0L)
  expect_identical(
    unique(direct$evidence_status),
    "historical_fig5a_label_inferred"
  )
  expect_identical(
    unique(direct$runtime_source_spec_status),
    "audited_source_identity_ambiguous"
  )
  expect_identical(
    unique(direct$runtime_admission_status),
    "not_admitted_source_identity_ambiguous"
  )
  expect_identical(
    all(grepl(
      "Mime1.*historical Fig.5A.*historical_fig5a_label_inferred",
      direct$runtime_source_spec_reference
    )),
    TRUE
  )
  expect_identical(
    all(grepl(
      "does not establish which implementation generated",
      direct$runtime_source_spec_difference
    )),
    TRUE
  )

  ridge <- direct[direct$pipeline_label == "Ridge", , drop = FALSE]
  lasso <- direct[direct$pipeline_label == "Lasso", , drop = FALSE]
  enet <- direct[grepl("^Enet", direct$pipeline_label), , drop = FALSE]
  expect_match(ridge$runtime_source_spec_difference, "omits alpha")
  expect_match(ridge$runtime_source_spec_difference, "default alpha=1")
  expect_match(ridge$runtime_source_spec_difference, "cross-validates alpha=0")
  expect_match(lasso$runtime_source_spec_difference, "fixed 10 folds")
  expect_match(lasso$runtime_source_spec_difference, "sample-size-dependent")
  expect_identical(
    all(grepl(
      "fixed 10 folds.*sample-size-dependent",
      enet$runtime_source_spec_difference
    )),
    TRUE
  )
  expect_identical(
    all(grepl(
      paste0(
        "declared inner folds, seed, measure, lambda rule, and ",
        "standardization"
      ),
      direct$runtime_source_spec_difference
    )),
    TRUE
  )
  expect_identical(unique(direct$adapter_status), "not_verified")
  expect_identical(any(direct$execution_status == "executable"), FALSE)
})

test_that("direct binary linear source identity remains unverified", {
  registry <- sb_candidate_registry()
  labels <- c(
    "LogisticRegression",
    "Ridge",
    "Lasso",
    sprintf("Enet[alpha=%.1f]", seq(0.1, 0.9, by = 0.1))
  )
  direct <- registry[
    registry$runtime_route_available &
      registry$runtime_task_type == "binary" &
      registry$model_type == "single" &
      registry$runtime_selector == "none" &
      registry$runtime_learner %in% c("logistic", "glmnet_logistic"), ,
    drop = FALSE
  ]

  expect_identical(nrow(direct), 12L)
  expect_identical(sort(direct$pipeline_label), sort(labels))
  expect_identical(anyDuplicated(direct$pipeline_label), 0L)
  expect_identical(
    unique(direct$evidence_status),
    "union_generated_candidate_endpoint_unverified"
  )
  expect_identical(all(direct$source_doi == ""), TRUE)
  expect_identical(all(direct$source_figure == ""), TRUE)
  expect_identical(
    unique(direct$runtime_source_spec_status),
    "audited_source_identity_unverified"
  )
  expect_identical(
    unique(direct$runtime_admission_status),
    "not_admitted_source_identity_unverified"
  )
  expect_identical(
    all(grepl(
      paste0(
        "Catalog construction reference:.*",
        "union_generated_candidate_endpoint_unverified.*",
        "source_mime_binary=no"
      ),
      direct$runtime_source_spec_reference
    )),
    TRUE
  )
  expect_identical(
    all(grepl(
      "Official Mime1.*neither LogisticRegression nor glmnet",
      direct$runtime_source_spec_reference
    )),
    TRUE
  )
  expect_identical(
    all(grepl(
      "does not establish that this operational implementation generated",
      direct$runtime_source_spec_difference
    )),
    TRUE
  )

  logistic <- direct[direct$runtime_learner == "logistic", , drop = FALSE]
  glmnet <- direct[direct$runtime_learner == "glmnet_logistic", , drop = FALSE]
  expect_identical(nrow(logistic), 1L)
  expect_identical(nrow(glmnet), 11L)
  expect_match(
    logistic$runtime_source_spec_difference,
    "unpenalized binomial-logit"
  )
  expect_match(logistic$runtime_source_spec_difference, "converged full-rank")
  expect_match(logistic$runtime_source_spec_difference, "warning or invalid fit")
  expect_identical(
    all(grepl(
      "label and alpha.*does not specify.*lambda selection",
      glmnet$runtime_source_spec_difference
    )),
    TRUE
  )
  expect_identical(
    all(grepl(
      paste0(
        "declared inner folds, seed, measure, lambda rule, and ",
        "standardization"
      ),
      glmnet$runtime_source_spec_difference
    )),
    TRUE
  )
  expect_identical(
    all(grepl(
      "inside each outer analysis fold.*without fallback",
      glmnet$runtime_source_spec_difference
    )),
    TRUE
  )
  expect_identical(unique(direct$adapter_status), "not_verified")
  expect_identical(any(direct$execution_status == "executable"), FALSE)
})

test_that("remaining direct survival source identities are explicit", {
  registry <- sb_candidate_registry()
  learner_map <- c(
    `survival-SVM` = "survivalsvm",
    CoxBoost = "coxboost",
    GBM = "gbm_cox",
    SuperPC = "superpc",
    plsRcox = "plsrcox",
    Ranger = "ranger_survival"
  )
  direct <- registry[
    registry$runtime_route_available &
      registry$runtime_task_type == "survival" &
      registry$model_type == "single" &
      registry$runtime_selector == "none" &
      registry$runtime_learner %in% unname(learner_map), ,
    drop = FALSE
  ]

  expect_identical(nrow(direct), 6L)
  expect_identical(sort(direct$pipeline_label), sort(names(learner_map)))
  expect_identical(anyDuplicated(direct$pipeline_label), 0L)
  expect_identical(
    unname(learner_map[direct$pipeline_label]),
    direct$runtime_learner
  )

  inferred_labels <- c("survival-SVM", "CoxBoost", "GBM", "SuperPC")
  inferred <- direct[direct$pipeline_label %in% inferred_labels, , drop = FALSE]
  unverified <- direct[
    direct$pipeline_label %in% c("plsRcox", "Ranger"), ,
    drop = FALSE
  ]
  expect_identical(nrow(inferred), 4L)
  expect_identical(
    unique(inferred$evidence_status),
    "historical_fig5a_label_inferred"
  )
  expect_identical(
    unique(inferred$runtime_source_spec_status),
    "audited_source_identity_ambiguous"
  )
  expect_identical(
    unique(inferred$runtime_admission_status),
    "not_admitted_source_identity_ambiguous"
  )
  expect_identical(
    unverified$evidence_status[match(
      c("plsRcox", "Ranger"),
      unverified$pipeline_label
    )],
    c(
      "published_101_117_compatibility_candidate",
      "table1_catalog_candidate_endpoint_unverified"
    )
  )
  expect_identical(
    unique(unverified$runtime_source_spec_status),
    "audited_source_identity_unverified"
  )
  expect_identical(
    unique(unverified$runtime_admission_status),
    "not_admitted_source_identity_unverified"
  )

  by_label <- function(label) {
    direct[direct$pipeline_label == label, , drop = FALSE]
  }
  coxboost <- by_label("CoxBoost")
  gbm <- by_label("GBM")
  superpc <- by_label("SuperPC")
  survival_svm <- by_label("survival-SVM")
  plsrcox <- by_label("plsRcox")
  ranger <- by_label("Ranger")

  expect_match(coxboost$runtime_source_spec_difference, "maxstepno=500")
  expect_match(coxboost$runtime_source_spec_difference, "maxstepno=200")
  expect_match(coxboost$runtime_source_spec_difference, "fixed step count")
  expect_match(gbm$runtime_source_spec_difference, "10,000")
  expect_match(gbm$runtime_source_spec_difference, "2,000")
  expect_match(gbm$runtime_source_spec_difference, "no early stopping")
  expect_match(superpc$runtime_source_spec_difference, "retries indefinitely")
  expect_match(superpc$runtime_source_spec_difference, "validation outcomes")
  expect_match(superpc$runtime_source_spec_difference, "fixed threshold")
  expect_match(survival_svm$runtime_source_spec_difference, "package defaults")
  expect_match(survival_svm$runtime_source_spec_difference, "raw predicted")
  expect_match(survival_svm$runtime_source_spec_difference, "Cox calibration")
  expect_match(plsrcox$runtime_source_spec_reference, "101/117 compatibility")
  expect_match(plsrcox$runtime_source_spec_difference, "cv.plsRcox")
  expect_match(plsrcox$runtime_source_spec_difference, "fixed component count")
  expect_match(ranger$runtime_source_spec_reference, "Table 1")
  expect_match(ranger$runtime_source_spec_difference, "randomForestSRC")
  expect_match(ranger$runtime_source_spec_difference, "not ranger")

  expect_identical(unique(direct$adapter_status), "not_verified")
  expect_identical(any(direct$execution_status == "executable"), FALSE)
})

test_that("Mime1 is only a reference when binary source identity is unverified", {
  registry <- sb_candidate_registry()
  labels <- c("NaiveBayes", "AdaBoost", "LogitBoost")
  direct <- registry[
    registry$runtime_route_available &
      registry$runtime_task_type == "binary" &
      registry$model_type == "single" &
      registry$runtime_selector == "none" &
      registry$pipeline_label %in% labels, ,
    drop = FALSE
  ]

  expect_identical(sort(direct$pipeline_label), sort(labels))
  expect_identical(anyDuplicated(direct$pipeline_label), 0L)
  expect_identical(
    unique(direct$runtime_source_spec_status),
    "audited_source_identity_unverified"
  )
  expect_identical(
    unique(direct$runtime_admission_status),
    "not_admitted_source_identity_unverified"
  )
  expect_identical(
    all(grepl(
      "Comparison reference only: Mime1.*source identity remains unverified",
      direct$runtime_source_spec_reference
    )),
    TRUE
  )

  naive_bayes <- direct[direct$pipeline_label == "NaiveBayes", , drop = FALSE]
  adaboost <- direct[direct$pipeline_label == "AdaBoost", , drop = FALSE]
  logitboost <- direct[direct$pipeline_label == "LogitBoost", , drop = FALSE]
  expect_match(naive_bayes$runtime_source_spec_difference, "usekernel=TRUE")
  expect_match(naive_bayes$runtime_source_spec_difference, "Gaussian")
  expect_match(adaboost$runtime_source_spec_difference, "Real AdaBoost")
  expect_match(adaboost$runtime_source_spec_difference, "fixed AdaBoost.M1")
  expect_match(logitboost$runtime_source_spec_difference, "11 through 101")
  expect_match(logitboost$runtime_source_spec_difference, "fixed iteration")
  expect_identical(
    all(grepl(
      "does not establish that Mime1 generated this candidate",
      direct$runtime_source_spec_difference
    )),
    TRUE
  )
  expect_identical(unique(direct$adapter_status), "not_verified")
  expect_identical(
    any(direct$execution_status == "executable"),
    FALSE
  )
})

test_that("candidate registry filters one task or route status", {
  all_candidates <- sb_candidate_registry()
  binary <- sb_candidate_registry(task_type = "binary")
  available <- sb_candidate_registry(route_status = "route_available")
  bad_task <- capture_sb_error(sb_candidate_registry(
    task_type = "diagnostic"
  ))
  bad_status <- capture_sb_error(sb_candidate_registry(
    route_status = "executable"
  ))

  expect_identical(nrow(binary), 584L)
  expect_identical(unique(binary$runtime_task_type), "binary")
  expect_identical(
    nrow(available),
    sum(all_candidates$runtime_route_status == "route_available")
  )
  expect_identical(unique(available$runtime_route_status), "route_available")
  expect_s3_class(bad_task, "stormbreaker_bad_candidate_registry_task")
  expect_s3_class(
    bad_status,
    "stormbreaker_bad_candidate_registry_route_status"
  )
})

# ---- survival-menu-pipelines ----

test_that("source-backed survival menu preserves registry pairs and MIME settings", {
  features <- paste0("gene_", seq_len(60L))
  menu <- sb_survival_menu_pipelines(
    selector_features = features,
    feature_set_id = "candidate60_test"
  )
  registry <- sb_pipeline_registry(task_type = "survival")
  expect_setequal(
    paste(menu$selector, menu$learner),
    paste(registry$selector, registry$learner)
  )
  expect_identical(anyDuplicated(menu$config_id), 0L)
  expect_identical(anyDuplicated(vapply(menu$pipeline, `[[`, "", "id")), 0L)
  reference <- sb_mime117_pipelines(
    inner_seed = 1L, selector_seed = 1L, learner_seed = 1L
  )
  mapping <- attr(menu, "mime117_mapping")
  expect_identical(nrow(mapping), 117L)
  expect_identical(mapping$official_model_label, reference$official_model_label)
  indices <- match(mapping$comparison_config_id, menu$config_id)
  expect_identical(anyNA(indices), FALSE)
  expect_identical(unname(menu$pipeline[indices]), unname(reference$pipeline))
  expect_identical(sum(menu$comparison_role == "MIME_structure_matched"), 117L)
  expect_identical(attr(menu, "registry_pairs"), nrow(registry))
  for (p in menu$pipeline) {
    if (p$selector == "fixed_features") {
      expect_identical(p$selector_params$features, features)
      expect_identical(p$selector_params$feature_set_id, "candidate60_test")
    }
    if (p$selector %in% c("lasso_nonzero", "coxboost_nonzero",
      "lasso_then_stepcox", "coxboost_then_stepcox", "enet_cox_nonzero")) {
      expect_identical(p$selector_params$inner_folds, 10L)
      expect_identical(p$selector_params$seed, 1L)
    }
    if (p$learner == "glmnet_cox") {
      expect_identical(p$learner_params$inner_folds, 10L)
      expect_identical(p$learner_params$inner_seed, 1L)
      expect_identical(p$learner_params$selection_rule, "lambda.min")
      expect_identical(p$learner_params$standardize, TRUE)
    }
    if (p$learner %in% c("coxboost_cv", "gbm_cox_cv", "superpc_cv")) {
      expect_identical(p$learner_params$inner_folds, 10L)
      expect_identical(p$learner_params$seed, 1L)
    }
    if (p$learner == "plsrcox_cv") {
      expect_identical(p$learner_params$inner_folds, 5L)
      expect_identical(p$learner_params$maximum_components, 10L)
    }
    if (p$learner == "gbm_cox_cv") {
      expect_identical(p$learner_params$maximum_trees, 10000L)
    }
    if (p$selector == "boruta_survival_roughfix") {
      expect_null(p$selector_params)
    }
  }
})

test_that("survival menu validates caller-declared features without data access", {
  invalid <- tryCatch(sb_survival_menu_pipelines(
    selector_features = c("g", "g"), feature_set_id = "duplicate"
  ), error = identity)
  expect_s3_class(invalid, "stormbreaker_bad_fixed_features")
  invalid <- tryCatch(sb_survival_menu_pipelines(
    selector_features = c("a", "b"), feature_set_id = "valid",
    selector_n_features = 0L
  ), error = identity)
  expect_s3_class(invalid, "error")
})
})

# ---- mime117-pipelines ----
local({
test_that("MIME117 generator returns 117 real admitted pipelines", {
  suite <- sb_mime117_pipelines(
    inner_folds = 3L,
    inner_seed = 8701L,
    selector_inner_folds = 3L,
    selector_seed = 8702L,
    rsf_nodesize = 5L
  )

  expect_s3_class(suite, "data.frame")
  expect_identical(nrow(suite), 117L)
  expect_identical(suite$official_order, seq_len(117L))
  expect_identical(anyDuplicated(suite$config_id), 0L)
  expect_identical(anyDuplicated(suite$official_model_label), 0L)
  expect_identical(
    anyDuplicated(vapply(suite$pipeline, `[[`, character(1), "id")),
    0L
  )
  expect_identical(
    as.integer(table(factor(
      suite$block,
      levels = c(
        "RSF", "Enet", "StepCox", "CoxBoost", "plsRcox", "SuperPC",
        "GBM", "survival-SVM", "Ridge", "Lasso"
      )
    ))),
    c(20L, 9L, 54L, 19L, 1L, 1L, 1L, 1L, 1L, 10L)
  )
  expect_identical(
    as.integer(table(suite$pipeline_type)),
    c(20L, 97L)
  )
  expect_true(all(vapply(
    suite$pipeline,
    inherits,
    logical(1),
    what = "stormbreaker_pipeline"
  )))

  lint_status <- vapply(suite$pipeline, function(pipeline) {
    sb_lint_pipeline(pipeline)$status
  }, character(1))
  expect_identical(unique(lint_status), "admitted")
})

test_that("MIME117 generator preserves variants and stage roles", {
  suite <- sb_mime117_pipelines(
    inner_folds = 3L,
    inner_seed = 8711L,
    selector_inner_folds = 3L,
    selector_seed = 8712L,
    rsf_nodesize = 5L
  )

  enet <- suite[grepl("Enet", suite$official_model_label, fixed = TRUE), ]
  stepcox <- suite[grepl("StepCox", suite$official_model_label, fixed = TRUE), ]
  expect_setequal(
    unique(enet$enet_alpha),
    seq(0.1, 0.9, by = 0.1)
  )
  expect_setequal(
    unique(stats::na.omit(stepcox$stepcox_direction)),
    c("both", "backward", "forward")
  )
  expect_identical(suite$official_model_label[[1L]], "RSF")
  expect_identical(suite$pipeline[[1L]]$selector, "none")
  expect_identical(suite$pipeline[[1L]]$learner, "rsf_survival")

  rsf_coxboost <- match("RSF + CoxBoost", suite$official_model_label)
  expect_identical(
    suite$pipeline[[rsf_coxboost]]$selector,
    "rsf_minimal_depth"
  )
  expect_identical(suite$pipeline[[rsf_coxboost]]$learner, "coxboost_cv")
  coxboost <- suite$pipeline[[match("CoxBoost", suite$official_model_label)]]
  plsrcox <- suite$pipeline[[match("plsRcox", suite$official_model_label)]]
  superpc <- suite$pipeline[[match("SuperPC", suite$official_model_label)]]
  gbm <- suite$pipeline[[match("GBM", suite$official_model_label)]]
  expect_identical(coxboost$learner, "coxboost_cv")
  expect_identical(coxboost$learner_params$start_penalty, 500)
  expect_identical(coxboost$learner_params$maximum_steps, 500L)
  expect_identical(plsrcox$learner, "plsrcox_cv")
  expect_identical(plsrcox$learner_params$maximum_components, 10L)
  expect_identical(plsrcox$learner_params$inner_folds, 5L)
  expect_identical(superpc$learner, "superpc_cv")
  expect_identical(superpc$learner_params$threshold_count, 20L)
  expect_identical(superpc$learner_params$tuning_components, 3L)
  expect_identical(superpc$learner_params$prediction_components, 1L)
  expect_identical(gbm$learner, "gbm_cox_cv")
  expect_identical(gbm$learner_params$maximum_trees, 10000L)
  expect_identical(gbm$learner_params$interaction_depth, 3L)
  expect_identical(gbm$learner_params$shrinkage, 0.001)

  lasso_stepcox <- match(
    "Lasso + StepCox[forward]",
    suite$official_model_label
  )
  expect_identical(
    suite$pipeline[[lasso_stepcox]]$selector,
    "lasso_then_stepcox"
  )
  expect_identical(suite$pipeline[[lasso_stepcox]]$learner, "coxph")
  expect_identical(
    suite$pipeline[[lasso_stepcox]]$selector_params$direction,
    "forward"
  )
  expect_identical(
    unique(suite$operational_fidelity),
    "fold_local_mime_training_reproduction_not_bug_for_bug"
  )
})
})

# ---- pipeline-registry ----
local({
registry_test_arguments <- function(row) {
  arguments <- list(
    selector = row$selector,
    learner = row$learner
  )
  if (identical(row$selector_parameter, "selector_n_features")) {
    arguments$selector_n_features <- 2L
  }
  if (identical(row$selector_parameter, "stepcox_direction")) {
    arguments$stepcox_direction <- "forward"
  }
  if (row$selector_parameter %in% c(
    "selector_rsf_nodesize_and_seed",
    "selector_rsf_nodesize_seed_and_stepcox_direction"
  )) {
    arguments$selector_rsf_nodesize <- 5L
    arguments$selector_seed <- 101L
  }
  if (row$selector_parameter %in% c(
    "selector_inner_folds_and_seed",
    "selector_inner_folds_seed_and_stepcox_direction"
  )) {
    arguments$selector_inner_folds <- 3L
    arguments$selector_seed <- 101L
  }
  if (row$selector_parameter %in% c(
    "selector_rsf_nodesize_seed_and_stepcox_direction",
    "selector_inner_folds_seed_and_stepcox_direction"
  )) {
    arguments$stepcox_direction <- "forward"
  }
  if (identical(
    row$selector_parameter,
    "selector_features_and_feature_set_id"
  )) {
    arguments$selector_features <- c("gene_a", "gene_b")
    arguments$feature_set_id <- "registry_pair"
  }
  if (identical(
    row$tuning_semantics,
    "lambda_selected_by_inner_cv"
  )) {
    arguments$inner_folds <- 3L
    arguments$inner_seed <- 101L
  }
  if (identical(row$tuning_semantics, "analysis_fold_internal_cv")) {
    arguments$learner_tuning_folds <- 3L
    arguments$learner_seed <- 101L
    if (identical(row$learner, "gbm_cox_cv")) {
      arguments$gbm_cv_max_trees <- 60L
    }
  }
  arguments
}

removed_python_learner_ids <- function() {
  c(
    "coxtime_fixed", "deephit_fixed", "deepsurv_fixed",
    "logistic_hazard_fixed", "pc_hazard_fixed"
  )
}

test_that("the public registry lists only admitted executable combinations", {
  registry <- sb_pipeline_registry()

  expect_s3_class(registry, "data.frame")
  expect_identical(
    names(registry),
    c(
      "source", "task_type", "selector", "learner",
      "admission_status", "selector_parameter", "tuning_semantics",
      "optional_dependency", "parameter_semantics"
    )
  )
  expect_setequal(
    unique(registry$source),
    c(
      "stormbreaker_builtin", "stormbreaker_abess_runtime",
      "stormbreaker_forest_selector_runtime",
      "stormbreaker_high_dimensional_runtime",
      "stormbreaker_binary_model_selector_runtime",
      "stormbreaker_survival_boosting_selector_runtime",
      "stormbreaker_pamr_selector_runtime",
      "stormbreaker_survival_weight_selector_runtime",
      "stormbreaker_mime117_runtime"
    )
  )
  expect_identical(unique(registry$admission_status), "admitted")
  expect_false(anyDuplicated(registry[c("selector", "learner")]) > 0L)
  expect_identical(
    nrow(sb_pipeline_registry(level = "learner")),
    176L
  )
  expect_false(any(registry$learner %in% removed_python_learner_ids()))
  for (learner in removed_python_learner_ids()) {
    expect_error(
      sb_pipeline(learner = learner),
      class = "stormbreaker_unsupported_learner"
    )
  }

  observed <- paste(
    registry$task_type,
    registry$selector,
    registry$learner,
    sep = "::"
  )
  survival_full_learners <- c(
    "coxph", "glmnet_cox", "catboost_cox",
    "catboost_aft_survival", "gbm_cox",
    sb_mboost_cox_learner_ids(),
    sb_mboost_aft_learner_ids(),
    sb_mboost_ranking_learner_ids(),
    sb_ranger_survival_learner_ids(), sb_rlt_survival_learner_ids(),
    sb_coxphf_learner_ids(), sb_bagged_knn_survival_learner_ids(),
    sb_ahaz_survival_learner_ids(), sb_bujar_linear_learner_ids(), sb_akritas_learner_ids(),
    "coxboost", "plsrcox",
    "superpc", sb_survivalsvm_learner_ids(),
    sb_xgboost_cox_learner_ids(),
    sb_xgboost_aft_learner_ids(),
    sb_survreg_learner_ids(),
    sb_flexsurv_parametric_learner_ids(),
    sb_flexsurvspline_learner_ids(),
    sb_mtlr_survival_learner_ids(),
    sb_pamr_survival_learner_ids(), "rpart_survival",
    "rsf_survival",
    "ctree_survival", "cforest_survival", "bart_survival",
    "aorsf_survival", sb_timereg_survival_learner_ids(),
    sb_ncvreg_survival_learner_ids(),
    "abess_cox_gic",
    sb_analysis_tuned_survival_learner_ids(),
    "bart_aft_survival"
  )
  survival_restricted_learners <- c(
    "pre_survival_rule_ensemble",
    "penaft_gehan_elastic_net"
  )
  survival_learners <- c(
    survival_full_learners,
    survival_restricted_learners
  )
  binary_full_learners <- c(
    "logistic", "glmnet_logistic", "gaussian_naive_bayes",
    sb_omics_specialist_binary_learner_ids(),
    sb_lightgbm_binary_learner_ids(),
    "catboost_binary", sb_ranger_binary_learner_ids(),
    sb_rlt_binary_learner_ids(), sb_wsrf_binary_learner_ids(),
    sb_random_glm_binary_learner_ids(),
    "adaboost_binary", "logitboost_binary",
    sb_classical_binary_learner_ids(),
    sb_specialized_binary_linear_learner_ids(),
    sb_xgboost_binary_learner_ids(), sb_ncvreg_binary_learner_ids(),
    "abess_binary_gic",
    sb_mclust_da_learner_ids(), sb_mixomics_da_learner_ids(),
    sb_hdda_binary_learner_ids()
  )
  binary_restricted_learners <- c(
    "pre_binary_rule_ensemble",
    sb_slope_binary_learner_ids(),
    "cancerclass_binary"
  )
  binary_learners <- c(binary_full_learners, binary_restricted_learners)
  paired_keys <- function(task, selectors, learners) {
    unlist(lapply(learners, function(learner) {
      paste(task, selectors, learner, sep = "::")
    }), use.names = FALSE)
  }
  expected <- c(
    paired_keys(
      "survival",
      c("none", "univ_cox_top_k"),
      survival_full_learners
    ),
    paste(
      "survival",
      "none",
      survival_restricted_learners,
      sep = "::"
    ),
    paired_keys(
      "binary",
      c("none", "univ_binary_auc_top_k"),
      binary_full_learners
    ),
    paste("binary", "none", binary_restricted_learners, sep = "::"),
    paste(
      "survival",
      "stepcox",
      survival_full_learners,
      sep = "::"
    ),
    paste(
      "survival",
      "rsf_minimal_depth",
      c(
        "coxboost_cv", "glmnet_cox", "gbm_cox_cv",
        "plsrcox_cv", "superpc_cv", sb_survivalsvm_learner_ids()
      ),
      sep = "::"
    ),
    paste(
      "survival",
      "coxboost_nonzero",
      c(
        "glmnet_cox", "gbm_cox_cv", "plsrcox_cv",
        "superpc_cv", sb_survivalsvm_learner_ids()
      ),
      sep = "::"
    ),
    paste(
      "survival",
      "lasso_nonzero",
      c(
        "coxboost_cv", "gbm_cox_cv", "plsrcox_cv",
        "rsf_survival", "superpc_cv", sb_survivalsvm_learner_ids()
      ),
      sep = "::"
    ),
    paste(
      "survival",
      c(
        "rsf_then_stepcox", "coxboost_then_stepcox",
        "lasso_then_stepcox"
      ),
      "coxph",
      sep = "::"
    ),
    paste(
      "binary",
      "abess_binary_gic_nonzero",
      setdiff(binary_full_learners, "abess_binary_gic"),
      sep = "::"
    ),
    paste(
      "survival",
      "abess_cox_gic_nonzero",
      setdiff(
        survival_full_learners,
        c("abess_cox_gic", sb_analysis_tuned_survival_learner_ids())
      ),
      sep = "::"
    ),
    paired_keys(
      "binary",
      c(
        "boruta_binary_roughfix", "varselrf_binary_oob_1se",
        "lasso_binary_nonzero", "stepglm_binary_aic",
        "glmboost_binary_nonzero", "rf_binary_importance",
        "enet_binary_nonzero", "gbm_binary_importance", "svm_rfe_binary", "gamboost_binary_nonzero", "xgboost_binary_importance", "plsrglm_binary_weights",
        setdiff(
          sb_high_dimensional_binary_selector_ids(),
          "stabsel_binary_lasso"
        )
      ),
      binary_full_learners
    ),
    paste(
      "binary",
      "stabsel_binary_lasso",
      setdiff(binary_full_learners, "glmnet_logistic"),
      sep = "::"
    ),
    paste(
      "survival",
      "boruta_survival_roughfix",
      setdiff(
        survival_full_learners,
        sb_analysis_tuned_survival_learner_ids()
      ),
      sep = "::"
    ),
    paste(
      "survival",
      "lpc_cox_nlogn",
      setdiff(
        survival_full_learners,
        sb_analysis_tuned_survival_learner_ids()
      ),
      sep = "::"
    ),
    paste(
      "survival",
      "stabsel_cox_lasso",
      setdiff(
        survival_full_learners,
        c(sb_analysis_tuned_survival_learner_ids(), "glmnet_cox")
      ),
      sep = "::"
    ),
    paste("survival", "fixed_features", survival_learners, sep = "::"),
    paste("binary", "fixed_features", binary_learners, sep = "::")
  )
  expected <- c(expected,
    paired_keys("binary", "pamr_binary_nonzero", binary_full_learners),
    paired_keys("survival", "pamr_survival_nonzero",
      setdiff(survival_full_learners, sb_analysis_tuned_survival_learner_ids())),
    paired_keys("survival", c("lasso_nonzero", "enet_cox_nonzero", "coxboost_nonzero", "rsf_minimal_depth"),
      setdiff(survival_full_learners, sb_analysis_tuned_survival_learner_ids())),
    paired_keys("survival", c("ridge_cox_weights", "plsrcox_weights", "svm_survival_weights"),
      setdiff(survival_full_learners, sb_analysis_tuned_survival_learner_ids())),
    paired_keys("survival", c("gbm_cox_importance", "glmboost_cox_nonzero", "superpc_survival_scores", "gamboost_cox_nonzero", "xgboost_cox_importance"),
      setdiff(survival_full_learners, sb_analysis_tuned_survival_learner_ids())),
    paired_keys("survival", "enet_cox_nonzero",
      unique(c("coxph", "glmnet_cox", "ranger_survival", "coxboost_cv",
        "gbm_cox_cv", "plsrcox_cv", "rsf_survival", "superpc_cv", sb_survivalsvm_learner_ids()))))
  expect_setequal(observed, expected)
})

test_that("every registry row resolves through the existing runtime", {
  registry <- sb_pipeline_registry()
  ids <- character(nrow(registry))

  for (index in seq_len(nrow(registry))) {
    row <- registry[index, , drop = FALSE]
    pipeline <- do.call(sb_pipeline, registry_test_arguments(row))
    runtime <- sb_resolve_pipeline_runtime(pipeline)

    expect_identical(pipeline$task_type, row$task_type)
    expect_identical(pipeline$selector, row$selector)
    expect_identical(pipeline$learner, row$learner)
    expect_identical(
      sb_validate_pipeline_contract(pipeline, row$task_type),
      pipeline
    )
    expect_identical(
      sb_validate_pipeline_runtime(runtime, pipeline),
      runtime
    )
    ids[[index]] <- pipeline$id
  }

  expect_identical(anyDuplicated(ids), 0L)
})

test_that("the registry can summarize one row per learner", {
  routes <- sb_pipeline_registry()
  learners <- sb_pipeline_registry(level = "learner")
  route_counts <- table(routes$learner)
  expected_selectors <- vapply(
    learners$learner,
    function(learner) {
      paste(routes$selector[routes$learner == learner], collapse = ", ")
    },
    character(1)
  )

  expect_identical(
    names(learners),
    c(
      "task_type", "learner", "pipeline_count", "selectors",
      "tuning_semantics", "optional_dependency", "parameter_semantics"
    )
  )
  expect_setequal(learners$learner, unique(routes$learner))
  expect_identical(
    learners$pipeline_count,
    unname(as.integer(route_counts[learners$learner]))
  )
  expect_identical(learners$selectors, unname(expected_selectors))

  binary <- sb_pipeline_registry(task_type = "binary", level = "learner")
  expect_identical(unique(binary$task_type), "binary")
  expect_setequal(
    binary$learner,
    unique(routes$learner[routes$task_type == "binary"])
  )
})

test_that("learner summaries can report installed R dependencies", {
  declared <- sb_pipeline_registry(level = "learner")
  checked <- sb_pipeline_registry(
    level = "learner",
    check_installed = TRUE
  )

  expect_identical(checked[names(declared)], declared)
  expect_identical(
    tail(names(checked), 4L),
    c(
      "r_packages", "python_modules", "missing_r_packages",
      "dependency_status"
    )
  )
  expect_identical(
    checked$dependency_status[checked$learner == "coxph"],
    "no_optional_dependency"
  )
  expect_identical(
    checked$dependency_status[
      checked$learner == "survreg_weibull"
    ],
    "r_packages_installed"
  )
  expect_identical(
    any(checked$learner %in% removed_python_learner_ids()),
    FALSE
  )
  expect_true(
    all(is.na(checked$python_modules) | checked$python_modules == "")
  )
  expect_identical(
    any(checked$dependency_status == "python_not_checked"),
    FALSE
  )
  expect_identical(
    setdiff(
      unique(checked$dependency_status),
      c(
        "no_optional_dependency", "r_packages_installed",
        "missing_r_packages", "python_not_checked"
      )
    ),
    character()
  )
})

test_that("fixed nonlinear learners are not printed as elastic nets", {
  binary <- utils::capture.output(print(sb_pipeline(
    learner = "bart_binary"
  )))
  survival <- utils::capture.output(print(sb_pipeline(
    learner = "bart_survival"
  )))

  expect_true(any(grepl("fixed before resampling", binary, fixed = TRUE)))
  expect_true(any(grepl("fixed before resampling", survival, fixed = TRUE)))
  expect_false(any(grepl("Elastic Net", c(binary, survival), fixed = TRUE)))
})

test_that("the linter admits a raw specification or constructed pipeline", {
  raw_report <- sb_lint_pipeline(list(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 2L
  ))
  object_report <- sb_lint_pipeline(
    sb_pipeline(learner = "logistic"),
    task_type = "binary"
  )

  expect_identical(
    names(raw_report),
    c(
      "status", "pipeline_id", "task_type", "selector", "learner",
      "issue_code", "message"
    )
  )
  expect_identical(raw_report$status, "admitted")
  expect_identical(
    raw_report$pipeline_id,
    "univ_binary_auc_top_k__logistic__k2"
  )
  expect_identical(raw_report$task_type, "binary")
  expect_true(is.na(raw_report$issue_code))
  expect_true(is.na(raw_report$message))

  expect_identical(object_report$status, "admitted")
  expect_identical(object_report$pipeline_id, "none__logistic")
})

test_that("the linter reports expected specification and contract errors", {
  missing_inner <- sb_lint_pipeline(list(learner = "glmnet_cox"))
  mismatch <- sb_lint_pipeline(
    list(learner = "logistic"),
    task_type = "survival"
  )
  unknown_argument <- sb_lint_pipeline(list(
    learner = "logistic",
    mystery = 1
  ))
  bad_container <- sb_lint_pipeline("none__coxph")
  bad_learner <- sb_lint_pipeline(list(learner = 123))
  bad_selector <- sb_lint_pipeline(list(
    selector = c("none", "bad")
  ))
  mutated <- sb_pipeline()
  mutated$risk_formula <- "changed"
  mutation <- sb_lint_pipeline(mutated)

  expect_identical(missing_inner$status, "rejected")
  expect_identical(
    missing_inner$issue_code,
    "stormbreaker_missing_inner_parameter"
  )
  expect_identical(
    mismatch$issue_code,
    "stormbreaker_pipeline_task_mismatch"
  )
  expect_identical(
    unknown_argument$issue_code,
    "stormbreaker_unknown_pipeline_argument"
  )
  expect_identical(
    bad_container$issue_code,
    "stormbreaker_bad_pipeline_spec"
  )
  expect_true(is.na(bad_container$selector))
  expect_true(is.na(bad_container$learner))
  expect_identical(bad_learner$selector, "none")
  expect_true(is.na(bad_learner$learner))
  expect_true(is.na(bad_selector$selector))
  expect_identical(bad_selector$learner, "coxph")
  expect_identical(
    mutation$issue_code,
    "stormbreaker_unresolved_pipeline"
  )
  expect_true(all(nzchar(c(
    missing_inner$message,
    mismatch$message,
    unknown_argument$message,
    bad_container$message,
    mutation$message
  ))))
})

test_that("runtime registry and adapter drift remains an internal error", {
  local_mocked_bindings(
    sb_resolve_pipeline_runtime = function(pipeline) {
      sb_abort(
        "injected registry drift",
        "stormbreaker_unresolved_learner_adapter"
      )
    }
  )

  error <- tryCatch(
    sb_lint_pipeline(list(learner = "logistic")),
    error = function(error) error
  )

  expect_s3_class(error, "stormbreaker_unresolved_learner_adapter")
  expect_identical(conditionMessage(error), "injected registry drift")
})

test_that("unexpected contract stormbreaker errors remain visible", {
  local_mocked_bindings(
    sb_validate_pipeline_contract = function(pipeline, task_type) {
      sb_abort(
        "injected internal contract failure",
        "stormbreaker_internal_contract_failure"
      )
    }
  )

  error <- tryCatch(
    sb_lint_pipeline(list(learner = "logistic")),
    error = function(error) error
  )

  expect_s3_class(error, "stormbreaker_internal_contract_failure")
  expect_identical(
    conditionMessage(error),
    "injected internal contract failure"
  )
})

test_that("unexpected constructor stormbreaker errors remain visible", {
  pipeline <- sb_pipeline(learner = "logistic")
  local_mocked_bindings(
    sb_pipeline = function(...) {
      sb_abort(
        "injected internal constructor failure",
        "stormbreaker_internal_constructor_failure"
      )
    }
  )

  error <- tryCatch(
    sb_lint_pipeline(pipeline),
    error = function(error) error
  )

  expect_s3_class(error, "stormbreaker_internal_constructor_failure")
  expect_identical(
    conditionMessage(error),
    "injected internal constructor failure"
  )
})

test_that("mutated pipeline reports use exact fields without defaults", {
  missing <- sb_pipeline(learner = "logistic")
  missing[["learner"]] <- NULL
  missing_report <- sb_lint_pipeline(missing)

  partial <- sb_pipeline(learner = "logistic")
  partial[["learner"]] <- NULL
  partial[["learners"]] <- "spoofed"
  partial_report <- sb_lint_pipeline(partial)

  expect_identical(missing_report$status, "rejected")
  expect_identical(
    missing_report$issue_code,
    "stormbreaker_unresolved_pipeline"
  )
  expect_identical(missing_report$selector, "none")
  expect_true(is.na(missing_report$learner))
  expect_identical(partial_report$status, "rejected")
  expect_identical(
    partial_report$issue_code,
    "stormbreaker_unresolved_pipeline"
  )
  expect_identical(partial_report$selector, "none")
  expect_true(is.na(partial_report$learner))
})

test_that("admission does not claim that optional packages are installed", {
  local_mocked_bindings(sb_lightgbm_available = function() FALSE)

  report <- sb_lint_pipeline(list(learner = "lightgbm_binary"))
  registry_row <- sb_pipeline_registry()
  registry_row <- registry_row[
    registry_row$learner == "lightgbm_binary" &
      registry_row$selector == "none", ,
    drop = FALSE
  ]

  expect_identical(report$status, "admitted")
  expect_identical(registry_row$optional_dependency, "lightgbm")
  expect_identical(registry_row$admission_status, "admitted")
})

test_that("the linter does not disguise unexpected internal errors", {
  local_mocked_bindings(
    sb_resolve_pipeline_runtime = function(pipeline) {
      stop("unexpected internal failure")
    }
  )

  error <- tryCatch(
    sb_lint_pipeline(list(learner = "logistic")),
    error = function(error) error
  )

  expect_s3_class(error, "simpleError")
  expect_false(inherits(error, "stormbreaker_error"))
  expect_identical(conditionMessage(error), "unexpected internal failure")

  pipeline <- sb_pipeline()
  local_mocked_bindings(
    sb_pipeline = function(...) {
      stop("unexpected reconstruction failure")
    }
  )
  reconstruction_error <- tryCatch(
    sb_lint_pipeline(pipeline),
    error = function(error) error
  )

  expect_s3_class(reconstruction_error, "simpleError")
  expect_false(inherits(reconstruction_error, "stormbreaker_error"))
  expect_identical(
    conditionMessage(reconstruction_error),
    "unexpected reconstruction failure"
  )
})
})

# ---- pipeline-specification ----
local({
# ---- pipeline-defensive-restraint ----

test_that("equivalent parameter representations remain executable", {
  p <- sb_pipeline(learner = "nnet_binary")
  canonical <- p
  p$learner_params$maxit <- as.numeric(p$learner_params$maxit)
  p$learner_params <- rev(p$learner_params)
  expect_identical(sb_validate_pipeline_contract(p, "binary"), canonical)
  chain <- sb_new_chain_pipeline(p, list(p))
  chain$learner_params <- rev(chain$learner_params)
  expect_identical(sb_validate_pipeline_contract(chain, "binary"),
    sb_new_chain_pipeline(canonical, list(canonical)))
  p$learner_params$maxit <- 200
  canonical$learner_params$maxit <- 200L
  expect_identical(sb_validate_pipeline_contract(p, "binary"), canonical)
})

test_that("equivalent parameters reach the native forest with usable types", {
  skip_if_not_installed("ranger")
  canonical <- sb_pipeline(learner = "ranger_binary", ranger_num_trees = 5L)
  p <- canonical
  p$learner_params$num_trees <- 5
  p$learner_params <- rev(p$learner_params)
  validated <- sb_validate_pipeline_contract(p, "binary")
  x <- data.frame(a = sin(1:30), b = cos(1:30))
  y <- rep(0:1, 15L)
  fitted <- sb_fit_ranger_binary(x, y, validated$learner_params)
  reference <- sb_fit_ranger_binary(x, y, canonical$learner_params)
  expect_identical(sb_predict_ranger_binary(fitted, x),
    sb_predict_ranger_binary(reference, x))
  chain <- sb_validate_pipeline_contract(sb_new_chain_pipeline(p, list(p)), "binary")
  expect_identical(chain$learner_params, canonical$learner_params)
  expect_identical(chain$selector_params$stages[[1L]]$learner_params,
    canonical$learner_params)
})

test_that("parameter equivalence does not hide scientific changes", {
  expect_identical(sb_pipeline_parameters_equal(
    list(a = 1L, b = list(x = 2L, y = 3L)),
    list(b = list(y = 3, x = 2), a = 1)), TRUE)
  expect_identical(sb_pipeline_parameters_equal(1, 1 + 1e-12), FALSE)
  expect_identical(sb_pipeline_parameters_equal(c(3L, 2L), c(2, 3)), FALSE)
  expect_identical(sb_pipeline_parameters_equal(list(a = 1L), list(b = 1)), FALSE)
  expect_identical(sb_pipeline_parameters_equal(list(1L, 2L), list(2, 1)), FALSE)
  expect_identical(sb_pipeline_parameters_equal(c(x = 1L), c(y = 1)), FALSE)
  expect_identical(sb_pipeline_parameters_equal(TRUE, 1L), FALSE)
  p <- sb_pipeline(learner = "nnet_binary")
  p$learner_params$decay <- p$learner_params$decay + 1e-12
  expect_s3_class(tryCatch(sb_validate_pipeline_contract(p, "binary"),
    error = identity), "stormbreaker_unresolved_pipeline")
})

test_that("chain annotations do not change the executable specification", {
  base <- sb_pipeline(learner = "logistic")
  stage <- sb_pipeline(learner = "logistic")
  chain <- sb_new_chain_pipeline(base, list(stage))
  chain$note <- "Clinician annotation"
  chain$chain_base$note <- "Development cohort"
  chain$selector_params$stages[[1L]]$note <- "First stage"
  expect_identical(sb_validate_pipeline_contract(chain, "binary"), chain)

  chain$positive_class <- 0L
  expect_identical(conditionMessage(tryCatch(
    sb_validate_pipeline_contract(chain, "binary"), error = identity)),
    "Selector chain does not match its executable specification.")
})

test_that("chain validation does not reconstruct an already built chain", {
  base <- sb_pipeline(learner = "logistic")
  chain <- sb_new_chain_pipeline(base, list(base))
  local_mocked_bindings(sb_new_chain_pipeline = function(...) {
    stop("constructor should not run during validation")
  })
  expect_identical(sb_validate_pipeline_contract(chain, "binary"), chain)
})

test_that("invalid chain stages and missing pipeline class still fail", {
  base <- sb_pipeline(learner = "logistic")
  chain <- sb_new_chain_pipeline(base, list(base))
  expect_s3_class(tryCatch(sb_validate_pipeline_contract(unclass(chain), "binary"),
    error = identity), "stormbreaker_bad_pipeline")
  empty <- chain
  empty$selector_params$stages <- list()
  expect_identical(conditionMessage(tryCatch(
    sb_validate_pipeline_contract(empty, "binary"), error = identity)),
    "A selector chain requires stages.")
  chain$selector_params$stages[[1L]] <- sb_pipeline()
  expect_s3_class(tryCatch(sb_validate_pipeline_contract(chain, "binary"),
    error = identity), "stormbreaker_pipeline_task_mismatch")
})

test_that("ordinary C50 warnings remain warnings and preserve the fit", {
  skip_if_not_installed("C50")
  native <- structure(list(size = 1L), class = "C5.0")
  local_mocked_bindings(C5.0 = function(...) {
    warning("informational native warning")
    native
  }, .package = "C50")
  seen <- character()
  result <- withCallingHandlers(sb_fit_c50_binary(
    data.frame(x = 1:4), factor(c(0, 0, 1, 1)), "c50_binary",
    sb_classical_binary_params("c50_binary")), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(result, native)
  expect_identical(seen, "informational native warning")
})

test_that("robust logistic still rejects unconverged native fits", {
  skip_if_not_installed("robustbase")
  local_mocked_bindings(glmrob = function(...) {
    list(converged = FALSE, coefficients = c(0, 1), iter = 1L)
  }, .package = "robustbase")
  result <- tryCatch(sb_fit_robust_logistic_mqle(
    data.frame(x = 1:4), factor(c(0, 0, 1, 1)),
    "robust_logistic_mqle_binary_fixed",
    sb_classical_binary_params("robust_logistic_mqle_binary_fixed")),
    error = identity)
  expect_s3_class(result, "stormbreaker_incomplete_robust_logistic_fit")
})

test_that("nnet warnings do not replace its convergence check", {
  skip_if_not_installed("nnet")
  local_mocked_bindings(nnet = function(...) {
    warning("native diagnostic")
    list(convergence = 1L)
  }, .package = "nnet")
  seen <- character()
  result <- tryCatch(withCallingHandlers(sb_fit_classical_binary(
    data.frame(x = 1:8), rep(c(0L, 1L), 4L), "nnet_binary",
    sb_classical_binary_params("nnet_binary")), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  expect_identical(seen, "native diagnostic")
  expect_s3_class(result, "stormbreaker_nnet_nonconvergence")
})

test_that("BART warnings do not bypass the posterior-tree check", {
  skip_if_not_installed("BART")
  local_mocked_bindings(pbart = function(...) {
    warning("native diagnostic")
    list()
  }, .package = "BART")
  seen <- character()
  result <- tryCatch(withCallingHandlers(sb_fit_bart_binary(
    matrix(1:8, ncol = 1L), rep(c(0L, 1L), 4L), "bart_binary",
    sb_classical_binary_params("bart_binary")), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  expect_identical(seen, "native diagnostic")
  expect_s3_class(result, "stormbreaker_incomplete_bart_binary_fit")
})

# ---- pipeline-preparation ----

test_that("one pipeline construction prepares its registry once", {
  native <- sb_pipeline_registry_data
  count <- 0L
  expected <- sb_pipeline("none", "coxph")
  actual <- with_mocked_bindings(
    sb_pipeline("none", "coxph"),
    sb_pipeline_registry_data = function() {
      count <<- count + 1L
      native()
    }
  )
  expect_identical(actual, expected)
  expect_identical(count, 1L)
})

test_that("registry helpers agree for local and standalone lookup", {
  registry <- sb_pipeline_registry_data()
  for (learner in c("coxph", "logistic", "glmnet_cox", "glmnet_logistic")) {
    expect_identical(
      sb_pipeline_registry_match("none", learner, registry),
      sb_pipeline_registry_match("none", learner)
    )
    expect_identical(
      sb_pipeline_uses_inner_cv(learner, registry),
      sb_pipeline_uses_inner_cv(learner)
    )
  }
  expect_identical(
    sb_pipeline_registry_match("unregistered", "coxph", registry),
    sb_pipeline_registry_match("unregistered", "coxph")
  )
})

# ---- pipeline ----

test_that("the minimum pipeline has one explicit selector and learner", {
  pipeline <- sb_pipeline()

  expect_s3_class(pipeline, "stormbreaker_pipeline")
  expect_identical(pipeline$id, "none__coxph")
  expect_identical(pipeline$task_type, "survival")
  expect_identical(pipeline$selector, "none")
  expect_identical(pipeline$learner, "coxph")
  expect_identical(pipeline$risk_type, "linear_predictor")
  expect_identical(pipeline$risk_reference, "zero")
  expect_identical(pipeline$risk_formula, "x_beta")
  expect_identical(pipeline$risk_direction, "higher_is_higher_hazard")
})

test_that("the minimum binary pipeline predicts class-one probability", {
  pipeline <- sb_pipeline(learner = "logistic")

  expect_s3_class(pipeline, "stormbreaker_pipeline")
  expect_identical(pipeline$id, "none__logistic")
  expect_identical(pipeline$task_type, "binary")
  expect_identical(pipeline$selector, "none")
  expect_identical(pipeline$learner, "logistic")
  expect_identical(pipeline$prediction_type, "probability")
  expect_identical(pipeline$positive_class, 1L)
  expect_identical(
    pipeline$prediction_direction,
    "higher_is_more_likely_positive"
  )
  expect_identical(pipeline$fallback, "forbidden")
})

test_that("nested glmnet logistic declares fixed alpha and inner resampling", {
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 5L,
    inner_seed = 1701L
  )

  expect_identical(pipeline$task_type, "binary")
  expect_identical(pipeline$selector, "none")
  expect_identical(pipeline$learner, "glmnet_logistic")
  expect_identical(pipeline$prediction_type, "probability")
  expect_identical(pipeline$positive_class, 1L)
  expect_identical(
    pipeline$prediction_direction,
    "higher_is_more_likely_positive"
  )
  expect_identical(
    pipeline$learner_params,
    list(
      alpha = 1,
      standardize = TRUE,
      inner_folds = 5L,
      inner_seed = 1701L,
      selection_rule = "lambda.min",
      measure = "deviance"
    )
  )
  expect_identical(
    pipeline$id,
    "none__glmnet_logistic__alpha1__inner5__seed1701"
  )
})

test_that("glmnet logistic IDs distinguish Ridge and Elastic Net", {
  ridge <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0,
    inner_folds = 3L,
    inner_seed = 101L
  )
  elastic <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0.4,
    inner_folds = 3L,
    inner_seed = 202L
  )

  expect_identical(
    ridge$id,
    "none__glmnet_logistic__alpha0__inner3__seed101"
  )
  expect_identical(
    elastic$id,
    "none__glmnet_logistic__alpha0.4__inner3__seed202"
  )
  expect_identical(ridge$learner_params$alpha, 0)
  expect_identical(elastic$learner_params$alpha, 0.4)
})

test_that("glmnet standardization is explicit without changing old defaults", {
  default <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 101L
  )
  explicit_default <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_standardize = TRUE
  )
  disabled <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_standardize = FALSE
  )

  expect_identical(explicit_default, default)
  expect_identical(disabled$learner_params$standardize, FALSE)
  expect_identical(
    disabled$id,
    paste0(
      "none__glmnet_logistic__alpha1__standardize0",
      "__inner3__seed101"
    )
  )
  expect_identical(
    sb_validate_pipeline_contract(disabled, "binary"),
    disabled
  )

  for (value in list(NA, 0, "FALSE", c(TRUE, FALSE))) {
    error <- capture_sb_error(sb_pipeline(
      learner = "glmnet_logistic",
      inner_folds = 3L,
      inner_seed = 101L,
      glmnet_standardize = value
    ))
    expect_s3_class(error, "stormbreaker_bad_glmnet_standardize")
  }
  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    glmnet_standardize = FALSE
  ))
  expect_s3_class(irrelevant, "stormbreaker_irrelevant_glmnet_parameter")
})

test_that("glmnet lambda rule is predeclared without changing old defaults", {
  default <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L
  )
  explicit_default <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_lambda_rule = "lambda.min"
  )
  one_se <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_lambda_rule = "lambda.1se"
  )
  one_se_without_standardization <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_standardize = FALSE,
    glmnet_lambda_rule = "lambda.1se"
  )
  named_default <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_lambda_rule = c(rule = "lambda.min")
  )
  named_one_se <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_lambda_rule = c(rule = "lambda.1se")
  )

  expect_identical(explicit_default, default)
  expect_identical(named_default, default)
  expect_identical(named_one_se, one_se)
  expect_false(identical(named_default$id, one_se$id))
  expect_identical(one_se$learner_params$selection_rule, "lambda.1se")
  expect_identical(
    one_se$id,
    "none__glmnet_cox__alpha1__lambda1se__inner3__seed101"
  )
  expect_identical(
    one_se_without_standardization$id,
    paste0(
      "none__glmnet_cox__alpha1__standardize0__lambda1se",
      "__inner3__seed101"
    )
  )
  expect_identical(sb_validate_pipeline_contract(one_se, "survival"), one_se)

  for (value in list(NA_character_, "min", 1, c("lambda.min", "lambda.1se"))) {
    error <- capture_sb_error(sb_pipeline(
      learner = "glmnet_cox",
      inner_folds = 3L,
      inner_seed = 101L,
      glmnet_lambda_rule = value
    ))
    expect_s3_class(error, "stormbreaker_bad_glmnet_lambda_rule")
  }
  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "coxph",
    glmnet_lambda_rule = "lambda.1se"
  ))
  expect_s3_class(irrelevant, "stormbreaker_irrelevant_glmnet_parameter")
})

test_that("glmnet tuning measure is predeclared and task specific", {
  default <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 101L
  )
  explicit_default <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_tuning_measure = "deviance"
  )
  auc <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_tuning_measure = c(metric = "auc")
  )
  concordance <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_tuning_measure = "C"
  )

  expect_identical(explicit_default, default)
  expect_identical(auc$learner_params$measure, "auc")
  expect_identical(concordance$learner_params$measure, "C")
  expect_identical(
    auc$id,
    "none__glmnet_logistic__alpha1__measureauc__inner3__seed101"
  )
  expect_identical(
    concordance$id,
    "none__glmnet_cox__alpha1__measureC__inner3__seed101"
  )
  expect_identical(sb_validate_pipeline_contract(auc, "binary"), auc)
  expect_identical(
    sb_validate_pipeline_contract(concordance, "survival"),
    concordance
  )

  binary_error <- capture_sb_error(sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_tuning_measure = "C"
  ))
  survival_error <- capture_sb_error(sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_tuning_measure = "auc"
  ))
  malformed_error <- capture_sb_error(sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_tuning_measure = c("deviance", "C")
  ))
  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "coxph",
    glmnet_tuning_measure = "C"
  ))
  expect_s3_class(binary_error, "stormbreaker_bad_glmnet_tuning_measure")
  expect_s3_class(survival_error, "stormbreaker_bad_glmnet_tuning_measure")
  expect_s3_class(malformed_error, "stormbreaker_bad_glmnet_tuning_measure")
  expect_s3_class(irrelevant, "stormbreaker_irrelevant_glmnet_parameter")
})

test_that("nested glmnet logistic requires safe declared components", {
  missing_error <- capture_sb_error(
    sb_pipeline(learner = "glmnet_logistic")
  )
  expect_s3_class(missing_error, "stormbreaker_missing_inner_parameter")

  selector_error <- capture_sb_error(sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "glmnet_logistic",
    selector_n_features = 2L,
    inner_folds = 3L,
    inner_seed = 101L
  ))
  expect_s3_class(
    selector_error,
    "stormbreaker_incompatible_binary_selector"
  )
})

test_that("the binary minimum rejects unresolved selection and tuning", {
  selector_error <- capture_sb_error(sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "logistic",
    selector_n_features = 1L
  ))
  expect_s3_class(
    selector_error,
    "stormbreaker_incompatible_binary_selector"
  )

  inner_error <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    inner_folds = 3L,
    inner_seed = 1L
  ))
  expect_s3_class(inner_error, "stormbreaker_irrelevant_inner_parameter")

  penalty_error <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    glmnet_alpha = 1
  ))
  expect_s3_class(
    penalty_error,
    "stormbreaker_irrelevant_glmnet_parameter"
  )
})

test_that("unsupported components stop instead of being substituted", {
  selector_error <- capture_sb_error(sb_pipeline(selector = "lasso"))
  expect_s3_class(selector_error, "stormbreaker_unsupported_selector")

  learner_error <- capture_sb_error(sb_pipeline(learner = "random_forest"))
  expect_s3_class(learner_error, "stormbreaker_unsupported_learner")
})

test_that("a fixed Top-K univariable Cox pipeline is explicit", {
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "coxph",
    selector_n_features = 2L
  )

  expect_s3_class(pipeline, "stormbreaker_pipeline")
  expect_identical(pipeline$id, "univ_cox_top_k__coxph__k2")
  expect_identical(pipeline$selector, "univ_cox_top_k")
  expect_identical(pipeline$selector_params, list(n_features = 2L))
  expect_identical(pipeline$learner, "coxph")
  expect_identical(pipeline$fallback, "forbidden")
})

test_that("Top-K selector parameters are required and predeclared", {
  missing_error <- capture_sb_error(
    sb_pipeline(selector = "univ_cox_top_k")
  )
  expect_s3_class(missing_error, "stormbreaker_missing_selector_parameter")

  invalid_error <- capture_sb_error(
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 1.5
    )
  )
  expect_s3_class(invalid_error, "stormbreaker_bad_integer")

  irrelevant_error <- capture_sb_error(
    sb_pipeline(selector = "none", selector_n_features = 2L)
  )
  expect_s3_class(
    irrelevant_error,
    "stormbreaker_irrelevant_selector_parameter"
  )
})

test_that("a nested glmnet Cox pipeline declares its inner resampling", {
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 5L,
    inner_seed = 1701L
  )

  expect_identical(pipeline$learner, "glmnet_cox")
  expect_identical(
    pipeline$learner_params,
    list(
      alpha = 1,
      standardize = TRUE,
      inner_folds = 5L,
      inner_seed = 1701L,
      selection_rule = "lambda.min",
      measure = "deviance"
    )
  )
  expect_identical(pipeline$fallback, "forbidden")
  expect_identical(
    pipeline$id,
    "none__glmnet_cox__alpha1__inner5__seed1701"
  )
})

test_that("glmnet pipeline IDs distinguish alpha and inner seeds", {
  first <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 101L
  )
  second <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0.9,
    inner_folds = 3L,
    inner_seed = 202L
  )

  expect_identical(
    first$id,
    "none__glmnet_cox__alpha0.5__inner3__seed101"
  )
  expect_identical(
    second$id,
    "none__glmnet_cox__alpha0.9__inner3__seed202"
  )
  expect_identical(first$learner_params$alpha, 0.5)
  expect_identical(second$learner_params$alpha, 0.9)
  expect_identical(identical(first$id, second$id), FALSE)
})

test_that("glmnet alpha is finite and bounded", {
  invalid_values <- list(-0.1, 1.1, NA_real_, Inf, "0.5", c(0.3, 0.7))

  for (value in invalid_values) {
    error <- capture_sb_error(sb_pipeline(
      learner = "glmnet_cox",
      glmnet_alpha = value,
      inner_folds = 3L,
      inner_seed = 101L
    ))
    expect_s3_class(error, "stormbreaker_bad_glmnet_alpha")
  }
})

test_that("alpha zero is an explicit Ridge Cox pipeline", {
  ridge <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0,
    inner_folds = 5L,
    inner_seed = 20260817L
  )

  expect_identical(ridge$learner_params$alpha, 0)
  expect_identical(
    ridge$id,
    "none__glmnet_cox__alpha0__inner5__seed20260817"
  )
  printed <- paste(capture.output(print(ridge)), collapse = "\n")
  expect_match(printed, "Penalty family: Ridge", fixed = TRUE)
})

test_that("glmnet alpha is rejected for ordinary Cox", {
  error <- capture_sb_error(sb_pipeline(glmnet_alpha = 0.5))

  expect_s3_class(error, "stormbreaker_irrelevant_glmnet_parameter")
})

test_that("nested glmnet Cox admits fold-local Top-K selection", {
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "glmnet_cox",
    selector_n_features = 2L,
    inner_folds = 5L,
    inner_seed = 1701L
  )

  expect_s3_class(pipeline, "stormbreaker_pipeline")
  expect_identical(
    pipeline$id,
    "univ_cox_top_k__glmnet_cox__k2__alpha1__inner5__seed1701"
  )

  too_small <- capture_sb_error(sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "glmnet_cox",
    selector_n_features = 1L,
    inner_folds = 5L,
    inner_seed = 1701L
  ))
  expect_s3_class(
    too_small,
    "stormbreaker_glmnet_too_few_selected_features"
  )
})

test_that("inner resampling arguments are required only for glmnet Cox", {
  missing_error <- capture_sb_error(sb_pipeline(learner = "glmnet_cox"))
  expect_s3_class(missing_error, "stormbreaker_missing_inner_parameter")

  irrelevant_error <- capture_sb_error(sb_pipeline(inner_folds = 5L))
  expect_s3_class(
    irrelevant_error,
    "stormbreaker_irrelevant_inner_parameter"
  )

  invalid_error <- capture_sb_error(sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 2L,
    inner_seed = 1701L
  ))
  expect_s3_class(invalid_error, "stormbreaker_bad_integer")
})

test_that("binary AUC Top-K is an explicit fold-local logistic pipeline", {
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 3L
  )

  expect_s3_class(pipeline, "stormbreaker_pipeline")
  expect_identical(
    pipeline$id,
    "univ_binary_auc_top_k__logistic__k3"
  )
  expect_identical(pipeline$task_type, "binary")
  expect_identical(pipeline$selector, "univ_binary_auc_top_k")
  expect_identical(pipeline$selector_params, list(n_features = 3L))
  expect_identical(pipeline$learner, "logistic")
  expect_identical(pipeline$prediction_type, "probability")
  expect_identical(pipeline$positive_class, 1L)
  expect_identical(pipeline$fallback, "forbidden")
})

test_that("binary AUC Top-K requires a predeclared positive K", {
  missing_error <- capture_sb_error(sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic"
  ))
  expect_s3_class(missing_error, "stormbreaker_missing_selector_parameter")

  invalid_error <- capture_sb_error(sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 0L
  ))
  expect_s3_class(invalid_error, "stormbreaker_bad_integer")
})

test_that("binary AUC Top-K is endpoint-specific and fold-local in glmnet", {
  survival_error <- capture_sb_error(sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "coxph",
    selector_n_features = 2L
  ))
  expect_s3_class(
    survival_error,
    "stormbreaker_incompatible_survival_selector"
  )

  nested <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "glmnet_logistic",
    selector_n_features = 2L,
    inner_folds = 3L,
    inner_seed = 101L
  )
  expect_identical(
    nested$id,
    paste0(
      "univ_binary_auc_top_k__glmnet_logistic__k2",
      "__alpha1__inner3__seed101"
    )
  )
})

test_that("the shared pipeline validator reconstructs every admitted family", {
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      learner = "coxph",
      selector_n_features = 2L
    ),
    sb_pipeline(
      learner = "glmnet_cox",
      glmnet_alpha = 0.5,
      inner_folds = 3L,
      inner_seed = 101L
    ),
    sb_pipeline(learner = "logistic"),
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "logistic",
      selector_n_features = 2L
    ),
    sb_pipeline(
      learner = "glmnet_logistic",
      glmnet_alpha = 0.5,
      inner_folds = 3L,
      inner_seed = 101L
    ),
    sb_pipeline(
      learner = "gaussian_naive_bayes",
      naive_bayes_var_smoothing = 1e-8
    ),
    sb_pipeline(
      learner = "lightgbm_binary",
      lightgbm_nrounds = 7L,
      lightgbm_learning_rate = 0.05,
      lightgbm_num_leaves = 9L,
      lightgbm_min_data_in_leaf = 3L,
      lightgbm_seed = 17L
    ),
    sb_pipeline(
      learner = "catboost_binary",
      catboost_iterations = 7L,
      catboost_learning_rate = 0.05,
      catboost_depth = 3L,
      catboost_l2_leaf_reg = 4,
      catboost_seed = 17L
    ),
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "ranger_binary",
      selector_n_features = 2L,
      ranger_num_trees = 7L,
      ranger_mtry = 1L,
      ranger_min_node_size = 10L,
      ranger_seed = 17L
    ),
    sb_pipeline(
      selector = "univ_cox_top_k",
      learner = "catboost_cox",
      selector_n_features = 2L,
      catboost_iterations = 7L,
      catboost_learning_rate = 0.05,
      catboost_depth = 3L,
      catboost_l2_leaf_reg = 4,
      catboost_seed = 17L
    ),
    sb_pipeline(
      selector = "univ_cox_top_k",
      learner = "ranger_survival",
      selector_n_features = 2L,
      ranger_num_trees = 7L,
      ranger_mtry = 1L,
      ranger_min_node_size = 3L,
      ranger_seed = 17L
    ),
    sb_pipeline(
      learner = "coxboost",
      coxboost_stepno = 7L,
      coxboost_penalty = 10
    ),
    sb_pipeline(
      learner = "plsrcox",
      plsrcox_components = 1L
    ),
    sb_pipeline(
      learner = "superpc",
      superpc_threshold = 0,
      superpc_components = 1L
    ),
    sb_pipeline(
      learner = "survivalsvm",
      survivalsvm_gamma_mu = 0.5
    )
  )

  for (pipeline in pipelines) {
    observed <- sb_validate_pipeline_contract(
      pipeline,
      pipeline$task_type
    )
    expect_identical(observed, pipeline)
  }
})

test_that("the shared pipeline validator rejects shortcuts and mutations", {
  list_error <- capture_sb_error(sb_validate_pipeline_contract(
    list(sb_pipeline()),
    "survival"
  ))
  expect_s3_class(list_error, "stormbreaker_bad_pipeline")

  mismatch_error <- capture_sb_error(sb_validate_pipeline_contract(
    sb_pipeline(learner = "logistic"),
    "survival"
  ))
  expect_s3_class(
    mismatch_error,
    "stormbreaker_pipeline_task_mismatch"
  )

  mutated <- sb_pipeline()
  mutated$risk_formula <- "different_formula"
  mutation_error <- capture_sb_error(sb_validate_pipeline_contract(
    mutated,
    "survival"
  ))
  expect_s3_class(mutation_error, "stormbreaker_unresolved_pipeline")
})

test_that("Gaussian Naive Bayes declares fixed variance smoothing", {
  pipeline <- sb_pipeline(learner = "gaussian_naive_bayes")

  expect_identical(pipeline$task_type, "binary")
  expect_identical(pipeline$selector, "none")
  expect_identical(pipeline$learner, "gaussian_naive_bayes")
  expect_identical(
    pipeline$id,
    "none__gaussian_naive_bayes__smooth1e-09"
  )
  expect_identical(
    pipeline$learner_params,
    list(
      var_smoothing = 1e-9,
      variance_estimator = "population",
      prior = "empirical"
    )
  )
  expect_identical(pipeline$prediction_type, "probability")
  expect_identical(pipeline$positive_class, 1L)
  expect_identical(pipeline$fallback, "forbidden")
  printed <- paste(capture.output(print(pipeline)), collapse = "\n")
  expect_match(printed, "Variance smoothing: 1e-09", fixed = TRUE)
})

test_that("Gaussian Naive Bayes IDs distinguish smoothing and selector", {
  no_smoothing <- sb_pipeline(
    learner = "gaussian_naive_bayes",
    naive_bayes_var_smoothing = 0
  )
  selected <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "gaussian_naive_bayes",
    selector_n_features = 3L,
    naive_bayes_var_smoothing = 1e-6
  )

  expect_identical(
    no_smoothing$id,
    "none__gaussian_naive_bayes__smooth0"
  )
  expect_identical(
    selected$id,
    "univ_binary_auc_top_k__gaussian_naive_bayes__k3__smooth1e-06"
  )
  expect_identical(selected$selector_params, list(n_features = 3L))
  expect_identical(selected$learner_params$var_smoothing, 1e-6)
})

test_that("Gaussian Naive Bayes IDs preserve adjacent smoothing values", {
  smoothing <- c(1.0000000000000001e-9, 1.0000000000000003e-9)
  pipelines <- lapply(smoothing, function(value) {
    sb_pipeline(
      learner = "gaussian_naive_bayes",
      naive_bayes_var_smoothing = value
    )
  })
  encoded <- sub("^.*__smooth", "", vapply(pipelines, `[[`, "", "id"))

  expect_false(identical(smoothing[[1L]], smoothing[[2L]]))
  expect_false(identical(pipelines[[1L]]$id, pipelines[[2L]]$id))
  expect_identical(as.numeric(encoded), smoothing)
})

test_that("Gaussian Naive Bayes parameters and selectors are explicit", {
  invalid_values <- list(-1e-9, NA_real_, Inf, "1e-9", c(1e-9, 1e-8))
  for (value in invalid_values) {
    error <- capture_sb_error(sb_pipeline(
      learner = "gaussian_naive_bayes",
      naive_bayes_var_smoothing = value
    ))
    expect_s3_class(
      error,
      "stormbreaker_bad_naive_bayes_var_smoothing"
    )
  }

  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    naive_bayes_var_smoothing = 1e-9
  ))
  expect_s3_class(
    irrelevant,
    "stormbreaker_irrelevant_naive_bayes_parameter"
  )

  survival_selector <- capture_sb_error(sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "gaussian_naive_bayes",
    selector_n_features = 1L
  ))
  expect_s3_class(
    survival_selector,
    "stormbreaker_incompatible_binary_selector"
  )
})

test_that("fixed LightGBM declares its complete binary configuration", {
  pipeline <- sb_pipeline(learner = "lightgbm_binary")

  expect_identical(pipeline$task_type, "binary")
  expect_identical(pipeline$selector, "none")
  expect_identical(pipeline$learner, "lightgbm_binary")
  expect_identical(
    pipeline$id,
    paste0(
      "none__lightgbm_binary__rounds100__lr0.1",
      "__leaves31__minleaf20__seed1"
    )
  )
  expect_identical(pipeline$learner_params$nrounds, 100L)
  expect_identical(pipeline$learner_params$learning_rate, 0.1)
  expect_identical(pipeline$learner_params$num_leaves, 31L)
  expect_identical(pipeline$learner_params$min_data_in_leaf, 20L)
  expect_identical(pipeline$learner_params$seed, 1L)
  expect_identical(pipeline$learner_params$objective, "binary")
  expect_identical(pipeline$learner_params$boosting, "gbdt")
  expect_identical(pipeline$learner_params$device_type, "cpu")
  expect_identical(pipeline$learner_params$num_threads, 1L)
  expect_identical(pipeline$learner_params$deterministic, TRUE)
  expect_identical(pipeline$learner_params$force_col_wise, TRUE)
  expect_identical(pipeline$learner_params$feature_fraction, 1)
  expect_identical(pipeline$learner_params$bagging_fraction, 1)
  expect_identical(pipeline$learner_params$bagging_freq, 0L)
  expect_identical(pipeline$prediction_type, "probability")
  expect_identical(pipeline$positive_class, 1L)
  expect_identical(pipeline$fallback, "forbidden")
  printed <- paste(capture.output(print(pipeline)), collapse = "\n")
  expect_match(printed, "Fixed LightGBM: 100 rounds", fixed = TRUE)
  expect_match(printed, "learning rate 0.1", fixed = TRUE)
})

test_that("fixed LightGBM IDs encode selector and declared parameters", {
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "lightgbm_binary",
    selector_n_features = 3L,
    lightgbm_nrounds = 12L,
    lightgbm_learning_rate = 0.025,
    lightgbm_num_leaves = 7L,
    lightgbm_min_data_in_leaf = 3L,
    lightgbm_seed = 99L
  )

  expect_identical(
    pipeline$id,
    paste0(
      "univ_binary_auc_top_k__lightgbm_binary__k3",
      "__rounds12__lr0.025__leaves7__minleaf3__seed99"
    )
  )
  expect_identical(pipeline$selector_params, list(n_features = 3L))
  expect_identical(pipeline$learner_params$nrounds, 12L)
  expect_identical(pipeline$learner_params$learning_rate, 0.025)
  expect_identical(pipeline$learner_params$num_leaves, 7L)
  expect_identical(pipeline$learner_params$min_data_in_leaf, 3L)
  expect_identical(pipeline$learner_params$seed, 99L)
})

test_that("fixed LightGBM IDs preserve adjacent learning rates", {
  rates <- c(0.10000000000000001, 0.10000000000000003)
  pipelines <- lapply(rates, function(rate) {
    sb_pipeline(
      learner = "lightgbm_binary",
      lightgbm_learning_rate = rate
    )
  })
  ids <- vapply(pipelines, `[[`, "", "id")
  encoded <- sub("^.*__lr", "", ids)
  encoded <- sub("__leaves.*$", "", encoded)

  expect_false(identical(rates[[1L]], rates[[2L]]))
  expect_false(identical(ids[[1L]], ids[[2L]]))
  expect_identical(as.numeric(encoded), rates)
})

test_that("fixed LightGBM parameters are validated and learner-specific", {
  integer_cases <- list(
    list(lightgbm_nrounds = 0L),
    list(lightgbm_num_leaves = 1L),
    list(lightgbm_min_data_in_leaf = -1L),
    list(lightgbm_seed = 0L),
    list(lightgbm_nrounds = 1.5)
  )
  for (arguments in integer_cases) {
    error <- capture_sb_error(do.call(
      sb_pipeline,
      c(list(learner = "lightgbm_binary"), arguments)
    ))
    expect_s3_class(error, "stormbreaker_bad_integer")
  }

  invalid_rates <- list(0, -0.1, NA_real_, Inf, "0.1", c(0.1, 0.2))
  for (rate in invalid_rates) {
    error <- capture_sb_error(sb_pipeline(
      learner = "lightgbm_binary",
      lightgbm_learning_rate = rate
    ))
    expect_s3_class(error, "stormbreaker_bad_lightgbm_learning_rate")
  }

  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    lightgbm_nrounds = 10L
  ))
  expect_s3_class(
    irrelevant,
    "stormbreaker_irrelevant_lightgbm_parameter"
  )

  survival_selector <- capture_sb_error(sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "lightgbm_binary",
    selector_n_features = 1L
  ))
  expect_s3_class(
    survival_selector,
    "stormbreaker_incompatible_binary_selector"
  )
})

test_that("fixed CatBoost declares its complete binary configuration", {
  pipeline <- sb_pipeline(learner = "catboost_binary")

  expect_identical(pipeline$task_type, "binary")
  expect_identical(pipeline$selector, "none")
  expect_identical(pipeline$learner, "catboost_binary")
  expect_identical(
    pipeline$id,
    paste0(
      "none__catboost_binary__iter100__lr0.1",
      "__depth6__l2reg3__seed1"
    )
  )
  expect_identical(pipeline$learner_params$iterations, 100L)
  expect_identical(pipeline$learner_params$learning_rate, 0.1)
  expect_identical(pipeline$learner_params$depth, 6L)
  expect_identical(pipeline$learner_params$l2_leaf_reg, 3)
  expect_identical(pipeline$learner_params$random_seed, 1L)
  expect_identical(pipeline$learner_params$loss_function, "Logloss")
  expect_identical(pipeline$learner_params$eval_metric, "Logloss")
  expect_identical(pipeline$learner_params$task_type, "CPU")
  expect_identical(pipeline$learner_params$thread_count, 1L)
  expect_identical(pipeline$learner_params$bootstrap_type, "Bayesian")
  expect_identical(pipeline$learner_params$bagging_temperature, 1)
  expect_identical(pipeline$learner_params$random_strength, 1)
  expect_identical(pipeline$learner_params$rsm, 1)
  expect_identical(pipeline$learner_params$grow_policy, "SymmetricTree")
  expect_identical(pipeline$learner_params$boosting_type, "Plain")
  expect_identical(pipeline$learner_params$use_best_model, FALSE)
  expect_identical(pipeline$learner_params$allow_writing_files, FALSE)
  expect_identical(pipeline$learner_params$save_snapshot, FALSE)
  expect_identical(pipeline$prediction_type, "probability")
  expect_identical(pipeline$positive_class, 1L)
  expect_identical(pipeline$fallback, "forbidden")
  printed <- paste(capture.output(print(pipeline)), collapse = "\n")
  expect_match(printed, "Fixed CatBoost: 100 iterations", fixed = TRUE)
  expect_match(printed, "learning rate 0.1", fixed = TRUE)
})

test_that("fixed CatBoost IDs encode selector and declared parameters", {
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "catboost_binary",
    selector_n_features = 3L,
    catboost_iterations = 12L,
    catboost_learning_rate = 0.025,
    catboost_depth = 4L,
    catboost_l2_leaf_reg = 2.5,
    catboost_seed = 99L
  )

  expect_identical(
    pipeline$id,
    paste0(
      "univ_binary_auc_top_k__catboost_binary__k3",
      "__iter12__lr0.025__depth4__l2reg2.5__seed99"
    )
  )
  expect_identical(pipeline$selector_params, list(n_features = 3L))
  expect_identical(pipeline$learner_params$iterations, 12L)
  expect_identical(pipeline$learner_params$learning_rate, 0.025)
  expect_identical(pipeline$learner_params$depth, 4L)
  expect_identical(pipeline$learner_params$l2_leaf_reg, 2.5)
  expect_identical(pipeline$learner_params$random_seed, 99L)
})

test_that("fixed CatBoost IDs preserve adjacent numeric settings", {
  rates <- c(0.10000000000000001, 0.10000000000000003)
  regularization <- c(3, 3.0000000000000004)
  pipelines <- Map(function(rate, l2) {
    sb_pipeline(
      learner = "catboost_binary",
      catboost_learning_rate = rate,
      catboost_l2_leaf_reg = l2
    )
  }, rates, regularization)
  ids <- vapply(pipelines, `[[`, "", "id")
  encoded_rates <- sub("^.*__lr", "", ids)
  encoded_rates <- sub("__depth.*$", "", encoded_rates)
  encoded_l2 <- sub("^.*__l2reg", "", ids)
  encoded_l2 <- sub("__seed.*$", "", encoded_l2)

  expect_false(identical(rates[[1L]], rates[[2L]]))
  expect_false(identical(regularization[[1L]], regularization[[2L]]))
  expect_false(identical(ids[[1L]], ids[[2L]]))
  expect_identical(as.numeric(encoded_rates), rates)
  expect_identical(as.numeric(encoded_l2), regularization)
})

test_that("fixed CatBoost parameters are validated and learner-specific", {
  integer_cases <- list(
    list(catboost_iterations = 0L),
    list(catboost_depth = 0L),
    list(catboost_seed = 0L),
    list(catboost_iterations = 1.5)
  )
  for (arguments in integer_cases) {
    error <- capture_sb_error(do.call(
      sb_pipeline,
      c(list(learner = "catboost_binary"), arguments)
    ))
    expect_s3_class(error, "stormbreaker_bad_integer")
  }

  depth_error <- capture_sb_error(sb_pipeline(
    learner = "catboost_binary",
    catboost_depth = 17L
  ))
  expect_s3_class(depth_error, "stormbreaker_bad_catboost_depth")

  invalid_rates <- list(0, -0.1, NA_real_, Inf, "0.1", c(0.1, 0.2))
  for (rate in invalid_rates) {
    error <- capture_sb_error(sb_pipeline(
      learner = "catboost_binary",
      catboost_learning_rate = rate
    ))
    expect_s3_class(error, "stormbreaker_bad_catboost_learning_rate")
  }

  invalid_l2 <- list(0, -1, NA_real_, Inf, "3", c(2, 3))
  for (l2 in invalid_l2) {
    error <- capture_sb_error(sb_pipeline(
      learner = "catboost_binary",
      catboost_l2_leaf_reg = l2
    ))
    expect_s3_class(error, "stormbreaker_bad_catboost_l2_leaf_reg")
  }

  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    catboost_iterations = 10L
  ))
  expect_s3_class(
    irrelevant,
    "stormbreaker_irrelevant_catboost_parameter"
  )

  survival_selector <- capture_sb_error(sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "catboost_binary",
    selector_n_features = 1L
  ))
  expect_s3_class(
    survival_selector,
    "stormbreaker_incompatible_binary_selector"
  )
})
})
