# ---- model-mboost ----

# ---- algorithm-boosting ----

sb_mboost_diagnostic_ids <- function() {
  c("gamboost_binary", "blackboost_binary",
    sb_mboost_cox_learner_ids(), sb_mboost_aft_learner_ids(),
    sb_mboost_ranking_learner_ids())
}

#' @rdname sb_gam_effects
#' @export
sb_boosting_diagnostics <- function(fit, type = "training_loss", max_pairs = 20L) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  learner <- fit$pipeline$learner
  if (learner == "bst_hinge_binary") return(sb_hinge_boosting_diagnostic(fit, type))
  if (learner %in% sb_classical_binary_boosting_learner_ids()) {
    return(sb_classical_boosting_diagnostic(fit, type))
  }
  if (learner %in% c("coxboost", "coxboost_cv")) {
    return(sb_coxboost_diagnostic(fit, type))
  }
  if (learner %in% c(sb_lightgbm_binary_learner_ids(), sb_catboost_diagnostic_ids())) {
    return(sb_light_cat_diagnostic(fit, type, max_pairs))
  }
  if (learner %in% c(sb_gbm_diagnostic_ids(), sb_xgboost_learner_ids())) {
    return(sb_gbm_xgboost_diagnostic(fit, type))
  }
  if (!learner %in% sb_mboost_diagnostic_ids()) sb_abort(
    "This learner has no supported native boosting diagnostic adapter.",
    "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("training_loss", "component_usage")) sb_abort(
    "Choose training_loss or component_usage.", "stormbreaker_unsupported_diagnostic")
  if (!requireNamespace("mboost", quietly = TRUE)) sb_abort(
    "Package mboost is required to read its saved diagnostic objects.",
    "stormbreaker_missing_dependency")
  state <- fit$learner_state
  engine <- unserialize(state$engine_raw)
  if (!inherits(engine, "mboost")) sb_abort("The original mboost engine is missing.",
    "stormbreaker_invalid_diagnostic")
  stop <- mboost::mstop(engine)
  if (!is.numeric(stop) || length(stop) != 1L || !is.finite(stop) ||
      stop < 1 || stop != floor(stop)) sb_abort("The original stopping iteration is invalid.",
    "stormbreaker_invalid_diagnostic")
  selections <- NULL
  if (type == "training_loss") {
    loss <- mboost::risk(engine)
    if (!is.numeric(loss) || length(loss) != stop + 1L || any(!is.finite(loss))) sb_abort(
      "The native in-bag loss sequence is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(iteration = seq.int(0L, stop), loss = unname(loss))
    measure <- paste("Native training objective:", engine$family@name)
  } else {
    if (startsWith(learner, "blackboost_")) sb_abort(
      "BlackBoost updates one tree base learner; component counts are not a variable-selection diagnostic.",
      "stormbreaker_unavailable_diagnostic")
    native_variables <- stats::variable.names(engine, usedonly = FALSE)
    selected <- mboost::selected(engine)
    if (!length(native_variables) || anyNA(native_variables) ||
        !is.numeric(selected) || length(selected) != stop || any(!is.finite(selected)) ||
        any(selected < 1 | selected > length(native_variables)) || any(selected != floor(selected))) sb_abort(
      "The original component selection sequence is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    variables <- unname(native_variables)
    if (!is.null(state$engine_features)) {
      index <- match(variables, state$engine_features)
      matched <- !is.na(index)
      variables[matched] <- state$features[index[matched]]
    }
    if (!all(variables %in% c(state$features, "(Intercept)"))) sb_abort(
      "A component cannot be mapped to its original input feature.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(component_id = seq_along(variables), feature = variables,
      native_component = names(native_variables),
      updates = tabulate(selected, nbins = length(variables)))
    selections <- data.frame(iteration = seq_along(selected), component_id = as.integer(selected),
      feature = variables[selected])
    measure <- "Native component update count"
  }
  structure(list(type = paste0("boosting_", type), table = tab, selections = selections,
    learner = learner, pipeline_id = fit$pipeline$id, stopping_iteration = as.integer(stop),
    measure = measure, source = "original saved mboost engine; no additional fitting",
    validation_scope = "in-bag development training, not held-out evaluation"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_boosting_diagnostic <- function(x) {
  iteration <- loss <- updates <- feature <- NULL
  if (x$type == "boosting_training_loss") {
    return(ggplot2::ggplot(x$table, ggplot2::aes(iteration, loss)) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::geom_point(data = x$table[c(1, nrow(x$table)), ], colour = "#2878A0", size = 1.5) +
      ggplot2::labs(title = "Boosting training loss", subtitle = x$learner,
        x = "Boosting iteration (zero = initial model)", y = "Native in-bag objective",
        caption = paste0(x$measure, ".\n",
          "Original scale and sign; not validation performance or a new stopping rule.")) + sb_viz_theme())
  }
  tab <- x$table
  tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$updates)])
  p <- ggplot2::ggplot(tab, ggplot2::aes(updates, feature)) +
    ggplot2::geom_col(fill = "#2878A0", width = .65) +
    ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
    ggplot2::labs(title = "Boosting component updates", subtitle = x$learner,
      x = "Number of iterations selecting this component", y = NULL,
      caption = "All original components retained, including zero counts and any intercept.\nUpdate frequency is not feature importance, effect magnitude or resampling stability.") + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 10, height = max(5, 2.5 + .28 * nrow(tab)))
  p
}

# ---- mboost-aft-fixed-adapters ----

sb_mboost_aft_learner_ids <- function() {
  c(
    "gamboost_weibull_aft",
    "gamboost_loglogistic_aft",
    "gamboost_lognormal_aft",
    "glmboost_weibull_aft",
    "glmboost_loglogistic_aft",
    "glmboost_lognormal_aft",
    "blackboost_weibull_aft",
    "blackboost_loglogistic_aft",
    "blackboost_lognormal_aft"
  )
}

sb_mboost_aft_params <- function(learner) {
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_mboost_aft_learner_ids()) {
    sb_abort(
      "The requested fixed mboost AFT learner is not supported.",
      "stormbreaker_unsupported_learner"
    )
  }

  distribution <- if (grepl("_weibull_", learner, fixed = TRUE)) {
    list(family = "Weibull", distribution = "weibull")
  } else if (grepl("_loglogistic_", learner, fixed = TRUE)) {
    list(family = "Loglog", distribution = "loglogistic")
  } else {
    list(family = "Lognormal", distribution = "lognormal")
  }
  common <- c(
    distribution,
    list(
      mstop = 100L,
      nu = 0.1,
      risk = "inbag",
      nuisance_range = c(0.1, 100)
    )
  )
  if (startsWith(learner, "gamboost_")) {
    return(c(common, list(
      algorithm = "gamboost",
      base_learner = "bbs",
      dfbase = 4
    )))
  }
  if (startsWith(learner, "glmboost_")) {
    return(c(common, list(
      algorithm = "glmboost",
      base_learner = "bolscw",
      center = TRUE
    )))
  }
  c(common, list(
    algorithm = "blackboost",
    base_learner = "btree",
    tree_controls = list(
      teststat = "quad",
      testtype = "Teststatistic",
      mincriterion = 0,
      minsplit = 10L,
      minbucket = 4L,
      maxdepth = 2L,
      saveinfo = FALSE
    ),
    seed = 1L
  ))
}

sb_mboost_aft_parameter_semantics <- function(learner) {
  params <- sb_mboost_aft_params(learner)
  paste0(
    params$algorithm, " ", params$distribution,
    " AFT likelihood, 100 boosting iterations, shrinkage 0.1, ",
    "positive scale-search interval [0.1,100], and ",
    params$base_learner, " base learner fixed before resampling"
  )
}

sb_mboost_aft_dependency <- function(learner) {
  params <- sb_mboost_aft_params(learner)
  packages <- if (identical(params$algorithm, "blackboost")) {
    c("mboost", "partykit")
  } else {
    "mboost"
  }
  list(
    packages = packages,
    available = all(vapply(
      packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    ))
  )
}

sb_require_mboost_aft_dependency <- function(learner) {
  dependency <- sb_mboost_aft_dependency(learner)
  if (!dependency$available) {
    missing <- dependency$packages[!vapply(
      dependency$packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )]
    sb_abort(
      paste0(
        "Package", if (length(missing) > 1L) "s" else "", " `",
        paste(missing, collapse = "`, `"), "` ",
        if (length(missing) > 1L) "are" else "is",
        " required for learner `", learner, "`."
      ),
      "stormbreaker_missing_mboost_aft_dependency"
    )
  }
  invisible(dependency)
}

sb_mboost_aft_internal_feature_names <- function(n_features) {
  paste0(".stormbreaker_mboost_aft_x", seq_len(n_features))
}

sb_mboost_aft_formula_data <- function(x, time, event) {
  engine_features <- sb_mboost_aft_internal_feature_names(ncol(x))
  engine_x <- x
  names(engine_x) <- engine_features
  formula_data <- sb_additional_survival_formula_data(
    engine_x,
    time,
    event
  )
  c(formula_data, list(engine_features = engine_features))
}

sb_mboost_aft_family <- function(params) {
  switch(params$family,
    Weibull = mboost::Weibull(params$nuisance_range),
    Loglog = mboost::Loglog(params$nuisance_range),
    Lognormal = mboost::Lognormal(params$nuisance_range)
  )
}

sb_fit_mboost_aft_engine <- function(formula_data, params) {
  if (identical(params$algorithm, "gamboost")) {
    withr::local_options(mboost_useMatrix = FALSE)
  }
  control <- mboost::boost_control(
    mstop = params$mstop,
    nu = params$nu,
    risk = params$risk,
    center = TRUE,
    trace = FALSE
  )
  fit_engine <- function() {
    family <- sb_mboost_aft_family(params)
    switch(params$algorithm,
      gamboost = mboost::gamboost(
        formula_data$formula,
        data = formula_data$data,
        family = family,
        control = control,
        baselearner = mboost::bbs,
        dfbase = params$dfbase
      ),
      glmboost = mboost::glmboost(
        formula_data$formula,
        data = formula_data$data,
        family = family,
        control = control,
        center = params$center
      ),
      blackboost = mboost::blackboost(
        formula_data$formula,
        data = formula_data$data,
        family = family,
        control = control,
        tree_controls = do.call(
          partykit::ctree_control,
          params$tree_controls
        )
      )
    )
  }
  fit <- function() {
    sb_classical_survival_stop_on_warning(
      fit_engine(),
      paste0("fixed ", params$algorithm, " ", params$distribution, " AFT learner"),
      "stormbreaker_mboost_aft_fit_warning"
    )
  }
  if ("seed" %in% names(params)) {
    return(withr::with_seed(params$seed, fit()))
  }
  fit()
}

sb_new_mboost_aft_model <- function(
  learner, features, engine_features, params, nuisance_scale, engine
) {
  structure(
    list(
      learner = learner,
      features = features,
      engine_features = engine_features,
      params = params,
      nuisance_scale = nuisance_scale,
      engine = engine
    ),
    class = c(
      paste0("stormbreaker_", learner),
      "stormbreaker_mboost_aft_fixed"
    )
  )
}

sb_fit_mboost_aft <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_mboost_aft_dependency(learner)
  formula_data <- sb_mboost_aft_formula_data(
    inputs$x,
    inputs$time,
    inputs$event
  )
  engine <- sb_fit_mboost_aft_engine(formula_data, params)
  nuisance_scale <- as.numeric(mboost::nuisance(engine))
  valid <- inherits(engine, params$algorithm) &&
    identical(as.integer(mboost::mstop(engine)), params$mstop) &&
    length(nuisance_scale) == 1L && !is.na(nuisance_scale) &&
    is.finite(nuisance_scale) &&
    nuisance_scale >= params$nuisance_range[[1L]] &&
    nuisance_scale <= params$nuisance_range[[2L]]
  if (!valid) {
    sb_abort(
      paste0("Learner `", learner, "` did not return its requested AFT fit."),
      "stormbreaker_incomplete_mboost_aft_fit"
    )
  }

  model <- sb_new_mboost_aft_model(
    learner,
    names(inputs$x),
    formula_data$engine_features,
    params,
    nuisance_scale,
    engine
  )
  sb_predict_mboost_aft(model, inputs$x)
  model
}

sb_validate_mboost_aft_model_x <- function(model, x) {
  valid_model <- inherits(model, "stormbreaker_mboost_aft_fixed") &&
    is.character(model$learner) && length(model$learner) == 1L &&
    !is.na(model$learner) &&
    model$learner %in% sb_mboost_aft_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.character(model$engine_features) &&
    length(model$engine_features) == length(model$features) &&
    !anyNA(model$engine_features) && !anyDuplicated(model$engine_features) &&
    is.numeric(model$nuisance_scale) &&
    length(model$nuisance_scale) == 1L &&
    !is.na(model$nuisance_scale) && is.finite(model$nuisance_scale) &&
    model$nuisance_scale > 0
  if (!valid_model) {
    sb_abort(
      "The supplied model is not a fitted mboost AFT learner.",
      "stormbreaker_bad_mboost_aft_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  sb_require_mboost_aft_dependency(model$learner)
  x
}

sb_mboost_aft_engine_x <- function(model, x) {
  engine_x <- x
  names(engine_x) <- model$engine_features
  engine_x
}

sb_mboost_aft_link <- function(model, x) {
  x <- sb_validate_mboost_aft_model_x(model, x)
  link <- withCallingHandlers(
    as.numeric(stats::predict(
      model$engine,
      newdata = sb_mboost_aft_engine_x(model, x),
      type = "link"
    )),
    warning = function(warning) {
      spline_extrapolation <- identical(
        model$params$algorithm,
        "gamboost"
      ) && identical(
        conditionMessage(warning),
        paste0(
          "Some 'x' values are beyond 'boundary.knots'; ",
          "Linear extrapolation used."
        )
      )
      if (spline_extrapolation) {
        invokeRestart("muffleWarning")
      }
      sb_abort(
        paste0(
          "The fixed ", model$learner, " prediction raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_mboost_aft_prediction_warning",
        parent = warning
      )
    }
  )
  if (length(link) != nrow(x) || anyNA(link) || any(!is.finite(link))) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid AFT scores."),
      "stormbreaker_invalid_mboost_aft_prediction"
    )
  }
  unname(link)
}

sb_predict_mboost_aft <- function(model, x) {
  -sb_mboost_aft_link(model, x)
}

sb_mboost_aft_survival_probability <- function(model, x, horizon) {
  link <- sb_mboost_aft_link(model, x)
  transformed <- (log(horizon) - link) / model$nuisance_scale
  probability <- switch(model$params$distribution,
    weibull = exp(-exp(transformed)),
    loglogistic = stats::plogis(-transformed),
    lognormal = stats::pnorm(transformed, lower.tail = FALSE)
  )
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid horizon survival."),
      "stormbreaker_invalid_mboost_aft_prediction"
    )
  }
  unname(probability)
}

sb_mboost_aft_direct_survival_source <- function(learner) {
  params <- sb_mboost_aft_params(learner)
  paste0(
    "mboost_", params$algorithm, "_",
    params$distribution, "_aft_direct"
  )
}

sb_mboost_aft_survival_source <- function(model) {
  sb_mboost_aft_direct_survival_source(model$learner)
}

sb_predict_mboost_aft_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  inputs <- sb_validate_additional_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  list(
    prediction_horizon = prediction_horizon,
    source = sb_mboost_aft_survival_source(model),
    ties = "not_applicable",
    probability = sb_mboost_aft_survival_probability(
      model,
      x_assessment,
      prediction_horizon
    )
  )
}

sb_extract_mboost_aft_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_mboost_aft_fixed")) {
    sb_bad_learner_state("The fitted mboost AFT learner state is invalid.")
  }
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_mboost_aft_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    survival_probability <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      source = probability$source,
      ties = probability$ties
    )
  }

  structure(
    list(
      learner = model$learner,
      features = model$features,
      engine_features = model$engine_features,
      params = model$params,
      nuisance_scale = model$nuisance_scale,
      risk_reference = paste0(
        "mboost_", model$params$algorithm, "_",
        model$params$distribution, "_aft_native"
      ),
      risk_formula = "negative_mboost_aft_log_time_link",
      model_serialization = "R_xdr_v3_mboost",
      engine_raw = serialize(
        model$engine,
        connection = NULL,
        ascii = FALSE,
        xdr = TRUE,
        version = 3L
      ),
      training_data_retention = "required_by_mboost_prediction_state",
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_mboost_aft_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_mboost_aft_model <- function(state, learner) {
  valid <- is.raw(state$engine_raw) && length(state$engine_raw) > 0L &&
    is.list(state$params) &&
    is.character(state$engine_features) &&
    length(state$engine_features) == length(state$features) &&
    !anyNA(state$engine_features) && !anyDuplicated(state$engine_features) &&
    is.numeric(state$nuisance_scale) &&
    length(state$nuisance_scale) == 1L &&
    !is.na(state$nuisance_scale) && is.finite(state$nuisance_scale) &&
    state$nuisance_scale > 0
  if (!valid) {
    sb_bad_learner_state("The saved mboost AFT learner state is invalid.")
  }
  sb_require_mboost_aft_dependency(learner)
  engine <- tryCatch(
    unserialize(state$engine_raw),
    error = function(error) {
      sb_abort(
        paste0("The saved `", learner, "` model could not be loaded."),
        "stormbreaker_invalid_mboost_aft_state",
        parent = error
      )
    }
  )
  if (!inherits(engine, state$params$algorithm) ||
    !sb_pipeline_parameters_equal(as.integer(mboost::mstop(engine)), state$params$mstop)) {
    sb_bad_learner_state("The saved mboost AFT learner engine is invalid.")
  }
  state$params$mstop <- as.integer(state$params$mstop)
  sb_new_mboost_aft_model(
    learner,
    state$features,
    state$engine_features,
    state$params,
    state$nuisance_scale,
    engine
  )
}

sb_predict_mboost_aft_state <- function(state, x) {
  learner <- if (is.list(state)) state$learner else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_mboost_aft_learner_ids()) {
    sb_bad_learner_state("The saved mboost AFT learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_restore_mboost_aft_model(state, learner)
  result <- list(risk = sb_predict_mboost_aft(model, x))
  probability <- state$survival_probability
  if (!is.null(probability)) {
    valid <- is.list(probability) &&
      is.numeric(probability$prediction_horizon) &&
      length(probability$prediction_horizon) == 1L &&
      !is.na(probability$prediction_horizon) &&
      is.finite(probability$prediction_horizon) &&
      probability$prediction_horizon > 0
    if (!valid) {
      sb_bad_learner_state("The saved mboost AFT horizon state is invalid.")
    }
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- sb_mboost_aft_survival_probability(
      model,
      x,
      probability$prediction_horizon
    )
  }
  result
}

sb_mboost_aft_learner_adapter <- function(
  learner, params = sb_mboost_aft_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_mboost_aft(x, time, event, learner, params)
    },
    predict = sb_predict_mboost_aft,
    predict_survival = sb_predict_mboost_aft_survival,
    extract_state = sb_extract_mboost_aft_state,
    predict_state = sb_predict_mboost_aft_state
  )
}

# ---- mboost-cox-fixed-adapters ----

sb_mboost_cox_learner_ids <- function() {
  c(
    "gamboost_cox",
    "glmboost_cox",
    "blackboost_cox"
  )
}

sb_mboost_cox_params <- function(learner) {
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_mboost_cox_learner_ids()) {
    sb_abort(
      "The requested fixed mboost Cox learner is not supported.",
      "stormbreaker_unsupported_learner"
    )
  }

  common <- list(
    mstop = 100L,
    nu = 0.1,
    risk = "inbag",
    family = "CoxPH"
  )
  specific <- switch(learner,
    gamboost_cox = list(
      base_learner = "bbs",
      dfbase = 4
    ),
    glmboost_cox = list(
      base_learner = "bolscw",
      center = TRUE
    ),
    blackboost_cox = list(
      base_learner = "btree",
      tree_controls = list(
        teststat = "quad",
        testtype = "Teststatistic",
        mincriterion = 0,
        minsplit = 10L,
        minbucket = 4L,
        maxdepth = 2L,
        saveinfo = FALSE
      )
    )
  )
  c(common, specific)
}

sb_mboost_cox_dependency <- function(learner) {
  packages <- if (identical(learner, "blackboost_cox")) {
    c("mboost", "partykit")
  } else {
    "mboost"
  }
  list(
    packages = packages,
    available = all(vapply(
      packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    ))
  )
}

sb_require_mboost_cox_dependency <- function(learner) {
  dependency <- sb_mboost_cox_dependency(learner)
  if (!dependency$available) {
    missing <- dependency$packages[!vapply(
      dependency$packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )]
    sb_abort(
      paste0(
        "Package", if (length(missing) > 1L) "s" else "", " `",
        paste(missing, collapse = "`, `"), "` ",
        if (length(missing) > 1L) "are" else "is",
        " required for learner `", learner, "`."
      ),
      "stormbreaker_missing_mboost_cox_dependency"
    )
  }
  invisible(dependency)
}

sb_new_mboost_cox_model <- function(learner, features, params, engine) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      engine = engine
    ),
    class = c(
      paste0("stormbreaker_", learner),
      "stormbreaker_mboost_cox_fixed"
    )
  )
}

sb_fit_mboost_cox <- function(x, time, event, learner, params) {
  inputs <- sb_validate_classical_survival_inputs(
    x,
    time,
    event,
    learner
  )
  sb_require_mboost_cox_dependency(learner)
  if (identical(learner, "gamboost_cox")) {
    # Native dense spline algebra avoids sparse-QR augmentation warnings;
    # the spline basis, degrees of freedom and warning checks are unchanged.
    withr::local_options(mboost_useMatrix = FALSE)
  }
  formula_data <- sb_additional_survival_formula_data(
    inputs$x,
    inputs$time,
    inputs$event
  )
  control <- mboost::boost_control(
    mstop = params$mstop,
    nu = params$nu,
    risk = params$risk,
    center = TRUE,
    trace = FALSE
  )

  engine <- sb_classical_survival_stop_on_warning(
    switch(learner,
      gamboost_cox = mboost::gamboost(
        formula_data$formula,
        data = formula_data$data,
        family = mboost::CoxPH(),
        control = control,
        baselearner = mboost::bbs,
        dfbase = params$dfbase
      ),
      glmboost_cox = mboost::glmboost(
        formula_data$formula,
        data = formula_data$data,
        family = mboost::CoxPH(),
        control = control,
        center = params$center
      ),
      blackboost_cox = mboost::blackboost(
        formula_data$formula,
        data = formula_data$data,
        family = mboost::CoxPH(),
        control = control,
        tree_controls = do.call(
          partykit::ctree_control,
          params$tree_controls
        )
      )
    ),
    paste0("fixed ", learner, " learner"),
    "stormbreaker_mboost_cox_fit_warning"
  )

  expected_class <- switch(learner,
    gamboost_cox = "gamboost",
    glmboost_cox = "glmboost",
    blackboost_cox = "blackboost"
  )
  if (!inherits(engine, expected_class) ||
    !identical(as.integer(mboost::mstop(engine)), params$mstop)) {
    sb_abort(
      paste0("Learner `", learner, "` did not return its requested fit."),
      "stormbreaker_incomplete_mboost_cox_fit"
    )
  }

  model <- sb_new_mboost_cox_model(
    learner,
    names(inputs$x),
    params,
    engine
  )
  sb_predict_mboost_cox(model, inputs$x)
  model
}

sb_validate_mboost_cox_model_x <- function(model, x) {
  valid_model <- inherits(model, "stormbreaker_mboost_cox_fixed") &&
    is.character(model$learner) && length(model$learner) == 1L &&
    !is.na(model$learner) &&
    model$learner %in% sb_mboost_cox_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features)
  if (!valid_model) {
    sb_abort(
      "The supplied model is not a fitted mboost Cox learner.",
      "stormbreaker_bad_mboost_cox_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  sb_require_mboost_cox_dependency(model$learner)
  x
}

sb_predict_mboost_cox <- function(model, x) {
  x <- sb_validate_mboost_cox_model_x(model, x)
  risk <- withCallingHandlers(
    as.numeric(stats::predict(
      model$engine,
      newdata = x,
      type = "link"
    )),
    warning = function(warning) {
      spline_extrapolation <- identical(
        model$learner,
        "gamboost_cox"
      ) && identical(
        conditionMessage(warning),
        paste0(
          "Some 'x' values are beyond 'boundary.knots'; ",
          "Linear extrapolation used."
        )
      )
      if (spline_extrapolation) {
        invokeRestart("muffleWarning")
      }
      sb_abort(
        paste0(
          "The fixed ", model$learner, " prediction raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_mboost_cox_prediction_warning",
        parent = warning
      )
    }
  )
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid risk scores."),
      "stormbreaker_invalid_mboost_cox_prediction"
    )
  }
  unname(risk)
}

sb_predict_mboost_cox_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_risk <- sb_predict_mboost_cox(model, inputs$x)
  assessment_risk <- sb_predict_mboost_cox(model, x_assessment)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = paste0(model$learner, "_breslow_native_link"),
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_extract_mboost_cox_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_mboost_cox_fixed")) {
    sb_bad_learner_state("The fitted mboost Cox learner state is invalid.")
  }

  baseline_hazard <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_mboost_cox_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    baseline_hazard <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      cumulative_baseline_hazard = probability$cumulative_baseline_hazard,
      baseline_survival = probability$baseline_survival,
      source = probability$source,
      ties = probability$ties
    )
  }

  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      risk_reference = "mboost_native",
      risk_formula = "mboost_coxph_link_prediction",
      model_serialization = "R_xdr_v3_mboost",
      engine_raw = serialize(
        model$engine,
        connection = NULL,
        ascii = FALSE,
        xdr = TRUE,
        version = 3L
      ),
      training_data_retention = "required_by_mboost_prediction_state",
      baseline_hazard = baseline_hazard
    ),
    class = c(
      "stormbreaker_mboost_cox_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_mboost_cox_model <- function(state, learner) {
  if (!is.raw(state$engine_raw) || length(state$engine_raw) == 0L ||
    !is.list(state$params)) {
    sb_bad_learner_state("The saved mboost Cox learner state is invalid.")
  }
  sb_require_mboost_cox_dependency(learner)
  engine <- tryCatch(
    unserialize(state$engine_raw),
    error = function(error) {
      sb_abort(
        paste0("The saved `", learner, "` model could not be loaded."),
        "stormbreaker_invalid_mboost_cox_state",
        parent = error
      )
    }
  )
  expected_class <- switch(learner,
    gamboost_cox = "gamboost",
    glmboost_cox = "glmboost",
    blackboost_cox = "blackboost"
  )
  if (!inherits(engine, expected_class)) {
    sb_bad_learner_state("The saved mboost Cox learner engine is invalid.")
  }
  sb_new_mboost_cox_model(
    learner,
    state$features,
    state$params,
    engine
  )
}

sb_predict_mboost_cox_state <- function(state, x) {
  learner <- if (is.list(state)) state$learner else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_mboost_cox_learner_ids()) {
    sb_bad_learner_state("The saved mboost Cox learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_restore_mboost_cox_model(state, learner)
  risk <- sb_predict_mboost_cox(model, x)
  result <- list(risk = risk)

  baseline <- state$baseline_hazard
  if (!is.null(baseline)) {
    valid_baseline <- is.list(baseline) &&
      is.numeric(baseline$prediction_horizon) &&
      length(baseline$prediction_horizon) == 1L &&
      !is.na(baseline$prediction_horizon) &&
      is.finite(baseline$prediction_horizon) &&
      baseline$prediction_horizon > 0
    if (!valid_baseline) {
      sb_bad_learner_state("The saved mboost Cox baseline state is invalid.")
    }
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      risk,
      baseline$cumulative_baseline_hazard
    )
  }
  result
}

sb_mboost_cox_learner_adapter <- function(
  learner, params = sb_mboost_cox_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_mboost_cox(x, time, event, learner, params)
    },
    predict = sb_predict_mboost_cox,
    predict_survival = sb_predict_mboost_cox_survival,
    extract_state = sb_extract_mboost_cox_state,
    predict_state = sb_predict_mboost_cox_state
  )
}

# ---- mboost-ranking-fixed-adapters ----

sb_mboost_ranking_learner_ids <- function() {
  c(
    "gamboost_gehan_aft",
    "gamboost_unoc",
    "glmboost_gehan_aft",
    "glmboost_unoc",
    "blackboost_gehan_aft",
    "blackboost_unoc"
  )
}

sb_mboost_ranking_params <- function(learner) {
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_mboost_ranking_learner_ids()) {
    sb_abort(
      "The requested fixed mboost ranking learner is not supported.",
      "stormbreaker_unsupported_learner"
    )
  }

  objective <- if (grepl("_gehan_", learner, fixed = TRUE)) {
    list(
      family = "Gehan",
      objective = "gehan_rank_aft",
      probability_method = "kaplan_meier_log_time_residual"
    )
  } else {
    list(
      family = "Cindex",
      objective = "uno_ipcw_concordance",
      sigma = 0.1,
      ipcw = 1,
      probability_method = "analysis_cox_breslow_calibration",
      calibration_ties = "breslow"
    )
  }
  common <- c(
    objective,
    list(
      mstop = 100L,
      nu = 0.1,
      risk = "inbag"
    )
  )
  if (startsWith(learner, "gamboost_")) {
    return(c(common, list(
      algorithm = "gamboost",
      base_learner = "bbs",
      dfbase = 4
    )))
  }
  if (startsWith(learner, "glmboost_")) {
    return(c(common, list(
      algorithm = "glmboost",
      base_learner = "bolscw",
      center = TRUE
    )))
  }
  c(common, list(
    algorithm = "blackboost",
    base_learner = "btree",
    tree_controls = list(
      teststat = "quad",
      testtype = "Teststatistic",
      mincriterion = 0,
      minsplit = 10L,
      minbucket = 4L,
      maxdepth = 2L,
      saveinfo = FALSE
    ),
    seed = 1L
  ))
}

sb_mboost_ranking_parameter_semantics <- function(learner) {
  params <- sb_mboost_ranking_params(learner)
  objective <- if (identical(params$family, "Gehan")) {
    paste0(
      "Gehan rank-based AFT loss with analysis-fold Kaplan-Meier ",
      "log-time residual survival"
    )
  } else {
    paste0(
      "Uno IPCW concordance loss with sigma 0.1 and fold-local IPCW, ",
      "plus analysis-fold Cox/Breslow horizon calibration"
    )
  }
  paste0(
    params$algorithm, " ", objective,
    ", 100 boosting iterations, shrinkage 0.1, and ",
    params$base_learner, " base learner fixed before resampling"
  )
}

sb_mboost_ranking_dependency <- function(learner) {
  params <- sb_mboost_ranking_params(learner)
  packages <- if (identical(params$algorithm, "blackboost")) {
    c("mboost", "partykit")
  } else {
    "mboost"
  }
  list(
    packages = packages,
    available = all(vapply(
      packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    ))
  )
}

sb_require_mboost_ranking_dependency <- function(learner) {
  dependency <- sb_mboost_ranking_dependency(learner)
  if (!dependency$available) {
    missing <- dependency$packages[!vapply(
      dependency$packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )]
    sb_abort(
      paste0(
        "Package", if (length(missing) > 1L) "s" else "", " `",
        paste(missing, collapse = "`, `"), "` ",
        if (length(missing) > 1L) "are" else "is",
        " required for learner `", learner, "`."
      ),
      "stormbreaker_missing_mboost_ranking_dependency"
    )
  }
  invisible(dependency)
}

sb_mboost_ranking_internal_feature_names <- function(n_features) {
  paste0(".stormbreaker_mboost_ranking_x", seq_len(n_features))
}

sb_mboost_ranking_formula_data <- function(x, time, event) {
  engine_features <- sb_mboost_ranking_internal_feature_names(ncol(x))
  engine_x <- x
  names(engine_x) <- engine_features
  formula_data <- sb_additional_survival_formula_data(
    engine_x,
    time,
    event
  )
  c(formula_data, list(engine_features = engine_features))
}

sb_mboost_ranking_family <- function(params) {
  if (identical(params$family, "Gehan")) {
    return(mboost::Gehan())
  }
  mboost::Cindex(sigma = params$sigma, ipcw = params$ipcw)
}

sb_fit_mboost_ranking_engine <- function(formula_data, params) {
  if (identical(params$algorithm, "gamboost")) {
    withr::local_options(mboost_useMatrix = FALSE)
  }
  control <- mboost::boost_control(
    mstop = params$mstop,
    nu = params$nu,
    risk = params$risk,
    center = TRUE,
    trace = FALSE
  )
  fit_engine <- function() {
    family <- sb_mboost_ranking_family(params)
    switch(params$algorithm,
      gamboost = mboost::gamboost(
        formula_data$formula,
        data = formula_data$data,
        family = family,
        control = control,
        baselearner = mboost::bbs,
        dfbase = params$dfbase
      ),
      glmboost = mboost::glmboost(
        formula_data$formula,
        data = formula_data$data,
        family = family,
        control = control,
        center = params$center
      ),
      blackboost = mboost::blackboost(
        formula_data$formula,
        data = formula_data$data,
        family = family,
        control = control,
        tree_controls = do.call(
          partykit::ctree_control,
          params$tree_controls
        )
      )
    )
  }
  fit <- function() {
    sb_classical_survival_stop_on_warning(
      fit_engine(),
      paste0("fixed ", params$algorithm, " ", params$objective, " learner"),
      "stormbreaker_mboost_ranking_fit_warning"
    )
  }
  if ("seed" %in% names(params)) {
    return(withr::with_seed(params$seed, fit()))
  }
  fit()
}

sb_new_mboost_ranking_model <- function(
  learner, features, engine_features, params, engine
) {
  structure(
    list(
      learner = learner,
      features = features,
      engine_features = engine_features,
      params = params,
      engine = engine
    ),
    class = c(
      paste0("stormbreaker_", learner),
      "stormbreaker_mboost_ranking_fixed"
    )
  )
}

sb_fit_mboost_ranking <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_mboost_ranking_dependency(learner)
  formula_data <- sb_mboost_ranking_formula_data(
    inputs$x,
    inputs$time,
    inputs$event
  )
  engine <- sb_fit_mboost_ranking_engine(formula_data, params)
  valid <- inherits(engine, params$algorithm) &&
    identical(as.integer(mboost::mstop(engine)), params$mstop)
  if (!valid) {
    sb_abort(
      paste0("Learner `", learner, "` did not return its requested fit."),
      "stormbreaker_incomplete_mboost_ranking_fit"
    )
  }

  model <- sb_new_mboost_ranking_model(
    learner,
    names(inputs$x),
    formula_data$engine_features,
    params,
    engine
  )
  sb_predict_mboost_ranking(model, inputs$x)
  model
}

sb_validate_mboost_ranking_model_x <- function(model, x) {
  valid_model <- inherits(model, "stormbreaker_mboost_ranking_fixed") &&
    is.character(model$learner) && length(model$learner) == 1L &&
    !is.na(model$learner) &&
    model$learner %in% sb_mboost_ranking_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.character(model$engine_features) &&
    length(model$engine_features) == length(model$features) &&
    !anyNA(model$engine_features) && !anyDuplicated(model$engine_features)
  if (!valid_model) {
    sb_abort(
      "The supplied model is not a fitted mboost ranking learner.",
      "stormbreaker_bad_mboost_ranking_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  sb_require_mboost_ranking_dependency(model$learner)
  x
}

sb_mboost_ranking_engine_x <- sb_mboost_aft_engine_x

sb_mboost_ranking_link <- function(model, x) {
  x <- sb_validate_mboost_ranking_model_x(model, x)
  link <- withCallingHandlers(
    as.numeric(stats::predict(
      model$engine,
      newdata = sb_mboost_ranking_engine_x(model, x),
      type = "link"
    )),
    warning = function(warning) {
      spline_extrapolation <- identical(
        model$params$algorithm,
        "gamboost"
      ) && identical(
        conditionMessage(warning),
        paste0(
          "Some 'x' values are beyond 'boundary.knots'; ",
          "Linear extrapolation used."
        )
      )
      if (spline_extrapolation) {
        invokeRestart("muffleWarning")
      }
      sb_abort(
        paste0(
          "The fixed ", model$learner, " prediction raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_mboost_ranking_prediction_warning",
        parent = warning
      )
    }
  )
  if (length(link) != nrow(x) || anyNA(link) || any(!is.finite(link))) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid rank scores."),
      "stormbreaker_invalid_mboost_ranking_prediction"
    )
  }
  unname(link)
}

sb_predict_mboost_ranking <- function(model, x) {
  -sb_mboost_ranking_link(model, x)
}

sb_fit_mboost_gehan_residual_survival <- function(time, event, link) {
  fit <- survival::survfit(
    survival::Surv(log(time) - link, event) ~ 1
  )
  curve <- list(
    time = as.numeric(fit$time),
    survival = as.numeric(fit$surv)
  )
  valid <- length(curve$time) > 0L &&
    length(curve$time) == length(curve$survival) &&
    !anyNA(curve$time) && all(is.finite(curve$time)) &&
    !is.unsorted(curve$time) &&
    !anyNA(curve$survival) && all(is.finite(curve$survival)) &&
    all(curve$survival >= 0 & curve$survival <= 1) &&
    all(diff(curve$survival) <= 0)
  if (!valid) {
    sb_abort(
      "The Gehan fit produced an invalid residual survival curve.",
      "stormbreaker_invalid_mboost_gehan_survival"
    )
  }
  curve
}

sb_predict_mboost_gehan_residual_survival <- function(
  curve, link, prediction_horizon
) {
  valid_curve <- is.list(curve) &&
    is.numeric(curve$time) && is.numeric(curve$survival) &&
    length(curve$time) > 0L &&
    length(curve$time) == length(curve$survival) &&
    !anyNA(curve$time) && all(is.finite(curve$time)) &&
    !is.unsorted(curve$time) &&
    !anyNA(curve$survival) && all(is.finite(curve$survival)) &&
    all(curve$survival >= 0 & curve$survival <= 1) &&
    all(diff(curve$survival) <= 0)
  if (!valid_curve) {
    sb_abort(
      "The saved Gehan residual survival curve is invalid.",
      "stormbreaker_invalid_mboost_gehan_survival"
    )
  }
  index <- findInterval(log(prediction_horizon) - link, curve$time)
  probability <- rep(1, length(index))
  use <- index > 0L
  probability[use] <- curve$survival[index[use]]
  if (anyNA(probability) || any(!is.finite(probability)) ||
    any(probability < 0 | probability > 1)) {
    sb_abort(
      "The Gehan fit produced invalid horizon survival probabilities.",
      "stormbreaker_invalid_mboost_gehan_survival"
    )
  }
  unname(probability)
}

sb_fit_mboost_unoc_calibration <- function(time, event, raw_risk, ties) {
  calibration <- sb_classical_survival_stop_on_warning(
    survival::coxph(
      survival::Surv(time, event) ~ raw_risk,
      data = data.frame(time = time, event = event, raw_risk = raw_risk),
      ties = ties,
      singular.ok = FALSE
    ),
    "mboost Uno-C Cox horizon calibration",
    "stormbreaker_mboost_unoc_calibration_warning"
  )
  coefficient <- unname(as.numeric(stats::coef(calibration)))
  if (length(coefficient) != 1L || is.na(coefficient) ||
    !is.finite(coefficient) || coefficient < 0) {
    sb_abort(
      paste0(
        "The mboost Uno-C horizon calibration did not preserve the ",
        "native higher-risk direction."
      ),
      "stormbreaker_invalid_mboost_unoc_calibration"
    )
  }
  coefficient
}

sb_mboost_ranking_direct_survival_source <- function(learner) {
  params <- sb_mboost_ranking_params(learner)
  if (!identical(params$family, "Gehan")) {
    return(NULL)
  }
  paste0("mboost_", params$algorithm, "_gehan_aft_residual_kaplan_meier")
}

sb_predict_mboost_ranking_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  inputs <- sb_validate_additional_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_link <- sb_mboost_ranking_link(model, inputs$x)
  assessment_link <- sb_mboost_ranking_link(model, x_assessment)

  if (identical(model$params$family, "Gehan")) {
    curve <- sb_fit_mboost_gehan_residual_survival(
      inputs$time,
      inputs$event,
      analysis_link
    )
    return(list(
      prediction_horizon = prediction_horizon,
      source = sb_mboost_ranking_direct_survival_source(model$learner),
      ties = "not_applicable",
      probability = sb_predict_mboost_gehan_residual_survival(
        curve,
        assessment_link,
        prediction_horizon
      ),
      residual_survival = curve
    ))
  }

  raw_analysis_risk <- -analysis_link
  raw_assessment_risk <- -assessment_link
  coefficient <- sb_fit_mboost_unoc_calibration(
    inputs$time,
    inputs$event,
    raw_analysis_risk,
    model$params$calibration_ties
  )
  calibrated_analysis_risk <- coefficient * raw_analysis_risk
  cumulative_hazard <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    calibrated_analysis_risk,
    prediction_horizon
  )
  result <- sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = paste0(model$learner, "_analysis_cox_breslow_calibration"),
    ties = model$params$calibration_ties,
    risk = coefficient * raw_assessment_risk
  )
  result$calibration_coefficient <- coefficient
  result
}

sb_extract_mboost_ranking_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_mboost_ranking_fixed")) {
    sb_bad_learner_state("The fitted mboost ranking learner state is invalid.")
  }

  residual_survival <- NULL
  calibration_coefficient <- NULL
  baseline_hazard <- NULL
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_mboost_ranking_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    if (identical(model$params$family, "Gehan")) {
      residual_survival <- probability$residual_survival
      survival_probability <- list(
        prediction_horizon = probability$prediction_horizon,
        time_unit = if ("time_unit" %in% names(context)) {
          context$time_unit
        } else {
          NA_character_
        },
        source = probability$source,
        ties = probability$ties
      )
    } else {
      calibration_coefficient <- probability$calibration_coefficient
      baseline_hazard <- list(
        prediction_horizon = probability$prediction_horizon,
        time_unit = if ("time_unit" %in% names(context)) {
          context$time_unit
        } else {
          NA_character_
        },
        cumulative_baseline_hazard = probability$cumulative_baseline_hazard,
        baseline_survival = probability$baseline_survival,
        source = probability$source,
        ties = probability$ties
      )
    }
  }

  structure(
    list(
      learner = model$learner,
      features = model$features,
      engine_features = model$engine_features,
      params = model$params,
      risk_reference = paste0(
        "mboost_", model$params$algorithm, "_",
        tolower(model$params$family), "_native"
      ),
      risk_formula = "negative_mboost_ranking_link",
      probability_method = model$params$probability_method,
      model_serialization = "R_xdr_v3_mboost",
      engine_raw = serialize(
        model$engine,
        connection = NULL,
        ascii = FALSE,
        xdr = TRUE,
        version = 3L
      ),
      training_data_retention = "required_by_mboost_prediction_state",
      residual_survival = residual_survival,
      calibration_coefficient = calibration_coefficient,
      baseline_hazard = baseline_hazard,
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_mboost_ranking_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_mboost_ranking_model <- function(state, learner) {
  valid <- is.raw(state$engine_raw) && length(state$engine_raw) > 0L &&
    is.list(state$params) &&
    is.character(state$engine_features) &&
    length(state$engine_features) == length(state$features) &&
    !anyNA(state$engine_features) && !anyDuplicated(state$engine_features)
  if (!valid) {
    sb_bad_learner_state("The saved mboost ranking learner state is invalid.")
  }
  sb_require_mboost_ranking_dependency(learner)
  engine <- tryCatch(
    unserialize(state$engine_raw),
    error = function(error) {
      sb_abort(
        paste0("The saved `", learner, "` model could not be loaded."),
        "stormbreaker_invalid_mboost_ranking_state",
        parent = error
      )
    }
  )
  if (!inherits(engine, state$params$algorithm) ||
    !sb_pipeline_parameters_equal(as.integer(mboost::mstop(engine)), state$params$mstop)) {
    sb_bad_learner_state("The saved mboost ranking learner engine is invalid.")
  }
  state$params$mstop <- as.integer(state$params$mstop)
  sb_new_mboost_ranking_model(
    learner,
    state$features,
    state$engine_features,
    state$params,
    engine
  )
}

sb_predict_mboost_ranking_state <- function(state, x) {
  learner <- if (is.list(state)) state$learner else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_mboost_ranking_learner_ids()) {
    sb_bad_learner_state("The saved mboost ranking learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_restore_mboost_ranking_model(state, learner)
  link <- sb_mboost_ranking_link(model, x)
  risk <- -link
  result <- list(risk = risk)

  if (identical(model$params$family, "Gehan")) {
    probability <- state$survival_probability
    if (!is.null(probability)) {
      horizon <- probability$prediction_horizon
      valid <- is.list(probability) && is.numeric(horizon) &&
        length(horizon) == 1L && !is.na(horizon) && is.finite(horizon) &&
        horizon > 0
      if (!valid) {
        sb_bad_learner_state("The saved mboost Gehan horizon state is invalid.")
      }
      result$prediction_horizon <- horizon
      result$time_unit <- probability$time_unit
      result$survival_probability <-
        sb_predict_mboost_gehan_residual_survival(
          state$residual_survival,
          link,
          horizon
        )
    }
    return(result)
  }

  baseline <- state$baseline_hazard
  if (!is.null(baseline)) {
    coefficient <- state$calibration_coefficient
    valid <- is.list(baseline) &&
      is.numeric(baseline$prediction_horizon) &&
      length(baseline$prediction_horizon) == 1L &&
      !is.na(baseline$prediction_horizon) &&
      is.finite(baseline$prediction_horizon) &&
      baseline$prediction_horizon > 0 &&
      is.numeric(coefficient) && length(coefficient) == 1L &&
      !is.na(coefficient) && is.finite(coefficient) && coefficient >= 0
    if (!valid) {
      sb_bad_learner_state("The saved mboost Uno-C horizon state is invalid.")
    }
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      coefficient * risk,
      baseline$cumulative_baseline_hazard
    )
  }
  result
}

sb_mboost_ranking_learner_adapter <- function(
  learner, params = sb_mboost_ranking_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_mboost_ranking(x, time, event, learner, params)
    },
    predict = sb_predict_mboost_ranking,
    predict_survival = sb_predict_mboost_ranking_survival,
    extract_state = sb_extract_mboost_ranking_state,
    predict_state = sb_predict_mboost_ranking_state
  )
}

# ---- model-gbm ----

# ---- algorithm-gbm-xgboost ----

sb_gbm_diagnostic_ids <- function() c("gbm_binary", "gbm_cox", "gbm_cox_cv")

sb_gbm_xgboost_diagnostic <- function(fit, type) {
  learner <- fit$pipeline$learner
  state <- fit$learner_state
  xgb <- learner %in% sb_xgboost_learner_ids()
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("training_loss", "cv_curve", "importance")) sb_abort(
    "This boosting engine supports training_loss, cv_curve or importance where originally available.",
    "stormbreaker_unsupported_diagnostic")
  if (xgb && type != "importance") sb_abort(
    "These XGBoost pipelines did not record an evaluation history; no loss or CV curve is reconstructed.",
    "stormbreaker_unavailable_diagnostic")
  selected <- NULL
  native <- NULL
  if (xgb) {
    model <- sb_load_xgboost_state_model(state, learner)
    native <- as.data.frame(xgboost::xgb.importance(feature_names = state$features, model = model$booster))
    linear <- identical(state$params$booster, "gblinear")
    columns <- if (linear) "Weight" else c("Gain", "Cover", "Frequency")
    if (nrow(native) && (!all(c("Feature", columns) %in% names(native)) ||
        anyDuplicated(native$Feature) || !all(native$Feature %in% state$features) ||
        any(!is.finite(as.matrix(native[columns]))))) sb_abort(
      "Native XGBoost feature scores are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(feature = state$features)
    for (column in columns) {
      values <- rep(0, length(state$features))
      if (nrow(native)) values[match(native$Feature, state$features)] <- native[[column]]
      tab[[column]] <- values
    }
    tab$value <- if (linear) tab$Weight else tab$Gain
    measure <- if (linear) "Original linear coefficient" else "Native split-gain share"
    scope <- if (linear) {
      scale <- if (learner %in% sb_xgboost_aft_learner_ids()) "log survival time" else
        if (learner %in% sb_xgboost_cox_learner_ids()) "log hazard" else "log odds for class 1"
      paste0("Coefficient scale: ", scale, ". Bias/base score excluded; not a complete prediction equation.")
    } else "Native gain, cover and split-frequency shares retained; not permutation importance or SHAP."
  } else if (type == "cv_curve") {
    if (!identical(learner, "gbm_cox_cv")) sb_abort(
      "This fixed GBM pipeline did not run inner CV.", "stormbreaker_unavailable_diagnostic")
    if (identical(state$tuning$search_status, "single_candidate_no_search")) sb_abort(
      "Only one tree count was allowed; no tree-count search or CV curve was computed.",
      "stormbreaker_unavailable_diagnostic")
    loss <- state$tuning$cv_error
    selected <- state$tuning$selected_trees
    if (!is.numeric(loss) || !length(loss) || any(!is.finite(loss)) ||
        length(selected) != 1L || !selected %in% seq_along(loss)) sb_abort(
      "The original GBM tuning sequence is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(iteration = seq_along(loss), loss = loss)
    measure <- "Native Cox inner-CV error"
    scope <- "Original inner-CV search, not validation of prior feature selection or untouched-cohort performance."
  } else {
    if (!requireNamespace("gbm", quietly = TRUE)) sb_abort("Package gbm is required.",
      "stormbreaker_missing_dependency")
    if (learner == "gbm_binary") {
      if (type == "training_loss") {
        if (is.null(fit$diagnostic_artifacts)) sb_abort(
          "Binary GBM training errors require capture_diagnostics when fitting.",
          "stormbreaker_missing_diagnostic_artifacts")
        engine <- fit$diagnostic_artifacts$learner$engine
      } else engine <- unserialize(state$engine_raw)
    } else engine <- unserialize(state$model_raw)
    if (!inherits(engine, "gbm")) sb_abort("The native GBM model is missing.",
      "stormbreaker_invalid_diagnostic")
    if (type == "training_loss") {
      loss <- engine$train.error
      if (!is.numeric(loss) || length(loss) != engine$n.trees || !length(loss) || any(!is.finite(loss))) sb_abort(
        "The original GBM training errors are incomplete or invalid.", "stormbreaker_invalid_diagnostic")
      tab <- data.frame(iteration = seq_along(loss), loss = loss)
      measure <- paste("Native training error:", engine$distribution$name)
      scope <- "Training sequence of the final fit; not the separate tuning search or held-out performance."
    } else {
      value <- gbm::relative.influence(engine, n.trees = engine$n.trees, scale. = FALSE, sort. = FALSE)
      if (length(value) != length(engine$var.names) || any(!is.finite(value)) ||
          !setequal(engine$var.names, state$features)) sb_abort(
        "Native GBM relative influence is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
      native <- data.frame(feature = engine$var.names, value = unname(value))
      tab <- native[match(state$features, native$feature), , drop = FALSE]
      rownames(tab) <- NULL
      measure <- "Native relative influence (unscaled)"
      scope <- "Original split-based influence; not permutation importance, a causal effect or SHAP."
    }
  }
  structure(list(type = paste0("boosted_tree_", type), table = tab, native_table = native,
    learner = learner, pipeline_id = fit$pipeline$id, selected_iteration = selected,
    measure = measure, interpretation = scope, source = "original fitted boosting model; no refitting"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_gbm_xgboost_diagnostic <- function(x) {
  value <- feature <- iteration <- loss <- NULL
  if (x$type == "boosted_tree_importance") {
    tab <- x$table
    tab$feature <- factor(tab$feature, levels = tab$feature[order(abs(tab$value))])
    p <- ggplot2::ggplot(tab, ggplot2::aes(value, feature)) +
      ggplot2::geom_col(fill = "#2878A0", width = .65) +
      ggplot2::geom_vline(xintercept = 0, colour = "#909090") +
      ggplot2::labs(title = x$measure, subtitle = x$learner, x = x$measure, y = NULL,
        caption = paste0(x$interpretation, "\nAll input features retained; zero-use features are not silently dropped.")) + sb_viz_theme()
    attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = max(5, 2.5 + .28 * nrow(tab)))
    return(p)
  }
  p <- ggplot2::ggplot(x$table, ggplot2::aes(iteration, loss)) +
    ggplot2::geom_point(colour = "#2878A0", size = 1) +
    ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
    ggplot2::labs(title = x$measure, subtitle = x$learner, x = "Number of trees", y = "Native error",
      caption = paste0(x$interpretation, "\nNo refitting, error-bar invention or alteration of the selected tree count.")) + sb_viz_theme()
  if (nrow(x$table) > 1L) p <- p + ggplot2::geom_line(colour = "#2878A0")
  if (!is.null(x$selected_iteration)) p <- p + ggplot2::geom_vline(
    xintercept = x$selected_iteration, linetype = 2, colour = "#D77935")
  p
}

# ---- gbm-cox-fixed-adapter ----

sb_gbm_available <- function() {
  requireNamespace("gbm", quietly = TRUE)
}

sb_validate_gbm_shrinkage <- function(shrinkage) {
  valid <- is.numeric(shrinkage) && length(shrinkage) == 1L &&
    !is.na(shrinkage) && is.finite(shrinkage) && shrinkage > 0
  if (!valid) {
    sb_abort(
      "`gbm_shrinkage` must be one finite positive number.",
      "stormbreaker_bad_gbm_shrinkage"
    )
  }
  as.numeric(shrinkage)
}

sb_validate_gbm_bag_fraction <- function(bag_fraction) {
  valid <- is.numeric(bag_fraction) && length(bag_fraction) == 1L &&
    !is.na(bag_fraction) && is.finite(bag_fraction) &&
    bag_fraction > 0 && bag_fraction <= 1
  if (!valid) {
    sb_abort(
      "`gbm_bag_fraction` must be one finite number in (0, 1].",
      "stormbreaker_bad_gbm_bag_fraction"
    )
  }
  as.numeric(bag_fraction)
}

sb_gbm_cox_params <- function(
    n_trees = 100L,
    interaction_depth = 1L,
    min_obs_node = 10L,
    shrinkage = 0.1,
    bag_fraction = 1,
    seed = 1L) {
  list(
    n_trees = sb_check_positive_integer(
      n_trees,
      "gbm_n_trees",
      minimum = 1L
    ),
    interaction_depth = sb_check_positive_integer(
      interaction_depth,
      "gbm_interaction_depth",
      minimum = 1L
    ),
    min_obs_node = sb_check_positive_integer(
      min_obs_node,
      "gbm_min_obs_node",
      minimum = 1L
    ),
    shrinkage = sb_validate_gbm_shrinkage(shrinkage),
    bag_fraction = sb_validate_gbm_bag_fraction(bag_fraction),
    seed = sb_check_positive_integer(seed, "gbm_seed", minimum = 1L),
    distribution = "coxph",
    keep_data = FALSE
  )
}

sb_gbm_cox_cv_params <- function(inner_folds = 10L, seed = 1L,
                                 maximum_trees = 10000L) {
  list(
    inner_folds = sb_check_positive_integer(
      inner_folds,
      "learner_tuning_folds",
      minimum = 3L
    ),
    seed = sb_check_positive_integer(seed, "learner_seed", minimum = 1L),
    maximum_trees = sb_check_positive_integer(
      maximum_trees,
      "gbm_cv_max_trees",
      minimum = 1L
    ),
    interaction_depth = 3L,
    min_obs_node = 10L,
    shrinkage = 0.001,
    bag_fraction = 0.5,
    distribution = "coxph",
    keep_data = FALSE,
    fit_interface = "formula"
  )
}

sb_fit_gbm_cox <- function(x, time, event, params) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time < 0)) {
    sb_abort(
      "`time` must contain non-negative values.",
      "stormbreaker_bad_adapter_time"
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  if (!any(event == 1L)) {
    sb_abort(
      "A fixed GBM-Cox analysis fold must contain at least one event.",
      "stormbreaker_gbm_cox_no_event"
    )
  }
  if (!sb_gbm_available()) {
    sb_abort(
      "Package `gbm` is required for the fixed GBM-Cox learner.",
      "stormbreaker_missing_gbm"
    )
  }

  booster <- withr::with_seed(
    params$seed,
    gbm::gbm.fit(
      x = as.matrix(x),
      y = survival::Surv(time, event),
      distribution = params$distribution,
      n.trees = params$n_trees,
      interaction.depth = params$interaction_depth,
      n.minobsinnode = params$min_obs_node,
      shrinkage = params$shrinkage,
      bag.fraction = params$bag_fraction,
      nTrain = nrow(x),
      keep.data = params$keep_data,
      verbose = FALSE,
      var.names = names(x),
      response.name = "Surv(time, event)"
    )
  )
  if (!inherits(booster, "gbm") ||
      !identical(as.integer(booster$n.trees), params$n_trees)) {
    sb_abort(
      "Fixed GBM-Cox did not fit the declared number of trees.",
      "stormbreaker_incomplete_gbm_fit"
    )
  }

  structure(
    list(
      features = names(x),
      params = params,
      n_trees = params$n_trees,
      booster = booster
    ),
    class = "stormbreaker_gbm_cox"
  )
}

sb_fit_gbm_cox_cv <- function(x, time, event, params) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time < 0)) {
    sb_abort(
      "`time` must contain non-negative values.",
      "stormbreaker_bad_adapter_time"
    )
  }
  if (!all(event %in% c(0, 1)) || !any(event == 1L)) {
    sb_abort(
      "A CV-tuned GBM-Cox analysis fold requires 0/1 events and an event.",
      "stormbreaker_gbm_cox_no_event"
    )
  }
  if (!sb_gbm_available()) {
    sb_abort(
      "Package `gbm` is required for the CV-tuned GBM-Cox learner.",
      "stormbreaker_missing_gbm"
    )
  }
  formula_data <- sb_additional_survival_formula_data(x, time, event)
  fit_gbm <- function(n_trees, cv_folds = params$inner_folds, search = NULL) {
    native <- gbm::gbm
    if (!is.null(search)) {
      reuse_cv <- function(...) {
        args <- list(...)
        ns <- asNamespace("gbm")
        i.train <- seq_len(args$nTrain)
        cv.group <- get("getCVgroup", ns)(args$distribution,
          args$class.stratify.cv, args$y, i.train, args$cv.folds, args$group)
        seeds <- as.integer(stats::runif(args$cv.folds, -(2^31 - 1), 2^31))
        # Serial native CV resets the seed for each fold. Its last fold's
        # terminal RNG is needed by the subsequent full-analysis booster.
        do_fold <- get("gbmDoFold", ns)
        last <- args[intersect(names(args), names(formals(do_fold)))]
        last$X <- args$cv.folds
        last$i.train <- i.train
        last$cv.group <- cv.group
        last$s <- seeds
        invisible(do.call(do_fold, last))
        list(error = search$cv.error[seq_len(args$n.trees)],
          predictions = search$cv.fitted)
      }
      environment(native) <- list2env(list(gbmCrossVal = reuse_cv),
        parent = environment(native))
    }
    booster <- NULL
    invisible(utils::capture.output(
      booster <- withr::with_seed(
        params$seed,
        native(
          formula = formula_data$formula,
          data = formula_data$data,
          distribution = params$distribution,
          n.trees = n_trees,
          interaction.depth = params$interaction_depth,
          n.minobsinnode = params$min_obs_node,
          shrinkage = params$shrinkage,
          bag.fraction = params$bag_fraction,
          cv.folds = cv_folds,
          n.cores = 1L,
          keep.data = params$keep_data,
          verbose = FALSE
        )
      )
    ))
    # Keep the same recorded native call for full and reused CV execution.
    booster$call[[1L]] <- quote(gbm::gbm)
    booster$m[[1L]] <- quote(gbm::gbm)
    booster
  }
  single_candidate <- params$maximum_trees == 1L
  # There is no tree-count choice with one candidate. Native gbm's one-row
  # CV reduction also fails; fit that candidate and record that no search ran.
  search_warned <- FALSE
  search <- if (single_candidate) NULL else withCallingHandlers(
    fit_gbm(params$maximum_trees),
    warning = function(w) search_warned <<- TRUE
  )
  best_trees <- if (single_candidate) 1L else as.integer(which.min(search$cv.error))
  if (length(best_trees) != 1L || is.na(best_trees) || best_trees < 1L ||
    best_trees > params$maximum_trees) {
    sb_abort(
      "CV-tuned GBM-Cox returned an invalid tree count.",
      "stormbreaker_invalid_gbm_cv"
    )
  }
  # Native gbm CV drops its error matrix when n.trees is one. The tree count
  # is already settled above; the final one-tree fit needs no further CV.
  booster <- if (best_trees == 1L) {
    fit_gbm(best_trees, cv_folds = 0L)
  } else if (search_warned) {
    fit_gbm(best_trees)
  } else if (best_trees == params$maximum_trees) {
    search
  } else {
    fit_gbm(best_trees, search = search)
  }
  if (!inherits(booster, "gbm") ||
    !identical(as.integer(booster$n.trees), best_trees)) {
    sb_abort(
      "CV-tuned GBM-Cox did not fit the selected number of trees.",
      "stormbreaker_incomplete_gbm_fit"
    )
  }
  tuning <- list(selected_trees = best_trees, cv_error = as.numeric(search$cv.error))
  if (single_candidate) tuning$search_status <- "single_candidate_no_search"
  structure(
    list(
      learner = "gbm_cox_cv",
      features = names(x),
      params = params,
      n_trees = best_trees,
      tuning = tuning,
      booster = booster
    ),
    class = "stormbreaker_gbm_cox"
  )
}

sb_predict_gbm_cox <- function(model, x) {
  if (!inherits(model, "stormbreaker_gbm_cox")) {
    sb_abort(
      "The supplied model is not a fitted fixed GBM-Cox learner.",
      "stormbreaker_bad_gbm_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed GBM-Cox schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (!sb_gbm_available()) {
    sb_abort(
      "Package `gbm` is required for fixed GBM-Cox prediction.",
      "stormbreaker_missing_gbm"
    )
  }

  risk <- as.numeric(stats::predict(
    model$booster,
    newdata = if (identical(model$params$fit_interface, "formula")) {
      x
    } else {
      as.matrix(x)
    },
    n.trees = model$n_trees,
    type = "link"
  ))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "Fixed GBM-Cox produced invalid log-hazard predictions.",
      "stormbreaker_invalid_gbm_prediction"
    )
  }
  unname(risk)
}

sb_predict_gbm_cox_survival <- function(
    model, x_analysis, time, event, x_assessment, prediction_horizon) {
  x_analysis <- sb_validate_predictor_frame(x_analysis, "x_analysis")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x_analysis))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x_analysis))
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    time
  )
  analysis_risk <- sb_predict_gbm_cox(model, x_analysis)
  assessment_risk <- sb_predict_gbm_cox(model, x_assessment)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    time,
    as.integer(event),
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = paste0(
      if ("learner" %in% names(model)) model$learner else "gbm_cox",
      "_breslow_log_hazard"
    ),
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_extract_gbm_cox_tree_model <- function(model) {
  if (!inherits(model, "stormbreaker_gbm_cox")) {
    sb_abort(
      "The supplied model is not a fitted fixed GBM-Cox learner.",
      "stormbreaker_bad_gbm_fit"
    )
  }
  if (!sb_gbm_available()) {
    sb_abort(
      "Package `gbm` is required to audit a GBM-Cox model.",
      "stormbreaker_missing_gbm"
    )
  }

  importance <- as.numeric(gbm::relative.influence(
    model$booster,
    n.trees = model$n_trees,
    scale. = FALSE,
    sort. = FALSE
  ))
  if (length(importance) != length(model$features) ||
      anyNA(importance) || any(!is.finite(importance))) {
    sb_abort(
      "Fixed GBM-Cox did not return one finite importance per feature.",
      "stormbreaker_invalid_gbm_importance"
    )
  }

  list(
    n_trees = model$n_trees,
    features = model$features,
    params = model$params,
    tuning = if ("tuning" %in% names(model)) model$tuning else NULL,
    importance_type = "relative_influence",
    model_serialization = "R_xdr_v3_gbm",
    model_raw = serialize(
      model$booster,
      connection = NULL,
      ascii = FALSE,
      xdr = TRUE,
      version = 3L
    ),
    feature_importance = data.frame(
      feature = model$features,
      importance = importance,
      stringsAsFactors = FALSE
    )
  )
}

sb_gbm_cox_learner_adapter <- function(
    n_trees = 100L,
    interaction_depth = 1L,
    min_obs_node = 10L,
    shrinkage = 0.1,
    bag_fraction = 1,
    seed = 1L) {
  params <- sb_gbm_cox_params(
    n_trees = n_trees,
    interaction_depth = interaction_depth,
    min_obs_node = min_obs_node,
    shrinkage = shrinkage,
    bag_fraction = bag_fraction,
    seed = seed
  )
  list(
    id = "gbm_cox",
    params = params,
    fit = function(x, time, event) {
      sb_fit_gbm_cox(x, time, event, params)
    },
    predict = sb_predict_gbm_cox,
    predict_survival = sb_predict_gbm_cox_survival,
    extract_tree_model = sb_extract_gbm_cox_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_gbm_cox_learner_state(model, context)
    },
    predict_state = sb_predict_gbm_cox_learner_state
  )
}

sb_gbm_cox_cv_learner_adapter <- function(inner_folds = 10L,
                                           seed = 1L,
                                           maximum_trees = 10000L) {
  params <- sb_gbm_cox_cv_params(
    inner_folds = inner_folds,
    seed = seed,
    maximum_trees = maximum_trees
  )
  list(
    id = "gbm_cox_cv",
    params = params,
    fit = function(x, time, event) {
      sb_fit_gbm_cox_cv(x, time, event, params)
    },
    predict = sb_predict_gbm_cox,
    predict_survival = sb_predict_gbm_cox_survival,
    extract_tree_model = sb_extract_gbm_cox_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_gbm_cox_learner_state(model, context)
    },
    predict_state = sb_predict_gbm_cox_learner_state
  )
}

# ---- model-xgboost ----

sb_xgboost_binary_learner_ids <- function() {
  c(
    "xgboost_binary",
    "xgboost_binary_dart",
    "xgboost_binary_linear"
  )
}

sb_xgboost_cox_learner_ids <- function() {
  c(
    "xgboost_cox",
    "xgboost_cox_dart",
    "xgboost_cox_linear"
  )
}

sb_xgboost_aft_learner_ids <- function() {
  c(
    "xgboost_aft",
    "xgboost_aft_dart",
    "xgboost_aft_linear"
  )
}

sb_xgboost_survival_learner_ids <- function() {
  c(sb_xgboost_cox_learner_ids(), sb_xgboost_aft_learner_ids())
}

sb_xgboost_learner_ids <- function() {
  c(sb_xgboost_binary_learner_ids(), sb_xgboost_survival_learner_ids())
}

sb_xgboost_booster <- function(learner) {
  if (learner %in% c(
    "xgboost_binary_dart",
    "xgboost_cox_dart",
    "xgboost_aft_dart"
  )) {
    return("dart")
  }
  if (learner %in% c(
    "xgboost_binary_linear",
    "xgboost_cox_linear",
    "xgboost_aft_linear"
  )) {
    return("gblinear")
  }
  "gbtree"
}

sb_xgboost_available <- function() {
  requireNamespace("xgboost", quietly = TRUE)
}

sb_xgboost_selector_gain <- function(model) {
  importance <- xgboost::xgb.importance(feature_names = model$features, model = model$booster)
  score <- numeric(length(model$features))
  if (nrow(importance)) score[match(importance$Feature, model$features)] <- importance$Gain
  score
}

sb_xgboost_params <- function(learner) {
  if (!learner %in% sb_xgboost_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed XGBoost learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  objective <- if (learner %in% sb_xgboost_binary_learner_ids()) {
    "binary:logistic"
  } else if (learner %in% sb_xgboost_cox_learner_ids()) {
    "survival:cox"
  } else {
    "survival:aft"
  }
  booster <- sb_xgboost_booster(learner)
  if (identical(booster, "gblinear")) {
    params <- list(
      nrounds = 100L,
      booster = booster,
      eta = 0.1,
      lambda = 1,
      alpha = 0,
      updater = "coord_descent",
      feature_selector = "cyclic",
      nthread = 1L,
      seed = 1L,
      objective = objective
    )
    if (learner %in% sb_xgboost_aft_learner_ids()) {
      params <- c(params, list(
        aft_loss_distribution = "normal",
        aft_loss_distribution_scale = 1
      ))
    }
    return(params)
  }
  tree_params <- list(
    nrounds = 100L,
    eta = 0.1,
    max_depth = 3L,
    min_child_weight = 1,
    subsample = 1,
    colsample_bytree = 1
  )
  if (identical(booster, "dart")) {
    params <- c(
      tree_params["nrounds"],
      list(booster = booster),
      tree_params[-1L],
      list(
        sample_type = "uniform",
        normalize_type = "tree",
        rate_drop = 0.1,
        skip_drop = 0,
        one_drop = FALSE,
        nthread = 1L,
        seed = 1L,
        objective = objective
      )
    )
    if (learner %in% sb_xgboost_aft_learner_ids()) {
      params <- c(params, list(
        aft_loss_distribution = "normal",
        aft_loss_distribution_scale = 1
      ))
    }
    return(params)
  }
  params <- c(
    tree_params,
    list(
      nthread = 1L,
      seed = 1L,
      objective = objective
    )
  )
  if (learner %in% sb_xgboost_aft_learner_ids()) {
    params <- c(params, list(
      aft_loss_distribution = "normal",
      aft_loss_distribution_scale = 1
    ))
  }
  params
}

sb_validate_xgboost_params <- function(learner, params) {
  expected <- names(sb_xgboost_params(learner))
  valid <- is.list(params) && identical(names(params), expected) &&
    is.numeric(params$nrounds) && length(params$nrounds) == 1L &&
    !is.na(params$nrounds) && params$nrounds >= 1 &&
    identical(params$nrounds, as.integer(params$nrounds)) &&
    is.numeric(params$nthread) && length(params$nthread) == 1L &&
    !is.na(params$nthread) && params$nthread >= 1 &&
    identical(params$nthread, as.integer(params$nthread)) &&
    is.numeric(params$seed) && length(params$seed) == 1L &&
    !is.na(params$seed) && params$seed >= 1 &&
    identical(params$seed, as.integer(params$seed)) &&
    is.character(params$objective) && length(params$objective) == 1L &&
    identical(
      params$objective,
      sb_xgboost_params(learner)$objective
    )
  booster <- sb_xgboost_booster(learner)
  if (booster %in% c("gbtree", "dart")) {
    numeric_controls <- unlist(params[c(
      "eta", "max_depth", "min_child_weight", "subsample",
      "colsample_bytree"
    )], use.names = FALSE)
    valid <- valid && !anyNA(numeric_controls) &&
      all(is.finite(numeric_controls)) && params$eta > 0 &&
      params$max_depth >= 1 && params$min_child_weight >= 0 &&
      params$subsample > 0 && params$subsample <= 1 &&
      params$colsample_bytree > 0 && params$colsample_bytree <= 1
  }
  if (identical(booster, "dart")) {
    dropout_controls <- unlist(
      params[c("rate_drop", "skip_drop")],
      use.names = FALSE
    )
    valid <- valid && identical(params$booster, booster) &&
      !anyNA(dropout_controls) && all(is.finite(dropout_controls)) &&
      all(dropout_controls >= 0 & dropout_controls <= 1) &&
      identical(params$sample_type, "uniform") &&
      identical(params$normalize_type, "tree") &&
      is.logical(params$one_drop) && length(params$one_drop) == 1L &&
      !is.na(params$one_drop)
  }
  if (identical(booster, "gblinear")) {
    linear_controls <- unlist(
      params[c("eta", "lambda", "alpha")],
      use.names = FALSE
    )
    valid <- valid && identical(params$booster, booster) &&
      !anyNA(linear_controls) && all(is.finite(linear_controls)) &&
      params$eta > 0 && params$lambda >= 0 && params$alpha >= 0 &&
      identical(params$updater, "coord_descent") &&
      identical(params$feature_selector, "cyclic")
  }
  if (learner %in% sb_xgboost_aft_learner_ids()) {
    valid <- valid &&
      identical(params$aft_loss_distribution, "normal") &&
      is.numeric(params$aft_loss_distribution_scale) &&
      length(params$aft_loss_distribution_scale) == 1L &&
      !is.na(params$aft_loss_distribution_scale) &&
      is.finite(params$aft_loss_distribution_scale) &&
      params$aft_loss_distribution_scale > 0
  }
  if (!valid) {
    sb_abort(
      "The fixed XGBoost parameter specification is invalid.",
      "stormbreaker_bad_xgboost_parameters"
    )
  }
  params$nrounds <- as.integer(params$nrounds)
  if (booster %in% c("gbtree", "dart")) {
    params$max_depth <- as.integer(params$max_depth)
  }
  params$nthread <- as.integer(params$nthread)
  params$seed <- as.integer(params$seed)
  params
}

sb_require_xgboost <- function() {
  if (!sb_xgboost_available()) {
    sb_abort(
      "Package `xgboost` is required for a fixed XGBoost learner.",
      "stormbreaker_missing_xgboost"
    )
  }
  invisible(TRUE)
}

sb_validate_xgboost_survival_inputs <- function(
  x, time, event, learner = "xgboost_cox"
) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  is_aft <- learner %in% sb_xgboost_aft_learner_ids()
  if (any(time <= 0)) {
    sb_abort(
      paste0(
        "A fixed XGBoost-", if (is_aft) "AFT" else "Cox",
        " analysis fold requires strictly positive times."
      ),
      if (is_aft) {
        "stormbreaker_xgboost_aft_nonpositive_time"
      } else {
        "stormbreaker_xgboost_cox_nonpositive_time"
      }
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  if (!any(event == 1L)) {
    sb_abort(
      paste0(
        "A fixed XGBoost-", if (is_aft) "AFT" else "Cox",
        " analysis fold must contain at least one event."
      ),
      if (is_aft) {
        "stormbreaker_xgboost_aft_no_event"
      } else {
        "stormbreaker_xgboost_cox_no_event"
      }
    )
  }
  list(x = x, time = time, event = as.integer(event))
}

sb_fit_xgboost <- function(x, endpoint, learner, params) {
  params <- sb_validate_xgboost_params(learner, params)
  if (learner %in% sb_xgboost_binary_learner_ids()) {
    x <- sb_validate_predictor_frame(x, "x")
    endpoint <- sb_validate_binary_outcome_vector(endpoint, nrow(x))
    if (!all(c(0L, 1L) %in% endpoint)) {
      sb_abort(
        "A fixed XGBoost binary analysis fold must contain both classes.",
        "stormbreaker_xgboost_binary_single_class"
      )
    }
    label <- as.numeric(endpoint)
  } else {
    inputs <- sb_validate_xgboost_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      learner
    )
    x <- inputs$x
    label <- if (learner %in% sb_xgboost_cox_learner_ids()) {
      ifelse(inputs$event == 1L, inputs$time, -inputs$time)
    } else {
      NULL
    }
  }
  sb_require_xgboost()

  dtrain <- if (learner %in% sb_xgboost_aft_learner_ids()) {
    matrix <- xgboost::xgb.DMatrix(
      as.matrix(x),
      nthread = params$nthread
    )
    xgboost::setinfo(matrix, "label_lower_bound", inputs$time)
    xgboost::setinfo(
      matrix,
      "label_upper_bound",
      ifelse(inputs$event == 1L, inputs$time, Inf)
    )
    matrix
  } else {
    xgboost::xgb.DMatrix(
      as.matrix(x),
      label = label,
      nthread = params$nthread
    )
  }
  booster_type <- sb_xgboost_booster(learner)
  engine_params <- if (identical(booster_type, "gblinear")) {
    list(
      objective = params$objective,
      booster = booster_type,
      eta = params$eta,
      lambda = params$lambda,
      alpha = params$alpha,
      updater = params$updater,
      feature_selector = params$feature_selector,
      nthread = params$nthread
    )
  } else {
    tree_params <- list(
      objective = params$objective,
      booster = booster_type,
      eta = params$eta,
      max_depth = params$max_depth,
      min_child_weight = params$min_child_weight,
      subsample = params$subsample,
      colsample_bytree = params$colsample_bytree,
      nthread = params$nthread
    )
    if (identical(booster_type, "dart")) {
      tree_params$sample_type <- params$sample_type
      tree_params$normalize_type <- params$normalize_type
      tree_params$rate_drop <- params$rate_drop
      tree_params$skip_drop <- params$skip_drop
      tree_params$one_drop <- params$one_drop
    }
    tree_params
  }
  if (learner %in% sb_xgboost_aft_learner_ids()) {
    engine_params$aft_loss_distribution <- params$aft_loss_distribution
    engine_params$aft_loss_distribution_scale <-
      params$aft_loss_distribution_scale
  }
  booster <- withr::with_seed(
    params$seed,
    xgboost::xgb.train(
      params = engine_params,
      data = dtrain,
      nrounds = params$nrounds,
      verbose = 0
    )
  )
  if (!inherits(booster, "xgb.Booster")) {
    sb_abort(
      "Fixed XGBoost did not return a usable booster.",
      "stormbreaker_incomplete_xgboost_fit"
    )
  }
  structure(
    list(
      learner = learner,
      features = names(x),
      params = params,
      booster = booster
    ),
    class = "stormbreaker_xgboost_fixed"
  )
}

sb_xgboost_prediction <- function(model, x) {
  if (!inherits(model, "stormbreaker_xgboost_fixed")) {
    sb_abort(
      "The supplied model is not a fitted fixed XGBoost learner.",
      "stormbreaker_bad_xgboost_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed XGBoost schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  sb_require_xgboost()
  dnew <- xgboost::xgb.DMatrix(
    as.matrix(x),
    nthread = model$params$nthread
  )
  value <- as.numeric(stats::predict(
    model$booster,
    dnew,
    outputmargin = model$learner %in% sb_xgboost_survival_learner_ids()
  ))
  if (model$learner %in% sb_xgboost_aft_learner_ids()) {
    value <- -value
  }
  valid <- length(value) == nrow(x) && !anyNA(value) && all(is.finite(value))
  if (model$learner %in% sb_xgboost_binary_learner_ids()) {
    valid <- valid && all(value >= 0 & value <= 1)
  }
  if (!valid) {
    sb_abort(
      "Fixed XGBoost produced invalid predictions.",
      "stormbreaker_invalid_xgboost_prediction"
    )
  }
  unname(value)
}

sb_predict_xgboost_cox_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  if (!model$learner %in% sb_xgboost_cox_learner_ids()) {
    sb_abort(
      "Survival probability requires the fixed XGBoost-Cox learner.",
      "stormbreaker_bad_xgboost_fit"
    )
  }
  inputs <- sb_validate_xgboost_survival_inputs(
    x_analysis, time, event, model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fixed XGBoost-Cox schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_risk <- sb_xgboost_prediction(model, inputs$x)
  assessment_risk <- sb_xgboost_prediction(model, x_assessment)
  baseline <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    baseline,
    source = "xgboost_cox_breslow_calibration",
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_xgboost_aft_survival_probability <- function(
  risk, prediction_horizon, scale
) {
  valid <- is.numeric(risk) && length(risk) > 0L && !anyNA(risk) &&
    all(is.finite(risk)) && is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L && !is.na(prediction_horizon) &&
    is.finite(prediction_horizon) && prediction_horizon > 0 &&
    is.numeric(scale) && length(scale) == 1L && !is.na(scale) &&
    is.finite(scale) && scale > 0
  if (!valid) {
    sb_abort(
      "Fixed XGBoost-AFT received invalid prediction inputs.",
      "stormbreaker_invalid_xgboost_aft_prediction"
    )
  }
  margin <- -risk
  probability <- stats::pnorm(
    (margin - log(prediction_horizon)) / scale
  )
  if (length(probability) != length(risk) || anyNA(probability) ||
    any(!is.finite(probability)) ||
    any(probability < 0 | probability > 1)) {
    sb_abort(
      "Fixed XGBoost-AFT produced invalid survival probabilities.",
      "stormbreaker_invalid_xgboost_aft_prediction"
    )
  }
  unname(probability)
}

sb_predict_xgboost_aft_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  if (!model$learner %in% sb_xgboost_aft_learner_ids()) {
    sb_abort(
      "Survival probability requires a fixed XGBoost-AFT learner.",
      "stormbreaker_bad_xgboost_fit"
    )
  }
  inputs <- sb_validate_xgboost_survival_inputs(
    x_analysis, time, event, model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fixed XGBoost-AFT schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  risk <- sb_xgboost_prediction(model, x_assessment)
  list(
    prediction_horizon = prediction_horizon,
    source = "xgboost_aft_normal_native",
    ties = "not_applicable",
    probability = sb_xgboost_aft_survival_probability(
      risk,
      prediction_horizon,
      model$params$aft_loss_distribution_scale
    )
  )
}

sb_xgboost_model_raw <- function(model) {
  sb_require_xgboost()
  raw <- xgboost::xgb.save.raw(model$booster, raw_format = "json")
  if (!is.raw(raw) || length(raw) == 0L) {
    sb_bad_learner_state("The fitted XGBoost booster could not be saved.")
  }
  raw
}

sb_extract_xgboost_learner_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_xgboost_fixed") ||
    !model$learner %in% sb_xgboost_learner_ids()) {
    sb_bad_learner_state("The fitted XGBoost learner state is invalid.")
  }
  baseline_hazard <- NULL
  if (model$learner %in% sb_xgboost_cox_learner_ids() &&
    "prediction_horizon" %in% names(context)) {
    probability <- sb_predict_xgboost_cox_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    baseline_hazard <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      cumulative_baseline_hazard = probability$cumulative_baseline_hazard,
      baseline_survival = probability$baseline_survival,
      source = probability$source,
      ties = probability$ties
    )
  }
  survival_probability <- NULL
  if (model$learner %in% sb_xgboost_aft_learner_ids() &&
    "prediction_horizon" %in% names(context)) {
    probability <- sb_predict_xgboost_aft_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x,
      context$prediction_horizon
    )
    survival_probability <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      source = probability$source,
      ties = probability$ties
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      model_raw = sb_xgboost_model_raw(model),
      risk_reference = if (model$learner %in% sb_xgboost_cox_learner_ids()) {
        "xgboost_native"
      } else if (model$learner %in% sb_xgboost_aft_learner_ids()) {
        "xgboost_aft_normal_native"
      } else {
        NULL
      },
      risk_formula = if (model$learner %in% sb_xgboost_cox_learner_ids()) {
        "xgboost_raw_margin"
      } else if (model$learner %in% sb_xgboost_aft_learner_ids()) {
        "negative_xgboost_aft_log_time_margin"
      } else {
        NULL
      },
      baseline_hazard = baseline_hazard,
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_xgboost_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_load_xgboost_state_model <- function(state, learner) {
  x <- state[["model_raw"]]
  if (!is.raw(x) || length(x) == 0L || !is.list(state[["params"]])) {
    sb_bad_learner_state("The saved XGBoost learner state is invalid.")
  }
  sb_require_xgboost()
  booster <- tryCatch(
    if ("as_booster" %in% names(formals(xgboost::xgb.load.raw))) {
      xgboost::xgb.load.raw(x, as_booster = TRUE)
    } else {
      xgboost::xgb.load.raw(x)
    },
    error = function(error) {
      sb_abort(
        paste0(
          "The saved XGBoost model could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_xgboost_state",
        parent = error
      )
    }
  )
  if (!inherits(booster, "xgb.Booster")) {
    sb_bad_learner_state("The saved XGBoost booster is invalid.")
  }
  structure(
    list(
      learner = learner,
      features = state[["features"]],
      params = state[["params"]],
      booster = booster
    ),
    class = "stormbreaker_xgboost_fixed"
  )
}

sb_predict_xgboost_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_xgboost_learner_ids()) {
    sb_bad_learner_state("The saved XGBoost learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_load_xgboost_state_model(state, learner)
  prediction <- sb_xgboost_prediction(model, x)
  if (learner %in% sb_xgboost_binary_learner_ids()) {
    return(list(probability = prediction))
  }

  result <- list(risk = prediction)
  baseline <- state[["baseline_hazard"]]
  if (!is.null(baseline)) {
    valid <- is.list(baseline) &&
      is.numeric(baseline$prediction_horizon) &&
      length(baseline$prediction_horizon) == 1L &&
      !is.na(baseline$prediction_horizon) &&
      is.finite(baseline$prediction_horizon) &&
      baseline$prediction_horizon > 0 &&
      is.character(baseline$time_unit) && length(baseline$time_unit) == 1L &&
      !is.na(baseline$time_unit) && nzchar(baseline$time_unit)
    if (!valid) {
      sb_bad_learner_state(
        "The saved XGBoost-Cox baseline-hazard state is invalid."
      )
    }
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      prediction,
      baseline$cumulative_baseline_hazard
    )
  } else if (!is.null(state[["survival_probability"]])) {
    probability <- state[["survival_probability"]]
    valid <- learner %in% sb_xgboost_aft_learner_ids() &&
      is.list(probability) &&
      is.numeric(probability$prediction_horizon) &&
      length(probability$prediction_horizon) == 1L &&
      !is.na(probability$prediction_horizon) &&
      is.finite(probability$prediction_horizon) &&
      probability$prediction_horizon > 0 &&
      is.character(probability$time_unit) &&
      length(probability$time_unit) == 1L &&
      !is.na(probability$time_unit) && nzchar(probability$time_unit) &&
      identical(probability$source, "xgboost_aft_normal_native") &&
      identical(probability$ties, "not_applicable")
    if (!valid) {
      sb_bad_learner_state(
        "The saved XGBoost-AFT survival-probability state is invalid."
      )
    }
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- sb_xgboost_aft_survival_probability(
      prediction,
      probability$prediction_horizon,
      model$params$aft_loss_distribution_scale
    )
  }
  result
}

sb_xgboost_learner_adapter <- function(
  learner,
  params = sb_xgboost_params(learner)
) {
  params <- sb_validate_xgboost_params(learner, params)
  fit <- function(x, endpoint, event = NULL) {
    if (learner %in% sb_xgboost_binary_learner_ids()) {
      return(sb_fit_xgboost(x, endpoint, learner, params))
    }
    sb_fit_xgboost(
      x,
      list(time = endpoint, event = event),
      learner,
      params
    )
  }
  adapter <- list(
    id = learner,
    params = params,
    fit = fit,
    predict = sb_xgboost_prediction,
    extract_state = function(model, context = list()) {
      sb_extract_xgboost_learner_state(model, context)
    },
    predict_state = sb_predict_xgboost_learner_state
  )
  if (learner %in% sb_xgboost_cox_learner_ids()) {
    adapter$predict_survival <- sb_predict_xgboost_cox_survival
  } else if (learner %in% sb_xgboost_aft_learner_ids()) {
    adapter$predict_survival <- sb_predict_xgboost_aft_survival
  }
  adapter
}

# ---- model-catboost ----

# ---- catboost-aft-survival-fixed-adapter ----

sb_catboost_aft_survival_params <- function(
    iterations = 100L,
    learning_rate = 0.1,
    depth = 6L,
    l2_leaf_reg = 3,
    seed = 1L) {
  list(
    iterations = sb_check_positive_integer(
      iterations,
      "catboost_iterations",
      minimum = 1L
    ),
    learning_rate = sb_validate_catboost_learning_rate(learning_rate),
    depth = sb_validate_catboost_depth(depth),
    l2_leaf_reg = sb_validate_catboost_l2_leaf_reg(l2_leaf_reg),
    random_seed = sb_check_positive_integer(
      seed,
      "catboost_seed",
      minimum = 1L
    ),
    loss_function = "SurvivalAft:dist=Normal;scale=1",
    eval_metric = "SurvivalAft:dist=Normal;scale=1",
    task_type = "CPU",
    thread_count = 1L,
    bootstrap_type = "Bayesian",
    bagging_temperature = 1,
    random_strength = 1,
    rsm = 1,
    grow_policy = "SymmetricTree",
    boosting_type = "Plain",
    use_best_model = FALSE,
    allow_writing_files = FALSE,
    save_snapshot = FALSE,
    logging_level = "Silent"
  )
}

sb_validate_catboost_aft_survival_inputs <- function(x, time, event) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time <= 0)) {
    sb_abort(
      paste0(
        "The fixed CatBoost AFT learner requires strictly positive ",
        "survival times."
      ),
      "stormbreaker_catboost_aft_nonpositive_time"
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  if (!any(event == 1L)) {
    sb_abort(
      "A fixed CatBoost AFT analysis fold must contain at least one event.",
      "stormbreaker_catboost_aft_no_event"
    )
  }
  list(x = x, time = time, event = as.integer(event))
}

sb_fit_catboost_aft_survival <- function(x, time, event, params) {
  inputs <- sb_validate_catboost_aft_survival_inputs(x, time, event)
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for the fixed CatBoost AFT learner.",
      "stormbreaker_missing_catboost"
    )
  }

  training_pool <- catboost::catboost.load_pool(
    data = as.matrix(inputs$x),
    label = cbind(
      inputs$time,
      ifelse(inputs$event == 1L, inputs$time, -1)
    )
  )
  booster <- catboost::catboost.train(
    learn_pool = training_pool,
    test_pool = NULL,
    params = params
  )
  tree_count <- as.integer(booster$tree_count)
  if (!identical(tree_count, params$iterations)) {
    sb_abort(
      "Fixed CatBoost AFT did not fit the declared number of iterations.",
      "stormbreaker_incomplete_catboost_fit"
    )
  }

  structure(
    list(
      features = names(inputs$x),
      params = params,
      iterations = params$iterations,
      distribution = "Normal",
      scale = 1,
      booster = booster
    ),
    class = "stormbreaker_catboost_aft_survival"
  )
}

sb_catboost_aft_raw_margin <- function(model, x) {
  if (!inherits(model, "stormbreaker_catboost_aft_survival")) {
    sb_abort(
      "The supplied model is not a fitted fixed CatBoost AFT learner.",
      "stormbreaker_bad_catboost_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed CatBoost AFT schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for fixed CatBoost AFT prediction.",
      "stormbreaker_missing_catboost"
    )
  }

  margin <- as.numeric(catboost::catboost.predict(
    model$booster,
    catboost::catboost.load_pool(as.matrix(x)),
    prediction_type = "RawFormulaVal",
    ntree_start = 0L,
    ntree_end = model$iterations,
    thread_count = 1L
  ))
  if (length(margin) != nrow(x) || anyNA(margin) || any(!is.finite(margin))) {
    sb_abort(
      "Fixed CatBoost AFT produced invalid log-time predictions.",
      "stormbreaker_invalid_catboost_prediction"
    )
  }
  unname(margin)
}

sb_predict_catboost_aft_survival <- function(model, x) {
  -sb_catboost_aft_raw_margin(model, x)
}

sb_catboost_aft_survival_probability <- function(
    risk, prediction_horizon, scale = 1) {
  valid <- is.numeric(risk) && length(risk) > 0L && !anyNA(risk) &&
    all(is.finite(risk)) && is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L && !is.na(prediction_horizon) &&
    is.finite(prediction_horizon) && prediction_horizon > 0 &&
    is.numeric(scale) && length(scale) == 1L && !is.na(scale) &&
    is.finite(scale) && scale > 0
  if (!valid) {
    sb_abort(
      "Fixed CatBoost AFT received invalid prediction inputs.",
      "stormbreaker_invalid_catboost_aft_prediction"
    )
  }
  probability <- stats::pnorm(
    (-risk - log(prediction_horizon)) / scale
  )
  if (anyNA(probability) || any(!is.finite(probability)) ||
    any(probability < 0 | probability > 1)) {
    sb_abort(
      "Fixed CatBoost AFT produced invalid survival probabilities.",
      "stormbreaker_invalid_catboost_aft_prediction"
    )
  }
  unname(probability)
}

sb_predict_catboost_aft_survival_survival <- function(
    model, x_analysis, time, event, x_assessment, prediction_horizon) {
  inputs <- sb_validate_catboost_aft_survival_inputs(
    x_analysis,
    time,
    event
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fixed CatBoost AFT schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  risk <- sb_predict_catboost_aft_survival(model, x_assessment)
  list(
    prediction_horizon = prediction_horizon,
    source = "catboost_aft_normal_native",
    ties = "not_applicable",
    probability = sb_catboost_aft_survival_probability(
      risk,
      prediction_horizon,
      model$scale
    )
  )
}

sb_extract_catboost_aft_survival_learner_state <- function(
    model, context = list()) {
  extracted <- sb_extract_catboost_tree_model(model)
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_catboost_aft_survival_survival(
      model,
      context[["x"]],
      context[["time"]],
      context[["event"]],
      context[["x"]][1L, , drop = FALSE],
      context[["prediction_horizon"]]
    )
    survival_probability <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context[["time_unit"]]
      } else {
        NA_character_
      },
      source = probability$source,
      ties = probability$ties
    )
  }
  structure(
    c(
      list(
        learner = "catboost_aft_survival",
        distribution = model$distribution,
        scale = model$scale,
        risk_reference = "catboost_aft_normal_native",
        risk_formula = "negative_catboost_aft_log_time_margin"
      ),
      extracted,
      list(survival_probability = survival_probability)
    ),
    class = c(
      "stormbreaker_catboost_aft_survival_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_catboost_aft_survival_learner_state <- function(state, x) {
  x <- sb_validate_learner_state_x(
    state,
    x,
    "catboost_aft_survival"
  )
  valid_state <- is.raw(state[["model_raw"]]) &&
    length(state[["model_raw"]]) > 0L &&
    is.numeric(state[["iterations"]]) &&
    length(state[["iterations"]]) == 1L &&
    !is.na(state[["iterations"]]) && state[["iterations"]] >= 1 &&
    identical(state[["distribution"]], "Normal") &&
    is.numeric(state[["scale"]]) && length(state[["scale"]]) == 1L &&
    !is.na(state[["scale"]]) && is.finite(state[["scale"]]) &&
    state[["scale"]] > 0
  if (!valid_state) {
    sb_bad_learner_state("The saved CatBoost AFT learner state is invalid.")
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for fixed CatBoost AFT prediction.",
      "stormbreaker_missing_catboost"
    )
  }
  booster <- tryCatch(
    unserialize(state[["model_raw"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved CatBoost AFT state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_catboost_state",
        parent = error
      )
    }
  )
  model <- structure(
    list(
      features = state[["features"]],
      params = state[["params"]],
      iterations = as.integer(state[["iterations"]]),
      distribution = state[["distribution"]],
      scale = state[["scale"]],
      booster = booster
    ),
    class = "stormbreaker_catboost_aft_survival"
  )
  result <- list(risk = sb_predict_catboost_aft_survival(model, x))
  probability <- state[["survival_probability"]]
  if (!is.null(probability)) {
    valid_probability <- is.list(probability) &&
      is.numeric(probability[["prediction_horizon"]]) &&
      length(probability[["prediction_horizon"]]) == 1L &&
      !is.na(probability[["prediction_horizon"]]) &&
      is.finite(probability[["prediction_horizon"]]) &&
      probability[["prediction_horizon"]] > 0 &&
      is.character(probability[["time_unit"]]) &&
      length(probability[["time_unit"]]) == 1L &&
      !is.na(probability[["time_unit"]]) &&
      nzchar(probability[["time_unit"]]) &&
      identical(probability[["source"]], "catboost_aft_normal_native") &&
      identical(probability[["ties"]], "not_applicable")
    if (!valid_probability) {
      sb_bad_learner_state(
        "The saved CatBoost AFT survival-probability state is invalid."
      )
    }
    result$prediction_horizon <- probability[["prediction_horizon"]]
    result$time_unit <- probability[["time_unit"]]
    result$survival_probability <- sb_catboost_aft_survival_probability(
      result$risk,
      probability[["prediction_horizon"]],
      model$scale
    )
  }
  result
}

sb_catboost_aft_survival_learner_adapter <- function(
    iterations = 100L,
    learning_rate = 0.1,
    depth = 6L,
    l2_leaf_reg = 3,
    seed = 1L) {
  params <- sb_catboost_aft_survival_params(
    iterations = iterations,
    learning_rate = learning_rate,
    depth = depth,
    l2_leaf_reg = l2_leaf_reg,
    seed = seed
  )
  list(
    id = "catboost_aft_survival",
    params = params,
    fit = function(x, time, event) {
      sb_fit_catboost_aft_survival(x, time, event, params)
    },
    predict = sb_predict_catboost_aft_survival,
    predict_survival = sb_predict_catboost_aft_survival_survival,
    extract_tree_model = sb_extract_catboost_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_catboost_aft_survival_learner_state(model, context)
    },
    predict_state = sb_predict_catboost_aft_survival_learner_state
  )
}

# ---- catboost-binary-fixed-adapter ----

sb_catboost_available <- function() {
  requireNamespace("catboost", quietly = TRUE)
}

sb_validate_catboost_learning_rate <- function(learning_rate) {
  valid <- is.numeric(learning_rate) && length(learning_rate) == 1L &&
    !is.na(learning_rate) && is.finite(learning_rate) && learning_rate > 0
  if (!valid) {
    sb_abort(
      "`catboost_learning_rate` must be one finite number greater than 0.",
      "stormbreaker_bad_catboost_learning_rate"
    )
  }
  as.numeric(learning_rate)
}

sb_validate_catboost_depth <- function(depth) {
  depth <- sb_check_positive_integer(depth, "catboost_depth", minimum = 1L)
  if (depth > 16L) {
    sb_abort(
      "`catboost_depth` must be an integer between 1 and 16 for CPU training.",
      "stormbreaker_bad_catboost_depth"
    )
  }
  depth
}

sb_validate_catboost_l2_leaf_reg <- function(l2_leaf_reg) {
  valid <- is.numeric(l2_leaf_reg) && length(l2_leaf_reg) == 1L &&
    !is.na(l2_leaf_reg) && is.finite(l2_leaf_reg) && l2_leaf_reg > 0
  if (!valid) {
    sb_abort(
      "`catboost_l2_leaf_reg` must be one finite number greater than 0.",
      "stormbreaker_bad_catboost_l2_leaf_reg"
    )
  }
  as.numeric(l2_leaf_reg)
}

sb_catboost_binary_params <- function(
    iterations = 100L,
    learning_rate = 0.1,
    depth = 6L,
    l2_leaf_reg = 3,
    seed = 1L) {
  list(
    iterations = sb_check_positive_integer(
      iterations,
      "catboost_iterations",
      minimum = 1L
    ),
    learning_rate = sb_validate_catboost_learning_rate(learning_rate),
    depth = sb_validate_catboost_depth(depth),
    l2_leaf_reg = sb_validate_catboost_l2_leaf_reg(l2_leaf_reg),
    random_seed = sb_check_positive_integer(
      seed,
      "catboost_seed",
      minimum = 1L
    ),
    loss_function = "Logloss",
    eval_metric = "Logloss",
    task_type = "CPU",
    thread_count = 1L,
    bootstrap_type = "Bayesian",
    bagging_temperature = 1,
    random_strength = 1,
    rsm = 1,
    grow_policy = "SymmetricTree",
    boosting_type = "Plain",
    use_best_model = FALSE,
    boost_from_average = TRUE,
    allow_writing_files = FALSE,
    save_snapshot = FALSE,
    logging_level = "Silent"
  )
}

sb_fit_catboost_binary <- function(x, outcome, params) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "A fixed CatBoost analysis fold must contain both classes.",
      "stormbreaker_catboost_single_class"
    )
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for the fixed CatBoost learner.",
      "stormbreaker_missing_catboost"
    )
  }

  training_pool <- catboost::catboost.load_pool(
    data = as.matrix(x),
    label = outcome
  )
  booster <- catboost::catboost.train(
    learn_pool = training_pool,
    test_pool = NULL,
    params = params
  )
  tree_count <- as.integer(booster$tree_count)
  if (!identical(tree_count, params$iterations)) {
    sb_abort(
      "Fixed CatBoost did not fit the declared number of iterations.",
      "stormbreaker_incomplete_catboost_fit"
    )
  }

  structure(
    list(
      features = names(x),
      params = params,
      iterations = params$iterations,
      booster = booster
    ),
    class = "stormbreaker_catboost_binary"
  )
}

sb_predict_catboost_binary <- function(model, x) {
  if (!inherits(model, "stormbreaker_catboost_binary")) {
    sb_abort(
      "The supplied model is not a fitted fixed CatBoost learner.",
      "stormbreaker_bad_catboost_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed CatBoost schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for fixed CatBoost prediction.",
      "stormbreaker_missing_catboost"
    )
  }

  prediction_pool <- catboost::catboost.load_pool(as.matrix(x))
  probability <- as.numeric(catboost::catboost.predict(
    model$booster,
    prediction_pool,
    prediction_type = "Probability",
    ntree_start = 0L,
    ntree_end = model$iterations,
    thread_count = 1L
  ))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      "Fixed CatBoost produced invalid class-one probabilities.",
      "stormbreaker_invalid_catboost_prediction"
    )
  }
  probability
}

sb_extract_catboost_tree_model <- function(model) {
  if (!inherits(model, c(
    "stormbreaker_catboost_binary",
    "stormbreaker_catboost_cox",
    "stormbreaker_catboost_aft_survival"
  ))) {
    sb_abort(
      "The supplied model is not a fitted fixed CatBoost learner.",
      "stormbreaker_bad_catboost_fit"
    )
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required to audit a CatBoost model.",
      "stormbreaker_missing_catboost"
    )
  }

  importance <- as.numeric(catboost::catboost.get_feature_importance(
    model$booster,
    type = "FeatureImportance",
    thread_count = 1L
  ))
  valid_importance <- length(importance) == length(model$features) &&
    !anyNA(importance) && all(is.finite(importance))
  if (!valid_importance) {
    sb_abort(
      "Fixed CatBoost did not return one finite importance per feature.",
      "stormbreaker_invalid_catboost_importance"
    )
  }

  list(
    iterations = model$iterations,
    tree_count = as.integer(model$booster$tree_count),
    features = model$features,
    params = model$params,
    model_serialization = "R_xdr_v3_catboost.Model",
    model_raw = serialize(
      model$booster,
      connection = NULL,
      ascii = FALSE,
      xdr = TRUE,
      version = 3L
    ),
    feature_importance = data.frame(
      feature = model$features,
      importance = importance,
      stringsAsFactors = FALSE
    )
  )
}

sb_catboost_binary_learner_adapter <- function(
    iterations = 100L,
    learning_rate = 0.1,
    depth = 6L,
    l2_leaf_reg = 3,
    seed = 1L) {
  params <- sb_catboost_binary_params(
    iterations = iterations,
    learning_rate = learning_rate,
    depth = depth,
    l2_leaf_reg = l2_leaf_reg,
    seed = seed
  )
  list(
    id = "catboost_binary",
    params = params,
    fit = function(x, outcome) {
      sb_fit_catboost_binary(x, outcome, params)
    },
    predict = sb_predict_catboost_binary,
    extract_tree_model = sb_extract_catboost_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_catboost_learner_state(model)
    },
    predict_state = sb_predict_catboost_learner_state
  )
}

# ---- catboost-cox-fixed-adapter ----

sb_catboost_cox_params <- function(
    iterations = 100L,
    learning_rate = 0.1,
    depth = 6L,
    l2_leaf_reg = 3,
    seed = 1L) {
  list(
    iterations = sb_check_positive_integer(
      iterations,
      "catboost_iterations",
      minimum = 1L
    ),
    learning_rate = sb_validate_catboost_learning_rate(learning_rate),
    depth = sb_validate_catboost_depth(depth),
    l2_leaf_reg = sb_validate_catboost_l2_leaf_reg(l2_leaf_reg),
    random_seed = sb_check_positive_integer(
      seed,
      "catboost_seed",
      minimum = 1L
    ),
    loss_function = "Cox",
    eval_metric = "Cox",
    task_type = "CPU",
    thread_count = 1L,
    bootstrap_type = "Bayesian",
    bagging_temperature = 1,
    random_strength = 1,
    rsm = 1,
    grow_policy = "SymmetricTree",
    boosting_type = "Plain",
    use_best_model = FALSE,
    allow_writing_files = FALSE,
    save_snapshot = FALSE,
    logging_level = "Silent"
  )
}

sb_fit_catboost_cox <- function(x, time, event, params) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time <= 0)) {
    sb_abort(
      paste0(
        "The fixed CatBoost Cox learner requires strictly positive ",
        "survival times."
      ),
      "stormbreaker_catboost_cox_nonpositive_time"
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  if (!any(event == 1L)) {
    sb_abort(
      "A fixed CatBoost Cox analysis fold must contain at least one event.",
      "stormbreaker_catboost_cox_no_event"
    )
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for the fixed CatBoost Cox learner.",
      "stormbreaker_missing_catboost"
    )
  }

  training_pool <- catboost::catboost.load_pool(
    data = as.matrix(x),
    label = ifelse(event == 1L, time, -time)
  )
  booster <- catboost::catboost.train(
    learn_pool = training_pool,
    test_pool = NULL,
    params = params
  )
  tree_count <- as.integer(booster$tree_count)
  if (!identical(tree_count, params$iterations)) {
    sb_abort(
      "Fixed CatBoost Cox did not fit the declared number of iterations.",
      "stormbreaker_incomplete_catboost_fit"
    )
  }

  structure(
    list(
      features = names(x),
      params = params,
      iterations = params$iterations,
      booster = booster
    ),
    class = "stormbreaker_catboost_cox"
  )
}

sb_predict_catboost_cox <- function(model, x) {
  if (!inherits(model, "stormbreaker_catboost_cox")) {
    sb_abort(
      "The supplied model is not a fitted fixed CatBoost Cox learner.",
      "stormbreaker_bad_catboost_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed CatBoost Cox schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for fixed CatBoost Cox prediction.",
      "stormbreaker_missing_catboost"
    )
  }

  prediction_pool <- catboost::catboost.load_pool(as.matrix(x))
  risk <- as.numeric(catboost::catboost.predict(
    model$booster,
    prediction_pool,
    prediction_type = "RawFormulaVal",
    ntree_start = 0L,
    ntree_end = model$iterations,
    thread_count = 1L
  ))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "Fixed CatBoost Cox produced invalid risk predictions.",
      "stormbreaker_invalid_catboost_prediction"
    )
  }
  unname(risk)
}

sb_predict_catboost_cox_survival <- function(
    model, x_analysis, time, event, x_assessment, prediction_horizon) {
  x_analysis <- sb_validate_predictor_frame(x_analysis, "x_analysis")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x_analysis))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x_analysis))
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    time
  )
  analysis_risk <- sb_predict_catboost_cox(model, x_analysis)
  assessment_risk <- sb_predict_catboost_cox(model, x_assessment)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    time,
    as.integer(event),
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = "catboost_cox_breslow_raw_margin",
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_catboost_cox_learner_adapter <- function(
    iterations = 100L,
    learning_rate = 0.1,
    depth = 6L,
    l2_leaf_reg = 3,
    seed = 1L) {
  params <- sb_catboost_cox_params(
    iterations = iterations,
    learning_rate = learning_rate,
    depth = depth,
    l2_leaf_reg = l2_leaf_reg,
    seed = seed
  )
  list(
    id = "catboost_cox",
    params = params,
    fit = function(x, time, event) {
      sb_fit_catboost_cox(x, time, event, params)
    },
    predict = sb_predict_catboost_cox,
    predict_survival = sb_predict_catboost_cox_survival,
    extract_tree_model = sb_extract_catboost_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_catboost_cox_learner_state(model, context)
    },
    predict_state = sb_predict_catboost_cox_learner_state
  )
}

# ---- model-lightgbm ----

# ---- algorithm-lightgbm-catboost ----

sb_catboost_diagnostic_ids <- function() {
  c("catboost_binary", "catboost_cox", "catboost_aft_survival")
}

sb_light_cat_diagnostic <- function(fit, type, max_pairs) {
  id <- fit$pipeline$learner
  state <- fit$learner_state
  light <- id %in% sb_lightgbm_binary_learner_ids()
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("importance", "tree_size", "interaction", "training_loss", "cv_curve")) sb_abort(
    "Choose a supported importance, tree_size or interaction diagnostic.",
    "stormbreaker_unsupported_diagnostic")
  if (type %in% c("training_loss", "cv_curve")) sb_abort(
    "No original evaluation history is retained by these pipelines; none is reconstructed.",
    "stormbreaker_unavailable_diagnostic")
  native <- NULL
  requested <- NULL
  display_limit <- NULL
  if (type == "importance") {
    tab <- state$feature_importance
    required <- if (light) c("gain", "cover", "frequency") else "importance"
    if (!is.data.frame(tab) || !all(c("feature", required) %in% names(tab)) ||
        anyDuplicated(tab$feature) || !setequal(tab$feature, state$features) ||
        !all(vapply(tab[required], function(z) is.numeric(z) && all(is.finite(z)), logical(1)))) sb_abort(
      "Saved native feature importance is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    native <- tab
    tab$value <- if (light) tab$gain else tab$importance
    measure <- if (light) "Native LightGBM gain (unscaled)" else "Native CatBoost FeatureImportance"
    interpretation <- if (light) "Saved gain, cover and split counts; not normalized XGBoost gain shares." else
      "Original CatBoost importance; not gain, SHAP, effect direction or a causal effect."
  } else if (type == "tree_size") {
    if (!light) sb_abort("This tree-size adapter is for LightGBM, not CatBoost.",
      "stormbreaker_unavailable_diagnostic")
    if (!requireNamespace("lightgbm", quietly = TRUE)) sb_abort("Package lightgbm is required.",
      "stormbreaker_missing_dependency")
    if (!requireNamespace("jsonlite", quietly = TRUE)) sb_abort("Package jsonlite is required to read the native LightGBM dump.",
      "stormbreaker_missing_dependency")
    booster <- lightgbm::lgb.load(model_str = state$model_string)
    native <- jsonlite::fromJSON(lightgbm::lgb.dump(booster), simplifyVector = FALSE)$tree_info
    leaves <- vapply(native, function(z) as.numeric(z$num_leaves), numeric(1))
    index <- vapply(native, function(z) as.numeric(z$tree_index), numeric(1))
    if (!length(leaves) || any(!is.finite(leaves)) || any(leaves < 1) ||
        any(leaves != floor(leaves)) || !identical(index, as.numeric(seq_along(leaves) - 1L))) sb_abort(
      "Native LightGBM tree counts or indices are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(tree = index + 1L, native_tree_index = index, terminal_nodes = as.integer(leaves))
    # Do not retain the full dump: the diagnostic needs counts, not every split.
    native <- NULL
    requested <- state$nrounds
    measure <- "Terminal nodes per original LightGBM tree"
    interpretation <- paste0("Actual saved trees: ", nrow(tab), "; configured rounds: ", requested,
      ". Fewer trees can be retained when no further valid split is available.")
  } else {
    if (light) sb_abort("This native interaction adapter is for CatBoost, not LightGBM.",
      "stormbreaker_unavailable_diagnostic")
    display_limit <- sb_check_positive_integer(max_pairs, "max_pairs", minimum = 1L)
    if (!requireNamespace("catboost", quietly = TRUE)) sb_abort("Package catboost is required.",
      "stormbreaker_missing_dependency")
    booster <- unserialize(state$model_raw)
    # With one input there are no distinct feature pairs. The native R wrapper
    # tries to name three columns on an empty result and fails before returning.
    native <- if (length(state$features) == 1L) {
      matrix(numeric(), nrow = 0L, ncol = 3L)
    } else {
      tryCatch(
        catboost::catboost.get_feature_importance(booster, type = "Interaction", thread_count = 1L),
        error = function(e) sb_abort(
          paste0("Native CatBoost interaction extraction failed: ", conditionMessage(e),
            ". No interaction scores were substituted."),
          "stormbreaker_unavailable_diagnostic", parent = e)
      )
    }
    if (is.null(native) || !length(native)) native <- matrix(numeric(), nrow = 0L, ncol = 3L)
    if (!is.matrix(native) || ncol(native) != 3L || !is.numeric(native) || any(!is.finite(native)) ||
        any(native[, 1:2, drop = FALSE] < 0) || any(native[, 1:2, drop = FALSE] >= length(state$features)) ||
        any(native[, 1:2, drop = FALSE] != floor(native[, 1:2, drop = FALSE])) ||
        any(native[, 1] == native[, 2])) sb_abort(
      "Native CatBoost interaction indices or scores are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(feature_1 = state$features[native[, 1] + 1L],
      feature_2 = state$features[native[, 2] + 1L], score = unname(native[, 3]),
      native_index_1 = unname(native[, 1]), native_index_2 = unname(native[, 2]))
    measure <- "Native CatBoost interaction strength"
    interpretation <- if (length(state$features) == 1L) {
      "Only one input feature: no distinct feature pairs exist; the native pair extractor was not called."
    } else {
      "Model-internal interaction strength, not SHAP interaction, clinical synergy or causality."
    }
  }
  structure(list(type = paste0("light_cat_", type), table = tab, native_table = native,
    learner = id, pipeline_id = fit$pipeline$id, measure = measure,
    interpretation = interpretation, requested_iterations = requested, max_pairs = display_limit,
    source = "original saved boosting model; no training, permutation or patient prediction"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_light_cat_diagnostic <- function(x) {
  if (x$type == "light_cat_importance") {
    z <- x
    z$type <- "boosted_tree_importance"
    return(sb_plot_gbm_xgboost_diagnostic(z))
  }
  if (x$type == "light_cat_tree_size") {
    z <- x
    z$type <- "forest_tree_size"
    return(sb_plot_forest_diagnostic(z) + ggplot2::labs(
      title = "LightGBM tree-size distribution", subtitle = x$learner,
      caption = paste0(x$interpretation, "\nAll actual trees counted; no artificial padding to the configured iteration count.")))
  }
  tab <- x$table
  tab <- tab[utils::head(order(tab$score, decreasing = TRUE), x$max_pairs), , drop = FALSE]
  pair_labels <- stats::setNames(paste(tab$feature_1, tab$feature_2, sep = " + "), seq_len(nrow(tab)))
  tab$pair <- factor(seq_len(nrow(tab)), levels = rev(seq_len(nrow(tab))))
  score <- pair <- NULL
  p <- ggplot2::ggplot(tab, ggplot2::aes(score, pair)) +
    ggplot2::geom_col(fill = "#2878A0", width = .65) +
    ggplot2::scale_y_discrete(labels = pair_labels) +
    ggplot2::labs(title = x$measure,
      subtitle = paste0(x$learner, ": showing ", nrow(tab), " of ", nrow(x$table), " reported pairs"),
      x = "Native interaction score", y = NULL,
      caption = paste0(x$interpretation, "\nComplete reported pairs remain in the table; unreported pairs are not assigned zero scores.")) + sb_viz_theme()
  if (!nrow(tab)) {
    p <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
      label = "No distinct feature pairs are available to display.") + ggplot2::theme_void() +
      ggplot2::labs(title = x$measure, subtitle = x$learner, caption = x$interpretation)
  }
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = max(5, 2.5 + .28 * nrow(tab)))
  p
}

# ---- lightgbm-binary-fixed-adapter ----

sb_lightgbm_binary_learner_ids <- function() {
  c(
    "lightgbm_binary",
    "lightgbm_binary_dart",
    "lightgbm_binary_goss",
    "lightgbm_binary_rf"
  )
}

sb_lightgbm_binary_variant <- function(learner) {
  variant <- switch(learner,
    lightgbm_binary = "gbdt",
    lightgbm_binary_dart = "dart",
    lightgbm_binary_goss = "goss",
    lightgbm_binary_rf = "rf",
    NULL
  )
  if (is.null(variant)) {
    sb_abort(
      paste0("Unknown fixed LightGBM learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  variant
}

sb_lightgbm_available <- function() {
  requireNamespace("lightgbm", quietly = TRUE)
}

sb_validate_lightgbm_learning_rate <- function(learning_rate) {
  valid <- is.numeric(learning_rate) && length(learning_rate) == 1L &&
    !is.na(learning_rate) && is.finite(learning_rate) && learning_rate > 0
  if (!valid) {
    sb_abort(
      "`lightgbm_learning_rate` must be one finite number greater than 0.",
      "stormbreaker_bad_lightgbm_learning_rate"
    )
  }
  as.numeric(learning_rate)
}

sb_lightgbm_binary_params <- function(
    nrounds = 100L,
    learning_rate = 0.1,
    num_leaves = 31L,
    min_data_in_leaf = 20L,
    seed = 1L,
    learner = "lightgbm_binary") {
  variant <- sb_lightgbm_binary_variant(learner)
  params <- list(
    nrounds = sb_check_positive_integer(
      nrounds,
      "lightgbm_nrounds",
      minimum = 1L
    ),
    learning_rate = sb_validate_lightgbm_learning_rate(learning_rate),
    num_leaves = sb_check_positive_integer(
      num_leaves,
      "lightgbm_num_leaves",
      minimum = 2L
    ),
    min_data_in_leaf = sb_check_positive_integer(
      min_data_in_leaf,
      "lightgbm_min_data_in_leaf",
      minimum = 0L
    ),
    seed = sb_check_positive_integer(
      seed,
      "lightgbm_seed",
      minimum = 1L
    ),
    objective = "binary",
    metric = "binary_logloss",
    boosting = if (identical(variant, "dart")) {
      "dart"
    } else if (identical(variant, "rf")) {
      "rf"
    } else {
      "gbdt"
    },
    device_type = "cpu",
    num_threads = 1L,
    deterministic = TRUE,
    force_col_wise = TRUE,
    feature_fraction = if (identical(variant, "rf")) 0.8 else 1,
    bagging_fraction = if (identical(variant, "rf")) 0.8 else 1,
    bagging_freq = if (identical(variant, "rf")) 1L else 0L,
    max_depth = -1L,
    lambda_l1 = 0,
    lambda_l2 = 0,
    boost_from_average = TRUE,
    verbosity = -1L
  )
  if (identical(variant, "dart")) {
    return(c(params, list(
      drop_rate = 0.1,
      max_drop = 50L,
      skip_drop = 0,
      uniform_drop = FALSE,
      xgboost_dart_mode = FALSE
    )))
  }
  if (identical(variant, "goss")) {
    return(c(params, list(
      data_sample_strategy = "goss",
      top_rate = 0.2,
      other_rate = 0.1
    )))
  }
  params
}

sb_fit_lightgbm_binary <- function(x, outcome, params, learner) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "A fixed LightGBM analysis fold must contain both classes.",
      "stormbreaker_lightgbm_single_class"
    )
  }
  if (!sb_lightgbm_available()) {
    sb_abort(
      "Package `lightgbm` is required for the fixed LightGBM learner.",
      "stormbreaker_missing_lightgbm"
    )
  }

  matrix_x <- as.matrix(x)
  training_params <- params[setdiff(names(params), "nrounds")]
  # Native alias avoids BiocGenerics format.list truncating binary_logloss.
  training_params$metric <- "binary"
  dataset <- lightgbm::lgb.Dataset(
    data = matrix_x,
    label = outcome,
    colnames = names(x),
    free_raw_data = TRUE
  )
  booster <- lightgbm::lgb.train(
    params = training_params,
    data = dataset,
    nrounds = params$nrounds,
    verbose = -1L,
    record = FALSE,
    early_stopping_rounds = NULL,
    serializable = TRUE
  )

  structure(
    list(
      learner = learner,
      features = names(x),
      params = params,
      nrounds = params$nrounds,
      booster = booster
    ),
    class = "stormbreaker_lightgbm_binary"
  )
}

sb_predict_lightgbm_binary <- function(model, x) {
  if (!inherits(model, "stormbreaker_lightgbm_binary")) {
    sb_abort(
      "The supplied model is not a fitted fixed LightGBM learner.",
      "stormbreaker_bad_lightgbm_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed LightGBM schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (!sb_lightgbm_available()) {
    sb_abort(
      "Package `lightgbm` is required for fixed LightGBM prediction.",
      "stormbreaker_missing_lightgbm"
    )
  }

  probability <- as.numeric(stats::predict(
    model$booster,
    as.matrix(x),
    type = "response",
    num_iteration = model$nrounds,
    params = list(num_threads = 1L)
  ))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      "Fixed LightGBM produced invalid class-one probabilities.",
      "stormbreaker_invalid_lightgbm_prediction"
    )
  }
  probability
}

sb_extract_lightgbm_tree_model <- function(model) {
  if (!inherits(model, "stormbreaker_lightgbm_binary")) {
    sb_abort(
      "The supplied model is not a fitted fixed LightGBM learner.",
      "stormbreaker_bad_lightgbm_fit"
    )
  }
  if (!sb_lightgbm_available()) {
    sb_abort(
      "Package `lightgbm` is required to audit a LightGBM model.",
      "stormbreaker_missing_lightgbm"
    )
  }

  model_string <- model$booster$save_model_to_string(
    num_iteration = model$nrounds
  )
  feature_line <- grep(
    "^feature_names=",
    strsplit(model_string, "\n", fixed = TRUE)[[1L]],
    value = TRUE
  )
  native_features <- strsplit(
    sub("^feature_names=", "", feature_line[[1L]]),
    " ",
    fixed = TRUE
  )[[1L]]

  observed <- if (!is.null(model$reused_native_importance) &&
    identical(model_string, model$reused_native_importance$model_string)) {
    model$reused_native_importance$observed
  } else as.data.frame(
    lightgbm::lgb.importance(model$booster, percentage = FALSE),
    stringsAsFactors = FALSE
  )
  importance <- data.frame(
    feature = model$features,
    gain = numeric(length(model$features)),
    cover = numeric(length(model$features)),
    frequency = integer(length(model$features)),
    stringsAsFactors = FALSE
  )
  if (nrow(observed) > 0L) {
    index <- match(observed$Feature, native_features)
    if (anyNA(index)) {
      sb_abort(
        "LightGBM importance names could not be mapped to the fitted features.",
        "stormbreaker_lightgbm_importance_mismatch"
      )
    }
    importance$gain[index] <- as.numeric(observed$Gain)
    importance$cover[index] <- as.numeric(observed$Cover)
    importance$frequency[index] <- as.integer(observed$Frequency)
  }

  list(
    nrounds = model$nrounds,
    features = model$features,
    params = model$params,
    model_string = model_string,
    feature_importance = importance
  )
}

sb_lightgbm_binary_learner_adapter <- function(
    nrounds = 100L,
    learning_rate = 0.1,
    num_leaves = 31L,
    min_data_in_leaf = 20L,
    seed = 1L,
    learner = "lightgbm_binary") {
  params <- sb_lightgbm_binary_params(
    nrounds = nrounds,
    learning_rate = learning_rate,
    num_leaves = num_leaves,
    min_data_in_leaf = min_data_in_leaf,
    seed = seed,
    learner = learner
  )
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_lightgbm_binary(x, outcome, params, learner)
    },
    predict = sb_predict_lightgbm_binary,
    extract_tree_model = sb_extract_lightgbm_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_lightgbm_learner_state(model)
    },
    predict_state = sb_predict_lightgbm_learner_state
  )
}
