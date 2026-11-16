#' Fixed signatures and feature mapping
#' @name sb_linear_signature
#' @description Apply fixed prediction formulas and explicit one-to-one feature mappings without refitting a model.
#' @rawRd \arguments{
#' \item{signature_id}{
#' One non-empty identifier for this signature.
#' }
#' \item{task_type}{
#' Exactly "binary" or "survival".
#' }
#' \item{coefficients}{
#' A non-empty finite numeric vector named by unique,
#' non-empty predictor column names.
#' }
#' \item{intercept}{
#' For a binary signature, one explicitly supplied finite
#' numeric intercept. For a survival signature, one finite numeric intercept;
#' \code{NULL} uses the conventional zero-reference linear predictor.
#' }
#' \item{source}{
#' Optional free-text citation or source identifier. It is stored
#' but not verified.
#' }
#' \item{input_description}{
#' Optional free-text description of the required
#' feature identifiers and preprocessing scale. It is stored but not applied.
#' }
#' \item{prediction_horizon}{
#' For survival only, an optional finite positive
#' horizon. It must be supplied together with time_unit and
#' baseline_survival.
#' }
#' \item{time_unit}{
#' For a survival horizon, exactly "day", "month", or "year".
#' }
#' \item{baseline_survival}{
#' For a survival horizon, the published baseline
#' survival probability at zero linear score, between zero and one inclusive.
#' }
#' \item{new_data}{
#' \strong{For \code{sb_apply_feature_map()}:} A non-empty data frame containing one unique patient ID and
#' every source feature named in \code{feature_map}. Other columns are ignored.
#' 
#' \strong{For \code{sb_predict_signature()}:} A non-empty data frame containing one unique patient ID and
#' every exactly named numeric signature feature.
#' }
#' \item{id_col}{
#' \strong{For \code{sb_apply_feature_map()}:} Name of the unique patient identifier column.
#' 
#' \strong{For \code{sb_predict_signature()}:} Name of the patient-ID column in new_data.
#' }
#' \item{feature_map}{
#' A non-empty data frame with exactly one column named
#' \code{source_feature} and one named \code{target_feature}. Both must contain unique,
#' non-empty character identifiers. Additional provenance columns are kept.
#' }
#' \item{target_features}{
#' Unique, non-empty target feature names. Their set must
#' exactly match \code{feature_map$target_feature}, and their order controls the
#' returned predictor order.
#' }
#' \item{spec}{
#' One \code{stormbreaker_linear_signature} object or one named list of
#' arguments accepted by \code{\link[=sb_linear_signature]{sb_linear_signature()}}.
#' }
#' \item{signature}{
#' A stormbreaker_linear_signature returned by
#' \code{\link[=sb_linear_signature]{sb_linear_signature()}}.
#' }
#' \item{records}{
#' A non-empty named list. Each element is one argument list or
#' \code{stormbreaker_linear_signature} accepted by \code{\link[=sb_lint_signature]{sb_lint_signature()}}. Names
#' must be unique and non-empty and are returned unchanged as \code{record_id}.
#' A single specification or signature object must still be wrapped, for
#' example \code{list(paper_a = specification)}.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_linear_signature()}}{
#' sb_linear_signature() records one explicit, fixed linear formula without
#' fitting, feature selection, tuning, or outcome access. Coefficient names are
#' the exact predictor column names required later by \code{\link[=sb_predict_signature]{sb_predict_signature()}}.
#' Their order is retained.
#' 
#' \strong{Returns:} A stormbreaker_linear_signature containing the exact formula,
#' prediction semantics, optional source metadata, and required feature order.
#' 
#' \strong{Details:} For a binary signature, the only supported interpretation is a published
#' logistic formula: probability(class 1) equals
#' plogis(intercept + sum(x_j * coefficient_j)).
#' The intercept must be supplied explicitly, including when its stated value
#' is zero; an omitted binary intercept cannot define an absolute probability.
#' 
#' For a survival signature, the linear score is returned as a risk value, with
#' larger values meaning higher predicted risk. If and only if a prediction
#' horizon, its time unit, and the baseline survival at zero linear score are
#' all supplied, survival probability is calculated as
#' S(h | x) = S0(h) ^ exp(intercept + sum(x_j * coefficient_j)).
#' 
#' The function deliberately does not infer gene identifiers, map probes,
#' transform or standardize expression, impute values, orient coefficients,
#' recover an omitted intercept or baseline survival, or verify a citation.
#' A caller-reviewed one-to-one map can be applied separately with
#' \code{\link[=sb_apply_feature_map]{sb_apply_feature_map()}}; signature prediction itself remains exact-name
#' only.
#' source and input_description are retained as free-text audit metadata only.
#' The caller remains responsible for establishing that predictor columns use
#' the exact identifiers and scale required by the original formula. A horizon
#' probability is not available when the publication gives only a relative-risk
#' score.
#' 
#' This object is a reproducible formula carrier, not a literature-evidence
#' registry entry and not proof that a signature was independently validated.
#' Like ordinary R objects it is mutable; "locked" means that prediction does
#' not learn from evaluation data, not that the object is cryptographically
#' immutable.
#' }
#' 
#' \subsection{\code{sb_apply_feature_map()}}{
#' \code{sb_apply_feature_map()} renames and reorders numeric new-patient features
#' using a caller-supplied mapping table. It performs no database lookup or
#' annotation inference. Only the requested patient ID and mapped source
#' columns are read; outcome and unrelated columns are omitted from the returned
#' data.
#' 
#' \strong{Returns:} A \code{stormbreaker_mapped_features} object containing the reordered
#' mapping table and an outcome-free \code{data} frame with the ID followed by the
#' mapped numeric target features.
#' 
#' \strong{Details:} \code{feature_map} must contain one unique \code{source_feature} for every unique
#' \code{target_feature}, and its target set must exactly equal \code{target_features}.
#' The order of \code{target_features} defines the returned predictor order. Extra
#' mapping-table columns are retained as caller-supplied provenance.
#' 
#' One-to-many and many-to-one mappings are deliberately rejected. The function
#' never chooses a first annotation, aggregates probes or transcripts, strips
#' identifier versions, changes case, normalizes expression, imputes values, or
#' fills an absent feature with zero. Those decisions require an explicit,
#' source-specific biological justification before this function is called.
#' 
#' A successful result can supply its \code{data} element directly to
#' \code{\link[=sb_predict_final]{sb_predict_final()}} or \code{\link[=sb_predict_signature]{sb_predict_signature()}}. Acceptance confirms only
#' that the supplied map is structurally one-to-one and executable on these
#' columns. It does not verify the species, annotation release, platform,
#' biological equivalence, expression scale, or cited source.
#' }
#' 
#' \subsection{\code{sb_lint_signature()}}{
#' \code{sb_lint_signature()} checks one \code{\link[=sb_linear_signature]{sb_linear_signature()}} object or one named
#' list of arguments for that constructor. It does not inspect patient data,
#' read outcomes, fit a model, verify a citation, or infer omitted formula or
#' preprocessing details. Expected specification problems are retained as a
#' one-row \code{rejected} report rather than disappearing when construction fails.
#' Unexpected internal errors remain errors.
#' 
#' \strong{Returns:} A one-row data frame containing admission status, safely reportable
#' signature metadata, output type, and an issue code and message for a
#' rejected specification.
#' 
#' \strong{Details:} \code{admitted} means only that STORMBREAKER can represent the stated formula.
#' It does not establish that the transcription is faithful to a publication,
#' that input identifiers or scales are compatible, that a cohort is suitable,
#' or that the signature has predictive or clinical validity.
#' }
#' 
#' \subsection{\code{sb_predict_signature()}}{
#' sb_predict_signature() applies the exact formula stored by
#' \code{\link[=sb_linear_signature]{sb_linear_signature()}} to a non-empty new-patient data frame. It reads only
#' the requested ID column and the signature's named numeric features. Unrelated
#' columns, including accidentally supplied endpoint columns, do not enter the
#' calculation. No coefficient, intercept, baseline survival, feature,
#' transformation, threshold, or model component is learned or changed.
#' 
#' \strong{Returns:} A successful or failed stormbreaker_final_prediction compatible
#' with the existing prediction-first external-evaluation workflow.
#' 
#' \strong{Details:} The result is a stormbreaker_final_prediction, so it can be passed directly
#' to \code{\link[=sb_evaluate_external]{sb_evaluate_external()}}, \code{\link[=sb_compare_external]{sb_compare_external()}}, and
#' \code{\link[=sb_bootstrap_external_comparison]{sb_bootstrap_external_comparison()}}. Binary signatures return the
#' probability of class 1. Survival signatures always return a linear risk
#' score and return a horizon survival probability only when the signature
#' explicitly stored all required horizon information.
#' 
#' Feature names are matched exactly and reordered to the coefficient order.
#' Missing, duplicated, non-numeric, missing-valued, or non-finite signature
#' inputs produce a stage-aware failed prediction. The function never drops a
#' feature, imputes a value, refits a formula, or substitutes another model.
#' }
#' 
#' \subsection{\code{sb_signature_inventory()}}{
#' \code{sb_signature_inventory()} applies \code{\link[=sb_lint_signature]{sb_lint_signature()}} sequentially to an
#' ordered, named list of literature records. Every input produces one output
#' row in the same order, including incomplete or malformed records. The list
#' names become \code{record_id}, a caller-supplied key that remains available even
#' when a rejected record has no usable \code{signature_id}.
#' 
#' \strong{Returns:} A data frame with one row per supplied record, in input order. It
#' adds \code{record_id} before the columns returned by \code{\link[=sb_lint_signature]{sb_lint_signature()}}.
#' 
#' \strong{Details:} This is an in-memory audit inventory, not a built-in or verified literature
#' registry. It does not search publications, verify citations, infer missing
#' formula details, inspect patient data, execute predictions, merge duplicate
#' signatures, rank records, or choose a model. Repeated \code{signature_id} values
#' are allowed because separate source records may describe the same formula;
#' unique non-empty \code{record_id} names keep those source records distinguishable.
#' Unexpected internal linter errors remain errors rather than being converted
#' into rejected literature records.
#' }
#' }
#' @rawRd \examples{
#' # sb_linear_signature
#' signature <- sb_linear_signature(
#'   "example_binary_signature",
#'   "binary",
#'   coefficients = c(gene_a = 0.7, gene_b = -0.4),
#'   intercept = -0.2,
#'   input_description = "Values prepared on the formula's stated scale"
#' )
#' signature
#' 
#' # sb_apply_feature_map
#' new_data <- data.frame(
#'   patient_id = c("p1", "p2"),
#'   ENSG_B = c(2, 3),
#'   ENSG_A = c(0.1, 0.4)
#' )
#' feature_map <- data.frame(
#'   source_feature = c("ENSG_B", "ENSG_A"),
#'   target_feature = c("gene_b", "gene_a"),
#'   annotation_release = "caller_declared"
#' )
#' sb_apply_feature_map(
#'   new_data,
#'   "patient_id",
#'   feature_map,
#'   c("gene_a", "gene_b")
#' )
#' 
#' # sb_lint_signature
#' sb_lint_signature(list(
#'   signature_id = "complete_binary_formula",
#'   task_type = "binary",
#'   coefficients = c(gene_a = 0.7),
#'   intercept = -0.2
#' ))
#' 
#' sb_lint_signature(list(
#'   signature_id = "intercept_not_reported",
#'   task_type = "binary",
#'   coefficients = c(gene_a = 0.7)
#' ))
#' 
#' # sb_predict_signature
#' signature <- sb_linear_signature(
#'   "example_binary_signature",
#'   "binary",
#'   coefficients = c(gene_a = 0.7, gene_b = -0.4),
#'   intercept = -0.2
#' )
#' new_data <- data.frame(
#'   patient_id = c("p1", "p2"),
#'   gene_a = c(0.1, 0.5),
#'   gene_b = c(-0.2, 0.4)
#' )
#' sb_predict_signature(signature, new_data, "patient_id")
#' 
#' # sb_signature_inventory
#' records <- list(
#'   paper_a = list(
#'     signature_id = "complete_formula",
#'     task_type = "binary",
#'     coefficients = c(gene_a = 0.7),
#'     intercept = -0.2,
#'     source = "Citation A"
#'   ),
#'   paper_b = list(
#'     signature_id = "intercept_not_reported",
#'     task_type = "binary",
#'     coefficients = c(gene_a = 0.5),
#'     source = "Citation B"
#'   )
#' )
#' sb_signature_inventory(records)
#' }
NULL

# ---- feature-specification ----

# ---- feature-mapping ----

sb_validate_feature_map_targets <- function(target_features) {
  valid <- is.character(target_features) && length(target_features) > 0L &&
    !anyNA(target_features) && all(nzchar(trimws(target_features))) &&
    !anyDuplicated(target_features)
  if (!valid) {
    sb_abort(
      "`target_features` must be unique, non-empty character names.",
      "stormbreaker_bad_feature_map_targets"
    )
  }
  target_features
}

sb_validate_explicit_feature_map <- function(feature_map) {
  required <- c("source_feature", "target_feature")
  if (!is.data.frame(feature_map) || nrow(feature_map) < 1L) {
    sb_abort(
      "`feature_map` must be a non-empty data frame.",
      "stormbreaker_bad_feature_map"
    )
  }
  positions <- lapply(required, function(column) {
    which(names(feature_map) == column)
  })
  if (any(lengths(positions) != 1L)) {
    sb_abort(
      paste0(
        "`feature_map` must contain exactly one `source_feature` and one ",
        "`target_feature` column."
      ),
      "stormbreaker_bad_feature_map"
    )
  }
  source <- feature_map[[positions[[1L]]]]
  target <- feature_map[[positions[[2L]]]]
  valid_values <- is.character(source) && is.character(target) &&
    !anyNA(source) && !anyNA(target) &&
    all(nzchar(trimws(source))) && all(nzchar(trimws(target)))
  if (!valid_values) {
    sb_abort(
      "Feature-map source and target identifiers must be non-empty strings.",
      "stormbreaker_bad_feature_map"
    )
  }
  if (anyDuplicated(source) || anyDuplicated(target)) {
    sb_abort(
      paste0(
        "`feature_map` must be one-to-one; resolve one-to-many or many-to-one ",
        "biology before applying it."
      ),
      "stormbreaker_non_bijective_feature_map"
    )
  }
  feature_map
}

#' @rdname sb_linear_signature
#' @export
sb_apply_feature_map <- function(
    new_data,
    id_col,
    feature_map,
    target_features) {
  patient_ids <- sb_extract_final_prediction_ids(new_data, id_col)
  target_features <- sb_validate_feature_map_targets(target_features)
  feature_map <- sb_validate_explicit_feature_map(feature_map)

  source <- feature_map$source_feature
  target <- feature_map$target_feature
  if (id_col %in% c(source, target)) {
    sb_abort(
      "The patient ID column cannot be a mapped source or target feature.",
      "stormbreaker_overlapping_feature_map_roles"
    )
  }
  missing_targets <- setdiff(target_features, target)
  extra_targets <- setdiff(target, target_features)
  if (length(missing_targets) > 0L || length(extra_targets) > 0L) {
    sb_abort(
      paste0(
        "The feature-map targets must exactly match `target_features`; ",
        "missing: ",
        if (length(missing_targets) == 0L) "none" else {
          paste(missing_targets, collapse = ", ")
        },
        "; extra: ",
        if (length(extra_targets) == 0L) "none" else {
          paste(extra_targets, collapse = ", ")
        },
        "."
      ),
      "stormbreaker_feature_map_target_mismatch"
    )
  }

  order <- match(target_features, target)
  ordered_map <- feature_map[order, , drop = FALSE]
  source <- ordered_map$source_feature
  source_counts <- vapply(source, function(feature) {
    sum(names(new_data) == feature)
  }, integer(1L))
  if (any(source_counts == 0L)) {
    sb_abort(
      paste0(
        "Prediction data are missing mapped source features: ",
        paste(source[source_counts == 0L], collapse = ", "),
        "."
      ),
      "stormbreaker_missing_feature_map_source"
    )
  }
  if (any(source_counts > 1L)) {
    sb_abort(
      paste0(
        "Prediction data duplicate mapped source features: ",
        paste(source[source_counts > 1L], collapse = ", "),
        "."
      ),
      "stormbreaker_duplicate_feature_map_source"
    )
  }

  source_positions <- vapply(source, function(feature) {
    which(names(new_data) == feature)
  }, integer(1L))
  x <- new_data[source_positions]
  names(x) <- target_features
  x <- sb_validate_predictor_frame(x, "mapped source features")
  id_position <- which(names(new_data) == id_col)
  mapped_data <- cbind(new_data[id_position], x)
  names(mapped_data) <- c(id_col, target_features)
  rownames(mapped_data) <- NULL
  rownames(ordered_map) <- NULL

  structure(
    list(
      status = "success",
      method = "explicit_one_to_one_feature_mapping",
      id_col = id_col,
      n = as.integer(length(patient_ids)),
      mapping = ordered_map,
      data = mapped_data
    ),
    class = "stormbreaker_mapped_features"
  )
}

#' @export
print.stormbreaker_mapped_features <- function(x, ...) {
  cat("<stormbreaker_mapped_features:", x$status, ">\n")
  cat("explicit one-to-one feature mapping\n")
  cat("Patients:", x$n, "| Mapped features:", nrow(x$mapping), "\n")
  cat("No feature was inferred, aggregated, imputed, or filled with zero.\n")
  invisible(x)
}

# ---- linear-signature ----

sb_signature_optional_text <- function(value, argument) {
  if (is.null(value)) {
    return(NULL)
  }
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    sb_abort(
      paste0(argument, " must be NULL or one non-empty string."),
      "stormbreaker_bad_signature_metadata"
    )
  }
  value
}

#' @rdname sb_linear_signature
#' @export
sb_linear_signature <- function(
    signature_id,
    task_type,
    coefficients,
    intercept = NULL,
    source = NULL,
    input_description = NULL,
    prediction_horizon = NULL,
    time_unit = NULL,
    baseline_survival = NULL) {
  if (!is.character(signature_id) || length(signature_id) != 1L ||
      is.na(signature_id) || !nzchar(trimws(signature_id))) {
    sb_abort(
      "signature_id must be one non-empty string.",
      "stormbreaker_bad_signature_id"
    )
  }
  if (!is.character(task_type) || length(task_type) != 1L ||
      is.na(task_type) || !task_type %in% c("binary", "survival")) {
    sb_abort(
      "task_type must be exactly \"binary\" or \"survival\".",
      "stormbreaker_bad_signature_task"
    )
  }
  coefficient_names <- names(coefficients)
  valid_coefficients <- is.numeric(coefficients) &&
    is.null(dim(coefficients)) && length(coefficients) > 0L &&
    !anyNA(coefficients) && all(is.finite(coefficients)) &&
    !is.null(coefficient_names) && !anyNA(coefficient_names) &&
    all(nzchar(trimws(coefficient_names))) &&
    anyDuplicated(coefficient_names) == 0L
  if (!valid_coefficients) {
    sb_abort(
      paste0(
        "coefficients must be a non-empty finite numeric vector named by ",
        "unique, non-empty predictor columns."
      ),
      "stormbreaker_bad_signature_coefficients"
    )
  }
  if (is.null(intercept)) {
    if (identical(task_type, "binary")) {
      sb_abort(
        paste0(
          "A binary logistic signature requires an explicitly stated ",
          "intercept; use intercept = 0 only when zero is the stated value."
        ),
        "stormbreaker_missing_signature_intercept"
      )
    }
    intercept <- 0
  }
  if (!is.numeric(intercept) || length(intercept) != 1L ||
      is.na(intercept) || !is.finite(intercept)) {
    sb_abort(
      "intercept must be one finite number.",
      "stormbreaker_bad_signature_intercept"
    )
  }
  source <- sb_signature_optional_text(source, "source")
  input_description <- sb_signature_optional_text(
    input_description,
    "input_description"
  )

  horizon_fields <- c(
    !is.null(prediction_horizon),
    !is.null(time_unit),
    !is.null(baseline_survival)
  )
  if (identical(task_type, "binary") && any(horizon_fields)) {
    sb_abort(
      "Survival-horizon fields cannot be attached to a binary signature.",
      "stormbreaker_bad_signature_horizon"
    )
  }
  if (identical(task_type, "survival") &&
      any(horizon_fields) && !all(horizon_fields)) {
    sb_abort(
      paste0(
        "A survival horizon requires prediction_horizon, time_unit, ",
        "and baseline_survival together."
      ),
      "stormbreaker_bad_signature_horizon"
    )
  }
  if (all(horizon_fields)) {
    prediction_horizon <- sb_validate_endpoint_horizon(prediction_horizon)
    time_unit <- sb_validate_endpoint_time_unit(time_unit)
    valid_baseline <- is.numeric(baseline_survival) &&
      length(baseline_survival) == 1L && !is.na(baseline_survival) &&
      is.finite(baseline_survival) && baseline_survival >= 0 &&
      baseline_survival <= 1
    if (!valid_baseline) {
      sb_abort(
        "baseline_survival must be one finite number between 0 and 1.",
        "stormbreaker_bad_signature_baseline_survival"
      )
    }
    baseline_survival <- as.numeric(baseline_survival)
  }

  prediction <- if (identical(task_type, "binary")) {
    list(
      output = "probability",
      positive_class = 1L,
      link = "logit"
    )
  } else {
    list(
      output = "risk",
      direction = "higher_linear_score_means_higher_risk"
    )
  }
  if (all(horizon_fields)) {
    prediction <- list(
      output = c("risk", "survival_probability"),
      direction = "higher_linear_score_means_higher_risk",
      prediction_horizon = prediction_horizon,
      time_unit = time_unit,
      baseline_survival = baseline_survival,
      formula = "baseline_survival^exp(linear_score)"
    )
  }
  result <- list(
    signature_id = signature_id,
    pipeline_id = paste0("signature__", signature_id),
    task_type = task_type,
    learner = "locked_linear_signature",
    features = coefficient_names,
    coefficients = coefficients,
    intercept = as.numeric(intercept),
    prediction = prediction,
    source = source,
    input_description = input_description
  )
  if (all(horizon_fields)) {
    result <- append(
      result,
      list(
        prediction_horizon = prediction_horizon,
        time_unit = time_unit,
        baseline_survival = baseline_survival
      ),
      after = match("intercept", names(result))
    )
  }
  structure(result, class = "stormbreaker_linear_signature")
}

sb_validate_linear_signature <- function(signature) {
  if (!inherits(signature, "stormbreaker_linear_signature") ||
      !is.list(signature)) {
    sb_abort(
      "signature must be returned by sb_linear_signature().",
      "stormbreaker_bad_linear_signature"
    )
  }
  sb_linear_signature(
    signature_id = signature[["signature_id"]],
    task_type = signature[["task_type"]],
    coefficients = signature[["coefficients"]],
    intercept = signature[["intercept"]],
    source = signature[["source"]],
    input_description = signature[["input_description"]],
    prediction_horizon = signature[["prediction_horizon"]],
    time_unit = signature[["time_unit"]],
    baseline_survival = signature[["baseline_survival"]]
  )
}

#' @export
print.stormbreaker_linear_signature <- function(x, ...) {
  cat("<stormbreaker_linear_signature>\n")
  cat("Signature:", x$signature_id, "| Task:", x$task_type, "\n")
  cat("Required numeric features:", length(x$features), "\n")
  if ("prediction_horizon" %in% names(x)) {
    cat(
      "Prediction horizon:", format(x$prediction_horizon),
      "[", x$time_unit, "]\n"
    )
  }
  if (!is.null(x$source)) {
    cat("Source metadata:", x$source, "\n")
  }
  invisible(x)
}

sb_signature_safe_metadata <- function(signature) {
  signature_id <- if (is.list(signature) &&
      sb_external_scalar_character(signature[["signature_id"]])) {
    signature[["signature_id"]]
  } else {
    NA_character_
  }
  task_type <- if (is.list(signature) &&
      sb_external_scalar_character(signature[["task_type"]]) &&
      signature[["task_type"]] %in% c("binary", "survival")) {
    signature[["task_type"]]
  } else {
    NA_character_
  }
  list(
    signature_id = signature_id,
    pipeline_id = if (is.na(signature_id)) {
      NA_character_
    } else {
      paste0("signature__", signature_id)
    },
    task_type = task_type,
    learner = "locked_linear_signature",
    prediction_horizon = if (is.list(signature)) {
      signature[["prediction_horizon"]]
    } else {
      NULL
    },
    time_unit = if (is.list(signature)) signature[["time_unit"]] else NULL
  )
}

sb_signature_prediction_rows <- function(
    metadata,
    patient_ids,
    status,
    error_class,
    error_message,
    values = NULL) {
  n <- length(patient_ids)
  rows <- data.frame(
    patient_id = patient_ids,
    pipeline_id = rep(metadata$pipeline_id, n),
    learner = rep(metadata$learner, n),
    stringsAsFactors = FALSE
  )
  if (identical(metadata$task_type, "binary")) {
    rows$positive_class <- rep.int(1L, n)
    rows$probability <- if (is.null(values)) {
      rep(NA_real_, n)
    } else {
      values$probability
    }
  } else if (identical(metadata$task_type, "survival")) {
    rows$risk <- if (is.null(values)) rep(NA_real_, n) else values$risk
    if (!is.null(metadata$prediction_horizon)) {
      rows$prediction_horizon <- rep(metadata$prediction_horizon, n)
      rows$time_unit <- rep(metadata$time_unit, n)
      rows$survival_probability <- if (is.null(values)) {
        rep(NA_real_, n)
      } else {
        values$survival_probability
      }
    }
  }
  rows$status <- rep(status, n)
  rows$error_class <- rep(error_class, n)
  rows$error_message <- rep(error_message, n)
  rows
}

sb_signature_prediction_failure <- function(
    signature,
    id_col,
    patient_ids,
    stage,
    error) {
  metadata <- sb_signature_safe_metadata(signature)
  error_class <- class(error)[[1L]]
  error_message <- conditionMessage(error)
  result <- list(
    status = "failed",
    task_type = metadata$task_type,
    pipeline_id = metadata$pipeline_id,
    learner = metadata$learner,
    id_col = if (is.character(id_col) && length(id_col) == 1L &&
        !is.na(id_col)) {
      id_col
    } else {
      NA_character_
    },
    signature = signature,
    predictions = sb_signature_prediction_rows(
      metadata,
      patient_ids,
      "failed",
      error_class,
      error_message
    ),
    error_stage = stage,
    error_class = error_class,
    error_message = error_message
  )
  structure(result, class = "stormbreaker_final_prediction")
}

sb_select_signature_prediction_features <- function(
    new_data,
    id_col,
    features) {
  if (id_col %in% features) {
    sb_abort(
      "The prediction ID column cannot also be a signature feature.",
      "stormbreaker_overlapping_signature_roles"
    )
  }
  counts <- vapply(features, function(feature) {
    sum(names(new_data) == feature)
  }, integer(1L))
  if (any(counts == 0L)) {
    sb_abort(
      paste0(
        "Prediction data are missing signature features: ",
        paste(features[counts == 0L], collapse = ", "),
        "."
      ),
      "stormbreaker_missing_signature_features"
    )
  }
  if (any(counts > 1L)) {
    sb_abort(
      paste0(
        "Prediction data duplicate signature features: ",
        paste(features[counts > 1L], collapse = ", "),
        "."
      ),
      "stormbreaker_duplicate_signature_feature_columns"
    )
  }
  positions <- vapply(features, function(feature) {
    which(names(new_data) == feature)
  }, integer(1L))
  x <- new_data[positions]
  names(x) <- features
  sb_validate_predictor_frame(x, "new_data signature features")
}

sb_signature_linear_score <- function(signature, x) {
  score <- rep(signature$intercept, nrow(x))
  for (feature in signature$features) {
    score <- score + x[[feature]] * signature$coefficients[[feature]]
  }
  if (anyNA(score) || any(!is.finite(score))) {
    sb_abort(
      "The fixed signature formula produced a non-finite linear score.",
      "stormbreaker_nonfinite_signature_score"
    )
  }
  unname(score)
}

sb_signature_prediction_values <- function(signature, x) {
  score <- sb_signature_linear_score(signature, x)
  if (identical(signature$task_type, "binary")) {
    probability <- stats::plogis(score)
    return(list(probability = unname(probability)))
  }
  values <- list(risk = score)
  if ("prediction_horizon" %in% names(signature)) {
    baseline <- signature$baseline_survival
    probability <- if (identical(baseline, 0)) {
      rep(0, length(score))
    } else {
      baseline^exp(score)
    }
    valid <- !anyNA(probability) && all(is.finite(probability)) &&
      all(probability >= 0 & probability <= 1)
    if (!valid) {
      sb_abort(
        "The fixed signature formula produced an invalid survival probability.",
        "stormbreaker_bad_signature_probability"
      )
    }
    values$survival_probability <- unname(probability)
  }
  values
}

sb_signature_prediction_provenance <- function() {
  list(
    package = "stormbreaker",
    package_version = as.character(utils::packageVersion("stormbreaker")),
    r_version = R.version.string,
    method = "locked_linear_signature_without_fitting"
  )
}

#' @rdname sb_linear_signature
#' @export
sb_predict_signature <- function(signature, new_data, id_col) {
  id_attempt <- sb_final_fit_attempt(
    sb_extract_final_prediction_ids(new_data, id_col)
  )
  if (!is.null(id_attempt$error)) {
    return(sb_signature_prediction_failure(
      signature,
      id_col,
      character(),
      "ids",
      id_attempt$error
    ))
  }
  patient_ids <- id_attempt$value

  signature_attempt <- sb_final_fit_attempt(
    sb_validate_linear_signature(signature)
  )
  if (!is.null(signature_attempt$error)) {
    return(sb_signature_prediction_failure(
      signature,
      id_col,
      patient_ids,
      "signature",
      signature_attempt$error
    ))
  }
  signature <- signature_attempt$value
  metadata <- sb_signature_safe_metadata(signature)

  feature_attempt <- sb_final_fit_attempt(
    sb_select_signature_prediction_features(
      new_data,
      id_col,
      signature$features
    )
  )
  if (!is.null(feature_attempt$error)) {
    return(sb_signature_prediction_failure(
      signature,
      id_col,
      patient_ids,
      "features",
      feature_attempt$error
    ))
  }
  value_attempt <- sb_final_fit_attempt(
    sb_signature_prediction_values(signature, feature_attempt$value)
  )
  if (!is.null(value_attempt$error)) {
    return(sb_signature_prediction_failure(
      signature,
      id_col,
      patient_ids,
      "formula",
      value_attempt$error
    ))
  }

  result <- list(
    status = "success",
    task_type = signature$task_type,
    pipeline_id = signature$pipeline_id,
    learner = signature$learner,
    id_col = id_col,
    signature = signature,
    predictions = sb_signature_prediction_rows(
      metadata,
      patient_ids,
      "success",
      NA_character_,
      NA_character_,
      value_attempt$value
    ),
    provenance = sb_signature_prediction_provenance(),
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )
  structure(result, class = "stormbreaker_final_prediction")
}

# ---- signature-inventory ----

#' @rdname sb_linear_signature
#' @export
sb_signature_inventory <- function(records) {
  record_names <- if (is.list(records)) names(records) else NULL
  core_signature_fields <- c("signature_id", "task_type", "coefficients")
  looks_like_one_raw_record <- is.list(records) &&
    !is.data.frame(records) &&
    !is.null(record_names) &&
    all(core_signature_fields %in% record_names) &&
    is.character(records[["signature_id"]]) &&
    length(records[["signature_id"]]) == 1L &&
    is.character(records[["task_type"]]) &&
    length(records[["task_type"]]) == 1L &&
    is.numeric(records[["coefficients"]]) &&
    is.null(dim(records[["coefficients"]]))
  looks_like_one_record <- inherits(
    records,
    "stormbreaker_linear_signature"
  ) || looks_like_one_raw_record
  if (looks_like_one_record) {
    sb_abort(
      paste0(
        "records must be an outer named list of records; wrap one ",
        "specification as list(record_id = specification)."
      ),
      "stormbreaker_bad_signature_inventory"
    )
  }
  valid <- is.list(records) && !is.data.frame(records) &&
    length(records) > 0L &&
    !is.null(record_names) && !anyNA(record_names) &&
    all(nzchar(trimws(record_names))) &&
    anyDuplicated(record_names) == 0L
  if (!valid) {
    sb_abort(
      paste0(
        "records must be a non-empty list named by unique, non-empty ",
        "record IDs."
      ),
      "stormbreaker_bad_signature_inventory"
    )
  }

  rows <- lapply(seq_along(records), function(index) {
    report <- sb_lint_signature(records[[index]])
    data.frame(
      record_id = names(records)[[index]],
      report,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

# ---- signature-linter ----

sb_signature_from_lint_spec <- function(spec) {
  if (inherits(spec, "stormbreaker_linear_signature")) {
    return(sb_validate_linear_signature(spec))
  }
  if (!is.list(spec)) {
    sb_abort(
      paste0(
        "spec must be one stormbreaker_linear_signature or one named ",
        "list of arguments for sb_linear_signature()."
      ),
      "stormbreaker_bad_signature_spec"
    )
  }
  if (length(spec) > 0L && (is.null(names(spec)) ||
      anyNA(names(spec)) || any(!nzchar(names(spec))) ||
      anyDuplicated(names(spec)) > 0L)) {
    sb_abort(
      "A signature argument list must have unique, non-empty names.",
      "stormbreaker_bad_signature_spec"
    )
  }
  unknown <- setdiff(names(spec), names(formals(sb_linear_signature)))
  if (length(unknown) > 0L) {
    sb_abort(
      paste0(
        "Unknown sb_linear_signature() argument",
        if (length(unknown) > 1L) "s" else "",
        ": ",
        paste0(unknown, collapse = ", "),
        "."
      ),
      "stormbreaker_unknown_signature_argument"
    )
  }
  required <- c("signature_id", "task_type", "coefficients")
  missing <- setdiff(required, names(spec))
  if (length(missing) > 0L) {
    sb_abort(
      paste0(
        "Missing required sb_linear_signature() argument",
        if (length(missing) > 1L) "s" else "",
        ": ",
        paste0(missing, collapse = ", "),
        "."
      ),
      "stormbreaker_missing_signature_argument"
    )
  }
  do.call(sb_linear_signature, spec)
}

sb_is_expected_signature_lint_error <- function(error) {
  inherits(error, c(
    "stormbreaker_bad_signature_spec",
    "stormbreaker_unknown_signature_argument",
    "stormbreaker_missing_signature_argument",
    "stormbreaker_bad_linear_signature",
    "stormbreaker_bad_signature_id",
    "stormbreaker_bad_signature_task",
    "stormbreaker_bad_signature_coefficients",
    "stormbreaker_missing_signature_intercept",
    "stormbreaker_bad_signature_intercept",
    "stormbreaker_bad_signature_metadata",
    "stormbreaker_bad_signature_horizon",
    "stormbreaker_bad_signature_baseline_survival",
    "stormbreaker_bad_prediction_horizon",
    "stormbreaker_bad_time_unit"
  ))
}

sb_signature_lint_report <- function(status, spec, signature = NULL,
                                     error = NULL) {
  source_object <- if (is.list(signature)) signature else spec
  safe_character <- function(name) {
    if (!is.list(source_object) || !name %in% names(source_object)) {
      return(NA_character_)
    }
    value <- source_object[[name]]
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(trimws(value))) {
      return(NA_character_)
    }
    value
  }
  coefficients <- if (is.list(source_object)) {
    source_object[["coefficients"]]
  } else {
    NULL
  }
  n_features <- if (is.numeric(coefficients) && is.null(dim(coefficients)) &&
      length(coefficients) > 0L) {
    as.integer(length(coefficients))
  } else {
    NA_integer_
  }
  output <- if (is.list(signature) && is.list(signature[["prediction"]])) {
    paste(signature[["prediction"]][["output"]], collapse = "+")
  } else {
    NA_character_
  }

  data.frame(
    status = status,
    signature_id = safe_character("signature_id"),
    task_type = safe_character("task_type"),
    n_features = n_features,
    output = output,
    source = safe_character("source"),
    issue_code = if (is.null(error)) NA_character_ else class(error)[[1L]],
    message = if (is.null(error)) NA_character_ else conditionMessage(error),
    stringsAsFactors = FALSE
  )
}

#' @rdname sb_linear_signature
#' @export
sb_lint_signature <- function(spec) {
  signature <- NULL
  result <- tryCatch(
    {
      signature <- sb_signature_from_lint_spec(spec)
      signature
    },
    stormbreaker_error = function(error) {
      if (!sb_is_expected_signature_lint_error(error)) {
        stop(error)
      }
      error
    }
  )

  if (inherits(result, "stormbreaker_error")) {
    return(sb_signature_lint_report(
      status = "rejected",
      spec = spec,
      signature = signature,
      error = result
    ))
  }
  sb_signature_lint_report(
    status = "admitted",
    spec = spec,
    signature = result
  )
}

# ---- constant-sum-representation ----

# Same training-only coordinates used by the formal-run adapters. This removes
# the known constant-sum direction, not low-variance principal components.
sb_constant_sum_projection <- function(x) {
  matrix_x <- as.matrix(x)
  p <- ncol(matrix_x)
  sums <- rowSums(matrix_x)
  projection <- list(features = names(x), basis = NULL, center = rep(0, p), row_sum = mean(sums))
  if (p > 1L && max(abs(sums - mean(sums))) < 1e-10) {
    basis <- stats::contr.helmert(p)
    projection$basis <- sweep(basis, 2L, sqrt(colSums(basis * basis)), "/")
    projection$center <- colMeans(matrix_x)
    reconstructed <- sweep(as.matrix(sb_constant_sum_project(x, projection)) %*%
      t(projection$basis), 2L, projection$center, "+")
    if (max(abs(reconstructed - matrix_x)) >= 1e-8) {
      sb_abort("Constant-sum coordinates did not reconstruct training inputs.",
        "stormbreaker_constant_sum_mismatch")
    }
  }
  projection
}

sb_constant_sum_project <- function(x, projection) {
  matrix_x <- as.matrix(x[, projection$features, drop = FALSE])
  if (is.null(projection$basis)) return(as.data.frame(matrix_x))
  if (max(abs(rowSums(matrix_x) - projection$row_sum)) > 1e-8) {
    sb_abort("Prediction input does not obey the fitted constant-sum representation.",
      "stormbreaker_constant_sum_mismatch")
  }
  transformed <- sweep(matrix_x, 2L, projection$center, "-") %*% projection$basis
  colnames(transformed) <- paste0("contrast_", seq_len(ncol(transformed)))
  as.data.frame(transformed)
}

sb_is_constant_sum_state <- function(state) {
  if (!is.list(state) || is.data.frame(state)) return(FALSE)
  identical(state$representation, "constant_sum_orthonormal_v1") ||
    identical(state$round2_representation, "constant_sum_orthonormal_v1") ||
    identical(state$variant, "brglm2_constant_sum_v1")
}

sb_constant_sum_learner_adapter <- function(adapter, constructor) {
  id <- adapter$id
  if (!id %in% c("brglm2_logistic", "lda_binary", "hda_binary", "fda_binary", "mda_binary", "loclda_binary")) {
    return(adapter)
  }
  native <- adapter
  make_native <- function(params) constructor(id, params)
  adapter$fit <- function(x, outcome) {
    projection <- sb_constant_sum_projection(x)
    if (is.null(projection$basis)) return(native$fit(x, outcome))
    inner_params <- native$params
    if (id == "loclda_binary") {
      inner_params$k_max <- sb_loclda_neighbors(inner_params, ncol(x), nrow(x), length(unique(outcome)))
    }
    inner <- make_native(inner_params)$fit(sb_constant_sum_project(x, projection), outcome)
    structure(list(learner = id, features = names(x), params = native$params,
      projection = projection, inner_params = inner_params, inner = inner),
      class = "stormbreaker_constant_sum_model")
  }
  adapter$predict <- function(model, x) {
    if (!inherits(model, "stormbreaker_constant_sum_model")) return(native$predict(model, x))
    make_native(model$inner_params)$predict(model$inner, sb_constant_sum_project(x, model$projection))
  }
  adapter$extract_state <- function(model, context = list()) {
    if (!inherits(model, "stormbreaker_constant_sum_model")) return(native$extract_state(model, context))
    structure(list(learner = id, features = model$features, params = model$params,
      representation = "constant_sum_orthonormal_v1", projection = model$projection,
      inner_params = model$inner_params,
      inner = make_native(model$inner_params)$extract_state(model$inner, context)),
      class = c("stormbreaker_constant_sum_state", "stormbreaker_learner_state"))
  }
  adapter$predict_state <- function(state, x) {
    if (!sb_is_constant_sum_state(state)) return(native$predict_state(state, x))
    params <- if (is.null(state$inner_params)) native$params else state$inner_params
    if (id == "hda_binary") loadNamespace("hda")
    make_native(params)$predict_state(sb_upgrade_legacy_names(state$inner),
      sb_constant_sum_project(x, state$projection))
  }
  if (!is.null(native$extract_coefficients)) {
    adapter$extract_coefficients <- function(model) {
      if (!inherits(model, "stormbreaker_constant_sum_model")) return(native$extract_coefficients(model))
      coefficients <- make_native(model$inner_params)$extract_coefficients(model$inner)
      beta <- as.vector(model$projection$basis %*% coefficients$coefficient[-1L])
      data.frame(term = c("(Intercept)", model$features),
        coefficient = c(coefficients$coefficient[[1L]] - sum(beta * model$projection$center), beta))
    }
  }
  adapter
}

sb_constant_sum_linear_parameters <- function(state) {
  inner <- sb_upgrade_legacy_names(state$inner)
  parameters <- if (inner$learner == "lda_binary") sb_discriminant_parameters(inner) else
    sb_extended_logistic_parameters(inner)
  if (is.null(state$projection$basis)) return(parameters)
  beta <- as.vector(state$projection$basis %*% parameters$coefficient)
  parameters$coefficient <- stats::setNames(beta, state$projection$features)
  parameters$intercept <- parameters$intercept - sum(beta * state$projection$center)
  if (!is.null(parameters$means)) {
    parameters$means <- sweep(parameters$means %*% t(state$projection$basis),
      2L, state$projection$center, "+")
    colnames(parameters$means) <- state$projection$features
  }
  if (!is.null(parameters$quadratic)) {
    parameters$quadratic <- matrix(0, length(beta), length(beta),
      dimnames = list(state$projection$features, state$projection$features))
  }
  parameters
}
