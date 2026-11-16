#' Data and endpoints
#' @name sb_binary_task
#' @description Prepare task data, join expression and clinical records, define endpoints and generate synthetic examples.
#' @rawRd \arguments{
#' \item{data}{
#' \strong{For \code{sb_binary_task()}, \code{sb_surv_task()}:} A data frame with one row per observation. Without \code{group_col},
#' each row is treated as an independent patient.
#' 
#' \strong{For \code{sb_multiclass_task()}:} Data frame, one row per observation.
#' 
#' \strong{For \code{sb_binary_endpoint()}:} A data frame containing at least the ID and outcome columns.
#' Unrequested columns are ignored and are not stored.
#' 
#' \strong{For \code{sb_surv_endpoint()}:} A data frame containing at least the ID, time, and event columns.
#' Unrequested columns are ignored and are not stored.
#' }
#' \item{id_col}{
#' \strong{For \code{sb_binary_task()}, \code{sb_surv_task()}:} Name of the unique observation identifier column.
#' 
#' \strong{For \code{sb_multiclass_task()}:} Unique observation identifier column.
#' 
#' \strong{For \code{sb_binary_endpoint()}, \code{sb_surv_endpoint()}:} Name of the unique patient identifier column.
#' 
#' \strong{For \code{sb_read_cohort()}:} ID column present in both files; defaults to \code{"patient_id"}.
#' }
#' \item{outcome_col}{
#' \strong{For \code{sb_binary_task()}:} Name of the numeric or logical outcome column, coded as
#' \code{0} or \code{1}.
#' 
#' \strong{For \code{sb_multiclass_task()}:} Column containing three or more observed classes.
#' 
#' \strong{For \code{sb_binary_endpoint()}:} Name of the numeric or logical \code{0/1} outcome column.
#' }
#' \item{predictors}{
#' \strong{For \code{sb_binary_task()}:} Character vector naming finite numeric predictor columns.
#' The reserved model term \code{"(Intercept)"} cannot be a predictor name.
#' 
#' \strong{For \code{sb_surv_task()}:} Character vector naming numeric predictor columns.
#' 
#' \strong{For \code{sb_multiclass_task()}:} Names of finite numeric predictor columns.
#' }
#' \item{decision_thresholds}{
#' \strong{For \code{sb_binary_task()}:} Optional unique, strictly increasing finite
#' probability thresholds between zero and one. These should represent
#' clinically meaningful action thresholds declared before model evaluation.
#' No default grid is supplied.
#' 
#' \strong{For \code{sb_surv_task()}:} Optional unique, strictly increasing finite
#' event-risk probabilities between zero and one. A prediction horizon is
#' required. These should be clinically meaningful action thresholds
#' declared before evaluation; no default grid is supplied.
#' }
#' \item{group_col}{
#' \strong{For \code{sb_binary_task()}:} Optional patient/group identifier column for repeated
#' observations. Row IDs remain unique. Generated evaluation folds, glmnet
#' learner tuning, and e1071 SVM probability calibration keep each group
#' together; supplied folds must do so too.
#' This does not change metric weighting or bootstrap sampling, and does not
#' configure resampling internal to other learner or selector engines.
#' 
#' \strong{For \code{sb_surv_task()}:} Optional patient/group identifier column for repeated
#' observations. Row IDs remain unique. Generated evaluation folds and glmnet
#' learner tuning keep each group together; supplied folds must do so too.
#' This does not change metric weighting or bootstrap sampling, and does not
#' configure resampling internal to other learner or selector engines.
#' 
#' \strong{For \code{sb_multiclass_task()}:} Optional patient/group column. Resampling keeps all rows
#' from a group together, including inner tuning folds.
#' }
#' \item{time_col}{
#' \strong{For \code{sb_surv_task()}:} Name of the non-negative survival-time column.
#' 
#' \strong{For \code{sb_surv_endpoint()}:} Name of the finite non-negative follow-up column.
#' }
#' \item{event_col}{
#' \strong{For \code{sb_surv_task()}:} Name of the event column, coded as \code{0} or \code{1}, where \code{1}
#' means the event was observed.
#' 
#' \strong{For \code{sb_surv_endpoint()}:} Name of the numeric or logical \code{0/1} event column, where
#' \code{1} means the event was observed.
#' }
#' \item{time_unit}{
#' Explicit unit label: \code{"day"}, \code{"month"}, or \code{"year"}.
#' }
#' \item{prediction_horizon}{
#' \strong{For \code{sb_surv_task()}:} Optional single positive time point for survival
#' probability prediction and IPCW Brier scoring. Its unit is \code{time_unit}.
#' 
#' \strong{For \code{sb_surv_endpoint()}:} Optional finite positive evaluation horizon in
#' \code{time_unit}. A final prediction with a saved horizon controls scoring.
#' }
#' \item{classes}{
#' Explicit class order. By default use observed factor levels,
#' or sorted character labels. No ordinal relationship is assumed.
#' }
#' \item{expression_file}{
#' Path to a UTF-8 comma-separated file with one row
#' per sample, an ID column and named gene columns. A combined file may also
#' contain clinical columns; assign predictor roles explicitly in the task.
#' }
#' \item{clinical_file}{
#' Optional path to a UTF-8 clinical CSV with exactly the
#' same unique sample IDs. NULL reads only the first file, suitable for
#' prediction without outcomes.
#' }
#' \item{n}{
#' Number of synthetic patients. Must be at least 10.
#' }
#' \item{seed}{
#' Positive integer controlling synthetic generation or the random split.
#' }
#' \item{task}{
#' For \code{sb_split_binary()}, a binary task with declared column roles.
#' }
#' \item{train_prop}{
#' Target proportion of independent patients or groups assigned to development,
#' strictly between zero and one; defaults to 0.7. Approximate with small strata.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_binary_task()}}{
#' \code{sb_binary_task()} records an explicit patient-level \code{0/1} outcome. Class
#' \code{1} is the positive class. The function never drops patients, imputes
#' values, converts labels, or guesses the positive class. Optional decision
#' thresholds are supplied by the researcher; the package does not derive or
#' optimize them from outcomes or model performance.
#' 
#' \strong{Returns:} A \code{stormbreaker_binary_task} object preserving the input rows and
#' recording outcome and predictor roles.
#' }
#' 
#' \subsection{\code{sb_split_binary()}}{
#' Split a binary cohort into development and an internal holdout. Strata are
#' based on group outcome profiles (only 0, only 1, or mixed). All records from
#' an individual or group stay on the same side. Both sides must contain both
#' classes; percentages can differ from the target after rounding or because
#' group sizes differ. This is not balancing or independent external validation.
#' \strong{Returns:} A list with a \code{development} task, an \code{assessment}
#' data frame, row-level \code{assignment}, aggregate \code{counts}, and split settings.
#' Passing the unchanged \code{assessment} frame to \code{finish_run()} labels
#' its plots as an internal holdout. Keep it out of screening and model selection.
#' }
#' 
#' \subsection{\code{sb_surv_task()}}{
#' \code{sb_surv_task()} is STORMBREAKER's data contract. It records column roles
#' and rejects ambiguous or incomplete endpoint data. The function never
#' removes patients, imputes values, converts time units, or chooses a
#' prediction horizon or decision threshold from the observed data.
#' 
#' \strong{Returns:} A \code{stormbreaker_surv_task} object. Its \code{data} element preserves the
#' input rows and its remaining elements record endpoint and predictor roles,
#' plus the prediction horizon and event-risk thresholds when declared.
#' }
#' 
#' \subsection{\code{sb_multiclass_task()}}{
#' Define a nominal multiclass classification task
#' 
#' \strong{Returns:} A \code{stormbreaker_multiclass_task} preserving rows and original labels.
#' }
#' 
#' \subsection{\code{sb_binary_endpoint()}}{
#' \code{sb_binary_endpoint()} records patient IDs and an observed \code{0/1} outcome
#' without accepting or retaining predictors. It is an outcome-only input for
#' \code{\link[=sb_evaluate_external]{sb_evaluate_external()}}, not a modeling task. One patient or one observed
#' class is allowed; metrics that require both classes remain unavailable.
#' 
#' \strong{Returns:} A \code{stormbreaker_binary_endpoint} containing only patient IDs,
#' outcomes, class counts, and endpoint-role metadata.
#' }
#' 
#' \subsection{\code{sb_surv_endpoint()}}{
#' \code{sb_surv_endpoint()} records patient IDs, observed follow-up, and event
#' status without accepting or retaining predictors. It is an outcome-only
#' input for \code{\link[=sb_evaluate_external]{sb_evaluate_external()}}, not a modeling task. One patient, all
#' events, or all censoring is allowed; metrics report when the endpoint does
#' not contain enough evaluable information.
#' 
#' \strong{Returns:} A \code{stormbreaker_surv_endpoint} containing only patient IDs, time,
#' event, counts, and time-semantics metadata.
#' }
#' 
#' \subsection{\code{sb_read_cohort()}}{
#' Read a sample-by-gene expression CSV and optionally join a clinical CSV.
#' IDs are preserved as text, including leading zeroes. Clinical rows are
#' aligned to expression order by ID, never by row position. Neither missing
#' measurements nor outcome labels are repaired or recoded.
#' 
#' \strong{Returns:} A data frame in expression-file row order. Pass explicit outcome,
#' time and predictor names to \code{\link[=sb_binary_task]{sb_binary_task()}} or \code{\link[=sb_surv_task]{sb_surv_task()}}.
#' 
#' \strong{Details:} Duplicate IDs or column names, missing/extra clinical IDs, and
#' overlapping non-ID columns are errors. Numeric columns are read as
#' numbers, while non-numeric labels remain text. The function does not
#' normalize expression, map probes, impute values or infer the endpoint.
#' Use study-specific processing before these modeling-ready CSV files.
#' }
#' 
#' \subsection{\code{sb_example_binary()}}{
#' This generator is only for examples and automated tests. It creates no
#' files and does not represent a clinical cohort. The outcome is balanced,
#' class \code{1} is positive, and all predictors are numeric.
#' 
#' \strong{Returns:} A data frame containing \code{patient_id}, \code{outcome}, and three numeric
#' predictor columns.
#' }
#' 
#' \subsection{\code{sb_example_survival()}}{
#' This generator is only for examples and automated tests. It creates no
#' files and does not represent a clinical cohort. \code{gene_a} influences the
#' simulated event hazard, while \code{gene_b} and \code{gene_c} are weaker noise
#' predictors. A common administrative censoring cutoff yields enough events
#' and censored patients for the package's five-fold example.
#' 
#' \strong{Returns:} A data frame containing \code{patient_id}, \code{time}, \code{event}, and three
#' numeric predictor columns. Time is expressed in arbitrary month-like
#' units for demonstration.
#' }
#' }
#' @rawRd \references{
#' Vickers AJ, Elkin EB (2006). Decision curve analysis: a novel
#' method for evaluating prediction models. \emph{Medical Decision Making},
#' 26(6), 565-574. \doi{10.1177/0272989X06295361}.
#' 
#' Vickers AJ, Cronin AM, Elkin EB, Gonen M (2008). Extensions to
#' decision curve analysis. \emph{BMC Medical Informatics and Decision Making},
#' 8, 53. \doi{10.1186/1472-6947-8-53}.
#' }
#' @rawRd \examples{
#' # sb_binary_task
#' d <- sb_example_binary(n = 30, seed = 1)
#' task <- sb_binary_task(
#'   d, "patient_id", "outcome",
#'   predictors = c("gene_a", "gene_b")
#' )
#' task
#' 
#' # sb_surv_task
#' d <- sb_example_survival(n = 30, seed = 1)
#' task <- sb_surv_task(
#'   d, "patient_id", "time", "event",
#'   predictors = c("gene_a", "gene_b"),
#'   time_unit = "month"
#' )
#' task
#' 
#' # sb_binary_endpoint
#' endpoint <- sb_binary_endpoint(
#'   data.frame(patient_id = c("e1", "e2"), outcome = c(0, 1)),
#'   "patient_id",
#'   "outcome"
#' )
#' endpoint
#' 
#' # sb_surv_endpoint
#' endpoint <- sb_surv_endpoint(
#'   data.frame(
#'     patient_id = c("e1", "e2"),
#'     time = c(3, 5),
#'     event = c(1, 0)
#'   ),
#'   "patient_id",
#'   "time",
#'   "event",
#'   "month"
#' )
#' endpoint
#' 
#' # sb_read_cohort
#' example_file <- system.file("extdata", "demo-binary-development.csv",
#'   package = "stormbreaker")
#' data <- sb_read_cohort(example_file)
#' sb_binary_task(data, "patient_id", "outcome", c("gene_a", "gene_b"))
#' 
#' # sb_example_binary
#' d <- sb_example_binary(n = 50, seed = 1)
#' head(d)
#' 
#' # sb_example_survival
#' d <- sb_example_survival(n = 50, seed = 1)
#' head(d)
#' }
NULL

# ---- task-data ----

# ---- read-cohort ----

sb_read_cohort_csv <- function(path, id_col, label) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !file.exists(path) || dir.exists(path)) {
    sb_abort(paste0("Supply an existing ", label, " CSV file."),
      "stormbreaker_bad_cohort_file")
  }
  data <- tryCatch(utils::read.csv(path, colClasses = "character",
    check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
    na.strings = c("", "NA"), fill = FALSE), error = function(error) {
      sb_abort(paste0("Cannot read ", label, " CSV: ", conditionMessage(error)),
        "stormbreaker_bad_cohort_file")
    })
  if (anyDuplicated(names(data))) {
    sb_abort(paste0(label, " CSV has duplicate columns: ",
      paste(unique(names(data)[duplicated(names(data))]), collapse = ", "), "."),
      "stormbreaker_duplicate_column_names")
  }
  if (!id_col %in% names(data)) {
    sb_abort(paste0(label, " CSV is missing ID column '", id_col, "'."),
      "stormbreaker_missing_endpoint_columns")
  }
  if (!nrow(data) || anyNA(data[[id_col]]) || any(!nzchar(trimws(data[[id_col]])))) {
    sb_abort(paste0(label, " CSV has no samples or contains missing/blank IDs."),
      "stormbreaker_bad_id")
  }
  duplicated_ids <- unique(data[[id_col]][duplicated(data[[id_col]])])
  if (length(duplicated_ids)) {
    sb_abort(paste0(label, " CSV contains duplicate IDs: ",
      paste(utils::head(duplicated_ids, 10L), collapse = ", "), "."),
      "stormbreaker_duplicate_id")
  }
  for (column in setdiff(names(data), id_col)) {
    data[[column]] <- utils::type.convert(data[[column]], as.is = TRUE)
  }
  data
}

#' @rdname sb_binary_task
#' @export
sb_read_cohort <- function(expression_file, clinical_file = NULL,
                           id_col = "patient_id") {
  sb_check_name(id_col, "id_col")
  expression <- sb_read_cohort_csv(expression_file, id_col, "Expression")
  if (is.null(clinical_file)) return(expression)
  clinical <- sb_read_cohort_csv(clinical_file, id_col, "Clinical")
  missing <- setdiff(expression[[id_col]], clinical[[id_col]])
  extra <- setdiff(clinical[[id_col]], expression[[id_col]])
  if (length(missing) || length(extra)) {
    show_ids <- function(x) {
      if (length(x)) paste(utils::head(x, 10L), collapse = ", ") else "none"
    }
    sb_abort(paste0("Sample IDs do not match. Missing clinical IDs (", length(missing),
      "): ", show_ids(missing), "; extra clinical IDs (", length(extra), "): ",
      show_ids(extra), ". No rows were dropped."), "stormbreaker_unmatched_cohort_ids")
  }
  overlap <- intersect(setdiff(names(expression), id_col), setdiff(names(clinical), id_col))
  if (length(overlap)) {
    sb_abort(paste0("Both CSV files contain non-ID columns: ",
      paste(overlap, collapse = ", "), ". Keep each column in one file."),
      "stormbreaker_overlapping_cohort_columns")
  }
  clinical <- clinical[match(expression[[id_col]], clinical[[id_col]]), , drop = FALSE]
  cbind(expression, clinical[setdiff(names(clinical), id_col)])
}

# ---- task-binary ----

sb_valid_decision_thresholds <- function(value) {
  is.numeric(value) && length(value) >= 1L &&
    !anyNA(value) && all(is.finite(value)) &&
    all(value > 0 & value < 1) &&
    !anyDuplicated(value) &&
    (length(value) == 1L || all(diff(value) > 0))
}

sb_validate_decision_thresholds <- function(value) {
  if (!sb_valid_decision_thresholds(value)) {
    sb_abort(
      paste0(
        "`decision_thresholds` must be unique, strictly increasing finite ",
        "probabilities between 0 and 1."
      ),
      "stormbreaker_bad_decision_thresholds"
    )
  }
  unname(as.numeric(value))
}

#' @rdname sb_binary_task
#' @export
sb_binary_task <- function(data, id_col, outcome_col, predictors,
                           decision_thresholds = NULL, group_col = NULL) {
  if (!is.data.frame(data)) {
    sb_abort("`data` must be a data.frame.", "stormbreaker_bad_data")
  }
  if (nrow(data) < 2L) {
    sb_abort(
      "`data` must contain at least two patients.",
      "stormbreaker_too_few_patients"
    )
  }
  if (anyDuplicated(names(data))) {
    sb_abort(
      "The data frame contains duplicate column names.",
      "stormbreaker_duplicate_column_names"
    )
  }

  sb_check_name(id_col, "id_col")
  sb_check_name(outcome_col, "outcome_col")
  if (identical(id_col, outcome_col)) {
    sb_abort(
      "ID and outcome must use different columns.",
      "stormbreaker_overlapping_endpoint_roles"
    )
  }
  required_columns <- c(id_col, outcome_col)
  missing_required <- setdiff(required_columns, names(data))
  if (length(missing_required) > 0L) {
    sb_abort(
      paste0(
        "Missing binary task columns: ",
        paste(missing_required, collapse = ", "), "."
      ),
      "stormbreaker_missing_endpoint_columns"
    )
  }

  id <- data[[id_col]]
  if (anyNA(id) || any(!nzchar(trimws(as.character(id))))) {
    sb_abort(
      "Patient identifiers contain missing or blank values.",
      "stormbreaker_bad_id"
    )
  }
  if (anyDuplicated(id)) {
    sb_abort(
      "Patient identifiers must be unique.",
      "stormbreaker_duplicate_id"
    )
  }

  outcome <- data[[outcome_col]]
  valid_outcome <- (is.numeric(outcome) || is.logical(outcome)) &&
    !anyNA(outcome) &&
    (!is.numeric(outcome) || all(is.finite(outcome))) &&
    all(outcome %in% c(0, 1))
  if (!valid_outcome) {
    sb_abort(
      "`outcome` must contain only finite 0/1 values.",
      "stormbreaker_bad_binary_outcome"
    )
  }
  if (sum(outcome == 1) < 1L) {
    sb_abort(
      "At least one positive patient with outcome=1 is required.",
      "stormbreaker_no_positive_class"
    )
  }
  if (sum(outcome == 0) < 1L) {
    sb_abort(
      "At least one negative patient with outcome=0 is required.",
      "stormbreaker_no_negative_class"
    )
  }

  if (!is.character(predictors) || length(predictors) < 1L ||
      anyNA(predictors) || any(!nzchar(predictors))) {
    sb_abort(
      "`predictors` must contain at least one non-empty column name.",
      "stormbreaker_bad_predictors"
    )
  }
  predictors <- unique(predictors)
  sb_validate_binary_predictor_names(predictors)
  missing_predictors <- setdiff(predictors, names(data))
  if (length(missing_predictors) > 0L) {
    sb_abort(
      paste0(
        "Missing predictor columns: ",
        paste(missing_predictors, collapse = ", "), "."
      ),
      "stormbreaker_missing_predictors"
    )
  }
  if (length(intersect(predictors, required_columns)) > 0L) {
    sb_abort(
      "Predictors cannot reuse the ID or outcome columns.",
      "stormbreaker_reserved_predictor"
    )
  }
  numeric_predictors <- vapply(data[predictors], is.numeric, logical(1))
  if (!all(numeric_predictors)) {
    sb_abort(
      paste0(
        "All minimum-core predictors must be numeric: ",
        paste(
          names(numeric_predictors)[!numeric_predictors],
          collapse = ", "
        ), "."
      ),
      "stormbreaker_bad_predictor_type"
    )
  }
  finite_predictors <- vapply(
    data[predictors],
    function(column) !anyNA(column) && all(is.finite(column)),
    logical(1)
  )
  if (!all(finite_predictors)) {
    sb_abort(
      paste0(
        "Predictors contain missing or non-finite values: ",
        paste(
          names(finite_predictors)[!finite_predictors],
          collapse = ", "
        ), "."
      ),
      "stormbreaker_bad_predictor_values"
    )
  }

  if (!is.null(decision_thresholds)) {
    decision_thresholds <- sb_validate_decision_thresholds(
      decision_thresholds
    )
  }
  group_col <- sb_validate_group_column(
    data, group_col, c(id_col, outcome_col), predictors
  )

  task <- list(
    data = data,
    id_col = id_col,
    outcome_col = outcome_col,
    predictors = predictors,
    positive_class = 1L,
    negative_class = 0L,
    n = nrow(data),
    n_positive = as.integer(sum(outcome == 1)),
    n_negative = as.integer(sum(outcome == 0)),
    independent_unit = "patient"
  )
  if (!is.null(decision_thresholds)) {
    task <- append(
      task,
      list(decision_thresholds = decision_thresholds),
      after = match("negative_class", names(task))
    )
  }
  if (!is.null(group_col)) task$group_col <- group_col
  class(task) <- "stormbreaker_binary_task"
  task
}

sb_validate_binary_task_contract <- function(task) {
  required_metadata_fields <- c(
    "id_col", "outcome_col", "predictors", "positive_class",
    "negative_class", "n", "n_positive", "n_negative",
    "independent_unit"
  )
  required_fields <- c("data", required_metadata_fields)
  if (!is.list(task) || length(setdiff(required_fields, names(task))) > 0L) {
    sb_abort(
      "The binary task is incomplete; recreate it with `sb_binary_task()`.",
      "stormbreaker_bad_task"
    )
  }
  has_decision_thresholds <- "decision_thresholds" %in% names(task)
  decision_thresholds <- if (has_decision_thresholds) {
    task$decision_thresholds
  } else {
    NULL
  }
  validated <- sb_binary_task(
    data = task$data,
    id_col = task$id_col,
    outcome_col = task$outcome_col,
    predictors = task$predictors,
    decision_thresholds = decision_thresholds,
    group_col = task$group_col
  )
  comparison_fields <- required_metadata_fields
  if (has_decision_thresholds) {
    comparison_fields <- c(comparison_fields, "decision_thresholds")
  }
  if ("group_col" %in% names(task)) {
    comparison_fields <- c(comparison_fields, "group_col")
  }
  metadata_ok <- identical(
    unclass(task)[comparison_fields],
    unclass(validated)[comparison_fields]
  )
  if (!metadata_ok) {
    sb_abort(
      "The binary task metadata changed; recreate it with `sb_binary_task()`.",
      "stormbreaker_mutated_task"
    )
  }
  validated
}

#' @export
print.stormbreaker_binary_task <- function(x, ...) {
  cat("<stormbreaker_binary_task>\n")
  cat(
    if (is.null(x$group_col)) "Patients:" else "Observations:", x$n,
    "| Positive (1):", x$n_positive,
    "| Negative (0):", x$n_negative, "\n"
  )
  cat("Outcome:", x$outcome_col, "| Predictors:", length(x$predictors), "\n")
  if (!is.null(x$group_col)) {
    cat("Groups:", length(unique(sb_task_group_ids(x))), "| Column:",
        x$group_col, "\n")
  }
  if ("decision_thresholds" %in% names(x)) {
    cat(
      "Decision thresholds:", length(x$decision_thresholds),
      "| Range:",
      paste(range(x$decision_thresholds), collapse = " to "), "\n"
    )
  }
  invisible(x)
}

# ---- task-multiclass ----

#' @rdname sb_binary_task
#' @export
sb_multiclass_task <- function(data, id_col, outcome_col, predictors,
                               classes = NULL, group_col = NULL) {
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data))) sb_abort(
    "Supply a nonempty data frame with unique column names.", "stormbreaker_bad_data")
  sb_check_name(id_col, "id_col")
  sb_check_name(outcome_col, "outcome_col")
  if (id_col == outcome_col || !all(c(id_col, outcome_col) %in% names(data))) sb_abort(
    "Supply distinct existing ID and outcome columns.", "stormbreaker_bad_multiclass_task")
  id <- as.character(data[[id_col]])
  if (anyNA(id) || any(!nzchar(trimws(id))) || anyDuplicated(id)) sb_abort(
    "Observation IDs must be unique and nonblank.", "stormbreaker_bad_id")
  outcome <- data[[outcome_col]]
  if (!is.atomic(outcome) || !is.null(dim(outcome)) || anyNA(outcome) ||
      any(!nzchar(trimws(as.character(outcome)))) ||
      (is.numeric(outcome) && any(!is.finite(outcome)))) sb_abort(
    "Class labels must be nonmissing and nonblank.", "stormbreaker_bad_multiclass_outcome")
  observed <- unique(as.character(outcome))
  if (is.null(classes)) classes <- if (is.factor(outcome))
    levels(outcome)[levels(outcome) %in% observed] else sort(observed)
  if (!is.character(classes) || length(classes) < 3L || anyNA(classes) ||
      anyDuplicated(classes) || !setequal(classes, observed)) sb_abort(
    "classes must name each observed class exactly once, with at least three classes.",
    "stormbreaker_bad_multiclass_outcome")
  if (!is.character(predictors) || !length(predictors) || anyNA(predictors) ||
      anyDuplicated(predictors) || !all(predictors %in% names(data)) ||
      any(predictors %in% c(id_col, outcome_col))) sb_abort(
    "Supply distinct predictor columns excluding ID and outcome.", "stormbreaker_bad_predictors")
  sb_validate_binary_predictor_names(predictors)
  sb_validate_predictor_frame(data[predictors], "predictors")
  group_col <- sb_validate_group_column(data, group_col, c(id_col, outcome_col), predictors)
  structure(list(data = data, id_col = id_col, outcome_col = outcome_col,
    predictors = predictors, classes = classes, group_col = group_col, n = nrow(data)),
    class = "stormbreaker_multiclass_task")
}

sb_validate_multiclass_task <- function(task) {
  if (!inherits(task, "stormbreaker_multiclass_task")) sb_abort(
    "Create the task with sb_multiclass_task().", "stormbreaker_bad_multiclass_task")
  sb_multiclass_task(task$data, task$id_col, task$outcome_col, task$predictors,
    task$classes, task$group_col)
}

sb_multiclass_folds <- function(task, folds, seed) {
  folds <- sb_check_positive_integer(folds, "folds", minimum = 2L)
  seed <- sb_check_positive_integer(seed, "seed", minimum = 1L)
  y <- factor(task$data[[task$outcome_col]], levels = task$classes)
  groups <- sb_task_group_ids(task)
  if (is.null(groups) || !anyDuplicated(groups)) {
    if (any(table(y) < folds)) sb_abort(
      "Each class needs at least folds observations.", "stormbreaker_infeasible_folds")
    assignment <- integer(task$n)
    withr::with_seed(seed, {
      for (level in levels(y)) {
        rows <- which(y == level)
        assignment[rows[sample.int(length(rows))]] <- rep(seq_len(folds), length.out = length(rows))
      }
    })
  } else {
    group_index <- match(groups, unique(groups))
    counts <- unclass(table(factor(group_index), y))
    if (any(colSums(counts > 0) < folds)) sb_abort(
      "Each class needs at least folds distinct groups.", "stormbreaker_infeasible_folds")
    totals <- colSums(counts)
    assigned <- integer(nrow(counts))
    balance <- matrix(0, folds, ncol(counts))
    withr::with_seed(seed, {
      order_groups <- sample.int(nrow(counts))
      order_groups <- order_groups[order(-rowSums(counts[order_groups, , drop = FALSE] > 0),
        -rowSums(counts[order_groups, , drop = FALSE]))]
      for (g in order_groups) {
        cost <- drop(balance %*% (counts[g, ] / totals^2))
        candidates <- which(cost == min(cost))
        chosen <- candidates[sample.int(length(candidates), 1L)]
        assigned[g] <- chosen
        balance[chosen, ] <- balance[chosen, ] + counts[g, ]
      }
    })
    assignment <- assigned[group_index]
  }
  counts <- table(factor(assignment, levels = seq_len(folds)), y)
  if (any(counts == 0)) sb_abort(
    "Generated folds do not all contain every class. Reduce folds or revise grouping.",
    "stormbreaker_infeasible_folds")
  as.integer(assignment)
}

# ---- task-survival ----

#' @rdname sb_binary_task
#' @export
sb_surv_task <- function(data, id_col, time_col, event_col, predictors,
                         time_unit, prediction_horizon = NULL,
                         decision_thresholds = NULL, group_col = NULL) {
  if (!is.data.frame(data)) {
    sb_abort("`data` must be a data.frame.", "stormbreaker_bad_data")
  }
  if (nrow(data) < 2L) {
    sb_abort(
      "`data` must contain at least two patients.",
      "stormbreaker_too_few_patients"
    )
  }
  if (anyDuplicated(names(data))) {
    sb_abort(
      "The data frame contains duplicate column names.",
      "stormbreaker_duplicate_column_names"
    )
  }

  sb_check_name(id_col, "id_col")
  sb_check_name(time_col, "time_col")
  sb_check_name(event_col, "event_col")
  endpoint_cols <- c(id_col, time_col, event_col)
  if (anyDuplicated(endpoint_cols)) {
    sb_abort(
      "ID, time, and event must use three different columns.",
      "stormbreaker_overlapping_endpoint_roles"
    )
  }

  missing_endpoint <- setdiff(endpoint_cols, names(data))
  if (length(missing_endpoint) > 0L) {
    sb_abort(
      paste0(
        "Missing endpoint columns: ",
        paste(missing_endpoint, collapse = ", "),
        "."
      ),
      "stormbreaker_missing_endpoint_columns"
    )
  }

  allowed_units <- c("day", "month", "year")
  if (!is.character(time_unit) || length(time_unit) != 1L ||
      is.na(time_unit) || !time_unit %in% allowed_units) {
    sb_abort(
      "`time_unit` must be exactly 'day', 'month', or 'year'.",
      "stormbreaker_bad_time_unit"
    )
  }

  if (!is.null(prediction_horizon)) {
    valid_horizon <- is.numeric(prediction_horizon) &&
      length(prediction_horizon) == 1L &&
      !is.na(prediction_horizon) &&
      is.finite(prediction_horizon) &&
      prediction_horizon > 0
    if (!valid_horizon) {
      sb_abort(
        "`prediction_horizon` must be one finite positive number.",
        "stormbreaker_bad_prediction_horizon"
      )
    }
    prediction_horizon <- as.numeric(prediction_horizon)
  }
  if (!is.null(decision_thresholds) && is.null(prediction_horizon)) {
    sb_abort(
      paste0(
        "`decision_thresholds` require a declared ",
        "`prediction_horizon`."
      ),
      "stormbreaker_decision_thresholds_require_prediction_horizon"
    )
  }
  if (!is.null(decision_thresholds)) {
    decision_thresholds <- sb_validate_decision_thresholds(
      decision_thresholds
    )
  }

  id <- data[[id_col]]
  if (anyNA(id) || any(!nzchar(trimws(as.character(id))))) {
    sb_abort(
      "Patient identifiers contain missing or blank values.",
      "stormbreaker_bad_id"
    )
  }
  if (anyDuplicated(id)) {
    sb_abort(
      "Patient identifiers must be unique.",
      "stormbreaker_duplicate_id"
    )
  }

  time <- data[[time_col]]
  if (!is.numeric(time) || anyNA(time) || any(!is.finite(time)) ||
      any(time < 0)) {
    sb_abort(
      "Survival time must contain finite, non-negative numeric values.",
      "stormbreaker_bad_time"
    )
  }

  event <- data[[event_col]]
  if ((!is.numeric(event) && !is.logical(event)) || anyNA(event) ||
      (is.numeric(event) && any(!is.finite(event))) ||
      !all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only finite 0/1 values.",
      "stormbreaker_bad_event"
    )
  }
  if (sum(event == 1) < 1L) {
    sb_abort(
      "At least one observed event is required.",
      "stormbreaker_no_events"
    )
  }
  if (sum(event == 0) < 1L) {
    sb_abort(
      "At least one censored patient is required.",
      "stormbreaker_no_censoring"
    )
  }

  if (!is.character(predictors) || length(predictors) < 1L ||
      anyNA(predictors) || any(!nzchar(predictors))) {
    sb_abort(
      "`predictors` must contain at least one non-empty column name.",
      "stormbreaker_bad_predictors"
    )
  }
  predictors <- unique(predictors)
  missing_predictors <- setdiff(predictors, names(data))
  if (length(missing_predictors) > 0L) {
    sb_abort(
      paste0(
        "Missing predictor columns: ",
        paste(missing_predictors, collapse = ", "),
        "."
      ),
      "stormbreaker_missing_predictors"
    )
  }
  if (length(intersect(predictors, endpoint_cols)) > 0L) {
    sb_abort(
      "Predictors cannot reuse the ID, time, or event columns.",
      "stormbreaker_reserved_predictor"
    )
  }
  numeric_predictors <- vapply(data[predictors], is.numeric, logical(1))
  if (!all(numeric_predictors)) {
    sb_abort(
      paste0(
        "All minimum-core predictors must be numeric: ",
        paste(names(numeric_predictors)[!numeric_predictors],
              collapse = ", "),
        "."
      ),
      "stormbreaker_bad_predictor_type"
    )
  }
  finite_predictors <- vapply(
    data[predictors],
    function(x) !anyNA(x) && all(is.finite(x)),
    logical(1)
  )
  if (!all(finite_predictors)) {
    sb_abort(
      paste0(
        "Predictors contain missing or non-finite values: ",
        paste(names(finite_predictors)[!finite_predictors],
              collapse = ", "),
        "."
      ),
      "stormbreaker_bad_predictor_values"
    )
  }

  group_col <- sb_validate_group_column(
    data, group_col, c(id_col, time_col, event_col), predictors
  )
  task <- list(
    data = data,
    id_col = id_col,
    time_col = time_col,
    event_col = event_col,
    predictors = predictors,
    time_unit = time_unit,
    n = nrow(data),
    n_events = as.integer(sum(event == 1)),
    n_censored = as.integer(sum(event == 0)),
    independent_unit = "patient"
  )
  if (!is.null(prediction_horizon)) {
    task <- append(
      task,
      list(prediction_horizon = prediction_horizon),
      after = match("time_unit", names(task))
    )
  }
  if (!is.null(decision_thresholds)) {
    task <- append(
      task,
      list(decision_thresholds = decision_thresholds),
      after = match("prediction_horizon", names(task))
    )
  }
  if (!is.null(group_col)) task$group_col <- group_col
  class(task) <- "stormbreaker_surv_task"
  task
}

sb_validate_task_contract <- function(task) {
  metadata_fields <- c(
    "id_col", "time_col", "event_col", "predictors", "time_unit", "n",
    "n_events", "n_censored", "independent_unit"
  )
  if ("group_col" %in% names(task)) {
    metadata_fields <- c(metadata_fields, "group_col")
  }
  if ("prediction_horizon" %in% names(task)) {
    metadata_fields <- append(
      metadata_fields,
      "prediction_horizon",
      after = match("time_unit", metadata_fields)
    )
  }
  if ("decision_thresholds" %in% names(task)) {
    if (!"prediction_horizon" %in% names(task)) {
      sb_abort(
        paste0(
          "`decision_thresholds` require a declared ",
          "`prediction_horizon`."
        ),
        "stormbreaker_decision_thresholds_require_prediction_horizon"
      )
    }
    metadata_fields <- append(
      metadata_fields,
      "decision_thresholds",
      after = match("prediction_horizon", metadata_fields)
    )
  }
  required_fields <- c("data", metadata_fields)
  if (!is.list(task) || length(setdiff(required_fields, names(task))) > 0L) {
    sb_abort(
      "The survival task is incomplete; recreate it with `sb_surv_task()`.",
      "stormbreaker_bad_task"
    )
  }
  prediction_horizon <- if ("prediction_horizon" %in% names(task)) {
    task[["prediction_horizon"]]
  } else {
    NULL
  }
  decision_thresholds <- if ("decision_thresholds" %in% names(task)) {
    task[["decision_thresholds"]]
  } else {
    NULL
  }

  validated <- sb_surv_task(
    data = task$data,
    id_col = task$id_col,
    time_col = task$time_col,
    event_col = task$event_col,
    predictors = task$predictors,
    time_unit = task$time_unit,
    prediction_horizon = prediction_horizon,
    decision_thresholds = decision_thresholds,
    group_col = task$group_col
  )
  metadata_ok <- identical(
    unclass(task)[metadata_fields],
    unclass(validated)[metadata_fields]
  )
  if (!metadata_ok) {
    sb_abort(
      "The survival task metadata changed; recreate it with `sb_surv_task()`.",
      "stormbreaker_mutated_task"
    )
  }

  validated
}

#' @export
print.stormbreaker_surv_task <- function(x, ...) {
  cat("<stormbreaker_surv_task>\n")
  cat(if (is.null(x$group_col)) "Patients:" else "Observations:",
      x$n, "| Events:", x$n_events,
      "| Censored:", x$n_censored, "\n")
  if (!is.null(x$group_col)) {
    cat("Groups:", length(unique(sb_task_group_ids(x))), "| Column:",
        x$group_col, "\n")
  }
  cat("Time:", x$time_col, "[", x$time_unit, "]",
      "| Event:", x$event_col, "\n")
  prediction_horizon <- if ("prediction_horizon" %in% names(x)) {
    x[["prediction_horizon"]]
  } else {
    NULL
  }
  if (!is.null(prediction_horizon)) {
    cat(
      "Prediction horizon:", format(prediction_horizon),
      "[", x$time_unit, "]\n"
    )
  }
  if ("decision_thresholds" %in% names(x)) {
    cat(
      "Decision thresholds:", length(x$decision_thresholds),
      "| Event risk by horizon\n"
    )
  }
  cat("Predictors:", length(x$predictors), "\n")
  invisible(x)
}

# ---- example-data ----

# ---- example-binary-data ----

#' @rdname sb_binary_task
#' @export
sb_example_binary <- function(n = 60L, seed = 20260818L) {
  n <- sb_check_positive_integer(n, "n", minimum = 10L)
  seed <- sb_check_positive_integer(seed, "seed", minimum = 1L)

  withr::with_seed(seed, {
    outcome <- sample(rep(c(0L, 1L), length.out = n))
    data.frame(
      patient_id = sprintf("B%03d", seq_len(n)),
      outcome = outcome,
      gene_a = stats::rnorm(n, mean = 0.5 * outcome),
      gene_b = stats::rnorm(n),
      gene_c = stats::rnorm(n),
      check.names = FALSE
    )
  })
}

# ---- example-data ----

#' @rdname sb_binary_task
#' @export
sb_example_survival <- function(n = 60L, seed = 20260812L) {
  n <- sb_check_positive_integer(n, "n", minimum = 10L)
  seed <- sb_check_positive_integer(seed, "seed", minimum = 1L)

  withr::with_seed(seed, {
    gene_a <- stats::rnorm(n)
    gene_b <- stats::rnorm(n)
    gene_c <- stats::rnorm(n)
    latent_rate <- exp(0.55 * gene_a - 0.15 * gene_b)
    event_time <- stats::rexp(n, rate = 0.055 * latent_rate) + 0.1
    target_events <- floor(n / 2)
    ordered_time <- sort(event_time)
    censor_time <- mean(ordered_time[c(target_events, target_events + 1L)])
    event <- as.integer(event_time <= censor_time)
    observed_time <- pmin(event_time, censor_time)

    data.frame(
      patient_id = sprintf("P%03d", seq_len(n)),
      time = round(observed_time, 3),
      event = event,
      gene_a = gene_a,
      gene_b = gene_b,
      gene_c = gene_c,
      check.names = FALSE
    )
  })
}
