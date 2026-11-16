# ---- resampling ----

sb_validate_group_ids <- function(group_ids, n) {
  valid <- is.atomic(group_ids) && is.null(dim(group_ids)) &&
    length(group_ids) == n && !anyNA(group_ids) &&
    all(nzchar(trimws(as.character(group_ids)))) &&
    (!is.numeric(group_ids) || all(is.finite(group_ids)))
  if (!valid) {
    sb_abort(
      "Group identifiers must be non-missing, non-blank values, one per row.",
      "stormbreaker_bad_group_ids"
    )
  }
  unname(as.character(group_ids))
}

sb_validate_group_column <- function(data, group_col, endpoint_cols,
                                      predictors) {
  if (is.null(group_col)) return(NULL)
  sb_check_name(group_col, "group_col")
  if (!group_col %in% names(data)) {
    sb_abort("The group column is missing.", "stormbreaker_missing_group_column")
  }
  if (group_col %in% c(endpoint_cols, predictors)) {
    sb_abort(
      "The group column cannot also be an ID, outcome, or predictor column.",
      "stormbreaker_overlapping_group_roles"
    )
  }
  sb_validate_group_ids(data[[group_col]], nrow(data))
  group_col
}

sb_task_group_ids <- function(task, rows = seq_len(task$n)) {
  if (is.null(task$group_col)) return(NULL)
  as.character(task$data[[task$group_col]][rows])
}

#' @rdname sb_binary_task
#' @export
sb_split_binary <- function(task, train_prop = 0.7, seed = 1L) {
  if (!inherits(task, "stormbreaker_binary_task")) {
    sb_abort("Supply a binary task.", "stormbreaker_bad_task")
  }
  task <- sb_validate_binary_task_contract(task)
  if (!is.numeric(train_prop) || length(train_prop) != 1L ||
      !is.finite(train_prop) || train_prop <= 0 || train_prop >= 1) {
    sb_abort("train_prop must lie strictly between 0 and 1.",
      "stormbreaker_bad_split_proportion")
  }
  seed <- sb_check_positive_integer(seed, "seed", minimum = 1L)
  ids <- as.character(task$data[[task$id_col]])
  group <- sb_task_group_ids(task)
  if (is.null(group)) group <- ids
  group_index <- match(group, unique(group))
  y <- as.integer(task$data[[task$outcome_col]])
  n_groups <- max(group_index)
  has_0 <- tabulate(group_index[y == 0L], n_groups) > 0L
  has_1 <- tabulate(group_index[y == 1L], n_groups) > 0L
  strata <- split(seq_len(n_groups), factor(has_0 + 2L * has_1, levels = 1:3))
  n <- unname(lengths(strata))
  at_least <- as.integer(n >= 2L)
  target <- pmax(at_least, pmin(n - at_least, floor(train_prop * n + 0.5)))

  # Unsplittable singleton strata can go to either side. Among the feasible
  # allocations choose the one closest to the target number of groups.
  choices <- lapply(seq_along(n), function(i) if (n[i] == 1L) 0:1 else target[i])
  candidates <- as.matrix(expand.grid(choices, KEEP.OUT.ATTRS = FALSE))
  has_both <- function(counts) {
    (counts[, 1L] + counts[, 3L] > 0) & (counts[, 2L] + counts[, 3L] > 0)
  }
  feasible <- has_both(candidates) &
    has_both(matrix(n, nrow(candidates), 3L, byrow = TRUE) - candidates)
  if (!any(feasible)) {
    sb_abort("Both partitions need both outcomes in disjoint patient groups.",
      "stormbreaker_insufficient_split_groups")
  }
  candidates <- candidates[feasible, , drop = FALSE]
  counts <- candidates[which.min(abs(rowSums(candidates) - train_prop * n_groups)), ]
  development_groups <- withr::with_seed(seed, unlist(Map(function(g, size) {
    if (size == 0) return(integer())
    g[sample.int(length(g), size)]
  }, strata, as.integer(counts)), use.names = FALSE))
  is_development <- group_index %in% development_groups
  role <- ifelse(is_development, "development", "assessment")
  development <- sb_binary_task(task$data[is_development, , drop = FALSE],
    task$id_col, task$outcome_col, task$predictors,
    decision_thresholds = task$decision_thresholds, group_col = task$group_col)
  summary <- do.call(rbind, lapply(c("development", "assessment"), function(part) {
    rows <- role == part
    data.frame(role = part, n_groups = length(unique(group[rows])),
      n_rows = sum(rows), n_negative = sum(y[rows] == 0L),
      n_positive = sum(y[rows] == 1L))
  }))
  assessment <- task$data[!is_development, , drop = FALSE]
  attr(assessment, "stormbreaker_assessment_role") <- "internal_holdout"
  list(development = development, assessment = assessment,
    assignment = data.frame(observation_id = ids, group_id = group,
      outcome = y, role = role),
    counts = summary, train_prop = train_prop, seed = seed)
}

sb_make_group_stratified_fold_id <- function(event, group_ids, folds, seed,
                                             error_class,
                                             require_both_classes = TRUE) {
  group_index <- match(group_ids, unique(group_ids))
  counts <- cbind(
    tabulate(group_index[event == 0L], max(group_index)),
    tabulate(group_index[event == 1L], max(group_index))
  )
  if (require_both_classes && any(colSums(counts > 0L) < folds)) {
    sb_abort(
      paste0("Cannot create ", folds, " grouped stratified folds: each ",
             "outcome must occur in at least ", folds, " distinct groups."),
      error_class
    )
  }
  if (!require_both_classes && nrow(counts) < folds) {
    sb_abort(paste0("Cannot create ", folds, " grouped folds from fewer than ",
                    folds, " distinct groups."), error_class)
  }
  group_fold <- integer(nrow(counts))
  fold_counts <- matrix(0, nrow = folds, ncol = 2L)
  totals <- colSums(counts)
  withr::with_seed(seed, {
    order_groups <- sample.int(nrow(counts))
    # Mixed-outcome groups first; then larger groups. This preserves both
    # outcomes in every fold whenever each occurs in at least K groups.
    order_groups <- order_groups[order(
      rowSums(counts[order_groups, , drop = FALSE] > 0L) != 2L,
      -rowSums(counts[order_groups, , drop = FALSE])
    )]
    for (group in order_groups) {
      # Minimize the increase in squared per-class fold proportions.
      cost <- drop(fold_counts %*% (counts[group, ] / totals^2))
      candidates <- which(cost == min(cost))
      if (!require_both_classes && any(rowSums(fold_counts) == 0L)) {
        # Probability calibration permits single-class training subproblems,
        # unlike per-fold AUC/deviance tuning. Still keep all K folds nonempty.
        candidates <- which(rowSums(fold_counts) == 0L)
        candidates <- candidates[cost[candidates] == min(cost[candidates])]
      }
      chosen <- candidates[[sample.int(length(candidates), 1L)]]
      group_fold[group] <- chosen
      fold_counts[chosen, ] <- fold_counts[chosen, ] + counts[group, ]
    }
  })
  group_fold[group_index]
}

sb_build_outer_folds <- function(task, assignment, source, seed) {
  id <- as.character(task$data[[task$id_col]])
  n_folds <- as.integer(length(unique(assignment)))
  splits <- lapply(seq_len(n_folds), function(fold_id) {
    assessment_rows <- which(assignment == fold_id)
    analysis_rows <- which(assignment != fold_id)
    list(
      fold = fold_id,
      analysis_rows = analysis_rows,
      assessment_rows = assessment_rows,
      analysis_ids = id[analysis_rows],
      assessment_ids = id[assessment_rows]
    )
  })

  structure(
    list(
      assignment = as.integer(assignment),
      splits = splits,
      source = source,
      seed = as.integer(seed),
      n_folds = n_folds
    ),
    class = "stormbreaker_folds"
  )
}

sb_make_stratified_folds <- function(task, folds, seed) {
  folds <- sb_check_positive_integer(folds, "folds", minimum = 2L)
  seed <- sb_check_positive_integer(seed, "seed", minimum = 1L)
  stratum_col <- if (inherits(task, "stormbreaker_binary_task")) {
    task$outcome_col
  } else {
    task$event_col
  }
  stratum_label <- if (inherits(task, "stormbreaker_binary_task")) {
    "outcome"
  } else {
    "event"
  }
  stratum_value <- as.integer(task$data[[stratum_col]])
  group_ids <- sb_task_group_ids(task)
  if (!is.null(group_ids) && anyDuplicated(group_ids)) {
    assignment <- sb_make_group_stratified_fold_id(
      stratum_value, group_ids, folds, seed, "stormbreaker_infeasible_folds"
    )
    return(sb_build_outer_folds(
      task, assignment, "generated_group_stratified", seed
    ))
  }
  strata_counts <- table(factor(stratum_value, levels = c(0, 1)))

  if (task$n < folds || any(strata_counts < folds)) {
    sb_abort(
      paste0(
        "Cannot create ", folds,
        " ", stratum_label,
        "-stratified folds: each of ", stratum_label, "=0 and ",
        stratum_label, "=1 must have ",
        "at least ", folds, " patients."
      ),
      "stormbreaker_infeasible_folds"
    )
  }

  assignment <- integer(task$n)
  withr::with_seed(seed, {
    for (stratum in c(0, 1)) {
      rows <- which(stratum_value == stratum)
      rows <- sample(rows, length(rows), replace = FALSE)
      assignment[rows] <- rep(seq_len(folds), length.out = length(rows))
    }
  })

  sb_build_outer_folds(
    task = task,
    assignment = assignment,
    source = "generated_stratified",
    seed = seed
  )
}

sb_make_preassigned_folds <- function(task, assignment) {
  id <- as.character(task$data[[task$id_col]])
  if (!is.numeric(assignment)) {
    sb_abort(
      "A preassigned `folds` vector must be numeric.",
      "stormbreaker_bad_fold_assignment"
    )
  }
  if (length(assignment) != task$n) {
    sb_abort(
      paste0(
        "A preassigned `folds` vector must contain exactly one value for ",
        "every task patient ID."
      ),
      "stormbreaker_bad_fold_assignment"
    )
  }
  assignment_names <- names(assignment)
  if (is.null(assignment_names)) {
    sb_abort(
      "A preassigned `folds` vector must be named with task patient IDs.",
      "stormbreaker_bad_fold_assignment"
    )
  }
  if (anyNA(assignment_names) || any(!nzchar(assignment_names))) {
    sb_abort(
      "Every preassigned `folds` name must be a non-empty patient ID.",
      "stormbreaker_bad_fold_assignment"
    )
  }
  if (anyDuplicated(assignment_names) > 0L) {
    sb_abort(
      "Preassigned `folds` names must be unique patient IDs.",
      "stormbreaker_bad_fold_assignment"
    )
  }
  if (!setequal(assignment_names, id)) {
    sb_abort(
      paste0(
        "Preassigned `folds` names must match the complete set of task ",
        "patient IDs."
      ),
      "stormbreaker_bad_fold_assignment"
    )
  }
  valid_values <- !anyNA(assignment) && all(is.finite(assignment)) &&
    all(assignment > 0) && all(assignment <= .Machine$integer.max) &&
    all(assignment == floor(assignment))
  if (!valid_values) {
    sb_abort(
      paste0(
        "Preassigned `folds` values must be finite positive whole numbers ",
        "within the R integer range."
      ),
      "stormbreaker_bad_fold_assignment"
    )
  }

  aligned <- as.integer(assignment[match(id, assignment_names)])
  labels <- sort(unique(aligned))
  if (length(labels) < 2L) {
    sb_abort(
      "A preassigned `folds` vector must contain at least two folds.",
      "stormbreaker_bad_fold_assignment"
    )
  }
  if (!identical(labels, seq_len(length(labels)))) {
    sb_abort(
      "Preassigned `folds` labels must be the consecutive integers 1 through K.",
      "stormbreaker_bad_fold_assignment"
    )
  }

  group_ids <- sb_task_group_ids(task)
  if (!is.null(group_ids) && any(vapply(split(aligned, group_ids), function(x) {
    length(unique(x)) != 1L
  }, logical(1)))) {
    sb_abort(
      "Preassigned folds split observations from the same group across folds.",
      "stormbreaker_group_fold_overlap"
    )
  }

  sb_build_outer_folds(
    task = task,
    assignment = aligned,
    source = "preassigned_patient_id",
    seed = NA_integer_
  )
}

sb_make_outer_folds <- function(task, folds, seed) {
  if (length(folds) <= 1L) {
    return(sb_make_stratified_folds(task, folds, seed))
  }
  sb_make_preassigned_folds(task, folds)
}
