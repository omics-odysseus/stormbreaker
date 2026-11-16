#' Cohort meta-analysis
#' @name sb_meta_external_discrimination
#' @description Synthesize compatible cohort-level estimates and paired model differences.
#' @rawRd \arguments{
#' \item{bootstraps}{
#' A named list of at least two
#' \code{stormbreaker_external_bootstrap} objects. Names must be unique, non-blank
#' cohort IDs. All objects must describe the same task, pipeline ID, and
#' learner label.
#' }
#' \item{confidence_level}{
#' \strong{For \code{sb_meta_external_discrimination()}, \code{sb_meta_external_calibration_slope()}, \code{sb_meta_external_mean_calibration()}:} Confidence level strictly between zero and one for
#' the pooled confidence interval and between-cohort prediction interval.
#' Defaults to 0.95.
#' 
#' \strong{For \code{sb_meta_external_discrimination_difference()}, \code{sb_meta_external_prediction_error_difference()}:} Confidence level strictly between 0 and 1. Defaults
#' to 0.95.
#' }
#' \item{metric}{
#' \strong{For \code{sb_meta_external_discrimination()}:} Discrimination estimand. \code{NULL} selects \code{"roc_auc"} for a
#' binary task and \code{"harrell_c"} for a survival task. Binary tasks otherwise
#' accept only \code{"roc_auc"}; survival tasks accept \code{"harrell_c"},
#' \code{"time_auc"}, or \code{"uno_c"}. Time-specific choices require one common
#' retained time unit and prediction horizon.
#' 
#' \strong{For \code{sb_meta_external_discrimination_difference()}:} \code{NULL} for the task-specific default, or one explicit
#' discrimination metric. Binary tasks allow only \code{"roc_auc"}; survival tasks
#' allow \code{"harrell_c"}, \code{"time_auc"}, or \code{"uno_c"}.
#' 
#' \strong{For \code{sb_meta_external_prediction_error_difference()}:} \code{NULL} for the task-specific default, or one explicit
#' prediction-error metric. Binary tasks allow \code{"brier_score"} or
#' \code{"log_loss"}; survival tasks allow only \code{"ipcw_brier_score"}.
#' }
#' \item{comparisons}{
#' A named list of at least two
#' \code{stormbreaker_external_bootstrap_comparison} objects returned by
#' \code{\link[=sb_bootstrap_external_comparison]{sb_bootstrap_external_comparison()}}. Outer names are treated as cohort
#' IDs and must be unique and non-blank.
#' }
#' \item{candidate_pipeline_id}{
#' Exact pipeline ID of one candidate retained in
#' every source comparison. The reference and candidate must have been
#' declared independently of external outcomes.
#' }
#' \item{cohorts}{
#' Named list of surv_risk_plots results for independent cohorts.
#' Use the same saved model and training cutoff in every cohort. Do not include
#' a pooled cohort alongside its constituent cohorts; IDs must not overlap.
#' }
#' \item{type}{
#' group compares High with Low. continuous reports HR per one unit
#' of the saved model score, without per-cohort standardization.
#' }
#' \item{x}{
#' A saved external discrimination, mean-calibration,
#' calibration-slope, discrimination-difference or prediction-error-difference
#' Meta-analysis object.
#' }
#' \item{show_prediction_interval}{
#' Logical; show the saved between-cohort
#' prediction interval on a separate row, distinct from the pooled CI.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_meta_external_discrimination()}}{
#' \code{sb_meta_external_discrimination()} performs one deliberately narrow
#' performance meta-analysis for the same locked prediction model evaluated in
#' multiple independent external cohorts. The default is ROC AUC for binary
#' tasks and horizon-free Harrell concordance for survival tasks. A survival
#' caller may instead explicitly request fixed-horizon cumulative/dynamic
#' \code{time_auc} or truncated \code{uno_c}; the function never chooses among these
#' estimands automatically. It does not pool calibration, prediction error, or
#' clinical utility.
#' 
#' \strong{Returns:} A \code{stormbreaker_external_discrimination_meta} containing cohort
#' point estimates and bootstrap standard errors, the random-effects summary,
#' explicit \code{metric_target}, method labels, software provenance,
#' prediction-interval caution, and the unchanged named source \code{bootstraps}.
#' No p value, ranking, or winner is returned.
#' 
#' \strong{Details:} Fixed-horizon \code{time_auc} and \code{uno_c} require every source bootstrap to retain
#' a complete \code{prediction_target} and all cohorts to share the exact same time
#' unit and prediction horizon. A horizon may have been supplied by the endpoint
#' or inherited from the locked prediction; only the effective target must
#' match. Harrell C remains horizon-free even when its source bootstrap also
#' contains horizon-specific metrics. The three survival discrimination
#' definitions are never pooled together.
#' 
#' Each input is a completed \code{\link[=sb_bootstrap_external]{sb_bootstrap_external()}} result. The cohort point
#' estimate is transformed with \code{qlogis()}, and its within-cohort standard error
#' is the sample standard deviation of all ordinary patient-bootstrap estimates
#' on that transformed scale. A random-effects model is then fitted with REML
#' and Knapp-Hartung inference through \code{metafor::rma.uni()}. The pooled estimate,
#' confidence interval, and prediction interval are back-transformed with
#' \code{plogis()}. Between-cohort variance and its square root remain on the logit
#' discrimination scale.
#' 
#' Every requested target bootstrap replicate must occur exactly once, be
#' successful and finite, and lie strictly between zero and one. The function
#' neither drops unavailable replicates nor applies boundary clipping or
#' continuity corrections. Failures of unrelated bootstrap metrics do not
#' invalidate an otherwise complete discrimination target.
#' 
#' At least two cohorts are required because one cohort cannot identify
#' between-cohort heterogeneity. Estimates of heterogeneity and prediction
#' intervals remain especially uncertain with few cohorts, and frequentist
#' coverage can be poor in some small or imbalanced meta-analyses. The returned
#' caution is therefore part of the result rather than a model-selection rule.
#' 
#' Matching \code{pipeline_id} and learner labels are necessary bookkeeping only.
#' They do not prove that cohorts are independent, use common eligibility and
#' outcome definitions, or received the identical fitted state, coefficient
#' formula, preprocessing, and predictor orientation. Those design conditions
#' must be established by the caller. The bootstrap uncertainty is conditional
#' on each locked model and prediction set and does not recover development-
#' fitting, feature-selection, tuning, or model-selection uncertainty.
#' }
#' 
#' \subsection{\code{sb_meta_external_discrimination_difference()}}{
#' \code{sb_meta_external_discrimination_difference()} synthesizes one direct,
#' paired comparison between a predeclared candidate and reference prediction
#' across independent external cohorts. Each source must come from
#' \code{\link[=sb_bootstrap_external_comparison]{sb_bootstrap_external_comparison()}}, which applies the same patient draw to
#' both models. The cohort effect is the raw discrimination difference
#' \code{candidate - reference}; its within-cohort standard error is the sample
#' standard deviation of all corresponding paired bootstrap differences.
#' 
#' \strong{Returns:} A \code{stormbreaker_external_discrimination_difference_meta} containing
#' the reference and candidate IDs, explicit metric target, cohort point
#' differences and paired-bootstrap standard errors, random-effects summary,
#' method labels, software provenance, cautions, and the unchanged named
#' source \code{comparisons}. No p value, ranking, winner, or model is returned.
#' 
#' \strong{Details:} Binary tasks use ROC AUC. Survival tasks default to horizon-free Harrell
#' concordance and may explicitly request cumulative/dynamic \code{time_auc} or
#' truncated \code{uno_c}. Time-specific metrics require every source comparison to
#' retain the exact same prediction horizon and time unit. The function never
#' mixes discrimination definitions or chooses a metric automatically beyond
#' the documented task-specific default.
#' 
#' Differences are meta-analysed on their raw scale with REML and Knapp-Hartung
#' inference through \code{metafor::rma.uni()}. Point estimates and every requested
#' target-metric replicate must be present exactly once, successful, and
#' finite. Each saved difference must exactly equal its candidate estimate
#' minus its reference estimate, and the synthesis uses that recomputed value.
#' Failed replicates are not discarded, differences are not clipped, and
#' confidence or prediction limits outside the mathematical interval \verb{[-1, 1]}
#' remain visible as a warning about the normal random-effects model.
#' 
#' A positive difference indicates higher candidate discrimination only. It is
#' not evidence of better calibration, prediction error, clinical utility, or
#' overall usefulness. The function does not report a p value, rank models,
#' select a winner, make indirect or network comparisons, or use external
#' outcomes to choose the reference or candidate. Matching pipeline labels do
#' not prove identical fitted states across cohorts, and cohort independence
#' remains the study team's responsibility.
#' }
#' 
#' \subsection{\code{sb_meta_external_calibration_slope()}}{
#' \code{sb_meta_external_calibration_slope()} pools one calibration-slope estimate
#' for the same locked prediction model evaluated in multiple independent
#' external cohorts. It accepts binary results and fixed-horizon right-censored
#' survival results from \code{\link[=sb_bootstrap_external]{sb_bootstrap_external()}}. Binary and survival tasks
#' are analysed separately and are never mixed.
#' 
#' \strong{Returns:} A \code{stormbreaker_external_calibration_slope_meta} containing cohort
#' slope estimates and bootstrap standard errors, the random-effects summary,
#' explicit time target, method labels, software provenance, cautions, and the
#' unchanged named source \code{bootstraps}. No p value, ranking, winner, or
#' recalibrated model is returned.
#' 
#' \strong{Details:} For binary outcomes, the source slope is from logistic recalibration of the
#' observed outcome on the logit predicted probability with a freely estimated
#' intercept. For survival outcomes, the source slope is the coefficient from
#' an Efron Cox recalibration model after administrative censoring at the locked
#' horizon, using \code{log(-log(S(h)))} as the transformed prediction. Every
#' survival cohort must retain the exact same time unit and horizon.
#' 
#' Calibration slopes are meta-analysed on their original, unbounded scale.
#' The within-cohort standard error is the sample standard deviation of every
#' ordinary patient-bootstrap slope. A random-effects model is fitted with REML
#' and Knapp-Hartung inference through \code{metafor::rma.uni()}. Finite negative
#' slopes and negative interval bounds are retained rather than clipped. A
#' slope of one is the ideal scaling target under the stated recalibration
#' model; a slope below one commonly indicates predictions that are too
#' extreme, whereas a slope above one commonly indicates insufficient spread.
#' A negative slope indicates inverse ordering under that recalibration model,
#' not merely excessive spread.
#' 
#' The point estimate and every requested target bootstrap replicate must occur
#' exactly once and be successful and finite. The function does not drop failed
#' replicates, transform slopes, or substitute another calibration measure.
#' Failures of unrelated bootstrap metrics do not invalidate an otherwise
#' complete calibration-slope target.
#' 
#' At least two cohorts are required. Heterogeneity and prediction intervals
#' remain especially uncertain with few cohorts. A pooled slope near one does
#' not establish mean or moderate calibration and can coexist with a poor O:E
#' ratio or calibration curve. This is a univariate synthesis and does not use
#' within-cohort correlations with discrimination or other calibration
#' measures. Matching pipeline and learner labels are bookkeeping only: the
#' caller must establish cohort independence, common eligibility and outcome
#' definitions, and an identical fitted model, preprocessing, predictor
#' orientation, and prediction horizon. Bootstrap uncertainty is conditional
#' on each locked model and prediction set and does not recover development-
#' fitting or model-selection uncertainty.
#' }
#' 
#' \subsection{\code{sb_meta_external_mean_calibration()}}{
#' \code{sb_meta_external_mean_calibration()} pools one observed-to-expected event
#' ratio (O:E) for the same locked prediction model evaluated in multiple
#' independent external cohorts. It accepts binary results and fixed-horizon
#' right-censored survival results from \code{\link[=sb_bootstrap_external]{sb_bootstrap_external()}}. It does not
#' pool calibration slopes, calibration curves, prediction error, or clinical
#' utility.
#' 
#' \strong{Returns:} A \code{stormbreaker_external_mean_calibration_meta} containing cohort
#' O:E estimates and bootstrap standard errors, the random-effects summary,
#' explicit time target, method labels, software provenance, cautions, and the
#' unchanged named source \code{bootstraps}. No p value, ranking, winner, or
#' recalibrated model is returned.
#' 
#' \strong{Details:} For binary outcomes, O:E is the observed event proportion divided by the
#' mean predicted event probability. For survival outcomes, it is the
#' Kaplan-Meier observed event probability divided by the mean predicted event
#' probability at the locked horizon. Every survival cohort must retain the
#' exact same time unit and horizon. O:E equal to one indicates agreement in
#' mean risk; values below one indicate average overprediction, and values above
#' one indicate average underprediction. O:E alone cannot reveal miscalibration
#' within the prediction range.
#' 
#' Each cohort estimate is transformed with \code{log()}, and its within-cohort
#' standard error is the sample standard deviation of all ordinary
#' patient-bootstrap estimates on that log scale. A random-effects model is
#' fitted with REML and Knapp-Hartung inference through \code{metafor::rma.uni()}.
#' The pooled estimate, confidence interval, and prediction interval are
#' exponentiated back to the ratio scale. Between-cohort variance and its square
#' root remain on the log O:E scale.
#' 
#' The point estimate and every requested target bootstrap replicate must occur
#' exactly once, be successful, finite, and strictly positive. The function
#' neither drops unavailable replicates nor applies a continuity correction.
#' This means that a cohort with zero observed events in the point estimate or
#' any bootstrap draw is not eligible for this log-ratio meta-analysis.
#' Failures of unrelated bootstrap metrics do not invalidate an otherwise
#' complete O:E target.
#' 
#' At least two cohorts are required. Heterogeneity and prediction intervals
#' remain especially uncertain with few cohorts, and a pooled O:E near one can
#' hide opposite directions of miscalibration across cohorts. Matching pipeline
#' and learner labels are bookkeeping only: the caller must establish cohort
#' independence, common eligibility and outcome definitions, and an identical
#' fitted model, preprocessing, predictor orientation, and prediction horizon.
#' Bootstrap uncertainty is conditional on each locked model and prediction set
#' and does not recover development-fitting or model-selection uncertainty.
#' }
#' 
#' \subsection{\code{sb_meta_external_prediction_error_difference()}}{
#' \code{sb_meta_external_prediction_error_difference()} synthesizes one direct,
#' paired comparison between a predeclared candidate and reference prediction
#' across independent external cohorts. Each source must come from
#' \code{\link[=sb_bootstrap_external_comparison]{sb_bootstrap_external_comparison()}}, which applies the same patient draw to
#' both models. The cohort effect is the raw prediction-error difference
#' \code{candidate - reference}; its within-cohort standard error is the sample
#' standard deviation of all corresponding paired bootstrap differences.
#' 
#' \strong{Returns:} A \code{stormbreaker_external_prediction_error_difference_meta}
#' containing the reference and candidate IDs, explicit metric target,
#' cohort point differences and paired-bootstrap standard errors,
#' random-effects summary, method labels, software provenance, cautions, and
#' the unchanged named source \code{comparisons}. No p value, ranking, winner, or
#' model is returned.
#' 
#' \strong{Details:} Binary tasks default to Brier score and may explicitly request log loss.
#' Survival tasks use the fixed-horizon IPCW Brier score and therefore require
#' every cohort to retain the exact same prediction horizon and time unit.
#' Lower values are better for every supported metric, so a negative difference
#' means lower prediction error for the candidate.
#' 
#' Differences are meta-analysed on their raw scale with REML and Knapp-Hartung
#' inference through \code{metafor::rma.uni()}. Point estimates and every requested
#' target-metric replicate must be present exactly once, successful, and
#' finite. Each saved difference must exactly equal its candidate estimate
#' minus its reference estimate, and the synthesis uses that recomputed value.
#' Failed replicates are not discarded and confidence or prediction limits are
#' not clipped.
#' 
#' A lower prediction error is not by itself evidence of clinical utility, nor
#' does it identify whether a difference arose from calibration,
#' discrimination, or both. The function does not report a p value, rank
#' models, select a winner, make indirect or network comparisons, or use
#' external outcomes to choose the reference, candidate, or metric. Prediction
#' error depends on cohort outcome incidence and case mix; the survival metric
#' additionally relies on the package's stated marginal censoring assumptions.
#' Matching pipeline labels do not prove identical fitted states across
#' cohorts, and cohort independence remains the study team's responsibility.
#' }
#' 
#' \subsection{\code{surv_meta_risk()}}{
#' Pool independent cohorts' risk-score Cox hazard ratios
#' 
#' \strong{Returns:} Cohort estimates, fixed and random effects, heterogeneity, and a
#' forest plot. Failed cohort estimates are returned and not silently pooled.
#' 
#' \strong{Details:} Univariate Cox models use Efron ties and Wald intervals. Meta-analysis
#' uses metafor inverse-variance fixed effects and DerSimonian-Laird random
#' effects with normal intervals. At least two evaluable independent cohorts
#' are required. Matching IDs can detect overlap, not prove independence.
#' Training-cohort associations are apparent, not independent validation.
#' }
#' 
#' \subsection{\code{sb_plot_meta_forest()}}{
#' Plot saved external Meta-analysis results as a forest plot
#' 
#' \strong{Returns:} A ggplot object with cohort points, saved cohort bootstrap CIs when
#' present, the saved pooled CI and optionally the saved prediction interval.
#' 
#' \strong{Details:} No Meta model is refitted and no cohort CI is inferred from its
#' standard error. Cohort interval levels are those in their source bootstrap
#' objects and can differ from the pooled level; each row is labelled.
#' Missing source intervals are explicit while cohort point estimates remain
#' visible. Pooled and prediction endpoints are copied unchanged, including
#' values beyond usual metric bounds. A prediction interval describes
#' between-cohort variation, not an individual patient's prediction.
#' The original small-cohort caution remains in the caption. Differences
#' retain candidate-minus-reference direction; calibration target lines are
#' one, and difference target lines are zero.
#' }
#' }
#' @rawRd \references{
#' Snell KIE, Ensor J, Debray TPA, Moons KGM, Riley RD (2018). Meta-analysis of
#' prediction model performance across multiple studies: Which scale helps
#' ensure between-study normality for the C-statistic and calibration
#' measures? \emph{Statistical Methods in Medical Research}, 27(11), 3505-3522.
#' \doi{10.1177/0962280217705678}.
#' 
#' Debray TPA, Damen JAAG, Snell KIE, et al. (2017). A guide to systematic
#' review and meta-analysis of prediction model performance. \emph{BMJ}, 356,
#' i6460. \doi{10.1136/bmj.i6460}.
#' 
#' Debray TPA, Damen JAAG, Snell KIE, et al. (2019). A framework for
#' meta-analysis of prediction model studies with binary and time-to-event
#' outcomes. \emph{Statistical Methods in Medical Research}, 28(9), 2768-2786.
#' \doi{10.1177/0962280218785504}.
#' 
#' Partlett C, Riley RD (2017). Random effects meta-analysis: Coverage
#' performance of 95\% confidence and prediction intervals following REML
#' estimation. \emph{Statistics in Medicine}, 36(2), 301-317.
#' \doi{10.1002/sim.7140}.
#' 
#' Haile SR, Guerra B, Soriano JB, Puhan MA, 3CIA collaboration (2017).
#' Multiple Score Comparison: a network meta-analysis approach to comparison
#' and external validation of prognostic scores. \emph{BMC Medical Research
#' Methodology}, 17, 172. \doi{10.1186/s12874-017-0433-2}
#' 
#' Debray TPA, Damen JAAG, Snell KIE, et al. (2017). A guide to systematic
#' review and meta-analysis of prediction model performance. \emph{BMJ}, 356,
#' i6460. \doi{10.1136/bmj.i6460}
#' 
#' Snell KIE, Ensor J, Debray TPA, Moons KGM, Riley RD (2018). Meta-analysis of
#' prediction model performance across multiple studies: Which scale helps
#' ensure between-study normality for the C-statistic and calibration
#' measures? \emph{Statistical Methods in Medical Research}, 27(11), 3505-3522.
#' \doi{10.1177/0962280217705678}.
#' 
#' Debray TPA, Damen JAAG, Snell KIE, et al. (2019). A framework for
#' meta-analysis of prediction model studies with binary and time-to-event
#' outcomes. \emph{Statistical Methods in Medical Research}, 28(9), 2768-2786.
#' \doi{10.1177/0962280218785504}.
#' 
#' Snell KIE, Hua H, Debray TPA, et al. (2016). Multivariate meta-analysis of
#' individual participant data helped externally validate the performance and
#' implementation of a prediction model. \emph{Journal of Clinical Epidemiology},
#' 69, 40-50. \doi{10.1016/j.jclinepi.2015.05.009}.
#' 
#' Stevens RJ, Poppe KK (2020). Validation of clinical prediction models: what
#' does the "calibration slope" really measure? \emph{Journal of Clinical
#' Epidemiology}, 118, 93-99. \doi{10.1016/j.jclinepi.2019.09.016}.
#' 
#' Debray TPA, Damen JAAG, Snell KIE, et al. (2019). A framework for
#' meta-analysis of prediction model studies with binary and time-to-event
#' outcomes. \emph{Statistical Methods in Medical Research}, 28(9), 2768-2786.
#' \doi{10.1177/0962280218785504}.
#' 
#' Snell KIE, Ensor J, Debray TPA, Moons KGM, Riley RD (2018). Meta-analysis of
#' prediction model performance across multiple studies: Which scale helps
#' ensure between-study normality for the C-statistic and calibration
#' measures? \emph{Statistical Methods in Medical Research}, 27(11), 3505-3522.
#' \doi{10.1177/0962280217705678}.
#' 
#' Gneiting T, Raftery AE (2007). Strictly proper scoring rules, prediction,
#' and estimation. \emph{Journal of the American Statistical Association},
#' 102(477), 359-378. \doi{10.1198/016214506000001437}
#' 
#' Debray TPA, Damen JAAG, Snell KIE, et al. (2017). A guide to systematic
#' review and meta-analysis of prediction model performance. \emph{BMJ}, 356,
#' i6460. \doi{10.1136/bmj.i6460}
#' 
#' McLernon DJ, Giardiello D, Van Calster B, et al. (2023). Assessing
#' performance and clinical usefulness in prediction models with survival
#' outcomes: practical guidance for Cox proportional hazards models.
#' \emph{Annals of Internal Medicine}, 176(1), 105-114.
#' \doi{10.7326/M22-0844}
#' }
#' @rawRd \examples{
#' # sb_meta_external_discrimination
#' \dontrun{
#' pooled <- sb_meta_external_discrimination(list(
#'   hospital_a = bootstrap_a,
#'   hospital_b = bootstrap_b
#' ))
#' pooled$summary
#' sb_meta_external_discrimination(
#'   list(hospital_a = bootstrap_a, hospital_b = bootstrap_b),
#'   metric = "time_auc"
#' )
#' }
#' 
#' # sb_meta_external_discrimination_difference
#' \dontrun{
#' pooled_difference <- sb_meta_external_discrimination_difference(
#'   list(hospital_a = paired_a, hospital_b = paired_b),
#'   candidate_pipeline_id = "candidate_model"
#' )
#' pooled_difference$summary
#' }
#' 
#' # sb_meta_external_calibration_slope
#' \dontrun{
#' pooled <- sb_meta_external_calibration_slope(list(
#'   hospital_a = bootstrap_a,
#'   hospital_b = bootstrap_b
#' ))
#' pooled$summary
#' }
#' 
#' # sb_meta_external_mean_calibration
#' \dontrun{
#' pooled <- sb_meta_external_mean_calibration(list(
#'   hospital_a = bootstrap_a,
#'   hospital_b = bootstrap_b
#' ))
#' pooled$summary
#' }
#' 
#' # sb_meta_external_prediction_error_difference
#' \dontrun{
#' pooled_difference <- sb_meta_external_prediction_error_difference(
#'   list(hospital_a = paired_a, hospital_b = paired_b),
#'   candidate_pipeline_id = "candidate_model"
#' )
#' pooled_difference$summary
#' }
#' }
NULL

# ---- meta-calibration ----

# ---- meta-external-calibration-slope ----

sb_external_calibration_slope_meta_member <- function(
    bootstrap,
    cohort_id) {
  if (!inherits(bootstrap, "stormbreaker_external_bootstrap") ||
      !is.list(bootstrap)) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` is not a result from `sb_bootstrap_external()`."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  task_type <- bootstrap[["task_type"]]
  pipeline_id <- bootstrap[["pipeline_id"]]
  learner <- bootstrap[["learner"]]
  n <- bootstrap[["n"]]
  reps <- bootstrap[["reps"]]
  valid_metadata <- sb_external_scalar_character(task_type) &&
    task_type %in% c("binary", "survival") &&
    sb_external_scalar_character(pipeline_id) &&
    sb_external_scalar_character(learner) &&
    is.numeric(n) && length(n) == 1L && !is.na(n) && is.finite(n) &&
    n >= 1 && n <= .Machine$integer.max && n == as.integer(n) &&
    is.numeric(reps) && length(reps) == 1L && !is.na(reps) &&
    is.finite(reps) && reps >= 2 && reps <= .Machine$integer.max &&
    reps == as.integer(reps)
  if (!valid_metadata) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has incomplete external-bootstrap metadata."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  metric_target <- sb_external_calibration_meta_target(
    bootstrap,
    cohort_id,
    task_type,
    "survival calibration slope"
  )
  point <- bootstrap[["point_estimates"]]
  replicates <- bootstrap[["replicates"]]
  valid_tables <- is.data.frame(point) &&
    all(c("metric", "estimate", "status") %in% names(point)) &&
    is.character(point$metric) && is.numeric(point$estimate) &&
    is.character(point$status) &&
    is.data.frame(replicates) &&
    all(c("replicate", "metric", "estimate", "status") %in%
      names(replicates)) &&
    is.numeric(replicates$replicate) &&
    is.character(replicates$metric) && is.numeric(replicates$estimate) &&
    is.character(replicates$status)
  if (!valid_tables) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` does not retain usable bootstrap metric tables."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  metric <- "calibration_slope"
  point_rows <- point[
    !is.na(point$metric) & point$metric == metric,
    ,
    drop = FALSE
  ]
  replicate_rows <- replicates[
    !is.na(replicates$metric) & replicates$metric == metric,
    ,
    drop = FALSE
  ]
  expected_replicates <- seq_len(as.integer(reps))
  replicate_ids_complete <- nrow(replicate_rows) == as.integer(reps) &&
    !anyNA(replicate_rows$replicate) &&
    all(is.finite(replicate_rows$replicate)) &&
    all(replicate_rows$replicate >= 1 &
      replicate_rows$replicate <= as.integer(reps)) &&
    all(replicate_rows$replicate == floor(replicate_rows$replicate)) &&
    identical(
      sort(as.integer(replicate_rows$replicate)),
      expected_replicates
    )
  point_available <- nrow(point_rows) == 1L &&
    identical(point_rows$status[[1L]], "success") &&
    length(point_rows$estimate) == 1L &&
    !is.na(point_rows$estimate[[1L]]) &&
    is.finite(point_rows$estimate[[1L]])
  replicates_complete <- replicate_ids_complete &&
    !anyNA(replicate_rows$status) &&
    all(replicate_rows$status == "success") &&
    !anyNA(replicate_rows$estimate) &&
    all(is.finite(replicate_rows$estimate))
  if (!point_available || !replicates_complete) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` requires one finite successful `", metric,
        "` point estimate and all ", as.integer(reps),
        " successful finite bootstrap estimates, numbered exactly once from ",
        "1 through `reps`. No failed replicate is dropped and no slope is ",
        "transformed or clipped."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  standard_error <- stats::sd(replicate_rows$estimate)
  if (!is.finite(standard_error) || standard_error <= 0) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has zero or non-finite bootstrap uncertainty for `", metric,
        "` on the original scale."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  list(
    task_type = task_type,
    pipeline_id = pipeline_id,
    learner = learner,
    n = as.integer(n),
    reps = as.integer(reps),
    point_estimate = as.numeric(point_rows$estimate[[1L]]),
    standard_error = as.numeric(standard_error),
    prediction_horizon = metric_target$prediction_horizon,
    time_unit = metric_target$time_unit
  )
}

#' @rdname sb_meta_external_discrimination
#' @export
sb_meta_external_calibration_slope <- function(
    bootstraps,
    confidence_level = 0.95) {
  valid_confidence <- is.numeric(confidence_level) &&
    length(confidence_level) == 1L && !is.na(confidence_level) &&
    is.finite(confidence_level) && confidence_level > 0 &&
    confidence_level < 1
  if (!valid_confidence) {
    sb_abort(
      "`confidence_level` must be one finite number strictly between 0 and 1.",
      "stormbreaker_bad_confidence_level"
    )
  }

  valid_input <- is.list(bootstraps) &&
    !inherits(bootstraps, "stormbreaker_external_bootstrap")
  if (!valid_input) {
    sb_abort(
      paste0(
        "`bootstraps` must be a named list of results from ",
        "`sb_bootstrap_external()`."
      ),
      "stormbreaker_bad_external_meta"
    )
  }
  if (length(bootstraps) < 2L) {
    sb_abort(
      "At least two independent external cohorts are required.",
      "stormbreaker_insufficient_external_cohorts"
    )
  }

  cohort_ids <- names(bootstraps)
  valid_ids <- !is.null(cohort_ids) &&
    length(cohort_ids) == length(bootstraps) &&
    !anyNA(cohort_ids) && all(nzchar(trimws(cohort_ids))) &&
    anyDuplicated(cohort_ids) == 0L
  if (!valid_ids) {
    sb_abort(
      "Outer list names must be unique, non-blank cohort IDs.",
      "stormbreaker_bad_external_meta"
    )
  }

  members <- lapply(seq_along(bootstraps), function(index) {
    sb_external_calibration_slope_meta_member(
      bootstraps[[index]],
      cohort_ids[[index]]
    )
  })
  task_types <- vapply(members, `[[`, character(1L), "task_type")
  if (length(unique(task_types)) != 1L) {
    sb_abort(
      "All cohorts must use the same task type.",
      "stormbreaker_external_meta_task_mismatch"
    )
  }
  pipeline_ids <- vapply(members, `[[`, character(1L), "pipeline_id")
  learners <- vapply(members, `[[`, character(1L), "learner")
  if (length(unique(pipeline_ids)) != 1L ||
      length(unique(learners)) != 1L) {
    sb_abort(
      paste0(
        "All cohorts must retain the same pipeline ID and learner label. ",
        "These labels do not prove identical fitted state or formula."
      ),
      "stormbreaker_external_meta_model_mismatch"
    )
  }

  task_type <- task_types[[1L]]
  metric_target <- list(
    has_horizon = identical(task_type, "survival"),
    prediction_horizon = NULL,
    time_unit = NULL
  )
  if (identical(task_type, "survival")) {
    horizons <- vapply(members, `[[`, numeric(1L), "prediction_horizon")
    time_units <- vapply(members, `[[`, character(1L), "time_unit")
    if (!all(horizons == horizons[[1L]]) ||
        !all(time_units == time_units[[1L]])) {
      sb_abort(
        paste0(
          "All survival cohorts must use the same time unit and prediction ",
          "horizon to meta-analyse calibration slope."
        ),
        "stormbreaker_external_meta_target_mismatch"
      )
    }
    metric_target$prediction_horizon <- horizons[[1L]]
    metric_target$time_unit <- time_units[[1L]]
  }

  estimand <- if (identical(task_type, "binary")) {
    paste0(
      "Random-effects mean calibration slope for one locked binary model ",
      "across an assumed distribution of independent exchangeable cohorts"
    )
  } else {
    paste0(
      "Random-effects mean calibration slope for one locked survival model ",
      "at ", format(metric_target$prediction_horizon), " ",
      metric_target$time_unit,
      " across an assumed distribution of independent exchangeable cohorts"
    )
  }

  if (!requireNamespace("metafor", quietly = TRUE)) {
    sb_abort(
      paste0(
        "Package `metafor` is required for external calibration-slope ",
        "meta-analysis. Install it with `install.packages(\"metafor\")`."
      ),
      "stormbreaker_missing_optional_dependency"
    )
  }

  cohorts <- data.frame(
    cohort_id = cohort_ids,
    n_patients = vapply(members, `[[`, integer(1L), "n"),
    point_estimate = vapply(members, `[[`, numeric(1L), "point_estimate"),
    standard_error = vapply(members, `[[`, numeric(1L), "standard_error"),
    n_bootstrap = vapply(members, `[[`, integer(1L), "reps"),
    stringsAsFactors = FALSE
  )
  fit <- metafor::rma.uni(
    yi = cohorts$point_estimate,
    sei = cohorts$standard_error,
    method = "REML",
    test = "knha",
    level = 100 * confidence_level,
    slab = cohort_ids
  )
  prediction <- stats::predict(fit)
  k <- as.integer(nrow(cohorts))
  summary <- data.frame(
    metric = "calibration_slope",
    k = k,
    estimate = as.numeric(fit$b),
    lower = as.numeric(fit$ci.lb),
    upper = as.numeric(fit$ci.ub),
    prediction_lower = as.numeric(prediction$pi.lb),
    prediction_upper = as.numeric(prediction$pi.ub),
    tau2 = as.numeric(fit$tau2),
    tau = sqrt(as.numeric(fit$tau2)),
    i2_percent = as.numeric(fit$I2),
    stringsAsFactors = FALSE
  )
  caution <- paste0(
    if (k < 10L) {
      paste0(
        "Only ", k, " cohorts were supplied; heterogeneity and the ",
        "between-cohort prediction interval are especially uncertain. "
      )
    } else {
      paste0(
        "The prediction interval assumes independent, exchangeable cohorts ",
        "and a correctly specified random-effects model. "
      )
    },
    "A pooled slope near 1 does not establish mean or moderate calibration; ",
    "inspect O:E and calibration curves separately."
  )

  structure(
    list(
      status = "success",
      task_type = task_type,
      pipeline_id = pipeline_ids[[1L]],
      learner = learners[[1L]],
      metric = "calibration_slope",
      metric_target = metric_target,
      confidence_level = as.numeric(confidence_level),
      estimand = estimand,
      interpretation = list(
        target = 1,
        below_one = "predictions_too_extreme_under_recalibration_model",
        above_one = "predictions_too_narrow_under_recalibration_model"
      ),
      method = list(
        model = "random_effects",
        effect_scale = "original_calibration_slope",
        within_cohort_standard_error =
          "sample_sd_of_patient_bootstrap_slope_replicates",
        tau2_estimator = "REML",
        inference = "Knapp-Hartung",
        back_transformation = "none"
      ),
      provenance = sb_external_meta_provenance(),
      cohorts = cohorts,
      summary = summary,
      prediction_interval_caution = caution,
      bootstraps = bootstraps
    ),
    class = "stormbreaker_external_calibration_slope_meta"
  )
}

#' @export
print.stormbreaker_external_calibration_slope_meta <- function(x, ...) {
  confidence <- 100 * x$confidence_level
  summary <- x$summary[1L, , drop = FALSE]
  cat("<stormbreaker_external_calibration_slope_meta:", x$status, ">\n")
  cat("fixed-model original-scale calibration-slope meta-analysis\n")
  cat(
    "Metric: calibration slope | Cohorts:", summary$k,
    "| Task:", x$task_type, "\n"
  )
  if (isTRUE(x$metric_target$has_horizon)) {
    cat(
      "Target:", format(x$metric_target$prediction_horizon),
      x$metric_target$time_unit, "\n"
    )
  }
  cat(
    "Summary slope:", format(summary$estimate),
    paste0("| ", format(confidence), "% CI:"),
    format(summary$lower), "to", format(summary$upper), "\n"
  )
  cat(
    paste0(format(confidence), "% prediction interval:"),
    format(summary$prediction_lower), "to",
    format(summary$prediction_upper), "\n"
  )
  cat(
    "Interpretation: target 1; a slope near 1 does not establish mean or ",
    "moderate calibration.\n",
    sep = ""
  )
  cat(
    "Scope: one univariate performance synthesis; no model ranking or ",
    "recalibration.\n",
    sep = ""
  )
  cat("Caution:", x$prediction_interval_caution, "\n")
  invisible(x)
}

# ---- meta-external-mean-calibration ----

sb_external_calibration_meta_target <- function(
    bootstrap,
    cohort_id,
    task_type,
    estimand_label) {
  if (identical(task_type, "binary")) {
    return(list(prediction_horizon = NULL, time_unit = NULL))
  }
  target <- bootstrap[["prediction_target"]]
  valid_target <- is.list(target) &&
    identical(target[["target_status"]], "complete") &&
    identical(target[["has_horizon"]], TRUE) &&
    is.numeric(target[["prediction_horizon"]]) &&
    length(target[["prediction_horizon"]]) == 1L &&
    !is.na(target[["prediction_horizon"]]) &&
    is.finite(target[["prediction_horizon"]]) &&
    target[["prediction_horizon"]] > 0 &&
    sb_external_scalar_character(target[["time_unit"]]) &&
    target[["time_unit"]] %in% c("day", "month", "year")
  if (!valid_target) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` requires a complete locked horizon and ",
        "time unit to meta-analyse ", estimand_label, "."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }
  list(
    prediction_horizon = as.numeric(target[["prediction_horizon"]]),
    time_unit = target[["time_unit"]]
  )
}

sb_external_mean_calibration_meta_member <- function(
    bootstrap,
    cohort_id) {
  if (!inherits(bootstrap, "stormbreaker_external_bootstrap") ||
      !is.list(bootstrap)) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` is not a result from `sb_bootstrap_external()`."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  task_type <- bootstrap[["task_type"]]
  pipeline_id <- bootstrap[["pipeline_id"]]
  learner <- bootstrap[["learner"]]
  n <- bootstrap[["n"]]
  reps <- bootstrap[["reps"]]
  valid_metadata <- sb_external_scalar_character(task_type) &&
    task_type %in% c("binary", "survival") &&
    sb_external_scalar_character(pipeline_id) &&
    sb_external_scalar_character(learner) &&
    is.numeric(n) && length(n) == 1L && !is.na(n) && is.finite(n) &&
    n >= 1 && n <= .Machine$integer.max && n == as.integer(n) &&
    is.numeric(reps) && length(reps) == 1L && !is.na(reps) &&
    is.finite(reps) && reps >= 2 && reps <= .Machine$integer.max &&
    reps == as.integer(reps)
  if (!valid_metadata) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has incomplete external-bootstrap metadata."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  metric_target <- sb_external_calibration_meta_target(
    bootstrap,
    cohort_id,
    task_type,
    "survival mean calibration"
  )
  point <- bootstrap[["point_estimates"]]
  replicates <- bootstrap[["replicates"]]
  valid_tables <- is.data.frame(point) &&
    all(c("metric", "estimate", "status") %in% names(point)) &&
    is.character(point$metric) && is.numeric(point$estimate) &&
    is.character(point$status) &&
    is.data.frame(replicates) &&
    all(c("replicate", "metric", "estimate", "status") %in%
      names(replicates)) &&
    is.numeric(replicates$replicate) &&
    is.character(replicates$metric) && is.numeric(replicates$estimate) &&
    is.character(replicates$status)
  if (!valid_tables) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` does not retain usable bootstrap metric tables."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  metric <- "mean_calibration_ratio"
  point_rows <- point[
    !is.na(point$metric) & point$metric == metric,
    ,
    drop = FALSE
  ]
  replicate_rows <- replicates[
    !is.na(replicates$metric) & replicates$metric == metric,
    ,
    drop = FALSE
  ]
  expected_replicates <- seq_len(as.integer(reps))
  replicate_ids_complete <- nrow(replicate_rows) == as.integer(reps) &&
    !anyNA(replicate_rows$replicate) &&
    all(is.finite(replicate_rows$replicate)) &&
    all(replicate_rows$replicate >= 1 &
      replicate_rows$replicate <= as.integer(reps)) &&
    all(replicate_rows$replicate == floor(replicate_rows$replicate)) &&
    identical(
      sort(as.integer(replicate_rows$replicate)),
      expected_replicates
    )
  point_available <- nrow(point_rows) == 1L &&
    identical(point_rows$status[[1L]], "success") &&
    length(point_rows$estimate) == 1L &&
    !is.na(point_rows$estimate[[1L]]) &&
    is.finite(point_rows$estimate[[1L]]) &&
    point_rows$estimate[[1L]] > 0
  replicates_complete <- replicate_ids_complete &&
    !anyNA(replicate_rows$status) &&
    all(replicate_rows$status == "success") &&
    !anyNA(replicate_rows$estimate) &&
    all(is.finite(replicate_rows$estimate)) &&
    all(replicate_rows$estimate > 0)
  if (!point_available || !replicates_complete) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` requires one finite successful `", metric,
        "` point estimate and all ", as.integer(reps),
        " successful bootstrap estimates, numbered exactly once from 1 ",
        "through `reps`, and strictly greater than zero. No failed replicate ",
        "is dropped and no continuity correction is applied."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  point_estimate <- point_rows$estimate[[1L]]
  log_replicates <- log(replicate_rows$estimate)
  log_se <- stats::sd(log_replicates)
  if (!is.finite(log_se) || log_se <= 0) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has zero or non-finite bootstrap uncertainty for `", metric,
        "` on the log scale."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  list(
    task_type = task_type,
    pipeline_id = pipeline_id,
    learner = learner,
    n = as.integer(n),
    reps = as.integer(reps),
    point_estimate = as.numeric(point_estimate),
    log_estimate = log(point_estimate),
    log_se = as.numeric(log_se),
    prediction_horizon = metric_target$prediction_horizon,
    time_unit = metric_target$time_unit
  )
}

#' @rdname sb_meta_external_discrimination
#' @export
sb_meta_external_mean_calibration <- function(
    bootstraps,
    confidence_level = 0.95) {
  valid_confidence <- is.numeric(confidence_level) &&
    length(confidence_level) == 1L && !is.na(confidence_level) &&
    is.finite(confidence_level) && confidence_level > 0 &&
    confidence_level < 1
  if (!valid_confidence) {
    sb_abort(
      "`confidence_level` must be one finite number strictly between 0 and 1.",
      "stormbreaker_bad_confidence_level"
    )
  }

  valid_input <- is.list(bootstraps) &&
    !inherits(bootstraps, "stormbreaker_external_bootstrap")
  if (!valid_input) {
    sb_abort(
      paste0(
        "`bootstraps` must be a named list of results from ",
        "`sb_bootstrap_external()`."
      ),
      "stormbreaker_bad_external_meta"
    )
  }
  if (length(bootstraps) < 2L) {
    sb_abort(
      "At least two independent external cohorts are required.",
      "stormbreaker_insufficient_external_cohorts"
    )
  }

  cohort_ids <- names(bootstraps)
  valid_ids <- !is.null(cohort_ids) &&
    length(cohort_ids) == length(bootstraps) &&
    !anyNA(cohort_ids) && all(nzchar(trimws(cohort_ids))) &&
    anyDuplicated(cohort_ids) == 0L
  if (!valid_ids) {
    sb_abort(
      "Outer list names must be unique, non-blank cohort IDs.",
      "stormbreaker_bad_external_meta"
    )
  }

  members <- lapply(seq_along(bootstraps), function(index) {
    sb_external_mean_calibration_meta_member(
      bootstraps[[index]],
      cohort_ids[[index]]
    )
  })
  task_types <- vapply(members, `[[`, character(1L), "task_type")
  if (length(unique(task_types)) != 1L) {
    sb_abort(
      "All cohorts must use the same task type.",
      "stormbreaker_external_meta_task_mismatch"
    )
  }
  pipeline_ids <- vapply(members, `[[`, character(1L), "pipeline_id")
  learners <- vapply(members, `[[`, character(1L), "learner")
  if (length(unique(pipeline_ids)) != 1L ||
      length(unique(learners)) != 1L) {
    sb_abort(
      paste0(
        "All cohorts must retain the same pipeline ID and learner label. ",
        "These labels do not prove identical fitted state or formula."
      ),
      "stormbreaker_external_meta_model_mismatch"
    )
  }

  task_type <- task_types[[1L]]
  metric_target <- list(
    has_horizon = identical(task_type, "survival"),
    prediction_horizon = NULL,
    time_unit = NULL
  )
  if (identical(task_type, "survival")) {
    horizons <- vapply(members, `[[`, numeric(1L), "prediction_horizon")
    time_units <- vapply(members, `[[`, character(1L), "time_unit")
    if (!all(horizons == horizons[[1L]]) ||
        !all(time_units == time_units[[1L]])) {
      sb_abort(
        paste0(
          "All survival cohorts must use the same time unit and prediction ",
          "horizon to meta-analyse mean calibration."
        ),
        "stormbreaker_external_meta_target_mismatch"
      )
    }
    metric_target$prediction_horizon <- horizons[[1L]]
    metric_target$time_unit <- time_units[[1L]]
  }

  estimand <- if (identical(task_type, "binary")) {
    paste0(
      "Exponentiated random-effects mean of log observed-to-expected event ",
      "ratio for one locked binary model across an assumed distribution of ",
      "independent exchangeable cohorts"
    )
  } else {
    paste0(
      "Exponentiated random-effects mean of log observed-to-expected event ",
      "ratio for one locked survival model at ",
      format(metric_target$prediction_horizon), " ", metric_target$time_unit,
      " across an assumed distribution of independent exchangeable cohorts"
    )
  }

  if (!requireNamespace("metafor", quietly = TRUE)) {
    sb_abort(
      paste0(
        "Package `metafor` is required for external mean-calibration ",
        "meta-analysis. Install it with `install.packages(\"metafor\")`."
      ),
      "stormbreaker_missing_optional_dependency"
    )
  }

  cohorts <- data.frame(
    cohort_id = cohort_ids,
    n_patients = vapply(members, `[[`, integer(1L), "n"),
    point_estimate = vapply(members, `[[`, numeric(1L), "point_estimate"),
    log_estimate = vapply(members, `[[`, numeric(1L), "log_estimate"),
    log_se = vapply(members, `[[`, numeric(1L), "log_se"),
    n_bootstrap = vapply(members, `[[`, integer(1L), "reps"),
    stringsAsFactors = FALSE
  )
  fit <- metafor::rma.uni(
    yi = cohorts$log_estimate,
    sei = cohorts$log_se,
    method = "REML",
    test = "knha",
    level = 100 * confidence_level,
    slab = cohort_ids
  )
  prediction <- stats::predict(fit)
  k <- as.integer(nrow(cohorts))
  summary <- data.frame(
    metric = "mean_calibration_ratio",
    k = k,
    estimate = exp(as.numeric(fit$b)),
    lower = exp(as.numeric(fit$ci.lb)),
    upper = exp(as.numeric(fit$ci.ub)),
    prediction_lower = exp(as.numeric(prediction$pi.lb)),
    prediction_upper = exp(as.numeric(prediction$pi.ub)),
    tau2_log = as.numeric(fit$tau2),
    tau_log = sqrt(as.numeric(fit$tau2)),
    i2_percent = as.numeric(fit$I2),
    stringsAsFactors = FALSE
  )
  caution <- paste0(
    if (k < 10L) {
      paste0(
        "Only ", k, " cohorts were supplied; heterogeneity and the ",
        "between-cohort prediction interval are especially uncertain. "
      )
    } else {
      paste0(
        "The prediction interval assumes independent, exchangeable cohorts ",
        "and a correctly specified random-effects model. "
      )
    },
    "A pooled O:E near 1 can hide opposite miscalibration across cohorts or ",
    "within the prediction range."
  )

  structure(
    list(
      status = "success",
      task_type = task_type,
      pipeline_id = pipeline_ids[[1L]],
      learner = learners[[1L]],
      metric = "mean_calibration_ratio",
      metric_target = metric_target,
      confidence_level = as.numeric(confidence_level),
      estimand = estimand,
      interpretation = list(
        target = 1,
        below_one = "average_overprediction",
        above_one = "average_underprediction"
      ),
      method = list(
        model = "random_effects",
        effect_scale = "log_observed_expected_ratio",
        within_cohort_standard_error =
          "sample_sd_of_log_patient_bootstrap_replicates",
        tau2_estimator = "REML",
        inference = "Knapp-Hartung",
        back_transformation = "exponential"
      ),
      provenance = sb_external_meta_provenance(),
      cohorts = cohorts,
      summary = summary,
      prediction_interval_caution = caution,
      bootstraps = bootstraps
    ),
    class = "stormbreaker_external_mean_calibration_meta"
  )
}

#' @export
print.stormbreaker_external_mean_calibration_meta <- function(x, ...) {
  confidence <- 100 * x$confidence_level
  summary <- x$summary[1L, , drop = FALSE]
  cat("<stormbreaker_external_mean_calibration_meta:", x$status, ">\n")
  cat("fixed-model external mean calibration meta-analysis\n")
  cat("Metric: O:E | Cohorts:", summary$k, "| Task:", x$task_type, "\n")
  if (isTRUE(x$metric_target$has_horizon)) {
    cat(
      "Target:", format(x$metric_target$prediction_horizon),
      x$metric_target$time_unit, "\n"
    )
  }
  cat(
    "Summary O:E:", format(summary$estimate),
    paste0("| ", format(confidence), "% CI:"),
    format(summary$lower), "to", format(summary$upper), "\n"
  )
  cat(
    paste0(format(confidence), "% prediction interval:"),
    format(summary$prediction_lower), "to",
    format(summary$prediction_upper), "\n"
  )
  cat(
    "Interpretation: O:E below 1 means overprediction; above 1 means ",
    "underprediction.\n",
    sep = ""
  )
  cat(
    "Scope: mean calibration only; this does not pool calibration slope or ",
    "curves.\n",
    sep = ""
  )
  cat("Caution:", x$prediction_interval_caution, "\n")
  invisible(x)
}

# ---- meta-external-prediction-error-difference ----

#' @rdname sb_meta_external_discrimination
#' @export
sb_meta_external_prediction_error_difference <- function(
    comparisons,
    candidate_pipeline_id,
    confidence_level = 0.95,
    metric = NULL) {
  core <- sb_external_paired_difference_meta_core(
    comparisons = comparisons,
    candidate_pipeline_id = candidate_pipeline_id,
    confidence_level = confidence_level,
    metric = metric,
    metric_family = "prediction_error"
  )
  metric <- core$metric
  metric_target <- core$metric_target
  estimand <- switch(
    metric,
    brier_score = paste0(
      "Random-effects mean paired candidate-minus-reference binary Brier ",
      "score difference across an assumed distribution of independent ",
      "exchangeable cohorts"
    ),
    log_loss = paste0(
      "Random-effects mean paired candidate-minus-reference binary log-loss ",
      "difference across an assumed distribution of independent ",
      "exchangeable cohorts"
    ),
    ipcw_brier_score = paste0(
      "Random-effects mean paired candidate-minus-reference IPCW Brier ",
      "score difference at ", format(metric_target$prediction_horizon), " ",
      metric_target$time_unit,
      " across an assumed distribution of independent exchangeable cohorts"
    )
  )
  k <- core$summary$k[[1L]]
  caution <- paste0(
    if (k < 10L) {
      paste0(
        "Only ", k, " cohorts were supplied; heterogeneity and the ",
        "between-cohort prediction interval are especially uncertain. "
      )
    } else {
      paste0(
        "The prediction interval assumes independent, exchangeable cohorts ",
        "and a correctly specified random-effects model. "
      )
    },
    "Raw-scale normal-model intervals are not clipped. A negative difference ",
    "means lower candidate prediction error only, not established clinical ",
    "utility or uniformly better calibration and discrimination."
  )

  structure(
    list(
      status = "success",
      task_type = core$task_type,
      reference_pipeline_id = core$reference_pipeline_id,
      candidate_pipeline_id = core$candidate_pipeline_id,
      contrast = "candidate_minus_reference",
      metric = metric,
      metric_target = metric_target,
      confidence_level = core$confidence_level,
      estimand = estimand,
      interpretation = list(
        null = 0,
        negative = "lower_candidate_prediction_error",
        positive = "lower_reference_prediction_error"
      ),
      method = list(
        model = "random_effects",
        effect_scale = "raw_candidate_minus_reference_prediction_error",
        within_cohort_standard_error =
          "sample_sd_of_paired_patient_bootstrap_differences",
        tau2_estimator = "REML",
        inference = "Knapp-Hartung",
        back_transformation = "none"
      ),
      provenance = sb_external_meta_provenance(),
      cohorts = core$cohorts,
      summary = core$summary,
      prediction_interval_caution = caution,
      comparisons = core$comparisons
    ),
    class = "stormbreaker_external_prediction_error_difference_meta"
  )
}

#' @export
print.stormbreaker_external_prediction_error_difference_meta <- function(
    x,
    ...) {
  confidence <- 100 * x$confidence_level
  summary <- x$summary[1L, , drop = FALSE]
  cat(
    "<stormbreaker_external_prediction_error_difference_meta:",
    x$status,
    ">\n"
  )
  cat("paired external prediction error difference meta-analysis\n")
  cat(
    "Contrast:", x$candidate_pipeline_id, "-", x$reference_pipeline_id,
    "(candidate - reference)\n"
  )
  cat(
    "Metric:", x$metric,
    "| Cohorts:", summary$k,
    "| Task:", x$task_type, "\n"
  )
  if (isTRUE(x$metric_target$has_horizon)) {
    cat(
      "Target:", format(x$metric_target$prediction_horizon),
      x$metric_target$time_unit, "\n"
    )
  }
  cat(
    sprintf(
      "Mean difference: %.3f | %.1f%% CI: %.3f to %.3f\n",
      summary$estimate,
      confidence,
      summary$lower,
      summary$upper
    )
  )
  cat(
    sprintf(
      "Between-cohort %.1f%% prediction interval: %.3f to %.3f\n",
      confidence,
      summary$prediction_lower,
      summary$prediction_upper
    )
  )
  cat(
    sprintf(
      "tau^2: %.4f | tau: %.4f | I^2: %.1f%%\n",
      summary$tau2,
      summary$tau,
      summary$i2_percent
    )
  )
  cat(x$prediction_interval_caution, "\n")
  cat(
    paste0(
      "Source comparisons are retained; no ranking or model selection is ",
      "performed. Model identity and cohort independence remain caller ",
      "responsibilities.\n"
    )
  )
  invisible(x)
}

# ---- meta-discrimination ----

# ---- meta-external-common ----

sb_external_meta_provenance <- function() {
  list(
    package = "stormbreaker",
    package_version = as.character(utils::packageVersion("stormbreaker")),
    r_version = R.version.string,
    metafor_version = as.character(utils::packageVersion("metafor"))
  )
}

sb_external_prediction_error_meta_metric <- function(task_type, metric) {
  if (is.null(metric)) {
    return(if (identical(task_type, "binary")) {
      "brier_score"
    } else {
      "ipcw_brier_score"
    })
  }
  if (!sb_external_scalar_character(metric)) {
    sb_abort(
      "`metric` must be NULL or one non-empty metric name.",
      "stormbreaker_bad_external_meta_metric"
    )
  }
  allowed <- if (identical(task_type, "binary")) {
    c("brier_score", "log_loss")
  } else {
    "ipcw_brier_score"
  }
  if (!metric %in% allowed) {
    sb_abort(
      paste0(
        "For a ", task_type, " task, `metric` must be one of: ",
        paste0("`", allowed, "`", collapse = ", "), "."
      ),
      "stormbreaker_bad_external_meta_metric"
    )
  }
  metric
}

sb_external_prediction_error_meta_target <- function(
    comparison,
    cohort_id,
    metric) {
  if (!identical(metric, "ipcw_brier_score")) {
    return(list(prediction_horizon = NULL, time_unit = NULL))
  }
  target <- comparison[["prediction_target"]]
  valid_target <- is.list(target) &&
    identical(target[["target_status"]], "complete") &&
    identical(target[["has_horizon"]], TRUE) &&
    is.numeric(target[["prediction_horizon"]]) &&
    length(target[["prediction_horizon"]]) == 1L &&
    !is.na(target[["prediction_horizon"]]) &&
    is.finite(target[["prediction_horizon"]]) &&
    target[["prediction_horizon"]] > 0 &&
    sb_external_scalar_character(target[["time_unit"]]) &&
    target[["time_unit"]] %in% c("day", "month", "year")
  if (!valid_target) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` requires a complete locked horizon and ",
        "time unit to meta-analyse `", metric, "` differences."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }
  list(
    prediction_horizon = as.numeric(target[["prediction_horizon"]]),
    time_unit = target[["time_unit"]]
  )
}

sb_external_paired_metric_estimates_valid <- function(values, metric) {
  if (length(values) == 0L || anyNA(values) || !all(is.finite(values))) {
    return(FALSE)
  }
  if (metric %in% c(
    "roc_auc", "harrell_c", "time_auc", "uno_c", "brier_score"
  )) {
    return(all(values >= 0 & values <= 1))
  }
  all(values >= 0)
}

sb_external_paired_difference_estimates_valid <- function(values, metric) {
  if (length(values) == 0L || anyNA(values) || !all(is.finite(values))) {
    return(FALSE)
  }
  if (metric %in% c(
    "roc_auc", "harrell_c", "time_auc", "uno_c", "brier_score"
  )) {
    return(all(values >= -1 & values <= 1))
  }
  TRUE
}

sb_external_paired_difference_meta_member <- function(
    comparison,
    cohort_id,
    candidate_pipeline_id,
    metric = NULL,
    metric_family = c("discrimination", "prediction_error")) {
  metric_family <- match.arg(metric_family)
  if (!inherits(
    comparison,
    "stormbreaker_external_bootstrap_comparison"
  ) || !is.list(comparison)) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` is not a result from ",
        "`sb_bootstrap_external_comparison()`."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  task_type <- comparison[["task_type"]]
  reference_pipeline_id <- comparison[["reference_pipeline_id"]]
  n <- comparison[["n"]]
  reps <- comparison[["reps"]]
  valid_metadata <- sb_external_scalar_character(comparison[["status"]]) &&
    comparison[["status"]] %in% c("success", "incomplete") &&
    sb_external_scalar_character(task_type) &&
    task_type %in% c("binary", "survival") &&
    sb_external_scalar_character(reference_pipeline_id) &&
    is.numeric(n) && length(n) == 1L && !is.na(n) && is.finite(n) &&
    n >= 1 && n <= .Machine$integer.max && n == as.integer(n) &&
    is.numeric(reps) && length(reps) == 1L && !is.na(reps) &&
    is.finite(reps) && reps >= 2 && reps <= .Machine$integer.max &&
    reps == as.integer(reps)
  if (!valid_metadata) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has incomplete paired-bootstrap metadata."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  candidate_status <- comparison[["candidate_status"]]
  valid_candidate_table <- is.data.frame(candidate_status) &&
    all(c("candidate_pipeline_id", "status") %in% names(candidate_status)) &&
    is.character(candidate_status$candidate_pipeline_id) &&
    is.character(candidate_status$status)
  if (!valid_candidate_table) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` does not retain usable candidate-status metadata."
      ),
      "stormbreaker_bad_external_meta"
    )
  }
  candidate_rows <- candidate_status[
    !is.na(candidate_status$candidate_pipeline_id) &
      candidate_status$candidate_pipeline_id == candidate_pipeline_id,
    ,
    drop = FALSE
  ]
  candidate_available <- nrow(candidate_rows) == 1L &&
    identical(candidate_rows$status[[1L]], "success")
  if (!candidate_available) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` requires one successful candidate `",
        candidate_pipeline_id, "` in the paired bootstrap comparison."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  if (identical(metric_family, "discrimination")) {
    metric <- sb_external_discrimination_meta_metric(task_type, metric)
    metric_target <- sb_external_discrimination_meta_target(
      comparison,
      cohort_id,
      metric
    )
    expected_direction <- "higher_is_better"
  } else {
    metric <- sb_external_prediction_error_meta_metric(task_type, metric)
    metric_target <- sb_external_prediction_error_meta_target(
      comparison,
      cohort_id,
      metric
    )
    expected_direction <- "lower_is_better"
  }
  metrics <- comparison[["metrics"]]
  if (!is.character(metrics) || !metric %in% metrics) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` does not retain the requested paired `",
        metric, "` metric."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  point <- comparison[["point_differences"]]
  replicates <- comparison[["differences"]]
  required_point <- c(
    "reference_pipeline_id", "candidate_pipeline_id", "contrast", "metric",
    "metric_direction", "reference_estimate", "candidate_estimate",
    "estimate", "status"
  )
  required_replicates <- c("replicate", required_point)
  valid_tables <- is.data.frame(point) &&
    all(required_point %in% names(point)) &&
    is.character(point$reference_pipeline_id) &&
    is.character(point$candidate_pipeline_id) &&
    is.character(point$contrast) && is.character(point$metric) &&
    is.character(point$metric_direction) &&
    is.numeric(point$reference_estimate) &&
    is.numeric(point$candidate_estimate) && is.numeric(point$estimate) &&
    is.character(point$status) && is.data.frame(replicates) &&
    all(required_replicates %in% names(replicates)) &&
    is.numeric(replicates$replicate) &&
    is.character(replicates$reference_pipeline_id) &&
    is.character(replicates$candidate_pipeline_id) &&
    is.character(replicates$contrast) && is.character(replicates$metric) &&
    is.character(replicates$metric_direction) &&
    is.numeric(replicates$reference_estimate) &&
    is.numeric(replicates$candidate_estimate) &&
    is.numeric(replicates$estimate) && is.character(replicates$status)
  if (!valid_tables) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` does not retain usable paired-difference tables."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  point_rows <- point[
    !is.na(point$reference_pipeline_id) &
      point$reference_pipeline_id == reference_pipeline_id &
      !is.na(point$candidate_pipeline_id) &
      point$candidate_pipeline_id == candidate_pipeline_id &
      !is.na(point$metric) & point$metric == metric,
    ,
    drop = FALSE
  ]
  replicate_rows <- replicates[
    !is.na(replicates$reference_pipeline_id) &
      replicates$reference_pipeline_id == reference_pipeline_id &
      !is.na(replicates$candidate_pipeline_id) &
      replicates$candidate_pipeline_id == candidate_pipeline_id &
      !is.na(replicates$metric) & replicates$metric == metric,
    ,
    drop = FALSE
  ]
  point_difference <-
    point_rows$candidate_estimate - point_rows$reference_estimate
  replicate_differences <-
    replicate_rows$candidate_estimate - replicate_rows$reference_estimate
  point_arithmetic_consistent <- identical(
    unname(point_rows$estimate),
    unname(point_difference)
  )
  replicate_arithmetic_consistent <- identical(
    unname(replicate_rows$estimate),
    unname(replicate_differences)
  )
  expected_replicates <- seq_len(as.integer(reps))
  replicate_ids_complete <- nrow(replicate_rows) == as.integer(reps) &&
    !anyNA(replicate_rows$replicate) &&
    all(is.finite(replicate_rows$replicate)) &&
    all(replicate_rows$replicate >= 1 &
      replicate_rows$replicate <= as.integer(reps)) &&
    all(replicate_rows$replicate == floor(replicate_rows$replicate)) &&
    identical(
      sort(as.integer(replicate_rows$replicate)),
      expected_replicates
    )
  point_available <- nrow(point_rows) == 1L &&
    point_arithmetic_consistent &&
    identical(point_rows$status[[1L]], "success") &&
    identical(point_rows$contrast[[1L]], "candidate_minus_reference") &&
    identical(point_rows$metric_direction[[1L]], expected_direction) &&
    sb_external_paired_metric_estimates_valid(
      c(point_rows$reference_estimate, point_rows$candidate_estimate),
      metric
    ) &&
    sb_external_paired_difference_estimates_valid(
      point_difference,
      metric
    )
  replicates_complete <- replicate_ids_complete &&
    replicate_arithmetic_consistent &&
    !anyNA(replicate_rows$status) &&
    all(replicate_rows$status == "success") &&
    !anyNA(replicate_rows$contrast) &&
    all(replicate_rows$contrast == "candidate_minus_reference") &&
    !anyNA(replicate_rows$metric_direction) &&
    all(replicate_rows$metric_direction == expected_direction) &&
    sb_external_paired_metric_estimates_valid(
      c(
        replicate_rows$reference_estimate,
        replicate_rows$candidate_estimate
      ),
      metric
    ) &&
    sb_external_paired_difference_estimates_valid(
      replicate_differences,
      metric
    )
  if (!point_available || !replicates_complete) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` requires one finite successful paired `",
        metric, "` point difference and all ", as.integer(reps),
        " paired bootstrap differences, numbered exactly once from 1 through ",
        "`reps`. Every saved difference must equal candidate minus reference; ",
        "no failed replicate is dropped and no difference is clipped."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  standard_error <- stats::sd(replicate_differences)
  if (!is.finite(standard_error) || standard_error <= 0) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has zero or non-finite paired bootstrap uncertainty for `", metric,
        "`."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  list(
    task_type = task_type,
    reference_pipeline_id = reference_pipeline_id,
    candidate_pipeline_id = candidate_pipeline_id,
    metric = metric,
    n = as.integer(n),
    reps = as.integer(reps),
    reference_point_estimate = point_rows$reference_estimate[[1L]],
    candidate_point_estimate = point_rows$candidate_estimate[[1L]],
    point_difference = as.numeric(point_difference[[1L]]),
    standard_error = as.numeric(standard_error),
    prediction_horizon = metric_target$prediction_horizon,
    time_unit = metric_target$time_unit
  )
}

sb_external_paired_difference_meta_core <- function(
    comparisons,
    candidate_pipeline_id,
    confidence_level,
    metric,
    metric_family = c("discrimination", "prediction_error")) {
  metric_family <- match.arg(metric_family)
  valid_confidence <- is.numeric(confidence_level) &&
    length(confidence_level) == 1L && !is.na(confidence_level) &&
    is.finite(confidence_level) && confidence_level > 0 &&
    confidence_level < 1
  if (!valid_confidence) {
    sb_abort(
      "`confidence_level` must be one finite number strictly between 0 and 1.",
      "stormbreaker_bad_confidence_level"
    )
  }
  if (!sb_external_scalar_character(candidate_pipeline_id)) {
    sb_abort(
      "`candidate_pipeline_id` must be one non-empty pipeline ID.",
      "stormbreaker_bad_external_meta"
    )
  }

  valid_input <- is.list(comparisons) && !inherits(
    comparisons,
    "stormbreaker_external_bootstrap_comparison"
  )
  if (!valid_input) {
    sb_abort(
      paste0(
        "`comparisons` must be a named list of results from ",
        "`sb_bootstrap_external_comparison()`."
      ),
      "stormbreaker_bad_external_meta"
    )
  }
  if (length(comparisons) < 2L) {
    sb_abort(
      "At least two independent external cohorts are required.",
      "stormbreaker_insufficient_external_cohorts"
    )
  }

  cohort_ids <- names(comparisons)
  valid_ids <- !is.null(cohort_ids) &&
    length(cohort_ids) == length(comparisons) &&
    !anyNA(cohort_ids) && all(nzchar(trimws(cohort_ids))) &&
    anyDuplicated(cohort_ids) == 0L
  if (!valid_ids) {
    sb_abort(
      "Outer list names must be unique, non-blank cohort IDs.",
      "stormbreaker_bad_external_meta"
    )
  }

  members <- lapply(seq_along(comparisons), function(index) {
    sb_external_paired_difference_meta_member(
      comparisons[[index]],
      cohort_ids[[index]],
      candidate_pipeline_id,
      metric,
      metric_family
    )
  })
  task_types <- vapply(members, `[[`, character(1L), "task_type")
  if (length(unique(task_types)) != 1L) {
    sb_abort(
      "All cohorts must use the same task type.",
      "stormbreaker_external_meta_task_mismatch"
    )
  }
  reference_ids <- vapply(
    members,
    `[[`,
    character(1L),
    "reference_pipeline_id"
  )
  if (length(unique(reference_ids)) != 1L) {
    sb_abort(
      paste0(
        "All cohorts must retain the same reference and candidate pipeline ",
        "IDs. Matching labels do not prove identical fitted states."
      ),
      "stormbreaker_external_meta_model_mismatch"
    )
  }

  metric <- members[[1L]]$metric
  has_horizon <- metric %in% c("time_auc", "uno_c", "ipcw_brier_score")
  metric_target <- list(
    has_horizon = has_horizon,
    prediction_horizon = NULL,
    time_unit = NULL
  )
  if (has_horizon) {
    horizons <- vapply(members, `[[`, numeric(1L), "prediction_horizon")
    time_units <- vapply(members, `[[`, character(1L), "time_unit")
    same_horizon <- all(horizons == horizons[[1L]])
    same_unit <- all(time_units == time_units[[1L]])
    if (!same_horizon || !same_unit) {
      sb_abort(
        paste0(
          "All cohorts must use the same time unit and prediction horizon ",
          "to meta-analyse paired `", metric, "` differences."
        ),
        "stormbreaker_external_meta_target_mismatch"
      )
    }
    metric_target$prediction_horizon <- horizons[[1L]]
    metric_target$time_unit <- time_units[[1L]]
  }

  if (!requireNamespace("metafor", quietly = TRUE)) {
    sb_abort(
      paste0(
        "Package `metafor` is required for paired external difference ",
        "meta-analysis. Install it with `install.packages(\"metafor\")`."
      ),
      "stormbreaker_missing_optional_dependency"
    )
  }

  cohorts <- data.frame(
    cohort_id = cohort_ids,
    n_patients = vapply(members, `[[`, integer(1L), "n"),
    reference_point_estimate = vapply(
      members,
      `[[`,
      numeric(1L),
      "reference_point_estimate"
    ),
    candidate_point_estimate = vapply(
      members,
      `[[`,
      numeric(1L),
      "candidate_point_estimate"
    ),
    point_difference = vapply(
      members,
      `[[`,
      numeric(1L),
      "point_difference"
    ),
    standard_error = vapply(members, `[[`, numeric(1L), "standard_error"),
    n_bootstrap = vapply(members, `[[`, integer(1L), "reps"),
    stringsAsFactors = FALSE
  )
  fit <- metafor::rma.uni(
    yi = cohorts$point_difference,
    sei = cohorts$standard_error,
    method = "REML",
    test = "knha",
    level = 100 * confidence_level,
    slab = cohort_ids
  )
  prediction <- stats::predict(fit)
  summary <- data.frame(
    metric = metric,
    contrast = "candidate_minus_reference",
    k = as.integer(nrow(cohorts)),
    estimate = as.numeric(fit$b),
    lower = as.numeric(fit$ci.lb),
    upper = as.numeric(fit$ci.ub),
    prediction_lower = as.numeric(prediction$pi.lb),
    prediction_upper = as.numeric(prediction$pi.ub),
    tau2 = as.numeric(fit$tau2),
    tau = sqrt(as.numeric(fit$tau2)),
    i2_percent = as.numeric(fit$I2),
    stringsAsFactors = FALSE
  )

  list(
    task_type = task_types[[1L]],
    reference_pipeline_id = reference_ids[[1L]],
    candidate_pipeline_id = candidate_pipeline_id,
    metric = metric,
    metric_target = metric_target,
    confidence_level = as.numeric(confidence_level),
    cohorts = cohorts,
    summary = summary,
    comparisons = comparisons
  )
}

# ---- meta-external-discrimination-difference ----

sb_external_discrimination_difference_meta_member <- function(
    comparison,
    cohort_id,
    candidate_pipeline_id,
    metric = NULL) {
  sb_external_paired_difference_meta_member(
    comparison = comparison,
    cohort_id = cohort_id,
    candidate_pipeline_id = candidate_pipeline_id,
    metric = metric,
    metric_family = "discrimination"
  )
}

#' @rdname sb_meta_external_discrimination
#' @export
sb_meta_external_discrimination_difference <- function(
    comparisons,
    candidate_pipeline_id,
    confidence_level = 0.95,
    metric = NULL) {
  core <- sb_external_paired_difference_meta_core(
    comparisons = comparisons,
    candidate_pipeline_id = candidate_pipeline_id,
    confidence_level = confidence_level,
    metric = metric,
    metric_family = "discrimination"
  )
  metric <- core$metric
  metric_target <- core$metric_target
  estimand <- switch(
    metric,
    roc_auc = paste0(
      "Random-effects mean paired candidate-minus-reference ROC AUC ",
      "difference across an assumed distribution of independent ",
      "exchangeable cohorts"
    ),
    harrell_c = paste0(
      "Random-effects mean paired candidate-minus-reference horizon-free ",
      "Harrell concordance difference across an assumed distribution of ",
      "independent exchangeable cohorts"
    ),
    time_auc = paste0(
      "Random-effects mean paired candidate-minus-reference cumulative/",
      "dynamic time-dependent AUC difference at ",
      format(metric_target$prediction_horizon), " ", metric_target$time_unit,
      " across an assumed distribution of independent exchangeable cohorts"
    ),
    uno_c = paste0(
      "Random-effects mean paired candidate-minus-reference truncated Uno ",
      "concordance difference at ",
      format(metric_target$prediction_horizon), " ", metric_target$time_unit,
      " across an assumed distribution of independent exchangeable cohorts"
    )
  )
  k <- core$summary$k[[1L]]
  caution <- paste0(
    if (k < 10L) {
      paste0(
        "Only ", k, " cohorts were supplied; heterogeneity and the ",
        "between-cohort prediction interval are especially uncertain. "
      )
    } else {
      paste0(
        "The prediction interval assumes independent, exchangeable cohorts ",
        "and a correctly specified random-effects model. "
      )
    },
    "Raw-scale normal-model intervals are not clipped to [-1, 1]. A positive ",
    "difference means higher candidate discrimination only, not better ",
    "calibration, clinical utility, or overall performance."
  )

  structure(
    list(
      status = "success",
      task_type = core$task_type,
      reference_pipeline_id = core$reference_pipeline_id,
      candidate_pipeline_id = core$candidate_pipeline_id,
      contrast = "candidate_minus_reference",
      metric = metric,
      metric_target = metric_target,
      confidence_level = core$confidence_level,
      estimand = estimand,
      interpretation = list(
        null = 0,
        positive = "higher_candidate_discrimination",
        negative = "higher_reference_discrimination"
      ),
      method = list(
        model = "random_effects",
        effect_scale = "raw_candidate_minus_reference_discrimination",
        within_cohort_standard_error =
          "sample_sd_of_paired_patient_bootstrap_differences",
        tau2_estimator = "REML",
        inference = "Knapp-Hartung",
        back_transformation = "none"
      ),
      provenance = sb_external_meta_provenance(),
      cohorts = core$cohorts,
      summary = core$summary,
      prediction_interval_caution = caution,
      comparisons = core$comparisons
    ),
    class = "stormbreaker_external_discrimination_difference_meta"
  )
}

#' @export
print.stormbreaker_external_discrimination_difference_meta <- function(
    x,
    ...) {
  confidence <- 100 * x$confidence_level
  summary <- x$summary[1L, , drop = FALSE]
  cat(
    "<stormbreaker_external_discrimination_difference_meta:",
    x$status,
    ">\n"
  )
  cat("paired external discrimination difference meta-analysis\n")
  cat(
    "Contrast:", x$candidate_pipeline_id, "-", x$reference_pipeline_id,
    "(candidate - reference)\n"
  )
  cat(
    "Metric:", x$metric,
    "| Cohorts:", summary$k,
    "| Task:", x$task_type, "\n"
  )
  if (isTRUE(x$metric_target$has_horizon)) {
    cat(
      "Target:", format(x$metric_target$prediction_horizon),
      x$metric_target$time_unit, "\n"
    )
  }
  cat(
    sprintf(
      "Mean difference: %.3f | %.1f%% CI: %.3f to %.3f\n",
      summary$estimate,
      confidence,
      summary$lower,
      summary$upper
    )
  )
  cat(
    sprintf(
      "Between-cohort %.1f%% prediction interval: %.3f to %.3f\n",
      confidence,
      summary$prediction_lower,
      summary$prediction_upper
    )
  )
  cat(
    sprintf(
      "tau^2: %.4f | tau: %.4f | I^2: %.1f%%\n",
      summary$tau2,
      summary$tau,
      summary$i2_percent
    )
  )
  cat(x$prediction_interval_caution, "\n")
  cat(
    paste0(
      "Source comparisons are retained; no ranking or model selection is ",
      "performed. Model identity and cohort independence remain caller ",
      "responsibilities.\n"
    )
  )
  invisible(x)
}

# ---- meta-external-discrimination ----

sb_external_discrimination_meta_metric <- function(task_type, metric) {
  if (is.null(metric)) {
    return(if (identical(task_type, "binary")) "roc_auc" else "harrell_c")
  }
  if (!sb_external_scalar_character(metric)) {
    sb_abort(
      "`metric` must be NULL or one non-empty metric name.",
      "stormbreaker_bad_external_meta_metric"
    )
  }
  allowed <- if (identical(task_type, "binary")) {
    "roc_auc"
  } else {
    c("harrell_c", "time_auc", "uno_c")
  }
  if (!metric %in% allowed) {
    sb_abort(
      paste0(
        "For a ", task_type, " task, `metric` must be one of: ",
        paste0("`", allowed, "`", collapse = ", "), "."
      ),
      "stormbreaker_bad_external_meta_metric"
    )
  }
  metric
}

sb_external_discrimination_meta_target <- function(
    bootstrap,
    cohort_id,
    metric) {
  if (!metric %in% c("time_auc", "uno_c")) {
    return(list(prediction_horizon = NULL, time_unit = NULL))
  }
  target <- bootstrap[["prediction_target"]]
  valid_target <- is.list(target) &&
    identical(target[["target_status"]], "complete") &&
    identical(target[["has_horizon"]], TRUE) &&
    is.numeric(target[["prediction_horizon"]]) &&
    length(target[["prediction_horizon"]]) == 1L &&
    !is.na(target[["prediction_horizon"]]) &&
    is.finite(target[["prediction_horizon"]]) &&
    target[["prediction_horizon"]] > 0 &&
    sb_external_scalar_character(target[["time_unit"]]) &&
    target[["time_unit"]] %in% c("day", "month", "year")
  if (!valid_target) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` requires a complete locked horizon and ",
        "time unit to meta-analyse `", metric, "`."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }
  list(
    prediction_horizon = as.numeric(target[["prediction_horizon"]]),
    time_unit = target[["time_unit"]]
  )
}

sb_external_discrimination_meta_member <- function(
    bootstrap,
    cohort_id,
    metric = NULL) {
  if (!inherits(bootstrap, "stormbreaker_external_bootstrap") ||
      !is.list(bootstrap)) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` is not a result from `sb_bootstrap_external()`."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  task_type <- bootstrap[["task_type"]]
  pipeline_id <- bootstrap[["pipeline_id"]]
  learner <- bootstrap[["learner"]]
  n <- bootstrap[["n"]]
  reps <- bootstrap[["reps"]]
  valid_metadata <- sb_external_scalar_character(task_type) &&
    task_type %in% c("binary", "survival") &&
    sb_external_scalar_character(pipeline_id) &&
    sb_external_scalar_character(learner) &&
    is.numeric(n) && length(n) == 1L && !is.na(n) && is.finite(n) &&
    n >= 1 && n <= .Machine$integer.max && n == as.integer(n) &&
    is.numeric(reps) && length(reps) == 1L && !is.na(reps) &&
    is.finite(reps) && reps >= 2 && reps <= .Machine$integer.max &&
    reps == as.integer(reps)
  if (!valid_metadata) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has incomplete external-bootstrap metadata."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  metric <- sb_external_discrimination_meta_metric(task_type, metric)
  metric_target <- sb_external_discrimination_meta_target(
    bootstrap,
    cohort_id,
    metric
  )
  point <- bootstrap[["point_estimates"]]
  replicates <- bootstrap[["replicates"]]
  valid_tables <- is.data.frame(point) &&
    all(c("metric", "estimate", "status") %in% names(point)) &&
    is.character(point$metric) && is.numeric(point$estimate) &&
    is.character(point$status) &&
    is.data.frame(replicates) &&
    all(c("replicate", "metric", "estimate", "status") %in%
      names(replicates)) &&
    is.numeric(replicates$replicate) &&
    is.character(replicates$metric) && is.numeric(replicates$estimate) &&
    is.character(replicates$status)
  if (!valid_tables) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` does not retain usable bootstrap metric tables."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  point_rows <- point[!is.na(point$metric) & point$metric == metric, , drop = FALSE]
  replicate_rows <- replicates[
    !is.na(replicates$metric) & replicates$metric == metric,
    ,
    drop = FALSE
  ]
  expected_replicates <- seq_len(as.integer(reps))
  replicate_ids_complete <- nrow(replicate_rows) == as.integer(reps) &&
    !anyNA(replicate_rows$replicate) &&
    all(is.finite(replicate_rows$replicate)) &&
    all(replicate_rows$replicate >= 1 &
      replicate_rows$replicate <= as.integer(reps)) &&
    all(replicate_rows$replicate == floor(replicate_rows$replicate)) &&
    identical(
      sort(as.integer(replicate_rows$replicate)),
      expected_replicates
    )
  point_available <- nrow(point_rows) == 1L &&
    identical(point_rows$status[[1L]], "success") &&
    length(point_rows$estimate) == 1L &&
    !is.na(point_rows$estimate[[1L]]) &&
    is.finite(point_rows$estimate[[1L]]) &&
    point_rows$estimate[[1L]] > 0 && point_rows$estimate[[1L]] < 1
  replicates_complete <- replicate_ids_complete &&
    !anyNA(replicate_rows$status) &&
    all(replicate_rows$status == "success") &&
    !anyNA(replicate_rows$estimate) &&
    all(is.finite(replicate_rows$estimate)) &&
    all(replicate_rows$estimate > 0 & replicate_rows$estimate < 1)
  if (!point_available || !replicates_complete) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` requires one finite successful `", metric,
        "` point estimate and all ", as.integer(reps),
        " successful bootstrap estimates, numbered exactly once from 1 ",
        "through `reps`, and strictly between 0 and 1. No failed replicate ",
        "is dropped and no boundary estimate is clipped."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  point_estimate <- point_rows$estimate[[1L]]
  logit_replicates <- stats::qlogis(replicate_rows$estimate)
  logit_se <- stats::sd(logit_replicates)
  if (!is.finite(logit_se) || logit_se <= 0) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has zero or non-finite bootstrap uncertainty for `", metric,
        "` on the logit scale."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  list(
    task_type = task_type,
    pipeline_id = pipeline_id,
    learner = learner,
    metric = metric,
    n = as.integer(n),
    reps = as.integer(reps),
    point_estimate = as.numeric(point_estimate),
    logit_estimate = stats::qlogis(point_estimate),
    logit_se = as.numeric(logit_se),
    prediction_horizon = metric_target$prediction_horizon,
    time_unit = metric_target$time_unit
  )
}

#' @rdname sb_meta_external_discrimination
#' @export
sb_meta_external_discrimination <- function(
    bootstraps,
    confidence_level = 0.95,
    metric = NULL) {
  valid_confidence <- is.numeric(confidence_level) &&
    length(confidence_level) == 1L && !is.na(confidence_level) &&
    is.finite(confidence_level) && confidence_level > 0 &&
    confidence_level < 1
  if (!valid_confidence) {
    sb_abort(
      "`confidence_level` must be one finite number strictly between 0 and 1.",
      "stormbreaker_bad_confidence_level"
    )
  }

  valid_input <- is.list(bootstraps) &&
    !inherits(bootstraps, "stormbreaker_external_bootstrap")
  if (!valid_input) {
    sb_abort(
      paste0(
        "`bootstraps` must be a named list of results from ",
        "`sb_bootstrap_external()`."
      ),
      "stormbreaker_bad_external_meta"
    )
  }
  if (length(bootstraps) < 2L) {
    sb_abort(
      "At least two independent external cohorts are required.",
      "stormbreaker_insufficient_external_cohorts"
    )
  }

  cohort_ids <- names(bootstraps)
  valid_ids <- !is.null(cohort_ids) &&
    length(cohort_ids) == length(bootstraps) &&
    !anyNA(cohort_ids) && all(nzchar(trimws(cohort_ids))) &&
    anyDuplicated(cohort_ids) == 0L
  if (!valid_ids) {
    sb_abort(
      "Outer list names must be unique, non-blank cohort IDs.",
      "stormbreaker_bad_external_meta"
    )
  }

  members <- lapply(seq_along(bootstraps), function(index) {
    sb_external_discrimination_meta_member(
      bootstraps[[index]],
      cohort_ids[[index]],
      metric
    )
  })
  task_types <- vapply(members, `[[`, character(1L), "task_type")
  if (length(unique(task_types)) != 1L) {
    sb_abort(
      "All cohorts must use the same task type.",
      "stormbreaker_external_meta_task_mismatch"
    )
  }
  pipeline_ids <- vapply(members, `[[`, character(1L), "pipeline_id")
  learners <- vapply(members, `[[`, character(1L), "learner")
  if (length(unique(pipeline_ids)) != 1L || length(unique(learners)) != 1L) {
    sb_abort(
      paste0(
        "All cohorts must retain the same pipeline ID and learner label. ",
        "These labels do not prove identical fitted state or formula."
      ),
      "stormbreaker_external_meta_model_mismatch"
    )
  }

  metric <- members[[1L]]$metric
  has_horizon <- metric %in% c("time_auc", "uno_c")
  metric_target <- list(
    has_horizon = has_horizon,
    prediction_horizon = NULL,
    time_unit = NULL
  )
  if (has_horizon) {
    horizons <- vapply(members, `[[`, numeric(1L), "prediction_horizon")
    time_units <- vapply(members, `[[`, character(1L), "time_unit")
    same_horizon <- all(horizons == horizons[[1L]])
    same_unit <- all(time_units == time_units[[1L]])
    if (!same_horizon || !same_unit) {
      sb_abort(
        paste0(
          "All cohorts must use the same time unit and prediction horizon ",
          "to meta-analyse `", metric, "`."
        ),
        "stormbreaker_external_meta_target_mismatch"
      )
    }
    metric_target$prediction_horizon <- horizons[[1L]]
    metric_target$time_unit <- time_units[[1L]]
  }
  estimand <- switch(
    metric,
    roc_auc = paste0(
      "Back-transformed random-effects mean of logit ROC AUC for one ",
      "locked model across an assumed distribution of independent ",
      "exchangeable cohorts"
    ),
    harrell_c = paste0(
      "Back-transformed random-effects mean of logit horizon-free Harrell ",
      "concordance for one locked model across an assumed distribution of ",
      "independent exchangeable cohorts"
    ),
    time_auc = paste0(
      "Back-transformed random-effects mean of logit cumulative/dynamic ",
      "time-dependent AUC at ", format(metric_target$prediction_horizon),
      " ", metric_target$time_unit, " across an assumed distribution of ",
      "independent exchangeable cohorts"
    ),
    uno_c = paste0(
      "Back-transformed random-effects mean of logit truncated Uno ",
      "concordance at ",
      format(metric_target$prediction_horizon), " ", metric_target$time_unit,
      " across an assumed distribution of independent exchangeable cohorts"
    )
  )

  if (!requireNamespace("metafor", quietly = TRUE)) {
    sb_abort(
      paste0(
        "Package `metafor` is required for external discrimination ",
        "meta-analysis. Install it with `install.packages(\"metafor\")`."
      ),
      "stormbreaker_missing_optional_dependency"
    )
  }

  cohorts <- data.frame(
    cohort_id = cohort_ids,
    n_patients = vapply(members, `[[`, integer(1L), "n"),
    point_estimate = vapply(members, `[[`, numeric(1L), "point_estimate"),
    logit_estimate = vapply(members, `[[`, numeric(1L), "logit_estimate"),
    logit_se = vapply(members, `[[`, numeric(1L), "logit_se"),
    n_bootstrap = vapply(members, `[[`, integer(1L), "reps"),
    stringsAsFactors = FALSE
  )
  fit <- metafor::rma.uni(
    yi = cohorts$logit_estimate,
    sei = cohorts$logit_se,
    method = "REML",
    test = "knha",
    level = 100 * confidence_level,
    slab = cohort_ids
  )
  prediction <- stats::predict(fit)
  k <- as.integer(nrow(cohorts))
  summary <- data.frame(
    metric = metric,
    k = k,
    estimate = stats::plogis(as.numeric(fit$b)),
    lower = stats::plogis(as.numeric(fit$ci.lb)),
    upper = stats::plogis(as.numeric(fit$ci.ub)),
    prediction_lower = stats::plogis(as.numeric(prediction$pi.lb)),
    prediction_upper = stats::plogis(as.numeric(prediction$pi.ub)),
    tau2_logit = as.numeric(fit$tau2),
    tau_logit = sqrt(as.numeric(fit$tau2)),
    i2_percent = as.numeric(fit$I2),
    stringsAsFactors = FALSE
  )
  caution <- if (k < 10L) {
    paste0(
      "Only ", k, " cohorts were supplied; heterogeneity and the ",
      "between-cohort prediction interval are especially uncertain."
    )
  } else {
    paste0(
      "The prediction interval assumes independent, exchangeable cohorts and ",
      "a correctly specified random-effects model."
    )
  }

  structure(
    list(
      status = "success",
      task_type = task_types[[1L]],
      pipeline_id = pipeline_ids[[1L]],
      learner = learners[[1L]],
      metric = metric,
      metric_target = metric_target,
      confidence_level = as.numeric(confidence_level),
      estimand = estimand,
      method = list(
        model = "random_effects",
        effect_scale = "logit_discrimination_probability",
        within_cohort_standard_error =
          "sample_sd_of_logit_patient_bootstrap_replicates",
        tau2_estimator = "REML",
        inference = "Knapp-Hartung",
        back_transformation = "inverse_logit"
      ),
      provenance = sb_external_meta_provenance(),
      cohorts = cohorts,
      summary = summary,
      prediction_interval_caution = caution,
      bootstraps = bootstraps
    ),
    class = "stormbreaker_external_discrimination_meta"
  )
}

#' @export
print.stormbreaker_external_discrimination_meta <- function(x, ...) {
  confidence <- 100 * x$confidence_level
  summary <- x$summary[1L, , drop = FALSE]
  cat("<stormbreaker_external_discrimination_meta:", x$status, ">\n")
  cat("fixed-model external discrimination meta-analysis\n")
  cat(
    "Metric:", x$metric,
    "| Cohorts:", summary$k,
    "| Task:", x$task_type, "\n"
  )
  if (isTRUE(x$metric_target$has_horizon)) {
    cat(
      "Target:", format(x$metric_target$prediction_horizon),
      x$metric_target$time_unit, "\n"
    )
  }
  cat(
    "Random effects: REML on logit discrimination | ",
    "Knapp-Hartung inference\n",
    sep = ""
  )
  cat(
    "Pooled:", format(summary$estimate),
    paste0("| ", format(confidence), "% CI: ["),
    format(summary$lower), ",", format(summary$upper), "]\n"
  )
  cat(
    paste0(format(confidence), "% prediction interval: ["),
    format(summary$prediction_lower), ",",
    format(summary$prediction_upper), "]\n"
  )
  cat(
    "Tau (logit discrimination):", format(summary$tau_logit),
    "| I-squared:", paste0(format(summary$i2_percent), "%\n")
  )
  cat(x$prediction_interval_caution, "\n")
  cat(
    paste0(
      "Cohort independence and common locked-model identity are caller ",
      "responsibilities; no model ranking or selection is performed.\n"
    )
  )
  invisible(x)
}
