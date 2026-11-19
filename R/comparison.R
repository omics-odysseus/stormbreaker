#' Cross-validation and cohort comparisons
#' @name sb_compare_cv
#' @description Evaluate pipelines with shared folds, repeated resampling or held-out cohorts.
#' @rawRd \arguments{
#' \item{task}{
#' \strong{For \code{sb_compare_cv()}, \code{sb_compare_repeated_cv()}:} A task returned by \code{\link[=sb_surv_task]{sb_surv_task()}} or \code{\link[=sb_binary_task]{sb_binary_task()}}.
#' 
#' \strong{For \code{sb_run_cv()}:} A task returned by \code{\link[=sb_surv_task]{sb_surv_task()}} or \code{\link[=sb_binary_task]{sb_binary_task()}}. A
#' survival task may declare one study-defined \code{prediction_horizon} and
#' clinically meaningful event-risk \code{decision_thresholds}; a binary task may
#' likewise declare action thresholds. The package never chooses either from
#' observed data or model performance.
#' A separate \code{group_col} keeps repeated observations together in generated
#' folds and glmnet learner tuning; supplied folds must respect these groups.
#' Prediction rows and metric weights still refer to observations. Other
#' engines' internal resampling is not configured by this declaration.
#' }
#' \item{pipelines}{
#' \strong{For \code{sb_compare_cv()}:} A non-empty list of distinct pipeline specifications
#' returned by \code{\link[=sb_pipeline]{sb_pipeline()}}. Input order is preserved.
#' 
#' \strong{For \code{sb_compare_repeated_cv()}:} Either a non-empty list of distinct pipeline specifications
#' returned by \code{\link[=sb_pipeline]{sb_pipeline()}} or one stack returned by \code{\link[=sb_stack]{sb_stack()}}. Input
#' and stack-member order are preserved within every repetition. A stack
#' repetition retains its member-plus-stack comparison; it does not select
#' among those rows.
#' 
#' \strong{For \code{sb_compare_cohorts()}:} Explicit \code{\link[=sb_pipeline]{sb_pipeline()}} specifications. If supplied,
#' leave \code{models} NULL.
#' 
#' \strong{For \code{sb_run_internal_external_cv()}:} Either a non-empty list of distinct, task-matched pipeline
#' specifications returned by \code{\link[=sb_pipeline]{sb_pipeline()}} or one task-matched stack from
#' \code{\link[=sb_stack]{sb_stack()}}. Specifications are fixed before the cycles; the function
#' does not tune across held-out cohorts. A stack is evaluated as one modeling
#' strategy; its members are not emitted as standalone IECV candidates.
#' 
#' \strong{For \code{bi_run()}:} Optional uniquely named list from \code{\link[=bi_pipeline]{bi_pipeline()}} or
#' \code{\link[=surv_pipeline]{surv_pipeline()}}. Mix single learners and ordered selector chains in one
#' comparison with shared folds. Names are used by \code{\link[=sb_complete_run]{sb_complete_run()}}.
#' Do not also supply models or selector; stages is ignored for this list.
#' }
#' \item{folds}{
#' \strong{For \code{sb_compare_cv()}:} Either one positive whole-number shared outer-fold count or a
#' multi-element numeric vector named by the complete unique task patient IDs.
#' Supplied values must be consecutive positive whole-number labels
#' \verb{1, ..., K}, with \code{K >= 2}, and are preserved after exact ID alignment.
#' 
#' \strong{For \code{sb_run_cv()}:} Either one positive whole-number cross-validation fold count,
#' for which each event or outcome stratum must contain at least that many
#' patients, or a multi-element numeric vector named by the complete unique
#' task patient IDs. Supplied values must be finite positive whole numbers
#' with labels exactly \verb{1, ..., K} for \code{K >= 2}. Vector order is ignored and
#' the aligned assignment is preserved exactly. A length-one value, even if
#' named, retains scalar-count behavior.
#' 
#' \strong{For \code{sb_compare_repeated_cv()}:} Either one positive whole-number outer-fold count or a non-empty
#' list of explicit patient-ID fold assignments. Each explicit assignment
#' follows the \code{folds} vector contract of \code{\link[=sb_compare_cv]{sb_compare_cv()}}.
#' 
#' \strong{For \code{sb_compare_cohorts()}:} Shared development CV settings.
#' 
#' \strong{For \code{bi_run()}:} Outer resampling controls.
#' }
#' \item{seed}{
#' \strong{For \code{sb_compare_cv()}, \code{sb_run_cv()}:} Positive integer used to generate scalar-count outer folds. It
#' does not choose a supplied patient-ID assignment, whose effective outer
#' seed is recorded as \code{NA}.
#' 
#' \strong{For \code{sb_compare_cohorts()}:} Shared development CV settings.
#' 
#' \strong{For \code{bi_run()}:} Outer resampling controls.
#' }
#' \item{reuse_selectors}{
#' Reuse successful, warning-free fits of the seeded
#' Boruta, varSelRF, survival Lasso, CoxBoost and RSF minimal-depth selectors,
#' and deterministic StepCox, within this comparison when their complete
#' analysis inputs and settings match exactly. Other selectors are unchanged.
#' No fit is shared across different outer or inner analysis inputs.
#' }
#' \item{progress}{
#' \strong{For \code{sb_compare_cv()}:} Show model and outer-fold progress messages. Defaults to
#' \code{interactive()}; set TRUE to include progress in script logs.
#' 
#' \strong{For \code{sb_run_cv()}:} Show outer-fold start/completion messages; defaults to
#' \code{interactive()}. Set TRUE for script logs.
#' 
#' \strong{For \code{sb_compare_cohorts()}:} Show model, fold, final-fit and evaluation progress.
#' Defaults to \code{interactive()}; use TRUE for script logs.
#' }
#' \item{pipeline}{
#' A matching pipeline returned by \code{\link[=sb_pipeline]{sb_pipeline()}}. \code{NULL}
#' selects the minimum default for the supplied task.
#' }
#' \item{seeds}{
#' For a scalar fold count, one or more positive whole-number seeds,
#' one per repetition. Use \code{NULL} with a list of explicit fold assignments.
#' }
#' \item{workers}{
#' Either \code{1L} (the unchanged serial default) or \code{2L}. Two local
#' processes execute independent repetitions, preserving within-repetition
#' sharing and input order. Engines remain single-threaded; memory use increases.
#' This does not parallelize pipelines within one repetition. Worker warnings
#' are relayed in repetition order and failures are not replaced or retried.
#' Native BlackBoost-Cox, RSF and AdaBoost learners consume a sequential caller
#' RNG stream and currently require \code{workers = 1}; they are not silently reseeded.
#' }
#' \item{development}{
#' A task from \code{\link[=sb_binary_task]{sb_binary_task()}} or \code{\link[=sb_surv_task]{sb_surv_task()}}.
#' }
#' \item{cohorts}{
#' \strong{For \code{sb_compare_cohorts()}:} A non-empty named list of same-task validation, test or
#' external cohort tasks. Patient/group IDs must be distinct across cohorts.
#' 
#' \strong{For \code{sb_run_internal_external_cv()}:} A named list of at least two \code{\link[=sb_binary_task]{sb_binary_task()}} objects or at
#' least two \code{\link[=sb_surv_task]{sb_surv_task()}} objects. Names are unique cohort IDs. All tasks
#' must have one task type and a common predictor-name set. Survival targets
#' and any decision-threshold declarations must also match exactly.
#' }
#' \item{roles}{
#' A named character vector matching \code{cohorts}, containing
#' \code{"validation"}, \code{"test"} or \code{"external"}. NULL labels all as external.
#' }
#' \item{models}{
#' \strong{For \code{sb_compare_cohorts()}:} Model IDs for \code{\link[=sb_model_pipelines]{sb_model_pipelines()}}; NULL uses its full menu
#' when \code{pipelines} is also NULL.
#' 
#' \strong{For \code{bi_run()}:} Method names from \code{\link[=sb_models]{sb_models()}}. NULL requests the whole task
#' menu. For binary tasks this is one shared combination pool: prediction
#' models and feature methods are not chosen from separate menus.
#' }
#' \item{inner_folds}{
#' Inner CV settings for glmnet menu models.
#' }
#' \item{inner_seed}{
#' Inner CV settings for glmnet menu models.
#' }
#' \item{development_name}{
#' Display name for development CV.
#' }
#' \item{view}{
#' Open the results page after computation. Defaults to FALSE;
#' set TRUE to explicitly request a browser.
#' }
#' \item{report_file}{
#' Optional HTML destination. Saves the page even when
#' \code{view = FALSE}. A new path is required.
#' }
#' \item{data}{
#' Data frame or an existing task. Repeated patients require a task
#' created with its group column; automatic row IDs do not infer patients.
#' }
#' \item{outcome}{
#' Outcome column (binary 0/1, or multiclass labels).
#' }
#' \item{predictors}{
#' Feature column names. NULL uses all non-role columns.
#' }
#' \item{id}{
#' ID column. NULL generates row IDs for independent rows.
#' }
#' \item{stages}{
#' Defaults to 1 when models are specified. Otherwise binary
#' defaults to c(1,2), and survival retains c(1,2,3). Binary 1 fits each
#' predictive method; 2 enumerates every supported extractor/model pair
#' within the chosen pool; explicit 3 enumerates two ordered extractors
#' followed by a model, with repetitions allowed. It is not just presets.
#' Large pools can be expensive; three-stage expansion is never implicit
#' for binary tasks. Survival retains its fixed three-stage presets.
#' }
#' \item{selector}{
#' Legacy explicit selector specification, overriding stages.
#' For new binary routes use \code{\link[=bi_pipeline]{bi_pipeline()}} with ordered model names.
#' }
#' \item{plots}{
#' NULL (default) draws comparison overviews, but only performs a
#' final fit and detailed plotting when exactly one pipeline is requested.
#' TRUE explicitly fits and plots all candidates; FALSE skips all plots and
#' final fits. Use \code{\link[=sb_complete_run]{sb_complete_run()}} for a chosen binary/survival pipeline
#' after batch comparison. No browser opens; unavailable views stay visible.
#' }
#' \item{decision_thresholds}{
#' Optional risk probabilities for decision curves.
#' Set these when planning the analysis. NULL leaves decision curves unrequested.
#' For an existing task, its declared settings are used.
#' }
#' \item{time}{
#' Survival time and 0/1 event columns.
#' }
#' \item{event}{
#' Survival time and 0/1 event columns.
#' }
#' \item{time_unit}{
#' Unit of the supplied follow-up times, default month.
#' }
#' \item{prediction_horizon}{
#' Optional survival prediction time in time_unit.
#' Enables horizon-specific metrics and curves. NULL requests risk ranking only.
#' For an existing task, its declared horizon is used.
#' }
#' \item{expr}{
#' An expression containing repeated comparisons, CV runs or final
#' fits. Do not change package implementations during the expression.
#' }
#' \item{enabled}{
#' Logical; FALSE executes the expression without reuse, even
#' inside an enabled outer context.
#' }
#' \item{diagnostics}{
#' Logical; TRUE returns a list with \code{result} and \code{reuse}
#' (selector hits, misses and entries, plus separately named
#' \code{final_learner_hits}, \code{final_learner_misses} and \code{final_learner_entries});
#' analogous \verb{cv_learner_*} fields describe CV learner reuse;
#' \code{learner_bytes}, \code{learner_byte_limit}, \code{learner_evictions} and
#' \code{learner_oversized} report the bounded learner collection;
#' \code{preparation_hits} and \code{preparation_misses} describe native tuning reuse;
#' FALSE returns the expression value unchanged.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_compare_cv()}}{
#' \code{sb_compare_cv()} runs every supplied pipeline with the same survival or
#' binary task and outer-fold specification. A scalar count creates one shared
#' deterministic event- or outcome-stratified assignment. A multi-element
#' numeric vector named by every task patient ID supplies one exact assignment
#' to every pipeline; names, not vector position, bind patients to folds. The
#' package does not repair or rebalance a supplied plan, derive it from model
#' performance, or certify its design quality. The returned table copies each
#' run's fold-stratified Harrell C, or
#' its IPCW Brier score, fold-local Kaplan-Meier reference Brier, Brier skill
#' score, cumulative/dynamic time-dependent AUC, truncated Uno IPCW
#' concordance, pooled observed/expected mean calibration, calibration slope,
#' and moderate-calibration ICI/E50/E90 when the survival task declared a
#' prediction horizon, or binary
#' ROC AUC,
#' fold-stratified average precision, patient-level OOF Brier
#' score and log loss, their analysis-prevalence references, and their skill
#' scores, descriptive OOF calibration intercept and slope, a separate saved
#' binary moderate-calibration curve's ICI/E50/E90 summaries, an optional saved
#' task-appropriate decision curve, the saved fold-local selection-stability
#' status and mean pairwise Jaccard, and the separate glmnet nonzero active-set
#' stability audit. Pipelines without a learned selector are marked not
#' applicable for selector stability; Ridge and non-glmnet learners are marked
#' not applicable for active-set stability. The comparison does
#' not rank predictions from different folds against each other or recompute a
#' different metric. Weak or moderate calibration equations are not applied to
#' predictions. The binary moderate-calibration curve is not refitted here and
#' ICI/E50/E90 have no universal acceptable cutoff. Decision curves are copied
#' as complete per-pipeline objects; they are not flattened, integrated,
#' smoothed, ranked, or used to select a threshold.
#' Survival censoring distributions, reference survival, probabilities, losses,
#' within-fold time-AUC and Uno C pair quantities, both saved weak-calibration
#' summaries, the moderate-calibration curve, and its ICI/E50/E90 summaries are
#' likewise taken from the saved CV result and are never re-estimated by this
#' comparison.
#' The moderate curve is descriptive and is not applied to predictions. No
#' survival calibration intercept or universal calibration cutoff is reported.
#' 
#' \strong{Returns:} A \code{stormbreaker_cv_comparison} containing the task type, a
#' descriptive \code{comparison} data frame, the complete named \code{results} from
#' \code{\link[=sb_run_cv]{sb_run_cv()}}, the shared fold count, effective seed, source, aligned
#' patient-to-fold assignment, and, when declared, a named list of copied
#' \code{decision_curves}.
#' 
#' \strong{Details:} This is a descriptive comparison. It does not rank pipelines, select a
#' winner, estimate uncertainty, or turn development CV into external
#' validation. Both Jaccard summaries are descriptive: outer analysis folds
#' overlap, and no confidence interval or chance correction is inferred. A
#' not-evaluable glmnet active-set audit does not erase otherwise complete OOF
#' predictions or metrics. A
#' failed pipeline remains in the table with its complete audited CV result and
#' no replacement learner; later supplied pipelines are still run in their
#' original order. Task, fold, and pipeline-contract errors are validated
#' before this collection behavior applies.
#' Fixed LightGBM, CatBoost, and ranger pipelines are compared as their declared
#' configurations; this function does not reinterpret them as tuned algorithms.
#' }
#' 
#' \subsection{\code{sb_run_cv()}}{
#' \code{sb_run_cv()} dispatches from the supplied task to a survival or binary
#' pipeline. A scalar \code{folds} value creates deterministic folds stratified by
#' event or outcome. A multi-element numeric vector named by every task patient
#' ID instead supplies the exact outer-fold assignment; names, not vector
#' position, bind patients to folds. The package preserves a structurally valid
#' supplied plan without shuffling, stratifying, repairing, or rebalancing it.
#' Such acceptance does not certify design quality, independence,
#' representativeness, or transportability. A real study should declare or
#' import its plan independently of model performance.
#' 
#' \strong{Returns:} A \code{stormbreaker_cv_result} with:
#' \itemize{
#' \item \code{predictions}: patient-level outcomes, OOF risk or class-1 probability,
#' fold identity, and status; a survival task with a declared horizon also
#' includes OOF survival probability and its status;
#' \item \code{metrics}: fold-stratified Harrell concordance, or binary ROC AUC,
#' non-interpolated average precision, and patient-level OOF Brier score and
#' log loss with outer-analysis prevalence references and skill scores, plus
#' descriptive OOF weak-calibration intercept and slope plus a separate pooled
#' binary moderate-calibration curve and ICI/E50/E90, plus an optional
#' patient-count decision curve with threshold operating characteristics at
#' declared thresholds; per-fold audits and failures are retained, and no
#' confidence interval is inferred from fold variability; a survival task
#' with a declared horizon additionally reports
#' analysis-fold IPCW Brier score, its censoring curves, a fold-local event
#' Kaplan-Meier reference, the aggregate Brier skill score, and a
#' within-fold cumulative/dynamic time-dependent AUC and truncated Uno IPCW
#' concordance, plus a pooled
#' three-knot Cox moderate-calibration curve and ICI/E50/E90, a pooled OOF
#' observed/expected mean-calibration ratio and calibration slope, and, when
#' declared, a pooled OOF survival decision curve at the locked horizon;
#' \item \code{runs}: one audit row per fold;
#' \item \code{selection_stability}: feature selection counts and frequencies over
#' successful selector folds, plus all unique outer-analysis fold-pair
#' Jaccard audits; pipelines without a learned selector are not applicable;
#' \item \code{glmnet_active_set_stability}: for glmnet with \code{alpha > 0}, exact nonzero
#' coefficient counts and frequencies plus all outer-fold active-set Jaccard
#' audits; Ridge and non-glmnet learners are not applicable, and empty-union
#' fold pairs remain not evaluable;
#' \item \code{selections}: for supervised selectors, selected features and complete
#' per-feature selection audit tables for every fold; a fixed selector retains its
#' unchanged feature block without a score table;
#' \item \code{tuning}: for glmnet learners, the selected \code{lambda}, predeclared lambda
#' rule and standardization setting, lambda alignment, inner patient-to-fold
#' assignment, error curve, compact fold-local
#' selection audits when requested, and any tuning failure for every outer
#' fold;
#' \item \code{coefficients}: every fitted linear-survival coefficient, including
#' zero-valued glmnet and SuperPC coefficients, every binary logistic
#' coefficient including its intercept, or the cancerclass inner-OOF
#' calibration intercept and score slope; Gaussian Naive Bayes and the fixed
#' tree learners mark this audit as not applicable;
#' \item \code{baseline_hazards}: when a survival horizon is declared, each Cox fold's
#' \code{H0(t)}, \code{S0(t)}, engine-specific ties, source, failure, and any
#' calibration coefficient multiplying the saved risk before applying
#' \code{H0(t)}; ranger split-rule forest folds are transparently marked not
#' applicable because their probabilities come directly from the forest;
#' \item \code{distribution_parameters}: for Gaussian Naive Bayes, every outer-fold
#' class prior, smoothing epsilon, mean, variance, and standard deviation;
#' \item \code{tree_models}: for fixed LightGBM, fixed CatBoost (binary, Cox, or AFT), fixed
#' GBM-Cox, fixed ranger binary or survival split-rule forests, fixed AdaBoost.M1, or
#' fixed LogitBoost, every outer-fold serialized model, declared setting,
#' fitted feature name, and all-feature importance audit; binary boosting also
#' records requested and fitted iteration counts with split-count importance;
#' \item \code{folds}: exact analysis and assessment patient IDs;
#' \item \code{provenance}: the outer-fold source, effective count and seed, task-order
#' patient-to-fold assignment, and R and dependency versions;
#' \item \code{complete_oof}: whether every patient has one finite OOF risk or binary
#' probability;
#' \item \code{complete_survival_probability_oof}: when requested, whether every patient
#' has one finite OOF survival probability at the declared horizon.
#' }
#' 
#' \strong{Details:} For each fold, the function fits the requested pipeline on analysis patients
#' and predicts assessment patients exactly once. A glmnet Cox or logistic
#' pipeline selects \code{lambda} by stratified inner CV using only that outer
#' fold's analysis patients, then refits on all of them. With a matching Top-K
#' selector, selection is repeated separately in every inner analysis fold;
#' it is never learned once from the corresponding inner assessment outcomes.
#' Ahaz Lasso and stepwise-SCAD instead select one lambda by BIC on the current
#' outer analysis fold only; they do not inspect its assessment patients.
#' A \code{fixed_features} selector instead applies one predeclared column block in
#' the same order in every fold. It never reads an outcome. A glmnet learner
#' tunes lambda within that already fixed block without treating it as
#' fold-local supervised selection.
#' Ordinary Cox, glmnet Cox, fixed CoxBoost, plsRcox, SuperPC, and the calibrated
#' linear survival-SVM return the zero-reference linear predictor \code{x_beta};
#' additive and RBF survival-SVM return their negative native kernel scores times
#' their analysis-fold Cox calibration coefficients. Fixed CatBoost Cox returns
#' its native raw margin. Fixed CatBoost SurvivalAft and XGBoost AFT return the
#' negative native log-time margin; CatBoost uses a fixed Normal distribution
#' for direct horizon survival. Fixed GBM-Cox and XGBoost-Cox return their
#' native link-scale log hazards, while mboost AFT returns the negative native
#' distributional log-time link. Ahaz learners return their native additive-
#' hazard \code{x_beta}; their
#' declared-horizon probability is separately labelled analysis-fold
#' Cox/Breslow calibration. Bujar returns negative predicted log time.
#' Discrete-time survival BART summarizes its
#' posterior mean curve as negative stepwise restricted mean survival over the
#' analysis-fold time grid; its horizon survival remains a direct posterior
#' probability rather than a proportional-hazards calibration. The RLT
#' survival learner reconstructs native discrete hazard increments into a
#' survival curve and reports negative stepwise restricted mean survival over
#' its native event grid. Full Aalen
#' regression reports cumulative hazard at the last analysis-fold event time;
#' proportional-odds survival regression reports \code{x_beta}. Every binary
#' learner returns the
#' probability of class \code{1}. Gaussian Naive Bayes estimates empirical priors
#' and smoothed
#' class-feature Gaussian parameters separately in every outer analysis fold.
#' Cancerclass fits its native all-feature Welch/correlation score in each
#' analysis fold and learns its probability mapping only from five-fold inner
#' out-of-fold scores; it is not admitted after a learned upstream selector.
#' The four LVQ learners instead expose one declared logistic mapping of the
#' nearest class-specific squared-prototype distance margin. It preserves the
#' native nearest-prototype classification boundary but is not described as a
#' calibrated native posterior. Scaling, codebook initialization, and prototype
#' updates are repeated separately inside each analysis fold.
#' When glmnet predictor standardization is enabled in the pipeline, glmnet
#' learns it separately inside each inner analysis fit and repeats it for the
#' outer-analysis refit; it is not a package-wide preprocessing step.
#' The fixed tree learners use only predeclared settings and supply no
#' assessment patients for early stopping or best-model selection. A binary
#' AUC Top-K selector, when requested, is scored and fitted separately within
#' every outer analysis fold. Harrell
#' concordance, ROC AUC, and non-interpolated average precision use only saved
#' OOF predictions within the same outer assessment fold. Average precision is
#' aggregated by fold positive count without ranking patients scored by
#' different models. Binary Brier score averages patient-level squared errors,
#' while non-clipped log loss averages patient-level negative log likelihood
#' from the saved class-1 OOF probabilities. Each reference probability is
#' learned from the corresponding outer analysis fold's class-1 prevalence;
#' Brier and log-loss skill compare model and reference loss sums. Descriptive
#' calibration-in-the-large and calibration slope are fitted only to the pooled
#' saved OOF probabilities and outcomes. They assess weak calibration and never
#' alter an OOF prediction; boundary probabilities are not clipped, and a
#' separated recalibration fit is not reported as a finite slope. A separate
#' binary moderate-calibration curve fits one three-knot restricted natural
#' cubic spline logistic model to \code{qlogis()} of the pooled OOF probabilities.
#' Its ICI, E50, and E90 summarize the empirical distribution of absolute
#' differences between predicted and smoothed observed class-1 probability.
#' Exact zero or one is not clipped and leaves only that finite-logit curve
#' unavailable. A post-fit bounded signed-margin linear program also leaves a
#' complete or quasi-completely separated curve unavailable even if logistic
#' IRLS reports numerical convergence. The curve is not applied as recalibration
#' and has no universal acceptable cutoff. If the binary task declares decision
#' thresholds, the saved OOF probabilities also enter a patient-count decision
#' curve using \code{probability >= threshold}. Net benefit, treat-all, and
#' treat-none are reported at the declared thresholds. The same TP/FP/TN/FN
#' counts give sensitivity, specificity, positive predictive value, and
#' negative predictive value. A zero denominator remains \code{NA}; predictive
#' values describe the pooled OOF prevalence. The package does not derive a
#' grid, smooth or integrate the curve, or select a threshold. Counts are
#' additive over patients and do not compare predictions from different folds.
#' For a learned selector, the already
#' saved outer-analysis feature sets are also summarized by per-feature
#' selection frequency and every unique
#' fold-pair Jaccard similarity. The mean Jaccard is withheld unless every
#' requested selector fold succeeded. Because outer analysis folds overlap,
#' this is descriptive stability without uncertainty or chance correction.
#' A glmnet learner with \code{alpha > 0} separately summarizes the exact nonzero
#' predictor coefficients at each outer fold's selected \code{lambda}. This active-
#' set audit is distinct from any explicit Top-K selector. Ridge (\code{alpha = 0})
#' and non-glmnet learners are not applicable, no coefficient cutoff is added,
#' and a fold pair with two empty active sets is left not evaluable rather than
#' called perfectly stable. This descriptive audit does not invalidate
#' otherwise complete OOF predictions.
#' Every fitted Cox-like linear survival learner and each ahaz learner also
#' retain their outer-fold coefficients so the declared linear predictor can be
#' reconstructed as \code{x_beta}. Additive and RBF survival-SVM outer-fold results
#' do not retain their kernel basis or native weights. A later \code{sb_fit_final()}
#' stores the scaled kernel basis and native weights required for outcome-blind
#' prediction, plus the actual RBF kernel parameter when applicable, without
#' outcomes or the fitted engine. Fixed CatBoost Cox and SurvivalAft, GBM-Cox,
#' and XGBoost survival learners instead retain their serialized models and
#' fixed prediction metadata.
#' Fixed ranger binary and
#' survival split-rule forests retain their serialized forest and descriptive impurity
#' importance. Fixed AdaBoost.M1 and LogitBoost retain their serialized native
#' models, requested and fitted iteration counts, and descriptive split-count
#' importance. Ranger survival risk is mortality, the row sum of the cumulative
#' hazard function. When the survival task declares one
#' \code{prediction_horizon}, ordinary Cox retains an Efron baseline. Glmnet Cox,
#' fixed CoxBoost, plsRcox, SuperPC, calibrated survival-SVM, CatBoost
#' Cox, GBM-Cox, XGBoost-Cox, and ahaz use a Breslow baseline matched to the
#' saved risk score; ahaz first fits its explicitly labelled nonnegative
#' analysis-fold Cox calibration coefficient. Ranger, CatBoost SurvivalAft,
#' XGBoost AFT, mboost AFT,
#' Bujar, SurvReg, randomForestSRC, partykit
#' Ctree/CForest, survival BART, Aalen, and proportional odds instead read
#' survival probability directly from
#' their fitted distributions, residual Kaplan-Meier curve, or coefficient
#' processes and do not invent a Cox baseline. Aalen uses its native
#' \code{exp(-H(t | x))}; an out-of-range native
#' estimate is reported as a probability failure rather than clipped or
#' recalibrated. Survival Rpart
#' uses an analysis-fold Breslow calibration matched to its relative-rate
#' score. The
#' corresponding IPCW Brier
#' score estimates each censoring distribution from the matching outer
#' analysis fold, orders tied events before censoring, and uses the
#' left-continuous censoring survival at event and prediction times. Censoring
#' exactly at the horizon is treated as known event-free at that instant. A
#' no-covariate event Kaplan-Meier reference is fitted in each outer analysis
#' fold and scored with the same assessment patients and IPCW weights. Model
#' and reference weighted losses are each divided by all assessment patients;
#' their aggregate loss ratio supplies a descriptive survival Brier skill
#' score.
#' 
#' The same locked survival probabilities also produce a cumulative/dynamic
#' time-dependent AUC. The marker is \code{1 - S(t)}: observed events at or before
#' the horizon are cases, patients followed strictly beyond it are controls,
#' and patients censored at or before it are excluded. A case is weighted by
#' the inverse marginal censoring survival immediately before its event time;
#' exact marker ties receive half credit. Case-control pairs are formed only
#' inside the matching outer assessment fold, never across predictions from
#' different fitted models. The headline estimate divides summed successful
#' fold pair credit by the corresponding summed pair denominator and is
#' withheld unless every requested fold is evaluable. This uses the marginal
#' non-informative censoring assumption. Time-dependent AUC measures
#' discrimination, not calibration, transportability, clinical utility,
#' significance, or a decision threshold. These development-set OOF quantities
#' are not external validation.
#' 
#' The saved risk scores separately produce Uno's truncated IPCW concordance at
#' the declared horizon. An index patient must have an observed event strictly
#' before the horizon; patients observed later and patients censored at the same
#' recorded time are comparators, while events tied at the same time are not.
#' Each eligible event receives 1 / G(T-)^2 from the matching outer-analysis
#' censoring curve, and exact risk ties receive half credit. Comparisons stay
#' within the matching outer assessment fold. The headline estimate is the
#' ratio of summed fold weighted pair credit to summed fold weighted pair
#' denominator, not a mean of fold estimates, and is withheld unless every fold
#' is evaluable. This statistic uses risk scores rather than survival
#' probabilities and does not compare scores from different outer models.
#' 
#' When a horizon-enabled survival task also declares event-risk thresholds,
#' the same saved OOF probabilities enter a survival decision curve using
#' \code{1 - survival_probability >= threshold}. At each threshold, the selected
#' fraction is multiplied by the Kaplan-Meier event and event-free
#' probabilities among selected OOF patients to estimate true- and false-
#' positive rates. Treat-all uses the Kaplan-Meier survival of all OOF patients
#' and treat-none is zero. A nonempty group whose follow-up ends before the
#' horizon is unavailable unless its Kaplan-Meier survival has already reached
#' the absorbing value zero. This published estimator assumes independent
#' non-informative censoring in the relevant risk groups and censoring
#' independent of model predictors. It neither uses the IPCW Brier/time-AUC
#' weights nor refits a model. The curve is a point estimate at predeclared
#' thresholds; it does not choose a threshold or prove clinical effectiveness.
#' 
#' Horizon-enabled survival results additionally fit one pooled moderate-
#' calibration curve to the saved OOF probabilities. A three-knot restricted
#' natural cubic spline of \code{log(-log(S(t)))} enters an Efron Cox model using the
#' complete right-censored follow-up; the calibrated survival probability is
#' then read at the locked horizon. ICI is the mean patient-level absolute
#' difference between predicted and calibrated survival probability, while E50
#' and E90 are its type-7 50th and 90th percentiles. The curve is descriptive,
#' assumes proportional hazards and independent non-informative censoring, and
#' has no universal acceptable cutoff. It is not applied as recalibration,
#' does not clip or change a prediction, and is distinct from mean or weak
#' calibration.
#' 
#' The same horizon predictions separately report weak calibration. Mean
#' calibration is the pooled Kaplan-Meier observed event probability at the
#' horizon divided by mean predicted event probability (O:E); its target is 1.
#' The calibration slope is the coefficient from an Efron Cox model of the
#' endpoint administratively censored at the horizon on
#' \code{log(-log(survival_probability))}; its target is 1. Strictly later outcomes
#' cannot change that slope. The two components retain independent status:
#' exact zero/one or constant predictions can leave the slope unavailable
#' without erasing the mean ratio, and unsupported Kaplan-Meier follow-up can
#' leave the mean ratio unavailable without erasing an otherwise estimable
#' slope. No probability is clipped or recalibrated, and no survival calibration
#' intercept, p value, or universal cutoff is reported. The pooled OOF values
#' describe predictions from different outer fits; they are not external
#' validation or refit-uncertainty estimates. Any calibration failure remains
#' optional and does not erase otherwise valid OOF evidence.
#' Failed fits remain failed in the run, selection, prediction, fitted-state,
#' and metric audits; no replacement learner is used.
#' }
#' 
#' \subsection{\code{sb_compare_repeated_cv()}}{
#' \code{sb_compare_repeated_cv()} repeats either \code{\link[=sb_compare_cv]{sb_compare_cv()}} for an ordered
#' pipeline list or \code{\link[=sb_run_stack_cv]{sb_run_stack_cv()}} for one predeclared stack, without
#' changing either function's fitting, prediction, metric, or failure behavior.
#' Within a pipeline comparison every pipeline uses the same outer folds; within
#' a stack repetition every member and the stack do so. Each patient's
#' out-of-fold predictions remain inside that repetition's complete result.
#' 
#' \strong{Returns:} A \code{stormbreaker_repeated_cv_comparison} containing a long
#' repeat-by-pipeline \code{comparison} table, every complete
#' single-run comparison or stack result in \code{repetitions}, and a long
#' \code{fold_assignments} table keyed by repetition and patient ID.
#' 
#' \strong{Details:} There are two explicit resampling modes. When \code{folds} is one fold count,
#' each value in \code{seeds} generates one event- or outcome-stratified assignment.
#' When \code{folds} is a list, each element must be one numeric vector named by all
#' task patient IDs and \code{seeds} must be \code{NULL}. Explicit assignments are passed
#' unchanged to \code{\link[=sb_compare_cv]{sb_compare_cv()}}, which aligns them by patient ID and neither
#' repairs nor rebalances them.
#' 
#' The long \code{comparison} table contains one row per repetition and pipeline.
#' It copies the corresponding single-comparison rows and does not pool patient
#' predictions across repetitions. Repeated splits show sensitivity to outer-
#' fold assignment; they are not independent cohorts or a confidence interval.
#' This function calculates no cross-repeat mean, standard deviation, interval,
#' ranking, winner, or final-model choice. Duplicate seeds or fold assignments
#' are preserved and remain visible rather than silently replaced.
#' }
#' 
#' \subsection{\code{sb_compare_cohorts()}}{
#' Run shared-fold development CV, fit each declared pipeline on development
#' patients, and apply those same fitted states to every supplied cohort.
#' The result opens as a model-by-cohort matrix in an interactive R session.
#' 
#' \strong{Returns:} A \code{stormbreaker_cohort_comparison} with development comparison,
#' final fits, cohort predictions/evaluations, named result sets, roles,
#' and an optional \code{report_file}.
#' 
#' \strong{Details:} Cohort outcomes are used only for evaluation, after all models
#' have been fitted on development data. No winner is chosen and validation
#' or test patients are not added to the fitting data. To assess a model
#' selected using validation results, keep a separate test cohort.
#' }
#' 
#' \subsection{\code{sb_run_internal_external_cv()}}{
#' \code{sb_run_internal_external_cv()} performs internal-external cross-validation
#' (IECV) over a named collection of compatible cohorts. In each cycle, one
#' complete cohort is excluded, either every predeclared pipeline or one
#' predeclared stack is fitted using only patients from the remaining cohorts,
#' and endpoint-free predictions are made for the excluded cohort. Only after
#' prediction is its outcome-only endpoint constructed and passed to
#' \code{\link[=sb_compare_external]{sb_compare_external()}}. Every cohort is excluded exactly once.
#' 
#' \strong{Returns:} A \code{stormbreaker_internal_external_cv} containing overall execution
#' \code{status}, the validation design, cohort and strategy summaries, one named
#' \code{fold} per held-out cohort with its development IDs, final fits,
#' endpoint-free predictions, and complete external comparison, plus a long
#' cohort-by-pipeline \code{results} table. Failed fits and predictions remain
#' visible and do not trigger fallback; later predeclared pipelines and
#' held-out-cohort cycles continue. Cohort, task, and strategy-contract
#' errors are validated before this collection behavior applies.
#' 
#' \strong{Details:} This function evaluates the supplied modeling strategy using
#' pooled development rows from the remaining cohorts. Pooling adds no intercepts
#' or strata by cohort, random effects, or batch correction.
#' Cohorts must therefore declare the same numeric predictor-name set, time
#' target, and decision thresholds. Predictor order is normalized to the first
#' cohort, but feature intersection, measurement harmonization, eligibility,
#' outcome-definition harmonization, and outcome-blind upstream preprocessing
#' remain the caller's responsibility.
#' 
#' A different fitted model is evaluated in every held-out cohort, and with
#' three or more cohorts their development sets overlap. Accordingly, the
#' cycle estimates are not independent evaluations of one locked model and
#' must not be passed to STORMBREAKER's external-performance meta-analysis
#' functions. No metric is pooled, no pipeline is ranked or selected, and no
#' final all-cohort model is fitted. If cycle results are inspected and used to
#' revise a pipeline, they are method-development evidence rather than final
#' independent external validation; that claim requires untouched data. Two
#' cohorts are executable, although the method literature has tentatively
#' recommended at least four studies for a useful heterogeneity appraisal.
#' 
#' Patient IDs must be globally unique across cohorts. This detects direct
#' patient overlap between a development set and its held-out cohort. If two
#' independent cohorts reuse local labels such as \code{"P001"}, namespace those
#' labels before creating the tasks rather than allowing the function to hide
#' the ambiguity by prefixing them automatically.
#' }
#' 
#' \subsection{\code{bi_run()}}{
#' Run binary models with short names and fixed algorithm defaults
#' 
#' \strong{Returns:} List with CV comparison, requested pipelines, unavailable combinations
#' and final-fit figures. The issues table offers readable explanations while
#' preserving original reasons. No winner is automatically selected.
#' The result retains development data for \code{\link[=sb_complete_run]{sb_complete_run()}}; treat saved
#' run objects as patient-level data, not anonymous reports.
#' Printing shows a short summary; all existing list fields remain available.
#' Optional plotting failures are recorded in issues without discarding fits.
#' }
#' 
#' \subsection{\code{sb_with_selector_reuse()}}{
#' This helper shares successful, warning-free Boruta binary/survival and
#' varSelRF, the high-dimensional omics selectors, ABESS, PAMR, binary model
#' selectors, binary AUC/survival Cox Top-K, survival Lasso/Elastic Net,
#' Ridge/PLS-Cox/survival-SVM weight
#' selection, CoxBoost, RSF minimal-depth, StepCox and survival boosting fits when the
#' complete predictor frame (including row and column identities), endpoint,
#' selector settings (including requested diagnostic history) and RNG kind match.
#' Eligible seeded selectors preserve the caller's random state;
#' deterministic StepCox does not use random draws.
#' Sequential MIME selectors use the same mechanism for their individual
#' stages, allowing identical prefixes to be shared with a single stage.
#' A non-first SVM-RFE stage also matches the original analysis predictors and
#' upstream selector specifications, retaining prefix refitting in its inner folds.
#' Other selectors execute normally.
#' No input or fitted state is persisted.
#' 
#' \strong{Returns:} The value of \code{expr}, or a result/diagnostics list.
#' 
#' \strong{Details:} Reuse does not preselect genes on all patients. Different outer-analysis,
#' inner-analysis and final-training inputs receive separate fits. Final fits
#' may also reuse fixed survival-SVM, logistic, Cox, Gaussian naive Bayes and
#' prepared glmnet fits, fixed StepGLM/GLMBoost/PLS-GLM and the classical binary
#' family (including trees, forests, discriminant, kernel and prototype models)
#' and mboost Cox/AFT/ranking or mgcv GAM fits after any selector.
#' Analysis-CV learners reuse their complete search and fit after any legally
#' registered selector; this does not create additional combinations. Models are reused only
#' when ordered predictors, patient IDs, endpoint and learner settings match.
#' Tuned-learner reuse includes the complete native search and selected fit,
#' not a search performed on other data.
#' Matching final inputs are retained within the expression, including when
#' another input for the same learner is fitted between two matching calls.
#' CV and final learner retention share a 256 MiB per-expression payload budget.
#' Stored models are serialized; entry sizes include matching inputs and metadata.
#' Least-recently-used entries are released when needed. An oversized model is
#' returned normally but not retained. Evicted inputs are fitted normally if
#' needed again. This bounds retained learner entries, not total R memory:
#' active fitting, returned results, selectors and preparation have other costs.
#' Each route still constructs its own selection, portable state, horizon
#' baseline, diagnostics and provenance. CV fits of fixed survival-SVM, logistic,
#' Cox, Gaussian naive Bayes, prepared glmnet, fixed StepGLM/GLMBoost/PLS-GLM
#' and the classical binary, mboost Cox/AFT/ranking and mgcv GAM families
#' use a separate in-memory collection of matching
#' analysis inputs; different analysis sets are fitted separately. Other CV learners
#' execute normally; tuning is never transferred to different analysis data.
#' Glmnet and SVM probability preparation also share identical native searches, matching learner
#' controls, patients, endpoints, groups and full upstream selection context.
#' Identical outer-selected columns alone never identify a learned-selector search.
#' Each consumer retains its own inner-selection method description.
#' SVM model matching additionally includes its fitted probability calibration.
#' Native BlackBoost-Cox, RSF and AdaBoost fits also match the caller's initial random state and
#' replay its random-state transition; different starting states are not shared.
#' Mutable mboost model environments are isolated by serialization, so changing
#' one consumer's fitted model cannot change another consumer's model.
#' SurvReg, ranger, flexsurv/spline, timereg, ahaz, Firth Cox, ncvreg and ABESS fits,
#' omics-specialist binary, mclust, mixOmics and HDDA fits also share identical
#' analysis inputs in CV and final fitting. Their native searches, when present,
#' remain confined to those exact analysis data.
#' CatBoost, fixed GBM/CoxBoost/PLS-Cox/SuperPC, MTLR, PAMR survival,
#' Akritas, bagged survival kNN, Bujar, RLT, WSRF, AdaBoost/LogitBoost and
#' rpart/RSF/ctree/cforest/AORSF survival fits are also eligible. Eligibility
#' alone does not guarantee a hit: inputs and required random states must match.
#' SLOPE, randomGLM and cancerclass may share the complete analysis-local fit,
#' including native searches or out-of-fold calibration. This does not add new
#' selector combinations or transfer internal calibration to other patients.
#' BART survival/AFT, PRE, penAFT, CV-tuned CoxBoost/PLS-Cox/SuperPC/GBM-Cox
#' also match their complete analysis inputs. Failed or warning-bearing fits
#' are not retained; every affected route preserves its own failure reporting.
#' XGBoost and LightGBM fits use native model serialization for independent
#' booster handles. LightGBM retains original fit importance values because its
#' text model rounds split gains; predictions and reported importance stay unchanged.
#' Errors are never cached. Ordinary selector fits with warnings are not cached.
#' Final learner fits with warnings are also executed again, including warnings
#' retained as fit metadata from unselected PLS-Cox components.
#' The context is discarded on normal
#' return or error; nested calls have independent contexts.
#' }
#' }
#' @rawRd \references{
#' McLernon DJ, Giardiello D, Van Calster B, et al. (2023).
#' Assessing performance and clinical usefulness in prediction models with
#' survival outcomes: practical guidance for Cox proportional hazards models.
#' \emph{Annals of Internal Medicine}, 176(1), 105-114.
#' \doi{10.7326/M22-0844}.
#' 
#' Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ (2011). On the C-statistics
#' for evaluating overall adequacy of risk prediction procedures with
#' censored survival data. \emph{Statistics in Medicine}, 30(10), 1105-1117.
#' \doi{10.1002/sim.4154}.
#' 
#' Ishwaran H, Kogalur UB, Blackstone EH, Lauer MS (2008). Random survival
#' forests. \emph{The Annals of Applied Statistics}, 2(3), 841-860.
#' \doi{10.1214/08-AOAS169}.
#' 
#' Wright MN, Ziegler A (2017). ranger: A fast implementation of random forests
#' for high dimensional data in C++ and R. \emph{Journal of Statistical Software},
#' 77(1), 1-17. \doi{10.18637/jss.v077.i01}.
#' 
#' Royston P, Parmar MKB, Sylvester R (2004). Construction and
#' validation of a prognostic model across several studies, with an
#' application in superficial bladder cancer. \emph{Statistics in Medicine},
#' 23, 907-926. \doi{10.1002/sim.1691}.
#' 
#' Debray TPA, Moons KGM, Ahmed I, Koffijberg H, Riley RD (2013). A framework
#' for developing, implementing, and evaluating clinical prediction models
#' in an individual participant data meta-analysis. \emph{Statistics in Medicine},
#' 32, 3158-3180. \doi{10.1002/sim.5732}.
#' 
#' Steyerberg EW, Harrell FE Jr (2016). Prediction models need appropriate
#' internal, internal-external, and external validation. \emph{Journal of Clinical
#' Epidemiology}, 69, 245-247. \doi{10.1016/j.jclinepi.2015.04.005}.
#' }
#' @rawRd \examples{
#' # sb_compare_cv
#' d <- sb_example_survival(n = 50, seed = 1)
#' task <- sb_surv_task(
#'   d, "patient_id", "time", "event",
#'   c("gene_a", "gene_b"), "month"
#' )
#' compared <- sb_compare_cv(
#'   task,
#'   list(
#'     sb_pipeline(),
#'     sb_pipeline(
#'       selector = "univ_cox_top_k",
#'       selector_n_features = 1
#'     )
#'   ),
#'   folds = 5,
#'   seed = 2026
#' )
#' compared
#' outer_plan <- setNames(
#'   compared$fold_assignment$fold,
#'   compared$fold_assignment$patient_id
#' )
#' identical(unname(outer_plan), compared$fold_assignment$fold)
#' 
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   binary_data <- sb_example_binary(n = 60, seed = 1)
#'   binary_task <- sb_binary_task(
#'     binary_data, "patient_id", "outcome",
#'     c("gene_a", "gene_b", "gene_c")
#'   )
#'   sb_compare_cv(
#'     binary_task,
#'     list(
#'       sb_pipeline(learner = "logistic"),
#'       sb_pipeline(learner = "gaussian_naive_bayes"),
#'       sb_pipeline(
#'         learner = "glmnet_logistic",
#'         inner_folds = 3,
#'         inner_seed = 1701
#'       )
#'     ),
#'     folds = 3,
#'     seed = 2026
#'   )
#' }
#' 
#' if (requireNamespace("lightgbm", quietly = TRUE)) {
#'   binary_data <- sb_example_binary(n = 50, seed = 2)
#'   binary_task <- sb_binary_task(
#'     binary_data, "patient_id", "outcome",
#'     c("gene_a", "gene_b", "gene_c")
#'   )
#'   sb_compare_cv(
#'     binary_task,
#'     list(
#'       sb_pipeline(
#'         learner = "lightgbm_binary",
#'         lightgbm_nrounds = 10
#'       ),
#'       sb_pipeline(
#'         selector = "univ_binary_auc_top_k",
#'         learner = "lightgbm_binary",
#'         selector_n_features = 2,
#'         lightgbm_nrounds = 10
#'       )
#'     ),
#'     folds = 5,
#'     seed = 2026
#'   )
#' }
#' 
#' if (requireNamespace("catboost", quietly = TRUE)) {
#'   binary_data <- sb_example_binary(n = 50, seed = 3)
#'   binary_task <- sb_binary_task(
#'     binary_data, "patient_id", "outcome",
#'     c("gene_a", "gene_b", "gene_c")
#'   )
#'   sb_compare_cv(
#'     binary_task,
#'     list(
#'       sb_pipeline(
#'         learner = "catboost_binary",
#'         catboost_iterations = 10
#'       ),
#'       sb_pipeline(
#'         selector = "univ_binary_auc_top_k",
#'         learner = "catboost_binary",
#'         selector_n_features = 2,
#'         catboost_iterations = 10
#'       )
#'     ),
#'     folds = 5,
#'     seed = 2026
#'   )
#' }
#' 
#' # sb_run_cv
#' d <- sb_example_survival(n = 50, seed = 1)
#' task <- sb_surv_task(
#'   d, "patient_id", "time", "event",
#'   c("gene_a", "gene_b"), "month",
#'   prediction_horizon = 6
#' )
#' result <- sb_run_cv(task, folds = 5, seed = 2026)
#' result
#' outer_plan <- setNames(
#'   result$provenance$fold_assignment$fold,
#'   result$provenance$fold_assignment$patient_id
#' )
#' reused <- sb_run_cv(task, folds = rev(outer_plan), seed = 999)
#' identical(
#'   reused$provenance$fold_assignment$fold,
#'   result$provenance$fold_assignment$fold
#' )
#' 
#' binary_data <- sb_example_binary(n = 50, seed = 1)
#' binary_task <- sb_binary_task(
#'   binary_data, "patient_id", "outcome",
#'   c("gene_a", "gene_b")
#' )
#' sb_run_cv(binary_task, folds = 5, seed = 2026)
#' sb_run_cv(
#'   binary_task,
#'   sb_pipeline(learner = "gaussian_naive_bayes"),
#'   folds = 5,
#'   seed = 2026
#' )
#' if (requireNamespace("lightgbm", quietly = TRUE)) {
#'   sb_run_cv(
#'     binary_task,
#'     sb_pipeline(
#'       learner = "lightgbm_binary",
#'       lightgbm_nrounds = 10
#'     ),
#'     folds = 5,
#'     seed = 2026
#'   )
#' }
#' if (requireNamespace("catboost", quietly = TRUE)) {
#'   sb_run_cv(
#'     binary_task,
#'     sb_pipeline(
#'       learner = "catboost_binary",
#'       catboost_iterations = 10
#'     ),
#'     folds = 5,
#'     seed = 2026
#'   )
#' }
#' 
#' # sb_compare_repeated_cv
#' d <- sb_example_binary(n = 48, seed = 1)
#' task <- sb_binary_task(
#'   d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
#' )
#' repeated <- sb_compare_repeated_cv(
#'   task,
#'   list(
#'     sb_pipeline(learner = "logistic"),
#'     sb_pipeline(learner = "gaussian_naive_bayes")
#'   ),
#'   folds = 3,
#'   seeds = c(11, 12)
#' )
#' repeated
#' 
#' # sb_compare_cohorts
#' make_task <- function(prefix, seed, n) {
#'   d <- sb_example_binary(n, seed)
#'   d$patient_id <- paste0(prefix, d$patient_id)
#'   sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#' }
#' result <- sb_compare_cohorts(
#'   make_task("dev_", 1, 50),
#'   list(Validation = make_task("val_", 2, 30),
#'        Test = make_task("test_", 3, 30)),
#'   roles = c(Validation = "validation", Test = "test"),
#'   models = c("logistic", "gaussian_naive_bayes"),
#'   folds = 3, view = FALSE
#' )
#' result$cohorts$Test$comparison
#' 
#' # sb_run_internal_external_cv
#' a <- sb_example_binary(n = 24, seed = 1)
#' b <- sb_example_binary(n = 24, seed = 2)
#' a$patient_id <- paste0("a_", a$patient_id)
#' b$patient_id <- paste0("b_", b$patient_id)
#' predictors <- c("gene_a", "gene_b")
#' cohorts <- list(
#'   cohort_a = sb_binary_task(
#'     a, "patient_id", "outcome", predictors
#'   ),
#'   cohort_b = sb_binary_task(
#'     b, "patient_id", "outcome", predictors
#'   )
#' )
#' sb_run_internal_external_cv(
#'   cohorts,
#'   list(sb_pipeline(learner = "logistic"))
#' )
#' 
#' # bi_run
#' genes <- c("gene_a", "gene_b")
#' development <- bi_example(40, seed = 42)
#' binary <- bi_run(development, outcome = "outcome", id = "patient_id",
#'   predictors = genes, models = "Logistic", folds = 2, seed = 42,
#'   plots = FALSE)
#' result_table(binary$comparison, metric = "roc_auc")
#' 
#' follow_up <- surv_example(40, seed = 42)
#' survival <- surv_run(follow_up, time = "time", event = "event",
#'   id = "patient_id", predictors = genes, time_unit = "month",
#'   models = "Cox", folds = 2, seed = 42, plots = FALSE)
#' result_table(survival$comparison, metric = "harrell_c")
#' 
#' # sb_with_selector_reuse
#' if (requireNamespace("Boruta", quietly = TRUE) &&
#'     requireNamespace("ranger", quietly = TRUE)) {
#'   d <- sb_example_binary(n = 100, seed = 42)
#'   task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#'   shared <- sb_with_selector_reuse({
#'     lapply(c("logistic", "gaussian_naive_bayes"), function(learner) {
#'       pipeline <- sb_pipeline("boruta_binary_roughfix", learner)
#'       sb_run_cv(task, pipeline, folds = 3, seed = 42)
#'     })
#'   }, diagnostics = TRUE)
#'   shared$reuse
#' }
#' }
NULL

# ---- cv-comparison ----

# ---- compare-cv ----

#' @rdname sb_compare_cv
#' @export
sb_compare_cv <- function(task, pipelines, folds = 5L, seed = 1L,
                          reuse_selectors = TRUE, progress = interactive()) {
  if (inherits(pipelines, "stormbreaker_stack")) {
    pipelines <- sb_upgrade_legacy_names(pipelines)
  } else if (is.list(pipelines)) pipelines <- lapply(pipelines, sb_upgrade_legacy_names)

  sb_check_progress(progress)
  task_type <- if (inherits(task, "stormbreaker_binary_task")) {
    "binary"
  } else if (inherits(task, "stormbreaker_surv_task")) {
    "survival"
  } else {
    sb_abort(
      "`task` must be created by `sb_surv_task()` or `sb_binary_task()`.",
      "stormbreaker_bad_task"
    )
  }
  valid_collection <- is.list(pipelines) && length(pipelines) > 0L &&
    all(vapply(
      pipelines,
      function(pipeline) inherits(pipeline, "stormbreaker_pipeline"),
      logical(1)
    ))
  if (!valid_collection) {
    sb_abort(
      "`pipelines` must be a non-empty list created with `sb_pipeline()`.",
      "stormbreaker_bad_pipeline_collection"
    )
  }

  pipeline_ids <- vapply(
    pipelines,
    function(pipeline) pipeline$id,
    character(1)
  )
  if (anyDuplicated(pipeline_ids) > 0L) {
    sb_abort(
      "Every compared pipeline must have a distinct pipeline ID.",
      "stormbreaker_duplicate_pipeline_id"
    )
  }
  matching_task <- vapply(
    pipelines,
    function(pipeline) identical(pipeline$task_type, task_type),
    logical(1)
  )
  if (!all(matching_task)) {
    sb_abort(
      "Every compared pipeline must match the supplied task type.",
      "stormbreaker_pipeline_task_mismatch"
    )
  }

  results <- sb_with_selector_reuse(lapply(
    seq_along(pipelines),
    function(i) {
      pipeline <- pipelines[[i]]
      sb_progress_message(progress, "CV model ", i, "/", length(pipelines), ": ",
        sb_display_selector(pipeline$selector), " + ",
        sb_display_learner(pipeline$learner, pipeline$learner_params))
      result <- sb_run_cv(task, pipeline, folds = folds, seed = seed, progress = progress)
      sb_progress_message(progress, "CV model ", i, "/", length(pipelines),
        " complete: ", result$status)
      result
    }
  ), enabled = reuse_selectors)
  names(results) <- pipeline_ids

  rows <- lapply(seq_along(results), function(index) {
    result <- results[[index]]
    pipeline <- pipelines[[index]]
    if (identical(task_type, "binary")) {
      roc_auc <- result$metrics$roc_auc
      average_precision <- result$metrics$average_precision
      brier <- result$metrics$brier_score
      log_loss <- result$metrics$log_loss
      calibration <- result$metrics$calibration
      moderate_calibration <- result$metrics$moderate_calibration
      selection_stability <- result$selection_stability
      active_set_stability <- result$glmnet_active_set_stability
      data.frame(
        pipeline_id = pipeline$id,
        selector = pipeline$selector,
        learner = pipeline$learner,
        status = result$status,
        complete_oof = result$complete_oof,
        roc_auc_status = roc_auc$status,
        roc_auc = roc_auc$estimate,
        comparable_pairs = roc_auc$comparable_pairs,
        average_precision_status = average_precision$status,
        average_precision = average_precision$estimate,
        brier_status = brier$status,
        brier_score = brier$estimate,
        n_scored = brier$n_scored,
        reference_brier_status = brier$reference$status,
        reference_brier = brier$reference$estimate,
        brier_skill_status = brier$skill_score$status,
        brier_skill_score = brier$skill_score$estimate,
        log_loss_status = log_loss$status,
        log_loss = log_loss$estimate,
        reference_log_loss_status = log_loss$reference$status,
        reference_log_loss = log_loss$reference$estimate,
        log_loss_skill_status = log_loss$skill_score$status,
        log_loss_skill_score = log_loss$skill_score$estimate,
        calibration_status = calibration$status,
        calibration_intercept = calibration$calibration_intercept,
        calibration_slope = calibration$calibration_slope,
        moderate_calibration_status = moderate_calibration$status,
        moderate_calibration_ici = moderate_calibration$ici,
        moderate_calibration_e50 = moderate_calibration$e50,
        moderate_calibration_e90 = moderate_calibration$e90,
        selection_stability_status = selection_stability$status,
        mean_pairwise_jaccard = selection_stability$mean_pairwise_jaccard,
        n_successful_selection_folds = selection_stability$n_successful_folds,
        n_folds = as.integer(nrow(result$runs)),
        n_successful_folds = as.integer(
          sum(result$runs$status == "success")
        ),
        glmnet_active_set_stability_status = active_set_stability$status,
        glmnet_active_set_mean_pairwise_jaccard =
          active_set_stability$mean_pairwise_jaccard,
        n_successful_glmnet_active_set_folds =
          active_set_stability$n_successful_folds,
        stringsAsFactors = FALSE
      )
    } else {
      metric <- result$metrics$harrell_c
      selection_stability <- result$selection_stability
      active_set_stability <- result$glmnet_active_set_stability
      row <- data.frame(
        pipeline_id = pipeline$id,
        selector = pipeline$selector,
        learner = pipeline$learner,
        status = result$status,
        complete_oof = result$complete_oof,
        harrell_c_status = metric$status,
        harrell_c = metric$estimate,
        comparable_pairs = metric$comparable_pairs,
        selection_stability_status = selection_stability$status,
        mean_pairwise_jaccard = selection_stability$mean_pairwise_jaccard,
        n_successful_selection_folds = selection_stability$n_successful_folds,
        n_folds = as.integer(nrow(result$runs)),
        n_successful_folds = as.integer(
          sum(result$runs$status == "success")
        ),
        glmnet_active_set_stability_status = active_set_stability$status,
        glmnet_active_set_mean_pairwise_jaccard =
          active_set_stability$mean_pairwise_jaccard,
        n_successful_glmnet_active_set_folds =
          active_set_stability$n_successful_folds,
        stringsAsFactors = FALSE
      )
      if (!is.null(result$metrics$ipcw_brier)) {
        ipcw_brier <- result$metrics$ipcw_brier
        weak_calibration <- result$metrics$weak_calibration
        row <- data.frame(
          pipeline_id = row$pipeline_id,
          selector = row$selector,
          learner = row$learner,
          status = row$status,
          complete_oof = row$complete_oof,
          complete_survival_probability_oof =
            result$complete_survival_probability_oof,
          prediction_horizon = ipcw_brier$prediction_horizon,
          harrell_c_status = row$harrell_c_status,
          harrell_c = row$harrell_c,
          comparable_pairs = row$comparable_pairs,
          ipcw_brier_status = ipcw_brier$status,
          ipcw_brier = ipcw_brier$estimate,
          reference_ipcw_brier_status = ipcw_brier$reference$status,
          reference_ipcw_brier = ipcw_brier$reference$estimate,
          ipcw_brier_skill_status = ipcw_brier$skill_score$status,
          ipcw_brier_skill_score = ipcw_brier$skill_score$estimate,
          time_auc_status = result$metrics$time_auc$status,
          time_auc = result$metrics$time_auc$estimate,
          time_auc_pair_denominator =
            result$metrics$time_auc$pair_denominator,
          uno_c_status = result$metrics$uno_c$status,
          uno_c = result$metrics$uno_c$estimate,
          uno_c_weighted_pair_denominator =
            result$metrics$uno_c$weighted_pair_denominator,
          calibration_status = result$metrics$calibration$status,
          calibration_ici = result$metrics$calibration$ici,
          calibration_e50 = result$metrics$calibration$e50,
          calibration_e90 = result$metrics$calibration$e90,
          mean_calibration_status = weak_calibration$mean_calibration$status,
          mean_calibration_ratio =
            weak_calibration$mean_calibration_ratio,
          calibration_slope_status = weak_calibration$slope_fit$status,
          calibration_slope = weak_calibration$calibration_slope,
          selection_stability_status = row$selection_stability_status,
          mean_pairwise_jaccard = row$mean_pairwise_jaccard,
          n_successful_selection_folds =
            row$n_successful_selection_folds,
          n_folds = row$n_folds,
          n_successful_folds = row$n_successful_folds,
          glmnet_active_set_stability_status =
            row$glmnet_active_set_stability_status,
          glmnet_active_set_mean_pairwise_jaccard =
            row$glmnet_active_set_mean_pairwise_jaccard,
          n_successful_glmnet_active_set_folds =
            row$n_successful_glmnet_active_set_folds,
          stringsAsFactors = FALSE
        )
      }
      row
    }
  })
  comparison <- do.call(rbind, rows)
  rownames(comparison) <- NULL
  comparison$prediction_variation <- vapply(results, function(result) {
    sb_prediction_variation(result$predictions, task_type)$status
  }, character(1))
  metric_status <- if (identical(task_type, "binary")) {
    comparison$roc_auc_status == "success" &
      comparison$average_precision_status == "success" &
      comparison$brier_status == "success" &
      comparison$reference_brier_status == "success" &
      comparison$brier_skill_status == "success" &
      comparison$log_loss_status == "success" &
      comparison$reference_log_loss_status == "success" &
      comparison$log_loss_skill_status == "success" &
      comparison$calibration_status == "success" &
      comparison$selection_stability_status %in% c(
        "success", "not_applicable"
      )
  } else {
    survival_metric_status <- comparison$harrell_c_status == "success" &
      comparison$selection_stability_status %in% c(
        "success", "not_applicable"
      )
    if ("ipcw_brier_status" %in% names(comparison)) {
      survival_metric_status <- survival_metric_status &
        comparison$ipcw_brier_status == "success" &
        comparison$reference_ipcw_brier_status == "success" &
        comparison$ipcw_brier_skill_status == "success" &
        comparison$time_auc_status == "success" &
        comparison$uno_c_status == "success"
    }
    survival_metric_status
  }
  complete <- all(comparison$status == "success") &&
    all(metric_status)

  compared <- list(
    status = if (complete) "success" else "incomplete",
    task_type = task_type,
    comparison = comparison,
    results = results,
    n_patients = results[[1L]]$task$n,
    folds = results[[1L]]$provenance$folds,
    seed = results[[1L]]$provenance$seed,
    fold_source = results[[1L]]$provenance$fold_source,
    fold_assignment = results[[1L]]$provenance$fold_assignment
  )
  if ("decision_thresholds" %in% names(task)) {
    decision_curves <- lapply(
      results,
      function(result) result$metrics$decision_curve
    )
    names(decision_curves) <- names(results)
    compared$decision_curves <- decision_curves
  }
  structure(compared, class = "stormbreaker_cv_comparison")
}

#' @export
print.stormbreaker_cv_comparison <- function(x, ...) {
  cat("<stormbreaker_cv_comparison:", x$status, ">\n")
  cat(
    "Pipelines:", nrow(x$comparison),
    "| Patients:", x$n_patients,
    "| Shared outer folds:", x$folds, "\n"
  )
  fold_label <- if (identical(x$fold_source, "preassigned_patient_id")) {
    "preassigned patient-ID folds"
  } else {
    "generated stratified folds"
  }
  cat("Outer resampling:", fold_label, "\n")
  metric_columns <- if (identical(x$task_type, "binary")) {
    c(
      "roc_auc", "average_precision", "brier_score", "brier_skill_score",
      "log_loss", "log_loss_skill_score",
      "calibration_intercept", "calibration_slope",
      "moderate_calibration_ici", "moderate_calibration_e50",
      "moderate_calibration_e90",
      "comparable_pairs", "mean_pairwise_jaccard",
      "glmnet_active_set_mean_pairwise_jaccard"
    )
  } else {
    survival_columns <- c(
      "harrell_c", "comparable_pairs", "mean_pairwise_jaccard",
      "glmnet_active_set_mean_pairwise_jaccard"
    )
    if ("ipcw_brier" %in% names(x$comparison)) {
      survival_columns <- c(
        "harrell_c", "ipcw_brier", "reference_ipcw_brier",
        "ipcw_brier_skill_score", "time_auc", "uno_c", "calibration_ici",
        "calibration_e50", "calibration_e90", "mean_calibration_ratio",
        "calibration_slope",
        "time_auc_pair_denominator", "uno_c_weighted_pair_denominator",
        "comparable_pairs",
        "mean_pairwise_jaccard",
        "glmnet_active_set_mean_pairwise_jaccard"
      )
    }
    survival_columns
  }
  print(
    x$comparison[c(
      "pipeline_id", "status", metric_columns
    )],
    row.names = FALSE,
    digits = 3
  )
  cat("No pipeline is selected by this comparison.\n")
  invisible(x)
}

# ---- compare-repeated-cv ----

sb_repeated_cv_specifications <- function(folds, seeds) {
  if (is.list(folds)) {
    if (length(folds) == 0L) {
      sb_abort(
        "`folds` must contain at least one patient-ID fold assignment.",
        "stormbreaker_bad_repeated_folds"
      )
    }
    if (!is.null(seeds)) {
      sb_abort(
        paste0(
          "`seeds` must be `NULL` when `folds` is a list of explicit ",
          "patient-ID assignments."
        ),
        "stormbreaker_repeated_fold_seed_conflict"
      )
    }
    valid_plan <- vapply(folds, function(plan) {
      plan_names <- names(plan)
      is.numeric(plan) && is.null(dim(plan)) && length(plan) > 1L &&
        !is.null(plan_names) && !anyNA(plan_names) &&
        all(nzchar(plan_names)) && anyDuplicated(plan_names) == 0L
    }, logical(1L))
    if (!all(valid_plan)) {
      sb_abort(
        paste0(
          "Every `folds` list element must be one multi-patient numeric ",
          "vector named by unique, non-empty patient IDs."
        ),
        "stormbreaker_bad_repeated_folds"
      )
    }
    return(list(
      folds = folds,
      seeds = rep(NA_integer_, length(folds))
    ))
  }

  if (length(folds) != 1L) {
    sb_abort(
      paste0(
        "For repeated generated folds, `folds` must be one fold count. ",
        "Wrap explicit patient-ID assignments in a list."
      ),
      "stormbreaker_bad_repeated_folds"
    )
  }
  valid_seeds <- is.numeric(seeds) && length(seeds) > 0L &&
    !anyNA(seeds) && all(is.finite(seeds)) && all(seeds > 0) &&
    all(seeds <= .Machine$integer.max) && all(seeds == floor(seeds))
  if (!valid_seeds) {
    sb_abort(
      "`seeds` must contain one or more finite positive whole numbers.",
      "stormbreaker_bad_repeat_seeds"
    )
  }

  list(
    folds = rep(list(folds), length(seeds)),
    seeds = as.integer(seeds)
  )
}

sb_run_repeated_cv_iteration <- function(task, pipelines, folds, seed) {
  if (inherits(pipelines, "stormbreaker_stack")) {
    return(sb_run_stack_cv(task, pipelines, folds = folds, seed = seed))
  }
  sb_compare_cv(task, pipelines, folds = folds, seed = seed)
}

sb_repeated_cv_result_metadata <- function(result) {
  if (inherits(result, "stormbreaker_stack_cv_result")) {
    return(list(
      folds = result$provenance$folds,
      seed = result$provenance$seed,
      fold_source = result$provenance$fold_source,
      fold_assignment = result$provenance$fold_assignment,
      n_patients = result$task$n
    ))
  }
  list(
    folds = result$folds,
    seed = result$seed,
    fold_source = result$fold_source,
    fold_assignment = result$fold_assignment,
    n_patients = result$n_patients
  )
}

sb_learner_uses_caller_rng <- function(id) {
  id %in% c("blackboost_cox", "rsf_survival", "adaboost_binary")
}

sb_parallel_cv_repetitions <- function(task, pipelines, specifications) {
  members <- if (inherits(pipelines, "stormbreaker_stack")) pipelines$pipelines else pipelines
  if (is.list(members) && any(vapply(members, function(p) {
    inherits(p, "stormbreaker_pipeline") && sb_learner_uses_caller_rng(p$learner)
  }, logical(1)))) {
    sb_abort(paste0("These pipelines include a learner that consumes the caller's ",
      "sequential RNG stream; use workers = 1 to preserve its results."),
      "stormbreaker_parallel_sequential_rng")
  }
  threads <- c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1")
  cluster <- withr::with_envvar(threads, parallel::makePSOCKcluster(
    min(2L, length(specifications$folds))))
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  initialize <- function(paths, package_path, rng, kind) {
    .libPaths(paths)
    options(rf.cores = 1L, mc.cores = 1L,
      stormbreaker.ranger_cindex_threads = 1L)
    do.call(RNGkind, as.list(kind))
    if (is.null(rng)) {
      if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", rng, .GlobalEnv)
    if (file.exists(file.path(package_path, "Meta", "package.rds"))) {
      loadNamespace("stormbreaker", lib.loc = dirname(package_path))
    } else {
      pkgload::load_all(package_path, quiet = TRUE, helpers = FALSE, attach = FALSE)
    }
    invisible(NULL)
  }
  environment(initialize) <- baseenv()
  rng <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    get(".Random.seed", .GlobalEnv) else NULL
  parallel::clusterCall(cluster, initialize, .libPaths(),
    getNamespaceInfo(asNamespace("stormbreaker"), "path"), rng, RNGkind())
  run <- function(i, task, pipelines, specifications, rng) {
    if (!is.null(rng)) {
      assign(".Random.seed", rng, .GlobalEnv)
    } else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    warnings <- list()
    value <- tryCatch(withCallingHandlers({
      seed <- specifications$seeds[[i]]
      get("sb_run_repeated_cv_iteration", asNamespace("stormbreaker"))(
        task, pipelines, specifications$folds[[i]], if (is.na(seed)) 1L else seed)
    }, warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    list(value = value, warnings = warnings, failed = inherits(value, "error"))
  }
  environment(run) <- baseenv()
  outputs <- parallel::parLapplyLB(cluster, seq_along(specifications$folds), run,
    task = task, pipelines = pipelines, specifications = specifications, rng = rng)
  lapply(outputs, function(output) {
    for (condition in output$warnings) warning(condition)
    if (output$failed) stop(output$value)
    output$value
  })
}

#' @rdname sb_compare_cv
#' @export
sb_compare_repeated_cv <- function(task, pipelines, folds = 5L,
                                   seeds = 1:5, workers = 1L) {
  if (inherits(pipelines, "stormbreaker_stack")) {
    pipelines <- sb_upgrade_legacy_names(pipelines)
  } else if (is.list(pipelines)) pipelines <- lapply(pipelines, sb_upgrade_legacy_names)

  if (!is.numeric(workers) || length(workers) != 1L || is.na(workers) ||
      !workers %in% c(1L, 2L)) {
    sb_abort("`workers` must be 1 or 2.", "stormbreaker_bad_workers")
  }
  specifications <- sb_repeated_cv_specifications(folds, seeds)
  repetitions <- if (workers == 2L) {
    sb_parallel_cv_repetitions(task, pipelines, specifications)
  } else lapply(
    seq_along(specifications$folds),
    function(repeat_id) {
      seed <- specifications$seeds[[repeat_id]]
      sb_run_repeated_cv_iteration(
        task,
        pipelines,
        folds = specifications$folds[[repeat_id]],
        seed = if (is.na(seed)) 1L else seed
      )
    }
  )
  names(repetitions) <- paste0("repeat_", seq_along(repetitions))
  metadata <- lapply(repetitions, sb_repeated_cv_result_metadata)

  comparison_rows <- lapply(seq_along(repetitions), function(repeat_id) {
    repetition <- repetitions[[repeat_id]]
    repeat_metadata <- metadata[[repeat_id]]
    data.frame(
      repeat_id = as.integer(repeat_id),
      folds = repeat_metadata$folds,
      seed = repeat_metadata$seed,
      fold_source = repeat_metadata$fold_source,
      repetition$comparison,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  comparison <- do.call(rbind, comparison_rows)
  rownames(comparison) <- NULL

  fold_rows <- lapply(seq_along(repetitions), function(repeat_id) {
    data.frame(
      repeat_id = as.integer(repeat_id),
      metadata[[repeat_id]]$fold_assignment,
      stringsAsFactors = FALSE
    )
  })
  fold_assignments <- do.call(rbind, fold_rows)
  rownames(fold_assignments) <- NULL

  complete <- all(vapply(
    repetitions,
    function(repetition) identical(repetition$status, "success"),
    logical(1L)
  ))
  result <- list(
    status = if (complete) "success" else "incomplete",
    task_type = repetitions[[1L]]$task_type,
    comparison = comparison,
    repetitions = repetitions,
    n_patients = metadata[[1L]]$n_patients,
    n_repeats = as.integer(length(repetitions)),
    fold_assignments = fold_assignments
  )
  if (inherits(pipelines, "stormbreaker_stack")) {
    result$strategy_type <- "stack"
  }
  structure(result, class = "stormbreaker_repeated_cv_comparison")
}

#' @export
print.stormbreaker_repeated_cv_comparison <- function(x, ...) {
  cat("<stormbreaker_repeated_cv_comparison:", x$status, ">\n")
  count_label <- if (identical(x$strategy_type, "stack")) {
    "Member/stack rows:"
  } else {
    "Pipelines:"
  }
  cat(
    count_label, length(unique(x$comparison$pipeline_id)),
    "| Patients:", x$n_patients,
    "| Outer repetitions:", x$n_repeats, "\n"
  )
  metric_columns <- if (identical(x$task_type, "binary")) {
    c("roc_auc", "average_precision", "brier_score", "log_loss")
  } else {
    columns <- "harrell_c"
    if ("ipcw_brier" %in% names(x$comparison)) {
      columns <- c(columns, "ipcw_brier", "time_auc", "uno_c")
    }
    columns
  }
  print(
    x$comparison[c(
      "repeat_id", "pipeline_id", "status", metric_columns
    )],
    row.names = FALSE,
    digits = 3
  )
  cat(
    paste0(
      "No cross-repeat mean, confidence interval, ranking, or winner ",
      "is computed.\n"
    )
  )
  invisible(x)
}

# ---- cohort-comparison ----

# ---- compare-cohorts ----

#' @rdname sb_compare_cv
#' @export
sb_compare_cohorts <- function(development, cohorts, roles = NULL,
                               pipelines = NULL, models = NULL,
                               folds = 5L, seed = 1L, inner_folds = 5L,
                               inner_seed = 1701L, development_name = "Development",
                               view = FALSE, report_file = NULL,
                               progress = interactive()) {
  sb_check_progress(progress)
  sb_check_page_request(view, report_file)
  if (!is.character(development_name) || length(development_name) != 1L ||
      is.na(development_name) || !nzchar(trimws(development_name))) {
    sb_abort("Supply one non-empty development_name.", "stormbreaker_bad_cohort_roles")
  }
  if (!is.list(cohorts) || is.object(cohorts) || !length(cohorts) ||
      is.null(names(cohorts)) || development_name %in% names(cohorts)) {
    sb_abort("Supply named cohorts distinct from development_name.",
      "stormbreaker_bad_cohort_roles")
  }
  if (is.null(roles)) roles <- stats::setNames(rep("external", length(cohorts)), names(cohorts))
  if (!is.character(roles) || anyNA(roles) || is.null(names(roles)) ||
      anyDuplicated(names(roles)) || !setequal(names(roles), names(cohorts)) ||
      !all(roles %in% c("validation", "test", "external"))) {
    sb_abort("roles must name every cohort as validation, test or external.",
      "stormbreaker_bad_cohort_roles")
  }
  roles <- roles[names(cohorts)]
  validated <- sb_validate_iecv_cohorts(c(stats::setNames(list(development), development_name), cohorts))
  development <- validated$tasks[[1L]]
  cohorts <- validated$tasks[-1L]
  if (!is.null(pipelines) && !is.null(models)) {
    sb_abort("Supply pipelines or models, not both.", "stormbreaker_bad_model_choice")
  }
  if (is.null(pipelines)) {
    pipelines <- sb_model_pipelines(validated$task_type, models, inner_folds, inner_seed)
  }
  pipelines <- sb_validate_iecv_pipelines(pipelines, validated$task_type)
  started <- proc.time()[["elapsed"]]
  cv <- sb_compare_cv(development, pipelines, folds = folds, seed = seed, progress = progress)
  fits <- sb_with_selector_reuse(lapply(seq_along(pipelines), function(i) {
    pipeline <- pipelines[[i]]
    sb_progress_message(progress, "Final fit ", i, "/", length(pipelines), ": ",
      sb_display_selector(pipeline$selector), " + ",
      sb_display_learner(pipeline$learner, pipeline$learner_params))
    fit <- sb_fit_final(development, pipeline)
    sb_progress_message(progress, "Final fit ", i, "/", length(pipelines), " complete: ", fit$status)
    fit
  }))
  names(fits) <- names(pipelines)
  predictions <- lapply(seq_along(cohorts), function(i) {
    cohort <- cohorts[[i]]
    sb_progress_message(progress, "Predict ", names(cohorts)[i])
    new_data <- cohort$data[c(cohort$id_col, development$predictors)]
    lapply(fits, function(fit) sb_predict_final(fit, new_data, cohort$id_col))
  })
  names(predictions) <- names(cohorts)
  evaluations <- lapply(seq_along(cohorts), function(i) {
    sb_progress_message(progress, "Evaluate ", names(cohorts)[i])
    endpoint <- sb_iecv_endpoint(cohorts[[i]], validated$task_type)
    sb_compare_external(predictions[[i]], endpoint)
  })
  names(evaluations) <- names(cohorts)
  result <- list(
    task_type = validated$task_type, development = cv, fits = fits,
    predictions = predictions, cohorts = evaluations,
    results = c(stats::setNames(list(cv), development_name), evaluations),
    cohort_roles = c(stats::setNames("development_cv", development_name), roles),
    model_pipelines = pipelines,
    elapsed_seconds = unname(proc.time()[["elapsed"]] - started)
  )
  class(result) <- "stormbreaker_cohort_comparison"
  sb_attach_result_page(result, view, report_file)
}

#' @export
print.stormbreaker_cohort_comparison <- function(x, ...) {
  cat("<stormbreaker_cohort_comparison>\n")
  cat("Task:", x$task_type, "| Pipelines:", length(x$model_pipelines),
    "| Evaluation cohorts:", length(x$cohorts), "\n")
  cat("Results:", paste(names(x$results), collapse = ", "), "\n")
  if (!is.null(x$report_file)) cat("Results page:", x$report_file, "\n")
  invisible(x)
}

# ---- internal-external-cv ----

sb_iecv_optional_common <- function(tasks, field, description) {
  present <- vapply(tasks, function(task) field %in% names(task), logical(1L))
  if (!all(present == present[[1L]])) {
    sb_abort(
      paste0("All cohorts must use the same ", description, "."),
      "stormbreaker_internal_external_target_mismatch"
    )
  }
  if (!present[[1L]]) {
    return(invisible(NULL))
  }
  reference <- tasks[[1L]][[field]]
  same <- vapply(tasks, function(task) {
    identical(task[[field]], reference)
  }, logical(1L))
  if (!all(same)) {
    sb_abort(
      paste0("All cohorts must use the same ", description, "."),
      "stormbreaker_internal_external_target_mismatch"
    )
  }
  invisible(reference)
}

sb_validate_iecv_cohorts <- function(cohorts) {
  valid_collection <- is.list(cohorts) && length(cohorts) >= 2L &&
    !inherits(cohorts, c(
      "stormbreaker_binary_task", "stormbreaker_surv_task"
    ))
  cohort_ids <- names(cohorts)
  valid_ids <- !is.null(cohort_ids) &&
    length(cohort_ids) == length(cohorts) &&
    !anyNA(cohort_ids) && all(nzchar(trimws(cohort_ids))) &&
    anyDuplicated(cohort_ids) == 0L
  if (!valid_collection || !valid_ids) {
    sb_abort(
      paste0(
        "`cohorts` must be a named list of at least two tasks with unique, ",
        "non-blank cohort IDs."
      ),
      "stormbreaker_bad_internal_external_cohorts"
    )
  }

  task_types <- vapply(cohorts, function(task) {
    if (inherits(task, "stormbreaker_binary_task")) {
      return("binary")
    }
    if (inherits(task, "stormbreaker_surv_task")) {
      return("survival")
    }
    sb_abort(
      paste0(
        "Every cohort must be created by `sb_binary_task()` or ",
        "`sb_surv_task()`."
      ),
      "stormbreaker_bad_internal_external_cohorts"
    )
  }, character(1L))
  if (length(unique(task_types)) != 1L) {
    sb_abort(
      "All internal-external validation cohorts must use the same task type.",
      "stormbreaker_internal_external_task_mismatch"
    )
  }

  task_type <- task_types[[1L]]
  tasks <- lapply(cohorts, function(task) {
    if (identical(task_type, "binary")) {
      sb_validate_binary_task_contract(task)
    } else {
      sb_validate_task_contract(task)
    }
  })
  names(tasks) <- cohort_ids

  reference_predictors <- tasks[[1L]]$predictors
  same_predictors <- vapply(tasks, function(task) {
    length(task$predictors) == length(reference_predictors) &&
      setequal(task$predictors, reference_predictors)
  }, logical(1L))
  if (!all(same_predictors)) {
    sb_abort(
      paste0(
        "All cohorts must declare the same predictor-name set; ",
        "STORMBREAKER does not intersect or harmonize features."
      ),
      "stormbreaker_internal_external_predictor_mismatch"
    )
  }

  if (identical(task_type, "survival")) {
    time_units <- vapply(tasks, `[[`, character(1L), "time_unit")
    if (length(unique(time_units)) != 1L) {
      sb_abort(
        "All survival cohorts must use the same time unit.",
        "stormbreaker_internal_external_target_mismatch"
      )
    }
    sb_iecv_optional_common(
      tasks, "prediction_horizon", "prediction-horizon declaration"
    )
  }
  sb_iecv_optional_common(
    tasks, "decision_thresholds", "decision-threshold declaration"
  )
  grouped <- vapply(tasks, function(task) !is.null(task$group_col), logical(1))
  if (any(grouped) && !all(grouped)) {
    sb_abort(
      "All cohorts must declare grouping when any cohort declares grouping.",
      "stormbreaker_internal_external_target_mismatch"
    )
  }
  if (all(grouped)) {
    groups <- unlist(lapply(tasks, function(task) {
      unique(sb_task_group_ids(task))
    }), use.names = FALSE)
    if (anyDuplicated(groups)) {
      sb_abort(
        "Group IDs cannot occur in both development and held-out cohorts.",
        "stormbreaker_internal_external_patient_overlap"
      )
    }
  }

  patient_ids <- lapply(tasks, function(task) {
    as.character(task$data[[task$id_col]])
  })
  all_ids <- unlist(patient_ids, use.names = FALSE)
  repeated <- unique(all_ids[duplicated(all_ids)])
  if (length(repeated) > 0L) {
    examples <- paste(utils::head(repeated, 3L), collapse = ", ")
    sb_abort(
      paste0(
        "Patient IDs must be globally unique across cohorts so a patient ",
        "cannot enter both development and held-out data. Repeated IDs: ",
        examples, ". Namespace cohort-local IDs before creating the tasks."
      ),
      "stormbreaker_internal_external_patient_overlap"
    )
  }

  list(
    tasks = tasks,
    task_type = task_type,
    cohort_ids = cohort_ids,
    patient_ids = patient_ids,
    predictors = reference_predictors
  )
}

sb_validate_iecv_pipelines <- function(pipelines, task_type) {
  valid_collection <- is.list(pipelines) && length(pipelines) > 0L &&
    !inherits(pipelines, "stormbreaker_pipeline")
  if (!valid_collection) {
    sb_abort(
      paste0(
        "`pipelines` must be a non-empty list of explicit specifications ",
        "created by `sb_pipeline()`."
      ),
      "stormbreaker_bad_internal_external_pipelines"
    )
  }
  pipelines <- lapply(pipelines, function(pipeline) {
    sb_validate_pipeline_contract(pipeline, task_type)
  })
  pipeline_ids <- vapply(pipelines, `[[`, character(1L), "id")
  if (anyDuplicated(pipeline_ids)) {
    sb_abort(
      "Every internal-external validation pipeline must have a distinct ID.",
      "stormbreaker_duplicate_pipeline_id"
    )
  }
  names(pipelines) <- pipeline_ids
  pipelines
}

sb_validate_iecv_strategy <- function(pipelines, task_type) {
  if (inherits(pipelines, "stormbreaker_stack")) {
    stack <- sb_validate_stack_contract(pipelines, task_type)
    return(list(
      type = "stack",
      specification = stack,
      ids = stack$id
    ))
  }
  validated <- sb_validate_iecv_pipelines(pipelines, task_type)
  list(
    type = "pipelines",
    specification = validated,
    ids = names(validated)
  )
}

sb_fit_iecv_strategy <- function(task, strategy) {
  if (identical(strategy$type, "stack")) {
    return(stats::setNames(
      list(sb_fit_stack_final(task, strategy$specification)),
      strategy$ids
    ))
  }
  sb_with_selector_reuse(lapply(strategy$specification, function(pipeline) {
    sb_fit_final(task, pipeline)
  }))
}

sb_predict_iecv_strategy <- function(fits, strategy, new_data, id_col) {
  if (identical(strategy$type, "stack")) {
    return(stats::setNames(
      list(sb_predict_stack_final(fits[[1L]], new_data, id_col)),
      strategy$ids
    ))
  }
  lapply(fits, function(fit) {
    sb_predict_final(fit, new_data, id_col)
  })
}

sb_iecv_task_rows <- function(task, reference, task_type) {
  role_names <- if (identical(task_type, "binary")) {
    c(task$id_col, task$outcome_col)
  } else {
    c(task$id_col, task$time_col, task$event_col)
  }
  reference_roles <- if (identical(task_type, "binary")) {
    c(reference$id_col, reference$outcome_col)
  } else {
    c(reference$id_col, reference$time_col, reference$event_col)
  }
  if (!is.null(reference$group_col)) {
    role_names <- c(role_names, task$group_col)
    reference_roles <- c(reference_roles, reference$group_col)
  }
  rows <- task$data[c(role_names, reference$predictors)]
  names(rows)[seq_along(reference_roles)] <- reference_roles
  rows
}

sb_iecv_development_task <- function(tasks, reference, task_type) {
  rows <- lapply(tasks, function(task) {
    sb_iecv_task_rows(task, reference, task_type)
  })
  data <- do.call(rbind, rows)
  rownames(data) <- NULL

  decision_thresholds <- if ("decision_thresholds" %in% names(reference)) {
    reference$decision_thresholds
  } else {
    NULL
  }
  if (identical(task_type, "binary")) {
    return(sb_binary_task(
      data,
      reference$id_col,
      reference$outcome_col,
      reference$predictors,
      decision_thresholds = decision_thresholds,
      group_col = reference$group_col
    ))
  }

  prediction_horizon <- if ("prediction_horizon" %in% names(reference)) {
    reference$prediction_horizon
  } else {
    NULL
  }
  sb_surv_task(
    data,
    reference$id_col,
    reference$time_col,
    reference$event_col,
    reference$predictors,
    reference$time_unit,
    prediction_horizon = prediction_horizon,
    decision_thresholds = decision_thresholds,
    group_col = reference$group_col
  )
}

sb_iecv_endpoint <- function(task, task_type) {
  if (identical(task_type, "binary")) {
    return(sb_binary_endpoint(task$data, task$id_col, task$outcome_col))
  }
  horizon <- if ("prediction_horizon" %in% names(task)) {
    task$prediction_horizon
  } else {
    NULL
  }
  sb_surv_endpoint(
    task$data,
    task$id_col,
    task$time_col,
    task$event_col,
    task$time_unit,
    prediction_horizon = horizon
  )
}

sb_iecv_cohort_table <- function(tasks, task_type) {
  cohort_ids <- names(tasks)
  table <- data.frame(
    cohort_id = cohort_ids,
    n_patients = vapply(tasks, `[[`, integer(1L), "n"),
    stringsAsFactors = FALSE
  )
  if (identical(task_type, "binary")) {
    table$n_positive <- vapply(tasks, `[[`, integer(1L), "n_positive")
    table$n_negative <- vapply(tasks, `[[`, integer(1L), "n_negative")
    return(table)
  }
  table$n_events <- vapply(tasks, `[[`, integer(1L), "n_events")
  table$n_censored <- vapply(tasks, `[[`, integer(1L), "n_censored")
  table$time_unit <- vapply(tasks, `[[`, character(1L), "time_unit")
  if ("prediction_horizon" %in% names(tasks[[1L]])) {
    table$prediction_horizon <- vapply(
      tasks, `[[`, numeric(1L), "prediction_horizon"
    )
  }
  table
}

#' @rdname sb_compare_cv
#' @export
sb_run_internal_external_cv <- function(cohorts, pipelines) {
  cohort_contract <- sb_validate_iecv_cohorts(cohorts)
  tasks <- cohort_contract$tasks
  task_type <- cohort_contract$task_type
  cohort_ids <- cohort_contract$cohort_ids
  reference <- tasks[[1L]]
  strategy <- sb_validate_iecv_strategy(pipelines, task_type)
  pipeline_ids <- strategy$ids

  folds <- lapply(seq_along(tasks), function(heldout_index) {
    heldout_id <- cohort_ids[[heldout_index]]
    development_ids <- cohort_ids[-heldout_index]
    development_task <- sb_iecv_development_task(
      tasks[-heldout_index], reference, task_type
    )
    fits <- sb_fit_iecv_strategy(development_task, strategy)

    heldout <- tasks[[heldout_index]]
    predictor_data <- heldout$data[c(
      heldout$id_col, cohort_contract$predictors
    )]
    predictions <- sb_predict_iecv_strategy(
      fits,
      strategy,
      predictor_data,
      heldout$id_col
    )

    endpoint <- sb_iecv_endpoint(heldout, task_type)
    comparison <- sb_compare_external(predictions, endpoint)
    list(
      status = comparison$status,
      heldout_cohort_id = heldout_id,
      development_cohort_ids = development_ids,
      development_patient_ids = as.character(
        development_task$data[[development_task$id_col]]
      ),
      heldout_patient_ids = cohort_contract$patient_ids[[heldout_index]],
      fits = fits,
      predictions = predictions,
      comparison = comparison
    )
  })
  names(folds) <- cohort_ids

  results <- lapply(folds, function(fold) {
    table <- fold$comparison$comparison
    cbind(
      data.frame(
        heldout_cohort_id = rep(fold$heldout_cohort_id, nrow(table)),
        n_development_cohorts = rep(
          as.integer(length(fold$development_cohort_ids)), nrow(table)
        ),
        n_development_patients = rep(
          as.integer(length(fold$development_patient_ids)), nrow(table)
        ),
        n_heldout_patients = rep(
          as.integer(length(fold$heldout_patient_ids)), nrow(table)
        ),
        stringsAsFactors = FALSE
      ),
      table
    )
  })
  results <- do.call(rbind, results)
  rownames(results) <- NULL

  complete <- all(vapply(
    folds, function(fold) identical(fold$status, "success"), logical(1L)
  ))
  result <- list(
    status = if (complete) "success" else "incomplete",
    task_type = task_type,
    validation_scheme = "leave_one_cohort_out_internal_external_cv",
    design = list(
      heldout_unit = "cohort",
      development_rule = "all_other_cohorts",
      development_strategy =
        "pooled_remaining_cohort_patients_without_cohort_effects",
      heldout_outcomes = "evaluation_only",
      fold_training_sets_overlap = length(tasks) > 2L,
      final_independent_external_validation = FALSE
    ),
    pipeline_ids = pipeline_ids,
    cohorts = sb_iecv_cohort_table(tasks, task_type),
    results = results,
    folds = folds
  )
  if (identical(strategy$type, "stack")) {
    result$strategy_type <- "stack"
  }
  structure(result, class = "stormbreaker_internal_external_cv")
}

#' @export
print.stormbreaker_internal_external_cv <- function(x, ...) {
  cat("<stormbreaker_internal_external_cv:", x$status, ">\n")
  cat("leave-one-cohort-out internal-external cross-validation\n")
  count_label <- if (identical(x$strategy_type, "stack")) {
    "Stacks:"
  } else {
    "Pipelines:"
  }
  cat(
    "Cohorts:", nrow(x$cohorts),
    "|", count_label, length(x$pipeline_ids),
    "| Task:", x$task_type, "\n"
  )
  metric_columns <- if (identical(x$task_type, "binary")) {
    c("roc_auc", "brier_score")
  } else {
    c("harrell_c", "ipcw_brier", "time_auc", "uno_c")
  }
  columns <- c(
    "heldout_cohort_id", "pipeline_id", "status",
    intersect(metric_columns, names(x$results))
  )
  print(x$results[columns], row.names = FALSE, digits = 3)
  cat(
    paste0(
      "Each held-out cohort uses a different fitted model; no metric is ",
      "pooled, ranked, or used to select a winner.\n"
    )
  )
  cat(
    paste0(
      "This is method-development IECV, not final independent external ",
      "validation after pipeline choices are revised.\n"
    )
  )
  invisible(x)
}
