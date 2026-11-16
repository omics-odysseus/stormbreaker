#' Model menus
#' @name sb_models
#' @description Inspect available models and construct task-specific pipeline menus.
#' @rawRd \arguments{
#' \item{task}{
#' One of binary, survival or multiclass.
#' }
#' \item{task_type}{
#' \strong{For \code{sb_model_catalog()}, \code{sb_model_pipelines()}:} Explicitly \code{"binary"} or \code{"survival"}.
#' 
#' \strong{For \code{sb_pipeline_registry()}:} Optional task filter: \code{"survival"} or \code{"binary"}.
#' 
#' \strong{For \code{sb_candidate_registry()}:} Optional task filter: \code{NULL}, \code{"survival"}, or \code{"binary"}.
#' }
#' \item{models}{
#' Character model IDs from \code{\link[=sb_model_catalog]{sb_model_catalog()}}, in desired
#' order. \code{NULL} means all catalog entries for this task, including entries
#' whose optional dependency is missing. No winner is selected.
#' 
#' \strong{For \code{sb_binary_workflow_plan()}:} Binary menu names or model IDs, or
#' \code{NULL} to inspect all binary menu entries, including selector-only methods.
#' }
#' \item{selector}{
#' For \code{sb_binary_workflow_plan()}, one registered binary selector ID (or
#' short selector name from the binary menu), or \code{"none"}. Defaults to RF
#' importance screening. All compatible predictors use the same selector.
#' }
#' \item{selector_options}{
#' Options for the requested selector, as in \code{sb_pipeline()}. With the RF
#' selector, NULL means \code{list(top_n = 5L, spearman_cutoff = .70)} in this
#' planner; use an empty list for the original importance-threshold settings.
#' }
#' \item{include_unfiltered}{
#' Explicitly include registered \code{"none"} routes when the requested
#' selector-learner pair is not admitted. They are labeled in the plan and
#' do not inherit selector options. Missing packages and failed fits do not
#' trigger a replacement route.
#' }
#' \item{inner_folds}{
#' \strong{For \code{sb_model_pipelines()}:} Inner CV settings for glmnet configurations.
#' Other learners retain their documented \code{sb_pipeline()} defaults.
#' 
#' \strong{For \code{sb_binary_menu_pipelines()}:} Number of inner folds for the glmnet logistic learner,
#' at least 3. Defaults to 5. It does not override selector-native tuning.
#' 
#' \strong{For \code{sb_binary_workflow_plan()}:} Number of inner glmnet learner
#' folds, at least 3. Selector-side tuning is unchanged.
#' 
#' \strong{For \code{sb_mime117_pipelines()}:} Number of inner folds used by downstream Cox glmnet
#' learners.
#' }
#' \item{inner_seed}{
#' \strong{For \code{sb_model_pipelines()}:} Inner CV settings for glmnet configurations.
#' Other learners retain their documented \code{sb_pipeline()} defaults.
#' 
#' \strong{For \code{sb_binary_menu_pipelines()}:} Positive integer seed for glmnet learner tuning. Other
#' learners and selectors retain their saved native defaults.
#' 
#' \strong{For \code{sb_binary_workflow_plan()}:} Seed for inner glmnet learner tuning.
#' 
#' \strong{For \code{sb_mime117_pipelines()}:} Positive integer seed for downstream Cox glmnet folds.
#' }
#' \item{selector_features}{
#' \strong{For \code{sb_binary_menu_pipelines()}:} Ordered, unique names for fixed-feature routes.
#' 
#' \strong{For \code{sb_binary_workflow_plan()}:} Ordered, unique feature names when
#' \code{selector = "fixed_features"}; otherwise NULL.
#' 
#' \strong{For \code{sb_survival_menu_pipelines()}:} Ordered, unique names for the caller-declared
#' fixed-feature routes. Every fixed-feature pipeline uses this same vector.
#' }
#' \item{feature_set_id}{
#' Non-empty label for the declared feature block.
#' }
#' \item{selector_n_features}{
#' \strong{For \code{sb_binary_menu_pipelines()}:} Positive integer, at least 2, for Top-K screening.
#' Defaults to 10 and is not selected from outcomes or clipped to the data.
#' 
#' \strong{For \code{sb_binary_workflow_plan()}:} Positive integer, at least 2,
#' where the registered selector requires a Top-K size (defaults to 5).
#' RF Top-K is instead controlled by \code{selector_options}.
#' 
#' \strong{For \code{sb_survival_menu_pipelines()}:} Positive integer, at least 2, for Top-K selectors.
#' The default is 30. This number must be declared before evaluation; the
#' constructor does not choose it from outcomes or clip it to observed data.
#' }
#' \item{level}{
#' Return one row per executable \code{"route"} or one row per
#' \code{"learner"}.
#' }
#' \item{check_installed}{
#' Whether to add a side-effect-free check of installed
#' R packages. This uses package locations only: it does not load optional R
#' namespaces, fit a model, or claim that a method suits the supplied data.
#' }
#' \item{route_status}{
#' Optional exact runtime-route status returned by the
#' function. \code{NULL} retains every status.
#' }
#' \item{selector_inner_folds}{
#' Number of inner folds used by CoxBoost and Lasso
#' first-stage selectors.
#' }
#' \item{selector_seed}{
#' Positive integer seed used by learned first-stage
#' selectors.
#' }
#' \item{rsf_nodesize}{
#' Positive integer terminal-node size used by RSF
#' first-stage selection.
#' }
#' \item{learner_tuning_folds}{
#' Number of internal folds for the reference
#' CoxBoost, SuperPC, and GBM learner fits.
#' }
#' \item{plsrcox_tuning_folds}{
#' Number of internal folds for the reference
#' plsRcox learner fit. Mime1 uses the plsRcox default of 5.
#' }
#' \item{learner_seed}{
#' Positive integer seed used by analysis-fold tuned
#' downstream learners.
#' }
#' \item{gbm_cv_max_trees}{
#' Maximum number of trees searched by reference GBM.
#' A value of one fits the sole candidate without an inner tree-count search
#' or CV curve; outer performance evaluation is unchanged.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_models()}}{
#' List short model names for doctor-facing entrypoints
#' 
#' \strong{Returns:} A model menu. Short names map to existing native configurations.
#' Binary rows declare can_predict, can_extract_features and stage_selector.
#' }
#' 
#' \subsection{\code{sb_model_catalog()}}{
#' The menu contains every admitted learner with selector \code{none}, plus explicit
#' Ridge and alpha-0.5 elastic-net configurations for the glmnet learner.
#' Existing glmnet learner IDs denote Lasso (alpha 1). Model variants are not
#' an exhaustive hyperparameter grid. Dependencies are reported, not installed.
#' 
#' \strong{Returns:} A data frame with model IDs, learner IDs, glmnet alpha, registry
#' semantics and dependency availability. Missing dependencies stay visible.
#' }
#' 
#' \subsection{\code{sb_model_pipelines()}}{
#' Every pipeline uses selector \code{none}: the task's user-supplied predictors
#' enter fitting without an added screening stage. Intrinsic sparsity, subset
#' selection, component estimation and tree splits remain part of the selected
#' learner; this does not force every coefficient to be nonzero.
#' 
#' \strong{Returns:} A named list of \code{\link[=sb_pipeline]{sb_pipeline()}} specifications, ready for
#' \code{\link[=sb_compare_cv]{sb_compare_cv()}} or \code{\link[=sb_fit_final]{sb_fit_final()}}.
#' }
#' 
#' \subsection{\code{sb_binary_menu_pipelines()}}{
#' Creates one specification for every registered binary selector-learner pair.
#' The constructor does not read expression data, fit models, rank methods or
#' choose genes. All fixed-feature routes use the same caller-declared vector.
#' Learned selectors retain their existing native defaults, including their
#' own tuning folds, seeds and Boruta budgets. Top-K size and glmnet learner
#' tuning controls are explicitly recorded in the returned pipelines.
#' 
#' \strong{Returns:} A data frame with unique \code{config_id}, \code{selector}, \code{learner},
#' \code{mime_family}, \code{comparison_role} and a \code{pipeline} list-column. The
#' \code{mime7_family_mapping} attribute records the seven family references and
#' their non-equivalence to MIME tuning. \code{registry_pairs} records coverage.
#' Counts describe configurations, not independent algorithm families.
#' 
#' \strong{Details:} The seven MIME classification families are identified for comparison, not
#' reproduced by this constructor. These fixed STORMBREAKER configurations do
#' not implement MIME's caret tuning grids or repeated cross-validation; run
#' the original MIME methods separately for that benchmark. The mapping is
#' deliberately labelled \code{algorithm_family_only}.
#' 
#' The \code{none} and \code{fixed_features} routes remain separate specifications.
#' They may use identical inputs when the task's predictors equal the fixed
#' feature vector. Use the shared comparison/reuse interfaces rather than
#' treating such duplicate inputs as additional independent evidence.
#' }
#' 
#' \subsection{\code{sb_binary_workflow_plan()}}{
#' Prepare one route per requested prediction configuration under a common
#' upstream selector. Lasso, Ridge and elastic net retain distinct alpha
#' values. This inspects the registry and local R package locations but does
#' not fit, rank, draw, export or install anything. Selector-only methods
#' are listed separately, not treated as classifiers. Fixed-feature routes
#' additionally need caller-declared features and a feature-set ID.
#' 
#' \strong{Returns:} A list with \code{pipelines} (named specifications for
#' \code{bi_run(..., pipelines = plan$pipelines, plots = FALSE)}) and
#' \code{candidates} (model names, IDs, effective selector, availability,
#' plan status, reasons and constructed pipeline IDs). A ready specification
#' still needs evaluation; data-specific failures are recorded by the run.
#' }
#' 
#' \subsection{\code{sb_survival_menu_pipelines()}}{
#' Creates one specification for every registered survival selector-learner
#' pair and includes all 117 MIME reference variants. Identical complete
#' pipeline IDs are deduplicated. This constructor only creates specifications;
#' it does not read expression data, fit models, rank methods or select genes.
#' 
#' \strong{Returns:} A data frame with a \code{pipeline} list-column, unique \code{config_id},
#' \code{selector}, \code{learner}, \code{mime_labels} and \code{comparison_role} columns.
#' The \code{mime117_mapping} attribute maps the 117 reference labels to menu
#' configuration IDs. The \code{registry_pairs} attribute records registry coverage.
#' Menu size can increase when the package admits additional methods; counts
#' describe configurations, not independent algorithm families.
#' 
#' \strong{Details:} Common MIME settings are saved in the returned pipeline objects: 10 inner
#' folds for glmnet, learned Lasso/CoxBoost selectors, CoxBoost, GBM and SuperPC;
#' 5 folds for PLS-Cox; a common declared seed of 1; and up to 10000 GBM trees.
#' The reference suite retains all Elastic Net alphas and StepCox directions.
#' Other registered methods keep their native package defaults, including the
#' default Boruta budget. Fixed-component and fixed-step extensions remain
#' separate from MIME's tuned learners. The low-level \code{\link[=sb_pipeline]{sb_pipeline()}} API keeps
#' its explicit parameter requirements and is not changed by this preset.
#' 
#' Parameter alignment does not imply identical internal fold assignments or
#' bug-for-bug MIME implementation. It does not turn a training-selected input
#' panel into independent or unbiased validation. See \code{\link[=sb_mime117_pipelines]{sb_mime117_pipelines()}}
#' for the reference suite's implementation boundaries.
#' }
#' 
#' \subsection{\code{sb_pipeline_registry()}}{
#' Returns the selector and learner combinations that the installed package
#' knows how to construct and resolve. \code{admitted} means that STORMBREAKER has a
#' matching built-in specification and adapter route. It does not mean that an
#' optional package is installed, that a model will fit a particular dataset,
#' or that the statistical design is appropriate. \code{source} identifies package
#' code, not literature provenance. At learner level, dependency columns
#' describe the learner itself (the \code{selector = "none"} route); use route level
#' when dependencies of a learner plus a feature selector must both be checked.
#' 
#' \strong{Returns:} A data frame. At route level, each row is an admitted selector and
#' learner combination. At learner level, each row summarizes one learner,
#' its available selectors, pipeline count, tuning semantics, and optional
#' dependency. With \code{check_installed = TRUE}, four additional columns list R
#' packages, a reserved module field, missing R packages, and one dependency
#' status. The current runtime is R-only, so \code{python_modules} is empty for
#' every registered route. Admission does not imply that optional software is
#' installed or operational.
#' }
#' 
#' \subsection{\code{sb_candidate_registry()}}{
#' \code{sb_candidate_registry()} returns the bundled task-aware catalog of 1,175 retained
#' historical candidate specifications and appends a descriptive mapping to
#' STORMBREAKER's current runtime. The source columns are preserved; in
#' particular, historical adapter, execution, and admission fields are not
#' upgraded because a similarly named runtime route now exists.
#' 
#' \strong{Returns:} A data frame retaining every source catalog column and adding the
#' current runtime task, mapped selector and learner, parameters declared by
#' the label, route status and availability, configuration status, source-
#' specification audit status/reference/difference, missing components,
#' current non-admission status, and mapping reason.
#' 
#' \strong{Details:} \code{runtime_route_status = "route_available"} means only that the current
#' package has a task-matched selector and learner family route capable of
#' representing the catalog shape. It does not verify that the source used the
#' same operational definition, package, defaults, tuning procedure, or
#' parameters. Every currently mapped row therefore remains non-admitted. The
#' appended \code{runtime_source_spec_status} separates routes that have not yet
#' been audited from audits that found a name-level ambiguity or operational
#' semantics that differ; it never promotes a family match to an executable
#' candidate. The function never constructs, fits, ranks, tunes, or selects a
#' pipeline.
#' 
#' The first source audit covers available routes that contain the StepCox
#' selector. The historical catalog merges the label across sources with
#' different operational definitions: official Mime1 and the historical Fig.5A
#' reconstruction differ in their common prefilter, forward start, minimum-
#' feature, warning, truncation, and fallback behavior. STORMBREAKER's current
#' fold-local AIC selector is therefore not silently called an exact historical
#' implementation. These rows are marked
#' \code{audited_source_identity_ambiguous} and remain
#' \code{not_admitted_source_identity_ambiguous}. This records provenance; it does
#' not change the admitted StepCox runtime.
#' 
#' The second source audit uses official Mime1 only as an operational comparison
#' for the direct binary \code{NaiveBayes}, \code{AdaBoost}, and \code{LogitBoost} labels. Mime1
#' tunes each method by 5-fold, 10-repeat ROC. Its Naive Bayes route uses a
#' kernel model, its AdaBoost grid includes both AdaBoost.M1 and Real AdaBoost,
#' and its LogitBoost route searches iteration counts. STORMBREAKER instead
#' exposes fixed Gaussian Naive Bayes, fixed AdaBoost.M1, and fixed LogitBoost
#' methods without those searches. Because the catalog rows provide no source
#' DOI or figure and identify themselves as union-generated, this comparison
#' does not attribute them to Mime1. They are marked
#' \code{audited_source_identity_unverified} and remain
#' \code{not_admitted_source_identity_unverified}. The three current STORMBREAKER
#' methods remain available under their own explicit specifications.
#' 
#' The third source audit covers the 11 direct survival \code{Ridge}, \code{Lasso}, and
#' elastic-net labels. Their catalog evidence is only
#' \code{historical_fig5a_label_inferred}, and the two available historical
#' implementations differ: official Mime1 uses fixed 10-fold cross-validation,
#' while the Fig.5A reconstruction uses a sample-size-dependent fold count.
#' Mime1's Ridge code also fits \code{alpha = 0} but obtains \code{lambda.min} from a
#' separate cross-validation call that omits \code{alpha}; the reconstruction
#' cross-validates \code{alpha = 0}. STORMBREAKER requires explicit alpha, inner
#' folds, seed, measure, lambda rule, and standardization inside outer analysis
#' resampling. The rows are therefore marked
#' \code{audited_source_identity_ambiguous}, remain non-admitted, and do not change
#' the current glmnet Cox route.
#' 
#' The fourth source audit covers the 12 direct binary \code{LogisticRegression},
#' \code{Ridge}, \code{Lasso}, and elastic-net labels. The catalog generator created them
#' as a label union with \code{union_generated_candidate_endpoint_unverified}; its
#' provenance marks these families outside the Mime1 binary method set, and the
#' rows provide no DOI or figure. The labels do not specify ordinary-logistic
#' convergence semantics or glmnet lambda selection, folds, seed, objective,
#' rule, standardization, preprocessing, or failures. They are marked
#' \code{audited_source_identity_unverified} and remain non-admitted. The current
#' explicit logistic and nested glmnet-logistic routes are unchanged.
#' 
#' The fifth source audit covers nine direct fixed survival labels.
#' \code{survival-SVM}, \code{CoxBoost}, \code{GBM}, and \code{SuperPC} carry only
#' \code{historical_fig5a_label_inferred}; official Mime1 and the historical Fig.5A
#' reconstruction differ in tuning, retry, preprocessing, or score semantics,
#' so these four rows are \code{audited_source_identity_ambiguous}. Direct \code{plsRcox}
#' is a 101/117 compatibility addition outside the reconstructed 296
#' single-model set, and \code{Ranger} is supported only by a Table 1 family/package
#' label while Mime1 uses RSF/randomForestSRC. Those two rows are
#' \code{audited_source_identity_unverified}. All six remain non-admitted. The
#' current fixed learners keep their own explicit parameters and are not
#' presented as exact reproductions of the catalog rows. The same audit records
#' that official Mime1 contains no GamBoost, GlmBoost, or BlackBoost survival
#' learner; local reconstructions provide comparison implementations, including
#' one incorrect GamBoost-to-glmboost dispatch. STORMBREAKER therefore exposes
#' three distinct native mboost algorithms while keeping their catalog rows
#' source-unverified.
#' 
#' Seven direct binary labels (\code{LDA}, \code{QDA}, \code{KNN}, \code{DecisionTree}, \code{RF}, \code{SVM},
#' and \code{GBM}) currently have family-level runtime routes but have not yet
#' received source-to-specification audits. They therefore carry
#' \code{runtime_source_spec_status = "not_yet_audited"} together with
#' \code{runtime_admission_status = "not_admitted_source_to_spec_unverified"}.
#' Route availability for these labels is not a claim of exact historical
#' reproduction.
#' 
#' Fixed runtime implementations now also create descriptive family matches for
#' binary and survival \code{XGBoost}, \code{RSF}, \code{Rpart}, \code{Ctree}, \code{CForest}, and the
#' compatible \code{SurvReg} selector routes. Direct matches without a completed
#' source audit remain \code{not_yet_audited}; StepCox-containing matches retain the
#' StepCox source-identity decision. No such row is promoted to an executable
#' historical specification.
#' 
#' Unresolved multi-stage labels remain unresolved rather than having roles
#' inferred from their order. Methodologically blocked rows remain blocked.
#' Names such as RSF and RF are not silently equated with ranger, and binary GBM
#' is not equated with LightGBM or CatBoost.
#' }
#' 
#' \subsection{\code{sb_mime117_pipelines()}}{
#' Reconstructs the 117 method configurations generated by the published
#' Mime1 loop structure as real STORMBREAKER pipeline objects in the original
#' source order. Shared selectors and learners are implemented once and reused
#' across configurations. Calling this function constructs and validates
#' specifications only: it does not fit models, inspect outcomes, rank routes,
#' or choose a winner.
#' 
#' \strong{Returns:} A 117-row data frame in reference source order. The \code{pipeline}
#' list-column contains validated \code{stormbreaker_pipeline} objects; the other
#' columns describe the reference method label and its two stage roles.
#' 
#' \strong{Details:} Learned selectors are fold-local when the returned pipelines are evaluated.
#' The suite preserves the method structure but intentionally does not copy
#' same-cohort apparent evaluation, unbounded retries, parallel defaults, or a
#' known inconsistent Ridge-alpha path from the reference implementation.
#' }
#' }
#' @rawRd \examples{
#' # sb_model_catalog
#' catalog <- sb_model_catalog("binary")
#' head(catalog[c("model_id", "r_packages", "dependency_status")])
#' 
#' # sb_model_pipelines
#' sb_model_pipelines("binary", c("logistic", "gaussian_naive_bayes"))
#' sb_model_pipelines("survival", c("coxph", "survreg_weibull"))
#' 
#' # sb_pipeline_registry
#' registry <- sb_pipeline_registry()
#' registry[c("task_type", "selector", "learner")]
#' sb_pipeline_registry(task_type = "binary", level = "learner")
#' sb_pipeline_registry(level = "learner", check_installed = TRUE)
#' 
#' # sb_candidate_registry
#' candidates <- sb_candidate_registry()
#' table(candidates$runtime_task_type, candidates$runtime_route_status)
#' sb_candidate_registry(
#'   task_type = "survival",
#'   route_status = "route_available"
#' )
#' }
NULL

#' Complete modeling workflows
#' @name sb_complete_run
#' @description Compare models, select a pipeline, complete its evaluation and inspect recorded issues.
#' @rawRd \arguments{
#' \item{run}{
#' Result from \code{\link[=bi_run]{bi_run()}} or \code{\link[=surv_run]{surv_run()}}.
#' }
#' \item{pipeline}{
#' Exact name from names(run$pipelines), chosen using development
#' evidence, not external outcomes. No winner is selected automatically.
#' }
#' \item{external}{
#' Named list of independent cohort data frames. Use the same
#' ID, predictor and endpoint column names as development data.
#' }
#' \item{plots}{
#' Include available algorithm and external ROC plots.
#' }
#' \item{capture_diagnostics}{
#' Retain available native diagnostic data in a new final fit. By default
#' this is enabled when plotting or explaining. Set TRUE with
#' \code{plots = FALSE} to defer plotting until after the run. Any already
#' saved fit is reused, not refitted.
#' }
#' \item{explain}{
#' Optional external cohort name to explain. Development data
#' supply the SHAP background; external outcomes never enter explanation.
#' }
#' \item{reverse}{
#' \code{FALSE} (default) keeps the development-defined direction.
#' For binary outcomes, explicit \code{TRUE} uses cohort outcomes to flip
#' probabilities when AUC < 0.5. Survival outcomes require \code{FALSE}.
#' }
#' \item{task}{
#' A task from \code{\link[=sb_binary_task]{sb_binary_task()}} or \code{\link[=sb_surv_task]{sb_surv_task()}}, with the
#' user's gene set supplied as \code{predictors}.
#' }
#' \item{models}{
#' Defaults to \code{"recommended"}: logistic/Lasso/random forest for
#' binary tasks, or Cox/Lasso Cox/ranger for survival. These are fixed,
#' complementary starting models, not selected using observed performance.
#' Supply IDs from \code{\link[=sb_model_catalog]{sb_model_catalog()}} to choose models, or explicitly use
#' \code{NULL} for all models. Missing engines remain failures, not substitutions.
#' }
#' \item{folds}{
#' Shared outer resampling settings passed to \code{\link[=sb_compare_cv]{sb_compare_cv()}}.
#' }
#' \item{seed}{
#' Shared outer resampling settings passed to \code{\link[=sb_compare_cv]{sb_compare_cv()}}.
#' }
#' \item{inner_folds}{
#' Glmnet inner CV settings; see \code{\link[=sb_model_pipelines]{sb_model_pipelines()}}.
#' }
#' \item{inner_seed}{
#' Glmnet inner CV settings; see \code{\link[=sb_model_pipelines]{sb_model_pipelines()}}.
#' }
#' \item{view}{
#' Open the results matrix after comparison. Defaults to FALSE;
#' set TRUE to explicitly request a browser.
#' }
#' \item{report_file}{
#' Optional new HTML path; saves a page even with \code{view = FALSE}.
#' }
#' \item{x}{
#' Result of bi_run or surv_run.
#' }
#' \item{comparison}{
#' A result from \code{\link[=sb_compare_cohorts]{sb_compare_cohorts()}}. At least one cohort
#' must have the explicit role \code{"validation"}.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_complete_run()}}{
#' Continue a chosen pipeline without repeating model comparison
#' 
#' \strong{Returns:} A common bundle with comparison, fit, genes, predictions, metrics,
#' evaluations, plots and issues. Predictions contain patient identifiers.
#' \code{predictions} retains original model outputs; \code{evaluated_predictions}
#' matches the reported metrics and contains original probabilities and flip
#' flags for processed evaluations. \code{direction_correction} records each cohort.
#' External failures are reported individually; they never trigger refitting.
#' Survival risk_groups records counts at the locked development median even
#' when plots is FALSE. Single groups are reported in issues, not repaired by
#' choosing a new external threshold; prediction metrics remain available.
#' With survival plots, risk_groups also retains logrank_p, logrank_status
#' and logrank_error. Failed log-rank tests warn and appear in issues while
#' usable plots remain available. Without plots the test is not computed.
#' Printing shows a short summary without patient rows. Saved algorithm plots
#' are reused unless a new explanation is requested.
#' }
#' 
#' \subsection{\code{sb_compare_models()}}{
#' Users create a binary or survival task and explicitly supply \code{predictors}
#' to that task constructor. This wrapper does not infer the task from outcome
#' values, discover genes, add a selector, select a winning model or read
#' external data. It runs the chosen pipelines on the same outer folds.
#' Failed models remain in the underlying comparison rather than being replaced.
#' 
#' \strong{Returns:} A \code{stormbreaker_cv_comparison} with additional \code{model_pipelines},
#' \code{input_genes}, \code{requested_models} and total \code{elapsed_seconds}. Its table
#' also contains \code{model_id}. Total time includes menu construction and CV,
#' but not final fitting or external prediction.
#' }
#' 
#' \subsection{\code{sb_run_issues()}}{
#' Explain issues from a doctor-facing model run
#' 
#' \strong{Returns:} Data frame with model, stage, category, Chinese explanation and
#' original reason. Plot unavailability is separate from fitting failures.
#' Unknown failures remain unclassified; raw evidence is never replaced.
#' }
#' 
#' \subsection{\code{sb_select_model()}}{
#' Select the highest equally weighted mean validation ROC AUC (binary) or
#' Harrell C (survival). The returned object is an existing training-only fit,
#' ready for \code{\link[=sb_predict_final]{sb_predict_final()}} and \code{saveRDS()}. No model is refitted.
#' 
#' \strong{Returns:} A \code{stormbreaker_final_fit} with a \code{model_selection} audit containing
#' the validation cohort names, metric, per-cohort scores and full ranking.
#' 
#' \strong{Details:} A candidate is eligible only when its final fit and development CV
#' succeeded, and every validation cohort has a successful finite primary
#' metric. Ineligible candidates remain in the ranking. Ties preserve the
#' original pipeline declaration order. Test and external cohorts never
#' enter selection. Keep the final test cohort out of the comparison call
#' until the model is saved; select once and assess it on the held-out data.
#' Validation scores used to choose a model are not unbiased estimates of
#' the selected model's performance. The function does not choose a decision
#' threshold, recalibrate predictions or add validation patients to training.
#' Exactly constant predictions in any validation cohort, or in every
#' development assessment fold, are ineligible even when AUC/C is 0.5.
#' Variation is checked within folds, not across their different baselines.
#' Saved metrics and fitted states are not changed. The selection audit
#' retains prediction diagnostics and exclusion reasons. Missing prediction
#' evidence cannot establish eligibility; test/external predictions are not
#' inspected. This is not a near-constant tolerance or a performance cutoff.
#' }
#' }
#' @rawRd \examples{
#' # sb_complete_run
#' development <- bi_example(40, seed = 42)
#' validation <- bi_example(30, seed = 17)
#' run <- bi_run(development, outcome = "outcome", id = "patient_id",
#'   predictors = c("gene_a", "gene_b"), models = "Logistic", folds = 2,
#'   seed = 42, plots = FALSE)
#' finished <- finish_run(run, "Logistic",
#'   external = list(Validation = validation), plots = FALSE)
#' result_table(finished$evaluations$Validation, metric = "roc_auc")
#' fit <- finished$fit
#' prediction <- predict_final(fit, validation[c("patient_id", "gene_a", "gene_b")],
#'   "patient_id")
#' head(prediction$predictions)
#' 
#' # sb_compare_models
#' d <- sb_example_binary(n = 60, seed = 42)
#' task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#' result <- sb_compare_models(
#'   task, models = c("logistic", "gaussian_naive_bayes"), folds = 3, seed = 42
#' )
#' result$comparison[c("model_id", "status", "roc_auc")]
#' 
#' # sb_select_model
#' make_task <- function(prefix, seed) {
#'   d <- sb_example_binary(40, seed)
#'   d$patient_id <- paste0(prefix, d$patient_id)
#'   sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#' }
#' compared <- sb_compare_cohorts(
#'   make_task("train", 1), list(Validation = make_task("val", 2)),
#'   roles = c(Validation = "validation"),
#'   models = c("logistic", "gaussian_naive_bayes"), folds = 3, view = FALSE
#' )
#' model <- sb_select_model(compared)
#' sb_predict_final(model, sb_example_binary(20, 3), "patient_id")
#' }
NULL

#' Multiclass models
#' @name multi_run
#' @description Fit and compare native multiclass models, predict probabilities and evaluate their alignment.
#' @rawRd \arguments{
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
#' \item{models}{
#' Method names from \code{\link[=sb_models]{sb_models()}}. NULL requests the whole task
#' menu. For binary tasks this is one shared combination pool: prediction
#' models and feature methods are not chosen from separate menus.
#' }
#' \item{folds}{
#' \strong{For \code{multi_run()}:} Outer resampling controls.
#' 
#' \strong{For \code{sb_compare_multiclass()}:} Outer stratified/grouped fold count and random seed.
#' }
#' \item{seed}{
#' \strong{For \code{multi_run()}:} Outer resampling controls.
#' 
#' \strong{For \code{sb_compare_multiclass()}:} Outer stratified/grouped fold count and random seed.
#' }
#' \item{plots}{
#' NULL draws detailed figures for one method only; TRUE draws all
#' requested methods and FALSE skips detailed fits. For selected multiclass
#' methods after comparison, use \code{\link[=sb_analyze_methods]{sb_analyze_methods()}}.
#' }
#' \item{task}{
#' A \code{\link[=sb_multiclass_task]{sb_multiclass_task()}}.
#' }
#' \item{methods}{
#' Supported classifier names; defaults to Lasso, random forest
#' and radial SVM, a fixed small starting set rather than a data-selected winner. See
#' sb_method_catalog("multiclass"). Standalone selectors (including Boruta,
#' varSelRF and VSURF) are not classifiers and are not accepted here.
#' }
#' \item{inner_folds}{
#' Inner tuning settings, used only within each
#' outer training set. Glmnet standardization is fitted within resampling.
#' }
#' \item{inner_seed}{
#' Inner tuning settings, used only within each
#' outer training set. Glmnet standardization is fitted within resampling.
#' }
#' \item{fit}{
#' A saved fit at \code{analysis$results[[method]]$fit} from
#' \code{\link[=sb_analyze_methods]{sb_analyze_methods()}} with a multiclass task.
#' }
#' \item{newdata}{
#' Data frame containing all fitted predictor columns.
#' }
#' \item{truth}{
#' Observed labels, in the same row order as probability.
#' }
#' \item{probability}{
#' Numeric matrix with one named column per class.
#' }
#' \item{x}{
#' A result from \code{\link[=sb_multiclass_metrics]{sb_multiclass_metrics()}}.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{multi_run()}}{
#' Run existing multiclass learners by short name
#' 
#' \strong{Returns:} Existing multiclass CV comparison. Multiclass selector chains and
#' three-stage presets are not implemented by this entrypoint.
#' }
#' 
#' \subsection{\code{sb_compare_multiclass()}}{
#' Compare native multiclass models with shared outer folds
#' 
#' \strong{Returns:} A \code{stormbreaker_multiclass_comparison} with a comparison table,
#' shared fold assignment, row-aligned OOF probabilities, per-class metrics
#' and each outer fit's tuning/feature results. A failed fold fails that
#' method's overall evaluation; partial predictions are not scored as complete.
#' No winner, clinical action threshold, DCA or confusion plot is generated.
#' Candidate genes must be prespecified, not selected using held-out outcomes.
#' }
#' 
#' \subsection{\code{sb_predict_multiclass()}}{
#' Predict probabilities from a saved native multiclass fit
#' 
#' \strong{Returns:} A probability matrix with rows in input order and columns in the
#' declared class order. No label conversion, refitting or normalization occurs.
#' }
#' 
#' \subsection{\code{sb_multiclass_metrics()}}{
#' Evaluate aligned multiclass probabilities
#' 
#' \strong{Returns:} A list with overall and per-class metrics and a confusion table.
#' Balanced accuracy is unweighted mean class recall. Macro F1 is unweighted
#' mean of class F1; zero predicted counts yield precision NA and F1 zero if
#' the class is present. If any class is absent, macro recall/F1/AUC are NA
#' rather than silently excluding it. AUC is one-versus-rest, macro averaged.
#' Brier score is mean sum of squared class errors (range 0 to 2). Log loss
#' clips probabilities only for the logarithm at 1e-15. Ties in predicted
#' class use the first declared class. These are row-weighted point estimates.
#' }
#' 
#' \subsection{\code{sb_plot_multiclass_confusion()}}{
#' Explicitly plot a multiclass confusion matrix
#' 
#' \strong{Returns:} A ggplot with exactly one count per true/predicted class pair.
#' No model fitting occurs. This plot is never generated by default.
#' }
#' }
NULL

# ---- doctor-workflow ----

# ---- doctor-complete ----

#' @rdname sb_complete_run
#' @export
sb_complete_run <- function(run, pipeline, external = list(), plots = TRUE, explain = NULL,
                            reverse = FALSE, capture_diagnostics = isTRUE(plots) || !is.null(explain)) {
  if (!is.logical(reverse) || length(reverse) != 1L || is.na(reverse)) {
    sb_abort("`reverse` must be TRUE or FALSE.", "stormbreaker_bad_reverse")
  }
  if (missing(reverse) && identical(run$task_type, "survival")) reverse <- FALSE
  if (!is.logical(capture_diagnostics) || length(capture_diagnostics) != 1L || is.na(capture_diagnostics)) {
    sb_abort("capture_diagnostics must be TRUE or FALSE.", "stormbreaker_bad_visualization_argument")
  }
  if (!is.character(pipeline) || length(pipeline) != 1L || is.na(pipeline) ||
    !pipeline %in% names(run$pipelines)) {
    stop("Choose one exact name from names(run$pipelines).")
  }
  task <- run$development_task
  if (is.null(task)) stop("This older run does not retain development data; use a current bi_run/surv_run result.")
  if (!is.list(external) || (length(external) && (is.null(names(external)) ||
    anyNA(names(external)) || any(!nzchar(names(external))) || anyDuplicated(names(external))))) {
    stop("Supply external cohorts as a uniquely named list of data frames.")
  }
  if (!is.null(explain) && (length(explain) != 1L || is.na(explain) || !explain %in% names(external))) {
    stop("explain must name one supplied external cohort.")
  }
  fit <- run$details[[pipeline]]$fit
  if (is.null(fit)) {
    fit <- sb_fit_final(task, run$pipelines[[pipeline]], capture_diagnostics = capture_diagnostics)
  }
  result <- structure(
    list(
      pipeline = pipeline, task_type = run$task_type,
      comparison = run$comparison$comparison, fit = fit, genes = data.frame(),
      predictions = list(), evaluated_predictions = list(), direction_correction = list(),
      metrics = list(), evaluations = list(), plots = list(),
      explanation = NULL, risk_groups = list(), issues = list()
    ),
    class = c("stormbreaker_completed_run", "list")
  )
  if (!identical(fit$status, "success")) {
    result$issues$fit <- fit
    return(result)
  }
  state <- fit$selection$state
  stages <- if (identical(state$id, "ordered_chain")) state$stages else list(state)
  if (!length(stages) || is.null(stages[[1]]$selected_features)) {
    result$genes <- data.frame(stage = 0L, method = "Input features", gene = task$predictors)
  } else {
    result$genes <- do.call(rbind, lapply(seq_along(stages), function(i) {
      data.frame(stage = i, method = stages[[i]]$id, gene = stages[[i]]$selected_features)
    }))
  }
  if (isTRUE(plots) || !is.null(explain)) {
    if (isTRUE(plots)) result$plots$comparison <- run$plots
    figures <- if (is.null(explain)) run$details[[pipeline]]$figures else NULL
    if (is.null(figures)) {
      figures <- tryCatch(sb_plots(fit,
        new_data = if (!is.null(explain)) external[[explain]][c(task$id_col, task$predictors)] else NULL,
        id_col = task$id_col,
        background_data = if (!is.null(explain)) task$data[c(task$id_col, task$predictors)] else NULL
      ), error = identity)
    }
    if (inherits(figures, "error")) {
      result$issues$algorithm <- conditionMessage(figures)
    } else {
      if (isTRUE(plots)) result$plots$algorithm <- figures$plots
      result$explanation <- figures$shap
      result$issues$algorithm <- figures$unavailable
    }
  }
  reference <- if (run$task_type == "survival" && length(external)) {
    tryCatch(surv_risk_reference(fit, task), error = identity)
  } else {
    NULL
  }
  for (cohort in names(external)) {
    attempt <- tryCatch(
      {
        d <- external[[cohort]]
        prediction <- sb_predict_final(fit, d[c(task$id_col, task$predictors)], task$id_col)
        endpoint <- if (run$task_type == "binary") {
          sb_binary_endpoint(d, task$id_col, task$outcome_col)
        } else {
          sb_surv_endpoint(d, task$id_col, task$time_col, task$event_col, task$time_unit,
            prediction_horizon = task$prediction_horizon
          )
        }
        evaluation <- sb_evaluate_external(prediction, endpoint, reverse = reverse)
        if (identical(attr(d, "stormbreaker_assessment_role", exact = TRUE), "internal_holdout")) {
          evaluation$assessment_role <- "internal_holdout"
        }
        list(prediction = prediction, evaluation = evaluation)
      },
      error = identity
    )
    if (inherits(attempt, "error")) {
      result$issues$external[[cohort]] <- conditionMessage(attempt)
      next
    }
    result$predictions[[cohort]] <- attempt$prediction$predictions
    result$evaluations[[cohort]] <- attempt$evaluation
    result$evaluated_predictions[[cohort]] <- attempt$evaluation$scored_predictions
    result$direction_correction[[cohort]] <- attempt$evaluation$direction_correction
    result$metrics[[cohort]] <- attempt$evaluation$metrics
    if (!identical(attempt$evaluation$status, "success")) {
      result$issues$external[[cohort]] <- attempt$evaluation
    } else if (isTRUE(plots) && run$task_type == "binary") {
      figure <- tryCatch(sb_plot_roc(attempt$evaluation, cohort_label = cohort), error = identity)
      if (inherits(figure, "error")) {
        result$issues$plot[[cohort]] <- conditionMessage(figure)
      } else {
        result$plots$external[[cohort]] <- sb_compact_plot(figure)
      }
    } else if (run$task_type == "survival") {
      if (inherits(reference, "error")) {
        result$issues$risk[[cohort]] <- conditionMessage(reference)
        next
      }
      scores <- attempt$prediction$predictions$risk
      n_high <- sum(scores > reference$cutoff)
      n_low <- length(scores) - n_high
      status <- if (n_high > 0 && n_low > 0) "two_groups" else "single_group"
      result$risk_groups[[cohort]] <- data.frame(
        status = status,
        cutoff = reference$cutoff, n_low = n_low, n_high = n_high,
        logrank_p = NA_real_, logrank_status = "not_computed",
        logrank_error = NA_character_
      )
      if (status == "single_group") {
        result$issues$risk[[cohort]] <- sprintf(
          "single_group: locked development cutoff %.8g gives Low=%d, High=%d; high-versus-low HR and log-rank comparison are unavailable. The cutoff was not changed.",
          reference$cutoff, n_low, n_high
        )
      }
      if (!isTRUE(plots)) next
      risk <- tryCatch(
        {
          ext_task <- sb_surv_task(
            external[[cohort]], task$id_col, task$time_col,
            task$event_col, task$predictors, task$time_unit
          )
          surv_risk_plots(reference, ext_task)
        },
        error = identity
      )
      if (inherits(risk, "error")) {
        result$issues$plot[[cohort]] <- conditionMessage(risk)
      } else {
        result$plots$external[[cohort]] <- risk$plot
        result$risk_groups[[cohort]]$logrank_p <- risk$logrank_p
        result$risk_groups[[cohort]]$logrank_status <- risk$logrank_status
        result$risk_groups[[cohort]]$logrank_error <- risk$logrank_error
        if (identical(risk$logrank_status, "failed")) {
          result$issues$risk[[cohort]] <- paste0("Log-rank test failed: ", risk$logrank_error)
        }
      }
    }
  }
  result
}

# ---- doctor-entrypoints ----

#' @rdname sb_models
#' @export
sb_models <- function(task = c("binary", "survival", "multiclass")) {
  task <- match.arg(task)
  if (task == "multiclass") {
    menu <- sb_method_catalog(task)
    menu$name <- menu$method
    menu$name[menu$method == "random_forest"] <- "RF"
    menu$name[menu$method == "svm_radial"] <- "SVM"
    menu$name[menu$method == "elastic_net"] <- "ElasticNet"
    return(menu)
  }
  menu <- sb_model_catalog(task)
  menu$name <- gsub("_(binary|survival|fixed|logistic|cox)(?=_|$)", "", menu$model_id, perl = TRUE)
  aliases <- c(
    logistic = "Logistic", coxph = "Cox", glmnet_logistic = "Lasso", glmnet_cox = "Lasso",
    ridge_logistic = "Ridge", ridge_cox = "Ridge", elastic_net_logistic = "ElasticNet", elastic_net_cox = "ElasticNet",
    random_forest_binary = "RF", ranger_survival = "RF", gbm_binary = "GBM", gbm_cox = "GBM",
    gaussian_naive_bayes = "NaiveBayes", bernoulli_naive_bayes_binary = "BernoulliNB",
    nearest_centroid_binary = "NearestCentroid", svm_radial_binary = "SVM",
    survivalsvm = "SVM", xgboost_binary = "XGBoost", xgboost_cox = "XGBoost",
    coxboost_cv = "CoxBoost", coxboost = "CoxBoost_fixed", lda_binary = "LDA", qda_binary = "QDA",
    plsrcox = "PLSR", plsrcox_cv = "PLSR_cv"
  )
  hits <- match(menu$model_id, names(aliases))
  menu$name[!is.na(hits)] <- unname(aliases[hits[!is.na(hits)]])
  menu$name <- make.unique(menu$name, sep = "-")
  if (task == "binary") menu <- sb_binary_unified_menu(menu)
  menu
}

sb_doctor_model_ids <- function(models, menu) {
  if (is.null(models)) {
    return(menu$model_id)
  }
  if (!is.character(models) || !length(models) || anyNA(models)) stop("Choose model names from sb_models().")
  if ("plsrcox" %in% menu$model_id) {
    old <- c(plsrcox = "plsrcox", plsrcox_cv = "plsrcox_cv", mboost = "glmboost_cox")
    alias <- match(tolower(models), names(old))
    models[!is.na(alias)] <- unname(old[alias[!is.na(alias)]])
  }
  index <- match(tolower(models), tolower(menu$name))
  index[is.na(index)] <- match(sb_canonical_learner_id(models[is.na(index)]), menu$model_id)
  if (anyNA(index)) stop("Unknown model: ", paste(models[is.na(index)], collapse = ", "), ". See sb_models().")
  if (anyDuplicated(index)) stop("Choose each model once.")
  menu$model_id[index]
}

sb_doctor_pipeline <- function(model, selector, task) {
  if (task == "binary" && is.null(selector)) return(sb_binary_model_route(model))
  ids <- sb_doctor_model_ids(model, sb_models(task))
  if (length(ids) != 1L) stop("Choose one final learner for a pipeline.")
  base <- sb_model_pipelines(task, ids)[[1L]]
  sb_doctor_pipeline_from_base(base, selector, task)
}

sb_doctor_pipeline_from_base <- function(base, selector, task) {
  if (is.null(selector)) {
    return(base)
  }
  if (!inherits(selector, "stormbreaker_selector_chain_spec") || selector$task_type != task) {
    stop("Selector and model tasks must match.")
  }
  stages <- lapply(selector$ids, sb_doctor_stage, learner = base$learner, base = base)
  if (length(stages) == 1L) {
    return(stages[[1L]])
  }
  sb_new_chain_pipeline(base, stages)
}

#' @rdname sb_pipeline
#' @export
bi_pipeline <- function(model = "RF", selector = NULL) sb_doctor_pipeline(model, selector, "binary")

#' @rdname sb_pipeline
#' @export
surv_pipeline <- function(model = "Cox", selector = NULL) sb_doctor_pipeline(model, selector, "survival")

sb_doctor_pipeline_set <- function(task, models, stages, selector = NULL) {
  if (task == "binary") return(sb_binary_model_plan(models, stages, selector))
  menu <- sb_models(task)
  ids <- sb_doctor_model_ids(models, menu)
  names_chosen <- menu$name[match(ids, menu$model_id)]
  bases <- sb_model_pipelines(task, ids)
  names(bases) <- names_chosen
  if (!is.numeric(stages) || !length(stages) || anyNA(stages) || !all(stages %in% 1:3)) {
    stop("stages must contain 1, 2 or 3.")
  }
  pipelines <- list()
  skipped <- list()
  add <- function(model, selections = NULL) {
    label <- paste(c(if (!is.null(selections)) selections$methods, model), collapse = " + ")
    p <- tryCatch(sb_doctor_pipeline_from_base(bases[[model]], selections, task), error = identity)
    if (inherits(p, "error")) {
      skipped[[label]] <<- conditionMessage(p)
    } else {
      pipelines[[label]] <<- p
    }
  }
  if (!is.null(selector)) {
    for (m in names_chosen) add(m, selector)
  } else {
    if (1 %in% stages) for (m in names_chosen) add(m)
    if (2 %in% stages) {
      for (s in names(sb_doctor_selector_menu(task))) {
        for (m in names_chosen) add(m, sb_doctor_selector(s, task))
      }
    }
    if (3 %in% stages) {
      presets <- list(c("Boruta", "StepCox", "Cox"), c("Lasso", "StepCox", "CoxBoost"))
      for (p in presets) if (p[3L] %in% names_chosen) add(p[3L], sb_doctor_selector(p[1:2], task))
    }
  }
  list(pipelines = pipelines, unavailable = skipped)
}

sb_doctor_data_task <- function(data, task, id, predictors, outcome = NULL,
                                time = NULL, event = NULL, time_unit = "month",
                                prediction_horizon = NULL, decision_thresholds = NULL) {
  if (inherits(data, c("stormbreaker_binary_task", "stormbreaker_surv_task", "stormbreaker_multiclass_task"))) {
    return(data)
  }
  if (!is.data.frame(data)) stop("Supply a data frame or a task.")
  if (is.null(id)) {
    id <- ".sb_row_id"
    if (id %in% names(data)) stop("Specify the ID column explicitly.")
    data[[id]] <- paste0("row", seq_len(nrow(data)))
  }
  if (is.null(predictors)) predictors <- setdiff(names(data), c(id, outcome, time, event))
  switch(task,
    binary = sb_binary_task(data, id, outcome, predictors,
      decision_thresholds = decision_thresholds),
    survival = sb_surv_task(data, id, time, event, predictors, time_unit,
      prediction_horizon = prediction_horizon, decision_thresholds = decision_thresholds),
    multiclass = sb_multiclass_task(data, id, outcome, predictors)
  )
}

sb_doctor_run <- function(task, type, models, stages, selector, folds, seed, plots, pipelines = NULL) {
  if ((type == "binary" && !inherits(task, "stormbreaker_binary_task")) ||
    (type == "survival" && !inherits(task, "stormbreaker_surv_task"))) {
    stop("Task and entrypoint do not match.")
  }
  if (is.null(pipelines)) {
    plan <- sb_doctor_pipeline_set(type, models, stages, selector)
  } else {
    if (!is.null(models) || !is.null(selector)) stop("Use pipelines or models/selector, not both.")
    if (!is.list(pipelines) || !length(pipelines) || is.null(names(pipelines)) ||
      anyNA(names(pipelines)) || any(!nzchar(names(pipelines))) || anyDuplicated(names(pipelines))) {
      stop("Supply pipelines as a non-empty, uniquely named list.")
    }
    pipelines <- lapply(pipelines, sb_validate_pipeline_contract, task_type = type)
    plan <- list(pipelines = pipelines, unavailable = list())
  }
  if (type == "survival" && anyDuplicated(sb_task_group_ids(task))) {
    for (name in names(plan$pipelines)) {
      p <- plan$pipelines[[name]]
      stage_ids <- if (p$selector == "ordered_chain") {
        vapply(p$selector_params$stages, `[[`, character(1), "selector")
      } else {
        p$selector
      }
      if (p$learner %in% sb_analysis_tuned_survival_learner_ids() || "coxboost_nonzero" %in% stage_ids) {
        plan$unavailable[[name]] <- "Native survival CV does not accept patient groups; not run on repeated patients."
        plan$pipelines[[name]] <- NULL
      }
    }
  }
  if (!length(plan$pipelines)) {
    stop("No requested combination is supported. Choose models/selectors from the task menu.")
  }
  result <- sb_compare_cv(task, unname(plan$pipelines), folds = folds, seed = seed)
  overview <- list()
  plot_errors <- list()
  if (is.null(plots) || isTRUE(plots)) {
    display_models <- sb_clinician_models(sb_viz_leaves(result))
    for (name in c("performance", if (type == "binary") "roc")) {
      p <- tryCatch(if (name == "performance") {
        sb_plot_performance(
          result,
          if (type == "binary") "roc_auc" else "harrell_c"
        )
      } else {
        sb_plot_roc(result)
      }, error = identity)
      if (inherits(p, "error")) {
        plot_errors[[name]] <- conditionMessage(p)
      } else {
        p <- sb_clinician_clean_plot(sb_compact_plot(p), display_models)
        if (name == "performance") p <- p + ggplot2::labs(title = "Model comparison", x = "Model")
        overview[[name]] <- p
      }
    }
  }
  details <- list()
  if (isTRUE(plots) || (is.null(plots) && length(plan$pipelines) == 1L)) {
    sb_with_selector_reuse(for (name in names(plan$pipelines)) {
      fit <- sb_fit_final(task, plan$pipelines[[name]], capture_diagnostics = TRUE)
      figures <- NULL
      if (fit$status == "success") {
        figures <- tryCatch(sb_plots(fit), error = function(error) {
          list(plots = list(), unavailable = list(algorithm = conditionMessage(error)))
        })
      }
      details[[name]] <- list(fit = fit, figures = figures)
    })
  }
  output <- list(
    comparison = result, pipelines = plan$pipelines, unavailable = plan$unavailable,
    details = details, plots = overview, plot_errors = plot_errors, task_type = type,
    development_task = task
  )
  output$issues <- sb_run_issues(output)
  structure(output, class = c("stormbreaker_run", "list"))
}

#' @rdname sb_compare_cv
#' @export
bi_run <- function(data, outcome = NULL, predictors = NULL, id = NULL, models = NULL,
                   stages = if (is.null(models)) c(1L, 2L) else 1L,
                   selector = NULL, folds = 5L, seed = 1L, plots = NULL,
                   pipelines = NULL, decision_thresholds = NULL) {
  task <- sb_doctor_data_task(data, "binary", id, predictors, outcome,
    decision_thresholds = decision_thresholds)
  sb_doctor_run(task, "binary", models, stages, selector, folds, seed, plots, pipelines)
}

#' @rdname sb_compare_cv
#' @export
surv_run <- function(data, time = NULL, event = NULL, predictors = NULL, id = NULL, models = NULL,
                     stages = if (is.null(models)) c(1L, 2L, 3L) else 1L,
                     selector = NULL, folds = 5L, seed = 1L, plots = NULL,
                     time_unit = "month", pipelines = NULL,
                     prediction_horizon = NULL, decision_thresholds = NULL) {
  task <- sb_doctor_data_task(data, "survival", id, predictors, time = time,
    event = event, time_unit = time_unit, prediction_horizon = prediction_horizon,
    decision_thresholds = decision_thresholds)
  sb_doctor_run(task, "survival", models, stages, selector, folds, seed, plots, pipelines)
}

#' @rdname multi_run
#' @export
multi_run <- function(data, outcome = NULL, predictors = NULL, id = NULL,
                      models = NULL, folds = 5L, seed = 1L, plots = NULL) {
  task <- sb_doctor_data_task(data, "multiclass", id, predictors, outcome)
  menu <- sb_models("multiclass")
  menu <- menu[menu$supported & !menu$method %in% sb_multiclass_standalone_selectors(), ]
  if (is.null(models)) {
    methods <- menu$method
  } else {
    index <- match(tolower(models), tolower(menu$name))
    if (anyNA(index)) stop("Choose supported multiclass learner names from sb_models('multiclass').")
    methods <- menu$method[index]
  }
  list(
    comparison = sb_compare_multiclass(task, methods = methods, folds = folds, seed = seed),
    details = if (isTRUE(plots) || (is.null(plots) && length(methods) == 1L)) {
      sb_analyze_methods(task, methods)
    } else {
      NULL
    },
    methods = methods, task_type = "multiclass"
  )
}

# ---- doctor-issues ----

#' @rdname sb_complete_run
#' @export
sb_run_issues <- function(x) {
  if (!is.list(x) || !inherits(x$comparison, "stormbreaker_cv_comparison")) stop("Supply a bi_run or surv_run result.")
  rules <- list(
    dependency = c("Package .*required|Install |not installed|there is no package", "\u7f3a\u5c11\u7b97\u6cd5\u4f9d\u8d56\uff0c\u8bf7\u5b89\u88c5\u63d0\u793a\u7684\u8f6f\u4ef6\u5305\u540e\u518d\u8fd0\u884c\u3002"),
    features = c("no features|no predictors|no genes|empty.*select|retain 1 or less|No features met|nonempty set", "\u7b5b\u9009\u540e\u53d8\u91cf\u4e0d\u8db3\uff1b\u53ef\u51cf\u5c11\u7b5b\u9009\u6b65\u9aa4\uff0c\u4e0d\u80fd\u968f\u610f\u8865\u57fa\u56e0\u51d1\u6570\u3002"),
    evidence = c("not retained|not captured|not saved|capture_diagnostics|disabled importance|OOB.*disabled", "\u6a21\u578b\u6ca1\u6709\u4fdd\u5b58\u8fd9\u9879\u7ed8\u56fe\u6570\u636e\uff1b\u4e0d\u5f71\u54cd\u5df2\u6210\u529f\u751f\u6210\u7684\u5176\u4ed6\u7ed3\u679c\u3002"),
    data = c("single.class|both.*classes|0/1|non-finite|finite|missing values|schema|columns.*match|same.*columns", "\u68c0\u67e5\u5206\u7ec4\u3001\u7f3a\u5931\u503c\u53ca\u8bad\u7ec3\u4e0e\u9a8c\u8bc1\u6570\u636e\u7684\u5217\u540d\uff1b\u4e5f\u53ef\u80fd\u662f\u8be5\u7b97\u6cd5\u5728\u8fd9\u4efd\u6570\u636e\u4e0a\u65e0\u6cd5\u4f30\u8ba1\u3002"),
    unsupported = c("not admitted|not supported|not implemented|No native|nonlinear kernel|cannot.*group|No specialized|no single exact|not a variable.selection", "\u8fd9\u4e2a\u7ec4\u5408\u6216\u56fe\u5f62\u4e0d\u9002\u7528\uff1b\u4e0d\u4f1a\u7528\u53e6\u4e00\u79cd\u7b97\u6cd5\u6216\u865a\u6784\u56fe\u5f62\u4ee3\u66ff\u3002")
  )
  rows <- list()
  add <- function(model, stage, reason) {
    if (is.null(reason) || is.na(reason) || !nzchar(reason)) reason <- "No detailed reason was saved."
    category <- "other"
    explanation <- "\u9700\u8981\u67e5\u770b\u539f\u59cb\u539f\u56e0\uff0c\u6682\u4e0d\u80fd\u81ea\u52a8\u5224\u65ad\u3002"
    for (k in names(rules)) {
      if (grepl(rules[[k]][1], reason, ignore.case = TRUE)) {
        category <- k
        explanation <- rules[[k]][2]
        break
      }
    }
    rows[[length(rows) + 1L]] <<- data.frame(
      model = model, stage = stage,
      category = category, explanation = explanation, reason = reason
    )
  }
  failures <- sb_failure_table(x$comparison)
  labels <- stats::setNames(as.character(names(x$pipelines)), vapply(x$pipelines, function(p) p$id, character(1)))
  for (i in seq_len(nrow(failures))) {
    f <- failures[i, ]
    if (f$issue_type == "pipeline" && any(failures$pipeline_id == f$pipeline_id & failures$issue_type == "run")) next
    model <- unname(labels[f$pipeline_id])
    if (is.na(model)) model <- f$pipeline_id
    reason <- f$error_message
    if (is.na(reason) || !nzchar(reason)) reason <- paste(f$metric, f$status)
    add(model, if (is.na(f$fold)) f$stage else paste0(f$stage, " ", f$fold), reason)
  }
  for (n in names(x$unavailable)) add(n, "combination", x$unavailable[[n]])
  for (n in names(x$plot_errors)) add("", "overview_plot", x$plot_errors[[n]])
  for (n in names(x$details)) {
    d <- x$details[[n]]
    if (d$fit$status != "success") add(n, "final_fit", d$fit$error_message)
    for (v in names(d$figures$unavailable)) {
      issue <- d$figures$unavailable[[v]]
      if (is.data.frame(issue)) {
        for (i in which(issue$status != "available")) {
          add(n, paste("plot", issue$view[i]), issue$reason[i])
        }
      } else {
        add(n, paste("plot", v), issue)
      }
    }
    native <- d$figures$native$status
    if (is.null(d$figures$unavailable$native) && is.data.frame(native)) for (i in which(native$status != "available")) add(n, paste("plot", native$view[i]), native$reason[i])
  }
  if (!length(rows)) {
    return(data.frame(
      model = character(), stage = character(),
      category = character(), explanation = character(), reason = character()
    ))
  }
  do.call(rbind, rows)
}

# ---- doctor-summary ----

#' @export
print.stormbreaker_run <- function(x, ...) {
  comparison <- x$comparison$comparison
  successful <- sum(comparison$status == "success")
  cat("<stormbreaker:", x$task_type, "comparison>\n")
  cat(
    "Models:", nrow(comparison), "| CV successful:", successful,
    "| Other:", nrow(comparison) - successful, "\n"
  )
  cat("Final fits retained:", length(x$details), "| Issues:", nrow(x$issues), "\n")
  cat("Model names: names(run$pipelines)\n")
  cat("Full comparison: run$comparison$comparison | Issues: run$issues\n")
  cat("Next: choose a pipeline with sb_complete_run().\n")
  invisible(x)
}

#' @export
print.stormbreaker_completed_run <- function(x, ...) {
  successful <- sum(vapply(x$evaluations, \(evaluation) evaluation$status == "success", logical(1)))
  cat("<stormbreaker:", x$task_type, "completed run>\n")
  cat("Pipeline:", x$pipeline, "| Fit:", x$fit$status, "\n")
  cat(
    "External cohorts evaluated:", length(x$evaluations),
    "| Successful:", successful, "\n"
  )
  cat("Results: $genes, $predictions, $metrics, $plots\n")
  cat("Unavailable results and failures: $issues\n")
  cat("Export: sb_export_run(finished, directory = \"results-new\")\n")
  invisible(x)
}

# ---- start ----

#' @rdname stormbreaker-package
#' @export
sb_start <- function(view = interactive()) {
  if (!is.logical(view) || length(view) != 1L || is.na(view)) sb_abort(
    "view must be TRUE or FALSE.", "stormbreaker_bad_report_input")
  html <- paste(readLines(system.file("start.html", package = "stormbreaker",
    mustWork = TRUE), encoding = "UTF-8", warn = FALSE), collapse = "\n")
  examples <- new.env(parent = baseenv())
  sys.source(system.file("examples", "advanced-examples.R",
    package = "stormbreaker", mustWork = TRUE), envir = examples)
  for (action in c("compare", "validate")) {
    recipe <- examples$.stormbreaker_advanced_examples[[paste0("doctor-", action)]]
    expressions <- utils::head(as.list(body(recipe))[-1L], -1L)
    code <- paste(c("library(stormbreaker)", vapply(expressions,
      function(expr) paste(deparse(expr, width.cutoff = 80L), collapse = "\n"),
      character(1))), collapse = "\n")
    html <- gsub(paste0("@@", toupper(action), "@@"), sb_report_escape(code), html, fixed = TRUE)
  }
  file <- tempfile("stormbreaker-start-", fileext = ".html")
  writeLines(enc2utf8(html), file, useBytes = TRUE)
  if (view) utils::browseURL(file)
  invisible(normalizePath(file, winslash = "/"))
}
