#' Pipeline specifications
#' @name sb_pipeline
#' @description Specify feature-selection stages, learners and ordered modeling routes.
#' @rawRd \arguments{
#' \item{selector}{
#' \strong{For \code{sb_pipeline()}:} Feature selector identifier, including \code{"none"},
#' \code{"lasso_binary_nonzero"}, \code{"stepglm_binary_aic"},
#' \code{"glmboost_binary_nonzero"}, \code{"rf_binary_importance"},
#' \code{"gamboost_binary_nonzero"} (active native binomial bbs spline terms),
#' \code{"xgboost_binary_importance"}, \code{"xgboost_cox_importance"} (positive
#' native tree gain; 100 rounds, eta 0.1, depth 3, single thread, seed 1),
#' \code{"enet_binary_nonzero"}, \code{"gbm_binary_importance"}, \code{"enet_cox_nonzero"},
#' \code{"gbm_cox_importance"}, \code{"glmboost_cox_nonzero"}, \code{"superpc_survival_scores"},
#' \code{"gamboost_cox_nonzero"} (active native Cox bbs spline terms),
#' \code{"ridge_cox_weights"}, \code{"plsrcox_weights"}, \code{"svm_survival_weights"},
#' \code{"fixed_features"}, \code{"univ_cox_top_k"}, \code{"stepcox"},
#' \code{"univ_binary_auc_top_k"}, \code{"abess_binary_gic_nonzero"},
#' \code{"abess_cox_gic_nonzero"}, \code{"boruta_binary_roughfix"},
#' \code{"boruta_survival_roughfix"}, \code{"varselrf_binary_oob_1se"},
#' \code{"sis_binary_nlogn"}, \code{"vsurf_binary_prediction"},
#' \code{"mrmr_binary_nlogn"}, \code{"jmi_binary_nlogn"},
#' \code{"cmim_binary_nlogn"}, \code{"disr_binary_nlogn"},
#' \code{"jmim_binary_nlogn"}, \code{"lpc_binary_nlogn"},
#' \code{"lpc_cox_nlogn"}, \code{"stabsel_binary_lasso"},
#' \code{"stabsel_cox_lasso"}, or an admitted MIME survival selector.
#' 
#' \strong{For \code{bi_pipeline()}:} Legacy selector specification, retained for compatibility.
#' New binary code supplies an ordered model vector instead. Survival may
#' use \code{\link[=surv_selector]{surv_selector()}}.
#' }
#' \item{learner}{
#' Learner identifier: \code{"coxph"}, \code{"glmnet_cox"},
#' \code{"logistic"}, \code{"glmnet_logistic"}, \code{"gaussian_naive_bayes"}, or
#' \code{"lightgbm_binary"}, \code{"lightgbm_binary_dart"},
#' \code{"lightgbm_binary_goss"}, \code{"lightgbm_binary_rf"},
#' \code{"catboost_binary"}, \code{"cancerclass_binary"},
#' \code{"ranger_binary"}, \code{"ranger_extratrees_binary"},
#' \code{"ranger_hellinger_binary"}, \code{"adaboost_binary"},
#' \code{"logitboost_binary"}, \code{"lda_binary"}, \code{"qda_binary"},
#' \code{"kernel_naive_bayes_binary"}, \code{"rda_binary"},
#' \code{"knn_binary"}, \code{"lvq1_binary"},
#' \code{"nearest_centroid_binary"}, \code{"bernoulli_naive_bayes_binary"},
#' \code{"olvq1_binary"}, \code{"lvq21_binary"},
#' \code{"lvq3_binary"}, \code{"nnet_binary"},
#' \code{"decision_tree_binary"},
#' \code{"random_forest_binary"}, \code{"svm_linear_binary"},
#' \code{"svm_radial_binary"},
#' \code{"svm_radial_weighted_binary"},
#' \code{"svm_polynomial_binary"}, \code{"svm_sigmoid_binary"},
#' \code{"gbm_binary"}, \code{"bart_binary"}, \code{"c50_binary"},
#' \code{"pre_binary_rule_ensemble"},
#' \code{"kernel_weighted_knn_binary"},
#' \code{"gaussian_process_binary"}, \code{"pamr_binary"},
#' \code{"neuralnet_multilayer_binary"},
#' \code{"ctree_binary"}, \code{"cforest_binary"},
#' \code{"bagging_binary"}, \code{"gamboost_binary"},
#' \code{"blackboost_binary"}, \code{"loclda_binary"},
#' \code{"aorsf_binary"},
#' \code{"earth_binary"}, \code{"fda_binary"}, \code{"mda_binary"},
#' \code{"liblinear_l2_logistic_binary"},
#' \code{"liblinear_l1_logistic_binary"},
#' \code{"liblinear_l2_l2loss_svc_binary"},
#' \code{"liblinear_l2_l1loss_svc_binary"},
#' \code{"liblinear_l1_l2loss_svc_binary"}, \code{"evtree_binary"},
#' \code{"sda_binary"}, \code{"rotation_forest_binary"},
#' \code{"rrf_binary"}, \code{"pda_binary"}, \code{"slda_binary"},
#' \code{"tree_binary"}, \code{"c50_rules_binary"},
#' \code{"random_ferns_binary"}, \code{"robust_logistic_mqle"},
#' \code{"bayesglm_logistic_cauchy"}, \code{"bst_hinge_binary"},
#' \code{"hda_binary"}, \code{"lssvm_radial_binary"},
#' \code{"stepglm_binary"}, \code{"glmboost_binary"},
#' \code{"plsrglm_binary"}, \code{"plsda_binary"},
#' \code{"splsda_binary"}, \code{"ncvreg_mcp_binary_bic"},
#' \code{"ncvreg_scad_binary_bic"}, \code{"rlt_reinforcement_binary"},
#' \code{"wsrf_binary"}, \code{"random_glm_binary"},
#' \code{"hdda_binary_cattell"},
#' \code{"abess_binary_gic"},
#' \code{"edda_binary_bic"}, \code{"mclustda_binary_bic"},
#' \code{"xgboost_binary"}, \code{"xgboost_binary_dart"},
#' \code{"xgboost_binary_linear"}, \code{"catboost_cox"},
#' \code{"catboost_aft_survival"},
#' \code{"gbm_cox"},
#' \code{"gbm_cox_cv"},
#' \code{"xgboost_cox"}, \code{"xgboost_cox_dart"},
#' \code{"xgboost_cox_linear"}, \code{"xgboost_aft"},
#' \code{"xgboost_aft_dart"}, \code{"xgboost_aft_linear"},
#' \code{"survreg_weibull"},
#' \code{"survreg_exponential"}, \code{"survreg_lognormal"},
#' \code{"survreg_loglogistic"},
#' \code{"flexsurv_gamma"}, \code{"flexsurv_gompertz"},
#' \code{"flexsurv_gengamma"},
#' \code{"flexsurvspline_hazard_k1"},
#' \code{"flexsurvspline_odds_k1"},
#' \code{"flexsurvspline_normal_k1"},
#' \code{"rpart_survival"}, \code{"rsf_survival"},
#' \code{"ctree_survival"}, \code{"cforest_survival"},
#' \code{"bart_survival"}, \code{"aorsf_survival"},
#' \code{"bart_aft_survival"}, \code{"pre_survival_rule_ensemble"},
#' \code{"penaft_gehan_elastic_net"},
#' \code{"aalen_additive_survival"}, \code{"prop_odds_survival"},
#' \code{"pamr_survival"},
#' \code{"gamboost_cox"}, \code{"glmboost_cox"},
#' \code{"blackboost_cox"}, \code{"gamboost_weibull_aft"},
#' \code{"gamboost_loglogistic_aft"},
#' \code{"gamboost_lognormal_aft"}, \code{"glmboost_weibull_aft"},
#' \code{"glmboost_loglogistic_aft"},
#' \code{"glmboost_lognormal_aft"}, \code{"blackboost_weibull_aft"},
#' \code{"blackboost_loglogistic_aft"},
#' \code{"blackboost_lognormal_aft"},
#' \code{"gamboost_gehan_aft"}, \code{"gamboost_unoc"},
#' \code{"glmboost_gehan_aft"}, \code{"glmboost_unoc"},
#' \code{"blackboost_gehan_aft"}, \code{"blackboost_unoc"},
#' \code{"ranger_survival"},
#' \code{"ranger_extratrees_survival"},
#' \code{"ranger_cindex_survival"}, \code{"ranger_maxstat_survival"},
#' \code{"coxboost"}, \code{"coxboost_cv"}, \code{"plsrcox"}, \code{"plsrcox_cv"},
#' \code{"superpc"}, \code{"superpc_cv"}, or
#' \code{"survivalsvm"}, \code{"survivalsvm_vanbelle1"},
#' \code{"survivalsvm_vanbelle2"}, \code{"survivalsvm_hybrid"},
#' \code{"survivalsvm_additive"},
#' \code{"survivalsvm_vanbelle1_additive"},
#' \code{"survivalsvm_vanbelle2_additive"},
#' \code{"survivalsvm_hybrid_additive"},
#' \code{"survivalsvm_rbf"},
#' \code{"survivalsvm_vanbelle1_rbf"},
#' \code{"survivalsvm_vanbelle2_rbf"},
#' \code{"survivalsvm_hybrid_rbf"},
#' \code{"ncvreg_mcp_cox_bic"}, \code{"ncvreg_scad_cox_bic"},
#' \code{"rlt_reinforcement_survival"}, \code{"coxphf_firth"},
#' \code{"bagged_knn_survival"}, \code{"akritas_survival"}, \code{"ahaz_lin_ying"},
#' \code{"ahaz_lasso_bic"}, \code{"ahaz_sscad_bic"},
#' \code{"bujar_linear_boost"}, \code{"abess_cox_gic"}.
#' }
#' \item{selector_n_features}{
#' Fixed number of features selected inside each
#' analysis fold for a Top-K selector. It must be declared before resampling,
#' must be at least 2 when paired with glmnet, and must be \code{NULL} for
#' \code{"none"} or \code{"fixed_features"}.
#' }
#' \item{inner_folds}{
#' Number of inner folds used to select \code{lambda} for a
#' glmnet learner. It must be at least 3 and must be \code{NULL} for every
#' non-glmnet learner.
#' }
#' \item{inner_seed}{
#' Positive integer seed for inner-fold assignment for a
#' glmnet learner. It must be \code{NULL} for every non-glmnet learner.
#' }
#' \item{glmnet_alpha}{
#' Fixed mixing parameter between 0 and 1 for a glmnet
#' learner: 0 is Ridge, 1 is Lasso, and intermediate values are Elastic Net.
#' \code{NULL} uses the backward-compatible Lasso default of 1. It must be \code{NULL}
#' for every non-glmnet learner.
#' }
#' \item{naive_bayes_var_smoothing}{
#' Non-negative variance-smoothing factor for
#' Gaussian Naive Bayes. \code{NULL} uses \code{1e-9}. It must be \code{NULL} for every
#' other learner.
#' }
#' \item{lightgbm_nrounds}{
#' Fixed number of training rounds for a fixed
#' LightGBM learner. \code{NULL} uses 100.
#' }
#' \item{lightgbm_learning_rate}{
#' Fixed positive learning rate for a fixed
#' LightGBM learner. \code{NULL} uses 0.1.
#' }
#' \item{lightgbm_num_leaves}{
#' Fixed maximum number of leaves per tree for a
#' fixed LightGBM learner. \code{NULL} uses 31.
#' }
#' \item{lightgbm_min_data_in_leaf}{
#' Fixed non-negative minimum leaf size for a
#' fixed LightGBM learner. \code{NULL} uses 20.
#' }
#' \item{lightgbm_seed}{
#' Fixed positive LightGBM seed. \code{NULL} uses 1.
#' }
#' \item{catboost_iterations}{
#' Fixed number of boosting iterations for
#' \code{learner = "catboost_binary"}, \code{"catboost_cox"}, or
#' \code{"catboost_aft_survival"}. \code{NULL} uses 100.
#' }
#' \item{catboost_learning_rate}{
#' Fixed positive learning rate for
#' \code{learner = "catboost_binary"}, \code{"catboost_cox"}, or
#' \code{"catboost_aft_survival"}. \code{NULL} uses 0.1.
#' }
#' \item{catboost_depth}{
#' Fixed symmetric-tree depth from 1 through 16 for
#' \code{learner = "catboost_binary"}, \code{"catboost_cox"}, or
#' \code{"catboost_aft_survival"}. \code{NULL} uses 6.
#' }
#' \item{catboost_l2_leaf_reg}{
#' Fixed positive L2 leaf regularization for
#' \code{learner = "catboost_binary"}, \code{"catboost_cox"}, or
#' \code{"catboost_aft_survival"}. \code{NULL} uses 3.
#' }
#' \item{catboost_seed}{
#' Fixed positive CatBoost seed. \code{NULL} uses 1.
#' }
#' \item{ranger_num_trees}{
#' Fixed number of trees for a fixed ranger learner.
#' \code{NULL} uses 500.
#' }
#' \item{ranger_mtry}{
#' Optional fixed number of candidate split variables for
#' a fixed ranger learner. \code{NULL} uses the deterministic rule \code{floor(sqrt(p))}
#' within each fitted feature set.
#' }
#' \item{ranger_min_node_size}{
#' Fixed minimum split-node size for a fixed ranger
#' learner. \code{NULL} uses 10 for binary probability forests and 3 for survival
#' forests.
#' }
#' \item{ranger_seed}{
#' Fixed positive ranger seed. \code{NULL} uses 1.
#' }
#' \item{selector_features}{
#' Ordered, unique numeric predictor names for
#' \code{selector = "fixed_features"}. They are declared before resampling and are
#' never learned from outcomes. It must be \code{NULL} for other selectors.
#' }
#' \item{feature_set_id}{
#' One non-empty caller-supplied label that distinguishes
#' a fixed feature block in pipeline IDs and audit output. It is required for
#' \code{selector = "fixed_features"} and must be \code{NULL} otherwise.
#' }
#' \item{glmnet_standardize}{
#' Whether glmnet standardizes predictors inside each
#' fitted model. \code{NULL} uses the backward-compatible default \code{TRUE}. Set
#' \code{FALSE} only when the supplied predictor scale is already appropriate for
#' penalization. It must be \code{NULL} for every non-glmnet learner. This does not
#' add imputation, normalization, or any other preprocessing step.
#' }
#' \item{glmnet_lambda_rule}{
#' Rule used inside each glmnet inner cross-validation
#' to select \code{lambda}: \code{"lambda.min"} chooses the optimum mean inner-CV
#' target, while \code{"lambda.1se"} chooses the largest \code{lambda} within one
#' standard error of that optimum. \code{NULL} uses the backward-compatible
#' \code{"lambda.min"} default. The rule is declared before resampling and must be
#' \code{NULL} for every non-glmnet learner.
#' }
#' \item{glmnet_tuning_measure}{
#' Inner cross-validation target for a glmnet
#' learner. \code{NULL} or \code{"deviance"} preserves the existing likelihood-loss
#' target. Binary logistic pipelines may instead predeclare \code{"auc"} for ROC
#' AUC; Cox pipelines may instead predeclare \code{"C"} for Harrell concordance.
#' The target is fixed before resampling and is never chosen from outer-fold
#' or external performance. It must be \code{NULL} for every other learner.
#' }
#' \item{stepcox_direction}{
#' Stepwise Cox direction: \code{"forward"}, \code{"backward"},
#' or \code{"both"}. It must be declared for \code{selector = "stepcox"} and for a
#' MIME sequential selector ending in StepCox, and must be \code{NULL} otherwise.
#' Forward selection starts from a null Cox model; backward and both-direction
#' selection start from the full model.
#' }
#' \item{gbm_n_trees}{
#' Fixed number of boosting trees for
#' \code{learner = "gbm_cox"}. \code{NULL} uses 100.
#' }
#' \item{gbm_interaction_depth}{
#' Fixed maximum interaction depth for
#' \code{learner = "gbm_cox"}. \code{NULL} uses 1.
#' }
#' \item{gbm_min_obs_node}{
#' Fixed minimum number of observations in a terminal
#' node for \code{learner = "gbm_cox"}. \code{NULL} uses 10.
#' }
#' \item{gbm_shrinkage}{
#' Fixed positive GBM learning rate for
#' \code{learner = "gbm_cox"}. \code{NULL} uses 0.1.
#' }
#' \item{gbm_bag_fraction}{
#' Fixed GBM training-row fraction in \verb{(0, 1]} for each
#' tree. \code{NULL} uses 1, so the default does not subsample rows.
#' }
#' \item{gbm_seed}{
#' Fixed positive GBM seed. \code{NULL} uses 1.
#' }
#' \item{coxboost_stepno}{
#' Fixed CoxBoost boosting-step count. \code{NULL} uses 100.
#' }
#' \item{coxboost_penalty}{
#' Optional fixed positive CoxBoost penalty. \code{NULL}
#' uses CoxBoost's documented fold-local rule of nine times the analysis-fold
#' event count.
#' }
#' \item{plsrcox_components}{
#' Fixed number of plsRcox components. \code{NULL} uses 2.
#' }
#' \item{superpc_threshold}{
#' Fixed non-negative SuperPC feature-score threshold.
#' \code{NULL} uses 1; STORMBREAKER does not cross-validate this value.
#' }
#' \item{superpc_components}{
#' Fixed SuperPC component count from 1 through 3.
#' \code{NULL} uses 1.
#' }
#' \item{survivalsvm_gamma_mu}{
#' Fixed positive survival-SVM regularization value.
#' \code{NULL} uses 1. Each learner uses its registered linear, additive, or RBF
#' kernel and a fixed quadratic-programming solver. The hybrid formulation
#' applies this scalar equally to both native regularizers. RBF learners use
#' median squared distance of standardized analysis predictors as the native
#' kernel parameter, without searching external or assessment outcomes.
#' }
#' \item{adaboost_iterations}{
#' Fixed maximum AdaBoost.M1 iteration count for
#' \code{learner = "adaboost_binary"}. \code{NULL} uses 50. The base learner is
#' the classification tree used by \code{fastAdaboost}, with \code{cp = 0}; the native
#' implementation can finish earlier when boosting converges.
#' }
#' \item{logitboost_iterations}{
#' Fixed decision-stump iteration count for
#' \code{learner = "logitboost_binary"}. \code{NULL} uses 51.
#' }
#' \item{stepglm_direction}{
#' Stepwise logistic direction: \code{"forward"},
#' \code{"backward"}, or \code{"both"} for \code{learner = "stepglm_binary"}.
#' \code{NULL} uses \code{"forward"}. Forward selection starts from the null logistic
#' model; backward and both-direction selection start from the full model. It
#' must be \code{NULL} for every other learner.
#' }
#' \item{selector_inner_folds}{
#' Number of analysis-local folds for
#' \code{"coxboost_nonzero"}, \code{"lasso_nonzero"}, \code{"enet_cox_nonzero"}, or sequential StepCox
#' selectors. It must be at least 3 and must be \code{NULL} for other selectors.
#' }
#' \item{selector_seed}{
#' Positive integer seed for an RSF, CoxBoost, or Lasso
#' outcome-learned MIME selector. It must be \code{NULL} for other selectors.
#' }
#' \item{selector_rsf_nodesize}{
#' Positive terminal-node size for
#' \code{"rsf_minimal_depth"} and \code{"rsf_then_stepcox"}. Those selectors fit 1,000
#' logrank trees before native conservative minimal-depth selection. It must
#' be \code{NULL} for other selectors.
#' }
#' \item{learner_tuning_folds}{
#' Number of analysis-fold internal folds for
#' \code{coxboost_cv}, \code{plsrcox_cv}, \code{superpc_cv}, or \code{gbm_cox_cv}. \code{NULL} uses
#' the reference default: 5 for plsRcox and 10 for the other learners.
#' }
#' \item{learner_seed}{
#' Positive integer seed for an analysis-fold tuned
#' survival learner. \code{NULL} uses 1. It must be \code{NULL} for other learners.
#' }
#' \item{gbm_cv_max_trees}{
#' Maximum tree count searched by \code{gbm_cox_cv}.
#' With one allowed tree, fit that sole candidate without an inner tree-count
#' search; tuning records \code{single_candidate_no_search} and no CV curve.
#' \code{NULL} uses the reference value 10,000 and must be \code{NULL} otherwise.
#' }
#' \item{selector_options}{
#' Named list for supported selectors. Boruta accepts
#' \code{num_trees = 500}, \code{maximum_runs = 100} (at least 11), and \code{seed = 1}.
#' Smaller budgets change the selection specification, not merely
#' runtime. Declare them before resampling; nondefault budgets appear in IDs.
#' Binary model selectors accept the following options.
#' LASSO: \code{inner_folds = 5}, \code{seed = 1}, \code{lambda_rule = "lambda.1se"}
#' (also accepts \code{"lambda.min"}); selects nonzero coefficients.
#' Enet uses the same options plus \code{alpha = 0.5} (0 < alpha <= 1).
#' GBM: \code{seed = 1}, fixed 100 depth-one trees, shrinkage 0.1, minimum
#' node size 10 and bag fraction 1; retains positive relative influence.
#' This can retain all predictors and is not a significance test.
#' Stepglm: \code{direction = "both"} (also \code{"forward"}, \code{"backward"}); AIC.
#' glmBoost: \code{iterations = 100}, \code{shrinkage = 0.1}; nonzero coefficients.
#' RF: \code{trees = 500}, \code{seed = 1}, \code{importance = "permutation"}
#' (also \code{"gini"}), \code{threshold = 0}; unscaled importance strictly above
#' threshold by default. With \code{top_n}, rank by unscaled importance (ties
#' follow input feature order), take up to that many and drop constant features.
#' This replaces the threshold rule. With \code{spearman_cutoff} in (0, 1],
#' keep a candidate only if its absolute training Spearman correlation with
#' every higher-ranked retained feature does not exceed the cutoff. Do not
#' refill discarded slots. Selection is repeated inside training folds. Empty selections
#' fail without adding substitute features; glmnet needs at least two retained
#' features. Must be \code{NULL} for other selectors.
#' }
#' \item{model}{
#' For binary tasks, one name or an ordered character vector from
#' \code{\link[=sb_models]{sb_models()}}, for example c("Lasso", "RF", "GBM"). Intermediate
#' methods extract features; the last method predicts. Survival takes one
#' final model name. No model is fitted by this constructor.
#' }
#' \item{methods}{
#' \strong{For \code{bi_selector()}:} Ordered binary selector names from historical code.
#' 
#' \strong{For \code{surv_selector()}:} Ordered survival names such as c("Lasso", "StepCox").
#' Existing native defaults, selection rules and fold-local fitting are kept.
#' Binary tasks use ordered names directly in \code{\link[=bi_pipeline]{bi_pipeline()}}, not a separate
#' binary selector entrypoint.
#' }
#' \item{spec}{
#' One \code{stormbreaker_pipeline} object or one named list of
#' arguments accepted by \code{\link[=sb_pipeline]{sb_pipeline()}}.
#' }
#' \item{task_type}{
#' Optional target task, either \code{"survival"} or \code{"binary"}.
#' When supplied, the pipeline must match it.
#' }
#' }
#' @rawRd \value{See the individual entries in Function reference.}
#' @rawRd \section{Function reference}{
#' \subsection{\code{sb_pipeline()}}{
#' The minimum core supports Cox proportional hazards, Firth-corrected Cox,
#' glmnet Cox, CoxBoost,
#' plsRcox, SuperPC, twelve survival-SVM configurations across four native
#' formulations and linear, additive, and radial-basis-function kernels,
#' unpenalized, mean bias-reduced, Cauchy-prior Bayesian, and glmnet logistic
#' regression, componentwise linear hinge-loss boosting, Gaussian Naive Bayes,
#' a one-dimensional heteroscedastic discriminant transformation with a native
#' Gaussian naive-Bayes rule, LDA, QDA, deterministic KNN, four native
#' learning-vector-quantization update algorithms, classification trees, rule
#' sets, Gini, ExtraTrees, Hellinger, and weighted-subspace probability forests,
#' linear, radial,
#' and class-weighted radial SVMs,
#' stepwise logistic regression, component-wise gradient boosting,
#' partial least squares logistic regression, fixed AdaBoost.M1 and LogitBoost,
#' fixed-configuration LightGBM GBDT, DART, GOSS, and random-forest training,
#' CatBoost binary, Cox, and Normal SurvivalAft models, GBM, and XGBoost
#' gbtree, DART, and linear boosters, GBM-Cox,
#' GamBoost-Cox, GlmBoost-Cox, BlackBoost-Cox, nine mboost AFT likelihood-by-
#' base-learner specifications, six mboost Gehan/Uno-C ranking specifications,
#' XGBoost-Cox, XGBoost normal AFT,
#' ranger logrank, ExtraTrees, C-index, and maxstat survival forests, RLT
#' reinforcement classification and survival forests, bagged k-nearest-neighbor
#' survival prediction, Lin-Ying, Lasso, and SCAD additive-hazards regression,
#' Buckley-James componentwise linear boosting, Weibull,
#' exponential, log-normal, and log-logistic SurvReg, Rpart
#' survival trees, randomForestSRC survival forests, partykit conditional-inference survival
#' trees and forests, discrete-time survival BART, oblique survival forests,
#' Aalen additive-hazards regression, and semiparametric proportional-odds
#' survival regression. Fixed HDA leaves its class-specific covariance estimates
#' unregularized, so each must be nonsingular; ordinarily, selected predictors
#' must be fewer than the samples in every analysis-fold class and not perfectly
#' collinear. Pipeline task type is
#' inferred from the learner and must match the task supplied to \code{\link[=sb_run_cv]{sb_run_cv()}}.
#' The glmnet mixing parameter \code{alpha} is fixed before resampling; \code{alpha = 0}
#' gives Ridge, \code{alpha = 1} gives Lasso, and values strictly between 0 and 1
#' give Elastic Net. Glmnet's predictor standardization is declared in the
#' pipeline and is performed separately by each inner analysis fit and final
#' outer-analysis refit; it is not a package-wide normalization or imputation
#' step. \code{lambda} is selected separately inside every outer analysis fold by a
#' predeclared \code{"lambda.min"} or \code{"lambda.1se"} rule; the package does not
#' compare those rules or choose one from outer-fold performance. Outer
#' assessment patients never enter inner tuning. When a
#' matching supervised Top-K selector is requested with glmnet, that selector
#' is refitted separately in every inner analysis fold through glmnet's
#' function-valued \code{exclude} filter. The outer-analysis selector defines only
#' the master path and final outer fit. Glmnet Cox explicitly uses Breslow
#' tied-event handling across dependency versions. Ordinary Cox, glmnet Cox,
#' fixed CoxBoost, plsRcox, SuperPC, and the calibrated linear survival-SVM
#' report a portable zero-reference linear predictor \code{x_beta}. Additive and RBF
#' survival-SVM report their negative native kernel scores multiplied by their
#' analysis-fold Cox calibration coefficients. CatBoost Cox reports its native
#' raw margin. CatBoost SurvivalAft reports the negative native log-time margin
#' and reads horizon survival directly from its fixed Normal distribution.
#' GBM-Cox and XGBoost-Cox report native link-scale log hazards. XGBoost normal
#' AFT reports the negative native log-time margin;
#' mboost AFT reports the negative native distributional log-time link. Mboost
#' Gehan and Uno-C ranking learners report the negative native ranking link.
#' Ahaz learners report their native additive-hazard linear predictor; a
#' declared-horizon probability uses analysis-fold Cox/Breslow calibration
#' because the native additive-hazard formula need not be a proper probability.
#' For tied deaths only, ahaz fitting uses seed 10101 to permute small negative
#' offsets within each tied group. The offset window is bounded by 1e-8 times
#' the recorded time and one quarter of the preceding distinct-time gap;
#' one death retains its recorded time. Censored observations and all originally
#' distinct time orderings are preserved. This is a numerical approximation,
#' not additional clinical timing information. Calibration and evaluation use
#' the original endpoints. The fitted audit records the seed, affected counts
#' and maximum shift. Seed sensitivity should be assessed for reported models.
#' Bujar reports negative predicted log time and uses an analysis-fold residual
#' Kaplan-Meier curve for a declared horizon.
#' AFT-BART reports negative posterior mean log time and evaluates its native
#' normal-error posterior at a declared horizon. PRE supplies binary probability
#' or a Cox rule-ensemble link with analysis-fold Breslow survival. Its native
#' procedure first learns the rule library from the complete current analysis
#' fold, then uses five-fold coefficient-only CV on that fixed library to choose
#' \code{lambda.1se}. Outer assessment patients remain untouched, but rule induction
#' is not regenerated separately inside those native coefficient-CV folds.
#' Penalized
#' semiparametric AFT reports negative centered log time and uses an
#' analysis-fold residual Kaplan-Meier curve at a declared horizon.
#' MTLR fits coupled logistic tasks over a ceiling-square-root analysis-fold
#' time grid. It reports negative native mean predicted survival time as risk
#' and reads a declared-horizon probability directly from its fitted curve.
#' PAMR survival forms two censoring-aware groups inside the analysis fold and
#' reports the native early-survival-group posterior as risk. Declared-horizon
#' survival uses a Cox/Breslow calibration fitted on that same analysis fold.
#' Ranger and randomForestSRC survival forests report native mortality. The
#' RLT survival forest reconstructs survival from native discrete hazard
#' increments and reports negative step-function restricted mean survival over
#' the native event-time grid.
#' Each SurvReg AFT distribution reports negative fitted log time, Rpart
#' reports log relative event rate, and partykit learners report the terminal
#' or ensemble distribution's final cumulative hazard. Full Aalen regression
#' reports cumulative hazard at the last analysis-fold event time, while the
#' proportional-odds model reports its zero-reference log failure-odds ratio
#' \code{x_beta}. These risk scores are not survival probabilities.
#' StepCox is a supervised survival selector based on Cox partial-likelihood
#' AIC. Forward selection starts from the null Cox model; backward and both-
#' direction selection start from the full candidate model. Its direction is
#' declared before resampling. Selection is repeated in every outer analysis
#' fold and, when followed by glmnet Cox, in every inner analysis fold.
#' The binary AUC
#' Top-K selector ranks either expression direction inside each outer analysis
#' fold and is admitted with ordinary logistic regression, Gaussian Naive
#' Bayes, LDA, QDA, deterministic KNN, LVQ1, optimized LVQ1, LVQ2.1, LVQ3,
#' classification trees, C5.0 rule sets, random forests,
#' linear, radial, and class-weighted radial SVMs, stepwise logistic regression,
#' glmBoost, plsRglm, PLS-DA, sparse PLS-DA, fixed GBM, fixed LightGBM, fixed
#' CatBoost, fixed XGBoost,
#' fixed ranger split-rule forests, fixed AdaBoost.M1, fixed LogitBoost, Earth,
#' FDA, MDA,
#' LIBLINEAR, evolutionary trees, shrinkage discriminant analysis, rotation
#' forests, regularized random forests, or glmnet logistic regression. Every
#' registered survival learner is also admitted
#' with the fold-local univariable Cox Top-K selector. No
#' fixed tree learner uses
#' outer assessment patients for early stopping or tuning. Binary learners
#' report the probability of class \code{1}. Unsupported components raise an error
#' and are never replaced silently.
#' 
#' \strong{Returns:} A \code{stormbreaker_pipeline} specification.
#' 
#' \strong{Details:} The fixed classical binary routes expose one explicit starting
#' configuration per learner rather than an implicit search. LDA and QDA use
#' analysis-fold empirical class priors. KNN standardizes within the analysis
#' fold, uses \code{k = 5}, and resolves equal distances by analysis-row order. The
#' four LVQ learners standardize inside each analysis fold and use an equal-prior
#' target for five-neighbor native codebook initialization. Eligibility filtering
#' means the realized class-specific prototype counts need not be equal. Codebook
#' size follows the
#' native default with a minimum of two. LVQ1 uses 100 updates per prototype and
#' alpha 0.03; optimized LVQ1 uses 40 and 0.3. LVQ2.1 uses 100, 0.03, and window
#' 0.3; LVQ3 additionally fixes epsilon 0.1. All use local seed 1. Their declared
#' class-one value is \code{plogis(d0 - d1)}, where \code{d0} and \code{d1} are the nearest
#' squared prototype distances for classes zero and one. This preserves the
#' native nearest-prototype decision boundary but is not a calibrated posterior.
#' Their portable states retain scaling values, patient-derived prototype
#' vectors, and outcome-derived prototype class labels. They contain no patient
#' identifiers or complete analysis table, but a prototype can equal an
#' initialized standardized training row. The
#' Rpart classification tree uses \code{cp = 0.01}, \code{minsplit = 10}, \code{maxdepth = 30}, and
#' no internal cross-validation. The random forest uses 500 trees,
#' \code{floor(sqrt(p))} candidate variables, terminal-node size 1, and seed 1.
#' Linear and radial SVMs use cost 1, analysis-fold scaling, probability output,
#' and seed 1; the radial routes use gamma \code{1 / p}, and the weighted route uses
#' analysis-fold inverse-frequency class weights normalized to mean weight 1.
#' Stepwise logistic regression uses AIC with a declared forward, backward, or
#' both direction. glmBoost uses 100 component-wise boosting iterations with
#' shrinkage 0.1 and no centering. plsRglm requests at most two scaled logistic
#' PLS components. Binary GBM uses 100 trees, depth 1, minimum node size 10,
#' shrinkage 0.1, bag fraction 1, and seed 1.
#' Mean bias-reduced logistic regression uses the \code{brglm2} adjusted-score fit
#' and remains finite under separation when the design is identifiable. DWD
#' uses exponent 1, native analysis-fold penalty selection and scaling; its
#' class-one output is a fixed two-score softmax of the native margin, not a
#' calibrated posterior. K-TSP uses Wilcoxon filtering, disjoint pairs, at most
#' ten candidate pairs, and native K selection inside the analysis fold; its
#' output is the fraction of selected pairs voting for class 1.
#' Fixed LightGBM learners use the native binary objective, 100 rounds,
#' learning rate 0.1, 31 leaves, minimum leaf size 20, one CPU thread,
#' deterministic mode, and seed 1. The GBDT learner uses all rows and features;
#' DART additionally fixes dropout rate 0.1, maximum drops 50, zero skip
#' probability, uniform-drop off, and XGBoost-compatible dropout off. GOSS uses
#' top and other gradient fractions 0.2 and 0.1. The random-forest regime uses
#' row and feature fractions 0.8 and bags every round. These regimes return
#' native class-one probabilities and are not selected or tuned from resampling
#' performance.
#' PLS-DA and sparse PLS-DA use native \code{mixOmics} dummy-response models with an
#' analysis-fold component cap of two and analysis-fold predictor scaling. The
#' effective component count is \code{min(2, p, n - 1)}. Sparse PLS-DA keeps
#' \code{min(10, p)} predictors on each component. These are fixed rules, not tuned
#' values. Their saved state is the exact native linear score projection and
#' contains no fitted patient rows. Because the native two-column scores are not
#' posterior probabilities, the binary output is explicitly defined as the
#' softmax class-one value, \code{plogis(score_1 - score_0)}. It preserves the native
#' greater-score decision boundary but is not labelled a native posterior.
#' These learners differ from plsRglm, which fits a logistic PLS generalized
#' linear model rather than a dummy-matrix PLS discriminant model.
#' Conditional-inference classification trees and forests, bootstrap Rpart
#' bagging, smooth additive and conditional-tree Binomial boosting, localized
#' LDA, and oblique classification forests use fixed controls. Forest and
#' bagging sizes, boosting iterations, neighbor cap, sampling rules, seeds,
#' and single-thread settings are declared before resampling.
#' Earth Binomial adaptive regression splines, polynomial flexible
#' discriminant analysis, mixture discriminant analysis, LIBLINEAR logistic and
#' support-vector formulations, evolutionary trees, analytic shrinkage
#' discriminant analysis, rotation forests, and regularized random forests use
#' one declared starting specification and never search outer or external
#' outcomes.
#' LIBLINEAR exposes L2- and L1-regularized logistic regression, plus
#' L2-regularized L2-loss SVC, L2-regularized L1-loss SVC, and L1-regularized
#' L2-loss SVC. Every formulation fixes cost 1, bias 1, analysis-fold scaling,
#' no cost search, and its documented native termination tolerance. The
#' stochastic native solvers use local seed 1. Logistic learners return native
#' class-one probabilities. SVC learners return the fixed two-class softmax
#' \code{plogis(decision_1 - decision_0)} of native decision values, preserving the
#' native decision boundary without claiming a calibrated posterior. The
#' alternate primal/dual solvers for an already represented objective and the
#' multiclass Crammer-Singer solver are not counted as additional binary
#' learners.
#' Fixed XGBoost gbtree uses 100 rounds, learning rate 0.1, depth 3, child
#' weight 1, full row and column sampling, one thread, and seed 1. The DART
#' learner uses the same trees with uniform tree dropout, tree normalization,
#' dropout rate 0.1, skip probability 0, and no forced dropout. The linear
#' booster uses deterministic coordinate descent in cyclic feature order,
#' learning rate 0.1, L2 penalty 1, L1 penalty 0, one thread, and seed 1. Each
#' booster is exposed separately for binary logistic, Cox, or normal AFT loss
#' because its fitted representation and update algorithm differ; ordinary
#' changes to rounds or tree depth are not counted as new learners. XGBoost AFT
#' fixes distribution scale at 1, represents an event as \verb{[time, time]} and a
#' right-censored observation as \verb{[time, Inf]}, and obtains declared-horizon
#' survival directly from the fitted normal AFT distribution. No route tunes
#' these controls. Mboost AFT crosses Weibull, log-logistic, and log-normal
#' likelihoods with GamBoost spline, GlmBoost component-wise linear, and
#' BlackBoost conditional-tree base learners. Every specification fixes 100
#' boosting iterations, shrinkage 0.1, and nuisance-scale search interval
#' \verb{[0.1, 100]} before resampling. Its risk is the negative native log-time
#' link, and its declared-horizon survival comes directly from the fitted AFT
#' distribution without a Cox baseline. No route searches these controls.
#' Mboost Gehan and Uno-C each cross the same spline, component-wise linear,
#' and conditional-tree base learners with 100 iterations and shrinkage 0.1.
#' Gehan uses a rank-based AFT loss and obtains horizon survival from a
#' Kaplan-Meier estimate of its analysis-fold log-time residual distribution.
#' Uno-C estimates censoring weights inside each analysis fold with fixed
#' sigmoid scale 0.1; its horizon probability uses a one-score Cox calibration
#' and Breslow baseline fitted only in that analysis fold. Both native links
#' increase with longer survival and are negated for the reported risk score.
#' Fixed ranger learners use 500 trees, one thread, seed 1,
#' \code{floor(sqrt(p))} candidate variables, and minimum node size 10 for binary or
#' 3 for survival. Gini and Hellinger probability forests, and logrank, C-index,
#' and maxstat survival forests, use bootstrap sampling. ExtraTrees probability
#' and survival forests instead use the full analysis fold without replacement
#' and one random split per candidate variable. Maxstat fixes \code{alpha = 0.5} and
#' \code{minprop = 0.1}. No ranger route searches these controls.
#' Fixed RLT reinforcement forests use 100 trees, native internal-node split
#' threshold \code{nmin = 5},
#' axis-aligned splits, embedded-model stop size 1, one core, and seed 1; their
#' remaining split, resampling, and muting controls use native defaults. The
#' binary learner returns the native class-one probability. The survival
#' learner converts the native per-event hazard increments to \code{cumprod(1 - h)}
#' survival curves and uses negative left-step restricted mean survival as its
#' ranking score. A survival analysis fold with fewer than two distinct event
#' times is rejected immediately before native fitting because RLT 3.2.6 exits
#' the R process for that one-point grid; an ordinary R condition cannot catch
#' that compiled-code failure.
#' Firth Cox fixes penalty 0.5, Breslow ties, Wald summaries, and 100 maximum
#' iterations. It saves only its unstandardized coefficients and obtains any
#' declared-horizon probability from an analysis-fold Breslow baseline; native
#' convergence warnings remain typed fit failures. Bagged kNN survival uses 50
#' bootstrap learners, all selected features, Mahalanobis distance, constant
#' neighbor weights, \code{k = max(1, floor(n / 10))}, one core, and seed 1. Its
#' direct survival curve is averaged on the native training-time grid and its
#' ranking score is negative left-step restricted mean survival. Because new
#' predictions require nearest-neighbor lookup, its saved state necessarily
#' retains the development predictors, follow-up times, and event indicators.
#' The weighted-subspace binary forest uses 500 trees,
#' \code{floor(log2(p) + 1)} candidate variables, leaf node size 2, serial execution,
#' native vote probabilities, and seed 1. Its saved prediction engine removes
#' training outcomes and patient-level OOB indices.
#' Weibull, log-normal,
#' and log-logistic SurvReg
#' estimate scale in each analysis fold; exponential SurvReg fixes scale at 1.
#' Survival Rpart uses \code{cp = 0.01}, \code{minsplit = 20}, depth 30,
#' and no internal cross-validation. randomForestSRC uses 1,000 logrank trees,
#' node size 5, its package mtry rule, one thread, and seed 1. Ctree uses
#' \code{mincriterion = 0.95}; CForest uses 500 trees, \code{ceiling(sqrt(p))} variables,
#' full replacement sampling, and seed 1. Both partykit learners use minimum
#' split and terminal-node sizes 20 and 7. None of these fixed routes performs
#' an implicit hyperparameter search. The oblique survival forest uses 500
#' fixed-seed aorsf trees; its scalar risk is native mortality and its declared-
#' horizon survival probability is returned directly by the native forest.
#' The Aalen learner fits the full time-varying additive coefficient process and
#' uses its cumulative hazard at the last analysis-fold event as the required
#' scalar ranking score. Its declared-horizon survival is the native
#' \code{exp(-H(t | x))}; additive-hazards estimates can yield values above one or
#' nonmonotone curves, so an invalid requested probability fails explicitly and
#' is never clipped or silently recalibrated. Horizon-enabled final fitting
#' checks that native probability for every development row, without depending
#' on patient order. The proportional-odds learner uses
#' 20 profile iterations and returns \code{x_beta}; its direct survival is
#' \code{1 / (1 + G0(t) * exp(x_beta))}. Both retain the native event-time
#' coefficient grid needed for portable prediction, but no patient IDs,
#' censoring-time vector, event-indicator vector, or fitted engine. They fix
#' local seed 1 for native tie jitter and perform no hyperparameter search.
#' The separate \code{coxboost_cv}, \code{plsrcox_cv}, \code{superpc_cv}, and \code{gbm_cox_cv}
#' routes reproduce their reference internal searches inside each analysis
#' fold; assessment outcomes never enter those searches.
#' ABESS Binomial and Cox learners use a fixed sequence support path from zero
#' through \code{min(10, p)} and select one support size by GIC inside each analysis
#' fold. Their task-specific nonzero-coefficient selectors start at support two;
#' when paired with glmnet, selection is repeated inside every inner analysis
#' fold. ABESS selector-to-ABESS learner chains and combinations with separately
#' tuned survival learners are not admitted.
#' Boruta selectors use ranger normalized permutation importance, 500 trees,
#' one thread, a fixed seed, multiplicity-adjusted \code{p = 0.01} decisions, at most
#' 100 runs, and \code{TentativeRoughFix}. They support binary and survival outcomes.
#' The binary varSelRF selector uses 5,000 trees for initial importance, 2,000
#' trees per elimination iteration, a 20\% drop fraction, and its native one-
#' standard-error OOB subset rule. Forest selectors are refitted inside each
#' analysis fold and, when paired with glmnet, inside every inner analysis fold.
#' The survival Boruta selector is not admitted with separately tuned survival
#' learners whose internal analysis-CV cannot rerun an arbitrary upstream
#' selector.
#' The binary SIS selector uses the native non-iterative marginal-utility
#' ranking and retains \code{min(p, max(2, floor(n / log(n))))} predictors. It uses
#' the pre-penalization \code{sis.ix0} screen rather than the package's downstream
#' penalized support. Binary VSURF uses its native threshold, interpretation,
#' and prediction stages and returns the final prediction subset. Both are
#' refitted inside each outer analysis fold and, with glmnet, inside every inner
#' analysis fold. Binary and survival LPC smooth native two-class t scores or
#' Cox scores across expression-derived principal components, rank by absolute
#' LPC score, and retain the same n/log(n) screen size. LPC's native repeated
#' half-split threshold calculation requires at least 50 candidate predictors.
#' LPC is likewise refitted in every applicable outer and glmnet inner analysis
#' fold. The HDDA learner fits class-specific high-dimensional
#' subspaces with the native Cattell dimension rule, analysis-fold scaling, and
#' class-one posterior probabilities.
#' Native mgcv GAM learners fit one cubic shrinkage smooth (\code{bs = "cs"},
#' \code{k = 5}) per selected numeric predictor. The binary learner uses a logit
#' Binomial model; the survival learner uses mgcv's Cox PH family and its native
#' Peto tie approximation. Basis size and single-thread execution are fixed
#' before resampling, while smoothness is estimated by REML separately inside
#' each analysis fold. There is no outer- or external-outcome tuning.
#' Model-based discriminant learners use native mclust posterior probabilities
#' with empirical analysis-fold class proportions. EDDA fits one Gaussian
#' component per class and selects one common covariance structure by BIC.
#' MclustDA searches one through three components and the covariance structure
#' separately within each class. The univariate grid is \code{E}/\code{V}; the
#' multivariate grid is the 14 native EM covariance models. Any native
#' data-specific reduction of the declared EDDA or MclustDA grid is returned as
#' a typed fit failure. MclustDA hierarchical initialization fixes the native
#' subset cap at 2000, \code{hcModelName = "VVV"}, \code{hcUse = "SVD"}, and seed 1 in a
#' local call; it neither depends on package-wide mclust settings nor changes
#' the caller's random-number state. EDDA uses the known class labels directly
#' and has no stochastic mixture initialization.
#' Penalized discriminant analysis uses native \code{mda::gen.ridge} with lambda 1
#' and one discriminant direction. Stabilised LDA uses \code{ipred::slda} with its
#' native \code{q = NULL} PC-score dimension rule. Both use empirical analysis-fold
#' priors and native posterior probabilities. The classic \code{tree} learner fixes
#' deviance splitting, \code{mincut = 5}, \code{minsize = 10}, and \code{mindev = 0.01}; the
#' C5.0 rule learner fixes ten trials, rule output, no winnowing, disabled early
#' stopping, and seed 1. Neither performs an implicit resampling search.
#' The fixed survival-SVM family crosses the native regression, Van Belle 1,
#' Van Belle 2, and hybrid formulations with linear, additive, and RBF kernels.
#' The
#' two ranking formulations and the hybrid use the native \code{makediff3}
#' comparable-pair method. All twelve standardize within each analysis fold and
#' use the quadratic-programming solver and one declared positive
#' regularization value. For the hybrid model, that scalar is applied equally
#' to its two native regularizers. Linear scores are reconstructed exactly as
#' portable coefficients. Additive and RBF scores are reconstructed exactly
#' from the required scaled kernel basis, native row weights, and intercept;
#' their saved state retains no outcomes or fitted \code{survivalsvm} engine. RBF
#' state also retains the actual native kernel parameter, which follows the
#' median squared distance of standardized analysis predictors inside each fit. Every score is
#' sign-aligned and calibrated to log-hazard scale by a Cox fit using only the
#' same analysis fold. The matching Breslow baseline supplies declared-horizon
#' survival probabilities. No formulation is selected from outer or external
#' outcomes, and no kernel parameter is searched.
#' SLOPE binary regression searches a fixed q grid and coefficient path by
#' native five-fold deviance CV inside the current analysis fold, then saves
#' original-scale logistic coefficients. Because that native CV cannot refit an
#' arbitrary supervised upstream selector, SLOPE is admitted only with all
#' features or a predeclared fixed feature block. The randomGLM binary learner
#' averages 100 bootstrap logistic models over native random feature subspaces,
#' limits each forward search to three candidates, disables interactions, and
#' saves only the per-bag coefficients needed for prediction.
#' The fixed cancerclass learner reproduces its native Welch ranking and
#' correlation-distance \code{z} score over every supplied feature, then estimates
#' one logistic probability mapping from five-fold out-of-fold scores inside
#' the current analysis fold. It is admitted only with all features or a
#' predeclared fixed feature block, so arbitrary learned upstream selectors
#' cannot leak their inner assessment outcomes into that calibration.
#' 
#' Every learner can instead use \code{selector = "fixed_features"} to apply one
#' caller-declared numeric feature block in the stated order. This selector is
#' outcome-independent: it is not refitted or scored in outer or inner folds,
#' and feature-selection stability is not applicable. \code{feature_set_id}
#' distinguishes multiple predeclared blocks that use the same learner. A
#' glmnet learner with \code{alpha > 0} can still report its separate nonzero active-
#' set stability without treating the fixed block as learned.
#' 
#' Nearest-centroid binary classification uses Euclidean class means without
#' scaling or shrinkage and uniform priors. Squared-distance softmax scores
#' are not calibrated risks. Bernoulli NB binarizes at x > 0, with Laplace
#' alpha 1 and training-fold empirical priors; all-positive data lose their
#' continuous information under this rule.
#' }
#' 
#' \subsection{\code{bi_pipeline()}}{
#' Create an ordered modeling route with existing native defaults
#' 
#' \strong{Returns:} An existing-compatible pipeline; no model is fitted.
#' }
#' 
#' \subsection{\code{bi_selector()}}{
#' Internal compatibility constructor for saved binary workflows
#' 
#' \strong{Returns:} An internal selector specification. New code uses \code{\link[=bi_pipeline]{bi_pipeline()}}.
#' }
#' 
#' \subsection{\code{surv_selector()}}{
#' Choose ordered survival feature stages by short name
#' 
#' \strong{Returns:} A survival-stage specification for \code{\link[=surv_pipeline]{surv_pipeline()}}.
#' }
#' 
#' \subsection{\code{sb_lint_pipeline()}}{
#' Checks one existing pipeline object or one named list of arguments for
#' \code{\link[=sb_pipeline]{sb_pipeline()}}. The check reuses the constructor, the pipeline contract,
#' and the built-in adapter resolver; it does not fit a model or inspect data.
#' Expected STORMBREAKER specification errors are returned as \code{rejected}.
#' Unexpected internal errors are not converted and remain visible.
#' 
#' \strong{Returns:} A one-row data frame containing status, pipeline identity, task,
#' selector, learner, and an issue code and message for a rejected
#' specification.
#' 
#' \strong{Details:} \code{admitted} means that the specification has a registered adapter route. It
#' does not check optional-package availability, study design quality, data
#' suitability, literature provenance, or whether fitting will succeed.
#' }
#' }
#' @rawRd \examples{
#' # sb_pipeline
#' sb_pipeline()
#' sb_pipeline(learner = "logistic")
#' sb_pipeline(learner = "gaussian_naive_bayes")
#' sb_pipeline(learner = "lightgbm_binary")
#' sb_pipeline(learner = "lightgbm_binary_dart")
#' sb_pipeline(learner = "catboost_binary")
#' sb_pipeline(learner = "catboost_cox")
#' sb_pipeline(learner = "gbm_cox")
#' sb_pipeline(learner = "ranger_binary")
#' sb_pipeline(learner = "ranger_extratrees_binary")
#' sb_pipeline(learner = "adaboost_binary")
#' sb_pipeline(learner = "logitboost_binary")
#' sb_pipeline(learner = "lda_binary")
#' sb_pipeline(learner = "random_forest_binary")
#' sb_pipeline(learner = "svm_radial_binary")
#' sb_pipeline(learner = "xgboost_binary")
#' sb_pipeline(learner = "xgboost_binary_dart")
#' sb_pipeline(learner = "xgboost_cox_linear")
#' sb_pipeline(learner = "xgboost_aft")
#' sb_pipeline(learner = "ranger_survival")
#' sb_pipeline(learner = "ranger_maxstat_survival")
#' sb_pipeline(learner = "rsf_survival")
#' sb_pipeline(learner = "abess_binary_gic")
#' sb_pipeline(learner = "abess_cox_gic")
#' sb_pipeline(learner = "edda_binary_bic")
#' sb_pipeline(learner = "mclustda_binary_bic")
#' sb_pipeline(learner = "plsda_binary")
#' sb_pipeline(learner = "splsda_binary")
#' sb_pipeline(learner = "survivalsvm_vanbelle1")
#' sb_pipeline(learner = "survivalsvm_vanbelle2")
#' sb_pipeline(learner = "survivalsvm_hybrid")
#' sb_pipeline(learner = "survivalsvm_additive")
#' sb_pipeline(learner = "survivalsvm_hybrid_additive")
#' sb_pipeline(learner = "survivalsvm_rbf")
#' sb_pipeline(learner = "survivalsvm_vanbelle1_rbf")
#' sb_pipeline(learner = "survivalsvm_vanbelle2_rbf")
#' sb_pipeline(learner = "survivalsvm_hybrid_rbf")
#' sb_pipeline(
#'   selector = "fixed_features",
#'   learner = "logistic",
#'   selector_features = c("gene_a", "gene_b"),
#'   feature_set_id = "two_gene"
#' )
#' sb_pipeline(
#'   selector = "univ_binary_auc_top_k",
#'   learner = "logistic",
#'   selector_n_features = 2L
#' )
#' sb_pipeline(
#'   selector = "univ_binary_auc_top_k",
#'   learner = "glmnet_logistic",
#'   selector_n_features = 2L,
#'   glmnet_alpha = 0.5,
#'   inner_folds = 3,
#'   inner_seed = 1701
#' )
#' sb_pipeline(
#'   selector = "stepcox",
#'   learner = "coxph",
#'   stepcox_direction = "forward"
#' )
#' sb_pipeline(
#'   learner = "glmnet_cox",
#'   glmnet_alpha = 0.5,
#'   inner_folds = 3,
#'   inner_seed = 1701
#' )
#' 
#' # sb_lint_pipeline
#' sb_lint_pipeline(list(learner = "logistic"))
#' sb_lint_pipeline(list(learner = "glmnet_cox"))
#' }
NULL

# ---- pipeline ----

sb_format_roundtrip_double <- function(value) {
  for (digits in seq_len(17L)) {
    candidate <- sprintf(paste0("%.", digits, "g"), value)
    if (identical(as.numeric(candidate), value)) {
      return(candidate)
    }
  }
  sprintf("%.17g", value)
}

sb_is_expected_pipeline_constructor_error <- function(error) {
  inherits(error, c(
    "stormbreaker_unsupported_selector",
    "stormbreaker_unsupported_learner",
    "stormbreaker_incompatible_binary_selector",
    "stormbreaker_incompatible_survival_selector",
    "stormbreaker_unsafe_nested_selector",
    "stormbreaker_unadmitted_pipeline",
    "stormbreaker_irrelevant_selector_parameter",
    "stormbreaker_missing_selector_parameter",
    "stormbreaker_bad_fixed_features",
    "stormbreaker_bad_feature_set_id",
    "stormbreaker_bad_stepcox_direction",
    "stormbreaker_bad_stepglm_direction",
    "stormbreaker_bad_selector_options",
    "stormbreaker_irrelevant_stepglm_parameter",
    "stormbreaker_irrelevant_inner_parameter",
    "stormbreaker_irrelevant_glmnet_parameter",
    "stormbreaker_missing_inner_parameter",
    "stormbreaker_bad_glmnet_alpha",
    "stormbreaker_bad_glmnet_standardize",
    "stormbreaker_bad_glmnet_lambda_rule",
    "stormbreaker_bad_glmnet_tuning_measure",
    "stormbreaker_glmnet_too_few_selected_features",
    "stormbreaker_bad_integer",
    "stormbreaker_bad_naive_bayes_var_smoothing",
    "stormbreaker_irrelevant_naive_bayes_parameter",
    "stormbreaker_bad_lightgbm_learning_rate",
    "stormbreaker_irrelevant_lightgbm_parameter",
    "stormbreaker_bad_catboost_learning_rate",
    "stormbreaker_bad_catboost_depth",
    "stormbreaker_bad_catboost_l2_leaf_reg",
    "stormbreaker_irrelevant_catboost_parameter",
    "stormbreaker_bad_gbm_shrinkage",
    "stormbreaker_bad_gbm_bag_fraction",
    "stormbreaker_irrelevant_gbm_parameter",
    "stormbreaker_irrelevant_ranger_parameter",
    "stormbreaker_bad_coxboost_penalty",
    "stormbreaker_bad_superpc_threshold",
    "stormbreaker_bad_superpc_components",
    "stormbreaker_bad_survivalsvm_gamma_mu",
    "stormbreaker_irrelevant_classical_survival_parameter",
    "stormbreaker_irrelevant_classical_binary_boosting_parameter"
  ))
}

sb_pipeline_parameters_equal <- function(x, y) {
  if (is.list(x) && is.list(y) && !is.object(x) && !is.object(y)) {
    if (length(x) != length(y)) return(FALSE)
    nx <- names(x)
    ny <- names(y)
    if (!is.null(nx) && !is.null(ny) && all(nzchar(nx)) &&
      !anyDuplicated(nx) && !anyDuplicated(ny) && setequal(nx, ny)) {
      y[] <- y[match(nx, ny)]
      names(y) <- nx
    }
    return(identical(attributes(x), attributes(y)) && all(vapply(
      seq_along(x), function(i) sb_pipeline_parameters_equal(x[[i]], y[[i]]),
      logical(1)
    )))
  }
  # Integer storage is not a different model setting; values remain exact.
  if (is.integer(x) && !is.object(x)) storage.mode(x) <- "double"
  if (is.integer(y) && !is.object(y)) storage.mode(y) <- "double"
  identical(x, y)
}

sb_pipeline_execution_fields <- function(task_type) {
  c("id", "task_type", "selector", "selector_params", "learner",
    "learner_params", "fallback",
    if (identical(task_type, "binary")) {
      c("prediction_type", "positive_class", "prediction_direction")
    } else {
      c("risk_type", "risk_reference", "risk_formula", "risk_direction")
    })
}

sb_validate_pipeline_contract <- function(pipeline, task_type) {
  pipeline <- sb_upgrade_legacy_names(pipeline)
  if (is.list(pipeline)) sb_reject_retired_learner(pipeline$learner)

  if (!inherits(pipeline, "stormbreaker_pipeline")) {
    sb_abort(
      "`pipeline` must be created by `sb_pipeline()`.",
      "stormbreaker_bad_pipeline"
    )
  }
  if (!task_type %in% c("survival", "binary") ||
    !identical(pipeline[["task_type"]], task_type)) {
    sb_abort(
      if (identical(task_type, "binary")) {
        "A binary task requires a binary pipeline."
      } else {
        "A survival task requires a survival pipeline."
      },
      "stormbreaker_pipeline_task_mismatch"
    )
  }
  if (identical(pipeline[["selector"]], "ordered_chain")) {
    return(sb_validate_chain_pipeline(pipeline, task_type))
  }

  selector <- pipeline[["selector"]]
  learner <- pipeline[["learner"]]
  uses_inner_cv <- sb_pipeline_uses_inner_cv(learner)
  is_catboost <- learner %in% c(
    "catboost_binary",
    "catboost_cox",
    "catboost_aft_survival"
  )
  is_ranger <- learner %in% sb_ranger_learner_ids()
  is_gbm_fixed <- identical(learner, "gbm_cox")
  is_gbm <- learner %in% c("gbm_cox", "gbm_cox_cv")
  is_analysis_tuned <- learner %in% sb_analysis_tuned_survival_learner_ids()
  is_coxboost <- identical(learner, "coxboost")
  is_plsrcox <- identical(learner, "plsrcox")
  is_superpc <- identical(learner, "superpc")
  is_survivalsvm <- learner %in% sb_survivalsvm_learner_ids()
  selector_params <- if ("selector_params" %in% names(pipeline)) {
    pipeline[["selector_params"]]
  } else {
    NULL
  }
  learner_params <- if ("learner_params" %in% names(pipeline)) {
    pipeline[["learner_params"]]
  } else {
    NULL
  }
  parameter <- function(name) {
    if (!is.list(learner_params) || !name %in% names(learner_params)) {
      return(NULL)
    }
    learner_params[[name]]
  }
  reconstructed <- tryCatch(
    sb_pipeline(
      selector = selector,
      learner = learner,
      selector_options = if (selector %in% c(sb_binary_model_selector_ids(),
        sb_boruta_selector_ids())) selector_params else NULL,
      selector_n_features = if (sb_selector_is_learned(selector) &&
        is.list(selector_params) &&
        "n_features" %in% names(selector_params)) {
        selector_params[["n_features"]]
      } else {
        NULL
      },
      selector_features = if (identical(selector, "fixed_features") &&
        is.list(selector_params)) {
        selector_params[["features"]]
      } else {
        NULL
      },
      feature_set_id = if (identical(selector, "fixed_features") &&
        is.list(selector_params)) {
        selector_params[["feature_set_id"]]
      } else {
        NULL
      },
      stepcox_direction = if (selector %in% c(
        "stepcox",
        sb_mime_composite_selector_ids()
      ) &&
        is.list(selector_params)) {
        selector_params[["direction"]]
      } else {
        NULL
      },
      selector_inner_folds = if (selector %in%
        sb_mime_tuned_selector_ids() && is.list(selector_params)) {
        selector_params[["inner_folds"]]
      } else {
        NULL
      },
      selector_seed = if (selector %in% sb_mime_survival_selector_ids() &&
        is.list(selector_params)) {
        selector_params[["seed"]]
      } else {
        NULL
      },
      selector_rsf_nodesize = if (selector %in%
        sb_mime_rsf_selector_ids() && is.list(selector_params)) {
        selector_params[["nodesize"]]
      } else {
        NULL
      },
      learner_tuning_folds = if (is_analysis_tuned) {
        parameter("inner_folds")
      } else {
        NULL
      },
      learner_seed = if (is_analysis_tuned) {
        parameter("seed")
      } else {
        NULL
      },
      gbm_cv_max_trees = if (identical(learner, "gbm_cox_cv")) {
        parameter("maximum_trees")
      } else {
        NULL
      },
      inner_folds = if (uses_inner_cv) {
        parameter("inner_folds")
      } else {
        NULL
      },
      inner_seed = if (uses_inner_cv) {
        parameter("inner_seed")
      } else {
        NULL
      },
      glmnet_alpha = if (uses_inner_cv) {
        parameter("alpha")
      } else {
        NULL
      },
      glmnet_standardize = if (uses_inner_cv) {
        parameter("standardize")
      } else {
        NULL
      },
      glmnet_lambda_rule = if (uses_inner_cv) {
        parameter("selection_rule")
      } else {
        NULL
      },
      glmnet_tuning_measure = if (uses_inner_cv) {
        parameter("measure")
      } else {
        NULL
      },
      naive_bayes_var_smoothing = if (identical(
        learner,
        "gaussian_naive_bayes"
      )) {
        parameter("var_smoothing")
      } else {
        NULL
      },
      lightgbm_nrounds = if (
        learner %in% sb_lightgbm_binary_learner_ids()
      ) {
        parameter("nrounds")
      } else {
        NULL
      },
      lightgbm_learning_rate = if (
        learner %in% sb_lightgbm_binary_learner_ids()
      ) {
        parameter("learning_rate")
      } else {
        NULL
      },
      lightgbm_num_leaves = if (
        learner %in% sb_lightgbm_binary_learner_ids()
      ) {
        parameter("num_leaves")
      } else {
        NULL
      },
      lightgbm_min_data_in_leaf = if (
        learner %in% sb_lightgbm_binary_learner_ids()
      ) {
        parameter("min_data_in_leaf")
      } else {
        NULL
      },
      lightgbm_seed = if (
        learner %in% sb_lightgbm_binary_learner_ids()
      ) {
        parameter("seed")
      } else {
        NULL
      },
      catboost_iterations = if (is_catboost) {
        parameter("iterations")
      } else {
        NULL
      },
      catboost_learning_rate = if (is_catboost) {
        parameter("learning_rate")
      } else {
        NULL
      },
      catboost_depth = if (is_catboost) {
        parameter("depth")
      } else {
        NULL
      },
      catboost_l2_leaf_reg = if (is_catboost) {
        parameter("l2_leaf_reg")
      } else {
        NULL
      },
      catboost_seed = if (is_catboost) {
        parameter("random_seed")
      } else {
        NULL
      },
      gbm_n_trees = if (is_gbm_fixed) {
        parameter("n_trees")
      } else {
        NULL
      },
      gbm_interaction_depth = if (is_gbm_fixed) {
        parameter("interaction_depth")
      } else {
        NULL
      },
      gbm_min_obs_node = if (is_gbm_fixed) {
        parameter("min_obs_node")
      } else {
        NULL
      },
      gbm_shrinkage = if (is_gbm_fixed) {
        parameter("shrinkage")
      } else {
        NULL
      },
      gbm_bag_fraction = if (is_gbm_fixed) {
        parameter("bag_fraction")
      } else {
        NULL
      },
      gbm_seed = if (is_gbm_fixed) {
        parameter("seed")
      } else {
        NULL
      },
      ranger_num_trees = if (is_ranger) {
        parameter("num_trees")
      } else {
        NULL
      },
      ranger_mtry = if (is_ranger) {
        parameter("mtry")
      } else {
        NULL
      },
      ranger_min_node_size = if (is_ranger) {
        parameter("min_node_size")
      } else {
        NULL
      },
      ranger_seed = if (is_ranger) {
        parameter("seed")
      } else {
        NULL
      },
      coxboost_stepno = if (is_coxboost) {
        parameter("stepno")
      } else {
        NULL
      },
      coxboost_penalty = if (is_coxboost) {
        parameter("penalty")
      } else {
        NULL
      },
      plsrcox_components = if (is_plsrcox) {
        parameter("components")
      } else {
        NULL
      },
      superpc_threshold = if (is_superpc) {
        parameter("threshold")
      } else {
        NULL
      },
      superpc_components = if (is_superpc) {
        parameter("components")
      } else {
        NULL
      },
      survivalsvm_gamma_mu = if (is_survivalsvm) {
        parameter("gamma_mu")
      } else {
        NULL
      },
      adaboost_iterations = if (identical(
        learner,
        "adaboost_binary"
      )) {
        parameter("iterations")
      } else {
        NULL
      },
      logitboost_iterations = if (identical(
        learner,
        "logitboost_binary"
      )) {
        parameter("iterations")
      } else {
        NULL
      },
      stepglm_direction = if (identical(
        learner,
        "stepglm_binary"
      )) {
        parameter("direction")
      } else {
        NULL
      }
    ),
    stormbreaker_error = function(error) {
      if (!sb_is_expected_pipeline_constructor_error(error)) {
        stop(error)
      }
      NULL
    }
  )

  # Preserve the recorded iteration budget in pipelines saved before 0.3.34.
  if (!is.null(reconstructed) && identical(learner, "nnet_binary") &&
    sb_pipeline_parameters_equal(learner_params$maxit, 200L)) reconstructed$learner_params$maxit <- 200L

  # Historical pipelines retain their declared model, not today's defaults.
  if (!is.null(reconstructed)) {
    if (identical(learner, "neuralnet_multilayer_binary") &&
      sb_pipeline_parameters_equal(learner_params$stepmax, 100000L)) {
      reconstructed$learner_params$stepmax <- 100000L
    }
    if (identical(learner, "neuralnet_multilayer_binary") &&
      sb_pipeline_parameters_equal(learner_params$hidden, c(5L, 3L)) &&
      identical(learner_params$error_function, "ce")) {
      reconstructed$learner_params$hidden <- c(5L, 3L)
      reconstructed$learner_params$error_function <- "ce"
    }
    if (identical(learner, "fda_binary") && sb_pipeline_parameters_equal(learner_params$degree, 2L)) {
      reconstructed$learner_params$degree <- 2L
    }
    if (identical(learner, "loclda_binary") && is.null(learner_params$k_rule)) {
      reconstructed$learner_params$k_rule <- NULL
    }
    if (learner %in% sb_survivalsvm_rbf_learner_ids() &&
      identical(learner_params$kernel_parameter_rule, "native_inverse_analysis_feature_count")) {
      reconstructed$learner_params$kernel_parameter_rule <- "native_inverse_analysis_feature_count"
    }
  }

  contract_fields <- sb_pipeline_execution_fields(task_type)
  contract_ok <- !is.null(reconstructed) && all(vapply(
    contract_fields,
    function(field) {
      if (field %in% c("selector_params", "learner_params")) {
        sb_pipeline_parameters_equal(pipeline[[field]], reconstructed[[field]])
      } else {
        identical(pipeline[[field]], reconstructed[[field]])
      }
    },
    logical(1)
  ))
  if (!contract_ok) {
    sb_abort(
      if (identical(task_type, "binary")) {
        "The supplied binary pipeline is not executable by the minimum core."
      } else {
        "The supplied pipeline is not executable by the minimum core."
      },
      "stormbreaker_unresolved_pipeline"
    )
  }
  pipeline$selector_params <- reconstructed$selector_params
  pipeline$learner_params <- reconstructed$learner_params
  pipeline
}

#' @rdname sb_pipeline
#' @export
sb_pipeline <- function(selector = "none", learner = "coxph",
                        selector_n_features = NULL,
                        inner_folds = NULL, inner_seed = NULL,
                        glmnet_alpha = NULL,
                        naive_bayes_var_smoothing = NULL,
                        lightgbm_nrounds = NULL,
                        lightgbm_learning_rate = NULL,
                        lightgbm_num_leaves = NULL,
                        lightgbm_min_data_in_leaf = NULL,
                        lightgbm_seed = NULL,
                        catboost_iterations = NULL,
                        catboost_learning_rate = NULL,
                        catboost_depth = NULL,
                        catboost_l2_leaf_reg = NULL,
                        catboost_seed = NULL,
                        ranger_num_trees = NULL,
                        ranger_mtry = NULL,
                        ranger_min_node_size = NULL,
                        ranger_seed = NULL,
                        selector_features = NULL,
                        feature_set_id = NULL,
                        glmnet_standardize = NULL,
                        glmnet_lambda_rule = NULL,
                        glmnet_tuning_measure = NULL,
                        stepcox_direction = NULL,
                        gbm_n_trees = NULL,
                        gbm_interaction_depth = NULL,
                        gbm_min_obs_node = NULL,
                        gbm_shrinkage = NULL,
                        gbm_bag_fraction = NULL,
                        gbm_seed = NULL,
                        coxboost_stepno = NULL,
                        coxboost_penalty = NULL,
                        plsrcox_components = NULL,
                        superpc_threshold = NULL,
                        superpc_components = NULL,
                        survivalsvm_gamma_mu = NULL,
                        adaboost_iterations = NULL,
                        logitboost_iterations = NULL,
                        stepglm_direction = NULL,
                        selector_inner_folds = NULL,
                        selector_seed = NULL,
                        selector_rsf_nodesize = NULL,
                        learner_tuning_folds = NULL,
                        learner_seed = NULL,
                        gbm_cv_max_trees = NULL, selector_options = NULL) {
  learner <- sb_canonical_learner_id(learner)
  sb_reject_retired_learner(learner)

  registry <- sb_pipeline_registry_data()
  if (!is.character(selector) || length(selector) != 1L ||
    is.na(selector) || !selector %in% unique(registry$selector)) {
    sb_abort(
      "The requested pipeline selector is not registered.",
      "stormbreaker_unsupported_selector"
    )
  }
  if (!is.character(learner) || length(learner) != 1L ||
    is.na(learner) || !learner %in% unique(registry$learner)) {
    sb_abort(
      "The requested pipeline learner is not registered.",
      "stormbreaker_unsupported_learner"
    )
  }
  learner <- unname(as.character(learner))
  admission <- sb_pipeline_registry_match(selector, learner, registry)
  learner_rows <- registry[registry$learner == learner, , drop = FALSE]
  learner_task_type <- unique(learner_rows$task_type)
  selector_task_types <- unique(registry$task_type[
    registry$selector == selector
  ])
  uses_inner_cv <- sb_pipeline_uses_inner_cv(learner, registry)
  if (nrow(admission) == 0L &&
    !learner_task_type %in% selector_task_types) {
    if (identical(learner_task_type, "binary")) {
      sb_abort(
        paste0(
          "Binary pipelines cannot use the survival-specific ",
          "`", selector, "` selector."
        ),
        "stormbreaker_incompatible_binary_selector"
      )
    }
    sb_abort(
      paste0(
        "Survival pipelines cannot use the binary-specific ",
        "`univ_binary_auc_top_k` selector."
      ),
      "stormbreaker_incompatible_survival_selector"
    )
  }
  if (nrow(admission) == 0L && uses_inner_cv &&
    !identical(selector, "none")) {
    sb_abort(
      paste0(
        "This nested selector and glmnet combination is not admitted. ",
        "Only combinations that refit selection inside every inner ",
        "analysis fold are safe."
      ),
      "stormbreaker_unsafe_nested_selector"
    )
  }
  if (nrow(admission) != 1L) {
    sb_abort(
      paste0(
        "The selector and learner combination `", selector, " -> ",
        learner, "` is not admitted."
      ),
      "stormbreaker_unadmitted_pipeline"
    )
  }
  pipeline_task_type <- admission$task_type[[1L]]
  options <- list(
    selector_n_features = selector_n_features,
    inner_folds = inner_folds,
    inner_seed = inner_seed,
    glmnet_alpha = glmnet_alpha,
    naive_bayes_var_smoothing = naive_bayes_var_smoothing,
    lightgbm_nrounds = lightgbm_nrounds,
    lightgbm_learning_rate = lightgbm_learning_rate,
    lightgbm_num_leaves = lightgbm_num_leaves,
    lightgbm_min_data_in_leaf = lightgbm_min_data_in_leaf,
    lightgbm_seed = lightgbm_seed,
    catboost_iterations = catboost_iterations,
    catboost_learning_rate = catboost_learning_rate,
    catboost_depth = catboost_depth,
    catboost_l2_leaf_reg = catboost_l2_leaf_reg,
    catboost_seed = catboost_seed,
    ranger_num_trees = ranger_num_trees,
    ranger_mtry = ranger_mtry,
    ranger_min_node_size = ranger_min_node_size,
    ranger_seed = ranger_seed,
    selector_features = selector_features,
    feature_set_id = feature_set_id,
    glmnet_standardize = glmnet_standardize,
    glmnet_lambda_rule = glmnet_lambda_rule,
    glmnet_tuning_measure = glmnet_tuning_measure,
    stepcox_direction = stepcox_direction,
    gbm_n_trees = gbm_n_trees,
    gbm_interaction_depth = gbm_interaction_depth,
    gbm_min_obs_node = gbm_min_obs_node,
    gbm_shrinkage = gbm_shrinkage,
    gbm_bag_fraction = gbm_bag_fraction,
    gbm_seed = gbm_seed,
    coxboost_stepno = coxboost_stepno,
    coxboost_penalty = coxboost_penalty,
    plsrcox_components = plsrcox_components,
    superpc_threshold = superpc_threshold,
    superpc_components = superpc_components,
    survivalsvm_gamma_mu = survivalsvm_gamma_mu,
    adaboost_iterations = adaboost_iterations,
    logitboost_iterations = logitboost_iterations,
    stepglm_direction = stepglm_direction,
    selector_inner_folds = selector_inner_folds,
    selector_seed = selector_seed,
    selector_rsf_nodesize = selector_rsf_nodesize,
    learner_tuning_folds = learner_tuning_folds,
    learner_seed = learner_seed,
    gbm_cv_max_trees = gbm_cv_max_trees,
    selector_options = selector_options
  )
  selection <- sb_pipeline_selector_config(selector, admission, uses_inner_cv, options)
  selector_params <- selection$params
  learner_params <- sb_pipeline_tuning_params(learner, pipeline_task_type, uses_inner_cv, options)
  learner_params <- sb_pipeline_tree_params(learner, learner_params, options)
  learner_params <- sb_pipeline_fixed_params(learner, learner_params, options)
  base_id <- paste0(selector, "__", learner, selection$suffix)
  pipeline_id <- sb_pipeline_model_id(base_id, learner, learner_params, uses_inner_cv)
  pipeline <- c(list(
    id = pipeline_id,
    task_type = pipeline_task_type,
    selector = selector,
    learner = learner
  ), sb_pipeline_prediction_spec(pipeline_task_type, learner, learner_params))
  pipeline$fallback <- "forbidden"
  if (length(selector_params) > 0L) {
    pipeline$selector_params <- selector_params
  }
  if (length(learner_params) > 0L) {
    pipeline$learner_params <- learner_params
  }
  class(pipeline) <- "stormbreaker_pipeline"
  pipeline
}

#' @export
print.stormbreaker_pipeline <- function(x, ...) {
  cat("<stormbreaker_pipeline:", x$id, ">\n")
  cat("Task:", x$task_type, "\n")
  cat("Selector:", x$selector, "| Learner:", x$learner, "\n")
  if (!is.null(x$selector_params)) {
    if (identical(x$selector, "fixed_features")) {
      cat(
        "Predeclared feature set:", x$selector_params$feature_set_id,
        "| Features:", x$selector_params$n_features, "\n"
      )
    } else if (identical(x$selector, "stepcox")) {
      cat("StepCox direction:", x$selector_params$direction, "\n")
    } else if (x$selector %in% sb_mime_rsf_selector_ids()) {
      cat(
        "RSF selector: 1000 trees | node size",
        x$selector_params$nodesize,
        "| seed", x$selector_params$seed,
        if (x$selector %in% sb_mime_composite_selector_ids()) {
          paste0(" | StepCox ", x$selector_params$direction)
        } else {
          ""
        },
        "\n"
      )
    } else if (x$selector %in% sb_mime_tuned_selector_ids()) {
      cat(
        "Tuned selector:", x$selector,
        "| inner folds", x$selector_params$inner_folds,
        "| seed", x$selector_params$seed,
        if (x$selector %in% sb_mime_composite_selector_ids()) {
          paste0(" | StepCox ", x$selector_params$direction)
        } else {
          ""
        },
        "\n"
      )
    } else {
      if (x$selector %in% sb_binary_model_selector_ids()) {
        cat("Selector options:", paste(names(x$selector_params),
          unlist(x$selector_params), sep = "=", collapse = ", "), "\n")
      } else {
        cat("Selector features:", x$selector_params$n_features, "\n")
      }
    }
  }
  if (!is.null(x$learner_params)) {
    if (identical(x$learner, "gaussian_naive_bayes")) {
      cat(
        "Variance smoothing:",
        format(
          x$learner_params$var_smoothing,
          trim = TRUE,
          scientific = TRUE,
          digits = 15L
        ),
        "\n"
      )
    } else if (x$learner %in% sb_lightgbm_binary_learner_ids()) {
      cat(
        if (identical(x$learner, "lightgbm_binary")) {
          "Fixed LightGBM:"
        } else {
          paste0(
            "Fixed LightGBM ",
            toupper(sb_lightgbm_binary_variant(x$learner)),
            ":"
          )
        },
        x$learner_params$nrounds, "rounds |",
        "learning rate", sb_format_roundtrip_double(
          x$learner_params$learning_rate
        ), "| leaves", x$learner_params$num_leaves,
        "| min leaf", x$learner_params$min_data_in_leaf,
        "| seed", x$learner_params$seed, "\n"
      )
    } else if (x$learner %in% c(
      "catboost_binary",
      "catboost_cox",
      "catboost_aft_survival"
    )) {
      cat(
        "Fixed CatBoost:", x$learner_params$iterations, "iterations |",
        "learning rate", sb_format_roundtrip_double(
          x$learner_params$learning_rate
        ), "| depth", x$learner_params$depth,
        "| L2 leaf", sb_format_roundtrip_double(
          x$learner_params$l2_leaf_reg
        ), "| seed", x$learner_params$random_seed, "\n"
      )
    } else if (x$learner %in% sb_analysis_tuned_survival_learner_ids()) {
      cat(
        "Analysis-fold tuned learner:", x$learner,
        "| inner folds", x$learner_params$inner_folds,
        "| seed", x$learner_params$seed,
        if (identical(x$learner, "gbm_cox_cv")) {
          paste0(" | maximum trees ", x$learner_params$maximum_trees)
        } else {
          ""
        },
        "\n"
      )
    } else if (identical(x$learner, "gbm_cox")) {
      cat(
        "Fixed GBM-Cox:", x$learner_params$n_trees, "trees |",
        "depth", x$learner_params$interaction_depth,
        "| min node", x$learner_params$min_obs_node,
        "| shrinkage", sb_format_roundtrip_double(
          x$learner_params$shrinkage
        ), "| bag fraction", sb_format_roundtrip_double(
          x$learner_params$bag_fraction
        ), "| seed", x$learner_params$seed, "\n"
      )
    } else if (x$learner %in% sb_mboost_cox_learner_ids()) {
      cat(
        "Fixed mboost Cox:", x$learner_params$mstop,
        "iterations | shrinkage",
        sb_format_roundtrip_double(x$learner_params$nu),
        "| base learner", x$learner_params$base_learner, "\n"
      )
    } else if (x$learner %in% sb_mboost_aft_learner_ids()) {
      cat(
        "Fixed mboost AFT:", x$learner_params$mstop,
        "iterations |", x$learner_params$distribution,
        "| base learner", x$learner_params$base_learner, "\n"
      )
    } else if (x$learner %in% sb_mboost_ranking_learner_ids()) {
      cat(
        "Fixed mboost ranking:", x$learner_params$mstop,
        "iterations |", x$learner_params$objective,
        "| base learner", x$learner_params$base_learner, "\n"
      )
    } else if (x$learner %in% sb_rlt_learner_ids()) {
      cat(
        "Fixed RLT reinforcement forest:",
        x$learner_params$num_trees, "trees | native nmin",
        x$learner_params$native_nmin,
        "| axis-aligned split | seed", x$learner_params$seed, "\n"
      )
    } else if (x$learner %in% sb_coxphf_learner_ids()) {
      cat(
        "Fixed Firth Cox: penalty", x$learner_params$penalty,
        "| Breslow ties | maximum iterations",
        x$learner_params$max_iterations, "\n"
      )
    } else if (x$learner %in% sb_akritas_learner_ids()) {
      cat("Native Akritas: fixed empirical-CDF bandwidth", x$learner_params$lambda,
        "| negative step-RMST | conditional product-limit survival\n")
    } else if (x$learner %in% sb_bagged_knn_survival_learner_ids()) {
      cat(
        "Fixed bagged kNN survival:",
        x$learner_params$num_base_learners,
        "bootstrap learners | k=max(1,floor(n/10)) | seed",
        x$learner_params$seed, "\n"
      )
    } else if (x$learner %in% sb_ahaz_survival_learner_ids()) {
      cat("Additive hazards:", sb_ahaz_parameter_semantics(x$learner), "\n")
    } else if (x$learner %in% sb_bujar_linear_learner_ids()) {
      cat(
        "Fixed Buckley-James linear boosting:",
        x$learner_params$iter_bj, "BJ iterations |",
        x$learner_params$mstop, "boosting steps | seed",
        x$learner_params$seed, "\n"
      )
    } else if (x$learner %in% sb_wsrf_binary_learner_ids()) {
      cat(
        "Fixed weighted-subspace forest:",
        x$learner_params$num_trees,
        "trees | mtry=floor(log2(p)+1) | seed",
        x$learner_params$seed, "\n"
      )
    } else if (x$learner %in% sb_ranger_learner_ids()) {
      cat(
        if (x$learner %in% sb_ranger_survival_learner_ids()) {
          "Fixed ranger survival forest:"
        } else {
          "Fixed ranger binary probability forest:"
        },
        x$learner_params$num_trees, "trees | mtry",
        if (is.null(x$learner_params$mtry)) {
          "floor(sqrt(p))"
        } else {
          x$learner_params$mtry
        }, "| min node", x$learner_params$min_node_size,
        "| split rule", x$learner_params$splitrule,
        "| seed", x$learner_params$seed, "\n"
      )
    } else if (identical(x$learner, "coxboost")) {
      cat(
        "Fixed CoxBoost:", x$learner_params$stepno, "steps | penalty",
        if (is.null(x$learner_params$penalty)) {
          "9 x analysis events"
        } else {
          sb_format_roundtrip_double(x$learner_params$penalty)
        }, "\n"
      )
    } else if (identical(x$learner, "plsrcox")) {
      cat("Fixed plsRcox components:", x$learner_params$components, "\n")
    } else if (identical(x$learner, "superpc")) {
      cat(
        "Fixed SuperPC: threshold",
        sb_format_roundtrip_double(x$learner_params$threshold),
        "| components", x$learner_params$components, "\n"
      )
    } else if (x$learner %in% sb_pamr_survival_learner_ids()) {
      cat("Fixed PAMR survival: two censoring-aware groups | threshold 1\n")
    } else if (x$learner %in% sb_survivalsvm_learner_ids()) {
      kernel_label <- switch(x$learner_params$kernel,
        lin_kernel = "linear",
        add_kernel = "additive",
        rbf_kernel = "RBF"
      )
      cat(
        if (identical(x$learner_params$type, "regression")) {
          paste0(
            "Fixed ", kernel_label, " survival-SVM gamma:"
          )
        } else {
          paste0(
            "Fixed ", kernel_label, " survival-SVM ",
            x$learner_params$type,
            " gamma:"
          )
        },
        sb_format_roundtrip_double(x$learner_params$gamma_mu), "\n"
      )
    } else if (x$learner %in%
      sb_classical_binary_boosting_learner_ids()) {
      cat(
        if (identical(x$learner, "adaboost_binary")) {
          "Fixed AdaBoost.M1:"
        } else {
          "Fixed LogitBoost:"
        },
        x$learner_params$iterations,
        if (identical(x$learner, "adaboost_binary")) {
          "maximum iterations | base learner"
        } else {
          "iterations | base learner"
        },
        x$learner_params$base_learner,
        "\n"
      )
    } else if (identical(x$learner, "stepglm_binary")) {
      cat(
        "Fixed stepwise logistic: direction",
        x$learner_params$direction,
        "| criterion AIC\n"
      )
    } else if (identical(x$learner, "glmboost_binary")) {
      cat(
        "Fixed glmBoost:", x$learner_params$iterations,
        "iterations | shrinkage",
        sb_format_roundtrip_double(x$learner_params$shrinkage), "\n"
      )
    } else if (identical(x$learner, "plsrglm_binary")) {
      cat(
        "Fixed plsRglm components:",
        x$learner_params$components,
        "\n"
      )
    } else if (x$learner %in% sb_mgcv_gam_learner_ids()) {
      cat(
        "GAM cubic shrinkage basis: k=",
        x$learner_params$basis_dimension,
        " | smoothness estimated by REML in each analysis fold\n",
        sep = ""
      )
    } else if (x$learner %in% sb_mclust_da_learner_ids()) {
      cat(
        "Model-based discriminant analysis:",
        x$learner_params$model_type,
        "| covariance/components selected by analysis-fold BIC\n"
      )
    } else if (x$learner %in% sb_mixomics_da_learner_ids()) {
      cat(
        "Fixed", x$learner_params$method,
        "| component cap", x$learner_params$requested_components,
        "| analysis-fold scaling | fixed binary softmax mapping\n"
      )
    } else if (x$learner %in% sb_hdda_binary_learner_ids()) {
      cat(
        "HDDA class-specific subspaces | Cattell dimension rule |",
        "analysis-fold scaling | native class-one posterior\n"
      )
    } else if (x$learner %in% sb_ncvreg_learner_ids()) {
      cat(
        x$learner_params$penalty,
        "penalty | canonical gamma",
        sb_format_roundtrip_double(x$learner_params$gamma),
        "| lambda selected by analysis-fold BIC\n"
      )
    } else if (x$learner %in% sb_xgboost_learner_ids()) {
      cat(
        "Fixed XGBoost booster:",
        sb_xgboost_booster(x$learner),
        "| objective", x$learner_params$objective,
        "| rounds", x$learner_params$nrounds,
        "| one thread | seed", x$learner_params$seed,
        "\n"
      )
    } else if (x$learner %in% c(
      sb_classical_binary_learner_ids(),
      sb_additional_survival_learner_ids(),
      sb_flexsurv_parametric_learner_ids(),
      sb_flexsurvspline_learner_ids(),
      sb_timereg_survival_learner_ids(),
      sb_rule_aft_learner_ids()
    )) {
      cat("Learner controls: fixed before resampling\n")
    } else if (x$learner %in% c("glmnet_cox", "glmnet_logistic")) {
      penalty_family <- if (identical(x$learner_params$alpha, 0)) {
        "Ridge"
      } else if (identical(x$learner_params$alpha, 1)) {
        "Lasso"
      } else {
        "Elastic Net"
      }
      cat("Penalty family:", penalty_family, "\n")
      cat(
        "Alpha:", x$learner_params$alpha,
        "| Inner folds:", x$learner_params$inner_folds,
        "| Lambda rule:", x$learner_params$selection_rule, "\n"
      )
      cat(
        "Glmnet standardization:",
        if (isTRUE(x$learner_params$standardize)) "enabled" else "disabled",
        "\n"
      )
    } else {
      cat("Learner controls: fixed before resampling\n")
    }
  }
  if (identical(x$task_type, "binary")) {
    cat(
      "Prediction:", x$prediction_type,
      "of class", x$positive_class, "\n"
    )
  } else {
    cat("Risk direction:", x$risk_direction, "\n")
  }
  cat("Fallback:", x$fallback, "\n")
  invisible(x)
}
