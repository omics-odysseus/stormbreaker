# ---- zzz-short-api ----

#' @rdname sb_models
#' @export
models <- sb_models

#' @rdname sb_binary_task
#' @export
bi_task <- sb_binary_task

#' @rdname sb_binary_task
#' @export
surv_task <- sb_surv_task

#' @rdname sb_binary_task
#' @export
multi_task <- sb_multiclass_task

#' @rdname sb_binary_task
#' @export
bi_endpoint <- sb_binary_endpoint

#' @rdname sb_binary_task
#' @export
surv_endpoint <- sb_surv_endpoint

#' @rdname sb_compare_cv
#' @export
run_cv <- sb_run_cv

#' @rdname sb_compare_cv
#' @export
compare_cv <- sb_compare_cv

#' @rdname sb_fit_final
#' @export
fit_final <- sb_fit_final

#' @rdname sb_fit_final
#' @export
predict_final <- sb_predict_final

#' @rdname sb_evaluate_external
#' @export
eval_external <- sb_evaluate_external

#' @rdname sb_compare_external
#' @export
compare_external <- sb_compare_external

#' @rdname sb_complete_run
#' @export
finish_run <- sb_complete_run

#' @rdname sb_plot_calibration
#' @export
clinical_plots <- sb_clinician_views

#' @rdname sb_native_diagnostics
#' @export
plots <- sb_plots

#' @rdname sb_render_report
#' @export
result_table <- sb_result_table

#' @rdname sb_analyze_methods
#' @export
selected_genes <- sb_selected_features

#' @rdname sb_render_report
#' @export
export_run <- sb_export_run

#' @rdname sb_render_report
#' @export
clinical_report <- sb_render_clinician_report

#' @rdname sb_model_shap
#' @export
model_shap <- sb_model_shap

#' @rdname sb_model_coefficients
#' @export
model_coefs <- sb_model_coefficients

#' @rdname sb_binary_task
#' @export
bi_example <- sb_example_binary

#' @rdname sb_binary_task
#' @export
surv_example <- sb_example_survival

#' @rdname sb_plot_performance
#' @export
plot_roc <- sb_plot_roc

#' @rdname sb_plot_performance
#' @export
plot_pr <- sb_plot_pr

#' @rdname sb_plot_calibration
#' @export
plot_calib <- sb_plot_calibration

#' @rdname sb_plot_calibration
#' @export
plot_dca <- sb_plot_dca

#' @rdname sb_plot_performance
#' @export
plot_perf <- sb_plot_performance

#' @rdname sb_plot_performance
#' @export
plot_heatmap <- sb_plot_performance_heatmap

#' @rdname sb_plot_performance
#' @export
plot_confusion <- sb_plot_confusion

#' @rdname sb_plot_calibration
#' @export
plot_prob <- sb_plot_probability_distribution

#' @rdname sb_survival_curves
#' @export
plot_survival <- sb_plot_survival_curves

#' @rdname sb_survival_curves
#' @export
plot_nomogram <- sb_plot_nomogram

#' @rdname sb_model_coefficients
#' @export
plot_coefs <- sb_plot_coefficients

#' @rdname sb_model_shap
#' @export
plot_shap <- sb_plot_shap_beeswarm

#' @rdname sb_model_shap
#' @export
plot_shap_importance <- sb_plot_shap_importance

#' @rdname sb_model_shap
#' @export
plot_shap_dependence <- sb_plot_shap_dependence

#' @rdname sb_model_shap
#' @export
plot_shap_waterfall <- sb_plot_shap_waterfall

#' @rdname sb_plot_performance
#' @export
plot_rank <- sb_plot_model_ranking

#' @rdname sb_plot_performance
#' @export
plot_cohorts <- sb_plot_cohort_metrics

#' @rdname sb_plot_calibration
#' @export
plot_panel <- sb_plot_clinician_panel

#' @rdname sb_native_diagnostics
#' @export
plot_pages <- sb_plot_pages

#' @rdname sb_meta_external_discrimination
#' @export
plot_forest <- sb_plot_meta_forest

#' @rdname sb_selector_diagnostics
#' @export
plot_selector <- sb_plot_selector_diagnostic

#' @rdname sb_native_diagnostics
#' @export
plot_algorithm <- sb_plot_algorithm_diagnostic

# ---- zzzz-specialty-api ----

#' @rdname sb_survival_curves
#' @export
add_surv_curve <- sb_add_survival_curve

#' @rdname sb_model_shap
#' @export
additive_shap <- sb_additive_shap

#' @rdname sb_flexsurv_diagnostics
#' @export
aft_diag <- sb_aft_diagnostics

#' @rdname sb_flexsurv_diagnostics
#' @export
ahaz_diag <- sb_ahaz_diagnostics

#' @rdname sb_model_coefficients
#' @export
path_diag <- sb_algorithm_diagnostics

#' @rdname sb_forest_diagnostics
#' @export
aorsf_diag <- sb_aorsf_diagnostics

#' @rdname sb_knn_neighborhood
#' @export
bag_knn_diag <- sb_bagged_knn_diagnostics

#' @rdname sb_forest_diagnostics
#' @export
bagging_diag <- sb_bagging_diagnostics

#' @rdname sb_gam_effects
#' @export
bart_post <- sb_bart_posterior

#' @rdname sb_gam_effects
#' @export
boost_diag <- sb_boosting_diagnostics

#' @rdname sb_flexsurv_diagnostics
#' @export
bujar_diag <- sb_bujar_diagnostics

#' @rdname sb_forest_diagnostics
#' @export
c50_diag <- sb_c50_diagnostics

#' @rdname sb_pls_diagnostics
#' @export
cancerclass_diag <- sb_cancerclass_diagnostics

#' @rdname sb_forest_diagnostics
#' @export
cforest_diag <- sb_cforest_diagnostics

#' @rdname sb_forest_diagnostics
#' @export
class_tree_diag <- sb_classification_tree_diagnostics

#' @rdname sb_survival_curves
#' @export
cox_forest <- sb_cox_forest

#' @rdname sb_forest_diagnostics
#' @export
ctree_diag <- sb_ctree_structure

#' @rdname sb_pls_diagnostics
#' @export
da_diag <- sb_discriminant_diagnostics

#' @rdname sb_gam_effects
#' @export
earth_diag <- sb_earth_diagnostics

#' @rdname sb_forest_diagnostics
#' @export
evtree_diag <- sb_evtree_structure

#' @rdname sb_pls_diagnostics
#' @export
fda_diag <- sb_fda_diagnostics

#' @rdname sb_forest_diagnostics
#' @export
ferns_diag <- sb_ferns_diagnostics

#' @rdname sb_flexsurv_diagnostics
#' @export
flexsurv_diag <- sb_flexsurv_diagnostics

#' @rdname sb_forest_diagnostics
#' @export
forest_diag <- sb_forest_diagnostics

#' @rdname sb_gam_effects
#' @export
gam_diag <- sb_gam_effects

#' @rdname sb_knn_neighborhood
#' @export
gp_diag <- sb_gaussian_process_diagnostics

#' @rdname sb_analyze_methods
#' @export
gene_overlap <- sb_gene_overlap

#' @rdname sb_analyze_methods
#' @export
gene_rank <- sb_gene_selection_ranking

#' @rdname sb_pls_diagnostics
#' @export
highdim_diag <- sb_high_dimensional_diagnostics

#' @rdname sb_knn_neighborhood
#' @export
knn_diag <- sb_knn_neighborhood

#' @rdname sb_pls_diagnostics
#' @export
ktsp_diag <- sb_ktsp_diagnostics

#' @rdname sb_knn_neighborhood
#' @export
loclda_diag <- sb_loclda_neighborhood

#' @rdname sb_knn_neighborhood
#' @export
lvq_diag <- sb_lvq_diagnostics

#' @rdname sb_knn_neighborhood
#' @export
mclust_diag <- sb_mclust_da_diagnostics

#' @rdname sb_flexsurv_diagnostics
#' @export
mtlr_diag <- sb_mtlr_diagnostics

#' @rdname sb_knn_neighborhood
#' @export
nb_diag <- sb_naive_bayes_diagnostics

#' @rdname sb_native_diagnostics
#' @export
native_diag <- sb_native_diagnostics

#' @rdname sb_gam_effects
#' @export
neural_diag <- sb_neural_diagnostics

#' @rdname sb_survival_curves
#' @export
nomogram <- sb_nomogram

#' @rdname sb_survival_curves
#' @export
nom_values <- sb_nomogram_values

#' @rdname sb_pls_diagnostics
#' @export
pamr_diag <- sb_pamr_diagnostics

#' @rdname sb_flexsurv_diagnostics
#' @export
penaft_diag <- sb_penaft_diagnostics

#' @rdname sb_bootstrap_external
#' @export
plot_boot_calib <- sb_plot_bootstrap_calibration

#' @rdname sb_bootstrap_external
#' @export
plot_ext_diff <- sb_plot_external_difference

#' @rdname sb_bootstrap_external
#' @export
plot_ext_ci <- sb_plot_external_intervals

#' @rdname sb_analyze_methods
#' @export
plot_overlap <- sb_plot_gene_overlap

#' @rdname sb_analyze_methods
#' @export
plot_gene_rank <- sb_plot_gene_selection_ranking

#' @rdname multi_run
#' @export
plot_multi_cm <- sb_plot_multiclass_confusion

#' @rdname sb_plot_performance
#' @export
plot_pred_cm <- sb_plot_prediction_confusion

#' @rdname sb_plot_performance
#' @export
plot_pred_roc <- sb_plot_prediction_roc

#' @rdname sb_analyze_methods
#' @export
plot_freq <- sb_plot_selection_frequency

#' @rdname sb_analyze_methods
#' @export
plot_jaccard <- sb_plot_selection_jaccard

#' @rdname sb_analyze_methods
#' @export
plot_upset <- sb_plot_selection_upset

#' @rdname sb_selector_diagnostics
#' @export
plot_scores <- sb_plot_selector_scores

#' @rdname sb_plot_calibration
#' @export
plot_thresholds <- sb_plot_thresholds

#' @rdname sb_pls_diagnostics
#' @export
pls_diag <- sb_pls_diagnostics

#' @rdname sb_model_coefficients
#' @export
randomglm_diag <- sb_randomglm_diagnostics

#' @rdname sb_forest_diagnostics
#' @export
rlt_diag <- sb_rlt_diagnostics

#' @rdname sb_forest_diagnostics
#' @export
rotation_diag <- sb_rotation_forest_diagnostics

#' @rdname sb_gam_effects
#' @export
rule_diag <- sb_rule_ensemble

#' @rdname sb_selector_diagnostics
#' @export
selector_diag <- sb_selector_diagnostics

#' @rdname sb_pls_diagnostics
#' @export
shrinkage_diag <- sb_shrinkage_discriminant_diagnostics

#' @rdname sb_model_coefficients
#' @export
slope_diag <- sb_slope_diagnostics

#' @rdname sb_pls_diagnostics
#' @export
superpc_diag <- sb_superpc_diagnostics

#' @rdname sb_survival_curves
#' @export
surv_curves <- sb_survival_curves

#' @rdname sb_knn_neighborhood
#' @export
surv_svm_diag <- sb_survival_svm_diagnostics

#' @rdname sb_knn_neighborhood
#' @export
svm_diag <- sb_svm_diagnostics

#' @rdname sb_flexsurv_diagnostics
#' @export
timereg_diag <- sb_timereg_diagnostics

#' @rdname sb_model_shap
#' @export
tree_shap <- sb_tree_shap

#' @rdname sb_forest_diagnostics
#' @export
tree_diag <- sb_tree_structure

#' @rdname sb_forest_diagnostics
#' @export
wsrf_diag <- sb_wsrf_diagnostics

#' @rdname sb_survival_curves
#' @export
surv_roc <- surv_plot_roc

#' @rdname sb_survival_curves
#' @export
surv_risk <- surv_risk_plots

#' @rdname sb_survival_curves
#' @export
risk_ref <- surv_risk_reference

#' @rdname sb_analyze_methods
#' @export
analyze_methods <- sb_analyze_methods

#' @rdname sb_evaluate_external
#' @export
prob_dist <- sb_probability_distribution

#' @rdname sb_evaluate_external
#' @export
threshold_stats <- sb_threshold_metrics

#' @rdname multi_run
#' @export
multi_metrics <- sb_multiclass_metrics

#' @rdname sb_native_diagnostics
#' @export
plot_diag <- sb_plot_algorithm_diagnostic

#' @rdname sb_selector_diagnostics
#' @export
plot_select <- sb_plot_selector_diagnostic

#' @rdname sb_model_shap
#' @export
plot_shap_imp <- sb_plot_shap_importance

#' @rdname sb_model_shap
#' @export
plot_shap_dep <- sb_plot_shap_dependence

#' @rdname sb_model_shap
#' @export
plot_shap_one <- sb_plot_shap_waterfall

#' @rdname sb_plot_performance
#' @export
plot_cm <- sb_plot_confusion

#' @rdname sb_survival_curves
#' @export
plot_surv <- sb_plot_survival_curves
