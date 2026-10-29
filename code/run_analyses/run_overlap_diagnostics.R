# ---------------------------------------------------------------------------------------------------
# Rerun primary IPD analyses to include propensity vectors in output and create overlap diagnostics
# (redundant with main analyses but adding for final sensitivity analysis)
# ---------------------------------------------------------------------------------------------------

here::i_am("code/run_analyses/run_overlap_diagnostics.R")

#library(abxGrowth)
devtools::load_all("../packages/abxGrowth")

library(SuperLearner)

source(here::here("code/run_analyses/SL.wrappers.R"))

# Outcome 1, Missingness -- update for guideline_abx variable
sl.library.with.abx <- list(c("SL.glm", "SL.screen.abx.lt.min_prop"),
                            c("SL.glm", "SL.screen.guideline_abx.glmnet"),
                            c("SL.glm.spline.age.haz", "SL.screen.abx.lt.min_prop"),
                            "SL.step.forward.spline.age.haz.guideline_abx",
                            "SL.glmnet",
                            "SL.ranger",
                            "SL.earth",
                            "SL.xgboost")

# Propensity, outcome 2
sl.library.without.abx <- list(c("SL.glm", "SL.screen.abx.lt.min_prop"),
                               c("SL.glm", "screen.glmnet"),
                               c("SL.glm.spline.age.haz", "SL.screen.abx.lt.min_prop"),
                               "SL.step.forward.spline.age.haz",
                               "SL.glmnet",
                               "SL.ranger",
                               "SL.earth",
                               "SL.xgboost")

# IPD ANALYSIS ----------------------------------------------------------------

# main dataset uses MSD MAL-ED only
data <- readRDS(here::here("data/ipd_data/ipd_data_no_etiology_bin_abx.Rds"))

imp_data <- impute_covariates(data = data,
                              imp_covariates = c("sex",
                                                 "age",
                                                 "ses_quintile",
                                                 "enr_haz",
                                                 "edu_bin",
                                                 "site",
                                                 "I_followup_days",
                                                 "I_followup_days_x_followup_days",
                                                 "dysentery",
                                                 "any_vom",
                                                 "any_fev",
                                                 "dehyd_level",
                                                 "lsstools",
                                                 "duration_pre_enroll",
                                                 "shigella_new" ,
                                                 "rotavirus_new",
                                                 "adenovirus_new" ,
                                                 "etec_new" ,
                                                 "cryptosporidium_new"  ,
                                                 "astrovirus_new" ,
                                                 "norovirus_gii_new" ,
                                                 "tepec_new",
                                                 "campylobacter_new" ,
                                                 "sapovirus_new" ,
                                                 "giardia_new"  ,
                                                 "e_bieneusi_new"  ,
                                                 "eaec_new"))
one_hot_data <- one_hot_encode(imp_data,
                               laz_var_name = "final_haz",
                               covariate_list = c("sex",
                                                  "age",
                                                  "ses_quintile",
                                                  "enr_haz",
                                                  "I_followup_days",
                                                  "I_followup_days_x_followup_days",
                                                  "site",
                                                  "edu_bin",
                                                  "num_hh_lt5",
                                                  "imp_water",
                                                  "imp_sanit"),
                               abx_var_name = "guideline_abx",
                               site_var_name = "site",
                               site_interaction = TRUE,
                               severity_list = c("dysentery",
                                                 "any_vom",
                                                 "lsstools",
                                                 "dehyd_level",
                                                 "duration_pre_enroll",
                                                 "any_fev"),
                               age_var_name = "age")

# The following variables are missing in some studies
# water / sanitation - missing abcd - set to 0
# num children in hh - missing maled - set to 1
# education - missing in abcd - set to 0
# fever - missing in abcd - set to 0

ipd_results <- agaipw(data = one_hot_data$data,
                      laz_var_name = "final_haz",
                      abx_var_name = "guideline_abx", 
                      infection_var_name = "shig_attr",
                      followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
                      site_var_name = one_hot_data$site_var_names,
                      covariate_list = one_hot_data$covariate_list,
                      pathogen_quantity_list = c("shigella_new",
                                                 "rotavirus_new",
                                                 "adenovirus_new",
                                                 "etec_new",
                                                 "cryptosporidium_new",
                                                 "astrovirus_new" ,
                                                 "norovirus_gii_new",
                                                 "tepec_new" ,
                                                 "campylobacter_new",
                                                 "sapovirus_new",
                                                 "e_bieneusi_new",
                                                 "giardia_new",
                                                 "eaec_new"),
                      severity_list = one_hot_data$severity_list,
                      no_etiology_var_name = "no_etiology",
                      first_id_var_name = "first_id",
                      outcome_type = "gaussian",
                      sl.library.outcome = sl.library.with.abx,
                      sl.library.outcome.2 = sl.library.without.abx,
                      sl.library.treatment = sl.library.without.abx,
                      sl.library.infection = sl.library.without.abx,
                      sl.library.missingness = sl.library.with.abx,
                      v_folds = 5, 
                      return_models = TRUE,
                      msm = TRUE,
                      msm_var_name = "age",
                      msm_formula = "age")

# get rid of models but save vectors bc object huge
ipd_results$aipw_est$aipw_models$outcome_model_1a <- NULL
ipd_results$aipw_est$aipw_models$outcome_model_1b <- NULL
ipd_results$aipw_est$aipw_models$outcome_model_2a <- NULL
ipd_results$aipw_est$aipw_models$outcome_model_2b_list <- NULL
ipd_results$aipw_est$aipw_models$prop_model_1a_list <- NULL
ipd_results$aipw_est$aipw_models$prop_model_1b_list <- NULL
ipd_results$aipw_est$aipw_models$prop_model_2a_1 <- NULL
ipd_results$aipw_est$aipw_models$prop_model_2a_2 <- NULL
ipd_results$aipw_est$aipw_models$prop_model_3a <- NULL
ipd_results$aipw_est$aipw_models$prop_model_3b <- NULL
ipd_results$aipw_est$aipw_models$msm_model_list <- NULL
saveRDS(ipd_results, here::here("results/no_etiology/ipd_no_etiology_guideline_abx_with_prop_vec.Rds"))

# ------------------------------------------------------------------------------

# MSD only
ipd_data <- readRDS(here::here("data/ipd_data/ipd_data_case_control_msd_tac_or_culture_shig_bin_abx.Rds"))

# Impute missing values, cases and controls separate
case_ipd_data <- ipd_data[ipd_data$case == 1,]
control_ipd_data <- ipd_data[ipd_data$case == 0,]

imp_data_case <- impute_covariates(data = case_ipd_data,
                                   imp_covariates = c("sex",
                                                      "age",
                                                      "ses_quintile",
                                                      "enr_haz",
                                                      "edu_bin",
                                                      "site",
                                                      "I_followup_days",
                                                      "I_followup_days_x_followup_days",
                                                      "dysentery",
                                                      "any_vom",
                                                      "any_fev",
                                                      "dehyd_level",
                                                      "lsstools",
                                                      "duration_pre_enroll",
                                                      "shigella_new" ,
                                                      "rotavirus_new",
                                                      "adenovirus_new" ,
                                                      "etec_new" ,
                                                      "cryptosporidium_new"  ,
                                                      "astrovirus_new" ,
                                                      "norovirus_new" ,
                                                      "tepec_new",
                                                      "campylobacter_new" ,
                                                      "sapovirus_new" ,
                                                      "giardia_new"  ,
                                                      "e_bieneusi_new"  ,
                                                      "eaec_new"))

imp_data_control <- impute_covariates(data = control_ipd_data,
                                      imp_covariates = c("sex",
                                                         "age",
                                                         "ses_quintile",
                                                         "enr_haz",
                                                         "edu_bin",
                                                         "site",
                                                         "I_followup_days",
                                                         "I_followup_days_x_followup_days"))


imp_data_case <- imp_data_case[,colnames(imp_data_case) %in% colnames(ipd_data), drop = FALSE]
imp_data_control <- imp_data_control[,colnames(imp_data_control) %in% colnames(ipd_data), drop = FALSE]
imp_data <- rbind(imp_data_case, imp_data_control)

one_hot_data <- one_hot_encode(imp_data,
                               laz_var_name = "final_haz",
                               covariate_list = c("sex",
                                                  "age",
                                                  "ses_quintile",
                                                  "enr_haz",
                                                  "I_followup_days",
                                                  "I_followup_days_x_followup_days",
                                                  "site",
                                                  "edu_bin",
                                                  "num_hh_lt5",
                                                  "imp_water",
                                                  "imp_sanit"),
                               abx_var_name = "guideline_abx",
                               site_var_name = "site",
                               site_interaction = TRUE,
                               severity_list = c("dysentery",
                                                 "any_vom",
                                                 "lsstools",
                                                 "dehyd_level",
                                                 "duration_pre_enroll",
                                                 "any_fev"),
                               age_var_name = "age")

ipd_results_cc <- agaipw(data = one_hot_data$data,
                         laz_var_name = "final_haz",
                         abx_var_name = "guideline_abx", 
                         case_var_name = "case",
                         followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
                         site_var_name = one_hot_data$site_var_names,
                         covariate_list = one_hot_data$covariate_list,
                         pathogen_quantity_list = c("shigella_new",
                                                    "rotavirus_new",
                                                    "adenovirus_new",
                                                    "etec_new",
                                                    "cryptosporidium_new",
                                                    "astrovirus_new" ,
                                                    "norovirus_new",
                                                    "tepec_new" ,
                                                    "campylobacter_new",
                                                    "sapovirus_new",
                                                    "e_bieneusi_new",
                                                    "giardia_new",
                                                    "eaec_new"),
                         severity_list = one_hot_data$severity_list,
                         first_id_var_name = "first_id",
                         outcome_type = "gaussian",
                         sl.library.outcome.case = sl.library.with.abx,
                         sl.library.outcome.control = sl.library.without.abx,
                         sl.library.treatment = sl.library.without.abx,
                         sl.library.infection = sl.library.without.abx,
                         sl.library.missingness.case = sl.library.with.abx,
                         sl.library.missingness.control = sl.library.without.abx,
                         v_folds = 5, 
                         return_models = TRUE,
                         msm = TRUE,
                         msm_var_name = "age",
                         msm_formula = "age",
                         case_control = TRUE)

saveRDS(ipd_results_cc, "results/case_control/ipd_results_msd_case_control_tac_or_culture_guideline_abx_with_prop_vec.Rds")

