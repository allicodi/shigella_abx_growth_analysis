# -------------------------------------------------------------------------------
# Script to run IPD no etiology & case control analyses stratified by age
# -------------------------------------------------------------------------------

here::i_am("code/run_analyses/run_age_stratified.R")

#library(abxGrowth)
devtools::load_all("../packages/abxGrowth")

library(SuperLearner)

source(here::here("code/run_analyses/SL.wrappers.R"))

# main dataset uses MSD MAL-ED only
no_etiology <- readRDS(here::here("data/ipd_data/ipd_data_no_etiology.Rds"))
case_control <- readRDS(here::here("data/ipd_data/ipd_data_case_control_msd_tac_or_culture_shig.Rds"))

# Outcome 1, Missingness
sl.library.with.abx <- list(c("SL.glm", "SL.screen.abx.lt.min_prop"),
                            c("SL.glm", "SL.screen.abx.glmnet"),
                            c("SL.glm.spline.age.haz", "SL.screen.abx.lt.min_prop"),
                            "SL.step.forward.spline.age.haz.abx",
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

# ------------------------------------------------------------------------------
# No etiology
# ------------------------------------------------------------------------------

imp_data <- impute_covariates(data = no_etiology,
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
                               abx_var_name = "all_abx",
                               site_var_name = "site",
                               site_interaction = TRUE,
                               severity_list = c("dysentery",
                                                 "any_vom",
                                                 "lsstools",
                                                 "dehyd_level",
                                                 "duration_pre_enroll",
                                                 "any_fev"),
                               age_var_name = "age")

one_hot_data_0_11 <- one_hot_data$data[one_hot_data$data$age >= 0 & one_hot_data$data$age < 12, ]
one_hot_data_12_23 <- one_hot_data$data[one_hot_data$data$age >= 12 & one_hot_data$data$age < 24, ]
one_hot_data_24_59 <- one_hot_data$data[one_hot_data$data$age >= 24 & one_hot_data$data$age < 60, ]

# The following variables are missing in some studies
# water / sanitation - missing abcd - set to 0
# num children in hh - missing maled - set to 1
# education - missing in abcd - set to 0
# fever - missing in abcd - set to 0

ipd_results_0_11 <- agaipw(data = one_hot_data_0_11,
                      laz_var_name = "final_haz",
                      abx_var_name = "all_abx", 
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

ipd_results_12_23 <- agaipw(data = one_hot_data_12_23,
                           laz_var_name = "final_haz",
                           abx_var_name = "all_abx", 
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

ipd_results_24_59 <- agaipw(data = one_hot_data_24_59,
                           laz_var_name = "final_haz",
                           abx_var_name = "all_abx", 
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

ipd_results_0_11$aipw_est$aipw_models <- NULL
ipd_results_12_23$aipw_est$aipw_models <- NULL
ipd_results_24_59$aipw_est$aipw_models <- NULL

saveRDS(ipd_results_0_11, here::here("results/no_etiology/ipd_no_etiology_0_11_no_models.Rds"))
saveRDS(ipd_results_12_23, here::here("results/no_etiology/ipd_no_etiology_12_23_no_models.Rds"))
saveRDS(ipd_results_24_59, here::here("results/no_etiology/ipd_no_etiology_24_59_no_models.Rds"))

# ------------------------------------------------------------------------------
# Case-control
# ------------------------------------------------------------------------------

# Impute missing values, cases and controls separate
case_ipd_data <- case_control[case_control$case == 1,]
control_ipd_data <- case_control[case_control$case == 0,]

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


imp_data_case <- imp_data_case[,colnames(imp_data_case) %in% colnames(case_control), drop = FALSE]
imp_data_control <- imp_data_control[,colnames(imp_data_control) %in% colnames(case_control), drop = FALSE]
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
                               abx_var_name = "all_abx",
                               site_var_name = "site",
                               site_interaction = TRUE,
                               severity_list = c("dysentery",
                                                 "any_vom",
                                                 "lsstools",
                                                 "dehyd_level",
                                                 "duration_pre_enroll",
                                                 "any_fev"),
                               age_var_name = "age")

one_hot_data_0_11 <- one_hot_data$data[one_hot_data$data$age >= 0 & one_hot_data$data$age < 12,]
one_hot_data_12_23 <- one_hot_data$data[one_hot_data$data$age >= 12 & one_hot_data$data$age < 23,]
one_hot_data_24_59 <- one_hot_data$data[one_hot_data$data$age >= 24 & one_hot_data$data$age < 60,]

ipd_results_cc_0_11 <- agaipw(data = one_hot_data_0_11,
                         laz_var_name = "final_haz",
                         abx_var_name = "all_abx", 
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

ipd_results_cc_12_23 <- agaipw(data = one_hot_data_12_23,
                              laz_var_name = "final_haz",
                              abx_var_name = "all_abx", 
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

ipd_results_cc_24_59 <- agaipw(data = one_hot_data_24_59,
                              laz_var_name = "final_haz",
                              abx_var_name = "all_abx", 
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

ipd_results_cc_0_11$aipw_est$aipw_models <- NULL
ipd_results_cc_12_23$aipw_est$aipw_models <- NULL
ipd_results_cc_24_59$aipw_est$aipw_models <- NULL

saveRDS(ipd_results_cc_0_11, here::here("results/case_control/ipd_cc_0_11_no_models.Rds"))
saveRDS(ipd_results_cc_12_23, here::here("results/case_control/ipd_cc_12_23_no_models.Rds"))
saveRDS(ipd_results_cc_24_59, here::here("results/case_control/ipd_cc_24_59_no_models.Rds"))




