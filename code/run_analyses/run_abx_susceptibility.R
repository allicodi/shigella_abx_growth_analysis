# ---------------------------------------------------------------------------
# Susceptibility analyses
# ---------------------------------------------------------------------------

here::i_am("code/run_analyses/run_abx_susceptibility.R")

devtools::load_all("../packages/abxGrowth/")

library(SuperLearner)
library(tidyr)

source(here::here("code/run_analyses/SL.wrappers.R"))

efgh_data <- readRDS(here::here("data/efgh_data/efgh_data.Rds"))

# Make var that counts NA abx test as its own thing
efgh_data$ast_given_abx_numeric <- ifelse(is.na(efgh_data$ast_given_abx), 0, 
                                          ifelse(efgh_data$ast_given_abx == "S", 1, 2))

# Outcome 1, Missingness
sl.library.with.abx <- list(c("SL.glm", "SL.screen.abx.lt.min_prop"),
                            c("SL.glm", "SL.screen.abx.glmnet"),
                            c("SL.glm.spline.age.haz", "SL.screen.abx.lt.min_prop"),
                            "SL.step.forward.spline.age.haz.abx",
                            "SL.glmnet",
                            #"SL.biglasso",
                            "SL.ranger",
                            "SL.earth",
                            "SL.xgboost")

# Propensity, outcome 2
sl.library.without.abx <- list(
  c("SL.glm", "SL.screen.abx.lt.min_prop"),
  c("SL.glm", "screen.glmnet"),
  c("SL.glm.spline.age.haz", "SL.screen.abx.lt.min_prop"),
  "SL.step.forward.spline.age.haz",
  "SL.glmnet",
  #"SL.biglasso",
  "SL.ranger",
  "SL.earth",
  "SL.xgboost")

efgh_data_imp <- impute_covariates(data = efgh_data,
                                   site_var_name = "enroll_site",
                                   imp_covariates = c("sex",
                                                      "enr_age_months",
                                                      "enr_haz",
                                                      "final_quintile_site",
                                                      "imp_water",
                                                      "imp_toi",
                                                      "enroll_ai_num_child",
                                                      "moth_ed_bin",
                                                      "I_mo3_days",
                                                      "I_mo3_days_x_mo3_days",
                                                      "enroll_diar_blood",
                                                      "enroll_diar_vom_days",
                                                      "enroll_diar_vom_num",
                                                      "enroll_diar_fever_days",
                                                      "enroll_diar_fever",
                                                      "enroll_diar_loose_num",
                                                      "enroll_cond_dehyd",
                                                      "duration_pre_enroll",
                                                      "shigella_new",               
                                                      "rotavirus_new",              
                                                      "adenovirus_new",             
                                                      "ETEC_new",                   
                                                      "cryptosporidium_new",        
                                                      "astrovirus_new",             
                                                      "norovirus_gii_new",              
                                                      "c_jejuni_new",               
                                                      "tEPEC_new",                  
                                                      "sapovirus_new",              
                                                      "e_bieneusi_new",             
                                                      "giardia_new",                
                                                      "EAEC_new"))

efgh_data_imp$age <- efgh_data_imp$enr_age_months
efgh_data_imp <- efgh_data_imp[which(!is.na(efgh_data_imp$tac_shigella_attributable)),]

one_hot_efgh <- one_hot_encode(data = efgh_data_imp,
                               laz_var_name = "mo3_haz",
                               covariate_list = c("sex",
                                                  "age",
                                                  "enr_haz",
                                                  "final_quintile_site",
                                                  "imp_water",
                                                  "imp_toi",
                                                  "enroll_ai_num_child",
                                                  "moth_ed_bin",
                                                  "I_mo3_days",
                                                  "I_mo3_days_x_mo3_days"),
                               abx_var_name = "all_abx",
                               site_var_name = "enroll_site",
                               site_interaction = TRUE,
                               severity_list = c("enroll_diar_blood",
                                                 "enroll_diar_vom_days",
                                                 "enroll_diar_vom_num",
                                                 "enroll_diar_fever_days",
                                                 "enroll_diar_fever",
                                                 "enroll_diar_loose_num",
                                                 "enroll_cond_dehyd",
                                                 "duration_pre_enroll"),
                               age_var_name = "age")

## Main results: EFGH ------------------------------------------------------------
# ITERATION 1: Normal. Save results + pull no etiology results, guideline recommended drug

efgh_results <- agaipw(data = one_hot_efgh$data,
                       laz_var_name = "mo3_haz",
                       abx_var_name = "all_abx", 
                       infection_var_name = "positive_tac_or_culture",
                       site_var_name = one_hot_efgh$site_var_names,
                       followup_var_names = c("I_mo3_days",
                                              "I_mo3_days_x_mo3_days"),
                       covariate_list = one_hot_efgh$covariate_list,
                       pathogen_quantity_list = c("shigella_new",               
                                                  "rotavirus_new",              
                                                  "adenovirus_new",             
                                                  "ETEC_new",                   
                                                  "cryptosporidium_new",        
                                                  "astrovirus_new",             
                                                  "norovirus_gii_new",              
                                                  "c_jejuni_new",               
                                                  "tEPEC_new",                  
                                                  "sapovirus_new",              
                                                  "e_bieneusi_new",             
                                                  "giardia_new",                
                                                  "EAEC_new"),
                       severity_list = one_hot_efgh$severity_list,
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

saveRDS(efgh_results, here::here("results/no_etiology/susceptibility/no_etiology.Rds"))

# ITERATION 2: Subset to infection definition to children with resistant Shigella

# ITERATION 2: keep shigella attr same, no_etiology = I(resistant to at least one WHO approved thing == 1) -- abx var = (0,1,2) where (I/R, S, NA to the drug you got)
#               then take the 'no etiology' aka resistant shigella output, no etiology + abx = resistant
#               compare to no etiology, WHO approved (iteration 1)

# abx resistance variable (ast_given_abx_numeric)
# 0 = NA, 1 = S, 2 = I/R

# Rename all_abx so i don't have to make a different version of the wrappers...
one_hot_efgh$data$all_abx_original <- one_hot_efgh$data$all_abx
one_hot_efgh$data$all_abx <- one_hot_efgh$data$ast_given_abx_numeric

efgh_results_IR <- agaipw(data = one_hot_efgh$data,
                          laz_var_name = "mo3_haz",
                          abx_var_name = "all_abx", 
                          infection_var_name = "positive_tac_or_culture",
                          site_var_name = one_hot_efgh$site_var_names,
                          followup_var_names = c("I_mo3_days",
                                                 "I_mo3_days_x_mo3_days"),
                          covariate_list = one_hot_efgh$covariate_list,
                          pathogen_quantity_list = c("shigella_new",               
                                                     "rotavirus_new",              
                                                     "adenovirus_new",             
                                                     "ETEC_new",                   
                                                     "cryptosporidium_new",        
                                                     "astrovirus_new",             
                                                     "norovirus_gii_new",              
                                                     "c_jejuni_new",               
                                                     "tEPEC_new",                  
                                                     "sapovirus_new",              
                                                     "e_bieneusi_new",             
                                                     "giardia_new",                
                                                     "EAEC_new"),
                          severity_list = one_hot_efgh$severity_list,
                          no_etiology_var_name = "resistant_WHO_approve",
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
                          msm_formula = "age",
                          parsimonious_propensity = FALSE)

saveRDS(efgh_results_IR, here::here("results/no_etiology/susceptibility/resistant_results.Rds"))

# ITERATION 3: Subset to infection definition to children with susceptible

# keep shigella attr same, no_etiology = I(susceptible to at least one WHO approved thing == 1) -- abx var = (0,1,2) where (I/R, S, NA to the drug you got)

efgh_results_S <- agaipw(data = one_hot_efgh$data,
                         laz_var_name = "mo3_haz",
                         abx_var_name = "all_abx", 
                         infection_var_name = "positive_tac_or_culture",
                         site_var_name = one_hot_efgh$site_var_names,
                         followup_var_names = c("I_mo3_days",
                                                "I_mo3_days_x_mo3_days"),
                         covariate_list = one_hot_efgh$covariate_list,
                         pathogen_quantity_list = c("shigella_new",               
                                                    "rotavirus_new",              
                                                    "adenovirus_new",             
                                                    "ETEC_new",                   
                                                    "cryptosporidium_new",        
                                                    "astrovirus_new",             
                                                    "norovirus_gii_new",              
                                                    "c_jejuni_new",               
                                                    "tEPEC_new",                  
                                                    "sapovirus_new",              
                                                    "e_bieneusi_new",             
                                                    "giardia_new",                
                                                    "EAEC_new"),
                         severity_list = one_hot_efgh$severity_list,
                         no_etiology_var_name = "susceptible_WHO_approve",
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
                         msm_formula = "age",
                         parsimonious_propensity = FALSE)

saveRDS(efgh_results_S, here::here("results/no_etiology/susceptibility/susceptible_results.Rds"))

