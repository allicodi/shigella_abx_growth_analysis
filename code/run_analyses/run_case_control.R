# ---------------------------------------------------------------------------
# Script to run no-etiology case-only analysis for meta-analysis
# ---------------------------------------------------------------------------

here::i_am("code/run_analyses/run_case_control.R")

#library(abxGrowth)
devtools::load_all("../packages/abxGrowth")

library(SuperLearner)

source(here::here("code/run_analyses/SL.wrappers.R"))

# Load data ----------------------------------------------

## PRIMARY ANALYSIS: Case definition = TAC or culture attributable Shigella
gems_data_tac_culture <- readRDS(here::here("data/gems_data/gems_case_control_tac_or_culture_shig.Rds"))
vida_data_tac_culture <- readRDS(here::here("data/vida_data/vida_case_control_tac_or_culture_shig.Rds"))
maled_data_tac_culture <- readRDS(here::here("data/maled_data/maled_case_control_tac_or_culture_shig.Rds"))

# MAL-ED MSD only
maled_msd_tac_culture <- maled_data_tac_culture
maled_data_msd_caseids <- unique(maled_msd_tac_culture$case_id[which(maled_msd_tac_culture$MSD == 1)])
maled_msd_tac_culture <- maled_msd_tac_culture[which(maled_msd_tac_culture$case_id %in% maled_data_msd_caseids),]

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

# GEMS ANALYSIS ------------------------------------------------------------

## Primary: TAC or culture definition
# Impute missing values, cases and controls separate
gems_case_data <- gems_data_tac_culture[gems_data_tac_culture$case == 1,]
gems_control_data <- gems_data_tac_culture[gems_data_tac_culture$case == 0,]

# Impute missing covariates
gems_case_data_imp <- impute_covariates(gems_case_data,
                                        imp_covariates = c("sex",
                                                           "age",
                                                           "ses_quintile",
                                                           "safe_water",
                                                           "safe_sanit",
                                                           "enr_haz",
                                                           "prim_caregiver_edu_bin",
                                                           "num_hh_lt5",
                                                           "shigella_new",
                                                           "rotavirus_new",
                                                           "adenovirus_new",
                                                           "st_etec_new",
                                                           "lt_etec_new",
                                                           "crypto_new",
                                                           "astro_new" ,
                                                           "noro_new",
                                                           "tepec_new" ,
                                                           "campy_new",
                                                           "sapo_new",
                                                           "e_bieneusi_new",
                                                           "giardia_new",
                                                           "EAEC_new", 
                                                           "dysentery",
                                                           "vomit",
                                                           "fever",
                                                           "lsstools",
                                                           "who_dehyd",
                                                           "duration_pre_enroll"))

gems_control_data_imp <- impute_covariates(gems_control_data,
                                           imp_covariates = c("sex",
                                                              "age",
                                                              "ses_quintile",
                                                              "safe_water",
                                                              "safe_sanit",
                                                              "enr_haz",
                                                              "prim_caregiver_edu_bin",
                                                              "num_hh_lt5"))

gems_case_data_imp <- gems_case_data_imp[,colnames(gems_case_data_imp) %in% colnames(gems_data_tac_culture), drop = FALSE]
gems_control_data_imp <- gems_control_data_imp[,colnames(gems_control_data_imp) %in% colnames(gems_data_tac_culture), drop = FALSE]
gems_data_imp <- rbind(gems_case_data_imp, gems_control_data_imp)

one_hot_gems <- one_hot_encode(gems_data_imp,
                               laz_var_name = "hazd60",
                               covariate_list = c("sex",
                                                  "age",
                                                  "ses_quintile",
                                                  "safe_water",
                                                  "safe_sanit",
                                                  "enr_haz",
                                                  "prim_caregiver_edu_bin",
                                                  "num_hh_lt5",
                                                  "I_followup_days",
                                                  "I_followup_days_x_followup_days",
                                                  "site"),
                               abx_var_name = "all_abx",
                               site_var_name = "site",
                               site_interaction = TRUE,
                               severity_list = c("dysentery",
                                                 "vomit",
                                                 "fever",
                                                 "lsstools",
                                                 "who_dehyd",
                                                 "duration_pre_enroll"),
                               age_var_name = "age")

gems_results <- agaipw(data = one_hot_gems$data,
                       laz_var_name = "hazd60",
                       abx_var_name = "all_abx", 
                       case_var_name = "case",
                       followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
                       site_var_name = one_hot_gems$site_var_names,
                       covariate_list = one_hot_gems$covariate_list,
                       pathogen_quantity_list = c("shigella_new",
                                                  "rotavirus_new",
                                                  "adenovirus_new",
                                                  "st_etec_new",
                                                  "lt_etec_new",
                                                  "crypto_new",
                                                  "astro_new" ,
                                                  "noro_new",
                                                  "tepec_new" ,
                                                  "campy_new",
                                                  "sapo_new",
                                                  "e_bieneusi_new",
                                                  "giardia_new",
                                                  "EAEC_new"),
                       severity_list = one_hot_gems$severity_list,
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

saveRDS(gems_results, here::here("results/case_control/gems_tac_or_culture.Rds"))
gems_results$aipw_est$aipw_models <- NULL
saveRDS(gems_results, here::here("results/case_control/gems_tac_or_culture_no_models.Rds"))

# VIDA ANALYSIS ------------------------------------------------------------

# Impute missing values, cases and controls separate
vida_case_data <- vida_data_tac_culture[vida_data_tac_culture$case == 1,]
vida_control_data <- vida_data_tac_culture[vida_data_tac_culture$case == 0,]

vida_case_data_imp <- impute_covariates(data = vida_case_data,
                                        site_var_name = "site",
                                        imp_by_site = TRUE,
                                        imp_covariates = c("sex",
                                                           "agemchild",
                                                           "ses_quintile",
                                                           "safe_water",
                                                           "safe_sanit",
                                                           "enr_haz",
                                                           "site",
                                                           "education_bin",
                                                           "num_hh_lt5",
                                                           "shigella_new",                   
                                                           "rotavirus_new",                  
                                                           "st_etec_new",                    
                                                           "crypto_new",                     
                                                           "adeno_new",                      
                                                           "astro_new",                      
                                                           "noro_new",                       
                                                           "tepec_new",                      
                                                           "campy_new",                      
                                                           "sapo_new",                       
                                                           "giardia_new",                    
                                                           "eaec_new",
                                                           "lsstools",
                                                           "vom_days",
                                                           "vom_freq",
                                                           "fever",
                                                           "dehydr",
                                                           "dysentery",
                                                           "duration_pre_enroll"))

vida_control_data_imp <- impute_covariates(data = vida_control_data,
                                           site_var_name = "site",
                                           imp_by_site = TRUE,
                                           imp_covariates = c("sex",
                                                              "agemchild",
                                                              "ses_quintile",
                                                              "safe_water",
                                                              "safe_sanit",
                                                              "enr_haz",
                                                              "site",
                                                              "education_bin",
                                                              "num_hh_lt5"))

vida_case_data_imp <- vida_case_data_imp[,colnames(vida_case_data_imp) %in% colnames(vida_data_tac_culture), drop = FALSE]
vida_control_data_imp <- vida_control_data_imp[,colnames(vida_control_data_imp) %in% colnames(vida_data_tac_culture), drop = FALSE]
vida_data_imp <- rbind(vida_case_data_imp, vida_control_data_imp)

vida_data_imp$age <- vida_data_imp$agemchild

one_hot_vida <- one_hot_encode(data = vida_data_imp,
                               laz_var_name = "hazd60",
                               covariate_list = c("sex",
                                                  "age",
                                                  "ses_quintile",
                                                  "safe_water",
                                                  "safe_sanit",
                                                  "enr_haz",
                                                  "site",
                                                  "education_bin",
                                                  "num_hh_lt5",
                                                  "I_followup_days",
                                                  "I_followup_days_x_followup_days"),
                               abx_var_name = "all_abx",
                               site_var_name = "site",
                               site_interaction = TRUE,
                               severity_list = c("lsstools",
                                                 "vom_days",
                                                 "vom_freq",
                                                 "fever",
                                                 "dehydr",
                                                 "dysentery",
                                                 "duration_pre_enroll"),
                               age_var_name = "age")

vida_results <- agaipw(data = one_hot_vida$data,
                       laz_var_name = "hazd60",
                       abx_var_name = "all_abx", 
                       case_var_name = "case",
                       followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
                       site_var_name = one_hot_vida$site_var_names,
                       covariate_list = one_hot_vida$covariate_list,
                       pathogen_quantity_list = c("shigella_new",                   
                                                  "rotavirus_new",                  
                                                  "st_etec_new",                    
                                                  "crypto_new",                     
                                                  "adeno_new",                      
                                                  "astro_new",                      
                                                  "noro_new",                       
                                                  "tepec_new",                      
                                                  "campy_new",                      
                                                  "sapo_new",                       
                                                  "giardia_new",                    
                                                  "eaec_new"),
                       severity_list = one_hot_vida$severity_list,
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

saveRDS(vida_results, here::here("results/case_control/vida_tac_or_culture.Rds"))
vida_results$aipw_est$aipw_models <- NULL
saveRDS(vida_results, here::here("results/case_control/vida_tac_or_culture_no_models.Rds"))

# MAL-ED ANALYSIS ----------------------------------------------------------

# Impute missing values, cases and controls separate
maled_case_data <- maled_data_tac_culture[maled_data_tac_culture$case == 1,]
maled_control_data <- maled_data_tac_culture[maled_data_tac_culture$case == 0,]

maled_case_data_imp <- impute_covariates(data = maled_case_data,
                                         site_var_name = "site",
                                         imp_by_site = TRUE,
                                         imp_covariates = c("site",
                                                            "sex",
                                                            "agemonths",
                                                            "baseline_haz",
                                                            "mated_bin", 
                                                            "drinkimp",
                                                            "sanitimp",
                                                            "wami_quintile",
                                                            "adenovirus_40_41_new",                  
                                                            "aeromonas_new",                         
                                                            "astrovirus_new",                        
                                                            "campylobacter_pan_new",                 
                                                            "cryptosporidium_new",                   
                                                            "cyclospora_new",                        
                                                            "e_histolytica_new",                     
                                                            "isospora_new",                          
                                                            "norovirus_new",                         
                                                            "rotavirus_new",                         
                                                            "salmonella_new",                        
                                                            "sapovirus_new",  
                                                            "shigella_new",
                                                            "st_etec_new",                           
                                                            "tEPEC_new",                             
                                                            "v_cholerae_new",                        
                                                            "ETEC_new",                              
                                                            "e_bieneusi_new",                        
                                                            "eaec_new",
                                                            "dysentery",
                                                            "fever",
                                                            "fever_days",
                                                            "dehyd",
                                                            "lsstools",
                                                            "daysvomit"))

maled_control_data_imp <- impute_covariates(data = maled_control_data,
                                            site_var_name = "site",
                                            imp_by_site = TRUE,
                                            imp_covariates = c("site",
                                                               "sex",
                                                               "agemonths",
                                                               "baseline_haz",
                                                               "wami_quintile",
                                                               "mated_bin", 
                                                               "drinkimp",
                                                               "sanitimp"))

maled_case_data_imp <- maled_case_data_imp[,colnames(maled_case_data_imp) %in% colnames(maled_data_tac_culture), drop = FALSE]
maled_control_data_imp <- maled_control_data_imp[,colnames(maled_control_data_imp) %in% colnames(maled_data_tac_culture), drop = FALSE]
maled_data_imp <- rbind(maled_case_data_imp, maled_control_data_imp)

# Add "age" and "enr_haz" to match wrappers (go back and rename in raw data later)
maled_data_imp$age <- maled_data_imp$agemonths
maled_data_imp$enr_haz <- maled_data_imp$baseline_haz

one_hot_maled <- one_hot_encode(data = maled_data_imp,
                                laz_var_name = "month3_haz",
                                covariate_list = c("site",
                                                   "sex",
                                                   "age",
                                                   "enr_haz",
                                                   "mated_bin", 
                                                   "drinkimp",
                                                   "sanitimp",
                                                   "wami_quintile",
                                                   "I_followup_days",
                                                   "I_followup_days_x_followup_days"),
                                abx_var_name = "all_abx",
                                site_var_name = "site",
                                site_interaction = TRUE,
                                severity_list = c("dysentery",
                                                  "fever",
                                                  "fever_days",
                                                  "dehyd",
                                                  "lsstools",
                                                  "daysvomit",
                                                  "duration_pre_abx"),
                                age_var_name = "age")

maled_results <- agaipw(data = one_hot_maled$data,
                        laz_var_name = "month3_haz",
                        abx_var_name = "all_abx", 
                        case_var_name = "case",
                        site_var_name = one_hot_maled$site_var_names,
                        followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
                        covariate_list = one_hot_maled$covariate_list,
                        pathogen_quantity_list = c("shigella_new",            
                                                   "adenovirus_40_41_new",                  
                                                   "aeromonas_new",                         
                                                   "astrovirus_new",                        
                                                   "campylobacter_pan_new",                 
                                                   "cryptosporidium_new",                   
                                                   "cyclospora_new",                        
                                                   "e_histolytica_new",                     
                                                   "isospora_new",                          
                                                   "norovirus_new",                         
                                                   "rotavirus_new",                         
                                                   "salmonella_new",                        
                                                   "sapovirus_new",                         
                                                   "st_etec_new",                           
                                                   "tEPEC_new",                             
                                                   "v_cholerae_new",                        
                                                   "ETEC_new",                              
                                                   "e_bieneusi_new",                        
                                                   "eaec_new"),
                        severity_list = one_hot_maled$severity_list,
                        no_etiology_var_name = "no_etiology",
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
                        case_control = TRUE,
                        first_id_var_name = "first_id")

saveRDS(maled_results, here::here("results/case_control/maled_tac_or_culture.Rds"))
maled_results$aipw_est$aipw_models <- NULL
saveRDS(maled_results, here::here("results/case_control/maled_tac_or_culture_no_models.Rds"))

#### MSD ONLY ####

# Impute missing values, cases and controls separate
maled_case_data <- maled_msd_tac_culture[maled_msd_tac_culture$case == 1,]
maled_control_data <- maled_msd_tac_culture[maled_msd_tac_culture$case == 0,]

maled_case_data_imp <- impute_covariates(data = maled_case_data,
                                         site_var_name = "site",
                                         imp_by_site = TRUE,
                                         imp_covariates = c("site",
                                                            "sex",
                                                            "agemonths",
                                                            "baseline_haz",
                                                            "mated_bin", 
                                                            "drinkimp",
                                                            "sanitimp",
                                                            "wami_quintile",
                                                            "adenovirus_40_41_new",                  
                                                            "aeromonas_new",                         
                                                            "astrovirus_new",                        
                                                            "campylobacter_pan_new",                 
                                                            "cryptosporidium_new",                   
                                                            "cyclospora_new",                        
                                                            "e_histolytica_new",                     
                                                            "isospora_new",                          
                                                            "norovirus_new",                         
                                                            "rotavirus_new",                         
                                                            "salmonella_new",                        
                                                            "sapovirus_new",  
                                                            "shigella_new",
                                                            "st_etec_new",                           
                                                            "tEPEC_new",                             
                                                            "v_cholerae_new",                        
                                                            "ETEC_new",                              
                                                            "e_bieneusi_new",                        
                                                            "eaec_new",
                                                            "dysentery",
                                                            "fever",
                                                            "fever_days",
                                                            "dehyd",
                                                            "lsstools",
                                                            "daysvomit"))

maled_control_data_imp <- impute_covariates(data = maled_control_data,
                                            site_var_name = "site",
                                            imp_by_site = TRUE,
                                            imp_covariates = c("site",
                                                               "sex",
                                                               "agemonths",
                                                               "baseline_haz",
                                                               "wami_quintile",
                                                               "mated_bin", 
                                                               "drinkimp",
                                                               "sanitimp"))

maled_case_data_imp <- maled_case_data_imp[,colnames(maled_case_data_imp) %in% colnames(maled_msd_tac_culture), drop = FALSE]
maled_control_data_imp <- maled_control_data_imp[,colnames(maled_control_data_imp) %in% colnames(maled_msd_tac_culture), drop = FALSE]
maled_data_imp <- rbind(maled_case_data_imp, maled_control_data_imp)

# Add "age" and "enr_haz" to match wrappers (go back and rename in raw data later)
maled_data_imp$age <- maled_data_imp$agemonths
maled_data_imp$enr_haz <- maled_data_imp$baseline_haz

one_hot_maled <- one_hot_encode(data = maled_data_imp,
                                laz_var_name = "month3_haz",
                                covariate_list = c("site",
                                                   "sex",
                                                   "age",
                                                   "enr_haz",
                                                   "mated_bin", 
                                                   "drinkimp",
                                                   "sanitimp",
                                                   "wami_quintile",
                                                   "I_followup_days",
                                                   "I_followup_days_x_followup_days"),
                                abx_var_name = "all_abx",
                                site_var_name = "site",
                                site_interaction = TRUE,
                                severity_list = c("dysentery",
                                                  "fever",
                                                  "fever_days",
                                                  "dehyd",
                                                  "lsstools",
                                                  "daysvomit",
                                                  "duration_pre_abx"),
                                age_var_name = "age")

maled_results <- agaipw(data = one_hot_maled$data,
                        laz_var_name = "month3_haz",
                        abx_var_name = "all_abx", 
                        case_var_name = "case",
                        site_var_name = one_hot_maled$site_var_names,
                        followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
                        covariate_list = one_hot_maled$covariate_list,
                        pathogen_quantity_list = c("shigella_new",            
                                                   "adenovirus_40_41_new",                  
                                                   "aeromonas_new",                         
                                                   "astrovirus_new",                        
                                                   "campylobacter_pan_new",                 
                                                   "cryptosporidium_new",                   
                                                   "cyclospora_new",                        
                                                   "e_histolytica_new",                     
                                                   "isospora_new",                          
                                                   "norovirus_new",                         
                                                   "rotavirus_new",                         
                                                   "salmonella_new",                        
                                                   "sapovirus_new",                         
                                                   "st_etec_new",                           
                                                   "tEPEC_new",                             
                                                   "v_cholerae_new",                        
                                                   "ETEC_new",                              
                                                   "e_bieneusi_new",                        
                                                   "eaec_new"),
                        severity_list = one_hot_maled$severity_list,
                        no_etiology_var_name = "no_etiology",
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
                        case_control = TRUE,
                        first_id_var_name = "first_id")

saveRDS(maled_results, here::here("results/case_control/maled_msd_tac_or_culture_MSD.Rds"))
maled_results$aipw_est$aipw_models <- NULL
saveRDS(maled_results, here::here("results/case_control/maled_msd_tac_or_culture_MSD_no_models.Rds"))

# IPD ANALYSIS -------------------------------------------------------------

# MSD only
ipd_data <- readRDS(here::here("data/ipd_data/ipd_data_case_control_msd_tac_or_culture_shig.Rds"))

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

ipd_results_cc <- agaipw(data = one_hot_data$data,
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

saveRDS(ipd_results_cc, "results/case_control/ipd_results_msd_case_control_tac_or_culture.Rds")
ipd_results_cc$aipw_est$aipw_models <- NULL
saveRDS(ipd_results_cc, "results/case_control/ipd_results_msd_case_control_no_models_tac_or_culture.Rds")


# All cases (include LSD from MAL-ED) --------------------------------------------------------
ipd_data_all <- readRDS(here::here("data/ipd_data/ipd_data_case_control_tac_or_culture_shig.Rds"))

# Impute missing values, cases and controls separate
case_ipd_data <- ipd_data_all[ipd_data_all$case == 1,]
control_ipd_data <- ipd_data_all[ipd_data_all$case == 0,]

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


imp_data_case <- imp_data_case[,colnames(imp_data_case) %in% colnames(ipd_data_all), drop = FALSE]
imp_data_control <- imp_data_control[,colnames(imp_data_control) %in% colnames(ipd_data_all), drop = FALSE]
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

ipd_results_cc <- agaipw(data = one_hot_data$data,
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

saveRDS(ipd_results_cc, "results/case_control/ipd_results_case_control_tac_or_culture.Rds")
ipd_results_cc$aipw_est$aipw_models <- NULL
saveRDS(ipd_results_cc, "results/case_control/ipd_results_case_control_no_models_tac_or_culture.Rds")


