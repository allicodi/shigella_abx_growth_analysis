# ------------------------------------------------------------------------------
# Script to add EED analyses based on Nature Communications review #1 comments
# ------------------------------------------------------------------------------
# 1. exposure = shigella cases + shigella detected under 30 or 35 vs no etiology + no shigella detected below threshold for exposure
# study specific estimates + overall estimate
# reuse existing pipeline, redefine exposure
# 
# 2. controls only -- shigella detected (30 or 35 cutoffs) vs no shigella detection
# 
# study specific estimates (for the studies that have it) + overall estimate
# use drtmle package rather than existing pipeline
# ------------------------------------------------------------------------------

here::i_am("code/run_analyses/run_EED.R")

#library(abxGrowth)
devtools::load_all("../packages/abxGrowth")

library(SuperLearner)
source(here::here("code/run_analyses/SL.wrappers.R"))

# Comparison to other diarrhea
gems_data <- readRDS(here::here("data/gems_data/gems_data.Rds"))
# maled_data <- readRDS(here::here("data/maled_data/maled_data.Rds"))
maled_data_MSD <- readRDS(here::here("data/maled_data/maled_data_shig_MSD_only.Rds"))
vida_data <- readRDS(here::here("data/vida_data/vida_data.Rds"))
efgh_data <- readRDS(here::here("data/efgh_data/efgh_data.Rds"))
abcd_data <- readRDS(here::here("data/abcd_data/abcd_data.Rds"))

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


sl.library.abcd.abx <- list("SL.glm", "SL.mean")

# --------------------------------------------------------------------------------
# 1a. Exposure = Diarrhea cases with any shigella detected; Comparator = Diarrhea cases with no Shigella at all
# --------------------------------------------------------------------------------

## ABCD: shigella_bin = shigella detected; comparator = no shigella detected (!shigella_bin)

abcd_data_select <- abcd_data[,colnames(abcd_data) %in% c("pid",
                                                          "site",
                                                          "dy1_ant_sex",
                                                          "agemchild",
                                                          "lfazscore",
                                                          "an_ses_quintile",
                                                          "an_tothhlt5",
                                                          "I_an_d90_timing",
                                                          "I_an_d90_timing_x_an_d90_timing",
                                                          "rotavirus_new",                  
                                                          "norovirus_gii_new",                  
                                                          "adenovirus_new",                 
                                                          "astrovirus_new",                 
                                                          "sapovirus_new",                  
                                                          "st_etec_new",                    
                                                          "shigella_new",                   
                                                          "campylobacter_new",              
                                                          "tepec_new",                      
                                                          "v_cholerae_new",                 
                                                          "salmonella_new",
                                                          "an_grp_01",
                                                          "shigella_likely",
                                                          "lazd90",
                                                          "no_etiology",
                                                          "shigella_bin")]

abcd_data_imp <- impute_covariates(abcd_data_select, 
                                   imp_covariates = c("site",
                                                      "dy1_ant_sex",
                                                      "agemchild",
                                                      "lfazscore",
                                                      "an_ses_quintile",
                                                      "an_tothhlt5",
                                                      "I_an_d90_timing",
                                                      "I_an_d90_timing_x_an_d90_timing",
                                                      "rotavirus_new",                  
                                                      "norovirus_gii_new",                  
                                                      "adenovirus_new",                 
                                                      "astrovirus_new",                 
                                                      "sapovirus_new",                  
                                                      "st_etec_new",                    
                                                      "shigella_new",                   
                                                      "campylobacter_new",              
                                                      "tepec_new",                      
                                                      "v_cholerae_new",                 
                                                      "salmonella_new"),
                                   site_var_name = "site")

abcd_data_imp$age <- abcd_data_imp$agemchild
abcd_data_imp$enr_haz <- abcd_data_imp$lfazscore
abcd_data_imp$guideline_abx <- factor(abcd_data_imp$an_grp_01, levels = 0:1, labels = c("No azithromycin", "Azithromycin"))

one_hot_abcd <- one_hot_encode(data = abcd_data_imp,
                               laz_var_name = "lazd90",
                               covariate_list = c("site",
                                                  "dy1_ant_sex",
                                                  "age",
                                                  "enr_haz",
                                                  "an_ses_quintile",
                                                  "an_tothhlt5",
                                                  "I_an_d90_timing",
                                                  "I_an_d90_timing_x_an_d90_timing"),
                               abx_var_name = "guideline_abx",
                               site_var_name = "site",
                               site_interaction = TRUE,
                               age_var_name = "age")

# "no_etiology" comparator in package = !shigella_bin
one_hot_abcd$data$shigella_not_detected <- as.numeric(!one_hot_abcd$data$shigella_bin)

abcd_results <- agaipw(data = one_hot_abcd$data,
                       laz_var_name = "lazd90",
                       abx_var_name = "guideline_abx",
                       infection_var_name = "shigella_bin",
                       site_var_name = one_hot_abcd$site_var_names,
                       followup_var_names = c("I_an_d90_timing","I_an_d90_timing_x_an_d90_timing"),
                       covariate_list = one_hot_abcd$covariate_list,
                       pathogen_quantity_list = c("rotavirus_new",                  
                                                  "norovirus_gii_new",                  
                                                  "adenovirus_new",                 
                                                  "astrovirus_new",                 
                                                  "sapovirus_new",                  
                                                  "st_etec_new",                    
                                                  "shigella_new",                   
                                                  "campylobacter_new",              
                                                  "tepec_new",                      
                                                  "v_cholerae_new",                 
                                                  "salmonella_new"),
                       no_etiology_var_name = "shigella_not_detected",
                       first_id_var_name = "pid",
                       outcome_type = "gaussian",
                       sl.library.outcome = sl.library.with.abx,
                       sl.library.outcome.2 = sl.library.without.abx,
                       sl.library.treatment = sl.library.abcd.abx,
                       sl.library.infection = sl.library.without.abx,
                       sl.library.missingness = sl.library.with.abx,
                       v_folds = 5, 
                       return_models = TRUE,
                       msm = TRUE,
                       msm_var_name = "age",
                       msm_formula = "age", 
                       all_other_diarrhea = TRUE)

abcd_results$aipw_est$aipw_models <- NULL
saveRDS(abcd_results, here::here("results/eed/case_only/abcd_results_35_no_models.Rds"))

## EFGH: create shigella_bin

efgh_data$shigella_bin <- ifelse(efgh_data$shigella_ct < 35, 1, 0)
efgh_data$shigella_not_detected <- as.numeric(!efgh_data$shigella_bin)

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
                                                      "EAEC_new",
                                                      "shigella_bin",
                                                      "shigella_not_detected"))

efgh_data_imp$age <- efgh_data_imp$enr_age_months
efgh_data_imp <- efgh_data_imp[which(!is.na(efgh_data_imp$tac_shigella_attributable)),]
efgh_data_imp$guideline_abx <- efgh_data_imp$who_rec_abx

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
                               abx_var_name = "guideline_abx",
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

efgh_results <- agaipw(data = one_hot_efgh$data,
                       laz_var_name = "mo3_haz",
                       abx_var_name = "guideline_abx", 
                       infection_var_name = "shigella_bin",
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
                       no_etiology_var_name = "shigella_not_detected",
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
                       all_other_diarrhea = TRUE)

efgh_results$aipw_est$aipw_models <- NULL
saveRDS(efgh_results, here::here("results/eed/case_only/efgh_results_35_no_models.Rds"))

## GEMS: create shigella_bin

gems_data$shigella_bin <- ifelse(gems_data$shigella_new > 0, 1, 0)
gems_data$shigella_not_detected <- as.numeric(!gems_data$shigella_bin)

gems_data$guideline_abx <- gems_data$who_abx

# Impute missing covariates
gems_data_imp <- impute_covariates(gems_data,
                                   imp_covariates = c("sex",
                                                      "age",
                                                      "ses_quintile",
                                                      "safe_water",
                                                      "safe_sanit",
                                                      "enr_haz",
                                                      "prim_caregiver_edu_bin",
                                                      "num_hh_lt5",
                                                      "rotavirus_attributable",
                                                      "adenovirus_attributable",
                                                      "aeromonas_attributable",
                                                      "astro_attributable",
                                                      "cryptosporidium_attributable",
                                                      "cyclospora_attributable",
                                                      "e_histolytica_attributable",
                                                      "isospora_attributable",
                                                      "noro_attributable",
                                                      "salmonella_attributable",
                                                      "sapovirus_attributable",
                                                      "st_etec_attributable",
                                                      "tepec_attributable",
                                                      "v_cholerae_attributable",
                                                      "EAEC_attributable",
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
                               abx_var_name = "guideline_abx",
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
                       abx_var_name = "guideline_abx", 
                       infection_var_name = "shigella_bin",
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
                       no_etiology_var_name = "shigella_not_detected",
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
                       all_other_diarrhea = TRUE)

gems_results$aipw_est$aipw_models <- NULL
saveRDS(gems_results, here::here("results/eed/case_only/gems_results_35_no_models.Rds"))

## VIDA

vida_data$shigella_bin <- ifelse(vida_data$shigella_new > 0, 1, 0)
vida_data$shigella_not_detected <- as.numeric(!vida_data$shigella_bin)

vida_data_imp <- impute_covariates(data = vida_data,
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
                                                      "I_followup_days",
                                                      "I_followup_days_x_followup_days",
                                                      "shigella_new",                   
                                                      "rotavirus_new",                  
                                                      "st_etec_new",                    
                                                      "crypto_new",                     
                                                      "adeno_new",                      
                                                      "astro_new",                      
                                                      "noro_gii_new",                       
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

vida_data_imp$age <- vida_data_imp$agemchild
vida_data_imp$guideline_abx <- vida_data_imp$who_abx

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
                               abx_var_name = "guideline_abx",
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
                       abx_var_name = "guideline_abx", 
                       infection_var_name = "shigella_bin",
                       followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
                       site_var_name = one_hot_vida$site_var_names,
                       covariate_list = one_hot_vida$covariate_list,
                       pathogen_quantity_list = c("shigella_new",                   
                                                  "rotavirus_new",                  
                                                  "st_etec_new",                    
                                                  "crypto_new",                     
                                                  "adeno_new",                      
                                                  "astro_new",                      
                                                  "noro_gii_new",                       
                                                  "tepec_new",                      
                                                  "campy_new",                      
                                                  "sapo_new",                       
                                                  "giardia_new",                    
                                                  "eaec_new"),
                       no_etiology_var_name = "shigella_not_detected",
                       first_id_var_name = "first_id",
                       severity_list = one_hot_vida$severity_list,
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
                       all_other_diarrhea = TRUE)

vida_results$aipw_est$aipw_models <- NULL
saveRDS(vida_results, here::here("results/eed/case_only/vida_results_35_no_models.Rds"))

## MAL-ED

maled_data_MSD$shigella_bin <- ifelse(maled_data_MSD$shigella_new > 0, 1, 0)
maled_data_MSD$shigella_not_detected <- as.numeric(!maled_data_MSD$shigella_bin)

maled_data_imp_MSD <- impute_covariates(data = maled_data_MSD,
                                        site_var_name = "site",
                                        imp_by_site = TRUE,
                                        imp_covariates = c("site",
                                                           "sex",
                                                           "agemonths",
                                                           "baseline_haz",
                                                           "wami_quintile",
                                                           "mated_bin", 
                                                           "drinkimp",
                                                           "sanitimp",
                                                           "adenovirus_40_41_new", 
                                                           "shigella_new",
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
                                                           "eaec_new",
                                                           "dysentery_g",
                                                           "fever_g",
                                                           "fever_days_g",
                                                           "dehyd_g",
                                                           "lsstools_g",
                                                           "daysvomit_g"))

# Add "age" and "enr_haz" to match wrappers (go back and rename in raw data later)
# and guideline_abx
maled_data_imp_MSD$age <- maled_data_imp_MSD$agemonths
maled_data_imp_MSD$enr_haz <- maled_data_imp_MSD$baseline_haz
maled_data_imp_MSD$guideline_abx <- maled_data_imp_MSD$who_abx

one_hot_maled_MSD <- one_hot_encode(data = maled_data_imp_MSD,
                                    laz_var_name = "month3_haz",
                                    covariate_list = c("site",
                                                       "sex",
                                                       "age",
                                                       "wami_quintile",
                                                       "enr_haz",
                                                       "mated_bin", 
                                                       "drinkimp",
                                                       "sanitimp",
                                                       "I_followup_days",
                                                       "I_followup_days_x_followup_days"),
                                    abx_var_name = "guideline_abx",
                                    site_var_name = "site",
                                    site_interaction = TRUE,
                                    severity_list = c("dysentery_g",
                                                      "fever_g",
                                                      "fever_days_g",
                                                      "dehyd_g",
                                                      "lsstools_g",
                                                      "daysvomit_g",
                                                      "duration_pre_abx_g"),
                                    age_var_name = "age")

maled_results_MSD <- agaipw(data = one_hot_maled_MSD$data,
                            laz_var_name = "month3_haz",
                            abx_var_name = "guideline_abx", 
                            infection_var_name = "shigella_bin",
                            site_var_name = one_hot_maled_MSD$site_var_names,
                            followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
                            covariate_list = one_hot_maled_MSD$covariate_list,
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
                            severity_list = one_hot_maled_MSD$severity_list,
                            no_etiology_var_name = "shigella_not_detected",
                            outcome_type = "gaussian",
                            sl.library.outcome = sl.library.with.abx,
                            sl.library.outcome.2 = sl.library.without.abx,
                            sl.library.treatment = sl.library.without.abx,
                            sl.library.infection = sl.library.without.abx,
                            sl.library.missingness = c("SL.mean"), #  mean in MSD only version bc only one missing obs
                            v_folds = 5, 
                            return_models = TRUE,
                            msm = TRUE,
                            msm_var_name = "age",
                            msm_formula = "age",
                            first_id_var_name = "pid", 
                            all_other_diarrhea = TRUE)

maled_results_MSD$aipw_est$aipw_models <- NULL
saveRDS(maled_results_MSD, here::here("results/eed/case_only/maled_msd_results_35_no_models.Rds"))

## IPD

# main dataset uses MSD MAL-ED only
data <- readRDS(here::here("data/ipd_data/ipd_data_no_etiology_bin_abx.Rds"))

data$shigella_bin <- ifelse(data$shigella_new > 0 , 1, 0)
data$shigella_not_detected <- as.numeric(!data$shigella_bin)

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
                      infection_var_name = "shigella_bin",
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
                      no_etiology_var_name = "shigella_not_detected",
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
                      all_other_diarrhea = TRUE)

ipd_results$aipw_est$aipw_models <- NULL
saveRDS(ipd_results, here::here("results/eed/case_only/ipd_results_35_no_models.Rds"))

# -----------------------------------------------------------------------------

# 1b. Exposure = Diarrhea cases with shigella TAC < 30 detected; Comparator = Diarrhea cases with no Shigella at all
# Hold off for now?? cutoff is ~30 anyways for studies that have it so do we really care abt this?

# -----------------------------------------------------------------------------

# 2a. Exposure = Controls with any Shigella detected; Comparator = Controls with no Shigella detected

recreate_augmentation_matrix <- function(first_id, aug_matrix){
  aug_id_matrix <- data.frame(first_id = first_id, 
                              eif = aug_matrix)
  
  aug_id_matrix <- aggregate(. ~ first_id, data = aug_id_matrix, FUN = sum)
  
  # scale by nrow(unique kids) / nrow(full data)
  new_augmentation <- aug_id_matrix[,2:ncol(aug_id_matrix)] * (nrow(aug_id_matrix) / length(first_id))
  
  return(new_augmentation)
}

## PRIMARY ANALYSIS: Case definition = TAC or culture attributable Shigella
gems_controls <- readRDS(here::here("data/gems_data/gems_case_control_tac_or_culture_shig.Rds")) %>%
  filter(case == 0)
vida_controls <- readRDS(here::here("data/vida_data/vida_case_control_tac_or_culture_shig.Rds")) %>%
  filter(case == 0)
maled_data_tac_culture <- readRDS(here::here("data/maled_data/maled_case_control_tac_or_culture_shig.Rds")) 

# MAL-ED MSD only
maled_msd_tac_culture <- maled_data_tac_culture
maled_data_msd_caseids <- unique(maled_msd_tac_culture$case_id[which(maled_msd_tac_culture$MSD == 1)])
maled_msd_tac_culture <- maled_msd_tac_culture[which(maled_msd_tac_culture$case_id %in% maled_data_msd_caseids),]
maled_controls <- maled_msd_tac_culture %>%
  filter(case == 0)

## GEMS

gems_controls$shigella_bin <- ifelse(gems_controls$shigella_new > 0, 1, 0)
gems_controls$shigella_levels <- ifelse(gems_controls$shigella_new == 0, 0, 
                                        ifelse(gems_controls$shigella_new < (35 - 30) / 3.322, 1, 2))

gems_control_data_imp <- impute_covariates(gems_controls,
                                           imp_covariates = c("sex",
                                                              "age",
                                                              "ses_quintile",
                                                              "safe_water",
                                                              "safe_sanit",
                                                              "enr_haz",
                                                              "prim_caregiver_edu_bin",
                                                              "num_hh_lt5"))

one_hot_gems_controls <- one_hot_encode(gems_control_data_imp,
                               laz_var_name = "hazd60",
                               covariate_list_control = c("sex",
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
                               abx_var_name = "who_abx",
                               site_var_name = "site",
                               site_interaction = TRUE,
                               age_var_name = "age")

gems_control_res <- drtmle(
  Y = one_hot_gems_controls$data$hazd60,
  A = one_hot_gems_controls$data$shigella_bin,
  W = one_hot_gems_controls$data[,colnames(one_hot_gems_controls$data) %in% one_hot_gems_controls$covariate_list_control],
  a_0 = c(0,1),
  family = gaussian(),
  stratify = TRUE,
  avg_over = "SL",
  SL_Q = sl.library.without.abx,
  SL_g = sl.library.without.abx,
  SL_Qr = sl.library.without.abx,
  SL_gr = sl.library.without.abx,
  guard = NULL,
  targeted_se = FALSE
)

# Unadjusted for multiple enrollments
# gems_control_ci <- ci(gems_control_res, est = "aiptw", contrast = c(-1, 1))

# Adjust for multiple enrollments
gems_new_aug <- recreate_augmentation_matrix(first_id = one_hot_gems_controls$data$first_id,
                                             aug_matrix = gems_control_res$ic_drtmle$ic)

# Create CI adjusting for multiple enrollments
gems_control_res_adjusted <- gems_control_res
gems_control_res_adjusted$aiptw$cov <- stats::cov(gems_new_aug) / dim(gems_new_aug)[1]
gems_control_adjusted_ci <- ci(gems_control_res_adjusted, est = "aiptw", contrast = c(-1, 1))

saveRDS(gems_control_res_adjusted, here::here("results/eed/control_only/gems_control_res_35.Rds"))
saveRDS(gems_control_adjusted_ci, here::here("results/eed/control_only/gems_control_ci_35.Rds"))

# VIDA

vida_controls$shigella_bin <- ifelse(vida_controls$shigella_new > 0, 1, 0)
vida_controls$shigella_levels <- ifelse(vida_controls$shigella_new == 0, 0, 
                                        ifelse(vida_controls$shigella_new < (35 - 30) / 3.322, 1, 2))

vida_controls$age <- vida_controls$agemchild

vida_control_data_imp <- impute_covariates(vida_controls,
                                           imp_covariates = c("sex",
                                                              "age",
                                                              "ses_quintile",
                                                              "safe_water",
                                                              "safe_sanit",
                                                              "enr_haz",
                                                              "site",
                                                              "education_bin",
                                                              "num_hh_lt5"))

one_hot_vida_controls <- one_hot_encode(vida_control_data_imp,
                                        laz_var_name = "hazd60",
                                        covariate_list_control = c("sex",
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
                                        abx_var_name = "who_abx",
                                        site_var_name = "site",
                                        site_interaction = TRUE,
                                        age_var_name = "age")

vida_control_res <- drtmle(
  Y = one_hot_vida_controls$data$hazd60,
  A = one_hot_vida_controls$data$shigella_bin,
  W = one_hot_vida_controls$data[,colnames(one_hot_vida_controls$data) %in% one_hot_vida_controls$covariate_list_control],
  a_0 = c(0,1),
  family = gaussian(),
  stratify = TRUE,
  avg_over = "SL",
  SL_Q = sl.library.without.abx,
  SL_g = sl.library.without.abx,
  SL_Qr = sl.library.without.abx,
  SL_gr = sl.library.without.abx,
  guard = NULL,
  targeted_se = FALSE
)

# vida_control_ci <- ci(vida_control_res, est = "aiptw", contrast = c(-1, 1))

# Adjust for multiple enrollments
vida_new_aug <- recreate_augmentation_matrix(first_id = one_hot_vida_controls$data$first_id,
                                             aug_matrix = vida_control_res$ic_drtmle$ic)

# Create CI adjusting for multiple enrollments
vida_control_res_adjusted <- vida_control_res
vida_control_res_adjusted$aiptw$cov <- stats::cov(vida_new_aug) / dim(vida_new_aug)[1]
vida_control_adjusted_ci <- ci(vida_control_res_adjusted, est = "aiptw", contrast = c(-1, 1))

saveRDS(vida_control_res_adjusted, here::here("results/eed/control_only/vida_control_res_35.Rds"))
saveRDS(vida_control_adjusted_ci, here::here("results/eed/control_only/vida_control_ci_35.Rds"))

## MAL-ED

maled_controls$shigella_bin <- ifelse(maled_controls$shigella_new > 0, 1, 0)
maled_controls$shigella_levels <- ifelse(maled_controls$shigella_new == 0, 0, 
                                        ifelse(maled_controls$shigella_new < (35 - 30) / 3.322, 1, 2))

maled_controls$age <- maled_controls$agemonths
maled_controls$enr_haz <- maled_controls$baseline_haz

maled_control_data_imp <- impute_covariates(data = maled_controls,
                                            site_var_name = "site",
                                            imp_by_site = TRUE,
                                            imp_covariates = c("site",
                                                               "sex",
                                                               "age",
                                                               "enr_haz",
                                                               "wami_quintile",
                                                               "mated_bin", 
                                                               "drinkimp",
                                                               "sanitimp"))

one_hot_maled_controls <- one_hot_encode(data = maled_control_data_imp,
                                laz_var_name = "month3_haz",
                                covariate_list_control = c("site",
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
                                age_var_name = "age")

maled_control_res <- drtmle(
  Y = one_hot_maled_controls$data$month3_haz,
  A = one_hot_maled_controls$data$shigella_bin,
  W = one_hot_maled_controls$data[,colnames(one_hot_maled_controls$data) %in% one_hot_maled_controls$covariate_list_control],
  a_0 = c(0,1),
  family = gaussian(),
  stratify = TRUE,
  avg_over = "SL",
  SL_Q = sl.library.without.abx,
  SL_g = sl.library.without.abx,
  SL_Qr = sl.library.without.abx,
  SL_gr = sl.library.without.abx,
  guard = NULL,
  targeted_se = FALSE
)

# Unadjusted
# maled_control_ci <- ci(maled_control_res, est = "aiptw", contrast = c(-1, 1))

# Adjust for multiple enrollments
maled_new_aug <- recreate_augmentation_matrix(first_id = one_hot_maled_controls$data$first_id,
                                             aug_matrix = maled_control_res$ic_drtmle$ic)

# Create CI adjusting for multiple enrollments
maled_control_res_adjusted <- maled_control_res
maled_control_res_adjusted$aiptw$cov <- stats::cov(maled_new_aug) / dim(maled_new_aug)[1]
maled_control_adjusted_ci <- ci(maled_control_res_adjusted, est = "aiptw", contrast = c(-1, 1))

saveRDS(maled_control_res_adjusted, here::here("results/eed/control_only/maled_control_res_35.Rds"))
saveRDS(maled_control_adjusted_ci, here::here("results/eed/control_only/maled_control_ci_35.Rds"))

## IPD 

# MSD only
ipd_data <- readRDS(here::here("data/ipd_data/ipd_data_case_control_msd_tac_or_culture_shig.Rds"))

ipd_controls <- ipd_data[ipd_data$case == 0,]

ipd_controls$shigella_bin <- ifelse(ipd_controls$shigella_new > 0, 1, 0)
ipd_controls$shigella_levels <- ifelse(ipd_controls$shigella_new == 0, 0, 
                                   ifelse(ipd_controls$shigella_new < (35 - 30) / 3.322, 1, 2))


imp_data_control <- impute_covariates(data = ipd_controls,
                                      imp_covariates = c("sex",
                                                         "age",
                                                         "ses_quintile",
                                                         "enr_haz",
                                                         "edu_bin",
                                                         "site",
                                                         "I_followup_days",
                                                         "I_followup_days_x_followup_days"))

one_hot_data <- one_hot_encode(imp_data_control,
                               laz_var_name = "final_haz",
                               covariate_list_control = c("sex",
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
                               age_var_name = "age")

ipd_control_res <- drtmle(
  Y = one_hot_data$data$final_haz,
  A = one_hot_data$data$shigella_bin,
  W = one_hot_data$data[,colnames(one_hot_data$data) %in% one_hot_data$covariate_list_control],
  a_0 = c(0,1),
  family = gaussian(),
  stratify = TRUE,
  avg_over = "SL",
  SL_Q = sl.library.without.abx,
  SL_g = sl.library.without.abx,
  SL_Qr = sl.library.without.abx,
  SL_gr = sl.library.without.abx,
  guard = NULL,
  targeted_se = FALSE
)

# Unadjusted
# ipd_control_ci <- ci(ipd_control_res, est = "aiptw", contrast = c(-1, 1))

# Adjust for multiple enrollments
ipd_new_aug <- recreate_augmentation_matrix(first_id = one_hot_data$data$first_id,
                                              aug_matrix = ipd_control_res$ic_drtmle$ic)

# Create CI adjusting for multiple enrollments
ipd_control_res_adjusted <- ipd_control_res
ipd_control_res_adjusted$aiptw$cov <- stats::cov(ipd_new_aug) / dim(ipd_new_aug)[1]
ipd_control_adjusted_ci <- ci(ipd_control_res_adjusted, est = "aiptw", contrast = c(-1, 1))

saveRDS(ipd_control_res_adjusted, here::here("results/eed/control_only/ipd_control_res_35.Rds"))
saveRDS(ipd_control_adjusted_ci, here::here("results/eed/control_only/ipd_control_ci_35.Rds"))

# -----------------------------------------------------------------------------
# 2b. Exposure = Controls with Shigella TAC < 30 detected; Comparator = Controls with no Shigella detected

gems_control_res_30 <- drtmle(
  Y = one_hot_gems_controls$data$hazd60,
  A = one_hot_gems_controls$data$shigella_levels,
  W = one_hot_gems_controls$data[,colnames(one_hot_gems_controls$data) %in% one_hot_gems_controls$covariate_list_control],
  a_0 = c(0,1,2),
  family = gaussian(),
  stratify = TRUE,
  avg_over = "SL",
  SL_Q = sl.library.without.abx,
  SL_g = sl.library.without.abx,
  SL_Qr = sl.library.without.abx,
  SL_gr = sl.library.without.abx,
  guard = NULL,
  targeted_se = FALSE
)

vida_control_res_30 <- drtmle(
  Y = one_hot_vida_controls$data$hazd60,
  A = one_hot_vida_controls$data$shigella_levels,
  W = one_hot_vida_controls$data[,colnames(one_hot_vida_controls$data) %in% one_hot_vida_controls$covariate_list_control],
  a_0 = c(0,1,2),
  family = gaussian(),
  stratify = TRUE,
  avg_over = "SL",
  SL_Q = sl.library.without.abx,
  SL_g = sl.library.without.abx,
  SL_Qr = sl.library.without.abx,
  SL_gr = sl.library.without.abx,
  guard = NULL,
  targeted_se = FALSE
)

maled_control_res_30 <- drtmle(
  Y = one_hot_maled_controls$data$month3_haz,
  A = one_hot_maled_controls$data$shigella_levels,
  W = one_hot_maled_controls$data[,colnames(one_hot_maled_controls$data) %in% one_hot_maled_controls$covariate_list_control],
  a_0 = c(0,1,2),
  family = gaussian(),
  stratify = TRUE,
  avg_over = "SL",
  SL_Q = sl.library.without.abx,
  SL_g = sl.library.without.abx,
  SL_Qr = sl.library.without.abx,
  SL_gr = sl.library.without.abx,
  guard = NULL,
  targeted_se = FALSE
)

ipd_control_res_30 <- drtmle(
  Y = one_hot_data$data$final_haz,
  A = one_hot_data$data$shigella_levels,
  W = one_hot_data$data[,colnames(one_hot_data$data) %in% one_hot_data$covariate_list_control],
  a_0 = c(0,1,2),
  family = gaussian(),
  stratify = TRUE,
  avg_over = "SL",
  SL_Q = sl.library.without.abx,
  SL_g = sl.library.without.abx,
  SL_Qr = sl.library.without.abx,
  SL_gr = sl.library.without.abx,
  guard = NULL,
  targeted_se = FALSE
)

# -----------------------------------------------------------------------------
# Adjust for multiple enrollments and calculate CIs
# -----------------------------------------------------------------------------

control_res_30 <- list(
  gems = gems_control_res_30,
  vida = vida_control_res_30,
  maled = maled_control_res_30,
  ipd = ipd_control_res_30
)

control_data_30 <- list(
  gems = one_hot_gems_controls$data,
  vida = one_hot_vida_controls$data,
  maled = one_hot_maled_controls$data,
  ipd = one_hot_data$data
)

for (study in names(control_res_30)) {
  
  # Adjust for multiple enrollments
  new_aug <- recreate_augmentation_matrix(
    first_id = control_data_30[[study]]$first_id,
    aug_matrix = control_res_30[[study]]$ic_drtmle$ic
  )
  
  # Create adjusted covariance matrix
  res_adjusted <- control_res_30[[study]]
  res_adjusted$aiptw$cov <- stats::cov(new_aug) / nrow(new_aug)
  
  # Calculate CI for Ct < 30 vs no Shigella detection
  adjusted_ci <- ci(
    res_adjusted,
    est = "aiptw",
    contrast = c(-1, 0, 1)
  )
  
  # Save adjusted results
  saveRDS(
    res_adjusted,
    here::here(paste0(
      "results/eed/control_only/",
      study, "_control_res_30.Rds"
    ))
  )
  
  # Save adjusted CI
  saveRDS(
    adjusted_ci,
    here::here(paste0(
      "results/eed/control_only/",
      study, "_control_ci_30.Rds"
    ))
  )
}