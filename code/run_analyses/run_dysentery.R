# ------------------------------------------------------------
# Script to run dysentery sub-infection analysis
# ------------------------------------------------------------

here::i_am("code/run_analyses/run_dysentery.R")

devtools::load_all("../packages/abxGrowth")

library(SuperLearner)

source(here::here("code/run_analyses/SL.wrappers.R"))

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

sl.library.abcd.abx <- list(c("SL.glm", "SL.mean")) 

########################################################################################
# Main IPD results
########################################################################################

## 1. No etiology comparison -------------------------------------------------

# main dataset uses MSD MAL-ED only
data <- readRDS(here::here("data/ipd_data/ipd_data_no_etiology.Rds"))

# temp troubleshoot
# data <- data[data$study %in% c("MALED", "GEMS", "VIDA"),]

# NOTE - 15 people missing dysentery
# Exclude
data <- data[-which(is.na(data$dysentery)),]

imp_data <- impute_covariates(data = data,
                              imp_covariates = c("sex",
                                                 "age",
                                                 "ses_quintile",
                                                 "enr_haz",
                                                 "edu_bin",
                                                 "site",
                                                 "I_followup_days",
                                                 "I_followup_days_x_followup_days",
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
                               severity_list = c("any_vom", # removed dysentery here
                                                 "lsstools",
                                                 "dehyd_level",
                                                 "duration_pre_enroll",
                                                 "any_fev"),
                               age_var_name = "age")

# Create sub-infection variable (0 = no shigella, 1 = no dysentery, 2 = dysentery)
one_hot_data$data$sub_infection <- ifelse(one_hot_data$data$shig_attr == 0, 0,
                                          ifelse(one_hot_data$data$dysentery == 0, 1, 2))

# The following variables are missing in some studies
# water / sanitation - missing abcd - set to 0
# num children in hh - missing maled - set to 1
# education - missing in abcd - set to 0
# fever - missing in abcd - set to 0

dysentery_res <- agaipw(
  data = one_hot_data$data,
  laz_var_name = "final_haz",
  abx_var_name = "all_abx", 
  infection_var_name = "shig_attr",
  subinfection_var_name = "sub_infection",
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
  return_models = FALSE,
  msm = FALSE,
  subinfection_analysis = TRUE
)

saveRDS(dysentery_res, here::here("results/no_etiology/ipd_no_etiology_dysentery.Rds"))

### STUDY SPECIFIC RESULTS ###

efgh_data <- one_hot_data$data[which(one_hot_data$data$study_EFGH == 1), ]
gems_data <- one_hot_data$data[which(one_hot_data$data$study_GEMS == 1), ]
vida_data <- one_hot_data$data[which(one_hot_data$data$study_VIDA == 1), ]
maled_data <- one_hot_data$data[which(one_hot_data$data$study_MALED == 1), ]

efgh_res <- agaipw(
  data = efgh_data,
  laz_var_name = "final_haz",
  abx_var_name = "all_abx", 
  infection_var_name = "shig_attr",
  subinfection_var_name = "sub_infection",
  followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
  site_var_name = one_hot_data$site_var_names[grepl("_efgh", one_hot_data$site_var_names)],
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
  return_models = FALSE,
  msm = FALSE,
  subinfection_analysis = TRUE
)

saveRDS(efgh_res, here::here("results/no_etiology/efgh_no_etiology_dysentery.Rds"))

gems_res <- agaipw(
  data = gems_data,
  laz_var_name = "final_haz",
  abx_var_name = "all_abx", 
  infection_var_name = "shig_attr",
  subinfection_var_name = "sub_infection",
  followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
  site_var_name = one_hot_data$site_var_names[grepl("_gems", one_hot_data$site_var_names)],
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
  return_models = FALSE,
  msm = FALSE,
  subinfection_analysis = TRUE
)

saveRDS(gems_res, here::here("results/no_etiology/gems_no_etiology_dysentery.Rds"))

vida_res <- agaipw(
  data = vida_data,
  laz_var_name = "final_haz",
  abx_var_name = "all_abx", 
  infection_var_name = "shig_attr",
  subinfection_var_name = "sub_infection",
  followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
  site_var_name = one_hot_data$site_var_names[grepl("_vida", one_hot_data$site_var_names)],
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
  return_models = FALSE,
  msm = FALSE,
  subinfection_analysis = TRUE
)

saveRDS(vida_res, here::here("results/no_etiology/vida_no_etiology_dysentery.Rds"))

maled_res <- agaipw(
  data = maled_data,
  laz_var_name = "final_haz",
  abx_var_name = "all_abx", 
  infection_var_name = "shig_attr",
  subinfection_var_name = "sub_infection",
  followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
  site_var_name = one_hot_data$site_var_names[grepl("_maled", one_hot_data$site_var_names)],
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
  # sl.library.missingness = sl.library.with.abx,
  sl.library.missingness = c("SL.mean"),
  v_folds = 3, 
  return_models = FALSE,
  msm = FALSE,
  subinfection_analysis = TRUE
)

saveRDS(maled_res, here::here("results/no_etiology/maled_no_etiology_dysentery.Rds"))


## Try age stratification
data_0_11 <- one_hot_data$data[which(one_hot_data$data$age < 12), ]
data_12_23 <- one_hot_data$data[which(one_hot_data$data$age >= 12 & one_hot_data$data$age < 24), ]
data_24_59 <- one_hot_data$data[which(one_hot_data$data$age >= 24), ]

res_0_11 <- agaipw(
  data = data_0_11,
  laz_var_name = "final_haz",
  abx_var_name = "all_abx", 
  infection_var_name = "shig_attr",
  subinfection_var_name = "sub_infection",
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
  return_models = FALSE,
  msm = FALSE,
  subinfection_analysis = TRUE
)

saveRDS(res_0_11, here::here("results/no_etiology/ipd_no_etiology_dysentery_0_11.Rds"))

res_12_23 <- agaipw(
  data = data_12_23,
  laz_var_name = "final_haz",
  abx_var_name = "all_abx", 
  infection_var_name = "shig_attr",
  subinfection_var_name = "sub_infection",
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
  return_models = FALSE,
  msm = FALSE,
  subinfection_analysis = TRUE
)

saveRDS(res_12_23, here::here("results/no_etiology/ipd_no_etiology_dysentery_12_23.Rds"))

res_24_59 <- agaipw(
  data = data_24_59,
  laz_var_name = "final_haz",
  abx_var_name = "all_abx", 
  infection_var_name = "shig_attr",
  subinfection_var_name = "sub_infection",
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
  return_models = FALSE,
  msm = FALSE,
  subinfection_analysis = TRUE
)

saveRDS(res_24_59, here::here("results/no_etiology/ipd_no_etiology_dysentery_24_59.Rds"))


# 2. Case control comparison -------------------------------------------------

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
                               severity_list = c("any_vom", # removed dysentery
                                                 "lsstools",
                                                 "dehyd_level",
                                                 "duration_pre_enroll",
                                                 "any_fev"),
                               age_var_name = "age")

# Create sub-infection variable (0 = control (not case), 1 = no dysentery, 2 = dysentery)
one_hot_data$data$sub_infection <- ifelse(one_hot_data$data$case == 0, 0,
                                          ifelse(one_hot_data$data$dysentery == 0, 1, 2))

dysentery_cc_res <- agaipw(data = one_hot_data$data,
                         laz_var_name = "final_haz",
                         abx_var_name = "all_abx", 
                         case_var_name = "case",
                         subinfection_var_name = "sub_infection",
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
                         return_models = FALSE,
                         msm = FALSE,
                         case_control = TRUE, 
                         subinfection_analysis = TRUE)

saveRDS(dysentery_cc_res, here::here("results/case_control/ipd_cc_dysentery.Rds"))

### Study specific resu;ts

gems_data <- one_hot_data$data[which(one_hot_data$data$study_GEMS == 1), ]
vida_data <- one_hot_data$data[which(one_hot_data$data$study_VIDA == 1), ]
maled_data <- one_hot_data$data[which(one_hot_data$data$study_MALED == 1), ]

gems_cc_res <- agaipw(data = gems_data,
                       laz_var_name = "final_haz",
                       abx_var_name = "all_abx", 
                       case_var_name = "case",
                       subinfection_var_name = "sub_infection",
                       followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
                       site_var_name = one_hot_data$site_var_names[grepl("_gems", one_hot_data$site_var_names)],
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
                       return_models = FALSE,
                       msm = FALSE,
                       case_control = TRUE, 
                       subinfection_analysis = TRUE)

saveRDS(gems_cc_res, here::here("results/case_control/gems_cc_dysentery.Rds"))

vida_cc_res <- agaipw(data = vida_data,
                      laz_var_name = "final_haz",
                      abx_var_name = "all_abx", 
                      case_var_name = "case",
                      subinfection_var_name = "sub_infection",
                      followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
                      site_var_name = one_hot_data$site_var_names[grepl("_vida", one_hot_data$site_var_names)],
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
                      return_models = FALSE,
                      msm = FALSE,
                      case_control = TRUE, 
                      subinfection_analysis = TRUE)


maled_cc_res <- agaipw(data = maled_data,
                      laz_var_name = "final_haz",
                      abx_var_name = "all_abx", 
                      case_var_name = "case",
                      subinfection_var_name = "sub_infection",
                      followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days"),
                      site_var_name = one_hot_data$site_var_names[grepl("_maled", one_hot_data$site_var_names)],
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
                      return_models = FALSE,
                      msm = FALSE,
                      case_control = TRUE, 
                      subinfection_analysis = TRUE)

##########
# Age stratified 

data_0_11 <- one_hot_data$data[which(one_hot_data$data$age < 12), ]
data_12_23 <- one_hot_data$data[which(one_hot_data$data$age >= 12 & one_hot_data$data$age < 24), ]
data_24_59 <- one_hot_data$data[which(one_hot_data$data$age >= 24), ]

cc_res_0_11 <- agaipw(data = data_0_11,
                      laz_var_name = "final_haz",
                      abx_var_name = "all_abx", 
                      case_var_name = "case",
                      subinfection_var_name = "sub_infection",
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
                      return_models = FALSE,
                      msm = FALSE,
                      case_control = TRUE, 
                      subinfection_analysis = TRUE)

saveRDS(cc_res_0_11, here::here("results/case_control/ipd_cc_dysentery_0_11.Rds"))

cc_res_12_23 <- agaipw(data = data_12_23,
                      laz_var_name = "final_haz",
                      abx_var_name = "all_abx", 
                      case_var_name = "case",
                      subinfection_var_name = "sub_infection",
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
                      return_models = FALSE,
                      msm = FALSE,
                      case_control = TRUE, 
                      subinfection_analysis = TRUE)

saveRDS(cc_res_12_23, here::here("results/case_control/ipd_cc_dysentery_12_23.Rds"))

cc_res_24_59 <- agaipw(data = data_24_59,
                      laz_var_name = "final_haz",
                      abx_var_name = "all_abx", 
                      case_var_name = "case",
                      subinfection_var_name = "sub_infection",
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
                      return_models = FALSE,
                      msm = FALSE,
                      case_control = TRUE, 
                      subinfection_analysis = TRUE)

saveRDS(cc_res_24_59, here::here("results/case_control/ipd_cc_dysentery_24_59.Rds"))
