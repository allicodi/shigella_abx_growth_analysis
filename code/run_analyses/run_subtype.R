# ------------------------------------------------------------
# Script to run subtype sub-infection analysis
# ------------------------------------------------------------

here::i_am("code/run_analyses/run_subtype.R")

# For developmental version locally
#devtools::load_all("../packages/abxGrowth/")

# From GitHub (installed using devtools::install_github("allicodi/abxGrowth"))
library(abxGrowth)

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

## 1. No etiology comparison

# main dataset uses MSD MAL-ED only
og_data <- readRDS(here::here("data/ipd_data/ipd_data_no_etiology.Rds"))

# Remove people who are TAC positive but no serotyping results (all of ABCD, a handful of MAL-ED, ~800 VIDA + ~800 GEMS)
data <- og_data[which(og_data$study != "ABCD"),]
data <- data[-which(data$shig_attr == 1 & (is.na(data$shigella_flex) | is.na(data$shigella_sonnei))), ]

# if both, remove (32 obs, mostly EFGH)
data <- data[-which(data$shigella_flex == 1 & data$shigella_sonnei == 1),]

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
                                                 "dysentery",
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
                               severity_list = c("any_vom", 
                                                 "dysentery",
                                                 "lsstools",
                                                 "dehyd_level",
                                                 "duration_pre_enroll",
                                                 "any_fev"),
                               age_var_name = "age")

# Create sub-infection variable (0 = no shigella, 1 = flexneri, 2 = sonnei, 3 = neither)
one_hot_data$data$sub_infection <- ifelse(one_hot_data$data$shig_attr == 0, 0,
                                          ifelse(one_hot_data$data$shigella_flex == 1, 1,
                                                 ifelse(one_hot_data$data$shigella_sonnei == 1, 2, 3)))

# The following variables are missing in some studies
# water / sanitation - missing abcd - set to 0
# num children in hh - missing maled - set to 1
# education - missing in abcd - set to 0
# fever - missing in abcd - set to 0

serotype_res <- agaipw(
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
  subinfection_analysis = TRUE, 
  complete_subinfection = FALSE 
)


saveRDS(serotype_res, here::here("results/no_etiology/ipd_no_etiology_serotype.Rds"))

# 2. Case Control -----------------------------------------------------------------

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

# Remove people who are cases but no serotyping results (a handful of MAL-ED, ~800 GEMS, ~800 VIDA)
data <- ipd_data[-which(ipd_data$case == 1 & (is.na(ipd_data$shigella_flex) | is.na(ipd_data$shigella_sonnei))), ]

# if both, remove
data <- data[-which(data$shigella_flex == 1 & data$shigella_sonnei == 1),]


# Impute missing values, cases and controls separate
case_ipd_data <- data[data$case == 1,]
control_ipd_data <- data[data$case == 0,]

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
                               severity_list = c("any_vom", 
                                                 "dysentery",
                                                 "lsstools",
                                                 "dehyd_level",
                                                 "duration_pre_enroll",
                                                 "any_fev"),
                               age_var_name = "age")

# make subinfection variable (0 = control, 1 = flex, 2 = sonnei, 3 = other)
one_hot_data$data$sub_infection <- ifelse(one_hot_data$data$case == 0, 0,
                                          ifelse(one_hot_data$data$shigella_flex == 1, 1,
                                                 ifelse(one_hot_data$data$shigella_sonnei == 1, 2, 3)))

subtype_cc_res <- agaipw(data = one_hot_data$data,
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
                           subinfection_analysis = TRUE, 
                           complete_subinfection = FALSE)

saveRDS(subtype_cc_res, here::here("results/case_control/ipd_cc_serotype.Rds"))
