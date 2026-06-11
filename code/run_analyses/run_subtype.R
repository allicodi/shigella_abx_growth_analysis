# ------------------------------------------------------------
# Script to run subtype sub-infection analysis
# ------------------------------------------------------------

here::i_am("code/run_analyses/run_subtype.R")

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

## 1. No etiology comparison

# main dataset uses MSD MAL-ED only
og_data <- readRDS(here::here("data/ipd_data/ipd_data_no_etiology.Rds"))

# subset data to those with serotyping results
# ABCD - none - exclude bc used old serotyping genes that were not very good
# GEMS - culture serotyping - no NAs bc people who did not have culture positive shigella were still coded as 0
#                           - **fill in culture negative gems with NAs so not mixing with TAC pos, culture neg**
# VIDA - culture serotyping - only in people with culture shigella, people culture negative coded as NA
# MAL-ED - TAC serotyping - some people have shigella = no (both tac and culture), but S_flexneri yes? from new dataset. 
#                         - after much investigating i think bc they were coded no shigella via AFE but still had enough to be retested/typed
#                         - similarly there are people who do have shigella attributable via tac but did not get retested
# EFGH - TAC serotyping - NA if no shigella tac attr but all positives have serotype results

# notes/questions:
# removing ABCD altogether but probably could leave the shigella negative people in for the no etiology estimate?

# a few people have both shigella and flex. i think this is ok broadly (shig flex vs no etiology, shig sonnei vs no etiology, never flex vs sonnei)
# but the way we coded it you have to fall into one category... 
# 23 efgh, 4 gems, 5 vida
# remove for now? can hack around it later

# MAL-ED weirdness with shigella attributable no but still subtyped. i think this is fine too? just noting 

# feels weird to be mixing tac and culture definitions (ex. VIDA there are a bunch of tac positive people who are excluded bc only culture positive has results & coded as NA,
# vs GEMS tac positive were coded as 0 so being treated as third category?? but i don't know that we want to throw away the tac pos, culture negs, especially given we are using 
# tac for maled & efgh... so not sure if it'd make more sense to adjust for TAC missingness?? that's more work might have to wait for david for, or could i could try to mimic what we 
# did for pathogen missingness in dsldensify analysis but idk )

# Remove people who are TAC positive but no serotyping results (all of ABCD, a handful of MAL-ED, ~800 VIDA + ~800 GEMS)
data <- og_data[which(og_data$study != "ABCD"),]
data <- data[-which(data$shig_attr == 1 & (is.na(data$shigella_flex) | is.na(data$shigella_sonnei))), ]

# temp if both, remove (32 obs, mostly EFGH). probably should handle this though.
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

# temp if both, remove. probably should handle this though.
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
