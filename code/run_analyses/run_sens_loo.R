# ------------------------------------------------------------------
# Run leave-one-out study sensitivity analyses of IPD data
# ------------------------------------------------------------------

here::i_am("code/run_analyses/run_sens_loo.R")

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

data <- readRDS(here::here("data/ipd_data/ipd_data_no_etiology_bin_abx.Rds"))

# ------------------------------------------------------------------
# Run leave-one-study-out analyses - no etiology
# ------------------------------------------------------------------

studies <- c("ABCD", "MALED", "GEMS", "VIDA", "EFGH")

# Check study variable and values before running
table(data$study)

for (excluded_study in studies) {
  
  message("Running IPD analysis excluding: ", excluded_study)
  
  # Exclude one study
  loo_data <- data[data$study != excluded_study, , drop = FALSE]
  
  # Impute covariates
  imp_data <- impute_covariates(
    data = loo_data,
    imp_covariates = c(
      "sex", "age", "ses_quintile", "enr_haz",
      "edu_bin", "site", "I_followup_days",
      "I_followup_days_x_followup_days",
      "dysentery", "any_vom", "any_fev",
      "dehyd_level", "lsstools", "duration_pre_enroll",
      "shigella_new", "rotavirus_new",
      "adenovirus_new", "etec_new",
      "cryptosporidium_new", "astrovirus_new",
      "norovirus_gii_new", "tepec_new",
      "campylobacter_new", "sapovirus_new",
      "giardia_new", "e_bieneusi_new", "eaec_new"
    )
  )
  
  # One-hot encode covariates
  one_hot_data <- one_hot_encode(
    imp_data,
    laz_var_name = "final_haz",
    covariate_list = c(
      "sex", "age", "ses_quintile", "enr_haz",
      "I_followup_days",
      "I_followup_days_x_followup_days",
      "site", "edu_bin", "num_hh_lt5",
      "imp_water", "imp_sanit"
    ),
    abx_var_name = "guideline_abx",
    site_var_name = "site",
    site_interaction = TRUE,
    severity_list = c(
      "dysentery", "any_vom", "lsstools",
      "dehyd_level", "duration_pre_enroll", "any_fev"
    ),
    age_var_name = "age"
  )
  
  # Run IPD analysis
  ipd_results <- agaipw(
    data = one_hot_data$data,
    laz_var_name = "final_haz",
    abx_var_name = "guideline_abx",
    infection_var_name = "shig_attr",
    followup_var_names = c(
      "I_followup_days",
      "I_followup_days_x_followup_days"
    ),
    site_var_name = one_hot_data$site_var_names,
    covariate_list = one_hot_data$covariate_list,
    pathogen_quantity_list = c(
      "shigella_new", "rotavirus_new",
      "adenovirus_new", "etec_new",
      "cryptosporidium_new", "astrovirus_new",
      "norovirus_gii_new", "tepec_new",
      "campylobacter_new", "sapovirus_new",
      "e_bieneusi_new", "giardia_new", "eaec_new"
    ),
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
    msm_formula = "age"
  )
  
  # Remove models to reduce file size
  ipd_results$aipw_est$aipw_models <- NULL
  
  # Save results for each excluded study
  study_name <- gsub("-", "_", tolower(excluded_study))
  
  saveRDS(
    ipd_results,
    here::here(
      "results/no_etiology",
      paste0(
        "ipd_loo_", study_name,
        "_no_etiology_guideline_abx_no_models.Rds"
      )
    )
  )
  
  message("Completed: ", excluded_study)
  
  rm(loo_data, imp_data, one_hot_data, ipd_results)
  gc()
}

# ------------------------------------------------------------------
# Run leave-one-study-out analyses - case-control
# ------------------------------------------------------------------

ipd_data <- readRDS(here::here(
  "data/ipd_data/ipd_data_case_control_msd_tac_or_culture_shig_bin_abx.Rds"
))

studies <- c("MALED", "GEMS", "VIDA")

table(ipd_data$study)

for (excluded_study in studies) {
  
  message("Running case-control IPD excluding: ", excluded_study)
  
  # Exclude one study
  loo_data <- ipd_data[
    ipd_data$study != excluded_study, ,
    drop = FALSE
  ]
  
  # Impute cases and controls separately
  case_ipd_data <- loo_data[loo_data$case == 1, , drop = FALSE]
  control_ipd_data <- loo_data[loo_data$case == 0, , drop = FALSE]
  
  imp_data_case <- impute_covariates(
    data = case_ipd_data,
    imp_covariates = c(
      "sex", "age", "ses_quintile", "enr_haz",
      "edu_bin", "site", "I_followup_days",
      "I_followup_days_x_followup_days",
      "dysentery", "any_vom", "any_fev",
      "dehyd_level", "lsstools", "duration_pre_enroll",
      "shigella_new", "rotavirus_new",
      "adenovirus_new", "etec_new",
      "cryptosporidium_new", "astrovirus_new",
      "norovirus_new", "tepec_new",
      "campylobacter_new", "sapovirus_new",
      "giardia_new", "e_bieneusi_new", "eaec_new"
    )
  )
  
  imp_data_control <- impute_covariates(
    data = control_ipd_data,
    imp_covariates = c(
      "sex", "age", "ses_quintile", "enr_haz",
      "edu_bin", "site", "I_followup_days",
      "I_followup_days_x_followup_days"
    )
  )
  
  # Combine imputed cases and controls
  imp_data_case <- imp_data_case[
    , colnames(imp_data_case) %in% colnames(ipd_data),
    drop = FALSE
  ]
  
  imp_data_control <- imp_data_control[
    , colnames(imp_data_control) %in% colnames(ipd_data),
    drop = FALSE
  ]
  
  imp_data <- rbind(imp_data_case, imp_data_control)
  
  # One-hot encode covariates
  one_hot_data <- one_hot_encode(
    imp_data,
    laz_var_name = "final_haz",
    covariate_list = c(
      "sex", "age", "ses_quintile", "enr_haz",
      "I_followup_days",
      "I_followup_days_x_followup_days",
      "site", "edu_bin", "num_hh_lt5",
      "imp_water", "imp_sanit"
    ),
    abx_var_name = "guideline_abx",
    site_var_name = "site",
    site_interaction = TRUE,
    severity_list = c(
      "dysentery", "any_vom", "lsstools",
      "dehyd_level", "duration_pre_enroll", "any_fev"
    ),
    age_var_name = "age"
  )
  
  # Run case-control IPD analysis
  ipd_results_cc <- agaipw(
    data = one_hot_data$data,
    laz_var_name = "final_haz",
    abx_var_name = "guideline_abx",
    case_var_name = "case",
    followup_var_names = c(
      "I_followup_days",
      "I_followup_days_x_followup_days"
    ),
    site_var_name = one_hot_data$site_var_names,
    covariate_list = one_hot_data$covariate_list,
    pathogen_quantity_list = c(
      "shigella_new", "rotavirus_new",
      "adenovirus_new", "etec_new",
      "cryptosporidium_new", "astrovirus_new",
      "norovirus_new", "tepec_new",
      "campylobacter_new", "sapovirus_new",
      "e_bieneusi_new", "giardia_new", "eaec_new"
    ),
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
    case_control = TRUE
  )
  
  # Remove models to reduce file size
  ipd_results_cc$aipw_est$aipw_models <- NULL
  
  # Save results for each excluded study
  study_name <- tolower(excluded_study)
  
  saveRDS(
    ipd_results_cc,
    here::here(
      "results/case_control",
      paste0(
        "ipd_loo_", study_name,
        "_msd_case_control_tac_or_culture_guideline_abx_no_models.Rds"
      )
    )
  )
  
  message("Completed: ", excluded_study)
  
  rm(
    loo_data, case_ipd_data, control_ipd_data,
    imp_data_case, imp_data_control, imp_data,
    one_hot_data, ipd_results_cc
  )
  gc()
}
