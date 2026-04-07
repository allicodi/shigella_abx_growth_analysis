# -------------------------------------------------------------------------
# Script to run analysis stratified by severity for case-control obs abx
# -------------------------------------------------------------------------

here::i_am("code/run_analyses/run_severity_stratified.R")

devtools::load_all("~/Documents/shigella_projects/packages/abxGrowth")

library(SuperLearner)

source(here::here("code/run_analyses/SL.wrappers.R"))

# msd only, lsd only, if both false include all (no subset)
msd <- FALSE
lsd <- FALSE

# change option for age stratified
age_stratified <- TRUE

age_0_11 <- FALSE
age_12_23 <- FALSE
age_24_59 <- TRUE

if(msd){
  data <- readRDS("data/ipd_data/ipd_data_case_control_msd_tac_or_culture_shig.Rds")
} else if(lsd){
  data <- readRDS("data/ipd_data/ipd_data_case_control_lsd_tac_or_culture_shig.Rds")
} else {
  data <- readRDS("data/ipd_data/ipd_data_case_control_tac_or_culture_shig.Rds")
}

if(age_stratified){
  if(age_0_11){
    # get cases in this age range, control is allowed to be older if matched 
    caseids_0_11 <- unique(data$case_id[which(data$age >= 0 & data$age < 12 & data$case == 1)])
    data <- data[which(data$case_id %in% caseids_0_11),]
    
    #data <- data[data$age >= 0 & data$age < 12,]
  } else if(age_12_23){
    
    caseids_12_23 <- unique(data$case_id[which(data$age >= 12 & data$age < 24 & data$case == 1) ])
    data <- data[which(data$case_id %in% caseids_12_23),]
    
    #data <- data[data$age >= 12 & data$age < 24,]
  } else{
    
    caseids_24_60 <- unique(data$case_id[which(data$age >= 24 & data$age < 60 & data$case == 1)])
    data <- data[which(data$case_id %in% caseids_24_60),]
    
    #data <- data[data$age >= 24 & data$age < 60,]
  }
}

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

case_data <- data[data$case == 1,]
control_data <- data[data$case == 0,]

# Impute cases and controls sep
imp_data_case <- impute_covariates(data = case_data,
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

imp_data_control <- impute_covariates(data = control_data,
                                      imp_covariates = c("sex",
                                                         "age",
                                                         "ses_quintile",
                                                         "enr_haz",
                                                         "edu_bin",
                                                         "site",
                                                         "I_followup_days",
                                                         "I_followup_days_x_followup_days"))


imp_data_case <- imp_data_case[,colnames(imp_data_case) %in% colnames(data), drop = FALSE]
imp_data_control <- imp_data_control[,colnames(imp_data_control) %in% colnames(data), drop = FALSE]
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

# SET PARAMETERS
data = one_hot_data$data
laz_var_name = "final_haz"
abx_var_name = "all_abx"
case_var_name = "case"
followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days")
site_var_name = one_hot_data$site_var_names
covariate_list = one_hot_data$covariate_list
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
                           "eaec_new")
severity_list = one_hot_data$severity_list
first_id_var_name = "first_id"
outcome_type = "gaussian"
sl.library.outcome.case = sl.library.with.abx
sl.library.outcome.control = sl.library.without.abx
sl.library.treatment = sl.library.without.abx
sl.library.infection = sl.library.without.abx
sl.library.missingness.case = sl.library.with.abx
sl.library.missingness.control = sl.library.without.abx
v_folds = 5
return_models = TRUE
msm = TRUE
msm_var_name = "age"
msm_formula = "age"
case_control = TRUE
seed = 12345

# ------------------------------------------------------------
# STEP 0: Create subsets of data for model fitting
# ------------------------------------------------------------

set.seed(seed)

# Rename abx levels if spaces in them
abx_levels <- levels(factor(data[[abx_var_name]]))
abx_levels_new <- ifelse(is.na(abx_levels), NA, gsub("[ /]", "_", abx_levels))
data[[abx_var_name]] <- factor(data[[abx_var_name]], levels = abx_levels, labels = abx_levels_new)

# get case vs control data
case_data <- data[data[[case_var_name]] == 1,]
control_data <- data[data[[case_var_name]] == 0,]

case_data_idx <- which(data[[case_var_name]] == 1)
control_data_idx <- which(data[[case_var_name]] == 0)

# Case data prep
I_Y_case <- ifelse(is.na(case_data[[laz_var_name]]), 1, 0)
Y_case <- case_data[[laz_var_name]]
attributes(Y_case) <- NULL
covariates_case <- case_data[, covariate_list, drop = FALSE]
severity_case <- case_data[, severity_list, drop = FALSE]
pathogen_q_case <- case_data[, pathogen_quantity_list, drop = FALSE]
abx_case <- case_data[, abx_var_name, drop = FALSE]

# Complete case data (excluding missing outcome)
case_data_complete <- case_data[!is.na(case_data[[laz_var_name]]), ]

Y_case_complete <- case_data_complete[[laz_var_name]]
attributes(Y_case_complete) <- NULL
covariates_case_complete <- case_data_complete[, covariate_list, drop = FALSE]
severity_case_complete <- case_data_complete[, severity_list, drop = FALSE]
pathogen_q_case_complete <- case_data_complete[, pathogen_quantity_list, drop = FALSE]
abx_case_complete <- case_data_complete[, abx_var_name, drop = FALSE]

# Control data prep
I_Y_control <- ifelse(is.na(control_data[[laz_var_name]]), 1, 0)
Y_control <- control_data[[laz_var_name]]
attributes(Y_control) <- NULL
if(is.null(covariate_list)){
  covariate_list <- covariate_list
}
covariates_control <- control_data[, covariate_list, drop = FALSE]

# Complete control data (excluding missing outcome)
control_data_complete <- control_data[!is.na(control_data[[laz_var_name]]), ]

Y_control_complete <- control_data_complete[[laz_var_name]]
attributes(Y_control_complete) <- NULL
covariates_control_complete <- control_data_complete[, covariate_list, drop = FALSE]

# ------------------------------------------------------------
# STEP 1: Fit & predict from outcome models
# ------------------------------------------------------------

## Model 1a: Outcome model in cases
outcome_model_1a <- SuperLearner::SuperLearner(Y = Y_case_complete, 
                                               X = data.frame(abx_case_complete,
                                                              covariates_case_complete, 
                                                              severity_case_complete, 
                                                              pathogen_q_case_complete),
                                               family = outcome_type, 
                                               SL.library = sl.library.outcome.case,
                                               cvControl = list(V = v_folds))

## Model 1b: Outcome model in controls
outcome_model_1b <- SuperLearner::SuperLearner(Y = Y_control_complete, 
                                               X = data.frame(covariates_control_complete),
                                               family = outcome_type, 
                                               SL.library = sl.library.outcome.control,
                                               cvControl = list(V = v_folds))

## Predict outcomes

# For each level of abx, predict setting abx = x, infection = 1 & abx = x, infection = 0
abx_levels <- unique(case_data[[abx_var_name]])[!is.na(unique(case_data[[abx_var_name]]))]

outcome_vectors_1a <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
outcome_vectors_1b <- data.frame(matrix(ncol = 1, nrow = nrow(data)))

colnames(outcome_vectors_1a) <- paste0("abx_observed")
colnames(outcome_vectors_1b) <- paste0("control")

outcome_vectors_1a[case_data_idx, 1] <- stats::predict(outcome_model_1a, newdata = case_data[, c(abx_var_name,
                                                                                                 covariate_list,
                                                                                                 severity_list,
                                                                                                 pathogen_quantity_list)], type = "response")$pred

# Replace NA controls with 0
outcome_vectors_1a[is.na(outcome_vectors_1a)] <- 0

outcome_vectors_1b[, 1] <- stats::predict(outcome_model_1b, newdata = data[,covariate_list], type = "response")$pred

# ------------------------------------------------------------
# STEP 2: Fit & predict from propensity models
# ------------------------------------------------------------

# If followup_days related variables, remove from missingness models
if(!any(is.na(followup_var_names))){
  covariate_list_no_followup <- covariate_list[!(covariate_list %in% followup_var_names)]
} else{
  covariate_list_no_followup <- covariate_list
}

prop_vectors_2a <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
prop_vectors_3a <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
prop_vectors_3b <- data.frame(matrix(ncol = 1, nrow = nrow(data)))

colnames(prop_vectors_2a) <- c("case")
colnames(prop_vectors_3a) <- paste0("abx_observed")
colnames(prop_vectors_3b) <- c("control")

## Part 2: Propensity models for shigella (or other infection) attribution

# 2a_1 = Shigella Attributable ~ BL Cov
prop_model_2a <- SuperLearner::SuperLearner(Y = data[[case_var_name]],
                                            X = data[, covariate_list, drop = FALSE],
                                            family = stats::binomial(),
                                            SL.library = sl.library.infection, 
                                            cvControl = list(V = v_folds))
tmp_pred_2a <- prop_model_2a$SL.pred
prop_vectors_2a[,1] <- tmp_pred_2a

## Part 3: Propensity models for missingness

covariates_case_no_site <- case_data[, covariate_list_no_followup, drop = FALSE]
covariates_control_no_site <- control_data[, covariate_list_no_followup, drop = FALSE]

## Missingness model in cases
prop_model_3a <- SuperLearner::SuperLearner(Y = I_Y_case,
                                            X = data.frame(abx_case,
                                                           covariates_case_no_site,
                                                           severity_case,
                                                           pathogen_q_case),
                                            family = stats::binomial(),
                                            SL.library = sl.library.missingness.case,
                                            cvControl = list(V = v_folds))

## Missingness model in controls
prop_model_3b <- SuperLearner::SuperLearner(Y = I_Y_control,
                                            X = data.frame(covariates_control_no_site),
                                            family = stats::binomial(),
                                            SL.library = sl.library.missingness.control,
                                            cvControl = list(V = v_folds))

# Predict setting each abx level

prop_vectors_3a[case_data_idx,1] <- prop_model_3a$SL.pred
# Fill in NAs with 0
prop_vectors_3a[is.na(prop_vectors_3a)] <- 0

prop_vectors_3b[control_data_idx,1] <- prop_model_3b$SL.pred
prop_vectors_3b[is.na(prop_vectors_3b)] <- 0

# -------------------------------------------------
# STEP 3: AIPW estimates and confidence intervals
# -------------------------------------------------

## Plug-in estimates
plug_ins_case <- mean(outcome_vectors_1a[case_data_idx,])
plug_ins_control <- mean(outcome_vectors_1b[case_data_idx,])

# 'individual level effect' save to do MSM stuff later
outcome_vectors_difference <- outcome_vectors_1a - outcome_vectors_1b

## EIF for bias corrections
case_eifs <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
control_eifs <- data.frame(matrix(ncol = 1, nrow = nrow(data)))

colnames(case_eifs) <- paste0("case_eif_observed")
colnames(control_eifs) <- paste0("control_eif_observed")

# Truncate any large propensity scores
ps_trunc_level <- 0.01
if(!is.na(ps_trunc_level)){
  #  any observation that is < ps_trunc_level or > 1-ps_trunc_level should be changed to ps_trunc_level or 1 - ps_trunc_level
  
  truncate_ps <- function(mat, ps_trunc_level){
    mat[mat < ps_trunc_level] <- ps_trunc_level
    mat[mat > 1- ps_trunc_level] <- 1- ps_trunc_level
    return(mat)
  }
  
  prop_vectors_2a <- truncate_ps(prop_vectors_2a, ps_trunc_level)
  prop_vectors_3a <- truncate_ps(prop_vectors_3a, ps_trunc_level)
  prop_vectors_3b <- truncate_ps(prop_vectors_3b, ps_trunc_level)
}

# 1 - Bias correction for case, abx level = a

I_Case <- data[[case_var_name]]
P_Case <- mean(prop_vectors_2a[,1])

I_Delta_0 <- as.numeric(!is.na(data[[laz_var_name]])) # Indicator NOT missing
#P_Delta_0__Case_all <- 1 - prop_vectors_3a[,i]
P_Delta_0__Case_all <- 1 - prop_vectors_3a[,1]

obs_outcome <- ifelse(is.na(data[[laz_var_name]]), 0, data[[laz_var_name]])  
Qbar_Case_Abx_a_Covariates <- outcome_vectors_1a[,1]

eif_case_vec <- (I_Case / P_Case) * (I_Delta_0 / P_Delta_0__Case_all) * (obs_outcome - Qbar_Case_Abx_a_Covariates) +
  (I_Case / P_Case) * (Qbar_Case_Abx_a_Covariates - plug_ins_case[1])

# correct for / 0 with P_Abx_a__Case_Covariates, should be 0'd out
eif_case_vec <- ifelse(is.nan(eif_case_vec), 0, eif_case_vec)

case_eifs[,1] <- eif_case_vec


# 2 - Bias correction for control
I_control <- ifelse(data[[case_var_name]] == 1, 0, 1)
P_control__Covariates <- 1 - prop_vectors_2a[,1]

I_Delta_0 <- as.numeric(!is.na(data[[laz_var_name]])) # Indicator NOT missing
I_Delta_0__Control_all <- 1 - prop_vectors_3b[,1]

P_Case__all <- prop_vectors_2a[,1]

Qbar_Control_Covariates <- outcome_vectors_1b[,1]

eif_control_vec <- (I_control / P_control__Covariates) * (I_Delta_0 / I_Delta_0__Control_all) * (P_Case__all / P_Case) * (obs_outcome - Qbar_Control_Covariates) +
  (I_Case / P_Case) * (Qbar_Control_Covariates - plug_ins_control)

control_eif <- eif_control_vec

# Get AIPWs
aipw_case <- plug_ins_case + colMeans(case_eifs)
aipw_control <- plug_ins_control + colMeans(control_eif)  

eif_matrix <- cbind(case_eifs, control_eif)

# Get id for each participant and recreate EIFs based on this if present
if(!is.null(first_id_var_name)){
  first_id_eif_matrix <- cbind(data.frame(first_id = data[[first_id_var_name]]), eif_matrix)
  first_id_eif_matrix <- aggregate(. ~ first_id, data = first_id_eif_matrix, FUN = sum)
  scaled_matrix <- first_id_eif_matrix[,-c(1)] * (nrow(first_id_eif_matrix) / nrow(eif_matrix))
  
}else{
  scaled_matrix <- eif_matrix
}

aipws_effect <- vector("numeric", length = 1)
eifs_effect <- data.frame(matrix(ncol = 1, nrow = nrow(scaled_matrix)))

names(aipws_effect) <- paste0("effect_observed")
colnames(eifs_effect) <- paste0("effect_observed")

aipw_effect <- aipw_case[1] - aipw_control[1]
aipws_effect[1] <- aipw_effect

idx_1 <- 1
idx_2 <- 1 + 1

gradient <- rep(0, 2)

gradient[idx_1] <- 1
gradient[idx_2] <- -1

gradient <- matrix(gradient, ncol = 1)
eif_effect <- as.numeric(as.matrix(scaled_matrix) %*% gradient)
eifs_effect[,1] <- eif_effect

results_df <- data.frame(abx_levels = "observed",
                         abx_level_case = aipw_case,
                         abx_level_control = rep(aipw_control, 1),
                         effect_inf_abx_level = aipws_effect)

eif_matrix_scaled <- cbind(scaled_matrix, eifs_effect)
cov_matrix <- stats::cov(eif_matrix_scaled)
eif_hat <- sqrt( diag(cov_matrix) / nrow(eif_matrix_scaled) )

results_object <- list(results_df = results_df,
                       plug_ins_case = plug_ins_case,
                       plug_ins_control = plug_ins_control,
                       eif_matrix = eif_matrix_scaled,
                       se = eif_hat)

class(results_object) <- "aipw_case_control_observed"


aipw_models <- list(outcome_model_1a = outcome_model_1a,
                    outcome_model_1b = outcome_model_1b,
                    prop_model_2a = prop_model_2a,
                    prop_model_3a = prop_model_3a,
                    prop_model_3b = prop_model_3b)

results <- list(results_object = results_object,
                aipw_models = aipw_models)

class(results) <- "agaipw_res"

results$aipw_models <- NULL

if(!age_stratified){
  if(msd){
    saveRDS(results, here::here("results/case_control/ipd_no_abx_msd.Rds"))
  } else if(lsd){
    saveRDS(results, here::here("results/case_control/ipd_no_abx_lsd.Rds"))
  } else{
    saveRDS(results, here::here("results/case_control/ipd_no_abx.Rds"))
  }
} else{
  if(age_0_11){
    saveRDS(results, here::here("results/case_control/ipd_no_abx_0_11.Rds"))
  } else if (age_12_23){
    saveRDS(results, here::here("results/case_control/ipd_no_abx_12_23.Rds"))
  } else{
    saveRDS(results, here::here("results/case_control/ipd_no_abx_24_59.Rds"))
  }
}

