# ------------------------------------------------
# Table 1 for IPD data
# ------------------------------------------------

library(tidyverse)
library(labelled)
library(gtsummary)
library(flextable)
library(officer)

here::i_am("misc/table1_ipd.R")

ipd_data <- readRDS("ipd_data/ipd_data_no_etiology.Rds")
ipd_case_control <- readRDS("ipd_data/ipd_data_case_control.Rds")

# CHECK THAT CASE == SHIG_ATTR
#table(ipd_case_control$study, ipd_case_control$case)
#table(ipd_data$study, ipd_data$shig_attr)
# N CASE == N SHIG ATTR

# Get N by study (including people diarrhea other etiology)
N_by_study <- ipd_data %>%
  group_by(study) %>%
  summarise(N = n(),
            N_kids = length(unique(first_id)))

N_by_study_cc <- ipd_case_control %>%
  group_by(study) %>%
  summarise(N = n(),
            N_kids = length(unique(first_id)))

control <- ipd_case_control[ipd_case_control$case == 0,] 

N_by_study_control <- control %>%
  group_by(study) %>%
  summarise(N = n(),
            N_kids = length(unique(first_id)))

# Total kids = N_kids in N_by_study + N_kids in N_by_study_control
sum(N_by_study$N_kids) + sum(N_by_study_control$N_kids)

# Subset case control data to controls only to merge in with no etiology data
# Remove differing columns
# Add N and control label
# Fill in NAs with 0s where applicable
control_subset <- ipd_case_control %>%
  mutate(all_abx = dplyr::case_when(
    all_abx == "Maybe effective antibiotics" ~ "Possibly effective antibiotics",
    all_abx == "WHO recommended antibiotics" ~ "Guideline recommended antibiotics",
    TRUE ~ all_abx  # retain existing values
  )) %>%
  mutate(all_abx = if_else(study == "MALED", "No or ineffective antibiotics", all_abx)) %>%
  filter(case == 0) %>%
  select(-case, -case_id, -shig_attr) %>%
  left_join(N_by_study_cc) %>%
  mutate(shig_status = "Control") %>%
  mutate(study_label = paste0(study, " (N = ", N, ")")) 
  
# Subset IPD data to only people with shigella or no etiology 
# (eliminate people who don't fall into either category)
# Add label and N
ipd_data_subset <- ipd_data %>%
  filter(!(shig_attr == 0 & no_etiology == 0)) %>%
  mutate(shig_status = if_else(shig_attr == 1, "Shigella", "No etiology")) %>%
  left_join(N_by_study, by = "study") %>%
  mutate(study_label = paste0(study, " (N = ", N, ")")) %>%
  select(-no_etiology, -culture_shig_attr, -shig_attr)

combo_ipd <- rbind(control_subset, ipd_data_subset)

combo_ipd$all_abx <- factor(combo_ipd$all_abx, levels = c("No or ineffective antibiotics",
                                                          "Possibly effective antibiotics",
                                                          "Guideline recommended antibiotics"))

# Set variables that were missing in some studies back to NA
# from IPD data prep:
# water / sanitation - missing abcd - set to 0
# num children in hh - missing maled - set to 1
# education - missing in abcd - set to 0
# fever - missing in abcd - set to 0

combo_ipd <- combo_ipd %>%
  mutate(imp_water = if_else(study == "ABCD", NA, imp_water),
          imp_sanit = if_else(study == "ABCD", NA, imp_water),
          num_hh_lt5 = if_else(study == "MALED", NA, num_hh_lt5),
          edu_bin = if_else(study == "ABCD", NA, edu_bin),
          any_fev = if_else(study == "ABCD", NA, any_fev))

# Create variables for detected, age group
combo_ipd <- combo_ipd %>%
  mutate(shigella_bin = if_else(shigella_new > 0 , 1, 0),
         rotavirus_bin = if_else(rotavirus_new > 0, 1, 0),
         adenovirus_bin = if_else(adenovirus_new > 0, 1, 0),
         etec_bin = if_else(etec_new > 0, 1, 0),
         cryptosporidium_bin = if_else(cryptosporidium_new > 0, 1, 0),
         astrovirus_bin = if_else(astrovirus_new > 0, 1, 0),
         norovirus_bin = if_else(norovirus_new > 0, 1, 0),
         tepec_bin = if_else(tepec_new > 0, 1, 0),
         campylobacter_bin = if_else(campylobacter_new > 0, 1, 0),
         sapovirus_bin = if_else(sapovirus_new > 0, 1, 0),
         giardia_bin = if_else(giardia_new > 0, 1, 0),
         e_bieneusi_bin = if_else(e_bieneusi_new > 0, 1, 0),
         eaec_bin = if_else(eaec_new > 0, 1, 0)) %>%
  set_variable_labels(shigella_bin = 'Shigella detected',
                      rotavirus_bin = "Rotavirus detected",
                      adenovirus_bin = "Adenovirus detected",
                      etec_bin = "ETEC detected",
                      cryptosporidium_bin = "Cryptosporidium detected",
                      astrovirus_bin = "Astrovirus detected",
                      norovirus_bin = "Norovirus detected",
                      tepec_bin = "tEPEC detected",
                      campylobacter_bin = "Campylobacter detected",
                      sapovirus_bin = "Sapovirus detected",
                      giardia_bin = "Giardia detected",
                      e_bieneusi_bin = "E bieneusi detected",
                      eaec_bin = "EAEC detected") %>%
  mutate(age_bin = if_else(age >= 0 & age < 12, "0-11 months",
                           if_else(age >= 12 & age <24, "12-23 months","24-59 months")),
         water_bin = if_else(imp_water < 0.5, 0, 1),
         sanit_bin = if_else(imp_sanit < 0.5, 0, 1)) %>%
  set_variable_labels(enr_haz = "Enrollment HAZ",
                      age_bin = "Age",
                      num_hh_lt5 = "Number of children in\nhousehold age <5 years",
                      edu_bin = "Primary caregiver education > primary school",
                      water_bin = "Improved water",
                      sanit_bin = "Improved sanitation",
                      sex = "Sex",
                      all_abx = "Antibiotic type",
                      dysentery = "Dysentery",
                      any_vom = "Any vomitting",
                      any_fev = "Any fever",
                      lsstools = "Number of loose stools",
                      dehyd_level = "Dehydration",
                      duration_pre_enroll = "Duration of episode before enrollment")

# Make new var for study x type
combo_ipd <- combo_ipd %>%
  mutate(study_shig = paste(study, shig_status, sep = "_"))

# Factor to ensure ordered
combo_ipd$study_shig <- factor(combo_ipd$study_shig, 
                          levels = c("GEMS_Shigella", "GEMS_No etiology", "GEMS_Control",
                                     "MALED_Shigella", "MALED_No etiology", "MALED_Control", 
                                     "VIDA_Shigella", "VIDA_No etiology", "VIDA_Control", 
                                     "ABCD_Shigella", "ABCD_No etiology",  
                                     "EFGH_Shigella", "EFGH_No etiology"))

table_1 <- tbl_summary(
  data = combo_ipd,
  by = study_shig,  # stratify by combined variable
  include = c("sex", "age_bin", "enr_haz", "num_hh_lt5", "edu_bin", "water_bin", "sanit_bin",       # COVARIATES
              "dysentery",  "any_vom", "any_fev", "lsstools", "dehyd_level", "duration_pre_enroll", # SEVERITY
              #"shigella_bin", "rotavirus_bin", "adenovirus_bin", "etec_bin", "cryptosporidium_bin", # PATHOGEN DETECTED
              #"astrovirus_bin", "norovirus_bin", "tepec_bin", "campylobacter_bin", "sapovirus_bin", 
              #"giardia_bin", "e_bieneusi_bin", "eaec_bin",
              "all_abx"),
  statistic = list(
    all_continuous() ~ "{mean} ({sd})",  # fallback
    enr_haz ~ "{mean} ({sd})",
    num_hh_lt5 ~ "{median} ({p25}, {p75})"
  ),
  missing = "ifany",
) %>%
  modify_header(label = "**Covariate**",
                c("stat_1", "stat_4", "stat_7", "stat_10", "stat_12") ~ "**Shigella**\n(N = {n})",
                c("stat_2", "stat_5", "stat_8", "stat_11", "stat_13") ~ "**No etiology**\n(N = {n})",
                c("stat_3", "stat_6", "stat_9") ~ "**Non-diarrheal control**\n(N = {n})") %>%
  modify_footnote_header("N episodes", columns = all_stat_cols()) %>%
  modify_footnote_header(
    footnote = "Controls in MALED are periods without diarrhea among 1715 children",
    columns = "stat_6",
    replace = FALSE
  ) %>%
  modify_footnote_body(
    footnote = "n (%)",
    columns = "label",
    rows = variable %in% c("sex", "age_bin","edu_bin",
                           "water_bin", "sanit_bin", "dysentery",
                           "any_vom", "any_fev","lsstools",
                           "dehyd_level", "all_abx") & row_type == "label",
    replace = FALSE
  ) %>%
  modify_footnote_body(
    footnote = "Mean (SD)",
    columns = "label",
    rows = variable %in% c("enr_haz","duration_pre_enroll") & row_type == "label",
    replace = FALSE
  ) %>%
  modify_footnote_body(
    footnote = "Median (IQR)",
    columns = "label",
    rows = variable %in% c("num_hh_lt5") & row_type == "label",
    replace = FALSE
  )  %>%
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**GEMS**",
                         c("stat_4", "stat_5", "stat_6") ~ "**MALED**",
                         c("stat_7", "stat_8", "stat_9") ~ "**VIDA**",
                         c("stat_10", "stat_11") ~ "**ABCD**",
                         c("stat_12", "stat_13") ~ "**EFGH**") 

# Flex table output to word
table_1_flex <- table_1 %>%
  as_flex_table() %>%
  fontsize(size = 7, part = "all") %>%
  width(width = 0.5) 

# Create and save Word document
read_docx() %>%
  body_add_flextable(table_1_flex) %>%
  body_end_section_landscape() %>%  
  print(target = here::here("results/final/figures/table1_summary.docx"))

# EXPORT AS WORD AND MANUALLY GET RID OF MISSING ROW WHERE APPLICABLE
############################################################################

# Supplementary table for pathogens
table_supp <- tbl_summary(
  data = combo_ipd,
  by = study_shig,  # stratify by combined variable
  include = c("shigella_bin", "rotavirus_bin", "adenovirus_bin", "etec_bin", "cryptosporidium_bin", # PATHOGEN DETECTED
              "astrovirus_bin", "norovirus_bin", "tepec_bin", "campylobacter_bin", "sapovirus_bin", 
              "giardia_bin", "e_bieneusi_bin", "eaec_bin"),
  # statistic = list(
  #   all_continuous() ~ "{mean} ({sd})",  # fallback
  #   enr_haz ~ "{mean} ({sd})",
  #   num_hh_lt5 ~ "{median} ({p25}, {p75})"
  # ),
  missing = "ifany",missing_text = "Missing",
) %>%
  modify_header(label = "**Covariate**",
                c("stat_1", "stat_4", "stat_7", "stat_10", "stat_12") ~ "**Shigella**\n(N = {n})",
                c("stat_2", "stat_5", "stat_8", "stat_11", "stat_13") ~ "**No etiology**\n(N = {n})",
                c("stat_3", "stat_6", "stat_9") ~ "**Non-diarrheal control**\n(N = {n})") %>%
  #modify_footnote_header("N episodes", columns = all_stat_cols(), replace = FALSE) %>%
  modify_footnote_header(
    footnote = "Controls in MALED are periods without diarrhea among 1715 children",
    columns = "stat_6",
    replace = FALSE
  ) %>%
  modify_footnote_body(
    footnote = "n (%)",
    columns = "label",
    rows = variable %in% c("sex", "age_bin","edu_bin",
                           "water_bin", "sanit_bin", "dysentery",
                           "any_vom", "any_fev","lsstools",
                           "dehyd_level", "all_abx") & row_type == "label",
    replace = FALSE
  ) %>%
  modify_footnote_body(
    footnote = "Mean (SD)",
    columns = "label",
    rows = variable %in% c("enr_haz","duration_pre_enroll") & row_type == "label",
    replace = FALSE
  ) %>%
  modify_footnote_body(
    footnote = "Median (IQR)",
    columns = "label",
    rows = variable %in% c("num_hh_lt5") & row_type == "label",
    replace = FALSE
  )  %>%
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**GEMS**",
                         c("stat_4", "stat_5", "stat_6") ~ "**MALED**",
                         c("stat_7", "stat_8", "stat_9") ~ "**VIDA**",
                         c("stat_10", "stat_11") ~ "**ABCD**",
                         c("stat_12", "stat_13") ~ "**EFGH**") 

# Flex table output to word
table_supp_flex <- table_supp %>%
  as_flex_table() %>%
  fontsize(size = 6.5, part = "all") %>%
  width(width = 0.5) 

# Create and save Word document
read_docx() %>%
  body_add_flextable(table_supp_flex) %>%
  body_end_section_landscape() %>%  
  print(target = here::here("results/final/figures/table_supp_summary.docx"))



############################################################################
# table1_noetiology_sev <- table1(~ dysentery + any_vom + any_fev + lsstools + 
#                                   dehyd_level + duration_pre_enroll | study_label*shig_status,
#                                 data = ipd_data_subset, 
#                                 overall = FALSE,
#                                 rowlabelhead = c("Severity"))
#                                 
# table1_noetiology_pathq <- table1(
#   ~ shigella_bin + rotavirus_bin + adenovirus_bin + etec_bin + cryptosporidium_bin + 
#     astrovirus_bin + norovirus_bin + tepec_bin + campylobacter_bin + sapovirus_bin + 
#     giardia_bin + e_bieneusi_bin + eaec_bin | study_label*shig_status, 
#   data = ipd_data_subset, 
#   overall = FALSE,
#   rowlabelhead = "Number of episodes with pathogen detected", 
#   render.continuous = function(x, ...) {
#     n <- sum(!is.na(x))
#     n1 <- sum(x == 1, na.rm = TRUE)
#     pct <- 100 * n1 / n
#     sprintf("%d (%.1f%%)", n1, pct)
#   }
# )                     
#   
# # ---------------------------------------------------------------------
# # Case control
# # ---------------------------------------------------------------------
# 
# ipd_cc_data <- readRDS("results/aipw/ipd_data_case_control.Rds")
# 
# # Get N by study (including people diarrhea other etiology)
# N_by_study_cc <- ipd_cc_data %>%
#   group_by(study) %>%
#   summarise(N = n())
# 
# # Only people with shigella or no etiology
# ipd_cc_data <- ipd_cc_data %>%
#   left_join(N_by_study_cc, by = "study") %>%
#   mutate(study_label = paste0(study, " (N = ", N, ")")) %>%
#   mutate(case_status = if_else(case == 1, "Case", "Control")) %>%
#   mutate(shigella_bin = if_else(shigella_new > 0 , 1, 0),
#          rotavirus_bin = if_else(rotavirus_new > 0, 1, 0),
#          adenovirus_bin = if_else(adenovirus_new > 0, 1, 0),
#          etec_bin = if_else(etec_new > 0, 1, 0),
#          cryptosporidium_bin = if_else(cryptosporidium_new > 0, 1, 0),
#          astrovirus_bin = if_else(astrovirus_new > 0, 1, 0),
#          norovirus_bin = if_else(norovirus_new > 0, 1, 0),
#          tepec_bin = if_else(tepec_new > 0, 1, 0),
#          campylobacter_bin = if_else(campylobacter_new > 0, 1, 0),
#          sapovirus_bin = if_else(sapovirus_new > 0, 1, 0),
#          giardia_bin = if_else(giardia_new > 0, 1, 0),
#          e_bieneusi_bin = if_else(e_bieneusi_new > 0, 1, 0),
#          eaec_bin = if_else(eaec_new > 0, 1, 0)) %>%
#   set_variable_labels(shigella_bin = 'Shigella detected',
#                       rotavirus_bin = "Rotavirus detected",
#                       adenovirus_bin = "Adenovirus detected",
#                       etec_bin = "ETEC detected",
#                       cryptosporidium_bin = "Cryptosporidium detected",
#                       astrovirus_bin = "Astrovirus detected",
#                       norovirus_bin = "Norovirus detected",
#                       tepec_bin = "tEPEC detected",
#                       campylobacter_bin = "Campylobacter detected",
#                       sapovirus_bin = "Sapovirus detected",
#                       giardia_bin = "Giardia detected",
#                       e_bieneusi_bin = "E bieneusi detected",
#                       eaec_bin = "EAEC detected")
# 
# table1_cc_cov <- table1(~ sex + age + ses_quintile + enr_haz + num_hh_lt5 +
#                                   edu_bin + imp_water + imp_sanit | study_label*case_status, 
#                                 data = ipd_cc_data, 
#                                 overall = FALSE,
#                                 rowlabelhead = c("Covariates"))
# 
# table1_cc_sev <- table1(~ dysentery + any_vom + any_fev + lsstools + 
#                                   dehyd_level + duration_pre_enroll | study_label*case_status,
#                                 data = ipd_cc_data, 
#                                 overall = FALSE,
#                                 rowlabelhead = c("Severity"))
# 
# table1_cc_pathq <- table1(
#   ~ shigella_bin + rotavirus_bin + adenovirus_bin + etec_bin + cryptosporidium_bin + 
#     astrovirus_bin + norovirus_bin + tepec_bin + campylobacter_bin + sapovirus_bin + 
#     giardia_bin + e_bieneusi_bin + eaec_bin | study_label*case_status, 
#   data = ipd_cc_data, 
#   overall = FALSE,
#   rowlabelhead = "Number of episodes with pathogen detected", 
#   render.continuous = function(x, ...) {
#     n <- sum(!is.na(x))
#     n1 <- sum(x == 1, na.rm = TRUE)
#     pct <- 100 * n1 / n
#     sprintf("%d (%.1f%%)", n1, pct)
#   }
# )        
