# ---------------------------------------------------------------
# Script to prep IPD data given individual prepped datasets
# ---------------------------------------------------------------

here::i_am("code/prep_data/prep_ipd.R")

library(tidyverse)
library(labelled)

#####################################################################
#                           CASE ONLY                               #
#####################################################################

full <- FALSE # FULL = TRUE all diarrhea cases (co-etiology meta-analysis); FALSE shigella tac or culture cases (original meta-analysis)
MALED_MSD <- FALSE

if(full){
  # Read in clean datasets from other studies
  abcd_data <- readRDS(here::here("data/abcd_data/abcd_data_full.Rds"))
  efgh_data <- readRDS(here::here("data/efgh_data/efgh_data_full.Rds"))
  gems_data <- readRDS(here::here("data/gems_data/gems_data_full.Rds"))
  vida_data <- readRDS(here::here("data/vida_data/vida_data_full.Rds"))
  maled_data <- readRDS(here::here("data/maled_data/maled_data_full.Rds"))
} else{
  # Read in clean datasets from other studies
  abcd_data <- readRDS(here::here("data/abcd_data/abcd_data.Rds"))
  efgh_data <- readRDS(here::here("data/efgh_data/efgh_data.Rds"))
  gems_data <- readRDS(here::here("data/gems_data/gems_data.Rds"))
  vida_data <- readRDS(here::here("data/vida_data/vida_data.Rds"))
  if(MALED_MSD){
    maled_data <- readRDS(here::here("data/maled_data/maled_data_shig_MSD_only.Rds"))
  } else{
    maled_data <- readRDS(here::here("data/maled_data/maled_data.Rds"))
  }
}


# Add column with study name to each dataset
abcd_data$study <- "ABCD"
efgh_data$study <- "EFGH"
gems_data$study <- "GEMS"
vida_data$study <- "VIDA"
maled_data$study <- "MALED"

# ------------------------------------------------------------------------------
#                                 COVARIATES
# ------------------------------------------------------------------------------

################################################################################
# 1. Sex (sex) - Male = 0, Female = 1
# 2. Age in months (age)
# 3. SES quintile (ses_quintile)
# 4. Baseline HAZ (enr_haz)
# 5. Baseline WAZ (enr_waz)
# 5. Days between episode and follow-up (followup_days) -- I_followup_days & I_followup_days_x_followup_days
# 6. Site (site) - but make it a factor for each study site individually?

# Missing in some (fill in with 0)
# 7. Water + sanitation score (missing in ABCD) (imp_water, imp_sanit)
# 8. Number of children in household < 5 (missing in MAL-ED) (num_hh_lt5)
# 9. Primary caregiver education (missing in ABCD) (edu_bin) 
################################################################################

# Rename & Reformat

# ABCD
abcd_data <- abcd_data %>%
  mutate(imp_water = 0,
         imp_sanit = 0,
         edu_bin = 0) %>%
  rename('sex' = dy1_ant_sex,
         'age' = agemchild,
         'ses_quintile' = an_ses_quintile,
         'enr_haz' = lfazscore,
         'enr_waz' = wfazscore,
         'I_followup_days' = I_an_d90_timing,
         'I_followup_days_x_followup_days' = I_an_d90_timing_x_an_d90_timing,
         'num_hh_lt5' = an_tothhlt5) 

abcd_data$sex <- ifelse(abcd_data$sex == "Female", 1, 0)
levels(abcd_data$site) <- c("Bangladesh_abcd", 
                            "Kenya_abcd", 
                            "Malawi_abcd", 
                            "Mali_abcd", 
                            "India_abcd", 
                            "Tanzania_abcd", 
                            "Pakistan_abcd")

# EFGH
efgh_data <- efgh_data %>%
  rename('age' = enr_age_months,
         'ses_quintile' = final_quintile_site,
         'site' = enroll_site,
         'I_followup_days' = I_mo3_days,
         'I_followup_days_x_followup_days' = I_mo3_days_x_mo3_days,
         'imp_sanit' = imp_toi,
         'num_hh_lt5' = enroll_ai_num_child,
         'edu_bin' = moth_ed_bin)

efgh_data$sex <- ifelse(efgh_data$sex == "Female", 1, 0)
levels(efgh_data$site) <- c("Bangladesh_efgh",
                            "Kenya_efgh",
                            "Malawi_efgh",
                            "Mali_efgh",
                            "Pakistan_efgh",
                            "Peru_efgh",
                            "The_Gambia_efgh")

# GEMS 
# naming conventions all match
gems_data$sex <- ifelse(gems_data$sex == "female", 1, 0)
levels(gems_data$site) <- c("The_Gambia_gems",
                            "Mali_gems",
                            "Mozambique_gems",
                            "Kenya_gems",
                            "India_gems",
                            "Bangladesh_gems",
                            "Pakistan_gems")

gems_data$edu_bin <- gems_data$prim_caregiver_edu_bin

gems_data$imp_water <- ifelse(gems_data$safe_water %in% c("Safely managed", "Basic"), 1, 0)
gems_data$imp_sanit <- ifelse(gems_data$safe_sanit %in% c("Safely manaed and basic"), 1, 0)

# VIDA
vida_data <- vida_data %>%
  rename('age' = agemchild,
         'edu_bin' = education_bin)

vida_data$sex <- ifelse(vida_data$sex == "Female", 1, 0)
levels(vida_data$site) <- c("The_Gambia_vida",
                            "Mali_vida",
                            "Kenya_vida")

vida_data$imp_water <- ifelse(vida_data$safe_water %in% c("Safely managed", "Basic"), 1, 0)
vida_data$imp_sanit <- ifelse(vida_data$safe_sanit %in% c("Safely managed and basic"), 1, 0)

# MAL-ED
maled_data <- maled_data %>%
  rename('age' = agemonths,
         'ses_quintile' = wami_quintile,
         'enr_haz' = baseline_haz,
         'enr_waz' = baseline_waz,
         'imp_water' = drinkimp,
         'imp_sanit' = sanitimp,
         'edu_bin' = mated_bin)

maled_data$sex <- ifelse(maled_data$sex == "female", 1, 0)
levels(maled_data$site) <- c("Bangladesh_maled",
                             "Brazil_maled",
                             "India_maled",
                             "Nepal_maled",
                             "Peru_maled",
                             "SouthAfrica_maled",
                             "Tanzania_maled")

maled_data$num_hh_lt5 <- 1

# ------------------------------------------------------------------------------
#                                 SEVERITY
# ------------------------------------------------------------------------------

################################################################################
# 1. Dysentery (0 = no, 1 = yes) (dysentery)
# 2. Vomiting during episode (0 = no, 1 = yes) (any_vom)
# 3. Dehydration level (Factor with levels No dehydration, Some dehydration, Severe dehydration) (dehyd_level)
# 4. Max number of loose stools (categorical) (Factor with levels 6 or less, 7 to 10, Over 10) (except vida is really 3 to 5, 6 to 10, over 10) (lsstools)
# 5. Duration of diarrhea pre-abx (duration_pre_enroll)

# 6. Fever during episode (0 = no, 1 = yes) (any_fev)
# 7. GEMS-like MSD (0 = no, 1, = yes) (MSD)
################################################################################

# ABCD
abcd_data <- abcd_data %>%
  mutate(dysentery = 0,
         any_vom = if_else(dy1_scrn_vomitall == "Yes", 1, 0),
         any_fev = 0,
         lsstools = if_else(dy1_scrn_lstools <=6, "6 or less per day",
                            if_else(dy1_scrn_lstools >= 7 & dy1_scrn_lstools <= 10, "7 to 10 per day", "Over 10 per day")),
         duration_pre_enroll = dy1_scrn_diardays + 1) %>%
  rename('dehyd_level' = dy1_scrn_dehydr, 
         'MSD' = gems_msd)

# EFGH
efgh_data <- efgh_data %>%
  mutate(any_vom = if_else(enroll_diar_vom_num > 0, 1, 0),
         lsstools = if_else(enroll_diar_loose_num <= 6, "6 or less per day",
                            if_else(enroll_diar_loose_num >= 7 & enroll_diar_loose_num <= 10, "7 to 10 per day", "Over 10 per day")),
         MSD = if_else(gems_msd == "Less-severe", 0, 1)) %>%
  rename('dysentery' = enroll_diar_blood,
         'any_fev' = enroll_diar_fever,
         'dehyd_level' = enroll_cond_dehyd)

# GEMS
gems_data <- gems_data %>%
  mutate(MSD = 1) %>%
  rename('any_vom' = vomit,
         'any_fev' = fever,
         'dehyd_level' = who_dehyd)

# VIDA
vida_data <- vida_data %>%
  mutate(any_vom = if_else(vom_days > 0, 1, 0),
         dehyd_level = if_else(dehydr == "None", "No dehydration", dehydr),
         MSD = 1) %>%
  rename('any_fev' = fever)

vida_data$lsstools <- ifelse(vida_data$lsstools == "3 stools" | vida_data$lsstools == "4 to 5 stools", "6 or less per day",
                             ifelse(vida_data$lsstools == "6 to 10 stools", "7 to 10 per day", "Over 10 per day"))

# MAL-ED
maled_data <- maled_data %>%
  mutate(any_vom = if_else(daysvomit > 0, 1, 0),
         dehyd_level = if_else(dehyd == "None", "No dehydration", dehyd),
         lsstools = if_else(lsstools <=6, "6 or less per day",
                            if_else(lsstools >= 7 & lsstools <= 10, "7 to 10 per day", "Over 10 per day"))) %>%
  rename('any_fev' = fever,
         'duration_pre_enroll' = duration_pre_abx)

# ------------------------------------------------------------------------------
#                             PATHOGEN QUANTITY
# ------------------------------------------------------------------------------

################################################################################
# Final list / naming conventions:

# "shigella_new"
# "rotavirus_new"
# "adenovirus_new"
# "st_etec_new"
# "etec_new"
# "cryptosporidium_new" 
# "astrovirus_new"
# "norovirus_new"
# "tepec_new"
# "campylobacter_new"
# "sapovirus_new"
# "giardia_new"
# "e_bieneusi_new"
# "eaec_new"

# 'salmonella_new'
# 'c_jejuni_coli_new'

# 9/19/25 add attributable (for WHO presentation figures)

# rotavirus_tac_attr
# norovirus_tac_attr
# adenovirus_tac_attr
# sapovirus_tac_attr
# astrovirus_tac_attr
# etec_tac_attr
# tepec_tac_attr
# cryptosporidium_tac_attr

################################################################################

# ABCD using conventions above

abcd_data <- abcd_data %>%
  rename('shigella_tac_attr' = shigella_likely,
         'rotavirus_tac_attr' = rotavirus_likely,
         'norovirus_gii_tac_attr' = norovirus_gii_likely,
         'adenovirus_tac_attr' = adenovirus_likely,
         'sapovirus_tac_attr' = sapovirus_likely,
         'astrovirus_tac_attr' = astrovirus_likely,
         'st_etec_tac_attr' = st_etec_likely, # look to see if we have lt etec in raw data
         'tepec_tac_attr' = tepec_likely,
         'cryptosporidium_tac_attr' = cryptosporidium_likely,
         'v_cholerae_tac_attr' = v_cholerae_likely,
         'salmonella_tac_attr' = salmonella_likely,
         'c_jejuni_coli_tac_attr' = c_jejuni_likely,
         'c_jejuni_coli_new' = c_jejuni_new) %>%
  rowwise() %>%
  mutate(
    co_etiology = if_else(any(c_across(ends_with("_tac_attr")) == 1), 1, 0)
  ) %>%
  ungroup()

# EFGH
# is it okay to set campylobacter_new = c_jejuni_new?

efgh_data <- efgh_data %>%
  rename(# Pathogen quantities
    'etec_new' = ETEC_new,
    'tepec_new' = tEPEC_new,
    'eaec_new' = EAEC_new,
    'campylobacter_new' = c_jejuni_new,
    # Pathogen attributable
    'shigella_tac_attr' = tac_shigella_attributable,
    'rotavirus_tac_attr' = rotavirus_attributable,
    'norovirus_gii_tac_attr' = norovirus_gii_attributable,
    'adenovirus_tac_attr' = adenovirus_40_41_attributable,
    'sapovirus_tac_attr' = sapovirus_attributable,
    'astrovirus_tac_attr' = astrovirus_attributable,
    'st_etec_tac_attr' = ST.ETEC_attributable, # look to see if we have lt etec in raw data
    'tepec_tac_attr' = tEPEC_attributable,
    'cryptosporidium_tac_attr' = cryptosporidium_attributable,
    'v_cholerae_tac_attr' = v_cholerae_attributable,
    'salmonella_tac_attr' = salmonella_attributable,
    'c_jejuni_coli_tac_attr' = c_jejuni_coli_attributable) %>%
  rowwise() %>%
  mutate(
    c_jejuni_coli_new = campylobacter_new, # same as general campylobacter
    co_etiology = if_else(any(c_across(ends_with("_tac_attr")) == 1), 1, 0)
  ) %>%
  ungroup()

# GEMS
# in code defined ETEC attributable and detected as ST or LT - not actually a TAC quantity for ETEC so can't make a transformed ETEC

gems_data <- gems_data %>%
  rename('cryptosporidium_new' = crypto_new,
         'astrovirus_new' = astro_new,
         'norovirus_gii_new' = noro_new,
         'campylobacter_new' = campy_new,
         'sapovirus_new' = sapo_new,
         'eaec_new' = EAEC_new,
         'shigella_tac_attr' = shigella_attributable_tac,
         'rotavirus_tac_attr' = rotavirus_attributable,
         'norovirus_gii_tac_attr' = noro_attributable,
         'adenovirus_tac_attr' = adenovirus_attributable,
         'sapovirus_tac_attr' = sapovirus_attributable,
         'astrovirus_tac_attr' = astro_attributable,
         'etec_tac_attr' = etec_attributable, 
         'tepec_tac_attr' = tepec_attributable,
         'cryptosporidium_tac_attr' = cryptosporidium_attributable,
         'v_cholerae_tac_attr' = v_cholerae_attributable,
         'c_jejuni_coli_tac_attr' = c_jejuni_coli_attributable,
         'salmonella_tac_attr' = salmonella_attributable,
         'st_etec_tac_attr' = st_etec_attributable) %>%
  rowwise() %>%
  mutate(
    co_etiology = if_else(any(c_across(ends_with("_tac_attr")) == 1), 1, 0)
  ) %>%
  ungroup()

# VIDA
# in data there's a binary etec alone and various TAC ETEC but unclear which if any are okay
vida_data <- vida_data %>%
  mutate(sapovirus_tac_attr = if_else(35 - (sapo_new * 3.322 ) < 21.9, 1, 0)) %>%
  rename('adenovirus_new' = adeno_new,
         'astrovirus_new' = astro_new,
         'norovirus_gii_new' = noro_gii_new,
         'campylobacter_new' = campy_new,
         'sapovirus_new' = sapo_new,
         'cryptosporidium_new' = crypto_new,
         'c_jejuni_coli_new' = campy_j_new,
         'shigella_tac_attr' = tac_shig,
         'rotavirus_tac_attr' = tac_rota,
         'norovirus_gii_tac_attr' = tac_norovirus_gii,
         'adenovirus_tac_attr' = tac_adeno_4041,
         'astrovirus_tac_attr' = tac_astro,
         'etec_tac_attr' = tac_st_etec, # look to see if we have lt etec in raw data
         'tepec_tac_attr' = tac_tepec,
         'cryptosporidium_tac_attr' = tac_crypto,
         'v_cholerae_tac_attr' = tac_vchol,
         'salmonella_tac_attr' = tac_salm,
         'c_jejuni_coli_tac_attr' = tac_campyj,
         'st_etec_tac_attr' = tac_st_etec) %>%
  rowwise() %>%
  mutate(
    campylobacter_new = c_jejuni_coli_new,
    co_etiology = if_else(any(c_across(ends_with("_tac_attr")) == 1), 1, 0)
  ) %>%
  ungroup()

# MAL-ED
maled_data <- maled_data %>%
  rename('adenovirus_new' = adenovirus_40_41_new,
         'etec_new' = ETEC_new,
         'tepec_new' = tEPEC_new,
         'campylobacter_new' = campylobacter_pan_new,
         'shigella_tac_attr' = shigella_attributable_tac,
         'rotavirus_tac_attr' = rotavirus_attributable,
         'norovirus_gii_tac_attr' = noro_gii_attributable,
         'adenovirus_tac_attr' = adenovirus_attributable,
         'sapovirus_tac_attr' = sapo_attributable,
         'astrovirus_tac_attr' = astro_attributable,
         'etec_tac_attr' = st_etec_attributable, # look to see if we have lt etec in raw data
         'tepec_tac_attr' = tepec_attributable,
         'cryptosporidium_tac_attr' = crypto_attributable,
         'v_cholerae_tac_attr' = v_cholerae_attributable,
         'salmonella_tac_attr' = salmonella_attributable,
         'c_jejuni_coli_tac_attr' = campylobacter_jejuni_coli_attributable) %>%
  rowwise() %>%
  mutate(
    st_etec_tac_attr = etec_tac_attr, 
    co_etiology = if_else(any(c_across(ends_with("_tac_attr")) == 1), 1, 0)
  ) %>%
  ungroup()

# Check no attributable etiology variable pathogens to include
# also v cholera, e histolytica, salmonella, isospora, aeromonas....
# ^^ it's ok if slightly different def of no_etiology

# ------------------------------------------------------------------------------
#                             SHIGELLA ATTR
# ------------------------------------------------------------------------------

################################################################################
# 1. TAC or culture attributable shigella (shig_attr)
# 2. Culture only Shigella (culture_shig_attr) (missing for ABCD)
################################################################################

# ABCD - note this is TAC only
abcd_data <- abcd_data %>%
  mutate('shig_attr' = shigella_tac_attr) %>%
  mutate(culture_shig_attr = rep(NA, nrow(abcd_data)))

# EFGH
efgh_data <- efgh_data %>%
  rename('shig_attr' = positive_tac_or_culture,
         'culture_shig_attr' = culture_shigella_positive)

# GEMS
gems_data <- gems_data %>%
  rename('shig_attr' = shigella_attributable,
         'culture_shig_attr' = shigella_attributable_culture)

# VIDA
vida_data <- vida_data %>%
  rename('shig_attr' = shigella_tac_or_culture,
         'culture_shig_attr' = shigella_culture_positive)

# MAL-ED
maled_data <- maled_data %>%
  rename('shig_attr' = shigella_attributable,
         'culture_shig_attr' = culture_shigella)

# ------------------------------------------------------------------------------
#                                    ANTIBIOTICS
# ------------------------------------------------------------------------------

################################################################################
# 1. Type of antibiotic received (all_abx) 
#.    (factor 3 levels - 0 = No or ineffective antibiotcs, 1 = Maybe effective antibiotics, 2 = WHO recommended antibiotics)
################################################################################

# ABCD
abcd_data$all_abx <- ifelse(abcd_data$an_grp_01 == 1, 2, 0)

# EFGH 
efgh_data$all_abx <- ifelse(efgh_data$all_abx == "No or ineffective antibiotics", 0,
                            ifelse(efgh_data$all_abx == "Possibly effective antibiotics", 1, 2))

# GEMS
gems_data$all_abx <- ifelse(gems_data$all_abx == "No or ineffective antibiotics", 0,
                            ifelse(gems_data$all_abx == "Possibly effective antibiotics", 1, 2))


# VIDA
vida_data$all_abx <- ifelse(vida_data$all_abx == "No or ineffective antibiotics", 0,
                            ifelse(vida_data$all_abx == "Possibly effective antibiotics", 1, 2))

# MAL-ED
maled_data$all_abx <- ifelse(maled_data$all_abx == "No or ineffective antibiotics", 0,
                             ifelse(maled_data$all_abx == "Possibly effective antibiotics", 1, 2))

# ------------------------------------------------------------------------------
#                                    OUTCOME
# ------------------------------------------------------------------------------

################################################################################
# 1. Followup HAZ (final_haz)
################################################################################

# ABCD
abcd_data <- abcd_data %>%
  rename('final_haz' = lazd90)

# EFGH
efgh_data <- efgh_data %>%
  rename('final_haz' = mo3_haz)

# GEMS
gems_data <- gems_data %>%
  rename('final_haz' = hazd60)

# VIDA
vida_data <- vida_data %>%
  rename('final_haz' = hazd60)

# MAL-ED
maled_data <- maled_data %>%
  rename('final_haz' = month3_haz)

# ------------------------------------------------------------------------------
#                             PATIENT ID VARIABLES
# ------------------------------------------------------------------------------

################################################################################
# 1. Unique patient identifier (first_id)
# 2. Identifier for episode (= first unique patient identifier) (child_id)
################################################################################

# ABCD - no re-enrollment
abcd_data <- abcd_data %>%
  mutate(child_id = pid) %>%
  rename('first_id' = pid)

# EFGH - pid is identifier for episode (fixed in IPD data 1/27/26)
efgh_data <- efgh_data %>%
  mutate(child_id = pid)

# GEMS - done
# VIDA - done
# MAL-ED - done

# ------------------------------------------------------------------------------
#                             Select and combine
# ------------------------------------------------------------------------------

# "shigella_new"
# "rotavirus_new"
# "adenovirus_new"
# "st_etec_new"
# "etec_new"
# "cryptosporidium_new" 
# "astrovirus_new"
# "norovirus_gii_new"
# "tepec_new"
# "campylobacter_new"
# "sapovirus_new"
# "giardia_new"
# "e_bieneusi_new"
# "eaec_new"

# ABCD 
abcd_data <- abcd_data %>%
  select(study,
         first_id,
         child_id,
         shig_attr,
         shigella_tac_attr, 
         culture_shig_attr,
         all_abx,
         sex,
         age,
         ses_quintile,
         enr_haz,
         site,
         edu_bin,
         num_hh_lt5,
         imp_water,
         imp_sanit,
         I_followup_days,
         I_followup_days_x_followup_days,
         dysentery,
         any_vom,
         any_fev,
         dehyd_level,
         lsstools,
         duration_pre_enroll,
         shigella_new,
         rotavirus_new,
         adenovirus_new,
         etec_new,
         st_etec_new,
         cryptosporidium_new,
         astrovirus_new,
         norovirus_gii_new,
         tepec_new,
         campylobacter_new,
         sapovirus_new,
         giardia_new,
         e_bieneusi_new,
         eaec_new,
         v_cholerae_new,
         salmonella_new,
         c_jejuni_coli_new,
         rotavirus_tac_attr,
         norovirus_gii_tac_attr,
         adenovirus_tac_attr, 
         sapovirus_tac_attr,
         astrovirus_tac_attr,
         # etec_tac_attr, 
         tepec_tac_attr,
         cryptosporidium_tac_attr, 
         v_cholerae_tac_attr, 
         salmonella_tac_attr,
         c_jejuni_coli_tac_attr,
         st_etec_tac_attr,
         co_etiology, 
         no_etiology,
         final_haz, 
         MSD)

# EFGH
efgh_data <- efgh_data %>%
  select(study,
         first_id,
         child_id,
         shig_attr,
         shigella_tac_attr,
         culture_shig_attr,
         all_abx,
         sex,
         age,
         ses_quintile,
         enr_haz,
         site,
         edu_bin,
         num_hh_lt5,
         imp_water,
         imp_sanit,
         I_followup_days,
         I_followup_days_x_followup_days,
         dysentery,
         any_vom,
         any_fev,
         dehyd_level,
         lsstools,
         duration_pre_enroll,
         shigella_new,
         rotavirus_new,
         adenovirus_new,
         etec_new,
         st_etec_new,
         cryptosporidium_new,
         astrovirus_new,
         norovirus_gii_new,
         tepec_new,
         campylobacter_new,
         sapovirus_new,
         giardia_new,
         e_bieneusi_new,
         eaec_new,
         v_cholerae_new,
         salmonella_new,
         c_jejuni_coli_new,
         rotavirus_tac_attr,
         norovirus_gii_tac_attr,
         adenovirus_tac_attr, 
         sapovirus_tac_attr,
         astrovirus_tac_attr,
         # etec_tac_attr, 
         tepec_tac_attr,
         cryptosporidium_tac_attr, 
         v_cholerae_tac_attr, 
         salmonella_tac_attr,
         c_jejuni_coli_tac_attr,
         st_etec_tac_attr,
         co_etiology, 
         no_etiology,
         final_haz, 
         MSD)

# GEMS
gems_data <- gems_data %>%
  select(study,
         first_id,
         child_id,
         shig_attr,
         shigella_tac_attr,
         culture_shig_attr,
         all_abx,
         sex,
         age,
         ses_quintile,
         enr_haz,
         site,
         edu_bin,
         num_hh_lt5,
         imp_water,
         imp_sanit,
         I_followup_days,
         I_followup_days_x_followup_days,
         dysentery,
         any_vom,
         any_fev,
         dehyd_level,
         lsstools,
         duration_pre_enroll,
         shigella_new,
         rotavirus_new,
         adenovirus_new,
         etec_new,
         st_etec_new,
         cryptosporidium_new,
         astrovirus_new,
         norovirus_gii_new,
         tepec_new,
         campylobacter_new,
         sapovirus_new,
         giardia_new,
         e_bieneusi_new,
         eaec_new,
         v_cholerae_new,
         salmonella_new,
         c_jejuni_coli_new,
         rotavirus_tac_attr,
         norovirus_gii_tac_attr,
         adenovirus_tac_attr, 
         sapovirus_tac_attr,
         astrovirus_tac_attr,
         # etec_tac_attr, 
         tepec_tac_attr,
         cryptosporidium_tac_attr, 
         v_cholerae_tac_attr, 
         salmonella_tac_attr,
         c_jejuni_coli_tac_attr,
         st_etec_tac_attr,
         co_etiology, 
         no_etiology,
         final_haz, 
         MSD)

# VIDA
vida_data <- vida_data %>%
  select(study,
         first_id,
         child_id,
         shig_attr,
         shigella_tac_attr,
         culture_shig_attr,
         all_abx,
         sex,
         age,
         ses_quintile,
         enr_haz,
         site,
         edu_bin,
         num_hh_lt5,
         imp_water,
         imp_sanit,
         I_followup_days,
         I_followup_days_x_followup_days,
         dysentery,
         any_vom,
         any_fev,
         dehyd_level,
         lsstools,
         duration_pre_enroll,
         shigella_new,
         rotavirus_new,
         adenovirus_new,
         etec_new,
         st_etec_new,
         cryptosporidium_new,
         astrovirus_new,
         norovirus_gii_new,
         tepec_new,
         campylobacter_new,
         sapovirus_new,
         giardia_new,
         e_bieneusi_new,
         eaec_new,
         v_cholerae_new,
         salmonella_new,
         c_jejuni_coli_new,
         rotavirus_tac_attr,
         norovirus_gii_tac_attr,
         adenovirus_tac_attr, 
         sapovirus_tac_attr,
         astrovirus_tac_attr,
         # etec_tac_attr, 
         tepec_tac_attr,
         cryptosporidium_tac_attr, 
         v_cholerae_tac_attr, 
         salmonella_tac_attr,
         c_jejuni_coli_tac_attr,
         st_etec_tac_attr,
         co_etiology, 
         no_etiology,
         final_haz, 
         MSD)

# MAL-ED
maled_data <- maled_data %>%
  ungroup() %>%
  select(study,
         first_id,
         child_id,
         shig_attr,
         shigella_tac_attr,
         culture_shig_attr,
         all_abx,
         sex,
         age,
         ses_quintile,
         enr_haz,
         site,
         edu_bin,
         num_hh_lt5,
         imp_water,
         imp_sanit,
         I_followup_days,
         I_followup_days_x_followup_days,
         dysentery,
         any_vom,
         any_fev,
         dehyd_level,
         lsstools,
         duration_pre_enroll,
         shigella_new,
         rotavirus_new,
         adenovirus_new,
         etec_new,
         st_etec_new,
         cryptosporidium_new,
         astrovirus_new,
         norovirus_gii_new,
         tepec_new,
         campylobacter_new,
         sapovirus_new,
         giardia_new,
         e_bieneusi_new,
         eaec_new,
         v_cholerae_new,
         salmonella_new,
         c_jejuni_coli_new,
         rotavirus_tac_attr,
         norovirus_gii_tac_attr,
         adenovirus_tac_attr, 
         sapovirus_tac_attr,
         astrovirus_tac_attr,
         # etec_tac_attr, 
         tepec_tac_attr,
         cryptosporidium_tac_attr, 
         v_cholerae_tac_attr, 
         salmonella_tac_attr,
         c_jejuni_coli_tac_attr,
         st_etec_tac_attr,
         co_etiology, 
         no_etiology,
         final_haz, 
         MSD)

# Combine all
combo_data <- rbind(abcd_data,
                    efgh_data,
                    gems_data,
                    vida_data,
                    maled_data)

combo_data$sex <- factor(combo_data$sex, levels = 0:1, labels = c("Male", "Female"))
combo_data$lsstools <- factor(combo_data$lsstools, levels = c("6 or less per day",
                                                              "7 to 10 per day",
                                                              "Over 10 per day"))

combo_data$all_abx <- factor(combo_data$all_abx, levels = 0:2, labels = c("No or ineffective antibiotics",
                                                                          "Possibly effective antibiotics",
                                                                          "Guideline recommended antibiotics"))
combo_data <- combo_data %>%
  set_variable_labels(study = "Study name",
                      first_id = "Identifier for child",
                      child_id = "Identifier for episode",
                      shig_attr = "TAC or culture Shigella attributable diarrhea",
                      shigella_tac_attr = "TAC attributable Shigella",
                      culture_shig_attr = "Culture Shigella attributable diarrhea",
                      all_abx = "Type of antibiotic received",
                      sex = "Sex",
                      age = "Age (months)",
                      ses_quintile = "Socioeconomic quintile",
                      enr_haz = "HAZ at enrollment",
                      site = "Site",
                      imp_water = "Improved water",
                      imp_sanit = "Improved sanitation",
                      num_hh_lt5 = "Number of children in household < age 5",
                      edu_bin = "Primary caregiver education > primary school",
                      I_followup_days = "Indicator followup days not missing",
                      I_followup_days_x_followup_days = "Indicator followup days not missing * followup days",
                      dysentery = "Dysentery",
                      any_vom = "Any vomitting",
                      any_fev = "Any fever",
                      dehyd_level = "Level of dehydration",
                      lsstools = "Number of loose stools",
                      duration_pre_enroll = "Duration of diarrhea before enrollment",
                      shigella_new = "TAC Shigella quantity",
                      rotavirus_new = "TAC rotavirus quantity",
                      adenovirus_new = "TAC adenovirus quantity",
                      etec_new = "TAC ETEC quantity",
                      cryptosporidium_new = "TAC cryptosporidium quantity",
                      astrovirus_new = "TAC astrovirus quantity",
                      norovirus_gii_new = "TAC norovirus GII quantity",
                      tepec_new = "TAC tEPEC quantity",
                      campylobacter_new = "TAC campylobacter quantity",
                      sapovirus_new = "TAC sapovirus quantity",
                      giardia_new = "TAC giardia quantity",
                      e_bieneusi_new = "TAC E Bieneusi quantity",
                      eaec_new = "TAC EAEC quantity",
                      v_cholerae_new = "TAC V Cholerae quantity",
                      st_etec_new = "TAC ST ETEC quantity",
                      salmonella_new = "TAC Salmonella quantity",
                      co_etiology = "Shigella diarrhea with co-etiology",
                      no_etiology = "Diarrhea with no known etiology",
                      final_haz = "HAZ at followup (day 60 or day 90)",
                      MSD = "GEMS definition of MSD (or as close as possible)")

if(full){
  saveRDS(combo_data, here::here("data/ipd_data/ipd_data_no_etiology_full.Rds"))
} else{
  if(MALED_MSD){
    saveRDS(combo_data, here::here("data/ipd_data/ipd_data_no_etiology.Rds"))
  } else{
    saveRDS(combo_data, here::here("data/ipd_data/ipd_data_no_etiology_inc_LSD.Rds"))
  }
}

#####################################################################
#                           CASE-CONTROL.                           #
#####################################################################

case_def <- "tac_or_culture_shig" # case definitions: tac_or_culture_shig, tac_shig, culture_shig, all_diar

gems_data <- readRDS(here::here(paste0("data/gems_data/gems_case_control_", case_def, ".Rds")))
maled_data <- readRDS(here::here(paste0("data/maled_data/maled_case_control_", case_def, ".Rds")))
vida_data <- readRDS(here::here(paste0("data/vida_data/vida_case_control_", case_def, ".Rds")))

# Add column with study name to each dataset
gems_data$study <- "GEMS"
vida_data$study <- "VIDA"
maled_data$study <- "MALED"

msd <- FALSE # subset MAL-ED to MSD cases
lsd <- FALSE

# MAL-ED remove LSD cases & matched controls 
# maled_LSD_cases <- maled_data$case_sid[which(maled_data$case == 1 & maled_data$MSD == 0)]
# maled_data <- maled_data[-which(maled_data$case_sid %in% maled_LSD_cases),]

# ------------------------------------------------------------------------------
#                                 COVARIATES
# ------------------------------------------------------------------------------

################################################################################
# 1. Sex (sex) - Male = 0, Female = 1
# 2. Age in months (age)
# 3. SES quintile (ses_quintile)
# 4. Baseline HAZ (enr_haz)
# 5. Days between episode and follow-up (followup_days) -- I_followup_days & I_followup_days_x_followup_days
# 6. Site (site) - but make it a factor for each study site individually?

# Missing in some (fill in with 0)
# 7. Water + sanitation score (missing in ABCD) (imp_water, imp_sanit)
# 8. Number of children in household < 5 (missing in MAL-ED) (num_hh_lt5)
# 9. Primary caregiver education (missing in ABCD) (edu_bin) 
################################################################################

# GEMS 
# naming conventions all match
gems_data$sex <- ifelse(gems_data$sex == "female", 1, 0)
levels(gems_data$site) <- c("The_Gambia_gems",
                            "Mali_gems",
                            "Mozambique_gems",
                            "Kenya_gems",
                            "India_gems",
                            "Bangladesh_gems",
                            "Pakistan_gems")

gems_data$edu_bin <- gems_data$prim_caregiver_edu_bin

gems_data$imp_water <- ifelse(gems_data$safe_water %in% c("Safely managed", "Basic"), 1, 0)
gems_data$imp_sanit <- ifelse(gems_data$safe_sanit %in% c("Safely manaed and basic"), 1, 0)

# VIDA
vida_data <- vida_data %>%
  rename('age' = agemchild,
         'edu_bin' = education_bin)

vida_data$sex <- ifelse(vida_data$sex == "Female", 1, 0)
levels(vida_data$site) <- c("The_Gambia_vida",
                            "Mali_vida",
                            "Kenya_vida")

vida_data$imp_water <- ifelse(vida_data$safe_water %in% c("Safely managed", "Basic"), 1, 0)
vida_data$imp_sanit <- ifelse(vida_data$safe_sanit %in% c("Safely managed and basic"), 1, 0)

# MAL-ED
maled_data <- maled_data %>%
  rename('age' = agemonths,
         'ses_quintile' = wami_quintile,
         'enr_haz' = baseline_haz,
         'imp_water' = drinkimp,
         'imp_sanit' = sanitimp,
         'edu_bin' = mated_bin)

maled_data$sex <- ifelse(maled_data$sex == "female", 1, 0)
levels(maled_data$site) <- c("Bangladesh_maled",
                             "Brazil_maled",
                             "India_maled",
                             "Nepal_maled",
                             "Peru_maled",
                             "SouthAfrica_maled",
                             "Tanzania_maled")

maled_data$num_hh_lt5 <- 1

# ------------------------------------------------------------------------------
#                                 SEVERITY
# ------------------------------------------------------------------------------

################################################################################
# 1. Dysentery (0 = no, 1 = yes) (dysentery)
# 2. Vomiting during episode (0 = no, 1 = yes) (any_vom)
# 3. Dehydration level (Factor with levels No dehydration, Some dehydration, Severe dehydration) (dehyd_level)
# 4. Max number of loose stools (categorical) (Factor with levels 6 or less, 7 to 10, Over 10) (except vida is really 3 to 5, 6 to 10, over 10) (lsstools)
# 5. Duration of diarrhea pre-abx (duration_pre_enroll)

# 6. Fever during episode (0 = no, 1 = yes) (any_fev)
################################################################################

# GEMS
gems_data <- gems_data %>%
  rename('any_vom' = vomit,
         'any_fev' = fever,
         'dehyd_level' = who_dehyd)

# VIDA
vida_data <- vida_data %>%
  mutate(any_vom = if_else(vom_days > 0, 1, 0),
         dehyd_level = if_else(dehydr == "None", "No dehydration", dehydr)) %>%
  rename('any_fev' = fever)

vida_data$lsstools <- ifelse(vida_data$lsstools == "3 stools" | vida_data$lsstools == "4 to 5 stools", "6 or less per day",
                             ifelse(vida_data$lsstools == "6 to 10 stools", "7 to 10 per day", "Over 10 per day"))

# MAL-ED
maled_data <- maled_data %>%
  mutate(any_vom = if_else(daysvomit > 0, 1, 0),
         dehyd_level = if_else(dehyd == "None", "No dehydration", dehyd),
         lsstools = if_else(lsstools <=6, "6 or less per day",
                            if_else(lsstools >= 7 & lsstools <= 10, "7 to 10 per day", "Over 10 per day"))) %>%
  rename('any_fev' = fever,
         'duration_pre_enroll' = duration_pre_abx)

# ------------------------------------------------------------------------------
#                             PATHOGEN QUANTITY
# ------------------------------------------------------------------------------

################################################################################
# Final list / naming conventions:

# "shigella_new"
# "rotavirus_new"
# "adenovirus_new"
# "st_etec_new"
# "etec_new"
# "cryptosporidium_new" 
# "astrovirus_new"
# "norovirus_new"
# "tepec_new"
# "campylobacter_new"
# "sapovirus_new"
# "giardia_new"
# "e_bieneusi_new"
# "eaec_new"
################################################################################

# GEMS
# in code defined ETEC attributable and detected as ST or LT - not actually a TAC quantity for ETEC so can't make a transformed ETEC

gems_data <- gems_data %>%
  rename('cryptosporidium_new' = crypto_new,
         'astrovirus_new' = astro_new,
         'norovirus_new' = noro_new,
         'campylobacter_new' = campy_new,
         'sapovirus_new' = sapo_new,
         'eaec_new' = EAEC_new)

# VIDA
# in data there's a binary etec alone and various TAC ETEC but unclear which if any are okay
vida_data <- vida_data %>%
  rename('adenovirus_new' = adeno_new,
         'astrovirus_new' = astro_new,
         'norovirus_new' = noro_new,
         'campylobacter_new' = campy_new,
         'sapovirus_new' = sapo_new,
         'cryptosporidium_new' = crypto_new)

# MAL-ED
maled_data <- maled_data %>%
  rename('adenovirus_new' = adenovirus_40_41_new,
         'etec_new' = ETEC_new,
         'tepec_new' = tEPEC_new,
         'campylobacter_new' = campylobacter_pan_new)

# Check no attributable etiology variable pathogens to include
# also v cholera, e histolytica, salmonella, isospora, aeromonas....
# ^^ it's ok if slightly different def of no_etiology

# ------------------------------------------------------------------------------
#                             SHIGELLA ATTR
# ------------------------------------------------------------------------------

################################################################################
# 1. TAC or culture attributable shigella (shig_attr)
# 2. culture attributable shigella

# Note we use 'case' variable instead but should already be 'case' for all
# this is just getting a uniform variable for tac_or_culture shigella and culture shigella
################################################################################

# GEMS
gems_data <- gems_data %>%
  rename('shig_attr' = shigella_attributable,
         'culture_shig_attr' = shigella_culture_pos,
         'tac_shig_attr' = shigella_attributable_tac)

# VIDA
vida_data <- vida_data %>%
  rename('shig_attr' = shigella_tac_or_culture,
         'culture_shig_attr' = shigella_culture_positive,
         'tac_shig_attr' = tac_shig)

# MAL-ED
maled_data <- maled_data %>%
  rename('shig_attr' = shigella_attributable,
         'culture_shig_attr' = culture_shigella,
         'tac_shig_attr' = tac_shigella_attributable)

################################################################################
# 1. Type of antibiotic received (all_abx) 
#.    (factor 3 levels - 0 = No or ineffective antibiotcs, 1 = Maybe effective antibiotics, 2 = WHO recommended antibiotics)
################################################################################

# GEMS
gems_data$all_abx <- ifelse(gems_data$all_abx == "No/Ineffective abx", 0,
                            ifelse(gems_data$all_abx == "Maybe effective abx", 1, 2))


# VIDA
vida_data$all_abx <- ifelse(vida_data$all_abx == "Ineffective or no abx", 0,
                            ifelse(vida_data$all_abx == "Maybe effective abx", 1, 2))

# MAL-ED
maled_data$all_abx <- ifelse(maled_data$all_abx == "No or ineffective antibiotics", 0,
                             ifelse(maled_data$all_abx == "Possibly effective antibiotics", 1, 2))
# ------------------------------------------------------------------------------
#                                    OUTCOME
# ------------------------------------------------------------------------------

################################################################################
# 1. Followup HAZ (final_haz)
################################################################################

# GEMS
gems_data <- gems_data %>%
  rename('final_haz' = hazd60)

# VIDA
vida_data <- vida_data %>%
  rename('final_haz' = hazd60)

# MAL-ED
maled_data <- maled_data %>%
  rename('final_haz' = month3_haz)

# ------------------------------------------------------------------------------
#                             PATIENT ID VARIABLES
# ------------------------------------------------------------------------------

################################################################################
# 1. Unique patient identifier (first_id)
# 2. Identifier for episode (= first unique patient identifier) (child_id)
# 3. Identifier for case (case_id, if case case_id = child_id)
################################################################################

# GEMS - done - make character
gems_data$child_id <- as.character(gems_data$child_id)
gems_data$first_id <- as.character(gems_data$first_id)
gems_data$case_id <- as.character(gems_data$case_id)

# VIDA - done - make character
vida_data$child_id <- as.character(vida_data$child_id)
vida_data$first_id <- as.character(vida_data$first_id)
vida_data$case_id <- as.character(vida_data$case_id)

# MAL-ED- done

# ------------------------------------------------------------------------------
#                             Select and combine
# ------------------------------------------------------------------------------

# "shigella_new"
# "rotavirus_new"
# "adenovirus_new"
# "st_etec_new"
# "etec_new"
# "cryptosporidium_new" 
# "astrovirus_new"
# "norovirus_new"
# "tepec_new"
# "campylobacter_new"
# "sapovirus_new"
# "giardia_new"
# "e_bieneusi_new"
# "eaec_new"

# GEMS
gems_data <- gems_data %>%
  select(study,
         first_id,
         child_id,
         case_id, 
         case, 
         shig_attr,
         culture_shig_attr,
         tac_shig_attr,
         all_abx,
         sex,
         age,
         ses_quintile,
         enr_haz,
         site,
         edu_bin,
         num_hh_lt5,
         imp_water,
         imp_sanit,
         I_followup_days,
         I_followup_days_x_followup_days,
         dysentery,
         any_vom,
         any_fev,
         dehyd_level,
         lsstools,
         duration_pre_enroll,
         shigella_new,
         rotavirus_new,
         adenovirus_new,
         etec_new,
         cryptosporidium_new,
         astrovirus_new,
         norovirus_new,
         tepec_new,
         campylobacter_new,
         sapovirus_new,
         giardia_new,
         e_bieneusi_new,
         eaec_new,
         final_haz)

# VIDA
vida_data <- vida_data %>%
  select(study,
         first_id,
         child_id,
         case_id, 
         case, 
         shig_attr,
         culture_shig_attr,
         tac_shig_attr,
         all_abx,
         sex,
         age,
         ses_quintile,
         enr_haz,
         site,
         edu_bin,
         num_hh_lt5,
         imp_water,
         imp_sanit,
         I_followup_days,
         I_followup_days_x_followup_days,
         dysentery,
         any_vom,
         any_fev,
         dehyd_level,
         lsstools,
         duration_pre_enroll,
         shigella_new,
         rotavirus_new,
         adenovirus_new,
         etec_new,
         cryptosporidium_new,
         astrovirus_new,
         norovirus_new,
         tepec_new,
         campylobacter_new,
         sapovirus_new,
         giardia_new,
         e_bieneusi_new,
         eaec_new,
         final_haz)

#  SUBSET TO MSD MAL-ED
if(msd){
  # Get IDs corresponding to MSD cases
  maled_data_msd_caseids <- unique(maled_data$case_id[which(maled_data$MSD == 1)])
  
  # Subset to MSD cases and their matched controls
  maled_data <- maled_data[which(maled_data$case_id %in% maled_data_msd_caseids),]
}
if(lsd){
  # Get IDs corresponding to MSD cases and their matched controls
  maled_data_msd_caseids <- unique(maled_data$case_id[which(maled_data$MSD == 1)])
  
  # Remove from data 
  # (can't just do MSD = 0 because then will still pull MSD cases because matched controls have MSD 0)
  maled_data <- subset(maled_data, !(case_id %in% maled_data_msd_caseids))
}

# MAL-ED
maled_data <- maled_data %>%
  ungroup() %>%
  select(study,
         first_id,
         child_id,
         case_id, 
         case, 
         shig_attr,
         culture_shig_attr,
         tac_shig_attr,
         all_abx,
         sex,
         age,
         ses_quintile,
         enr_haz,
         site,
         edu_bin,
         num_hh_lt5,
         imp_water,
         imp_sanit,
         I_followup_days,
         I_followup_days_x_followup_days,
         dysentery,
         any_vom,
         any_fev,
         dehyd_level,
         lsstools,
         duration_pre_enroll,
         shigella_new,
         rotavirus_new,
         adenovirus_new,
         etec_new,
         cryptosporidium_new,
         astrovirus_new,
         norovirus_new,
         tepec_new,
         campylobacter_new,
         sapovirus_new,
         giardia_new,
         e_bieneusi_new,
         eaec_new,
         final_haz)

# Combine all
if(!lsd){
  combo_data <- rbind(gems_data,
                      vida_data,
                      maled_data)
} else {
  combo_data <- maled_data
}


combo_data$sex <- factor(combo_data$sex, levels = 0:1, labels = c("Male", "Female"))
combo_data$lsstools <- factor(combo_data$lsstools, levels = c("6 or less per day",
                                                              "7 to 10 per day",
                                                              "Over 10 per day"))

combo_data$all_abx <- factor(combo_data$all_abx, levels = 0:2, labels = c("No or ineffective antibiotics",
                                                                          "Maybe effective antibiotics",
                                                                          "WHO recommended antibiotics"))

combo_data$dehyd_level <- factor(combo_data$dehyd_level, levels = c("No dehydration", "Some dehydration", "Severe dehydration"))

combo_data <- combo_data %>%
  set_variable_labels(study = "Study name",
                      first_id = "Identifier for child",
                      child_id = "Identifier for episode",
                      case_id = "Identifier for case episode",
                      case = "Case (=1 case, =0 control)",
                      shig_attr = "TAC or culture Shigella attributable diarrhea",
                      tac_shig_attr = "TAC attributable Shigella diarrhea",
                      culture_shig_attr = "Culture attributable Shigella diarrhea",
                      all_abx = "Type of antibiotic received",
                      sex = "Sex",
                      age = "Age (months)",
                      ses_quintile = "Socioeconomic quintile",
                      enr_haz = "HAZ at enrollment",
                      site = "Site",
                      imp_water = "Improved water",
                      imp_sanit = "Improved sanitation",
                      num_hh_lt5 = "Number of children in household < age 5",
                      edu_bin = "Primary caregiver education > primary school",
                      I_followup_days = "Indicator followup days not missing",
                      I_followup_days_x_followup_days = "Indicator followup days not missing * followup days",
                      dysentery = "Dysentery",
                      any_vom = "Any vomitting",
                      any_fev = "Any fever",
                      dehyd_level = "Level of dehydration",
                      lsstools = "Number of loose stools",
                      duration_pre_enroll = "Duration of diarrhea before enrollment",
                      shigella_new = "TAC Shigella quantity",
                      rotavirus_new = "TAC rotavirus quantity",
                      adenovirus_new = "TAC adenovirus quantity",
                      etec_new = "TAC ETEC quantity",
                      cryptosporidium_new = "TAC cryptosporidium quantity",
                      astrovirus_new = "TAC astrovirus quantity",
                      norovirus_new = "TAC norovirus quantity",
                      tepec_new = "TAC tEPEC quantity",
                      campylobacter_new = "TAC campylobacter quantity",
                      sapovirus_new = "TAC sapovirus quantity",
                      giardia_new = "TAC giardia quantity",
                      e_bieneusi_new = "TAC E Bieneusi quantity",
                      eaec_new = "TAC EAEC quantity",
                      final_haz = "HAZ at followup (day 60 or day 90)")

if(msd){
  saveRDS(combo_data, here::here(paste0("data/ipd_data/ipd_data_case_control_msd_", case_def,".Rds")))
} else if(lsd){
  saveRDS(combo_data, here::here(paste0("data/ipd_data/ipd_data_case_control_lsd_", case_def,".Rds")))
} else {
  saveRDS(combo_data, here::here(paste0("data/ipd_data/ipd_data_case_control_", case_def, ".Rds")))
}

