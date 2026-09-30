# ---------------------------------------------------------------
# Script to create cleaned ABCD dataset
# ---------------------------------------------------------------

here::i_am("code/prep_data/prep_abcd.R")

library(tidyverse)
library(labelled)
library(haven)
library(corrr)
library(factoextra)

prep_abcd <- function(){
  # only use for vomitting variable
  abcd_data_otr <- read.csv(here::here("data/abcd_data/raw_data/ABCD_OTR_26dec2023.csv"))
  abcd_data_otr <- abcd_data_otr[, colnames(abcd_data_otr) %in% c("pid", "dy1_scrn_vomitall")]
  abcd_data_otr$pid <- as.character(abcd_data_otr$pid)
  
  # raw data + merge in vomitting var
  full_abcd_data <- read_dta(here::here("data/abcd_data/raw_data/ABCD_TAC_withnewvars_23Mar2022.dta"))
  full_abcd_data <- full_abcd_data %>%
    left_join(abcd_data_otr, by = "pid")
  
  # NEW 9/18/25 - make GEMS MSD variable
  # at least 1 of the following
  # sunken eyes (redcap dictionary says have for deaths?? but doesn't seem like overall)
  # loss of skin turgor (recap dictionary says dy1_scrn_sknrtrnnorm but also can't find in my rawest data)
  # IV rehydratoin (redcap dictionary says d\1_scrn_stabilneed but can't find)
  # dysentery (exclusion criteria here)
  # hospitalized (redcap dictioniary an_hosp10_yn first 10 days, closest to during episode??)
  
  # for now, use: 
  # hosp within 10 days and/or severe dehydration
  full_abcd_data$gems_msd <- ifelse(full_abcd_data$an_hosp10_yn == 1 | full_abcd_data$dy1_scrn_dehydr > 1, 1, 0)
  
  # Get new vars from full data
  sub_abcd_data <- full_abcd_data[,c("pid", 
                                     "sample_type", # 1 = stool, 2 = swab
                                     "an_d90_timing", 
                                     
                                     "lazd90",
                                     "wlzd90",
                                     "wazd90",
                                     
                                     "day2diar",
                                     "day3diar",
                                     "agemchild" ,
                                     "month_en" ,
                                     "site" ,
                                     "dy1_ant_sex" ,
                                     "an_ses_quintile" ,
                                     "avemuac" ,
                                     
                                     "lfazscore" ,
                                     "wflzscore",
                                     "wfazscore",
                                     
                                     "an_tothhlt5" ,
                                     "an_grp_01" ,
                                     
                                     "dy1_scrn_vomitall",
                                     "dy1_scrn_dehydr",
                                     "dy1_scrn_lstools",
                                     "dy1_scrn_diardays",
                                     
                                     # Add all pathogens from raw data
                                     # do not have an adjusted version for these
                                     "campylobacter_pan",
                                     "norovirus",
                                     "lt_etec",
                                     "etec",
                                     "eaec", 
                                     "giardia", 
                                     "e_bieneusi",
                                     
                                     # adjusted for swab samples
                                     "adenovirus_40_41_adjust",
                                     "astrovirus_adjust",
                                     "campylobacter_jejuni_coli_adjust",
                                     "cryptosporidium_adjust",
                                     "e_histolytica_adjust",
                                     "norovirus_gii_adjust",
                                     "aeromonas_adjust",
                                     "isospora_adjust",
                                     "cyclospora_adjust",
                                     "v_cholerae_adjust",
                                     "rotavirus_adjust",
                                     "salmonella_adjust",               
                                     "sapovirus_adjust",
                                     "shigella_eiec_adjust",
                                     "st_etec_adjust",
                                     "tepec_adjust",  
                                     
                                     "wlzdiff",
                                     "wazdiff",
                                     "lazdiff" ,
                                     
                                     "an_meanlen_d1",
                                     "an_meanwt_d1",
                                     "an_meanlen_d90",
                                     "an_meanwt_d90",
                                     "an_hosp90death", # add hospitalisation or death by D90
                                     "an_hosp90_yn", # hospitalization by 90
                                     "an_death90",   # death by 90
                                     "gems_msd")]
  
  sub_abcd_data$pid <- as.numeric(sub_abcd_data$pid)
  
  # make difference vars for length and weight
  sub_abcd_data$lendiff <- sub_abcd_data$an_meanlen_d90 - sub_abcd_data$an_meanlen_d1
  sub_abcd_data$wtdiff <- sub_abcd_data$an_meanwt_d90 - sub_abcd_data$an_meanwt_d1
  
  sub_abcd_data <- sub_abcd_data %>%
    set_variable_labels(lendiff = "Difference in average day 90 length and day 1 length",
                        wtdiff = "Difference in average day 90 weight and day 1 weight")
  
  ### NEW 2/10/26: redo pathogen quantity variables starting at raw data ###
  
  # transformation equivalent to log(2^(35 - sub_abcd_data$shigella_eiec), base = 10)
  
  # shigella cutoff 28.7
  sub_abcd_data$shigella_likely <- ifelse(sub_abcd_data$shigella_eiec_adjust < 28.7, 1, 0)
  sub_abcd_data$shigella_new <- (35 - sub_abcd_data$shigella_eiec_adjust) / 3.322
  
  # adenovirus cutoff 23.6
  sub_abcd_data$adenovirus_likely <- ifelse(sub_abcd_data$adenovirus_40_41_adjust < 23.6, 1, 0)
  sub_abcd_data$adenovirus_new <- (35 - sub_abcd_data$adenovirus_40_41_adjust) / 3.322
  
  # astrovirus cutoff 25.2
  sub_abcd_data$astrovirus_likely <- ifelse(sub_abcd_data$astrovirus_adjust < 25.2, 1, 0)
  sub_abcd_data$astrovirus_new <- (35 - sub_abcd_data$astrovirus_adjust) / 3.322
  
  # campy jejuni coli cutoff 16.3
  sub_abcd_data$c_jejuni_likely <- ifelse(sub_abcd_data$campylobacter_jejuni_coli_adjust < 16.3, 1, 0)
  sub_abcd_data$c_jejuni_new <- (35 - sub_abcd_data$campylobacter_jejuni_coli_adjust) / 3.322
  
  # cryptosporidium cutoff 25.6
  sub_abcd_data$cryptosporidium_likely <- ifelse(sub_abcd_data$cryptosporidium_adjust < 25.6, 1, 0)
  sub_abcd_data$cryptosporidium_new <- (35 - sub_abcd_data$cryptosporidium_adjust) / 3.322
  
  # cyclospora cutoff 33.8
  sub_abcd_data$cyclospora_likely <- ifelse(sub_abcd_data$cyclospora_adjust < 33.8, 1, 0)
  sub_abcd_data$cyclospora_new <- (35 - sub_abcd_data$cyclospora_adjust) / 3.322
  
  # e hist cutoff 26.6
  sub_abcd_data$e_hist_likely <- ifelse(sub_abcd_data$e_histolytica_adjust < 26.6, 1, 0)
  sub_abcd_data$e_hist_new <- (35 - sub_abcd_data$e_histolytica_adjust) / 3.322
  
  # isospora cutoff 32.6
  sub_abcd_data$isospora_likely <- ifelse(sub_abcd_data$isospora_adjust < 32.6, 1, 0)
  sub_abcd_data$isospora_new <- (35 - sub_abcd_data$isospora_adjust) / 3.322
  
  # norovirus cutoff 25.7
  sub_abcd_data$norovirus_gii_likely <- ifelse(as.numeric(sub_abcd_data$norovirus_gii_adjust) < 25.7, 1, 0)
  sub_abcd_data$norovirus_gii_new <- (35 - as.numeric(sub_abcd_data$norovirus_gii_adjust)) / 3.322
  
  # rotavirus cutoff 32
  sub_abcd_data$rotavirus_likely <- ifelse(sub_abcd_data$rotavirus_adjust < 32, 1, 0)
  sub_abcd_data$rotavirus_new <- (35 - sub_abcd_data$rotavirus_adjust) / 3.322
  
  # salmonella cutoff 31.9
  sub_abcd_data$salmonella_likely <- ifelse(sub_abcd_data$salmonella_adjust < 31.9, 1, 0)
  sub_abcd_data$salmonella_new <- (35 - sub_abcd_data$salmonella_adjust) / 3.322
  
  # sapovirus cutoff 22.5
  sub_abcd_data$sapovirus_likely <- ifelse(sub_abcd_data$sapovirus_adjust < 22.5, 1, 0)
  sub_abcd_data$sapovirus_new <- (35 - sub_abcd_data$sapovirus_adjust) / 3.322
  
  # st-etec cutoff 25.4
  sub_abcd_data$st_etec_likely <- ifelse(sub_abcd_data$st_etec_adjust < 25.4, 1, 0)
  sub_abcd_data$st_etec_new <- (35 - sub_abcd_data$st_etec_adjust) / 3.322
  
  # tepec cutoff 18.1
  sub_abcd_data$tepec_likely <- ifelse(sub_abcd_data$tepec_adjust < 18.1, 1, 0)
  sub_abcd_data$tepec_new <- (35 - sub_abcd_data$tepec_adjust) / 3.322
  
  # v cholerae cutoff 32.6
  sub_abcd_data$v_cholerae_likely <- ifelse(sub_abcd_data$v_cholerae_adjust < 32.6, 1, 0)
  sub_abcd_data$v_cholerae_new <- (35 - sub_abcd_data$v_cholerae_adjust) / 3.322
  
  # pathogens of interest without explicit cutoffs (unadjusted)
  
  # LT etec
  sub_abcd_data$lt_etec_adjust <- ifelse(sub_abcd_data$sample_type == 2 & sub_abcd_data$lt_etec < 35, sub_abcd_data$lt_etec - 2.78, sub_abcd_data$lt_etec)
  sub_abcd_data$lt_etec_new <- (35 - sub_abcd_data$lt_etec_adjust) / 3.322
  
  # ETEC (general) 
  sub_abcd_data$etec_new <- ifelse(sub_abcd_data$lt_etec_new > sub_abcd_data$st_etec_new, sub_abcd_data$lt_etec_new, sub_abcd_data$st_etec_new)
  
  # campylobacter (campylobacter_pan) note we want this instead of jejuni and coli but cutoff seems to be for the other one??
  # unclear if adjustment is for this either but applying it 
  sub_abcd_data$campylobacter_adjust <- ifelse(sub_abcd_data$sample_type == 2 & as.numeric(sub_abcd_data$campylobacter_pan) < 35, as.numeric(sub_abcd_data$campylobacter_pan) - 0.81, as.numeric(sub_abcd_data$campylobacter_pan))
  sub_abcd_data$campylobacter_new <- (35 - as.numeric(sub_abcd_data$campylobacter_adjust)) / 3.322
  
  # E Bienusi - no adjustment
  sub_abcd_data$e_bieneusi_new <- (35 - as.numeric(sub_abcd_data$e_bieneusi)) / 3.322
  
  # giardia
  sub_abcd_data$giardia_adjust <- ifelse(sub_abcd_data$sample_type == 2 & as.numeric(sub_abcd_data$giardia) < 35, as.numeric(sub_abcd_data$giardia) - 3.88, as.numeric(sub_abcd_data$giardia))
  sub_abcd_data$giardia_new <- (35 - as.numeric(sub_abcd_data$giardia)) / 3.322
  
  # eaec
  sub_abcd_data$eaec_adjust <- ifelse(sub_abcd_data$sample_type == 2 & sub_abcd_data$eaec < 35, sub_abcd_data$eaec - 0.95, sub_abcd_data$eaec)
  sub_abcd_data$eaec_new <- (35 - as.numeric(sub_abcd_data$eaec)) / 3.322
  
  # Make no_etiology variable using the pathogens of interest we have likely etiology for
  sub_abcd_data$no_etiology <- ifelse(rowSums(sub_abcd_data[,c("shigella_likely",
                                                               "adenovirus_likely",
                                                               "astrovirus_likely",
                                                               # "campylobacter_likely",
                                                               "c_jejuni_likely",
                                                               "cryptosporidium_likely",
                                                               "cyclospora_likely",
                                                               "e_hist_likely",
                                                               "isospora_likely",
                                                               "norovirus_gii_likely",
                                                               "rotavirus_likely",
                                                               "salmonella_likely",
                                                               "sapovirus_likely",
                                                               "st_etec_likely",
                                                               "tepec_likely",
                                                               "v_cholerae_likely")], na.rm = TRUE) == 0, 1, 0)
  
  # also make the _bin variables that were once in the OTR dataset
  # List of pathogens to make _bin variables for
  pathogens <- c(
    "shigella", "adenovirus", "astrovirus", "c_jejuni", "cryptosporidium",
    "cyclospora", "e_hist", "isospora", "norovirus_gii", "rotavirus",
    "salmonella", "sapovirus", "st_etec", "tepec", "v_cholerae",
    "etec", "campylobacter", "e_bieneusi", "giardia", "eaec"
  )
  
  for (p in pathogens) {
    new_var <- paste0(p, "_new")
    bin_var <- paste0(p, "_bin")
    
    # If _new exists, use _new > 0; otherwise fallback to raw < 35
    if (new_var %in% names(sub_abcd_data)) {
      sub_abcd_data[[bin_var]] <- ifelse(sub_abcd_data[[new_var]] > 0, 1, 0)
    } else {
      warning(paste("No variable found for", p))
    }
  }
  
  # drop existing _new and _bin variables from abcd_data otr dataset
  # abcd_data <- abcd_data[, !grepl("_new$", names(abcd_data))]
  # abcd_data <- abcd_data[, !grepl("_bin$", names(abcd_data))]
  # 
  # # join raw data with new tac variables to OTR data
  # abcd_data <- left_join(abcd_data, sub_abcd_data, by = "pid")
  abcd_data <- sub_abcd_data
  
  abcd_data$I_an_d90_timing <- ifelse(is.na(abcd_data$an_d90_timing), 0, 1)
  abcd_data$I_an_d90_timing_x_an_d90_timing <- ifelse(is.na(abcd_data$an_d90_timing), 0, abcd_data$an_d90_timing)
  
  # eliminate laz and waz >6 or <-6
  abcd_data$lazd90 <- ifelse(abcd_data$lazd90 > 6 | abcd_data$lazd90 < -6, NA, abcd_data$lazd90)
  abcd_data$wlzd90 <- ifelse(abcd_data$wlzd90 > 6 | abcd_data$wlzd90 < -6, NA, abcd_data$wlzd90)
  abcd_data$wazd90 <- ifelse(abcd_data$wazd90 > 6 | abcd_data$wazd90 < -6, NA, abcd_data$wazd90)
  
  abcd_data$lfazscore <- ifelse(abcd_data$lfazscore > 6 | abcd_data$lfazscore < -6, NA, abcd_data$lfazscore)
  abcd_data$wfazscore <- ifelse(abcd_data$wfazscore > 6 | abcd_data$wfazscore < -6, NA, abcd_data$wfazscore)
  abcd_data$wflzscore <- ifelse(abcd_data$wflzscore > 6 | abcd_data$wflzscore < -6, NA, abcd_data$wflzscore)
  
  abcd_data$lazdiff <- ifelse(is.na(abcd_data$lazd90) | is.na(abcd_data$lfazscore), NA, abcd_data$lazdiff)
  
  # add stunting var based on lazd90
  abcd_data$stunting_bin <- ifelse(abcd_data$lazd90 < -2, 1, 0)
  
  # SAME PREP FROM OTR ANALYSIS
  abcd_data$dy1_scrn_vomitall <- factor(abcd_data$dy1_scrn_vomitall, levels=c(1,2), 
                                        labels=c("Yes", "No"))
  
  abcd_data$dy1_scrn_dehydr <- factor(abcd_data$dy1_scrn_dehydr, levels=c(1,2,3), 
                                      labels=c("No dehydration", "Some dehydration", "Severe dehydration"))
  
  abcd_data$site <- factor(abcd_data$site, levels=c(2:8),
                           labels=c("Bangladesh", "Kenya", "Malawi", "Mali", "India", "Tanzania", "Pakistan"))
  
  abcd_data$dy1_ant_sex <- factor(abcd_data$dy1_ant_sex, levels=c(1,2),
                                  labels=c("Male", "Female"))
  
  abcd_data$an_ses_quintile <- factor(abcd_data$an_ses_quintile, levels=c(1:5),
                                      labels = c("1st quintile of SES",
                                                 "2nd quintile of SES",
                                                 "3rd quintile of SES",
                                                 "4th quintile of SES",
                                                 "5th quintile of SES"))
  
  abcd_data$month_en <- factor(abcd_data$month_en, levels = c(1:12),
                               labels = c("January",
                                          "February",
                                          "March",
                                          "April",
                                          "May",
                                          "June",
                                          "July",
                                          "August",
                                          "September",
                                          "October",
                                          "November",
                                          "December"))
  
  # USE PATHOGEN CUTOFFS FROM ABCD SUPPLEMENT TO MAKE LIKELY VARS BY TAC
  # "Median value of GEMS/MAL-ED cut offs ("likely etiology" cutoffs used in present study)
  # Transformation for "shigella_new" was (35 - CTValue) / 3.322
  # Transform back to original scale
  
  # added (1/21/25)
  # Variable for detected (rather than likely/attributable) == pathogen_bin 
  # already in cleaned data, if not in cleaned data make var
  
  # shigella cutoff 28.7
  # transform back to original scale (though now we want the new scale ones anyways...)
  # abcd_data$shigella_tac <- 35 - 3.322*abcd_data$shigella_new 
  # abcd_data$shigella_likely <- ifelse(abcd_data$shigella_tac < 28.7, 1, 0)
  # 
  # # adenovirus cutoff 23.6
  # abcd_data$adenovirus_tac <- 35 - 3.322*abcd_data$adenovirus_new
  # abcd_data$adenovirus_likely <- ifelse(abcd_data$adenovirus_tac < 23.6, 1, 0)
  # 
  # # No aeromonas cutoff
  # 
  # # astrovirus cutoff 25.2
  # abcd_data$astrovirus_tac <- 35 - 3.322*abcd_data$astrovirus_new
  # abcd_data$astrovirus_likely <- ifelse(abcd_data$astrovirus_tac < 25.2,1, 0)
  # 
  # # campy jejuni coli cutoff 16.3
  # abcd_data$campylobacter_tac <- abcd_data$campylobacter_jejuni_coli
  # abcd_data$campylobacter_likely <- ifelse(abcd_data$campylobacter_tac < 16.3, 1, 0)
  # 
  # # cryptosporidium cutoff 25.6
  # abcd_data$cryptosporidium_tac <- 35 - 3.322*abcd_data$cryptosporidium_new
  # abcd_data$cryptosporidium_likely <- ifelse(abcd_data$cryptosporidium_tac < 25.6, 1, 0)
  # 
  # # cyclospora cutoff 33.8
  # abcd_data$cyclospora_tac <- abcd_data$cyclospora
  # abcd_data$cyclospora_likely <- ifelse(abcd_data$cyclospora_tac < 33.8, 1, 0)
  # 
  # # e hist cutoff 26.6
  # abcd_data$e_histolytica_tac <- abcd_data$e_histolytica
  # abcd_data$e_hist_likely <- ifelse(abcd_data$e_histolytica_tac < 26.6, 1, 0)
  # 
  # # isospora cutoff 32.6
  # abcd_data$isospora_tac <- abcd_data$isospora
  # abcd_data$isospora_likely <- ifelse(abcd_data$isospora < 32.6, 1, 0)
  # 
  # # norovirus cutoff 25.7
  # abcd_data$norovirus_tac <- 35 - 3.322*abcd_data$norovirus_new
  # abcd_data$norovirus_likely <- ifelse(abcd_data$norovirus_tac < 25.7, 1, 0)
  # 
  # # rotavirus cutoff 32
  # abcd_data$rotavirus_tac <- 35 - 3.322*abcd_data$rotavirus_new
  # abcd_data$rotavirus_likely <- ifelse(abcd_data$rotavirus_tac < 32.0, 1, 0)
  # 
  # # salmonella cutoff 31.9
  # abcd_data$salmonella_tac <- 35 - 3.322*abcd_data$salmonella_new
  # abcd_data$salmonella_likely <- ifelse(abcd_data$salmonella_tac < 31.9, 1, 0)
  # 
  # # sapovirus cutoff 22.5
  # abcd_data$sapovirus_tac <- abcd_data$sapovirus
  # abcd_data$sapovirus_likely <- ifelse(abcd_data$sapovirus_tac < 22.5, 1, 0)
  # 
  # # st-etec cutoff 25.4
  # abcd_data$st_etec_tac <- 35 - 3.322*abcd_data$st_etec_new
  # abcd_data$st_etec_likely <- ifelse(abcd_data$st_etec_tac < 25.4, 1, 0)
  # 
  # #tepec cutoff 18.1
  # abcd_data$tepec_tac <- 35 - 3.322*abcd_data$tepec_new
  # abcd_data$tepec_likely <- ifelse(abcd_data$tepec_tac < 18.1, 1, 0)
  # 
  # # v cholerae cutoff 32.6
  # abcd_data$v_cholerae_tac <- 35 - 3.322*abcd_data$v_cholerae_new
  # abcd_data$v_cholerae_likely <- ifelse(abcd_data$v_cholerae_tac < 32.6, 1, 0)
  # 
  # # ETEC (general)
  # abcd_data$etec_tac <- abcd_data$etec
  # abcd_data$etec_bin <- ifelse(abcd_data$etec_tac < 35, 1, 0)
  # 
  # # campylobacter (campylobacter_pan) note we want this instead of jejuni and coli but cutoff seems to be for the other one??
  # abcd_data$campylobacter_tac <- abcd_data$campylobacter_pan
  # 
  # # E Bienusi
  # abcd_data$e_bieneusi_tac <- abcd_data$e_bieneusi
  # abcd_data$e_bieneusi_bin <- ifelse(abcd_data$e_bieneusi_tac < 35, 1, 0)
  # 
  # # giardia
  # abcd_data$giardia_tac <- abcd_data$giardia
  # abcd_data$giardia_bin <- ifelse(abcd_data$giardia_tac < 35, 1, 0)
  # 
  # #eaec
  # abcd_data$eaec_tac <- abcd_data$eaec
  # abcd_data$eaec_bin <- ifelse(abcd_data$eaec_tac < 35, 1, 0)
  # 
  # # MAKE TAC TRANSFORMED VARIABLES FOR 
  # # etec
  # abcd_data$eaec_new <- (35 - abcd_data$eaec) / 3.322
  # abcd_data$etec_new <- (35 - as.numeric(abcd_data$etec)) / 3.322
  # abcd_data$giardia_new <- (35 - as.numeric(abcd_data$giardia)) / 3.322
  # abcd_data$e_bieneusi_new <- (35 - as.numeric(abcd_data$e_bieneusi_tac)) / 3.322
  # abcd_data$c_jejuni_coli_new <- (35 - as.numeric(abcd_data$campylobacter_jejuni_coli)) / 3.322
  # 
  # # New 1/30/25
  # # Make no_etiology variable using the pathogens of interest we have likely etiology for
  # abcd_data$no_etiology <- ifelse(rowSums(abcd_data[,c("shigella_likely",
  #                                                       "adenovirus_likely",
  #                                                        "astrovirus_likely",
  #                                                        "campylobacter_likely",
  #                                                        "cryptosporidium_likely",
  #                                                        "cyclospora_likely",
  #                                                        "e_hist_likely",
  #                                                        "isospora_likely",
  #                                                        "norovirus_likely",
  #                                                        "rotavirus_likely",
  #                                                        "salmonella_likely",
  #                                                        "sapovirus_likely",
  #                                                        "st_etec_likely",
  #                                                        "tepec_likely",
  #                                                        "v_cholerae_likely")], na.rm = TRUE) == 0, 1, 0)
  
  # Missing likely for lt etec, campy, sapo, e bienusi, giardia, eaec
  
  # temp manual fix till talk w. liz + sara:
  # one missing shigella_eiec but has shigella_flex and sonnei undetected so add as 0
  abcd_data$shigella_new[is.na(abcd_data$shigella_new)] <- 0
  abcd_data$shigella_bin[is.na(abcd_data$shigella_bin)] <- 0
  abcd_data$shigella_likely[is.na(abcd_data$shigella_likely)] <- 0
  
  # Label variables for tables
  abcd_data <- abcd_data %>%
    set_variable_labels(agemchild = "Age (months)",
                        month_en = "Enrollment month",
                        site = "Site",
                        dy1_ant_sex = "Sex",
                        an_ses_quintile = "Socioeconomic status quintile",
                        avemuac = "Average mid-upper arm circumference at baseline",
                        lfazscore = "Length-for-age z-score at baseline",
                        an_tothhlt5 = "Number of household members under age 5",
                        an_grp_01 = "Azithromycin",
                        lazd90 = "Length-for-age z-score at day 90",
                        wlzd90 = "Weight-for-length z-score at day 90",
                        day2diar = "3+ watery stools on day 2",
                        day3diar = "3+ watery stools on day 3",
                        lazdiff = "Difference in length-for-age z-score enrollment to day 90",
                        stunting_bin = "Stunting (day 90 LAZ < -2)",
                        an_hosp90death = "Indicates hospitalisation or death by D90",
                        an_death90 = "Indicates death by D90",
                        an_hosp90_yn = "Indicates hospitalisation by D90",
                        cryptosporidium_likely = "Cryptosporidium attribution likely",
                        shigella_likely = "Shigella attribution likely",
                        adenovirus_likely = "Adenovirus attribution likely",
                        rotavirus_likely = "Rotavirus attribution likely",
                        st_etec_likely = "ST Etec attribution likely",
                        astrovirus_likely = "Astrovirus attribution likely",
                        norovirus_gii_likely = "Norovirus attribution GII likely",
                        tepec_likely = "T Epec attribution likely",
                        no_etiology = "Diarrhea not attributable to other pathogens of interest",
                        shigella_bin = "Shigella detected",
                        rotavirus_bin = "Rotavirus detected",
                        adenovirus_bin = "Adenovirus detected",
                        etec_bin = "ETEC detected",
                        cryptosporidium_bin = "Cryprosporidium detected",
                        astrovirus_bin = "Astrovirus detected",
                        norovirus_gii_bin = "Norovirus GII detected",
                        tepec_bin = "TEPEC detected",
                        campylobacter_bin = "Campylobacter detected",
                        sapovirus_bin = "Sapovirus detected",
                        e_bieneusi_bin = "E bienusi detected",
                        giardia_bin = "Giardia detected",
                        eaec_bin = "EAEC detected",
                        rotavirus_new = "Rotavirus TAC scaled",
                        norovirus_gii_new = "Norovirus GII TAC scaled",
                        adenovirus_new = "Adenovirus TAC scaled",
                        astrovirus_new = "Astrovirus TAC scaled",
                        sapovirus_new = "Sapovirus TAC scaled",
                        st_etec_new = "ST ETEC TAC scaled",
                        shigella_new = "Shigella TAC scaled",
                        campylobacter_new = "Campylobacter TAC scaled",
                        tepec_new = "tEPEC TAC scaled",
                        v_cholerae_new = "V Cholerae TAC scaled",
                        salmonella_new = "Salmonella TAC scaled", 
                        etec_new = "ETEC TAC scaled",
                        eaec_new = "EAEC TAC scaled",
                        e_bieneusi_new ="E Bieneusi TAC scaled",
                        c_jejuni_new = "C jejuni coli scaled",
                        giardia_new = "Giardia TAC scaled",
                        an_d90_timing = "Day of follow-up visit",
                        I_an_d90_timing = "Indicator day of follow-up visit not missing",
                        I_an_d90_timing_x_an_d90_timing = "Indicator day of follow-up visit not missing * day",
                        gems_msd = "GEMS-definition-like MSD (hosp within 10 days and/or severe dehydration")
  
  return(abcd_data)
  
}

abcd_data_full <- prep_abcd()
abcd_data <- abcd_data_full[which(!is.na(abcd_data_full$shigella_new)), , drop = FALSE] # note no missing TAC so identical, keeping for consistency with other scripts

saveRDS(abcd_data, here::here("data/abcd_data/abcd_data.Rds"))
saveRDS(abcd_data_full, here::here("data/abcd_data/abcd_data_full.Rds"))
