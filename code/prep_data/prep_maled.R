# ---------------------------------------------------------------
# Script to create cleaned MAL-ED datasets
# ---------------------------------------------------------------

here::i_am("code/prep_data/prep_maled.R")

library(tidyverse)
library(labelled)
library(haven)
library(corrr)
library(factoextra)

#---------------------------------------------------------------
# CASE-ONLY ANALYSIS
# --------------------------------------------------------------

#' Function to create clean MAL-ED case dataset
#' 
#' By default, this will return all diarrhea episodes (modified for co-etiology 
#' meta-analysis). For the original meta-analysis, we define case as TAC 
#' or culture attributable Shigella so subset the dataset from this function to those who
#' have Shigella TAC results.
prep_mal_ed<- function(){
  
  # Load all data 
  zscore_data <- read.csv(here::here("data/maled_data/raw_data/zscores_all.csv"))
  tac_data <- read.csv(here::here("data/maled_data/raw_data/Tac_Sep2018.csv"))
  maled_bl <- read.csv(here::here("data/maled_data/raw_data/maled_baseline_all.csv"))
  diarrhea_data <- read.csv(here::here("data/maled_data/raw_data/diarrhea.csv"))
  maled_full <- read.csv(here::here("data/maled_data/raw_data/maled_full.csv"))
  
  # for culture data
  micro_data <- read.csv(here::here("data/maled_data/raw_data/micro_x.csv"))
  
  # Get wami quintile by site based on full data (not just diarrhea episodees)
  maled_bl <- maled_bl %>%
    group_by(Country_ID) %>%
    mutate(wami_quintile = ntile(wamiimp, 5))
  
  # Subset to kids with stooltype == "D1" to get diarrhea episodes
  tac_data <- tac_data[tac_data$stooltype == "D1",]
  
  # Get date of sample collection (use as episode date)
  tac_data$date <- as.POSIXct(tac_data$date, format = "%m/%d/%Y", tz = "UTC")
  
  # Get attribution by AFE > 0.5
  tac_data$shigella_attributable <- ifelse(tac_data$shigella_eiec_afe > 0.5, 1, 0)
  tac_data$adenovirus_attributable <- ifelse(tac_data$adenovirus_40_41_afe > 0.5, 1, 0)
  tac_data$aeromonas_attributable <- ifelse(tac_data$aeromonas_afe > 0.5, 1, 0)
  tac_data$astro_attributable <- ifelse(tac_data$astrovirus_afe > 0.5, 1, 0)
  tac_data$campylobacter_jejuni_coli_attributable <- ifelse(tac_data$campylobacter_jejuni_coli_afe > 0.5, 1, 0)
  tac_data$crypto_attributable <- ifelse(tac_data$cryptosporidium_afe > 0.5, 1, 0)
  tac_data$cyclospora_attributable <- ifelse(tac_data$cyclospora_afe > 0.5, 1, 0)
  tac_data$e_histolytica_attributable <- ifelse(tac_data$e_histolytica_afe > 0.5, 1, 0)
  tac_data$isospora_attributable <- ifelse(tac_data$isospora_afe > 0.5, 1, 0)
  tac_data$noro_attributable <- ifelse(tac_data$norovirus_afe > 0.5, 1, 0)
  tac_data$noro_gii_attributable <- ifelse(tac_data$norovirus_gii_afe > 0.5, 1, 0)
  tac_data$rotavirus_attributable <- ifelse(tac_data$rotavirus_afe > 0.5, 1, 0)
  tac_data$salmonella_attributable <- ifelse(tac_data$salmonella_afe > 0.5, 1, 0)
  tac_data$sapo_attributable <- ifelse(tac_data$sapovirus_afe > 0.5, 1, 0)
  tac_data$st_etec_attributable <- ifelse(tac_data$ST_ETEC_afe > 0.5, 1, 0)
  tac_data$tepec_attributable <- ifelse(tac_data$tEPEC_afe > 0.5, 1, 0)
  tac_data$v_cholerae_attributable <- ifelse(tac_data$v_cholerae_afe > 0.5, 1, 0)
  
  tac_data$ETEC_attributable <- ifelse(tac_data$ETEC_afe > 0.5, 1, 0)
  tac_data$e_bieneusi_attributable <- ifelse(tac_data$e_bieneusi_afe > 0.5, 1, 0)
  tac_data$giardia_attributable <- ifelse(tac_data$giardia_afe > 0.5, 1, 0)
  tac_data$eaec_attributable <- ifelse(tac_data$EAEC_afe > 0.5, 1, 0)
  
  # Get re-scaled pathogen quantities
  tac_data$shigella_new <- (35 - tac_data$shigella_eiec) / 3.322
  tac_data$adenovirus_40_41_new <- (35 - tac_data$adenovirus_40_41) / 3.322
  tac_data$aeromonas_new <- (35 - tac_data$aeromonas) / 3.322
  tac_data$astrovirus_new <- (35 - tac_data$astrovirus) / 3.322
  tac_data$campylobacter_pan_new <- (35 - tac_data$campylobacter_pan) / 3.322
  tac_data$c_jejuni_coli_new <- (35 - tac_data$campylobacter_jejuni_coli) / 3.322
  tac_data$cryptosporidium_new <- (35 - tac_data$cryptosporidium) / 3.322
  tac_data$cyclospora_new <- (35 - tac_data$cyclospora) / 3.322
  tac_data$e_histolytica_new <- (35 - tac_data$e_histolytica) / 3.322
  tac_data$isospora_new <- (35 - tac_data$isospora) / 3.322
  tac_data$norovirus_new <- (35 - tac_data$norovirus) / 3.322
  tac_data$norovirus_gii_new <- (35 - tac_data$norovirus_gii) / 3.322
  tac_data$rotavirus_new <- (35 - tac_data$rotavirus) / 3.322
  tac_data$salmonella_new <- (35 - tac_data$salmonella) / 3.322
  tac_data$sapovirus_new <- (35 - tac_data$sapovirus) / 3.322
  tac_data$st_etec_new <- (35 - tac_data$ST_ETEC) / 3.322
  tac_data$tEPEC_new <- (35 - tac_data$tEPEC) / 3.322
  tac_data$v_cholerae_new <- (35 - tac_data$v_cholerae) / 3.322
  tac_data$ETEC_new <- (35 - tac_data$ETEC) / 3.322
  tac_data$e_bieneusi_new <- (35 - tac_data$e_bieneusi) / 3.322
  tac_data$eaec_new <- (35 - tac_data$EAEC) / 3.322
  tac_data$giardia_new <- (35 - tac_data$giardia) / 3.322
  
  
  # Get culture attributable Shigella
  # Subset to kids with stooltype == "D1" to get diarrhea episodes (by culture)
  micro_diar <- micro_data[micro_data$stooltype == "D1",]
  
  # Get date of stool sample (use as episode date)
  micro_diar$episode_date <- as.POSIXct(micro_diar$srfdate, format = "%d%b%y:%H:%M:%S", tz = "UTC")
  
  micro_diar <- micro_diar %>%
    select(pid,
           episode_date,
           shigella) %>%
    rename(date = 'episode_date',
           culture_shigella = 'shigella') %>%
    #deduplicate, if there is NA or 0/1 keep 0/1, if it is 0 and 1 keep 1
    mutate(
      culture_priority = case_when(
        is.na(culture_shigella) ~ -1,         # lowest priority
        culture_shigella == 0 ~ 1,            # medium
        culture_shigella == 1 ~ 2             # highest
      )
    ) %>%
    arrange(pid, date, desc(culture_priority)) %>%
    distinct(pid, date, .keep_all = TRUE) %>%
    select(-culture_priority)
  
  #join culture info 
  tac_data <- left_join(tac_data, micro_diar, by = c("pid", "date"))
  
  tac_data$shigella_attributable_tac <- tac_data$shigella_attributable
  
  # Drop rows without TAC results 
  # NO LONGER DOING THIS NOW CO-ETIOLOGY
  # tac_data <- tac_data[-which(is.na(tac_data$shigella_attributable_tac)),]
  
  tac_data$shigella_attributable <- ifelse(
    !is.na(tac_data$culture_shigella) & tac_data$culture_shigella == 1, 
    1, 
    tac_data$shigella_attributable_tac
  )
  
  tac_data$no_etiology <- ifelse(rowSums(tac_data[,c("shigella_attributable",
                                                     "adenovirus_attributable",
                                                     "aeromonas_attributable",
                                                     "astro_attributable",
                                                     "campylobacter_jejuni_coli_attributable",
                                                     "crypto_attributable",
                                                     "cyclospora_attributable",
                                                     "e_histolytica_attributable",
                                                     "isospora_attributable",
                                                     "noro_attributable",
                                                     "rotavirus_attributable",
                                                     "salmonella_attributable",
                                                     "sapo_attributable",
                                                     "st_etec_attributable",
                                                     "tepec_attributable",
                                                     "v_cholerae_attributable")], na.rm = TRUE) == 0, 1, 0)
  
  # for other pathogens, use any tac < 35
  tac_data$rota_detect <- ifelse(tac_data$rotavirus < 35, 1, 0)
  tac_data$adeno_detect <- ifelse(tac_data$adenovirus_40_41 < 35, 1, 0)
  tac_data$ETEC_detect <- ifelse(tac_data$ETEC < 35, 1, 0)
  tac_data$crypto_detect <- ifelse(tac_data$cryptosporidium < 35,1,0)
  tac_data$astro_detect <- ifelse(tac_data$astrovirus < 35, 1, 0)
  tac_data$noro_detect <- ifelse(tac_data$norovirus < 35, 1, 0)
  tac_data$tepec_detect <- ifelse(tac_data$tEPEC < 35, 1, 0)
  tac_data$campy_detect <- ifelse(tac_data$campylobacter_pan < 35, 1, 0)
  tac_data$sapo_detect <- ifelse(tac_data$sapovirus < 35, 1, 0)
  tac_data$e_bieneusi_detect <- ifelse(tac_data$e_bieneusi < 35, 1, 0)
  tac_data$giardia_detect <- ifelse(tac_data$giardia < 35, 1, 0)
  tac_data$EAEC_detect <- ifelse(tac_data$EAEC < 35, 1, 0)
  
  # Get initial abx treatment variables
  tac_data$any_abx <- tac_data$abxtrt 
  tac_data$who_abx <- ifelse(tac_data$macrotrt == 1 | tac_data$fluorotrt == 1, 1, 0)
  
  tac_data$maybe_eff_abx <- ifelse(tac_data$cephalotrt == 1 | tac_data$sulfontrt == 1 | tac_data$tetratrt == 1 | 
                                     tac_data$othertrt == 1, 1, 0)
  
  tac_data$ineff_abx <- ifelse(tac_data$peniciltrt == 1 |
                                 tac_data$metrontrt == 1 |
                                 tac_data$unknowtrt == 1, 1, 0)
  
  tac_data$no_abx <- ifelse(tac_data$who_abx == 0 & tac_data$maybe_eff_abx == 0 & tac_data$ineff_abx == 0, 1, 0)
  
  tac_data$ineff_abx <- ifelse(tac_data$ineff_abx == 1 & (tac_data$who_abx == 1 | tac_data$maybe_eff_abx == 1), 0, tac_data$ineff_abx)
  tac_data$maybe_eff_abx <- ifelse(tac_data$maybe_eff_abx == 1 & tac_data$who_abx == 1, 0, tac_data$maybe_eff_abx)
  
  tac_data$all_abx <- ifelse(tac_data$no_abx == 1 | tac_data$ineff_abx == 1, 0,
                             ifelse(tac_data$maybe_eff_abx == 1, 1, 2))
  
  tac_data$all_abx <- factor(tac_data$all_abx, levels = 0:2, labels = c("No or ineffective antibiotics",
                                                                        "Possibly effective antibiotics",
                                                                        "Guideline recommended antibiotics"))
  
  # Select relevant variables from TAC dataset
  tac_data <- tac_data %>%
    select(pid,
           sid,
           country_id,
           date,
           agedays,
           shigella_attributable,
           shigella_attributable_tac,
           culture_shigella,
           adenovirus_attributable,
           aeromonas_attributable,
           astro_attributable,
           campylobacter_jejuni_coli_attributable,
           crypto_attributable,
           cyclospora_attributable,
           e_histolytica_attributable,
           isospora_attributable,
           noro_attributable,
           noro_gii_attributable,
           rotavirus_attributable,
           salmonella_attributable,
           sapo_attributable,
           st_etec_attributable,
           tepec_attributable,
           v_cholerae_attributable,
           ETEC_attributable,
           e_bieneusi_attributable,
           giardia_attributable,
           eaec_attributable,
           no_etiology,
           rota_detect,
           adeno_detect,
           ETEC_detect,
           crypto_detect,
           astro_detect,
           noro_detect,
           tepec_detect,
           campy_detect,
           sapo_detect,
           e_bieneusi_detect,
           giardia_detect,
           EAEC_detect,
           shigella_new,
           adenovirus_40_41_new,
           aeromonas_new,
           astrovirus_new,
           campylobacter_pan_new,
           c_jejuni_coli_new,
           cryptosporidium_new,
           cyclospora_new,
           e_histolytica_new,
           isospora_new,
           norovirus_new,
           norovirus_gii_new,
           rotavirus_new,
           salmonella_new,
           sapovirus_new,
           st_etec_new,
           tEPEC_new,
           v_cholerae_new,
           ETEC_new,
           e_bieneusi_new,
           eaec_new,
           giardia_new,
           any_abx,
           who_abx,
           maybe_eff_abx, 
           ineff_abx,
           no_abx,
           all_abx,
           prop_ebf30)
  
  # define daily who/any abx variables
  maled_full$date <- as.POSIXct(maled_full$date, format = "%d%b%Y", tz = "UTC")
  
  maled_full$who_abx <- ifelse(maled_full$safmacrolide == 1 | maled_full$saffluoro == 1, 1, 0)
  
  maled_full$maybe_eff_abx <- ifelse(maled_full$safcephalo == 1 | maled_full$safsulfon == 1 |
                                       maled_full$safother == 1, 1, 0)
  
  maled_full$ineff_abx <- ifelse(maled_full$safpenicillin == 1 |
                                   maled_full$safmetron == 1 | 
                                   maled_full$saftetracycl == 1 | 
                                   maled_full$safunknown == 1, 1, 0)
  
  maled_full$no_abx <- ifelse(maled_full$who_abx == 0 & maled_full$maybe_eff_abx == 0 & maled_full$ineff_abx == 0 ,1, 0)
  
  maled_full$ineff_abx <- ifelse(maled_full$ineff_abx == 1 & (maled_full$who_abx == 1 | maled_full$maybe_eff_abx ==1), 0, maled_full$ineff_abx)
  maled_full$maybe_eff_abx <- ifelse(maled_full$who_abx == 1 & maled_full$maybe_eff_abx == 1, 0, maled_full$maybe_eff_abx)
  
  maled_full$any_abx <- ifelse(maled_full$safpenicillin == 1 |
                                 maled_full$safcephalo == 1 | 
                                 maled_full$safsulfon == 1 |
                                 maled_full$safmacrolide == 1 |
                                 maled_full$saftetracycl == 1 |
                                 maled_full$saffluoro == 1 |
                                 maled_full$safunknown == 1 | 
                                 maled_full$safmetron == 1 | 
                                 maled_full$safother == 1, 1, 0)
  
  maled_full$all_abx <- ifelse(maled_full$no_abx == 1 | maled_full$ineff_abx == 1, 0,
                               ifelse(maled_full$maybe_eff_abx == 1, 1, 2))
  
  maled_full$all_abx <- factor(maled_full$all_abx, levels = 0:2, labels = c("No or ineffective antibiotics",
                                                                            "Possibly effective antibiotics",
                                                                            "Guideline recommended antibiotics"))
  
  # remake severity variables such that it's severity before antibiotics
  severity_df <- lapply(1:nrow(tac_data), function(i){
    row <- tac_data[i, ]
    
    episode_info <- maled_full[which(maled_full$stooldiaage == row$agedays &
                                       maled_full$Pid == row$pid),]
    
    if(nrow(episode_info) == 0){
      return(data.frame(maxb = NA,
                        fever = NA,
                        fever_days = NA,
                        maxls = NA,
                        sumvom = NA,
                        maxdehyd = NA,
                        alri = NA,
                        safcough = NA,
                        safshb = NA,
                        fstab = NA,
                        duration_pre_abx = 999)) # flag to remove row
    }
    
    # If no antibiotics, return overall info
    if(sum(episode_info$who_abx) == 0 & sum(episode_info$any_abx) == 0 & sum(episode_info$maybe_eff_abx) == 0){
      return(data.frame(maxb = episode_info$maxb[1],
                        fever = episode_info$fever[1],
                        fever_days = sum(episode_info$saffev, na.rm = TRUE),
                        maxls = episode_info$maxls[1],
                        sumvom = episode_info$sumvom[1],
                        maxdehyd = episode_info$maxdehyd[1],
                        alri = max(episode_info$alri),
                        safcough = max(episode_info$safcough),
                        safshb = max(episode_info$safshb),
                        fstab = max(episode_info$fstab),
                        duration_pre_abx = nrow(episode_info))) # returning length of episode 
    } else{
      # duration of episode prior to and including day they got antibiotics 
      
      # would like to specify for type coded in tac_data, for example if someone
      # gets ineffective abx then WHO abx later in episode, count prior to WHO abx
      # (added & all_abx == row$all_abx)
      # but there are some instances where safmacrolide (daily macrolide var) == 0 for whole episode 
      # even though macrotrt (episode treated with macrolides) == 1
      # so leaving as min(age[any_abx == 1]) for now
      
      pre_abx <- episode_info %>%
        mutate(first_abx = min(age[any_abx == 1])) %>% # & all_abx == row$all_abx])) %>%
        filter(age <= first_abx)
      
      # duration of episode prior to and including day they got antibiotics = nrow(pre_abx)
      
      return(data.frame(maxb = max(pre_abx$safblood, na.rm = TRUE),
                        fever = max(pre_abx$saffev, na.rm = TRUE), 
                        fever_days = sum(pre_abx$saffev, na.rm = TRUE),
                        maxls = max(pre_abx$safnumls, na.rm = TRUE),
                        sumvom = sum(pre_abx$safvom, na.rm = TRUE),
                        maxdehyd = max(pre_abx$safdehyd, na.rm = TRUE),
                        alri = max(pre_abx$alri, na.rm = TRUE),
                        safcough = max(pre_abx$safcough, na.rm = TRUE),
                        safshb = max(pre_abx$safshb, na.rm = TRUE),
                        fstab = max(episode_info$fstab, na.rm = TRUE),
                        duration_pre_abx = nrow(pre_abx)))
    }
    
  })
  
  severity_df <- do.call(rbind, severity_df)
  # severity_df[severity_df == -Inf] <- NA
  
  tac_data <- cbind(tac_data, severity_df)
  
  # Remove 999s (not applicable now with different merge)
  # tac_data <- tac_data[-which(tac_data$duration_pre_abx == 999),]
  
  # If any_abx = 1 and fstab = 0, received abx before episode began
  # Mark duration_pre_abx = 0
  tac_data$duration_pre_abx <- ifelse(tac_data$any_abx == 1 & tac_data$fstab == 0, 0, tac_data$duration_pre_abx)
  
  # Get dates of z-score measurements
  zscore_data$date <- strptime(zscore_data$date, format = "%d%b%Y", tz = "UTC")
  
  # Get baseline growth (HAZ at or within one month before episode date) & month 3 growth (HAZ closest to 3 months after episode +- 45 days)
  baseline_and_month3_df <- lapply(1:nrow(tac_data), function(i, zscore_data, tac_data){
    
    x <- tac_data[i,]
    
    # Subset zscore_data for the same participant
    sub_zscore <- zscore_data[zscore_data$pid == x$pid, ]
    
    # get all dates prior to episode date AND within 75 days
    baseline_candidates <- sub_zscore[
      sub_zscore$date <= x$date &
        abs(sub_zscore$date - x$date) <= 75,
    ]
    
    if (nrow(baseline_candidates) > 0) {
      
      # order by closeness to episode date
      baseline_candidates <- baseline_candidates[
        order(abs(baseline_candidates$date - x$date)),
      ]
      
      # initialize as NA
      baseline_date   <- NA
      baseline_haz    <- NA
      baseline_waz    <- NA
      baseline_whz    <- NA
      baseline_weight <- NA
      baseline_length <- NA
      
      # loop through candidates and take first non-missing HAZ
      for (j in 1:nrow(baseline_candidates)) {
        
        row <- baseline_candidates[j, ]
        
        haz <- row$haz
        if (is.na(haz)) haz <- row$zhei
        if (is.na(haz)) haz <- row$zheiorig
        
        if (!is.na(haz)) {
          
          baseline_date   <- row$date
          baseline_haz    <- haz
          baseline_waz    <- row$zwei
          baseline_whz    <- row$whz
          baseline_weight <- row$weight
          baseline_length <- row$length
          
          # WHZ fallback
          if (is.na(baseline_whz)) baseline_whz <- row$zwfl
          if (is.na(baseline_whz)) baseline_whz <- row$zwflorig
          
          break
        }
      }
      
    } else {
      baseline_date <- NA
      baseline_haz <- NA
      baseline_waz <- NA
      baseline_whz <- NA
      baseline_weight <- NA
      baseline_length <- NA
    }
    
    
    # # get all dates that are prior to episode date
    # baseline_dates <- sub_zscore$date[sub_zscore$date <= x$date]
    # 
    # # get baseline date closest to episode date and corresponding HAZ
    # # also add zwfl (weight for length zscore) or whz, use whz default
    # # zwei (weight for age zscore)
    # # weight
    # # length
    # if (length(baseline_dates) > 0) {
    #   baseline_date <- baseline_dates[which.min(abs(baseline_dates - x$date))]
    #   baseline_haz <- sub_zscore$haz[sub_zscore$date == baseline_date] #check-- looks like zlen (derived within tac codebook) == baseline_date so yay
    #   
    #   # NEW for describing tanzania
    #   baseline_waz <- sub_zscore$zwei[sub_zscore$date == baseline_date]
    #   baseline_whz <- sub_zscore$whz[sub_zscore$date == baseline_date]
    #   baseline_weight <- sub_zscore$weight[sub_zscore$date == baseline_date]
    #   baseline_length <- sub_zscore$length[sub_zscore$date == baseline_date]
    #   
    #   # If baseline_haz missing using HAZ, try using zhei, zheiorig
    #   if(is.na(baseline_haz)){
    #     baseline_haz <- sub_zscore$zhei[sub_zscore$date == baseline_date]
    #   }
    #   if(is.na(baseline_haz)){
    #     baseline_haz <- sub_zscore$zheiorig[sub_zscore$date == baseline_date]
    #   }
    #   
    #   # if baseline_whz missing using WHZ, try using zwfl, zwflorig
    #   if(is.na(baseline_whz)){
    #     baseline_whz <- sub_zscore$zwfl[sub_zscore$date == baseline_date]
    #   }
    #   if(is.na(baseline_whz)){
    #     baseline_whz <- sub_zscore$zwflorig[sub_zscore$date == baseline_date]
    #   }
    #   
    #   # Check to make sure date is within 75 days of baseline measurement
    #   if(is.na(baseline_date) | abs(x$date - baseline_date) > 75){
    #     baseline_date <- NA
    #     baseline_haz <- NA
    #   }
    #   
    #   # QUESTION if baseline HAZ missing in measurement closest to episode but there's another one within 75 day window that is present,
    #   # should i use that?
    #   
    # } else {
    #   baseline_date <- NA
    #   baseline_haz <- NA  # No baseline date found
    # }
    
    # get month3 date closest to 90 days post episode and corresponding HAZ
    target_month3_date <- x$date + days(90)
    sub_zscore <- sub_zscore[sub_zscore$date > baseline_date,]
    
    if(nrow(sub_zscore) == 0){
      month3_date <- NA
      month3_haz <- NA
    } else{
      month3_date <- sub_zscore$date[which.min(abs(sub_zscore$date - target_month3_date))]
      month3_haz <- sub_zscore$haz[sub_zscore$date == month3_date]
      
      # NEW for describing tanzania
      month3_waz <- sub_zscore$zwei[sub_zscore$date == month3_date]
      month3_whz <- sub_zscore$whz[sub_zscore$date == month3_date]
      month3_weight <- sub_zscore$weight[sub_zscore$date == month3_date]
      month3_length <- sub_zscore$length[sub_zscore$date == month3_date]
      
      # Check to make sure date is within 45 days of 3mo followup (updated from 75 days)
      if(length(month3_date) == 0 || abs(month3_date - target_month3_date) > 45){
        month3_date <- NA
        month3_haz <- NA
        
        month3_waz <- NA
        month3_whz <- NA
        month3_weight <- NA
        month3_length <- NA
      }
      
    }
    
    return(data.frame(baseline_haz = baseline_haz,
                      baseline_date = baseline_date,
                      baseline_waz = baseline_waz,
                      baseline_whz = baseline_whz,
                      baseline_weight = baseline_weight,
                      baseline_length = baseline_length,
                      month3_haz = month3_haz,
                      month3_date = month3_date,
                      month3_waz = month3_waz, 
                      month3_whz = month3_whz,
                      month3_weight = month3_weight, 
                      month3_length = month3_length ))
  }, zscore_data = zscore_data, tac_data = tac_data)
  
  baseline_and_month3_df <- do.call(rbind, baseline_and_month3_df) 
  final_df <- cbind(tac_data, baseline_and_month3_df)
  
  # rename covariates
  final_df <- final_df %>%
    rename(
      "episode_date" = date,
      "dysentery" = maxb,
      "lsstools" = maxls,
      "dehyd" = maxdehyd,
      "daysvomit" = sumvom,
      "cough" = safcough,
      "shortbreath" = safshb)
  
  # select covariates from bl data
  maled_bl <- maled_bl %>%
    select(Pid,
           CAFSEX,
           #mated,
           ageexbfimp, 
           Country_ID,
           incomeabovemed, #note not seeing this in the dictionary
           edimp,          #continuous maternal education
           incomemean,     #income? not in dictionary but liz said to use
           wamiimp,        #continuous version of SES score
           wami_quintile,
           drinkimp,
           sanitimp) %>%    # WAMI quintile by site
    rename("pid" = Pid,
           "sex" = CAFSEX,
           #"mated_bin" = mated,
           "site" = Country_ID,
           "maxagebf" = ageexbfimp,
           "mated_cont" = edimp,
           "income" = incomemean,
           "ses_wami" = wamiimp)
  
  maled_bl$mated_bin <- ifelse(maled_bl$mated_cont >= 6, 1, 0)
  
  maled_bl$wami_quintile <- factor(maled_bl$wami_quintile, levels = 1:5, labels = c("1st quintile of SES",
                                                                                    "2nd quintile of SES",
                                                                                    "3rd quintile of SES",
                                                                                    "4th quintile of SES",
                                                                                    "5th quintile of SES"))
  
  maled_bl$site <- factor(maled_bl$site, 
                          levels = c("BGD",
                                     "BRF",
                                     "INV",
                                     "NEB",
                                     "PEL",
                                     "PKN",
                                     "SAV",
                                     "TZH"),
                          labels = c("Bangladesh",
                                     "Brazil",
                                     "India",
                                     "Nepal",
                                     "Peru",
                                     "Pakistan",
                                     "South Africa",
                                     "Tanzania"))
  maled_bl$sex <- factor(maled_bl$sex, levels = c(1,2), labels = c("male", "female"))
  maled_bl$incomeabovemed <- factor(maled_bl$incomeabovemed, levels = c(0,1), labels = c("Income below country median",
                                                                                         "Income above country median"))
  
  # join covariates into final_df
  final_df <- left_join(final_df, maled_bl, by = "pid")
  
  # get rid of pakistan
  # Get rid of Pakistan and drop unused factor levels
  final_df <- final_df[which(final_df$site != "Pakistan"),]
  final_df$site <- droplevels(final_df$site)
  
  final_df$followup_days <- as.numeric(difftime(final_df$month3_date,final_df$baseline_date , units = "days"))
  final_df$I_followup_days <- ifelse(is.na(final_df$followup_days), 0, 1)
  final_df$I_followup_days_x_followup_days <- ifelse(is.na(final_df$followup_days), 0, final_df$followup_days)
  
  final_df$agemonths <- round(final_df$agedays / 30.44, 1)
  
  final_df$dehyd <- factor(final_df$dehyd, levels = c(0,1,2), labels = c("None", "Some dehydration", "Severe dehydration"))
  
  # Get rid of extreme HAZ observations
  final_df$month3_haz <- ifelse(final_df$month3_haz < -6 | final_df$month3_haz > 6, NA, final_df$month3_haz)
  final_df$baseline_haz <- ifelse(final_df$baseline_haz < -6 | final_df$baseline_haz > 6, NA, final_df$baseline_haz)
  
  final_df$hazdiff <- final_df$month3_haz - final_df$baseline_haz
  final_df$wazdiff <- final_df$month3_waz - final_df$baseline_waz
  final_df$wlzdiff <- final_df$month3_whz - final_df$baseline_whz
  final_df$lendiff <- final_df$month3_length - final_df$baseline_length
  final_df$wtdiff <- final_df$month3_weight - final_df$baseline_weight
  
  # add child ID as pid for sake of bootstrap make sure grabbing all episodes?? 
  #final_df$child_id <- final_df$pid
  
  # Add in GEMS definition of MSD from diarrhea data
  
  sub_diarrhea_data <- diarrhea_data[diarrhea_data$Pid %in% final_df$pid,]
  names(sub_diarrhea_data)[names(sub_diarrhea_data) == 'Pid'] <- "pid"
  names(sub_diarrhea_data)[names(sub_diarrhea_data) == 'age'] <- "agedays"
  names(sub_diarrhea_data)[names(sub_diarrhea_data) == 'gemsdef'] <- "MSD"
  names(sub_diarrhea_data)[names(sub_diarrhea_data) == 'incidabxtrt'] <- "initiated_abx_during_episode"
  
  final_df <- left_join(final_df, 
                        sub_diarrhea_data[,c("pid", "agedays", "MSD", "initiated_abx_during_episode")], 
                        by = c("pid" = "pid", 
                               "agedays" = "agedays"))
  
  # remove people who had antibiotics but it was not during episode
  # No longer doing this?? 3/20/25
  # check w/ liz incidabxtrt vs fstab
  # added fstab, marking duration as 0 but leaving them in
  # final_df <- final_df[-which(final_df$initiated_abx_during_episode == 0 & (final_df$any_abx == 1 | final_df$who_abx == 1 | final_df$maybe_eff_abx == 1)),]
  
  # Use sample ID as first ID (same convention as case control data)
  final_df <- final_df %>%
    arrange(pid, agedays) %>%
    group_by(pid) %>%
    mutate(first_id = sid[1]) %>%
    mutate(child_id = sid)
  
  final_df <- final_df %>%
    select(pid,
           sid, 
           first_id,
           child_id,
           episode_date,
           agedays,
           agemonths,
           shigella_attributable,
           shigella_attributable_tac,
           culture_shigella,
           adenovirus_attributable,
           aeromonas_attributable,
           astro_attributable,
           campylobacter_jejuni_coli_attributable,
           crypto_attributable,
           cyclospora_attributable,
           e_histolytica_attributable,
           isospora_attributable,
           noro_attributable,
           noro_gii_attributable,
           rotavirus_attributable,
           salmonella_attributable,
           sapo_attributable,
           st_etec_attributable,
           tepec_attributable,
           v_cholerae_attributable,
           ETEC_attributable,
           e_bieneusi_attributable,
           giardia_attributable,
           eaec_attributable,
           shigella_new,
           adenovirus_40_41_new,
           aeromonas_new,
           astrovirus_new,
           campylobacter_pan_new,
           c_jejuni_coli_new,
           cryptosporidium_new,
           cyclospora_new,
           e_histolytica_new,
           isospora_new,
           norovirus_new,
           norovirus_gii_new,
           rotavirus_new,
           salmonella_new,
           sapovirus_new,
           st_etec_new,
           tEPEC_new,
           v_cholerae_new,
           ETEC_new,
           e_bieneusi_new,
           eaec_new,
           giardia_new,
           no_etiology,
           rota_detect,
           adeno_detect,
           ETEC_detect,
           crypto_detect,
           astro_detect,
           noro_detect,
           tepec_detect,
           campy_detect,
           sapo_detect,
           e_bieneusi_detect,
           giardia_detect,
           EAEC_detect,
           any_abx,
           who_abx,
           maybe_eff_abx,
           ineff_abx,
           no_abx,
           all_abx,
           duration_pre_abx,
           dysentery,
           fever,
           fever_days,
           dehyd,
           lsstools,
           daysvomit,
           cough,
           shortbreath,
           alri,
           income,
           incomeabovemed,
           mated_cont,
           mated_bin,
           ses_wami,
           wami_quintile,
           drinkimp,
           sanitimp,
           baseline_haz,
           baseline_date,
           month3_haz,
           month3_date,
           baseline_waz ,
           baseline_whz ,
           baseline_weight, 
           baseline_length ,
           month3_waz ,
           month3_whz,
           month3_weight,
           month3_length,
           sex, 
           maxagebf,
           prop_ebf30,
           site,
           followup_days,
           I_followup_days,
           I_followup_days_x_followup_days,
           hazdiff, 
           wazdiff,
           wlzdiff,
           lendiff,
           wtdiff,
           MSD) %>%
    set_variable_labels(pid = "Participant ID",
                        episode_date = "Date of diarrhea episode",
                        agedays = "Age at sample collection (days)",
                        agemonths = "Age at sample collection (months)",
                        shigella_attributable = "Shigella attributable (AFE > 0.5 or culture)",
                        shigella_attributable_tac = "Shigella attributable (AFE > 0.5)",
                        culture_shigella = "Culture Shigella positive",
                        any_abx = "Received any antibiotics",
                        who_abx = "Received WHO approved antibiotics",
                        maybe_eff_abx = "Recieved maybe effective antibiotics", 
                        ineff_abx = "Recieved ineffective antibiotics",
                        no_abx = "Did not receive antibiotics",
                        all_abx = "Antibiotic type received",
                        baseline_haz = "HAZ at baseline (before & closest to episode)",
                        baseline_date = "Date of baseline HAZ measurement",
                        month3_haz = "HAZ at three months (after & closest to 90 days post-episode)",
                        month3_date = "Date of month three HAZ measurement",
                        dysentery = "Dysentery",
                        lsstools = "Max number of loose stools during episode",
                        dehyd = "Maximum severity of dehydration during diarrhea episode",
                        fever = "Reported fever during episode",
                        fever_days = "Days reported fever during episode",
                        daysvomit = "Days vommitted during episode",
                        cough = "Maternal report of cough",
                        shortbreath = "Maternal report of shortness of breath",
                        alri = "ALRI definition met",
                        rotavirus_attributable = "Rotavirus attributable (AFE > 0.5)",
                        crypto_attributable = "Cryptosporidium attributable (AFE > 0.5)",
                        adenovirus_attributable = "Adenovirus attributable (AFE > 0.5)",
                        ETEC_attributable = "ETEC attributable (AFE >0.5)",
                        astro_attributable = "Astrovirus attributable (AFE > 0.5)",
                        noro_attributable = "Norovirus attributable (AFE > 0.5)",
                        tepec_attributable = "tEPEC attributable (AFE > 0.5)",
                        sapo_attributable = "Sapovirus attributable (AFE > 0.5)",
                        e_bieneusi_attributable = "E bieneusi attributable (AFE > 0.5)",
                        giardia_attributable = "Giardia attributable (AFE > 0.5)",
                        eaec_attributable = "EAEC attributable (AFE > 0.5)",
                        no_etiology = "No attributable etiology",
                        rota_detect = "Rotavirus detected",
                        adeno_detect = "Adenovirus detected",
                        ETEC_detect = "ETEC detected",
                        crypto_detect = "Cryptosporidium detected",
                        astro_detect = "Astrovirus detected",
                        noro_detect = "Norovirus detected",
                        tepec_detect = "tEPEC detected",
                        campy_detect = "Campylobacter detected",
                        sapo_detect = "Sapovirus detected",
                        e_bieneusi_detect = "E Bieneusi detected",
                        giardia_detect = "Giardia detected",
                        EAEC_detect = "EAEC detected",
                        sex = "Sex",
                        mated_cont = "Years of maternal education",
                        mated_bin = "Mother completed >=6 years of school",
                        ses_wami = "WAMI Socioeconomic Status Score",
                        wami_quintile = "WAMI quintile (by site)",
                        drinkimp = "Improved drinking water",
                        sanitimp = "Improved sanitation",
                        maxagebf = "Max age of breastfeeding (imputed mean for country if missing)", #note could only find imputed version, could remove imputed values if needed. also concerned this is > age at episode in many cases. prop var better
                        prop_ebf30 = "Proportion of days of exclusive breastfeeding of the 30 days prior to episode",
                        site = "Site",
                        incomeabovemed = "Income above country median",
                        income = "Mean income",
                        followup_days = "Days between baseline HAZ and month 3 HAZ measurement",
                        hazdiff = "Difference between month 3 and baseline HAZ",
                        MSD = "Moderate to severe diarrhea (by GEMS definition)")
  
  return(final_df)
  
}

maled_data_full <- prep_mal_ed()
maled_data <- maled_data_full[which(!is.na(maled_data_full$shigella_attributable_tac)), , drop = FALSE] # subset to shigella TAC available 
maled_data_MSD <- maled_data[which(maled_data$MSD == 1), , drop = FALSE] # NEW version also excluding LSD

saveRDS(maled_data_full, here::here("data/maled_data/maled_data_full.Rds"))
saveRDS(maled_data, here::here("data/maled_data/maled_data.Rds"))
saveRDS(maled_data_MSD, here::here("data/maled_data/maled_data_shig_MSD_only.Rds"))

# -------------------------------------------------------------------
# CASE-CONTROL ANALYSIS
# -------------------------------------------------------------------

# GEMS inclusion criteria for controls:
# 1) resides in deomographic surveillance system area
# 2) age +-2mo for cases 0-11 mo
# age +-4mo for cases 12-59 mo
# ^^may not exceed boundaries of case
# 3) sex
# 4) enrolled within 14 days of presentation of case
# 5) no diarrhea in previous 7 days

# For MAL-ED:
# for each shigella diarrhea stool (case):
#   pull non-diarrhea stool same site, age, sex, date +- 14days
#   check if had diarrhea in previous 7 days

prep_maled_case_control <- function(case_def = "tac_or_culture_shig_diar", max_controls = NA){
  
  # set seed for sampling controls if not using all
  if(!is.na(max_controls)) set.seed(12345)
  
  # Load all data 
  zscore_data <- read.csv(here::here("data/maled_data/raw_data/zscores_all.csv"))
  tac_data <- read.csv(here::here("data/maled_data/raw_data/Tac_Sep2018.csv"))
  maled_bl <- read.csv(here::here("data/maled_data/raw_data/maled_baseline_all.csv"))
  diarrhea_data <- read.csv(here::here("data/maled_data/raw_data/diarrhea.csv"))
  maled_full <- read.csv(here::here("data/maled_data/raw_data/maled_full.csv"))
  
  # for culture data
  micro_data <- read.csv(here::here("data/maled_data/raw_data/micro_x.csv"))
  
  # Get wami quintile by site based on full data (not just diarrhea episodees)
  maled_bl <- maled_bl %>%
    group_by(Country_ID) %>%
    mutate(wami_quintile = ntile(wamiimp, 5))
  
  # defile daily abx variables
  maled_full$who_abx <- ifelse(maled_full$safmacrolide == 1 | maled_full$saffluoro == 1, 1, 0)
  maled_full$maybe_eff_abx <- ifelse(maled_full$safcephalo == 1 | maled_full$safsulfon == 1 | maled_full$saftetracycl == 1 | 
                                       maled_full$safother == 1, 1, 0)
  
  maled_full$ineff_abx <- ifelse(maled_full$safpenicillin == 1 | maled_full$safmetron == 1 | maled_full$safunknown == 1, 1, 0)
  maled_full$no_abx <- ifelse(maled_full$who_abx == 0 & maled_full$maybe_eff_abx == 0 & maled_full$ineff_abx == 0, 1, 0)
  
  maled_full$ineff_abx <- ifelse(maled_full$ineff_abx == 1 & (maled_full$who_abx == 1 | maled_full$maybe_eff_abx == 1), 0, maled_full$ineff_abx)
  maled_full$maybe_eff_abx <- ifelse(maled_full$maybe_eff_abx == 1 & maled_full$who_abx == 1, 0, maled_full$maybe_eff_abx)
  
  maled_full$all_abx <- ifelse(maled_full$no_abx == 1 | maled_full$ineff_abx == 1, 0,
                               ifelse(maled_full$maybe_eff_abx == 1, 1, 2))
  
  maled_full$all_abx <- factor(maled_full$all_abx, levels = 0:2, labels = c("No or ineffective antibiotics",
                                                                            "Possibly effective antibiotics",
                                                                            "Guideline recommended antibiotics"))
  
  maled_full$any_abx <- ifelse(maled_full$safpenicillin == 1 |
                                 maled_full$safcephalo == 1 | 
                                 maled_full$safsulfon == 1 |
                                 maled_full$safmacrolide == 1 |
                                 maled_full$saftetracycl == 1 |
                                 maled_full$saffluoro == 1 |
                                 maled_full$safunknown == 1 | 
                                 maled_full$safmetron == 1 | 
                                 maled_full$safother == 1, 1, 0)
  
  #### Define Shigella cases (tac & culture, culture only)
  
  ## TAC attributable
  tac_data$tac_shigella_attributable <- ifelse(tac_data$shigella_eiec_afe > 0.5, 1, 0)
  
  ## Culture attributable
  
  # Get date of stool sample (use as episode date)
  micro_data$episode_date <- as.POSIXct(micro_data$srfdate, format = "%d%b%y:%H:%M:%S", tz = "UTC")
  
  micro_diar <- micro_data %>%
    select(pid,
           episode_date,
           shigella) %>%
    rename(date = 'episode_date',
           culture_shigella = 'shigella') %>%
    #deduplicate, if there is NA or 0/1 keep 0/1, if it is 0 and 1 keep 1
    mutate(
      culture_priority = case_when(
        is.na(culture_shigella) ~ -1,         # lowest priority
        culture_shigella == 0 ~ 1,            # medium
        culture_shigella == 1 ~ 2             # highest
      )
    ) %>%
    arrange(pid, date, desc(culture_priority)) %>%
    distinct(pid, date, .keep_all = TRUE) %>%
    select(-culture_priority)
  
  #join culture info 
  # Get date of sample collection (use as episode date)
  tac_data$date <- as.POSIXct(tac_data$date, format = "%m/%d/%Y", tz = "UTC")
  tac_data <- left_join(tac_data, micro_diar, by = c("pid", "date"))
  
  # Drop rows without TAC Shigella results 
  # NO LONGER DOING THIS FOR CO-ETIOLOGY META_ANALYSIS
  # tac_data <- tac_data[-which(is.na(tac_data$tac_shigella_attributable)),]
  
  tac_data$shigella_attributable <- ifelse(
    !is.na(tac_data$culture_shigella) & tac_data$culture_shigella == 1, 
    1, 
    tac_data$tac_shigella_attributable
  )
  
  # ----- Identify cases and controls -----
  
  if(case_def == "tac_or_culture_shig_diar"){
    
    # Case = diarrhea attributable to Shigella via tac or culture
    tac_data$case <- ifelse(tac_data$stooltype == "D1" & tac_data$shigella_attributable == 1, 1, 0)
    
    # Get rid of non-shigella diarrhea
    tac_data_shig_only <- tac_data[-which(tac_data$case == 0 & tac_data$stooltype == "D1"),]
    
    # Get cases
    case_data <- tac_data_shig_only[which(tac_data_shig_only$case == 1),]
    
    case_data$case_pid <- case_data$pid
    case_data$case_sid <- case_data$sid
    
    # Dataset of eligible controls
    all_controls <- tac_data_shig_only[tac_data_shig_only$stooltype == "M1",]
    
  } else if(case_def == "tac_shig_diar"){
    
    # Case = diarrhea attributable to Shigella via TAC
    tac_data$case <- ifelse(tac_data$stooltype == "D1" & tac_data$tac_shigella_attributable == 1, 1, 0)
    
    # Get rid of non-shigella diarrhea
    tac_data_shig_only <- tac_data[-which(tac_data$case == 0 & tac_data$stooltype == "D1"),]
    
    # Get cases
    case_data <- tac_data_shig_only[which(tac_data_shig_only$case == 1),]
    
    case_data$case_pid <- case_data$pid
    case_data$case_sid <- case_data$sid
    
    # Dataset of eligible controls
    all_controls <- tac_data_shig_only[tac_data_shig_only$stooltype == "M1",]
    
  }else if(case_def == "culture_shig_diar"){
    # Case = diarrhea attributable to Shigella via culture
    # AND non-NA tac results
    tac_data$case <- ifelse(tac_data$stooltype == "D1" & tac_data$culture_shigella == 1 & !is.na(tac_data$tac_shigella_attributable), 1, 0)
    tac_data_shig_only <- tac_data[-which(tac_data$case == 0 & tac_data$stooltype == "D1"),]
    
    # Get cases
    case_data <- tac_data_shig_only[which(tac_data_shig_only$case == 1),]
    
    case_data$case_pid <- case_data$pid
    case_data$case_sid <- case_data$sid
    
    # Dataset of eligible controls
    all_controls <- tac_data_shig_only[tac_data_shig_only$stooltype == "M1",]
    
  } else if(case_def == "all_diar"){
    # Case = diarrhea in general
    
    tac_data$case <- ifelse(tac_data$stooltype == "D1", 1, 0)
    
    # Get cases
    case_data <- tac_data[which(tac_data$case == 1),]
    
    case_data$case_pid <- case_data$pid
    case_data$case_sid <- case_data$sid
    
    # Dataset of eligible controls
    all_controls <- tac_data[tac_data$stooltype == "M1",]
  } else{
    stop("Invalid case definition")
  }
  
  # Make control data
  
  # person ID = pid
  # sample ID = sid
  
  matched_controls <- lapply(1:nrow(case_data), function(i, case_data, all_controls, tac_data, max_controls){
    row <- case_data[i,]
    
    # get age range to match depending on case age
    # age range = 0-11 months (1-364 days) --> +- 2mo (30.44*2 mo= 61 days)
    # age range = 12+ months (365 days +) --> +- 4mo (30.44*4 mo = 122 days)
    if(row$agedays < 365){
      min_age <- max(0, row$agedays - 61)
      max_age <- min(row$agedays + 61, 364)
    } else{
      min_age <- max(365, row$agedays - 122)
      max_age <- row$agedays + 122
    }
    
    # match sex, site, time, age, not their own control
    matching_controls <- all_controls %>%
      filter(cafsex == row$cafsex) %>%
      filter(country_id == row$country_id) %>%
      filter(date < (row$date + days(15)) & date > row$date - days(15)) %>%
      filter(agedays >= min_age & agedays <= max_age) %>%
      filter(pid != row$pid) %>%
      group_by(pid) %>%
      slice_max(order_by = date, n = 1) %>% # Keep the latest sample per individual
      ungroup()
    
    if(nrow(matching_controls) > 0){
      # for each matching control, make sure no diarrhea 7 days prior
      control_eligible <- rep(TRUE, nrow(matching_controls))
      for(j in 1:nrow(matching_controls)){
        control_row <- matching_controls[j,]
        tac_data_match <- tac_data %>%
          filter(pid == control_row$pid) %>%
          filter(date <= control_row$date & date > control_row$date - days(7))
        
        if(any(tac_data_match$stooltype == "D1")){
          control_eligible[j] <- FALSE
        } 
      }
      
      # eliminate ineligible controls
      matching_controls <- matching_controls[control_eligible,]
      
      # if nrow(matching_controls > max_controls), take max_controls num of controls
      if(!is.na(max_controls) & nrow(matching_controls) > max_controls){
        ctrl_samp <- sample(1:nrow(matching_controls), max_controls)
        matching_controls <- matching_controls[ctrl_samp,]
      }
      
      matching_controls$case_pid <- row$pid
      matching_controls$case_sid <- row$sid
      matching_controls$no_match <- FALSE
      
    } else{
      # No matching controls
      matching_controls[1,] <- NA
      matching_controls$case_pid <- row$pid
      matching_controls$case_sid <- row$sid
      matching_controls$no_match <- TRUE
    }
    
    return(matching_controls)
    
  }, case_data = case_data, all_controls = all_controls, tac_data = tac_data, max_controls = max_controls)
  
  matched_controls <- do.call(rbind, matched_controls)
  
  # identify cases with no matches and remove from data (rare, only occurring in all case data (not shig case))
  bad_case_sids <- matched_controls$case_sid[matched_controls$no_match == TRUE]
  all_tac_cc <- rbind(case_data, matched_controls[,colnames(matched_controls) != "no_match"])
  if(length(bad_case_sids) > 0){
    all_tac_cc <- all_tac_cc[-which(all_tac_cc$case_sid %in% bad_case_sids), ]
  }
  
  # Get attribution by AFE > 0.5
  all_tac_cc$adenovirus_attributable <- ifelse(all_tac_cc$adenovirus_40_41_afe > 0.5, 1, 0)
  all_tac_cc$aeromonas_attributable <- ifelse(all_tac_cc$aeromonas_afe > 0.5, 1, 0)
  all_tac_cc$astro_attributable <- ifelse(all_tac_cc$astrovirus_afe > 0.5, 1, 0)
  all_tac_cc$campylobacter_jejuni_coli_attributable <- ifelse(all_tac_cc$campylobacter_jejuni_coli_afe > 0.5, 1, 0)
  all_tac_cc$crypto_attributable <- ifelse(all_tac_cc$cryptosporidium_afe > 0.5, 1, 0)
  all_tac_cc$cyclospora_attributable <- ifelse(all_tac_cc$cyclospora_afe > 0.5, 1, 0)
  all_tac_cc$e_histolytica_attributable <- ifelse(all_tac_cc$e_histolytica_afe > 0.5, 1, 0)
  all_tac_cc$isospora_attributable <- ifelse(all_tac_cc$isospora_afe > 0.5, 1, 0)
  all_tac_cc$noro_attributable <- ifelse(all_tac_cc$norovirus_afe > 0.5, 1, 0)
  all_tac_cc$rotavirus_attributable <- ifelse(all_tac_cc$rotavirus_afe > 0.5, 1, 0)
  all_tac_cc$salmonella_attributable <- ifelse(all_tac_cc$salmonella_afe > 0.5, 1, 0)
  all_tac_cc$sapo_attributable <- ifelse(all_tac_cc$sapovirus_afe > 0.5, 1, 0)
  all_tac_cc$st_etec_attributable <- ifelse(all_tac_cc$ST_ETEC_afe > 0.5, 1, 0)
  all_tac_cc$tepec_attributable <- ifelse(all_tac_cc$tEPEC_afe > 0.5, 1, 0)
  all_tac_cc$v_cholerae_attributable <- ifelse(all_tac_cc$v_cholerae_afe > 0.5, 1, 0)
  
  all_tac_cc$ETEC_attributable <- ifelse(all_tac_cc$ETEC_afe > 0.5, 1, 0)
  all_tac_cc$e_bieneusi_attributable <- ifelse(all_tac_cc$e_bieneusi_afe > 0.5, 1, 0)
  all_tac_cc$giardia_attributable <- ifelse(all_tac_cc$giardia_afe > 0.5, 1, 0)
  all_tac_cc$eaec_attributable <- ifelse(all_tac_cc$EAEC_afe > 0.5, 1, 0)
  
  all_tac_cc$no_etiology <- ifelse(rowSums(all_tac_cc[,c("shigella_attributable",
                                                         "adenovirus_attributable",
                                                         "aeromonas_attributable",
                                                         "astro_attributable",
                                                         "campylobacter_jejuni_coli_attributable",
                                                         "crypto_attributable",
                                                         "cyclospora_attributable",
                                                         "e_histolytica_attributable",
                                                         "isospora_attributable",
                                                         "noro_attributable",
                                                         "rotavirus_attributable",
                                                         "salmonella_attributable",
                                                         "sapo_attributable",
                                                         "st_etec_attributable",
                                                         "tepec_attributable",
                                                         "v_cholerae_attributable")], na.rm = TRUE) == 0, 1, 0)
  
  # for other pathogens, use any tac < 35
  all_tac_cc$rota_detect <- ifelse(all_tac_cc$rotavirus < 35, 1, 0)
  all_tac_cc$adeno_detect <- ifelse(all_tac_cc$adenovirus_40_41 < 35, 1, 0)
  all_tac_cc$ETEC_detect <- ifelse(all_tac_cc$ETEC < 35, 1, 0)
  all_tac_cc$crypto_detect <- ifelse(all_tac_cc$cryptosporidium < 35,1,0)
  all_tac_cc$astro_detect <- ifelse(all_tac_cc$astrovirus < 35, 1, 0)
  all_tac_cc$noro_detect <- ifelse(all_tac_cc$norovirus < 35, 1, 0)
  all_tac_cc$tepec_detect <- ifelse(all_tac_cc$tEPEC < 35, 1, 0)
  all_tac_cc$campy_detect <- ifelse(all_tac_cc$campylobacter_pan < 35, 1, 0)
  all_tac_cc$sapo_detect <- ifelse(all_tac_cc$sapovirus < 35, 1, 0)
  all_tac_cc$e_bieneusi_detect <- ifelse(all_tac_cc$e_bieneusi < 35, 1, 0)
  all_tac_cc$giardia_detect <- ifelse(all_tac_cc$giardia < 35, 1, 0)
  all_tac_cc$EAEC_detect <- ifelse(all_tac_cc$EAEC < 35, 1, 0)
  
  
  # Get re-scaled pathogen quantities
  all_tac_cc$shigella_new <- (35 - all_tac_cc$shigella_eiec) / 3.322
  all_tac_cc$adenovirus_40_41_new <- (35 - all_tac_cc$adenovirus_40_41) / 3.322
  all_tac_cc$aeromonas_new <- (35 - all_tac_cc$aeromonas) / 3.322
  all_tac_cc$astrovirus_new <- (35 - all_tac_cc$astrovirus) / 3.322
  all_tac_cc$campylobacter_pan_new <- (35 - all_tac_cc$campylobacter_pan) / 3.322
  all_tac_cc$cryptosporidium_new <- (35 - all_tac_cc$cryptosporidium) / 3.322
  all_tac_cc$cyclospora_new <- (35 - all_tac_cc$cyclospora) / 3.322
  all_tac_cc$e_histolytica_new <- (35 - all_tac_cc$e_histolytica) / 3.322
  all_tac_cc$isospora_new <- (35 - all_tac_cc$isospora) / 3.322
  all_tac_cc$norovirus_new <- (35 - all_tac_cc$norovirus) / 3.322
  all_tac_cc$rotavirus_new <- (35 - all_tac_cc$rotavirus) / 3.322
  all_tac_cc$salmonella_new <- (35 - all_tac_cc$salmonella) / 3.322
  all_tac_cc$sapovirus_new <- (35 - all_tac_cc$sapovirus) / 3.322
  all_tac_cc$st_etec_new <- (35 - all_tac_cc$ST_ETEC) / 3.322
  all_tac_cc$tEPEC_new <- (35 - all_tac_cc$tEPEC) / 3.322
  all_tac_cc$v_cholerae_new <- (35 - all_tac_cc$v_cholerae) / 3.322
  all_tac_cc$ETEC_new <- (35 - all_tac_cc$ETEC) / 3.322
  all_tac_cc$e_bieneusi_new <- (35 - all_tac_cc$e_bieneusi) / 3.322
  all_tac_cc$eaec_new <- (35 - all_tac_cc$EAEC) / 3.322
  all_tac_cc$giardia_new <- (35 - all_tac_cc$giardia) / 3.322
  
  
  # Get initial abx treatment variables
  all_tac_cc$any_abx <- all_tac_cc$abxtrt 
  all_tac_cc$who_abx <- ifelse(all_tac_cc$macrotrt == 1 | all_tac_cc$fluorotrt == 1, 1, 0)
  
  all_tac_cc$maybe_eff_abx <- ifelse(all_tac_cc$cephalotrt == 1 | all_tac_cc$sulfontrt == 1 | all_tac_cc$tetratrt == 1 | 
                                       all_tac_cc$othertrt == 1, 1, 0)
  
  all_tac_cc$ineff_abx <- ifelse(all_tac_cc$peniciltrt == 1 |
                                   all_tac_cc$metrontrt == 1 |
                                   all_tac_cc$unknowtrt == 1, 1, 0)
  
  all_tac_cc$no_abx <- ifelse(all_tac_cc$who_abx == 0 & all_tac_cc$maybe_eff_abx == 0 & all_tac_cc$ineff_abx == 0, 1, 0)
  
  all_tac_cc$ineff_abx <- ifelse(all_tac_cc$ineff_abx == 1 & (all_tac_cc$who_abx == 1 | all_tac_cc$maybe_eff_abx == 1), 0, all_tac_cc$ineff_abx)
  all_tac_cc$maybe_eff_abx <- ifelse(all_tac_cc$maybe_eff_abx == 1 & all_tac_cc$who_abx == 1, 0, all_tac_cc$maybe_eff_abx)
  
  all_tac_cc$all_abx <- ifelse(all_tac_cc$no_abx == 1 | all_tac_cc$ineff_abx == 1, 0,
                               ifelse(all_tac_cc$maybe_eff_abx == 1, 1, 2))
  
  all_tac_cc$all_abx <- factor(all_tac_cc$all_abx, levels = 0:2, labels = c("No or ineffective antibiotics",
                                                                            "Possibly effective antibiotics",
                                                                            "Guideline recommended antibiotics"))
  
  # Drop controls with abx 0-15 days before sample (healthy controls only)
  # all_tac_cc <- all_tac_cc[-which(all_tac_cc$case == 0 & all_tac_cc$abx15 == 1),]
  
  # Select relevant variables from TAC dataset
  all_tac_cc <- all_tac_cc %>%
    select(pid,
           sid,
           case_pid,
           case_sid,
           case,
           country_id,
           date,
           agedays,
           shigella_attributable,
           tac_shigella_attributable,
           culture_shigella,
           adenovirus_attributable,
           aeromonas_attributable,
           astro_attributable,
           campylobacter_jejuni_coli_attributable,
           crypto_attributable,
           cyclospora_attributable,
           e_histolytica_attributable,
           isospora_attributable,
           noro_attributable,
           rotavirus_attributable,
           salmonella_attributable,
           sapo_attributable,
           st_etec_attributable,
           tepec_attributable,
           v_cholerae_attributable,
           ETEC_attributable,
           e_bieneusi_attributable,
           giardia_attributable,
           eaec_attributable,
           no_etiology,
           rota_detect,
           adeno_detect,
           ETEC_detect,
           crypto_detect,
           astro_detect,
           noro_detect,
           tepec_detect,
           campy_detect,
           sapo_detect,
           e_bieneusi_detect,
           giardia_detect,
           EAEC_detect,
           shigella_new,
           adenovirus_40_41_new,
           aeromonas_new,
           astrovirus_new,
           campylobacter_pan_new,
           cryptosporidium_new,
           cyclospora_new,
           e_histolytica_new,
           isospora_new,
           norovirus_new,
           rotavirus_new,
           salmonella_new,
           sapovirus_new,
           st_etec_new,
           tEPEC_new,
           v_cholerae_new,
           ETEC_new,
           e_bieneusi_new,
           eaec_new,
           giardia_new,
           # now also include quantities
           rotavirus,
           adenovirus_40_41,
           ETEC,
           cryptosporidium,
           astrovirus,
           norovirus,
           tEPEC,
           campylobacter_pan,
           sapovirus,
           e_bieneusi,
           giardia,
           EAEC,
           any_abx,
           who_abx,
           maybe_eff_abx, 
           ineff_abx,
           no_abx,
           all_abx,
           prop_ebf30) %>%
    rename("tac_rotavirus" = rotavirus,
           "tac_adenovirus" = adenovirus_40_41,
           "tac_etec" = ETEC,
           "tac_crypto" = cryptosporidium,
           "tac_astrovirus" = astrovirus,
           "tac_norovirus" = norovirus,
           "tac_tEPEC" = tEPEC,
           "tac_campylobacter_pan" = campylobacter_pan,
           "tac_sapovirus" = sapovirus,
           "tac_e_bieneusi" = e_bieneusi,
           "tac_giardia" = giardia,
           "tac_EAEC" = EAEC)
  
  maled_full$date <- as.POSIXct(maled_full$date, format = "%d%b%Y", tz = "UTC")
  
  severity_pre_abx <- function(i){
    row <- all_tac_cc[i, ]
    
    if(row$case == 1){
      # CASE
      episode_info <- maled_full[which(maled_full$stooldiaage == row$agedays &
                                         maled_full$Pid == row$pid),]
      
      if(nrow(episode_info) == 0){
        return(data.frame(maxb = NA,
                          fever = NA,
                          fever_days = NA,
                          maxls = NA,
                          sumvom = NA,
                          maxdehyd = NA,
                          alri = NA,
                          safcough = NA,
                          safshb = NA,
                          fstab = NA,
                          duration_pre_abx = 999)) # flag to remove row
      }
      
      # If no antibiotics, return overall info
      if(sum(episode_info$who_abx) == 0 & sum(episode_info$any_abx) == 0 & sum(episode_info$maybe_eff_abx) == 0){
        return(data.frame(maxb = episode_info$maxb[1],
                          fever = episode_info$fever[1],
                          fever_days = sum(episode_info$saffev, na.rm = TRUE),
                          maxls = episode_info$maxls[1],
                          sumvom = episode_info$sumvom[1],
                          maxdehyd = episode_info$maxdehyd[1],
                          alri = max(episode_info$alri),
                          safcough = max(episode_info$safcough),
                          safshb = max(episode_info$safshb),
                          fstab = max(episode_info$fstab),
                          duration_pre_abx = nrow(episode_info))) # returning length of episode 
      } else{
        
        # would like to specify for type coded in tac_data, for example if someone
        # gets ineffective abx then WHO abx later in episode, count prior to WHO abx
        # (added & all_abx == row$all_abx)
        # but there are some instances where safmacrolide (daily macrolide var) == 0 for whole episode 
        # even though macrotrt (episode treated with macrolides) == 1
        # so leaving as min(age[any_abx == 1]) for now
        
        pre_abx <- episode_info %>%
          mutate(first_abx = min(age[any_abx == 1])) %>% # & all_abx == row$all_abx])) %>%
          filter(age <= first_abx)
        
        return(data.frame(maxb = max(pre_abx$safblood, na.rm = TRUE),
                          fever = max(pre_abx$saffev, na.rm = TRUE), 
                          fever_days = sum(pre_abx$saffev, na.rm = TRUE),
                          maxls = max(pre_abx$safnumls, na.rm = TRUE),
                          sumvom = sum(pre_abx$safvom, na.rm = TRUE),
                          maxdehyd = max(pre_abx$safdehyd, na.rm = TRUE),
                          alri = max(pre_abx$alri, na.rm = TRUE),
                          safcough = max(pre_abx$safcough, na.rm = TRUE),
                          safshb = max(pre_abx$safshb, na.rm = TRUE),
                          fstab = max(episode_info$fstab, na.rm = TRUE),
                          duration_pre_abx = nrow(pre_abx)))
      }
    } else{
      # CONTROL
      
      # check to make sure not taking abx on sample date (could be on abx for something else)
      full_row <- maled_full[which(maled_full$Pid == row$pid & maled_full$date == row$date),]
      
      # full_row$any_abx == 1 if any of the saf daily abx vars == 1
      # nrow(full_row) == 0 if no entry in longitudinal data for that day
      if((nrow(full_row) == 0) || full_row$any_abx == 1) {
        # 999 to indicate drop row
        return(data.frame(maxb = 999,
                          fever = 999, 
                          fever_days = 999,
                          maxls = 999,
                          sumvom = 999,
                          maxdehyd = 999,
                          alri = 999,
                          safcough =999,
                          safshb = 999,
                          fstab = 999,
                          duration_pre_abx = 999))
      } else {
        # NA because not adjusting for severity in controls
        return(data.frame(maxb = NA,
                          fever = NA, 
                          fever_days = NA,
                          maxls = NA,
                          sumvom = NA,
                          maxdehyd = NA,
                          alri = NA,
                          safcough =NA,
                          safshb = NA,
                          fstab = NA,
                          duration_pre_abx = NA))
      }
      
      
    }
  }
  
  severity_df <- lapply(1:nrow(all_tac_cc), severity_pre_abx)
  severity_df <- do.call(rbind, severity_df)
  
  all_tac_cc <- cbind(all_tac_cc, severity_df)
  
  # Drop any controls taking abx on day of sample (indicated by 999 in severity cols)
  all_tac_cc <- all_tac_cc[-which(all_tac_cc$case == 0 & all_tac_cc$maxb == 999),]
  
  # drop diarrhea episodes coded with 999 in duration (not applicable now that switched severity merge)
  # all_tac_cc <- all_tac_cc[-which(all_tac_cc$case == 1 & all_tac_cc$duration_pre_abx == 999),]
  
  # If any_abx = 1 and fstab = 0, received abx before episode began
  # Mark duration_pre_abx = 0
  all_tac_cc$duration_pre_abx <- ifelse(all_tac_cc$any_abx == 1 & all_tac_cc$fstab == 0, 0, all_tac_cc$duration_pre_abx)
  
  # Get dates of z-score measurements
  zscore_data$date <- strptime(zscore_data$date, format = "%d%b%Y", tz = "UTC")
  
  # Get baseline growth (HAZ at or within one month before episode date) & month 3 growth (HAZ closest to 3 months after episode)
  baseline_and_month3_df <- lapply(1:nrow(all_tac_cc), function(i, zscore_data, tac_data){
    
    x <- tac_data[i,]
    
    # Subset zscore_data for the same participant
    sub_zscore <- zscore_data[zscore_data$pid == x$pid, ]
    
    # get all dates prior to episode date AND within 75 days
    baseline_candidates <- sub_zscore[
      sub_zscore$date <= x$date &
        abs(sub_zscore$date - x$date) <= 75,
    ]
    
    if (nrow(baseline_candidates) > 0) {
      
      # order by closeness to episode date
      baseline_candidates <- baseline_candidates[
        order(abs(baseline_candidates$date - x$date)),
      ]
      
      # initialize as NA
      baseline_date   <- NA
      baseline_haz    <- NA
      baseline_waz    <- NA
      baseline_whz    <- NA
      baseline_weight <- NA
      baseline_length <- NA
      
      # loop through candidates and take first non-missing HAZ
      for (j in 1:nrow(baseline_candidates)) {
        
        row <- baseline_candidates[j, ]
        
        haz <- row$haz
        if (is.na(haz)) haz <- row$zhei
        if (is.na(haz)) haz <- row$zheiorig
        
        if (!is.na(haz)) {
          
          baseline_date   <- row$date
          baseline_haz    <- haz
          baseline_waz    <- row$zwei
          baseline_whz    <- row$whz
          baseline_weight <- row$weight
          baseline_length <- row$length
          
          # WHZ fallback
          if (is.na(baseline_whz)) baseline_whz <- row$zwfl
          if (is.na(baseline_whz)) baseline_whz <- row$zwflorig
          
          break
        }
      }
      
    } else {
      baseline_date <- NA
      baseline_haz <- NA
      baseline_waz <- NA
      baseline_whz <- NA
      baseline_weight <- NA
      baseline_length <- NA
    }
    
    # get month3 date closest to 90 days post episode and corresponding HAZ
    target_month3_date <- x$date + days(90)
    sub_zscore <- sub_zscore[sub_zscore$date > baseline_date,]
    
    if(nrow(sub_zscore) == 0){
      month3_date <- NA
      month3_haz <- NA
      
      month3_waz <- NA
      month3_whz <- NA
      month3_weight <- NA
      month3_length <- NA
    } else{
      month3_date <- sub_zscore$date[which.min(abs(sub_zscore$date - target_month3_date))]
      month3_haz <- sub_zscore$haz[sub_zscore$date == month3_date]
      
      # NEW for describing tanzania
      month3_waz <- sub_zscore$zwei[sub_zscore$date == month3_date]
      month3_whz <- sub_zscore$whz[sub_zscore$date == month3_date]
      month3_weight <- sub_zscore$weight[sub_zscore$date == month3_date]
      month3_length <- sub_zscore$length[sub_zscore$date == month3_date]
      
      # Check to make sure date is within 45 days of 3mo followup (updated from 75 days)
      if(length(month3_date) == 0 || abs(month3_date - target_month3_date) > 45){
        month3_date <- NA
        month3_haz <- NA
        
        month3_waz <- NA
        month3_whz <- NA
        month3_weight <- NA
        month3_length <- NA
      }
      
    }
    
    return(data.frame(baseline_haz = baseline_haz,
                      baseline_date = baseline_date,
                      baseline_waz = baseline_waz,
                      baseline_whz = baseline_whz,
                      baseline_weight = baseline_weight,
                      baseline_length = baseline_length,
                      month3_haz = month3_haz,
                      month3_date = month3_date,
                      month3_waz = month3_waz, 
                      month3_whz = month3_whz,
                      month3_weight = month3_weight, 
                      month3_length = month3_length ))
    
  }, zscore_data = zscore_data, tac_data = all_tac_cc)
  
  
  baseline_and_month3_df <- do.call(rbind, baseline_and_month3_df) 
  final_df <- cbind(all_tac_cc, baseline_and_month3_df)
  
  # rename covariates
  final_df <- final_df %>%
    rename(
      "episode_date" = date,
      "dysentery" = maxb,
      "lsstools" = maxls,
      "dehyd" = maxdehyd,
      "daysvomit" = sumvom,
      "cough" = safcough,
      "shortbreath" = safshb)
  
  # select covariates from bl data
  maled_bl <- maled_bl %>%
    select(Pid,
           CAFSEX,
           #mated,
           ageexbfimp, 
           Country_ID,
           incomeabovemed, #note not seeing this in the dictionary
           edimp,          #continuous maternal education
           incomemean,     #income? not in dictionary but liz said to use
           wamiimp,        #continuous version of SES score
           wami_quintile,
           drinkimp,
           sanitimp) %>%    # WAMI quintile by site
    rename("pid" = Pid,
           "sex" = CAFSEX,
           #"mated_bin" = mated,
           "site" = Country_ID,
           "maxagebf" = ageexbfimp,
           "mated_cont" = edimp,
           "income" = incomemean,
           "ses_wami" = wamiimp)
  
  maled_bl$mated_bin <- ifelse(maled_bl$mated_cont >= 6, 1, 0)
  
  maled_bl$wami_quintile <- factor(maled_bl$wami_quintile, levels = 1:5, labels = c("1st quintile of SES",
                                                                                    "2nd quintile of SES",
                                                                                    "3rd quintile of SES",
                                                                                    "4th quintile of SES",
                                                                                    "5th quintile of SES"))
  
  maled_bl$site <- factor(maled_bl$site, 
                          levels = c("BGD",
                                     "BRF",
                                     "INV",
                                     "NEB",
                                     "PEL",
                                     "PKN",
                                     "SAV",
                                     "TZH"),
                          labels = c("Bangladesh",
                                     "Brazil",
                                     "India",
                                     "Nepal",
                                     "Peru",
                                     "Pakistan",
                                     "South Africa",
                                     "Tanzania"))
  maled_bl$sex <- factor(maled_bl$sex, levels = c(1,2), labels = c("male", "female"))
  maled_bl$incomeabovemed <- factor(maled_bl$incomeabovemed, levels = c(0,1), labels = c("Income below country median",
                                                                                         "Income above country median"))
  
  # join covariates into final_df
  final_df <- left_join(final_df, maled_bl, by = "pid")
  
  # get rid of pakistan
  # Get rid of Pakistan and drop unused factor levels
  final_df <- final_df[which(final_df$site != "Pakistan"),]
  final_df$site <- droplevels(final_df$site)
  
  final_df$followup_days <- as.numeric(difftime(final_df$month3_date,final_df$baseline_date , units = "days"))
  final_df$I_followup_days <- ifelse(is.na(final_df$followup_days), 0, 1)
  final_df$I_followup_days_x_followup_days <- ifelse(is.na(final_df$followup_days), 0, final_df$followup_days)
  
  final_df$agemonths <- round(final_df$agedays / 30.44, 1)
  
  final_df$dehyd <- factor(final_df$dehyd, levels = c(0,1,2), labels = c("None", "Some dehydration", "Severe dehydration"))
  
  # Get rid of extreme HAZ observations
  final_df$month3_haz <- ifelse(final_df$month3_haz < -6 | final_df$month3_haz > 6, NA, final_df$month3_haz)
  final_df$baseline_haz <- ifelse(final_df$baseline_haz < -6 | final_df$baseline_haz > 6, NA, final_df$baseline_haz)
  
  final_df$hazdiff <- final_df$month3_haz - final_df$baseline_haz
  final_df$wazdiff <- final_df$month3_waz - final_df$baseline_waz
  final_df$wlzdiff <- final_df$month3_whz - final_df$baseline_whz
  final_df$lendiff <- final_df$month3_length - final_df$baseline_length
  final_df$wtdiff <- final_df$month3_weight - final_df$baseline_weight
  
  # Add in GEMS definition of MSD from diarrhea data
  
  sub_diarrhea_data <- diarrhea_data[diarrhea_data$Pid %in% final_df$pid,]
  names(sub_diarrhea_data)[names(sub_diarrhea_data) == 'Pid'] <- "pid"
  names(sub_diarrhea_data)[names(sub_diarrhea_data) == 'age'] <- "agedays"
  names(sub_diarrhea_data)[names(sub_diarrhea_data) == 'gemsdef'] <- "MSD"
  names(sub_diarrhea_data)[names(sub_diarrhea_data) == 'incidabxtrt'] <- "initiated_abx_during_episode"
  
  final_df <- left_join(final_df, 
                        sub_diarrhea_data[,c("pid", "agedays", "MSD", "initiated_abx_during_episode")], 
                        by = c("pid" = "pid", 
                               "agedays" = "agedays"))
  
  final_df <- final_df %>%
    select(pid,
           sid,
           case_pid,
           case_sid,
           case,
           episode_date,
           agedays,
           agemonths,
           shigella_attributable,
           tac_shigella_attributable,
           culture_shigella,
           adenovirus_attributable,
           aeromonas_attributable,
           astro_attributable,
           campylobacter_jejuni_coli_attributable,
           crypto_attributable,
           cyclospora_attributable,
           e_histolytica_attributable,
           isospora_attributable,
           noro_attributable,
           rotavirus_attributable,
           salmonella_attributable,
           sapo_attributable,
           st_etec_attributable,
           tepec_attributable,
           v_cholerae_attributable,
           ETEC_attributable,
           e_bieneusi_attributable,
           giardia_attributable,
           eaec_attributable,
           no_etiology,
           rota_detect,
           adeno_detect,
           ETEC_detect,
           crypto_detect,
           astro_detect,
           noro_detect,
           tepec_detect,
           campy_detect,
           sapo_detect,
           e_bieneusi_detect,
           giardia_detect,
           EAEC_detect,
           tac_rotavirus,
           tac_adenovirus,
           tac_etec,
           tac_crypto,
           tac_astrovirus,
           tac_norovirus,
           tac_tEPEC,
           tac_campylobacter_pan,
           tac_sapovirus,
           tac_e_bieneusi,
           tac_giardia,
           tac_EAEC,
           shigella_new,
           adenovirus_40_41_new,
           aeromonas_new,
           astrovirus_new,
           campylobacter_pan_new,
           cryptosporidium_new,
           cyclospora_new,
           e_histolytica_new,
           isospora_new,
           norovirus_new,
           rotavirus_new,
           salmonella_new,
           sapovirus_new,
           st_etec_new,
           giardia_new,
           tEPEC_new,
           v_cholerae_new,
           ETEC_new,
           e_bieneusi_new,
           eaec_new,
           no_etiology,
           any_abx,
           who_abx,
           maybe_eff_abx,
           ineff_abx,
           no_abx,
           all_abx,
           duration_pre_abx,
           dysentery,
           fever,
           fever_days,
           dehyd,
           lsstools,
           daysvomit,
           cough,
           shortbreath,
           alri,
           income,
           incomeabovemed,
           mated_cont,
           mated_bin,
           ses_wami,
           wami_quintile,
           drinkimp,
           sanitimp,
           baseline_haz,
           baseline_date,
           month3_haz,
           month3_date,
           baseline_waz ,
           baseline_whz ,
           baseline_weight, 
           baseline_length ,
           month3_waz ,
           month3_whz,
           month3_weight,
           month3_length,
           sex, 
           maxagebf,
           prop_ebf30,
           site,
           followup_days,
           I_followup_days,
           I_followup_days_x_followup_days,
           hazdiff, 
           wazdiff,
           wlzdiff,
           lendiff,
           wtdiff,
           MSD) %>%
    set_variable_labels(pid = "Participant ID",
                        sid = "Sample ID",
                        case_pid = "Participant ID of case",
                        case_sid = "Sample ID of case",
                        case = "Case",
                        episode_date = "Date of diarrhea episode",
                        agedays = "Age at sample collection (days)",
                        agemonths = "Age at sample collection (months)",
                        shigella_attributable = "Shigella attributable (AFE > 0.5)",
                        any_abx = "Received any antibiotics",
                        who_abx = "Received WHO approved antibiotics",
                        maybe_eff_abx = "Recieved maybe effective antibiotics",
                        ineff_abx = "Recieved ineffective or no antibiotics",
                        baseline_haz = "HAZ at baseline (before & closest to episode)",
                        baseline_date = "Date of baseline HAZ measurement",
                        month3_haz = "HAZ at three months (after & closest to 90 days post-episode)",
                        month3_date = "Date of month three HAZ measurement",
                        dysentery = "Dysentery",
                        lsstools = "Max number of loose stools during episode",
                        dehyd = "Maximum severity of dehydration during diarrhea episode",
                        fever = "Reported fever during episode",
                        fever_days = "Days reported fever during episode",
                        daysvomit = "Days vommitted during episode",
                        cough = "Maternal report of cough",
                        shortbreath = "Maternal report of shortness of breath",
                        alri = "ALRI definition met",
                        rotavirus_attributable = "Rotavirus attributable (AFE > 0.5)",
                        crypto_attributable = "Cryptosporidium attributable (AFE > 0.5)",
                        adenovirus_attributable = "Adenovirus attributable (AFE > 0.5)",
                        ETEC_attributable = "ETEC attributable (AFE >0.5)",
                        astro_attributable = "Astrovirus attributable (AFE > 0.5)",
                        noro_attributable = "Norovirus attributable (AFE > 0.5)",
                        tepec_attributable = "tEPEC attributable (AFE > 0.5)",
                        sapo_attributable = "Sapovirus attributable (AFE > 0.5)",
                        e_bieneusi_attributable = "E bieneusi attributable (AFE > 0.5)",
                        giardia_attributable = "Giardia attributable (AFE > 0.5)",
                        eaec_attributable = "EAEC attributable (AFE > 0.5)",
                        no_etiology = "No other attributable etiology",
                        rota_detect = "Rotavirus detected",
                        adeno_detect = "Adenovirus detected",
                        ETEC_detect = "ETEC detected",
                        crypto_detect = "Cryptosporidium detected",
                        astro_detect = "Astrovirus detected",
                        noro_detect = "Norovirus detected",
                        tepec_detect = "tEPEC detected",
                        campy_detect = "Campylobacter detected",
                        sapo_detect = "Sapovirus detected",
                        e_bieneusi_detect = "E Bieneusi detected",
                        giardia_detect = "Giardia detected",
                        EAEC_detect = "EAEC detected",
                        sex = "Sex",
                        mated_cont = "Years of maternal education",
                        mated_bin = "Mother completed >=6 years of school",
                        ses_wami = "WAMI Socioeconomic Status Score",
                        wami_quintile = "WAMI quintile (by site)",
                        drinkimp = "Improved drinking water",
                        sanitimp = "Improved sanitation",
                        maxagebf = "Max age of breastfeeding (imputed mean for country if missing)", #note could only find imputed version, could remove imputed values if needed. also concerned this is > age at episode in many cases. prop var better
                        prop_ebf30 = "Proportion of days of exclusive breastfeeding of the 30 days prior to episode",
                        site = "Site",
                        incomeabovemed = "Income above country median",
                        income = "Mean income",
                        followup_days = "Days between baseline HAZ and month 3 HAZ measurement",
                        hazdiff = "Difference between month 3 and baseline HAZ",
                        MSD = "Moderate to severe diarrhea (by GEMS definition)")
  
  # Add same variables as VIDA/GEMS for bootstrap
  
  # first_id = associated with the child
  # case_id = associated with the case
  # child_id = associated with the episode 
  
  # when first sample is case then first = case = child
  final_df <- final_df %>%
    arrange(pid, agedays) %>%
    group_by(pid) %>%
    mutate(first_id = sid[1]) %>%
    mutate(case_id = case_sid,
           child_id = sid)
  
  return(final_df)
  
}

maled_case_control_data_culture_or_tac <- prep_maled_case_control(case_def = "tac_or_culture_shig_diar")
saveRDS(maled_case_control_data_culture_or_tac, here::here("data/maled_data/maled_case_control_tac_or_culture_shig.Rds"))

maled_case_control_data_tac <- prep_maled_case_control(case_def = "tac_shig_diar")
saveRDS(maled_case_control_data_tac, here::here("data/maled_data/maled_case_control_tac_shig.Rds"))

maled_case_control_data_culture <- prep_maled_case_control(case_def = "culture_shig_diar")
saveRDS(maled_case_control_data_culture, here::here("data/maled_data/maled_case_control_culture_shig.Rds"))

maled_case_control_data_all <- prep_maled_case_control(case_def = "all_diar", max_controls = 3)
saveRDS(maled_case_control_data_all, here::here("data/maled_data/maled_case_control_all.Rds"))

# ---------------------------------------------------------------------------

# Prep monthly growth effect data for longitudinal analysis (figure of the year)


prep_maled_case_control_monthx <- function(case_def = "tac_or_culture_shig_diar", max_controls = NA, month_x_days = 90){
  
  # set seed for sampling controls if not using all
  if(!is.na(max_controls)) set.seed(12345)
  
  # Load all data 
  zscore_data <- read.csv(here::here("data/maled_data/raw_data/zscores_all.csv"))
  tac_data <- read.csv(here::here("data/maled_data/raw_data/Tac_Sep2018.csv"))
  maled_bl <- read.csv(here::here("data/maled_data/raw_data/maled_baseline_all.csv"))
  diarrhea_data <- read.csv(here::here("data/maled_data/raw_data/diarrhea.csv"))
  maled_full <- read.csv(here::here("data/maled_data/raw_data/maled_full.csv"))
  
  # for culture data
  micro_data <- read.csv(here::here("data/maled_data/raw_data/micro_x.csv"))
  
  # Get wami quintile by site based on full data (not just diarrhea episodees)
  maled_bl <- maled_bl %>%
    group_by(Country_ID) %>%
    mutate(wami_quintile = ntile(wamiimp, 5))
  
  # defile daily abx variables
  maled_full$who_abx <- ifelse(maled_full$safmacrolide == 1 | maled_full$saffluoro == 1, 1, 0)
  maled_full$maybe_eff_abx <- ifelse(maled_full$safcephalo == 1 | maled_full$safsulfon == 1 | maled_full$saftetracycl == 1 | 
                                       maled_full$safother == 1, 1, 0)
  
  maled_full$ineff_abx <- ifelse(maled_full$safpenicillin == 1 | maled_full$safmetron == 1 | maled_full$safunknown == 1, 1, 0)
  maled_full$no_abx <- ifelse(maled_full$who_abx == 0 & maled_full$maybe_eff_abx == 0 & maled_full$ineff_abx == 0, 1, 0)
  
  maled_full$ineff_abx <- ifelse(maled_full$ineff_abx == 1 & (maled_full$who_abx == 1 | maled_full$maybe_eff_abx == 1), 0, maled_full$ineff_abx)
  maled_full$maybe_eff_abx <- ifelse(maled_full$maybe_eff_abx == 1 & maled_full$who_abx == 1, 0, maled_full$maybe_eff_abx)
  
  maled_full$all_abx <- ifelse(maled_full$no_abx == 1 | maled_full$ineff_abx == 1, 0,
                               ifelse(maled_full$maybe_eff_abx == 1, 1, 2))
  
  maled_full$all_abx <- factor(maled_full$all_abx, levels = 0:2, labels = c("No or ineffective antibiotics",
                                                                            "Possibly effective antibiotics",
                                                                            "Guideline recommended antibiotics"))
  
  maled_full$any_abx <- ifelse(maled_full$safpenicillin == 1 |
                                 maled_full$safcephalo == 1 | 
                                 maled_full$safsulfon == 1 |
                                 maled_full$safmacrolide == 1 |
                                 maled_full$saftetracycl == 1 |
                                 maled_full$saffluoro == 1 |
                                 maled_full$safunknown == 1 | 
                                 maled_full$safmetron == 1 | 
                                 maled_full$safother == 1, 1, 0)
  
  #### Define Shigella cases (tac & culture, culture only)
  
  ## TAC attributable
  tac_data$tac_shigella_attributable <- ifelse(tac_data$shigella_eiec_afe > 0.5, 1, 0)
  
  ## Culture attributable
  
  # Get date of stool sample (use as episode date)
  micro_data$episode_date <- as.POSIXct(micro_data$srfdate, format = "%d%b%y:%H:%M:%S", tz = "UTC")
  
  micro_diar <- micro_data %>%
    select(pid,
           episode_date,
           shigella) %>%
    rename(date = 'episode_date',
           culture_shigella = 'shigella') %>%
    #deduplicate, if there is NA or 0/1 keep 0/1, if it is 0 and 1 keep 1
    mutate(
      culture_priority = case_when(
        is.na(culture_shigella) ~ -1,         # lowest priority
        culture_shigella == 0 ~ 1,            # medium
        culture_shigella == 1 ~ 2             # highest
      )
    ) %>%
    arrange(pid, date, desc(culture_priority)) %>%
    distinct(pid, date, .keep_all = TRUE) %>%
    select(-culture_priority)
  
  #join culture info 
  # Get date of sample collection (use as episode date)
  tac_data$date <- as.POSIXct(tac_data$date, format = "%m/%d/%Y", tz = "UTC")
  tac_data <- left_join(tac_data, micro_diar, by = c("pid", "date"))
  
  # Drop rows without TAC Shigella results 
  # NO LONGER DOING THIS FOR CO-ETIOLOGY META_ANALYSIS
  # tac_data <- tac_data[-which(is.na(tac_data$tac_shigella_attributable)),]
  
  tac_data$shigella_attributable <- ifelse(
    !is.na(tac_data$culture_shigella) & tac_data$culture_shigella == 1, 
    1, 
    tac_data$tac_shigella_attributable
  )
  
  # ----- Identify cases and controls -----
  
  if(case_def == "tac_or_culture_shig_diar"){
    
    # Case = diarrhea attributable to Shigella via tac or culture
    tac_data$case <- ifelse(tac_data$stooltype == "D1" & tac_data$shigella_attributable == 1, 1, 0)
    
    # Get rid of non-shigella diarrhea
    tac_data_shig_only <- tac_data[-which(tac_data$case == 0 & tac_data$stooltype == "D1"),]
    
    # Get cases
    case_data <- tac_data_shig_only[which(tac_data_shig_only$case == 1),]
    
    case_data$case_pid <- case_data$pid
    case_data$case_sid <- case_data$sid
    
    # Dataset of eligible controls
    all_controls <- tac_data_shig_only[tac_data_shig_only$stooltype == "M1",]
    
  } else if(case_def == "tac_shig_diar"){
    
    # Case = diarrhea attributable to Shigella via TAC
    tac_data$case <- ifelse(tac_data$stooltype == "D1" & tac_data$tac_shigella_attributable == 1, 1, 0)
    
    # Get rid of non-shigella diarrhea
    tac_data_shig_only <- tac_data[-which(tac_data$case == 0 & tac_data$stooltype == "D1"),]
    
    # Get cases
    case_data <- tac_data_shig_only[which(tac_data_shig_only$case == 1),]
    
    case_data$case_pid <- case_data$pid
    case_data$case_sid <- case_data$sid
    
    # Dataset of eligible controls
    all_controls <- tac_data_shig_only[tac_data_shig_only$stooltype == "M1",]
    
  }else if(case_def == "culture_shig_diar"){
    # Case = diarrhea attributable to Shigella via culture
    # AND non-NA tac results
    tac_data$case <- ifelse(tac_data$stooltype == "D1" & tac_data$culture_shigella == 1 & !is.na(tac_data$tac_shigella_attributable), 1, 0)
    tac_data_shig_only <- tac_data[-which(tac_data$case == 0 & tac_data$stooltype == "D1"),]
    
    # Get cases
    case_data <- tac_data_shig_only[which(tac_data_shig_only$case == 1),]
    
    case_data$case_pid <- case_data$pid
    case_data$case_sid <- case_data$sid
    
    # Dataset of eligible controls
    all_controls <- tac_data_shig_only[tac_data_shig_only$stooltype == "M1",]
    
  } else if(case_def == "all_diar"){
    # Case = diarrhea in general
    
    tac_data$case <- ifelse(tac_data$stooltype == "D1", 1, 0)
    
    # Get cases
    case_data <- tac_data[which(tac_data$case == 1),]
    
    case_data$case_pid <- case_data$pid
    case_data$case_sid <- case_data$sid
    
    # Dataset of eligible controls
    all_controls <- tac_data[tac_data$stooltype == "M1",]
  } else{
    stop("Invalid case definition")
  }
  
  # Make control data
  
  # person ID = pid
  # sample ID = sid
  
  matched_controls <- lapply(1:nrow(case_data), function(i, case_data, all_controls, tac_data, max_controls){
    row <- case_data[i,]
    
    # get age range to match depending on case age
    # age range = 0-11 months (1-364 days) --> +- 2mo (30.44*2 mo= 61 days)
    # age range = 12+ months (365 days +) --> +- 4mo (30.44*4 mo = 122 days)
    if(row$agedays < 365){
      min_age <- max(0, row$agedays - 61)
      max_age <- min(row$agedays + 61, 364)
    } else{
      min_age <- max(365, row$agedays - 122)
      max_age <- row$agedays + 122
    }
    
    # match sex, site, time, age, not their own control
    matching_controls <- all_controls %>%
      filter(cafsex == row$cafsex) %>%
      filter(country_id == row$country_id) %>%
      filter(date < (row$date + days(15)) & date > row$date - days(15)) %>%
      filter(agedays >= min_age & agedays <= max_age) %>%
      filter(pid != row$pid) %>%
      group_by(pid) %>%
      slice_max(order_by = date, n = 1) %>% # Keep the latest sample per individual
      ungroup()
    
    if(nrow(matching_controls) > 0){
      # for each matching control, make sure no diarrhea 7 days prior
      control_eligible <- rep(TRUE, nrow(matching_controls))
      for(j in 1:nrow(matching_controls)){
        control_row <- matching_controls[j,]
        tac_data_match <- tac_data %>%
          filter(pid == control_row$pid) %>%
          filter(date <= control_row$date & date > control_row$date - days(7))
        
        if(any(tac_data_match$stooltype == "D1")){
          control_eligible[j] <- FALSE
        } 
      }
      
      # eliminate ineligible controls
      matching_controls <- matching_controls[control_eligible,]
      
      # if nrow(matching_controls > max_controls), take max_controls num of controls
      if(!is.na(max_controls) & nrow(matching_controls) > max_controls){
        ctrl_samp <- sample(1:nrow(matching_controls), max_controls)
        matching_controls <- matching_controls[ctrl_samp,]
      }
      
      matching_controls$case_pid <- row$pid
      matching_controls$case_sid <- row$sid
      matching_controls$no_match <- FALSE
      
    } else{
      # No matching controls
      matching_controls[1,] <- NA
      matching_controls$case_pid <- row$pid
      matching_controls$case_sid <- row$sid
      matching_controls$no_match <- TRUE
    }
    
    return(matching_controls)
    
  }, case_data = case_data, all_controls = all_controls, tac_data = tac_data, max_controls = max_controls)
  
  matched_controls <- do.call(rbind, matched_controls)
  
  # identify cases with no matches and remove from data (rare, only occurring in all case data (not shig case))
  bad_case_sids <- matched_controls$case_sid[matched_controls$no_match == TRUE]
  all_tac_cc <- rbind(case_data, matched_controls[,colnames(matched_controls) != "no_match"])
  if(length(bad_case_sids) > 0){
    all_tac_cc <- all_tac_cc[-which(all_tac_cc$case_sid %in% bad_case_sids), ]
  }
  
  # Get attribution by AFE > 0.5
  all_tac_cc$adenovirus_attributable <- ifelse(all_tac_cc$adenovirus_40_41_afe > 0.5, 1, 0)
  all_tac_cc$aeromonas_attributable <- ifelse(all_tac_cc$aeromonas_afe > 0.5, 1, 0)
  all_tac_cc$astro_attributable <- ifelse(all_tac_cc$astrovirus_afe > 0.5, 1, 0)
  all_tac_cc$campylobacter_jejuni_coli_attributable <- ifelse(all_tac_cc$campylobacter_jejuni_coli_afe > 0.5, 1, 0)
  all_tac_cc$crypto_attributable <- ifelse(all_tac_cc$cryptosporidium_afe > 0.5, 1, 0)
  all_tac_cc$cyclospora_attributable <- ifelse(all_tac_cc$cyclospora_afe > 0.5, 1, 0)
  all_tac_cc$e_histolytica_attributable <- ifelse(all_tac_cc$e_histolytica_afe > 0.5, 1, 0)
  all_tac_cc$isospora_attributable <- ifelse(all_tac_cc$isospora_afe > 0.5, 1, 0)
  all_tac_cc$noro_attributable <- ifelse(all_tac_cc$norovirus_afe > 0.5, 1, 0)
  all_tac_cc$rotavirus_attributable <- ifelse(all_tac_cc$rotavirus_afe > 0.5, 1, 0)
  all_tac_cc$salmonella_attributable <- ifelse(all_tac_cc$salmonella_afe > 0.5, 1, 0)
  all_tac_cc$sapo_attributable <- ifelse(all_tac_cc$sapovirus_afe > 0.5, 1, 0)
  all_tac_cc$st_etec_attributable <- ifelse(all_tac_cc$ST_ETEC_afe > 0.5, 1, 0)
  all_tac_cc$tepec_attributable <- ifelse(all_tac_cc$tEPEC_afe > 0.5, 1, 0)
  all_tac_cc$v_cholerae_attributable <- ifelse(all_tac_cc$v_cholerae_afe > 0.5, 1, 0)
  
  all_tac_cc$ETEC_attributable <- ifelse(all_tac_cc$ETEC_afe > 0.5, 1, 0)
  all_tac_cc$e_bieneusi_attributable <- ifelse(all_tac_cc$e_bieneusi_afe > 0.5, 1, 0)
  all_tac_cc$giardia_attributable <- ifelse(all_tac_cc$giardia_afe > 0.5, 1, 0)
  all_tac_cc$eaec_attributable <- ifelse(all_tac_cc$EAEC_afe > 0.5, 1, 0)
  
  all_tac_cc$no_etiology <- ifelse(rowSums(all_tac_cc[,c("shigella_attributable",
                                                         "adenovirus_attributable",
                                                         "aeromonas_attributable",
                                                         "astro_attributable",
                                                         "campylobacter_jejuni_coli_attributable",
                                                         "crypto_attributable",
                                                         "cyclospora_attributable",
                                                         "e_histolytica_attributable",
                                                         "isospora_attributable",
                                                         "noro_attributable",
                                                         "rotavirus_attributable",
                                                         "salmonella_attributable",
                                                         "sapo_attributable",
                                                         "st_etec_attributable",
                                                         "tepec_attributable",
                                                         "v_cholerae_attributable")], na.rm = TRUE) == 0, 1, 0)
  
  # for other pathogens, use any tac < 35
  all_tac_cc$rota_detect <- ifelse(all_tac_cc$rotavirus < 35, 1, 0)
  all_tac_cc$adeno_detect <- ifelse(all_tac_cc$adenovirus_40_41 < 35, 1, 0)
  all_tac_cc$ETEC_detect <- ifelse(all_tac_cc$ETEC < 35, 1, 0)
  all_tac_cc$crypto_detect <- ifelse(all_tac_cc$cryptosporidium < 35,1,0)
  all_tac_cc$astro_detect <- ifelse(all_tac_cc$astrovirus < 35, 1, 0)
  all_tac_cc$noro_detect <- ifelse(all_tac_cc$norovirus < 35, 1, 0)
  all_tac_cc$tepec_detect <- ifelse(all_tac_cc$tEPEC < 35, 1, 0)
  all_tac_cc$campy_detect <- ifelse(all_tac_cc$campylobacter_pan < 35, 1, 0)
  all_tac_cc$sapo_detect <- ifelse(all_tac_cc$sapovirus < 35, 1, 0)
  all_tac_cc$e_bieneusi_detect <- ifelse(all_tac_cc$e_bieneusi < 35, 1, 0)
  all_tac_cc$giardia_detect <- ifelse(all_tac_cc$giardia < 35, 1, 0)
  all_tac_cc$EAEC_detect <- ifelse(all_tac_cc$EAEC < 35, 1, 0)
  
  
  # Get re-scaled pathogen quantities
  all_tac_cc$shigella_new <- (35 - all_tac_cc$shigella_eiec) / 3.322
  all_tac_cc$adenovirus_40_41_new <- (35 - all_tac_cc$adenovirus_40_41) / 3.322
  all_tac_cc$aeromonas_new <- (35 - all_tac_cc$aeromonas) / 3.322
  all_tac_cc$astrovirus_new <- (35 - all_tac_cc$astrovirus) / 3.322
  all_tac_cc$campylobacter_pan_new <- (35 - all_tac_cc$campylobacter_pan) / 3.322
  all_tac_cc$cryptosporidium_new <- (35 - all_tac_cc$cryptosporidium) / 3.322
  all_tac_cc$cyclospora_new <- (35 - all_tac_cc$cyclospora) / 3.322
  all_tac_cc$e_histolytica_new <- (35 - all_tac_cc$e_histolytica) / 3.322
  all_tac_cc$isospora_new <- (35 - all_tac_cc$isospora) / 3.322
  all_tac_cc$norovirus_new <- (35 - all_tac_cc$norovirus) / 3.322
  all_tac_cc$rotavirus_new <- (35 - all_tac_cc$rotavirus) / 3.322
  all_tac_cc$salmonella_new <- (35 - all_tac_cc$salmonella) / 3.322
  all_tac_cc$sapovirus_new <- (35 - all_tac_cc$sapovirus) / 3.322
  all_tac_cc$st_etec_new <- (35 - all_tac_cc$ST_ETEC) / 3.322
  all_tac_cc$tEPEC_new <- (35 - all_tac_cc$tEPEC) / 3.322
  all_tac_cc$v_cholerae_new <- (35 - all_tac_cc$v_cholerae) / 3.322
  all_tac_cc$ETEC_new <- (35 - all_tac_cc$ETEC) / 3.322
  all_tac_cc$e_bieneusi_new <- (35 - all_tac_cc$e_bieneusi) / 3.322
  all_tac_cc$eaec_new <- (35 - all_tac_cc$EAEC) / 3.322
  all_tac_cc$giardia_new <- (35 - all_tac_cc$giardia) / 3.322
  
  
  # Get initial abx treatment variables
  all_tac_cc$any_abx <- all_tac_cc$abxtrt 
  all_tac_cc$who_abx <- ifelse(all_tac_cc$macrotrt == 1 | all_tac_cc$fluorotrt == 1, 1, 0)
  
  all_tac_cc$maybe_eff_abx <- ifelse(all_tac_cc$cephalotrt == 1 | all_tac_cc$sulfontrt == 1 | all_tac_cc$tetratrt == 1 | 
                                       all_tac_cc$othertrt == 1, 1, 0)
  
  all_tac_cc$ineff_abx <- ifelse(all_tac_cc$peniciltrt == 1 |
                                   all_tac_cc$metrontrt == 1 |
                                   all_tac_cc$unknowtrt == 1, 1, 0)
  
  all_tac_cc$no_abx <- ifelse(all_tac_cc$who_abx == 0 & all_tac_cc$maybe_eff_abx == 0 & all_tac_cc$ineff_abx == 0, 1, 0)
  
  all_tac_cc$ineff_abx <- ifelse(all_tac_cc$ineff_abx == 1 & (all_tac_cc$who_abx == 1 | all_tac_cc$maybe_eff_abx == 1), 0, all_tac_cc$ineff_abx)
  all_tac_cc$maybe_eff_abx <- ifelse(all_tac_cc$maybe_eff_abx == 1 & all_tac_cc$who_abx == 1, 0, all_tac_cc$maybe_eff_abx)
  
  all_tac_cc$all_abx <- ifelse(all_tac_cc$no_abx == 1 | all_tac_cc$ineff_abx == 1, 0,
                               ifelse(all_tac_cc$maybe_eff_abx == 1, 1, 2))
  
  all_tac_cc$all_abx <- factor(all_tac_cc$all_abx, levels = 0:2, labels = c("No or ineffective antibiotics",
                                                                            "Possibly effective antibiotics",
                                                                            "Guideline recommended antibiotics"))
  
  # Drop controls with abx 0-15 days before sample (healthy controls only)
  # all_tac_cc <- all_tac_cc[-which(all_tac_cc$case == 0 & all_tac_cc$abx15 == 1),]
  
  # Select relevant variables from TAC dataset
  all_tac_cc <- all_tac_cc %>%
    select(pid,
           sid,
           case_pid,
           case_sid,
           case,
           country_id,
           date,
           agedays,
           shigella_attributable,
           tac_shigella_attributable,
           culture_shigella,
           adenovirus_attributable,
           aeromonas_attributable,
           astro_attributable,
           campylobacter_jejuni_coli_attributable,
           crypto_attributable,
           cyclospora_attributable,
           e_histolytica_attributable,
           isospora_attributable,
           noro_attributable,
           rotavirus_attributable,
           salmonella_attributable,
           sapo_attributable,
           st_etec_attributable,
           tepec_attributable,
           v_cholerae_attributable,
           ETEC_attributable,
           e_bieneusi_attributable,
           giardia_attributable,
           eaec_attributable,
           no_etiology,
           rota_detect,
           adeno_detect,
           ETEC_detect,
           crypto_detect,
           astro_detect,
           noro_detect,
           tepec_detect,
           campy_detect,
           sapo_detect,
           e_bieneusi_detect,
           giardia_detect,
           EAEC_detect,
           shigella_new,
           adenovirus_40_41_new,
           aeromonas_new,
           astrovirus_new,
           campylobacter_pan_new,
           cryptosporidium_new,
           cyclospora_new,
           e_histolytica_new,
           isospora_new,
           norovirus_new,
           rotavirus_new,
           salmonella_new,
           sapovirus_new,
           st_etec_new,
           tEPEC_new,
           v_cholerae_new,
           ETEC_new,
           e_bieneusi_new,
           eaec_new,
           giardia_new,
           # now also include quantities
           rotavirus,
           adenovirus_40_41,
           ETEC,
           cryptosporidium,
           astrovirus,
           norovirus,
           tEPEC,
           campylobacter_pan,
           sapovirus,
           e_bieneusi,
           giardia,
           EAEC,
           any_abx,
           who_abx,
           maybe_eff_abx, 
           ineff_abx,
           no_abx,
           all_abx,
           prop_ebf30) %>%
    rename("tac_rotavirus" = rotavirus,
           "tac_adenovirus" = adenovirus_40_41,
           "tac_etec" = ETEC,
           "tac_crypto" = cryptosporidium,
           "tac_astrovirus" = astrovirus,
           "tac_norovirus" = norovirus,
           "tac_tEPEC" = tEPEC,
           "tac_campylobacter_pan" = campylobacter_pan,
           "tac_sapovirus" = sapovirus,
           "tac_e_bieneusi" = e_bieneusi,
           "tac_giardia" = giardia,
           "tac_EAEC" = EAEC)
  
  maled_full$date <- as.POSIXct(maled_full$date, format = "%d%b%Y", tz = "UTC")
  
  severity_pre_abx <- function(i){
    row <- all_tac_cc[i, ]
    
    if(row$case == 1){
      # CASE
      episode_info <- maled_full[which(maled_full$stooldiaage == row$agedays &
                                         maled_full$Pid == row$pid),]
      
      if(nrow(episode_info) == 0){
        return(data.frame(maxb = NA,
                          fever = NA,
                          fever_days = NA,
                          maxls = NA,
                          sumvom = NA,
                          maxdehyd = NA,
                          alri = NA,
                          safcough = NA,
                          safshb = NA,
                          fstab = NA,
                          duration_pre_abx = 999)) # flag to remove row
      }
      
      # If no antibiotics, return overall info
      if(sum(episode_info$who_abx) == 0 & sum(episode_info$any_abx) == 0 & sum(episode_info$maybe_eff_abx) == 0){
        return(data.frame(maxb = episode_info$maxb[1],
                          fever = episode_info$fever[1],
                          fever_days = sum(episode_info$saffev, na.rm = TRUE),
                          maxls = episode_info$maxls[1],
                          sumvom = episode_info$sumvom[1],
                          maxdehyd = episode_info$maxdehyd[1],
                          alri = max(episode_info$alri),
                          safcough = max(episode_info$safcough),
                          safshb = max(episode_info$safshb),
                          fstab = max(episode_info$fstab),
                          duration_pre_abx = nrow(episode_info))) # returning length of episode 
      } else{
        
        # would like to specify for type coded in tac_data, for example if someone
        # gets ineffective abx then WHO abx later in episode, count prior to WHO abx
        # (added & all_abx == row$all_abx)
        # but there are some instances where safmacrolide (daily macrolide var) == 0 for whole episode 
        # even though macrotrt (episode treated with macrolides) == 1
        # so leaving as min(age[any_abx == 1]) for now
        
        pre_abx <- episode_info %>%
          mutate(first_abx = min(age[any_abx == 1])) %>% # & all_abx == row$all_abx])) %>%
          filter(age <= first_abx)
        
        return(data.frame(maxb = max(pre_abx$safblood, na.rm = TRUE),
                          fever = max(pre_abx$saffev, na.rm = TRUE), 
                          fever_days = sum(pre_abx$saffev, na.rm = TRUE),
                          maxls = max(pre_abx$safnumls, na.rm = TRUE),
                          sumvom = sum(pre_abx$safvom, na.rm = TRUE),
                          maxdehyd = max(pre_abx$safdehyd, na.rm = TRUE),
                          alri = max(pre_abx$alri, na.rm = TRUE),
                          safcough = max(pre_abx$safcough, na.rm = TRUE),
                          safshb = max(pre_abx$safshb, na.rm = TRUE),
                          fstab = max(episode_info$fstab, na.rm = TRUE),
                          duration_pre_abx = nrow(pre_abx)))
      }
    } else{
      # CONTROL
      
      # check to make sure not taking abx on sample date (could be on abx for something else)
      full_row <- maled_full[which(maled_full$Pid == row$pid & maled_full$date == row$date),]
      
      # full_row$any_abx == 1 if any of the saf daily abx vars == 1
      # nrow(full_row) == 0 if no entry in longitudinal data for that day
      if((nrow(full_row) == 0) || full_row$any_abx == 1) {
        # 999 to indicate drop row
        return(data.frame(maxb = 999,
                          fever = 999, 
                          fever_days = 999,
                          maxls = 999,
                          sumvom = 999,
                          maxdehyd = 999,
                          alri = 999,
                          safcough =999,
                          safshb = 999,
                          fstab = 999,
                          duration_pre_abx = 999))
      } else {
        # NA because not adjusting for severity in controls
        return(data.frame(maxb = NA,
                          fever = NA, 
                          fever_days = NA,
                          maxls = NA,
                          sumvom = NA,
                          maxdehyd = NA,
                          alri = NA,
                          safcough =NA,
                          safshb = NA,
                          fstab = NA,
                          duration_pre_abx = NA))
      }
      
      
    }
  }
  
  severity_df <- lapply(1:nrow(all_tac_cc), severity_pre_abx)
  severity_df <- do.call(rbind, severity_df)
  
  all_tac_cc <- cbind(all_tac_cc, severity_df)
  
  # Drop any controls taking abx on day of sample (indicated by 999 in severity cols)
  all_tac_cc <- all_tac_cc[-which(all_tac_cc$case == 0 & all_tac_cc$maxb == 999),]
  
  # drop diarrhea episodes coded with 999 in duration (not applicable now that switched severity merge)
  # all_tac_cc <- all_tac_cc[-which(all_tac_cc$case == 1 & all_tac_cc$duration_pre_abx == 999),]
  
  # If any_abx = 1 and fstab = 0, received abx before episode began
  # Mark duration_pre_abx = 0
  all_tac_cc$duration_pre_abx <- ifelse(all_tac_cc$any_abx == 1 & all_tac_cc$fstab == 0, 0, all_tac_cc$duration_pre_abx)
  
  # Get dates of z-score measurements
  zscore_data$date <- strptime(zscore_data$date, format = "%d%b%Y", tz = "UTC")
  
  # Get baseline growth (HAZ at or within one month before episode date) & month 3 growth (HAZ closest to 3 months after episode)
  baseline_and_monthx_df <- lapply(1:nrow(all_tac_cc), function(i, zscore_data, tac_data){
    
    x <- tac_data[i,]
    
    # Subset zscore_data for the same participant
    sub_zscore <- zscore_data[zscore_data$pid == x$pid, ]
    
    # get all dates prior to episode date AND within 75 days
    baseline_candidates <- sub_zscore[
      sub_zscore$date <= x$date &
        abs(sub_zscore$date - x$date) <= 75,
    ]
    
    if (nrow(baseline_candidates) > 0) {
      
      # order by closeness to episode date
      baseline_candidates <- baseline_candidates[
        order(abs(baseline_candidates$date - x$date)),
      ]
      
      # initialize as NA
      baseline_date   <- NA
      baseline_haz    <- NA
      baseline_waz    <- NA
      baseline_whz    <- NA
      baseline_weight <- NA
      baseline_length <- NA
      
      # loop through candidates and take first non-missing HAZ
      for (j in 1:nrow(baseline_candidates)) {
        
        row <- baseline_candidates[j, ]
        
        haz <- row$haz
        if (is.na(haz)) haz <- row$zhei
        if (is.na(haz)) haz <- row$zheiorig
        
        if (!is.na(haz)) {
          
          baseline_date   <- row$date
          baseline_haz    <- haz
          baseline_waz    <- row$zwei
          baseline_whz    <- row$whz
          baseline_weight <- row$weight
          baseline_length <- row$length
          
          # WHZ fallback
          if (is.na(baseline_whz)) baseline_whz <- row$zwfl
          if (is.na(baseline_whz)) baseline_whz <- row$zwflorig
          
          break
        }
      }
      
    } else {
      baseline_date <- NA
      baseline_haz <- NA
      baseline_waz <- NA
      baseline_whz <- NA
      baseline_weight <- NA
      baseline_length <- NA
    }
    
    # get monthx date closest to 90 days post episode and corresponding HAZ
    # month 3 date = 90 days by default, dynamic arg for longitudinal analysis
    target_monthx_date <- x$date + days(month_x_days)
    sub_zscore <- sub_zscore[sub_zscore$date > baseline_date,]
    
    if(nrow(sub_zscore) == 0){
      monthx_date <- NA
      monthx_haz <- NA
      
      monthx_waz <- NA
      monthx_whz <- NA
      monthx_weight <- NA
      monthx_length <- NA
    } else{
      monthx_date <- sub_zscore$date[which.min(abs(sub_zscore$date - target_monthx_date))]
      monthx_haz <- sub_zscore$haz[sub_zscore$date == monthx_date]
      
      # NEW for describing tanzania
      monthx_waz <- sub_zscore$zwei[sub_zscore$date == monthx_date]
      monthx_whz <- sub_zscore$whz[sub_zscore$date == monthx_date]
      monthx_weight <- sub_zscore$weight[sub_zscore$date == monthx_date]
      monthx_length <- sub_zscore$length[sub_zscore$date == monthx_date]
      
      # Check to make sure date is within 45 days of 3mo followup (updated from 75 days)
      if(length(monthx_date) == 0 || abs(monthx_date - target_monthx_date) > 45){
        monthx_date <- NA
        monthx_haz <- NA
        
        monthx_waz <- NA
        monthx_whz <- NA
        monthx_weight <- NA
        monthx_length <- NA
      }
      
    }
    
    return(data.frame(baseline_haz = baseline_haz,
                      baseline_date = baseline_date,
                      baseline_waz = baseline_waz,
                      baseline_whz = baseline_whz,
                      baseline_weight = baseline_weight,
                      baseline_length = baseline_length,
                      monthx_haz = monthx_haz,
                      monthx_date = monthx_date,
                      monthx_waz = monthx_waz, 
                      monthx_whz = monthx_whz,
                      monthx_weight = monthx_weight, 
                      monthx_length = monthx_length ))
    
  }, zscore_data = zscore_data, tac_data = all_tac_cc)
  
  
  baseline_and_monthx_df <- do.call(rbind, baseline_and_monthx_df) 
  final_df <- cbind(all_tac_cc, baseline_and_monthx_df)
  
  # rename covariates
  final_df <- final_df %>%
    rename(
      "episode_date" = date,
      "dysentery" = maxb,
      "lsstools" = maxls,
      "dehyd" = maxdehyd,
      "daysvomit" = sumvom,
      "cough" = safcough,
      "shortbreath" = safshb)
  
  # select covariates from bl data
  maled_bl <- maled_bl %>%
    select(Pid,
           CAFSEX,
           #mated,
           ageexbfimp, 
           Country_ID,
           incomeabovemed, #note not seeing this in the dictionary
           edimp,          #continuous maternal education
           incomemean,     #income? not in dictionary but liz said to use
           wamiimp,        #continuous version of SES score
           wami_quintile,
           drinkimp,
           sanitimp) %>%    # WAMI quintile by site
    rename("pid" = Pid,
           "sex" = CAFSEX,
           #"mated_bin" = mated,
           "site" = Country_ID,
           "maxagebf" = ageexbfimp,
           "mated_cont" = edimp,
           "income" = incomemean,
           "ses_wami" = wamiimp)
  
  maled_bl$mated_bin <- ifelse(maled_bl$mated_cont >= 6, 1, 0)
  
  maled_bl$wami_quintile <- factor(maled_bl$wami_quintile, levels = 1:5, labels = c("1st quintile of SES",
                                                                                    "2nd quintile of SES",
                                                                                    "3rd quintile of SES",
                                                                                    "4th quintile of SES",
                                                                                    "5th quintile of SES"))
  
  maled_bl$site <- factor(maled_bl$site, 
                          levels = c("BGD",
                                     "BRF",
                                     "INV",
                                     "NEB",
                                     "PEL",
                                     "PKN",
                                     "SAV",
                                     "TZH"),
                          labels = c("Bangladesh",
                                     "Brazil",
                                     "India",
                                     "Nepal",
                                     "Peru",
                                     "Pakistan",
                                     "South Africa",
                                     "Tanzania"))
  maled_bl$sex <- factor(maled_bl$sex, levels = c(1,2), labels = c("male", "female"))
  maled_bl$incomeabovemed <- factor(maled_bl$incomeabovemed, levels = c(0,1), labels = c("Income below country median",
                                                                                         "Income above country median"))
  
  # join covariates into final_df
  final_df <- left_join(final_df, maled_bl, by = "pid")
  
  # get rid of pakistan
  # Get rid of Pakistan and drop unused factor levels
  final_df <- final_df[which(final_df$site != "Pakistan"),]
  final_df$site <- droplevels(final_df$site)
  
  final_df$followup_days <- as.numeric(difftime(final_df$monthx_date,final_df$baseline_date , units = "days"))
  final_df$I_followup_days <- ifelse(is.na(final_df$followup_days), 0, 1)
  final_df$I_followup_days_x_followup_days <- ifelse(is.na(final_df$followup_days), 0, final_df$followup_days)
  
  final_df$agemonths <- round(final_df$agedays / 30.44, 1)
  
  final_df$dehyd <- factor(final_df$dehyd, levels = c(0,1,2), labels = c("None", "Some dehydration", "Severe dehydration"))
  
  # Get rid of extreme HAZ observations
  final_df$monthx_haz <- ifelse(final_df$monthx_haz < -6 | final_df$monthx_haz > 6, NA, final_df$monthx_haz)
  final_df$baseline_haz <- ifelse(final_df$baseline_haz < -6 | final_df$baseline_haz > 6, NA, final_df$baseline_haz)
  
  final_df$hazdiff <- final_df$monthx_haz - final_df$baseline_haz
  final_df$wazdiff <- final_df$monthx_waz - final_df$baseline_waz
  final_df$wlzdiff <- final_df$monthx_whz - final_df$baseline_whz
  final_df$lendiff <- final_df$monthx_length - final_df$baseline_length
  final_df$wtdiff <- final_df$monthx_weight - final_df$baseline_weight
  
  # Add in GEMS definition of MSD from diarrhea data
  
  sub_diarrhea_data <- diarrhea_data[diarrhea_data$Pid %in% final_df$pid,]
  names(sub_diarrhea_data)[names(sub_diarrhea_data) == 'Pid'] <- "pid"
  names(sub_diarrhea_data)[names(sub_diarrhea_data) == 'age'] <- "agedays"
  names(sub_diarrhea_data)[names(sub_diarrhea_data) == 'gemsdef'] <- "MSD"
  names(sub_diarrhea_data)[names(sub_diarrhea_data) == 'incidabxtrt'] <- "initiated_abx_during_episode"
  
  final_df <- left_join(final_df, 
                        sub_diarrhea_data[,c("pid", "agedays", "MSD", "initiated_abx_during_episode")], 
                        by = c("pid" = "pid", 
                               "agedays" = "agedays"))
  
  final_df <- final_df %>%
    select(pid,
           sid,
           case_pid,
           case_sid,
           case,
           episode_date,
           agedays,
           agemonths,
           shigella_attributable,
           tac_shigella_attributable,
           culture_shigella,
           adenovirus_attributable,
           aeromonas_attributable,
           astro_attributable,
           campylobacter_jejuni_coli_attributable,
           crypto_attributable,
           cyclospora_attributable,
           e_histolytica_attributable,
           isospora_attributable,
           noro_attributable,
           rotavirus_attributable,
           salmonella_attributable,
           sapo_attributable,
           st_etec_attributable,
           tepec_attributable,
           v_cholerae_attributable,
           ETEC_attributable,
           e_bieneusi_attributable,
           giardia_attributable,
           eaec_attributable,
           no_etiology,
           rota_detect,
           adeno_detect,
           ETEC_detect,
           crypto_detect,
           astro_detect,
           noro_detect,
           tepec_detect,
           campy_detect,
           sapo_detect,
           e_bieneusi_detect,
           giardia_detect,
           EAEC_detect,
           tac_rotavirus,
           tac_adenovirus,
           tac_etec,
           tac_crypto,
           tac_astrovirus,
           tac_norovirus,
           tac_tEPEC,
           tac_campylobacter_pan,
           tac_sapovirus,
           tac_e_bieneusi,
           tac_giardia,
           tac_EAEC,
           shigella_new,
           adenovirus_40_41_new,
           aeromonas_new,
           astrovirus_new,
           campylobacter_pan_new,
           cryptosporidium_new,
           cyclospora_new,
           e_histolytica_new,
           isospora_new,
           norovirus_new,
           rotavirus_new,
           salmonella_new,
           sapovirus_new,
           st_etec_new,
           giardia_new,
           tEPEC_new,
           v_cholerae_new,
           ETEC_new,
           e_bieneusi_new,
           eaec_new,
           no_etiology,
           any_abx,
           who_abx,
           maybe_eff_abx,
           ineff_abx,
           no_abx,
           all_abx,
           duration_pre_abx,
           dysentery,
           fever,
           fever_days,
           dehyd,
           lsstools,
           daysvomit,
           cough,
           shortbreath,
           alri,
           income,
           incomeabovemed,
           mated_cont,
           mated_bin,
           ses_wami,
           wami_quintile,
           drinkimp,
           sanitimp,
           baseline_haz,
           baseline_date,
           monthx_haz,
           monthx_date,
           baseline_waz ,
           baseline_whz ,
           baseline_weight, 
           baseline_length ,
           monthx_waz ,
           monthx_whz,
           monthx_weight,
           monthx_length,
           sex, 
           maxagebf,
           prop_ebf30,
           site,
           followup_days,
           I_followup_days,
           I_followup_days_x_followup_days,
           hazdiff, 
           wazdiff,
           wlzdiff,
           lendiff,
           wtdiff,
           MSD) %>%
    set_variable_labels(pid = "Participant ID",
                        sid = "Sample ID",
                        case_pid = "Participant ID of case",
                        case_sid = "Sample ID of case",
                        case = "Case",
                        episode_date = "Date of diarrhea episode",
                        agedays = "Age at sample collection (days)",
                        agemonths = "Age at sample collection (months)",
                        shigella_attributable = "Shigella attributable (AFE > 0.5)",
                        any_abx = "Received any antibiotics",
                        who_abx = "Received WHO approved antibiotics",
                        maybe_eff_abx = "Recieved maybe effective antibiotics",
                        ineff_abx = "Recieved ineffective or no antibiotics",
                        baseline_haz = "HAZ at baseline (before & closest to episode)",
                        baseline_date = "Date of baseline HAZ measurement",
                        monthx_haz = "HAZ at three months (after & closest to 90 days post-episode)",
                        monthx_date = "Date of month three HAZ measurement",
                        dysentery = "Dysentery",
                        lsstools = "Max number of loose stools during episode",
                        dehyd = "Maximum severity of dehydration during diarrhea episode",
                        fever = "Reported fever during episode",
                        fever_days = "Days reported fever during episode",
                        daysvomit = "Days vommitted during episode",
                        cough = "Maternal report of cough",
                        shortbreath = "Maternal report of shortness of breath",
                        alri = "ALRI definition met",
                        rotavirus_attributable = "Rotavirus attributable (AFE > 0.5)",
                        crypto_attributable = "Cryptosporidium attributable (AFE > 0.5)",
                        adenovirus_attributable = "Adenovirus attributable (AFE > 0.5)",
                        ETEC_attributable = "ETEC attributable (AFE >0.5)",
                        astro_attributable = "Astrovirus attributable (AFE > 0.5)",
                        noro_attributable = "Norovirus attributable (AFE > 0.5)",
                        tepec_attributable = "tEPEC attributable (AFE > 0.5)",
                        sapo_attributable = "Sapovirus attributable (AFE > 0.5)",
                        e_bieneusi_attributable = "E bieneusi attributable (AFE > 0.5)",
                        giardia_attributable = "Giardia attributable (AFE > 0.5)",
                        eaec_attributable = "EAEC attributable (AFE > 0.5)",
                        no_etiology = "No other attributable etiology",
                        rota_detect = "Rotavirus detected",
                        adeno_detect = "Adenovirus detected",
                        ETEC_detect = "ETEC detected",
                        crypto_detect = "Cryptosporidium detected",
                        astro_detect = "Astrovirus detected",
                        noro_detect = "Norovirus detected",
                        tepec_detect = "tEPEC detected",
                        campy_detect = "Campylobacter detected",
                        sapo_detect = "Sapovirus detected",
                        e_bieneusi_detect = "E Bieneusi detected",
                        giardia_detect = "Giardia detected",
                        EAEC_detect = "EAEC detected",
                        sex = "Sex",
                        mated_cont = "Years of maternal education",
                        mated_bin = "Mother completed >=6 years of school",
                        ses_wami = "WAMI Socioeconomic Status Score",
                        wami_quintile = "WAMI quintile (by site)",
                        drinkimp = "Improved drinking water",
                        sanitimp = "Improved sanitation",
                        maxagebf = "Max age of breastfeeding (imputed mean for country if missing)", #note could only find imputed version, could remove imputed values if needed. also concerned this is > age at episode in many cases. prop var better
                        prop_ebf30 = "Proportion of days of exclusive breastfeeding of the 30 days prior to episode",
                        site = "Site",
                        incomeabovemed = "Income above country median",
                        income = "Mean income",
                        followup_days = "Days between baseline HAZ and month 3 HAZ measurement",
                        hazdiff = "Difference between month 3 and baseline HAZ",
                        MSD = "Moderate to severe diarrhea (by GEMS definition)")
  
  # Add same variables as VIDA/GEMS for bootstrap
  
  # first_id = associated with the child
  # case_id = associated with the case
  # child_id = associated with the episode 
  
  # when first sample is case then first = case = child
  final_df <- final_df %>%
    arrange(pid, agedays) %>%
    group_by(pid) %>%
    mutate(first_id = sid[1]) %>%
    mutate(case_id = case_sid,
           child_id = sid)
  
  return(final_df)
  
}


maled_1mo <-prep_maled_case_control_monthx(month_x_days = 30)
saveRDS(maled_1mo, here::here("data/maled_data/longitudinal/maled_1mo.Rds"))

maled_2mo <-prep_maled_case_control_monthx(month_x_days = 60)
saveRDS(maled_2mo, here::here("data/maled_data/longitudinal/maled_2mo.Rds"))

maled_3mo <-prep_maled_case_control_monthx(month_x_days = 90)
saveRDS(maled_3mo, here::here("data/maled_data/longitudinal/maled_3mo.Rds"))

maled_4mo <- prep_maled_case_control_monthx(month_x_days = 120)
saveRDS(maled_4mo, here::here("data/maled_data/longitudinal/maled_4mo.Rds"))

maled_5mo <- prep_maled_case_control_monthx(month_x_days = 150)
saveRDS(maled_5mo, here::here("data/maled_data/longitudinal/maled_5mo.Rds"))

maled_6mo <- prep_maled_case_control_monthx(month_x_days = 180)
saveRDS(maled_6mo, here::here("data/maled_data/longitudinal/maled_6mo.Rds"))

maled_7mo <- prep_maled_case_control_monthx(month_x_days = 210)
saveRDS(maled_7mo, here::here("data/maled_data/longitudinal/maled_7mo.Rds"))

maled_8mo <- prep_maled_case_control_monthx(month_x_days = 240)
saveRDS(maled_8mo, here::here("data/maled_data/longitudinal/maled_8mo.Rds"))

maled_9mo <- prep_maled_case_control_monthx(month_x_days = 270)
saveRDS(maled_9mo, here::here("data/maled_data/longitudinal/maled_9mo.Rds")) 

maled_10mo <- prep_maled_case_control_monthx(month_x_days = 300)
saveRDS(maled_10mo, here::here("data/maled_data/longitudinal/maled_10mo.Rds"))

maled_11mo <- prep_maled_case_control_monthx(month_x_days = 330)
saveRDS(maled_11mo, here::here("data/maled_data/longitudinal/maled_11mo.Rds"))

maled_12mo <- prep_maled_case_control_monthx(month_x_days = 365)
saveRDS(maled_12mo, here::here("data/maled_data/longitudinal/maled_12mo.Rds"))



