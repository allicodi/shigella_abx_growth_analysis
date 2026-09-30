# ---------------------------------------------------------------
# Script to create cleaned VIDA datasets
# ---------------------------------------------------------------

here::i_am("code/prep_data/prep_vida.R")

library(tidyverse)
library(labelled)
library(haven)
library(corrr)
library(factoextra)

#---------------------------------------------------------------
# CASE-ONLY ANALYSIS
# --------------------------------------------------------------

#' Function to create clean VIDA case dataset
#' 
#' By default, this will return all diarrhea episodes (modified for co-etiology 
#' meta-analysis). For the original meta-analysis, we define case as TAC 
#' or culture attributable Shigella so subset the dataset from this function to those who
#' have Shigella TAC results.
prep_vida <- function(){
  
  data <- read.csv(here::here("data/vida_data/raw_data/VIDA_Rogawski.csv"))
  data_w_dss <- readxl::read_xlsx(here::here("data/vida_data/raw_data/VIDA_Rogawski_withDSScount and alt.xlsx"))
  
  # Just merge in dss_alt because excel formatting/var types is weird and messing things up, don't want fix
  data_w_dss <- data_w_dss[,c("CHILDID", "DSS_alt")]
  data <- left_join(data, data_w_dss, by = "CHILDID")
  
  # make first ID based on DSS_alt
  # group by dss alt, arrange by enroll_date, the first CHILDID == first id
  data <- data %>%
    group_by(DSS_alt) %>%
    arrange(DSS_alt, ENROLL_DATE) %>%
    mutate(first_id = first(CHILDID)) %>%
    ungroup()
  
  # -------------------------------------------------
  # Make PCA SES score in cases and controls by site 
  # -------------------------------------------------
  data$fridge <- ifelse(!is.na(data$F4A_HOUSE_FRIDGE), 
                        data$F4A_HOUSE_FRIDGE,
                        data$F7_HOUSE_FRIDGE)
  data$tv <- ifelse(!is.na(data$F4A_HOUSE_TELE),
                    data$F4A_HOUSE_TELE,
                    data$F7_HOUSE_TELE)
  data$electricity <- ifelse(!is.na(data$F4A_HOUSE_ELEC),
                             data$F4A_HOUSE_ELEC,
                             data$F7_HOUSE_ELEC)
  data$motorcycle <- ifelse(!is.na(data$F4A_HOUSE_SCOOT),
                            data$F4A_HOUSE_SCOOT,
                            data$F7_HOUSE_SCOOT)
  data$radio <- ifelse(!is.na(data$F4A_HOUSE_RADIO),
                       data$F4A_HOUSE_RADIO,
                       data$F7_HOUSE_RADIO)
  data$bike <- ifelse(!is.na(data$F4A_HOUSE_BIKE),
                      data$F4A_HOUSE_BIKE,
                      data$F7_HOUSE_BIKE)
  data$car <- ifelse(!is.na(data$F4A_HOUSE_CAR),
                     data$F4A_HOUSE_CAR,
                     data$F7_HOUSE_CAR)
  data$boat <- ifelse(!is.na(data$F4A_HOUSE_BOAT),
                      data$F4A_HOUSE_BOAT,
                      data$F7_HOUSE_BOAT)
  data$phone <- ifelse(!is.na(data$F4A_HOUSE_PHONE),
                       data$F4A_HOUSE_PHONE,
                       data$F7_HOUSE_PHONE)
  data$cart <- ifelse(!is.na(data$F4A_HOUSE_CART),
                      data$F4A_HOUSE_CART,
                      data$F7_HOUSE_CART)
  data$agland <- ifelse(!is.na(data$F4A_HOUSE_AGLAND),
                        data$F4A_HOUSE_AGLAND,
                        data$F7_HOUSE_AGLAND)
  
  cols_for_pca <- c("fridge", "tv", "electricity", "motorcycle",
                    "radio", "bike", "car", "boat", "phone", "cart", "agland")
  
  ses_pca <- prcomp(data[complete.cases(data[,cols_for_pca]), cols_for_pca], 
                    center = TRUE, scale = TRUE)
  
  pca_score <- as.matrix(data[, cols_for_pca]) %*% ses_pca$rotation[, 1]
  
  data$ses_score <- pca_score[, 1]   
  
  data <- data %>%
    group_by(SITE) %>%
    mutate(ses_quintile = ntile(ses_score, 5)) 
  
  # ----------------------------------------------------------------------------
  
  # Subset to cases only
  case_data <- data[data$VIDA_CASE == 1 | data$VIDA_CASE == 2,]
  
  # Covariates:
  # sex, age, site, bl_haz, edu, bf week prior, days, additional pathogens, kids under 5 in hh
  # SES = improved drinking, improved sanitation, electricity, tv, fridge
  
  # Severity: 
  # dysentery, fever, dehydration, loose stools, vomit
  case_data <- case_data %>%
    rename('child_id' = CHILDID,
           'first_id' = first_id,
           'sex' = GENDER,
           'agemchild' = AGE,
           'agedchild' = BASE_AGE_DAYS,
           'site' = SITE,
           'enr_haz' = BASE_HAZ,
           'enr_waz' = BASE_WAZ,
           'hazd60' = END_HAZ,
           'hazdiff' = DELTA_HAZ,
           'education' = F4A_PRIM_SCHL,
           'followup_days' = DURDAYS,
           'num_hh_lt5' = F4A_YNG_CHILDREN,
           'bf_wk_before' = F4A_PRI_BMILK,
           'lsstools' = F4A_DAILY_MAX,
           'vom_days' = F4A_DAYS_VOMIT,
           'vom_freq' = F4A_FREQ_VOMIT,
           'fever' = F4A_DRH_FEVER,
           'dehydr' = WHO_DEHYD,
           'dysentery' = DYS_IND,
           'safe_water' = IMPROV_WATER,     
           'safe_sanit' = IMPROV_SANIT)
  
  # remove don't know --> NA
  case_data$fever <- ifelse(case_data$fever == 9, NA, case_data$fever)
  
  case_data$lsstools <- factor(case_data$lsstools, levels = 1:4, labels = c("3 stools",
                                                                            "4 to 5 stools",
                                                                            "6 to 10 stools",
                                                                            "More than 10 stools"))
  # Get rid of extreme HAZ observations
  case_data$enr_haz <- ifelse(case_data$enr_haz < -6 | case_data$enr_haz > 6, NA, case_data$enr_haz)
  case_data$enr_waz <- ifelse(case_data$enr_waz < -6 | case_data$enr_waz > 6, NA, case_data$enr_waz)
  case_data$hazd60 <- ifelse(case_data$hazd60 < -6 | case_data$hazd60 > 6, NA, case_data$hazd60)
  
  case_data$hazdiff <- ifelse(is.na(case_data$enr_haz) | is.na(case_data$hazd60), NA, case_data$hazdiff)
  
  case_data$vom_days <- ifelse(is.na(case_data$vom_days), 0, case_data$vom_days)
  case_data$vom_freq <- ifelse(is.na(case_data$vom_freq), 0, case_data$vom_freq)
  
  case_data$I_followup_days <- ifelse(is.na(case_data$followup_days), 0 , 1)
  case_data$I_followup_days_x_followup_days <- ifelse(is.na(case_data$followup_days), 0 , case_data$followup_days)
  
  case_data$duration_pre_enroll <- case_data$F4A_DRH_DAYS
  
  # Pathogens
  case_data$shigella_culture_positive <- ifelse(rowSums(case_data[, c("F16_SHIG_SPP", 
                                                                      "SHIG_BOYDII", 
                                                                      "SHIG_DYSENT", 
                                                                      "SHIG_FLEX", 
                                                                      "SHIG_SONNEI")] == 1, na.rm = TRUE) > 0, 1, 0)
  
  case_data$rotavirus <- case_data$F18_RES_ROTAVIRUS
  case_data$st_etec <- case_data$ETEC_ST
  case_data$crypto <- case_data$F18_RES_CRYPTOSPOR
  case_data$adeno <- case_data$ADENO_4041
  case_data$etec <- case_data$ETEC_ALL
  case_data$astro <- case_data$F19_ASTRO_VIRUS
  case_data$noro <- case_data$NORO_ANY
  case_data$noro_gii <- case_data$F19_NORO_GII
  case_data$tepec <- case_data$TEPEC
  case_data$campy <- case_data$CAMPY_ANY
  case_data$campy_j <- case_data$F16_CAMPY_JEJUNI
  case_data$sapo <- case_data$F19_SAPO_VIRUS
  case_data$giardia <- case_data$F18_RES_GIARDIA
  case_data$eaec <- case_data$EAEC
  
  # Create 'new' scaled variables
  case_data$shigella_new <- (35 - case_data$TAC_SHIGELLA_EIEC) / 3.322
  case_data$rotavirus_new <- (35 - case_data$TAC_ROTAVIRUS) / 3.322
  case_data$st_etec_new <- (35 - case_data$TAC_ST_ETEC) / 3.322
  case_data$lt_etec_new <- (35 - case_data$TAC_LT_ETEC) / 3.322
  # etec max of st and lt
  case_data$etec_new <- ifelse(case_data$st_etec_new > case_data$lt_etec_new, case_data$st_etec_new, case_data$lt_etec_new)
  case_data$crypto_new <- (35 - case_data$TAC_CRYPTO) / 3.322
  case_data$adeno_new <- (35 - case_data$TAC_ADENO4041) / 3.322
  # case_data$etec <- (35 - case_data$TAC_ST_ETEC) / 3.322
  case_data$astro_new <- (35 - case_data$TAC_ASTROVIRUS) / 3.322
  # case_data$noro_new <- (35 - case_data$TAC_NORO) / 3.322
  case_data$noro_gii_new <- (35 - case_data$TAC_NORO_GII) / 3.322
  case_data$tepec_new <- (35 - case_data$TAC_TEPEC) / 3.322
  case_data$campy_new <- (35 - case_data$TAC_CAMPY_ANY) / 3.322
  case_data$campy_j_new <- (35 - case_data$TAC_CAMPY) / 3.322
  case_data$sapo_new <- (35 - case_data$TAC_SAPOVIRUS) / 3.322
  case_data$giardia_new <- (35 - case_data$TAC_GIARDIA) / 3.322
  case_data$eaec_new <- (35 - case_data$TAC_EAEC) / 3.322
  case_data$e_bieneusi_new <- (35 - case_data$TAC_E_BIENEUSI) / 3.322
  case_data$v_cholerae_new <- (35 - case_data$TAC_V_CHOLERAE) / 3.322
  
  case_data$salmonella_new <- (35 - case_data$TAC_SALMONELLA) / 3.322
  
  # TAC using cutoffs from GEMS - https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(16)31529-X/abstract
  case_data$tac_shig <- ifelse(case_data$TAC_SHIGELLA_EIEC < 27.9, 1, 0)
  
  case_data$tac_adeno_4041 <- ifelse(case_data$TAC_ADENO4041 < 22.7, 1, 0)
  case_data$tac_astro <- ifelse(case_data$TAC_ASTROVIRUS < 22.2, 1, 0)
  case_data$tac_campyj <- ifelse(case_data$TAC_CAMPY < 15.4, 1, 0)
  case_data$tac_crypto <- ifelse(case_data$TAC_CRYPTO < 24.0, 1, 0)
  case_data$tac_cyclosporidium <- ifelse(case_data$TAC_CYCLO < 29.6, 1, 0)
  case_data$tac_ehist <- ifelse(case_data$TAC_E_HISTOLYTICA < 32.8, 1, 0)
  case_data$tac_norovirus_gii <- ifelse(case_data$TAC_NORO_GII < 23.4, 1, 0)
  case_data$tac_rota <- ifelse(case_data$TAC_ROTAVIRUS < 32.6, 1, 0)
  case_data$tac_salm <- ifelse(case_data$TAC_SALMONELLA < 30.7, 1, 0)
  case_data$tac_st_etec <- ifelse(case_data$TAC_ST_ETEC < 22.8, 1, 0)
  case_data$tac_tepec <- ifelse(case_data$TAC_TEPEC < 16.0, 1, 0)
  case_data$tac_vchol <- ifelse(case_data$TAC_V_CHOLERAE < 33.8, 1, 0)
  
  # get rid of people who don't have TAC results
  # stopped doing for co-etiology meta-analysi
  # case_data <- case_data[-which(is.na(case_data$tac_shig)),]
  
  case_data$shigella_tac_or_culture <- ifelse(
    rowSums(case_data[, c("shigella_culture_positive", "tac_shig")] == 1, na.rm = TRUE) > 0,
    1, 0
  )
  
  case_data$no_etiology <- ifelse(rowSums(case_data[,c("shigella_tac_or_culture",
                                                       "tac_adeno_4041",
                                                       "tac_astro",
                                                       "tac_campyj",
                                                       "tac_crypto",
                                                       "tac_cyclosporidium",
                                                       "tac_ehist",
                                                       "tac_norovirus_gii",
                                                       "tac_rota",
                                                       "tac_salm",
                                                       "tac_st_etec",
                                                       "tac_tepec",
                                                       "tac_vchol")], na.rm = TRUE) == 0, 1, 0)
  
  # detected
  case_data$rota_detected <- ifelse(case_data$TAC_ROTAVIRUS < 35, 1, 0)
  case_data$adeno_4041_detected <- ifelse(case_data$TAC_ADENO4041 < 35, 1, 0)
  case_data$etec_detected <- ifelse(case_data$TAC_ST_ETEC < 35 | case_data$TAC_LT_ETEC < 35, 1, 0)
  case_data$crypto_detected <- ifelse(case_data$TAC_CRYPTO < 35, 1, 0)
  case_data$astro_detected <- ifelse(case_data$TAC_ASTROVIRUS < 35, 1, 0)
  case_data$noro_detected <- ifelse(case_data$TAC_NORO_GI < 35 | case_data$TAC_NORO_GII < 35, 1, 0)
  case_data$tepec_detected <- ifelse(case_data$TAC_TEPEC < 35, 1, 0)
  case_data$campy_detected <- ifelse(case_data$TAC_CAMPY_ANY < 35, 1, 0)
  case_data$sapo_detected <- ifelse(case_data$TAC_SAPOVIRUS < 35, 1, 0)
  case_data$e_bieneusi_detected <- ifelse(case_data$TAC_E_BIENEUSI < 35, 1, 0)
  case_data$giardia_detected <- ifelse(case_data$TAC_GIARDIA < 35, 1, 0)
  case_data$eaec_detected <- ifelse(case_data$TAC_EAEC < 35, 1, 0)
  
  # Antibiotics
  # WHO approved = azithromycin, ciprofloxacin, ceftriaxone, pivmecillinam
  case_data$who_abx <- ifelse(
    rowSums(cbind(
      case_data$F4B_TRT_GIVE_AZI, case_data$F4B_TRT_PRES_AZI,
      case_data$F4B_TRT_GIVE_CPNR, case_data$F4B_TRT_PRES_CPNR,
      case_data$F4B_TRT_GIVE_CEF, case_data$F4B_TRT_PRES_CEF,
      case_data$F4B_TRT_GIVE_SLPY, case_data$F4B_TRT_PRES_SLPY
    ) == 1, na.rm = TRUE) > 0 
      # No longer including hometrt as of 6/11/26
      # |
      # case_data$F4A_HOMETRT_AB_SPEC %in% c(
      #   "CEFTRIAXONE", "CEFTRIAXONE 500mg", "ciprofloxacine", "CECTRIAXONE"
      # )
    , 1, 0
  )
  
  # New 9/4/25 - get abx given so can make ast_given_abx variable ------
  
  # try without azithro
  # case_data$azithro <- ifelse(
  #   rowSums(cbind(
  #     case_data$F4B_TRT_GIVE_AZI, case_data$F4B_TRT_PRES_AZI
  #   ) == 1, na.rm = TRUE) > 0, 1, 0
  # )
  
  case_data$cipro <- ifelse(
    rowSums(cbind(
      case_data$F4B_TRT_GIVE_CPNR, case_data$F4B_TRT_PRES_CPNR
    ) == 1, na.rm = TRUE) > 0 #|
      #case_data$F4A_HOMETRT_AB_SPEC %in% c(
      #  "ciprofloxacine"
      #)
  , 1, 0
  )
  
  case_data$ceft <- ifelse(
    rowSums(cbind(
      case_data$F4B_TRT_GIVE_CEF, case_data$F4B_TRT_PRES_CEF
    ) == 1, na.rm = TRUE) > 0, 
      #|
      #case_data$F4A_HOMETRT_AB_SPEC %in% c(
      #  "CEFTRIAXONE", "CEFTRIAXONE 500mg", "CECTRIAXONE"
      #),
    1, 0
  )
  
  # --------------------------------------------------------------------
  
  case_data$maybe_eff_abx <- ifelse(
    rowSums(cbind(
      case_data$F4B_TRT_GIVE_AMOX, case_data$F4B_TRT_PRES_AMOX,
      case_data$F4B_TRT_GIVE_AMPI, case_data$F4B_TRT_PRES_AMPI,
      case_data$F4B_TRT_GIVE_CXL, case_data$F4B_TRT_PRES_CXL,
      case_data$F4B_TRT_GIVE_GENT, case_data$F4B_TRT_PRES_GENT,
      case_data$F4B_TRT_GIVE_CHLOR, case_data$F4B_TRT_PRES_CHLOR,
      case_data$F4B_TRT_GIVE_ERY, case_data$F4B_TRT_PRES_ERY,
      case_data$F4B_TRT_GIVE_MACR, case_data$F4B_TRT_PRES_MACR,
      case_data$F4B_TRT_GIVE_NALID, case_data$F4B_TRT_PRES_NALID,
      case_data$F4B_TRT_GIVE_OTHR, case_data$F4B_TRT_PRES_OTHR
    ) == 1, na.rm = TRUE) > 0 # |
      # case_data$F4A_HOMETRT_AB_SPEC %in% c(
      #   "COTRIMOXAZOLE", "AMOXICILLIN", "AMXICILLINE+METRO", "AMOXI + ACIDE CLAVUL",
      #   "amoxil", "amoxyl", "AMOXYL", "CHLORAMPHENICOL", "TETRACYCLINE", "Amoxicillin",
      #   "AMOXACILLIN", "COTRIMOXAZOLLE", "COTROMOXAZOLE", "COTRIMOXAZOL", "COTIRMOXAZOLE",
      #   "COTROMOXAZOL", "SULFAMIDE", "SULFAMIDES", "SULFADIME", "sulfamide", "AMOXICILLINE",
      #   "ERYTHROMYCINE", "ERYTROMYCINE", "COTRIMXAZOLE", "NIFLUROXAZIDE", "doxy",
      #   "Doxycylline", "SEPTRIN, FLAGYL", "SEPTRIN", "SEPTRIN SYRUP", "SEPTRIN TABLET", "SEPTRIN  TABLET"
      # )
    , 1, 0
  )
  
  case_data$ineff_abx <- ifelse(
    rowSums(cbind(case_data$F4B_TRT_GIVE_PEN, case_data$F4B_TRT_PRES_PEN) == 1, na.rm = TRUE) > 0, # |
      # case_data$F4A_HOMETRT_AB_SPEC %in% c("METRONIDAZOLE","PHARMACY", "A SYRUP AND TABLETS", 
      #                                      "cloxacillin" ,"Metronidazole" , "PARACETAMOL","CEFADROXIL",
      #                                      "METRONODAZOLE", "ENTAMIZOLE","OREX", "FLAGYL", "ORACEFAL",
      #                                      "METRO PERFUSION","METRONIDAZOL","METRONIDAZOLE SYRUP" ,"flaggyl"),
    1, 0
  )
  
  # take highest of abx
  case_data$ineff_abx <- ifelse(case_data$ineff_abx == 1 & (case_data$who_abx == 1 | case_data$maybe_eff_abx == 1), 0, case_data$ineff_abx)
  case_data$maybe_eff_abx <- ifelse(case_data$who_abx == 1 & case_data$maybe_eff_abx == 1, 0, case_data$maybe_eff_abx)
  case_data$no_abx <- ifelse(case_data$ineff_abx == 0 & case_data$maybe_eff_abx == 0 & case_data$who_abx == 0, 1, 0)
  
  case_data$all_abx <- ifelse(case_data$ineff_abx == 1 | case_data$no_abx == 1, 0, 
                              ifelse(case_data$maybe_eff_abx == 1, 1, 2))
  
  case_data$all_abx <- factor(case_data$all_abx, levels = 0:2, labels = c("No or ineffective antibiotics",
                                                                          "Possibly effective antibiotics",
                                                                          "Guideline recommended antibiotics"))
  
  case_data$any_abx <- ifelse(
    rowSums(case_data[, c( # "F4A_HOMETRT_AB",
                          "F4B_TRT_PRES_CXL",
                          "F4B_TRT_GIVE_CXL",
                          "F4B_TRT_PRES_GENT",
                          "F4B_TRT_GIVE_GENT",
                          "F4B_TRT_PRES_CHLOR",
                          "F4B_TRT_GIVE_CHLOR",
                          "F4B_TRT_PRES_ERY",
                          "F4B_TRT_GIVE_ERY",
                          "F4B_TRT_PRES_AZI",
                          "F4B_TRT_GIVE_AZI",
                          "F4B_TRT_PRES_MACR",
                          "F4B_TRT_GIVE_MACR",
                          "F4B_TRT_PRES_PEN",
                          "F4B_TRT_GIVE_PEN",
                          "F4B_TRT_PRES_AMOX",
                          "F4B_TRT_GIVE_AMOX",
                          "F4B_TRT_PRES_CEF",
                          "F4B_TRT_GIVE_CEF",
                          "F4B_TRT_PRES_CIP",
                          "F4B_TRT_GIVE_CIP",
                          "F4B_TRT_PRES_AMPI",
                          "F4B_TRT_GIVE_AMPI",
                          "F4B_TRT_PRES_NALID",
                          "F4B_TRT_GIVE_NALID",
                          "F4B_TRT_PRES_CPNR",
                          "F4B_TRT_GIVE_CPNR",
                          "F4B_TRT_PRES_SLPY",
                          "F4B_TRT_GIVE_SLPY",
                          "F4B_TRT_PRES_FLAG",
                          "F4B_TRT_GIVE_FLAG",
                          "F4B_TRT_PRES_OTHR",
                          "F4B_TRT_GIVE_OTHR")] == 1, na.rm = TRUE) > 0, 1, 0)
  
  
  # New 9/4/25 -- incorporate Shigella resistance data ----
  
  # if any I/R, classify as I/R
  # if all S, classify as S
  
  # Ceftriaxone
  case_data$ceft_S_IR <- apply(
    case_data[, c("CEFTR_RIS_COLONY1", "CEFTR_RIS_COLONY2",
                  "CEFTR_RIS_COLONY3", "CEFTR_RIS_COLONY4")],
    1,
    function(x) {
      if (all(is.na(x))) {
        return(NA)
      } else if (any(x %in% c("I", "R"), na.rm = TRUE)) {
        return("I/R")
      } else if (all(x == "S", na.rm = TRUE)) {
        return("S")
      } else {
        # mixed case: some S, some NA, but no I/R
        return("S")   # or return(NA) if you want stricter handling
      }
    }
  )
  
  # Ciprofloxacin
  case_data$cipro_S_IR <- apply(
    case_data[, c("CIPRO_RIS_COLONY1", "CIPRO_RIS_COLONY2",
                  "CIPRO_RIS_COLONY3", "CIPRO_RIS_COLONY4")],
    1,
    function(x) {
      if (all(is.na(x))) {
        return(NA)
      } else if (any(x %in% c("I", "R"), na.rm = TRUE)) {
        return("I/R")
      } else if (all(x == "S", na.rm = TRUE)) {
        return("S")
      } else {
        # mixed case: some S, some NA, but no I/R
        return("S")  
      }
    }
  )
  
  
  # Azithromycin
  # case_data$azithro_S_IR <- apply(
  #   case_data[, c("AZITHRO_FLEX_COLONY1", "AZITHRO_FLEX_COLONY2",
  #                 "AZITHRO_FLEX_COLONY3", "AZITHRO_FLEX_COLONY4",
  #                 "AZITHRO_SONNEI_COLONY1", "AZITHRO_SONNEI_COLONY2",
  #                 "AZITHRO_SONNEI_COLONY3", "AZITHRO_SONNEI_COLONY4")],
  #   1,
  #   function(x) {
  #     if (all(is.na(x))) {
  #       return(NA)
  #     } else if (any(x %in% c("NWT"), na.rm = TRUE)) {
  #       return("I/R")
  #     } else if (all(x == "WT", na.rm = TRUE)) {
  #       return("S")
  #     } else {
  #       # mixed case: some WT, some NA, but no NWT
  #       return("S")   
  #     }
  #   }
  # )
  
  # create overall resistant and susceptible variables w same logic as EFGH
  # but revisit in meeting friday bc seems weird?? 
  case_data$resistant_WHO_approve <- apply(
    #case_data[, c("azithro_S_IR", "cipro_S_IR", "ceft_S_IR")],
    case_data[, c("cipro_S_IR", "ceft_S_IR")],
    1,
    function(x) if (any(x == "I/R", na.rm = TRUE)) 1 else 0
  )
  
  case_data$susceptible_WHO_approve <- apply(
    #case_data[, c("azithro_S_IR", "cipro_S_IR", "ceft_S_IR")],
    case_data[, c("cipro_S_IR", "ceft_S_IR")],
    1,
    function(x) if (any(x == "S", na.rm = TRUE)) 1 else 0
  )
  
  # Now create AST_given_abx variable for if they are resistant to an antibiotic they received
  # if susceptible for any abx received, mark as S
  
  # variables for S_IR if given drug
  #case_data$azithro_given_S_IR <- ifelse(case_data$azithro == 1, case_data$azithro_S_IR, NA)
  case_data$cipro_given_S_IR <- ifelse(case_data$cipro == 1, case_data$cipro_S_IR, NA)
  case_data$ceft_given_S_IR <- ifelse(case_data$ceft == 1, case_data$ceft_S_IR, NA)
  
  case_data$ast_given_abx <- apply(
    #case_data[, c("azithro_given_S_IR", "cipro_given_S_IR", "ceft_given_S_IR")],
    case_data[, c("cipro_given_S_IR", "ceft_given_S_IR")],
    1,
    function(x) {
      if (all(is.na(x))) {
        return(NA)
      } else if (any(x == "S", na.rm = TRUE)) {
        return("S")   # S dominates
      } else if (all(x == "I/R", na.rm = TRUE)) {
        return("I/R")
      } 
    }
  )
  
  # -------------------------------------------------------
  
  case_data$site <- factor(case_data$site, levels = c(1,2,3), labels = c("The Gambia", "Mali", "Kenya"))
  
  case_data$education <- ifelse(case_data$education == 7, NA, case_data$education)
  case_data$education <- factor(case_data$education, levels = 1:6, labels = c("None", "Less than primary", "Completed primary",
                                                                              "Completed secondary", "Post secondary", "Religious only"))
  
  case_data$education_bin <- ifelse(case_data$education %in% c("None", "Less than primary", "Completed primary", "Religious only"), 0, 1)
  
  case_data$safe_water <- factor(case_data$safe_water, levels = 1:5, labels = c("Safely managed",
                                                                                "Basic",
                                                                                "Limited",
                                                                                "Unimproved",
                                                                                "Surface water"))
  
  case_data$safe_sanit <- factor(case_data$safe_sanit, levels = 1:4, labels = c("Safely managed and basic",
                                                                                "Limited",
                                                                                "Unimproved",
                                                                                "Open defication"))
  
  case_data$sex <- factor(case_data$sex, levels = 1:2, labels = c("Male", "Female"))
  
  case_data$dehydr <- factor(case_data$dehydr, levels = c(0,1,2), labels = c("None", "Some dehydration", "Severe dehydration"))
  
  case_data$ses_quintile <- factor(case_data$ses_quintile, levels = 1:5, labels = c("1st quintile of SES",
                                                                                    "2nd quintile of SES",
                                                                                    "3rd quintile of SES",
                                                                                    "4th quintile of SES",
                                                                                    "5th quintile of SES"))
  
  case_data <- case_data %>%
    select(child_id,
           first_id,
           sex, 
           agemchild,
           agedchild,
           site,
           enr_haz,
           enr_waz,
           hazd60,
           hazdiff,
           education,
           education_bin,
           followup_days,
           I_followup_days,
           I_followup_days_x_followup_days,
           num_hh_lt5,
           bf_wk_before,
           lsstools,
           vom_days,
           vom_freq,
           fever, 
           dehydr,
           dysentery,
           duration_pre_enroll,
           safe_water,
           safe_sanit,
           ses_quintile,
           shigella_culture_positive,
           tac_shig,
           shigella_tac_or_culture,
           SHIG_FLEX,
           SHIG_SONNEI,
           rotavirus,
           st_etec,
           crypto,
           adeno,
           etec,
           astro,
           noro_gii,
           tepec,
           campy,
           campy_j,
           sapo,
           giardia,
           eaec,
           tac_adeno_4041,
           tac_astro,
           tac_campyj,
           tac_crypto,
           tac_cyclosporidium,
           tac_ehist,
           tac_norovirus_gii,
           tac_rota,
           tac_salm,
           tac_st_etec,
           tac_tepec,
           tac_vchol,
           no_etiology,
           rota_detected,
           adeno_4041_detected,
           etec_detected,
           crypto_detected,
           noro_detected,
           tepec_detected,
           campy_detected,
           sapo_detected,
           e_bieneusi_detected,
           giardia_detected,
           eaec_detected,
           shigella_new,
           rotavirus_new,
           st_etec_new,
           lt_etec_new, 
           etec_new,
           crypto_new,
           adeno_new,
           astro_new,
           noro_gii_new,
           tepec_new,
           campy_new,
           campy_j_new,
           sapo_new,
           giardia_new,
           eaec_new,
           e_bieneusi_new,
           v_cholerae_new,
           salmonella_new,
           who_abx,
           maybe_eff_abx,
           ineff_abx,
           no_abx,
           all_abx,
           any_abx,
           #azithro,
           cipro,
           ceft,
           #azithro_S_IR,
           cipro_S_IR,
           ceft_S_IR,
           ast_given_abx,
           resistant_WHO_approve,
           susceptible_WHO_approve) %>%
    rename("shig_flex" = SHIG_FLEX,
           "shig_sonnei" = SHIG_SONNEI) %>%
    set_variable_labels(child_id = "Child ID",
                        first_id = "First child ID",
                        sex = "Sex", 
                        agemchild = "Age (months)",
                        agedchild = "Age (days)",
                        site = "Study site",
                        enr_haz = "Baseline HAZ",
                        enr_waz = "Baseline WAZ",
                        hazd60 = "Day 60 HAZ",
                        hazdiff = "Change in HAZ enrollment - day 60",
                        education = "Primary caregiver education",
                        education_bin = "Primary caregiver education secondary school or greater",
                        followup_days = "Days between enrollment and follow-up",
                        num_hh_lt5 = "Number of children in household under age 5",
                        bf_wk_before = "Child breastfed week before illness",
                        lsstools = "Maximum number of loose stools worst day of episode",
                        vom_days = "Number of days child vomited",
                        vom_freq = "Times vomited on worst day",
                        fever = "Fever >38C", 
                        dehydr = "WHO dehydrated",
                        dysentery = "Dysentery",
                        safe_water = "JMP Improved Water" ,
                        safe_sanit = "JMP Improved Sanitation",
                        ses_quintile = "SES quintile by site",
                        shigella_culture_positive = "Shigella via culture",
                        tac_shig = "TAC attributable Shigella",
                        shigella_tac_or_culture = "TAC or culture positive Shigella",
                        rotavirus = "Rotavirus",
                        st_etec = "ST ETEC",
                        crypto = "Cryptosporidium",
                        adeno = "Adenovirus 40/41",
                        tac_rota = "TAC attributable rotavirus",
                        tac_st_etec = "TAC attributable ST ETEC",
                        tac_crypto = "TAC attributable cryptosporidium",
                        tac_adeno_4041 = "TAC attributable adenovirus 40/41",
                        who_abx = "WHO approved antibiotics",
                        maybe_eff_abx = "Maybe effective antibiotics",
                        ineff_abx = "Ineffective antibiotics",
                        no_abx = "No antibiotics",
                        all_abx = "Antibiotic treatment received",
                        any_abx = "Any antibiotics",
                        #azithro = "Azithromycin prescribed",
                        cipro = "Ciprofloxacin prescribed",
                        ceft = "Ceftriaxone prescribed",
                        #azithro_S_IR = "Shigella resistance to azithromycin",
                        cipro_S_IR = "Shigella resistance to ciprofloxacin",
                        ceft_S_IR = "Shigella resistance to ceftriaxone",
                        ast_given_abx = "Resistance to antibiotic prescribed",
                        resistant_WHO_approve = "Resistant to WHO approved abx",
                        susceptible_WHO_approve = "Susceptible to WHO approved abx")
  
  return(case_data)
  
}

vida_data_full <- prep_vida()
vida_data <- vida_data_full[which(!is.na(vida_data_full$tac_shig)), , drop = FALSE]

saveRDS(vida_data_full, here::here("data/vida_data/vida_data_full.Rds"))
saveRDS(vida_data, here::here("data/vida_data/vida_data.Rds"))

# -------------------------------------------------------------------
# CASE-CONTROL ANALYSIS
# -------------------------------------------------------------------

prep_vida_case_control <- function(case_def = "tac_or_culture_shig_diar"){
  
  data <- read.csv(here::here("data/vida_data/raw_data/VIDA_Rogawski.csv"))
  data_w_dss <- readxl::read_xlsx(here::here("data/vida_data/raw_data/VIDA_Rogawski_withDSScount and alt.xlsx"))
  
  # Just merge in dss_alt because excel formatting/var types is weird and messing things up, don't want fix
  data_w_dss <- data_w_dss[,c("CHILDID", "DSS_alt")]
  data <- left_join(data, data_w_dss, by = "CHILDID")
  
  # make first ID based on DSS_alt
  # group by dss alt, arrange by enroll_date, the first CHILDID == first id
  data <- data %>%
    group_by(DSS_alt) %>%
    arrange(DSS_alt, ENROLL_DATE) %>%
    mutate(first_id = first(CHILDID)) %>%
    ungroup()
  
  # -------------------------------------------------
  # Make PCA SES score in cases and controls by site 
  # -------------------------------------------------
  data$fridge <- ifelse(!is.na(data$F4A_HOUSE_FRIDGE), 
                        data$F4A_HOUSE_FRIDGE,
                        data$F7_HOUSE_FRIDGE)
  data$tv <- ifelse(!is.na(data$F4A_HOUSE_TELE),
                    data$F4A_HOUSE_TELE,
                    data$F7_HOUSE_TELE)
  data$electricity <- ifelse(!is.na(data$F4A_HOUSE_ELEC),
                             data$F4A_HOUSE_ELEC,
                             data$F7_HOUSE_ELEC)
  data$motorcycle <- ifelse(!is.na(data$F4A_HOUSE_SCOOT),
                            data$F4A_HOUSE_SCOOT,
                            data$F7_HOUSE_SCOOT)
  data$radio <- ifelse(!is.na(data$F4A_HOUSE_RADIO),
                       data$F4A_HOUSE_RADIO,
                       data$F7_HOUSE_RADIO)
  data$bike <- ifelse(!is.na(data$F4A_HOUSE_BIKE),
                      data$F4A_HOUSE_BIKE,
                      data$F7_HOUSE_BIKE)
  data$car <- ifelse(!is.na(data$F4A_HOUSE_CAR),
                     data$F4A_HOUSE_CAR,
                     data$F7_HOUSE_CAR)
  data$boat <- ifelse(!is.na(data$F4A_HOUSE_BOAT),
                      data$F4A_HOUSE_BOAT,
                      data$F7_HOUSE_BOAT)
  data$phone <- ifelse(!is.na(data$F4A_HOUSE_PHONE),
                       data$F4A_HOUSE_PHONE,
                       data$F7_HOUSE_PHONE)
  data$cart <- ifelse(!is.na(data$F4A_HOUSE_CART),
                      data$F4A_HOUSE_CART,
                      data$F7_HOUSE_CART)
  data$agland <- ifelse(!is.na(data$F4A_HOUSE_AGLAND),
                        data$F4A_HOUSE_AGLAND,
                        data$F7_HOUSE_AGLAND)
  
  cols_for_pca <- c("fridge", "tv", "electricity", "motorcycle",
                    "radio", "bike", "car", "boat", "phone", "cart", "agland")
  
  ses_pca <- prcomp(data[complete.cases(data[,cols_for_pca]), cols_for_pca], 
                    center = TRUE, scale = TRUE)
  
  pca_score <- as.matrix(data[, cols_for_pca]) %*% ses_pca$rotation[, 1]
  
  data$ses_score <- pca_score[, 1]   
  
  data <- data %>%
    group_by(SITE) %>%
    mutate(ses_quintile = ntile(ses_score, 5)) 
  
  data <- data %>%
    rename('child_id' = CHILDID,
           'case_id' = CASEID,
           'first_id' = first_id,
           'sex' = GENDER,
           'agemchild' = AGE,
           'agedchild' = BASE_AGE_DAYS,
           'site' = SITE,
           'enr_haz' = BASE_HAZ,
           'hazd60' = END_HAZ,
           'enr_whz' = BASE_WHZ,
           'whzd60' = END_WHZ,
           'enr_waz' = BASE_WAZ,
           'wazd60' = END_WAZ,
           'hazdiff' = DELTA_HAZ,
           'whzdiff' = DELTA_WHZ,
           'wazdiff' = DELTA_WAZ,
           'education' = F4A_PRIM_SCHL,
           'followup_days' = DURDAYS,
           'num_hh_lt5' = F4A_YNG_CHILDREN,
           'bf_wk_before' = F4A_PRI_BMILK,
           'lsstools' = F4A_DAILY_MAX,
           'vom_days' = F4A_DAYS_VOMIT,
           'vom_freq' = F4A_FREQ_VOMIT,
           'fever' = F4A_DRH_FEVER,
           'dehydr' = WHO_DEHYD,
           'dysentery' = DYS_IND,
           'safe_water' = IMPROV_WATER,     # water, sanitation, electricity proxy for SES
           'safe_sanit' = IMPROV_SANIT)
  
  data$vom_days <- ifelse(data$F4A_ANY_VOMIT == 1, data$vom_days, data$F4A_ANY_VOMIT)
  data$vom_freq <- ifelse(data$F4A_ANY_VOMIT == 1, data$vom_freq, data$F4A_ANY_VOMIT)
  
  # vomiting not unified well
  # cases has categories for 1, 2-4, 5 or more but controls is just 3 or more
  data$mult_vom <- ifelse(data$vom_freq >= 2, 1, 0)
  data$mult_vom <- ifelse(!is.na(data$F7_VOMIT), data$F7_VOMIT, data$mult_vom)
  
  # add in comparable vars for controls if applicable
  data$education <- ifelse(is.na(data$education), data$F7_PRIM_SCHL, data$education)
  data$num_hh_lt5 <- ifelse(is.na(data$num_hh_lt5), data$F7_YNG_CHILDREN, data$num_hh_lt5)
  data$bf_wk_before <- ifelse(is.na(data$bf_wk_before), data$F7_CUR_BMILK, data$bf_wk_before) # this one not exactly the same
  data$dysentery <- ifelse(is.na(data$dysentery), data$F7_BLOOD, data$dysentery)
  data$fever <- ifelse(is.na(data$fever), data$F7_FEVER, data$fever)
  
  data$I_followup_days <- ifelse(is.na(data$followup_days), 0 , 1)
  data$I_followup_days_x_followup_days <- ifelse(is.na(data$followup_days), 0 , data$followup_days)
  
  data$duration_pre_enroll <- data$F4A_DRH_DAYS
  
  # Duration post-enrollment -- from diarrhea dictionary
  # added 9/4/26
  # check with liz
  # only include people who fully completed? memory aid -- 1
  # data has 0 through 14, pic in manuscript has 1 - 14, so verify dates w pre-enroll 
  # pre-enroll includes enrollment date
  data$duration_post_enroll <- ifelse(!is.na(data$F9_DRH_LAST) & data$F9_MEMORY_AID == 1, data$F9_DRH_LAST, NA)
  
  # Death -- in facility, at 60 day visit, overall
  # already made in gems as follows:
  # data_full$death = overall (in facility + 60 day visit)
  # data_full$died60 = 60 day visit
  # data_full$diedhosp = died in hospital before discharge
  
  data$death <- data$DEATH_IND
  data$diedhosp <- ifelse(data$F4B_OUTCOME == 5, 1, 0)
  data$died60 <- ifelse(data$death == 1 & data$diedhosp == 0, 1, 0) # don't have F5 var but based on dictionary can work backwards
  
  # can't make visit health facility during follow up bc missing F5 vars
  # Define Shigella
  data$tac_shig <- ifelse(data$TAC_SHIGELLA_EIEC < 27.9, 1, 0)
  
  # Only look at people who have tac results
  # stopped doing for co-etiology meta-analysis
  # data <- data[-which(is.na(data$tac_shig)),]
  
  data$shigella_culture_positive <- ifelse(rowSums(data[, c("F16_SHIG_SPP", 
                                                            "SHIG_BOYDII", 
                                                            "SHIG_DYSENT", 
                                                            "SHIG_FLEX", 
                                                            "SHIG_SONNEI")] == 1, na.rm = TRUE) > 0, 1, 0)
  
  data$shigella_tac_or_culture <- ifelse(
    rowSums(data[, c("shigella_culture_positive", "tac_shig")] == 1, na.rm = TRUE) > 0,
    1, 0
  )
  
  if(case_def == "tac_or_culture_shig_diar"){
    
    # Only look at people who have tac results
    # stopped doing for co-etiology meta-analysis
    data <- data[-which(is.na(data$tac_shig)),]
    
    data$case <- ifelse(data$TYPE == 1, 1, 0)
    
    # eliminate people who are diarrhea cases but negative for shigella (or NA result)
    data <- data[-which(data$case == 1 & (data$shigella_tac_or_culture == 0 | is.na(data$shigella_tac_or_culture))),]
    
    
    case_ids <- unique(data$child_id[which(data$case == 1)])
    data <- data[-which(data$case == 0 & !(data$case_id %in% case_ids)),]
    
    
  } else if(case_def == "tac_shig"){
    
    # Only look at people who have tac results
    # stopped doing for co-etiology meta-analysis
    data <- data[-which(is.na(data$tac_shig)),]
    
    data$case <- ifelse(data$TYPE == 1, 1, 0)
    
    data <- data[-which(data$case == 1 & (data$tac_shig == 0 | is.na(data$tac_shig))),]
    
    case_ids <- unique(data$child_id[which(data$case == 1)])
    data <- data[-which(data$case == 0 & !(data$case_id %in% case_ids)),]
    
  }else if(case_def == "culture_shig_diar"){
    
    # Only look at people who have tac results
    # stopped doing for co-etiology meta-analysis
    data <- data[-which(is.na(data$tac_shig)),]
    
    # we want case definition to be culture positive, but still subset to people who have 
    # shigella tac results ^^
    
    data$case <- ifelse(data$TYPE == 1, 1, 0)
    
    data <- data[-which(data$case == 1 & (data$shigella_culture_positive == 0 | is.na(data$shigella_culture_positive))),]
    
    case_ids <- unique(data$child_id[which(data$case == 1)])
    data <- data[-which(data$case == 0 & !(data$case_id %in% case_ids)),]
    
  } else{
    data$case <- ifelse(data$TYPE == 1, 1, 0)
    
    # all diarrhea cases and controls as is
  }
  
  # Covariates:
  # sex, age, site, bl_haz, edu, bf week prior, days, additional pathogens, kids under 5 in hh
  # SES = improved drinking, improved sanitation, electricity
  
  # Severity: 
  # dysentery, fever, dehydration, loose stools, vomit
  
  # Get rid of extreme HAZ observations
  data$enr_haz <- ifelse(data$enr_haz < -6 | data$enr_haz > 6, NA, data$enr_haz)
  data$hazd60 <- ifelse(data$hazd60 < -6 | data$hazd60 > 6, NA, data$hazd60)
  
  data$enr_whz <- ifelse(data$enr_whz < -6 | data$enr_whz > 6, NA, data$enr_whz)
  data$whzd60 <- ifelse(data$whzd60 < -6 | data$whzd60 > 6, NA, data$whzd60)
  
  data$enr_waz <- ifelse(data$enr_waz < -6 | data$enr_waz > 6, NA, data$enr_waz)
  data$wazd60 <- ifelse(data$wazd60 < -6 | data$wazd60 > 6, NA, data$wazd60)
  
  data$hazdiff <- ifelse(is.na(data$enr_haz) | is.na(data$hazd60), NA, data$hazdiff)
  data$whzdiff <- ifelse(is.na(data$enr_whz) | is.na(data$whzd60), NA, data$whzdiff)
  data$wazdiff <- ifelse(is.na(data$enr_waz) | is.na(data$wazd60), NA, data$wazdiff)
  
  data$vom_days <- ifelse(is.na(data$vom_days), 0, data$vom_days)
  data$vom_freq <- ifelse(is.na(data$vom_freq), 0, data$vom_freq)
  
  # Other Pathogens
  
  data$rotavirus <- data$F18_RES_ROTAVIRUS
  data$st_etec <- data$ETEC_ST
  data$crypto <- data$F18_RES_CRYPTOSPOR
  data$adeno <- data$ADENO_4041
  data$etec <- data$ETEC_ALL
  data$astro <- data$F19_ASTRO_VIRUS
  data$noro <- data$NORO_ANY
  data$tepec <- data$TEPEC
  data$campy <- data$CAMPY_ANY
  data$sapo <- data$F19_SAPO_VIRUS
  data$giardia <- data$F18_RES_GIARDIA
  data$eaec <- data$EAEC
  
  # TAC using cutoffs from GEMS - https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(16)31529-X/abstract
  data$tac_st_etec <- ifelse(data$TAC_ST_ETEC < 22.8, 1, 0)
  data$tac_rota <- ifelse(data$TAC_ROTAVIRUS < 32.6, 1, 0)
  data$tac_crypto <- ifelse(data$TAC_CRYPTO < 24.0, 1, 0)
  data$tac_adeno_4041 <- ifelse(data$TAC_ADENO4041 < 22.7, 1, 0)
  data$tac_astro <- ifelse(data$TAC_ASTROVIRUS < 22.2, 1, 0)
  data$tac_campyj <- ifelse(data$TAC_CAMPY < 15.4, 1, 0)
  data$tac_cyclosporidium <- ifelse(data$TAC_CYCLO < 29.6, 1, 0)
  data$tac_ehist <- ifelse(data$TAC_E_HISTOLYTICA < 32.8, 1, 0)
  data$tac_norovirus_gii <- ifelse(data$TAC_NORO_GII < 23.4, 1, 0)
  data$tac_salm <- ifelse(data$TAC_SALMONELLA < 30.7, 1, 0)
  data$tac_tepec <- ifelse(data$TAC_TEPEC < 16.0, 1, 0)
  data$tac_vchol <- ifelse(data$TAC_V_CHOLERAE < 33.8, 1, 0)
  
  # detected
  data$rota_detected <- ifelse(data$TAC_ROTAVIRUS < 35, 1, 0)
  data$adeno_4041_detected <- ifelse(data$TAC_ADENO4041 < 35, 1, 0)
  data$etec_detected <- ifelse(data$TAC_ST_ETEC < 35 | data$TAC_LT_ETEC < 35, 1, 0)
  data$crypto_detected <- ifelse(data$TAC_CRYPTO < 35, 1, 0)
  data$astro_detected <- ifelse(data$TAC_ASTROVIRUS < 35, 1, 0)
  data$noro_detected <- ifelse(data$TAC_NORO_GI < 35 | data$TAC_NORO_GII < 35, 1, 0)
  data$tepec_detected <- ifelse(data$TAC_TEPEC < 35, 1, 0)
  data$campy_detected <- ifelse(data$TAC_CAMPY_ANY < 35, 1, 0)
  data$sapo_detected <- ifelse(data$TAC_SAPOVIRUS < 35, 1, 0)
  data$e_bieneusi_detected <- ifelse(data$TAC_E_BIENEUSI < 35, 1, 0)
  data$giardia_detected <- ifelse(data$TAC_GIARDIA < 35, 1, 0)
  data$eaec_detected <- ifelse(data$TAC_EAEC < 35, 1, 0)
  
  # Create 'new' scaled variables
  data$shigella_new <- (35 - data$TAC_SHIGELLA_EIEC) / 3.322
  data$rotavirus_new <- (35 - data$TAC_ROTAVIRUS) / 3.322
  data$st_etec_new <- (35 - data$TAC_ST_ETEC) / 3.322
  data$lt_etec_new <- (35 - data$TAC_LT_ETEC) / 3.322
  data$etec_new <- ifelse(data$st_etec_new > data$lt_etec_new, data$st_etec_new, data$lt_etec_new)
  data$crypto_new <- (35 - data$TAC_CRYPTO) / 3.322
  data$adeno_new <- (35 - data$TAC_ADENO4041) / 3.322
  # data$etec <- (35 - data$TAC_ST_ETEC) / 3.322
  data$astro_new <- (35 - data$TAC_ASTROVIRUS) / 3.322
  data$noro_gii_new <- (35 - data$TAC_NORO_GII) / 3.322
  data$tepec_new <- (35 - data$TAC_TEPEC) / 3.322
  data$campy_new <- (35 - data$TAC_CAMPY_ANY) / 3.322
  data$campy_j_new <- (35 - data$TAC_CAMPY) / 3.322
  data$sapo_new <- (35 - data$TAC_SAPOVIRUS) / 3.322
  data$giardia_new <- (35 - data$TAC_GIARDIA) / 3.322
  data$eaec_new <- (35 - data$TAC_EAEC) / 3.322
  data$e_bieneusi_new <- (35 - data$TAC_E_BIENEUSI) / 3.322
  
  data$v_cholerae_new <- (35 - data$TAC_V_CHOLERAE) / 3.322
  data$salmonella_new <- (35 - data$TAC_SALMONELLA) / 3.322
  
  # Antibiotics
  # WHO approved = azithromycin, ciprofloxacin, ceftriaxone, pivmecillinam
  data$who_abx <- ifelse(
    rowSums(cbind(
      data$F4B_TRT_GIVE_AZI, data$F4B_TRT_PRES_AZI,
      data$F4B_TRT_GIVE_CPNR, data$F4B_TRT_PRES_CPNR,
      data$F4B_TRT_GIVE_CEF, data$F4B_TRT_PRES_CEF,
      data$F4B_TRT_GIVE_SLPY, data$F4B_TRT_PRES_SLPY
    ) == 1, na.rm = TRUE) > 0, # |
      # data$F4A_HOMETRT_AB_SPEC %in% c(
      #   "CEFTRIAXONE", "CEFTRIAXONE 500mg", "ciprofloxacine", "CECTRIAXONE"
      # ), 
    1, 0
  )
  
  data$azithro <- ifelse(
    rowSums(cbind(
      data$F4B_TRT_GIVE_AZI, data$F4B_TRT_PRES_AZI
    ) == 1, na.rm = TRUE) > 0,
    1, 0
  )
  
  data$cipro <- ifelse(
    rowSums(cbind(
      data$F4B_TRT_GIVE_CPNR, data$F4B_TRT_PRES_CPNR
    ) == 1, na.rm = TRUE) > 0, 
    1, 0
  )
  
  data$ceft <- ifelse(
    rowSums(cbind(
      data$F4B_TRT_GIVE_CEF, data$F4B_TRT_PRES_CEF
    ) == 1, na.rm = TRUE) > 0,
    1, 0
  )
  
  data$maybe_eff_abx <- ifelse(
    rowSums(cbind(
      data$F4B_TRT_GIVE_AMOX, data$F4B_TRT_PRES_AMOX,
      data$F4B_TRT_GIVE_AMPI, data$F4B_TRT_PRES_AMPI,
      data$F4B_TRT_GIVE_CXL, data$F4B_TRT_PRES_CXL,
      data$F4B_TRT_GIVE_GENT, data$F4B_TRT_PRES_GENT,
      data$F4B_TRT_GIVE_CHLOR, data$F4B_TRT_PRES_CHLOR,
      data$F4B_TRT_GIVE_ERY, data$F4B_TRT_PRES_ERY,
      data$F4B_TRT_GIVE_MACR, data$F4B_TRT_PRES_MACR,
      data$F4B_TRT_GIVE_NALID, data$F4B_TRT_PRES_NALID,
      data$F4B_TRT_GIVE_OTHR, data$F4B_TRT_PRES_OTHR
    ) == 1, na.rm = TRUE) > 0, # |
      # data$F4A_HOMETRT_AB_SPEC %in% c(
      #   "COTRIMOXAZOLE", "AMOXICILLIN", "AMXICILLINE+METRO", "AMOXI + ACIDE CLAVUL",
      #   "amoxil", "amoxyl", "AMOXYL", "CHLORAMPHENICOL", "TETRACYCLINE", "Amoxicillin",
      #   "AMOXACILLIN", "COTRIMOXAZOLLE", "COTROMOXAZOLE", "COTRIMOXAZOL", "COTIRMOXAZOLE",
      #   "COTROMOXAZOL", "SULFAMIDE", "SULFAMIDES", "SULFADIME", "sulfamide", "AMOXICILLINE",
      #   "ERYTHROMYCINE", "ERYTROMYCINE", "COTRIMXAZOLE", "NIFLUROXAZIDE", "doxy",
      #   "Doxycylline", "SEPTRIN, FLAGYL", "SEPTRIN", "SEPTRIN SYRUP", "SEPTRIN TABLET"
      # ), 
    1, 0
  )
  
  data$ineff_abx <- ifelse(
    rowSums(cbind(data$F4B_TRT_GIVE_PEN, data$F4B_TRT_PRES_PEN) == 1, na.rm = TRUE) > 0, # |
      # data$F4A_HOMETRT_AB_SPEC %in% c("METRONIDAZOLE","PHARMACY", "A SYRUP AND TABLETS", 
      #                                 "cloxacillin" ,"Metronidazole" , "PARACETAMOL","CEFADROXIL",
      #                                 "METRONODAZOLE", "ENTAMIZOLE","OREX", "FLAGYL", "ORACEFAL",
      #                                 "METRO PERFUSION","METRONIDAZOL","METRONIDAZOLE SYRUP" ,"flaggyl"),
    1, 0
  )
  
  # take highest of abx
  data$ineff_abx <- ifelse(data$ineff_abx == 1 & (data$who_abx == 1 | data$maybe_eff_abx == 1), 0, data$ineff_abx)
  data$maybe_eff_abx <- ifelse(data$who_abx == 1 & data$maybe_eff_abx == 1, 0, data$maybe_eff_abx)
  data$no_abx <- ifelse(data$ineff_abx == 0 & data$maybe_eff_abx == 0 & data$who_abx == 0, 1, 0)
  
  data$all_abx <- ifelse(data$ineff_abx == 1 | data$no_abx == 1, 0, 
                         ifelse(data$maybe_eff_abx == 1, 1, 2))
  
  data$all_abx <- factor(data$all_abx, levels = 0:2, labels = c("Ineffective or no abx", "Maybe effective abx", "WHO approved abx"))
  
  data$any_abx <- ifelse(
    rowSums(data[, c(# "F4A_HOMETRT_AB",
                     "F4B_TRT_PRES_CXL",
                     "F4B_TRT_GIVE_CXL",
                     "F4B_TRT_PRES_GENT",
                     "F4B_TRT_GIVE_GENT",
                     "F4B_TRT_PRES_CHLOR",
                     "F4B_TRT_GIVE_CHLOR",
                     "F4B_TRT_PRES_ERY",
                     "F4B_TRT_GIVE_ERY",
                     "F4B_TRT_PRES_AZI",
                     "F4B_TRT_GIVE_AZI",
                     "F4B_TRT_PRES_MACR",
                     "F4B_TRT_GIVE_MACR",
                     "F4B_TRT_PRES_PEN",
                     "F4B_TRT_GIVE_PEN",
                     "F4B_TRT_PRES_AMOX",
                     "F4B_TRT_GIVE_AMOX",
                     "F4B_TRT_PRES_CEF",
                     "F4B_TRT_GIVE_CEF",
                     "F4B_TRT_PRES_CIP",
                     "F4B_TRT_GIVE_CIP",
                     "F4B_TRT_PRES_AMPI",
                     "F4B_TRT_GIVE_AMPI",
                     "F4B_TRT_PRES_NALID",
                     "F4B_TRT_GIVE_NALID",
                     "F4B_TRT_PRES_CPNR",
                     "F4B_TRT_GIVE_CPNR",
                     "F4B_TRT_PRES_SLPY",
                     "F4B_TRT_GIVE_SLPY",
                     "F4B_TRT_PRES_FLAG",
                     "F4B_TRT_GIVE_FLAG",
                     "F4B_TRT_PRES_OTHR",
                     "F4B_TRT_GIVE_OTHR")] == 1, na.rm = TRUE) > 0, 1, 0)
  
  # also add abx for controls so we can eliminate controls who are taking abx
  
  # Any antibiotics for controls
  data$any_abx <- ifelse(
    rowSums(
      cbind(
        data$F7_MED_GENT, data$F7_MED_COTR, data$F7_MED_CHLOR,
        data$F7_MED_ERYTH, data$F7_MED_AZITH, data$F7_MED_OMACR,
        data$F7_MED_PENI, data$F7_MED_AMOXY, data$F7_MED_AMPI,
        data$F7_MED_NALID, data$F7_MED_CIPRO, data$F7_MED_SELE,
        data$F7_MED_OTHERANT
      ) == 1, na.rm = TRUE
    ) > 0, 1, data$any_abx
  )
  
  # get rid of controls with abx, dysentery
  data <- data[-which(data$case == 0 & data$any_abx == 1),]
  data <- data[-which(data$case == 0 & data$dysentery == 1),]
  
  data$site <- factor(data$site, levels = c(1,2,3), labels = c("The Gambia", "Mali", "Kenya"))
  
  data$education <- ifelse(data$education == 7, NA, data$education)
  data$education <- factor(data$education, levels = 1:6, labels = c("None", "Less than primary", "Completed primary",
                                                                    "Completed secondary", "Post secondary", "Religious only"))
  
  data$education_bin <- ifelse(data$education %in% c("None", "Less than primary", "Completed primary", "Religious only"), 0, 1)
  
  data$safe_water <- factor(data$safe_water, levels = 1:5, labels = c("Safely managed",
                                                                      "Basic",
                                                                      "Limited",
                                                                      "Unimproved",
                                                                      "Surface water"))
  
  data$safe_sanit <- factor(data$safe_sanit, levels = 1:4, labels = c("Safely managed and basic",
                                                                      "Limited",
                                                                      "Unimproved",
                                                                      "Open defication"))
  
  data$sex <- factor(data$sex, levels = 1:2, labels = c("Male", "Female"))
  
  data$dehydr <- factor(data$dehydr, levels = c(0,1,2), labels = c("None", "Some dehydration", "Severe dehydration"))
  
  data$vom_freq <- factor(data$vom_freq, levels = c(0,1,2,3), labels = c("None", "One", "Two to four", "Five or more"))
  
  data$ses_quintile <- factor(data$ses_quintile, levels = 1:5, labels = c("1st quintile of SES",
                                                                          "2nd quintile of SES",
                                                                          "3rd quintile of SES",
                                                                          "4th quintile of SES",
                                                                          "5th quintile of SES"))
  data <- data %>%
    select(child_id,
           first_id,
           case_id,
           case,
           shigella_tac_or_culture,
           sex, 
           agemchild,
           agedchild,
           site,
           ses_quintile,
           enr_haz,
           hazd60,
           hazdiff,
           enr_whz,
           whzd60,
           whzdiff,
           enr_waz,
           wazd60,
           wazdiff,
           education,
           education_bin,
           followup_days,
           I_followup_days,
           I_followup_days_x_followup_days,
           num_hh_lt5,
           bf_wk_before,
           lsstools,
           vom_days,
           vom_freq,
           mult_vom,
           fever, 
           dehydr,
           dysentery,
           duration_pre_enroll,
           duration_post_enroll,
           safe_water,
           safe_sanit,
           shigella_culture_positive,
           SHIG_FLEX,
           SHIG_SONNEI,
           rotavirus,
           st_etec,
           crypto,
           adeno,
           etec,
           astro,
           noro,
           tepec,
           campy,
           sapo,
           giardia,
           eaec,
           tac_shig,
           tac_rota,
           tac_st_etec,
           tac_crypto,
           tac_adeno_4041,
           tac_norovirus_gii,
           tac_astro,
           tac_tepec,
           tac_crypto,
           tac_vchol,
           tac_salm,
           tac_campyj,
           rota_detected,
           adeno_4041_detected,
           etec_detected,
           crypto_detected,
           astro_detected,
           noro_detected,
           tepec_detected,
           campy_detected,
           sapo_detected,
           e_bieneusi_detected,
           giardia_detected,
           eaec_detected,
           who_abx,
           shigella_new,
           rotavirus_new,
           st_etec_new,
           crypto_new,
           adeno_new,
           astro_new,
           noro_gii_new,
           tepec_new,
           campy_new,
           campy_j_new,
           sapo_new,
           etec_new,
           giardia_new,
           e_bieneusi_new,
           eaec_new,
           v_cholerae_new,
           salmonella_new,
           death,
           died60,
           diedhosp,
           maybe_eff_abx,
           ineff_abx,
           no_abx,
           all_abx,
           any_abx,
           azithro,
           cipro,
           ceft) %>%
    rename("shig_flex" = SHIG_FLEX,
           "shig_sonnei" = SHIG_SONNEI) %>%
    set_variable_labels(child_id = "Child ID",
                        first_id = "First ID",
                        case_id = "Case ID",
                        case = "Case",
                        sex = "Sex", 
                        agemchild = "Age (months)",
                        agedchild = "Age (days)",
                        site = "Study site",
                        enr_haz = "Baseline HAZ",
                        hazd60 = "Day 60 HAZ",
                        hazdiff = "Change in HAZ enrollment - day 60",
                        education = "Primary caregiver education",
                        education_bin = "Primary caregiver education secondary school or greater ",
                        ses_quintile = "SES quintile",
                        followup_days = "Days between enrollment and follow-up",
                        num_hh_lt5 = "Number of children in household under age 5",
                        bf_wk_before = "Case: Child breastfed week before illness | Control: Normal diet includes breast milk",
                        lsstools = "Maximum number of loose stools worst day of episode",
                        vom_days = "Number of days child vomited",
                        vom_freq = "Times vomited on worst day",
                        mult_vom = "Cases: Vomitted more than once | Controls: Vomitted more than twice",
                        fever = "Cases: Fever >38C | Control: Fever >38C week before", 
                        dehydr = "WHO dehydrated",
                        dysentery = "Cases: Dysentery | Control: NA (exclude controls with dysentery)",
                        safe_water = "JMP Improved Water" ,
                        safe_sanit = "JMP Improved Sanitation",
                        ses_quintile = 'SES quintile',
                        shigella_culture_positive = "Shigella culture",
                        rotavirus = "Rotavirus",
                        st_etec = "ST ETEC",
                        crypto = "Cryptosporidium",
                        adeno = "Adenovirus 40/41",
                        etec = "ETEC",
                        astro = "Astrovirus",
                        noro = "Norovirus",
                        tepec = "tEPEC",
                        campy = "Campylobacter",
                        sapo = "Sapovirus",
                        giardia = "Giardia",
                        eaec = "EAEC",
                        tac_shig = "TAC attributable Shigella",
                        tac_rota = "TAC attributable rotavirus",
                        tac_st_etec = "TAC attributable ST ETEC",
                        tac_crypto = "TAC attributable cryptosporidium",
                        tac_adeno_4041 = "TAC attributable adenovirus 40/41",
                        rota_detected = "Rotavirus detected",
                        adeno_4041_detected = "Adenovirus detected",
                        etec_detected = "ETEC detected",
                        crypto_detected = "Cryptosporidium detected",
                        astro_detected = "Astrovirus detected",
                        noro_detected = "Norovirus detected",
                        tepec_detected = "tEPEC detected",
                        campy_detected = "Campylobacter detected",
                        sapo_detected = "Sapovirus detected",
                        e_bieneusi_detected = "E bieneusi detected",
                        giardia_detected = "Giardia detected",
                        eaec_detected = "EAEC detected",
                        who_abx = "WHO approved antibiotics",
                        maybe_eff_abx = "Maybe effective antibiotics",
                        ineff_abx = "Ineffective or no antibiotics",
                        any_abx = "Any antibiotics",
                        azithro = "Received azithromycin",
                        cipro = "Received ciprofloxacin",
                        ceft = "Received ceftriaxone")
  
  return(data)
  
}

vida_case_control_tac_or_culture_shig <- prep_vida_case_control("tac_or_culture_shig_diar")
vida_case_control_tac_shig <- prep_vida_case_control("tac_shig")
vida_case_control_culture_shig <- prep_vida_case_control("culture_shig_diar")
vida_case_control_all_diar <- prep_vida_case_control("all_diar")

saveRDS(vida_case_control_tac_or_culture_shig, here::here("data/vida_data/vida_case_control_tac_or_culture_shig.Rds"))
saveRDS(vida_case_control_tac_shig, here::here("data/vida_data/vida_case_control_tac_shig.Rds"))
saveRDS(vida_case_control_culture_shig, here::here("data/vida_data/vida_case_control_culture_shig.Rds"))
saveRDS(vida_case_control_all_diar, here::here("data/vida_data/vida_case_control_all_diar.Rds"))

