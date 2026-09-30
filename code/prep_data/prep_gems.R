# ---------------------------------------------------------------
# Script to create cleaned GEMS datasets
# ---------------------------------------------------------------

here::i_am("code/prep_data/prep_gems.R")

library(tidyverse)
library(labelled)
library(haven)
library(corrr)
library(factoextra)

# --------------------------------------------------------------------
# CASE-ONLY DATA PREP 
# --------------------------------------------------------------------

#' Function to prep GEMS case-only data
#' 
#' By default, this will return all diarrhea episodes (modified for co-etiology 
#' meta-analysis). For the original meta-analysis, we define case as TAC 
#' or culture attributable Shigella so subset the dataset from this function to those who
#' have Shigella TAC results.
prep_gems <- function(){
  
  # Load data
  case_control_data <- read.csv(here::here("data/gems_data/raw_data/GEMS-1 case-control dataset.csv"))
  load(here::here("data/gems_data/raw_data/GEMS with AFEs.RData")) # called 'data'
  amr <- haven::read_sas(here::here("data/gems_data/raw_data/amr.sas7bdat")) 
  
  # -------------------------------------------------
  # Make PCA SES score in cases and controls by site 
  # -------------------------------------------------
  case_control_data$fridge <- ifelse(!is.na(case_control_data$F4A_HOUSE_FRIDGE), 
                                     case_control_data$F4A_HOUSE_FRIDGE,
                                     case_control_data$F7_HOUSE_FRIDGE)
  case_control_data$tv <- ifelse(!is.na(case_control_data$F4A_HOUSE_TELE),
                                 case_control_data$F4A_HOUSE_TELE,
                                 case_control_data$F7_HOUSE_TELE)
  case_control_data$electricity <- ifelse(!is.na(case_control_data$F4A_HOUSE_ELEC),
                                          case_control_data$F4A_HOUSE_ELEC,
                                          case_control_data$F7_HOUSE_ELEC)
  case_control_data$motorcycle <- ifelse(!is.na(case_control_data$F4A_HOUSE_SCOOT),
                                         case_control_data$F4A_HOUSE_SCOOT,
                                         case_control_data$F7_HOUSE_SCOOT)
  case_control_data$radio <- ifelse(!is.na(case_control_data$F4A_HOUSE_RADIO),
                                    case_control_data$F4A_HOUSE_RADIO,
                                    case_control_data$F7_HOUSE_RADIO)
  case_control_data$bike <- ifelse(!is.na(case_control_data$F4A_HOUSE_BIKE),
                                   case_control_data$F4A_HOUSE_BIKE,
                                   case_control_data$F7_HOUSE_BIKE)
  case_control_data$car <- ifelse(!is.na(case_control_data$F4A_HOUSE_CAR),
                                  case_control_data$F4A_HOUSE_CAR,
                                  case_control_data$F7_HOUSE_CAR)
  case_control_data$boat <- ifelse(!is.na(case_control_data$F4A_HOUSE_BOAT),
                                   case_control_data$F4A_HOUSE_BOAT,
                                   case_control_data$F7_HOUSE_BOAT)
  case_control_data$phone <- ifelse(!is.na(case_control_data$F4A_HOUSE_PHONE),
                                    case_control_data$F4A_HOUSE_PHONE,
                                    case_control_data$F7_HOUSE_PHONE)
  case_control_data$cart <- ifelse(!is.na(case_control_data$F4A_HOUSE_CART),
                                   case_control_data$F4A_HOUSE_CART,
                                   case_control_data$F7_HOUSE_CART)
  case_control_data$agland <- ifelse(!is.na(case_control_data$F4A_HOUSE_AGLAND),
                                     case_control_data$F4A_HOUSE_AGLAND,
                                     case_control_data$F7_HOUSE_AGLAND)
  
  cols_for_pca <- c("fridge", "tv", "electricity", "motorcycle",
                    "radio", "bike", "car", "boat", "phone", "cart", "agland")
  
  ses_pca <- prcomp(case_control_data[complete.cases(case_control_data[,cols_for_pca]), cols_for_pca], 
                    center = TRUE, scale = TRUE)
  
  pca_score <- as.matrix(case_control_data[, cols_for_pca]) %*% ses_pca$rotation[, 1]
  
  case_control_data$ses_score <- pca_score[, 1] 
  
  case_control_data <- case_control_data %>%
    group_by(SITE) %>%
    mutate(ses_quintile = ntile(desc(ses_score), 5)) # flip so lower = worse, higher = best
  
  # -------------------------------------------------------------------------
  
  # subset to cases 
  afe_cases <- data[data$case.control == 1,]
  case_data <- case_control_data[case_control_data$Type == "Case",]
  
  # with TAC data only
  # note some people only in one or the other dataset -- look into later? 
  # case_data_full <- merge(case_data, afe_cases, 
  #                         by.x = "FIRSTID", 
  #                         by.y = "childid")
  
  # case_control_data CHILDID all end in 1 and do not match data, get rid of
  case_data$pid <- as.numeric(substr(as.character(case_data$CHILDID), 
                                     1, 
                                     nchar(as.character(case_data$CHILDID)) - 1))
  
  case_data$child_id <- case_data$pid
  
  case_data$first_id <- as.numeric(case_data$FIRSTID)
  
  case_data_full <- left_join(case_data, afe_cases,
                              by = c("child_id" = "childid"),
                              keep = FALSE)
  
  ## OUTCOME - Change in HAZ enrollment --> day 60
  # case_data_full$hazdiff <- case_data_full$F4B_MED_HAZ - case_data_full$F5_HAZ
  # day 60 = F5_HAZ
  # enrollment = F4B_MED_HAZ
  
  # Eliminate HAZ > 6 or < -6
  case_data_full$F5_HAZ <- ifelse(case_data_full$F5_HAZ > 6 | case_data_full$F5_HAZ  < -6, NA, case_data_full$F5_HAZ)
  case_data_full$F4B_MED_HAZ <- ifelse(case_data_full$F4B_MED_HAZ > 6 | case_data_full$F4B_MED_HAZ  < -6, NA, case_data_full$F4B_MED_HAZ)
  
  # Same with WAZ
  case_data_full$F4B_WAZ <- ifelse(case_data_full$F4B_WAZ < -6 | case_data_full$F4B_WAZ > 6, NA, case_data_full$F4B_WAZ)
  
  case_data_full$hazdiff <- case_data_full$F5_HAZ - case_data_full$F4B_MED_HAZ
  case_data_full$hazd60 <- case_data_full$F5_HAZ
  
  ## SHIGELLA - Shigella attributable diarrhea by AFE
  case_data_full$shigella_attributable_tac <- ifelse(case_data_full$shigella_eiec_afe > 0.5, 1, 0)
  
  ## Shigella attributable diarrhea by culture
  case_data_full$shigella_attributable_culture <- ifelse(
    rowSums(case_data[, c(
      "F16CORR_S_DYSENT", "F16CORR_S_FLEXNERI", "F16CORR_S_BOYDII", 
      "F16CORR_S_SONNEI", "F16CORR_S_NONTYPE", "F16CORR_SHIG_1A", 
      "F16CORR_SHIG_1B", "F16CORR_SHIG_2A", "F16CORR_SHIG_2B", 
      "F16CORR_SHIG_3A", "F16CORR_SHIG_3B", "F16CORR_SHIG_4A", 
      "F16CORR_SHIG_4B", "F16CORR_SHIG_4C", "F16CORR_SHIG_5A", 
      "F16CORR_SHIG_5B", "F16CORR_SHIG_6", "F16CORR_SHIG_X", 
      "F16CORR_SHIG_Y", "F16CORR_SHIG_NONTYP"
    )] == 1, na.rm = TRUE) > 0, 1, 0)
  
  # If NA tac results, make sure it is NA
  case_data_full$shigella_attributable <- ifelse(
    case_data_full$shigella_attributable_tac == 1 | case_data_full$shigella_attributable_culture == 1,
    1,
    0
  )
  
  case_data_full$shigella_attributable <- ifelse(is.na(case_data_full$shigella_attributable_tac), NA, case_data_full$shigella_attributable)
  
  # Get rid of people without Shigella TAC results
  # ELIMINATED FOR CO-ETIOLOGY ANALYSIS 2/11/26
  # case_data_full <- case_data_full[!is.na(case_data_full$shigella_attributable), ]
  
  ## ANTIBIOTICS - vars for any abx, who approved abx
  case_data_full$any_abx <- 0
  
  case_data_full$who_abx <- 0
  case_data_full$azithro <- 0
  case_data_full$maybe_eff_abx <- 0
  case_data_full$ineff_abx <- 0
  case_data_full$no_abx <- 0
  
  who_approved_names <- c("ARITHROMYCIN", 
                          "ASITHROMYCIN",
                          "ASITHROMYCIN, ERYTHROMYCIN",
                          "AZIMOX & AZITHRRO",
                          "AZITHROMYCIN",
                          "AZITHROMYCIN, ERYTHROMYCIN",
                          "AZITHROMYEIN",
                          "ASITHROMYCINN, NALIDIXICACID",
                          "AZITHTOMYCIN",
                          "CEFTRIAXONE",
                          "CEFTRIAXONE,AMOXICILLINE",
                          "CEFTRAXONE",
                          "CEFTRIAZONE",
                          "CIFROFLEXACIN",
                          "CIPROFLACIN",
                          "CIPROFLAXACIN",
                          "CIPROFLAXZCIN",
                          "CIPROFLEXACIN",
                          "CIPROFLOXACIN",
                          "CIPROFLOXAEIN",
                          "CIPROFLOXCIN",
                          "CIPROFLOXICIN",
                          "CIPROFLOZACIN",
                          "CIPROFLXACIN",
                          "CIPROFOLOXACIN",
                          "CIPROXIN",
                          "PIVMECILLINAM",  # note not seeing pivmecillinam in data
                          "TAB CIPROFLOXACIN",
                          "INJ CEFTRIAXONE.",
                          "INJ. CEFTRIAXONE (250MG)",
                          "INJ. CEFTRIAXONE 250MG",
                          "INJ. CEFTRIAXONE",
                          "CRFTRIAZONE",
                          "INJ CEFXONE") 
  
  azithro_names <- c("ARITHROMYCIN", 
                     "ASITHROMYCIN",
                     "ASITHROMYCIN, ERYTHROMYCIN",
                     "AZIMOX & AZITHRRO",
                     "AZITHROMYCIN",
                     "AZITHROMYCIN, ERYTHROMYCIN",
                     "AZITHROMYEIN",
                     "ASITHROMYCINN, NALIDIXICACID",
                     "AZITHTOMYCIN")
  
  maybe_effective_abx <- c("AMIKACIN", 
                           "AMOXICILLIN", 
                           "AMPICILLIN", 
                           "AUGMENTIN/CO-AMOXICLAV", 
                           "CEFACOR", 
                           "CEFIXIME", 
                           "CEFPODOXIME", 
                           "CEFTRAZIDIME", 
                           "CEFTRIAXONE", 
                           "CEFTRUOXIME", 
                           "CHLORAMPHENICOL", 
                           "CLARITHROMYCIN", 
                           "COTRIMOXAZOLE", 
                           "DOXYCYCLINE", 
                           "ERYTHROMYCIN", 
                           "FURAZOLIDONE", 
                           "GENTAMICIN", 
                           "LEVOFLOXACIN", 
                           "MEROPENEM", 
                           "NALIDIXIC ACID", 
                           "NIFUROXAZIDE", 
                           "NITROFURANTOIN", 
                           "OTHER MACROLIDES", 
                           "STREPTOMYCIN", 
                           "TETRACYCLINE", 
                           "THIAMPHENICOL",
                           "TETRACYCLINE CAP",
                           "TERACYCLINE",
                           "TETRACYCLINE CAPS" ,                                       
                           "TETRACYCLINE CAPSUL" ,
                           "TETRACYCLINE CAPSULE",
                           "AMOXACILLAN SYRUP",
                           "AMOXLYL" ,
                           "AMOXYL" ,
                           "AMOXYL SYRUP" ,
                           "CEFIXINE,COTRIMOXAZOL",
                           "COTRIMOXOZOL",
                           "AMOXICILINE",
                           "CO TRIMOXAZOLE",
                           "CCOTRIMOXAZOLE",
                           "CONTRIMOXAZOLE",
                           "CORTIMOXAZOLE",
                           "COTRINOKAZOLE",
                           "COTRIMOXAZOLE, METRONIDAZOLE", #check bc one effective, one ineffective
                           "AMPICILINE",
                           "BIODROXIL/ COTRIMOXAZOLE",
                           "COTRIMOXAZ0LE",
                           "COTRI",
                           "AMOXICYLLINE, METRONIDAZOLE",
                           "AMOXACILLIN",
                           "AMOXICILLINE SP",
                           "ADVENT",
                           "AMOXIL SYRUP" ,
                           "AMOXIL TAB.",
                           "BACITRACIN (COTRINOXAZOLE)",
                           "AMOSICILLIN",
                           "AMOXACILLY",                                             
                           "AMOXICILLIN,METRONIDAZOLE",
                           "AMOXICILLINE",
                           "AMOXOCILLAN",                                              
                           "AMOXOCILLEN",
                           "AMOXOCILLIN",                                              
                           "AMOXYCILLAN",
                           "AMOXYCILLIN",                                              
                           "AMOXYCILLIN & CLAVULANIC ACID",
                           "AMOXICILLINE, COTRIMOXAZOLE",
                           "BIODROXYL,AMOXICILLINE",
                           "CAPSULES AMOX",
                           "PARACETAMOL, AMOXICILLIN",
                           "PRINCIMOX",
                           "CHLORANPHENICOL",
                           "CHLORAMPHENICOLE",
                           "CLORAMPHENICOLE",
                           "CLORAFENICOL",
                           "CLOTRIMOXAZOLE",
                           "CO-TREMOXAZOLE" ,
                           "CO-TRIMOXAZOLE",
                           "COTRIMOXAZOL",
                           "CO-TRIMOXAZOLE, METRONIDAZOLE",
                           "CORTRIMOXAZOLE",
                           "SEPTRIN",                     # brand name?
                           "SEPTRIN (COTRIMOXAZOLE" ,
                           "SEPTRIN SYRP.",
                           "SEPTRIN SYRUP",
                           "SEPTRIN TAB." ,
                           "SEPTRIN TABS" ,
                           "SEPTRINE SYRUP",
                           "SEPTRINE SYRUP.",
                           "SYTUP SEPTRINE",  
                           "TRIMETHOPRIM & SULFAHENAZOLE",
                           "TRIMETHOPRIME SULFAMETHAZOLE" ,
                           "UCLAPRIM",
                           "SYP AMOXYCILLIN",
                           "SYP. COTRIMOXAZOLE",
                           "SYP AMOXYCILLIN",
                           "ERTHROMYCIN",
                           "ERTHROMYCIN",
                           "EROIHROMYCIN",
                           "ERYTHROMYCINE",
                           "ERYTHRONYCIN",
                           "ERCYTHCOMYSH",
                           "ERYTHNOMYCIN",                                             
                           "ERYTHOMYCINE" ,
                           "ERYTHROMYCIN, METRONIDAZOLE",
                           "ERYTHROMYCM",
                           "ERYTHROMYUME",
                           "ERYTROMICINE",
                           "ERYTROMIXYNE",
                           "SYP. ERYTHROMYCINE",
                           "CORTRIMOXAZOLE, ERYTHROMYCIN" ,
                           "COTRINOXAZOLE",
                           "COTRTMOXAZOLE",
                           "CORIMOXAZOLE",
                           "COTAMOXAZOLE" ,
                           "COTAMOXAZOLE",                                             
                           "COTRIMASCAZOLE",
                           "COTRIMOSCAZOLE",
                           "COTRIMOTAZOLE" ,
                           "COTRIMOXAZOLE & METRONIDAZOLE",
                           "COTRIMOXAZOLE SYRUP",
                           "COTRIMOXOZOLE",   
                           "COFRIMUXAZOLE",
                           "SEPRINE COTRIMOXASOLE",
                           "COMMOXAZOLE" ,
                           "COTRIMZAZOLE SYRUP",
                           "NALIDEXIC ACID",
                           "NADIDIXIC ACID",
                           "NALEDIXIC ACID",
                           "NALIDENIC ACID",
                           "NALIDENIE ACID",
                           "NALIDERIC ACID",
                           "NALIDEXIE ACID",
                           "NALIDIXIE ACID",
                           "NALIDXIC ACID",
                           "NALIDENE ACID",
                           "NALIDOXIC ACID",
                           "SYP. ENTAMEZOLE (SYP. METRONIDAZOL- OLE & NALIDIXIC ACID)",
                           "SYP. NEGRAM" ,
                           "METRONIDAZOLE, NALIIDIXIC ACID" ,
                           "SYP NALIDIX ACID",
                           "CEFIXINE",
                           "CEFIXIME SYRUP",
                           "CEFEXINE FRIHYDRATE",
                           "CEFIXIME SODIUM.",
                           "SYP. SEFEXIME & I/M INJ.",
                           "INJ CEFIXIME",
                           "SYP CEFIXME",
                           "MAXIMA (CEFIXIME)",
                           "TAXIM-O",
                           "DOXY",
                           "DOXY CYCLINE",
                           "HICANCIL",
                           "INJ-GENTAMICIN",
                           "O FLOXACIN",
                           "OFLAXACIN",
                           "XINTOF - M (OFLOXACIN & METRONIDAZOLE)" ,
                           "OFLOX & ORNIDAZOLE",
                           "OFLOXAC",
                           "OFLOXACAN",
                           "XINTOF-M",
                           "OFLOXACIN",
                           "OFLOXACIN & METRONIDAZOLE",
                           "OFLOXACIN & ORINIDAZOLE",
                           "OFLOXACIN & ORNIDAZNE",
                           "OFLOXACIN & ORNIDAZOLE",
                           "OFLOXACIN AMIKACIN",
                           "OFLOXACIN AND METRONIDAZOLE",
                           "OFLOXACIN AND ORNIDAZOLE",
                           "OFLOXACIN SUSPENSION",
                           "OFLOXACIN-HETROGYL(METRONIDA)",
                           "OFLOXACIN, ORNIDAZOLE",
                           "OFLOXACINE & ORNIDAZOLE",
                           "OFLOXACINS & METRONIDAZOLE",
                           "OFOXACIN & ORNIDAZOLE",
                           "SYR OFLONAC (OFLOXACIN)" ,
                           "SUSP. OFLOXACIN",
                           "SYR OFLONAC (OFLOXACIN)"  ,
                           "02 SUSPENSION (OFLOXACIN ORNIDAZOLE)",
                           "02 (OOFLOXACIN)",
                           "SUSP. ZENFLOX"  ,
                           "SYR TIDFLOX",
                           "ZENFLOX",
                           "ORNIDAZOLE SUSP. & OFLOXACIN SUSP.",
                           "SEPTIN",
                           "SEPTRAN",
                           "SEPTRIN FLAGYL",
                           "SEPTRIN- AMOXIL",
                           "SEPTRINE",
                           "SYP SEPTRAN" ,
                           "SEPTRIN AND FLAGYL",
                           "SEDTRIN",
                           "SEPTRUM,FLAGYL SYRUP",
                           "SEPRIN SYRUP",
                           "SULFAMETHOXAZOLE",
                           "BACTRIM 200/40MG",
                           "BACTRIM SP",
                           "BACTRIM 200/40MG",
                           "BACTRIM SP",
                           "WALLAMYCIN",
                           "SYR. WALAMYCIN",
                           "NOFLAIACIN",
                           "NORFLOXACIN",
                           "NORFLOXACIN & METRONIDAZOLE",
                           "NORFLOXACIN, METRONIDAZOLE",
                           "NORFLOXACIN, TINDAZOLE",
                           "NORFLOXACIN,METRONIDAZOLE",
                           "NAROFLOXACIN"  ,
                           "NORMET",
                           "NORMET SYR.",
                           "SYR. NORMET.",
                           "SYRUP NOMMET",
                           "DIARYL" ,
                           "NI FUROXAZIDE",
                           "AMPILOX",
                           "AMPICILLINE" ,
                           "AMOXILLINE",
                           "CEFOTAXIME",
                           "CEFOTAMIME",
                           "CEFOTAXIME",
                           "CEFOTAXIME(250MG)",
                           "INJ-CEFTAXIME",
                           "INJ. CEFOTAXIME",
                           "INJ. CEFOTAXIME 500MG",
                           "INJ/CEFOTAXIME" ,
                           "CETRIMAXAZOLE SYRUP",
                           "CURAM",
                           "COLISTIN",
                           "FOSOMIN(FOSFOMYCIN)",
                           "INJ CEFTROXIME" ,
                           "MONOCEF",
                           "NOR-METROGL.",
                           "OSPAMOX",
                           "SUL FAGMANIDINE",
                           "SULEAMIDE",
                           "SULFAGUANIDINE",
                           "SULFAGUANIDINE CP.",
                           "SULFAMIDE",
                           "CEFPODOXINE PROXELIL",
                           "CEFPODOXIME PROCEFIL",
                           "CAFOTAXIME")
  
  # this contains abx we know are ineffective. 
  # the ineffective variable will also contain people with write-ins that appear to be typos, "unknown", and no abx
  ineffective_abx <- c("CEFADROXIL", 
                       "CLINDAMYCIN", 
                       "CLOXACILLIN", 
                       "ENTAMIZOLE", 
                       "FLUCLOXACILLIN", 
                       "METRONIDAZOLE (FLAGYL)", 
                       "MUPIROCIN", 
                       "PENICILLIN", 
                       "PYRAZINAMIDE",
                       "FLAGYL",
                       "FLAGYL SYRUP",
                       "SRP FLAGYL" ,
                       "SYP. FLAGYL",
                       "FLAGYL, AMOXYL",
                       "FLAGYL, METRONDAZOLE",
                       "METRONIDAZOLE",
                       "METRONIDAZOLE SYRUP",
                       "METRONIDAZOL",
                       "METRCORIDAZOLE",
                       "METROINDAZOLE",
                       "METROMDAZOLE",
                       "METROMIDAZOLE",
                       "METROMIDAZOLE SYRU[",
                       "METRONDAZOLE SYP.",
                       "METRONEDAZOLE",
                       "METRONIDAJOLE",
                       "METRONIDAZOLE & COTIMOXAZOL",
                       "METRONIDAZOLE & DILOXANIDE FUROATE",
                       "METRONIDAZOLE AND FURAZOLIDONE",
                       "METRONIDAZOLE SUSPENSION",
                       "METRONIDAZOLE SYR",
                       "NONFLOXACIN & METROGYL SYRUP.",
                       "SYR. METROGYL",
                       "DELOXACIN & METRONIDAZOLE",
                       "CEFADROLIL",
                       "FAGYL",
                       "FLAGGYL TABLETS",                                          
                       "FLAGLY",
                       "FLAGYL SYP",
                       "FLAGYL TABLET",
                       "FLAGYL TABLETS.",                                          
                       "FLAGYL TABS",
                       "SYP FLAGYL",
                       "SYP FLEGYL",
                       "FLAYGILE" ,
                       "FLAGYL, CEPHRADINE, EEFAZOLIN",
                       "FLAGYL, METRONIDAZOLE"   ,
                       "SYP METROINIDAZOLE",
                       "SYP METRONIDAZOLE" ,
                       "SYP METRONIDAZOL",
                       "SYP. METRONIDAZOLE",
                       "SYRUP METRONIDAZOLE",
                       "MEHROGYL",
                       "METHONIDAZOLE" ,
                       "METOQUIDAZOLE",
                       "METRCOAIDAZOLE",                                           
                       "METRCORIDARZODE" ,
                       "METROGYL",                                                
                       "METROGYL, OFLOXACIN" ,
                       "METRONIDAZOLE, CEFIXIME",
                       "METRONIDOZAL",
                       "METRONIDOZOL",
                       "METRONIDOZOLE",
                       "METRONISAZOLE",
                       "METROZINE",
                       "METRUNIDAZEL",
                       "NORFLEX-METRONIDAZOLE",
                       "BIOBROXIL",
                       "BIODROXIL",
                       "BOIDROXIL",
                       "ENTAMEZOLE",
                       "GRAMAGYL",
                       "PARACETAMOL",
                       "SMECTA" ,
                       "FLUCAZOL",
                       "ENTEROQUINOL",
                       "CEFALEXINE",
                       "CEFAZOLIN SODIUM/IV" ,
                       "CEFLEX",
                       "CEPHALEXIN",
                       "SYP KEFLEX CEPHALEXIN",
                       "SPORNIDEX",
                       "CEPHESEIN",
                       "DELOXACINE, ORNIDAZOLE,",
                       "NITAZOXAMIDE",
                       "NITAZOXANIDE",
                       "SYNOBE",
                       "SYP(ENTAMIZOLE)",
                       "CEPHRADINE",
                       "(NOT ABLE TO)RED AND BLUE CAPSELS.",
                       "UNKNOWN",
                       "UNKNOW CAPSUL",
                       "UNKNOWN CAPSULES",
                       "UNKNOWN CAPSUL",
                       "SURUP",
                       "CATMAX",
                       "NETRO",
                       "GANDIDA",
                       "MOTHER NOT ABLE TO",
                       "NOT KNOWN",
                       "WHITE CAPSULE",
                       "CAPSULES UNKNOWN",
                       "UNKNOWN CAPSULE",
                       "OREX",
                       "GANIDA",
                       "BRISTOPEN",
                       "SEPHINE",
                       "TARICIN OZ",
                       "SYR. DISFLO-OZ",    
                       "SYR. 02",            
                       "DIDF",
                       "GENERAL",
                       "FUROYONE",
                       "ALQUNOL",
                       "02",
                       "OF-M",
                       "PEPHRODIN",
                       "ERYTHROGIA",
                       "FAYTHRONYCIN",
                       "NAME NOT KNOWN",
                       "DK",                
                       "1/2 INJ DK",
                       "INJ DK.",
                       "INJ. DK",
                       "IV INJ. NAME DK",
                       "INJ. EXCEF.",
                       "I.M. A/BIOTICS",    
                       "UNKNOWN.",
                       "DONT KNOW")
  
  # ----- Antibiotics pre-enroll -----
  
  # No longer including pre-enrollment antibiotics 6/11/26
  
  # Any abx = hometrt abx marked 1
  # case_data_full$any_abx <- ifelse(case_data_full$F4A_HOMETRT_AB == 1, 1, case_data_full$any_abx)
  # 
  # # WHO recommended abx
  # case_data_full$who_abx <- ifelse(case_data_full$F4A_HOMETRT_AB_SPEC.x %in% who_approved_names, 1, case_data_full$who_abx)
  # 
  # # azithro
  # case_data_full$azithro <- ifelse(case_data_full$F4A_HOMETRT_AB_SPEC.x %in% azithro_names, 1, case_data_full$azithro)
  # 
  # # Maybe effective antibiotics
  # case_data_full$maybe_eff_abx <- ifelse(case_data_full$F4A_HOMETRT_AB_SPEC.x %in% maybe_effective_abx, 1, case_data_full$maybe_eff_abx)
  # 
  # # Ineffective abx
  # case_data_full$ineff_abx <- ifelse(case_data_full$F4A_HOMETRT_AB_SPEC.x %in% ineffective_abx, 1, case_data_full$ineff_abx)
  
  # ----- Antibiotics hospital or home -----
  
  # cotrimoxazole - TRT_GIVE_CXL | TRT_PRES_CXL
  # gentamycin - TRT_GIVE_GENT | TRT_PRES_GENT
  # chloramphenicol - TRT_GIVE_CHLOR | TRT_PRES_CHLOR
  # erythromycin - TRT_GIVE_ERY | TRT_PRES_ERY
  # azithromycin - TRT_GIVE_AZI | TRT_PRES_AZI
  # other macrolides - TRT_GIVE_MACR | TRT_PRES_MACR
  # penicillin - TRT_GIVE_PEN | TRT_PRES_PEN
  # amoxycillin - TRT_GIVE_AMOX | TRT_PRES_AMOX
  # ampicillin - TRT_GIVE_AMPI | TRT_PRES_AMPI
  # nalidixic acid - TRT_PRES_NALID | TRT_GIVE_NALID
  # ciprofloxacin - TRT_PRES_CPNR | TRT_GIVE_CPNR
  # selexid/pivmecillinam - TRT_PRES_SLPY | TRT_GIVE_SLPY
  # other antibiotic - TRT_PRES_OTHR | TRT_GIVE_OTHR
  
  case_data_full$any_abx <- ifelse(case_data_full$F4B_TRT_GIVE_CXL == 1 | case_data_full$F4B_TRT_PRES_CXL.x == 1 |
                                     case_data_full$F4B_TRT_GIVE_GENT == 1 | case_data_full$F4B_TRT_PRES_GENT.x == 1 |
                                     case_data_full$F4B_TRT_GIVE_CHLOR == 1 | case_data_full$F4B_TRT_PRES_CHLOR.x == 1 |
                                     case_data_full$F4B_TRT_GIVE_ERY == 1 | case_data_full$F4B_TRT_PRES_ERY.x == 1 |
                                     case_data_full$F4B_TRT_GIVE_AZI == 1 | case_data_full$F4B_TRT_PRES_AZI.x == 1 |
                                     case_data_full$F4B_TRT_GIVE_MACR == 1 | case_data_full$F4B_TRT_PRES_MACR.x == 1 |
                                     case_data_full$F4B_TRT_GIVE_PEN == 1 | case_data_full$F4B_TRT_PRES_PEN.x == 1 |
                                     case_data_full$F4B_TRT_GIVE_AMOX == 1 | case_data_full$F4B_TRT_PRES_AMOX.x == 1 |
                                     case_data_full$F4B_TRT_GIVE_AMPI == 1 | case_data_full$F4B_TRT_PRES_AMPI.x == 1 |
                                     case_data_full$F4B_TRT_GIVE_NALID == 1 | case_data_full$F4B_TRT_PRES_NALID.x == 1 |
                                     case_data_full$F4B_TRT_GIVE_CPNR == 1 | case_data_full$F4B_TRT_PRES_CPNR.x == 1 |
                                     case_data_full$F4B_TRT_GIVE_SLPY == 1 | case_data_full$F4B_TRT_PRES_SLPY.x == 1 |
                                     case_data_full$F4B_TRT_GIVE_OTHR == 1 | case_data_full$F4B_TRT_PRES_OTHR.x == 1,
                                   1, case_data_full$any_abx)
  
  # note ceft does not have its own category
  case_data_full$who_abx <- ifelse(case_data_full$F4B_TRT_GIVE_AZI == 1 | case_data_full$F4B_TRT_PRES_AZI.x == 1 |
                                     case_data_full$F4B_TRT_GIVE_CPNR == 1 | case_data_full$F4B_TRT_PRES_CPNR.x == 1 |
                                     case_data_full$F4B_TRT_GIVE_SLPY == 1 | case_data_full$F4B_TRT_PRES_SLPY.x == 1,
                                   1, case_data_full$who_abx)
  
  case_data_full$azithro <- ifelse(case_data_full$F4B_TRT_GIVE_AZI == 1 | case_data_full$F4B_TRT_PRES_AZI.x == 1,
                                   1, case_data_full$azithro)
  # Maybe effective: 
  case_data_full$maybe_eff_abx <- ifelse(case_data_full$F4B_TRT_GIVE_CXL == 1 | case_data_full$F4B_TRT_PRES_CXL.x == 1 |
                                           case_data_full$F4B_TRT_GIVE_GENT == 1 | case_data_full$F4B_TRT_PRES_GENT.x == 1 |
                                           case_data_full$F4B_TRT_GIVE_CHLOR == 1 | case_data_full$F4B_TRT_PRES_CHLOR.x == 1 |
                                           case_data_full$F4B_TRT_GIVE_ERY == 1 | case_data_full$F4B_TRT_PRES_ERY.x == 1 |
                                           case_data_full$F4B_TRT_GIVE_MACR == 1 | case_data_full$F4B_TRT_PRES_MACR.x == 1 |
                                           case_data_full$F4B_TRT_GIVE_AMOX == 1 | case_data_full$F4B_TRT_PRES_AMOX.x == 1 |
                                           case_data_full$F4B_TRT_GIVE_AMPI == 1 | case_data_full$F4B_TRT_PRES_AMPI.x == 1 |
                                           case_data_full$F4B_TRT_GIVE_NALID == 1 | case_data_full$F4B_TRT_PRES_NALID.x == 1 |
                                           case_data_full$F4B_TRT_GIVE_OTHR == 1 | case_data_full$F4B_TRT_PRES_OTHR.x == 1,
                                         1, case_data_full$maybe_eff_abx)
  
  # Ineffective:
  case_data_full$ineff_abx <- ifelse(case_data_full$F4B_TRT_GIVE_PEN == 1 | case_data_full$F4B_TRT_PRES_PEN.x == 1, 1, case_data_full$ineff_abx)
  
  # make sure only fall into one category, taking highest
  case_data_full$ineff_abx <- ifelse(case_data_full$ineff_abx == 1 & (case_data_full$maybe_eff_abx == 1 | case_data_full$who_abx == 1), 0, case_data_full$ineff_abx)
  case_data_full$maybe_eff_abx <- ifelse(case_data_full$who_abx == 1 & case_data_full$maybe_eff_abx == 1, 0, case_data_full$maybe_eff_abx)
  
  # if ineff = 0, maybe = 0, and who = 0, categorize as no abx
  case_data_full$no_abx <- ifelse(case_data_full$who_abx == 0 & case_data_full$maybe_eff_abx == 0 & case_data_full$ineff_abx == 0, 1, case_data_full$no_abx)
  
  # Make 'all abx' variable so levels all in one
  case_data_full$all_abx <- ifelse(case_data_full$no_abx == 1 | case_data_full$ineff_abx == 1, 0, 
                                   ifelse(case_data_full$maybe_eff_abx == 1, 1, 2))
  case_data_full$all_abx <- factor(case_data_full$all_abx, levels = 0:2, labels = c("No or ineffective antibiotics",
                                                                                    "Possibly effective antibiotics",
                                                                                    "Guideline recommended antibiotics"))
  
  ### NEW 9/5/25 -- adding AMR data
  
  amr <- amr %>%
    filter(Pathotype == "Shigella") %>%
    mutate(SID = as.numeric(SID)) %>%
    select("SID", "mphA") %>%
    mutate(mphA = if_else(mphA == 0, "S", "I/R")) %>%
    rename('S_IR_azithro' = mphA)
  
  case_data_full <- left_join(case_data_full, amr, by = "SID") 
  
  case_data_full$azithro_given_S_IR <- ifelse(case_data_full$azithro == 1, case_data_full$S_IR_azithro, NA)
  case_data_full$ast_given_abx <- case_data_full$azithro_given_S_IR
  
  # create overall resistant and susceptible variables w same logic as EFGH
  # but revisit in meeting friday bc seems weird?? 
  case_data_full$resistant_WHO_approve <- apply(
    #case_data[, c("azithro_S_IR", "cipro_S_IR", "ceft_S_IR")],
    case_data_full[, c("S_IR_azithro")],
    1,
    function(x) if (any(x == "I/R", na.rm = TRUE)) 1 else 0
  )
  
  case_data_full$susceptible_WHO_approve <- apply(
    #case_data[, c("azithro_S_IR", "cipro_S_IR", "ceft_S_IR")],
    case_data_full[, c("S_IR_azithro")],
    1,
    function(x) if (any(x == "S", na.rm = TRUE)) 1 else 0
  )
  
  
  ## COVARIATES & SEVERITY VARS
  
  # create variables for other pathogen attribution
  case_data_full$adenovirus_attributable <- ifelse(case_data_full$adenovirus_40_41_afe > 0.5, 1, 0)
  case_data_full$aeromonas_attributable <- ifelse(case_data_full$aeromonas_afe > 0.5, 1, 0)
  case_data_full$astro_attributable <- ifelse(case_data_full$astrovirus_afe > 0.5, 1, 0)
  case_data_full$cryptosporidium_attributable <- ifelse(case_data_full$cryptosporidium_afe > 0.5, 1, 0)
  case_data_full$cyclospora_attributable <- ifelse(case_data_full$cyclospora_afe > 0.5, 1, 0)
  case_data_full$e_histolytica_attributable <- ifelse(case_data_full$e_histolytica_afe > 0.5, 1, 0)
  case_data_full$isospora_attributable <- ifelse(case_data_full$isospora_afe > 0.5, 1, 0)
  case_data_full$noro_attributable <- ifelse(case_data_full$norovirus_gii_afe > 0.5, 1, 0)
  case_data_full$rotavirus_attributable <- ifelse(case_data_full$rotavirus_afe > 0.5, 1, 0)
  case_data_full$salmonella_attributable <- ifelse(case_data_full$salmonella_afe > 0.5, 1, 0)
  case_data_full$sapovirus_attributable <- ifelse(case_data_full$sapovirus_afe > 0.5, 1, 0)
  case_data_full$st_etec_attributable <- ifelse(case_data_full$ST_ETEC_afe > 0.5, 1, 0)
  case_data_full$etec_attributable <- ifelse(case_data_full$LT_ETEC_afe > 0.5 | case_data_full$ST_ETEC_afe > 0.5, 1, 0)
  case_data_full$tepec_attributable <- ifelse(case_data_full$TEPEC_afe > 0.5, 1, 0)
  case_data_full$v_cholerae_attributable <- ifelse(case_data_full$v_cholerae_afe > 0.5, 1, 0)
  case_data_full$c_jejuni_coli_attributable <- ifelse(case_data_full$c_jejuni_coli_afe > 0.5, 1, 0)
  
  case_data_full$EAEC_attributable <- ifelse(case_data_full$EAEC_afe > 0.5, 1, 0)
  
  # Indicator for diarrhea with no known etiology
  case_data_full$no_etiology <- ifelse(rowSums(case_data_full[,c("shigella_attributable",
                                                                 "adenovirus_attributable",
                                                                 "aeromonas_attributable",
                                                                 "astro_attributable",
                                                                 "cryptosporidium_attributable",
                                                                 "cyclospora_attributable",
                                                                 "e_histolytica_attributable",
                                                                 "isospora_attributable",
                                                                 "noro_attributable",
                                                                 "rotavirus_attributable",
                                                                 "salmonella_attributable",
                                                                 "sapovirus_attributable",
                                                                 "st_etec_attributable",
                                                                 "tepec_attributable",
                                                                 "v_cholerae_attributable",
                                                                 "EAEC_attributable",
                                                                 "c_jejuni_coli_attributable")], na.rm = TRUE) == 0, 1, 0)
  
  
  # Note for 'detected' also using AFE because don't all seem to be on a scale to 35
  # note this is a different list than above 
  case_data_full$rota_detected <- ifelse(case_data_full$rotavirus_afe != 0, 1, 0)
  case_data_full$adeno_detected <- ifelse(case_data_full$adenovirus_40_41_afe != 0, 1, 0)
  case_data_full$etec_detected <- ifelse(case_data_full$ST_ETEC_afe != 0 | case_data_full$LT_ETEC_afe != 0 , 1, 0)
  case_data_full$crypto_detected <- ifelse(case_data_full$cryptosporidium_afe != 0 , 1, 0)
  case_data_full$astro_detected <- ifelse(case_data_full$astrovirus_afe != 0, 1, 0)
  case_data_full$noro_detected <- ifelse(case_data_full$norovirus_gii_afe != 0, 1, 0)
  case_data_full$tepec_detected <- ifelse(case_data_full$TEPEC_afe != 0, 1, 0)
  case_data_full$campy_detected <- ifelse(case_data_full$campylobacter_pan < 35, 1, 0) # NOTE not using AFE here
  case_data_full$sapo_detected <- ifelse(case_data_full$sapovirus_afe != 0, 1, 0)
  case_data_full$e_bieneusi_detected <- ifelse(case_data_full$e_bieneusi < 35, 1, 0)
  case_data_full$giardia_detected <- ifelse(case_data_full$giardia < 35 , 1, 0)
  case_data_full$EAEC_detected <- ifelse(case_data_full$EAEC_afe != 0, 1, 0)
  
  # 'New' variables = scaled TAC
  case_data_full$shigella_new <- case_data_full$shigella_eiec
  case_data_full$rotavirus_new <- case_data_full$rotavirus
  case_data_full$adenovirus_new <- case_data_full$adenovirus_40_41
  # ETEC new is max of st and lt etec (min tac quantity = max of rescaled tac quantity)
  case_data_full$etec_new <- ifelse(case_data_full$ST_ETEC > case_data_full$LT_ETEC, case_data_full$ST_ETEC, case_data_full$LT_ETEC)
  case_data_full$st_etec_new <- case_data_full$ST_ETEC
  case_data_full$lt_etec_new <- case_data_full$LT_ETEC
  case_data_full$crypto_new <- case_data_full$cryptosporidium
  case_data_full$astro_new <- case_data_full$astrovirus
  case_data_full$noro_new <- case_data_full$norovirus_gii
  case_data_full$tepec_new <- case_data_full$TEPEC
  # campy no scaled
  case_data_full$campy_new <- ifelse(case_data_full$campylobacter_pan > 35, 35, case_data_full$campylobacter_pan)
  case_data_full$campy_new <- (35 - case_data_full$campy_new) / 3.322
  case_data_full$sapo_new <- case_data_full$sapovirus
  # e bieneusi, giardia, eaec no scaled
  case_data_full$e_bieneusi_new <- ifelse(case_data_full$e_bieneusi > 35, 35, case_data_full$e_bieneusi)
  case_data_full$e_bieneusi_new <- (35 - case_data_full$e_bieneusi_new) / 3.322
  case_data_full$giardia_new <- ifelse(case_data_full$giardia> 35, 35, case_data_full$giardia)
  case_data_full$giardia_new <- (35 - case_data_full$giardia_new) / 3.322
  case_data_full$EAEC_new <-  case_data_full$EAEC.y
  
  case_data_full$v_cholerae_new <- case_data_full$v_cholerae
  
  case_data_full$salmonella_new <- case_data_full$salmonella
  case_data_full$c_jejuni_coli_new <- case_data_full$c_jejuni_coli
  
  # # other path cols
  # other_path <- c("STEC_afe", "EAEC_afe", "TEPEC_afe", "LT_ETEC_afe",
  #                       "h_pylori_afe", "aeromonas_afe", "sapovirus_afe",
  #                       "astrovirus_afe", "salmonella_afe", "isospora_afe",
  #                       "c_jejuni_coli_afe", "norovirus_gii_afe", "v_cholerae_afe",
  #                       "cyclospora_afe", "e_histolytica_afe")
  # 
  # # Calculate the 'other_path_attributable' column
  # case_data_full$other_path_attributable <- ifelse(
  #   rowSums(case_data_full[, other_path] > 0.5, na.rm = TRUE) > 0, 1, 0
  # )
  
  # create variable for days between enrollment and 60-day follow-up
  case_data_full$followupdate <- as.Date(case_data_full$F5_DATE)
  case_data_full$followup_days <- case_data_full$followupdate - case_data_full$enrolldate
  case_data_full$I_followup_days <- ifelse(is.na(case_data_full$followup_days), 0, 1)
  case_data_full$I_followup_days_x_followup_days <- ifelse(is.na(case_data_full$followup_days), 0, case_data_full$followup_days)
  
  # get other variables of interest
  
  # covariates:
  # sex (F2_GENDER - 0 boy, 1 girl)
  # age (F2_age - age in months)
  # HAZ (F4B_MED_HAZ - height-for-age z-score at enrollment)
  # site  (SITE - site - 1-7 but not seeing dictionary)
  # primary_caretaker_edu (F4A_PRIM_SCHL - ordinal 1-5, 6 religious only, 7 NA)
  # child_under_60mo_hh (F4A_YNG_CHILDREN - how many children younger than 60 months live in household)
  # breastfed - (F4A_BREASTFED - is child currently breastfed - 0 no, 1 partial, 2 exclusive)
  
  # --- NEW 12/9/24 FOR SES TO MATCH VIDA: ---
  # safe_water
  # safe_sanit
  # electricity
  
  # Make joint monitoring programme definition of improved water
  # based on data dictionary for VIDA
  case_data_full <- case_data_full %>%
    mutate(
      MS_WATER = coalesce(F4A_MS_WATER, F7_MS_WATER),
      TIME_WATER = coalesce(F4A_TIME_WATER, F7_TIME_WATER),
      WATER_AVAIL = coalesce(F4A_WATER_AVAIL, F7_WATER_AVAIL),
      MS_SPEC = coalesce(F4A_MS_SPEC, F7_MS_SPEC) 
    ) %>%
    mutate(
      safe_water = case_when(
        (MS_WATER %in% c(1, 2, 17, 7, 8, 9, 10, 11, 15, 16, 3) | grepl("tap|pipe", MS_SPEC, ignore.case = TRUE)) &
          TIME_WATER %in% c(1, 2) & WATER_AVAIL %in% c(1, 2) ~ 1,
        (MS_WATER %in% c(1, 2, 9, 15) | grepl("pipe", MS_SPEC, ignore.case = TRUE)) &
          is.na(TIME_WATER) & WATER_AVAIL %in% c(1, 2) ~ 1,
        (MS_WATER %in% c(1, 2, 17, 7, 8, 9, 10, 11, 15, 16, 3) | grepl("tap|pipe", MS_SPEC, ignore.case = TRUE)) &
          TIME_WATER %in% c(1, 2) & WATER_AVAIL %in% c(3, 4) ~ 2,
        (MS_WATER %in% c(1, 2, 9, 15) | grepl("pipe", MS_SPEC, ignore.case = TRUE)) &
          is.na(TIME_WATER) & WATER_AVAIL %in% c(3, 4) ~ 2,
        (MS_WATER %in% c(1, 2, 17, 7, 8, 9, 10, 11, 15, 16, 3) | grepl("tap|pipe", MS_SPEC, ignore.case = TRUE)) &
          TIME_WATER %in% c(3, 4, 5) ~ 3,
        MS_WATER %in% c(4, 5, 12) | MS_SPEC == "open" ~ 4,
        MS_WATER %in% c(6, 13, 14, 18, 19) ~ 5,
        MS_WATER == 18 & !grepl("tap|pipe|open", MS_SPEC, ignore.case = TRUE) | is.na(TIME_WATER) | is.na(WATER_AVAIL) ~ NA_real_,
        TRUE ~ NA_real_ # Default case
      )
    )
  
  # GEMS -- > VIDA
  # GEMS X = VIDA Y
  # GEMS 1 = VIDA 1
  # 2 = 4
  # 3 = 5 or 6
  # 4 = 2 
  # 5 = 9
  # 6 = 10 (specify)
  # 7 = 4
  # HANGING, HANGING LATRINE, HANGING TOILET = 8
  
  # note this is closest to VIDA, don't have all the same options
  case_data_full <- case_data_full %>%
    mutate(
      MAIN_WASTE = coalesce(F4A_FAC_WASTE, F7_FAC_WASTE),
      MAIN_SPEC = coalesce(F4A_FAC_SPEC, F7_FAC_SPEC),
      SHARE_FAC = coalesce(F4A_SHARE_FAC, F7_SHARE_FAC),
      safe_sanit = case_when(
        MAIN_WASTE %in% c(1, 4, 7, 3) & SHARE_FAC == 0 ~ 1,
        MAIN_WASTE %in% c(1, 4, 7, 3) & SHARE_FAC >= 1 ~ 2,
        MAIN_WASTE %in% c(3) | SHARE_FAC %in% c("HANGING", "HANGING LATRINE", "HANGING TOILET") ~ 3,
        MAIN_WASTE %in% c(5, 6) ~ 4,
        is.na(MAIN_WASTE) | is.na(SHARE_FAC) ~ NA_real_,
        TRUE ~ NA_real_ 
      )
    )
  
  case_data_full$safe_water <- factor(case_data_full$safe_water, levels = 1:5, labels = c("Safely managed",
                                                                                          "Basic",
                                                                                          "Limited",
                                                                                          "Unimproved",
                                                                                          "Surface water"))
  
  case_data_full$safe_sanit <- factor(case_data_full$safe_sanit, levels = 1:4, labels = c("Safely managed and basic",
                                                                                          "Limited",
                                                                                          "Unimproved",
                                                                                          "Open defication"))
  
  # severity: 
  # - dysentery (F4B_OUTCOME_DYS - left hospital with dysentery)
  #             (F4A_DRH_STOOLS = 4 - describe diarrhea as bloody)
  #             (F4A_DRH_BLOOD = 1 - blood in stool since illness began, 9 NA)
  #
  # - vomitting more than 3x per day (F4A_DRH_VOMIT = 1 - vomit 3 or more times / day since illness began)
  # - fever (F4A_DRH_FEVER = 1 - fever 38C or higher since illness began )
  # - loose stools (categorical number of loose stools on worst day )
  # - dehydration: define same way as VIDA
  
  # based on vida definition
  case_data_full$SEVERE_DEHYD <- ifelse(case_data_full$F4B_MENTAL == 2, 1, 0)
  case_data_full$SEVERE_DEHYD <- ifelse(case_data_full$F4B_EYES == 1, case_data_full$SEVERE_DEHYD + 1, case_data_full$SEVERE_DEHYD)
  case_data_full$SEVERE_DEHYD <- ifelse(case_data_full$F4A_DRH_LESSDRINK == 1 | case_data_full$F4A_DRH_UNDRINK == 1, case_data_full$SEVERE_DEHYD + 1, case_data_full$SEVERE_DEHYD)
  case_data_full$SEVERE_DEHYD <- ifelse(case_data_full$F4B_SKIN == 2, case_data_full$SEVERE_DEHYD + 1, case_data_full$SEVERE_DEHYD)
  
  case_data_full$MOD_DEHYD <- ifelse(case_data_full$F4B_MENTAL == 1, 1, 0)
  case_data_full$MOD_DEHYD <- ifelse(case_data_full$F4B_EYES == 1, case_data_full$MOD_DEHYD + 1, case_data_full$MOD_DEHYD)
  case_data_full$MOD_DEHYD <- ifelse(case_data_full$F4A_DRH_THIRST == 1, case_data_full$MOD_DEHYD + 1, case_data_full$MOD_DEHYD)
  case_data_full$MOD_DEHYD <- ifelse(case_data_full$F4B_SKIN == 1, case_data_full$MOD_DEHYD + 1, case_data_full$MOD_DEHYD)
  
  case_data_full$who_dehyd <- ifelse(case_data_full$SEVERE_DEHYD >= 2, 2, 0)
  case_data_full$who_dehyd <- ifelse(case_data_full$F4B_EYES == 0 & 
                                       case_data_full$SEVERE_DEHYD == 0 & 
                                       case_data_full$MOD_DEHYD >= 2, 1, case_data_full$who_dehyd)
  case_data_full$who_dehyd <- ifelse(case_data_full$F4B_EYES == 0 & 
                                       case_data_full$SEVERE_DEHYD == 1 & 
                                       case_data_full$MOD_DEHYD >= 1, 1, case_data_full$who_dehyd)
  case_data_full$who_dehyd <- ifelse(case_data_full$F4B_EYES == 1 & 
                                       case_data_full$SEVERE_DEHYD == 1 & 
                                       case_data_full$MOD_DEHYD >= 2, 1, case_data_full$who_dehyd)
  
  case_data_full$who_dehyd <- factor(case_data_full$who_dehyd, levels = 0:2, labels = c("No dehydration",
                                                                                        "Some dehydration",
                                                                                        "Severe dehydration"))
  
  case_data_full$dysentery <- ifelse(case_data_full$F4B_OUTCOME_DYS == 1 | 
                                       case_data_full$F4A_DRH_STOOLS == 4 | 
                                       case_data_full$F4A_DRH_BLOOD == 1, 1, 0)
  
  case_data_full$lsstools <- factor(case_data_full$F4A_MAX_STOOLS, levels = 1:3, labels = c("6 or less per day",
                                                                                            "7 to 10 per day",
                                                                                            "Over 10 per day"))
  
  # Duration before enrollment
  case_data_full$duration_pre_enroll <- case_data_full$F4A_DRH_DAYS
  
  case_data_select <- case_data_full %>%
    select(pid,
           first_id,
           child_id,
           F2_GENDER,
           F2_AGE,
           F4B_MED_HAZ,
           F4B_WAZ,
           SITE,
           F4A_PRIM_SCHL,
           F4A_YNG_CHILDREN,
           F4A_BREASTFED,
           dysentery,
           lsstools,
           who_dehyd,
           F4A_DRH_VOMIT,
           F4A_DRH_FEVER,
           duration_pre_enroll,
           safe_sanit,
           safe_water, 
           ses_quintile,
           shigella_attributable,
           shigella_attributable_tac,
           shigella_attributable_culture,
           Shig_flex, 
           Shig_sonnei,
           rotavirus_attributable,
           etec_attributable,
           adenovirus_attributable,
           aeromonas_attributable,
           astro_attributable,
           cryptosporidium_attributable,
           cyclospora_attributable,
           e_histolytica_attributable,
           isospora_attributable,
           noro_attributable,
           rotavirus_attributable,
           salmonella_attributable,
           sapovirus_attributable,
           st_etec_attributable,
           tepec_attributable,
           v_cholerae_attributable,
           c_jejuni_coli_attributable,
           EAEC_attributable,
           no_etiology,
           rota_detected,
           adeno_detected,
           etec_detected,
           crypto_detected,
           astro_detected,
           noro_detected,
           tepec_detected,
           campy_detected,
           sapo_detected,
           e_bieneusi_detected,
           giardia_detected,
           EAEC_detected,
           shigella_eiec_afe,
           adenovirus_40_41_afe,
           aeromonas_afe,
           astrovirus_afe,
           cryptosporidium_afe,
           cyclospora_afe,
           e_histolytica_afe ,
           isospora_afe,
           norovirus_gii_afe,
           rotavirus_afe,
           salmonella_afe,
           sapovirus_afe,
           ST_ETEC_afe,
           LT_ETEC_afe,
           TEPEC_afe,
           shigella_new,
           rotavirus_new,
           adenovirus_new,
           st_etec_new,
           lt_etec_new,
           etec_new,
           crypto_new,
           astro_new,
           noro_new,
           tepec_new,
           campy_new,
           sapo_new,
           e_bieneusi_new,
           giardia_new,
           EAEC_new,
           v_cholerae_new,
           salmonella_new,
           c_jejuni_coli_new,
           v_cholerae_afe,
           followup_days,
           I_followup_days,
           I_followup_days_x_followup_days,
           hazdiff,
           hazd60,
           any_abx,
           who_abx,
           maybe_eff_abx,
           ineff_abx,
           no_abx,
           all_abx,
           ast_given_abx, 
           S_IR_azithro,
           resistant_WHO_approve,
           susceptible_WHO_approve
    ) %>%
    rename("sex" = F2_GENDER,
           "age" = F2_AGE,
           "enr_haz" = F4B_MED_HAZ,
           "enr_waz" = F4B_WAZ,
           "site" = SITE,
           "prim_caregiver_edu" = F4A_PRIM_SCHL,
           "num_hh_lt5" = F4A_YNG_CHILDREN,
           "breastfed" = F4A_BREASTFED,
           "vomit" = F4A_DRH_VOMIT,
           "fever" = F4A_DRH_FEVER,
           "shig_flex" = Shig_flex,
           "shig_sonnei" = Shig_sonnei) %>%
    mutate(fever = if_else(fever == 9, NA, fever)) # 9s for fever are NAs
  
  # if culture negative, shig_flex and shig_sonnei should be NA (we do not know types of shig tac positive who were cultuer neg)
  case_data_select <- case_data_select %>%
    mutate(shig_flex = if_else(shigella_attributable_culture == 0, NA, shig_flex),
           shig_sonnei = if_else(shigella_attributable_culture == 0, NA, shig_sonnei))
  
  # Transform factors as needed
  
  # Sites: (key from matching IDs in gems1_form03_clean)
  # 1 - The Gambia
  # 2 - Mali
  # 3 - Mozambique
  # 4 - Kenya
  # 5 - India
  # 6 - Bangladesh
  # 7 - Pakistan
  
  case_data_select$site <- factor(case_data_select$site, levels = 1:7, labels = c("The Gambia",
                                                                                  "Mali",
                                                                                  "Mozambique",
                                                                                  "Kenya",
                                                                                  "India",
                                                                                  "Bangladesh",
                                                                                  "Pakistan"))
  
  case_data_select$sex <- factor(case_data_select$sex, levels = c(0,1), labels = c("male", "female"))
  
  case_data_select$prim_caregiver_edu <- ifelse(case_data_select$prim_caregiver_edu == 7, NA, case_data_select$prim_caregiver_edu)
  case_data_select$prim_caregiver_edu <- factor(case_data_select$prim_caregiver_edu, levels = 1:6, labels = c("No formal schooling",
                                                                                                              "Less than primary",
                                                                                                              "Completed primary",
                                                                                                              "Completed secondary",
                                                                                                              "Post secondary",
                                                                                                              "Religious education only"))
  
  case_data_select$prim_caregiver_edu_bin <- ifelse(case_data_select$prim_caregiver_edu %in% c("No formal schooling",
                                                                                               "Less than primary",
                                                                                               "Completed primary",
                                                                                               "Religious education only"), 0, 1)
  
  
  case_data_select$breastfed <- factor(case_data_select$breastfed, levels = c(0,1,2), labels = c("no",
                                                                                                 "partial",
                                                                                                 "exclusive"))
  
  case_data_select$ses_quintile <- factor(case_data_select$ses_quintile, levels = 1:5, labels = c("1st quintile of SES",
                                                                                                  "2nd quintile of SES",
                                                                                                  "3rd quintile of SES",
                                                                                                  "4th quintile of SES",
                                                                                                  "5th quintile of SES"))
  
  # LABEL EVERYTHING FOR DESCRIPTIVE TABLES FUNCTION
  case_data_select <- case_data_select %>%
    set_variable_labels(pid = "Participant ID",
                        child_id = "Child ID",
                        first_id = "First ID",
                        sex = "Sex",
                        age = "Age (months)",
                        enr_haz = "HAZ at enrollment",
                        enr_waz = "WAZ at enrollment",
                        site = "Enrollment site",
                        prim_caregiver_edu = "Education level of primary caregiver secondary school or greater",
                        prim_caregiver_edu_bin = "Primary caregiver education level (binary)",
                        num_hh_lt5 = "Number of children in household <5 years old",
                        breastfed = "Is child currently breastfed?",
                        dysentery = "Dysentery",
                        who_dehyd = "WHO defined dehydration",
                        lsstools = "Maximum number of loose stools in 24hr period",
                        vomit = "Vomitted >3x/day since illness began",
                        fever = "Fever >=38C since illness began",
                        shigella_attributable = "Shigella attributable diarrhea via TAC (AFE>0.5) or culture",
                        shigella_attributable_tac = "Shigella attributable diarrhea via TAC (AFE>0.5)",
                        shigella_attributable_culture = "Shigella attributable diarrhea via culture",
                        shig_flex = "S. flexneri",
                        shig_sonnei = "S. sonnei",
                        rotavirus_attributable = "Rotavirus attributable diarrhea (AFE>0.5)",
                        etec_attributable = "ST/LT ETEC attributable diarrhea (AFE>0.5)",
                        cryptosporidium_attributable = "Cryptosporidium attributable diarrhea (AFE>0.5)",
                        adenovirus_attributable = "Adenovirus attributable diarrhea (AFE>0.5)",
                        astro_attributable = "Astrovirus attributable diarrhea (AFE>0.5)",
                        noro_attributable = "Norovirus attributable diarrhea (AFE>0.5)",
                        tepec_attributable = "TEPEC attributable diarrhea (AFE>0.5)", 
                        sapovirus_attributable = "Sapovirus attributable diarrhea (AFE>0.5)",
                        EAEC_attributable = "EAEC attributable diarrhea (AFE>0.5)",
                        no_etiology = "Diarrhea not attributable to any pathogens (by AFE)",
                        rota_detected = "Rotavirus detected",
                        adeno_detected = "Adenovirus detected",
                        etec_detected = "ETEC detected",
                        crypto_detected = "Cryptosporidium detected",
                        astro_detected = "Astrovirus detected",
                        noro_detected = "Norovirus detected",
                        tepec_detected = "tEPEC detected",
                        campy_detected = "Campylobacter detected",
                        sapo_detected = "Sapovirus detected",
                        e_bieneusi_detected = "E bieneusi detected",
                        giardia_detected = "Giardia detected",
                        EAEC_detected = "EAEC detected",
                        hazdiff = "Difference in HAZ enrollment - day 60",
                        hazd60 = "HAZ on day 60 post-enrollment",
                        followup_days = "Days between enrollment and follow-up",
                        any_abx = "Received any antibiotics",
                        who_abx = "Received WHO approved antibiotics",
                        maybe_eff_abx = "Received maybe effective antibiotics",
                        all_abx = "Type of antibiotics recieved",
                        ineff_abx = "Recieved ineffective antibiotics",
                        no_abx = "Did not receive antibiotics",
                        safe_water = "JMP Improved Water" ,
                        safe_sanit = "JMP Improved Sanitation",
                        ses_quintile = "SES quintile by site")
  
  return(case_data_select)
  
}

# Create and save analysis-ready datasets 

gems_data_full <- prep_gems() # Includes all cases
gems_data  <- gems_data_full[which(!is.na(gems_data_full$shigella_attributable_tac)), , drop = FALSE] # Removes cases missing Shigella TAC

saveRDS(gems_data, here::here("data/gems_data/gems_data.Rds"))
saveRDS(gems_data_full, here::here("data/gems_data/gems_data_full.Rds"))

# --------------------------------------------------------------------
# CASE-CONTROL DATA PREP 
# --------------------------------------------------------------------

#' Function to create case control dataset for gems
#' 
#' Potential case vs control definitions:
#' - All diarrhea cases vs controls (all_diar)
#' - TAC or culture Shigella attributable diarrhea vs no diarrhea (tac_or_culture_shig_diar)
#' - TAC shigella diarrhea vs no diarrhea (tac_shig_diar)
#' - Culture Shigella attributable diarrhea vs no diarrhea (culture_shig_diar)
prep_gems_case_control <- function(case_def = "tac_or_culture_shig_diar"){
  
  # Load data
  case_control_data <- read.csv(here::here("data/gems_data/raw_data/GEMS-1 case-control dataset.csv"))
  load(here::here("data/gems_data/raw_data/GEMS with AFEs.RData")) # called 'data'
  amr <- haven::read_sas(here::here("data/gems_data/raw_data/amr.sas7bdat")) 
  
  # -------------------------------------------------
  # Make PCA SES score in cases and controls by site 
  # -------------------------------------------------
  case_control_data$fridge <- ifelse(!is.na(case_control_data$F4A_HOUSE_FRIDGE), 
                                     case_control_data$F4A_HOUSE_FRIDGE,
                                     case_control_data$F7_HOUSE_FRIDGE)
  case_control_data$tv <- ifelse(!is.na(case_control_data$F4A_HOUSE_TELE),
                                 case_control_data$F4A_HOUSE_TELE,
                                 case_control_data$F7_HOUSE_TELE)
  case_control_data$electricity <- ifelse(!is.na(case_control_data$F4A_HOUSE_ELEC),
                                          case_control_data$F4A_HOUSE_ELEC,
                                          case_control_data$F7_HOUSE_ELEC)
  case_control_data$motorcycle <- ifelse(!is.na(case_control_data$F4A_HOUSE_SCOOT),
                                         case_control_data$F4A_HOUSE_SCOOT,
                                         case_control_data$F7_HOUSE_SCOOT)
  case_control_data$radio <- ifelse(!is.na(case_control_data$F4A_HOUSE_RADIO),
                                    case_control_data$F4A_HOUSE_RADIO,
                                    case_control_data$F7_HOUSE_RADIO)
  case_control_data$bike <- ifelse(!is.na(case_control_data$F4A_HOUSE_BIKE),
                                   case_control_data$F4A_HOUSE_BIKE,
                                   case_control_data$F7_HOUSE_BIKE)
  case_control_data$car <- ifelse(!is.na(case_control_data$F4A_HOUSE_CAR),
                                  case_control_data$F4A_HOUSE_CAR,
                                  case_control_data$F7_HOUSE_CAR)
  case_control_data$boat <- ifelse(!is.na(case_control_data$F4A_HOUSE_BOAT),
                                   case_control_data$F4A_HOUSE_BOAT,
                                   case_control_data$F7_HOUSE_BOAT)
  case_control_data$phone <- ifelse(!is.na(case_control_data$F4A_HOUSE_PHONE),
                                    case_control_data$F4A_HOUSE_PHONE,
                                    case_control_data$F7_HOUSE_PHONE)
  case_control_data$cart <- ifelse(!is.na(case_control_data$F4A_HOUSE_CART),
                                   case_control_data$F4A_HOUSE_CART,
                                   case_control_data$F7_HOUSE_CART)
  case_control_data$agland <- ifelse(!is.na(case_control_data$F4A_HOUSE_AGLAND),
                                     case_control_data$F4A_HOUSE_AGLAND,
                                     case_control_data$F7_HOUSE_AGLAND)
  
  cols_for_pca <- c("fridge", "tv", "electricity", "motorcycle",
                    "radio", "bike", "car", "boat", "phone", "cart", "agland")
  
  ses_pca <- prcomp(case_control_data[complete.cases(case_control_data[,cols_for_pca]), cols_for_pca], 
                    center = TRUE, scale = TRUE)
  
  pca_score <- as.matrix(case_control_data[, cols_for_pca]) %*% ses_pca$rotation[, 1]
  
  case_control_data$ses_score <- pca_score[, 1] 
  
  case_control_data <- case_control_data %>%
    group_by(SITE) %>%
    mutate(ses_quintile = ntile(desc(ses_score), 5)) # flip so lower = worse, higher = best
  
  case_control_data$ses_quintile <- factor(case_control_data$ses_quintile,  levels = 1:5, labels = c("1st quintile of SES",
                                                                                                     "2nd quintile of SES",
                                                                                                     "3rd quintile of SES",
                                                                                                     "4th quintile of SES",
                                                                                                     "5th quintile of SES"))
  
  # case_control_data CHILDID all end in 1 and do not match data, get rid of
  case_control_data$child_id <- as.numeric(substr(as.character(case_control_data$CHILDID), 
                                                  1, 
                                                  nchar(as.character(case_control_data$CHILDID)) - 1))
  
  case_control_data$case_id <- as.numeric(substr(as.character(case_control_data$CASEID), 
                                                 1, 
                                                 nchar(as.character(case_control_data$CASEID)) - 1))
  
  
  case_control_data$first_id <- case_control_data$FIRSTID
  
  # merge on child_id minus the 1
  data_full <- left_join(case_control_data, data, 
                         by = c("child_id" = "childid"), 
                         keep = FALSE)
  
  # Rename transformed quantities to end with '_new'
  data_full <- data_full %>%
    rename('shigella_new' = shigella_eiec,
           'rotavirus_new' = rotavirus,
           'adenovirus_new' = adenovirus_40_41,
           'st_etec_new' = ST_ETEC,
           'lt_etec_new' = LT_ETEC,
           'crypto_new' = cryptosporidium,
           'astro_new' = astrovirus,
           'noro_new' = norovirus_gii,
           'tepec_new' = TEPEC,
           'sapo_new' = sapovirus)
  
  data_full$etec_new <- ifelse(data_full$st_etec_new > data_full$lt_etec_new, data_full$st_etec_new, data_full$lt_etec_new)
  data_full$campy_new <- ifelse(data_full$campylobacter_pan > 35, 35, data_full$campylobacter_pan)
  data_full$campy_new <- (35 - data_full$campy_new) / 3.322
  data_full$e_bieneusi_new <- ifelse(data_full$e_bieneusi > 35, 35, data_full$e_bieneusi)
  data_full$e_bieneusi_new <- (35 - data_full$e_bieneusi_new) / 3.322
  data_full$giardia_new <- ifelse(data_full$giardia> 35, 35, data_full$giardia)
  data_full$giardia_new <- (35 - data_full$giardia_new) / 3.322
  data_full$EAEC_new <-  data_full$EAEC.y
  data_full$v_cholerae_new <- data_full$v_cholerae
  data_full$salmonella_new <- data_full$salmonella
  data_full$c_jejuni_coli_new <- data_full$c_jejuni_coli
  
  # create variables for other pathogen attribution & detected
  data_full$adenovirus_attributable <- ifelse(data_full$adenovirus_40_41_afe > 0.5, 1, 0)
  data_full$aeromonas_attributable <- ifelse(data_full$aeromonas_afe > 0.5, 1, 0)
  data_full$astro_attributable <- ifelse(data_full$astrovirus_afe > 0.5, 1, 0)
  data_full$cryptosporidium_attributable <- ifelse(data_full$cryptosporidium_afe > 0.5, 1, 0)
  data_full$cyclospora_attributable <- ifelse(data_full$cyclospora_afe > 0.5, 1, 0)
  data_full$e_histolytica_attributable <- ifelse(data_full$e_histolytica_afe > 0.5, 1, 0)
  data_full$isospora_attributable <- ifelse(data_full$isospora_afe > 0.5, 1, 0)
  data_full$noro_attributable <- ifelse(data_full$norovirus_gii_afe > 0.5, 1, 0)
  data_full$rotavirus_attributable <- ifelse(data_full$rotavirus_afe > 0.5, 1, 0)
  data_full$salmonella_attributable <- ifelse(data_full$salmonella_afe > 0.5, 1, 0)
  data_full$sapovirus_attributable <- ifelse(data_full$sapovirus_afe > 0.5, 1, 0)
  data_full$st_etec_attributable <- ifelse(data_full$ST_ETEC_afe > 0.5, 1, 0)
  data_full$etec_attributable <- ifelse(data_full$LT_ETEC_afe > 0.5 | data_full$ST_ETEC_afe > 0.5, 1, 0)
  data_full$tepec_attributable <- ifelse(data_full$TEPEC_afe > 0.5, 1, 0)
  data_full$v_cholerae_attributable <- ifelse(data_full$v_cholerae_afe > 0.5, 1, 0)
  data_full$c_jejuni_coli_attributable <- ifelse(data_full$c_jejuni_coli_afe > 0.5, 1, 0)
  data_full$EAEC_attributable <- ifelse(data_full$EAEC_afe > 0.5, 1, 0)
  
  # note this is a different list than above 
  data_full$rota_detected <- ifelse(data_full$rotavirus_afe != 0, 1, 0)
  data_full$adeno_detected <- ifelse(data_full$adenovirus_40_41_afe != 0, 1, 0)
  data_full$etec_detected <- ifelse(data_full$ST_ETEC_afe != 0 | data_full$LT_ETEC_afe != 0 , 1, 0)
  data_full$crypto_detected <- ifelse(data_full$cryptosporidium_afe != 0 , 1, 0)
  data_full$astro_detected <- ifelse(data_full$astrovirus_afe != 0, 1, 0)
  data_full$noro_detected <- ifelse(data_full$norovirus_gii_afe != 0, 1, 0)
  data_full$tepec_detected <- ifelse(data_full$TEPEC_afe != 0, 1, 0)
  data_full$campy_detected <- ifelse(data_full$campylobacter_pan < 35, 1, 0) # NOTE not using AFE here
  data_full$sapo_detected <- ifelse(data_full$sapovirus_afe != 0, 1, 0)
  data_full$e_bieneusi_detected <- ifelse(data_full$e_bieneusi < 35, 1, 0)
  data_full$giardia_detected <- ifelse(data_full$giardia < 35 , 1, 0)
  data_full$EAEC_detected <- ifelse(data_full$EAEC_afe != 0, 1, 0)
  
  # culture pathogen attribution
  data_full$rotavirus_elisa <- data_full$F18_RES_ROTAVIRUS
  data_full$crypto_elisa <- data_full$F18_RES_CRYPTOSPOR
  data_full$ecoli_culture <- ifelse(data_full$F16CORR_ECOLI == 1 & data_full$F17CORR_RESULT_ESTA == 1, 1, 0)
  data_full$adenovirus_elisa <- data_full$F18_RES_ADENO4041
  
  data_full$etec_pcr <- ifelse(data_full$ETEC_ALL == 1, 1, 0)
  data_full$astro_pcr <- ifelse(data_full$F19_ASTRO_VIRUS == 1, 1, 0)
  data_full$noro_pcr <- ifelse(data_full$NOROVirus == 1, 1, 0)
  data_full$tEPEC_culture <- ifelse(data_full$tEPEC == 1, 1, 0)
  data_full$campy_culture <- ifelse(data_full$Campy == 1, 1, 0)
  data_full$sapo_pcr <- ifelse(data_full$F19_SAPO_VIRUS == 1, 1, 0)
  data_full$EAEC_culture <- ifelse(data_full$EAEC.x == 1, 1, 0)
  data_full$giardia_elisa <- ifelse(data_full$F18_RES_GIARDIA == 1, 1, 0)
  
  ## SHIGELLA - Shigella attributable diarrhea by AFE
  data_full$shigella_attributable_tac <- ifelse(data_full$shigella_eiec_afe > 0.5, 1, 0)
  
  ## SHIGELLA - Define Shigella culture positive
  data_full$shigella_culture_pos <- ifelse(
    rowSums(data_full[, c(
      "F16CORR_S_DYSENT", "F16CORR_S_FLEXNERI", "F16CORR_S_BOYDII", 
      "F16CORR_S_SONNEI", "F16CORR_S_NONTYPE", "F16CORR_SHIG_1A", 
      "F16CORR_SHIG_1B", "F16CORR_SHIG_2A", "F16CORR_SHIG_2B", 
      "F16CORR_SHIG_3A", "F16CORR_SHIG_3B", "F16CORR_SHIG_4A", 
      "F16CORR_SHIG_4B", "F16CORR_SHIG_4C", "F16CORR_SHIG_5A", 
      "F16CORR_SHIG_5B", "F16CORR_SHIG_6", "F16CORR_SHIG_X", 
      "F16CORR_SHIG_Y", "F16CORR_SHIG_NONTYP"
    )] == 1, na.rm = TRUE) > 0, 1, 0)
  
  
  # If NA tac results, make sure it is NA (we only want to use people who have AFE data)
  data_full$shigella_attributable <- ifelse(
    data_full$shigella_attributable_tac == 1 | data_full$shigella_culture_pos == 1,
    1,
    0
  )
  
  data_full$shigella_attributable <- ifelse(is.na(data_full$shigella_attributable_tac), NA, data_full$shigella_attributable)
  
  # CASES -- define cases vs controls
  if(case_def == "tac_or_culture_shig_diar"){
    
    # Type == Case means they have diarrhea & TAC or culture attributable Shigella
    data_full$case <- ifelse(data_full$Type == "Case", 1, 0)
    
    # Drop people who have diarrhea (are cases) but it's not TAC attributable.
    data_full <- data_full[-which(data_full$case == 1 & (data_full$shigella_attributable == 0 | is.na(data_full$shigella_attributable))),]
    
    case_ids <- unique(data_full$child_id[which(data_full$case == 1)])
    
    # only keep the matched controls for the subset of cases with TAC attributable Shigella
    # note the matching is not always 1:1
    data_full <- data_full[-which(data_full$case == 0 & !(data_full$case_id %in% case_ids)),]
    
    
  } else if(case_def == "tac_shig_diar"){
    # Type == Case means they have diarrhea & TAC attributable Shigella
    data_full$case <- ifelse(data_full$Type == "Case", 1, 0)
    
    # Drop people who have diarrhea (are cases) but it's not TAC attributable.
    data_full <- data_full[-which(data_full$case == 1 & (data_full$shigella_attributable_tac == 0 | is.na(data_full$shigella_attributable_tac))),]
    
    
    case_ids <- unique(data_full$child_id[which(data_full$case == 1)])
    
    # only keep the matched controls for the subset of cases with TAC attributable Shigella
    # note the matching is not always 1:1
    data_full <- data_full[-which(data_full$case == 0 & !(data_full$case_id %in% case_ids)),]
    
  } else if(case_def == "culture_shig_diar"){
    
    # Type == Case means they have diarrhea
    data_full$case <- ifelse(data_full$Type == "Case", 1, 0)
    
    # Drop people who have diarrhea but it's not culture attributable (OR IF SHIGELLA ATTRIBUTABLE TAC IS NA)
    data_full <- data_full[-which(data_full$case == 1 & (data_full$shigella_culture_pos == 0 | is.na(data_full$shigella_attributable_tac))),]
    
    case_ids <- unique(data_full$child_id[which(data_full$case == 1)])
    
    # only keep the matched controls for the subset of cases with TAC attributable Shigella
    # note the matching is not always 1:1
    data_full <- data_full[-which(data_full$case == 0 & !(data_full$case_id %in% case_ids)),]
    
    
  } else if(case_def == "all_diar"){
    
    # Keep all diarrhea cases and controls
    
    # Type == Case means they have diarrhea
    data_full$case <- ifelse(data_full$Type == "Case", 1, 0)
    
    # case_ids <- unique(data_full$child_id[which(data_full$case == 1)])
    
    # only keep the matched controls for the subset of cases with TAC attributable Shigella
    # note the matching is not always 1:1
    # data_full <- data_full[-which(data_full$case == 0 & !(data_full$case_id %in% case_ids)),]
    
    
  } else{
    stop("Case definition not supported")
  }
  
  ## OUTCOME - Change in HAZ enrollment --> day 60
  # data_full$hazdiff <- data_full$F4B_MED_HAZ - data_full$F5_HAZ
  # day 60 = F5_HAZ for both case and control
  # enrollment = F4B_MED_HAZ for cases, F7 for controls
  
  # eliminate any extreme obs first
  data_full$F5_HAZ <- ifelse(data_full$F5_HAZ < -6 | data_full$F5_HAZ > 6, NA, data_full$F5_HAZ)
  data_full$F7_MED_HAZ <- ifelse(data_full$F7_MED_HAZ < -6 | data_full$F7_MED_HAZ > 6, NA, data_full$F7_MED_HAZ)
  data_full$F4B_MED_HAZ <- ifelse(data_full$F4B_MED_HAZ < -6 | data_full$F4B_MED_HAZ > 6, NA, data_full$F4B_MED_HAZ)
  
  data_full$F5_WHZ <- ifelse(data_full$F5_WHZ < -6 | data_full$F5_WHZ > 6, NA, data_full$F5_WHZ)
  data_full$F7_MED_WHZ <- ifelse(data_full$F7_MED_WHZ < -6 | data_full$F7_MED_WHZ > 6, NA, data_full$F7_MED_WHZ)
  data_full$F4B_MED_WHZ <- ifelse(data_full$F4B_MED_WHZ < -6 | data_full$F4B_MED_WHZ > 6, NA, data_full$F4B_MED_WHZ)
  
  data_full$F5_WAZ <- ifelse(data_full$F5_WAZ < -6 | data_full$F5_WAZ > 6, NA, data_full$F5_WAZ)
  # WAZ does not have the med variables
  #data_full$F7_MED_WAZ <- ifelse(data_full$F7_MED_WAZ < -6 | data_full$F7_MED_WAZ > 6, NA, data_full$F7_MED_WAZ)
  #data_full$F4B_MED_WAZ <- ifelse(data_full$F4B_MED_WAZ < -6 | data_full$F4B_MED_WAZ > 6, NA, data_full$F4B_MED_WAZ)
  data_full$F7_WAZ <- ifelse(data_full$F7_WAZ < -6 | data_full$F7_WAZ > 6, NA, data_full$F7_WAZ)
  data_full$F4B_WAZ <- ifelse(data_full$F4B_WAZ < -6 | data_full$F4B_WAZ > 6, NA, data_full$F4B_WAZ)
  
  data_full$hazdiff <- ifelse(data_full$case == 1, 
                              data_full$F5_HAZ - data_full$F4B_MED_HAZ,
                              data_full$F5_HAZ - data_full$F7_MED_HAZ)
  
  data_full$whzdiff <- ifelse(data_full$case == 1,
                              data_full$F5_WHZ - data_full$F4B_MED_WHZ,
                              data_full$F5_WHZ - data_full$F7_MED_WHZ)
  
  data_full$wazdiff <- ifelse(data_full$case == 1, 
                              data_full$F5_WAZ - data_full$F4B_WAZ,
                              data_full$F5_WAZ - data_full$F7_WAZ)
  
  data_full$hazd60 <- data_full$F5_HAZ
  data_full$whzd60 <- data_full$F5_WHZ
  data_full$wazd60 <- data_full$F5_WAZ
  
  ## ANTIBIOTICS - vars for any abx, who approved abx
  data_full$any_abx <- 0
  
  data_full$who_abx <- 0
  data_full$maybe_eff_abx <- 0
  data_full$ineff_abx <- 0
  data_full$no_abx <- 0
  data_full$azithro <- 0
  data_full$cipro <- 0 

  who_approved_names <- c("ARITHROMYCIN", 
                          "ASITHROMYCIN",
                          "ASITHROMYCIN, ERYTHROMYCIN",
                          "AZIMOX & AZITHRRO",
                          "AZITHROMYCIN",
                          "AZITHROMYCIN, ERYTHROMYCIN",
                          "AZITHROMYEIN",
                          "ASITHROMYCINN, NALIDIXICACID",
                          "AZITHTOMYCIN",
                          "CEFTRIAXONE",
                          "CEFTRIAXONE,AMOXICILLINE",
                          "CEFTRAXONE",
                          "CEFTRIAZONE",
                          "CIFROFLEXACIN",
                          "CIPROFLACIN",
                          "CIPROFLAXACIN",
                          "CIPROFLAXZCIN",
                          "CIPROFLEXACIN",
                          "CIPROFLOXACIN",
                          "CIPROFLOXAEIN",
                          "CIPROFLOXCIN",
                          "CIPROFLOXICIN",
                          "CIPROFLOZACIN",
                          "CIPROFLXACIN",
                          "CIPROFOLOXACIN",
                          "CIPROXIN",
                          "PIVMECILLINAM",  # note not seeing pivmecillinam in data
                          "TAB CIPROFLOXACIN",
                          "INJ CEFTRIAXONE.",
                          "INJ. CEFTRIAXONE (250MG)",
                          "INJ. CEFTRIAXONE 250MG",
                          "INJ. CEFTRIAXONE",
                          "CRFTRIAZONE",
                          "INJ CEFXONE") 
  
  maybe_effective_abx <- c("AMIKACIN", 
                           "AMOXICILLIN", 
                           "AMPICILLIN", 
                           "AUGMENTIN/CO-AMOXICLAV", 
                           "CEFACOR", 
                           "CEFIXIME", 
                           "CEFPODOXIME", 
                           "CEFTRAZIDIME", 
                           "CEFTRIAXONE", 
                           "CEFTRUOXIME", 
                           "CHLORAMPHENICOL", 
                           "CLARITHROMYCIN", 
                           "COTRIMOXAZOLE", 
                           "DOXYCYCLINE", 
                           "ERYTHROMYCIN", 
                           "FURAZOLIDONE", 
                           "GENTAMICIN", 
                           "LEVOFLOXACIN", 
                           "MEROPENEM", 
                           "NALIDIXIC ACID", 
                           "NIFUROXAZIDE", 
                           "NITROFURANTOIN", 
                           "OTHER MACROLIDES", 
                           "STREPTOMYCIN", 
                           "TETRACYCLINE", 
                           "THIAMPHENICOL",
                           "TETRACYCLINE CAP",
                           "TERACYCLINE",
                           "TETRACYCLINE CAPS" ,                                       
                           "TETRACYCLINE CAPSUL" ,
                           "TETRACYCLINE CAPSULE",
                           "AMOXACILLAN SYRUP",
                           "AMOXLYL" ,
                           "AMOXYL" ,
                           "AMOXYL SYRUP" ,
                           "CEFIXINE,COTRIMOXAZOL",
                           "COTRIMOXOZOL",
                           "AMOXICILINE",
                           "CO TRIMOXAZOLE",
                           "CCOTRIMOXAZOLE",
                           "CONTRIMOXAZOLE",
                           "CORTIMOXAZOLE",
                           "COTRINOKAZOLE",
                           "COTRIMOXAZOLE, METRONIDAZOLE", #check bc one effective, one ineffective
                           "AMPICILINE",
                           "BIODROXIL/ COTRIMOXAZOLE",
                           "COTRIMOXAZ0LE",
                           "COTRI",
                           "AMOXICYLLINE, METRONIDAZOLE",
                           "AMOXACILLIN",
                           "AMOXICILLINE SP",
                           "ADVENT",
                           "AMOXIL SYRUP" ,
                           "AMOXIL TAB.",
                           "BACITRACIN (COTRINOXAZOLE)",
                           "AMOSICILLIN",
                           "AMOXACILLY",                                             
                           "AMOXICILLIN,METRONIDAZOLE",
                           "AMOXICILLINE",
                           "AMOXOCILLAN",                                              
                           "AMOXOCILLEN",
                           "AMOXOCILLIN",                                              
                           "AMOXYCILLAN",
                           "AMOXYCILLIN",                                              
                           "AMOXYCILLIN & CLAVULANIC ACID",
                           "AMOXICILLINE, COTRIMOXAZOLE",
                           "BIODROXYL,AMOXICILLINE",
                           "CAPSULES AMOX",
                           "PARACETAMOL, AMOXICILLIN",
                           "PRINCIMOX",
                           "CHLORANPHENICOL",
                           "CHLORAMPHENICOLE",
                           "CLORAMPHENICOLE",
                           "CLORAFENICOL",
                           "CLOTRIMOXAZOLE",
                           "CO-TREMOXAZOLE" ,
                           "CO-TRIMOXAZOLE",
                           "COTRIMOXAZOL",
                           "CO-TRIMOXAZOLE, METRONIDAZOLE",
                           "CORTRIMOXAZOLE",
                           "SEPTRIN",                     # brand name?
                           "SEPTRIN (COTRIMOXAZOLE" ,
                           "SEPTRIN SYRP.",
                           "SEPTRIN SYRUP",
                           "SEPTRIN TAB." ,
                           "SEPTRIN TABS" ,
                           "SEPTRINE SYRUP",
                           "SEPTRINE SYRUP.",
                           "SYTUP SEPTRINE",  
                           "TRIMETHOPRIM & SULFAHENAZOLE",
                           "TRIMETHOPRIME SULFAMETHAZOLE" ,
                           "UCLAPRIM",
                           "SYP AMOXYCILLIN",
                           "SYP. COTRIMOXAZOLE",
                           "SYP AMOXYCILLIN",
                           "ERTHROMYCIN",
                           "ERTHROMYCIN",
                           "EROIHROMYCIN",
                           "ERYTHROMYCINE",
                           "ERYTHRONYCIN",
                           "ERCYTHCOMYSH",
                           "ERYTHNOMYCIN",                                             
                           "ERYTHOMYCINE" ,
                           "ERYTHROMYCIN, METRONIDAZOLE",
                           "ERYTHROMYCM",
                           "ERYTHROMYUME",
                           "ERYTROMICINE",
                           "ERYTROMIXYNE",
                           "SYP. ERYTHROMYCINE",
                           "CORTRIMOXAZOLE, ERYTHROMYCIN" ,
                           "COTRINOXAZOLE",
                           "COTRTMOXAZOLE",
                           "CORIMOXAZOLE",
                           "COTAMOXAZOLE" ,
                           "COTAMOXAZOLE",                                             
                           "COTRIMASCAZOLE",
                           "COTRIMOSCAZOLE",
                           "COTRIMOTAZOLE" ,
                           "COTRIMOXAZOLE & METRONIDAZOLE",
                           "COTRIMOXAZOLE SYRUP",
                           "COTRIMOXOZOLE",   
                           "COFRIMUXAZOLE",
                           "SEPRINE COTRIMOXASOLE",
                           "COMMOXAZOLE" ,
                           "COTRIMZAZOLE SYRUP",
                           "NALIDEXIC ACID",
                           "NADIDIXIC ACID",
                           "NALEDIXIC ACID",
                           "NALIDENIC ACID",
                           "NALIDENIE ACID",
                           "NALIDERIC ACID",
                           "NALIDEXIE ACID",
                           "NALIDIXIE ACID",
                           "NALIDXIC ACID",
                           "NALIDENE ACID",
                           "NALIDOXIC ACID",
                           "SYP. ENTAMEZOLE (SYP. METRONIDAZOL- OLE & NALIDIXIC ACID)",
                           "SYP. NEGRAM" ,
                           "METRONIDAZOLE, NALIIDIXIC ACID" ,
                           "SYP NALIDIX ACID",
                           "CEFIXINE",
                           "CEFIXIME SYRUP",
                           "CEFEXINE FRIHYDRATE",
                           "CEFIXIME SODIUM.",
                           "SYP. SEFEXIME & I/M INJ.",
                           "INJ CEFIXIME",
                           "SYP CEFIXME",
                           "MAXIMA (CEFIXIME)",
                           "TAXIM-O",
                           "DOXY",
                           "DOXY CYCLINE",
                           "HICANCIL",
                           "INJ-GENTAMICIN",
                           "O FLOXACIN",
                           "OFLAXACIN",
                           "XINTOF - M (OFLOXACIN & METRONIDAZOLE)" ,
                           "OFLOX & ORNIDAZOLE",
                           "OFLOXAC",
                           "OFLOXACAN",
                           "XINTOF-M",
                           "OFLOXACIN",
                           "OFLOXACIN & METRONIDAZOLE",
                           "OFLOXACIN & ORINIDAZOLE",
                           "OFLOXACIN & ORNIDAZNE",
                           "OFLOXACIN & ORNIDAZOLE",
                           "OFLOXACIN AMIKACIN",
                           "OFLOXACIN AND METRONIDAZOLE",
                           "OFLOXACIN AND ORNIDAZOLE",
                           "OFLOXACIN SUSPENSION",
                           "OFLOXACIN-HETROGYL(METRONIDA)",
                           "OFLOXACIN, ORNIDAZOLE",
                           "OFLOXACINE & ORNIDAZOLE",
                           "OFLOXACINS & METRONIDAZOLE",
                           "OFOXACIN & ORNIDAZOLE",
                           "SYR OFLONAC (OFLOXACIN)" ,
                           "SUSP. OFLOXACIN",
                           "SYR OFLONAC (OFLOXACIN)"  ,
                           "02 SUSPENSION (OFLOXACIN ORNIDAZOLE)",
                           "02 (OOFLOXACIN)",
                           "SUSP. ZENFLOX"  ,
                           "SYR TIDFLOX",
                           "ZENFLOX",
                           "ORNIDAZOLE SUSP. & OFLOXACIN SUSP.",
                           "SEPTIN",
                           "SEPTRAN",
                           "SEPTRIN FLAGYL",
                           "SEPTRIN- AMOXIL",
                           "SEPTRINE",
                           "SYP SEPTRAN" ,
                           "SEPTRIN AND FLAGYL",
                           "SEDTRIN",
                           "SEPTRUM,FLAGYL SYRUP",
                           "SEPRIN SYRUP",
                           "SULFAMETHOXAZOLE",
                           "BACTRIM 200/40MG",
                           "BACTRIM SP",
                           "BACTRIM 200/40MG",
                           "BACTRIM SP",
                           "WALLAMYCIN",
                           "SYR. WALAMYCIN",
                           "NOFLAIACIN",
                           "NORFLOXACIN",
                           "NORFLOXACIN & METRONIDAZOLE",
                           "NORFLOXACIN, METRONIDAZOLE",
                           "NORFLOXACIN, TINDAZOLE",
                           "NORFLOXACIN,METRONIDAZOLE",
                           "NAROFLOXACIN"  ,
                           "NORMET",
                           "NORMET SYR.",
                           "SYR. NORMET.",
                           "SYRUP NOMMET",
                           "DIARYL" ,
                           "NI FUROXAZIDE",
                           "AMPILOX",
                           "AMPICILLINE" ,
                           "AMOXILLINE",
                           "CEFOTAXIME",
                           "CEFOTAMIME",
                           "CEFOTAXIME",
                           "CEFOTAXIME(250MG)",
                           "INJ-CEFTAXIME",
                           "INJ. CEFOTAXIME",
                           "INJ. CEFOTAXIME 500MG",
                           "INJ/CEFOTAXIME" ,
                           "CETRIMAXAZOLE SYRUP",
                           "CURAM",
                           "COLISTIN",
                           "FOSOMIN(FOSFOMYCIN)",
                           "INJ CEFTROXIME" ,
                           "MONOCEF",
                           "NOR-METROGL.",
                           "OSPAMOX",
                           "SUL FAGMANIDINE",
                           "SULEAMIDE",
                           "SULFAGUANIDINE",
                           "SULFAGUANIDINE CP.",
                           "SULFAMIDE",
                           "CEFPODOXINE PROXELIL",
                           "CEFPODOXIME PROCEFIL",
                           "CAFOTAXIME")
  
  # this contains abx we know are ineffective. 
  # the ineffective variable will also contain people with write-ins that appear to be typos, "unknown", and no abx
  ineffective_abx <- c("CEFADROXIL", 
                       "CLINDAMYCIN", 
                       "CLOXACILLIN", 
                       "ENTAMIZOLE", 
                       "FLUCLOXACILLIN", 
                       "METRONIDAZOLE (FLAGYL)", 
                       "MUPIROCIN", 
                       "PENICILLIN", 
                       "PYRAZINAMIDE",
                       "FLAGYL",
                       "FLAGYL SYRUP",
                       "SRP FLAGYL" ,
                       "SYP. FLAGYL",
                       "FLAGYL, AMOXYL",
                       "FLAGYL, METRONDAZOLE",
                       "METRONIDAZOLE",
                       "METRONIDAZOLE SYRUP",
                       "METRONIDAZOL",
                       "METRCORIDAZOLE",
                       "METROINDAZOLE",
                       "METROMDAZOLE",
                       "METROMIDAZOLE",
                       "METROMIDAZOLE SYRU[",
                       "METRONDAZOLE SYP.",
                       "METRONEDAZOLE",
                       "METRONIDAJOLE",
                       "METRONIDAZOLE & COTIMOXAZOL",
                       "METRONIDAZOLE & DILOXANIDE FUROATE",
                       "METRONIDAZOLE AND FURAZOLIDONE",
                       "METRONIDAZOLE SUSPENSION",
                       "METRONIDAZOLE SYR",
                       "NONFLOXACIN & METROGYL SYRUP.",
                       "SYR. METROGYL",
                       "DELOXACIN & METRONIDAZOLE",
                       "CEFADROLIL",
                       "FAGYL",
                       "FLAGGYL TABLETS",                                          
                       "FLAGLY",
                       "FLAGYL SYP",
                       "FLAGYL TABLET",
                       "FLAGYL TABLETS.",                                          
                       "FLAGYL TABS",
                       "SYP FLAGYL",
                       "SYP FLEGYL",
                       "FLAYGILE" ,
                       "FLAGYL, CEPHRADINE, EEFAZOLIN",
                       "FLAGYL, METRONIDAZOLE"   ,
                       "SYP METROINIDAZOLE",
                       "SYP METRONIDAZOLE" ,
                       "SYP METRONIDAZOL",
                       "SYP. METRONIDAZOLE",
                       "SYRUP METRONIDAZOLE",
                       "MEHROGYL",
                       "METHONIDAZOLE" ,
                       "METOQUIDAZOLE",
                       "METRCOAIDAZOLE",                                           
                       "METRCORIDARZODE" ,
                       "METROGYL",                                                
                       "METROGYL, OFLOXACIN" ,
                       "METRONIDAZOLE, CEFIXIME",
                       "METRONIDOZAL",
                       "METRONIDOZOL",
                       "METRONIDOZOLE",
                       "METRONISAZOLE",
                       "METROZINE",
                       "METRUNIDAZEL",
                       "NORFLEX-METRONIDAZOLE",
                       "BIOBROXIL",
                       "BIODROXIL",
                       "BOIDROXIL",
                       "ENTAMEZOLE",
                       "GRAMAGYL",
                       "PARACETAMOL",
                       "SMECTA" ,
                       "FLUCAZOL",
                       "ENTEROQUINOL",
                       "CEFALEXINE",
                       "CEFAZOLIN SODIUM/IV" ,
                       "CEFLEX",
                       "CEPHALEXIN",
                       "SYP KEFLEX CEPHALEXIN",
                       "SPORNIDEX",
                       "CEPHESEIN",
                       "DELOXACINE, ORNIDAZOLE,",
                       "NITAZOXAMIDE",
                       "NITAZOXANIDE",
                       "SYNOBE",
                       "SYP(ENTAMIZOLE)",
                       "CEPHRADINE",
                       "(NOT ABLE TO)RED AND BLUE CAPSELS.",
                       "UNKNOWN",
                       "UNKNOW CAPSUL",
                       "UNKNOWN CAPSULES",
                       "UNKNOWN CAPSUL",
                       "SURUP",
                       "CATMAX",
                       "NETRO",
                       "GANDIDA",
                       "GANIDA")
  
  azithro_names <- c("ARITHROMYCIN", 
                     "ASITHROMYCIN",
                     "ASITHROMYCIN, ERYTHROMYCIN",
                     "AZIMOX & AZITHRRO",
                     "AZITHROMYCIN",
                     "AZITHROMYCIN, ERYTHROMYCIN",
                     "AZITHROMYEIN",
                     "ASITHROMYCINN, NALIDIXICACID",
                     "AZITHTOMYCIN")
  
  # ----- Antibiotics pre-enroll -----
  
  # No longer including hometrt antibiotics 6/11/26
  
  # Any abx = hometrt abx marked 1
  # data_full$any_abx <- ifelse(data_full$F4A_HOMETRT_AB == 1, 1, data_full$any_abx)
  # 
  # # WHO recommended abx
  # data_full$who_abx <- ifelse(data_full$F4A_HOMETRT_AB_SPEC.x %in% who_approved_names, 1, data_full$who_abx)
  # 
  # # Maybe effective antibiotics
  # data_full$maybe_eff_abx <- ifelse(data_full$F4A_HOMETRT_AB_SPEC.x %in% maybe_effective_abx, 1, data_full$maybe_eff_abx)
  # 
  # # Ineffective abx
  # data_full$ineff_abx <- ifelse(data_full$F4A_HOMETRT_AB_SPEC.x %in% ineffective_abx, 1, data_full$ineff_abx)
  # 
  # # azithro
  # data_full$azithro <- ifelse(data_full$F4A_HOMETRT_AB_SPEC.x %in% azithro_names, 1, data_full$azithro)
  
  # Antibiotics hospital or home
  
  # cotrimoxazole - TRT_GIVE_CXL | TRT_PRES_CXL
  # gentamycin - TRT_GIVE_GENT | TRT_PRES_GENT
  # chloramphenicol - TRT_GIVE_CHLOR | TRT_PRES_CHLOR
  # erythromycin - TRT_GIVE_ERY | TRT_PRES_ERY
  # azithromycin - TRT_GIVE_AZI | TRT_PRES_AZI
  # other macrolides - TRT_GIVE_MACR | TRT_PRES_MACR
  # penicillin - TRT_GIVE_PEN | TRT_PRES_PEN
  # amoxycillin - TRT_GIVE_AMOX | TRT_PRES_AMOX
  # ampicillin - TRT_GIVE_AMPI | TRT_PRES_AMPI
  # nalidixic acid - TRT_PRES_NALID | TRT_GIVE_NALID
  # ciprofloxacin - TRT_PRES_CPNR | TRT_GIVE_CPNR
  # selexid/pivmecillinam - TRT_PRES_SLPY | TRT_GIVE_SLPY
  # other antibiotic - TRT_PRES_OTHR | TRT_GIVE_OTHR
  
  # Any antibiotics for cases
  data_full$any_abx <- ifelse(
    rowSums(
      cbind(
        data_full$F4B_TRT_GIVE_CXL, data_full$F4B_TRT_PRES_CXL.x,
        data_full$F4B_TRT_GIVE_GENT, data_full$F4B_TRT_PRES_GENT.x,
        data_full$F4B_TRT_GIVE_CHLOR, data_full$F4B_TRT_PRES_CHLOR.x,
        data_full$F4B_TRT_GIVE_ERY, data_full$F4B_TRT_PRES_ERY.x,
        data_full$F4B_TRT_GIVE_AZI, data_full$F4B_TRT_PRES_AZI.x,
        data_full$F4B_TRT_GIVE_MACR, data_full$F4B_TRT_PRES_MACR.x,
        data_full$F4B_TRT_GIVE_PEN, data_full$F4B_TRT_PRES_PEN.x,
        data_full$F4B_TRT_GIVE_AMOX, data_full$F4B_TRT_PRES_AMOX.x,
        data_full$F4B_TRT_GIVE_AMPI, data_full$F4B_TRT_PRES_AMPI.x,
        data_full$F4B_TRT_GIVE_NALID, data_full$F4B_TRT_PRES_NALID.x,
        data_full$F4B_TRT_GIVE_CPNR, data_full$F4B_TRT_PRES_CPNR.x,
        data_full$F4B_TRT_GIVE_SLPY, data_full$F4B_TRT_PRES_SLPY.x,
        data_full$F4B_TRT_GIVE_OTHR, data_full$F4B_TRT_PRES_OTHR.x
      ) == 1, na.rm = TRUE
    ) > 0, 1, data_full$any_abx
  )
  
  # Any antibiotics for controls
  data_full$any_abx <- ifelse(
    rowSums(
      cbind(
        data_full$F7_MED_GENT, data_full$F7_MED_COTR, data_full$F7_MED_CHLOR,
        data_full$F7_MED_ERYTH, data_full$F7_MED_AZITH, data_full$F7_MED_OMACR,
        data_full$F7_MED_PENI, data_full$F7_MED_AMOXY, data_full$F7_MED_AMPI,
        data_full$F7_MED_NALID, data_full$F7_MED_CIPRO, data_full$F7_MED_SELE,
        data_full$F7_MED_OTHERANT
      ) == 1, na.rm = TRUE
    ) > 0, 1, data_full$any_abx
  )
  
  # note ceft does not have its own category
  # WHO antibiotics for cases
  data_full$who_abx <- ifelse(
    rowSums(
      cbind(
        data_full$F4B_TRT_GIVE_AZI, data_full$F4B_TRT_PRES_AZI.x,
        data_full$F4B_TRT_GIVE_CPNR, data_full$F4B_TRT_PRES_CPNR.x,
        data_full$F4B_TRT_GIVE_SLPY, data_full$F4B_TRT_PRES_SLPY.x
      ) == 1, na.rm = TRUE
    ) > 0, 1, data_full$who_abx
  )
  
  # WHO antibiotics for controls
  data_full$who_abx <- ifelse(
    rowSums(
      cbind(
        data_full$F7_MED_AZITH,
        data_full$F7_MED_CIPRO,
        data_full$F7_MED_SELE
      ) == 1, na.rm = TRUE
    ) > 0, 1, data_full$who_abx
  )
  
  # Maybe effective antibiotics for cases
  data_full$maybe_eff_abx <- ifelse(
    rowSums(
      cbind(
        data_full$F4B_TRT_GIVE_CXL, data_full$F4B_TRT_PRES_CXL.x,
        data_full$F4B_TRT_GIVE_GENT, data_full$F4B_TRT_PRES_GENT.x,
        data_full$F4B_TRT_GIVE_CHLOR , data_full$F4B_TRT_PRES_CHLOR.x,
        data_full$F4B_TRT_GIVE_ERY, data_full$F4B_TRT_PRES_ERY.x,
        data_full$F4B_TRT_GIVE_MACR ,data_full$F4B_TRT_PRES_MACR.x,
        data_full$F4B_TRT_GIVE_AMOX, data_full$F4B_TRT_PRES_AMOX.x,
        data_full$F4B_TRT_GIVE_AMPI ,data_full$F4B_TRT_PRES_AMPI.x,
        data_full$F4B_TRT_GIVE_NALID,data_full$F4B_TRT_PRES_NALID.x,
        data_full$F4B_TRT_GIVE_OTHR,data_full$F4B_TRT_PRES_OTHR.x
      ) == 1, na.rm = TRUE
    ) > 0, 1, data_full$maybe_eff_abx
  )
  
  # Ineffective:
  data_full$ineff_abx <- ifelse(
    rowSums(
      cbind(
        data_full$F4B_TRT_GIVE_PEN, data_full$F4B_TRT_PRES_PEN.x
      ) == 1, na.rm = TRUE
    ) > 0, 1, data_full$ineff_abx)
  
  # make sure only fall into one category, taking highest
  data_full$ineff_abx <- ifelse(data_full$ineff_abx == 1 & (data_full$maybe_eff_abx == 1 | data_full$who_abx == 1), 0, data_full$ineff_abx)
  data_full$maybe_eff_abx <- ifelse(data_full$who_abx == 1 & data_full$maybe_eff_abx == 1, 0, data_full$maybe_eff_abx)
  
  # if ineff = 0, maybe = 0, and who = 0, categorize as no abx
  data_full$no_abx <- ifelse(data_full$who_abx == 0 & data_full$maybe_eff_abx == 0 & data_full$ineff_abx == 0, 1, data_full$no_abx)
  
  # Make 'all abx' variable so levels all in one
  data_full$all_abx <- ifelse(data_full$no_abx == 1 | data_full$ineff_abx == 1, 0, 
                              ifelse(data_full$maybe_eff_abx == 1, 1, 2))
  data_full$all_abx <- factor(data_full$all_abx, levels = 0:2, labels = c("No/Ineffective abx", "Maybe effective abx", "WHO approved abx"))
  
  # azithro specifically
  data_full$azithro <- ifelse(data_full$F4B_TRT_GIVE_AZI == 1 | data_full$F4B_TRT_PRES_AZI.x == 1, 
                                   1, data_full$azithro)
  
  # cipro
  data_full$cipro <- ifelse(
    data_full$F4B_TRT_GIVE_CPNR == 1 | data_full$F4B_TRT_PRES_CPNR.x == 1, 1, data_full$cipro
  )
  
  # ceft not explicitly in GEMS
  
  # Drop controls with abx (healthy controls only)
  data_full <- data_full[-which(data_full$case == 0 & data_full$any_abx == 1),]
  
  ## COVARIATES & SEVERITY VARS
  
  # create variable for days between enrollment and 60-day follow-up
  data_full$followupdate <- as.Date(data_full$F5_DATE)
  data_full$followup_days <- data_full$followupdate - as.Date(data_full$ENROLLDATE)
  data_full$I_followup_days <- ifelse(is.na(data_full$followup_days), 0, 1)
  data_full$I_followup_days_x_followup_days <- ifelse(is.na(data_full$followup_days), 0, data_full$followup_days)
  
  # get other variables of interest
  
  # covariates:
  # sex (F2_GENDER - 0 boy, 1 girl for cases, F6 for controls)
  data_full$sex <- ifelse(data_full$case == 1, 
                          data_full$F2_GENDER,
                          data_full$F6_GENDER)
  
  # age (F2_age - age in months cases, F6 controls)
  data_full$age <- ifelse(data_full$case == 1,
                          data_full$F2_AGE,
                          data_full$F6_AGE)
  
  # HAZ (F4B_MED_HAZ - height-for-age z-score at enrollment)
  data_full$enr_haz <- ifelse(data_full$case == 1,
                              data_full$F4B_MED_HAZ,
                              data_full$F7_MED_HAZ)
  
  data_full$enr_whz <- ifelse(data_full$case == 1,
                              data_full$F4B_MED_WHZ,
                              data_full$F7_MED_WHZ)
  
  data_full$enr_waz <- ifelse(data_full$case == 1,
                              data_full$F4B_WAZ,
                              data_full$F7_WAZ)
  
  # site  (SITE - site - 1-7 but not seeing dictionary)
  
  # primary_caretaker_edu (F4A_PRIM_SCHL - ordinal 1-5, 6 religious only, 7 NA)
  data_full$prim_caregiver_edu <- ifelse(data_full$case == 1, 
                                         data_full$F4A_PRIM_SCHL,
                                         data_full$F7_PRIM_SCHL)
  
  # num_hh_lt5 (F4A_YNG_CHILDREN - how many children younger than 60 months live in household)
  data_full$num_hh_lt5 <- ifelse(data_full$case == 1,
                                 data_full$F4A_YNG_CHILDREN, 
                                 data_full$F7_YNG_CHILDRN)
  
  # breastfed - (F4A_BREASTFED - is child currently breastfed - 0 no, 1 partial, 2 exclusive)
  data_full$breastfed <- ifelse(data_full$case == 1,
                                data_full$F4A_BREASTFED,
                                data_full$F7_BREASTFED)
  
  ### NEW 9/5/25 -- adding AMR data
  
  amr <- amr %>%
    filter(Pathotype == "Shigella") %>%
    mutate(SID = as.numeric(SID)) %>%
    select("SID", "mphA") %>%
    mutate(mphA = if_else(mphA == 0, "S", "I/R")) %>%
    rename('S_IR_azithro' = mphA)
  
  data_full <- left_join(data_full, amr, by = "SID") 
  
  data_full$azithro_given_S_IR <- ifelse(data_full$azithro == 1, data_full$S_IR_azithro, NA)
  data_full$ast_given_abx <- data_full$azithro_given_S_IR
  
  # create overall resistant and susceptible variables w same logic as EFGH
  data_full$resistant_WHO_approve <- apply(
    data_full[, c("S_IR_azithro")],
    1,
    function(x) if (any(x == "I/R", na.rm = TRUE)) 1 else 0
  )
  
  data_full$susceptible_WHO_approve <- apply(
    data_full[, c("S_IR_azithro")],
    1,
    function(x) if (any(x == "S", na.rm = TRUE)) 1 else 0
  )
  
  
  # --- NEW 12/9/24 FOR SES TO MATCH VIDA: ---
  # safe_water
  # safe_sanit
  
  # Make joint monitoring programme definition of improved water
  # based on data dictionary for VIDA
  data_full <- data_full %>%
    mutate(
      MS_WATER = coalesce(F4A_MS_WATER, F7_MS_WATER),
      TIME_WATER = coalesce(F4A_TIME_WATER, F7_TIME_WATER),
      WATER_AVAIL = coalesce(F4A_WATER_AVAIL, F7_WATER_AVAIL),
      MS_SPEC = coalesce(F4A_MS_SPEC, F7_MS_SPEC) 
    ) %>%
    mutate(
      safe_water = case_when(
        (MS_WATER %in% c(1, 2, 17, 7, 8, 9, 10, 11, 15, 16, 3) | grepl("tap|pipe", MS_SPEC, ignore.case = TRUE)) &
          TIME_WATER %in% c(1, 2) & WATER_AVAIL %in% c(1, 2) ~ 1,
        (MS_WATER %in% c(1, 2, 9, 15) | grepl("pipe", MS_SPEC, ignore.case = TRUE)) &
          is.na(TIME_WATER) & WATER_AVAIL %in% c(1, 2) ~ 1,
        (MS_WATER %in% c(1, 2, 17, 7, 8, 9, 10, 11, 15, 16, 3) | grepl("tap|pipe", MS_SPEC, ignore.case = TRUE)) &
          TIME_WATER %in% c(1, 2) & WATER_AVAIL %in% c(3, 4) ~ 2,
        (MS_WATER %in% c(1, 2, 9, 15) | grepl("pipe", MS_SPEC, ignore.case = TRUE)) &
          is.na(TIME_WATER) & WATER_AVAIL %in% c(3, 4) ~ 2,
        (MS_WATER %in% c(1, 2, 17, 7, 8, 9, 10, 11, 15, 16, 3) | grepl("tap|pipe", MS_SPEC, ignore.case = TRUE)) &
          TIME_WATER %in% c(3, 4, 5) ~ 3,
        MS_WATER %in% c(4, 5, 12) | MS_SPEC == "open" ~ 4,
        MS_WATER %in% c(6, 13, 14, 18, 19) ~ 5,
        MS_WATER == 18 & !grepl("tap|pipe|open", MS_SPEC, ignore.case = TRUE) | is.na(TIME_WATER) | is.na(WATER_AVAIL) ~ NA_real_,
        TRUE ~ NA_real_ # Default case
      )
    )
  
  # GEMS -- > VIDA
  # GEMS X = VIDA Y
  # GEMS 1 = VIDA 1
  # 2 = 4
  # 3 = 5 or 6
  # 4 = 2 
  # 5 = 9
  # 6 = 10 (specify)
  # 7 = 4
  # HANGING, HANGING LATRINE, HANGING TOILET = 8
  
  # note this is closest to VIDA, don't have all the same options
  data_full <- data_full %>%
    mutate(
      MAIN_WASTE = coalesce(F4A_FAC_WASTE, F7_FAC_WASTE),
      MAIN_SPEC = coalesce(F4A_FAC_SPEC, F7_FAC_SPEC),
      SHARE_FAC = coalesce(F4A_SHARE_FAC, F7_SHARE_FAC),
      safe_sanit = case_when(
        MAIN_WASTE %in% c(1, 4, 7, 3) & SHARE_FAC == 0 ~ 1,
        MAIN_WASTE %in% c(1, 4, 7, 3) & SHARE_FAC >= 1 ~ 2,
        MAIN_WASTE %in% c(3) | SHARE_FAC %in% c("HANGING", "HANGING LATRINE", "HANGING TOILET") ~ 3,
        MAIN_WASTE %in% c(5, 6) ~ 4,
        is.na(MAIN_WASTE) | is.na(SHARE_FAC) ~ NA_real_,
        TRUE ~ NA_real_ 
      )
    )
  
  data_full$safe_water <- factor(data_full$safe_water, levels = 1:5, labels = c("Safely managed",
                                                                                "Basic",
                                                                                "Limited",
                                                                                "Unimproved",
                                                                                "Surface water"))
  
  data_full$safe_sanit <- factor(data_full$safe_sanit, levels = 1:4, labels = c("Safely managed and basic",
                                                                                "Limited",
                                                                                "Unimproved",
                                                                                "Open defication"))
  
  # Severity variables
  # Cases only -- ignore comparable variables for controls
  
  
  # - dysentery (F4B_OUTCOME_DYS - left hospital with dysentery)
  #             (F4A_DRH_STOOLS = 4 - describe diarrhea as bloody)
  #             (F4A_DRH_BLOOD = 1 - blood in stool since illness began, 9 NA)
  
  data_full$dysentery <- ifelse(data_full$F4B_OUTCOME_DYS == 1 | 
                                  data_full$F4A_DRH_STOOLS == 4 | 
                                  data_full$F4A_DRH_BLOOD == 1, 1, 0)
  
  # control -- has had blood in stool within last 7 days
  data_full$dysentery <- ifelse(!is.na(data_full$F7_BLOOD), data_full$F7_BLOOD, data_full$dysentery)
  
  # NEW exclude controls with dysentery
  data_full <- data_full[-which(data_full$case == 0 & data_full$dysentery == 1),]
  
  # - vomitting more than 3x per day (F4A_DRH_VOMIT = 1 - vomit 3 or more times / day since illness began)
  # for controls vomit more than 3x per day in last seven days
  data_full$vomit <- ifelse(data_full$case == 1,
                            data_full$F4A_DRH_VOMIT,
                            data_full$F7_VOMIT)
  
  # - fever (F4A_DRH_FEVER = 1 - fever 38C or higher since illness began )
  # for controls at least 38 C in last seven days
  data_full$fever <- ifelse(data_full$case == 1,
                            data_full$F4A_DRH_FEVER,
                            data_full$F7_FEVER)
  
  data_full$fever <- ifelse(data_full$fever == 9, NA, data_full$fever)
  
  # based on vida definition
  data_full$SEVERE_DEHYD <- ifelse(data_full$F4B_MENTAL == 2, 1, 0)
  data_full$SEVERE_DEHYD <- ifelse(data_full$F4B_EYES == 1, data_full$SEVERE_DEHYD + 1, data_full$SEVERE_DEHYD)
  data_full$SEVERE_DEHYD <- ifelse(data_full$F4A_DRH_LESSDRINK == 1 | data_full$F4A_DRH_UNDRINK == 1, data_full$SEVERE_DEHYD + 1, data_full$SEVERE_DEHYD)
  data_full$SEVERE_DEHYD <- ifelse(data_full$F4B_SKIN == 2, data_full$SEVERE_DEHYD + 1, data_full$SEVERE_DEHYD)
  
  data_full$MOD_DEHYD <- ifelse(data_full$F4B_MENTAL == 1, 1, 0)
  data_full$MOD_DEHYD <- ifelse(data_full$F4B_EYES == 1, data_full$MOD_DEHYD + 1, data_full$MOD_DEHYD)
  data_full$MOD_DEHYD <- ifelse(data_full$F4A_DRH_THIRST == 1, data_full$MOD_DEHYD + 1, data_full$MOD_DEHYD)
  data_full$MOD_DEHYD <- ifelse(data_full$F4B_SKIN == 1, data_full$MOD_DEHYD + 1, data_full$MOD_DEHYD)
  
  data_full$who_dehyd <- ifelse(data_full$SEVERE_DEHYD >= 2, 2, 0)
  data_full$who_dehyd <- ifelse(data_full$F4B_EYES == 0 & 
                                  data_full$SEVERE_DEHYD == 0 & 
                                  data_full$MOD_DEHYD >= 2, 1, data_full$who_dehyd)
  data_full$who_dehyd <- ifelse(data_full$F4B_EYES == 0 & 
                                  data_full$SEVERE_DEHYD == 1 & 
                                  data_full$MOD_DEHYD >= 1, 1, data_full$who_dehyd)
  data_full$who_dehyd <- ifelse(data_full$F4B_EYES == 1 & 
                                  data_full$SEVERE_DEHYD == 1 & 
                                  data_full$MOD_DEHYD >= 2, 1, data_full$who_dehyd)
  
  data_full$who_dehyd <- factor(data_full$who_dehyd, levels = 0:2, labels = c("No dehydration",
                                                                              "Some dehydration",
                                                                              "Severe dehydration"))
  
  data_full$lsstools <- factor(data_full$F4A_MAX_STOOLS, levels = 1:3, labels = c("6 or less per day",
                                                                                  "7 to 10 per day",
                                                                                  "Over 10 per day"))
  
  # Duration before enrollment
  data_full$duration_pre_enroll <- data_full$F4A_DRH_DAYS
  
  # Duration post-enrollment -- from diarrhea dictionary
  # added 9/4/26
  # check with liz
  # only include people who fully completed? memory aid -- 1
  # data has 0 through 14, pic in manuscript has 1 - 14, so verify dates w pre-enroll 
  # pre-enroll includes enrollment date
  data_full$duration_post_enroll <- ifelse(!is.na(data_full$F9_DRH_LAST) & data_full$F9_MEMORY_AID == 1, data_full$F9_DRH_LAST, NA)
  
  # Death -- in facility, at 60 day visit, overall
  # data_full$death = overall (in facility + 60 day visit)
  # data_full$died60 = 60 day visit
  # data_full$diedhosp = died in hospital before discharge
  
  # visit health facility during follow up
  # not quite same as rehosp 
  # 7899 total == 1
  data_full$visithf60 <- ifelse(data_full$F5_STATUS == 1,
                                as.integer(
                                  rowSums(
                                    data_full[, c(
                                      "F5_EXP_DRH_VISIT",
                                      "F5_EXP_COU_VISIT",
                                      "F5_EXP_DYS_VISIT",
                                      "F5_EXP_FEVER_VISIT",
                                      "F5_EXP_OTHR_VISIT",
                                      "F5_EXP_OTHR2_VISIT"
                                    )] == 1,
                                    na.rm = TRUE
                                  ) > 0 
                                ), NA)
  
  data_select <- data_full %>%
    select(child_id,
           case_id,
           first_id,
           case,
           sex,
           age,
           enr_haz,
           enr_whz,
           enr_waz,
           SITE,
           prim_caregiver_edu,
           num_hh_lt5,
           breastfed,
           dysentery,
           vomit,
           fever,
           duration_pre_enroll,
           duration_post_enroll,
           who_dehyd,
           lsstools,
           safe_sanit,
           safe_water, 
           ses_quintile,
           shigella_culture_pos,
           shigella_attributable_tac, 
           shigella_attributable,
           Shig_flex, 
           Shig_sonnei,
           rotavirus_attributable,
           noro_attributable,
           adenovirus_attributable,
           sapovirus_attributable,
           astro_attributable,
           etec_attributable,
           tepec_attributable,
           cryptosporidium_attributable,
           v_cholerae_attributable,
           c_jejuni_coli_attributable,
           salmonella_attributable,
           st_etec_attributable,
           rota_detected,
           adeno_detected,
           etec_detected,
           crypto_detected,
           astro_detected,
           noro_detected,
           tepec_detected,
           campy_detected,
           sapo_detected,
           e_bieneusi_detected,
           giardia_detected,
           EAEC_detected,
           shigella_culture_pos,
           rotavirus_elisa,
           adenovirus_elisa, 
           crypto_elisa, 
           etec_pcr,
           astro_pcr,
           noro_pcr,
           tEPEC_culture,
           campy_culture,
           sapo_pcr,
           EAEC_culture,
           giardia_elisa,
           adenovirus_new,
           astro_new,
           crypto_new,
           noro_new,
           rotavirus_new,
           sapo_new,
           shigella_new,
           st_etec_new,
           lt_etec_new,
           etec_new,
           tepec_new,
           campy_new,
           e_bieneusi_new,
           giardia_new,
           EAEC_new,
           v_cholerae_new,
           salmonella_new,
           c_jejuni_coli_new,
           followup_days,
           I_followup_days,
           I_followup_days_x_followup_days,
           hazdiff,
           whzdiff,
           hazd60,
           whzd60,
           wazd60,
           death,
           died60, 
           diedhosp,
           visithf60,
           any_abx,
           who_abx,
           maybe_eff_abx,
           ineff_abx, 
           all_abx,
           azithro,
           cipro) %>%
    rename("site" = SITE,
           "shig_flex" = Shig_flex,
           "shig_sonnei" = Shig_sonnei)
  
  # if culture negative, shig_flex and shig_sonnei should be NA (we do not know types of shig tac positive who were cultuer neg)
  data_select <- data_select %>%
    mutate(shig_flex = if_else(shigella_culture_pos == 0, NA, shig_flex),
           shig_sonnei = if_else(shigella_culture_pos == 0, NA, shig_sonnei))
  
  # Transform factors as needed
  
  # Sites: (key from matching IDs in gems1_form03_clean)
  # 1 - The Gambia
  # 2 - Mali
  # 3 - Mozambique
  # 4 - Kenya
  # 5 - India
  # 6 - Bangladesh
  # 7 - Pakistan
  
  data_select$site <- factor(data_select$site, levels = 1:7, labels = c("The Gambia",
                                                                        "Mali",
                                                                        "Mozambique",
                                                                        "Kenya",
                                                                        "India",
                                                                        "Bangladesh",
                                                                        "Pakistan"))
  
  data_select$sex <- factor(data_select$sex, levels = c(0,1), labels = c("male", "female"))
  
  data_select$prim_caregiver_edu <- ifelse(data_select$prim_caregiver_edu == 7, NA, data_select$prim_caregiver_edu)
  data_select$prim_caregiver_edu <- factor(data_select$prim_caregiver_edu, levels = 1:6, labels = c("No formal schooling",
                                                                                                    "Less than primary",
                                                                                                    "Completed primary",
                                                                                                    "Completed secondary",
                                                                                                    "Post secondary",
                                                                                                    "Religious education only"))
  
  data_select$prim_caregiver_edu_bin <- ifelse(data_select$prim_caregiver_edu %in% c("No formal schooling",
                                                                                     "Less than primary",
                                                                                     "Completed primary",
                                                                                     "Religious education only"), 0, 1)
  
  data_select$breastfed <- factor(data_select$breastfed, levels = c(0,1,2), labels = c("no",
                                                                                       "partial",
                                                                                       "exclusive"))
  
  
  # LABEL EVERYTHING FOR DESCRIPTIVE TABLES FUNCTION
  data_select <- data_select %>%
    set_variable_labels(case = "Case",
                        child_id = "Child ID",
                        case_id = "Case ID",
                        first_id = "First ID",
                        sex = "Sex",
                        age = "Age (months)",
                        enr_haz = "HAZ at enrollment",
                        enr_whz = "WHZ at enrollment",
                        site = "Enrollment site",
                        ses_quintile = "SES quintile",
                        prim_caregiver_edu = "Education level of primary caregiver",
                        prim_caregiver_edu_bin = "Primary caregiver education secondary school or greater",
                        num_hh_lt5 = "Number of children in household <5 years old",
                        breastfed = "Is child currently breastfed?",
                        dysentery = "Cases: Dysentery during episode | Controls: NA (exclude controls with dysentery)",
                        vomit = "Cases: Vomitted >3x/day since illness began | Controls: Vomitted >3x/day in last seven days",
                        fever = "Cases: Fever >=38C since illness began | Controls: >=38C in last seven days",
                        who_dehyd = "Cases: WHO defined dehydration | Controls: NA",
                        lsstools = "Cases: Number of loose stools on worst day",
                        shigella_culture_pos = "Shigella attributable via culture",
                        shigella_attributable_tac = "Shigella attributable diarrhea (AFE>0.5)",
                        shigella_attributable = "Shigella attributable diarrhea (AFE>0.5 or culture)",
                        rotavirus_attributable = "Rotavirus attributable diarrhea (AFE>0.5)",
                        st_etec_attributable = "ST ETEC attributable diarrhea (AFE>0.5)",
                        cryptosporidium_attributable = "Cryptosporidium attributable diarrhea (AFE>0.5)",
                        adenovirus_attributable = "Adenovirus attributable diarrhea (AFE>0.5)",
                        rotavirus_elisa = "ELISA immunoassay rotavirus",
                        crypto_elisa = "ELISA immunoassay cryptosporidium",
                        adenovirus_elisa = "ELISA adenovirus 40/41",
                        rota_detected = "Rotavirus detected",
                        adeno_detected = "Adenovirus 40/41 detected",
                        etec_detected = "ETEC detected",
                        crypto_detected = "Cryptosporidium detected",
                        astro_detected = "Astrovirus detected",
                        noro_detected = "Norovirus detected",
                        tepec_detected = "tEPEC detected",
                        campy_detected = "Campylobacter detected",
                        sapo_detected = "Sapovirus detected",
                        e_bieneusi_detected = "E bieneusi detected",
                        giardia_detected = "Giardia detected",
                        EAEC_detected = "EAEC detected",
                        shigella_culture_pos = "Shigella culture positive",
                        rotavirus_elisa = "Rotavirus ELISA immunoassay positive",
                        adenovirus_elisa = "Adenovirus 40/41 ELISA immunoassay positive",
                        crypto_elisa = "Cryptosporidium ELISA immunoassay positive",
                        etec_pcr = "ETEC",
                        astro_pcr = "Astrovirus",
                        noro_pcr = "Norovirus",
                        tEPEC_culture = "tEPEC",
                        campy_culture = "Campylobacter",
                        sapo_pcr = "Sapovirus",
                        EAEC_culture = "EAEC",
                        giardia_elisa = "Giardia",
                        hazdiff = "Difference in HAZ enrollment - day 60",
                        hazd60 = "HAZ at day 60",
                        followup_days = "Days between enrollment and follow-up",
                        any_abx = "Received any antibiotics",
                        who_abx = "Received WHO approved antibiotics",
                        maybe_eff_abx = "Recieved maybe effective antibiotics",
                        ineff_abx = "Recieved ineffective or no antibioitcs",
                        all_abx = "Type of antibiotics received",
                        azithro = "Received azithromycin",
                        cipro = "Received ciprofloxacin",
                        safe_water = "JMP Improved Water" ,
                        safe_sanit = "JMP Improved Sanitation")
  
  return(data_select)
  
}

# Potential case vs control definitions:
# - All diarrhea cases (all_diar)
# - TAC or culture Shigella attributable diarrhea vs no diarrhea (tac_or_culture_shig_diar)
# - TAC only Shigella diarrhea (tac_shig_diar)
# - Culture Shigella attributable diarrhea vs no diarrhea (culture_shig_diar)

gems_case_control_tac_or_culture <- prep_gems_case_control(case_def = "tac_or_culture_shig_diar")
gems_case_control_tac <- prep_gems_case_control(case_def = "tac_shig_diar")
gems_case_control_culture <- prep_gems_case_control(case_def = "culture_shig_diar")
gems_case_control_all <- prep_gems_case_control(case_def = "all_diar")

saveRDS(gems_case_control_tac_or_culture, here::here("data/gems_data/gems_case_control_tac_or_culture_shig.Rds"))
saveRDS(gems_case_control_tac, here::here("data/gems_data/gems_case_control_tac.Rds"))
saveRDS(gems_case_control_culture, here::here("data/gems_data/gems_case_control_culture.Rds"))
saveRDS(gems_case_control_all, here::here("data/gems_data/gems_case_control_all_diar.Rds"))
