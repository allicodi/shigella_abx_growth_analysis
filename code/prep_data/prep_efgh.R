# ---------------------------------------------------------------
# Script to create cleaned EFGH datasets
# ---------------------------------------------------------------

here::i_am("code/prep_data/prep_efgh.R")

library(tidyverse)
library(labelled)
library(haven)
library(corrr)
library(factoextra)

#' Function to help with cleaning antibiotic related data 
#' (adapted from UW medical management script)
clean_abx <- function(){
  
  # Four potential antibiotic timings -- preenrollment, at facility, for home, light touch 
  
  # 1. Pre-enrollment
  preenroll_abx <- readRDS(here::here("data/efgh_data/raw_data/DCS_04_enrollment.Rds")) %>%
    select(pid, enroll_hs_beh_tr___9, enroll_hs_beh_tr_ant___1:enroll_hs_beh_tr_oth_3, enroll_hs_beh_tr_oth1, enroll_hs_beh_tr_oth2) %>%
    rename(amox_admin = enroll_hs_beh_tr_ant___1,
           az_admin = enroll_hs_beh_tr_ant___2, 
           amp_admin = enroll_hs_beh_tr_ant___4, 
           aug_admin = enroll_hs_beh_tr_ant___5, 
           ceft_admin = enroll_hs_beh_tr_ant___6,  # ceftriaxone
           cip_admin = enroll_hs_beh_tr_ant___7,  # ciprofloxacin
           cef_admin = enroll_hs_beh_tr_ant___8,  # Cefuroxime
           cefix_admin = enroll_hs_beh_tr_ant___24, # cefixime
           clar_admin = enroll_hs_beh_tr_ant___9, # Clarithromycin
           clin_admin = enroll_hs_beh_tr_ant___10,  # Clindamycin
           chlo_admin = enroll_hs_beh_tr_ant___11,  # Chloramphenicol
           cot_admin = enroll_hs_beh_tr_ant___12,  # Cotrimoxazole/Septrin/trimethoprim-sulfamethoxazole
           dox_admin = enroll_hs_beh_tr_ant___13, 
           ery_admin = enroll_hs_beh_tr_ant___14, 
           gen_admin = enroll_hs_beh_tr_ant___15, 
           met_admin = enroll_hs_beh_tr___9,
           pen_admin = enroll_hs_beh_tr_ant___16, 
           pyr_admin = enroll_hs_beh_tr_ant___17, 
           str_admin = enroll_hs_beh_tr_ant___18, 
           tet_admin = enroll_hs_beh_tr_ant___19, 
           unk_admin = enroll_hs_beh_tr_ant___20,
           oth1_admin = enroll_hs_beh_tr_ant___21, 
           oth2_admin = enroll_hs_beh_tr_ant___22, 
           oth3_admin = enroll_hs_beh_tr_ant___23, 
           oth1_name = enroll_hs_beh_tr_oth_1,  
           oth2_name = enroll_hs_beh_tr_oth_2,  
           oth3_name = enroll_hs_beh_tr_oth_3,
           oth4_name = enroll_hs_beh_tr_oth1,
           oth5_name = enroll_hs_beh_tr_oth2) %>% 
    rowwise() %>%
    mutate(unique_abx_given = sum(across(ends_with("_admin")), na.rm = T)) %>%
    filter(unique_abx_given>=1 | !is.na(oth4_name) | !is.na(oth5_name)) %>%
    select(-unique_abx_given) %>%
    ungroup() %>%
    mutate_all(as.character) %>%
    # <- reshape data to be long rather than wide
    pivot_longer(cols = -pid, names_to = "abx", values_to = "value") %>%
    separate(abx, into = c("antibiotic", "variable"), sep = "_", extra = "merge") %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    # filter to only rows with administered antibiotics
    filter(admin==1 | !is.na(name)) %>%
    mutate(drug = antibiotic,
           drug_type = "antibiotic",
           treatment = "Received pre-enrollment",
           drug_route = NA,
           dose = NA,
           frequency = NA,
           indication = NA,
           dt = NA) %>%
    mutate(drug = case_when(drug=="amox" ~ "amoxicillin",
                            drug=="amp" ~ "ampicillin",
                            drug=="aug" ~ "augmentin",
                            drug=="az" ~ "azithromycin",
                            drug=="cef" ~ "cefuroxime",
                            drug=="cefix" ~ "cefixime",
                            drug=="ceft" ~ "ceftriaxone",
                            drug=="chlor" ~ "chloramphenicol",
                            drug=="cip" ~ "ciprofloxacin", 
                            drug=="clar" ~ "clarithromycin",
                            drug=="cot" ~ "cotrim",
                            drug=="dox" ~ "doxycycline",
                            drug=="ery" ~ "erythromycin",
                            drug=="gen" ~ "gentamicin",
                            drug=="met" ~ "metronidazole",
                            drug=="pen" ~ "penicillin",
                            drug=="tet" ~ "tetracycline",
                            drug=="unk" ~ "unknown",
                            TRUE ~ drug
    )) %>%
    select(pid, treatment, drug, drug_type, name, drug_route, dose, frequency, indication, dt)
  
  
  ############### CLEAN WRITE-INS ############### 
  
  preenroll_abx <- preenroll_abx %>%
    filter(!(name %in% c("paracetamol","Paracetamol","PARACETAMOL","Parcentamol","Pracentamol",
                         "ANALGESICS(PARACETAMOL)", "calpol", "Calpol", "artemether lumefantrine", 
                         "Carfenol", "CLORURO DE SODIO 9%", "Cough syrup", "Décoction",
                         "BISMUTOL", "cafenol", "Cafenol", "Cafenol tablet", "CLORFENAMINA", "CORTIMET 9 MG",
                         "Diclofenac", "Diclofenac IM", "DIMEHINDRINATO", "DIMENHIDRATO", "dolar",
                         "DONPERERIDONE", "FEXOFENADINE", "gravinate", "GRAVOL", "HERBAL", "HERBAL MEDICATION",
                         "HERBAL MEDICINE", "HIOSCINA", "hydrosac ", "inj ringer 250ml", "Inj ringer 600ml.",
                         "DIMEHIDRINATO", "DIMENHIDRINATO", "ENTEROGERMINA", "IVERMECTINA",
                         "iv rehydration ", "KETOTIFEN", "KETOTIFEN ", "MON BÉBÉ ", "Montelukast sodium", "multivitamin syrup",
                         "Multivitamin syrup", "MWARUBAINI", "nitazoxanide", "Nitazoxanide", "ONDENSTERON ",
                         "METAMIZOL", "ONDANSETRON", "ONDANSETRON ", "Pnanodol", "panadol", "RANITIDINA",
                         "Paracétamol", "Paracetamol ", "Paracetamol and Promethazine syrup", "Paracetamol syrup", 
                         "Paracetamol tablet", "Paracetamol tablets", "PCM SYRUP", "PERACITAMOL", "piriton", 
                         "powder syrup ", "probiotic ", "Pumpkin leaves", "REHIDRATACION INTRAVENOSA", "REHIDRATACION ENDOVENOSA",
                         "REHIDRATACION EV.", "syphilis chlorpheniramine", "Vitamin c syrup", "ONDARSATORN", "TIEMONIUM METHYLSULPHAT",
                         "pcm", "syp paracetamole", "SOLUCION SALINA", "solución salina normal", "paracetamol syrup",
                         "syb paracetamol,syp domel", "syp paracetemol", "syphilis calpol", "Tetracycline eye ointment",
                         "Tisane ", "Tisane", "TISANE", "TISANE ", "traditional concoction", "NITAZOXANIDE", 
                         "SOLUCIÓN SALINA", "Smecta ","homeopathic medicine ","multivitamin","promethazine syrup",
                         "Vitamin B complex"))) %>%
    mutate(drug = ifelse(drug=="oth1" & name=="AMIKACINA", "amikacin", drug),
           drug = ifelse(drug=="oth4" & name=="Amoxil", "amoxicillin", drug),
           drug = ifelse(drug=="oth1" & name=="AMPICLOX", "ampicillin/cloxacillin", drug),
           drug = ifelse(drug=="oth4" & name=="syp Loratadine,syp Cefaclor", "cefaclor", drug),
           drug = ifelse(drug %in% c("oth1", "oth4") & name %in% c("CEFIXIME", "syp.cefixime", "syp cefixime"), "cefixime", drug),
           drug = ifelse(drug=="oth1" & name=="Ceftizidine i/v inj.", "ceftazidime", drug),
           drug = ifelse(drug=="oth4" & name=="inj ceftrixone 250mg", "ceftriaxone", drug),
           drug = ifelse(drug %in% c("oth1", "oth4") & name %in% c("CIPROFLOXACINO", "ciprofloxacin", "CIPROFLOXACIN", "ciprofloxacin "), "ciprofloxacin", drug),
           drug = ifelse(drug %in% c("oth1", "oth4") & name %in% c("BACTRIL", "BACTRIM", "SULFAMETAZOL", "SULFAMETAZOL + TRIMETROPINA", "SULFAMETOXAZOL",
                                                                   "SULFAMETOXAZOL.", "SULFAMETROXAZOL", "sulphamethoxazole/Trimethoprim", "ROPSIL", "ROSTRIN",
                                                                   "SULFAME", "BACTEROL", "Cot rime", "syp sulphamethoxazole", "syp septran DS"), "cotrim", drug),
           drug = ifelse(drug=="oth1" & name=="Doxycline", "doxycycline", drug),
           drug = ifelse(drug %in% c("oth1", "oth2", "oth3", "oth4") & name %in% c("FURAZOLIDONA", "FUROZOLIDONA", "FURAZOLIDONA"), "furazolidone", drug),
           drug = ifelse(drug %in% c("oth1", "oth4") & name %in% c("metronidazole ", "Metronidazole ", "METRONIDAZOLE ", "Entamozole"), "metronidazole", drug),
           drug = ifelse(drug %in% c("oth1") & name=="DIAREN", "nifuroxazide", drug)) %>%
    select(-name)
  
  # 2. At facility
  facility_abx <- readRDS(here::here("data/efgh_data/raw_data/DCS_04a_medical_management.Rds"))  %>% 
    select(pid, med_ant_choice___1:med_ant_oth3_ind_oth) %>% 
    rename(med_ant_amox_admin = med_ant_choice___1,
           med_ant_az_admin = med_ant_choice___2, 
           med_ant_amp_admin = med_ant_choice___3, 
           med_ant_aug_admin = med_ant_choice___4, 
           med_ant_cefix_admin = med_ant_choice___5, 
           med_ant_ceft_admin = med_ant_choice___6, 
           med_ant_cip_admin = med_ant_choice___7, 
           med_ant_cef_admin = med_ant_choice___8, 
           med_ant_clin_admin = med_ant_choice___9, 
           med_ant_chlo_admin = med_ant_choice___10, 
           med_ant_cot_admin = med_ant_choice___11, 
           med_ant_dox_admin = med_ant_choice___12, 
           med_ant_ery_admin = med_ant_choice___13, 
           med_ant_gen_admin = med_ant_choice___14, 
           med_ant_met_admin = med_ant_choice___15, 
           med_ant_pen_admin = med_ant_choice___16, 
           med_ant_pyr_admin = med_ant_choice___17, 
           med_ant_str_admin = med_ant_choice___18, 
           med_ant_tet_admin = med_ant_choice___19, 
           med_ant_oth1_admin = med_ant_choice___20, 
           med_ant_oth2_admin = med_ant_choice___21, 
           med_ant_oth3_admin = med_ant_choice___22, 
           med_ant_clar_admin = med_ant_choice___23, 
           med_ant_oth1_name = med_ant_oth1,  
           med_ant_oth2_name = med_ant_oth2,  
           med_ant_oth3_name = med_ant_oth3, 
           med_ant_cip_complete = med_ant_cipro_complete) %>% 
    rowwise() %>%
    mutate(unique_abx_given = sum(across(ends_with("_admin")))) %>%
    filter(unique_abx_given>=1) %>%
    select(-unique_abx_given) %>%
    ungroup() %>%
    filter(pid != 4108711) %>% # JUST UNTIL MALI RESPONDS
    mutate_all(as.character) %>%
    # <- reshape data to be long rather than wide
    pivot_longer(cols = -pid, names_to = "abx", values_to = "value") %>%
    separate(abx, into = c("prefix1", "prefix2", "antibiotic", "variable"), sep = "_", extra = "merge") %>%
    select(-prefix1,-prefix2) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    # filter to only rows with administered antibiotics
    filter(admin==1) %>%
    mutate(drug_route = case_when(rt ==1 ~ "IV", 
                                  rt== 2 ~ "oral", 
                                  rt==3 ~ "IM", 
                                  rt==4 ~ "other"),
           dose = as.numeric(dos),
           dose_known = as.numeric(dosknown),
           treatment = "Administered",
           frequency = fq,
           indication = ind,
           drug = antibiotic,
           drug_type = "antibiotic") %>%
    select(pid, treatment,drug, drug_type, drug_route, rt_oth, dose_known, dose, frequency, fq_oth, 
           indication, ind_oth, complete, name, dt, dt_com)
  
  ############### CLEAN WRITE-INS ############### 
  
  facility_abx <- facility_abx %>%
    filter(!(name %in% c("Syp.calpol"))) %>%
    mutate(drug = ifelse(drug=="oth1" & name==" MEROPENAM","meropenam",drug),
           drug = ifelse(drug=="oth1" & name=="AMIKACINA", "amikacin", drug),
           drug = ifelse(drug=="oth1" & name=="BENZYLPENICILLIN", "benzylpenicillin", drug),
           drug = ifelse(drug=="oth1" & (name=="cloxacillin" | name=="Cloxacillin"),"cloxacillin", drug),
           drug = ifelse(drug=="oth3" & name=="floxapen ", "flucloxacillin", drug),
           drug = ifelse(drug=="oth1" & name=="FLUCLOXACILLINE", "flucloxacillin", drug),
           drug = ifelse(drug=="oth1" & (name=="Gentamycen eye drops" | name=="Gentamycin Ear Drops"), "gentamicin",drug),
           drug = ifelse(drug=="oth2" & name=="GENTAMYCIN", "gentamicin", drug),
           drug = ifelse(drug=="oth1" & name=="LEVOFLOXACIN", "levofloxacin", drug),
           drug = ifelse(drug=="oth1" & (name=="Metronidazole" | name=="Metronidazole " |
                                           name=="METRONIDAZOLE " | name=="Métronidazole "), "metronidazole", drug),
           drug = ifelse(drug=="oth1" & name=="PENICILINA PROCAINICA", "procain penicillin", drug)) %>%
    mutate(drug_route = ifelse(drug_route=="other", rt_oth, drug_route)) %>%
    select(-rt_oth,-name)
  
  # 3. For home
  
  home_abx <- readRDS(here::here("data/efgh_data/raw_data/DCS_04a_medical_management.Rds"))  %>% 
    select(pid, med_ant_choice_hm___1:med_ant_oth3_ind_hm_oth) %>% 
    rename(med_ant_amox_presc = med_ant_choice_hm___1,
           med_ant_az_presc = med_ant_choice_hm___2, 
           med_ant_amp_presc = med_ant_choice_hm___3, 
           med_ant_aug_presc = med_ant_choice_hm___4, 
           med_ant_cefix_presc = med_ant_choice_hm___5, 
           med_ant_ceft_presc = med_ant_choice_hm___6, 
           med_ant_cip_presc = med_ant_choice_hm___7, 
           med_ant_cef_presc = med_ant_choice_hm___8, 
           med_ant_clin_presc = med_ant_choice_hm___9, 
           med_ant_chlor_presc = med_ant_choice_hm___10, 
           med_ant_cot_presc = med_ant_choice_hm___11, 
           med_ant_dox_presc = med_ant_choice_hm___12, 
           med_ant_ery_presc = med_ant_choice_hm___13, 
           med_ant_gen_presc = med_ant_choice_hm___14, 
           med_ant_pen_presc = med_ant_choice_hm___15, 
           med_ant_pyr_presc = med_ant_choice_hm___16, 
           med_ant_str_presc = med_ant_choice_hm___17, 
           med_ant_tet_presc = med_ant_choice_hm___18, 
           med_ant_oth1_presc = med_ant_choice_hm___19, 
           med_ant_oth2_presc = med_ant_choice_hm___20, 
           med_ant_oth3_presc = med_ant_choice_hm___21, 
           med_ant_met_presc = med_ant_choice_hm___22,
           med_ant_clar_presc = med_ant_choice_hm___23,
           med_ant_oth1_name = med_ant_oth1_hm,  
           med_ant_oth2_name = med_ant_oth2_hm,  
           med_ant_oth3_name = med_ant_oth_hm) %>% 
    rowwise() %>%
    mutate(unique_abx_given = sum(across(ends_with("_presc")))) %>%
    filter(unique_abx_given>=1) %>%
    filter(pid != 4108711) %>% # JUST UNTIL MALI RESPONDS
    select(-unique_abx_given) %>%
    ungroup() %>%
    mutate(across(everything(), as.character)) %>%
    # <- reshape data to be long rather than wide
    pivot_longer(cols = -pid, names_to = "abx", values_to = "value") %>%
    separate(abx, into = c("prefix1", "prefix2", "antibiotic", "variable"), sep = "_", extra = "merge") %>%
    select(-prefix1,-prefix2) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    # filter to only rows with administered antibiotics
    filter(presc==1) %>%
    mutate(dose = as.numeric(dos_hm),
           dose_known = as.numeric(hm_dosknown),
           treatment = "Prescribed",
           frequency = fq_hm,
           indication = ind_hm,
           drug = antibiotic,
           drug_type = "antibiotic",
           duration = pres_dur_hm) %>%
    select(pid, treatment, drug, drug_type, dose_known, dose, frequency, fq_oth_hm, duration, 
           indication, ind_hm_oth, name)
  
  ############### CLEAN WRITE-INS ############### 
  
  home_abx <- home_abx %>%
    mutate(drug=ifelse(drug %in% c("oth1", "oth2", "oth3"), "oth", drug)) %>%
    mutate(drug = ifelse(drug=="oth" & name=="Amikacin", "amikacin", drug),
           drug = ifelse(drug=="oth" & (name=="ampiclox" | name=="Ampiclox" |
                                          name=="AMPICLOX" | name=="Ampiclox syrup" |
                                          name=="Ampiclox syrup (ampicillin-cloxacillin)" |
                                          name=="AMPICLOX SYRUP"), "ampiclox", drug),
           drug = ifelse(drug=="oth" & (name=="Biodroxil" | name=="BIODROXIL" |
                                          name=="Biodroxil " | name=="Fedrox"), "cefadroxil", drug),
           drug = ifelse(drug=="oth" & name=="Cefpodoxime", "cefpodoxime", drug),
           drug = ifelse(drug=="oth" & (name=="cloxacillin" | name=="Cloxacillin" |
                                          name=="cloxacillin " | name=="Cloxacillin " |
                                          name=="cloxacillin syrup" | name=="Cloxacillin syrup" |
                                          name=="CLOXACILLIN SYRUP"), "cloxacillin", drug),
           drug = ifelse(drug=="oth" & (name=="Entamizole" | name=="ENTAMIZOLE"), "entamizole", drug),
           drug = ifelse(drug=="oth" & name=="ERYTROMICINE", "erythromycin", drug),
           drug = ifelse(drug=="oth" & (name=="Flagyl" | name=="Flagyl " |
                                          name=="METRANIDAZOLE" | name=="METRANIDAZOLE " |
                                          name=="METRANIDAZOLLE" | name=="METRODINAZOL" |
                                          name=="METRONIDAZOL" | name=="Metronidazol " |
                                          name=="metronidazole" | name=="Metronidazole" |
                                          name=="METRONIDAZOLE" | name=="Metronidazole " |
                                          name=="METRONIDAZOLE " | name=="Métronidazole " |
                                          name=="Metronidazole  " | name=="METRONIDAZOLE  " |
                                          name=="Metrononidazole"), "metronidazole", drug),
           drug = ifelse(drug=="oth" & (name=="floxapen" | name=="Floxapen" |
                                          name=="FLOXAPEN" | name=="floxapen syrup"), "floxapen", drug),
           drug = ifelse(drug=="oth" & (name=="FLUCLOXACACILLIN" | name=="flucloxacillin" |
                                          name=="Flucloxacillin" | name=="FLUCLOXACILLIN" |
                                          name=="FLUCLOXACILLIN SYRUP" | name=="FLUCLOXACILLINE" |
                                          name=="SYRUP FLOXAPEN"),"flucloxacillin", drug),
           drug = ifelse(drug=="oth" & (name=="Furazolidona" | name=="FURAZOLIDONA" |
                                          name=="FUROZOLIDONA"), "furazolidone", drug),
           drug = ifelse(drug=="oth" & (name=="Gentamicin ear drop" | name=="Gentamicin eye drop" |
                                          name=="Gentamycin ear drops" | name=="Gentamycin eye drop" |
                                          name=="Gentamycin eye drops"), "gentamicin", drug),
           drug = ifelse(drug=="oth" & name=="MUPIROCINA", "mupirocin", drug),
           drug = ifelse(drug=="oth" & name=="NITROFURANTOINA", "nitrofurantoin", drug),
           drug = ifelse(drug=="oth" & name=="Orelox", "cefpodoxime", drug),
           drug = ifelse(drug=="oth" & (name=="PIVMECILLINAM" | name=="PIVMECILLINAM "),
                         "pivmecillinam", drug),
           drug = ifelse(drug=="oth" & name=="SYRUP AMPICLOX", "ampiclox", drug),
           drug = ifelse(drug=="oth" & name=="TETRACYCLINE EYE OINTMENT", "tetracycline", drug),
           drug = ifelse(drug=="oth" & (name=="Thiobactin " | name=="Thiobactin"), "thiamphenicol", drug)) %>%
    filter(!(name %in% c("Triple action cream","XTRADERM CREAM","Teething syrup","syrup piriton",
                         "SYRUP FLUGONE DM","SRUP PIRITON","Smecta ","probiotic","Rigix",
                         "PREDNISOLONE","Ondensetron","Nystatn oral","Kalamine lotion","AMBROXOL",
                         "CLOTRIMAZOLE CREAM","Flucazol","MULTIVITAMIN", "tabs Albendazole"))) %>%
    select(-name)
  
  # 4. Light touch 
  unwell <- readRDS(here::here("data/efgh_data/raw_data/DCS_08_unwell_child.Rds"))
  
  light <- unwell %>%
    select(pid,uv_purpose,uv_treat_pres___1:uv_oth2_pres) %>%
    filter(uv_purpose==3) %>%
    select(-uv_purpose) %>%
    rename(uv_cip_light = uv_treat_pres___1,
           uv_az_light = uv_treat_pres___2,
           uv_cot_light = uv_treat_pres___3,
           uv_ceft_light = uv_treat_pres___4,
           uv_cefix_light = uv_treat_pres___5,
           uv_oth1_light = uv_treat_pres___6,
           uv_oth2_light = uv_treat_pres___7) %>% 
    rowwise() %>%
    mutate(unique_abx_given = sum(across(ends_with("_light")))) %>%
    filter(unique_abx_given>=1) %>%
    select(-unique_abx_given) %>%
    ungroup() %>%
    mutate(across(everything(), as.character)) %>%
    # <- reshape data to be long rather than wide
    pivot_longer(cols = -pid, names_to = "abx", values_to = "value") %>%
    separate(abx, into = c("prefix1", "antibiotic", "variable"), sep = "_", extra = "merge") %>%
    select(-prefix1) %>%
    #############
  # TEMPORARY FILTER - SHOULD BE FIXED 7/11/2024
  #############
  # filter(pid!=7110154) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
    filter(light==1) %>%
    mutate(treatment = "Light Contact",
           drug = antibiotic,
           drug_type = "antibiotic",
           unit = ifelse(unit==1, "ml", 
                         ifelse(unit==2, "mg", unit)),
           # dose = paste0(dosage,unit),
           dose = as.numeric(dosage),
           indication = "shigella") %>%
    select(pid, treatment, drug, drug_type, pres, dose, spec, indication)
  
  light <- light %>%
    mutate(drug = ifelse(drug=="oth1" & spec=="meropenam", "meropenam", drug),
           drug = ifelse(drug=="oth1" & spec=="METRONIDAZOL", "metronidazole", drug),
           drug = ifelse(drug=="oth1" & (spec=="PIVMECILLINAM" | spec=="PIVMECILLINAM "),
                         "pivmecillinam", drug),
           drug = ifelse(drug=="oth1" & spec=="Trimethoprim/Sulfamethoxazole", "cot", drug)) %>%
    filter(!(spec %in% c("cough syp","Vogalène"))) %>%
    select(-spec, -pres) %>%
    mutate(drug = case_when(drug=="az" ~ "azithromycin",
                            drug=="ceft" ~ "ceftriaxone",
                            drug=="cip" ~ "ciprofloxacin",
                            drug=="cot" ~ "cotrim",
                            TRUE ~ drug))
  
  # FINAL - Combine datasets
  
  all_abx <- bind_rows(facility_abx,home_abx,light) %>%
    mutate(indication = case_when(indication==1 ~ "diarrhea",
                                  indication==2 ~ "dysentery",
                                  indication==3 ~ "sepsis",
                                  indication==4 ~ "lower respiratory infection",
                                  indication==5 ~ "malnutrition",
                                  indication==6 ~ "other",
                                  indication=="shigella" ~ "shigella")) %>%
    mutate(frequency = case_when(frequency==1 ~ "1x/day",
                                 frequency==2 ~ "2x/day",
                                 frequency==3 ~ "3x/day",
                                 frequency==4 ~ "4x/day",
                                 frequency==5 ~ "other")) %>%
    mutate(indication = ifelse(indication=="other" & !is.na(ind_oth), ind_oth, indication),
           indication = ifelse(indication=="other" & !is.na(ind_hm_oth), ind_hm_oth, indication),
           frequency = ifelse(frequency=="other" & !is.na(fq_oth), fq_oth, frequency),
           frequency = ifelse(frequency=="other" & !is.na(fq_oth_hm), fq_oth_hm, frequency)) %>%
    select(-fq_oth, -ind_hm_oth, -ind_oth, -fq_oth_hm, -complete, -dt_com, -duration, -dose_known)  %>%
    mutate(drug = case_when(drug=="amox" ~ "amoxicillin",
                            drug=="amp" ~ "ampicillin",
                            drug=="aug" ~ "augmentin",
                            drug=="az" ~ "azithromycin",
                            drug=="cef" ~ "cefuroxime",
                            drug=="cefix" ~ "cefixime",
                            drug=="ceft" ~ "ceftriaxone",
                            drug=="chlor" ~ "chloramphenicol",
                            drug=="cip" ~ "ciprofloxacin", 
                            drug=="clar" ~ "clarithromycin",
                            drug=="cot" ~ "cotrim",
                            drug=="ery" ~ "erythromycin",
                            drug=="gen" ~ "gentamicin",
                            drug=="met" ~ "metronidazole",
                            drug=="pen" ~ "penicillin",
                            drug=="tet" ~ "tetracycline",
                            TRUE ~ drug
    )) %>%
    mutate(dose = as.character(dose))
  
  # note- pre-enrollment does not have indication info
  all_abx <- all_abx %>%
    bind_rows(preenroll_abx)
  
  # --------------------------------------------------------------------------
  # Adding susceptibility data (5/28/25)
  # --------------------------------------------------------------------------
  
  ast_data <- readRDS(here::here("data/efgh_data/raw_data/DCS_16_rectal_swab_results.Rds")) %>%
    select(pid,
           sw_date,
           swi_ast_amp_int,
           swi_ast_az_int,
           swi_ast_ceft_int,
           swi_ast_cipro_int,
           swi_ast_nal_int,
           swi_ast_piv_int,
           swi_ast_tri_int) %>%
    mutate(pid = as.character(pid)) %>%
    mutate(swi_ast_az_int = if_else(swi_ast_az_int == 1, "S",
                                    if_else(swi_ast_az_int == 2 | swi_ast_az_int == 3, "I/R", NA)),
           swi_ast_cipro_int = if_else(swi_ast_cipro_int == 1, "S",
                                       if_else(swi_ast_cipro_int == 2 | swi_ast_cipro_int == 3, "I/R", NA)),
           swi_ast_ceft_int = if_else(swi_ast_ceft_int == 1, "S",
                                      if_else(swi_ast_ceft_int == 2 | swi_ast_ceft_int == 3, "I/R", NA)),
           swi_ast_amp_int = if_else(swi_ast_amp_int == 1, "S",
                                     if_else(swi_ast_amp_int == 2 | swi_ast_amp_int == 3, "I/R", NA)),
           swi_ast_nal_int = if_else(swi_ast_nal_int == 1, "S",
                                     if_else(swi_ast_nal_int == 2 | swi_ast_nal_int == 3, "I/R", NA)),
           swi_ast_piv_int = if_else(swi_ast_piv_int == 1, "S",
                                     if_else(swi_ast_piv_int == 2 | swi_ast_piv_int == 3, "I/R", NA)),
           swi_ast_tri_int = if_else(swi_ast_tri_int == 1, "S",
                                     if_else(swi_ast_tri_int == 2 | swi_ast_tri_int == 3, "I/R", NA))) %>%
    group_by(pid) %>%
    #group by PID, if there are multiple entries, only make it "S" if all are S
    summarise(
      sw_date = sw_date[1],
      swi_ast_az_int = case_when(
        all(is.na(swi_ast_az_int)) ~ NA_character_,
        all(swi_ast_az_int == "S", na.rm = TRUE) ~ "S",
        any(swi_ast_az_int == "I/R", na.rm = TRUE) ~ "I/R",
        TRUE ~ NA_character_
      ),
      swi_ast_cipro_int = case_when(
        all(is.na(swi_ast_cipro_int)) ~ NA_character_,
        all(swi_ast_cipro_int == "S", na.rm = TRUE) ~ "S",
        any(swi_ast_cipro_int == "I/R", na.rm = TRUE) ~ "I/R",
        TRUE ~ NA_character_
      ),
      swi_ast_ceft_int = case_when(
        all(is.na(swi_ast_ceft_int)) ~ NA_character_,
        all(swi_ast_ceft_int == "S", na.rm = TRUE) ~ "S",
        any(swi_ast_ceft_int == "I/R", na.rm = TRUE) ~ "I/R",
        TRUE ~ NA_character_
      ),
      swi_ast_nal_int = case_when(
        all(is.na(swi_ast_nal_int)) ~ NA_character_,
        all(swi_ast_nal_int == "S", na.rm = TRUE) ~ "S",
        any(swi_ast_nal_int == "I/R", na.rm = TRUE) ~ "I/R",
        TRUE ~ NA_character_
      ),
      swi_ast_amp_int = case_when(
        all(is.na(swi_ast_amp_int)) ~ NA_character_,
        all(swi_ast_amp_int == "S", na.rm = TRUE) ~ "S",
        any(swi_ast_amp_int == "I/R", na.rm = TRUE) ~ "I/R",
        TRUE ~ NA_character_
      ),
      swi_ast_piv_int = case_when(
        all(is.na(swi_ast_piv_int)) ~ NA_character_,
        all(swi_ast_piv_int == "S", na.rm = TRUE) ~ "S",
        any(swi_ast_piv_int == "I/R", na.rm = TRUE) ~ "I/R",
        TRUE ~ NA_character_
      ),
      swi_ast_tri_int = case_when(
        all(is.na(swi_ast_tri_int)) ~ NA_character_,
        all(swi_ast_tri_int == "S", na.rm = TRUE) ~ "S",
        any(swi_ast_tri_int == "I/R", na.rm = TRUE) ~ "I/R",
        TRUE ~ NA_character_
      )
    )
  
  abx_ast_data <- left_join(all_abx, ast_data) %>%
    mutate(ast_given_abx = case_when(
      drug == "ciprofloxacin" ~ swi_ast_cipro_int,
      drug == "azithromycin" ~ swi_ast_az_int,
      drug == "ceftriaxone" ~ swi_ast_ceft_int,
      drug == "ampicillin" ~ swi_ast_amp_int,
      drug == "pivmecillinam" ~ swi_ast_piv_int,
      TRUE ~ NA_character_
    )) 
  
  # --------------------------------------------------------------------------
  
  final_dataset <- abx_ast_data %>%
    mutate(any_abx = 1, 
           who_rec_abx = ifelse(drug=="azithromycin" | drug=="ciprofloxacin" | 
                                  drug=="ceftriaxone" | drug=="pivmecillinam", 1, 0),
           # ADDED 3/3/25 FROM UPDATED SCRIPT
           maybe_eff_abx = ifelse(drug %in% c("amoxicillin","ampicillin","ampicillin/cloxacillin",
                                              "cefixime","cefuroxime","chloramphenicol",
                                              "cotrim","doxycycline","erythromycin",
                                              "tetracycline","clarithromycin","gentamicin",
                                              "ampiclox","augmentin","amikacin",
                                              "cefaclor","cefpodoxime","ceftazidime",
                                              "furazolidone","levofloxacin","meropenam",
                                              "nifuroxazide","nitrofurantoin",
                                              "thiamphenicol"), 1, 0),
           ineff_abx = ifelse(drug %in% c("metronidazole","penicillin",
                                          "procain penicillin","unknown",
                                          "benzylpenicillin","cefadroxil",
                                          "cloxacillin","entamizole","flucloxacillin",
                                          "mupirocin"), 1, 0),
           anti_diarrhea_abx = ifelse(indication %in% c("Acute Gasroentritis", "Acute Gasroentritis",
                                                        "acute gastroentritis", "Blood visible on rectal swabs.",
                                                        "diarrhea", "dysentery", "Dysentery and Urinary tract infection ",
                                                        "Enteric fever","Enteric Fever","ENTERIC FEVER",
                                                        "Enteric fever ","gastritis","GASTRITIS","GASTRO-ENTERITIS",
                                                        "GASTRO ENTERITIS","gastroenteritis","Gastroenteritis",
                                                        "GASTROENTERITIS","gastroenteritis ","Gastroenteritis ",
                                                        "Gastroenteritis/ Upper Respiratory Tract Infection",
                                                        "gastroentritis","Hyperthermia, diarrhea and DD ",
                                                        "suspected cholera", "Suspicion of enteric fever",
                                                        "upper respiratory tract infection  and diarrhea",
                                                        "UTI(Hematuria) we had no ciprofloxacin so started cefixime due to diarrhea also",
                                                        "Hyperthermia,  diarrhea ","shigella"), 1, 0),
           azithromycin = ifelse(drug == "azithromycin", 1, 0),
           ciprofloxacin = ifelse(drug == "ciprofloxacin", 1, 0),
           ceftriaxone = ifelse(drug == "ceftriaxone", 1, 0),
           pivmecillinam = ifelse(drug == "pivmecillinam", 1, 0))%>%
    # TODO- CHECK ON THIS ABT RECEIVED PRE-ENROLLMENT
    mutate(ineff_abx = ifelse(treatment =="Received pre-enrollment", 0, ineff_abx),
           maybe_eff_abx = ifelse(treatment == "Received pre-enrollment", 0, maybe_eff_abx),
           who_rec_abx = ifelse(treatment=="Received pre-enrollment", 0, who_rec_abx),
           anti_diarrhea_abx = ifelse(treatment=="Received pre-enrollment", 0, anti_diarrhea_abx),
           azithromycin = ifelse(treatment=="Received pre-enrollment", 0, azithromycin),
           ciprofloxacin = ifelse(treatment=="Received pre-enrollment", 0, ciprofloxacin),
           ceftriaxone = ifelse(treatment=="Received pre-enrollment", 0, ceftriaxone),
           pivmecillinam = ifelse(treatment=="Received pre-enrollment", 0, pivmecillinam)) %>%
    group_by(pid) %>%
    summarise(who_rec_abx = max(who_rec_abx),
              maybe_eff_abx = max(maybe_eff_abx),
              ineff_abx = max(ineff_abx),
              anti_diarrhea_abx = max(anti_diarrhea_abx),
              azithromycin = max(azithromycin),
              ciprofloxacin = max(ciprofloxacin),
              ceftriaxone = max(ceftriaxone),
              pivmecillinam = max(pivmecillinam),
              # if given multiple abx and any are susceptible, mark as susceptible
              ast_given_abx = ifelse(any(ast_given_abx == "S"), "S", ast_given_abx))#,
  # resistant_WHO_approve = max(resistant_WHO_approve),
  # susceptible_WHO_approve = max(susceptible_WHO_approve))
  
  # WHO, Maybe, and Ineff should be mutually exclusive - if received multiple go with highest
  final_dataset$ineff_abx <- ifelse(final_dataset$ineff_abx == 1 & (final_dataset$maybe_eff_abx == 1 | final_dataset$who_rec_abx == 1), 0, final_dataset$ineff_abx)
  final_dataset$maybe_eff_abx <- ifelse(final_dataset$who_rec_abx == 1 & final_dataset$maybe_eff_abx == 1, 0, final_dataset$maybe_eff_abx)
  final_dataset$no_abx <- ifelse(final_dataset$who_rec_abx == 0 &
                                   final_dataset$maybe_eff_abx == 0 &
                                   final_dataset$ineff_abx == 0, 1, 0)
  
  # return final dataset of pids who received abx
  return(final_dataset)
  
}

#' Function to create death/rehospitalization outcome
#' Re-used from Allison O. OTR project
#' 
create_hosp90death <- function(){
  DCS_follow <- readRDS(here::here("data/efgh_data/raw_data/DCS_07_follow_up.rds")) %>%
    select(pid, foll_visit, foll_new_dia_adm, foll_new_ill_adm,
           foll_new_dia_adv_adm, foll_new_ill_adv_adm)
  #Consider sensitivity analysis to include recommended for hospitalization vs. admission
  #n=10 recommended for admission but not admitted for new diarrheal illness
  #n=25 recommended for admission but not admitted for new non-diarrheal illness
  
  # Transform dataset to wide data
  DCS_follow <- DCS_follow %>%
    mutate(foll_visit = case_when(
      foll_visit == 1 ~ "1mo",
      foll_visit == 2 ~ "3mo"
    ))
  
  DCS_follow <- pivot_wider(DCS_follow, 
                            names_from = foll_visit, 
                            values_from = c(foll_new_dia_adm, foll_new_ill_adm, foll_new_dia_adv_adm, foll_new_ill_adv_adm))
  
  DCS_hosp_record <- readRDS(here::here("data/efgh_data/raw_data/DCS_09_hospital_record_abstraction.rds")) %>%
    select(pid, hosp_ward, hosp_efgh_out)
  
  DCS_unwell <- readRDS(here::here("data/efgh_data/raw_data/DCS_08_unwell_child.rds")) %>%
    select(pid, uv_outcome, uv_outcome_oth)
  
  DCS_mortality <- readRDS(here::here("data/efgh_data/raw_data/DCS_10a_mortality.rds")) %>%
    select(pid, mort_date)
  
  # ----- Combine datasets for distinct observations
  dfs <- list(DCS_follow, DCS_hosp_record, DCS_unwell, DCS_mortality)
  
  hosp90death <- Reduce(function(x,y) merge(x, y, by = "pid", all =TRUE), dfs)
  
  hosp90death <- hosp90death %>%
    distinct()
  
  # ----- Create hierarchy to prioritize relevant uv_outcome observations
  hosp90death <- hosp90death %>%
    group_by(pid) %>%
    mutate(priority_uv_outcome = min(uv_outcome)) %>%
    ungroup()
  
  hosp90death <- hosp90death %>%
    distinct( pid, foll_new_dia_adm_1mo, foll_new_dia_adm_3mo, 
              foll_new_ill_adm_1mo, foll_new_ill_adm_3mo,
              foll_new_dia_adv_adm_1mo, foll_new_dia_adv_adm_3mo,
              foll_new_ill_adv_adm_1mo, foll_new_ill_adv_adm_3mo, mort_date,
              priority_uv_outcome, .keep_all = TRUE) %>%
    select(-uv_outcome, -uv_outcome_oth)
  
  # ----- Identify those who completed second follow-up (for whom outcome is known)
  follow <- readRDS(here::here("data/efgh_data/raw_data/DCS_07_follow_up.rds")) %>%
    select(pid, foll_visit) %>%
    filter(foll_visit==2)
  
  pid_complete <- unique(follow$pid)
  
  # ----- Construct outcome
  hosp90death <- hosp90death %>%
    mutate(hosp90 = NA,
           hosp90 = ifelse(priority_uv_outcome==1 |
                             hosp_ward==1 | #short-stay wards excluded from outcome
                             foll_new_dia_adm_1mo==1 |
                             foll_new_dia_adm_3mo==1 |
                             foll_new_ill_adm_1mo==1 |
                             foll_new_ill_adm_3mo==1, 1, hosp90),
           
           death90 = NA,
           death90 = ifelse(hosp_efgh_out==3 |
                              is.na(mort_date)==FALSE, 1, death90))
  
  
  hosp90death <- hosp90death %>%
    mutate(an_hosp90death = NA,
           an_hosp90death = ifelse(hosp90==1 | death90==1, 1, an_hosp90death),
           
           an_hosp90death = if_else(is.na(an_hosp90death)==TRUE &
                                      pid %in% pid_complete, 0, an_hosp90death),
           # allison c. added for hosp90 and death90 alone 9/1/26
           hosp90 = if_else(is.na(hosp90)==TRUE &
                                      pid %in% pid_complete, 0, hosp90),
           death90 = if_else(is.na(death90)==TRUE &
                                      pid %in% pid_complete, 0, death90))
  
  outcomes_final <- hosp90death %>%
    select(pid, an_hosp90death, death90, hosp90)
  
  return(outcomes_final)
}


#' Function to create clean EFGH case dataset
#' 
#' By default, this will return all diarrhea episodes (modified for co-etiology 
#' meta-analysis). For the original meta-analysis, we define case as TAC 
#' or culture attributable Shigella so subset the dataset from this function to those who
#' have Shigella TAC results.
prep_efgh <- function(){
  
  enroll_data <- readRDS(here::here("data/efgh_data/raw_data/DCS_04_enrollment.Rds")) %>%
    select(pid, 
           enroll_site, 
           enroll_diar_blood,
           enroll_diar_vom_days,
           enroll_diar_fever,
           enroll_diar_fever_days,
           enroll_diar_vom_num,
           enroll_diar_loose_num,
           enroll_cond_dehyd,
           enroll_cg_moth_ed, 
           enroll_ai_num_child,
           enroll_ai_wt,
           enroll_ai_toi,
           enroll_date) %>%
    mutate(across(everything(), ~ replace_na(., 0)))
  
  enroll_data$enroll_diar_blood <- ifelse(enroll_data$enroll_diar_blood == 3, NA, enroll_data$enroll_diar_blood)
  enroll_data$enroll_diar_blood <- ifelse(enroll_data$enroll_diar_blood == 1, 1, 0)
  
  enroll_data$enroll_date <- as.Date(enroll_data$enroll_date)
  
  # EFGH Water type
  # enroll_data$enroll_ai_wt <- factor(enroll_data$enroll_ai_wt, levels = 1:18,
  #                                    labels = c("Piped water into dwelling",        #improved
  #                                               "Piped into compound/yard/plot",    #improved
  #                                               "Piped to neighbor",                #improved
  #                                               "Public tap/standpipe",             #improved
  #                                               "Tube well or borehold",            #improved
  #                                               "Protected well",                   #improved
  #                                               "Unprotected well",                 #unimproved
  #                                               "Water from spring (protected)",    #improved
  #                                               "Water from spring (unprotected)",  #unimproved
  #                                               "Rainwater",                        #unimproved
  #                                               "Tanker-truck",                     #improved
  #                                               "Cart with small tank/drum",        #improved
  #                                               "Water kiosk (filtered)",           #improved
  #                                               "Water kiosk (unfiltered)",         #unimproved
  #                                               "Bottled water",                    #improved
  #                                               "Sachet water",                     #improved
  #                                               "Surface water (river/dam/lake/pond/canal/irrigation channel)", #unimproved
  #                                               "Other"))                           #improved
  
  enroll_data$imp_water <- ifelse(enroll_data$enroll_ai_wt %in% c(1,2,3,4,5,6,8,11,12,13,15,16,18), 1, 0)
  
  # MAL-ED definitions: (trying to match closely)
  # piped into dwelling = 1
  # piped to yard / plot = 2
  # public stand/pipe = 3
  # tube well or borehole = 4
  # protected well = 5
  # unprotected well = 6
  # surface water = 7
  # other = 8
  
  # 1,2,3,4,5,8 = improved
  # 6,7 = unimproved
  
  # Toilet type
  # enroll_data$enroll_ai_toi <- factor(enroll_data$enroll_ai_toi, levels = 1:12,
  #                                     labels = c("Flushed or pour flush to piped sewer system", #improved
  #                                                "Flushed or pour flush to septic tank",        #improved
  #                                                "Flushed or pour flush to somewhere else",     #improved
  #                                                "Flushed or pour flush, don't know where",     #improved
  #                                                "Ventilated improved pit latrine",             #improved
  #                                                "Pit latrine with slab floor",                 #improved
  #                                                "Pit latrine without slab floor",              #improved
  #                                                "Composting toilet",                           #improved
  #                                                "Bucket toilet",                               #unimproved
  #                                                "Hanging toilet/latrine",                      #improved
  #                                                "Open defecation/bush/field",                  #unimproved
  #                                                "Other"))                                      #improved
  enroll_data$imp_toi <- ifelse(enroll_data$enroll_ai_toi %in% c(1:8,10,12), 1, 0)
  
  # MAL-ED definitions:
  # no facility/bush/bucket/field = 1
  # pit no flush = 2
  # flushed to sewer = 3
  # flush to septic = 4
  # flush to pit = 5
  # flush elsewhere = 6
  # other = 7
  
  # 1 = unimproved
  # 2,3,4,5,6,7 = improved
  
  enroll_data$enroll_cond_dehyd <- factor(enroll_data$enroll_cond_dehyd,
                                          levels = 1:3,
                                          labels = c("Severe dehydration",
                                                     "Some dehydration",
                                                     "No dehydration"))
  
  enroll_data$enroll_site <- factor(enroll_data$enroll_site,
                                    levels = 1:7,
                                    labels = c("Bangladesh",
                                               "Kenya",
                                               "Malawi",
                                               "Mali",
                                               "Pakistan",
                                               "Peru",
                                               "The Gambia"))
  
  # NOTE also seems like UW analysis excluded Koranic school only but leaving for now
  # 8 = NA, 6 = declined to answer 
  enroll_data$enroll_cg_moth_ed <- ifelse(enroll_data$enroll_cg_moth_ed == 8 | enroll_data$enroll_cg_moth_ed == 6, NA, enroll_data$enroll_cg_moth_ed)
  enroll_data$enroll_cg_moth_ed <- factor(enroll_data$enroll_cg_moth_ed,
                                          levels = c(1:5, 7),
                                          labels = c("None",
                                                     "Less than primary school",
                                                     "Primary school only",
                                                     "Some secondary school",
                                                     "Secondary school or greater",
                                                     "Koranic school only"))
  
  # binary education to try to match <6 years, >= 6 years
  enroll_data$moth_ed_bin <- ifelse(
    is.na(enroll_data$enroll_cg_moth_ed), 
    NA, 
    ifelse(
      enroll_data$enroll_cg_moth_ed %in% c(
        "None", 
        "Less than primary school", 
        "Primary school only",
        "Koranic school only"
      ), 
      0, 
      1
    )
  )
  
  # Duration of diarrhea pre enrollment
  screen_data <- readRDS(here::here("data/efgh_data/raw_data/DCS_02_screening.Rds")) %>%
    select(serial_id, 
           scr_diar_days) %>%
    rename('duration_pre_enroll' = scr_diar_days)
  
  preenroll_data <- readRDS(here::here('data/efgh_data/raw_data/DCS_03_preenrollment.Rds')) %>%
    select(serial_id, 
           pid)
  
  screen_data <- left_join(screen_data, preenroll_data, by = 'serial_id')   %>%
    select(-serial_id) %>%
    drop_na()
  
  laz_data <- readRDS(here::here("data/efgh_data/raw_data/DCS_anthro.Rds")) %>%
    select(pid,
           sex,
           enr_age_months,
           enr_haz,
           enr_whz, 
           enr_waz, 
           enr_weight_before_kg,
           enr_lenhei_cm, 
           wk4_haz,
           mo3_lenhei,
           mo3_wt,
           mo3_haz,
           mo3_waz,
           mo3_whz,
           wk4_visit_dt,
           mo3_visit_dt)
  
  laz_data$enr_haz <- ifelse(laz_data$enr_haz > 6 | laz_data$enr_haz < -6, NA, laz_data$enr_haz)
  laz_data$wk4_haz <- ifelse(laz_data$wk4_haz > 6 | laz_data$wk4_haz < -6, NA, laz_data$wk4_haz)
  laz_data$mo3_haz <- ifelse(laz_data$mo3_haz > 6 | laz_data$mo3_haz < -6, NA, laz_data$mo3_haz)
  
  laz_data$enr_waz <- ifelse(laz_data$enr_waz > 6 | laz_data$enr_waz < -6, NA, laz_data$enr_waz)
  laz_data$mo3_waz <- ifelse(laz_data$mo3_waz > 6 | laz_data$mo3_waz < -6, NA, laz_data$mo3_waz)
  
  laz_data$enr_whz <- ifelse(laz_data$enr_whz > 6 | laz_data$enr_whz < -6, NA, laz_data$enr_whz)
  laz_data$mo3_whz <- ifelse(laz_data$mo3_whz > 6 | laz_data$mo3_whz < -6, NA, laz_data$mo3_whz)
  
  laz_data$wk4_visit_dt <- as.Date(laz_data$wk4_visit_dt)
  laz_data$mo3_visit_dt <- as.Date(laz_data$mo3_visit_dt)
  
  # NEW 11/11/24 -- add in breastfeeding variable 
  # using derivation code from UW, binary version of variable
  breast_data <- readRDS(here::here("data/efgh_data/raw_data/DCS_04_enrollment.Rds")) %>%
    left_join(laz_data %>% select(pid,enr_age_months)) %>%
    select(pid,enroll_hh_no_56,enroll_hh_breast,enroll_hh_yes_fluid,enroll_hh_yes_food, enr_age_months) %>%
    mutate(enroll_hh_yes_fluid = ifelse(enroll_hh_yes_fluid==99, NA, enroll_hh_yes_fluid),
           enroll_hh_yes_food = ifelse(enroll_hh_yes_food==99, NA, enroll_hh_yes_food)) %>%
    rowwise() %>%
    mutate(age_exclusive_breastfeeding = suppressWarnings(min(enroll_hh_yes_fluid,enroll_hh_yes_food,na.rm = T)),
           age_exclusive_breastfeeding = ifelse(age_exclusive_breastfeeding=="Inf", NA, age_exclusive_breastfeeding),
           age_exclusive_breastfeeding = ifelse(enroll_hh_breast==2, 0, age_exclusive_breastfeeding)) %>%
    ungroup() %>%
    mutate(exclusive_breastfeeding_bin = ifelse(age_exclusive_breastfeeding<6, 0,
                                                ifelse(age_exclusive_breastfeeding>=6, 1, NA)),
           exclusive_breastfeeding_bin = ifelse(enroll_hh_no_56==1 & is.na(exclusive_breastfeeding_bin), 1, exclusive_breastfeeding_bin)) %>%
    select(pid, exclusive_breastfeeding_bin)
  
  wealth_data <- readRDS(here::here("data/efgh_data/raw_data/DCS_wealth_index.Rds")) %>%
    left_join(data.frame(pid = enroll_data$pid,
                         site = enroll_data$enroll_site), by = "pid") %>%
    group_by(site) %>%
    mutate(final_quintile_site = ntile(final_index, 5)) %>%
    ungroup() %>%
    select(-c("site"))
  
  wealth_data$final_quintile_site <- factor(wealth_data$final_quintile_site, levels = 1:5, labels = c("1st quintile of SES",
                                                                                                      "2nd quintile of SES",
                                                                                                      "3rd quintile of SES",
                                                                                                      "4th quintile of SES",
                                                                                                      "5th quintile of SES"))
  
  # Use helper function to get abx data
  # Returns dataset containing PIDs who received (any) abx, who abx, and abx indicated for diarrhea
  # Uses same cleaning method from UW analysis but missing light-touch information
  abx_data <- clean_abx()
  abx_data$pid <- as.numeric(abx_data$pid)
  
  shigella_data <- readRDS(here::here("data/efgh_data/raw_data/DCS_shigella.Rds")) %>%
    select("pid",
           "positive_tac_or_culture",
           "culture_shigella_positive",
           "tac_shigella_detected",
           "tac_shigella_attributable")
  
  # SHIGELLA ONLY ONE WITH CULTURE, rest tac only
  # some people have both rectal and whole stool, does not look like whole stool used to define *_bin variables so drop 
  
  # NEW VERSION USING TAC PROCESSED FROM DATA LOCK
  tac_data <- readRDS(here::here("data/efgh_data/raw_data/tac_processed.Rds")) %>%
    filter(shigella_sample_type == "Rectal") %>%
    select("pid",
           "shigella_attributable",
           "shigella_serotype",
           "adenovirus_40_41_attributable",
           "aeromonas_attributable",
           "astrovirus_attributable",
           "c_jejuni_coli_attributable",
           "cryptosporidium_attributable",
           "cyclospora_attributable",
           "e_histolytica_attributable",
           "norovirus_gii_attributable",
           "rotavirus_attributable",
           "salmonella_attributable",
           "sapovirus_attributable",
           "ST.ETEC_attributable",
           "tEPEC_attributable",
           "v_cholerae_attributable",
           "isospora_attributable",
           "rotavirus_ct",
           "shigella_ct",
           "adenovirus_40_41_ct",
           "ETEC_ct",
           "cryptosporidium_ct",
           "astrovirus_ct",
           "norovirus_gi_ct",
           "norovirus_gii_ct",
           "c_jejuni_coli_ct",
           "tEPEC_ct",
           "sapovirus_ct",
           "e_bieneusi_ct",
           "giardia_ct",
           "EAEC_ct",
           "ST.ETEC_ct",
           "v_cholerae_ct",
           "salmonella_ct") %>%
    # Binary shigella flex and shigella sonnei
    mutate(
      shigella_flex = ifelse(
        is.na(shigella_serotype),
        NA_real_,
        as.numeric(grepl("S\\.flexneri", shigella_serotype))
      ),
      shigella_sonnei = ifelse(
        is.na(shigella_serotype),
        NA_real_,
        as.numeric(grepl("S\\.sonnei", shigella_serotype))
      )) %>%
    left_join(shigella_data[,c("pid", 
                               "tac_shigella_attributable",
                               "positive_tac_or_culture")], by = "pid") %>% 
    mutate(
      # Transform TAC to scale starting at 0; (35 - CTValue) / 3.322
      shigella_new = (35 - shigella_ct) / 3.322,
      rotavirus_new = (35 - rotavirus_ct) / 3.322,
      adenovirus_new = (35 - adenovirus_40_41_ct) / 3.322,
      ETEC_new = (35 - ETEC_ct) / 3.322,
      st_etec_new = (35 - ST.ETEC_ct) / 3.322,
      cryptosporidium_new = (35 - cryptosporidium_ct) / 3.322,
      astrovirus_new = (35 - astrovirus_ct) / 3.322,
      norovirus_gii_new = (35 - norovirus_gii_ct) / 3.322,
      c_jejuni_new = (35 - c_jejuni_coli_ct) / 3.322,
      tEPEC_new = (35 - tEPEC_ct) / 3.322,
      sapovirus_new = (35 - sapovirus_ct) / 3.322,
      e_bieneusi_new = (35 - e_bieneusi_ct) / 3.322,
      giardia_new = (35 - giardia_ct) / 3.322,
      EAEC_new = (35 - EAEC_ct) / 3.322,
      v_cholerae_new = (35 - v_cholerae_ct) / 3.322,
      salmonella_new = (35 - salmonella_ct) / 3.322
    ) %>%
    mutate(
      bacteria_attr = if_else(
        coalesce(shigella_attributable, 0) == 1 |
          coalesce(ST.ETEC_attributable, 0) == 1 |
          coalesce(c_jejuni_coli_attributable, 0) == 1 |
          coalesce(salmonella_attributable, 0) == 1 |
          coalesce(tEPEC_attributable, 0) == 1 |
          coalesce(v_cholerae_attributable, 0) == 1 |
          coalesce(aeromonas_attributable, 0) == 1 |
          coalesce(cryptosporidium_attributable, 0) == 1 |
          coalesce(cyclospora_attributable, 0) == 1 |
          coalesce(e_histolytica_attributable, 0) == 1 |
          coalesce(isospora_attributable, 0) == 1,
        1, 0
      ),
      viral_attr = if_else(
        coalesce(rotavirus_attributable, 0) == 1 |
          coalesce(adenovirus_40_41_attributable, 0) == 1 |
          coalesce(astrovirus_attributable, 0) == 1 |
          coalesce(norovirus_gii_attributable, 0) == 1 |
          coalesce(sapovirus_attributable, 0) == 1,
        1, 0
      )
    )
  
  # Susceptibility data general (all people, not just those who got abx)
  ast_data <- readRDS(here::here("data/efgh_data/raw_data/DCS_16_rectal_swab_results.Rds")) %>%
    select(pid,
           swi_ast_az_int,
           swi_ast_ceft_int,
           swi_ast_cipro_int,
           swi_ast_piv_int) %>%
    mutate(swi_ast_az_int = if_else(swi_ast_az_int == 1, "S",
                                    if_else(swi_ast_az_int == 2 | swi_ast_az_int == 3, "I/R", NA)),
           swi_ast_cipro_int = if_else(swi_ast_cipro_int == 1, "S",
                                       if_else(swi_ast_cipro_int == 2 | swi_ast_cipro_int == 3, "I/R", NA)),
           swi_ast_ceft_int = if_else(swi_ast_ceft_int == 1, "S",
                                      if_else(swi_ast_ceft_int == 2 | swi_ast_ceft_int == 3, "I/R", NA)),
           swi_ast_piv_int = if_else(swi_ast_piv_int == 1, "S",
                                     if_else(swi_ast_piv_int == 2 | swi_ast_piv_int == 3, "I/R", NA))) %>%
    group_by(pid) %>%
    #group by PID, if there are multiple entries, only make it "S" if all are S
    summarise(
      swi_ast_az_int = case_when(
        all(is.na(swi_ast_az_int)) ~ NA_character_,
        all(swi_ast_az_int == "S", na.rm = TRUE) ~ "S",
        any(swi_ast_az_int == "I/R", na.rm = TRUE) ~ "I/R",
        TRUE ~ NA_character_
      ),
      swi_ast_cipro_int = case_when(
        all(is.na(swi_ast_cipro_int)) ~ NA_character_,
        all(swi_ast_cipro_int == "S", na.rm = TRUE) ~ "S",
        any(swi_ast_cipro_int == "I/R", na.rm = TRUE) ~ "I/R",
        TRUE ~ NA_character_
      ),
      swi_ast_ceft_int = case_when(
        all(is.na(swi_ast_ceft_int)) ~ NA_character_,
        all(swi_ast_ceft_int == "S", na.rm = TRUE) ~ "S",
        any(swi_ast_ceft_int == "I/R", na.rm = TRUE) ~ "I/R",
        TRUE ~ NA_character_
      ),
      swi_ast_piv_int = case_when(
        all(is.na(swi_ast_piv_int)) ~ NA_character_,
        all(swi_ast_piv_int == "S", na.rm = TRUE) ~ "S",
        any(swi_ast_piv_int == "I/R", na.rm = TRUE) ~ "I/R",
        TRUE ~ NA_character_
      )
    ) %>%
    rowwise() %>%
    mutate(
      resistant_WHO_approve = if_else(
        any(c_across(c(swi_ast_az_int, swi_ast_cipro_int, swi_ast_ceft_int, swi_ast_piv_int)) == "I/R", na.rm = TRUE),
        1, 0
      ),
      susceptible_WHO_approve = if_else(
        any(c_across(c(swi_ast_az_int, swi_ast_cipro_int, swi_ast_ceft_int, swi_ast_piv_int)) == "S", na.rm = TRUE),
        1, 0
      )
    ) %>%
    ungroup()
  
  # Find diarrhea with no etiology
  tac_data$no_etiology <- ifelse(rowSums(tac_data[, c(
    "positive_tac_or_culture", # based on TAC or culture definition
    "rotavirus_attributable",
    "ST.ETEC_attributable",
    "cryptosporidium_attributable",
    "adenovirus_40_41_attributable",
    "astrovirus_attributable",
    "c_jejuni_coli_attributable",
    "cyclospora_attributable",
    "e_histolytica_attributable",
    "norovirus_gii_attributable",
    "salmonella_attributable",
    "sapovirus_attributable",
    "tEPEC_attributable",
    "v_cholerae_attributable",
    "isospora_attributable",
    "aeromonas_attributable"
  )], na.rm = TRUE) == 0, 1, 0)
  
  # remove before remerge
  tac_data <- tac_data %>%
    select(-tac_shigella_attributable, -positive_tac_or_culture, -shigella_attributable)
  
  # Get severity information (GEMS_MSD)
  # ADDED 9/1/26 -- get duration of episode after enrollment by taking episode diarrhea length - diar_days_at_enroll
  sev_data <- readRDS(here::here("data/efgh_data/raw_data/DCS_severity_scores.Rds")) %>%
    mutate("duration_post_enroll" = episode_diarrhea_length - diar_days_at_enroll) %>%
    select("pid", "gems_msd", "duration_post_enroll")
  
  combo_data <- left_join(tac_data, laz_data, by = "pid") %>%
    merge(wealth_data, by = "pid") %>%
    merge(enroll_data, by = "pid") %>%
    left_join(abx_data, by = "pid") %>%
    merge(shigella_data, by = "pid") %>%
    left_join(breast_data, by = "pid") %>%
    left_join(screen_data, by = 'pid') %>%
    left_join(ast_data, by = 'pid') %>%
    left_join(sev_data, by = 'pid')
  
  # Drop people with NA for tac_attributable_shigella
  # NOW KEEP THESE PEOPLE FOR CO-ETIOLOGY ANALYSIS
  # combo_data <- combo_data[-which(is.na(combo_data$tac_shigella_attributable)),]
  
  # if was not in abx dataset, did not receive abx
  combo_data$no_abx <- ifelse(is.na(combo_data$no_abx), 1, combo_data$no_abx)
  combo_data$who_rec_abx <- ifelse(is.na(combo_data$who_rec_abx), 0, combo_data$who_rec_abx)
  combo_data$maybe_eff_abx <- ifelse(is.na(combo_data$maybe_eff_abx), 0, combo_data$maybe_eff_abx)
  combo_data$ineff_abx <- ifelse(is.na(combo_data$ineff_abx), 0, combo_data$ineff_abx)
  combo_data$anti_diarrhea_abx <- ifelse(is.na(combo_data$anti_diarrhea_abx), 0, combo_data$anti_diarrhea_abx)
  combo_data$azithromycin <- ifelse(is.na(combo_data$azithromycin), 0, combo_data$azithromycin)
  combo_data$ciprofloxacin <- ifelse(is.na(combo_data$ciprofloxacin), 0, combo_data$ciprofloxacin)
  combo_data$ceftriaxone <- ifelse(is.na(combo_data$ceftriaxone), 0, combo_data$ceftriaxone)
  combo_data$pivmecillinam <- ifelse(is.na(combo_data$pivmecillinam), 0, combo_data$pivmecillinam)
  
  combo_data$all_abx <- ifelse(combo_data$no_abx == 1 | combo_data$ineff_abx == 1, 0, 
                               ifelse(combo_data$maybe_eff_abx == 1, 1, 2))
  
  combo_data$all_abx <- factor(combo_data$all_abx, levels = 0:2, labels = c("No or ineffective antibiotics",
                                                                            "Possibly effective antibiotics",
                                                                            "Guideline recommended antibiotics"))
  
  # Vars for delta wk4 and delta mo3
  combo_data$enr_wk4_delta_haz <- combo_data$wk4_haz - combo_data$enr_haz
  combo_data$enr_mo3_delta_haz <- combo_data$mo3_haz - combo_data$enr_haz
  
  # Vars for other deltas
  combo_data$enr_mo3_delta_waz <- combo_data$mo3_waz - combo_data$enr_waz
  combo_data$enr_mo3_delta_whz <- combo_data$mo3_whz - combo_data$enr_whz
  combo_data$enr_mo3_delta_lenhei <- combo_data$mo3_lenhei - combo_data$enr_lenhei_cm
  combo_data$enr_mo3_delta_weight <- combo_data$mo3_wt - combo_data$enr_weight_before_kg
  
  # NEW 12/2/24 - add covariate for days to follow-up 
  combo_data$wk4_days <- as.Date(combo_data$wk4_visit_dt) - as.Date(combo_data$enroll_date)
  combo_data$mo3_days <- as.Date(combo_data$mo3_visit_dt) - as.Date(combo_data$enroll_date)
  
  # Add indicator that follow-up days is present
  combo_data$I_mo3_days <- ifelse(is.na(combo_data$mo3_days), 0, 1)
  combo_data$I_mo3_days_x_mo3_days <- ifelse(is.na(combo_data$I_mo3_days * combo_data$mo3_days), 
                                             0, 
                                             combo_data$I_mo3_days * combo_data$mo3_days)
  
  # Add death and rehospitalization variable - reusing code from OTR project
  death_rehosp <- create_hosp90death()
  # note there are people missing from death_rehosp that appear in combo_data
  # these people are also missing other outcomes (lost to follow up)
  combo_data <- combo_data %>%
    left_join(death_rehosp, by = "pid")
  
  # Add labels to dataset for nice display in gtsummary tables
  combo_data <- combo_data %>%
    select(pid,
           enroll_site,
           enr_age_months, 
           sex, 
           enr_haz,
           enr_whz,
           enr_waz,
           enr_weight_before_kg,
           enr_lenhei_cm, 
           wk4_haz,
           wk4_visit_dt, 
           mo3_haz,
           mo3_waz, 
           mo3_whz, 
           mo3_wt, 
           mo3_lenhei,
           mo3_wt,
           mo3_visit_dt,  
           enr_wk4_delta_haz,
           enr_mo3_delta_haz,
           enr_mo3_delta_waz,
           enr_mo3_delta_whz,
           enr_mo3_delta_lenhei,
           enr_mo3_delta_weight,
           wk4_days,                     
           mo3_days,
           I_mo3_days,
           I_mo3_days_x_mo3_days,
           death90,
           hosp90,
           an_hosp90death, 
           shigella_serotype,
           adenovirus_40_41_attributable,
           aeromonas_attributable,
           astrovirus_attributable,
           c_jejuni_coli_attributable,
           cryptosporidium_attributable,
           cyclospora_attributable,
           e_histolytica_attributable,
           norovirus_gii_attributable,
           rotavirus_attributable,
           salmonella_attributable,
           sapovirus_attributable,
           ST.ETEC_attributable,
           tEPEC_attributable,
           v_cholerae_attributable,
           isospora_attributable,
           rotavirus_ct,
           shigella_ct,
           adenovirus_40_41_ct,
           ETEC_ct,
           cryptosporidium_ct,
           astrovirus_ct,
           norovirus_gi_ct,
           norovirus_gii_ct,
           c_jejuni_coli_ct,
           tEPEC_ct,
           sapovirus_ct,
           e_bieneusi_ct,
           giardia_ct,
           EAEC_ct,
           ST.ETEC_ct,
           v_cholerae_ct,
           salmonella_ct,
           shigella_flex,
           shigella_sonnei,
           shigella_new,
           rotavirus_new,
           adenovirus_new,
           ETEC_new,
           st_etec_new,
           cryptosporidium_new,
           astrovirus_new,
           norovirus_gii_new,
           c_jejuni_new,
           tEPEC_new,
           sapovirus_new,
           e_bieneusi_new,
           giardia_new,
           EAEC_new,
           v_cholerae_new,
           salmonella_new,
           bacteria_attr,
           viral_attr,
           no_etiology,
           final_quintile,
           final_quintile_site,
           exclusive_breastfeeding_bin,
           enroll_cg_moth_ed,
           enroll_ai_num_child,
           enroll_ai_wt,
           enroll_ai_toi,
           enroll_date,
           imp_water,
           imp_toi,
           moth_ed_bin,
           duration_pre_enroll,
           duration_post_enroll,
           gems_msd,
           enroll_diar_blood,
           enroll_diar_vom_days,
           enroll_diar_fever,
           enroll_diar_fever_days,
           enroll_diar_vom_num,
           enroll_diar_loose_num,
           enroll_cond_dehyd,
           all_abx,
           who_rec_abx,
           maybe_eff_abx,                
           ineff_abx, 
           anti_diarrhea_abx,
           azithromycin,
           ciprofloxacin,
           ceftriaxone,
           pivmecillinam,
           no_abx, 
           ast_given_abx,
           swi_ast_az_int,
           swi_ast_cipro_int,
           swi_ast_ceft_int,
           swi_ast_piv_int,
           resistant_WHO_approve,
           susceptible_WHO_approve,
           ast_given_abx,
           positive_tac_or_culture,
           culture_shigella_positive,
           tac_shigella_detected,
           tac_shigella_attributable
           ) %>%
    set_variable_labels(pid = "Participant ID (for given episode)",
                        sex = "Sex",
                        enr_age_months = "Age at enrollment (months)",
                        enr_haz = "HAZ at enrollment",
                        enr_whz = "WHZ at enrollment",
                        enr_waz = "WAZ at enrollment",
                        enr_weight_before_kg = "Mean pre-hydration weight at enrollment",
                        enr_lenhei_cm = "Mean length at enrollment",
                        wk4_haz = "Week 4 HAZ",
                        mo3_haz = "Month 3 HAZ",
                        wk4_visit_dt = "Date of week 4 visit",
                        mo3_visit_dt = "Date of month 3 visit",
                        mo3_waz = "Month 3 WAZ",
                        mo3_whz = "Month 3 WHZ",
                        mo3_wt = "Month 3 weight",
                        mo3_lenhei = "Month 3 length",
                        mo3_visit_dt = "Month 3 visit date",
                        an_hosp90death = "Month 3 rehospitalization or death",
                        death90 = "Death by day 90",
                        hosp90 = "Hospitalization by day 90",
                        wk4_days = "Days between enrollment and week 4 visit",
                        mo3_days = "Days between enrollment and month 3 visit",
                        duration_post_enroll = "Duration of diarrhea after enrollment (>=3 watery stools in 24hr period)",
                        final_quintile = "SES Quintile",
                        final_quintile_site = "SES Quintile by site",
                        enroll_site = "Enrollment site",
                        enroll_cg_moth_ed = "Mother's education",
                        moth_ed_bin = "Mother's education more than primary school",
                        anti_diarrhea_abx = "Received antibiotics with diarrheal indication",
                        who_rec_abx = "Received WHO approved antibiotics",
                        maybe_eff_abx = "Maybe effective antibiotics",
                        ineff_abx = "Ineffective antibiotics",
                        no_abx = "No antibiotics",
                        all_abx = "Type of antibiotics received",
                        azithromycin = "Azithromycin",
                        ciprofloxacin = "Ciprofloxacin",
                        ceftriaxone = "Ceftriaxone",
                        pivmecillinam = "Pivmecillinam",
                        ast_given_abx = "S or I/R for given antibiotic(s)",
                        resistant_WHO_approve = "Resistant to one or more WHO approved abx",
                        susceptible_WHO_approve = "Susceptible for one or more WHO approved abx",
                        positive_tac_or_culture = "Shigella positive PCR or culture",
                        culture_shigella_positive = "Shigella positive culture",
                        tac_shigella_detected = "Shigella detected by PCR",
                        tac_shigella_attributable = "Shigella attributable by PCR",
                        enroll_diar_blood = "Dysentery pre-enrollment",
                        enroll_diar_loose_num = "Max number of loose stools on worst day pre-enrollment",
                        enroll_diar_vom_days = "Number of days vomiting pre-enrollment",
                        enroll_diar_vom_num = "Number of vomiting episodes on worst vomiting day pre-enrollment",
                        enroll_diar_fever = "Caregiver reported fever",
                        enroll_diar_fever_days = "Number of days of fever pre-enrollment",
                        enroll_ai_num_child = "Number of children in household under age 5",
                        duration_pre_enroll = "Diarrhea duration prior to enrollment",
                        imp_water = "Improved water",
                        imp_toi = "Improved toilets",
                        enroll_cond_dehyd = "Dehydration assessment",
                        rotavirus_attributable = "Rotavirus attributable",
                        ST.ETEC_attributable = "ST ETEC attributable",
                        cryptosporidium_attributable = "Cryptosporidium attributable",
                        adenovirus_40_41_attributable = "Adenovirus attributable",
                        astrovirus_attributable = "Astrovirus attributable",
                        c_jejuni_coli_attributable = "C jejuni attributable",
                        cyclospora_attributable = "Cyclospora attributable",
                        e_histolytica_attributable = "E histolytica attributable",
                        norovirus_gii_attributable = "Norovirus GII attributable",
                        salmonella_attributable = "Salmonella attributable",
                        sapovirus_attributable = "Sapovirus attributable",
                        tEPEC_attributable = "T EPEC attributable",
                        v_cholerae_attributable = "V Cholerae attributable",
                        isospora_attributable = "Cytoisospora belli attributable",
                        aeromonas_attributable = "Aeromonas attributable",
                        no_etiology = "Not TAC attributable to any pathogens",
                        shigella_new = "Shigella TAC scaled",
                        rotavirus_new = "Rotavirus TAC scaled",
                        adenovirus_new = "Adenovirus TAC scaled",
                        ETEC_new = "ETEC TAC scaled",
                        cryptosporidium_new = "Cryptosporidium TAC scaled",
                        astrovirus_new = "Astrovirus TAC scaled",
                        norovirus_gii_new = "Norovirus GII TAC scaled",
                        c_jejuni_new = "C Jejuni TAC scaled",
                        salmonella_new = "Salmonella scaled",
                        tEPEC_new = "tEPEC TAC scaled",
                        sapovirus_new = "Sapovirus TAC scaled",
                        e_bieneusi_new = "E Bieneusi TAC scaled",
                        giardia_new = "Giardia TAC scaled",
                        EAEC_new = 'EAEC TAC scaled',
                        st_etec_new = "ST ETEC scaled",
                        v_cholerae_new = "V. Cholerae scaled",
                        exclusive_breastfeeding_bin = "Exclusively breastfed for >= 6 months",
                        enr_wk4_delta_haz = "Change in HAZ - Week 4",
                        enr_mo3_delta_haz = "Change in HAZ - Month 3",
                        enr_mo3_delta_waz = "Change in WAZ - Month 3",
                        enr_mo3_delta_whz = "Change in WHZ - Month 3",
                        enr_mo3_delta_lenhei = "Change in length - Month 3",
                        enr_mo3_delta_weight = "Change in weight - Month 3",
                        wk4_days = "Days between enrollment and week 4 follow-up visit",
                        mo3_days = "Days between enrollment and month 3 follow-up visit",
                        I_mo3_days = "Indicator mo3_days not missing",
                        I_mo3_days_x_mo3_days = "Indicator mo3_days not missing * mo3_days",
                        gems_msd = "GEMS definition of moderate to severe diarrhea",
                        shigella_serotype = "Shgiella serotype",
                        shigella_flex = "Shigella flexneri serotype",
                        shigella_sonnei = "Shigella sonnei serotype", 
                        enroll_date = "Enrollment date",
                        swi_ast_az_int = "Azithromycin resistant",
                        swi_ast_cipro_int = "Ciprofloxacin resistant",
                        swi_ast_ceft_int = "Ceftriaxone resistant",
                        swi_ast_piv_int = "Pivmecillinam resistant",
                        enroll_ai_toi = "Toilet type",
                        enroll_ai_wt = "Water type")
  
  # ADD IN CHILD ID 
  # first ID is the first pid for a given child
  # pid = child id = associated with the episode
  child_id_data <- readRDS(here::here("data/efgh_data/raw_data/DCS_multiple_enrollments.Rds"))
  child_id_final <- data.frame(first_id = child_id_data$pid,
                               pid = child_id_data$pid1)
  
  
  child_id_final <- rbind(child_id_final,
                          data.frame(first_id = child_id_data$pid,
                                     pid = child_id_data$pid2) %>% drop_na())
  
  child_id_final <- rbind(child_id_final,
                          data.frame(first_id = child_id_data$pid,
                                     pid  = child_id_data$pid3) %>% drop_na())
  
  child_id_final <- rbind(child_id_final,
                          data.frame(first_id = child_id_data$pid,
                                     pid = child_id_data$pid4) %>% drop_na())
  
  combo_data <- left_join(combo_data, child_id_final, by = "pid") 
    
  # combo_data$child_id <- combo_data$first_id
  # fix line 9/1/26, had already fixed in IPD
  combo_data$child_id <- combo_data$pid
  
  combo_data <- combo_data %>%
    set_variable_labels(first_id = "PID from first enrollment (identifies unique child)",
                        child_id = "PID of episode (= pid, naming convention to match other studies)")
  
  return(combo_data)
}

efgh_data_full <- prep_efgh()
efgh_data <- efgh_data_full[which(!is.na(efgh_data_full$tac_shigella_attributable)), , drop = FALSE] # drop NA tac results

saveRDS(efgh_data, here::here("data/efgh_data/efgh_data.Rds"))
saveRDS(efgh_data_full, here::here("data/efgh_data/efgh_data_full.Rds"))
