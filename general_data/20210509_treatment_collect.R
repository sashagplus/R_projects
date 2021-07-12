library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library("rpart")
library(rpart.plot)
library(readxl)    
library(caret)

#library(gpio)

rm(list=ls())


#0. Read stard lines to g+ lines mapping
#gp = genetikaplus
patient_to_sample=read.table(str_c("/Users/sashakugel/gplus_dropbox/",
                                   "Genetika+ Dropbox/Genetika+SharedDrive/",
                                   "01_Protocol_Development/01_10_Data_Science/",
                                   "2020_GP_Patient_Codes-StarD.tsv"),
                             header=TRUE,
                             stringsAsFactors = F,
                             sep="\t")

patient_to_sample %<>% mutate(GP...CODE.clean=str_remove(GP...CODE, "\\*"))

patient_to_sample %<>% mutate(GP...CODE.numeric=as.numeric(str_remove(GP...CODE.clean, "LCL-")))

colnames(patient_to_sample) = c("id.stard", "id.ruidCellLine", "id.gpCode.raw", "id.gp", "id.gp.numeric")

#1. Read response data, and summarise

read_level=function(filename)
{
  levelx=read.csv(filename)
  
  levelx %<>% 
    mutate_at(vars(starts_with("Study.med.code.")), ~ as.character(.x)) %<>%
    mutate_at(vars(starts_with("Prescribed.med.code.")), ~ as.character(.x)) %<>%
  mutate_at(vars(starts_with("Prescribed.med.daily.dose.")), ~ as.character(.x)) %<>%
    mutate_at(vars(starts_with("Med.type")), ~ as.character(.x))
  
  
   return(levelx)
}

#level 1
pathname=stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                        "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/",
                        "DataforPatientHistory/DepressionDataNIH/",
                        "StarD/STAR_DClinicalData/ASCII Files for viewing/")

cc_levels_filenames=data.frame(levels=c("1","2", "2A", "3","4" ),
                               filename=c("Level1/Visits/CC_1.csv",
                                          "Level2/Visits/CC_2.csv",
                                          "Level2a/Visits/CC_2A.csv", 
                                          "Level3/Visits/CC_3.csv", 
                                          "Level4/Visits/CC_4.csv"))
alllevels.full=data.frame()
for(cc in 1:nrow(cc_levels_filenames))
{
  level_cc=read_level(stringr::str_c(pathname, cc_levels_filenames$filename[cc]))
  alllevels.full %<>% bind_rows(level_cc)
}


#Write combined output to csv for other team memebrs use
# write.csv(alllevels.full, 
#           "/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/Exploratory_Outputs/2021_05_13_CC_all_levels_combined.csv",
#           row.names = FALSE)
# 
# write.csv(alllevels.full %>% filter(Patient.ID %in% patient_to_sample$id.stard), 
#           "/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/Exploratory_Outputs/2021_05_13_CC_all_levels_combined_GPPatientsOnly.csv",
#           row.names = FALSE)


#alllevels=bind_rows(level1, level2, level2a, level3, level4)

alllevels=alllevels.full  %>% select(Patient.ID, 
                                     Treatment.level, 
                                     Week.in.level, 
                                     starts_with("Study.med.code."), 
                                     QIDS.C.percent.improvement..transcribed., 
                                     QIDS.C.beginning.score..transcribed., 
                                     QIDS.C.current.score..transcribed.) 

#Filter for genetikaplus patients
alllevels %<>% filter(Patient.ID %in% patient_to_sample$id.stard)

#Filter single drug therapies
# alllevels %<>% filter(Study.med.code.2==-2) %>% 
#                 filter(Study.med.code.1!=-2) 

#alllevels %<>% ungroup %<>% select(-Study.med.code.2, -Study.med.code.3, -Study.med.code.4) %<>% distinct()

alllevels %<>% mutate_at(vars(starts_with("Study.med")), 
                         ~ ifelse(.x %in% c("-2", "-3", "-5", "0"), NA, .x) )

treatment_counts=alllevels %>% group_by_at(vars(starts_with("Study.med"))) %>% 
  tally() 


treatment_counts.processd=treatment_counts %>% tibble::rowid_to_column("treatment") %>% 
  gather(key="Study", value="med", -n, -treatment ) %>% 
  group_by(treatment, n) %>% 
  summarise(treatment_string=toString(sort(unique(med[!is.na(med)]))), .groups = "keep") %>% 
  right_join(treatment_counts %>% tibble::rowid_to_column("treatment"), by=c("treatment", "n")) %>% 
  group_by(treatment_string) %>% 
  mutate(grp_indx=cur_group_id()) %>%
  ungroup() %>% 
  rowwise %>% 
  mutate(num_of_drugs=sum(!is.na(c_across(starts_with("Study.med")))))   





alllevels %<>% left_join(treatment_counts.processd ,
                         by=c("Study.med.code.1", "Study.med.code.2", "Study.med.code.3", "Study.med.code.4")) 
alllevels %<>% group_by(Patient.ID, Treatment.level, treatment_string) %>% mutate(max_week=max(Week.in.level))



#Mark the last week of each level
#alllevels %<>% group_by(Patient.ID, Treatment.level) %>% mutate(max_week=max(Week.in.level))

#In patients selected by GP, single drug treatments didn't have more then 2 recorded drugs per level.
#alllevels %>% group_by(Patient.ID, Treatment.level) %>% summarise(meds=n_distinct(Study.med.code.1)) %>% View

#Filter last week !!ADD IF WEEK>IN>LEVEL SINGLE LINE THEN TAKE IT
alllevels %<>% filter(Week.in.level != 0) %>% filter((!Week.in.level %in% c(0.1, 2)) | max_week==Week.in.level) 
#alllevels %<>% select(-max_week)

#' required output:
#' 1. Patient_ID
#' 2. Treatment Level
#' 3. Drug
#' 3. latest_week_in_level
#' 4. qids_start
#' 5. qids_current_week
#' 6. qids_improvement
#' 7. response_imp50p
#' 8. response_imp50p_bl10p
#' 9. response_imp50p_bl20p
#' 10. remission_p50imp_below6current


single_drugs= alllevels %>% group_by(Patient.ID, treatment_string ) %>% 
                            filter(Week.in.level == max(Week.in.level)) %>% 
                            slice(which.max(QIDS.C.percent.improvement..transcribed.)) #==max(QIDS.C.percent.improvement..transcribed.)) %>%
                            #distinct()

single_drugs %<>% select(Patient.ID,
                         Treatment.level,
                         latest_week_in_level=max_week,
                         Treatment_single_drug=treatment_string,
                         qids_start=QIDS.C.beginning.score..transcribed. ,
                         qids_current=QIDS.C.current.score..transcribed. ,
                         qids_improvment=QIDS.C.percent.improvement..transcribed.)


single_drugs %<>% mutate(response_imp50p=ifelse(qids_improvment>=50, "RESP", "NONRESP"),
                         remission_p50imp_below6current=ifelse(qids_improvment>=50 & qids_current<=6, "REMS", "NONREMS")) %<>% 
                  mutate(response_imp50p_bl10p=ifelse(qids_improvment<55 & qids_improvment>45 , 
                                                      stringr::str_c(response_imp50p, "BL"), response_imp50p),
                         response_imp50p_bl20p=ifelse(qids_improvment<60 & qids_improvment>40 , 
                                                      stringr::str_c(response_imp50p, "BL"), response_imp50p))


single_drugs %<>% left_join(patient_to_sample %>% select(id.stard, id.gp.numeric ), 
                            by=c("Patient.ID"="id.stard")) %<>% 
                  relocate(id.gp.numeric )

write.csv(single_drugs,
          stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                          "Genetika+SharedDrive/01_Protocol_Development/",
                          "01_10_Data_Science/Exploratory_Outputs/responsivness/",
                          "2021_06_14_responsiveness_labeling_qids_based_latest_allTreatments.csv"),
          row.names = F)


