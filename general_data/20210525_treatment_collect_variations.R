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
#alllevels %<>% filter(Patient.ID %in% patient_to_sample$id.stard)


alllevels %<>% mutate_at(vars(starts_with("Study.med")), 
                         ~ ifelse(.x %in% c("-2", "-3", "-5", "0"), NA, .x) )
#Explore treatments:
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


#Mark the last week of each level
alllevels %<>% group_by(Patient.ID, Treatment.level, treatment_string) %>% mutate(max_week=max(Week.in.level))

#In patients selected by GP, single drug treatments didn't have more then 2 recorded drugs per level.
#alllevels %>% group_by(Patient.ID, Treatment.level) %>% summarise(meds=n_distinct(Study.med.code.1)) %>% View

#Filter last week
alllevels %<>% filter(Week.in.level != 0) %>% filter((!Week.in.level %in% c(0.1, 2)) | max_week==Week.in.level) 
#alllevels %<>% select(-max_week)

#' required output:
#' 1. Patient_ID
#' 2. Treatment Level
#' 3. Treatment String
#' 3. latest_week_in_level
#' 4. qids_start
#' 5. qids_current_week
#' 6. qids_improvement
#' 7. response_imp50p
#' 8. response_imp50p_bl10p
#' 9. response_imp50p_bl20p
#' 10. remission_p50imp_below6current


treatment_summary= alllevels %>% group_by(Patient.ID,Treatment.level,  treatment_string ) %>% 
                            filter(Week.in.level == max(Week.in.level)) %>% 
                            slice(which.max(QIDS.C.percent.improvement..transcribed.)) #==max(QIDS.C.percent.improvement..transcribed.)) %>%
                            #distinct()

treatment_summary %<>% select(Patient.ID,
                         Treatment.level,
                         Treatment=treatment_string,
                         latest_week_in_level=max_week,
                         qids_start=QIDS.C.beginning.score..transcribed. ,
                         qids_current=QIDS.C.current.score..transcribed. ,
                         qids_improvment=QIDS.C.percent.improvement..transcribed.)


treatment_summary %<>% mutate(response_imp50p=ifelse(qids_improvment>=50, "RESP", "NONRESP"),
                         remission_p50imp_below6current=ifelse(qids_improvment>=50 & qids_current<=6, "REMS", "NONREMS")) %<>% 
                  mutate(response_imp50p_bl10p=ifelse(qids_improvment<55 & qids_improvment>45 , 
                                                      stringr::str_c(response_imp50p, "BL"), response_imp50p),
                         response_imp50p_bl20p=ifelse(qids_improvment<60 & qids_improvment>40 , 
                                                      stringr::str_c(response_imp50p, "BL"), response_imp50p))


treatment_summary %<>% left_join(patient_to_sample %>% select(id.stard, id.gp.numeric, id.ruidCellLine ), 
                            by=c("Patient.ID"="id.stard")) %<>% 
                  relocate(id.gp.numeric, gp.id.ruidCellLine=id.ruidCellLine)


celllines_path=stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/", 
                                "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/",
                                "DataforPatientHistory/DepressionDataNIH/StarD/",
                                "STAR_DClinicalData/ASCII Files for viewing/")

celllines_filenames=data.frame(files=c("Level1/cellline.csv",
                      "Level2/cellline-.csv",
                      "Level2a/cellline.csv",
                      "Level3/cellline.csv",
                      "Level4/cellline.csv"), 
                      Treatment.level=c("Level 1",
                                        "Level 2", 
                                        "Level 2A", 
                                        "Level 3", 
                                        "Level 4"))
  

celllines=data.frame()
for(cc in 1:nrow(celllines_filenames))
{
  cc.file=read.csv(stringr::str_c(celllines_path, celllines_filenames$files[cc]))
  cc.file %<>% mutate(Treatment.level=celllines_filenames$Treatment.level[cc])
  celllines %<>% bind_rows(cc.file)
}
  
treatment_summary %<>% left_join(celllines, by=c("Patient.ID", "Treatment.level")) %<>% 
                        relocate("cell.line", .after = gp.id.ruidCellLine) %<>%
                        rename(stard.rui.cell.line=cell.line) %<>%  
                        mutate(is.rui.identical=ifelse(!is.na(gp.id.ruidCellLine),gp.id.ruidCellLine== stard.rui.cell.line, NA)) %<>%
                         relocate("is.rui.identical", .after = gp.id.ruidCellLine) 

treatment_summary %<>% rename(stard.Patient.ID=Patient.ID)

# write.csv(treatment_summary,
#           stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
#                           "Genetika+SharedDrive/01_Protocol_Development/",
#                           "01_10_Data_Science/Exploratory_Outputs/responsivness/",
#                           "2021_05_25_responsiveness_labeling_qids_based_latest_allSTARD_includeCelllines.csv"),
#           row.names = F)

# #May 30 st
# #Nadav asked for the following:
# # Patients that satisfy :
# # 1. Responders to CIT+BUP (level 2)
# # 2. Did not move to level 3
# # 3. Were definite NRs in level 1 (i.e., non-borderliner nonresponders)
# # 4. Treatment lasted at least 9 weeks or more.
# level1_nonresp_ids=treatment_summary %>% filter(Treatment.level=="Level 1" & 
#                                                   qids_improvment<40 & 
#                                                   Treatment!="") %>% 
#                                           pull(stard.Patient.ID) %>% unique()
# 
# desired_ids=treatment_summary %>%  group_by(stard.Patient.ID) %>% 
#                                     mutate(max_level=max(Treatment.level)) %>% 
#                                     filter(Treatment=="Bupropion, Citalopram (Celexa)" & 
#                                         latest_week_in_level>=9 & 
#                                         max_level<"Level 3" &
#                                         stard.Patient.ID %in% level1_nonresp_ids) %>% 
#                                       pull(stard.Patient.ID)



# #Nadav asked for the following:
# # table with  :
# # 1. stard patient id
# # 2. genetika plus patient id
# # 3. cellline rui id
# # 4. For every level 1, 2, 2A, 3, 4
# #  a. "Treatment.level"
#     # b. "latest_week_in_level"
#     # c. "Treatment"
#     # d. "response_imp50p"
#     # e. "remission_p50imp_below6current"
# 
# ts=treatment_summary %>% ungroup %>% 
#   #filter(!is.na(id.gp.numeric)) %>% 
#   select(stard.Patient.ID, id.gp.numeric,gp.id.ruidCellLine, Treatment.level, Treatment, latest_week_in_level, response_imp50p, remission_p50imp_below6current) %>% 
#   distinct() %>% 
#   ungroup %>% 
#   group_by(stard.Patient.ID, id.gp.numeric, gp.id.ruidCellLine , Treatment.level) %>% 
#   mutate(gid=row_number()) %>% 
#   relocate(gid, .after = id.gp.numeric)
# 
# ts %<>%  mutate(t1=Treatment.level) %<>% 
#   pivot_wider(names_from = t1, values_from =c("Treatment.level", "Treatment", "latest_week_in_level", "response_imp50p", "remission_p50imp_below6current") , names_sep="__") 
# 
# order_of_columns=c("Treatment.level", "latest_week_in_level", "Treatment", "response_imp50p", "remission_p50imp_below6current")
# 
# names_of_columns=data.frame(org=colnames(ts)[-c(1:4)] ) %>% 
#   separate(org, into=c("first", "second"), sep="__", remove = F) %>% 
#   rowwise %>% 
#   mutate(newname=str_c(second, "__", first)) %>% 
#   arrange(second, match(first, order_of_columns))
# 
# ts %<>% select(stard.Patient.ID, id.gp.numeric, gp.id.ruidCellLine, names_of_columns$org) %>% 
#   arrange(stard.Patient.ID) 
# 
# write.csv(ts,
#           stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
#                           "Genetika+SharedDrive/01_Protocol_Development/",
#                           "01_10_Data_Science/Exploratory_Outputs/responsivness/",
#                           "2021_05_30_responsiveness_labelingByColumn_GPpatients.csv"),
#           row.names = F)
# write.csv(ts,
#           stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
#                          "Genetika+SharedDrive/01_Protocol_Development/",
#                          "01_10_Data_Science/Exploratory_Outputs/responsivness/",
#                          "2021_05_31_responsiveness_labelingByColumn_ALLpatients.csv"),
#           row.names = F)















