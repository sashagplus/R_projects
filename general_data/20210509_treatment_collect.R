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

# read_excel_allsheets <- function(filename, tibble = FALSE) {
#   # I prefer straight data.frames
#   # but if you like tidyverse tibbles (the default with read_excel)
#   # then just pass tibble = TRUE
#   sheets <- readxl::excel_sheets(filename)
#   x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
#   if(!tibble) x <- lapply(x, as.data.frame)
#   names(x) <- sheets
#   x
# }

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
  
  #levelx %<>% select(Patient.ID, Treatment.level, Week.in.level, starts_with("Study.med.code."), QIDS.C.percent.improvement..transcribed.) 
  
  levelx %<>% 
    mutate_at(vars(starts_with("Study.med.code.")), ~ as.character(.x)) %<>%
    mutate_at(vars(starts_with("Prescribed.med.code.")), ~ as.character(.x)) %<>%
  mutate_at(vars(starts_with("Prescribed.med.daily.dose.")), ~ as.character(.x)) %<>%
    mutate_at(vars(starts_with("Med.type")), ~ as.character(.x))
  
  
   return(levelx)
}

#level 1
pathname="/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/DataforPatientHistory/DepressionDataNIH/StarD/STAR_DClinicalData/ASCII Files for viewing/"

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

alllevels=alllevels.full  %>% select(Patient.ID, Treatment.level, Week.in.level, starts_with("Study.med.code."), QIDS.C.percent.improvement..transcribed., QIDS.C.beginning.score..transcribed., QIDS.C.current.score..transcribed.) 
alllevels %<>% filter(Patient.ID %in% patient_to_sample$id.stard)

alllevels %<>% filter(Study.med.code.2==-2) %>% 
                filter(Study.med.code.1!=-2) 

alllevels %<>% ungroup %<>% select(-Study.med.code.2, -Study.med.code.3, -Study.med.code.4) %<>% distinct()

alllevels %<>% group_by(Patient.ID, Treatment.level) %>% mutate(max_week=max(Week.in.level))

alllevels %<>% filter(Week.in.level != 0) %>% filter((!Week.in.level %in% c(0.1, 2)) | max_week==Week.in.level) 

alllevels %<>% select(-max_week)

single_drugs= alllevels %>% group_by_at(vars(Patient.ID, starts_with("Study.med.code."))) %>% 
                            filter(QIDS.C.percent.improvement..transcribed. == max(QIDS.C.percent.improvement..transcribed.)) %>% 
                            #filter(Study.med.code.2==-2) %>% 
                            #filter(Study.med.code.1!=-2) %>% 
                            #select(-Treatment.level, -Week.in.level) %>%
                            distinct()

single_drugs.byWeek= alllevels %>% group_by_at(vars(Patient.ID, starts_with("Study.med.code."))) %>% 
                            filter(Week.in.level == max(Week.in.level)) %>% 
                            filter(QIDS.C.percent.improvement..transcribed.==max(QIDS.C.percent.improvement..transcribed.)) %>%
                            #filter(Study.med.code.2==-2) %>% 
                            #filter(Study.med.code.1!=-2) %>% 
                            #select(-Treatment.level, -Week.in.level) %>%
                            distinct()

#single_drugs %<>% mutate_at(vars(Treatment.level), ~ as.numeric(str_remove(str_remove(.x, "Level "), "A"))) 
#single_drugs %<>% group_by(Patient.ID, Study.med.code.1) %>% 
#  filter(Treatment.level == max(Treatment.level)) 

#single_drugs %<>% ungroup %<>% select(-Study.med.code.2, -Study.med.code.3, -Study.med.code.4, -Week.in.level) %<>% distinct()

#single_drugs %<>% group_by(Patient.ID, Treatment.level, Study.med.code.1) %<>% 
#                  summarise(QIDS.C.percent.improvement..transcribed.=max(QIDS.C.percent.improvement..transcribed., na.rm=TRUE), .groups = "keep")

single_drugs %<>% distinct() %<>% spread(key=Study.med.code.1, value=QIDS.C.percent.improvement..transcribed.)




write.csv(single_drugs, 
          "/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/20210509_treatment_qids_levels.csv",
          row.names = F)

write.csv(single_drugs.byWeek, 
          "/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/20210512_treatment_qids_levels_selecLatest.csv",
          row.names = F)


