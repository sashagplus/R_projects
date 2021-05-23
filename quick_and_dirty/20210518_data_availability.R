library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library("rpart")
library(rpart.plot)
library(readxl)    
library(caret)

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
#Temporarly - reading the file Sari prepared. Later will be the file I curate from Claudia's extraction
responses=read.csv(str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                          "Genetika+SharedDrive/01_Protocol_Development/",
                          "01_10_Data_Science/20210509_treatment_qids_levels.csv"))
responses %<>% select(-Treatment.level, -Week.in.level)
responses %<>% gather(key="Treat", value = "qids", -Patient.ID, na.rm = T)

responses %<>% mutate(Treat.abbrv=case_when(Treat=="Bupropion" ~ "BUP",
                                            Treat=="Citalopram..Celexa." ~ "CIT",
                                            Treat=="Mirtazapine" ~ "MIRT",
                                            Treat=="Nortriptyline" ~ "NTP",
                                            Treat=="Sertraline" ~ "STL",
                                            Treat=="Tranylclypromine" ~ "TCL",
                                            Treat=="Venlafaxine" ~ "VLF"))

responses %<>% mutate(response = ifelse(qids>=50, "R", "NR"))

responses.byWeek=read.csv(str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/20210512_treatment_qids_levels_selecLatest.csv"))
responses.byWeek %<>% select(-Treatment.level, -Week.in.level, Patient.ID, Treat=Study.med.code.1, qids.byMaxWeek=QIDS.C.percent.improvement..transcribed.)
#responses.byWeek %<>% gather(key="Study.med.code.1", value = "QIDS.C.percent.improvement..transcribed.", -Patient.ID, na.rm = T)

responses.byWeek %<>% mutate(Treat.abbrv=case_when(Treat=="Bupropion" ~ "BUP",
                                            Treat=="Citalopram (Celexa)" ~ "CIT",
                                            Treat=="Mirtazapine" ~ "MIRT",
                                            Treat=="Nortriptyline" ~ "NTP",
                                            Treat=="Sertraline" ~ "STL",
                                            Treat=="Tranylclypromine" ~ "TCL",
                                            Treat=="Venlafaxine" ~ "VLF"))


responses.byWeek %<>% select(-Treat)
responses.byWeek %<>% mutate(response = ifelse(qids.byMaxWeek>=50, "R", "NR"))

 print(1)                                        

 
responsSelected=responses 
 
responsSelected %<>% select(Patient.ID, Treat.abbrv, response) %<>% distinct

#3. Read rnaseq data and summarise
rnaseq_data=read.csv(stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                                      "Genetika+SharedDrive/01_Protocol_Development/",
                                      "01_06_RNA-seq/Analysis/raw_data/database/",
                                      "rnaseq_sample_log_analysis.csv"))
rnaseq_data %<>% mutate(CellLIne.numeric=as.numeric(str_remove(CellLIne, "L")))

rnaseq_data %<>% mutate(Treatment.abbrv=case_when(Treatment=="bupropion" ~ "BUP",
                                                  Treatment=="citalopram" ~ "CIT",
                                                  Treatment=="mirtazapine" ~ "MIRT",
                                                  Treatment=="nortriptyline" ~ "NTP",
                                                  Treatment=="vehicle" ~ "VEH"))

rnaseq_data.LineTreat= rnaseq_data %>% select(CellLIne.numeric, Treatment.abbrv) %>% distinct() %>% mutate(rnaseq=1)

#rnaseq_data.LineTreat %<>% filter(Treatment.abbrv!="VEH") %>% left_join(patient_to_sample %>% select(id.stard, id.gp.numeric), by=c("CellLIne.numeric"="id.gp.numeric")) 
#4. imaging data 

imaging_data=read.csv(stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                                     "01_Protocol_Development/01_05_Imaging/", 
                                     "01_Experiments_and_Results/04_Differentiation_Experiments/",
                                     "20210505_File_Location_well_index_per_line_including2019Data.csv"))

imaging_data %<>% filter(Line != 19)
colnames(imaging_data) = c( "DF", "Line", "Well", "Treat", "Days", "BUP.Responder",
                            "folder.spines", "folder.coloc", 
                            "file.coloc.roi.xlsx", "file.coloc.psd.csv", 
                            "file.coloc.syn.csv", "pixel.to.um.ratio")

imaging_data %<>% filter(Days >= 7)

imaging_data.LineTreat=imaging_data %>% select(Line, Treat) %>% distinct() %>% mutate(imaging=1)

#imaging_data.LineTreat %<>% filter(Treat!="VEH") %>% left_join(patient_to_sample %>% select(id.stard, id.gp.numeric), by=c("Line"="id.gp.numeric")) 


pts=patient_to_sample %>% select(id.stard, id.gp.numeric) %>% 
  left_join(responsSelected, by=c("id.stard"="Patient.ID")) %>% 
  left_join(rnaseq_data.LineTreat, by=c(c("id.gp.numeric"="CellLIne.numeric", "Treat.abbrv"="Treatment.abbrv"))) %>% 
  left_join(imaging_data.LineTreat, by=c("id.gp.numeric"="Line", "Treat.abbrv"="Treat")) 



pts %<>% gather(key="AnalysisType", value="cnt", -id.stard, -id.gp.numeric, -Treat.abbrv, -response) %>% drop_na(cnt) 

data_avaialble=pts %>% group_by(AnalysisType, Treat.abbrv, response) %>% summarise(sumcnt=sum(cnt), .groups = "keep") %>% spread(key=Treat.abbrv, value=sumcnt)


library(flextable)

#explore R/NR ratios

flextable::flextable(data_avaialble) %>% 
  merge_v(~AnalysisType) %>%
  set_header_labels(AnalysisType="Analysis Type", response="R / NR") %>%
  align( i = NULL, j = NULL, align = "center", part="all") 





