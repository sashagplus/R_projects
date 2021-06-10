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
responses=read.csv(stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                                  "Genetika+SharedDrive/01_Protocol_Development/",
                                  "01_10_Data_Science/Exploratory_Outputs/responsivness/",
                                  "2021_05_24_responsiveness_labeling_qids_based_latest.csv"))

responses %<>% mutate(Treat.abbrv=case_when(Treatment_single_drug=="Bupropion" ~ "BUP",
                                            Treatment_single_drug=="Citalopram (Celexa)" ~ "CIT",
                                            Treatment_single_drug=="Mirtazapine" ~ "MIRT",
                                            Treatment_single_drug=="Nortriptyline" ~ "NTP",
                                            Treatment_single_drug=="Sertraline" ~ "STL",
                                            Treatment_single_drug=="Tranylclypromine" ~ "TCL",
                                            Treatment_single_drug=="Venlafaxine" ~ "VLF"))


 
responsSelected=responses 
 
responsSelected %<>% select(Patient.ID, Treat.abbrv, response_imp50p, remission_p50imp_below6current) %<>% distinct

#3. Read rnaseq data and summarise
rnaseq_data=read.csv(stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                                      "Genetika+SharedDrive/01_Protocol_Development/",
                                      "01_06_RNA-seq/Analysis/raw_data/database/",
                                      "rnaseq_sample_log.csv"))
rnaseq_data %<>% mutate(CellLIne.numeric=as.numeric(str_remove(CellLIne, "L")))

rnaseq_data %<>% mutate(Treatment.abbrv=case_when(Treatment=="bupropion" ~ "BUP",
                                                  Treatment=="citalopram" ~ "CIT",
                                                  Treatment=="mirtazapine" ~ "MIRT",
                                                  Treatment=="nortriptyline" ~ "NTP",
                                                  Treatment=="vehicle" ~ "VEH"))

rnaseq_data.LineTreat= rnaseq_data %>% select(CellLIne.numeric, Treatment.abbrv) %>% distinct() %>% mutate(rnaseq=1)

#rnaseq_data.LineTreat %<>% filter(Treatment.abbrv!="VEH") %>% left_join(patient_to_sample %>% select(id.stard, id.gp.numeric), by=c("CellLIne.numeric"="id.gp.numeric")) 
#4. imaging data 

imaging_data=gpio::read_excel_allsheets(stringr::str_c("/Users/sashakugel/gplus_dropbox/",
                                      "Genetika+ Dropbox/Genetika+SharedDrive/",
                                      "01_Protocol_Development/01_10_Data_Science/",
                                      "DATABASE/imaging/",
                                      "imaging_data_file_map_all_data_updated_weekly.xlsx"))
imaging_data = imaging_data$Sheet1

imaging_data %<>% filter(Line != 19) #not enough neurons; low fields count
colnames(imaging_data) = c( "DF", "Line", "Well", "Treat", "Days", # "BUP.Responder",
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

pts %<>% mutate(rnaseq=ifelse(rnaseq==1, "rnaseq", ""),
                imaging=ifelse(imaging==1, "imaging", "")) 

pts %<>% filter(!is.na(rnaseq) | !is.na(imaging))

totals_by_tretment_and_analysis=pts %>% group_by(rnaseq, imaging, Treat.abbrv) %>% 
  tally() %>% 
  spread(key=Treat.abbrv, value=n) %>% 
  ungroup %>% bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))


#pts %<>% gather(key="AnalysisType", value="cnt", -id.stard, -id.gp.numeric, -Treat.abbrv, -response) %>% drop_na(cnt) 

#data_avaialble=pts %>% group_by(AnalysisType, Treat.abbrv, response) %>% summarise(sumcnt=sum(cnt), .groups = "keep") %>% spread(key=Treat.abbrv, value=sumcnt)

data_avaialble_response=pts %>% group_by(rnaseq, imaging, Treat.abbrv, response_imp50p) %>% 
              tally() %>% 
              mutate(analysis_type=ifelse((!is.na(rnaseq) & !is.na(imaging)), "both", ifelse(!is.na(rnaseq), rnaseq, imaging))) %>% 
              ungroup %>% 
              select(-rnaseq, -imaging) %>% 
              relocate(analysis_type) %>% 
              mutate(treat_response=stringr::str_c(Treat.abbrv, "_", response_imp50p), .keep = c( "unused") ) %>% 
              spread(key=treat_response, value=n)


data_avaialble_remission=pts %>% group_by(rnaseq, imaging, Treat.abbrv, remission_p50imp_below6current) %>% 
  tally() %>% 
  mutate(analysis_type=ifelse((!is.na(rnaseq) & !is.na(imaging)), "both", ifelse(!is.na(rnaseq), rnaseq, imaging))) %>% 
  ungroup %>% 
  select(-rnaseq, -imaging) %>% 
  relocate(analysis_type) %>% 
  mutate(treat_response=stringr::str_c(Treat.abbrv, "_", remission_p50imp_below6current), .keep = c( "unused") ) %>% 
  spread(key=treat_response, value=n)

library(flextable)

#explore R/NR ratios

# flextable::flextable(data_avaialble) %>% 
#   merge_v(~analysis_type) %>%
#   set_header_labels(analysis_type="Analysis Type", response="R / NR") %>%
#   align( i = NULL, j = NULL, align = "center", part="all") 

flextable::flextable(data_avaialble_response) %>% 
  merge_v(~analysis_type) %>%
  set_header_labels(analysis_type="Analysis Type", BUP_NONRESP="NR", BUP_RESP="R", CIT_NONRESP="NR", CIT_RESP="R", MIRT_NONRESP="NR",MIRT_RESP="R", NTP_NONRESP="NR", NTP_RESP="R") %>%
  add_header_row(values=c("", "BUP", "CIT", "MIRT", "NTP"), colwidths = c(1, 2,2, 2,2)) %>%
  add_header_row(values=c("RESPONSE"), colwidths = c(9)) %>%
  align( i = NULL, j = NULL, align = "center", part="all") %>% vline(j = c(1,3, 5,7,9), border = officer::fp_border(color="black", width = 1))


flextable::flextable(data_avaialble_remission) %>% 
  merge_v(~analysis_type) %>%
  set_header_labels(analysis_type="Analysis Type", BUP_NONREMS="NR", BUP_REMS="R", CIT_NONREMS="NR", CIT_REMS="R", MIRT_NONREMS="NR",MIRT_REMS="R", NTP_NONREMS="NR", NTP_REMS="R") %>%
  add_header_row(values=c("", "BUP", "CIT", "MIRT", "NTP"), colwidths = c(1, 2,2, 2,2)) %>%
  add_header_row(values=c("REMISSION"), colwidths = c(9)) %>%
  align( i = NULL, j = NULL, align = "center", part="all") %>% vline(j = c(1,3, 5,7,9), border = officer::fp_border(color="black", width = 1))



remission_data=data_avaialble_remission %>% tibble::column_to_rownames("analysis_type") %>% 
                                            t %>% 
                                            as.data.frame() %>% 
                                            tibble::rownames_to_column("temp") %>%  
                                            tidyr::separate(col="temp", into=c("Treatment", "Remission"), sep="_") %>% 
                                            select(Treatment, Remission, remission_congruent=both, remission_rnaseq=rnaseq,  remission_imaging=imaging ) %>%
  rowwise() %>%
                                            mutate(remission_Total=sum(remission_congruent, remission_rnaseq, remission_imaging, na.rm = T))

response_data=data_avaialble_response %>% tibble::column_to_rownames("analysis_type") %>% 
                                          t %>% 
                                          as.data.frame() %>% 
                                          tibble::rownames_to_column("temp") %>%  
                                          tidyr::separate(col="temp", into=c("Treatment", "Response"), sep="_") %>% 
                                          select(Treatment, Response, response_congruent=both, response_rnaseq=rnaseq,  response_imaging=imaging ) %>%
  rowwise() %>%
mutate(response_Total=sum(response_congruent, response_rnaseq, response_imaging, na.rm = T))

response_data %<>% group_by(Treatment) %>% summarise_at(.vars=vars(response_congruent,response_rnaseq,response_imaging,response_Total),
                                                       ~ sum(.x, na.rm=T)) %>% 
  mutate(Response="Total", .after=Treatment) %>% rowwise %>% mutate(response_rnaseq=response_rnaseq+response_congruent,
                                                                    response_imaging=response_imaging+response_congruent) %>% 
  bind_rows(response_data) %>% arrange(match(Treatment, c("BUP", "MIRT", "NTP", "CIT")), 
                                       match(Response, c("NONRESP", "RESP", "Total")))

remission_data %<>% group_by(Treatment) %>% summarise_at(.vars=vars(remission_congruent,remission_rnaseq,remission_imaging,remission_Total),
                                                        ~ sum(.x, na.rm=T)) %>% 
  mutate(Remission="Total", .after=Treatment) %>% rowwise %>% mutate(remission_rnaseq=remission_rnaseq+remission_congruent,
                                                                     remission_imaging=remission_imaging+remission_congruent) %>% 
  bind_rows(remission_data) %>% arrange(match(Treatment, c("BUP", "MIRT", "NTP", "CIT")), 
                                       match(Remission, c("NONRESP", "RESP", "Total")))


alldata=response_data %>% mutate(temp=ifelse(Response=="RESP", "REMS", ifelse(Response=="NONRESP" ,"NONREMS", "Total"))) %>% 
  right_join(remission_data, by=c("Treatment", "temp"="Remission")) %>% 
  select(-temp)

alldata %<>% arrange(match(Treatment, c("BUP", "MIRT", "NTP", "CIT")), 
                     match(Response, c("NONRESP", "RESP", "Total")))

 # alldata %<>% rowwise() %>% mutate(Total=sum(response_congruent, response_rnaseq, response_imaging, na.rm = T)) %>%
 #             relocate("Total", .before=response_congruent)

ft=
  flextable::flextable(alldata) %>% 
#  merge_v(~Treatment) %>%
  set_header_labels(Treatment="", 
                    Response="Response/ Remission Label", 
                    response_congruent="congruent", 
                    response_rnaseq="RNASeq",  
                    response_imaging="Imaging",
                    response_Total="Total",
                    remission_congruent="congruent", 
                    remission_rnaseq="RNASeq",  
                    remission_imaging="Imaging",
                    remission_Total="Total") %>%
  vline(j = c(5,9), border = officer::fp_border(color="grey", width = 1), part="all") %>%
  add_header_row(values=c("","","Response", "Remission"), colwidths = c(1, 1,4, 4)) %>%
  border_outer(part="all", border = officer::fp_border(color="red", width = 2) ) %>%
  add_header_row(values=c("Genetika+ Data Availability"), colwidths = c(10)) %>%
  align( i = NULL, j = NULL, align = "center", part="all") %>% 
  align( i = NULL, j = 2, align = "left", part="body") %>% 
  vline(j = c(2,6), border = officer::fp_border(color="black", width = 1), part="all") %>%
  hline(i=c(3,6, 9), border=officer::fp_border(color="lightgrey", width = 1, style="solid")) %>% 
  hline(i=c(2,5, 8), j=c(2:10), border=officer::fp_border(color="lightgrey", width = 1, style="dashed")) %>% 
  bold(  j = c(6, 10), bold = TRUE, part = "all") %>% 
  color(i = c(3,6,9), j = c(3,7), color="grey", part = "body", source = j) %>%
  merge_v(~Treatment) %>%
  fix_border_issues()
  
#save_as_pptx(ft, path="/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/Exploratory_Outputs/2021_05_30_data_availability.pptx")
print(ft)


