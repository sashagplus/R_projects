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

responses.full=responses

responses.byWeek=read.csv(str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/20210512_treatment_qids_levels_selecLatest.csv"))
responses.byWeek %<>% select(-Treatment.level, -Week.in.level, Patient.ID, Treat=Study.med.code.1, 
                             qids.latest.imp=QIDS.C.percent.improvement..transcribed.,
                             qids.latest.beginning=QIDS.C.beginning.score..transcribed.,
                             qids.latest.current=QIDS.C.current.score..transcribed.)
#responses.byWeek %<>% gather(key="Study.med.code.1", value = "QIDS.C.percent.improvement..transcribed.", -Patient.ID, na.rm = T)

responses.byWeek %<>% mutate(Treat.abbrv=case_when(Treat=="Bupropion" ~ "BUP",
                                            Treat=="Citalopram (Celexa)" ~ "CIT",
                                            Treat=="Mirtazapine" ~ "MIRT",
                                            Treat=="Nortriptyline" ~ "NTP",
                                            Treat=="Sertraline" ~ "STL",
                                            Treat=="Tranylclypromine" ~ "TCL",
                                            Treat=="Venlafaxine" ~ "VLF"))


responses.byWeek %<>% select(-Treat)

responses.byWeek.full=responses.byWeek


 print(1)                                        

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

#3a read imaging data
imaging_data=read.csv(str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                                "Genetika+SharedDrive/01_Protocol_Development/",
                                "01_05_Imaging/01_Experiments_and_Results/04_Differentiation_Experiments/",
                                "Line_summary_treatments.csv"))

imaging_data %<>% select(Line=L, everything()) %>% select(-DF, -DF.1, -L.1) 
imaging_data %<>% gather(key="Treat", value="response.sari", -Line ) %<>%
                  mutate(response.sari=ifelse(response.sari=="", NA, response.sari)) %<>%
                  drop_na(response.sari)

imaging_data %<>% mutate(Treat.abbrv=case_when(Treat=="BUP_RESPONSE" ~ "BUP",
                                               Treat=="CIT_RESPONSE" ~ "CIT",
                                               Treat=="MIRT_RESPONSE" ~ "MIRT",
                                               Treat=="NTP_RESPONSE" ~ "NTP",
                                               Treat=="vehicle" ~ "VEH"))


imaging_data %<>% filter(Treat.abbrv=="CIT")
#4. Cross all data
patient_to_sample %<>%  select(id.stard, id.gp.numeric) %<>%
                        left_join(responses, by=c("id.stard"="Patient.ID"))

  patient_to_sample %<>%  left_join(responses.byWeek %>% filter(Treat.abbrv=="CIT"), by=c("id.stard"="Patient.ID", "Treat.abbrv"))

patient_to_sample %<>% mutate(response=ifelse(qids >=50, "R", "NR"), 
                              response.byMaxWeek=ifelse(qids.latest.imp >=50 & qids.latest.current<=6, "R", "NR"))
patient_to_sample %<>% mutate(borderline=ifelse(qids>=40 & qids <=60, TRUE, FALSE),
                              borderline.byMaxWeek=ifelse(qids.latest.imp>=40 & qids <=60, TRUE, FALSE))
patient_to_sample %<>% filter(Treat.abbrv=="CIT")

patient_to_sample %<>% select(id.stard, 
                             id.gp.numeric, 
                             Treat.abbrv.auto=Treat.abbrv, 
                             qids.maxImpInLevel.auto=qids, 
                             response.maxImpInLevel.auto=response, 
                             qids.latest.imp.auto=qids.latest.imp,
                             qids.latest.beginning.auto=qids.latest.beginning,
                             qids.latest.current.auto=qids.latest.current, 
                             REMISSION.latest.auto=response.byMaxWeek, 
                             borderline.auto=borderline) 



#Read Claudias classification:
mat=gpio::read_excel_allsheets("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/Exploratory_Outputs/responsivness/2021_05_12_responsivness_CIT_compare_Claudia_and_PreviousClassification.xlsx")

mat=mat$Sheet1

mat %<>% mutate(Line=as.numeric(str_remove(Line, "L")))
mat %<>% mutate(response.claudia=case_when(`Clauda response`=="BLNR" ~ "NR",
                                           `Clauda response`=="BR" ~ "R",
                                           TRUE ~ `Clauda response`))
mat %<>% mutate(borderline.claudia=case_when(`Clauda response`=="BLNR" ~ TRUE,
                                           `Clauda response`=="BR" ~ TRUE,
                                           TRUE ~ FALSE))

mat %<>% select(Line, raw.claudia=`Clauda response`,  qids.claudia=`Claudia qids`, response.claudia, borderline.claudia)


rnaseq_data %<>% filter(Treatment.abbrv=="CIT") %>% 
  select(CellLIne.numeric, Treatment.abbrv, Remission) %>% 
  distinct() %>% 
  mutate(Remission.Yishai=ifelse(Remission=="NRES", "NR", "R")) %>%
  select(-Remission)


  
mat %<>% left_join(patient_to_sample , by=c("Line"="id.gp.numeric"))
mat %<>% left_join(rnaseq_data , by=c("Line"="CellLIne.numeric"))
mat %<>% left_join(imaging_data %>% select(Line, response.sari) %>% distinct , by=c("Line"))


#!! Sari has sample 10 twice, once NR and once BL NR
mat %<>% filter(!(Line==10 & response.sari=="NR"))

mat %<>% distinct()

mat %<>% mutate(response.sari.binary=ifelse(response.sari=="R", "R", stringr::str_extract(response.sari, "NR")))

library(flextable)

flxtbl2x2=function(df)
{
  #creates formatted table of 2X2 cross table
  rname=colnames(df)[1]
  cname=colnames(df)[2]
  df %<>% spread(key=2, value=3)
  flextable::flextable(df %>% mutate(header=stringr::str_to_title(stringr::str_replace_all(rname, "\\.", "\n")), .before=!!rname) %>% 
                         rename(val=!!rname)) %>% 
    merge_v(~header) %>%
    set_header_labels(header="", val="", NR="NR", R="R") %>%
    add_header_row(values=c("", stringr::str_to_title(stringr::str_replace_all(cname, "\\.", "\n"))), colwidths = c(2,2)) %>%
    align( i = NULL, j = NULL, align = "center", part="all") %>%
    print
}



flxtbl2x2(df=mat %>% select(response.claudia, response.maxImpInLevel.auto) %>% table() %>% as.data.frame())


flxtbl2x2(df=mat %>% select(response.claudia, REMISSION.latest.auto) %>% table() %>% as.data.frame())

flxtbl2x2(df=mat %>% select(response.claudia, Remission.Yishai) %>% table() %>% as.data.frame())

flxtbl2x2(df=mat %>% select(Remission.Yishai, response.maxImpInLevel.auto) %>% table() %>% as.data.frame())

flxtbl2x2(df=mat %>% select(Remission.Yishai, REMISSION.latest.auto) %>% table() %>% as.data.frame())

flxtbl2x2(df=mat %>% select(Remission.Yishai, REMISSION.latest.auto) %>% table() %>% as.data.frame())

flxtbl2x2(df=mat %>% select(Remission.Yishai, response.sari.binary) %>% table() %>% as.data.frame())

flxtbl2x2(df=mat %>% select(REMISSION.latest.auto, Remission.Yishai) %>% table() %>% as.data.frame())
flxtbl2x2(df=mat %>% select(REMISSION.latest.auto, response.sari.binary) %>% table() %>% as.data.frame())

# write.csv(mat, 
#           "/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/Exploratory_Outputs/2021_05_13_comparing_response_calssifications.csv",
#           row.names = FALSE)




#explore R/NR ratios
responses.full %<>% distinct()

responses.full %<>% mutate(response=ifelse(qids>=50, "R", "NR"))
responses.byWeek.full %<>% mutate(response=ifelse(qids.latest.imp>=50, "R", "NR"))

flextable::flextable(as.data.frame( responses.byWeek.full %>% select(Treat.abbrv, response) %>% table) %>% spread(key=response, value=Freq)) 
flextable::flextable(as.data.frame( responses.full %>% select(Treat.abbrv, response) %>% table) %>% spread(key=response, value=Freq)) 


flextable::flextable(mat %>% filter(id.stard==3586)) %>% 
  delete_part( part = "header") %>% 
  add_header_row(values=stringr::str_to_title( stringr::str_replace_all(colnames(mat), "\\.", "\n"))) %>% 
  align( i = NULL, j = NULL, align = "center", part="all") %>% 
  hline( i = NULL, j = NULL, border = NULL, part = "header") %>% 
  color(i=1, j=c(2,3,10, 12, 13, 16,18), color = "red")


flextable::flextable(mat %>% filter(id.stard==1804)) %>% 
  delete_part( part = "header") %>% 
  add_header_row(values=stringr::str_to_title( stringr::str_replace_all(colnames(mat), "\\.", "\n"))) %>% 
  align( i = NULL, j = NULL, align = "center", part="all") %>% 
  hline( i = NULL, j = NULL, border = NULL, part = "header") %>% 
  color(i=1, j=c(2,3,10, 12, 13, 16,18), color = "red")


