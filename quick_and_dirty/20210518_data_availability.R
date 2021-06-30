library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library("rpart")
library(rpart.plot)
library(readxl)    
library(caret)
library(flextable)

rm(list=ls())

extract_date=Sys.Date()
#return the closest PREVIOUS weekday to curr_date
#if weekdays(curr_date)=="Monday" then return the date of previous weekday
getDateOfLatestWeekday=function(curr_date=Sys.Date(),
                       weekday="Monday")
{
  d=as.Date(curr_date, "%Y%m%d");
  if(weekdays(d)==weekday)
  {
    d=d-1
  }
  prev.days <- seq(d-7,d,by='day')
  return(prev.days[weekdays(prev.days)==weekday])
}


current_treatments=c("BUP", "MIRT", "CIT", "NTP", "VN")

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
                                  "2021_06_14_responsiveness_labeling_qids_based_latest_allTreatments.csv"))

responses %<>% mutate(Treat.abbrv=case_when(Treatment_single_drug=="Bupropion" ~ "BUP",
                                            Treatment_single_drug=="Citalopram (Celexa)" ~ "CIT",
                                            Treatment_single_drug=="Mirtazapine" ~ "MIRT",
                                            Treatment_single_drug=="Nortriptyline" ~ "NTP",
                                            Treatment_single_drug=="Sertraline" ~ "STL",
                                            Treatment_single_drug=="Tranylclypromine" ~ "TCL",
                                            Treatment_single_drug=="Venlafaxine" ~ "VN"))


 
responsSelected=responses 
 
responsSelected %<>% select(Patient.ID, Treat.abbrv, response_imp50p, remission_p50imp_below6current) %<>% distinct

#3. Read rnaseq data and summarise
rnaseq_data=read.csv(stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                                      "Genetika+SharedDrive/01_Protocol_Development/",
                                      "01_06_RNA-seq/Analysis/raw_data/database/",
                                      "rnaseq_sample_log_dateYYYYMMDD.csv"))
rnaseq_data %<>% mutate(CellLIne.numeric=as.numeric(str_remove(CellLIne, "L")))

rnaseq_data %<>% mutate(Treatment.abbrv=case_when(Treatment=="bupropion" ~ "BUP",
                                                  Treatment=="citalopram" ~ "CIT",
                                                  Treatment=="mirtazapine" ~ "MIRT",
                                                  Treatment=="nortriptyline" ~ "NTP",
                                                  Treatment=="vehicle" ~ "VEH"))
#calculate which are old, and which were added in the latest week
rnaseq_data %<>% dplyr::mutate(date_added=as.Date(as.character(date_added), "%Y%m%d")) %<>% 
                dplyr::mutate(isNew=ifelse(date_added>=getDateOfLatestWeekday(curr_date = extract_date), TRUE, FALSE)) 

rnaseq_data.LineTreat= rnaseq_data %>% select(CellLIne.numeric, Treatment.abbrv, isNew) %>% distinct() %>% mutate(rnaseq=1)

#rnaseq_data.LineTreat %<>% filter(Treatment.abbrv!="VEH") %>% left_join(patient_to_sample %>% select(id.stard, id.gp.numeric), by=c("CellLIne.numeric"="id.gp.numeric")) 
#4. imaging data 

imaging_data=gpio::read_excel_allsheets(stringr::str_c("/Users/sashakugel/gplus_dropbox/",
                                      "Genetika+ Dropbox/Genetika+SharedDrive/",
                                      "01_Protocol_Development/01_10_Data_Science/",
                                      "DATABASE/imaging/",
                                      "imaging_data_file_map_all_data_updated_weekly_dateYYYYMMDD.xlsx"))
imaging_data = imaging_data$Sheet1

imaging_data %<>% filter(Line != 19) #not enough neurons; low fields count
colnames(imaging_data) = c( "date_added",
                            "DF", "Line", "Well", "Treat", "Days", # "BUP.Responder",
                            "folder.spines", "folder.coloc", 
                            "file.coloc.roi.xlsx", "file.coloc.psd.csv", 
                            "file.coloc.syn.csv", "pixel.to.um.ratio")

imaging_data %<>% filter(Days >= 7)

#calculate which are old, and which were added in the latest week
imaging_data %<>% dplyr::mutate(date_added=as.Date(as.character(date_added), "%Y%m%d")) %<>% 
                  dplyr::mutate(isNew=ifelse(date_added>=getDateOfLatestWeekday(curr_date = extract_date), TRUE, FALSE)) 

imaging_data.LineTreat=imaging_data %>% select(Line, Treat, isNew) %>% distinct() %>% mutate(imaging=1)

#imaging_data.LineTreat %<>% filter(Treat!="VEH") %>% left_join(patient_to_sample %>% select(id.stard, id.gp.numeric), by=c("Line"="id.gp.numeric")) 


pts=merge(patient_to_sample %>% select(id.stard, id.gp.numeric), data.frame(Treat.abbrv=current_treatments))
pts %<>% filter(!(id.gp.numeric==92 & Treat.abbrv=="MIRT")) #Why?
pts %<>%  left_join(rnaseq_data.LineTreat, by=c(c("id.gp.numeric"="CellLIne.numeric", "Treat.abbrv"="Treatment.abbrv"))) 

pts %<>%  left_join(imaging_data.LineTreat, by=c("id.gp.numeric"="Line", "Treat.abbrv"="Treat"),
                    suffix = c(".rnaseq", ".imaging")) 

pts %<>% mutate(rnaseq_oon=ifelse(rnaseq==1, ifelse(isNew.rnaseq==TRUE, "rnaseq_new", "rnaseq_old")  , ""),
                imaging_oon=ifelse(imaging==1, ifelse(isNew.imaging==TRUE, "imaging_new", "imaging_old"), "")) %<>%
  dplyr::mutate(rnaseq=ifelse(rnaseq==1, "rnaseq", ""),
                imaging=ifelse(imaging==1, "imaging", ""))
pts %<>% left_join(responsSelected, by=c("id.stard"="Patient.ID", "Treat.abbrv")) 

pts %<>% filter(!is.na(rnaseq) | !is.na(imaging))

#!!!TEMP
#pts[16:21, "rnaseq_oon"]="rnaseq_new"; pts[18:22, "imaging_oon"]="imaging_new"

data_avaialble_response=pts %>% group_by(Treat.abbrv, response_imp50p, rnaseq_oon, imaging_oon) %>% tally() %>%
                                dplyr::mutate_at(.vars=dplyr::vars("rnaseq_oon", "imaging_oon"), ~ifelse(is.na(.x), "", .x)) %>% dplyr::ungroup() %>%
                                dplyr::mutate(analysis_type_str=stringr::str_c(rnaseq_oon, "XXX", imaging_oon), .keep="unused") %>%
                                tidyr::spread(key="analysis_type_str", value=n)

data_avaialble_response %<>% dplyr::mutate(both=rowSums(dplyr::select(dplyr::select(., dplyr::starts_with("rnaseq_")), dplyr::contains("imaging")), na.rm = T)) %>%
                            dplyr::mutate(imaging=rowSums(dplyr::select(., dplyr::starts_with("XXXimaging")), na.rm = T),
                                          rnaseq=rowSums(dplyr::select(., dplyr::ends_with("XXX")), na.rm = T)) 

#calculate totals per treatment and response:
data_avaialble_response %<>% dplyr::mutate(total=rowSums(dplyr::select(., both, imaging, rnaseq), na.rm = T))

treatment.sum=data_avaialble_response %>% dplyr::group_by(Treat.abbrv) %>% 
                                          dplyr::summarise_at(.vars=dplyr::vars(both, imaging, rnaseq, total), ~sum(.x, na.rm=T)) %>% 
                                          dplyr::ungroup()

#data_avaialble_response %<>% dplyr::select(-rnaseq_oldXXXimaging_old,  -rnaseq_oldXXX, -XXXimaging_old)

data_avaialble_response %<>% dplyr::mutate(both_="",
                                           rnaseq_="",
                                           imaging_="")


# define brackets
# 1. both, going over all options ([rna new,img new], [rna old, img new], [rna new, img old] )
if ("rnaseq_newXXXimaging_new" %in% colnames(data_avaialble_response)) data_avaialble_response %<>% 
                                                    dplyr::mutate(both_=ifelse(!is.na(rnaseq_newXXXimaging_new), 
                                                                               stringr::str_c(both_, "+", rnaseq_newXXXimaging_new, " both"), 
                                                                               ""))

if ("rnaseq_oldXXXimaging_new" %in% colnames(data_avaialble_response)) data_avaialble_response %<>% 
                                                      dplyr::mutate(both_=ifelse(!is.na(rnaseq_oldXXXimaging_new), 
                                                                                 stringr::str_c(both_, 
                                                                                                ifelse(both_!="", ",", "")
                                                                                                , "+", 
                                                                                                rnaseq_oldXXXimaging_new, 
                                                                                                " img"), 
                                                                                 both_))

if ("rnaseq_newXXXimaging_old" %in% colnames(data_avaialble_response)) data_avaialble_response %<>% 
                                                      dplyr::mutate(both_=ifelse(!is.na(rnaseq_newXXXimaging_old), 
                                                                                 stringr::str_c(both_, 
                                                                                                ifelse(both_!="", ",", "")
                                                                                                , "+", 
                                                                                                rnaseq_newXXXimaging_old, 
                                                                                                " rna"), 
                                                                                 both_))

# 2. rna
if ("rnaseq_newXXX" %in% colnames(data_avaialble_response)) data_avaialble_response %<>% 
                                                dplyr::mutate(rnaseq_=ifelse(!is.na(rnaseq_newXXX), 
                                                                           stringr::str_c(rnaseq_, 
                                                                                          ifelse(rnaseq_!="", ",", "")
                                                                                          , "+", 
                                                                                          rnaseq_newXXX, 
                                                                                          " rna"), 
                                                                           rnaseq_))

# 3. imaging
if ("XXXimaging_new" %in% colnames(data_avaialble_response)) data_avaialble_response %<>% 
                                                dplyr::mutate(imaging_=ifelse(!is.na(XXXimaging_new), 
                                                                             stringr::str_c(imaging_, 
                                                                                            ifelse(imaging_!="", ",", "")
                                                                                            , "+", 
                                                                                            XXXimaging_new, 
                                                                                            " img"), 
                                                                             imaging_))


data_avaialble_response %<>% dplyr::select(-dplyr::contains("XXX"))

data_avaialble_response %<>% dplyr::mutate(both_str=as.character(ifelse(both_!="", stringr::str_c(both, "(", both_, ")"), both)),
                                           imaging_str=as.character(ifelse(imaging_!="", stringr::str_c(imaging, "(", imaging_, ")"), imaging)),
                                           rnaseq_str=as.character(ifelse(rnaseq_!="", stringr::str_c(rnaseq, "(", rnaseq_, ")"), rnaseq)))

data_avaialble_response %<>% dplyr::select(-both, -imaging, -rnaseq, -both_, -rnaseq_, -imaging_)

print(1)

alldata=treatment.sum %>% dplyr::mutate(both_str=as.character(both),imaging_str=as.character(imaging), rnaseq_str=as.character(rnaseq), .keep="unused") %>% 
                          dplyr::mutate(response_imp50p="Total") %>% dplyr::bind_rows(data_avaialble_response) %>% 
                          dplyr::rename(Treatment=Treat.abbrv, Response=response_imp50p) %>%
                          dplyr::relocate(Response, .after = Treatment) 

print(2)
alldata %<>% dplyr::arrange(Treatment, 
                     match(Response, c("NONRESP", "RESP", "Total"))) %>%
              dplyr::relocate(total, .after = dplyr::last_col())


ft_response=
  flextable::flextable(alldata ) %>% 
  #  merge_v(~Treatment) %>%
  set_header_labels(Treatment="Treatment", 
                    Response="Response Label", 
                    both_str="congruent", 
                    rnaseq_str="RNASeq",  
                    imaging_str="Imaging",
                    total="Total") %>%
  vline(j = c(5), border = officer::fp_border(color="grey", width = 1), part="all") %>%
  add_header_row(values=c("","","Response"), colwidths = c(1, 1,4)) %>%
  border_outer(part="all", border = officer::fp_border(color="red", width = 2) ) %>%
  add_header_row(values=c("Genetika+ Data Availability"), colwidths = c(6)) %>%
  align( i = NULL, j = NULL, align = "center", part="all") %>% 
  align( i = NULL, j = 2, align = "left", part="body") %>% 
  vline(j = c(2), border = officer::fp_border(color="black", width = 1), part="all") %>%
  hline(i=c(3,6,9,12), border=officer::fp_border(color="lightgrey", width = 1, style="solid")) %>% 
  hline(i=c(2,5,8,13), j=c(2:6), border=officer::fp_border(color="lightgrey", width = 1, style="dashed")) %>% 
  bold(  j = c(6), bold = TRUE, part = "all") %>% 
  color(i = c(3,6,9,12), j = c(3), color="grey", part = "body", source = j) %>%
  autofit( add_w = 0.1, add_h = 0.1, part = c("body", "header")) %>%
  merge_v(~Treatment) %>%
  fix_border_issues()




print(ft_response)


