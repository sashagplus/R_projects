#cross processed lines from claudia with all lines and treatments in G+
library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)

rm(list=ls())

claud=gpio::read_excel_allsheets("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/2021_05_31_Table of analysed lines_Claudia.xlsx")$Data;

drug_abbrvs=data.frame(drug_name=c("Bupropion", "Citalopram (Celexa)", "Mirtazapine", "Nortriptyline", "Sertraline", "Tranylclypromine", "Venlafaxine"),
                       drug_abbrv=c("bup", "cit", "mirt", "ntp", "sert", "trp", "ven"))

status=claud[3:nrow(claud), 1:13];
colnames(status)=claud[2,1:13];
status %<>% gather(key="analysis", value="st", -Line) %>%
            filter(!is.na(st)) %>%
            mutate(gp_id=as.numeric(str_replace(str_replace(Line, "\\*", ""), "LCL-", ""))) %>%
            mutate(treatment=str_replace(str_replace(analysis, "rna_", ""), "img_", "")) %>%
            mutate(analysis_type=str_replace(str_replace(analysis, treatment, ""), "_", "")) ;

status %<>% select(-Line, -analysis) %>% spread(key=analysis_type, value=st, fill="")

resp_sql=str_c("select r2.id_gp_numeric, r2.patient_id , r2.treatment_single_drug, r2.response_imp50p, r2.treatment_level
                  from gpd.responsivness r2 ")

resp=gpio::dbGetQuery(resp_sql);
resp %<>% left_join(drug_abbrvs, by=c("treatment_single_drug"="drug_name"));
resp %<>% left_join(status, by=c("id_gp_numeric"="gp_id", "drug_abbrv"="treatment"))

resp %<>% select(-drug_abbrv) %>%
          rename(treatment=treatment_single_drug, stard_patient_id=patient_id, gp_patient_id=id_gp_numeric)

resp %<>% arrange(gp_patient_id, treatment_level, treatment)

resp %<>% relocate(treatment_level, .after=stard_patient_id)

colnames(resp)=str_c(str_c("a", seq(1:ncol(resp)),"_"), colnames(resp))
#caution!
# gpdb::db.management.addTableFromDataFrame(df = resp,
#                                           newTableName = "lines_status_vs_treatments_20210706",
#                                           schemaName = "pgp")









