library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library("rpart")
library(rpart.plot)
library(readxl)    



data_folder=stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Sasha Kugel/",
                           "Clinical Data/ICD10/2021-code-descriptions-tabular-order/")


icd10=read.table(stringr::str_c(data_folder,
                                        "icd10cm_order_2021.txt"), 
                 sep="\t", stringsAsFactors = F, quote = "")


icd10 %<>% separate(V1, into=c("prefix", "long_description"), sep=76)
icd10 %<>% separate(prefix, into=c("prefix", "short_description"), sep=16)
icd10 %<>% separate(prefix, into=c("prefix", "hipaaheader"), sep=14)
icd10 %<>% separate(prefix, into=c("order_number", "icd_10_cm"), sep=6)

icd10 %<>% mutate_all(~ str_trim(., side=c("both")))

write.csv(icd10 %>% filter (str_detect(shortdescription, " depress" ) | str_detect(longdescription, " depress" )),
          stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                         "01_Protocol_Development/01_10_Data_Science/DataforPatientHistory/ICD_andOther_codes/",
                         "icd10_depress_codes.csv"),
          quote = TRUE,
          row.names = FALSE)


