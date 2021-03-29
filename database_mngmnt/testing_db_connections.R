#' creator name: Sasha 
#' created date: March 21, 2021
#' Objective: Test data upload and download from postgresql 
#' Input: none
#' Output: noneß 
#' 
library(RPostgreSQL)
library(DBI)
rm(list=ls())
patient_to_sample=read.table(str_c("/Users/sashakugel/gplus_dropbox/",
                                  "Genetika+ Dropbox/Genetika+SharedDrive/",
                                  "01_Protocol_Development/01_10_Data_Science/",
                                  "2020_GP_Patient_Codes-StarD.tsv"),
                             header=TRUE,
                             stringsAsFactors = F,
                             sep="\t")




patient_to_sample %<>% mutate(GP...CODE.clean=str_remove(GP...CODE, "\\*"))

patient_to_sample %<>% mutate(GP...CODE.numeric=as.numeric(str_remove(GP...CODE.clean, "LCL-")))

db = "postgres"

host_db = "localhost" #i.e. # i.e. 'ec2-54-83-201-96.compute-1.amazonaws.com'  

db_port <- "5432"  # or any other port specified by the DBA

db_user <- "postgres"  

db_password <- "genetikapostgre";

con <- dbConnect(RPostgreSQL::PostgreSQL(), 
                 dbname = db, 
                 host=host_db, 
                 port=db_port, 
                 user=db_user, 
                 password=db_password,
                 options="-c search_path=gplus_beta")  


dbWriteTable(con, "patients_stardID_to_gpId", patient_to_sample)


#dbRemoveTable(con, "patients_stardID_to_gpId")

tester=dbReadTable(con, "patients_stardID_to_gpId")

print(str_c("Are equal? ", all_equal(patient_to_sample, tester)))
#data frames maybe of different formats


dbDisconnect(con) 
