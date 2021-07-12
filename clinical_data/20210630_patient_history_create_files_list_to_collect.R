library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)

find_output_file_path="/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/DataforPatientHistory/DepressionDataNIH/StarD/STAR_DClinicalData/GenetikaVersions/20210630_stard_patient_history_ascii_viewing_files_list.txt"

f=read.csv(find_output_file_path, header = F)

colnames(f)="filepath"

f %<>% mutate(path=dirname(filepath),
              filename=basename(filepath))

f %<>% mutate(to = strsplit(filepath, "/")) %>%
  unnest(to) %>%
  group_by(filepath) %>%
  mutate(row = row_number()) %>%
  spread(row, to) 
break()
colnames(f)[2:ncol(f)]=c("initpath","first", "second", "third")

f %<>% select(-initpath)

