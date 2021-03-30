#Research the current SNPs in streamlit, their source and competability with future
library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)


setwd("GenetikaPlus_Projects_R/genetics/")
data_folder="/Users/sashakugel/data/"


pharma=read.table(stringr::str_c(data_folder,
                                 "StahlsGuide/pharmaGKB/annotations/",
                                 "clinical_ann_metadata.tsv"),
                  sep="\t",
                  stringsAsFactors=FALSE,
                  quote = "",
                  header=TRUE
)

pharma %<>% filter(!Level.of.Evidence %in% c("3","4"))

pharma %>% filter(str_detect(Related.Diseases,regex("Depress", ignore_case = TRUE) )) %>% View


