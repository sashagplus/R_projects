library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library("rpart")
library(rpart.plot)

rm(list=ls())

#read data
# rnaseq_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
#                       "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/rnaseq/primary_models/",
#                       "20210315_rnaseq_Bup_7d_premodeling_data.csv")
# 
# rnaseq=read.csv(rnaseq_filename)
# 
# spines_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
#                       "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
#                       "04_Differentiation_Experiments/Data Science/Spines_modeling/",
#                       "20210315_spines_Bup_7d_maxDF_premodeling_data.csv")
# 
# spines=read.csv(spines_filename)
# spines %<>% separate(name, into=c("DF", "Line", "Well", "Days", "fn_Well", "fn_Field"))
# 
# 
# coloc_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
#                      "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
#                      "04_Differentiation_Experiments/Data Science/Colocolization_modeling/",
#                      "20210315_coloc_Bup_7d_maxDF_premodeling_data.csv")
# 
# coloc=read.csv(coloc_filename)
# 
# 
# # coloc %<>% select(FIELD_NO, Line, BUP.Responder, response, mean_prediction)
# # rnaseq %<>% select(Line, response, mean_prediction)
# # spines %<>% select(Line, BUP.Responder, response, mean_prediction)
# # 
#  rnaseq %<>% mutate(Line=as.numeric(gsub('L', '', Line)))
# 
# alllines=full_join(rnaseq %>% select(Line=Line, BUP.Responder=group, RNASeq=Line) %>% distinct ,
#                     full_join(spines %>% select(Line=Line, spines=Line, BUP.Responder) %>% distinct %>% mutate(Line=as.numeric(Line)), 
#                               coloc %>% select(BUP.Responder, Line=Line, coloc=Line) %>% distinct , 
#                               by=c("Line", "BUP.Responder")),
#                     by=c("Line", "BUP.Responder"))
#  
# alllines %<>% left_join(rnaseq, by=c("Line", "BUP.Responder"="group"))  
#  
# alllines %<>% left_join(spines %>% select(-DF, -Well, -Days, -fn_Field, -fn_Well) %>% mutate(Line=as.numeric(Line)), 
#                         by=c("Line", "BUP.Responder")) 
# 
# alllines %<>% left_join(coloc %>% select(-FIELD_NO, -DF) , 
#                         by=c("Line", "BUP.Responder")) 


rnaseq_pred_filename=str_c("//Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                      "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/rnaseq/primary_models/",
                      "20210314_rnaseq_Bup_7d_modeling_predictions.csv")

rnaseq_pred=read.csv(rnaseq_pred_filename)

spines_pred_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                      "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
                      "04_Differentiation_Experiments/Data Science/Spines_modeling/",
                      "20210314_spines_Bup_7d_maxDF_modeling_predictions.csv")

spines_pred=read.csv(spines_pred_filename)

coloc_pred_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                     "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
                     "04_Differentiation_Experiments/Data Science/Colocolization_modeling/",
                     "20210314_coloc_Bup_7d_maxDF_modeling_predictions.csv")

coloc_pred=read.csv(coloc_pred_filename)


coloc_pred %<>% select(FIELD_NO, Line, BUP.Responder, response, mean_prediction)
rnaseq_pred %<>% select(Line, response, mean_prediction)
spines_pred %<>% select(Line, BUP.Responder, response, mean_prediction)

rnaseq_pred %<>% mutate(Line=as.numeric(gsub('L', '', Line)))
# 
# 
# alllines %<>% select(Line, BUP.Responder, RNASeq, spines, coloc)
# 
# 
# 
# 
# alllines %<>% left_join(rnaseq_pred %>% select(rnaseq_pred=mean_prediction, everything()), by=c("Line", BUP.Responder="response" ))
# 
# alllines %<>% left_join(spines_pred %>% select(spines_pred=mean_prediction, Line, response), by=c("Line", BUP.Responder="response" ))
# alllines %<>% left_join(coloc_pred %>% select(coloc_pred=mean_prediction, Line, response), by=c("Line", BUP.Responder="response" ))
alllines=full_join(rnaseq_pred %>% select(Line=Line, response=response, RNASeq=mean_prediction) %>% distinct ,
                   full_join(spines_pred %>% select(Line=Line, spines=mean_prediction, response) %>% distinct %>% mutate(Line=as.numeric(Line)), 
                             coloc_pred %>% select(response, Line=Line, coloc=mean_prediction) %>% distinct , 
                             by=c("Line", "response")),
                   by=c("Line", "response"))




alllines %<>% mutate(RNASeq = round(RNASeq))
alllines %<>% mutate(spines = round(spines))
alllines %<>% mutate(coloc = round(coloc))

alllines %<>% rowwise %<>% mutate(overall_pred=ifelse(all(!is.na(.data$RNASeq), !is.na(.data$spines), !is.na(.data$coloc)), 
                                                      round((.data$RNASeq+ .data$spines+ .data$coloc)/3), 
                                         ifelse(!is.na(.data$RNASeq),.data$RNASeq, .data$coloc )))

print("why?")


accuracies=data.frame()

#for(i in 1:nrow(coloc))
for(i in 1:20)
{ 
  train=c(sample(which(alllines$response==0), length(which(alllines$response==0))*0.5) ,
          sample(which(alllines$response==1), length(which(alllines$response==1))*0.5) )
  test=c(1:nrow(alllines))[-train]
  
  cm.train <- confusionMatrix(factor(alllines$response[train]), 
                              factor(alllines$overall_pred[train]))
  
  cm.test <- confusionMatrix(factor(alllines$response[test]), 
                              factor(alllines$overall_pred[test]))
  
  accuracies %<>% bind_rows(data.frame(nn=i,
                                       train=cm.train$overall["Accuracy"], 
                                       test=cm.test$overall["Accuracy"]))



}


