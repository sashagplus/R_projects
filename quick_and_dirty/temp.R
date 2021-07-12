library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library("rpart")
library(rpart.plot)
library(readxl)    
library(caret)

rm(list=ls())

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

data_folder=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                  "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/",
                  "04_Differentiation_Experiments/Compiled_Results/Colocalization/20210311_Bup_by_Field_for_Sasha/")


coloc=read_excel_allsheets(str_c(data_folder, "20210311_Bup_by_Field_DF42to74.xlsx"))
#coloc1=read_excel_allsheets(str_c(data_folder, "20210311_Bup_by_Field_DF42to74.xlsx"))

coloc=coloc$Sheet1
#coloc1=coloc1$Sheet1

colnames(coloc)=str_replace_all(colnames(coloc),"'", "")
#colnames(coloc1)=str_replace_all(colnames(coloc1),"'", "")

#filtering
coloc %<>% filter(TREATMENT=="BUP" & DAYS==7) 
coloc %<>% group_by(Line, `BUP Responder`) %>% filter(DF==max(DF)) 

coloc.identifiers=coloc %>% select(FIELD_NO, DF, Line, `BUP Responder`,  TREATMENT, DAYS, `SN/SR` )
#coloc %<>% filter(`SN/SR`=="SN")
#coloc %<>% select(-FIELD_NO, -DF, -`SN/SR`) %>% group_by(Line, `BUP Responder`) %>% summarise_if(is.numeric, median) 
coloc %<>% ungroup %<>% select(-Line, -FIELD_NO, -DF, -`SN/SR`, -`pix/um`, -TREATMENT, -DAYS)
#coloc %<>% tibble::column_to_rownames("Line")


coloc.bupResponder=coloc$`BUP Responder`
coloc %<>% mutate(`BUP Responder`=ifelse(`BUP Responder`=="NR", 0,1))

coloc %<>% mutate_all(as.numeric)
#coloc = lapply(coloc, is.numeric)

ttests=coloc %>%
  summarise_each(funs(t.test(.[`BUP Responder` == 0], .[`BUP Responder` == 1])$p.value), vars = -'BUP Responder')

ttests=as.data.frame(t(ttests))
#ttests %<>% tibble::rownames_to_column("ensembleID")
#colnames(ttests)[2]="pval"

coloc %<>% ungroup()
coloc.identifiers %<>% ungroup()

# output_filename_predata=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
#                               "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
#                               "04_Differentiation_Experiments/Data Science/Colocolization_modeling/",
#                               "20210315_coloc_Bup_7d_maxDF_premodeling_data.csv")
# 
# write.csv(bind_cols(coloc.identifiers %>% select(FIELD_NO, DF, Line), coloc), 
#           output_filename_predata,
#           quote = TRUE,
#           row.names = FALSE)



set.seed(NULL)
accuracies=data.frame()
predictions=data.frame(response=coloc$`BUP Responder`)
tests=data.frame(response=coloc$`BUP Responder`)
models=list()
#for(i in 1:nrow(coloc))
for(i in 1:20)
{ 
  print(i)
  train=c(sample(which(coloc$`BUP Responder`==0), 74*0.8) ,
          sample(which(coloc$`BUP Responder`==1), 57*0.8) )
  test=c(1:nrow(coloc))[-train]
  
  tt=data.frame(response=coloc$`BUP Responder`)
  tt$trainORtest=""
  tt[train,]="train"
  tt[test,]="test"
  colnames(tt)[2]=str_c("trainORtest_", i)
  tests %<>% bind_cols(tt %>% select(-response))
  # train=c(1:nrow(coloc))[-i]
  # test=c(i)
  
  
  #mylogit <- glm(`BUP Responder` ~ ., data = coloc[train,], family = "binomial")
  
  xgb.model <- xgboost(data = as.matrix(coloc[train,-1]),
                       label = coloc[train,] %>% pull(`BUP Responder`),
                       max.depth = 4, #prms$max.depth[p],
                       eta = 0.1, #prms$eta[p],
                       nthread = 5, #prms$nthreads[p],
                       nrounds = 25, #prms$nrounds[p] , #Why 63? xgb.cv showed several times for it to be best fit
                       eval_metric = "error",
                       objective = "binary:logistic",
                       verbose = 0)
  
  predicted <- predict(xgb.model, as.matrix(coloc[,-1]));
  
  #predicted=predict(mylogit, newdata = coloc[,-1], type = "response")
  
  pp=data.frame(predicted)
  colnames(pp)=str_c("prediction_", ncol(predictions))
  
  models[[str_c("prediction_", ncol(predictions))]]=xgb.model
  
  predictions %<>% bind_cols(pp)
  
  #colnames(predictions)[ncol(predictions)]=str_c("prediction_", ncol(predictions))
  
  #predicted=round(predict(mylogit, newdata = coloc[,-1], type = "response"))
  predicted=round(predict(xgb.model, newdata = as.matrix(coloc[,-1])))
  
  
  cm.train <- confusionMatrix(factor(coloc$`BUP Responder`[train]), 
                              factor(predicted[train]))
  #print("Train")
  #print(cm.train$overall)
  
  #cm.test <- ifelse(predicted[test]==coloc$`BUP Responder`[test], 100, 0);
  
  cm.test <- confusionMatrix(factor(coloc$`BUP Responder`[test]),
                             factor(predicted[test]))
  
  # accuracies %<>% bind_rows(data.frame(nn=i,
  #                                      train=cm.train$overall["Accuracy"], 
  #                                      test=cm.test))
  accuracies %<>% bind_rows(data.frame(nn=i,
                                       train=cm.train$overall["Accuracy"], 
                                       test=cm.test$overall["Accuracy"],
                                       test.sensetivity=cm.test$byClass[["Sensitivity"]],
                                       test.specificity=cm.test$byClass[["Specificity"]]))
}

predictions =bind_cols(coloc.identifiers, predictions)

predictions %<>% group_by_at(vars(!starts_with("prediction_"))) %>% 
  rowwise() %>%
  mutate(mean_prediction=mean(c_across( starts_with("prediction_"))))



# output_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
#                       "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
#                       "04_Differentiation_Experiments/Data Science/Colocolization_modeling/",
#                       "20210314_coloc_Bup_7d_maxDF_modeling_predictions.csv")
# 
# write.csv(predictions, output_filename,
#           quote = TRUE,
#           row.names = FALSE)



print(str_c("mean accuracy:",  mean(accuracies$test)))
print(str_c("mean sensetivity:",  mean(accuracies$test.sensetivity)))
print(str_c("mean specificity:",  mean(accuracies$test.specificity)))


# #output predictions
# 
# pp=predictions %>% select(FIELD_NO, DF, Line, `BUP Responder`, TREATMENT, DAYS, `SN/SR`, response, prediction_6, prediction_16)
# pp %<>% bind_cols(tests %>% select(trainORtest_6, trainORtest_16))
# write.csv(pp, "/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/Exploratory_Outputs/neurobiology_imaging/Coloc/20210426_predictions_coloc_models.csv", row.names = F)
# 
# 
# 
# 
# ## CIT TEMP!!!
# cit=read_excel_allsheets(str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/04_Differentiation_Experiments/Compiled_Results/Colocalization/20210427_by_field_for_Sasha/20210427_CIT_by_Field.xlsx"))
# #coloc1=read_excel_allsheets(str_c(data_folder, "20210311_Bup_by_Field_DF42to74.xlsx"))
# 
# cit=cit$Sheet1
# #coloc1=coloc1$Sheet1
# 
# colnames(cit)=str_replace_all(colnames(cit),"'", "")
# #colnames(coloc1)=str_replace_all(colnames(coloc1),"'", "")
# 
# #filtering
# cit %<>% filter((TREATMENT %in% c("CIT")) & DAYS==7) 
# cit %<>% group_by(Line, `BUP Responder`) %>% filter(DF==max(DF)) 
# 
# cit.identifiers=cit %>% select(FIELD_NO, DF, Line, `BUP Responder`,  TREATMENT, DAYS, `SN/SR` )
# 
# cit %<>% ungroup %<>% select(-Line, -FIELD_NO, -DF, -`SN/SR`, -`pix/um`, -TREATMENT, -DAYS)
# #coloc %<>% tibble::column_to_rownames("Line")
# 
# 
# cit.bupResponder=cit$`BUP Responder`
# cit %<>% rowwise() %<>% mutate(response=ifelse(`BUP Responder`=="NR", 0,1)) %>% relocate(response) %>% select(-`BUP Responder`)
# 
# cit %<>% mutate_all(as.numeric)
# 
# cit %<>% ungroup()
# cit.identifiers %<>% ungroup()
# 
# 
# predicted_6 <- data.frame(response=cit$response, prediction_6=predict(models[["prediction_6"]], as.matrix(cit[,-1])))
# predicted_16 <-data.frame(prediction_16=predict(models[["prediction_16"]], as.matrix(cit[,-1])))
# 
# cit.predictions =bind_cols(cit.identifiers, predicted_6)
# cit.predictions %<>% bind_cols(predicted_16)
# 
# print(confusionMatrix(factor(round(cit.predictions$prediction_6)),
#                       factor(cit.predictions$response), positive = "1"))
# 
# print(confusionMatrix(factor(round(cit.predictions$prediction_16)),
#                       factor(cit.predictions$response), positive = "1"))
# 
# cit.predictions %<>% mutate(pred6_bin=round(prediction_6))
# cit.predictions %<>% mutate(pred16_bin=round(prediction_16))
# 
# borderline=c(10, 32, 27)
# 
# 
# cit.predictions %<>% mutate(response_borderline=ifelse(Line %in% borderline, 1, response))
# 
# print(confusionMatrix(factor(round(cit.predictions$prediction_6)),
#                       factor(cit.predictions$response_borderline), positive = "1"))
# 
# print(confusionMatrix(factor(round(cit.predictions$prediction_16)),
#                       factor(cit.predictions$response_borderline), positive = "1"))
# 
# 
# write.csv(cit.predictions, "/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/Exploratory_Outputs/neurobiology_imaging/Coloc/20210427_predictions_coloc_models_CIT.csv", row.names = F)
# 
# 
# #NTP
# ntp=read_excel_allsheets(str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/04_Differentiation_Experiments/Compiled_Results/Colocalization/20210427_by_field_for_Sasha/20210427_NTP_by_Field.xlsx"))
# #coloc1=read_excel_allsheets(str_c(data_folder, "20210311_Bup_by_Field_DF42to74.xlsx"))
# 
# ntp=ntp$Sheet1
# #coloc1=coloc1$Sheet1
# 
# colnames(ntp)=str_replace_all(colnames(ntp),"'", "")
# #colnames(coloc1)=str_replace_all(colnames(coloc1),"'", "")
# 
# #filtering
# ntp %<>% filter((TREATMENT %in% c("NTP")) & DAYS==7) 
# ntp %<>% group_by(Line, `BUP Responder`) %>% filter(DF==max(DF)) 
# 
# ntp.identifiers=ntp %>% select(FIELD_NO, DF, Line, `BUP Responder`,  TREATMENT, DAYS, `SN/SR` )
# 
# ntp %<>% ungroup %<>% select(-Line, -FIELD_NO, -DF, -`SN/SR`, -`pix/um`, -TREATMENT, -DAYS)
# #coloc %<>% tibble::column_to_rownames("Line")
# 
# 
# ntp.bupResponder=ntp$`BUP Responder`
# ntp %<>% rowwise() %<>% mutate(response=ifelse(`BUP Responder`=="NR", 0,1)) %>% relocate(response) %>% select(-`BUP Responder`)
# 
# ntp %<>% mutate_all(as.numeric)
# 
# ntp %<>% ungroup()
# ntp.identifiers %<>% ungroup()
# 
# 
# ntp_predicted_6 <- data.frame(response=ntp$response, prediction_6=predict(models[["prediction_6"]], as.matrix(ntp[,-1])))
# ntp_predicted_16 <-data.frame(prediction_16=predict(models[["prediction_16"]], as.matrix(ntp[,-1])))
# 
# ntp.predictions =bind_cols(ntp.identifiers, ntp_predicted_6)
# ntp.predictions %<>% bind_cols(ntp_predicted_16)
# 
# print(confusionMatrix(factor(round(ntp.predictions$prediction_6)),
#                       factor(ntp.predictions$response), positive = "1"))
# 
# print(confusionMatrix(factor(round(ntp.predictions$prediction_16)),
#                       factor(ntp.predictions$response), positive = "1"))
# 
# ntp.predictions %<>% mutate(pred6_bin=round(prediction_6))
# ntp.predictions %<>% mutate(pred16_bin=round(prediction_16))
# 
# 
# write.csv(ntp.predictions, "/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/Exploratory_Outputs/neurobiology_imaging/Coloc/20210427_predictions_coloc_models_NTP.csv", row.names = F)
# 
