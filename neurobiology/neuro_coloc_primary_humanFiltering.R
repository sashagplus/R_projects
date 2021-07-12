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
                  "04_Differentiation_Experiments/Compiled_Results/Colocalization/20210427_by_field_for_Sasha/")


coloc=read_excel_allsheets(str_c(data_folder, "20210428_BUP_CIT_NTP_by_Field.xlsx"))
#coloc1=read_excel_allsheets(str_c(data_folder, "20210311_Bup_by_Field_DF42to74.xlsx"))

coloc=coloc$Sheet1
#coloc1=coloc1$Sheet1

colnames(coloc)=str_replace_all(colnames(coloc),"'", "")
#colnames(coloc1)=str_replace_all(colnames(coloc1),"'", "")

coloc %<>% mutate(TREATMENT=str_replace_all(TREATMENT, "'", "")) %<>%
          mutate(TREATMENT=ifelse(TREATMENT=="CIT_7D", "CIT", TREATMENT)) %<>%
          mutate(TREATMENT=ifelse(TREATMENT=="MIRT_7D", "MIRT", TREATMENT))


#filtering
#coloc %<>% filter(! TREATMENT %in% c("MIRT", "VEH", "NTP"))
coloc %<>% filter(! TREATMENT %in% c( "VEH", "MIRT"))
#coloc %<>% filter(TREATMENT=="BUP" & DAYS==7) 
coloc %<>% group_by(Line, `BUP Responder`) %>% filter(DF==max(DF)) %<>% ungroup()

set.seed(1)


# train <- caret::createDataPartition(as.matrix(coloc %>%  group_by(TREATMENT, `BUP Responder`) %>% mutate(group_ind=group_indices()) %>% pull(group_ind)), 
#                                          p = .80, 
#                                          list = FALSE, 
#                                          times = 1)
# test=c(1:nrow(coloc))[-train]


coloc.identifiers=coloc %>% select(FIELD_NO, DF, Line, `BUP Responder`,  TREATMENT, DAYS, `SN/SR` )
#coloc %<>% filter(`SN/SR`=="SN")
#coloc %<>% select(-FIELD_NO, -DF, -`SN/SR`) %>% group_by(Line, `BUP Responder`) %>% summarise_if(is.numeric, median) 
coloc %<>% ungroup %<>% select(-Line, -FIELD_NO, -DF, -`SN/SR`, -`pix/um`, -TREATMENT, -DAYS)
#coloc %<>% tibble::column_to_rownames("Line")


coloc.bupResponder=coloc$`BUP Responder`
coloc %<>% mutate(`BUP Responder`=ifelse(`BUP Responder`=="NR", 0,1))

coloc %<>% mutate_all(as.numeric)
#coloc = lapply(coloc, is.numeric)

# ttests=coloc %>%
#   summarise_each(funs(t.test(.[`BUP Responder` == 0], .[`BUP Responder` == 1])$p.value), vars = -'BUP Responder')
# 
# ttests=as.data.frame(t(ttests))
#ttests %<>% tibble::rownames_to_column("ensembleID")
#colnames(ttests)[2]="pval"



coloc %<>% ungroup()
coloc.identifiers %<>% ungroup()

# output_filename_predata=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
#                               "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
#                               "04_Differentiation_Experiments/Data Science/Colocolization_modeling/",
#                               "20210427_coloc_Bup_CIT_NTP_maxDF_premodeling_data.csv")
# 
# write.csv(bind_cols(coloc.identifiers %>% select(FIELD_NO, DF, Line), coloc), 
#           output_filename_predata,
#           quote = TRUE,
#           row.names = FALSE)




accuracies=data.frame()
predictions=data.frame(response=coloc$`BUP Responder`)
tests=data.frame(response=coloc$`BUP Responder`)
models=list()
predictions_perTreatment=data.frame()
#for(i in 1:nrow(coloc))
for(i in 1:20)
{ 
  print(i)
  train <- caret::createDataPartition(as.matrix(coloc.identifiers %>%  group_by(TREATMENT, `BUP Responder`) %>% mutate(group_ind=group_indices()) %>% pull(group_ind)), 
                                      p = .80, 
                                      list = FALSE, 
                                      times = 1)
  test=c(1:nrow(coloc))[-train]
  
  #keep the train/test partitioning
  tt=data.frame(response=coloc$`BUP Responder`)
  tt$trainORtest=""
  tt[train,]="train"
  tt[test,]="test"
  colnames(tt)[2]="trainORtest"
  
  
  #Model
  xgb.model <- xgboost(data = as.matrix(coloc[train,-1]),
                       label = coloc[train,] %>% pull(`BUP Responder`),
                       max.depth = 4, #prms$max.depth[p],
                       eta = 0.1, #prms$eta[p],
                       nthread= 5, #prms$nthreads[p],
                       nrounds = 25, #prms$nrounds[p] , #Why 63? xgb.cv showed several times for it to be best fit
                       eval_metric = "error",
                       objective = "binary:logistic",
                       verbose = 0)

  predicted <- predict(xgb.model, as.matrix(coloc[,-1]));
  
  
  #Predict
  pp=data.frame(predicted_probability=predicted)
  pp %<>% mutate(predicted_class=ifelse(round(predicted_probability)==1, "R", "NR"))
  
  #recording predictions per treatment
  predictions_perTreatment %<>% bind_rows(bind_cols(coloc.identifiers %>% select(TREATMENT, `BUP Responder`), 
                                                    pp %>% select(predicted_class),
                                                    tt %>% select(-1)) %>%
                                            ftable() %>% 
                                            as.data.frame() %>% ungroup() %>%
                                            mutate(temp=str_c(trainORtest, "_", predicted_class), .keep = c( "unused")) %>% 
                                            spread(key=temp, value=Freq) %>% 
                                            mutate(model=i, .before=c("TREATMENT"))
  )
  
  
  
  pp %<>% rename_all(~ str_c(.x, "_", i))
  tt %<>% rename_all(~ str_c(.x, "_", i))
  
  models[[str_c("prediction_", ncol(predictions))]]=xgb.model
  
  predictions %<>% bind_cols(pp, tt %>% select(-1))
 
  predicted=round(predict(xgb.model, newdata = as.matrix(coloc[,-1])))
  
  
  #recording model performance
  cm.train <- confusionMatrix(factor(coloc$`BUP Responder`[train]), 
                              factor(predicted[train]))
  
  cm.test <- confusionMatrix(factor(coloc$`BUP Responder`[test]),
                             factor(predicted[test]))
  
  accuracies %<>% bind_rows(bind_cols( data.frame(nn=i,
                                       train.accuracy=cm.train$overall["Accuracy"], 
                                       test.accuracy=cm.test$overall["Accuracy"],
                                        test.sensetivity=cm.test$byClass[["Sensitivity"]],
                                        test.specificity=cm.test$byClass[["Specificity"]]),
                                      as.data.frame(cm.train$table) %>% 
                                        mutate_at(.vars=vars(c("Prediction", "Reference")), ~as.character(ifelse(.x==1, "R", "NR"))) %>% 
                                        mutate(temp=str_c("prediction", Prediction, "_", "actual",Reference), .keep=c("unused")) %>% 
                                        spread(key=temp, value=Freq) %>% 
                                        rename_all(~ str_c("train_", .x)),
                                      as.data.frame(cm.test$table) %>% 
                                        mutate_at(.vars=vars(c("Prediction", "Reference")), ~as.character(ifelse(.x==1, "R", "NR"))) %>% 
                                        mutate(temp=str_c("prediction", Prediction, "_", "actual",Reference), .keep=c("unused")) %>% 
                                        spread(key=temp, value=Freq) %>% 
                                        rename_all(~ str_c("test_", .x))
                                      )
  )
}

predictions =bind_cols(coloc.identifiers, predictions)

# predictions %<>% group_by_at(vars(!starts_with("prediction_"))) %>% 
#                 rowwise() %>%
#                  mutate(mean_prediction=mean(c_across( starts_with("prediction_"))))



predictions_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                      "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
                      "04_Differentiation_Experiments/Data Science/Colocolization_modeling/",
                      "20210428_coloc_Bup_CIT_NTP_7d_maxDF_modeling_predictions_RvsNR.csv")

write.csv(predictions, predictions_filename,
          quote = TRUE,
          row.names = FALSE)

accuracies_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                      "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
                      "04_Differentiation_Experiments/Data Science/Colocolization_modeling/",
                      "20210428_coloc_Bup_CIT_NTP_7d_maxDF_modelSummary.csv")

write.csv(accuracies, accuracies_filename,
          quote = TRUE,
          row.names = FALSE)

predictions_perTreatment_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                          "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
                          "04_Differentiation_Experiments/Data Science/Colocolization_modeling/",
                          "20210428_coloc_Bup_CIT_NTP_7d_maxDF_predictions_perTreatment.csv")

write.csv(predictions_perTreatment, predictions_perTreatment_filename,
          quote = TRUE,
          row.names = FALSE)


print(str_c("mean accuracy:",  mean(accuracies$test.accuracy)))
print(str_c("mean sensetivity:",  mean(accuracies$test.sensetivity)))
print(str_c("mean specificity:",  mean(accuracies$test.specificity)))


#model 8 and model 18

print(confusionMatrix(factor(round(predictions %>% filter(trainORtest_8=="test") %>% pull(predicted_probability_8))),
                      factor(predictions %>% filter(trainORtest_8=="test") %>% pull(response)),
                             positive = "1"))








break()

## CIT TEMP!!!
cit=coloc[which(coloc.identifiers$TREATMENT=="VEH", )]
cit.identifiers=cit %>% select(FIELD_NO, DF, Line, `BUP Responder`,  TREATMENT, DAYS, `SN/SR` )

cit %<>% ungroup %<>% select(-Line, -FIELD_NO, -DF, -`SN/SR`, -`pix/um`, -TREATMENT, -DAYS)
#coloc %<>% tibble::column_to_rownames("Line")


cit.bupResponder=cit$`BUP Responder`
cit %<>% rowwise() %<>% mutate(response=ifelse(`BUP Responder`=="NR", 0,1)) %>% relocate(response) %>% select(-`BUP Responder`)

cit %<>% mutate_all(as.numeric)

cit %<>% ungroup()
cit.identifiers %<>% ungroup()


predicted_6 <- data.frame(response=cit$response, prediction_6=predict(models[["prediction_6"]], as.matrix(cit[,-1])))
predicted_16 <-data.frame(prediction_16=predict(models[["prediction_16"]], as.matrix(cit[,-1])))

cit.predictions =bind_cols(cit.identifiers, predicted_6)
cit.predictions %<>% bind_cols(predicted_16)

print(confusionMatrix(factor(round(cit.predictions$prediction_6)),
                      factor(cit.predictions$response), positive = "1"))

print(confusionMatrix(factor(round(cit.predictions$prediction_16)),
                      factor(cit.predictions$response), positive = "1"))

cit.predictions %<>% mutate(pred6_bin=round(prediction_6))
cit.predictions %<>% mutate(pred16_bin=round(prediction_16))

borderline=c(10, 32, 27)


cit.predictions %<>% mutate(response_borderline=ifelse(Line %in% borderline, 1, response))

print(confusionMatrix(factor(round(cit.predictions$prediction_6)),
                      factor(cit.predictions$response_borderline), positive = "1"))

print(confusionMatrix(factor(round(cit.predictions$prediction_16)),
                      factor(cit.predictions$response_borderline), positive = "1"))


write.csv(cit.predictions, "/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/Exploratory_Outputs/neurobiology_imaging/Coloc/20210427_predictions_coloc_models_CIT.csv", row.names = F)


#NTP
ntp=read_excel_allsheets(str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/04_Differentiation_Experiments/Compiled_Results/Colocalization/20210427_by_field_for_Sasha/20210427_NTP_by_Field.xlsx"))
#coloc1=read_excel_allsheets(str_c(data_folder, "20210311_Bup_by_Field_DF42to74.xlsx"))

ntp=ntp$Sheet1
#coloc1=coloc1$Sheet1

colnames(ntp)=str_replace_all(colnames(ntp),"'", "")
#colnames(coloc1)=str_replace_all(colnames(coloc1),"'", "")

#filtering
ntp %<>% filter((TREATMENT %in% c("NTP")) & DAYS==7) 
ntp %<>% group_by(Line, `BUP Responder`) %>% filter(DF==max(DF)) 

ntp.identifiers=ntp %>% select(FIELD_NO, DF, Line, `BUP Responder`,  TREATMENT, DAYS, `SN/SR` )

ntp %<>% ungroup %<>% select(-Line, -FIELD_NO, -DF, -`SN/SR`, -`pix/um`, -TREATMENT, -DAYS)
#coloc %<>% tibble::column_to_rownames("Line")


ntp.bupResponder=ntp$`BUP Responder`
ntp %<>% rowwise() %<>% mutate(response=ifelse(`BUP Responder`=="NR", 0,1)) %>% relocate(response) %>% select(-`BUP Responder`)

ntp %<>% mutate_all(as.numeric)

ntp %<>% ungroup()
ntp.identifiers %<>% ungroup()


ntp_predicted_6 <- data.frame(response=ntp$response, prediction_6=predict(models[["prediction_6"]], as.matrix(ntp[,-1])))
ntp_predicted_16 <-data.frame(prediction_16=predict(models[["prediction_16"]], as.matrix(ntp[,-1])))

ntp.predictions =bind_cols(ntp.identifiers, ntp_predicted_6)
ntp.predictions %<>% bind_cols(ntp_predicted_16)

print(confusionMatrix(factor(round(ntp.predictions$prediction_6)),
                      factor(ntp.predictions$response), positive = "1"))

print(confusionMatrix(factor(round(ntp.predictions$prediction_16)),
                      factor(ntp.predictions$response), positive = "1"))

ntp.predictions %<>% mutate(pred6_bin=round(prediction_6))
ntp.predictions %<>% mutate(pred16_bin=round(prediction_16))


write.csv(ntp.predictions, "/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/Exploratory_Outputs/neurobiology_imaging/Coloc/20210427_predictions_coloc_models_NTP.csv", row.names = F)



