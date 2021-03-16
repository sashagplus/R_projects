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

output_filename_predata=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                              "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
                              "04_Differentiation_Experiments/Data Science/Colocolization_modeling/",
                              "20210315_coloc_Bup_7d_maxDF_premodeling_data.csv")

write.csv(bind_cols(coloc.identifiers %>% select(FIELD_NO, DF, Line), coloc), 
          output_filename_predata,
          quote = TRUE,
          row.names = FALSE)




accuracies=data.frame()
predictions=data.frame(response=coloc$`BUP Responder`)
#for(i in 1:nrow(coloc))
for(i in 1:20)
{ 
  print(i)
  train=c(sample(which(coloc$`BUP Responder`==0), 74*0.8) ,
          sample(which(coloc$`BUP Responder`==1), 57*0.8) )
  test=c(1:nrow(coloc))[-train]
  
  # train=c(1:nrow(coloc))[-i]
  # test=c(i)
  
  
  mylogit <- glm(`BUP Responder` ~ ., data = coloc[train,], family = "binomial")
  
  predicted=predict(mylogit, newdata = coloc[,-1], type = "response")
  
  pp=data.frame(predicted)
  colnames(pp)=str_c("prediction_", ncol(predictions))
  
  predictions %<>% bind_cols(pp)
  
  #colnames(predictions)[ncol(predictions)]=str_c("prediction_", ncol(predictions))
  
  predicted=round(predict(mylogit, newdata = coloc[,-1], type = "response"))
  
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
                                       test=cm.test$overall["Accuracy"]))
}

predictions =bind_cols(coloc.identifiers, predictions)

predictions %<>% group_by_at(vars(!starts_with("prediction_"))) %>% 
                rowwise() %>%
                 mutate(mean_prediction=mean(c_across( starts_with("prediction_"))))



output_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                      "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
                      "04_Differentiation_Experiments/Data Science/Colocolization_modeling/",
                      "20210314_coloc_Bup_7d_maxDF_modeling_predictions.csv")

write.csv(predictions, output_filename,
          quote = TRUE,
          row.names = FALSE)



print(mean(accuracies$test))








