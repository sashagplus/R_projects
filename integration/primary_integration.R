library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library("rpart")
library(rpart.plot)

rm(list=ls())


rnaseq_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                      "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/rnaseq/primary_models/",
                      "20210315_rnaseq_Bup_7d_premodeling_data.csv")

rnaseq=read.csv(rnaseq_filename)

spines_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                      "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
                      "04_Differentiation_Experiments/Data Science/Spines_modeling/",
                      "20210315_spines_Bup_7d_maxDF_premodeling_data.csv")

spines=read.csv(spines_filename)

coloc_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                     "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
                     "04_Differentiation_Experiments/Data Science/Colocolization_modeling/",
                     "20210315_coloc_Bup_7d_maxDF_premodeling_data.csv")

coloc=read.csv(coloc_filename)


# coloc %<>% select(FIELD_NO, Line, BUP.Responder, response, mean_prediction)
# rnaseq %<>% select(Line, response, mean_prediction)
# spines %<>% select(Line, BUP.Responder, response, mean_prediction)
# 
# rnaseq %<>% mutate(Line=as.numeric(gsub('L', '', Line)))

# rname=data.frame(rnum=1:nrow(rnaseq), BUP=rnaseq$group)
# cname=data.frame(cnum=1:nrow(coloc), BUP=coloc$BUP.Responder)
# sname=data.frame(snum=1:nrow(spines), BUP=spines$BUP.Responder)

responsive=merge(merge(rnaseq %>% filter(group==1) %>% select(-group),
                 coloc %>% filter(BUP.Responder==1) %>% select(-BUP.Responder)),
                 spines %>% filter(BUP.Responder==1) %>% select(-BUP.Responder))

responsive %<>% mutate(BUP=1) %>% select(BUP, everything())


nonresponsive=merge(merge(rnaseq %>% filter(group==0) %>% select(-group),
                       coloc %>% filter(BUP.Responder==0) %>% select(-BUP.Responder)),
                 spines %>% filter(BUP.Responder==0) %>% select(-BUP.Responder))

nonresponsive %<>% mutate(BUP=0) %>% select(BUP, everything())



simulated=bind_rows(responsive, nonresponsive)


#simulated$BUP = sample(simulated$BUP)
accuracies.xgb=data.frame()
#predictions=data.frame(response=dndrt_data$BUP.Responder)
for(i in 1:10)
{  
  train=c(sample(which(simulated$BUP==0), length(which(simulated$BUP==0))*0.5) ,
          sample(which(simulated$BUP==1), length(which(simulated$BUP==1))*0.5) )
  test=c(1:nrow(simulated))[-train]
  
  # train=c(1:nrow(dndrt_data))[-1]
  # test=c(i)
  
  xgb.model <- xgboost(data = as.matrix(simulated[train,-1]),
                       label = simulated[train,1],
                       max.depth = 2, #prms$max.depth[p],
                       eta = 0.2, #prms$eta[p],
                       nthread = 5, #prms$nthreads[p],
                       nrounds = 30, #prms$nrounds[p] , #Why 63? xgb.cv showed several times for it to be best fit
                       eval_metric = "error",
                       objective = "binary:logistic",
                       verbose = 0)
  
  #predicted <- predict(xgb.model, as.matrix(simulated[,-1]));
  
  #pp=data.frame(predicted)
  #colnames(pp)=str_c("prediction_", ncol(predictions))
  
  #predictions %<>% bind_cols(pp)
  
  predicted <- round(predict(xgb.model, as.matrix(simulated[,-1])));
  
  cm.train <- confusionMatrix(factor(simulated$BUP[train], levels = unique(simulated$BUP)), 
                              factor(predicted[train],levels = unique(simulated$BUP)))
  #print("Train")
  #print(cm.train$overall)
  
  cm.test <- confusionMatrix(factor(simulated$BUP[test], levels = unique(simulated$BUP)),
                                        factor(predicted[test], levels = unique(simulated$BUP)))
  
  accuracies.xgb %<>% bind_rows(data.frame(nn=i,
                                       train=cm.train$overall["Accuracy"], 
                                       test=cm.test$overall["Accuracy"]))
}
  
print(mean(accuracies.xgb$test))
