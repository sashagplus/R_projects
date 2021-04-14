library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library("rpart")
library(rpart.plot)

rm(list=ls())

#read data
rnaseq_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                      "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/rnaseq/primary_models/",
                      "20210315_rnaseq_Bup_7d_premodeling_data.csv")

rnaseq=read.csv(rnaseq_filename)

spines_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                      "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/",
                      "04_Differentiation_Experiments/Data Science/Spines_modeling/",
                      "20210315_spines_Bup_7d_maxDF_premodeling_data.csv")

spines=read.csv(spines_filename)
spines %<>% separate(name, into=c("DF", "Line", "Well", "Days", "fn_Well", "fn_Field"))

spines %<>% select(-DF,-Well,-Days, -fn_Well, -fn_Field)

spines %<>% group_by(Line, BUP.Responder) %>% summarise_all(median)


coloc_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                     "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/",
                     "04_Differentiation_Experiments/Data Science/Colocolization_modeling/",
                     "20210315_coloc_Bup_7d_maxDF_premodeling_data.csv")

coloc=read.csv(coloc_filename)

coloc %<>% select(-FIELD_NO,-DF)

coloc %<>% group_by(Line, BUP.Responder) %>% summarise_all(median)


# coloc %<>% select(FIELD_NO, Line, BUP.Responder, response, mean_prediction)
# rnaseq %<>% select(Line, response, mean_prediction)
# spines %<>% select(Line, BUP.Responder, response, mean_prediction)
#
 rnaseq %<>% mutate(Line=as.numeric(gsub('L', '', Line)))

#look into data availability
 alllines_names=full_join(full_join(spines %>% select(Line, BUP.Responder) %>% mutate(Spines="V")  %>% mutate(Line=as.numeric(Line)),
                              coloc %>% select(Line, BUP.Responder) %>% mutate(Coloc="V") ,
                              by=c("Line", "BUP.Responder")),
                    rnaseq %>% select(Line) %>% mutate(RNASeq="V") ,
                    by=c("Line"))
 
 alllines_names %<>% arrange(Line) %>% 
                mutate(all_present=ifelse(Spines=="V" & Coloc=="V" & RNASeq=="V", "V", "")) %>%
              arrange(all_present)
write.csv(alllines_names, 
          stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                          "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/",
                          "Integration/Data/data_availability_gplus_7d_bupropion.csv"),
          quote = TRUE,
          row.names = FALSE)
break()
 

alllines=full_join(rnaseq  ,
                    full_join(spines %>% mutate(Line=as.numeric(Line)),
                              coloc  ,
                              by=c("Line", "BUP.Responder")),
                    by=c("Line", "group"="BUP.Responder"))


alllines %<>% drop_na() 

xgb.model <- xgboost(data = as.matrix(alllines[,-c(1:2)]),
                     label = alllines[,2],
                     max.depth = 2, #prms$max.depth[p],
                     eta = 0.2, #prms$eta[p],
                     nthread = 5, #prms$nthreads[p],
                     nrounds = 30, #prms$nrounds[p] , #Why 63? xgb.cv showed several times for it to be best fit
                     eval_metric = "error",
                     objective = "binary:logistic",
                     verbose = 0)
mylogit <- glm(`group` ~ ., data = alllines[,-1], family = "binomial")

predicted <- predict(mylogit, (alllines[,-1]));


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
library(ggplot2)
alllines.pca <- prcomp(alllines[,-c(1:2)])
summary(alllines.pca)

autoplot(prcomp(alllines[,-c(1:2)]), data=alllines[,-1], col="group")



