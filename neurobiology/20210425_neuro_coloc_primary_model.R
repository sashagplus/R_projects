# Analyze and model colocolization data
library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library("rpart")
library(rpart.plot)
library(readxl)    
library(caret)
rm(list=ls())


#reading input
files_treatToVehicle=read.csv(
          str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science",
                "/Exploratory_Outputs/neurobiology_imaging/Coloc/",
                "20210425_coloc_treatment_to_vehicle_files_mapping.csv"))

rois=read.csv(str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science",
                "/Exploratory_Outputs/neurobiology_imaging/Coloc/",
                "20210425_coloc_rois_collected_data.csv"))

feature_inspect=read.csv(str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                               "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/",
                               "Exploratory_Outputs/neurobiology_imaging/Coloc/",
                               "20210421_coloc_feature_inspection_for_modeling_SN.csv"))

samples_response=read.csv(str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                                "Genetika+SharedDrive/01_Protocol_Development/",
                                "01_05_Imaging/01_Experiments_and_Results/04_Differentiation_Experiments/",
                                "Line_summary_treatments.csv"))

samples_response %<>% rename(Line=L) %>% 
  gather(key="Responder_class", value="Responder_value", -DF, -Line) %>%
  rowwise() %>% 
  mutate(Responder_class=str_split(Responder_class, "_")[[1]][1]) 

feature_inspect %<>% filter(!feature_name %in% c("folder.coloc", "image_nROI_Length", "PSD", "SYN")) 
feature_inspect %<>% filter(to_model=="Yes" | is_identifier==TRUE) 
feature_inspect %<>% mutate(feature_name=str_replace_all(feature_name, " ", ".")) 

feature_inspect %<>% filter(file_type=="ROI")


#1. create for all features "to_model" a normalized field to vehicle
# if interesting - check correlation between all features
#2. Model with evertyhing


#1. create for all features "to_model" a normalized field to vehicle
#rois %>% group_by(DF, Line, Treat, Field) %>%  tally() %>% pull(Treat) %>% table()
#rois %>% group_by(DF, Line, Treat) %>%  tally() %>% pull(Treat) %>% table()
#Not enough data to work by lines, thu will be working by field.

#Each treatment field will be normalized towards the median of the VEH
features_to_model = feature_inspect %>% filter(is_identifier==FALSE) %>% pull(feature_name)
features_identifiers = feature_inspect %>% filter(is_identifier==TRUE) %>% pull(feature_name)

#Remove dendrites having 0 or NA length
rois %<>% filter(!(is.na(Dendrite.length.micron) & is.na(Total.dendrite.length.uM.DANA))) 

#rois %>% filter(!is.na(Dendrite.length.micron) & !is.na(Total.dendrite.length.uM.DANA)) %>% View
rois %<>% mutate(Dendrite.length.micron=replace_na(Dendrite.length.micron, 0), 
                Total.dendrite.length.uM.DANA=replace_na(Total.dendrite.length.uM.DANA, 0)) %>%
        mutate(Dendrite.length.micron=Dendrite.length.micron + Total.dendrite.length.uM.DANA) %>%
  filter(Dendrite.length.micron!=0)

print(rois %>% nrow)
print(table(rois$Treat))


rois.vehicle = rois %>% filter(Treat=="VEH")
#fix Days 73 into 7 and 3
rois.vehicle %<>% tibble::rowid_to_column("rowid") 
rois.vehicle= rois.vehicle %>%
  select(rowid, Days) %>%
  mutate(Days=as.character(Days)) %>% 
  mutate(Days=ifelse(Days=="73", "7 3", Days)) %>% 
  separate(Days, into=c("seven", "three"), sep=" ", remove = F, fill="left") %>% 
  gather(key, value, -rowid, -Days, na.rm = T) %>%
  mutate(Days=as.integer(value)) %>% 
  select(-key, -value) %>%
  right_join(rois.vehicle %>% select(-Days), by=c("rowid")) %>%
  select(-rowid) %>%
  as.data.frame()


rois %<>% gather(key="feature_to_model", value="value", all_of(features_to_model)) 
rois.vehicle %<>% gather(key="feature_to_model", value="value", all_of(features_to_model)) 

rois.vehicle %<>% group_by(DF, Line, Days, BUP.Responder, roi_filename, feature_to_model) %<>%
                  summarise(median.value=median(value, na.rm=T), .groups = "keep")


rois.treatment = rois %>% filter(Treat!="VEH")
rois.treatment %<>% left_join(rois.vehicle, 
                              by=c("DF",  "Line", "Days", "BUP.Responder", "roi_filename", "feature_to_model")) 


rois.treatment %<>% bind_rows(rois.treatment %>% mutate(feature_to_model=str_c(feature_to_model, "_subtract"), 
                                                       value=value-median.value)) %<>%
                  bind_rows(rois.treatment %>% mutate(feature_to_model=str_c(feature_to_model, "_fold"), 
                                           value=ifelse(median.value==0, value, value/median.value)))


rois.treatment %<>% select(-median.value) %<>%  spread(key=feature_to_model, value=value)

rois.treatment %<>% left_join(samples_response, by=c("DF", "Line", "Treat"="Responder_class")) %>% relocate(Responder_value) %>% select(-BUP.Responder) 


rois.treatment %<>% rowwise() %>% mutate(Responder_value_rnr= str_trim(str_split(Responder_value, "(_)")[[1]][1]),
                                        Responder_value_bl=str_remove(str_split(Responder_value, "\\(")[[1]][2], "\\)")) %>% 
                                        relocate(Responder_value_rnr, Responder_value_bl) 
rois.treatment %<>% mutate(response=ifelse(Responder_value=="R", 1, 0)) %>% relocate(response)
rois.treatment %<>% group_by(Treat) %>% mutate(treat_id=cur_group_id()-1) %<>% ungroup() %>% relocate(treat_id)



rois.treatment.identifiers=rois.treatment %>% select(treat_id:ImageNumber)

features_to_model=colnames(rois.treatment)[which(!colnames(rois.treatment) %in% colnames(rois.treatment.identifiers))]






modelmat=rois.treatment %>% select(treat_id, !!features_to_model)

# #all 4 classes
# trainIndex <- caret::createDataPartition(as.matrix(modelmat[,1]), 
#                                          p = .75, 
#                                          list = FALSE, 
#                                          times = 1)
# testIndex=c(1:nrow(modelmat))[-trainIndex]
# 
# #Selecting relevant features
# gdata.modeling=model.matrix(treat_id~. , modelmat)[,-1]
# 
# cv_output <- cv.glmnet(gdata.modeling[trainIndex,], 
#                        modelmat$treat_id[trainIndex],
#                        alpha = 1, 
#                        family="multinomial",
#                        nfolds = 10)
# 
# best_lam=cv_output$lambda.min
# 
# lasso_best <- glmnet(gdata.modeling[trainIndex,], 
#                      modelmat$treat_id[trainIndex], 
#                      alpha = 1, 
#                      lambda = best_lam,
#                      family="multinomial"
# )
# 
# pred <- predict(lasso_best, 
#                 s = best_lam, 
#                 newx = gdata.modeling,
#                 type="class")
# 
# pred.df=data.frame(res.original=modelmat$treat_id,
#                    res.predicted=as.numeric(pred),
#                    stringsAsFactors = F)
# 
# 
# #C. Identify informative features using LASSO coefficients 
# cfcnt=coef(lasso_best,s=best_lam ,exact=TRUE) 
# features_selected=c()
# for(i in 1:4)
# {
#   fs_i=cfcnt[[i]] %>% as.matrix() %>% as.data.frame() %>% tibble::rownames_to_column("feature") %>% rename(prob= 2) %>% filter(prob==0) %>% pull(feature)
#   features_selected=c(features_selected, fs_i)
# }
# 
# features_selected=unique(features_selected)
# 
# xgb.model <- xgboost(data = as.matrix(gdata.modeling[trainIndex,features_selected]),
#                      label = modelmat$treat_id[trainIndex],
#                      max.depth = 4, #prms$max.depth[p],
#                      eta = 0.1, #prms$eta[p],
#                      nthread = 5, #prms$nthreads[p],
#                      nrounds = 25, #prms$nrounds[p] , #Why 63? xgb.cv showed several times for it to be best fit
#                      eval_metric = "merror",
#                      objective = "multi:softmax",
#                      verbose = 0,
#                      num_class=4)
# 
# predicted <- predict(xgb.model, as.matrix(as.matrix(gdata.modeling[trainIndex,features_selected])))
# 
# cm.train <- confusionMatrix(factor(modelmat$treat_id[trainIndex]), 
#                             factor(predicted[trainIndex]))
# 
# 
# cm.test <- confusionMatrix(factor(modelmat$treat_id[testIndex]), 
#                            factor(predicted[testIndex]))
# 



set.seed(1)
#BUP vs CIT
#Selecting relevant features
modelmat.backup=modelmat
modelmat %<>% filter(treat_id %in% c(0, 1))

ttests=modelmat %>%
  summarise_each(funs(t.test(.[treat_id == 0], .[treat_id == 1])$p.value), vars = -treat_id)

ttests=as.data.frame(t(ttests))
ttests %<>% tibble::rownames_to_column("feature") %>% 
  rename(pval=2) %>%
  mutate(pval=round(pval, 3)) 

features_selected=ttests %>% filter(pval<=0.005) %>% pull(feature)

trainIndex <- caret::createDataPartition(as.matrix(modelmat[,1]), 
                                         p = .75, 
                                         list = FALSE, 
                                         times = 1)
testIndex=c(1:nrow(modelmat))[-trainIndex]

gdata.modeling=model.matrix(treat_id~. , modelmat)[,-1]

# cv_output <- cv.glmnet(gdata.modeling[trainIndex,], 
#                        modelmat$treat_id[trainIndex],
#                        alpha = 1, 
#                        family="binomial",
#                        nfolds = 10)
# 
# best_lam=cv_output$lambda.min
# 
# lasso_best <- glmnet(gdata.modeling[trainIndex,], 
#                      modelmat$treat_id[trainIndex], 
#                      alpha = 1, 
#                      lambda = best_lam,
#                      family="multinomial"
# )
# 
# pred <- predict(lasso_best, 
#                 s = best_lam, 
#                 newx = gdata.modeling,
#                 type="class")
# 
# pred.df=data.frame(res.original=modelmat$treat_id,
#                    res.predicted=as.numeric(pred),
#                    stringsAsFactors = F)


# #C. Identify informative features using LASSO coefficients 
# cfcnt=coef(lasso_best,s=best_lam ,exact=TRUE) 
# features_selected=c()
# for(i in 1:2)
# {
#   fs_i=cfcnt[[i]] %>% as.matrix() %>% as.data.frame() %>% tibble::rownames_to_column("feature") %>% rename(prob= 2) %>% filter(prob==0) %>% pull(feature)
#   features_selected=c(features_selected, fs_i)
# }
# 
# features_selected=unique(features_selected)
# accuracies=data.frame()
# predictions=data.frame(response=coloc$`BUP Responder`)
# tests=data.frame(response=coloc$`BUP Responder`)
# models=list()
#for(i in 1:nrow(coloc))
#for(i in 1:20)
{ 
  
  xgb.model <- xgboost(data = as.matrix(gdata.modeling[trainIndex,features_selected]),
                       label = modelmat$treat_id[trainIndex],
                       max.depth = 4, #prms$max.depth[p],
                       eta = 0.1, #prms$eta[p],
                       nthread = 5, #prms$nthreads[p],
                       nrounds = 25, #prms$nrounds[p] , #Why 63? xgb.cv showed several times for it to be best fit
                       eval_metric = "error",
                       objective = "binary:logistic",
                       verbose = 0)
  
  predicted <- predict(xgb.model, as.matrix(as.matrix(gdata.modeling[,features_selected])))
  
  cm.train <- confusionMatrix(factor(modelmat$treat_id[trainIndex]), 
                              factor(round(predicted[trainIndex])))
  
  
  cm.test <- confusionMatrix(factor(modelmat$treat_id[testIndex]), 
                             factor(round(predicted[testIndex])))
  
  
  
  print(cm.train)
  print(cm.test)

}

