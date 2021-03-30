#Genetika+ - exercise

library(dplyr)
library(magrittr)
library(caret)
library(glmnet)
library(xgboost)

rm(list=ls())
set.seed(1)

#read input files function
#
#1. read filefirst 
#2. common minimal processing: 
#    a. unify features format tolower
#    b. mark feature names with origin
#    c. remove date and crcid, since non informative
read_csv=function(filename, 
                  filepath=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                                "01_Protocol_Development/01_10_Data_Science/DataforPatientHistory/DepressionDataNIH/",
                                "StarD/STAR_DClinicalData/ASCII Files raw data/"),
                  withDate=FALSE)
{
     
     readdata = read.csv(file = stringr::str_c(filepath, filename),
                          header = T,
                          stringsAsFactors = F);
     
     names(readdata) %<>% tolower;
     
     readdata %<>% dplyr::select(-dplyr::matches("crcid"))
     if(!withDate) 
     {
         readdata %<>% dplyr::select(-dplyr::matches("date"))
     }
     
     
     sffx=stringr::str_c(".",
                         stringr::str_split(stringr::str_split(basename(filename),
                                            pattern = "\\.")[[1]][1],
                         pattern="_")[[1]][1])
     
     readdata %<>% dplyr::rename_at(dplyr::vars(-c(.data$id)),
                                    ~ stringr::str_c(., sffx))
     
     return(readdata)
}

# Read data and primary processing: #######################################
#a. Read depandable data, and bind to a single df -> gdata
#b. Read response data, validate and bind to gdata
####
     files_names=c(
                   "CC_qids.csv",
                    "Enrollment/CRS.csv",
                   "Enrollment/DM.csv",
                   "Enrollment/EL.csv",
                   "Enrollment/HRSD.csv",
                   "Enrollment/MHX.csv",
                   "Enrollment/PDS.csv",
                   "Enrollment/PHX.csv"
     )
     
#collect dependable data
     for(file_counter in c(2:length(files_names)))
     {
          
          mat=read_csv(filename = files_names[file_counter])
        
          if(exists("gdata"))
          {
               gdata %<>% dplyr::inner_join(mat,
                                             by=c("id"))
          }
          else
          {
               gdata=mat
          }
          
          print(stringr::str_c(file_counter, ": ", 
                               files_names[file_counter],
                              " read and binded."))
     }
     
#collect reponse data, and validate 
     # cc_datae=read_csv(filename = files_names[1],
     #                  filepath="/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Sasha Kugel/Clinical Data/exercise helpers/")

     cc_data=read_csv(filename = "CC.csv",
                      filepath=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                                     "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/",
                                     "DataforPatientHistory/DepressionDataNIH/StarD/STAR_DClinicalData/ASCII Files raw data/",
                                     "Level2/Visits/"),
                     withDate = TRUE)
     cc_dataA=read_csv(filename = "CC.csv",
                      filepath=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                                     "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/",
                                     "DataforPatientHistory/DepressionDataNIH/StarD/STAR_DClinicalData/ASCII Files raw data/",
                                     "Level2a/Visits/"),
                      withDate = TRUE)
     cc_dataA$stds5.CC=as.character(cc_dataA$stds5.CC)
     cc_dataA$stds6.CC=as.character(cc_dataA$stds6.CC)
     
     cc_data %<>% bind_rows(cc_dataA)
     
     cc_data %<>% filter_all(any_vars(!is.na(.) & .!=""))
     
     cc_copy=cc_data
     
     {
         #Seeking for gold standard
         #gs= cc_copy %>% group_by(id) %>% filter(date.CC==max(date.CC)) 
         gs= cc_copy %>% group_by(id) %>% filter(qcimp_r.CC==max(qcimp_r.CC)) 
         # cc_data %<>% dplyr::group_by(.data$id) %<>%
         #     dplyr::summarise(ave_imp=mean(.data$qcimp_r.CC, na.rm = T),
         #                      max_imp=max(.data$qcimp_r.CC, na.rm = T)) %<>%
         #     dplyr::ungroup() %<>%
         #     as.data.frame()
         
         #gs %<>% dplyr::mutate(id_name=stringr::str_c("patient_", .data$id)) 
         #gs %<>% filter(id_name %in% rownames(gdata)) %<>% ungroup
         
         gs %<>% dplyr::mutate(response=ifelse(.data$qcimp_r.CC>50, 1,0))
         #gs %<>% tibble::column_to_rownames("id_name") 
         #gs %<>% dplyr::select(-.data$id) 
         
         gs %<>% select(id, response, starts_with("stmd"))
         
         gs %<>% tidyr::gather(key="trtmnt_key", value="code", -id, -response)
         
        citalopram_patients=gs %>% filter(code==103) %>% pull(id) %>% unique
     }
     
     
     cc_data%<>% select(-.data$date.CC)
     # cc_data %<>% dplyr::group_by(.data$id) %<>%
     #                dplyr::summarise(ave_imp=mean(.data$qids.percent.improvement.CC, na.rm = T),
     #                                 max_imp=max(.data$qids.percent.improvement.CC, na.rm = T)) %<>%
     #                dplyr::ungroup() %<>% 
     #                as.data.frame()
     # 
     
     cc_data %<>% dplyr::group_by(.data$id) %<>%
         dplyr::summarise(ave_imp=mean(.data$qcimp_r.CC, na.rm = T),
                          max_imp=max(.data$qcimp_r.CC, na.rm = T)) %<>%
         dplyr::ungroup() %<>%
         as.data.frame()
     
               # Validating calculated fields
               # test_cc=read_csv(filename = "8_CC_qids_calculations.csv")
               # test_cc %<>% dplyr::right_join(cc_data,
               #                                by="id",
               #                                suffix=c(".test", ".cc_data"))
               # test_cc %>% dplyr::filter(round(.data$ave_imp.test,1) !=
               #                              round(.data$ave_imp.cc_data,1)) %>% 
               #              View
               # test_cc %>% dplyr::filter(round(.data$max_imp.test,1) != 
               #                             round(.data$max_imp.cc_data,1)) %>% 
               #              View
               # #test passes; both views are empty
               
#calculate target field: response 
     #1 - responsive patients, QIDS % IMP >50
     #0 - else (unresponsive)
     cc_data %<>% dplyr::mutate(response=ifelse(.data$max_imp>50, 1,0))
     #cc_data %<>% dplyr::mutate(response=ifelse(.data$max_imp>5, 1,0))
 break()    
     gdata = dplyr::inner_join(cc_data,
                               gdata,
                               by="id")
    
     print(stringr::str_c("1: ", 
                          files_names[1],
                          " target feature read and formated."))



# Preprocessing data: #######################################
#A. Impute/remove feature with NAs; data will not hold NAs.
#B. Remove single value features - noninformative.
#C. Explore/impute/rescale negative values
#D. Continuous variables - standartise with z-score.
#E. Dichotomous and categorical values - standartize into 0-1 range.
####
     
 #A. check missing values (per column) and process:
     nas=gdata %>% dplyr::summarise_all(~ sum(is.na(.))) %>% 
                   t %>% 
                    as.data.frame() %>%
                    tibble::rownames_to_column("feature_name")
     colnames(nas)[2]="nas_count"
     nas %<>% dplyr::filter(nas_count > 0)
     
     #1. Spanish imputation:
     # The summed rule: 	if(Spanish DM, HRSD and PHX are NA, then -2,
     #                      else max(Spanish DM, HRSD and PHX, na.rm=T)
     gdata %<>% dplyr::mutate(spanish = ifelse((is.na(.data$spanish.DM) & 
                                                is.na(.data$spanish.HRSD) & 
                                                is.na(.data$spanish.PDS)), 
                                               -2, 
                                               pmax(.data$spanish.DM,
                                                    .data$spanish.HRSD, 
                                                    .data$spanish.PDS, 
                                                    na.rm=T )))
     gdata %<>% dplyr::select(-.data$spanish.DM,
                              -.data$spanish.HRSD, 
                              -.data$spanish.PDS)
     
     #2. Removing ncmdd.EL, ncage.EL and hdtot.HRSD - see rational in WORD document
     gdata %<>% dplyr::select(-.data$ncmdd.EL, 
                              -.data$ncage.EL, 
                              -.data$hdtot.HRSD)
     
     #3. Filling NA in axi_oth.PHX, axii_oth.PHX with -3, 
     #since we don't know why the field was left NA 
     #(no exclusion in other fields can be infered)
     gdata %<>% dplyr::mutate_at(dplyr::vars(.data$axi_oth.PHX, .data$axii_oth.PHX), 
                                 ~ -3)
     #4. Attempt to predict missing values in UNKNOWN (unknown.EL) according to other race columns:
     #white.EL, black.EL, asian.EL, amind.EL, hawai.EL
                    #exploring:
                    # xx=gdata %>% dplyr::select(id, 
                    #                     white.EL, 
                    #                     black.EL, 
                    #                     asian.EL, 
                    #                     amind.EL, 
                    #                     hawai.EL, 
                    #                     unknown.EL) %>%
                    #                mutate(s=white.EL+ black.EL+ asian.EL+ amind.EL+ hawai.EL) 
                    # xx %>% filter(is.na(unknown.EL) & s==0) %>% View   #returns empty set
     #4 CONCLUSION: unknown.EL will be remmoved.
     gdata %<>% dplyr::select(-.data$unknown.EL)   
     
     
#B. Remove feature having single value; non informative.
     single_value_features=gdata %>% dplyr::summarise_all(dplyr::n_distinct) %>% 
                                        t %>% 
                                        as.data.frame() %>% 
                                        tibble::rownames_to_column("feature") %>% 
                                        filter(.data$V1==1) %>%
                                        pull(.data$feature)
     gdata %<>% dplyr::select(-tidyselect::all_of(single_value_features))
       
     
 #C. Explore negative values, remove 
     
     #a. removing 2 patients having QIDS % <0
     gdata %<>% dplyr::filter(.data$max_imp >= 0)
     gdata %<>% dplyr::select(-.data$ave_imp, 
                              -.data$max_imp)
     
     #b. Remove features having over 20% non applicable values 
     #    (otherwise linearity is damaged by negative values)
     # if not stated otherwise, the non-appplicable value is -2
     negative_features<-names(gdata)[sapply(gdata, function(x) min(x))<0]
     
     feature_characteristics=read.csv(stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                                                     "Sasha Kugel/Clinical Data/exercise helpers/",
                                                       "features_charatesitics.csv"),
                                      header = T,
                                      stringsAsFactors = F)
     
     negative_features=feature_characteristics %>% 
                              dplyr::filter(.data$feature %in% negative_features)
     negative_features %<>% dplyr::mutate(nonapplicable_values=
                                               ifelse(is.na(.data$nonapplicable_values) &
                                                        is.na(.data$legit_negative) &
                                                        .data$type!="cont", 
                                                       -2,
                                                       nonapplicable_values))
     fsummary=gdata %>% dplyr::summarise_at(dplyr::vars(negative_features %>% 
                                                             dplyr::pull(.data$feature)),
                                             ~ sum(. <0)) %>% 
                         t %>% 
                         as.data.frame() %>% 
                         tibble::rownames_to_column("feature") %>% 
                         dplyr::rename(cnt_neg=.data$V1) 
     
     negative_features %<>% dplyr::left_join(fsummary,
                                             by="feature")
     
     fsummary=gdata %>% dplyr::summarise_at(dplyr::vars(negative_features %>% 
                                                             dplyr::filter(.data$notes=="negative") %>% 
                                                             dplyr::pull(.data$feature)),
                                            ~sum(.!=-2)) %>% 
                         t %>% 
                         as.data.frame() %>% 
                         tibble::rownames_to_column("feature") %>% 
                         dplyr::rename(cnt_neg_correct_by=.data$V1) 
                    
     negative_features %<>% dplyr::left_join(fsummary,
                                      by="feature")
     
     negative_features %<>% dplyr::mutate(cnt_neg_corrected=ifelse(!is.na(.data$cnt_neg_correct_by),
                                                                   .data$cnt_neg-.data$cnt_neg_correct_by,
                                                                   .data$cnt_neg))
     negative_features %<>% dplyr::mutate(cnt_neg_prcnt=.data$cnt_neg_corrected/nrow(gdata))
     
     negative_features=negative_features %>% dplyr::filter(.data$cnt_neg_prcnt>=0.2) %>% 
                                             dplyr::pull(.data$feature)
     #remove features having over 20% nonapplicable values
     gdata %<>% dplyr::select(-tidyselect::all_of(negative_features))
     
     
     
     #update feature_characteristics df
     feature_characteristics %<>% dplyr::filter(.data$feature %in% colnames(gdata))
     
     

#D. normalize continuous variables with z score
     gdata %<>% dplyr::mutate_at(dplyr::vars(feature_characteristics %>% 
                                    dplyr::filter(.data$type=="cont" & 
                                                       .data$subtype!="identifier") %>% 
                                    dplyr::pull(.data$feature)),
                          ~ (.-mean(.))/sd(.))
     
#E. standartize range of dichotomous and categorical values
     gdata %<>% dplyr::mutate_at(dplyr::vars(feature_characteristics %>% 
                                    dplyr::filter(.data$type!="cont") %>% 
                                    dplyr::pull(.data$feature)),
                          ~ (.-min(.))/(max(.)-min(.)))
     
     print("Data preprocessed.")


# Data exploration: #######################################
#Simple histogram overview

     # pdf(stringr::str_c("C:/Users/owner/Dropbox/",
     #                        "Genetika+ExerciseConfidentialFile/",
     #                        "plots/histograms.pdf"))
     # for(feature in colnames(gdata)[-1])
     # {
     #      p=hist(gdata[,feature], 
     #             main=stringr::str_c("Histogram for ", feature),
     #             xlab=feature,
     #             labels = TRUE)
     #      print(p)
     # }
     # dev.off()


# Modeling #######################################
# Feature selection and classification/prediction modeling.
#A. Partition the data and prep
#B. Use LASSO model to perform primary classification (on TRAIN)
#C. Identify informative features using LASSO coefficients 
#D. Train xgboost model using features from C (on TRAIN)
#E. Predict on TRAIN, TEST and report.
####
#A. Partition the data and prep
     #filter only citalopram patients
     gdata %<>% filter(id %in% citalopram_patients)
     
     gdata %<>% dplyr::mutate(id_name=stringr::str_c("patient_", .data$id)) 
     gdata %<>% tibble::column_to_rownames("id_name") 
     gdata %<>% dplyr::select(-.data$id) 
     
     response_vector=gdata %>% dplyr::pull(.data$response)
     names(response_vector) = rownames(gdata)
          
     #Partition the data into to trainining and testing sets, 
     #balanced towards the smaller group
     trainIndex <- caret::createDataPartition(response_vector, 
                                              p = .75, 
                                              list = FALSE, 
                                              times = 1)
     numToBalance=min(table(response_vector[trainIndex]))
     
     balancedTrainIndex=data.frame(inx=trainIndex,
                                    resp=response_vector[trainIndex]) %>%
                         dplyr::filter(.data$resp==1) %>%
                         dplyr::pull(.data$Resample1) %>%
                         sample(numToBalance)
     balancedTrainIndex=c(balancedTrainIndex,
                          data.frame(inx=trainIndex,
                                     resp=response_vector[trainIndex]) %>%
                              dplyr::filter(.data$resp==0) %>%
                              dplyr::pull(.data$Resample1)
                         )
     

#B. Use LASSO model to perform primary classification (on TRAIN)
     gdata.modeling=model.matrix(response~. , gdata)[,-1]
     
     cv_output <- cv.glmnet(gdata.modeling[balancedTrainIndex,], 
                            response_vector[balancedTrainIndex],
                            alpha = 1, 
                            family="binomial",
                            nfolds = 10)
     
     best_lam=cv_output$lambda.min
     
     lasso_best <- glmnet(gdata.modeling[balancedTrainIndex,], 
                          response_vector[balancedTrainIndex], 
                          alpha = 1, 
                          lambda = best_lam,
                          family="binomial"
                          )
     
     pred <- predict(lasso_best, 
                     s = best_lam, 
                     newx = gdata.modeling,
                     type="class")
     
     pred.df=data.frame(res.original=response_vector,
                        res.predicted=as.numeric(pred),
                        stringsAsFactors = F)
     
     print(confusionMatrix(table(pred.df[balancedTrainIndex,])))
     
#C. Identify informative features using LASSO coefficients 
     cfcnt=coef(lasso_best,s=best_lam ,exact=TRUE) %>% as.matrix() 
     cfcnt %<>% as.data.frame() %<>% 
          tibble::rownames_to_column("feature")
     colnames(cfcnt)[2]="prob"
     cfcnt %<>% dplyr::filter(.data$prob!=0)
     
     subgdata=gdata %>% dplyr::select(.data$response, 
                                      tidyselect::any_of(cfcnt %>% 
                                                         dplyr::pull(.data$feature)))
     
     subgdata.modeling=model.matrix(response~. , subgdata)[,-1]


#D. Train xgboost model using features from C (on TRAIN)

               # #Get best number of rounds for the task, based on TRAIN data
               # params <- list(booster = "gbtree", 
               #                objective = "binary:logistic", 
               #                eta=0.1, 
               #                gamma=0, 
               #                max_depth=6, 
               #                min_child_weight=1, 
               #                subsample=1, 
               #                colsample_bytree=1)
               # 
               # dtrain=xgb.DMatrix(data = subgdata.modeling[balancedTrainIndex,],
               #                    label = response_vector[balancedTrainIndex])
               # 
               # xgbcv <- xgb.cv( params = params,
               #                     data = dtrain,
               #                     nrounds = 100,
               #                     nfold = 6,
               #                     showsd = T,
               #                     stratified = T,
               #                     print_every_n =  10,
               #                     early_stopping_rounds = 20,
               #                     maximize = F)

     #Since we know that the testing set is highly skewed (unbalanced), 
     #we will train the model by considering the nonresponsive group be less represented
     #Despite the balanced TRAIN set.
     x_scale_pos_weight = round(sqrt(max(table(response_vector[-balancedTrainIndex])) / 
                                          min(table(response_vector[-balancedTrainIndex]))))
     
     bstSparse <- xgboost(data = subgdata.modeling[balancedTrainIndex,],
                          label = response_vector[balancedTrainIndex],
                          max.depth = 4,
                          eta = 0.1,
                          nthread = 3,
                          nrounds = 20 , #Why 63? xgb.cv showed several times for it to be best fit
                          scale_pos_weight = x_scale_pos_weight,
                          eval_metric = "error",
                          objective = "binary:logistic")
 
#E. Predict on TRAIN, TEST and report.    
     xgb.pred <- predict(bstSparse, subgdata.modeling);
     
     xgb.pred.df=data.frame(res.original=response_vector,
                            res.predicted=as.numeric(xgb.pred),
                            stringsAsFactors = F);
     
     xgb.pred.df %<>% dplyr::mutate(res.predicted.binary=ifelse(res.predicted>=0.5, 1, 0));
     
     print(confusionMatrix(table(xgb.pred.df[balancedTrainIndex,-2])));
     print(confusionMatrix(table(xgb.pred.df[-balancedTrainIndex,-2])))
     
     

     {
        #gold standard experts opinion
         #gs %<>% filter(id %in% citalopram_patients)
         gs %<>% mutate(is_cit=ifelse(id %in% citalopram_patients, 1, 0))
         gs %<>% ungroup %>%  select(id, response, is_cit) %>% distinct 
         
         print("Experts opinion")
         print(confusionMatrix(table(gs %>% select(-id))))
     }



