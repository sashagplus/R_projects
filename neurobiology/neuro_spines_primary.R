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


data_folder=stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                           "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/04_Differentiation_Experiments/")


files_location=read.csv(stringr::str_c(data_folder,
                                       "20210228_File_Location_well_index_per_line.csv"))


#bupoprion only + vehicle
files_location %<>% filter(Treat %in% c("BUP", "VEH"))

# 
# 1. replace in files_location all \\ with / and remove prefix C:\\
# 2. for every spine folder collect all spine files, read them and count how many sheets are in those files 
files_location %<>% mutate(folder.spines=str_replace_all(str_replace(folder.spines, "C:", ""), "\\\\", "/"))
files_location %<>% mutate(folder.spines=str_replace(folder.spines, "Owner/", "sashakugel/gplus_dropbox/"))

all_spine_files=list()
spine_files=data.frame();
for (well in 1:nrow(files_location))
{
  well_files_list=list.files(str_c(files_location$folder.spines[well], "/"), pattern = str_c("(", files_location$Well[well] ,").*\\.xlsx$"))
  if(length(well_files_list)==0) {well_files_list=NA}
  well_files=data.frame(Well=files_location$Well[well],
                        folder.spines=files_location$folder.spines[well],
                        filenames=well_files_list)

  spine_files %<>% bind_rows(well_files)
}

files_location %<>% left_join(spine_files, by=c("Well", "folder.spines"))
files_location %<>% select(DF, Line, Well, Treat, Days, BUP.Responder, filenames, everything())

write.csv(files_location,
          stringr::str_c(data_folder,
                         "20210228_File_Location_well_index_per_line_SPINES.csv"),
          quote = TRUE,
          row.names = FALSE)


#Collect data for dendrites (sometimes more then 1 in every file)
#Transpose all avaiable data to have a single enrty for every dendrtie
#For spines: pick only numeric columns, and summarise them.
# dndrt_data=data.frame()
# for (field in 1:nrow(files_location))
# {
#   if(!is.na(files_location$filenames[field]))
#   {  
#     spines_file=read_excel_allsheets(filename = str_c(files_location$folder.spines[field], "/",
#                                                     spine_files$filenames[field]))
#     
#     #Trees
#     treespines=spines_file$`Tree Spines - Dendrites` %>% select(Type, starts_with("Density"))
#     colnames(treespines)[-1]=c(1:(ncol(treespines)-1))
#       
#     treespines %<>% gather(key="Centrifugal", value="Density(1/µm)", -Type) %>%
#                   mutate(Type.Density=str_c(Type, ".Density")) %>%  
#                   select(-Type)  %>% 
#                   spread(key=Type.Density, value = `Density(1/µm)`) %>% 
#                   mutate(Centrifugal=as.numeric(Centrifugal))
#     
#     
#     
#     
#     dendrite=data.frame(c(files_location[field,1:8]), 
#                         spines_file$`Segment - Dendrites` %>% mutate(`Z Angle`=as.numeric(`Z Angle`)))
#     
#     dendrite %<>% left_join(treespines, by=c("Centrifugal"))
#    
#     spines=spines_file$`Spine Details - Automatic (Back`
#     spines=spines[1:(nrow(spines)-2),]
#       
#     spines %<>% group_by(Tree, `Branch Order`, `Spine Type`)
#       spines %<>% select_if(is.numeric)
#       
#       spines %<>% group_by(Tree, `Branch Order`, `Spine Type`) %>%
#                   summarise(across(everything(), 
#                               list(mean = mean, min=min, max=max, median=median), 
#                               .names = "{.col}.spines.{.fn}"), .groups = "keep") 
#       
#       spines %<>% ungroup () %>% 
#                 gather(key="key", value="value", -Tree, -`Branch Order`, -`Spine Type`) %>% 
#                 mutate(new_key=str_c(`Spine Type`, ".", key)) %>% 
#                 select(-key, -`Spine Type`) %>% 
#                 spread(key="new_key", value="value") 
#       
#       dendrite %<>% left_join(spines, by=c("Tree"="Tree", "Centrifugal"="Branch Order"))
#     
#       dndrt_data %<>% bind_rows(dendrite)
#     }
# }
# dndrt_data1=dndrt_data
# #for now, taking only Tree=1 and Centrifugal =1 
# dndrt_data %<>% filter(  Tree==1 & Centrifugal==1) 
# 
# 
# dndrt_data %<>% select(-Tree, -Centrifugal)
# 
# dndrt_data %<>% filter(Treat=="BUP" & Days==7)
# dndrt_data %<>% group_by(Line, BUP.Responder) %>% filter(DF==max(DF)) 
# 
# output_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
#                       "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
#                       "04_Differentiation_Experiments/Data Science/Spines_modeling/",
#                       "20210314_spines_Bup_7d_maxDF_Tree1_Centrifugal1_preprocessed_data.csv")
# 
# write.csv(dndrt_data, output_filename,
#           quote = TRUE,
#           row.names = FALSE)


processed_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                      "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/",
                      "04_Differentiation_Experiments/Data Science/Spines_modeling/",
                      "20210314_spines_Bup_7d_maxDF_Tree1_Centrifugal1_preprocessed_data.csv")

dndrt_data=read.csv(processed_filename)


dndrt_data %<>% mutate(name=str_c(DF, "_", Line, "_", Well, "_", Days, "_", 
                                  str_remove(filenames, ".xlsx")))


dndrt_data %<>% select(-DF, -Line, -Well, -Days, -filenames) %>%
                select(name, everything())

dndrt_data_info=dndrt_data %>% select(name, Treat, BUP.Responder, folder.spines)

dndrt_data %<>% select(-folder.spines)
dndrt_data %<>% tibble::column_to_rownames ("name")


dndrt_data %<>%
  select_if(~ !any(is.na(.)))

dndrt_data = dndrt_data %>% select ( BUP.Responder) %>% 
            bind_cols(dndrt_data %>% select_if(~is.numeric(.)))



dndrt_data %<>% mutate(BUP.Responder=ifelse(BUP.Responder=="NR", 0,1))

# accuracies=data.frame()
# for(i in 1:20)
# {  
#   train=c(sample(which(dndrt_data$BUP.Responder==0), 40*0.875) ,
#           sample(which(dndrt_data$BUP.Responder==1), 32*0.875) )
#   test=c(1:nrow(dndrt_data))[-train]
#   
#   mylogit <- glm(BUP.Responder ~ ., data = dndrt_data[train,], family = "binomial")
#   
#   predicted=round(predict(mylogit, newdata = dndrt_data[,-1], type = "response"))
#   
#   cm.train <- confusionMatrix(factor(dndrt_data$BUP.Responder[train]), 
#                               factor(predicted[train]))
#   print("Train")
#   print(cm.train$overall)
#   
#   cm.test <- confusionMatrix(factor(dndrt_data$BUP.Responder[test]), 
#                               factor(predicted[test]))
#   
#   accuracies %<>% bind_rows(data.frame(nn=i,
#                                       train=cm.train$overall["Accuracy"], 
#                                        test=cm.test$overall["Accuracy"]))
#   }
# 
# 
# print(1)
# bupresponder=ifelse(dndrt_data$BUP.Responder==0,
#                      "NR", "R")
# 
# bupresponder.color=ifelse(dndrt_data$BUP.Responder==0,
#                     "blue", "red")
# treat.color=ifelse(dndrt_data_info$BUP.Responder==0,
#                    "blue", "red")
# 
# output_folder="/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/Exploratory_Outputs/neurobiology_imaging/Spines/"

# pdf(str_c(output_folder, "11032021_spines_boxplots.pdf"))
# 
# for(i in colnames(dndrt_data)[-1]) {
#   
#   boxplot(dndrt_data[,i]~bupresponder, ylab = i, xlab = "BUP.Responder")
#   
#   
# }
# dev.off()


# dndrt_data.pca <- prcomp(dndrt_data[,-1])
# summary(dndrt_data.pca)
# 
# bupresponder.num=dndrt_data$BUP.Responder
# dndrt_data$BUP.Responder=bupresponder
#autoplot(prcomp(dndrt_data[,-1]), data=dndrt_data, col="BUP.Responder")

#library(rgl)
#plot3d(dndrt_data.pca$x[,1], dndrt_data.pca$x[,2], dndrt_data.pca$x[,3], col = rainbow(1000))


# fit <- rpart(BUP.Responder ~ ., data = dndrt_data[train,], method = 'class')
# rpart.plot(fit, extra = 106)
# 
# predicted=predict(fit, dndrt_data[test,], type = 'class')
# print(2)
# accuracies.dt=data.frame()
# for(i in 1:20)
# {  
#   train=c(sample(which(dndrt_data$BUP.Responder==0), 40*0.875) ,
#           sample(which(dndrt_data$BUP.Responder==1), 32*0.875) )
#   test=c(1:nrow(dndrt_data))[-train]
# 
#   dt.fit <- rpart(BUP.Responder ~ ., data = dndrt_data[train,], method = 'class')
#   
#   predicted=predict(dt.fit, dndrt_data, type = 'class')
#   
#   cm.train <- confusionMatrix(factor(dndrt_data$BUP.Responder[train]), 
#                               factor(predicted[train]))
#   print("Train")
#   print(cm.train$overall)
#   
#   cm.test <- confusionMatrix(factor(dndrt_data$BUP.Responder[test]), 
#                              factor(predicted[test]))
#   
#   accuracies.dt %<>% bind_rows(data.frame(nn=i,
#                                        train=cm.train$overall["Accuracy"], 
#                                        test=cm.test$overall["Accuracy"]))
#   
# }


output_filename_predata=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                              "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
                              "04_Differentiation_Experiments/Data Science/Spines_modeling/",
                              "20210315_spines_Bup_7d_maxDF_premodeling_data.csv")



write.csv(bind_cols(dndrt_data_info %>% select(-BUP.Responder, -folder.spines, -Treat), dndrt_data), 
          output_filename_predata,
          quote = TRUE,
          row.names = FALSE)




library(xgboost)

#set.seed(0)
set.seed(NULL)
# prms=merge(seq(1,10,1), seq(0.1,0.9, 0.1))
# colnames(prms)=c("max.depth", "eta")
# prms %<>% full_join(data.frame(nthreads=c(1:10)), by=as.character())
# prms %<>% full_join(data.frame(nrounds=c(1:30)), by=as.character())
# prms$accuracy.train=NA
# prms$accuracy.test=NA

#mean.accuracies.prms.xgb=data.frame()
# for(p in 21091:nrow(prms))
# {
  #if(p%%1000==0) print(p)
  
  accuracies.xgb=data.frame()
  predictions=data.frame(response=dndrt_data$BUP.Responder)
  for(i in 1:20)
  {  
    train=c(sample(which(dndrt_data$BUP.Responder==0), 40*0.875) ,
            sample(which(dndrt_data$BUP.Responder==1), 32*0.875) )
    test=c(1:nrow(dndrt_data))[-train]
  
    # train=c(1:nrow(dndrt_data))[-1]
    # test=c(i)
    
    xgb.model <- xgboost(data = as.matrix(dndrt_data[train,-1]),
                         label = dndrt_data[train,1],
                         max.depth = 2, #prms$max.depth[p],
                         eta = 0.2, #prms$eta[p],
                         nthread = 5, #prms$nthreads[p],
                         nrounds = 30, #prms$nrounds[p] , #Why 63? xgb.cv showed several times for it to be best fit
                         eval_metric = "error",
                         objective = "binary:logistic",
                         verbose = 0)
    
    predicted <- predict(xgb.model, as.matrix(dndrt_data[,-1]));
    
    pp=data.frame(predicted)
    colnames(pp)=str_c("prediction_", ncol(predictions))
    
    predictions %<>% bind_cols(pp)
    
    predicted <- round(predict(xgb.model, as.matrix(dndrt_data[,-1])));
    
    cm.train <- confusionMatrix(factor(dndrt_data$BUP.Responder[train], levels = unique(dndrt_data$BUP.Responder)), 
                                factor(predicted[train],levels = unique(dndrt_data$BUP.Responder)))
    #print("Train")
    #print(cm.train$overall)
    
    cm.test <- cm.test <- confusionMatrix(factor(dndrt_data$BUP.Responder[test], levels = unique(dndrt_data$BUP.Responder)),
                                          factor(predicted[test], levels = unique(dndrt_data$BUP.Responder)))
  
    
  #  cm.test <- ifelse(predicted[test]==dndrt_data$BUP.Responder[test], 100, 0);
    
    accuracies.xgb %<>% bind_rows(data.frame(nn=i,
                                            train=cm.train$overall["Accuracy"],
                                            test=cm.test$overall["Accuracy"]))
    # accuracies.xgb %<>% bind_rows(data.frame(nn=i,
    #                                          train=cm.train$overall["Accuracy"], 
    #                                          test=cm.test))
  }
  #prms$accuracy.train[p]=mean(accuracies.xgb$train)
  #prms$accuracy.test[p]=mean(accuracies.xgb$test)
#}

print(mean(accuracies.xgb$test))

predictions =bind_cols(dndrt_data_info, predictions)

predictions %<>% separate(name, into=c("DF", "Line", "Well", "Day", "Well_filename", "Field_filename")) %<>% 
  unite("filename", Well_filename, Field_filename, sep = "_")

predictions %<>% group_by_at(vars(!starts_with("prediction_"))) %>% 
                  rowwise() %>%
                  mutate(mean_prediction=mean(c_across( starts_with("prediction_"))))

output_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                      "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/", 
                      "04_Differentiation_Experiments/Data Science/Spines_modeling/",
                      "20210314_spines_Bup_7d_maxDF_modeling_predictions.csv")

write.csv(predictions, output_filename,
          quote = TRUE,
          row.names = FALSE)

  
# parameters_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
#                       "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/",
#                       "04_Differentiation_Experiments/Data Science/Spines_modeling/",
#                       "spines_xgb_params_tweaking.csv")
# 
# write.csv(prms, parameters_filename,
#           quote = TRUE,
#           row.names = FALSE)



