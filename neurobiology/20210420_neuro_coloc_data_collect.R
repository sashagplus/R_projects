library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library("rpart")
library(rpart.plot)
library(readxl)    
library(caret)

library(gpio)

rm(list=ls())

# read_excel_allsheets <- function(filename, tibble = FALSE) {
#   # I prefer straight data.frames
#   # but if you like tidyverse tibbles (the default with read_excel)
#   # then just pass tibble = TRUE
#   sheets <- readxl::excel_sheets(filename)
#   x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
#   if(!tibble) x <- lapply(x, as.data.frame)
#   names(x) <- sheets
#   x
# }
# 

data_folder=stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
                           "01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/04_Differentiation_Experiments/")


files_location=read.csv(stringr::str_c(data_folder,
                                       "20210228_File_Location_well_index_per_line_including2019Data.csv"))

files_location %<>% select(DF, Line, Well, Treat, Days, BUP.Responder, folder.coloc, 
                           Coloc.file...Image_n_ROI_Length...xls, Coloc.file...PSD..csv, Coloc.file...SYN..csv) %>% distinct() 

colnames(files_location) = c( "DF", "Line", "Well", "Treat", "Days", "BUP.Responder",
                              "folder.coloc", "image_nROI_Length", "PSD", "SYN")

files_location %<>% mutate(folder.coloc=str_replace_all(str_replace(folder.coloc, "C:", ""), "\\\\", "/"))
files_location %<>% mutate(folder.coloc=str_replace(folder.coloc, "Owner/", "sashakugel/gplus_dropbox/"))

#TEMP Fix
files_location %<>% mutate(folder.coloc=ifelse(folder.coloc=="/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/04_Differentiation_Experiments/2020_07_Df63/Coloc/Diff63_L23_Image.xlsx
", "/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/04_Differentiation_Experiments/2020_07_Df63/Coloc/Diff63_L23_Image.xlsx
/DF63_L23", folder.coloc))

files_location %<>% mutate(roi_filepath=str_c(folder.coloc, "/", image_nROI_Length, ".xlsx"),
                           psd_filepath=str_c(folder.coloc, "/", PSD, ".csv"),
                           syn_filepath=str_c(folder.coloc, "/", SYN, ".csv"))

files_location %<>% select(-folder.coloc, -image_nROI_Length, -PSD, -SYN)

#TEMP: remove df==62 with line 10, wrong well labeling. Line 10 also repeats in DF71
files_location %<>% filter(DF!=62)



vehicles=files_location %>% filter(Treat=="VEH")
treatments=files_location %>% filter(Treat!="VEH")

vehicles %<>% mutate(Days=as.character(Days)) %>% 
  mutate(Days=ifelse(Days=="73", "7 3", Days)) %>% 
  separate(Days, into=c("seven", "three"), sep=" ", remove = F, fill="left") %>% 
  gather(key, value, -DF, -Line, -Well, -Treat, -Days, -BUP.Responder, -roi_filepath, -psd_filepath, -syn_filepath, na.rm = T) %>% 
  mutate(Days=as.integer(value)) %>% select(-key, -value)

vehicles %<>% select(-Well) %<>% distinct()
treatments %<>% select(-Well) %<>% distinct()

files_treatToVehicle=treatments %>% left_join(vehicles, 
                          by=c("Line", "DF", "Days", "BUP.Responder"), 
                          suffix = c(".treatment", ".vehicle")) %>%
  as.data.frame() 

# write.csv(files_treatToVehicle,
#          str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
#          "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science",
#          "/Exploratory_Outputs/neurobiology_imaging/Coloc/",
#          "20210425_coloc_treatment_to_vehicle_files_mapping.csv"),
#          row.names = F)




all_roi_files=list()
df_cnames=NULL;
distinct_roi_file_names=files_location %>% select(roi_filepath, psd_filepath, syn_filepath) %>% distinct()
rois=data.frame();
psds=data.frame();
syns=data.frame();
for (rownum in 1:nrow(distinct_roi_file_names))
{
  print(rownum)
  
  roi_filename=distinct_roi_file_names$roi_filepath[rownum]
  if(file.exists(roi_filename))
  {
    roi_file=gpio::read_excel_allsheets(filename = roi_filename)
    
    all_roi_files[[rownum]] = roi_file
    
    sheetname = names(roi_file)[grep("diff", tolower(names(roi_file)))]
    
    roi_datasheet=roi_file[[sheetname]] 
    
    if (!"Dendrite length pix" %in% colnames(roi_datasheet))
    {
      missing_column=(roi_datasheet %>% select(starts_with("...")) %>% colnames)[1]
      roi_datasheet %<>% rename(`Dendrite length pix`=!!missing_column) 
    }
    
    
    
    if(is.null(df_cnames)) 
    {
      df_cnames=data.frame(cnames=c(colnames(roi_datasheet), "tester"), second=1)
      colnames(df_cnames)[2]="file_1"
    }
    else
    {
      df_cnames %<>% left_join(data.frame(cnames=colnames(roi_datasheet), second=1), by=c("cnames"))
      colnames(df_cnames)[ncol(df_cnames)]=str_c("file_", rownum)
    }
    roi_datasheet %<>% mutate(`Dendrite length pix`=as.numeric(ifelse(`Dendrite length pix`=="-", NA, `Dendrite length pix`)))
    
    roi_datasheet %<>% mutate(roi_filename=roi_filename) %>%
                      relocate(roi_filename)
    rois %<>% bind_rows(roi_datasheet)
  }
  else
  {
    print(str_c(roi_filename, "\n\n", "IS MISSING"))
  }
    
  #PSD
  psd_filename=distinct_roi_file_names$psd_filepath[rownum]
  if(file.exists(psd_filename))
  {
    #psd_file=read_excel_allsheets(filename = psd_filename)
    psd_datasheet=read.csv(psd_filename)
    
    #sheetname = names(psd_file)[grep("diff", tolower(names(psd_file)))]
    
    #psd_datasheet=psd_file[[sheetname]] 
    
    psd_datasheet %<>% mutate(psd_filename=psd_filename,
                              roi_filename=roi_filename) %>%
                      relocate(psd_filename, roi_filename)
    psds %<>% bind_rows(psd_datasheet)
  }
  
  #PSD
  syn_filename=distinct_roi_file_names$syn_filepath[rownum]
  if(file.exists(syn_filename))
  {
    #syn_file=read_excel_allsheets(filename = syn_filename)
    syn_datasheet=read.csv(syn_filename)
    #sheetname = names(psd_file)[grep("diff", tolower(names(psd_file)))]
    
    #psd_datasheet=psd_file[[sheetname]] 
    
    syn_datasheet %<>% mutate(syn_filename=syn_filename,
                              roi_filename=roi_filename) %>%
                      relocate( syn_filename, roi_filename)
    syns %<>% bind_rows(syn_datasheet)
  }
}

break()
rois %<>% rowwise() %>%  mutate(Well=str_remove(str_split(.data$FileName_PSD95, "_")[[1]][2], "well"),
                               Field=str_split(.data$FileName_PSD95, "_")[[1]][3], .after=roi_filename) 


#join rois with files_location
rois_to_files_location=files_location %>% right_join(rois, by=c("roi_filepath"="roi_filename", "Well")) 

#which rois didnt have a match in files_location
rois_not_in_files=rois_to_files_location %>% filter(is.na(DF))

#which df+line +well didn't appear in rois
files_not_in_rois=files_location %>% anti_join(rois_to_files_location, by=c("DF", "Line", "Well", "roi_filepath")) 

write.csv(files_not_in_rois,
stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/",
"01_Protocol_Development/01_05_Imaging/01_Experiments_and_Results/04_Differentiation_Experiments/20210428_missing_wells_rois_sasha.csv")
, row.names = F)


#Inspecting features
feature_inspect=melt(lapply(rois, class)) %>%  
                  select(feature_name=L1, feature_type=value) %>% 
                  mutate(to_model=ifelse(feature_type=="character", "No", "Yes")) %>%
                  mutate(file_type="ROI", .before = feature_name)

feature_inspect %<>% bind_rows(melt(lapply(psds, class)) %>%  
                               select(feature_name=L1, feature_type=value) %>% 
                                mutate(to_model=ifelse(feature_type=="character", "No", "Yes")) %>%
                                mutate(file_type="PSD", .before = feature_name)
                        )

feature_inspect %<>% bind_rows(melt(lapply(syns, class)) %>%  
                                 select(feature_name=L1, feature_type=value) %>% 
                                 mutate(to_model=ifelse(feature_type=="character", "No", "Yes")) %>%
                                 mutate(file_type="SYN", .before = feature_name)
                      )

identifier_fields=c("folder.coloc", "image_nROI_Length", "roi_filename", "Well", "Field", "Group_Index", "ImageNumber")

feature_inspect %<>% mutate(is_identifier=ifelse(feature_name %in% identifier_fields, TRUE, FALSE))

feature_inspect_backup=feature_inspect
#Sari went over the features, and now we can exclude those that are not necessary  
feature_inspect=read.csv(str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                               "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/",
                               "Exploratory_Outputs/neurobiology_imaging/Coloc/",
                               "20210421_coloc_feature_inspection_for_modeling_SN.csv"))

feature_inspect %<>% filter(!feature_name %in% c("folder.coloc", "image_nROI_Length", "PSD", "SYN")) 


feature_inspect %<>% filter(to_model=="Yes" | is_identifier==TRUE) 

#select only modelable features
rois %<>% select(one_of(feature_inspect %>% filter(file_type=="ROI") %>% pull(feature_name))) 
psds %<>% select(one_of(feature_inspect %>% filter(file_type=="PSD") %>% pull(feature_name))) 
syns %<>% select(one_of(feature_inspect %>% filter(file_type=="SYN") %>% pull(feature_name))) 


rois_backup=rois

#match data to classifiers
rois %<>% left_join(files_location %>% 
                     select(-psd_filepath, -syn_filepath) %>%
                    rename(roi_filename=roi_filepath), by=c("roi_filename", "Well")) %>% 
          relocate(files_location %>% 
                     select(-psd_filepath, -syn_filepath) %>%
                     rename(roi_filename=roi_filepath) %>% colnames()) 

print("writing output")
write.csv(rois,
          str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science",
                "/Exploratory_Outputs/neurobiology_imaging/Coloc/",
                "20210425_coloc_rois_collected_data.csv"),
          row.names = F)


#field will be considered replicate for now

break()
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
set.seed(1)
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
                                           test=cm.test$overall["Accuracy"],
                                           test.sensetivity=cm.test$byClass[["Sensitivity"]],
                                           test.specificity=cm.test$byClass[["Specificity"]]))
  # accuracies.xgb %<>% bind_rows(data.frame(nn=i,
  #                                          train=cm.train$overall["Accuracy"], 
  #                                          test=cm.test))
}
#prms$accuracy.train[p]=mean(accuracies.xgb$train)
#prms$accuracy.test[p]=mean(accuracies.xgb$test)
#}

print(str_c("mean accuracy:",  mean(accuracies.xgb$test)))
print(str_c("mean sensetivity:",  mean(accuracies.xgb$test.sensetivity)))
print(str_c("mean specificity:",  mean(accuracies.xgb$test.specificity)))


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



