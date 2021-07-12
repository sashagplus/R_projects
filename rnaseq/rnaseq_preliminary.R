#setwd("~/GenetikaPlus_Projects_R/rnaseq")
library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library("rpart")
library(rpart.plot)

rm(list=ls())

data_folder=stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_06_RNA-seq/Analysis/Analysis full/")

data_files=data.frame(type=c("exp", "grp", "exp", "grp", "exp", "grp"),
                      filename=c("Final_tables/2021_02_00/2021_02_22_7d_Bup_cpm_exclude_L38.csv",
                                 "2021_02_11_newbatch_groups_ex_L38.csv",
                                 "cit_exp_2021_02_22_7d_Citalopram_cpm.csv",
                                 "cit_grp_2021_02_14_7d_citalopram_groups.csv",
                                 "mirt_exp_2021_02_22_7d_Mirt_cpm.csv",
                                 "mirt_grp_2021_02_14_7d_mirtazapine_groups.csv"))

#read data
bup_grp=read.csv(file=stringr::str_c(data_folder,data_files$filename[2]))
bup_grp %<>% distinct

bup_exp=read.csv(file=stringr::str_c(data_folder,data_files$filename[1]))

sample_profiles=data.frame(lname=colnames(bup_exp)[-c(1:3)], exp_col_index=(1:ncol(bup_exp))[-c(1:3)] )
sample_profiles %<>% separate(col=.data$lname, into=c("sample", "rest"), sep="\\.", extra = "merge", remove=F)
sample_profiles %<>% mutate(treatment=case_when(grepl("bup", rest, ignore.case = TRUE) ~ "Bup",
                                                grepl("Treated", rest, ignore.case = TRUE) ~"Bup",
                                                TRUE ~ "Vehicle"))


sample_profiles %<>% dplyr::left_join(bup_grp, 
                                      by=c("sample", "treatment"))
sample_profiles %<>% mutate(nname=stringr::str_c(sample, "_",treatment))
sample_profiles %<>% group_by(sample, treatment ) 

bup_exp %<>% gather(key="lname", value="cnts", -ensemblID, -gene_ID, -data)

bup_exp %<>% left_join(sample_profiles %>% select(lname, sample, treatment, group, nname),
                       by=c("lname"))
bup_exp %<>% group_by(ensemblID, gene_ID, data, sample, treatment, group, nname) %<>%
  summarise(med_cnt=median(cnts), .groups = "keep")

bup_exp %<>% ungroup() %<>% select(-nname) %<>% spread(key=treatment, value = med_cnt) %<>% as.data.frame()

bup_exp %<>%rowwise() %<>% mutate(fc=ifelse(Vehicle==0,
                                            ifelse(Bup==0, 1, Bup),
                                            Bup/Vehicle))

bup_exp %<>% mutate(logfc=log(fc,2))

cut_bup_exp = bup_exp%>% filter(logfc<=log(0.7,2) | logfc >=log(1.3,2))


mat=bup_exp %>% ungroup %>%select(ensemblID, sample, group, logfc) %>% spread (key=ensemblID, value=logfc)


#mat=cut_bup_exp %>% ungroup %>%select(ensemblID, sample, group, logfc) %>% spread (key=ensemblID, value=logfc)

mat %<>%
  mutate(across(starts_with("ENS"), ~ ifelse(is.infinite(.x), min(is.finite(.x))*1.5, .x)))

mat %<>% as.data.frame() %<>% tibble::column_to_rownames("sample")


ttests=mat %>%
  summarise_each(funs(t.test(.[group == "Non-Res"], .[group == "Responders"])$p.value), vars = -group)

ttests=as.data.frame(t(ttests))
ttests %<>% tibble::rownames_to_column("ensembleID")
colnames(ttests)[2]="pval"

# nas=colSums(is.na(mat))
# submat=mat[,names(which(nas==0))] %>% as.data.frame()
# for ( i in 3:ncol(submat)) {submat[is.infinite(submat[,i]),i]=min(submat[is.finite(submat[,i]), i])*1.5}

sig_features=ttests %>% filter(pval<=0.01) %>% pull(ensembleID)
#mat %<>% select(group, !!sig_features)

# genes=bup_exp %>% select(ensemblID, gene_ID, data) %>% distinct
# genes %<>% mutate(labels=stringr::str_c(ensemblID, "  ", gene_ID, " (", data, ")"))
# 
# pdf("boxplots.pdf")
# for(i in sig_features) {
#   ll= genes %>% filter(ensemblID==i) %>% pull(labels)
#   boxplot(mat[,i]~mat$group, ylab = ll)
#
# }
# dev.off()

output_filename_processed=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                              "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/rnaseq/primary_models/",
                              "20210315_rnaseq_Bup_7d_processed.Rdata")

 save(mat, file=output_filename_processed)

 
 set.seed(0)
load(output_filename_processed)

#break()
 
library(caret)
library(xgboost)
# balance partition
mat.group=mat$group
resp_partition = split(sample(which(mat.group=="Responders")),1:8)
nonresp_partition = split(sample(which(mat.group=="Non-Res")),1:8)

mat %<>% mutate(group=ifelse(group=="Responders", 1, 0))

mat1=mat


# output_filename_predata=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
#                               "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/rnaseq/primary_models/",
#                               "20210315_rnaseq_Bup_7d_premodeling_data.csv")
# 
# write.csv(mat %>% tibble::rownames_to_column("Line"), output_filename_predata,
#           quote = TRUE,
#           row.names = FALSE)



#set.seed(0)
#set.seed(NULL)
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
mat.group=mat$group
#resp_partition = split(sample(which(mat.group==1)),1:4)
#nonresp_partition = split(sample(which(mat.group==0)),1:4)
accuracies.xgb=data.frame()
predictions=data.frame(Line=row.names(mat) ,response=mat$group)
models=list()

for(i in 1:8)
{  
  
  exclude_samples=c(resp_partition[[i]], nonresp_partition[[i]])
  
  xgb.model <- xgboost(data = as.matrix(mat[-exclude_samples,sig_features]),
                       label = mat[-exclude_samples,1],
                       max.depth = 2, #prms$max.depth[p],
                       eta = 0.2, #prms$eta[p],
                       nthread = 5, #prms$nthreads[p],
                       nrounds = 30, #prms$nrounds[p] , #Why 63? xgb.cv showed several times for it to be best fit
                       eval_metric = "error",
                       objective = "binary:logistic",
                       verbose = 0)

  
  #mylogit <- glm(group ~ ., data = mat[-exclude_samples, c("group",sig_features)], family = "binomial")
  
  #predicted <- predict(xgb.model, as.matrix(mat[,-1]));
  predicted <- round(predict(xgb.model, as.matrix(mat[,sig_features])));
  
  
  #predicted=round(predict(mylogit, newdata = mat[,sig_features], type = "response"))
  
  pp=data.frame(predicted)
  colnames(pp)=str_c("prediction_", ncol(predictions))
  
  predictions %<>% bind_cols(pp)
  
  
  cm.train <- confusionMatrix(factor(mat$group[-exclude_samples], levels = unique(mat$group)), 
                              factor(predicted[-exclude_samples],levels = unique(mat$group)))
  #print("Train")
  #print(cm.train$overall)
  
  cm.test <- cm.test <- confusionMatrix(factor(mat$group[exclude_samples], levels = unique(mat$group)),
                                        factor(predicted[exclude_samples], levels = unique(mat$group)))
  
  
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

predictions %<>% mutate(sum_pred=rowSums(across(starts_with("prediction_")))) 
predictions %<>% mutate(pred_using_all=ifelse(sum_pred<=2, 0, 1))
print(table(predictions$pred_using_all, predictions$response))

print(str_c("mean accuracy:",  mean(accuracies.xgb$test, na.rm = T)))
print(str_c("mean sensetivity:",  mean(accuracies.xgb$test.sensetivity, na.rm = T)))
print(str_c("mean specificity:",  mean(accuracies.xgb$test.specificity, na.rm = T)))


{#modeling full set + simulated data
  xgb.model <- xgboost(data = as.matrix(mat[,sig_features]),
                       label = mat[,1],
                       max.depth = 2, #prms$max.depth[p],
                       eta = 0.2, #prms$eta[p],
                       nthread = 5, #prms$nthreads[p],
                       nrounds = 30, #prms$nrounds[p] , #Why 63? xgb.cv showed several times for it to be best fit
                       eval_metric = "error",
                       objective = "binary:logistic",
                       verbose = 0)
  
  
  numtosim_1=800
  numtosim_0=1200
  smat=bind_rows(data.frame(group=rep(1, numtosim_1)),
                 data.frame(group=rep(0, numtosim_0)));  
  
  for(s in sig_features)
  {
    df1=data.frame(runif(numtosim_1, 
                     min=min(range(mat[which(mat[,1]==1), s]))-2*sd(mat[which(mat[,1]==0), s]), 
                     max=max(range(mat[which(mat[,1]==1), s]))+2*sd(mat[which(mat[,1]==0), s])))
    colnames(df1)[1]=s
    
    
    df0=data.frame(runif(numtosim_0, 
                         min=min(range(mat[which(mat[,1]==0), s]))-2*sd(mat[which(mat[,1]==1), s]), 
                         max=max(range(mat[which(mat[,1]==0), s]))+2*sd(mat[which(mat[,1]==1), s])))
    colnames(df0)[1]=s
    
    smat %<>% bind_cols(bind_rows(df1, df0))
  }

  spredicted <- round(predict(xgb.model, as.matrix(smat[,sig_features])));
  
  confusionMatrix(factor(smat$group, levels = unique(smat$group)), 
                  factor(spredicted,levels = unique(smat$group)))
}


#predictions =bind_cols(dndrt_data_info, predictions)

# predictions %<>% separate(name, into=c("DF", "Line", "Well", "Day", "Well_filename", "Field_filename")) %<>% 
#   unite("filename", Well_filename, Field_filename, sep = "_")

# predictions %<>% group_by_at(vars(!starts_with("prediction_"))) %>% 
#   rowwise() %>%
#   mutate(mean_prediction=mean(c_across( starts_with("prediction_"))))
# 
# output_filename=str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
#                       "Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/rnaseq/primary_models/",
#                       "20210314_rnaseq_Bup_7d_modeling_predictions.csv")
# 
# write.csv(predictions, output_filename,
#           quote = TRUE,
#           row.names = FALSE)
