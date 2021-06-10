#cross reference 2020_04_07 streamlit output with what patients actually got
library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library(pROC)

rm(list=ls())

#functions:

calculate_score=function(df, 
                         new_score_name,
                         w_history =0, 
                         w_ge = 0, 
                         w_imaging = 0,
                         w_genetic =0,
                         w_side_effect = 0){
  
  df %<>% mutate(!!new_score_name := (w_history * .data$score_history + 
                                        w_ge * .data$score_gene_expression +
                                        w_imaging * .data$score_imaging   +
                                        w_genetic * .data$genetic_rank +
                                        w_side_effect * .data$side_effect_score) );
  return(df);
}

#twoColor == TRUE, then 2 color. Otherwise 3 color
calculate_category=function(df, 
                         score_column_name,
                         category_column_name,
                         twoColor = TRUE){
  
  if (twoColor)
  {
    df %<>% mutate(!!category_column_name := ifelse(!!as.name(score_column_name)>=0.5, 
                                                 "green", 
                                                 "red" ));
  }
  else
  {
    df %<>% mutate(!!category_column_name := ifelse(!!as.name(score_column_name) >= 1 , 
                                                  "green", 
                                                  ifelse(!!as.name(score_column_name) < 0, 
                                                         "red", 
                                                         "yellow")))
  }
  return(df);
}

timestamp()

#Processing: 
#read data
predicted=read.csv("/Users/sashakugel/Documents/python_projects/PycharmProjects/remission-prediction-streamlit/data/20210423_predicted_integrated_jupyter.csv",
                   stringsAsFactors = F) 
predicted %<>% select(-X)


actual=read.csv("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/api_development/streamlit/2021_04_01_data_for_streamlit_runs/20210406_treatmen_outcomes_ca.csv",
                stringsAsFactors = F) 
actual %<>% select(-X)


merged=actual %>% left_join(predicted, by=c("Patient.ID"="ID", "Main_Drug"="Treated_with"))


merged %<>% filter(Main_Drug==drugs)

#Scores to calculate: all possible combinations, with weight 1
# integrative score with weights (3,4,4,3,0.5)

score_types=c( "score_history", "score_gene_expression", "score_imaging", "genetic_rank", "side_effect_score")
score_abbreviation=c("h", "ge", "i", "g","se")
score_weights=seq(0,2, 0.25)
#score_weights=c(0.5, 3, 4)

feature_weight=merge(data.frame(score_types=score_types, score_abbreviation=score_abbreviation), score_weights) %>% rename( weight=y) %>% arrange(score_types)
feature_weight %<>% mutate(score_name=str_c(score_abbreviation, "_", weight))

feature_weight %<>% group_by(score_abbreviation)
df_list=feature_weight %>% group_split()

score_perms=expand.grid(df_list[[1]]$score_name, df_list[[2]]$score_name, df_list[[3]]$score_name, df_list[[4]]$score_name, df_list[[5]]$score_name) 
colnames(score_perms)=unlist(lapply(df_list, function(x) unique(x$score_types)))

score_perms %<>% mutate_all(~as.numeric(str_remove_all(as.character(.x), paste(str_c(score_abbreviation, "_"), collapse = "|"))))
  

score_perms %<>% mutate(score_num=1:n())
score_perms=score_perms[,c("score_num", score_types )] 
score_perms %<>% rowwise() %<>% mutate(score_name=str_c(c("score", str_c(score_abbreviation, 
                                                          c(score_history, score_gene_expression, score_imaging, genetic_rank, side_effect_score))), collapse="_"),
                                       category_name=str_c(c("category", str_c(score_abbreviation, 
                                                                         c(score_history, score_gene_expression, score_imaging, genetic_rank, side_effect_score))), collapse="_"),
                        .after=score_num)


print("moving to calculating scores")
score_perms$auc=0
for(sp in 1:nrow(score_perms))
{
  if((sp %% 1000)==0) {print(sp)}
  merged=calculate_score(merged, 
                              new_score_name = score_perms$score_name[sp], 
                              w_history = score_perms$score_history[sp],
                              w_ge = score_perms$score_gene_expression[sp],
                              w_imaging = score_perms$score_imaging[sp],
                              w_genetic = score_perms$genetic_rank[sp],
                              w_side_effect = score_perms$side_effect_score[sp])

  merged=calculate_category(merged, 
                            score_column_name = score_perms$score_name[sp], 
                            category_column_name = score_perms$category_name[sp])
    
  
  roc_obj <- roc(merged$Response_Main_Drug, merged[,score_perms$score_name[sp]], quiet = TRUE)
  score_perms$auc[sp]=auc(roc_obj)
}

# for(sp in 1:nrow(score_perms))
# {
# }

write.csv(score_perms,
          str_c("/Users/sashakugel/gplus_dropbox/",
                "Genetika+ Dropbox/Genetika+SharedDrive/",
                "01_Protocol_Development/01_10_Data_Science/",
                "api_development/streamlit/",
                "2021_04_01_data_for_streamlit_runs/",
                "20210423_predicted_auc_cutoff_weights0t2by0.25.csv"),
          row.names = F)

timestamp()
break()





# score_perms=data.frame()
# for(cc in 1:length(score_types))
# {
#   lfw=feature_weight_indices
#   colnames(lfw)=str_c(colnames(feature_weight_indices), "_" ,cc)
#   
#   if(cc>1)
#   {
#     lfw=crossing(score_perms, lfw)
#     .vars <- rlang::syms(colnames(lfw))
#     lfw %<>%  group_by(id = row_number()) %>% 
#       gather(col, value, -id) %>% 
#       filter(!any(duplicated(value))) %>% 
#       spread(col, value) %>%
#       ungroup() %>%
#       select(-id) %>% 
#       rowwise() %>%
#       mutate(grp = paste(sort(c(!!!.vars)), collapse = "_")) %>%
#       group_by(grp) %>%
#       slice(1) %>%
#       ungroup() %>%
#       select(-grp)  
#     
#     
#   }
#   score_perms %<>% bind_rows(lfw)
# }
# 
# score_perms %<>% mutate(score_num=1:n(), .before=fwi_1)
# 
# score_perms %<>% gather(key="key", val="val", -score_num) %>% 
#                 drop_na(val) %>% 
#                 left_join(feature_weight, by=c("val"="ind")) 
# 
# score_perms %<>% select(-val, -key)
# 
# dup_fwi=duplicated(score_perms %>% select(score_num, feature)) | duplicated(score_perms %>% select(score_num, feature), fromLast = TRUE)
# 
# score_perms=score_perms[-which(dup_fwi),]
# 
# score_perms %<>% spread(key=feature, value = weight, fill=0) 
# 
# score_perms %<>% select(-score_num) %>% distinct() 
# 
# 
#   cc_combn=combn(score_types,cc) %>% t
# colnames(cc_combn)=str_c("feature_", 1:ncol(cc_combn))
# 
# ww_combn=combn(score_weights,min(cc, length(score_weights))) %>% t
# colnames(ww_combn)=str_c("w_feature_", 1:ncol(ww_combn))
# 
#   
# score_perms %<>% mutate(weights_num=1:n(), .before=feature_1)
# score_perms %<>% mutate(weights_num=str_c("e",weights_num))
# score_perms %<>% gather(key="feature", val="non", -weights_num) %>% 
#   drop_na(non) %>% 
#   select(-feature) 
# 
# score_perms %<>% left_join(expand_grid(non=score_types,weights=seq(0,5, 0.1)), by=c("non")) 
# 
# score_perms %<>% mutate(weights_num=str_c(weights_num, ".w", weights)) 
# 
# score_perms %<>% spread(key=non, val=weights, fill=0) 


predicted %<>% mutate(score_weighted=3*score_history + 
                                      4*(score_gene_expression+score_imaging)   +
                                      3*genetic_rank +
                                      0.5*side_effect_score) %>%
              mutate(score_weighted_category=ifelse(score_weighted>=1, "green", ifelse(score_weighted<0, "red", "yellow")))

predicted %<>% mutate(score_genetics_only=0*score_history + 
                                          0*(score_gene_expression+score_imaging)   +
                                          1*genetic_rank +
                                          0*side_effect_score) %>%
                mutate(score_genetics_only_category=ifelse(score_genetics_only>=1, "green", ifelse(score_genetics_only<0, "red", "yellow")))


predicted %<>% mutate(score_genetics_history=1*score_history + 
                                          0*(score_gene_expression+score_imaging)   +
                                          1*genetic_rank +
                                          1*side_effect_score) %>%
                mutate(score_genetics_history_category=ifelse(score_genetics_history>=1, "green", ifelse(score_genetics_history<0, "red", "yellow")))

predicted %<>% mutate(score_history_only=1*score_history + 
                                              0*(score_gene_expression+score_imaging)   +
                                              0*genetic_rank +
                                              0*side_effect_score)  %>%
                mutate(score_history_only_category=ifelse(score_history_only>=1, "green", ifelse(score_history_only<0, "red", "yellow")))

predicted %<>% mutate(score_gene_expression_only=0*score_history + 
                                    1*score_gene_expression +
                                    0*score_imaging   +
                                    0*genetic_rank +
                                    0*side_effect_score)  %>%
              mutate(score_gene_expression_only_category=ifelse(score_gene_expression_only>=1, "green", ifelse(score_gene_expression_only<0, "red", "yellow")))

predicted %<>% mutate(score_imaging_only=0*score_history + 
                                                  0*score_gene_expression +
                                                  1*score_imaging   +
                                                  0*genetic_rank +
                                                  0*side_effect_score)  %>%
                            mutate(score_imaging_only_category=ifelse(score_imaging_only>=1, "green", ifelse(score_imaging_only<0, "red", "yellow")))




summed_3color=merged %>% group_by(Response_Main_Drug, score_category) %>% tally() %>% rename(category=score_category, gplus_weighted_count=n)

summed_3color %<>% left_join(merged %>% group_by(Response_Main_Drug, score_genetics_only_category) %>% tally() %>% rename(category=score_genetics_only_category, genetics_only_count=n),
                      by=c("Response_Main_Drug", "category"))

summed_3color %<>% left_join(merged %>% group_by(Response_Main_Drug, score_genetics_history_category) %>% tally() %>% rename(category=score_genetics_history_category, genetics_history_only_count=n),
                      by=c("Response_Main_Drug", "category"))

summed_3color %<>% left_join(merged %>% group_by(Response_Main_Drug, score_history_only_category) %>% tally() %>% rename(category=score_history_only_category, history_only_count=n),
                      by=c("Response_Main_Drug", "category"))

summed_3color %<>% left_join(merged %>% group_by(Response_Main_Drug, score_gene_expression_only_category) %>% tally() %>% rename(category=score_gene_expression_only_category, gene_expression_only_count=n),
                             by=c("Response_Main_Drug", "category"))

summed_3color %<>% left_join(merged %>% group_by(Response_Main_Drug, score_imaging_only_category) %>% tally() %>% rename(category=score_imaging_only_category, imaging_only_count=n),
                             by=c("Response_Main_Drug", "category"))

summed_3color %<>% gather(key="model", value="count", -Response_Main_Drug, -category) %>% 
                    spread(key=category, value=.data$count) %>% 
                    select(model, everything()) %>% 
                    arrange(model)

summed_3color %<>% rowwise() %>% 
                    mutate(total=sum(green, red, yellow, na.rm = T)) %>% 
                    mutate(true_positive=ifelse(Response_Main_Drug=="R", (green/total)*100, 0)) %>% 
                    mutate(false_positive=ifelse(Response_Main_Drug=="NR", (green/total)*100, 0)) %>%
                    mutate(true_negative=ifelse(Response_Main_Drug=="NR", (red/total)*100, 0)) %>%
                    mutate(false_negative=ifelse(Response_Main_Drug=="R", (red/total)*100, 0)) 





#2color
cutoff=0.5
predicted %<>% mutate(score_weighted_category=ifelse(score_weighted>=cutoff, "green", "red"))

predicted %<>% mutate(score_genetics_only_category=ifelse(score_genetics_only>=cutoff, "green", "red"))

predicted %<>% mutate(score_genetics_history_category=ifelse(score_genetics_history>=cutoff, "green", "red"))

predicted %<>% mutate(score_history_only_category=ifelse(score_history_only>=cutoff, "green", "red"))

predicted %<>% mutate(score_gene_expression_only_category=ifelse(score_gene_expression_only>=cutoff, "green", "red"))

predicted %<>% mutate(score_imaging_only_category=ifelse(score_imaging_only>=cutoff, "green", "red"))

merged=actual %>% left_join(predicted, by=c("Patient.ID"="ID", "Main_Drug"="Treated_with"))
merged %<>% filter(Main_Drug==drugs)

summed_2color=merged %>% group_by(Response_Main_Drug, score_weighted_category) %>% tally() %>% rename(category=score_weighted_category, gplus_weighted_count=n)

summed_2color %<>% left_join(merged %>% group_by(Response_Main_Drug, score_genetics_only_category) %>% tally() %>% rename(category=score_genetics_only_category, genetics_only_count=n),
                             by=c("Response_Main_Drug", "category"))

summed_2color %<>% left_join(merged %>% group_by(Response_Main_Drug, score_genetics_history_category) %>% tally() %>% rename(category=score_genetics_history_category, genetics_history_only_count=n),
                             by=c("Response_Main_Drug", "category"))

summed_2color %<>% left_join(merged %>% group_by(Response_Main_Drug, score_history_only_category) %>% tally() %>% rename(category=score_history_only_category, history_only_count=n),
                             by=c("Response_Main_Drug", "category"))

summed_2color %<>% left_join(merged %>% group_by(Response_Main_Drug, score_gene_expression_only_category) %>% tally() %>% rename(category=score_gene_expression_only_category, gene_expression_only_count=n),
                             by=c("Response_Main_Drug", "category"))

summed_2color %<>% left_join(merged %>% group_by(Response_Main_Drug, score_imaging_only_category) %>% tally() %>% rename(category=score_imaging_only_category, imaging_only_count=n),
                             by=c("Response_Main_Drug", "category"))


summed_2color %<>% gather(key="model", value="count", -Response_Main_Drug, -category) %>% 
                    spread(key=category, value=.data$count) %>% 
                    select(model, everything()) %>% 
                    arrange(model) 

summed_2color %<>% rowwise() %>% 
                  mutate(total=sum(green, red, na.rm = T)) %>% 
                  mutate(true_positive=ifelse(Response_Main_Drug=="R", (green/total)*100, 0)) %>% 
                  mutate(false_positive=ifelse(Response_Main_Drug=="NR", (green/total)*100, 0)) %>%
                  mutate(true_negative=ifelse(Response_Main_Drug=="NR", (red/total)*100, 0)) %>%
                  mutate(false_negative=ifelse(Response_Main_Drug=="R", (red/total)*100, 0)) 

print(merged %>% select(Response_Main_Drug, score_category) %>% table())
ss=sum(merged %>% select(Response_Main_Drug, score_category) %>% table())
print(merged %>% select(Response_Main_Drug, score_category) %>% table()/ss)

write.csv(summed_3color, 
          "/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/api_development/streamlit/2021_04_01_data_for_streamlit_runs/20210420_integration_comparison_3color.csv",
          row.names = F)

write.csv(summed_2color, 
          "/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/api_development/streamlit/2021_04_01_data_for_streamlit_runs/20210420_integration_comparison_2color.csv",
          row.names = F)

#merged %>% select(Treated_With=Main_Drug, Response_to_Treated=Response_Main_Drug, score_color) %>% table












