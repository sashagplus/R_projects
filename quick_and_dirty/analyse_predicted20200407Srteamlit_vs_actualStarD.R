#cross reference 2020_04_07 streamlit output with what patients actually got
library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)


rm(list=ls())

predicted=read.csv("/Users/sashakugel/Documents/python_projects/PycharmProjects/remission-prediction-streamlit/data/20210413_predicted_integrated_jupyter.csv",
                   stringsAsFactors = F)

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

predicted=calculate_score(predicted, new_score_name = "ggg", 
                          w_history = 3,
                          w_ge = 4,
                          w_imaging = 4,
                          w_genetic = 3,
                          w_side_effect = 0.5)

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


actual=read.csv("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/api_development/streamlit/2021_04_01_data_for_streamlit_runs/20210406_treatmen_outcomes_ca.csv",
                stringsAsFactors = F)


merged=actual %>% left_join(predicted, by=c("Patient.ID"="ID", "Main_Drug"="Treated_with"))


merged %<>% filter(Main_Drug==drugs)

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












