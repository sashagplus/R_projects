#cross reference 2020_04_07 streamlit output with what patients actually got

rm(list=ls())

predicted=read.csv("/Users/sashakugel/Documents/python_projects/PycharmProjects/remission-prediction-streamlit/data/20210413_predicted_integrated_jupyter.csv",
                   stringsAsFactors = F)


actual=read.csv("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/api_development/streamlit/2021_04_01_data_for_streamlit_runs/20210406_treatmen_outcomes_ca.csv",
                stringsAsFactors = F)


merged=actual %>% left_join(predicted, by=c("Patient.ID"="ID", "Main_Drug"="Treated_with"))


merged %<>% filter(Main_Drug==drugs)


print(merged %>% select(Response_Main_Drug, score_category) %>% table())
ss=sum(merged %>% select(Response_Main_Drug, score_category) %>% table())
print(merged %>% select(Response_Main_Drug, score_category) %>% table()/ss)

#merged %>% select(Treated_With=Main_Drug, Response_to_Treated=Response_Main_Drug, score_color) %>% table












