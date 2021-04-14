

rm(list=ls())

# full_history=read.csv("/Users/sashakugel/Documents/python_projects/PycharmProjects/remission-prediction-streamlit/data/stard_with_biological_data_full.csv",
#                   check.names = F)
# cnames_history=colnames(full_history)
# 
# full_history=read.csv("/Users/sashakugel/Documents/python_projects/PycharmProjects/remission-prediction-streamlit/data/stard_with_biological_data_full.csv",
#                   stringsAsFactors =F)
# 

full_bio=read.csv("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/api_development/streamlit/2021_04_01_data_for_streamlit_runs/stard_with_biological_data_revised_2021_SN_Spines_preprocessed.csv",
                  check.names  = F)
cnames_history=colnames(full_bio)
full_bio=read.csv("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/api_development/streamlit/2021_04_01_data_for_streamlit_runs/stard_with_biological_data_revised_2021_SN_Spines_preprocessed.csv",
                  stringsAsFactors = F)


DATA_DIR="/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/DataforPatientHistory/DepressionDataNIH/StarD/STAR_DClinicalData/ASCII Files raw data/"

crs = read.csv(str_c(DATA_DIR, 'Enrollment/CRS.csv')) %>% filter(ID %in% unique(full_bio$ID))
dm = read.csv(str_c(DATA_DIR , 'Enrollment/DM.csv')) %>% filter(ID %in% unique(full_bio$ID))
el = read.csv(str_c(DATA_DIR , 'Enrollment/EL.csv')) %>% filter(ID %in% unique(full_bio$ID))
hrsd = read.csv(str_c(DATA_DIR , 'Enrollment/HRSD.csv')) %>% filter(ID %in% unique(full_bio$ID))
mhx = read.csv(str_c(DATA_DIR , 'Enrollment/MHX.csv')) %>% filter(ID %in% unique(full_bio$ID))
pds = read.csv(str_c(DATA_DIR , 'Enrollment/PDS.csv')) %>% filter(ID %in% unique(full_bio$ID))
phx = read.csv(str_c(DATA_DIR , 'Enrollment/PHX.csv')) %>% filter(ID %in% unique(full_bio$ID))
cc = read.csv(str_c(DATA_DIR , 'Level1/Visits/CC.csv')) %>% filter(ID %in% unique(full_bio$ID))
qs = read.csv(str_c(DATA_DIR , 'Level1/Visits/QS.csv')) %>% filter(ID %in% unique(full_bio$ID) & DATE==0)

medication_columns = c('PSMD1', 'PSMD2', 'PSMD3', 'PSMD4', 'PSMD5', 'PSMD6', 'PSMD7', 'PSMD8')

mhx %<>% mutate(ZOL_EV=ifelse(PSMD1==107 | 
                                PSMD2==107 | 
                                PSMD3==107 | 
                                PSMD4==107 | 
                                PSMD5==107 | 
                                PSMD6==107 | 
                                PSMD7==107 | 
                                PSMD8==107, 1, 0))



#find the features from those. Gad built his model on those. 
# dfc=data.frame(feature=colnames(full_history), crs=NA,
#                                               dm=NA,
#                                               el=NA,
#                                               hrsd=NA,
#                                               mhx=NA,
#                                               pds=NA,
#                                               phx=NA,
#                                               cc=NA,
#                                               qs=NA)


df=crs %>% left_join(dm, by="ID")
df %<>% left_join(el, by="ID")
df %<>% left_join(hrsd, by="ID")
df %<>% left_join(mhx, by="ID")
df %<>% left_join(pds, by="ID")
df %<>% left_join(phx, by="ID")
df %<>% left_join(cc, by="ID")
df %<>% left_join(qs, by="ID")

            

# for(cn in 1:nrow(dfc))
# {
#   for(cc in 2:ncol(dfc))
#   {
#     if(dfc$feature[cn] %in% colnames(eval(parse(text=colnames(dfc)[cc])))) dfc[cn, cc]=TRUE
#   }
# }

clinical_history_model_features=c("QSTOT", "EMPL", "SAGIT", "SENGY", "BLACK", "HDTOT_R", "SMDSD", "SCHOOL", "HINSG", "HENGY", "HSANX", "TESHK", "HEMIN", "TRWIT", "TERMD", "WHITE", "FRLNE", "FRCAR", "PHACH", "HSUIC", "HMDSD", "ANAVD", "ZOL_EV", "epino", "SSOIN") 

df %<>% select(c(ID, clinical_history_model_features))

#history_data=data.frame()

full_bio%<>% select(-clinical_history_model_features)
full_bio %<>% left_join(df, by="ID")

full_bio %<>% select(ID, ID.1, Intermal.ID, Treatment, clinical_history_model_features, everything())
colnames(full_bio)=cnames_history

write.csv(full_bio,
          stringr::str_c("/Users/sashakugel/Documents/",
                         "python_projects/PycharmProjects/",
                         "remission-prediction-streamlit/",
                         "data/2021_04_01_data_for_streamlit_runs_full_new.csv"),
          quote = TRUE,
          row.names = FALSE)





