#map streamlit file input patients to side effects

rm(list=ls())

data_folder="/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/DataforPatientHistory/"

#1. Read side effects file
side_effects = read.csv(stringr::str_c(data_folder, "StahlsGuide/2020_12_08_Sorted_side_effects_table.csv"))



biodata=read.csv(stringr::str_c("/Users/sashakugel/Documents/",
                                "python_projects/PycharmProjects/",
                                "remission-prediction-streamlit/",
                                "data/2021_04_01_data_for_streamlit_runs_full_new.csv"), 
                check.names = F,
                 stringsAsFactors = F)
cnames=colnames(biodata)

biodata=read.csv(stringr::str_c("/Users/sashakugel/Documents/",
                                "python_projects/PycharmProjects/",
                                "remission-prediction-streamlit/",
                                "data/2021_04_01_data_for_streamlit_runs_full_new.csv"), 
                 stringsAsFactors = F)


#sedata = side effect data
sedata=read.csv(stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                              "Genetika+SharedDrive/01_Protocol_Development/",
                              "01_10_Data_Science/DataforPatientHistory/",
                              "StahlsGuide/20210404_stard_patients_side_effects_fields_data.csv"),
                stringsAsFactors = F)


sedata %<>% filter(id %in% unique(biodata$ID))

sedata %<>% group_by(id, condition) %>% filter(level==min(level)) 
sedata %<>% mutate(side_effect_present=ifelse(side_effect_present==TRUE, 1, 0))

sedata %<>% group_by(id, level, condition) %>% 
  summarise(side_effect_present=sum(side_effect_present, na.rm=T), .groups = "keep") %>% 
  mutate(side_effect_present=ifelse(side_effect_present>0, 1, 0)) 

side_effects %<>% select(-Question) %>% filter(row_number()<1)

sedata %<>% mutate(condition=case_when(condition == "Sexually Active"  ~ "SexuallyActive",
                                       condition == "Insomnia"  ~ "Insomniac",
                                       condition == "Obesity"  ~ "Obese",
                                        TRUE ~ condition))

sedata  %<>% tidyr::spread(key=condition,
                           value=side_effect_present)                   
sedata %<>% ungroup %<>% select(-level)                     

output= read.csv(text=str_c(colnames(side_effects), collapse = ","))                     
output %<>% bind_rows( sedata) %<>% select(ID=id, everything())


biodata %<>% left_join(output, by="ID")


#filter and sum  #remove filter later
biodata %<>% mutate(Treatment=ifelse(str_detect(.data$Treatment, "BUP"), "Bupropion", Treatment))
biodata %<>% mutate(Treatment=ifelse(str_detect(.data$Treatment, "CIT") |
                                       str_detect(.data$Treatment, "CIR"), "Citalopram", Treatment))
biodata %<>% mutate(Treatment=ifelse(str_detect(.data$Treatment, "MIRT"), "Mirtazapine", Treatment))
biodata %<>% mutate(Treatment=ifelse(str_detect(.data$Treatment, "NTP"), "Nortriptyline", Treatment))

features_to_sum=c("B.TPI1P2", "B.LIN37", "B.KLC3", "B.CXCL11", "T.TPI1P2", "T.LIN37", "T.KLC3", "T.CXCL11", "B.MinDistPSD_SYN_puncta", "T.MinDistPSD_SYN_puncta", "B.ThinSpineLength", "T.ThinSpineLength", "B.ThinSpine_per_Length", "T.ThinSpine_per_Length")

bio_biodata = biodata %>% select(ID, Treatment, all_of(features_to_sum))

bio_biodata %<>% mutate(TPI1P2=abs(T.TPI1P2-B.TPI1P2), 
                        LIN37=abs(T.LIN37- B.LIN37),
                        KLC3=abs(T.KLC3-B.KLC3),
                        CXCL11=abs(T.CXCL11-B.CXCL11),
                        MinDistPSD_SYN_puncta=abs(T.MinDistPSD_SYN_puncta-B.MinDistPSD_SYN_puncta),
                        ThinSpineLength=abs(T.ThinSpineLength-B.ThinSpineLength),
                        ThinSpine_per_Length=abs(T.ThinSpine_per_Length-B.ThinSpine_per_Length))

bio_biodata %<>% group_by(ID, Treatment) %>% summarise_all(median, na.rm=T)


biodata %<>% select(-ID.1, -Intermal.ID,  -ID.2, -all_of(features_to_sum))

biodata %<>% distinct()

biodata %<>% left_join(bio_biodata, by=c("ID", "Treatment"))



#add genetics data

genetic_rules=read.csv(stringr::str_c("/Users/sashakugel/Documents/python_projects/",
                                      "PycharmProjects/remission-prediction-streamlit/",
                                      "data/20210407_genetic_rules.csv"),
                       header = T, stringsAsFactors = F)



variables_to_transform= unique(genetic_rules$Variable)
  
genetics=read.table(stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                                "Genetika+SharedDrive/01_Protocol_Development/",
                                "01_10_Data_Science/DataforPatientHistory/",
                                "DepressionDataNIH/StarD/dp-s18-gen-tg-dataset6/",
                                "STARD_study18_genotypes/Hamilton_STARD_PK_Genotypes_082107.txt"), 
                  stringsAsFactors = F, sep="\t",
                  check.names = F, 
                  header=T)

gene_variables=colnames(genetics)[-1]

for(vtt in variables_to_transform)
{
  genetics %<>% mutate("{vtt}__1" := str_c(!!as.name(str_c(vtt, "_1")), !!as.name(str_c(vtt, "_2"))))
  genetics %<>% mutate("{vtt}__2" := str_c(!!as.name(str_c(vtt, "_2")), !!as.name(str_c(vtt, "_1")))) 
}



# genetics=read.table(stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
#                                    "Genetika+SharedDrive/01_Protocol_Development/",
#                                    "01_10_Data_Science/DataforPatientHistory/",
#                                    "DepressionDataNIH/StarD/dp-s18-gen-tg-dataset6/",
#                                    "STARD_study18_genotypes/Hamilton_STARD_PK_Genotypes_082107.txt"), 
#                     stringsAsFactors = F, sep="\t",
#                     header=T)
# 
ids_converter=read.table(stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/",
                                        "Genetika+SharedDrive/01_Protocol_Development/",
                                        "01_10_Data_Science/DataforPatientHistory/DepressionDataNIH/",
                                        "StarD/dp-s18-star-d-id-files/STAR_D_IDs_by_cellline.txt"), 
                         stringsAsFactors = F, sep="\t",
                         header=T)


genetics %<>% left_join(ids_converter, by=c("CellID"="cell.line"))  

genetics %<>% select(CellID, Patient.ID, starts_with(str_c(variables_to_transform, "__"))) 

genetics %<>% filter(Patient.ID %in% biodata$ID) 

genetics %<>% filter(CellID!="04C29695")

genetics %<>% select(-CellID)

biodata %<>% left_join(genetics, by=c("ID"="Patient.ID"))


cnames = cnames[which(! cnames %in% c( "Intermal.ID"))]
cnames=cnames[-which(cnames=="ID")[-1]]

biodata %<>% relocate(all_of(colnames(output)[-1]), .after = last_col()) 

colnames(biodata)[1:length(cnames)]=cnames
write.csv(biodata,
          stringr::str_c("/Users/sashakugel/Documents/",
                         "python_projects/PycharmProjects/",
                         "remission-prediction-streamlit/",
                         "data/2021_04_01_data_for_streamlit_runs_full_new_withSideEffects_wGenetics.csv"),
          quote = TRUE,
          row.names = FALSE)


print("Done.")








