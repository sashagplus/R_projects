library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)

rm(list=ls())
data_folder="/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/DataforPatientHistory/"

#1. Read side effects file
side_effects = read.csv(stringr::str_c(data_folder, "StahlsGuide/2020_12_08_Sorted_side_effects_table.csv"))
side_effects=side_effects[-1,]
side_effects %<>% tidyr::gather(key="symptom", value="indication", -Question)
side_effects %<>% filter(indication!="")

stard_levels=c("Enrollment", "level1", "level2", "level2a", "level3", "level4")
locations=read.csv(stringr::str_c(data_folder, "StahlsGuide/2021_04_04_mapping_side_effects_to_starD_features.csv"),
                   stringsAsFactors = F)


#2. Read patient reports files
# Read according to conditions data frame
# locations %<>% filter(level.or.enrollement!=""|
#                          information.file!=""|
#                          folder!=""|
#                          file!="")

locations %<>% 
  mutate(across(where(is.character), str_trim))

locations %<>% mutate(form.field.code=tolower(form.field.code))

locations %<>% filter(Use=="Yes") 


reports_folder="/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/DataforPatientHistory/DepressionDataNIH/StarD/STAR_DClinicalData/ASCII Files raw data/"
patient_reports=data.frame()

for(loc in 1:nrow(locations))
{
  for (stg in stard_levels)
  {
    if(stg=="Enrollment" & locations$level.or.enrollement[loc]=="Enrollment")
    {
      read_data=read.csv(stringr:::str_c(reports_folder, locations$level.or.enrollement[loc], "/", locations$file[loc]))
      read_data %<>% rename_all(tolower)
      read_data %<>% select(any_of(c("id","level")), !!locations$form.field.code[loc])
      read_data %<>% mutate(condition=locations$condition[loc],
                            file_level=stg,
                            field=locations$form.field.code[loc])
      read_data %<>% mutate(side_effect_present=ifelse(!!sym(locations$form.field.code[loc]) >= locations$cutoff[loc], TRUE, FALSE))
      read_data %<>%  select(-!!sym(locations$form.field.code[loc]))
      patient_reports %<>% bind_rows(read_data)
    }
    
    if(stg!="Enrollment" & locations$level.or.enrollement[loc]=="Level")
    {
      read_data=read.csv(stringr::str_c(reports_folder,
                                         stg, "/",
                                        locations$folder[loc], "/", 
                                        locations$file[loc]))
      read_data %<>% rename_all(tolower)
      read_data %<>% select(id,level, !!locations$form.field.code[loc])
      read_data %<>% mutate(condition=locations$condition[loc],
                            file_level=stg,
                            field=locations$form.field.code[loc])
      read_data %<>% mutate(side_effect_present=ifelse(!!sym(locations$form.field.code[loc]) >= locations$cutoff[loc], TRUE, FALSE))
      read_data %<>%  select(-!!sym(locations$form.field.code[loc]))
      patient_reports %<>% bind_rows(read_data)
    }
  }
}


write.csv(patient_reports,
          stringr::str_c(data_folder,
                         "StahlsGuide/20210404_stard_patients_side_effects_fields_data.csv"),
          quote = TRUE,
          row.names = FALSE)


print("Done.")










