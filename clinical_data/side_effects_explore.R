library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)

data_folder="/Users/sashakugel/data/"

#1. Read side effects file
side_effects = read.csv(stringr::str_c(data_folder, "StahlsGuide/2020_12_08_Sorted_side_effects_table.csv"))
side_effects=side_effects[-1,]
side_effects %<>% tidyr::gather(key="symptom", value="indication", -Question)
side_effects %<>% filter(indication!="")

stard_levels=c("Enrollment", "level1", "level2", "level2a", "level3", "level4")
conditions=read.csv(stringr::str_c(data_folder, "StahlsGuide/mapping side effects to starD features.csv"))
#2. Read patient reports files
# Read according to conditions data frame
conditions %<>% filter(Enrollment_file_loc!=""|
                         Enrollment_field!=""|
                         Level_file_loc!=""|
                         Level_field!="")
reports_folder="STAR_DClinicalData/STAR_DClinicalData/ASCII Files raw data/"
patient_reports=data.frame()
for(cnd in 1:nrow(conditions))
{
  for (stg in stard_levels)
  {
    if(stg=="Enrollment" & conditions[cnd, "Enrollment_file_loc"]!="")
    {
      read_data=read.csv(stringr:::str_c(data_folder, reports_folder, conditions[cnd, "Enrollment_file_loc"]))
      read_data %<>% select(any_of(c("ID","LEVEL")), !!conditions[cnd, "Enrollment_field"])
      read_data %<>% mutate(condition=conditions[cnd, "condition"],
                            file_level=stg,
                            field=conditions[cnd, "Enrollment_field"],
                            .before=ID)
      patient_reports %<>% bind_rows(read_data)
    }
    
    if(stg!="Enrollment" & conditions[cnd, "Level_file_loc"]!="")
    {
      read_data=read.csv(stringr:::str_c(data_folder,
                                         reports_folder,
                                         stg, "/",
                                         conditions[cnd, "Level_file_loc"]))
      read_data %<>% select(ID,LEVEL, !!conditions[cnd, "Level_field"])
      read_data %<>% mutate(condition=conditions[cnd, "condition"],
                            file_level=stg,
                            field=conditions[cnd, "Level_field"],
                            .before=ID)
      patient_reports %<>% bind_rows(read_data)
    }
  }
}

patient_reports %<>% tidyr::gather(key="field_",
                                   value="val",
                                   -condition,
                                   -file_level,
                                   -field,
                                   -ID,
                                   -LEVEL) %<>%
  filter(field==field_)

patient_reports %<>% mutate(condition_present = case_when(
  condition=="Obesity" & val>1 ~ TRUE,
  condition=="Fatigue" & val>0 ~ TRUE,
  condition=="Liver malfunction" & val>0 ~ TRUE,
  condition=="Sexually Active" & val>0 ~ TRUE,
  condition=="Anxiety" & val>0 ~ TRUE,
  condition=="Insomnia" & val>0 ~ TRUE,
  TRUE ~ FALSE
))

#3. Read patient treatment and outcome files.
#read CC files, and summarise properly (group by)

trtmnt_folder="STAR_DClinicalData/STAR_DClinicalData/ASCII Files for viewing/"
patients_trtmnt_files=c("CC_1.csv", "CC_2.csv", "CC_2A.csv", "CC_3.csv", "CC_4.csv")
patient_visits=data.frame()
for (stg in 2:6)
{
  
  cc=read.csv(stringr:::str_c(data_folder,
                              trtmnt_folder,
                              stard_levels[stg], "/Visits/CC_",
                              toupper(stringr::str_remove(stard_levels[stg], "level")),
                              ".csv"))
  cc$Study.med.code.3 = as.character(cc$Study.med.code.3)
  cc$Study.med.code.2 = as.character(cc$Study.med.code.2)
  cc$Prescribed.med.code.6 = as.character(cc$Prescribed.med.code.6)
  cc$Prescribed.med.daily.dose.5 = as.character(cc$Prescribed.med.daily.dose.5)
  cc$Prescribed.med.daily.dose.6 = as.character(cc$Prescribed.med.daily.dose.6)
  cc$Med.type = as.character(cc$Med.type)
  patient_visits %<>% bind_rows(cc)
  
}

#Filter entries with prescription
patients_prescribed=patient_visits %>% filter(Study.meds.prescribed=="Yes") %>%
  select(Patient.ID, Treatment.level, starts_with("Prescribed.med.code."))

patients_prescribed %<>% tidyr::gather(key="code_number", value="med", -Patient.ID, -Treatment.level)
patients_prescribed %<>% select(-code_number)
patients_prescribed %<>% distinct()
patients_prescribed %<>% mutate(med=ifelse(med=="Citalopram (Celexa)", "Citalopram", med))
patients_prescribed %<>% mutate(med=ifelse(med=="Tranylclypromine", "Tranylcypromine", med))
patients_prescribed %<>% mutate(Treatment.level=tolower(str_remove_all(Treatment.level, " ")))

patients_remission=patient_visits %>% filter(Patient.in.remission=="Yes") %>%
  select(Patient.ID, Treatment.level, Patient.in.remission)
patients_remission %<>% distinct()
patients_remission %<>% mutate(Treatment.level=tolower(str_remove_all(Treatment.level, " ")))


patients_prescribed %<>% filter(.data$med %in% side_effects$Question)

#cross reference
mat = patients_prescribed %>% left_join(patient_reports %>% select(condition, file_level, ID, condition_present),
                                        by=c("Patient.ID"="ID", "Treatment.level"="file_level"))

mat %<>% left_join(side_effects, by=c("med"="Question", "condition"="symptom"))

mat %<>% left_join(patients_remission, by=c("Patient.ID"="Patient.ID", "Treatment.level"="Treatment.level"))

mat %<>% mutate(Patient.in.remission=ifelse(is.na(Patient.in.remission), "No", Patient.in.remission))

#If Patient X went into Remission with drug D, how many times would we recommend drug D to such patients?













