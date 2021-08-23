

prevprovided=gpio::storage.read.csv(stringr::str_c("01_Protocol_Development/01_10_Data_Science/DATABASE/",
                                                  "base_files/2020_GP_Patient_Codes-StarD.csv")) %>% rename_all(~ str_c("gplus_", .x))


level_cellines=gpio::storage.read.csv(stringr::str_c("01_Protocol_Development/01_10_Data_Science/",
                                                     "DataforPatientHistory/DepressionDataNIH/",
                                                     "StarD/STAR_D_patient_history_clean/ASCII Files raw data/",
                                                     "Level1/cellline.csv")) %>% rename_all(~ str_c("level_", .x))



stard_ids=gpio::storage.read.csv(stringr::str_c("01_Protocol_Development/01_10_Data_Science/",
                                                     "DataforPatientHistory/DepressionDataNIH/",
                                                     "StarD/STAR_D_id_files/",
                                                     "dist302.study18_Auxiliary.csv")) %>% rename_all(~ str_c("stardids_", .x))



stard_genetics=gpio::storage.read.csv(stringr::str_c("01_Protocol_Development/01_10_Data_Science/DataforPatientHistory/",
                                                                   "DepressionDataNIH/StarD/STAR_D_genetics_set2_clean/",
                                                                   "STARD_Discovery_Genotypes/",
                                                                   "SD Name to CellID Converter For STARD.txt"),
                                                    sep="\t") %>% rename_all(~ str_c("stardGen_", .x))




















