library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)

rm(list=ls())

study_design=gpio::dbGetQuery("select * from gpd.imaging_study_design")
study_design %<>% select(-file_id,
                          -spines_directory, -coloc_directory, 
                          -coloc_roi_xlsx, -coloc_psd_csv, -coloc_syn_csv)

roi=gpio::dbGetQuery("select * from gpd.imaging_coloc_roi")
#character columns, except well and field, are file information; no analyzable
roi %<>% select(-filename_dapi, -filename_map2, -filename_psd95, 
                -filename_synapsin, -pathname_dapi, -pathname_map2, 
                -pathname_psd95, -pathname_synapsin, -url_dapi, -url_map2, 
                -url_psd95, -url_synapsin     )


psd=gpio::dbGetQuery("select * from gpd.imaging_coloc_psd")
syn=gpio::dbGetQuery("select * from gpd.imaging_coloc_syn")
#normalize microns and pixels
#read file, select the proper features and normalize
norm_sql=stringr::str_c('select ipsam.database_colnames_psd,ipsam.database_colnames_syn  , ipsac."Correct_by_deviding_by_pix_to_um" , ipsac."action" 
                          from gpd.imaging_psd_syn_action_map ipsam left join
                          	   gpd.imaging_psd_syn_action_codes ipsac 
                          	    on ipsam."Correct_by_deviding_by_pix_to_um" =ipsac."Correct_by_deviding_by_pix_to_um" 
                          ')
norm_info=gpio::dbGetQuery(norm_sql)


#normalize 
psd %<>% left_join(study_design %>% select(id, coloc_um_to_pixel_ratio), 
                   by=c("study_design_id"="id")) %>% 
          mutate_at(vars(norm_info %>% filter(action=="divide by ration") %>% pull(database_colnames_psd)), 
                    ~.x*coloc_um_to_pixel_ratio) %>% 
          mutate_at(vars(norm_info %>% filter(action=="divide by ration square") %>% pull(database_colnames_psd)), 
                    ~.x*(coloc_um_to_pixel_ratio^2)) 

syn %<>% left_join(study_design %>% select(id, coloc_um_to_pixel_ratio), 
                   by=c("study_design_id"="id")) %>% 
          mutate_at(vars(norm_info %>% filter(action=="divide by ration") %>% pull(database_colnames_syn)), 
                    ~.x*coloc_um_to_pixel_ratio) %>% 
          mutate_at(vars(norm_info %>% filter(action=="divide by ration square") %>% pull(database_colnames_syn)), 
                    ~.x*(coloc_um_to_pixel_ratio^2)) 

spines=gpio::dbGetQuery("select * from gpd.imaging_spines")
#dnd_planar_angle, dnd_max_angle, spn_contact_area_squaremicrom - should be numeric
spines %<>% mutate_at(vars(dnd_planar_angle, dnd_max_angle, spn_contact_area_squaremicrom),
                      ~as.numeric(.x))
#dnd_base_coordinates discarded for now
#spn_assigned_type, spn_is_a_2d_spine - are single value at the moment; remove
spines %<>% select(-dnd_base_coordinates, -spn_assigned_type, -spn_is_a_2d_spine)

#calculate 3D distance between spn_base_coordinate and spn_head_position
spines %<>% mutate_at(vars(spn_base_coordinate, spn_head_position), ~stringr::str_remove(stringr::str_remove(.x, "\\("), "\\)")) %>% 
  separate(spn_base_coordinate, ", ", into = c("x.base", "y.base","z.base")) %>% 
  separate(spn_head_position, ", ", into = c("x.head", "y.head","z.head")) %>% 
  mutate_at(vars(x.base, y.base, z.base, x.head, y.head, z.head), ~as.numeric(.x)) %>% 
  mutate(spn_head_base_distance=sqrt((x.head-x.base)^2+(y.head-y.base)^2+(z.head-z.base)^2), .keep="unused") 

#transform dnd_terminal_type spn_spine_type spn_is_attached into numeric values for modeling
#dnd_terminal_type: Normal = 1, Branch = 2,  Generated = 3
#spn_spine_type: filopodia = 1  mushroom =2          stubby=3      thin=4  none=5 other-6
#spn_is_attached yes=1 no=0
spines %<>% mutate(dnd_terminal_type=case_when(dnd_terminal_type=="Normal" ~ "1",
                                              dnd_terminal_type=="Branch" ~ "2",
                                              dnd_terminal_type=="Generated" ~ "3"),
                  spn_spine_type=case_when(spn_spine_type=="filopodia" ~ "1",
                                           spn_spine_type=="mushroom" ~ "2",
                                           spn_spine_type=="stubby" ~ "3",
                                           spn_spine_type=="thin" ~ "4",
                                           spn_spine_type=="none" ~ "5",
                                           TRUE  ~ "6"),
                  spn_is_attached=case_when(spn_is_attached=="yes" ~ "1",
                                            spn_is_attached=="no" ~ "0")) %>%
          mutate_at(vars(dnd_terminal_type, spn_spine_type, spn_is_attached), ~as.numeric(.x)) 
      


#Remove all 3 days treatment
study_design %<>% filter(days!=3) %<>% select(-days, -coloc_um_to_pixel_ratio, -well_original)
roi %<>% filter(study_design_id %in% study_design$id)
psd %<>% filter(study_design_id %in% study_design$id)
syn %<>% filter(study_design_id %in% study_design$id)
spines %<>% filter(study_design_id %in% study_design$id)

#join with study design data and upload to database, for future analytics and modeling
#also load the coding used in the different tables
#processed data schema: pgp

#join with study design
join_relocate_writeDB=function(df, sd, new_table_name)
{
  df %<>% left_join(sd, by=c("study_design_id"="id", "well")) %>% 
    relocate(colnames(sd %>% select(-id)), .after = study_design_id)   
  
  sts=gpdb::db.management.addTableFromDataFrame(df = df, newTableName = new_table_name, schemaName = "pgp")
  print(sts)
  return(sts)
}

join_relocate_writeDB(df = roi, sd = study_design, new_table_name = "imaging_roi_processed_sevendays")
join_relocate_writeDB(df = psd, sd = study_design, new_table_name = "imaging_psd_processed_sevendays")
join_relocate_writeDB(df = syn, sd = study_design, new_table_name = "imaging_syn_processed_sevendays")
join_relocate_writeDB(df = spines, sd = study_design, new_table_name = "imaging_spines_processed_sevendays")


categorical_coding=data.frame(database_colnames=rep("dnd_terminal_type", 3),
                              chr_value=c("Normal", "Branch", "Generated"),
                              int_value=c(1:3)) %>%
                  bind_rows(data.frame(database_colnames=rep("spn_spine_type", 6),
                                       chr_value=c("filopodia", "mushroom", "stubby", "thin", "none", "<other>"),
                                       int_value=c(1:6))) %>%
                  bind_rows(data.frame(database_colnames=rep("spn_is_attached", 2),
                                       chr_value=c("yes", "no"),
                                       int_value=c(1,0))) %>%
                  mutate(table_name="imaging_spines_processed_sevendays")


gpdb::db.management.addTableFromDataFrame(df = categorical_coding, newTableName = "imaging_sevendays_categorical_encoding", schemaName = "pgp")





