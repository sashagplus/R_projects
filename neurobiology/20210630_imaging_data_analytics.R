library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library(Polychrome)
library(ggplot2)

rm(list=ls())

# load response
responses=gpio::dbGetQuery("select * from gpd.responsivness")

responses %<>% mutate(treat.abbrv=case_when(treatment_single_drug=="Bupropion" ~ "BUP",
                                            treatment_single_drug=="Citalopram (Celexa)" ~ "CIT",
                                            treatment_single_drug=="Mirtazapine" ~ "MIRT",
                                            treatment_single_drug=="Nortriptyline" ~ "NTP",
                                            treatment_single_drug=="Sertraline" ~ "STL",
                                            treatment_single_drug=="Tranylclypromine" ~ "TCL",
                                            treatment_single_drug=="Venlafaxine" ~ "VN"))

#################### PSD ####################
psd=gpio::dbGetQuery("select * from pgp.imaging_psd_processed_sevendays")

#Currently removing 92 mirt; we don't have the response for that.
psd %<>% filter(!(line==92 & treatment=="MIRT"))
#1. explore study design
plines=psd %>% select(df, line, treatment) %>% distinct %>% left_join(responses %>% select(id_gp_numeric,
                                                                        treat.abbrv,
                                                                        response_imp50p),
                                                    by=c("line"="id_gp_numeric", "treatment"="treat.abbrv")) %>%
              mutate(response_imp50p=ifelse(is.na(response_imp50p), "", response_imp50p))


#test if every sample has vehicle
#plines %>% select(-response_imp50p) %>% mutate(n=1) %>% spread(key=treatment, value=n) %>% arrange(line) %>% filter(is.na(VEH))

psd.grp=psd %>% select(-study_design_id, -file_id, -well, -field, -imagenumber, -objectnumber) %>%
  group_by(df, line, treatment) %>% 
  summarise_all(~median(.x, na.rm = T)) 

psd.grp %<>% left_join(plines, by=c("df", "line", "treatment")) %>%
              relocate(response_imp50p, .after=treatment)
  
pallete.psd.line=palette36.colors(n = length(unique(psd.grp$line)))

#to choose the colors play with the following: 
#pallete.psd.treatment <- createPalette(length(unique(psd.grp$treatment)), c("#010101", "#ff0000"), M=50); swatch(pallete.psd.treatment)

pallete.psd.treatment <- c("#89657C", "#F16965", "#4DE868", "#5551EB", "#E14DE2")

#pallete.psd.treatment=sky.colors(n = length(unique(psd.grp$treatment)))

psd.grp %<>% group_by(line) %>% 
  mutate(color.line=pallete.psd.line[cur_group_id()]) %>% 
  ungroup()

psd.grp %<>% group_by(treatment) %>% 
  mutate(color.treatment=pallete.psd.treatment[cur_group_id()]) %>% 
  ungroup()

psd.grp %<>% arrange(line, df, treatment)

psd.grp %<>% mutate(label=str_c(line, ":", df,":", case_when(treatment=="BUP" ~ "B",
                                                             treatment=="CIT" ~ "C",
                                                             treatment=="MIRT" ~ "M",
                                                             treatment=="NTP" ~ "N",
                                                             treatment=="VEH" ~"V")))




filter_from_string = function(xdf, ...) {
  filter(xdf, !!!parse_exprs(paste(..., sep = ";")))
}

pca_and_plotgg=function(mat,
                        treatment_filter=NULL,
                        PC1_filter=NULL,
                        PC2_filter=NULL,
                        label_column=NULL,
                        shape_column="treatment", 
                        linetype_column="treatment", 
                        group_column="treatment", 
                        color_column="treatment",
                        gtitle="")
{
  
  if (!is.null(treatment_filter)) 
  {
    noVEH=treatment_filter[-which(treatment_filter=="VEH")]
    indcs=mat %>% filter(treatment %in% noVEH) %>% select(df, line) %>% distinct
    mat=indcs %>% left_join(mat, by=c("df", "line")) %>% filter(treatment %in% treatment_filter)
  }
  
  identifier_features=c("line", "df","treatment", "color.line", "color.treatment", "label", "response_imp50p")
  
  mat.pca=bind_cols(mat %>% select(any_of(identifier_features)), 
                    as.data.frame(prcomp(mat %>% select(-any_of(identifier_features)))$x))
  
  if (!is.null(PC1_filter)) mat.pca=filter_from_string(xdf=mat.pca, PC1_filter) 
  if (!is.null(PC2_filter)) mat.pca=filter_from_string(xdf=mat.pca, PC2_filter) 
   
  #mat.pca %<>% mutate_at(vars(treatment, line), ~as.factor(.x))
  
  if(group_column=="treatment")
  {
    
    mat.pca %<>% arrange(treatment)
    color_manual=mat.pca %>% arrange(treatment) %>% pull(color.treatment) %>% unique
  }  else
  {#line
    
    geom_line_color=mat.pca %>% arrange(treatment) %>% pull(color.treatment) %>% unique
    mat.pca %<>% arrange(line)
    color_manual=mat.pca %>% arrange(line) %>% pull(color.line) %>% unique
  }

  
  
  g=ggplot(mat.pca, 
           aes(x = PC1, y = PC2, 
                      shape=factor(get(shape_column)), linetype=get(linetype_column), 
                      group=get(group_column), color=factor(get(color_column)))) + 
    geom_point(size=4) + 
    scale_color_manual(values = color_manual)+
    ggtitle(gtitle) + 
    #scale_fill_manual(values=geom_line_color)+
    labs(shape=str_to_title(shape_column), 
         linetype=str_to_title(linetype_column), 
         color=str_to_title(color_column))
  
  if(group_column=="treatment")
  {
    g=g+geom_line(aes(group=treatment)) 
  } else
  {#line
    g=g+geom_line(aes(group=treatment), color="grey") 
  }
  if(!is.null(label_column)) g=g+geom_text(aes(label=get(label_column)),hjust=0, vjust=0) 
  
  g=g+theme(plot.title=element_text(face="bold", hjust = 0.5, size=23),
                   axis.title.x=element_text(face="bold", size=18),
                   axis.title.y=element_text(face="bold", size=18))
  return(g)
  
}

#plots
output_directory="/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_10_Data_Science/Exploratory_Outputs/neurobiology_imaging/2021_07_04_overall_imaging_explore"
break;
print("creating plots into pdf:" )

# pdf(file=str_c(output_directory, "/", "2021_06_04_psd_pca_exploratories.pdf"),
#       width = 15, height=11)
#All treatments
pca_and_plotgg(mat = psd.grp,gtitle = "PSD PCA all lines colored by Treatment")

pca_and_plotgg(mat = psd.grp,gtitle = "PSD PCA all lines colored by Treatment 2K<PC1, -1K<PC2<2.5K",
                PC1_filter = c("PC1 >= 2000"),
              PC2_filter = c("PC2 >= -1000", "PC2 <= 2500")
               )
#BUP
pca_and_plotgg(mat = psd.grp,treatment_filter = c("BUP", "VEH"),
               group_column = "line", color_column = "line",
               gtitle = "PSD PCA:  BUP vs VEH lines colored by Line")


pca_and_plotgg(mat = psd.grp,treatment_filter = c("BUP", "VEH"),
               group_column = "line", color_column = "line",
               PC1_filter = c("PC1 >= 0"),
               gtitle = "PSD PCA:  BUP vs VEH lines colored by Line PC1>0")

#CIT
pca_and_plotgg(mat = psd.grp,treatment_filter = c("CIT", "VEH"),
               group_column = "line", color_column = "line",
               gtitle = "PSD PCA:  CIT vs VEH lines colored by Line")

pca_and_plotgg(mat = psd.grp,treatment_filter = c("CIT", "VEH"),
               group_column = "line", color_column = "line",
               PC1_filter = c("PC1 <= 2500", "PC1 >=-2500"),
               PC2_filter = c("PC2 >= 0"),
               gtitle = "PSD PCA:  CIT vs VEH lines colored by Line -2.5K<PC1<2.5K, PC2>0")

#NTP
pca_and_plotgg(mat = psd.grp,treatment_filter = c("NTP", "VEH"),
               group_column = "line", color_column = "line", label_column = "label",
               gtitle = "PSD PCA:  NTP vs VEH lines colored by Line")


#MIRT
pca_and_plotgg(mat = psd.grp,treatment_filter = c("MIRT", "VEH"),
               group_column = "line", color_column = "line",
               gtitle = "PSD PCA:  MIRT vs VEH lines colored by Line")

pca_and_plotgg(mat = psd.grp,treatment_filter = c("MIRT", "VEH"),
               group_column = "line", color_column = "line",
               PC1_filter = c("PC1 <= 1500"),
               PC2_filter = c("PC2 <= 0"),
               gtitle = "PSD PCA:  MIRT vs VEH lines colored by Line PC1<1.5K, PC2<0")

# dev.off()
