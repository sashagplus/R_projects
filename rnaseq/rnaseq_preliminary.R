#setwd("~/GenetikaPlus_Projects_R/rnaseq")
library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library("rpart")
library(rpart.plot)



data_folder=stringr::str_c("/Users/sashakugel/gplus_dropbox/Genetika+ Dropbox/Genetika+SharedDrive/01_Protocol_Development/01_06_RNA-seq/Analysis/Analysis full/")

data_files=data.frame(type=c("exp", "grp", "exp", "grp", "exp", "grp"),
                      filename=c("Final_tables/2021_02_00/2021_02_22_7d_Bup_cpm_exclude_L38.csv",
                                 "2021_02_11_newbatch_groups_ex_L38.csv",
                                 "cit_exp_2021_02_22_7d_Citalopram_cpm.csv",
                                 "cit_grp_2021_02_14_7d_citalopram_groups.csv",
                                 "mirt_exp_2021_02_22_7d_Mirt_cpm.csv",
                                 "mirt_grp_2021_02_14_7d_mirtazapine_groups.csv"))

#read data
bup_grp=read.csv(file=stringr::str_c(data_folder,data_files$filename[2]))
bup_grp %<>% distinct

bup_exp=read.csv(file=stringr::str_c(data_folder,data_files$filename[1]))

sample_profiles=data.frame(lname=colnames(bup_exp)[-c(1:3)], exp_col_index=(1:ncol(bup_exp))[-c(1:3)] )
sample_profiles %<>% separate(col=.data$lname, into=c("sample", "rest"), sep="\\.", extra = "merge", remove=F)
sample_profiles %<>% mutate(treatment=case_when(grepl("bup", rest, ignore.case = TRUE) ~ "Bup",
                                                grepl("Treated", rest, ignore.case = TRUE) ~"Bup",
                                                TRUE ~ "Vehicle"))


sample_profiles %<>% dplyr::left_join(bup_grp, 
                                      by=c("sample", "treatment"))
sample_profiles %<>% mutate(nname=stringr::str_c(sample, "_",treatment))
sample_profiles %<>% group_by(sample, treatment ) 

bup_exp %<>% gather(key="lname", value="cnts", -ensemblID, -gene_ID, -data)

bup_exp %<>% left_join(sample_profiles %>% select(lname, sample, treatment, group, nname),
                       by=c("lname"))
bup_exp %<>% group_by(ensemblID, gene_ID, data, sample, treatment, group, nname) %<>% 
  summarise(med_cnt=median(cnts), .groups = "keep")

bup_exp %<>% ungroup() %<>% select(-nname) %<>% spread(key=treatment, value = med_cnt) %<>% as.data.frame()

bup_exp %<>%rowwise() %<>% mutate(fc=ifelse(Vehicle==0, 
                                            ifelse(Bup==0, 1, Bup), 
                                            Bup/Vehicle))

bup_exp %<>% mutate(logfc=log(fc,2))

cut_bup_exp = bup_exp%>% filter(logfc<=log(0.7,2) | logfc >=log(1.3,2)) 


mat=bup_exp %>% ungroup %>%select(ensemblID, sample, group, logfc) %>% spread (key=ensemblID, value=logfc)


#mat=cut_bup_exp %>% ungroup %>%select(ensemblID, sample, group, logfc) %>% spread (key=ensemblID, value=logfc)

mat %<>%
  mutate(across(starts_with("ENS"), ~ ifelse(is.infinite(.x), min(is.finite(.x))*1.5, .x)))

mat %<>% as.data.frame() %<>% tibble::column_to_rownames("sample")

ttests=mat %>%
  summarise_each(funs(t.test(.[group == "Non-Res"], .[group == "Responders"])$p.value), vars = -group)

ttests=as.data.frame(t(ttests))
ttests %<>% tibble::rownames_to_column("ensembleID")
colnames(ttests)[2]="pval"

# nas=colSums(is.na(mat))
# submat=mat[,names(which(nas==0))] %>% as.data.frame()
# for ( i in 3:ncol(submat)) {submat[is.infinite(submat[,i]),i]=min(submat[is.finite(submat[,i]), i])*1.5}

sig_features=ttests %>% filter(pval<=0.005) %>% pull(ensembleID)

genes=bup_exp %>% select(ensemblID, gene_ID, data) %>% distinct
genes %<>% mutate(labels=stringr::str_c(ensemblID, "  ", gene_ID, " (", data, ")"))

pdf("boxplots.pdf")
for(i in sig_features) {
  ll= genes %>% filter(ensemblID==i) %>% pull(labels)
  boxplot(mat[,i]~mat$group, ylab = ll)
  
}
dev.off()



# balance partition
mat.group=mat$group
resp_partition = split(sample(which(mat.group=="Responders")),1:8)
nonresp_partition = split(sample(which(mat.group=="Non-Res")),1:8)

mat %<>% mutate(group=ifelse(group=="Responders", 1, 0))
accs=c()
for(i in 1:8)
{
  exclude_samples=c(resp_partition[[i]], nonresp_partition[[i]])
  
  mylogit <- glm(group ~ ., data = mat[-exclude_samples, c("group",sig_features)], family = "binomial")
  
  predicted=round(predict(mylogit, newdata = mat[exclude_samples,], type = "response"))
  accs=c(accs, length(which(mat$group[exclude_samples]==predicted))/length(exclude_samples))
}


submat = mat[,c("group",sig_features)]