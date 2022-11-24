
# Objective: Explore BP PRS variants with potential confounders (e.g. BMI, obsesity, smoking)
# From PhenoScanner: http://www.phenoscanner.medschl.cam.ac.uk/
# Author: Xiaonan Liu
# Date: 18Nov2021

library(data.table)
library(yaml)
library(dplyr)
library(tidyverse)
library(phenoscanner)

config = yaml.load_file("config.yml")

# Read files
BP <- read.csv(file.path(config$Evangelou2018$inputs,"Betas.csv"))

##
#res8<-phenoscanner(snpquery = "rs10418305",pvalue=5e-8,proxies = "EUR",r2=0.8,build = 37)

# 1. Extract list of SNPs and transform in the right format (i.e. add "chr" to ID not starting with "rs")
BP_SNPs=BP%>%
  mutate(SNP_query=ifelse(str_detect(rsid,"rs"),rsid,paste0("chr",rsid)))%>%pull(SNP_query)

# Note: The reason why I put below in our own function is to add the ability to filter for traits of interest
# between batches to save memory.

# phenoscanner() actually run 10 SNPs for one query so even though one specify batch size as 100, it still
# splits them into size of 10 blocks.

Run_PhenoScanner<-function(full_list=BP_SNPs,batch_size=100,
                           trait_ID=c("EFO_0004340","EFO_0001073","EFO_0006527"),
                           pvalue){
  # Split total BP SNPs into batch size of 100
  
  batches=split(1:length(full_list), ceiling(seq_along(1:length(full_list))/batch_size))
  
  pheno_subsets = list()
  
  for (i in 1:length(batches)){
    
    # Run phenoscanner
    res<-phenoscanner(snpquery = full_list[batches[[i]]],pvalue=pvalue,proxies = "EUR",r2=0.8,build = 37)
    
    # For assoc result for each SNP, we only keep trait of BMI, obesity and smoking 
    pheno_subsets[[i]]<-res$result%>%
      filter(efo %in% trait_ID & ancestry=="European")
    
    Sys.sleep(5)
    
  }
  
  pheno_subsets_df = do.call(rbind, pheno_subsets)
  
  return(pheno_subsets_df)
  
}



pheno_subsets_df<-Run_PhenoScanner(pvalue=5e-8)

write.csv(pheno_subsets_df,file=file.path(config$Evangelou2018$outputs,"BmiObsSmok_PhenoScan_5e-8.csv"))

bmi_obs_smok<-pheno_subsets_df%>%
  ## Keep only important columns
  select(snp,ref_hg19_coordinates,ref_a1,ref_a2,efo)%>%
  # select unique rsid
  distinct(.keep_all = TRUE)%>%
  mutate(hg19_coordinates=str_sub(ref_hg19_coordinates,4,-1))

# For checking
#bmi_obs_smok%>%filter(hg19_coordinates%in%inner$V1)%>%View()


# Export chr:pos so we can use --exclude in PLINK 

write.table(bmi_obs_smok[c("hg19_coordinates")],file = file.path(config$Evangelou2018$outputs,"BmiObsSmok_5e-8.txt"),
            col.names = FALSE,row.names = FALSE,quote = FALSE)


## Run for different threshold
pheno_subsets_df<-Run_PhenoScanner(pvalue=5e-6)

write.csv(pheno_subsets_df,file=file.path(config$Evangelou2018$outputs,"BmiObsSmok_PhenoScan_5e-6.csv"))

bmi_obs_smok<-pheno_subsets_df%>%
  select(snp,ref_hg19_coordinates,ref_a1,ref_a2,efo)%>%
  # select unique rsid
  distinct(.keep_all = TRUE)%>%
  mutate(hg19_coordinates=str_sub(ref_hg19_coordinates,4,-1))

write.table(bmi_obs_smok[c("hg19_coordinates")],file = file.path(config$Evangelou2018$outputs,"BmiObsSmok_5e-6.txt"),
            col.names = FALSE,row.names = FALSE,quote = FALSE)

## Check if BMI/Obesity assoc SNPs are on FTO gene

FTO_check<-phenoscanner(snpquery = bmi_obs_smok$ref_hg19_coordinates,pvalue=5e-8,proxies = "EUR",build = 37)

unique(FTO_check$snps$hgnc)

###########################################################
# Below is my try on searching variants assoc with BMI, obesity and smoking first then compare with BP variants
# Note: The downside of this approach is that it doesn't seem to take LD into account.

#bmi_pheno<-read.table(file.path(config$Evangelou2018$inputs,'Body_mass_index_PhenoScanner_GWAS.tsv'),sep='\t',header=TRUE)

#obs_pheno<-read.table(file.path(config$Evangelou2018$inputs,'Obesity_PhenoScanner_GWAS.tsv'),sep='\t',header=TRUE)

#smok_pheno<-read.table(file.path(config$Evangelou2018$inputs,'Smoking_status_PhenoScanner_GWAS.tsv'),sep='\t',header=TRUE)

## 1. Filter on EFO (https://www.ebi.ac.uk/efo/)

#bmi_pheno_subset<-bmi_pheno%>%
 # filter(efo=="EFO_0004340" & ancestry=="European")%>%
 # select(rsid,hg19_coordinates,a1,a2)%>%
  #mutate(BMI=1)%>%
  # remove duplicates SNPs
  #distinct(.keep_all = TRUE)


#obs_pheno_subset<-obs_pheno%>%
#  filter(efo=="EFO_0001073" & ancestry=="European")%>%
#  select(rsid,hg19_coordinates,a1,a2)%>%
#  mutate(obs=1)%>%
  # remove duplicates SNPs
 # distinct(.keep_all = TRUE)

#smok_pheno_subset<-smok_pheno%>%
#  filter(efo=="EFO_0006527" & ancestry=="European")%>%
#  select(rsid,hg19_coordinates,a1,a2)%>%
#  mutate(smok=1)%>%
  # remove duplicates SNPs
#  distinct(.keep_all = TRUE)


# After checking no multi-allelic SNPs in any of this pheno_subset above, and all the reference allele matched across 3,
# we can full join them to obtain the full list of SNPs associated with BMI or obesity or smoking.

#bmi_obs_smok<-full_join(bmi_pheno_subset,obs_pheno_subset,by=c("rsid","hg19_coordinates","a1","a2"))%>%
#  full_join(.,smok_pheno_subset,by=c("rsid","hg19_coordinates","a1","a2"))%>%
#  mutate(hg19_coordinates=str_sub(hg19_coordinates,4,-1))
                              

# Inner join with Evangelou2018 summary stats

#inner1<-inner_join(bmi_obs_smok,BP,by=c("hg19_coordinates"="chr_pos","a1"="effect_allele","a2"="reference_allele"))

#inner2<-inner_join(bmi_obs_smok,BP,by=c("hg19_coordinates"="chr_pos","a1"="reference_allele","a2"="effect_allele"))

#inner<-union(inner1,inner2)


# Export chr:pos so we can use --exclude in PLINK 

#write.table(inner[c("hg19_coordinates")],file = file.path(config$Evangelou2018$outputs,"bmi_obs_smok.txt"),
 #           col.names = FALSE,row.names = FALSE,quote = FALSE)

