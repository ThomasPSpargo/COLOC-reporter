library(data.table)
library(tidyverse)

path <- "/Users/tom/OneDrive - King's College London/PhD/PhD project/COLOC/FilesForManuscript"

###read in the lava results
lavaRes<- tibble(fread(file.path(path,"LAVAall_writeup.tsv")))

#Calculate smallest rg significance thresholds and the number of traits with a gwas significant snp
lavaRes <- lavaRes %>%
  mutate(min_pval_is_sig_5e.08_phen1=min_snp_pval_phen1<5e-08,
         min_pval_is_sig_5e.08_phen2=min_snp_pval_phen2<5e-08,
         sigthresh = case_when(p<0.05/605 ~ "Bonf_sig",
                               p.fdr<0.05 ~ "FDR_sig",
                               p<0.05 ~ "Nominal_sig",
                               TRUE ~ "nonsig"),
         ntraits_w_gwasSig_snp = case_when(min_pval_is_sig_5e.08_phen1 & min_pval_is_sig_5e.08_phen2 ~ 2,
                                           min_pval_is_sig_5e.08_phen1 | min_pval_is_sig_5e.08_phen2 ~ 1,
                                           TRUE ~ 0)
  )

#Crosstabulate
signif_crosstab<-xtabs(~ntraits_w_gwasSig_snp+sigthresh,lavaRes)

#Write out
write.csv(signif_crosstab,file.path(path,"signif_crosstab.csv"))


