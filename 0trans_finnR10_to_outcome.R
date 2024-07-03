

rm(list = ls())
library(TwoSampleMR)
library(tidyverse)
library(data.table)
setwd("M:/metab1400_BC_MR/")
trait <- 'finngen_R10_C3_BREAST_EXALLC.gz'

finn <- fread(trait,data.table = F)
head(finn)

finndata0 <- format_data(finn,type = 'outcome',snp_col = "rsids", beta_col = "beta", se_col = "sebeta", 
                        eaf_col = "af_alt", effect_allele_col = "alt", other_allele_col = "ref", 
                        pval_col = "pval", chr_col = "#chrom", pos_col = "pos")
finndata <- finndata0 %>% separate_rows(SNP, sep = ",")

finn_info <- fread('Finn_R10_manifest.csv',data.table = F)

trait_row <- finn_info[grepl(trait,finn_info$path_https),]

finndata$ncase.outcome <- trait_row$num_cases
finndata$ncontrol.outcome <- trait_row$num_controls
finndata$samplesize.outcome <- trait_row$num_cases + trait_row$num_controls

outcome <- trait_row$phenotype
finndata$outcome <- outcome

save(finndata,file = '0finnR10_BC_as_outcome.rdata')


