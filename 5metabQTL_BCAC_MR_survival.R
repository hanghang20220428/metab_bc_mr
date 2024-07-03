

rm(list = ls())
library(TwoSampleMR)
library(tidyverse)
library(data.table)
setwd("M:/metab1400_BC_MR/")
load("./0survival_data_BCAC/BCAC_survival_2021_meta.rdata")

#读取数据
alltrait <- fread('1400metab_trait.csv',data.table = F)

AllOR <- data.frame()

setwd('M:/metab1400_BC_MR/all_metab_1400_noLD_asExp_0.1_100kb')
expfile <- list.files(pattern = 'rdata')
for (i in expfile) {
  setwd('M:/metab1400_BC_MR/all_metab_1400_noLD_asExp_0.1_100kb')
  print(which(expfile==i))
  
  trait <-sapply(strsplit(i,split = '.rda'), '[',1)  
  trait0 <- gsub('~',':',trait)
  
  load(i)
  out <- BCAC_survival_2021_meta %>% dplyr::filter(SNP %in% exp$SNP)
  if (is.null(out)) {
    next              }
  else{
    repeat{dat <- try(TwoSampleMR::harmonise_data(exposure_dat = exp, outcome_dat = out,
                                                  action = 2))
    if(!('try-error' %in% class(dat))){
      break }}
    dat <- subset(dat, dat$mr_keep == TRUE)
    if (nrow(dat) == 0) {
      next            }
    else{
      res <- GagnonMR::primary_MR_analysis(dat = dat)
      res$exposure <- trait0 
      res$outcome <- "BCAC"
      OR <- generate_odds_ratios(res)
      if (res$pval>0.05 | res$pval=='NaN') {
        OR <- cbind(OR,het_Qval=NA,pleio_pval=NA,
                    direction=NA)
        AllOR <- rbind(AllOR,OR)
        next           }
      else{
        setwd("M:/metab1400_BC_MR/")
        if (OR$method == 'Inverse variance weighted') {
          het <- mr_heterogeneity(dat,method_list = 'mr_ivw')
          pleio <- mr_pleiotropy_test(dat)
          direct <- directionality_test(dat)
          OR <- cbind(OR,het_Qval=het$Q_pval,pleio_pval=pleio$pval,
                      direction=direct$correct_causal_direction)
          save(exp,out,dat,res,OR,file = paste0('./6metabQTL_BCAC_survival_MR/',trait,
                                                "_mrResult.rdata"))
          AllOR <- rbind(AllOR,OR)

          single <- mr_leaveoneout(dat)
          mr_leaveoneout_plot(single)
          ggsave(filename = paste0('./6metabQTL_BCAC_survival_MR/',trait,'_leaveoneout_plot.pdf'),width = 6,height = 5)
          

          mr_scatter_plot(res,dat)
          ggsave(filename = paste0('./6metabQTL_BCAC_survival_MR/',trait,'_mr_scatter_plot.pdf'),width = 6,height = 5)
          

          res_single <- mr_singlesnp(dat,all_method = c("mr_ivw")) 
          mr_forest_plot(res_single)
          ggsave(filename = paste0('./6metabQTL_BCAC_survival_MR/',trait,'_mr_forest_plot.pdf'),width = 5,height = 6)
          

          mr_funnel_plot(res_single)
          ggsave(filename = paste0('./6metabQTL_BCAC_survival_MR/',trait,'_mr_funnel_plot2.pdf'),width = 6,height = 4.5)
          
        }
        else{
          direct <- directionality_test(dat)
          OR <- cbind(OR,het_Qval=NA,pleio_pval=NA,
                      direction=direct$correct_causal_direction)
          save(exp,out,dat,res,OR,file = paste0('./6metabQTL_BCAC_survival_MR/',trait,
                                                "_mrResult.rdata"))
          AllOR <- rbind(AllOR,OR)
        }
      }
    }
  }
}

setwd("M:/metab1400_BC_MR/")


AllOR$FDR <- p.adjust(AllOR$pval,method = "BH")


AllOR$BFR <- p.adjust(AllOR$pval,method = "bonferroni")

save(AllOR,file = 'all_metabQTL_mrres_BCAC_survival.rdata')
write.csv(AllOR,file = 'all_metabQTL_mrres_BCAC_survival.csv')

AllOR_NOhete1 <- subset(AllOR,FDR<0.05 & het_Qval>0.05 & pleio_pval>0.05 & direction==TRUE) 
AllOR_NOhete2 <- subset(AllOR,FDR<0.05 & is.na(het_Qval) & direction==TRUE)
AllOR_NOhete_FDR <- rbind(AllOR_NOhete1,AllOR_NOhete2)
save(AllOR_NOhete_FDR,file = 'all_metabQTL_mrres_BCAC_survival_NOhete_FDR.rdata')
write.csv(AllOR_NOhete_FDR,file = 'all_metabQTL_mrres_BCAC_survival_NOhete_FDR.csv')

AllOR_NOhete11 <- subset(AllOR,BFR<0.05 & het_Qval>0.05 & pleio_pval>0.05 & direction==TRUE) 
AllOR_NOhete22 <- subset(AllOR,BFR<0.05 & is.na(het_Qval) & direction==TRUE)
AllOR_NOhete_BFR <- rbind(AllOR_NOhete11,AllOR_NOhete22)
save(AllOR_NOhete_BFR,file = 'all_metabQTL_mrres_BCAC_survival_NOhete_BFR.rdata')
write.csv(AllOR_NOhete_BFR,file = 'all_metabQTL_mrres_BCAC_survival_NOhete_BFR.csv')






