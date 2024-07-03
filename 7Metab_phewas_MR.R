


rm(list = ls())
library(TwoSampleMR)
library(tidyverse)
library(data.table)
 
setwd("H:/FinnR9_all_outcome_rdata")
outfile <- list.files(pattern = 'rdata')
metab <- c('3,5-dichloro-2,6-dihydroxybenzoic acid levels','Carnitine C14 levels','Epiandrosterone sulfate levels',
           'Glyco-beta-muricholate levels','N4-acetylcytidine levels')

AllOR <- data.frame()
for (j in metab) {
  setwd("M:/metab1400_BC_MR/all_metab_1400_noLD_asExp_0.1_100kb")
  load(paste0(j,'.rdata'))
for (i in outfile) {
  setwd("H:/FinnR9_all_outcome_rdata")
  load(i)
  trait <-sapply(strsplit(i,split = '_finndata'), '[',1)  
  
  out <- finndata %>% dplyr::filter(SNP %in% exp$SNP)
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
      res$exposure <- j
      res$outcome <- trait
      OR <- generate_odds_ratios(res)
      
      setwd("M:/metab1400_BC_MR/7Metab_phewas_MR")
      if (OR$method == 'Inverse variance weighted') {
        het <- mr_heterogeneity(dat,method_list = 'mr_ivw')
        pleio <- mr_pleiotropy_test(dat)
        direct <- directionality_test(dat)
        OR <- cbind(OR,het_Qval=het$Q_pval,pleio_pval=pleio$pval,
                    direction=direct$correct_causal_direction)
        save(exp,out,dat,res,OR,file = paste0(j,'——',trait,
                                              "_mrResult.rdata"))
        AllOR <- rbind(AllOR,OR)

        single <- mr_leaveoneout(dat)
        mr_leaveoneout_plot(single)
        ggsave(filename = paste0(j,'——',trait,'_leaveoneout_plot.pdf'),width = 6,height = 5)
        

        mr_scatter_plot(res,dat)
        ggsave(filename = paste0(j,'——',trait,'_mr_scatter_plot.pdf'),width = 6,height = 5)
        

        res_single <- mr_singlesnp(dat,all_method = c("mr_ivw")) 
        mr_forest_plot(res_single)
        ggsave(filename = paste0(j,'——',trait,'_mr_forest_plot.pdf'),width = 5,height = 6)
        

        mr_funnel_plot(res_single)
        ggsave(filename = paste0(j,'——',trait,'_mr_funnel_plot2.pdf'),width = 6,height = 4.5)
        
      }
      else{
        direct <- directionality_test(dat)
        OR <- cbind(OR,het_Qval=NA,pleio_pval=NA,
                    direction=direct$correct_causal_direction)
        save(exp,out,dat,res,OR,file = paste0(j,'——',trait,
                                              "_mrResult.rdata"))
        AllOR <- rbind(AllOR,OR)
      }
    }
  }
}
}
table(AllOR$exposure)
setwd("M:/metab1400_BC_MR/")


AllOR$FDR <- p.adjust(AllOR$pval,method = "BH")
AllOR$BFR <- p.adjust(AllOR$pval,method = "bonferroni")

save(AllOR,file = '7Metab_phewas_MR.rdata')
write.csv(AllOR,file = '7Metab_phewas_MR.csv')

AllOR_NOhete1 <- subset(AllOR,het_Qval>0.05 & pleio_pval>0.05 & direction==TRUE) 
AllOR_NOhete2 <- subset(AllOR,is.na(het_Qval) & direction==TRUE)
AllOR_NOhete <- rbind(AllOR_NOhete1,AllOR_NOhete2)
save(AllOR_NOhete,file = '7Metab_phewas_MR_NOhete.rdata')

AllOR_NOhete1 <- subset(AllOR,FDR<0.05 & het_Qval>0.05 & pleio_pval>0.05 & direction==TRUE) 
AllOR_NOhete2 <- subset(AllOR,FDR<0.05 & is.na(het_Qval) & direction==TRUE)
AllOR_NOhete <- rbind(AllOR_NOhete1,AllOR_NOhete2)
save(AllOR_NOhete,file = '7Metab_phewas_MR_FDR_NOhete.rdata')
fwrite(AllOR_NOhete,file = '7Metab_phewas_MR_FDR_NOhete.csv')

AllOR_NOhete1 <- subset(AllOR,BFR<0.05 & het_Qval>0.05 & pleio_pval>0.05 & direction==TRUE) 
AllOR_NOhete2 <- subset(AllOR,BFR<0.05 & is.na(het_Qval) & direction==TRUE)
AllOR_NOhete <- rbind(AllOR_NOhete1,AllOR_NOhete2)
save(AllOR_NOhete,file = '7Metab_phewas_MR_BFR_NOhete.rdata')
write.csv(AllOR_NOhete,file = '7Metab_phewas_MR_BFR_NOhete.csv')





