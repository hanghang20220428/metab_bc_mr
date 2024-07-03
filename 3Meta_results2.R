

rm(list = ls())
library(meta)
library(tidyverse)
library(export)
library(data.table)
setwd('M:/metab1400_BC_MR')

out_plots <- function(filename,pic_width=5,pic_height=7){
  graph2png(file=filename,width=pic_width,height=pic_height)
  graph2pdf(file=filename,width=pic_width,height=pic_height)
}

load('1all_metabQTL_mrres_Finn_NOhete_FDR.rdata')
AllOR_Finn <- AllOR_NOhete_FDR

load('2.1all_metabQTL_mrres_BCAC_NOhete_FDR.rdata')
AllOR_BCAC <- AllOR_NOhete_FDR

load('2.5all_metabQTL_mrres_BCAC_iCOGS_NOhete_FDR.rdata')
AllOR_BCAC_iCOGS <- AllOR_NOhete_FDR

same_exp <- Reduce(intersect, list(AllOR_Finn$exposure, AllOR_BCAC$exposure,AllOR_BCAC_iCOGS$exposure))
save(same_exp,file = '3intersect_metab_finn_BCAC.rdata')

setwd('M:/metab1400_BC_MR/3Meta_results2/')
for (i in 1:length(same_exp)) {
  trait0 <- same_exp[i]
  trait <- gsub('[:]','~',trait0)
  
  sig_Finn1 <- AllOR_Finn[AllOR_Finn$exposure==trait0,]
  sig_Finn1$study <- 'FinnGen'
  sig_Finn1 <- sig_Finn1 %>% dplyr::select(exposure,study,b,se)
  
  sig_BCAC1 <- AllOR_BCAC[AllOR_BCAC$exposure==trait0,]
  sig_BCAC1$study <- 'BCAC_OncoArray'
  sig_BCAC1 <- sig_BCAC1 %>% dplyr::select(exposure,study,b,se)
  
  sig_BCAC1_iCOGS <- AllOR_BCAC_iCOGS[AllOR_BCAC_iCOGS$exposure==trait0,]
  sig_BCAC1_iCOGS$study <- 'BCAC_iCOGS'
  sig_BCAC1_iCOGS <- sig_BCAC1_iCOGS %>% dplyr::select(exposure,study,b,se)
  
  sig_merge <- rbind(sig_Finn1,sig_BCAC1,sig_BCAC1_iCOGS)
  
  sig_meta <- meta::metagen(TE = b, seTE = se, data = sig_merge, studlab = study, sm = "OR",
                            common = T,random = T)
  
  settings.meta("JAMA")
  meta::forest(sig_meta)
  out_plots(filename=paste0('Meta_BC_',trait) ,pic_width=8,pic_height=6)
  
  inf <- meta::metainf(sig_meta)
  forest(inf)
  out_plots(filename=paste0('Inf_BC_',trait),pic_width=8,pic_height=6)
}


