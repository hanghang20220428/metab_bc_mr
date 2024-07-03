

rm(list = ls())
library(tidyverse)
library(cowplot)
setwd('M:/metab1400_BC_MR')

load('1all_metabQTL_mrres_Finn_NOhete_FDR.rdata')
AllOR_Finn <- AllOR_NOhete_FDR

load('2.1all_metabQTL_mrres_BCAC_NOhete_FDR.rdata')
AllOR_BCAC <- AllOR_NOhete_FDR

load('2.5all_metabQTL_mrres_BCAC_iCOGS_NOhete_FDR.rdata')
AllOR_BCAC_iCOGS <- AllOR_NOhete_FDR

same_exp <- Reduce(intersect, list(AllOR_Finn$exposure, AllOR_BCAC$exposure,AllOR_BCAC_iCOGS$exposure))

sig_Finn <- AllOR_Finn[AllOR_Finn$exposure %in% same_exp,]
sig_Finn$study <- 'FinnGen'
colnames(sig_Finn)
sig_Finn <- sig_Finn %>% dplyr::select(outcome,exposure,study,or,or_lci95,or_uci95,FDR)

sig_BCAC <- AllOR_BCAC[AllOR_BCAC$exposure %in% same_exp,]
sig_BCAC$study <- 'BCAC_OncoArray'
sig_BCAC <- sig_BCAC %>% dplyr::select(outcome,exposure,study,or,or_lci95,or_uci95,FDR)

sig_BCAC_iCOGS <- AllOR_BCAC_iCOGS[AllOR_BCAC_iCOGS$exposure %in% same_exp,]
sig_BCAC_iCOGS$study <- 'BCAC_iCOGS'
sig_BCAC_iCOGS <- sig_BCAC_iCOGS %>% dplyr::select(outcome,exposure,study,or,or_lci95,or_uci95,FDR)

sig_metab <- rbind(sig_Finn,sig_BCAC,sig_BCAC_iCOGS)
colnames(sig_metab)

ggplot(sig_metab, aes(y=exposure, x=or, colour=study)) +
  geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95), height=.2) +
  geom_point(size=2)+
  scale_color_manual(values=c("#4197d8","#9E0142", "#f8c120"))+
  geom_vline(xintercept=1, linetype='longdash') +
  theme_minimal_hgrid(10, rel_small = 1) +
  facet_wrap(~study, ncol=1)+
  labs(y = "", x = "OR (95%CI)",
       title="MR results between significant metabolites and BC risk")+
  theme(legend.position = "none",
        strip.text = element_text(face = 'bold'),
        axis.title.x=element_text(size=12,face = 'bold'),
        axis.title.y=element_text(size=12,face = 'bold'))

ggsave(filename = '4.2forest_Allmetab_MR.pdf',width=8,height=5)










