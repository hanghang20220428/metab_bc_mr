

rm(list = ls())
library(tidyverse)
setwd('M:/metab1400_BC_MR')
load("./all_metabQTL_mrres_BCAC_survival_NOhete_BFR.rdata")

AllOR_NOhete_BFR <- AllOR_NOhete_BFR %>% dplyr::select(exposure,b,BFR)


ggplot(AllOR_NOhete_BFR,aes(b,forcats::fct_reorder(exposure, b))) + 
  geom_segment(aes(xend=0, yend = exposure)) +
  geom_point(aes(color=BFR, size = abs(b))) +
  scale_color_gradient(low = "#1A9993FF",high = "#9E0142") +
  scale_size_continuous(range=c(5, 6)) +
  theme_bw() + 
  xlab("Beta") +
  ylab(NULL) + 
  ggtitle("Significant MR results between metabolites and breast cancer survival")
ggsave(filename = 'Dot_Allmetab_survival_MR.pdf',width=9,height=4)










