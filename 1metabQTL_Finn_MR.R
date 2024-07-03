

rm(list = ls())
library(TwoSampleMR)
library(tidyverse)
library(data.table)
setwd("M:/metab1400_BC_MR/")
load("M:/metab1400_BC_MR/0finnR10_BC_as_outcome.rdata")

setwd('D:/all_metab_1400')
all_metab <- list.files(pattern = 'gz')

alltrait <- fread('1400metab_trait.csv',data.table = F)

AllOR <- data.frame()
Wrongfiles <- data.frame()

for (i in all_metab) {

  print(which(all_metab==i))
  

  phecode<- sapply(strsplit(i,split = "_"),'[',1)
  
  trait0 <- alltrait$Trait[which(alltrait$GCST==phecode)]
  

  trait <- gsub('[*]','',trait0)
  trait <- gsub('[/]','_',trait)
  trait <- gsub('[:]','~',trait)
  
  setwd('D:/all_metab_1400')
  
  data0 <- try(fread(i,data.table = F, integer64 = "numeric"))
  head(data0)
  if('try-error' %in% class(data0)){
    Wrongfiles <- rbind(Wrongfiles,i)
    next
  }else{

    data1 <- data0 %>% dplyr::select("variant_id","effect_allele","other_allele","beta","standard_error",
                                     "effect_allele_frequency","p_value","chromosome","base_pair_location")
    colnames(data1) <- c("SNP","effect_allele","other_allele","beta","se",
                         "eaf",'p',"chr","pos")
    data1$samplesize <- alltrait$Sample_size[which(alltrait$GCST==phecode)]
    data2 <- na.omit(data1)
    

    exp_dat1 <- data2 %>% dplyr::select("SNP","p")
    

    colnames(exp_dat1)[1] <- 'rsid'
    colnames(exp_dat1)[2] <- 'pval'
    

    exp_dat2<- try(ieugwasr::ld_clump_local(exp_dat1,clump_p=5e-08,clump_r2=0.1,clump_kb=100,
                                            bfile="d:/1kg.v3/EUR",#欧洲人参考
                                            plink_bin="d:/plink_win64_20230116/plink.exe"))#也可以用plinkR路径
    if('try-error' %in% class(exp_dat2)){
      next
    }else{
      if (nrow(exp_dat2)>0) {
        

        sameSNP <- intersect(data2$SNP,exp_dat2$rsid)
        data <- data2[data2$SNP %in% sameSNP,] %>% 
          dplyr::arrange(p) %>% 
          dplyr::distinct(SNP,.keep_all = T)
        

        exp <- TwoSampleMR::format_data(data,type = "exposure",  snp_col = "SNP",
                                        beta_col = "beta",se_col = "se",
                                        eaf_col = "eaf",effect_allele_col = "effect_allele",
                                        other_allele_col = "other_allele",
                                        pval_col = "p",samplesize_col = "samplesize",
                                        chr_col = "chr",pos_col = "pos")

        exp$R2<-exp$beta.exposure*exp$beta.exposure*2*(exp$eaf.exposure)*(1-exp$eaf.exposure)
        exp$f<-(exp$samplesize.exposure-2)*exp$R2/(1-exp$R2)
        exp$f

        exp <- exp[exp$f>10,]
        
        if (is.null(exp)) {
          next            }
        else{
          out <- finndata %>% dplyr::filter(SNP %in% exp$SNP)
          if (is.null(out)) {
            next              }
          else{
            outSNP <- finndata$SNP[finndata$pval.outcome<1e-5] 
            
            exp <- exp %>% dplyr::filter(!(SNP %in% outSNP))
            if (is.null(exp)) {
              next            }
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
              res$outcome <- "FinnGen"
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
                  save(exp,out,dat,res,OR,file = paste0('./metabQTL_Finn_MR/',trait,
                                                        "_mrResult.rdata"))
                  AllOR <- rbind(AllOR,OR)

                  single <- mr_leaveoneout(dat)
                  mr_leaveoneout_plot(single)
                  ggsave(filename = paste0('./metabQTL_Finn_MR/',trait,'_leaveoneout_plot.pdf'),width = 6,height = 5)
                  
                  mr_scatter_plot(res,dat)
                  ggsave(filename = paste0('./metabQTL_Finn_MR/',trait,'_mr_scatter_plot.pdf'),width = 6,height = 5)
                  
                  res_single <- mr_singlesnp(dat,all_method = c("mr_ivw")) 
                  mr_forest_plot(res_single)
                  ggsave(filename = paste0('./metabQTL_Finn_MR/',trait,'_mr_forest_plot.pdf'),width = 5,height = 6)
                  
                  mr_funnel_plot(res_single)
                  ggsave(filename = paste0('./metabQTL_Finn_MR/',trait,'_mr_funnel_plot2.pdf'),width = 6,height = 4.5)
                  
                }
                else{
                  direct <- directionality_test(dat)
                  OR <- cbind(OR,het_Qval=NA,pleio_pval=NA,
                              direction=direct$correct_causal_direction)
                  save(exp,out,dat,res,OR,file = paste0('./metabQTL_Finn_MR/',trait,
                                                        "_mrResult.rdata"))
                  AllOR <- rbind(AllOR,OR)
                }
              }
            }
          }
        }
      }
    }
  }
}
}
setwd("M:/metab1400_BC_MR/")

AllOR$FDR <- p.adjust(AllOR$pval,method = "BH")

AllOR$BFR <- p.adjust(AllOR$pval,method = "bonferroni")

save(AllOR,file = '1all_metabQTL_mrres_Finn.rdata')
write.csv(AllOR,file = '1all_metabQTL_mrres_Finn.csv')

AllOR_NOhete1 <- subset(AllOR,pval<0.05 & het_Qval>0.05 & pleio_pval>0.05 & direction==TRUE) 
AllOR_NOhete2 <- subset(AllOR,pval<0.05 & is.na(het_Qval) & direction==TRUE)
AllOR_NOhete_pval <- rbind(AllOR_NOhete1,AllOR_NOhete2)
save(AllOR_NOhete_pval,file = 'all_metabQTL_mrres_Finn_NOhete_pval.rdata')
write.csv(AllOR_NOhete_pval,file = 'all_metabQTL_mrres_Finn_NOhete_pval.csv')

AllOR_NOhete1 <- subset(AllOR,FDR<0.05 & het_Qval>0.05 & pleio_pval>0.05 & direction==TRUE) 
AllOR_NOhete2 <- subset(AllOR,FDR<0.05 & is.na(het_Qval) & direction==TRUE)
AllOR_NOhete_FDR <- rbind(AllOR_NOhete1,AllOR_NOhete2)
save(AllOR_NOhete_FDR,file = 'all_metabQTL_mrres_Finn_NOhete_FDR.rdata')
write.csv(AllOR_NOhete_FDR,file = 'all_metabQTL_mrres_Finn_NOhete_FDR.csv')

AllOR_NOhete11 <- subset(AllOR,BFR<0.05 & het_Qval>0.05 & pleio_pval>0.05 & direction==TRUE) 
AllOR_NOhete22 <- subset(AllOR,BFR<0.05 & is.na(het_Qval) & direction==TRUE)
AllOR_NOhete_BFR <- rbind(AllOR_NOhete11,AllOR_NOhete22)
save(AllOR_NOhete_BFR,file = 'all_metabQTL_mrres_Finn_NOhete_BFR.rdata')
write.csv(AllOR_NOhete_BFR,file = 'all_metabQTL_mrres_Finn_NOhete_BFR.csv')

AllOR$label = ifelse(AllOR$BFR < 0.05,AllOR$exposure, NA)
AllOR$Effect <- ifelse(AllOR$FDR<0.05 & AllOR$b>0,"Positive",
                       ifelse(AllOR$FDR<0.05 & AllOR$b<0,"Negative","Insignificant"))

ggplot(AllOR,aes(x = b, y = -log10(FDR),color=Effect)) +
  geom_point(aes(size = abs(b)), alpha = 0.9) +
  scale_color_manual(values=c("#B3B3B3","#5E4FA2", "#f8c120"))+
  geom_vline(xintercept = 0, linetype = 2)+
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 6.5))+
  labs(x = "Beta (effect size)", 
       y = parse(text = "-log[10]*(FDR)"),
       title = 'All MR results between proteins and BC risk (FinnGen)')+
  ggrepel::geom_label_repel(aes(label = label),size = 3,
                            color="black",box.padding = unit(0.4, "lines"), 
                            segment.color = "black",
                            segment.size = 0.4)
ggsave(filename = '1volcano_plot_Finn.pdf',width = 6.5,height = 4.5) 




