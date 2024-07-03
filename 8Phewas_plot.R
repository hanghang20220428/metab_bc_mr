

rm(list = ls())
library(tidyverse)
library(data.table)
library(tibble)
library(circlize)
library(ComplexHeatmap)
setwd("M:/metab1400_BC_MR/")

load("M:/metab1400_BC_MR/7Metab_phewas_MR_BFR_NOhete.rdata")
finninfo <- fread('Finn_R9_data.csv',data.table = F)
head(AllOR_NOhete)
head(finninfo)
data0 <- merge(AllOR_NOhete,finninfo,by.x = 'outcome',by.y = 'phenocode')
head(data0)

data <- data0 %>% dplyr::select(exposure,name,b,category)
data$category <- sapply(strsplit(data$category,split = '[(]'),'[',1 )
data$exposure <- as.factor(data$exposure)
data$category <- as.factor(data$category)
levels(data$exposure) 
levels(data$category)

data <- pivot_wider(data,names_from=exposure,values_from = b) %>% column_to_rownames('name') %>% 
  dplyr::select(-category,category)

levels(data$category)

col_n <- c("#a6cee3", "#1f78b4", "#b2df8a","#33a02c", "#fb9a99",
           "#e31a1c","#fdbf6f", "#ff7f00", "#cab2d6", "#62B197",
           "#E16E6D","#9392BE","#D0E7ED","#D5E4A8","#6a3d9a",
           "#4197d8","#E79397", "#f8c120","#CC88B0","#a1d5b9",
           '#083356','#ce1554', "#FFFFBF","#3288BD", "#e1abbc",
           "#6a73cf","#edd064","#0eb0c8")
names(col_n) <- unique(data$category)
                    
col_fun = list(col_1 = colorRamp2(c(-1, 0, 1), 
                                  c("#6a3d9a", "white","#e31a1c")),
               col_2 = colorRamp2(c(-1, 0, 1), 
                                  c("#6a3d9a", "white", "#e31a1c")),
               col_3 = colorRamp2(c(-1, 0, 1), 
                                  c("#6a3d9a", "white", "#e31a1c")),
               col_4 = colorRamp2(c(-1, 0, 1), 
                                  c("#6a3d9a", "white", "#e31a1c")),
               col_5 = colorRamp2(c(-1, 0, 1), 
                                  c("#6a3d9a", "white", "#e31a1c")),
               col_6 = col_n
)

#画图
circos.clear()
if(T){
pdf("8Phewas_plot.pdf", height = 30, width = 30)
circos.par$gap.degree <- 60
circos.par$start.degree <- 30
circos.par$track.margin <- c(0.001, 0.001)

for (i in 1:6) {

  data_tmp <- as.matrix(data[,i])
  if (i == 1) {
    rownames(data_tmp) <- rownames(data)
  }
  colnames(data_tmp) <- colnames(data)[i]

  if (i<6) {
    circos.heatmap(data_tmp, 
                   col = col_fun[[i]],
                   rownames.side = "outside", 
                   rownames.cex = 0.8,cluster = F,
                   cell.border = "white",show.sector.labels = T,
                   track.height = 0.025)
  } else {
    circos.heatmap(data_tmp, 
                   col = col_fun[[i]],
                   rownames.side = "outside", 
                   cluster = T,
                   track.height = 0.02)
  }
}

lgd1 <- Legend(title = "Beta", border = "black", grid_height = unit(6, "mm"),
               legend_width = unit(25, "mm"),
               at = c(-1, 0, 1), title_position = "topcenter",
               col_fun = col_fun[[1]], direction = "horizontal")
pd <- packLegend(lgd1)
draw(pd, x = unit(0.6, "npc"), y = unit(0.7, "npc"))

m=0.8
n=1.2
p=23
circos.track(track.index = get.current.track.index(),
             panel.fun = function(x, y) {
               if(CELL_META$sector.numeric.index == 1) { # the last sector
                 circos.rect(CELL_META$cell.xlim[2] + convert_x(p, "mm"), -0.5,
                             CELL_META$cell.xlim[2] + convert_x(p, "mm"), 11,border=NA)
                 
                 circos.text(CELL_META$cell.xlim[2] + convert_x(p, "mm"), m,
                             colnames(data)[6], cex = 1, facing = "inside")
                 circos.text(CELL_META$cell.xlim[2] + convert_x(p, "mm"), m+n,
                             colnames(data)[5], cex = 1, facing = "inside")
                 circos.text(CELL_META$cell.xlim[2] + convert_x(p, "mm"), m+2*n,
                             colnames(data)[4], cex = 1, facing = "inside")
                 circos.text(CELL_META$cell.xlim[2] + convert_x(p, "mm"), m+3*n,
                             colnames(data)[3], cex = 1, facing = "inside")
                 circos.text(CELL_META$cell.xlim[2] + convert_x(p, "mm"), m+4*n,
                             colnames(data)[2], cex = 1, facing = "inside")
                 circos.text(CELL_META$cell.xlim[2] + convert_x(p, "mm"), m+5*n,
                             colnames(data)[1], cex = 1, facing = "inside")
                 
               }
             }, bg.border = NA)

lgd = Legend(labels = names(col_n),
             legend_gp = gpar(fill = col_n), title = "Category", 
             ncol = 1, row_gap = unit(1, "mm"))

draw(lgd, x = unit(0.515, "npc"), y = unit(0.5, "npc"))

dev.off()
circos.clear()
} 
                  
   
                    
