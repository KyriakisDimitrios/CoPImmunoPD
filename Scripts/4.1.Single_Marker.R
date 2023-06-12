options(future.globals.maxSize = 1000000 * 1024^2)
set.seed(2422012)

# Single cell libraries
library(Seurat)
# Rest libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(scCustomize)
library(Nebulosa)
library(dittoSeq)
library(RColorBrewer)
library(viridis )
source("Scripts/Functions.R")


setwd("C:/Users/kyria.000/Documents/PhD/Projects/Feng2023/CoPImmunoPD/")
resultdir = "Result/Final_Figures/Single_Markers/"
dir.create(resultdir)

Subbbb <- readRDS("Subset_raarranged.rds")
DefaultAssay(Subbbb) <- "RNA"

features <- c("GZMA","GZMB","PRF1","FCGR3A","IFNG","GZMK","GNLY","CCL5")
DefaultAssay(Subbbb) <- "RNA"
for(feature in features){
  p1 <-Plot_Density_Custom(Subbbb, features = feature)+ggtitle(feature)+xlab("")
  ggsave(plot=p1,filename=paste0(resultdir,"Density_",feature,".pdf"))
}

p1 <- DotPlot(Subbbb,assay = "RNA", features = features,,group.by = "Combined")+RotatedAxis()
p1 
ggsave(plot=p1,filename=paste0(resultdir,"DotPlot_Individual_Features.pdf"))
ggsave(plot=p1+coord_flip(),filename=paste0(resultdir,"DotPlot_Individual_Features_flipped.pdf"))
