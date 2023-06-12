options(future.globals.maxSize = 1000000 * 1024^2)
set.seed(2422012)

# Single cell libraries
library(Seurat)
# Rest libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(scCustomize)
library(dittoSeq)
library(RColorBrewer)
library(viridis )
library(stringr)

source("Scripts/Functions.R")
setwd("C:/Users/kyria.000/Documents/PhD/Projects/Feng2023/CoPImmunoPD/")

resultdirviolin = "./Result/Feng_Figures/Violins_Heatmaps/"


Subbbb <- readRDS("Subset_raarranged.rds")
DefaultAssay(Subbbb) <- "SCT"

# ================================== Calculate Average Expression ==========================================
avg_SCT_cl <- AverageExpression(Subbbb,group.by="Combined",return.seurat = T)
avg_SCT_cl$Cluster <- as.vector(unlist(lapply(colnames(avg_SCT_cl),function(x){str_split(x,pattern="_")[[1]][1]})))
avg_SCT_cl$Condition <- unlist(lapply(colnames(avg_SCT_cl),function(x){str_split(x,pattern="_")[[1]][2]}))
avg_SCT_cl$Cluster <- factor(avg_SCT_cl$Cluster,
    levels=c(" CD45RO-CCR7+", " CD45RO+CCR7+", " CD45RO+CCR7-", " CD45RO-CCR7-"))
avg_SCT_cl$Cluster
# -----------------------------------------------------------------------------------------------------------



# =========================
features <- c("GZMA", "GZMB", "GZMH", "GZMM", "KLRG1", "PRF1", "GNLY", "NKG7", "KLRD1")
gene_title <- paste(features,collapse="_")
Idents(Subbbb) <- "CCellType"
p <- VlnPlot(Subbbb, features,stack=T,split.by="Condition",cols=c("black","red"),flip =T)
pdf(paste0(resultdirviolin,"Violin_",gene_title,".pdf"))
plot(p)
dev.off()
p 

p <- dittoHeatmap(avg_SCT_cl,features,
    annot.colors = c(dittoColors(1)[seq_len(4)],c("black","red")),
    annot.by = c("Cluster","Condition"),cluster_cols=F)
pdf(paste0(resultdirviolin,"Heatmap_",gene_title,".pdf"))
print(p)
dev.off()
p 
# ------------------------


# =========================
features <- c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "GNLY", "NKG7", "KLRD1")
gene_title <- paste(features,collapse="_")
Idents(Subbbb) <- "CCellType"
VlnPlot(Subbbb, features,stack=T,split.by="Condition",cols=c("black","red"),flip =T)

p <- VlnPlot(Subbbb, features,stack=T,split.by="Condition",cols=c("black","red"),flip =T)
pdf(paste0(resultdirviolin,"Violin_",gene_title,".pdf"))
plot(p)
dev.off()
p 
# ------------------------


# =========================
# removed HLA-DRA, HLA-DQA1, PTPRM
graft <- c("PRF1","GZMB","KLRD1","HLA-A","HLA-DRB1","HLA-DQB1")
antigen <- c("CD74","CD8A","KLRD1","HLA-A","HLA-DRB1","TAPBP","HLA-DQB1")
celladhension <- c("SPN","CD2","CD8A","ITGB2","HLA-A","ITGAL","HLA-DRB1","HLA-DQB1")
allograft <- c("PRF1","GZMB","HLA-A","HLA-DRB1","HLA-DQB1")
autoimmune <- c("PRF1","GZMB","HLA-A","HLA-DRB1","HLA-DQB1")
nkdeat <- c("PPP3CC","ITGB2","PRF1","GZMB","KLRD1","LCP2","HLA-A","ITGAL")
features <- unique(c(graft,antigen,celladhension,allograft,autoimmune,nkdeat))

Idents(Subbbb) <- "CCellType"
p <- VlnPlot(Subbbb, features,stack=T,split.by="Condition",cols=c("black","red"),flip =T)
pdf(paste0(resultdirviolin,"Violin_selected_basedon_pathways.pdf"))
plot(p)
dev.off()
p 

p <- dittoHeatmap(avg_SCT_cl,features,
    annot.colors = c(dittoColors(1)[seq_len(4)],c("black","red")),
    annot.by = c("Cluster","Condition"),cluster_cols=F)
pdf(paste0(resultdirviolin,"Heatmap_selected_basedon_pathways.pdf"))
print(p)
dev.off()
p 
# ------------------------





# =========================
# HLA-DRA, HLA-DQA1, HLA-DQB1, ADCY9, PTPRM
graft <- c("HLA-DPB1","GZMB","HLA-C","KLRD1","HLA-A","HLA-DRB1","HLA-DPA1")
allograft <- c("HLA-DPB1","GZMB","HLA-C","HLA-A","HLA-DRB1","HLA-DPA1")
typeI <- c("HLA-DPB1","GZMB","HLA-C","HLA-A","HLA-DRB1","HLA-DPA1")
hummant <- c("ITGB2","NFATC3","HLA-C","HLA-A","ITGAL","PPP3CC","HLA-DPB1","PRKACB","HLA-DRB1","HLA-DPA1")
antigen <- c("CD74","HLA-DPB1","HLA-C","KLRD1","HLA-A","HLA-DRB1","HLA-DPA1")
autimmune <- c("HLA-DPB1","HLA-C","GZMB","HLA-A","HLA-DRB1","HLA-DPA1")
celladjension <- c("ITGB2","HLA-DPB1","HLA-C","HLA-A","ITGAL","HLA-DRB1","HLA-DPA1")
NKcell <- c("VAV3","PPP3CC","ITGB2","PLCG2","HLA-C","GZMB","KLRD1","HLA-A","ITGAL","HCST")

features <- unique(c(graft,allograft,typeI,hummant,
    antigen,autimmune,celladjension,NKcell))

Idents(Subbbb) <- "CCellType"
p <- VlnPlot(Subbbb, features,stack=T,split.by="Condition",cols=c("black","red"),flip =T)

pdf(paste0(resultdirviolin,"Violin_selected_basedon_pathways2.pdf"))
print(p)
dev.off()
p 
p<- dittoHeatmap(avg_SCT_cl,features,
    annot.colors = c(dittoColors(1)[seq_len(4)],c("black","red")),
    annot.by = c("Cluster","Condition"),cluster_cols=F)
pdf(paste0(resultdirviolin,"Heatmap_selected_basedon_pathways2.pdf"))
print(p)
dev.off()
p 
# ------------------------
