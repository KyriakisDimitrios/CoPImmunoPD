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

source("Scripts/Functions.R")
setwd("C:/Users/kyria.000/Documents/PhD/Projects/Feng2023/CoPImmunoPD/")
Subbbb <- readRDS("Subset_raarranged.rds")
DefaultAssay(Subbbb) <- "RNA"


resultdir = "Result/Final_Figures/Combined_Markers/"
dir.create(resultdir)


Subbbb$CCR7cond <- as.vector(Subbbb$CCellType)
Subbbb$CCR7cond[Subbbb$CCellType %in% c("CD45RO+CCR7-","CD45RO-CCR7-") ] <- "CCR7-"
Subbbb$CCR7cond[Subbbb$CCellType %in% c("CD45RO-CCR7+","CD45RO+CCR7+")] <- "CCR7+"
DimPlot(Subbbb,group.by = "CCR7cond")
Subbbb$Combined <- factor(Subbbb$Combined ,levels=c("HC: CD45RO-CCR7+",
                                                    "PD: CD45RO-CCR7+", 
                                                    "HC: CD45RO+CCR7+","PD: CD45RO+CCR7+", "HC: CD45RO+CCR7-",
                                                    "PD: CD45RO+CCR7-", "HC: CD45RO-CCR7-",
                                                    "PD: CD45RO-CCR7-")) 

if(!file.exists("subset_ccr7_up.rds") | !file.exists("subset_ccr7_down.rds")){
  subset_ccr7_up <- subset(x = Subbbb,subset= CCR7cond=="CCR7+" ,downsample = 900)
  subset_ccr7_down <- subset(x = Subbbb,subset= CCR7cond=="CCR7-" ,downsample = 3000)
  saveRDS(subset_ccr7_up,"subset_ccr7_up.rds")
  saveRDS(subset_ccr7_down,"subset_ccr7_down.rds")
}else{
  subset_ccr7_up<- readRDS("subset_ccr7_up.rds")
  subset_ccr7_down <- readRDS("subset_ccr7_down.rds")
}



DefaultAssay(Subbbb) <- "RNA"
# ---------------------------------------- IFNG+GZMK+ -----------------------------------------
features <- c("GZMK","IFNG")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)



# ---------------------------------------- GZMA_GZMB_PRF1_IFNG -----------------------------------------
features <- c("GZMA","GZMB","PRF1","IFNG")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)



# ---------------------------------------- GZMA_GZMB_PRF1 -----------------------------------------------
features <- c("GZMA","GZMB","PRF1")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)


# ---------------------------------------- GZMA_GZMB -----------------------------------------------
features <- c("GZMA","GZMB")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)



# ---------------------------------------- "GZMA","GZMB","PRF1","FCGR3A" -----------------------------------------------
features <-  c("GZMA","GZMB","PRF1","FCGR3A")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)

# ---------------------------------------- "GZMA" "PRF1" -----------------------------------------------
features <- c("GZMA","PRF1")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)

# ---------------------------------------- "GZMB" "PRF1" -----------------------------------------------
features <- c("GZMB","PRF1")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)


# ---------------------------------------- "GZMB" "GZMK" -----------------------------------------------
features <- c("GZMB","GZMK")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)



# ---------------------------------------- "GZMA" "GZMK" -----------------------------------------------
features <- c("GZMA","GZMK")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)


# ---------------------------------------- "PRF1" "GZMK" -----------------------------------------------
features <- c("PRF1","GZMK")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)

# ---------------------------------------- "CCL5" "GZMK" -----------------------------------------------
features <- c("CCL5","GZMK")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)


# ---------------------------------------- "GZMA" "GZMB" -----------------------------------------------
features <- c("GZMA","GZMB")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)

# ---------------------------------------- "GZMA" "GZMB" PRF1 -----------------------------------------------
features <- c("GZMA","GZMB","PRF1")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)


# ---------------------------------------- "CCL5" "GZMK -----------------------------------------------
features <- c("CCL5","GZMK")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)

# ---------------------------------------- "IFNG" "GZMK -----------------------------------------------
features <- c("IFNG","GZMK")
title_genes1 <- paste0(paste(features,collapse="+"),"+")
title_genes2 <- paste(features,collapse="_")
p1 <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = features)+ggtitle(title_genes1)+xlab("")
pdf(paste0(resultdir,"Overall_Density_",title_genes2,".pdf"))
p1
dev.off()
Overall_dens_Count_Joint_plot(Seurat_CCR7pos=subset_ccr7_up,Seurat_CCR7neg=subset_ccr7_down,features=features,resultdir=resultdir)
data_df_res_all<- Combined_Percent_Expressing(seurat_object=Subbbb,features=features,group_by = "Condition",split_by = "CellType",threshold = 0)
data_df_res <- data_df_res_all$final_df
expression_info <- data_df_res_all$expression_info
as.data.frame(data_df_res)
write.csv(as.data.frame(data_df_res),paste0(resultdir,"Percentage_Exp_",title_genes2,".csv"))
expression_info_vec <- as.vector(expression_info[[title_genes1]])
names(expression_info_vec) <- rownames(x = expression_info)
Subbbb <- AddMetaData(
  object = Subbbb,
  metadata = expression_info_vec,
  col.name = title_genes2
)


# ============================================================================================================
# ============================================ DOT PLOTS =====================================================
# ============================================================================================================
source("Scripts/Functions.R")
features <- c("GZMA_GZMB","GZMA_PRF1","GZMB_PRF1",
              "GZMA_GZMB_PRF1",
              "GZMA_GZMB_PRF1_IFNG")
p1 <- myDotPlot(Subbbb,assay = "RNA",
                cols = c("lightgrey", '#14a84d'),
                features = features, 
                group.by = "Combined")+
  RotatedAxis()+
  scale_colour_gradientn(colours =c("lightgrey", '#14a84d'),breaks =  seq(0,50,5))+
  scale_size(breaks = seq(0,30,5),range = c(1, 8))
p1
ggsave(plot=p1,filename=paste0(resultdir,"DotPlot_Joint_1.pdf"))
ggsave(plot=p1+coord_flip(),filename=paste0(resultdir,"DotPlot_Joint_1_flip.pdf"))


features <- c("GZMA_GZMB","GZMA_PRF1","GZMB_PRF1",
              "GZMA_GZMB_PRF1",
              "GZMA_GZMB_PRF1_IFNG",
              "GZMA_GZMB_PRF1_FCGR3A")
p1 <- myDotPlot(Subbbb,assay = "RNA",
                cols = c("lightgrey", '#14a84d'),
                features = features, 
                group.by = "Combined")+
  RotatedAxis()+
  scale_colour_gradientn(colours =c("lightgrey", '#14a84d'),breaks =  seq(0,50,5))+
  scale_size(breaks = seq(0,30,5),range = c(1, 8))
p1
ggsave(plot=p1,filename=paste0(resultdir,"DotPlot_Joint_2.pdf"))
ggsave(plot=p1+coord_flip(),filename=paste0(resultdir,"DotPlot_Joint_2_flip.pdf"))



features <- c("GZMA_GZMK","GZMB_GZMK",
              "PRF1_GZMK",
              "CCL5_GZMK",
              "IFNG_GZMK")
p1 <- myDotPlot(Subbbb,assay = "RNA",
                cols = c("lightgrey", '#14a84d'),
                features = features, 
                group.by = "Combined")+
  RotatedAxis()+
  scale_colour_gradientn(colours =c("lightgrey", '#14a84d'),breaks =  seq(0,50,5))+
  scale_size(breaks = seq(0,30,5),range = c(1, 8))
p1
ggsave(plot=p1,filename=paste0(resultdir,"DotPlot_Joint_3.pdf"))
ggsave(plot=p1+coord_flip(),filename=paste0(resultdir,"DotPlot_Joint_3_flip.pdf"))



