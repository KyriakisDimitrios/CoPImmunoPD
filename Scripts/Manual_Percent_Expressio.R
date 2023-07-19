
require(ggplot2)
require(Seurat)
library(scCustomize)

setwd("C:/Users/kyria.000/Documents/PhD/Projects/Feng2023/CoPImmunoPD")
# Load Object
Subbbb <- readRDS("Subset_raarranged.rds")

Subbbb


table(Subbbb$CellType,Subbbb$Condition)
median(Subbbb$nFeature_RNA)
median(Subbbb$nCount_RNA)


# Select genes
Subbbb$Combined <- paste0(Subbbb$Condition,"_",Subbbb$CellType)

genes_to_plot <- c("TIGIT","EOMES","TOX","HAVCR2","MAF","LAG3","CTLA4")
percent_express_RNA <- Percent_Expressing(seurat_object = Subbbb,
    assay="RNA",
    features = genes_to_plot,
    group_by = "Combined")

percent_express_SCT <- Percent_Expressing(seurat_object = Subbbb,
    assay="SCT",
    features = genes_to_plot,
    group_by = "Combined")

write.csv(as.data.frame(percent_express_RNA),"Percent_Expression_extra_markers.csv")


percent_express_RNA

percent_express_SCT

VlnPlot(Subbbb,group.by="Combined",features=genes_to_plot,stack=T)



