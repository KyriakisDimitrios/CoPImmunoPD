require(Seurat)

# Load Object
Subbbb <- readRDS("Subset_raarranged.rds")

# Select genes
genes_to_plot <- c("GZMA", "GZMB")

title_genes <- paste0(paste(features,collapse="_"))
p <- VlnPlot(Subbbb, features=genes_to_plot,stack=T,group.by="CCellType",split.by="Condition",cols=c("black","red"),flip =T)

# Save Plot as pdf
ggsave(plot=p,filename=paste0(resultdir,"Violin_Plot",title_genes,".pdf"))
# Save Plot as png
ggsave(plot=p,filename=paste0(resultdir,"Violin_Plot",title_genes,".png"))