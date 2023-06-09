if (!require("Seurat")) install.packages("Seurat")
require(Seurat)

# Load Object
Subbbb <- readRDS("Subset_raarranged.rds")

# Select genes
genes_to_plot <- c("GZMA", "GZMB")

title_genes <- paste0(paste(features,collapse="_"))
p <- VlnPlot(Subbbb, features=genes_to_plot,stack=T,group.by="CCellType",split.by="Condition",cols=c("black","red"),flip =T)

# Save Plot as pdf
# In case you select to many genes maybe you need to adjust height or width of the saved plot
# e.g ggsave(plot=p,filename=paste0(resultdir,"Violin_Plot",title_genes,".pdf"), width=12, height=10)
ggsave(plot=p,filename=paste0(resultdir,"Violin_Plot",title_genes,".pdf"))
# Save Plot as png
ggsave(plot=p,filename=paste0(resultdir,"Violin_Plot",title_genes,".png"))