# ID - 9-4133000034498
# https://g.co/fi/r/VP23UN

options(future.globals.maxSize = 1000000 * 1024^2)
set.seed(2422012)

# Single cell libraries
library(Seurat)
library(sctransform)
library(rliger)
library(SeuratWrappers)
library(conos)
library(scater)
library(scDblFinder)
library(scran)
library(sctransform)
library(scry)
library(SoupX)

# Rest libraries
library(BiocParallel)
library(ggplot2)
library(dplyr)
library(cowplot)

library(scCustomize)
library(stringr)
library(dittoSeq)
library(DESeq2)
library(enrichR)


setwd("C:/Users/kyria.000/Documents/PhD/Projects/Feng2023/")
Subbbb <- readRDS("Subset_raarranged.rds")

table(Subbbb$Condition,Subbbb$CellType)


library(RColorBrewer)
library(viridis )
Subbbb$pseudo_order <- -Subbbb@reductions$umap@cell.embeddings[,2]
FeaturePlot(Subbbb,"pseudo_order") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(Subbbb,"pseudo_order") + scale_colour_gradientn(colours = viridis(256, option = "B"))



# ================================== Slingshot ORDERING  =======================================

library(slingshot)
sce <- as.SingleCellExperiment(Subbbb, assay = "SCT")
reducedDim(sce, "PCA", withDimnames=TRUE) <- Subbbb[['pca']]@cell.embeddings
reducedDim(sce, "UMAP", withDimnames=TRUE) <- Subbbb[['umap']]@cell.embeddings
sce.sling <- slingshot(sce, reducedDim='PCA',clusterLabels="CCellType")


embedded <- embedCurves(sce.sling, "UMAP")
embedded <- slingCurves(embedded)[[1]] # only 1 path.
embedded <- data.frame(embedded$s[embedded$ord,])
embedded

plotUMAP(sce.sling, colour_by="slingPseudotime_1") +
  geom_path(data=embedded, aes(x=UMAP_1, y=UMAP_2), size=1.2)

Subbbb$Pseudotime  <- sce.sling$slingPseudotime_1

library(sm)

library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)
# Plot

table(Subbbb$Condition)


## An example with non-unique breaks:
x <- seq(0,10,1)
b <- seq(0,10,2)
x
.bincode(x, b, TRUE)
.bincode(x, b, FALSE)
.bincode(x, b, TRUE, TRUE)
.bincode(x, b, FALSE, TRUE)


Subbbb$Binned <- factor(.bincode(Subbbb$Pseudotime, breaks=seq(0,100,10), TRUE, TRUE),levels=c(1:10))
levels(Subbbb$Binned) <- c(" [0-10]","(10-20]","(20-30]","(30-40]","(40-50]",
                           "(50-60]","(60-70]","(70-80]","(80-90]","(90-100]")
Subbbb$Pseudotime[which (Subbbb$Pseudotime<1)[1:10]]
Subbbb$Binned[which (Subbbb$Pseudotime<1)[1:10]]


Subbbb$Binned20 <- factor(.bincode(Subbbb$Pseudotime, breaks=seq(0,100,5), TRUE, TRUE))
levels(Subbbb$Binned20) <- c(" [0-05]", "(05-10]",
                             "(10-15]","(15-20]",
                             "(20-25]","(25-30]",
                             "(30-35]","(35-40]",
                             "(40-45]","(45-50]",
                             "(50-55]","(55-60]",
                             "(60-65]","(65-70]",
                             "(70-75]","(75-80]",
                             "(80-85]","(85-90]",
                             "(90-95]","(95-100]")


pt_20 <- prop.table(table(Subbbb$Binned20,Subbbb$Condition),margin = 2)
write.ftable(ftable(pt_20),file = "Result/Pseudotime/Prob_Binned_20.csv", sep = ",",
             quote = FALSE)


pt_10 <- prop.table(table(Subbbb$Binned,Subbbb$Condition),margin = 2)
write.ftable(ftable(pt_10),file = "Result/Pseudotime/Prob_Binned_10.csv", sep = ",",
             quote = FALSE)



Subbbb$Pseudotime[which (Subbbb$Pseudotime<1)[1:10]]
Subbbb$Binned20[which (Subbbb$Pseudotime<1)[1:10]]


Subbbb$revBinned10 <- factor(Subbbb$Binned, levels= rev(levels(Subbbb$Binned))) 
Subbbb$revBinned20 <- factor(Subbbb$Binned20, levels= rev(levels(Subbbb$Binned20))) 



pdf("Result/Pseudotime/Bar_Stack_Pseudotime.pdf",width=12)
dittoBarPlot(Subbbb , "revBinned10", 
             group.by = "Condition",
             legend.show=TRUE,retain.factor.levels = TRUE,x.labels.rotate = TRUE
)+scale_fill_viridis_d(option = "plasma",direction = -1)+coord_flip()+ggtitle("")+ theme_cowplot()

dittoBarPlot(Subbbb , "revBinned20", 
             group.by = "Condition",
             legend.show=TRUE,retain.factor.levels = TRUE,x.labels.rotate = TRUE
)+scale_fill_viridis_d(option = "plasma",direction = -1)+coord_flip()+ggtitle("")+ theme_cowplot()

dittoBarPlot(Subbbb , "Condition", 
             group.by = "Binned",
             color.panel =c("black","red")
)+ theme_cowplot()+RotatedAxis()+ggtitle("")
dittoBarPlot(Subbbb , "Condition", 
             group.by = "Binned20",
             color.panel =c("black","red")
)+RotatedAxis()+ggtitle("")
dev.off()






dittoBarPlot(Subbbb , "revBinned", 
             group.by = "Condition",x.labels.rotate = TRUE,retain.factor.levels = TRUE,
)+coord_flip()+ scale_fill_viridis_d(option = "plasma",direction = -1)


ggplot(Subbbb@meta.data, aes(x = Pseudotime,color=Condition))  +
  geom_histogram(binwidth=1.5)


pdf("Result/Pseudotime/Ridge_Pseudotime.pdf",width=9,height=7)
ggplot(Subbbb@meta.data, aes(x = Pseudotime, y = CCellType, fill = ..x..))  +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.0) +
  facet_wrap(Subbbb$Condition ~ .,ncol=1,scale="free_y") +
  scale_fill_viridis(name = "Pseudotime", option = "C") +
  labs(title = 'Pseudotime')+theme_cowplot()
ggplot(Subbbb@meta.data, aes(x = Pseudotime, y = CCellType, fill = ..x..))  +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  facet_wrap(Subbbb$Condition ~ .,ncol=1,scale="free_y") +
  scale_fill_viridis(name = "Pseudotime", option = "C") +
  labs(title = 'Pseudotime')+theme_cowplot()
dev.off()

pdf("Result/Pseudotime/Ridge_Pseudotime_Overlays.pdf",width=12)
ggplot(Subbbb@meta.data, aes(x=Pseudotime, fill=Condition, colour = Condition)) +
  geom_density(alpha=.1,bw=1.5)+#y+
  scale_fill_manual(values=c('black','red'))+
  scale_colour_manual(values=c('black','red'))+
  facet_wrap(Subbbb$CCellType ~ .,ncol=2,scale="free_y") +
  labs(title = 'Pseudotime')+theme_ggprism_mod()+ggtitle("")

ggplot(Subbbb@meta.data, aes(x = Pseudotime, y = Condition, fill = Condition))  +
  geom_density_ridges_gradient(rel_min_height = 0.0)+#y+
  scale_fill_manual(values=c('black','red'))+
  facet_wrap(Subbbb$CCellType ~ .,ncol=2,scale="free_y") +
  labs(title = 'Pseudotime')+theme_ggprism_mod()+ggtitle("")

ggplot(Subbbb@meta.data, aes(x = Pseudotime, y = Condition, fill = Condition))  +
  geom_density_ridges_gradient(rel_min_height = 0.01)+#y+
  scale_fill_manual(values=c('black','red'))+
  facet_wrap(Subbbb$CCellType ~ .,ncol=2,scale="free_y") +
  labs(title = 'Pseudotime')+theme_ggprism_mod()+ggtitle("")
dev.off()








ggplot(Subbbb@meta.data, aes(x = Pseudotime, y = Condition, fill = Condition))  +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01)+
  scale_fill_manual(values=c('black','red'))+
  facet_wrap(Subbbb$CCellType ~ .,ncol=2,scale="free_y") +
  labs(title = 'Pseudotime')+theme_cowplot()


ggplot(Subbbb@meta.data, aes(x = Pseudotime, y = CCellType, fill = ..x..))  +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  facet_wrap(Subbbb$Condition ~ .,ncol=1) +
  scale_fill_viridis(name = "Temp. [F]", option = "C") +
  labs(title = 'Pseudotime') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )

RidgePlot(Subbbb,"Pseudotime",
          group.by = "CCellType",fill.by="feature")+
  NoLegend()+
  scale_fill_viridis(name = "Temp. [F]", option = "C") +
  facet_wrap(as.factor(Subbbb$Condition),ncol=1)



ggplot(data=Subbbb@meta.data, aes(x=Pseudotime, group=Condition,fill=Condition)) +
  geom_density(adjust=1.5, alpha=.4) +facet_wrap(as.factor(Subbbb$CCellType))


FeaturePlot(Subbbb,c("Pseudotime"),split.by = "Condition") & scale_colour_gradientn(colours = viridis(256, option = "B"))


p <- FeaturePlot(Subbbb,c("Pseudotime")) + scale_colour_gradientn(colours = viridis(256, option = "B"))
pdf("Result/Pseudotime/Umap_Pseudotime.pdf")
p 
dev.off()
FeaturePlot(Subbbb,c("pseudo_order","Pseudotime")) + scale_colour_gradientn(colours = viridis(256, option = "B"))


p <- RidgePlot(Subbbb,"Pseudotime",
               group.by = "Condition",
               cols = c('black','red'))+
  NoLegend()+
  facet_wrap(Subbbb$CCellType,scale="free_y")

pdf("Result/Pseudotime/Ridge_Pseudotime.pdf")
RidgePlot(Subbbb,"Pseudotime",
          group.by = "CCellType")+
  NoLegend()+
  facet_wrap(as.factor(Subbbb$Condition),scale="free_y",ncol=1)
p
dev.off()


# ================================================================================================



# ======= UMAP ORDERING 
RidgePlot(Subbbb,"pseudo_order",
          group.by = "Condition",
          cols = c('black','red'))+
  NoLegend()+
  facet_wrap(Subbbb$CCellType,scale="free")
RidgePlot(Subbbb,"pseudo_order",
          group.by = "CCellType")+
  NoLegend()+
  facet_wrap(as.factor(Subbbb$Condition),scale="free",ncol=1)
# ===================== 

# ============================== Calcualte geometric mean ================================================
unique(Subbbb$CCellType)

mean1_hc <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO-CCR7+" & Subbbb$Condition=="HC"  & Subbbb$Pseudotime >0])))
mean1_pd <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO-CCR7+" & Subbbb$Condition=="PD"  & Subbbb$Pseudotime >0])))

mean2_hc <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO+CCR7+" & Subbbb$Condition=="HC"  & Subbbb$Pseudotime >0])))
mean2_pd <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO+CCR7+" & Subbbb$Condition=="PD" & Subbbb$Pseudotime >0])))

mean3_hc <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO+CCR7-" & Subbbb$Condition=="HC" & Subbbb$Pseudotime >0])))
mean3_pd <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO+CCR7-" & Subbbb$Condition=="PD" & Subbbb$Pseudotime >0])))

mean4_hc <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO-CCR7-" & Subbbb$Condition=="HC" & Subbbb$Pseudotime >0])))
mean4_pd <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO-CCR7-" & Subbbb$Condition=="PD" & Subbbb$Pseudotime >0])))

Subbbb$Pseudotime <- Subbbb$Pseudotime+0.01

mean1_hc <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO-CCR7+" & Subbbb$Condition=="HC"])))
mean1_pd <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO-CCR7+" & Subbbb$Condition=="PD" ])))

mean2_hc <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO+CCR7+" & Subbbb$Condition=="HC" ])))
mean2_pd <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO+CCR7+" & Subbbb$Condition=="PD" ])))

mean3_hc <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO+CCR7-" & Subbbb$Condition=="HC" ])))
mean3_pd <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO+CCR7-" & Subbbb$Condition=="PD"])))

mean4_hc <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO-CCR7-" & Subbbb$Condition=="HC" ])))
mean4_pd <- exp(mean(log(Subbbb$Pseudotime[Subbbb$CCellType=="CD45RO-CCR7-" & Subbbb$Condition=="PD" ])))


gm_df <- data.frame(CellType = unique(Subbbb$Combined),
           GM = c(mean2_hc,mean1_hc,mean4_hc,mean3_hc,
             mean2_pd,mean1_pd,mean4_pd,mean3_pd),
           Condition = c(rep("HC",4),rep("PD",4)), 
           CellTypes = rep(c("CD45RO+CCR7+","CD45RO-CCR7+","CD45RO-CCR7-","CD45RO+CCR7-"),2)
           )
gm_df$CellTypes <- factor(gm_df$CellTypes,levels=levels(Subbbb$CCellType))

pdf("Result/Pseudotime/GM_Pseudotime.pdf")
ggplot(gm_df,aes(x=CellTypes,y=GM,color=Condition,fill=Condition,group=Condition))+
  geom_point()+
  geom_line()+theme_cowplot()+RotatedAxis()+scale_colour_manual(values =c('black','red'))
dev.off()
# ===========================================================================================================



pdens <- Plot_Density_Joint_Only(seurat_object = Subbbb, features = c("GZMA","GZMB","PRF1","FCGR3A"),viridis_palette="inferno")



# ==================== DEGs on overall Pseudotime ============================
sce$Pseudotime <- Subbbb$Pseudotime
pseudo <- TSCAN::testPseudotime(sce, pseudotime=sce$Pseudotime)
pseudo$SYMBOL <- rowData(sce)$SYMBOL
pseudo_sign <- pseudo[pseudo$p.value < 0.01,]
up.right <- pseudo_sign[pseudo_sign$logFC > 0,]
up.right_ord <- up.right[order(abs(up.right$logFC),decreasing = T),]
up.right_ord$SYMBOL <- rownames(up.right_ord)
up.left <- pseudo_sign[pseudo_sign$logFC < 0,]
up.left_ord <- up.left[order(abs(up.left$logFC),decreasing = T),]
up.left_ord$SYMBOL <- rownames(up.left_ord)
best_left <- head(up.left_ord$SYMBOL, 20)
best_right <- head(up.right_ord$SYMBOL, 20)
library(viridis)
p <- plotHeatmap(sce,features=c(best_left,best_right), 
            order_columns_by='Pseudotime',
            center=TRUE,column_annotation_colors= list(Pseudotime=viridis(256, option = "B")))
pdf("Result/Pseudotime/Heatmap_top20_Pseudotime.pdf")
p
dev.off()
# ==================== DEGs on HC Pseudotime ============================
sce_hc <- subset(sce, , Condition=="HC")
pseudo <- TSCAN::testPseudotime(sce_hc, pseudotime=sce_hc$Pseudotime)
pseudo$SYMBOL <- rowData(sce_hc)$SYMBOL
pseudo_sign_hc <- pseudo[pseudo$p.value < 0.01,]
up.right <- pseudo_sign_hc[pseudo_sign_hc$logFC > 0,]
up.right_ord <- up.right[order(abs(up.right$logFC),decreasing = T),]
up.right_ord$SYMBOL <- rownames(up.right_ord)
up.left <- pseudo_sign_hc[pseudo_sign_hc$logFC < 0,]
up.left_ord <- up.left[order(abs(up.left$logFC),decreasing = T),]
up.left_ord$SYMBOL <- rownames(up.left_ord)
best_left <- head(up.left_ord$SYMBOL, 20)
best_right <- head(up.right_ord$SYMBOL, 20)
library(viridis)
p <- plotHeatmap(sce_hc,features=c(best_left,best_right), 
                 order_columns_by='Pseudotime',
                 center=TRUE,column_annotation_colors= list(Pseudotime=viridis(256, option = "B")))
p
pdf("Result/Pseudotime/Heatmap_top20_Pseudotime_HC.pdf")
p
dev.off()
# ==================== DEGs on HC Pseudotime ============================
sce_pd <- subset(sce, , Condition=="PD")
pseudo <- TSCAN::testPseudotime(sce_pd, pseudotime=sce_pd$Pseudotime)
pseudo$SYMBOL <- rowData(sce_pd)$SYMBOL
pseudo_sign_pd <- pseudo[pseudo$p.value < 0.01,]
up.right <- pseudo_sign_pd[pseudo_sign_pd$logFC > 0,]
up.right_ord <- up.right[order(abs(up.right$logFC),decreasing = T),]
up.right_ord$SYMBOL <- rownames(up.right_ord)
up.left <- pseudo_sign_pd[pseudo_sign_pd$logFC < 0,]
up.left_ord <- up.left[order(abs(up.left$logFC),decreasing = T),]
up.left_ord$SYMBOL <- rownames(up.left_ord)
best_left <- head(up.left_ord$SYMBOL, 20)
best_right <- head(up.right_ord$SYMBOL, 20)
p <- plotHeatmap(sce_pd,features=c(best_left,best_right), 
                 order_columns_by='Pseudotime',
                 center=TRUE,column_annotation_colors= list(Pseudotime=viridis(256, option = "B")))
p
pdf("Result/Pseudotime/Heatmap_top20_Pseudotime_PD.pdf")
p
dev.off()

rownames(pseudo_sign_pd)[ ]
# ================================= VENNN COMMON ==============================
library(ggvenn)
library(dplyr)

pseudo_sign_hc <- as.data.frame(pseudo_sign_hc)
pseudo_sign_pd <- as.data.frame(pseudo_sign_pd)


x <- list( 
  HC = rownames(pseudo_sign_hc %>% filter(  FDR <0.05)),
  PD = rownames(pseudo_sign_pd %>% filter(FDR <0.05))
)
x

p <- ggvenn(
  x, 
  fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
  stroke_size = 0.5, set_name_size = 4
)
p




unique_degs_pd <- x$PD[! x$PD %in% x$HC]
pseudo_degs_pd <- as.data.frame(enrichR::enrichr(genes = unique_degs_pd, databases = "Reactome_2022"))
head(pseudo_degs_pd)
pseudo_degs_pd$log10pval <- -log10(x = pseudo_degs_pd[, paste("Reactome_2022", sep = ".", "P.value")])
pseudo_degs_pd$term <- pseudo_degs_pd$Reactome_2022.Term
pseudo_degs_pd$term <- factor(pseudo_degs_pd$term,levels=pseudo_degs_pd$term[order(pseudo_degs_pd$log10pval)])
pseudo_degs_pd <- pseudo_degs_pd[order(pseudo_degs_pd$log10pval,decreasing = T),]

unique_degs_hc <- x$HC[! x$HC %in% x$PD]
pseudo_degs_hc <- as.data.frame(enrichR::enrichr(genes = unique_degs_hc, databases = "Reactome_2022"))
head(pseudo_degs_hc)
pseudo_degs_hc$log10pval <- -log10(x = pseudo_degs_hc[, paste("Reactome_2022", sep = ".", "P.value")])
pseudo_degs_hc$term <- pseudo_degs_hc$Reactome_2022.Term
pseudo_degs_hc$term <- factor(pseudo_degs_hc$term,levels=pseudo_degs_hc$term[order(pseudo_degs_hc$log10pval)])
pseudo_degs_hc <- pseudo_degs_hc[order(pseudo_degs_hc$log10pval,decreasing = T),]

text_col <- 'term'


source("Functions.R")
p <- plot_enrichR_wrap_fun(
  ident.1="PD",
  ident.2="HC",
  pos.er = pseudo_degs_pd,
  neg.er = pseudo_degs_hc,
  enrich.database="Reactome_2022" )


print(p)




ggplot(data = pseudo_degs_pd[1:20,], aes_string(x = text_col, y = "log10pval"))+
  geom_bar(stat = "identity", fill = "dodgerblue") +
  coord_flip() + xlab("Pathway") +
  scale_fill_manual(values = c("black","blue"), drop = FALSE) +
  ylab("-log10(pval)") +
  ggtitle(paste('Reactome_2022', 'PD', sep = "_", "positive markers")) +
  theme_classic() +
  geom_text(aes_string(label = text_col, y = 0),
            size = 4,
            color = "black",
            position = position_dodge(1),
            hjust = 0)+
  theme(axis.title.y= element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())








pseudo_degs <- as.data.frame(enrichR::enrichr(genes = c(best_left,best_right), databases = "Reactome_2022"))
head(pseudo_degs)
pseudo_degs$log10pval <- -log10(x = pseudo_degs[, paste("Reactome_2022", sep = ".", "P.value")])
pseudo_degs$term <- factor(pseudo_degs$term,levels=pseudo_degs$term[order(pseudo_degs$log10pval)])
pseudo_degs <- pseudo_degs[order(pseudo_degs$log10pval),]


pseudo_degs$term <- factor(pseudo_degs$Reactome_2022.Term,levels=pseudo_degs$Reactome_2022.Term[order(pseudo_degs$log10pval)])

pseudo_degs[1:20, ]

ggplot(data = pseudo_degs[1:20, ], aes_string(x = 'term', y = "log10pval"))+
  geom_bar(stat = "identity", fill = "dodgerblue") +
  coord_flip() + xlab("Pathway") +
  scale_fill_manual(values = c("black","blue"), drop = FALSE) +
  ylab("-log10(pval)") +
  theme_classic() +
  geom_text(aes_string(label ='term' , y = 0),
            size = 4,
            color = "black",
            position = position_dodge(1),
            hjust = 0)+
  theme(axis.title.y= element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())








unique_degs_pd <- x$PD[! x$PD %in% x$HC]
pseudo_degs_pd <- as.data.frame(enrichR::enrichr(genes = unique_degs_pd, databases = "KEGG_2021_Human"))
head(pseudo_degs_pd)
pseudo_degs_pd$log10pval <- -log10(x = pseudo_degs_pd[, paste("KEGG_2021_Human", sep = ".", "P.value")])
pseudo_degs_pd$term <- pseudo_degs_pd$KEGG_2021_Human.Term
pseudo_degs_pd$term <- factor(pseudo_degs_pd$term,levels=pseudo_degs_pd$term[order(pseudo_degs_pd$log10pval)])
pseudo_degs_pd <- pseudo_degs_pd[order(pseudo_degs_pd$log10pval,decreasing = T),]

unique_degs_hc <- x$HC[! x$HC %in% x$PD]
pseudo_degs_hc <- as.data.frame(enrichR::enrichr(genes = unique_degs_hc, databases = "KEGG_2021_Human"))
head(pseudo_degs_hc)
pseudo_degs_hc$log10pval <- -log10(x = pseudo_degs_hc[, paste("KEGG_2021_Human", sep = ".", "P.value")])
pseudo_degs_hc$term <- pseudo_degs_hc$KEGG_2021_Human.Term
pseudo_degs_hc$term <- factor(pseudo_degs_hc$term,levels=pseudo_degs_hc$term[order(pseudo_degs_hc$log10pval)])
pseudo_degs_hc <- pseudo_degs_hc[order(pseudo_degs_hc$log10pval,decreasing = T),]

text_col <- 'term'


source("Functions.R")
p <- plot_enrichR_wrap_fun(
  ident.1="PD",
  ident.2="HC",
  pos.er = pseudo_degs_pd,
  neg.er = pseudo_degs_hc,
  enrich.database="KEGG_2021_Human" )
print(p)








plotExpression(sce, features=best)

on.first.path <- !is.na(sce$TSCAN.first)





Subbbb <- FindVariableFeatures(Subbbb,assay = "integrated",
                               nfeatures=100, do.plot = FALSE)

var_genes <- Subbbb@assays$integrated@var.features

library(tradeSeq)

nonna.pseudo <- pathStat(sce.sling)
not.on.path <- is.na(nonna.pseudo)
nonna.pseudo[not.on.path] <- 0
cell.weights <- !not.on.path
storage.mode(cell.weights) <- "numeric"



fit <- fitGAM(Subbbb@assays$SCT@data[var_genes,], 
              pseudotime=Subbbb$Pseudotime)+ scale_colour_gradientn(colours = viridis(256, option = "B"))





