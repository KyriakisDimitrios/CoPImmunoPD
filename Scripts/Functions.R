

gene_score_calc <- function(dataframe_obj){
    dataframe_obj <- dataframe_obj%>% mutate(
                cluster = ifelse(avg_log2FC >0 , "PD" ,"HC"), 
                enrichment.ratio = ifelse(pct.1>pct.2, pct.1 / (pct.2 + .000001), pct.2 / (pct.1 + .000001)), 
                diff.pct = pct.1 - pct.2,
                tstat= -log10(p_val_adj + 1e-320) * sign(avg_log2FC) ,
                gene.score = avg_log2FC * enrichment.ratio) %>% 
                dplyr::select(gene,tstat, gene.score,p_val, p_val_adj, diff.pct, enrichment.ratio,
                    pct.1, pct.2, avg_log2FC,cluster) %>%
                mutate(across(where(is.numeric), .f = ~ round(.x, 5)))
    return(dataframe_obj)
}





DEenrichRPlot_Man <- function(
    object,
    ident.1 = NULL,
    ident.2 = NULL,
    balanced = TRUE,
    logfc.threshold.return_all=0,
    logfc.threshold = 0.25,
    p.val.cutoff = 0.05,
    assay = NULL,
    max.genes,
    test.use = 'wilcox',
    cols = NULL,
    enrich.database = NULL,
    num.pathway = 10,
    return.gene.list = FALSE,
    ...
) {
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay

    all.markers <- FindMarkers(
        object = object,
        ident.1 = ident.1,
        ident.2 = ident.2,
        only.pos = FALSE,
        logfc.threshold = logfc.threshold.return_all,
        test.use = test.use,
        assay = assay
    )

    pos.markers <- all.markers[all.markers[, 2] > logfc.threshold & all.markers[, 1] < p.val.cutoff, , drop = FALSE]
    pos.markers.list <- rownames(x = pos.markers)[1:min(max.genes, nrow(x = pos.markers))]
    pos.er <- enrichR::enrichr(genes = pos.markers.list, databases = enrich.database)
    pos.er <- do.call(what = cbind, args = pos.er)
    pos.er$log10pval <- -log10(x = pos.er[, paste(enrich.database, sep = ".", "P.value")])
    pos.er$term <- pos.er[, paste(enrich.database, sep = ".", "Term")]
    pos.er <- pos.er[1:num.pathway, ]
    pos.er$term <- factor(x = pos.er$term, levels = pos.er$term[order(pos.er$log10pval)])
    gene.list <- list(pos = pos.er)

    neg.markers <- all.markers[all.markers[, 2] < -logfc.threshold & all.markers[, 1] < p.val.cutoff, , drop = FALSE]
    neg.markers.list <- rownames(x = neg.markers)[1:min(max.genes, nrow(x = neg.markers))]
    neg.er <- enrichR::enrichr(genes = neg.markers.list, databases = enrich.database)
    neg.er <- do.call(what = cbind, args = neg.er)
    neg.er$log10pval <- -log10(x = neg.er[, paste(enrich.database, sep = ".", "P.value")])
    neg.er$term <- neg.er[, paste(enrich.database, sep = ".", "Term")]
    neg.er <- neg.er[1:num.pathway, ]
    neg.er$term <- factor(x = neg.er$term, levels = neg.er$term[order(neg.er$log10pval)])

    if(isTRUE(length(neg.er$term) == 0) & isTRUE(length(pos.er == 0))){
        stop("No positive or negative marker genes identified")
    }

    else{
        if(isTRUE(length(neg.er$term) == 0)){
        gene.list <- list(pos = pos.er)
        }
        else{
          gene.list <- list(pos = pos.er, neg = neg.er)
        }
    }


    if(nrow(pos.markers) == 0){
        message("No positive markers to plot")

        if (isTRUE(x = balanced)) {

        p2 <- ggplot(data = neg.er, aes_string(x = "term", y = "log10pval")) +
            geom_bar(stat = "identity", fill = "indianred2") +
            coord_flip() + xlab("Pathway") +
            scale_fill_manual(values = cols, drop = FALSE) +
            ylab("-log10(pval)") +
            ggtitle(paste(enrich.database, ident.1, sep = "_", "negative markers")) +
            theme_classic() +
            geom_text(aes_string(label = "term", y = 0),
                    size = 5,
                    color = "black",
                    position = position_dodge(1),
                    hjust = 0)+
            theme(axis.title.y= element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank())
        p <- p2

        }
        else{
        stop("Nothing to plot")
        }
    }

    p1 <- ggplot(data = pos.er, aes_string(x = "term", y = "log10pval"))+
        geom_bar(stat = "identity", fill = "dodgerblue") +
        coord_flip() + xlab("Pathway") +
        scale_fill_manual(values = c("black","blue"), drop = FALSE) +
        ylab("-log10(pval)") +
        ggtitle(paste(enrich.database, ident.1, sep = "_", "positive markers")) +
        theme_classic() +
        geom_text(aes_string(label = "term", y = 0),
                size = 5,
                color = "black",
                position = position_dodge(1),
                hjust = 0)+
        theme(axis.title.y= element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())

    p2 <- ggplot(data = neg.er, aes_string(x = "term", y = "log10pval")) +
        geom_bar(stat = "identity", fill = "indianred2") +
        coord_flip() + xlab("Pathway") +
        ylab("-log10(pval)") +
        ggtitle(paste(enrich.database, ident.1, sep = "_", "negative markers")) +
        theme_classic() +
        geom_text(aes_string(label = "term", y = 0),
                  size = 5,
                  color = "black",
                  position = position_dodge(1),
                  hjust = 0)+
        theme(axis.title.y= element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

    

    p <- p1+p2




    all.markers$gene <- rownames(all.markers)
    all.markers$cluster <- unlist(lapply(abs(all.markers$avg_log2FC),function(x){if (  x > 1 ) { "PD" } else { "HC" }}))

    all.markers  <- all.markers  %>% mutate(enrichment.ratio = pct.1 / (pct.2 + .000001),
                    diff.pct = pct.1 - pct.2,
                    gene.score = avg_log2FC * enrichment.ratio)
    all.markers  <- all.markers  %>% 
        dplyr::select(gene, gene.score, p_val_adj, diff.pct, enrichment.ratio,
                        pct.1, pct.2, avg_log2FC,cluster) %>%
            mutate(across(where(is.numeric), .f = ~ round(.x, 5)))
    all.markers  <- all.markers  %>% group_by(cluster) %>%
                arrange(desc(gene.score), .by_group = TRUE)

    return(list(plot=p,
        DEG_mat = all.markers ,
        enrich_pos_table=pos.er,
        enrich_neg_table=neg.er))
}



enrichR_wrap_fun_pos <- function(all.markers,
        logfc.threshold=0.25,
        p.val.cutoff=0.05,
        num.pathway = 10,max.genes =500,
        enrich.database)
    {
    pos.markers <- all.markers[all.markers[,"avg_log2FC"] > logfc.threshold & all.markers[,"p_val"] < p.val.cutoff, , drop = FALSE]
    pos.markers.list <- rownames(x = pos.markers)[1:min(max.genes, nrow(x = pos.markers))]
    pos.er <- enrichR::enrichr(genes = pos.markers.list, databases = enrich.database)
    pos.er <- do.call(what = cbind, args = pos.er)
    pos.er$log10pval <- -log10(x = pos.er[, paste(enrich.database, sep = ".", "P.value")])
    pos.er$term <- pos.er[, paste(enrich.database, sep = ".", "Term")]
    pos.er$term <- factor(x = pos.er$term, levels = pos.er$term[order(pos.er$log10pval)])
    pos.er[,2] <- paste0("ratio ",pos.er[,2])
    return(list(pos_mat = pos.er,
                pos.markers=pos.markers))
}


enrichR_wrap_fun_neg<- function(all.markers,
        logfc.threshold=0.25,
        p.val.cutoff=0.05,
        num.pathway = 10,max.genes =500,
        enrich.database)
    {
    neg.markers <- all.markers[all.markers[,"avg_log2FC"] < -logfc.threshold & all.markers[,"p_val"] < p.val.cutoff, , drop = FALSE]
    neg.markers.list <- rownames(x = neg.markers)[1:min(max.genes, nrow(x = neg.markers))]
    neg.er <- enrichR::enrichr(genes = neg.markers.list, databases = enrich.database)
    neg.er <- do.call(what = cbind, args = neg.er)
    neg.er$log10pval <- -log10(x = neg.er[, paste(enrich.database, sep = ".", "P.value")])
    neg.er$term <- neg.er[, paste(enrich.database, sep = ".", "Term")]
    neg.er$term <- factor(x = neg.er$term, levels = neg.er$term[order(neg.er$log10pval)])
    neg.er[,2] <- paste0("ratio ",neg.er[,2])
    return(list(neg_mat = neg.er,
                neg.markers=neg.markers))
}



enrichR_wrap_fun <- function(all.markers,
        logfc.threshold=0.25,
        p.val.cutoff=0.05,
        num.pathway = 10,max.genes =500,
        enrich.database)
    {
    pos.markers <- all.markers[all.markers[,"avg_log2FC"] > logfc.threshold & all.markers[,"p_val"] < p.val.cutoff, , drop = FALSE]
    pos.markers.list <- rownames(x = pos.markers)[1:min(max.genes, nrow(x = pos.markers))]
    pos.er <- enrichR::enrichr(genes = pos.markers.list, databases = enrich.database)
    pos.er <- do.call(what = cbind, args = pos.er)
    pos.er$log10pval <- -log10(x = pos.er[, paste(enrich.database, sep = ".", "P.value")])
    pos.er$term <- pos.er[, paste(enrich.database, sep = ".", "Term")]
    pos.er$term <- factor(x = pos.er$term, levels = pos.er$term[order(pos.er$log10pval)])
    pos.er[,2] <- paste0("ratio ",pos.er[,2])

    neg.markers <- all.markers[all.markers[,"avg_log2FC"] < -logfc.threshold & all.markers[,"p_val"] < p.val.cutoff, , drop = FALSE]
    neg.markers.list <- rownames(x = neg.markers)[1:min(max.genes, nrow(x = neg.markers))]
    neg.er <- enrichR::enrichr(genes = neg.markers.list, databases = enrich.database)
    neg.er <- do.call(what = cbind, args = neg.er)
    neg.er$log10pval <- -log10(x = neg.er[, paste(enrich.database, sep = ".", "P.value")])
    neg.er$term <- neg.er[, paste(enrich.database, sep = ".", "Term")]
    neg.er$term <- factor(x = neg.er$term, levels = neg.er$term[order(neg.er$log10pval)])
    neg.er[,2] <- paste0("ratio ",neg.er[,2])

    return(list(pos_mat = pos.er,
                pos.markers=pos.markers,
                neg_mat = neg.er,
                neg.markers=neg.markers))
}











plot_enrichR_wrap_fun <- function(ident.1,ident.2,pos.er,neg.er,enrich.database,num.pathway=10){
    text_col <-  "term"
    # if(enrich.database=="Reactome_2022"){
    #     text_col <- "Term_text"
    # }
    pos.er <- pos.er[1:num.pathway, ]
    neg.er <- neg.er[1:num.pathway, ]

    pos.er$term <- factor(pos.er$term,levels=pos.er$term[order(pos.er$log10pval)])
    neg.er$term <- factor(neg.er$term,levels=neg.er$term[order(neg.er$log10pval)])

    gene.list <- list(pos = pos.er)
    if(isTRUE(length(neg.er$term) == 0) & isTRUE(length(pos.er == 0))){
        stop("No positive or negative marker genes identified")
    }
    else{
        if(isTRUE(length(neg.er$term) == 0)){
            gene.list <- list(pos = pos.er)
        }
        else{
          gene.list <- list(pos = pos.er, neg = neg.er)
        }
    }
    if(nrow(pos.er) == 0){
        message("No positive markers to plot")

        if (isTRUE(x = balanced)) {

            p2 <- ggplot(data = neg.er, aes_string(x = text_col, y = "log10pval")) +
                geom_bar(stat = "identity", fill = "indianred2") +
                coord_flip() + xlab("Pathway") +
                scale_fill_manual(values = cols, drop = FALSE) +
                ylab("-log10(pval)") +
                ggtitle(paste(enrich.database, ident.1, sep = "_", "negative markers")) +
                theme_classic() +
                geom_text(aes_string(label = text_col, y = 0),
                        size = 4,
                        color = "black",
                        position = position_dodge(1),
                        hjust = 0)+
                theme(axis.title.y= element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank())
            p <- p2
        }
        else{
        stop("Nothing to plot")
        }
    }

    p1 <- ggplot(data = pos.er, aes_string(x = text_col, y = "log10pval"))+
        geom_bar(stat = "identity", fill = "dodgerblue") +
        coord_flip() + xlab("Pathway") +
        scale_fill_manual(values = c("black","blue"), drop = FALSE) +
        ylab("-log10(pval)") +
        ggtitle(paste(enrich.database, ident.1, sep = "_", "positive markers")) +
        theme_classic() +
        geom_text(aes_string(label = text_col, y = 0),
                size = 4,
                color = "black",
                position = position_dodge(1),
                hjust = 0)+
        theme(axis.title.y= element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())

    p2 <- ggplot(data = neg.er, aes_string(x = text_col, y = "log10pval")) +
        geom_bar(stat = "identity", fill = "indianred2") +
        coord_flip() + xlab("Pathway") +
        ylab("-log10(pval)") +
        ggtitle(paste(enrich.database, ident.1, sep = "_", "negative markers")) +
        theme_classic() +
        geom_text(aes_string(label = text_col, y = 0),
                  size = 4,
                  color = "black",
                  position = position_dodge(1),
                  hjust = 0)+
        theme(axis.title.y= element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    
    p <- p1+p2
}








entrezMapper <- function(DEGenes, geneIDs, species="mmu") {  
    require(tidyverse)

    # Kegg convention
    if (species == "mmu") {
        OrgDb <- org.Mm.eg.db
        require(org.Mm.eg.db)
    } else if (species == "hsa") {
        OrgDb <- org.Hs.eg.db
        require(org.Hs.eg.db)
    }
    keytypes(OrgDb) # UNIPROT, ALIAS etc
    
    idsFromENSEMBL <- bitr(geneIDs$ensembl.id,
                            fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL", "GENENAME"),
                            OrgDb=OrgDb, drop = F)
    
    idsFromSYMBOL <- bitr(geneIDs$gene.name,
                            fromType="SYMBOL", toType=c("ENTREZID", "ENSEMBL", "GENENAME"),
                            OrgDb=OrgDb, drop = F)

    entrezID1 <- list()
    entrezID2 <- list()
    entrezID <- list()

    mismatched <- 0
    # ENSEMBL FOR
    for (symbol_g in DEGenes$gene.name ) {
        inds <- which(idsFromENSEMBL$ENSEMBL == DEGenes$ensembl.id[DEGenes$gene.name == symbol_g])
        entrezID1[[symbol_g]] <- idsFromENSEMBL$ENTREZID[inds]
        if (length(entrezID1[[symbol_g]]) == 0){
            entrezID1[[symbol_g]] <- NA
        }
        entrezID[[symbol_g]] <- entrezID1[[symbol_g]]
    }
    
    # SYMBOL FOR
    for (symbol_g in DEGenes$gene.name ) {
        inds <- which(idsFromSYMBOL$SYMBOL == DEGenes$gene.name[DEGenes$gene.name == symbol_g])
        entrezID2[[symbol_g]] <- idsFromSYMBOL$ENTREZID[inds]
        if (length(entrezID2[[symbol_g]]) == 0){
            entrezID2[[symbol_g]] <- NA
        }
        # print(symbol_g)
        # print(entrezID1[[symbol_g]])
        # print(entrezID2[[symbol_g]])
        if(length(entrezID2[[symbol_g]]) >1 ){
            entrezID2[[symbol_g]] <- entrezID2[[symbol_g]][1]
        }
        if(length(entrezID1[[symbol_g]]) >1 ){
            entrezID1[[symbol_g]] <- entrezID1[[symbol_g]][1]
        }

        if ( !is.na(entrezID1[[symbol_g]]) && !is.na(entrezID2[[symbol_g]]) ) {
            if (entrezID1[[symbol_g]] != entrezID2[[symbol_g]]) {
                # print("WARNING: entrezID mismatch!!!!!!!!")
                mismatched <- mismatched + 1
                entrezID[[symbol_g]] <- c(entrezID2[[symbol_g]])
            }
        }


        if (all(is.na(entrezID1[[symbol_g]]))) {
            entrezID[[symbol_g]] <- entrezID2[[symbol_g]]
        }
        if(length(entrezID[[symbol_g]])==2){
            entrezID[[symbol_g]] <- entrezID[[symbol_g]][1]
        }
    }

    
    print(table(is.na(entrezID)))
    DEGenes$entrezID <- as.character(entrezID)
    return(DEGenes)
}


dotplot_enrich_dim <- function(pos_mat , neg_mat,n_paths=10){
    sub_pos_mat <- pos_mat[1:n_paths,]
    ord_sub_pos_mat <- sub_pos_mat
    ord_sub_pos_mat$term <- factor(sub_pos_mat$term,levels=rev(sub_pos_mat$term))
    ord_sub_pos_mat["P.value"] <- ord_sub_pos_mat[colnames(ord_sub_pos_mat)[grep("P.value", colnames(ord_sub_pos_mat))][1]]
    ord_sub_pos_mat["Odds.Ratio"] <- ord_sub_pos_mat[colnames(ord_sub_pos_mat)[grep("Odds.Ratio", colnames(ord_sub_pos_mat))][1]]


    sub_neg_mat <- neg_mat[1:n_paths,]
    ord_sub_neg_mat <- sub_neg_mat
    ord_sub_neg_mat$term <- factor(sub_neg_mat$term,levels=rev(sub_neg_mat$term))
    ord_sub_neg_mat["P.value"] <- ord_sub_neg_mat[colnames(ord_sub_neg_mat)[grep("P.value", colnames(ord_sub_neg_mat))][1]]
    ord_sub_neg_mat["Odds.Ratio"] <- ord_sub_neg_mat[colnames(ord_sub_neg_mat)[grep("Odds.Ratio", colnames(ord_sub_neg_mat))][1]]

    
    p <- ggplot(ord_sub_pos_mat,
            aes(    x=log10pval,
                    y=term,
                    color=P.value,
                    size=Odds.Ratio
                    ),
            title= enrich.database)+
            geom_point() +
            scale_color_continuous(low="red", high="blue",guide=guide_colorbar(reverse=TRUE))+ 
            ylab(NULL)+
            DOSE::theme_dose(12)+
    ggplot(ord_sub_neg_mat,
            aes(    x=log10pval,
                    y=term,
                    color=P.value,
                    size=Odds.Ratio
                    ),
            title= enrich.database)+
            geom_point() +
            scale_color_continuous(low="red", high="blue",guide=guide_colorbar(reverse=TRUE))+ 
            ylab(NULL)+
            DOSE::theme_dose(12)
    return(p)
}



PercentAbove_Seurat <- function(x, threshold) {
  return((length(x = x[x > threshold]) / length(x = x))*100)
}



Combined_Percent_Expressing <- function(
  seurat_object,
  features,
  threshold = 0,
  group_by = NULL,
  split_by = NULL,
  entire_object = FALSE,
  slot =NULL,
  assay = NULL
) {
    # Check Seurat
    Is_Seurat(seurat_object = seurat_object)

    # set assay (if null set to active assay)
    assay <- assay %||% DefaultAssay(object = seurat_object)

    # Check features exist in object
    features_list <- Gene_Present(data = seurat_object, gene_list = features, print_msg = FALSE, case_check = TRUE, seurat_assay = assay)[[1]]

    # Check group_by is in object
    if (!is.null(x = group_by)) {
        possible_groups <- colnames(seurat_object@meta.data)
        if (!group_by %in% possible_groups) {
            cli_abort("Grouping variable {.val {group_by}} was not found in Seurat Object.")
        }
    }

    # Check split_by is in object
    if (!is.null(x = split_by)) {
        possible_groups <- colnames(seurat_object@meta.data)
        if (!split_by %in% possible_groups) {
            cli_abort("Splitting variable {.val {split_by}} was not found in Seurat Object.")
        }
    }

    # Pull Expression Info
    title_genes <- paste0(paste(features,collapse="+"),"+")
    cells <- unlist(x = CellsByIdentities(object = seurat_object, idents = NULL))
    expression_info_org <- FetchData(object = seurat_object, vars = features_list, cells = cells, slot = slot)
    expression_info <- expression_info_org
    expression_info[expression_info>=1] <- 1
    expression_info[title_genes] <- ifelse(rowSums(expression_info)==length(features),1,0)
    expression_info[paste0("Joint_",title_genes)] <-rowSums(expression_info_org)

    features_list <- colnames(expression_info)
    # Add grouping variable
    if (entire_object) {
        expression_info$id <- "All_Cells"
    } else {
        expression_info$id <- if (is.null(x = group_by)) {
            Idents(object = seurat_object)[cells, drop = TRUE]
        } else {
            seurat_object[[group_by, drop = TRUE]][cells, drop = TRUE]
        }
    }
    if (!is.factor(x = expression_info$id)) {
        expression_info$id <- factor(x = expression_info$id)
    }
    id.levels <- levels(x = expression_info$id)
    expression_info$id <- as.vector(x = expression_info$id)

    # Split data if split.by is true
    if (!is.null(x = split_by)) {
        splits <- seurat_object[[split_by, drop = TRUE]][cells, drop = TRUE]
        expression_info$id <- paste(expression_info$id, splits, sep = '_')
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }

    # Calculate percent expressing
    percent_expressing <- lapply(
        X = unique(x = expression_info$id),
        FUN = function(ident) {
            data.use <- expression_info[expression_info$id == ident, 1:(ncol(x = expression_info) - 1), drop = FALSE]
            pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove_Seurat, threshold = threshold)
            return(list(pct.exp = pct.exp))
        }
    )
    names(x = percent_expressing) <- unique(x = expression_info$id)

    # Convert & return data.frame
    row_dim_names <- features_list
    col_dim_names <- names(percent_expressing)
    mat_dims <- list(row_dim_names, col_dim_names)
    final_df <- data.frame(matrix(unlist(percent_expressing), nrow = length(features_list), byrow = FALSE, dimnames = mat_dims), stringsAsFactors = FALSE)
    expression_info <- expression_info[colnames(seurat_object),]

    return(list(final_df=final_df,expression_info=expression_info))
}


Is_Seurat <- function(
  seurat_object
) {
  if (!inherits(what = "Seurat", x = seurat_object)) {
    cli_abort(message = "{.code seurat_object} provided is not an object of class: Seurat.")
  }
}


`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}


`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}



Overall_dens_Count_Joint_plot <- function(Seurat_CCR7pos,Seurat_CCR7neg,features,resultdir){
    title_genes_1 <- paste0(paste(features,collapse="+"),"+")
    title_genes_2 <- paste(features,collapse="_")
    
    # ================================== Positive CCR7 =================================
    data_df_res_all <- Combined_Percent_Expressing(
        seurat_object=Seurat_CCR7pos,
        features= features,
        threshold = 0,
        group_by = "CellType",
        split_by = "Condition",
        entire_object = FALSE,
        slot = "counts",
        assay = "RNA"
    )
    data_df_res <- data_df_res_all$final_df
    expression_info <- data_df_res_all$expression_info[colnames(Seurat_CCR7pos),]
    expression_info_vec <- as.vector(expression_info[[title_genes_1]])
    names(expression_info_vec) <- rownames(x = expression_info)
    Seurat_CCR7pos <- AddMetaData(
        object = Seurat_CCR7pos,
        metadata = expression_info_vec,
        col.name = title_genes_2
    )

    # ================================== NEGATTIVE CCR7 =================================
    data_df_res_all <- Combined_Percent_Expressing(
        seurat_object=Seurat_CCR7neg,
        features= features,
        threshold = 0,
        group_by = "CellType",
        split_by = "Condition",
        entire_object = FALSE,
        slot = "counts",
        assay = "RNA"
    )
    data_df_res <- data_df_res_all$final_df
    expression_info <- data_df_res_all$expression_info[colnames(Seurat_CCR7neg),]
    expression_info_vec <- as.vector(expression_info[[title_genes_1]])
    names(expression_info_vec) <- rownames(x = expression_info)
    Seurat_CCR7neg <- AddMetaData(
        object = Seurat_CCR7neg,
        metadata = expression_info_vec,
        col.name = title_genes_2
    )



    UP_VECTOR <-   FetchData(object = Seurat_CCR7pos, vars = title_genes_2)
    down_VECTOR <-   FetchData(object = Seurat_CCR7neg, vars = title_genes_2)
    if(sum(UP_VECTOR)==0){
        p2.1 <- FeaturePlot(Seurat_CCR7pos,  features=title_genes_2,pt.size=0.2,order=T,split.by = "Combined", cols = c("gray","gray")) +patchwork::plot_layout(ncol = 4, nrow = 1)
    }else{
        p2.1 <- FeaturePlot(Seurat_CCR7pos,  features=title_genes_2,pt.size=0.2,order=T,split.by = "Combined") +patchwork::plot_layout(ncol = 4, nrow = 1)
    }
    if(sum(down_VECTOR)==0){
        p2.2 <- FeaturePlot(Seurat_CCR7neg,  features=title_genes_2,pt.size=0.2,order=T,split.by = "Combined", cols = c("gray","gray")) +patchwork::plot_layout(ncol = 4, nrow = 1)
    }else{
        p2.2 <- FeaturePlot(Seurat_CCR7neg,  features=title_genes_2,pt.size=0.2,order=T,split.by = "Combined") +patchwork::plot_layout(ncol = 4, nrow = 1)
    }

    p2 <- p2.1/p2.2 
    # &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Blues")),limits = range(c(0, 1)))
    #)scale_colour_continuous(limits = range(c(0, 1)))
    ggsave(plot=p2,filename=paste0(resultdir,"Split_Counts_",title_genes_2,".pdf"),width=12,height=8)
    return(p2)
}







myDotPlot <- function(
  object,
  assay = NULL,
  features,
  cols = c("grey90", "#C51B7D"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = 'radius',
  scale.min = NA,
  scale.max = NA
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))

  data.features <- FetchData(object = object, vars = features, cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(
      what = rbind,
      args = lapply(X = data.plot, FUN = unlist)
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  } else if (ngroup < 5 & scale) {
    warning(
      "Scaling data with a low number of groups may produce misleading results",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'pct.exp']
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log1p(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(
      X = as.character(x = data.plot$id),
      FUN = gsub,
      FUN.VALUE = character(length = 1L),
      pattern =  paste0(
        '^((',
        paste(sort(x = levels(x = object), decreasing = TRUE), collapse = '|'),
        ')_)'
      ),
      replacement = '',
      USE.NAMES = FALSE
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )
  }
  color.by <- ifelse(test = split.colors, yes = 'colors', no = 'pct.exp')
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
      x = 'Features',
      y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
    ) +
    theme_cowplot()
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(
      facets = ~feature.groups,
      scales = "free_x",
      space = "free_x",
      switch = "y"
    ) + theme(
      panel.spacing = unit(x = 1, units = "lines"),
      strip.background = element_blank()
    )
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  } else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = 'Percent Expressed'))
  }
  return(plot)
}