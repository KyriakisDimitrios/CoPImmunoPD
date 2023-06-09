

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




    return(final_df)
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