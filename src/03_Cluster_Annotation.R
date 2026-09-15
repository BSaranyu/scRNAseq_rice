library(Seurat)
library(dplyr)
library(tidyr)
library(data.table)

DefaultAssay(mergedDat) <- "RNA"
mergedDat <- NormalizeData(mergedDat, normalization.method = "LogNormalize",scale.factor = 1e6)
mergedDat[["RNA"]] <- JoinLayers(mergedDat[["RNA"]])

marker <- as.data.frame(fread('data/ReportedMarker.csv',header = T))
rownames(marker) <- 1:dim(marker)[1]
colnames(marker)[5] <- "gene" 

marker1 <- marker[marker$Organ %in% c("Leaf","Stem"),]
marker2 <- marker1 %>% separate_rows(Tissue,sep = ",")
marker2 <- marker2[marker2$Tissue != "leaf", ]
marker2 <- marker2[marker2$Tissue != "flag leaf", ]

rna_data <- GetAssayData(mergedDat, assay = "RNA", layer = "data")

hvg_genes <- VariableFeatures(mergedDat)
hvg_markers <- intersect(hvg_genes, unique(marker1$gene))
extra_genes <- marker1$gene[marker1$Tissue == "mesocotyl"]
final_genes <- unique(c(hvg_markers, extra_genes))
final_genes <- intersect(final_genes, rownames(rna_data))
marker1_filtered <- marker1[marker1$gene %in% final_genes, ]
expr_data <- as.data.frame(rna_data[final_genes, ])

cluster_ids <- mergedDat$clusters
n_clusters <- length(unique(Idents(mergedDat)))
weight_data <- matrix(NA, nrow = length(final_genes), ncol = n_clusters)
marker1_hvg <- marker1[marker1$gene %in% hvg_genes, ]

rownames(weight_data) <- final_genes
colnames(weight_data) <- paste0("Cluster ", sort(unique(Idents(mergedDat))))

for (i in 0:20) {
  cells_in_cluster <- names(cluster_ids[cluster_ids == i])
  
  for (j in seq_along(final_genes)) {
    gene <- final_genes[j]
    if (gene %in% rownames(expr_data)) {
      expr_vals <- as.numeric(expr_data[gene, cells_in_cluster])
      expr_vals <- expm1(expr_vals)
      avg_expr <- mean(expr_vals, na.rm = TRUE)
      weight_data[j, i + 1] <- log2(avg_expr + 1)
    } else {
      weight_data[j, i + 1] <- NA
    }
    
  }
}

gene_variance <- apply(weight_data, 1, var, na.rm = TRUE)
weight_data_with_var <- cbind(weight_data, Variance = gene_variance)
weight_data_with_var <- cbind(weight_data_with_var, gene = rownames(weight_data_with_var))
weight_data <- merge(weight_data_with_var, marker1_filtered[, c("gene", "Organ", "Tissue")], by = "gene")
weight_data <- weight_data %>% filter(!if_all(2:(n_clusters + 2), is.na))

expr_scale_data <- as.data.frame(t(scale(t(expr_data))))
expr_scale_data <- merge(marker2,expr_scale_data,by.x='gene',by.y="row.names")
expr_scale_data <- merge(weight_data[,c("gene","Variance")],expr_scale_data,by="gene")
expr_scale_data$Variance <- as.numeric(expr_scale_data$Variance)
expr_scale_data[,colnames(rna_data)] <- expr_scale_data[,colnames(rna_data)]*expr_scale_data$Variance

na_rows <- apply(expr_scale_data[, colnames(rna_data)], 1, function(x) all(is.na(x))) 
expr_scale_data <- expr_scale_data[!na_rows, ] 

fun <- function(x) {
  tissue_name <- unique(expr_scale_data$Tissue)[x]
  tmp <- expr_scale_data[expr_scale_data$Tissue == tissue_name, colnames(rna_data)]
  
  tmp1 <- as.data.frame(t(as.data.frame(colSums(tmp, na.rm = TRUE))))
  rownames(tmp1) <- tissue_name
  return(tmp1)
}

MICI_out <- do.call('rbind', parallel::mclapply(
  1:length(unique(expr_scale_data$Tissue)),
  function(x) { fun(x) },
  mc.cores = 1
))

MICI_out$cell_type <- as.character(rownames(MICI_out))

cell_type_col <- MICI_out$cell_type
expr_mat <- MICI_out[, -which(colnames(MICI_out) == "cell_type")]
MICI_result <- sapply(expr_mat, function(x) {
  idx <- which.max(x)  
  cell_type_col[idx]
})

MICI_result <- data.frame(cell_type = MICI_result)
colnames(MICI_result) <- c("identified_ct")
umap_coords <- Embeddings(mergedDat, reduction = "umap")

plt <- as.data.frame(umap_coords)
plt$cell <- rownames(umap_coords)
plt$conditions <- mergedDat$conditions
plt$cluster <- mergedDat$clusters

MICI_result <- merge(plt,MICI_result,by.x="cell",by.y="row.names")
MICI_result$identified_ct <- as.character(MICI_result$identified_ct)
tmp <- as.data.frame(t(MICI_out[,-dim(MICI_out)[2]]))
max_values <- apply(tmp, 1, max, na.rm = TRUE) 
unct <- rownames(tmp[max_values < 1,]) # adjust the threshold to 1 (default = 2)
MICI_result[MICI_result$cell %in% unct,]$identified_ct <- "Unknown" 

dat1 <- as.data.frame(table(MICI_result$identified_ct,MICI_result$cluster))
colnames(dat1) <- c("cell_type","cluster","number_of_cell")
dat2 <- as.data.frame(table(MICI_result$cluster))
colnames(dat2) <- c("cluster","number_of_cluster")
dat3 <- merge(dat1,dat2,by="cluster")
dat3$ratio <- signif(dat3$number_of_cell/dat3$number_of_cluster,3)
dat4 <- dat3 %>% group_by(cluster) %>% top_n(1,ratio) %>% as.data.frame()
dat5 <- reshape2::dcast(dat3,cluster~cell_type,value.var = "ratio")
cluster_assign <- merge(dat4[,c("cluster","cell_type")],dat5,by="cluster")

cellType <- setNames(cluster_assign$cell_type, levels(mergedDat))
mergedDat <- RenameIdents(mergedDat, cellType)
mergedDat$celltype <- Idents(mergedDat)
