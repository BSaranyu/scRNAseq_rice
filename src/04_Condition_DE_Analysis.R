library(Seurat)
library(dplyr)
library(clusterProfiler)
library(biomaRt)
library(AnnotationHub)

mergedDat <- readRDS(file = "mergedDat.rds")
find_DEGs_between_conditions <- function(seurat_obj, cell_type, condt1, condt2, 
                                         logfc_threshold = 0, min_pct = 0.1, 
                                         p_adj_threshold = 10000, test_method = "wilcox") {
  
  subset_obj <- subset(seurat_obj, idents = cell_type)
  Idents(subset_obj) <- subset_obj@meta.data$conditions
  
  degs <- FindMarkers(
    subset_obj,
    ident.1 = condt1,
    ident.2 = condt2,
    logfc.threshold = logfc_threshold,
    min.pct = min_pct,
    test.use = test_method
  )
  degs_sig <- degs[degs$p_val_adj < p_adj_threshold, ]
  degs_sig$gene <- rownames(degs_sig)
  
  return(degs_sig)
}
cell_types <- levels(Idents(mergedDat))
DEGs_list_byCellType <- list()
for (cell in cell_types) {
  DEGs_list_byCellType[[cell]] <- find_DEGs_between_conditions(
    seurat_obj = mergedDat,
    cell_type = cell,
    condt1 = "RD2",
    condt2 = "RN1"
  )
}
add_significance <- function(degs_df, 
                             logfc_threshold = 1, 
                             p_adj_threshold = 0.05,
                             up_label = "Up in Drought",
                             down_label = "Down in Drought",
                             ns_label = "Not significant") {
  
  degs_df$significance <- ns_label
  degs_df$significance[
    degs_df$avg_log2FC > logfc_threshold & 
      degs_df$p_val_adj < p_adj_threshold
  ] <- up_label
  degs_df$significance[
    degs_df$avg_log2FC < -logfc_threshold & 
      degs_df$p_val_adj < p_adj_threshold
  ] <- down_label
  
  return(degs_df)
}
DEGs_list_byCellType <- lapply(
  DEGs_list_byCellType,
  add_significance
)
ms_DEGs    <- DEGs_list_byCellType[["leaf mesophyll"]]
mspre_DEGs <- DEGs_list_byCellType[["Mesophyll precursor"]]
ep_DEGs    <- DEGs_list_byCellType[["leaf epidermis"]]
mst_DEGs   <- DEGs_list_byCellType[["mestome sheath"]]
init_DEGs  <- DEGs_list_byCellType[["leaf initial cell"]]
phm_DEGs   <- DEGs_list_byCellType[["Phloem"]]
ukn_DEGs   <- DEGs_list_byCellType[["Unkown"]]

# GO Enrichment Analysis
hub <- AnnotationHub()
rice <- hub[["AH114586"]]
ensembl_oryzaSativaJ <- useMart(biomart = "plants_mart", dataset = "osativa_eg_gene", host = "https://plants.ensembl.org")
gene_symb <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   mart = ensembl_oryzaSativaJ)
colnames(gene_symb)[colnames(gene_symb) == "ensembl_gene_id"] <- "RAP"
colnames(gene_symb)[colnames(gene_symb) == "external_gene_name"] <- "Gene.ID"
IDtable <- read.csv("data/riceIDtable.csv")

convert_to_RAP <- function(genes, gene_symb) {  
  is_RAP <- stringr::str_detect(genes, "^Os\\d{2}g\\d{7}$")  
  genes_df <- data.frame(gene_input = genes, stringsAsFactors = FALSE) %>%
    dplyr::mutate(RAP_mine = ifelse(is_RAP, gene_input, NA)) %>%
    dplyr::left_join(gene_symb, by = c("gene_input" = "Gene.ID")) %>%
    dplyr::mutate(RAP_final = dplyr::coalesce(RAP_mine, RAP)) %>%
    dplyr::select(RAP_final) %>%
    dplyr::distinct() %>%
    na.omit()  
  return(genes_df$RAP_final)
}

run_GO_enrichment <- function(rap_genes, IDtable, OrgDb, ontology = "BP") {
  genes_eid <- IDtable[match(rap_genes, IDtable$rapdb), "entrezgene"]
  genes_eid <- as.character(genes_eid[!is.na(genes_eid)])  
  go_result <- enrichGO(
    gene = genes_eid,
    OrgDb = OrgDb,
    ont = ontology,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500
  )  
  return(as.data.frame(go_result))
}

process_celltype_GO <- function(DEG_df,
                                gene_symb,
                                IDtable,
                                OrgDb,
                                celltype_name = NULL) {
 
  DEG_up <- DEG_df %>% dplyr::filter(significance == "Up in Drought")
  DEG_down <- DEG_df %>% dplyr::filter(significance == "Down in Drought")
  
  rap_up <- convert_to_RAP(DEG_up$gene, gene_symb)
  rap_down <- convert_to_RAP(DEG_down$gene, gene_symb)
  
  go_up <- run_GO_enrichment(rap_up, IDtable, OrgDb)
  go_down <- run_GO_enrichment(rap_down, IDtable, OrgDb)
  
  return(list(
    celltype = celltype_name,
    up_GO = go_up,
    down_GO = go_down,
    up_genes = rap_up,
    down_genes = rap_down
  ))
}
DEG_list <- list(
  Mesophyll = ms_DEGs,
  Mesophyll_precursor = mspre_DEGs,
  Epidermis = ep_DEGs,
  Mestome = mst_DEGs,
  Initial_cell = init_DEGs,
  Phloem = phm_DEGs,
  Unknown = ukn_DEGs
)

GO_results <- lapply(names(DEG_list), function(ct) { 
  process_celltype_GO(
    DEG_df = DEG_list[[ct]],
    gene_symb = gene_symb,
    IDtable = IDtable,
    OrgDb = rice,
    celltype_name = ct
  )
})
names(GO_results) <- names(DEG_list)

Mesophyll_go_up <- GO_results$Mesophyll$up_GO
Mesophyll_precursor_go_up <- GO_results$Mesophyll_precursor$up_GO
Epidermis_go_up <- GO_results$Epidermis$up_GO
Mestome_go_up <- GO_results$Mestome$up_GO
Initial_cell_go_up <- GO_results$Initial_cell$up_GO
Phloem_go_up <- GO_results$Phloem$up_GO
Unknown_go_up <- GO_results$Unknown$up_GO

Mesophyll_go_down <- GO_results$Mesophyll$down_GO
Mesophyll_precursor_go_down <- GO_results$Mesophyll_precursor$down_GO
Epidermis_go_down <- GO_results$Epidermis$down_GO
Mestome_go_down <- GO_results$Mestome$down_GO
Initial_cell_go_down <- GO_results$Initial_cell$down_GO
Phloem_go_down <- GO_results$Phloem$down_GO
Unknown_go_down <- GO_results$Unknown$down_GO
