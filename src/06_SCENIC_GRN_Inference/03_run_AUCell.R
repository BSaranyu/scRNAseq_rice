library(dplyr)
library(Seurat)
library(AUCell)
library(biomaRt)
library(AnnotationHub)
library(clusterProfiler)

normal_df <- read.csv("normal_mesophyll_regulon_Q90.csv", stringsAsFactors = FALSE)
drought_df <- read.csv("drought_mesophyll_regulon_Q90.csv", stringsAsFactors = FALSE)
normal_epidermis_df <- read.csv("normal_epidermis_regulon_Q90.csv", stringsAsFactors = FALSE)
drought_epidermis_df <- read.csv("drought_epidermis_regulon_Q90.csv", stringsAsFactors = FALSE)
make_regulon_list <- function(df) {
  df %>%
    filter(grepl("directAnnotation", TF_highConf)) %>%
    dplyr::select(motif, targetGenesList) %>%
    group_by(motif) %>%
    summarise(all_targets = paste(targetGenesList, collapse = ";"), .groups = "drop") %>%
    mutate(all_targets = strsplit(all_targets, ";")) %>%
    deframe()
}
regulons_normal <- make_regulon_list(normal_df)
regulons_drought <- make_regulon_list(drought_df)
regulons_normal_epidermis <- make_regulon_list(normal_epidermis_df)
regulons_drought_epidermis <- make_regulon_list(drought_epidermis_df)

mesophyll <- readRDS("leaf_mesophyll.rds")
epidermis <- readRDS("leaf_epidermis.rds")
n_mesophyll <- subset(mesophyll, subset = conditions == "RN1")
d_mesophyll <- subset(mesophyll, subset = conditions == "RD2")
n_epidermis <- subset(epidermis, subset = conditions == "RN1")
d_epidermis <- subset(epidermis, subset = conditions == "RD2")
n_exprMatrix <- GetAssayData(n_mesophyll, assay = "RNA", layer = "data")
d_exprMatrix <- GetAssayData(d_mesophyll, assay = "RNA", layer = "data")
n_exprMatrix_epidermis <- GetAssayData(n_epidermis, assay = "RNA", layer = "data")
d_exprMatrix_epidermis <- GetAssayData(d_epidermis, assay = "RNA", layer = "data")

rankings_normal <- AUCell_buildRankings(n_exprMatrix, plotStats = FALSE)
auc_normal <- AUCell_calcAUC(regulons_normal, rankings_normal)
auc_normal_mat <- getAUC(auc_normal)
rankings_drought <- AUCell_buildRankings(d_exprMatrix, plotStats = FALSE)
auc_drought <- AUCell_calcAUC(regulons_drought, rankings_drought)
auc_drought_mat <- getAUC(auc_drought)
rankings_normal_epidermis <- AUCell_buildRankings(n_exprMatrix_epidermis, plotStats = FALSE)
auc_normal_epidermis <- AUCell_calcAUC(regulons_normal_epidermis, rankings_normal_epidermis)
auc_normal_epidermis_mat <- getAUC(auc_normal_epidermis)
rankings_drought_epidermis <- AUCell_buildRankings(d_exprMatrix_epidermis, plotStats = FALSE)
auc_drought_epidermis <- AUCell_calcAUC(regulons_drought_epidermis, rankings_drought_epidermis)
auc_drought_epidermis_mat <- getAUC(auc_drought_epidermis)

# GO Analysis
ensembl_oryzaSativaJ <- useMart(biomart = "plants_mart", dataset = "osativa_eg_gene", host = "https://plants.ensembl.org")
gene_symb <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl_oryzaSativaJ)
names(gene_symb)[names(gene_symb) == "ensembl_gene_id"] <- "RAP"
gene_symb$external_gene_name[is.na(gene_symb$external_gene_name) | gene_symb$external_gene_name == ""] <- gene_symb$RAP[is.na(gene_symb$external_gene_name) | gene_symb$external_gene_name == ""]
gene_to_RAP <- setNames(gene_symb$RAP, gene_symb$external_gene_name)
IDtable <- read.csv("riceIDtable.csv", stringsAsFactors = FALSE)
rap_to_entrez <- setNames(as.character(IDtable$entrezgene), IDtable$rapdb)
hub <- AnnotationHub()
rice <- hub[["AH114586"]]
convert_regulon_to_entrez <- function(tf_name, regulon_list, gene_to_RAP, rap_to_entrez, include_tf = TRUE) {
  if (!tf_name %in% names(regulon_list)) return(character(0))
  genes <- regulon_list[[tf_name]]
  if (include_tf) {
    genes <- c(tf_name, genes)
  }
  rap_ids <- gene_to_RAP[genes]
  rap_ids <- unname(rap_ids)
  rap_ids <- rap_ids[!is.na(rap_ids) & rap_ids != ""]
  entrez_ids <- rap_to_entrez[rap_ids]
  entrez_ids <- unname(entrez_ids)
  entrez_ids <- entrez_ids[!is.na(entrez_ids) & entrez_ids != ""]
  unique(entrez_ids)
}
run_go_for_regulon <- function(tf_name, regulon_list, gene_to_RAP, rap_to_entrez, OrgDb,
                               ont = "BP", p_cutoff = 0.05, q_cutoff = 0.2,
                               minGSSize = 10, maxGSSize = 500) {
  entrez_ids <- convert_regulon_to_entrez(
    tf_name = tf_name,
    regulon_list = regulon_list,
    gene_to_RAP = gene_to_RAP,
    rap_to_entrez = rap_to_entrez,
    include_tf = TRUE
  )
  if (length(entrez_ids) < minGSSize) return(NULL) 
  go_res <- tryCatch({
    enrichGO(
      gene = entrez_ids,
      OrgDb = OrgDb,
      ont = ont,
      pvalueCutoff = p_cutoff,
      pAdjustMethod = "BH",
      qvalueCutoff = q_cutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      readable = FALSE
    )
  }, error = function(e) NULL) 
  if (is.null(go_res) || nrow(as.data.frame(go_res)) == 0) return(NULL) 
  go_df <- as.data.frame(go_res) %>%
    filter(!is.na(Description), !is.na(p.adjust), !is.na(Count)) %>%
    mutate(TF = tf_name)
  go_df
}
run_go_all_regulons <- function(regulon_list, celltype, condition,
                                gene_to_RAP, rap_to_entrez, OrgDb,
                                ont = "BP", p_cutoff = 0.05, q_cutoff = 0.2,
                                minGSSize = 10, maxGSSize = 500) { 
  tf_names <- names(regulon_list)
  results <- purrr::map(tf_names, function(tf) {
    res <- run_go_for_regulon(
      tf_name = tf,
      regulon_list = regulon_list,
      gene_to_RAP = gene_to_RAP,
      rap_to_entrez = rap_to_entrez,
      OrgDb = OrgDb,
      ont = ont,
      p_cutoff = p_cutoff,
      q_cutoff = q_cutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    ) 
    if (!is.null(res)) {
      res$Celltype <- celltype
      res$Condition <- condition
    }
    res
  })
  bind_rows(results)
}
# Normal mesophyll
go_mesophyll_normal <- run_go_all_regulons(
  regulon_list = regulons_normal,
  celltype = "Mesophyll",
  condition = "Normal",
  gene_to_RAP = gene_to_RAP,
  rap_to_entrez = rap_to_entrez,
  OrgDb = rice
)
# Drought mesophyll
go_mesophyll_drought <- run_go_all_regulons(
  regulon_list = regulons_drought,
  celltype = "Mesophyll",
  condition = "Drought",
  gene_to_RAP = gene_to_RAP,
  rap_to_entrez = rap_to_entrez,
  OrgDb = rice
)
# Normal epidermis
go_epidermis_normal <- run_go_all_regulons(
  regulon_list = regulons_normal_epidermis,
  celltype = "Epidermis",
  condition = "Normal",
  gene_to_RAP = gene_to_RAP,
  rap_to_entrez = rap_to_entrez,
  OrgDb = rice
)
# Drought epidermis
go_epidermis_drought <- run_go_all_regulons(
  regulon_list = regulons_drought_epidermis,
  celltype = "Epidermis",
  condition = "Drought",
  gene_to_RAP = gene_to_RAP,
  rap_to_entrez = rap_to_entrez,
  OrgDb = rice
)
go_all_results <- bind_rows(go_mesophyll_normal, go_mesophyll_drought, go_epidermis_normal, go_epidermis_drought)
