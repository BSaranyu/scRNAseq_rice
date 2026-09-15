library(dplyr)
library(readr)
library(data.table)
library(glue)
library(org.Osativa.eg.db)
library(biomaRt)
library(RcisTarget)

# Specify cell type and condition
cell_type <- "mesophyll"
condition <- "normal"

grn <- read_tsv(glue("output/{cell_type}/{cell_type}_{condition}_GRNBoost2.tsv"))
quantile_cutoff <- quantile(grn$importance, 0.9)
grn_filtered <- grn[grn$importance >= quantile_cutoff, ]
geneLists <- grn_filtered %>%
  group_by(TF) %>%
  summarise(targets = list(unique(target))) %>%
  deframe()
motifRankings <- importRankings(glue("data/Ranking/{cell_type}/{condition}_{cell_type}_cisTarget_rankings.feather"), indexCol = "motif")
motifAnnotation <- fread("data/motifAnnotations_rice.tsv")
genRap <- motifAnnotation$TF_RAP
rap_to_msu <- AnnotationDbi::select(
  org.Osativa.eg.db,
  keys = genRap,
  columns = c("GID"),   # MSU ID column
  keytype = "RAP"
)
ensembl_oryzaSativaJ <- useEnsemblGenomes(biomart = "plants_mart", dataset = "osativa_eg_gene")
gene_symb <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl_oryzaSativaJ)
names(gene_symb)[names(gene_symb) == "ensembl_gene_id"] <- "RAP"
combined_df <- merge(rap_to_msu, gene_symb, by = "RAP", all.x = TRUE)
combined_df$GID[is.na(combined_df$GID) | combined_df$GID == ""] <- combined_df$RAP[is.na(combined_df$GID) | combined_df$GID == ""]
combined_df$GID[combined_df$GID == ""] <- combined_df$RAP[combined_df$GID == ""]
combined_df$external_gene_name[is.na(combined_df$external_gene_name) | combined_df$external_gene_name == ""] <- combined_df$RAP[is.na(combined_df$external_gene_name) | combined_df$external_gene_name == ""]
gene_map <- setNames(combined_df$external_gene_name, combined_df$RAP)
converted_genes <- gene_map[genRap]
motifAnnotation$TF_id <- make.unique(ifelse(!is.na(converted_genes), converted_genes, genRap))
setnames(motifAnnotation, c("motif_id", "TF","TF_id"), c("motif", "TF_MSU","TF"))
motifAnnotation[, directAnnotation := TF_highConf %in% c("High", "Medium","Low")]
motifAnnotation[, inferred_Orthology := FALSE]
motifAnnotation[, inferred_MotifSimil := FALSE]
motifAnnotation[, annotationSource := ifelse(
  directAnnotation == TRUE,
  "directAnnotation",
  "inferredBy_MotifSimilarity"
)]
motifAnnotation[, annotationSource := factor(annotationSource, levels = c(
  "directAnnotation", 
  "inferredBy_Orthology", 
  "inferredBy_MotifSimilarity", 
  "inferredBy_MotifSimilarity_n_Orthology"
))]
setkeyv(motifAnnotation, c("motif", "TF"))
motifs_AUC <- calcAUC(geneLists, motifRankings)
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                                           motifAnnot = motifAnnotation,
                                           motifAnnot_highConfCat = c("directAnnotation"),
                                           motifAnnot_lowConfCat = character(0),
                                           idColumn = "TF")
motifEnrichmentTable$TF_match <- motifEnrichmentTable$geneSet == motifEnrichmentTable$TF
colnames(motifEnrichmentTable)[2] <- "motif"
regulonTargetsInfo <- addSignificantGenes(
  resultsTable = motifEnrichmentTable,
  rankings = motifRankings,
  geneSets = geneLists,
  method = "aprox"
)
strongLinks <- regulonTargetsInfo %>%
  filter(NES > 3, nEnrGenes >= 50)
strongLinks$targetGenesList <- strsplit(strongLinks$enrichedGenes, ";")
tfTargets <- setNames(strongLinks$targetGenesList, strongLinks$TF_highConf)
strongLinks$targetGenesList <- sapply(strongLinks$targetGenesList, function(x) paste(x, collapse = ";"))
