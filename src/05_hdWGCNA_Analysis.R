library(Seurat)
library(hdWGCNA)
library(dplyr)
library(patchwork)
library(biomaRt)
library(stringr)
library(AnnotationHub)
library(clusterProfiler)

# Co-expression network constrction
mergedDat <- readRDS(file = "mergedDat.rds")
mergedDat <- SetupForWGCNA(
  mergedDat,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "coexpression network of rice leaves" 
)
mergedDat <- MetacellsByGroups(
  seurat_obj = mergedDat,
  group.by = c("celltype", "conditions"),
  reduction = 'harmony',
  k = 25,
  max_shared = 10,
  ident.group = 'celltype',
  assay = "RNA"
)
mergedDat <- NormalizeMetacells(mergedDat)
mergedDat <- ScaleMetacells(mergedDat, features=VariableFeatures(mergedDat))
mergedDat <- RunPCAMetacells(mergedDat, features=VariableFeatures(mergedDat))
DefaultAssay(mergedDat) <- "RNA"
mergedDat <- SetDatExpr(
  mergedDat,
  group_name = c(mergedDat$celltype %>% unique()),
  group.by='celltype',
  assay = 'RNA',
  layer = 'data' 
)
mergedDat <- TestSoftPowers(
  mergedDat,
  networkType = 'signed'
)
plot_list <- PlotSoftPowers(mergedDat)
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(mergedDat)
mergedDat <- ConstructNetwork(mergedDat, tom_name = 'rice leaves')
TOM <- GetTOM(mergedDat)
mergedDat <- ScaleData(mergedDat, features=VariableFeatures(mergedDat))
mergedDat$celltype_condition <- paste(mergedDat$celltype, mergedDat$conditions, sep="_")
hMEs <- GetMEs(mergedDat) # harmonized module eigengenes
MEs <- GetMEs(mergedDat, harmonized=FALSE) # module eigengenes
mergedDat <- ModuleConnectivity(mergedDat, group.by = 'celltype', group_name = c(mergedDat$celltype %>% unique()))
modules <- GetModules(mergedDat) %>% subset(module != 'grey')
head(modules[,1:6])
hub_genes <- GetHubGenes(mergedDat, n_hubs = 10)
mergedDat <- ModuleExprScore(
  mergedDat,
  n_genes = 25,
  method='UCell'
)
mergedDat <- RunHarmonyMetacells(mergedDat, group.by.vars='conditions')
mergedDat <- RunUMAPMetacells(mergedDat, reduction='harmony', dims=1:15)

plot_list <- ModuleFeaturePlot(
  mergedDat,
  features='hMEs', 
  order=TRUE
)
wrap_plots(plot_list, ncol=6)

# GO Enrichment Analysis
ensembl_oryzaSativaJ <- useMart(biomart = "plants_mart", dataset = "osativa_eg_gene", host = "https://plants.ensembl.org")
gene_symb <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl_oryzaSativaJ)
colnames(gene_symb) <- c("RAP", "Gene.ID")
convert_to_RAP <- function(genes) {
  is_RAP <- stringr::str_detect(genes, "^Os\\d{2}g\\d{7}$")
  genes_df <- data.frame(gene_input = genes, stringsAsFactors = FALSE) %>%
    mutate(RAP_mine = ifelse(is_RAP, gene_input, NA)) %>%
    left_join(gene_symb, by = c("gene_input" = "Gene.ID")) %>%
    mutate(RAP_final = coalesce(RAP_mine, RAP)) %>%
    dplyr::select(gene_input, RAP_final) %>%
    distinct()
  na.omit(genes_df)
}
IDtable <- read.csv("riceIDtable.csv")

hub <- AnnotationHub()
rice <- hub[["AH114586"]]

modules <- GetModules(mergedDat)
module_list <- c("black", "yellow", "pink", "cyan", "tan", "midnightblue", "blue",
                 "green","turquoise","red","greenyellow","magenta","brown","salmon",
                 "purple","lightcyan")
GO_results_modules <- list()
for (mod in module_list) {
  
  module_genes <- modules %>%
    filter(module == mod) %>%
    pull(gene_name)
  
  rap_df <- convert_to_RAP(module_genes)
  gene_rap <- rap_df$RAP_final
  
  genes_eid <- IDtable[match(gene_rap, IDtable$rapdb), "entrezgene"]
  genes_eid <- as.character(genes_eid[!is.na(genes_eid)])
  
  if (length(genes_eid) < 10) next
  
  go_res <- enrichGO(
    gene = genes_eid,
    OrgDb = rice,
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500
  )
  
  go_df <- as.data.frame(go_res) %>%
    filter(p.adjust < 0.05) %>%
    mutate(Module = mod,
           logP = -log10(p.adjust))
  
  GO_results_modules[[mod]] <- go_df
}
