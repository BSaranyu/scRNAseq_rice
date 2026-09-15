library(Seurat)
library(SeuratObject)
library(dplyr)
library(patchwork)
library(DoubletFinder)
library(ggplot2)
library(data.table)
library(tidyverse)
library(readxl)

rn1.cts <- Read10X_h5(filename = "filtered_feature_bc_matrix_rn1.h5")

r_mt <- rownames(rn1.cts)
mt_table <- read_excel("data/mt_cp_geneList.xlsx", sheet = "mt") %>% as.data.frame()
mt_genes <- mt_table$external_gene_name
mt_matched <- r_mt %in% mt_genes
r_mt[mt_matched] <- gsub("^gene-","MT-",r_mt[mt_matched])
rownames(rn1.cts) <- r_mt

r_cp <- rownames(rn1.cts)
cp_table <- read_excel("data/mt_cp_geneList.xlsx", sheet = "cp") %>% as.data.frame()
cp_genes <- cp_table$external_gene_name
cp_matched <- r_cp %in% cp_genes
r_cp[cp_matched] <- gsub("^gene-","PT-",r_cp[cp_matched])
rownames(rn1.cts) <- r_cp

rn1.seurat.object <- CreateSeuratObject(counts = rn1.cts, project = "RN1", min.cells = 5, min.features = 300)
rn1.seurat.object[["percent.mt"]] <- PercentageFeatureSet(rn1.seurat.object, pattern = "^MT-")
rn1.seurat.object[["percent.cp"]] <- PercentageFeatureSet(rn1.seurat.object, pattern = "^PT-")

VlnPlot(rn1.seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cp"), ncol = 4, pt.size=0)
plot1 <- FeatureScatter(rn1.seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(rn1.seurat.object, feature1 = "nCount_RNA", feature2 = "percent.cp")
plot3 <- FeatureScatter(rn1.seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 + plot3

rn1.filtered <- subset(rn1.seurat.object, subset = nFeature_RNA > 300 & 
                         percent.mt < 2 & percent.cp < 5 & nCount_RNA > 500 & nCount_RNA < 25000)

VlnPlot(rn1.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cp"), ncol = 4, pt.size=0)
plot1 <- FeatureScatter(rn1.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(rn1.filtered, feature1 = "nCount_RNA", feature2 = "percent.cp")
plot3 <- FeatureScatter(rn1.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 + plot3

rn1.filtered <- rn1.filtered %>% 
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(features = VariableFeatures(object = rn1.filtered), npcs = 100) %>%
  RunUMAP(dims = 1:75) %>%
  FindNeighbors(dims = 1:75) %>%
  FindClusters(resolution = 0.3)

plotHV <- VariableFeaturePlot(rn1.filtered)
plotPCA <- DimPlot(rn1.filtered, reduction = "pca") + NoLegend() + ggtitle("PCA of Normal") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
plotPC_n <- ElbowPlot(rn1.filtered, ndims = 100) + NoLegend() + ggtitle("Elbow Plot of the First 100 PCs (Normal)") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
plotRN1 <- UMAPPlot(object = rn1.filtered, label=T) + ggtitle('RN1')

sweep.res.rn1 <- paramSweep(rn1.filtered, PCs = 1:50, sct = FALSE)
sweep.stat.rn1 <- summarizeSweep(sweep.res.rn1, GT = FALSE)
bcmvn.rn1 <- find.pK(sweep.stat.rn1)
ggplot(bcmvn.rn1, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  ggtitle('RN1')

pK.rn1 <- bcmvn.rn1 %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK.rn1 <- as.numeric(as.character(pK.rn1[[1]]))
ann.rn1 <- rn1.filtered@meta.data$seurat_clusters
homotypic.rn1 <- modelHomotypic(ann.rn1)
nExp.rn1 <- round(0.08*nrow(rn1.filtered@meta.data))
nExp.rn1.adj <- round(nExp.rn1*(1-homotypic.rn1))
rn1.filtered <- doubletFinder(rn1.filtered, PCs = 1:50, pN = 0.25, pK = pK.rn1, nExp = nExp.rn1.adj, reuse.pANN = FALSE, sct = FALSE)

names(rn1.filtered@meta.data)
DimPlot(rn1.filtered, reduction = 'umap', group.by = 'DF.classifications_0.25_0.29_643') +
  ggtitle('RN1') + 
  xlab("UMAP 1") + 
  ylab("UMAP 2")

rn1.singlet <- subset(rn1.filtered, subset = DF.classifications_0.25_0.29_643 == "Singlet")

rd2.cts <- Read10X_h5(filename = "filtered_feature_bc_matrix_rd2.h5")

r_mt <- rownames(rd2.cts)
mt_table <- read_excel("data/mt_cp_geneList.xlsx", sheet = "mt") %>% as.data.frame()
mt_genes <- mt_table$external_gene_name
mt_matched <- r_mt %in% mt_genes
r_mt[mt_matched] <- gsub("^gene-","MT-",r_mt[mt_matched])
rownames(rd2.cts) <- r_mt

r_cp <- rownames(rd2.cts)
cp_table <- read_excel("data/mt_cp_geneList.xlsx", sheet = "cp") %>% as.data.frame()
cp_genes <- cp_table$external_gene_name
cp_matched <- r_cp %in% cp_genes
r_cp[cp_matched] <- gsub("^gene-","PT-",r_cp[cp_matched])
rownames(rd2.cts) <- r_cp

rd2.seurat.object <- CreateSeuratObject(counts = rd2.cts, project = "RD2", min.cells = 5, min.features = 300)
rd2.seurat.object[["percent.mt"]] <- PercentageFeatureSet(rd2.seurat.object, pattern = "^MT-")
rd2.seurat.object[["percent.cp"]] <- PercentageFeatureSet(rd2.seurat.object, pattern = "^PT-")
VlnPlot(rd2.seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cp"), ncol = 4, pt.size=0)
plot1 <- FeatureScatter(rd2.seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(rd2.seurat.object, feature1 = "nCount_RNA", feature2 = "percent.cp")
plot3 <- FeatureScatter(rd2.seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 + plot3

rd2.filtered <- subset(rd2.seurat.object, subset = nFeature_RNA > 300 &  
                         percent.mt < 2 & percent.cp < 5 & nCount_RNA > 500 & nCount_RNA < 25000)
VlnPlot(rd2.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cp"), ncol = 4, pt.size=0)
plot1 <- FeatureScatter(rd2.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(rd2.filtered, feature1 = "nCount_RNA", feature2 = "percent.cp")
plot3 <- FeatureScatter(rd2.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 + plot3

rd2.filtered <- rd2.filtered %>% 
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(features = VariableFeatures(object = rd2.filtered), npcs = 100) %>%
  RunUMAP(dims = 1:75) %>%
  FindNeighbors(dims = 1:75) %>%
  FindClusters(resolution = 0.3)
plotHV2 <- VariableFeaturePlot(rd2.filtered)
plotPCA2 <- DimPlot(rd2.filtered, reduction = "pca") + NoLegend() + ggtitle("PCA of Drought") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
plotPC2_n <- ElbowPlot(rd2.filtered, ndims = 100) + NoLegend() + ggtitle("Elbow Plot of the First 100 PCs (Drought)") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
plotRD2 <- UMAPPlot(object = rd2.filtered, label=T) + ggtitle('RD2')

sweep.res.rd2 <- paramSweep(rd2.filtered, PCs = 1:50, sct = FALSE)
sweep.stat.rd2 <- summarizeSweep(sweep.res.rd2, GT = FALSE)
bcmvn.rd2 <- find.pK(sweep.stat.rd2)
ggplot(bcmvn.rd2, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  ggtitle('RD2')

pK.rd2 <- bcmvn.rd2 %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK.rd2 <- as.numeric(as.character(pK.rd2[[1]])) # pk = 0.21
ann.rd2 <- rd2.filtered@meta.data$seurat_clusters
homotypic.rd2 <- modelHomotypic(ann.rd2)
nExp.rd2 <- round(0.08*nrow(rd2.filtered@meta.data))
nExp.rd2.adj <- round(nExp.rd2*(1-homotypic.rd2)) # nExp = 934

rd2.filtered <- doubletFinder(rd2.filtered, PCs = 1:50, pN = 0.25, pK = pK.rd2, nExp = nExp.rd2.adj, reuse.pANN = FALSE, sct = FALSE)

names(rd2.filtered@meta.data)
DimPlot(rd2.filtered, reduction = 'umap', group.by = 'DF.classifications_0.25_0.26_1251') +
  ggtitle('RD2') + 
  xlab("UMAP 1") + 
  ylab("UMAP 2")
rd2.singlet <- subset(rd2.filtered, subset = DF.classifications_0.25_0.26_1251 == "Singlet")
