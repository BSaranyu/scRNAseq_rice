library(Seurat)
library(SeuratObject)
library(dplyr)
library(patchwork)
library(DoubletFinder)
library(harmony)
library(ggplot2)
library(data.table)
library(tidyverse)
library(readxl)

mergedDat <- merge(rn1.singlet, y = rd2.singlet, project = "Rice Leaf", merge.data = TRUE)
names(mergedDat@meta.data)[names(mergedDat@meta.data) == "orig.ident"] <- "conditions"

mergedDat <- SCTransform(mergedDat)
mergedDat <- RunPCA(object = mergedDat, feature = VariableFeatures(object = mergedDat), npcs = 100)
mergedDat <- RunHarmony(object = mergedDat, group.by.vars = 'conditions', dim = 1:75)
mergedDat <- RunUMAP(object = mergedDat, reduction = "harmony", dim = 1:75)
mergedDat <- FindNeighbors(object = mergedDat, dims = 1:75)
mergedDat <- FindClusters(object = mergedDat, resolution = 0.70)
names(mergedDat@meta.data)[names(mergedDat@meta.data) == "seurat_clusters"] <- "clusters"

clusterHarmony <- DimPlot(mergedDat, reduction = 'umap', group.by = 'clusters', label = TRUE)
conditionHarmony <- DimPlot(mergedDat, reduction = 'umap', group.by = 'conditions') 
clusterHarmony|conditionHarmony
