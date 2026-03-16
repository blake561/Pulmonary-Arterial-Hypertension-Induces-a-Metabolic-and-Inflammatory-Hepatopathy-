##Install Packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
BiocManager::install("DESeq2")

##Load Packages

library(Seurat)
library(ggplot2)
library(SeuratData)
library(EnhancedVolcano)
library(DESeq2)
library(BiocManager)
library(tidyr)
library(gridExtra)
library(tidyverse)

##Load Data (From Integration and Clustering Sheet)

total.mtx_deg <- readRDS("/users/1/blake561/Fontan_total.mtx_final")

Fontan <- total.mtx_deg

#Filtering Out ATAC Data

rna_counts <- GetAssayData(total.mtx_deg, assay = "RNA", slot = "counts")
gene_features <- !grepl(":", rownames(rna_counts))
gene_only_counts <- rna_counts[gene_features, ]
filtered_obj <- CreateSeuratObject(counts = gene_only_counts, meta.data = total.mtx_deg@meta.data)
filtered_counts <- GetAssayData(filtered_obj, assay = "RNA", slot = "counts")
filtered_metadata <- filtered_obj@meta.data
filtered_obj_new <- CreateSeuratObject(counts = gene_only_counts, meta.data = filtered_obj@meta.data)

saveRDS(filtered_obj_new, file ="/users/1/blake561/Fontan_total.mtx_filteredNoATAC")

filtered_obj_new <- readRDS('/users/1/blake561/Fontan_total.mtx_filteredNoATAC')

#Separating By Cell Type for DEG Analysis

Fontan$subcelltype <- Idents(Fontan)

Hepatocyte <- subset(Fontan, subset = subcelltype == "Hepatocyte")

Endothelial <- subset(Fontan, subset = subcelltype == "Endothelial Cell")

HSC <- subset(Fontan, subset = subcelltype == "HSC/mFB")

Macrophage <- subset(Fontan, subset = subcelltype == "Macrophage")

##UMAPs by Cell Types##

#Not Separated by Clusters

DimPlot(Hepatocyte, reduction = "umap", split.by = "orig.ident", cols =c('Hepatocyte' = 'chartreuse3'), raster = FALSE) 

DimPlot(Endothelial, reduction = "umap", split.by = "orig.ident", cols =c("Endothelial Cell" = 'deepskyblue'), raster = FALSE) 

DimPlot(Macrophage, reduction = "umap", split.by = "orig.ident", cols =c('Macrophage' = 'violet'), raster = FALSE) 

DimPlot(HSC, reduction = "umap", split.by = "orig.ident", cols =c('HSC/mFB' = 'slateblue1'), raster = FALSE) 

#DEG Analysis

#Hepatocytes

Idents(Hepatocyte) <- Hepatocyte$orig.ident

hep_deg <- FindMarkers(
  object = Hepatocyte,
  ident.1 = "Fontan",
  ident.2 = "Fontan Control",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25 
)

hep_deg <- hep_deg[
  (abs(hep_deg$avg_log2FC) >= 0.25) &
    (hep_deg$p_val_adj < 0.05) &
    !grepl(":", rownames(hep_deg)),
]

write.csv(hep_deg, "hep_Fontan_FC_0.25.csv", row.names = TRUE)

#Endothelial 

Idents(Endothelial) <- Endothelial$orig.ident

endo_deg <- FindMarkers(
  object = Endothelial,
  ident.1 = "Fontan",
  ident.2 = "Fontan Control",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25 
)

endo_deg <- endo_deg[
  (abs(endo_deg$avg_log2FC) >= 0.25) &
    (endo_deg$p_val_adj < 0.05) &
    !grepl(":", rownames(endo_deg)),
]

write.csv(endo_deg, "endo_Fontan_FC_0.25.csv", row.names = TRUE)

#HSC

Idents(HSC) <- HSC$orig.ident

hsc_deg <- FindMarkers(
  object = HSC,
  ident.1 = "Fontan",
  ident.2 = "Fontan Control",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25 
)

hsc_deg <- hsc_deg[
  (abs(hsc_deg$avg_log2FC) >= 0.25) &
    (hsc_deg$p_val_adj < 0.05) &
    !grepl(":", rownames(hsc_deg)),
]

write.csv(hsc_deg, "hsc_Fontan_FC_0.25.csv", row.names = TRUE)

#Macrophage

Idents(Macrophage) <- Macrophage$orig.ident

mac_deg <- FindMarkers(
  object = Macrophage,
  ident.1 = "Fontan",
  ident.2 = "Fontan Control",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25 
)

mac_deg <- mac_deg[
  (abs(mac_deg$avg_log2FC) >= 0.25) &
    (mac_deg$p_val_adj < 0.05) &
    !grepl(":", rownames(mac_deg)),
]

write.csv(mac_deg, "mac_Fontan_FC_0.25.csv", row.names = TRUE)

# #Creating Matrix for DEG
# 
# cts <- AggregateExpression(filtered_obj,
#                            group.by = c("annotation", "sample_num"),
#                            assays = "RNA",
#                            slot = "counts",
#                            return.seurat = FALSE)
# 
# cts <- cts$RNA
# cts.t <- t(cts)
# cts.t <- as.data.frame(cts.t)
# splitRows <- gsub('_.*', '', rownames(cts.t))
# 
# cts.split <- split.data.frame(cts.t,
#                               f = factor(splitRows))
# 
# cts.split.modified <- lapply(cts.split, function(x){
#   rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
#   t(x)
#   
# })
# 
# 
# 

## Differential Analysis - Cell Types

# #Hepatocyte DEGs (Whole Cluster)
# 
# counts_HP <- cts.split.modified$Hepatocyte
# colData <- data.frame(samples = colnames(counts_HP))
# 
# colData$condition <- c('Control','Control','Fontan','Fontan','Fontan','Fontan')
# 
# colData <- column_to_rownames(colData, var = 'samples')
# 
# dds_HP <- DESeqDataSetFromMatrix(countData = counts_HP,
#                                  colData = colData,
#                                  design = ~condition)
# 
# keep <- rowSums(counts(dds_HP)) >=10
# dds_HP <- dds_HP[keep,]
# 
# dds_HP <- DESeq(dds_HP)
# 
# #These comparisons are c("condition", "numerator", "denominator")
# 
# Fontan_vs_Control_HP_res <- results(dds_HP, contrast =c("condition", "Fontan", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in Fontan Relative to Control 
# 
# #Saving File for Pathway Analysis
# 
# write.csv(Fontan_vs_Control_HP_res, file = "MB250403 Pseudobulk HP DE Fontan_vs_Control")
# d1 <- as.data.frame(Fontan_vs_Control_HP_res)
# d2 <- subset(d1, subset = padj <= 0.05 &
#                abs(log2FoldChange) >= 0.25) #padj is FDR
# write.csv(d2, file = "Fontan_vs_Control Hepatocyte Pathway Analysis.csv")
# 
# #Volcano Plot
# 
# EnhancedVolcano(Fontan_vs_Control_HP_res,
#                 lab = rownames(Fontan_vs_Control_HP_res),
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 pCutoff = 0.05,
#                 FCcutoff = 0.5,
#                 title = "Fontan vs Control Hepatocytes")
# 
# ##Macrophage Comparision## 
# 
# counts_MP <- cts.split.modified$Macrophage
# colData <- data.frame(samples = colnames(counts_MP))
# 
# colData$condition <- c('Control','Control','Fontan','Fontan','Fontan','Fontan')
# #View(colData)
# 
# colData <- column_to_rownames(colData, var = 'samples')
# 
# dds_MP <- DESeqDataSetFromMatrix(countData = counts_MP,
#                                  colData = colData,
#                                  design = ~condition)
# 
# keep <- rowSums(counts(dds_MP)) >=10
# dds_MP <- dds_MP[keep,]
# 
# dds_MP <- DESeq(dds_MP)
# resultsNames(dds_MP)
# 
# #These comparisons are c("condition", "numerator", "denominator")
# 
# Fontan_vs_Control_MP_res <- results(dds_MP, contrast =c("condition", "Fontan", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in Fontan Relative to Control 
# view(Fontan_vs_Control_MP_res)
# 
# #Saving File
# 
# write.csv(Fontan_vs_Control_MP_res, file = "MB250403 Pseudobulk MP DE Fontan_vs_Control")
# m1 <- as.data.frame(Fontan_vs_Control_MP_res)
# m2 <- subset(m1, subset = padj <= 0.05 &
#                abs(log2FoldChange) >= 0.25) #padj is FDR
# write.csv(m2, file = "Fontan_vs_Control Macrophage Pathway Analysis.csv")
# 
# #Volcano Plot
# 
# EnhancedVolcano(Fontan_vs_Control_MP_res,
#                 lab = rownames(Fontan_vs_Control_MP_res),
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 pCutoff = 0.05,
#                 FCcutoff = 0.5,
#                 title = "Fontan vs Control Macrophage")
# 
# ##Endothelial DEGs (Whole Cluster)
# 
# counts_EC <- cts.split.modified$Endothelial
# colData <- data.frame(samples = colnames(counts_EC))
# 
# colData$condition <- c('Control','Control','Fontan','Fontan','Fontan','Fontan')
# 
# colData <- column_to_rownames(colData, var = 'samples')
# 
# dds_EC <- DESeqDataSetFromMatrix(countData = counts_EC,
#                                  colData = colData,
#                                  design = ~condition)
# 
# keep <- rowSums(counts(dds_EC)) >=10
# dds_EC <- dds_EC[keep,]
# 
# dds_EC <- DESeq(dds_EC)
# resultsNames(dds_EC)
# rownames(results(dds_EC))
# 
# #Comparison 
# 
# Fontan_vs_Control_EC_res <- results(dds_EC, contrast =c("condition", "Fontan", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in Fontan Relative to Control 
# 
# #Saving File
# 
# write.csv(Fontan_vs_Control_EC_res, file = "MB250403 Pseudobulk EC DE Fontan_vs_Control")
# e1 <- as.data.frame(Fontan_vs_Control_EC_res)
# e2 <- subset(e1, subset = padj <= 0.05 &
#                abs(log2FoldChange) >= 0.25) #padj is FDR
# write.csv(e2, file = "Fontan_vs_Control Endothelial Pathway Analysis.csv")
# 
# #Volcano Plot
# 
# EnhancedVolcano(Fontan_vs_Control_EC_res,
#                 lab = rownames(Fontan_vs_Control_EC_res),
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 pCutoff = 0.05,
#                 FCcutoff = 0.5,
#                 title = "Fontan vs Control Endothelial")
# 
# ##HSC/mFB DEGs (Whole Cluster)
# 
# counts_FB <- cts.split.modified$`HSC/mFB`
# colData <- data.frame(samples = colnames(counts_FB))
# 
# colData$condition <- c('Control','Control','Fontan','Fontan','Fontan','Fontan')
# 
# colData <- column_to_rownames(colData, var = 'samples')
# 
# dds_FB <- DESeqDataSetFromMatrix(countData = counts_FB,
#                                  colData = colData,
#                                  design = ~condition)
# 
# keep <- rowSums(counts(dds_FB)) >=10
# dds_FB <- dds_FB[keep,]
# 
# dds_FB <- DESeq(dds_FB)
# resultsNames(dds_FB)
# 
# #Comparison 
# 
# Fontan_vs_Control_FB_res <- results(dds_FB, contrast =c("condition", "Fontan", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in Fontan Relative to Control 
# 
# #Saving File
# 
# write.csv(Fontan_vs_Control_FB_res, file = "MB250403 Pseudobulk FB DE Fontan_vs_Control")
# f1 <- as.data.frame(Fontan_vs_Control_FB_res)
# f2 <- subset(f1, subset = padj <= 0.05 &
#                abs(log2FoldChange) >= 0.25) #padj is FDR
# write.csv(f2, file = "Fontan_vs_Control HSC mFB Pathway Analysis.csv")
# 
# #Volcano Plot
# 
# EnhancedVolcano(Fontan_vs_Control_FB_res,
#                 lab = rownames(Fontan_vs_Control_FB_res),
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 pCutoff = 0.05,
#                 FCcutoff = 0.5,
#                 title = "Fontan vs Control HSC/mFB")
# 
# #Lymphocyte Comparision 
# 
# ##Lymphocyte DEGs (Whole Cluster)
# 
# counts_LC <- cts.split.modified$Lymphocyte
# colData <- data.frame(samples = colnames(counts_LC))
# 
# colData$condition <- c('Control','Control','Fontan','Fontan','Fontan','Fontan')
# 
# colData <- column_to_rownames(colData, var = 'samples')
# 
# dds_LC <- DESeqDataSetFromMatrix(countData = counts_LC,
#                                  colData = colData,
#                                  design = ~condition)
# 
# keep <- rowSums(counts(dds_LC)) >=10
# dds_LC <- dds_LC[keep,]
# 
# dds_LC <- DESeq(dds_LC)
# resultsNames(dds_LC)
# 
# #Comparison 
# 
# Fontan_vs_Control_LC_res <- results(dds_LC, contrast =c("condition", "Fontan", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in Fontan Relative to Control 
# #view(Fontan_vs_Control_LC_res)
# 
# #Saving File
# 
# write.csv(Fontan_vs_Control_LC_res, file = "MB250403 Pseudobulk LC DE Fontan_vs_Control")
# l1 <- as.data.frame(Fontan_vs_Control_LC_res)
# l2 <- subset(l1, subset = padj <= 0.05 &
#                abs(log2FoldChange) >= 0.25) #padj is FDR
# write.csv(l2, file = "Fontan_vs_Control Lymphocyte Pathway Analysis.csv")
# 
# #Volcano Plot
# 
# EnhancedVolcano(Fontan_vs_Control_LC_res,
#                 lab = rownames(Fontan_vs_Control_LC_res),
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 pCutoff = 0.05,
#                 FCcutoff = 0.5,
#                 title = "Fontan vs Control Lymphocyte")
# 
# ##Cholangiocytes DEGs (Whole Cluster)
# 
# counts_CL <- cts.split.modified$Cholangiocyte
# colData <- data.frame(samples = colnames(counts_CL))
# 
# colData$condition <- c('Control','Control','Fontan','Fontan','Fontan','Fontan')
# 
# colData <- column_to_rownames(colData, var = 'samples')
# 
# dds_CL <- DESeqDataSetFromMatrix(countData = counts_CL,
#                                  colData = colData,
#                                  design = ~condition)
# 
# keep <- rowSums(counts(dds_CL)) >=10
# dds_CL <- dds_CL[keep,]
# 
# dds_CL <- DESeq(dds_CL)
# resultsNames(dds_CL)
# 
# #Comparison 
# 
# Fontan_vs_Control_CL_res <- results(dds_CL, contrast =c("condition", "Fontan", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in Fontan Relative to Control 
# #view(Fontan_vs_Control_CL_res)
# 
# #Saving File
# 
# write.csv(Fontan_vs_Control_CL_res, file = "MB250403 Pseudobulk CL DE Fontan_vs_Control")
# c1 <- as.data.frame(Fontan_vs_Control_CL_res)
# c2 <- subset(c1, subset = padj <= 0.05 &
#                abs(log2FoldChange) >= 0.25) #padj is FDR
# write.csv(c2, file = "Fontan_vs_Control Cholangiocyte Pathway Analysis.csv")
# 
# #Volcano Plot
# 
# EnhancedVolcano(Fontan_vs_Control_CL_res,
#                 lab = rownames(Fontan_vs_Control_CL_res),
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 pCutoff = 0.05,
#                 FCcutoff = 0.5,
#                 title = "Fontan vs Control Cholangiocyte")
# 
# 
# #Wilcoxon Test
# 
# cells_hepatocytes <- rownames(filtered_obj@meta.data)[filtered_obj@meta.data$`cell_identity` == "Hepatocyte"]
# hepatocytes <- subset(filtered_obj, cells = cells_hepatocytes)
# 
# Idents(hepatocytes) <- hepatocytes@meta.data$orig.ident
# hep_deg <- FindMarkers(
#   object = hepatocytes,
#   ident.1 = "Fontan",
#   ident.2 = "Fontan Control",
#   test.use = "wilcox",
#   logfc.threshold = 1.5,
#   min.pct = 0.25 #Original run was 10% instead of 25%
# )
# 
# hep_deg <- hep_deg[abs(hep_deg$avg_log2FC) > 1.5 & hep_deg$p_val_adj < 0.00001, ]
# 
# 
# cells_hsc_orig <- rownames(total.mtx_deg@meta.data)[total.mtx_deg@meta.data$`cell_identity` == "HSC/mFB"]
# hsc_orig <- subset(total.mtx_deg, cells = cells_hsc_orig)
# 
# 
# cells_hsc <- rownames(filtered_obj@meta.data)[filtered_obj@meta.data$`cell_identity` == "HSC/mFB"]
# hsc <- subset(filtered_obj, cells = cells_hsc)
# 
# Idents(hsc) <- hsc@meta.data$orig.ident
# 
# 
# Idents(hsc) <- hsc@meta.data$orig.ident
# hsc_deg <- FindMarkers(
#   object = hsc,
#   ident.1 = "Fontan",
#   ident.2 = "Fontan Control",
#   test.use = "wilcox",
#   logfc.threshold = 1,
#   min.pct = 0.25 #Original run was 10% instead of 25%
# )
# 
# hsc_deg <- hsc_deg[abs(hsc_deg$avg_log2FC) > 1.5 & hsc_deg$p_val_adj < 0.00001, ]
# 
# write.csv(hsc_deg, "hsc_Fontan_0.25.csv", row.names = TRUE)
# 
# #Fontan Methods Test 
# 
# cells_hsc <- rownames(filtered_obj_new@meta.data)[filtered_obj_new@meta.data$cell_identity == "HSC/mFB"]
# hsc_new <- subset(filtered_obj_new, cells = cells_hsc)
# 
# # Normalize new object:
# hsc_new <- NormalizeData(hsc_new, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# 
# deg_results <- FindMarkers(
#   object = hsc,
#   ident.1 = "Fontan",
#   ident.2 = "Fontan Control",
#   test.use = "LR",                      # Use the likelihood ratio test
#   logfc.threshold = log2(1.5),          # Only test genes with |log2 fold-change| > log2(1.5)
#   min.pct = 0.25
# )
# 
# 
# deg_results <- deg_results[deg_results$p_val < 1e-5, ]
# 
# deg_results <- deg_results[abs(deg_results$avg_log2FC) > log2(1.5), ]
# 
# write.csv(deg_results, "HSC_DEG_Fontan_vs_Control_testMB250409.csv", row.names = TRUE)
# 
