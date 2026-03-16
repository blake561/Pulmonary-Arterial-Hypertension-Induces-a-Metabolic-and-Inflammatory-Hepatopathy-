
#Load Packages

library(BiocManager)
library(tidyverse)
library(Seurat)
library(gridExtra)
library(Azimuth)
library(harmony)
library(BiocManager)
library(multtest)
library(metap)

library(Seurat)
library(ggplot2)
library(SeuratData)
library(EnhancedVolcano)
library(DESeq2)
library(BiocManager)
library(tidyr)
library(gridExtra)
library(tidyverse)

#### Load Data ####

total.mtx <- readRDS(file = "/users/1/blake561/total_matrix_NASH_filt")

# Check it
ncol(total.mtx)
table(gsub("_.*", "", colnames(total.mtx)))

#Remove Column Artifacts from Doublet Finder
total.mtx$pANN_0.25_0.02_2201 <- NULL #N1
total.mtx$DF.classifications_0.25_0.02_2201 <- NULL
total.mtx$pANN_0.25_0.005_1868 <- NULL #N2
total.mtx$DF.classifications_0.25_0.005_1868 <- NULL
total.mtx$pANN_0.25_0.11_523 <- NULL #N3
total.mtx$DF.classifications_0.25_0.11_523 <- NULL
total.mtx$pANN_0.25_0.005_8453 <- NULL #CN1
total.mtx$DF.classifications_0.25_0.005_8453 <- NULL
total.mtx$pANN_0.25_0.005_919 <- NULL #CN2
total.mtx$DF.classifications_0.25_0.005_919 <- NULL
total.mtx$pANN_0.25_0.17_411 <- NULL #CN3
total.mtx$DF.classifications_0.25_0.17_411 <- NULL

total.mtx_filtered <- total.mtx 

#### Standard Workflow / Pre-processing #### 
total.mtx_filtered <- NormalizeData(object = total.mtx)
total.mtx_filtered <- FindVariableFeatures(object = total.mtx_filtered)
total.mtx_filtered <- ScaleData(object = total.mtx_filtered)
total.mtx_filtered <- RunPCA(object = total.mtx_filtered)

print(total.mtx_filtered[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(total.mtx_filtered, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(total.mtx_filtered)+
  ggtitle("Total Matrix Elbow Plot")

total.mtx_filtered <- FindNeighbors(object = total.mtx_filtered, dims = 1:12)
total.mtx_filtered <- FindClusters(object = total.mtx_filtered, dims = 1:12)
total.mtx_filtered <- RunUMAP(object = total.mtx_filtered, dims = 1:12)

DimPlot(total.mtx_filtered, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)
DimPlot(total.mtx_filtered, reduction = "umap", group.by = "orig.ident", raster = FALSE, label = TRUE)

# Create sample column
total.mtx_filtered$sample_id <- rownames(total.mtx_filtered@meta.data)
total.mtx_filtered@meta.data <- separate(total.mtx_filtered@meta.data, col = "sample_id", into = c("sample", "barcode"),
                                         sep = "_")


p1 <- DimPlot(total.mtx_filtered, reduction = "umap", group.by = 'sample')
p2 <- DimPlot(total.mtx_filtered, reduction = "umap")
grid.arrange(p1,p2, ncol = 2)

# Run Harmony
total.mtx_harmony <- total.mtx_filtered %>%
  RunHarmony(group.by.vars = 'sample', plot_convergence = FALSE)

# Check harmony was added
total.mtx_harmony@reductions
total.mtx_harmony.embed <- Embeddings(total.mtx_harmony, "harmony")
total.mtx_harmony.embed[1:10,1:10]

# Run downstream - ADD ASSIGNMENT HERE and fix dims
total.mtx_harmony <- total.mtx_harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:12) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:12) %>%
  FindClusters(resolution = 0.5)

# Visualize
p3 <- DimPlot(total.mtx_harmony, reduction = "umap", group.by = "sample") +
  ggtitle("Integrated With Harmony")
p4 <- DimPlot(total.mtx_harmony, reduction = "umap") +
  ggtitle("Clusters")

p3 + p4

saveRDS(total.mtx_harmony, file = "NASH_total.mtx_harmony.rds")
total.mtx_harmony <- readRDS(file = "NASH_total.mtx_harmony.rds")

#### Cluster Identification & Annotation #### 

total.mtx_harmony <- FindClusters(total.mtx_harmony, resolution = c(0.1,0.2,0.25,0.3,0.5))
total.mtx_harmony <- FindClusters(total.mtx_harmony, resolution = c(0.25))

DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE, raster = FALSE) 
DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.2", label = TRUE, raster = FALSE)
DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.25", label = TRUE, raster = FALSE) 
DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.3", label = TRUE, raster = FALSE) 
DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.5", label = TRUE, raster = FALSE) 

Idents(total.mtx_harmony) <- total.mtx_harmony$RNA_snn_res.0.25 

total.mtx_harmony <- JoinLayers(total.mtx_harmony)
total.mtx_harmony <- RunAzimuth(total.mtx_harmony, reference = "/users/1/blake561/Liver Reference")

#UMAPs for Control vs. Fontan

Idents(total.mtx_harmony) <- "predicted.celltype.l1"
p5 <- DimPlot(total.mtx_harmony, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, raster = FALSE)
p6 <- DimPlot(total.mtx_harmony, reduction = "umap", group.by = "predicted.celltype.l2", label = TRUE, raster = FALSE)

p5+p6

#UMAPs for Control vs. Fontan

p5 <- DimPlot(total.mtx_harmony, reduction = "umap", group.by = "seurat_clusters", split.by ='orig.ident',label = TRUE, raster = FALSE)

#Checking Confidence in Predictions

table(Idents(total.mtx_harmony), total.mtx_harmony$predicted.celltype.l1.score < 0.5) #Most clusters are confident predictions, plasma cells and erythrocytes are a bit uncertain

Idents(total.mtx_harmony) <- "predicted.celltype.l2"
table(Idents(total.mtx_harmony), total.mtx_harmony$predicted.celltype.l2.score < 0.5) #Heterogeneity for T/NK cells and myeloid cells, rest are homogeneous
table(Idents(total.mtx_harmony), total.mtx_harmony$predicted.celltype.l2) #Heterogeneity for T/NK cells and myeloid cells, rest are homogeneous


total.mtx_harmony$condition <- total.mtx_harmony$orig.ident

#Overlapped on Same Graph 

DimPlot(total.mtx_harmony, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE) 

#Set Resolution of Your Choice as Seurat Cluster Numbers and Identity as Seurat Clusters

total.mtx_harmony$seurat_clusters <- total.mtx_harmony$RNA_snn_res.0.25 
Idents(total.mtx_harmony) <- total.mtx_harmony$seurat_clusters

DimPlot(total.mtx_harmony, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE) 

### FindConservedMarkers() / Begin Annotation for this Dataset

c1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 0, only.pos = TRUE, grouping.var = 'orig.ident')
c1_genes <- rownames(head(c1_markers,20))
FeaturePlot(total.mtx_harmony, features = c1_genes, min.cutoff = "q10")

c2_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 1, only.pos = TRUE, grouping.var = 'orig.ident')
c2_genes <- rownames(head(c2_markers,20))
FeaturePlot(total.mtx_harmony, features = c2_genes, min.cutoff = "q10")

c3_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 2, only.pos = TRUE, grouping.var = 'orig.ident')
c3_genes <- rownames(head(c3_markers,20))
FeaturePlot(total.mtx_harmony, features = c3_genes, min.cutoff = "q10")

c4_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 3, only.pos = TRUE, grouping.var = 'orig.ident')
c4_genes <- rownames(head(c4_markers,20))
FeaturePlot(total.mtx_harmony, features = c4_genes, min.cutoff = "q10")

c5_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 4, only.pos = TRUE, grouping.var = 'orig.ident')
c5_genes <- rownames(head(c5_markers,20))
FeaturePlot(total.mtx_harmony, features = c5_genes, min.cutoff = "q10")

c6_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 5, only.pos = TRUE, grouping.var = 'orig.ident')
c6_genes <- rownames(head(c6_markers,20))
FeaturePlot(total.mtx_harmony, features = c6_genes, min.cutoff = "q10")

c7_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 6, only.pos = TRUE, grouping.var = 'orig.ident')
c7_genes <- rownames(head(c7_markers,20))
FeaturePlot(total.mtx_harmony, features = c7_genes, min.cutoff = "q10")

c8_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 7, only.pos = TRUE, grouping.var = 'orig.ident') #Not present in PAH, only Control
c8_genes <- rownames(head(c8_markers,20))
FeaturePlot(total.mtx_harmony, features = c8_genes, min.cutoff = "q10")

c9_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 8, only.pos = TRUE, grouping.var = 'orig.ident') #Not present in PAH, only Control
c9_genes <- rownames(head(c9_markers,20))
FeaturePlot(total.mtx_harmony, features = c9_genes, min.cutoff = "q10")

c10_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 9, only.pos = TRUE, grouping.var = 'orig.ident')
c10_genes <- rownames(head(c10_markers,20))
FeaturePlot(total.mtx_harmony, features = c10_genes, min.cutoff = "q10")

c11_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 10, only.pos = TRUE, grouping.var = 'orig.ident')
c11_genes <- rownames(head(c11_markers,20))
FeaturePlot(total.mtx_harmony, features = c11_genes, min.cutoff = "q10")

c12_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 11, only.pos = TRUE, grouping.var = 'orig.ident')
c12_genes <- rownames(head(c12_markers,20))
FeaturePlot(total.mtx_harmony, features = c12_genes, min.cutoff = "q10")

c13_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 12, only.pos = TRUE, grouping.var = 'orig.ident')
c13_genes <- rownames(head(c13_markers,20))
FeaturePlot(total.mtx_harmony, features = c13_genes, min.cutoff = "q10")

c14_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 13, only.pos = TRUE, grouping.var = 'orig.ident')
c14_genes <- rownames(head(c14_markers,20))
FeaturePlot(total.mtx_harmony, features = h14_genes, min.cutoff = "q10")

#Manual Annotation of Clusters and Removing Noise

Idents(total.mtx_harmony) <- total.mtx_harmony$seurat_clusters #Reassinging Numerical Labels for Renaming

total.mtx_harmony <- RenameIdents(total.mtx_harmony, "0" = "Hepatocyte")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "1" = "Hepatocyte")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "2" = "Hepatocyte")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "3" = "Endothelial Cell")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "4" = "Macrophage")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "5" = "Lymphocyte")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "6" = "Cholangiocyte")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "7" = "HSC/mFB")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "8" = "Endothelial Cell")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "10" = "HSC/mFB")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "11" = "Hepatocyte")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "12" = "Endothelial Cell")

total.mtx_harmony_filt <- total.mtx_harmony

total.mtx_cluster <- subset(total.mtx_harmony_filt, idents = c('9'), invert = TRUE) #Excluding Junk Clusters

#Confirm with UMAPs

total.mtx_cluster$subcelltype <- Idents(total.mtx_cluster) 

DimPlot(total.mtx_cluster, reduction = "umap", group.by = "subcelltype", split.by = 'orig.ident', label = TRUE, raster = FALSE) 

total.mtx_cluster$sample_num <- paste0(total.mtx_cluster$subcelltype, total.mtx_cluster$sample)

total.mtx_deg <- total.mtx_cluster

total.mtx_deg$annotation <- Idents(total.mtx_deg)

total.mtx_deg$sample_num <- paste0(total.mtx_deg$subcelltype, total.mtx_deg$sample)

###

# Check what columns you have
head(total.mtx_cluster@meta.data)

# Create abundance table by sample
abundance_by_sample <- total.mtx_cluster@meta.data %>%
  group_by(sample, subcelltype) %>%
  summarise(count = n()) %>%
  group_by(sample) %>%
  mutate(total_cells = sum(count),
         relative_abundance = count / total_cells * 100) %>%
  ungroup()

# View it
print(abundance_by_sample)

# Wide format (samples as columns, cell types as rows)
abundance_wide_sample <- abundance_by_sample %>%
  select(sample, subcelltype, relative_abundance) %>%
  pivot_wider(names_from = sample, values_from = relative_abundance, values_fill = 0)

# Create condition column (CN = Control, N = NASH)
total.mtx_cluster$condition <- ifelse(grepl("^CN", total.mtx_cluster$sample), "Control", "NASH")

# Abundance by condition (group level)
abundance_by_condition <- total.mtx_cluster@meta.data %>%
  group_by(condition, subcelltype) %>%
  summarise(count = n()) %>%
  group_by(condition) %>%
  mutate(total_cells = sum(count),
         relative_abundance = count / total_cells * 100) %>%
  ungroup()

# Wide format for condition
abundance_wide_condition <- abundance_by_condition %>%
  select(condition, subcelltype, relative_abundance) %>%
  pivot_wider(names_from = condition, values_from = relative_abundance, values_fill = 0)

# Save CSVs
write.csv(abundance_by_sample, "NASH_celltype_abundance_by_sample_long.csv", row.names = FALSE)
write.csv(abundance_wide_sample, "NASH_celltype_abundance_by_sample_wide.csv", row.names = FALSE)
write.csv(abundance_by_condition, "NASH_celltype_abundance_by_condition_long.csv", row.names = FALSE)
write.csv(abundance_wide_condition, "NASH_celltype_abundance_by_condition_wide.csv", row.names = FALSE)

# Combined summary table
abundance_combined <- abundance_by_sample %>%
  select(sample, subcelltype, relative_abundance) %>%
  pivot_wider(names_from = sample, values_from = relative_abundance, values_fill = 0) %>%
  left_join(abundance_wide_condition, by = "subcelltype")

write.csv(abundance_combined, "NASH_celltype_abundance_combined.csv", row.names = FALSE)

saveRDS(total.mtx_cluster, file = "NASH_total.mtx_final")


