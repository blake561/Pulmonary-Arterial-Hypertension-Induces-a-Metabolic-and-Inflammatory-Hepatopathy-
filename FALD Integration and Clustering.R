#Library / Install Packages
install.packages("tidyverse")
install.packages("Seurat")
install.packages("gridExtra")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TFBSTools")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
install.packages("harmony")
install.packages("metap")

install.packages("BiocManager")

#Seurat Disk/Azimuth Dependencies and Installation

install.packages("hdf5r")
install.packages("cli")
install.packages("crayon")
install.packages("Matrix")
install.packages("R6")
install.packages("rlang")
install.packages("stringi")
install.packages("withr")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github('satijalab/azimuth', ref = 'master')

#Metap Dependencies

install.packages("lattice")
install.packages("TFisher")
BiocManager::install('multtest')
install.packages("metap")
BiocManager::install('multtest')
install.packages("https://cran.r-project.org/src/contrib/Archive/metap/metap_1.0.tar.gz", repos=NULL, type="source") #Newer version of metap requires qqconf, this is an older version to circumvent that requirement! 

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

#### Load Data ####

total.mtx <- readRDS(file = "/users/1/blake561/total_matrix_fontan")

#Remove Column Artifacts from Doublet Finder

total.mtx$pANN_0.25_0.005_304 <- NULL #F1
total.mtx$DF.classifications_0.25_0.005_304 <- NULL
total.mtx$pANN_0.25_0.005_182 <- NULL #F2
total.mtx$DF.classifications_0.25_0.005_182 <- NULL
total.mtx$pANN_0.25_0.005_764 <- NULL #F3
total.mtx$DF.classifications_0.25_0.005_764 <- NULL
total.mtx$pANN_0.25_0.005_101 <- NULL #F4
total.mtx$DF.classifications_0.25_0.005_101 <- NULL
total.mtx$pANN_0.25_0.005_337 <- NULL #CF1
total.mtx$DF.classifications_0.25_0.005_337 <- NULL
total.mtx$pANN_0.25_0.01_1807 <- NULL #CF2
total.mtx$DF.classifications_0.25_0.01_1807 <- NULL


view(total.mtx@meta.data)

VlnPlot(total.mtx_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

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

#### Performing Integration via Harmony #### 

total.mtx_harmony <- total.mtx_filtered %>%
  RunHarmony(group.by.vars = 'sample', plot_convergence = FALSE)

total.mtx_harmony@reductions

total.mtx_harmony.embed <- Embeddings(total.mtx_harmony, "harmony")
total.mtx_harmony.embed[1:10,1:10]

total.mtx_harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:12)%>%
  FindNeighbors(reduction = 'harmony', dims = 1:12)%>%
  FindClusters(resolution = 0.5)

p3 <- DimPlot(total.mtx_harmony, reduction = "umap", group.by = "sample")+
  ggtitle("Integrated With Harmony")
p4 <- DimPlot(total.mtx_harmony, reduction = "umap")

grid.arrange(p1,p3, ncol = 2)

#### Cluster Identification & Annotation #### 

total.mtx_harmony <- FindClusters(total.mtx_harmony, resolution = c(0.1,0.25,0.3,0.35))

DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE, raster = FALSE) 
DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.25", label = TRUE, raster = FALSE) 
DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.3", label = TRUE, raster = FALSE) 
DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.35", label = TRUE, raster = FALSE) 

#Load in Human Liver Azimuth Reference - Need to Download Files from Azimuth Database and Upload to MSI in Reference Folder

Idents(total.mtx_harmony) <- total.mtx_harmony$RNA_snn_res.0.4 #Do print afterwards to confirm 23 clusters, because that is what we see at resolution of 0.7

liverref <- LoadReference("/users/1/blake561/Reference")

total.mtx_harmony <- JoinLayers(total.mtx_harmony)
total.mtx_harmony <- RunAzimuth(total.mtx_harmony, reference = "/users/1/blake561/Reference")

Idents(total.mtx_harmony) <- "predicted.celltype.l1"
p5 <- DimPlot(total.mtx_harmony, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, raster = FALSE)

#Save the File/Checkpoint

saveRDS(total.mtx_harmony, file = "total.mtx_harmony_fontan")
total.mtx_harmony <- readRDS("/users/1/blake561/total.mtx_harmony_fontan")

#UMAPs for Control vs. Fontan

#Grouping by Condition

total.mtx_harmony$condition <- total.mtx_harmony$orig.ident

total.mtx_fontan <- subset(total.mtx_harmony, condition == "Fontan")
p1 <- DimPlot(total.mtx_fontan, reduction = "umap", group.by = "condition", raster = FALSE, label = TRUE)

total.mtx_FC <- subset(total.mtx_harmony, condition == "Fontan Control")
p2 <- DimPlot(total.mtx_FC, reduction = "umap", group.by = "condition", raster = FALSE, label = TRUE)

#Grouping by Cell Type 

#Basic Prediction

total.mtx_fontan <- subset(total.mtx_harmony, condition == "Fontan")
p1 <- DimPlot(total.mtx_fontan, reduction = "umap", group.by = "predicted.celltype.l1", raster = FALSE, label = TRUE)

total.mtx_FC <- subset(total.mtx_harmony, condition == "Fontan Control")
p2 <- DimPlot(total.mtx_FC, reduction = "umap", group.by = "predicted.celltype.l1", raster = FALSE, label = TRUE)

#Complex Prediction

p1 <- DimPlot(total.mtx_fontan, reduction = "umap", group.by = "predicted.celltype.l2", raster = FALSE, label = TRUE)

p2 <- DimPlot(total.mtx_FC, reduction = "umap", group.by = "predicted.celltype.l2", raster = FALSE, label = TRUE)

#Overlapped on Same Graph 

DimPlot(total.mtx_harmony, reduction = "umap", group.by = "orig.ident", label = TRUE, raster = FALSE) 

#Checking Confidence in Predictions

table(Idents(total.mtx_harmony), total.mtx_harmony$predicted.celltype.l1.score < 0.5) #Most clusters are confident predictions, plasma cells and erythrocytes are a bit uncertain

Idents(total.mtx_harmony) <- "predicted.celltype.l2"
table(Idents(total.mtx_harmony), total.mtx_harmony$predicted.celltype.l2.score < 0.5) #Heterogeneity for T/NK cells and myeloid cells, rest are homogeneous
table(Idents(total.mtx_harmony), total.mtx_harmony$predicted.celltype.l2) #Heterogeneity for T/NK cells and myeloid cells, rest are homogeneous


total.mtx_harmony$condition <- total.mtx_harmony$orig.ident

#Set Resolution of Your Choice as Seurat Cluster Numbers and Identity as Seurat Clusters

total.mtx_harmony$seurat_clusters <- total.mtx_harmony$RNA_snn_res.0.35 #Resolution = 0.35
Idents(total.mtx_harmony) <- total.mtx_harmony$seurat_clusters

DimPlot(total.mtx_harmony, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)


#Hepatocytes

h1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 13, only.pos = TRUE, grouping.var = 'orig.ident')
h1_genes <- rownames(head(h1_markers,20))
#FeaturePlot(total.mtx_harmony, features = h1_genes, min.cutoff = "q10")

h2_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 0, only.pos = TRUE, grouping.var = 'orig.ident')
h2_genes <- rownames(head(h2_markers,20))
#FeaturePlot(total.mtx_harmony, features = h2_genes, min.cutoff = "q10")

h3_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 14, only.pos = TRUE, grouping.var = 'orig.ident')
h3_genes <- rownames(head(h3_markers,20))
#FeaturePlot(total.mtx_harmony, features = h3_genes, min.cutoff = "q10")

h4_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 1, only.pos = TRUE, grouping.var = 'orig.ident')
h4_genes <- rownames(head(h4_markers,20))
#FeaturePlot(total.mtx_harmony, features = h4_genes, min.cutoff = "q10")

h5_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 2, only.pos = TRUE, grouping.var = 'orig.ident')
h5_genes <- rownames(head(h5_markers,20))
#FeaturePlot(total.mtx_harmony, features = h5_genes, min.cutoff = "q10")

h6_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 3, only.pos = TRUE, grouping.var = 'orig.ident')
h6_genes <- rownames(head(h6_markers,20))
#FeaturePlot(total.mtx_harmony, features = h6_genes, min.cutoff = "q10")

h7_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 4, only.pos = TRUE, grouping.var = 'orig.ident')
h7_genes <- rownames(head(h7_markers,20))
#FeaturePlot(total.mtx_harmony, features = h7_genes, min.cutoff = "q10")

h8_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 6, only.pos = TRUE, grouping.var = 'orig.ident') #Not present in PAH, only Control
h8_genes <- rownames(head(h8_markers,20))
#FeaturePlot(total.mtx_harmony, features = h8_genes, min.cutoff = "q10")

#Endothelial Cells

e1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 9, only.pos = TRUE, grouping.var = 'orig.ident')
e1_genes <- rownames(head(e1_markers,20))
#FeaturePlot(total.mtx_harmony, features = e1_genes, min.cutoff = "q20")

e2_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 7, only.pos = TRUE, grouping.var = 'orig.ident')
e2_genes <- rownames(head(e2_markers,20))
#FeaturePlot(total.mtx_harmony, features = e2_genes, min.cutoff = "q10")

#HSC

hsc1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 11, only.pos = TRUE, grouping.var = 'orig.ident')
hsc1_genes <- rownames(head(hsc1_markers,20))
#FeaturePlot(total.mtx_harmony, features = hsc1_genes, min.cutoff = "q20")

#Myeloid

m1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 5, only.pos = TRUE, grouping.var = 'orig.ident')
m1_genes <- rownames(head(m1_markers,20))
#FeaturePlot(total.mtx_harmony, features = m1_genes, min.cutoff = "q10")

m2_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 12, only.pos = TRUE, grouping.var = 'orig.ident')
m2_genes <- rownames(head(m2_markers,20))
#FeaturePlot(total.mtx_harmony, features = m2_genes, min.cutoff = "q10")

#T/NK

t1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 8, only.pos = TRUE, grouping.var = 'orig.ident')
t1_genes <- rownames(head(t1_markers,20))
#FeaturePlot(total.mtx_harmony, features = t1_genes, min.cutoff = "q10")

#Cholangiocyte

c1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 10, only.pos = TRUE, grouping.var = 'orig.ident')
c1_genes <- rownames(head(c1_markers,20))
#FeaturePlot(total.mtx_harmony, features = c1_genes, min.cutoff = "q10")

#Marker Genes Plots for Confirmation

Idents(total.mtx_harmony) <- "predicted.celltype.l1"

DotPlot(total.mtx_harmony, features = h1_genes, cols = c('blue', 'red'))
DotPlot(total.mtx_harmony, features = h2_genes, cols = c('blue', 'red'))
DotPlot(total.mtx_harmony, features = h3_genes, cols = c('blue', 'red'))
DotPlot(total.mtx_harmony, features = h4_genes, cols = c('blue', 'red'))
DotPlot(total.mtx_harmony, features = h5_genes, cols = c('blue', 'red')) #Maybe Check Top 50
DotPlot(total.mtx_harmony, features = h6_genes, cols = c('blue', 'red'))
DotPlot(total.mtx_harmony, features = h7_genes, cols = c('blue', 'red'))
DotPlot(total.mtx_harmony, features = h8_genes, cols = c('blue', 'red'))

DotPlot(total.mtx_harmony, features = e1_genes, cols = c('blue', 'red'))
DotPlot(total.mtx_harmony, features = e2_genes, cols = c('blue', 'red'))

DotPlot(total.mtx_harmony, features = hsc1_genes, cols = c('blue', 'red'))

DotPlot(total.mtx_harmony, features = m1_genes, cols = c('blue', 'red'))
DotPlot(total.mtx_harmony, features = m2_genes, cols = c('blue', 'red'))

DotPlot(total.mtx_harmony, features = t1_genes, cols = c('blue', 'red'))

DotPlot(total.mtx_harmony, features = c1_genes, cols = c('blue', 'red'))

#Manual Annotation of harmonys and Removing Noise

total.mtx_harmony$seurat_clusters <- total.mtx_harmony$RNA_snn_res.0.35

Idents(total.mtx_harmony) <- total.mtx_harmony$seurat_clusters #Reassinging Numerical Labels for Renaming

total.mtx_harmony <- RenameIdents(total.mtx_harmony, "0" = "Hepatocyte 1")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "1" = "Hepatocyte 2")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "2" = "Hepatocyte 3")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "3" = "Hepatocyte 4")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "4" = "Hepatocyte 5")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "5" = "Macrophage 1")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "6" = "Hepatocyte 6")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "7" = "Endothelial 1")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "8" = "Lymphocyte 1")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "9" = "Endothelial 2")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "10" = "Cholangiocyte 1")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "11" = "HSC/mFB 1")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "12" = "Macrophage 2")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "13" = "Endothelial 3")

total.mtx_harmony$subcelltype <- Idents(total.mtx_harmony) 

total.mtx_final <- subset(total.mtx_harmony, idents = c('14'), invert = TRUE) #Excluding Junk clusters within the filter

#Creating General Cell Labeled filtereds (No filtered #)

total.mtx_final <- RenameIdents(total.mtx_final,"Hepatocyte 1" = 'Hepatocyte')
total.mtx_final <- RenameIdents(total.mtx_final, "Hepatocyte 2" = 'Hepatocyte')
total.mtx_final <- RenameIdents(total.mtx_final, "Hepatocyte 3" = 'Hepatocyte')
total.mtx_final <- RenameIdents(total.mtx_final, "Hepatocyte 4" = 'Hepatocyte')
total.mtx_final <- RenameIdents(total.mtx_final, "Hepatocyte 5" = 'Hepatocyte')
total.mtx_final <- RenameIdents(total.mtx_final, "Hepatocyte 6" = 'Hepatocyte')
total.mtx_final <- RenameIdents(total.mtx_final, "Macrophage 1" = 'Macrophage')
total.mtx_final <- RenameIdents(total.mtx_final, "Endothelial 1" = "Endothelial Cell")
total.mtx_final <- RenameIdents(total.mtx_final, "Lymphocyte 1" = "Lymphocyte")
total.mtx_final <- RenameIdents(total.mtx_final, "Endothelial 2" = "Endothelial Cell")
total.mtx_final <- RenameIdents(total.mtx_final, "Cholangiocyte 1" = "Cholangiocyte")
total.mtx_final <- RenameIdents(total.mtx_final, "HSC/mFB 1" = "HSC/mFB")
total.mtx_final <- RenameIdents(total.mtx_final, "Macrophage 2" = "Macrophage")
total.mtx_final <- RenameIdents(total.mtx_final, "Endothelial 3" = "Endothelial Cell")

total.mtx_final$cell_identity <- Idents(total.mtx_final) 

total.mtx_final$orig.ident <- factor(total.mtx_final$orig.ident, levels = c('Fontan Control', 'Fontan'))
DimPlot(total.mtx_final, reduction = "umap", group.by = "cell_identity", split.by ='orig.ident',label = TRUE, raster = FALSE) 

DimPlot(total.mtx_final, reduction = "umap", group.by = "cell_identity", label = TRUE, raster = FALSE) 

#Cell Type Proportions 

#Proportion Bar Plot for Bottom 

cell_table <- table(total.mtx_final$orig.ident, total.mtx_final$subcelltype)
cell_data <- as.data.frame(cell_table)
write_csv(cell_data, file ='Fontan Cleaned Cell Subcelltype Counts - MB250403.csv')

sample_table_main <- as.data.frame(table(total.mtx_final$sample, total.mtx_final$cell_identity))
write_csv(sample_table_main, file = "Fontan Main Celltype Counts - MB250403")

#Save File for DEG Analysis
total.mtx_final$sample_num <- paste0(total.mtx_final$cell_identity, total.mtx_final$sample)

saveRDS(total.mtx_final, file = "Fontan_total.mtx_final")
