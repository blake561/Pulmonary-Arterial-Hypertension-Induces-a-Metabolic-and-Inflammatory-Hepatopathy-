#Required/Downloaded Packages for Pipeline - Put here
install.packages("tidyverse")
install.packages("Seurat")
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))
install.packages("magrittr")
install.packages("parallel")
install.packages("Matrix")
install.packages("fields")
install.packages("KernSmooth")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')


if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")



# Load Library
library(Seurat)
library(tidyverse)
library(DoubletFinder)
library(ggplot2)
library(KernSmooth)
library(fields)


#### Load Data Into R ####

# Control Samples 

C1_mtx <- ReadMtx(mtx = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_02_GEX_FL/filtered_feature_bc_matrix/matrix.mtx.gz",
                  features = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_02_GEX_FL/filtered_feature_bc_matrix/features.tsv.gz",
                  cells = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_02_GEX_FL/filtered_feature_bc_matrix/barcodes.tsv.gz")

C1_seurat_obj <- CreateSeuratObject(counts = C1_mtx, project = "Control", min.cells = 3, min.features = 200)

C2_mtx <- ReadMtx(mtx = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_04_GEX_FL/filtered_feature_bc_matrix/matrix.mtx.gz",
                  features = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_04_GEX_FL/filtered_feature_bc_matrix/features.tsv.gz",
                  cells = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_04_GEX_FL/filtered_feature_bc_matrix/barcodes.tsv.gz")

C2_seurat_obj <- CreateSeuratObject(counts = C2_mtx, project = "Control", min.cells = 3, min.features = 200)

#C3_mtx <- ReadMtx(mtx = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_06_GEX_FL/filtered_feature_bc_matrix/matrix.mtx.gz",
                  #features = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_06_GEX_FL/filtered_feature_bc_matrix/features.tsv.gz",
                  #cells = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_06_GEX_FL/filtered_feature_bc_matrix/barcodes.tsv.gz")

#C3_seurat_obj <- CreateSeuratObject(counts = C3_mtx, project = "PAH", min.cells = 3, min.features = 200)

C4_mtx <- ReadMtx(mtx = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_08_GEX_FL/filtered_feature_bc_matrix/matrix.mtx.gz",
                  features = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_08_GEX_FL/filtered_feature_bc_matrix/features.tsv.gz",
                  cells = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_08_GEX_FL/filtered_feature_bc_matrix/barcodes.tsv.gz")

C4_seurat_obj <- CreateSeuratObject(counts = C4_mtx, project = "Control", min.cells = 3, min.features = 200)

C5_mtx <- ReadMtx(mtx = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_10_GEX_FL/filtered_feature_bc_matrix/matrix.mtx.gz",
                  features = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_10_GEX_FL/filtered_feature_bc_matrix/features.tsv.gz",
                  cells = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_10_GEX_FL/filtered_feature_bc_matrix/barcodes.tsv.gz")

C5_seurat_obj <- CreateSeuratObject(counts = C5_mtx, project = "Control", min.cells = 3, min.features = 200)

#PAH Samples 

P1_mtx <- ReadMtx(mtx = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_01_GEX_FL/filtered_feature_bc_matrix/matrix.mtx.gz",
                  features = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_01_GEX_FL/filtered_feature_bc_matrix/features.tsv.gz",
                  cells = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_01_GEX_FL/filtered_feature_bc_matrix/barcodes.tsv.gz")

P1_seurat_obj <- CreateSeuratObject(counts = P1_mtx, project = "PAH", min.cells = 3, min.features = 200)

P2_mtx <- ReadMtx(mtx = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_03_GEX_FL/filtered_feature_bc_matrix/matrix.mtx.gz",
                  features = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_03_GEX_FL/filtered_feature_bc_matrix/features.tsv.gz",
                  cells = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_03_GEX_FL/filtered_feature_bc_matrix/barcodes.tsv.gz")

P2_seurat_obj <- CreateSeuratObject(counts = P2_mtx, project = "PAH", min.cells = 3, min.features = 200)

P3_mtx <- ReadMtx(mtx = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_05_GEX_FL/filtered_feature_bc_matrix/matrix.mtx.gz",
                  features = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_05_GEX_FL/filtered_feature_bc_matrix/features.tsv.gz",
                  cells = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_05_GEX_FL/filtered_feature_bc_matrix/barcodes.tsv.gz")

P3_seurat_obj <- CreateSeuratObject(counts = P3_mtx, project = "PAH", min.cells = 3, min.features = 200)

P4_mtx <- ReadMtx(mtx = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_07_GEX_FL/filtered_feature_bc_matrix/matrix.mtx.gz",
                  features = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_07_GEX_FL/filtered_feature_bc_matrix/features.tsv.gz",
                  cells = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_07_GEX_FL/filtered_feature_bc_matrix/barcodes.tsv.gz")

P4_seurat_obj <- CreateSeuratObject(counts = P4_mtx, project = "PAH", min.cells = 3, min.features = 200)

P5_mtx <- ReadMtx(mtx = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_09_GEX_FL/filtered_feature_bc_matrix/matrix.mtx.gz",
                  features = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_09_GEX_FL/filtered_feature_bc_matrix/features.tsv.gz",
                  cells = "/home/prin0088/data_release/umgc/data_delivery/2024-q4/20241204_AV224201_Aviti_Run_140_12_04_24/Prins_Project_015_GEX/Analysis/cellranger_MB241112_09_GEX_FL/filtered_feature_bc_matrix/barcodes.tsv.gz")

P5_seurat_obj <- CreateSeuratObject(counts = P5_mtx, project = "PAH", min.cells = 3, min.features = 200)


#### Pre-Processing ####

# Quality Control - Visualization of Data

#Control 

C1_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(C1_seurat_obj, pattern = "^MT")
C2_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(C2_seurat_obj, pattern = "^MT")
C4_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(C4_seurat_obj, pattern = "^MT")
C5_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(C5_seurat_obj, pattern = "^MT")

#PAH Samples# 

P1_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(P1_seurat_obj, pattern = "^MT")
P2_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(P2_seurat_obj, pattern = "^MT")
P3_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(P3_seurat_obj, pattern = "^MT")
P4_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(P4_seurat_obj, pattern = "^MT")
P5_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(P5_seurat_obj, pattern = "^MT")

##Normalization

#Control Samples

#Control 1

C1_seurat_filt <- subset(C1_seurat_obj, subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 8000 &
                           nCount_RNA < 750000 &
                           percent.mt < 5)
C1_seurat_filt <- NormalizeData(C1_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

#Control 2

C2_seurat_filt <- subset(C2_seurat_obj, subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 7000 &
                           nCount_RNA < 55000 &
                           percent.mt < 5)
C2_seurat_filt <- NormalizeData(C2_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

#Control 4

C4_seurat_filt <- subset(C4_seurat_obj, subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 8500 &
                           nCount_RNA < 80000 &
                           percent.mt < 5)
C4_seurat_filt <- NormalizeData(C4_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

#Control 5

C5_seurat_filt <- subset(C5_seurat_obj, subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 7500 &
                           nCount_RNA < 55000 &
                           percent.mt < 5)
C5_seurat_filt <- NormalizeData(C5_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

#PAH Samples

#PAH 1

P1_seurat_filt <- subset(P1_seurat_obj, subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 8000 &
                           nCount_RNA < 750000 &
                           percent.mt < 5)
P1_seurat_filt <- NormalizeData(P1_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

#PAH 2

P2_seurat_filt <- subset(P2_seurat_obj, subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 7000 &
                           nCount_RNA < 55000 &
                           percent.mt < 5)
P2_seurat_filt <- NormalizeData(P2_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

#PAH 3

P3_seurat_filt <- subset(P3_seurat_obj, subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 6500 &
                           nCount_RNA < 42000 &
                           percent.mt < 5)
P3_seurat_filt <- NormalizeData(P3_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

#PAH 4

P4_seurat_filt <- subset(P4_seurat_obj, subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 8500 &
                           nCount_RNA < 80000 &
                           percent.mt < 5)
P4_seurat_filt <- NormalizeData(P4_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

#PAH 5

P5_seurat_filt <- subset(P5_seurat_obj, subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 7500 &
                           nCount_RNA < 55000 &
                           percent.mt < 5)
P5_seurat_filt <- NormalizeData(P5_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

##Find Variable Features##

#Control Samples#

#Control 1

C1_seurat_filt <- FindVariableFeatures(C1_seurat_filt, selection.method = "vst", nfeatures = 2000)

#Control 2

C2_seurat_filt <- FindVariableFeatures(C2_seurat_filt, selection.method = "vst", nfeatures = 2000)

#Control 3

#C3_seurat_filt <- FindVariableFeatures(C3_seurat_filt, selection.method = "vst", nfeatures = 2000)

#Control 4

C4_seurat_filt <- FindVariableFeatures(C4_seurat_filt, selection.method = "vst", nfeatures = 2000)

#Control 5

C5_seurat_filt <- FindVariableFeatures(C5_seurat_filt, selection.method = "vst", nfeatures = 2000)

#PAH Samples

#PAH 1

P1_seurat_filt <- FindVariableFeatures(P1_seurat_filt, selection.method = "vst", nfeatures = 2000)

#PAH 2

P2_seurat_filt <- FindVariableFeatures(P2_seurat_filt, selection.method = "vst", nfeatures = 2000)

#PAH 3

P3_seurat_filt <- FindVariableFeatures(P3_seurat_filt, selection.method = "vst", nfeatures = 2000)

#PAH 4

P4_seurat_filt <- FindVariableFeatures(P4_seurat_filt, selection.method = "vst", nfeatures = 2000)

#PAH 5

P5_seurat_filt <- FindVariableFeatures(P5_seurat_filt, selection.method = "vst", nfeatures = 2000)

## Scale Data##

#Control Samples#

C1_allgenes <- rownames(C1_seurat_filt)
C1_seurat_filt <- ScaleData(C1_seurat_filt, features = C1_allgenes)

C2_allgenes <- rownames(C2_seurat_filt)
C2_seurat_filt <- ScaleData(C2_seurat_filt, features = C2_allgenes)

#C3_allgenes <- rownames(C3_seurat_filt)
#C3_seurat_filt <- ScaleData(C3_seurat_filt, features = C3_allgenes)

C4_allgenes <- rownames(C4_seurat_filt)
C4_seurat_filt <- ScaleData(C4_seurat_filt, features = C4_allgenes)

C5_allgenes <- rownames(C5_seurat_filt)
C5_seurat_filt <- ScaleData(C5_seurat_filt, features = C5_allgenes)

#PAH Samples#

P1_allgenes <- rownames(P1_seurat_filt)
P1_seurat_filt <- ScaleData(P1_seurat_filt, features = P1_allgenes)

P2_allgenes <- rownames(P2_seurat_filt)
P2_seurat_filt <- ScaleData(P2_seurat_filt, features = P2_allgenes)

P3_allgenes <- rownames(P3_seurat_filt)
P3_seurat_filt <- ScaleData(P3_seurat_filt, features = P3_allgenes)

P4_allgenes <- rownames(P4_seurat_filt)
P4_seurat_filt <- ScaleData(P4_seurat_filt, features = P4_allgenes)

P5_allgenes <- rownames(P5_seurat_filt)
P5_seurat_filt <- ScaleData(P5_seurat_filt, features = P5_allgenes)

#### Principal Component Analysis (PCA) - Linear Dimension Reduction ####

#Control Samples#

#Control 1

C1_seurat_filt <- RunPCA(C1_seurat_filt, features = VariableFeatures(object = C1_seurat_filt))
#print(C1_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5) 
#DimHeatmap(C1_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
#ElbowPlot(C1_seurat_filt)+
  #ggtitle("C1 Elbow Plot")

#Control 2

C2_seurat_filt <- RunPCA(C2_seurat_filt, features = VariableFeatures(object = C2_seurat_filt))
#print(C2_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5)
#DimHeatmap(C2_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
#ElbowPlot(C2_seurat_filt)+
  #ggtitle("C2 Elbow Plot")

#Control 3

#C3_seurat_filt <- RunPCA(C3_seurat_filt, features = VariableFeatures(object = C3_seurat_filt))
#print(C3_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5)
#DimHeatmap(C3_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
#ElbowPlot(C3_seurat_filt)+
  #ggtitle("C3 Elbow Plot")

#Control 4

C4_seurat_filt <- RunPCA(C4_seurat_filt, features = VariableFeatures(object = C4_seurat_filt))
#print(C4_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5)
#DimHeatmap(C4_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
#ElbowPlot(C4_seurat_filt)+
  #ggtitle("C4 Elbow Plot")

#Control 5

C5_seurat_filt <- RunPCA(C5_seurat_filt, features = VariableFeatures(object = C5_seurat_filt))
#print(C5_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5)
#DimHeatmap(C5_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
#ElbowPlot(C5_seurat_filt)+
  #ggtitle("C5 Elbow Plot")

#PAH Samples#

#PAH 1

P1_seurat_filt <- RunPCA(P1_seurat_filt, features = VariableFeatures(object = P1_seurat_filt))
#print(P1_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5)
#DimHeatmap(P1_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
#ElbowPlot(P1_seurat_filt)+
  #ggtitle("P1 Elbow Plot")

#PAH 2

P2_seurat_filt <- RunPCA(P2_seurat_filt, features = VariableFeatures(object = P2_seurat_filt))
#print(P2_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5)
#DimHeatmap(P2_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
#ElbowPlot(P2_seurat_filt)+
  #ggtitle("P2 Elbow Plot")

#PAH 3

P3_seurat_filt <- RunPCA(P3_seurat_filt, features = VariableFeatures(object = P3_seurat_filt))
#print(P3_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5)
#DimHeatmap(P3_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
#ElbowPlot(P3_seurat_filt)+
  #ggtitle("P3 Elbow Plot")

#PAH 4

P4_seurat_filt <- RunPCA(P4_seurat_filt, features = VariableFeatures(object = P4_seurat_filt))
#print(P4_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5)
#DimHeatmap(P4_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
#ElbowPlot(P4_seurat_filt)+
  #ggtitle("P4 Elbow Plot")

#PAH 5

P5_seurat_filt <- RunPCA(P5_seurat_filt, features = VariableFeatures(object = P5_seurat_filt))
#print(P5_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5)
#DimHeatmap(P5_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
#ElbowPlot(P5_seurat_filt)+
  #ggtitle("P5 Elbow Plot")

#### Clustering ####

#Control Samples#

#Control 1

C1_seurat_filt <- FindNeighbors(C1_seurat_filt, dims = 1:18)
C1_seurat_filt <- FindClusters(object = C1_seurat_filt)
C1_seurat_filt <- RunUMAP(C1_seurat_filt, dims = 1:18)

#DimPlot(C1_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

#To view what to group by, do view(C1_seurat_filt @meta.data); should just do seurat_clusters because RNA_snn looked the same w/ clusters

#Control 2

C2_seurat_filt <- FindNeighbors(C2_seurat_filt, dims = 1:18)
C2_seurat_filt <- FindClusters(object = C2_seurat_filt)
C2_seurat_filt <- RunUMAP(C2_seurat_filt, dims = 1:18)

#DimPlot(C2_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

#Control 3

#C3_seurat_filt <- FindNeighbors(C3_seurat_filt, dims = 1:19)
#C3_seurat_filt <- FindClusters(object = C3_seurat_filt)
#C3_seurat_filt <- RunUMAP(C3_seurat_filt, dims = 1:19)

#DimPlot(C3_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

#Control 4

C4_seurat_filt <- FindNeighbors(C4_seurat_filt, dims = 1:19)
C4_seurat_filt <- FindClusters(object = C4_seurat_filt)
C4_seurat_filt <- RunUMAP(C4_seurat_filt, dims = 1:19)

#DimPlot(C4_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

#Control 5

C5_seurat_filt <- FindNeighbors(C5_seurat_filt, dims = 1:18)
C5_seurat_filt <- FindClusters(object = C5_seurat_filt)
C5_seurat_filt <- RunUMAP(C5_seurat_filt, dims = 1:18)

#DimPlot(C5_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

#PAH Samples#

#PAH 1

P1_seurat_filt <- FindNeighbors(P1_seurat_filt, dims = 1:18)
P1_seurat_filt <- FindClusters(object = P1_seurat_filt)
P1_seurat_filt <- RunUMAP(P1_seurat_filt, dims = 1:18)

#DimPlot(P1_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

#PAH 2

P2_seurat_filt <- FindNeighbors(P2_seurat_filt, dims = 1:19)
P2_seurat_filt <- FindClusters(object = P2_seurat_filt)
P2_seurat_filt <- RunUMAP(P2_seurat_filt, dims = 1:19)

#DimPlot(P2_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

#PAH 3

P3_seurat_filt <- FindNeighbors(P3_seurat_filt, dims = 1:18)
P3_seurat_filt <- FindClusters(object = P3_seurat_filt)
P3_seurat_filt <- RunUMAP(P3_seurat_filt, dims = 1:18)

#DimPlot(P3_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

#PAH 4

P4_seurat_filt <- FindNeighbors(P4_seurat_filt, dims = 1:17)
P4_seurat_filt <- FindClusters(object = P4_seurat_filt)
P4_seurat_filt <- RunUMAP(P4_seurat_filt, dims = 1:17)

#Changed Dimensions from [1:18] to [1:17] Due to Large Number of Clusters Initially

#DimPlot(P4_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

#PAH 5

P5_seurat_filt <- FindNeighbors(P5_seurat_filt, dims = 1:15)
P5_seurat_filt <- FindClusters(object = P5_seurat_filt)
P5_seurat_filt <- RunUMAP(P5_seurat_filt, dims = 1:15)

#Changed Dimensions from [1:18] to [1:15] Due to Large Number of Clusters Initially

#DimPlot(P5_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

#Combining Samples Into Groups#

#Control Samples#

#Contotal_matrix <- merge(C1_seurat_filt, 
                         #y = list(C2_seurat_filt, C3_seurat_filt, C4_seurat_filt, C5_seurat_filt),
                         #add.cell.ids = c("C1", "C2", "C3", "C4", "C5"), 
                         #project = "Control")

#Con_seurat_filt <- FindVariableFeatures(Contotal_matrix, selection.method = "vst", nfeatures = 2000)
#Con_allgenes <- rownames(Con_seurat_filt)
#Con_seurat_filt <- ScaleData(Con_seurat_filt, features = Con_allgenes)
#Con_seurat_filt <- RunPCA(Con_seurat_filt, features = VariableFeatures(object = Con_seurat_filt))


#print(Con_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5)
#DimHeatmap(Con_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
#ElbowPlot(Con_seurat_filt)+
#ggtitle("Control All Samples Elbow Plot")

#Con_seurat_filt <- FindNeighbors(Con_seurat_filt, dims = 1:8) #Was Initially [1:17] - 34 Clusters, 
#Con_seurat_filt <- FindClusters(Con_seurat_filt)
#Con_seurat_filt <- RunUMAP(Con_seurat_filt, dims = 1:8)

#DimPlot(Con_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

#PAH Samples#

#PAHtotal_matrix <- merge(P1_seurat_filt, 
                         #y = list(P2_seurat_filt, P3_seurat_filt, P4_seurat_filt, P5_seurat_filt),
                         #add.cell.ids = c("P1", "P2", "P3", "P4", "P5"), 
                         #project = "PAH")

#PAH_seurat_filt <- FindVariableFeatures(PAHtotal_matrix, selection.method = "vst", nfeatures = 2000)
#PAH_allgenes <- rownames(PAH_seurat_filt)
#PAH_seurat_filt <- ScaleData(PAH_seurat_filt, features = PAH_allgenes)
#PAH_seurat_filt <- RunPCA(PAH_seurat_filt, features = VariableFeatures(object = PAH_seurat_filt))


#print(PAH_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5)
#DimHeatmap(PAH_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
#ElbowPlot(PAH_seurat_filt)+
#ggtitle("PAH All Samples Elbow Plot")

#PAH_seurat_filt <- FindNeighbors(PAH_seurat_filt, dims = 1:11) 
#PAH_seurat_filt <- FindClusters(PAH_seurat_filt)
#PAH_seurat_filt <- RunUMAP(PAH_seurat_filt, dims = 1:11)

#DimPlot(PAH_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

##Doublet Finder##

#C1 ---
C1_sweep_res_list <- paramSweep(C1_seurat_filt, PCs = 1:18, sct = FALSE)
C1_sweep_stats <- summarizeSweep(C1_sweep_res_list, GT = FALSE)
C1_find_pK <- find.pK(C1_sweep_stats)

#Visualization
ggplot(C1_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
C1_pK <- C1_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
C1_pK <- as.numeric(as.character(C1_pK[[1]]))

## Homotypic Doublet Proportion Estimate
C1_annotations <- C1_seurat_filt@meta.data$seurat_clusters
C1_homotypic_prop <- modelHomotypic(C1_annotations)
C1_nExp.poi <- round(0.136*nrow(C1_seurat_filt@meta.data)) # 17,000 cells recovered, %0.8 for every 1000 cells recovered (i.e 25,000 --> calculation for red number, 25 * 0.8 = 20% --> 0.2; really was 24,600 but round to nearest whole number)
C1_nEXP.poi.adj <- round(C1_nExp.poi *(1-C1_homotypic_prop))

C1_seurat_filt <- doubletFinder(C1_seurat_filt, PCs = 1:18, pN = 0.25, pK = C1_pK, nExp = C1_nEXP.poi.adj,
                                  reuse.pANN = FALSE, sct = FALSE)
DimPlot(C1_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_2486") #Change the green name by finding the DF.classifcations file name by typing view(C1_seurat_filt@meta.data)

C1_seurat_filt <- subset(C1_seurat_filt, subset = DF.classifications_0.25_0.005_2486 == "Singlet") #Change this DF.classifcations as well! 
saveRDS(C1_seurat_filt, file = "C1_doubremoved")

#C2 ---
C2_sweep_res_list <- paramSweep(C2_seurat_filt, PCs = 1:18, sct = FALSE)
C2_sweep_stats <- summarizeSweep(C2_sweep_res_list, GT = FALSE)
C2_find_pK <- find.pK(C2_sweep_stats)

#Visualization
ggplot(C2_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
C2_pK <- C2_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
C2_pK <- as.numeric(as.character(C2_pK[[1]]))

## Homotypic Doublet Proportion Estimate
C2_annotations <- C2_seurat_filt@meta.data$seurat_clusters
C2_homotypic_prop <- modelHomotypic(C2_annotations)
C2_nExp.poi <- round(0.08*nrow(C2_seurat_filt@meta.data)) # 10,000 cells recovered, %0.8 for every 1000 cells recovered
C2_nEXP.poi.adj <- round(C2_nExp.poi *(1-C2_homotypic_prop))

C2_seurat_filt <- doubletFinder(C2_seurat_filt, PCs = 1:18, pN = 0.25, pK = C2_pK, nExp = C2_nEXP.poi.adj,
                                reuse.pANN = FALSE, sct = FALSE)
DimPlot(C2_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_702")

C2_seurat_filt <- subset(C2_seurat_filt, subset = DF.classifications_0.25_0.005_702 == "Singlet")
saveRDS(C2_seurat_filt, file = "C2_doubremoved")

#C4 ---
C4_sweep_res_list <- paramSweep(C4_seurat_filt, PCs = 1:19, sct = FALSE)
C4_sweep_stats <- summarizeSweep(C4_sweep_res_list, GT = FALSE)
C4_find_pK <- find.pK(C4_sweep_stats)

#Visualization
ggplot(C4_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
C4_pK <- C4_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
C4_pK <- as.numeric(as.character(C4_pK[[1]]))

## Homotypic Doublet Proportion Estimate
C4_annotations <- C4_seurat_filt@meta.data$seurat_clusters
C4_homotypic_prop <- modelHomotypic(C4_annotations)
C4_nExp.poi <- round(0.136*nrow(C4_seurat_filt@meta.data)) # 17,000 cells recovered, %0.8 for every 1000 cells recovered
C4_nEXP.poi.adj <- round(C4_nExp.poi *(1-C4_homotypic_prop))

C4_seurat_filt <- doubletFinder(C4_seurat_filt, PCs = 1:19, pN = 0.25, pK = C4_pK, nExp = C4_nEXP.poi.adj,
                                reuse.pANN = FALSE, sct = FALSE)
DimPlot(C4_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_2104")

C4_seurat_filt <- subset(C4_seurat_filt, subset = DF.classifications_0.25_0.005_2104 == "Singlet")
saveRDS(C4_seurat_filt, file = "C4_doubremoved")

#C5 ---
C5_sweep_res_list <- paramSweep(C5_seurat_filt, PCs = 1:18, sct = FALSE)
C5_sweep_stats <- summarizeSweep(C5_sweep_res_list, GT = FALSE)
C5_find_pK <- find.pK(C5_sweep_stats)

#Visualization
ggplot(C5_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
C5_pK <- C5_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
C5_pK <- as.numeric(as.character(C5_pK[[1]]))

## Homotypic Doublet Proportion Estimate
C5_annotations <- C5_seurat_filt@meta.data$seurat_clusters
C5_homotypic_prop <- modelHomotypic(C5_annotations)
C5_nExp.poi <- round(0.104*nrow(C5_seurat_filt@meta.data)) # 24,000 cells recovered, %0.8 for every 1000 cells recovered
C5_nEXP.poi.adj <- round(C5_nExp.poi *(1-C5_homotypic_prop))

C5_seurat_filt <- doubletFinder(C5_seurat_filt, PCs = 1:18, pN = 0.25, pK = C5_pK, nExp = C5_nEXP.poi.adj,
                                reuse.pANN = FALSE, sct = FALSE)
DimPlot(C5_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.23_1220")

C5_seurat_filt <- subset(C5_seurat_filt, subset = DF.classifications_0.25_0.23_1220 == "Singlet")
saveRDS(C5_seurat_filt, file = "C5_doubremoved")

##PAH Samples

#P1 ---
P1_sweep_res_list <- paramSweep(P1_seurat_filt, PCs = 1:18, sct = FALSE)
P1_sweep_stats <- summarizeSweep(P1_sweep_res_list, GT = FALSE)
P1_find_pK <- find.pK(P1_sweep_stats)

#Visualization
ggplot(P1_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
P1_pK <- P1_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
P1_pK <- as.numeric(as.character(P1_pK[[1]]))

## Homotypic Doublet Proportion Estimate
P1_annotations <- P1_seurat_filt@meta.data$seurat_clusters
P1_homotypic_prop <- modelHomotypic(P1_annotations)
P1_nExp.poi <- round(0.088*nrow(P1_seurat_filt@meta.data)) # 11,000 cells recovered, %0.8 for every 1000 cells recovered
P1_nEXP.poi.adj <- round(P1_nExp.poi *(1-P1_homotypic_prop))

P1_seurat_filt <- doubletFinder(P1_seurat_filt, PCs = 1:18, pN = 0.25, pK = P1_pK, nExp = P1_nEXP.poi.adj,
                                reuse.pANN = FALSE, sct = FALSE)
DimPlot(P1_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_876")

P1_seurat_filt <- subset(P1_seurat_filt, subset = DF.classifications_0.25_0.005_876 == "Singlet")
saveRDS(P1_seurat_filt, file = "P1_doubremoved")

#P2 ---
P2_sweep_res_list <- paramSweep(P2_seurat_filt, PCs = 1:19, sct = FALSE)
P2_sweep_stats <- summarizeSweep(P2_sweep_res_list, GT = FALSE)
P2_find_pK <- find.pK(P2_sweep_stats)

#Visualization
ggplot(P2_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
P2_pK <- P2_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
P2_pK <- as.numeric(as.character(P2_pK[[1]]))

## Homotypic Doublet Proportion Estimate
P2_annotations <- P2_seurat_filt@meta.data$seurat_clusters
P2_homotypic_prop <- modelHomotypic(P2_annotations)
P2_nExp.poi <- round(0.128*nrow(P2_seurat_filt@meta.data)) # 16,000 cells recovered, %0.8 for every 1000 cells recovered
P2_nEXP.poi.adj <- round(P2_nExp.poi *(1-P2_homotypic_prop))

P2_seurat_filt <- doubletFinder(P2_seurat_filt, PCs = 1:19, pN = 0.25, pK = P2_pK, nExp = P2_nEXP.poi.adj,
                                reuse.pANN = FALSE, sct = FALSE)
DimPlot(P2_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_1823")

P2_seurat_filt <- subset(P2_seurat_filt, subset = DF.classifications_0.25_0.005_1823 == "Singlet")
saveRDS(P2_seurat_filt, file = "P2_doubremoved")

#P3 ---
P3_sweep_res_list <- paramSweep(P3_seurat_filt, PCs = 1:18, sct = FALSE)
P3_sweep_stats <- summarizeSweep(P3_sweep_res_list, GT = FALSE)
P3_find_pK <- find.pK(P3_sweep_stats)

#Visualization
ggplot(P3_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
P3_pK <- P3_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
P3_pK <- as.numeric(as.character(P3_pK[[1]]))

## Homotypic Doublet Proportion Estimate
P3_annotations <- P3_seurat_filt@meta.data$seurat_clusters
P3_homotypic_prop <- modelHomotypic(P3_annotations)
P3_nExp.poi <- round(0.144*nrow(P3_seurat_filt@meta.data)) # 18,000 cells recovered, %0.8 for every 1000 cells recovered
P3_nEXP.poi.adj <- round(P3_nExp.poi *(1-P3_homotypic_prop))

P3_seurat_filt <- doubletFinder(P3_seurat_filt, PCs = 1:18, pN = 0.25, pK = P3_pK, nExp = P3_nEXP.poi.adj,
                                reuse.pANN = FALSE, sct = FALSE)
DimPlot(P3_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_2332")

P3_seurat_filt <- subset(P3_seurat_filt, subset = DF.classifications_0.25_0.005_2332 == "Singlet")
saveRDS(P3_seurat_filt, file = "P3_doubremoved")

#P4 ---
P4_sweep_res_list <- paramSweep(P4_seurat_filt, PCs = 1:17, sct = FALSE)
P4_sweep_stats <- summarizeSweep(P4_sweep_res_list, GT = FALSE)
P4_find_pK <- find.pK(P4_sweep_stats)

#Visualization
ggplot(P4_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
P4_pK <- P4_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
P4_pK <- as.numeric(as.character(P4_pK[[1]]))

## Homotypic Doublet Proportion Estimate
P4_annotations <- P4_seurat_filt@meta.data$seurat_clusters
P4_homotypic_prop <- modelHomotypic(P4_annotations)
P4_nExp.poi <- round(0.136*nrow(P4_seurat_filt@meta.data)) # 24,000 cells recovered, %0.8 for every 1000 cells recovered
P4_nEXP.poi.adj <- round(P4_nExp.poi *(1-P4_homotypic_prop))

P4_seurat_filt <- doubletFinder(P4_seurat_filt, PCs = 1:17, pN = 0.25, pK = P4_pK, nExp = P4_nEXP.poi.adj,
                                reuse.pANN = FALSE, sct = FALSE)
DimPlot(P4_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_2071")

P4_seurat_filt <- subset(P4_seurat_filt, subset = DF.classifications_0.25_0.005_2071 == "Singlet")
saveRDS(P4_seurat_filt, file = "P4_doubremoved")

#P5 ---
P5_sweep_res_list <- paramSweep(P5_seurat_filt, PCs = 1:15, sct = FALSE)
P5_sweep_stats <- summarizeSweep(P5_sweep_res_list, GT = FALSE)
P5_find_pK <- find.pK(P5_sweep_stats)

#Visualization
ggplot(P5_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
P5_pK <- P5_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
P5_pK <- as.numeric(as.character(P5_pK[[1]]))

## Homotypic Doublet Proportion Estimate
P5_annotations <- P5_seurat_filt@meta.data$seurat_clusters
P5_homotypic_prop <- modelHomotypic(P5_annotations)
P5_nExp.poi <- round(0.168*nrow(P5_seurat_filt@meta.data)) # 21,000 cells recovered, %0.8 for every 1000 cells recovered
P5_nEXP.poi.adj <- round(P5_nExp.poi *(1-P5_homotypic_prop))

P5_seurat_filt <- doubletFinder(P5_seurat_filt, PCs = 1:15, pN = 0.25, pK = P5_pK, nExp = P5_nEXP.poi.adj,
                                reuse.pANN = FALSE, sct = FALSE)
DimPlot(P5_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_3283")

P5_seurat_filt <- subset(P5_seurat_filt, subset = DF.classifications_0.25_0.005_3283 == "Singlet")
saveRDS(P5_seurat_filt, file = "P5_doubremoved")


##### SAVE ALL DATA TO NEW File (TO integrate)

total_matrix <- merge(C1_seurat_filt, y = list(C2_seurat_filt, C4_seurat_filt, C5_seurat_filt, 
                                                    P1_seurat_filt, P2_seurat_filt, P3_seurat_filt, P4_seurat_filt, 
                                                    P5_seurat_filt), 
                         add.cell.ids = c("C1","C2", "C4", "C5", "P1", "P2", "P3", 
                                          "P4", "P5"), project = ("All.Samples"))
saveRDS(total_matrix, file = "total_matrix")

