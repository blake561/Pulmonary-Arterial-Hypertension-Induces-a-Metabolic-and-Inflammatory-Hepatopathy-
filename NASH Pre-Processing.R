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

library(Matrix)
library(Seurat)
library(tidyverse)
library(DoubletFinder)
library(ggplot2)
library(KernSmooth)
library(fields)


#### Load Data Into R ####

##NASH Livers

#Non-Rerun Samples

N1_mtx <- ReadMtx(mtx = "/users/1/blake561/NASH Files/GSM6556466_humanNASH_9_matrix.mtx.gz",
                  features = "/users/1/blake561/NASH Files/GSM6556466_humanNASH_9_features.tsv.gz",
                  cells = "/users/1/blake561/NASH Files/GSM6556466_humanNASH_9_barcodes.tsv.gz")

N1_seurat_obj <- CreateSeuratObject(counts = N1_mtx, project = "NASH", min.cells = 3, min.features = 200)

N2_mtx <- ReadMtx(mtx = "/users/1/blake561/NASH Files/GSM6556465_humanNASH_8_matrix.mtx.gz",
                  features = "/users/1/blake561/NASH Files/GSM6556465_humanNASH_8_features.tsv.gz",
                  cells = "/users/1/blake561/NASH Files/GSM6556465_humanNASH_8_barcodes.tsv.gz")

N2_seurat_obj <- CreateSeuratObject(counts = N2_mtx, project = "NASH", min.cells = 3, min.features = 200)

N3_mtx <- ReadMtx(mtx = "/users/1/blake561/NASH Files/GSM6556464_humanNASH_7_matrix.mtx.gz",
                  features = "/users/1/blake561/NASH Files/GSM6556464_humanNASH_7_features.tsv.gz",
                  cells = "/users/1/blake561/NASH Files/GSM6556464_humanNASH_7_barcodes.tsv.gz")

N3_seurat_obj <- CreateSeuratObject(counts = N3_mtx, project = "NASH", min.cells = 3, min.features = 200)

N4_mtx <- ReadMtx(mtx = "/users/1/blake561/NASH Files/GSM6556463_humanNASH_6_matrix.mtx.gz",
                  features = "/users/1/blake561/NASH Files/GSM6556463_humanNASH_6_features.tsv.gz",
                  cells = "/users/1/blake561/NASH Files/GSM6556463_humanNASH_6_barcodes.tsv.gz")

N4_seurat_obj <- CreateSeuratObject(counts = N4_mtx, project = "NASH", min.cells = 3, min.features = 200)

#NASH Controls

CN1_mtx <- ReadMtx(mtx = "/users/1/blake561/NASH Files/GSM6556449_humanCTRL_1_matrix.mtx.gz",
                   features = "/users/1/blake561/NASH Files/GSM6556449_humanCTRL_1_features.tsv.gz",
                   cells = "/users/1/blake561/NASH Files/GSM6556449_humanCTRL_1_barcodes.tsv.gz")

CN1_seurat_obj <- CreateSeuratObject(counts = CN1_mtx, project = "NASH Control", min.cells = 3, min.features = 200)

CN2_mtx <- ReadMtx(mtx = "/users/1/blake561/NASH Files/GSM6556450_humanCTRL_2_matrix.mtx.gz",
                   features = "/users/1/blake561/NASH Files/GSM6556450_humanCTRL_2_features.tsv.gz",
                   cells = "/users/1/blake561/NASH Files/GSM6556450_humanCTRL_2_barcodes.tsv.gz")

CN2_seurat_obj <- CreateSeuratObject(counts = CN2_mtx, project = "NASH Control", min.cells = 3, min.features = 200)

CN3_mtx <- ReadMtx(mtx = "/users/1/blake561/NASH Files/GSM6556451_humanCTRL_3_matrix.mtx.gz",
                   features = "/users/1/blake561/NASH Files/GSM6556451_humanCTRL_3_features.tsv.gz",
                   cells = "/users/1/blake561/NASH Files/GSM6556451_humanCTRL_3_barcodes.tsv.gz")

CN3_seurat_obj <- CreateSeuratObject(counts = CN3_mtx, project = "NASH Control", min.cells = 3, min.features = 200)

#### Pre-Processing ####

# Quality Control - Visualization of Data

#NASH 

N1_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(N1_seurat_obj, pattern = "^MT")
N2_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(N2_seurat_obj, pattern = "^MT")
N3_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(N3_seurat_obj, pattern = "^MT")
N4_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(N4_seurat_obj, pattern = "^MT")

#NASH Controls

CN1_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(CN1_seurat_obj, pattern = "^MT")
CN2_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(CN2_seurat_obj, pattern = "^MT")
CN3_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(CN3_seurat_obj, pattern = "^MT")

##Normalization

#NASH 

#N1

N1_seurat_filt <- subset(N1_seurat_obj, subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 8000 &
                           percent.mt < 10) #NASH samples had more mitochondrial contamination relative to Fontan samples and PAH samples
N1_seurat_filt <- NormalizeData(N1_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

VlnPlot(N1_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#N2

N2_seurat_filt <- subset(N2_seurat_obj, subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 8000 &
                           percent.mt < 10)
N2_seurat_filt <- NormalizeData(N2_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

VlnPlot(N2_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#N3

N3_seurat_filt <- subset(N3_seurat_obj, subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 8000 &
                           percent.mt < 10)
N3_seurat_filt <- NormalizeData(N3_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

VlnPlot(N3_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#N4

N4_seurat_filt <- subset(N4_seurat_obj, subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 8000 &
                           percent.mt < 10)
N4_seurat_filt <- NormalizeData(N4_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

VlnPlot(N4_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#NASH 

#N1

CN1_seurat_filt <- subset(CN1_seurat_obj, subset = nFeature_RNA > 200 & 
                            nFeature_RNA < 8000 &
                            percent.mt < 10) #more mitochondrial DNA compared to other samples 
CN1_seurat_filt <- NormalizeData(CN1_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

VlnPlot(CN1_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#N2

CN2_seurat_filt <- subset(CN2_seurat_obj, subset = nFeature_RNA > 200 & 
                            nFeature_RNA < 8000 &
                            percent.mt < 10)
CN2_seurat_filt <- NormalizeData(CN2_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

VlnPlot(CN2_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#N3

CN3_seurat_filt <- subset(CN3_seurat_obj, subset = nFeature_RNA > 200 & 
                            nFeature_RNA < 8000 &
                            percent.mt < 10)
CN3_seurat_filt <- NormalizeData(CN3_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

VlnPlot(CN3_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

##Find Variable Features##

#NASH

N1_seurat_filt <- FindVariableFeatures(N1_seurat_filt, selection.method = "vst", nfeatures = 2000)
N2_seurat_filt <- FindVariableFeatures(N2_seurat_filt, selection.method = "vst", nfeatures = 2000)
N3_seurat_filt <- FindVariableFeatures(N3_seurat_filt, selection.method = "vst", nfeatures = 2000)
N4_seurat_filt <- FindVariableFeatures(N4_seurat_filt, selection.method = "vst", nfeatures = 2000)

#NASH Controls 

CN1_seurat_filt <- FindVariableFeatures(CN1_seurat_filt, selection.method = "vst", nfeatures = 2000)
CN2_seurat_filt <- FindVariableFeatures(CN2_seurat_filt, selection.method = "vst", nfeatures = 2000)
CN3_seurat_filt <- FindVariableFeatures(CN3_seurat_filt, selection.method = "vst", nfeatures = 2000)

##Scale Data##

#NASH 

N1_allgenes <- rownames(N1_seurat_filt)
N1_seurat_filt <- ScaleData(N1_seurat_filt, features = N1_allgenes)

N2_allgenes <- rownames(N2_seurat_filt)
N2_seurat_filt <- ScaleData(N2_seurat_filt, features = N2_allgenes)

N3_allgenes <- rownames(N3_seurat_filt)
N3_seurat_filt <- ScaleData(N3_seurat_filt, features = N3_allgenes)

N4_allgenes <- rownames(N4_seurat_filt)
N4_seurat_filt <- ScaleData(N4_seurat_filt, features = N4_allgenes)

#NASH Controls 

CN1_allgenes <- rownames(CN1_seurat_filt)
CN1_seurat_filt <- ScaleData(CN1_seurat_filt, features = CN1_allgenes)

CN2_allgenes <- rownames(CN2_seurat_filt)
CN2_seurat_filt <- ScaleData(CN2_seurat_filt, features = CN2_allgenes)

CN3_allgenes <- rownames(CN3_seurat_filt)
CN3_seurat_filt <- ScaleData(CN3_seurat_filt, features = CN3_allgenes)

#### Principal Component Analysis (PCA) - Linear Dimension Reduction ####

#NASH 

N1_seurat_filt <- RunPCA(N1_seurat_filt, features = VariableFeatures(object = N1_seurat_filt))
print(N1_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5) 
DimHeatmap(N1_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(N1_seurat_filt)+
ggtitle("N1 Elbow Plot")

N2_seurat_filt <- RunPCA(N2_seurat_filt, features = VariableFeatures(object = N2_seurat_filt))
print(N2_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5) 
DimHeatmap(N2_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(N2_seurat_filt)+
ggtitle("N2 Elbow Plot")

N3_seurat_filt <- RunPCA(N3_seurat_filt, features = VariableFeatures(object = N3_seurat_filt))
print(N3_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5) 
DimHeatmap(N3_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(N3_seurat_filt)+
ggtitle("N3 Elbow Plot")

N4_seurat_filt <- RunPCA(N4_seurat_filt, features = VariableFeatures(object = N4_seurat_filt))
print(N4_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5) 
DimHeatmap(N4_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(N4_seurat_filt)+
ggtitle("N4 Elbow Plot")

#NASH Controls 

CN1_seurat_filt <- RunPCA(CN1_seurat_filt, features = VariableFeatures(object = CN1_seurat_filt))
print(CN1_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5) 
DimHeatmap(CN1_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(CN1_seurat_filt)+
ggtitle("CN1 Elbow Plot")

CN2_seurat_filt <- RunPCA(CN2_seurat_filt, features = VariableFeatures(object = CN2_seurat_filt))
print(CN2_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5) 
DimHeatmap(CN2_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(CN2_seurat_filt)+
ggtitle("CN2 Elbow Plot")

CN3_seurat_filt <- RunPCA(CN3_seurat_filt, features = VariableFeatures(object = CN3_seurat_filt))
print(CN3_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5) 
DimHeatmap(CN3_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(CN3_seurat_filt)+
ggtitle("CN3 Elbow Plot")

#NASH 

N1_seurat_filt <- FindNeighbors(N1_seurat_filt, dims = 1:6)
N1_seurat_filt <- FindClusters(object = N1_seurat_filt)
N1_seurat_filt <- RunUMAP(N1_seurat_filt, dims = 1:6)

DimPlot(N1_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

N2_seurat_filt <- FindNeighbors(N2_seurat_filt, dims = 1:7)
N2_seurat_filt <- FindClusters(object = N2_seurat_filt)
N2_seurat_filt <- RunUMAP(N2_seurat_filt, dims = 1:7)

DimPlot(N2_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

N3_seurat_filt <- FindNeighbors(N3_seurat_filt, dims = 1:6)
N3_seurat_filt <- FindClusters(object = N3_seurat_filt)
N3_seurat_filt <- RunUMAP(N3_seurat_filt, dims = 1:6)

DimPlot(N3_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

N4_seurat_filt <- FindNeighbors(N4_seurat_filt, dims = 1:5)
N4_seurat_filt <- FindClusters(object = N4_seurat_filt)
N4_seurat_filt <- RunUMAP(N4_seurat_filt, dims = 1:5)

DimPlot(N4_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

CN1_seurat_filt <- FindNeighbors(CN1_seurat_filt, dims = 1:5)
CN1_seurat_filt <- FindClusters(object = CN1_seurat_filt)
CN1_seurat_filt <- RunUMAP(CN1_seurat_filt, dims = 1:5)

DimPlot(CN1_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

CN2_seurat_filt <- FindNeighbors(CN2_seurat_filt, dims = 1:7)
CN2_seurat_filt <- FindClusters(object = CN2_seurat_filt)
CN2_seurat_filt <- RunUMAP(CN2_seurat_filt, dims = 1:7)

DimPlot(CN2_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

CN3_seurat_filt <- FindNeighbors(CN3_seurat_filt, dims = 1:5)
CN3_seurat_filt <- FindClusters(object = CN3_seurat_filt)
CN3_seurat_filt <- RunUMAP(CN3_seurat_filt, dims = 1:5)

DimPlot(CN3_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

#N1
N1_sweep_res_list <- paramSweep(N1_seurat_filt, PCs = 1:6, sct = FALSE)
N1_sweep_stats <- summarizeSweep(N1_sweep_res_list, GT = FALSE)
N1_find_pK <- find.pK(N1_sweep_stats)

#Visualization
ggplot(N1_find_pK, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  ggtitle("N1 pK Sweep")

#Store optimal pK Value
N1_pK <- N1_find_pK %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
N1_pK <- as.numeric(as.character(N1_pK[[1]]))
cat("Optimal pK for N1:", N1_pK, "\n")

# Doublet rate estimation
N1_cell_count <- nrow(N1_seurat_filt@meta.data)
N1_doublet_rate <- (N1_cell_count / 1000) * 0.008  # 0.8% per 1000 cells
cat("N1 cell count:", N1_cell_count, "\n")
cat("Estimated doublet rate:", round(N1_doublet_rate * 100, 1), "%\n")

# Homotypic doublet adjustment
N1_annotations <- N1_seurat_filt@meta.data$seurat_clusters
N1_homotypic_prop <- modelHomotypic(N1_annotations)
N1_nExp.poi <- round(N1_doublet_rate * N1_cell_count)
N1_nExp.poi.adj <- round(N1_nExp.poi * (1 - N1_homotypic_prop))
cat("Expected doublets (adjusted):", N1_nExp.poi.adj, "\n")

# Run DoubletFinder
N1_seurat_filt <- doubletFinder(N1_seurat_filt, 
                                PCs = 1:6, 
                                pN = 0.25, 
                                pK = N1_pK, 
                                nExp = N1_nExp.poi.adj,
                                reuse.pANN = FALSE, 
                                sct = FALSE)

# Dynamically get the DF column names (instead of hardcoding)
df_cols <- grep("^DF.classifications", colnames(N1_seurat_filt@meta.data), value = TRUE)
cat("DoubletFinder classification column:", df_cols, "\n")

# Visualize doublets
DimPlot(N1_seurat_filt, reduction = 'umap', group.by = df_cols[1])

# Check doublet counts before removing
table(N1_seurat_filt@meta.data[[df_cols[1]]])

# Remove doublets
N1_seurat_filt <- subset(N1_seurat_filt, subset = !!sym(df_cols[1]) == "Singlet")

# Verify removal
cat("Cells after doublet removal:", ncol(N1_seurat_filt), "\n")

# Save
saveRDS(N1_seurat_filt, file = "NASH1_doubremoved.rds")

# ==============================================================================
# N2
# ==============================================================================
N2_sweep_res_list <- paramSweep(N2_seurat_filt, PCs = 1:7, sct = FALSE)
N2_sweep_stats <- summarizeSweep(N2_sweep_res_list, GT = FALSE)
N2_find_pK <- find.pK(N2_sweep_stats)

# Visualization
ggplot(N2_find_pK, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  ggtitle("N2 pK Sweep")

# Store optimal pK Value
N2_pK <- N2_find_pK %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
N2_pK <- as.numeric(as.character(N2_pK[[1]]))
cat("Optimal pK for N2:", N2_pK, "\n")

# Dynamic doublet rate calculation
N2_cell_count <- nrow(N2_seurat_filt@meta.data)
N2_doublet_rate <- (N2_cell_count / 1000) * 0.008
cat("N2 cell count:", N2_cell_count, "\n")
cat("Estimated doublet rate:", round(N2_doublet_rate * 100, 1), "%\n")

# Homotypic doublet adjustment
N2_annotations <- N2_seurat_filt@meta.data$seurat_clusters
N2_homotypic_prop <- modelHomotypic(N2_annotations)
N2_nExp.poi <- round(N2_doublet_rate * N2_cell_count)
N2_nExp.poi.adj <- round(N2_nExp.poi * (1 - N2_homotypic_prop))
cat("Expected doublets (adjusted):", N2_nExp.poi.adj, "\n")

# Run DoubletFinder
N2_seurat_filt <- doubletFinder(N2_seurat_filt, 
                                PCs = 1:7, 
                                pN = 0.25, 
                                pK = N2_pK, 
                                nExp = N2_nExp.poi.adj,
                                reuse.pANN = FALSE, 
                                sct = FALSE)

# Dynamically get DF column names
N2_df_cols <- grep("^DF.classifications", colnames(N2_seurat_filt@meta.data), value = TRUE)
cat("DoubletFinder classification column:", N2_df_cols, "\n")

# Visualize doublets
DimPlot(N2_seurat_filt, reduction = 'umap', group.by = N2_df_cols[1])

# Check doublet counts
table(N2_seurat_filt@meta.data[[N2_df_cols[1]]])

# Remove doublets
N2_seurat_filt <- subset(N2_seurat_filt, subset = !!sym(N2_df_cols[1]) == "Singlet")
cat("Cells after doublet removal:", ncol(N2_seurat_filt), "\n")

# Save
saveRDS(N2_seurat_filt, file = "N2_doubremoved.rds")


# ==============================================================================
# N3
# ==============================================================================
N3_sweep_res_list <- paramSweep(N3_seurat_filt, PCs = 1:6, sct = FALSE)
N3_sweep_stats <- summarizeSweep(N3_sweep_res_list, GT = FALSE)
N3_find_pK <- find.pK(N3_sweep_stats)

# Visualization
ggplot(N3_find_pK, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  ggtitle("N3 pK Sweep")

# Store optimal pK Value
N3_pK <- N3_find_pK %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
N3_pK <- as.numeric(as.character(N3_pK[[1]]))
cat("Optimal pK for N3:", N3_pK, "\n")

# Dynamic doublet rate calculation
N3_cell_count <- nrow(N3_seurat_filt@meta.data)
N3_doublet_rate <- (N3_cell_count / 1000) * 0.008
cat("N3 cell count:", N3_cell_count, "\n")
cat("Estimated doublet rate:", round(N3_doublet_rate * 100, 1), "%\n")

# Homotypic doublet adjustment
N3_annotations <- N3_seurat_filt@meta.data$seurat_clusters
N3_homotypic_prop <- modelHomotypic(N3_annotations)
N3_nExp.poi <- round(N3_doublet_rate * N3_cell_count)
N3_nExp.poi.adj <- round(N3_nExp.poi * (1 - N3_homotypic_prop))
cat("Expected doublets (adjusted):", N3_nExp.poi.adj, "\n")

# Run DoubletFinder
N3_seurat_filt <- doubletFinder(N3_seurat_filt, 
                                PCs = 1:6, 
                                pN = 0.25, 
                                pK = N3_pK, 
                                nExp = N3_nExp.poi.adj,
                                reuse.pANN = FALSE, 
                                sct = FALSE)

# Dynamically get DF column names
N3_df_cols <- grep("^DF.classifications", colnames(N3_seurat_filt@meta.data), value = TRUE)
cat("DoubletFinder classification column:", N3_df_cols, "\n")

# Visualize doublets
DimPlot(N3_seurat_filt, reduction = 'umap', group.by = N3_df_cols[1])

# Check doublet counts
table(N3_seurat_filt@meta.data[[N3_df_cols[1]]])

# Remove doublets
N3_seurat_filt <- subset(N3_seurat_filt, subset = !!sym(N3_df_cols[1]) == "Singlet")
cat("Cells after doublet removal:", ncol(N3_seurat_filt), "\n")

# Save
saveRDS(N3_seurat_filt, file = "NASH3_doubremoved.rds")


# ==============================================================================
# N4
# ==============================================================================
N4_sweep_res_list <- paramSweep(N4_seurat_filt, PCs = 1:5, sct = FALSE)
N4_sweep_stats <- summarizeSweep(N4_sweep_res_list, GT = FALSE)
N4_find_pK <- find.pK(N4_sweep_stats)

# Visualization
ggplot(N4_find_pK, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  ggtitle("N4 pK Sweep")

# Store optimal pK Value
N4_pK <- N4_find_pK %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
N4_pK <- as.numeric(as.character(N4_pK[[1]]))
cat("Optimal pK for N4:", N4_pK, "\n")

# Dynamic doublet rate calculation
N4_cell_count <- nrow(N4_seurat_filt@meta.data)
N4_doublet_rate <- (N4_cell_count / 1000) * 0.008
cat("N4 cell count:", N4_cell_count, "\n")
cat("Estimated doublet rate:", round(N4_doublet_rate * 100, 1), "%\n")

# Homotypic doublet adjustment
N4_annotations <- N4_seurat_filt@meta.data$seurat_clusters
N4_homotypic_prop <- modelHomotypic(N4_annotations)
N4_nExp.poi <- round(N4_doublet_rate * N4_cell_count)
N4_nExp.poi.adj <- round(N4_nExp.poi * (1 - N4_homotypic_prop))
cat("Expected doublets (adjusted):", N4_nExp.poi.adj, "\n")

# Run DoubletFinder
N4_seurat_filt <- doubletFinder(N4_seurat_filt, 
                                PCs = 1:5, 
                                pN = 0.25, 
                                pK = N4_pK, 
                                nExp = N4_nExp.poi.adj,
                                reuse.pANN = FALSE, 
                                sct = FALSE)

# Dynamically get DF column names
N4_df_cols <- grep("^DF.classifications", colnames(N4_seurat_filt@meta.data), value = TRUE)
cat("DoubletFinder classification column:", N4_df_cols, "\n")

# Visualize doublets
DimPlot(N4_seurat_filt, reduction = 'umap', group.by = N4_df_cols[1])

# Check doublet counts
table(N4_seurat_filt@meta.data[[N4_df_cols[1]]])

# Remove doublets
N4_seurat_filt <- subset(N4_seurat_filt, subset = !!sym(N4_df_cols[1]) == "Singlet")
cat("Cells after doublet removal:", ncol(N4_seurat_filt), "\n")

# Save
saveRDS(N4_seurat_filt, file = "NASH4_doubremoved.rds")


# ==============================================================================
# CN1
# ==============================================================================
CN1_sweep_res_list <- paramSweep(CN1_seurat_filt, PCs = 1:5, sct = FALSE)
CN1_sweep_stats <- summarizeSweep(CN1_sweep_res_list, GT = FALSE)
CN1_find_pK <- find.pK(CN1_sweep_stats)

# Visualization
ggplot(CN1_find_pK, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  ggtitle("CN1 pK Sweep")

# Store optimal pK Value
CN1_pK <- CN1_find_pK %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
CN1_pK <- as.numeric(as.character(CN1_pK[[1]]))
cat("Optimal pK for CN1:", CN1_pK, "\n")

# Dynamic doublet rate calculation
CN1_cell_count <- nrow(CN1_seurat_filt@meta.data)
CN1_doublet_rate <- (CN1_cell_count / 1000) * 0.008
cat("CN1 cell count:", CN1_cell_count, "\n")
cat("Estimated doublet rate:", round(CN1_doublet_rate * 100, 1), "%\n")

# Homotypic doublet adjustment
CN1_annotations <- CN1_seurat_filt@meta.data$seurat_clusters
CN1_homotypic_prop <- modelHomotypic(CN1_annotations)
CN1_nExp.poi <- round(CN1_doublet_rate * CN1_cell_count)
CN1_nExp.poi.adj <- round(CN1_nExp.poi * (1 - CN1_homotypic_prop))
cat("Expected doublets (adjusted):", CN1_nExp.poi.adj, "\n")

# Run DoubletFinder
CN1_seurat_filt <- doubletFinder(CN1_seurat_filt, 
                                 PCs = 1:5, 
                                 pN = 0.25, 
                                 pK = CN1_pK, 
                                 nExp = CN1_nExp.poi.adj,
                                 reuse.pANN = FALSE, 
                                 sct = FALSE)

# Dynamically get DF column names
CN1_df_cols <- grep("^DF.classifications", colnames(CN1_seurat_filt@meta.data), value = TRUE)
cat("DoubletFinder classification column:", CN1_df_cols, "\n")

# Visualize doublets
DimPlot(CN1_seurat_filt, reduction = 'umap', group.by = CN1_df_cols[1])

# Check doublet counts
table(CN1_seurat_filt@meta.data[[CN1_df_cols[1]]])

# Remove doublets
CN1_seurat_filt <- subset(CN1_seurat_filt, subset = !!sym(CN1_df_cols[1]) == "Singlet")
cat("Cells after doublet removal:", ncol(CN1_seurat_filt), "\n")

# Save
saveRDS(CN1_seurat_filt, file = "CNASH1_doubremoved.rds")


# ==============================================================================
# CN2
# ==============================================================================
CN2_sweep_res_list <- paramSweep(CN2_seurat_filt, PCs = 1:7, sct = FALSE)
CN2_sweep_stats <- summarizeSweep(CN2_sweep_res_list, GT = FALSE)
CN2_find_pK <- find.pK(CN2_sweep_stats)

# Visualization
ggplot(CN2_find_pK, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  ggtitle("CN2 pK Sweep")

# Store optimal pK Value
CN2_pK <- CN2_find_pK %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
CN2_pK <- as.numeric(as.character(CN2_pK[[1]]))
cat("Optimal pK for CN2:", CN2_pK, "\n")

# Dynamic doublet rate calculation
CN2_cell_count <- nrow(CN2_seurat_filt@meta.data)
CN2_doublet_rate <- (CN2_cell_count / 1000) * 0.008
cat("CN2 cell count:", CN2_cell_count, "\n")
cat("Estimated doublet rate:", round(CN2_doublet_rate * 100, 1), "%\n")

# Homotypic doublet adjustment
CN2_annotations <- CN2_seurat_filt@meta.data$seurat_clusters
CN2_homotypic_prop <- modelHomotypic(CN2_annotations)
CN2_nExp.poi <- round(CN2_doublet_rate * CN2_cell_count)
CN2_nExp.poi.adj <- round(CN2_nExp.poi * (1 - CN2_homotypic_prop))
cat("Expected doublets (adjusted):", CN2_nExp.poi.adj, "\n")

# Run DoubletFinder
CN2_seurat_filt <- doubletFinder(CN2_seurat_filt, 
                                 PCs = 1:7, 
                                 pN = 0.25, 
                                 pK = CN2_pK, 
                                 nExp = CN2_nExp.poi.adj,
                                 reuse.pANN = FALSE, 
                                 sct = FALSE)

# Dynamically get DF column names
CN2_df_cols <- grep("^DF.classifications", colnames(CN2_seurat_filt@meta.data), value = TRUE)
cat("DoubletFinder classification column:", CN2_df_cols, "\n")

# Visualize doublets
DimPlot(CN2_seurat_filt, reduction = 'umap', group.by = CN2_df_cols[1])

# Check doublet counts
table(CN2_seurat_filt@meta.data[[CN2_df_cols[1]]])

# Remove doublets
CN2_seurat_filt <- subset(CN2_seurat_filt, subset = !!sym(CN2_df_cols[1]) == "Singlet")
cat("Cells after doublet removal:", ncol(CN2_seurat_filt), "\n")

# Save
saveRDS(CN2_seurat_filt, file = "CNASH2_doubremoved.rds")


# ==============================================================================
# CN3
# ==============================================================================
CN3_sweep_res_list <- paramSweep(CN3_seurat_filt, PCs = 1:6, sct = FALSE)
CN3_sweep_stats <- summarizeSweep(CN3_sweep_res_list, GT = FALSE)
CN3_find_pK <- find.pK(CN3_sweep_stats)

# Visualization
ggplot(CN3_find_pK, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  ggtitle("CN3 pK Sweep")

# Store optimal pK Value
CN3_pK <- CN3_find_pK %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
CN3_pK <- as.numeric(as.character(CN3_pK[[1]]))
cat("Optimal pK for CN3:", CN3_pK, "\n")

# Dynamic doublet rate calculation
CN3_cell_count <- nrow(CN3_seurat_filt@meta.data)
CN3_doublet_rate <- (CN3_cell_count / 1000) * 0.008
cat("CN3 cell count:", CN3_cell_count, "\n")
cat("Estimated doublet rate:", round(CN3_doublet_rate * 100, 1), "%\n")

# Homotypic doublet adjustment
CN3_annotations <- CN3_seurat_filt@meta.data$seurat_clusters
CN3_homotypic_prop <- modelHomotypic(CN3_annotations)
CN3_nExp.poi <- round(CN3_doublet_rate * CN3_cell_count)
CN3_nExp.poi.adj <- round(CN3_nExp.poi * (1 - CN3_homotypic_prop))
cat("Expected doublets (adjusted):", CN3_nExp.poi.adj, "\n")

# Run DoubletFinder
CN3_seurat_filt <- doubletFinder(CN3_seurat_filt, 
                                 PCs = 1:6, 
                                 pN = 0.25, 
                                 pK = CN3_pK, 
                                 nExp = CN3_nExp.poi.adj,
                                 reuse.pANN = FALSE, 
                                 sct = FALSE)

# Dynamically get DF column names
CN3_df_cols <- grep("^DF.classifications", colnames(CN3_seurat_filt@meta.data), value = TRUE)
cat("DoubletFinder classification column:", CN3_df_cols, "\n")

# Visualize doublets
DimPlot(CN3_seurat_filt, reduction = 'umap', group.by = CN3_df_cols[1])

# Check doublet counts
table(CN3_seurat_filt@meta.data[[CN3_df_cols[1]]])

# Remove doublets
CN3_seurat_filt <- subset(CN3_seurat_filt, subset = !!sym(CN3_df_cols[1]) == "Singlet")
cat("Cells after doublet removal:", ncol(CN3_seurat_filt), "\n")

# Save
saveRDS(CN3_seurat_filt, file = "CNASH3_doubremoved.rds")

N1 <- readRDS(file = "/projects/standard/prin0088/shared/NASH1_doubremoved.rds")
N2 <- readRDS(file = "/projects/standard/prin0088/shared/N2_doubremoved.rds")
N3 <- readRDS(file = "/projects/standard/prin0088/shared/NASH3_doubremoved.rds")
N4 <- readRDS(file = "/projects/standard/prin0088/shared/NASH4_doubremoved.rds")
CN1 <- readRDS(file = "/projects/standard/prin0088/shared/CNASH1_doubremoved.rds")
CN2 <- readRDS(file = "/projects/standard/prin0088/shared/CNASH2_doubremoved.rds")
CN3 <- readRDS(file = "/projects/standard/prin0088/shared/CNASH3_doubremoved.rds")

# Check cell counts before merging
cat("N1:", ncol(N1), "\n")
cat("N2:", ncol(N2), "\n")
cat("N3:", ncol(N3), "\n")
cat("N4:", ncol(N4), "\n")
cat("CN1:", ncol(CN1), "\n")
cat("CN2:", ncol(CN2), "\n")
cat("CN3:", ncol(CN3), "\n")

# Merge with correct variable names
new_total_matrix_NASH_filt <- merge(N1, y = list(N2, N3, N4, CN1, CN2, CN3), 
                                    add.cell.ids = c("N1","N2", "N3", "N4","CN1","CN2","CN3"), 
                                    project = "All.Samples")

cat("Total merged:", ncol(new_total_matrix_NASH_filt), "\n")

saveRDS(new_total_matrix_NASH_filt, file = 'total_matrix_NASH_filt')


