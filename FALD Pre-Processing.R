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

##Fontan Livers
F1_mtx <- ReadMtx(mtx = "/users/1/blake561/Fontan Files/GSM6997746_GEX_Fontan1_matrix.mtx.gz",
                  features = "/users/1/blake561/Fontan Files/GSM6997746_GEX_Fontan1_features.tsv.gz",
                  cells = "/users/1/blake561/Fontan Files/GSM6997746_GEX_Fontan1_barcodes.tsv.gz")

F1_seurat_obj <- CreateSeuratObject(counts = F1_mtx, project = "Fontan", min.cells = 3, min.features = 200)

F2_mtx <- ReadMtx(mtx = "/users/1/blake561/Fontan Files/GSM6997748_GEX_Fontan2_matrix.mtx.gz",
                  features = "/users/1/blake561/Fontan Files/GSM6997748_GEX_Fontan2_features.tsv.gz",
                  cells = "/users/1/blake561/Fontan Files/GSM6997748_GEX_Fontan2_barcodes.tsv.gz")

F2_seurat_obj <- CreateSeuratObject(counts = F2_mtx, project = "Fontan", min.cells = 3, min.features = 200)

F3_mtx <- ReadMtx(mtx = "/users/1/blake561/Fontan Files/GSM6997750_GEX_Fontan3_matrix.mtx.gz",
                  features = "/users/1/blake561/Fontan Files/GSM6997750_GEX_Fontan3_features.tsv.gz",
                  cells = "/users/1/blake561/Fontan Files/GSM6997750_GEX_Fontan3_barcodes.tsv.gz")

F3_seurat_obj <- CreateSeuratObject(counts = F3_mtx, project = "Fontan", min.cells = 3, min.features = 200)

F4_mtx <- ReadMtx(mtx = "/users/1/blake561/Fontan Files/GSM6997752_GEX_Fontan4_matrix.mtx.gz",
                  features = "/users/1/blake561/Fontan Files/GSM6997752_GEX_Fontan4_features.tsv.gz",
                  cells = "/users/1/blake561/Fontan Files/GSM6997752_GEX_Fontan4_barcodes.tsv.gz")

F4_seurat_obj <- CreateSeuratObject(counts = F4_mtx, project = "Fontan", min.cells = 3, min.features = 200)

##Fontan Controls

CF1_mtx <- ReadMtx(mtx = "/users/1/blake561/Fontan Files/GSM6997741_GEX_Ctrl_matrix.mtx.gz",
                   features = "/users/1/blake561/Fontan Files/GSM6997741_GEX_Ctrl_features.tsv.gz",
                   cells = "/users/1/blake561/Fontan Files/GSM6997741_GEX_Ctrl_barcodes.tsv.gz")

CF1_seurat_obj <- CreateSeuratObject(counts = CF1_mtx, project = "Fontan Control", min.cells = 3, min.features = 200)

CF2_mtx <- ReadMtx(mtx = "/users/1/blake561/Fontan Files/GSM6997744_GEX_Ctrl3_matrix.mtx.gz",
                   features = "/users/1/blake561/Fontan Files/GSM6997744_GEX_Ctrl3_features.tsv.gz",
                   cells = "/users/1/blake561/Fontan Files/GSM6997744_GEX_Ctrl3_barcodes.tsv.gz")

CF2_seurat_obj <- CreateSeuratObject(counts = CF2_mtx, project = "Fontan Control", min.cells = 3, min.features = 200)


#Fontan

F1_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(F1_seurat_obj, pattern = "^MT")
F2_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(F2_seurat_obj, pattern = "^MT")
F3_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(F3_seurat_obj, pattern = "^MT")
F4_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(F4_seurat_obj, pattern = "^MT")

#Fontan Controls

CF1_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(CF1_seurat_obj, pattern = "^MT")
CF2_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(CF2_seurat_obj, pattern = "^MT")

##Normalization

#Fontan 

F1_seurat_filt <- subset(F1_seurat_obj, subset = nFeature_RNA > 800 & #Filtering parameters from manuscript methods. 
                           nFeature_RNA < 20000 &
                           nCount_RNA < 25000 &
                           percent.mt < 10)
F1_seurat_filt <- NormalizeData(F1_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

F2_seurat_filt <- subset(F2_seurat_obj, subset = nFeature_RNA > 800 & 
                           nFeature_RNA < 20000 &
                           nCount_RNA < 25000 &
                           percent.mt < 10)
F2_seurat_filt <- NormalizeData(F2_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

F3_seurat_filt <- subset(F3_seurat_obj, subset = nFeature_RNA > 800 & #Altered Filtering Bc Sample is Much Messier than the Others
                           nFeature_RNA < 20000 &
                           nCount_RNA < 25000 &
                           percent.mt < 10)
F3_seurat_filt <- NormalizeData(F3_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

F4_seurat_filt <- subset(F4_seurat_obj, subset = nFeature_RNA > 800 & 
                           nFeature_RNA < 20000 &
                           nCount_RNA < 25000 &
                           percent.mt < 10)
F4_seurat_filt <- NormalizeData(F4_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

#Fontan Controls

CF1_seurat_filt <- subset(CF1_seurat_obj, subset = nFeature_RNA > 800 & 
                            nFeature_RNA < 20000 &
                            nCount_RNA < 25000 &
                            percent.mt < 10)
CF1_seurat_filt <- NormalizeData(CF1_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)

CF2_seurat_filt <- subset(CF2_seurat_obj, subset = nFeature_RNA > 800 & 
                            nFeature_RNA < 20000 &
                            nCount_RNA < 25000 &
                            percent.mt < 10) 
CF2_seurat_filt <- NormalizeData(CF2_seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)


VlnPlot(F1_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(F2_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(F3_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(F4_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(CF1_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CF2_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



ggplot(F3_seurat_filt@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point() + xlab("Count Depth") + ylab("Number of genes")

##Find Variable Features##

#Fontan 

F1_seurat_filt <- FindVariableFeatures(F1_seurat_filt, selection.method = "vst", nfeatures = 2000)
F2_seurat_filt <- FindVariableFeatures(F2_seurat_filt, selection.method = "vst", nfeatures = 2000)
F3_seurat_filt <- FindVariableFeatures(F3_seurat_filt, selection.method = "vst", nfeatures = 2000)
F4_seurat_filt <- FindVariableFeatures(F4_seurat_filt, selection.method = "vst", nfeatures = 2000)

#Fontan Controls

CF1_seurat_filt <- FindVariableFeatures(CF1_seurat_filt, selection.method = "vst", nfeatures = 2000)
CF2_seurat_filt <- FindVariableFeatures(CF2_seurat_filt, selection.method = "vst", nfeatures = 2000)

##Scale Data##

#Fontan 

F1_allgenes <- rownames(F1_seurat_filt)
F1_seurat_filt <- ScaleData(F1_seurat_filt, features = F1_allgenes)

F2_allgenes <- rownames(F2_seurat_filt)
F2_seurat_filt <- ScaleData(F2_seurat_filt, features = F2_allgenes)

F3_allgenes <- rownames(F3_seurat_filt)
F3_seurat_filt <- ScaleData(F3_seurat_filt, features = F3_allgenes)

F4_allgenes <- rownames(F4_seurat_filt)
F4_seurat_filt <- ScaleData(F4_seurat_filt, features = F4_allgenes)

#Fontan Controls

CF1_allgenes <- rownames(CF1_seurat_filt)
CF1_seurat_filt <- ScaleData(CF1_seurat_filt, features = CF1_allgenes)

CF2_allgenes <- rownames(CF2_seurat_filt)
CF2_seurat_filt <- ScaleData(CF2_seurat_filt, features = CF2_allgenes)


#### Principal Component Analysis (PCA) - Linear Dimension Reduction ####

#Fontan 

F1_seurat_filt <- RunPCA(F1_seurat_filt, features = VariableFeatures(object = F1_seurat_filt))
print(F1_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5) 
DimHeatmap(F1_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(F1_seurat_filt)+
  ggtitle("F1 Elbow Plot")

F2_seurat_filt <- RunPCA(F2_seurat_filt, features = VariableFeatures(object = F2_seurat_filt))
print(F2_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5) 
DimHeatmap(F2_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(F2_seurat_filt)+
  ggtitle("F2 Elbow Plot")

F3_seurat_filt <- RunPCA(F3_seurat_filt, features = VariableFeatures(object = F3_seurat_filt))
print(F3_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5) 
DimHeatmap(F3_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(F3_seurat_filt)+
  ggtitle("F3 Elbow Plot")

F4_seurat_filt <- RunPCA(F4_seurat_filt, features = VariableFeatures(object = F4_seurat_filt))
print(F4_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5) 
DimHeatmap(F4_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(F4_seurat_filt)+
  ggtitle("F4 Elbow Plot")

#Fontan Controls

CF1_seurat_filt <- RunPCA(CF1_seurat_filt, features = VariableFeatures(object = CF1_seurat_filt))
print(CF1_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5) 
DimHeatmap(CF1_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(CF1_seurat_filt)+
ggtitle("CF1 Elbow Plot")

CF2_seurat_filt <- RunPCA(CF2_seurat_filt, features = VariableFeatures(object = CF2_seurat_filt))
print(CF2_seurat_filt[["pca"]], dims = 1:5, nfeatures = 5) 
DimHeatmap(CF2_seurat_filt, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(CF2_seurat_filt)+
ggtitle("CF2 Elbow Plot")

#### Clustering ####

#Fontan 

F1_seurat_filt <- FindNeighbors(F1_seurat_filt, dims = 1:10) #14
F1_seurat_filt <- FindClusters(object = F1_seurat_filt)
F1_seurat_filt <- RunUMAP(F1_seurat_filt, dims = 1:10)

DimPlot(F1_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

F2_seurat_filt <- FindNeighbors(F2_seurat_filt, dims = 1:8) #10 #16
F2_seurat_filt <- FindClusters(object = F2_seurat_filt)
F2_seurat_filt <- RunUMAP(F2_seurat_filt, dims = 1:8)

DimPlot(F2_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

F3_seurat_filt <- FindNeighbors(F3_seurat_filt, dims = 1:8) #13
F3_seurat_filt <- FindClusters(object = F3_seurat_filt)
F3_seurat_filt <- RunUMAP(F3_seurat_filt, dims = 1:8)

DimPlot(F3_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

F4_seurat_filt <- FindNeighbors(F4_seurat_filt, dims = 1:9) #14
F4_seurat_filt <- FindClusters(object = F4_seurat_filt)
F4_seurat_filt <- RunUMAP(F4_seurat_filt, dims = 1:9)

DimPlot(F4_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)

#Fontan Controls

CF1_seurat_filt <- FindNeighbors(CF1_seurat_filt, dims = 1:6) #12
CF1_seurat_filt <- FindClusters(object = CF1_seurat_filt)
CF1_seurat_filt <- RunUMAP(CF1_seurat_filt, dims = 1:6)

DimPlot(CF1_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)


CF2_seurat_filt <- FindNeighbors(CF2_seurat_filt, dims = 1:10) #12
CF2_seurat_filt <- FindClusters(object = CF2_seurat_filt)
CF2_seurat_filt <- RunUMAP(CF2_seurat_filt, dims = 1:10)

DimPlot(CF2_seurat_filt, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)


##Doublet Finder##

#Fontan Samples

#F1

F1_sweep_res_list <- paramSweep(F1_seurat_filt, PCs = 1:10, sct = FALSE)
F1_sweep_stats <- summarizeSweep(F1_sweep_res_list, GT = FALSE)
F1_find_pK <- find.pK(F1_sweep_stats)

#Visualization
ggplot(F1_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
F1_pK <- F1_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
F1_pK <- as.numeric(as.character(F1_pK[[1]]))

F1_annotations <- F1_seurat_filt@meta.data$seurat_clusters
F1_homotypic_prop <- modelHomotypic(F1_annotations)
F1_nExp.poi <- round(0.051984*nrow(F1_seurat_filt@meta.data)) # 6228 cells recovered, %0.8 for every 1000 cells recovered (i.e 25,000 --> calculation for red number, 25 * 0.8 = 20% --> 0.2; really was 24,600 but round to nearest whole number)
F1_nEXP.poi.adj <- round(F1_nExp.poi *(1-F1_homotypic_prop))

F1_seurat_filt <- doubletFinder(F1_seurat_filt, PCs = 1:10, pN = 0.25, pK = F1_pK, nExp = F1_nEXP.poi.adj,
                                reuse.pANN = FALSE, sct = FALSE)
DimPlot(F1_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_304") #Change the green name by finding the DF.classifcations file name by typing view(C1_seurat_filt@meta.data)

F1_seurat_filt <- subset(F1_seurat_filt, subset = DF.classifications_0.25_0.005_304 == "Singlet") #Change this DF.classifcations as well! 
saveRDS(F1_seurat_filt, file = "F1_doubremoved")


#F2

F2_sweep_res_list <- paramSweep(F2_seurat_filt, PCs = 1:8, sct = FALSE)
F2_sweep_stats <- summarizeSweep(F2_sweep_res_list, GT = FALSE)
F2_find_pK <- find.pK(F2_sweep_stats)

#Visualization
ggplot(F2_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
F2_pK <- F2_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
F2_pK <- as.numeric(as.character(F2_pK[[1]]))

F2_annotations <- F2_seurat_filt@meta.data$seurat_clusters
F2_homotypic_prop <- modelHomotypic(F2_annotations)
F2_nExp.poi <- round(0.039712*nrow(F2_seurat_filt@meta.data)) # 3928 cells recovered, %0.8 for every 1000 cells recovered (i.e 25,000 --> calculation for red number, 25 * 0.8 = 20% --> 0.2; really was 24,600 but round to nearest whole number)
F2_nEXP.poi.adj <- round(F2_nExp.poi *(1-F2_homotypic_prop))

F2_seurat_filt <- doubletFinder(F2_seurat_filt, PCs = 1:8, pN = 0.25, pK = F2_pK, nExp = F2_nEXP.poi.adj,
                                reuse.pANN = FALSE, sct = FALSE)
DimPlot(F2_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_182") #Change the green name by finding the DF.classifcations file name by typing view(C1_seurat_filt@meta.data)

F2_seurat_filt <- subset(F2_seurat_filt, subset = DF.classifications_0.25_0.01_138 == "Singlet") #Change this DF.classifcations as well! 
saveRDS(F2_seurat_filt, file = "F2_doubremoved")

#F3

F3_sweep_res_list <- paramSweep(F3_seurat_filt, PCs = 1:8, sct = FALSE)
F3_sweep_stats <- summarizeSweep(F3_sweep_res_list, GT = FALSE)
F3_find_pK <- find.pK(F3_sweep_stats)

#Visualization
ggplot(F3_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
F3_pK <- F3_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
F3_pK <- as.numeric(as.character(F3_pK[[1]]))

F3_annotations <- F3_seurat_filt@meta.data$seurat_clusters
F3_homotypic_prop <- modelHomotypic(F3_annotations)
F3_nExp.poi <- round(0.082256*nrow(F3_seurat_filt@meta.data)) # 3928 cells recovered, %0.8 for every 1000 cells recovered (i.e 25,000 --> calculation for red number, 25 * 0.8 = 20% --> 0.2; really was 24,600 but round to nearest whole number)
F3_nEXP.poi.adj <- round(F3_nExp.poi *(1-F3_homotypic_prop))

F3_seurat_filt <- doubletFinder(F3_seurat_filt, PCs = 1:8, pN = 0.25, pK = F3_pK, nExp = F3_nEXP.poi.adj,
                                reuse.pANN = FALSE, sct = FALSE)
DimPlot(F3_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_764") #Change the green name by finding the DF.classifcations file name by typing view(C1_seurat_filt@meta.data)

F3_seurat_filt <- subset(F3_seurat_filt, subset = DF.classifications_0.25_0.005_764 == "Singlet") #Change this DF.classifcations as well! 
saveRDS(F3_seurat_filt, file = "F3_doubremoved")

#F4

F4_sweep_res_list <- paramSweep(F4_seurat_filt, PCs = 1:9, sct = FALSE)
F4_sweep_stats <- summarizeSweep(F4_sweep_res_list, GT = FALSE)
F4_find_pK <- find.pK(F4_sweep_stats)

#Visualization
ggplot(F4_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
F4_pK <- F4_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
F4_pK <- as.numeric(as.character(F4_pK[[1]]))

F4_annotations <- F4_seurat_filt@meta.data$seurat_clusters
F4_homotypic_prop <- modelHomotypic(F4_annotations)
F4_nExp.poi <- round(0.0259*nrow(F4_seurat_filt@meta.data)) # 3928 cells recovered, %0.8 for every 1000 cells recovered (i.e 25,000 --> calculation for red number, 25 * 0.8 = 20% --> 0.2; really was 24,600 but round to nearest whole number)
F4_nEXP.poi.adj <- round(F4_nExp.poi *(1-F4_homotypic_prop))

F4_seurat_filt <- doubletFinder(F4_seurat_filt, PCs = 1:9, pN = 0.25, pK = F4_pK, nExp = F4_nEXP.poi.adj,
                                reuse.pANN = FALSE, sct = FALSE)
DimPlot(F4_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_101") #Change the green name by finding the DF.classifcations file name by typing view(C1_seurat_filt@meta.data)

F4_seurat_filt <- subset(F4_seurat_filt, subset = DF.classifications_0.25_0.005_101 == "Singlet") #Change this DF.classifcations as well! 
saveRDS(F4_seurat_filt, file = "F4_doubremoved")

#Fontan Controls 

#CF1

CF1_sweep_res_list <- paramSweep(CF1_seurat_filt, PCs = 1:6, sct = FALSE)
CF1_sweep_stats <- summarizeSweep(CF1_sweep_res_list, GT = FALSE)
CF1_find_pK <- find.pK(CF1_sweep_stats)

#Visualization
ggplot(CF1_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
CF1_pK <- CF1_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
CF1_pK <- as.numeric(as.character(CF1_pK[[1]]))

CF1_annotations <- CF1_seurat_filt@meta.data$seurat_clusters
CF1_homotypic_prop <- modelHomotypic(CF1_annotations)
CF1_nExp.poi <- round(0.054248*nrow(CF1_seurat_filt@meta.data)) # 6228 cells recovered, %0.8 for every 1000 cells recovered (i.e 25,000 --> calculation for red number, 25 * 0.8 = 20% --> 0.2; really was 24,600 but round to nearest whole number)
CF1_nEXP.poi.adj <- round(CF1_nExp.poi *(1-CF1_homotypic_prop))

CF1_seurat_filt <- doubletFinder(CF1_seurat_filt, PCs = 1:6, pN = 0.25, pK = CF1_pK, nExp = CF1_nEXP.poi.adj,
                                 reuse.pANN = FALSE, sct = FALSE)
DimPlot(CF1_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_337") #Change the green name by finding the DF.classifcations file name by typing view(C1_seurat_filt@meta.data)

CF1_seurat_filt <- subset(CF1_seurat_filt, subset =  DF.classifications_0.25_0.005_337 == "Singlet") #Change this DF.classifcations as well! 
saveRDS(CF1_seurat_filt, file = "CF1_doubremoved")

#CF2

CF2_sweep_res_list <- paramSweep(CF2_seurat_filt, PCs = 1:10, sct = FALSE)
CF2_sweep_stats <- summarizeSweep(CF2_sweep_res_list, GT = FALSE)
CF2_find_pK <- find.pK(CF2_sweep_stats)

#Visualization
ggplot(CF2_find_pK, aes(pK, BCmetric, group = 1))+
  geom_point()+
  geom_line()

#Store Value
CF2_pK <- CF2_find_pK%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
CF2_pK <- as.numeric(as.character(CF2_pK[[1]]))

CF2_annotations <- CF2_seurat_filt@meta.data$seurat_clusters
CF2_homotypic_prop <- modelHomotypic(CF2_annotations)
CF2_nExp.poi <- round(0.127112*nrow(CF2_seurat_filt@meta.data)) # 6228 cells recovered, %0.8 for every 1000 cells recovered (i.e 25,000 --> calculation for red number, 25 * 0.8 = 20% --> 0.2; really was 24,600 but round to nearest whole number)
CF2_nEXP.poi.adj <- round(CF2_nExp.poi *(1-CF2_homotypic_prop))

CF2_seurat_filt <- doubletFinder(CF2_seurat_filt, PCs = 1:10, pN = 0.25, pK = CF2_pK, nExp = CF2_nEXP.poi.adj,
                                 reuse.pANN = FALSE, sct = FALSE)
DimPlot(CF2_seurat_filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.01_1807") #Change the green name by finding the DF.classifcations file name by typing view(C1_seurat_filt@meta.data)

CF2_seurat_filt <- subset(CF2_seurat_filt, subset = DF.classifications_0.25_0.01_1807 == "Singlet") #Change this DF.classifcations as well! 
saveRDS(CF2_seurat_filt, file = "CF2_doubremoved")

#Combining for Integration


#Fontan

total_matrix_fontan <- merge(F1_seurat_filt, y = list(F2_seurat_filt, F3_seurat_filt, F4_seurat_filt, 
                                                      CF1_seurat_filt, CF2_seurat_filt), 
                             add.cell.ids = c("F1","F2", "F3", "F4", "CF1", "CF2"), project = ("All.Samples"))
saveRDS(total_matrix_fontan, file = 'total_matrix_fontan')

CF2_seurat_filt <- readRDS(file="/users/1/blake561/CF2_doubremoved")

