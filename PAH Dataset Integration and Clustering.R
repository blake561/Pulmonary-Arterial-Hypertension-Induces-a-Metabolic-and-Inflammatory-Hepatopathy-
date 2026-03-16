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

total.mtx <- readRDS(file = "/users/1/blake561/total_matrix")


#Remove Column Artifacts from Doublet Finder

total.mtx$pANN_0.25_0.005_2027 <- NULL #C1
total.mtx$DF.classifications_0.25_0.005_2027 <- NULL
total.mtx$pANN_0.25_0.005_702 <- NULL #C2
total.mtx$DF.classifications_0.25_0.005_702 <- NULL
total.mtx$pANN_0.25_0.005_2104 <- NULL #C4
total.mtx$DF.classifications_0.25_0.005_2104 <- NULL
total.mtx$pANN_0.25_0.23_1220 <- NULL #C5
total.mtx$DF.classifications_0.25_0.23_1220 <- NULL
total.mtx$pANN_0.25_0.005_876 <- NULL #P1
total.mtx$DF.classifications_0.25_0.005_876 <- NULL
total.mtx$pANN_0.25_0.005_1823 <- NULL #P2
total.mtx$DF.classifications_0.25_0.005_1823 <- NULL
total.mtx$pANN_0.25_0.005_2332 <- NULL #P3
total.mtx$ DF.classifications_0.25_0.005_2332 <- NULL
total.mtx$pANN_0.25_0.005_2071 <- NULL #P4
total.mtx$DF.classifications_0.25_0.005_2071 <- NULL
total.mtx$pANN_0.25_0.005_3283 <- NULL #P5
total.mtx$DF.classifications_0.25_0.005_3283 <- NULL


total.mtx_filtered <- subset(total.mtx, subset = nCount_RNA > 500 &
                               nFeature_RNA > 200 &
                               percent.mt < 5)


#### Standard Workflow / Pre-processing #### 
total.mtx_filtered <- NormalizeData(object = total.mtx_filtered)
total.mtx_filtered <- FindVariableFeatures(object = total.mtx_filtered)
total.mtx_filtered <- ScaleData(object = total.mtx_filtered)
total.mtx_filtered <- RunPCA(object = total.mtx_filtered)

#print(total.mtx_filtered[["pca"]], dims = 1:5, nfeatures = 5)
#DimHeatmap(total.mtx_filtered, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(total.mtx_filtered)+
ggtitle("Total Matrix Elbow Plot")


#ElbowPlot(total.mtx_filtered)
total.mtx_filtered <- FindNeighbors(object = total.mtx_filtered, dims = 1:15)
total.mtx_filtered <- FindClusters(object = total.mtx_filtered, dims = 1:15)
total.mtx_filtered <- RunUMAP(object = total.mtx_filtered, dims = 1:15)


#DimPlot(total.mtx_filtered, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)
#DimPlot(total.mtx_filtered, reduction = "umap", group.by = "orig.ident", raster = FALSE, label = TRUE)

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
  RunUMAP(reduction = 'harmony', dims = 1:15)%>%
  FindNeighbors(reduction = 'harmony', dims = 1:15)%>%
  FindClusters(resolution = 0.5)

p3 <- DimPlot(total.mtx_harmony, reduction = "umap", group.by = "sample")+
  ggtitle("Integrated With Harmony")
p4 <- DimPlot(total.mtx_harmony, reduction = "umap")

grid.arrange(p1,p3, ncol = 2)

#### Cluster Identification & Annotation #### 

total.mtx_harmony <- FindClusters(total.mtx_harmony, resolution = c(0.1,0.3,0.4,0.5,0.7,0.9,1))

# Visualize

#DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE, raster = FALSE) # 9 Clusters
#DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.3", label = TRUE, raster = FALSE) # 16 Clusters
#DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.4", label = TRUE, raster = FALSE) # 18 Clusters
#DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.5", label = TRUE, raster = FALSE) # 21 Clusters
#DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.7", label = TRUE, raster = FALSE) # 23 Clusters
#DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.0.9", label = TRUE, raster = FALSE) # 27 Clusters
#DimPlot(total.mtx_harmony, reduction = "umap", group.by = "RNA_snn_res.1", label = TRUE, raster = FALSE) # 31 Clusters

#Load in Human Liver Azimuth Reference - Need to Download Files from Azimuth Database and Upload to MSI in Reference Folder

Idents(total.mtx_harmony) <- total.mtx_harmony$RNA_snn_res.0.7 #Do print afterwards to confirm 23 clusters, because that is what we see at resolution of 0.7

liverref <- LoadReference("/users/1/blake561/Reference")

total.mtx_harmony <- JoinLayers(total.mtx_harmony)
total.mtx_harmony <- RunAzimuth(total.mtx_harmony, reference = "/users/1/blake561/Reference")

Idents(total.mtx_harmony) <- "predicted.celltype.l1"
p5 <- DimPlot(total.mtx_harmony, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, raster = FALSE)
p6 <- DimPlot(total.mtx_harmony, reduction = "umap", group.by = "predicted.celltype.l2", label = TRUE, raster = FALSE)

#Save the File/Checkpoint

saveRDS(total.mtx_harmony, file = "total.mtx_harmony")
total.mtx_harmony <- readRDS("/users/1/blake561/total.mtx_harmony")

#UMAPs for Control vs. PAH

#Grouping by Condition

total.mtx_harmony$condition <- total.mtx_harmony$orig.ident

total.mtx_PAH <- subset(total.mtx_harmony, condition == "PAH")
p1 <- DimPlot(total.mtx_PAH, reduction = "umap", group.by = "condition", raster = FALSE, label = TRUE)

total.mtx_Control <- subset(total.mtx_harmony, condition == "Control")
p2 <- DimPlot(total.mtx_Control, reduction = "umap", group.by = "condition", raster = FALSE, label = TRUE)


#Grouping by Cell Type 

#Basic Prediction

total.mtx_PAH <- subset(total.mtx_harmony, condition == "PAH")
p1 <- DimPlot(total.mtx_PAH, reduction = "umap", group.by = "predicted.celltype.l1", raster = FALSE, label = TRUE)

total.mtx_Control <- subset(total.mtx_harmony, condition == "Control")
p2 <- DimPlot(total.mtx_Control, reduction = "umap", group.by = "predicted.celltype.l1", raster = FALSE, label = TRUE)

#Complex Prediction

p1 <- DimPlot(total.mtx_PAH, reduction = "umap", group.by = "predicted.celltype.l2", raster = FALSE, label = TRUE)

p2 <- DimPlot(total.mtx_Control, reduction = "umap", group.by = "predicted.celltype.l2", raster = FALSE, label = TRUE)


#Overlapped on Same Graph 

DimPlot(total.mtx_harmony, reduction = "umap", group.by = "orig.ident", label = TRUE, raster = FALSE) 

#Checking Confidence in Predictions

table(Idents(total.mtx_harmony), total.mtx_harmony$predicted.celltype.l1.score < 0.5) #Most clusters are confident predictions, plasma cells and erythrocytes are a bit uncertain

Idents(total.mtx_harmony) <- "predicted.celltype.l2"
table(Idents(total.mtx_harmony), total.mtx_harmony$predicted.celltype.l2.score < 0.5) #Heterogeneity for T/NK cells and myeloid cells, rest are homogeneous
table(Idents(total.mtx_harmony), total.mtx_harmony$predicted.celltype.l2) #Heterogeneity for T/NK cells and myeloid cells, rest are homogeneous


total.mtx_harmony$condition <- total.mtx_harmony$orig.ident

#Set Resolution of Your Choice as Seurat Cluster Numbers and Identity as Seurat Clusters

total.mtx_harmony$seurat_clusters <- total.mtx_harmony$RNA_snn_res.0.7 #Resolution = 0.7
Idents(total.mtx_harmony) <- total.mtx_harmony$seurat_clusters

### FindConservedMarkers() / Begin Annotation for this Dataset

#Hepatocytes

h1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 0, only.pos = TRUE, grouping.var = 'orig.ident')
h1_genes <- rownames(head(h1_markers,20))
FeaturePlot(total.mtx_harmony, features = h1_genes, min.cutoff = "q10")

h2_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 3, only.pos = TRUE, grouping.var = 'orig.ident')
h2_genes <- rownames(head(h2_markers,20))
FeaturePlot(total.mtx_harmony, features = h2_genes, min.cutoff = "q10")

h3_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 4, only.pos = TRUE, grouping.var = 'orig.ident')
h3_genes <- rownames(head(h3_markers,20))
FeaturePlot(total.mtx_harmony, features = h3_genes, min.cutoff = "q10")

h4_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 5, only.pos = TRUE, grouping.var = 'orig.ident')
h4_genes <- rownames(head(h4_markers,20))
FeaturePlot(total.mtx_harmony, features = h4_genes, min.cutoff = "q10")

h5_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 6, only.pos = TRUE, grouping.var = 'orig.ident')
h5_genes <- rownames(head(h5_markers,20))
FeaturePlot(total.mtx_harmony, features = h5_genes, min.cutoff = "q10")

h6_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 7, only.pos = TRUE, grouping.var = 'orig.ident')
h6_genes <- rownames(head(h6_markers,20))
FeaturePlot(total.mtx_harmony, features = h6_genes, min.cutoff = "q10")

h7_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 12, only.pos = TRUE, grouping.var = 'orig.ident')
h7_genes <- rownames(head(h7_markers,20))
FeaturePlot(total.mtx_harmony, features = h7_genes, min.cutoff = "q10")

h8_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 13, only.pos = TRUE, grouping.var = 'orig.ident') #Not present in PAH, only Control
h8_genes <- rownames(head(h8_markers,20))
FeaturePlot(total.mtx_harmony, features = h8_genes, min.cutoff = "q10")

h9_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 16, only.pos = TRUE, grouping.var = 'orig.ident') #Not present in PAH, only Control
h9_genes <- rownames(head(h9_markers,20))
FeaturePlot(total.mtx_harmony, features = h9_genes, min.cutoff = "q10")

h10_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 17, only.pos = TRUE, grouping.var = 'orig.ident')
h10_genes <- rownames(head(h10_markers,20))
FeaturePlot(total.mtx_harmony, features = h10_genes, min.cutoff = "q10")

#Endothelial

e1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 1, only.pos = TRUE, grouping.var = 'orig.ident')
e1_genes <- rownames(head(e1_markers,20))
FeaturePlot(total.mtx_harmony, features = e1_genes, min.cutoff = "q20")

e2_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 15, only.pos = TRUE, grouping.var = 'orig.ident')
e2_genes <- rownames(head(e2_markers,20))
FeaturePlot(total.mtx_harmony, features = e2_genes, min.cutoff = "q10")

e3_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 23, only.pos = TRUE, grouping.var = 'orig.ident')
e3_genes <- rownames(head(e3_markers,20))
FeaturePlot(total.mtx_harmony, features = e3_genes, min.cutoff = "q10")

#Myeloid

m1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 9, only.pos = TRUE, grouping.var = 'orig.ident')
m1_genes <- rownames(head(m1_markers,24))
FeaturePlot(total.mtx_harmony, features = m1_genes, min.cutoff = "q10")

#Cholangiocyte

c1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 14, only.pos = TRUE, grouping.var = 'orig.ident')
c1_genes <- rownames(head(c1_markers,20))
FeaturePlot(total.mtx_harmony, features = c1_genes, min.cutoff = "q10")


#T/NK

t1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 8, only.pos = TRUE, grouping.var = 'orig.ident')
t1_genes <- rownames(head(t1_markers,20))
FeaturePlot(total.mtx_harmony, features = t1_genes, min.cutoff = "q10")

t2_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 20, only.pos = TRUE, grouping.var = 'orig.ident')
t2_genes <- rownames(head(t2_markers,20))
FeaturePlot(total.mtx_harmony, features = t2_genes, min.cutoff = "q10")

t3_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 21, only.pos = TRUE, grouping.var = 'orig.ident')
t3_genes <- rownames(head(t3_markers,20))
FeaturePlot(total.mtx_harmony, features = t3_genes, min.cutoff = "q10")

#Erythrocyte/Plasma 

ep1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 18, only.pos = TRUE, grouping.var = 'orig.ident')
ep1_genes <- rownames(head(ep1_markers,20))
FeaturePlot(total.mtx_harmony, features = ep1_genes, min.cutoff = "q10")

#B Cells

b1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 19, only.pos = TRUE, grouping.var = 'orig.ident')
b1_genes <- rownames(head(b1_markers,20))
FeaturePlot(total.mtx_harmony, features = b1_genes, min.cutoff = "q10")

#HSC/mFB

hm1_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 2, only.pos = TRUE, grouping.var = 'orig.ident')
hm1_genes <- rownames(head(hm1_markers,24))
FeaturePlot(total.mtx_harmony, features = hm1_genes, min.cutoff = "q10")

hm2_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 10, only.pos = TRUE, grouping.var = 'orig.ident')
hm2_genes <- rownames(head(hm2_markers,24))
FeaturePlot(total.mtx_harmony, features = hm2_genes, min.cutoff = "q10")

hm3_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 11, only.pos = TRUE, grouping.var = 'orig.ident')
hm3_genes <- rownames(head(hm3_markers,24))
FeaturePlot(total.mtx_harmony, features = hm3_genes, min.cutoff = "q10")

hm4_markers <- FindConservedMarkers(total.mtx_harmony, ident.1 = 22, only.pos = TRUE, grouping.var = 'orig.ident')
hm4_genes <- rownames(head(hm4_markers,20))
FeaturePlot(total.mtx_harmony, features = hm4_genes, min.cutoff = "q10")


# Marker Genes - Dot Plots

#Hepatocytes

DotPlot(total.mtx_harmony, features = c('SLC39A14', 'CP', 'CPS1', 'SAA2', 'NOS1AP', 'CFHR4', 'SLC13A5', 'C4BPB', 'C8B', 'ABCC2', 'PCK1', 'C5','AL136456.1','GPAM'), cols = c("blue", "red"))

DotPlot(total.mtx_harmony, features = c('ALDOB', 'AL391117.1', 'APOB', 'ERRFI1', 'FGL1', 'C3', 'TF', 'SDS', 'CYP3A5', 'SLC7A2', 'SLCO1B3', 'ADH4', 'SLC38A4', 'CFH'), cols = c("blue", "red"))

DotPlot(total.mtx_harmony, features = c('ELOVL6', 'ALDH1A2', 'GPAM', 'ACSM2B', 'EBNA1BP2', 'TENM2', 'AC007262.2', 'CPB2', 'BHMT', 'PLG', 'CYP2C9', 'ACSM2A',"ADH4"), cols = c("blue", "red"))


#Endothelial Cells
DotPlot(total.mtx_cluster, features = c('PTPRB', 'ST6GALNAC3', 'NRG3', 'STAB2', 'BMPER', 'ROBO2', 'FLT1', 'LDB2', 'EGFL7'), cols = c("blue", "red"))

#Myeloid Cells
DotPlot(total.mtx_harmony, features = c('GRK3', 'DMXL2','TBXAS1','CD163','EPB41L3','MCTP1','PDE4B','SAT1','CTSS'), cols = c("blue", "red"))

#Cholangiocytes
DotPlot(total.mtx_harmony, features = c('CTNND2','ANXA4', "CASC15", 'PKHD1', 'BICC1', 'DCDC2', 'RALYL', 'FGFR2','CFTR','SLC12A2'), cols = c("blue", "red"))

#T/NK Cells
DotPlot(total.mtx_harmony, features = c('PARP8', 'SKAP1', 'PRKCH', 'FYN', 'STAT4', 'CD247', 'MCTP2','CD96'), cols = c("blue", "red"))

#Erythrocyte/Plasma Cells
DotPlot(total.mtx_harmony, features = c('CA1', 'PIP5K1B', 'ANK1'), cols = c("blue", "red"))

#B Cells

DotPlot(total.mtx_harmony, features = c('BANK1','BLK','ADAM28','MS4A1','FCRL5',"WDFY4",'PAX5','EBF1'), cols = c("blue", "red"))

#HSC/mFB
DotPlot(total.mtx_harmony, features = c('ITGA9', 'CCBE1', 'ZFPM2', 'ADAMTSL1','ADAMTS12','PDGFRA','COL25A1'), cols = c("blue", "red"))


#Manual Annotation of Clusters and Removing Noise

Idents(total.mtx_harmony) <- total.mtx_harmony$seurat_clusters #Reassinging Numerical Labels for Renaming

total.mtx_harmony <- RenameIdents(total.mtx_harmony, "0" = "Hepatocyte 1")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "1" = "Endothelial 1")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "2" = "HSC/mFB 1")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "3" = "Hepatocyte 3")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "4" = "Hepatocyte 6")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "5" = "Hepatocyte 7")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "6" = "Hepatocyte 5")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "7" = "Hepatocyte 2")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "8" = "Lymphocyte 1")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "9" = "Macrophage")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "10" = "HSC/mFB 2")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "11" = "HSC/mFB 3")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "12" = "Hepatocyte 4")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "13" = "Hepatocyte 9")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "14" = "Cholangiocyte")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "15" = 'Endothelial 2')
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "16" = "Hepatocyte 10")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "17" = 'Hepatocyte 8')
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "18" = "Plasma/Erythrocyte")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "19" = "B Cell")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "20" = "Lymphocyte 2")
total.mtx_harmony <- RenameIdents(total.mtx_harmony, "21" = "Lymphocyte 3")

total.mtx_harmony$subcelltype <- Idents(total.mtx_harmony) 

#Creating Subset Without Junk Clusters

total.mtx_cluster <- subset(total.mtx_harmony, idents = c('22','23'), invert = TRUE) #Excluding Junk Clusters


# Subsetting by Cell Type
HP_list <- c('Hepatocyte 1', 'Hepatocyte 2', 'Hepatocyte 3', 'Hepatocyte 4', 'Hepatocyte 5',
             'Hepatocyte 6', 'Hepatocyte 7', 'Hepatocyte 8', 'Hepatocyte 9', 'Hepatocyte 10')

HP_cluster.mtx <- subset(total.mtx_cluster, subcelltype %in% HP_list)
Idents(HP_cluster.mtx) <- HP_cluster.mtx$subcelltype

DimPlot(HP_cluster.mtx, reduction = "umap", group.by = "subcelltype", raster = FALSE, label = TRUE)
DimPlot(HP_cluster.mtx, reduction = "umap", group.by = "condition", raster = FALSE, label = TRUE)

#Marker Genes for Hepatocyte Clusters

#PAH Clusters

HP1_markers <- FindConservedMarkers(HP_cluster.mtx, ident.1 = "Hepatocyte 1", only.pos = TRUE, grouping.var = "condition")
HP1_markers_list <- as.data.frame(HP1_markers)
HP1_markers_list$genes <- rownames(HP1_markers_list) 
write_csv(HP1_markers_list, file = "~/Hepatocyte Cluster 1 Markers - MB250225.csv")

HP2_markers <- FindConservedMarkers(HP_cluster.mtx, ident.1 = "Hepatocyte 2", only.pos = TRUE, grouping.var = "condition")
HP2_markers_list <- as.data.frame(HP2_markers)
HP2_markers_list$genes <- rownames(HP2_markers_list) 
write_csv(HP2_markers_list, file = "~/Hepatocyte Cluster 2 Markers - MB250225.csv")

HP3_markers <- FindConservedMarkers(HP_cluster.mtx, ident.1 = "Hepatocyte 3", only.pos = TRUE, grouping.var = "condition")
HP3_markers_list <- as.data.frame(HP3_markers)
HP3_markers_list$genes <- rownames(HP3_markers_list) 
write_csv(HP3_markers_list, file = "~/Hepatocyte Cluster 3 Markers - MB250225.csv")

HP4_markers <- FindConservedMarkers(HP_cluster.mtx, ident.1 = "Hepatocyte 4", only.pos = TRUE, grouping.var = "condition")
HP4_markers_list <- as.data.frame(HP4_markers)
HP4_markers_list$genes <- rownames(HP4_markers_list) 
write_csv(HP4_markers_list, file = "~/Hepatocyte Cluster 4 Markers - MB250225.csv")

HP7_markers <- FindConservedMarkers(HP_cluster.mtx, ident.1 = "Hepatocyte 7", only.pos = TRUE, grouping.var = "condition")
HP7_markers_list <- as.data.frame(HP7_markers)
HP7_markers_list$genes <- rownames(HP7_markers_list) 
write_csv(HP7_markers_list, file = "~/Hepatocyte Cluster 7 Markers - MB250225.csv")

#Control Clusters

HP5_markers <- FindConservedMarkers(HP_cluster.mtx, ident.1 = "Hepatocyte 5", only.pos = TRUE, grouping.var = "condition")
HP5_markers_list <- as.data.frame(HP5_markers)
HP5_markers_list$genes <- rownames(HP5_markers_list) 
write_csv(HP5_markers_list, file = "~/Hepatocyte Cluster 5 Markers - MB250225.csv")

HP6_markers <- FindConservedMarkers(HP_cluster.mtx, ident.1 = "Hepatocyte 6", only.pos = TRUE, grouping.var = "condition")
HP6_markers_list <- as.data.frame(HP6_markers)
HP6_markers_list$genes <- rownames(HP6_markers_list) 
write_csv(HP6_markers_list, file = "~/Hepatocyte Cluster 6 Markers - MB250225.csv")

HP8_markers <- FindConservedMarkers(HP_cluster.mtx, ident.1 = "Hepatocyte 8", only.pos = TRUE, grouping.var = "condition")
HP8_markers_list <- as.data.frame(HP8_markers)
HP8_markers_list$genes <- rownames(HP8_markers_list) 
write_csv(HP8_markers_list, file = "~/Hepatocyte Cluster 8 Markers - MB250225.csv")

HP9_markers <- FindConservedMarkers(HP_cluster.mtx, ident.1 = "Hepatocyte 9", only.pos = TRUE, grouping.var = "condition")
HP9_markers_list <- as.data.frame(HP9_markers)
HP9_markers_list$genes <- rownames(HP9_markers_list) 
write_csv(HP9_markers_list, file = "~/Hepatocyte Cluster 9 Markers - MB250225.csv")

HP10_markers <- FindConservedMarkers(HP_cluster.mtx, ident.1 = "Hepatocyte 10", only.pos = TRUE, grouping.var = "condition")
HP10_markers_list <- as.data.frame(HP10_markers)
HP10_markers_list$genes <- rownames(HP10_markers_list) 
write_csv(HP10_markers_list, file = "~/Hepatocyte Cluster 10 Markers - MB250225.csv")






EC_list <- c('Endothelial 1', 'Endothelial 2')

EC_cluster.mtx <- subset(total.mtx_cluster, cell_identity %in% EC_list)

HSC_list <- c('HSC/mFB 1', 'HSC/mFB 2','HSC/mFB 3')

HSC_cluster.mtx <- subset(total.mtx_cluster, cell_identity %in% HSC_list)


DimPlot(total.mtx_cluster, reduction = "umap", group.by = "cell_identity", label = TRUE, raster = FALSE) 

#Creating General Cell Labeled Clusters (No Cluster #)

total.mtx_cluster <- RenameIdents(total.mtx_cluster, "Hepatocyte 1" = 'Hepatocyte')
total.mtx_cluster <- RenameIdents(total.mtx_cluster, "Endothelial 1" = 'Endothelial')
total.mtx_cluster <- RenameIdents(total.mtx_cluster, "HSC/mFB 1" = 'HSC/mFB')
total.mtx_cluster <- RenameIdents(total.mtx_cluster, "Hepatocyte 3" = 'Hepatocyte')
total.mtx_cluster <- RenameIdents(total.mtx_cluster, "Hepatocyte 6" = 'Hepatocyte')
total.mtx_cluster <- RenameIdents(total.mtx_cluster, "Hepatocyte 7" = 'Hepatocyte')
total.mtx_cluster <- RenameIdents(total.mtx_cluster,  "Hepatocyte 5" = 'Hepatocyte')
total.mtx_cluster <- RenameIdents(total.mtx_cluster,  "Hepatocyte 2" = 'Hepatocyte')
total.mtx_cluster <- RenameIdents(total.mtx_cluster, "Lymphocyte 1" = 'Lymphocyte')
total.mtx_cluster <- RenameIdents(total.mtx_cluster,  "Macrophage" = 'Macrophage')
total.mtx_cluster <- RenameIdents(total.mtx_cluster, "HSC/mFB 2" = 'HSC/mFB')
total.mtx_cluster <- RenameIdents(total.mtx_cluster,  "HSC/mFB 3" = 'HSC/mFB')
total.mtx_cluster <- RenameIdents(total.mtx_cluster,  "Hepatocyte 4" = 'Hepatocyte')
total.mtx_cluster <- RenameIdents(total.mtx_cluster,  "Hepatocyte 9" = 'Hepatocyte')
total.mtx_cluster <- RenameIdents(total.mtx_cluster, "Cholangiocyte" = 'Cholangiocyte')
total.mtx_cluster <- RenameIdents(total.mtx_cluster,  'Endothelial 2' = 'Endothelial')
total.mtx_cluster <- RenameIdents(total.mtx_cluster,  "Hepatocyte 10" = 'Hepatocyte')
total.mtx_cluster <- RenameIdents(total.mtx_cluster,  'Hepatocyte 8' = 'Hepatocyte')
total.mtx_cluster <- RenameIdents(total.mtx_cluster,  "Plasma/Erythrocyte" = 'Plasma')
total.mtx_cluster <- RenameIdents(total.mtx_cluster,  "B Cell" = 'B Cell')
total.mtx_cluster <- RenameIdents(total.mtx_cluster,  "Lymphocyte 2" = 'Lymphocyte')
total.mtx_cluster <- RenameIdents(total.mtx_cluster,  "Lymphocyte 3" = 'Lymphocyte')

total.mtx_cluster$cell_identity <- Idents(total.mtx_cluster) 

#Creating Control and PAH Separate UMAPs

DimPlot(total.mtx_cluster, reduction = "umap", group.by = "cell_identity", split.by ='condition',label = TRUE, raster = FALSE) 

#Proportion Bar Plot for Bottom 

cell_table <- table(total.mtx_cluster$condition, total.mtx_cluster$subcelltype)
cell_data <- as.data.frame(cell_table)
write_csv(cell_data, file ='Cell Subcelltype Counts - MB250303.csv')

sample_table <- as.data.frame(table(total.mtx_cluster$sample, total.mtx_cluster$subcelltype))
write_csv(sample_table, file = "Subcelltype Sample Data - MB250303.csv")

sample_table_main <- as.data.frame(table(total.mtx_cluster$sample, total.mtx_cluster$cell_identity))
write_csv(sample_table_main, file = "Main Celltype Counts - MB250303")

ggplot(cell_data, aes(x = Var1, fill = Var2))+
  geom_bar(method = "fill")


props <- read.csv(file= "/users/1/blake561/Cell Prop Data MB250123.csv")


ggplot(props, aes(x = Condition, y = Proportion, fill = Cell.Type))+
  geom_col(position = "fill", width = 0.4)+
  scale_fill_manual(values = c('Hepatocyte' = 'orchid', 'Endothelial' = 'slateblue1', 
                               'HSC/mFB' = 'deeppink1', 'Macrophage' = 'chartreuse3', 
                               'Lymphocyte' = 'dodgerblue', 'B Cell' = 'yellow2', 
                               'Cholangiocyte' = 'salmon', 'Plasma' = 'darkorange1'))



#Save File for DEG Analysis
total.mtx_cluster$sample_num <- paste0(total.mtx_cluster$cell_identity, total.mtx_cluster$sample)

saveRDS(total.mtx_cluster, file = "total.mtx_final")

total.mtx_filtered <- readRDS(file = "/users/1/blake561/total.mtx_final")


#Feature plots - see specific genes from expression based on the graphs
FeaturePlot(total.mtx_final, features = ('STAT3'), min.cutoff = "q10", raster = FALSE)
FeaturePlot(PAH, features = ('IL6'), min.cutoff = "q10", raster = FALSE) # NO BUENO
FeaturePlot(total.mtx_final, features = c('INHBA','INHBB'), min.cutoff = "q10", raster = FALSE)
FeaturePlot(total.mtx_final, features = ('HIF1A'), min.cutoff = "q10", raster = FALSE) # ENDO - HIGH
FeaturePlot(total.mtx_final, features = ('GDF15'), min.cutoff = "q10", raster = FALSE) # ENDO - HIGH
FeaturePlot(total.mtx_final, features = ('BMP10'), min.cutoff = "q10", raster = FALSE) # ENDO - HIGH

#Further Separating Macrophage 

# Check the elbow plot for your macrophage subset
ElbowPlot(mac_subset, ndims = 15)

# Look at PC loadings and variance explained
print(mac_subset[["pca"]], dims = 1:20, nfeatures = 5)

# Check how much variance is explained by each PC
pca_var <- mac_subset[["pca"]]@stdev^2
pca_var_percent <- round(pca_var / sum(pca_var) * 100, 2)
print(pca_var_percent[1:20])

# Compare with 1:12 if you want to be slightly more inclusive
mac_subset <- FindNeighbors(mac_subset, reduction = "pca", dims = 1:12)
mac_subset <- FindClusters(mac_subset, resolution = 0.05)
mac_subset <- RunUMAP(mac_subset, reduction = "pca", dims = 1:12)

# Visualize the result
DimPlot(mac_subset, group.by = "seurat_clusters", label = TRUE)
DimPlot(mac_subset, group.by = "seurat_clusters", split.by = 'condition', label = TRUE)

#Resident vs Recruited Macrophage Marker Genes#

Idents(mac_subset) <- "RNA_snn_res.0.02"

# Find markers for all clusters
mac_markers <- FindAllMarkers(mac_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For resident macrophage markers
FeaturePlot(mac_subset, features = c("VSIG4", "TIMD4", "MARCO", "CD163", "MRC1"))

# For recruited macrophage markers
FeaturePlot(mac_subset, features = c('CCR2', 'S100A8', 'TNF','TREM1','S100A4'))

mac_subset <- AddModuleScore(mac_subset, 
                             features = list(c("VSIG4", "TIMD4", "MARCO", "CD163", "MRC1")), 
                             name = "Resident_Score")

mac_subset <- AddModuleScore(mac_subset, 
                             features = list(c('CCR2', 'S100A8','S100A4')), 
                             name = "Recruited_Score")


# Visualize the scores
FeaturePlot(mac_subset, features = c("Resident_Score1", "Recruited_Score1"), ncol = 2)

# Use the normalized data slot instead of scaled data (cluster size influence)
DoHeatmap(mac_subset, 
          features = c("VSIG4", "TIMD4", "MARCO", "CD163", "MRC1",'CCR2', 'S100A8', 'TNF','TREM1','S100A4'),
          group.by = "seurat_clusters",
          slot = "data") + 
  scale_fill_viridis_c()

#Reassigning Cluster Identities 

mac_subset <- RenameIdents(mac_subset,
                           "0" = "Recruited Macrophage",    
                           "1" = "Resident Kupffer",  
                           "2" = "Resident Kupffer"
)

# Save the annotations
mac_subset$macrophage_subtype <- Idents(mac_subset)

# Split by condition and create separate plots
control_mac <- subset(mac_subset, subset = orig.ident == "Control")
pah_mac <- subset(mac_subset, subset = orig.ident == "PAH")

# Create UMAPs
control_umap <- DimPlot(control_mac, group.by = "macrophage_subtype", 
                        label = TRUE, raster = FALSE) +
  ggtitle("Control - Macrophage Subtypes") +
  theme_classic()

pah_umap <- DimPlot(pah_mac, group.by = "macrophage_subtype", 
                    label = TRUE, raster = FALSE) +
  ggtitle("PAH - Macrophage Subtypes") +
  theme_classic()

print(control_umap)
print(pah_umap)



#Save PAH Macrophage Subset 

saveRDS(mac_subset, file = "PAH Macrophage Subset")

mac_subset <- readRDS(file = "/users/1/blake561/PAH Macrophage Subset")

#Markers for Clusters 

Idents(mac_subset) <- "RNA_snn_res.0.05"

basic_markers <- FindAllMarkers(mac_subset,
                                min.pct = 0.05,
                                logfc.threshold = 0.25,
                                only.pos = TRUE)

print(paste("Basic markers found:", nrow(basic_markers)))
if(nrow(basic_markers) > 0) {
  print(head(basic_markers, 20))
}

write.csv(basic_markers, "PAH_macrophage_basic_markers.csv", row.names = FALSE)