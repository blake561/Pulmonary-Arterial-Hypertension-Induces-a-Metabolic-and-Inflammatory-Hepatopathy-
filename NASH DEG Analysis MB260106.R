##Load Packages
library(Seurat)
library(ggplot2)
library(EnhancedVolcano)
library(tidyverse)

##Load Data
NASH <- readRDS("/users/1/blake561/NASH_total.mtx_final")  

# Create condition column
total.mtx_deg$condition <- ifelse(grepl("^CN", total.mtx_deg$sample), "Control", "NASH")

# JoinLayers if using Seurat v5
total.mtx_deg <- JoinLayers(total.mtx_deg)

# Confirm subcelltype is set
Idents(total.mtx_deg) <- total.mtx_deg$subcelltype
table(total.mtx_deg$subcelltype, total.mtx_deg$condition)

##=============================================================================
## Subset by Cell Type
##=============================================================================

Hepatocyte <- subset(total.mtx_deg, subset = subcelltype == "Hepatocyte")
Endothelial <- subset(total.mtx_deg, subset = subcelltype == "Endothelial Cell")
HSC <- subset(total.mtx_deg, subset = subcelltype == "HSC/mFB")
Macrophage <- subset(total.mtx_deg, subset = subcelltype == "Macrophage")
Cholangiocyte <- subset(total.mtx_deg, subset = subcelltype == "Cholangiocyte")
Lymphocyte <- subset(total.mtx_deg, subset = subcelltype == "Lymphocyte")

##=============================================================================
## UMAPs by Cell Type
##=============================================================================

DimPlot(Hepatocyte, reduction = "umap", split.by = "condition", cols = c('Hepatocyte' = 'chartreuse3'), raster = FALSE) + ggtitle("Hepatocyte")
DimPlot(Endothelial, reduction = "umap", split.by = "condition", cols = c("Endothelial Cell" = 'deepskyblue'), raster = FALSE) + ggtitle("Endothelial Cell")
DimPlot(HSC, reduction = "umap", split.by = "condition", cols = c('HSC/mFB' = 'slateblue1'), raster = FALSE) + ggtitle("HSC/mFB")
DimPlot(Macrophage, reduction = "umap", split.by = "condition", cols = c('Macrophage' = 'violet'), raster = FALSE) + ggtitle("Macrophage")
DimPlot(Cholangiocyte, reduction = "umap", split.by = "condition", cols = c('Cholangiocyte' = 'orange'), raster = FALSE) + ggtitle("Cholangiocyte")
DimPlot(Lymphocyte, reduction = "umap", split.by = "condition", cols = c('Lymphocyte' = 'firebrick'), raster = FALSE) + ggtitle("Lymphocyte")

##=============================================================================
## DEG Analysis
##=============================================================================

# Hepatocyte
Idents(Hepatocyte) <- Hepatocyte$condition
hep_deg <- FindMarkers(
  object = Hepatocyte,
  ident.1 = "NASH",
  ident.2 = "Control",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25
)
hep_deg <- hep_deg[
  (abs(hep_deg$avg_log2FC) >= 0.25) &
    (hep_deg$p_val_adj < 0.05),
]
write.csv(hep_deg, "hep_NASH_FC_0.25.csv", row.names = TRUE)

# Endothelial Cell
Idents(Endothelial) <- Endothelial$condition
endo_deg <- FindMarkers(
  object = Endothelial,
  ident.1 = "NASH",
  ident.2 = "Control",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25
)
endo_deg <- endo_deg[
  (abs(endo_deg$avg_log2FC) >= 0.25) &
    (endo_deg$p_val_adj < 0.05),
]
write.csv(endo_deg, "endo_NASH_FC_0.25.csv", row.names = TRUE)

# HSC/mFB
Idents(HSC) <- HSC$condition
hsc_deg <- FindMarkers(
  object = HSC,
  ident.1 = "NASH",
  ident.2 = "Control",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25
)
hsc_deg <- hsc_deg[
  (abs(hsc_deg$avg_log2FC) >= 0.25) &
    (hsc_deg$p_val_adj < 0.05),
]
write.csv(hsc_deg, "hsc_NASH_FC_0.25.csv", row.names = TRUE)

# Macrophage
Idents(Macrophage) <- Macrophage$condition
mac_deg <- FindMarkers(
  object = Macrophage,
  ident.1 = "NASH",
  ident.2 = "Control",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25
)
mac_deg <- mac_deg[
  (abs(mac_deg$avg_log2FC) >= 0.25) &
    (mac_deg$p_val_adj < 0.05),
]
write.csv(mac_deg, "mac_NASH_FC_0.25.csv", row.names = TRUE)

# Cholangiocyte
Idents(Cholangiocyte) <- Cholangiocyte$condition
chol_deg <- FindMarkers(
  object = Cholangiocyte,
  ident.1 = "NASH",
  ident.2 = "Control",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25
)
chol_deg <- chol_deg[
  (abs(chol_deg$avg_log2FC) >= 0.25) &
    (chol_deg$p_val_adj < 0.05),
]
write.csv(chol_deg, "chol_NASH_FC_0.25.csv", row.names = TRUE)

# Lymphocyte
Idents(Lymphocyte) <- Lymphocyte$condition
lymph_deg <- FindMarkers(
  object = Lymphocyte,
  ident.1 = "NASH",
  ident.2 = "Control",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25
)
lymph_deg <- lymph_deg[
  (abs(lymph_deg$avg_log2FC) >= 0.25) &
    (lymph_deg$p_val_adj < 0.05),
]
write.csv(lymph_deg, "lymph_NASH_FC_0.25.csv", row.names = TRUE)

##=============================================================================
## Summary
##=============================================================================

cat("DEG Summary:\n")
cat("Hepatocyte:", nrow(hep_deg), "DEGs\n")
cat("Endothelial:", nrow(endo_deg), "DEGs\n")
cat("HSC/mFB:", nrow(hsc_deg), "DEGs\n")
cat("Macrophage:", nrow(mac_deg), "DEGs\n")
cat("Cholangiocyte:", nrow(chol_deg), "DEGs\n")
cat("Lymphocyte:", nrow(lymph_deg), "DEGs\n")