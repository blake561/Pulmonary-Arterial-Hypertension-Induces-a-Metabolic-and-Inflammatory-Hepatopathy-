# Load Packages
library(BiocManager)
library(tidyverse)
library(Seurat)
library(gridExtra)
library(Azimuth)
library(harmony)
library(multtest)
library(metap)
library(ggplot2)
library(SeuratData)
library(EnhancedVolcano)
library(DESeq2)
library(tidyr)
library(rlang)
library(devtools)
library(NMF)
library(circlize)
library(ComplexHeatmap)
library(CellChat)
library(ggalluvial)
library(matrixStats)

#Load File

PAH <- readRDS(file = '/projects/standard/prin0088/shared/Madelyn Liver and Kidney Projects/PAH_total.mtx_final')
NASH <- readRDS("/users/1/blake561/NASH_total.mtx_final")  
Fontan <- readRDS(file = '/projects/standard/prin0088/shared/Madelyn Liver and Kidney Projects/Fontan_total.mtx_final')

#Separating By Cell Type for DEG Analysis


#UMAPs for Figures

#PAH

#Hepatocyte

# Create a new metadata column that combines cell identity and condition (for unique color mapping)
Hepatocyte$ColorByCondition <- Hepatocyte$condition

# Custom colors to map per panel
condition_colors <- c("Control" = "salmon", "PAH" = "steelblue")

# Now plot with split.by, and color based on your new metadata column
DimPlot(
  Hepatocyte,
  reduction = "umap",
  split.by  = "condition",
  group.by  = "ColorByCondition",  # This applies color within each split
  cols      = condition_colors,
  raster    = FALSE
)

#Endothelial 

# Create a new metadata column that combines cell identity and condition (for unique color mapping)
Endothelial$ColorByCondition <- Endothelial$condition

# Custom colors to map per panel
condition_colors <- c("Control" = "salmon", "PAH" = "steelblue")

# Now plot with split.by, and color based on your new metadata column
DimPlot(
  Endothelial,
  reduction = "umap",
  split.by  = "condition",
  group.by  = "ColorByCondition",  # This applies color within each split
  cols      = condition_colors,
  raster    = FALSE
)

#HSC

# Create a new metadata column that combines cell identity and condition (for unique color mapping)
HSC$ColorByCondition <- HSC$condition

# Custom colors to map per panel
condition_colors <- c("Control" = "salmon", "PAH" = "steelblue")

# Now plot with split.by, and color based on your new metadata column
DimPlot(
  HSC,
  reduction = "umap",
  split.by  = "condition",
  group.by  = "ColorByCondition",  # This applies color within each split
  cols      = condition_colors,
  raster    = FALSE
)

#Macrophage

# Create a new metadata column that combines cell identity and condition (for unique color mapping)
Macrophage$ColorByCondition <- Macrophage$condition

# Custom colors to map per panel
condition_colors <- c("Control" = "salmon", "PAH" = "steelblue")

# Now plot with split.by, and color based on your new metadata column
DimPlot(
  Macrophage,
  reduction = "umap",
  split.by  = "condition",
  group.by  = "ColorByCondition",  # This applies color within each split
  cols      = condition_colors,
  raster    = FALSE
)

#NASH

# ── NASH: Hepatocyte ───────────────────────────────────────────────

NASH_Hepatocyte$ColorByCondition <- NASH_Hepatocyte$condition

# Custom colors to map per panel
condition_colors <- c("Control" = "darkgoldenrod2", "NASH" = "orchid4")

# Now plot with split.by, and color based on your new metadata column
DimPlot(
  NASH_Hepatocyte,
  reduction = "umap",
  split.by  = "condition",
  group.by  = "ColorByCondition",  # This applies color within each split
  cols      = condition_colors,
  raster    = FALSE
)

# ── NASH: Endothelial ───────────────────────────────────────────────

NASH_Endothelial$ColorByCondition <- NASH_Endothelial$condition

# Custom colors to map per panel
condition_colors <- c("Control" = "darkgoldenrod2", "NASH" = "orchid4")

# Now plot with split.by, and color based on your new metadata column
DimPlot(
  NASH_Endothelial,
  reduction = "umap",
  split.by  = "condition",
  group.by  = "ColorByCondition",  # This applies color within each split
  cols      = condition_colors,
  raster    = FALSE
)

# ── NASH: HSC ───────────────────────────────────────────────────────

NASH_HSC$ColorByCondition <- NASH_HSC$condition

# Custom colors to map per panel
condition_colors <- c("Control" = "darkgoldenrod2", "NASH" = "orchid4")

# Now plot with split.by, and color based on your new metadata column
DimPlot(
  NASH_HSC,
  reduction = "umap",
  split.by  = "condition",
  group.by  = "ColorByCondition",  # This applies color within each split
  cols      = condition_colors,
  raster    = FALSE
)

# ── NASH: Macrophage ────────────────────────────────────────────────

NASH_Macrophage$ColorByCondition <- NASH_Macrophage$condition

# Custom colors to map per panel
condition_colors <- c("Control" = "darkgoldenrod2", "NASH" = "orchid4")

# Now plot with split.by, and color based on your new metadata column
DimPlot(
  NASH_Macrophage,
  reduction = "umap",
  split.by  = "condition",
  group.by  = "ColorByCondition",  # This applies color within each split
  cols      = condition_colors,
  raster    = FALSE
)

#Fontan 

# ── Fontan: Hepatocyte ───────────────────────────────────────────────

Fontan_Hepatocyte$ColorByCondition <- Fontan_Hepatocyte$orig.ident

# Custom colors to map per panel
condition_colors <- c("Fontan Control" = "sienna", "Fontan" = "olivedrab")

# Now plot with split.by, and color based on your new metadata column
DimPlot(
  Fontan_Hepatocyte,
  reduction = "umap",
  split.by  = "orig.ident",
  group.by  = "ColorByCondition",  # This applies color within each split
  cols      = condition_colors,
  raster    = FALSE
)

# ── Fontan: Endothelial ───────────────────────────────────────────────

Fontan_Endothelial$ColorByorig.ident <- Fontan_Endothelial$orig.ident

# Custom colors to map per panel
orig.ident_colors <- c("Fontan Control" = "sienna", "Fontan" = "olivedrab")

# Now plot with split.by, and color based on your new metadata column
DimPlot(
  Fontan_Endothelial,
  reduction = "umap",
  split.by  = "orig.ident",
  group.by  = "ColorByorig.ident",  # This applies color within each split
  cols      = orig.ident_colors,
  raster    = FALSE
)

# ── Fontan: HSC ───────────────────────────────────────────────────────

Fontan_HSC$ColorByorig.ident <- Fontan_HSC$orig.ident

# Custom colors to map per panel
orig.ident_colors <- c("Fontan Control" = "sienna", "Fontan" = "olivedrab")

# Now plot with split.by, and color based on your new metadata column
DimPlot(
  Fontan_HSC,
  reduction = "umap",
  split.by  = "orig.ident",
  group.by  = "ColorByorig.ident",  # This applies color within each split
  cols      = orig.ident_colors,
  raster    = FALSE
)

# ── Fontan: Macrophage ────────────────────────────────────────────────

Fontan_Macrophage$ColorByorig.ident <- Fontan_Macrophage$orig.ident

# Custom colors to map per panel
orig.ident_colors <- c("Fontan Control" = "sienna", "Fontan" = "olivedrab")

# Now plot with split.by, and color based on your new metadata column
DimPlot(
  Fontan_Macrophage,
  reduction = "umap",
  split.by  = "orig.ident",
  group.by  = "ColorByorig.ident",  # This applies color within each split
  cols      = orig.ident_colors,
  raster    = FALSE
)

#Proportions 

##Cell Type Proportions##

cell_table <- table(Fontan$orig.ident, Fontan$subcelltype)
cell_data <- as.data.frame(cell_table)
write_csv(cell_data, file ='Fontan Subcelltype Counts - MB250422.csv')

Fontan_props <- read.csv(file= "/users/1/blake561/Fontan Subcelltype Counts - MB250422 Revised.csv") #I need to edit this file to match the input format

ggplot(Fontan_props, aes(x = Condition, y = Proportion, fill = Cell.Type))+
  geom_col(position = "fill", width = 0.4)+
  scale_fill_manual(values = c(
    'Hepatocyte'    = '#33CC33',
    'Endothelial Cell'   = '#33CCCC',
    'HSC/mFB'       = '#CC99FF',
    'Macrophage'    = '#FF66FF',
    'Lymphocyte'    = 'salmon',
    'B Cell'        = '#FFCC00',
    'Cholangiocyte' = 'dodgerblue',
    'Plasma Cell'        = '#CC0000'
  ))

NASH_props <- read.csv(file= "/users/1/blake561/NASH Subcelltype Counts - MB250415 Revised.csv") #I need to edit this file to match the input format

ggplot(NASH_props, aes(x = Condition, y = Proportion, fill = Cell.Type))+
  geom_col(position = "fill", width = 0.4)+
  scale_fill_manual(values = c(
    'Hepatocyte'    = '#33CC33',
    'Endothelial Cell'   = '#33CCCC',
    'HSC/mFB'       = '#CC99FF',
    'Macrophage'    = '#FF66FF',
    'Lymphocyte'    = 'salmon',
    'B Cell'        = '#FFCC00',
    'Cholangiocyte' = 'dodgerblue',
    'Plasma Cell'        = '#CC0000'
  ))

##UMAPS##

#PAH

# Define your custom color palette
values <- c(
  'Hepatocyte'    = '#33CC33',
  'Endothelial Cell'   = '#33CCCC',
  'HSC/mFB'       = '#CC99FF',
  'Macrophage'    = '#FF66FF',
  'Lymphocyte'    = 'salmon',
  'B Cell'        = '#FFCC00',
  'Cholangiocyte' = 'dodgerblue',
  'Plasma Cell'        = '#CC0000'
)

# Generate the UMAP plot with your color mapping

#Overall

DimPlot(
  PAH,
  reduction = "umap",
  group.by = "cell_identity",
  label = TRUE,
  raster = FALSE,
  cols = values
)

#Split by Condition

DimPlot(
  PAH,
  reduction = "umap",
  group.by = "cell_identity",
  split.by = 'condition',
  label = TRUE,
  raster = FALSE,
  cols = values
)

#Fontan 

#Overall

DimPlot(
  Fontan,
  reduction = "umap",
  group.by = "cell_identity",
  label = TRUE,
  raster = FALSE,
  cols = values
)

#Split by Condition

DimPlot(
  Fontan,
  reduction = "umap",
  group.by = "cell_identity",
  split.by = 'orig.ident',
  label = TRUE,
  raster = FALSE,
  cols = values
)

#NASH 

#Fontan 

#Overall

DimPlot(
  NASH,
  reduction = "umap",
  group.by = "cell_identity",
  label = TRUE,
  raster = FALSE,
  cols = values
)

#Split by Condition

DimPlot(
  NASH,
  reduction = "umap",
  group.by = "cell_identity",
  split.by = 'condition',
  label = TRUE,
  raster = FALSE,
  cols = values
)


##Marker Gene Plots for NASH and FALD##

#Fontan

figure_order <- c(
  'MYOM1', 'CRP', 'CP', 'NNMT', 'ALDH1A2',                         # Hepatocyte 
  'PTPRB', 'ST6GALNAC3', 'STAB2', 'ROBO2', 'FLT1',                                 # Endothelial
  'CCBE1', 'ADAMTSL1', 'ADAMTS12', 'PDGFRA', 'COL25A1',                            # HSC
  'GRK3', 'TBXAS1', 'MCTP1', 'EPB41L3','RAB31',                                   # Macrophage
  'MCTP2', 'CD247', 'CD96', 'PYHIN1','PRKCQ',                                                      # Lymphocyte
  'CTNND2', 'ANXA4', 'PKHD1', 'DCDC2', 'RALYL'                                     # Cholangio
)

Idents(Fontan) <- Fontan$cell_identity

# Reorder and clean identities
Idents(Fontan) <- factor(Idents(Fontan), 
                         levels = c("Hepatocyte", "Endothelial Cell", "HSC/mFB", 
                                    "Macrophage", "Lymphocyte", "Cholangiocyte"))

# Subset and drop unused identities
Fontan_subset <- subset(Fontan, downsample = 1000)
Fontan_subset <- Fontan_subset[, !is.na(Idents(Fontan_subset))]
Idents(Fontan_subset) <- droplevels(Idents(Fontan_subset))

# Now make the heatmap
DoHeatmap(
  Fontan_subset, 
  features = figure_order, 
  group.colors = c(
    'Hepatocyte'        = '#33CC33',
    'Endothelial Cell'  = '#33CCCC',
    'HSC/mFB'           = '#CC99FF',
    'Macrophage'        = '#FF66FF',
    'Lymphocyte'        = 'salmon',
    'Cholangiocyte'     = 'dodgerblue'
  ), 
  group.bar.height = 0.05
)

#NASH 

Idents(NASH) <- NASH$seurat_clusters #Reassinging Numerical Labels for Renaming

NASH <- RenameIdents(NASH, "0" = "Hepatocyte")
NASH <- RenameIdents(NASH, "1" = "Endothelial Cell")
NASH <- RenameIdents(NASH, "2" = "Macrophage")
NASH <- RenameIdents(NASH, "3" = "Lymphocyte")
NASH <- RenameIdents(NASH, "4" = "Endothelial Cell")
NASH <- RenameIdents(NASH, "6" = "hep/mFB")
NASH <- RenameIdents(NASH, "7" = "Macrophage")
NASH <- RenameIdents(NASH, "8" = "Plasma Cell")
NASH <- RenameIdents(NASH, "9" = "Cholangiocyte")
NASH <- RenameIdents(NASH, "10" = "B Cell")
NASH <- RenameIdents(NASH, "11" = "Lymphocyte")

NASH$subcelltype <- Idents(NASH)

control_ids <- c("CON")  
NASH$condition <- ifelse(NASH$orig.ident %in% control_ids, "Control", "NASH")


figure_order <- c(
  'GAS2', 'AGMO', 'LURAP1L', 'ARHGEF26','SMOC1',                         # Hepatocyte (filtered)
  'PTPRB', 'ST6GALNAC3', 'STAB2', 'ROBO2', 'FLT1',                                 # Endothelial
  'LAMA2', 'CCBE1', 'ADGRB3', 'PDGFRA', 'COL25A1',                            # HSC
  'TLR2', 'TBXAS1', 'CD163', 'NDST3','RAB31',                                   # Macrophage
  'CARD11', 'CD247', 'SYTL3', 'PARP8','SKAP1', #Lymphocyte
  'BCL11A', 'BLK', 'FCRL5','FCRL1','POU2AF1', #B Cell
  'CTNND2', 'ANXA4', 'PKHD1', 'DCDC2', 'RALYL',   # Cholangio
  'DIAPH3', 'POLQ') #Plasma

Idents(NASH) <- NASH$cell_identity


# Reorder and clean identities
Idents(NASH) <- factor(Idents(NASH), 
                         levels = c("Hepatocyte", "Endothelial Cell", "HSC/mFB", 
                                    "Macrophage", "Lymphocyte", "B Cell", "Cholangiocyte", "Plasma Cell"))

# Subset and drop unused identities
NASH_subset <- subset(NASH, downsample = 1000)
NASH_subset <- NASH_subset[, !is.na(Idents(NASH_subset))]
Idents(NASH_subset) <- droplevels(Idents(NASH_subset))

# Now make the heatmap
DoHeatmap(
  NASH_subset, 
  features = figure_order, 
  group.colors = c(
    'Hepatocyte'    = '#33CC33',
    'Endothelial Cell'   = '#33CCCC',
    'HSC/mFB'       = '#CC99FF',
    'Macrophage'    = '#FF66FF',
    'Lymphocyte'    = 'salmon',
    'B Cell'        = '#FFCC00',
    'Cholangiocyte' = 'dodgerblue',
    'Plasma Cell'        = '#CC0000'
  ),
  group.bar.height = 0.05
)

# hep_markers <- FindConservedMarkers(NASH, ident.1 = "Hepatocyte", only.pos = TRUE, grouping.var = "orig.ident")
# hep_genes <- rownames(head(hep_markers, 20))
# FeaturePlot(NASH, features = hep_genes, min.cutoff = "q10")
# 
# endo_markers <- FindConservedMarkers(NASH, ident.1 = "Endothelial Cell", only.pos = TRUE, grouping.var = "orig.ident")
# endo_genes <- rownames(head(endo_markers, 20))
# FeaturePlot(NASH, features = endo_genes, min.cutoff = "q10")

# hsc_markers <- FindConservedMarkers(NASH, ident.1 = "HSC/mFB", only.pos = TRUE, grouping.var = "orig.ident")
# hsc_genes <- rownames(head(hsc_markers, 20))
# FeaturePlot(NASH, features = hsc_genes, min.cutoff = "q10")

# mac_markers <- FindConservedMarkers(NASH, ident.1 = "Macrophage", only.pos = TRUE, grouping.var = "orig.ident")
# mac_genes <- rownames(head(mac_markers, 20))
# FeaturePlot(NASH, features = mac_genes, min.cutoff = "q10")
# 
# lymph_markers <- FindConservedMarkers(NASH, ident.1 = "Lymphocyte", only.pos = TRUE, grouping.var = "orig.ident")
# lymph_genes <- rownames(head(lymph_markers, 20))
# FeaturePlot(NASH, features = lymph_genes, min.cutoff = "q10")

# bcell_markers <- FindConservedMarkers(NASH, ident.1 = "B Cell", only.pos = TRUE, grouping.var = "orig.ident")
# bcell_genes <- rownames(head(bcell_markers, 20))
# FeaturePlot(NASH, features = bcell_genes, min.cutoff = "q10")

plasma_markers <- FindConservedMarkers(NASH, ident.1 = "Plasma Cell", only.pos = TRUE, grouping.var = "orig.ident")
plasma_genes <- rownames(head(plasma_markers, 20))
FeaturePlot(NASH, features = plasma_genes, min.cutoff = "q10")

#Individual Cell Type Abundances for PAH, NASH, and FALD Datasets 

# ==============================================================================
# CELL TYPE RELATIVE ABUNDANCE PER SAMPLE - PAH, NASH, FALD LIVER DATASETS
# ==============================================================================

library(Seurat)
library(dplyr)
library(tidyr)

# ==============================================================================
# 1. PAH DATASET
# ==============================================================================

cat("=== PAH Dataset ===\n")

# Calculate abundance per individual sample
pah_abundance <- PAH@meta.data %>%
  group_by(sample, cell_identity) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(sample) %>%
  mutate(
    total_cells = sum(count),
    relative_abundance = (count / total_cells) * 100
  ) %>%
  ungroup()

# Add condition info
pah_condition_map <- PAH@meta.data %>%
  select(sample, condition) %>%
  distinct()

pah_abundance <- pah_abundance %>%
  left_join(pah_condition_map, by = "sample") %>%
  select(
    Sample = sample,
    Condition = condition,
    CellType = cell_identity,
    Count = count,
    TotalCells = total_cells,
    RelativeAbundance = relative_abundance
  )

# Create wide format (samples as rows, cell types as columns)
pah_wide <- pah_abundance %>%
  select(Sample, Condition, CellType, RelativeAbundance) %>%
  pivot_wider(
    names_from = CellType,
    values_from = RelativeAbundance,
    values_fill = 0
  )

# Export
write.csv(pah_abundance, "PAH_cell_type_abundance_long.csv", row.names = FALSE)
write.csv(pah_wide, "PAH_cell_type_abundance_wide.csv", row.names = FALSE)

cat("\nPAH abundance (long format):\n")
print(pah_abundance)
cat("\nPAH abundance (wide format):\n")
print(pah_wide)

# ==============================================================================
# 2. NASH DATASET
# ==============================================================================
# Cell types stored in: subcelltype
# Samples stored in: orig.ident (CON, NASH1, NASH2, NASH3)

cat("\n\n=== NASH Dataset ===\n")

# Set identities to subcelltype
Idents(NASH) <- "subcelltype"

# Calculate abundance per sample
nash_abundance <- NASH@meta.data %>%
  group_by(orig.ident, subcelltype) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(orig.ident) %>%
  mutate(
    total_cells = sum(count),
    relative_abundance = (count / total_cells) * 100
  ) %>%
  ungroup() %>%
  mutate(
    # Add condition based on sample name
    Condition = ifelse(orig.ident == "CON", "Control", "NASH")
  ) %>%
  select(
    Sample = orig.ident,
    Condition,
    CellType = subcelltype,
    Count = count,
    TotalCells = total_cells,
    RelativeAbundance = relative_abundance
  )

# Create wide format
nash_wide <- nash_abundance %>%
  select(Sample, Condition, CellType, RelativeAbundance) %>%
  pivot_wider(
    names_from = CellType,
    values_from = RelativeAbundance,
    values_fill = 0
  )

# Export
write.csv(nash_abundance, "NASH_cell_type_abundance_long.csv", row.names = FALSE)
write.csv(nash_wide, "NASH_cell_type_abundance_wide.csv", row.names = FALSE)

cat("\nNASH abundance (long format):\n")
print(nash_abundance)
cat("\nNASH abundance (wide format):\n")
print(nash_wide)

# ==============================================================================
# 3. FALD/FONTAN DATASET
# ==============================================================================
# Cell types stored in: cell_identity
# Samples stored in: sample
# Condition stored in: orig.ident (Fontan vs Control)

cat("\n\n=== FALD/Fontan Dataset ===\n")

# Calculate abundance per individual sample
fald_abundance <- Fontan@meta.data %>%
  group_by(sample, cell_identity) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(sample) %>%
  mutate(
    total_cells = sum(count),
    relative_abundance = (count / total_cells) * 100
  ) %>%
  ungroup()

# Add condition info from orig.ident
fald_condition_map <- Fontan@meta.data %>%
  select(sample, orig.ident) %>%
  distinct() %>%
  mutate(Condition = ifelse(orig.ident == "Fontan", "FALD", "Control"))

fald_abundance <- fald_abundance %>%
  left_join(fald_condition_map, by = "sample") %>%
  select(
    Sample = sample,
    Condition,
    CellType = cell_identity,
    Count = count,
    TotalCells = total_cells,
    RelativeAbundance = relative_abundance
  )

# Create wide format
fald_wide <- fald_abundance %>%
  select(Sample, Condition, CellType, RelativeAbundance) %>%
  pivot_wider(
    names_from = CellType,
    values_from = RelativeAbundance,
    values_fill = 0
  )

# Export
write.csv(fald_abundance, "FALD_cell_type_abundance_long.csv", row.names = FALSE)
write.csv(fald_wide, "FALD_cell_type_abundance_wide.csv", row.names = FALSE)

cat("\nFALD abundance (long format):\n")
print(fald_abundance)
cat("\nFALD abundance (wide format):\n")
print(fald_wide)

# ==============================================================================
# 4. COMBINED SUMMARY (OPTIONAL)
# ==============================================================================

cat("\n\n=== Creating Combined Dataset ===\n")

# Add dataset identifier and combine
pah_abundance$Dataset <- "PAH"
nash_abundance$Dataset <- "NASH"
fald_abundance$Dataset <- "FALD"

combined_abundance <- bind_rows(pah_abundance, nash_abundance, fald_abundance) %>%
  select(Dataset, Sample, Condition, CellType, Count, TotalCells, RelativeAbundance)

write.csv(combined_abundance, "Combined_liver_cell_type_abundance.csv", row.names = FALSE)

cat("Combined dataset exported!\n")
cat("\nSummary of samples per dataset:\n")
combined_abundance %>%
  select(Dataset, Sample, Condition) %>%
  distinct() %>%
  print()

cat("\n=== All files exported successfully! ===\n")
cat("Files created:\n")
cat("- PAH_cell_type_abundance_long.csv\n")
cat("- PAH_cell_type_abundance_wide.csv\n")
cat("- NASH_cell_type_abundance_long.csv\n")
cat("- NASH_cell_type_abundance_wide.csv\n")
cat("- FALD_cell_type_abundance_long.csv\n")
cat("- FALD_cell_type_abundance_wide.csv\n")
cat("- Combined_liver_cell_type_abundance.csv\n")

colnames(NASH_orig@meta.data)
head(NASH_orig@meta.data, 20)


# Check all metadata columns
colnames(NASH@meta.data)

# Check total cells in the object
ncol(NASH)

# Check what's in orig.ident
table(NASH$orig.ident)

# Look at the first few rows of metadata to see structure
head(NASH@meta.data, 20)

# Check subcelltype counts
table(NASH$subcelltype)

# If there's a 'sample' column, check it
if("sample" %in% colnames(NASH@meta.data)) {
  table(NASH$sample)
}

# Check for any column that might contain sample info
# Look for columns with a small number of unique values (likely sample IDs)
sapply(NASH@meta.data, function(x) length(unique(x)))

# Find a barcode that appears in both NASH1 and NASH2
test_barcode <- intersect(barcodes_nash1, barcodes_nash2)[1]
cat("Testing barcode:", test_barcode, "\n")

# Get full cell names
cell_nash1 <- paste0(test_barcode, "-1_2")
cell_nash2 <- paste0(test_barcode, "-1_3")

# Extract expression for top variable genes
var_genes <- head(VariableFeatures(NASH), 20)

# Compare expression
expr_nash1 <- as.numeric(NASH@assays$RNA@counts[var_genes, cell_nash1])
expr_nash2 <- as.numeric(NASH@assays$RNA@counts[var_genes, cell_nash2])

# Correlation - if same cell, should be r = 1.0
cor(expr_nash1, expr_nash2)

# Direct comparison
data.frame(
  gene = var_genes,
  NASH1 = expr_nash1,
  NASH2 = expr_nash2,
  identical = expr_nash1 == expr_nash2
)

# First, check the structure of the original object
cat("=== NASH_orig Structure ===\n")
colnames(NASH_orig@meta.data)
table(NASH_orig$orig.ident)

# Check total cells
ncol(NASH_orig)

# Look at barcode structure
head(colnames(NASH_orig), 10)

# Check barcode suffixes
table(substr(colnames(NASH_orig), nchar(colnames(NASH_orig))-1, nchar(colnames(NASH_orig))))

# Extract barcodes for each group (adjust based on what orig.ident shows)
barcodes_orig_nash1 <- gsub("-1_2$", "", colnames(NASH_orig)[NASH_orig$orig.ident == "NASH1"])
barcodes_orig_nash2 <- gsub("-1_3$", "", colnames(NASH_orig)[NASH_orig$orig.ident == "NASH2"])
barcodes_orig_nash3 <- gsub("-1_4$", "", colnames(NASH_orig)[NASH_orig$orig.ident == "NASH3"])

# Check overlap
cat("\n=== Barcode Overlap ===\n")
cat("NASH1 cells:", length(barcodes_orig_nash1), "\n")
cat("NASH2 cells:", length(barcodes_orig_nash2), "\n")
cat("NASH3 cells:", length(barcodes_orig_nash3), "\n")
cat("Overlap NASH1-NASH2:", length(intersect(barcodes_orig_nash1, barcodes_orig_nash2)), "\n")
cat("Overlap %:", length(intersect(barcodes_orig_nash1, barcodes_orig_nash2)) / length(barcodes_orig_nash1) * 100, "%\n")

# Expression correlation test
cat("\n=== Expression Correlation Test ===\n")
test_barcode_orig <- intersect(barcodes_orig_nash1, barcodes_orig_nash2)[1]
cat("Testing barcode:", test_barcode_orig, "\n")

cell_orig_nash1 <- paste0(test_barcode_orig, "-1_2")
cell_orig_nash2 <- paste0(test_barcode_orig, "-1_3")

# Get variable features (or use first 20 genes if none set)
if(length(VariableFeatures(NASH_orig)) > 0) {
  test_genes <- head(VariableFeatures(NASH_orig), 20)
} else {
  test_genes <- head(rownames(NASH_orig), 20)
}

expr_orig_nash1 <- as.numeric(NASH_orig@assays$RNA@counts[test_genes, cell_orig_nash1])
expr_orig_nash2 <- as.numeric(NASH_orig@assays$RNA@counts[test_genes, cell_orig_nash2])

cat("Correlation:", cor(expr_orig_nash1, expr_orig_nash2), "\n")

# Direct comparison
data.frame(
  gene = test_genes,
  NASH1 = expr_orig_nash1,
  NASH2 = expr_orig_nash2,
  identical = expr_orig_nash1 == expr_orig_nash2
)

table(PAH$cell_identity)

library(tidyverse)

# Calculate HSC/mFB abundance by sample
hsc_abundance <- PAH@meta.data %>%
  group_by(sample) %>%
  summarise(
    total_cells = n(),
    HSC_mFB_count = sum(cell_identity == "HSC/mFB"),
    HSC_mFB_percent = HSC_mFB_count / total_cells * 100
  ) %>%
  # Add condition info
  
  mutate(condition = ifelse(grepl("^C|Control|CON", sample, ignore.case = TRUE), "Control", "PAH"))

# View results
print(hsc_abundance)

# Summary by condition
hsc_abundance %>%
  group_by(condition) %>%
  summarise(
    n_samples = n(),
    mean_HSC_percent = mean(HSC_mFB_percent),
    sd_HSC_percent = sd(HSC_mFB_percent)
  )

# Export
write.csv(hsc_abundance, "PAH_HSC_mFB_abundance_by_sample.csv", row.names = FALSE)
