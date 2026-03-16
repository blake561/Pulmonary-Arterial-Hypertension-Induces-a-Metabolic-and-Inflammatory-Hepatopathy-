# ==============================================================================
# NASH REVISED DATASET - VISUALIZATION CODE
# ==============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# ==============================================================================
# STEP 1: SET UP IDENTITIES AND CONDITION
# ==============================================================================

# Set cell type identity
Idents(NASH) <- NASH$subcelltype

# Verify condition column exists
table(NASH$condition)
table(NASH$subcelltype)

# ==============================================================================
# STEP 2: OVERALL UMAP WITH CELL TYPE COLORS
# ==============================================================================

# Define custom color palette for cell types
celltype_colors <- c(
  'Hepatocyte'       = '#33CC33',
  'Endothelial Cell' = '#33CCCC',
  'HSC/mFB'          = '#CC99FF',
  'Macrophage'       = '#FF66FF',
  'Lymphocyte'       = 'salmon',
  'Cholangiocyte'    = 'dodgerblue'
)

# Overall UMAP
DimPlot(
  NASH,
  reduction = "umap",
  group.by = "subcelltype",
  label = FALSE,
  raster = FALSE,
  cols = celltype_colors
) + ggtitle("NASH Dataset - All Cell Types")

# ==============================================================================
# STEP 3: UMAP SPLIT BY CONDITION
# ==============================================================================

DimPlot(
  NASH,
  reduction = "umap",
  group.by = "subcelltype",
  split.by = "condition",
  label = FALSE,
  raster = FALSE,
  cols = celltype_colors
) + ggtitle("NASH Dataset - Split by Condition")

# ==============================================================================
# STEP 4: CELL-TYPE SPECIFIC UMAPS
# ==============================================================================

# Subset each cell type
NASH_Hepatocyte <- subset(NASH, subset = subcelltype == "Hepatocyte")
NASH_Endothelial <- subset(NASH, subset = subcelltype == "Endothelial Cell")
NASH_HSC <- subset(NASH, subset = subcelltype == "HSC/mFB")
NASH_Macrophage <- subset(NASH, subset = subcelltype == "Macrophage")

# Define condition colors
condition_colors <- c("Control" = "darkgoldenrod2", "NASH" = "orchid4")

# ── NASH: Hepatocyte ───────────────────────────────────────────────
NASH_Hepatocyte$ColorByCondition <- NASH_Hepatocyte$condition

DimPlot(
  NASH_Hepatocyte,
  reduction = "umap",
  split.by = "condition",
  group.by = "ColorByCondition",
  cols = condition_colors,
  raster = FALSE
) + ggtitle("Hepatocytes")

# ── NASH: Endothelial ───────────────────────────────────────────────
NASH_Endothelial$ColorByCondition <- NASH_Endothelial$condition

DimPlot(
  NASH_Endothelial,
  reduction = "umap",
  split.by = "condition",
  group.by = "ColorByCondition",
  cols = condition_colors,
  raster = FALSE
) + ggtitle("Endothelial Cells")

# ── NASH: HSC/mFB ───────────────────────────────────────────────────
NASH_HSC$ColorByCondition <- NASH_HSC$condition

DimPlot(
  NASH_HSC,
  reduction = "umap",
  split.by = "condition",
  group.by = "ColorByCondition",
  cols = condition_colors,
  raster = FALSE
) + ggtitle("HSC/mFB")

# ── NASH: Macrophage ────────────────────────────────────────────────
NASH_Macrophage$ColorByCondition <- NASH_Macrophage$condition

DimPlot(
  NASH_Macrophage,
  reduction = "umap",
  split.by = "condition",
  group.by = "ColorByCondition",
  cols = condition_colors,
  raster = FALSE
) + ggtitle("Macrophages")

# ==============================================================================
# STEP 5: RELATIVE ABUNDANCE BAR PLOT
# ==============================================================================

# Calculate cell counts by condition and cell type
cell_counts <- NASH@meta.data %>%
  group_by(condition, subcelltype) %>%
  summarise(count = n(), .groups = "drop")

# Calculate proportions within each condition
cell_props <- cell_counts %>%
  group_by(condition) %>%
  mutate(
    total = sum(count),
    proportion = count / total
  ) %>%
  ungroup()

# View the proportions
print(cell_props)

# Export to CSV
write.csv(cell_props, "NASH_Revised_Subcelltype_Proportions.csv", row.names = FALSE)

# Create stacked bar plot
ggplot(cell_props, aes(x = condition, y = proportion, fill = subcelltype)) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = c(
    'Hepatocyte'       = '#33CC33',
    'Endothelial Cell' = '#33CCCC',
    'HSC/mFB'          = '#CC99FF',
    'Macrophage'       = '#FF66FF',
    'Lymphocyte'       = 'salmon',
    'Cholangiocyte'    = 'dodgerblue'
  )) +
  labs(
    x = "Condition",
    y = "Proportion",
    fill = "Cell Type",
    title = "Cell Type Proportions in NASH Dataset"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12)
  )

# Also create per-sample proportions for more detailed view
cell_counts_sample <- NASH@meta.data %>%
  group_by(sample, condition, subcelltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(
    total = sum(count),
    proportion = count / total
  ) %>%
  ungroup()

# Export sample-level proportions
write.csv(cell_counts_sample, "NASH_Revised_Subcelltype_Proportions_BySample.csv", row.names = FALSE)

# Sample-level bar plot
ggplot(cell_counts_sample, aes(x = sample, y = proportion, fill = subcelltype)) +
  geom_col(position = "fill", width = 0.7) +
  scale_fill_manual(values = c(
    'Hepatocyte'       = '#33CC33',
    'Endothelial Cell' = '#33CCCC',
    'HSC/mFB'          = '#CC99FF',
    'Macrophage'       = '#FF66FF',
    'Lymphocyte'       = 'salmon',
    'Cholangiocyte'    = 'dodgerblue'
  )) +
  labs(
    x = "Sample",
    y = "Proportion",
    fill = "Cell Type",
    title = "Cell Type Proportions by Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ==============================================================================
# STEP 6: MARKER GENE HEATMAP
# ==============================================================================

# Define marker genes for each cell type
figure_order <- c(
  # Hepatocyte
  'GAS2', 'AGMO', 'CPS1',
  'CYP3A5', 'ACSM2A',
  # Endothelial
  'PTPRB', 'ST6GALNAC3', 'STAB2', 'ROBO2', 'FLT1',
  # HSC/mFB
  'LAMA2', 'CCBE1', 'ADGRB3', 'PDGFRA', 'COL25A1',
  # Macrophage
  'TLR2', 'TBXAS1', 'CD163', 'NDST3', 'RAB31',
  # Lymphocyte
  'CARD11', 'CD247', 'SYTL3', 'PARP8', 'SKAP1',
  # Cholangiocyte
  'CTNND2', 'ANXA4', 'PKHD1', 'DCDC2', 'RALYL'
)

# Set identity to subcelltype
Idents(NASH) <- NASH$subcelltype

# Reorder identities for heatmap
Idents(NASH) <- factor(Idents(NASH), 
                       levels = c("Hepatocyte", "Endothelial Cell", "HSC/mFB", 
                                  "Macrophage", "Lymphocyte", "Cholangiocyte"))

# Downsample for visualization
set.seed(42)
NASH_subset <- subset(NASH, downsample = 1000)
NASH_subset <- NASH_subset[, !is.na(Idents(NASH_subset))]
Idents(NASH_subset) <- droplevels(Idents(NASH_subset))

# Create heatmap
DoHeatmap(
  NASH_subset, 
  features = figure_order, 
  group.colors = c(
    'Hepatocyte'       = '#33CC33',
    'Endothelial Cell' = '#33CCCC',
    'HSC/mFB'          = '#CC99FF',
    'Macrophage'       = '#FF66FF',
    'Lymphocyte'       = 'salmon',
    'Cholangiocyte'    = 'dodgerblue'
  ), 
  group.bar.height = 0.05
) + ggtitle("NASH Dataset - Marker Gene Expression")

# ==============================================================================
# STEP 7: SAVE PLOTS (OPTIONAL)
# ==============================================================================

# Overall UMAP
pdf("NASH_Revised_UMAP_Overall.pdf", width = 10, height = 8)
DimPlot(NASH, reduction = "umap", group.by = "subcelltype", label = TRUE, 
        raster = FALSE, cols = celltype_colors) + ggtitle("NASH Dataset - All Cell Types")
dev.off()

# Split UMAP
pdf("NASH_Revised_UMAP_SplitByCondition.pdf", width = 14, height = 6)
DimPlot(NASH, reduction = "umap", group.by = "subcelltype", split.by = "condition",
        label = TRUE, raster = FALSE, cols = celltype_colors) + ggtitle("NASH Dataset - Split by Condition")
dev.off()

# Proportion bar plot
pdf("NASH_Revised_CellType_Proportions.pdf", width = 8, height = 6)
ggplot(cell_props, aes(x = condition, y = proportion, fill = subcelltype)) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = celltype_colors) +
  labs(x = "Condition", y = "Proportion", fill = "Cell Type", 
       title = "Cell Type Proportions in NASH Dataset") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# Heatmap
pdf("NASH_Revised_MarkerGene_Heatmap.pdf", width = 12, height = 10)
DoHeatmap(NASH_subset, features = figure_order, 
          group.colors = celltype_colors, group.bar.height = 0.05) + 
  ggtitle("NASH Dataset - Marker Gene Expression")
dev.off()

cat("\n=== VISUALIZATION COMPLETE ===\n")
cat("Generated:\n")
cat("1. Overall UMAP with cell type colors\n")
cat("2. UMAP split by condition\n")
cat("3. Cell-type specific UMAPs (Control vs NASH)\n")
cat("4. Cell type proportion bar plots (by condition and by sample)\n")
cat("5. Marker gene heatmap\n")
cat("6. CSV files with proportion data\n")