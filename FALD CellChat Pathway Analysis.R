# ==============================================================================
# FONTAN INDIVIDUAL SAMPLE CELLCHAT ANALYSIS - COMPLETE WORKFLOW
# ==============================================================================

# Load required libraries
library(CellChat)
library(Seurat)
library(ggplot2)
library(dplyr)

# Load ggrepel for enhanced plots
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
library(ggrepel)

# ==============================================================================
# STEP 1: SUBSET INDIVIDUAL SAMPLES
# ==============================================================================

# Subset each sample using the "sample" column
CF1 <- subset(Fontan, subset = sample == "CF1")
CF2 <- subset(Fontan, subset = sample == "CF2")
F1 <- subset(Fontan, subset = sample == "F1") 
F2 <- subset(Fontan, subset = sample == "F2")
F3 <- subset(Fontan, subset = sample == "F3")
F4 <- subset(Fontan, subset = sample == "F4")

# Verify the subsets
cat("=== SAMPLE VERIFICATION ===\n")
cat("CF1:", ncol(CF1), "cells\n")
cat("CF2:", ncol(CF2), "cells\n")
cat("F1:", ncol(F1), "cells\n") 
cat("F2:", ncol(F2), "cells\n")
cat("F3:", ncol(F3), "cells\n")
cat("F4:", ncol(F4), "cells\n")

# Check cell types in each sample
cat("\nCell types in CF1:\n")
print(table(CF1$cell_identity))

cat("\nCell types in F1:\n")
print(table(F1$cell_identity))

# ==============================================================================
# STEP 2: CONVERT INDIVIDUAL SEURAT OBJECTS TO CELLCHAT
# ==============================================================================

# Function to create CellChat from Seurat object
create_cellchat_from_seurat <- function(seurat_obj, sample_name) {
  
  cat("Processing sample:", sample_name, "\n")
  cat("Number of cells:", ncol(seurat_obj), "\n")
  
  # Extract expression data (use normalized data) - updated for newer Seurat
  expression_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  
  # Extract metadata and add cell type column
  cell_metadata <- seurat_obj@meta.data
  
  # The cell types are already in cell_identity column
  cat("Available metadata columns:", colnames(cell_metadata), "\n")
  
  # Create CellChat object using the cell_identity column
  cellchat <- createCellChat(object = expression_data, meta = cell_metadata, group.by = "cell_identity")
  
  # Add a samples column to avoid warnings
  cellchat@meta$samples <- sample_name
  
  # Set cell type identity
  cellchat <- setIdent(cellchat, ident.use = "cell_identity")
  
  # Set database (human data)
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB
  
  # Process the CellChat object
  cat("Processing CellChat pipeline for", sample_name, "...\n")
  
  # Subset data to use only expressed genes
  cellchat <- subsetData(cellchat)
  
  # Identify overexpressed genes and interactions
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Project gene expression onto PPI network (optional step, skip if function doesn't exist)
  tryCatch({
    cellchat <- projectData(cellchat, PPI.human)
    cat("Successfully projected data onto PPI network\n")
  }, error = function(e) {
    cat("Note: projectData function not available, skipping PPI projection step\n")
    cat("This is fine - CellChat will still work without this step\n")
  })
  
  # Compute communication probabilities
  cellchat <- computeCommunProb(cellchat)
  
  # Filter out communications with few cells
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # Compute pathway-level communication
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Aggregate cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  
  cat("Completed processing for", sample_name, "\n")
  cat("Number of pathways:", length(cellchat@netP$pathways), "\n\n")
  
  return(cellchat)
}

# ==============================================================================
# STEP 3: PROCESS ALL SAMPLES
# ==============================================================================

# Create list of Seurat objects
seurat_objects <- list(
  "CF1" = CF1,
  "CF2" = CF2,
  "F1" = F1,
  "F2" = F2, 
  "F3" = F3,
  "F4" = F4
)

# Process each sample through CellChat
cellchat_objects <- list()

for(sample_name in names(seurat_objects)) {
  cellchat_objects[[sample_name]] <- create_cellchat_from_seurat(
    seurat_objects[[sample_name]], 
    sample_name
  )
}

# ==============================================================================
# STEP 4: ORGANIZE INTO CONTROL AND FONTAN GROUPS
# ==============================================================================

# Separate into control and Fontan groups
control_cellchats <- list(
  "CF1" = cellchat_objects[["CF1"]],
  "CF2" = cellchat_objects[["CF2"]]
)

fontan_cellchats <- list(
  "F1" = cellchat_objects[["F1"]],
  "F2" = cellchat_objects[["F2"]],
  "F3" = cellchat_objects[["F3"]],
  "F4" = cellchat_objects[["F4"]]
)

# Verify the groups
cat("=== GROUP VERIFICATION ===\n")
cat("Control samples:", names(control_cellchats), "\n")
cat("Fontan samples:", names(fontan_cellchats), "\n")

# Check cell types in each sample (should be consistent)
cat("\nCell types in each sample:\n")
for(sample_name in names(cellchat_objects)) {
  cell_types <- levels(cellchat_objects[[sample_name]]@idents)
  cat(sample_name, ":", paste(cell_types, collapse = ", "), "\n")
}

# ==============================================================================
# STEP 5: STATISTICAL ANALYSIS FUNCTIONS
# ==============================================================================

# Function to calculate pathway strengths for a specific cell type
calculate_celltype_pathways <- function(cellchat_obj, cell_types) {
  if (is.null(cellchat_obj@netP$prob)) return(NULL)
  
  prob_array <- cellchat_obj@netP$prob
  pathway_names <- dimnames(prob_array)[[3]]
  all_celltypes <- rownames(prob_array[,,1])
  
  # Filter to cell types that actually exist
  existing_cells <- intersect(cell_types, all_celltypes)
  if (length(existing_cells) == 0) return(NULL)
  
  pathway_strengths <- setNames(rep(0, length(pathway_names)), pathway_names)
  
  for (pathway in pathway_names) {
    pathway_matrix <- prob_array[,,pathway]
    # Sum all communications involving these cells (both sending and receiving)
    outgoing <- sum(pathway_matrix[existing_cells, ], na.rm = TRUE)
    incoming <- sum(pathway_matrix[, existing_cells], na.rm = TRUE)
    pathway_strengths[pathway] <- outgoing + incoming
  }
  
  return(pathway_strengths)
}

# Main statistical analysis function (for n=2 vs n=4)
analyze_pathways_with_stats <- function(control_list, fontan_list, cell_types, cell_name) {
  
  cat("\n=== Statistical Analysis for", cell_name, "===\n")
  
  # Calculate pathway strengths for each sample
  control_strengths <- lapply(control_list, function(obj) {
    calculate_celltype_pathways(obj, cell_types)
  })
  
  fontan_strengths <- lapply(fontan_list, function(obj) {
    calculate_celltype_pathways(obj, cell_types)
  })
  
  # Remove NULL results
  control_strengths <- control_strengths[!sapply(control_strengths, is.null)]
  fontan_strengths <- fontan_strengths[!sapply(fontan_strengths, is.null)]
  
  if (length(control_strengths) == 0 || length(fontan_strengths) == 0) {
    cat("No valid data for", cell_name, "\n")
    return(NULL)
  }
  
  # Get common pathways across all samples
  all_pathway_names <- Reduce(intersect, lapply(c(control_strengths, fontan_strengths), names))
  
  # Create matrices for statistical testing
  control_matrix <- do.call(cbind, lapply(control_strengths, function(x) x[all_pathway_names]))
  fontan_matrix <- do.call(cbind, lapply(fontan_strengths, function(x) x[all_pathway_names]))
  
  # Perform statistical tests for each pathway
  results_list <- list()
  
  for (pathway in all_pathway_names) {
    control_vals <- control_matrix[pathway, ]
    fontan_vals <- fontan_matrix[pathway, ]
    
    # Calculate basic statistics
    control_mean <- mean(control_vals, na.rm = TRUE)
    fontan_mean <- mean(fontan_vals, na.rm = TRUE)
    control_sd <- sd(control_vals, na.rm = TRUE)
    fontan_sd <- sd(fontan_vals, na.rm = TRUE)
    
    # Perform statistical test
    test_result <- tryCatch({
      wilcox.test(control_vals, fontan_vals, exact = FALSE)
    }, error = function(e) {
      t.test(control_vals, fontan_vals)
    })
    
    # Calculate effect size (Cohen's d)
    pooled_sd <- sqrt(((length(control_vals)-1) * control_sd^2 + (length(fontan_vals)-1) * fontan_sd^2) / 
                        (length(control_vals) + length(fontan_vals) - 2))
    cohens_d <- (fontan_mean - control_mean) / pooled_sd
    
    results_list[[pathway]] <- data.frame(
      Pathway = pathway,
      Control_Mean = round(control_mean, 4),
      Control_SD = round(control_sd, 4),
      Fontan_Mean = round(fontan_mean, 4),
      Fontan_SD = round(fontan_sd, 4),
      Difference = round(fontan_mean - control_mean, 4),
      P_Value = test_result$p.value,
      Cohens_D = round(cohens_d, 3),
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results and apply multiple testing correction
  results_df <- do.call(rbind, results_list)
  results_df$FDR <- p.adjust(results_df$P_Value, method = "fdr")
  
  # Add significance categories
  results_df$Significance <- with(results_df, ifelse(
    FDR < 0.001, "***",
    ifelse(FDR < 0.01, "**",
           ifelse(FDR < 0.05, "*",
                  ifelse(FDR < 0.1, ".", "")))))
  
  # Sort by FDR
  results_df <- results_df[order(results_df$FDR), ]
  
  # Print results
  cat("\nTop 10 pathways for", cell_name, ":\n")
  print(head(results_df[, c("Pathway", "Control_Mean", "Fontan_Mean", "Difference", 
                            "P_Value", "FDR", "Cohens_D", "Significance")], 10))
  
  # ENHANCED ANALYSIS: Create publication-quality plots
  cat("\n--- Enhanced Analysis for", cell_name, "---\n")
  enhanced_plots <- create_enhanced_pathway_plots(results_df, cell_name)
  
  return(results_df)
}

# ==============================================================================
# STEP 6: ENHANCED VISUALIZATION FUNCTIONS
# ==============================================================================

create_enhanced_pathway_plots <- function(results_df, cell_name) {
  
  # Filter out pathways where both conditions are near zero
  filtered_results <- results_df[results_df$Control_Mean > 0.001 | results_df$Fontan_Mean > 0.001, ]
  
  # Calculate Log2 fold change
  filtered_results$Log2_Fold_Change <- log2((filtered_results$Fontan_Mean + 0.001) / 
                                              (filtered_results$Control_Mean + 0.001))
  
  # 1. TOP 5 PATHWAYS BAR PLOT
  top_5 <- head(filtered_results[order(-abs(filtered_results$Log2_Fold_Change)), ], 5)
  
  barplot_top5 <- ggplot(top_5, aes(x = reorder(Pathway, Log2_Fold_Change), y = Log2_Fold_Change)) +
    geom_col(aes(fill = Log2_Fold_Change > 0), width = 0.7) +
    geom_text(aes(label = paste0(round(Log2_Fold_Change, 2), "x")), 
              hjust = ifelse(top_5$Log2_Fold_Change > 0, -0.1, 1.1), 
              size = 4) +
    scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                      labels = c("TRUE" = "Increased in Fontan", "FALSE" = "Decreased in Fontan")) +
    coord_flip() +
    labs(
      title = paste("Top 5 Most Altered", cell_name, "Pathways in Fontan"),
      subtitle = "Ranked by Log2 Fold Change Magnitude",
      x = "Pathway",
      y = "Log2 Fold Change (Fontan/Control)",
      fill = "Direction"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  # 2. TOP 10 PATHWAYS BAR PLOT
  top_10 <- head(filtered_results[order(-abs(filtered_results$Log2_Fold_Change)), ], 10)
  
  barplot_top10 <- ggplot(top_10, aes(x = reorder(Pathway, Log2_Fold_Change), y = Log2_Fold_Change)) +
    geom_col(aes(fill = Log2_Fold_Change > 0), width = 0.7) +
    geom_text(aes(label = paste0(round(Log2_Fold_Change, 2), "x")), 
              hjust = ifelse(top_10$Log2_Fold_Change > 0, -0.1, 1.1), 
              size = 3) +
    scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                      labels = c("TRUE" = "Increased in Fontan", "FALSE" = "Decreased in Fontan")) +
    coord_flip() +
    labs(
      title = paste("Top 10 Most Altered", cell_name, "Pathways in Fontan"),
      subtitle = "Ranked by Log2 Fold Change Magnitude",
      x = "Pathway",
      y = "Log2 Fold Change (Fontan/Control)",
      fill = "Direction"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  # Print plots
  cat("\n--- Top 5 Pathways ---\n")
  print(barplot_top5)
  
  cat("\n--- Top 10 Pathways ---\n")
  print(barplot_top10)
  
  # Create summary tables
  summary_table_5 <- top_5[, c("Pathway", "Control_Mean", "Fontan_Mean", "Fontan_SD", 
                               "Log2_Fold_Change", "FDR", "Significance")]
  
  summary_table_10 <- top_10[, c("Pathway", "Control_Mean", "Fontan_Mean", "Fontan_SD", 
                                 "Log2_Fold_Change", "FDR", "Significance")]
  
  cat("\n=== TOP 5 ALTERED PATHWAYS SUMMARY ===\n")
  cat("Cell Type:", cell_name, "\n\n")
  print(summary_table_5)
  
  cat("\n=== TOP 10 ALTERED PATHWAYS SUMMARY ===\n")
  cat("Cell Type:", cell_name, "\n\n")
  print(summary_table_10)
  
  return(list(
    barplot_top5 = barplot_top5,
    barplot_top10 = barplot_top10,
    summary_table_5 = summary_table_5,
    summary_table_10 = summary_table_10
  ))
}

# Function for cross-cell-type comparison (showing fold changes regardless of significance)
create_cross_celltype_comparison <- function(all_results) {
  
  # Combine top pathways from all cell types
  combined_data <- data.frame()
  
  for(cell_type in names(all_results)) {
    if(!is.null(all_results[[cell_type]])) {
      df <- all_results[[cell_type]]
      df$Log2_Fold_Change <- log2((df$Fontan_Mean + 0.001) / (df$Control_Mean + 0.001))
      df$Cell_Type <- cell_type
      
      # Get top 10 pathways by fold change magnitude (regardless of significance)
      df_top <- head(df[order(-abs(df$Log2_Fold_Change)), ], 10)
      
      combined_data <- rbind(combined_data, df_top)
    }
  }
  
  if(nrow(combined_data) == 0) {
    cat("No pathway data found across cell types\n")
    return(NULL)
  }
  
  # Create heatmap
  heatmap_plot <- ggplot(combined_data, aes(x = Cell_Type, y = reorder(Pathway, Log2_Fold_Change))) +
    geom_tile(aes(fill = Log2_Fold_Change), color = "white", linewidth = 0.5) +
    geom_text(aes(label = paste0(round(Log2_Fold_Change, 1), "x")), 
              color = "black", fontface = "bold", size = 2.5) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", 
      midpoint = 0, name = "Log2 FC\n(Fontan/Control)"
    ) +
    labs(
      title = "Pathway Alterations Across Cell Types in Fontan",
      subtitle = "Color = Log2 Fold Change, Text = Fold Change Value",
      x = "Cell Type", 
      y = "Pathway"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9)
    )
  
  print(heatmap_plot)
  return(heatmap_plot)
}

# ==============================================================================
# STEP 7: RUN ANALYSIS FOR ALL CELL TYPES
# ==============================================================================

# Define cell types to analyze - FOCUS ON 4 MAIN LIVER CELL TYPES
# Check what cell types you actually have:
cat("Available cell types:", levels(factor(CF1$cell_identity)), "\n")

# Define the 4 cell types you want to analyze
cell_types_to_analyze <- list(
  "Hepatocytes" = "Hepatocyte",
  "Endothelial_Cells" = "Endothelial Cell",
  "HSCs" = "HSC/mFB",
  "Macrophages" = "Macrophage"
)

# Store results
all_results <- list()

# Run analysis for each cell type
for (cell_name in names(cell_types_to_analyze)) {
  all_results[[cell_name]] <- analyze_pathways_with_stats(
    control_list = control_cellchats,
    fontan_list = fontan_cellchats,
    cell_types = cell_types_to_analyze[[cell_name]],
    cell_name = cell_name
  )
}

# ==============================================================================
# STEP 8: SAVE RESULTS
# ==============================================================================

# Save the CellChat objects and all results for future use
save(cellchat_objects, control_cellchats, fontan_cellchats, all_results,
     file = "fontan_individual_cellchat_analysis.RData")
