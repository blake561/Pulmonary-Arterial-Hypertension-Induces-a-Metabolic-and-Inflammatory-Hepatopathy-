# Load required libraries
library(CellChat)
library(Seurat)
library(ggplot2)
library(dplyr)

# ==============================================================================
# STEP 1: CONVERT INDIVIDUAL SEURAT OBJECTS TO CELLCHAT
# ==============================================================================

# Function to create CellChat from Seurat object
create_cellchat_from_seurat <- function(seurat_obj, sample_name) {
  
  cat("Processing sample:", sample_name, "\n")
  cat("Number of cells:", ncol(seurat_obj), "\n")
  
  # Extract expression data (use normalized data) - fix for newer Seurat versions
  expression_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  
  # Extract metadata
  cell_metadata <- seurat_obj@meta.data
  
  # Check what cell type column you have (adjust as needed)
  cat("Available metadata columns:", colnames(cell_metadata), "\n")
  
  # Create CellChat object with explicit group.by parameter and suppress samples warning
  cellchat <- createCellChat(object = expression_data, meta = cell_metadata, group.by = "cell_identity")
  
  # Add a samples column to avoid the warning (though it doesn't affect functionality)
  cellchat@meta$samples <- sample_name
  
  # Set cell type identity (already set via group.by above, but keeping for clarity)
  cellchat <- setIdent(cellchat, ident.use = "cell_identity")
  
  # Set database (choose human or mouse)
  CellChatDB <- CellChatDB.human  # or CellChatDB.mouse
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
    cellchat <- projectData(cellchat, PPI.human)  # or PPI.mouse
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
# STEP 2: PROCESS ALL SAMPLES
# ==============================================================================

# Create list of Seurat objects
seurat_objects <- list(
  "C1" = C1,
  "C2" = C2, 
  "C4" = C4,
  "C5" = C5,
  "P1" = P1,
  "P2" = P2,
  "P3" = P3,
  "P4" = P4,
  "P5" = P5
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
# STEP 3: ORGANIZE INTO CONTROL AND PAH GROUPS
# ==============================================================================

# Separate into control and PAH groups
control_cellchats <- list(
  "C1" = cellchat_objects[["C1"]],
  "C2" = cellchat_objects[["C2"]],
  "C4" = cellchat_objects[["C4"]],
  "C5" = cellchat_objects[["C5"]]
)

pah_cellchats <- list(
  "P1" = cellchat_objects[["P1"]],
  "P2" = cellchat_objects[["P2"]],
  "P3" = cellchat_objects[["P3"]],
  "P4" = cellchat_objects[["P4"]],
  "P5" = cellchat_objects[["P5"]]
)

# Verify the groups
cat("=== VERIFICATION ===\n")
cat("Control samples:", names(control_cellchats), "\n")
cat("PAH samples:", names(pah_cellchats), "\n")

# Check cell types in each sample (should be consistent)
cat("\nCell types in each sample:\n")
for(sample_name in names(cellchat_objects)) {
  cell_types <- levels(cellchat_objects[[sample_name]]@idents)
  cat(sample_name, ":", paste(cell_types, collapse = ", "), "\n")
}

# ==============================================================================
# STEP 4: STATISTICAL ANALYSIS FUNCTIONS
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

# Main statistical analysis function
analyze_pathways_with_stats <- function(control_list, pah_list, cell_types, cell_name) {
  
  cat("\n=== Statistical Analysis for", cell_name, "===\n")
  
  # Calculate pathway strengths for each sample
  control_strengths <- lapply(control_list, function(obj) {
    calculate_celltype_pathways(obj, cell_types)
  })
  
  pah_strengths <- lapply(pah_list, function(obj) {
    calculate_celltype_pathways(obj, cell_types)
  })
  
  # Remove NULL results
  control_strengths <- control_strengths[!sapply(control_strengths, is.null)]
  pah_strengths <- pah_strengths[!sapply(pah_strengths, is.null)]
  
  if (length(control_strengths) == 0 || length(pah_strengths) == 0) {
    cat("No valid data for", cell_name, "\n")
    return(NULL)
  }
  
  # Get common pathways across all samples
  all_pathway_names <- Reduce(intersect, lapply(c(control_strengths, pah_strengths), names))
  
  # Create matrices for statistical testing
  control_matrix <- do.call(cbind, lapply(control_strengths, function(x) x[all_pathway_names]))
  pah_matrix <- do.call(cbind, lapply(pah_strengths, function(x) x[all_pathway_names]))
  
  # Perform statistical tests for each pathway
  results_list <- list()
  
  for (pathway in all_pathway_names) {
    control_vals <- control_matrix[pathway, ]
    pah_vals <- pah_matrix[pathway, ]
    
    # Calculate basic statistics
    control_mean <- mean(control_vals, na.rm = TRUE)
    pah_mean <- mean(pah_vals, na.rm = TRUE)
    control_sd <- sd(control_vals, na.rm = TRUE)
    pah_sd <- sd(pah_vals, na.rm = TRUE)
    
    # Perform statistical test
    test_result <- tryCatch({
      wilcox.test(control_vals, pah_vals, exact = FALSE)
    }, error = function(e) {
      t.test(control_vals, pah_vals)
    })
    
    # Calculate effect size (Cohen's d)
    pooled_sd <- sqrt(((length(control_vals)-1) * control_sd^2 + (length(pah_vals)-1) * pah_sd^2) / 
                        (length(control_vals) + length(pah_vals) - 2))
    cohens_d <- (pah_mean - control_mean) / pooled_sd
    
    results_list[[pathway]] <- data.frame(
      Pathway = pathway,
      Control_Mean = round(control_mean, 4),
      Control_SD = round(control_sd, 4),
      PAH_Mean = round(pah_mean, 4),
      PAH_SD = round(pah_sd, 4),
      Difference = round(pah_mean - control_mean, 4),
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
  cat("\nTop 5 pathways for", cell_name, ":\n")
  print(head(results_df[, c("Pathway", "Control_Mean", "PAH_Mean", "Difference", 
                            "P_Value", "FDR", "Cohens_D", "Significance")], 5))
  
  # Create visualization - top 5 pathways
  top_pathways <- head(results_df, 5)
  
  plot_df <- data.frame(
    pathway = factor(top_pathways$Pathway, levels = rev(top_pathways$Pathway)),
    difference = top_pathways$Difference,
    significance = top_pathways$Significance,
    fdr = top_pathways$FDR
  )
  
  p <- ggplot(plot_df, aes(x = pathway, y = difference)) +
    geom_col(aes(fill = difference > 0), width = 0.7) +
    geom_text(aes(label = significance), hjust = ifelse(plot_df$difference > 0, -0.1, 1.1), 
              size = 4, fontface = "bold") +
    scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
                      labels = c("FALSE" = paste("Decreased in PAH", cell_name), 
                                 "TRUE" = paste("Increased in PAH", cell_name))) +
    coord_flip() +
    labs(title = paste("Top 5", cell_name, "Pathways (PAH vs Control)"),
         subtitle = "*** p<0.001, ** p<0.01, * p<0.05, . p<0.1",
         x = "Pathway",
         y = "Mean Difference (PAH - Control)",
         fill = "Direction") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  print(p)
  
  return(results_df)
}

# ==============================================================================
# STEP 5: RUN ANALYSIS FOR ALL CELL TYPES
# ==============================================================================

# Define cell types to analyze (focused on main liver cell types)
cell_types_to_analyze <- list(
  "Hepatocytes" = "Hepatocyte",
  "Endothelial" = "Endothelial",
  "HSCs" = "HSC/mFB", 
  "Macrophages" = "Macrophage"
)

# Store results
all_results <- list()

# Run analysis for each cell type
for (cell_name in names(cell_types_to_analyze)) {
  all_results[[cell_name]] <- analyze_pathways_with_stats(
    control_list = control_cellchats,
    pah_list = pah_cellchats,
    cell_types = cell_types_to_analyze[[cell_name]],
    cell_name = cell_name
  )
}

# Save the CellChat objects for future use
save(cellchat_objects, control_cellchats, pah_cellchats, all_results,
     file = "individual_cellchat_analysis.RData")

