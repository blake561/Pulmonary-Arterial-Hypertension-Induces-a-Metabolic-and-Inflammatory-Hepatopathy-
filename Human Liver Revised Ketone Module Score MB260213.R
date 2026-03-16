# ══════════════════════════════════════════════════════════════════════════════
# HEPATOCYTE MODULE SCORING - KETONE METABOLISM
# Comparing PAH, NASH, and FALD liver datasets
# ══════════════════════════════════════════════════════════════════════════════

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ══════════════════════════════════════════════════════════════════════════════
# 1. LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════

PAH <- readRDS(file = '/projects/standard/prin0088/shared/Madelyn Liver and Kidney Projects/PAH_total.mtx_final')
NASH <- readRDS("/users/1/blake561/NASH_total.mtx_final")  
FALD <- readRDS(file = '/projects/standard/prin0088/shared/Madelyn Liver and Kidney Projects/Fontan_total.mtx_final')

# ══════════════════════════════════════════════════════════════════════════════
# 2. PREPARE HEPATOCYTE SUBSETS
# ══════════════════════════════════════════════════════════════════════════════

# PAH - uses cell_identity
DefaultAssay(PAH) <- "RNA"
Idents(PAH) <- "cell_identity"
Hepatocyte <- subset(PAH, idents = "Hepatocyte")

# NASH - uses subcelltype  
DefaultAssay(NASH) <- "RNA"
Idents(NASH) <- "subcelltype"
hep_NASH <- subset(NASH, idents = "Hepatocyte")

# Fontan/FALD - uses cell_identity
DefaultAssay(FALD) <- "RNA"
Idents(FALD) <- "cell_identity"
Fontan_Hepatocyte <- subset(FALD, idents = "Hepatocyte")

# ══════════════════════════════════════════════════════════════════════════════
# 3. DEFINE KETONE METABOLISM GENE SET
# ══════════════════════════════════════════════════════════════════════════════

Ketone <- c(
  "ACAT1",   # acetyl-CoA acetyltransferase 1
  "ACAT2",   # acetyl-CoA acetyltransferase 2
  "BDH1",    # 3-hydroxybutyrate dehydrogenase, type 1
  "BDH2",    # 3-hydroxybutyrate dehydrogenase, type 2
  "HMGCL",   # 3-hydroxymethyl-3-methylglutaryl-CoA lyase
  "HMGCS1",  # 3-hydroxy-3-methylglutaryl-CoA synthase 1 (soluble)
  "HMGCS2",  # 3-hydroxy-3-methylglutaryl-CoA synthase 2 (mitochondrial)
  "OXCT1",   # 3-oxoacid CoA transferase 1
  "OXCT2"    # 3-oxoacid CoA transferase 2
)

# ══════════════════════════════════════════════════════════════════════════════
# 4. PLOTTING FUNCTION WITH MEDIAN + IQR ERROR BARS
# ══════════════════════════════════════════════════════════════════════════════

plot_module_score <- function(seurat_obj, gene_list, pathway_name, 
                              display_name, dataset_name, 
                              ctrl_color, disease_color,
                              condition_col = "condition",
                              ctrl_label = "Control",
                              disease_label = NULL) {
  
  if (is.null(disease_label)) {
    disease_label <- dataset_name
  }
  
  # Add module score
  seurat_obj <- AddModuleScore(seurat_obj, list(gene_list), 
                               name = paste0(pathway_name, "_"), assay = "RNA")
  score_col <- paste0(pathway_name, "_1")
  
  # Extract and format data
  df <- seurat_obj@meta.data %>%
    transmute(
      score = .data[[score_col]],
      condition = ifelse(.data[[condition_col]] == ctrl_label, "Control", disease_label)
    ) %>%
    filter(!is.na(condition)) %>%
    mutate(condition = factor(condition, levels = c("Control", disease_label)))
  
  # Statistical test
  vals_ctrl <- df$score[df$condition == "Control"]
  vals_dis <- df$score[df$condition == disease_label]
  
  if (length(vals_ctrl) >= 3 && length(vals_ctrl) <= 5000 &&
      length(vals_dis) >= 3 && length(vals_dis) <= 5000) {
    sw_ctrl <- shapiro.test(vals_ctrl)$p.value
    sw_dis <- shapiro.test(vals_dis)$p.value
    normal <- sw_ctrl > 0.05 && sw_dis > 0.05
  } else {
    normal <- FALSE
  }
  
  if (normal) {
    p_val <- t.test(vals_ctrl, vals_dis)$p.value
  } else {
    p_val <- wilcox.test(vals_ctrl, vals_dis)$p.value
  }
  
  if (p_val < 0.001) {
    p_label <- "p < 0.001"
  } else {
    p_label <- paste0("p = ", formatC(p_val, format = "e", digits = 2))
  }
  
  # Calculate statistics: median, IQR for error bars
  stats_df <- df %>%
    group_by(condition) %>%
    summarise(
      median_score = median(score, na.rm = TRUE),
      min_score = min(score, na.rm = TRUE),
      max_score = max(score, na.rm = TRUE),
      q25 = quantile(score, 0.25, na.rm = TRUE),
      q75 = quantile(score, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  y_max <- max(df$score, na.rm = TRUE)
  y_min <- min(df$score, na.rm = TRUE)
  y_range <- y_max - y_min
  
  colors <- setNames(c(ctrl_color, disease_color), c("Control", disease_label))
  
  # Create plot
  p <- ggplot(df, aes(x = condition, y = score, fill = condition)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.9) +
    
    # Median crossbar
    stat_summary(fun = median, geom = "crossbar", width = 0.5, 
                 color = "black", fatten = 1, linewidth = 0.5) +
    
    # Error bars showing IQR
    geom_errorbar(data = stats_df,
                  aes(x = condition, y = median_score, 
                      ymin = q25, ymax = q75),
                  width = 0.15, linewidth = 0.8, color = "black",
                  inherit.aes = FALSE) +
    
    # Median value labels
    geom_text(data = stats_df,
              aes(x = condition, y = y_max + y_range * 0.08,
                  label = sprintf("%.3f", median_score)),
              inherit.aes = FALSE, size = 3.5, fontface = "plain") +
    
    # P-value annotation
    annotate("text", x = 1.5, y = y_max + y_range * 0.18,
             label = p_label, size = 3.5) +
    
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    labs(
      y = "Module Score",
      title = display_name
    )
  
  return(p)
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. GENERATE KETONE METABOLISM PLOTS
# ══════════════════════════════════════════════════════════════════════════════

cat("Processing: Ketone Metabolism\n")

# PAH
p_PAH <- plot_module_score(
  seurat_obj = Hepatocyte,
  gene_list = Ketone,
  pathway_name = "Ketone",
  display_name = "Ketone Metabolism",
  dataset_name = "PAH",
  ctrl_color = "#FA8072",
  disease_color = "#4682B4",
  condition_col = "condition",
  ctrl_label = "Control",
  disease_label = "PAH"
)

# NASH
p_NASH <- plot_module_score(
  seurat_obj = hep_NASH,
  gene_list = Ketone,
  pathway_name = "Ketone",
  display_name = "Ketone Metabolism",
  dataset_name = "NASH",
  ctrl_color = "#DAA520",
  disease_color = "#9932CC",
  condition_col = "condition",
  ctrl_label = "Control",
  disease_label = "NASH"
)

# FALD
p_FALD <- plot_module_score(
  seurat_obj = Fontan_Hepatocyte,
  gene_list = Ketone,
  pathway_name = "Ketone",
  display_name = "Ketone Metabolism",
  dataset_name = "FALD",
  ctrl_color = "#A0522D",
  disease_color = "#6B8E23",
  condition_col = "orig.ident",
  ctrl_label = "Fontan Control",
  disease_label = "FALD"
)

# ══════════════════════════════════════════════════════════════════════════════
# 6. ARRANGE & SAVE - 1x3 LAYOUT (PAH | NASH | FALD)
# ══════════════════════════════════════════════════════════════════════════════

final_figure <- p_PAH | p_NASH | p_FALD

ggsave("Hepatocyte_Ketone_Metabolism.pdf", final_figure, width = 10, height = 4, dpi = 300)
ggsave("Hepatocyte_Ketone_Metabolism.png", final_figure, width = 10, height = 4, dpi = 300)

cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
cat("DONE: Ketone Metabolism module scores saved!\n")
cat("  - Hepatocyte_Ketone_Metabolism.pdf\n")
cat("  - Hepatocyte_Ketone_Metabolism.png\n")
cat("═══════════════════════════════════════════════════════════════════════════════\n")