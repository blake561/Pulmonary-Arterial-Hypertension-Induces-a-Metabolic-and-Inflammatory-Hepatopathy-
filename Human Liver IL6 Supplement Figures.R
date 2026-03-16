

hsc_PAH <- subset(PAH, idents ='HSC/mFB')
hsc_NASH <- subset(NASH, idents ='hep/mFB')
Fontan_HSC <- subset(Fontan, idents='HSC/mFB')

install.packages("ggtext")

# ── Load Required Libraries ───────────────────────────────────────────────
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggtext)  # For subscript formatting

# ── Function to Compute Pseudobulk Log2 Fold Change ───────────────────────
compute_log2fc_pseudobulk <- function(seurat_obj, group_var, control_label, case_label, gene, disease_name) {
  DefaultAssay(seurat_obj) <- "RNA"
  counts <- GetAssayData(seurat_obj, slot = "counts")[gene, ]
  meta <- seurat_obj@meta.data
  
  df <- data.frame(
    group = meta[[group_var]],
    counts = as.numeric(counts)
  ) %>%
    filter(group %in% c(control_label, case_label))
  
  # Pseudobulk sum
  pseudobulk <- df %>%
    group_by(group) %>%
    summarise(total = sum(counts), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = total, values_fill = 0)
  
  log2fc <- log2((pseudobulk[[case_label]] + 1) / (pseudobulk[[control_label]] + 1))
  
  # Wilcoxon p-value (skip if one group is all zeros or empty)
  counts_case <- df$counts[df$group == case_label]
  counts_ctrl <- df$counts[df$group == control_label]
  
  p_val <- tryCatch({
    if (length(unique(counts_case)) > 1 || length(unique(counts_ctrl)) > 1) {
      wilcox.test(counts_case, counts_ctrl)$p.value
    } else {
      NA
    }
  }, error = function(e) NA)
  
  tibble(
    disease = disease_name,
    log2FC = log2fc,
    p = p_val,
    FDR = ifelse(is.na(p_val), NA, p.adjust(p_val, method = "fdr"))
  )
}

# ── Run for Each Disease Group: HSCs only ─────────────────────────────────
il6_fc_df <- bind_rows(
  compute_log2fc_pseudobulk(hsc_PAH,    "condition",    "Control",       "PAH",           "IL6", "PAH"),
  compute_log2fc_pseudobulk(hsc_NASH,   "condition",    "Control",       "NASH",          "IL6", "NASH"),
  compute_log2fc_pseudobulk(Fontan_HSC, "orig.ident",   "Fontan Control","Fontan",        "IL6", "Fontan")
) %>%
  mutate(
    disease = factor(disease, levels = c("PAH", "NASH", "Fontan")),
    fdr_label = ifelse(FDR < 0.001, 
                       "P<sub>adj</sub> < 0.001", 
                       paste0("P<sub>adj</sub> = ", signif(FDR, 3)))
  )

# ── Define color palette (matching module score plots) ────────────────────
disease_colors <- c(
  "PAH"    = "steelblue",
  "NASH"   = "orchid4",
  "Fontan" = "olivedrab"
)

# ── Plot ───────────────────────────────────────────────────────────────────
ggplot(il6_fc_df, aes(x = disease, y = log2FC, fill = disease)) +
  geom_col(color = "black", linewidth = 0.5, width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_richtext(
    aes(label = fdr_label, y = log2FC + ifelse(log2FC > 0, 0.3, -0.3)),
    size = 4,
    fontface = "bold",
    fill = NA,
    label.color = NA
  ) +
  scale_fill_manual(values = disease_colors) +
  labs(
    title = "IL6 Log2 Fold Change in HSCs (Disease vs. Control)",
    x = "",
    y = "Log2 Fold Change"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# ── Load Required Libraries ─────────────────────────────────────────────────────
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(ggtext)

# ── Function to Compute True Pseudobulk Log2 Fold Change per Cell Type ─────────
compute_log2fc_all_subtypes <- function(seurat_obj, group_var, control_label, case_label, gene) {
  DefaultAssay(seurat_obj) <- "RNA"
  counts <- GetAssayData(seurat_obj, slot = "counts")[gene, ]
  meta <- seurat_obj@meta.data
  
  df <- data.frame(
    cell_identity = meta$cell_identity,
    group = meta[[group_var]],
    counts = as.numeric(counts)
  ) %>%
    filter(group %in% c(control_label, case_label))
  
  # Pseudobulk sums
  fc_df <- df %>%
    group_by(cell_identity, group) %>%
    summarise(total = sum(counts), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = total, values_fill = 0) %>%
    mutate(log2fc = log2((.data[[case_label]] + 1) / (.data[[control_label]] + 1)))
  
  # P-value computation
  pvals <- df %>%
    group_by(cell_identity) %>%
    summarise(
      p = tryCatch(
        wilcox.test(counts[group == case_label], counts[group == control_label])$p.value,
        error = function(e) NA
      ),
      .groups = "drop"
    ) %>%
    mutate(fdr = p.adjust(p, method = "fdr"))
  
  left_join(fc_df, pvals, by = "cell_identity")
}

# ── Run for IL6: PAH vs Control ────────────────────────────────────────────────
log2fc_result <- compute_log2fc_all_subtypes(
  seurat_obj     = PAH,
  group_var      = "condition",
  control_label  = "Control",
  case_label     = "PAH",
  gene           = "IL6"
) %>%
  mutate(
    # Reorder cell_identity by log2fc (descending)
    cell_identity = factor(cell_identity, levels = cell_identity[order(-log2fc)]),
    # Create label for P_adj
    fdr_label = case_when(
      is.na(fdr) ~ "P<sub>adj</sub> = NaN",
      fdr < 0.001 ~ "P<sub>adj</sub> < 0.001",
      TRUE ~ paste0("P<sub>adj</sub> = ", signif(fdr, 2))
    ),
    # Determine fill color
    direction = ifelse(log2fc > 0, "positive", "negative")
  )

# ── Plot ────────────────────────────────────────────────────────────────────────
ggplot(log2fc_result, aes(x = cell_identity, y = log2fc, fill = direction)) +
  geom_col(color = "black", linewidth = 0.5, width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_richtext(
    aes(label = fdr_label, y = log2fc + ifelse(log2fc > 0, 0.3, -0.3)),
    size = 3.5,
    fontface = "bold",
    fill = NA,
    label.color = NA
  ) +
  scale_fill_manual(values = c("positive" = "#E74C3C", "negative" = "#4A90E2")) +
  labs(
    title = "IL6 Log<sub>2</sub> Fold Change (PAH vs. Control)",
    x = "Cell Type",
    y = "Log<sub>2</sub> Fold Change"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.title.y = element_markdown(),
    plot.title = element_markdown(hjust = 0.5, face = "bold")
  )