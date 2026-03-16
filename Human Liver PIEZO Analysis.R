# ══════════════════════════════════════════════════════════════════════════════
# PIEZO1/2 EXPRESSION IN HSC/mFB - PAH vs. CONTROL
# ══════════════════════════════════════════════════════════════════════════════

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# ── Subset HSC/mFB from PAH dataset ───────────────────────────────────────────
Idents(PAH) <- "cell_identity"
HSC_PAH <- subset(PAH, idents = "HSC/mFB")

# ── Confirm genes are present ──────────────────────────────────────────────────
piezo_genes <- c("PIEZO1", "PIEZO2")
piezo_found <- intersect(piezo_genes, rownames(HSC_PAH))
cat("Piezo genes found:", piezo_found, "\n")

# ── Add condition label ────────────────────────────────────────────────────────
# Adjust this pattern to match your sample naming convention
HSC_PAH$condition <- ifelse(grepl("^C|Control", HSC_PAH$sample, ignore.case = TRUE),
                            "Control", "PAH")
HSC_PAH$condition <- factor(HSC_PAH$condition, levels = c("Control", "PAH"))

# ══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION
# ══════════════════════════════════════════════════════════════════════════════

# ── 1. Violin plots ────────────────────────────────────────────────────────────
p_vln <- VlnPlot(
  HSC_PAH,
  features = piezo_found,
  group.by = "condition",
  cols = c("Control" = "#4DAFED", "PAH" = "#E84B36"),
  pt.size = 0,
  ncol = 2
) & theme(plot.title = element_text(face = "bold"))

# ── 2. Feature plots on UMAP ──────────────────────────────────────────────────
p_feat <- FeaturePlot(
  HSC_PAH,
  features = piezo_found,
  ncol = 2,
  order = TRUE,
  cols = c("lightgrey", "darkred")
)

# ── 3. Dot plot (expression level + % expressing) ─────────────────────────────
Idents(HSC_PAH) <- "condition"
p_dot <- DotPlot(
  HSC_PAH,
  features = piezo_found,
  cols = c("lightgrey", "darkred")
) + 
  RotatedAxis() +
  ggtitle("PIEZO1/2 - HSC/mFB: Control vs. PAH") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

print(p_vln)
print(p_feat)
print(p_dot)

# ══════════════════════════════════════════════════════════════════════════════
# PSEUDOBULK STATS - t-test per gene
# ══════════════════════════════════════════════════════════════════════════════

# Pseudobulk by sample
pseudo_HSC <- AggregateExpression(HSC_PAH, return.seurat = TRUE, group.by = "sample")

# Assign condition
pseudo_HSC$condition <- ifelse(grepl("^C|Control", rownames(pseudo_HSC@meta.data), ignore.case = TRUE),
                               "Control", "PAH")

# Extract normalized expression for PIEZO1/2
expr_df <- as.data.frame(t(GetAssayData(pseudo_HSC, assay = "RNA", layer = "data")[piezo_found, , drop = FALSE]))
expr_df$condition <- pseudo_HSC$condition
expr_df$sample <- rownames(expr_df)

cat("\n=== PSEUDOBULK EXPRESSION SUMMARY ===\n")
print(expr_df)

# t-tests
for (gene in piezo_found) {
  ctrl_vals <- expr_df[expr_df$condition == "Control", gene]
  pah_vals  <- expr_df[expr_df$condition == "PAH", gene]
  
  cat(paste0("\n── ", gene, " ──\n"))
  cat("Control mean:", round(mean(ctrl_vals), 4), "\n")
  cat("PAH mean:    ", round(mean(pah_vals), 4), "\n")
  
  if (length(unique(c(ctrl_vals, pah_vals))) > 1) {
    tt <- t.test(pah_vals, ctrl_vals)
    cat("t-test p-value:", round(tt$p.value, 4), "\n")
    cat("Fold change (PAH/Control):", round(mean(pah_vals) / mean(ctrl_vals), 3), "\n")
  } else {
    cat("Insufficient variance for t-test\n")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# PSEUDOBULK SINGLE GENE EXPRESSION - PIEZO1 in HSC/mFB
# ══════════════════════════════════════════════════════════════════════════════

library(Seurat)
library(dplyr)
library(ggplot2)

# ── Create pseudobulk HSC object ───────────────────────────────────────────────
Idents(PAH) <- "cell_identity"
HSC <- subset(PAH, idents = "HSC/mFB")
pseudo_HSC <- AggregateExpression(HSC, return.seurat = TRUE, group.by = "sample")

# ── Assign condition ───────────────────────────────────────────────────────────
# Check order first!
cat("Sample order:\n")
print(rownames(pseudo_HSC@meta.data))

pseudo_HSC$condition <- ifelse(grepl("^C|Control", rownames(pseudo_HSC@meta.data), ignore.case = TRUE),
                               "Control", "PAH")
cat("\nCondition assignment:\n")
print(data.frame(sample = rownames(pseudo_HSC@meta.data), condition = pseudo_HSC$condition))

# ── Extract PIEZO1 expression ──────────────────────────────────────────────────
PIEZO1_expr <- FetchData(pseudo_HSC, vars = "PIEZO1", layer = "data")

HSC_PIEZO1 <- data.frame(
  sample    = rownames(pseudo_HSC@meta.data),
  condition = factor(pseudo_HSC$condition, levels = c("Control", "PAH")),
  PIEZO1_expression = PIEZO1_expr$PIEZO1
)

print(HSC_PIEZO1)

# ── Statistics ─────────────────────────────────────────────────────────────────
t_PIEZO1 <- t.test(PIEZO1_expression ~ condition, data = HSC_PIEZO1)
cat("\nPIEZO1 t-test p-value:", signif(t_PIEZO1$p.value, 3), "\n")

PIEZO1_summary <- HSC_PIEZO1 %>%
  group_by(condition) %>%
  summarise(
    n      = n(),
    mean   = mean(PIEZO1_expression),
    sd     = sd(PIEZO1_expression),
    se     = sd / sqrt(n),
    .groups = "drop"
  )

print(PIEZO1_summary)

# ── Plot ───────────────────────────────────────────────────────────────────────
ggplot(HSC_PIEZO1, aes(x = condition, y = PIEZO1_expression, fill = condition)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, colour = "black") +
  geom_errorbar(data = PIEZO1_summary,
                aes(x = condition, y = mean, ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(size = 3, shape = 21, fill = "white", colour = "black") +
  annotate("text", x = 1.5,
           y = max(HSC_PIEZO1$PIEZO1_expression) * 1.15,
           label = paste0("p = ", signif(t_PIEZO1$p.value, 3)),
           size = 4.5) +
  scale_fill_manual(values = c("Control" = "#4DAFED", "PAH" = "#E84B36")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  labs(
    y     = "PIEZO1 Expression (Pseudobulk)",
    title = "HSC/mFB: PIEZO1 Expression (PAH vs. Control)"
  )

# ── Export ─────────────────────────────────────────────────────────────────────
write.csv(HSC_PIEZO1, "PAH_HSC_PIEZO1_Pseudobulk_Expression.csv", row.names = FALSE)