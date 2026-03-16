
##Install Packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
BiocManager::install("DESeq2")

##Load Packages

library(Seurat)
library(ggplot2)
library(SeuratData)
library(EnhancedVolcano)
library(DESeq2)
library(BiocManager)
library(tidyr)
library(gridExtra)
library(tidyverse)

##Load Data (From Integration and Clustering Sheet)

total.mtx_deg <- readRDS("/users/1/blake561/total.mtx_final")

total.mtx_deg$annotation <- Idents(total.mtx_deg)

total.mtx_deg$sample_num <- paste0(total.mtx_deg$cell_identity, total.mtx_deg$sample)


#Creating Matrix for DEG

cts <- AggregateExpression(total.mtx_deg,
                           group.by = c("annotation", "sample_num"),
                           assays = "RNA",
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$RNA
cts.t <- t(cts)
cts.t <- as.data.frame(cts.t)
splitRows <- gsub('_.*', '', rownames(cts.t))

cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

## Differential Analysis - Cell Types

#Hepatocyte DEGs (Whole Cluster)

counts_HP <- cts.split.modified$Hepatocyte
colData <- data.frame(samples = colnames(counts_HP))

colData$condition <- c('Control','Control','Control','Control','PAH','PAH','PAH','PAH','PAH')

colData <- column_to_rownames(colData, var = 'samples')

dds_HP <- DESeqDataSetFromMatrix(countData = counts_HP,
                                 colData = colData,
                                 design = ~condition)

keep <- rowSums(counts(dds_HP)) >=10
dds_HP <- dds_HP[keep,]

dds_HP <- DESeq(dds_HP)

#Hepatocyte Comparision 

#These comparisons are c("condition", "numerator", "denominator")

PAH_vs_Control_HP_res <- results(dds_HP, contrast =c("condition", "PAH", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in PAH Relative to Control 

#Saving File for Pathway Analysis

write.csv(PAH_vs_Control_HP_res, file = "MB250115 Pseudobulk HP DE PAH_vs_Control")
d1 <- as.data.frame(PAH_vs_Control_HP_res)
d2 <- subset(d1, subset = padj <= 0.05 &
               abs(log2FoldChange) >= 0.25) #padj is FDR
write.csv(d2, file = "PAH_vs_Control Hepatocyte Pathway Analysis.csv")

#Volcano Plot

EnhancedVolcano(PAH_vs_Control_HP_res,
                lab = rownames(PAH_vs_Control_HP_res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = "PAH vs Control Hepatocytes")

#Macrophage Comparision 

counts_MP <- cts.split.modified$Macrophage
colData <- data.frame(samples = colnames(counts_MP))

colData$condition <- c('Control','Control','Control','Control','PAH','PAH','PAH','PAH','PAH')
#View(colData)

colData <- column_to_rownames(colData, var = 'samples')

dds_MP <- DESeqDataSetFromMatrix(countData = counts_MP,
                                 colData = colData,
                                 design = ~condition)

keep <- rowSums(counts(dds_MP)) >=10
dds_MP <- dds_MP[keep,]

dds_MP <- DESeq(dds_MP)
resultsNames(dds_MP)

#These comparisons are c("condition", "numerator", "denominator")

PAH_vs_Control_MP_res <- results(dds_MP, contrast =c("condition", "PAH", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in PAH Relative to Control 
view(PAH_vs_Control_MP_res)

#Saving File for Pathway Analysis

##Macrophage DEGs (Whole Cluster)

counts_MP <- cts.split.modified$Macrophage
colData <- data.frame(samples = colnames(counts_MP))

colData$condition <- c('Control','Control','Control','Control','PAH','PAH','PAH','PAH','PAH')
#View(colData)

colData <- column_to_rownames(colData, var = 'samples')

dds_MP <- DESeqDataSetFromMatrix(countData = counts_MP,
                                 colData = colData,
                                 design = ~condition)

keep <- rowSums(counts(dds_MP)) >=10
dds_MP <- dds_MP[keep,]

dds_MP <- DESeq(dds_MP)
resultsNames(dds_MP)

#Comparison 

PAH_vs_Control_MP_res <- results(dds_MP, contrast =c("condition", "PAH", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in PAH Relative to Control 
view(PAH_vs_Control_MP_res)

#Saving File

write.csv(PAH_vs_Control_MP_res, file = "MB250123 Pseudobulk MP DE PAH_vs_Control")
m1 <- as.data.frame(PAH_vs_Control_MP_res)
m2 <- subset(m1, subset = padj <= 0.05 &
               abs(log2FoldChange) >= 0.25) #padj is FDR
write.csv(m2, file = "PAH_vs_Control Macrophage Pathway Analysis.csv")

#Volcano Plot

EnhancedVolcano(PAH_vs_Control_MP_res,
                lab = rownames(PAH_vs_Control_MP_res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = "PAH vs Control Macrophage")

#Lymphocyte Comparision 

##Lymphocyte DEGs (Whole Cluster)

counts_LC <- cts.split.modified$Lymphocyte
colData <- data.frame(samples = colnames(counts_LC))

colData$condition <- c('Control','Control','Control','Control','PAH','PAH','PAH','PAH','PAH')

colData <- column_to_rownames(colData, var = 'samples')

dds_LC <- DESeqDataSetFromMatrix(countData = counts_LC,
                                 colData = colData,
                                 design = ~condition)

keep <- rowSums(counts(dds_LC)) >=10
dds_LC <- dds_LC[keep,]

dds_LC <- DESeq(dds_LC)
resultsNames(dds_LC)

#Comparison 

PAH_vs_Control_LC_res <- results(dds_LC, contrast =c("condition", "PAH", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in PAH Relative to Control 
view(PAH_vs_Control_LC_res)

#Saving File

write.csv(PAH_vs_Control_LC_res, file = "MB250115 Pseudobulk LC DE PAH_vs_Control")
l1 <- as.data.frame(PAH_vs_Control_LC_res)
l2 <- subset(l1, subset = padj <= 0.05 &
               abs(log2FoldChange) >= 0.25) #padj is FDR
write.csv(l2, file = "PAH_vs_Control Lymphocyte Pathway Analysis.csv")

#Volcano Plot

EnhancedVolcano(PAH_vs_Control_LC_res,
                lab = rownames(PAH_vs_Control_LC_res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = "PAH vs Control Lymphocyte")

##Cholangiocytes DEGs (Whole Cluster)

counts_CL <- cts.split.modified$Cholangiocyte
colData <- data.frame(samples = colnames(counts_CL))

colData$condition <- c('Control','Control','Control','Control','PAH','PAH','PAH','PAH','PAH')

colData <- column_to_rownames(colData, var = 'samples')

dds_CL <- DESeqDataSetFromMatrix(countData = counts_CL,
                                 colData = colData,
                                 design = ~condition)

keep <- rowSums(counts(dds_CL)) >=10
dds_CL <- dds_CL[keep,]

dds_CL <- DESeq(dds_CL)
resultsNames(dds_CL)

#Comparison 

PAH_vs_Control_CL_res <- results(dds_CL, contrast =c("condition", "PAH", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in PAH Relative to Control 
#view(PAH_vs_Control_CL_res)

#Saving File

write.csv(PAH_vs_Control_CL_res, file = "MB250115 Pseudobulk CL DE PAH_vs_Control")
c1 <- as.data.frame(PAH_vs_Control_CL_res)
c2 <- subset(c1, subset = padj <= 0.05 &
               abs(log2FoldChange) >= 0.25) #padj is FDR
write.csv(c2, file = "PAH_vs_Control Cholangiocyte Pathway Analysis.csv")

#Volcano Plot

EnhancedVolcano(PAH_vs_Control_CL_res,
                lab = rownames(PAH_vs_Control_CL_res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = "PAH vs Control Cholangiocyte")

##Endothelial DEGs (Whole Cluster)

counts_EC <- cts.split.modified$Endothelial
colData <- data.frame(samples = colnames(counts_EC))

colData$condition <- c('Control','Control','Control','Control','PAH','PAH','PAH','PAH','PAH')

colData <- column_to_rownames(colData, var = 'samples')

dds_EC <- DESeqDataSetFromMatrix(countData = counts_EC,
                                 colData = colData,
                                 design = ~condition)

keep <- rowSums(counts(dds_EC)) >=10
dds_EC <- dds_EC[keep,]

dds_EC <- DESeq(dds_EC)
resultsNames(dds_EC)

#Comparison 

PAH_vs_Control_EC_res <- results(dds_EC, contrast =c("condition", "PAH", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in PAH Relative to Control 

#Saving File

write.csv(PAH_vs_Control_EC_res, file = "MB250115 Pseudobulk EC DE PAH_vs_Control")
e1 <- as.data.frame(PAH_vs_Control_EC_res)
e2 <- subset(e1, subset = padj <= 0.05 &
               abs(log2FoldChange) >= 0.25) #padj is FDR
write.csv(e2, file = "PAH_vs_Control Endothelial Pathway Analysis.csv")

#Volcano Plot

EnhancedVolcano(PAH_vs_Control_EC_res,
                lab = rownames(PAH_vs_Control_EC_res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = "PAH vs Control Endothelial")

##HSC/mFB DEGs (Whole Cluster)

counts_FB <- cts.split.modified$`HSC/mFB`
colData <- data.frame(samples = colnames(counts_FB))

colData$condition <- c('Control','Control','Control','Control','PAH','PAH','PAH','PAH','PAH')

colData <- column_to_rownames(colData, var = 'samples')

dds_FB <- DESeqDataSetFromMatrix(countData = counts_FB,
                                 colData = colData,
                                 design = ~condition)

keep <- rowSums(counts(dds_FB)) >=10
dds_FB <- dds_FB[keep,]

dds_FB <- DESeq(dds_FB)
resultsNames(dds_FB)

#Comparison 

PAH_vs_Control_FB_res <- results(dds_FB, contrast =c("condition", "PAH", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in PAH Relative to Control 

#Saving File

write.csv(PAH_vs_Control_FB_res, file = "MB250115 Pseudobulk FB DE PAH_vs_Control")
f1 <- as.data.frame(PAH_vs_Control_FB_res)
f2 <- subset(f1, subset = padj <= 0.05 &
               abs(log2FoldChange) >= 0.25) #padj is FDR
write.csv(f2, file = "PAH_vs_Control HSC mFB Pathway Analysis.csv")

#Volcano Plot

EnhancedVolcano(PAH_vs_Control_FB_res,
                lab = rownames(PAH_vs_Control_FB_res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = "PAH vs Control HSC/mFB")

#Plasma/Erythrocyte DEGs (Whole Cluster)

counts_PE <- cts.split.modified$`Plasma/Erythrocyte`
colData <- data.frame(samples = colnames(counts_PE))

view(coldata)

colData$condition <- c('Control','Control','Control','PAH','PAH','PAH','PAH','PAH') #Only 3 Control Samples Have P/E - Not Sure Why

colData <- column_to_rownames(colData, var = 'samples')

dds_PE <- DESeqDataSetFromMatrix(countData = counts_PE,
                                 colData = colData,
                                 design = ~condition)

keep <- rowSums(counts(dds_PE)) >=10
dds_PE <- dds_PE[keep,]

dds_PE <- DESeq(dds_PE)
resultsNames(dds_PE)

#Comparison 

PAH_vs_Control_PE_res <- results(dds_PE, contrast =c("condition", "PAH", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in PAH Relative to Control 
view(PAH_vs_Control_PE_res)

#Saving File

write.csv(PAH_vs_Control_PE_res, file = "MB250115 Pseudobulk PE DE PAH_vs_Control")
p1 <- as.data.frame(PAH_vs_Control_PE_res)
p2 <- subset(p1, subset = padj <= 0.05 &
               abs(log2FoldChange) >= 0.25) #padj is FDR
write.csv(p2, file = "PAH_vs_Control Plasma Erythrocyte Pathway Analysis.csv")

#Volcano Plot

EnhancedVolcano(PAH_vs_Control_PE_res,
                lab = rownames(PAH_vs_Control_PE_res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = "PAH vs Control Plasma Erythrocyte")

#Plasma/Erythrocyte DEGs (Whole Cluster)

counts_BC <- cts.split.modified$`B Cell`
colData <- data.frame(samples = colnames(counts_BC))

colData$condition <- c('Control','Control','Control','Control','PAH','PAH','PAH','PAH','PAH')

colData <- column_to_rownames(colData, var = 'samples')

dds_BC <- DESeqDataSetFromMatrix(countData = counts_BC,
                                 colData = colData,
                                 design = ~condition)

keep <- rowSums(counts(dds_BC)) >=10
dds_BC <- dds_BC[keep,]

dds_BC <- DESeq(dds_BC)
resultsNames(dds_BC)

#Comparison 

PAH_vs_Control_BC_res <- results(dds_BC, contrast =c("condition", "PAH", "Control"), cooksCutoff = FALSE, independentFiltering = FALSE) #Different in PAH Relative to Control 
view(PAH_vs_Control_BC_res)

#Saving File

write.csv(PAH_vs_Control_BC_res, file = "MB250115 Pseudobulk BC DE PAH_vs_Control")
b1 <- as.data.frame(PAH_vs_Control_BC_res)
b2 <- subset(b1, subset = padj <= 0.05 &
               abs(log2FoldChange) >= 0.25) #padj is FDR
write.csv(b2, file = "PAH_vs_Control B Cell Pathway Analysis.csv")

#Volcano Plot

EnhancedVolcano(PAH_vs_Control_BC_res,
                lab = rownames(PAH_vs_Control_BC_res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = "PAH vs Control B Cell")




#Combining all Veh vs Con DEGs into one graph - need to pseudobulk and run DESeq2 first

#Load in Files

PAH_vs_Control_HP_res <- read.csv(file='/users/1/blake561/MB250115 Pseudobulk HP DE PAH_vs_Control')
PAH_vs_Control_EC_res <- read.csv(file='/users/1/blake561/MB250115 Pseudobulk EC DE PAH_vs_Control')
PAH_vs_Control_FB_res <- read.csv(file='/users/1/blake561/MB250115 Pseudobulk FB DE PAH_vs_Control')
PAH_vs_Control_MP_res <- read.csv(file='/users/1/blake561/MB250123 Pseudobulk MP DE PAH_vs_Control')
PAH_vs_Control_CL_res <- read.csv(file='/users/1/blake561/MB250115 Pseudobulk CL DE PAH_vs_Control')
PAH_vs_Control_LC_res <- read.csv(file='/users/1/blake561/MB250115 Pseudobulk LC DE PAH_vs_Control')
PAH_vs_Control_BC_res <- read.csv(file='/users/1/blake561/MB250115 Pseudobulk BC DE PAH_vs_Control')
PAH_vs_Control_PE_res <- read.csv(file='/users/1/blake561/MB250115 Pseudobulk PE DE PAH_vs_Control')

#turn deseq objects into data frames
cv_hep <- as.data.frame(PAH_vs_Control_HP_res)
cv_endo <- as.data.frame(PAH_vs_Control_EC_res)
cv_hsc <- as.data.frame(PAH_vs_Control_FB_res)
cv_mac <- as.data.frame(PAH_vs_Control_MP_res)
cv_cholang <- as.data.frame(PAH_vs_Control_CL_res)
cv_lymph <- as.data.frame(PAH_vs_Control_LC_res)
cv_b <- as.data.frame(PAH_vs_Control_BC_res)
cv_plasma <- as.data.frame(PAH_vs_Control_PE_res)

#put numbers in front to set left to right number order you want
cv_celltypes_combined <- rbind(cbind(cv_hep, group = "1 Hepatocytes"), cbind(cv_endo, group = "2 Endothelial"), cbind(cv_hsc, group = "3 HSC/mFB"),
                               cbind(cv_mac, group = "4 Macrophage"), cbind(cv_cholang, group = "5 Cholangiocyte"), cbind(cv_lymph, group = "6 Lymphocyte"), 
                               cbind(cv_b, group = "7 B Cell"), cbind(cv_plasma, group = "8 Plasma"))

view(cv_celltypes_combined)

padj.cutoff <- 0.05
lfc.cutoff <- 0.25
threshold <- cv_celltypes_combined$padj < padj.cutoff & abs(cv_celltypes_combined$log2FoldChange) > lfc.cutoff
cv_celltypes_combined$threshold <- threshold


cv_celltypes_combined2 <- subset(cv_celltypes_combined, cv_celltypes_combined$threshold == "TRUE" | cv_celltypes_combined$threshold == "FALSE")
view(cv_celltypes_combined2)

cv_celltypes_combined3 <- subset(cv_celltypes_combined, cv_celltypes_combined$threshold == "TRUE")

ggplot(cv_celltypes_combined2) +
  aes(x=group, y=log2FoldChange, color = threshold, levels = c()) + geom_jitter(alpha = 0.4) +
  scale_color_manual(name = "threshold",
                     values = c("TRUE" = "red1", "FALSE" = "grey90")) +
  ggtitle("PAH vs Control no NAs") + theme_classic() + 
  ylim (-7,7) + guides(y=guide_axis(minor.ticks = TRUE)) + 
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 



ggplot(cv_celltypes_combined3) +
  aes(x=group, y=log2FoldChange, color = threshold, levels = c()) + geom_jitter(alpha = 0.3) +
  scale_color_manual(name = "threshold",
                     values = c("TRUE" = "red1")) +
  ggtitle("Veh vs control Sig only") + theme_classic() + 
  xlab("Cell Type") + 
  ylab("Log2FC") +
  theme(legend.position = "right",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


#changing colors so pos and neg are different
pos_threshold <- cv_celltypes_combined2$padj < padj.cutoff & cv_celltypes_combined2$log2FoldChange > lfc.cutoff
cv_celltypes_combined2$pos_threshold <- pos_threshold

neg_threshold <- cv_celltypes_combined2$padj < padj.cutoff & cv_celltypes_combined2$log2FoldChange < -(lfc.cutoff)
cv_celltypes_combined2$neg_threshold <- neg_threshold

level <- cv_celltypes_combined2 %>%
  mutate(ifelse(cv_celltypes_combined2$pos_threshold == "TRUE" & cv_celltypes_combined2$neg_threshold == "FALSE", "level 1",
                ifelse(cv_celltypes_combined2$pos_threshold == "FALSE" & cv_celltypes_combined2$neg_threshold == "FALSE", "level 2",
                       ifelse(cv_celltypes_combined2$pos_threshold == "FALSE" & cv_celltypes_combined2$neg_threshold == "TRUE", "level 3", NA))))

cv_celltypes_combined2$level <- level

#FINAL graph used for figure
ggplot(cv_celltypes_combined2) +
  aes(x=group, y=log2FoldChange, color = level$`ifelse(...)`, levels = c()) + geom_jitter(alpha = 0.4) +
  scale_color_manual(name = "level",
                     values = c("level 1" = "red1", "level 2" = "grey90", "level 3" = "mediumblue")) +
  ggtitle("Veh vs control no NAs") + theme_classic() + geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  ylim (-7,7) + guides(y=guide_axis(minor.ticks = TRUE)) + 
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 





ggplot(cv_celltypes_combined2) +
  aes(x=group, y=log2FoldChange, color = level$`ifelse(...)`, levels = c()) + geom_jitter(alpha = 0.4) +
  scale_color_manual(name = "level",
                     values = c("level 1" = "red1", "level 2" = "grey90", "level 3" = "mediumblue")) +
  ggtitle("Veh vs control no NAs") + theme_classic() + geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  ylim (-7,7) + guides(y=guide_axis(minor.ticks = TRUE)) + 
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 


# Ketogenesis Genes
genes_of_interest <- c("HMGCS2", "HMGCL", "ACAT1", "BDH1", "SCP2", "CPT1A", "OXCT1")

#Hepatocytes 

d1 <- as.data.frame(PAH_vs_Control_HP_res)
d2 <- subset(d1, rownames(d1) %in% genes_of_interest)

d2$gene <- rownames(d2)  # Add gene names as a column

# Create the bar plot with padj values as text labels
ggplot(d2, aes(x = gene, y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ifelse(padj < 0.05, sprintf("padj = %.2e", padj), "")), 
            vjust = -0.5, size = 3, color = "black") +  # Add padj values as labels
  scale_fill_manual(values = c("blue", "red"), labels = c("Down", "Up")) +
  labs(
    title = "PAH vs Control: Hepatocytes - Key Ketogenic Genes",
    x = "Gene",
    y = "Log2 Fold Change",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Endothelial Cells

e1 <- as.data.frame(PAH_vs_Control_EC_res)
e2 <- subset(e1, rownames(e1) %in% genes_of_interest)

e2$gene <- rownames(d2)  # Add gene names as a column

ggplot(e2, aes(x = gene, y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ifelse(padj < 0.05, sprintf("padj = %.2e", padj), "")), 
            vjust = -0.5, size = 3, color = "black") +  # Add padj values as labels
  scale_fill_manual(values = c("blue", "red"), labels = c("Down", "Up")) +
  labs(
    title = "PAH vs Control: Endothelial Cells - Key Ketogenic Genes",
    x = "Gene",
    y = "Log2 Fold Change",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#HSC/mFB

f1 <- as.data.frame(PAH_vs_Control_FB_res)
f2 <- subset(f1, rownames(1) %in% genes_of_interest)

f2$gene <- rownames(d2)  # Add gene names as a column

ggplot(f2, aes(x = gene, y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ifelse(padj < 0.05, sprintf("padj = %.2e", padj), "")), 
            vjust = -0.5, size = 3, color = "black") +  # Add padj values as labels
  scale_fill_manual(values = c("blue", "red"), labels = c("Down", "Up")) +
  labs(
    title = "PAH vs Control: HSC/mFB - Key Ketogenic Genes",
    x = "Gene",
    y = "Log2 Fold Change",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Vasoactive Genes
genes_of_interest <- c("BMP9", "BMP10", "ENG", "EDN1",'GDF2','BMPR2','FLT1','GDF15') #Vasoactive Genes

genes_of_interest <- c("MKI67","CCNB1",'AFP','JAG1','YAP1','TF','CYP3A4','CYP2C9','CYP1A2') 

genes_of_interest <- c("GCK","HSD17B13",'FASN','ME1','THRSP','CYP17A1','CYP2A7','CYP4F22') #Hepatocyte Maturation Markers


#Hepatocytes 

d1 <- as.data.frame(PAH_vs_Control_HP_res)
d2 <- subset(d1, rownames(d1) %in% genes_of_interest)

d2$gene <- rownames(d2)  # Add gene names as a column

# Create the bar plot with padj values as text labels
ggplot(d2, aes(x = gene, y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ifelse(padj < 0.05, sprintf("padj = %.2e", padj), "")), 
            vjust = -0.5, size = 3, color = "black") +  # Add padj values as labels
  scale_fill_manual(values = c("blue", "red"), labels = c("Down", "Up")) +
  labs(
    title = "PAH vs Control: Hepatocytes - Maturation Markers",
    x = "Gene",
    y = "Log2 Fold Change",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Endothelial Cells

e1 <- as.data.frame(PAH_vs_Control_EC_res)
e2 <- subset(e1, rownames(e1) %in% genes_of_interest)

e2$gene <- rownames(e2)  # Add gene names as a column

ggplot(e2, aes(x = gene, y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ifelse(padj < 0.05, sprintf("padj = %.2e", padj), "")), 
            vjust = -0.5, size = 3, color = "black") +  # Add padj values as labels
  scale_fill_manual(values = c("red", "blue"), labels = c("Up", 'Down')) +
  labs(
    title = "PAH vs Control: Endothelial Cells - Vasoactive Genes",
    x = "Gene",
    y = "Log2 Fold Change",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#HSC/mFB

f1 <- as.data.frame(PAH_vs_Control_FB_res)
f2 <- subset(f1, rownames(f1) %in% genes_of_interest)

f2$gene <- rownames(f2)  # Add gene names as a column

ggplot(f2, aes(x = gene, y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ifelse(padj < 0.05, sprintf("padj = %.2e", padj), "")), 
            vjust = -0.5, size = 3, color = "black") +  # Add padj values as labels
  scale_fill_manual(values = c("blue", "red"), labels = c("Down", "Up")) +
  labs(
    title = "PAH vs Control: HSC/mFB - Vasoactive Genes",
    x = "Gene",
    y = "Log2 Fold Change",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


