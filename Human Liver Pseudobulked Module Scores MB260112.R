#Load File

PAH <- readRDS(file = '/projects/standard/prin0088/shared/Madelyn Liver and Kidney Projects/PAH_total.mtx_final')

# ═══════════════════════════════════════════════════════════════════════════════
# PSEUDOBULK MODULE SCORE ANALYSIS - PAH LIVER DATASET
# Hepatocytes: Cytochrome P450 | HSC/mFB: PI3K-Akt & HIF-1
# ═══════════════════════════════════════════════════════════════════════════════

library(Seurat)
library(tidyverse)

# ═══════════════════════════════════════════════════════════════════════════════
# DEFINE SHARED GENE SETS
# ═══════════════════════════════════════════════════════════════════════════════

# ── JAK-STAT Signaling Gene Set ────────────────────────────────────────────────
JAK_STAT <- c(
  "AKT3", "SPRY3", "SPRY1", "SPRY2", "STAM2", "IRF9", "PIAS3", "IL24", "CISH",
  "IL22RA2", "SOCS4", "CNTF", "CNTFR", "CREBBP", "CSF2", "CSF2RA", "CSF2RB",
  "CSF3", "CSF3R", "CSH1", "CTF1", "IL23R", "SPRED1", "IFNLR1", "SPRED2",
  "EP300", "EPO", "EPOR", "AKT1", "AKT2", "CLCF1", "PIK3R5", "CBLC", "GH1",
  "GH2", "GHR", "IFNL2", "IFNL3", "IFNL1", "GRB2", "IL19", "SOCS7", "IFNE",
  "IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNA10",
  "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA21", "IFNAR1", "IFNAR2",
  "IFNB1", "IFNG", "IFNGR1", "IFNGR2", "IFNW1", "IL2", "IL2RA", "IL2RB",
  "IL2RG", "IL3", "IL3RA", "IL4", "IL4R", "IL5", "IL5RA", "IL6", "IL6R",
  "IL6ST", "IL7", "IL7R", "IL9", "IL9R", "IL10", "IL10RA", "IL10RB", "IL11",
  "IL11RA", "IL12A", "IL12B", "IL12RB1", "IL12RB2", "IL13", "IL13RA1",
  "IL13RA2", "IL15", "IL15RA", "JAK1", "JAK2", "JAK3", "LEP", "LEPR", "LIF",
  "LIFR", "MPL", "MYC", "OSM", "IL20", "IL21R", "IL22", "IL23A", "PIAS4",
  "PIK3CA", "PIK3CB", "PIM1", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2",
  "IL20RA", "IL20RB", "IL26", "PRL", "PRLR", "IFNK", "PTPN6", "PTPN11",
  "IL22RA1", "IL21", "CCND1", "BCL2L1", "CRLF2", "SOS1", "SOS2", "STAT1",
  "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6", "TPO", "TYK2",
  "STAM", "SPRY4", "PIK3R3", "TSLP", "PIAS1", "SOCS1", "CBL", "CBLB",
  "SOCS2", "CCND2", "CCND3", "SOCS3", "PIAS2", "OSMR", "SOCS5"
)

# ═══════════════════════════════════════════════════════════════════════════════
# 1) HEPATOCYTES - CYTOCHROME P450 OXIDATION
# ═══════════════════════════════════════════════════════════════════════════════

# ── Subset Hepatocytes ─────────────────────────────────────────────────────────
Hepatocytes <- subset(PAH, cell_identity == "Hepatocyte")

# ── Create Pseudobulk Object ───────────────────────────────────────────────────
pseudo_Hep <- AggregateExpression(Hepatocytes, return.seurat = TRUE, group.by = "sample")

# ── Assign condition ───────────────────────────────────────────────────────────
pseudo_Hep$condition <- ifelse(grepl("^C|Control", rownames(pseudo_Hep@meta.data), ignore.case = TRUE), 
                               "Control", "PAH")

# ── Cytochrome P450 Gene Set ───────────────────────────────────────────────────
CYP450 <- c(
  "CYP46A1", "CYP4F8", "CYP2U1", "CYP2R1", "CYP4F22", "CYB5A", "CYP1A1", "CYP1A2", "CYP1B1",
  "CYP2A6", "CYP2A7", "CYP3A7", "CYP2A13", "CYP2B6", "CYP2C19", "CYP2C8", "CYP2C9", "CYP2C18",
  "CYP2D6", "CYP2E1", "CYP2F1", "CYP2J2", "CYP3A4", "CYP3A5", "CYP4A11", "CYP4B1", "CYP7A1",
  "CYP8B1", "CYP11A1", "CYP11B1", "CYP11B2", "CYP17A1", "CYP19A1", "CYP24A1", "CYP26A1",
  "CYP27A1", "CYP27B1", "CYP51A1", "CYB5R3", "CYP4Z1", "CYP2G1P", "CYP4X1", "CYP4A22", "CYP4V2",
  "CYP2S1", "CYP27C1", "CYP26C1", "CYP4F3", "CYB5R4", "CYP39A1", "CYB5R2", "CYB5R1", "POR",
  "CYP2W1", "CYP26B1", "CYP20A1", "CYP4F11", "CYP3A43", "CYP4F12", "CYB5B", "CYP4F2", "CYP7B1"
)

CYP450 <- intersect(CYP450, rownames(pseudo_Hep))
cat("Hepatocyte CYP450 genes found:", length(CYP450), "\n")

# ── Add Module Score ───────────────────────────────────────────────────────────
pseudo_Hep <- AddModuleScore(pseudo_Hep, list(CYP450), name = "CYP450_", assay = "RNA")

# ── Build Data Frame ───────────────────────────────────────────────────────────
Hep_cyp450 <- pseudo_Hep@meta.data %>%
  mutate(sample = rownames(pseudo_Hep@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Hepatocyte",
    module = "CYP450_Oxidation",
    score = CYP450_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

# ── Statistical Test ───────────────────────────────────────────────────────────
t_hep <- t.test(score ~ condition, data = Hep_cyp450)
cat("\nHepatocyte CYP450 - p-value:", signif(t_hep$p.value, 3), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# HEPATOCYTES - ADD HIF-1, AEROBIC GLYCOLYSIS, & OXPHOS
# ═══════════════════════════════════════════════════════════════════════════════

# ── HIF-1 Signaling Gene Set ───────────────────────────────────────────────────
HIF1 <- c(
  "AKT3", "CDKN1A", "CDKN1B", "EGLN2", "EGLN3", "CREBBP", "CYBB", "LDHAL6A", "EDN1",
  "EGF", "EGFR", "EIF4E", "EIF4EBP1", "ENO1", "ENO2", "ENO3", "EP300", "EPAS1", "EPO",
  "ERBB2", "AKT1", "AKT2", "ALDOA", "ALDOB", "ALDOC", "FLT1", "MTOR", "EIF4E1B",
  "GAPDH", "ANGPT1", "ANGPT2", "MKNK2", "HIF1A", "HK1", "HK2", "HK3", "HMOX1", "IFNG",
  "IFNGR1", "IFNGR2", "IGF1", "IGF1R", "IL6", "IL6R", "INS", "INSR", "LDHA", "LDHB",
  "LDHC", "ARNT", "LTBR", "NFKB1", "NOS2", "NOS3", "NPPA", "SERPINE1", "ANGPT4",
  "PDHA1", "PDHA2", "PDHB", "PDK1", "PFKFB3", "PFKL", "PFKM", "PFKP", "PGK1", "PGK2",
  "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2", "PLCG1", "PLCG2", "EGLN1",
  "PRKCA", "PRKCB", "PRKCG", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "BCL2", "RELA",
  "RPS6", "RPS6KB1", "RPS6KB2", "SLC2A1", "STAT3", "ELOC", "ELOB", "TEK", "TF", "TFRC",
  "TIMP1", "TLR4", "VEGFA", "VHL", "HKDC1", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G",
  "CUL2", "PIK3R3", "MKNK1", "LDHAL6B", "EIF4E2", "RBX1"
)

# ── Aerobic Glycolysis Gene Set (user-provided) ────────────────────────────────
AeroGlycolysis <- c(
  "ENO1", "ALDOA", "GAPDH", "GPI", "HK1", "LDHA", 
  "PFKM", "PGAM2", "PGK1", "PKM", "SLC2A1", "TPI1"
)

# ── OXPHOS Gene Set (user-provided) ────────────────────────────────────────────
OxPhos <- c(
  "COX6CP3", "COX17", "TCIRG1", "ATP5PD", "ATP5MG", "UQCR11", "COX6B2", "NDUFA11", 
  "ATP6V1G3", "COX4I1", "COX5B", "COX6A1", "COX6A2", "COX6B1", "COX6C", "COX7A1", 
  "COX7A2", "COX7B", "COX7C", "COX8A", "COX10", "COX11", "COX15", "CYC1", 
  "ATP6V0E2", "COX7B2", "ATP6V0A2", "ATP6V0D2", "ATP6V1C2", "PPA2", "UQCRQ", 
  "UQCR10", "COX8C", "NDUFS7", "ATP5MC1P5", "UQCRHL", "MT-ATP6", "MT-ATP8", 
  "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", 
  "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6", "NDUFA1", "NDUFA2", "NDUFA3", 
  "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7", "NDUFA8", "NDUFA9", "NDUFA10", 
  "NDUFAB1", "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", 
  "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFC1", "NDUFC2", "NDUFS1", 
  "NDUFS2", "NDUFS3", "NDUFV1", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS8", 
  "NDUFV2", "NDUFV3", "ATP12A", "ATP4A", "ATP4B", "ATP5F1A", "ATP5F1B", 
  "ATP6V0A4", "ATP5F1C", "ATP5F1D", "ATP6V1D", "ATP5F1E", "ATP5PB", "ATP5MC1", 
  "ATP6V1H", "ATP5MC2", "ATP5MC3", "ATP5ME", "ATP5PF", "ATP6V1A", "ATP6V1B1", 
  "ATP6V1B2", "ATP6V0C", "ATP6V1C1", "ATP6V1E1", "ATP6V0B", "ATP6V1G2", 
  "ATP6V0A1", "ATP6AP1", "ATP5PO", "PPA1", "NDUFA4L2", "SDHA", "SDHB", "SDHC", 
  "SDHD", "LHPP", "UQCR10P1", "UQCRB", "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRH", 
  "COX4I2", "ATP6V0E1", "ATP6V1E2", "ATP6V0D1", "COX7A2L", "ATP6V1F", "COX5A", 
  "ATP6V1G1", "ATP5MF"
)

# ── Intersect with available genes ─────────────────────────────────────────────
HIF1_hep <- intersect(HIF1, rownames(pseudo_Hep))
AeroGlycolysis_hep <- intersect(AeroGlycolysis, rownames(pseudo_Hep))
OxPhos_hep <- intersect(OxPhos, rownames(pseudo_Hep))

cat("Hepatocyte HIF-1 genes found:", length(HIF1_hep), "\n")
cat("Hepatocyte Aerobic Glycolysis genes found:", length(AeroGlycolysis_hep), "\n")
cat("Hepatocyte OxPhos genes found:", length(OxPhos_hep), "\n")

# ── Add Module Scores ──────────────────────────────────────────────────────────
pseudo_Hep <- AddModuleScore(pseudo_Hep, list(HIF1_hep), name = "HIF1_", assay = "RNA")
pseudo_Hep <- AddModuleScore(pseudo_Hep, list(AeroGlycolysis_hep), name = "AeroGlycolysis_", assay = "RNA")
pseudo_Hep <- AddModuleScore(pseudo_Hep, list(OxPhos_hep), name = "OxPhos_", assay = "RNA")

# ── Build Data Frames ──────────────────────────────────────────────────────────
Hep_hif1 <- pseudo_Hep@meta.data %>%
  mutate(sample = rownames(pseudo_Hep@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Hepatocyte",
    module = "HIF1_Signaling",
    score = HIF1_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

Hep_aeroglycolysis <- pseudo_Hep@meta.data %>%
  mutate(sample = rownames(pseudo_Hep@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Hepatocyte",
    module = "Aerobic_Glycolysis",
    score = AeroGlycolysis_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

Hep_oxphos <- pseudo_Hep@meta.data %>%
  mutate(sample = rownames(pseudo_Hep@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Hepatocyte",
    module = "OxPhos",
    score = OxPhos_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

# ── Statistical Tests ──────────────────────────────────────────────────────────
t_hep_hif1 <- t.test(score ~ condition, data = Hep_hif1)
t_hep_aeroglycolysis <- t.test(score ~ condition, data = Hep_aeroglycolysis)
t_hep_oxphos <- t.test(score ~ condition, data = Hep_oxphos)

cat("\nHepatocyte HIF-1 - p-value:", signif(t_hep_hif1$p.value, 3), "\n")
cat("Hepatocyte Aerobic Glycolysis - p-value:", signif(t_hep_aeroglycolysis$p.value, 3), "\n")
cat("Hepatocyte OxPhos - p-value:", signif(t_hep_oxphos$p.value, 3), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# HEPATOCYTES - GLUCONEOGENESIS (add after OxPhos section)
# ═══════════════════════════════════════════════════════════════════════════════

# ── Gluconeogenesis Gene Set ───────────────────────────────────────────────────
Gluconeogenesis <- c(
  "PCK1",
  "ENO1", "ENO2", "ENO3", "ENO4",
  "PGAM1", "PGAM2",
  "PGK1", "PGK2",
  "GAPDH",
  "TPI1",
  "ALDOA", "ALDOB", "ALDOC",
  "FBP1", "FBP2"
)

Gluconeogenesis_hep <- intersect(Gluconeogenesis, rownames(pseudo_Hep))
cat("Hepatocyte Gluconeogenesis genes found:", length(Gluconeogenesis_hep), "\n")

# ── Add Module Score ───────────────────────────────────────────────────────────
pseudo_Hep <- AddModuleScore(pseudo_Hep, list(Gluconeogenesis_hep), name = "Gluconeogenesis_", assay = "RNA")

# ── Build Data Frame ───────────────────────────────────────────────────────────
Hep_gluconeo <- pseudo_Hep@meta.data %>%
  mutate(sample = rownames(pseudo_Hep@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Hepatocyte",
    module = "Gluconeogenesis",
    score = Gluconeogenesis_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

# ── Statistical Test ───────────────────────────────────────────────────────────
t_hep_gluconeo <- t.test(score ~ condition, data = Hep_gluconeo)
cat("Hepatocyte Gluconeogenesis - p-value:", signif(t_hep_gluconeo$p.value, 3), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# HEPATOCYTES - FATTY ACID METABOLISM
# ═══════════════════════════════════════════════════════════════════════════════

# ── Fatty Acid Metabolism Gene Set ─────────────────────────────────────────────
FAM <- c(
  "ACAA1", "ACAA2", "ACADL", "ACADM", "ACADS", "ACADSB", "ACADVL",
  "ACAT1", "ACAT2", "ACOX1", "ACOX3", "ACSL1", "ACSL3", "ACSL4",
  "ACSL5", "ACSL6", "ADH1A", "ADH1B", "ADH1C", "ADH4", "ADH5",
  "ADH6", "ADH7", "ADHFE1", "ALDH1A3", "ALDH1B1", "ALDH2",
  "ALDH3A1", "ALDH3A2", "ALDH7A1", "ALDH9A1", "CPT1A", "CPT1B",
  "CPT1C", "CPT2", "CYP4A11", "CYP4A22", "ECHS1", "EHHADH", "GCDH",
  "HADH", "HADHA", "HADHB", "HSD17B10", "HSD17B4"
)

FAM_hep <- intersect(FAM, rownames(pseudo_Hep))
cat("Hepatocyte Fatty Acid Metabolism genes found:", length(FAM_hep), "\n")

# ── Add Module Score ───────────────────────────────────────────────────────────
pseudo_Hep <- AddModuleScore(pseudo_Hep, list(FAM_hep), name = "FAM_", assay = "RNA")

# ── Build Data Frame ───────────────────────────────────────────────────────────
Hep_fam <- pseudo_Hep@meta.data %>%
  mutate(sample = rownames(pseudo_Hep@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Hepatocyte",
    module = "Fatty_Acid_Metabolism",
    score = FAM_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

# ── Statistical Test ───────────────────────────────────────────────────────────
t_hep_fam <- t.test(score ~ condition, data = Hep_fam)
cat("Hepatocyte Fatty Acid Metabolism - p-value:", signif(t_hep_fam$p.value, 3), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# HEPATOCYTES - JAK-STAT SIGNALING
# ═══════════════════════════════════════════════════════════════════════════════

# ── Intersect with available genes ─────────────────────────────────────────────
JAK_STAT_hep <- intersect(JAK_STAT, rownames(pseudo_Hep))
cat("Hepatocyte JAK-STAT genes found:", length(JAK_STAT_hep), "\n")

# ── Add Module Score ───────────────────────────────────────────────────────────
pseudo_Hep <- AddModuleScore(pseudo_Hep, list(JAK_STAT_hep), name = "JAK_STAT_", assay = "RNA")

# ── Build Data Frame ───────────────────────────────────────────────────────────
Hep_jakstat <- pseudo_Hep@meta.data %>%
  mutate(sample = rownames(pseudo_Hep@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Hepatocyte",
    module = "JAK_STAT",
    score = JAK_STAT_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

# ── Statistical Test ───────────────────────────────────────────────────────────
t_hep_jakstat <- t.test(score ~ condition, data = Hep_jakstat)
cat("Hepatocyte JAK-STAT - p-value:", signif(t_hep_jakstat$p.value, 3), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 2) HSC/mFB - PI3K-AKT & HIF-1 SIGNALING
# ═══════════════════════════════════════════════════════════════════════════════

# ── Subset HSC/mFB ─────────────────────────────────────────────────────────────
HSC <- subset(PAH, cell_identity == "HSC/mFB")

# ── Create Pseudobulk Object ───────────────────────────────────────────────────
pseudo_HSC <- AggregateExpression(HSC, return.seurat = TRUE, group.by = "sample")

# ── Assign condition ───────────────────────────────────────────────────────────
pseudo_HSC$condition <- ifelse(grepl("^C|Control", rownames(pseudo_HSC@meta.data), ignore.case = TRUE), 
                               "Control", "PAH")

# ── PI3K-Akt Gene Set ──────────────────────────────────────────────────────────
PI3K_AKT <- c(
  "AKT3", "BCL2L11", "SGK2", "LPAR6", "CDK2", "CDK4", "CDK6", 
  "CDKN1A", "CDKN1B", "LAMC3", "CREB3", "GNB5", "YWHAQ", "CHAD", 
  "CDC37", "CHRM1", "CHRM2", "CHUK", "THEM4", "PIK3AP1", 
  "COL1A1", "COL1A2", "COL2A1", "COL4A1", "COL4A2", "COL4A3", 
  "COL4A4", "COL4A5", "COL4A6", "COL6A1", "COL6A2", "COL6A3", "COL9A1", "COL9A2", 
  "COL9A3", "COMP", "COL6A6", "CREB1", "ATF2", "ATF6B", "CSF1", "CSF1R", "CSF3", 
  "CSF3R", "CSH1", "CSH2", "CSHL1", "PIK3R6", "CREB3L4", "LPAR1", "EFNA1", "EFNA2", 
  "EFNA3", "EFNA4", "EFNA5", "EGF", "EGFR", "EPHA2", "EIF4B", "EIF4E", "EIF4EBP1", 
  "CRTC2", "EPO", "EPOR", "ERBB2", "ERBB3", "ERBB4", "EREG", "AKT1", "AKT2", 
  "F2R", "FGF1", "FGF2", "FGF3", "FGF4", "FGF5", "FGF6", "FGF7", "FGF8", "FGF9", 
  "FGF10", "FGFR1", "FGFR3", "FGFR2", "FGFR4", "VEGFD", "LAMB4", "ITGA11", 
  "PHLPP2", "FOXO3", "FLT1", "FLT3", "FLT3LG", "PHLPP1", "FLT4", "FN1", "PIK3R5", 
  "LPAR3", "SGK3", "MTOR", "EIF4E1B", "G6PC1", "COL6A5", "FGF20", "FGF21", "GDNF", 
  "GH1", "GH2", "GHR", "FGF22", "GNB1", "GNB2", "GNB3", "GNG3", "GNG4", "GNG5", 
  "GNG7", "GNG10", "GNG11", "GNGT1", "GNGT2", "PPP2R3B", "ANGPT1", "LAMA1", 
  "LPAR4", "ANGPT2", "GRB2", "GSK3B", "PKN3", "GYS1", "GYS2", "HGF", "NR4A1", 
  "HRAS", "HSP90AA1", "HSP90AB1", "TNC", "IBSP", "IFNA1", "IFNA2", "IFNA4", 
  "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNA10", "IFNA13", "IFNA14", "IFNA16", 
  "IFNA17", "IFNA21", "IFNAR1", "IFNAR2", "IFNB1", "IGF1", "IGF1R", "IGF2", 
  "IKBKB", "IL2", "IL2RA", "FASLG", "IL2RB", "IL2RG", "IL3", "IL3RA", "IL4", 
  "IL4R", "IL6", "IL6R", "IL7", "IL7R", "INS", "INSR", "ITGA6", "IRS1", "ITGA1", 
  "ITGA2", "ITGA2B", "ITGA3", "ITGA4", "ITGA5", "ITGA7", "ITGA9", "ITGAV", 
  "ITGB1", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8", "JAK1", "JAK2", 
  "JAK3", "AREG", "KDR", "KIT", "KRAS", "LAMA2", "LAMA3", "LAMA4", "LAMA5", 
  "LAMB1", "LAMB2", "LAMB3", "LAMC1", "LAMC2", "MCL1", "MDM2", "MET", "KITLG", 
  "MTCP1", "MYB", "MYC", "ATF4", "NFKB1", "NGF", "NGFR", "NOS3", "NRAS", "NRTN", 
  "NTF3", "NTF4", "NTRK1", "NTRK2", "OSM", "PCK1", "PCK2", "ANGPT4", "PDGFA", 
  "PDGFB", "PDGFRA", "PDGFRB", "PDPK1", "GNG13", "PGF", "PIK3CA", "PIK3CB", 
  "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "GNG2", "DDIT4", "PPP2R3C", "PPP2CA", 
  "PPP2CB", "PPP2R1A", "PPP2R1B", "PPP2R2A", "PPP2R2B", "PPP2R2C", "PPP2R3A", 
  "PPP2R5A", "PPP2R5B", "PPP2R5C", "PPP2R5D", "PPP2R5E", "PRKAA1", "PRKAA2", 
  "PRKCA", "PPP2R2D", "PKN1", "PKN2", "MAPK1", "MAPK3", "GNG12", "PDGFC", 
  "MAP2K1", "MAP2K2", "PRL", "PRLR", "PSPN", "RELN", "LPAR5", "BAD", "PTEN", 
  "PTK2", "RPTOR", "G6PC2", "RAC1", "RAF1", "RBL2", "GNB4", "CCND1", "BCL2", 
  "RELA", "RET", "BCL2L1", "RHEB", "RPS6", "RPS6KB1", "RPS6KB2", "RXRA", "BDNF", 
  "TNN", "MLST8", "SGK1", "CREB3L2", "SOS1", "SOS2", "SPP1", "BRCA1", "STK11", 
  "SYK", "TEK", "TGFA", "THBS1", "THBS2", "THBS3", "THBS4", "TLR2", "TLR4", 
  "TNR", "TNXB", "TP53", "HSP90B1", "TSC1", "TSC2", "VEGFA", "VEGFB", "VEGFC", 
  "VTN", "VWF", "YWHAB", "YWHAE", "YWHAG", "YWHAH", "YWHAZ", "PDGFD", "FGF23", 
  "TCL1A", "CASP9", "CREB3L3", "PIK3R3", "ITGA10", "ITGA8", "IKBKG", "FGF18", 
  "FGF17", "FGF16", "CCND2", "CCND3", "CCNE1", "ARTN", "CREB3L1", "CCNE2", 
  "LPAR2", "OSMR", "MAGI1", "G6PC3", "CD19", "GNG8", "EIF4E2", "CREB5", "TCL1B", 
  "MAGI2", "FGF19"
)

PI3K_AKT <- intersect(PI3K_AKT, rownames(pseudo_HSC))
cat("HSC PI3K-Akt genes found:", length(PI3K_AKT), "\n")

# ── HIF-1 Signaling Gene Set ───────────────────────────────────────────────────
HIF1 <- c(
  "AKT3", "CDKN1A", "CDKN1B", "EGLN2", "EGLN3", "CREBBP", "CYBB", "LDHAL6A", "EDN1",
  "EGF", "EGFR", "EIF4E", "EIF4EBP1", "ENO1", "ENO2", "ENO3", "EP300", "EPAS1", "EPO",
  "ERBB2", "AKT1", "AKT2", "ALDOA", "ALDOB", "ALDOC", "FLT1", "MTOR", "EIF4E1B",
  "GAPDH", "ANGPT1", "ANGPT2", "MKNK2", "HIF1A", "HK1", "HK2", "HK3", "HMOX1", "IFNG",
  "IFNGR1", "IFNGR2", "IGF1", "IGF1R", "IL6", "IL6R", "INS", "INSR", "LDHA", "LDHB",
  "LDHC", "ARNT", "LTBR", "NFKB1", "NOS2", "NOS3", "NPPA", "SERPINE1", "ANGPT4",
  "PDHA1", "PDHA2", "PDHB", "PDK1", "PFKFB3", "PFKL", "PFKM", "PFKP", "PGK1", "PGK2",
  "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2", "PLCG1", "PLCG2", "EGLN1",
  "PRKCA", "PRKCB", "PRKCG", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "BCL2", "RELA",
  "RPS6", "RPS6KB1", "RPS6KB2", "SLC2A1", "STAT3", "ELOC", "ELOB", "TEK", "TF", "TFRC",
  "TIMP1", "TLR4", "VEGFA", "VHL", "HKDC1", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G",
  "CUL2", "PIK3R3", "MKNK1", "LDHAL6B", "EIF4E2", "RBX1"
)

HIF1 <- intersect(HIF1, rownames(pseudo_HSC))
cat("HSC HIF-1 genes found:", length(HIF1), "\n")

# ── Add Module Scores ──────────────────────────────────────────────────────────
pseudo_HSC <- AddModuleScore(pseudo_HSC, list(PI3K_AKT), name = "PI3K_AKT_", assay = "RNA")
pseudo_HSC <- AddModuleScore(pseudo_HSC, list(HIF1), name = "HIF1_", assay = "RNA")

# ── Build Data Frames ──────────────────────────────────────────────────────────
HSC_pi3k <- pseudo_HSC@meta.data %>%
  mutate(sample = rownames(pseudo_HSC@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "HSC_mFB",
    module = "PI3K_Akt",
    score = PI3K_AKT_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

HSC_hif1 <- pseudo_HSC@meta.data %>%
  mutate(sample = rownames(pseudo_HSC@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "HSC_mFB",
    module = "HIF1_Signaling",
    score = HIF1_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

# ── Statistical Tests ──────────────────────────────────────────────────────────
t_hsc_pi3k <- t.test(score ~ condition, data = HSC_pi3k)
t_hsc_hif1 <- t.test(score ~ condition, data = HSC_hif1)
cat("HSC PI3K-Akt - p-value:", signif(t_hsc_pi3k$p.value, 3), "\n")
cat("HSC HIF-1 - p-value:", signif(t_hsc_hif1$p.value, 3), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# HSC/mFB - JAK-STAT SIGNALING
# ═══════════════════════════════════════════════════════════════════════════════

# ── Intersect with available genes ─────────────────────────────────────────────
JAK_STAT_hsc <- intersect(JAK_STAT, rownames(pseudo_HSC))
cat("HSC JAK-STAT genes found:", length(JAK_STAT_hsc), "\n")

# ── Add Module Score ───────────────────────────────────────────────────────────
pseudo_HSC <- AddModuleScore(pseudo_HSC, list(JAK_STAT_hsc), name = "JAK_STAT_", assay = "RNA")

# ── Build Data Frame ───────────────────────────────────────────────────────────
HSC_jakstat <- pseudo_HSC@meta.data %>%
  mutate(sample = rownames(pseudo_HSC@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "HSC_mFB",
    module = "JAK_STAT",
    score = JAK_STAT_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

# ── Statistical Test ───────────────────────────────────────────────────────────
t_hsc_jakstat <- t.test(score ~ condition, data = HSC_jakstat)
cat("HSC JAK-STAT - p-value:", signif(t_hsc_jakstat$p.value, 3), "\n")


# ═══════════════════════════════════════════════════════════════════════════════
# 3) COMBINE RESULTS & EXPORT
# ═══════════════════════════════════════════════════════════════════════════════

# ── Combine individual sample scores ───────────────────────────────────────────
combined_scores <- bind_rows(Hep_cyp450, HSC_pi3k, HSC_hif1) %>%
  arrange(cell_type, module, condition, sample)

# ── Create summary statistics ──────────────────────────────────────────────────
summary_stats <- combined_scores %>%
  group_by(cell_type, module, condition) %>%
  summarise(
    n = n(),
    mean = mean(score),
    median = median(score),
    sd = sd(score),
    se = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(p_value = case_when(
    module == "CYP450_Oxidation" ~ t_hep$p.value,
    module == "PI3K_Akt" ~ t_hsc_pi3k$p.value,
    module == "HIF1_Signaling" ~ t_hsc_hif1$p.value
  ))

# ── Export ─────────────────────────────────────────────────────────────────────
write.csv(combined_scores, "PAH_Liver_Pseudobulk_Module_Scores_by_Sample.csv", row.names = FALSE)
write.csv(summary_stats, "PAH_Liver_Pseudobulk_Module_Summary.csv", row.names = FALSE)

cat("\n✓ Exported: PAH_Liver_Pseudobulk_Module_Scores_by_Sample.csv")
cat("\n✓ Exported: PAH_Liver_Pseudobulk_Module_Summary.csv\n")

print(combined_scores)
print(summary_stats)

# ═══════════════════════════════════════════════════════════════════════════════
# PSEUDOBULK SINGLE GENE EXPRESSION - IL6 in HSC/mFB
# ═══════════════════════════════════════════════════════════════════════════════

# ── Extract IL6 Expression ─────────────────────────────────────────────────────
# Get normalized expression from the RNA assay
IL6_expr <- FetchData(pseudo_HSC, vars = "IL6", layer = "data")

# Build data frame
HSC_IL6 <- data.frame(
  sample = rownames(pseudo_HSC@meta.data),
  condition = factor(pseudo_HSC$condition, levels = c("Control", "PAH")),
  IL6_expression = IL6_expr$IL6
)

print(HSC_IL6)

# ── Statistical Test ───────────────────────────────────────────────────────────
t_IL6 <- t.test(IL6_expression ~ condition, data = HSC_IL6)
cat("\nHSC IL6 expression - p-value:", signif(t_IL6$p.value, 3), "\n")

# ── Summary Stats ──────────────────────────────────────────────────────────────
IL6_summary <- HSC_IL6 %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean = mean(IL6_expression),
    median = median(IL6_expression),
    sd = sd(IL6_expression),
    se = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(p_value = t_IL6$p.value)

print(IL6_summary)

# ── Visualization ──────────────────────────────────────────────────────────────
ggplot(HSC_IL6, aes(x = condition, y = IL6_expression, fill = condition)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, colour = "black") +
  geom_errorbar(data = IL6_summary,
                aes(x = condition, y = mean, ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(size = 3, shape = 21, fill = "white", colour = "black") +
  annotate("text", x = 1.5,
           y = max(HSC_IL6$IL6_expression) * 1.1,
           label = paste0("p = ", signif(t_IL6$p.value, 3)),
           size = 4) +
  scale_fill_manual(values = c("Control" = "salmon", "PAH" = "steelblue")) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  labs(
    y = "IL6 Expression (Pseudobulk)",
    title = "HSC/mFB: IL6 Expression (PAH vs Control)"
  )

# ── Export ─────────────────────────────────────────────────────────────────────
write.csv(HSC_IL6, "PAH_HSC_IL6_Pseudobulk_Expression.csv", row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# PSEUDOBULK SINGLE GENE EXPRESSION - HIF1A in Hepatocytes
# ═══════════════════════════════════════════════════════════════════════════════

# ── Extract HIF1A Expression ───────────────────────────────────────────────────
# Get normalized expression from the RNA assay
HIF1A_expr <- FetchData(pseudo_Hep, vars = "HIF1A", layer = "data")

# Build data frame
Hep_HIF1A <- data.frame(
  sample = rownames(pseudo_Hep@meta.data),
  condition = factor(pseudo_Hep$condition, levels = c("Control", "PAH")),
  HIF1A_expression = HIF1A_expr$HIF1A
)

print(Hep_HIF1A)

# ── Statistical Test ───────────────────────────────────────────────────────────
t_HIF1A <- t.test(HIF1A_expression ~ condition, data = Hep_HIF1A)
cat("\nHepatocyte HIF1A expression - p-value:", signif(t_HIF1A$p.value, 3), "\n")

# ── Summary Stats ──────────────────────────────────────────────────────────────
HIF1A_summary <- Hep_HIF1A %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean = mean(HIF1A_expression),
    median = median(HIF1A_expression),
    sd = sd(HIF1A_expression),
    se = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(p_value = t_HIF1A$p.value)

print(HIF1A_summary)

# ── Visualization ──────────────────────────────────────────────────────────────
ggplot(Hep_HIF1A, aes(x = condition, y = HIF1A_expression, fill = condition)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, colour = "black") +
  geom_errorbar(data = HIF1A_summary,
                aes(x = condition, y = mean, ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(size = 3, shape = 21, fill = "white", colour = "black") +
  annotate("text", x = 1.5,
           y = max(Hep_HIF1A$HIF1A_expression) * 1.1,
           label = paste0("p = ", signif(t_HIF1A$p.value, 3)),
           size = 4) +
  scale_fill_manual(values = c("Control" = "salmon", "PAH" = "steelblue")) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  labs(
    y = "HIF1A Expression (Pseudobulk)",
    title = "Hepatocyte: HIF1A Expression (PAH vs Control)"
  )

# ── Export ─────────────────────────────────────────────────────────────────────
write.csv(Hep_HIF1A, "PAH_Hepatocyte_HIF1A_Pseudobulk_Expression.csv", row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# HSC/mFB - ECM-RECEPTOR INTERACTION (add to existing HSC analysis)
# ═══════════════════════════════════════════════════════════════════════════════

# ── ECM-Receptor Interaction Gene Set (KEGG) ───────────────────────────────────
ECM_RECEPTOR <- c(
  "AGRN", "CD36", "CD44", "CD47", "COL1A1", "COL1A2", "COL2A1", "COL3A1", 
  "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6", "COL5A1", 
  "COL5A2", "COL5A3", "COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6",
  "COL9A1", "COL9A2", "COL9A3", "COL11A1", "COL11A2", "COMP", "DAG1",
  "DMD", "DTNA", "DTNB", "FN1", "FNDC1", "GP1BA", "GP1BB", "GP5", "GP6", 
  "GP9", "HMMR", "HSPG2", "IBSP", "ITGA1", "ITGA2", "ITGA2B", "ITGA3", 
  "ITGA4", "ITGA5", "ITGA6", "ITGA7", "ITGA8", "ITGA9", "ITGA10", "ITGA11", 
  "ITGAV", "ITGB1", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8",
  "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMB3", 
  "LAMB4", "LAMC1", "LAMC2", "LAMC3", "NPNT", "RELN", "SDC1", "SDC2", 
  "SDC3", "SDC4", "SGCA", "SGCB", "SGCD", "SGCE", "SGCG", "SGCZ", "SNTA1", 
  "SNTB1", "SNTB2", "SNTG1", "SNTG2", "SPP1", "SSPN", "SV2A", "SV2B", 
  "SV2C", "THBS1", "THBS2", "THBS3", "THBS4", "TNC", "TNN", "TNR", "TNXB", 
  "UTRN", "VTN", "VWF"
)

ECM_RECEPTOR <- intersect(ECM_RECEPTOR, rownames(pseudo_HSC))
cat("HSC ECM-Receptor genes found:", length(ECM_RECEPTOR), "\n")

# ── Add Module Score ───────────────────────────────────────────────────────────
pseudo_HSC <- AddModuleScore(pseudo_HSC, list(ECM_RECEPTOR), name = "ECM_Receptor_", assay = "RNA")

# ── Build Data Frame ───────────────────────────────────────────────────────────
HSC_ecm <- pseudo_HSC@meta.data %>%
  mutate(sample = rownames(pseudo_HSC@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "HSC_mFB",
    module = "ECM_Receptor",
    score = ECM_Receptor_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

# ── Statistical Test ───────────────────────────────────────────────────────────
t_hsc_ecm <- t.test(score ~ condition, data = HSC_ecm)
cat("HSC ECM-Receptor - p-value:", signif(t_hsc_ecm$p.value, 3), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# ENDOTHELIAL CELLS - HIF-1, PI3K-AKT, OXPHOS, GLYCOLYSIS, & PPP
# ═══════════════════════════════════════════════════════════════════════════════

# ── Subset Endothelial Cells ───────────────────────────────────────────────────
Endothelial <- subset(PAH, cell_identity == "Endothelial")

# ── Create Pseudobulk Object ───────────────────────────────────────────────────
pseudo_Endo <- AggregateExpression(Endothelial, return.seurat = TRUE, group.by = "sample")

# ── Assign condition ───────────────────────────────────────────────────────────
pseudo_Endo$condition <- ifelse(grepl("^C|Control", rownames(pseudo_Endo@meta.data), ignore.case = TRUE), 
                                "Control", "PAH")

# ── HIF-1 Signaling Gene Set ───────────────────────────────────────────────────
HIF1 <- c(
  "AKT3", "CDKN1A", "CDKN1B", "EGLN2", "EGLN3", "CREBBP", "CYBB", "LDHAL6A", "EDN1",
  "EGF", "EGFR", "EIF4E", "EIF4EBP1", "ENO1", "ENO2", "ENO3", "EP300", "EPAS1", "EPO",
  "ERBB2", "AKT1", "AKT2", "ALDOA", "ALDOB", "ALDOC", "FLT1", "MTOR", "EIF4E1B",
  "GAPDH", "ANGPT1", "ANGPT2", "MKNK2", "HIF1A", "HK1", "HK2", "HK3", "HMOX1", "IFNG",
  "IFNGR1", "IFNGR2", "IGF1", "IGF1R", "IL6", "IL6R", "INS", "INSR", "LDHA", "LDHB",
  "LDHC", "ARNT", "LTBR", "NFKB1", "NOS2", "NOS3", "NPPA", "SERPINE1", "ANGPT4",
  "PDHA1", "PDHA2", "PDHB", "PDK1", "PFKFB3", "PFKL", "PFKM", "PFKP", "PGK1", "PGK2",
  "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2", "PLCG1", "PLCG2", "EGLN1",
  "PRKCA", "PRKCB", "PRKCG", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "BCL2", "RELA",
  "RPS6", "RPS6KB1", "RPS6KB2", "SLC2A1", "STAT3", "ELOC", "ELOB", "TEK", "TF", "TFRC",
  "TIMP1", "TLR4", "VEGFA", "VHL", "HKDC1", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G",
  "CUL2", "PIK3R3", "MKNK1", "LDHAL6B", "EIF4E2", "RBX1"
)

# ── PI3K-Akt Signaling Gene Set ────────────────────────────────────────────────
PI3K_AKT <- c(
  "AKT3", "BCL2L11", "SGK2", "LPAR6", "CDK2", "CDK4", "CDK6", 
  "CDKN1A", "CDKN1B", "LAMC3", "CREB3", "GNB5", "YWHAQ", "CHAD", 
  "CDC37", "CHRM1", "CHRM2", "CHUK", "THEM4", "PIK3AP1", 
  "COL1A1", "COL1A2", "COL2A1", "COL4A1", "COL4A2", "COL4A3", 
  "COL4A4", "COL4A5", "COL4A6", "COL6A1", "COL6A2", "COL6A3", "COL9A1", "COL9A2", 
  "COL9A3", "COMP", "COL6A6", "CREB1", "ATF2", "ATF6B", "CSF1", "CSF1R", "CSF3", 
  "CSF3R", "CSH1", "CSH2", "CSHL1", "PIK3R6", "CREB3L4", "LPAR1", "EFNA1", "EFNA2", 
  "EFNA3", "EFNA4", "EFNA5", "EGF", "EGFR", "EPHA2", "EIF4B", "EIF4E", "EIF4EBP1", 
  "CRTC2", "EPO", "EPOR", "ERBB2", "ERBB3", "ERBB4", "EREG", "AKT1", "AKT2", 
  "F2R", "FGF1", "FGF2", "FGF3", "FGF4", "FGF5", "FGF6", "FGF7", "FGF8", "FGF9", 
  "FGF10", "FGFR1", "FGFR3", "FGFR2", "FGFR4", "VEGFD", "LAMB4", "ITGA11", 
  "PHLPP2", "FOXO3", "FLT1", "FLT3", "FLT3LG", "PHLPP1", "FLT4", "FN1", "PIK3R5", 
  "LPAR3", "SGK3", "MTOR", "EIF4E1B", "G6PC1", "COL6A5", "FGF20", "FGF21", "GDNF", 
  "GH1", "GH2", "GHR", "FGF22", "GNB1", "GNB2", "GNB3", "GNG3", "GNG4", "GNG5", 
  "GNG7", "GNG10", "GNG11", "GNGT1", "GNGT2", "PPP2R3B", "ANGPT1", "LAMA1", 
  "LPAR4", "ANGPT2", "GRB2", "GSK3B", "PKN3", "GYS1", "GYS2", "HGF", "NR4A1", 
  "HRAS", "HSP90AA1", "HSP90AB1", "TNC", "IBSP", "IFNA1", "IFNA2", "IFNA4", 
  "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNA10", "IFNA13", "IFNA14", "IFNA16", 
  "IFNA17", "IFNA21", "IFNAR1", "IFNAR2", "IFNB1", "IGF1", "IGF1R", "IGF2", 
  "IKBKB", "IL2", "IL2RA", "FASLG", "IL2RB", "IL2RG", "IL3", "IL3RA", "IL4", 
  "IL4R", "IL6", "IL6R", "IL7", "IL7R", "INS", "INSR", "ITGA6", "IRS1", "ITGA1", 
  "ITGA2", "ITGA2B", "ITGA3", "ITGA4", "ITGA5", "ITGA7", "ITGA9", "ITGAV", 
  "ITGB1", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8", "JAK1", "JAK2", 
  "JAK3", "AREG", "KDR", "KIT", "KRAS", "LAMA2", "LAMA3", "LAMA4", "LAMA5", 
  "LAMB1", "LAMB2", "LAMB3", "LAMC1", "LAMC2", "MCL1", "MDM2", "MET", "KITLG", 
  "MTCP1", "MYB", "MYC", "ATF4", "NFKB1", "NGF", "NGFR", "NOS3", "NRAS", "NRTN", 
  "NTF3", "NTF4", "NTRK1", "NTRK2", "OSM", "PCK1", "PCK2", "ANGPT4", "PDGFA", 
  "PDGFB", "PDGFRA", "PDGFRB", "PDPK1", "GNG13", "PGF", "PIK3CA", "PIK3CB", 
  "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "GNG2", "DDIT4", "PPP2R3C", "PPP2CA", 
  "PPP2CB", "PPP2R1A", "PPP2R1B", "PPP2R2A", "PPP2R2B", "PPP2R2C", "PPP2R3A", 
  "PPP2R5A", "PPP2R5B", "PPP2R5C", "PPP2R5D", "PPP2R5E", "PRKAA1", "PRKAA2", 
  "PRKCA", "PPP2R2D", "PKN1", "PKN2", "MAPK1", "MAPK3", "GNG12", "PDGFC", 
  "MAP2K1", "MAP2K2", "PRL", "PRLR", "PSPN", "RELN", "LPAR5", "BAD", "PTEN", 
  "PTK2", "RPTOR", "G6PC2", "RAC1", "RAF1", "RBL2", "GNB4", "CCND1", "BCL2", 
  "RELA", "RET", "BCL2L1", "RHEB", "RPS6", "RPS6KB1", "RPS6KB2", "RXRA", "BDNF", 
  "TNN", "MLST8", "SGK1", "CREB3L2", "SOS1", "SOS2", "SPP1", "BRCA1", "STK11", 
  "SYK", "TEK", "TGFA", "THBS1", "THBS2", "THBS3", "THBS4", "TLR2", "TLR4", 
  "TNR", "TNXB", "TP53", "HSP90B1", "TSC1", "TSC2", "VEGFA", "VEGFB", "VEGFC", 
  "VTN", "VWF", "YWHAB", "YWHAE", "YWHAG", "YWHAH", "YWHAZ", "PDGFD", "FGF23", 
  "TCL1A", "CASP9", "CREB3L3", "PIK3R3", "ITGA10", "ITGA8", "IKBKG", "FGF18", 
  "FGF17", "FGF16", "CCND2", "CCND3", "CCNE1", "ARTN", "CREB3L1", "CCNE2", 
  "LPAR2", "OSMR", "MAGI1", "G6PC3", "CD19", "GNG8", "EIF4E2", "CREB5", "TCL1B", 
  "MAGI2", "FGF19"
)

# ── OXPHOS Gene Set (KEGG Oxidative Phosphorylation) ───────────────────────────
OXPHOS <- c(
  # Complex I (NADH dehydrogenase)
  "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7", "NDUFA8",
  "NDUFA9", "NDUFA10", "NDUFA11", "NDUFA12", "NDUFA13", "NDUFAB1", "NDUFB1", 
  "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9",
  "NDUFB10", "NDUFB11", "NDUFC1", "NDUFC2", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4",
  "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2", "NDUFV3",
  # Complex II (Succinate dehydrogenase)
  "SDHA", "SDHB", "SDHC", "SDHD",
  # Complex III (Cytochrome bc1)
  "UQCR10", "UQCR11", "UQCRB", "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRH", "UQCRQ",
  "CYC1", "CYCS",
  # Complex IV (Cytochrome c oxidase)
  "COX4I1", "COX4I2", "COX5A", "COX5B", "COX6A1", "COX6A2", "COX6B1", "COX6B2",
  "COX6C", "COX7A1", "COX7A2", "COX7A2L", "COX7B", "COX7B2", "COX7C", "COX8A",
  "COX8C", "MT-CO1", "MT-CO2", "MT-CO3",
  # Complex V (ATP synthase)
  "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E", "ATP5MC1", "ATP5MC2",
  "ATP5MC3", "ATP5ME", "ATP5MF", "ATP5MG", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO",
  "MT-ATP6", "MT-ATP8",
  # Alternative names for some ATP synthase subunits
  "ATP5A1", "ATP5B", "ATP5C1", "ATP5D", "ATP5E", "ATP5G1", "ATP5G2", "ATP5G3",
  "ATP5H", "ATP5I", "ATP5J", "ATP5J2", "ATP5L", "ATP5O"
)

# ── Aerobic Glycolysis Gene Set ────────────────────────────────────────────────
AeroGlycolysis <- c(
  "ENO1", "ALDOA", "GAPDH", "GPI", "HK1", "LDHA", 
  "PFKM", "PGAM2", "PGK1", "PKM", "SLC2A1", "TPI1"
)

# ── Pentose Phosphate Pathway Gene Set ─────────────────────────────────────────
PPP <- c(
  "FBP1", "PRPS1L1", "ALDOA", "ALDOB", "RPIA", "ALDOC", "G6PD", "PGLS", "GPI", "DERA",
  "PFKL", "PFKM", "PFKP", "PGD", "PGM1", "PGM2", "PRPS1", "PRPS2", "RPE", "RBKS",
  "TALDO1", "TKT", "RPEL1", "TKTL1", "TKTL2", "FBP2", "H6PD"
)

# ── Intersect with available genes ─────────────────────────────────────────────
HIF1_endo <- intersect(HIF1, rownames(pseudo_Endo))
PI3K_AKT_endo <- intersect(PI3K_AKT, rownames(pseudo_Endo))
OXPHOS_endo <- intersect(OXPHOS, rownames(pseudo_Endo))
AeroGlycolysis_endo <- intersect(AeroGlycolysis, rownames(pseudo_Endo))
PPP_endo <- intersect(PPP, rownames(pseudo_Endo))

cat("Endothelial HIF-1 genes found:", length(HIF1_endo), "\n")
cat("Endothelial PI3K-Akt genes found:", length(PI3K_AKT_endo), "\n")
cat("Endothelial OXPHOS genes found:", length(OXPHOS_endo), "\n")
cat("Endothelial Aerobic Glycolysis genes found:", length(AeroGlycolysis_endo), "\n")
cat("Endothelial PPP genes found:", length(PPP_endo), "\n")

# ── Add Module Scores ──────────────────────────────────────────────────────────
pseudo_Endo <- AddModuleScore(pseudo_Endo, list(HIF1_endo), name = "HIF1_", assay = "RNA")
pseudo_Endo <- AddModuleScore(pseudo_Endo, list(PI3K_AKT_endo), name = "PI3K_AKT_", assay = "RNA")
pseudo_Endo <- AddModuleScore(pseudo_Endo, list(OXPHOS_endo), name = "OXPHOS_", assay = "RNA")
pseudo_Endo <- AddModuleScore(pseudo_Endo, list(AeroGlycolysis_endo), name = "AeroGlycolysis_", assay = "RNA")
pseudo_Endo <- AddModuleScore(pseudo_Endo, list(PPP_endo), name = "PPP_", assay = "RNA")

# ── Build Data Frames ──────────────────────────────────────────────────────────
Endo_hif1 <- pseudo_Endo@meta.data %>%
  mutate(sample = rownames(pseudo_Endo@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Endothelial",
    module = "HIF1_Signaling",
    score = HIF1_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

Endo_pi3k <- pseudo_Endo@meta.data %>%
  mutate(sample = rownames(pseudo_Endo@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Endothelial",
    module = "PI3K_Akt",
    score = PI3K_AKT_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

Endo_oxphos <- pseudo_Endo@meta.data %>%
  mutate(sample = rownames(pseudo_Endo@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Endothelial",
    module = "OXPHOS",
    score = OXPHOS_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

Endo_aeroglycolysis <- pseudo_Endo@meta.data %>%
  mutate(sample = rownames(pseudo_Endo@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Endothelial",
    module = "Aerobic_Glycolysis",
    score = AeroGlycolysis_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

Endo_ppp <- pseudo_Endo@meta.data %>%
  mutate(sample = rownames(pseudo_Endo@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Endothelial",
    module = "PPP",
    score = PPP_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

# ── Statistical Tests ──────────────────────────────────────────────────────────
t_endo_hif1 <- t.test(score ~ condition, data = Endo_hif1)
t_endo_pi3k <- t.test(score ~ condition, data = Endo_pi3k)
t_endo_oxphos <- t.test(score ~ condition, data = Endo_oxphos)
t_endo_aeroglycolysis <- t.test(score ~ condition, data = Endo_aeroglycolysis)
t_endo_ppp <- t.test(score ~ condition, data = Endo_ppp)

cat("\nEndothelial HIF-1 - p-value:", signif(t_endo_hif1$p.value, 3), "\n")
cat("Endothelial PI3K-Akt - p-value:", signif(t_endo_pi3k$p.value, 3), "\n")
cat("Endothelial OXPHOS - p-value:", signif(t_endo_oxphos$p.value, 3), "\n")
cat("Endothelial Aerobic Glycolysis - p-value:", signif(t_endo_aeroglycolysis$p.value, 3), "\n")
cat("Endothelial PPP - p-value:", signif(t_endo_ppp$p.value, 3), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# ENDOTHELIAL CELLS - JAK-STAT SIGNALING
# ═══════════════════════════════════════════════════════════════════════════════

# ── Intersect with available genes ─────────────────────────────────────────────
JAK_STAT_endo <- intersect(JAK_STAT, rownames(pseudo_Endo))
cat("Endothelial JAK-STAT genes found:", length(JAK_STAT_endo), "\n")

# ── Add Module Score ───────────────────────────────────────────────────────────
pseudo_Endo <- AddModuleScore(pseudo_Endo, list(JAK_STAT_endo), name = "JAK_STAT_", assay = "RNA")

# ── Build Data Frame ───────────────────────────────────────────────────────────
Endo_jakstat <- pseudo_Endo@meta.data %>%
  mutate(sample = rownames(pseudo_Endo@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Endothelial",
    module = "JAK_STAT",
    score = JAK_STAT_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

# ── Statistical Test ───────────────────────────────────────────────────────────
t_endo_jakstat <- t.test(score ~ condition, data = Endo_jakstat)
cat("Endothelial JAK-STAT - p-value:", signif(t_endo_jakstat$p.value, 3), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# MACROPHAGES - COMPLEMENT & COAGULATION CASCADE + HIF-1 SIGNALING
# ═══════════════════════════════════════════════════════════════════════════════

# ── Subset Macrophages ─────────────────────────────────────────────────────────
Macrophage <- subset(PAH, cell_identity == "Macrophage")

# ── Create Pseudobulk Object ───────────────────────────────────────────────────
pseudo_Mac <- AggregateExpression(Macrophage, return.seurat = TRUE, group.by = "sample")

# ── Assign condition ───────────────────────────────────────────────────────────
pseudo_Mac$condition <- ifelse(grepl("^C|Control", rownames(pseudo_Mac@meta.data), ignore.case = TRUE), 
                               "Control", "PAH")

# ── Complement & Coagulation Cascade Gene Set ──────────────────────────────────
Complement <- c(
  "C4B", "PROCR", "MASP2", "CFHR4", "CFHR3", "VSIG4", "CLU", "CPB2",
  "CR1", "CR1L", "CR2", "CD55", "CFD", "A2M", "F2", "F2R", "F2RL2",
  "F3", "F5", "F7", "F8", "F9", "F10", "F11", "F12", "F13A1", "F13B",
  "FGA", "FGB", "FGG", "SERPIND1", "CFH", "CFHR1", "CFHR2", "CFI",
  "ITGAM", "ITGAX", "ITGB2", "KLKB1", "KNG1", "MBL2", "CD46",
  "SERPINC1", "SERPINE1", "SERPINB2", "SERPINA5", "SERPINA1", "SERPINE2",
  "PLAT", "PLAU", "PLAUR", "PLG", "SERPINF2", "PROC", "PROS1", "MASP1",
  "BDKRB1", "BDKRB2", "CFB", "TFPI", "THBD", "SERPING1",
  "C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C3", "C3AR1",
  "C4A", "C4B", "C4BPA", "C4BPB", "C5", "C5AR1", "C6", "C7",
  "C8A", "C8B", "C8G", "C9", "VTN", "VWF", "CFHR5", "F2RL3", "CD59"
)

Complement_mac <- intersect(Complement, rownames(pseudo_Mac))
cat("Macrophage Complement genes found:", length(Complement_mac), "\n")

# ── Add Module Score ───────────────────────────────────────────────────────────
pseudo_Mac <- AddModuleScore(pseudo_Mac, list(Complement_mac), name = "Complement_", assay = "RNA")

# ── Build Data Frame ───────────────────────────────────────────────────────────
Mac_complement <- pseudo_Mac@meta.data %>%
  mutate(sample = rownames(pseudo_Mac@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Macrophage",
    module = "Complement_Coagulation",
    score = Complement_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

# ── Statistical Test ───────────────────────────────────────────────────────────
t_mac_complement <- t.test(score ~ condition, data = Mac_complement)
cat("Macrophage Complement - p-value:", signif(t_mac_complement$p.value, 3), "\n")

# ── HIF-1 Signaling Gene Set (redefine for macrophage context) ─────────────────
HIF1_mac_full <- c(
  "AKT3", "CDKN1A", "CDKN1B", "EGLN2", "EGLN3", "CREBBP", "CYBB", "LDHAL6A", "EDN1",
  "EGF", "EGFR", "EIF4E", "EIF4EBP1", "ENO1", "ENO2", "ENO3", "EP300", "EPAS1", "EPO",
  "ERBB2", "AKT1", "AKT2", "ALDOA", "ALDOB", "ALDOC", "FLT1", "MTOR", "EIF4E1B",
  "GAPDH", "ANGPT1", "ANGPT2", "MKNK2", "HIF1A", "HK1", "HK2", "HK3", "HMOX1", "IFNG",
  "IFNGR1", "IFNGR2", "IGF1", "IGF1R", "IL6", "IL6R", "INS", "INSR", "LDHA", "LDHB",
  "LDHC", "ARNT", "LTBR", "NFKB1", "NOS2", "NOS3", "NPPA", "SERPINE1", "ANGPT4",
  "PDHA1", "PDHA2", "PDHB", "PDK1", "PFKFB3", "PFKL", "PFKM", "PFKP", "PGK1", "PGK2",
  "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2", "PLCG1", "PLCG2", "EGLN1",
  "PRKCA", "PRKCB", "PRKCG", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "BCL2", "RELA",
  "RPS6", "RPS6KB1", "RPS6KB2", "SLC2A1", "STAT3", "ELOC", "ELOB", "TEK", "TF", "TFRC",
  "TIMP1", "TLR4", "VEGFA", "VHL", "HKDC1", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G",
  "CUL2", "PIK3R3", "MKNK1", "LDHAL6B", "EIF4E2", "RBX1"
)

HIF1_mac <- intersect(HIF1_mac_full, rownames(pseudo_Mac))
cat("Macrophage HIF-1 genes found:", length(HIF1_mac), "\n")

# ── Add Module Score ───────────────────────────────────────────────────────────
pseudo_Mac <- AddModuleScore(pseudo_Mac, list(HIF1_mac), name = "HIF1_", assay = "RNA")

# ── Build Data Frame ───────────────────────────────────────────────────────────
Mac_hif1 <- pseudo_Mac@meta.data %>%
  mutate(sample = rownames(pseudo_Mac@meta.data)) %>%
  transmute(
    sample = sample,
    cell_type = "Macrophage",
    module = "HIF1_Signaling",
    score = HIF1_1,
    condition = factor(condition, levels = c("Control", "PAH"))
  ) %>%
  filter(!is.na(condition))

# ── Statistical Test ───────────────────────────────────────────────────────────
t_mac_hif1 <- t.test(score ~ condition, data = Mac_hif1)
cat("Macrophage HIF-1 - p-value:", signif(t_mac_hif1$p.value, 3), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# UPDATE COMBINED RESULTS
# ═══════════════════════════════════════════════════════════════════════════════

combined_scores <- bind_rows(
  Hep_cyp450, Hep_hif1, Hep_aeroglycolysis, Hep_oxphos, Hep_gluconeo, Hep_fam, Hep_jakstat,
  HSC_pi3k, HSC_hif1, HSC_ecm, HSC_jakstat,
  Endo_hif1, Endo_pi3k, Endo_oxphos, Endo_aeroglycolysis, Endo_ppp, Endo_jakstat,
  Mac_complement, Mac_hif1
) %>%
  arrange(cell_type, module, condition, sample)

summary_stats <- combined_scores %>%
  group_by(cell_type, module, condition) %>%
  summarise(
    n = n(),
    mean = mean(score),
    median = median(score),
    sd = sd(score),
    se = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(p_value = case_when(
    cell_type == "Hepatocyte" & module == "CYP450_Oxidation" ~ t_hep$p.value,
    cell_type == "Hepatocyte" & module == "HIF1_Signaling" ~ t_hep_hif1$p.value,
    cell_type == "Hepatocyte" & module == "Aerobic_Glycolysis" ~ t_hep_aeroglycolysis$p.value,
    cell_type == "Hepatocyte" & module == "OxPhos" ~ t_hep_oxphos$p.value,
    cell_type == "Hepatocyte" & module == "Gluconeogenesis" ~ t_hep_gluconeo$p.value,
    cell_type == "Hepatocyte" & module == "Fatty_Acid_Metabolism" ~ t_hep_fam$p.value,
    cell_type == "Hepatocyte" & module == "JAK_STAT" ~ t_hep_jakstat$p.value,
    cell_type == "HSC_mFB" & module == "PI3K_Akt" ~ t_hsc_pi3k$p.value,
    cell_type == "HSC_mFB" & module == "HIF1_Signaling" ~ t_hsc_hif1$p.value,
    cell_type == "HSC_mFB" & module == "ECM_Receptor" ~ t_hsc_ecm$p.value,
    cell_type == "HSC_mFB" & module == "JAK_STAT" ~ t_hsc_jakstat$p.value,
    cell_type == "Endothelial" & module == "HIF1_Signaling" ~ t_endo_hif1$p.value,
    cell_type == "Endothelial" & module == "PI3K_Akt" ~ t_endo_pi3k$p.value,
    cell_type == "Endothelial" & module == "OXPHOS" ~ t_endo_oxphos$p.value,
    cell_type == "Endothelial" & module == "Aerobic_Glycolysis" ~ t_endo_aeroglycolysis$p.value,
    cell_type == "Endothelial" & module == "PPP" ~ t_endo_ppp$p.value,
    cell_type == "Endothelial" & module == "JAK_STAT" ~ t_endo_jakstat$p.value,
    cell_type == "Macrophage" & module == "Complement_Coagulation" ~ t_mac_complement$p.value,
    cell_type == "Macrophage" & module == "HIF1_Signaling" ~ t_mac_hif1$p.value
  ))

# ── Export ─────────────────────────────────────────────────────────────────────
write.csv(combined_scores, "PAH_Liver_Pseudobulk_Module_Scores_by_Sample.csv", row.names = FALSE)
write.csv(summary_stats, "PAH_Liver_Pseudobulk_Module_Summary.csv", row.names = FALSE)

cat("\n✓ Exported: PAH_Liver_Pseudobulk_Module_Scores_by_Sample.csv")
cat("\n✓ Exported: PAH_Liver_Pseudobulk_Module_Summary.csv\n")

print(combined_scores)
print(summary_stats)