# Configuration for TNBC analysis.R
# Gene marker sets, thresholds, and cluster labels used throughout the pipeline.

cell_state_markers <- list(
    Proliferative = c("MKI67", "PCNA", "TOP2A", "CCNB1", "AURKA", "CDK1", "CENPF"),
    Epithelial_like = c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "GATA3"),
    Basal_like = c("KRT5", "KRT14", "ACTA2", "MYLK", "SNAI2"),
    Mesenchymal_like = c("VIM", "CDH2", "FN1", "SNAI1", "TWIST1", "ZEB1", "ZEB2"),
    Stressed = c("HSPA1A", "HSPA1B", "DNAJA1", "HSP90AA1", "HSPB1", "XBP1", "HERPUD1"),
    Cycling = c("MKI67", "PCNA", "TOP2A", "CCNB1", "AURKA", "BIRC5"),
    CD8_T_cells = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7"),
    CD4_T_cells = c("CD4", "IL7R", "FOXP3", "CTLA4", "ICOS"),
    Tregs = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"),
    B_cells = c("CD19", "MS4A1", "CD79A", "IGHM", "IGKC"),
    Macrophages = c("CD68", "CD163", "CSF1R", "MRC1", "MSR1"),
    TAMs = c("CD68", "CD163", "MRC1", "APOE", "C1QA", "C1QB", "TREM2"),
    Dendritic_cells = c("CD1C", "FCER1A", "CLEC9A", "XCR1", "LAMP3"),
    NK_cells = c("NKG7", "GNLY", "KLRD1", "NCR1", "FCGR3A"),
    Myeloid_cells = c("CD14", "LYZ", "AIF1", "FCN1", "S100A8", "S100A9"),
    Fibroblasts = c("PDGFRA", "PDGFRB", "COL1A1", "COL1A2", "FAP", "SFRP2"),
    CAFs = c("FAP", "SFRP2", "COL10A1", "MMP11", "POSTN", "ACTA2"),
    Pericytes = c("RGS5", "MCAM", "PDGFRB", "NOTCH3"),
    Endothelial = c("PECAM1", "VWF", "CD34", "KDR", "CLDN5")
)

antigen_presentation_genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1",
                                 "HLA-DPB1", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9")

lineage_markers <- list(
    Tumor = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "GATA3"),
    Immune = c("PTPRC", "CD3D", "CD3E", "CD8A", "CD4", "CD19", "CD68", "CD14", "NKG7"),
    Stromal = c("PDGFRA", "PDGFRB", "COL1A1", "COL1A2", "FAP", "VIM", "ACTA2")
)

tumor_markers <- c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "GATA3")
immune_markers <- c("PTPRC", "CD3D", "CD3E", "CD8A", "CD4", "CD19", "CD68", "CD14", "NKG7")
stromal_markers <- c("PDGFRA", "PDGFRB", "COL1A1", "COL1A2", "FAP", "PECAM1", "VWF", "CD34")

candidate_markers <- c("SSR4", "XBP1", "HERPUD1", "HSP90B1", "SEC11C", "FKBP11",
                       "EPCAM", "KRT7", "KRT19", "VIM", "CD14", "MGP", "TGFBR2", "HLA-DRA", "B2M")

stromal_validation_markers <- c("PDGFRA", "PDGFRB", "COL1A1", "COL1A2",
                                 "FAP", "VIM", "ACTA2", "SFRP2")

immune_clusters <- c("T_cells", "Myeloid_Cells", "Activated_Immune_Cells", "B_Cells_Plasma")
caf_cluster <- "CAF"

# NOTE: this cluster -> label mapping reflects the identities found for THIS
# specific dataset/run at resolution = 0.2 (10 clusters, numbered 0-9 by
# Seurat's FindClusters). It is not a general-purpose mapping: if you rerun
# clustering on different data, a different Seurat/R version, or a different
# resolution, cluster numbers can come out in a different order and this
# mapping will need to be redone by inspecting marker genes per cluster.
TNBC.cluster.ids <- c(
    "0" = "Tumor_Epithelial",
    "1" = "T_cells",
    "2" = "Tumor_proliferative",
    "3" = "Myeloid_Cells",
    "4" = "Metabolic_tumor",
    "5" = "ER_stressed_AP_suppressed_plasma_cells",
    "6" = "CAF",
    "7" = "Activated_Immune_Cells",
    "8" = "Endothelial_Cells",
    "9" = "B_Cells_Plasma"
)

cluster_order_corrected <- c(
    "ER_stressed_AP_suppressed_plasma_cells",
    "CAF",
    "Tumor_Epithelial",
    "Tumor_proliferative",
    "Metabolic_tumor",
    "Endothelial_Cells",
    "Activated_Immune_Cells",
    "T_cells",
    "Myeloid_Cells",
    "B_Cells_Plasma"
)
