# TNBC SINGLE-CELL RNA SEQUENCE ANALYSIS - COMPLETE VERSION

# Step 1: Setup and packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

packages_needed <- c("Seurat", "ggplot2", "dplyr", "tidyr", "pheatmap",
                     "RColorBrewer", "gridExtra", "cowplot", "patchwork", "viridis", "ggpubr")

for(pkg in packages_needed) {
    if (!require(pkg, character.only = TRUE)) {
        if(pkg == "Seurat") {
            BiocManager::install(pkg)
        } else {
            install.packages(pkg)
        }
        library(pkg, character.only = TRUE)
    }
}

# Step 2: Set working directory and data loading
setwd("C:/Users/user/OneDrive/Desktop/me/uni/y4 s2/Research/Seurat Object")

if(file.exists("SeuratObject_TNBC.rds")) {
    tnbc_obj <- readRDS("SeuratObject_TNBC.rds")
    tnbc_obj <- UpdateSeuratObject(tnbc_obj)
    cat("TNBC object loaded successfully!\n")
    cat("Dimensions:", dim(tnbc_obj), "\n")
    cat("Number of cells:", ncol(tnbc_obj), "\n")
} else {
    stop("SeuratObject_TNBC.rds not found in working directory")
}

dir.create("analysis_output_tnbc", showWarnings = FALSE)
setwd("analysis_output_tnbc")

# Step 3: Define gene sets
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

cat("Total cell state markers:", length(unlist(cell_state_markers)), "\n")
cat("Total antigen presentation genes:", length(antigen_presentation_genes), "\n")

# Step 4: Calculating cell state and AP scores
expr_matrix <- GetAssayData(tnbc_obj, layer = "data")

for(state in names(cell_state_markers)) {
    genes <- intersect(cell_state_markers[[state]], rownames(expr_matrix))
    if(length(genes) >= 2) {
        tnbc_obj[[paste0("Score_", state)]] <- colMeans(expr_matrix[genes, ], na.rm = TRUE)
        cat("  Added", state, "-", length(genes), "genes\n")
    }
}

ap_genes <- intersect(antigen_presentation_genes, rownames(expr_matrix))
tnbc_obj$AP_score <- colMeans(expr_matrix[ap_genes, ], na.rm = TRUE)
tnbc_obj$AP_suppressed <- tnbc_obj$AP_score < quantile(tnbc_obj$AP_score, 0.25, na.rm = TRUE)
cat("  Added AP score -", length(ap_genes), "genes\n")
cat("  AP-suppressed cells:", sum(tnbc_obj$AP_suppressed), 
    "(", round(mean(tnbc_obj$AP_suppressed)*100, 1), "%)\n")

# Step 5: Dimension reduction and clustering
tnbc_obj <- FindVariableFeatures(tnbc_obj)
tnbc_obj <- ScaleData(tnbc_obj)
tnbc_obj <- RunPCA(tnbc_obj, features = VariableFeatures(tnbc_obj), npcs = 50)
tnbc_obj <- RunUMAP(tnbc_obj, dims = 1:10)
tnbc_obj <- FindNeighbors(tnbc_obj, dims = 1:10)
tnbc_obj <- FindClusters(tnbc_obj, resolution = 0.2)

cat("Number of clusters:", length(unique(tnbc_obj$seurat_clusters)), "\n")
print(table(tnbc_obj$seurat_clusters))

# Step 6: Identify suppressor cluster and marker genes
cluster_ranking <- tnbc_obj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(
        n_cells = n(),
        mean_AP = round(mean(AP_score, na.rm = TRUE), 3),
        pct_suppressed = round(mean(AP_suppressed, na.rm = TRUE) * 100, 1)
    ) %>%
    arrange(desc(pct_suppressed))

print(cluster_ranking)
top_cluster <- cluster_ranking$seurat_clusters[1]
cat("\n>>> Cluster", top_cluster, "has highest AP suppression at",
    cluster_ranking$pct_suppressed[1], "%\n")

all_markers <- FindAllMarkers(tnbc_obj, only.pos = TRUE,
                               min.pct = 0.25, logfc.threshold = 0.25)

top10_markers_per_cluster <- all_markers %>%
    group_by(cluster) %>%
    slice_head(n = 10) %>%
    arrange(cluster)
write.csv(top10_markers_per_cluster, "top10_markers_per_cluster.csv", row.names = FALSE)

top_cluster_markers <- all_markers %>%
    filter(cluster == top_cluster) %>%
    slice_head(n = 20)
print(top_cluster_markers)

# Step 7: Cluster characterization
get_cluster_scores <- function(obj, gene_list) {
    results <- data.frame()
    for(cl in unique(obj$seurat_clusters)) {
        cells_in_cl <- obj$seurat_clusters == cl
        row_data <- data.frame(seurat_clusters = cl)
        for(set_name in names(gene_list)) {
            genes <- intersect(gene_list[[set_name]], rownames(obj))
            if(length(genes) > 0) {
                expr <- GetAssayData(obj, layer = "data")[genes, cells_in_cl, drop = FALSE]
                avg_score <- mean(colMeans(expr, na.rm = TRUE), na.rm = TRUE)
                row_data[[set_name]] <- round(avg_score, 3)
            } else {
                row_data[[set_name]] <- NA
            }
        }
        results <- rbind(results, row_data)
    }
    return(results)
}

lineage_markers <- list(
    Tumor = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "GATA3"),
    Immune = c("PTPRC", "CD3D", "CD3E", "CD8A", "CD4", "CD19", "CD68", "CD14", "NKG7"),
    Stromal = c("PDGFRA", "PDGFRB", "COL1A1", "COL1A2", "FAP", "VIM", "ACTA2")
)

lineage_scores <- get_cluster_scores(tnbc_obj, lineage_markers)
functional_scores <- get_cluster_scores(tnbc_obj, cell_state_markers[1:6])

cluster_characterization <- cluster_ranking %>%
    left_join(lineage_scores, by = "seurat_clusters") %>%
    left_join(functional_scores, by = "seurat_clusters") %>%
    arrange(desc(pct_suppressed))

cluster_characterization <- cluster_characterization %>%
    rowwise() %>%
    mutate(
        dominant_lineage = case_when(
            Tumor > Immune & Tumor > Stromal ~ "Tumor",
            Immune > Tumor & Immune > Stromal ~ "Immune",
            Stromal > Tumor & Stromal > Immune ~ "Stromal",
            TRUE ~ "Mixed"
        ),
        dominant_state = case_when(
            Proliferative >= Epithelial_like & Proliferative >= Basal_like & 
            Proliferative >= Mesenchymal_like & Proliferative >= Stressed & 
            Proliferative >= Cycling ~ "Proliferative",
            Epithelial_like >= Proliferative & Epithelial_like >= Basal_like & 
            Epithelial_like >= Mesenchymal_like & Epithelial_like >= Stressed & 
            Epithelial_like >= Cycling ~ "Epithelial_like",
            Basal_like >= Proliferative & Basal_like >= Epithelial_like & 
            Basal_like >= Mesenchymal_like & Basal_like >= Stressed & 
            Basal_like >= Cycling ~ "Basal_like",
            Mesenchymal_like >= Proliferative & Mesenchymal_like >= Epithelial_like & 
            Mesenchymal_like >= Basal_like & Mesenchymal_like >= Stressed & 
            Mesenchymal_like >= Cycling ~ "Mesenchymal_like",
            Stressed >= Proliferative & Stressed >= Epithelial_like & 
            Stressed >= Basal_like & Stressed >= Mesenchymal_like & 
            Stressed >= Cycling ~ "Stressed",
            TRUE ~ "Cycling"
        )
    ) %>%
    ungroup()

write.csv(cluster_characterization, "cluster_characterization.csv", row.names = FALSE)

# CLUSTER NAMES 
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

tnbc_obj <- RenameIdents(tnbc_obj, TNBC.cluster.ids)
tnbc_obj$cluster_name_corrected <- Idents(tnbc_obj)

# Cluster order for visualization (ordered by AP suppression)
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

tnbc_obj$cluster_name_ordered <- factor(tnbc_obj$cluster_name_corrected, 
                                         levels = cluster_order_corrected)

# Verify name changes
print(table(tnbc_obj$cluster_name_corrected))


# Step 8: Correlation analysis within suppressor cluster
cells_in_top <- tnbc_obj$seurat_clusters == top_cluster
ap_suppressed_in_top <- tnbc_obj$AP_suppressed[cells_in_top]
expr_matrix_top <- GetAssayData(tnbc_obj, layer = "data")[, cells_in_top]

correlations <- data.frame()
for(gene in rownames(expr_matrix_top)) {
    gene_expr <- expr_matrix_top[gene, ]
    if(sd(gene_expr) > 0) {
        cor_test <- cor.test(gene_expr, as.numeric(ap_suppressed_in_top), method = "spearman")
        correlations <- rbind(correlations, data.frame(
            gene = gene,
            correlation = cor_test$estimate,
            p_value = cor_test$p.value
        ))
    }
}
correlations <- correlations[order(-abs(correlations$correlation)), ]

top_cluster_genes <- head(all_markers[all_markers$cluster == top_cluster, "gene"], 100)
top_corr_genes <- head(correlations$gene[correlations$correlation > 0], 50)
candidate_genes <- intersect(top_cluster_genes, top_corr_genes)
cat("\nCandidate genes associated with AP suppression in ER-stressed cells:\n")
print(candidate_genes)

# Step 9: Cell type scores
tumor_markers <- c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "GATA3")
immune_markers <- c("PTPRC", "CD3D", "CD3E", "CD8A", "CD4", "CD19", "CD68", "CD14", "NKG7")
stromal_markers <- c("PDGFRA", "PDGFRB", "COL1A1", "COL1A2", "FAP", "PECAM1", "VWF", "CD34")

safe_colMeans <- function(mat, genes) {
    genes_present <- intersect(genes, rownames(mat))
    if(length(genes_present) == 0) return(rep(NA, ncol(mat)))
    sub_mat <- mat[genes_present, , drop = FALSE]
    colMeans(sub_mat, na.rm = TRUE)
}

tnbc_obj$score_tumor <- safe_colMeans(expr_matrix, intersect(tumor_markers, rownames(tnbc_obj)))
tnbc_obj$score_immune <- safe_colMeans(expr_matrix, intersect(immune_markers, rownames(tnbc_obj)))
tnbc_obj$score_stromal <- safe_colMeans(expr_matrix, intersect(stromal_markers, rownames(tnbc_obj)))

tnbc_obj$cell_type <- case_when(
    tnbc_obj$score_tumor > tnbc_obj$score_immune & tnbc_obj$score_tumor > tnbc_obj$score_stromal ~ "Tumor",
    tnbc_obj$score_immune > tnbc_obj$score_tumor & tnbc_obj$score_immune > tnbc_obj$score_stromal ~ "Immune",
    tnbc_obj$score_stromal > tnbc_obj$score_tumor & tnbc_obj$score_stromal > tnbc_obj$score_immune ~ "Stromal",
    TRUE ~ "Mixed"
)

# Step 10: Patient analysis - AP score across patients
ap_by_patient <- tnbc_obj@meta.data %>%
    group_by(group) %>%
    summarise(
        n_cells = n(),
        mean_AP = round(mean(AP_score, na.rm = TRUE), 3),
        sd_AP = round(sd(AP_score, na.rm = TRUE), 3),
        pct_suppressed = round(mean(AP_suppressed, na.rm = TRUE) * 100, 1)
    ) %>%
    arrange(desc(mean_AP))
print(ap_by_patient)

p_ap_patient <- ggplot(ap_by_patient, aes(x = reorder(group, -mean_AP), 
                                           y = mean_AP, fill = mean_AP)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = mean_AP - sd_AP, ymax = mean_AP + sd_AP), width = 0.2) +
    geom_text(aes(label = round(mean_AP, 2)), vjust = -0.5, size = 3) +
    scale_fill_gradient(low = "steelblue", high = "red", name = "Mean AP Score") +
    theme_minimal(base_size = 10) +
    labs(title = "Antigen Presentation Score by Patient",
         x = "Patient", y = "Mean AP Score ± SD") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("AP_Score_by_Patient.png", p_ap_patient, width = 12, height = 6, dpi = 300)

# Step 11: UMAP split by patient
if("umap" %in% Reductions(tnbc_obj)) {
    p_split_wide <- DimPlot(tnbc_obj, reduction = "umap", pt.size = 0.2, 
                             split.by = "group", ncol = 4) +
        ggtitle("UMAP of TNBC Cells Split by Patient") +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
              strip.text = element_text(face = "bold", size = 8))
    ggsave("UMAP_Split_by_Patient_Wide.png", p_split_wide, width = 20, height = 16, dpi = 300)
}

# Step 12: Immune cell infiltration percentage per patient
immune_clusters <- c("T_cells", "Myeloid_Cells", "Activated_Immune_Cells", "B_Cells_Plasma")

immune_infiltration <- tnbc_obj@meta.data %>%
    group_by(group) %>%
    summarise(
        total_cells = n(),
        immune_cells = sum(cluster_name_corrected %in% immune_clusters, na.rm = TRUE),
        pct_immune = round(immune_cells / total_cells * 100, 2)
    ) %>%
    arrange(desc(pct_immune))
print(immune_infiltration)

p_immune <- ggplot(immune_infiltration, aes(x = reorder(group, -pct_immune), 
                                             y = pct_immune, fill = pct_immune)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", pct_immune)), vjust = -0.5, size = 3) +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "% Immune Cells") +
    theme_minimal(base_size = 10) +
    labs(title = "Immune Cell Infiltration by Patient",
         x = "Patient", y = "Percentage of Immune Cells") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("Immune_Infiltration_by_Patient.png", p_immune, width = 12, height = 6, dpi = 300)

# Step 13: Cancer-Associated Fibroblasts percentage per patient
caf_cluster <- "CAF"

caf_by_patient <- tnbc_obj@meta.data %>%
    group_by(group) %>%
    summarise(
        total_cells = n(),
        caf_cells = sum(cluster_name_corrected == caf_cluster, na.rm = TRUE),
        pct_caf = round(caf_cells / total_cells * 100, 2)
    ) %>%
    arrange(desc(pct_caf))
print(caf_by_patient)

p_caf <- ggplot(caf_by_patient, aes(x = reorder(group, -pct_caf), 
                                     y = pct_caf, fill = pct_caf)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", pct_caf)), vjust = -0.5, size = 3) +
    scale_fill_gradient(low = "lightgreen", high = "darkgreen", name = "% CAFs") +
    theme_minimal(base_size = 10) +
    labs(title = "Cancer-Associated Fibroblasts (CAFs) by Patient",
         x = "Patient", y = "Percentage of CAFs") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("CAFs_by_Patient.png", p_caf, width = 12, height = 6, dpi = 300)

# Step 14: Correlation between immune infiltration and AP suppression
correlation_data <- immune_infiltration %>%
    left_join(ap_by_patient, by = "group") %>%
    select(group, pct_immune, mean_AP, pct_suppressed)
print(correlation_data)

if(nrow(correlation_data) >= 3) {
    cor_test <- cor.test(correlation_data$pct_immune, 
                          correlation_data$pct_suppressed, 
                          method = "pearson")
    
    cat("\nPearson correlation (Immune Infiltration vs % AP-Suppressed):\n")
    cat("  r =", round(cor_test$estimate, 3), "\n")
    cat("  p-value =", format(cor_test$p.value, scientific = TRUE), "\n")
    
    p_corr <- ggplot(correlation_data, aes(x = pct_immune, y = pct_suppressed)) +
        geom_point(size = 4, alpha = 0.7, color = "steelblue") +
        geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink") +
        geom_text(aes(label = group), vjust = -1, hjust = 0.5, size = 3) +
        annotate("text", x = max(correlation_data$pct_immune) * 0.7, 
                 y = max(correlation_data$pct_suppressed) * 0.9,
                 label = paste("r =", round(cor_test$estimate, 3),
                               "\np =", format(cor_test$p.value, digits = 3)),
                 size = 4, hjust = 0) +
        theme_minimal(base_size = 12) +
        labs(title = "Correlation: Immune Infiltration vs AP Suppression",
             x = "Immune Cell Infiltration (%)", 
             y = "AP-Suppressed Cells (%)") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave("Correlation_Immune_vs_AP_Suppression.png", p_corr, 
           width = 8, height = 6, dpi = 300)
}

# Step 15: Bar plot showing percentage of each cluster per patient
cluster_patient_pct <- tnbc_obj@meta.data %>%
    group_by(group, cluster_name_ordered) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(pct = n_cells / sum(n_cells) * 100)

ordered_clusters <- levels(tnbc_obj$cluster_name_ordered)

p_cluster_patient <- ggplot(cluster_patient_pct, 
                            aes(x = group, y = pct, fill = factor(cluster_name_ordered, levels = ordered_clusters))) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_viridis_d(name = "Cluster") +
    theme_minimal(base_size = 10) +
    labs(title = "Cluster Composition by Patient",
         x = "Patient", y = "Percentage of Cells (%)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "right")
ggsave("Cluster_Composition_by_Patient.png", p_cluster_patient, 
       width = 14, height = 8, dpi = 300)

# Step 16: Heatmap with actual cluster names
heatmap_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_head(n = 5) %>%
    pull(gene) %>%
    unique()
heatmap_markers <- intersect(heatmap_markers, rownames(tnbc_obj))

if(length(heatmap_markers) > 0) {
    avg_exp <- AverageExpression(tnbc_obj, features = heatmap_markers, 
                                  group.by = "cluster_name_ordered")$RNA
    avg_exp_scaled <- t(scale(t(avg_exp)))
    
    pheatmap(avg_exp_scaled,
             main = "Top 5 Marker Genes per Cluster (Named Clusters)",
             color = colorRampPalette(c("blue", "white", "red"))(100),
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             border_color = NA,
             filename = "Heatmap_Named_Clusters.png",
             width = 12, height = 10)
    cat("Heatmap with named clusters saved successfully!\n")
}

# Step 17: Visualizations (UMAP, t-SNE, AP score plots)
if(!"umap" %in% Reductions(tnbc_obj)) {
    cat("ERROR: UMAP not found. Please run RunUMAP() first.\n")
} else {

# 17.1: UMAP colored by cluster name
p_umap_named <- DimPlot(tnbc_obj, reduction = "umap", group.by = "cluster_name_ordered", 
                         label = TRUE, repel = TRUE, label.size = 5) +
    ggtitle("A. TNBC Cell Clusters (Ordered by Suppression)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("1_UMAP_by_Cluster_Named.png", p_umap_named, width = 14, height = 10, dpi = 300)

# 17.2: t-SNE colored by cluster
if(!"tsne" %in% Reductions(tnbc_obj)) {
    cat("  Running t-SNE for visualization...\n")
    tnbc_obj <- RunTSNE(tnbc_obj, dims = 1:10, perplexity = 30, check_duplicates = FALSE)
}
p_tsne <- DimPlot(tnbc_obj, reduction = "tsne", group.by = "cluster_name_ordered", label = TRUE, repel = TRUE, label.size = 4) +
    ggtitle("B. t-SNE of TNBC Cells (Colored by Cluster)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("2_tSNE_by_Cluster.png", p_tsne, width = 12, height = 10, dpi = 300)

# 17.4: AP suppression status on UMAP
tnbc_obj$AP_status <- ifelse(tnbc_obj$AP_suppressed, "AP-Suppressed", "Normal AP")
p_ap_status <- DimPlot(tnbc_obj, reduction = "umap", group.by = "AP_status",
                        cols = c("AP-Suppressed" = "red", "Normal AP" = "lightgray")) +
    ggtitle("D. AP-Suppressed Cells on UMAP") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("4_AP_Status_UMAP.png", p_ap_status, width = 10, height = 8, dpi = 300)

# 17.5a: Bar plot of AP suppression by cluster (NAMED)
cluster_ranking <- cluster_ranking %>%
    arrange(desc(pct_suppressed)) %>%
    mutate(rank = 1:n())

# Updated cluster name mapping for barplot
cluster_name_mapping <- c(
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

cluster_ranking$cluster_name <- cluster_name_mapping[as.character(cluster_ranking$seurat_clusters)]

p_ap_bar <- ggplot(cluster_ranking, aes(x = reorder(cluster_name, -pct_suppressed),
                                         y = pct_suppressed, 
                                         fill = pct_suppressed)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", pct_suppressed)), 
              vjust = -0.5, size = 4) +
    scale_fill_gradient(low = "steelblue", high = "red", name = "% Suppressed") +
    theme_minimal(base_size = 12) +
    labs(title = "E. Antigen Presentation Suppression by Cell Type",
         x = "Cell Type", 
         y = "% of Cells with Low AP Score") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
ggsave("5_AP_Suppression_Barplot_Named.png", p_ap_bar, width = 12, height = 7, dpi = 300)

# 17.5b: Bar plot of AP suppression by RANK
p_ap_bar_rank <- ggplot(cluster_ranking, aes(x = factor(rank), y = pct_suppressed, fill = pct_suppressed)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", pct_suppressed)), vjust = -0.5, size = 5) +
    scale_fill_gradient(low = "steelblue", high = "red", name = "% Suppressed") +
    scale_x_discrete(labels = as.character(1:10)) +
    theme_minimal(base_size = 14) +
    labs(title = "Antigen Presentation Suppression by Rank",
         x = "Rank (1 = Highest Suppression)", y = "% of Cells with Low AP Score") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(size = 12))
ggsave("5_AP_Suppression_Barplot_Ranked.png", p_ap_bar_rank, width = 10, height = 6, dpi = 300)

# 17.6: Top cluster highlighted on UMAP
tnbc_obj$highlight <- ifelse(tnbc_obj$seurat_clusters == top_cluster, 
                              paste("Cluster", top_cluster, "(Highest Suppression)"), 
                              "Other Clusters")
highlight_colors <- c("red", "lightgray")
names(highlight_colors) <- c(paste("Cluster", top_cluster, "(Highest Suppression)"), "Other Clusters")
p_highlight <- DimPlot(tnbc_obj, reduction = "umap", group.by = "highlight", cols = highlight_colors) +
    ggtitle(paste("F. Cluster", top_cluster, "- Primary Suppressor Population")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("6_Top_Cluster_Highlighted.png", p_highlight, width = 10, height = 8, dpi = 300)

# 17.7a: Heatmap of top markers
if(length(heatmap_markers) > 0) {
    avg_exp <- AverageExpression(tnbc_obj, features = heatmap_markers, group.by = "cluster_name_ordered")$RNA
    avg_exp_scaled <- t(scale(t(avg_exp)))
    pheatmap(avg_exp_scaled,
             main = "G. Top 5 Marker Genes per Cluster",
             color = colorRampPalette(c("blue", "white", "red"))(100),
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             border_color = NA,
             filename = "7_Marker_Genes_Heatmap.png",
             width = 12, height = 10)
}

# 17.7b: Candidate marker heatmap
candidate_markers <- c("SSR4", "XBP1", "HERPUD1", "HSP90B1", "SEC11C", "FKBP11", 
                       "EPCAM", "KRT7", "KRT19", "VIM", "CD14", "MGP", "TGFBR2", "HLA-DRA", "B2M")
candidate_markers <- intersect(candidate_markers, rownames(tnbc_obj))

if(length(candidate_markers) > 0) {
    avg_exp <- AverageExpression(tnbc_obj, features = candidate_markers, group.by = "cluster_name_ordered")$RNA
    avg_exp_scaled <- t(scale(t(avg_exp)))
    pheatmap(avg_exp_scaled,
             main = "Candidate Marker Genes",
             color = colorRampPalette(c("blue", "white", "red"))(100),
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             border_color = NA,
             filename = "8_Candidate_Markers_Heatmap.png",
             width = 10, height = 8)
}

# 17.8: Feature plots for top candidate markers
if(length(candidate_genes) >= 1) {
    p_candidate <- FeaturePlot(tnbc_obj, features = candidate_genes[1], reduction = "umap") +
        scale_color_gradientn(colors = c("lightgray", "yellow", "red")) +
        ggtitle(paste("Expression of", candidate_genes[1])) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave("10_Top_Candidate_Expression.png", p_candidate, width = 10, height = 8, dpi = 300)
}

if(length(candidate_genes) >= 2) {
    p_second <- FeaturePlot(tnbc_obj, features = candidate_genes[2], reduction = "umap") +
        scale_color_gradientn(colors = c("lightgray", "yellow", "red")) +
        ggtitle(paste("Expression of", candidate_genes[2])) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave("11_Second_Candidate_Expression.png", p_second, width = 10, height = 8, dpi = 300)
}

} # End UMAP check

# Step 18: Correlation analyses
# 18.1: Correlation: AP suppression % vs CAFs % per patient
if(exists("caf_by_patient") && exists("ap_by_patient")) {
    correlation_caf_ap <- caf_by_patient %>%
        left_join(ap_by_patient, by = "group") %>%
        select(group, pct_caf, pct_suppressed, mean_AP)
    
    if(nrow(correlation_caf_ap) >= 3) {
        cor_test_caf <- cor.test(correlation_caf_ap$pct_caf, 
                                   correlation_caf_ap$pct_suppressed, 
                                   method = "pearson")
        
        p_corr_caf_ap <- ggplot(correlation_caf_ap, aes(x = pct_caf, y = pct_suppressed)) +
            geom_point(size = 5, alpha = 0.7, color = "darkgreen") +
            geom_smooth(method = "lm", se = TRUE, color = "red", fill = "lightgreen") +
            geom_text(aes(label = group), vjust = -1.2, hjust = 0.5, size = 4) +
            annotate("text", x = max(correlation_caf_ap$pct_caf) * 0.7, 
                     y = max(correlation_caf_ap$pct_suppressed) * 0.9,
                     label = paste("r =", round(cor_test_caf$estimate, 3),
                                   "\np =", format(cor_test_caf$p.value, digits = 3)),
                     size = 5, hjust = 0, fontface = "bold") +
            theme_minimal(base_size = 14) +
            labs(title = "Correlation: CAFs vs AP Suppression",
                 x = "CAFs (%)", 
                 y = "AP-Suppressed Cells (%)") +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
        ggsave("8_Correlation_CAFs_vs_AP_Suppression.png", p_corr_caf_ap, width = 9, height = 7, dpi = 300)
    }
}

# 18.2: Correlation: Immune infiltration % vs CAFs %
if(exists("immune_infiltration") && exists("caf_by_patient")) {
    correlation_immune_caf <- immune_infiltration %>%
        left_join(caf_by_patient, by = "group") %>%
        select(group, pct_immune, pct_caf)
    
    if(nrow(correlation_immune_caf) >= 3) {
        cor_test_immune_caf <- cor.test(correlation_immune_caf$pct_immune, 
                                          correlation_immune_caf$pct_caf, 
                                          method = "pearson")
        
        p_corr_immune_caf <- ggplot(correlation_immune_caf, aes(x = pct_immune, y = pct_caf)) +
            geom_point(size = 5, alpha = 0.7, color = "purple") +
            geom_smooth(method = "lm", se = TRUE, color = "red", fill = "plum") +
            geom_text(aes(label = group), vjust = -1.2, hjust = 0.5, size = 4) +
            annotate("text", x = max(correlation_immune_caf$pct_immune) * 0.7, 
                     y = max(correlation_immune_caf$pct_caf) * 0.9,
                     label = paste("r =", round(cor_test_immune_caf$estimate, 3),
                                   "\np =", format(cor_test_immune_caf$p.value, digits = 3)),
                     size = 5, hjust = 0, fontface = "bold") +
            theme_minimal(base_size = 14) +
            labs(title = "Correlation: Immune Infiltration vs CAFs",
                 x = "Immune Cell Infiltration (%)", 
                 y = "CAFs (%)") +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
        ggsave("9_Correlation_Immune_vs_CAFs.png", p_corr_immune_caf, width = 9, height = 7, dpi = 300)
    }
}

# 18.3: Correlation: CAFs vs Immune cells
if(exists("caf_by_patient") && exists("immune_infiltration")) {
    correlation_caf_immune <- caf_by_patient %>%
        left_join(immune_infiltration, by = "group") %>%
        select(group, pct_caf, pct_immune)
    
    if(nrow(correlation_caf_immune) >= 3) {
        cor_test_caf_immune <- cor.test(correlation_caf_immune$pct_caf, 
                                         correlation_caf_immune$pct_immune, 
                                         method = "pearson")
        
        p_corr_caf_immune <- ggplot(correlation_caf_immune, aes(x = pct_caf, y = pct_immune)) +
            geom_point(size = 5, alpha = 0.7, color = "orange") +
            geom_smooth(method = "lm", se = TRUE, color = "red", fill = "yellow") +
            geom_text(aes(label = group), vjust = -1.2, hjust = 0.5, size = 4) +
            annotate("text", x = max(correlation_caf_immune$pct_caf) * 0.7, 
                     y = max(correlation_caf_immune$pct_immune) * 0.9,
                     label = paste("r =", round(cor_test_caf_immune$estimate, 3),
                                   "\np =", format(cor_test_caf_immune$p.value, digits = 3)),
                     size = 5, hjust = 0, fontface = "bold") +
            theme_minimal(base_size = 14) +
            labs(title = "Correlation: CAFs vs Immune Cells",
                 x = "CAFs (%)", 
                 y = "Immune Cells (%)") +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
        ggsave("9b_Correlation_CAF_vs_Immune.png", p_corr_caf_immune, width = 9, height = 7, dpi = 300)
        
        cat("\n=== CAF vs Immune Correlation ===\n")
        cat("r =", round(cor_test_caf_immune$estimate, 3), "\n")
        cat("p-value =", format(cor_test_caf_immune$p.value, scientific = TRUE), "\n")
    }
}

# Step 19: marker validation for Cluster 5
stromal_validation_markers <- c("PDGFRA", "PDGFRB", "COL1A1", "COL1A2", 
                                 "FAP", "VIM", "ACTA2", "SFRP2")
stromal_markers_present <- intersect(stromal_validation_markers, rownames(tnbc_obj))
cat("  Stromal markers present:", paste(stromal_markers_present, collapse = ", "), "\n")

if(length(stromal_markers_present) >= 2) {
    n_markers <- min(length(stromal_markers_present), 6)
    
    p_stromal_features <- FeaturePlot(tnbc_obj, features = stromal_markers_present[1:n_markers], 
                                        reduction = "umap", ncol = 3) +
        scale_color_gradientn(colors = c("lightgray", "yellow", "red")) &
        theme_minimal() &
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
    
    ggsave("10_Stromal_Markers_Validation.png", p_stromal_features, width = 15, height = 12, dpi = 300)
    cat("  Saved: 10_Stromal_Markers_Validation.png\n")
    
    cluster5_cells_corrected <- tnbc_obj$cluster_name_corrected == "ER_stressed_AP_suppressed_plasma_cells"
    
    stromal_scores <- data.frame()
    for(marker in stromal_markers_present) {
        marker_expr <- FetchData(tnbc_obj, marker)[,1]
        cluster5_mean <- mean(marker_expr[cluster5_cells_corrected], na.rm = TRUE)
        other_mean <- mean(marker_expr[!cluster5_cells_corrected], na.rm = TRUE)
        stromal_scores <- rbind(stromal_scores, data.frame(
            Marker = marker,
            Cluster5_Mean = round(cluster5_mean, 4),
            Other_Clusters_Mean = round(other_mean, 4),
            Fold_Change = round(cluster5_mean / (other_mean + 0.0001), 2)
        ))
    }
    
    cat("\n  Stromal marker expression in Cluster 5 vs other clusters:\n")
    cat("  (LOW expression confirms Cluster 5 is NOT stromal!)\n")
    print(stromal_scores)
    
    write.csv(stromal_scores, "Cluster5_Stromal_Marker_Validation.csv", row.names = FALSE)
    cat("  Saved: Cluster5_Stromal_Marker_Validation.csv\n")
} else {
    cat("  Insufficient stromal markers found for validation\n")
}

# Step 20: Feature plot validations
DefaultAssay(tnbc_obj) <- "RNA"

# ER stress markers (Cluster 5 specific)
FeaturePlot(tnbc_obj, features = c("XBP1"))
ggsave("FeaturePlot_XBP1.png", width = 10, height = 8, dpi = 300)

FeaturePlot(tnbc_obj, features = c("HERPUD1"))
ggsave("FeaturePlot_HERPUD1.png", width = 10, height = 8, dpi = 300)

FeaturePlot(tnbc_obj, features = c("HSP90B1"))
ggsave("FeaturePlot_HSP90B1.png", width = 10, height = 8, dpi = 300)

FeaturePlot(tnbc_obj, features = c("SSR4"))
ggsave("FeaturePlot_SSR4.png", width = 10, height = 8, dpi = 300)

# Epithelial marker
FeaturePlot(tnbc_obj, features = c("EPCAM"))
ggsave("FeaturePlot_EPCAM.png", width = 10, height = 8, dpi = 300)

# Pan-immune marker (CD45)
FeaturePlot(tnbc_obj, features = c("PTPRC"))
ggsave("FeaturePlot_PTPRC.png", width = 10, height = 8, dpi = 300)

# T cell markers (Cluster 1 - T_cells)
FeaturePlot(tnbc_obj, features = c("CD3D", "CD3E", "CD8A", "CD4"))
ggsave("Validation_Tcell_Markers.png", width = 14, height = 10, dpi = 300)

# Proliferation marker
FeaturePlot(tnbc_obj, features = c("MKI67"))
ggsave("FeaturePlot_MKI67.png", width = 10, height = 8, dpi = 300)

# B cell marker
FeaturePlot(tnbc_obj, features = c("CD19"))
ggsave("FeaturePlot_CD19.png", width = 10, height = 8, dpi = 300)

# Myeloid marker
FeaturePlot(tnbc_obj, features = c("CD14"))
ggsave("FeaturePlot_CD14.png", width = 10, height = 8, dpi = 300)

# Mesenchymal marker
FeaturePlot(tnbc_obj, features = c("VIM"))
ggsave("FeaturePlot_VIM.png", width = 10, height = 8, dpi = 300)

# CAF markers (Cluster 6)
FeaturePlot(tnbc_obj, features = c("FAP", "SFRP2", "COL1A1", "POSTN"))
ggsave("Validation_CAF_Markers.png", width = 14, height = 10, dpi = 300)

# Endothelial markers
FeaturePlot(tnbc_obj, features = c("PECAM1"))
ggsave("FeaturePlot_PECAM1.png", width = 10, height = 8, dpi = 300)

FeaturePlot(tnbc_obj, features = c("VWF"))
ggsave("FeaturePlot_VWF.png", width = 10, height = 8, dpi = 300)

# B cell/Plasma cell markers
FeaturePlot(tnbc_obj, features = c("MS4A1"))
ggsave("FeaturePlot_CD20.png", width = 10, height = 8, dpi = 300)

FeaturePlot(tnbc_obj, features = c("IGKC"))
ggsave("FeaturePlot_IGKC.png", width = 10, height = 8, dpi = 300)

FeaturePlot(tnbc_obj, features = c("JCHAIN"))
ggsave("FeaturePlot_JCHAIN.png", width = 10, height = 8, dpi = 300)

# Combined feature plot for EPCAM and PTPRC (CD45) in one image
p_combined_features <- FeaturePlot(tnbc_obj, features = c("EPCAM", "PTPRC"), 
                                    reduction = "umap", ncol = 2) &
    scale_color_gradientn(colors = c("lightgray", "yellow", "red")) &
    theme_minimal() &
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
ggsave("FeaturePlot_EPCAM_PTPRC_Combined.png", p_combined_features, width = 14, height = 6, dpi = 300)
cat("✓ Saved combined feature plot: FeaturePlot_EPCAM_PTPRC_Combined.png\n")

# Additional validation plots for specific clusters
# Tumor_Epithelial markers (Cluster 0)
FeaturePlot(tnbc_obj, features = c("EPCAM", "KRT8", "KRT18", "KRT19"))
ggsave("Validation_Tumor_Epithelial_Markers.png", width = 14, height = 10, dpi = 300)

# Tumor_proliferative markers (Cluster 2)
FeaturePlot(tnbc_obj, features = c("MKI67", "TOP2A", "PCNA", "CCNB1"))
ggsave("Validation_Tumor_Proliferative_Markers.png", width = 14, height = 10, dpi = 300)

# Metabolic_tumor markers (Cluster 4)
FeaturePlot(tnbc_obj, features = c("HMGCS1", "SLC2A1", "LDHA", "PKM"))
ggsave("Validation_Metabolic_Tumor_Markers.png", width = 14, height = 10, dpi = 300)

# Step 21: CANDIDATE MARKERS FOR TUMOR_EPITHELIAL (Cluster 0) AND CAFs (Cluster 6) 
# Get cluster numbers
cluster_numbers <- tnbc_obj@meta.data %>%
    select(seurat_clusters, cluster_name_corrected) %>%
    distinct()

tumor_epithelial_cluster_num <- cluster_numbers$seurat_clusters[cluster_numbers$cluster_name_corrected == "Tumor_Epithelial"][1]
caf_cluster_num <- cluster_numbers$seurat_clusters[cluster_numbers$cluster_name_corrected == "CAF"][1]

# Candidate markers for Tumor_Epithelial (Cluster 0)
if(!is.na(tumor_epithelial_cluster_num)) {
    cat("\n--- Candidate Markers for Tumor_Epithelial (Original Cluster", tumor_epithelial_cluster_num, ") ---\n")
    
    tumor_epithelial_candidate_markers <- all_markers %>%
        filter(cluster == tumor_epithelial_cluster_num) %>%
        arrange(desc(avg_log2FC)) %>%
        head(20)
    
    print(tumor_epithelial_candidate_markers[, c("gene", "avg_log2FC", "pct.1", "p_val_adj")])
    write.csv(tumor_epithelial_candidate_markers, "Candidate_Markers_Tumor_Epithelial.csv", row.names = FALSE)
    cat("✓ Saved: Candidate_Markers_Tumor_Epithelial.csv\n")
    
    # Feature plot for top Tumor_Epithelial markers
    top_tumor_epithelial <- head(tumor_epithelial_candidate_markers$gene, 6)
    p_tumor_epithelial_features <- FeaturePlot(tnbc_obj, features = top_tumor_epithelial, reduction = "umap", ncol = 3) &
        scale_color_gradientn(colors = c("lightgray", "yellow", "red")) &
        theme_minimal()
    ggsave("Top_Tumor_Epithelial_Markers_FeaturePlot.png", p_tumor_epithelial_features, width = 15, height = 12, dpi = 300)
    cat("✓ Saved: Top_Tumor_Epithelial_Markers_FeaturePlot.png\n")
}

# Candidate markers for CAFs (Cluster 6)
if(!is.na(caf_cluster_num)) {
    cat("\n--- Candidate Markers for CAFs (Original Cluster", caf_cluster_num, ") ---\n")
    
    caf_candidate_markers <- all_markers %>%
        filter(cluster == caf_cluster_num) %>%
        arrange(desc(avg_log2FC)) %>%
        head(20)
    
    print(caf_candidate_markers[, c("gene", "avg_log2FC", "pct.1", "p_val_adj")])
    write.csv(caf_candidate_markers, "Candidate_Markers_CAFs.csv", row.names = FALSE)
    cat("✓ Saved: Candidate_Markers_CAFs.csv\n")
    
    # Feature plot for top CAF markers
    top_caf <- head(caf_candidate_markers$gene, 6)
    p_caf_features <- FeaturePlot(tnbc_obj, features = top_caf, reduction = "umap", ncol = 3) &
        scale_color_gradientn(colors = c("lightgray", "yellow", "red")) &
        theme_minimal()
    ggsave("Top_CAF_Markers_FeaturePlot.png", p_caf_features, width = 15, height = 12, dpi = 300)
    cat("✓ Saved: Top_CAF_Markers_FeaturePlot.png\n")
}

# Step 22: FINAL SUMMARY

top_cluster_original <- as.character(cluster_ranking$seurat_clusters[1])
top_pct_original <- cluster_ranking$pct_suppressed[1]
top_n_cells <- cluster_ranking$n_cells[1]

# Get lineage and state 
top_lineage <- unique(cluster_characterization$dominant_lineage[cluster_characterization$seurat_clusters == top_cluster_original])
top_state <- unique(cluster_characterization$dominant_state[cluster_characterization$seurat_clusters == top_cluster_original])

cat("1. PRIMARY SUPPRESSOR CLUSTER:\n")
cat("   - Name: ER_stressed_AP_suppressed_plasma_cells\n")
cat("   - AP suppression rate:", cluster_ranking$pct_suppressed[1], "%\n")
cat("   - Number of cells:", format(top_n_cells, big.mark = ","), "\n")
cat("   - Dominant functional state:", top_state, "\n")
cat("   - Top markers: XBP1, HERPUD1, HSP90B1, FKBP11, SSR4 (ER stress genes)\n\n")


cat("2. CORRELATION RESULTS (based on PTPRC/CD45 feature plots):\n")
if(exists("cor_test")) {
    cat("   - Immune Infiltration vs AP Suppression: r =", round(cor_test$estimate, 3),
        " (p =", format(cor_test$p.value, digits = 3), ")\n")
}
if(exists("cor_test_caf")) {
    cat("   - CAFs vs AP Suppression: r =", round(cor_test_caf$estimate, 3), 
        " (p =", format(cor_test_caf$p.value, digits = 3), ")\n")
}
if(exists("cor_test_caf_immune")) {
    cat("   - CAFs vs Immune Cells: r =", round(cor_test_caf_immune$estimate, 3),
        " (p =", format(cor_test_caf_immune$p.value, digits = 3), ")\n\n")
}

cat("3. CANDIDATE MARKERS:\n")
if(exists("tumor_epithelial_candidate_markers") && nrow(tumor_epithelial_candidate_markers) > 0) {
    cat("   Tumor_Epithelial (Cluster 0):\n")
    for(i in 1:min(5, nrow(tumor_epithelial_candidate_markers))) {
        cat("     -", tumor_epithelial_candidate_markers$gene[i], 
            "(log2FC =", round(tumor_epithelial_candidate_markers$avg_log2FC[i], 2), ")\n")
    }
}
if(exists("caf_candidate_markers") && nrow(caf_candidate_markers) > 0) {
    cat("\n   CAFs (Cluster 6):\n")
    for(i in 1:min(5, nrow(caf_candidate_markers))) {
        cat("     -", caf_candidate_markers$gene[i], 
            "(log2FC =", round(caf_candidate_markers$avg_log2FC[i], 2), ")\n")
    }
}

cat("\n4. PATIENT ANALYSIS:\n")
cat("   - Total patients analyzed:", length(unique(tnbc_obj$group)), "\n")
cat("   - Patients:", paste(unique(tnbc_obj$group), collapse = ", "), "\n")

cat("\n5. FILES SAVED:\n")
cat("   - UMAPs: 1_UMAP_by_Cluster_Named.png, 2_tSNE_by_Cluster.png, UMAP_Split_by_Patient_Wide.png\n")
cat("   - Bar plots: 5_AP_Suppression_Barplot_Named.png, 5_AP_Suppression_Barplot_Ranked.png\n")
cat("   - Heatmaps: 7_Marker_Genes_Heatmap.png, 8_Candidate_Markers_Heatmap.png, Heatmap_Named_Clusters.png\n")
cat("   - Feature plots: FeaturePlot_XBP1.png, FeaturePlot_HERPUD1.png, FeaturePlot_SSR4.png\n")
cat("   - Combined feature plot: FeaturePlot_EPCAM_PTPRC_Combined.png\n")
cat("   - Correlation plots: 8_Correlation_CAFs_vs_AP_Suppression.png, 9_Correlation_Immune_vs_CAFs.png, 9b_Correlation_CAF_vs_Immune.png\n")
cat("   - Patient analysis: AP_Score_by_Patient.png, Immune_Infiltration_by_Patient.png, CAFs_by_Patient.png, Cluster_Composition_by_Patient.png\n")
cat("   - Candidate markers: Candidate_Markers_Tumor_Epithelial.csv, Candidate_Markers_CAFs.csv, Top_Tumor_Epithelial_Markers_FeaturePlot.png, Top_CAF_Markers_FeaturePlot.png\n")
cat("   - Validation plots: Validation_Tumor_Epithelial_Markers.png, Validation_Tumor_Proliferative_Markers.png, Validation_Metabolic_Tumor_Markers.png, Validation_Tcell_Markers.png, Validation_CAF_Markers.png\n")
cat("   - Data files: top10_markers_per_cluster.csv, cluster_characterization.csv, Cluster5_Stromal_Marker_Validation.csv\n")
