# ==============================================================================
# File: manual_annotation.R
# Description: performing manual annotation by using the infromation of the
# provided gene panel
# ==============================================================================

# Setting up environment =======================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation
# Set seed
set.seed(42)

target_dir <- "/path/to/file/"
.libPaths(target_dir)

# Loading relevant libraries  ==================================================
library(tidyverse)
library(celldex)
library(SingleR)
library(Seurat)
library(scRNAseq)
library(scran)
library(viridis)
library(pheatmap)
library(Matrix)
library(HPAanalyze)
library(hpar)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(SimBu)
library(gridExtra)
library(ggplot2)
library(CellChat)
library(gmp)
library(NMF)
library(ggalluvial)
library(reticulate)
library(grid)
library(gridBase)
library(arrow)
library(classInt)
library(sf)
library(remotes)
library(RColorBrewer)
library(scales)
library(fields)
library(reshape2)
library(scales)
library(ComplexHeatmap)
library(circlize)

# Main  ========================================================================

# load seurat object and split it
seu <- LoadSeuratRds("/path/to/file/")

# split the seu object
patients <- unique(seu$patient_id)

seurat_list <- list()
for (patient in patients) {
  seurat_list[[patient]] <- subset(seu, subset = patient_id == patient)
}

cell_genes <- list(
  t_cells = c("CD3E", "CD3G", "CD4", "CD8A", "CD8B", "GZMB", "IL4", "IL5", "RORC", "GATA3", "KLRG1", "RGS1", "SLAMF6", "MKI67", "FOXP3", "TCF7", "IRF1", "BCL6", "CMKLR1", "ENTPD1"),
  mhc_molecules = c("HLA.DQA1", "HLA.DRB1", "HLA.E", "HLA.DRA", "HLA.B"),
  trm_cells = c("TGFB1", "CD69", "ITGAE", "ITGA1", "CXCR6", "EOMES", "IL15"),
  keratinocytes = c("KRT14", "KRT20", "KRT13", "KRT9", "KRT12", "KRT1", "KRT8", "KRT2", "KRT6A", "KRT7", "KRT19", "KRT18", "KRT3", "KRT4", "KRT10"),
  b_cells = c("CD22", "CD79A", "IGLL1", "FCER2", "MS4A1", "CD72", "IGLL5", "CD19", "TCL1A"),
  checkpoints = c("LAG3", "CTLA4", "CD27", "CD86", "REL", "HAVCR2", "CD200", "ICOS", "ICOSLG", "CD28", "PDCD1", "CD274", "CD96", "CD40", "TIGIT", "PDCD1LG2", "IDO1", "TNFRSF9"),
  nk_cells = c("NKG7", "GZMK", "IL12B", "IL12RB2", "GNLY", "GZMA", "CD226", "GZMH", "KLRB1", "NCR1"),
  oncogenes_tumor_suppressor_genes = c("DNAJC12", "CTSC", "CLEC2B", "ALDOC", "MBOAT7", "CAMK1", "ANXA1", "THBS1", "VWF", "CDC45", "BLVRB", "CD44", "MCAM", "STAT3", "SCARB2", "STAP1", "VPS37B", "TP53", "ENDOD1", "CA2", "S100A1", "S100P", "VAMP5", "MAP2K1", "GIMAP7", "SAMD9L", "RRS1", "BASP1", "POGLUT3", "PDLIM1", "RORA", "HMGB1", "PLSCR1", "STAT1", "SIRPG", "PTPRC", "CTNNB1", "TAP1", "SLC39A8", "HMGN5", "PSMB10", "SPATS2L", "XAF1", "PRG4", "SEMA4A", "S100A8", "S100A9", "GBP5", "GBP1", "SH2D1A", "RUNX2", "CEP55", "ADAMTS13", "VEGFA", "SNX30", "AGER", "AUH", "USP18", "FBXO6", "CASK", "RGCC", "DDX58", "MRPS26", "PAX8", "MX1", "SIRT1", "PIK3C2B", "FCGR3A", "RUNX1", "ANO6", "CASP4", "ANXA2", "NELL2", "GBP4", "GBP2", "CMPK2", "PAEP", "RXRA", "BRAF", "LEF1", "PALLD", "DDX60L", "CASP1", "LIMA1", "CDH16", "SCPEP1", "SPINT1", "SPINT2", "MYO1F", "RASGRP4", "PDPN", "SERPINA1", "PRF1"),
  cell_adhesion_migration = c("ICAM1", "SELL", "ITGB7", "VCAM1", "PECAM1", "ITGA4", "LMCD1", "MADCAM1"),
  macrophages = c("CD68", "LYVE1", "BCAT1", "IL1B", "CD276", "MAOA", "CD74", "ARG1", "TREM2", "CD83", "SAMD9", "PSTPIP2", "CD80", "CD163", "ITGAM", "MRC1", "NOS2"),
  chemokines = c("CCL22", "CCR4", "CCL2", "CCR4", "CCR5", "CCL13", "CXCL12", "CXCR4", "CCL19", "CCR7", "CXCL9", "CXCL10", "CXCL11", "CXCR3", "CXCL13", "CXCR5", "CXCL5", "CXCR6", "CCL11", "CCR5", "CXCL8", "XCL1", "XCR1", "CCL20", "CCR6", "CCL28", "XCL1", "XCR1", "CXCL1", "CXCL14", "CXCR4", "CXCL17", "CCL5"),
  dendritic_cells = c("BATF3", "IRF4", "FLT3", "LAMP3", "ITGAX", "XCR1", "CLEC9A", "CD1C", "RSAD2"),
  cytokines = c("IL23A", "IFNG", "IL21", "IL23R", "IL10", "IL36G", "IL36A", "IL7", "IL18", "IL4", "IL5", "IL3", "TNFSF14", "IL12A", "IL17F", "IL17A", "IL22", "IL1RN", "IL27", "IFNB1", "IL15", "IL6", "IL33", "TNF", "IL32"),
  others = c("NTAN1", "CD38", "WDFY1", "HERC6", "TLR4", "SOX10", "TNFRSF8", "SMAD2"),
  inf_inducable_genes = c("OAS3", "LGALS3BP", "IFI44L", "IFIT5", "IFIT1", "IFIT3", "IFI44", "HELZ2", "TRANK1", "HERC5"),
  complement = c("C5", "CR2")
)

Idents(seu) <- "seurat_clusters"

seu <- subset(seu, features = rownames(seu)[!grepl("^Blank", rownames(seu))])
all_markers <- FindAllMarkers(seu)
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 5) %>%
  arrange(cluster, p_val) %>%
  group_by(cluster) %>%
  slice_max(order_by = abs(avg_log2FC), n = 5) %>%
  arrange(cluster, desc(abs(avg_log2FC)))

unique_genes <- unique(top_markers$gene)

# Create a binary matrix indicating whether a gene is a marker for a cluster
marker_matrix <- matrix(0, nrow = length(unique_genes), ncol = length(unique(top_markers$cluster)))
rownames(marker_matrix) <- unique_genes
colnames(marker_matrix) <- unique(top_markers$cluster)

for (i in 1:nrow(top_markers)) {
  gene <- top_markers$gene[i]
  cluster <- top_markers$cluster[i]
  marker_matrix[gene, cluster] <- 1
}

marker_matrix_t <- t(marker_matrix)

seu <- subset(seu, features = rownames(seu)[!grepl("^Blank", rownames(seu))])
all_markers <- FindAllMarkers(seu)
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

unique_genes <- unique(top_markers$gene)

# Create a binary matrix indicating whether a gene is a marker for a cluster
marker_matrix <- matrix(0, nrow = length(unique_genes), ncol = length(unique(top_markers$cluster)))
rownames(marker_matrix) <- unique_genes
colnames(marker_matrix) <- unique(top_markers$cluster)

for (i in 1:nrow(top_markers)) {
  gene <- top_markers$gene[i]
  cluster <- top_markers$cluster[i]
  marker_matrix[gene, cluster] <- 1
}

# Create a color mapping for each cell type
cell_types <- names(cell_genes)
colors <- rainbow(length(cell_types))
names(colors) <- cell_types

# Create a matrix indicating the cell type for each gene
cell_type_matrix <- matrix(NA, nrow = nrow(marker_matrix), ncol = ncol(marker_matrix))
rownames(cell_type_matrix) <- rownames(marker_matrix)
colnames(cell_type_matrix) <- colnames(marker_matrix)

for (cell_type in cell_types) {
  genes <- cell_genes[[cell_type]]
  cell_type_matrix[rownames(cell_type_matrix) %in% genes, ] <- cell_type
}

# Define colors for absent (0) and present (1) with cell type colors
col_fun <- function(x) {
  ifelse(is.na(x), "white", colors[x])
}

cell_type_matrix_t <- t(cell_type_matrix)

marker_matrix_t <- t(marker_matrix)

# Define a color mapping for cell types
cell_types <- unique(as.vector(cell_type_matrix_t))
colors <- rainbow(length(cell_types))
names(colors) <- cell_types

# remove KRT16, CD3D
genes_to_remove <- c("STAT4", "CD3D", "KRT16", "MX2", "S1PR1")

marker_matrix_t <- marker_matrix_t[, !colnames(marker_matrix_t) %in% genes_to_remove]
cell_type_matrix_t <- cell_type_matrix_t[, !colnames(cell_type_matrix_t) %in% genes_to_remove]

# Convert the data frames back to matrices
marker_matrix_t <- as.matrix(marker_matrix_t)
cell_type_matrix_t <- as.matrix(cell_type_matrix_t)

cell_types_factor <- factor(unique(unlist(cell_type_matrix_t[1,])))

numeric_cell_type_matrix <- apply(cell_type_matrix_t, c(1, 2), function(x) as.numeric(cell_types_factor[match(x, levels(cell_types_factor))]))

numeric_cell_type_matrix[marker_matrix_t != 1] <- 0

new_matrix <- cell_type_matrix_t
new_matrix[marker_matrix_t == 0] <- "NA"

df <- as.data.frame(as.table(new_matrix))

colnames(df) <- c("Cluster", "Gene", "CellType")

df <- df %>% mutate(CellType = ifelse(CellType == "None", "NA", CellType))

# other kind of plot
filtered_df <- df %>% filter(CellType != "NA")

# Calculate the number of genes of each cell type for each cluster
result_df <- filtered_df %>%
  group_by(Cluster, CellType) %>%
  summarise(GeneCount = n()) %>%
  ungroup()

# Check the structure of the new data frame
str(result_df)

heatmap_data <- result_df %>%
  spread(key = Cluster, value = GeneCount, fill = 0)

# Convert the data frame to a matrix for the heatmap
heatmap_matrix <- as.matrix(heatmap_data[,-1])
rownames(heatmap_matrix) <- heatmap_data$CellType

cell_types <- c("b_cells", "dendritic_cells", "keratinocytes", "macrophages", "nk_cells", "t_cells", "trm_cells")
heatmap_matrix_just_cell_types <- heatmap_matrix[rownames(heatmap_matrix) %in% cell_types, ]

# Create the heatmap
pdf(paste0("/path/to/file/"),
    width = 7,
    height = 3)
ggplot(melt(heatmap_matrix_just_cell_types), aes(x = factor(Var2), y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "Cluster", y = "Cell Type", fill = "Marker Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank())
dev.off()

# assigning the known cluster - cell type associations
manual_annotations <- rep(NA, ncol(seu))

# Assign annotations based on Seurat clusters
manual_annotations[seu$seurat_clusters == 0] <- "Macrophage"
manual_annotations[seu$seurat_clusters %in% c(1, 2)] <- "Keratinocytes"
manual_annotations[seu$seurat_clusters == 5] <- "T_cells"

# Add the manual_annotations vector as a new column to the metadata
seu <- AddMetaData(seu, metadata = manual_annotations, col.name = "manual_annotations")

# Verify the new column
head(seu@meta.data)

SaveSeuratRds(seu, "/path/to/file/")

pdf(file = "/path/to/file/",
    width = 10,
    height = 7)
Idents(seu) <- "manual_annotations"
umap_plot <- DimPlot(seu, reduction = "umap")
umap_plot + labs(x = "UMAP 1", y = "UMAP 2")
dev.off()

# individual p-values
marker_matrix <- heatmap_matrix_just_cell_types

# Step 1: Assign each cluster the cell type with the highest value
most_probable_cell_types <- apply(marker_matrix, 2, function(x) rownames(marker_matrix)[which.max(x)])

# Step 2: Perform a statistical test for each cluster to calculate the p-value
p_values <- sapply(1:ncol(marker_matrix), function(cluster) {
  cell_type <- most_probable_cell_types[cluster]
  cell_type_index <- which(rownames(marker_matrix) == cell_type)
  
  # Create a contingency table for the current cluster
  contingency_table <- matrix(c(marker_matrix[cell_type_index, cluster], sum(marker_matrix[cell_type_index, ]) - marker_matrix[cell_type_index, cluster],
                                sum(marker_matrix[, cluster]) - marker_matrix[cell_type_index, cluster], sum(marker_matrix) - sum(marker_matrix[cell_type_index, ]) - sum(marker_matrix[, cluster]) + marker_matrix[cell_type_index, cluster]),
                              nrow = 2)
  
  # Perform Fisher's exact test
  fisher_test <- fisher.test(contingency_table)
  
  # Return the p-value
  return(fisher_test$p.value)
})

# Print the results
results <- data.frame(Cluster = colnames(marker_matrix), Most_Probable_Cell_Type = most_probable_cell_types, P_Value = p_values)
print(results)


# p-value when using the highest log2FC as the assignment
filtered_cell_genes <- cell_genes[cell_types]

# Step 1: Assign each cluster to the cell type with the highest avg_log2FC for its marker genes
assign_clusters <- function(all_markers, filtered_cell_genes) {
  cluster_assignments <- data.frame(Cluster = unique(all_markers$cluster), Cell_Type = NA, Max_avg_log2FC = NA)
  
  for (cluster in unique(all_markers$cluster)) {
    cluster_markers <- all_markers %>% filter(cluster == !!cluster)
    max_avg_log2FC <- -Inf
    best_cell_type <- NA
    
    for (cell_type in names(filtered_cell_genes)) {
      cell_type_genes <- filtered_cell_genes[[cell_type]]
      avg_log2FC <- mean(cluster_markers %>% filter(gene %in% cell_type_genes) %>% pull(avg_log2FC), na.rm = TRUE)
      
      if (!is.na(avg_log2FC) && avg_log2FC > max_avg_log2FC) {
        max_avg_log2FC <- avg_log2FC
        best_cell_type <- cell_type
      }
    }
    
    cluster_assignments <- cluster_assignments %>% mutate(Cell_Type = ifelse(Cluster == cluster, best_cell_type, Cell_Type),
                                                          Max_avg_log2FC = ifelse(Cluster == cluster, max_avg_log2FC, Max_avg_log2FC))
  }
  
  return(cluster_assignments)
}

cluster_assignments <- assign_clusters(all_markers, filtered_cell_genes)

# Step 2: Perform a statistical test for each cluster to calculate the p-value
calculate_p_values <- function(all_markers, cluster_assignments, filtered_cell_genes) {
  p_values <- sapply(1:nrow(cluster_assignments), function(i) {
    cluster <- cluster_assignments$Cluster[i]
    cell_type <- cluster_assignments$Cell_Type[i]
    cell_type_genes <- filtered_cell_genes[[cell_type]]
    
    cluster_markers <- all_markers %>% filter(cluster == !!cluster)
    positive_genes <- sum(cluster_markers$gene %in% cell_type_genes & cluster_markers$avg_log2FC > 0)
    negative_genes <- sum(cluster_markers$gene %in% cell_type_genes & cluster_markers$avg_log2FC <= 0)
    other_genes <- nrow(cluster_markers) - positive_genes - negative_genes
    
    contingency_table <- matrix(c(positive_genes, negative_genes, other_genes, nrow(all_markers) - positive_genes - negative_genes - other_genes), nrow = 2)
    fisher_test <- fisher.test(contingency_table)
    
    return(fisher_test$p.value)
  })
  
  cluster_assignments <- cluster_assignments %>% mutate(P_Value = p_values)
  return(cluster_assignments)
}

cluster_assignments <- calculate_p_values(all_markers, cluster_assignments, filtered_cell_genes)

# Print the results
print(cluster_assignments)

png(paste0("/path/to/file/"),
    width = 12,
    height = 6,
    units = "in",
    res = 600)
ht <- Heatmap(
  matrix = numeric_cell_type_matrix,  # Use the modified numeric matrix
  name = "Cell Types",  # Legend name
  col = unname(color_mapping),  # Use the color mapping
  rect_gp = gpar(col = "white"),  # Set cell border color
  cluster_rows = FALSE,  # Disable row clustering
  cluster_columns = FALSE,  # Disable column clustering
  show_row_names = TRUE,  # Show row names
  show_column_names = TRUE,  # Show column names
  row_names_gp = gpar(fontsize = 10),  # Customize row name font size
  column_names_gp = gpar(fontsize = 10)
)
dev.off()


Idents(seu) <- "seurat_clusters"
# umap
DimPlot(seu)

png(paste0("/path/to/file/"),
    width = 12,
    height = 8,
    units = "in",
    res = 600)
FeaturePlot(seu, cell_genes$b_cells)
dev.off()

b_marker <- c("CD74", "HLA.DRA", "HLA.B", "HLA.E", "CTNNB1")
png(paste0("/path/to/file/"),
    width = 8,
    height = 8,
    units = "in",
    res = 600)
FeaturePlot(seu, b_marker)
dev.off()


