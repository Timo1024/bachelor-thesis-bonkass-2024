# ==============================================================================
# File: analyze_cell_atlas_annotations.R
# Description: Analysis of the annotations inferred by SingleR.
# Genes associated to specific cell types described in the gene panel were
# used to assess the annotation quality
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
library(ggrepel)
library(ComplexHeatmap)

# Main  ========================================================================

# load seurat object
seu <- LoadSeuratRds("/path/to/file/")

Idents(seu) <- "SingleR_HCA"

# plot each cell type in a umap
idents <- unique(Idents(seu))
plots <- list()
for (ident in idents) {
  cells.highlight <- WhichCells(seu, idents = ident)
  p <- DimPlot(seu, 
               reduction = "umap", 
               cells.highlight = cells.highlight,
               sizes.highlight = 0.5,
               cols.highlight = "#3333dd",
               pt.size = 0.5) +
    ggtitle(ident) +
    theme(legend.position = "none")
  plots[[ident]] <- p
}
combined_plot <- cowplot::plot_grid(plotlist = plots, ncol = 3)

png(paste0("/path/to/file/"),
    width = 13,
    height = 7,
    units = "in",
    res = 300)
print(combined_plot)
dev.off()

plot_umap_by_celltype <- function(seu, celltype, responder_col = "responder_status", celltype_col = "SingleR_HCA") {
  plot_color <- rep("grey", nrow(seu@meta.data))
  plot_color[seu@meta.data[[celltype_col]] == celltype] <- seu@meta.data[[responder_col]][seu@meta.data[[celltype_col]] == celltype]
  seu <- AddMetaData(seu, metadata = plot_color, col.name = "plot_color")
  seu$plot_color <- factor(seu$plot_color, levels = c("responder", "non_responder", "grey"))
  DimPlot(seu, reduction = "umap", group.by = "plot_color", order = c("responder", "non_responder", "grey")) +
    scale_color_manual(values = c("responder" = "blue", "non_responder" = "red", "grey" = "#cccccc"),
                       labels = c("responder" = "Responder", "non_responder" = "Non Responder", "grey" = "Other")) +
    ggtitle(celltype) +
    labs(x = "UMAP 1", y = "UMAP 2")
}
cell_types <- unique(seu$SingleR_HCA)
umap_plots <- lapply(cell_types, function(celltype) {
  plot_umap_by_celltype(seu, celltype)
})
combined_plot <- wrap_plots(umap_plots, ncol = 3) +
  plot_layout(guides = "collect")

png(paste0("/path/to/file/"),
    width = 13,
    height = 7,
    units = "in",
    res = 300)
print(combined_plot)
dev.off()

# just for macs
plot_umap_by_celltype <- function(seu, celltype = "Macrophage", responder_col = "responder_status", celltype_col = "SingleR_HCA") {
  plot_color <- rep("grey", nrow(seu@meta.data))
  plot_color[seu@meta.data[[celltype_col]] == celltype] <- seu@meta.data[[responder_col]][seu@meta.data[[celltype_col]] == celltype]
  seu <- AddMetaData(seu, metadata = plot_color, col.name = "plot_color")
  seu$plot_color <- factor(seu$plot_color, levels = c("responder", "non_responder", "grey"))
  DimPlot(seu, reduction = "umap", group.by = "plot_color", order = c("responder", "non_responder", "grey")) +
    scale_color_manual(values = c("responder" = "blue", "non_responder" = "red", "grey" = "#cccccc"),
                       labels = c("responder" = "Responder", "non_responder" = "Non Responder", "grey" = "Other")) +
    ggtitle(celltype) +
    labs(x = "UMAP 1", y = "UMAP 2")
}
p <- plot_umap_by_celltype(seu)
png(paste0("/path/to/file/"),
    width = 8,
    height = 5,
    units = "in",
    res = 300)
print(p)
dev.off()

# visualize amount of cells per cell_type (for responders and non responders)
# for general annotations
metadata <- seu@meta.data
cell_counts_patients <- metadata %>%
  group_by(SingleR_HCA, patient_id) %>%
  summarise(count = n()) %>%
  ungroup()

# Calculate total counts for each responder status
total_counts_patients <- metadata %>%
  group_by(patient_id) %>%
  summarise(total = n())

# Join total counts with cell counts
cell_counts_patients <- cell_counts_patients %>%
  left_join(total_counts_patients, by = "patient_id")
  # left_join(total_counts, by = "responder_status")

# Calculate proportions
cell_counts_patients <- cell_counts_patients %>%
  mutate(proportion = count / total)

cell_counts_responders <- metadata %>%
  group_by(SingleR_HCA, responder_status) %>%
  summarise(count = n()) %>%
  ungroup()

# Calculate total counts for each responder status
total_counts_responders <- metadata %>%
  group_by(responder_status) %>%
  summarise(total = n())

# Join total counts with cell counts
cell_counts_responders <- cell_counts_responders %>%
  left_join(total_counts_responders, by = "responder_status")

# Calculate proportions
cell_counts_responders <- cell_counts_responders %>%
  mutate(proportion = count / total)

cell_counts_patients <- cell_counts_patients %>%
  mutate(group = paste0("patient ", patient_id)) %>%
  select(-patient_id)
cell_counts_responders <- cell_counts_responders %>%
  mutate(group = responder_status) %>%
  select(-responder_status) %>%
  mutate(group = ifelse(group == "non_responder", "non responder", group))

head(cell_counts_patients)
head(cell_counts_responders)

# merge dfs
cell_counts <- rbind(cell_counts_patients, cell_counts_responders)

# prepare for heatmap
heatmap_data <- reshape2::dcast(cell_counts, SingleR_HCA ~ group, value.var = "proportion")
rownames(heatmap_data) <- heatmap_data$SingleR_HCA
heatmap_data$SingleR_HCA <- NULL

breaks = breaks = c(seq(0, 0.2, length.out = 50), seq(0.2, 1, length.out = 50))
colors = viridis(100, option = "viridis")
col_fun = circlize::colorRamp2(breaks, colors)

heatmap_data <- heatmap_data %>%
  select(`patient 0`, `patient 1`, `patient 2`, `patient 3`, responder, `non responder`)

# heatmap
png(paste0("/path/to/file/"),
    width = 7,
    height = 7,
    units = "in",
    res = 300)
Heatmap(
  as.matrix(heatmap_data),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  col = col_fun,
  # column_title = "Responder Status",
  column_title = "Groups",
  name = "Proportion",
  column_title_side = "bottom",
  column_names_rot = 45,
  column_names_centered = TRUE,
  row_names_gp = gpar(fontsize = 10), 
  row_names_max_width = unit(4, "cm"),
  column_split = c(rep(1, 4), rep(2, 2)),
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- heatmap_data[i, j]
    text_color <- ifelse(value < 0.2, "white", "black")
    grid.text(sprintf("%.2f%%", value * 100), x, y, gp = gpar(fontsize = 8, col = text_color))
  }
)
dev.off()

# visualize amount of cells per cell_type (for responders and non responders)
# for fine annotations
metadata <- seu@meta.data
cell_counts_patients <- metadata %>%
  group_by(SingleR_HCA_fine, patient_id) %>%
  summarise(count = n()) %>%
  ungroup()

# Calculate total counts for each responder status
total_counts_patients <- metadata %>%
  group_by(patient_id) %>%
  summarise(total = n())

# Join total counts with cell counts
cell_counts_patients <- cell_counts_patients %>%
  left_join(total_counts_patients, by = "patient_id")

# Calculate proportions
cell_counts_patients <- cell_counts_patients %>%
  mutate(proportion = count / total)

cell_counts_responders <- metadata %>%
  group_by(SingleR_HCA_fine, responder_status) %>%
  summarise(count = n()) %>%
  ungroup()

# Calculate total counts for each responder status
total_counts_responders <- metadata %>%
  group_by(responder_status) %>%
  summarise(total = n())

# Join total counts with cell counts
cell_counts_responders <- cell_counts_responders %>%
  left_join(total_counts_responders, by = "responder_status")

# Calculate proportions
cell_counts_responders <- cell_counts_responders %>%
  mutate(proportion = count / total)

cell_counts_patients <- cell_counts_patients %>%
  mutate(group = paste0("patient ", patient_id)) %>%
  select(-patient_id)
cell_counts_responders <- cell_counts_responders %>%
  mutate(group = responder_status) %>%
  select(-responder_status) %>%
  mutate(group = ifelse(group == "non_responder", "non responder", group))

head(cell_counts_patients)
head(cell_counts_responders)

# merge dfs
cell_counts <- rbind(cell_counts_patients, cell_counts_responders)

# prepare for heatmap
heatmap_data <- reshape2::dcast(cell_counts, SingleR_HCA_fine ~ group, value.var = "proportion")
rownames(heatmap_data) <- heatmap_data$SingleR_HCA_fine
heatmap_data$SingleR_HCA_fine <- NULL

breaks = c(seq(0, 0.04, length.out = 100))
# breaks = breaks = c(seq(0, 0.01, length.out = 50), seq(0.01, 1, length.out = 50))
colors = viridis(100, option = "viridis")
col_fun = circlize::colorRamp2(breaks, colors)

heatmap_data <- heatmap_data %>%
  select(`patient 0`, `patient 1`, `patient 2`, `patient 3`, responder, `non responder`)

# remove all rows where the rowname includes "Keratinocytes"
heatmap_data_filtered <- heatmap_data[!grepl("Keratinocytes", rownames(heatmap_data)), ]

# heatmap
png(paste0("/path/to/file/"),
    width = 10,
    height = 12,
    units = "in",
    res = 300)
Heatmap(
  as.matrix(heatmap_data_filtered),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  col = col_fun,
  column_title = "Groups",
  name = "Proportion",
  column_title_side = "bottom",
  column_names_rot = 45,
  column_names_centered = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_split = c(rep(1, 4), rep(2, 2)),
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- heatmap_data_filtered[i, j]
    text_color <- ifelse(value < 0.01, "white", "black")
    grid.text(sprintf("%.2f%%", value * 100), x, y, gp = gpar(fontsize = 8, col = text_color))
  },
  row_names_max_width = unit(8, "cm"),
  heatmap_legend_param = list(
    at = c(0, 0.01, 0.02, 0.03, 0.04),
    labels = c("0%", "1%", "2%", "3%", "4%")
  )
)
dev.off()

# heatmap
heatmap_data <- reshape2::dcast(cell_counts, SingleR_HCA ~ patient_id, value.var = "proportion")
rownames(heatmap_data) <- heatmap_data$SingleR_HCA
heatmap_data$SingleR_HCA <- NULL

breaks = breaks = c(seq(0, 0.2, length.out = 50), seq(0.2, 1, length.out = 50))
colors = viridis(100, option = "viridis")
col_fun = circlize::colorRamp2(breaks, colors)

png(paste0("/path/to/file/"),
    width = 7,
    height = 7,
    units = "in",
    res = 300)
Heatmap(
  as.matrix(heatmap_data),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  col = col_fun,
  column_title = "Patient IDs",
  name = "Proportion",
  column_title_side = "bottom",
  column_names_rot = 0,
  row_names_gp = gpar(fontsize = 10), 
  row_names_max_width = unit(4, "cm")
)
dev.off()

heatmap_data_responder_status <- heatmap_data
heatmap_data_patients <- heatmap_data

heatmap_data <- merge(heatmap_data_responder_status, heatmap_data_patients, by = "row.names")
rownames(heatmap_data) <- heatmap_data$Row.names
heatmap_data$Row.names <- NULL

# scatter plot
heatmap_data[is.na(heatmap_data)] <- 0
heatmap_data$Cell_Type <- rownames(heatmap_data)
heatmap_data_long <- melt(heatmap_data, id.vars = "Cell_Type", variable.name = "Column", value.name = "Number")
# Create scatter plot
ggplot(heatmap_data_long, aes(x = reorder(Cell_Type, Number), y = Number, color = Column)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Cell Type", y = "Number", color = "Column")


# make pie chart
heatmap_data[is.na(heatmap_data)] <- 0
heatmap_data$Cell_Type <- rownames(heatmap_data)
heatmap_data_long <- melt(heatmap_data, id.vars = "Cell_Type", variable.name = "Column", value.name = "Number")

png(paste0("/path/to/file/"),
    width = 7,
    height = 7,
    # width = 10,
    # height = 15,
    units = "in",
    res = 300)
ggplot(heatmap_data_long, aes(x="", y=Number, fill=Cell_Type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  facet_wrap(~ Column, scales = "free") +
  theme_void() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  scale_fill_brewer(palette="Dark2")
dev.off()

##############################################################################
# compare Gene Expression Between Responders and Non-Responders per Cell Type #
##############################################################################

# Extract the metadata
metadata <- seu@meta.data

# Get unique cell types
cell_types <- unique(metadata$SingleR_HCA)

# Initialize a list to store results
results_list <- list()

# Loop through each cell type
for (cell_type in cell_types) {
  # Create a subset of the Seurat object for the current cell type
  subset_seurat <- subset(seu, subset = SingleR_HCA == cell_type)
  
  Idents(subset_seurat) <- "responder_status"
  
  # Perform differential expression analysis between responders and non-responders
  markers <- FindMarkers(subset_seurat, 
                         ident.1 = "responder", 
                         ident.2 = "non_responder", 
                         logfc.threshold = 0.25) # Adjust logfc.threshold as needed
  
  # Add cell type information to the results
  markers$cell_type <- cell_type
  
  # Store results in the list
  results_list[[cell_type]] <- markers
}

# Combine results into a single data frame
combined_results <- bind_rows(results_list, .id = "cell_type")

top_markers <- combined_results %>%
  top_n(20, avg_log2FC) # Select top 20 markers

ggplot(filter(combined_results, cell_type == "T_cells"), 
       aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Volcano Plot for Specific Cell Type",
       x = "Log Fold Change",
       y = "-Log10 Adjusted P-value")

combined_results <- combined_results %>%
  mutate(significant = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 1, "Significant", "Not Significant"))
combined_results <- combined_results %>%
  mutate(regulation = case_when(
    p_val_adj < 0.05 & avg_log2FC > 1 ~ "Upregulated",
    p_val_adj < 0.05 & avg_log2FC < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
))

# remove rows which start with "Blank"
rows_to_remove <- grep("^Blank", rownames(combined_results), value = TRUE)
combined_results <- combined_results[!rownames(combined_results) %in% rows_to_remove, ]

top_genes <- combined_results %>%
  filter(p_val_adj < 0.05) %>%  # Ensure only significant genes are considered
  arrange(desc(abs(avg_log2FC))) %>%  # Sort by absolute fold change
  slice_head(n = 30)  # Select the top 10

gene_names <- rownames(top_genes)
cleaned_gene_names <- gsub("\\.\\.\\.\\d+$", "", gene_names)
top_genes$gene <- cleaned_gene_names

png(paste0("/path/to/file/"),
    width = 13,
    height = 6,
    units = "in",
    res = 300)
ggplot(combined_results, aes(x = avg_log2FC, y = pmin(-log10(p_val_adj), 100), color = regulation)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +  # Add a significance threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +  # Add fold change threshold lines
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +  # Color based on regulation
  facet_wrap(~ cell_type) +
  geom_text_repel(data = top_genes, aes(label = gene), color = "black", size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 100) +  # Add text labels
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes Across Cell Types",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "right")
dev.off()

# big dotplot
gene_information <- read.csv("/path/to/file/", header = TRUE, sep = ";")

cell_types <- c("t_cells", "trm_cells", "b_cells", "keratinocytes", "nk_cells", "macrophages", "dendritic_cells")
gene_information_cells <- gene_information[gene_information$type %in% cell_types,]
gene_information_cells$features.plot <- gene_information_cells$gene
gene_information_cells$label <- gene_information_cells$type
data.anno <- gene_information_cells[,c(4,5)]

data.usage <- DotPlot(seu, features = gene_information_cells$gene, group.by = "SingleR_HCA")$data

df.plot <- plyr::join(data.usage, data.anno)

df.plot$id <- factor(df.plot$id, levels = sort(levels(df.plot$id), decreasing = T))

# rename cell types
df.plot$label <- gsub("t_cells", "t-cells", df.plot$label)
df.plot$label <- gsub("trm_cells", "trm-cells", df.plot$label)
df.plot$label <- gsub("b_cells", "b-cells", df.plot$label)
df.plot$label <- gsub("nk_cells", "nk-cells", df.plot$label)
df.plot$label <- gsub("dendritic_cells", "dendritic cells", df.plot$label)

df.plot$id <- gsub("T_cells", "t-cells", df.plot$id)
df.plot$id <- gsub("NK_cell", "nk-cells", df.plot$id)
df.plot$id <- gsub("Macrophage", "macrophages", df.plot$id)
df.plot$id <- gsub("Keratinocytes", "keratinocytes", df.plot$id)
df.plot$id <- gsub("DC", "dendritic cells", df.plot$id)
df.plot$id <- gsub("B_cell", "b-cells", df.plot$id)

# convert back to factor
df.plot$id <- as.factor(df.plot$id)

pdf(paste0("/path/to/file/"),
    width = 18,
    height = 5)
ggplot(df.plot,aes(x=features.plot,y =  as.numeric(id),size = pct.exp, color = avg.exp.scaled))+
  geom_point() + 
  scale_size("% detected", range = c(0.5,6)) +
  scale_color_gradientn(colours = viridis::viridis(20),
                        guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                        name = "Average\nexpression") +
  cowplot::theme_cowplot() + 
  ylab("Annotated Cell types") + xlab("Markers") + theme_bw() +
  scale_y_continuous(breaks = 1:length(levels(df.plot$id)),labels = levels(df.plot$id))+
  facet_grid(~label, scales="free_x",space = "free") + 
  theme_classic() +
  theme(
    axis.text.x = element_text(size=10, angle=90, hjust=1, vjust = 0.5, color="black", margin = margin(t=10)),#face="bold"),
    axis.text.y = element_text(size=12, color="black"),
    axis.title.x = element_text(size=14,colour = 'black',vjust = -0.8,hjust = 0.5),
    axis.title.y = element_text(size=14,colour = 'black',vjust = 1.4,hjust = 0.5),
    
    axis.ticks.y = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_blank(),
    
    panel.spacing=unit(5, "mm"),
    strip.text.x = element_text(size=12, color = "black",
                                vjust = 0.5,margin = margin(b = 10,t=3)),
    strip.background = element_rect(colour="white", fill="white",size = 1),
    plot.title = element_text(hjust = 0.5)
    
  ) +
  labs(title = "")
dev.off()

top_markers <- filter(combined_results, cell_type == "T_cells") %>%
  top_n(20, avg_log2FC) # Select top 20 markers

# Plot heatmap
DoHeatmap(seu, features = top_markers$gene)

# make nice plot
gene_names <- rownames(combined_results)
cleaned_gene_names <- gsub("\\.\\.\\.\\d+$", "", gene_names)
combined_results$gene <- cleaned_gene_names

top_genes <- combined_results %>%
  group_by(cell_type) %>%
  top_n(10, avg_log2FC) # Select top 10 genes per cell type; adjust as needed

# Get top gene names
top_genes_list <- top_genes$gene

unique_top_genes_list <- unique(top_genes_list)
top_genes_seu <- subset(seu, features = unique_top_genes_list)

# Extract the expression data for top genes
expression_data <- GetAssayData(top_genes_seu, slot = "data")

valid_genes <- unique_top_genes_list[unique_top_genes_list %in% rownames(expression_data)]

expression_data_subset <- expression_data[valid_genes, ]

expression_data_dense <- as.matrix(expression_data_subset)

pheatmap(expression_data_dense, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         show_rownames = TRUE, 
         show_colnames = TRUE,
         annotation_col = top_genes_seu@meta.data[, c("SingleR_HCA", "responder_status")],
         main = "Heatmap of Top Differentially Expressed Genes")

# Create a heatmap using pheatmap
pheatmap(expression_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = TRUE,
         main = "Heatmap of Top Differentially Expressed Genes",
         annotation_col = top_genes_seu@meta.data[, c("SingleR_HCA", "responder_status")])









