# ==============================================================================
# File: check_cluster_quality.R
# Description: this script checks the quality of the clustering by calculating
# the silhouette widths
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

ptm = Sys.time()

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
library(patchwork)
library(circlize)

# Main  ========================================================================

seu <- LoadSeuratRds("/path/to/file/")

Idents(seu) <- "seurat_clusters"

resolutions <- seq(0.5, 1.5, by = 0.1)

silhouette_scores <- c()

for (resolution in resolutions) {
  
  seu <- FindClusters(seu, resolution = resolution)

  # Convert factor to numeric
  ident_numeric <- as.numeric(as.factor(Idents(seu)))

  # Calculate distance matrix
  dist_matrix <- dist(seu@reductions$pca@cell.embeddings)

  # Calculate silhouette scores
  silhouette_scores <- silhouette(ident_numeric, dist_matrix)

  # Convert silhouette scores to a data frame
  silhouette_df <- data.frame(
    cluster = silhouette_scores[, 1],
    neighbor = silhouette_scores[, 2],
    silhouette_width = silhouette_scores[, 3]
  )

  write.csv(silhouette_df, paste0("/path/to/", resolution, "/file/"), row.names = FALSE)
  # silhouette_df <- as.data.frame(read_csv(paste0("/path/to/", resolution, "/file/"), show_col_types = FALSE))

  silhouette_df <- silhouette_df %>%
    arrange(cluster, desc(silhouette_width))
  silhouette_df$cell <- seq_along(silhouette_df$silhouette_width)

  combined_silhouette_score <- mean(silhouette_df$silhouette_width)
  
  print(resolution)
  print(combined_silhouette_score)
  print("")
  
  silhouette_scores <- c(silhouette_scores, combined_silhouette_score)

  pdf(paste0("/path/to/", resolution, "/file/"),
      width = 13,
      height = 6)
  plot(ggplot(silhouette_df) +
    geom_col(aes(x = cell, y = silhouette_width, group = cluster, fill = factor(cluster)), position = position_dodge(width = 0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = viridis::viridis(length(unique(silhouette_df$cluster)))) +
    xlab("Cells") +
    ylab("Silhouette Width") +
    labs(title = "", fill = "Clusters") +
    theme_minimal())
  dev.off()
  
}

print(max(silhouette_scores))





