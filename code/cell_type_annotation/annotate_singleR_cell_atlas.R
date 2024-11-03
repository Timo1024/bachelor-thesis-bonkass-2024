# ==============================================================================
# File: annotate_singleR_cell_atlas.R
# Description: performing reference annotation using the human primary cell
# atlas as reference when applying SingleR
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

# Main  ========================================================================

# load seurat object
seu <- LoadSeuratRds("/path/to/file/")

# load primary cell atlas
ref <- celldex::HumanPrimaryCellAtlasData()
ref <- ref[,grepl('Keratinocytes|^B_cell|T_cells|DC|Macrophage|NK_cell', ref$label.main)]

# make singleR annotations
norm_counts <- LayerData(seu, assay = "RNA", layer = 'data')

ct_ann <- SingleR(test = norm_counts,
                  ref = ref,
                  # labels = ref$label.fine,
                  labels = ref$label.main,
                  de.method = 'wilcox')

# make plot of annotation quality
png(paste0("/path/to/file/"),
    width = 7,
    height = 8,
    units = "in",
    res = 300)
print(plotScoreHeatmap(ct_ann))
dev.off()

png(paste0("/path/to/file/"),
    width = 7,
    height = 8,
    units = "in",
    res = 300)
print(plotDeltaDistribution(ct_ann, ncol = 3, dots.on.top = FALSE))
dev.off()

# Add to seurat object
seu <- AddMetaData(seu, ct_ann$pruned.labels, col.name = 'SingleR_HCA')

seu <- SetIdent(seu, value = "SingleR_HCA")
png(paste0("/path/to/file/"),
    width = 8,
    height = 7,
    units = "in",
    res = 300)
print(DimPlot(seu, label = T , repel = T, label.size = 3) + NoLegend())
dev.off()

pdf(file = "/path/to/file/",
    width = 10,
    height = 7)
Idents(seu) <- "SingleR_HCA"
umap_plot <- DimPlot(seu, reduction = "umap")
umap_plot + labs(x = "UMAP 1", y = "UMAP 2")
dev.off()

# save seurat object
saveRDS(seu, file = paste0("/path/to/file/"))








