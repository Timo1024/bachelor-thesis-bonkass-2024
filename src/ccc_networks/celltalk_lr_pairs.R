# ==============================================================================
# File: celltalk_lr_pairs.R
# Description: getting all ligand-receptor pairs which can be found since they
# need to be included in the 300 gene panel.
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
library(celltalker)
library(igraph)
library(ggraph)
library(dplyr)
library(tidyr)
library(reshape2)
library(ComplexHeatmap)

# Main =========================================================================

seu <- LoadSeuratRds("/path/to/file/")

possible_interactions <- ramilowski_pairs[ramilowski_pairs$ligand %in% rownames(seu) & ramilowski_pairs$receptor %in% rownames(seu),]

print(possible_interactions$pair)















