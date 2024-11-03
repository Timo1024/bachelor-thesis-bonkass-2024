# ==============================================================================
# File: batch_correction.R
# Description: This script performs batch correction, centering and
# normalization, filtering and doublet removal
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
library(remotes)
library(RColorBrewer)
library(scales)
library(fields)
library(reshape2)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(scDblFinder)
library(patchwork)

# Main  ========================================================================

# add new samples here
patient_list <- list(
  "1" = c("responder", "/path/to/file/"),
  "2" = c("responder", "/path/to/file/"),
  "3" = c("responder", "/path/to/file/"),
  "4" = c("responder", "/path/to/file/"),
  "5" = c("responder", "/path/to/file/"),
  "6" = c("non_responder", "/path/to/file/"),
  "7" = c("non_responder", "/path/to/file/"),
  "8" = c("non_responder", "/path/to/file/"),
  "9" = c("non_responder", "/path/to/file/"),
  "10" = c("non_responder", "/path/to/file/")
)

# Function to process and create Seurat object for each patient
process_patient <- function(patient_id, responder_status, path, min.features = 6) {
  file_path <- path
  cell_by_gene <- read.csv(file_path, numerals = c("no.loss"))
  
  rownames(cell_by_gene) <- cell_by_gene$cell
  cell_by_gene <- cell_by_gene[,-1]
  
  cell_by_gene_matrix <- as.matrix(cell_by_gene)
  cell_by_gene_matrix <- t(cell_by_gene_matrix)
  
  rownames(cell_by_gene_matrix) <- colnames(cell_by_gene)
  colnames(cell_by_gene_matrix) <- rownames(cell_by_gene)
  
  cell_by_gene_matrix_converted <- as(cell_by_gene_matrix, "dgCMatrix")
  
  seu_spacial <- CreateSeuratObject(counts = cell_by_gene_matrix_converted, project = paste0("spacial_region_", patient_id), min.cells = 0, min.features = min.features)
  
  seu_spacial <- subset(seu_spacial, features = rownames(seu_spacial)[!grepl("^Blank\\.\\d+$", rownames(seu_spacial))])
  
  print(length(rownames(seu_spacial)))
  
  # Add metadata
  seu_spacial$patient_id <- patient_id
  seu_spacial$responder_status <- responder_status
  
  return(seu_spacial)
}

seu_spacial_list <- list()
for (patient in names(patient_list)) {
  # Extract the response status and file path
  response_status <- patient_list[[patient]][1]
  file_path <- patient_list[[patient]][2]
  
  seu_ <- process_patient(patient_id = patient, responder_status = response_status, path = file_path)
  
  seu_ <- NormalizeData(seu_)
  seu_ <- FindVariableFeatures(seu_, selection.method = "vst", nfeatures = 2000)
  
  seu_spacial_list <- append(seu_spacial_list, seu_)
  
}

anchors <- FindIntegrationAnchors(object.list = seu_spacial_list)

integrated <- IntegrateData(anchorset = anchors)

integrated <- ScaleData(integrated, features = rownames(integrated))
integrated <- RunPCA(integrated, features = VariableFeatures(object = integrated))
integrated <- FindNeighbors(integrated, dims = 1:12)
integrated <- FindClusters(integrated, resolution = 0.5)
integrated <- RunUMAP(integrated, dims = 1:12)

integrated_joined <- JoinLayers(integrated, assay = "RNA")

# Function to detect and remove doublets from a Seurat object
detect_and_remove_doublets <- function(seurat_obj) {
  counts_data <- seurat_obj@assays$RNA@layers$counts
  sce_obj <- SingleCellExperiment(assays = list(counts = counts_data))
  
  sce_obj <- scDblFinder(sce_obj)
  
  # Add doublet information back to Seurat object
  seurat_obj$doublet_score <- sce_obj$scDblFinder.score
  seurat_obj$doublet_class <- sce_obj$scDblFinder.class
  
  # Filter out doublets
  seurat_obj <- subset(seurat_obj, subset = doublet_class == "singlet")
  
  return(seurat_obj)
}

# Detect and remove doublets for each Seurat object
integrated_joined <- detect_and_remove_doublets(integrated_joined)

# get top variable features
# Split the Seurat object by patient_id
seurat_list <- SplitObject(integrated_joined, split.by = "patient_id")

# Initialize an empty list to store genes present in less than 10 cells for each patient
gene_counts_list <- list()

# Loop through each patient and count the number of cells expressing each gene
for (patient in names(seurat_list)) {
  # Get the counts matrix for the patient
  counts_matrix <- seurat_list[[patient]]@assays$RNA@layers$counts

  # Calculate the number of cells in which each gene is expressed
  gene_counts <- rowSums(counts_matrix > 0)
  names(gene_counts) <- rownames(seurat_list[[patient]])

  gene_counts_list[[patient]] <- gene_counts
}

# Combine the gene counts into a matrix
gene_counts_matrix <- do.call(cbind, gene_counts_list)

# Calculate the mean for each gene
gene_means <- rowMeans(gene_counts_matrix)

# Set a threshold for the mean count
threshold <- 5

# Get the list of genes with a mean count below the threshold
genes_above_threshold <- names(gene_means[gene_means >= threshold])

# Print the list of genes
print(genes_above_threshold)

integrated_joined_filtered <- subset(integrated_joined, features = genes_above_threshold)

SaveSeuratRds(integrated_joined_filtered, file = "/path/to/file/")
