# ==============================================================================
# File: prepare_for_misty.R
# Description: prepare data so that it can be used in python as input for MISTy
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
library(Seurat)
library(SeuratHelper)
library(anndata)
library(dplyr)

# Main =========================================================================

seu <- LoadSeuratRds("/path/to/file/")

patients <- unique(seu$patient_id)

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

for (patient in patients) {
  seu_subset <- subset(seu, subset = patient_id == patient)
  
  gene_expression_matrix <- as.matrix(seu_subset@assays$RNA$data)
  gene_expression_matrix_t <- t(gene_expression_matrix)
  
  cell_metadata <- seu_subset@meta.data
  
  gene_metadata <- data.frame(gene_name = rownames(seu_subset@assays$RNA$data))
  
  metadata_sample <- read.csv(paste0(patient_list[[patient]][2], "cell_metadata.csv"), numerals = c("no.loss"))
  
  cell_metadata$id <- rownames(cell_metadata)
  metadata_sample$id <- metadata_sample$EntityID
  
  # merge metadata dfs
  cell_metadata <- merge(cell_metadata, metadata_sample, by = "id")
  rownames(cell_metadata) <- cell_metadata$id
  
  cell_metadata_for_adata <- cell_metadata %>%
    select(id, center_x, center_y, seurat_clusters, SingleR_HCA, DAPI_raw, DAPI_high_pass, PolyT_raw, PolyT_high_pass, solidity, anisotropy, transcript_count, perimeter_area_ratio, volume)
  
  write.csv(cell_metadata_for_adata, paste0("/path/to/file/", patient, "/path/to/file/"))
  write.csv(gene_expression_matrix_t, paste0("/path/to/file/", patient, "/path/to/file/"), row.names = TRUE)
  
}
