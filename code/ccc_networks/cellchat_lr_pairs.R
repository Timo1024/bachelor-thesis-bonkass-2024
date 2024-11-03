# ==============================================================================
# File: cellchat_lr_pairs.R
# Description: getting all ligand-receptor pairs which can be found since they
# need to be included in the 300 gene panel.
# Additionally looking at S100A8/A9 - TLR4 interaction
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

# Main =========================================================================

seu <- LoadSeuratRds("/path/to/file/")

# cellchat
CellChatDB <- CellChatDB.human

possible_interactions <- CellChatDB$interaction[CellChatDB$interaction$ligand %in% rownames(seu) & CellChatDB$interaction$receptor %in% rownames(seu),]

CellChatDB$interaction[
  CellChatDB$interaction$ligand == "S100A8" | 
    CellChatDB$interaction$receptor == "S100A8", ]

# celltalker
data("ramilowski_pairs")

s100a8_tlr4_pairs <- ramilowski_pairs[
  ramilowski_pairs$ligand == "S100A9" & ramilowski_pairs$receptor == "TLR4", ]

# Display the result
print(s100a8_tlr4_pairs)

# load transcript data
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
patients <- unique(seu$patient_id)

transcripts_list <- list()
for (patient in patients) {
  transcripts_list[[patient]] <- read.csv(paste0(patient_list[[patient]][2], "/path/to/file/"), numerals = c("no.loss"))
}

df <- transcripts_list[["6"]]

# Filter the data frame to include only rows with the gene "S100A8"
s100a8_df <- df[df$gene %in% c("S100A8", "S100A9", "TLR4"), ]

# Plot the locations of "S100A8" genes
ggplot(s100a8_df, aes(x = as.numeric(global_x), y = as.numeric(global_y), color = gene)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  labs(title = "Spatial Location of S100A8 Genes",
       x = "X Coordinate",
       y = "Y Coordinate")

# get normalized S100A8/A9 and TLR4 counts
# Initialize an empty data frame to store the results
results_df <- data.frame(
  sample = character(),
  S100A8 = numeric(),
  S100A9 = numeric(),
  TLR4 = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each sample
for (i in 1:length(transcripts_list)) {
  # Get the current sample data frame
  sample_df <- transcripts_list[[i]]
  
  # Calculate the total number of rows in the sample
  total_rows <- nrow(sample_df)
  
  # Calculate the normalized counts for each gene
  s100a8_count <- sum(sample_df$gene == "S100A8") / total_rows
  s100a9_count <- sum(sample_df$gene == "S100A9") / total_rows
  tlr4_count <- sum(sample_df$gene == "TLR4") / total_rows
  
  # Add the results to the results data frame
  results_df <- rbind(results_df, data.frame(
    sample = paste("Sample", i),
    S100A8 = s100a8_count,
    S100A9 = s100a9_count,
    TLR4 = tlr4_count
  ))
}

# Display the results
print(results_df)

ggplot(results_df, aes(x = sample, y = TLR4)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_minimal() +
  labs(title = "Normalized S100A8 Counts Across Samples",
       x = "Sample",
       y = "Normalized S100A8 Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

