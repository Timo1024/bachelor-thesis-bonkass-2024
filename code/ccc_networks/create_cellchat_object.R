# ==============================================================================
# File: create_cellchat_object.R
# Description: creating the CellChat object to infer CCCs
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

metadata <- list()
for (patient in patients) {
  metadata[[patient]] <- read.csv(paste0(patient_list[[patient]][2], "/path/to/file/"), numerals = c("no.loss"))
}

# split seurat object into the four patients
seurat_list <- list()
for (patient in patients) {
  seurat_list[[patient]] <- subset(seu, subset = patient_id == patient)
}

cellchat_list <- list()
for (patient in patients) {
  
  # Prepare input data for CellChat analysis
  data.input <- GetAssayData(seurat_list[[patient]], layer = "data", assay = "RNA")
  
  Idents(seurat_list[[patient]]) <- "SingleR_HCA"
  meta = data.frame(labels = Idents(seurat_list[[patient]]), row.names = names(Idents(seurat_list[[patient]])))
  
  # get location of all cell ids
  rownames(metadata[[patient]]) <- metadata[[patient]]$EntityID
  metadata[[patient]] <- metadata[[patient]] %>%
    mutate(x = as.numeric(center_x),
           y = as.numeric(center_y)) %>%
    select(x, y)
  # get just the rows which are also in seu
  metadata_subset <- metadata[[patient]][rownames(metadata[[patient]]) %in% names(seurat_list[[patient]]@active.ident), ]
  
  spatial.locs <- as.matrix(metadata_subset)
  
  # create scale.factors df
  ratio <- 1 # as the x and y are already microns
  
  cell_distances <- computeCellDistance(spatial.locs)
  cell_distances_vector <- as.numeric(cell_distances)
  min_distance <- min(cell_distances_vector, na.rm = TRUE)
  tol <- min_distance / 2
  
  scale.factors <- data.frame(ratio = ratio, tol = tol)
  
  # create cellchat object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                             datatype = "spatial", coordinates = spatial.locs, spatial.factors = scale.factors,
                             do.sparse = TRUE)

  # Set the ligand-receptor interaction database
  CellChatDB <- CellChatDB.human
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  
  # preprocessing
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat)
  future::plan("multisession", workers = 4) 
  cellchat <- identifyOverExpressedGenes(cellchat, thresh.p = 0.01, min.cells = 2)
  cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = FALSE)
  
  # communication probability and infer cellular communication network
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, # trim = 0.1,
                                distance.use = TRUE, interaction.range = 1000, scale.distance = 1, # interaction.range = 250, scale.distance = 0.01, 
                                contact.range = 1000) # contact.range = 10
  
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10) # min.cells = 10)
  
  # cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat) #, thresh = 0)
  
  # Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  
  cellchat_list[[patient]] <- cellchat
}

# save cellchat objects
for (patient in patients) {
  # Construct the filename
  filename <- paste0("/path/to/file/", patient, "/path/to/file/")

  # Save the cellchat object to a file
  saveRDS(cellchat_list[[patient]], file = filename)
}

