# ==============================================================================
# File: analyze_interactions_misty.R
# Description: Compare results of MISTy CCC analysis and plot them
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
library(metap)
library(RColorBrewer)

# Main =========================================================================

# load seurat object
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

# usable_patients
usable_patients <- c("1", "2", "7", "8")

# lr interactions
# load the files for all patients
misty_output_list <- list()
for (patient in usable_patients) {
  misty_output_list[[patient]] <- read.csv(paste0("/path/to/file/", patient, "/path/to/file/"))
}

# get unique genes
unique_genes <- c()
for (patient in usable_patients) {
  unique_genes <- unique(c(unique_genes, misty_output_list[[patient]]$target, misty_output_list[[patient]]$predictor))
}

expanded_dfs <- list()
for (patient in usable_patients) {
  current_df <- misty_output_list[[patient]]
  
  all_pairs <- expand.grid(target = unique_genes, predictor = unique_genes)
  
  all_pairs$view <- "extra"
  all_pairs$importances <- NA
  
  merged_df <- merge(all_pairs, current_df, 
                     by = c("target", "predictor", "view"), 
                     all.x = TRUE)
  
  expanded_dfs[[patient]] <- merged_df
  
  expanded_dfs[[patient]] <- expanded_dfs[[patient]][, c("target", "predictor", "view", "importances.y")]
  
  # Rename the 'importances.y' column to 'importances'
  colnames(expanded_dfs[[patient]])[colnames(expanded_dfs[[patient]]) == "importances.y"] <- "importances"
}

# calculate difference nd p-value
# Combine the data frames for responders and non-responders
responders <- rbind(expanded_dfs[["1"]], expanded_dfs[["2"]])
non_responders <- rbind(expanded_dfs[["7"]], expanded_dfs[["8"]])

# Calculate the mean importance for each gene-gene pair in responders and non-responders
mean_importance_responders <- aggregate(importances ~ target + predictor, data = responders, FUN = mean, na.rm = TRUE)
mean_importance_non_responders <- aggregate(importances ~ target + predictor, data = non_responders, FUN = mean, na.rm = TRUE)

# Merge the mean importance data frames
mean_importance <- merge(mean_importance_responders, mean_importance_non_responders, by = c("target", "predictor"), suffixes = c("_responders", "_non_responders"))

# Perform a t-test to calculate the p-value for the difference in means
mean_importance$p_value <- apply(mean_importance, 1, function(row) {
  target <- row["target"]
  predictor <- row["predictor"]
  
  # Extract the importances for the current gene-gene pair
  importances_responders <- responders$importances[responders$target == target & responders$predictor == predictor]
  importances_non_responders <- non_responders$importances[non_responders$target == target & non_responders$predictor == predictor]
  
  # Check if there are enough observations for the t-test
  if (sum(!is.na(importances_responders)) < 2 || sum(!is.na(importances_non_responders)) < 2) {
    return(1)  # Set p-value to 1 (maximal) if not enough non-NA observations
  } else {
    # Perform the t-test
    t_test <- t.test(importances_responders, importances_non_responders, na.rm = TRUE)
    return(t_test$p.value)
  }
})

mean_importance$difference <- mean_importance$importances_responders - mean_importance$importances_non_responders

# View the result
print(mean_importance)

# dotplot lr pairs
dot_plot <- ggplot(mean_importance, aes(x = target, y = predictor, color = difference, size = -log10(p_value))) +
  geom_point() +
  scale_color_viridis(option = "D") +
  scale_size_continuous(range = c(1, 10)) +
  theme_minimal() +
  labs(title = "", x = "Target (Receptor)", y = "Predictor (Ligand)", size = "-log10(p-value)", color = "Change of \naverage importance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(dot_plot)

# Print the plot
pdf(paste0("/path/to/file/"),
    width = 9,
    height = 11)
print(dot_plot)
dev.off()










