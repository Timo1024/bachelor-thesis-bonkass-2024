# ==============================================================================
# File: cell_counts_and_pca_qc_check.R
# Description: this script performs quality control by inspecting the 
# amount of cells and transcripts and looking at the misidentification control
# words. Additionally the quality of the PCA is assessed by looking at the 
# amount of variance explained by each principal component
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

seu <- LoadSeuratRds("/path/to/file/")

sum(seu$nCount_RNA)
length(seu$nCount_RNA)

patients <- unique(seu$patient_id)
patients <- as.character(patients)
seurat_list <- list()

for (patient_ in patients) {
  patient_ <- as.character(patient_)
  print(patient_)
  seurat_list[[patient_]] <- subset(seu, subset = patient_id == patient_)
}

metadata <- list()
for (patient_ in patients) {
  metadata[[patient]] <- read.csv(paste0(patient_list[[patient]][2], "filename.csv"))
}

transcripts_list <- list()
for (patient_ in patients) {
  transcripts_list[[patient]] <- read.csv(paste0(patient_list[[patient]][2], "filename.csv"))
}

amount_transcripts <- 0
amount_negative_ones <- 0
for (patient in patients) {
  amount_transcripts <- amount_transcripts + nrow(transcripts_list[[patient]])
  amount_negative_ones <- amount_negative_ones + nrow(transcripts_list[[patient]][transcripts_list[[patient]]$cell_id == -1,])
}

print(amount_negative_ones/amount_transcripts)

cell_by_gene <- list()
for (patient_ in patients) {
  cell_by_gene[[patient]] <- read.csv(paste0(patient_list[[patient]][2], "filename.csv"))
}

cell_by_gene_all <- do.call(rbind, cell_by_gene)
cell_by_gene_all$sum <- rowSums(cell_by_gene_all %>% select(-cell))
print(mean(cell_by_gene_all$sum))

# amount of cells
cell_amount_list <- list()

for (patient in patients) {
  cell_amount_list[[patient]] <- c(
    nrow(metadata[[patient]]),
    length(seurat_list[[patient]]$nCount_RNA)
  )
}

cell_amount_df <- as.data.frame(cell_amount_list)
colnames(cell_amount_df) <- patients
rownames(cell_amount_df) <- c("before filtering", "after filtering")

cell_amount_df_sum <- cell_amount_df
cell_amount_df_sum$sum <- rowSums(cell_amount_df_sum)

cell_amount_df_thesis <- as.data.frame(t(cell_amount_df)) %>%
  mutate(percentage = (`after filtering` / `before filtering`) * 100)

cell_amount_df_long <- cell_amount_df %>%
  rownames_to_column(var = "Status") %>%
  pivot_longer(cols = -Status, names_to = "Patient", values_to = "Value") %>%
  mutate(Patient = paste("Patient", Patient))

# amount of transcripts
transcripts_amount_list <- list()

for (patient in patients) {
  transcripts_amount_list[[patient]] <- c(
    nrow(transcripts_list[[patient]]),
    sum(seurat_list[[patient]]$nCount_RNA)
  )
}

transcripts_amount_df <- as.data.frame(transcripts_amount_list)
colnames(transcripts_amount_df) <- patients
rownames(transcripts_amount_df) <- c("before filtering", "after filtering")

transcripts_amount_df_thesis <- as.data.frame(t(transcripts_amount_df)) %>%
  mutate(percentage = (`after filtering` / `before filtering`) * 100)

transcripts_amount_df_long <- transcripts_amount_df %>%
  rownames_to_column(var = "Status") %>%
  pivot_longer(cols = -Status, names_to = "Patient", values_to = "Value") %>%
  mutate(Patient = paste("Patient", Patient))

# plots
pdf(file = "/path/to/file/",
    width = 5,
    height = 5)
cell_amount_df_long$Status <- factor(cell_amount_df_long$Status, levels = c("before filtering", "after filtering"))
ggplot(cell_amount_df_long, aes(x = Patient, y = Value, fill = Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Patient",
       y = "Amount of cells") +
  scale_fill_manual(values = c("before filtering" = "#00bfc4", "after filtering" = "#f8766d"),
                    labels = c("Before \nfiltering", "After \nfiltering")) +
  theme_minimal()+
  theme(legend.spacing = unit(2, "cm"),  # Adjust the spacing between legend items
        legend.key.size = unit(1.5, "lines"),
        legend.title = element_blank())
dev.off()

pdf(file = "/path/to/file/",
    width = 5,
    height = 5)
transcripts_amount_df_long$Status <- factor(transcripts_amount_df_long$Status, levels = c("before filtering", "after filtering"))
ggplot(transcripts_amount_df_long, aes(x = Patient, y = Value, fill = Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Patient",
       y = "Amount of transcripts") +
  scale_fill_manual(values = c("before filtering" = "#00bfc4", "after filtering" = "#f8766d"),
                    labels = c("Before \nfiltering", "After \nfiltering")) +
  theme_minimal()+
  theme(legend.spacing = unit(2, "cm"),  # Adjust the spacing between legend items
        legend.key.size = unit(1.5, "lines"),
        legend.title = element_blank())
dev.off()


# quality control using 15 misidentification control words 

cell_by_gene_all_sums <- cell_by_gene_all
cell_by_gene_all_sums$negative_sums <- rowSums(cell_by_gene_all_sums[,302:316])

sum(cell_by_gene_all_sums$negative_sums) / sum(cell_by_gene_all_sums$sum)

# check quality of PCA
stdev <- seu[["pca"]]@stdev
variance_explained <- stdev^2 / sum(stdev^2)

variance_df <- data.frame(PC = 1:length(variance_explained), VarianceExplained = variance_explained)
print(variance_df)

pdf(file = "/path/to/file/",
    width = 10,
    height = 5)
ggplot(variance_df, aes(x = PC, y = VarianceExplained)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 12, linetype = "dashed", col = "grey") +
  labs(x = "Principal Components", y = "Variance Explained") +
  theme_minimal()
dev.off()

# cumulative sum of 12th PC
cumulative_variance <- cumsum(variance_df$VarianceExplained)
print(cumulative_variance[12])





