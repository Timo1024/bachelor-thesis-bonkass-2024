# ==============================================================================
# File: manual_annotations_comparison_singleR_HCA.R
# Description: comparison of the manual and reference annotation approach.
# A new meta data column is created which included the annotations which are
# the same for both methods.
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
library(ggalluvial)
library(caret)

# Main =========================================================================

seu <- LoadSeuratRds("/path/to/file/")

alluvial_data <- data.frame(
  SingleR_HCA = seu$SingleR_HCA,
  Manual = seu$manual_annotations
)

ggplot(alluvial_data,
       aes(axis1 = SingleR_HCA, axis2 = Manual)) +
  geom_alluvium(aes(fill = SingleR_HCA)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  labs(title = "Comparison of Annotations",
       x = "Annotations",
       y = "Count")

# confusion matrix
conf_matrix <- table(seu$SingleR_HCA, seu$manual_annotations)
conf_matrix <- as.matrix(conf_matrix)

# scale and log transform
log_transformed_matrix <- log1p(conf_matrix)
scaled_matrix <- scale(log_transformed_matrix, center = TRUE, scale = TRUE)

df <- as.data.frame(as.table(scaled_matrix))
colnames(df) <- c("SingleR", "Manual", "Value")
df$highlight <- FALSE
df[c(3, 10, 18),]$highlight <- TRUE

pdf(file = "/path/to/file/",
    width = 5,
    height = 4)
ggplot(df, aes(x = Manual, y = SingleR, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_tile(data = df[df$highlight == TRUE, ], aes(x = Manual, y = SingleR), color = "black", size = 1.5, fill = NA) +
  theme_minimal() +
  labs(title = "",
       x = "Manual",
       y = "SingleR",
       fill = "Scaled \nValue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

# new column for combined_annotations
seu$manual_annotations[seu$manual_annotations == "T-cells"] <- "T_cells"
seu$SingleR_HCA[seu$SingleR_HCA == "Macrophage"] <- "Macrophages"
seu$combined_annotations <- ifelse(seu$manual_annotations == seu$SingleR_HCA, seu$manual_annotations, NA)

SaveSeuratRds(seu, "/path/to/file/")

# create plots for how much of which cell type is present in which patient
annotations_df <- data.frame(patient = seu$patient_id, annotation = seu$combined_annotations)
annotations_df <- na.omit(annotations_df)

# Count the number of cells for each cell type per patient
cell_counts <- annotations_df %>%
  group_by(patient, annotation) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the total number of cells per patient
total_counts <- annotations_df %>%
  group_by(patient) %>%
  summarise(total = n(), .groups = 'drop')

normalized_counts <- cell_counts %>%
  left_join(total_counts, by = "patient") %>%
  mutate(normalized_count = count / total)

normalized_counts <- normalized_counts %>%
  complete(patient, annotation, fill = list(normalized_count = 0))

normalized_counts$patient <- factor(paste("Sample", normalized_counts$patient), levels = paste("Sample", as.character(1:10)))

normalized_counts <- normalized_counts %>%
  mutate(responder_status = ifelse(patient %in% paste("Sample", 1:5), "Responder", "Non-responder"))

# Plot the counts, splitting the plot between patients 1-5 and 6-10
pdf(file = "/path/to/file/",
    width = 10,
    height = 7)
ggplot(normalized_counts, aes(x = annotation, y = normalized_count, fill = annotation)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ patient, scales = "free_x", ncol = 5) +
  theme_minimal() +
  labs(title = "",
       x = "Cell type",
       y = "Relative counts",
       fill = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
dev.off()

# p-values for different proportions
# Separate the data for responders and non-responders
responders <- normalized_counts %>% filter(responder_status == "Responder")
non_responders <- normalized_counts %>% filter(responder_status == "Non-responder")

# Perform the Wilcoxon rank-sum test for each cell type
keratinocytes_p_value <- wilcox.test(responders %>% filter(annotation == "Keratinocytes") %>% pull(normalized_count),
                                     non_responders %>% filter(annotation == "Keratinocytes") %>% pull(normalized_count))$p.value

macrophages_p_value <- wilcox.test(responders %>% filter(annotation == "Macrophages") %>% pull(normalized_count),
                                   non_responders %>% filter(annotation == "Macrophages") %>% pull(normalized_count))$p.value

t_cells_p_value <- wilcox.test(responders %>% filter(annotation == "T_cells") %>% pull(normalized_count),
                               non_responders %>% filter(annotation == "T_cells") %>% pull(normalized_count))$p.value

# Print the p-values
cat("P-value for Keratinocytes:", keratinocytes_p_value, "\n")
cat("P-value for Macrophages:", macrophages_p_value, "\n")
cat("P-value for T_cells:", t_cells_p_value, "\n")

mean_counts <- normalized_counts %>%
  group_by(annotation, responder_status) %>%
  summarise(mean_normalized_count = mean(normalized_count), .groups = 'drop')

# Calculate the difference between the mean normalized counts
mean_counts_diff <- mean_counts %>%
  spread(responder_status, mean_normalized_count) %>%
  mutate(difference = Responder - `Non-responder`)

# Print the differences
print(mean_counts_diff)

# plot
mean_counts_diff <- data.frame(
  annotation = c("Keratinocytes", "Macrophages", "T_cells"),
  Non_responder = c(0.909, 0.0646, 0.0264),
  Responder = c(0.901, 0.0458, 0.0537),
  difference = c(-0.00842, -0.0189, 0.0273),
  p_value = c(0.84, 0.69, 0.31)
)

# Calculate the relative difference
mean_counts_diff <- mean_counts_diff %>%
  mutate(relative_difference = difference / ((Non_responder + Responder) / 2))

# Create the plot
pdf(file = "/path/to/file/",
    width = 7,
    height = 4)
ggplot(mean_counts_diff, aes(x = annotation, y = relative_difference, fill = annotation)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste("p =", round(p_value, 3)), y = ifelse(relative_difference > 0, -0.02, 0.0)), vjust = ifelse(mean_counts_diff$relative_difference > 0, 1, -1)) +
  labs(title = "",
       x = "Cell type",
       y = "Relative difference in normalized count",
       fill = "Cell type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
dev.off()
