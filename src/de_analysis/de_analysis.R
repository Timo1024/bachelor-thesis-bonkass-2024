# ==============================================================================
# File: de_analysis.R
# Description: differential gene expression analysis by plotting an MA plot
# and normalizing
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
library(sf)
library(remotes)
library(RColorBrewer)
library(scales)
library(fields)
library(reshape2)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)

# Main =========================================================================

seu <- LoadSeuratRds("/path/to/file/")

Idents(seu) <- "responder_status"
cells_responder <- WhichCells(seu, idents = "responder")
cells_non_responder <- WhichCells(seu, idents = "non_responder")
de_results <- FindMarkers(seu, ident.1 = "responder", ident.2 = "non_responder", min.pct = 0, thresh.use = 0, logfc.threshold = 0, test.use = "t")

colnames(de_results)[which(names(de_results) == "pct.1")] <- "prevalence_responders"
colnames(de_results)[which(names(de_results) == "pct.2")] <- "prevalence_non_responders"

plot(de_results$p_val_adj)

prevalence_responders <- de_results$prevalence_responders
prevalence_non_responders <- de_results$prevalence_non_responders

de_results$p_val_adj_pseudocount <- de_results$p_val_adj
de_results$p_val_adj_pseudocount[de_results$p_val_adj_pseudocount == 0] <- .Machine$double.xmin

# Combine results into a data frame
results_df <- data.frame(
  gene = rownames(de_results),
  log2FC = de_results$avg_log2FC,
  neg_log10_p_value_adjusted = -log10(de_results$p_val_adj_pseudocount),
  p_value_adjusted = de_results$p_val_adj,
  prevalence_responders = prevalence_responders,
  prevalence_non_responders = prevalence_non_responders
)

# making MA plot
avg_expression <- AggregateExpression(seu, features = rownames(seu), group.by = "responder_status")
avg_expression <- as.data.frame(avg_expression)

# Calculate M (log2 fold change) and A (average expression)
avg_expression$M <- with(avg_expression, log2(RNA.responder / RNA.non.responder))
avg_expression$A <- with(avg_expression, 0.5 * (log2(RNA.responder) + log2(RNA.non.responder)))

# Create the MA plot
ggplot(avg_expression, aes(x = A, y = M)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(x = "Average Expression (A)", y = "Log2 Fold Change (M)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey")  # Line for log2FC 0

# lowess normalization
# Perform LOWESS regression
lowess_fit <- lowess(avg_expression$A, avg_expression$M)

# Create a data frame with the fitted values
lowess_df <- data.frame(A = lowess_fit$x, M_fitted = lowess_fit$y)
lowess_df$gene <- rownames(avg_expression)
avg_expression$gene <- rownames(avg_expression)

avg_expression <- avg_expression %>%
  select(-A)

# Merge the fitted values with the original data
avg_expression <- merge(avg_expression, lowess_df, by = "gene", all.x = TRUE)

# Normalize the data based on the LOWESS fit
avg_expression$M_normalized <- avg_expression$M - avg_expression$M_fitted

# original data with lowess regression line
pdf("/path/to/file/",
    width = 10,
    height = 6)
ggplot(avg_expression, aes(x = A, y = M)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(x = "Average Expression (A)", y = "Log2 Fold Change (M)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey")
dev.off()

# Create the MA plot with the normalized data
ggplot(avg_expression, aes(x = A, y = M_normalized)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(x = "Average Expression (A)", y = "Normalized Log2 Fold Change (M)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey")  # Line for average expression 0

# merge dfs
all_infos_de <- merge(results_df, avg_expression, by = "gene")

log2FC_thresh <- 2

all_infos_de$color <- ifelse(abs(all_infos_de$M_normalized) > log2FC_thresh , ifelse(all_infos_de$M_normalized > log2FC_thresh, "#e87d72", "#56bcc2"), "grey")
all_infos_de$mark_text <- ifelse(abs(all_infos_de$M_normalized) > log2FC_thresh & all_infos_de$neg_log10_p_value_adjusted > 50, "yes", "no")

pdf("/path/to/file/",
    width = 10,
    height = 6)
ggplot(all_infos_de, aes(x = A, y = M_normalized)) +
  geom_point(alpha = 0.6, aes(color = color, size = neg_log10_p_value_adjusted)) +
  scale_size("-log10(p-value)", range = c(0,3)) +
  scale_color_manual(
    name = "",
    values = c("#4433dd", "#dd4422", "grey"),
    labels = c("Downregulated", "Upregulated", "Not significant")
  ) +
  theme_minimal() +
  labs(x = "Average Expression (A)", y = "Normalized Log2 Fold Change (M)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = log2FC_thresh, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log2FC_thresh, linetype = "dashed", color = "grey") +
  geom_text_repel(
    aes(label = ifelse(mark_text != "no", as.character(gene), '')),
    size = 3,
    max.overlaps = 100
  ) +
  annotate("text", x = 16, y = 2.3, label = paste0("Log2FC > ", log2FC_thresh), size = 4, color = "grey") +
  annotate("text", x = 16, y = -2.3, label = paste0("Log2FC < -", log2FC_thresh), size = 4, color = "grey")
dev.off()
