# ==============================================================================
# File: cellchat_getting_communications.R
# Description: inferring the CCC networks and visualizing them as a circos plot
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
library(ComplexHeatmap)
library(RColorBrewer)
library(reshape2)

# Main =========================================================================

# load cellchat objects
seu <- LoadSeuratRds("/path/to/file/")

patients <- unique(seu$patient_id)

cellchat_list <- list()
for (patient in patients) {
  # Construct the filename
  filename <- paste0("/path/to/file/", patient, "/path/to/file/")

  # Load the cellchat object from the file and store it in the list
  cellchat_list[[patient]] <- readRDS(file = filename)
}

circos_plot <- function(ligand_receptor_frame,
                        cell_group_colors,
                        patient,
                        ligand_color="blue",
                        receptor_color="red",
                        cex_outer=0.5,
                        cex_inner=0.4) {
  
  # Reformat data
  part1 <- ligand_receptor_frame %>%
    mutate(lig = sapply(strsplit(interaction, split = "_"), function(x) x[[1]])) %>%
    mutate(rec = sapply(strsplit(interaction, split = "_"), function(x) x[[2]])) %>%
    select(cell_type1, lig, value) %>%
    distinct() %>%
    mutate(type = "lig")
  
  part2 <- ligand_receptor_frame %>%
    mutate(lig = sapply(strsplit(interaction, split = "_"), function(x) x[[1]])) %>%
    mutate(rec = sapply(strsplit(interaction, split = "_"), function(x) x[[2]])) %>%
    select(cell_type2, rec, value) %>%
    distinct() %>%
    mutate(type = "rec")
  
  colnames(part1) <- colnames(part2) <- c("classes", "lig.rec", "value", "type")
  
  part12 <- rbind(part1, part2) %>%
    group_by(classes) %>%
    group_split() %>%
    lapply(function(x) {
      x %>%
        mutate(ordered.lig.rec = paste(type, lig.rec, sep = "_")) %>%
        mutate(ranges = as.numeric(as.factor(ordered.lig.rec))) %>%
        select(-ordered.lig.rec)
    }) %>%
    do.call(rbind, .)
  
  to.join <- ligand_receptor_frame %>%
    mutate(lig = sapply(strsplit(interaction, split = "_"), function(x) x[[1]])) %>%
    mutate(rec = sapply(strsplit(interaction, split = "_"), function(x) x[[2]])) %>%
    select(cell_type1, cell_type2, lig, rec, value)
  
  colnames(to.join)[1:3] <- c("classes", "to.class", "lig.rec")
  
  part3 <- part12
  
  joined <- left_join(part3, to.join, by = c("classes", "lig.rec"))
  joined$to.rec <- NA
  joined$value <- joined$value.x
  joined <- joined %>%
    select(-value.x, -value.y)
  
  for (i in 1:nrow(joined)) {
    sub.group <- joined[i,]
    sub.joined <- joined %>% filter(classes == sub.group$to.class)
    joined$to.rec[i] <- sub.joined[match(sub.group$rec, sub.joined$lig.rec), "ranges"] %>% pull()
  }
  
  final.construct <- joined
  
  # Repair single class
  single.class <- final.construct %>%
    group_by(classes) %>%
    summarize(max_range = max(ranges)) %>%
    filter(max_range == 1) %>%
    pull(classes)
  
  if (length(single.class) > 0) {
    for (i in 1:length(single.class)) {
      row.add <- final.construct[final.construct$classes == single.class[i],][1,]
      row.add$ranges <- 2
      final.construct <- rbind(final.construct, row.add)
      final.construct <- final.construct %>%
        arrange(classes)
    }
  }
  
  final.construct <- final.construct %>%
    arrange(classes, ranges)

  circos.clear()
  circos.par(gap.degree = 10, track.margin = c(0, 0.2))
  circos.initialize(factors = final.construct$classes, x = final.construct$ranges)
  
  # Outer track: Cell type labels
  circos.track(
    ylim = c(0, 1),
    track.height = 0.1,
    panel.fun = function(x, y) {
      circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1], CELL_META$cell.xlim[2], CELL_META$cell.ylim[2],
                  # col = "white")
                  col = cell_group_colors[CELL_META$sector.numeric.index])
      circos.text(CELL_META$xcenter, y = 2.5, CELL_META$sector.index, facing = "downward", cex = cex_outer)
    }
  )
  
  # Inner track: Ligand and receptor labels
  circos.track(
    ylim = c(0, 1),
    track.height = 0.05,
    bg.border = "white"
  )
  
  final.construct2 <- final.construct %>%
    select(classes, lig.rec, ranges, type) %>%
    distinct() %>%
    arrange(classes, ranges)
  
  ref.tab <- unname(table(final.construct2$classes))
  sec.multi <- (ref.tab - 1) / ref.tab
  names(sec.multi) <- names(table(final.construct2$classes))
  
  int.types.list <- final.construct2 %>%
    group_split(classes)
  
  names(int.types.list) <- sapply(int.types.list, function(x) x$classes[1])
  
  int.types.list.multi <- int.types.list[!names(int.types.list) %in% single.class]
  int.types.list.individ <- int.types.list[names(int.types.list) %in% single.class]
  
  for (i in 1:length(int.types.list.multi)) {
    for (a in 1:nrow(int.types.list.multi[[i]])) {
      if (a == 1) {
        sec.multi.use <- sec.multi[names(sec.multi) == int.types.list.multi[[i]]$classes[1]]
        circos.rect(1, 0, 1 + sec.multi.use * a, 1, sector.index = int.types.list.multi[[i]]$classes[a],
                    col = ifelse(int.types.list.multi[[i]]$type[a] == "lig", ligand_color, receptor_color), track.index = 2)
        circos.text(1 + sec.multi.use * a / 2, 4, sector.index = int.types.list.multi[[i]]$classes[a],
                    labels = int.types.list.multi[[i]]$lig.rec[a], track.index = 2, facing = "downward", cex = cex_inner)
      } else {
        sec.multi.use <- sec.multi[names(sec.multi) == int.types.list.multi[[i]]$classes[1]]
        circos.rect(1 + sec.multi.use * (a - 1), 0, 1 + sec.multi.use * a, 1, sector.index = int.types.list.multi[[i]]$classes[a],
                    col = ifelse(int.types.list.multi[[i]]$type[a] == "lig", ligand_color, receptor_color), track.index = 2)
        circos.text(1 + sec.multi.use * a - sec.multi.use / 2, 4, sector.index = int.types.list.multi[[i]]$classes[a],
                    labels = int.types.list.multi[[i]]$lig.rec[a], track.index = 2, facing = "downward", cex = cex_inner)
      }
    }
  }
  
  if (length(int.types.list.individ) > 0) {
    for (i in 1:length(int.types.list.individ)) {
      circos.rect(1, 0, 2, 1, sector.index = int.types.list.individ[[i]]$classes[1],
                  col = ifelse(int.types.list.individ[[i]]$type[1] == "lig", ligand_color, receptor_color), track.index = 2)
      circos.text(1.5, 4, sector.index = int.types.list.individ[[i]]$classes[1],
                  labels = int.types.list.individ[[i]]$lig.rec[1], track.index = 2, facing = "downward", cex = cex_inner)
    }
  }
  
  # Draw links
  final.construct3 <- joined %>%
    select(classes, lig.rec, ranges, to.class, to.rec, value) %>%
    distinct()
  
  split.construct <- final.construct3 %>%
    split(.$classes)
  
  final.construct3 <- lapply(split.construct, function(x) {
    class.length <- length(unique(x$ranges))
    if (class.length == 1) {
      x[,"ranges"] <- 1.5
      x
    } else {
      x
    }
  }) %>%
    do.call(rbind, .)
  
  int.types.list <- final.construct3 %>%
    group_split(classes)
  
  names(int.types.list) <- sapply(int.types.list, function(x) x$classes[1])
  
  colorramp <- viridis(100)
  
  for (i in 1:length(int.types.list)) {
    for (a in 1:nrow(int.types.list[[i]])) {
      target <- which(!is.na(match(names(int.types.list), int.types.list[[i]]$to.class[[a]])))
      
      all_colors <- colors()
      random_color <- sample(all_colors, 1)
      
      current_color <- colorramp[ceiling(int.types.list[[i]]$value[a] * 100)]
      
      scale <- int.types.list[[i]]$value[a]
      
      print(scale)
      
      if (length(target) > 0) {
        if (!int.types.list[[i]]$to.class[[a]] %in% single.class) {
          circos.link(int.types.list[[i]]$classes[a], 1 + sec.multi[i] * int.types.list[[i]]$ranges[a] - sec.multi[i] / 2,
                      int.types.list[[i]]$to.class[[a]], 1 + sec.multi[target] * int.types.list[[i]]$to.rec[a] - sec.multi[target] / 2,
                      0.43, 0.43, directional = 0,
                      col = "black",
                      # arr.width = 0.2,
                      # arr.length = 0.2,
                      lwd = min(c(scale * 5, 5)))
        } else {
          circos.link(int.types.list[[i]]$classes[a], 1 + sec.multi[i] * int.types.list[[i]]$ranges[a] - sec.multi[i] / 2,
                      int.types.list[[i]]$to.class[[a]], 1.5,
                      0.43, 0.43, directional = 0,
                      col = "black",
                      # arr.width = 0.2,
                      # arr.length = 0.2,
                      lwd = min(c(scale * 5, 5)))
        }
      }
    }
  }
}

{
  current_patient <- "10"
  
  pdf(paste0("/path/to/file/", current_patient, "/path/to/file/"),
      width = 6, # 14
      height = 6)
  # Set up the multi-panel plot layout for 4 plots in a single row
  par(mfrow = c(1, 1))
  
  # no interactions in 4
  for (patient in c(current_patient)) {
    
    cc <- cellchat_list[[patient]]
  
    interaction_data <- subsetCommunication(cc)
    
    interaction_data$source <- gsub("cluster", "", interaction_data$source)
    interaction_data$target <- gsub("cluster", "", interaction_data$target)
    
    interaction_data$ligand <- sapply(interaction_data$ligand, function(x) {
      if(x == "PECAM1") {
        return(paste0("L-", x))
      } else {
        return(x)
      }
    })
    interaction_data$receptor <- sapply(interaction_data$receptor, function(x) {
      if(x == "PECAM1") {
        return(paste0("R-", x))
      } else {
        return(x)
      }
    })
    
    ligand_receptor_frame_2 <- as.tibble(data.frame(
      cell_type1 = as.character(interaction_data$source),
      cell_type2 = as.character(interaction_data$target),
      interaction = as.character(paste0(interaction_data$ligand, "_", interaction_data$receptor)),
      value = interaction_data$prob * 300
    ))
    
    cell_group_colors <- c(
      "Keratinovytes" = "lightgrey",
      "B-cell" = "lightblue",
      "DC" = "lightgreen",
      "T-cells" = "lightcoral",
      "NK-cell" = "lightgoldenrod",
      "Macrophage" = "lightpink"
    )
    
    circos_plot(ligand_receptor_frame_2, cell_group_colors, patient)
    
  }
  dev.off()
}

table(subsetCommunication(cellchat_list[["1"]])$annotation)

table(subsetCommunication(cellchat_list[["1"]]))

head(cellchat_list[["2"]]@net$weight)
head(cellchat_list[["2"]]@net$pval)

# Function to add a pseudocount and scale a matrix, handling constant columns
scale_matrix_with_pseudocount <- function(x, reference_order) {
  x <- x[reference_order, reference_order]
  return(x)
}

reference_order <- rownames(cellchat_list[["1"]]@net$weight)

# Add a pseudocount and scale the weight data frames for each sample
scaled_weight_list <- lapply(cellchat_list[1:10], function(x) scale_matrix_with_pseudocount(x@net$weight, reference_order))

# Combine the scaled data frames into a matrix
scaled_combined_matrix <- do.call(cbind, scaled_weight_list)

# Create a vector indicating the sample each column belongs to
sample_labels <- rep(paste0("Sample_", 1:10), each = ncol(cellchat_list[[1]]@net$weight))

# Convert sample_labels to a factor with levels in the correct order
sample_labels <- factor(sample_labels, levels = paste0("Sample_", 1:10))

# Create a vector indicating the cell type each column belongs to
cell_types <- rep(colnames(cellchat_list[[1]]@net$weight), times = 10)

# Create the annotation bar
column_annotation <- HeatmapAnnotation(
  CellType = cell_types,
  col = list(`Target cell type` = c("Keratinocytes" = "red", "Macrophages" = "blue", "DC" = "green", 
                          "T_cells" = "purple", "NK_cell" = "orange", "B_cell" = "yellow"))
)

set2_colors <- brewer.pal(n = 6, name = "Set2")
names(set2_colors) <- colnames(cellchat_list[[1]]@net$weight)
column_annotation <- HeatmapAnnotation(
  `Target cell type` = cell_types,
  col = list(`Target cell type` = set2_colors)
)

# Create the heatmap with the viridis color palette
ht <- Heatmap(scaled_combined_matrix, 
        name = "Interaction Weights", 
        column_split = sample_labels,
        top_annotation = column_annotation,
        col = viridis(256),
        cluster_rows = TRUE, 
        cluster_columns = FALSE, 
        show_column_names = FALSE,
        show_row_names = TRUE,
        heatmap_legend_param = list(title = "Weight"),
        row_title = "Source (Ligands)")

pdf(paste0("/path/to/file/"),
    width = 14, # 14
    height = 3)
draw(ht, column_title = "Target (Receptor)")
dev.off()

# differences
# Function to extract values at the same position from all matrices
extract_values_at_position <- function(weight_list, row, col) {
  sapply(weight_list, function(x) x[row, col])
}

# Get the dimensions of the matrices
n_rows <- nrow(scaled_weight_list[[1]])
n_cols <- ncol(scaled_weight_list[[1]])

# Initialize matrices to store the changes and p-values
change_matrix <- matrix(0, n_rows, n_cols)
p_value_matrix <- matrix(0, n_rows, n_cols)

# Loop through each element in the matrix
for (i in 1:n_rows) {
  for (j in 1:n_cols) {
    # Extract values at the same position from all matrices
    values <- extract_values_at_position(scaled_weight_list, i, j)
    
    # Split the values into responders and non-responders
    responders_values <- values[1:5]
    non_responders_values <- values[6:10]
    
    # Calculate the change in mean interaction strength
    change_matrix[i, j] <- mean(responders_values) - mean(non_responders_values)
    
    # Perform a t-test to get the p-value
    t_test_result <- t.test(responders_values, non_responders_values)
    p_value_matrix[i, j] <- t_test_result$p.value
  }
}

# Set row and column names for the matrices
rownames(change_matrix) <- rownames(scaled_weight_list[[1]])
colnames(change_matrix) <- colnames(scaled_weight_list[[1]])
rownames(p_value_matrix) <- rownames(scaled_weight_list[[1]])
colnames(p_value_matrix) <- colnames(scaled_weight_list[[1]])

# Print the results
cat("Change in mean interaction strength for each pair:\n")
print(change_matrix)
cat("\nP-values for each pair:\n")
print(p_value_matrix)

# Convert matrices to data frames for plotting
change_matrix_scaled <- scale(change_matrix, center = FALSE)
change_df <- melt(change_matrix_scaled)
colnames(change_df) <- c("Row", "Column", "Change")

p_value_df <- melt(p_value_matrix)
colnames(p_value_df) <- c("Row", "Column", "PValue")

# Merge the data frames
merged_df <- merge(change_df, p_value_df, by = c("Row", "Column"))

merged_df$PValue[is.na(merged_df$PValue)] <- 0
merged_df$Change[is.na(merged_df$Change)] <- 0

dot_plot <- ggplot(merged_df, aes(x = Column, y = Row, color = Change, size = -log10(PValue))) +
  geom_point() +
  scale_color_viridis(option = "E") +
  # scale_color_viridis(option = "viridis") +
  scale_size_continuous(range = c(1, 10)) +
  theme_minimal() +
  labs(title = "", x = "Target (Receptor)", y = "Source (Ligand)", size = "-log10(p-value)", color = "Change of \naverage strength")

# Print the plot
pdf(paste0("/path/to/file/"),
    width = 9, # 14
    height = 4)
print(dot_plot)
dev.off()


