# ==============================================================================
# File: celltalk_compare_samples.R
# Description: Infer CCCs by using celltalker and plot them
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

setwd("/path/to/file/")

install.packages("devtools")
devtools::install()

library(celltalkerModified)

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

# replace all _ with - in SingleR_HCA and SingleR_HCA_fine metadata columns
seu$SingleR_HCA <- gsub("_", "-", seu$SingleR_HCA)

ct <- celltalkerModified::celltalk(
  input_object = seu,
  metadata_grouping = "SingleR_HCA",
  ligand_receptor_pairs = ramilowski_pairs,
  number_cells_required = 5,
  min_expression = 10,
  max_expression = 100000,
  scramble_times = 10
)

compared <- celltalkerModified::compare_group_interactions(
  seu,
  ct,
  "patient_id",
  "responder_status",
  "SingleR_HCA"
)

significant_interactions <- lapply(compared, function(model) {
  summary_model <- summary(model)
  p_value <- summary_model$coefficients[2, 4] # Extracting the p-value for the sample_groupsresponder coefficient
  if (p_value < 0.05) {
    return(model)
  } else {
    return(NULL)
  }
})

# Remove NULL values from the list
significant_interactions <- Filter(Negate(is.null), significant_interactions)

# Print the names of significant interactions
names(significant_interactions)

# make a df out of significant interactions
# Create an empty data frame to store the results
results <- data.frame(Ligand = character(),
                      Receptor = character(),
                      CellType1 = character(),
                      CellType2 = character(),
                      pValue = numeric(),
                      Difference = numeric(),
                      stringsAsFactors = FALSE)

# List of your models
models <- significant_interactions

# Function to extract information from a model
extract_info <- function(model, model_name) {
  # Split the model name at underscores
  parts <- unlist(strsplit(model_name, "_"))
  celltype1 <- parts[1]
  ligand <- parts[2]
  celltype2 <- parts[3]
  receptor <- parts[4]
  
  summary_model <- summary(model)
  p_value <- summary_model$coefficients[2, 4]
  difference <- summary_model$coefficients[2, 1]
  
  return(data.frame(Ligand = ligand,
                    Receptor = receptor,
                    CellType1 = celltype1,
                    CellType2 = celltype2,
                    pValue = p_value,
                    Difference = difference,
                    stringsAsFactors = FALSE))
}

# Extract information from each model and add it to the results data frame
for (model_name in names(models)) {
  results <- rbind(results, extract_info(models[[model_name]], model_name))
}

# View the results
print(results)

# combine same cell type pairs
# Function to combine p-values and calculate mean difference
combine_pvalues_and_mean_diff <- function(data) {
  combined_results <- data.frame(CellType1 = character(),
                                 CellType2 = character(),
                                 CombinedPValue = numeric(),
                                 MeanDifference = numeric(),
                                 stringsAsFactors = FALSE)
  
  unique_pairs <- unique(data[, c("CellType1", "CellType2")])
  
  for (i in 1:nrow(unique_pairs)) {
    pair <- unique_pairs[i, ]
    subset_data <- subset(data, CellType1 == pair$CellType1 & CellType2 == pair$CellType2)
    
    if (nrow(subset_data) == 1) {
      combined_pvalue <- subset_data$pValue
    } else {
      combined_pvalue <- sumlog(subset_data$pValue)$p
    }
    
    # combined_pvalue <- sumlog(subset_data$pValue)$p
    mean_difference <- mean(subset_data$Difference)
    
    combined_results <- rbind(combined_results, data.frame(CellType1 = pair$CellType1,
                                                           CellType2 = pair$CellType2,
                                                           CombinedPValue = combined_pvalue,
                                                           MeanDifference = mean_difference,
                                                           stringsAsFactors = FALSE))
  }
  
  return(combined_results)
}

# Combine p-values and calculate mean differences
combined_results <- combine_pvalues_and_mean_diff(results)

# View the results
print(combined_results)

order <- c("DC", "Keratinocytes", "T-cells", "NK-cell", "Macrophages", "B-cell")
combined_results$CellType1 <- factor(combined_results$CellType1, levels = order)
combined_results$CellType2 <- factor(combined_results$CellType2, levels = order)

# plot
dot_plot <- ggplot(combined_results, aes(x = CellType1, y = CellType2, color = MeanDifference, size = -log10(CombinedPValue))) +
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

# chord diagram
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
                    labels = int.types.list.multi[[i]]$lig.rec[a], track.index = 2, facing = "clockwise", niceFacing = TRUE, cex = cex_inner)
      } else {
        sec.multi.use <- sec.multi[names(sec.multi) == int.types.list.multi[[i]]$classes[1]]
        circos.rect(1 + sec.multi.use * (a - 1), 0, 1 + sec.multi.use * a, 1, sector.index = int.types.list.multi[[i]]$classes[a],
                    col = ifelse(int.types.list.multi[[i]]$type[a] == "lig", ligand_color, receptor_color), track.index = 2)
        circos.text(1 + sec.multi.use * a - sec.multi.use / 2, 4, sector.index = int.types.list.multi[[i]]$classes[a],
                    labels = int.types.list.multi[[i]]$lig.rec[a], track.index = 2, facing = "clockwise", niceFacing = TRUE, cex = cex_inner)
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
      
      scale <- int.types.list[[i]]$value[a]
      
      # Set color based on the sign of the value
      current_color <- ifelse(scale < 0, "red", "blue")
      
      # Set line width based on the absolute value of the scale (value)
      line_width <- abs(scale) * 5  # Adjust the multiplier as needed to control the thickness
      
      if (length(target) > 0) {
        if (!int.types.list[[i]]$to.class[[a]] %in% single.class) {
          circos.link(int.types.list[[i]]$classes[a], 1 + sec.multi[i] * int.types.list[[i]]$ranges[a] - sec.multi[i] / 2,
                      int.types.list[[i]]$to.class[[a]], 1 + sec.multi[target] * int.types.list[[i]]$to.rec[a] - sec.multi[target] / 2,
                      0.43, 0.43, directional = 0,
                      col = current_color,    # Use the current color based on the value's sign
                      lwd = min(c(line_width, 5)))  # Use line width based on the absolute value
        } else {
          circos.link(int.types.list[[i]]$classes[a], 1 + sec.multi[i] * int.types.list[[i]]$ranges[a] - sec.multi[i] / 2,
                      int.types.list[[i]]$to.class[[a]], 1.5,
                      0.43, 0.43, directional = 0,
                      col = current_color,    # Use the current color based on the value's sign
                      lwd = min(c(line_width, 5)))  # Use line width based on the absolute value
        }
      }
    }
  }
}

lr_frame <- results %>%
  # Create interaction column by concatenating Ligand and Receptor
  mutate(interaction = paste("L-", Ligand, "_R-", Receptor, sep = "")) %>%
  # Select only the relevant columns and rename them
  select(cell_type1 = CellType1, cell_type2 = CellType2, interaction, value = Difference) %>%
  # Convert to a tibble
  as_tibble()

# View the result
lr_frame

color_scheme <- brewer.pal(6, "Set3")

# Assign colors to cell types
cell_group_colors <- c("B-cell" = color_scheme[1], 
                       "DC" = color_scheme[2], 
                       "T-cells" = color_scheme[3], 
                       "Keratinocytes" = color_scheme[4], 
                       "NK-cell" = color_scheme[5], 
                       "Macrophages" = color_scheme[6])
pdf(paste0("/path/to/file/"),
    width = 6, # 14
    height = 6)
circos_plot(lr_frame, cell_group_colors, "Patient1")
dev.off()











