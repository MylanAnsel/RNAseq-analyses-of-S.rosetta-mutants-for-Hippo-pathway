# Mylan Ansel, Institut Pasteur, 2025

# This script makes volcanot plots with colors to differentiate DEG that pass the filtering
# steps from the others.

## -----------------------------------------------------------------------------------------------------------------------------------------
# Load required libraries
library(dplyr)
library(ggplot2)
library(stringr)
## -----------------------------------------------------------------------------------------------------------------------------------------

# Define output directory
output_dir <- "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/VolcanoPlot"

# Ensure the output directory exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define input
input_conditions_paths <- list(
  # WT as reference
  Ins_WT_Yorkie = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/YorkievsWT.complete.txt",
  Ins_WT_Warts = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/WartsvsWT.complete.txt",
  Ins_WT_Hippo = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/HippovsWT.complete.txt",
  Strict_WT_Warts = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/filtered_deg_lists/Warts_insertion_WT_filtered.txt",
  Strict_WT_Yorkie = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/filtered_deg_lists/Yorkie_insertion_WT_filtered.txt",
  Strict_WT_Hippo = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/filtered_deg_lists/Hippo_insertion_WT_filtered.txt"
)

# Define thresholds for DEG filtering and padj cap for plot
padj_threshold <- 0.05
log2foldchange_threshold <- 0
padj_cap_value <- 1e-10

## -----------------------------------------------------------------------------------------------------------------------------------------

#' Processes a DEG file, filters by significance, and adds metadata.
#'
#' @param file_path Path to the DEG table.
#' @param cols Columns to select from the input table.
#' @param padj_thresh Adjusted p-value threshold for filtering.
#' @param log2fc_thresh Log2 fold change threshold for filtering.
#' @param df_name_prefix Prefix used to derive comparison and type (e.g., "Del_WT_YorkieUP").
#' @return A processed dataframe with gene_id, comparison, type, species columns.
process_deg_data <- function(file_path, cols, padj_thresh, log2fc_thresh, df_name_prefix) {
  df <- read.table(file_path, header = TRUE, sep = "\t")[, cols]
  colnames(df)[1]<-"Id"
  df <- df %>%
    mutate(significant = ifelse(padj < padj_thresh & abs(log2FoldChange) > log2fc_thresh,"yes","no"))%>%
    mutate(Id = sub("gene:", "", Id))
  
  type_val <- "Deletion" # Default
  if (startsWith(df_name_prefix, "Ins")) {
    type_val <- "all"
  } else {
    type_val <- "filter"
  }
  
  comparison_val <- ""
  if (grepl("_WT_", df_name_prefix)) {
    comparison_val <- paste0(df_name_prefix %>% stringr::str_remove("^(Strict|Ins)_WT_"), "_WT")
  } else {
    warning(paste("Could not determine comparison type for:", df_name_prefix))
    comparison_val <- df_name_prefix # Fallback
  }
  processed_df <- df %>%
    dplyr::rename(gene_id = Id) %>%
    dplyr::mutate(
      comparison = comparison_val,
      type = type_val,
      species = "Salpingoeca_rosetta"
    )
  return(processed_df)
}

#' Loads and combines all relevant DEG tables into a single dataframe for comprehensive filtering.
#'
#' @param conditions_list A named list of file paths for each condition.
#' @param padj_thresh Adjusted p-value threshold.
#' @param log2fc_thresh Log2 fold change threshold.
#' @return A single dataframe containing all processed DEGs.
load_all_deg_data_for_filtering <- function(conditions_list, padj_thresh, log2fc_thresh) {
  processed_dfs_list <- list()
  for (cond_name in names(conditions_list)) {
    cols_to_select <- if (startsWith(cond_name, "Ins")) {
      c(1, 112, 113, 116)
    } else if (startsWith(cond_name, "Strict")) {
      c(1:4)
    } else {
      stop(paste("Unknown condition prefix for:", cond_name))
    }
    
    # Process all files
    processed_dfs_list[[cond_name]] <- process_deg_data(
      file_path = conditions_list[[cond_name]],
      cols = cols_to_select,
      padj_thresh = padj_thresh,
      log2fc_thresh = log2foldchange_threshold,
      df_name_prefix = cond_name
    )
  }
  combined_results <- bind_rows(processed_dfs_list)
  return(combined_results)
}

# --- Data Loading and Preprocessing ---

all_raw_deg_data_combined <- load_all_deg_data_for_filtering(
  conditions_list = input_conditions_paths,
  padj_thresh = padj_threshold,
  log2fc_thresh = log2foldchange_threshold
)

# Remove duplicate rows based on key columns (keeping both 'all' and 'filter' initially for a gene if they exist)
# We will handle the prioritization in the color assignment step
all_raw_deg_data_combined_unique <- all_raw_deg_data_combined %>%
  distinct(gene_id, comparison, type, .keep_all = TRUE)

# Remove rows where padj is NA
all_raw_deg_data_no_na <- all_raw_deg_data_combined_unique %>%
  filter(!is.na(padj))

# Group by gene_id and comparison to determine a single color for each gene within a comparison
all_raw_deg_data_combined_for_plot <- all_raw_deg_data_no_na %>%
  group_by(gene_id, comparison) %>%
  mutate(
    # Check if this gene-comparison pair contains any significant entry
    is_significant_in_group = any(significant == "yes"),
    # Check if this gene-comparison pair contains any "filter" type entry
    is_filtered_in_group = any(type == "filter")
  ) %>%
  ungroup() %>%
  mutate(
    point_color = case_when(
      # If not significant in any entry for this gene-comparison group, it's "not_significant"
      !is_significant_in_group ~ "not_significant",
      # If significant AND present in a "filter" type entry, it's "significant_and_filtered" (highest priority)
      is_significant_in_group & is_filtered_in_group ~ "significant_and_filtered",
      # Otherwise, if significant but ONLY present in "all" types, it's "significant_all"
      is_significant_in_group & !is_filtered_in_group ~ "significant_all",
      TRUE ~ "other" # Fallback, should not be reached with the above logic
    )
  ) %>%
  # Now that the final color is determined for each gene-comparison pair,
  # we can ensure only one row per gene-comparison pair is kept for plotting.
  # This prevents a gene from being plotted twice with different colors.
  distinct(gene_id, comparison, .keep_all = TRUE) %>%
  mutate(
    point_color = factor(point_color, levels = c("not_significant", "significant_all", "significant_and_filtered"))
  )


# Define the custom color palette for the plots
custom_colors <- c(
  "not_significant" = "darkgray",
  "significant_all" = "#ff7f50",
  "significant_and_filtered" = "darkblue"
)

# Get unique comparisons to iterate through for plotting
unique_comparisons <- unique(all_raw_deg_data_combined_for_plot$comparison)

# --- Generate and Save Volcano Plots ---
for (comp in unique_comparisons) {
  # Filter data for the current comparison
  plot_data <- all_raw_deg_data_combined_for_plot %>%
    filter(comparison == comp)
  
  # Create the volcano plot using ggplot2
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = pmin(-log10(padj), -log10(padj_cap_value)),color = point_color))+
    # Add data points with transparency and size
    geom_point(alpha = 0.7, size = 1.5) +
    # Manually set colors and legend labels for the `point_color` factor
    scale_color_manual(values = custom_colors,
                       labels = c("Not significant", "Significant (All)", "Significant (Filtered)"),
                       name = "Gene Status") + # Customize legend title and labels
    geom_hline(yintercept = -log10(padj_threshold), linetype = "longdash", color = "black",linewidth=0.5) +
    geom_vline(xintercept = 0, linetype = "longdash", color = "black",linewidth=0.5) +
    geom_hline(yintercept = -log10(padj_cap_value), linetype = "dotted", color = "black",linewidth=0.5) +
    # Set plot title and axis labels
    labs(
      title = paste("Volcano Plot for", comp),
      x = expression(log[2]~"Fold Change"), # Using expression for proper log2 notation
      y = expression(-log[10]~"Adjusted p-value") # Using expression for proper log10 notation
    ) +
    # Use a minimal theme for a clean look
    theme_minimal() +
    # Customize theme elements for better readability and aesthetics
    theme(
      axis.text.x = element_text(size = 10, angle = 0, vjust = 0, hjust = 0.4, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title =  element_text(size = 13, color = "black"),
      axis.line = element_line(color="black"),
      legend.position = "right", # Position legend on the right
      legend.title = element_text(size = 12, face = "bold"), # Adjust legend title font
      legend.text = element_text(size = 11), # Adjust legend item text font
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") # Add some margin around the plot
    )
  
  print(p)
  # Define file path for saving the plot
  file_name <- paste0("Volcano_Plot_", comp, ".svg")
  output_file_path <- file.path(output_dir, file_name)
  # Save the plot to the specified output directory
  ggsave(  output_file_path <- file.path(output_dir, file_name), plot = p, width = 5, height = 5)
  message(paste("Saved volcano plot for", comp, "to:", output_file_path))
}