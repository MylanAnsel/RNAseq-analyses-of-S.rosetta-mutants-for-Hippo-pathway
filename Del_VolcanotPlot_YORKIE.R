# Mylan Ansel, Institut Pasteur, 2025

# This script makes a volcanot plots for YorkieKO deletion in which DEG shared with
# consistent in WartsKO are also highlighted.
# It also makes the comparison between log2FoldChange of Yorkie genes which are consistent DEG in Warts-KO
# and all the other ones


## -----------------------------------------------------------------------------------------------------------------------------------------
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(ggrepel)
## -----------------------------------------------------------------------------------------------------------------------------------------
# Define output directory
output_dir <- "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses"

# Ensure the output directory exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define input
input_conditions_paths <- list(
  # WT as reference
  Del_WT_Yorkie = "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/tables/YCvsSr.complete.txt",
  Del_WT_Warts = "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/tables/WBvsSr.complete.txt",
  Strict_WT_Warts = "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Heatmap/filtered_deg_lists/insertion_intersected_direction/Warts_strict_filtered_DEGs_intersect_direction.txt",
  Strict_WT_Yorkie = "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Heatmap/filtered_deg_lists/Yorkie_strict_filtered_DEGs.txt"
)

# Define thresholds for DEG filtering
padj_threshold <- 0.05
log2foldchange_threshold <- 0
padj_cap_value <- 1e-1

## -----------------------------------------------------------------------------------------------------------------------------------------

# Processes a DEG file, filters by significance, and adds metadata.
#
# @param file_path Path to the DEG table.
# @param cols Columns to select from the input table.
# @param padj_thresh Adjusted p-value threshold for filtering.
# @param log2fc_thresh Log2 fold change threshold for filtering.
# @param df_name_prefix Prefix used to derive comparison and type (e.g., "Del_WT_YorkieUP").
# @return A processed dataframe with gene_id, comparison, type, species columns.
process_deg_data <- function(file_path, cols, padj_thresh, log2fc_thresh, df_name_prefix) {
  df <- read.table(file_path, header = TRUE, sep = "\t")[, cols]
  colnames(df)[1]<-"Id"
  df <- df %>%
    mutate(significant = ifelse(padj < padj_thresh & abs(log2FoldChange) > log2fc_thresh,"yes","no"))%>%
    mutate(Id = sub("gene:", "", Id))
  
  type_val <- "Deletion" # Default
  if (startsWith(df_name_prefix, "Del")) {
    type_val <- "all"
  } else {
    type_val <- "filter"
  }
  
  comparison_val <- ""
  if (grepl("_WT_", df_name_prefix)) {
    comparison_val <- paste0(df_name_prefix %>% stringr::str_remove("^(Strict|Del)_WT_"), "_WT")
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

# Loads and combines all relevant DEG tables into a single dataframe for comprehensive filtering.
#
# @param conditions_list A named list of file paths for each condition.
# @param padj_thresh Adjusted p-value threshold.
# @param log2fc_thresh Log2 fold change threshold.
# @return A single dataframe containing all processed DEGs.
load_all_deg_data_for_filtering <- function(conditions_list, padj_thresh, log2fc_thresh) {
  processed_dfs_list <- list()
  for (cond_name in names(conditions_list)) {
    cols_to_select <- if (startsWith(cond_name, "Del")) {
      c(1, 38, 39, 42)
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

# Group by gene_id and comparison to determine relevant flags for each gene within a comparison
all_raw_deg_data_for_color_determination <- all_raw_deg_data_no_na %>%
  group_by(gene_id, comparison) %>%
  mutate(
    # Check if this gene-comparison pair contains any significant entry
    is_significant_in_group = any(significant == "yes"),
    # Check if this gene-comparison pair contains any "filter" type entry
    is_filtered_in_group = any(type == "filter")
  ) %>%
  ungroup()

# --- Identify genes significant AND filtered in BOTH Yorkie_WT and Warts_WT ---

# Identify genes that are significant_and_filtered in Yorkie_WT
yorkie_sig_filtered_genes <- all_raw_deg_data_for_color_determination %>%
  filter(comparison == "Yorkie_WT", is_significant_in_group == TRUE, is_filtered_in_group == TRUE) %>%
  pull(gene_id) %>%
  unique()

# Identify genes that are significant_and_filtered in Warts_WT
warts_sig_filtered_genes <- all_raw_deg_data_for_color_determination %>%
  filter(comparison == "Warts_WT", is_significant_in_group == TRUE, is_filtered_in_group == TRUE) %>%
  pull(gene_id) %>%
  unique()

# Find genes that are significant_and_filtered in BOTH Yorkie_WT and Warts_WT
shared_sig_filtered_genes <- intersect(yorkie_sig_filtered_genes, warts_sig_filtered_genes)

# Identify genes that are significant in Warts_WT (any type), but not in the "both" category
warts_only_significant_genes_for_yorkie_plot <- setdiff(warts_sig_filtered_genes, shared_sig_filtered_genes)

# --- Assign final point_color, prioritizing the "both" condition for Yorkie_WT plot ---
all_raw_deg_data_combined_for_plot <- all_raw_deg_data_for_color_determination %>%
  mutate(
    point_color = case_when(
      # Highest priority: significant_and_filtered in both AND current comparison is Yorkie_WT
      (comparison == "Yorkie_WT" & gene_id %in% shared_sig_filtered_genes) ~ "significant_and_filtered_both",
      # Next priority: significant in Warts_WT (but not "both") AND current comparison is Yorkie_WT
      (comparison == "Yorkie_WT" & gene_id %in% warts_only_significant_genes_for_yorkie_plot) ~ "significant_in_warts_only",
      # Standard significant_and_filtered for the current comparison (if not caught by above)
      is_significant_in_group & is_filtered_in_group ~ "significant_and_filtered",
      # All others are not significant
      TRUE ~ "not_significant"
    ),
    # Adjust point size for specific categories
    point_size_raw = case_when(
      point_color == "significant_and_filtered_both" ~ 3,
      point_color == "significant_in_warts_only" ~ 2.5,
      point_color == "significant_and_filtered" ~ 1.5,# Slightly larger for Warts only
      TRUE ~ 1
    )
  ) %>%
  # Now that the final color and raw size are determined for each gene-comparison pair,
  # we can ensure only one row per gene-comparison pair is kept for plotting.
  # This prevents a gene from being plotted twice with different colors.
  distinct(gene_id, comparison, .keep_all = TRUE) %>%
  mutate(
    # Re-factor with the new level. Ensure order is logical for legend.
    point_color = factor(point_color, levels = c("not_significant", "significant_and_filtered", "significant_in_warts_only", "significant_and_filtered_both")),
    # Convert point_size to a factor to be used with scale_size_manual
    point_size = factor(point_size_raw, levels = c("1","1.5", "2.5", "3"))
  )


# Define the custom color palette for the plots
custom_colors <- c(
  "not_significant" = "darkgray",
  "significant_and_filtered" = "darkgray",
  "significant_in_warts_only" = "blue", # New color for Warts only
  "significant_and_filtered_both" = "darkblue" # Color for the combined category
)

# Get unique comparisons to iterate through for plotting
unique_comparisons <- "Yorkie_WT"

# --- Generate and Save Volcano Plots ---
for (comp in unique_comparisons) {
  # Filter data for the current comparison
  plot_data <- all_raw_deg_data_combined_for_plot %>%
    filter(comparison == comp)
  
  # Filter for genes to label (only for Yorkie_WT and only the "both" category)
  label_data <- plot_data %>%
    filter(comparison == "Yorkie_WT", point_color == "significant_and_filtered_both")
  
  # Create the volcano plot using ggplot2
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = pmin(-log10(padj), -log10(padj_cap_value)), color = point_color, size = point_size, alpha = point_size)) +
    # Add data points with transparency and size mapped to `point_size`
    geom_point() +
    # Manually set colors and legend labels for the `point_color` factor
    scale_color_manual(values = custom_colors,
                       labels = c("Not significant", "Significant (Filtered)", "Significant (Warts Only)", "Significant (Yorkie & Warts Filtered)"), # Updated labels
                       name = "Gene Status") + # Customize legend title and labels
    # Manually set sizes for the `point_size` aesthetic
    scale_size_manual(values = c("1" = 0.5, "1.5" = 0.5,"2.5" = 1, "3" = 1), guide = "none") + # Hide size legend
    # Manually set alpha for the `point_size` aesthetic
    scale_alpha_manual(values = c("1" = 0.2,"1.5" = 0.2, "2.5" = 0.9, "3" = 0.9), guide = "none") + # Adjusted alpha values
    
    geom_hline(yintercept = -log10(padj_threshold), linetype = "longdash", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "longdash", color = "black", linewidth = 0.5) +
    geom_hline(yintercept = -log10(padj_cap_value), linetype = "dotted", color = "black", linewidth = 0.5) +
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
  file_name <- paste0("Volcano_Plot_Yorkie_shared_Warts.svg")
  output_file_path <- file.path(output_dir, file_name)
  # Save the plot to the specified output directory
  #ggsave(output_file_path, plot = p, width = 10, height = 5)
}

# Add labels only for the "Yorkie_WT" plot and for the specific "both" genes
if (comp == "Yorkie_WT" && nrow(label_data) > 0) {
  p <- p + geom_text_repel(data = label_data,
                           aes(label = gene_id),
                           box.padding = 0.5,
                           point.padding = 0.5,
                           segment.color = 'grey50',
                           size = 3.5, # Adjust label text size
                           color = "black") # Set label text color
}


# --- Generate Histograms for Yorkie_WT log2FoldChange ---

# Data for Yorkie_WT (all genes)
yorkie_all_genes_data <- all_raw_deg_data_combined_for_plot %>%
  filter(comparison == "Yorkie_WT" & (point_color == "significant_and_filtered" | point_color == "not_significant"))%>%
  mutate(condition="yorkie")%>%
  dplyr::select(log2FoldChange,condition) 

# Get gene_ids from Strict_WT_Warts
strict_warts_genes <- all_raw_deg_data_combined_for_plot %>%
  filter(comparison == "Yorkie_WT" & (point_color == "significant_in_warts_only" | point_color == "significant_and_filtered_both")) %>%
  mutate(condition="warts")%>%
  dplyr::select(log2FoldChange,condition) 


# Merge the two datasets by gene_id
superplot_data <- rbind(yorkie_all_genes_data, strict_warts_genes)%>%
  mutate(log2FoldChange=as.numeric(log2FoldChange),
         condition=as.factor(condition))

Averages <- superplot_data %>%
  group_by(condition) %>% 
  summarize(across(log2FoldChange, mean))

p_violin <- ggplot(superplot_data, aes(x = condition, y = log2FoldChange)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) + # Add horizontal line at 0
  geom_point(aes(color = condition),size=1,position = position_jitterdodge(jitter.width = 0.9, dodge.width = 1),alpha=0.8)+
  geom_beeswarm(data=Averages,aes(fill=condition),shape=23, size=4,cex=3) +
  
  stat_compare_means(data=superplot_data, aes(x=condition,y=log2FoldChange),
                     comparisons = list(c("warts",  "yorkie")), method="wilcox.test", 
                     paired=FALSE,size=7, shape = 18,cex=2,label.x = 1.2, label.y = c(1.5),
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) +
  
  labs(title = "",x = "",y = expression(log[2]~"Fold Change")) +
  scale_color_manual(values = c("yorkie" = "darkgray", "warts" = "darkblue")) + 
  scale_fill_manual(values = c("yorkie" = "darkgray", "warts" = "darkblue")) + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    axis.line = element_line(color="black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )+
 ylim(-2.2,2)

print(p_violin)

# Define file path for saving the violin plot
file_name_violinplot <- "Dotplot_YorkieKO-all_vs_YorkieKO-Warts_log2FC.svg"
output_file_path_violinplot <- file.path(output_dir, file_name_violinplot)
ggsave(output_file_path_violinplot, plot = p_violin, width = 4, height = 6)

