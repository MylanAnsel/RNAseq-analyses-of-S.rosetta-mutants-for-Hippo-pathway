# Mylan Ansel, Institut Pasteur, 2025

# This script save DEG subsets based on their origin, direction (UP/DOWN), and filter used 
# It also generates Heatmaps with protein long name extracted from Uniprot using API.

## -----------------------------------------------------------------------------------------------------------------------------------------
# Load required libraries
library(dplyr)
library(ggplot2)
library(viridis)
library(cluster)
library(readr)
library(stringr)
library(text2vec)
library(tidyverse)
library(ggupset)
library(stringdist)
library(ComplexUpset)
library(UpSetR) # For the separated Upset plots
library(ComplexHeatmap) # For heatmaps
library(circlize) # For heatmap color customization
library(httr)
library(jsonlite)
library(readxl)
library(data.table)
library(DESeq2)
## -----------------------------------------------------------------------------------------------------------------------------------------
#
# Define output directory
output_dir <- "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap"

input_conditions_paths <- list(
  # WT as reference
  Ins_WT_HippoUP = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/HippovsWT.up.txt",
  Ins_WT_YorkieUP = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/YorkievsWT.up.txt",
  Ins_WT_WartsUP = "W:/mansel/RNAseq/For_publication/Ko_Insertion/output_R_analyses/tables/WartsvsWT.up.txt",
  
  Ins_WT_HippoDOWN = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/HippovsWT.down.txt",
  Ins_WT_YorkieDOWN = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/YorkievsWT.down.txt",
  Ins_WT_WartsDOWN = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/WartsvsWT.down.txt",
  Ins_WT_Neg2UP = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/Neg2vsWT.down.txt", # Neg2 vs WT, so Neg2 is the condition here
  Ins_WT_Neg2DOWN = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/Neg2vsWT.up.txt", # Neg2 vs WT, so Neg2 is the condition here
  
  # Neg2 as reference
  Ins_Neg2_HippoDOWN = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/Neg2vsHippo.up.txt", # Hippo vs Neg2, so Hippo is the condition
  Ins_Neg2_YorkieDOWN = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/YorkievsNeg2.down.txt",
  Ins_Neg2_WartsDOWN = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/WartsvsNeg2.down.txt",
  Ins_Neg2_HippoUP = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/Neg2vsHippo.down.txt", # Hippo vs Neg2, so Hippo is the condition
  Ins_Neg2_YorkieUP = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/YorkievsNeg2.up.txt",
  Ins_Neg2_WartsUP = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/tables/WartsvsNeg2.up.txt"
  
)

# Define thresholds for DEG filtering
padj_threshold <- 0.05
log2foldchange_threshold <- 0

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
  df <- read.table(file_path, header = TRUE, sep = "\t")[, cols] %>%
    filter(padj < padj_thresh & abs(log2FoldChange) > log2fc_thresh) %>%
    mutate(Id = sub("gene:", "", Id))
  
  type_val <- "insertion" # Default
  if (startsWith(df_name_prefix, "Del")) {
    type_val <- "deletion"
  }
  
  comparison_val <- ""
  if (df_name_prefix %in% c("Ins_WT_Neg2UP", "Ins_WT_Neg2DOWN")) {
    comparison_val <- "WT_vs_Neg2" # For WT_vs_Neg2, comparison itself doesn't include UP/DOWN
  } else if (grepl("_WT_", df_name_prefix)) {
    comparison_val <- paste0(df_name_prefix %>% stringr::str_remove("^(Ins|Del)_WT_"), "_WT")
  } else if (grepl("_Neg2_", df_name_prefix)) {
    comparison_val <- paste0(df_name_prefix %>% stringr::str_remove("^(Ins|Del)_Neg2_"), "_Neg2")
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

## -----------------------------------------------------------------------------------------------------------------------------------------

#' Loads and combines all relevant DEG tables into a single dataframe for comprehensive filtering.
#'
#' @param conditions_list A named list of file paths for each condition.
#' @param padj_thresh Adjusted p-value threshold.
#' @param log2fc_thresh Log2 fold change threshold.
#' @return A single dataframe containing all processed DEGs needed for suspicious DEG detection.
load_all_deg_data_for_filtering <- function(conditions_list, padj_thresh, log2fc_thresh) {
  processed_dfs_list <- list()
  for (cond_name in names(conditions_list)) {
    cols_to_select <- if (startsWith(cond_name, "Del")) {
      c(1, 38, 39, 42)
    } else if (startsWith(cond_name, "Ins")) {
      c(1, 112, 113, 116)
    } else {
      stop(paste("Unknown condition prefix for:", cond_name))
    }
    
    # Process all files relevant for suspicious DEG detection
    # This includes all KO_vs_WT, KO_vs_Neg2, and WT_vs_Neg2 comparisons
    if (grepl("^(Ins|Del)_(WT|Neg2)_(Yorkie|Warts|Hippo)(UP|DOWN)$", cond_name) ||
        grepl("Ins_WT_Neg2(UP|DOWN)", cond_name)) {
      processed_dfs_list[[cond_name]] <- process_deg_data(
        file_path = conditions_list[[cond_name]],
        cols = cols_to_select,
        padj_thresh = padj_thresh,
        log2fc_thresh = log2foldchange_threshold,
        df_name_prefix = cond_name
      )
    }
  }
  combined_results <- bind_rows(processed_dfs_list)
  return(combined_results)
}

## -----------------------------------------------------------------------------------------------------------------------------------------

#' Prepares data for Upset plots from a combined DEG dataframe.
#' Adds 'genotype' and 'direction' columns.
#'
#' @param combined_df A dataframe with processed DEG data.
#' @return A list containing:
#'      - significant_genes_df: Filtered and mutated dataframe for Upset (now with genotype and direction columns).
#'      - binary_matrix: A binary matrix indicating gene set membership.
#'      - all_set_names: A vector of all unique set names.
prepare_upset_data <- function(combined_df) {
  significant_genes_df <- combined_df %>%
    mutate(
      # Derive direction based on comparison name or log2FoldChange for WT_vs_Neg2
      direction = case_when(
        grepl("UP", comparison) ~ "UP",
        grepl("DOWN", comparison) ~ "DOWN",
        comparison == "WT_vs_Neg2" & log2FoldChange > 0 ~ "UP", # For WT_vs_Neg2, determine direction from log2FoldChange
        comparison == "WT_vs_Neg2" & log2FoldChange < 0 ~ "DOWN",
        TRUE ~ NA_character_ # Should not happen for valid DEG
      ),
      # Extract genotype from comparison string
      genotype = case_when(
        grepl("Yorkie", comparison) ~ "Yorkie",
        grepl("Warts", comparison) ~ "Warts",
        grepl("Hippo", comparison) ~ "Hippo",
        comparison == "WT_vs_Neg2" ~ "Neg2_vs_WT", # Assign a sensible genotype for this comparison
        TRUE ~ NA_character_
      ),
      # Create a unique condition_type for each set in Upset plot
      condition_type = case_when(
        genotype == "Neg2_vs_WT" ~ paste0(comparison, "_", direction), # WT_vs_Neg2_UP or WT_vs_Neg2_DOWN
        TRUE ~ paste(comparison, type, sep = "_") # e.g., YorkieUP_WT_insertion
      )
    )
  
  
  all_unique_genes <- unique(significant_genes_df$gene_id)
  plot_data <- data.frame(gene_id = all_unique_genes)
  all_set_names <- unique(significant_genes_df$condition_type)
  
  for (set_name in all_set_names) {
    plot_data[[set_name]] <- plot_data$gene_id %in%
      significant_genes_df$gene_id[significant_genes_df$condition_type == set_name]
  }
  
  binary_matrix <- significant_genes_df %>%
    mutate(value = TRUE) %>%
    pivot_wider(
      id_cols = gene_id,
      names_from = condition_type,
      values_from = value,
      values_fill = FALSE
    )
  
  return(list(
    significant_genes_df = significant_genes_df,
    binary_matrix = binary_matrix,
    all_set_names = all_set_names
  ))
}

#' Filters and saves specific DEG lists based on consistency and comparisons,
#' including "suspicious" DEG, filtered consistent DEG, and "low-filtered" DEG.
#'
#' @param combined_df The comprehensive DEG dataframe containing all KO_vs_WT, KO_vs_Neg2,
#'                    and WT_vs_Neg2 comparisons.
#' @param output_dir Directory to save the output files.
#' @return A dataframe with all DEGs after removing suspicious DEGs.
save_filtered_deg_lists <- function(combined_df, output_dir) {
  
  # Ensure the output directories exist
  homer_output_dir <- file.path(output_dir, "homer")
  deg_output_dir <- file.path(output_dir, "filtered_deg_lists")
  low_filtered_deg_output_dir <- file.path(output_dir, "low_filtered_deg_lists") # New directory for low_filtered files
  
  if (!dir.exists(homer_output_dir)) {
    dir.create(homer_output_dir, recursive = TRUE)
  }
  if (!dir.exists(deg_output_dir)) {
    dir.create(deg_output_dir, recursive = TRUE)
  }
  if (!dir.exists(low_filtered_deg_output_dir)) { # Create new directory
    dir.create(low_filtered_deg_output_dir, recursive = TRUE)
  }
  
  message("\n--- Calculating Suspicious DEGs ---")
  
  # 1. Get all DEGs from WT_vs_Neg2
  wt_neg2_deg_genes <- combined_df %>%
    filter(comparison == "WT_vs_Neg2") %>%
    pull(gene_id) %>% unique()
  
  # 2. Identify KO_vs_WT genes that are NOT in KO_vs_Neg2 (for the same direction and type)
  # Filter for KO_vs_WT and KO_vs_Neg2 comparisons
  ko_deg_comparisons <- combined_df %>%
    filter(genotype %in% c("Hippo", "Warts", "Yorkie")) # Already filtered for valid genotypes by import function
  
  suspicious_ko_wt_genes <- c()
  for (genotype_val in c("Hippo", "Warts", "Yorkie")) {
    for (direction in c("UP", "DOWN")) {
      for (type_val in c("insertion", "deletion")) {
        if (genotype_val == "Hippo" && type_val == "deletion") {
          next # Skip deletion for Hippo, as it's not present
        }
        
        # Construct the specific comparison names
        comp_wt_full <- paste0(genotype_val, direction, "_WT")
        comp_neg2_full <- paste0(genotype_val, direction, "_Neg2")
        
        # Get genes for KO_vs_WT for this specific direction/genotype/type
        genes_in_ko_wt_specific <- ko_deg_comparisons %>%
          filter(comparison == comp_wt_full, type == type_val) %>%
          pull(gene_id)
        
        # Get genes for KO_vs_Neg2 for this specific direction/genotype/type
        genes_in_ko_neg2_specific <- ko_deg_comparisons %>%
          filter(comparison == comp_neg2_full, type == "insertion") %>% # Only insertion because no Neg2 deletion
          pull(gene_id)
        
        # Genes in KO_vs_WT but NOT in KO_vs_Neg2 (for the same direction and type)
        absent_from_neg2 <- setdiff(genes_in_ko_wt_specific, genes_in_ko_neg2_specific)
        suspicious_ko_wt_genes <- c(suspicious_ko_wt_genes, absent_from_neg2)
      }
    }
  }
  
  # Combine all suspicious genes for the main filter
  suspicious_deg_genes <- unique(c(wt_neg2_deg_genes, suspicious_ko_wt_genes))
  
  if (length(suspicious_deg_genes) > 0) {
    write.table(suspicious_deg_genes, file.path(deg_output_dir, "Suspicious_DEG_genes.txt"),
                sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    message(paste0("Saved ", length(suspicious_deg_genes), " Suspicious_DEG_genes.txt"))
  } else {
    message("No Suspicious DEGs identified.")
  }
  
  # Filter the original combined_df to remove *all* suspicious DEGs
  filtered_combined_df <- combined_df %>%
    filter(!(gene_id %in% suspicious_deg_genes))
  
  # Create a "low filtered" dataframe by removing ONLY WT_vs_Neg2 DEGs
  low_filtered_combined_df <- combined_df %>%
    filter(!(gene_id %in% wt_neg2_deg_genes))
  
  message("\n--- Saving all DEG lists after removing Suspicious DEGs ---")
  
  # Define the core genotypes for iteration
  core_genotypes <- c("Hippo", "Warts", "Yorkie")
  
  # Loop through genotypes, references, directions, and types to save filtered all-DEG lists
  for (genotype in core_genotypes) {
    for (ref_type in c("WT", "Neg2")) {
      for (direction in c("UP", "DOWN")) {
        for (del_or_ins in c("deletion", "insertion")) {
          # Handle Hippo specific case (only insertion)
          if (genotype == "Hippo" && del_or_ins == "deletion") {
            next
          }
          
          # Construct the specific comparison name as it appears in the 'comparison' column
          current_comparison_name <- paste0(genotype, direction, "_", ref_type) # e.g., YorkieUP_WT
          
          # Save the fully filtered DEG lists
          deg_subset <- filtered_combined_df %>%
            filter(comparison == current_comparison_name, type == del_or_ins)
          
          if (nrow(deg_subset) > 0) {
            filename_all_deg <- file.path(deg_output_dir, paste0("DEG_", genotype, "_", direction, "_", ref_type, "_", del_or_ins, "_filtered.txt"))
            write.table(deg_subset, filename_all_deg, row.names = FALSE,quote = FALSE,sep = "\t")
            message("Saved ", basename(filename_all_deg))
          } else {
            message("No filtered DEG found for ", genotype, " ", direction, " ", ref_type, " (", del_or_ins, ")")
          }
          
          # Save the low-filtered DEG lists
          low_deg_subset <- low_filtered_combined_df %>%
            filter(comparison == current_comparison_name, type == del_or_ins)
          
          if (nrow(low_deg_subset) > 0) {
            filename_low_deg <- file.path(low_filtered_deg_output_dir, paste0("DEG_", genotype, "_", direction, "_", ref_type, "_", del_or_ins, "_low_filtered.txt"))
            write.table(low_deg_subset, filename_low_deg, row.names = FALSE, quote = FALSE, sep = "\t")
            message("Saved ", basename(filename_low_deg))
          } else {
            message("No low-filtered DEG found for ", genotype, " ", direction, " ", ref_type, " (", del_or_ins, ")")
          }
        }
      }
    }
  }
  
  # Save the WT_vs_Neg2 DEG list separately (as it's not a KO genotype and only 'insertion' type)
  # This file should remain the same as the "suspicious" WT_vs_Neg2 genes are filtered out
  # in both filtered_combined_df and low_filtered_combined_df.
  # For the "low_filtered" context, this file should be the original WT_vs_Neg2 DEG.
  wt_neg2_deg_original <- combined_df %>%
    filter(comparison == "WT_vs_Neg2")
  
  if (nrow(wt_neg2_deg_original) > 0) {
    # Save for the main filtered directory (it should technically be empty if WT_vs_Neg2 were all suspicious)
    # The original script saves it from combined_df, so keeping that logic for filtered.
    write.table(wt_neg2_deg_original, file.path(deg_output_dir, "DEG_WT_vs_Neg2_filtered.txt"), row.names = FALSE,quote = FALSE,sep = "\t")
    message("Saved DEG_WT_vs_Neg2_filtered.txt (using all suspicious filtering for consistency)")
    
    # Save the original WT_vs_Neg2 for the low_filtered directory (as it's the basis of the low filter)
    write.table(wt_neg2_deg_original, file.path(low_filtered_deg_output_dir, "DEG_WT_vs_Neg2_low_filtered.txt"), row.names = FALSE,quote = FALSE,sep = "\t")
    message("Saved DEG_WT_vs_Neg2_low_filtered.txt")
  } else {
    message("No DEG found for WT vs Neg2.")
  }
  
  message("\n--- Saving Consistent DEGs (shared between deletion and insertion, after removing Suspicious DEGs) ---")
  
  # Consistent DEGs (shared between insertion and deletion) after removing suspicious DEGs (full filter)
  for (genotype in core_genotypes) {
    if (genotype == "Hippo" | genotype == "Yorkie") { # no 'consistent' definition for Hippo and Yorkie
      next
    }
    
    for (ref_type in c("WT", "Neg2")) { # Loop through reference types
      current_comparison_base <- paste0(genotype, "_", ref_type) # e.g., Warts_WT
      
      consistent_genes_current_per_comp <- filtered_combined_df %>% # Uses fully filtered df
        filter(genotype == !!genotype, endsWith(comparison, paste0("_", ref_type))) %>% # Filter by ref_type
        group_by(gene_id, comparison) %>%
        summarize(
          has_deletion = any(type == "deletion"),
          has_insertion = any(type == "insertion"),
          .groups = "drop"
        ) %>%
        filter(has_deletion & has_insertion) %>%
        pull(gene_id) %>% unique()
      
      if (length(consistent_genes_current_per_comp) > 0) {
        consistent_deg_df <- filtered_combined_df %>% # Uses fully filtered df
          filter(gene_id %in% consistent_genes_current_per_comp,
                 genotype == !!genotype,
                 endsWith(comparison, paste0("_", ref_type))) %>% # Filter by ref_type
          distinct(gene_id, comparison, type, .keep_all = TRUE) # Ensure unique rows
        
        filename_consistent_deg <- file.path(deg_output_dir, paste0("Consistent_DEG_", genotype, "_", ref_type, "_filtered.txt"))
        write.table(consistent_deg_df, filename_consistent_deg, row.names = FALSE,quote = FALSE,sep = "\t")
        message("Saved ", basename(filename_consistent_deg))
        
        # Also save for Homer2
        write.table(consistent_genes_current_per_comp, file.path(homer_output_dir, paste0(genotype, "_", ref_type, "_filtered.txt")),
                    sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
        message("Saved homer/", genotype, "_", ref_type, "_filtered.txt (consistent, after removing suspicious DEGs)")
      } else {
        message("No consistent DEGs found for ", genotype, " ", ref_type, " after removing suspicious DEGs.")
      }
    }
  }
  
  message("\n--- Saving Consistent DEGs (shared between deletion and insertion, after removing ONLY WT_vs_Neg2 DEGs) ---")
  
  # Consistent DEGs (shared between insertion and deletion) after removing ONLY WT_vs_Neg2 DEGs (low filter)
  for (genotype in core_genotypes) {
    if (genotype == "Hippo" | genotype == "Yorkie") { # no 'consistent' definition for Hippo and Yorkie
      next
    }
    
    for (ref_type in c("WT", "Neg2")) { # Loop through reference types
      current_comparison_base <- paste0(genotype, "_", ref_type) # e.g., Warts_WT
      
      consistent_genes_current_per_comp_low <- low_filtered_combined_df %>% # Uses low filtered df
        filter(genotype == !!genotype, endsWith(comparison, paste0("_", ref_type))) %>% # Filter by ref_type
        group_by(gene_id, comparison) %>%
        summarize(
          has_deletion = any(type == "deletion"),
          has_insertion = any(type == "insertion"),
          .groups = "drop"
        ) %>%
        filter(has_deletion & has_insertion) %>%
        pull(gene_id) %>% unique()
      
      if (length(consistent_genes_current_per_comp_low) > 0) {
        consistent_deg_df_low <- low_filtered_combined_df %>% # Uses low filtered df
          filter(gene_id %in% consistent_genes_current_per_comp_low,
                 genotype == !!genotype,
                 endsWith(comparison, paste0("_", ref_type))) %>% # Filter by ref_type
          distinct(gene_id, comparison, type, .keep_all = TRUE)
        
        filename_consistent_deg_low <- file.path(low_filtered_deg_output_dir, paste0("Consistent_DEG_", genotype, "_", ref_type, "_low_filtered.txt"))
        write.table(consistent_deg_df_low, filename_consistent_deg_low, row.names = FALSE,quote = FALSE,sep = "\t")
        message("Saved ", basename(filename_consistent_deg_low))
      } else {
        message("No consistent DEGs found for ", genotype, " ", ref_type, " after removing only WT_vs_Neg2 DEGs.")
      }
    }
  }
  
  
  # For Warts, save only insertion genes after removing suspicious DEGs (full filter)
  for (ref_type in c("WT", "Neg2")) {
    Warts_insertion_genes_filtered <- filtered_combined_df %>%
      filter(genotype == "Warts" & type == "insertion" & endsWith(comparison, paste0("_", ref_type))) %>%
      pull(gene_id) %>% unique()
    
    if(length(Warts_insertion_genes_filtered) > 0) {
      Warts_insertion_df <- filtered_combined_df %>%
        filter(gene_id %in% Warts_insertion_genes_filtered,
               genotype == "Warts",
               type == "insertion",
               endsWith(comparison, paste0("_", ref_type))) %>%
        distinct(gene_id, comparison, type, .keep_all = TRUE)
      
      filename_Warts_insertion <- file.path(deg_output_dir, paste0("Warts_insertion_", ref_type, "_filtered.txt"))
      write.table(Warts_insertion_df, filename_Warts_insertion, row.names = FALSE,quote = FALSE,sep = "\t")
      message("Saved ", basename(filename_Warts_insertion))
      
      write.table(Warts_insertion_genes_filtered, file.path(homer_output_dir, paste0("Warts_", ref_type, "_filtered.txt")),
                  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
      message("Saved homer/Warts_", ref_type, "_filtered.txt (insertion only, after removing suspicious DEGs)")
    } else {
      message("No Warts insertion DEGs found for ", ref_type, " after removing suspicious DEGs.")
    }
    
    # For Warts, save only insertion genes after removing ONLY WT_vs_Neg2 DEGs (low filter)
    Warts_insertion_genes_low_filtered <- low_filtered_combined_df %>%
      filter(genotype == "Warts" & type == "insertion" & endsWith(comparison, paste0("_", ref_type))) %>%
      pull(gene_id) %>% unique()
    
    if(length(Warts_insertion_genes_low_filtered) > 0) {
      Warts_insertion_df_low <- low_filtered_combined_df %>%
        filter(gene_id %in% Warts_insertion_genes_low_filtered,
               genotype == "Warts",
               type == "insertion",
               endsWith(comparison, paste0("_", ref_type))) %>%
        distinct(gene_id, comparison, type, .keep_all = TRUE)
      
      filename_Warts_insertion_low <- file.path(low_filtered_deg_output_dir, paste0("Warts_insertion_", ref_type, "_low_filtered.txt"))
      write.table(Warts_insertion_df_low, filename_Warts_insertion_low, row.names = FALSE,quote = FALSE,sep = "\t")
      message("Saved ", basename(filename_Warts_insertion_low))
      
      write.table(Warts_insertion_genes_low_filtered, file.path(homer_output_dir, paste0("Warts_", ref_type, "_low_filtered.txt")),
                  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
      message("Saved homer/Warts_", ref_type, "_low_filtered.txt (insertion only, after removing suspicious DEGs)")
    } else {
      message("No Warts insertion DEGs found for ", ref_type, " after removing only WT_vs_Neg2 DEGs.")
    }
  }
  
  
  # For Yorkie, save only insertion genes after removing suspicious DEGs (full filter)
  for (ref_type in c("WT", "Neg2")) {
    Yorkie_insertion_genes_filtered <- filtered_combined_df %>%
      filter(genotype == "Yorkie" & type == "insertion" & endsWith(comparison, paste0("_", ref_type))) %>%
      pull(gene_id) %>% unique()
    
    if(length(Yorkie_insertion_genes_filtered) > 0) {
      Yorkie_insertion_df <- filtered_combined_df %>%
        filter(gene_id %in% Yorkie_insertion_genes_filtered,
               genotype == "Yorkie",
               type == "insertion",
               endsWith(comparison, paste0("_", ref_type))) %>%
        distinct(gene_id, comparison, type, .keep_all = TRUE)
      
      filename_Yorkie_insertion <- file.path(deg_output_dir, paste0("Yorkie_insertion_", ref_type, "_filtered.txt"))
      write.table(Yorkie_insertion_df, filename_Yorkie_insertion, row.names = FALSE,quote = FALSE,sep = "\t")
      message("Saved ", basename(filename_Yorkie_insertion))
      
      write.table(Yorkie_insertion_genes_filtered, file.path(homer_output_dir, paste0("Yorkie_", ref_type, "_filtered.txt")),
                  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
      message("Saved homer/Yorkie_", ref_type, "_filtered.txt (insertion only, after removing suspicious DEGs)")
    } else {
      message("No Yorkie insertion DEGs found for ", ref_type, " after removing suspicious DEGs.")
    }
    
    # For Yorkie, save only insertion genes after removing ONLY WT_vs_Neg2 DEGs (low filter)
    Yorkie_insertion_genes_low_filtered <- low_filtered_combined_df %>%
      filter(genotype == "Yorkie" & type == "insertion" & endsWith(comparison, paste0("_", ref_type))) %>%
      pull(gene_id) %>% unique()
    
    if(length(Yorkie_insertion_genes_low_filtered) > 0) {
      Yorkie_insertion_df_low <- low_filtered_combined_df %>%
        filter(gene_id %in% Yorkie_insertion_genes_low_filtered,
               genotype == "Yorkie",
               type == "insertion",
               endsWith(comparison, paste0("_", ref_type))) %>%
        distinct(gene_id, comparison, type, .keep_all = TRUE)
      
      filename_Yorkie_insertion_low <- file.path(low_filtered_deg_output_dir, paste0("Yorkie_insertion_", ref_type, "_low_filtered.txt"))
      write.table(Yorkie_insertion_df_low, filename_Yorkie_insertion_low, row.names = FALSE,quote = FALSE,sep = "\t")
      message("Saved ", basename(filename_Yorkie_insertion_low))
      
      write.table(Yorkie_insertion_genes_low_filtered, file.path(homer_output_dir, paste0("Yorkie_", ref_type, "_low_filtered.txt")),
                  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
      message("Saved homer/Yorkie_", ref_type, "_low_filtered.txt (insertion only, after removing suspicious DEGs)")
    } else {
      message("No Yorkie insertion DEGs found for ", ref_type, " after removing only WT_vs_Neg2 DEGs.")
    }
  }
  
  
  # For Hippo, save only insertion genes (as it has no deletion for consistency) after removing suspicious DEGs (full filter)
  for (ref_type in c("WT", "Neg2")) {
    hippo_insertion_genes_filtered <- filtered_combined_df %>%
      filter(genotype == "Hippo" & type == "insertion" & endsWith(comparison, paste0("_", ref_type))) %>%
      pull(gene_id) %>% unique()
    
    if(length(hippo_insertion_genes_filtered) > 0) {
      hippo_insertion_df <- filtered_combined_df %>%
        filter(gene_id %in% hippo_insertion_genes_filtered,
               genotype == "Hippo",
               type == "insertion",
               endsWith(comparison, paste0("_", ref_type))) %>%
        distinct(gene_id, comparison, type, .keep_all = TRUE)
      
      filename_hippo_insertion <- file.path(deg_output_dir, paste0("Hippo_insertion_", ref_type, "_filtered.txt"))
      write.table(hippo_insertion_df, filename_hippo_insertion, row.names = FALSE,quote = FALSE,sep = "\t")
      message("Saved ", basename(filename_hippo_insertion))
      
      write.table(hippo_insertion_genes_filtered, file.path(homer_output_dir, paste0("Hippo_", ref_type, "_filtered.txt")),
                  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
      message("Saved homer/Hippo_", ref_type, "_filtered.txt (insertion only, after removing suspicious DEGs)")
    } else {
      message("No Hippo insertion DEGs found for ", ref_type, " after removing suspicious DEGs.")
    }
    
    # For Hippo, save only insertion genes (as it has no deletion for consistency) after removing ONLY WT_vs_Neg2 DEGs (low filter)
    hippo_insertion_genes_low_filtered <- low_filtered_combined_df %>%
      filter(genotype == "Hippo" & type == "insertion" & endsWith(comparison, paste0("_", ref_type))) %>%
      pull(gene_id) %>% unique()
    
    if(length(hippo_insertion_genes_low_filtered) > 0) {
      hippo_insertion_df_low <- low_filtered_combined_df %>%
        filter(gene_id %in% hippo_insertion_genes_low_filtered,
               genotype == "Hippo",
               type == "insertion",
               endsWith(comparison, paste0("_", ref_type))) %>%
        distinct(gene_id, comparison, type, .keep_all = TRUE)
      
      filename_hippo_insertion_low <- file.path(low_filtered_deg_output_dir, paste0("Hippo_insertion_", ref_type, "_low_filtered.txt"))
      write.table(hippo_insertion_df_low, filename_hippo_insertion_low, row.names = FALSE,quote = FALSE,sep = "\t")
      message("Saved ", basename(filename_hippo_insertion_low))
      
      write.table(hippo_insertion_genes_low_filtered, file.path(homer_output_dir, paste0("Hippo_", ref_type, "_low_filtered.txt")),
                  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
      message("Saved homer/Hippo_", ref_type, "_low_filtered.txt (insertion only, after removing suspicious DEGs)")
    } else {
      message("No Hippo insertion DEGs found for ", ref_type, " after removing only WT_vs_Neg2 DEGs.")
    }
  }
  
  return(filtered_combined_df) # Return the fully filtered dataframe for plotting
}
## -----------------------------------------------------------------------------------------------------------------------------------------
# --- Main Execution Loop ---

# 1. Load ALL relevant raw DEG data needed for the *entire* filtering process ONCE.
# This 'all_raw_deg_data_combined' will contain KO_vs_WT, KO_vs_Neg2, and WT_vs_Neg2 comparisons.
all_raw_deg_data_combined <- load_all_deg_data_for_filtering(
  conditions_list = input_conditions_paths,
  padj_thresh = padj_threshold,
  log2fc_thresh = log2foldchange_threshold
)

# Check if any data was loaded at all
if (nrow(all_raw_deg_data_combined) == 0) {
  stop("No DEG data loaded from input_conditions_paths. Please check paths and file contents.")
}

# 2. Add 'genotype' and 'direction' columns to the comprehensive raw data.
# This step is critical before calculating suspicious DEGs.
# The 'prepare_upset_data' function is adapted to also prepare this intermediate dataframe.
# We call it here to get the 'genotype' and 'direction' columns which are crucial for `save_filtered_deg_lists`.
temp_prepared_all_data <- prepare_upset_data(all_raw_deg_data_combined)$significant_genes_df

# 3. Perform filtering (suspicious DEG removal) and save all filtered lists.
# This function returns the fully filtered dataframe which will be used for plotting.
filtered_deg_for_plotting <- save_filtered_deg_lists(
  combined_df = temp_prepared_all_data, # Pass the comprehensive data with genotype/direction
  output_dir = output_dir
)

## -----------------------------------------------------------------------------------------------------------------------------------------
#' Retrieves protein information and sequences from UniProtKB or UniParc.
#'
#' This function queries the UniProt API to convert a list of input IDs
#' (e.g., gene IDs) into UniProt primary IDs, gene symbols, RefSeq cross-references,
#' protein names, and sequences. If an ID is not found in UniProtKB, it attempts
#' to search UniParc.
#'
#' @param input_ids A character vector of IDs to convert.
#' @param Organism The NCBI Taxonomy ID for the organism (e.g., "946362" for Salpingoeca rosetta).
#' @param Gene_name Logical; if TRUE, the input IDs are treated as gene names/symbols for the query.
#' @return A list containing three data frames:
#'    - `Gene_primary_Results`: Data from UniProtKB for found IDs.
#'    - `Uniparc_Results`: Data from UniParc for IDs not found in UniProtKB but found in UniParc.
#'    - `Not_Found_Ids`: A character vector of input IDs that were not found in either database.
get_conversion_ids <- function(input_ids, Organism, Gene_name) {
  results_df <- data.frame(
    Input_ID = character(), Gene_primary_ID = character(),
    Symbol_ID = character(), Xref_refseq_ID = character(),
    sequence = character(), Protein_name = character(), stringsAsFactors = FALSE
  )
  
  uniparc_df <- data.frame(
    Input_ID = character(), Uniparc_ID = character(), sequence = character(), stringsAsFactors = FALSE
  )
  
  not_found_ids <- c()
  
  for (Input_id in input_ids) {
    # Try UniProtKB first
    base_url_uniprotkb <- "https://rest.uniprot.org/uniprotkb/search"
    query_field <- if (Gene_name) "?query=gene_exact:" else "?query="
    query_url_uniprotkb <- paste0(
      base_url_uniprotkb, query_field,  Input_id,
      "+AND+organism_id:", Organism,
      "&fields=sequence,gene_primary,accession,xref_refseq,gene_orf,protein_name"
    )
    
    response_uniprotkb <- GET(query_url_uniprotkb, accept_json())
    
    if (status_code(response_uniprotkb) == 200) {
      response_content_uniprotkb <- httr::content(response_uniprotkb, as = "text")
      response_json_uniprotkb <- fromJSON(response_content_uniprotkb, flatten = TRUE)
      
      if (length(response_json_uniprotkb$results) > 0) {
        first_result <- response_json_uniprotkb$results[1, ]
        
        seq <- first_result$sequence.value %||% NA_character_
        gene_primary <- first_result$primaryAccession %||% NA_character_
        
        symbol <- NA_character_
        if (!is.null(first_result$genes) && length(first_result$genes[[1]]$geneName.value) > 0) {
          symbol <- first_result$genes[[1]]$geneName.value[1]
        } else if (!is.null(first_result$genes) && length(first_result$genes[[1]]$orfNames[[1]]$value) > 0 && Organism == "946362") {
          symbol <- first_result$genes[[1]]$orfNames[[1]]$value[1]
        }
        
        refseq <- NA_character_
        if (!is.null(first_result$uniProtKBCrossReferences) && length(first_result$uniProtKBCrossReferences[[1]]$id) > 0) {
          refseq_xrefs <- first_result$uniProtKBCrossReferences[[1]] %>%
            filter(database == "RefSeq")
          if (nrow(refseq_xrefs) > 0) {
            refseq <- refseq_xrefs$id[1]
          }
        }
        
        protein_name <- first_result$proteinDescription.recommendedName.fullName.value %||% NA_character_
        
        results_df <- rbind(results_df, data.frame(
          Input_ID = Input_id, Gene_primary_ID = gene_primary,
          Symbol_ID = symbol, Xref_refseq_ID = refseq, Protein_name = protein_name,
          sequence = seq, stringsAsFactors = FALSE
        ))
      } else {
        message(paste("No results found in UniProtKB for:", Input_id, ". Testing UniParc..."))
        base_url_uniparc <- "https://rest.uniprot.org/uniparc/search"
        query_url_uniparc <- paste0(
          base_url_uniparc, "?query=", Input_id,
          "+AND+organism_id:", Organism,
          "&fields=sequence,length,uniparc_id"
        )
        
        response_uniparc <- GET(query_url_uniparc, accept_json())
        
        if (status_code(response_uniparc) == 200) {
          response_content_uniparc <- httr::content(response_uniparc, as = "text")
          response_json_uniparc <- fromJSON(response_content_uniparc, flatten = TRUE)
          
          if (length(response_json_uniparc$results) > 0) {
            longest_idx <- which.max(as.numeric(response_json_uniparc$results$sequence.length))
            seq_uniparc <- response_json_uniparc$results$sequence.value[[longest_idx]] %||% NA_character_
            id_uniparc <- response_json_uniparc$results$uniParcId[[longest_idx]] %||% NA_character_
            
            uniparc_df <- rbind(uniparc_df, data.frame(
              Input_ID = Input_id, Uniparc_ID = id_uniparc, sequence = seq_uniparc, stringsAsFactors = FALSE
            ))
          } else {
            message(paste("No results found in UniParc either for:", Input_id))
            not_found_ids <- c(not_found_ids, Input_id)
          }
        } else {
          warning(paste("Error retrieving data from UniParc for", Input_id, ": HTTP Status", status_code(response_uniparc)))
          not_found_ids <- c(not_found_ids, Input_id)
        }
      }
    } else {
      warning(paste("Error retrieving data from UniProtKB for", Input_id, ": HTTP Status", status_code(response_uniprotkb), ". Testing UniParc..."))
      base_url_uniparc <- "https://rest.uniprot.org/uniparc/search"
      query_url_uniparc <- paste0(
        base_url_uniparc, "?query=", Input_id,
        "+AND+organism_id:", Organism,
        "&fields=sequence,length,uniparc_id"
      )
      
      response_uniparc <- GET(query_url_uniparc, accept_json())
      
      if (status_code(response_uniparc) == 200) {
        response_content_uniparc <- httr::content(response_uniparc, as = "text")
        response_json_uniparc <- fromJSON(response_content_uniparc, flatten = TRUE)
        
        if (length(response_json_uniparc$results) > 0) {
          longest_idx <- which.max(as.numeric(response_json_uniparc$results$sequence.length))
          seq_uniparc <- response_json_uniparc$results$sequence.value[[longest_idx]] %||% NA_character_
          id_uniparc <- response_json_uniparc$results$uniParcId[[longest_idx]] %||% NA_character_
          
          uniparc_df <- rbind(uniparc_df, data.frame(
            Input_ID = Input_id, Uniparc_ID = id_uniparc, sequence = seq_uniparc, stringsAsFactors = FALSE
          ))
        } else {
          message(paste("No results found in UniParc either for:", Input_id))
          not_found_ids <- c(not_found_ids, Input_id)
        }
      } else {
        warning(paste("Error retrieving data from UniParc for", Input_id, ": HTTP Status", status_code(response_uniparc)))
        not_found_ids <- c(not_found_ids, Input_id)
      }
    }
  }
  return(list(
    Gene_primary_Results = results_df,
    Uniparc_Results = uniparc_df,
    Not_Found_Ids = not_found_ids
  ))
}


#' Helper function to handle NULL values gracefully
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    y
  } else {
    x
  }
}

#' Writes protein sequences to a FASTA file.
#'
#' @param results_df A dataframe containing Input_ID, Gene_primary_ID (optional),
#'    Uniparc_ID (optional), and sequence.
#' @param output_file The path to the output FASTA file.
#' @param organism_abbr A short abbreviation for the organism (e.g., "Sros").
write_fasta <- function(results_df, output_file, organism_abbr) {
  fileConn <- file(output_file, open = "w")
  for (i in 1:nrow(results_df)) {
    header_id <- ifelse(!is.na(results_df$Gene_primary_ID[i]),
                        results_df$Gene_primary_ID[i],
                        results_df$Uniparc_ID[i]
    )
    sequence <- ifelse(!is.na(results_df$sequence.x[i]),
                       results_df$sequence.x[i],
                       results_df$sequence.y[i]
    )
    
    if (!is.na(header_id) && !is.na(sequence)) {
      fasta_header <- paste0(">", organism_abbr, "|", results_df$Input_ID[i], "_", header_id)
      writeLines(fasta_header, fileConn)
      writeLines(sequence, fileConn)
    } else {
      warning(paste("Skipping FASTA entry for Input_ID:", results_df$Input_ID[i], "due to missing ID or sequence."))
    }
  }
  close(fileConn)
}

#' Writes a list of IDs that were not found to a text file.
#'
#' @param not_found_ids A character vector of IDs that could not be converted.
#' @param output_file The path to the output text file.
write_not_found_ids <- function(not_found_ids, output_file) {
  if (length(not_found_ids) > 0) {
    writeLines(not_found_ids, output_file)
  } else {
    message(paste("No not-found IDs to write to:", output_file))
  }
}

## -----------------------------------------------------------------------------------------------------------------------------------------

#' Extracts top DEGs and prepares data for heatmap.
#'
#' @param conditions_list A named list of file paths to complete DEG tables.
#' @param num_top_genes Number of top genes to select (by padj).
#' @param padj_thresh Adjusted p-value threshold for initial filtering.
#' @param log2fc_thresh Log2 fold change threshold for initial filtering.
#' @param columns_to_select Columns to select from each complete DEG table.
#' @return A list containing `list_top_genes` (vector of gene IDs) and `combined_top_genes_df` (dataframe of top genes).
prepare_heatmap_data <- function(conditions_list, num_top_genes, padj_thresh, log2fc_thresh, columns_to_select) {
  all_top_genes_list <- list()
  for (condition in names(conditions_list)) {
    data <- read.table(conditions_list[[condition]], header = TRUE, sep = "\t")[, columns_to_select]
    colnames(data)[1] <- "Id"
    top_genes <- data %>%
      filter(padj < padj_thresh & abs(log2FoldChange) > log2fc_thresh) %>%
      distinct(Id, .keep_all = TRUE) %>%
      arrange(padj) %>%
      head(num_top_genes) # Number must be lower of equal to total number of gene present
    all_top_genes_list[[condition]] <- top_genes
  }
  
  # Determine the maximum length among the 'Id' columns
  max_len <- max(sapply(all_top_genes_list, function(x) length(x$Id)))
  
  # Pad the 'Id' columns with NA values to match the maximum length
  padded_ids_list <- lapply(all_top_genes_list, function(x) {
    ids <- x$Id
    length(ids) <- max_len # This pads with NA if 'ids' is shorter than max_len
    return(ids)
  })
  
  # Combine the padded 'Id' columns into a dataframe
  combined_top_genes_df <- as.data.frame(padded_ids_list)
  
  colnames(combined_top_genes_df) <- names(all_top_genes_list)
  
  # Use `na.omit` before `unique()` to ensure that only actual gene IDs are considered for the list_top_genes
  list_top_genes <- combined_top_genes_df %>%
    pivot_longer(cols = everything(), values_to = "Combined") %>%
    dplyr::select(Combined) %>%
    na.omit() %>% # Remove NA entries
    unique() %>%
    pull(Combined)
  
  return(list(list_top_genes = list_top_genes, combined_top_genes_df = combined_top_genes_df))
}


#' Generates and saves a ComplexHeatmap.
#'
#' @param mat The expression matrix (scaled).
#' @param anno The annotation dataframe for columns.
#' @param genotype_colors Named vector of colors for genotypes.
#' @param output_file Path to save the SVG plot.
#' @param row_km Number of row clusters (k-means).
#' @param column_km Number of column clusters (k-means).
#' @param row_font_size Font size for row names.
#' @param column_font_size Font size for column names.
#' @param heatmap_title Title for the heatmap legend.
#' @param heatmap_legend_direction Direction of the heatmap legend.
#' @param annotation_legend_direction Direction of the annotation legend.
plot_expression_heatmap <- function(mat, anno, genotype_colors, output_file,
                                    row_km, column_km, row_font_size, column_font_size,
                                    heatmap_title = "Scaled Expression", heatmap_legend_direction = "vertical",
                                    annotation_legend_direction = "vertical",width = width, height = height) {
  # 'primary_results' is expected to be available in the global environment
  # or passed as an argument if its scope is limited.
  # For now, assuming it's loaded globally as in the original script.
  label_map <- setNames(
    primary_results$combine,
    primary_results$V1
  )
  
  mat <- mat[, rownames(anno)]
  rownames(mat) <- gsub("^gene:", "", rownames(mat))
  rownames(mat) <- ifelse(rownames(mat) %in% names(label_map),
                          label_map[rownames(mat)],
                          rownames(mat)
  )
  
  # Scale the matrix by rows
  scaled_mat <- t(scale(t(mat)))
  
  # Define the color function
  col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("darkblue", "white", "gold"))
  
  # Define the column annotation
  col_annotation <- HeatmapAnnotation(
    genotype = anno$genotype,
    col = list(genotype = genotype_colors),
    border = TRUE,
    gp = gpar(col = "darkgray", lwd = 0.5),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
  
  cluster_rows_flag <- row_km > 0
  
  heatmap <- Heatmap(
    scaled_mat,
    name = "Expression",
    col = col_fun,
    cluster_rows = cluster_rows_flag,
    cluster_columns = column_km > 0,
    row_order = if (!cluster_rows_flag) rownames(scaled_mat) else NULL,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    row_km = row_km,
    column_km = column_km,
    row_names_gp = gpar(fontsize = row_font_size),
    column_names_gp = gpar(fontsize = column_font_size),
    rect_gp = gpar(col = "darkgray", lwd = 0.5),
    border = TRUE,
    top_annotation = col_annotation,
    row_names_max_width = unit(20, "cm"),
    heatmap_legend_param = list(
      title = heatmap_title,
      legend_height = unit(4, "cm"),
      direction = heatmap_legend_direction
    )
  )
  annotation_legend <- Legend(
    title = "Genotype",
    at = names(genotype_colors),
    labels = names(genotype_colors),
    legend_gp = gpar(fill = genotype_colors),
    direction = annotation_legend_direction
  )
  
  svg(output_file, width = width, height = height) # Adjust dimensions as needed
  draw(
    heatmap,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    annotation_legend_list = list(annotation_legend)
  )
  dev.off()
}


##-------------------- KO post strict filter vs WT (subset) -----------------------------------##

output_dir <- "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/"
# You might want to create the directory if it doesn't exist:
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

conditions_ins_heatmap <- list(
  WT_Warts = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/filtered_deg_lists/Warts_insertion_WT_filtered.txt",
  WT_Yorkie = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/filtered_deg_lists/Yorkie_insertion_WT_filtered.txt",
  WT_Hippo = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/filtered_deg_lists/Hippo_insertion_WT_filtered.txt"
)

dds <- readRDS("W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Ins_out_deseq2_dds.rds")

# Prepare data for deletion heatmap
heatmap_data_ins <- prepare_heatmap_data(
  conditions_list = conditions_ins_heatmap,
  num_top_genes = 50, # Original script selected 30
  padj_thresh = 0.05,
  log2fc_thresh = 0,
  columns_to_select = c(1:4) # Columns relevant for deletion data
)

# Filter to only keep some interesting genes
# Vector of genes you want to keep
genes_to_keep <- c(
  "PTSG_11789", "PTSG_10931", "PTSG_05399", "PTSG_09711", "PTSG_04943",
  "PTSG_12667", "PTSG_03144", "PTSG_01648", "PTSG_07368", "PTSG_03426",
  "PTSG_05995", "PTSG_11757", "PTSG_02708", "PTSG_08232", "PTSG_07978",
  "PTSG_09125", "PTSG_09127", "PTSG_04127", "PTSG_04961", "PTSG_01647"
)

# Filter rows where any column matches a gene in genes_to_keep
list_top_genes_ins <- heatmap_data_ins$list_top_genes
list_top_genes_ins <- genes_to_keep[genes_to_keep %in% list_top_genes_ins]
list_top_genes_ins <- gsub("PT","gene:PT",list_top_genes_ins)

combined_top_genes_ins <- heatmap_data_ins$combined_top_genes_df%>%
  filter(WT_Warts %in% genes_to_keep |
           WT_Yorkie %in% genes_to_keep |
           WT_Hippo %in% genes_to_keep)
combined_top_genes_ins <- gsub("PT","gene:PT",combined_top_genes_ins)

# Save all unique DEG (top 60) for deletion heatmap
write.table(sub("gene:", "", list_top_genes_ins),
            file = file.path(output_dir, "KO_WT_small_filtered-DEG.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)

# Run Uniprot API function to get uniprot long name with the gene present on the heatmap
input_ids_for_uniprot <- read.table(
  file.path(output_dir, "KO_WT_small_filtered-DEG.txt"),
  sep = "\t", header = FALSE
)$V1

organism_ncbi_id <- "946362" # Salpingoeca rosetta
organism_abbr_name <- "Sros"
input_is_gene_name <- FALSE # Set to TRUE if input_ids are gene symbols, FALSE if they are protein/gene IDs

output_conversion_prefix <- file.path(output_dir, "KO_WT_small_filtered-DEG_uniprot")

### Run UniProt Conversion ###
message(paste("Starting UniProt conversion for", length(input_ids_for_uniprot), "genes..."))
conversion_results <- get_conversion_ids(
  input_ids = input_ids_for_uniprot,
  Organism = organism_ncbi_id,
  Gene_name = input_is_gene_name
)

primary_results_uniprot <- conversion_results$Gene_primary_Results
uniparc_results_uniprot <- conversion_results$Uniparc_Results
not_found_ids_uniprot <- conversion_results$Not_Found_Ids

# Merge results from UniProtKB and UniParc
# Ensure column names for 'sequence' are distinct before merging if they came from different sources
# Renaming here to avoid conflict: sequence.x for uniprotkb, sequence.y for uniparc
# The write_fasta function will then pick the correct one
founded_results_merged <- merge(
  primary_results_uniprot,
  uniparc_results_uniprot,
  by = "Input_ID",
  all = TRUE,
  suffixes = c(".x", ".y") # Use .x for primary_results_uniprot and .y for uniparc_results_uniprot
)

# Save the merged conversion table
write.table(founded_results_merged,
            paste0(output_conversion_prefix, "_conversion_table.txt"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
)

primary_results <- read.table(paste0(output_conversion_prefix, "_conversion_table.txt"),sep="\t",quote = "")%>%
  mutate(combine = ifelse(is.na(V5), V1, paste0(V1, " (", V5, ")")))%>%
  dplyr::select(V1, combine)

# Extract normalized counts
vsd_ins <- vst(dds)
transformed_counts_ins <- assay(vsd_ins)
mat_ins <- transformed_counts_ins[list_top_genes_ins, ]

# Prepare annotation for insertion heatmap 
anno_ins <- as.data.frame(colData(dds)[, "genotype", drop = FALSE])
cond_ins <- c("WT", "Yorkie","Warts","Hippo") 
anno_ins <- anno_ins[anno_ins$genotype %in% cond_ins, , drop = FALSE]

# Modify genotype labels (without "algo" specific labels)
anno_ins$genotype <- gsub("WT", "WT", anno_ins$genotype)
anno_ins$genotype <- gsub("Warts", "Warts KO", anno_ins$genotype)
anno_ins$genotype <- gsub("Yorkie", "Yorkie KO", anno_ins$genotype)
anno_ins$genotype <- gsub("Hippo", "Hippo KO", anno_ins$genotype)
anno_ins$genotype <- as.factor(anno_ins$genotype)
anno_ins$genotype <- droplevels(anno_ins$genotype)

# Define custom colors for deletion genotypes
genotype_colors_ins <- c(
  "WT" = "blue",
  "Yorkie KO" = "orange",
  "Warts KO" = "red",
  "Hippo KO" = "#e632c2"
)

# Plot deletion heatmap
plot_expression_heatmap(
  mat = mat_ins,
  anno = anno_ins,
  genotype_colors = genotype_colors_ins,
  output_file = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/Heatmap_KO_WT_small_filtered_DEG.svg",
  row_km = 0,
  column_km = 0,
  row_font_size = 10,
  column_font_size = 10,
  width = 16, 
  height = 8
)

##-------------------- KO post strict filter vs WT (all) -----------------------------------##

output_dir <- "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/"
# You might want to create the directory if it doesn't exist:
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

conditions_ins_heatmap <- list(
  WT_Warts = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/filtered_deg_lists/Warts_insertion_WT_filtered.txt",
  WT_Yorkie = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/filtered_deg_lists/Yorkie_insertion_WT_filtered.txt",
  WT_Hippo = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/filtered_deg_lists/Hippo_insertion_WT_filtered.txt"
)

dds <- readRDS("W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Ins_out_deseq2_dds.rds")

# Prepare data for deletion heatmap
heatmap_data_ins <- prepare_heatmap_data(
  conditions_list = conditions_ins_heatmap,
  num_top_genes = 150, 
  padj_thresh = 0.05,
  log2fc_thresh = 0,
  columns_to_select = c(1:4) # Columns relevant for deletion data
)

list_top_genes_ins <- heatmap_data_ins$list_top_genes
list_top_genes_ins <- gsub("PT","gene:PT",list_top_genes_ins)
combined_top_genes_ins <- heatmap_data_ins$combined_top_genes_df
combined_top_genes_ins <- gsub("PT","gene:PT",combined_top_genes_ins)

combined_top_genes_ins <- heatmap_data_ins$combined_top_genes_df%>%
  filter(WT_Warts %in% genes_to_keep |
           WT_Yorkie %in% genes_to_keep |
           WT_Hippo %in% genes_to_keep)
combined_top_genes_ins <- gsub("PT","gene:PT",combined_top_genes_ins)

# Save all unique DEG (top 60) for deletion heatmap
write.table(sub("gene:", "", list_top_genes_ins),
            file = file.path(output_dir, "KO_WT_filtered-DEG.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)

# Run Uniprot API function to get uniprot long name with the gene present on the heatmap
input_ids_for_uniprot <- read.table(
  file.path(output_dir, "KO_WT_filtered-DEG.txt"),
  sep = "\t", header = FALSE
)$V1

organism_ncbi_id <- "946362" # Salpingoeca rosetta
organism_abbr_name <- "Sros"
input_is_gene_name <- FALSE # Set to TRUE if input_ids are gene symbols, FALSE if they are protein/gene IDs

output_conversion_prefix <- file.path(output_dir, "KO_WT_filtered-DEG_uniprot")

### Run UniProt Conversion ###
message(paste("Starting UniProt conversion for", length(input_ids_for_uniprot), "genes..."))
conversion_results <- get_conversion_ids(
  input_ids = input_ids_for_uniprot,
  Organism = organism_ncbi_id,
  Gene_name = input_is_gene_name
)

primary_results_uniprot <- conversion_results$Gene_primary_Results
uniparc_results_uniprot <- conversion_results$Uniparc_Results
not_found_ids_uniprot <- conversion_results$Not_Found_Ids

# Merge results from UniProtKB and UniParc
# Ensure column names for 'sequence' are distinct before merging if they came from different sources
# Renaming here to avoid conflict: sequence.x for uniprotkb, sequence.y for uniparc
# The write_fasta function will then pick the correct one
founded_results_merged <- merge(
  primary_results_uniprot,
  uniparc_results_uniprot,
  by = "Input_ID",
  all = TRUE,
  suffixes = c(".x", ".y") # Use .x for primary_results_uniprot and .y for uniparc_results_uniprot
)

# Save the merged conversion table
write.table(founded_results_merged,
            paste0(output_conversion_prefix, "_conversion_table.txt"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
)

primary_results <- read.table(paste0(output_conversion_prefix, "_conversion_table.txt"),sep="\t",quote = "")%>%
  mutate(combine = ifelse(is.na(V5), V1, paste0(V1, " (", V5, ")")))%>%
  dplyr::select(V1, combine)

# Extract normalized counts
vsd_ins <- vst(dds)
transformed_counts_ins <- assay(vsd_ins)
mat_ins <- transformed_counts_ins[list_top_genes_ins, ]

# Prepare annotation for insertion heatmap 
anno_ins <- as.data.frame(colData(dds)[, "genotype", drop = FALSE])
cond_ins <- c("WT", "Yorkie","Warts","Hippo") 
anno_ins <- anno_ins[anno_ins$genotype %in% cond_ins, , drop = FALSE]

# Modify genotype labels (without "algo" specific labels)
anno_ins$genotype <- gsub("WT", "WT", anno_ins$genotype)
anno_ins$genotype <- gsub("Warts", "Warts KO", anno_ins$genotype)
anno_ins$genotype <- gsub("Yorkie", "Yorkie KO", anno_ins$genotype)
anno_ins$genotype <- gsub("Hippo", "Hippo KO", anno_ins$genotype)
anno_ins$genotype <- as.factor(anno_ins$genotype)
anno_ins$genotype <- droplevels(anno_ins$genotype)

# Define custom colors for deletion genotypes
genotype_colors_ins <- c(
  "WT" = "blue",
  "Yorkie KO" = "orange",
  "Warts KO" = "red",
  "Hippo KO" = "#e632c2"
)

# Plot deletion heatmap
plot_expression_heatmap(
  mat = mat_ins,
  anno = anno_ins,
  genotype_colors = genotype_colors_ins,
  output_file = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/Heatmap_KO_WT_filtered_DEG.svg",
  row_km = 2,
  column_km = 3,
  row_font_size = 10,
  column_font_size = 10,
  width = 16, 
  height = 35
)

##-------------------- Neg2 vs WT -----------------------------------##

output_dir <- "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/"
# You might want to create the directory if it doesn't exist:
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

conditions_ins_heatmap <- list(
  WT_Neg2 = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/filtered_deg_lists/DEG_WT_vs_Neg2_filtered.txt"
)

dds <- readRDS("W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Ins_out_deseq2_dds.rds")

# Prepare data for deletion heatmap
heatmap_data_ins <- prepare_heatmap_data(
  conditions_list = conditions_ins_heatmap,
  num_top_genes = 50, # Original script selected 30
  padj_thresh = 0.05,
  log2fc_thresh = 0,
  columns_to_select = c(1:4) # Columns relevant for deletion data
)

list_top_genes_ins <- heatmap_data_ins$list_top_genes
list_top_genes_ins <- gsub("PT","gene:PT",list_top_genes_ins)
combined_top_genes_ins <- heatmap_data_ins$combined_top_genes_df
combined_top_genes_ins <- gsub("PT","gene:PT",combined_top_genes_ins)

# Save all unique DEG (top 60) for deletion heatmap
write.table(sub("gene:", "", list_top_genes_ins),
            file = file.path(output_dir, "Neg2_WT-DEG.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)

# Run Uniprot API function to get uniprot long name with the gene present on the heatmap
input_ids_for_uniprot <- read.table(
  file.path(output_dir, "Neg2_WT-DEG.txt"),
  sep = "\t", header = FALSE
)$V1

organism_ncbi_id <- "946362" # Salpingoeca rosetta
organism_abbr_name <- "Sros"
input_is_gene_name <- FALSE # Set to TRUE if input_ids are gene symbols, FALSE if they are protein/gene IDs

output_conversion_prefix <- file.path(output_dir, "Neg2_WT-DEG_uniprot")

### Run UniProt Conversion ###
message(paste("Starting UniProt conversion for", length(input_ids_for_uniprot), "genes..."))
conversion_results <- get_conversion_ids(
  input_ids = input_ids_for_uniprot,
  Organism = organism_ncbi_id,
  Gene_name = input_is_gene_name
)

primary_results_uniprot <- conversion_results$Gene_primary_Results
uniparc_results_uniprot <- conversion_results$Uniparc_Results
not_found_ids_uniprot <- conversion_results$Not_Found_Ids

# Merge results from UniProtKB and UniParc
# Ensure column names for 'sequence' are distinct before merging if they came from different sources
# Renaming here to avoid conflict: sequence.x for uniprotkb, sequence.y for uniparc
# The write_fasta function will then pick the correct one
founded_results_merged <- merge(
  primary_results_uniprot,
  uniparc_results_uniprot,
  by = "Input_ID",
  all = TRUE,
  suffixes = c(".x", ".y") # Use .x for primary_results_uniprot and .y for uniparc_results_uniprot
)

# Save the merged conversion table
write.table(founded_results_merged,
            paste0(output_conversion_prefix, "_conversion_table.txt"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
)

primary_results <- read.table(paste0(output_conversion_prefix, "_conversion_table.txt"),sep="\t",quote = "")%>%
  mutate(combine = ifelse(is.na(V5), V1, paste0(V1, " (", V5, ")")))%>%
  dplyr::select(V1, combine)

# Extract normalized counts
vsd_ins <- vst(dds)
transformed_counts_ins <- assay(vsd_ins)
mat_ins <- transformed_counts_ins[list_top_genes_ins, ]

# Prepare annotation for insertion heatmap 
anno_ins <- as.data.frame(colData(dds)[, "genotype", drop = FALSE])
cond_ins <- c("WT", "Neg2") 
anno_ins <- anno_ins[anno_ins$genotype %in% cond_ins, , drop = FALSE]

# Modify genotype labels (without "algo" specific labels)
anno_ins$genotype <- as.factor(anno_ins$genotype)
anno_ins$genotype <- droplevels(anno_ins$genotype)

# Define custom colors for deletion genotypes
genotype_colors_ins <- c(
  "WT" = "blue",
  "Neg2" = "#86c5da"
)

# Plot deletion heatmap
plot_expression_heatmap(
  mat = mat_ins,
  anno = anno_ins,
  genotype_colors = genotype_colors_ins,
  output_file = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/Heatmap_Neg_WT_DEG.svg",
  row_km = 3,
  column_km = 2,
  row_font_size = 10,
  column_font_size = 10,
  height=12,
  width=10
)
