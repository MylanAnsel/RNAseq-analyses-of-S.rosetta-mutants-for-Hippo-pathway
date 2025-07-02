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
output_dir <- "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Heatmap"

input_conditions_paths <- list(
  # WT as reference
  Del_WT_YorkieUP = "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/tables/YCvsSr.up.txt",
  Del_WT_WartsUP = "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/tables/WBvsSr.up.txt",
  Del_WT_YorkieDOWN = "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/tables/YCvsSr.down.txt",
  Del_WT_WartsDOWN = "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/tables/WBvsSr.down.txt"
)

filter <- list(
  strict = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/filtered_deg_lists/Suspicious_DEG_genes.txt"
)

# Define the insertion filter file paths
insertion_filter <- list(
  Ins_strict_WT_Warts = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/filtered_deg_lists/Warts_insertion_WT_filtered.txt",
  Ins_strict_WT_Yorkie = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Heatmap/filtered_deg_lists/Yorkie_insertion_WT_filtered.txt"
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
  
  type_val <- "Deletion" # Default
  if (startsWith(df_name_prefix, "Del")) {
    type_val <- "deletion"
  }
  
  comparison_val <- ""
  if (grepl("_WT_", df_name_prefix)) {
    comparison_val <- paste0(df_name_prefix %>% stringr::str_remove("^(Ins|Del)_WT_"), "_WT")
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
#' This function does not filter by primary reference or comparison type during loading.
#'
#' @param conditions_list A named list of file paths for each condition.
#' @param padj_thresh Adjusted p-value threshold.
#' @param log2fc_thresh Log2 fold change threshold.
#' @return A single dataframe containing all processed DEGs.
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

## -----------------------------------------------------------------------------------------------------------------------------------------

#' Prepares data for Upset plots by adding genotype and direction.
#'
#' @param combined_deg_df A dataframe of combined DEG data.
#' @return A list containing the dataframe with significant genes and a list of all genes.
prepare_upset_data <- function(combined_deg_df) {
  # Add genotype and direction columns based on comparison name
  significant_genes_df <- combined_deg_df %>%
    mutate(
      genotype = case_when(
        grepl("Yorkie", comparison) ~ "Yorkie",
        grepl("Warts", comparison) ~ "Warts",
        grepl("Hippo", comparison) ~ "Hippo",
        TRUE ~ NA_character_
      ),
      direction = case_when(
        grepl("UP", comparison) ~ "UP",
        grepl("DOWN", comparison) ~ "DOWN",
        TRUE ~ NA_character_
      )
    ) %>%
    # Filter out rows where genotype or direction couldn't be determined
    filter(!is.na(genotype) & !is.na(direction))
  
  all_genes <- unique(significant_genes_df$gene_id)
  return(list(significant_genes_df = significant_genes_df, all_genes = all_genes))
}

## -----------------------------------------------------------------------------------------------------------------------------------------

#' Filters and saves specific DEG lists based on external filter tables.
#'
#' @param combined_df The comprehensive DEG dataframe.
#' @param output_dir Directory to save the output files.
#' @param filter_list A named list containing paths to 'strict' filter files.
#' @return A dataframe with DEGs after applying the strict filter.
save_filtered_deg_lists <- function(combined_df, output_dir, filter_list) {
  
  # Ensure the output directories exist
  homer_output_dir <- file.path(output_dir, "homer")
  deg_output_dir <- file.path(output_dir, "filtered_deg_lists")
  
  if (!dir.exists(homer_output_dir)) {
    dir.create(homer_output_dir, recursive = TRUE)
  }
  if (!dir.exists(deg_output_dir)) {
    dir.create(deg_output_dir, recursive = TRUE)
  }
  
  message("\n--- Loading Filter Tables ---")

  # Load strict filter genes (single column, no header)
  if (file.exists(filter_list$strict)) {
    strict_filter_genes_df <- read.table(filter_list$strict, header = FALSE, sep = "\t")
    strict_filter_genes <- unique(strict_filter_genes_df[[1]]) # Access the first column
    message(paste0("Loaded ", length(strict_filter_genes), " genes from strict filter: ", filter_list$strict))
  } else {
    strict_filter_genes <- c()
    warning(paste("Strict filter file not found:", filter_list$strict))
  }
  
  
  message("\n--- Applying Filters and Saving DEG lists ---")
  
  # Define the genotypes for iteration
  target_genotypes <- c("Warts", "Yorkie")
  
  for (genotype in target_genotypes) {
    # Filter for strict_filtered DEGs (must be absent from strict_filter_genes)
    strict_filtered_deg_df <- low_filtered_deg_df %>% # Start from low_filtered to apply strict on top
      filter(!(gene_id %in% strict_filter_genes))
    
    # Save strict_filtered DEGs for the current genotype
    if (nrow(strict_filtered_deg_df) > 0) {
      filename_strict_filtered <- file.path(deg_output_dir, paste0(genotype, "_strict_filtered_DEGs.txt"))
      write.table(strict_filtered_deg_df, filename_strict_filtered, row.names = FALSE, quote = FALSE, sep = "\t")
      message("Saved ", basename(filename_strict_filtered))
    } else {
      message("No strict-filtered DEGs found for ", genotype)
    }
    
    # Save gene_ids for Homer (strict_filtered)
    strict_filtered_genes_homer <- unique(strict_filtered_deg_df$gene_id)
    if (length(strict_filtered_genes_homer) > 0) {
      filename_strict_homer <- file.path(homer_output_dir, paste0(genotype, "_strict_filtered_genes_for_homer.txt"))
      write.table(strict_filtered_genes_homer, filename_strict_homer, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
      message("Saved ", basename(filename_strict_homer), " for Homer.")
    } else {
      message("No strict-filtered genes for Homer for ", genotype)
    }
  }
  
  return(strict_filtered_deg_df) # Return the strictly filtered dataframe
}

## -----------------------------------------------------------------------------------------------------------------------------------------
#' Filters deletion DEG tables by intersecting with insertion DEG tables,
#' ensuring genes are present and have the same direction (UP/DOWN),
#' and saves gene_id lists for Homer.
#'
#' @param deletion_deg_dir Directory containing the deletion DEG tables
#'        ("Warts_low_filtered_DEGs.txt", "Yorkie_low_filtered_DEGs.txt",
#'         "Warts_strict_filtered_DEGs.txt", "Yorkie_strict_filtered_DEGs.txt").
#' @param insertion_filter A named list of paths to the insertion DEG tables.
#' @param output_base_dir The base output directory where the 'filtered_deg_lists'
#'        and 'homer' subdirectories will be created for the intersected results.
#' @return None. Saves the filtered DEG tables and Homer gene lists to the output directory.
filter_by_insertion_and_direction <- function(deletion_deg_dir, insertion_filter, output_base_dir) {
  
  # Define specific output directories
  output_deg_dir <- file.path(output_base_dir, "filtered_deg_lists", "insertion_intersected_direction")
  output_homer_dir <- file.path(output_base_dir, "homer", "insertion_intersected_direction")
  
  # Create output directories if they don't exist
  if (!dir.exists(output_deg_dir)) {
    dir.create(output_deg_dir, recursive = TRUE)
  }
  if (!dir.exists(output_homer_dir)) {
    dir.create(output_homer_dir, recursive = TRUE)
  }
  
  # Helper function to load DEG data and select relevant columns (gene_id, direction)
  load_deg_data_for_intersection <- function(file_path) {
    if (file.exists(file_path)) {
      deg_table <- read.table(file_path, header = TRUE, sep = "\t")
      # Select gene_id and direction for intersection
      return(deg_table %>% dplyr::select(gene_id, direction))
    } else {
      warning(paste("File not found:", file_path))
      return(data.frame(gene_id = character(0), direction = character(0))) # Return empty dataframe
    }
  }
  
  message("\n--- Loading Insertion DEG Data for Intersection ---")
  # Load insertion DEG dataframes
  ins_strict_warts_df <- load_deg_data_for_intersection(insertion_filter$Ins_strict_WT_Warts)
  ins_strict_yorkie_df <- load_deg_data_for_intersection(insertion_filter$Ins_strict_WT_Yorkie)
  
  message("\n--- Filtering Deletion DEGs by Insertion Data and Direction ---")
  
  # Define a generic processing function to avoid repetition
  process_deg_pair <- function(deletion_file_path, insertion_df) {
    if (file.exists(deletion_file_path)) {
      deletion_deg_df <- read.table(deletion_file_path, header = TRUE, sep = "\t")
      
      # Perform an inner join on gene_id and direction
      # This ensures both presence and matching direction
      intersected_deg_df <- inner_join(deletion_deg_df, insertion_df, by = c("gene_id", "direction"))
      
      # --- Save full DEG table ---
      output_filename_base <- basename(deletion_file_path) %>% tools::file_path_sans_ext()
      output_full_deg_file <- file.path(output_deg_dir, paste0(output_filename_base, "_intersect_direction.txt"))
      write.table(intersected_deg_df, output_full_deg_file, row.names = FALSE, quote = FALSE, sep = "\t")
      message(paste("Saved", nrow(intersected_deg_df), "intersected (by gene_id and direction) DEGs to", output_full_deg_file))
      
      # --- Save gene_id list for Homer ---
      intersected_gene_ids <- unique(intersected_deg_df$gene_id)
      if (length(intersected_gene_ids) > 0) {
        output_homer_file <- file.path(output_homer_dir, paste0(output_filename_base, "_intersect_direction_for_homer.txt"))
        write.table(intersected_gene_ids, output_homer_file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
        message(paste("Saved", length(intersected_gene_ids), "gene IDs for Homer to", output_homer_file))
      } else {
        message(paste("No intersected gene IDs found for Homer from", deletion_file_path))
      }
      
    } else {
      warning(paste("Deletion DEG file not found:", deletion_file_path))
    }
  }

  # Process Warts strict filtered
  process_deg_pair(
    deletion_file_path = file.path(deletion_deg_dir, "Warts_strict_filtered_DEGs.txt"),
    insertion_df = ins_strict_warts_df
  )
  
  # Process Yorkie strict filtered
  process_deg_pair(
    deletion_file_path = file.path(deletion_deg_dir, "Yorkie_strict_filtered_DEGs.txt"),
    insertion_df = ins_strict_yorkie_df
  )
}


## -----------------------------------------------------------------------------------------------------------------------------------------
# --- Main Execution Loop ---

# 1. Load ALL relevant raw DEG data needed for the *entire* filtering process ONCE.
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
temp_prepared_all_data <- prepare_upset_data(all_raw_deg_data_combined)$significant_genes_df

# 3. Perform filtering based on the provided low and strict filter files.
filtered_deg_final <- save_filtered_deg_lists(
  combined_df = temp_prepared_all_data, # Pass the comprehensive data with genotype/direction
  output_dir = output_dir,
  filter_list = filter
)

# 4. Save DEG which are consistent between insertion and deletion KO
output_base_dir <- "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Heatmap"

# The directory where the *previously* filtered deletion DEG files are located
deletion_deg_directory <- file.path(output_base_dir, "filtered_deg_lists")

# Call the function
filter_by_insertion_and_direction(
  deletion_deg_dir = deletion_deg_directory,
  insertion_filter = insertion_filter,
  output_base_dir = output_base_dir
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



##-------------------- KO post strict filter vs WT -----------------------------------##

output_dir <- "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Heatmap/"
# You might want to create the directory if it doesn't exist:
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

conditions_Del_heatmap <- list(
  WT_Warts = "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Heatmap/filtered_deg_lists/Warts_strict_filtered_DEGs.txt",
  WT_Yorkie = "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Heatmap/filtered_deg_lists/Yorkie_strict_filtered_DEGs.txt"
)

dds <- readRDS("W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Del_out_deseq2_dds.rds")

# Prepare data for deletion heatmap
heatmap_data_del <- prepare_heatmap_data(
  conditions_list = conditions_Del_heatmap,
  num_top_genes = 100, # Original script selected 30
  padj_thresh = 0.05,
  log2fc_thresh = 0,
  columns_to_select = c(1:4) # Columns relevant for deletion data
)

list_top_genes_del <- heatmap_data_del$list_top_genes
list_top_genes_del <- gsub("PT","gene:PT",list_top_genes_del)
combined_top_genes_del <- heatmap_data_del$combined_top_genes_df
combined_top_genes_del <- gsub("PT","gene:PT",combined_top_genes_del)

# Save all unique DEG (top 60) for deletion heatmap
write.table(sub("gene:", "", list_top_genes_del),
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
vsd_del <- vst(dds)
transformed_counts_del <- assay(vsd_del)
mat_del <- transformed_counts_del[list_top_genes_del, ]

# Prepare annotation for Deletion heatmap 
anno_del <- as.data.frame(colData(dds)[, "genotype", drop = FALSE])
cond_del <- c("Sr", "YC","WB") 
anno_del <- anno_del[anno_del$genotype %in% cond_del, , drop = FALSE]

# Modify genotype labels (without "algo" specific labels)
anno_del$genotype <- gsub("Sr", "WT", anno_del$genotype)
anno_del$genotype <- gsub("WB", "Warts KO", anno_del$genotype)
anno_del$genotype <- gsub("YC", "Yorkie KO", anno_del$genotype)
anno_del$genotype <- as.factor(anno_del$genotype)
anno_del$genotype <- droplevels(anno_del$genotype)

# Define custom colors for deletion genotypes
genotype_colors_del <- c(
  "WT" = "darkblue",
  "Yorkie KO" = "darkorange",
  "Warts KO" = "darkred"
)

# Plot deletion heatmap
plot_expression_heatmap(
  mat = mat_del,
  anno = anno_del,
  genotype_colors = genotype_colors_del,
  output_file = "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Heatmap/Heatmap_KO_WT_filtered_DEG.svg",
  row_km = 4,
  column_km = 3,
  row_font_size = 10,
  column_font_size = 10,
  width = 16, 
  height = 30
)

##-------------------- KO post strict filter vs WT + Consistent DEG -----------------------------------##

output_dir <- "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Heatmap/"
# You might want to create the directory if it doesn't exist:
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

conditions_Del_heatmap <- list(
  WT_Warts = "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Heatmap/filtered_deg_lists/insertion_intersected_direction/Warts_strict_filtered_DEGs_intersect_direction.txt"
)

dds <- readRDS("W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Del_out_deseq2_dds.rds")

# Prepare data for deletion heatmap
heatmap_data_del <- prepare_heatmap_data(
  conditions_list = conditions_Del_heatmap,
  num_top_genes = 90, # Original script selected 30
  padj_thresh = 0.05,
  log2fc_thresh = 0,
  columns_to_select = c(1:4) # Columns relevant for deletion data
)

# Decide of gene order for the Heatmap
genes_to_keep <- c(
  "PTSG_03144","PTSG_04943", "PTSG_09711","PTSG_10931","PTSG_11789","PTSG_12667",
  "PTSG_01648","PTSG_03426","PTSG_05995","PTSG_07368",
  "PTSG_01647","PTSG_04961","PTSG_09125","PTSG_09127",
  "PTSG_12062","PTSG_00031","PTSG_07323","PTSG_01855","PTSG_13168","PTSG_05312",
  "PTSG_09328","PTSG_07640"
)

# Filter rows where any column matches a gene in genes_to_keep
list_top_genes_del <- heatmap_data_del$list_top_genes
list_top_genes_del <- genes_to_keep[genes_to_keep %in% list_top_genes_del]
list_top_genes_del <- gsub("PT","gene:PT",list_top_genes_del)

combined_top_genes_del <- heatmap_data_del$combined_top_genes_df%>%
  filter(WT_Warts %in% genes_to_keep )
combined_top_genes_del <- gsub("PT","gene:PT",combined_top_genes_del)

# Save all unique DEG (top 60) for deletion heatmap
write.table(sub("gene:", "", list_top_genes_del),
            file = file.path(output_dir, "KO_WT_consistent_filtered-DEG.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)

# Run Uniprot API function to get uniprot long name with the gene present on the heatmap
input_ids_for_uniprot <- read.table(
  file.path(output_dir, "KO_WT_consistent_filtered-DEG.txt"),
  sep = "\t", header = FALSE
)$V1

organism_ncbi_id <- "946362" # Salpingoeca rosetta
organism_abbr_name <- "Sros"
input_is_gene_name <- FALSE # Set to TRUE if input_ids are gene symbols, FALSE if they are protein/gene IDs

output_conversion_prefix <- file.path(output_dir, "KO_WT_consistent_filtered-DEG_uniprot")

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
vsd_del <- vst(dds)
transformed_counts_del <- assay(vsd_del)
mat_del <- transformed_counts_del[list_top_genes_del, ]

# Prepare annotation for Deletion heatmap 
anno_del <- as.data.frame(colData(dds)[, "genotype", drop = FALSE])
cond_del <- c("Sr", "YC","WB") 
anno_del <- anno_del[anno_del$genotype %in% cond_del, , drop = FALSE]

# Modify genotype labels (without "algo" specific labels)
anno_del$genotype <- gsub("Sr", "WT", anno_del$genotype)
anno_del$genotype <- gsub("WB", "Warts KO", anno_del$genotype)
anno_del$genotype <- gsub("YC", "Yorkie KO", anno_del$genotype)
anno_del$genotype <- as.factor(anno_del$genotype)
anno_del$genotype <- droplevels(anno_del$genotype)

# Define custom colors for deletion genotypes
genotype_colors_del <- c(
  "WT" = "darkblue",
  "Yorkie KO" = "darkorange",
  "Warts KO" = "darkred"
)

# Plot deletion heatmap
plot_expression_heatmap(
  mat = mat_del,
  anno = anno_del,
  genotype_colors = genotype_colors_del,
  output_file = "W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Heatmap/Heatmap_KO_WT_consistent_filtered_DEG.svg",
  row_km = 0,
  column_km = 0,
  row_font_size = 10,
  column_font_size = 10,
  width = 10, 
  height = 8
)

