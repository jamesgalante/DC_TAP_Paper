# Script: resize_and_merge_dc_tap_elements.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/resize_and_merge_dc_tap_elements.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
})

message("Loading input files")
results_with_element_gene_pair_categories_modified <- read_tsv(snakemake@input$results_with_element_gene_pair_categories_modified)


### RESIZE AND MERGE FUNCTIONS ================================================

message("Defining functions")
# Function to resize and merge elements in crispr data
resize_crispr_elements <- function(crispr, size = 500) {
  
  # Extract unique elements and create GRanges object
  elements <- distinct(select(crispr, targeting_chr_hg38, targeting_start_hg38, targeting_end_hg38, intended_target_name_hg38))
  elements <- makeGRangesFromDataFrame(elements, seqnames.field = "targeting_chr_hg38",
                                       start.field = "targeting_start_hg38", end.field = "targeting_end_hg38",
                                       keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
  
  # Resize elements and merge overlaps
  resize_elements <- width(elements) < size
  elements[resize_elements] <- resize(elements[resize_elements], width = size, fix = "center")
  
  # Get merged element uids for each resized element
  elements_merged <- reduce(elements, with.revmap = TRUE)
  merged_uids <- lapply(elements_merged$revmap, FUN = function(x) elements[x]$intended_target_name_hg38)
  elements_merged$merged_uid <- merged_uids
  
  # Convert to table in long format
  elements_merged_df <- elements_merged %>% 
    as.data.frame() %>% 
    select(targeting_chr_hg38 = seqnames, targeting_start_hg38 = start, targeting_end_hg38 = end, merged_uid) %>% 
    unnest(cols = merged_uid)
  
  # Add resized element coordinates to crispr results based on original merged element uids
  crispr_merged <- crispr %>% 
    left_join(
      elements_merged_df %>% 
        dplyr::rename( 
          resized_merged_targeting_chr_hg38 = targeting_chr_hg38, 
          resized_merged_targeting_start_hg38 = targeting_start_hg38, 
          resized_merged_targeting_end_hg38 = targeting_end_hg38 
        ), 
      by = c("intended_target_name_hg38" = "merged_uid")
    )
  
  # For every gene - merged element pair, summarize the Regulated column
  crispr_merged <- crispr_merged %>% 
    mutate(resized_merged_element_gene_pair_identifier_hg38 = paste0(gene_symbol, "|", resized_merged_targeting_chr_hg38, ":", resized_merged_targeting_start_hg38, "-", resized_merged_targeting_end_hg38)) %>%
  
  return(crispr_merged)
}


### DATA PROCESSING ===========================================================

message("Resizing and merging elements")
# Merge elements for DC-TAP-seq datasets
k562_dc_tap <- resize_crispr_elements(filter(results_with_element_gene_pair_categories_modified, cell_type == "K562"), size = 500)
wtc11_dc_tap <- resize_crispr_elements(filter(results_with_element_gene_pair_categories_modified, cell_type == "WTC11"), size = 500)

# Combine with other CRISPR data to create new combined file
formatted_resized_and_merged_dc_tap_output <- bind_rows(k562_dc_tap, wtc11_dc_tap)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(formatted_resized_and_merged_dc_tap_output, snakemake@output$resized_and_merged_input_for_chromatin_categorization_pipeline)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)



# # Find merged elements in the K562 data
# k562_merged_elements <- k562_dc_tap %>%
#   # Group by merged coordinates
#   group_by(resized_merged_targeting_chr_hg38, 
#            resized_merged_targeting_start_hg38, 
#            resized_merged_targeting_end_hg38) %>%
#   # Count unique original elements for each merged location
#   summarize(
#     num_original_elements = n_distinct(intended_target_name_hg38),
#     original_elements = paste(unique(intended_target_name_hg38), collapse = ", "),
#     .groups = "drop"
#   ) %>%
#   # Keep only cases with multiple original elements merged
#   filter(num_original_elements > 1) %>%
#   # Sort by chromosome and position
#   arrange(resized_merged_targeting_chr_hg38, resized_merged_targeting_start_hg38)
# 
# # Get detailed information about the merged elements
# k562_merged_details <- k562_dc_tap %>%
#   # Find rows belonging to merged elements
#   inner_join(
#     k562_merged_elements %>% 
#       select(resized_merged_targeting_chr_hg38, 
#              resized_merged_targeting_start_hg38, 
#              resized_merged_targeting_end_hg38),
#     by = c("resized_merged_targeting_chr_hg38", 
#            "resized_merged_targeting_start_hg38", 
#            "resized_merged_targeting_end_hg38")
#   ) %>%
#   # Select all columns including the annotation columns
#   select(
#     original_element = intended_target_name_hg38,
#     original_chr = targeting_chr_hg38,
#     original_start = targeting_start_hg38,
#     original_end = targeting_end_hg38,
#     merged_chr = resized_merged_targeting_chr_hg38,
#     merged_start = resized_merged_targeting_start_hg38,
#     merged_end = resized_merged_targeting_end_hg38,
#     gene_symbol,
#     significant,
#     DistalElement_Gene,
#     Random_DistalElement_Gene,
#     Positive_Control_DistalElement_Gene,
#     selfPromoter,
#     DistalPromoter_Gene
#   ) %>%
#   # Sort by merged coordinates and original element
#   arrange(merged_chr, merged_start, original_element)
# 
# # Get unique pairs while preserving all columns
# distinct_pairs <- k562_merged_details %>%
#   distinct(original_element, merged_chr, merged_start, merged_end, .keep_all = TRUE)
# 
# See all merged pairs where the individual elements for a given gene were differentially significant
# blah <- k562_merged_details %>%
#   # Create merged_element identifier
#   mutate(merged_element = paste0(gene_symbol, "|", merged_chr, ":", merged_start, "-", merged_end)) %>%
#   # Group by merged_element
#   group_by(merged_element) %>%
#   # Filter to keep only merged elements with different significant values
#   filter(n_distinct(significant) > 1) %>%
#   # Ungroup for further operations
#   ungroup()