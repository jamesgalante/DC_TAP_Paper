# Script: modify_specific_pairs_in_final_file.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/modify_specific_pairs_in_final_file.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
suppressPackageStartupMessages(
  library(tidyverse)
)

message("Loading input files")
summarized_categories <- read_tsv(snakemake@input$summarized_categories)


### MODIFYING =================================================================

# Call design_file_target_name = RTN4_TSS_8 an enhancer as it's currently labelled as a promoter, but is not an active TSS in IGV
# Change element_location to "distal"
# Change DistalElement_Gene to TRUE
  # Keep Random_DistalElement_Gene FALSE as this element was not chosen randomly and was initially intended to be a positive control

# "CCDC26|chr8:130594299-130594600" can be converted to Positive Control DistalElement Gene
# element_gene_pair_identifier_hg38 == "CCDC26|chr8:130594299-130594600"
# DistalElement_Gene = TRUE
# Positive_Control_DistalElement_Gene = TRUE

Formatted_DC_TAP_Seq_Results <- summarized_categories %>%
  mutate(
    # Modify RTN4_TSS_8 to be an enhancer
    element_location = ifelse(design_file_target_name == "RTN4_TSS_8", "distal", element_location),
    
    # Set DistalElement_Gene to TRUE for RTN4_TSS_8
    DistalElement_Gene = ifelse(design_file_target_name == "RTN4_TSS_8", TRUE, DistalElement_Gene),
    
    # Modify CCDC26 entry (overlaps target gene intron - but of an isoform that's not expressed)
    DistalElement_Gene = ifelse(element_gene_pair_identifier_hg19 == "CCDC26|chr8:130594299-130594600", TRUE, DistalElement_Gene)
  )

### CALCULATE SUMMARY STATISTICS ==============================================

# For each category, calculate the number of Significant pairs (upregulated v. downregulated) - calculate the the number of Non-significant pairs (high power v. underpowered)
get_significant_pair_stats <- function(df, category) {
  df %>%
    filter(!!sym(category), significant == T) %>%
    mutate(Downregulated = pct_change_effect_size < 0) %>%
    select(significant, Downregulated) %>%
    table()
}
get_nonsignificant_pair_stats <- function(df, category) {
  df %>%
    filter(!!sym(category), significant == F) %>%
    mutate(Underpowered = ifelse(power_at_effect_size_15 < 0.8, TRUE, FALSE)) %>%
    select(significant, Underpowered) %>%
    table()
}

# List of categories to summarize
categories <- c(
  "DistalElement_Gene",
  "selfPromoter",
  "DistalPromoter_Gene",
  "Positive_Control_DistalElement_Gene",
  "Random_DistalElement_Gene"
)

# Function to build summary table for a given cell type
get_summary_table_by_celltype <- function(filter_cell_type) {
  df_cell <- Formatted_DC_TAP_Seq_Results %>% filter(cell_type == filter_cell_type)
  
  map_dfr(categories, function(cat) {
    sig_tbl <- get_significant_pair_stats(df_cell, cat)
    nonsig_tbl <- get_nonsignificant_pair_stats(df_cell, cat)
    
    # Safely retrieve table values or 0 if not present
    downreg_false <- if ("TRUE" %in% rownames(sig_tbl) && "FALSE" %in% colnames(sig_tbl)) {
      sig_tbl["TRUE", "FALSE"]
    } else 0
    
    downreg_true <- if ("TRUE" %in% rownames(sig_tbl) && "TRUE" %in% colnames(sig_tbl)) {
      sig_tbl["TRUE", "TRUE"]
    } else 0
    
    underpowered_false <- if ("FALSE" %in% rownames(nonsig_tbl) && "FALSE" %in% colnames(nonsig_tbl)) {
      nonsig_tbl["FALSE", "FALSE"]
    } else 0
    
    underpowered_true <- if ("FALSE" %in% rownames(nonsig_tbl) && "TRUE" %in% colnames(nonsig_tbl)) {
      nonsig_tbl["FALSE", "TRUE"]
    } else 0
    
    tibble(
      Category = cat,
      Upregulated = downreg_false,
      Downregulated = downreg_true,
      WellPowered = underpowered_false,
      UnderPowered = underpowered_true
    )
  })
}

# Generate summary tables for each cell type
summary_K562 <- get_summary_table_by_celltype("K562")
summary_WTC11 <- get_summary_table_by_celltype("WTC11")

# Print the summary tables
cat("Summary table for K562:\n")
print(summary_K562)

cat("\nSummary table for WTC11:\n")
print(summary_WTC11)


### REMOVE UNDERPOWERED PAIRS =================================================

# We wanted to get the statistics on how many pairs are underpowered, but we don't want this in the final table
# Modify the categories to remove underpowered nonsignificant pairs
Formatted_DC_TAP_Seq_Results <- Formatted_DC_TAP_Seq_Results %>%
  mutate(
    across(
      c(DistalElement_Gene, selfPromoter, DistalPromoter_Gene, Positive_Control_DistalElement_Gene, Random_DistalElement_Gene),
      ~ case_when(
        significant == FALSE & power_at_effect_size_15 < 0.8 ~ FALSE,
        TRUE ~ .x
      )
    )
  )


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(Formatted_DC_TAP_Seq_Results, snakemake@output$Formatted_DC_TAP_Seq_Results)
write_tsv(summary_K562, snakemake@output$summary_K562)
write_tsv(summary_WTC11, snakemake@output$summary_WTC11)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)