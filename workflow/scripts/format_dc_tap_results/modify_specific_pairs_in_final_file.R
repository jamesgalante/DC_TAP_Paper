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
element_gene_pairs <- read_tsv(snakemake@input$results_with_element_gene_pair_categories)

# Read in the SCEPTRE discovery results with confidence intervals
k562_sceptre_results_with_CIs <- readRDS(snakemake@input$discovery_results_w_CIs[[1]]) %>% mutate(cell_type = "K562")
wtc11_sceptre_results_with_CIs <- readRDS(snakemake@input$discovery_results_w_CIs[[2]]) %>% mutate(cell_type = "WTC11")
sceptre_results_with_CIs <- rbind(k562_sceptre_results_with_CIs, wtc11_sceptre_results_with_CIs)


### MODIFY SPECIFIC PAIRS =====================================================

message("Modifying specific element-gene pairs")
modified_pairs <- element_gene_pairs %>%
  mutate(
    # RTN4_TSS_8: Change from promoter to enhancer
    element_location = ifelse(design_file_target_name == "RTN4_TSS_8", "distal", element_location),
    
    # CCDC26: Fix overlapping gene isoform that's not expressed
    DistalElement_Gene = ifelse(element_gene_pair_identifier_hg19 == "CCDC26|chr8:130594299-130594600", 
                                TRUE, DistalElement_Gene),
    Positive_Control_DistalElement_Gene = ifelse(element_gene_pair_identifier_hg19 == "CCDC26|chr8:130594299-130594600", 
                                                 TRUE, Positive_Control_DistalElement_Gene)
  )


### CALCULATE SUMMARY STATISTICS ==============================================

message("Calculating summary statistics by cell type")

# Helper functions for summarizing pair stats
summarize_pairs <- function(df, category, cell_type) {
  # Get significant pair stats (upregulated vs downregulated)
  sig_stats <- df %>%
    filter(!!sym(category), significant == TRUE, cell_type == !!cell_type) %>%
    summarize(
      Upregulated = sum(pct_change_effect_size >= 0),
      Downregulated = sum(pct_change_effect_size < 0)
    )
  
  # Get non-significant pair stats (well-powered vs underpowered)
  nonsig_stats <- df %>%
    filter(!!sym(category), significant == FALSE, cell_type == !!cell_type) %>%
    summarize(
      WellPowered = sum(power_at_effect_size_15 >= 0.8),
      UnderPowered = sum(power_at_effect_size_15 < 0.8)
    )
  
  # Combine the stats
  tibble(
    Category = category,
    Upregulated = sig_stats$Upregulated,
    Downregulated = sig_stats$Downregulated,
    WellPowered = nonsig_stats$WellPowered,
    UnderPowered = nonsig_stats$UnderPowered
  )
}

# Categories to summarize
categories <- c(
  "DistalElement_Gene",
  "selfPromoter",
  "DistalPromoter_Gene",
  "Positive_Control_DistalElement_Gene",
  "Random_DistalElement_Gene"
)

# Generate summary tables for each cell type
summary_K562 <- map_dfr(categories, ~summarize_pairs(modified_pairs, .x, "K562"))
summary_WTC11 <- map_dfr(categories, ~summarize_pairs(modified_pairs, .x, "WTC11"))

# Display summary tables for verification
message("Summary statistics for K562:")
print(summary_K562)
message("\nSummary statistics for WTC11:")
print(summary_WTC11)


### REMOVE UNDERPOWERED PAIRS =================================================

message("Removing underpowered non-significant pairs")
final_pairs <- modified_pairs %>%
  mutate(across(
    all_of(categories),
    ~if_else(significant == FALSE & power_at_effect_size_15 < 0.8, FALSE, .x)
  ))

### ADD IN CONFIDENCE INTERVALS ===============================================

# Add in the k562 and wtc11 results by gene_id + design_file_target_name + cell_type
final_pairs <- final_pairs %>%
  left_join(
    sceptre_results_with_CIs %>% 
      select(response_id, grna_target, se_fold_change, cell_type), 
    by = c("gene_id" = "response_id", "design_file_target_name" = "grna_target", "cell_type")
  ) %>%
  mutate(
    standard_error_pct_change_effect_size = se_fold_change * 100
  ) %>%
  select(
    -se_fold_change
  ) %>%
  relocate(
    standard_error_pct_change_effect_size, .after = pct_change_effect_size
  )


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(final_pairs, snakemake@output$results_with_element_gene_pair_categories_modified)
write_tsv(summary_K562, snakemake@output$summary_K562)
write_tsv(summary_WTC11, snakemake@output$summary_WTC11)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)