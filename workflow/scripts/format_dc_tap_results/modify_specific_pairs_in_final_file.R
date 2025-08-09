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


### ADD IN CONFIDENCE INTERVALS ===============================================

# Add in the k562 and wtc11 results by gene_id + design_file_target_name + cell_type
final_pairs <- modified_pairs %>%
  left_join(
    sceptre_results_with_CIs %>% 
      select(response_id, grna_target, fold_change_effect_size = fold_change, se_fold_change, cell_type), 
    by = c("gene_id" = "response_id", "design_file_target_name" = "grna_target", "cell_type")
  ) %>%
  # Calculate all effect size statistics
  mutate(
    # Standard errors
    standard_error_fold_change = se_fold_change,
    standard_error_pct_change = se_fold_change * 100,
    standard_error_log_2_FC = se_fold_change / (fold_change_effect_size * log(2)),
    
    # Confidence intervals for fold change
    lower_CI_95_fold_change = fold_change_effect_size - 1.96 * standard_error_fold_change,
    upper_CI_95_fold_change = fold_change_effect_size + 1.96 * standard_error_fold_change,
    
    # Confidence intervals for log2 fold change
    lower_CI_95_log_2_FC = log_2_FC_effect_size - 1.96 * standard_error_log_2_FC,
    upper_CI_95_log_2_FC = log_2_FC_effect_size + 1.96 * standard_error_log_2_FC,
    
    # Confidence intervals for percent change
    lower_CI_95_pct_change = pct_change_effect_size - 1.96 * standard_error_pct_change,
    upper_CI_95_pct_change = pct_change_effect_size + 1.96 * standard_error_pct_change
  ) %>%
  # Remove original SE column after calculations
  select(-se_fold_change) %>%
  # Rearrange columns in the specified order
  relocate(
    fold_change_effect_size, log_2_FC_effect_size, pct_change_effect_size,
    standard_error_fold_change, standard_error_log_2_FC, standard_error_pct_change,
    lower_CI_95_fold_change, upper_CI_95_fold_change,
    lower_CI_95_log_2_FC, upper_CI_95_log_2_FC,
    lower_CI_95_pct_change, upper_CI_95_pct_change, .after = cell_type
  )


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(final_pairs, snakemake@output$results_with_element_gene_pair_categories_modified)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)