# Script: summary_of_element_gene_categories.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/summary_of_element_gene_categories.rda"))
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
results_wo_pos_controls <- read_tsv(snakemake@input$results_wo_pos_controls)


### CALCULATE SUMMARY STATISTICS ==============================================

message("Calculating summary statistics by cell type")

# Helper functions for summarizing pair stats (updated to include wo_pos_controls)
# NOTE: Still only using effect size 15 for summary statistics as requested
summarize_pairs <- function(df, category, cell_type) {
  # Get significant pair stats (upregulated vs downregulated) - original data
  sig_stats <- df %>%
    filter(!!sym(category), significant == TRUE, cell_type == !!cell_type) %>%
    summarize(
      Upregulated = sum(pct_change_effect_size >= 0, na.rm = TRUE),
      Downregulated = sum(pct_change_effect_size < 0, na.rm = TRUE)
    )
  
  # Get non-significant pair stats (well-powered vs underpowered) - original data
  nonsig_stats <- df %>%
    filter(!!sym(category), significant == FALSE, cell_type == !!cell_type) %>%
    summarize(
      WellPowered = sum(power_at_effect_size_15 >= 0.8, na.rm = TRUE),
      UnderPowered = sum(power_at_effect_size_15 < 0.8, na.rm = TRUE)
    )
  
  # Get non-significant pair stats for wo_pos_controls data
  nonsig_stats_wo_pos <- df %>%
    filter(!!sym(category), significant_wo_pos_controls == FALSE, cell_type == !!cell_type) %>%
    summarize(
      WellPowered_wo_pos_controls = sum(power_at_effect_size_15_wo_pos_controls >= 0.8, na.rm = TRUE),
      UnderPowered_wo_pos_controls = sum(power_at_effect_size_15_wo_pos_controls < 0.8, na.rm = TRUE)
    )
  
  # Get significant pair stats for wo_pos_controls data (custom FDR)
  sig_stats_wo_pos <- df %>%
    filter(!!sym(category), significant_wo_pos_controls == TRUE, cell_type == !!cell_type) %>%
    summarize(
      Upregulated_wo_pos_controls = sum(pct_change_effect_size >= 0, na.rm = TRUE),
      Downregulated_wo_pos_controls = sum(pct_change_effect_size < 0, na.rm = TRUE)
    )
  
  # Combine the stats
  result <- tibble(
    Category = category,
    Upregulated = sig_stats$Upregulated,
    Downregulated = sig_stats$Downregulated,
    WellPowered = nonsig_stats$WellPowered,
    UnderPowered = nonsig_stats$UnderPowered,
    Upregulated_wo_pos_controls = sig_stats_wo_pos$Upregulated_wo_pos_controls,
    Downregulated_wo_pos_controls = sig_stats_wo_pos$Downregulated_wo_pos_controls,
    WellPowered_wo_pos_controls = nonsig_stats_wo_pos$WellPowered_wo_pos_controls,
    UnderPowered_wo_pos_controls = nonsig_stats_wo_pos$UnderPowered_wo_pos_controls
  )
  
  # Set wo_pos_controls columns to NA for excluded categories
  if(category %in% c("selfPromoter", "Positive_Control_DistalElement_Gene", "Positive_Control_selfPromoter")) {
    result$Upregulated_wo_pos_controls <- NA_integer_
    result$Downregulated_wo_pos_controls <- NA_integer_
    result$WellPowered_wo_pos_controls <- NA_integer_
    result$UnderPowered_wo_pos_controls <- NA_integer_
  }
  
  return(result)
}

# Categories to summarize
categories <- c(
  "DistalElement_Gene",
  "selfPromoter",
  "DistalPromoter_Gene",
  "Positive_Control_DistalElement_Gene",
  "Positive_Control_selfPromoter", 
  "Random_DistalElement_Gene"
)

# Generate summary tables for each cell type
summary_K562 <- map_dfr(categories, ~summarize_pairs(results_wo_pos_controls, .x, "K562")) %>% mutate(cell_type = "K562")
summary_WTC11 <- map_dfr(categories, ~summarize_pairs(results_wo_pos_controls, .x, "WTC11")) %>% mutate(cell_type = "WTC11")


### COMBINING THE TABLES ======================================================

full_summary <- rbind(summary_K562, summary_WTC11)
colnames(full_summary) <- c("Pair_Category", "Significant_Upregulated", "Significant_Downregulated", "Non_Significant_Well_Powered", "Non_Significant_Under_Powered", "Significant_Upregulated_no_Pos_Control", "Significant_Downregulated_no_Pos_Control", "Non_Significant_Well_Powered_no_Pos_Control", "Non_Significant_Under_Powered_no_Pos_Control", "Cell_Type")


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(full_summary, snakemake@output$summary_of_element_gene_categories_supplementary_table, col_names = T)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)