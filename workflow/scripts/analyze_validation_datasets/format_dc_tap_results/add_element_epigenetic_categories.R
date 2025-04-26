# Script: add_element_epigenetic_categories.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/add_element_epigenetic_categories.rda"))
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
resized_Formatted_DC_TAP_Seq_Results <- read_tsv(snakemake@input$resized_Formatted_DC_TAP_Seq_Results)
categorized_data <- read_tsv(snakemake@input$categorized_data)


### ADDING CATEGORIES =========================================================

# Add the categories to the resized_Formatted data
Formatted_DC_TAP_Seq_Results_w_Categories <- resized_Formatted_DC_TAP_Seq_Results %>%
  left_join(
    categorized_data %>%
      select(element_gene_pair_identifier_hg38, element_category, element_category_simple, element_category_with_dnase, CellType) %>%
      dplyr::rename(
        cell_type = CellType
      ),
    by = c("element_gene_pair_identifier_hg38", "cell_type")
  )


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(Formatted_DC_TAP_Seq_Results_w_Categories, snakemake@output$Formatted_DC_TAP_Seq_Results_w_Categories)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)