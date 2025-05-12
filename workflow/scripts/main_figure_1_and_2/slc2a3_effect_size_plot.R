# Script: slc2a3_effect_size_plot.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/slc2a3_effect_size_plot.rda"))
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
results_with_element_gene_pair_categories_modified <- read_tsv(snakemake@input$results_with_element_gene_pair_categories_modified)
wtc11_calibration_check_results <- readRDS(snakemake@input$wtc11_calibration_check_results)


### CALCULATE CIs FOR NEG CONTROLS ============================================

# Filter for the pair of interest
results_with_element_gene_pair_categories_modified %>% 
  filter(cell_type == "WTC11", element_gene_pair_identifier_hg38 == "SLC2A3|chr12:7960676-7960977")


# For calibration check - 15 negative control guides are tested against the gene of interest: "ENSG00000059804"
# Thus a confidence interval isn't calculated on all negative control pairs at the same time - I have each neg ctrl guide tested against "ENSG00000059804"
# But I don't have a standard error.
# I think will have to manually calculate mu and y from the final sceptre object of the singleton analysis (is this different...)
# Have to email Eugeen


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)