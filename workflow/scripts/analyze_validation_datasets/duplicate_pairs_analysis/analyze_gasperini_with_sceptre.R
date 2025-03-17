# Script: analyze_gasperini_with_sceptre.R

### SETUP =====================================================================

# Saving image for debugging
save.image("RDA_objects/analyze_gasperini_with_sceptre.rda")
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

# Load in the necessary packages
message("Loading packages")
suppressPackageStartupMessages({
  library(devtools)
  library(tidyverse)
  library(Matrix)
})

# Download Sceptre
message("Installing Sceptre")
devtools::install_github("katsevich-lab/sceptre")
message("Sceptre Installation Complete")
message("Loading Sceptre")
library(sceptre)
message("Sceptre Loading Complete")

# Load input files
message("Loading input files")
raw_counts <- readRDS(snakemake@input$raw_counts)
binarized_guide_counts <- readRDS(snakemake@input$binarized_guide_counts)
response_id_target_pairs <- read_tsv(snakemake@input$response_id_target_pairs)
guide_target_pairs <- read_tsv(snakemake@input$guide_target_pairs)


### CREATE SCEPTRE OBJECT =====================================================

# Because the pre-calculated grna_n_nonzero and grna_n_umis covariates are not independent (due to a binary grna_matrix) we must remove one from the formula
new_formula <- formula(~ log(response_n_nonzero) + log(response_n_umis) + log(grna_n_nonzero))

# There is also an issue where grna_n_nonzero is 0, so the log(grna_n_nonzero) == -Inf ; To fix this, we must remove these cells
filtered_grna_matrix <- binarized_guide_counts[, colSums(binarized_guide_counts) != 0]
# Let's filter the raw_counts matrix as well, so the column names are matched
filtered_response_matrix <- raw_counts[, colSums(binarized_guide_counts) != 0]

# Initialize sceptre object
message("Initializing Sceptre Object")
sceptre_object <- import_data(
  response_matrix = filtered_response_matrix,
  grna_matrix = filtered_grna_matrix,
  grna_target_data_frame = guide_target_pairs,
  moi = "high",
  response_names = rownames(raw_counts)
)

# Set analysis parameters
message("Setting analysis parameters")
sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = response_id_target_pairs,
  side = "both",
  grna_integration_strategy = "union",
  control_group = "complement",
  multiple_testing_method = "BH",
  formula_object = new_formula
)


### RUN SCEPTRE ===============================================================

# Assign grnas based on a threshold of 1 (as the gasperini data has already been assigned guides)
# For the sake of this analysis, this is okay, but if the results are weird, we may want to run mixture assignment
sceptre_object <- assign_grnas(
  sceptre_object = sceptre_object,
  method = "thresholding",
  threshold = 1
)

# Run quality control with default settings
sceptre_object <- run_qc(sceptre_object = sceptre_object)

# Calibration check can be skipped as there are no negative controls
# Power check can be skipped as there are no negative or positive controls (i filtered out pos controls in prior script)

# Run discovery analysis
sceptre_object <- run_discovery_analysis(sceptre_object = sceptre_object)


### SAVE OUTPUT ===============================================================

# Write all sceptre outputs to an output directory
write_outputs_to_directory(
  sceptre_object = sceptre_object, 
  directory = dirname(snakemake@output$discovery_results)
)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

