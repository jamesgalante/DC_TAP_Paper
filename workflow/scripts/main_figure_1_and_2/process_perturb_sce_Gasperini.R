# Script: process_perturb_sce_Gasperini.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/process_perturb_sce_Gasperini.rda"))
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
  library(SingleCellExperiment)
})

message("Loading input files")
perturb_sce <- readRDS(snakemake@input$perturb_sce)


### EXTRACT DATA FROM PERTURB_SCE =============================================

# Extract the raw counts matrix
raw_counts <- assay(perturb_sce, "counts")

# Extract the binarized guide counts matrix
binarized_guide_counts <- assay(altExp(perturb_sce, "grna_perts"), "perts")

# Extract guide information
guide_targets <- as.data.frame(rowData(altExp(perturb_sce, "grna_perts"))) %>% 
  filter(target_type == "enh")
rownames(guide_targets) <- NULL


### FILTER BASED ON ZERO COUNTS ===============================================

# Let's filter the raw_counts matrix, so the column names are matched
raw_counts <- raw_counts[, colSums(binarized_guide_counts) != 0]

# There is an issue where grna_n_nonzero is 0, so the log(grna_n_nonzero) == -Inf ; To fix this, we must remove these cells
binarized_guide_counts <- binarized_guide_counts[, colSums(binarized_guide_counts) != 0]


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
saveRDS(raw_counts, snakemake@output$raw_counts)
saveRDS(binarized_guide_counts, snakemake@output$binarized_guide_counts)
write_tsv(guide_targets, snakemake@output$guide_targets)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)