# Script: compare_replicates.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/compare_replicates.rda"))
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



###  Doing Thing =========================================================

# See maddie code about sampling procedure for splits - because did a lot of splits

# take original k562 outputs - split by lane (suffix of barcode)
# run sceptre on each split
# Compare effect sizes between splits


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)