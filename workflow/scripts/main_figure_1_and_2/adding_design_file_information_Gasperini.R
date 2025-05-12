# Script: adding_design_file_information_Gasperini.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/adding_design_file_information_Gasperini.rda"))
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
  library(cowplot)
  library(GenomicRanges)
})

message("Loading input files")
gasperini <- read_tsv(snakemake@input$gasperini)

# Load in the guide_targets files
guide_targets <- read_tsv(snakemake@input$guide_targets)


### ADD TSS CONTROL INFORMATION ===============================================

# Combine the guide_targets and the gasperini by their coordinates (hg19) for each target
# Combine the guide targets
guide_targets_all <- guide_targets %>%
  group_by(target_name) %>%
  dplyr::slice(1)

# Parse gasperini: extract coordinate information
combined_parsed <- gasperini %>%
  # Split 'name' at the "|" and discard the gene part (we already have measuredGeneSymbol)
  separate(name, into = c(NA, "hg19_target_coords"), sep = "\\|", remove = FALSE) %>%
  # Remove trailing ":." from the coordinates and extract components
  mutate(
    hg19_target_coords = str_remove(hg19_target_coords, ":\\.$"),
    hg19_target_chr   = str_extract(hg19_target_coords, "chr[^:]+"),
    hg19_target_start = as.numeric(str_extract(hg19_target_coords, "(?<=:)[0-9]+")),
    hg19_target_end   = as.numeric(str_extract(hg19_target_coords, "(?<=-)[0-9]+"))
  )

# Join the parsed combined data with the combined TSS controls by matching coordinates.
combined_joined <- combined_parsed %>%
  left_join(
    guide_targets_all %>% select(target_name, target_type, target_chr, target_start, target_end),
    by = c("hg19_target_chr" = "target_chr",
           "hg19_target_start" = "target_start",
           "hg19_target_end" = "target_end")
  )


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(combined_joined, snakemake@output$results_with_design_file_features)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)