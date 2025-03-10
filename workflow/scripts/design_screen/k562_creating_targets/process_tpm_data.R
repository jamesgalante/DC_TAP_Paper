# Script: process_tpm_data.R
### SETUP =====================================================================
# Saving image for debugging
if (!file.exists("RDA_objects/k562_creating_targets")) { dir.create("RDA_objects/k562_creating_targets", recursive = TRUE) }
save.image(paste0("RDA_objects/k562_creating_targets/process_tpm_data.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")

### LOADING FILES =============================================================
# Required packages and functions
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

# Load input files
message("Loading input files...")
tpm_file <- snakemake@input$tpm_file
tpm <- read_csv(tpm_file, col_types = cols(.default = col_character(), tpm = col_double()))

# Get parameters from Snakemake
tpm_threshold <- snakemake@params$tpm_threshold

### TPM FILTERING ============================================================
message(paste("Filtering genes with TPM >=", tpm_threshold))
# Count number and percentage of genes above TPM threshold
n_above_threshold <- sum(tpm$tpm >= tpm_threshold)
message(paste("Number of genes above", tpm_threshold, "TPM:", n_above_threshold))

# Filter genes above TPM threshold
tpm_filtered <- tpm %>% filter(tpm >= tpm_threshold)

### VISUALIZATION ============================================================
message("Creating TPM distribution plot...")
tpm_plot <- ggplot(tpm, aes(tpm)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = tpm_threshold, color = "red") +
  labs(
    x = "Transcripts-per-million (TPM)", 
    title = paste0("Genes above ", tpm_threshold, " TPM: ", n_above_threshold)
  ) +
  scale_x_log10() +
  theme_bw()

### SAVE OUTPUT ==============================================================
message("Saving output files")
saveRDS(tpm_filtered, file = snakemake@output$tpm_filtered)
ggsave(snakemake@output$tpm_plot, plot = tpm_plot, width = 7, height = 5)

### CLEAN UP =================================================================
message("Closing log file")
sink()
sink(type = "message")
close(log)