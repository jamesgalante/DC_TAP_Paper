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
k562_summary <- read_tsv(snakemake@input$k562_summary) %>% mutate(cell_type = "K562")
wtc11_summary <- read_tsv(snakemake@input$wtc11_summary) %>% mutate(cell_type = "WTC11")


### COMBINING THE TABLES ======================================================

full_summary <- rbind(k562_summary, wtc11_summary)
colnames(full_summary) <- c("Significant_Upregulated", "Significant_Downregulated", "Non_Significant_Well_Powered", "Non_Significant_Under_Powered", "Cell_Type")


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(full_summary, snakemake@output$summary_of_element_gene_categories_supplementary_table)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)