# Script: create_gasperini_sceptre_inputs.R

### SETUP =====================================================================

# Saving image for debugging
save.image("RDA_objects/create_gasperini_sceptre_inputs.rda")
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

# required packages
message("Loading packages")
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(rtracklayer)
  library(readr)
  library(tidyverse)
})

# Load input files
message("Loading input files")
perturb_sce <- readRDS(snakemake@input$perturb_sce)
crispr_pipeline_output <- read_tsv(snakemake@input$crispr_pipeline_output)


### CREATE SCEPTRE INPUT FILES ================================================

# Extract the raw counts file from the perturb_sce rds object
raw_counts <- assay(perturb_sce, "counts")

# Extract the binarized guide counts matrix
binarized_guide_counts <- assay(altExp(perturb_sce, "grna_perts"), "perts")

# Create the response id - target pairs file
response_id_target_pairs <- crispr_pipeline_output %>% select(perturbation, gene, target_type)

# Create the guide - target pairs file
guide_target_pairs <- as.data.frame(rowData(altExp(perturb_sce, "grna_perts"))) %>% dplyr::select(name, target_name)
rownames(guide_target_pairs) <- NULL


### FILTERING =================================================================

# Are all genes in the response id - target pairs file in the original raw_counts matrix
all(response_id_target_pairs$gene %in% rownames(raw_counts))
# yes

#### ==========

# Are all perturbations in the response id - target pairs file in the guide - target pairs file
all(response_id_target_pairs$perturbation %in% guide_target_pairs$target_name)
# no

# Which perturbations are in the response-target file but not in the guide-target file
unique(response_id_target_pairs %>% filter(!perturbation %in% guide_target_pairs$target_name) %>% select(perturbation))
# chr8:101912497-101912997

# Let's remove this from the response_id_target_pairs file 
# If Sceptre tries to test a perturbation-gene pair that doesn't exist, it will throw an error
# Let's also remove all positive controls from the response_id_target_pairs file and subset the response_id and target columns
# We can also rename the columns to "grna_target" and "response_id" for Sceptre
response_id_target_pairs <- response_id_target_pairs %>% 
  filter(perturbation != "chr8:101912497-101912997") %>%
  filter(target_type == "enh") %>% 
  select(perturbation, gene) %>%
  dplyr::rename(grna_target = perturbation, response_id = gene)

#### ==========
  
# Are all guides in guide_target_pairs unique
length(unique(guide_target_pairs$name)) == nrow(guide_target_pairs)
# yes

# We can also rename the guide - target dataframe for Sceptre
guide_target_pairs <- guide_target_pairs %>%
  dplyr::rename(grna_id = "name", grna_target = "target_name")

### SAVE OUTPUT ===============================================================

# Save the raw counts
saveRDS(raw_counts, snakemake@output$raw_counts)

# Save the binarized guide counts
saveRDS(binarized_guide_counts, snakemake@output$binarized_guide_counts)

# Save the response id - target pairs file
write_tsv(response_id_target_pairs, snakemake@output$response_id_target_pairs)

# Save the guide - target pairs file
write_tsv(guide_target_pairs, snakemake@output$guide_target_pairs)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

