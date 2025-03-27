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
  source(file.path(snakemake@scriptdir, "../../process_validation_datasets/sceptre_setup/gene_target_pairing_functions.R"))
})

# Load input files
message("Loading input files")
perturb_sce <- readRDS(snakemake@input$perturb_sce)
crispr_pipeline_output <- read_tsv(snakemake@input$crispr_pipeline_output)
annot <- import(snakemake@input$annot)


### CREATE SCEPTRE INPUT FILES ================================================

# Extract the raw counts file from the perturb_sce rds object
raw_counts <- assay(perturb_sce, "counts")

# Extract the binarized guide counts matrix
binarized_guide_counts <- assay(altExp(perturb_sce, "grna_perts"), "perts")

# Create the response id - target pairs file
response_id_target_pairs <- crispr_pipeline_output %>% select(perturbation, gene, target_type)

# Create the guide - target pairs file and filter out positive controls
guide_targets <- as.data.frame(rowData(altExp(perturb_sce, "grna_perts"))) %>% filter(target_type == "enh")
rownames(guide_targets) <- NULL
guide_target_pairs <- guide_targets %>% dplyr::select(name, target_name)


### CUSTOM RESPONSE_ID_TARGET_PAIRS ===========================================

# One option is to only use pairs that were tested in the original gasperini analysis
# Another option is to remake the pairs to be tested with a larger distance threshold (2Mb) rather than 1Mb

# Need to have some form of "guide_targets" to put into pipeline
output <- find_genes_near_targets(guide_targets = guide_targets, annotation_file = annot, gene_ids = rownames(raw_counts), max_distance = 2e6)
response_id_target_pairs <- output[[1]]
report <- output[[2]]


### FILTERING =================================================================

# Are all genes in the response id - target pairs file in the original raw_counts matrix
all(response_id_target_pairs$gene %in% rownames(raw_counts))
# yes

#### ==========

# Are all perturbations in the response id - target pairs file in the guide - target pairs file
all(response_id_target_pairs$perturbation %in% guide_target_pairs$target_name)
# no

# Which perturbations are in the response-target file but not in the guide-target file
unique(response_id_target_pairs %>% filter(!grna_group %in% guide_target_pairs$target_name) %>% select(grna_group))

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

