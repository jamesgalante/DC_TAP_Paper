# Script: create_sceptre_diffex_input_Gasperini.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/create_sceptre_diffex_input_Gasperini.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading packages")
suppressPackageStartupMessages({
  library(sceptre)
  library(data.table)
  library(rtracklayer)
  library(tidyverse)
  library(GenomicRanges)
  source(file.path(snakemake@scriptdir, "../process_validation_datasets/sceptre_setup/gene_target_pairing_functions.R"))
})

# Load input files
message("Loading input files")
annot <- import(snakemake@input$annot)
raw_counts <- readRDS(snakemake@input$raw_counts)
binarized_guide_counts <- readRDS(snakemake@input$binarized_guide_counts)
guide_targets <- read_tsv(snakemake@input$guide_targets)


### CREATE GENE_GRNA_GROUP_PAIRS ==============================================

# Create the gene-gRNA pairs using the find_genes_near_targets function
output <- find_genes_near_targets(
  guide_targets = guide_targets, 
  annotation_file = annot, 
  gene_ids = rownames(raw_counts), 
  max_distance = 2e6
)
gene_gRNA_group_pairs <- output[[1]]
distances_table <- output[[2]]

# Modify for SCEPTRE input specifications
gene_gRNA_group_pairs <- gene_gRNA_group_pairs %>%
  dplyr::select(response_id, grna_target = grna_group)


### CREATE GRNA_GROUPS_TABLE ==================================================

# Create the grna_groups_table which matches guide names to target names
gRNA_groups_table <- guide_targets %>% 
  dplyr::select(name, target_name) %>%
  dplyr::rename(grna_id = "name", grna_target = "target_name")


### CREATE THE METADATA FILE ==================================================

# Create the metadata indicating the batch for each cell barcode
cell_barcode <- colnames(raw_counts)

# Extract prep batch from cell barcode suffixes
# The first digit after the underscore indicates the prep batch (1 or 2)
suffixes <- sub(".*-\\d+_", "", cell_barcode)  # Extract the suffix
cell_batches <- as.factor(str_extract(suffixes, "^\\d+"))

# Create the metadata dataframe
metadata <- data.frame(row.names = cell_barcode, batch = cell_batches)


### CREATE SCEPTRE OBJECT =====================================================

# Create sceptre_object
sceptre_object <- import_data(
  response_matrix = raw_counts,
  grna_matrix = binarized_guide_counts,
  grna_target_data_frame = gRNA_groups_table,
  moi = "high",
  extra_covariates = metadata
)

# Because the pre-calculated grna_n_nonzero and grna_n_umis covariates are not independent (due to a binary grna_matrix) we must remove one from the formula
new_formula <- formula(~ log(response_n_nonzero) + log(response_n_umis) + log(grna_n_nonzero) + batch)

# Set analysis parameters
sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = gene_gRNA_group_pairs,
  side = "both",
  grna_integration_strategy = "union",
  formula_object = new_formula
)

print(sceptre_object)


### SAVE OUTPUT ===============================================================

message("Saving output files")
saveRDS(gene_gRNA_group_pairs, snakemake@output$gene_gRNA_group_pairs)
saveRDS(gRNA_groups_table, snakemake@output$gRNA_groups_table)
saveRDS(metadata, snakemake@output$metadata)
saveRDS(sceptre_object, snakemake@output$sceptre_diffex_input)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)