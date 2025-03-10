# Script: analyze_gasperini_with_sceptre.R

### SETUP =====================================================================

# Saving image for debugging
save.image("RDA_objects/compare_gasperini_Sceptre_and_MAST.rda")
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
  library(tidyverse)
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(rtracklayer)
})

# Load input files
message("Loading input files")
gasperini_sceptre_results <- readRDS(snakemake@input$gasperini_sceptre_results)
combined_training <- read_tsv(snakemake@input$combined_training) %>% filter(Dataset == "Gasperini2019")
annot <- import(snakemake@input$annotation_file)


### FORMAT GENE MAPPING =======================================================

# Get associated gene symbol for each ensemble id
# extract gene locus annotations and only retain autosomes
genes <- annot[annot$type == "gene"]
genes <- genes[seqnames(genes) %in% c(paste0("chr", 1:22), "chrX")]

# Format gene mapping from ENSEMBLE ID to Gene Symbol
gene_mapping <- genes %>%
  as.data.frame() %>%
  select(gene_name, gene_id) %>%
  mutate(gene_id = sub("\\..+", "", gene_id)) %>%
  filter(gene_id %in% gasperini_sceptre_results$response_id) %>%
  dplyr::rename(response_id = "gene_id")

# Check for duplicates in gene_symbol to ENSEMBL ID mapping
duplicates <- gene_mapping %>%
  group_by(gene_name) %>%
  filter(n() > 1) %>%
  ungroup()

# Check if any of the gene names are represented in the Combined Training data
if(any(duplicates$gene_name %in% combined_training$measuredGeneSymbol)){
  warning("Some gene symbols map to multiple ENSEMBL IDs. These will result in multiple rows in the merged dataframe.")
  print(duplicates)
}


### ADD ENSEMBLE ID TO TRAINING ===============================================

# Add response_id to combined_training by mapping measuredGeneSymbol to gene_name
combined_training_mapped <- combined_training %>%
  left_join(gene_mapping, by = c("measuredGeneSymbol" = "gene_name"))

# Check for any NAs in response_id
missing_response_ids <- combined_training_mapped %>%
  filter(is.na(response_id)) %>%
  distinct(measuredGeneSymbol)

if(nrow(missing_response_ids) > 0){
  warning("Some measuredGeneSymbols in combined_training did not have a corresponding ENSEMBL ID in perturb_sce.")
  print(missing_response_ids)
}


### MERGED AND FORMAT =========================================================

# Extract pert_name using regex
combined_training_mapped <- combined_training_mapped %>%
  mutate(grna_target = gsub("^[^|]*\\|([^:]*:[^:]*-[^:]*)[:].*", "\\1", name))

# Merge the datasets
merged_df_all <- combined_training_mapped %>%
  left_join(gasperini_sceptre_results, by = c("grna_target" = "grna_target", "response_id" = "response_id"))

# Convert log2_fold_change and EffectSize to Percent Change
merged_df_all <- merged_df_all %>%
  mutate(
    Sceptre_pctChange = (2^log_2_fold_change - 1) * 100,
    MAST_pctChange = EffectSize * 100
  )


### SAVE OUTPUT ===============================================================

# Write the merged object for later plotting
saveRDS(merged_df_all, snakemake@output$gasperini_MAST_and_Sceptre)

### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

