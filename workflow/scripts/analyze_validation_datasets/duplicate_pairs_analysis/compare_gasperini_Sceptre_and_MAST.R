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
# There are none


### ADD ENSEMBLE ID TO TRAINING ===============================================

# Add response_id to combined_training by mapping measuredGeneSymbol to gene_name
combined_training_mapped <- combined_training %>%
  left_join(gene_mapping, by = c("measuredGeneSymbol" = "gene_name"))

# Check for any NAs in response_id
missing_response_ids <- combined_training_mapped %>%
  filter(is.na(response_id)) %>%
  distinct(measuredGeneSymbol)
# There are none


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


### CREATE THE SCEPTRE VERSION ONLY ===========================================

# When merging with the combined_training, we're essentially deleting any sceptre pairs that don't match
# To avoid this and keep all Sceptre pairs, let's just merge the gasperini sceptre results with selected columns from merged_df_all so that we can get the hg38 conversion
# And let's convert the Log2FC to pctChange
gasperini_sceptre_results_w_symbol <- gasperini_sceptre_results %>%
  left_join(merged_df_all %>% 
              select(chrom, chromStart, chromEnd, grna_target) %>%
              group_by(grna_target) %>%
              dplyr::slice(1), 
            by = "grna_target") %>%
  left_join(gene_mapping, by = "response_id") %>%
  mutate(Sceptre_pctChange = (2^log_2_fold_change - 1) * 100) %>%
  filter(!is.na(chrom), !is.na(Sceptre_pctChange))


### SAVE OUTPUT ===============================================================

# Write the merged object for later plotting
saveRDS(merged_df_all, snakemake@output$gasperini_MAST_and_Sceptre)
saveRDS(gasperini_sceptre_results_w_symbol, snakemake@output$gasperini_sceptre_results_w_symbol)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

