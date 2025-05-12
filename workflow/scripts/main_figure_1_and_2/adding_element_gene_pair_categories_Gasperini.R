# Script: adding_element_gene_pair_categories_Gasperini.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/adding_element_gene_pair_categories_Gasperini.rda"))
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
combined_joined <- read_tsv(snakemake@input$results_with_design_file_and_genomic_feature_overlaps)
create_ensemble_encode_input <-read_tsv(snakemake@input$create_ensemble_encode_input)
guide_targets <- read_tsv(snakemake@input$guide_targets)
discovery_results <- readRDS(snakemake@input$discovery_results)


### PREPARE BASE DATA =======================================================

message("Preparing base data")
# Handle NA values in ValidConnection and create element location classifier
combined_joined <- combined_joined %>%
  mutate(
    ValidConnection = ifelse(is.na(ValidConnection), "", ValidConnection),
    distal_or_promoter = ifelse(
      is.na(abc_tss_overlap) & is.na(gencode_promoter_overlap), 
      "distal", "promoter"
    )
  )


### ADD ELEMENT-GENE PAIR CATEGORIES =========================================

message("Adding element-gene pair categories")
# Define filter flags for distal elements
ValidConnection_Flags <- c("distance <1000", "overlaps target gene exon", "overlaps target gene intron")

# Add all category flags in a single step with detailed descriptions
combined_joined_w_categories <- combined_joined %>%
  mutate(
    # DistalElement-Gene: 
    # True if the element is a distal regulatory element (not overlapping any promoter) 
    # AND it passes all ENCODE filters (not overlapping target gene exon/intron, not too close to TSS)
    # These are candidate enhancers that are not promoters and don't overlap their target genes
    DistalElement_Gene = !str_detect(
      ValidConnection, 
      str_c(ValidConnection_Flags, collapse = "|")
    ) & (distal_or_promoter == "distal"),
    
    # DistalPromoter-Gene:
    # True if the element overlaps any promoter (in ABC, GENCODE, or defined TSS controls)
    # BUT does NOT overlap the promoter of its tested target gene
    # These are elements that might function as promoters for one gene while regulating another
    DistalPromoter_Gene = (distal_or_promoter == "promoter") & !measuredGene_in_abc_tss & !measuredGene_in_gencode_promoter,
    
    # SelfPromoter:
    # True if the element overlaps the promoter of its tested target gene
    # (in ABC, GENCODE)
    # These elements represent the gene's own promoter regulating the gene
    selfPromoter = measuredGene_in_abc_tss | measuredGene_in_gencode_promoter
  )


### PREPARE GUIDE INFORMATION ================================================

message("Preparing guide information")
# Compile guide information for each target
target_names_w_guides <- guide_targets %>% 
    group_by(target_name) %>% 
    summarise(all_names = paste(name, collapse = ","), .groups = "drop") %>%
    mutate(cell_type = "K562")


### PREPARE ENCODE DATA =======================================================

message("Preparing ENCODE data for joining")
# Pre-process the encode data to include all needed columns in one dataframe
encode_data_for_join <- create_ensemble_encode_input %>%
  # Select all needed columns for a single join
  select(
    chrom, chromStart, chromEnd, 
    measuredGeneSymbol, Reference, 
    measuredEnsemblID, pValue, distToTSS
  )


### CREATE FINAL OUTPUT =====================================================

message("Creating final output")
# Process data to create final output with a single encode join
final_output <- combined_joined_w_categories %>%
  # Add all ENCODE data in a single join
  left_join(
    encode_data_for_join,
    by = c(
      "chrom" = "chrom", 
      "chromStart" = "chromStart", 
      "chromEnd" = "chromEnd", 
      "measuredGeneSymbol" = "measuredGeneSymbol", 
      "Reference" = "Reference"
    )
  ) %>%
  # Add guide information
  left_join(
    target_names_w_guides,
    by = c("target_name", "CellType" = "cell_type")
  ) %>%
  # Add in the se_fold_change
  left_join(
    discovery_results %>% select(response_id, grna_target, p_value, fold_change, se_fold_change),
    by = c("target_name" = "grna_target", "measuredEnsemblID" = "response_id")
  ) %>%
  # Calculate percent change, create identifiers, and calculate standard errors and confidence intervals
  mutate(
    # Effect size calculations
    pct_change_effect_size = EffectSize * 100,
    log_2_FC_effect_size = log2(fold_change),
    fold_change_effect_size = fold_change,
    
    # Standard errors
    standard_error_pct_change = se_fold_change * 100,
    standard_error_log_2_FC = se_fold_change / (fold_change_effect_size * log(2)),
    standard_error_fold_change = se_fold_change,
    
    # Confidence intervals for fold change
    lower_CI_95_fold_change = fold_change_effect_size - 1.96 * standard_error_fold_change,
    upper_CI_95_fold_change = fold_change_effect_size + 1.96 * standard_error_fold_change,
    
    # Confidence intervals for log2 fold change
    lower_CI_95_log_2_FC = log_2_FC_effect_size - 1.96 * standard_error_log_2_FC,
    upper_CI_95_log_2_FC = log_2_FC_effect_size + 1.96 * standard_error_log_2_FC,
    
    # Confidence intervals for percent change
    lower_CI_95_pct_change = pct_change_effect_size - 1.96 * standard_error_pct_change,
    upper_CI_95_pct_change = pct_change_effect_size + 1.96 * standard_error_pct_change,
    
    # Identifiers
    intended_target_name_hg38 = paste0(chrom, ":", chromStart, "-", chromEnd),
    element_gene_pair_identifier_hg38 = paste0(measuredGeneSymbol, "|", intended_target_name_hg38),
    element_gene_pair_identifier_hg19 = paste0(measuredGeneSymbol, "|", hg19_target_coords),
    distance_to_gencode_gene_TSS = abs(distToTSS)
  ) %>%
  # Select and organize final columns
  select(
    # HG38 target coordinates
    intended_target_name_hg38, 
    targeting_chr_hg38 = chrom, 
    targeting_start_hg38 = chromStart, 
    targeting_end_hg38 = chromEnd,
    
    # HG19 target coordinates
    intended_target_name_hg19 = hg19_target_coords, 
    targeting_chr_hg19 = hg19_target_chr, 
    targeting_start_hg19 = hg19_target_start, 
    targeting_end_hg19 = hg19_target_end,
    
    # Target names
    element_gene_pair_identifier_hg38,
    element_gene_pair_identifier_hg19,
    
    # Gene and cell type information
    gene_symbol = measuredGeneSymbol, 
    gene_id = measuredEnsemblID, 
    cell_type = CellType,
    
    # Effect size calculations and confidence intervals
    fold_change_effect_size,
    log_2_FC_effect_size, 
    pct_change_effect_size,
    
    # Standard errors
    standard_error_fold_change,
    standard_error_log_2_FC,
    standard_error_pct_change,
    
    # Confidence intervals
    lower_CI_95_fold_change, 
    upper_CI_95_fold_change,
    lower_CI_95_log_2_FC, 
    upper_CI_95_log_2_FC,
    lower_CI_95_pct_change, 
    upper_CI_95_pct_change,
    
    # Statistical measures
    sceptre_p_value = pValue, 
    sceptre_adj_p_value = pValueAdjusted, 
    significant = Significant,
    
    # Distance to TSS
    distance_to_gencode_gene_TSS, 
    distance_to_abc_canonical_TSS,
    
    # TSS coordinates HG38
    chrTSS_hg38 = chrTSS, 
    startTSS_hg38 = startTSS, 
    endTSS_hg38 = endTSS,
    
    # Element feature overlap
    element_location = distal_or_promoter,
    gencode_promoter_overlap, 
    abc_tss_overlap,
    gencode_protein_coding_gene_body_overlap = protein_coding_gene_body_overlap,
    
    # Element-Gene category
    DistalElement_Gene, 
    DistalPromoter_Gene, 
    selfPromoter, 
    
    # Power Simulation Results
    power_at_effect_size_10 = PowerAtEffectSize10, 
    power_at_effect_size_15 = PowerAtEffectSize15, 
    power_at_effect_size_20 = PowerAtEffectSize20, 
    power_at_effect_size_25 = PowerAtEffectSize25, 
    power_at_effect_size_50 = PowerAtEffectSize50,
    
    # Screen Design Information
    design_file_target_name = target_name, 
    design_file_type = target_type, 
    guide_ids = all_names
  )


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(final_output, snakemake@output$results_with_element_gene_pair_categories)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)