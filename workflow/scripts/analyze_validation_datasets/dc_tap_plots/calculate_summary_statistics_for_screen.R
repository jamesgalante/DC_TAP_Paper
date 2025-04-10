# Script: calculate_summary_statistics_for_screen.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/calculate_summary_statistics_for_screen.rda"))
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
combined_validation <- read_tsv(snakemake@input$combined_validation) %>% filter(pred_id == "ENCODE_rE2G")

# Load in the guide_targets files
k562_guide_targets <- read_tsv(snakemake@input$guide_targets[[1]])
wtc11_guide_targets <- read_tsv(snakemake@input$guide_targets[[2]])

# Load in the training dataset
encode_training_dataset <- read_tsv(snakemake@input$encode_training_dataset)

# Load in the create_ensemble_encode input files
create_ensemble_encode_input <- bind_rows(lapply(snakemake@input$create_ensemble_encode_input, read_tsv))


### ADD TSS CONTROL INFORMATION ===============================================

# Combine the guide_targets and the combined_validation by their coordinates (hg19) for each target
# Combine the guide targets
guide_targets_all <- bind_rows(
  k562_guide_targets,
  wtc11_guide_targets
) %>%
  group_by(target_name) %>%
  dplyr::slice(1)

# Parse combined_validation: extract coordinate information
combined_parsed <- combined_validation %>%
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

# For all TSS Controls, get the target gene name - separate by "_", scrap anything that's just a number, "TSS", "AND", or has "chr" in it
combined_joined <- combined_joined %>%
  mutate(TSS_control_gene = case_when(
    target_type %in% c("tss_pos", "tss_random") ~ 
      sapply(strsplit(target_name, "_"), function(x) 
        paste(unique(x[!x %in% c("TSS", "AND") & !grepl("^\\d+$", x) & !grepl("^chr", x)]), 
              collapse = "_")),
    TRUE ~ NA_character_
  ))


### FIGURING OUT DISTAL ELEMENT POSITIVE CONTROL TARGET GENES =================

# There are a few distal element positive controls. We want to figure out what those controls' intended targets were
# Get all positive control distal elements from K562 because these aren't labelled with their tested control gene
positive_control_distal_elements_k562 <- combined_joined %>% filter(target_type == "DE", Reference == "K562_DC_TAP_Seq")

# Import list of all control genes used in K562 screen
control_genes <- c("TMEM98", "RASSF7", "RAE1", "PLP2", "PIM2", "GATA1", "UBR5", "LYL1", "DNASE2", "PPIF", "CEP104", "CTSD", "NFE2", "PHF20L1", "LRRC47", "RRP8", 
                   "TIMM10B", "RBM38", "LRRCC1", "CD164", "MYC", "PIM1", "FAM83A", "FADS1", "ALAS2", "EPB41", "TMEM41B", "JUNB", "UROS", "SVIP", "HBD", "CCDC26", "SMIM1")

# Retrieve all DEs that are paired with a control_gene
t <- positive_control_distal_elements_k562 %>% filter(measuredGeneSymbol %in% control_genes)

# For all DE targets that are paired with one gene after filtering, this is the intended control gene - save these
easy_des <- names((t %>% pull(target_name) %>% table())[t %>% pull(target_name) %>% table() == 1])

# Create a map for these easy DE targets
easy_map_df <- t %>% 
  filter(target_name %in% easy_des) %>% 
  select(target_name, measuredGeneSymbol) %>%
  dplyr::rename(de_assigned_gene = measuredGeneSymbol)

# After looking at IGV, I've estimated the other DE positive controls
manual_mappings <- list(
  # chr1:3691430-3691731 is closest to SMIM1
  "chr1:3691430-3691731" = "SMIM1",
  
  # These are all closer to HBD
  "chr11:5297067-5297368" = "HBD",
  "chr11:5301767-5302068" = "HBD",
  "chr11:5305872-5306173" = "HBD",
  "chr11:5309369-5309670" = "HBD",
  
  # These are closest to JUNB
  "chr19:12895376-12895677" = "JUNB",
  "chr19:12900763-12901064" = "JUNB",
  
  # This one is closest to LYL1
  "chr19:13215352-13215653" = "LYL1",
  
  # This one is closest to RBM38
  "chr20:55990343-55990644" = "RBM38",
  
  # These are MYC
  "chr8:128911048-128911349" = "MYC", 
  "chr8:128972548-128972849" = "MYC",
  
  # These are CCDC26
  "chr8:130594299-130594600" = "CCDC26",
  "chr8:130701786-130702087" = "CCDC26",
  "chr8:130705120-130705421" = "CCDC26",
  
  # These are closest to GATA1
  "chrX:48641339-48641640" = "GATA1",
  "chrX:48658921-48659222" = "GATA1",
  
  # These are closest to PIM2
  "chrX:48797830-48798131" = "PIM2",
  "chrX:48799012-48799313" = "PIM2",
  
  # These are closest to PLP2
  "chrX:49005070-49005371" = "PLP2",
  "chrX:49022767-49023068" = "PLP2",
  "chrX:49023638-49023939" = "PLP2"
)

# Create a mapping for these other DE positive control target - Gene pairs
manual_mappings_df <- tibble(
  target_name = names(manual_mappings),
  de_assigned_gene = unlist(manual_mappings)
)

# Combine both dataframes
de_gene_map_df <- bind_rows(easy_map_df, manual_mappings_df)

# Add the intended positive control gene for each distal element target
combined_joined <- combined_joined %>%
  left_join(de_gene_map_df, by = "target_name") %>%
  mutate(de_assigned_gene = case_when(
    target_type == "DE" ~ de_assigned_gene,
    TRUE ~ NA
  )) %>%
  # Now add the WTC11 DE targets - which have the intended positive control gene indicated in the target_name
  mutate(de_assigned_gene = case_when(
    Reference == "WTC11_DC_TAP_Seq" & target_type == "DE" ~ sapply(strsplit(target_name, "_"), function(x) 
      paste(unique(x[x != "DE" & !grepl("^\\d+$", x)]), 
            collapse = "_")),
    TRUE ~ de_assigned_gene
  )) %>%
  # For the two de_assigned_genes that are OCT4 for WTC11, change these to `POU5F1`
  mutate(de_assigned_gene = ifelse(ExperimentCellType == "WTC11" & measuredGeneSymbol == "POU5F1" & target_type == "DE", "POU5F1", de_assigned_gene))


### LABELLING EACH PAIR BY CATEGORY ===========================================

# E-G table: We want to label each pair as T/F with regard to each of the following categories
  # Category: `DistalElement-Gene`
  # Category: `DistalPromoter-Gene`
  # Category: `SelfPromoter`
  # Category: `Positive Control DistalElement-Gene`
  # Category: `Random DistalElement-Gene`
  # Category: `Random Validation DistalElement-Gene`

# Define vector of ENCODE filter flags
ValidConnection_Flags <- c("TSS targeting guide\\(s\\)", "distance <1000", "overlaps target gene exon", "overlaps potential promoter", "overlaps target gene intron")

# Category: `DistalElement-Gene`
# These are all pairs that pass ENCODE filters - non-significant pairs that pass ENCODE filters, but are underpowered are included in this set
combined_joined_w_categories <- combined_joined %>%
  mutate(DistalElement_Gene = ifelse(!str_detect(ValidConnection, str_c(ValidConnection_Flags, collapse = "|")), TRUE, FALSE))                 

# There are two cases that we want to account for
  # When there's "overlaps potential promoter" - what promoter is it overlapping -> is that promoter the same as the measuredGeneSymbol?
  # "TSS targeting guide(s)" - In this case is the measuredGeneSymbol == TSS_control_gene?
# When the ValidConnection parameter includes, "overlaps potential promoter", extract all promoters (if there are multiple) and put into list
combined_joined_w_categories <- combined_joined_w_categories %>%
  rowwise() %>%  # Process each row individually
  mutate(PromoterGeneList = if (str_detect(ValidConnection, "overlaps potential promoter:")) {
    # Extract the substring, split on "|" and trim each element,
    # then wrap in a list so the column remains a list-column.
    list(trimws(unlist(str_split(
      str_extract(ValidConnection, "(?<=overlaps potential promoter:)[^;]+"),
      "\\|"
    ))))
  } else {
    list(NA_character_)
  }) %>%
  ungroup()

# For each target that is overlapping a promoter or has TSS targeting guide(s), let's define a DistalPromoter_Gene column for when that target is tested against a gene that it isn't overlapping
# We prioritize the results returned by "overlaps potential promoter" for deciding if it's selfPromoter or DistalPromoter_Gene. See edge cases for when these don't align: https://docs.google.com/presentation/d/10NFOre0GwunuMbL625IgyWAsmY2uE3rBhf7SY4LbRZM/edit?usp=sharing
# There's only one instance (if we prioritize "overlap") where we have to manually change DistalPromoter_Gene to selfPromoter: name == "FAM83A|chr8:124192121-124192621:."
combined_joined_w_categories <- combined_joined_w_categories %>%
  # Create a helper column that picks the promoter type
  mutate(
    promoter_type = case_when(
      str_detect(ValidConnection, "overlaps potential promoter:") ~ "overlap",
      str_detect(ValidConnection, "TSS targeting") ~ "tss",
      TRUE ~ "none"
    ),
    # When both conditions are present, use "overlap" as the primary deciding factor
    promoter_type = if_else(
      str_detect(ValidConnection, "overlaps potential promoter:") & str_detect(ValidConnection, "TSS targeting"),
      "overlap",
      promoter_type
    )
  ) %>%
  mutate(
    selfPromoter = case_when(
      promoter_type == "overlap" ~ map2_lgl(measuredGeneSymbol, PromoterGeneList, ~ .x %in% .y),
      promoter_type == "tss" ~ (measuredGeneSymbol == TSS_control_gene),
      TRUE ~ FALSE
    ),
    DistalPromoter_Gene = case_when(
      promoter_type == "overlap" ~ map2_lgl(measuredGeneSymbol, PromoterGeneList, ~ !(.x %in% .y)),
      promoter_type == "tss" ~ (measuredGeneSymbol != TSS_control_gene),
      TRUE ~ FALSE
    )
  ) %>%
  # Change that one value: for the specific FAM83A target, set DistalPromoter_Gene to FALSE and selfPromoter to TRUE.
  mutate(
    DistalPromoter_Gene = if_else(name == "FAM83A|chr8:124192121-124192621:.", FALSE, DistalPromoter_Gene),
    selfPromoter = if_else(name == "FAM83A|chr8:124192121-124192621:.", TRUE, selfPromoter)
  )

# Category: `Positive Control DistalElement-Gene`
combined_joined_w_categories <- combined_joined_w_categories %>%
  mutate(Positive_Control_DistalElement_Gene = ifelse(target_type == "DE" & de_assigned_gene == measuredGeneSymbol, TRUE, FALSE))

# Category: `Random DistalElement-Gene`
# For the random set, we can define that as every pair where targets are from the random 2Mb loci (except that one instance of NANOG DE pos control being within a 2Mb region in wtc11)
# This should just be all "enh"
combined_joined_w_categories <- combined_joined_w_categories %>%
  mutate(Random_DistalElement_Gene = ifelse(target_type == "enh" & DistalElement_Gene == TRUE, TRUE, FALSE))

# Category: `Random Validation DistalElement-Gene`
# Random element pairs that are also DistalElement pairs
# Create GRanges objects for overlap detection (without merging)
random_gr <- combined_joined_w_categories %>%
  mutate(random_pair_id = row_number()) %>%
  filter(ExperimentCellType == "K562", Random_DistalElement_Gene == TRUE) %>%
  { GRanges(
    seqnames = paste0(.$chrom, ":", .$measuredGeneSymbol),
    ranges = IRanges(start = .$chromStart, end = .$chromEnd),
    random_pair_id = .$random_pair_id
  ) }

training_gr <- encode_training_dataset %>%
  filter(pred_id == "ENCODE_rE2G") %>%
  { GRanges(
    seqnames = paste0(.$chrom, ":", .$measuredGeneSymbol),
    ranges = IRanges(start = .$chromStart, end = .$chromEnd),
  ) }

# Find overlaps (since the seqnames include the gene name, overlapping hits should match on the gene as well)
ovl_hits <- findOverlaps(random_gr, training_gr)

# Get the list of random_pair_ids with matching overlaps
matching_overlap_ids <- random_gr$random_pair_id[queryHits(ovl_hits)]

# Mark as TRUE for those NOT overlapping training (i.e. valid for Random Validation)
combined_joined_w_categories <- combined_joined_w_categories %>%
  mutate(random_pair_id = row_number()) %>%
  mutate(Random_Validation_DistalElement_Gene = case_when(
    ExperimentCellType == "K562" ~ random_pair_id %in% setdiff(random_gr$random_pair_id, matching_overlap_ids),
    ExperimentCellType == "WTC11" ~ Random_DistalElement_Gene
  )) %>%
  select(-random_pair_id)


### CALCULATE SUMMARY STATISTICS ==============================================

# For each category, calculate the number of Significant pairs (upregulated v. downregulated) - calculate the the number of Non-significant pairs (high power v. underpowered)
get_significant_pair_stats <- function(df, category) {
  df %>%
    filter(!!sym(category), Significant == T) %>%
    mutate(Downregulated = EffectSize < 0) %>%
    select(Significant, Downregulated) %>%
    table()
}
get_nonsignificant_pair_stats <- function(df, category) {
  df %>%
    filter(!!sym(category), Significant == F) %>%
    mutate(Underpowered = str_detect(ValidConnection, "< 80% power at 15% effect size")) %>%
    select(Significant, Underpowered) %>%
    table()
}

# List of categories to summarize
categories <- c(
  "DistalElement_Gene",
  "selfPromoter",
  "DistalPromoter_Gene",
  "Positive_Control_DistalElement_Gene",
  "Random_DistalElement_Gene",
  "Random_Validation_DistalElement_Gene"
)

# Function to build summary table for a given cell type
get_summary_table_by_celltype <- function(cell_type) {
  df_cell <- combined_joined_w_categories %>% filter(ExperimentCellType == cell_type)
  
  map_dfr(categories, function(cat) {
    sig_tbl <- get_significant_pair_stats(df_cell, cat)
    nonsig_tbl <- get_nonsignificant_pair_stats(df_cell, cat)
    
    # Safely retrieve table values or 0 if not present
    downreg_false <- if ("TRUE" %in% rownames(sig_tbl) && "FALSE" %in% colnames(sig_tbl)) {
      sig_tbl["TRUE", "FALSE"]
    } else 0
    
    downreg_true <- if ("TRUE" %in% rownames(sig_tbl) && "TRUE" %in% colnames(sig_tbl)) {
      sig_tbl["TRUE", "TRUE"]
    } else 0
    
    underpowered_false <- if ("FALSE" %in% rownames(nonsig_tbl) && "FALSE" %in% colnames(nonsig_tbl)) {
      nonsig_tbl["FALSE", "FALSE"]
    } else 0
    
    underpowered_true <- if ("FALSE" %in% rownames(nonsig_tbl) && "TRUE" %in% colnames(nonsig_tbl)) {
      nonsig_tbl["FALSE", "TRUE"]
    } else 0
    
    tibble(
      Category = cat,
      Upregulated = downreg_false,
      Downregulated = downreg_true,
      WellPowered = underpowered_false,
      UnderPowered = underpowered_true
    )
  }) %>% mutate(ExperimentCellType = cell_type)
}

# Generate summary tables for each cell type
summary_K562 <- get_summary_table_by_celltype("K562")
summary_WTC11 <- get_summary_table_by_celltype("WTC11")

# Print the summary tables
cat("Summary table for K562:\n")
print(summary_K562)

cat("\nSummary table for WTC11:\n")
print(summary_WTC11)


### REMOVE UNDERPOWERED PAIRS =================================================

# We wanted to get the statistics on how many pairs are underpowered, but we don't want this in the final table
# Modify the categories to remove underpowered nonsignificant pairs
combined_joined_w_categories_fixed <- combined_joined_w_categories %>%
  mutate(
    DistalElement_Gene = case_when(
      Significant == FALSE & 
        str_detect(ValidConnection, "< 80% power at 15% effect size") ~ FALSE,
      TRUE ~ DistalElement_Gene
    ),
    selfPromoter = case_when(
      Significant == FALSE & 
        str_detect(ValidConnection, "< 80% power at 15% effect size") ~ FALSE,
      TRUE ~ selfPromoter
    ),
    DistalPromoter_Gene = case_when(
      Significant == FALSE & 
        str_detect(ValidConnection, "< 80% power at 15% effect size") ~ FALSE,
      TRUE ~ DistalPromoter_Gene
    ),
    Positive_Control_DistalElement_Gene = case_when(
      Significant == FALSE & 
        str_detect(ValidConnection, "< 80% power at 15% effect size") ~ FALSE,
      TRUE ~ Positive_Control_DistalElement_Gene
    ),
    Random_DistalElement_Gene = case_when(
      Significant == FALSE & 
        str_detect(ValidConnection, "< 80% power at 15% effect size") ~ FALSE,
      TRUE ~ Random_DistalElement_Gene
    ),
    Random_Validation_DistalElement_Gene = case_when(
      Significant == FALSE & 
        str_detect(ValidConnection, "< 80% power at 15% effect size") ~ FALSE,
      TRUE ~ Random_Validation_DistalElement_Gene
    )
  )


### CREATE SUMMARIZED COMBINED JOINED W CATEGORIES ============================

# Add the guide names for the guide targets to create column 2 of the IGVF formatted file
target_names_w_guides <- rbind(
  k562_guide_targets %>% 
    filter(!target_type %in% c("safe_targeting", "negative_control")) %>%
    group_by(target_name) %>% 
    summarise(all_names = paste(name, collapse = ","), .groups = "drop") %>%
    mutate(cell_type = "K562"),
  wtc11_guide_targets %>%
    filter(!target_type %in% c("safe_targeting", "negative_control")) %>%
    group_by(target_name) %>% 
    summarise(all_names = paste(name, collapse = ","), .groups = "drop") %>%
    mutate(cell_type = "WTC11")
)

# Create IGVF formatted file
igvf_formatted_file <- combined_joined_w_categories_fixed %>%
  left_join(
    create_ensemble_encode_input %>% 
      select(chrom, chromStart, chromEnd, measuredGeneSymbol, Reference, measuredEnsemblID, pValue), 
    by = c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "Reference")) %>%
  left_join(
    target_names_w_guides,
    by = c("target_name", "ExperimentCellType" = "cell_type")
  ) %>%
  mutate(
    intended_target_name = paste0(chrom, ":", chromStart, "-", chromEnd),
    `guide_id(s)` = as.character(all_names),
    targeting_chr = as.character(chrom),
    targeting_start = as.integer(chromStart),
    targeting_end = as.integer(chromEnd),
    type = as.character(target_type),
    gene_id = as.character(measuredEnsemblID),
    gene_symbol = as.character(measuredGeneSymbol),
    sceptre_log2_fc = as.numeric(log2(EffectSize + 1)),
    sceptre_p_value = as.numeric(pValue),
    sceptre_adj_p_value = as.numeric(pValueAdjusted),
    significant = as.logical(Significant),
    sample_term_name = as.character(ExperimentCellType),
    sample_term_id = as.character(NA_character_),
    sample_summary_short = as.character(ExperimentCellType),
    power_at_effect_size_10 = as.numeric(PowerAtEffectSize10),
    power_at_effect_size_15 = as.numeric(PowerAtEffectSize15),
    power_at_effect_size_20 = as.numeric(PowerAtEffectSize20),
    power_at_effect_size_25 = as.numeric(PowerAtEffectSize25),
    power_at_effect_size_50 = as.numeric(PowerAtEffectSize50),
    notes = as.character(NA_character_)
  ) %>%
  select(intended_target_name, `guide_id(s)`, targeting_chr, targeting_start, targeting_end, type, gene_id, gene_symbol, sceptre_log2_fc,
         sceptre_p_value, sceptre_adj_p_value, significant, sample_term_name, sample_term_id, sample_summary_short, power_at_effect_size_10, 
         power_at_effect_size_15, power_at_effect_size_20, power_at_effect_size_25, power_at_effect_size_50, notes)

# Add distance to TSS, the categories, original target names, and hg19 coordinates
summarized_categories <- igvf_formatted_file %>%
  left_join(
    create_ensemble_encode_input %>%
      mutate(CellType = str_extract(Reference, "^[^_]+")) %>%  # Extract string before first "_"
      select(chrom, chromStart, chromEnd, measuredGeneSymbol, CellType, distToTSS),
    by = c("targeting_chr" = "chrom", "targeting_start" = "chromStart", "targeting_end" = "chromEnd", 
           "gene_symbol" = "measuredGeneSymbol", "sample_term_name" = "CellType")
  ) %>%
  left_join(
    combined_joined_w_categories_fixed %>%
      select(chrom, chromStart, chromEnd, measuredGeneSymbol, ExperimentCellType, 
             hg19_target_coords, EffectSize, chrTSS, startTSS, endTSS, hg19_target_chr, hg19_target_start, 
             hg19_target_end, target_name, DistalElement_Gene, selfPromoter, DistalPromoter_Gene, ValidConnection, 
             Positive_Control_DistalElement_Gene, Random_DistalElement_Gene, Random_Validation_DistalElement_Gene),
    by = c("targeting_chr" = "chrom", "targeting_start" = "chromStart", "targeting_end" = "chromEnd", 
           "gene_symbol" = "measuredGeneSymbol", "sample_term_name" = "ExperimentCellType")
  ) %>%
  dplyr::rename(
    intended_target_name_hg38 = intended_target_name,
    targeting_chr_hg38 = targeting_chr,
    targeting_start_hg38 = targeting_start,
    targeting_end_hg38 = targeting_end,
    intended_target_name_hg19 = hg19_target_coords,
    log2FC_EffectSize = sceptre_log2_fc,
    pctChange_EffectSize = EffectSize,
    chrTSS_hg38 = chrTSS,
    startTSS_hg38 = startTSS,
    endTSS_hg38 = endTSS,
    targeting_chr_hg19 = hg19_target_chr,
    targeting_start_hg19 = hg19_target_start,
    targeting_end_hg19 = hg19_target_end,
    design_file_target_name = target_name
  )


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(combined_joined_w_categories_fixed, snakemake@output$combined_joined_w_categories)
write_tsv(summarized_categories, snakemake@output$summarized_categories)
write_tsv(igvf_formatted_file, snakemake@output$igvf_formatted_file)
write_tsv(summary_K562, snakemake@output$summary_K562)
write_tsv(summary_WTC11, snakemake@output$summary_WTC11)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)
