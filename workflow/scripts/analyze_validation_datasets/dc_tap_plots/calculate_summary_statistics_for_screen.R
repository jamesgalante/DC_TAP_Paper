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

# Load in abc canonical TSS file
abc_canonical_tss <- read_tsv(snakemake@input$abc_canonical_tss)


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


### ADDING ABC TSS OVERLAP COLUMN =============================================

# Subset for only protein_coding genes
abc_canonical_tss <- abc_canonical_tss %>% filter(gene_type == "protein_coding")

# Convert tibble to GRanges object
abc_canonical_tss_gr <- GRanges(
  seqnames = abc_canonical_tss$`#chr`,
  ranges = IRanges(start = abc_canonical_tss$start, end = abc_canonical_tss$end),
  strand = abc_canonical_tss$strand,
  name = abc_canonical_tss$name,
  score = abc_canonical_tss$score,
  Ensembl_ID = abc_canonical_tss$Ensembl_ID,
  gene_type = abc_canonical_tss$gene_type
)

# Resize to 1000bp (from 500bp) while fixing the center point
abc_canonical_tss_extended_gr <- resize(abc_canonical_tss_gr, width = 1000, fix = "center")

# Create a GRanges object with all elements
combined_joined_elements_gr <- GRanges(
  seqnames = combined_joined$chrom,
  ranges = IRanges(start = combined_joined$chromStart, end = combined_joined$chromEnd),
  strand = "*",
  name = combined_joined$name
)

# Overlap the abc promoters with all elements
overlaps <- findOverlaps(abc_canonical_tss_extended_gr, combined_joined_elements_gr)

# Create a mapping from elements to TSS genes
element_to_tss_map <- data.frame(
  element_idx = subjectHits(overlaps),
  tss_idx = queryHits(overlaps)
)

# Get the gene names for each TSS
element_to_tss_map$tss_gene <- abc_canonical_tss_gr$name[element_to_tss_map$tss_idx]

# Group by element and create a list of overlapping TSS genes
element_tss_list <- element_to_tss_map %>%
  group_by(element_idx) %>%
  summarize(overlapping_tss = list(tss_gene))

# Create a named vector for easy lookup
element_names <- combined_joined$name
names(element_names) <- seq_along(element_names)

# Add the column to combined_joined
combined_joined <- combined_joined %>%
  mutate(
    abc_tss_overlap = NA,
    element_idx = row_number()
  ) %>%
  left_join(element_tss_list, by = c("element_idx" = "element_idx")) %>%
  mutate(
    abc_tss_overlap = case_when(
      !is.na(overlapping_tss) ~ map_chr(overlapping_tss, ~ paste(.x, collapse = ",")),
      is.na(overlapping_tss) ~ NA_character_
    )
  ) %>%
  select(-element_idx, -overlapping_tss)


### ADD DISTANCE TO ABC TSS COLUMN ============================================

# Load in all genes (including lincRNAs etc)
abc_canonical_tss_all_genes <- read_tsv(snakemake@input$abc_canonical_tss)

# Check which measured genes are not in the ABC TSS list
measured_genes <- unique(combined_joined$measuredGeneSymbol)
abc_genes <- unique(abc_canonical_tss_all_genes$name)
missing_genes <- setdiff(measured_genes, abc_genes)
print(missing_genes)

# Do these genes have other names that are in our abc file? - looked at GeneCards for any aliases
# "TENT5B": "FAM46B" %in% abc_genes
# "KHDC4": "KIAA0907" %in% abc_genes
# "ZNRD2": "SSSCA1" %in% abc_genes
# "ANTKMT": "FAM173A" %in% abc_genes
# "JPT2": "HN1L" %in% abc_genes
# "ATP5MC1": "ATP5G1" %in% abc_genes
# "FDX2": "FDX1L" %in% abc_genes
# "GATD3A": "C21orf33" %in% abc_genes
# "RBIS": "C8orf59" %in% abc_genes

# Define alternative names mapping
alternative_names <- c(
  "TENT5B" = "FAM46B",
  "KHDC4" = "KIAA0907", 
  "ZNRD2" = "SSSCA1",
  "ANTKMT" = "FAM173A", 
  "JPT2" = "HN1L",
  "ATP5MC1" = "ATP5G1", 
  "FDX2" = "FDX1L",
  "GATD3A" = "C21orf33", 
  "RBIS" = "C8orf59"
)

# Create a lookup table for gene TSS centers
tss_centers <- abc_canonical_tss_all_genes %>%
  mutate(tss_center = (start + end) / 2) %>%
  select(`#chr`, name, tss_center)

# Calculate element centers and join with TSS data
combined_joined <- combined_joined %>%
  mutate(
    element_center = (chromStart + chromEnd) / 2,
    # Use alternative name for joining if it exists, otherwise use original name
    join_gene_name = ifelse(
      measuredGeneSymbol %in% names(alternative_names),
      alternative_names[measuredGeneSymbol],
      measuredGeneSymbol
    )
  ) %>%
  left_join(
    tss_centers %>% select(`#chr`, name, tss_center), 
    by = c("chrom" = "#chr", "join_gene_name" = "name")
  ) %>%
  mutate(
    # Calculate distance only when chromosomes match and TSS exists
    distance_to_abc_canonical_TSS = if_else(
      !is.na(tss_center), 
      abs(element_center - tss_center),
      NA_real_
    )
  ) %>%
  select(-element_center, -tss_center, -join_gene_name)

# Quick summary
message("Elements with TSS distance: ", sum(!is.na(combined_joined$distance_to_abc_canonical_TSS)) / nrow(combined_joined) * 100, "%")

### ADD COLUMN FOR OVERLAP WITH GENE BODY =====================================

# Create a column for a general flag indicating overlap with a protein coding gene body
combined_joined <- combined_joined %>%
  rowwise() %>%  
  mutate(
    # First create the list column as before
    protein_coding_gene_body_overlap = case_when(
      is.na(ValidConnection) ~ list(NA_character_),
      str_detect(ValidConnection, "overlaps potential protein coding gene body:") ~ {
        # Extract the substring, split on "|" and trim each element
        genes <- str_extract(ValidConnection, "(?<=overlaps potential protein coding gene body:)[^;]+")
        list(trimws(unlist(str_split(genes, "\\|"))))
      },
      TRUE ~ list(NA_character_)
    ),
    # Then create a string version using ALL elements in the list
    protein_coding_gene_body_overlap = case_when(
      all(is.na(unlist(protein_coding_gene_body_overlap))) ~ NA_character_,
      TRUE ~ paste(unlist(protein_coding_gene_body_overlap), collapse = ",")
    )
  ) %>%
  ungroup()


### ADD COLUMN FOR PROMOTER OVERLAP ===========================================

# Create the gencode_promoter_overlap column
combined_joined <- combined_joined %>%
  rowwise() %>%  
  mutate(
    # First create the list column as before
    gencode_promoter_overlap = case_when(
      is.na(ValidConnection) ~ list(NA_character_),
      str_detect(ValidConnection, "overlaps potential promoter:") ~ {
        # Extract the substring, split on "|" and trim each element
        genes <- str_extract(ValidConnection, "(?<=overlaps potential promoter:)[^;]+")
        list(trimws(unlist(str_split(genes, "\\|"))))
      },
      TRUE ~ list(NA_character_)
    ),
    # Then create a string version using ALL elements in the list
    gencode_promoter_overlap = case_when(
      all(is.na(unlist(gencode_promoter_overlap))) ~ NA_character_,
      TRUE ~ paste(unlist(gencode_promoter_overlap), collapse = ",")
    )
  ) %>%
  ungroup() 

# Add three columns indicating if the measuredGeneSymbol is in "TSS_control_gene", "abc_tss_overlap" or "gencode_promoter_overlap" T/F
combined_joined <- combined_joined %>%
  # First convert empty strings to NA
  mutate(
    abc_tss_overlap = if_else(abc_tss_overlap == "", NA_character_, abc_tss_overlap)
  ) %>%
  mutate(
    # First, convert comma-separated strings back to lists where needed
    abc_tss_overlap_list = case_when(
      is.na(abc_tss_overlap) ~ list(NA_character_),
      TRUE ~ lapply(strsplit(abc_tss_overlap, ","), function(x) x)
    ),
    gencode_promoter_overlap_list = case_when(
      is.na(gencode_promoter_overlap) ~ list(NA_character_),
      TRUE ~ lapply(strsplit(gencode_promoter_overlap, ","), function(x) x)
    )
  ) %>%
  rowwise() %>%
  mutate(
    # Now check if measuredGeneSymbol is in each list
    measuredGene_in_TSS_control = !is.na(TSS_control_gene) && measuredGeneSymbol == TSS_control_gene,
    
    # For abc_tss_overlap - check if measuredGeneSymbol is in the list as an exact match
    measuredGene_in_abc_tss = if(all(is.na(unlist(abc_tss_overlap_list)))) {
      FALSE
    } else {
      measuredGeneSymbol %in% unlist(abc_tss_overlap_list)
    },
    
    # For gencode_promoter_overlap - check if measuredGeneSymbol is in the list as an exact match
    measuredGene_in_gencode_promoter = if(all(is.na(unlist(gencode_promoter_overlap_list)))) {
      FALSE
    } else {
      measuredGeneSymbol %in% unlist(gencode_promoter_overlap_list)
    }
  ) %>%
  ungroup() %>%
  # Optionally remove the temporary list columns if not needed further
  select(-abc_tss_overlap_list, -gencode_promoter_overlap_list) %>%
  # Create a T/F column if the abc or gencode promoter overlap lists are TRUE
  mutate(abc_or_gencode_overlap = measuredGene_in_abc_tss | measuredGene_in_gencode_promoter)


### LABELLING EACH PAIR BY CATEGORY ===========================================

# E-G table: We want to label each pair as T/F with regard to each of the following categories
  # Category: `DistalElement-Gene`
  # Category: `DistalPromoter-Gene`
  # Category: `SelfPromoter`
  # Category: `Positive Control DistalElement-Gene`
  # Category: `Random DistalElement-Gene`
  # Category: `Random Validation DistalElement-Gene`

# Define vector of ENCODE filter flags
ValidConnection_Flags <- c("distance <1000", "overlaps target gene exon", "overlaps target gene intron") # Because the other flags are already processed and can be deduced by separate columns

# Edit the ValidConnection column to be "" when NA
combined_joined <- combined_joined %>%
  mutate(ValidConnection = ifelse(is.na(ValidConnection), "", ValidConnection))

# Create a new column indicating if the element overlaps a promoter defined by either abc_tss_overlap, gencode_promoter_overlap, or TSS_control_gene
combined_joined <- combined_joined %>%
  mutate(distal_or_promoter = ifelse(is.na(abc_tss_overlap) & is.na(gencode_promoter_overlap) & is.na(TSS_control_gene), "distal", "promoter"))

# DistalElement-Gene
  # filter out if element-gene "distance <1000", 
  # filter out if "abc_tss_overlap" "gencode_promoter_overlap" or "TSS_control_gene" are non-empty (meaning the element overlaps a promoter in abc or gencode OR is labelled as a TSS control gene) 
  # filter out if element "overlaps target gene intron", "overlaps target gene exon"

# Category: `DistalElement-Gene`
# These are all pairs that pass ENCODE filters - non-significant pairs that pass ENCODE filters, but are underpowered are included in this set
combined_joined_w_categories <- combined_joined %>%
  mutate(DistalElement_Gene = 
           ifelse(
             !str_detect(ValidConnection, str_c(ValidConnection_Flags, collapse = "|")) & (distal_or_promoter == "distal"), TRUE, FALSE
           )
         )                 

# IF "abc_tss_overlap" "gencode_promoter_overlap" or "TSS_control_gene" are non-empty (meaning the element overlaps with a promoter in abc or gencode OR is labelled as a TSS control gene)
  # DistalPromoter-Gene
    # filter out if "measuredGene_in_abc_tss", "measuredGene_in_gencode_promoter", OR "measuredGene_in_TSS_control"
  # selfPromoter
    # use if "measuredGene_in_abc_tss", "measuredGene_in_gencode_promoter", OR "measuredGene_in_TSS_control"

# Category: `DistalPromoter-Gene`
# Elements that overlap promoters (in ABC, GENCODE, or TSS controls) but not the measured gene's promoter
combined_joined_w_categories <- combined_joined_w_categories %>%
  mutate(DistalPromoter_Gene = 
           ifelse(
             distal_or_promoter == "promoter" & !measuredGene_in_abc_tss & !measuredGene_in_gencode_promoter & !measuredGene_in_TSS_control, TRUE, FALSE
           )
        )

# Category: `selfPromoter`
# Elements that overlap the measured gene's own promoter
combined_joined_w_categories <- combined_joined_w_categories %>%
  mutate(selfPromoter = 
           ifelse(
             measuredGene_in_abc_tss | measuredGene_in_gencode_promoter | measuredGene_in_TSS_control, 
             TRUE, FALSE
           )
        )

# Positive Control DistalElement-Gene
  # If the target_type == "DE" AND the intended positive control gene == the measuredGeneSymbol

# Category: `Positive Control DistalElement-Gene`
# Distal elements designed as positive controls paired with their intended target genes
combined_joined_w_categories <- combined_joined_w_categories %>%
  mutate(Positive_Control_DistalElement_Gene = 
           ifelse(
             target_type == "DE" & measuredGeneSymbol == de_assigned_gene & DistalElement_Gene, 
             TRUE, FALSE
           )
  )

# Random DistalElement-Gene
  # If the target_type == "enh" AND "DistalElement-Gene" == TRUE

# Category: `Random DistalElement-Gene`
# Random enhancers that are distal elements (not promoters)
combined_joined_w_categories <- combined_joined_w_categories %>%
  mutate(Random_DistalElement_Gene = 
           ifelse(
             target_type == "enh" & DistalElement_Gene == TRUE, 
             TRUE, FALSE
           )
  )

# Random Validation DistalElement-Gene
  # If Random DistalElement-Gene == TRUE AND element-gene pair is NOT in the training dataset

# Category: `Random Validation DistalElement-Gene`
# Random distal elements that are not in the training dataset
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
    mutate(Underpowered = ifelse(str_detect(ValidConnection, "< 80% power at 15% effect size"), TRUE, FALSE)) %>%
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
    sample_term_id = ifelse(ExperimentCellType == "K562", "EFO:0002067", "EFO:0009747"),
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
  select(-sample_term_id, -sample_summary_short, -notes) %>%
  left_join(
    combined_joined_w_categories_fixed %>%
      select(chrom, chromStart, chromEnd, measuredGeneSymbol, ExperimentCellType, 
             hg19_target_coords, EffectSize, chrTSS, startTSS, endTSS, distance_to_abc_canonical_TSS, hg19_target_chr, hg19_target_start, 
             hg19_target_end, target_name, distal_or_promoter, DistalElement_Gene, selfPromoter, DistalPromoter_Gene, TSS_control_gene,
             Positive_Control_DistalElement_Gene, Random_DistalElement_Gene, Random_Validation_DistalElement_Gene, de_assigned_gene,
             protein_coding_gene_body_overlap, gencode_promoter_overlap, abc_tss_overlap),
    by = c("targeting_chr" = "chrom", "targeting_start" = "chromStart", "targeting_end" = "chromEnd", 
           "gene_symbol" = "measuredGeneSymbol", "sample_term_name" = "ExperimentCellType")
  ) %>%
  mutate(pctEffectSize = EffectSize * 100) %>%
  dplyr::rename(
    intended_target_name_hg38 = intended_target_name,
    targeting_chr_hg38 = targeting_chr,
    targeting_start_hg38 = targeting_start,
    targeting_end_hg38 = targeting_end,
    intended_target_name_hg19 = hg19_target_coords,
    log_2_FC_effect_size = sceptre_log2_fc,
    pct_change_effect_size = pctEffectSize,
    chrTSS_hg38 = chrTSS,
    startTSS_hg38 = startTSS,
    endTSS_hg38 = endTSS,
    targeting_chr_hg19 = hg19_target_chr,
    targeting_start_hg19 = hg19_target_start,
    targeting_end_hg19 = hg19_target_end,
    design_file_target_name = target_name,
    distance_to_gencode_gene_TSS = distToTSS,
    element_location = distal_or_promoter,
    intended_positive_control_distal_element_target_gene = de_assigned_gene,
    intended_positive_control_target_gene = TSS_control_gene,
    cell_type = sample_term_name,
    gencode_protein_coding_gene_body_overlap = protein_coding_gene_body_overlap,
    design_file_type = type
  ) %>%
  mutate(
    element_gene_pair_identifier_hg38 = paste0(gene_symbol, "|", intended_target_name_hg38),
    element_gene_pair_identifier_hg19 = paste0(gene_symbol, "|", intended_target_name_hg19)
  ) %>%
  select(
    # HG38 target coordinates
    intended_target_name_hg38, 
    targeting_chr_hg38, 
    targeting_start_hg38, 
    targeting_end_hg38,
    
    # HG19 target coordinates
    intended_target_name_hg19, 
    targeting_chr_hg19, 
    targeting_start_hg19, 
    targeting_end_hg19,
    
    # Target names
    element_gene_pair_identifier_hg38,
    element_gene_pair_identifier_hg19,
    
    # Gene and cell type information
    gene_symbol, 
    gene_id, 
    cell_type,
    
    # Effect size and statistical measures
    log_2_FC_effect_size, 
    pct_change_effect_size, 
    sceptre_p_value, 
    sceptre_adj_p_value, 
    significant,
    
    # Distance to TSS
    distance_to_gencode_gene_TSS, 
    distance_to_abc_canonical_TSS,
    
    # TSS coordinates HG38
    chrTSS_hg38, 
    startTSS_hg38, 
    endTSS_hg38,
    
    # Element feature overlap
    element_location,
    gencode_promoter_overlap, 
    abc_tss_overlap,
    gencode_protein_coding_gene_body_overlap,
    
    # Element-Gene category
    DistalElement_Gene, 
    DistalPromoter_Gene, 
    selfPromoter, 
    Positive_Control_DistalElement_Gene, 
    Random_DistalElement_Gene, 
    Random_Validation_DistalElement_Gene,
    
    # Power Simulation Results
    power_at_effect_size_10, 
    power_at_effect_size_15, 
    power_at_effect_size_20, 
    power_at_effect_size_25, 
    power_at_effect_size_50,
    
    # Screen Design Information
    intended_positive_control_distal_element_target_gene, 
    intended_positive_control_target_gene,
    design_file_target_name, 
    design_file_type, 
    `guide_id(s)`
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
