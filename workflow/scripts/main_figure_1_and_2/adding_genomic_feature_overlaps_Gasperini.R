# Script: adding_genomic_feature_overlaps_Gasperini.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/adding_genomic_feature_overlaps_Gasperini.rda"))
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
combined_joined <- read_tsv(snakemake@input$results_with_design_file_features)
abc_canonical_tss <- read_tsv(snakemake@input$abc_canonical_tss)


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

# Create a lookup table for gene TSS centers
tss_centers <- abc_canonical_tss_all_genes %>%
  mutate(tss_center = (start + end) / 2) %>%
  select(`#chr`, name, tss_center)

# Calculate element centers and join with TSS data
combined_joined <- combined_joined %>%
  mutate(
    element_center = (chromStart + chromEnd) / 2
  ) %>%
  left_join(
    tss_centers %>% select(`#chr`, name, tss_center), 
    by = c("chrom" = "#chr", "measuredGeneSymbol" = "name")
  ) %>%
  mutate(
    # Calculate distance only when chromosomes match and TSS exists
    distance_to_abc_canonical_TSS = if_else(
      !is.na(tss_center), 
      abs(element_center - tss_center),
      NA_real_
    )
  ) %>%
  select(-element_center, -tss_center)

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


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(combined_joined, snakemake@output$results_with_design_file_and_genomic_feature_overlaps)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)