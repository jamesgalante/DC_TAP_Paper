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
})

message("Loading input files")
combined_validation <- read_tsv(snakemake@input$combined_validation)
combined_training <- read_tsv(snakemake@input$combined_training)

# Load in the guide_targets files
k562_guide_targets <- read_tsv(snakemake@input$guide_targets[[1]])
wtc11_guide_targets <- read_tsv(snakemake@input$guide_targets[[2]])


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

# There are a few distal element positive controls. We want to figure out what those controls intended targets were
positive_control_distal_elements_k562 <- combined_joined %>% 
  filter(target_type == "DE", category == "K562 DC TAP Seq") # Get all K562 distal elements because these aren't named by what gene they target
# control_genes <- k562_gene_table %>% filter(type == "control") %>% pull(gene) %>% unique()
control_genes <- c("TMEM98", "RASSF7", "RAE1", "PLP2", "PIM2", "GATA1", "UBR5", "LYL1", 
                   "DNASE2", "PPIF", "CEP104", "CTSD", "NFE2", "PHF20L1", "LRRC47", 
                   "RRP8", "TIMM10B", "RBM38", "LRRCC1", "CD164", "MYC", "PIM1", 
                   "FAM83A", "FADS1", "ALAS2", "EPB41", "TMEM41B", "JUNB", "UROS", 
                   "SVIP", "HBD", "CCDC26", "SMIM1")
t <- positive_control_distal_elements_k562 %>% filter(measuredGeneSymbol %in% control_genes)
# Can see how many of these target gene pairs are unique? - as in which targets are only paired with one gene
t %>% pull(target_name) %>% table() # Most have more than one gene they're paired with

# There are 38 unique K562 targets that are DEs
# 2 dropped out at power analysis because Diffex was NA (they didn't pass sceptre QC) `setdiff(all_de_controls, combined_joined %>% filter(target_type == "DE") %>% pull(target_name) %>% unique())`
# So now for each of these 40 targets, i have to figure out what the original intended gene symbol was to be the control

# 17 DEs are only paired with on of the "control" genes
easy_des <- names((t %>% pull(target_name) %>% table())[t %>% pull(target_name) %>% table() == 1])
# So for these 17 DEs, that control gene must be the gene that it was originally designed against

# For the others, I guess i can go through each of these groups of targets
# chr1:3691430-3691731
# This guy is paired with three TSSs: SMIM, LRRC47, CEP104
# However it's closest to SMIM

# chr11:5297067-5297368    chr11:5301767-5302068    chr11:5305872-5306173    chr11:5309369-5309670
# These guys are all paired with 3 TSSs: HBD, TIMM10B, RRP8
# However they're all closer to HBD

# chr19:12895376-12895677  chr19:12900763-12901064  
# These two guys are both paired with: JUNB, DNASE, LYL1
# However they're closest to JUNB

# chr19:13215352-13215653
# This guys is paired with: JUNB, DNASE, LYL1
# However it's closest to LYL1

# chr20:55990343-55990644
# This guy is paired with RAE1, RBM38
# However it's closest to RBM38

# chr8:128911048-128911349 chr8:128972548-128972849 
# These two guys are paired with MYC, CCDC26
# Seems like these are MYC though

# chr8:130594299-130594600 chr8:130701786-130702087 chr8:130705120-130705421
# These guys are paired with MYC, CCDC26
# Seems like these are CCDC26 though

# chrX:48641339-48641640   chrX:48658921-48659222   
# These guys are paired with GATA1, PIM2, PLP2
# Closest to GATA1

# chrX:48797830-48798131   chrX:48799012-48799313   
# These guys are paired with GATA1, PIM2, PLP2
# Closest to PIM2

# chrX:49005070-49005371   chrX:49022767-49023068   chrX:49023638-49023939
# These guys are paired with GATA1, PIM2, PLP2
# Closest to PLP2

# These pairs make the most logical sense, but I would have to confirm with the screen design - for now though, we can use this
# Let's make a map from all DE targets to the genes that they're supposed to be positive controls against

# We already have the easy des, let's manually do the others
easy_map_df <- t %>% 
  filter(target_name %in% easy_des) %>% 
  select(target_name, measuredGeneSymbol) %>%
  dplyr::rename(de_assigned_gene = measuredGeneSymbol)

# Now add the manual mappings based on your analysis
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

# Combine the easy mappings with the manual mappings
manual_mappings_df <- tibble(
  target_name = names(manual_mappings),
  de_assigned_gene = unlist(manual_mappings)
)

# Combine both dataframes
de_gene_map_df <- bind_rows(easy_map_df, manual_mappings_df)

# Let's now just add the intended positive control gene for each distal element target
combined_joined <- combined_joined %>%
  left_join(de_gene_map_df, by = "target_name") %>%
  mutate(de_assigned_gene = case_when(
    target_type == "DE" ~ de_assigned_gene,
    TRUE ~ NA
  )) %>%
  # Now add the WTC11 DE targets
  mutate(de_assigned_gene = case_when(
    category == "WTC11 DC TAP Seq" & target_type == "DE" ~ sapply(strsplit(target_name, "_"), function(x) 
      paste(unique(x[x != "DE" & !grepl("^\\d+$", x)]), 
            collapse = "_")),
    TRUE ~ de_assigned_gene
  ))


### CREATING THE SUMMARY STATS TABLE

# E-G table:  Include the following columns:
# Category: `DistalElement-Gene`
  # Filtered for ENCODE stuff 
  # Remove "TSS targeting guide(s)", "distance <1000", "overlaps target gene exon", "overlaps potential promoter", "overlaps target gene intron"
  # Don't have to filter out "< 80% power at 15% effect size" - because we want stats on this
# Category: `DistalPromoter-Gene`
  # All promoter-overlapping or TSS targets when paired with a gene not self
# Category: `SelfPromoter`
  # All promoter-overlapping or TSS targets when paired with themselves
# Category: `Positive Control DistalElement-Gene`
  # All "DE"s when paired with the gene that they're supposed to be tested against
# Category: `Random DistalElement-Gene`
  # From all `DistalElement-Gene` pairs, remove the non-random pairs (For the random set, we can define that as every pair where targets are from the random 2Mb loci (except that one instance of NANOG DE pos control being within a 2Mb region in wtc11))
# I don't think this should exist
# Category: `Random Validation DistalElement-Gene`
  # Pairs that should be in the validation dataset


# Escape parentheses for regex matching
ValidConnection_Flags <- c(
  "TSS targeting guide\\(s\\)", 
  "distance <1000", 
  "overlaps target gene exon", 
  "overlaps potential promoter", 
  "overlaps target gene intron"
)

# Category: `DistalElement-Gene`
combined_joined_w_categories <- combined_joined %>%
  mutate(DistalElement_Gene = case_when(
    Significant == TRUE & !str_detect(ValidConnection, str_c(ValidConnection_Flags, collapse = "|")) ~ TRUE,
    Significant == FALSE & !str_detect(ValidConnection, str_c(ValidConnection_Flags, collapse = "|")) ~ TRUE,
    TRUE ~ FALSE
  ))                 

# There are two cases that we want to account for
  # When there's "overlaps potential promoter"
    # what promoter is it overlapping -> is that promoter the same as the measuredGeneSymbol
  # "TSS targeting guide(s)"
    # In this case is the measuredGeneSymbol == TSS_control_gene

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
# We prioritize the results returned by "overlaps potential promoter" for deciding if it's selfPromoter or DistalPromoter_Gene
combined_joined_w_categories <- combined_joined_w_categories %>%
  # Create a helper column that picks the promoter type
  mutate(
    promoter_type = case_when(
      str_detect(ValidConnection, "overlaps potential promoter:") ~ "overlap",
      str_detect(ValidConnection, "TSS targeting") ~ "tss",
      TRUE ~ "none"
    ),
    # When both conditions are present, you could decide to use "overlap"
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
  )

# Category: `Positive Control DistalElement-Gene`
combined_joined_w_categories <- combined_joined_w_categories %>%
  mutate(Positive_Control_DistalElement_Gene = case_when(
    target_type == "DE" & de_assigned_gene == measuredGeneSymbol ~ TRUE,
    TRUE ~ FALSE
  ))

# Category: `Random DistalElement-Gene`
# Need the locus for each target
# For the random set, we can define that as every pair where targets are from the random 2Mb loci (except that one instance of NANOG DE pos control being within a 2Mb region in wtc11)
# This should just be all "enh" and all "tss_random" right? - yes
combined_joined_w_categories <- combined_joined_w_categories %>%
  mutate(Random_DistalElement_Gene = case_when(
    target_type %in% c("tss_random", "enh") ~ TRUE,
    TRUE ~ FALSE
  ))

# Category: `Random Validation DistalElement-Gene`
# Random element pairs that are also DistalElement pairs
combined_joined_w_categories <- combined_joined_w_categories %>%
  mutate(Random_Validation_DistalElement_Gene = ifelse(Random_DistalElement_Gene == TRUE & DistalElement_Gene == TRUE, TRUE, FALSE))



### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
# saveRDS(combined_joined, snakemake@output$unfiltered_validation_with_guide_targets)

# Save a text file with all of these summary statistics
# Save a bunch of plots which display all of these numbers


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)
