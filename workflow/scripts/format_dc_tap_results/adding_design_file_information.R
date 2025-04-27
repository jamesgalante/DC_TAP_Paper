
# Script: adding_design_file_information.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/adding_design_file_information.rda"))
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
combined_validation <- read_tsv(snakemake@input$combined_validation)

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
results_with_design_file_features <- combined_joined %>%
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
  mutate(de_assigned_gene = ifelse(CellType == "WTC11" & measuredGeneSymbol == "POU5F1" & target_type == "DE", "POU5F1", de_assigned_gene))


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(results_with_design_file_features, snakemake@output$results_with_design_file_features)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)