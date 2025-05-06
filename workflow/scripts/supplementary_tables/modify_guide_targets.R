# Script: modify_guide_targets.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/modify_guide_targets.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
suppressPackageStartupMessages(
  library(tidyverse)
)

message("Loading input files")
# Load in the guide_targets files
k562_guide_targets <- read_tsv(snakemake@input$guide_targets[[1]])
wtc11_guide_targets <- read_tsv(snakemake@input$guide_targets[[2]])

# Load in the design files
k562_design_file <- read_tsv(snakemake@input$k562_design_file)
wtc11_design_file <- read_tsv(snakemake@input$wtc11_design_file)


### K562 ======================================================================

# Replace the spacer column, which has the PAM sequence
k562_guide_targets_supp_table <- k562_guide_targets %>%
  mutate(spacer = k562_design_file$GuideSequenceMinusG) %>%
  dplyr::rename(GuideSequenceMinusG = spacer)


### WTC11 =====================================================================

# Replace the sgOpti name with CROP for compatibility with guide_targets file
wtc11_design_file$name <- gsub("WTC11_Random_Screen_sgOpti-DC", "WTC11_Random_Screen_Crop", wtc11_design_file$name)

# Add in 10 guides from the Jones group - obtained from previous TF-Perturb-Seq experiment
jones_group_guides <- data.frame(
  name = c(
    "WTC11_Random_Screen_Crop_SOX2-1553", "WTC11_Random_Screen_Crop_TCF4", "WTC11_Random_Screen_Crop_YBX1-821",
    "WTC11_Random_Screen_Crop_YBX1-823", "WTC11_Random_Screen_Crop_OCT4", "WTC11_Random_Screen_Crop_NANOG",
    "WTC11_Random_Screen_Crop_SOX2", "WTC11_Random_Screen_Crop_BAG3", "WTC11_Random_Screen_Crop_ROCK1",
    "WTC11_Random_Screen_Crop_GSK3B"
  ),
  GuideSequenceMinusG = c(
    "CGGCAGGATCGGCCAGAGG", "CATCACCATGGACTCCCCCG", "CCATAGAGTACCCGAGAGGA",
    "CCGCCGGCCGCCATAGAGAC", "GGTGAAATGAGGGCTTGCGA", "CCAGCAGAACGTTAAAATCC",
    "CCCTGACAGCCCCCGTCACA", "TTCATAAAGGTGCCCGGCGC", "CGGGGCGCGGACGCTCGGAA",
    "GGGTGGCTCGGAGATGCGAC"
  )
)

# Replace the spacer column, which has the PAM sequence
wtc11_guide_targets_supp_table <- wtc11_guide_targets %>%
  left_join(
    wtc11_design_file %>%
      select(name, GuideSequenceMinusG),
    by = "name"
  ) %>%
  mutate(spacer = GuideSequenceMinusG) %>%
  select(-GuideSequenceMinusG) %>%
  dplyr::rename(GuideSequenceMinusG = spacer) %>%
  # Update GuideSequenceMinusG for Jones group guides
  rows_update(jones_group_guides, by = "name") %>%
  # Move GuideSequenceMinusG after strand
  relocate(GuideSequenceMinusG, .after = strand)


# Note: There are 6 guide sequences that are duplicated, and thus have two named guides with the same sequence. 
# Because the duplicates have the same target, this doesn't change the outcome of the screen
# But, because the screen was run with each individual guide, these duplicates are not removed from the table
# The names with duplicated sequences:
  # WTC11_Random_Screen_Crop_14742 WTC11_Random_Screen_Crop_14743 
  # WTC11_Random_Screen_Crop_7289 WTC11_Random_Screen_Crop_ROCK1 
  # WTC11_Random_Screen_Crop_14744 WTC11_Random_Screen_Crop_14745 
  # WTC11_Random_Screen_Crop_14746 WTC11_Random_Screen_Crop_14747 
  # WTC11_Random_Screen_Crop_2422 WTC11_Random_Screen_Crop_BAG3 
  # WTC11_Random_Screen_Crop_14748 WTC11_Random_Screen_Crop_14749 
# Two of these duplicate pairs are because 10 guides were added from the Jones group


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(k562_guide_targets_supp_table, snakemake@output$k562_guide_targets_supp_table)
write_tsv(wtc11_guide_targets_supp_table, snakemake@output$wtc11_guide_targets_supp_table)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)