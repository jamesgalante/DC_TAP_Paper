# Script: modify_specific_pairs_in_final_file.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/modify_specific_pairs_in_final_file.rda"))
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
summarized_categories <- read_tsv(snakemake@input$summarized_categories)


### MODIFYING =================================================================

# Call design_file_target_name = RTN4_TSS_8 an enhancer as it's currently labelled as a promoter, but is not an active TSS in IGV
# Change element_location to "distal"
# Change DistalElement_Gene to TRUE
  # Keep Random_DistalElement_Gene FALSE as this element was not chosen randomly and was initially intended to be a positive control

# "CCDC26|chr8:130594299-130594600" can be converted to Positive Control DistalElement Gene
# element_gene_pair_identifier_hg38 == "CCDC26|chr8:130594299-130594600"
# DistalElement_Gene = TRUE
# Positive_Control_DistalElement_Gene = TRUE

Formatted_DC_TAP_Seq_Results <- summarized_categories %>%
  mutate(
    # Modify RTN4_TSS_8 to be an enhancer
    element_location = ifelse(design_file_target_name == "RTN4_TSS_8", "distal", element_location),
    
    # Set DistalElement_Gene to TRUE for RTN4_TSS_8
    DistalElement_Gene = ifelse(design_file_target_name == "RTN4_TSS_8", TRUE, DistalElement_Gene),
    
    # Modify CCDC26 entry
    DistalElement_Gene = ifelse(element_gene_pair_identifier_hg38 == "CCDC26|chr8:130594299-130594600", TRUE, DistalElement_Gene),
    Positive_Control_DistalElement_Gene = ifelse(element_gene_pair_identifier_hg38 == "CCDC26|chr8:130594299-130594600", TRUE, Positive_Control_DistalElement_Gene)
  )


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(Formatted_DC_TAP_Seq_Results, snakemake@output$Formatted_DC_TAP_Seq_Results)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)