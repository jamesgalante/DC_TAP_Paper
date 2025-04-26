# Script: resize_dc_tap_elements.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/resize_dc_tap_elements.rda"))
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
  library(GenomicRanges)
})

message("Loading input files")
Formatted_DC_TAP_Seq_Results <- read_tsv(snakemake@input$Formatted_DC_TAP_Seq_Results)


### RESIZE DC TAP ELEMENTS ====================================================

# For overlapping DC TAP elements with different tracks, we want to resize the elements to 500 bp
message("Resizing DC TAP elements to 500bp")

# Function to resize genomic elements using GenomicRanges
resize_genomic_elements <- function(df, size = 500) {
  # Create GRanges object from the input dataframe
  gr <- GRanges(
    seqnames = df$targeting_chr_hg38,
    ranges = IRanges(
      start = df$targeting_start_hg38,
      end = df$targeting_end_hg38
    )
  )
  
  # Identify elements that need resizing (smaller than target size)
  resize_elements <- width(gr) < size
  
  # Resize elements that need it, keeping the center fixed
  gr[resize_elements] <- resize(gr[resize_elements], width = size, fix = "center")
  
  # Create new columns with resized coordinates
  df$resized_targeting_chr_hg38 <- df$targeting_chr_hg38
  df$resized_targeting_start_hg38 <- start(gr)
  df$resized_targeting_end_hg38 <- end(gr)
  df$resized_intended_target_name_hg38 <- paste0(
    df$targeting_chr_hg38, ":",
    df$resized_targeting_start_hg38, "-",
    df$resized_targeting_end_hg38
  )
  
  # Reorder columns to put resized columns at the beginning
  df <- df %>%
    relocate(
      resized_intended_target_name_hg38,
      resized_targeting_chr_hg38,
      resized_targeting_start_hg38,
      resized_targeting_end_hg38
    )
  
  return(df)
}

# Apply the resize function to the dataframe
resized_Formatted_DC_TAP_Seq_Results <- resize_genomic_elements(Formatted_DC_TAP_Seq_Results, size = 500)

message("Resizing complete")


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(resized_Formatted_DC_TAP_Seq_Results, snakemake@output$resized_Formatted_DC_TAP_Seq_Results)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)