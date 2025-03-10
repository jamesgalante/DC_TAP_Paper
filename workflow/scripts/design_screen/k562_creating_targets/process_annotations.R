# Script: process_annotations.R

### SETUP =====================================================================

# Saving image for debugging
if (!file.exists("RDA_objects/k562_creating_targets")) { dir.create("RDA_objects/k562_creating_targets", recursive = TRUE) }
save.image(paste0("RDA_objects/k562_creating_targets/process_annotations.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")

### LOADING FILES =============================================================

# Required packages and functions
suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
})

# Load input files
message("Loading input files...")
annotation_file <- import(snakemake@input$annotation_file)


### PROCESS ANNOTATIONS ======================================================
message("Processing annotations...")

# Add gene_base_id by removing version number from Ensembl IDs
annotations$gene_base_id <- sub("\\..*", "", annotations$gene_id)

# Filter for protein-coding genes on main chromosomes
filtered_annotations <- annotations[
  annotations$type == "gene" &
    seqnames(annotations) %in% paste0("chr", c(1:22, "X")) &
    annotations$gene_type == "protein_coding"
]

# Add TSS position based on gene strand
filtered_annotations$tss <- ifelse(
  as.character(strand(filtered_annotations)) == "+",
  start(filtered_annotations),
  end(filtered_annotations)
)

# Create a GRanges object for the TSS
filtered_annotations$tss_range <- GRanges(
  seqnames = seqnames(filtered_annotations),
  ranges = IRanges(
    start = filtered_annotations$tss,
    end = filtered_annotations$tss + 1
  ),
  strand = strand(filtered_annotations)
)

message(paste("Number of protein-coding genes retained:", length(filtered_annotations)))

### SAVE OUTPUT ==============================================================
message("Saving output files")
saveRDS(filtered_annotations, file = snakemake@output$processed_annotation)

### CLEAN UP =================================================================
message("Closing log file")
sink()
sink(type = "message")
close(log)