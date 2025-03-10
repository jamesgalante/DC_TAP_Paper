# Script: generate_outputs.R
### SETUP =====================================================================
# Saving image for debugging
if (!file.exists("RDA_objects/k562_creating_targets")) { dir.create("RDA_objects/k562_creating_targets", recursive = TRUE) }
save.image(paste0("RDA_objects/k562_creating_targets/generate_outputs.rda"))
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
  library(GenomicRanges)
})

# Load input files
message("Loading input files...")
sampled_loci <- readRDS(snakemake@input$sampled_loci)
processed_annotation <- readRDS(snakemake@input$processed_annotation)

### CREATE LOCI GENE INFO ====================================================
message("Creating loci summary file...")
# Generate loci summary with selected genes
loci_gene_info <- data.frame(
  locus_central_gene = sampled_loci$locus,
  locus_chr = sampled_loci$chr,
  locus_start = sampled_loci$locus_start,
  locus_end = sampled_loci$locus_end,
  total_dhs = sampled_loci$dhs,
  total_genes = sampled_loci$genes,
  genes_above_tpm = sampled_loci$genes_above_tpms,
  gene_ids = sapply(sampled_loci$gene_ids, function(x) paste(x, collapse = ";")),
  gene_names = sapply(sampled_loci$gene_names, function(x) paste(x, collapse = ";")),
  genes_above_tpm_ids = sapply(sampled_loci$genes_above_tpm_ids, function(x) paste(x, collapse = ";"))
)

### CREATE CANDIDATE PEAKS BED ===============================================
message("Creating candidate peaks file...")
# Create a BED file with all candidate peaks for the screen
candidate_peaks <- data.frame()

for (i in 1:nrow(sampled_loci)) {
  # Get list of peaks for this locus
  peak_ids <- unlist(strsplit(sampled_loci$peak_sample[i], ";"))
  
  if (length(peak_ids) > 0) {
    for (peak_id in peak_ids) {
      # Parse peak coordinates
      peak_parts <- strsplit(peak_id, ":")[[1]]
      chr <- peak_parts[1]
      
      coord_parts <- strsplit(peak_parts[2], "_")[[1]]
      start <- as.numeric(coord_parts[1])
      end <- as.numeric(coord_parts[2])
      
      # Get corresponding read count and quantile
      locus_dhs_ids <- sampled_loci$locus_dhs_ids[[i]]
      peak_index <- match(peak_id, locus_dhs_ids)
      
      if (!is.na(peak_index)) {
        read_count <- sampled_loci$dhs_reads[[i]][peak_index]
        quantile <- sampled_loci$dhs_quants[[i]][peak_index]
        
        # Get TSS distance
        tss_distance <- sampled_loci$dhs_nearest_tss_distance[[i]][peak_index]
        nearest_gene <- sampled_loci$dhs_nearest_tss[[i]][peak_index]
        
        # Add to candidate peaks data frame
        peak_row <- data.frame(
          chr = chr,
          start = start,
          end = end,
          peak_id = peak_id,
          locus = sampled_loci$locus[i],
          read_count = read_count,
          quantile = quantile,
          nearest_gene = nearest_gene,
          tss_distance = tss_distance
        )
        
        candidate_peaks <- rbind(candidate_peaks, peak_row)
      }
    }
  }
}

### CREATE SELECTED LOCI BED =================================================
message("Creating selected loci BED file...")
# Create a BED file with the selected loci
selected_loci <- data.frame(
  chr = sampled_loci$chr,
  start = sampled_loci$locus_start,
  end = sampled_loci$locus_end,
  name = sampled_loci$locus,
  score = 0,
  strand = ".",
  total_dhs = sampled_loci$dhs,
  sampled_dhs = sapply(strsplit(sampled_loci$peak_sample, ";"), length),
  genes_above_tpm = sampled_loci$genes_above_tpms
)

# Ensure all data frames have the right column types
loci_gene_info <- as.data.frame(loci_gene_info, stringsAsFactors = FALSE)
candidate_peaks <- as.data.frame(candidate_peaks, stringsAsFactors = FALSE)
selected_loci <- as.data.frame(selected_loci, stringsAsFactors = FALSE)

### SAVE OUTPUT ==============================================================
message("Saving output files")
write_tsv(loci_gene_info, file = snakemake@output$loci_gene_info)
write_tsv(candidate_peaks, file = snakemake@output$candidate_peaks)
write_tsv(selected_loci, file = snakemake@output$selected_loci)

### CLEAN UP =================================================================
message("Closing log file")
sink()
sink(type = "message")
close(log)