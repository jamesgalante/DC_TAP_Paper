# Script: sample_loci.R
### SETUP =====================================================================
# Saving image for debugging
if (!file.exists("RDA_objects/k562_creating_targets")) { dir.create("RDA_objects/k562_creating_targets", recursive = TRUE) }
save.image(paste0("RDA_objects/k562_creating_targets/sample_loci.rda"))
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
  library(cowplot)
})

# Load input files
message("Loading input files...")
usable_loci <- readRDS(snakemake@input$filtered_loci)

# Get parameters
num_loci <- snakemake@params$num_loci
seed_value <- snakemake@params$seed

### PREPARE FOR SAMPLING =====================================================
message("Preparing loci for sampling...")
message(paste("Will sample", num_loci, "loci using seed", seed_value))

# Add chromosomal information to loci if not present
if (!"chr" %in% colnames(usable_loci)) {
  usable_loci$chr <- NA
  for (i in 1:nrow(usable_loci)) {
    # Extract chromosome from the first DHS ID
    if (length(usable_loci$locus_dhs_ids[[i]]) > 0) {
      chr <- strsplit(usable_loci$locus_dhs_ids[[i]][1], ":")[[1]][1]
      usable_loci$chr[i] <- chr
    }
  }
}

### SAMPLE LOCI =============================================================
# Set seed for reproducibility
set.seed(seed_value)

message("Sampling loci...")
# Sample loci based on criteria
sampled_loci <- usable_loci %>%
  sample_n(size = num_loci)

### SAMPLE PEAKS WITHIN LOCI ================================================
message("Sampling peaks within each locus...")
# Sample DHS peaks from each locus
sampled_loci$peak_sample <- NA
sampled_loci$sampled_stats <- NA
sampled_loci$sampled_quants <- NA
sampled_loci$sampled_reads <- NA
sampled_loci$sampled_nearest_tss_distance <- NA

for (i in 1:nrow(sampled_loci)) {
  distances <- sampled_loci$dhs_nearest_tss_distance[[i]]
  quant_list <- sampled_loci$dhs_quants[[i]]
  peak_list <- sampled_loci$locus_dhs_ids[[i]]
  reads_list <- sampled_loci$dhs_reads[[i]]
  
  # Get peaks >100kb from TSS
  g100_kb_peak_indexs <- which(distances > 100000)
  g100_kb_sampled_indexs <- sample(g100_kb_peak_indexs, size = min(10, length(g100_kb_peak_indexs)))
  g100kb_sampled_peaks <- peak_list[g100_kb_sampled_indexs]
  
  # Count how many of the sampled distant peaks are in top 75% quantile
  q75_sampled_already <- sum(quant_list[g100_kb_sampled_indexs] == "75_1")
  
  # Get peaks <100kb but >100bp from TSS
  l100_kb_peak_indexs <- which(distances <= 100000 & distances > 100)
  
  if (length(l100_kb_peak_indexs) > 0) {
    l100kb_peaks <- peak_list[l100_kb_peak_indexs]
    quants_l100kb_peaks <- quant_list[l100_kb_peak_indexs]
    
    # Top 75% quantile peaks
    top75_l100kb <- l100kb_peaks[quants_l100kb_peaks == "75_1"]
    # Other quantile peaks
    bottom75_l100kb <- l100kb_peaks[quants_l100kb_peaks != "75_1"]
    
    # Sample from top 75% quantile peaks
    if (length(top75_l100kb) >= (20 - q75_sampled_already)) {
      top75_l100kb_sampled <- sample(top75_l100kb, size = 20 - q75_sampled_already)
    } else {
      top75_l100kb_sampled <- top75_l100kb
    }
    
    # Sample from other quantile peaks
    remaining_slots <- 40 - length(top75_l100kb_sampled) - length(g100kb_sampled_peaks)
    if (length(bottom75_l100kb) >= remaining_slots) {
      bottom75_l100kb_sampled <- sample(bottom75_l100kb, size = remaining_slots)
    } else {
      bottom75_l100kb_sampled <- bottom75_l100kb
    }
    
    # Combine all sampled peaks
    full_sample <- unique(c(g100kb_sampled_peaks, top75_l100kb_sampled, bottom75_l100kb_sampled))
  } else {
    # If no peaks in the desired distance range, just use what we have from >100kb
    full_sample <- g100kb_sampled_peaks
  }
  
  # Create statistics for the sampled peaks
  sampled_peak_indexs <- match(full_sample, peak_list)
  sampled_distances <- distances[sampled_peak_indexs]
  sampled_quants <- quant_list[sampled_peak_indexs]
  sampled_reads <- reads_list[sampled_peak_indexs]
  
  g100 <- sum(sampled_distances > 100000)
  l100 <- sum(sampled_distances <= 100000 & sampled_distances > 1000)
  q75 <- sum(sampled_quants == "75_1")
  b75 <- sum(sampled_quants != "75_1")
  
  # Store sampled peaks and their statistics
  sampled_loci$peak_sample[i] <- paste(full_sample, collapse = ";")
  sampled_loci$sampled_stats[i] <- paste(c(g100, l100, q75, b75), collapse = ";")
  sampled_loci$sampled_quants[i] <- paste(sampled_quants, collapse = ";")
  sampled_loci$sampled_reads[i] <- paste(sampled_reads, collapse = ";")
  sampled_loci$sampled_nearest_tss_distance[i] <- paste(sampled_distances, collapse = ";")
}

### CHECK FOR OVERLAPPING LOCI ===============================================
message("Checking for overlapping loci...")
# Create GRanges object for sampled loci
sampled_loci_gr <- GRanges(
  seqnames = sampled_loci$chr,
  ranges = IRanges(start = sampled_loci$locus_start, end = sampled_loci$locus_end),
  locus = sampled_loci$locus
)

# Check for overlapping loci
overlaps <- findOverlaps(sampled_loci_gr, sampled_loci_gr, ignore.self = TRUE)
if (length(overlaps) > 0) {
  message("Warning: Found overlapping loci in the sample!")
  print(overlaps)
} else {
  message("No overlapping loci found in the sample.")
}

### PARSE STATISTICS ========================================================
# Parse sampled stats
sampled_loci$g100_peaks <- sapply(strsplit(sampled_loci$sampled_stats, ";"), function(x) as.numeric(x[1]))
sampled_loci$l100_peaks <- sapply(strsplit(sampled_loci$sampled_stats, ";"), function(x) as.numeric(x[2]))
sampled_loci$q75_peaks <- sapply(strsplit(sampled_loci$sampled_stats, ";"), function(x) as.numeric(x[3]))
sampled_loci$b75_peaks <- sapply(strsplit(sampled_loci$sampled_stats, ";"), function(x) as.numeric(x[4]))

### VISUALIZATION ============================================================
message("Creating sampling visualization plots...")
# Create summary plots
g100_plot <- ggplot(sampled_loci, aes(x = g100_peaks)) +
  geom_histogram(bins = 10) +
  labs(title = "Distribution of peaks >100kb from TSS in sampled loci",
       x = "Number of peaks >100kb from TSS",
       y = "Count") +
  theme_bw()

q75_plot <- ggplot(sampled_loci, aes(x = q75_peaks)) +
  geom_histogram(bins = 10) +
  labs(title = "Distribution of peaks in top 75% quantile in sampled loci",
       x = "Number of peaks in top 75% quantile",
       y = "Count") +
  theme_bw()

scatter_plot <- ggplot(sampled_loci, aes(x = q75_peaks, y = g100_peaks)) +
  geom_point() +
  labs(title = "Distribution of sampled loci",
       x = "Peaks in top 75% quantile",
       y = "Peaks >100kb from TSS") +
  theme_bw()

# Plot chromosome distribution
chr_plot <- ggplot(sampled_loci, aes(x = chr)) +
  geom_bar() +
  labs(title = "Chromosome distribution of sampled loci",
       x = "Chromosome",
       y = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine plots
combined_plot <- plot_grid(
  g100_plot, q75_plot, scatter_plot, chr_plot,
  ncol = 2,
  labels = c("A", "B", "C", "D")
)

### SAVE OUTPUT ==============================================================
message("Saving output files")
saveRDS(sampled_loci, file = snakemake@output$sampled_loci)
ggsave(snakemake@output$sampling_plots, plot = combined_plot, width = 12, height = 10)

### CLEAN UP =================================================================
message("Closing log file")
sink()
sink(type = "message")
close(log)