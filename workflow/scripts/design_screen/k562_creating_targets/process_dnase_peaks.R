# Script: process_dnase_peaks.R
### SETUP =====================================================================
# Saving image for debugging
if (!file.exists("RDA_objects/k562_creating_targets")) { dir.create("RDA_objects/k562_creating_targets", recursive = TRUE) }
save.image(paste0("RDA_objects/k562_creating_targets/process_dnase_peaks.rda"))
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
  library(ggplot2)
  library(GenomicRanges)
  library(cowplot)
})

# Load input files
message("Loading input files...")
dnase_peaks_file <- snakemake@input$dnase_peaks
message(paste("Loading DNase-seq peaks from:", dnase_peaks_file))

### IMPORT DNASE PEAKS =======================================================
# Column names in a narrowPeaks bed file
peak_colnames <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", 
                   "signalValue", "pValue", "qValue", "peak", "reads")

# Column types in a narrowPeaks bed file with added read counts
peak_cols <- cols(
  chrom = col_character(),
  chromStart = col_double(),
  chromEnd = col_double(),
  name = col_character(),
  score = col_double(),
  strand = col_character(),
  signalValue = col_double(),
  pValue = col_double(),
  qValue = col_double(),
  peak = col_double(),
  reads = col_double()
)

message("Reading DNase-seq peaks...")
# Import ENCODE DNase-seq peaks
dhs <- read_tsv(dnase_peaks_file, col_names = peak_colnames, col_types = peak_cols)

### PROCESS PEAKS ============================================================
message("Processing DNase-seq peaks...")
# Normalize reads for peak width by calculating RPKM
dhs <- dhs %>% 
  mutate(
    peak_width = chromEnd - chromStart,
    reads_rpkm = reads / (sum(reads) / 1e6) / (peak_width / 1000),
    id = paste0(chrom, ":", chromStart, "_", chromEnd)
  )

# Calculate quantiles for peak read counts
read_quantiles <- quantile(dhs$reads, probs = seq(0, 1, 0.25))
message("DNase-seq read count quantiles:")
print(read_quantiles)

# Assign quantile categories to each peak
dhs <- dhs %>%
  mutate(quant = case_when(
    reads < read_quantiles[2] ~ "0_25",
    reads >= read_quantiles[2] & reads < read_quantiles[3] ~ "25_50",
    reads >= read_quantiles[3] & reads < read_quantiles[4] ~ "50_75",
    reads >= read_quantiles[4] ~ "75_1"
  ))

### VISUALIZATION ============================================================
message("Creating DNase-seq peak distribution plots...")
reads_hist <- ggplot(dhs, aes(x = reads)) +
  geom_histogram(bins = 50, fill = "orange") +
  labs(title = paste(nrow(dhs), "K562 DHS"), 
       x = "DNase-seq reads in DHS") +
  theme_bw()

reads_ecdf <- ggplot(dhs, aes(x = reads)) + 
  stat_ecdf(geom = "step") +
  labs(title = "K562 DHS", 
       x = "DNase-seq reads in DHS", 
       y = "cumulative fraction") +
  theme_bw()

rpkm_ecdf <- ggplot(dhs, aes(x = reads_rpkm)) + 
  stat_ecdf(geom = "step") +
  labs(title = "RPKM normalized K562 DHS", 
       x = "RPKM normalized DNase-seq reads in DHS", 
       y = "cumulative fraction") +
  theme_bw()

# Combine plots
combined_plot <- plot_grid(
  reads_hist, reads_ecdf, rpkm_ecdf, 
  ncol = 1, 
  labels = c("A", "B", "C")
)

### CREATE GENOMIC RANGES ====================================================
# Create GRanges object for DNase peaks
dhs_gr <- makeGRangesFromDataFrame(dhs, keep.extra.columns = TRUE)

message(paste("Total number of DNase-seq peaks:", length(dhs_gr)))
message(paste("Peak quantile breakdown:", paste(table(dhs$quant), collapse = ", ")))

### SAVE OUTPUT ==============================================================
message("Saving output files")
saveRDS(dhs_gr, file = snakemake@output$processed_peaks)
cowplot::save_plot(snakemake@output$peak_distribution, combined_plot, 
                   base_width = 8, base_height = 12)

### CLEAN UP =================================================================
message("Closing log file")
sink()
sink(type = "message")
close(log)