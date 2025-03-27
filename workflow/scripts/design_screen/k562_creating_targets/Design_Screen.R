############################################################
# Exploration ENCODE CRISPRi Screen Design
#
# This script explores candidate loci for an ENCODE enhancer 
# screen. It performs the following major tasks:
#   - Loads TPM data (from Gasperini et al. 2019) to filter genes above a threshold (e.g. 50 TPM).
#   - Loads genome annotations (GENCODE v26 and v29) and processes them.
#   - Loads and normalizes DNase-seq (DHS) peak data.
#   - Computes locus statistics for candidate gene loci (e.g., number of DHS peaks,
#     genes, TPM values, distances from peaks to TSS).
#   - Randomly samples sets of loci to assess variability in locus composition.
#   - Produces various summary plots and writes output files.
#
# Inputs (file paths and parameters) are defined in the "Input Variables" section.
# Outputs (tables, BED files, plots) are written to the specified paths.
############################################################

###############################
# Input Variables and Settings
###############################

# --- FILE PATHS ---
tpm_filepath            <- "resources/design_screen/k562_creating_targets/tpm.csv"
gencode_v26_filepath    <- "results/design_screen/gencode.v26lift37.annotation.gtf.gz"
gencode_v29_filepath    <- "results/design_screen/gencode.v29.annotation.gtf.gz"
dhs_filepath            <- "results/design_screen/ENCFF185XRG_w_ENCFF325RTP_q30_sorted.txt"
rep6_loci_filepath      <- "resources/design_screen/k562_creating_targets/loci_only_bed_version_rep6.bed"

# --- OUTPUT PATHS (example outputs) ---
locus_stats_output      <- "results/design_screen/locus_stats_full_rep6.txt"
rep_summary_output      <- "rep_summary_tab.txt"

# --- PARAMETERS ---
tpm_threshold           <- 50      # Minimum TPM to consider a gene expressed.
locus_width             <- 2e6     # Locus window width (e.g., 2Mb).
n_loci_sample           <- 25      # Number of loci to randomly sample.
n_sample_reps           <- 100     # Number of repetitions for sampling.
seed_value              <- 211005  # Seed for random sampling reproducibility.
n_workers               <- 5       # Number of workers for parallel processing.

###############################
# Setup: Load Libraries
###############################

library(tidyverse)
library(here)
library(rtracklayer)
library(BiocParallel)
library(cowplot)

###############################
# Loading TPM Data
###############################

# Load the TPM data (assumes TPM values in a column named "tpm" and gene IDs in "gene").
tpm <- read_csv(tpm_filepath,
                col_types = cols(.default = col_character(), tpm = col_double()))
n_above_50 <- sum(tpm$tpm >= tpm_threshold)
cat("Number of genes above", tpm_threshold, "TPM:", n_above_50, "\n")

# Plot TPM distribution
p_tpm <- ggplot(tpm, aes(tpm)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = tpm_threshold, color = "red") +
  labs(x = "Transcripts-per-million (TPM)",
       title = paste("Genes above", tpm_threshold, "TPM:", n_above_50)) +
  scale_x_log10() +
  theme_bw()
print(p_tpm)

###############################
# Loading Genome Annotations
###############################

# Import GENCODE annotations for v26
v26 <- rtracklayer::import(gencode_v26_filepath, format = "gtf")

# Remove version numbers from gene IDs in v26.
v26$gene_base_id <- sub("\\..*", "", v26$gene_id)

###############################
# Loading and Processing DNase-seq (DHS) Data
###############################

# Define column names and types for the DHS file (narrowPeak format).
peak_colnames <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
                   "signalValue", "pValue", "qValue", "peak", "reads")
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

# Import DHS data.
dhs <- read_tsv(dhs_filepath, col_names = peak_colnames, col_types = peak_cols)

# Normalize reads to calculate RPKM.
dhs <- dhs %>% 
  mutate(peak_width = chromEnd - chromStart,
         reads_rpkm = reads / (sum(reads) / 1e6) / (peak_width / 1000),
         id = paste0(chrom, ":", chromStart, "_", chromEnd))

# Calculate quantiles and assign quantile labels.
q_thresholds <- seq(from = 0.0, to = 1, by = 0.25)
quants <- sapply(q_thresholds, function(q) quantile(as.numeric(dhs$reads), q))
dhs$quant <- sapply(dhs$reads, function(read) {
  if (read < quants[2]) {
    "0_25"
  } else if (read < quants[3]) {
    "25_50"
  } else if (read < quants[4]) {
    "50_75"
  } else {
    "75_1"
  }
})

# Plot DHS reads distribution.
p_dhs <- ggplot(dhs, aes(x = reads)) +
  geom_histogram(bins = 50, fill = "red") +
  labs(title = paste(nrow(dhs), "K562 DHS"), x = "DNase-seq reads in DHS") +
  theme_bw()
print(p_dhs)

###############################
# Define Locus Statistics Function
###############################

# This function calculates statistics for a given gene locus.
compute_locus_stats <- function(gene, genes_annot, peaks, tpm_data, 
                                tpm_threshold = 50, locus_width = 2e6) {
  
  # Create a locus window centered at the gene's TSS.
  locus <- resize(gene, width = locus_width, fix = "start")
  locusstart <- start(locus)
  locusend <- end(locus)
  
  # Find all genes and peaks that overlap the locus.
  locus_genes <- subsetByOverlaps(genes_annot, locus, ignore.strand = TRUE)
  locus_dhs <- subsetByOverlaps(peaks, locus, ignore.strand = TRUE)
  
  # Compute distances from DHS peaks to nearest TSS.
  x  <- distanceToNearest(ranges(locus_dhs), locus_genes$tss)
  locus_distance_nearest_tss_to_peak <- mcols(x)[,1]
  locus_nearest_tss <- locus_genes$gene_base_id[subjectHits(x)]
  
  # Identify genes with TPM data.
  locus_gene_ids <- locus_genes$gene_base_id 
  genes_with_tpm_data <- intersect(locus_gene_ids, tpm_data$gene)
  genes_without_tmp_data <- setdiff(locus_gene_ids, tpm_data$gene)
  
  # Get TPM values for genes.
  locus_tpm <- tpm_data %>% filter(gene %in% genes_with_tpm_data)
  tpm_values <- locus_tpm$tpm
  
  # Identify genes above the TPM threshold.
  genes_above_tpm <- locus_tpm %>% filter(tpm >= tpm_threshold) %>% pull(gene)
  
  # Create and return a tibble of locus statistics.
  output <- tibble(
    locus = gene$gene_base_id,
    locus_start = locusstart,
    locus_end = locusend,
    dhs = length(locus_dhs),
    locus_dhs_ids = list(locus_dhs$id),
    dhs_reads = list(locus_dhs$reads),
    dhs_quants = list(locus_dhs$quant),
    dhs_nearest_tss = list(locus_nearest_tss),
    dhs_nearest_tss_distance = list(locus_distance_nearest_tss_to_peak),
    genes = length(locus_genes),
    gene_ids = list(locus_gene_ids),
    gene_names = list(locus_genes$gene_name),
    gene_start = list(start(locus_genes)),
    genes_above_tpms = length(genes_above_tpm),
    genes_above_tpm = length(genes_above_tpm),
    genes_above_tpm_ids = list(genes_above_tpm),
    genes_with_tpm_data = length(genes_with_tpm_data),
    tpm_values = list(tpm_values),
    genes_without_tmp_data = length(genes_without_tmp_data)
  )
  return(output)
}

###############################
# Compute Locus Statistics for Candidate Genes
###############################

# Filter TPM data to genes with TPM >= threshold.
tpm_genes <- tpm %>% filter(tpm >= tpm_threshold)

# Select protein-coding genes on autosomes and X from v26.
genes <- v26[v26$type == "gene" &
               seqnames(v26) %in% paste0("chr", c(1:22, "X")) &
               v26$gene_type == "protein_coding"]

# Add a TSS column (strand-aware, here simplified by using the gene start).
genes$tss <- IRanges(start(resize(genes, 1)), start(resize(genes, 1)) + 1)

# Subset annotations for genes with TPM data.
tpm_genes_annot <- genes[genes$gene_base_id %in% tpm_genes$gene,]
tpm_genes <- tpm_genes %>% filter(gene %in% tpm_genes_annot$gene_base_id)
tpm_genes_annot <- split(tpm_genes_annot, f = tpm_genes_annot$gene_base_id)

# Compute locus statistics in parallel.
locus_stats <- lapply(tpm_genes_annot, FUN = compute_locus_stats,
                        genes_annot = genes,
                        peaks = makeGRangesFromDataFrame(dhs, keep.extra.columns = TRUE),
                        tpm_data = tpm, tpm_threshold = tpm_threshold, locus_width = locus_width)
locus_stats <- bind_rows(locus_stats)

###############################
# Plot Locus Statistics
###############################

# Define colors for plotting.
colors <- c("DHS" = "red", "Genes" = "blue", "Genes above 50 TPM" = "darkblue")

# Plot histograms of locus statistics.
p_locus <- locus_stats %>% 
  select(locus, DHS = dhs, Genes = genes, `Genes above 50 TPM` = genes_above_tpm) %>% 
  pivot_longer(cols = -locus, names_to = "stat", values_to = "value") %>% 
  ggplot(aes(x = value, fill = stat)) +
  facet_wrap(~stat, scales = "free") +
  geom_histogram(bins = 20) +
  labs(title = paste("Locus stats for all candidate loci (", nrow(locus_stats), ")"),
       x = "Count", y = "Number of loci") +
  scale_fill_manual(values = colors) +
  theme_bw()
print(p_locus)

###############################
# Random Sampling of Loci
###############################

# Read rep6 loci file (if needed for filtering).
rep6_loci <- read.table(rep6_loci_filepath, header = FALSE, sep = "\t")
colnames(rep6_loci) <- c("chr", "start", "end", "seqnames")

# Filter locus statistics to those in rep6 (example).
locus_stats_rep6 <- locus_stats %>% filter(locus %in% rep6_loci$seqnames)
# Convert all list columns to a character string where each element is collapsed by ";"
locus_stats_rep6_flat <- locus_stats_rep6 %>%
  mutate(across(where(is.list), ~ sapply(., function(x) paste(x, collapse = ";"))))
# Now write out the flattened tibble to a file.
write.table(locus_stats_rep6_flat, locus_stats_output, quote = FALSE, sep = "\t", row.names = FALSE)

# Randomly sample n_loci_sample loci, repeated n_sample_reps times.
set.seed(seed_value)
sample_reps <- replicate(sample_n(locus_stats, size = n_loci_sample), n = n_sample_reps, simplify = FALSE)
sample_reps <- bind_rows(sample_reps, .id = "rep") %>% mutate(rep = factor(rep, levels = unique(rep)))

# (Additional processing such as checking overlaps, sampling DHS peaks, and writing output files 
#  would continue here.)

###############################
# Additional Analyses and Plotting
###############################

# Example: Violin plot of locus stats across sampling repetitions.
p_violin <- sample_reps %>% 
  select(rep, locus, DHS = dhs, Genes = genes, `Genes above 50 TPM` = genes_above_tpm) %>% 
  pivot_longer(cols = -c(rep, locus), names_to = "stat", values_to = "value") %>% 
  ggplot(aes(x = rep, y = value, fill = stat)) +
  facet_wrap(~stat, scales = "free", ncol = 1) +
  geom_violin() +
  geom_jitter(width = 0.15, shape = 21, size = 1.2) +
  stat_summary(fun = median, geom = "point", size = 2.5, fill = "firebrick2", shape = 23) +
  labs(title = "Locus stats across repetitive sampling",
       x = "Repetition", y = "Count") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(legend.position = "none")
print(p_violin)

# Further downstream analyses, file outputs, and plotting code would follow here.

############################################################
# End of Script
############################################################