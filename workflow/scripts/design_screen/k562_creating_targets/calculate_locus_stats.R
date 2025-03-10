# Script: calculate_locus_stats.R
### SETUP =====================================================================
# Saving image for debugging
if (!file.exists("RDA_objects/k562_creating_targets")) { dir.create("RDA_objects/k562_creating_targets", recursive = TRUE) }
save.image(paste0("RDA_objects/k562_creating_targets/calculate_locus_stats.rda"))
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
  library(BiocParallel)
})

# Load input files
message("Loading input files...")
tpm_filtered <- readRDS(snakemake@input$tpm_filtered)
genes <- readRDS(snakemake@input$processed_annotation)
dhs_gr <- readRDS(snakemake@input$processed_peaks)

# Get parameters from Snakemake
locus_width <- snakemake@params$locus_width

### DEFINE LOCUS STATS FUNCTION ==============================================
message("Defining compute_locus_stats function...")
compute_locus_stats <- function(gene, genes_annot, peaks, tpm_data, tpm_threshold = 50,
                                locus_width = 2e6) {
  
  # Create locus window, centered on TSS of gene
  locus <- resize(gene, width = locus_width, fix = "start")
  locusstart <- as.numeric(start(locus))
  locusend <- as.numeric(end(locus))
  
  # Get all genes and peaks within that locus window
  locus_genes <- subsetByOverlaps(genes_annot, locus, ignore.strand = TRUE)
  locus_dhs <- subsetByOverlaps(peaks, locus, ignore.strand = TRUE)
  
  # Calculate distance from each peak to the nearest TSS
  x <- distanceToNearest(locus_dhs, locus_genes$tss_range)
  locus_distance_nearest_tss_to_peak <- mcols(x)[,1]
  locus_nearest_tss <- locus_genes$gene_base_id[subjectHits(x)]
  
  # Get genes with tpm data
  locus_gene_start <- start(locus_genes)
  locus_gene_ids <- locus_genes$gene_name
  locus_genes_ids <- locus_genes$gene_base_id
  genes_with_tpm_data <- intersect(locus_genes_ids, tpm_data$gene)
  genes_without_tmp_data <- setdiff(locus_genes_ids, tpm_data$gene)
  
  # Get tpm values for genes if available
  locus_tpm <- filter(tpm_data, gene %in% genes_with_tpm_data)
  tpm_values <- pull(locus_tpm, tpm)
  
  # Get genes above tpm threshold
  genes_above_tpm <- locus_tpm %>% 
    filter(tpm >= tpm_threshold) %>% 
    pull(gene)
  
  # Create output data frame
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
    genes = length(locus_genes_ids),
    gene_ids = list(locus_genes_ids),
    gene_names = list(locus_gene_ids),
    gene_start = list(locus_gene_start),
    genes_above_tpms = length(genes_above_tpm),
    genes_above_tpm_ids = list(genes_above_tpm),
    genes_with_tpm_data = length(genes_with_tpm_data),
    tpm_values = list(tpm_values),
    genes_without_tmp_data = length(genes_without_tmp_data)
  )
  
  return(output)
}

### PREPARE ANALYSIS =========================================================
message("Setting up for locus statistics computation...")
# Extract annotations for tpm filtered genes
tpm_genes_annot <- genes[genes$gene_base_id %in% tpm_filtered$gene]

message(paste("Number of candidate genes:", length(tpm_genes_annot)))

# Split tpm filtered gene annotations into GRangesList, one gene per element
tpm_genes_annot_list <- split(tpm_genes_annot, f = tpm_genes_annot$gene_base_id)

### COMPUTE LOCUS STATISTICS =================================================
# Register backend for parallel computing
register(MulticoreParam(workers = min(parallel::detectCores() - 1, 4)))

# Compute locus statistics for each potential target gene
message("Computing locus statistics in parallel (this may take some time)...")
locus_stats <- bplapply(tpm_genes_annot_list, FUN = compute_locus_stats, 
                        genes_annot = genes, peaks = dhs_gr, 
                        tpm_data = tpm_filtered, tpm_threshold = 50, 
                        locus_width = locus_width)

# Combine into one data frame
locus_stats_df <- bind_rows(locus_stats)

message(paste("Computed statistics for", nrow(locus_stats_df), "loci"))

### VISUALIZATION ============================================================
message("Creating locus statistics plots...")
# Colors for locus stats
colors <- c("DHS" = "#cc7000", "Genes" = "#4d9eff", "Genes above 50 TPM" = "#0052cc")

# Plot locus stats for all candidate loci
locus_hist_plot <- locus_stats_df %>% 
  select(locus, DHS = dhs, Genes = genes, `Genes above 50 TPM` = genes_above_tpms) %>% 
  pivot_longer(cols = -locus, names_to = "stat", values_to = "value") %>% 
  ggplot(aes(x = value, fill = stat)) +
  facet_wrap(~stat, scale = "free") +
  geom_histogram(bins = 20) +
  labs(title = paste0("Locus stats all candidate loci (", nrow(locus_stats_df), ")"),
       x = "Number of DHS/Genes/Genes above 50 TPM",
       y = "Number of loci") +
  scale_fill_manual(values = colors) +
  theme_bw()

# Extract DNase-seq reads in DHS
dnase_reads <- locus_stats_df %>% 
  select(locus, dhs_reads) %>% 
  unnest(cols = dhs_reads)

# Plot DNase-seq reads distribution across all DHS within loci
dnase_hist_plot <- ggplot(dnase_reads, aes(x = dhs_reads)) +
  geom_histogram(bins = 50, fill = "#cc7000") +
  labs(title = paste(nrow(dnase_reads), "K562 DHS within", nrow(locus_stats_df), "loci"),
       x = "DNase-seq reads in DHS") +
  theme_bw()

# Combine plots
combined_plot <- plot_grid(locus_hist_plot, dnase_hist_plot, 
                           ncol = 1, 
                           labels = c("A", "B"))

### SAVE OUTPUT ==============================================================
message("Saving output files")
saveRDS(locus_stats_df, file = snakemake@output$locus_stats)
ggsave(snakemake@output$locus_plots, plot = combined_plot, width = 10, height = 8)

### CLEAN UP =================================================================
message("Closing log file")
sink()
sink(type = "message")
close(log)