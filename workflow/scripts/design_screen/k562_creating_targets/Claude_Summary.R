#!/usr/bin/env Rscript

#####################################################################
# ENCODE CRISPRi Screen Design
# Organized from original_evvie_script.Rmd
#####################################################################

###############
# Load Libraries
###############
library(tidyverse)
library(here)
library(rtracklayer)
library(BiocParallel)
library(cowplot)
library(gridExtra)

###############
# Configuration Parameters
###############

# Random seed for reproducibility
RANDOM_SEED <- 69167

# Thresholds for filtering
TPM_THRESHOLD <- 50
LOCUS_WIDTH <- 2e6
DHS_THRESHOLD <- 20
KB100_THRESHOLD <- 10
Q75_THRESHOLD <- 20
NUM_LOCI <- 25

# File paths for input data
TPM_FILE <- "/Users/ejagoda/Documents/HGRM/EP_benchmarking/Encode_Crispri_Screen/tpm.csv"
DHS_FILE <- "/Users/ejagoda/Documents/HGRM/EP_benchmarking/Encode_Crispri_Screen/final_round_design/ENCFF185XRG_w_ENCFF325RTP_q30_sorted.txt"
PROMOTER_CLASSIFICATIONS_FILE <- "prom_gw_prediction.csv"
K562_CONTROL_GENES_FILE <- "/Users/ejagoda/Documents/HGRM/Encode_Crispri_Screen/K562_Positive_Controls_update.txt"

# Gencode annotation files
V26LIFT37_FILE <- "/Users/ejagoda/Documents/HGRM/EP_benchmarking/Encode_Crispri_Screen/gencode.v26lift37.annotation.gtf.gz"
V29_FILE <- "/Users/ejagoda/Documents/HGRM/Encode_Crispri_Screen/wct11/gencode.v29.annotation.gtf.gz"

# Output directories
OUTPUT_DIR <- "/Users/ejagoda/Documents/HGRM/Encode_Crispri_Screen/final_round_design/Potential_loci/"
RANDOM_DESIGN_DIR <- "/Users/ejagoda/Documents/HGRM/Encode_Crispri_Screen/totally_random_design/"

# URLs to GENCODE annotations (commented out, using local files instead)
# v26lift37_url <- "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh37_mapping/gencode.v26lift37.annotation.gtf.gz"
# v29_url <- "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz"






###############
# Input and Output Paths
###############

# Base directories
BASE_DIR <- "/Users/ejagoda/Documents/HGRM/EP_benchmarking/Encode_Crispri_Screen"
RANDOM_DESIGN_DIR <- "/Users/ejagoda/Documents/HGRM/Encode_Crispri_Screen/totally_random_design"
FINAL_DESIGN_DIR <- "/Users/ejagoda/Documents/HGRM/Encode_Crispri_Screen/final_round_design"

# Input files
INPUT_TPM_FILE <- file.path(BASE_DIR, "tpm.csv")
INPUT_DHS_FILE <- file.path(FINAL_DESIGN_DIR, "ENCFF185XRG_w_ENCFF325RTP_q30_sorted.txt")
INPUT_PROMOTER_CLASSIFICATIONS_FILE <- "prom_gw_prediction.csv"  # Relative path
INPUT_K562_CONTROL_GENES_FILE <- "/Users/ejagoda/Documents/HGRM/Encode_Crispri_Screen/K562_Positive_Controls_update.txt"
INPUT_V26LIFT37_FILE <- file.path(BASE_DIR, "gencode.v26lift37.annotation.gtf.gz")
INPUT_V29_FILE <- file.path(BASE_DIR, "wct11/gencode.v29.annotation.gtf.gz")
INPUT_REP6_LOCI_BED <- file.path(RANDOM_DESIGN_DIR, "potential_loci/loci_only_bed_version_rep6.bed")
INPUT_REP6_GENES_TAPSEQ <- file.path(RANDOM_DESIGN_DIR, "K562.rep6.GenesForTAPseq_wV29_genenames.txt")

# Output directories
OUTPUT_POTENTIAL_LOCI_DIR <- file.path(FINAL_DESIGN_DIR, "Potential_loci")
OUTPUT_RANDOM_DESIGN_DIR <- file.path(RANDOM_DESIGN_DIR, "potential_loci")

# Locus stats outputs
OUTPUT_REP6_LOCUS_INFO <- file.path(OUTPUT_RANDOM_DESIGN_DIR, "rep_6_locus_info_genes_>50tpm.txt")
OUTPUT_REP6_WITH_3GENES <- file.path(OUTPUT_RANDOM_DESIGN_DIR, "rep_6_loci_with_3_other_>50tpm_genes.txt")
OUTPUT_REP6_WITH_TPM_THRESHOLDS <- file.path(OUTPUT_RANDOM_DESIGN_DIR, "rep_6_loci_with_3_other_>50tpm_genes_w_g_40_w_g_30.txt")

# Gene list outputs
OUTPUT_GENE_LIST_TPM_TAB <- file.path(OUTPUT_RANDOM_DESIGN_DIR, "K562_rep_6_genes_g30_and_controls_tpm_tab_take2.txt")
OUTPUT_NEGATIVE_CONTROLS <- file.path(RANDOM_DESIGN_DIR, "k562_negative_controls_random_set.txt")
OUTPUT_FINAL_GENE_LIST <- file.path(RANDOM_DESIGN_DIR, "K562.rep6.GenesForTAPseq_wV29_genenames_controls_fixed.txt")

# Sample iteration output templates (use sprintf or paste0 with iteration number)
OUTPUT_TEMPLATE_LOCI_BED <- file.path(OUTPUT_POTENTIAL_LOCI_DIR, "loci_only_bed_version_rep%d.bed")
OUTPUT_TEMPLATE_PEAKS_TXT <- file.path(OUTPUT_POTENTIAL_LOCI_DIR, "sampled_peaks_rep%d.txt")
OUTPUT_TEMPLATE_LOCI_TXT <- file.path(OUTPUT_POTENTIAL_LOCI_DIR, "sampled_loci_rep%d.txt")

# Best loci bed file
OUTPUT_BEST_LOCI_BED <- file.path(FINAL_DESIGN_DIR, "best_loci_loci_only_for_bedtools_intersect.bed")

# Create output directories if they don't exist
dir.create(OUTPUT_POTENTIAL_LOCI_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_RANDOM_DESIGN_DIR, recursive = TRUE, showWarnings = FALSE)






###############
# Load and Process Data
###############

# Load TPM data
tpm <- read_csv(TPM_FILE, col_types = cols(.default = col_character(), tpm = col_double()))

# Import genome annotations
v26 <- rtracklayer::import(V26LIFT37_FILE)
v29 <- rtracklayer::import(V29_FILE)

# Get Ensembl gene ids without version for v26
v26$gene_base_id <- sub("\\..*", "", v26$gene_id)

# Column names in a narrowPeaks bed file
peak_colnames <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue",
                   "pValue", "qValue", "peak", "reads")

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

# Import ENCODE DNase-seq peaks
dhs <- read_tsv(DHS_FILE, col_names = peak_colnames, col_types = peak_cols)

# Normalize reads for peak width by calculating RPKM
dhs <- dhs %>% 
  mutate(peak_width = chromEnd - chromStart,
         reads_rpkm = reads / (sum(reads) / 1e6) / (peak_width / 1000))

dhs$id = paste0(dhs$chrom,":",dhs$chromStart,"_",dhs$chromEnd)

# Assign quantile categories to DHS peaks
q_thresholds = seq(from = 0.0, to = 1, by = 0.25)
quants = c()
for (q in q_thresholds){
  quant = quantile(as.numeric(paste0(dhs$reads)), q)
  quants = c(quants, quant)
}

dhs$quant = "x"
for (i in 1:nrow(dhs)){
  read = as.numeric(paste0(dhs[i,"reads"]))
  if (read < quants[2]){
    dhs[i,"quant"] = "0_25"
  }
  else if (read >= quants[2] & read < quants[3]){
    dhs[i,"quant"] = "25_50"
  }
  else if (read >= quants[3] & read < quants[4]){
    dhs[i,"quant"] = "50_75"
  }
  else if (read >= quants[4] & read < quants[5]){
    dhs[i,"quant"] = "75_1"
  }
}

# Create GRanges object for DNase peaks
dhs_gr <- makeGRangesFromDataFrame(dhs, keep.extra.columns = TRUE)

###############
# Function Definitions
###############

# Function to count all genes and peaks within a certain distance of a target gene
compute_locus_stats <- function(gene, genes_annot, peaks, tpm_data, tpm_threshold = 50,
                                locus_width = 2e6) {
  
  # Create locus window, centered on TSS of gene
  locus <- resize(gene, width = locus_width, fix = "start")
  locusstart = as.numeric(paste0(start(locus)))
  locusend = as.numeric(paste0(end(locus)))
  
  # Get all genes and peaks within that locus window
  locus_genes <- subsetByOverlaps(genes_annot, locus, ignore.strand = TRUE)
  locus_dhs <- subsetByOverlaps(peaks, locus, ignore.strand = TRUE)
  
  # Calculate distance from each peak to nearest TSS
  x = distanceToNearest(ranges(locus_dhs), locus_genes$tss)
  locus_distance_nearest_tss_to_peak <- mcols(x)[,1]
  locus_nearest_tss = locus_genes$gene_base_id[subjectHits(x)]
  
  # Get genes with tpm data
  locus_gene_start <- start(locus_genes)
  locus_gene_ids <- locus_genes$gene_name 
  locus_genes <- locus_genes$gene_base_id
  genes_with_tpm_data <- intersect(locus_genes, tpm$gene)
  genes_without_tmp_data <- setdiff(locus_genes, tpm$gene)
  
  # Get tpm values for genes if available
  locus_tpm <- filter(tpm_data, gene %in% genes_with_tpm_data)
  tpm_values <- pull(locus_tpm, tpm)
  
  # Get genes above tpm threshold
  genes_above_tpm <- locus_tpm %>% 
    filter(tpm >= tpm_threshold) %>% 
    pull(gene)
  
  # Create output data frame with all the stats
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
    gene_ids = list(locus_genes),
    gene_names = list(locus_gene_ids),
    gene_start = list(locus_gene_start),
    genes_above_tpms = length(genes_above_tpm),
    genes_above_tpm = length(genes_above_tpm),
    genes_above_tpm_ids = list(genes_above_tpm),
    genes_with_tpm_data = length(genes_with_tpm_data),
    tpm_values = list(tpm_values),
    genes_without_tmp_data = length(genes_without_tmp_data)
  )
  
  return(output)
}

###############
# Process and Filter Genes
###############

# Get all genes above TPM threshold
tpm_genes <- filter(tpm, tpm >= 50)

# Only retain protein-coding gene locus annotations on autosomes and chromosome X
genes <- v26[v26$type == "gene" &
               seqnames(v26) %in% paste0("chr", c(1:22, "X")) &
               v26$gene_type == "protein_coding"]

# Add the strand oriented start position as tss
genes$tss = IRanges(start(resize(genes,1)), start(resize(genes,1)) + 1) 

# Extract annotations for tpm filtered genes
tpm_genes_annot <- genes[genes$gene_base_id %in% tpm_genes$gene]
tpm_genes <- filter(tpm_genes, gene %in% tpm_genes_annot$gene_base_id)

# Split tpm filtered gene annotations into GRangesList, one gene per element
tpm_genes_annot <- split(tpm_genes_annot, f = tpm_genes_annot$gene_base_id)

# Register backend for parallel computing
register(MulticoreParam(workers = 5))

# Compute locus statistics for each potential target gene
locus_stats <- bplapply(tpm_genes_annot, FUN = compute_locus_stats, genes_annot = genes,
                        peaks = dhs_gr, tpm_data = tpm, tpm_threshold = 50, locus_width = 2e6)

# Combine into one data frame
locus_stats <- bind_rows(locus_stats)
locus_stats_df <- data.frame(locus_stats)

###############
# Analyze DHS Peaks and TSS Distances
###############

# Calculate statistics for DHS peaks by quantile
locus_stats_df$q75counts <- "x"
locus_stats_df$bottom75 <- "x"
locus_stats_df$q75counts_p <- "x"
locus_stats_df$bottom75_p <- "x"

for (i in 1:nrow(locus_stats_df)){
  quant_list <- locus_stats_df[i,"dhs_quants"][[1]]
  if (is.na(as.numeric(table(quant_list)["75_1"]))){
    locus_stats_df[i,"q75counts"] <- 0
  }
  else{
    locus_stats_df[i,"q75counts"] <- as.numeric(table(quant_list)["75_1"])
    locus_stats_df[i,"q75counts_p"] <- as.numeric(table(quant_list)["75_1"])/length(quant_list)
  }
  locus_stats_df[i,"bottom75"] <- sum(table(quant_list)) - as.numeric(table(quant_list)["75_1"])
  locus_stats_df[i,"bottom75_p"] <- as.numeric(paste0(locus_stats_df[i,"bottom75"]))/length(quant_list)
}

# Calculate TSS distance statistics for DHS peaks
locus_stats_df$g100kb <- "x"
locus_stats_df$l100kb <- "x"
locus_stats_df$g100kb_p <- "x"
locus_stats_df$l100kb_p <- "x"

for (i in 1:nrow(locus_stats_df)){
  distances <- locus_stats_df$dhs_nearest_tss_distance[[i]]
  if (length(distances) == 0){
    locus_stats_df[i,"g100kb"] <- 0
    locus_stats_df[i,"l100kb"] <- 0
  }
  else{
    g100kb_count <- length(which(distances > 100000))
    l100kb_count <- length(which(distances <= 100000 & distances > 1000)) # Want at least 1kb from promoter
    locus_stats_df[i,"g100kb"] <- g100kb_count
    locus_stats_df[i,"l100kb"] <- l100kb_count
  }
  locus_stats_df[i,"g100kb_p"] <- as.numeric(paste0(locus_stats_df[i,"g100kb"]))/length(distances)
  locus_stats_df[i,"l100kb_p"] <- as.numeric(paste0(locus_stats_df[i,"l100kb"]))/length(distances)
}

###############
# Filter Loci Based on Criteria
###############

# Get the potentially usable loci based on DHS quantiles and TSS distances
usable_loci <- locus_stats_df[
  as.numeric(paste0(locus_stats_df$q75counts)) >= DHS_THRESHOLD & 
    as.numeric(paste0(locus_stats_df$bottom75)) >= DHS_THRESHOLD & 
    as.numeric(paste0(locus_stats_df$g100kb)) >= KB100_THRESHOLD & 
    as.numeric(paste0(locus_stats_df$l100kb)) >= 30, 
]

###############
# Sample Peaks from Usable Loci
###############

# For each usable locus, sample peaks based on distance and quantile criteria
usable_loci$peak_sample <- "x"
usable_loci$peak_sample_length <- "x"

for (i in 1:nrow(usable_loci)){
  set.seed(100)
  distances <- usable_loci$dhs_nearest_tss_distance[[i]]
  quant_list <- usable_loci$dhs_quants[[i]]
  peak_list <- usable_loci$locus_dhs_ids[[i]]
  
  # Sample peaks > 100kb from TSS
  g100_kb_peak_indexs <- which(distances > 100000)
  g100_kb_sampled_indexs <- sample(g100_kb_peak_indexs, size = min(10, length(g100_kb_peak_indexs)))
  g100kb_sampled_peaks <- peak_list[g100_kb_sampled_indexs]
  q75_sampled_already <- length(which(quant_list[g100_kb_sampled_indexs] == "75_1"))
  
  # Sample peaks <= 100kb from TSS
  l100_kb_peak_indexs <- which(distances <= 100000 & distances > 100)
  l100kb_peaks <- peak_list[l100_kb_peak_indexs]
  quants_l100kb_peaks <- quant_list[l100_kb_peak_indexs]
  
  # Sample top 75% quantile peaks
  top75_l100kb <- l100kb_peaks[which(quants_l100kb_peaks == "75_1")]
  bottom75_l100kb <- l100kb_peaks[which(quants_l100kb_peaks != "75_1")]
  
  if (length(top75_l100kb) >= (20 - q75_sampled_already)){
    top75_l100kb_sampled <- sample(top75_l100kb, size = 20 - q75_sampled_already)
  } else {
    top75_l100kb_sampled <- top75_l100kb
  }
  
  # Sample bottom 75% quantile peaks
  if (length(bottom75_l100kb) >= (40 - (length(top75_l100kb_sampled) + q75_sampled_already) - (10-q75_sampled_already))){
    bottom75_l100kb_sampled <- sample(bottom75_l100kb, size = (40 - (length(top75_l100kb_sampled) + q75_sampled_already) - (10-q75_sampled_already)))
  } else {
    bottom75_l100kb_sampled <- bottom75_l100kb
  }
  
  # Combine all sampled peaks
  full_sample <- unique(c(g100kb_sampled_peaks, top75_l100kb_sampled, bottom75_l100kb_sampled))
  usable_loci[i,"peak_sample"] <- paste0(full_sample, collapse = ";")
  usable_loci[i,"peak_sample_length"] <- length(full_sample)
}

###############
# Calculate Stats for Sampled Peaks
###############

usable_loci$sampled_stats <- "x"
usable_loci$sampled_quants <- "x" 
usable_loci$sampled_nearest_tss_distance <- "x"
usable_loci$sampled_reads <- "x"

for (i in 1:nrow(usable_loci)){
  distances <- usable_loci$dhs_nearest_tss_distance[[i]]
  quant_list <- usable_loci$dhs_quants[[i]]
  reads_list <- usable_loci$dhs_reads[[i]]
  all_locus_peak_list <- usable_loci$locus_dhs_ids[[i]]
  sampled_peak_list <- str_split(usable_loci$peak_sample[[i]], pattern = ";")[[1]]
  
  sampled_peak_indexs <- c()
  for (j in 1:length(sampled_peak_list)){
    peak <- sampled_peak_list[j]
    peak_index <- match(peak, all_locus_peak_list)
    sampled_peak_indexs <- c(sampled_peak_indexs, peak_index)
  }
  
  # Calculate statistics for sampled peaks
  sampled_distances <- distances[sampled_peak_indexs]
  g100 <- length(which(sampled_distances > 100000))
  l100 <- length(which(sampled_distances <= 100000 & sampled_distances > 1000))
  sampled_quants <- quant_list[sampled_peak_indexs]
  sampled_reads <- reads_list[sampled_peak_indexs]
  q75 <- length(which(sampled_quants == "75_1"))
  b75 <- length(which(sampled_quants != "75_1"))
  
  usable_loci[i,"sampled_stats"] <- paste0(c(g100, l100, q75, b75), collapse = ";")
  usable_loci[i,"sampled_quants"] <- paste0(sampled_quants, collapse = ";")
  usable_loci[i,"sampled_reads"] <- paste0(sampled_reads, collapse = ";")
  usable_loci[i,"sampled_nearest_tss_distance"] <- paste0(sampled_distances, collapse = ";")
}

###############
# Filter for Best Loci
###############

# Select loci with ideal stats (10 >100kb, 30 <=100kb, 20 top75%, 20 bottom75%)
best_loci <- usable_loci[usable_loci$sampled_stats == "10;30;20;20", ]

# Add chromosome information
best_loci$chr <- "x"
for (i in 1:nrow(best_loci)){
  chr <- str_split(best_loci[i,"locus_dhs_ids"][[1]][1], ":")[[1]][1]
  best_loci[i,"chr"] <- chr
}

###############
# Sample Loci and Generate Outputs
###############

# Create BED format output for best loci
bed_version <- best_loci[, c("chr", "locus_start", "locus_end", "locus")]

# Run multiple iterations of loci sampling
for (k in 0:100) {
  set.seed(RANDOM_SEED + k)
  
  # Sample 25 loci
  sample_best_loci_indexes <- sample(1:nrow(best_loci), size = NUM_LOCI)
  best_loci_sample_1 <- best_loci[sample_best_loci_indexes, ]
  
  # Output BED file for sampled loci
  write.table(
    bed_version[bed_version$locus %in% best_loci_sample_1$locus, ],
    paste0(OUTPUT_DIR, "loci_only_bed_version_rep", k, ".bed"),
    sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  
  # Prepare gene information for output
  best_loci_sample_1$gene_ids_print <- "x"
  best_loci_sample_1$gene_names_print <- "x"
  
  for (i in 1:nrow(best_loci_sample_1)) {
    gene_ids <- best_loci_sample_1[i,"gene_ids"][[1]]
    best_loci_sample_1[i,"gene_ids_print"] <- paste0(gene_ids, collapse = ";")
    
    gene_names <- best_loci_sample_1[i,"gene_names"][[1]]
    best_loci_sample_1[i,"gene_names_print"] <- paste0(gene_names, collapse = ";")
  }
  
  # Output tab-delimited file with peak information
  outfile_name <- paste0(OUTPUT_DIR, "sampled_peaks_rep", k, ".txt")
  write(
    paste("locus", "peak", "peak_reads", "peak_quant", "peak_nearest_tss_distance", "locus_genes_id", "locus_gene_names", sep = "\t"),
    outfile_name
  )
  
  for (i in 1:nrow(best_loci_sample_1)) {
    peaks <- str_split(best_loci_sample_1[i,"peak_sample"], pattern = ";")[[1]]
    sampled_reads <- str_split(best_loci_sample_1[i,"sampled_reads"], pattern = ";")[[1]]
    sampled_quants <- str_split(best_loci_sample_1[i,"sampled_quants"], pattern = ";")[[1]]
    sampled_nearest_tss_distance <- str_split(best_loci_sample_1[i,"sampled_nearest_tss_distance"], pattern = ";")[[1]]
    
    for (j in 1:length(peaks)) {
      write(
        paste(
          best_loci_sample_1[i,"locus"],
          peaks[j],
          sampled_reads[j],
          sampled_quants[j],
          sampled_nearest_tss_distance[j],
          best_loci_sample_1[i,"gene_ids_print"],
          best_loci_sample_1[i,"gene_names_print"],
          sep = "\t"
        ),
        outfile_name, append = TRUE
      )
    }
  }
  
  # Output summarized loci information
  best_loci_sample_1_printv1 <- data.frame(
    cbind(
      best_loci_sample_1$locus,
      best_loci_sample_1$dhs,
      best_loci_sample_1$gene_ids_print,
      best_loci_sample_1$gene_names_print,
      best_loci_sample_1$peak_sample,
      best_loci_sample_1$sampled_reads,
      best_loci_sample_1$sampled_quants,
      best_loci_sample_1$sampled_nearest_tss_distance
    )
  )
  
  colnames(best_loci_sample_1_printv1) <- c(
    "locus_central_gene",
    "total_dhs",
    "locus_gene_ids",
    "locus_gene_names",
    "sampled_peaks",
    "sampled_peak_DHSreads",
    "sampled_peak_quants",
    "sampled_peaks_distance_to_nearest_tss"
  )
  
  write.table(
    best_loci_sample_1_printv1,
    paste0(OUTPUT_DIR, "sampled_loci_rep", k, ".txt"),
    quote = FALSE, sep = "\t", row.names = FALSE
  )
}

###############
# Process Specific Rep 6
###############

# Process specific rep 6 loci (if available)
rep_6_loci <- read.table(
  paste0(RANDOM_DESIGN_DIR, "potential_loci/loci_only_bed_version_rep6.bed"),
  header = FALSE, sep = '\t'
)
colnames(rep_6_loci) <- c("chr", "start", "end", "seqnames")

locus_stats_rep6 <- locus_stats_df[locus_stats_df$locus %in% rep_6_loci$seqnames, ]
locus_stats_rep6_small <- locus_stats_rep6[, c("locus", "genes_above_tpms", "genes_above_tpm_ids")]

# Format genes above TPM
locus_stats_rep6_small$genes_above_tpms_list <- sapply(
  locus_stats_rep6_small$genes_above_tpm_ids, 
  function(x) paste(x, collapse = ";")
)
locus_stats_rep6_small$genes_above_tpm_ids <- NULL

# Save to file
write.table(
  locus_stats_rep6_small,
  paste0(RANDOM_DESIGN_DIR, "potential_loci/rep_6_locus_info_genes_>50tpm.txt"),
  quote = FALSE, sep = '\t', row.names = FALSE
)

# Additional processing for rep 6
rep6_tab <- read.table(
  paste0(RANDOM_DESIGN_DIR, "potential_loci/rep_6_locus_info_genes_>50tpm.txt"),
  header = TRUE, sep = "\t"
)

# Clean up gene list format
good_tpm_list <- gsub("c[(]", "", rep6_tab$genes_above_tpms_list)
good_tpm_listt <- gsub("[)]", "", good_tpm_list)
rep6_tab$genes_above_tpms_list_good <- good_tpm_listt
rep6_tab$genes_above_tpms_list <- NULL

# Add gene names and sample additional genes
set.seed(20)
rep6_tab$additional_3genes <- "x"
rep6_tab$locus_gene_name <- "x"
v26_df <- data.frame(v26)

for (i in 1:nrow(rep6_tab)) {
  locus_gene <- rep6_tab[i, "locus"]    
  locus_gene_name <- unique(v26_df[v26_df$gene_base_id == locus_gene, "gene_name"])
  rep6_tab[i, "locus_gene_name"] <- locus_gene_name
  
  if (rep6_tab[i, "genes_above_tpms"] != 1) {
    gene_list <- str_split(rep6_tab$genes_above_tpms_list_good[i], pattern = ", ")[[1]]
    other_genes_list <- gene_list[which(grepl(gene_list, pattern = locus_gene) == FALSE)]
    other_genes_list_gene_name <- c()
    
    for (gene in other_genes_list) {
      gene_name <- unique(v26_df[v26_df$gene_base_id == gene, "gene_name"])
      other_genes_list_gene_name <- c(other_genes_list_gene_name, gene_name)
    }
    
    if (length(other_genes_list) <= 3) {
      rep6_tab[i, "additional_3genes"] <- paste(other_genes_list_gene_name, collapse = ";")
    } else {
      sampled_other_genes <- sample(other_genes_list_gene_name, 3)
      rep6_tab[i, "additional_3genes"] <- paste(sampled_other_genes, collapse = ";")
    }
  }
}

# Save rep 6 with additional genes
write.table(
  rep6_tab,
  paste0(RANDOM_DESIGN_DIR, "potential_loci/rep_6_loci_with_3_other_>50tpm_genes.txt"),
  quote = FALSE, sep = '\t', row.names = FALSE
)

# Get TPM statistics for different thresholds
rep6_tab$n_genes_g_40tpm <- "x"
rep6_tab$n_genes_g_30tpm <- "x"
rep6_tab$genes_g_40tpm <- "x"
rep6_tab$genes_g_30tpm <- "x"

for (i in 1:nrow(rep6_tab)) {
  locus_gene <- rep6_tab[i, "locus"]
  locus_gene_name <- unique(v26_df[v26_df$gene_base_id == locus_gene, "gene_name"])
  locus_stats_tab <- locus_stats_df[locus_stats_df$locus == locus_gene, ]
  gene_list <- locus_stats_tab[1, "gene_names"][[1]]
  gene_tpm_tab <- data.frame(tpm[tpm$gene_name %in% gene_list, ])
  
  rep6_tab[i, "n_genes_g_40tpm"] <- nrow(gene_tpm_tab[gene_tpm_tab$tpm >= 40, ])
  rep6_tab[i, "genes_g_40tpm"] <- paste0(gene_tpm_tab[gene_tpm_tab$tpm >= 40, "gene_name"], collapse = ";")
  rep6_tab[i, "n_genes_g_30tpm"] <- nrow(gene_tpm_tab[gene_tpm_tab$tpm >= 30, ])
  rep6_tab[i, "genes_g_30tpm"] <- paste0(gene_tpm_tab[gene_tpm_tab$tpm >= 30, "gene_name"], collapse = ";")
}

# Save with additional TPM thresholds
write.table(
  rep6_tab,
  paste0(RANDOM_DESIGN_DIR, "potential_loci/rep_6_loci_with_3_other_>50tpm_genes_w_g_40_w_g_30.txt"),
  quote = FALSE, sep = '\t', row.names = FALSE
)

###############
# Select Genes for TAPseq
###############

# Extract genes for experiment
genes_k <- c()
for (i in 1:nrow(rep6_tab)) {
  genes_k <- c(genes_k, rep6_tab[i, "locus_gene_name"])
  ad_genes <- rep6_tab[i, "genes_g_30tpm"] 
  if (ad_genes != "x") {
    other_genes <- str_split(ad_genes, pattern = ";")[[1]] 
    genes_k <- c(genes_k, other_genes)
  }
}

# Get TPM data for selected genes
exper_tab <- data.frame(tpm[tpm$gene_name %in% genes_k, ])
exper_tab$type <- "x"
for (i in 1:nrow(exper_tab)) {
  if (exper_tab[i, "gene_name"] %in% rep6_tab$locus_gene_name) {
    exper_tab[i, "type"] <- "central_gene"
  } else {
    exper_tab[i, "type"] <- "locus_gene"
  }
}

# Process control genes
k562_control_genes <- read.table(K562_CONTROL_GENES_FILE, header = FALSE, sep = '\t')
k562_control_genes$V2 <- NULL
colnames(k562_control_genes) <- "gene"
k562_control_genes_w_tpm <- tpm[tpm$gene_name %in% k562_control_genes$gene, ]
k562_control_genes_w_tpm$type <- "control"

# Combine experimental and control genes
full_k562_tpm_tab <- data.frame(rbind(exper_tab, k562_control_genes_w_tpm))

# Output combined gene table
write.table(
  full_k562_tpm_tab,
  paste0(RANDOM_DESIGN_DIR, "potential_loci/K562_rep_6_genes_g30_and_controls_tpm_tab_take2.txt"),
  quote = FALSE, sep = "\t", row.names = FALSE
)

###############
# Create Negative Controls
###############

# Filter for protein-coding genes in v29
v29_protein <- v29[v29$gene_type == "protein_coding" & 
                     v29$type == "gene" &
                     seqnames(v29) %in% paste0("chr", c(1:22, "X")), ]

# Select negative controls in specific TPM range
tpm_subset <- tpm[tpm$tpm >= 30 & 
                    tpm$tpm <= 150 & 
                    tpm$gene_name %in% full_k562_tpm_tab$gene_name == FALSE & 
                    tpm$gene_name %in% v29_protein$gene_name, ]

# Sample 50 negative controls
set.seed(100)
tpm_subset_50 <- tpm_subset[sample(1:nrow(tpm_subset), size = 50), ]

# Add type and use flags
tpm_subset_50$type <- "neg_control"
tpm_subset_50$use <- "yes"

# Read file with genes for TAPseq (if available)
# This assumes the file exists, otherwise this part would need to be modified
full_set_use <- read.table(
  paste0(RANDOM_DESIGN_DIR, "K562.rep6.GenesForTAPseq_wV29_genenames.txt"),
  header = TRUE, sep = '\t'
)

# Add v29 gene names
full_set_use$gene_name_v29 <- "x"
for (i in 1:nrow(full_set_use)) {
  ensemble_name <- full_set_use[i, "gene"]
  name_v29 <- v29[v29$gene_base_id == ensemble_name, "gene_name"][1]
  full_set_use[i, "gene_name_v29"] <- name_v29
}

# Find and replace non-protein-coding genes
# In the original script, RAB30-AS1 was identified as a lincRNA and replaced with SNAPIN
non_protein <- full_set_use[full_set_use$gene_name_v29 %in% setdiff(full_set_use$gene_name_v29, v29_protein$gene_name), ]

# Add replacement gene (SNAPIN)
tpm_add <- tpm[tpm$gene_name == "SNAPIN", ]
tpm_add$type <- "neg_control"
tpm_add$use <- "yes"
tpm_add$gene_name_v29 <- tpm_add$gene_name

# Create final gene set for TAPseq
full_set_use2 <- rbind(full_set_use, tpm_add)
full_set_use3 <- full_set_use2[full_set_use2$gene_name != "RAB30-AS1", ]

# Output final gene set
write.table(
  full_set_use3,
  paste0(RANDOM_DESIGN_DIR, "K562.rep6.GenesForTAPseq_wV29_genenames_controls_fixed.txt"),
  quote = FALSE, sep = '\t', row.names = FALSE
)

###############
# Visualization Functions (Used Throughout the Workflow)
###############

# These plotting functions can be called at various points in the workflow

# Plot TPM distribution
plot_tpm_distribution <- function() {
  # Number and percentage of genes above TPM threshold
  n_above_tpm <- sum(tpm$tpm >= TPM_THRESHOLD)
  
  # Plot TPM distribution
  ggplot(tpm, aes(tpm)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = TPM_THRESHOLD, color = "red") +
    labs(x = "Transcripts-per-million (TPM)", 
         title = paste0("Genes above ", TPM_THRESHOLD, " TPM: ", n_above_tpm)) +
    scale_x_log10() +
    theme_bw()
}

# Plot DNase reads distribution
plot_dhs_distribution <- function() {
  ggplot(dhs, aes(x = reads)) +
    geom_histogram(bins = 50, fill = hsv(h = 0.1, s = 0.8, alpha = 1)) +
    labs(title = paste(nrow(dhs), "K562 DHS"), x = "DNase-seq reads in DHS") +
    theme_bw()
}

# Plot DNase reads cumulative distribution
plot_dhs_cumulative <- function() {
  ggplot(dhs, aes(x = reads)) + 
    stat_ecdf(geom = "step") +
    labs(title = "K562 DHS", x = "DNase-seq reads in DHS", y = "cumulative fraction") +
    theme_bw()
}

# Plot RPKM normalized DNase reads cumulative distribution
plot_dhs_rpkm_cumulative <- function() {
  ggplot(dhs, aes(x = reads_rpkm)) + 
    stat_ecdf(geom = "step") +
    labs(title = "rpkm normalized K562 DHS", 
         x = "rpkm normalized DNase-seq reads in DHS", 
         y = "cumulative fraction") +
    theme_bw()
}

# Plot locus stats
plot_locus_stats <- function() {
  # Colors for locus stats
  colors <- c("DHS" = hsv(h = 0.1, s = 0.8, alpha = 1),
              "Genes" = hsv(h = 0.62, s = 0.3, alpha = 1), 
              "Genes above 50 TPM" = hsv(h = 0.62, s = 0.8, alpha = 1))
  
  # Plot locus stats for all candidate loci
  locus_stats %>% 
    dplyr::select(locus, DHS = dhs, Genes = genes, `Genes above 50 TPM` = genes_above_tpm) %>% 
    pivot_longer(cols = -locus, names_to = "stat", values_to = "value") %>% 
    ggplot(., aes(x = value, fill = stat)) +
    facet_wrap(~stat, scale = "free") +
    geom_histogram(bins = 20) +
    labs(title = paste0("Locus stats all candidate loci (", nrow(locus_stats), ")"),
         x = "Number of DHS/Genes/Genes above 50 TPM",
         y = "Number of loci") +
    scale_fill_manual(values = colors) +
    theme_bw()
}

# Plot all DNase reads in loci
plot_all_dnase_reads_in_loci <- function() {
  # Extract DNase-seq reads in DHS
  dnase_reads <- locus_stats %>% 
    dplyr::select(locus, dhs_reads) %>% 
    unnest(cols = dhs_reads)
  
  # Plot DNase-seq reads distribution across all DHS within loci
  ggplot(dnase_reads, aes(x = dhs_reads)) +
    geom_histogram(bins = 50, fill = hsv(h = 0.1, s = 0.8, alpha = 1)) +
    labs(title = paste(nrow(dnase_reads), "K562 DHS within", nrow(locus_stats), "loci"),
         x = "DNase-seq reads in DHS") +
    theme_bw()
}

# Compare usable vs non-usable loci
plot_usable_vs_all <- function() {
  # Add flag for best loci
  locus_stats_df$isbest <- 0
  for (i in 1:nrow(locus_stats_df)) {
    if (paste0(locus_stats_df[i, "locus"]) %in% paste0(best_loci$locus)) {
      locus_stats_df[i, "isbest"] <- 1
    }
  }
  
  # Sort by percentile of q75 peaks
  locus_stats_df_odhs <- locus_stats_df[order(as.numeric(paste0(locus_stats_df$q75counts_p)), decreasing = FALSE), ]
  locus_stats_df_odhs$dhs_rank <- 1:nrow(locus_stats_df_odhs)
  
  # Plot the distribution of q75 peaks percentage
  plot(
    locus_stats_df_odhs$dhs_rank, 
    as.numeric(paste0(locus_stats_df_odhs$q75counts_p)), 
    pch = 16, cex = 0.5, 
    main = "All TPM > 50 loci", 
    ylab = "% q75 peaks", 
    xlab = "rank"
  )
  legend("topleft", fill = c("black", "red"), legend = c("all", "usable"))
  
  # Highlight usable loci
  best <- locus_stats_df_odhs[locus_stats_df_odhs$isbest == 1, ]
  points(
    best$dhs_rank, 
    as.numeric(paste0(best$q75counts_p)), 
    pch = 16, 
    cex = 0.5, 
    col = "red"
  )
  
  # Boxplot comparisons
  par(mfrow = c(2, 2))
  boxplot(
    as.numeric(paste0(locus_stats_df$q75counts_p)) ~ locus_stats_df$isbest, 
    ylab = "% DHS q75", 
    xlab = "usable"
  )
  
  boxplot(
    as.numeric(paste0(locus_stats_df$g100kb_p)) ~ locus_stats_df$isbest, 
    ylab = "% >100kb tss", 
    xlab = "usable"
  )
  
  boxplot(
    as.numeric(paste0(locus_stats_df$dhs)) ~ locus_stats_df$isbest, 
    ylab = "n_peaks", 
    xlab = "usable"
  )
  
  boxplot(
    as.numeric(paste0(locus_stats_df$genes)) ~ locus_stats_df$isbest, 
    ylab = "n_genes", 
    xlab = "usable"
  )
  par(mfrow = c(1, 1))
}

# Plot gene expression
plot_gene_expression <- function() {
  # Sort by TPM
  full_k562_tpm_tab_o <- full_k562_tpm_tab[order(full_k562_tpm_tab$tpm), ]
  full_k562_tpm_tab_o$rank <- 1:nrow(full_k562_tpm_tab_o)
  
  # Plot all genes
  p1 <- ggplot(full_k562_tpm_tab_o, aes(x = as.numeric(paste0(rank)), y = as.numeric(paste0(tpm)))) +
    geom_point() +
    geom_text(
      data = full_k562_tpm_tab_o[full_k562_tpm_tab_o$tpm > 1000, ],
      aes(label = paste0(gene_name, "_", type)),
      hjust = 1.2,
      check_overlap = FALSE
    ) +
    ylab("TPM") +
    xlab("rank") +
    theme_bw()
  
  # Plot high expression genes
  pz <- ggplot(full_k562_tpm_tab_o[full_k562_tpm_tab_o$tpm > 200, ], 
               aes(x = as.numeric(paste0(rank)), y = as.numeric(paste0(tpm)))) +
    geom_point() +
    geom_text(
      data = full_k562_tpm_tab_o[full_k562_tpm_tab_o$tpm > 1000, ],
      aes(label = paste0(gene_name, "_", type)),
      hjust = 1.2
    ) +
    ylab("TPM") +
    xlab("rank") +
    theme_bw()
  
  grid.arrange(p1, pz, nrow = 1)
}

###############
# Final Output Message
###############

cat("\nScript completed successfully!\n")
cat("Summary of outputs:\n")
cat("- Number of total candidate loci:", nrow(locus_stats_df), "\n")
cat("- Number of usable loci:", nrow(usable_loci), "\n")
cat("- Number of best loci:", nrow(best_loci), "\n")
cat("- Generated", k+1, "iterations of sampled loci\n")
cat("- Selected", length(unique(genes_k)), "unique genes for targeting\n")