# Script: igv_file_creation.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/igv_file_creation.rda"))
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
  library(rtracklayer)
  library(GenomicRanges)
  library(R.utils)
})

message("Loading input files")
results_w_categories <- read_tsv(snakemake@input$combined_joined_w_categories)

# guide targets - removing negative_control and safe_targeting
k562_guide_targets <- read_tsv(snakemake@input$k562_guide_targets) %>% filter(!target_type %in% c("negative_control", "safe_targeting"))
wtc11_guide_targets <- read_tsv(snakemake@input$wtc11_guide_targets) %>% filter(!target_type %in% c("negative_control", "safe_targeting"))

# Split up WTC11 and K562
wtc11 <- results_w_categories %>% filter(ExperimentCellType == "WTC11")
k562 <- results_w_categories %>% filter(ExperimentCellType == "K562")

# Load in the chain file for liftOver
chain_hg19_to_hg38 <- import.chain(snakemake@input$hg19ToHg38_chain_file)

# Load in the annotation file
annot <- import(snakemake@input$annot)

# Load in the TAP gene type files
k562_gene_table <- read_tsv(snakemake@input$k562_gene_table)
wtc11_gene_table <- read_tsv(snakemake@input$wtc11_gene_table)


### LIFTOVER GUIDE_TARGETS GUIDES =============================================

# --- Process K562 guide targets ---
# Create a GRanges object from the K562 guide_targets (assumed to be in hg19)
k562_guides_gr <- GRanges(
  seqnames = k562_guide_targets$chr,
  ranges = IRanges(start = k562_guide_targets$start, end = k562_guide_targets$end),
  strand = k562_guide_targets$strand
)

# LiftOver to hg38
k562_guides_hg38 <- liftOver(k562_guides_gr, chain_hg19_to_hg38)

# For each guide, if there are multiple mappings, pick the first one. - this only happens for one example: K562_Random_Screen_Crop_8994
k562_guides_hg38 <- GRangesList(lapply(k562_guides_hg38, function(gr) {
  if (length(gr) > 0) gr[1] else gr
}))

# unlist the liftOver results
k562_guides_hg38 <- unlist(k562_guides_hg38)

# Add the new coordinates to the original K562 guide targets tibble
k562_guide_targets_hg38 <- k562_guide_targets %>%
  mutate(
    guide_chr_hg38 = as.character(seqnames(k562_guides_hg38)),
    guide_start_hg38 = start(k562_guides_hg38),
    guide_end_hg38 = end(k562_guides_hg38),
    strand_hg38 = as.character(strand(k562_guides_hg38))
  )

# --- Process WTC11 guide targets ---
# Create a GRanges object from the WTC11 guide_targets (assumed to be in hg19)
wtc11_guides_gr <- GRanges(
  seqnames = wtc11_guide_targets$chr,
  ranges = IRanges(start = wtc11_guide_targets$start, end = wtc11_guide_targets$end),
  strand = wtc11_guide_targets$strand
)

# LiftOver to hg38
wtc11_guides_hg38 <- liftOver(wtc11_guides_gr, chain_hg19_to_hg38)
wtc11_guides_hg38 <- unlist(wtc11_guides_hg38)

# Add the new coordinates to the original WTC11 guide targets tibble
wtc11_guide_targets_hg38 <- wtc11_guide_targets %>%
  mutate(
    guide_chr_hg38 = as.character(seqnames(wtc11_guides_hg38)),
    guide_start_hg38 = start(wtc11_guides_hg38),
    guide_end_hg38 = end(wtc11_guides_hg38),
    strand_hg38 = as.character(strand(wtc11_guides_hg38))
  )


### MERGE GUIDE_TARGETS WITH RESULTS ==========================================

# The results already has the hg38 coordinates of each target, so just join these together
k562_guide_targets_hg38_full <- k562_guide_targets_hg38 %>%
  left_join(
    k562 %>% group_by(target_name) %>% dplyr::slice(1) %>% select(chrom, chromStart, chromEnd, target_name), 
    by = "target_name"
  ) %>%
  dplyr::rename(
    target_chr_hg38 = chrom, 
    target_start_hg38 = chromStart, 
    target_end_hg38 = chromEnd
  )

wtc11_guide_targets_hg38_full <- wtc11_guide_targets_hg38 %>%
  left_join(
    wtc11 %>% group_by(target_name) %>% dplyr::slice(1) %>% select(chrom, chromStart, chromEnd, target_name), 
    by = "target_name"
  ) %>%
  dplyr::rename(
    target_chr_hg38 = chrom, 
    target_start_hg38 = chromStart, 
    target_end_hg38 = chromEnd
  )


### DEFINE COLOR FUNCTIONS ====================================================

# Define color mapping for target types (modify as needed)
target_type_colors <- list(
  "tss_pos" = "#0d95a1",    # Teal
  "tss_random" = "#bb0e3d",  # Red
  "enh"     = "#2757ea",    # Blue
  "DE"      = "#ff9999"     # Light pink
)
get_color <- function(target_type) {
  if (target_type %in% names(target_type_colors)) {
    return(target_type_colors[[target_type]])
  } else {
    return("#000000") # Default black
  }
}
# Convert hex color (e.g., "#0d95a1") to "R,G,B" string
hex_to_rgb <- function(hex_color) {
  hex_color <- gsub("^#", "", hex_color)
  r <- strtoi(substr(hex_color, 1, 2), base = 16)
  g <- strtoi(substr(hex_color, 3, 4), base = 16)
  b <- strtoi(substr(hex_color, 5, 6), base = 16)
  paste(r, g, b, sep = ",")
}

# For each target-gene pair (indicated by an arc in bedpe file) 
  # - we want a track (arc) colored by (in order of precedence)
    # "selfPromoter", "Green"
    # "DistalPromoter_Gene", "Red"
    # "Positive_Control_DistalElement_Gene", "Pink"
    # "Random_Validation_DistalElement_Gene", - random element that make it past encode filters "Blue"
    # "Random_DistalElement_Gene" - random elements that don't make it past encode filters "Orange"
    # "DistalElement_Gene" - all elements that pass ENCODE filters (not necessarily random) "Black"
    # All other pairs tested - I don't think there should be any here "Purple"

# For BEDPE arcs, assign color by category precedence:
# selfPromoter (Green), DistalPromoter_Gene (Red), Positive_Control_DistalElement_Gene (Pink),
# Random_Validation_DistalElement_Gene (Blue), Random_DistalElement_Gene (Orange),
# DistalElement_Gene (Black), else Purple.
get_arc_color <- function(row) {
  if (isTRUE(row$selfPromoter)) {
    return(hex_to_rgb("#00FF00"))  # Green
  } else if (isTRUE(row$DistalPromoter_Gene)) {
    return(hex_to_rgb("#FF0000"))  # Red
  } else if (isTRUE(row$Positive_Control_DistalElement_Gene)) {
    return(hex_to_rgb("#FFC0CB"))  # Pink
  } else if (isTRUE(row$Random_Validation_DistalElement_Gene)) {
    return(hex_to_rgb("#0000FF"))  # Blue
  } else if (isTRUE(row$Random_DistalElement_Gene)) {
    return(hex_to_rgb("#FFA500"))  # Orange
  } else if (isTRUE(row$DistalElement_Gene)) {
    return(hex_to_rgb("#000000"))  # Black
  } else {
    return(hex_to_rgb("#800080"))  # Purple
  }
}

### CREATE FILES ==============================================================

# Helper functions to create BED dataframes
create_guides_bed <- function(df) {
  df %>%
    filter(!is.na(guide_start_hg38), !is.na(guide_end_hg38)) %>%
    mutate(
      chr = guide_chr_hg38,
      start_0based = guide_start_hg38 - 1,
      end = guide_end_hg38,
      name = name,
      score = 1000,
      strand = strand_hg38,
      thickStart = start_0based,
      thickEnd = end,
      itemRgb = map_chr(target_type, ~ hex_to_rgb(get_color(.x)))
    ) %>%
    select(chr, start_0based, end, name, score, strand, thickStart, thickEnd, itemRgb) %>%
    arrange(chr, start_0based)
}

create_targets_bed <- function(df) {
  df %>%
    select(target_chr_hg38, target_start_hg38, target_end_hg38, target_name, target_type) %>%
    distinct() %>%
    filter(!is.na(target_chr_hg38), !is.na(target_start_hg38), !is.na(target_end_hg38)) %>%
    mutate(
      chr = target_chr_hg38,
      start_0based = target_start_hg38 - 1,
      end = target_end_hg38,
      name = target_name,
      score = 1000,
      strand = "+",
      thickStart = start_0based,
      thickEnd = end,
      itemRgb = map_chr(target_type, ~ hex_to_rgb(get_color(.x)))
    ) %>%
    select(chr, start_0based, end, name, score, strand, thickStart, thickEnd, itemRgb) %>%
    arrange(chr, start_0based)
}

# Split merged guide targets by cell type
k562_all <- k562_guide_targets_hg38_full
wtc11_all <- wtc11_guide_targets_hg38_full

# Create BED data for K562
k562_guides_bed <- create_guides_bed(k562_all)
k562_targets_bed <- create_targets_bed(k562_all)

# Create BED data for WTC11
wtc11_guides_bed <- create_guides_bed(wtc11_all)
wtc11_targets_bed <- create_targets_bed(wtc11_all)


### CREATE ARCS BEDPE FILES FROM RESULTS ======================================
# We now create arcs using the original results table (which has the categories).
# Each arc is drawn between the hg38 target region and the TSS of the gene it is paired with.
# Note: In the results table, "chrom", "chromStart", "chromEnd" refer to the target region;
# "chrTSS", "startTSS", "endTSS" refer to the TSS coordinates.
create_arcs_bedpe_results <- function(df) {
  df %>%
    filter(!is.na(chrom), !is.na(chromStart), !is.na(chromEnd),
           !is.na(chrTSS), !is.na(startTSS), !is.na(endTSS)) %>%
    mutate(
      chrom1 = chrom,
      start1 = chromStart - 1,
      end1 = chromEnd,
      chrom2 = chrTSS,
      start2 = startTSS - 1,
      end2 = endTSS,
      name = name,
      score = 1000,
      strand1 = "+",
      strand2 = "+"
    ) %>%
    rowwise() %>%
    mutate(arcColor = get_arc_color(pick(everything()))) %>%
    ungroup() %>%
    select(chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2, arcColor)
}

# There are a few rows where the chrTSS, startTSS, and endTSS aren't filled in
# > k562[which(is.na(k562$chrTSS)),] %>% pull(measuredGeneSymbol) %>% unique()
# [1] "KHDC4"   "ZNRD2"   "ANTKMT"  "JPT2"    "ATP5MC1" "GATD3A"  "RBIS"   
# > wtc11[which(is.na(wtc11$chrTSS)),] %>% pull(measuredGeneSymbol) %>% unique()
# [1] "TENT5B" "FDX2"
# Let's just add those results manually


# UPDATE MISSING TSS COORDINATES ==============================================

# Identify missing genes from both datasets
missing_k562 <- k562 %>% filter(is.na(chrTSS)) %>% pull(measuredGeneSymbol) %>% unique()
missing_wtc11 <- wtc11 %>% filter(is.na(chrTSS)) %>% pull(measuredGeneSymbol) %>% unique()
missing_genes <- unique(c(missing_k562, missing_wtc11))

# Generate TSS annotations using promoters() which is strand-aware:
genes <- annot[annot$type == "gene"]
genes_reduced <- unlist(reduce(split(genes, f = genes$gene_name)))
# Use promoters() to get a 2-bp region at the TSS (5' end) regardless of strand.
tss <- promoters(genes_reduced, upstream = 0, downstream = 2)

# Convert TSS GRanges to a data frame
tss_df <- data.frame(
  chrTSS = as.character(seqnames(tss)),
  startTSS = start(tss),
  endTSS = end(tss),
  measuredGeneSymbol = names(tss),
  stringsAsFactors = FALSE
)

# Filter the TSS data for missing genes only
annot_missing <- tss_df %>% filter(measuredGeneSymbol %in% missing_genes)

# Update k562: join and replace only missing TSS values
k562 <- k562 %>% 
  left_join(annot_missing, by = "measuredGeneSymbol", suffix = c("", ".new")) %>%
  mutate(
    chrTSS   = if_else(is.na(chrTSS), chrTSS.new, chrTSS),
    startTSS = if_else(is.na(startTSS), startTSS.new, startTSS),
    endTSS   = if_else(is.na(endTSS), endTSS.new, endTSS)
  ) %>%
  select(-chrTSS.new, -startTSS.new, -endTSS.new)

# Update wtc11 similarly
wtc11 <- wtc11 %>% 
  left_join(annot_missing, by = "measuredGeneSymbol", suffix = c("", ".new")) %>%
  mutate(
    chrTSS   = if_else(is.na(chrTSS), chrTSS.new, chrTSS),
    startTSS = if_else(is.na(startTSS), startTSS.new, startTSS),
    endTSS   = if_else(is.na(endTSS), endTSS.new, endTSS)
  ) %>%
  select(-chrTSS.new, -startTSS.new, -endTSS.new)

# Create arcs for K562 and WTC11 using the results table
k562_arcs_bedpe <- create_arcs_bedpe_results(k562)
wtc11_arcs_bedpe <- create_arcs_bedpe_results(wtc11)


### CREATING TAP GENE TRACKS ==================================================

# We want a track of all genes that were read out with TAP Seq
k562_genes <- k562_gene_table$gene
wtc11_genes <- wtc11_gene_table$gene

# Filter the TSS data for K562 genes only
k562_genes_tss <- tss_df %>% filter(measuredGeneSymbol %in% k562_genes)
wtc11_genes_tss <- tss_df %>% filter(measuredGeneSymbol %in% wtc11_genes)

# There are a couple genes that aren't included
# > k562_gene_table %>% filter(!gene %in% k562_genes_tss$measuredGeneSymbol) %>% pull(gene)
# [1] "FAM173A" "SSSCA1"  "C8orf59"
# > wtc11_gene_table %>% filter(!gene %in% wtc11_genes_tss$measuredGeneSymbol) %>% pull(gene)
# [1] "FAM46B" "FAM64A"

# Perhaps these genes can be found by another name
# K562: "FAM173A" = "ANTKMT", "SSSCA1" = "ZNRD2, "C8orf59" = "RBIS"
# WTC11: "FAM46B" = "TENT5B", "FAM64A" = "PIMREG"

# Define rename mappings for missing genes
k562_rename <- c("FAM173A" = "ANTKMT", "SSSCA1" = "ZNRD2", "C8orf59" = "RBIS")
wtc11_rename <- c("FAM46B" = "TENT5B", "FAM64A" = "PIMREG")

# Update gene names in gene tables
k562_gene_table_updated <- k562_gene_table %>%
  mutate(gene = recode(gene, !!!k562_rename))
wtc11_gene_table_updated <- wtc11_gene_table %>%
  mutate(gene = recode(gene, !!!wtc11_rename))

# Use the updated gene lists
k562_genes_updated <- k562_gene_table_updated$gene
wtc11_genes_updated <- wtc11_gene_table_updated$gene

# Filter TSS data using the updated gene names
k562_genes_tss <- tss_df %>% filter(measuredGeneSymbol %in% k562_genes_updated)
wtc11_genes_tss <- tss_df %>% filter(measuredGeneSymbol %in% wtc11_genes_updated)

# Convert the TSS data to BED format (0-based start) for track files
k562_genes_tss_bed <- k562_genes_tss %>%
  mutate(
    chrom = chrTSS,
    start = startTSS - 1,  # convert to 0-based
    end = endTSS,
    name = measuredGeneSymbol,
    score = 1000,
    strand = "+"
  ) %>%
  select(chrom, start, end, name, score, strand) %>%
  arrange(chrom, start)

wtc11_genes_tss_bed <- wtc11_genes_tss %>%
  mutate(
    chrom = chrTSS,
    start = startTSS - 1,
    end = endTSS,
    name = measuredGeneSymbol,
    score = 1000,
    strand = "+"
  ) %>%
  select(chrom, start, end, name, score, strand) %>%
  arrange(chrom, start)


### SIGNIFICANCE ARCS =========================================================

# Create a function to generate BEDPE for significant arcs only
create_significance_arcs_bedpe <- function(df) {
  df %>%
    filter(Significant == TRUE,
           !is.na(chrom), !is.na(chromStart), !is.na(chromEnd),
           !is.na(chrTSS), !is.na(startTSS), !is.na(endTSS)) %>%
    mutate(
      chrom1 = chrom,
      start1 = chromStart - 1,
      end1 = chromEnd,
      chrom2 = chrTSS,
      start2 = startTSS - 1,
      end2 = endTSS,
      name = name,
      score = 1000,
      strand1 = "+",
      strand2 = "+"
    ) %>%
    rowwise() %>%
    mutate(arcColor = if_else(EffectSize < 0,
                              hex_to_rgb("#FF0000"),  # red for downregulated
                              hex_to_rgb("#00FF00"))) %>%  # green for upregulated
    ungroup() %>%
    select(chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2, arcColor)
}

# Generate significance-only arcs for K562 and WTC11
k562_signif_arcs_bedpe <- create_significance_arcs_bedpe(k562)
wtc11_signif_arcs_bedpe <- create_significance_arcs_bedpe(wtc11)


### SAVE OUTPUT FILES =========================================================

message("Saving Output")

# For K562
write_lines('track name="K562 Guide Targets (hg38)" description="Guide Target Regions (hg38) - K562" visibility=2 itemRgb="On"', snakemake@output$k562_guides)
write_tsv(k562_guides_bed, snakemake@output$k562_guides, append = TRUE, col_names = FALSE)

write_lines('track name="K562 Target Regions (hg38)" description="Unique Target Regions (hg38) - K562" visibility=2 itemRgb="On"', snakemake@output$k562_targets)
write_tsv(k562_targets_bed, snakemake@output$k562_targets, append = TRUE, col_names = FALSE)

cat('chr1\tx1\tx2\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tcolor\n', file = snakemake@output$k562_arcs)
write_tsv(k562_arcs_bedpe, snakemake@output$k562_arcs, append = TRUE, col_names = FALSE)

# For WTC11
write_lines('track name="WTC11 Guide Targets (hg38)" description="Guide Target Regions (hg38) - WTC11" visibility=2 itemRgb="On"', snakemake@output$wtc11_guides)
write_tsv(wtc11_guides_bed, snakemake@output$wtc11_guides, append = TRUE, col_names = FALSE)

write_lines('track name="WTC11 Target Regions (hg38)" description="Unique Target Regions (hg38) - WTC11" visibility=2 itemRgb="On"', snakemake@output$wtc11_targets)
write_tsv(wtc11_targets_bed, snakemake@output$wtc11_targets, append = TRUE, col_names = FALSE)

cat('chr1\tx1\tx2\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tcolor\n', file = snakemake@output$wtc11_arcs)
write_tsv(wtc11_arcs_bedpe, snakemake@output$wtc11_arcs, append = TRUE, col_names = FALSE)

# Save combined color maps (from both cell types)
all_guides_bed <- bind_rows(create_guides_bed(k562_guide_targets_hg38_full), create_guides_bed(wtc11_guide_targets_hg38_full))
all_targets_bed <- bind_rows(create_targets_bed(k562_guide_targets_hg38_full), create_targets_bed(wtc11_guide_targets_hg38_full))
guide_color_map <- all_guides_bed %>% select(name, itemRgb) %>% distinct()
target_color_map <- all_targets_bed %>% select(name, itemRgb) %>% distinct()

write_tsv(guide_color_map, snakemake@output$guide_color_map)
write_tsv(target_color_map, snakemake@output$target_color_map)

# Save the TSS track files
# For K562:
track_header_k562 <- 'track name="K562 TAP Seq Genes" description="TSS regions for TAP Seq genes in K562" visibility=2 itemRgb="On"'
write_lines(track_header_k562, snakemake@output$k562_TAP_seq_genes)
write_tsv(k562_genes_tss_bed, snakemake@output$k562_TAP_seq_genes, append = TRUE, col_names = FALSE)

# For WTC11:
track_header_wtc11 <- 'track name="WTC11 TAP Seq Genes" description="TSS regions for TAP Seq genes in WTC11" visibility=2 itemRgb="On"'
write_lines(track_header_wtc11, snakemake@output$wtc11_TAP_seq_genes)
write_tsv(wtc11_genes_tss_bed, snakemake@output$wtc11_TAP_seq_genes, append = TRUE, col_names = FALSE)

# Save the significance arcs using Snakemake output objects
cat('chr1\tx1\tx2\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tcolor\n', file = snakemake@output$k562_signif_arcs)
write_tsv(k562_signif_arcs_bedpe, snakemake@output$k562_signif_arcs, append = TRUE, col_names = FALSE)

cat('chr1\tx1\tx2\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tcolor\n', file = snakemake@output$wtc11_signif_arcs)
write_tsv(wtc11_signif_arcs_bedpe, snakemake@output$wtc11_signif_arcs, append = TRUE, col_names = FALSE)



### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)