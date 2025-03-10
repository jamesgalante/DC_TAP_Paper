## Reformat CRISPR enhancer screen pipeline output to ENCODE format

# Saving image for debugging
save.image(paste0("RDA_objects/create_encode_dataset_", snakemake@wildcards$sample, ".rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(rtracklayer)
})

# load input data ----------------------------------------------------------------------------------

# column types main results input file
results_cols <- cols(
  .default = col_double(),  # for power columns, which can vary in names
  perturbation = col_character(),
  gene = col_character(),
  logFC = col_double(),
  ci_high = col_double(),
  ci_low = col_double(),
  pvalue = col_double(),
  pval_adj = col_double(),
  pert_chr = col_character(),
  pert_start = col_integer(),
  pert_end = col_integer(),
  gene_chr = col_character(),
  gene_tss = col_integer(),
  gene_strand = col_character(),
  dist_to_tss = col_double(),
  pert_level = col_character(),
  target_type = col_character(),
  cells = col_double(),
  avg_expr = col_double(),
  disp_outlier_deseq2 = col_logical()
)

# load differential expression results
results <- read_tsv(snakemake@input$results, col_types = results_cols, progress = FALSE)

# load genome annotations
annot <- import(snakemake@input$annot)

# filter out any transcripts that should be ignored
annot <- annot[!annot$transcript_id %in% snakemake@params$ignore_txs]

# remove any version numbers from gene ids
annot$gene_id <- sub("\\..+$", "", annot$gene_id)

# column types in guide targets file
guide_targets_cols <- cols(
  chr = col_character(),
  start = col_double(),
  end = col_double(),
  name = col_character(),
  strand = col_character(),
  spacer = col_character(),
  target_chr = col_character(),
  target_start = col_double(),
  target_end = col_double(),
  target_name = col_character(),
  target_strand = col_character(),
  target_type = col_character()
)

# load guide targets file
guide_targets <- read_tsv(snakemake@input$guide_targets, col_types = guide_targets_cols,
                          progress = FALSE)

# add gene names respectively ensembl ids to results -----------------------------------------------

# create gene name - ensembl id table
gene_ids <- mcols(annot)[, c("gene_id", "gene_name")] %>% 
  as.data.frame() %>% 
  distinct() %>% 
  dplyr::rename(measuredEnsemblID = gene_id, measuredGeneSymbol = gene_name)

# get number of pairs in results to later check whether gene name - ensembl id matching was 1-to-1
results_pairs <- nrow(results)

# add missing gene id to data and rename columns according to ENCODE format
if (snakemake@params$gene_ids == "ensembl") {
  
  results <- results %>% 
    dplyr::rename(measuredEnsemblID = gene) %>% 
    left_join(gene_ids, "measuredEnsemblID")
  
} else if (snakemake@params$gene_ids == "symbol") {
  
  results <- results %>% 
    dplyr::rename(measuredGeneSymbol = gene) %>% 
    left_join(gene_ids, "measuredGeneSymbol")
  
} else {
  stop("Invalid 'gene_ids' snakemake parameter.", call. = FALSE)
}

# verify that gene name - ensembl id matching was 1-to-1
if (nrow(results) != results_pairs) {
  stop("Gene name - Ensembl ID matching not 1-to-1!", call. = FALSE)
}

# identify and tag promoter targeting guides, and intronic and promoter enhancers ------------------

# create GRangesList with all annotated exons per annotated gene and transcripts
exons <- annot[annot$type == "exon"]
genes <- split(exons, f = exons$gene_id)
txs <- split(exons, f = exons$transcript_id)

# get gene locus coordinates
gene_loci <- unlist(range(genes))

# get all promoters for annotated transcripts
promoters <- unlist(promoters(range(txs), upstream = snakemake@params$tss_min_dist))

# get gene id for every transcript id
gene_ids_per_tx_id <- mcols(exons)[, c("transcript_id", "gene_id")] %>% 
  as.data.frame() %>% 
  distinct() %>% 
  deframe()

# add gene id to promoters metadata columns
mcols(promoters) <- DataFrame(feature_id = gene_ids_per_tx_id[names(promoters)], type = "promoter")

# add metadata columns to exons and gene loci
mcols(exons) <- DataFrame(feature_id = exons$gene_id, type = "exon")
mcols(gene_loci) <- DataFrame(feature_id = names(gene_loci), type = "gene")

# combine exons, genes and promoters into one feature GRanges object
features <- c(exons, gene_loci, promoters)
names(features) <- NULL

# extract perturbation coordinates
pert_coords <- results %>% 
  select(pert_chr, pert_start, pert_end, perturbation) %>% 
  distinct() %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)

# overlap perturbations with features
feature_ovl <- findOverlapPairs(pert_coords, features, ignore.strand = TRUE)

# create table with overlapping features per perturbation in long format
feature_ovl <- tibble(perturbation = first(feature_ovl)$perturbation,
                      feature_id = second(feature_ovl)$feature_id,
                      type = second(feature_ovl)$type)

# add tested genes for perturbations overlapping any features
feature_ovl <- feature_ovl %>% 
  distinct() %>% 
  left_join(select(results, perturbation, measuredEnsemblID), by = "perturbation")

# Add gene symbols to feature_ovl (solely for adding the gene symbols to the promoter overlap line) - added by James
feature_ovl <- feature_ovl %>% 
  left_join(gene_ids, by = c("feature_id" = "measuredEnsemblID"))

# summarize feature overlaps
feature_ovl_summary <- feature_ovl %>% 
  group_by(perturbation, measuredEnsemblID, type) %>% 
  summarize(ovl = n() > 0,
            ovl_target = unique(measuredEnsemblID) %in% feature_id,
            promoter_gene = if (type[1] == "promoter") paste(unique(measuredGeneSymbol), collapse = "|") else NA_character_, # added by James
            .groups = "drop")

# reformat overlap summary
flag_summary <- feature_ovl_summary %>% 
  select(-promoter_gene) %>%
  pivot_longer(cols = c(ovl, ovl_target), names_to = "ovl_type", values_to = "ovl") %>% 
  unite(type, ovl_type, type, sep = "_") %>% 
  pivot_wider(names_from = type, values_from = ovl, values_fill = FALSE)

# Separately, aggregate promoter_gene values for promoter rows - added by James
promoter_summary <- feature_ovl_summary %>%
  filter(type == "promoter") %>%
  group_by(perturbation, measuredEnsemblID) %>%
  summarize(promoter_gene = paste(unique(promoter_gene[!is.na(promoter_gene)]), collapse = ";"),
            .groups = "drop")

# Join the flag summary and promoter summary so each key appears only once - added by James
feature_ovl_summary_final <- left_join(flag_summary, promoter_summary, by = c("perturbation", "measuredEnsemblID"))

# infer whether perturbation overlaps intron
feature_ovl_summary_final <- feature_ovl_summary_final %>% 
  mutate(ovl_intron = if_else((ovl_gene & !ovl_exon), true = TRUE, false = FALSE),
         ovl_target_intron = if_else((ovl_intron & ovl_target_gene), true = TRUE, false = FALSE))

# add selected overlaps to results
results_tmp <- feature_ovl_summary_final %>% 
  select(perturbation, measuredEnsemblID, ovl_target_exon, ovl_target_intron, ovl_promoter, promoter_gene) %>% 
  left_join(results, ., by = c("perturbation", "measuredEnsemblID"))

# create ValidConnection column based on overlaps
results <- results_tmp %>%
  rowwise() %>%
  mutate(ValidConnection = {
    conn <- paste(c(
      if (isTRUE(ovl_target_exon)) "overlaps target gene exon" else NULL,
      if (isTRUE(ovl_target_intron)) "overlaps target gene intron" else NULL,
      if (isTRUE(ovl_promoter)) paste("overlaps potential promoter:", promoter_gene) else NULL
    ), collapse = "; ")
    if (nchar(conn) == 0) "TRUE" else conn
  }) %>%
  ungroup()

# add TSS targeting flag to known TSS controls (if TSS control pattern is provided)
if (!is.null(snakemake@params$tss_ctrl_tag)) {
  results <- results %>% 
    mutate(ValidConnection = if_else(target_type %in% snakemake@params$tss_ctrl_tag,
                                     true = paste(ValidConnection, "TSS targeting guide(s)", sep="; "), 
                                     false = ValidConnection))
}

# add gRNA sequences -------------------------------------------------------------------------------

# get whether perturbation id refers to guide or target name
if (snakemake@wildcards$strategy == "perGRNA") {
  pert_id <- "name"
} else {
  pert_id <- "target_name"
}

# get all gRNA sequences per perturbation and concatenate to one string
grna_seqs <- guide_targets %>% 
  select(perturbation = all_of(pert_id), spacer) %>% 
  group_by(perturbation) %>% 
  summarize(guideSeq = paste(spacer, collapse = ";"))

# add to data
results <- left_join(results, grna_seqs, by = "perturbation")

# add NA for 'guideSpacerSeq'
results <- mutate(results, guideSpacerSeq = NA_character_)

# create ENCODE style output -----------------------------------------------------------------------

# add significance column and gene TSS coordinates
output <- results %>% 
  mutate(Significant = pval_adj < snakemake@params$padj_threshold, 
         startTSS = gene_tss, endTSS = gene_tss + 1, strandGene = gene_strand)

# add additional meta data columns in ENCODE style to create output
output <- output %>% 
  mutate(strandPerturbationTarget = ".",
         PerturbationTargetID = paste0(pert_chr, ":", pert_start, "-", pert_end, ":",
                                       strandPerturbationTarget),
         name = paste(measuredGeneSymbol, PerturbationTargetID, sep = "|"),
         Notes = NA_character_,
         Reference = snakemake@params$reference)

# rename columns containing simulated power
colnames(output) <- sub("^power_effect_size_", "PowerAtEffectSize", colnames(output))

# reformat to ENCODE style (last 3 columns are not ENCODE format and need to be stripped for upload)
output <- output %>% 
  select(chrom = pert_chr, chromStart = pert_start, chromEnd = pert_end, name, EffectSize = logFC,
         strandPerturbationTarget, PerturbationTargetID, chrTSS = gene_chr, startTSS,
         endTSS, strandGene, EffectSize95ConfidenceIntervalLow = ci_low,
         EffectSize95ConfidenceIntervalHigh = ci_high, measuredGeneSymbol, measuredEnsemblID,
         guideSpacerSeq, guideSeq, Significant, pValue = pvalue, pValueAdjusted = pval_adj,
         starts_with("PowerAtEffectSize"), ValidConnection, Notes, Reference,
         distToTSS = dist_to_tss, avgGeneExpr = avg_expr, nPertCells = cells)

# make sure output is sorted according to genomic coordinates
output <- arrange(output, chrom, chromStart, chromEnd, measuredGeneSymbol)

# write to .tsv file
write_tsv(output, file = snakemake@output[[1]])
