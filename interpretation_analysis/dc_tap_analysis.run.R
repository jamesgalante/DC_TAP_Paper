suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggdist)
  library(cowplot)
})

source("dc_tap_analysis.functions.R")


### --- DEFINE FILE PATHS AND PARAMS ---

# download from: https://github.com/EngreitzLab/ENCODE_rE2G/blob/dev/resources/external_features/gene_promoter_class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.tsv
promoter_class_path <- "ENCODE_rE2G/resources/external_features/gene_promoter_class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.tsv"

# download from: https://github.com/EngreitzLab/ENCODE_rE2G/blob/dev/reference/CollapsedGeneBounds.hg38.TSS500bp.bed"
abc_genes_path <- "ENCODE_rE2G/reference/CollapsedGeneBounds.hg38.TSS500bp.bed"

## inputs for combining all annotations
# parameters for categorization
remove_promoters <- TRUE
quantiles <- c(0.1, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.75, 0.85, 0.9, 0.95)
cell_types <- c("K562", "GM12878", "HCT116", "Jurkat", "WTC11")

# results from chrom-annotate pipeline
chrom_annotate_results <- "chrom_annotate_results"
enh_lists <- c(file.path(chrom_annotate_results, cell_types, "EnhancerList.extended.tsv"))
names(enh_lists) <- cell_types

crispr_annot_files <- c(other_crispr = file.path(chrom_annotate_results, "CRISPR_data", "other_crispr_screens.tsv"),
    dc_tap_gasperini = file.path(chrom_annotate_results, "CRISPR_data", "dc_tap_gasperini.random_distal_element_pairs.tsv"))

# results from direct effect rate annotations
direct_effect_results <- "direct_effect_results"
crispr_annot_files <- c(other_crispr = file.path(direct_effect_results, "other_crispr_screens.direct_effect.tsv.gz")
    dc_tap_gasperini = file.path(direct_effect_results, "dc_tap_gasperini.random_distal_element_pairs.direct_effect.tsv.gz"))

# results from encode_re2g pipeline to get all candidate element-gene pairs
e2g_base <- "ENCODE_rE2G_results"
e2g_res <- c(GM12878 = "2025_0227_validation_new_inputs", HCT116 = "2025_0227_validation_new_inputs", K562 = "2025_0227_validation_new_inputs",
    Jurkat = "2025_0227_validation_new_inputs", WTC11 = "2025_0227_validation_new_inputs")
e2g_files <- lapply(cell_types, function(ct) file.path(e2g_base, e2g_res[ct], paste0(ct, "_H3K27ac_megamap"),
                                                "dhs_h3k27ac_megamap", "encode_e2g_predictions.tsv.gz")) %>% unlist() %>% setNames(cell_types)

## categorized pairs file paths 
recat_dir <- file.path("categorized_pairs")
crispr_encode_categorized_file <- file.path(recat_dir, "other_crispr_screens.categorized.tsv.gz")
crispr_sceptre_categorized_file <- file.path(recat_dir, "dc_tap_gasperini.random_distal_element_pairs.categorized.tsv.gz")
thresholds_file <- file.path(recat_dir, "thresholds.tsv")

## model benchmarking results from CRISPR_comparison pipeline
merged_dc_tap_benchmarking_file <- file.path("CRISPR_benchmarks/results/DC_TAP/expt_pred_merged_annot.txt.gz")
e2g_threshold <- 0.201

### --- ANNOTATE DATA WITH CHROMATIN CATEGORIES ---
enh <- lapply(enh_lists, read_enh_lists, remove_promoters) %>% 
    rbindlist(idcol = "cell_type") %>% as.data.frame() %>%
    mutate(CTCF.H3K27ac.ratio = (CTCF.RPM) / (H3K27ac.RPM + 0.001)) %>% 
    replace(is.na(.), 0)
message("Read in candidate elements!")

thresholds <- get_category_thresholds(enh, quantiles) %>% arrange(cell_type, feature)
fwrite(thresholds, thresholds_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
message("Saved thresholds!")

dc_tap_gasperini <- read_sceptre_crispr_data(crispr_annot_files[["dc_tap_gasperini"]], direct_effect_files[["dc_tap_gasperini"]])
other_crispr <- read_sceptre_crispr_data(crispr_annot_files[["other_crispr"]], direct_effect_files[["other_crispr"]])
message("Read in and formatted CRISPR data!")

ubiq_expr_genes <- fread(promoter_class_path, sep = "\t") %>% filter(is_ubiquitous_uniform %in% c("True", TRUE)) %>% pull(TargetGene)
crispr_sceptre <- dc_tap_gasperini %>%
        categorize_elements(thresholds, H3K27ac_q_high = 0.9, H3K27ac_q_low = 0.5) %>% 
        mutate(ubiq_category = ifelse(measuredGeneSymbol %in% ubiq_expr_genes, "Ubiq. expr. gene", "Other gene"))
fwrite(crispr_sceptre, crispr_sceptre_categorized_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
other_crispr <- other_crispr %>% 
        categorize_elements(thresholds, H3K27ac_q_high = 0.9, H3K27ac_q_low = 0.5) %>% 
        mutate(ubiq_category = ifelse(measuredGeneSymbol %in% ubiq_expr_genes, "Ubiq. expr. gene", "Other gene"))
fwrite(crispr_encode, crispr_encode_categorized_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
message("Categorized CRISPR data!")

enh_cat <- enh %>% 
    select(chr, start, end, cell_type, H3K27me3_peak_overlap, CTCF_peak_overlap, H3K27ac_peak_overlap, H3K27ac.RPM.expandedRegion) %>% 
    categorize_elements(thresholds, H3K27ac_q_high = 0.9, H3K27ac_q_low = 0.5)

## summarize by category
enh_pairs_cat <- annotate_genomewide_pairs(enh_cat, e2g_files, cell_types, remove_promoters, 2e6)
gw_pairs <- enh_pairs_cat[[1]]
fwrite(gw_pairs, file.path(recat_dir, "all_gw_pairs.categorized.summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
message("Categorized genome-wide pairs!")

### --- PLOT CHROMATIN CATEGORY PROPORTIONS ---
crispr_sceptre_use <- crispr_sceptre %>%
        mutate(data_category = ifelse(Dataset == "Gasperini2019", "training", "validation")) %>% 
        mutate(direct_vs_indirect = ifelse(EffectSize < 0, direct_vs_indirect_negative, direct_vs_indirect_positive))

# DC-TAP by cell type (Fig. S5)
z <- plot_category_proportion_by_dataset_sceptre(gw_pairs, crispr_sceptre_use, "element_category",
        direct_effect_weighted = TRUE, direct_effect_threshold = 0)
    ggsave(file.path(recat_dir, "proportion_by_cell_type.element_category.de_weighted.for_dc_tap.pdf"), z, height = 6, width = 8)

# DC-TAP combined (Fig. 5b)
crispr_dc_combined <- crispr_sceptre %>%
    mutate(Dataset = ifelse(Dataset %in% c("K562_DC_TAP", "WTC11_DC_TAP"), "DC_TAP", Dataset)) %>% 
    mutate(cell_type = ifelse(Dataset == "DC_TAP", "K562_WTC11", cell_type))
gw_k562 <- gw_pairs %>% filter(cell_type == "K562")    
gw_pairs_combined <- gw_pairs %>% mutate(cell_type = ifelse(cell_type %in% c("K562", "WTC11"), "K562_WTC11", cell_type)) %>% 
    rbind(gw_k562)
z <- plot_category_proportion_by_dataset_sceptre(gw_pairs_combined, crispr_dc_combined, "element_category",
    direct_effect_weighted = TRUE, direct_effect_threshold = 0)
ggsave(file.path(recat_dir, "proportion_by_cell_type.element_category.de_weighted.combine_dc_tap.pdf"), z, height = 6, width = 7)


### -- COMPUTE RESULTS BY PROMOTER CLASS (Fig. 5c,d) ---
this_results_dir <- "housekeeping_genes"
crispr_sceptre_use <- crispr_sceptre %>%
    mutate(Dataset = ifelse(Dataset %in% c("K562_DC_TAP", "WTC11_DC_TAP"), "DC_TAP", Dataset)) %>% 
    mutate(EffectSize = EffectSize * 100) %>% 
    mutate(direct_vs_indirect = ifelse(EffectSize < 0, direct_vs_indirect_negative, direct_vs_indirect_positive))
group_var <- "ubiq_category"
for (use_allpower in c(TRUE, FALSE)){
    if (use_allpower) {
        out_dir <- file.path(this_results_dir, "v3_combinedDC_TAP_all_metrics_comparison_allPower"); dir.create(out_dir, showWarnings = FALSE)
        out_dir2 <- file.path(this_results_dir, "v3_combinedDC_TAP_all_metrics_significance_allPower"); dir.create(out_dir2, showWarnings = FALSE)
    } else {
        out_dir <- file.path(this_results_dir, "v3_combinedDC_TAP_all_metrics_comparison_wellPowered"); dir.create(out_dir, showWarnings = FALSE)
        out_dir2 <- file.path(this_results_dir, "v3_combinedDC_TAP_all_metrics_significance_wellPowered"); dir.create(out_dir2, showWarnings = FALSE)
    }

    # direct effect weighted
    out_prefix <- paste0(out_dir, "/weighted_"); out_prefix2 <- paste0(out_dir2, "/weighted_")
    plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
        direct_effect_weighted = TRUE, direct_effect_threshold = NULL)
    compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
        direct_effect_weighted = TRUE, direct_effect_threshold = NULL)

    # direct effect filtered
    out_prefix <- paste0(out_dir, "/filter50_"); out_prefix2 <- paste0(out_dir2, "/filter50_")
    plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
        direct_effect_weighted = FALSE, direct_effect_threshold = 0.5)
    compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
        direct_effect_weighted = FALSE, direct_effect_threshold = 0.5)

    # no adjustments for p(direct)
    out_prefix <- paste0(out_dir, "/"); out_prefix2 <- paste0(out_dir2, "/")
    plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
        direct_effect_weighted = FALSE, direct_effect_threshold = 0)
    compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
        direct_effect_weighted = FALSE, direct_effect_threshold = 0)
}

### --- COMPARE TO ENCODE-rE2G (Fig. 6) ---
crispr_annot <- fread(merged_dc_tap_benchmarking_file) %>% 
        mutate(Dataset = ifelse(Dataset %in% c("K562_DC_TAP", "WTC11_DC_TAP"), "DC_TAP", Dataset)) %>% 
        mutate(CTCF_category = ifelse(element_category == "CTCF element", "CTCF element", "Other element")) %>% 
        filter(Dataset == "DC_TAP")
crispr_annot_direct <- filter(crispr_annot, direct_vs_indirect_negative >= 0.5)

this_out_dir = paste0(this_results_dir, "/score_by_category"); dir.create(this_out_dir, showWarnings = FALSE)
out_prefix <- paste0(this_out_dir, "/")
out_prefix_direct <- paste0(this_out_dir, "/filter50pct_")

plot_score_by_category(crispr_annot_direct, "element_category", "Regulated", e2g_threshold, out_prefix_direct)
plot_score_by_category(crispr_annot_direct, "enhancerness", "Regulated", e2g_threshold, out_prefix_direct)
plot_score_by_category(crispr_annot_direct, "CTCF_category", "Regulated", e2g_threshold, out_prefix_direct)
plot_score_by_category(crispr_annot_direct, "ubiq_category", "Regulated", e2g_threshold, out_prefix_direct)
plot_score_by_category(crispr_annot_direct, "EffectSize_5pct", "Regulated", e2g_threshold, out_prefix_direct)
plot_score_by_category(crispr_annot, "pDirect_90pct", "Regulated", e2g_threshold, out_prefix)


### --- ESTIMATE UNDETECTED VS DETECTED SIGNIFICANT POSITIVES ---
this_out <- file.path(this_results_dir, "detected_positives"); dir.create(this_out, showWarnings = FALSE)
out_prefix <- paste0(this_out, "/5_10_")
es_bins <- c("<=5%", "(5%, 10%]", "(10%, 15%]", "(15%, 25%]", "(25%, 50%]", "(50%, Inf)")
bin_min <- c(0, 5, 10, 15, 25, 50); names(bin_min) <- es_bins
bin_max <- c(5, 10, 15, 25, 50, Inf); names(bin_max) <- es_bins        
plot_positives_by_effect_size(crispr_sceptre, es_bins, bin_min, bin_max, "Regulated", out_prefix)

