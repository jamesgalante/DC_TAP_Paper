
# Rule to get numbers for the paper
rule calculate_summary_statistics_for_screen:
  input:
    combined_validation = "results/benchmark_validation_datasets/crispr_benchmarking/expt_pred_merged_annot/validation_expt_pred_merged_annot.txt.gz",
    guide_targets = expand("results/process_validation_datasets/{sample}/guide_targets.tsv", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    encode_training_dataset = "resources/analyze_validation_datasets/expt_pred_merged_annot/training_expt_pred_merged_annot.txt.gz",
    create_ensemble_encode_input = expand("results/benchmark_validation_datasets/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_{sample}_0.13gStd_Sceptre_perCRE_0.8pwrAt15effect_GRCh38.tsv.gz", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    abc_canonical_tss = "resources/benchmark_validation_datasets/crispr_benchmarking/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed"
  output:
    combined_joined_w_categories = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/combined_joined_w_categories.tsv",
    summarized_categories = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/summarized_categories.tsv",
    igvf_formatted_file = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/igvf_formatted_file.tsv"
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/calculate_summary_statistics_for_screen.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/format_dc_tap_results/calculate_summary_statistics_for_screen.R"
    
# Rule to deal with specific pairs
rule modify_specific_pairs_in_final_file:
  input:
    summarized_categories = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/summarized_categories.tsv"
  output:
    Formatted_DC_TAP_Seq_Results = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/Formatted_DC_TAP_Seq_Results.tsv",
    summary_K562 = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/summary_K562.tsv",
    summary_WTC11 = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/summary_WTC11.tsv"
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/modify_specific_pairs_in_final_file.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/format_dc_tap_results/modify_specific_pairs_in_final_file.R"
    
# Resize DC TAP elements for epigenetic categories overlap
rule resize_dc_tap_elements:
  input:
    Formatted_DC_TAP_Seq_Results = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/Formatted_DC_TAP_Seq_Results.tsv"
  output:
    resized_Formatted_DC_TAP_Seq_Results = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/resized_Formatted_DC_TAP_Seq_Results.tsv"
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/resize_dc_tap_elements.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/format_dc_tap_results/resize_dc_tap_elements.R"

# Add Maya's epigenetic categories for each pair
rule add_element_epigenetic_categories:
  input:
    resized_Formatted_DC_TAP_Seq_Results = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/resized_Formatted_DC_TAP_Seq_Results.tsv",
    categorized_data = "resources/analyze_validation_datasets/formatting_dc_tap_results/all_DC_TAP.h3k27me3_quantile_50.ratio_quantile_35.h3k27ac_quantiles_60_90.WTC11_as_WTC11.tsv"
  output:
    Formatted_DC_TAP_Seq_Results_w_Categories = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/Formatted_DC_TAP_Seq_Results_w_Categories.tsv"
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/add_element_epigenetic_categories.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/format_dc_tap_results/add_element_epigenetic_categories.R"
