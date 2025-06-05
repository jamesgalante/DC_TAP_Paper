
# This script is for creating the main figure 3 panels

# Effect sizes of significant hits in different perturb-seq screens
rule compare_effect_sizes_to_other_screens:
  input:
    gasperini_results = "results/main_figure_1_and_2/duplicate_pairs_analysis/results_with_element_gene_pair_categories.tsv", 
    dc_tap_results = "results/formatted_dc_tap_results/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements.tsv",
    klann = "resources/main_figure_3/ENCODE_Klann_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz",
    morrisv1 = "resources/main_figure_3/ENCODE_Morrisv1_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz",
    morrisv2 = "resources/main_figure_3/ENCODE_Morrisv2_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz",
    xie = "resources/main_figure_3/ENCODE_Xie_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz"
  output:
    effect_size_comparison = "results/main_figure_3/compare_effect_sizes_to_other_screens/effect_size_comparison.pdf"
  log: "results/main_figure_3/logs/compare_effect_sizes_to_other_screens.log"
  conda:
    "../envs/all_packages.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/main_figure_3/compare_effect_sizes_to_other_screens.R"
    
# Comparing power for all the datasets
rule power_to_detect_change_for_all_datasets:
  input:
    gasperini_results = "results/main_figure_1_and_2/duplicate_pairs_analysis/results_with_element_gene_pair_categories.tsv", 
    dc_tap_results = "results/formatted_dc_tap_results/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements.tsv",
    klann = "resources/main_figure_3/ENCODE_Klann_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz",
    morrisv1 = "resources/main_figure_3/ENCODE_Morrisv1_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz",
    morrisv2 = "resources/main_figure_3/ENCODE_Morrisv2_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz",
    xie = "resources/main_figure_3/ENCODE_Xie_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz"
  output:
    power_for_all_datasets_by_proportion = "results/main_figure_3/power_to_detect_change_for_all_datasets/power_for_all_datasets_by_proportion.pdf",
    power_for_all_datasets_by_count = "results/main_figure_3/power_to_detect_change_for_all_datasets/power_for_all_datasets_by_count.pdf"
  log: "results/main_figure_3/logs/power_to_detect_change_for_all_datasets.log"
  conda:
    "../envs/all_packages.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/main_figure_3/power_to_detect_change_for_all_datasets.R"
