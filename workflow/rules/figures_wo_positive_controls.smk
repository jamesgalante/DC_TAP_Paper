
# These rules pertain to any figure that was also made with significance called without positive controls included in FDR correction

# Create distance by effect size plots
rule distance_by_effect_size_plots_FDR_wo_pos_controls:
  input:
    screen_results = "results/formatted_dc_tap_results/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements.tsv"
  output:
    k562_distance_by_es = "results/figures_wo_positive_controls/k562_distance_by_es.pdf",
    wtc11_distance_by_es = "results/figures_wo_positive_controls/wtc11_distance_by_es.pdf"
  log: "results/figures_wo_positive_controls/logs/distance_by_effect_size_plots.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/figures_wo_positive_controls/distance_by_effect_size_plots.R"

# Effect Size Plots
rule effect_size_plots_FDR_wo_pos_controls:
  input:
    screen_results = "results/formatted_dc_tap_results/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements.tsv"
  output:
    effect_size_boxplot = "results/figures_wo_positive_controls/effect_size_boxplot.pdf",
    positive_control_self_promoters_effect_size_boxplot = "results/figures_wo_positive_controls/positive_control_self_promoters_effect_size_boxplot.pdf"
  log: "results/figures_wo_positive_controls/logs/effect_size_plots.log"
  conda:
    "../envs/all_packages.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/figures_wo_positive_controls/effect_size_plots.R"

# Effect sizes of significant hits in different perturb-seq screens
rule compare_effect_sizes_to_other_screens_FDR_wo_pos_controls:
  input:
    gasperini_results = "results/main_figure_1_and_2/duplicate_pairs_analysis/results_with_element_gene_pair_categories.tsv",
    dc_tap_results = "results/formatted_dc_tap_results/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements.tsv",
    klann = "resources/main_figure_3/ENCODE_Klann_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz",
    morrisv1 = "resources/main_figure_3/ENCODE_Morrisv1_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz",
    morrisv2 = "resources/main_figure_3/ENCODE_Morrisv2_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz",
    xie = "resources/main_figure_3/ENCODE_Xie_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz"
  output:
    effect_size_comparison = "results/figures_wo_positive_controls/compare_effect_sizes_to_other_screens/effect_size_comparison.pdf"
  log: "results/figures_wo_positive_controls/logs/compare_effect_sizes_to_other_screens.log"
  conda:
    "../envs/all_packages.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/figures_wo_positive_controls/compare_effect_sizes_to_other_screens.R"
    
# Comparing power for all the datasets
rule power_to_detect_change_for_all_datasets_FDR_wo_pos_controls:
  input:
    gasperini_results = "results/main_figure_1_and_2/duplicate_pairs_analysis/results_with_element_gene_pair_categories.tsv",
    dc_tap_results = "results/formatted_dc_tap_results/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements.tsv",
    klann = "resources/main_figure_3/ENCODE_Klann_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz",
    morrisv1 = "resources/main_figure_3/ENCODE_Morrisv1_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz",
    morrisv2 = "resources/main_figure_3/ENCODE_Morrisv2_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz",
    xie = "resources/main_figure_3/ENCODE_Xie_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz"
  output:
    power_for_all_datasets_by_proportion = "results/figures_wo_positive_controls/power_to_detect_change_for_all_datasets/power_for_all_datasets_by_proportion.pdf"
  log: "results/figures_wo_positive_controls/logs/power_to_detect_change_for_all_datasets.log"
  conda:
    "../envs/all_packages.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/figures_wo_positive_controls/power_to_detect_change_for_all_datasets.R"
