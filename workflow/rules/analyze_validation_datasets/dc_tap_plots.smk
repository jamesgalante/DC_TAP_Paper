
# Create distance by effect size plots
rule distance_by_effect_size_plots:
  input:
    discovery_results = expand("results/process_validation_datasets/{sample}/differential_expression/results_run_discovery_analysis.rds", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    final_sceptre_object = expand("results/process_validation_datasets/{sample}/differential_expression/final_sceptre_object.rds", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    distances = expand("results/process_validation_datasets/{sample}/distances.tsv", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"])
  output:
    k562_distance_by_es = "results/analyze_validation_datasets/dc_tap_plots/distance_by_effect_size_plots/k562_distance_by_es.pdf",
    wtc11_distance_by_es = "results/analyze_validation_datasets/dc_tap_plots/distance_by_effect_size_plots/wtc11_distance_by_es.pdf"
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/distance_by_effect_size_plots.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/dc_tap_plots/distance_by_effect_size_plots.R"
    
rule power_simulation_plots:
  input:
    k562_power_output = "results/process_validation_datasets/K562_DC_TAP_Seq/power_analysis/output_0.13gStd_Sceptre_perCRE.tsv",
    wtc11_power_output = "results/process_validation_datasets/WTC11_DC_TAP_Seq/power_analysis/output_0.13gStd_Sceptre_perCRE.tsv",
    gasperini_power_output = "resources/analyze_validation_datasets/power_simulation_plots/output_0.13gStd_MAST_perCRE.tsv.gz",
    k562_tpm_file = "resources/analyze_validation_datasets/power_simulation_plots/TPM_from_NC_cells_only_P16_10x.txt"
  output:
    power_sim_results_line_plot = "results/analyze_validation_datasets/dc_tap_plots/power_plots/power_sim_results_line_plot.pdf",
    power_sim_results_bar_plot = "results/analyze_validation_datasets/dc_tap_plots/power_plots/power_sim_results_bar_plot.pdf"
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/power_simulation_plots.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/dc_tap_plots/power_simulation_plots.R"
    
rule effect_size_plots:
  input:
    experiment_summary_table = "results/analyze_validation_datasets/dc_tap_plots/summarized_categories.tsv"
  output:
    effect_size_boxplot = "results/analyze_validation_datasets/dc_tap_plots/effect_size_plots/effect_size_boxplot.pdf"
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/effect_size_plots.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/dc_tap_plots/effect_size_plots.R" 
    
# Create a rule to calculate the effect size distribution of screens
# basically showing that the median effect size is like 15% or something or wahtever it is (maybe this changes by TPM)
# And then comparing the number of cells needed to detect a 15% change between different screens (perturb seq style v dc tap)
      
# Create a rule ranking (on xaxis) the encode score (~40 bins) and on the y axis the effect size 
# Plot this with all significant pairs and without significant pairs
  # Compare dc tap with other screens demonstrating that there are less false negatives


    
# rule choosing_enhancer_category_thresholds:
#   input:
#     # Input the bed files and such
#   output:
#     # Plots for paper on threshold decisions
#   log: "results/analyze_validation_datasets/dc_tap_plots/logs/choosing_enhancer_category_thresholds.log"
#   conda:
#     "../../envs/analyze_crispr_screen.yml"
#   resources:
#     mem = "8G",
#     time = "1:00:00"
#   script:
#     "../../scripts/analyze_validation_datasets/dc_tap_plots/choosing_enhancer_category_thresholds.R"

# Create the enhancer categories based off chosen thresholds
# rule create_enhancer_categories:
#   input:
#     combined_validation = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_validation_expt_pred_merged_annot.txt",
#     combined_training = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_training_expt_pred_merged_annot.txt"
#   output:
#     labelled_combined_validation = "results/analyze_validation_datasets/dc_tap_plots/expt_pred_merged_annot/labelled_combined_validation_expt_pred_merged_annot.txt",
#     labelled_combined_training = "results/analyze_validation_datasets/dc_tap_plots/expt_pred_merged_annot/labelled_combined_training_expt_pred_merged_annot.txt",
#   log: "results/analyze_validation_datasets/dc_tap_plots/logs/create_enhancer_categories_and_plot.log"
#   conda:
#     "../../envs/analyze_crispr_screen.yml"
#   resources:
#     mem = "8G",
#     time = "1:00:00"
#   script:
#     "../../scripts/analyze_validation_datasets/dc_tap_plots/create_enhancer_categories_and_plot.R"

# Rule to create the proportion dc tap proportion plots without the "control" genes
# rule dc_tap_proportions_plots:
#   input:
#     labelled_combined_validation = "results/analyze_validation_datasets/dc_tap_plots/expt_pred_merged_annot/labelled_combined_validation_expt_pred_merged_annot.txt",
#     labelled_combined_training = "results/analyze_validation_datasets/dc_tap_plots/expt_pred_merged_annot/labelled_combined_training_expt_pred_merged_annot.txt",
#     guide_targets = expand("results/process_validation_datasets/{sample}/guide_targets.tsv", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"])
#   output:
#     proportions_plot_all = "results/analyze_validation_datasets/dc_tap_plots/proportions_plots/proportions_plot_all.pdf",
#     proportions_plot_pos = "results/analyze_validation_datasets/dc_tap_plots/proportions_plots/proportions_plot_pos.pdf",
#     fold_change_plot = "results/analyze_validation_datasets/dc_tap_plots/proportions_plots/fold_change_plot.pdf",
#     chi_squared_results = "results/analyze_validation_datasets/dc_tap_plots/proportions_plots/chi_squared_results.tsv"
#   log: "results/analyze_validation_datasets/dc_tap_plots/logs/dc_tap_proportions_plots.log"
#   conda:
#     "../../envs/analyze_crispr_screen.yml"
#   resources:
#     mem = "32G",
#     time = "2:00:00"
#   script:
#     "../../scripts/analyze_validation_datasets/dc_tap_plots/dc_tap_proportions_plots.R"
  

