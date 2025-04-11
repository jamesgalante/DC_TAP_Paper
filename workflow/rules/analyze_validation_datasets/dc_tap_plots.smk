

# These plots take the CellRanger output and plot exploratory statistics
rule plot_screen_results:
  input:
    # What is number of unique genes per cell - by batch and cell type
    # What is number of gene UMIs per cell - by batch and cell type
    # For each gene (because ~300 and ~200) - what are the number of UMIs per cell - by batch and cell type
      # histogram of gene expression
    
    # What is the number of unique guides per cell - by batch and cell type
    # What is the number of guide UMIs per cell - by batch and cell type
      # histogram of guide expression
      
    # Would be cool to have a proportion plot of all UMIs (guide and gene) and then a proportion plot by specific gene or specific guide
    # Can plot for each cell, the proportion of reads to each gene or guide and then like categorize based on perturbation - can see visual change... idk

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
    
    

# Create a rule ranking (on xaxis) the encode score (~40 bins) and on the y axis the effect size 
# Plot this with all significant pairs and without significant pairs
  # Compare dc tap with other screens demonstrating that there are less false negatives


    


